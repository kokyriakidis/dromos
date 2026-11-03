"""Placed lines staged by bucket while the input streams, for the final sort: each bucket (BUCKET_SHIFT's reference
bases of a contig) keeps its lines, compressed by zstd a FRAME at a time. A line is staged with its start u64, end
u64 and length u32; a slice hands its lines over with the contig's rank before them (HEADER bytes in all)."""

from std.ffi import OwnedDLHandle, external_call
from std.memory import unsafe_memcpy

from dromos.bgzf import load_zstd, zstd_compress, zstd_decompress
from dromos.output import decimal_digits, put_padded
from dromos.stream import imm_untracked, mut_untracked

comptime BUCKET_SHIFT = 20
"""A bucket is 2^20 reference bases of a contig (pgcat's BUCKET)."""
comptime FRAME = 1 << 20
"""A bucket's lines are compressed once this many bytes of them gather."""
comptime ZSTD_LEVEL = 1
comptime HEADER = 24
"""Bytes before each line a slice keeps: contig rank u32, start u64, end u64, line length u32."""


@always_inline
def get32(p: ImmPointer[UInt8, _], at: Int) -> Int:
    return (
        Int(p[unsafe_offset=at])
        | Int(p[unsafe_offset=at + 1]) << 8
        | Int(p[unsafe_offset=at + 2]) << 16
        | Int(p[unsafe_offset=at + 3]) << 24
    )


@always_inline
def get64(p: ImmPointer[UInt8, _], at: Int) -> Int:
    return get32(p, at) | get32(p, at + 4) << 32


@always_inline
def store(p: MutPointer[UInt8, _], v: Int, bytes: Int):
    for k in range(bytes):
        p[unsafe_offset=k] = UInt8((v >> (8 * k)) & 0xFF)


def put_record(mut out: List[UInt8], at: Int, rank: Int, start: Int, end: Int):
    """Fills the HEADER bytes at `at` of `out` before the line that follows them to the end of `out`."""
    var h = out.unsafe_ptr().unsafe_offset(at)
    store(h, rank, 4)
    store(h.unsafe_offset(4), start, 8)
    store(h.unsafe_offset(12), end, 8)
    store(h.unsafe_offset(20), len(out) - at - HEADER, 4)


struct Bucket(Copyable, Movable):
    """A bucket's lines while the GAM streams: zstd frames (and each one's size inflated), then those not yet
    compressed. Each line is kept with its start u64, end u64 and length u32."""

    var frames: List[List[UInt8]]
    var sizes: List[Int]
    var pending: List[UInt8]
    var records: Int

    def __init__(out self):
        self.frames = List[List[UInt8]]()
        self.sizes = List[Int]()
        self.pending = List[UInt8]()
        self.records = 0


@fieldwise_init
struct Rec(TrivialRegisterPassable):
    """A line of a bucket being sorted: its key, where it is in the bucket's bytes, and its output length."""

    var start: Int
    var end: Int
    var off: Int
    var len: Int
    var out_len: Int


struct Slot(Movable):
    """A bucket of the window being written: its lines inflated, sorted, and where its text goes in the window."""

    var rank: Int
    var bucket: Int
    var raw: List[UInt8]
    var recs: List[Rec]
    var text_len: Int
    var at: Int

    def __init__(out self):
        self.rank = 0
        self.bucket = 0
        self.raw = List[UInt8]()
        self.recs = List[Rec]()
        self.text_len = 0
        self.at = 0


struct Staging(Movable):
    """The buckets, by contig rank then bucket; the contigs' names by rank."""

    var buckets: List[List[Bucket]]
    var names: List[String]
    var zstd: OwnedDLHandle
    var shares: Int
    """How many tasks the buckets are shared among when lines are handed over (`take_share`)."""

    def __init__(out self, var names: List[String], shares: Int) raises:
        self.buckets = List[List[Bucket]](length=len(names), fill=List[Bucket]())
        self.names = names^
        self.zstd = load_zstd()
        self.shares = shares

    def make_room(mut self, recs: ImmPointer[UInt8, _], total: Int):
        """Makes the buckets the records `recs[:total]` need (on one thread, before they are handed over)."""
        var at = 0
        while at < total:
            var r = get32(recs, at)
            var b = get64(recs, at + 4) >> BUCKET_SHIFT
            if len(self.buckets[r]) <= b:
                self.buckets[r].resize(b + 1, Bucket())
            at += HEADER + get32(recs, at + 20)

    def take_share(mut self, share: Int, recs: ImmPointer[UInt8, _], total: Int) raises:
        """Takes the records of `recs[:total]` in share `share`'s buckets (bucket b of rank r: (r << 12 | b) mod
        shares), compressing a bucket's lines as a frame once FRAME bytes of them gather. Shares run at once."""
        var at = 0
        while at < total:
            var r = get32(recs, at)
            var b = get64(recs, at + 4) >> BUCKET_SHIFT
            var length = get32(recs, at + 20)
            if (r << 12 | b) % self.shares == share:
                var bucket = mut_untracked(Pointer(to=self.buckets[r][b]))
                bucket[].pending.extend(Span(unsafe_ptr=recs.unsafe_offset(at + 4), length=HEADER - 4 + length))
                bucket[].records += 1
                if len(bucket[].pending) >= FRAME:
                    bucket[].frames.append(
                        zstd_compress(self.zstd, bucket[].pending.unsafe_ptr(), len(bucket[].pending), ZSTD_LEVEL)
                    )
                    bucket[].sizes.append(len(bucket[].pending))
                    bucket[].pending = List[UInt8]()
            at += HEADER + length

    def order(self) -> List[Tuple[Int, Int]]:
        """The buckets with lines, in output order: by contig rank, then position."""
        var order = List[Tuple[Int, Int]]()
        for r in range(len(self.buckets)):
            for b in range(len(self.buckets[r])):
                if self.buckets[r][b].records > 0:
                    order.append((r, b))
        return order^

    def sort(mut self, mut slot: Slot) raises:
        """Slot `slot`'s bucket inflated into one buffer and its lines sorted (pgcat's order: start, end, line), each
        with its output length; the bucket's memory freed."""
        ref bucket = self.buckets[slot.rank][slot.bucket]
        var total = len(bucket.pending)
        for s in bucket.sizes:
            total += s
        slot.raw.resize(unsafe_uninit_length=total)
        var at = 0
        for f in range(len(bucket.frames)):
            zstd_decompress(self.zstd, Span(bucket.frames[f]), slot.raw.unsafe_ptr().unsafe_offset(at), bucket.sizes[f])
            at += bucket.sizes[f]
        if len(bucket.pending) > 0:
            unsafe_memcpy(
                dest=slot.raw.unsafe_ptr().unsafe_offset(at), src=bucket.pending.unsafe_ptr(), count=len(bucket.pending)
            )
        slot.recs = List[Rec](capacity=bucket.records)
        bucket.frames = List[List[UInt8]]()
        bucket.sizes = List[Int]()
        bucket.pending = List[UInt8]()
        var name_len = self.names[slot.rank].byte_length()
        var p = imm_untracked(slot.raw.unsafe_ptr())
        at = 0
        slot.text_len = 0
        while at < total:
            var s = get64(p, at)
            var e = get64(p, at + 8)
            var length = get32(p, at + 16)
            # contig \t start+1 \t end \t line \n
            var out_len = name_len + 1 + decimal_digits(s + 1) + 1 + decimal_digits(e) + 1 + length + 1
            slot.recs.append(Rec(s, e, at + 20, length, out_len))
            slot.text_len += out_len
            at += 20 + length

        def less(a: Rec, b: Rec) {p} -> Bool:
            if a.start != b.start:
                return a.start < b.start
            if a.end != b.end:
                return a.end < b.end
            var c = external_call["memcmp", Int32](p.unsafe_offset(a.off), p.unsafe_offset(b.off), min(a.len, b.len))
            if c != 0:
                return c < 0
            return a.len < b.len

        sort(Span(slot.recs), less)

    def lay_out(self, slot: Slot, dst: MutPointer[UInt8, _]):
        """Slot `slot`'s lines, as `contig  start+1  end  line`, at `dst`."""
        var name = self.names[slot.rank].as_bytes()
        var src = slot.raw.unsafe_ptr()
        var at = slot.at
        for rec in slot.recs:
            var o = dst.unsafe_offset(at)
            unsafe_memcpy(dest=o, src=name.unsafe_ptr(), count=len(name))
            var k = len(name)
            o[unsafe_offset=k] = UInt8(ord("\t"))
            k += 1
            k += put_padded(o.unsafe_offset(k), rec.start + 1, 0)
            o[unsafe_offset=k] = UInt8(ord("\t"))
            k += 1
            k += put_padded(o.unsafe_offset(k), rec.end, 0)
            o[unsafe_offset=k] = UInt8(ord("\t"))
            k += 1
            unsafe_memcpy(dest=o.unsafe_offset(k), src=src.unsafe_offset(rec.off), count=rec.len)
            k += rec.len
            o[unsafe_offset=k] = UInt8(ord("\n"))
            at += rec.out_len
