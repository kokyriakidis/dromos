"""`dromos index`: a GAM made searchable by reference position, as pgphase reads it, in one pass and one process. It
writes the file `pgcat index-gaf --slim --ref-map` writes from the GAM's GAF (the same lines, blocks and index).

Each alignment is placed on the reference by the nodes of its walk that lie on the reference sample's paths
(refmap.mojo), written as `contig  start+1  end  <its 12 GAF columns>`, sorted by contig name, start, end and line,
bgzipped in blocks of exactly BLOCK bytes and tabix-indexed. Header lines `#max_span <contig> <bases>` give the most
graph bases any alignment on each contig spans.

While the GAM streams, the placed lines gather by bucket (BUCKET reference bases of a contig), each bucket's lines
compressed by zstd a FRAME at a time beside the GAM's own tasks. Then a window of buckets at a time is inflated and
sorted, each bucket on a thread of its own, laid out, and its whole blocks deflated in parallel and written in
order, the index built as they go.
"""

from std.ffi import OwnedDLHandle, external_call
from std.memory import unsafe_memcpy, unsafe_memmove
from std.os import getenv
from std.sys import stderr
from std.time import perf_counter_ns

from dromos.bgzf import BLOCK, Deflater, eof_block, load_libdeflate, load_zstd, zstd_compress, zstd_decompress
from dromos.gam import Alignment, gaf_line
from dromos.output import put_padded
from dromos.parallel import ParallelTask, parallel_for
from dromos.refmap import NO_ANCHOR, Anchor, Placed, RefMap, reason_name
from dromos.stream import MessageSink, imm_untracked, mut_untracked
from dromos.tabix import TbiBuilder

comptime BUCKET_SHIFT = 20
"""A bucket is 2^20 reference bases of a contig (pgcat's BUCKET)."""
comptime FRAME = 1 << 20
"""A bucket's lines are compressed once this many bytes of them gather."""
comptime ZSTD_LEVEL = 1
comptime HEADER = 24
"""Bytes before each line a slice keeps: contig rank u32, start u64, end u64, line length u32."""
comptime BLOCKS_PER_TASK = 16

# Why an alignment has no place, counted per slice: input lines, placed, unmapped, unplaced (then by reason, in
# refmap.mojo's order) and malformed.
comptime C_LINES = 0
comptime C_PLACED = 1
comptime C_UNMAPPED = 2
comptime C_UNPLACED = 3
comptime C_WHY = 4
comptime C_MALFORMED = 8
comptime N_COUNTS = 9

comptime _SKIP = 0
comptime _PLACED = 1
comptime _UNMAPPED = 2
comptime _UNPLACED = 3
comptime _MALFORMED = 4


@always_inline
def _get32(p: ImmPointer[UInt8, _], at: Int) -> Int:
    return (
        Int(p[unsafe_offset=at])
        | Int(p[unsafe_offset=at + 1]) << 8
        | Int(p[unsafe_offset=at + 2]) << 16
        | Int(p[unsafe_offset=at + 3]) << 24
    )


@always_inline
def _get64(p: ImmPointer[UInt8, _], at: Int) -> Int:
    return _get32(p, at) | _get32(p, at + 4) << 32


@always_inline
def _put(mut out: List[UInt8], v: Int, bytes: Int):
    for k in range(bytes):
        out.append(UInt8((v >> (8 * k)) & 0xFF))


def _parse_u64(p: ImmPointer[UInt8, _], at: Int, end: Int) -> Int:
    """The decimal number `p[at:end]`, or -1 (pgcat's parse_u64: 1 to 19 digits, nothing else)."""
    if end <= at or end - at > 19:
        return -1
    var v = 0
    for i in range(at, end):
        var c = Int(p[unsafe_offset=i]) - ord("0")
        if c < 0 or c > 9:
            return -1
        v = v * 10 + c
    return v if v >= 0 else -1


def _parse_walk(p: ImmPointer[UInt8, _], at: Int, end: Int, mut walk: List[Int]) -> Bool:
    """A GAF walk (`>12<34`) as steps `node << 1 | reverse`; false if it is not one."""
    walk.clear()
    var i = at
    while i < end:
        var c = p[unsafe_offset=i]
        var reverse: Int
        if c == UInt8(ord(">")):
            reverse = 0
        elif c == UInt8(ord("<")):
            reverse = 1
        else:
            return False
        i += 1
        var s = i
        while i < end and p[unsafe_offset=i] != UInt8(ord(">")) and p[unsafe_offset=i] != UInt8(ord("<")):
            i += 1
        var node = _parse_u64(p, s, i)
        if node < 0:
            return False
        walk.append(node << 1 | reverse)
    return True


struct _Slice(Movable):
    """A slice's own state: its scratch, what it placed this step, and its counts over the whole GAM."""

    var alignment: Alignment
    var line: List[UInt8]
    var walk: List[Int]
    var anchors: List[Anchor]
    var recs: List[UInt8]
    var counts: List[Int]
    var span: List[Int]

    def __init__(out self, contigs: Int):
        self.alignment = Alignment()
        self.line = List[UInt8]()
        self.walk = List[Int]()
        self.anchors = List[Anchor]()
        self.recs = List[UInt8]()
        self.counts = List[Int](length=N_COUNTS, fill=0)
        self.span = List[Int](length=contigs, fill=0)


struct _Bucket(Copyable, Movable):
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


struct _Job(Movable):
    """A FRAME of a bucket's lines, to compress."""

    var rank: Int
    var bucket: Int
    var raw: List[UInt8]
    var size: Int
    var frame: List[UInt8]

    def __init__(out self, rank: Int, bucket: Int, var raw: List[UInt8]):
        self.rank = rank
        self.bucket = bucket
        self.size = len(raw)
        self.raw = raw^
        self.frame = List[UInt8]()


@fieldwise_init
struct _Rec(TrivialRegisterPassable):
    """A line of a bucket being sorted: its key, where it is in the bucket's bytes, and its output length."""

    var start: Int
    var end: Int
    var off: Int
    var len: Int
    var out_len: Int


struct _Slot(Movable):
    """A bucket of the window being written: its lines inflated, sorted, and where its text goes in the window."""

    var rank: Int
    var bucket: Int
    var raw: List[UInt8]
    var recs: List[_Rec]
    var text_len: Int
    var at: Int

    def __init__(out self):
        self.rank = 0
        self.bucket = 0
        self.raw = List[UInt8]()
        self.recs = List[_Rec]()
        self.text_len = 0
        self.at = 0


@always_inline
def _digits(v: Int) -> Int:
    var n = 1
    var bound = 10
    while v >= bound and n < 19:
        n += 1
        bound *= 10
    return n


struct IndexSink(MessageSink, Movable):
    """GAM messages placed and gathered by bucket (see the module's notes)."""

    var map: RefMap
    var zstd: OwnedDLHandle
    var rank: List[Int]
    """Each contig's rank: contigs in name order."""
    var by_rank: List[Int]
    var states: List[_Slice]
    var buckets: List[List[_Bucket]]
    """By contig rank, then bucket."""
    var compressing: List[_Job]
    """Frames the running step compresses, beside its slices."""
    var queued: List[_Job]
    """Frames for the step after it."""

    def __init__(out self, var map: RefMap, threads: Int) raises:
        var n = len(map.contigs)
        # Output order is by contig name (bytes), then position, then the line itself.
        self.by_rank = List[Int](capacity=n)
        for c in range(n):
            var k = len(self.by_rank)
            self.by_rank.append(c)
            while k > 0 and _name_less(map.contigs[self.by_rank[k]], map.contigs[self.by_rank[k - 1]]):
                self.by_rank.swap_elements(k, k - 1)
                k -= 1
        self.rank = List[Int](length=n, fill=0)
        for r in range(n):
            self.rank[self.by_rank[r]] = r
        self.map = map^
        self.zstd = load_zstd()
        self.states = List[_Slice]()
        for _ in range(4 * threads):
            self.states.append(_Slice(n))
        self.buckets = List[List[_Bucket]](length=n, fill=List[_Bucket]())
        self.compressing = List[_Job]()
        self.queued = List[_Job]()

    def slices(self) -> Int:
        return len(self.states)

    def take(
        mut self,
        slice: Int,
        data: ImmPointer[UInt8, ImmUntrackedOrigin],
        messages: ImmPointer[Int, ImmUntrackedOrigin],
        first: Int,
        last: Int,
    ) raises:
        var st = mut_untracked(Pointer(to=self.states[slice]))
        st[].recs.clear()
        var placed = Placed(0, 0, 0)
        for m in range(first, last):
            st[].alignment.decode(data, messages[unsafe_offset=2 * m], messages[unsafe_offset=2 * m + 1])
            st[].line.clear()
            gaf_line[False](data, st[].alignment, self.map, st[].line)
            var p = st[].line.unsafe_ptr().as_imm()
            var n = len(st[].line) - 1
            var why = 0
            var aligned = 0
            var kind = self._classify(p, n, st[].walk, st[].anchors, placed, why, aligned)
            if kind == _SKIP:
                continue
            st[].counts[C_LINES] += 1
            if kind == _PLACED:
                st[].counts[C_PLACED] += 1
                var r = self.rank[placed.contig]
                st[].span[r] = max(st[].span[r], aligned)
                _put(st[].recs, r, 4)
                _put(st[].recs, placed.start, 8)
                _put(st[].recs, placed.end, 8)
                _put(st[].recs, n, 4)
                st[].recs.extend(Span(unsafe_ptr=p, length=n))
            elif kind == _UNMAPPED:
                st[].counts[C_UNMAPPED] += 1
            elif kind == _UNPLACED:
                st[].counts[C_UNPLACED] += 1
                st[].counts[C_WHY + why - NO_ANCHOR] += 1
            else:
                st[].counts[C_MALFORMED] += 1

    def _classify(
        self,
        p: ImmPointer[UInt8, _],
        n: Int,
        mut walk: List[Int],
        mut anchors: List[Anchor],
        mut placed: Placed,
        mut why: Int,
        mut aligned: Int,
    ) -> Int:
        """pgcat's classify, for one GAF line `p[:n]` (no newline)."""
        if n == 0 or p[unsafe_offset=0] == UInt8(ord("#")) or p[unsafe_offset=0] == UInt8(ord("@")):
            return _SKIP
        # The 12 columns' starts, and the end of the 12th.
        var starts = SIMD[DType.int64, 16](0)
        var ends = SIMD[DType.int64, 16](0)
        var field = 0
        var s = 0
        for i in range(n + 1):
            if i == n or p[unsafe_offset=i] == UInt8(ord("\t")):
                starts[field] = Int64(s)
                ends[field] = Int64(i)
                field += 1
                s = i + 1
                if field == 12:
                    break
        if field < 12:
            return _MALFORMED
        var strand_at = Int(starts[4])
        var strand_len = Int(ends[4]) - strand_at
        if strand_len == 1 and p[unsafe_offset=strand_at] == UInt8(ord("*")):
            return _UNMAPPED
        if strand_len != 1 or (
            p[unsafe_offset=strand_at] != UInt8(ord("+")) and p[unsafe_offset=strand_at] != UInt8(ord("-"))
        ):
            return _MALFORMED
        var ts = _parse_u64(p, Int(starts[7]), Int(ends[7]))
        var te = _parse_u64(p, Int(starts[8]), Int(ends[8]))
        if ts < 0 or te < 0:
            return _MALFORMED
        if not _parse_walk(p, Int(starts[5]), Int(ends[5]), walk):
            return _MALFORMED
        why = self.map.project(walk, ts, te, anchors, placed)
        if why != 0:
            return _UNPLACED
        aligned = te - ts
        return _PLACED

    def side_tasks(self) -> Int:
        return len(self.compressing)

    def side_task(mut self, index: Int) raises:
        var job = mut_untracked(Pointer(to=self.compressing[index]))
        job[].frame = zstd_compress(self.zstd, job[].raw.unsafe_ptr(), len(job[].raw), ZSTD_LEVEL)
        job[].raw = List[UInt8]()

    def after_step(mut self, slices: Int) raises:
        while len(self.compressing) > 0:
            var job = self.compressing.pop()
            var rank = job.rank
            var b = job.bucket
            var size = job.size
            var frame = List[UInt8]()
            swap(frame, job.frame)
            self.buckets[rank][b].frames.append(frame^)
            self.buckets[rank][b].sizes.append(size)
        for s in range(slices):
            var recs = imm_untracked(self.states[s].recs.unsafe_ptr())
            var total = len(self.states[s].recs)
            var at = 0
            while at < total:
                var r = _get32(recs, at)
                var b = _get64(recs, at + 4) >> BUCKET_SHIFT
                var length = _get32(recs, at + 20)
                if len(self.buckets[r]) <= b:
                    self.buckets[r].resize(b + 1, _Bucket())
                ref bucket = self.buckets[r][b]
                bucket.pending.extend(Span(unsafe_ptr=recs.unsafe_offset(at + 4), length=HEADER - 4 + length))
                bucket.records += 1
                if len(bucket.pending) >= FRAME:
                    var raw = List[UInt8]()
                    swap(raw, bucket.pending)
                    self.queued.append(_Job(r, b, raw^))
                at += HEADER + length
        swap(self.compressing, self.queued)

    def finish(mut self, out_path: String, level: Int, threads: Int) raises:
        """Writes the index: `out_path` and `out_path`.tbi, each through a temporary name."""
        var t0 = Int(perf_counter_ns())
        var n = len(self.map.contigs)
        var counts = List[Int](length=N_COUNTS, fill=0)
        var max_span = List[Int](length=n, fill=0)
        for si in range(len(self.states)):
            ref st = self.states[si]
            for k in range(N_COUNTS):
                counts[k] += st.counts[k]
            for r in range(n):
                max_span[r] = max(max_span[r], st.span[r])
        var why = String()
        for k in range(4):
            if k > 0:
                why += ", "
            why += String(t"{reason_name(NO_ANCHOR + k)} {counts[C_WHY + k]}")
        print(
            (
                t"[dromos] {counts[C_LINES]} alignments: {counts[C_PLACED]} placed, {counts[C_UNMAPPED]} unmapped,"
                t" {counts[C_UNPLACED]} without a unique reference interval ({why}), {counts[C_MALFORMED]} malformed"
            ),
            file=stderr,
        )
        # The buckets with lines, in output order.
        var order = List[Tuple[Int, Int]]()
        for r in range(n):
            for b in range(len(self.buckets[r])):
                if self.buckets[r][b].records > 0:
                    order.append((r, b))
        var lib = load_libdeflate()
        var w = _Writer(out_path, level)
        # The longest alignment on each contig, for pgphase: a site whose every allele spans more graph bases cannot
        # be traversed whole by any read. Every contig is listed, 0 for one without alignments.
        var header = String()
        for r in range(n):
            header += String(t"#max_span\t{self.map.contigs[self.by_rank[r]]}\t{max_span[r]}\n")
        w.header(header)
        var window = 4 * threads
        var slots = List[_Slot]()
        for _ in range(window):
            slots.append(_Slot())
        var bad = List[Int](length=1, fill=0)
        var gather_ns = 0
        var layout_ns = 0
        var index_ns = 0
        var serial_ns = 0
        var start = 0
        while start < len(order):
            var count = min(window, len(order) - start)
            for k in range(count):
                slots[k].rank = order[start + k][0]
                slots[k].bucket = order[start + k][1]
            var gather = _WindowTask(
                0,
                mut_untracked(Pointer(to=self)),
                mut_untracked(slots.unsafe_ptr()),
                mut_untracked(Pointer(to=w)),
                mut_untracked(bad.unsafe_ptr()),
            )
            var g0 = Int(perf_counter_ns())
            parallel_for(gather, count, threads)
            var g1 = Int(perf_counter_ns())
            gather_ns += g1 - g0
            if bad[0] != 0:
                raise Error("a bucket's lines did not inflate")
            # Where each bucket's text goes, after what the last window left of its last block.
            var at = len(w.buf)
            for k in range(count):
                slots[k].at = at
                at += slots[k].text_len
            w.buf.resize(unsafe_uninit_length=at)
            var lay_out = _WindowTask(
                1,
                mut_untracked(Pointer(to=self)),
                mut_untracked(slots.unsafe_ptr()),
                mut_untracked(Pointer(to=w)),
                mut_untracked(bad.unsafe_ptr()),
            )
            var l0 = Int(perf_counter_ns())
            parallel_for(lay_out, count, threads)
            var l1 = Int(perf_counter_ns())
            layout_ns += l1 - l0
            for k in range(count):
                var contig = self.by_rank[slots[k].rank]
                var end = w.base + slots[k].at
                for rec in slots[k].recs:
                    end += rec.out_len
                    w.pending.append(_Tab(contig, rec.start, rec.end, end))
                slots[k].raw = List[UInt8]()
                slots[k].recs = List[_Rec]()
            var i0 = Int(perf_counter_ns())
            w.write_blocks(lib, len(w.buf) // BLOCK, threads)
            var i1 = Int(perf_counter_ns())
            w.index_written(self.map.contigs)
            index_ns += Int(perf_counter_ns()) - i1
            serial_ns += i0 - l1
            start += count
        w.write_blocks(lib, (len(w.buf) + BLOCK - 1) // BLOCK, threads)
        w.index_written(self.map.contigs)
        w.close(lib)
        var t1 = Int(perf_counter_ns())
        print(
            (
                t"[dromos] {out_path} and {out_path}.tbi: {counts[C_PLACED]} alignments ({Float64(t1 - t0) * 1e-9} s"
                t" sorting and writing)"
            ),
            file=stderr,
        )
        if Bool(getenv("DROMOS_TIMING")):
            print(
                (
                    t"dromos timing (s): inflate+sort {Float64(gather_ns) * 1e-9}, lay out {Float64(layout_ns) * 1e-9},"
                    t" deflate {Float64(w.deflate_ns) * 1e-9}, write {Float64(w.write_ns) * 1e-9}, tabix"
                    t" {Float64(index_ns) * 1e-9}, between {Float64(serial_ns) * 1e-9}"
                ),
                file=stderr,
            )

    def gather(mut self, mut slot: _Slot) raises:
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
        slot.recs = List[_Rec](capacity=bucket.records)
        bucket.frames = List[List[UInt8]]()
        bucket.sizes = List[Int]()
        bucket.pending = List[UInt8]()
        var name_len = self.map.contigs[self.by_rank[slot.rank]].byte_length()
        var p = imm_untracked(slot.raw.unsafe_ptr())
        at = 0
        slot.text_len = 0
        while at < total:
            var s = _get64(p, at)
            var e = _get64(p, at + 8)
            var length = _get32(p, at + 16)
            # contig \t start+1 \t end \t line \n
            var out_len = name_len + 1 + _digits(s + 1) + 1 + _digits(e) + 1 + length + 1
            slot.recs.append(_Rec(s, e, at + 20, length, out_len))
            slot.text_len += out_len
            at += 20 + length

        def less(a: _Rec, b: _Rec) {p} -> Bool:
            if a.start != b.start:
                return a.start < b.start
            if a.end != b.end:
                return a.end < b.end
            var c = external_call["memcmp", Int32](p.unsafe_offset(a.off), p.unsafe_offset(b.off), min(a.len, b.len))
            if c != 0:
                return c < 0
            return a.len < b.len

        sort(Span(slot.recs), less)

    def lay_out(self, slot: _Slot, dst: MutPointer[UInt8, _]):
        """Slot `slot`'s lines, as `contig  start+1  end  line`, at `dst`."""
        var name = self.map.contigs[self.by_rank[slot.rank]].as_bytes()
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


@always_inline
def _name_less(a: String, b: String) -> Bool:
    """`a` before `b` by their bytes (Rust's ordering of byte strings)."""
    var x = a.as_bytes()
    var y = b.as_bytes()
    for i in range(min(len(x), len(y))):
        if x[i] != y[i]:
            return x[i] < y[i]
    return len(x) < len(y)


@fieldwise_init
struct _Tab(TrivialRegisterPassable):
    """A line written, waiting for the block holding its end: its contig, interval, and its end in the stream."""

    var contig: Int
    var beg: Int
    var end: Int
    var stream_end: Int


struct _Writer(Movable):
    """The BGZF output, cut into blocks of exactly BLOCK bytes, so the virtual offset of any byte follows from its
    stream offset and the blocks' addresses (pgcat's IndexedWriter): lines are indexed once the block holding their
    end is written."""

    var path: String
    var tmp: String
    var file: FileHandle
    var level: Int
    var buf: List[UInt8]
    var base: Int
    """The stream offset of buf[0], a block boundary."""
    var written: Int
    """Compressed bytes written: the address of the next block."""
    var addrs: List[Int]
    """Every block's address, by block number."""
    var pending: List[_Tab]
    var head: Int
    var index: TbiBuilder
    var outs: List[List[UInt8]]
    var deflate_ns: Int
    var write_ns: Int

    def __init__(out self, path: String, level: Int) raises:
        self.path = path
        self.tmp = path + ".tmp"
        self.file = open(self.tmp, "w")
        self.level = level
        self.buf = List[UInt8]()
        self.base = 0
        self.written = 0
        self.addrs = List[Int]()
        self.pending = List[_Tab]()
        self.head = 0
        self.index = TbiBuilder()
        self.outs = List[List[UInt8]]()
        self.deflate_ns = 0
        self.write_ns = 0

    def header(mut self, text: String) raises:
        """Header lines (each starting with '#', with its newline), before any line: tabix skips them, and the index
        starts after them."""
        self.buf.extend(text.as_bytes())
        if len(self.buf) >= BLOCK:
            raise Error("header longer than a BGZF block")
        # Block 0 starts at address 0: the first line starts at (0, length).
        self.index.set_first_offset(len(self.buf))

    def write_blocks(mut self, lib: OwnedDLHandle, blocks: Int, threads: Int) raises:
        """Deflates and writes the first `blocks` blocks of buf (the last may be partial, when finishing), keeping
        the rest."""
        if blocks == 0:
            return
        var take = min(blocks * BLOCK, len(self.buf))
        var tasks = (blocks + BLOCKS_PER_TASK - 1) // BLOCKS_PER_TASK
        while len(self.outs) < tasks:
            self.outs.append(List[UInt8]())
        var failed = List[Int](length=1, fill=0)
        var task = _DeflateTask(
            Pointer(to=lib).as_imm(),
            self.level,
            imm_untracked(self.buf.unsafe_ptr()),
            take,
            mut_untracked(self.outs.unsafe_ptr()),
            mut_untracked(failed.unsafe_ptr()),
        )
        var d0 = Int(perf_counter_ns())
        parallel_for(task, tasks, threads)
        var d1 = Int(perf_counter_ns())
        self.deflate_ns += d1 - d0
        if failed[0] != 0:
            raise Error("BGZF: a block did not deflate")
        for t in range(tasks):
            var p = self.outs[t].unsafe_ptr().as_imm()
            var at = 0
            while at < len(self.outs[t]):
                self.addrs.append(self.written)
                var size = (Int(p[unsafe_offset=at + 16]) | Int(p[unsafe_offset=at + 17]) << 8) + 1
                self.written += size
                at += size
            self.file.write_bytes(Span(self.outs[t]))
        var rest = len(self.buf) - take
        if rest > 0:
            unsafe_memmove(
                dest=mut_untracked(self.buf.unsafe_ptr()),
                src=imm_untracked(self.buf.unsafe_ptr()).unsafe_offset(take),
                count=rest,
            )
        self.buf.resize(unsafe_uninit_length=rest)
        self.base += take
        self.write_ns += Int(perf_counter_ns()) - d1

    def _voffset(self, s: Int) -> Int:
        """The virtual offset of stream offset `s`, in a written block or at the start of the next: a block's end
        reads as the start of the next block, as htslib reports it after reading a line."""
        if s == self.base:
            return self.written << 16
        return self.addrs[s // BLOCK] << 16 | (s % BLOCK)

    def index_written(mut self, contigs: List[String]) raises:
        """Indexes the lines whose ends are written."""
        while self.head < len(self.pending):
            var t = self.pending[self.head]
            if t.stream_end > self.base:
                break
            self.index.push(t.contig, contigs[t.contig], t.beg, t.end, self._voffset(t.stream_end))
            self.head += 1
        if self.head == len(self.pending):
            self.pending.clear()
            self.head = 0

    def close(mut self, lib: OwnedDLHandle) raises:
        """The EOF block, the index, and both files under their names."""
        var final_offset = self.written << 16
        var eof = List[UInt8]()
        eof_block(eof)
        self.file.write_bytes(Span(eof))
        self.file.close()
        self.index.finish(final_offset)
        var tbi = self.path + ".tbi"
        self.index.write(lib, tbi + ".tmp")
        _rename(self.tmp, self.path)
        _rename(tbi + ".tmp", tbi)


def _rename(src: String, dst: String) raises:
    var a = List[UInt8](src.as_bytes())
    a.append(0)
    var b = List[UInt8](dst.as_bytes())
    b.append(0)
    if external_call["rename", Int32](a.unsafe_ptr(), b.unsafe_ptr()) != 0:
        raise Error(t"cannot rename {src} to {dst}")


@fieldwise_init
struct _DeflateTask[lib_origin: ImmOrigin](ParallelTask):
    """Task t deflates blocks t * BLOCKS_PER_TASK onward of `data[:length]` into outs[t]."""

    var lib: ImmPointer[OwnedDLHandle, Self.lib_origin]
    var level: Int
    var data: ImmPointer[UInt8, ImmUntrackedOrigin]
    var length: Int
    var outs: MutPointer[List[UInt8], MutUntrackedOrigin]
    var failed: MutPointer[Int, MutUntrackedOrigin]

    def run(self, index: Int):
        ref out = self.outs[unsafe_offset=index]
        out.clear()
        try:
            var d = Deflater(self.lib[], self.level)
            var at = index * BLOCKS_PER_TASK * BLOCK
            var end = min(self.length, at + BLOCKS_PER_TASK * BLOCK)
            while at < end:
                var n = min(BLOCK, end - at)
                d.block(self.lib[], self.data.unsafe_offset(at), n, out)
                at += n
            d.free(self.lib[])
        except:
            self.failed[] = 1


@fieldwise_init
struct _WindowTask(ParallelTask):
    """Stage 0: slot k's bucket inflated and sorted; stage 1: its lines laid out in the writer's buffer."""

    var stage: Int
    var sink: MutPointer[IndexSink, MutUntrackedOrigin]
    var slots: MutPointer[_Slot, MutUntrackedOrigin]
    var writer: MutPointer[_Writer, MutUntrackedOrigin]
    var bad: MutPointer[Int, MutUntrackedOrigin]

    def run(self, index: Int):
        if self.stage == 0:
            try:
                self.sink[].gather(self.slots[unsafe_offset=index])
            except:
                self.bad[] = 1
            return
        self.sink[].lay_out(self.slots[unsafe_offset=index], self.writer[].buf.unsafe_ptr())
