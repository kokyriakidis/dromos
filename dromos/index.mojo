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
from dromos.gam import Alignment, slim_line
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


def _append_slim(mut out: List[UInt8], p: ImmPointer[UInt8, _], n: Int):
    """pgcat's append_line with slim: the GAF line `p[:n]`'s first 12 columns, then of the rest only its SLIM_TAGS
    (rq:i:, the rescue's MAPQ; pu:Z:, why a line had no place; rs:Z:, where the rescue put it), in their order."""
    var field = 0
    var s = 0
    for i in range(n + 1):
        if i < n and p[unsafe_offset=i] != UInt8(ord("\t")):
            continue
        var length = i - s
        if field < 12:
            if field > 0:
                out.append(UInt8(ord("\t")))
            out.extend(Span(unsafe_ptr=p.unsafe_offset(s), length=length))
        elif length >= 5 and p[unsafe_offset=s + 2] == UInt8(ord(":")) and p[unsafe_offset=s + 4] == UInt8(ord(":")):
            var a = p[unsafe_offset=s]
            var b = p[unsafe_offset=s + 1]
            var c = p[unsafe_offset=s + 3]
            var keep = (
                (a == UInt8(ord("r")) and b == UInt8(ord("q")) and c == UInt8(ord("i")))
                or (a == UInt8(ord("p")) and b == UInt8(ord("u")) and c == UInt8(ord("Z")))
                or (a == UInt8(ord("r")) and b == UInt8(ord("s")) and c == UInt8(ord("Z")))
            )
            if keep:
                out.append(UInt8(ord("\t")))
                out.extend(Span(unsafe_ptr=p.unsafe_offset(s), length=length))
        field += 1
        s = i + 1


@always_inline
def _store(p: MutPointer[UInt8, _], v: Int, bytes: Int):
    for k in range(bytes):
        p[unsafe_offset=k] = UInt8((v >> (8 * k)) & 0xFF)


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
    """A slice's own state: its scratch, what it placed this step and the last (being gathered into buckets), and
    its counts over the whole GAM."""

    var alignment: Alignment
    var walk: List[Int]
    var anchors: List[Anchor]
    var recs: List[UInt8]
    var prev: List[UInt8]
    var counts: List[Int]
    var span: List[Int]

    def __init__(out self, contigs: Int):
        self.alignment = Alignment()
        self.walk = List[Int]()
        self.anchors = List[Anchor]()
        self.recs = List[UInt8]()
        self.prev = List[UInt8]()
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
    var gatherers: Int
    """Side tasks gathering the last step's lines into buckets, each its share of the buckets."""
    var gathering: Bool
    """Whether the last step placed lines, for the next to gather."""
    var lines: Bool
    """Whether the messages are a GAF's lines (else a GAM's alignments)."""

    def __init__(out self, var map: RefMap, threads: Int, lines: Bool = False) raises:
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
        self.gatherers = threads
        self.gathering = False
        self.lines = lines

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
        if self.lines:
            self._take_lines(st[], data, messages, first, last)
            return
        for m in range(first, last):
            st[].alignment.decode[slim=True](data, messages[unsafe_offset=2 * m], messages[unsafe_offset=2 * m + 1])
            # The line goes straight after room for its record header, kept if it is placed.
            var at = len(st[].recs)
            st[].recs.resize(unsafe_uninit_length=at + HEADER)
            var path = slim_line(data, st[].alignment, self.map, st[].recs, st[].walk)
            var n = len(st[].recs) - at - HEADER - 1
            var why = 0
            var aligned = 0
            var kind = self._place(
                data,
                st[].alignment,
                st[].recs.unsafe_ptr().as_imm().unsafe_offset(at + HEADER),
                n,
                path,
                st[],
                placed,
                why,
                aligned,
            )
            if kind == _PLACED:
                var r = self.rank[placed.contig]
                st[].counts[C_LINES] += 1
                st[].counts[C_PLACED] += 1
                st[].span[r] = max(st[].span[r], aligned)
                var h = st[].recs.unsafe_ptr().unsafe_offset(at)
                _store(h, r, 4)
                _store(h.unsafe_offset(4), placed.start, 8)
                _store(h.unsafe_offset(12), placed.end, 8)
                _store(h.unsafe_offset(20), n, 4)
                st[].recs.resize(unsafe_uninit_length=at + HEADER + n)
                continue
            st[].recs.resize(unsafe_uninit_length=at)
            if kind == _SKIP:
                continue
            st[].counts[C_LINES] += 1
            if kind == _UNMAPPED:
                st[].counts[C_UNMAPPED] += 1
            elif kind == _UNPLACED:
                st[].counts[C_UNPLACED] += 1
                st[].counts[C_WHY + why - NO_ANCHOR] += 1
            else:
                st[].counts[C_MALFORMED] += 1

    def _take_lines(
        self,
        mut st: _Slice,
        data: ImmPointer[UInt8, ImmUntrackedOrigin],
        messages: ImmPointer[Int, ImmUntrackedOrigin],
        first: Int,
        last: Int,
    ):
        """Lines `first ..< last` of a GAF, read as pgcat reads them: a line's trailing carriage return left out, the
        line classified as text, and a placed one kept slim (its 12 columns and SLIM_TAGS)."""
        var placed = Placed(0, 0, 0)
        for m in range(first, last):
            var start = messages[unsafe_offset=2 * m]
            var n = messages[unsafe_offset=2 * m + 1] - start
            var p = data.unsafe_offset(start)
            if n > 0 and p[unsafe_offset=n - 1] == UInt8(ord("\r")):
                n -= 1
            var why = 0
            var aligned = 0
            var kind = self._classify(p, n, st.walk, st.anchors, placed, why, aligned)
            if kind == _SKIP:
                continue
            st.counts[C_LINES] += 1
            if kind == _PLACED:
                var r = self.rank[placed.contig]
                st.counts[C_PLACED] += 1
                st.span[r] = max(st.span[r], aligned)
                var at = len(st.recs)
                st.recs.resize(unsafe_uninit_length=at + HEADER)
                _append_slim(st.recs, p, n)
                var h = st.recs.unsafe_ptr().unsafe_offset(at)
                _store(h, r, 4)
                _store(h.unsafe_offset(4), placed.start, 8)
                _store(h.unsafe_offset(12), placed.end, 8)
                _store(h.unsafe_offset(20), len(st.recs) - at - HEADER, 4)
            elif kind == _UNMAPPED:
                st.counts[C_UNMAPPED] += 1
            elif kind == _UNPLACED:
                st.counts[C_UNPLACED] += 1
                st.counts[C_WHY + why - NO_ANCHOR] += 1
            else:
                st.counts[C_MALFORMED] += 1

    def _place(
        self,
        data: ImmPointer[UInt8, _],
        a: Alignment,
        p: ImmPointer[UInt8, _],
        n: Int,
        path: Tuple[Int, Int],
        mut st: _Slice,
        mut placed: Placed,
        mut why: Int,
        mut aligned: Int,
    ) -> Int:
        """pgcat's classify of the GAF line `p[:n]` (no newline) slim_line wrote for alignment `a`, from what it
        gave (the walk, the path start and end) rather than by reading the line again. A name with a tab in it moves
        the line's columns: that line is read as text, as pgcat reads it."""
        for k in range(a.name_length):
            if data[unsafe_offset=a.name_at + k] == UInt8(ord("\t")):
                return self._classify(p, n, st.walk, st.anchors, placed, why, aligned)
        if a.name_length > 0:
            var c = data[unsafe_offset=a.name_at]
            if c == UInt8(ord("#")) or c == UInt8(ord("@")):
                return _SKIP
        if not a.has_path or len(a.maps) == 0:
            return _UNMAPPED
        var ts = path[0]
        var te = path[1]
        # A negative number is no number to pgcat.
        if te < 0:
            return _MALFORMED
        for step in st.walk:
            if step < 0:
                return _MALFORMED
        why = self.map.project(st.walk, ts, te, st.anchors, placed)
        if why != 0:
            return _UNPLACED
        aligned = te - ts
        return _PLACED

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
        return self.gatherers if self.gathering else 0

    def side_task(mut self, index: Int) raises:
        """Gathers the last step's lines of this task's buckets (bucket b of rank r: (r << 12 | b) mod gatherers),
        compressing a bucket's lines as a frame once FRAME bytes of them gather."""
        for s in range(len(self.states)):
            var recs = imm_untracked(self.states[s].prev.unsafe_ptr())
            var total = len(self.states[s].prev)
            var at = 0
            while at < total:
                var r = _get32(recs, at)
                var b = _get64(recs, at + 4) >> BUCKET_SHIFT
                var length = _get32(recs, at + 20)
                if (r << 12 | b) % self.gatherers == index:
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

    def after_step(mut self, slices: Int) raises:
        """The step's lines become the next step's to gather; the buckets they need are made here, on one thread."""
        self.gathering = False
        for s in range(len(self.states)):
            self.states[s].prev.clear()
            if s < slices:
                swap(self.states[s].recs, self.states[s].prev)
            else:
                self.states[s].recs.clear()
            var recs = imm_untracked(self.states[s].prev.unsafe_ptr())
            var total = len(self.states[s].prev)
            self.gathering = self.gathering or total > 0
            var at = 0
            while at < total:
                var r = _get32(recs, at)
                var b = _get64(recs, at + 4) >> BUCKET_SHIFT
                if len(self.buckets[r]) <= b:
                    self.buckets[r].resize(b + 1, _Bucket())
                at += HEADER + _get32(recs, at + 20)

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
        # Windows of buckets, three deep: each step inflates and sorts window j, lays out window j - 1, deflates the
        # whole blocks of window j - 2 and writes window j - 3, all at once. Two slot sets, text buffers and output
        # sets take turns.
        var window = 4 * threads
        var windows = (len(order) + window - 1) // window
        var slots = List[List[_Slot]]()
        for _ in range(2):
            var group = List[_Slot]()
            for _ in range(window):
                group.append(_Slot())
            slots.append(group^)
        # Each window's first stream offset (its lines'), then the end of the last; each window's buckets; how many
        # outputs each output set holds.
        var win_start = List[Int](length=windows + 1, fill=0)
        var win_count = List[Int](length=windows, fill=0)
        var deflated = List[Int](length=2, fill=0)
        var bad = List[Int](length=1, fill=0)
        var end = len(w.texts[0])
        var step_ns = 0
        var serial_ns = 0
        var t_loop = Int(perf_counter_ns())
        for j in range(windows + 3):
            var s0 = Int(perf_counter_ns())
            # Window j's buckets, to inflate and sort.
            var gathers = 0
            if j < windows:
                gathers = min(window, len(order) - j * window)
                win_count[j] = gathers
                for k in range(gathers):
                    slots[j % 2][k].rank = order[j * window + k][0]
                    slots[j % 2][k].bucket = order[j * window + k][1]
            # Window j - 1, sorted, laid out after what window j - 2 leaves of its last block.
            var lay = j - 1
            var layouts = 0
            if lay >= 0 and lay < windows:
                layouts = win_count[lay]
                var x = lay % 2
                var start = end
                var base = (start // BLOCK) * BLOCK
                var at = start - base
                for k in range(layouts):
                    slots[x][k].at = at
                    at += slots[x][k].text_len
                w.texts[x].resize(unsafe_uninit_length=at)
                if lay > 0 and start > base:
                    var y = 1 - x
                    unsafe_memcpy(
                        dest=mut_untracked(w.texts[x].unsafe_ptr()),
                        src=imm_untracked(w.texts[y].unsafe_ptr()).unsafe_offset(base - w.text_base[y]),
                        count=start - base,
                    )
                w.text_base[x] = base
                win_start[lay] = start
                end = base + at
                win_start[lay + 1] = end
            # Window j - 2's whole blocks, to deflate.
            var dw = j - 2
            var deflates = 0
            var deflate_len = 0
            if dw >= 0 and dw < windows:
                var x = dw % 2
                deflate_len = (win_start[dw + 1] // BLOCK) * BLOCK - w.text_base[x]
                var blocks = deflate_len // BLOCK
                deflates = (blocks + BLOCKS_PER_TASK - 1) // BLOCKS_PER_TASK
                while len(w.outs[x]) < deflates:
                    w.outs[x].append(List[UInt8]())
                deflated[x] = deflates
            # Window j - 3's blocks, to write.
            var wr = j - 3
            var writes = 1 if wr >= 0 and wr < windows else 0
            var task = _FinishTask(
                mut_untracked(Pointer(to=self)),
                Pointer(to=lib).as_imm(),
                level,
                writes,
                mut_untracked(Pointer(to=w.file)),
                mut_untracked(w.outs[(j + 1) % 2].unsafe_ptr()),
                deflated[(j + 1) % 2] if writes == 1 else 0,
                gathers,
                mut_untracked(slots[j % 2].unsafe_ptr()),
                layouts,
                mut_untracked(slots[(j + 1) % 2].unsafe_ptr()),
                mut_untracked(w.texts[(j + 1) % 2].unsafe_ptr()),
                deflates,
                imm_untracked(w.texts[j % 2].unsafe_ptr()),
                deflate_len,
                mut_untracked(w.outs[j % 2].unsafe_ptr()),
                mut_untracked(bad.unsafe_ptr()),
            )
            var s1 = Int(perf_counter_ns())
            parallel_for(task, writes + gathers + layouts + deflates, threads)
            var s2 = Int(perf_counter_ns())
            if bad[0] == 1:
                raise Error("a bucket's lines did not inflate")
            if bad[0] == 2:
                raise Error("BGZF: a block did not deflate")
            # Window j - 1's lines wait for their blocks; window j - 2's blocks have their addresses, and the lines
            # ending in them are indexed.
            if layouts > 0:
                var x = lay % 2
                for k in range(layouts):
                    var contig = self.by_rank[slots[x][k].rank]
                    var stream_end = w.text_base[x] + slots[x][k].at
                    for rec in slots[x][k].recs:
                        stream_end += rec.out_len
                        w.pending.append(_Tab(contig, rec.start, rec.end, stream_end))
                    slots[x][k].raw = List[UInt8]()
                    slots[x][k].recs = List[_Rec]()
            if deflates > 0:
                w.addressed(dw % 2, deflates, w.text_base[dw % 2] + deflate_len)
                w.index_written(self.map.contigs)
            step_ns += s2 - s1
            serial_ns += (s1 - s0) + (Int(perf_counter_ns()) - s2)
        # The last window's partial block (or the header alone), then the end.
        var last = (windows - 1) % 2 if windows > 0 else 0
        var t_tail = Int(perf_counter_ns())
        w.finish_tail(lib, last, end)
        w.index_written(self.map.contigs)
        var t_close = Int(perf_counter_ns())
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
                    t"dromos timing (s): before {Float64(t_loop - t0) * 1e-9}, sort and write steps"
                    t" {Float64(step_ns) * 1e-9}, between {Float64(serial_ns) * 1e-9}, tail"
                    t" {Float64(t_close - t_tail) * 1e-9}, close {Float64(t1 - t_close) * 1e-9}"
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
    end has its address."""

    var path: String
    var tmp: String
    var file: FileHandle
    var level: Int
    var texts: List[List[UInt8]]
    """Two windows' text, taking turns: the one being laid out, the one being deflated."""
    var text_base: List[Int]
    """The stream offset of each text's first byte, a block boundary."""
    var outs: List[List[List[UInt8]]]
    """Two windows' deflated blocks, taking turns: the ones being deflated, the ones being written."""
    var base: Int
    """The stream bytes in blocks with addresses."""
    var written: Int
    """The compressed bytes of those blocks: the address of the next block."""
    var addrs: List[Int]
    """Every block's address, by block number."""
    var pending: List[_Tab]
    var head: Int
    var index: TbiBuilder

    def __init__(out self, path: String, level: Int) raises:
        self.path = path
        self.tmp = path + ".tmp"
        self.file = open(self.tmp, "w")
        self.level = level
        self.texts = List[List[UInt8]](length=2, fill=List[UInt8]())
        self.text_base = List[Int](length=2, fill=0)
        self.outs = List[List[List[UInt8]]](length=2, fill=List[List[UInt8]]())
        self.base = 0
        self.written = 0
        self.addrs = List[Int]()
        self.pending = List[_Tab]()
        self.head = 0
        self.index = TbiBuilder()

    def header(mut self, text: String) raises:
        """Header lines (each starting with '#', with its newline), before any line: tabix skips them, and the index
        starts after them."""
        self.texts[0].extend(text.as_bytes())
        if len(self.texts[0]) >= BLOCK:
            raise Error("header longer than a BGZF block")
        # Block 0 starts at address 0: the first line starts at (0, length).
        self.index.set_first_offset(len(self.texts[0]))

    def addressed(mut self, x: Int, outputs: Int, through: Int):
        """The blocks in output set `x`'s first `outputs` outputs, deflated from the stream up to `through`, get
        their addresses."""
        for t in range(outputs):
            var p = self.outs[x][t].unsafe_ptr().as_imm()
            var at = 0
            while at < len(self.outs[x][t]):
                self.addrs.append(self.written)
                var size = (Int(p[unsafe_offset=at + 16]) | Int(p[unsafe_offset=at + 17]) << 8) + 1
                self.written += size
                at += size
        self.base = through

    def finish_tail(mut self, lib: OwnedDLHandle, x: Int, end: Int) raises:
        """Deflates and writes what text `x` holds past the last whole block, up to stream offset `end`."""
        var rest = end - self.base
        if rest <= 0:
            return
        self.outs[x].clear()
        self.outs[x].append(List[UInt8]())
        var d = Deflater(lib, self.level)
        d.block(lib, self.texts[x].unsafe_ptr().unsafe_offset(self.base - self.text_base[x]), rest, self.outs[x][0])
        d.free(lib)
        self.file.write_bytes(Span(self.outs[x][0]))
        self.addressed(x, 1, end)

    def _voffset(self, s: Int) -> Int:
        """The virtual offset of stream offset `s`, in a block with an address or at the start of the next: a
        block's end reads as the start of the next block, as htslib reports it after reading a line."""
        if s == self.base:
            return self.written << 16
        return self.addrs[s // BLOCK] << 16 | (s % BLOCK)

    def index_written(mut self, contigs: List[String]) raises:
        """Indexes the lines whose ends lie in blocks with addresses."""
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
struct _FinishTask[lib_origin: ImmOrigin](ParallelTask):
    """One step of the final phase: `writes` (0 or 1) writes a window's deflated blocks in order, `gathers` inflate
    and sort a window's buckets, `layouts` lay out the window before it, and `deflates` deflate the whole blocks of
    the window before that, BLOCKS_PER_TASK blocks each."""

    var sink: MutPointer[IndexSink, MutUntrackedOrigin]
    var lib: ImmPointer[OwnedDLHandle, Self.lib_origin]
    var level: Int
    var writes: Int
    var file: MutPointer[FileHandle, MutUntrackedOrigin]
    var write_outs: MutPointer[List[UInt8], MutUntrackedOrigin]
    var write_count: Int
    var gathers: Int
    var gather_slots: MutPointer[_Slot, MutUntrackedOrigin]
    var layouts: Int
    var layout_slots: MutPointer[_Slot, MutUntrackedOrigin]
    var layout_dst: MutPointer[UInt8, MutUntrackedOrigin]
    var deflates: Int
    var deflate_src: ImmPointer[UInt8, ImmUntrackedOrigin]
    var deflate_len: Int
    var deflate_outs: MutPointer[List[UInt8], MutUntrackedOrigin]
    var bad: MutPointer[Int, MutUntrackedOrigin]

    def run(self, index: Int):
        var i = index
        if i < self.writes:
            for t in range(self.write_count):
                self.file[].write_bytes(Span(self.write_outs[unsafe_offset=t]))
            return
        i -= self.writes
        if i < self.gathers:
            try:
                self.sink[].gather(self.gather_slots[unsafe_offset=i])
            except:
                self.bad[] = 1
            return
        i -= self.gathers
        if i < self.layouts:
            self.sink[].lay_out(self.layout_slots[unsafe_offset=i], self.layout_dst)
            return
        i -= self.layouts
        ref out = self.deflate_outs[unsafe_offset=i]
        out.clear()
        try:
            var d = Deflater(self.lib[], self.level)
            var at = i * BLOCKS_PER_TASK * BLOCK
            var stop = min(self.deflate_len, at + BLOCKS_PER_TASK * BLOCK)
            while at < stop:
                var n = min(BLOCK, stop - at)
                d.block(self.lib[], self.deflate_src.unsafe_offset(at), n, out)
                at += n
            d.free(self.lib[])
        except:
            self.bad[] = 2
