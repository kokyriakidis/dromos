"""`dromos index`: a GAM made searchable by reference position, as pgphase reads it, in one pass and one process. It
writes the file `pgcat index-gaf --slim --ref-map` writes from the GAM's GAF (the same lines, blocks and index).

Each alignment is placed on the reference by the nodes of its walk that lie on the reference sample's paths
(refmap.mojo), written as `contig  start+1  end  <its 12 GAF columns>`, sorted by contig name, start, end and line,
bgzipped in blocks of exactly BLOCK bytes and tabix-indexed. Header lines `#max_span <contig> <bases>` give the most
graph bases any alignment on each contig spans.

While the input streams, each slice places its alignments (or classifies its GAF lines, gaf_text.mojo) and keeps
the placed ones' lines; the step after, side tasks hand them to the staging, by bucket (staging.mojo). At the end
the buckets are sorted and written, the index built as they go (sort_write.mojo).
"""

from std.ffi import external_call
from std.sys import stderr
from std.time import perf_counter_ns

from dromos.gaf_text import MALFORMED, PLACED, SKIP, UNMAPPED, UNPLACED, append_slim, classify
from dromos.gam import Alignment, slim_line
from dromos.refmap import NO_ANCHOR, Anchor, Placed, RefMap, reason_name
from dromos.sort_write import write_sorted
from dromos.staging import HEADER, Staging, put_record
from dromos.stream import MessageSink, imm_untracked, mut_untracked

# Why an alignment has no place, counted per slice: input lines, placed, unmapped, unplaced (then by reason, in
# refmap.mojo's order) and malformed.
comptime C_LINES = 0
comptime C_PLACED = 1
comptime C_UNMAPPED = 2
comptime C_UNPLACED = 3
comptime C_WHY = 4
comptime C_MALFORMED = 8
comptime N_COUNTS = 9


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


struct IndexSink(MessageSink, Movable):
    """GAM messages placed and gathered by bucket (see the module's notes)."""

    var map: RefMap
    var rank: List[Int]
    """Each contig's rank: contigs in name order."""
    var by_rank: List[Int]
    var states: List[_Slice]
    var staging: Staging
    var gatherers: Int
    """Side tasks handing the last step's lines over, each its share of the buckets."""
    var gathering: Bool
    """Whether the last step placed lines, for the next to gather."""
    var lines: Bool
    """Whether the messages are a GAF's lines (else a GAM's alignments)."""

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
        var names = List[String](capacity=n)
        for r in range(n):
            names.append(map.contigs[self.by_rank[r]])
        self.map = map^
        self.states = List[_Slice]()
        for _ in range(4 * threads):
            self.states.append(_Slice(n))
        self.staging = Staging(names^, threads)
        self.gatherers = threads
        self.gathering = False
        self.lines = False

    def begin(mut self, lines: Bool) raises:
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
            if kind == PLACED:
                var r = self.rank[placed.contig]
                st[].counts[C_LINES] += 1
                st[].counts[C_PLACED] += 1
                st[].span[r] = max(st[].span[r], aligned)
                st[].recs.resize(unsafe_uninit_length=at + HEADER + n)
                put_record(st[].recs, at, r, placed.start, placed.end)
                continue
            st[].recs.resize(unsafe_uninit_length=at)
            if kind == SKIP:
                continue
            st[].counts[C_LINES] += 1
            if kind == UNMAPPED:
                st[].counts[C_UNMAPPED] += 1
            elif kind == UNPLACED:
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
        """The lines of pieces `first ..< last` of a GAF (frame_lines), read as pgcat reads them: a line's trailing
        carriage return left out, the line classified as text, and a placed one kept slim (its 12 columns and
        SLIM_TAGS)."""
        var placed = Placed(0, 0, 0)
        for m in range(first, last):
            var at = messages[unsafe_offset=2 * m]
            var stop = messages[unsafe_offset=2 * m + 1]
            while at < stop:
                var hit = external_call["memchr", Int](data.unsafe_offset(at), Int32(ord("\n")), stop - at)
                var end = stop if hit == 0 else hit - Int(data)
                self._take_line(st, data.unsafe_offset(at), end - at, placed)
                at = end + 1

    @always_inline
    def _take_line(self, mut st: _Slice, p: ImmPointer[UInt8, _], length: Int, mut placed: Placed):
        """One GAF line, `p[:length]` without its newline."""
        var n = length
        if n > 0 and p[unsafe_offset=n - 1] == UInt8(ord("\r")):
            n -= 1
        var why = 0
        var aligned = 0
        var kind = classify(self.map, p, n, st.walk, st.anchors, placed, why, aligned)
        if kind == SKIP:
            return
        st.counts[C_LINES] += 1
        if kind == PLACED:
            var r = self.rank[placed.contig]
            st.counts[C_PLACED] += 1
            st.span[r] = max(st.span[r], aligned)
            var at = len(st.recs)
            st.recs.resize(unsafe_uninit_length=at + HEADER)
            append_slim(st.recs, p, n)
            put_record(st.recs, at, r, placed.start, placed.end)
        elif kind == UNMAPPED:
            st.counts[C_UNMAPPED] += 1
        elif kind == UNPLACED:
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
                return classify(self.map, p, n, st.walk, st.anchors, placed, why, aligned)
        if a.name_length > 0:
            var c = data[unsafe_offset=a.name_at]
            if c == UInt8(ord("#")) or c == UInt8(ord("@")):
                return SKIP
        if not a.has_path or len(a.maps) == 0:
            return UNMAPPED
        var ts = path[0]
        var te = path[1]
        # A negative number is no number to pgcat.
        if te < 0:
            return MALFORMED
        for step in st.walk:
            if step < 0:
                return MALFORMED
        why = self.map.project(st.walk, ts, te, st.anchors, placed)
        if why != 0:
            return UNPLACED
        aligned = te - ts
        return PLACED

    def side_tasks(self) -> Int:
        return self.gatherers if self.gathering else 0

    def side_task(mut self, index: Int) raises:
        """Hands the last step's lines in share `index`'s buckets over to the staging."""
        for s in range(len(self.states)):
            self.staging.take_share(index, imm_untracked(self.states[s].prev.unsafe_ptr()), len(self.states[s].prev))

    def after_step(mut self, slices: Int) raises:
        """The step's lines become the next step's to hand over; the buckets they need are made here, on one thread."""
        self.gathering = False
        for s in range(len(self.states)):
            self.states[s].prev.clear()
            if s < slices:
                swap(self.states[s].recs, self.states[s].prev)
            else:
                self.states[s].recs.clear()
            self.staging.make_room(imm_untracked(self.states[s].prev.unsafe_ptr()), len(self.states[s].prev))
            self.gathering = self.gathering or len(self.states[s].prev) > 0

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
        write_sorted(self.staging, max_span, out_path, level, threads)
        var t1 = Int(perf_counter_ns())
        print(
            (
                t"[dromos] {out_path} and {out_path}.tbi: {counts[C_PLACED]} alignments ({Float64(t1 - t0) * 1e-9} s"
                t" sorting and writing)"
            ),
            file=stderr,
        )


@always_inline
def _name_less(a: String, b: String) -> Bool:
    """`a` before `b` by their bytes (Rust's ordering of byte strings)."""
    var x = a.as_bytes()
    var y = b.as_bytes()
    for i in range(min(len(x), len(y))):
        if x[i] != y[i]:
            return x[i] < y[i]
    return len(x) < len(y)
