"""`dromos surject`: reads' graph alignments as alignments to one reference contig (SAM), the records `pgcatalog
surject` writes for them, in the input's order (pgcatalog's surject.rs `surject_line` and `run`, ported line for
line; linearize.mojo ports the alignment it makes of each read).

Each read is rebuilt from its GAF line: the walk's node letters and the `cs` string (qualities from `bq`). A read
base aligned to a node the contig visits once takes that node's position; the read's orientation is that of most such
bases, and its placed stretches (blocks) are linearized. Each end is then clipped where the read diverges from its
own graph alignment. Reads with no base placed on the contig are left out.

A GAM's alignments are first written as the GAF lines `vg convert -G` writes (gam.mojo's `gaf_line`), so a GAM and
its GAF give the same records; an alignment none of whose nodes the contig visits once has no record and is passed
over from its slim decode."""

from std.atomic import Atomic
from std.ffi import external_call
from std.memory import Allocation, Layout, alloc, dealloc
from std.sys import stderr

from dromos.bam import RecordWriter
from dromos.bam_sort import BamStaging
from dromos.contig import Contig, Reference
from dromos.gaf_text import parse_unsigned
from dromos.gam import Alignment, gaf_line, put_bytes, put_int, put_text
from dromos.linearize import (
    CLIP_MATCH,
    CLIP_MISMATCH,
    S,
    Block,
    CigarOp,
    Window,
    clip_cigar,
    clip_gap,
    clip_lengths,
    linearize,
    upper,
)
from dromos.nodes import NodeTable, comp
from dromos.output import write_all
from dromos.parallel import Background, ParallelTask, parallel_for
from dromos.stream import MessageSink, imm_untracked, mut_untracked


@fieldwise_init
struct _Run(TrivialRegisterPassable):
    """Read letters placed on the contig: their offset in the read (walk orientation), how many, the reference
    position of the first, and +1 or -1 as the next goes along the reference."""

    var q: Int
    var n: Int
    var r: Int
    var sense: Int


struct _Scratch(Movable):
    """An alignment worker's lists, kept between reads."""

    var fields: List[Int]
    var walk: List[Int]
    var lens: List[Int]
    var read: List[UInt8]
    var score: List[Int]
    var runs: List[_Run]
    var blocks: List[Block]
    var merged: List[Block]
    var cigar: List[CigarOp]
    var hits: List[Int]

    def __init__(out self):
        self.fields = List[Int]()
        self.walk = List[Int]()
        self.lens = List[Int]()
        self.read = List[UInt8]()
        self.score = List[Int]()
        self.runs = List[_Run]()
        self.blocks = List[Block]()
        self.merged = List[Block]()
        self.cigar = List[CigarOp]()
        self.hits = List[Int]()


comptime _TAB = UInt8(ord("\t"))
comptime _NEWLINE = UInt8(ord("\n"))
comptime _COLON = UInt8(ord(":"))
comptime _STAR = UInt8(ord("*"))
comptime _PLUS = UInt8(ord("+"))
comptime _MINUS = UInt8(ord("-"))
comptime _FORWARD = UInt8(ord(">"))
comptime _REVERSE = UInt8(ord("<"))
comptime _HASH = UInt8(ord("#"))
comptime _QUAL = UInt8(ord("I"))
comptime _UPPER_A = UInt8(ord("A"))
comptime _UPPER_Z = UInt8(ord("Z"))
comptime _ZERO = UInt8(ord("0"))
comptime _NINE = UInt8(ord("9"))


@always_inline
def _alpha(c: UInt8) -> Bool:
    var u = upper(c)
    return u >= _UPPER_A and u <= _UPPER_Z


@always_inline
def _digit(c: UInt8) -> Bool:
    return c >= _ZERO and c <= _NINE


def _is(p: ImmPointer[UInt8, _], fields: List[Int], f: Int, text: StaticString) -> Bool:
    """Whether field `f` is `text`."""
    var at = fields[2 * f]
    var t = text.as_bytes()
    if fields[2 * f + 1] - at != len(t):
        return False
    for k in range(len(t)):
        if p[unsafe_offset=at + k] != t[k]:
            return False
    return True


@always_inline
def _number(p: ImmPointer[UInt8, _], fields: List[Int], f: Int, limit: Int = Int.MAX) -> Int:
    """Field `f` as an unsigned number (`parse_unsigned`), or -1."""
    return parse_unsigned(p, fields[2 * f], fields[2 * f + 1], limit)


def _tag(p: ImmPointer[UInt8, _], fields: List[Int], key: StaticString) -> Int:
    """The first field starting with `key` (cs:Z:, bq:Z:), or -1."""
    var t = key.as_bytes()
    for f in range(len(fields) // 2):
        var at = fields[2 * f]
        if fields[2 * f + 1] - at < len(t):
            continue
        var hit = True
        for k in range(len(t)):
            hit = hit and p[unsafe_offset=at + k] == t[k]
        if hit:
            return f
    return -1


def _parse_walk(p: ImmPointer[UInt8, _], fields: List[Int], f: Int, mut walk: List[Int]) -> Bool:
    """Field `f`, a GAF walk, as steps `node << 1 | reverse`; false unless it is all > or < each followed by a node
    id."""
    walk.clear()
    var i = fields[2 * f]
    var end = fields[2 * f + 1]
    while i < end:
        var c = p[unsafe_offset=i]
        if c != _FORWARD and c != _REVERSE:
            return False
        i += 1
        var st = i
        while i < end and _digit(p[unsafe_offset=i]):
            i += 1
        var node = parse_unsigned(p, st, i)
        if node < 0:
            return False
        walk.append(node << 1 | Int(c == _REVERSE))
    return True


def _split(p: ImmPointer[UInt8, _], n: Int, mut fields: List[Int]):
    """The tab-separated fields of `p[:n]`, each one's start and end."""
    fields.clear()
    var at = 0
    while True:
        var hit = external_call["memchr", Int](p.unsafe_offset(at), Int32(ord("\t")), n - at) if at < n else 0
        var e = n if hit == 0 else hit - Int(p)
        fields.append(at)
        fields.append(e)
        if e >= n:
            return
        at = e + 1


@always_inline
def _first_column(p: ImmPointer[UInt8, _], fields: List[Int]) -> Int:
    """The GAF's first column among the fields: 0, or 3 when column 5 is not a strand (a line with three columns
    before it, as pgcatalog reads them)."""
    var nf = len(fields) // 2
    if nf > 4 and (_is(p, fields, 4, "+") or _is(p, fields, 4, "-")):
        return 0
    return 3


@always_inline
def _place(contig: Contig, step: Int, length: Int, off: Int) -> Tuple[Int, Int]:
    """The reference position of letter `off` of walk step `step` (of `length` letters) and its direction (+1 or -1
    along the reference), or direction 0 if the contig does not visit its node once."""
    var pl = contig.placed(step >> 1)
    if pl < 0:
        return (0, 0)
    var rev = step & 1 == 1
    var node_rev = pl & 1 == 1
    var start = pl >> 1
    var fo = length - 1 - off if rev else off
    var r = start + (length - 1 - fo) if node_rev else start + fo
    return (r, 1 if rev == node_rev else -1)


struct Surjector(Movable):
    """What `surject_line` needs: the contigs, the node letters, and the options; each contig's name in the SAM and
    its catalog's windows."""

    var reference: Reference
    var nodes: NodeTable
    var sam_names: List[String]
    var graph: Bool
    var windows: List[List[Window]]
    var records: RecordWriter
    var complements: List[UInt8]
    """`comp` of every byte."""

    def __init__(
        out self,
        var reference: Reference,
        var nodes: NodeTable,
        var sam_names: List[String],
        graph: Bool,
        var windows: List[List[Window]],
    ):
        self.reference = reference^
        self.nodes = nodes^
        self.sam_names = sam_names^
        self.graph = graph
        self.windows = windows^
        self.records = RecordWriter()
        self.complements = List[UInt8](capacity=256)
        for b in range(256):
            self.complements.append(comp(UInt8(b)))

    def touches(self, p: ImmPointer[UInt8, _], n: Int, mut fields: List[Int], mut walk: List[Int]) -> Bool:
        """Whether GAF line `p[:n]` (no newline) may place a base on a contig: its walk holds a node a contig
        visits once (`surject_line` places none otherwise)."""
        _split(p, n, fields)
        var k = _first_column(p, fields)
        if len(fields) // 2 < k + 12 or not _parse_walk(p, fields, k + 5, walk):
            return False
        for step in walk:
            if self.reference.placed_anywhere(step >> 1):
                return True
        return False

    def surject_read[
        bam: Bool
    ](self, p: ImmPointer[UInt8, _], n: Int, mut s: _Scratch, mut out: List[UInt8]) raises -> Int:
        """The records of GAF line `p[:n]` (no newline) appended to `out`, one on each contig it places a base on
        (in contig order), each as `surject_line` makes it; how many."""
        if len(self.reference.contigs) == 1:
            return Int(self.surject_line[bam](p, n, 0, s, out))
        _split(p, n, s.fields)
        var k = _first_column(p, s.fields)
        if len(s.fields) // 2 < k + 12 or not _parse_walk(p, s.fields, k + 5, s.walk):
            return 0
        s.hits.clear()
        for step in s.walk:
            self.reference.add_contigs_of(step >> 1, s.hits)
        var made = 0
        for i in range(len(s.hits)):
            var c = s.hits[i]
            made += Int(self.surject_line[bam](p, n, c, s, out, parsed=True))
        return made

    def surject_line[
        bam: Bool
    ](
        self, p: ImmPointer[UInt8, _], n: Int, c: Int, mut s: _Scratch, mut out: List[UInt8], parsed: Bool = False
    ) raises -> Bool:
        """The record on contig `c` of GAF line `p[:n]` (no newline) appended to `out`, as SAM text (with its
        newline) or as BAM; false, nothing appended, if it places no base on the contig. `parsed`: the line's fields
        and walk are in `s` already (`surject_read` found its contigs from them)."""
        ref contig = self.reference.contigs[c]
        if not parsed:
            _split(p, n, s.fields)
        var k = _first_column(p, s.fields)
        if len(s.fields) // 2 < k + 12 or _is(p, s.fields, k + 5, "*"):
            return False
        var minus = _is(p, s.fields, k + 4, "-")
        var qlen = _number(p, s.fields, k + 1)
        var qstart = _number(p, s.fields, k + 2)
        var qend = _number(p, s.fields, k + 3)
        var pstart = _number(p, s.fields, k + 7)
        if qlen < 0 or qstart < 0 or qend < 0 or pstart < 0:
            return False
        var mapq = _number(p, s.fields, k + 11, 0xFFFF_FFFF)
        if mapq < 0:
            mapq = 255
        var cs_field = _tag(p, s.fields, "cs:Z:")
        if cs_field < 0:
            return False
        var bq_field = _tag(p, s.fields, "bq:Z:")
        if not parsed and not _parse_walk(p, s.fields, k + 5, s.walk):
            return False
        # No base can be placed unless a node of the walk is the contig's.
        var any_placed = False
        for step in s.walk:
            any_placed = any_placed or contig.placed(step >> 1) >= 0
        if not any_placed:
            return False
        var walk_len = len(s.walk)
        s.lens.clear()
        for step in s.walk:
            s.lens.append(self.nodes.length(step >> 1))
        var w = 0
        var off = pstart
        while w < walk_len and off >= s.lens[w]:
            off -= s.lens[w]
            w += 1
        # Read letters in walk orientation, placed runs, and each letter's clipping score against the walk (an
        # indel's on the letter before it).
        s.read.clear()
        s.runs.clear()
        s.score.clear()
        var i = s.fields[2 * cs_field] + 5
        var cs_end = s.fields[2 * cs_field + 1]
        while i < cs_end:
            var op = p[unsafe_offset=i]
            i += 1
            var st = i
            if op == _COLON:
                while i < cs_end and _digit(p[unsafe_offset=i]):
                    i += 1
                var left = parse_unsigned(p, st, i)
                if left < 0:
                    return False
                while left > 0:
                    if w >= walk_len:
                        return False
                    var step = s.walk[w]
                    var length = s.lens[w]
                    if length == 0:
                        # No letters: a node the graph lacks.
                        return False
                    var letters = self.nodes.start(step >> 1)
                    var take = min(left, length - off)
                    var at = len(s.read)
                    s.read.resize(unsafe_uninit_length=at + take)
                    var dst = s.read.unsafe_ptr().unsafe_offset(at)
                    if step & 1 == 1:
                        self.nodes.copy_reverse(dst, letters + length - off - 1, take)
                    else:
                        self.nodes.copy_forward(dst, letters + off, take)
                    s.score.resize(at + take, CLIP_MATCH)
                    var pl = _place(contig, step, length, off)
                    if pl[1] != 0:
                        s.runs.append(_Run(len(s.read) - take, take, pl[0], pl[1]))
                    off += take
                    left -= take
                    if off == length:
                        w += 1
                        off = 0
            elif op == _STAR:
                if i + 2 > cs_end or w >= walk_len:
                    return False
                s.read.append(upper(p[unsafe_offset=i + 1]))
                i += 2
                s.score.append(CLIP_MISMATCH)
                var pl = _place(contig, s.walk[w], s.lens[w], off)
                if pl[1] != 0:
                    s.runs.append(_Run(len(s.read) - 1, 1, pl[0], pl[1]))
                off += 1
                if off == s.lens[w]:
                    w += 1
                    off = 0
            elif op == _PLUS:
                while i < cs_end and _alpha(p[unsafe_offset=i]):
                    i += 1
                for t in range(st, i):
                    s.read.append(upper(p[unsafe_offset=t]))
                var count = i - st
                for _ in range(count):
                    s.score.append(0)
                if len(s.score) >= count + 1:
                    s.score[len(s.score) - count - 1] += clip_gap(count)
            elif op == _MINUS:
                while i < cs_end and _alpha(p[unsafe_offset=i]):
                    i += 1
                var left = i - st
                if len(s.score) > 0:
                    s.score[len(s.score) - 1] += clip_gap(left)
                while left > 0:
                    if w >= walk_len or s.lens[w] == 0:
                        return False
                    var take = min(left, s.lens[w] - off)
                    off += take
                    left -= take
                    if off == s.lens[w]:
                        w += 1
                        off = 0
            else:
                return False
        var fwd_bases = 0
        var rev_bases = 0
        for r in s.runs:
            if r.sense == 1:
                fwd_bases += r.n
            else:
                rev_bases += r.n
        if fwd_bases + rev_bases == 0:
            return False
        var fwd = fwd_bases >= rev_bases
        var rl = len(s.read)
        s.blocks.clear()
        for r in s.runs:
            if (r.sense == 1) == fwd:
                s.blocks.append(Block(r.q, r.r, r.n) if fwd else Block(rl - (r.q + r.n), r.r + 1 - r.n, r.n))

        def less(a: Block, b: Block) -> Bool:
            return a.less(b)

        sort(Span(s.blocks), less)
        s.merged.clear()
        for b in s.blocks:
            var last = len(s.merged) - 1
            if last >= 0 and s.merged[last].q + s.merged[last].n == b.q and s.merged[last].r + s.merged[last].n == b.r:
                s.merged[last].n += b.n
                continue
            s.merged.append(b)
        # The read in the BAM's orientation: a GAF read on `-` is the path's letters reverse-complemented.
        if not fwd:
            var r = s.read.unsafe_ptr()
            var lo = 0
            var hi = rl - 1
            while lo < hi:
                var a = r[unsafe_offset=lo]
                r[unsafe_offset=lo] = self.complements[Int(r[unsafe_offset=hi])]
                r[unsafe_offset=hi] = self.complements[Int(a)]
                lo += 1
                hi -= 1
            if lo == hi:
                r[unsafe_offset=lo] = self.complements[Int(r[unsafe_offset=lo])]
        var pos = 0
        if not linearize(Span(s.read), s.merged, Span(contig.seq), self.graph, self.windows[c], pos, s.cigar):
            return False
        # Clipping on the read's alignment to the walk, in the BAM's orientation.
        if not fwd:
            s.score.reverse()
        var clips = clip_lengths(s.score)
        var nc = len(s.cigar)
        var lead = s.cigar[0].length if nc > 0 and s.cigar[0].op == S else 0
        var trail = s.cigar[nc - 1].length if nc > 1 and s.cigar[nc - 1].op == S else 0
        var clipped = clip_cigar(pos, s.cigar, max(clips[0], lead), max(clips[1], trail))
        # bq covers the whole read; the aligned part (qstart..qend, on the read as sequenced) is what was rebuilt.
        # bq runs along the read as sequenced; the walk's letters along the read reversed when it is on `-`.
        var reverse_qual = minus == fwd
        var bq_at = s.fields[2 * bq_field] + 5 if bq_field >= 0 else 0
        var bq_len = s.fields[2 * bq_field + 1] - bq_at if bq_field >= 0 else -1
        var qual_at = -1
        if bq_len == qlen and qend <= qlen and qend >= qstart and qend - qstart == rl:
            qual_at = bq_at + qstart
        comptime if bam:
            self.records.put(
                out,
                p.unsafe_offset(s.fields[2 * k]),
                s.fields[2 * k + 1] - s.fields[2 * k],
                0 if fwd else 16,
                c,
                clipped[0],
                mapq,
                clipped[1],
                Span(s.read),
                p,
                qual_at,
                reverse_qual,
            )
            return True
        put_bytes(out, p, s.fields[2 * k], s.fields[2 * k + 1] - s.fields[2 * k])
        if fwd:
            put_text(out, "\t0\t")
        else:
            put_text(out, "\t16\t")
        out.extend(self.sam_names[c].as_bytes())
        out.append(_TAB)
        put_int(out, clipped[0] + 1)
        out.append(_TAB)
        put_int(out, mapq)
        out.append(_TAB)
        for o in clipped[1]:
            put_int(out, o.length)
            out.append(o.op)
        put_text(out, "\t*\t0\t0\t")
        out.extend(Span(s.read))
        out.append(_TAB)
        if qual_at < 0:
            for _ in range(rl):
                out.append(_QUAL)
        elif reverse_qual:
            for t in range(qual_at + rl - 1, qual_at - 1, -1):
                out.append(p[unsafe_offset=t])
        else:
            put_bytes(out, p, qual_at, rl)
        out.append(_NEWLINE)
        return True


comptime SLICES_PER_THREAD = 16
"""Slices of a chunk's messages: each copies its share of the reads to align into a piece."""
comptime IN_FLIGHT = 512 << 20
"""The lines of pieces published and not yet written out, at most (but for one step's more): reading waits for the
aligners past it."""
comptime RING_STEPS = 8
"""Steps' pieces the ring holds."""


struct _Slice(Movable):
    """What a slice of a chunk needs as it gathers reads, and its count of alignments read."""

    var alignment: Alignment
    var fields: List[Int]
    var walk: List[Int]
    var alignments: Int

    def __init__(out self):
        self.alignment = Alignment()
        self.fields = List[Int]()
        self.walk = List[Int]()
        self.alignments = 0


struct _Piece(Copyable, Movable):
    """A slice's reads of one chunk: the GAF lines of those that touch a contig (newline-ended), then their
    records."""

    var lines: List[UInt8]
    var text: List[UInt8]
    var bytes: Int
    """Its lines' bytes, as published."""
    var written: Int
    """Reads with a record."""
    var records: Int

    def __init__(out self):
        self.lines = List[UInt8]()
        self.text = List[UInt8]()
        self.bytes = 0
        self.written = 0
        self.records = 0


comptime _PUBLISHED = 0
"""Pieces published: piece g (a sequence number) lies in ring slot g mod the ring's size."""
comptime _NEXT = 1
"""The next piece an aligner takes."""
comptime _FINISHED = 2
"""1 once no more pieces will be published; 2 to stop at once (after an error)."""
comptime _FAILED = 3
comptime _CONTROLS = 4


@fieldwise_init
struct _Aligners[bam: Bool](Movable, ParallelTask):
    """The aligners, on threads of their own for the whole run: each takes the next piece published, makes its
    records and marks it done, until no more will be."""

    var sink: MutPointer[SurjectSink[Self.bam], MutUntrackedOrigin]

    def run(self, index: Int):
        self.sink[].align_loop(index)


@fieldwise_init
struct _Writers[bam: Bool](ParallelTask):
    """The writer side tasks, run on their own (when reading waits, and at the end)."""

    var sink: MutPointer[SurjectSink[Self.bam], MutUntrackedOrigin]
    var failed: MutPointer[Bool, MutUntrackedOrigin]

    def run(self, index: Int):
        try:
            self.sink[].side_task(index)
        except:
            self.failed[] = True


struct SurjectSink[bam: Bool](MessageSink, Movable):
    """Each slice gathers the GAF lines of its reads that touch a contig (a GAM's alignments written as GAF) into a
    piece of the ring; after the step the pieces are published, in order. Aligners on threads of their own take the
    pieces one at a time as they are published and make their records, so no thread waits at the end of a batch for
    the slowest read of it. In each step, the pieces done since the last (a prefix, in order) are written out: as
    SAM on stdout by one side task, or, with `bam`, handed to the staging buckets by `shares` side tasks. Reading
    waits when the pieces not yet written hold IN_FLIGHT bytes of lines."""

    var s: Surjector
    var states: List[_Slice]
    var workers: List[_Scratch]
    var ring: List[_Piece]
    var ring_ptr: MutPointer[_Piece, MutUntrackedOrigin]
    var done: Allocation[Atomic[Int64]]
    """Per ring slot: 1 + the sequence number of the last piece done there."""
    var done_ptr: MutPointer[Atomic[Int64], MutUntrackedOrigin]
    var control: Allocation[Atomic[Int64]]
    var control_ptr: MutPointer[Atomic[Int64], MutUntrackedOrigin]
    var aligners: List[Background[_Aligners[Self.bam]]]
    var threads: Int
    var produced: Int
    """The sequence number of the first piece the next step's slices fill."""
    var write_from: Int
    var write_to: Int
    """The pieces the next step writes out."""
    var in_flight: Int
    """Lines of the pieces published and not yet written."""
    var lines: Bool
    var written: Int
    """Reads with a record, over the pieces written."""
    var records: Int
    var failed: Bool
    var staging: BamStaging
    var shares: Int

    def __init__(out self, var s: Surjector, threads: Int, header: String, var staging: BamStaging):
        """`header` is written on stdout first (SAM); BAM records go to `staging`."""
        self.s = s^
        self.states = List[_Slice]()
        for _ in range(SLICES_PER_THREAD * threads):
            self.states.append(_Slice())
        self.workers = List[_Scratch]()
        for _ in range(threads):
            self.workers.append(_Scratch())
        var slots = RING_STEPS * len(self.states)
        self.ring = List[_Piece](length=slots, fill=_Piece())
        self.ring_ptr = mut_untracked(self.ring.unsafe_ptr())
        self.done = alloc(Layout[Atomic[Int64]](length=slots))
        self.done_ptr = mut_untracked(self.done.unsafe_ptr())
        for i in range(slots):
            self.done_ptr.unsafe_offset(i).unsafe_write(Atomic[Int64](0))
        self.control = alloc(Layout[Atomic[Int64]](length=_CONTROLS))
        self.control_ptr = mut_untracked(self.control.unsafe_ptr())
        for i in range(_CONTROLS):
            self.control_ptr.unsafe_offset(i).unsafe_write(Atomic[Int64](0))
        self.aligners = List[Background[_Aligners[Self.bam]]]()
        self.threads = threads
        self.produced = 0
        self.write_from = 0
        self.write_to = 0
        self.in_flight = 0
        self.lines = False
        self.written = 0
        self.records = 0
        self.failed = False
        self.shares = staging.shares
        self.staging = staging^
        comptime if not Self.bam:
            write_all(1, header.unsafe_ptr(), header.byte_length())

    def __deinit__(deinit self):
        # The aligners are joined by `finish`; any still running are stopped and joined before the ring goes.
        self.control_ptr.unsafe_offset(_FINISHED)[].store(2)
        for k in range(len(self.aligners)):
            self.aligners[k].join()
        dealloc(self.done^)
        dealloc(self.control^)

    @always_inline
    def _slot(self, g: Int) -> MutPointer[_Piece, MutUntrackedOrigin]:
        return self.ring_ptr.unsafe_offset(g % len(self.ring))

    def begin(mut self, lines: Bool) raises:
        self.lines = lines
        self.aligners.append(
            Background(_Aligners[Self.bam](mut_untracked(Pointer(to=self))), self.threads, self.threads)
        )

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
        var piece = self._slot(self.produced + slice)
        piece[].lines.clear()
        piece[].text.clear()
        piece[].written = 0
        piece[].records = 0
        for m in range(first, last):
            var at = messages[unsafe_offset=2 * m]
            var stop = messages[unsafe_offset=2 * m + 1]
            if self.lines:
                # A piece of whole lines: each read as Rust's lines() gives it (no newline, nor a carriage return
                # before it), empty lines and comments passed over.
                while at < stop:
                    var hit = external_call["memchr", Int](data.unsafe_offset(at), Int32(ord("\n")), stop - at)
                    var end = stop if hit == 0 else hit - Int(data)
                    var n = end - at
                    if n > 0 and data[unsafe_offset=at + n - 1] == UInt8(ord("\r")):
                        n -= 1
                    if n > 0 and data[unsafe_offset=at] != _HASH:
                        st[].alignments += 1
                        if self.s.touches(data.unsafe_offset(at), n, st[].fields, st[].walk):
                            put_bytes(piece[].lines, data, at, n)
                            piece[].lines.append(_NEWLINE)
                    at = end + 1
                continue
            ref a = st[].alignment
            a.decode[slim=True](data, at, stop)
            if a.name_length > 0 and data[unsafe_offset=a.name_at] == _HASH:
                continue
            st[].alignments += 1
            # No base can be placed unless a node of the alignment is a contig's.
            var any_placed = False
            for map in a.maps:
                any_placed = any_placed or self.s.reference.placed_anywhere(map.node)
            if any_placed:
                a.decode(data, at, stop)
                gaf_line(data, a, self.s.nodes, piece[].lines)

    def align_loop(mut self, index: Int):
        """Aligner `index`: the next piece published, one at a time, until no more will be."""
        var scratch = mut_untracked(Pointer(to=self.workers[index]))
        ref published = self.control_ptr.unsafe_offset(_PUBLISHED)[]
        ref finished = self.control_ptr.unsafe_offset(_FINISHED)[]
        while True:
            var g = Int(self.control_ptr.unsafe_offset(_NEXT)[].fetch_add(1))
            while g >= Int(published.load()):
                if finished.load() != 0 and g >= Int(published.load()):
                    return
                _ = external_call["usleep", Int32](UInt32(200))
            if finished.load() == 2:
                return
            var piece = self._slot(g)
            try:
                var p = piece[].lines.unsafe_ptr().as_imm()
                var n = len(piece[].lines)
                var at = 0
                while at < n:
                    var hit = external_call["memchr", Int](p.unsafe_offset(at), Int32(ord("\n")), n - at)
                    var end = hit - Int(p)
                    var made = self.s.surject_read[Self.bam](p.unsafe_offset(at), end - at, scratch[], piece[].text)
                    if made > 0:
                        piece[].written += 1
                        piece[].records += made
                    at = end + 1
            except:
                self.control_ptr.unsafe_offset(_FAILED)[].store(1)
            piece[].lines = List[UInt8]()
            self.done_ptr.unsafe_offset(g % len(self.ring))[].store(Int64(g + 1))

    def _done_through(self, start: Int) -> Int:
        """The first piece from `start` on not yet done (or not yet published)."""
        var g = start
        var published = Int(self.control_ptr.unsafe_offset(_PUBLISHED)[].load())
        while g < published and Int(self.done_ptr.unsafe_offset(g % len(self.ring))[].load()) == g + 1:
            g += 1
        return g

    def _check(self) raises:
        if self.control_ptr.unsafe_offset(_FAILED)[].load() != 0:
            raise Error("a read could not be aligned")

    def side_tasks(self) -> Int:
        if self.write_to == self.write_from:
            return 0
        comptime if Self.bam:
            return self.shares
        return 1

    def side_task(mut self, index: Int) raises:
        comptime if Self.bam:
            # Share `index` of the staging buckets takes its records of every piece, in order.
            for g in range(self.write_from, self.write_to):
                var piece = self._slot(g)
                self.staging.take_share(index, piece[].text.unsafe_ptr(), len(piece[].text))
        else:
            for g in range(self.write_from, self.write_to):
                var piece = self._slot(g)
                write_all(1, piece[].text.unsafe_ptr(), len(piece[].text))

    def _wrote(mut self) raises:
        """Accounts for the pieces just written out (on one thread): counted, their records freed, their slots
        free; then the next pieces to write, those done since."""
        comptime if Self.bam:
            if self.write_to > self.write_from:
                self.staging.after_shares()
        for g in range(self.write_from, self.write_to):
            var piece = self._slot(g)
            self.written += piece[].written
            self.records += piece[].records
            self.in_flight -= piece[].bytes
            piece[].text = List[UInt8]()
        self.write_from = self.write_to
        self.write_to = self._done_through(self.write_from)
        comptime if Self.bam:
            for g in range(self.write_from, self.write_to):
                var piece = self._slot(g)
                self.staging.make_room(piece[].text.unsafe_ptr(), len(piece[].text))

    def _write_now(mut self) raises:
        """Writes out the pieces due, on `threads` threads, now (not beside a step)."""
        var task = _Writers[Self.bam](mut_untracked(Pointer(to=self)), mut_untracked(Pointer(to=self.failed)))
        parallel_for(task, self.side_tasks(), self.threads)
        if self.failed:
            raise Error("records could not be written")

    def after_step(mut self, slices: Int) raises:
        self._check()
        # The pieces the slices filled are published, in order.
        for k in range(slices):
            var piece = self._slot(self.produced + k)
            piece[].bytes = len(piece[].lines)
            self.in_flight += piece[].bytes
        self.produced += slices
        self.control_ptr.unsafe_offset(_PUBLISHED)[].store(Int64(self.produced))
        # The pieces written in the step just run, and those to write next.
        self._wrote()
        # Room for the next step's pieces: slots written out, and lines in flight under the bound. Done pieces are
        # written out here while waiting, so the wait always ends.
        while self.produced + len(self.states) > self.write_from + len(self.ring) or self.in_flight >= IN_FLIGHT:
            self._check()
            if self.write_to > self.write_from:
                self._write_now()
                self._wrote()
            else:
                _ = external_call["usleep", Int32](UInt32(500))
                self.write_to = self._done_through(self.write_from)
                comptime if Self.bam:
                    for g in range(self.write_from, self.write_to):
                        var piece = self._slot(g)
                        self.staging.make_room(piece[].text.unsafe_ptr(), len(piece[].text))

    def finish(mut self, threads: Int, complete: Bool) raises:
        """Once the input is `complete`, waits for the aligners and writes out what is left. After an error (not
        `complete`), the aligners stop, and SAM records made so far are written out."""
        if not complete:
            self.control_ptr.unsafe_offset(_FINISHED)[].store(2)
            for k in range(len(self.aligners)):
                self.aligners[k].join()
            self.aligners.clear()
            comptime if not Self.bam:
                if self.write_to > self.write_from:
                    self._write_now()
            return
        self.control_ptr.unsafe_offset(_FINISHED)[].store(1)
        # Pieces written in the step that ran last are accounted for by after_step; the rest now.
        while self.write_from < self.produced:
            self._check()
            if self.write_to > self.write_from:
                self._write_now()
                self._wrote()
            else:
                _ = external_call["usleep", Int32](UInt32(500))
                self.write_to = self._done_through(self.write_from)
                comptime if Self.bam:
                    for g in range(self.write_from, self.write_to):
                        var piece = self._slot(g)
                        self.staging.make_room(piece[].text.unsafe_ptr(), len(piece[].text))
        for k in range(len(self.aligners)):
            self.aligners[k].join()
        self.aligners.clear()
        self._check()

    def report(self, what: String):
        """The counts, on stderr: alignments read, and those with a record on `what` (the contig, or the contigs)."""
        var n = 0
        for t in range(len(self.states)):
            n += self.states[t].alignments
        var contigs = len(self.s.reference.contigs)
        if contigs == 1:
            print(
                (
                    t"[dromos] surject: {n} alignments, {self.written} written, {n - self.written} without a base on"
                    t" {what}"
                ),
                file=stderr,
            )
            return
        print(
            (
                t"[dromos] surject: {n} alignments, {self.written} with a base on {what} ({self.records} records),"
                t" {n - self.written} without"
            ),
            file=stderr,
        )
