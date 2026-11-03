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
from std.sys import stderr

from dromos.contig import Contig, comp
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
from dromos.nodes import NodeTable
from dromos.output import write_all
from dromos.parallel import ParallelTask, parallel_for
from dromos.stream import MessageSink, mut_untracked


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


struct Surjector(Movable):
    """What `surject_line` needs: the contig, the node letters, and the options."""

    var contig: Contig
    var nodes: NodeTable
    var sam_contig: String
    var graph: Bool
    var windows: List[Window]

    def __init__(
        out self, var contig: Contig, var nodes: NodeTable, sam_contig: String, graph: Bool, var windows: List[Window]
    ):
        self.contig = contig^
        self.nodes = nodes^
        self.sam_contig = sam_contig
        self.graph = graph
        self.windows = windows^

    @always_inline
    def _place(self, step: Int, length: Int, off: Int) -> Tuple[Int, Int]:
        """The reference position of letter `off` of walk step `step` (of `length` letters) and its direction (+1 or
        -1 along the reference), or direction 0 if the contig does not visit its node once."""
        var pl = self.contig.placed(step >> 1)
        if pl < 0:
            return (0, 0)
        var rev = step & 1 == 1
        var node_rev = pl & 1 == 1
        var start = pl >> 1
        var fo = length - 1 - off if rev else off
        var r = start + (length - 1 - fo) if node_rev else start + fo
        return (r, 1 if rev == node_rev else -1)

    def touches(self, p: ImmPointer[UInt8, _], n: Int, mut fields: List[Int], mut walk: List[Int]) -> Bool:
        """Whether GAF line `p[:n]` (no newline) may place a base on the contig: its walk holds a node the contig
        visits once (`surject_line` places none otherwise)."""
        _split(p, n, fields)
        var k = _first_column(p, fields)
        if len(fields) // 2 < k + 12 or not _parse_walk(p, fields, k + 5, walk):
            return False
        for step in walk:
            if self.contig.placed(step >> 1) >= 0:
                return True
        return False

    def surject_line(self, p: ImmPointer[UInt8, _], n: Int, mut s: _Scratch, mut out: List[UInt8]) -> Bool:
        """The SAM record (with its newline) of GAF line `p[:n]` (no newline) appended to `out`; false, nothing
        appended, if it places no base on the contig."""
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
        if not _parse_walk(p, s.fields, k + 5, s.walk):
            return False
        # No base can be placed unless a node of the walk is the contig's.
        var any_placed = False
        for step in s.walk:
            any_placed = any_placed or self.contig.placed(step >> 1) >= 0
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
                    var letters = self.nodes.bases(step >> 1)
                    var take = min(left, length - off)
                    var at = len(s.read)
                    s.read.resize(unsafe_uninit_length=at + take)
                    var dst = s.read.unsafe_ptr().unsafe_offset(at)
                    if step & 1 == 1:
                        var last = length - off - 1
                        for t in range(take):
                            dst[unsafe_offset=t] = comp(letters[unsafe_offset=last - t])
                    else:
                        for t in range(take):
                            dst[unsafe_offset=t] = upper(letters[unsafe_offset=off + t])
                    s.score.resize(at + take, CLIP_MATCH)
                    var pl = self._place(step, length, off)
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
                var pl = self._place(s.walk[w], s.lens[w], off)
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
            s.read.reverse()
            for t in range(rl):
                s.read[t] = comp(s.read[t])
        var pos = 0
        if not linearize(Span(s.read), s.merged, Span(self.contig.seq), self.graph, self.windows, pos, s.cigar):
            return False
        # Clipping on the read's alignment to the walk, in the BAM's orientation.
        if not fwd:
            s.score.reverse()
        var clips = clip_lengths(s.score)
        var nc = len(s.cigar)
        var lead = s.cigar[0].length if nc > 0 and s.cigar[0].op == S else 0
        var trail = s.cigar[nc - 1].length if nc > 1 and s.cigar[nc - 1].op == S else 0
        var clipped = clip_cigar(pos, s.cigar, max(clips[0], lead), max(clips[1], trail))
        put_bytes(out, p, s.fields[2 * k], s.fields[2 * k + 1] - s.fields[2 * k])
        if fwd:
            put_text(out, "\t0\t")
        else:
            put_text(out, "\t16\t")
        out.extend(self.sam_contig.as_bytes())
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
        # bq covers the whole read; the aligned part (qstart..qend, on the read as sequenced) is what was rebuilt.
        # bq runs along the read as sequenced; the walk's letters along the read reversed when it is on `-`.
        var reverse_qual = minus == fwd
        var bq_at = s.fields[2 * bq_field] + 5 if bq_field >= 0 else 0
        var bq_len = s.fields[2 * bq_field + 1] - bq_at if bq_field >= 0 else -1
        if bq_len == qlen and qend <= qlen and qend >= qstart and qend - qstart == rl:
            if reverse_qual:
                for t in range(bq_at + qend - 1, bq_at + qstart - 1, -1):
                    out.append(p[unsafe_offset=t])
            else:
                put_bytes(out, p, bq_at + qstart, rl)
        else:
            for _ in range(rl):
                out.append(_QUAL)
        out.append(_NEWLINE)
        return True


comptime SLICES_PER_THREAD = 16
"""Slices of a chunk's messages: each copies its share of the reads to align into a piece of the batch."""
comptime BATCH_BYTES = 256 << 20
"""The lines a batch gathers before it is aligned. A read's alignment costs from nothing to a second, and a step waits
for its slowest task: aligned a few chunks at a time, the threads wait at the end of a batch, not of every chunk."""


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


struct _Piece(Movable):
    """A slice's reads of one chunk: the GAF lines of those that touch the contig (newline-ended), then their
    records."""

    var lines: List[UInt8]
    var text: List[UInt8]
    var written: Int

    def __init__(out self):
        self.lines = List[UInt8]()
        self.text = List[UInt8]()
        self.written = 0


struct _Batch(Movable):
    """Pieces in input order: `used` of them filled, `bytes` of lines in all."""

    var pieces: List[_Piece]
    var used: Int
    var bytes: Int

    def __init__(out self):
        self.pieces = List[_Piece]()
        self.used = 0
        self.bytes = 0


comptime _FREE = -1


@fieldwise_init
struct _Align(ParallelTask):
    """The sink's alignment workers, run once more at the end for the last batch."""

    var sink: MutPointer[SurjectSink, MutUntrackedOrigin]
    var failed: MutPointer[Bool, MutUntrackedOrigin]

    def run(self, index: Int):
        try:
            self.sink[].side_task(index)
        except:
            self.failed[] = True


struct SurjectSink(MessageSink, Movable):
    """Each slice gathers the GAF lines of its reads that touch the contig (a GAM's alignments written as GAF) into a
    piece of the filling batch. A full batch (or the last) is aligned in the next step by side tasks, one a thread,
    each taking pieces in turn and writing their records; in the step after, one more side task writes the records
    out in order (on stdout). Three batches take turns: filling, aligning, writing."""

    var s: Surjector
    var states: List[_Slice]
    var workers: List[_Scratch]
    var batches: List[_Batch]
    var filling: Int
    var aligning: Int
    """The batch the next step aligns, or _FREE."""
    var writing: Int
    """The batch the next step writes out, or _FREE."""
    var next_piece: Atomic[Int64]
    var lines: Bool
    var written: Int
    """Records made, over the batches aligned."""
    var failed: Bool

    def __init__(out self, var s: Surjector, threads: Int):
        self.s = s^
        self.states = List[_Slice]()
        for _ in range(SLICES_PER_THREAD * threads):
            self.states.append(_Slice())
        self.workers = List[_Scratch]()
        for _ in range(threads):
            self.workers.append(_Scratch())
        self.batches = List[_Batch]()
        for _ in range(3):
            self.batches.append(_Batch())
        self.filling = 0
        self.aligning = _FREE
        self.writing = _FREE
        self.next_piece = Atomic[Int64](0)
        self.lines = False
        self.written = 0
        self.failed = False
        self._room()
        var header = List[UInt8]()
        put_text(header, "@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:")
        header.extend(self.s.sam_contig.as_bytes())
        put_text(header, "\tLN:")
        put_int(header, len(self.s.contig.seq))
        put_text(header, "\n@PG\tID:dromos-surject\tPN:dromos\tVN:0.1.0\n")
        write_all(1, header.unsafe_ptr(), len(header))

    def _room(mut self):
        """Pieces enough in the filling batch for the next step's slices."""
        ref b = self.batches[self.filling]
        while len(b.pieces) < b.used + len(self.states):
            b.pieces.append(_Piece())

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
        ref b = self.batches[self.filling]
        var piece = mut_untracked(Pointer(to=b.pieces[b.used + slice]))
        piece[].lines.clear()
        piece[].text.clear()
        piece[].written = 0
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
            # No base can be placed unless a node of the alignment is the contig's.
            var any_placed = False
            for map in a.maps:
                any_placed = any_placed or self.s.contig.placed(map.node) >= 0
            if any_placed:
                a.decode(data, at, stop)
                gaf_line(data, a, self.s.nodes, piece[].lines)

    def side_tasks(self) -> Int:
        return (len(self.workers) if self.aligning != _FREE else 0) + (1 if self.writing != _FREE else 0)

    def side_task(mut self, index: Int) raises:
        if index >= len(self.workers) or self.aligning == _FREE:
            self._write(self.writing)
            return
        # A worker: the batch's pieces, one at a time, until none is left.
        var scratch = mut_untracked(Pointer(to=self.workers[index]))
        ref b = self.batches[self.aligning]
        while True:
            var k = Int(self.next_piece.fetch_add(1))
            if k >= b.used:
                return
            var piece = mut_untracked(Pointer(to=b.pieces[k]))
            var p = piece[].lines.unsafe_ptr().as_imm()
            var n = len(piece[].lines)
            var at = 0
            while at < n:
                var hit = external_call["memchr", Int](p.unsafe_offset(at), Int32(ord("\n")), n - at)
                var end = hit - Int(p)
                if self.s.surject_line(p.unsafe_offset(at), end - at, scratch[], piece[].text):
                    piece[].written += 1
                at = end + 1

    def _write(mut self, batch: Int):
        ref b = self.batches[batch]
        for k in range(b.used):
            write_all(1, b.pieces[k].text.unsafe_ptr(), len(b.pieces[k].text))

    def after_step(mut self, slices: Int) raises:
        # The batch aligned in the step just run is the next to write; the one written is free.
        var free = self.writing
        self.writing = self.aligning
        self.aligning = _FREE
        if self.writing != _FREE:
            ref done = self.batches[self.writing]
            for k in range(done.used):
                self.written += done.pieces[k].written
        if free != _FREE:
            self.batches[free].used = 0
            self.batches[free].bytes = 0
        ref b = self.batches[self.filling]
        for k in range(b.used, b.used + slices):
            b.bytes += len(b.pieces[k].lines)
        b.used += slices
        # A full batch is aligned next (the last, whatever its size, by `finish`).
        if b.bytes >= BATCH_BYTES:
            self.aligning = self.filling
            self.next_piece.store(0)
            for k in range(3):
                if k != self.aligning and k != self.writing:
                    self.filling = k
                    break
        self._room()

    def finish(mut self, contig: String, threads: Int, complete: Bool) raises:
        """Writes what was aligned and, once the input is `complete`, aligns and writes the last batch; then the
        counts. Also after an error (not `complete`), so the records made so far are out."""
        if self.writing != _FREE:
            self._write(self.writing)
            self.writing = _FREE
        if complete and self.aligning == _FREE and self.batches[self.filling].used > 0:
            self.aligning = self.filling
            self.next_piece.store(0)
        if complete and self.aligning != _FREE:
            var task = _Align(mut_untracked(Pointer(to=self)), mut_untracked(Pointer(to=self.failed)))
            parallel_for(task, len(self.workers), threads)
            if self.failed:
                raise Error("a read could not be aligned")
            self.after_step(0)
            self._write(self.writing)
            self.writing = _FREE
        var n = 0
        for t in range(len(self.states)):
            n += self.states[t].alignments
        print(
            t"[dromos] surject: {n} alignments, {self.written} written, {n - self.written} without a base on {contig}",
            file=stderr,
        )
