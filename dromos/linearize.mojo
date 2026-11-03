"""A read's blocks on a reference contig as one alignment to it (pgcatalog's surject.rs, ported line for line so
`dromos surject` writes the CIGARs pgcatalog writes): the banded affine-gap alignment (`Banded`, `align_global`,
`align_end`), the chain of blocks (`chain`), their joining (`linearize`), and the clipping of read ends that diverge
from the read's own graph alignment (`clip_lengths`, `clip_cigar`).

Between two blocks, a gap as long in the read as on the reference (a substitution bubble) is aligned base for base,
one with nothing on one side is an insertion or a deletion, and any other (an indel bubble) is aligned globally
(affine gaps, in a band around the diagonal as wide as the stretches' length difference and BAND_SLACK more; past
MAX_CELLS band cells, an insertion then a deletion). The read's bases before the first block and after the last are
aligned to the reference beyond them, its far end free, and kept if at least END_MATCHES of them match; else
soft-clipped. With `graph` placement a stretch of more than SMALL_BUBBLE bases between blocks is an insertion and a
deletion, and such a read end is clipped. A catalog window both of whose anchors the chain places is aligned as an
allele, at ASM_COSTS.
"""

from std.ffi import external_call

comptime M = UInt8(ord("M"))
comptime I = UInt8(ord("I"))
comptime D = UInt8(ord("D"))
comptime S = UInt8(ord("S"))
comptime LOWER_A = UInt8(ord("a"))
comptime LOWER_Z = UInt8(ord("z"))

comptime MAX_CELLS = 60_000_000
"""Band cells of an alignment, at most."""
comptime BAND_SLACK = 128
"""The band beyond the two stretches' length difference."""
comptime END_BAND = 2_000
"""The band of a read end's alignment beyond BAND_SLACK, at most."""
comptime SMALL_BUBBLE = 50
"""The longest stretch off the reference graph placement still aligns."""
comptime END_MATCHES = 0.9
"""The share of a read end's bases its alignment must match to be kept rather than clipped."""
comptime INF = 536870911
"""i32::MAX / 4, as pgcatalog's band starts every cell."""


@fieldwise_init
struct Costs(TrivialRegisterPassable):
    """Alignment costs (a match costs nothing)."""

    var mismatch: Int
    var open: Int
    var extend: Int


comptime READ_COSTS = Costs(5, 6, 3)
"""A read against the reference, gaps cheap enough to follow HiFi indels."""
comptime ASM_COSTS = Costs(40, 78, 7)
"""An allele against the reference as the truth sets spell it (minimap2's asm5 scores as costs)."""


@fieldwise_init
struct CigarOp(Equatable, TrivialRegisterPassable):
    var op: UInt8
    var length: Int


@fieldwise_init
struct Block(Equatable, TrivialRegisterPassable):
    """A placed stretch of a read: read offset, reference position (0-based), length."""

    var q: Int
    var r: Int
    var n: Int

    @always_inline
    def less(self, other: Self) -> Bool:
        """Tuple order: read offset, reference position, length."""
        if self.q != other.q:
            return self.q < other.q
        if self.r != other.r:
            return self.r < other.r
        return self.n < other.n


@fieldwise_init
struct Window(Equatable, TrivialRegisterPassable):
    """A catalog window: reference bases [lo, hi), 0-based."""

    var lo: Int
    var hi: Int


@always_inline
def push(mut cigar: List[CigarOp], op: UInt8, n: Int):
    """Appends `n` of `op`, joining a run of the same op."""
    if n == 0:
        return
    var k = len(cigar)
    if k > 0 and cigar[k - 1].op == op:
        cigar[k - 1].length += n
        return
    cigar.append(CigarOp(op, n))


@always_inline
def upper(c: UInt8) -> UInt8:
    """An ASCII letter in capitals; any other byte as it is."""
    return c - 32 if c >= LOWER_A and c <= LOWER_Z else c


@always_inline
def same_letter(a: UInt8, b: UInt8) -> Bool:
    """Equal, ignoring ASCII case."""
    return upper(a) == upper(b)


struct Banded(Movable):
    """The affine-gap alignment of a read stretch `q` against a reference stretch `r` from both starts, in a band
    around the diagonal (|j - i| <= half), as one traceback byte per cell (cell (i, j) at `i * w + j - i + half`):
    bits 0-1 where the best score came from (0 a match or mismatch, 1 a deletion, 2 an insertion), bit 2 a deletion
    opened there, bit 3 an insertion opened there. The last row's best scores are kept for a free reference end."""

    var half: Int
    var w: Int
    var tb: List[UInt8]
    var last: List[Int32]

    def __init__(out self):
        self.half = 0
        self.w = 0
        self.tb = List[UInt8]()
        self.last = List[Int32]()

    def fill(mut self, q: Span[UInt8, _], r: Span[UInt8, _], half: Int, c: Costs) -> Bool:
        """Fills the band; false (nothing filled) past MAX_CELLS cells."""
        var n = len(q)
        var m = len(r)
        var w = 2 * half + 1
        if (n + 1) * w > MAX_CELLS:
            return False
        self.half = half
        self.w = w
        self.tb.clear()
        self.tb.resize((n + 1) * w, 0)
        # Rows by band offset d = j - i + half: (i, j - 1) is d - 1 in the same row, (i - 1, j) is d + 1 in the row
        # before, (i - 1, j - 1) is d in the row before. H and F are kept two rows at a time (row i at (i & 1) * w),
        # E one. A row reads only cells the row before wrote, so none is reset (pgcatalog resets each to INF).
        var hs = List[Int32](length=2 * w, fill=Int32(INF))
        var fs = List[Int32](length=2 * w, fill=Int32(INF))
        var es = List[Int32](length=w, fill=Int32(INF))
        var qu = List[UInt8](capacity=n)
        for x in q:
            qu.append(upper(x))
        var ru = List[UInt8](capacity=m)
        for x in r:
            ru.append(upper(x))
        var rp = ru.unsafe_ptr()
        var tb = self.tb.unsafe_ptr()
        var e_row = es.unsafe_ptr()
        var mismatch = Int32(c.mismatch)
        var extend = Int32(c.extend)
        var open_cost = Int32(c.open + c.extend)
        # Row 0: (0, 0) scores 0; the cells after it are a deletion from it.
        var h0 = hs.unsafe_ptr()
        var f0 = fs.unsafe_ptr()
        h0[unsafe_offset=half] = 0
        e_row[unsafe_offset=half] = Int32(INF)
        f0[unsafe_offset=half] = Int32(INF)
        for j in range(1, min(half, m) + 1):
            var d = j + half
            var opened = h0[unsafe_offset=d - 1] + open_cost
            var extended = e_row[unsafe_offset=d - 1] + extend
            var e = min(opened, extended)
            var bits = UInt8(4) if opened <= extended else UInt8(0)
            # diag and f are INF: e is taken when it is at most INF (always, past (0, 0)).
            var h: Int32
            if Int32(INF) <= e:
                h = Int32(INF)
            else:
                bits |= 1
                h = e
            h0[unsafe_offset=d] = h
            e_row[unsafe_offset=d] = e
            f0[unsafe_offset=d] = Int32(INF)
            tb[unsafe_offset=d] = bits
        for i in range(1, n + 1):
            var h_cur = hs.unsafe_ptr().unsafe_offset((i & 1) * w)
            var h_prev = hs.unsafe_ptr().unsafe_offset(((i + 1) & 1) * w)
            var f_cur = fs.unsafe_ptr().unsafe_offset((i & 1) * w)
            var f_prev = fs.unsafe_ptr().unsafe_offset(((i + 1) & 1) * w)
            var row = tb.unsafe_offset(i * w)
            var lo = max(i - half, 0)
            var hi = min(i + half, m)
            if lo > hi:
                # Past the reference's end (a free end's band): this row and the rest have no cells.
                break
            var qi = qu[i - 1]
            # The row's first cell has no left neighbour in the band (E is INF); from j = 0 no diagonal either.
            var d = lo + half - i
            var f = Int32(INF)
            var bits = UInt8(0)
            if d + 1 < w:
                var opened = h_prev[unsafe_offset=d + 1] + open_cost
                var extended = f_prev[unsafe_offset=d + 1] + extend
                f = min(opened, extended)
                bits = UInt8(8) if opened <= extended else UInt8(0)
            var diag = Int32(INF)
            if lo > 0:
                diag = h_prev[unsafe_offset=d] + (0 if qi == rp[unsafe_offset=lo - 1] else mismatch)
            var e_left = Int32(INF)
            var h_left: Int32
            if diag <= e_left and diag <= f:
                h_left = diag
            elif e_left <= f:
                bits |= 1
                h_left = e_left
            else:
                bits |= 2
                h_left = f
            h_cur[unsafe_offset=d] = h_left
            e_row[unsafe_offset=d] = e_left
            f_cur[unsafe_offset=d] = f
            row[unsafe_offset=d] = bits
            # The rest: the last band cell (d = w - 1) has no cell above it in the band (F is INF).
            var full = min(hi, i + half - 1)
            for j in range(lo + 1, hi + 1):
                d = j + half - i
                var opened_e = h_left + open_cost
                var extended_e = e_left + extend
                var e = min(opened_e, extended_e)
                bits = UInt8(4) if opened_e <= extended_e else UInt8(0)
                f = Int32(INF)
                if j <= full:
                    var opened_f = h_prev[unsafe_offset=d + 1] + open_cost
                    var extended_f = f_prev[unsafe_offset=d + 1] + extend
                    f = min(opened_f, extended_f)
                    bits |= UInt8(8) if opened_f <= extended_f else UInt8(0)
                diag = h_prev[unsafe_offset=d] + (0 if qi == rp[unsafe_offset=j - 1] else mismatch)
                var h: Int32
                if diag <= e and diag <= f:
                    h = diag
                elif e <= f:
                    bits |= 1
                    h = e
                else:
                    bits |= 2
                    h = f
                h_cur[unsafe_offset=d] = h
                e_row[unsafe_offset=d] = e
                f_cur[unsafe_offset=d] = f
                row[unsafe_offset=d] = bits
                h_left = h
                e_left = e
        self.last = List[Int32](Span(hs)[(n & 1) * w : (n & 1) * w + w])
        return True

    @always_inline
    def last_score(self, n: Int, j: Int) -> Optional[Int]:
        """The best score at (n, j) of the last row, if j is in the band."""
        var d = j - n + self.half
        if d >= 0 and d < self.w:
            return Int(self.last[d])
        return None

    def trace(self, i0: Int, j0: Int) -> List[CigarOp]:
        """The operations from (0, 0) to (i0, j0); on ties a match or mismatch, so gaps go as far left as the costs
        allow."""
        var i = i0
        var j = j0
        var state = 0  # 0 best, 1 in a deletion, 2 in an insertion
        var rev = List[UInt8](capacity=i + j)
        while i > 0 or j > 0:
            var d = j - i + self.half
            var bits = self.tb[i * self.w + d] if d >= 0 and d < self.w else UInt8(0)
            if state == 0:
                var src = bits & 3
                if src == 0:
                    rev.append(M)
                    i -= 1
                    j -= 1
                elif src == 1:
                    state = 1
                else:
                    state = 2
            elif state == 1:
                rev.append(D)
                j -= 1
                if bits & 4 != 0:
                    state = 0
            else:
                rev.append(I)
                i -= 1
                if bits & 8 != 0:
                    state = 0
        var out = List[CigarOp]()
        for k in range(len(rev) - 1, -1, -1):
            push(out, rev[k], 1)
        return out^


@always_inline
def half_band(n: Int, m: Int) -> Int:
    """The band half-width an alignment of an `n`-letter and an `m`-letter stretch needs: their difference, and
    BAND_SLACK more for indels that cancel out."""
    return abs(n - m) + BAND_SLACK


def align_global(q: Span[UInt8, _], r: Span[UInt8, _], c: Costs = READ_COSTS) -> List[CigarOp]:
    """The global alignment of `q` (read) against `r` (reference) with affine gaps, as CIGAR operations; past
    MAX_CELLS band cells, an insertion then a deletion."""
    var n = len(q)
    var m = len(r)
    if n > 0 and m > 0:
        var b = Banded()
        if b.fill(q, r, half_band(n, m), c):
            return b.trace(n, m)
    var out = List[CigarOp]()
    push(out, I, n)
    push(out, D, m)
    return out^


def align_end(q: Span[UInt8, _], r: Span[UInt8, _], mut ops: List[CigarOp]) -> Int:
    """`q` aligned from its start against `r` from its start, `r`'s far end free: the CIGAR into `ops` and how many
    letters of `r` it uses, or -1 if it matches fewer than END_MATCHES of `q`'s letters (or is too large to align)."""
    ops.clear()
    var n = len(q)
    var m = len(r)
    if n == 0:
        return 0
    if m == 0:
        return -1
    # The far end is free, so the band need only hold the indels a read end may have, not the reference's slack.
    var b = Banded()
    if not b.fill(q, r, min(abs(n - m), END_BAND) + BAND_SLACK, READ_COSTS):
        return -1
    var end = -1
    var best = 0
    for j in range(1, m + 1):
        if case .Some(s) = b.last_score(n, j):
            if end < 0 or s < best:
                best = s
                end = j
    if end < 0:
        return -1
    ops = b.trace(n, end)
    var matches = 0
    var i = 0
    var j = 0
    for o in ops:
        for _ in range(o.length):
            if o.op == M:
                matches += Int(same_letter(q[i], r[j]))
                i += 1
                j += 1
            elif o.op == I:
                i += 1
            else:
                j += 1
    if Float64(matches) >= END_MATCHES * Float64(n):
        return end
    return -1


comptime CHAIN_LOOKBACK = 64
"""Predecessors a block's chaining looks back over."""
comptime CHAIN_MAX_STEP = 100_000
"""The farthest step in read or reference a join may cross."""


@always_inline
def join_cost(d: Int) -> Float64:
    """The cost of joining two blocks whose read and reference steps differ by `d` bases: minimap2's 0.01 w |d| +
    0.5 log2 |d| at w 5. libm's log2, as Rust's f64::log2 calls it, so scores tie where pgcatalog's do."""
    if d == 0:
        return 0.0
    var x = Float64(d)
    return 0.05 * x + 0.5 * external_call["log2", Float64](x)


def chain(blocks: List[Block]) -> List[Block]:
    """The best chain of blocks (sorted by read offset): each later than the one before in both read and reference,
    with the most bases less the join costs. The chain's blocks, in order."""
    var n = len(blocks)
    if n <= 1:
        return blocks.copy()
    var score = List[Float64](length=n, fill=0.0)
    var prev = List[Int](length=n, fill=-1)
    for j in range(n):
        var bj = blocks[j]
        score[j] = Float64(bj.n)
        for i in range(max(j - CHAIN_LOOKBACK, 0), j):
            var bi = blocks[i]
            var qe = bi.q + bi.n
            var re = bi.r + bi.n
            if qe > bj.q or re > bj.r:
                continue
            var dq = bj.q - qe
            var dr = bj.r - re
            if dq > CHAIN_MAX_STEP or dr > CHAIN_MAX_STEP:
                continue
            var s = score[i] + Float64(bj.n) - join_cost(abs(dq - dr))
            if s > score[j]:
                score[j] = s
                prev[j] = i
    # The best score; on ties the first.
    var j = 0
    for k in range(1, n):
        if score[k] > score[j]:
            j = k
    var out = List[Block]()
    out.append(blocks[j])
    while prev[j] >= 0:
        j = prev[j]
        out.append(blocks[j])
    out.reverse()
    return out^


# End clipping, as BWA-MEM scores a local alignment, on the read's alignment to the graph (its cs edits): a base
# matching the walk +1, a mismatch -4, an indel -(6 + its length) (free beyond SMALL_BUBBLE: a structural variant),
# and an end clipped only where that gains more than CLIP_PENALTY.
comptime CLIP_MATCH = 1
comptime CLIP_MISMATCH = -4
comptime CLIP_GAP_OPEN = -6
comptime CLIP_PENALTY = 5


@always_inline
def clip_gap(n: Int) -> Int:
    """The score of an indel of `n` bases at the clipping scores."""
    return 0 if n > SMALL_BUBBLE else CLIP_GAP_OPEN - n


def clip_lengths(score: List[Int]) -> Tuple[Int, Int]:
    """Read bases to clip at each end, from per-base scores (an indel's score on the base before it): the cut before
    a matching base with the lowest prefix score, and after one with the lowest suffix score, where that gains more
    than CLIP_PENALTY."""
    var n = len(score)
    var p = score.unsafe_ptr()
    var total = 0
    var best = 0
    var left = 0
    for i in range(n):
        var v = p[unsafe_offset=i]
        if v == CLIP_MATCH and total < best:
            best = total
            left = i
        total += v
    if -best <= CLIP_PENALTY:
        left = 0
    total = 0
    best = 0
    var right = n
    for i in range(n - 1, left - 1, -1):
        var v = p[unsafe_offset=i]
        if v == CLIP_MATCH and total < best:
            best = total
            right = i + 1
        total += v
    if -best <= CLIP_PENALTY:
        right = n
    return (left, n - right)


def _remerge(runs: List[CigarOp]) -> List[CigarOp]:
    var merged = List[CigarOp]()
    for o in runs:
        push(merged, o.op, o.length)
    return merged^


def clip_cigar(start: Int, cigar: List[CigarOp], head: Int, tail: Int) -> Tuple[Int, List[CigarOp]]:
    """The alignment (start, CIGAR) with `head` read bases soft-clipped at its start and `tail` at its end (soft clips
    already there included in the counts); a gap left at either end goes with the clip."""
    if head == 0 and tail == 0:
        return (start, cigar.copy())
    var qlen = 0
    for o in cigar:
        if o.op != D:
            qlen += o.length
    if head + tail >= qlen:
        return (start, cigar.copy())
    # One op per read base or deleted base, cut, then run-length again.
    var q = 0
    var pos = start
    var runs = List[CigarOp]()
    for o in cigar:
        var op = o.op
        var read_base = op != D
        for _ in range(o.length):
            var in_head = q < head
            var in_tail = q >= qlen - tail
            if in_head or (in_tail and read_base):
                if read_base:
                    push(runs, S, 1)
                if in_head and (op == M or op == D):
                    pos += 1
            elif not (in_tail and op == D):
                # (a deletion after the last kept base is dropped)
                push(runs, op, 1)
            if read_base:
                q += 1
    # A clip ends at an aligned base: gaps beside a soft clip join it (insertions) or move the start (deletions).
    while True:
        var k = 0
        while k < len(runs) and runs[k].op == S:
            k += 1
        if k >= len(runs):
            break
        var o = runs[k]
        if o.op == I:
            _ = runs.pop(k)
            if k > 0:
                runs[k - 1].length += o.length
            else:
                runs.insert(0, CigarOp(S, o.length))
        elif o.op == D:
            _ = runs.pop(k)
            pos += o.length
        else:
            break
        runs = _remerge(runs)
    while True:
        var k = len(runs) - 1
        while k >= 0 and runs[k].op == S:
            k -= 1
        if k < 0:
            break
        var o = runs[k]
        if o.op == I:
            _ = runs.pop(k)
            if k < len(runs):
                runs[k].length += o.length
            else:
                runs.append(CigarOp(S, o.length))
        elif o.op == D:
            _ = runs.pop(k)
        else:
            break
        runs = _remerge(runs)
    return (pos, runs^)


def _reversed(s: Span[UInt8, _]) -> List[UInt8]:
    var out = List[UInt8](capacity=len(s))
    for k in range(len(s) - 1, -1, -1):
        out.append(s[k])
    return out^


def _covers(kept: List[Block], r: Int) -> Bool:
    for b in kept:
        if b.r <= r and r < b.r + b.n:
            return True
    return False


def linearize(
    query: Span[UInt8, _],
    blocks: List[Block],
    reference: Span[UInt8, _],
    graph: Bool,
    windows: List[Window],
    mut start: Int,
    mut cigar: List[CigarOp],
) -> Bool:
    """The placed blocks of read `query` (the BAM's orientation, sorted by read offset) as an alignment to `reference`
    (the contig's letters): its first reference position (0-based) into `start` and CIGAR into `cigar`; false if no
    block. `graph`: graph placement. The read is aligned at each window (sorted, disjoint) it holds both anchors of
    (the reference bases just before and just after it) as its allele, at ASM_COSTS."""
    cigar.clear()
    var ref_len = len(reference)
    var kept = chain(blocks)
    # The windows whose two anchors the chain places: their bases taken out of the blocks, so each is one gap.
    var spanned = List[Window]()
    if len(kept) > 0:
        var lo = kept[0].r
        var hi = kept[len(kept) - 1].r + kept[len(kept) - 1].n
        # The first window not wholly before lo (windows[k].hi < lo for every earlier k).
        var a = 0
        var b = len(windows)
        while a < b:
            var mid = (a + b) // 2
            if windows[mid].hi < lo:
                a = mid + 1
            else:
                b = mid
        for k in range(a, len(windows)):
            var w = windows[k]
            if w.lo > hi:
                break
            if w.lo > 0 and _covers(kept, w.lo - 1) and _covers(kept, w.hi):
                spanned.append(w)
    if len(spanned) > 0:
        var cut = List[Block](capacity=len(kept) + len(spanned))
        for blk in kept:
            var parts = List[Block]()
            parts.append(blk)
            for w in spanned:
                var next = List[Block]()
                for p in parts:
                    var e = p.r + p.n
                    if e <= w.lo or p.r >= w.hi:
                        next.append(p)
                        continue
                    if p.r < w.lo:
                        next.append(Block(p.q, p.r, w.lo - p.r))
                    if e > w.hi:
                        next.append(Block(p.q + (w.hi - p.r), w.hi, e - w.hi))
                parts = next^
            cut.extend(parts^)
        kept = cut^
    if len(kept) == 0:
        return False
    var first = kept[0]
    start = first.r
    var ops = List[CigarOp]()
    # The head: the read's first bases, reversed, against the reference before the first block, reversed.
    var head = first.q
    var reach = head + max(head // 5, 50)
    var r0 = max(first.r - reach, 0)
    var used = -1
    if head > 0 and first.r > 0 and not (graph and head > SMALL_BUBBLE):
        used = align_end(_reversed(query[:head]), _reversed(reference[r0 : first.r]), ops)
    if used >= 0:
        start = first.r - used
        for k in range(len(ops) - 1, -1, -1):
            push(cigar, ops[k].op, ops[k].length)
    else:
        push(cigar, S, head)
    push(cigar, M, first.n)
    for k in range(1, len(kept)):
        var a = kept[k - 1]
        var b = kept[k]
        var qa = a.q + a.n
        var ra = a.r + a.n
        var gq = b.q - qa
        var gr = b.r - ra
        if Window(ra, b.r) in spanned:
            for o in align_global(query[qa : b.q], reference[ra : b.r], ASM_COSTS):
                push(cigar, o.op, o.length)
        elif graph and max(gq, gr) > SMALL_BUBBLE:
            push(cigar, I, gq)
            push(cigar, D, gr)
        elif gq == gr:
            push(cigar, M, gq)
        elif gr == 0:
            push(cigar, I, gq)
        elif gq == 0:
            push(cigar, D, gr)
        else:
            for o in align_global(query[qa : b.q], reference[ra : b.r]):
                push(cigar, o.op, o.length)
        push(cigar, M, b.n)
    var last = kept[len(kept) - 1]
    var tq = last.q + last.n
    var tr = last.r + last.n
    var tail = len(query) - tq
    var tail_reach = min(tail + max(tail // 5, 50), max(ref_len - tr, 0))
    used = -1
    if tail > 0 and tail_reach > 0 and not (graph and tail > SMALL_BUBBLE):
        used = align_end(query[tq:], reference[tr : tr + tail_reach], ops)
    if used >= 0:
        for o in ops:
            push(cigar, o.op, o.length)
    else:
        push(cigar, S, tail)
    # An end alignment may start or finish with a deletion: not a CIGAR's first or last operation.
    var lead = 0
    while lead < len(cigar) and cigar[lead].op == D:
        start += cigar[lead].length
        lead += 1
    if lead > 0:
        var rest = List[CigarOp](capacity=len(cigar) - lead)
        for k in range(lead, len(cigar)):
            rest.append(cigar[k])
        cigar = rest^
    while len(cigar) > 0 and cigar[len(cigar) - 1].op == D:
        _ = cigar.pop()
    return True
