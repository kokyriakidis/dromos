"""A GBZ graph as the snarl finder reads it: its nodes and edges split into weakly connected components, each a compact
bidirected graph (`Bigraph`): nodes by rank (ascending id), their lengths, and every handle's successors in one CSR.

A handle is `2 * rank + reverse` within its component; a visit (node id and orientation, as the GBWT encodes it) is
`id << 1 | reverse`. The edges come from the GBWT's records through gbz-ffi (`GbzParts.node_edges`), read on every
thread by ranges of node ids.
"""

from std.memory import Layout, alloc, dealloc
from std.sys import stderr
from std.time import perf_counter_ns

from dromos.flat import Flat
from dromos.graph import GbzParts
from dromos.parallel import ParallelTask, parallel_for

comptime NO_RANK = UInt32.MAX
"""A rank, handle or component for none."""
comptime ID_SPAN = 1 << 16
"""Node ids one read task takes."""


struct Bigraph(Movable):
    """One weakly connected component of a bidirected graph: `ids[r]` is node rank r's id (ascending), `lengths[r]`
    its sequence length, and handle h's successors are `targets[offsets[h] ..< offsets[h + 1]]` (handles: an edge
    a -> b is in a's list, and its other reading, flip(b) -> flip(a), in flip(b)'s)."""

    var ids: Flat[UInt64]
    var lengths: Flat[UInt32]
    var offsets: Flat[UInt32]
    var targets: Flat[UInt32]

    def __init__(out self):
        self.ids = Flat[UInt64]()
        self.lengths = Flat[UInt32]()
        self.offsets = Flat[UInt32]()
        self.targets = Flat[UInt32]()

    def __init__(
        out self, var ids: Flat[UInt64], var lengths: Flat[UInt32], var offsets: Flat[UInt32], var targets: Flat[UInt32]
    ):
        self.ids = ids^
        self.lengths = lengths^
        self.offsets = offsets^
        self.targets = targets^

    @always_inline
    def nodes(self) -> Int:
        return len(self.ids)

    @always_inline
    def visit(self, handle: UInt32) -> UInt64:
        """A handle as a visit: `id << 1 | reverse`."""
        return self.ids[Int(handle >> 1)] << 1 | UInt64(handle & 1)


struct _IdRange(Movable):
    """What one read task found for node ids `first ..< first + count`: per id, 2 successor counts (forward then
    reverse; absent ids have none) and its length (NO_RANK when absent), and the successors as global handles
    (`2 * (id - min id) + reverse`)."""

    var starts: Flat[UInt32]
    """Handle `first_handle + k`'s successors are `succ[starts[k] ..< starts[k + 1]]`."""
    var succ: Flat[UInt32]
    var lengths: Flat[UInt32]
    var cross: Flat[UInt32]
    """Edges (x, y) to a higher id outside the range, as pairs of id offsets: joined after the ranges' own."""
    var runs: Flat[UInt32]
    """The range's nodes as runs of one component: (component, nodes, first rank) triples."""
    var roots: Int

    def __init__(out self):
        self.starts = Flat[UInt32]()
        self.succ = Flat[UInt32]()
        self.lengths = Flat[UInt32]()
        self.cross = Flat[UInt32]()
        self.runs = Flat[UInt32]()
        self.roots = 0


@fieldwise_init
struct _ReadTask(ParallelTask):
    var parts: ImmPointer[GbzParts, ImmUntrackedOrigin]
    var ranges: MutPointer[_IdRange, MutUntrackedOrigin]
    var lo: Int
    var hi: Int
    var bad: MutPointer[Int, MutUntrackedOrigin]

    def run(self, index: Int):
        var first = self.lo + index * ID_SPAN
        var count = min(ID_SPAN, self.hi + 1 - first)
        var r = _IdRange()
        try:
            var raw = List[UInt64]()
            var counts = List[Int64]()
            r.lengths = Flat[UInt32](unsafe_uninit_length=count)
            self.parts[].node_edges(first, count, raw, counts, r.lengths.unsafe_mut_ptr())
            r.starts = Flat[UInt32](unsafe_uninit_length=2 * count + 1)
            r.succ = Flat[UInt32](unsafe_uninit_length=len(raw))
            var at = 0
            var rawp = raw.unsafe_ptr()
            var countp = counts.unsafe_ptr()
            for k in range(2 * count):
                r.starts[k] = UInt32(at)
                for _ in range(max(Int(countp[unsafe_offset=k]), 0)):
                    var v = rawp[unsafe_offset=at]
                    r.succ[at] = UInt32(Int(v >> 1) - self.lo) << 1 | UInt32(v & 1)
                    at += 1
            r.starts[2 * count] = UInt32(at)
        except:
            self.bad[] = 1
        self.ranges[unsafe_offset=index] = r^


@always_inline
def _find(parent: MutPointer[UInt32, _], x: UInt32) -> UInt32:
    """The root of `x`, halving the path on the way."""
    var i = x
    while True:
        var p = parent[unsafe_offset=Int(i)]
        if p == i:
            return i
        var g = parent[unsafe_offset=Int(p)]
        parent[unsafe_offset=Int(i)] = g
        i = g


@fieldwise_init
struct _BuildTask(ParallelTask):
    """Component `order[index]`'s Bigraph, from the read ranges."""

    var ranges: ImmPointer[_IdRange, ImmUntrackedOrigin]
    var members: ImmPointer[UInt32, ImmUntrackedOrigin]
    """Every component's nodes (offsets from the min id), component by component, ascending."""
    var member_starts: ImmPointer[Int, ImmUntrackedOrigin]
    var rank: ImmPointer[UInt32, ImmUntrackedOrigin]
    var lo: Int
    var out: MutPointer[Bigraph, MutUntrackedOrigin]

    def run(self, index: Int):
        var a = self.member_starts[unsafe_offset=index]
        var b = self.member_starts[unsafe_offset=index + 1]
        var n = b - a
        var ids = Flat[UInt64](unsafe_uninit_length=n)
        var lengths = Flat[UInt32](unsafe_uninit_length=n)
        var offsets = Flat[UInt32](unsafe_uninit_length=2 * n + 1)
        var total = 0
        for r in range(n):
            var x = Int(self.members[unsafe_offset=a + r])
            var range_ = self.ranges.unsafe_offset(x // ID_SPAN)
            var k = x % ID_SPAN
            ids[r] = UInt64(self.lo + x)
            lengths[r] = range_[].lengths[k]
            for o in range(2):
                offsets[2 * r + o] = UInt32(total)
                total += Int(range_[].starts[2 * k + o + 1] - range_[].starts[2 * k + o])
        offsets[2 * n] = UInt32(total)
        var targets = Flat[UInt32](unsafe_uninit_length=total)
        var at = 0
        for r in range(n):
            var x = Int(self.members[unsafe_offset=a + r])
            var range_ = self.ranges.unsafe_offset(x // ID_SPAN)
            var k = x % ID_SPAN
            for j in range(Int(range_[].starts[2 * k]), Int(range_[].starts[2 * k + 2])):
                var g = range_[].succ[j]
                targets[at] = self.rank[unsafe_offset=Int(g >> 1)] << 1 | (g & 1)
                at += 1
        self.out[unsafe_offset=index] = Bigraph(ids^, lengths^, offsets^, targets^)


def rss_gb() -> Float64:
    """This process's resident memory (GB), from /proc/self/status; 0 where there is none."""
    try:
        with open("/proc/self/status", "r") as f:
            var text = f.read()
            var at = text.find("VmRSS:")
            if at < 0:
                return 0
            var rest = text[byte = at + 6 :]
            var kb = Int(String(rest.split("kB")[0].strip()))
            return Float64(kb) / 1048576.0
    except:
        return 0


def _log_time(verbose: Bool, what: StaticString, t0: Int) -> Int:
    var t = perf_counter_ns()
    if verbose:
        print(t"[bigraph] {what}: {Float64(t - t0) / 1e9} s, rss {rss_gb()} GB", file=stderr)
    return t


struct GraphEdges(Movable):
    """A graph's nodes and edges as read from the GBZ, by ranges of node ids (what `components` splits)."""

    var ranges: List[_IdRange]
    var lo: Int
    var hi: Int
    var count: Int

    def __init__(out self):
        self.ranges = List[_IdRange]()
        self.lo = 0
        self.hi = -1
        self.count = 0


def read_edges(parts: GbzParts, threads: Int, verbose: Bool = False) raises -> GraphEdges:
    """Every node's length and both orientations' successors, read from the GBZ's GBWT on `threads` threads."""
    var t = perf_counter_ns()
    var out = GraphEdges()
    var lo, hi = parts.id_bounds()
    if hi == 0:
        return out^
    var span = hi + 1 - lo
    if 2 * span >= Int(NO_RANK):
        raise Error(t"node ids {lo}..{hi} span too many for 32-bit handles")
    var nranges = (span + ID_SPAN - 1) // ID_SPAN
    out.ranges = List[_IdRange](capacity=nranges)
    for _ in range(nranges):
        out.ranges.append(_IdRange())
    var bad = 0
    var read = _ReadTask(
        Pointer(to=parts).unsafe_origin_cast[ImmUntrackedOrigin](),
        out.ranges.unsafe_ptr().unsafe_origin_cast[MutUntrackedOrigin](),
        lo,
        hi,
        Pointer(to=bad).unsafe_origin_cast[MutUntrackedOrigin](),
    )
    parallel_for(read, nranges, threads)
    _ = len(out.ranges)
    if bad != 0:
        raise Error("reading the graph's edges failed")
    out.lo = lo
    out.hi = hi
    var count = 0
    for q in range(nranges):
        ref lengths = out.ranges[q].lengths
        for k in range(len(lengths)):
            if lengths[k] != NO_RANK:
                count += 1
    out.count = count
    _ = _log_time(verbose, "edges read", t)
    return out^


@always_inline
def _root(parent: ImmPointer[UInt32, _], x: UInt32) -> UInt32:
    """The root of `x`, without changing the forest (other threads read it at the same time)."""
    var i = x
    while True:
        var p = parent[unsafe_offset=Int(i)]
        if p == i:
            return i
        i = p


@fieldwise_init
struct _LocalJoin(ParallelTask):
    """Range `index`'s edges within the range joined (each range writes only its own ids' entries); edges leaving
    it kept for later. The smaller root wins, so a root is its tree's smallest id."""

    var ranges: MutPointer[_IdRange, MutUntrackedOrigin]
    var parent: MutPointer[UInt32, MutUntrackedOrigin]

    def run(self, index: Int):
        var r = self.ranges.unsafe_offset(index)
        var base = index * ID_SPAN
        var m = len(r[].lengths)
        var limit = UInt32(base + m)
        for k in range(m):
            var x = UInt32(base + k)
            self.parent[unsafe_offset=Int(x)] = x
        for k in range(m):
            var x = UInt32(base + k)
            for j in range(Int(r[].starts[2 * k]), Int(r[].starts[2 * k + 2])):
                var y = r[].succ[j] >> 1
                if y <= x:
                    continue
                if y >= limit:
                    r[].cross.append(x)
                    r[].cross.append(y)
                    continue
                var rx = _find(self.parent, x)
                var ry = _find(self.parent, y)
                if rx < ry:
                    self.parent[unsafe_offset=Int(ry)] = rx
                elif ry < rx:
                    self.parent[unsafe_offset=Int(rx)] = ry


@fieldwise_init
struct _Roots(ParallelTask):
    """Each id's root (NO_RANK for ids the graph lacks), and how many roots the range has."""

    var ranges: MutPointer[_IdRange, MutUntrackedOrigin]
    var parent: ImmPointer[UInt32, ImmUntrackedOrigin]
    var root: MutPointer[UInt32, MutUntrackedOrigin]

    def run(self, index: Int):
        var r = self.ranges.unsafe_offset(index)
        var base = index * ID_SPAN
        var roots = 0
        for k in range(len(r[].lengths)):
            var x = UInt32(base + k)
            if r[].lengths[k] == NO_RANK:
                self.root[unsafe_offset=Int(x)] = NO_RANK
                continue
            var rt = _root(self.parent, x)
            self.root[unsafe_offset=Int(x)] = rt
            if rt == x:
                roots += 1
        r[].roots = roots


@fieldwise_init
struct _Number(ParallelTask):
    """Each root's component number (from the range's first), written at the root in `number`."""

    var ranges: ImmPointer[_IdRange, ImmUntrackedOrigin]
    var first: ImmPointer[Int, ImmUntrackedOrigin]
    var root: ImmPointer[UInt32, ImmUntrackedOrigin]
    var number: MutPointer[UInt32, MutUntrackedOrigin]

    def run(self, index: Int):
        var base = index * ID_SPAN
        var c = self.first[unsafe_offset=index]
        for k in range(len(self.ranges[unsafe_offset=index].lengths)):
            var x = base + k
            if self.root[unsafe_offset=x] == UInt32(x):
                self.number[unsafe_offset=x] = UInt32(c)
                c += 1


@fieldwise_init
struct _Label(ParallelTask):
    """Each id's component (in place of its root) and the range's runs of one component."""

    var ranges: MutPointer[_IdRange, MutUntrackedOrigin]
    var number: ImmPointer[UInt32, ImmUntrackedOrigin]
    var comp: MutPointer[UInt32, MutUntrackedOrigin]

    def run(self, index: Int):
        var r = self.ranges.unsafe_offset(index)
        var base = index * ID_SPAN
        var last = NO_RANK
        for k in range(len(r[].lengths)):
            var x = base + k
            var rt = self.comp[unsafe_offset=x]
            if rt == NO_RANK:
                continue
            var c = self.number[unsafe_offset=Int(rt)]
            self.comp[unsafe_offset=x] = c
            if c == last:
                var n = len(r[].runs)
                r[].runs[n - 2] += 1
            else:
                r[].runs.append(c)
                r[].runs.append(1)
                r[].runs.append(0)
                last = c


@fieldwise_init
struct _Rank(ParallelTask):
    """Each node's rank in its component, and its place in the members of all components."""

    var ranges: ImmPointer[_IdRange, ImmUntrackedOrigin]
    var comp: ImmPointer[UInt32, ImmUntrackedOrigin]
    var member_starts: ImmPointer[Int, ImmUntrackedOrigin]
    var rank: MutPointer[UInt32, MutUntrackedOrigin]
    var members: MutPointer[UInt32, MutUntrackedOrigin]

    def run(self, index: Int):
        var r = self.ranges.unsafe_offset(index)
        var base = index * ID_SPAN
        var run = 0
        var left = 0
        var next_rank = UInt32(0)
        for k in range(len(r[].lengths)):
            var x = base + k
            var c = self.comp[unsafe_offset=x]
            if c == NO_RANK:
                continue
            if left == 0:
                left = Int(r[].runs[run + 1])
                next_rank = r[].runs[run + 2]
                run += 3
            self.rank[unsafe_offset=x] = next_rank
            self.members[unsafe_offset=self.member_starts[unsafe_offset=Int(c)] + Int(next_rank)] = UInt32(x)
            next_rank += 1
            left -= 1


struct Components(Movable, Sized):
    """A graph's weakly connected components, each a Bigraph (`graphs`, in order of their smallest node id), and the
    component of each node id (`of`)."""

    var graphs: List[Bigraph]
    var label: Flat[UInt32]
    """Node id `lo + i`'s component, NO_RANK for an id the graph lacks."""
    var lo: Int

    def __init__(out self):
        self.graphs = List[Bigraph]()
        self.label = Flat[UInt32]()
        self.lo = 0

    def __init__(out self, var graphs: List[Bigraph], var label: Flat[UInt32], lo: Int):
        self.graphs = graphs^
        self.label = label^
        self.lo = lo

    def __len__(self) -> Int:
        return len(self.graphs)

    def of(self, id: Int) -> Optional[Int]:
        """The component holding node `id`; none if the graph lacks it."""
        var i = id - self.lo
        if i < 0 or i >= len(self.label) or self.label[i] == NO_RANK:
            return None
        return Int(self.label[i])


def components(var g: GraphEdges, threads: Int, verbose: Bool = False) raises -> Components:
    """The graph's weakly connected components, each a Bigraph, in order of their smallest node id, and each node's
    component. Every step but joining the edges between id ranges runs on `threads` threads."""
    if g.count == 0:
        return Components()
    var t = perf_counter_ns()
    var lo = g.lo
    var span = g.hi + 1 - lo
    var nranges = len(g.ranges)
    var rp = g.ranges.unsafe_ptr().unsafe_origin_cast[MutUntrackedOrigin]()
    var parent = Flat[UInt32](unsafe_uninit_length=span)
    var pp = parent.unsafe_mut_ptr().unsafe_origin_cast[MutUntrackedOrigin]()
    var local = _LocalJoin(rp, pp)
    parallel_for(local, nranges, threads)
    for q in range(nranges):
        ref cross = g.ranges[q].cross
        for j in range(0, len(cross), 2):
            var rx = _find(pp, cross[j])
            var ry = _find(pp, cross[j + 1])
            if rx < ry:
                pp[unsafe_offset=Int(ry)] = rx
            elif ry < rx:
                pp[unsafe_offset=Int(rx)] = ry
        cross = Flat[UInt32]()
    t = _log_time(verbose, "components joined", t)
    # Roots (each component's smallest id), numbered in order.
    var comp = Flat[UInt32](unsafe_uninit_length=span)
    var cp = comp.unsafe_mut_ptr().unsafe_origin_cast[MutUntrackedOrigin]()
    var roots = _Roots(rp, pp.as_imm().unsafe_origin_cast[ImmUntrackedOrigin](), cp)
    parallel_for(roots, nranges, threads)
    var first = Flat[Int](unsafe_uninit_length=nranges + 1)
    var ncomp = 0
    for q in range(nranges):
        first[q] = ncomp
        ncomp += g.ranges[q].roots
    var number = _Number(
        rp.as_imm().unsafe_origin_cast[ImmUntrackedOrigin](),
        first.unsafe_ptr().unsafe_origin_cast[ImmUntrackedOrigin](),
        cp.as_imm().unsafe_origin_cast[ImmUntrackedOrigin](),
        pp,
    )
    parallel_for(number, nranges, threads)
    var label = _Label(rp, pp.as_imm().unsafe_origin_cast[ImmUntrackedOrigin](), cp)
    parallel_for(label, nranges, threads)
    # Component sizes, and each run's first rank.
    var sizes = Flat[Int](length=ncomp, fill=0)
    for q in range(nranges):
        ref runs = g.ranges[q].runs
        for j in range(0, len(runs), 3):
            var c = Int(runs[j])
            runs[j + 2] = UInt32(sizes[c])
            sizes[c] += Int(runs[j + 1])
    var member_starts = Flat[Int](length=ncomp + 1, fill=0)
    for c in range(ncomp):
        member_starts[c + 1] = member_starts[c] + sizes[c]
    var members = Flat[UInt32](unsafe_uninit_length=g.count)
    var rank_task = _Rank(
        rp.as_imm().unsafe_origin_cast[ImmUntrackedOrigin](),
        cp.as_imm().unsafe_origin_cast[ImmUntrackedOrigin](),
        member_starts.unsafe_ptr().unsafe_origin_cast[ImmUntrackedOrigin](),
        pp,
        members.unsafe_mut_ptr().unsafe_origin_cast[MutUntrackedOrigin](),
    )
    parallel_for(rank_task, nranges, threads)
    _ = len(first) + len(sizes)
    t = _log_time(verbose, "components numbered", t)
    var out = List[Bigraph](capacity=ncomp)
    for _ in range(ncomp):
        out.append(Bigraph())
    var build = _BuildTask(
        rp.as_imm().unsafe_origin_cast[ImmUntrackedOrigin](),
        members.unsafe_ptr().unsafe_origin_cast[ImmUntrackedOrigin](),
        member_starts.unsafe_ptr().unsafe_origin_cast[ImmUntrackedOrigin](),
        pp.as_imm().unsafe_origin_cast[ImmUntrackedOrigin](),
        lo,
        out.unsafe_ptr().unsafe_origin_cast[MutUntrackedOrigin](),
    )
    parallel_for(build, ncomp, threads)
    # The tasks reached these through untracked pointers: uses after them keep them alive until they are done.
    _ = len(g.ranges) + len(members) + len(member_starts) + len(parent)
    _ = _log_time(verbose, "components built", t)
    return Components(out^, comp^, lo)


def gbz_components(parts: GbzParts, threads: Int, verbose: Bool = False) raises -> Components:
    """The graph's weakly connected components, each a Bigraph, in order of their smallest node id."""
    return components(read_edges(parts, threads, verbose), threads, verbose)
