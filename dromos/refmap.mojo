"""Alignments placed on the reference: pgcat's reference map, read from the graph's cache (cache.mojo: a reference
sample's PATHS, FIRST, PLACEMENTS and MORE sections, with the nodes' LENGTHS), and its placement rule (pgcat's
gaf.rs, `RefMap::project`, itself pggaf's walk_ref_interval), ported line for line so dromos places every alignment
where `pgcat index-gaf` does.
"""

from dromos.cache import (
    FIRST,
    LENGTHS,
    MORE,
    MORE_AT,
    MORE_NODES,
    NONE,
    PATH_BITS,
    PATHS,
    PLACEMENTS,
    START_SHIFT,
    CacheFile,
    u32_at,
    u64_at,
)
from std.math import ceildiv

from dromos.nodes import Nodes
from dromos.parallel import ParallelTask, parallel_for
from dromos.stream import imm_untracked, mut_untracked, untracked

# Why an alignment with a walk has no place (pgcat's Unplaced, in its order); 0 is placed.
comptime NO_ANCHOR = 1
"""No node of the aligned part of the walk lies on the reference."""
comptime MISSING_NODE = 2
"""A node of the walk is not in the graph."""
comptime AMBIGUOUS = 3
"""The anchors fit the reference at more than one place."""
comptime INCONSISTENT = 4
"""The anchors fit no single path in one orientation and order (an inversion, a rearrangement)."""


def reason_name(why: Int) -> StaticString:
    """A reason's name, as pgcat writes it."""
    if why == NO_ANCHOR:
        return "no_anchor"
    if why == MISSING_NODE:
        return "missing_node"
    if why == AMBIGUOUS:
        return "ambiguous"
    return "inconsistent"


@fieldwise_init
struct Anchor(TrivialRegisterPassable):
    """A walk node on the reference, and the part of it the alignment covers."""

    var node: Int
    var reverse: Bool
    var local_start: Int
    var local_end: Int


@fieldwise_init
struct Placement(TrivialRegisterPassable):
    """A node's place on a reference path."""

    var path: Int
    var step: Int
    var start: Int
    var reverse: Bool


@fieldwise_init
struct Placed(Equatable, TrivialRegisterPassable):
    """A reference interval: contig, start and end (0-based, end exclusive)."""

    var contig: Int
    var start: Int
    var end: Int


@fieldwise_init
struct ContigPath(Copyable):
    """A reference path of one contig: its index, its first placement (step 0's), and each step's node."""

    var path: Int
    var first: Int
    var nodes: List[Int]


comptime _NODES_PER_TASK = 1 << 22


@fieldwise_init
struct _StepsOf(ParallelTask):
    """`RefMap.set_steps` for each range of node ids, and once more for the further placements."""

    var map: ImmPointer[RefMap, ImmUntrackedOrigin]
    var out: MutPointer[List[List[ContigPath]], MutUntrackedOrigin]
    var slot: ImmPointer[List[Int], ImmUntrackedOrigin]
    var ranges: Int

    def run(self, index: Int):
        self.map[].set_steps(self.out[], self.slot[], index, self.ranges)


@always_inline
def _set_step(mut paths: List[ContigPath], k: Int, node: Int):
    """Placement `k` is of `node`: its step's node, if `k` is a step of one of `paths`."""
    for ref cp in paths:
        if k >= cp.first and k < cp.first + len(cp.nodes):
            cp.nodes[k - cp.first] = node
            return


struct RefMap(Movable, Nodes):
    """A reference sample's placements, mapped from the graph's cache: node lengths and placements are read in
    place."""

    var cache: CacheFile
    var sample: String
    var size: Int
    var node_len: ImmPointer[UInt32, ImmUntrackedOrigin]
    var first: ImmPointer[UInt32, ImmUntrackedOrigin]
    var placements: ImmPointer[UInt64, ImmUntrackedOrigin]
    var placement_count: Int
    var path_first: List[Int]
    """Each path's first placement, and the end of the last."""
    var more_nodes: ImmPointer[UInt64, ImmUntrackedOrigin]
    """The nodes with further placements, in id order."""
    var more_at: ImmPointer[UInt64, ImmUntrackedOrigin]
    """Where each one's further placements start in `more`."""
    var more: ImmPointer[UInt32, ImmUntrackedOrigin]
    var more_count: Int
    var has_more: List[UInt64]
    """A bit per node id: set when the node has further placements."""
    var path_contig: List[Int]
    var contigs: List[String]

    def __init__(out self, path: String, sample: String) raises:
        var cache = CacheFile(path)
        var index = cache.sample_index(sample)
        if not index:
            raise Error(t"{path}: no reference sample {sample} in the graph's cache")
        var s = index.value()
        var damaged = String(t"{path}: a damaged graph cache; make it again with dromos prepare")
        var lengths = cache.section(LENGTHS)
        var first = cache.section(FIRST, s)
        var placements = cache.section(PLACEMENTS, s)
        var paths = cache.section(PATHS, s)
        var more_nodes = cache.section(MORE_NODES, s)
        var more_at = cache.section(MORE_AT, s)
        var more = cache.section(MORE, s)
        for sec in [lengths, first, placements, paths, more_nodes, more_at]:
            cache.advise(sec)
        self.size = cache.count
        if first.length != 4 * self.size:
            raise Error(damaged)
        self.node_len = cache.ptr(lengths).unsafe_bitcast[UInt32]()
        self.first = cache.ptr(first).unsafe_bitcast[UInt32]()
        self.placements = cache.ptr(placements).unsafe_bitcast[UInt64]()
        self.placement_count = placements.length // 8
        # The paths, their contigs and the contigs' names.
        var p = cache.ptr(paths)
        var n = u64_at(p, 0) if paths.length >= 8 else -1
        var names_at = 8 + 8 * (n + 1) + (4 * n + 7) // 8 * 8
        if n < 0 or names_at + 8 > paths.length:
            raise Error(damaged)
        self.path_first = List[Int](capacity=n + 1)
        self.path_contig = List[Int](capacity=n)
        for i in range(n + 1):
            self.path_first.append(u64_at(p, 8 + 8 * i))
        for i in range(n):
            self.path_contig.append(u32_at(p, 8 + 8 * (n + 1) + 4 * i))
        var c = u64_at(p, names_at)
        var at = names_at + 8
        self.contigs = List[String](capacity=c)
        for _ in range(c):
            var length = u64_at(p, at) if at + 8 <= paths.length else paths.length
            if at + 8 + length > paths.length:
                raise Error(damaged)
            self.contigs.append(
                String(StringSlice(unsafe_from_utf8=Span(unsafe_ptr=p.unsafe_offset(at + 8), length=length)))
            )
            at += 8 + length
        for i in range(n):
            if self.path_contig[i] >= c or self.path_first[i] > self.path_first[i + 1]:
                raise Error(damaged)
        if self.path_first[n] != self.placement_count:
            raise Error(damaged)
        # Nodes with further placements.
        self.more_count = more_nodes.length // 8
        if more_at.length != 8 * (self.more_count + 1):
            raise Error(damaged)
        self.more_nodes = cache.ptr(more_nodes).unsafe_bitcast[UInt64]()
        self.more_at = cache.ptr(more_at).unsafe_bitcast[UInt64]()
        self.more = cache.ptr(more).unsafe_bitcast[UInt32]()
        if Int(self.more_at[unsafe_offset=self.more_count]) * 4 > more.length:
            raise Error(damaged)
        self.has_more = List[UInt64](length=(self.size + 63) // 64, fill=0)
        for i in range(self.more_count):
            var node = Int(self.more_nodes[unsafe_offset=i])
            if node >= self.size:
                raise Error(damaged)
            self.has_more[node >> 6] |= UInt64(1) << UInt64(node & 63)
        self.sample = sample
        self.cache = cache^

    @always_inline
    def _len(self, node: Int) -> Int:
        """Node `node`'s entry in LENGTHS (NONE for an id the graph lacks); `node` is below `size`."""
        return Int(self.node_len[unsafe_offset=node])

    @always_inline
    def length(self, node: Int) -> Int:
        """Node `node`'s length; 0 for an id the graph lacks."""
        if node < 0 or node >= self.size:
            return 0
        var length = self._len(node)
        return 0 if length == NONE else length

    def has_letters(self) -> Bool:
        return False

    def letter(self, node: Int, reverse: Bool, i: Int) -> UInt8:
        """The map holds no letters (`has_letters` is false)."""
        return UInt8(ord("N"))

    @always_inline
    def _first(self, node: Int) -> Int:
        """Node `node`'s first placement; NONE for none (FIRST holds it plus one, 0 for none)."""
        return Int(self.first[unsafe_offset=node] - 1)

    @always_inline
    def _placement(self, k: Int) -> Placement:
        var v = Int(self.placements[unsafe_offset=k])
        var path = (v >> 1) & ((1 << PATH_BITS) - 1)
        return Placement(path, k - self.path_first[path], v >> START_SHIFT, v & 1 == 1)

    @always_inline
    def _more(self, node: Int) -> Tuple[Int, Int]:
        """Where node `node`'s further placements start in `more`, and how many there are (0 for none)."""
        if (self.has_more[node >> 6] >> UInt64(node & 63)) & 1 == 0:
            return (0, 0)
        var lo = 0
        var hi = self.more_count
        while lo < hi:
            var mid = (lo + hi) // 2
            if Int(self.more_nodes[unsafe_offset=mid]) < node:
                lo = mid + 1
            else:
                hi = mid
        var at = Int(self.more_at[unsafe_offset=lo])
        return (at, Int(self.more_at[unsafe_offset=lo + 1]) - at)

    @always_inline
    def _placement_of(self, node: Int, first: Int, more: Tuple[Int, Int], j: Int) -> Placement:
        """Placement `j` of node `node` (its first, then its further ones)."""
        if j == 0:
            return self._placement(first)
        return self._placement(Int(self.more[unsafe_offset=more[0] + j - 1]))

    def contig_paths(self, name: String) -> List[ContigPath]:
        """The paths of contig `name` (none if the map has no such contig), each with its steps' nodes: placements
        are kept in path order and a path's placements in step order, so each path is a range of them, and the nodes
        come from a scan of every node's placements."""
        var out = List[ContigPath]()
        var contig = -1
        for c in range(len(self.contigs)):
            if self.contigs[c] == name:
                contig = c
        for path in range(len(self.path_contig)):
            if self.path_contig[path] == contig:
                var lo = self.path_first[path]
                var hi = self.path_first[path + 1]
                out.append(ContigPath(path, lo, List[Int](length=hi - lo, fill=-1)))
        if len(out) == 0:
            return out^
        for node in range(self.size):
            var k = self._first(node)
            if k != NONE:
                _set_step(out, k, node)
        for i in range(self.more_count):
            var node = Int(self.more_nodes[unsafe_offset=i])
            for j in range(Int(self.more_at[unsafe_offset=i]), Int(self.more_at[unsafe_offset=i + 1])):
                _set_step(out, Int(self.more[unsafe_offset=j]), node)
        return out^

    def paths_by_contig(self, threads: Int) -> List[List[ContigPath]]:
        """Every contig's paths (by contig number), each with its steps' nodes, from one scan of every node's
        placements (`contig_paths` for all of them)."""
        var out = List[List[ContigPath]](length=len(self.contigs), fill=List[ContigPath]())
        # Each path's place in its contig's list.
        var slot = List[Int](length=len(self.path_contig), fill=0)
        for path in range(len(self.path_contig)):
            ref list = out[self.path_contig[path]]
            slot[path] = len(list)
            var lo = self.path_first[path]
            list.append(ContigPath(path, lo, List[Int](length=self.path_first[path + 1] - lo, fill=-1)))

        # Each placement is one node's: tasks over ranges of node ids write distinct steps.
        var task = _StepsOf(
            untracked(Pointer(to=self)),
            mut_untracked(Pointer(to=out)),
            imm_untracked(Pointer(to=slot)),
            ceildiv(self.size, _NODES_PER_TASK),
        )
        parallel_for(task, task.ranges + 1, threads)
        # The tasks reach `slot` through an untracked pointer: kept alive until here.
        _ = slot^
        return out^

    def set_steps(self, mut out: List[List[ContigPath]], slot: List[Int], range_index: Int, ranges: Int):
        """Range `range_index` of the node ids' first placements, or (the one past the ranges) the further ones."""
        if range_index < ranges:
            for node in range(range_index * _NODES_PER_TASK, min(self.size, (range_index + 1) * _NODES_PER_TASK)):
                var k = self._first(node)
                if k != NONE:
                    self._set_step_of(out, slot, k, node)
            return
        for i in range(self.more_count):
            var node = Int(self.more_nodes[unsafe_offset=i])
            for j in range(Int(self.more_at[unsafe_offset=i]), Int(self.more_at[unsafe_offset=i + 1])):
                self._set_step_of(out, slot, Int(self.more[unsafe_offset=j]), node)

    @always_inline
    def _set_step_of(self, mut out: List[List[ContigPath]], slot: List[Int], k: Int, node: Int):
        """Placement `k` is of `node`: its step's node, in the path holding it (the first whose end is past it)."""
        var lo = 0
        var hi = len(self.path_contig)
        while lo < hi:
            var mid = (lo + hi) // 2
            if self.path_first[mid + 1] <= k:
                lo = mid + 1
            else:
                hi = mid
        if lo >= len(self.path_contig):
            return
        ref cp = out[self.path_contig[lo]][slot[lo]]
        if k >= cp.first and k < cp.first + len(cp.nodes):
            cp.nodes[k - cp.first] = node

    def first_placement_of(self, path: Int) -> Int:
        """Path `path`'s first placement (placements are in path order)."""
        return self.path_first[path]

    @always_inline
    def placement(self, k: Int) -> Placement:
        return self._placement(k)

    def anchors(self, walk: List[Int], ts: Int, te: Int, mut out: List[Anchor]) -> Int:
        """The walk's nodes (each `node << 1 | reverse`) that overlap its aligned part [ts, te) and lie on the
        reference, in walk order; 0, or why there are none."""
        out.clear()
        if len(walk) == 0 or te <= ts:
            return NO_ANCHOR
        var off = 0
        for step in walk:
            var node = step >> 1
            if node < 0 or node >= self.size:
                return MISSING_NODE
            var length = self._len(node)
            if length == NONE:
                return MISSING_NODE
            var next = off + length
            if length != 0 and next > ts and off < te:
                var local_start = ts - off if ts > off else 0
                var local_end = te - off if te < next else length
                if local_start < local_end and self._first(node) != NONE:
                    out.append(Anchor(node, step & 1 == 1, local_start, local_end))
            off = next
            if off >= te:
                break
        return NO_ANCHOR if len(out) == 0 else 0

    @always_inline
    def _interval(self, a: Anchor, p: Placement) -> Tuple[Int, Int]:
        """The reference bases an anchor's aligned part covers at placement `p`."""
        if a.reverse == p.reverse:
            return (p.start + a.local_start, p.start + a.local_end)
        var length = self._len(a.node)
        return (p.start + (length - a.local_end), p.start + (length - a.local_start))

    def _try_from(self, anchors: List[Anchor], p0: Placement, mut out: Placed) -> Int:
        """The interval the anchors give starting from the first one's placement `p0`: every later anchor must have
        exactly one placement on the same path, in the same relative orientation and strictly further along it."""
        var a0 = anchors[0]
        var forward = a0.reverse == p0.reverse
        var prev = p0.step
        var span = self._interval(a0, p0)
        var lo = span[0]
        var hi = span[1]
        for i in range(1, len(anchors)):
            var a = anchors[i]
            var first = self._first(a.node)
            var more = self._more(a.node)
            var found = False
            var at = Placement(0, 0, 0, False)
            for j in range(1 + more[1]):
                var p = self._placement_of(a.node, first, more, j)
                if p.path != p0.path or (a.reverse == p.reverse) != forward:
                    continue
                if (forward and p.step <= prev) or (not forward and p.step >= prev):
                    continue
                if found:
                    return AMBIGUOUS
                found = True
                at = p
            if not found:
                return INCONSISTENT
            prev = at.step
            var s = self._interval(a, at)
            lo = min(lo, s[0])
            hi = max(hi, s[1])
        if lo < hi:
            out = Placed(self.path_contig[p0.path], lo, hi)
            return 0
        return INCONSISTENT

    def project(self, walk: List[Int], ts: Int, te: Int, mut anchors: List[Anchor], mut out: Placed) -> Int:
        """The reference interval of the alignment of walk `walk` (each step `node << 1 | reverse`) over its part
        [ts, te), into `out`; 0, or why it has none. Two starts giving different intervals leave it ambiguous."""
        var why = self.anchors(walk, ts, te, anchors)
        if why != 0:
            return why
        var a0 = anchors[0]
        var first = self._first(a0.node)
        var more = self._more(a0.node)
        var resolved = False
        why = INCONSISTENT
        for j in range(1 + more[1]):
            var c = Placed(0, 0, 0)
            var got = self._try_from(anchors, self._placement_of(a0.node, first, more, j), c)
            if got == 0:
                if not resolved:
                    resolved = True
                    out = c
                elif out != c:
                    return AMBIGUOUS
            elif got == AMBIGUOUS:
                why = AMBIGUOUS
        return 0 if resolved else why
