"""A reference contig as `dromos surject` projects reads onto it (pgcatalog's surject.rs `load_contig`, over its
reference.rs `ReferencePath`), built from the graph's cache (its reference map and node table) instead of the graph: the contig's
letters, and each node the contig's path visits once, with where and whether reversed.

The contig's paths (fragments) are laid end to end in the order of their offsets (then of their paths), each step at
its placement's start; bases no fragment covers are N, and the contig is as long as its farthest fragment end."""

from dromos.linearize import upper
from dromos.nodes import NodeTable, comp
from dromos.parallel import ParallelTask, parallel_for
from dromos.refmap import ContigPath, RefMap


struct Contig(Movable):
    var chrom: String
    """`<sample>#0#<contig>`, the SAM's contig name unless another is given."""
    var seq: List[UInt8]
    var base: Int
    """The smallest node id in `place`."""
    var place: List[Int]
    """For node `base + i`: `start << 1 | reverse` if the contig visits it once, else -1."""

    def __init__(out self):
        """No contig: a place for one to be made in."""
        self.chrom = String()
        self.seq = List[UInt8]()
        self.base = 0
        self.place = List[Int]()

    def __init__(out self, map: RefMap, nodes: NodeTable, contig: String) raises:
        self = Contig(map, nodes, contig, map.contig_paths(contig))

    def __init__(out self, map: RefMap, nodes: NodeTable, contig: String, paths: List[ContigPath]) raises:
        """Contig `contig` from its paths (`RefMap.contig_paths`)."""
        self.chrom = String(t"{map.sample}#0#{contig}")
        if len(paths) == 0:
            raise Error(t"no reference path {self.chrom} in the graph")
        # Fragments in offset order (step 0's start), then path order.
        var order = List[Int]()
        for i in range(len(paths)):
            if len(paths[i].nodes) > 0:
                order.append(i)
        for a in range(1, len(order)):
            var k = a
            while k > 0:
                var x = map.placement(paths[order[k]].first).start
                var y = map.placement(paths[order[k - 1]].first).start
                if x > y or (x == y and paths[order[k]].path > paths[order[k - 1]].path):
                    break
                order.swap_elements(k, k - 1)
                k -= 1
        var length = 0
        var lo = map.size
        var hi = -1
        for i in order:
            ref cp = paths[i]
            for s in range(len(cp.nodes)):
                var node = cp.nodes[s]
                if node < 0:
                    raise Error(t"{self.chrom}: a step without a node in the reference map")
                lo = min(lo, node)
                hi = max(hi, node)
                length = max(length, map.placement(cp.first + s).start + map.length(node))
        self.seq = List[UInt8](length=length, fill=UInt8(ord("N")))
        self.base = lo
        var visits = List[Int](length=hi - lo + 1, fill=0)
        self.place = List[Int](length=hi - lo + 1, fill=-1)
        for i in order:
            ref cp = paths[i]
            for s in range(len(cp.nodes)):
                var node = cp.nodes[s]
                var p = map.placement(cp.first + s)
                var n = map.length(node)
                var letters = nodes.start(node)
                if p.reverse:
                    for k in range(n):
                        self.seq[p.start + k] = comp(nodes.at(letters + n - 1 - k))
                else:
                    for k in range(n):
                        self.seq[p.start + k] = nodes.at(letters + k)
                visits[node - lo] += 1
                self.place[node - lo] = p.start << 1 | Int(p.reverse)
        for k in range(len(visits)):
            if visits[k] != 1:
                self.place[k] = -1

    @always_inline
    def placed(self, node: Int) -> Int:
        """`start << 1 | reverse` of a node the contig visits once, else -1."""
        var k = node - self.base
        if k < 0 or k >= len(self.place):
            return -1
        return self.place[k]


@fieldwise_init
struct _MakeContigs[map_origin: ImmOrigin, nodes_origin: ImmOrigin](ParallelTask):
    """Contig `i` of `names` made into `out[i]` from its paths, a contig a task."""

    var map: ImmPointer[RefMap, Self.map_origin]
    var nodes: ImmPointer[NodeTable, Self.nodes_origin]
    var names: ImmPointer[String, ImmUntrackedOrigin]
    var paths: ImmPointer[List[ContigPath], ImmUntrackedOrigin]
    var out: MutPointer[Contig, MutUntrackedOrigin]
    var failed: MutPointer[Int, MutUntrackedOrigin]

    def run(self, index: Int):
        try:
            self.out[unsafe_offset=index] = Contig(
                self.map[], self.nodes[], self.names[unsafe_offset=index], self.paths[unsafe_offset=index]
            )
        except:
            self.failed[] = index + 1


struct Reference(Movable):
    """The contigs reads are surjected to (in reference order: a contig's number is its place here) and, for any
    node, the contigs it may lie on: the node id ranges of the contigs' once-visited nodes cut into segments, each
    with the contigs whose range covers it."""

    var contigs: List[Contig]
    var names: List[String]
    """Each contig's name without its sample (chr20)."""
    var seg_start: List[Int]
    """Each segment's first node id, then the end of the last."""
    var seg_first: List[Int]
    """Where each segment's contigs start in `seg_contigs`, then the end of the last."""
    var seg_contigs: List[Int]
    var lo: Int
    var hi: Int
    """The node ids of all segments: [lo, hi)."""

    def __init__(out self, map: RefMap, nodes: NodeTable, names: List[String], threads: Int) raises:
        """The contigs `names` (in that order), made at once on `threads` threads."""
        var paths = List[List[ContigPath]]()
        if len(names) == 1:
            # One contig: its paths alone, from one scan.
            paths.append(map.contig_paths(names[0]))
        else:
            var by_contig = map.paths_by_contig(threads)
            for name in names:
                var c = -1
                for k in range(len(map.contigs)):
                    if map.contigs[k] == name:
                        c = k
                if c < 0:
                    raise Error(t"no reference path {map.sample}#0#{name} in the graph")
                var taken = List[ContigPath]()
                swap(taken, by_contig[c])
                paths.append(taken^)
        self.contigs = List[Contig]()
        for _ in range(len(names)):
            self.contigs.append(Contig())
        self.names = names.copy()
        var failed = List[Int](length=1, fill=0)
        var task = _MakeContigs(
            Pointer(to=map),
            Pointer(to=nodes),
            self.names.unsafe_ptr().as_imm().unsafe_origin_cast[ImmUntrackedOrigin](),
            paths.unsafe_ptr().as_imm().unsafe_origin_cast[ImmUntrackedOrigin](),
            self.contigs.unsafe_ptr().unsafe_origin_cast[MutUntrackedOrigin](),
            failed.unsafe_ptr().unsafe_origin_cast[MutUntrackedOrigin](),
        )
        parallel_for(task, len(names), threads)
        # The paths and the flag are read by the tasks through untracked pointers: kept alive until here.
        _ = paths^
        if failed[0] > 0:
            # Made again on this thread, for its error.
            _ = Contig(map, nodes, names[failed[0] - 1])
            raise Error(t"{names[failed[0] - 1]}: its contig could not be made")
        # Segments: cut at every range's ends.
        var cuts = List[Int]()
        for k in range(len(self.contigs)):
            ref c = self.contigs[k]
            if len(c.place) > 0:
                cuts.append(c.base)
                cuts.append(c.base + len(c.place))
        sort(Span(cuts))
        self.seg_start = List[Int]()
        for x in cuts:
            if len(self.seg_start) == 0 or self.seg_start[len(self.seg_start) - 1] != x:
                self.seg_start.append(x)
        self.seg_first = List[Int]()
        self.seg_contigs = List[Int]()
        for s in range(len(self.seg_start) - 1):
            self.seg_first.append(len(self.seg_contigs))
            for k in range(len(self.contigs)):
                ref c = self.contigs[k]
                if len(c.place) > 0 and c.base <= self.seg_start[s] and c.base + len(c.place) > self.seg_start[s]:
                    self.seg_contigs.append(k)
        self.seg_first.append(len(self.seg_contigs))
        self.lo = self.seg_start[0] if len(self.seg_start) > 0 else 0
        self.hi = self.seg_start[len(self.seg_start) - 1] if len(self.seg_start) > 0 else 0

    @always_inline
    def _segment(self, node: Int) -> Int:
        """The segment holding node id `node` (within [lo, hi))."""
        var lo = 0
        var hi = len(self.seg_start) - 1
        while hi - lo > 1:
            var mid = (lo + hi) // 2
            if self.seg_start[mid] <= node:
                lo = mid
            else:
                hi = mid
        return lo

    @always_inline
    def placed_anywhere(self, node: Int) -> Bool:
        """Whether some contig visits node `node` once."""
        if node < self.lo or node >= self.hi:
            return False
        var s = self._segment(node)
        for i in range(self.seg_first[s], self.seg_first[s + 1]):
            if self.contigs[self.seg_contigs[i]].placed(node) >= 0:
                return True
        return False

    def add_contigs_of(self, node: Int, mut out: List[Int]):
        """Adds to `out` (kept sorted, without repeats) the contigs that visit node `node` once."""
        if node < self.lo or node >= self.hi:
            return
        var s = self._segment(node)
        for i in range(self.seg_first[s], self.seg_first[s + 1]):
            var k = self.seg_contigs[i]
            if self.contigs[k].placed(node) < 0:
                continue
            var at = len(out)
            while at > 0 and out[at - 1] > k:
                at -= 1
            if at > 0 and out[at - 1] == k:
                continue
            out.insert(at, k)
