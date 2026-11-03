"""A reference contig as `dromos surject` projects reads onto it (pgcatalog's surject.rs `load_contig`, over its
reference.rs `ReferencePath`), built from the graph's cache (its reference map and node table) instead of the graph: the contig's
letters, and each node the contig's path visits once, with where and whether reversed.

The contig's paths (fragments) are laid end to end in the order of their offsets (then of their paths), each step at
its placement's start; bases no fragment covers are N, and the contig is as long as its farthest fragment end."""

from dromos.linearize import upper
from dromos.nodes import NodeTable
from dromos.refmap import RefMap


comptime _A = UInt8(ord("A"))
comptime _C = UInt8(ord("C"))
comptime _G = UInt8(ord("G"))
comptime _T = UInt8(ord("T"))
comptime _N = UInt8(ord("N"))


@always_inline
def comp(b: UInt8) -> UInt8:
    """pgcatalog's complement: of A, C, G or T in either case, in capitals; N for any other letter."""
    var u = upper(b)
    if u == _A:
        return _T
    if u == _C:
        return _G
    if u == _G:
        return _C
    if u == _T:
        return _A
    return _N


struct Contig(Movable):
    var chrom: String
    """`<sample>#0#<contig>`, the SAM's contig name unless another is given."""
    var seq: List[UInt8]
    var base: Int
    """The smallest node id in `place`."""
    var place: List[Int]
    """For node `base + i`: `start << 1 | reverse` if the contig visits it once, else -1."""

    def __init__(out self, map: RefMap, nodes: NodeTable, contig: String) raises:
        var paths = map.contig_paths(contig)
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
        self.seq = List[UInt8](length=length, fill=_N)
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
