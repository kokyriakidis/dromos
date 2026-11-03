"""Alignments placed on the reference: pgcat's reference map, read whole, and its placement rule (pgcat's gaf.rs,
`RefMap::project`, itself pggaf's walk_ref_interval), ported line for line so dromos places every alignment where
`pgcat index-gaf` does.

The map (`pgcat ref-map` writes it once per graph) is little-endian:

    "PGCATRM1"                         magic
    u64 n, n bytes                     the reference sample
    u64 size                           node ids 0 .. size - 1
    u32 node_len[size]                 each node's length; NONE for an id the graph lacks
    u32 first[size]                    the node's first placement; NONE for a node off the reference
    u64 count, count x 17 bytes        placements: path u32, step u32, start u64, reverse u8
    u64 k, then k times:               nodes the reference visits more than once, by node id:
        u64 node, u64 m, u32 more[m]       the placements after its first
    u64 p, u32 path_contig[p]          each reference path's contig
    u64 c, then c times u64 n, n bytes the contigs' names
"""

from dromos.nodes import Nodes, read_file

comptime MAGIC = "PGCATRM1"
comptime NONE = 0xFFFF_FFFF
"""A node_len or first entry with no value."""
comptime PLACEMENT_BYTES = 17

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


@always_inline
def _u32(p: ImmPointer[UInt8, _], at: Int) -> Int:
    return (
        Int(p[unsafe_offset=at])
        | Int(p[unsafe_offset=at + 1]) << 8
        | Int(p[unsafe_offset=at + 2]) << 16
        | Int(p[unsafe_offset=at + 3]) << 24
    )


@always_inline
def _u64(p: ImmPointer[UInt8, _], at: Int) -> Int:
    return _u32(p, at) | _u32(p, at + 4) << 32


struct RefMap(Movable, Nodes):
    """A pgcat reference map read whole into memory: node lengths and placements are read in place."""

    var data: List[UInt8]
    var sample: String
    var size: Int
    var node_len_at: Int
    var first_at: Int
    var placements_at: Int
    var more_nodes: List[Int]
    """The nodes with further placements, in id order."""
    var more_at: List[Int]
    """Where each one's list of further placements starts in `data`."""
    var more_count: List[Int]
    """How many further placements each one has."""
    var has_more: List[UInt64]
    """A bit per node id: set when the node has further placements."""
    var path_contig: List[Int]
    var contigs: List[String]

    def __init__(out self, path: String) raises:
        self.data = read_file(path)
        var n = len(self.data)
        var p = self.data.unsafe_ptr()
        var bad = String(t"{path}: not a pgcat ref map (pgcat ref-map writes one)")
        var ok = n >= 16
        var magic = MAGIC.as_bytes()
        for i in range(8):
            ok = ok and self.data[i] == magic[i]
        if not ok:
            raise Error(bad)
        var at = 8
        var length = _u64(p, at)
        at += 8
        if length > n - at:
            raise Error(bad)
        self.sample = String(StringSlice(unsafe_from_utf8=Span(self.data)[at : at + length]))
        at += length
        if at + 8 > n:
            raise Error(bad)
        self.size = _u64(p, at)
        at += 8
        self.node_len_at = at
        self.first_at = at + 4 * self.size
        at = self.first_at + 4 * self.size
        if at + 8 > n:
            raise Error(bad)
        var count = _u64(p, at)
        at += 8
        self.placements_at = at
        at += PLACEMENT_BYTES * count
        if at + 8 > n:
            raise Error(bad)
        var k = _u64(p, at)
        at += 8
        self.more_nodes = List[Int](capacity=k)
        self.more_at = List[Int](capacity=k)
        self.more_count = List[Int](capacity=k)
        self.has_more = List[UInt64](length=(self.size + 63) // 64, fill=0)
        for _ in range(k):
            if at + 16 > n:
                raise Error(bad)
            var node = _u64(p, at)
            var m = _u64(p, at + 8)
            at += 16
            if node >= self.size or (len(self.more_nodes) > 0 and node <= self.more_nodes[len(self.more_nodes) - 1]):
                raise Error(bad)
            self.more_nodes.append(node)
            self.more_at.append(at)
            self.more_count.append(m)
            self.has_more[node >> 6] |= UInt64(1) << UInt64(node & 63)
            at += 4 * m
        if at + 8 > n:
            raise Error(bad)
        var paths = _u64(p, at)
        at += 8
        if at + 4 * paths + 8 > n:
            raise Error(bad)
        self.path_contig = List[Int](capacity=paths)
        for i in range(paths):
            self.path_contig.append(_u32(p, at + 4 * i))
        at += 4 * paths
        var c = _u64(p, at)
        at += 8
        self.contigs = List[String](capacity=c)
        for _ in range(c):
            if at + 8 > n:
                raise Error(bad)
            var name_length = _u64(p, at)
            at += 8
            if name_length > n - at:
                raise Error(bad)
            self.contigs.append(String(StringSlice(unsafe_from_utf8=Span(self.data)[at : at + name_length])))
            at += name_length
        for i in range(paths):
            if self.path_contig[i] >= c:
                raise Error(bad)

    @always_inline
    def _ptr(self) -> ImmPointer[UInt8, origin_of(self.data)]:
        return self.data.unsafe_ptr().as_imm()

    @always_inline
    def length(self, node: Int) -> Int:
        """Node `node`'s length; 0 for an id the graph lacks."""
        if node < 0 or node >= self.size:
            return 0
        var length = _u32(self._ptr(), self.node_len_at + 4 * node)
        return 0 if length == NONE else length

    def has_letters(self) -> Bool:
        return False

    def letter(self, node: Int, reverse: Bool, i: Int) -> UInt8:
        """The map holds no letters (`has_letters` is false)."""
        return UInt8(ord("N"))

    @always_inline
    def _first(self, node: Int) -> Int:
        return _u32(self._ptr(), self.first_at + 4 * node)

    @always_inline
    def _placement(self, k: Int) -> Placement:
        var p = self._ptr()
        var at = self.placements_at + PLACEMENT_BYTES * k
        return Placement(_u32(p, at), _u32(p, at + 4), _u64(p, at + 8), p[unsafe_offset=at + 16] != 0)

    @always_inline
    def _more(self, node: Int) -> Tuple[Int, Int]:
        """Where node `node`'s further placements are listed in `data`, and how many there are (0 for none)."""
        if (self.has_more[node >> 6] >> UInt64(node & 63)) & 1 == 0:
            return (0, 0)
        var lo = 0
        var hi = len(self.more_nodes)
        while lo < hi:
            var mid = (lo + hi) // 2
            if self.more_nodes[mid] < node:
                lo = mid + 1
            else:
                hi = mid
        return (self.more_at[lo], self.more_count[lo])

    @always_inline
    def _placement_of(self, node: Int, first: Int, more: Tuple[Int, Int], j: Int) -> Placement:
        """Placement `j` of node `node` (its first, then its further ones)."""
        if j == 0:
            return self._placement(first)
        return self._placement(_u32(self._ptr(), more[0] + 4 * (j - 1)))

    def anchors(self, walk: List[Int], ts: Int, te: Int, mut out: List[Anchor]) -> Int:
        """The walk's nodes (each `node << 1 | reverse`) that overlap its aligned part [ts, te) and lie on the
        reference, in walk order; 0, or why there are none."""
        out.clear()
        if len(walk) == 0 or te <= ts:
            return NO_ANCHOR
        var p = self._ptr()
        var off = 0
        for step in walk:
            var node = step >> 1
            if node < 0 or node >= self.size:
                return MISSING_NODE
            var length = _u32(p, self.node_len_at + 4 * node)
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
        var length = _u32(self._ptr(), self.node_len_at + 4 * a.node)
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
