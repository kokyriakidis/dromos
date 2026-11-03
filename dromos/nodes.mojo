"""A graph's nodes as dromos reads them: each node's length and its letters, indexed by node id, from the graph's
cache (cache.mojo: its NODES, LENGTHS, STARTS and LETTERS)."""

from std.math import iota

from dromos.cache import BLOCK, LENGTHS, LETTERS, NONE, STARTS, CacheFile
from dromos.linearize import upper


trait Nodes:
    """What GAF lines need of a graph's nodes: each node's length (0 for an id the graph lacks) and, for cs, its
    letters."""

    def length(self, node: Int) -> Int:
        ...

    def has_letters(self) -> Bool:
        ...

    def letter(self, node: Int, reverse: Bool, i: Int) -> UInt8:
        ...


struct NodeTable(Movable, Nodes):
    """The graph's node lengths and, with `letters`, letters, mapped from its cache: pages are read as they are first
    touched, by whichever thread touches them."""

    var cache: CacheFile
    var count: Int
    var lengths: ImmPointer[UInt32, ImmUntrackedOrigin]
    var starts: ImmPointer[UInt64, ImmUntrackedOrigin]
    var bases: ImmPointer[UInt8, ImmUntrackedOrigin]
    var wide: Bool
    """Letters as bytes, not 4-bit codes."""
    var table: SIMD[DType.uint8, 16]
    var with_letters: Bool
    var forward_pairs: List[UInt16]
    """For each byte of two 4-bit letters: the low one's letter then the high one's, in capitals."""
    var reverse_pairs: List[UInt16]
    """For each byte of two 4-bit letters: the complements (`comp`) of the high one's letter, then the low one's."""

    def __init__(out self, path: String, letters: Bool = True) raises:
        var cache = CacheFile(path)
        var lengths = cache.section(LENGTHS)
        var starts = cache.section(STARTS)
        var bases = cache.section(LETTERS)
        cache.advise(lengths)
        if letters:
            cache.advise(starts)
            cache.advise(bases)
        self.count = cache.count
        self.lengths = cache.ptr(lengths).unsafe_bitcast[UInt32]()
        self.starts = cache.ptr(starts).unsafe_bitcast[UInt64]()
        self.bases = cache.ptr(bases)
        self.wide = cache.width == 8
        self.table = cache.table
        self.with_letters = letters
        self.forward_pairs = List[UInt16](capacity=256)
        self.reverse_pairs = List[UInt16](capacity=256)
        for byte in range(256):
            var low = self.table[byte & 15]
            var high = self.table[byte >> 4]
            self.forward_pairs.append(UInt16(upper(low)) | UInt16(upper(high)) << 8)
            self.reverse_pairs.append(UInt16(comp(high)) | UInt16(comp(low)) << 8)
        self.cache = cache^

    @always_inline
    def length(self, node: Int) -> Int:
        """Node `node`'s length; 0 for an id the graph lacks."""
        if node < 0 or node >= self.count:
            return 0
        var n = Int(self.lengths[unsafe_offset=node])
        return 0 if n == NONE else n

    def has_letters(self) -> Bool:
        return self.with_letters

    @always_inline
    def start(self, node: Int) -> Int:
        """Where node `node`'s letters start (an index for `at`): its block's start and the lengths before it in the
        block."""
        var block = node - node % BLOCK
        var v = self.lengths.unsafe_offset(block).unsafe_load[width=BLOCK, alignment=4 * BLOCK]()
        var lane = iota[DType.uint32, BLOCK]()
        var keep = lane.lt(UInt32(node - block)) & v.ne(UInt32(NONE))
        var before = keep.select(v, 0).cast[DType.uint64]().reduce_add()
        return Int(self.starts[unsafe_offset=node // BLOCK] + before)

    @always_inline
    def at(self, index: Int) -> UInt8:
        """The letter at `index`."""
        if self.wide:
            return self.bases[unsafe_offset=index]
        var code = (self.bases[unsafe_offset=index >> 1] >> UInt8((index & 1) << 2)) & 15
        return self.table[Int(code)]

    def copy_forward(self, dst: MutPointer[UInt8, _], index: Int, count: Int):
        """Letters `index ..< index + count` (as `at` numbers them) at `dst`, in capitals: two a lookup."""
        if self.wide:
            for t in range(count):
                dst[unsafe_offset=t] = upper(self.bases[unsafe_offset=index + t])
            return
        var t = 0
        if index & 1 == 1 and count > 0:
            dst[unsafe_offset=0] = upper(self.at(index))
            t = 1
        var src = self.bases.unsafe_offset((index + t) >> 1)
        var pairs = (count - t) >> 1
        var out = dst.unsafe_offset(t)
        for k in range(pairs):
            out.unsafe_offset(2 * k).unsafe_bitcast[UInt16]()[] = self.forward_pairs[Int(src[unsafe_offset=k])]
        t += 2 * pairs
        if t < count:
            dst[unsafe_offset=t] = upper(self.at(index + t))

    def copy_reverse(self, dst: MutPointer[UInt8, _], last: Int, count: Int):
        """The complements (`comp`) of letters `last`, `last - 1`, ... (`count` of them) at `dst`: a stretch read on
        the reverse strand, two a lookup."""
        if self.wide:
            for t in range(count):
                dst[unsafe_offset=t] = comp(self.bases[unsafe_offset=last - t])
            return
        var t = 0
        if last & 1 == 0 and count > 0:
            dst[unsafe_offset=0] = comp(self.at(last))
            t = 1
        # Letter last - t is a byte's high one: it and the low one before it, a byte at a time, downwards.
        var top = (last - t) >> 1
        var pairs = (count - t) >> 1
        var out = dst.unsafe_offset(t)
        for k in range(pairs):
            out.unsafe_offset(2 * k).unsafe_bitcast[UInt16]()[] = self.reverse_pairs[
                Int(self.bases[unsafe_offset=top - k])
            ]
        t += 2 * pairs
        if t < count:
            dst[unsafe_offset=t] = comp(self.at(last - t))

    @always_inline
    def letter(self, node: Int, reverse: Bool, i: Int) -> UInt8:
        """Letter `i` of node `node` read in orientation `reverse` (its reverse complement when true); N past the
        node's end or for a node the graph lacks (an alignment to another graph)."""
        if i < 0 or i >= self.length(node):
            return UInt8(ord("N"))
        var s = self.start(node)
        if not reverse:
            return self.at(s + i)
        return complement(self.at(s + self.length(node) - 1 - i))


@always_inline
def complement(c: UInt8) -> UInt8:
    """A base's complement, case kept; any other letter as it is."""
    if c == UInt8(ord("A")):
        return UInt8(ord("T"))
    if c == UInt8(ord("C")):
        return UInt8(ord("G"))
    if c == UInt8(ord("G")):
        return UInt8(ord("C"))
    if c == UInt8(ord("T")):
        return UInt8(ord("A"))
    if c == UInt8(ord("a")):
        return UInt8(ord("t"))
    if c == UInt8(ord("c")):
        return UInt8(ord("g"))
    if c == UInt8(ord("g")):
        return UInt8(ord("c"))
    if c == UInt8(ord("t")):
        return UInt8(ord("a"))
    return c


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
