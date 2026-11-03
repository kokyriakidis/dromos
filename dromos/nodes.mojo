"""A graph's nodes as dromos reads them: each node's length and, optionally, its letters, indexed by node id.

The table is one file, little-endian, every section starting at a multiple of 8 bytes:

    "DROMOSN1"                  magic
    u64 count                   node ids 0 .. count - 1 (an id the graph lacks has length 0)
    u32 lengths[count]          then zero padding to a multiple of 8
    u64 has_letters             1 if the letters follow, else 0
    u64 offsets[count + 1]      node i's letters are bases[offsets[i] .. offsets[i + 1]], forward strand
    u8  bases[offsets[count]]

`dromos node-table` writes it from a GBZ graph (the file `pgcatalog node-table` writes).
"""

from dromos.bgzf import MappedFile

comptime MAGIC = "DROMOSN1"
"""The table's first eight bytes."""


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
    """A node table, mapped into memory: its pages are read as they are first touched, by whichever thread touches
    them (a table of letters is several GB; reading it whole first held every run up)."""

    var file: MappedFile
    var count: Int
    var lengths_at: Int
    var offsets_at: Int
    """The byte offset of `offsets`, or -1 when the table has no letters."""
    var bases_at: Int

    def __init__(out self, path: String) raises:
        self.file = MappedFile(path)
        var n = self.file.size
        var magic = MAGIC.as_bytes()
        var ok = n >= 24
        for i in range(min(8, n)):
            ok = ok and self.file.ptr[unsafe_offset=i] == magic[i]
        if not ok:
            raise Error(t"{path}: not a dromos node table")
        self.count = Int(_u64_at(self.file.ptr, 8))
        self.lengths_at = 16
        var at = self.lengths_at + 4 * self.count
        at += (8 - at % 8) % 8
        if at + 8 > n:
            raise Error(t"{path}: truncated node table")
        var has_letters = _u64_at(self.file.ptr, at) == 1
        at += 8
        self.offsets_at = -1
        self.bases_at = -1
        if has_letters:
            self.offsets_at = at
            self.bases_at = at + 8 * (self.count + 1)
            if self.bases_at > n or self.bases_at + Int(_u64_at(self.file.ptr, self.offsets_at + 8 * self.count)) > n:
                raise Error(t"{path}: truncated node letters")

    @always_inline
    def _u64(self, at: Int) -> UInt64:
        return self.file.ptr.unsafe_offset(at).unsafe_bitcast[UInt64]()[]

    @always_inline
    def length(self, node: Int) -> Int:
        """Node `node`'s length; 0 for an id the graph lacks."""
        if node < 0 or node >= self.count:
            return 0
        return Int(self.file.ptr.unsafe_offset(self.lengths_at).unsafe_bitcast[UInt32]()[unsafe_offset=node])

    def has_letters(self) -> Bool:
        return self.offsets_at >= 0

    @always_inline
    def bases(self, node: Int) -> ImmPointer[UInt8, ImmUntrackedOrigin]:
        """Node `node`'s letters on its forward strand (`length(node)` of them; the table must have letters)."""
        return self.file.ptr.unsafe_offset(self.bases_at + Int(self._u64(self.offsets_at + 8 * node)))

    @always_inline
    def letter(self, node: Int, reverse: Bool, i: Int) -> UInt8:
        """Letter `i` of node `node` read in orientation `reverse` (its reverse complement when true)."""
        var start = Int(self._u64(self.offsets_at + 8 * node))
        var length = Int(self._u64(self.offsets_at + 8 * (node + 1))) - start
        var bases = self.file.ptr.unsafe_offset(self.bases_at + start)
        if not reverse:
            return bases[unsafe_offset=i]
        return complement(bases[unsafe_offset=length - 1 - i])


@always_inline
def _u64_at(p: ImmPointer[UInt8, _], at: Int) -> UInt64:
    """The little-endian u64 at byte `at` of the table (a multiple of 8 in it)."""
    return p.unsafe_offset(at).unsafe_bitcast[UInt64]()[]


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
