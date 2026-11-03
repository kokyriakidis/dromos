"""Snarls from a giraffe distance index (`.dist`, libbdsg's SnarlDistanceIndex), read in place: the index is mapped,
never loaded, and a contig's snarls are found by walking its connected component's snarl tree.

The file: a 4-byte magic (67 A1 7C C6), then libbdsg's arena image (little-endian), whose first allocation is one
packed int vector (`length` entries of `width` bits, least significant bit first across 64-bit words). Entry 0 starts
the root record: [tag, 1023 ^ version, components N, node id slots M, min node id, depth], the N root children's
record offsets, then per node id its record's offset and its place there. Records start with a tag
(`type << 9 | flags`). Version 3 and 4 share one type numbering (4 changes only how vg fills the index); version 5
renumbers the types and is refused.

A chain record lists its children from entry +15: trivial-snarl records (runs of up to 255 chain nodes, each node
`id << 1 | reverse`) and snarl records between them, each record's size in the slot before it. A snarl stores no
boundary nodes: its start is the last node of the trivial-snarl record before it, its end the first node of the
record after it (the chain's first record when a looping chain ends with the snarl). A snarl's children (chains,
single nodes) are listed at the offset in its entry +7; a simple snarl's children are single nodes, held inline.

Orientation follows vg's SnarlManager::regularize, as `vg snarls` writes it: per chain, every consecutive pair of chain
nodes is one of vg's snarls (the edge-only ones too, which the index folds into its trivial-snarl records); when fewer
than half of them ascend by node id the chain, and every snarl in it, is flipped (start and end swapped and reversed).
Parents are named in their own chain's orientation. Checked against vg snarls on HPRC v2.1 (51,186,896 snarls, every
one the same).
"""

from dromos.bgzf import MappedFile
from dromos.snarls import NO_PARENT, Snarl

comptime MAGIC = SIMD[DType.uint8, 4](0x67, 0xA1, 0x7C, 0xC6)
comptime VERSION_SENTINEL = 1023
comptime READ_VERSIONS = "3 and 4"
"""The versions this reader takes."""
comptime ORIENTATION = (
    "vg SnarlManager::regularize: per chain, flipped when fewer than half of its consecutive node pairs ascend by id"
)

# Record types of versions 3 and 4 (tag >> 9).
comptime ROOT = 1
comptime NODE = 2
comptime DISTANCED_NODE = 3
comptime TRIVIAL_SNARL = 4
comptime DISTANCED_TRIVIAL_SNARL = 5
comptime SIMPLE_SNARL = 6
comptime DISTANCED_SIMPLE_SNARL = 7
comptime SNARL = 8
comptime DISTANCED_SNARL = 9
comptime OVERSIZED_SNARL = 10
comptime ROOT_SNARL = 11
comptime DISTANCED_ROOT_SNARL = 12
comptime CHAIN = 13
comptime DISTANCED_CHAIN = 14
comptime MULTICOMPONENT_CHAIN = 15


@always_inline
def _is_trivial(t: Int) -> Bool:
    return t == TRIVIAL_SNARL or t == DISTANCED_TRIVIAL_SNARL


@always_inline
def _is_simple(t: Int) -> Bool:
    return t == SIMPLE_SNARL or t == DISTANCED_SIMPLE_SNARL


@always_inline
def _is_snarl(t: Int) -> Bool:
    return t == SNARL or t == DISTANCED_SNARL or t == OVERSIZED_SNARL


@always_inline
def _is_chain(t: Int) -> Bool:
    return t == CHAIN or t == DISTANCED_CHAIN or t == MULTICOMPONENT_CHAIN


@always_inline
def _is_node(t: Int) -> Bool:
    return t == NODE or t == DISTANCED_NODE


@always_inline
def _is_root_snarl(t: Int) -> Bool:
    return t == ROOT_SNARL or t == DISTANCED_ROOT_SNARL


def _u64(p: ImmPointer[UInt8, _], at: Int) -> Int:
    var v = 0
    for k in range(8):
        v |= Int(p[unsafe_offset=at + k]) << (8 * k)
    return v


@always_inline
def _entry(words: ImmPointer[UInt64, _], width: Int, mask: UInt64, i: Int) -> Int:
    """Entry `i` of a packed int vector: bits [i * width, (i + 1) * width), least significant first."""
    var b = i * width
    var k = b >> 6
    var o = b & 63
    var v = words[unsafe_offset=k] >> UInt64(o)
    if o + width > 64:
        v |= words[unsafe_offset=k + 1] << UInt64(64 - o)
    return Int(v & mask)


@fieldwise_init
struct _Header(TrivialRegisterPassable):
    var words: ImmPointer[UInt64, ImmUntrackedOrigin]
    var length: Int
    var width: Int
    var mask: UInt64
    var version: Int
    var components: Int
    var id_slots: Int
    var min_id: Int


def _header(path: String, file: MappedFile) raises -> _Header:
    """The index's int vector found and its root record read (`file` borrowed for the whole call, so it stays mapped
    on every path)."""
    var p = file.ptr
    var n = file.size
    var bad = String(t"{path}: not a distance index (vg/libbdsg SnarlDistanceIndex)")
    if n < 64:
        raise Error(bad)
    for i in range(4):
        if p[unsafe_offset=i] != MAGIC[i]:
            raise Error(bad)
    # The arena header after the magic: its first allocation (a 16-byte pointer at byte 36: an offset from where it
    # lies) is the int vector's object.
    var rel = _u64(p, 36)
    if rel == Int.MAX or rel <= 0 or 36 + rel + 48 > n:
        raise Error(t"{path}: a damaged distance index (its arena header)")
    var obj = 36 + rel
    var length = _u64(p, obj)
    var width = _u64(p, obj + 8)
    var words = _u64(p, obj + 24)
    var data_rel = _u64(p, obj + 40)
    if data_rel == Int.MAX:
        raise Error(t"{path}: a damaged distance index (no data)")
    var data = obj + 40 + data_rel
    if (
        width < 1
        or width > 64
        or data < 0
        or data % 8 != 0
        or words < 0
        or data + 8 * words > n
        or words * 64 < length * width
        or length < 6
    ):
        raise Error(t"{path}: a damaged distance index (its int vector)")
    var w = ImmPointer[UInt8, ImmUntrackedOrigin](unsafe_from_address=Int(p) + data).unsafe_bitcast[UInt64]()
    var mask = UInt64.MAX if width == 64 else (UInt64(1) << UInt64(width)) - 1
    if _entry(w, width, mask, 0) >> 9 != ROOT:
        raise Error(t"{path}: a damaged distance index (no root record)")
    var version = _entry(w, width, mask, 1) ^ VERSION_SENTINEL
    if version < 3 or version > 64:
        raise Error(
            t"{path}: a distance index from before libbdsg's versioned format (2025); dromos reads versions"
            t" {READ_VERSIONS}. Build the catalog from vg snarls instead: dromos catalog --snarls-cache DIR (vg chunk"
            t" and vg snarls run for the contigs not cached there)"
        )
    if version > 4:
        raise Error(
            t"{path}: a distance index of version {version}; dromos reads versions {READ_VERSIONS} (version 5"
            t" renumbers its records). Build the catalog from vg snarls instead: dromos catalog --snarls-cache DIR (vg"
            t" chunk and vg snarls run for the contigs not cached there)"
        )
    var components = _entry(w, width, mask, 2)
    var id_slots = _entry(w, width, mask, 3)
    var min_id = _entry(w, width, mask, 4)
    if 6 + components + 2 * id_slots > length:
        raise Error(t"{path}: a damaged distance index (its root record)")
    return _Header(w, length, width, mask, version, components, id_slots, min_id)


struct DistIndex(Movable):
    """A `.dist` mapped into memory, its root record read."""

    var path: String
    var file: MappedFile
    var words: ImmPointer[UInt64, ImmUntrackedOrigin]
    var length: Int
    """Entries in the int vector."""
    var width: Int
    var mask: UInt64
    var version: Int
    var components: Int
    var id_slots: Int
    var min_id: Int

    def __init__(out self, path: String) raises:
        self.path = path
        var file = MappedFile(path, ahead=False)
        var h = _header(path, file)
        self.file = file^
        self.words = h.words
        self.length = h.length
        self.width = h.width
        self.mask = h.mask
        self.version = h.version
        self.components = h.components
        self.id_slots = h.id_slots
        self.min_id = h.min_id

    @always_inline
    def at(self, i: Int) -> Int:
        """Entry `i` of the int vector (callers keep `i` below `length`)."""
        return _entry(self.words, self.width, self.mask, i)

    def entry(self, i: Int) raises -> Int:
        if i < 0 or i >= self.length:
            raise Error(t"{self.path}: a damaged distance index (an offset past its end)")
        return self.at(i)

    @always_inline
    def kind(self, record: Int) raises -> Int:
        return self.entry(record) >> 9

    def max_id(self) -> Int:
        return self.min_id + self.id_slots - 1

    def node_record(self, id: Int) -> Int:
        """The record holding node `id` (0: the index has no such node)."""
        if id < self.min_id or id > self.max_id():
            return 0
        return self.at(6 + self.components + 2 * (id - self.min_id))

    def node_count(self) -> Int:
        """How many node ids the index holds."""
        var n = 0
        var base = 6 + self.components
        for i in range(self.id_slots):
            if self.at(base + 2 * i) != 0:
                n += 1
        return n

    def parent(self, record: Int) raises -> Int:
        """A record's parent (0: the root)."""
        var t = self.kind(record)
        if _is_trivial(t):
            return self.entry(record + 1)
        if _is_node(t) or _is_simple(t) or _is_snarl(t) or _is_chain(t):
            return self.entry(record + 2)
        if _is_root_snarl(t):
            return 0
        raise Error(t"{self.path}: a damaged distance index (a record of type {t})")

    def top(self, node: Int) raises -> Int:
        """The root child (a chain, a node or a root snarl) whose snarl tree holds node `node`."""
        var r = self.node_record(node)
        if r == 0:
            raise Error(t"{self.path}: node {node} is not in the distance index")
        var steps = 0
        while True:
            var up = self.parent(r)
            if up == 0:
                return r
            r = up
            steps += 1
            if steps > 1 << 20:
                raise Error(t"{self.path}: a damaged distance index (a cycle of parents)")

    def snarls(self, top: Int, mut out: List[Snarl]) raises:
        """The snarls under root child `top`, oriented as vg snarls writes them (the module's notes), into `out`."""
        var t = self.kind(top)
        if _is_root_snarl(t):
            raise Error(
                t"{self.path}: a root snarl (a component with no chain through it), which dromos does not read; build"
                t" the catalog from vg snarls instead (--snarls-cache DIR)"
            )
        if _is_chain(t):
            self._chain(top, NO_PARENT, NO_PARENT, out)
        elif not _is_node(t):
            raise Error(t"{self.path}: a damaged distance index (a root child of type {t})")

    def chain_words(self, chain: Int, mut words: List[UInt64], mut snarl_at: List[Int]) raises:
        """A chain's nodes in order (`id << 1 | reverse`), and before which of them each of its snarls lies (index
        into `words`, with the snarl's record): pairs (record, k) for a snarl between words[k - 1] and words[k] (k
        == len(words) when a looping chain ends with it)."""
        words.clear()
        snarl_at.clear()
        var first = chain + 15
        var last = self.entry(chain + 8) >> 2
        var p = first
        var steps = 0
        while True:
            var t = self.kind(p)
            if _is_trivial(t):
                var k = self.entry(p + 2)
                var stride = 2 if t == DISTANCED_TRIVIAL_SNARL else 1
                var at = p + 8 if t == DISTANCED_TRIVIAL_SNARL else p + 3
                if k < 1 or k > 255:
                    raise Error(t"{self.path}: a damaged distance index (a run of {k} chain nodes)")
                for j in range(k):
                    words.append(UInt64(self.entry(at + stride * j)))
            elif _is_snarl(t) or _is_simple(t):
                if len(words) == 0:
                    raise Error(t"{self.path}: a damaged distance index (a chain starting with a snarl)")
                snarl_at.append(p)
                snarl_at.append(len(words))
            else:
                raise Error(t"{self.path}: a damaged distance index (a chain child of type {t})")
            if p == last:
                break
            var size = self.entry(p - 1)
            if size <= 0:
                raise Error(t"{self.path}: a damaged distance index (a chain's sizes)")
            p += size + 2
            steps += 1
            if p > self.length or steps > self.length:
                raise Error(t"{self.path}: a damaged distance index (a chain runs past its end)")
        if len(words) == 0:
            raise Error(t"{self.path}: a damaged distance index (a chain with no nodes)")

    def _chain(self, chain: Int, parent_start: UInt64, parent_end: UInt64, mut out: List[Snarl]) raises:
        var words = List[UInt64]()
        var snarl_at = List[Int]()
        self.chain_words(chain, words, snarl_at)
        var looping = (self.entry(chain + 6) >> 1) == (self.entry(chain + 7) >> 1)
        # vg's vote: every consecutive pair of the chain's nodes is a snarl of its chain.
        var pairs = len(words) - 1
        var ascending = 0
        for k in range(len(words) - 1):
            if words[k] >> 1 <= words[k + 1] >> 1:
                ascending += 1
        if looping and len(words) > 1:
            pairs += 1
            if words[len(words) - 1] >> 1 <= words[0] >> 1:
                ascending += 1
        var flip = 2 * ascending < pairs
        for i in range(len(snarl_at) // 2):
            var record = snarl_at[2 * i]
            var k = snarl_at[2 * i + 1]
            var s = words[k - 1]
            var e: UInt64
            if k < len(words):
                e = words[k]
            elif looping:
                e = words[0]
            else:
                raise Error(t"{self.path}: a damaged distance index (a chain ending with a snarl)")
            if flip:
                var t = s
                s = e ^ 1
                e = t ^ 1
            out.append(Snarl(s, e, parent_start, parent_end))
            if _is_snarl(self.kind(record)):
                var n = self.entry(record + 1)
                var list = self.entry(record + 7)
                for c in range(n):
                    var child = self.entry(list + c)
                    var ct = self.kind(child)
                    if _is_chain(ct):
                        self._chain(child, s, e, out)
                    elif not _is_node(ct):
                        raise Error(t"{self.path}: a damaged distance index (a snarl child of type {ct})")
