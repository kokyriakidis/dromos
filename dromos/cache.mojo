"""The graph's cache: what dromos reads of a GBZ graph, laid out once for memory mapping (made by prepare.mojo, kept as
`G.gbz.dromos`). Commands map it and read only the sections they need: index the reference sample's placements and
the node lengths, gam-to-gaf the lengths and letters, surject all of these.

The file is little-endian. A preamble of PREAMBLE bytes, then sections, each at a multiple of ALIGN bytes, then the
table of sections:

    "DROMOSGC"                  magic
    u64 version                 VERSION; a cache of another version is made again
    u64 size                    the file's size
    u64 table_at, u64 sections  the table of sections: u32 kind, u32 sample, u64 offset, u64 length each
    u64 graph_size, u64 graph_mtime, u64 graph_hash
                                the graph it was made from (`GraphId`): a cache whose graph differs is stale
    u64 checksum                FNV-1a of the preamble (this field 0) and the table
    u64 n, n bytes              the graph's file name

The sections (kind; those of a reference sample carry its number in SAMPLES):

    NODES       u64 count (node ids 0 ..< count), u64 width (letter bits: 4 or 8), u64 letters (their index space),
                u8 table[16] (the letter of each 4-bit code)
    SAMPLES     u64 n, then n times u64 n, n bytes: the reference samples, in the graph's order
    LENGTHS     u32[count rounded up to BLOCK]: each node's length; NONE for an id the graph lacks (and the padding)
    STARTS      u64[count / BLOCK, rounded up]: where the letters of node BLOCK * i start; a node's start is its
                block's plus the lengths before it in the block
    LETTERS     every node's letters, forward strand, in id order: 4-bit codes two a byte (the first in the low half)
                or bytes; with 4 bits, each NODE_CHUNK of ids starts at an even index (a byte of its own)
    PATHS       u64 paths, u64 path_first[paths + 1] (each path's first placement), u32 path_contig[paths] (padded to
                8 bytes), u64 contigs, then each contig's name (u64 n, n bytes): the sample's paths in path id order,
                every phase, contigs numbered as first met
    FIRST       u32[count]: the node's first placement plus one; 0 for a node off the sample's paths
    PLACEMENTS  u64[placements]: start << START_SHIFT | path << 1 | reverse, in path order, each path's in step order
    MORE_NODES  u64[k]: the nodes the sample's paths visit more than once, in id order
    MORE_AT     u64[k + 1]: where each one's further placements start in MORE
    MORE        u32[]: the placements after each such node's first, in placement order
"""

from std.ffi import external_call
from std.os import stat
from std.sys.info import CompilationTarget

from dromos.bgzf import MappedFile, RawFile
from dromos.graph import GraphId
from dromos.output import OutFile, errno

comptime MAGIC = "DROMOSGC"
comptime VERSION = 1
"""The format's version: any change to the layout raises it."""
comptime PREAMBLE = 4096
comptime ALIGN = 4096
comptime NONE = 0xFFFF_FFFF
"""A LENGTHS or FIRST entry with no value."""
comptime BLOCK = 32
"""Nodes per entry of STARTS."""
comptime NODE_CHUNK = 1 << 20
"""Node ids per task when the cache is made (a multiple of BLOCK)."""
comptime START_SHIFT = 24
comptime PATH_BITS = 23
"""A placement holds a path number below 2^PATH_BITS and a start below 2^(64 - START_SHIFT)."""
comptime CODES = "ACGTNRYKMSWBDHVU"
"""The letters 4-bit codes stand for, in code order; a graph with another letter keeps bytes."""

comptime NODES = 1
comptime SAMPLES = 2
comptime LENGTHS = 3
comptime STARTS = 4
comptime LETTERS = 5
comptime PATHS = 10
comptime FIRST = 11
comptime PLACEMENTS = 12
comptime MORE_NODES = 13
comptime MORE_AT = 14
comptime MORE = 15

comptime _MAX_NAME = PREAMBLE - 80
comptime _SECTION_BYTES = 24


@fieldwise_init
struct Section(Copyable, TrivialRegisterPassable):
    var kind: Int
    var sample: Int
    var offset: Int
    var length: Int


@always_inline
def u64_at(p: ImmPointer[UInt8, _], at: Int) -> Int:
    return Int(p.unsafe_offset(at).unsafe_bitcast[UInt64]()[])


@always_inline
def u32_at(p: ImmPointer[UInt8, _], at: Int) -> Int:
    return Int(p.unsafe_offset(at).unsafe_bitcast[UInt32]()[])


def put_u64(mut out: List[UInt8], value: Int):
    for i in range(8):
        out.append(UInt8((value >> (8 * i)) & 0xFF))


def put_u32(mut out: List[UInt8], value: Int):
    for i in range(4):
        out.append(UInt8((value >> (8 * i)) & 0xFF))


def put_text(mut out: List[UInt8], text: String):
    put_u64(out, text.byte_length())
    out.extend(text.as_bytes())


def _fnv(h: UInt64, data: Span[UInt8, _]) -> UInt64:
    var x = h
    for b in data:
        x = (x ^ UInt64(b)) * 0x100000001B3
    return x


def _checksum(preamble: Span[UInt8, _], table: Span[UInt8, _]) -> UInt64:
    """The preamble's checksum: FNV-1a of its bytes, the checksum's own 8 as 0, then the table's."""
    var h = _fnv(UInt64(0xCBF29CE484222325), preamble[:64])
    h = _fnv(h, Span(List[UInt8](length=8, fill=0)))
    h = _fnv(h, preamble[72:])
    return _fnv(h, table)


# Why a file is not a usable cache, as `examine` finds it.
comptime CURRENT = 0
comptime MISSING = 1
comptime NOT_A_CACHE = 2
"""Not a dromos cache at all: never written over."""
comptime OTHER_VERSION = 3
comptime DAMAGED = 4
"""A dromos cache cut short or damaged."""


struct Header(Movable):
    """A cache's preamble and table of sections, read and checked."""

    var size: Int
    var graph: GraphId
    var name: String
    var sections: List[Section]

    def __init__(out self):
        self.size = 0
        self.graph = GraphId(0, 0, 0)
        self.name = String()
        self.sections = List[Section]()

    def find(self, kind: Int, sample: Int = 0) -> Optional[Section]:
        for s in self.sections:
            if s.kind == kind and s.sample == sample:
                return s
        return None


def read_header(file: RawFile, mut header: Header) -> Int:
    """Reads and checks a cache's preamble and table: CURRENT if they are whole and agree, else why not (the graph is
    not compared)."""
    var head = List[UInt8](length=PREAMBLE, fill=0)
    var got = file.read(head.unsafe_ptr(), PREAMBLE, 0)
    var magic = MAGIC.as_bytes()
    if got < 8 or Span(head)[:8] != magic:
        return NOT_A_CACHE
    if got < PREAMBLE:
        return DAMAGED
    var p = head.unsafe_ptr().as_imm()
    if u64_at(p, 8) != VERSION:
        return OTHER_VERSION
    var size = u64_at(p, 16)
    var table_at = u64_at(p, 24)
    var count = u64_at(p, 32)
    var name_length = u64_at(p, 72)
    if size != file.size or name_length > _MAX_NAME or count > 1 << 16:
        return DAMAGED
    if table_at < PREAMBLE or table_at + _SECTION_BYTES * count != size:
        return DAMAGED
    var table = List[UInt8](length=_SECTION_BYTES * count, fill=0)
    if file.read(table.unsafe_ptr(), len(table), table_at) != len(table):
        return DAMAGED
    if UInt64(u64_at(p, 64)) != _checksum(Span(head), Span(table)):
        return DAMAGED
    header.size = size
    header.graph = GraphId(u64_at(p, 40), u64_at(p, 48), UInt64(u64_at(p, 56)))
    header.name = String(StringSlice(unsafe_from_utf8=Span(head)[80 : 80 + name_length]))
    header.sections.clear()
    var t = table.unsafe_ptr().as_imm()
    for i in range(count):
        var at = _SECTION_BYTES * i
        var s = Section(u32_at(t, at), u32_at(t, at + 4), u64_at(t, at + 8), u64_at(t, at + 16))
        if s.offset < PREAMBLE or s.offset % ALIGN != 0 or s.length > table_at - s.offset:
            return DAMAGED
        header.sections.append(s)
    return CURRENT


def examine(path: String, mut header: Header) -> Int:
    """Whether `path` holds a cache (`read_header`), or MISSING."""
    try:
        _ = stat(path)
    except:
        return MISSING
    try:
        var file = RawFile(path)
        if not file.seekable:
            return NOT_A_CACHE
        return read_header(file, header)
    except:
        return NOT_A_CACHE


def _need(header: Header, kind: Int, path: String, sample: Int = 0) raises -> Section:
    var s = header.find(kind, sample)
    if not s:
        raise Error(t"{path}: a damaged graph cache (a section missing); make it again with dromos prepare")
    return s.value()


struct CacheFile(Movable):
    """A cache, mapped into memory (nothing read ahead until `advise`): its sections are read in place."""

    var file: MappedFile
    var path: String
    var header: Header
    var count: Int
    """Node ids 0 ..< count."""
    var width: Int
    var table: SIMD[DType.uint8, 16]
    var samples: List[String]

    def __init__(out self, path: String) raises:
        var header = Header()
        var why = read_header(RawFile(path), header)
        var damaged = String(t"{path}: a damaged graph cache; make it again with dromos prepare")
        if why == NOT_A_CACHE:
            raise Error(t"{path}: not a dromos graph cache")
        if why == OTHER_VERSION:
            raise Error(t"{path}: a graph cache of another dromos version; make it again with dromos prepare")
        if why != CURRENT:
            raise Error(damaged)
        var file = MappedFile(path, ahead=False)
        if file.size != header.size:
            raise Error(damaged)
        var nodes = _need(header, NODES, path)
        var p = file.ptr.unsafe_offset(nodes.offset)
        if nodes.length < 40:
            raise Error(damaged)
        self.count = u64_at(p, 0)
        self.width = u64_at(p, 8)
        self.table = p.unsafe_offset(24).unsafe_load[width=16]()
        var s = _need(header, SAMPLES, path)
        var q = file.ptr.unsafe_offset(s.offset)
        var n = u64_at(q, 0) if s.length >= 8 else 0
        self.samples = List[String](capacity=n)
        var at = 8
        for _ in range(n):
            var length = u64_at(q, at) if at + 8 <= s.length else s.length
            if at + 8 + length > s.length:
                raise Error(damaged)
            self.samples.append(
                String(StringSlice(unsafe_from_utf8=Span(unsafe_ptr=q.unsafe_offset(at + 8), length=length)))
            )
            at += 8 + length
        var rounded = (self.count + BLOCK - 1) // BLOCK * BLOCK
        if (
            _need(header, LENGTHS, path).length != 4 * rounded
            or _need(header, STARTS, path).length != 8 * (rounded // BLOCK)
            or (self.width != 4 and self.width != 8)
        ):
            raise Error(damaged)
        self.file = file^
        self.path = path
        self.header = header^

    def section(self, kind: Int, sample: Int = 0) raises -> Section:
        return _need(self.header, kind, self.path, sample)

    @always_inline
    def ptr(self, s: Section) -> ImmPointer[UInt8, ImmUntrackedOrigin]:
        return self.file.ptr.unsafe_offset(s.offset)

    def advise(self, s: Section):
        """Has the kernel read section `s` ahead."""
        self.file.advise(s.offset, s.length)

    def sample_index(self, name: String) -> Optional[Int]:
        for i in range(len(self.samples)):
            if self.samples[i] == name:
                return i
        return None


comptime _PROT_READ_WRITE: Int32 = 3
comptime _MAP_SHARED: Int32 = 1
comptime _MADV_DONTNEED: Int32 = 4
comptime _PAGE = 4096
comptime _EOPNOTSUPP: Int32 = 95


struct OutMap(Movable):
    """A region of a file mapped for writing: what is written there is the file's (buffered writes to one file take
    its lock in turn; a mapping's pages are written from every thread at once). Unmapped when it goes."""

    var ptr: MutPointer[UInt8, MutUntrackedOrigin]
    var length: Int
    var address: Int

    def __init__(out self, fd: Int32, offset: Int, length: Int) raises:
        self.length = length
        self.address = 0
        if length > 0:
            self.address = external_call["mmap", Int](Int(0), length, _PROT_READ_WRITE, _MAP_SHARED, fd, offset)
            if self.address == -1:
                raise Error(t"cannot map the graph cache for writing (errno {errno()})")
        self.ptr = MutPointer[UInt8, MutUntrackedOrigin](unsafe_from_address=self.address if length > 0 else 8)

    def release(self, offset: Int, length: Int):
        """Lets go of the whole pages of bytes `offset ..< offset + length`, written: they stay in the page cache, to
        be written out, but no longer in this process's memory."""
        var lo = (offset + _PAGE - 1) // _PAGE * _PAGE
        var hi = (offset + length) // _PAGE * _PAGE
        if hi > lo:
            _ = external_call["madvise", Int32](self.ptr.unsafe_offset(lo), hi - lo, _MADV_DONTNEED)

    def __deinit__(deinit self):
        if self.address != 0:
            _ = external_call["munmap", Int32](self.ptr, self.length)


struct CacheWriter(Movable):
    """A cache being written: sections laid out one after another, written at their offsets (or through a mapping)
    from any thread, then the table and the preamble (`finish`)."""

    var file: OutFile
    var at: Int
    var size: Int
    """The file's size so far, as `map` grew it."""
    var sections: List[Section]

    def __init__(out self, path: String) raises:
        self.file = OutFile(path, mappable=True)
        self.at = PREAMBLE
        self.size = 0
        self.sections = List[Section]()

    def map(mut self, offset: Int, length: Int) raises -> OutMap:
        """Bytes `offset ..< offset + length` of the file (a section's), mapped for writing; the file grows to hold
        them."""
        if offset + length > self.size:
            comptime truncate = "ftruncate64" if CompilationTarget.is_linux() else "ftruncate"
            # Room taken on the disk now: a full disk is an error here, not a fault when a mapped page is written.
            comptime if CompilationTarget.is_linux():
                var grow = offset + length - self.size
                if external_call["fallocate", Int32](self.file.fd, Int32(0), self.size, grow) != 0:
                    var e = errno()
                    if e != _EOPNOTSUPP:
                        raise Error(t"cannot write {self.file.path} (errno {e})")
            if external_call[truncate, Int32](self.file.fd, offset + length) != 0:
                raise Error(t"cannot write {self.file.path} (errno {errno()})")
            self.size = offset + length
        return OutMap(self.file.fd, offset, length)

    def reserve(mut self, kind: Int, sample: Int, length: Int) -> Int:
        """Room for a section of `length` bytes: its offset."""
        var offset = self.at
        self.sections.append(Section(kind, sample, offset, length))
        self.at = (offset + length + ALIGN - 1) // ALIGN * ALIGN
        return offset

    def drop_last(mut self):
        """Takes the last section reserved back."""
        var s = self.sections.pop()
        self.at = s.offset

    def add(mut self, kind: Int, sample: Int, data: Span[UInt8, _]) raises:
        """A section of `data`."""
        self.file.put(data, self.reserve(kind, sample, len(data)))

    def finish(mut self, graph: GraphId, name: String) raises:
        """Writes the table, then the preamble, and closes the file."""
        if name.byte_length() > _MAX_NAME:
            raise Error(t"a graph file name of more than {_MAX_NAME} bytes")
        var table = List[UInt8]()
        for s in self.sections:
            put_u32(table, s.kind)
            put_u32(table, s.sample)
            put_u64(table, s.offset)
            put_u64(table, s.length)
        var head = List[UInt8]()
        head.extend(MAGIC.as_bytes())
        put_u64(head, VERSION)
        put_u64(head, self.at + len(table))
        put_u64(head, self.at)
        put_u64(head, len(self.sections))
        put_u64(head, graph.size)
        put_u64(head, graph.mtime)
        put_u64(head, Int(graph.hash))
        put_u64(head, 0)
        put_text(head, name)
        head.resize(PREAMBLE, 0)
        var sum = _checksum(Span(head), Span(table))
        for i in range(8):
            head[64 + i] = UInt8((sum >> UInt64(8 * i)) & 0xFF)
        self.file.put(Span(table), self.at)
        self.file.put(Span(head), 0)
        self.file.close()
