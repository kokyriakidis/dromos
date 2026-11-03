"""A GBZ graph, read by the gbz crate (Rust, behind gbz-ffi's C interface: libdromos_gbz): its tags, path metadata,
path walks, node sequences, edges and bidirectional search. Everything above the crate is here, in Mojo: what
pgcatalog's graph.rs learns from a GBZ without loading it (reference samples and paths, the reference sample to use, a
fingerprint), ported line for line.

The library is found beside the dromos binary, then at $DROMOS_GBZ_LIB.
"""

from std.ffi import OwnedDLHandle, _DLCallable, external_call
from std.memory import unsafe_memcpy
from std.os import getenv, stat
from std.os.path import exists
from std.sys.info import CompilationTarget

from dromos.bgzf import RawFile

comptime ABI = 2
"""The C interface's version (gbz-ffi's `ABI`)."""
comptime LIBRARY = "libdromos_gbz" + (".dylib" if CompilationTarget.is_macos() else ".so")
comptime MAGIC = "GBZ "
"""A GBZ's first four bytes."""
comptime REFERENCE_SAMPLES_KEY = "reference_samples"
"""The GBWT tag naming the reference samples (gbz::REFERENCE_SAMPLES_KEY)."""
comptime NO_LENGTH = UInt64.MAX
"""A node length `GbzGraph.lengths` gives for an id the graph lacks."""


def _c_string(s: String) -> List[UInt8]:
    var c = List[UInt8](s.as_bytes())
    c.append(0)
    return c^


def _exe_dir() -> String:
    """The folder of the running program; empty if it cannot be told."""
    var buf = List[UInt8](length=4096, fill=0)
    var n = 0
    comptime if CompilationTarget.is_macos():
        var size = UInt32(4096)
        if external_call["_NSGetExecutablePath", Int32](buf.unsafe_ptr(), Pointer(to=size)) == 0:
            n = Int(external_call["strlen", Int](buf.unsafe_ptr()))
    else:
        var proc = _c_string("/proc/self/exe")
        n = max(external_call["readlink", Int](proc.unsafe_ptr(), buf.unsafe_ptr(), 4095), 0)
    var path = String(StringSlice(unsafe_from_utf8=Span(buf)[:n]))
    var slash = path.rfind("/")
    return String(path[byte=:slash]) if slash > 0 else String("/" if slash == 0 else "")


def is_gbz(path: String) -> Bool:
    """Whether the file at `path` is a GBZ, told by its first bytes (not its name)."""
    try:
        var f = RawFile(path)
        if not f.seekable:
            return False
        var head = List[UInt8](length=4, fill=0)
        return f.read(head.unsafe_ptr(), 4, 0) == 4 and Span(head) == MAGIC.as_bytes()
    except:
        return False


struct GbzLib(Movable):
    """libdromos_gbz, loaded: every call into the crate goes through one of these."""

    var lib: OwnedDLHandle

    def __init__(out self) raises:
        var names = List[String]()
        var dir = _exe_dir()
        if dir:
            names.append(dir + "/" + LIBRARY)
        var env = getenv("DROMOS_GBZ_LIB")
        if env:
            names.append(env)
        for name in names:
            if exists(name):
                self.lib = OwnedDLHandle(name)
                var abi = Int(self.lib.get_function[UInt32]("dgbz_abi")())
                if abi != ABI:
                    raise Error(
                        t"{name} is gbz-ffi version {abi}, this dromos needs {ABI}: rebuild with pixi run build"
                    )
                return
        raise Error(t"{LIBRARY} not found beside dromos ({dir}) or at DROMOS_GBZ_LIB: build it with pixi run build")

    def error(self) -> String:
        """The message of this thread's last failed call."""
        try:
            var length = 0
            var p = self.lib.get_function[Int]("dgbz_error")(Pointer(to=length))
            return String(_borrowed(p, length))
        except:
            return String("gbz-ffi failed")

    def check(self, status: Int32) raises:
        if status != 0:
            raise Error(self.error())

    def take_bytes(self, mut buf: _Buffer) raises -> List[UInt8]:
        """The bytes of an FFI buffer, which is freed."""
        var out = List[UInt8](unsafe_uninit_length=buf.len)
        if buf.len > 0:
            unsafe_memcpy(dest=out.unsafe_ptr(), src=_address[UInt8](buf.ptr), count=buf.len)
        self.lib.get_function[NoneType]("dgbz_bytes_free")(Pointer(to=buf))
        return out^

    def take_u64s(self, mut buf: _Buffer) raises -> List[UInt64]:
        """The u64s of an FFI buffer, which is freed."""
        var out = List[UInt64](unsafe_uninit_length=buf.len)
        if buf.len > 0:
            unsafe_memcpy(dest=out.unsafe_ptr(), src=_address[UInt64](buf.ptr), count=buf.len)
        self.lib.get_function[NoneType]("dgbz_u64s_free")(Pointer(to=buf))
        return out^


def _call[R: TrivialRegisterPassable](lib: GbzLib, name: StaticString) raises -> _DLCallable[R, ImmUntrackedOrigin]:
    """A function of the library, to call many times: valid while the library is loaded (its owner lives)."""
    return rebind[_DLCallable[R, ImmUntrackedOrigin]](lib.lib.get_function[R](name))


@fieldwise_init
struct _Buffer(TrivialRegisterPassable):
    """gbz-ffi's DgbzU64s and DgbzBytes: an array the library made, freed by it."""

    var ptr: Int
    var len: Int


@always_inline
def _address[T: AnyType](address: Int) -> ImmPointer[T, ImmUntrackedOrigin]:
    return ImmPointer[T, ImmUntrackedOrigin](unsafe_from_address=address)


def _borrowed(address: Int, length: Int) -> StringSlice[ImmUntrackedOrigin]:
    """The `length` bytes at `address` (none if it is 0), as text."""
    if address == 0 or length == 0:
        return StringSlice[ImmUntrackedOrigin]()
    return StringSlice(unsafe_from_utf8=Span(unsafe_ptr=_address[UInt8](address), length=length))


@fieldwise_init
struct PathName(Copyable, TrivialRegisterPassable):
    """A path's name in the metadata: sample, contig, phase (haplotype) and fragment (a fragment's first base)."""

    var sample: Int
    var contig: Int
    var phase: Int
    var fragment: Int


struct Metadata(Movable):
    """A GBZ's path metadata: read from the file on its own (`read_metadata`), or the loaded graph's (valid while the
    graph lives)."""

    var lib: GbzLib
    var ptr: Int
    var owned: Bool

    def __init__(out self, var lib: GbzLib, ptr: Int, owned: Bool):
        self.lib = lib^
        self.ptr = ptr
        self.owned = owned

    def __deinit__(deinit self):
        if self.owned:
            try:
                self.lib.lib.get_function[NoneType]("dgbz_metadata_free")(self.ptr)
            except:
                pass

    def paths(self) raises -> Int:
        return self.lib.lib.get_function[Int]("dgbz_md_paths")(self.ptr)

    def samples(self) raises -> Int:
        return self.lib.lib.get_function[Int]("dgbz_md_samples")(self.ptr)

    def contigs(self) raises -> Int:
        return self.lib.lib.get_function[Int]("dgbz_md_contigs")(self.ptr)

    def path_names(self) raises -> List[PathName]:
        """Every path's name, by path id."""
        var n = self.paths()
        var raw = List[UInt64](length=4 * n, fill=0)
        var got = self.lib.lib.get_function[Int]("dgbz_md_path_names")(self.ptr, 0, n, raw.unsafe_ptr())
        var out = List[PathName](capacity=got)
        for i in range(got):
            out.append(PathName(Int(raw[4 * i]), Int(raw[4 * i + 1]), Int(raw[4 * i + 2]), Int(raw[4 * i + 3])))
        return out^

    def sample_name(self, id: Int) raises -> String:
        """Sample `id`'s name, or its number when it has none (the crate's `sample_name`)."""
        var length = 0
        var p = self.lib.lib.get_function[Int]("dgbz_md_sample")(self.ptr, id, Pointer(to=length))
        return String(_borrowed(p, length)) if p != 0 else String(id)

    def contig_name(self, id: Int) raises -> String:
        """Contig `id`'s name, or its number when it has none (the crate's `contig_name`)."""
        var length = 0
        var p = self.lib.lib.get_function[Int]("dgbz_md_contig")(self.ptr, id, Pointer(to=length))
        return String(_borrowed(p, length)) if p != 0 else String(id)

    def sample_id(self, name: String) raises -> Optional[Int]:
        var id = self.lib.lib.get_function[Int]("dgbz_md_sample_id")(self.ptr, name.unsafe_ptr(), name.byte_length())
        return Optional(id) if id >= 0 else None

    def contig_id(self, name: String) raises -> Optional[Int]:
        var id = self.lib.lib.get_function[Int]("dgbz_md_contig_id")(self.ptr, name.unsafe_ptr(), name.byte_length())
        return Optional(id) if id >= 0 else None


@fieldwise_init
struct SearchState(Copyable, Equatable, TrivialRegisterPassable):
    """A state of bidirectional GBWT search (gbz-ffi's DgbzState): the last node each way, GBWT-encoded
    (`node << 1 | reverse`), and its range. `count` is how many haplotypes it holds."""

    var forward_node: UInt64
    var forward_start: UInt64
    var forward_end: UInt64
    var reverse_node: UInt64
    var reverse_start: UInt64
    var reverse_end: UInt64

    def __init__(out self):
        self = Self(0, 0, 0, 0, 0, 0)

    @always_inline
    def count(self) -> Int:
        return Int(self.forward_end - self.forward_start)

    @always_inline
    def node(self) -> Int:
        """The forward node's id."""
        return Int(self.forward_node >> 1)

    @always_inline
    def reverse(self) -> Bool:
        """Whether the forward node is visited in reverse."""
        return self.forward_node & 1 == 1


struct GbzGraph(Movable):
    """A GBZ graph, loaded whole by the crate. Read-only: any thread may call it."""

    var lib: GbzLib
    var ptr: Int

    def __init__(out self, path: String) raises:
        """Loads the graph at `path`; a failure's message is the crate's, after the path (as pgcatalog's)."""
        self.lib = GbzLib()
        self.ptr = 0
        var status = self.lib.lib.get_function[Int32]("dgbz_load")(
            path.unsafe_ptr(), path.byte_length(), Pointer(to=self.ptr)
        )
        if status != 0:
            raise Error(t"{path}: {self.lib.error()}")

    def __deinit__(deinit self):
        try:
            self.lib.lib.get_function[NoneType]("dgbz_free")(self.ptr)
        except:
            pass

    def metadata(self) raises -> Metadata:
        """The graph's path metadata, valid while the graph lives."""
        var md = self.lib.lib.get_function[Int]("dgbz_metadata")(self.ptr)
        if md == 0:
            raise Error("the GBZ has no path metadata")
        return Metadata(GbzLib(), md, False)

    def nodes(self) raises -> Int:
        return self.lib.lib.get_function[Int]("dgbz_nodes")(self.ptr)

    def min_node(self) raises -> Int:
        return self.lib.lib.get_function[Int]("dgbz_min_node")(self.ptr)

    def max_node(self) raises -> Int:
        return self.lib.lib.get_function[Int]("dgbz_max_node")(self.ptr)

    def has_node(self, node: Int) raises -> Bool:
        return self.lib.lib.get_function[Int32]("dgbz_has_node")(self.ptr, node) != 0

    def sequence(self, node: Int) raises -> Span[UInt8, ImmUntrackedOrigin]:
        """Node `node`'s letters on its forward strand, in the graph's memory; empty for a node the graph lacks."""
        var length = 0
        var p = self.lib.lib.get_function[Int]("dgbz_sequence")(self.ptr, node, Pointer(to=length))
        if p == 0:
            return Span[UInt8, ImmUntrackedOrigin]()
        return Span(unsafe_ptr=_address[UInt8](p), length=length)

    def lengths(self, first: Int, count: Int, dst: MutPointer[UInt64, _]) raises:
        """The lengths of nodes `first ..< first + count` into `dst`; NO_LENGTH for an id the graph lacks."""
        self.lib.lib.get_function[NoneType]("dgbz_sequence_lengths")(self.ptr, first, count, dst)

    def sequences(self, first: Int, count: Int, address: MutPointer[UInt64, _], length: MutPointer[UInt64, _]) raises:
        """The letters of nodes `first ..< first + count`, in the graph's memory: each one's address into `address`
        (0 for an id the graph lacks) and its length into `length`."""
        self.lib.lib.get_function[NoneType]("dgbz_sequences")(self.ptr, first, count, address, length)

    def paths(self) raises -> Int:
        return self.lib.lib.get_function[Int]("dgbz_paths")(self.ptr)

    def path(self, id: Int, reverse: Bool = False) raises -> List[UInt64]:
        """Path `id`'s visits (each `node << 1 | reverse`), walked in the GBWT; its reverse when `reverse`."""
        var buf = _Buffer(0, 0)
        self.lib.check(self.lib.lib.get_function[Int32]("dgbz_path")(self.ptr, id, Int32(reverse), Pointer(to=buf)))
        return self.lib.take_u64s(buf)

    def _edges(self, name: StaticString, node: Int, reverse: Bool, mut out: List[UInt64]) raises -> Bool:
        out.resize(max(len(out), 16), 0)
        var n = self.lib.lib.get_function[Int]("dgbz_" + name)(
            self.ptr, node, Int32(reverse), out.unsafe_ptr(), len(out)
        )
        if n > len(out):
            out.resize(n, 0)
            n = self.lib.lib.get_function[Int]("dgbz_" + name)(self.ptr, node, Int32(reverse), out.unsafe_ptr(), n)
        if n < 0:
            out.clear()
            return False
        out.resize(n, 0)
        return True

    def successors(self, node: Int, reverse: Bool, mut out: List[UInt64]) raises -> Bool:
        """The visits that may follow (node, reverse), each `node << 1 | reverse`, into `out`; false for a node the
        graph lacks."""
        return self._edges("successors", node, reverse, out)

    def predecessors(self, node: Int, reverse: Bool, mut out: List[UInt64]) raises -> Bool:
        """The visits that may precede (node, reverse), as `successors`."""
        return self._edges("predecessors", node, reverse, out)

    def search_state(self, node: Int, reverse: Bool) raises -> Optional[SearchState]:
        """The state of the haplotypes visiting (node, reverse); none if no haplotype does."""
        var s = SearchState()
        if self.lib.lib.get_function[Int32]("dgbz_search_state")(self.ptr, node, Int32(reverse), Pointer(to=s)) == 0:
            return None
        return s

    def _follow(self, name: StaticString, state: SearchState, mut out: List[SearchState]) raises -> Bool:
        out.resize(max(len(out), 16), SearchState())
        var n = self.lib.lib.get_function[Int]("dgbz_" + name)(self.ptr, Pointer(to=state), out.unsafe_ptr(), len(out))
        if n > len(out):
            out.resize(n, SearchState())
            n = self.lib.lib.get_function[Int]("dgbz_" + name)(self.ptr, Pointer(to=state), out.unsafe_ptr(), n)
        if n < 0:
            out.clear()
            return False
        out.resize(n, SearchState())
        return True

    def follow_forward(self, state: SearchState, mut out: List[SearchState]) raises -> Bool:
        """The non-empty states `state` extends to, one per successor its haplotypes take, into `out`; false when its
        node is not in the graph."""
        return self._follow("follow_forward", state, out)

    def follow_backward(self, state: SearchState, mut out: List[SearchState]) raises -> Bool:
        """As `follow_forward`, extending the state backward."""
        return self._follow("follow_backward", state, out)


trait NodeSource:
    """Where the graph's cache takes node lengths and letters from: the GBZ's parts, or nodes held in memory."""

    def lengths_into(
        self, first: Int, count: Int, dst: MutPointer[UInt32, MutUntrackedOrigin]
    ) raises -> Tuple[Int, Int]:
        """The lengths of nodes `first ..< first + count` into `dst`, 0xFFFF_FFFF for an id the graph lacks: their
        total, and the first id whose length a u32 cannot hold (-1 for none)."""
        ...

    def letters_into(
        self, first: Int, count: Int, dst: MutPointer[UInt8, MutUntrackedOrigin], capacity: Int
    ) raises -> Int:
        """The letters of nodes `first ..< first + count` the graph has, one after another, into `dst` (room for
        `capacity`): how many."""
        ...


struct GbzParts(Movable, NodeSource):
    """A GBZ's parts for reading every node and the reference paths fast (gbz-ffi's Parts): the GBWT, loaded by the
    crate, and the node sequences' pieces, read at the same time and unpacked as they are asked for. Node ids, real
    nodes, lengths and letters mean what they mean in `GbzGraph`. Read-only: any thread may call it."""

    var lib: GbzLib
    var ptr: Int
    var _follow: _DLCallable[Int, ImmUntrackedOrigin]
    var _edges: _DLCallable[Int, ImmUntrackedOrigin]
    var _search: _DLCallable[Int32, ImmUntrackedOrigin]
    """The search's calls, looked up once (a lookup costs more than a small call's work, and takes a lock)."""

    def __init__(out self, path: String) raises:
        self.lib = GbzLib()
        self.ptr = 0
        self._follow = _call[Int](self.lib, "dgbz_parts_follow_forward")
        self._edges = _call[Int](self.lib, "dgbz_parts_edges")
        self._search = _call[Int32](self.lib, "dgbz_parts_search_state")
        var status = self.lib.lib.get_function[Int32]("dgbz_parts_load")(
            path.unsafe_ptr(), path.byte_length(), Pointer(to=self.ptr)
        )
        if status != 0:
            raise Error(self.lib.error())

    def __deinit__(deinit self):
        try:
            self.lib.lib.get_function[NoneType]("dgbz_parts_free")(self.ptr)
        except:
            pass

    def metadata(self) raises -> Metadata:
        """The graph's path metadata, valid while the parts live."""
        var md = self.lib.lib.get_function[Int]("dgbz_parts_metadata")(self.ptr)
        if md == 0:
            raise Error("the GBZ has no path metadata")
        return Metadata(GbzLib(), md, False)

    def max_node(self) raises -> Int:
        return self.lib.lib.get_function[Int]("dgbz_parts_max_node")(self.ptr)

    def find_nodes(self, threads: Int) raises:
        """Finds which node ids the graph has (a scan of the GBWT's records) on `threads` threads, if not found yet:
        `lengths` and `letters` do it first otherwise, on every core."""
        self.lib.lib.get_function[NoneType]("dgbz_parts_find_nodes")(self.ptr, threads)

    def paths(self) raises -> Int:
        return self.lib.lib.get_function[Int]("dgbz_parts_paths")(self.ptr)

    def path(self, id: Int, reverse: Bool = False) raises -> List[UInt64]:
        """Path `id`'s visits (each `node << 1 | reverse`), walked in the GBWT; its reverse when `reverse`."""
        var buf = _Buffer(0, 0)
        self.lib.check(
            self.lib.lib.get_function[Int32]("dgbz_parts_path")(self.ptr, id, Int32(reverse), Pointer(to=buf))
        )
        return self.lib.take_u64s(buf)

    def lengths(self, first: Int, count: Int, dst: MutPointer[UInt64, _]) raises:
        """The lengths of nodes `first ..< first + count` into `dst`; NO_LENGTH for an id the graph lacks."""
        self.lib.lib.get_function[NoneType]("dgbz_parts_lengths")(self.ptr, first, count, dst)

    def lengths_u32(self, first: Int, count: Int, absent: UInt32, dst: MutPointer[UInt32, _]) raises -> Tuple[Int, Int]:
        """The lengths of nodes `first ..< first + count` as u32 into `dst` (`absent` for an id the graph lacks, a
        length cut to 32 bits): their total, and the first id whose length a u32 cannot hold (-1 for none)."""
        var sum = UInt64(0)
        var long = self.lib.lib.get_function[Int]("dgbz_parts_lengths_u32")(
            self.ptr, first, count, absent, dst, Pointer(to=sum)
        )
        return (Int(sum), long)

    def letters(self, first: Int, count: Int, dst: MutPointer[UInt8, _], capacity: Int) raises -> Int:
        """The letters of nodes `first ..< first + count` the graph has, one after another, into `dst` (room for
        `capacity`): how many."""
        var n = self.lib.lib.get_function[Int]("dgbz_parts_letters")(self.ptr, first, count, dst, capacity)
        if n == -2:
            raise Error("a node letter outside the graph's alphabet: a damaged GBZ")
        if n < 0:
            raise Error("node letters overflow their buffer")
        return n

    def lengths_into(
        self, first: Int, count: Int, dst: MutPointer[UInt32, MutUntrackedOrigin]
    ) raises -> Tuple[Int, Int]:
        return self.lengths_u32(first, count, UInt32(0xFFFF_FFFF), dst)

    def letters_into(
        self, first: Int, count: Int, dst: MutPointer[UInt8, MutUntrackedOrigin], capacity: Int
    ) raises -> Int:
        return self.letters(first, count, dst, capacity)

    def search_state(self, node: Int, reverse: Bool) raises -> Optional[SearchState]:
        """The state of the haplotypes visiting (node, reverse) (`GbzGraph.search_state`); none if no haplotype
        does."""
        var s = SearchState()
        var found = self._search(self.ptr, node, Int32(reverse), Pointer(to=s))
        if found == 0:
            return None
        return s

    def follow_forward(
        self, states: ImmPointer[SearchState, _], n: Int, mut out: List[SearchState], mut counts: List[Int64]
    ) raises:
        """For each of the `n` states at `states`, the non-empty states it extends to forward (`GbzGraph.
        follow_forward`), one state's after another into `out`, and how many into `counts` (-1 when its node is not in
        the graph): one call into the crate for many states."""
        counts.resize(unsafe_uninit_length=n)
        out.resize(unsafe_uninit_length=max(len(out), 4 * n + 16))
        var call = self._follow
        var done = 0
        var written = 0
        while done < n:
            var got = call(
                self.ptr,
                states.unsafe_offset(done),
                n - done,
                out.unsafe_ptr().unsafe_offset(written),
                len(out) - written,
                counts.unsafe_ptr().unsafe_offset(done),
            )
            for i in range(done, done + got):
                written += max(Int(counts[i]), 0)
            done += got
            if done < n:
                out.resize(unsafe_uninit_length=2 * len(out))
        out.resize(unsafe_uninit_length=written)

    def edges(
        self, visits: ImmPointer[UInt64, _], n: Int, predecessors: Bool, mut out: List[UInt64], mut counts: List[Int64]
    ) raises:
        """For each of the `n` visits at `visits`, its successors (`GbzGraph.successors`) or, with `predecessors`, its
        predecessors, one visit's after another into `out`, and how many into `counts` (-1 for a node the graph
        lacks): one call into the crate for many visits."""
        counts.resize(unsafe_uninit_length=n)
        out.resize(unsafe_uninit_length=max(len(out), 4 * n + 16))
        var call = self._edges
        var done = 0
        var written = 0
        while done < n:
            var got = call(
                self.ptr,
                visits.unsafe_offset(done),
                n - done,
                Int32(predecessors),
                out.unsafe_ptr().unsafe_offset(written),
                len(out) - written,
                counts.unsafe_ptr().unsafe_offset(done),
            )
            for i in range(done, done + got):
                written += max(Int(counts[i]), 0)
            done += got
            if done < n:
                out.resize(unsafe_uninit_length=2 * len(out))
        out.resize(unsafe_uninit_length=written)


def gbwt_tag(path: String, key: String) raises -> Optional[String]:
    """The GBWT's tag `key`, read from the start of the file (pgcatalog's graph.rs gbwt_tags)."""
    var lib = GbzLib()
    var buf = _Buffer(0, 0)
    var found = Int32(0)
    lib.check(
        lib.lib.get_function[Int32]("dgbz_file_gbwt_tag")(
            path.unsafe_ptr(),
            path.byte_length(),
            key.unsafe_ptr(),
            key.byte_length(),
            Pointer(to=buf),
            Pointer(to=found),
        )
    )
    if found == 0:
        return None
    var bytes = lib.take_bytes(buf)
    return String(StringSlice(unsafe_from_utf8=Span(bytes)))


def read_metadata(path: String) raises -> Metadata:
    """The GBZ's path metadata, the GBWT body before it read past rather than loaded: seconds, against a whole-graph
    load."""
    var lib = GbzLib()
    var md = 0
    lib.check(lib.lib.get_function[Int32]("dgbz_file_metadata")(path.unsafe_ptr(), path.byte_length(), Pointer(to=md)))
    if md == 0:
        raise Error(t"{path}: the GBZ has no path metadata")
    return Metadata(lib^, md, True)


def reference_samples(path: String) raises -> List[String]:
    """The reference samples a GBZ declares (the GBWT's `reference_samples` tag), read from the start of the file."""
    var out = List[String]()
    var tag = gbwt_tag(path, REFERENCE_SAMPLES_KEY)
    if tag:
        for s in tag.value().split(" "):
            if s.byte_length() > 0:
                out.append(String(s))
    return out^


def reference_paths(path: String, sample: String) raises -> List[Tuple[String, List[String]]]:
    """The sample's phase-0 paths by contig, contigs in name order, each contig's in fragment order, named as vg
    names them (`CHM13#0#chr6[31350872]` for a fragment starting at 31350872), from the GBWT's metadata."""
    var md = read_metadata(path)
    var out = List[Tuple[String, List[String]]]()
    var sample_id = md.sample_id(sample)
    if not sample_id:
        return out^
    var contigs = List[String]()
    var contig_index = Dict[String, Int]()
    var fragments = List[List[Tuple[Int, String]]]()
    for p in md.path_names():
        if p.sample != sample_id.value() or p.phase != 0:
            continue
        var contig = md.contig_name(p.contig)
        var name = String(t"{sample}#0#{contig}")
        if p.fragment != 0:
            name += String(t"[{p.fragment}]")
        if contig not in contig_index:
            contig_index[contig] = len(contigs)
            contigs.append(contig)
            fragments.append(List[Tuple[Int, String]]())
        fragments[contig_index[contig]].append((p.fragment, name^))
    var order = List[Int](capacity=len(contigs))
    for i in range(len(contigs)):
        order.append(i)

    def by_name(a: Int, b: Int) {contigs} -> Bool:
        return contigs[a] < contigs[b]

    def by_fragment(a: Tuple[Int, String], b: Tuple[Int, String]) -> Bool:
        return a[0] < b[0] or (a[0] == b[0] and a[1] < b[1])

    sort(Span(order), by_name)
    for c in order:
        sort(Span(fragments[c]), by_fragment)
        var names = List[String](capacity=len(fragments[c]))
        for f in fragments[c]:
            names.append(f[1])
        out.append((contigs[c], names^))
    return out^


def resolve_reference_sample(path: String, given: Optional[String]) raises -> String:
    """The reference sample: the one given, or the graph's only declared one."""
    if given:
        return given.value()
    var found = reference_samples(path)
    if len(found) == 1:
        return found[0]
    if len(found) == 0:
        raise Error(t"{path} declares no reference sample; give one with --reference-sample")
    raise Error(t"{path} declares several reference samples ({', '.join(found)}); choose one with --reference-sample")


comptime _FINGERPRINT_SPAN = 4 << 20


@fieldwise_init
struct GraphId(Equatable, TrivialRegisterPassable):
    """What tells a graph file from another, cheaply: its size, its modification time (seconds) and an FNV-1a hash of
    its first and last 4 MB (header, tags and index structure)."""

    var size: Int
    var mtime: Int
    var hash: UInt64


def graph_id(path: String) raises -> GraphId:
    """The graph file's `GraphId` (of the file a link names)."""
    var f = RawFile(path)
    var size = f.size
    var span = min(_FINGERPRINT_SPAN, size)
    var buf = List[UInt8](length=span, fill=0)
    var h = UInt64(0xCBF29CE484222325)
    for start in [0, max(size - _FINGERPRINT_SPAN, 0)]:
        if f.read(buf.unsafe_ptr(), span, start) != span:
            raise Error(t"{path}: cannot read")
        for b in buf:
            h = (h ^ UInt64(b)) * 0x100000001B3
    return GraphId(size, Int(stat(path).st_mtimespec.tv_sec), h)


def fingerprint(path: String) raises -> String:
    """An identity for the graph's content, cheap to compute: the file size and the hash of `graph_id`."""
    var id = graph_id(path)
    comptime digits = "0123456789abcdef"
    var hex = String()
    for i in range(16):
        hex += String(digits[byte=Int((id.hash >> UInt64(60 - 4 * i)) & 15)])
    return String(t"size={id.size} fnv1a={hex}")
