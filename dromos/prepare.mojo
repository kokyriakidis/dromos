"""The graph's cache made from a GBZ (`dromos prepare`, or a command's first use of the graph; the format in
cache.mojo), and found or made again by every command that reads it (`graph_cache`).

One load of the graph's parts serves every section: the reference samples' paths are walked on threads of their own
while the nodes are found and their lengths read on every core; then the letters are unpacked and packed by range on
threads of their own while every sample's placements are laid out, a path a task, first placements set by compare
and swap and further ones gathered. Each section is written through a mapping of the file by the tasks that make it
(buffered writes to one file take its lock in turn; a mapping's pages are written from every thread at once).
"""

from std.atomic import Atomic
from std.ffi import external_call
from std.os import getenv
from std.sys import stderr
from std.time import perf_counter_ns

from dromos.cache import (
    BLOCK,
    CODES,
    CURRENT,
    DAMAGED,
    FIRST,
    LENGTHS,
    LETTERS,
    MORE,
    MORE_AT,
    MORE_NODES,
    NODE_CHUNK,
    NODES,
    NONE,
    NOT_A_CACHE,
    OTHER_VERSION,
    PATH_BITS,
    PATHS,
    PLACEMENTS,
    SAMPLES,
    START_SHIFT,
    STARTS,
    CacheFile,
    CacheWriter,
    Header,
    OutMap,
    examine,
    put_text,
    put_u32,
    put_u64,
)
from dromos.graph import GbzParts, NodeSource, graph_id, is_gbz, reference_samples
from dromos.parallel import Background, ParallelTask, parallel_for
from dromos.stream import imm_untracked, mut_untracked, untracked


struct _Timer:
    """Phase times on stderr with DROMOS_TIMING=1."""

    var on: Bool
    var t: Int

    def __init__(out self):
        self.on = Bool(getenv("DROMOS_TIMING"))
        self.t = perf_counter_ns()

    def lap(mut self, what: StaticString):
        var now = perf_counter_ns()
        if self.on:
            print(t"dromos timing: {what} {Float64(now - self.t) * 1e-9} s", file=stderr)
        self.t = now


@always_inline
def _bytes[T: AnyType, o: Origin, //](p: Pointer[T, o]) -> Pointer[UInt8, o]:
    """`p`'s memory as bytes, keeping its owner's origin."""
    return p.unsafe_bitcast[UInt8]()


# The nodes: NODES, LENGTHS, STARTS and LETTERS.


@fieldwise_init
struct _LengthsTask[S: NodeSource, origin: Origin](ParallelTask):
    """Range `index` of NODE_CHUNK ids: its lengths, their total and its first id too long for a u32."""

    var graph: Pointer[Self.S, Self.origin]
    var count: Int
    var out: MutPointer[UInt32, MutUntrackedOrigin]
    var sums: MutPointer[Int, MutUntrackedOrigin]
    var too_long: MutPointer[Int, MutUntrackedOrigin]
    var bad: MutPointer[Int, MutUntrackedOrigin]

    def run(self, index: Int):
        var first = index * NODE_CHUNK
        var n = min(NODE_CHUNK, self.count - first)
        try:
            var got = self.graph[].lengths_into(first, n, self.out.unsafe_offset(first))
            self.sums[unsafe_offset=index] = got[0]
            self.too_long[unsafe_offset=index] = got[1]
        except:
            self.bad[] = 1


comptime _CODE_LETTERS = SIMD[DType.uint8, 16](65, 67, 71, 84, 78, 82, 89, 75, 77, 83, 87, 66, 68, 72, 86, 85)
"""CODES as bytes."""
comptime _PACK = 64
"""Letters packed a vector at a time."""


@always_inline
def _codes[n: Int](v: SIMD[DType.uint8, n]) -> SIMD[DType.uint8, n]:
    """Each letter's 4-bit code (its place in CODES); 255 for a letter without one."""
    var code = SIMD[DType.uint8, n](255)
    comptime for c in range(16):
        code = v.eq(_CODE_LETTERS[c]).select(SIMD[DType.uint8, n](c), code)
    return code


def pack_letters(src: ImmPointer[UInt8, _], size: Int, dst: MutPointer[UInt8, _]) -> Bool:
    """`size` letters as 4-bit codes, two a byte (the first in the low half), into `dst`; false if a letter has no
    code."""
    var over = SIMD[DType.uint8, _PACK // 2](0)
    var i = 0
    while i + _PACK <= size:
        var code = _codes[_PACK](src.unsafe_offset(i).unsafe_load[width=_PACK]())
        var halves = code.deinterleave()
        var both = halves[0] | (halves[1] << 4)
        over |= halves[0] | halves[1]
        dst.unsafe_offset(i // 2).unsafe_store(both)
        i += _PACK
    var rest = over.reduce_or()
    while i < size:
        var lo = _codes[1](SIMD[DType.uint8, 1](src[unsafe_offset=i]))[0]
        var hi = _codes[1](SIMD[DType.uint8, 1](src[unsafe_offset=i + 1]))[0] if i + 1 < size else UInt8(0)
        rest |= lo | hi
        dst[unsafe_offset=i // 2] = lo | (hi << 4)
        i += 2
    return rest < 16


@always_inline
def _block_sum(lengths: ImmPointer[UInt32, _], block: Int) -> Int:
    """The letters of the BLOCK nodes from `block` (a multiple of BLOCK), NONE taken as 0."""
    var v = lengths.unsafe_offset(block).unsafe_load[width=BLOCK]()
    return Int(v.ne(UInt32(NONE)).select(v, 0).cast[DType.uint64]().reduce_add())


@fieldwise_init
struct _LettersTask[S: NodeSource, origin: Origin](ParallelTask):
    """Worker `index` takes ranges of NODE_CHUNK ids in turn, keeping its buffers: a range's letters are unpacked,
    packed to 4-bit codes (or kept as bytes) into the range's place in LETTERS, with the starts of its blocks. A
    letter without a code sets `bad` to 2."""

    var graph: Pointer[Self.S, Self.origin]
    var count: Int
    var width: Int
    var lengths: ImmPointer[UInt32, ImmUntrackedOrigin]
    var sums: ImmPointer[Int, ImmUntrackedOrigin]
    var chunk_start: ImmPointer[Int, ImmUntrackedOrigin]
    """Each range's first letter index (even with 4 bits)."""
    var letters: ImmPointer[OutMap, ImmUntrackedOrigin]
    var starts: MutPointer[UInt64, MutUntrackedOrigin]
    var next: MutPointer[Atomic[Int64], MutUntrackedOrigin]
    var bad: MutPointer[Int, MutUntrackedOrigin]

    def run(self, worker: Int):
        var chunks = (self.count + NODE_CHUNK - 1) // NODE_CHUNK
        var buffer = List[UInt8]()
        while self.bad[] == 0:
            var index = Int(self.next[].fetch_add(1))
            if index >= chunks:
                return
            var first = index * NODE_CHUNK
            var n = min(NODE_CHUNK, self.count - first)
            var size = self.sums[unsafe_offset=index]
            var at = self.chunk_start[unsafe_offset=index]
            # Each block's start; the last range's blocks run past `count`, over the padding's NONE.
            var s = at
            for b in range(first // BLOCK, (first + n + BLOCK - 1) // BLOCK):
                self.starts[unsafe_offset=b] = UInt64(s)
                s += _block_sum(self.lengths, b * BLOCK)
            var into: MutPointer[UInt8, MutUntrackedOrigin]
            var offset = at if self.width == 8 else at // 2
            if self.width == 8:
                into = self.letters[].ptr.unsafe_offset(at)
            else:
                buffer.resize(unsafe_uninit_length=size)
                into = mut_untracked(buffer.unsafe_ptr())
            try:
                if self.graph[].letters_into(first, n, into, size) != size:
                    self.bad[] = 1
                    return
            except:
                self.bad[] = 1
                return
            if self.width == 4 and not pack_letters(into, size, self.letters[].ptr.unsafe_offset(offset)):
                self.bad[] = 2
                return
            self.letters[].release(offset, (size + 1) // 2 if self.width == 4 else size)


comptime NO_CODE = "a node letter without a 4-bit code"
"""The error `NodeWriter.finish` raises for such a letter: the cache is then made again with bytes."""


struct NodeWriter(Movable):
    """The node sections being written: the lengths read into their section (`read_lengths`), then the letters and
    block starts by `letters` tasks (on threads of their own, beside other work, or not), then NODES (`finish`)."""

    var count: Int
    var chunks: Int
    var width: Int
    var lengths_map: OutMap
    var starts_map: OutMap
    var letters_map: OutMap
    var sums: List[Int]
    var chunk_start: List[Int]
    var total: Int
    var state: List[Int64]
    """The letters tasks' next range and failure (0, 1, or 2 for a letter without a code), on the heap: the tasks
    reach them however this moves."""

    def __init__(out self, count: Int, width: Int) raises:
        self.count = count
        self.chunks = (count + NODE_CHUNK - 1) // NODE_CHUNK
        self.width = width
        self.lengths_map = OutMap(-1, 0, 0)
        self.starts_map = OutMap(-1, 0, 0)
        self.letters_map = OutMap(-1, 0, 0)
        self.sums = List[Int](length=self.chunks, fill=0)
        self.chunk_start = List[Int](capacity=self.chunks)
        self.total = 0
        self.state = List[Int64](length=2, fill=0)

    def lengths(self) -> ImmPointer[UInt32, ImmUntrackedOrigin]:
        """Every node's length, NONE for an id the graph lacks, padded to a multiple of BLOCK: in LENGTHS, valid while
        this lives."""
        return self.lengths_map.ptr.as_imm().unsafe_origin_cast[ImmUntrackedOrigin]().unsafe_bitcast[UInt32]()

    def read_lengths[S: NodeSource](mut self, mut w: CacheWriter, graph: S, threads: Int) raises:
        """Reads the lengths into LENGTHS, and lays out STARTS and LETTERS."""
        if self.count >= Int(NONE):
            raise Error(t"{self.count} node ids: more than a graph cache holds")
        var rounded = (self.count + BLOCK - 1) // BLOCK * BLOCK
        self.lengths_map = w.map(w.reserve(LENGTHS, 0, 4 * rounded), 4 * rounded)
        var lengths = self.lengths_map.ptr.unsafe_bitcast[UInt32]()
        for i in range(self.count, rounded):
            lengths[unsafe_offset=i] = UInt32(NONE)
        var too_long = List[Int](length=self.chunks, fill=-1)
        var bad = 0
        var task = _LengthsTask(
            Pointer(to=graph),
            self.count,
            lengths,
            mut_untracked(self.sums.unsafe_ptr()),
            mut_untracked(too_long.unsafe_ptr()),
            mut_untracked(Pointer(to=bad)),
        )
        parallel_for(task, self.chunks, threads)
        if bad != 0:
            raise Error("cannot read the graph's node lengths")
        for long in too_long:
            if long >= 0:
                raise Error(t"node {long}: more letters than a graph cache holds")
        self.starts_map = w.map(w.reserve(STARTS, 0, 8 * (rounded // BLOCK)), 8 * (rounded // BLOCK))
        for c in range(self.chunks):
            self.chunk_start.append(self.total)
            self.total += self.sums[c]
            if self.width == 4:
                self.total += self.total % 2
        var size = (self.total + 1) // 2 if self.width == 4 else self.total
        self.letters_map = w.map(w.reserve(LETTERS, 0, size), size)

    def letters[S: NodeSource, o: Origin](mut self, graph: Pointer[S, o]) -> _LettersTask[S, o]:
        """The letters task, its workers to run as `workers` says; this must outlive them."""
        var state = self.state.unsafe_ptr()
        return _LettersTask(
            graph,
            self.count,
            self.width,
            self.lengths(),
            untracked(self.sums.unsafe_ptr()),
            untracked(self.chunk_start.unsafe_ptr()),
            untracked(Pointer(to=self.letters_map)),
            self.starts_map.ptr.unsafe_bitcast[UInt64](),
            state.unsafe_origin_cast[MutUntrackedOrigin]().unsafe_bitcast[Atomic[Int64]](),
            state.unsafe_offset(1).unsafe_origin_cast[MutUntrackedOrigin]().unsafe_bitcast[Int](),
        )

    def workers(self, threads: Int) -> Int:
        return max(min(threads, self.chunks), 1)

    def finish(self, mut w: CacheWriter) raises:
        """NODES, once the letters tasks are done."""
        if self.state[1] == 2:
            raise Error(NO_CODE)
        if self.state[1] != 0:
            raise Error(t"cannot write {w.file.path}: the graph's node letters")
        var info = List[UInt8]()
        put_u64(info, self.count)
        put_u64(info, self.width)
        put_u64(info, self.total)
        info.extend(CODES.as_bytes())
        w.add(NODES, 0, Span(info))


def write_nodes[S: NodeSource](mut w: CacheWriter, graph: S, count: Int, width: Int, threads: Int) raises -> NodeWriter:
    """The node sections, from `graph`'s nodes 0 ..< count, letters `width` bits each (4 or 8), in turn: done, with
    their lengths."""
    var nodes = NodeWriter(count, width)
    nodes.read_lengths(w, graph, threads)
    var task = nodes.letters(Pointer(to=graph))
    parallel_for(task, nodes.workers(threads), threads)
    nodes.finish(w)
    return nodes^


# The reference samples: PATHS, FIRST, PLACEMENTS and MORE_* each.


struct ReferenceSample(Movable):
    """A reference sample's paths in path id order (every phase): their ids in the graph, walks (each step
    `node << 1 | reverse`), fragment offsets and contigs, and the contigs' names (numbered as first met)."""

    var name: String
    var ids: List[Int]
    var walks: List[List[UInt64]]
    var fragments: List[Int]
    var path_contig: List[Int]
    var contigs: List[String]

    def __init__(out self, name: String):
        self.name = name
        self.ids = List[Int]()
        self.walks = List[List[UInt64]]()
        self.fragments = List[Int]()
        self.path_contig = List[Int]()
        self.contigs = List[String]()


@fieldwise_init
struct _PlaceTask(ParallelTask):
    """A reference path's placements (job `index`: sample << 32 | path, longest paths first, so the longest does not
    start last), written into PLACEMENTS: its steps' positions run from its fragment's offset by the nodes' lengths.
    Each step's node takes the step as its first placement if it has none or a later one; the placement that is not
    first is kept (`node << 32 | placement`) as a further one. The walk is freed."""

    var jobs: ImmPointer[Int, ImmUntrackedOrigin]
    var samples: MutPointer[ReferenceSample, MutUntrackedOrigin]
    var bases: ImmPointer[List[Int], ImmUntrackedOrigin]
    """Each sample's paths' first placements."""
    var placements: ImmPointer[OutMap, ImmUntrackedOrigin]
    var firsts: ImmPointer[OutMap, ImmUntrackedOrigin]
    var lengths: ImmPointer[UInt32, ImmUntrackedOrigin]
    var more: MutPointer[List[List[UInt64]], MutUntrackedOrigin]
    var bad: MutPointer[Int, MutUntrackedOrigin]

    def run(self, index: Int):
        var job = self.jobs[unsafe_offset=index]
        var s = job >> 32
        var path = job & 0xFFFF_FFFF
        ref sample = self.samples[unsafe_offset=s]
        var first = self.firsts[unsafe_offset=s].ptr.unsafe_bitcast[Atomic[UInt32]]()
        var walk = List[UInt64]()
        swap(walk, sample.walks[path])
        var n = len(walk)
        var k = self.bases[unsafe_offset=s][path]
        var out = self.placements[unsafe_offset=s].ptr.unsafe_bitcast[UInt64]().unsafe_offset(k)
        var more = List[UInt64]()
        var pos = sample.fragments[path]
        var tag = UInt64(path << 1)
        for step in range(n):
            var v = walk[step]
            var node = Int(v >> 1)
            out[unsafe_offset=step] = UInt64(pos) << UInt64(START_SHIFT) | tag | (v & 1)
            # FIRST holds a placement plus one, 0 for none (a new file's zeros): one less is NONE for none.
            var cur = first[unsafe_offset=node].load()
            var placed = UInt32(k + step + 1)
            var later = placed
            while placed - 1 < cur - 1:
                if first[unsafe_offset=node].compare_exchange(cur, placed):
                    later = cur
                    break
            if later != 0:
                more.append(UInt64(node) << 32 | UInt64(later - 1))
            var length = self.lengths[unsafe_offset=node]
            if length == UInt32(NONE):
                self.bad[] = 2
                return
            pos += Int(length)
        if pos >= 1 << (64 - START_SHIFT):
            self.bad[] = 3
            return
        self.more[unsafe_offset=s][path] = more^
        self.placements[unsafe_offset=s].release(8 * k, 8 * n)


def write_placements(
    mut w: CacheWriter,
    mut samples: List[ReferenceSample],
    lengths: ImmPointer[UInt32, _],
    count: Int,
    threads: Int,
) raises:
    """Each reference sample's sections (sample i's numbered i), its paths' placements laid out all at once (their
    walks freed), then SAMPLES. `lengths`: every node's (NONE for an id the graph lacks)."""
    var n = len(samples)
    var bases = List[List[Int]](capacity=n)
    var placements = List[OutMap](capacity=n)
    var firsts = List[OutMap](capacity=n)
    var more = List[List[List[UInt64]]](capacity=n)
    var jobs = List[Int]()
    var sizes = List[Int]()
    for s in range(n):
        ref sp = samples[s]
        var paths = len(sp.walks)
        if paths >= 1 << PATH_BITS:
            raise Error(t"{paths} paths of reference sample {sp.name}: more than a graph cache holds")
        var b = List[Int](capacity=paths + 1)
        var total = 0
        for i in range(paths):
            b.append(total)
            total += len(sp.walks[i])
            jobs.append(s << 32 | i)
            sizes.append(len(sp.walks[i]))
        b.append(total)
        if total >= Int(NONE):
            raise Error(t"{total} steps on the paths of {sp.name}: more than a graph cache holds")
        var info = List[UInt8]()
        put_u64(info, paths)
        for x in b:
            put_u64(info, x)
        for c in sp.path_contig:
            put_u32(info, c)
        if paths % 2 == 1:
            put_u32(info, 0)
        put_u64(info, len(sp.contigs))
        for c in sp.contigs:
            put_text(info, c)
        w.add(PATHS, s, Span(info))
        placements.append(w.map(w.reserve(PLACEMENTS, s, 8 * total), 8 * total))
        firsts.append(w.map(w.reserve(FIRST, s, 4 * count), 4 * count))
        more.append(List[List[UInt64]](length=paths, fill=List[UInt64]()))
        bases.append(b^)
    var order = List[Int](capacity=len(jobs))
    for i in range(len(jobs)):
        order.append(i)

    def longer(a: Int, b: Int) {sizes} -> Bool:
        return sizes[a] > sizes[b]

    sort(Span(order), longer)
    var sorted = List[Int](capacity=len(jobs))
    for i in order:
        sorted.append(jobs[i])
    var bad = 0
    var task = _PlaceTask(
        imm_untracked(sorted.unsafe_ptr()),
        mut_untracked(samples.unsafe_ptr()),
        imm_untracked(bases.unsafe_ptr()),
        imm_untracked(placements.unsafe_ptr()),
        imm_untracked(firsts.unsafe_ptr()),
        untracked(lengths),
        mut_untracked(more.unsafe_ptr()),
        mut_untracked(Pointer(to=bad)),
    )
    parallel_for(task, len(sorted), threads)
    _ = len(sorted) + len(samples) + len(bases) + len(placements) + len(firsts) + len(more)
    if bad == 2:
        raise Error("a reference path visits a node the graph lacks")
    if bad == 3:
        raise Error("a reference path longer than a graph cache holds")
    if bad != 0:
        raise Error(t"cannot write {w.file.path}")
    firsts.clear()
    for s in range(n):
        # Further placements, by node; a node's in placement order.
        var pairs = List[UInt64]()
        for m in more[s]:
            pairs.extend(Span(m))
        sort(Span(pairs))
        var nodes = List[UInt64]()
        var at = List[UInt64]()
        var rest = List[UInt32](capacity=len(pairs))
        for i in range(len(pairs)):
            var node = pairs[i] >> 32
            if i == 0 or node != pairs[i - 1] >> 32:
                nodes.append(node)
                at.append(UInt64(len(rest)))
            rest.append(UInt32(pairs[i] & 0xFFFF_FFFF))
        at.append(UInt64(len(rest)))
        w.add(MORE_NODES, s, Span(unsafe_ptr=_bytes(nodes.unsafe_ptr()), length=8 * len(nodes)))
        w.add(MORE_AT, s, Span(unsafe_ptr=_bytes(at.unsafe_ptr()), length=8 * len(at)))
        w.add(MORE, s, Span(unsafe_ptr=_bytes(rest.unsafe_ptr()), length=4 * len(rest)))
    var info = List[UInt8]()
    put_u64(info, n)
    for i in range(n):
        put_text(info, samples[i].name)
    w.add(SAMPLES, 0, Span(info))


# From the GBZ.


@fieldwise_init
struct _WalkTask[origin: Origin](Movable, ParallelTask):
    """Path `ids[order[index]]` walked in the GBWT, into `walks[order[index]]`."""

    var graph: Pointer[GbzParts, Self.origin]
    var ids: ImmPointer[Int, ImmUntrackedOrigin]
    var order: ImmPointer[Int, ImmUntrackedOrigin]
    var walks: MutPointer[List[UInt64], MutUntrackedOrigin]
    var bad: MutPointer[Int, MutUntrackedOrigin]

    def run(self, index: Int):
        var j = self.order[unsafe_offset=index]
        try:
            self.walks[unsafe_offset=j] = self.graph[].path(self.ids[unsafe_offset=j])
        except:
            self.bad[] = 1


def make_cache(
    graph_path: String, path: String, samples: List[String], required: Optional[String], width: Int, threads: Int
) raises:
    """The graph's cache at `path`, letters `width` bits each, with the placements of each of `samples` that has paths
    (`required` must)."""
    var timer = _Timer()
    var id = graph_id(graph_path)
    var graph = GbzParts(graph_path)
    timer.lap("graph loaded")
    var md = graph.metadata()
    var names = md.path_names()
    var kept = List[ReferenceSample]()
    for s in samples:
        var sid = md.sample_id(s)
        var sp = ReferenceSample(s)
        if sid:
            var contig_ids = Dict[Int, Int]()
            for i in range(len(names)):
                var p = names[i]
                if p.sample != sid.value():
                    continue
                if p.contig not in contig_ids:
                    contig_ids[p.contig] = len(sp.contigs)
                    sp.contigs.append(md.contig_name(p.contig))
                sp.ids.append(i)
                sp.fragments.append(p.fragment)
                sp.path_contig.append(contig_ids[p.contig])
        if len(sp.ids) == 0:
            if required and required.value() == s:
                raise Error(t"no paths of sample {s} in the graph")
            continue
        kept.append(sp^)
    var ids = List[Int]()
    var firsts = List[Int]()
    var most = 0
    for i in range(len(kept)):
        firsts.append(len(ids))
        ids.extend(Span(kept[i].ids))
        most = max(most, len(kept[i].ids))
    # The samples' paths in turn: a graph's first paths are its longest contigs', each walked on one thread, and
    # every sample's start first.
    var order = List[Int](capacity=len(ids))
    for j in range(most):
        for i in range(len(kept)):
            if j < len(kept[i].ids):
                order.append(firsts[i] + j)
    var walks = List[List[UInt64]](length=len(ids), fill=List[UInt64]())
    var bad = 0
    # The paths are walked on threads of their own while the nodes are found and read.
    var walking = Background(
        _WalkTask(
            Pointer(to=graph),
            imm_untracked(ids.unsafe_ptr()),
            imm_untracked(order.unsafe_ptr()),
            mut_untracked(walks.unsafe_ptr()),
            mut_untracked(Pointer(to=bad)),
        ),
        len(ids),
        threads,
    )
    graph.find_nodes(threads)
    timer.lap("nodes found")
    var count = graph.max_node() + 1
    var w = CacheWriter(path)
    var nodes = NodeWriter(count, width)
    nodes.read_lengths(w, graph, threads)
    timer.lap("node lengths")
    # The letters on threads of their own too, while the placements are laid out.
    var lettering = Background(nodes.letters(Pointer(to=graph)), nodes.workers(threads), threads)
    walking.join()
    _ = len(ids) + len(order)
    if bad != 0:
        raise Error("unreadable reference path")
    timer.lap("reference paths walked")
    var at = 0
    for i in range(len(kept)):
        ref sp = kept[i]
        sp.walks = List[List[UInt64]](capacity=len(sp.ids))
        for _ in range(len(sp.ids)):
            var walk = List[UInt64]()
            swap(walk, walks[at])
            sp.walks.append(walk^)
            at += 1
    write_placements(w, kept, nodes.lengths(), count, threads)
    timer.lap("placements")
    lettering.join()
    _ = graph.max_node()
    nodes.finish(w)
    timer.lap("node letters")
    var slash = graph_path.rfind("/")
    w.finish(id, String(graph_path[byte = slash + 1 :]) if slash >= 0 else graph_path)
    timer.lap("cache written")


# Found, or made again.


def default_cache(graph: String) -> String:
    """Where a graph's cache is kept unless --cache says: beside it, `G.gbz.dromos`."""
    return graph + ".dromos"


def folder_of(path: String) -> String:
    """The folder a path names a file in ("." for a bare name)."""
    var slash = path.rfind("/")
    if slash < 0:
        return String(".")
    return String(path[byte = : max(slash, 1)])


def name_of(path: String) -> String:
    var slash = path.rfind("/")
    return String(path[byte = slash + 1 :]) if slash >= 0 else path


def _c(path: String) -> List[UInt8]:
    var c = List[UInt8](path.as_bytes())
    c.append(0)
    return c^


def _writable(folder: String) -> Bool:
    var c = _c(folder)
    return external_call["access", Int32](c.unsafe_ptr(), Int32(2)) == 0


def graph_cache(graph: String, cache: Optional[String], sample: Optional[String], threads: Int) raises -> String:
    """The path of `graph`'s cache (`cache`, or beside the graph), current and holding `sample`'s placements: made
    (through a temporary name in its folder, then renamed) when missing, of another dromos version, damaged, of an
    older state of the graph or without `sample`. A file that is not a dromos cache, or another graph's cache, is
    never written over."""
    var path = cache.or_else(default_cache(graph))
    if not is_gbz(graph):
        raise Error(t"{graph}: not a GBZ graph")
    var id = graph_id(graph)
    var header = Header()
    var why = examine(path, header)
    var keep = List[String]()
    var reason: String
    if why == CURRENT:
        if header.graph == id:
            if not sample:
                return path
            var held = CacheFile(path).samples.copy()
            for s in held:
                if s == sample.value():
                    return path
            keep = held^
            reason = String(t"adding reference sample {sample.value()}")
        elif header.name != name_of(graph):
            raise Error(t"{path} is the graph cache of {header.name}, not of {graph}: give another --cache")
        else:
            reason = "the graph has changed"
    elif why == NOT_A_CACHE:
        raise Error(t"{path} is not a dromos graph cache: remove it or give another --cache")
    elif why == OTHER_VERSION:
        reason = "made by another dromos version"
    elif why == DAMAGED:
        reason = "damaged"
    else:
        reason = "first use"
    var samples = reference_samples(graph)
    for s in keep:
        if s not in samples:
            samples.append(s)
    if sample and sample.value() not in samples:
        samples.append(sample.value())
    var folder = folder_of(path)
    if not _writable(folder):
        raise Error(
            t"the graph's cache {path} must be made ({reason}) and {folder} is not writable: give --cache FILE in a"
            t" writable folder, or have it made there with dromos prepare --graph {graph}"
        )
    print(t"[dromos] making the graph's cache ({reason}), kept for later runs: {path}", file=stderr)
    var part = path + ".part" + String(external_call["getpid", Int32]())
    var a = _c(part)
    try:
        try:
            make_cache(graph, part, samples, sample, 4, threads)
        except e:
            if String(e) != NO_CODE:
                raise e^
            make_cache(graph, part, samples, sample, 8, threads)
    except e:
        _ = external_call["unlink", Int32](a.unsafe_ptr())
        raise e^
    var b = _c(path)
    if external_call["rename", Int32](a.unsafe_ptr(), b.unsafe_ptr()) != 0:
        _ = external_call["unlink", Int32](a.unsafe_ptr())
        raise Error(t"cannot write {path}")
    return path
