"""The graph files dromos reads alignments against, made from a GBZ: the reference map (`dromos ref-map`, a port of
pgcatalog's refmap.rs, its format in refmap.mojo) and the node table (`dromos node-table`, a port of pgcatalog's
nodetable.rs, its format in nodes.mojo), each byte for byte as pgcatalog writes it. The crate reads the graph; the
files are laid out here, node ranges and paths on all threads, and written at their offsets.
"""

from std.os import getenv
from std.sys import stderr
from std.time import perf_counter_ns

from dromos.graph import GbzParts
from dromos.output import OutFile
from dromos.parallel import ParallelTask, parallel_for
from dromos.refmap import MAGIC as REFMAP_MAGIC, NONE
from dromos.nodes import MAGIC as NODES_MAGIC
from dromos.stream import imm_untracked, mut_untracked

comptime NODE_CHUNK = 1 << 20
"""Node ids per task."""


@always_inline
def _put_u64(mut out: List[UInt8], value: UInt64):
    for i in range(8):
        out.append(UInt8((value >> UInt64(8 * i)) & 0xFF))


@always_inline
def _put_u32(mut out: List[UInt8], value: UInt32):
    for i in range(4):
        out.append(UInt8((value >> UInt32(8 * i)) & 0xFF))


def _put_text(mut out: List[UInt8], text: String):
    _put_u64(out, UInt64(text.byte_length()))
    out.extend(text.as_bytes())


@always_inline
def _bytes[T: AnyType, o: Origin, //](p: Pointer[T, o]) -> Pointer[UInt8, o]:
    """`p`'s memory as bytes, keeping its owner's origin (so the owner lives until the bytes are written)."""
    return p.unsafe_bitcast[UInt8]()


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


@fieldwise_init
struct _LengthsTask[graph_origin: Origin](ParallelTask):
    """Node lengths for ranges of NODE_CHUNK ids, worker w taking ranges w, w + workers, ...: `absent` for an id the
    graph lacks, else the length as u32 (cut to 32 bits, as pgcatalog's `as u32`); each range's total letters, and
    its first id whose length a u32 cannot hold (-1 for none)."""

    var graph: Pointer[GbzParts, Self.graph_origin]
    var count: Int
    var absent: UInt32
    var out: MutPointer[UInt32, MutUntrackedOrigin]
    var sums: MutPointer[Int, MutUntrackedOrigin]
    var too_long: MutPointer[Int, MutUntrackedOrigin]
    var bad: MutPointer[Int, MutUntrackedOrigin]
    var workers: Int

    def run(self, worker: Int):
        var chunks = (self.count + NODE_CHUNK - 1) // NODE_CHUNK
        for index in range(worker, chunks, self.workers):
            var first = index * NODE_CHUNK
            var n = min(NODE_CHUNK, self.count - first)
            try:
                var got = self.graph[].lengths_u32(first, n, self.absent, self.out.unsafe_offset(first))
                self.sums[unsafe_offset=index] = got[0]
                self.too_long[unsafe_offset=index] = got[1]
            except:
                self.bad[] = 1
                return


def _lengths(graph: GbzParts, count: Int, absent: UInt32, threads: Int, mut sums: List[Int]) raises -> List[UInt32]:
    """Every node id's length (`absent` for one the graph lacks), the total letters of each range of NODE_CHUNK ids
    into `sums`; an error at the first id whose length a u32 cannot hold when `absent` is 0 (the node table's)."""
    var out = List[UInt32](unsafe_uninit_length=count)
    var chunks = (count + NODE_CHUNK - 1) // NODE_CHUNK
    sums = List[Int](length=chunks, fill=0)
    var too_long = List[Int](length=chunks, fill=-1)
    var bad = 0
    var task = _LengthsTask(
        Pointer(to=graph),
        count,
        absent,
        mut_untracked(out.unsafe_ptr()),
        mut_untracked(sums.unsafe_ptr()),
        mut_untracked(too_long.unsafe_ptr()),
        mut_untracked(Pointer(to=bad)),
        min(threads, chunks),
    )
    parallel_for(task, min(threads, chunks), threads)
    # The tasks reach these through untracked pointers: reads after them keep them alive until they are done.
    _ = len(sums) + len(too_long)
    if bad != 0:
        raise Error("cannot read the graph's node lengths")
    if absent == 0:
        for long in too_long:
            if long >= 0:
                var n = UInt64(0)
                graph.lengths(long, 1, Pointer(to=n))
                raise Error(t"node {long}: {n} letters, more than a node table holds")
    return out^


# The node table.


@fieldwise_init
struct _LettersTask[graph_origin: Origin](ParallelTask):
    """Ranges of NODE_CHUNK ids' offsets and letters, written at their places, worker w taking ranges w, w + workers,
    ... with its buffers kept."""

    var graph: Pointer[GbzParts, Self.graph_origin]
    var file: ImmPointer[OutFile, ImmUntrackedOrigin]
    var count: Int
    var lengths: ImmPointer[UInt32, ImmUntrackedOrigin]
    var starts: ImmPointer[Int, ImmUntrackedOrigin]
    """Each range's first letter in `bases`, and the end of the last."""
    var offsets_at: Int
    var bases_at: Int
    var bad: MutPointer[Int, MutUntrackedOrigin]
    var workers: Int

    def run(self, worker: Int):
        var offsets = List[UInt64](capacity=NODE_CHUNK + 1)
        var letters = List[UInt8]()
        var f = self.file
        var chunks = (self.count + NODE_CHUNK - 1) // NODE_CHUNK
        for index in range(worker, chunks, self.workers):
            var first = index * NODE_CHUNK
            var n = min(NODE_CHUNK, self.count - first)
            var last = first + n == self.count
            var at = self.starts[unsafe_offset=index]
            offsets.resize(unsafe_uninit_length=n + Int(last))
            for i in range(n):
                offsets[i] = UInt64(at)
                at += Int(self.lengths[unsafe_offset=first + i])
            if last:
                offsets[n] = UInt64(at)
            var ok = f[].write_at(_bytes(offsets.unsafe_ptr()), 8 * len(offsets), self.offsets_at + 8 * first)
            var size = self.starts[unsafe_offset=index + 1] - self.starts[unsafe_offset=index]
            letters.resize(unsafe_uninit_length=size)
            try:
                if self.graph[].letters(first, n, letters.unsafe_ptr(), size) != size:
                    ok = False
            except:
                ok = False
            ok = ok and f[].write_at(letters.unsafe_ptr(), size, self.bases_at + self.starts[unsafe_offset=index])
            if not ok:
                self.bad[] = 1
                return


def write_node_table(graph_path: String, path: String, letters: Bool, threads: Int) raises:
    """`pgcatalog node-table`: the graph's node lengths and, with `letters`, its letters, as one flat file (the
    format in nodes.mojo)."""
    var timer = _Timer()
    var graph = GbzParts(graph_path)
    timer.lap("graph loaded")
    graph.find_nodes(threads)
    timer.lap("nodes found")
    var count = graph.max_node() + 1
    var sums = List[Int]()
    var lengths = _lengths(graph, count, 0, threads, sums)
    timer.lap("node lengths")
    var file = OutFile(path)
    var head = List[UInt8]()
    head.extend(NODES_MAGIC.as_bytes())
    _put_u64(head, UInt64(count))
    file.put(Span(head), 0)
    if not file.write_at(_bytes(lengths.unsafe_ptr()), 4 * count, 16):
        raise Error(t"cannot write {path}")
    var at = 16 + 4 * count
    var tail = List[UInt8]()
    if count % 2 == 1:
        _put_u32(tail, 0)
    _put_u64(tail, UInt64(Int(letters)))
    file.put(Span(tail), at)
    at += len(tail)
    if letters:
        var starts = List[Int](capacity=len(sums) + 1)
        var total = 0
        for s in sums:
            starts.append(total)
            total += s
        starts.append(total)
        var bad = 0
        var task = _LettersTask(
            Pointer(to=graph),
            imm_untracked(Pointer(to=file)),
            count,
            imm_untracked(lengths.unsafe_ptr()),
            imm_untracked(starts.unsafe_ptr()),
            at,
            at + 8 * (count + 1),
            mut_untracked(Pointer(to=bad)),
            min(threads, len(sums)),
        )
        parallel_for(task, min(threads, len(sums)), threads)
        # The tasks reach these through untracked pointers: reads after them keep them alive until they are done.
        _ = len(lengths) + len(starts)
        if bad != 0:
            raise Error(t"cannot write {path}")
    file.close()
    timer.lap("node table written")


# The reference map.


@fieldwise_init
struct _WalkTask[graph_origin: Origin](ParallelTask):
    """Task 0 finds the graph's nodes (on `threads` threads of its own, overlapping the walks, which leave threads
    idle at their end); task i walks reference path i - 1 in the GBWT."""

    var graph: Pointer[GbzParts, Self.graph_origin]
    var ids: ImmPointer[Int, ImmUntrackedOrigin]
    var walks: MutPointer[List[UInt64], MutUntrackedOrigin]
    var bad: MutPointer[Int, MutUntrackedOrigin]
    var threads: Int

    def run(self, index: Int):
        try:
            if index == 0:
                self.graph[].find_nodes(self.threads)
            else:
                self.walks[unsafe_offset=index - 1] = self.graph[].path(self.ids[unsafe_offset=index - 1])
        except:
            self.bad[] = 1


comptime PLACEMENT_BYTES = 17


@fieldwise_init
struct _PlacementsTask(ParallelTask):
    """A reference path's placements (path, step, start, reverse; 17 bytes each), written at their place: its steps'
    positions run from its fragment's offset by the nodes' lengths."""

    var file: ImmPointer[OutFile, ImmUntrackedOrigin]
    var order: ImmPointer[Int, ImmUntrackedOrigin]
    """Paths longest first, so the longest does not start last."""
    var walks: ImmPointer[List[UInt64], ImmUntrackedOrigin]
    var fragments: ImmPointer[Int, ImmUntrackedOrigin]
    var bases: ImmPointer[Int, ImmUntrackedOrigin]
    """Each path's first placement."""
    var node_len: ImmPointer[UInt32, ImmUntrackedOrigin]
    var placements_at: Int
    var bad: MutPointer[Int, MutUntrackedOrigin]

    def run(self, index: Int):
        var path = self.order[unsafe_offset=index]
        ref walk = self.walks[unsafe_offset=path]
        var n = len(walk)
        var out = List[UInt8](unsafe_uninit_length=PLACEMENT_BYTES * n)
        var p = out.unsafe_ptr()
        var pos = UInt64(self.fragments[unsafe_offset=path])
        for step in range(n):
            var v = walk[step]
            var q = p.unsafe_offset(PLACEMENT_BYTES * step)
            q.unsafe_bitcast[UInt32]()[] = UInt32(path)
            q.unsafe_offset(4).unsafe_bitcast[UInt32]()[] = UInt32(step)
            q.unsafe_offset(8).unsafe_bitcast[UInt64]()[] = pos
            q[unsafe_offset=16] = UInt8(v & 1)
            pos += UInt64(self.node_len[unsafe_offset=Int(v >> 1)])
        var at = self.placements_at + PLACEMENT_BYTES * self.bases[unsafe_offset=path]
        if not self.file[].write_at(p, len(out), at):
            self.bad[] = 1


@fieldwise_init
struct _FirstTask(ParallelTask):
    """For the node ids in range `index` (of `span` ids): each one's first placement, and its further ones, listed
    as the map lists them (u64 node, u64 count, then each u32), by node id."""

    var walks: ImmPointer[List[UInt64], ImmUntrackedOrigin]
    var paths: Int
    var size: Int
    var span: Int
    var first: MutPointer[UInt32, MutUntrackedOrigin]
    var more: MutPointer[List[UInt8], MutUntrackedOrigin]
    var more_nodes: MutPointer[Int, MutUntrackedOrigin]

    def run(self, index: Int):
        var lo = index * self.span
        var hi = min(lo + self.span, self.size)
        for node in range(lo, hi):
            self.first[unsafe_offset=node] = UInt32(NONE)
        # (node, k) of each further placement, in placement order.
        var pairs = List[UInt64]()
        var k = 0
        var lo64 = UInt64(lo)
        var width = UInt64(hi - lo)
        for i in range(self.paths):
            ref walk = self.walks[unsafe_offset=i]
            for v in walk:
                var node = v >> 1
                if node - lo64 < width:
                    if self.first[unsafe_offset=Int(node)] == UInt32(NONE):
                        self.first[unsafe_offset=Int(node)] = UInt32(k)
                    else:
                        pairs.append(node << 32 | UInt64(k))
                k += 1
        # By node; a node's placements stay in placement order.
        sort(Span(pairs))
        var out = List[UInt8]()
        var nodes = 0
        var j = 0
        while j < len(pairs):
            var node = pairs[j] >> 32
            var e = j
            while e < len(pairs) and pairs[e] >> 32 == node:
                e += 1
            _put_u64(out, node)
            _put_u64(out, UInt64(e - j))
            for x in range(j, e):
                _put_u32(out, UInt32(pairs[x] & 0xFFFF_FFFF))
            nodes += 1
            j = e
        self.more[unsafe_offset=index] = out^
        self.more_nodes[unsafe_offset=index] = nodes


def write_ref_map(graph_path: String, sample: String, path: String, threads: Int) raises:
    """`pgcatalog ref-map`: where each node lies on `sample`'s paths, saved for dromos (the format in refmap.mojo)."""
    var timer = _Timer()
    var graph = GbzParts(graph_path)
    timer.lap("graph loaded")
    var md = graph.metadata()
    var sample_id = md.sample_id(sample)
    if not sample_id:
        raise Error(t"no paths of sample {sample} in the graph")
    var size = graph.max_node() + 1
    # The sample's paths in path id order, every phase; contigs numbered as first met.
    var ids = List[Int]()
    var fragments = List[Int]()
    var path_contig = List[UInt32]()
    var contigs = List[String]()
    var contig_ids = Dict[Int, Int]()
    var names = md.path_names()
    for id in range(len(names)):
        var p = names[id]
        if p.sample != sample_id.value():
            continue
        if p.contig not in contig_ids:
            contig_ids[p.contig] = len(contigs)
            contigs.append(md.contig_name(p.contig))
        ids.append(id)
        fragments.append(p.fragment)
        path_contig.append(UInt32(contig_ids[p.contig]))
    var paths = len(ids)
    var walks = List[List[UInt64]](length=paths, fill=List[UInt64]())
    var bad = 0
    var walk_task = _WalkTask(
        Pointer(to=graph),
        imm_untracked(ids.unsafe_ptr()),
        mut_untracked(walks.unsafe_ptr()),
        mut_untracked(Pointer(to=bad)),
        threads,
    )
    parallel_for(walk_task, paths + 1, threads)
    _ = len(ids)
    if bad != 0:
        raise Error("unreadable reference path")
    if paths == 0:
        raise Error(t"no paths of sample {sample} in the graph")
    timer.lap("nodes found, reference paths walked")
    var sums = List[Int]()
    var node_len = _lengths(graph, size, UInt32(NONE), threads, sums)
    timer.lap("node lengths")
    var bases = List[Int](capacity=paths)
    var count = 0
    for w in walks:
        bases.append(count)
        count += len(w)
    var order = List[Int](capacity=paths)
    for i in range(paths):
        order.append(i)

    def longer(a: Int, b: Int) {walks} -> Bool:
        return len(walks[a]) > len(walks[b])

    sort(Span(order), longer)
    var file = OutFile(path)
    var head = List[UInt8]()
    head.extend(REFMAP_MAGIC.as_bytes())
    _put_text(head, sample)
    _put_u64(head, UInt64(size))
    file.put(Span(head), 0)
    var node_len_at = len(head)
    if not file.write_at(_bytes(node_len.unsafe_ptr()), 4 * size, node_len_at):
        raise Error(t"cannot write {path}")
    var count_bytes = List[UInt8]()
    _put_u64(count_bytes, UInt64(count))
    var placements_at = node_len_at + 8 * size + 8
    file.put(Span(count_bytes), placements_at - 8)
    var place_task = _PlacementsTask(
        imm_untracked(Pointer(to=file)),
        imm_untracked(order.unsafe_ptr()),
        imm_untracked(walks.unsafe_ptr()),
        imm_untracked(fragments.unsafe_ptr()),
        imm_untracked(bases.unsafe_ptr()),
        imm_untracked(node_len.unsafe_ptr()),
        placements_at,
        mut_untracked(Pointer(to=bad)),
    )
    parallel_for(place_task, paths, threads)
    _ = len(order) + len(fragments) + len(bases) + len(node_len)
    if bad != 0:
        raise Error(t"cannot write {path}")
    timer.lap("placements")
    var first = List[UInt32](unsafe_uninit_length=size)
    var ranges = max(threads, 1)
    var span = (size + ranges - 1) // ranges
    var more = List[List[UInt8]](length=ranges, fill=List[UInt8]())
    var more_nodes = List[Int](length=ranges, fill=0)
    var first_task = _FirstTask(
        imm_untracked(walks.unsafe_ptr()),
        paths,
        size,
        span,
        mut_untracked(first.unsafe_ptr()),
        mut_untracked(more.unsafe_ptr()),
        mut_untracked(more_nodes.unsafe_ptr()),
    )
    parallel_for(first_task, ranges, threads)
    _ = len(walks)
    if not file.write_at(_bytes(first.unsafe_ptr()), 4 * size, node_len_at + 4 * size):
        raise Error(t"cannot write {path}")
    timer.lap("first and further placements")
    var at = placements_at + PLACEMENT_BYTES * count
    var k = 0
    for n in more_nodes:
        k += n
    var tail = List[UInt8]()
    _put_u64(tail, UInt64(k))
    for m in more:
        tail.extend(Span(m))
    _put_u64(tail, UInt64(paths))
    for c in path_contig:
        _put_u32(tail, c)
    _put_u64(tail, UInt64(len(contigs)))
    for c in contigs:
        _put_text(tail, c)
    file.put(Span(tail), at)
    file.close()
    timer.lap("reference map written")
