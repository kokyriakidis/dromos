"""A distance index read against its graph, and its snarls as a table (`dromos snarls --dist`).

dist.mojo reads a `.dist` alone. `check_graph` refuses another graph's: the index must hold the graph's node ids (the
same smallest and largest id and the same number of nodes), each path of the graph's reference samples must start
and end in one of its components (a path never leaves its component), and each step of every top-level chain with no
snarl between its two nodes must be one of the graph's edges. `dist_snarls` gives every component's snarls, oriented
as vg snarls writes them (dist.mojo's notes); `run_dist_snarls` writes them as the five columns of
`dromos snarls --vg-file`, the component being the index's root child (snarlfind.mojo's `write_records_tsv`).
"""

from dromos.dist import CHAIN, DistIndex
from dromos.graph import GbzParts, reference_samples
from dromos.parallel import ParallelTask, parallel_for
from dromos.snarlfind import write_records_tsv
from dromos.snarls import Snarl
from dromos.stream import imm_untracked, mut_untracked, untracked

comptime EDGE_BATCH = 1 << 16
"""Chain steps whose edges are looked up in one call into the crate."""


def check_nodes(d: DistIndex, parts: GbzParts) raises:
    """Refuses a distance index whose node ids are not the graph's: the same smallest and largest id and the same
    number of nodes (`parts` must have found its nodes, `GbzParts.find_nodes`)."""
    var g = parts.node_range()
    var count = d.node_count()
    if d.min_id != g[1] or d.max_id() != g[2] or count != g[0]:
        raise Error(
            t"{d.path} is not this graph's distance index: it holds {count} nodes, ids {d.min_id} to {d.max_id()};"
            t" the graph {g[0]} nodes, ids {g[1]} to {g[2]}"
        )


@fieldwise_init
struct _EndsTask(ParallelTask):
    """Path `ids[index]` walked: its first and last visits (0 for an empty path)."""

    var parts: ImmPointer[GbzParts, ImmUntrackedOrigin]
    var ids: ImmPointer[Int, ImmUntrackedOrigin]
    var firsts: MutPointer[UInt64, MutUntrackedOrigin]
    var lasts: MutPointer[UInt64, MutUntrackedOrigin]
    var errors: MutPointer[String, MutUntrackedOrigin]

    def run(self, index: Int):
        try:
            var w = self.parts[].path(self.ids[unsafe_offset=index])
            if len(w) > 0:
                self.firsts[unsafe_offset=index] = w[0]
                self.lasts[unsafe_offset=index] = w[len(w) - 1]
        except e:
            self.errors[unsafe_offset=index] = String(e)


def check_paths(d: DistIndex, parts: GbzParts, graph: String, threads: Int) raises:
    """Refuses a distance index in which a path of the graph's reference samples (every phase and fragment) starts in
    one component and ends in another, or visits a node the index lacks: the paths walked on `threads` threads."""
    var md = parts.metadata()
    var names = md.path_names()
    var samples = List[Int]()
    for s in reference_samples(graph):
        var id = md.sample_id(s)
        if id:
            samples.append(id.value())
    var ids = List[Int]()
    for id in range(len(names)):
        if names[id].sample in samples:
            ids.append(id)
    var firsts = List[UInt64](length=len(ids), fill=0)
    var lasts = List[UInt64](length=len(ids), fill=0)
    var errors = List[String](length=len(ids), fill=String())
    var task = _EndsTask(
        untracked(Pointer(to=parts)),
        imm_untracked(ids.unsafe_ptr()),
        mut_untracked(firsts.unsafe_ptr()),
        mut_untracked(lasts.unsafe_ptr()),
        mut_untracked(errors.unsafe_ptr()),
    )
    parallel_for(task, len(ids), threads)
    # The tasks reached these through untracked pointers: uses after them keep them alive until they are done.
    _ = len(ids) + len(firsts) + len(lasts) + len(errors)
    for k in range(len(ids)):
        if errors[k]:
            raise Error(t"walking the graph's path {ids[k]}: {errors[k]}")
        if firsts[k] == 0:
            continue
        var p = names[ids[k]]
        var name = String(t"{md.sample_name(p.sample)}#{p.phase}#{md.contig_name(p.contig)}")
        if p.fragment != 0:
            name += String(t"[{p.fragment}]")
        for v in [firsts[k], lasts[k]]:
            if d.node_record(Int(v >> 1)) == 0:
                raise Error(
                    t"{d.path} is not this graph's distance index: the graph's path {name} visits node {v >> 1},"
                    t" which it lacks"
                )
        if d.top(Int(firsts[k] >> 1)) != d.top(Int(lasts[k] >> 1)):
            raise Error(
                t"{d.path} is not this graph's distance index: the graph's path {name} starts in one of its"
                t" components and ends in another"
            )


def _check_edges(
    d: DistIndex,
    parts: GbzParts,
    component: Int,
    mut froms: List[UInt64],
    mut tos: List[UInt64],
    mut succ: List[UInt64],
    mut counts: List[Int64],
) raises:
    """Each `froms[j]` has an edge to `tos[j]`; the lists are emptied."""
    parts.edges(froms.unsafe_ptr(), len(froms), False, succ, counts)
    var at = 0
    for j in range(len(froms)):
        var n = Int(counts[j])
        if n < 0:
            raise Error(t"{d.path} is not this graph's distance index: node {froms[j] >> 1} is not in the graph")
        var found = False
        for x in range(n):
            found = found or succ[at + x] == tos[j]
        at += n
        if not found:
            raise Error(
                t"{d.path} is not this graph's distance index: the chain of its component {component} steps from"
                t" node {froms[j] >> 1} to {tos[j] >> 1}, which the graph does not join"
            )
    froms.clear()
    tos.clear()


def check_chain_edges(d: DistIndex, parts: GbzParts, component: Int) raises:
    """Refuses a distance index whose top-level chain of component `component` (if it is one) steps between nodes the
    graph has no edge between: each pair of consecutive chain nodes with no snarl between them is one of the graph's
    edges."""
    var top = d.root(component)
    if d.kind(top) < CHAIN:
        return
    var words = List[UInt64]()
    var snarl_at = List[Int]()
    d.chain_words(top, words, snarl_at)
    var skip = List[Bool](length=len(words) + 1, fill=False)
    for i in range(len(snarl_at) // 2):
        skip[snarl_at[2 * i + 1]] = True
    var froms = List[UInt64]()
    var tos = List[UInt64]()
    var succ = List[UInt64]()
    var counts = List[Int64]()
    for k in range(1, len(words)):
        if skip[k]:
            continue
        froms.append(words[k - 1])
        tos.append(words[k])
        if len(froms) == EDGE_BATCH:
            _check_edges(d, parts, component, froms, tos, succ, counts)
    if len(froms) > 0:
        _check_edges(d, parts, component, froms, tos, succ, counts)


@fieldwise_init
struct _ChainTask(ParallelTask):
    """Component `index`'s top-level chain checked against the graph's edges."""

    var d: ImmPointer[DistIndex, ImmUntrackedOrigin]
    var parts: ImmPointer[GbzParts, ImmUntrackedOrigin]
    var errors: MutPointer[String, MutUntrackedOrigin]

    def run(self, index: Int):
        try:
            check_chain_edges(self.d[], self.parts[], index)
        except e:
            self.errors[unsafe_offset=index] = String(e)


def check_graph(d: DistIndex, graph: String, threads: Int) raises:
    """Refuses a distance index of another graph than `graph` (the module's notes): its node ids, its reference paths'
    components and its top-level chains' edges, on `threads` threads."""
    var parts = GbzParts(graph)
    parts.find_nodes(threads)
    check_nodes(d, parts)
    check_paths(d, parts, graph, threads)
    var errors = List[String](length=d.components, fill=String())
    var task = _ChainTask(
        untracked(Pointer(to=d)), imm_untracked(Pointer(to=parts)), mut_untracked(errors.unsafe_ptr())
    )
    parallel_for(task, d.components, threads)
    _ = len(errors)
    _ = parts.ptr
    for e in errors:
        if e:
            raise Error(e)


@fieldwise_init
struct _SnarlsTask(ParallelTask):
    """Component `index`'s snarls."""

    var d: ImmPointer[DistIndex, ImmUntrackedOrigin]
    var out: MutPointer[List[Snarl], MutUntrackedOrigin]
    var errors: MutPointer[String, MutUntrackedOrigin]

    def run(self, index: Int):
        try:
            self.d[].snarls(self.d[].root(index), self.out[unsafe_offset=index])
        except e:
            self.errors[unsafe_offset=index] = String(e)


def dist_snarls(d: DistIndex, threads: Int) raises -> List[List[Snarl]]:
    """Every component's snarls (the index's root children, in its order; none for a lone node), oriented as vg snarls
    writes them, on `threads` threads."""
    var out = List[List[Snarl]](capacity=d.components)
    for _ in range(d.components):
        out.append(List[Snarl]())
    var errors = List[String](length=d.components, fill=String())
    var task = _SnarlsTask(
        untracked(Pointer(to=d)), mut_untracked(out.unsafe_ptr()), mut_untracked(errors.unsafe_ptr())
    )
    parallel_for(task, d.components, threads)
    _ = len(out) + len(errors)
    for k in range(len(errors)):
        if errors[k]:
            raise Error(t"component {k}: {errors[k]}")
    return out^


def run_dist_snarls(dist: String, graph: Optional[String], table: String, threads: Int) raises -> Tuple[Int, Int, Int]:
    """`dromos snarls --dist`: the distance index's snarls written to `table` as `dromos snarls --vg-file`'s columns,
    after the index is checked against `graph` when one is given. Returns its version, its components and the
    snarls written."""
    var d = DistIndex(dist)
    if graph:
        check_graph(d, graph.value(), threads)
    var comps = dist_snarls(d, threads)
    write_records_tsv(table, comps)
    var total = 0
    for c in comps:
        total += len(c)
    return (d.version, d.components, total)
