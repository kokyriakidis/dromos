"""vg snarls' snarls of a bidirected graph, found in Mojo: vg's IntegratedSnarlFinder (src/integrated_snarl_finder.cpp,
vg v1.70; MIT) on one weakly connected component, then what vg's SnarlManager and `vg snarls` make of its output
(src/snarls.cpp, src/subcommand/snarls_main.cpp): chains, their orientation (`regularize`) and the snarls written
(those with a node inside).

The finder, in vg's terms (handles "read into" the adjacency component at their end):
1. Adjacency components: an edge a -> b joins a's end and flip(b)'s end (MergedAdjacencyGraph's constructor).
2. Their 3-edge-connected components (Tsin's algorithm as vg implements it, algorithms/
   three_edge_connected_components.cpp) are merged: the cactus graph, whose edges are the graph's nodes.
3. The cactus's cycles (`cycles_in_cactus`): for each cycle edge reading into a vertex, the next edge out of it.
4. The bridge forest (cycles merged) and its longest leaf-to-leaf path (`longest_paths_in_forest`), rooted at the
   longest cycle when that is longer; each forest vertex's edge towards its deepest leaf.
5. The decomposition (`traverse_computed_decomposition`): rooted at that path or cycle; a snarl for each cactus
   vertex on a chain (cycle), bridge paths turned into cycles as they are met, cycles they cross pinched.

Each step is vg's, on flat arrays (ranks, CSR members, a union-find with member rings) in place of vg's hash maps;
vg's iteration orders (hash tables, its union-find's child sets) are not reproduced, so where vg breaks a tie by
them, this does by its own order.

What `vg snarls` (no -T) writes of them (snarls_main.cpp): every snarl with a node inside (`shallow_contents`), in
the orientation `SnarlManager::regularize` gives it: a chain whose snarls mostly run from a higher node id to a lower
is reversed, every snarl with it. vg's `parent` field holds the parent's orientation before `regularize` (it is copied
before `finish`, snarls.cpp:88, and `SnarlManager::flip` changes only the snarl itself, snarls.cpp:829-855), which
depends on vg's hash orders; here a parent is given as written (regularized).

Entry points: `find_snarls` (one component, a `SnarlTree`), `find_all_snarls` (a GBZ's every component, on every
thread), `snarls_of_nodes` (the snarls of the components holding given nodes, as snarls records: the catalog's
source), `write_snarls_tsv` (`dromos snarls`' table). Checked against vg snarls (tests/integration/compare.py,
compare_genome.sh, fuzz.py) and BubbleFinder (bfcheck.py).
"""

from std.os import getenv
from std.sys import stderr
from std.time import perf_counter_ns

from dromos.bigraph import Bigraph, Components, components, read_edges, rss_gb
from dromos.flat import Flat
from dromos.graph import GbzParts
from dromos.output import OutFile
from dromos.parallel import ParallelTask, parallel_for
from dromos.snarls import NO_PARENT, Snarl
from dromos.stream import imm_untracked, mut_untracked

comptime NONE = UInt32.MAX
"""No handle, vertex or snarl."""


# --- union-find -----------------------------------------------------------------------------------------------------


struct _UnionFind(Movable):
    """Union by size, path halving."""

    var parent: Flat[UInt32]
    var size: Flat[UInt32]

    def __init__(out self, n: Int):
        self.parent = Flat[UInt32](unsafe_uninit_length=n)
        self.size = Flat[UInt32](unsafe_uninit_length=n)
        for i in range(n):
            self.parent[i] = UInt32(i)
            self.size[i] = 1

    @always_inline
    def find(mut self, x: UInt32) -> UInt32:
        var i = x
        while True:
            var p = self.parent[Int(i)]
            if p == i:
                return i
            var g = self.parent[Int(p)]
            self.parent[Int(i)] = g
            i = g

    @always_inline
    def union(mut self, a: UInt32, b: UInt32) -> UInt32:
        """Joins a's and b's sets: the new root."""
        var ra = self.find(a)
        var rb = self.find(b)
        if ra == rb:
            return ra
        if self.size[Int(ra)] < self.size[Int(rb)]:
            ra, rb = rb, ra
        self.parent[Int(rb)] = ra
        self.size[Int(ra)] += self.size[Int(rb)]
        return ra


@inline(.never)
def _dense_labels(mut uf: _UnionFind, n: Int, mut labels: Flat[UInt32]) -> Int:
    """Each element's set as a number 0.. in order of the sets' first elements: how many sets."""
    labels.resize(unsafe_uninit_length=n)
    var root_label = Flat[UInt32](length=n, fill=NONE)
    var count = 0
    for i in range(n):
        var r = Int(uf.find(UInt32(i)))
        if root_label[r] == NONE:
            root_label[r] = UInt32(count)
            count += 1
        labels[i] = root_label[r]
    return count


@inline(.never)
def _group(labels: Flat[UInt32], groups: Int, mut starts: Flat[UInt32], mut members: Flat[UInt32]):
    """Elements by label (counting sort, ascending within a label): label g's are members[starts[g] ..< starts[g+1]]."""
    var n = len(labels)
    starts = Flat[UInt32](length=groups + 1, fill=0)
    for i in range(n):
        starts[Int(labels[i]) + 1] += 1
    for g in range(groups):
        starts[g + 1] += starts[g]
    members = Flat[UInt32](unsafe_uninit_length=n)
    var at = Flat[UInt32](unsafe_uninit_length=groups)
    for g in range(groups):
        at[g] = starts[g]
    for i in range(n):
        var g = Int(labels[i])
        members[Int(at[g])] = UInt32(i)
        at[g] += 1


# --- 3-edge-connected components (Tsin, as vg does it) ---------------------------------------------------------------


@inline(.never)
def _three_edge_connected(
    vertices: Int, vstart: Flat[UInt32], vmembers: Flat[UInt32], vid: Flat[UInt32], mut uf: _UnionFind
):
    """vg's `three_edge_connected_component_merges_dense` (Norouzi and Tsin 2014, with vg's handling of bridges and
    self loops) on the multigraph whose vertices are adjacency components and whose edges are the graph's nodes: a
    member handle h of vertex v leads to vid[h ^ 1]; a self loop is followed once, from its forward handle. Merges
    each 3-edge-connected component in `uf`."""
    var dfs_counter = Flat[UInt32](unsafe_uninit_length=vertices)
    var dfs_exit = Flat[UInt32](unsafe_uninit_length=vertices)
    var low_point = Flat[UInt32](unsafe_uninit_length=vertices)
    var degree = Flat[Int64](length=vertices, fill=0)
    var path_tail = Flat[UInt32](unsafe_uninit_length=vertices)
    var on_path = Flat[Bool](unsafe_uninit_length=vertices)
    var visited = Flat[Bool](length=vertices, fill=False)
    # The DFS stack: a vertex, the next member to look at (members are taken last to first, as vg pops its list),
    # whether the tree edge to the parent was seen, whether we came back from a recursion on the current member.
    var st_vertex = Flat[UInt32]()
    var st_cursor = Flat[UInt32]()
    var st_saw_parent = Flat[Bool]()
    var st_recursing = Flat[Bool]()

    @always_inline
    def neighbor(h: UInt32) {vid} -> UInt32:
        return vid[Int(h ^ 1)]

    @always_inline
    def skip(v: UInt32, h: UInt32) {vid} -> Bool:
        # A self loop's reverse handle: vg follows the loop once.
        return vid[Int(h ^ 1)] == v and (h & 1) == 1

    def absorb_path(into_: UInt32, start: UInt32, past_end: UInt32) {mut path_tail, mut on_path, mut degree, mut uf}:
        var into = into_
        var here = start
        while here != past_end:
            if here == NONE:
                break
            if on_path[Int(here)]:
                if into == NONE:
                    into = here
                else:
                    degree[Int(into)] = degree[Int(into)] + degree[Int(here)] - 2
                    _ = uf.union(into, here)
            here = path_tail[Int(here)]

    var counter = UInt32(1)
    for root in range(vertices):
        if visited[root]:
            continue
        st_vertex.append(UInt32(root))
        st_cursor.append(0)
        st_saw_parent.append(False)
        st_recursing.append(False)
        while len(st_vertex) > 0:
            var top = len(st_vertex) - 1
            var cur = st_vertex[top]
            var c = Int(cur)
            if not visited[c]:
                visited[c] = True
                dfs_counter[c] = counter
                counter += 1
                low_point[c] = dfs_counter[c]
                path_tail[c] = NONE
                on_path[c] = True
                st_cursor[top] = vstart[c + 1]
                continue
            # Skip members that are not edges to follow.
            var cursor = st_cursor[top]
            while not st_recursing[top] and cursor > vstart[c] and skip(cur, vmembers[Int(cursor) - 1]):
                cursor -= 1
            st_cursor[top] = cursor
            if cursor > vstart[c]:
                var nb = neighbor(vmembers[Int(cursor) - 1])
                var nbi = Int(nb)
                if not st_recursing[top]:
                    degree[c] += 1
                    if not visited[nbi]:
                        st_recursing[top] = True
                        st_vertex.append(nb)
                        st_cursor.append(0)
                        st_saw_parent.append(False)
                        st_recursing.append(False)
                        continue
                    if top > 0 and nb == st_vertex[top - 1] and not st_saw_parent[top]:
                        st_saw_parent[top] = True
                    elif dfs_counter[nbi] < dfs_counter[c]:
                        if dfs_counter[nbi] < low_point[c]:
                            absorb_path(NONE, cur, NONE)
                            low_point[c] = dfs_counter[nbi]
                            on_path[c] = True
                            path_tail[c] = NONE
                    elif dfs_counter[c] < dfs_counter[nbi]:
                        degree[c] -= 2
                        var replacement = cur
                        var candidate = path_tail[c]
                        while (
                            candidate != NONE
                            and dfs_counter[Int(candidate)] <= dfs_counter[nbi]
                            and dfs_exit[Int(candidate)] >= dfs_exit[nbi]
                        ):
                            replacement = candidate
                            candidate = path_tail[Int(replacement)]
                        absorb_path(NONE, cur, path_tail[Int(replacement)])
                        path_tail[c] = path_tail[Int(replacement)]
                    else:
                        degree[c] -= 1
                    st_cursor[top] = cursor - 1
                else:
                    if low_point[nbi] == dfs_counter[nbi]:
                        # Back across a bridge: hide it.
                        degree[nbi] -= 1
                        degree[c] -= 1
                    else:
                        if degree[nbi] == 2:
                            on_path[nbi] = False
                        if low_point[c] <= low_point[nbi]:
                            absorb_path(cur, nb, NONE)
                        else:
                            low_point[c] = low_point[nbi]
                            absorb_path(NONE, cur, NONE)
                            on_path[c] = True
                            path_tail[c] = nb
                    st_recursing[top] = False
                    st_cursor[top] = cursor - 1
            else:
                dfs_exit[c] = counter
                _ = st_vertex.pop()
                _ = st_cursor.pop()
                _ = st_saw_parent.pop()
                _ = st_recursing.pop()


# --- the result ------------------------------------------------------------------------------------------------------


struct SnarlTree(Movable, Sized):
    """A component's snarls as `vg snarls` writes them, in its order (each top-level snarl, then its subtree: a
    snarl, then its child chains' snarls, the chains last to first, each in order). For snarl i: `start[i]`,
    `end[i]` (visits, `id << 1 | reverse`, regularized), `parent[i]` (its index, NONE at the top level), `chain[i]`
    (its chain's number) and `rank[i]` (its place in the chain). Trivial snarls (no node inside, which vg snarls
    leaves out) are kept in chains' counts (`chain_length`) but not listed."""

    var start: List[UInt64]
    var end: List[UInt64]
    var parent: List[UInt32]
    var chain: List[UInt32]
    var rank: List[UInt32]
    var raw_start: List[UInt64]
    var raw_end: List[UInt64]
    """The snarl as the finder made it, before `regularize` (what vg writes as a child's `parent`)."""
    var chain_length: List[UInt32]
    """Every snarl of the chain, trivial ones included."""
    var chain_tied: List[Bool]
    """The chain's `regularize` vote was tied (as many snarls running down as up): vg keeps its unregularized
    direction, which depends on vg's hash orders; here the finder's is kept."""
    var chain_kind: List[UInt8]
    """Where the chain came from (a diagnostic): 0 a cycle of the cactus graph, 1 a bridge path, 2 a pinched cycle,
    3 the root cycle."""
    var chain_parent: List[UInt32]
    """The chain's parent snarl, NONE at the top level (as a snarl index here, or NONE if the parent is trivial,
    which cannot be)."""
    var stats: List[Int]
    """Per-step nanoseconds and counts, for profiling."""

    def __init__(out self):
        self.start = List[UInt64]()
        self.end = List[UInt64]()
        self.parent = List[UInt32]()
        self.chain = List[UInt32]()
        self.rank = List[UInt32]()
        self.raw_start = List[UInt64]()
        self.raw_end = List[UInt64]()
        self.chain_length = List[UInt32]()
        self.chain_parent = List[UInt32]()
        self.chain_kind = List[UInt8]()
        self.chain_tied = List[Bool]()
        self.stats = List[Int]()

    def __len__(self) -> Int:
        return len(self.start)


struct _Frame(Copyable, TrivialRegisterPassable):
    """A snarl or chain of the decomposition being made: its bounds, the record made for it, where its to-do list
    starts on the shared to-do stack."""

    var is_snarl: Bool
    var first: UInt32
    var second: UInt32
    var base: Int
    var seen: Bool
    var record: UInt32

    def __init__(out self, is_snarl: Bool, first: UInt32, second: UInt32, base: Int):
        self.is_snarl = is_snarl
        self.first = first
        self.second = second
        self.base = base
        self.seen = False
        self.record = NONE


# --- the finder ------------------------------------------------------------------------------------------------------


struct _Finder:
    """One component's decomposition, vg's steps on arrays. Handles are `2 * rank + reverse`."""

    var n: Int
    var handles: Int
    # The cactus graph: each handle's cactus vertex (static), the vertices' members.
    var cv: Flat[UInt32]
    var cactus_vertices: Int
    # Forest vertex of each cactus vertex (static after `_forest`).
    var fv_of_c: Flat[UInt32]
    var forest_vertices: Int
    # Next edge along a cycle, for an edge reading into a cycle vertex (vg's next_along_cycle), NONE if none.
    var nac: Flat[UInt32]
    # Where each cycle edge's `nac` came from (a diagnostic): 0 the cactus's cycles, 1 a bridge path, 2 a pinch.
    var nac_kind: Flat[UInt8]
    # The longest cycle: its length and an edge on it (NONE for none).
    var cycle_length: UInt64
    var cycle_edge: UInt32
    # The longest bridge forest path (edges) and its length; `has_path` false when the cycle won its component.
    var has_path: Bool
    var path_length: UInt64
    var path: Flat[UInt32]
    # Each forest vertex's edge towards its deepest leaf (vg's towards_deepest_leaf), NONE if none.
    var deepest: Flat[UInt32]
    # Scratch of `forest`: the forest vertices' members, and longest_paths_in_forest's DFS records.
    var fstart: Flat[UInt32]
    var fmembers: Flat[UInt32]
    var seen: Flat[Bool]
    var parent_edge: Flat[UInt32]
    var leaf_length: Flat[UInt64]
    var second: Flat[UInt32]
    var has_second: Flat[Bool]
    var sub_root: Flat[UInt32]
    var sub_length: Flat[UInt64]
    # The dynamic cactus union-find (vertices) with member rings over handles.
    var cuf: _UnionFind
    var ring: Flat[UInt32]
    var ring_head: Flat[UInt32]
    # The decomposition: snarls (raw finder orientation, local handles) and chains.
    var s_start: Flat[UInt32]
    var s_end: Flat[UInt32]
    var s_chain: Flat[UInt32]
    var s_next: Flat[UInt32]
    var c_parent: Flat[UInt32]
    var c_first: Flat[UInt32]
    var c_last: Flat[UInt32]
    var c_count: Flat[UInt32]
    var c_kind: Flat[UInt8]
    var visited: Flat[Bool]
    var visited_count: Int
    # Ties vg breaks by its hash orders (diagnostics): other cycles as long as the longest, children as deep as the
    # deepest or second deepest, leaf-leaf paths as long as the longest so far.
    var tie_cycles: Int
    var tie_leaves: Int
    var tie_paths: Int
    # How length ties are broken (a test switch, DROMOS_SNARL_TIES): unset, the first of equals wins (vg's rule, in this
    # code's order); `last`, the last; `mask:N`, the last at the k-th tie met if bit k of N is set.
    var ties_last: Bool
    var tie_mask: Int
    var ties_met: Int

    def __init__(out self, n: Int):
        self.n = n
        self.handles = 2 * n
        self.cv = Flat[UInt32]()
        self.cactus_vertices = 0
        self.fv_of_c = Flat[UInt32]()
        self.forest_vertices = 0
        self.nac = Flat[UInt32]()
        self.nac_kind = Flat[UInt8]()
        self.cycle_length = 0
        self.cycle_edge = NONE
        self.has_path = False
        self.path_length = 0
        self.path = Flat[UInt32]()
        self.deepest = Flat[UInt32]()
        self.fstart = Flat[UInt32]()
        self.fmembers = Flat[UInt32]()
        self.seen = Flat[Bool]()
        self.parent_edge = Flat[UInt32]()
        self.leaf_length = Flat[UInt64]()
        self.second = Flat[UInt32]()
        self.has_second = Flat[Bool]()
        self.sub_root = Flat[UInt32]()
        self.sub_length = Flat[UInt64]()
        self.cuf = _UnionFind(0)
        self.ring = Flat[UInt32]()
        self.ring_head = Flat[UInt32]()
        self.s_start = Flat[UInt32]()
        self.s_end = Flat[UInt32]()
        self.s_chain = Flat[UInt32]()
        self.s_next = Flat[UInt32]()
        self.c_parent = Flat[UInt32]()
        self.c_first = Flat[UInt32]()
        self.c_last = Flat[UInt32]()
        self.c_count = Flat[UInt32]()
        self.c_kind = Flat[UInt8]()
        self.visited = Flat[Bool]()
        self.visited_count = 0
        self.tie_cycles = 0
        self.tie_leaves = 0
        self.tie_paths = 0
        var mode = getenv("DROMOS_SNARL_TIES")
        self.ties_last = mode == "last"
        self.tie_mask = 0
        self.ties_met = 0
        if mode.startswith("mask:"):
            try:
                self.tie_mask = Int(String(mode[byte=5:]))
            except:
                pass

    def tie_last(mut self) -> Bool:
        """At a length tie: whether the later of the equals wins (see `ties_last`)."""
        var k = self.ties_met
        self.ties_met += 1
        return self.ties_last or ((self.tie_mask >> k) & 1) == 1

    @always_inline
    def fv(self, h: UInt32) -> UInt32:
        """The forest vertex handle h reads into."""
        return self.fv_of_c[Int(self.cv[Int(h)])]

    @always_inline
    def cfind(mut self, h: UInt32) -> UInt32:
        """The dynamic cactus vertex (its root) handle h reads into."""
        return self.cuf.find(self.cv[Int(h)])

    def cmerge(mut self, a: UInt32, b: UInt32):
        var ra = self.cfind(a)
        var rb = self.cfind(b)
        if ra == rb:
            return
        var x = self.ring_head[Int(ra)]
        var y = self.ring_head[Int(rb)]
        var t = self.ring[Int(x)]
        self.ring[Int(x)] = self.ring[Int(y)]
        self.ring[Int(y)] = t
        var r = self.cuf.union(ra, rb)
        self.ring_head[Int(r)] = x

    @inline(.never)
    def members(mut self, h: UInt32, mut out: Flat[UInt32]):
        """The handles reading into h's (dynamic) cactus vertex."""
        out.clear()
        var start = self.ring_head[Int(self.cfind(h))]
        var x = start
        while True:
            out.append(x)
            x = self.ring[Int(x)]
            if x == start:
                break

    @always_inline
    def length(self, g: Bigraph, h: UInt32) -> UInt64:
        return UInt64(g.lengths[Int(h >> 1)])

    @always_inline
    def visit_node(mut self, h: UInt32):
        var i = Int(h >> 1)
        if not self.visited[i]:
            self.visited[i] = True
            self.visited_count += 1

    # --- steps 1-2: the cactus graph ---

    @inline(.never)
    def cactus(mut self, g: Bigraph):
        var H = self.handles
        var uf = _UnionFind(H)
        for h in range(H):
            for k in range(Int(g.offsets[h]), Int(g.offsets[h + 1])):
                # An edge a -> b is in a's list and, as flip(b) -> flip(a), in flip(b)'s: both join a and flip(b).
                var other = g.targets[k] ^ 1
                if UInt32(h) <= other:
                    _ = uf.union(UInt32(h), other)
        var vid = Flat[UInt32]()
        var vertices = _dense_labels(uf, H, vid)
        _ = uf^
        var vstart = Flat[UInt32]()
        var vmembers = Flat[UInt32]()
        _group(vid, vertices, vstart, vmembers)
        var tuf = _UnionFind(vertices)
        _three_edge_connected(vertices, vstart, vmembers, vid, tuf)
        _ = vstart^
        _ = vmembers^
        var label = Flat[UInt32]()
        self.cactus_vertices = _dense_labels(tuf, vertices, label)
        self.cv = Flat[UInt32](unsafe_uninit_length=H)
        for h in range(H):
            self.cv[h] = label[Int(vid[h])]

    # --- step 3: cycles in the cactus ---

    @inline(.never)
    def cycles(mut self, g: Bigraph):
        var H = self.handles
        var C = self.cactus_vertices
        var cstart = Flat[UInt32]()
        var cmembers = Flat[UInt32]()
        _group(self.cv, C, cstart, cmembers)
        self.nac = Flat[UInt32](length=H, fill=NONE)
        self.nac_kind = Flat[UInt8](length=H, fill=0)
        var frame_of = Flat[UInt32](length=C, fill=NONE)
        var st_here = Flat[UInt32]()
        var st_cursor = Flat[UInt32]()
        var best = UInt64(0)
        var best_edge = NONE
        for root in range(C):
            if frame_of[root] != NONE:
                continue
            st_here.append(cmembers[Int(cstart[root])])
            st_cursor.append(NONE)
            while len(st_here) > 0:
                var top = len(st_here) - 1
                var here = st_here[top]
                var head = Int(self.cv[Int(here)])
                if st_cursor[top] == NONE:
                    frame_of[head] = UInt32(top)
                    st_cursor[top] = cstart[head + 1]
                # The next member to follow (last to first), not the edge we came in by (unless at the root).
                var cursor = st_cursor[top]
                while cursor > cstart[head] and top > 0 and cmembers[Int(cursor) - 1] == here:
                    cursor -= 1
                if cursor > cstart[head]:
                    st_cursor[top] = cursor - 1
                    var edge_into = cmembers[Int(cursor) - 1] ^ 1
                    var connected = Int(self.cv[Int(edge_into)])
                    if frame_of[connected] == NONE:
                        st_here.append(edge_into)
                        st_cursor.append(NONE)
                    elif UInt32(top) > frame_of[connected]:
                        # A back edge up the stack closes a cycle.
                        var length = self.length(g, edge_into)
                        var prev = edge_into
                        for i in range(Int(frame_of[connected]) + 1, len(st_here)):
                            length += self.length(g, st_here[i])
                            self.nac[Int(prev)] = st_here[i]
                            prev = st_here[i]
                        self.nac[Int(prev)] = edge_into
                        if length > best or (length == best and best > 0 and self.tie_last()):
                            best = length
                            best_edge = edge_into
                            self.tie_cycles = 0
                        elif length == best:
                            self.tie_cycles += 1
                else:
                    _ = st_here.pop()
                    _ = st_cursor.pop()
        self.cycle_length = best
        self.cycle_edge = best_edge

    # --- step 4: the bridge forest and its longest path ---

    @inline(.never)
    def forest(mut self, g: Bigraph):
        var H = self.handles
        var C = self.cactus_vertices
        var fuf = _UnionFind(C)
        for h in range(H):
            var nx = self.nac[h]
            if nx != NONE:
                _ = fuf.union(self.cv[h], self.cv[Int(nx)])
        self.forest_vertices = _dense_labels(fuf, C, self.fv_of_c)
        _ = fuf^
        var F = self.forest_vertices
        var flabel = Flat[UInt32](unsafe_uninit_length=H)
        for h in range(H):
            flabel[h] = self.fv(UInt32(h))
        _group(flabel, F, self.fstart, self.fmembers)
        _ = flabel^
        self.deepest = Flat[UInt32](length=F, fill=NONE)
        self.seen = Flat[Bool](length=F, fill=False)
        self.parent_edge = Flat[UInt32](unsafe_uninit_length=F)
        self.leaf_length = Flat[UInt64](length=F, fill=0)
        self.second = Flat[UInt32](length=F, fill=NONE)
        self.has_second = Flat[Bool](length=F, fill=False)
        self.sub_root = Flat[UInt32](unsafe_uninit_length=F)
        self.sub_length = Flat[UInt64](length=F, fill=0)
        if self.cycle_edge != NONE:
            self.try_root(g, self.cycle_edge, self.cycle_length)
        for f in range(F):
            self.try_root(g, self.fmembers[Int(self.fstart[f])], 0)
        # Only `deepest` (and the path) outlive this step.
        self.seen = Flat[Bool]()
        self.parent_edge = Flat[UInt32]()
        self.leaf_length = Flat[UInt64]()
        self.second = Flat[UInt32]()
        self.has_second = Flat[Bool]()
        self.sub_root = Flat[UInt32]()
        self.sub_length = Flat[UInt64]()
        self.fstart = Flat[UInt32]()
        self.fmembers = Flat[UInt32]()

    def mix(mut self, parent: UInt32, child_edge: UInt32, child_length: UInt64):
        """Child edge `child_edge` (its subtree's leaf path `child_length`) into parent's deepest and second deepest
        children, as vg does it."""
        var p = Int(parent)
        var pd = self.deepest[p]
        if pd != NONE and (
            self.leaf_length[Int(self.fv(pd))] == child_length
            or (self.has_second[p] and self.leaf_length[Int(self.fv(self.second[p]))] == child_length)
        ):
            self.tie_leaves += 1
        if pd == NONE:
            self.deepest[p] = child_edge
        elif self.leaf_length[Int(self.fv(pd))] < child_length or (
            self.leaf_length[Int(self.fv(pd))] == child_length and self.tie_last()
        ):
            self.second[p] = pd
            self.has_second[p] = True
            self.deepest[p] = child_edge
        elif not self.has_second[p]:
            self.second[p] = child_edge
            self.has_second[p] = True
        elif self.leaf_length[Int(self.fv(self.second[p]))] < child_length or (
            self.leaf_length[Int(self.fv(self.second[p]))] == child_length and self.tie_last()
        ):
            self.second[p] = child_edge

    @inline(.never)
    def try_root(mut self, g: Bigraph, traversal_root: UInt32, root_cycle_length: UInt64):
        """vg's `try_root` in `longest_paths_in_forest`: a DFS of the bridge tree from `traversal_root`'s forest
        vertex, each vertex's deepest and second deepest child, the longest leaf-leaf path, which (if at least
        `root_cycle_length`) is recorded and becomes the tree's root."""
        if self.seen[Int(self.fv(traversal_root))]:
            return
        var st_here = Flat[UInt32]()
        var st_cursor = Flat[UInt32]()
        st_here.append(traversal_root)
        st_cursor.append(NONE)
        while len(st_here) > 0:
            var top = len(st_here) - 1
            var here = st_here[top]
            var fh = self.fv(here)
            var f = Int(fh)
            if st_cursor[top] == NONE:
                self.seen[f] = True
                self.parent_edge[f] = here ^ 1
                self.sub_root[f] = fh
                st_cursor[top] = self.fstart[f + 1]
            var cursor = st_cursor[top]
            # Bridges only: members whose other side is in another forest vertex.
            while cursor > self.fstart[f] and self.fv(self.fmembers[Int(cursor) - 1] ^ 1) == fh:
                cursor -= 1
            if cursor > self.fstart[f]:
                st_cursor[top] = cursor - 1
                var edge_into = self.fmembers[Int(cursor) - 1] ^ 1
                if not self.seen[Int(self.fv(edge_into))]:
                    st_here.append(edge_into)
                    st_cursor.append(NONE)
                continue
            # Done with every child.
            var d = self.deepest[f]
            if top > 0:
                var ph = self.fv(st_here[top - 1])
                self.leaf_length[f] = self.length(g, here)
                if d != NONE:
                    self.leaf_length[f] += self.leaf_length[Int(self.fv(d))]
                self.mix(ph, here, self.leaf_length[f])
            if self.has_second[f] or top == 0:
                var here_length = UInt64(0)
                if d != NONE:
                    here_length += self.leaf_length[Int(self.fv(d))]
                if self.has_second[f]:
                    here_length += self.leaf_length[Int(self.fv(self.second[f]))]
                if self.sub_root[f] != fh and here_length == self.sub_length[f] and here_length > 0:
                    self.tie_paths += 1
                if (
                    self.sub_root[f] == fh
                    or here_length > self.sub_length[f]
                    or (here_length == self.sub_length[f] and here_length > 0 and self.tie_last())
                ):
                    self.sub_root[f] = fh
                    self.sub_length[f] = here_length
            if top > 0 and self.sub_length[f] > 0:
                var p = Int(self.fv(st_here[top - 1]))
                if self.sub_root[p] != UInt32(p) and self.sub_length[p] == self.sub_length[f]:
                    self.tie_paths += 1
                if (
                    self.sub_root[p] == UInt32(p)
                    or self.sub_length[p] < self.sub_length[f]
                    or (self.sub_length[p] == self.sub_length[f] and self.tie_last())
                ):
                    self.sub_root[p] = self.sub_root[f]
                    self.sub_length[p] = self.sub_length[f]
            if top == 0 and self.sub_length[f] >= root_cycle_length:
                self.record_path(g, fh, traversal_root)
            _ = st_here.pop()
            _ = st_cursor.pop()

    def record_path(mut self, g: Bigraph, root: UInt32, traversal_root: UInt32):
        """The tree's longest leaf-leaf path wins: record it and root the tree at its convergence point."""
        self.has_path = True
        self.path_length = self.sub_length[Int(root)]
        self.path.clear()
        var conv = self.sub_root[Int(root)]
        var c = Int(conv)
        if self.has_second[c]:
            var down = Flat[UInt32]()
            down.append(self.second[c])
            var x = self.deepest[Int(self.fv(self.second[c]))]
            while x != NONE:
                down.append(x)
                x = self.deepest[Int(self.fv(x))]
            for k in range(len(down) - 1, -1, -1):
                self.path.append(down[k] ^ 1)
        if self.deepest[c] != NONE:
            self.path.append(self.deepest[c])
            var x = self.deepest[Int(self.fv(self.deepest[c]))]
            while x != NONE:
                self.path.append(x)
                x = self.deepest[Int(self.fv(x))]
        # Re-root along the path from the convergence point to the traversal root.
        var up = Flat[UInt32]()
        var at = conv
        while at != root:
            var pe = self.parent_edge[Int(at)]
            up.append(pe)
            at = self.fv(pe)
        while len(up) > 0:
            var pce = up.pop()
            var child = Int(self.fv(pce))
            var parent = self.fv(pce ^ 1)
            var cd = self.deepest[child]
            if cd != NONE and self.fv(cd) == parent:
                if self.has_second[child]:
                    self.deepest[child] = self.second[child]
                    self.has_second[child] = False
                else:
                    self.deepest[child] = NONE
            self.leaf_length[child] = self.length(g, pce)
            if self.deepest[child] != NONE:
                self.leaf_length[child] += self.leaf_length[Int(self.fv(self.deepest[child]))]
            self.mix(parent, pce, self.leaf_length[child])
        if len(self.path) == 0:
            self.path.append(traversal_root)

    # --- step 5: the decomposition ---

    def new_chain(mut self, parent: UInt32) -> UInt32:
        var c = UInt32(len(self.c_parent))
        self.c_parent.append(parent)
        self.c_first.append(NONE)
        self.c_last.append(NONE)
        self.c_count.append(0)
        return c

    def new_snarl(mut self, start: UInt32, end: UInt32, chain: UInt32) -> UInt32:
        var s = UInt32(len(self.s_start))
        self.s_start.append(start)
        self.s_end.append(end)
        self.s_chain.append(chain)
        self.s_next.append(NONE)
        var c = Int(chain)
        if self.c_last[c] == NONE:
            self.c_first[c] = s
        else:
            self.s_next[Int(self.c_last[c])] = s
        self.c_last[c] = s
        self.c_count[c] += 1
        return s

    def cycle_path(mut self, start: UInt32, end_root: UInt32) raises -> Flat[UInt32]:
        """vg's `find_cycle_path_in_cactus`: the cycles (each by its edge into the vertex it shares with the one
        before) from start's cactus vertex to one that touches `end_root`."""
        var roots = Flat[UInt32]()
        var mem = Flat[UInt32]()
        self.members(start, mem)
        for i_m in range(len(mem)):
            var m = mem[i_m]
            if self.nac[Int(m)] != NONE:
                roots.append(m)
        var st_edge = Flat[UInt32]()
        var st_base = Flat[Int]()
        var st_expanded = Flat[Bool]()
        var children = Flat[UInt32]()
        for i_root in range(len(roots)):
            var root = roots[i_root]
            st_edge.append(root)
            st_base.append(len(children))
            st_expanded.append(False)
            while len(st_edge) > 0:
                var top = len(st_edge) - 1
                if not st_expanded[top]:
                    st_expanded[top] = True
                    var entry = st_edge[top]
                    var e = self.nac[Int(entry)]
                    while e != entry:
                        if self.cfind(e) == end_root:
                            return st_edge.copy()
                        self.members(e, mem)
                        for i_m in range(len(mem)):
                            var m = mem[i_m]
                            if m != e and self.nac[Int(m)] != NONE:
                                children.append(m)
                        e = self.nac[Int(e)]
                if len(children) > st_base[top]:
                    var child = children.pop()
                    st_edge.append(child)
                    st_base.append(len(children))
                    st_expanded.append(False)
                else:
                    _ = st_edge.pop()
                    _ = st_base.pop()
                    _ = st_expanded.pop()
        raise Error("snarlfind: could not find a cycle path in the cactus graph")

    @inline(.never)
    def decompose(mut self, g: Bigraph) raises:
        var H = self.handles
        var C = self.cactus_vertices
        self.cuf = _UnionFind(C)
        # Member rings, one per cactus vertex.
        self.ring = Flat[UInt32](unsafe_uninit_length=H)
        self.ring_head = Flat[UInt32](length=C, fill=NONE)
        var last = Flat[UInt32](length=C, fill=NONE)
        for h in range(H):
            var c = Int(self.cv[h])
            if self.ring_head[c] == NONE:
                self.ring_head[c] = UInt32(h)
            else:
                self.ring[Int(last[c])] = UInt32(h)
            last[c] = UInt32(h)
        for c in range(C):
            self.ring[Int(last[c])] = self.ring_head[c]
        _ = last^
        self.visited = Flat[Bool](length=self.n, fill=False)
        self.visited_count = 0
        var frames = List[_Frame]()
        var todo = Flat[UInt32]()
        var mem = Flat[UInt32]()

        var rooted = False
        while self.visited_count < self.n:
            if rooted:
                raise Error("snarlfind: nodes left out of the decomposition")
            rooted = True
            if self.cycle_edge == NONE or (self.has_path and self.cycle_length <= self.path_length):
                if not self.has_path:
                    raise Error("snarlfind: no root path or cycle")
                var first_edge = self.path[0]
                if self.visited[Int(first_edge >> 1)]:
                    continue
                if self.path_length == 0:
                    self.members(first_edge, mem)
                    for i_m in range(len(mem)):
                        var m = mem[i_m]
                        if (m & 1) == 0:
                            self.visit_node(m)
                else:
                    for i in range(1, len(self.path)):
                        self.deepest[Int(self.fv(self.path[i - 1]))] = self.path[i]
                    frames.append(_Frame(True, NONE, NONE, len(todo)))
                    todo.append(first_edge ^ 1)
                    self.members(first_edge ^ 1, mem)
                    for i_m in range(len(mem)):
                        var m = mem[i_m]
                        if m == first_edge ^ 1:
                            continue
                        if self.nac[Int(m)] != NONE:
                            todo.append(m)
                        elif (m & 1) == 0 and self.cfind(m ^ 1) == self.cfind(m):
                            self.visit_node(m)
            else:
                if self.visited[Int(self.cycle_edge >> 1)]:
                    continue
                frames.append(_Frame(True, NONE, NONE, len(todo)))
                frames.append(_Frame(False, self.cycle_edge, self.cycle_edge, len(todo)))

            while len(frames) > 0:
                var top = len(frames) - 1
                if top > 0 and not frames[top].seen:
                    frames[top].seen = True
                    var first = frames[top].first
                    var second = frames[top].second
                    if frames[top].is_snarl:
                        frames[top].record = self.new_snarl(first, second, frames[top - 1].record)
                        self.visit_node(first)
                        self.visit_node(second)
                        self.members(first, mem)
                        for i_m in range(len(mem)):
                            var m = mem[i_m]
                            if m == first or (m ^ 1) == second:
                                continue
                            if self.fv(m ^ 1) != self.fv(m):
                                todo.append(m)
                            elif self.nac[Int(m)] != NONE:
                                todo.append(m)
                            elif (m & 1) == 0 and self.cfind(m ^ 1) == self.cfind(m):
                                self.visit_node(m)
                    else:
                        var parent = frames[top - 1].record if top > 1 else NONE
                        frames[top].record = self.new_chain(parent)
                        self.c_kind.append(3 if first == second else self.nac_kind[Int(second)])
                        var region = len(todo)
                        var here = first
                        while True:
                            todo.append(here)
                            here = self.nac[Int(here)]
                            if here == second:
                                break
                        # Forward along the chain when popped.
                        var a = region
                        var b = len(todo) - 1
                        while a < b:
                            var t = todo[a]
                            todo[a] = todo[b]
                            todo[b] = t
                            a += 1
                            b -= 1
                if len(todo) > frames[top].base:
                    var task = todo.pop()
                    if frames[top].is_snarl:
                        var outgoing = self.nac[Int(task)]
                        if outgoing != NONE:
                            frames.append(_Frame(False, outgoing, task, len(todo)))
                        else:
                            self.bridge(task, todo, mem)
                    else:
                        frames.append(_Frame(True, task, self.nac[Int(task)], len(todo)))
                else:
                    _ = frames.pop()

    @inline(.never)
    def bridge(mut self, task: UInt32, mut todo: Flat[UInt32], mut mem: Flat[UInt32]) raises:
        """A bridge edge reading into the current snarl: its path to the deepest leaf becomes a cycle through the
        snarl (cycles it crosses pinched), and the leaf's cactus vertex joins the snarl's."""
        var edge = task ^ 1
        var head = self.cfind(edge)
        var deep = self.deepest[Int(self.fv(edge))]
        while deep != NONE:
            var next_back = self.cfind(deep ^ 1)
            if head != next_back:
                var cpath = self.cycle_path(edge, next_back)
                while len(cpath) > 0:
                    var through = cpath[len(cpath) - 1]
                    var end = through
                    while True:
                        end = self.nac[Int(end)]
                        if self.cfind(end) == self.cfind(deep ^ 1):
                            break
                    self.cmerge(through, deep ^ 1)
                    var t = self.nac[Int(through)]
                    self.nac[Int(through)] = self.nac[Int(end)]
                    self.nac[Int(end)] = t
                    self.nac_kind[Int(through)] = 2
                    self.nac_kind[Int(end)] = 2
                    if self.nac[Int(through)] == through:
                        self.nac[Int(through)] = NONE
                    if self.nac[Int(end)] == end:
                        self.nac[Int(end)] = NONE
                    _ = cpath.pop()
            self.nac[Int(edge)] = deep
            self.nac_kind[Int(edge)] = 1
            edge = deep
            head = self.cfind(edge)
            deep = self.deepest[Int(self.fv(edge))]
        if edge == task ^ 1:
            # A single bridge: a node held in the snarl, no chain.
            self.visit_node(edge)
        else:
            self.nac[Int(edge)] = task ^ 1
            self.nac_kind[Int(edge)] = 1
        self.members(edge, mem)
        for i_m in range(len(mem)):
            var m = mem[i_m]
            if self.nac[Int(m)] != NONE:
                todo.append(m)
            elif (m & 1) == 0 and self.cfind(m ^ 1) == self.cfind(m):
                self.visit_node(m)
        self.cmerge(edge, task)


# --- what vg snarls writes -------------------------------------------------------------------------------------------


@inline(.never)
def _trivial(g: Bigraph, start: UInt32, end: UInt32) -> Bool:
    """vg's `shallow_contents` finds no node: every neighbor inward of the start and of the end is a boundary node."""
    var a = start >> 1
    var b = end >> 1
    for k in range(Int(g.offsets[Int(start)]), Int(g.offsets[Int(start) + 1])):
        var x = g.targets[k] >> 1
        if x != a and x != b:
            return False
    for k in range(Int(g.offsets[Int(end ^ 1)]), Int(g.offsets[Int(end ^ 1) + 1])):
        var x = g.targets[k] >> 1
        if x != a and x != b:
            return False
    return True


def find_snarls(g: Bigraph) raises -> SnarlTree:
    """The snarls `vg snarls` writes for the component `g` (see the module notes)."""
    var out = SnarlTree()
    var n = g.nodes()
    if n == 0:
        return out^
    var t = perf_counter_ns()
    var f = _Finder(n)
    f.cactus(g)
    var t1 = perf_counter_ns()
    f.cycles(g)
    var t2 = perf_counter_ns()
    f.forest(g)
    var t3 = perf_counter_ns()
    f.decompose(g)
    var t4 = perf_counter_ns()
    # Regularize: per chain, a vote over all its snarls (trivial ones too) in the finder's orientation.
    var S = len(f.s_start)
    var nchains = len(f.c_parent)
    var flip = Flat[Bool](length=nchains, fill=False)
    out.chain_tied = List[Bool](length=nchains, fill=False)
    for c in range(nchains):
        var ascending = 0
        var s = f.c_first[c]
        while s != NONE:
            if g.ids[Int(f.s_start[Int(s)] >> 1)] <= g.ids[Int(f.s_end[Int(s)] >> 1)]:
                ascending += 1
            s = f.s_next[Int(s)]
        flip[c] = 2 * ascending < Int(f.c_count[c])
        out.chain_tied[c] = 2 * ascending == Int(f.c_count[c])
    # Chain members in final order, and each snarl's children chains.
    var members = Flat[UInt32](unsafe_uninit_length=S)
    var mstart = Flat[UInt32](unsafe_uninit_length=nchains + 1)
    var at = 0
    for c in range(nchains):
        mstart[c] = UInt32(at)
        var s = f.c_first[c]
        while s != NONE:
            members[at] = s
            at += 1
            s = f.s_next[Int(s)]
        if flip[c]:
            var a = Int(mstart[c])
            var b = at - 1
            while a < b:
                var x = members[a]
                members[a] = members[b]
                members[b] = x
                a += 1
                b -= 1
    mstart[nchains] = UInt32(at)
    var rank_in = Flat[UInt32](unsafe_uninit_length=S)
    for c in range(nchains):
        for k in range(Int(mstart[c]), Int(mstart[c + 1])):
            rank_in[Int(members[k])] = UInt32(k - Int(mstart[c]))
    var child_start = Flat[UInt32](length=S + 1, fill=0)
    for c in range(nchains):
        var p = f.c_parent[c]
        if p != NONE:
            child_start[Int(p) + 1] += 1
    for s in range(S):
        child_start[s + 1] += child_start[s]
    var child_chains = Flat[UInt32](unsafe_uninit_length=len(f.c_parent))
    var fill = child_start.copy()
    var roots = Flat[UInt32]()
    for c in range(nchains):
        var p = f.c_parent[c]
        if p == NONE:
            roots.append(UInt32(c))
        else:
            child_chains[Int(fill[Int(p)])] = UInt32(c)
            fill[Int(p)] += 1
    # Written order and indices of the non-trivial snarls.
    var trivial = Flat[Bool](unsafe_uninit_length=S)
    for s in range(S):
        trivial[s] = _trivial(g, f.s_start[s], f.s_end[s])
    var index = Flat[UInt32](length=S, fill=NONE)
    out.chain_length = List[UInt32](unsafe_uninit_length=nchains)
    for c in range(nchains):
        out.chain_length[c] = f.c_count[c]
    var stack = Flat[UInt32]()
    for i_rc in range(len(roots)):
        var rc = roots[i_rc]
        for k in range(Int(mstart[Int(rc)]), Int(mstart[Int(rc) + 1])):
            stack.append(members[k])
            while len(stack) > 0:
                var s = Int(stack.pop())
                if trivial[s]:
                    continue
                var a = f.s_start[s]
                var b = f.s_end[s]
                out.raw_start.append(g.visit(a))
                out.raw_end.append(g.visit(b))
                var c = Int(f.s_chain[s])
                if flip[c]:
                    var x = a
                    a = b ^ 1
                    b = x ^ 1
                index[s] = UInt32(len(out.start))
                out.start.append(g.visit(a))
                out.end.append(g.visit(b))
                var p = f.c_parent[c]
                out.parent.append(index[Int(p)] if p != NONE else NONE)
                out.chain.append(UInt32(c))
                out.rank.append(rank_in[s])
                for j in range(Int(child_start[s]), Int(child_start[s + 1])):
                    var cc = Int(child_chains[j])
                    for k2 in range(Int(mstart[cc + 1]) - 1, Int(mstart[cc]) - 1, -1):
                        stack.append(members[k2])
    out.chain_kind = f.c_kind.data.copy()
    out.chain_parent = List[UInt32](unsafe_uninit_length=nchains)
    for c in range(nchains):
        var p = f.c_parent[c]
        out.chain_parent[c] = index[Int(p)] if p != NONE else NONE
    var t5 = perf_counter_ns()
    out.stats = [
        t1 - t,
        t2 - t1,
        t3 - t2,
        t4 - t3,
        t5 - t4,
        S,
        nchains,
        f.cactus_vertices,
        f.forest_vertices,
        f.tie_cycles,
        f.tie_leaves,
        f.tie_paths,
    ]
    return out^


def as_snarls(t: SnarlTree) -> List[Snarl]:
    """The snarls as dromos's snarls records (`dromos.snarls.Snarl`, what a vg snarls file reads as): start, end and
    the parent's start and end (NO_PARENT at the top level), in vg snarls' order. The parent is given as it is
    written (regularized), where vg's file has its unregularized orientation (see the module notes): read it as a
    snarl in either direction (`dromos.snarls.canonical`)."""
    var out = List[Snarl](capacity=len(t))
    for i in range(len(t)):
        var p = t.parent[i]
        if p == NONE:
            out.append(Snarl(t.start[i], t.end[i], NO_PARENT, NO_PARENT))
        else:
            out.append(Snarl(t.start[i], t.end[i], t.start[Int(p)], t.end[Int(p)]))
    return out^


comptime TIES = (
    "vg IntegratedSnarlFinder's; where vg breaks a tie by its hash tables' order (cycles or leaf-to-leaf paths of"
    " equal length, a chain's tied orientation vote), the first in dromos's own order wins (the same on any threads)"
)
"""How ties vg leaves to its hash orders are broken (a catalog header's `snarls.ties`)."""


@fieldwise_init
struct _FindTask(ParallelTask):
    """Component `order[index]`'s snarls, as a SnarlTree or (`records`) as snarls records; its graph is freed once
    they are found."""

    var comps: MutPointer[Bigraph, MutUntrackedOrigin]
    var order: ImmPointer[Int, ImmUntrackedOrigin]
    var trees: MutPointer[SnarlTree, MutUntrackedOrigin]
    var records: MutPointer[List[Snarl], MutUntrackedOrigin]
    var as_records: Bool
    var errors: MutPointer[String, MutUntrackedOrigin]

    def run(self, index: Int):
        var k = self.order[unsafe_offset=index]
        var slot = self.comps.unsafe_offset(k)
        var g = slot.unsafe_take_pointee()
        slot.unsafe_write(Bigraph())
        try:
            var t = find_snarls(g)
            if self.as_records:
                self.records[unsafe_offset=k] = as_snarls(t)
            else:
                self.trees[unsafe_offset=k] = t^
        except e:
            self.errors[unsafe_offset=k] = String(t"component {k}: {e}")


def _find(
    mut comps: Components, wanted: List[Int], threads: Int, as_records: Bool
) raises -> Tuple[List[SnarlTree], List[List[Snarl]]]:
    """The snarls of the components `wanted`, the biggest first on `threads` threads, each component's graph freed
    once they are found: per component, a tree (or, `as_records`, its snarls records); empty for those not wanted."""
    var n = len(comps.graphs)
    var order = List[Int]()
    var seen = List[Bool](length=n, fill=False)
    for k in wanted:
        if not seen[k]:
            seen[k] = True
            order.append(k)
    ref graphs = comps.graphs

    def bigger(a: Int, b: Int) {graphs} -> Bool:
        return graphs[a].nodes() > graphs[b].nodes() or (graphs[a].nodes() == graphs[b].nodes() and a < b)

    sort(Span(order), bigger)
    var trees = List[SnarlTree](capacity=n if not as_records else 0)
    var records = List[List[Snarl]](capacity=n if as_records else 0)
    for _ in range(n):
        if as_records:
            records.append(List[Snarl]())
        else:
            trees.append(SnarlTree())
    var errors = List[String](length=n, fill=String())
    var task = _FindTask(
        mut_untracked(comps.graphs.unsafe_ptr()),
        imm_untracked(order.unsafe_ptr()),
        mut_untracked(trees.unsafe_ptr()),
        mut_untracked(records.unsafe_ptr()),
        as_records,
        mut_untracked(errors.unsafe_ptr()),
    )
    parallel_for(task, len(order), threads)
    # The tasks reached these through untracked pointers: uses after them keep them alive until they are done.
    _ = len(comps.graphs) + len(order) + len(trees) + len(records)
    for e in errors:
        if e:
            raise Error(t"finding snarls: {e}")
    return (trees^, records^)


def _log(verbose: Bool, what: StaticString, t0: Int) -> Int:
    var t = perf_counter_ns()
    if verbose:
        print(t"[snarls] {what}: {Float64(t - t0) / 1e9} s, rss {rss_gb()} GB", file=stderr)
    return t


def find_all_snarls(path: String, threads: Int, verbose: Bool = False) raises -> List[SnarlTree]:
    """Every weakly connected component's snarls of the GBZ graph at `path` (a tree per component, in order of their
    smallest node id), as `vg snarls` writes them, on `threads` threads. The GBZ is dropped once its edges are read."""
    var t = perf_counter_ns()
    var edges = read_edges(GbzParts(path), threads, verbose)
    t = _log(verbose, "graph loaded and edges read", t)
    var comps = components(edges^, threads, verbose)
    comps.label = Flat[UInt32]()
    t = _log(verbose, "components made", t)
    var wanted = List[Int](capacity=len(comps))
    for k in range(len(comps)):
        wanted.append(k)
    var found = _find(comps, wanted, threads, False)
    _ = _log(verbose, "snarls found", t)
    var trees = List[SnarlTree]()
    swap(trees, found[0])
    return trees^


def snarls_of_nodes(
    parts: GbzParts, groups: List[List[Int]], threads: Int, verbose: Bool = False
) raises -> List[List[Snarl]]:
    """For each group of node ids, the snarls of the graph's components that hold them (each such component's once,
    the components in order of their smallest node id), as `vg snarls` writes them for those components
    (`as_snarls`). Only those components' snarls are found. A node the graph lacks is refused."""
    var t = perf_counter_ns()
    var comps = components(read_edges(parts, threads, verbose), threads, verbose)
    t = _log(verbose, "edges read, components made", t)
    var of_group = List[List[Int]]()
    var wanted = List[Int]()
    for group in groups:
        var ks = List[Int]()
        for id in group:
            var k = comps.of(id)
            if not k:
                raise Error(t"node {id} is not in the graph")
            var known = False
            for x in ks:
                known = known or x == k.value()
            if not known:
                ks.append(k.value())
                wanted.append(k.value())
        sort(Span(ks))
        of_group.append(ks^)
    comps.label = Flat[UInt32]()
    var found = _find(comps, wanted, threads, True)
    _ = _log(verbose, "snarls found", t)
    ref records = found[1]
    var out = List[List[Snarl]](capacity=len(groups))
    for ks in of_group:
        var total = 0
        for k in ks:
            total += len(records[k])
        var snarls = List[Snarl](capacity=total)
        for k in ks:
            snarls.extend(Span(records[k]))
        out.append(snarls^)
    return out^


def _flush(mut f: OutFile, mut buf: String, mut at: Int) raises:
    f.put(buf.as_bytes(), at)
    at += buf.byte_length()
    buf = String()


def write_snarls_tsv(path: String, trees: List[SnarlTree], stats: Bool = False) raises:
    """The trees as a table, a line per snarl: `component start end parent_start parent_end chain rank chain_length
    chain_kind tied` (visits as `id << 1 | reverse`, `-` for no parent), each component's in vg snarls' order.
    With `stats`, a line per component on stderr: its sizes, each step's milliseconds and the ties met."""
    var f = OutFile(path)
    var at = 0
    var buf = String()
    for k in range(len(trees)):
        ref t = trees[k]
        for i in range(len(t)):
            buf.write(k, "\t", t.start[i], "\t", t.end[i], "\t")
            var p = t.parent[i]
            if p == NONE:
                buf.write("-\t-")
            else:
                buf.write(t.start[Int(p)], "\t", t.end[Int(p)])
            var c = Int(t.chain[i])
            buf.write(
                "\t",
                c,
                "\t",
                t.rank[i],
                "\t",
                t.chain_length[c],
                "\t",
                t.chain_kind[c],
                "\t",
                Int(t.chain_tied[c]),
                "\n",
            )
            if buf.byte_length() > (1 << 20):
                _flush(f, buf, at)
        if stats and len(t.stats) > 0:
            ref st = t.stats
            print(
                (
                    t"component {k}: snarls {len(t)} (all {st[5]}) chains {st[6]} cactus {st[7]} forest {st[8]}; ms"
                    t" cactus {st[0] // 1000000} cycles {st[1] // 1000000} forest {st[2] // 1000000} decompose"
                    t" {st[3] // 1000000} output {st[4] // 1000000}; ties cycles {st[9]} leaves {st[10]} paths"
                    t" {st[11]}"
                ),
                file=stderr,
            )
    _flush(f, buf, at)
    f.close()


def write_records_tsv(path: String, snarls: List[Snarl]) raises:
    """Snarls records (a vg snarls file's, `dromos.snarls.read_snarls`) as the first five columns of
    `write_snarls_tsv`, all in component 0, for comparing."""
    var f = OutFile(path)
    var at = 0
    var buf = String()
    _put_records(f, buf, at, 0, snarls)
    _flush(f, buf, at)
    f.close()


def write_records_tsv(path: String, components: List[List[Snarl]]) raises:
    """Each component's snarls records (a distance index's, `dromos.distgraph.dist_snarls`) as the first five
    columns of `write_snarls_tsv`, the component its index in `components`."""
    var f = OutFile(path)
    var at = 0
    var buf = String()
    for k in range(len(components)):
        _put_records(f, buf, at, k, components[k])
    _flush(f, buf, at)
    f.close()


def _put_records(mut f: OutFile, mut buf: String, mut at: Int, component: Int, snarls: List[Snarl]) raises:
    for s in snarls:
        buf.write(component, "\t", s.start, "\t", s.end, "\t")
        if s.parent_start == NO_PARENT:
            buf.write("-\t-\n")
        else:
            buf.write(s.parent_start, "\t", s.parent_end, "\n")
        if buf.byte_length() > (1 << 20):
            _flush(f, buf, at)
