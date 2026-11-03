"""One contig's records for the site catalog, or those of regions of it, from its snarls (pgcatalog's catalog.rs,
record.rs and reference.rs, ported): every snarl the reference path passes through, at every nesting level, with all
its alleles, as `vg deconstruct -a` writes them, without what needs haplotype names (genotypes, AC/AF/AN/NS). Their
VCF text (`dromos view`) is pgcatalog's byte for byte.

Pass 1 finds each snarl's reference traversal; snarls the reference does not pass through once have no record. The
rest are sorted by position, ancestors first at a tie, and built in pgcatalog's batches: a snarl's haplotype walks by
forward GBWT search, its alleles (the reference, then distinct sequences, most haplotypes first), then LV/PS from the
ancestors built so far; records are placed in blocks (pgcatalog's: BLOCK_BP windows, BLOCK_BYTES of text) and each
sealed block encoded as the catalog's columns (`encode_block`), in parallel. Each step builds one batch while the
blocks sealed so far are encoded; threads with nothing left in a step build later batches. Heavy snarls are searched
ahead by the pool's workers when they have nothing else to do, their haplotypes in parts; a record whose line would
be longer than STREAM_BYTES is encoded as a job of its own. Every parallel step is a job of the pool (pool.mojo), so
several contigs can be built at once on the same workers.

The crate does the search (a batch of states per call) and gives the node letters, a page of node ids at a time as
they are first needed.
"""

from std.atomic import Atomic
from std.ffi import external_call
from std.memory import Allocation, Layout, OwnedPointer, alloc, dealloc, unsafe_memcpy
from std.os import getenv
from std.sys import stderr
from std.time import perf_counter_ns

from dromos.graph import NO_LENGTH, GbzParts, Metadata, SearchState
from dromos.output import decimal_digits, put_padded
from dromos.parallel import ParallelTask, parallel_for
from dromos.pool import Pool, PoolTask
from dromos.sites import (
    BLOCK_BP,
    BLOCK_BYTES,
    BOUND,
    COLUMNS,
    EXC,
    FLAG_PAD,
    FLAG_PS,
    FLAG_REVERSED,
    FLAGS,
    LV,
    NODEID,
    NODELEN,
    NWALK,
    OWN_BLOCK_BYTES,
    PAD,
    POS,
    PS,
    SEQ,
    WALK,
    WLEN,
    Compressor,
    EncodedBlock,
    LetterPacker,
    complement,
    put_oriented,
    put_varint,
    upper,
    zigzag,
)
from dromos.snarls import Snarl, canonical, put_visit
from dromos.stream import imm_untracked, mut_untracked, untracked

comptime NONE = UInt32.MAX
comptime BATCH = 50_000
"""Snarls built per batch: enough to keep every thread busy, few enough that a batch of records fits in memory."""
comptime STREAM_BYTES = 64 << 20
"""A record whose VCF line would be longer is encoded on a thread of its own."""
comptime SEARCH_BATCH = 64
"""Search states extended per call into the crate."""
comptime PAGE_BITS = 16
comptime PAGE = 1 << PAGE_BITS
"""Node ids whose letters are unpacked together."""


# Node letters, a page at a time.


struct _Page(Copyable, Movable):
    var offsets: List[UInt32]
    """Node `first + i`'s letters are letters[offsets[i] ..< offsets[i + 1]] (none for an id the graph lacks)."""
    var letters: List[UInt8]

    def __init__(out self):
        self.offsets = List[UInt32]()
        self.letters = List[UInt8]()


struct Graph(Movable):
    """The GBZ's parts (the crate's GBWT and node letters), with node letters unpacked a page of PAGE ids at a time,
    by whichever thread first needs one (the others wait for it). Read-only otherwise: any thread may use it."""

    var parts: GbzParts
    var pages: List[_Page]
    var page_ptr: MutPointer[_Page, MutUntrackedOrigin]
    var states: Allocation[Atomic[Int64]]
    """Per page: 0 not unpacked, 1 being unpacked, 2 unpacked, 3 failed."""
    var state_ptr: MutPointer[Atomic[Int64], MutUntrackedOrigin]
    var page_count: Int
    var max_node: Int

    def __init__(out self, path: String, threads: Int) raises:
        self.parts = GbzParts(path)
        self.parts.find_nodes(threads)
        self.max_node = self.parts.max_node()
        self.page_count = (self.max_node + PAGE) >> PAGE_BITS
        self.pages = List[_Page](length=self.page_count, fill=_Page())
        self.page_ptr = mut_untracked(self.pages.unsafe_ptr())
        self.states = alloc(Layout[Atomic[Int64]](length=max(self.page_count, 1)))
        self.state_ptr = mut_untracked(self.states.unsafe_ptr())
        for i in range(self.page_count):
            self.state_ptr.unsafe_offset(i).unsafe_write(Atomic[Int64](0))

    def __deinit__(deinit self):
        dealloc(self.states^)

    def _unpack(self, page: Int) raises:
        var first = page << PAGE_BITS
        var count = min(PAGE, self.max_node + 1 - first)
        var lengths = List[UInt64](unsafe_uninit_length=count)
        self.parts.lengths(first, count, lengths.unsafe_ptr())
        ref p = self.page_ptr[unsafe_offset=page]
        p.offsets = List[UInt32](unsafe_uninit_length=count + 1)
        var total = 0
        for i in range(count):
            p.offsets[i] = UInt32(total)
            if lengths[i] != NO_LENGTH:
                total += Int(lengths[i])
            if total > Int(UInt32.MAX):
                raise Error(t"nodes {first} to {first + count - 1}: more letters than a page holds")
        p.offsets[count] = UInt32(total)
        p.letters = List[UInt8](unsafe_uninit_length=total)
        if self.parts.letters(first, count, p.letters.unsafe_ptr(), total) != total:
            raise Error(t"nodes {first} to {first + count - 1}: letters missing")

    def ensure(self, page: Int) raises:
        """Unpacks the page if no thread has; waits if another is."""
        ref state = self.state_ptr[unsafe_offset=page]
        var s = state.load()
        if s == 2:
            return
        var expected = Int64(0)
        if s == 0 and state.compare_exchange(expected, 1):
            try:
                self._unpack(page)
            except e:
                state.store(3)
                raise e^
            state.store(2)
            return
        while True:
            s = state.load()
            if s == 2:
                return
            if s == 3:
                raise Error("node letters could not be unpacked")
            _ = external_call["sched_yield", Int32]()

    @always_inline
    def letters(self, node: Int) raises -> Span[UInt8, MutUntrackedOrigin]:
        """Node `node`'s letters, forward strand; empty for an id the graph lacks."""
        var page = node >> PAGE_BITS
        if node < 0 or page >= self.page_count:
            return Span[UInt8, MutUntrackedOrigin]()
        self.ensure(page)
        ref p = self.page_ptr[unsafe_offset=page]
        var i = node & (PAGE - 1)
        var a = Int(p.offsets[i])
        return Span(unsafe_ptr=mut_untracked(p.letters.unsafe_ptr().unsafe_offset(a)), length=Int(p.offsets[i + 1]) - a)

    @always_inline
    def length(self, node: Int) raises -> Int:
        return len(self.letters(node))


@fieldwise_init
struct _PrefetchTask(ParallelTask):
    var graph: ImmPointer[Graph, ImmUntrackedOrigin]
    var first: Int
    var bad: MutPointer[Int, MutUntrackedOrigin]

    def run(self, index: Int):
        try:
            self.graph[].ensure(self.first + index)
        except:
            self.bad[] = 1


def prefetch(graph: Graph, lo: Int, hi: Int, threads: Int) raises:
    """Unpacks the pages of node ids `lo ..= hi` on `threads` threads."""
    if hi < lo:
        return
    var first = max(lo, 0) >> PAGE_BITS
    var last = min(hi >> PAGE_BITS, graph.page_count - 1)
    var bad = 0
    var task = _PrefetchTask(untracked(Pointer(to=graph)), first, mut_untracked(Pointer(to=bad)))
    parallel_for(task, last - first + 1, threads)
    if bad != 0:
        raise Error("node letters could not be unpacked")


# The reference path.


struct ReferencePath(Copyable, Movable):
    """A contig's reference path: its steps (visits), their positions, and the steps of each node, from every
    fragment of `<sample>#0#<contig>` laid end to end in offset order; a traversal never runs from one fragment into
    the next."""

    var chrom: String
    """`<sample>#0#<contig>`: the CHROM written in the catalog."""
    var steps: List[UInt64]
    var pos: List[Int]
    """0-based position of each step's first base on the contig."""
    var frag_beg: List[UInt32]
    var frag_end: List[UInt32]
    """For each step, the first step of its fragment and one past its last."""
    var lo: Int
    var first: List[UInt32]
    """For node `lo + i`, its first step (NONE if the path does not visit it)."""
    var next: List[UInt32]
    """For each step, the next step on the same node (NONE after the last)."""
    var length: Int
    """One past the contig's last reference base."""

    def __init__(out self, graph: Graph, md: Metadata, sample: String, contig: String, threads: Int) raises:
        self.chrom = String(t"{sample}#0#{contig}")
        var missing = String(t"no reference path {self.chrom} in the graph")
        var sample_id = md.sample_id(sample)
        var contig_id = md.contig_id(contig)
        if not sample_id or not contig_id:
            raise Error(missing)
        var fragments = List[Tuple[Int, Int]]()
        var names = md.path_names()
        for id in range(len(names)):
            var p = names[id]
            if p.sample == sample_id.value() and p.contig == contig_id.value() and p.phase == 0:
                fragments.append((p.fragment, id))
        if len(fragments) == 0:
            raise Error(missing)

        def by_offset(a: Tuple[Int, Int], b: Tuple[Int, Int]) -> Bool:
            return a[0] < b[0] or (a[0] == b[0] and a[1] < b[1])

        sort(Span(fragments), by_offset)
        self.steps = List[UInt64]()
        self.pos = List[Int]()
        self.frag_beg = List[UInt32]()
        self.frag_end = List[UInt32]()
        self.length = 0
        var walks = List[List[UInt64]]()
        var lo = Int.MAX
        var hi = -1
        for f in fragments:
            var w = graph.parts.path(f[1])
            for v in w:
                lo = min(lo, Int(v >> 1))
                hi = max(hi, Int(v >> 1))
            walks.append(w^)
        prefetch(graph, lo, hi, threads)
        for k in range(len(fragments)):
            var beg = UInt32(len(self.steps))
            var at = fragments[k][0]
            for v in walks[k]:
                self.pos.append(at)
                var n = graph.length(Int(v >> 1))
                if n == 0:
                    raise Error(t"{self.chrom}: a step on node {v >> 1}, which the graph lacks")
                at += n
                self.steps.append(v)
                self.frag_beg.append(beg)
            var end = UInt32(len(self.steps))
            for _ in range(len(self.frag_end), len(self.steps)):
                self.frag_end.append(end)
            self.length = max(self.length, at)
        self.lo = lo if hi >= 0 else 0
        self.first = List[UInt32](length=max(hi - self.lo + 1, 0), fill=NONE)
        self.next = List[UInt32](length=len(self.steps), fill=NONE)
        for k in range(len(self.steps) - 1, -1, -1):
            var i = Int(self.steps[k] >> 1) - self.lo
            self.next[k] = self.first[i]
            self.first[i] = UInt32(k)

    @always_inline
    def first_step(self, node: Int) -> UInt32:
        var i = node - self.lo
        if i < 0 or i >= len(self.first):
            return NONE
        return self.first[i]


@fieldwise_init
struct Traversal(Copyable, TrivialRegisterPassable):
    """A snarl's traversal by the reference path: from step `s0` to `e0`, forward or `backward` (its visits flipped)."""

    var s0: Int
    var e0: Int
    var backward: Bool


@always_inline
def _blocked(block: List[UInt64], node: UInt64) -> Bool:
    for b in block:
        if b == node:
            return True
    return False


def traversals(r: ReferencePath, sn: Snarl, block: List[UInt64], mut found: Traversal) -> Int:
    """vg's PathTraversalFinder on the reference path (pgcatalog's reference.rs `traversals`): how many traversals
    run from the snarl's start to its end, the last into `found`. From each visit of the start node the path is
    walked forward, and if that does not reach the end, backward; it may not step onto a node of `block` (just
    outside the snarl, the start's predecessors and the end's successors), or past the end of its fragment."""
    var k0 = r.first_step(Int(sn.start >> 1))
    if k0 == NONE or r.first_step(Int(sn.end >> 1)) == NONE:
        return 0
    var count = 0
    var end = sn.end
    var end_check = end ^ 1
    var s = k0
    while s != NONE:
        var s0 = Int(s)
        var fbeg = Int(r.frag_beg[s0])
        var fend = Int(r.frag_end[s0])
        var k = s0
        while True:
            var h = r.steps[k]
            if h == end or k + 1 >= fend:
                break
            var nxt = r.steps[k + 1]
            if nxt == end or not _blocked(block, nxt >> 1):
                k += 1
            else:
                break
        if r.steps[k] == end:
            found = Traversal(s0, k, False)
            count += 1
        else:
            k = s0
            while True:
                var h = r.steps[k] ^ 1
                if h == end or k == fbeg:
                    break
                var prev = r.steps[k - 1]
                if prev == end_check or not _blocked(block, prev >> 1):
                    k -= 1
                else:
                    break
            if r.steps[k] == end_check:
                found = Traversal(s0, k, True)
                count += 1
        s = r.next[s0]
    return count


# Records.


struct Record(Movable):
    """A snarl's record: its walks (the reference's, then one per distinct ALT sequence), and how its line is
    written. An allele's bases are its walk's inner nodes' letters, written as they are needed."""

    var present: Bool
    """False when every haplotype spells the reference: no record."""
    var built: Bool
    var pos: Int
    var start: UInt64
    var end: UInt64
    var visits: List[UInt64]
    var starts: List[Int]
    var ends: List[Int]
    """Walk k is visits[starts[k] ..< ends[k]] (the search's walks are kept where it found them)."""
    var lens: List[Int]
    """Each allele's length as written, padding included."""
    var pad: Int
    """The base written before every allele when they are not all one length, else -1."""
    var reversed: Bool
    """The reference reaches the snarl's end node first: alleles are written reverse-complemented and walks
    backwards."""
    var at_length: Int
    """The AT field's walks' bytes, separators excluded."""

    def __init__(out self):
        self.present = False
        self.built = False
        self.pos = 0
        self.start = 0
        self.end = 0
        self.visits = List[UInt64]()
        self.starts = List[Int]()
        self.ends = List[Int]()
        self.lens = List[Int]()
        self.pad = -1
        self.reversed = False
        self.at_length = 0

    @always_inline
    def walks(self) -> Int:
        return len(self.starts)

    def last(self) -> Int:
        """Last reference base the record covers, as tabix computes it for VCF."""
        return self.pos + max(self.lens[0], 1) - 1

    def min_allele(self) -> Int:
        """The shortest allele's bases, the padding base not counted."""
        var m = Int.MAX
        for l in self.lens:
            m = min(m, l)
        return m - Int(self.pad >= 0)

    def max_walk(self) -> Int:
        var m = 0
        for k in range(self.walks()):
            m = max(m, self.ends[k] - self.starts[k])
        return m

    def id_length(self) -> Int:
        return 2 + decimal_digits(Int(self.start >> 1)) + decimal_digits(Int(self.end >> 1))

    def line_length(self, chrom_length: Int, lv: Int, ps_length: Int) -> Int:
        """The line's length, newline excluded; `ps_length` 0 for no PS."""
        var alleles = 0
        for l in self.lens:
            alleles += l
        var n = self.walks()
        var length = chrom_length + 1 + decimal_digits(self.pos) + 1 + self.id_length() + 1
        length += alleles + (n - 2) + StaticString("\t\t60\t.\tAT=").byte_length() + self.at_length + (n - 1)
        length += StaticString(";LV=").byte_length() + decimal_digits(lv)
        if ps_length > 0:
            length += StaticString(";PS=").byte_length() + ps_length
        return length


# Building records.


@always_inline
def _hash_bytes(p: ImmPointer[UInt8, _], n: Int) -> UInt64:
    """A 64-bit hash of `p[:n]`, eight bytes at a time (equal hashes are confirmed byte by byte)."""
    var h = UInt64(0x9E3779B97F4A7C15) ^ UInt64(n)
    var i = 0
    while i + 8 <= n:
        var w = p.unsafe_offset(i).unsafe_bitcast[UInt64]().unsafe_load[alignment=1]()
        h = (h ^ w) * 0xFF51AFD7ED558CCD
        h ^= h >> 32
        i += 8
    var tail = UInt64(0)
    var shift = UInt64(0)
    while i < n:
        tail |= UInt64(p[unsafe_offset=i]) << shift
        shift += 8
        i += 1
    h = (h ^ tail) * 0xC4CEB9FE1A85EC53
    return h ^ (h >> 29)


def _text(graph: Graph, visits: Span[UInt64, _], mut out: List[UInt8]) raises:
    """A walk's bases between its boundary nodes, in the snarl's orientation (record.rs `walk_string`)."""
    out.clear()
    for i in range(1, len(visits) - 1):
        var v = visits[i]
        var s = graph.letters(Int(v >> 1))
        var at = len(out)
        out.resize(unsafe_uninit_length=at + len(s))
        put_oriented(out.unsafe_ptr().unsafe_offset(at), s, v & 1 == 1)


def _walk_length(graph: Graph, visits: Span[UInt64, _]) raises -> Int:
    var n = 0
    for i in range(1, len(visits) - 1):
        n += graph.length(Int(visits[i] >> 1))
    return n


struct WalkSet(Movable, Sized):
    """Haplotype walks, one after another, with how many haplotypes take each."""

    var visits: List[UInt64]
    var bounds: List[Int]
    """Walk k is visits[bounds[k] ..< bounds[k + 1]]."""
    var counts: List[Int]

    def __init__(out self):
        self.visits = List[UInt64]()
        self.bounds = List[Int]()
        self.bounds.append(0)
        self.counts = List[Int]()

    def clear(mut self):
        self.visits.clear()
        self.bounds.clear()
        self.bounds.append(0)
        self.counts.clear()

    @always_inline
    def __len__(self) -> Int:
        return len(self.counts)

    @always_inline
    def walk(self, k: Int) -> Span[UInt64, origin_of(self.visits)]:
        return Span(self.visits)[self.bounds[k] : self.bounds[k + 1]]


@always_inline
def _visits_less(a: Span[UInt64, _], b: Span[UInt64, _]) -> Bool:
    for i in range(min(len(a), len(b))):
        if a[i] != b[i]:
            return a[i] < b[i]
    return len(a) < len(b)


comptime SEGMENT = 32
"""Visits per segment of the search's paths."""


struct Searcher(Movable):
    """One thread's haplotype walk search, its buffers kept from snarl to snarl.

    The search's paths share prefixes: a path is a chain of segments, each up to SEGMENT visits in a row in `pool`
    and a parent (a segment and how many of its visits the path takes). A state's first extension continues its
    segment in place while it has room and no other path has continued it; others start segments of their own. A
    walk found is copied out a segment at a time, so the paths extended side by side (a batch of states per call into
    the crate) still read back in runs."""

    var pool: List[UInt64]
    var seg_len: List[Int]
    var seg_parent: List[Int]
    var seg_parent_len: List[Int]
    var stack_seg: List[Int]
    var stack_len: List[Int]
    var stack_state: List[SearchState]
    var batch_seg: List[Int]
    var batch_len: List[Int]
    var batch_state: List[SearchState]
    var next_states: List[SearchState]
    var counts: List[Int64]
    var chain: List[Int]

    def __init__(out self):
        self.pool = List[UInt64]()
        self.seg_len = List[Int]()
        self.seg_parent = List[Int]()
        self.seg_parent_len = List[Int]()
        self.stack_seg = List[Int]()
        self.stack_len = List[Int]()
        self.stack_state = List[SearchState]()
        self.batch_seg = List[Int]()
        self.batch_len = List[Int]()
        self.batch_state = List[SearchState]()
        self.next_states = List[SearchState]()
        self.counts = List[Int64]()
        self.chain = List[Int]()

    @always_inline
    def _segment(mut self, parent: Int, parent_len: Int, v: UInt64) -> Int:
        var seg = len(self.seg_len)
        self.pool.resize(unsafe_uninit_length=len(self.pool) + SEGMENT)
        self.pool[seg * SEGMENT] = v
        self.seg_len.append(1)
        self.seg_parent.append(parent)
        self.seg_parent_len.append(parent_len)
        return seg

    def _copy_out(mut self, seg: Int, length: Int, last: UInt64, mut out: WalkSet):
        """The path ending at visit `length` of `seg`, then `last`, as a walk."""
        self.chain.clear()
        var s = seg
        var n = length
        while s >= 0:
            self.chain.append(s)
            self.chain.append(n)
            n = self.seg_parent_len[s]
            s = self.seg_parent[s]
        var k = len(self.chain) - 2
        while k >= 0:
            var at = self.chain[k] * SEGMENT
            out.visits.extend(Span(self.pool)[at : at + self.chain[k + 1]])
            k -= 2
        out.visits.append(last)
        out.bounds.append(len(out.visits))

    def search(mut self, graph: Graph, first: SearchState, start: UInt64, end: UInt64, mut out: WalkSet) raises:
        """vg's list_haplotypes (record.rs `haplotype_walks`) from the haplotypes of `first` (a state at `start`):
        the distinct walks they take from `start` to `end`, with how many take each, by forward GBWT search,
        appended to `out`. A walk that steps onto a boundary node backwards has left the snarl and is dropped.
        States are extended SEARCH_BATCH at a time (the walks come in another order than pgcatalog's: the same walks
        and counts)."""
        var start_back = start ^ 1
        var end_back = end ^ 1
        self.pool.clear()
        self.seg_len.clear()
        self.seg_parent.clear()
        self.seg_parent_len.clear()
        _ = self._segment(-1, 0, start)
        self.stack_seg.clear()
        self.stack_len.clear()
        self.stack_state.clear()
        self.stack_seg.append(0)
        self.stack_len.append(1)
        self.stack_state.append(first)
        while len(self.stack_seg) > 0:
            var take = min(SEARCH_BATCH, len(self.stack_seg))
            self.batch_seg.clear()
            self.batch_len.clear()
            self.batch_state.clear()
            for _ in range(take):
                self.batch_seg.append(self.stack_seg.pop())
                self.batch_len.append(self.stack_len.pop())
                self.batch_state.append(self.stack_state.pop())
            graph.parts.follow_forward(self.batch_state.unsafe_ptr(), take, self.next_states, self.counts)
            var j = 0
            for b in range(take):
                var seg = self.batch_seg[b]
                var length = self.batch_len[b]
                for _ in range(max(Int(self.counts[b]), 0)):
                    var nxt = self.next_states[j]
                    j += 1
                    var v = nxt.forward_node
                    if v == end_back or v == start_back:
                        continue
                    if v == end:
                        self._copy_out(seg, length, v, out)
                        out.counts.append(nxt.count())
                        continue
                    if length == self.seg_len[seg] and length < SEGMENT:
                        # The path's first extension: its segment continues in place.
                        self.pool[seg * SEGMENT + length] = v
                        self.seg_len[seg] += 1
                        self.stack_seg.append(seg)
                        self.stack_len.append(length + 1)
                    else:
                        self.stack_seg.append(self._segment(seg, length, v))
                        self.stack_len.append(1)
                    self.stack_state.append(nxt)


struct Builder(Movable):
    """One thread's record building and block encoding, its buffers kept from snarl to snarl and block to block."""

    var searcher: Searcher
    var encoder: Encoder
    var found: WalkSet
    var order: List[Int]
    var text_a: List[UInt8]
    var text_b: List[UInt8]

    def __init__(out self) raises:
        self.searcher = Searcher()
        self.encoder = Encoder()
        self.found = WalkSet()
        self.order = List[Int]()
        self.text_a = List[UInt8]()
        self.text_b = List[UInt8]()

    def build(mut self, graph: Graph, r: ReferencePath, sn: Snarl, t: Traversal, mut out: Record) raises:
        """The record of a snarl the reference passes through (record.rs `build_record`); not present when every
        haplotype spells the reference."""
        self.found.clear()
        var first = graph.parts.search_state(Int(sn.start >> 1), sn.start & 1 == 1)
        if first:
            self.searcher.search(graph, first.value(), sn.start, sn.end, self.found)
        self.finish(graph, r, sn, t, out)

    def merge(mut self, parts: Span[WalkSet, _]):
        """The walks of searches of parts of a snarl's haplotypes as one search's: equal walks one, their counts
        summed, into `found`."""
        self.found.clear()
        var all = WalkSet()
        for i in range(len(parts)):
            ref p = parts[i]
            for k in range(len(p)):
                all.visits.extend(p.walk(k))
                all.bounds.append(len(all.visits))
                all.counts.append(p.counts[k])
        self.order.clear()
        for k in range(len(all)):
            self.order.append(k)
        var visits = imm_untracked(all.visits.unsafe_ptr())
        var bounds = imm_untracked(all.bounds.unsafe_ptr())

        def less(x: Int, y: Int) {visits, bounds} -> Bool:
            var ax = bounds[unsafe_offset=x]
            var ay = bounds[unsafe_offset=y]
            return _visits_less(
                Span(unsafe_ptr=visits.unsafe_offset(ax), length=bounds[unsafe_offset=x + 1] - ax),
                Span(unsafe_ptr=visits.unsafe_offset(ay), length=bounds[unsafe_offset=y + 1] - ay),
            )

        sort(Span(self.order), less)
        var i = 0
        while i < len(self.order):
            var k = self.order[i]
            var count = all.counts[k]
            var j = i + 1
            while j < len(self.order) and all.walk(self.order[j]) == all.walk(k):
                count += all.counts[self.order[j]]
                j += 1
            self.found.visits.extend(all.walk(k))
            self.found.bounds.append(len(self.found.visits))
            self.found.counts.append(count)
            i = j
        _ = len(all)

    def finish(mut self, graph: Graph, r: ReferencePath, sn: Snarl, t: Traversal, mut out: Record) raises:
        """The record from the snarl's haplotype walks, in `found`."""
        out.present = False
        # The reference's walk, in the snarl's orientation.
        var ref_walk = List[UInt64]()
        if t.backward:
            for k in range(t.s0, t.e0 - 1, -1):
                ref_walk.append(r.steps[k] ^ 1)
        else:
            for k in range(t.s0, t.e0 + 1):
                ref_walk.append(r.steps[k])
        var start_pos = r.pos[t.s0]
        var end_pos = r.pos[t.e0]
        # As vg: VCF position and padding base come from the boundary node the reference reaches first.
        var use_start = start_pos < end_pos
        var first_path = r.steps[t.s0 if use_start else t.e0]
        var first_path_pos = start_pos if use_start else end_pos
        var first_trav = ref_walk[0] if use_start else ref_walk[len(ref_walk) - 1]
        var same = ((first_trav & 1) == (first_path & 1)) == use_start
        var fp = graph.letters(Int(first_path >> 1))
        if len(fp) == 0:
            raise Error(t"node {first_path >> 1}: not in the graph")
        var fp_reverse = first_path & 1 == 1
        var prev_char: UInt8
        if same:
            first_path_pos += len(fp)
            # The oriented sequence's last letter.
            prev_char = complement(fp[0]) if fp_reverse else upper(fp[len(fp) - 1])
        else:
            # The complement of its first letter.
            var first_letter = complement(fp[len(fp) - 1]) if fp_reverse else upper(fp[0])
            prev_char = complement(first_letter)
        var pos = first_path_pos + 1

        # Alleles: the reference, then distinct sequences, most haplotypes first.
        ref found = self.found
        self.order.clear()
        var ref_span = Span(ref_walk)
        for w in range(len(found)):
            if found.walk(w) != ref_span:
                self.order.append(w)
        if len(self.order) == 0:
            return
        var found_ptr = imm_untracked(found.visits.unsafe_ptr())
        var bounds_ptr = imm_untracked(found.bounds.unsafe_ptr())
        var counts_ptr = imm_untracked(found.counts.unsafe_ptr())

        def before(x: Int, y: Int) {found_ptr, bounds_ptr, counts_ptr} -> Bool:
            var cx = counts_ptr[unsafe_offset=x]
            var cy = counts_ptr[unsafe_offset=y]
            if cx != cy:
                return cx > cy
            var ax = bounds_ptr[unsafe_offset=x]
            var ay = bounds_ptr[unsafe_offset=y]
            return _visits_less(
                Span(unsafe_ptr=found_ptr.unsafe_offset(ax), length=bounds_ptr[unsafe_offset=x + 1] - ax),
                Span(unsafe_ptr=found_ptr.unsafe_offset(ay), length=bounds_ptr[unsafe_offset=y + 1] - ay),
            )

        sort(Span(self.order), before)
        var ref_len = _walk_length(graph, ref_span)
        # Kept walks: index into `found` (-1: the reference), length, hash once needed. Equal sequences have equal
        # lengths, so only walks as long as a kept one are hashed, and equal hashes confirmed byte by byte.
        var kept = List[Int]()
        var kept_len = List[Int]()
        var kept_hash = List[UInt64]()
        var kept_hashed = List[Bool]()
        kept.append(-1)
        kept_len.append(ref_len)
        kept_hash.append(0)
        kept_hashed.append(False)
        var substitution = True
        for w in self.order:
            var walk = found.walk(w)
            var length = _walk_length(graph, walk)
            var hashed = False
            var h = UInt64(0)
            var seen = False
            for k in range(len(kept)):
                if kept_len[k] != length:
                    continue
                if not hashed:
                    _text(graph, walk, self.text_a)
                    h = _hash_bytes(self.text_a.unsafe_ptr(), len(self.text_a))
                    hashed = True
                var other = ref_span if kept[k] < 0 else found.walk(kept[k])
                if not kept_hashed[k]:
                    _text(graph, other, self.text_b)
                    kept_hash[k] = _hash_bytes(self.text_b.unsafe_ptr(), len(self.text_b))
                    kept_hashed[k] = True
                if kept_hash[k] != h:
                    continue
                _text(graph, other, self.text_b)
                if Span(self.text_a) == Span(self.text_b):
                    seen = True
                    break
            if seen:
                continue
            substitution = substitution and length == ref_len
            kept.append(w)
            kept_len.append(length)
            kept_hash.append(h)
            kept_hashed.append(hashed)
        if len(kept) < 2:
            return
        out.present = True
        out.start = sn.start
        out.end = sn.end
        out.reversed = not use_start
        # The walks stay where the search put them, the reference's after them: the record takes the buffer.
        out.starts.clear()
        out.ends.clear()
        out.lens.clear()
        out.starts.append(len(found.visits))
        out.ends.append(len(found.visits) + len(ref_walk))
        for k in range(1, len(kept)):
            out.starts.append(found.bounds[kept[k]])
            out.ends.append(found.bounds[kept[k] + 1])
        found.visits.extend(Span(ref_walk))
        out.visits.clear()
        swap(out.visits, found.visits)
        out.pad = -1 if substitution else Int(prev_char)
        for l in kept_len:
            out.lens.append(l + Int(not substitution))
        if not substitution:
            pos -= 1
        out.pos = pos
        var at = 0
        for k in range(len(out.starts)):
            for i in range(out.starts[k], out.ends[k]):
                at += 1 + decimal_digits(Int(out.visits[i] >> 1))
        out.at_length = at


# A contig's catalog.


@fieldwise_init
struct RefTrav(Copyable, TrivialRegisterPassable):
    """A snarl's single reference traversal, and its sort key: the 0-based position of its leftmost boundary node,
    then its depth in the snarl tree."""

    var snarl: Int
    var t: Traversal
    var key: Int
    var depth: Int


@fieldwise_init
struct Summary(Copyable, Movable):
    var snarls: Int
    var on_reference: Int
    var records: Int
    var blocks: Int
    var seconds: Float64


comptime TRAVERSE_SPAN = 4096
"""Snarls per pass-1 task."""
comptime NO_TRAVERSAL = -1


@fieldwise_init
struct _TraverseTask(PoolTask):
    """Pass 1 on snarls `index * TRAVERSE_SPAN ..`: each one's reference traversal, if it has exactly one (and, with
    regions, can overlap one), into `found` (NO_TRAVERSAL's s0 otherwise)."""

    var graph: ImmPointer[Graph, ImmUntrackedOrigin]
    var reference: ImmPointer[ReferencePath, ImmUntrackedOrigin]
    var snarls: ImmPointer[Snarl, ImmUntrackedOrigin]
    var count: Int
    var found: MutPointer[Traversal, MutUntrackedOrigin]
    var regions: ImmPointer[Tuple[Int, Int], ImmUntrackedOrigin]
    var region_count: Int
    var bad: MutPointer[Int, MutUntrackedOrigin]

    def run(self, index: Int, worker: Int):
        ref r = self.reference[]
        ref g = self.graph[]
        var lo = index * TRAVERSE_SPAN
        var hi = min(lo + TRAVERSE_SPAN, self.count)
        var cand = List[Int]()
        var starts = List[UInt64]()
        var ends = List[UInt64]()
        for i in range(lo, hi):
            self.found[unsafe_offset=i] = Traversal(NO_TRAVERSAL, 0, False)
            var sn = self.snarls[unsafe_offset=i]
            if r.first_step(Int(sn.start >> 1)) != NONE and r.first_step(Int(sn.end >> 1)) != NONE:
                cand.append(i)
                starts.append(sn.start)
                ends.append(sn.end)
        var preds = List[UInt64]()
        var pred_counts = List[Int64]()
        var succs = List[UInt64]()
        var succ_counts = List[Int64]()
        try:
            g.parts.edges(starts.unsafe_ptr(), len(starts), True, preds, pred_counts)
            g.parts.edges(ends.unsafe_ptr(), len(ends), False, succs, succ_counts)
        except:
            self.bad[] = 1
            return
        var block = List[UInt64]()
        var pa = 0
        var sa = 0
        for j in range(len(cand)):
            var i = cand[j]
            var sn = self.snarls[unsafe_offset=i]
            block.clear()
            var pc = max(Int(pred_counts[j]), 0)
            for k in range(pc):
                block.append(preds[pa + k] >> 1)
            pa += pc
            var sc = max(Int(succ_counts[j]), 0)
            for k in range(sc):
                block.append(succs[sa + k] >> 1)
            sa += sc
            var t = Traversal(NO_TRAVERSAL, 0, False)
            if traversals(r, sn, block, t) != 1:
                continue
            if self.region_count > 0:
                # The record lies within the boundary nodes' span, or starts one base before it; exact overlap is
                # checked on the record itself.
                var left = t.s0 if r.pos[t.s0] <= r.pos[t.e0] else t.e0
                var right = t.e0 if left == t.s0 else t.s0
                var span_end: Int
                try:
                    span_end = r.pos[right] + g.length(Int(r.steps[right] >> 1))
                except:
                    self.bad[] = 1
                    return
                var near = False
                for q in range(self.region_count):
                    var g = self.regions[unsafe_offset=q]
                    near = near or not (r.pos[left] > g[1] or span_end + 1 < g[0])
                if not near:
                    continue
            self.found[unsafe_offset=i] = t


comptime LOOKAHEAD = 8
"""Batches after the one a step builds whose records fill the step's end."""
comptime HEAVY_SPAN = 50_000
comptime HEAVY_NODES = 2_000
"""Snarls spanning at least HEAVY_SPAN reference bases, or HEAVY_NODES node ids (a graph's ids run along it, so the
nodes between a snarl's boundaries), are heavy: their haplotypes are searched ahead of their turn, by helper threads
from the start of the contig, split into parts (each a range of the haplotypes) that threads take one at a time."""

comptime _PART = 0
comptime _HEAVY = 1
comptime _BUILD = 2
comptime _UNIT = 3


@fieldwise_init
struct _Job(Copyable, TrivialRegisterPassable):
    """A step's task: a part of a heavy snarl's haplotypes, the next no thread has taken (`a`: the heavy snarl), the
    heavy snarl's record awaited (`a`), another record built (`a`: the batch's record), or a sealed unit encoded
    (`a`: the unit)."""

    var kind: Int
    var a: Int


@fieldwise_init
struct _Item(Copyable, TrivialRegisterPassable):
    """A record to place: its index in the batch, LV, PS (the parent snarl's visits, if `has_ps`), and its VCF line's
    length (which decides its block)."""

    var record: Int
    var lv: Int
    var ps_start: UInt64
    var ps_end: UInt64
    var has_ps: Bool
    var length: Int


comptime _TAKEN = 0
comptime _LEFT = 1
comptime _DONE = 2


struct _Heavy(Movable):
    """A heavy snarl: its parts (ranges of the start state's haplotypes) and their walks, the record built from them,
    and counters: parts taken, parts left to search (the thread finishing the last builds the record), done."""

    var refs_at: Int
    var record: Record
    var parts: List[WalkSet]
    var ranges: List[Tuple[Int, Int]]
    var state: Allocation[Atomic[Int64]]
    var state_ptr: MutPointer[Atomic[Int64], MutUntrackedOrigin]

    def __init__(out self, refs_at: Int, var ranges: List[Tuple[Int, Int]]):
        self.refs_at = refs_at
        self.record = Record()
        self.parts = List[WalkSet]()
        for _ in range(len(ranges)):
            self.parts.append(WalkSet())
        self.state = alloc(Layout[Atomic[Int64]](length=3))
        self.state_ptr = mut_untracked(self.state.unsafe_ptr())
        self.state_ptr.unsafe_offset(_TAKEN).unsafe_write(Atomic[Int64](0))
        self.state_ptr.unsafe_offset(_LEFT).unsafe_write(Atomic[Int64](Int64(len(ranges))))
        self.state_ptr.unsafe_offset(_DONE).unsafe_write(Atomic[Int64](0))
        self.ranges = ranges^

    def __deinit__(deinit self):
        dealloc(self.state^)


@fieldwise_init
struct _Shared(Movable):
    """What a contig's tasks read: the graph, the reference, the snarls and their traversals, the heavy snarls."""

    var graph: ImmPointer[Graph, ImmUntrackedOrigin]
    var reference: ImmPointer[ReferencePath, ImmUntrackedOrigin]
    var snarls: ImmPointer[Snarl, ImmUntrackedOrigin]
    var refs: ImmPointer[RefTrav, ImmUntrackedOrigin]
    var heavy: MutPointer[_Heavy, MutUntrackedOrigin]
    var bad: MutPointer[Int, MutUntrackedOrigin]


def _search_part(sh: _Shared, mut b: Builder, h: Int, p: Int) raises:
    """Part `p` of heavy snarl `h`'s haplotypes searched; the last part's thread builds the record."""
    ref heavy = sh.heavy[unsafe_offset=h]
    var rt = sh.refs[unsafe_offset=heavy.refs_at]
    var sn = sh.snarls[unsafe_offset=rt.snarl]
    ref g = sh.graph[]
    var first = g.parts.search_state(Int(sn.start >> 1), sn.start & 1 == 1)
    if first:
        # A range of the start state's haplotypes: the forward range is what the search extends, and each haplotype
        # follows its own walk, so the parts' walks and counts add up to the whole's.
        var s = first.value()
        var range = heavy.ranges[p]
        s.forward_end = s.forward_start + UInt64(range[1])
        s.forward_start += UInt64(range[0])
        b.searcher.search(g, s, sn.start, sn.end, heavy.parts[p])
    if heavy.state_ptr[unsafe_offset=_LEFT].fetch_sub(1) == 1:
        b.merge(Span(heavy.parts))
        b.finish(g, sh.reference[], sn, rt.t, heavy.record)
        heavy.parts = List[WalkSet]()
        heavy.state_ptr[unsafe_offset=_DONE].store(1)


def _take_part(sh: _Shared, mut b: Builder, h: Int) raises -> Bool:
    """Searches the next part of heavy snarl `h` no thread has taken; false when none is left."""
    ref heavy = sh.heavy[unsafe_offset=h]
    var p = Int(heavy.state_ptr[unsafe_offset=_TAKEN].fetch_add(1))
    if p >= len(heavy.ranges):
        return False
    _search_part(sh, b, h, p)
    return True


struct _HelperTask(Movable, PoolTask):
    """The heavy snarls' parts, in position order, searched ahead of the steps that need them: a background job of
    the pool, an index a part (of the heavy snarl whose parts it falls among: `first[h] ..< first[h + 1]`), each
    taking the next part of its snarl not yet taken (a step may have taken them first)."""

    var sh: _Shared
    var first: List[Int]
    var builders: MutPointer[Builder, MutUntrackedOrigin]
    """One per pool worker."""

    def __init__(out self, var sh: _Shared, var first: List[Int], builders: MutPointer[Builder, MutUntrackedOrigin]):
        self.sh = sh^
        self.first = first^
        self.builders = builders

    def run(self, index: Int, worker: Int):
        var a = 0
        var b = len(self.first) - 1
        while b - a > 1:
            var m = (a + b) // 2
            if self.first[m] <= index:
                a = m
            else:
                b = m
        try:
            _ = _take_part(self.sh, self.builders[unsafe_offset=worker], a)
        except:
            self.sh.bad[] = 1


# Blocks: records placed in them, then encoded as the catalog's columns.


struct Placed(Movable):
    """A record as a block holds it: the record, its LV, and its PS (the nearest ancestor with a record), if any."""

    var record: Record
    var lv: Int
    var ps_start: UInt64
    var ps_end: UInt64
    var has_ps: Bool

    def __init__(out self, var record: Record, lv: Int, ps_start: UInt64, ps_end: UInt64, has_ps: Bool):
        self.record = record^
        self.lv = lv
        self.ps_start = ps_start
        self.ps_end = ps_end
        self.has_ps = has_ps


@always_inline
def _at_visit(r: Record, k: Int, i: Int) -> UInt64:
    """Visit `i` of walk `k` as the AT field writes it: when the record is reversed, the walk backwards, each visit
    flipped."""
    return r.visits[r.ends[k] - 1 - i] ^ 1 if r.reversed else r.visits[r.starts[k] + i]


struct IdIndex(Movable):
    """Node ids to their index in a block's node table: open addressing (linear probing) over a power-of-two table
    kept at most half full, its slots' keys id + 1 (0: empty)."""

    var keys: List[UInt64]
    var values: List[UInt32]
    var used: Int

    def __init__(out self):
        self.keys = List[UInt64](length=1024, fill=0)
        self.values = List[UInt32](length=1024, fill=0)
        self.used = 0

    def clear(mut self):
        for i in range(len(self.keys)):
            self.keys[i] = 0
        self.used = 0

    @always_inline
    def _slot(self, key: UInt64) -> Int:
        var h = key * 0x9E3779B97F4A7C15
        return Int(h >> 32) & (len(self.keys) - 1)

    def _grow(mut self):
        var old_keys = self.keys^
        var old_values = self.values^
        self.keys = List[UInt64](length=2 * len(old_keys), fill=0)
        self.values = List[UInt32](length=2 * len(old_keys), fill=0)
        self.used = 0
        for i in range(len(old_keys)):
            if old_keys[i] != 0:
                _ = self.insert(old_keys[i] - 1, old_values[i])

    @always_inline
    def insert(mut self, id: UInt64, value: UInt32) -> Bool:
        """Adds `id` with `value` unless it is there; true if it was added."""
        if 2 * (self.used + 1) > len(self.keys):
            self._grow()
        var key = id + 1
        var mask = len(self.keys) - 1
        var at = self._slot(key)
        while True:
            var k = self.keys[at]
            if k == key:
                return False
            if k == 0:
                self.keys[at] = key
                self.values[at] = value
                self.used += 1
                return True
            at = (at + 1) & mask

    @always_inline
    def set(mut self, id: UInt64, value: UInt32):
        """`id`'s value (it must be there)."""
        var key = id + 1
        var mask = len(self.keys) - 1
        var at = self._slot(key)
        while self.keys[at] != key:
            at = (at + 1) & mask
        self.values[at] = value

    @always_inline
    def get(self, id: UInt64) -> Int:
        """`id`'s value (it must be there)."""
        var key = id + 1
        var mask = len(self.keys) - 1
        var at = self._slot(key)
        while self.keys[at] != key:
            at = (at + 1) & mask
        return Int(self.values[at])


struct Encoder(Movable):
    """One thread's block encoding: its column buffers, node-table scratch and zstd context, kept from block to
    block."""

    var cols: List[List[UInt8]]
    var ids: List[UInt64]
    var index: IdIndex
    var compressor: Compressor

    def __init__(out self) raises:
        self.cols = List[List[UInt8]](length=COLUMNS, fill=List[UInt8]())
        self.ids = List[UInt64]()
        self.index = IdIndex()
        self.compressor = Compressor()


def encode_block(
    graph: Graph, items: Span[Placed, _], mut e: Encoder, level: Int, threads: Int, mut out: EncodedBlock
) raises:
    """The records of a block (in file order) as the catalog's columns, each compressed (docs/catalog-format.md):
    walks as the AT field writes them, their inner visits as indexes into the block's node table (its nodes' ids,
    lengths and letters)."""
    for ref c in e.cols:
        c.clear()
    var prev_pos = items[0].record.pos
    var prev_first = 0
    for i in range(len(items)):
        ref it = items[i]
        ref r = it.record
        put_varint(e.cols[POS], UInt64(r.pos - prev_pos))
        prev_pos = r.pos
        var n0 = r.ends[0] - r.starts[0]
        var first = _at_visit(r, 0, 0)
        var last = _at_visit(r, 0, n0 - 1)
        var ok = (last ^ 1 == r.start and first ^ 1 == r.end) if r.reversed else (
            first == r.start and last == r.end
        )
        if not ok:
            raise Error(t"the record of snarl {r.start}..{r.end}: its first walk does not run between its ID's nodes")
        var flags = UInt8(0)
        if r.reversed:
            flags |= FLAG_REVERSED
        if r.pad >= 0:
            flags |= FLAG_PAD
            e.cols[PAD].append(UInt8(r.pad))
        if it.has_ps:
            flags |= FLAG_PS
        e.cols[FLAGS].append(flags)
        put_varint(e.cols[BOUND], zigzag(Int(first) - prev_first))
        put_varint(e.cols[BOUND], zigzag(Int(last) - Int(first)))
        prev_first = Int(first)
        put_varint(e.cols[NWALK], UInt64(r.walks()))
        for k in range(r.walks()):
            var m = r.ends[k] - r.starts[k]
            if m < 2 or _at_visit(r, k, 0) != first or _at_visit(r, k, m - 1) != last:
                raise Error(t"the record of snarl {r.start}..{r.end}: a walk does not share its first walk's ends")
            put_varint(e.cols[WLEN], UInt64(m - 2))
        put_varint(e.cols[LV], UInt64(it.lv))
        if it.has_ps:
            put_varint(e.cols[PS], zigzag(Int(it.ps_start) - Int(r.start)))
            put_varint(e.cols[PS], zigzag(Int(it.ps_end) - Int(it.ps_start)))
    # The node table: the distinct inner nodes, by id.
    e.ids.clear()
    e.index.clear()
    for i in range(len(items)):
        ref r = items[i].record
        for k in range(r.walks()):
            for j in range(r.starts[k] + 1, r.ends[k] - 1):
                var id = r.visits[j] >> 1
                if e.index.insert(id, 0):
                    e.ids.append(id)
    sort(e.ids)
    for i in range(len(e.ids)):
        e.index.set(e.ids[i], UInt32(i))
    # WALK: each inner visit as its table index (a zigzag difference from the visit before, across walks and
    # records) and its orientation.
    var prev = 0
    for i in range(len(items)):
        ref r = items[i].record
        for k in range(r.walks()):
            var m = r.ends[k] - r.starts[k]
            for x in range(1, m - 1):
                var v = _at_visit(r, k, x)
                var idx = e.index.get(v >> 1)
                put_varint(e.cols[WALK], zigzag(idx - prev) << 1 | (v & 1))
                prev = idx
    # The table: ids as differences, lengths, letters (LetterPacker).
    var prev_id = UInt64(0)
    var packer = LetterPacker()
    var seq = List[UInt8]()
    var exc = List[UInt8]()
    swap(seq, e.cols[SEQ])
    swap(exc, e.cols[EXC])
    for id in e.ids:
        put_varint(e.cols[NODEID], id - prev_id)
        prev_id = id
        var s = graph.letters(Int(id))
        put_varint(e.cols[NODELEN], UInt64(len(s)))
        packer.add(s, seq, exc)
    packer.finish(seq)
    swap(seq, e.cols[SEQ])
    swap(exc, e.cols[EXC])
    out.data.clear()
    for c in range(COLUMNS):
        out.sizes[2 * c] = e.compressor.column(Span(e.cols[c]), level, threads, out.data)
        out.sizes[2 * c + 1] = len(e.cols[c])


struct Unit(Movable):
    """A block as it is filled: its records, then (once sealed) encoded. A giant's record is encoded elsewhere
    (_GiantTask), its block set from there."""

    var window: Int
    var items: List[Placed]
    var text: Int
    """The records' VCF text: their lines and newlines."""
    var first: Int
    var last: Int
    var max_end: Int
    var min_allele: Int
    var max_walk: Int
    var count: Int
    var sealed: Bool
    """Takes no further records: an own block (a record over OWN_BLOCK_BYTES), a giant, or followed by another."""
    var giant: Bool
    var done: Bool
    """Its block is encoded."""
    var block: EncodedBlock

    def __init__(out self, window: Int, pos: Int, sealed: Bool):
        self.window = window
        self.items = List[Placed]()
        self.text = 0
        self.first = pos
        self.last = pos
        self.max_end = 0
        self.min_allele = Int.MAX
        self.max_walk = 0
        self.count = 0
        self.sealed = sealed
        self.giant = False
        self.done = False
        self.block = EncodedBlock()

    def add(mut self, r: Record, length: Int):
        self.text += length + 1
        self.last = r.pos
        self.max_end = max(self.max_end, r.last())
        self.min_allele = min(self.min_allele, r.min_allele())
        self.max_walk = max(self.max_walk, r.max_walk())
        self.count += 1

    def encode(mut self, graph: Graph, mut e: Encoder, level: Int, threads: Int) raises:
        encode_block(graph, Span(self.items), e, level, threads, self.block)
        self.finish()

    def finish(mut self):
        """Its block's index fields; its records freed."""
        self.block.first_pos = self.first
        self.block.last_pos = self.last
        self.block.max_end = self.max_end
        self.block.records = self.count
        self.block.min_allele = self.min_allele
        self.block.max_walk = self.max_walk
        self.items = List[Placed]()
        self.done = True


struct BlockWriter(Movable):
    """A contig's blocks as its records arrive in POS order: `place` puts each in a block (a unit), or `place_giant`
    in one of its own; the caller encodes the sealed units (any thread, each its own) and `take_done` keeps the
    units done, in order. Units are numbered from 0 in file order; `units` holds those not yet taken, from
    `base`."""

    var units: List[Unit]
    var base: Int
    var blocks: List[EncodedBlock]

    def __init__(out self):
        self.units = List[Unit]()
        self.base = 0
        self.blocks = List[EncodedBlock]()

    def _seal_last(mut self):
        if len(self.units) > 0:
            self.units[len(self.units) - 1].sealed = True

    def place(mut self, var item: Placed, length: Int) -> Int:
        """A record whose line is `length` bytes (newline excluded) into its block; returns the unit. A record
        starts a new block when it has one of its own, or the block filling is sealed, of another BLOCK_BP window,
        or BLOCK_BYTES long."""
        var pos = item.record.pos
        var window = pos // BLOCK_BP
        var own = length > OWN_BLOCK_BYTES
        var n = len(self.units)
        var new_block = own or n == 0
        if not new_block:
            ref b = self.units[n - 1]
            new_block = b.sealed or b.window != window or b.text >= BLOCK_BYTES
        if new_block:
            self._seal_last()
            self.units.append(Unit(window, pos, own))
        ref b = self.units[len(self.units) - 1]
        b.add(item.record, length)
        b.items.append(item^)
        return self.base + len(self.units) - 1

    def place_giant(mut self, r: Record, length: Int) -> Int:
        """A unit for a record encoded elsewhere, after every unit so far (which it seals)."""
        self._seal_last()
        var u = Unit(r.pos // BLOCK_BP, r.pos, True)
        u.giant = True
        u.add(r, length)
        self.units.append(u^)
        return self.base + len(self.units) - 1

    def unit(mut self, id: Int) -> MutPointer[Unit, MutUntrackedOrigin]:
        """Unit `id`: valid until the next `place`, `place_giant` or `take_done`."""
        return mut_untracked(self.units.unsafe_ptr().unsafe_offset(id - self.base))

    def seal(mut self):
        """No further records: the block filling is sealed too."""
        self._seal_last()

    def take_done(mut self):
        """Keeps the units done, up to the first not done, in order."""
        var n = 0
        while n < len(self.units) and self.units[n].done:
            n += 1
        if n == 0:
            return
        var rest = List[Unit]()
        for i in range(len(self.units)):
            if i < n:
                var b = EncodedBlock()
                swap(b, self.units[i].block)
                self.blocks.append(b^)
        while len(self.units) > n:
            rest.append(self.units.pop(n))
        self.units = rest^
        self.base += n


struct _GiantOut(Movable):
    """A giant record's block, made on a thread of its own: its block, and 1 once made (2 if that failed)."""

    var block: EncodedBlock
    var state: Atomic[Int64]
    var unit: Int
    var error: String

    def __init__(out self, unit: Int):
        self.block = EncodedBlock()
        self.state = Atomic[Int64](0)
        self.unit = unit
        self.error = String()


struct _GiantTask(Movable, PoolTask):
    """A record whose line is longer than STREAM_BYTES, encoded on a thread of its own (zstd's own threads compress
    its long columns)."""

    var graph: ImmPointer[Graph, ImmUntrackedOrigin]
    var item: List[Placed]
    var level: Int
    var threads: Int
    var result: MutPointer[_GiantOut, MutUntrackedOrigin]

    def __init__(
        out self,
        graph: ImmPointer[Graph, ImmUntrackedOrigin],
        var item: Placed,
        level: Int,
        threads: Int,
        result: MutPointer[_GiantOut, MutUntrackedOrigin],
    ):
        self.graph = graph
        self.item = List[Placed]()
        self.item.append(item^)
        self.level = level
        self.threads = threads
        self.result = result

    def run(self, index: Int, worker: Int):
        try:
            var e = Encoder()
            encode_block(self.graph[], Span(self.item), e, self.level, self.threads, self.result[].block)
            self.result[].state.store(1)
        except e:
            self.result[].error = String(e)
            self.result[].state.store(2)


struct _Batch(Movable):
    """A batch of snarls: their records, then, once built, the items to place, in file order."""

    var lo: Int
    var hi: Int
    var records: List[Record]
    var items: List[_Item]

    def __init__(out self, lo: Int, hi: Int):
        self.lo = lo
        self.hi = hi
        self.records = List[Record]()
        for _ in range(hi - lo):
            self.records.append(Record())
        self.items = List[_Item]()


struct _StepTask(PoolTask):
    """One step's jobs, taken in turn by `workers` workers (each with its own builder): the batch's records built
    (`cur`), and the sealed blocks encoded."""

    var sh: _Shared
    var jobs: ImmPointer[_Job, ImmUntrackedOrigin]
    var count: Int
    var next: MutPointer[Atomic[Int64], MutUntrackedOrigin]
    var cur: MutPointer[_Batch, MutUntrackedOrigin]
    var writer: MutPointer[BlockWriter, MutUntrackedOrigin]
    var builders: MutPointer[Builder, MutUntrackedOrigin]
    var level: Int
    var threads: Int
    var fill_batch: ImmPointer[Int, ImmUntrackedOrigin]
    var fill: ImmPointer[Int, ImmUntrackedOrigin]
    """Records of the next batches to build while a step's last jobs finish: each one's batch (its address) and
    index there."""
    var fill_count: Int
    var fill_next: MutPointer[Atomic[Int64], MutUntrackedOrigin]
    var left: MutPointer[Atomic[Int64], MutUntrackedOrigin]
    """The step's jobs not finished."""

    def __init__(
        out self,
        var sh: _Shared,
        jobs: ImmPointer[_Job, ImmUntrackedOrigin],
        count: Int,
        next: MutPointer[Atomic[Int64], MutUntrackedOrigin],
        cur: MutPointer[_Batch, MutUntrackedOrigin],
        writer: MutPointer[BlockWriter, MutUntrackedOrigin],
        builders: MutPointer[Builder, MutUntrackedOrigin],
        level: Int,
        threads: Int,
        fill_batch: ImmPointer[Int, ImmUntrackedOrigin],
        fill: ImmPointer[Int, ImmUntrackedOrigin],
        fill_count: Int,
        fill_next: MutPointer[Atomic[Int64], MutUntrackedOrigin],
        left: MutPointer[Atomic[Int64], MutUntrackedOrigin],
    ):
        self.fill_batch = fill_batch
        self.fill = fill
        self.fill_count = fill_count
        self.fill_next = fill_next
        self.left = left
        self.sh = sh^
        self.jobs = jobs
        self.count = count
        self.next = next
        self.cur = cur
        self.writer = writer
        self.builders = builders
        self.level = level
        self.threads = threads

    def run(self, seat: Int, worker: Int):
        ref b = self.builders[unsafe_offset=seat]
        while True:
            var i = Int(self.next[].fetch_add(1))
            if i >= self.count:
                break
            try:
                self._job(b, self.jobs[unsafe_offset=i])
            except:
                self.sh.bad[] = 1
            _ = self.left[].fetch_sub(1)
        # Every job is taken: while the last ones run, records of the next batch (the cheapest first).
        while self.left[].load() > 0:
            var f = Int(self.fill_next[].fetch_add(1))
            if f >= self.fill_count:
                break
            ref batch = MutPointer[_Batch, MutUntrackedOrigin](unsafe_from_address=self.fill_batch[unsafe_offset=f])[]
            var k = self.fill[unsafe_offset=f]
            var rt = self.sh.refs[unsafe_offset=batch.lo + k]
            try:
                b.build(
                    self.sh.graph[], self.sh.reference[], self.sh.snarls[unsafe_offset=rt.snarl], rt.t, batch.records[k]
                )
                batch.records[k].built = True
            except:
                self.sh.bad[] = 1

    def _job(self, mut b: Builder, job: _Job) raises:
        ref g = self.sh.graph[]
        if job.kind == _PART:
            _ = _take_part(self.sh, b, job.a)
        elif job.kind == _HEAVY:
            ref heavy = self.sh.heavy[unsafe_offset=job.a]
            var spins = 0
            while heavy.state_ptr[unsafe_offset=_DONE].load() == 0:
                if self.sh.bad[] != 0:
                    raise Error("a heavy snarl's search failed")
                spins += 1
                if spins < 64:
                    _ = external_call["sched_yield", Int32]()
                else:
                    _ = external_call["usleep", Int32](UInt32(200))
        elif job.kind == _BUILD:
            ref batch = self.cur[]
            var rt = self.sh.refs[unsafe_offset=batch.lo + job.a]
            b.build(g, self.sh.reference[], self.sh.snarls[unsafe_offset=rt.snarl], rt.t, batch.records[job.a])
        else:
            var u = self.writer[].unit(job.a)
            u[].encode(g, b.encoder, self.level, self.threads)


struct _Phases:
    """A contig's phase times, summed, on stderr with DROMOS_TIMING=1."""

    var on: Bool
    var t: Int
    var seconds: List[Float64]

    def __init__(out self):
        self.on = Bool(getenv("DROMOS_TIMING"))
        self.t = perf_counter_ns()
        self.seconds = List[Float64](length=9, fill=0.0)

    def lap(mut self, phase: Int):
        var now = perf_counter_ns()
        self.seconds[phase] += Float64(now - self.t) * 1e-9
        self.t = now

    def report(self, chrom: String):
        if not self.on:
            return
        var names: List[String] = [
            "snarls read",
            "parents",
            "pass 1",
            "sorted, heavy snarls",
            "steps planned and placed",
            "steps",
            "blocks kept",
            "finished",
            "prepared, waiting",
        ]
        var text = String(t"dromos timing: {chrom}:")
        for i in range(len(names)):
            text += String(t" {names[i]} {self.seconds[i]} s;")
        print(text, file=stderr)


def _bytes_less(a: Span[UInt8, _], b: Span[UInt8, _]) -> Bool:
    for i in range(min(len(a), len(b))):
        if a[i] != b[i]:
            return a[i] < b[i]
    return len(a) < len(b)


def _id_key(r: Record, mut out: List[UInt8]):
    """The record's ID and a tab: lines at one position are ordered by it (their bytes' order)."""
    out.clear()
    put_visit(out, r.start)
    put_visit(out, r.end)
    out.append(UInt8(ord("\t")))


def _visit_length(v: UInt64) -> Int:
    return 1 + decimal_digits(Int(v >> 1))


def _place(
    mut batch: _Batch,
    mut emitted: List[Bool],
    parent: List[UInt32],
    snarls: List[Snarl],
    refs: List[RefTrav],
    regions: List[Tuple[Int, Int]],
    chrom: String,
) -> Int:
    """A built batch's records to write, in file order: LV / PS from the ancestors built so far, as vg's
    update_nesting_info_tags (the ancestors that have a record count, the nearest names the parent). A record
    outside every region (if any are given) is left out. Returns how many."""
    for k in range(batch.lo, batch.hi):
        emitted[refs[k].snarl] = batch.records[k - batch.lo].present
    batch.items.clear()
    for k in range(batch.lo, batch.hi):
        ref r = batch.records[k - batch.lo]
        if not r.present:
            continue
        if len(regions) > 0:
            var inside = False
            for g in regions:
                inside = inside or (r.pos <= g[1] and r.last() >= g[0])
            if not inside:
                continue
        var lv = 0
        var ps_start = UInt64(0)
        var ps_end = UInt64(0)
        var has_ps = False
        var c = parent[refs[k].snarl]
        while c != NONE:
            var pi = Int(c)
            if emitted[pi]:
                lv += 1
                if not has_ps:
                    ps_start = snarls[pi].start
                    ps_end = snarls[pi].end
                    has_ps = True
            c = parent[pi]
        var ps_length = _visit_length(ps_start) + _visit_length(ps_end) if has_ps else 0
        var length = r.line_length(chrom.byte_length(), lv, ps_length)
        batch.items.append(_Item(k - batch.lo, lv, ps_start, ps_end, has_ps, length))
    var records_ptr = imm_untracked(batch.records.unsafe_ptr())
    var id_a = List[UInt8]()
    var id_b = List[UInt8]()

    def item_less(a: _Item, b: _Item) {records_ptr, mut id_a, mut id_b} -> Bool:
        ref ra = records_ptr[unsafe_offset=a.record]
        ref rb = records_ptr[unsafe_offset=b.record]
        if ra.pos != rb.pos:
            return ra.pos < rb.pos
        # Order within a position is the lines' byte order, which the ID (unique, then a tab) decides.
        _id_key(ra, id_a)
        _id_key(rb, id_b)
        return _bytes_less(Span(id_a), Span(id_b))

    sort(Span(batch.items), item_less)
    return len(batch.items)


@always_inline
def _slot(a: UInt64, b: UInt64, size: Int) -> Int:
    var h = (a * 0x9E3779B97F4A7C15) ^ (b * 0xC2B2AE3D27D4EB4F)
    h ^= h >> 29
    return Int(h) & (size - 1)


def _parts(span: Int, nodes: Int, haplotypes: Int, threads: Int) -> List[Tuple[Int, Int]]:
    """A heavy snarl's parts: ranges of its start state's haplotypes, more for a wider snarl."""
    var n = max(min(max(span // HEAVY_SPAN, nodes // HEAVY_NODES), 16 * threads, haplotypes), 1)
    var out = List[Tuple[Int, Int]]()
    for p in range(n):
        out.append((haplotypes * p // n, haplotypes * (p + 1) // n))
    return out^


struct ContigPlan(Movable):
    """A selection's catalog, prepared: its snarls (and parents), their reference traversals in position order, the
    batches, and the heavy snarls, which helper threads search from here on. Everything the helpers reach is on the
    heap, so the plan may move."""

    var reference: OwnedPointer[ReferencePath]
    var snarls: List[Snarl]
    var parent: List[UInt32]
    var refs: List[RefTrav]
    var bounds: List[Int]
    """Batch b is refs[bounds[b] ..< bounds[b + 1]]."""
    var heavy: List[_Heavy]
    var heavy_of: List[Int]
    var widths: List[Int]
    """Per reference traversal, the node ids between its snarl's boundaries."""
    var regions: List[Tuple[Int, Int]]
    var pool: MutPointer[Pool, MutUntrackedOrigin]
    var rank: Int
    """The contig's place among those built at once: the pool runs its jobs before those of contigs ranked after."""
    var helper_builders: List[Builder]
    var helpers: OwnedPointer[_HelperTask]
    var helper_job: Int
    """The helpers' job in the pool (-1 once stopped)."""
    var bad: List[Int]
    """1 once a task failed."""
    var timing: _Phases
    var t0: Int

    def __init__(
        out self,
        graph: Graph,
        var reference: ReferencePath,
        var snarls: List[Snarl],
        regions: List[Tuple[Int, Int]],
        threads: Int,
        pool: MutPointer[Pool, MutUntrackedOrigin],
        rank: Int,
    ) raises:
        self.t0 = perf_counter_ns()
        self.timing = _Phases()
        self.reference = OwnedPointer(reference^)
        self.regions = regions.copy()
        self.pool = pool
        self.rank = rank
        self.helper_job = -1
        self.bad = List[Int](length=1, fill=0)
        self.heavy = List[_Heavy]()
        self.snarls = snarls^
        self.timing.lap(0)
        var rp = untracked(self.reference.ptr())
        ref r = rp[]
        var n = len(self.snarls)

        # Each snarl's parent, by its key in either direction (the last snarl with a key wins, as a map built in
        # order): an open-addressing table of the keys.
        var size = 1
        while size < 2 * n + 2:
            size <<= 1
        var slot_a = List[UInt64](length=size, fill=UInt64.MAX)
        var slot_b = List[UInt64](length=size, fill=UInt64.MAX)
        var slot_i = List[UInt32](length=size, fill=NONE)
        for i in range(n):
            var c = canonical(self.snarls[i].start, self.snarls[i].end)
            var at = _slot(c[0], c[1], size)
            while slot_i[at] != NONE and (slot_a[at] != c[0] or slot_b[at] != c[1]):
                at = (at + 1) & (size - 1)
            slot_a[at] = c[0]
            slot_b[at] = c[1]
            slot_i[at] = UInt32(i)
        self.parent = List[UInt32](length=n, fill=NONE)
        for i in range(n):
            var sn = self.snarls[i]
            if not sn.has_parent():
                continue
            var c = canonical(sn.parent_start, sn.parent_end)
            var at = _slot(c[0], c[1], size)
            while slot_i[at] != NONE:
                if slot_a[at] == c[0] and slot_b[at] == c[1]:
                    self.parent[i] = slot_i[at]
                    break
                at = (at + 1) & (size - 1)
        self.timing.lap(1)

        # Pass 1: reference traversals, kept when there is exactly one (a snarl the reference crosses twice would
        # need one record per crossing) and, with regions, when it can overlap one.
        var found = List[Traversal](unsafe_uninit_length=n)
        var pass1 = _TraverseTask(
            untracked(Pointer(to=graph)),
            rp,
            imm_untracked(self.snarls.unsafe_ptr()),
            n,
            mut_untracked(found.unsafe_ptr()),
            imm_untracked(self.regions.unsafe_ptr()),
            len(self.regions),
            mut_untracked(self.bad.unsafe_ptr()),
        )
        pool[].run(pass1, (n + TRAVERSE_SPAN - 1) // TRAVERSE_SPAN, 0, rank)
        if self.bad[0] != 0:
            raise Error(t"{r.chrom}: the graph could not be searched")
        self.refs = List[RefTrav]()
        for i in range(n):
            var t = found[i]
            if t.s0 == NO_TRAVERSAL:
                continue
            var depth = 0
            var c = self.parent[i]
            while c != NONE:
                c = self.parent[Int(c)]
                depth += 1
            self.refs.append(RefTrav(i, t, min(r.pos[t.s0], r.pos[t.e0]), depth))
        self.timing.lap(2)

        def ref_less(a: RefTrav, b: RefTrav) -> Bool:
            if a.key != b.key:
                return a.key < b.key
            if a.depth != b.depth:
                return a.depth < b.depth
            return a.snarl < b.snarl

        sort(Span(self.refs), ref_less)

        # Batches: a batch ends where the next key is at least two bases on, so records (POS is the key or one base
        # later) never interleave across batches.
        self.bounds = List[Int]()
        self.bounds.append(0)
        while self.bounds[len(self.bounds) - 1] < len(self.refs):
            var lo = self.bounds[len(self.bounds) - 1]
            var hi = min(lo + BATCH, len(self.refs))
            while hi < len(self.refs) and self.refs[hi].key < self.refs[hi - 1].key + 2:
                hi += 1
            self.bounds.append(hi)

        # Heavy snarls, in position order, and the helpers that search them ahead.
        self.heavy_of = List[Int](length=len(self.refs), fill=-1)
        self.widths = List[Int](capacity=len(self.refs))
        for k in range(len(self.refs)):
            var t = self.refs[k].t
            var sn = self.snarls[self.refs[k].snarl]
            var span = abs(r.pos[t.s0] - r.pos[t.e0])
            var nodes = abs(Int(sn.end >> 1) - Int(sn.start >> 1))
            self.widths.append(nodes)
            if span >= HEAVY_SPAN or nodes >= HEAVY_NODES:
                var first = graph.parts.search_state(Int(sn.start >> 1), sn.start & 1 == 1)
                var haplotypes = first.value().count() if first else 0
                self.heavy_of[k] = len(self.heavy)
                self.heavy.append(_Heavy(k, _parts(span, nodes, haplotypes, threads)))
        var first = List[Int]()
        first.append(0)
        for h in range(len(self.heavy)):
            first.append(first[len(first) - 1] + len(self.heavy[h].ranges))
        var parts = first[len(first) - 1]
        self.helper_builders = List[Builder]()
        for _ in range(pool[].threads if parts > 0 else 0):
            self.helper_builders.append(Builder())
        self.helpers = OwnedPointer(
            _HelperTask(
                _Shared(
                    untracked(Pointer(to=graph)),
                    rp,
                    untracked(self.snarls.unsafe_ptr()),
                    untracked(self.refs.unsafe_ptr()),
                    mut_untracked(self.heavy.unsafe_ptr()),
                    mut_untracked(self.bad.unsafe_ptr()),
                ),
                first^,
                mut_untracked(self.helper_builders.unsafe_ptr()),
            )
        )
        if parts > 0:
            self.helper_job = pool[].submit(self.helpers[], parts, 1, rank)
        self.timing.lap(3)

    def shared(mut self, graph: Graph) -> _Shared:
        return _Shared(
            untracked(Pointer(to=graph)),
            untracked(self.reference.ptr()),
            untracked(self.snarls.unsafe_ptr()),
            untracked(self.refs.unsafe_ptr()),
            mut_untracked(self.heavy.unsafe_ptr()),
            mut_untracked(self.bad.unsafe_ptr()),
        )

    def stop(mut self):
        """Stops the helpers (after the parts they are searching) and waits for them."""
        if self.helper_job >= 0:
            self.pool[].cancel(self.helper_job)
            self.helper_job = -1

    def __deinit__(deinit self):
        if self.helper_job >= 0:
            self.pool[].cancel(self.helper_job)


struct ContigTail(Movable):
    """A catalog whose records are all placed: the giant records still encoding on threads of their own, and the
    blocks left to keep."""

    var plan: ContigPlan
    var writer: BlockWriter
    var giants: List[OwnedPointer[_GiantOut]]
    var giant_tasks: List[OwnedPointer[_GiantTask]]
    var giant_jobs: List[Int]
    """The giants' jobs in the pool."""
    var records: Int

    def __init__(
        out self,
        var plan: ContigPlan,
        var writer: BlockWriter,
        var giants: List[OwnedPointer[_GiantOut]],
        var giant_tasks: List[OwnedPointer[_GiantTask]],
        var giant_jobs: List[Int],
        records: Int,
    ):
        self.plan = plan^
        self.writer = writer^
        self.giants = giants^
        self.giant_tasks = giant_tasks^
        self.giant_jobs = giant_jobs^
        self.records = records


def _giants_done(mut giants: List[OwnedPointer[_GiantOut]], mut w: BlockWriter):
    """The blocks of the giants made since, to their units."""
    for i in range(len(giants)):
        if giants[i][].state.load() == 1 and giants[i][].unit >= w.base:
            var u = w.unit(giants[i][].unit)
            if not u[].done:
                swap(u[].block, giants[i][].block)
                u[].finish()


def run_steps(graph: Graph, var plan: ContigPlan, level: Int, threads: Int) raises -> ContigTail:
    """A prepared catalog's records built, placed and encoded, in steps: a step's tasks build one batch's records
    and, at the same time, encode the blocks sealed so far; between steps the batch just built is placed in blocks
    (one thread) and the blocks done are kept, in order. Once every task of a step is taken, the threads build
    records of the next LOOKAHEAD batches until its last tasks finish, so no step waits idle on one long task. Heavy
    snarls' parts no helper has taken are searched in the step that needs them; giant records are encoded on threads
    of their own."""
    plan.timing.lap(8)
    var sh = plan.shared(graph)
    var batches = len(plan.bounds) - 1
    ref r = plan.reference[]
    var chrom = r.chrom
    var w = BlockWriter()
    var emitted = List[Bool](length=len(plan.snarls), fill=False)
    var builders = List[Builder]()
    for _ in range(threads):
        builders.append(Builder())
    var giants = List[OwnedPointer[_GiantOut]]()
    var giant_tasks = List[OwnedPointer[_GiantTask]]()
    var giant_jobs = List[Int]()
    var records_written = 0
    # The batch built and up to LOOKAHEAD after it, whose records fill the steps' ends.
    var cur = OwnedPointer(_Batch(plan.bounds[0], plan.bounds[1]) if batches > 0 else _Batch(0, 0))
    var futures = List[OwnedPointer[_Batch]]()
    var made = min(batches, 1)
    for step in range(batches + 1):
        while len(futures) < LOOKAHEAD and made < batches:
            futures.append(OwnedPointer(_Batch(plan.bounds[made], plan.bounds[made + 1])))
            made += 1
        if step == batches:
            w.seal()
        # The jobs: heavy snarls' parts no helper has taken first, then the sealed blocks, the other records, and
        # last the heavy snarls' records awaited (a helper may still be searching a part).
        var jobs = List[_Job]()
        for k in range(cur[].lo, cur[].hi):
            var h = plan.heavy_of[k]
            if h >= 0:
                var taken = Int(plan.heavy[h].state_ptr[unsafe_offset=_TAKEN].load())
                for _ in range(len(plan.heavy[h].ranges) - taken):
                    jobs.append(_Job(_PART, h))
        for u in range(len(w.units)):
            ref unit = w.units[u]
            if unit.sealed and not unit.done and not unit.giant:
                jobs.append(_Job(_UNIT, w.base + u))
        # Records widest first (the node ids between a snarl's boundaries, as the heavy rule counts them): the
        # costly ones start early and the step ends on cheap ones.
        var first_build = len(jobs)
        for k in range(cur[].lo, cur[].hi):
            if plan.heavy_of[k] < 0 and not cur[].records[k - cur[].lo].built:
                jobs.append(_Job(_BUILD, k - cur[].lo))
        var widths = imm_untracked(plan.widths.unsafe_ptr().unsafe_offset(cur[].lo))

        def wider_job(a: _Job, b: _Job) {widths} -> Bool:
            var wa = widths[unsafe_offset=a.a]
            var wb = widths[unsafe_offset=b.a]
            return wa > wb or (wa == wb and a.a < b.a)

        sort(Span(jobs)[first_build:], wider_job)
        for k in range(cur[].lo, cur[].hi):
            if plan.heavy_of[k] >= 0:
                jobs.append(_Job(_HEAVY, plan.heavy_of[k]))
        # The next batches' other records not built, batch by batch, each's cheapest first, to fill the step's end.
        var fill = List[Int]()
        var fill_batch = List[Int]()
        for f in range(len(futures)):
            var address = Int(futures[f].ptr())
            ref fb = futures[f][]
            var first = len(fill)
            for k in range(fb.lo, fb.hi):
                if plan.heavy_of[k] < 0 and not fb.records[k - fb.lo].built:
                    fill.append(k - fb.lo)
                    fill_batch.append(address)
            var fb_widths = imm_untracked(plan.widths.unsafe_ptr().unsafe_offset(fb.lo))

            def narrower(a: Int, b: Int) {fb_widths} -> Bool:
                var wa = fb_widths[unsafe_offset=a]
                var wb = fb_widths[unsafe_offset=b]
                return wa < wb or (wa == wb and a < b)

            sort(Span(fill)[first:], narrower)
        var fill_next = Atomic[Int64](0)
        var left = Atomic[Int64](Int64(len(jobs)))
        var next = Atomic[Int64](0)
        var task = _StepTask(
            _Shared(sh.graph, sh.reference, sh.snarls, sh.refs, sh.heavy, sh.bad),
            imm_untracked(jobs.unsafe_ptr()),
            len(jobs),
            mut_untracked(Pointer(to=next)),
            mut_untracked(cur.ptr()),
            mut_untracked(Pointer(to=w)),
            mut_untracked(builders.unsafe_ptr()),
            level,
            threads,
            imm_untracked(fill_batch.unsafe_ptr()),
            imm_untracked(fill.unsafe_ptr()),
            len(fill),
            mut_untracked(Pointer(to=fill_next)),
            mut_untracked(Pointer(to=left)),
        )
        plan.timing.lap(4)
        plan.pool[].run(task, min(threads, len(jobs)), 0, plan.rank)
        plan.timing.lap(5)
        # The tasks reach these through untracked pointers: reads after them keep them alive until they are done.
        _ = len(jobs) + len(builders) + Int(next.load()) + chrom.byte_length()
        _ = len(fill) + len(fill_batch) + Int(fill_next.load()) + Int(left.load()) + len(futures)
        if plan.bad[0] != 0:
            plan.stop()
            raise Error(t"{chrom}: a record could not be built or encoded")
        _giants_done(giants, w)
        w.take_done()
        plan.timing.lap(6)
        if step < batches:
            for k in range(cur[].lo, cur[].hi):
                if plan.heavy_of[k] >= 0:
                    swap(cur[].records[k - cur[].lo], plan.heavy[plan.heavy_of[k]].record)
            records_written += _place(cur[], emitted, plan.parent, plan.snarls, plan.refs, plan.regions, chrom)
            for ref it in cur[].items:
                var record = Record()
                swap(record, cur[].records[it.record])
                var item = Placed(record^, it.lv, it.ps_start, it.ps_end, it.has_ps)
                if it.length <= STREAM_BYTES:
                    _ = w.place(item^, it.length)
                    continue
                # A giant record leaves the batch for a thread of its own.
                var unit = w.place_giant(item.record, it.length)
                giants.append(OwnedPointer(_GiantOut(unit)))
                giant_tasks.append(
                    OwnedPointer(
                        _GiantTask(sh.graph, item^, level, threads, mut_untracked(giants[len(giants) - 1].ptr()))
                    )
                )
                giant_jobs.append(plan.pool[].submit(giant_tasks[len(giant_tasks) - 1][], 1, 0, plan.rank))
        cur = futures.pop(0) if len(futures) > 0 else OwnedPointer(_Batch(0, 0))
        plan.timing.lap(4)
    return ContigTail(plan^, w^, giants^, giant_tasks^, giant_jobs^, records_written)


def finish_catalog(var tail: ContigTail) raises -> Tuple[Summary, List[EncodedBlock]]:
    """The giants still encoding awaited, the last blocks kept, the helpers stopped: the contig's blocks."""
    ref w = tail.writer
    for job in tail.giant_jobs:
        tail.plan.pool[].wait(job)
    for i in range(len(tail.giants)):
        if tail.giants[i][].state.load() != 1:
            raise Error(
                t"{tail.plan.reference[].chrom}: a giant record could not be encoded: {tail.giants[i][].error}"
            )
    _giants_done(tail.giants, w)
    w.take_done()
    if len(w.units) > 0:
        raise Error(t"{tail.plan.reference[].chrom}: blocks left unencoded")
    tail.plan.stop()
    tail.plan.timing.lap(7)
    tail.plan.timing.report(tail.plan.reference[].chrom)
    var blocks = List[EncodedBlock]()
    swap(blocks, w.blocks)
    var summary = Summary(
        len(tail.plan.snarls),
        len(tail.plan.refs),
        tail.records,
        len(blocks),
        Float64(perf_counter_ns() - tail.plan.t0) * 1e-9,
    )
    return (summary^, blocks^)


def build_catalog(
    graph: Graph,
    reference: ReferencePath,
    var snarls: List[Snarl],
    regions: List[Tuple[Int, Int]],
    level: Int,
    threads: Int,
    pool: MutPointer[Pool, MutUntrackedOrigin],
) raises -> Tuple[Summary, List[EncodedBlock]]:
    """The blocks of the reference's contig, or of `regions` of it (1-based, inclusive: only the records whose
    [POS, POS + len(REF) - 1] overlaps one, as `tabix contig:start-end` would return them): prepared, built and
    finished in turn."""
    var plan = ContigPlan(graph, reference.copy(), snarls^, regions, threads, pool, 0)
    return finish_catalog(run_steps(graph, plan^, level, threads))
