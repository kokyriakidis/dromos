"""One contig's site catalog, or one region of it, from its snarls (pgcatalog's catalog.rs, record.rs and
reference.rs, ported): every snarl the reference path passes through, at every nesting level, with all its alleles,
as `vg deconstruct -a` writes them, without what needs haplotype names (genotypes, AC/AF/AN/NS).

Pass 1 finds each snarl's reference traversal; snarls the reference does not pass through once have no record. The
rest are sorted by position, ancestors first at a tie, and built in batches on every thread: a snarl's haplotype walks
by forward GBWT search, its alleles (the reference, then distinct sequences, most haplotypes first), then LV/PS from
the ancestors built so far, and the lines written straight into the catalog's blocks, which are compressed in
parallel. A record whose line is longer than STREAM_BYTES is compressed as it is written, in a block of its own.

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
from dromos.parallel import Background, ParallelTask, parallel_for
from dromos.pgsites import FrameStream, Writer
from dromos.snarls import Snarl, canonical, put_visit, read_snarls
from dromos.stream import imm_untracked, mut_untracked, untracked

comptime NONE = UInt32.MAX
comptime BATCH = 50_000
"""Snarls built per batch: enough to keep every thread busy, few enough that a batch of records fits in memory."""
comptime STREAM_BYTES = 64 << 20
"""Records with longer lines are not formatted in memory: they are compressed from their walks as written."""
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


# Letters as alleles are written: upper case, or reverse-complemented.


@always_inline
def _upper(c: UInt8) -> UInt8:
    return c - 32 if c >= UInt8(ord("a")) and c <= UInt8(ord("z")) else c


@always_inline
def _complement(c: UInt8) -> UInt8:
    """The crate's COMPLEMENT: A, C, G, T in either case to their complement in capitals, any other letter to N."""
    var u = _upper(c)
    if u == UInt8(ord("A")):
        return UInt8(ord("T"))
    if u == UInt8(ord("C")):
        return UInt8(ord("G"))
    if u == UInt8(ord("G")):
        return UInt8(ord("C"))
    if u == UInt8(ord("T")):
        return UInt8(ord("A"))
    return UInt8(ord("N"))


comptime _W = 16
"""Letters converted at a time."""


@always_inline
def _upper_lanes(v: SIMD[DType.uint8, _W]) -> SIMD[DType.uint8, _W]:
    return (v.ge(UInt8(ord("a"))) & v.le(UInt8(ord("z")))).select(v - 32, v)


@always_inline
def _complement_lanes(v: SIMD[DType.uint8, _W]) -> SIMD[DType.uint8, _W]:
    comptime L = SIMD[DType.uint8, _W]
    var u = _upper_lanes(v)
    var out = L(UInt8(ord("N")))
    out = u.eq(UInt8(ord("A"))).select(L(UInt8(ord("T"))), out)
    out = u.eq(UInt8(ord("C"))).select(L(UInt8(ord("G"))), out)
    out = u.eq(UInt8(ord("G"))).select(L(UInt8(ord("C"))), out)
    return u.eq(UInt8(ord("T"))).select(L(UInt8(ord("A"))), out)


@always_inline
def put_oriented(dst: MutPointer[UInt8, _], src: Span[UInt8, _], reverse: Bool):
    """A node's letters as visited (pgcatalog's oriented_seq): upper case forward, reverse-complemented backward."""
    var n = len(src)
    var s = src.unsafe_ptr()
    var i = 0
    if reverse:
        while i + _W <= n:
            dst.unsafe_store(i, _complement_lanes(s.unsafe_load[width=_W](n - i - _W).reversed()))
            i += _W
        while i < n:
            dst[unsafe_offset=i] = _complement(s[unsafe_offset=n - 1 - i])
            i += 1
    else:
        while i + _W <= n:
            dst.unsafe_store(i, _upper_lanes(s.unsafe_load[width=_W](i)))
            i += _W
        while i < n:
            dst[unsafe_offset=i] = _upper(s[unsafe_offset=i])
            i += 1


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


trait LineSink:
    """Where a record's line is written."""

    def room(mut self) raises -> Tuple[MutPointer[UInt8, MutUntrackedOrigin], Int]:
        """Memory for the next bytes: where, and how many (at least one)."""
        ...

    def grew(mut self, n: Int) raises:
        """`n` bytes were written where `room` said."""
        ...


struct MemorySink(LineSink):
    """A line written into memory that has room for it."""

    var at: MutPointer[UInt8, MutUntrackedOrigin]

    def __init__(out self, at: MutPointer[UInt8, MutUntrackedOrigin]):
        self.at = at

    @always_inline
    def room(mut self) raises -> Tuple[MutPointer[UInt8, MutUntrackedOrigin], Int]:
        return (self.at, Int.MAX)

    @always_inline
    def grew(mut self, n: Int) raises:
        self.at = self.at.unsafe_offset(n)


struct StreamSink(LineSink):
    """A line compressed as it is written."""

    var stream: FrameStream

    def __init__(out self, var stream: FrameStream):
        self.stream = stream^

    def room(mut self) raises -> Tuple[MutPointer[UInt8, MutUntrackedOrigin], Int]:
        return (self.stream.end_ptr(), self.stream.room())

    def grew(mut self, n: Int) raises:
        self.stream.grew(n)


@always_inline
def _put[S: LineSink](mut out: S, src: Span[UInt8, _]) raises:
    var done = 0
    while done < len(src):
        var r = out.room()
        var n = min(r[1], len(src) - done)
        unsafe_memcpy(dest=r[0], src=src.unsafe_ptr().unsafe_offset(done), count=n)
        out.grew(n)
        done += n


@always_inline
def _put_byte[S: LineSink](mut out: S, b: UInt8) raises:
    var r = out.room()
    r[0][] = b
    out.grew(1)


@always_inline
def _put_number[S: LineSink](mut out: S, value: Int) raises:
    var digits = Array[UInt8, 24](fill=0)
    var n = put_padded(digits.unsafe_ptr(), value, 0)
    _put(out, Span(digits)[:n])


@always_inline
def _put_letters[S: LineSink](mut out: S, src: Span[UInt8, _], reverse: Bool) raises:
    """A node's letters as visited, piece by piece as the sink has room."""
    var n = len(src)
    var done = 0
    while done < n:
        var r = out.room()
        var k = min(r[1], n - done)
        if reverse:
            put_oriented(r[0], src[n - done - k : n - done], True)
        else:
            put_oriented(r[0], src[done : done + k], False)
        out.grew(k)
        done += k


def _put_allele[S: LineSink](graph: Graph, r: Record, k: Int, mut out: S) raises:
    """Allele k: the padding base, then its walk's inner nodes' letters (reverse-complemented when the record is)."""
    if r.pad >= 0:
        _put_byte(out, UInt8(r.pad))
    var a = r.starts[k]
    var b = r.ends[k]
    if b - a < 2:
        return
    if r.reversed:
        # The bases reverse-complemented: the inner nodes backwards, each in the opposite orientation.
        for i in range(b - 2, a, -1):
            var v = r.visits[i] ^ 1
            _put_letters(out, graph.letters(Int(v >> 1)), v & 1 == 1)
    else:
        for i in range(a + 1, b - 1):
            var v = r.visits[i]
            _put_letters(out, graph.letters(Int(v >> 1)), v & 1 == 1)


def _put_walk[S: LineSink](r: Record, k: Int, mut out: S) raises:
    """Walk k as the AT field writes it (backwards when the record is reversed)."""
    var a = r.starts[k]
    var b = r.ends[k]
    for i in range(b - a):
        var v = r.visits[b - 1 - i] if r.reversed else r.visits[a + i]
        _put_byte(out, UInt8(ord("<")) if (v & 1 == 1) != r.reversed else UInt8(ord(">")))
        _put_number(out, Int(v >> 1))


def _walk_text_length(r: Record, k: Int) -> Int:
    var n = 0
    for i in range(r.starts[k], r.ends[k]):
        n += 1 + decimal_digits(Int(r.visits[i] >> 1))
    return n


def write_line[S: LineSink](graph: Graph, r: Record, chrom: String, lv: Int, ps: Span[UInt8, _], mut out: S) raises:
    """The record's VCF line, without its newline (record.rs `write_line`)."""
    _put(out, chrom.as_bytes())
    _put_byte(out, UInt8(ord("\t")))
    _put_number(out, r.pos)
    _put_byte(out, UInt8(ord("\t")))
    var id = List[UInt8]()
    put_visit(id, r.start)
    put_visit(id, r.end)
    _put(out, Span(id))
    _put_byte(out, UInt8(ord("\t")))
    for k in range(r.walks()):
        if k > 1:
            _put_byte(out, UInt8(ord(",")))
        _put_allele(graph, r, k, out)
        if k == 0:
            _put_byte(out, UInt8(ord("\t")))
    _put(out, "\t60\t.\tAT=".as_bytes())
    for k in range(r.walks()):
        if k > 0:
            _put_byte(out, UInt8(ord(",")))
        _put_walk(r, k, out)
    _put(out, ";LV=".as_bytes())
    _put_number(out, lv)
    if len(ps) > 0:
        _put(out, ";PS=".as_bytes())
        _put(out, ps)


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
    """One thread's record building, its buffers kept from snarl to snarl."""

    var searcher: Searcher
    var found: WalkSet
    var order: List[Int]
    var text_a: List[UInt8]
    var text_b: List[UInt8]

    def __init__(out self):
        self.searcher = Searcher()
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
            prev_char = _complement(fp[0]) if fp_reverse else _upper(fp[len(fp) - 1])
        else:
            # The complement of its first letter.
            var first_letter = _complement(fp[len(fp) - 1]) if fp_reverse else _upper(fp[0])
            prev_char = _complement(first_letter)
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


def vcf_header(chrom: String, length: Int) -> String:
    return String(
        "##fileformat=VCFv4.2\n",
        '##INFO=<ID=AT,Number=R,Type=String,Description="Allele Traversal as path in graph">\n',
        '##INFO=<ID=LV,Number=1,Type=Integer,Description="Level in the snarl tree (0=top level)">\n',
        '##INFO=<ID=PS,Number=1,Type=String,Description="ID of variant corresponding to parent snarl">\n',
        "##contig=<ID=",
        chrom,
        ",length=",
        length,
        ">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
    )


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
struct _TraverseTask(ParallelTask):
    """Pass 1 on snarls `index * TRAVERSE_SPAN ..`: each one's reference traversal, if it has exactly one (and, for
    a region, can overlap it), into `found` (NO_TRAVERSAL's s0 otherwise)."""

    var graph: ImmPointer[Graph, ImmUntrackedOrigin]
    var reference: ImmPointer[ReferencePath, ImmUntrackedOrigin]
    var snarls: ImmPointer[Snarl, ImmUntrackedOrigin]
    var count: Int
    var found: MutPointer[Traversal, MutUntrackedOrigin]
    var region: Bool
    var beg: Int
    var end: Int
    var bad: MutPointer[Int, MutUntrackedOrigin]

    def run(self, index: Int):
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
            if self.region:
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
                if r.pos[left] > self.end or span_end + 1 < self.beg:
                    continue
            self.found[unsafe_offset=i] = t


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
    heavy snarl's record awaited (`a`), another record built (`a`: the batch's record), or a unit's lines written and
    the unit compressed if sealed (`a`: the unit)."""

    var kind: Int
    var a: Int


@fieldwise_init
struct _Item(Copyable, TrivialRegisterPassable):
    """A record to write: its index in the batch, LV, PS (bytes in the batch's PS text; length 0 for none), its
    line's length, and its unit and offset there (-1 for a giant)."""

    var record: Int
    var lv: Int
    var ps_at: Int
    var ps_length: Int
    var length: Int
    var unit: Int
    var offset: Int


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


struct _HelperTask(ParallelTask):
    """Helper threads: the heavy snarls' parts, in position order, ahead of the steps that need them."""

    var sh: _Shared
    var count: Int
    var next: MutPointer[Atomic[Int64], MutUntrackedOrigin]
    var stop: MutPointer[Atomic[Int64], MutUntrackedOrigin]
    var builders: MutPointer[Builder, MutUntrackedOrigin]

    def __init__(
        out self,
        var sh: _Shared,
        count: Int,
        next: MutPointer[Atomic[Int64], MutUntrackedOrigin],
        stop: MutPointer[Atomic[Int64], MutUntrackedOrigin],
        builders: MutPointer[Builder, MutUntrackedOrigin],
    ):
        self.sh = sh^
        self.count = count
        self.next = next
        self.stop = stop
        self.builders = builders

    def run(self, worker: Int):
        ref b = self.builders[unsafe_offset=worker]
        while self.stop[].load() == 0:
            var h = Int(self.next[].fetch_add(1))
            if h >= self.count:
                return
            try:
                while _take_part(self.sh, b, h):
                    pass
            except:
                self.sh.bad[] = 1
                return


struct _GiantOut(Movable):
    """A giant record's block, made on a thread of its own: its frame, and 1 once made (2 if that failed)."""

    var frame: List[UInt8]
    var state: Atomic[Int64]
    var unit: Int
    var error: String

    def __init__(out self, unit: Int):
        self.frame = List[UInt8]()
        self.state = Atomic[Int64](0)
        self.unit = unit
        self.error = String()


struct _GiantTask(ParallelTask):
    """A record whose line is longer than STREAM_BYTES, compressed as it is written (zstd's own threads compress)."""

    var graph: ImmPointer[Graph, ImmUntrackedOrigin]
    var record: Record
    var chrom: String
    var lv: Int
    var ps: List[UInt8]
    var length: Int
    var level: Int
    var threads: Int
    var result: MutPointer[_GiantOut, MutUntrackedOrigin]

    def __init__(
        out self,
        graph: ImmPointer[Graph, ImmUntrackedOrigin],
        var record: Record,
        chrom: String,
        lv: Int,
        var ps: List[UInt8],
        length: Int,
        level: Int,
        threads: Int,
        result: MutPointer[_GiantOut, MutUntrackedOrigin],
    ):
        self.graph = graph
        self.record = record^
        self.chrom = chrom
        self.lv = lv
        self.ps = ps^
        self.length = length
        self.level = level
        self.threads = threads
        self.result = result

    def run(self, index: Int):
        try:
            var sink = StreamSink(FrameStream(self.length + 1, self.level, self.threads))
            write_giant(self.graph[], self.record, self.chrom, self.lv, Span(self.ps), sink, self.threads)
            var newline = sink.room()
            newline[0][] = UInt8(ord("\n"))
            sink.grew(1)
            self.result[].frame = sink.stream.finish()
            self.result[].state.store(1)
        except e:
            self.result[].error = String(e)
            self.result[].state.store(2)


comptime GIANT_WINDOW = 256 << 20
"""The most of a giant record's line made at once."""


@fieldwise_init
struct _PieceTask(ParallelTask):
    """Pieces `first + index` of a giant record's line (allele k for piece k, then walk k), each into its buffer."""

    var graph: ImmPointer[Graph, ImmUntrackedOrigin]
    var record: ImmPointer[Record, ImmUntrackedOrigin]
    var first: Int
    var buffers: MutPointer[List[UInt8], MutUntrackedOrigin]
    var bad: MutPointer[Int, MutUntrackedOrigin]

    def run(self, index: Int):
        var piece = self.first + index
        ref r = self.record[]
        var n = r.walks()
        var k = piece if piece < n else piece - n
        var length = r.lens[k] if piece < n else _walk_text_length(r, k)
        ref buf = self.buffers[unsafe_offset=index]
        buf = List[UInt8](unsafe_uninit_length=length)
        var sink = MemorySink(mut_untracked(buf.unsafe_ptr()))
        try:
            if piece < n:
                _put_allele(self.graph[], r, k, sink)
            else:
                _put_walk(r, k, sink)
        except:
            self.bad[] = 1
            return
        if Int(sink.at) - Int(buf.unsafe_ptr()) != length:
            self.bad[] = 1


def write_giant(
    graph: Graph, r: Record, chrom: String, lv: Int, ps: Span[UInt8, _], mut out: StreamSink, threads: Int
) raises:
    """`write_line` for a giant record: its alleles and walks made on `threads` threads, a window of them at a time,
    and handed to zstd in order."""
    _put(out, chrom.as_bytes())
    _put_byte(out, UInt8(ord("\t")))
    _put_number(out, r.pos)
    _put_byte(out, UInt8(ord("\t")))
    var id = List[UInt8]()
    put_visit(id, r.start)
    put_visit(id, r.end)
    _put(out, Span(id))
    _put_byte(out, UInt8(ord("\t")))
    var n = r.walks()
    var pieces = 2 * n
    var first = 0
    while first < pieces:
        var last = first
        var bytes = 0
        while last < pieces and (last == first or bytes < GIANT_WINDOW):
            bytes += r.lens[last] if last < n else _walk_text_length(r, last - n)
            last += 1
        var buffers = List[List[UInt8]](length=last - first, fill=List[UInt8]())
        var bad = 0
        var task = _PieceTask(
            untracked(Pointer(to=graph)),
            untracked(Pointer(to=r)),
            first,
            mut_untracked(buffers.unsafe_ptr()),
            mut_untracked(Pointer(to=bad)),
        )
        parallel_for(task, last - first, threads)
        if bad != 0:
            raise Error("a giant record's line could not be made")
        for i in range(first, last):
            # The separator before the piece: the REF ends at a tab, ALTs are split by commas, then the AT field.
            if i == 1:
                _put_byte(out, UInt8(ord("\t")))
            elif i > 1 and i < n:
                _put_byte(out, UInt8(ord(",")))
            elif i == n:
                if n == 1:
                    _put_byte(out, UInt8(ord("\t")))
                _put(out, "\t60\t.\tAT=".as_bytes())
            elif i > n:
                _put_byte(out, UInt8(ord(",")))
            ref buf = buffers[i - first]
            out.stream.feed(buf.unsafe_ptr(), len(buf))
        first = last
    _put(out, ";LV=".as_bytes())
    _put_number(out, lv)
    if len(ps) > 0:
        _put(out, ";PS=".as_bytes())
        _put(out, ps)


struct _Batch(Movable):
    """A batch of snarls: their records, then, once built, the items to write and where."""

    var lo: Int
    var hi: Int
    var records: List[Record]
    var items: List[_Item]
    var ps_text: List[UInt8]

    def __init__(out self, lo: Int, hi: Int):
        self.lo = lo
        self.hi = hi
        self.records = List[Record]()
        for _ in range(hi - lo):
            self.records.append(Record())
        self.items = List[_Item]()
        self.ps_text = List[UInt8]()


struct _StepTask(ParallelTask):
    """One step's jobs, taken in turn by `workers` workers (each with its own builder): the batch's records built
    (`cur`), the previous batch's lines written and its sealed units compressed (`prev`)."""

    var sh: _Shared
    var jobs: ImmPointer[_Job, ImmUntrackedOrigin]
    var count: Int
    var next: MutPointer[Atomic[Int64], MutUntrackedOrigin]
    var cur: MutPointer[_Batch, MutUntrackedOrigin]
    var prev: ImmPointer[_Batch, ImmUntrackedOrigin]
    var unit_lines: ImmPointer[List[Int], ImmUntrackedOrigin]
    """Per unit (from the writer's base), the previous batch's items placed in it."""
    var writer: MutPointer[Writer, MutUntrackedOrigin]
    var builders: MutPointer[Builder, MutUntrackedOrigin]
    var chrom: ImmPointer[String, ImmUntrackedOrigin]
    var level: Int
    var threads: Int

    def __init__(
        out self,
        var sh: _Shared,
        jobs: ImmPointer[_Job, ImmUntrackedOrigin],
        count: Int,
        next: MutPointer[Atomic[Int64], MutUntrackedOrigin],
        cur: MutPointer[_Batch, MutUntrackedOrigin],
        prev: ImmPointer[_Batch, ImmUntrackedOrigin],
        unit_lines: ImmPointer[List[Int], ImmUntrackedOrigin],
        writer: MutPointer[Writer, MutUntrackedOrigin],
        builders: MutPointer[Builder, MutUntrackedOrigin],
        chrom: ImmPointer[String, ImmUntrackedOrigin],
        level: Int,
        threads: Int,
    ):
        self.sh = sh^
        self.jobs = jobs
        self.count = count
        self.next = next
        self.cur = cur
        self.prev = prev
        self.unit_lines = unit_lines
        self.writer = writer
        self.builders = builders
        self.chrom = chrom
        self.level = level
        self.threads = threads

    def run(self, worker: Int):
        ref b = self.builders[unsafe_offset=worker]
        while True:
            var i = Int(self.next[].fetch_add(1))
            if i >= self.count:
                return
            try:
                self._job(b, self.jobs[unsafe_offset=i])
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
            ref w = self.writer[]
            var u = w.unit(job.a)
            ref prev = self.prev[]
            for k in self.unit_lines[unsafe_offset=job.a - w.base]:
                var it = prev.items[k]
                var sink = MemorySink(w.at(job.a, it.offset))
                var start = Int(sink.at)
                write_line(
                    g,
                    prev.records[it.record],
                    self.chrom[],
                    it.lv,
                    Span(prev.ps_text)[it.ps_at : it.ps_at + it.ps_length],
                    sink,
                )
                if Int(sink.at) - start != it.length:
                    raise Error("a record's line is not the length it was measured")
            if u[].sealed and not u[].done:
                u[].compress(self.level, self.threads)


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
            "blocks written",
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


def _place(
    mut batch: _Batch,
    mut emitted: List[Bool],
    parent: List[UInt32],
    snarls: List[Snarl],
    refs: List[RefTrav],
    region: Optional[Tuple[Int, Int]],
    chrom: String,
    mut w: Writer,
) -> Int:
    """A built batch's records to write, in file order, each placed in its unit: LV / PS from the ancestors built so
    far, as vg's update_nesting_info_tags (the ancestors that have a record count, the nearest names the parent).
    Returns how many."""
    for k in range(batch.lo, batch.hi):
        emitted[refs[k].snarl] = batch.records[k - batch.lo].present
    batch.items.clear()
    batch.ps_text.clear()
    for k in range(batch.lo, batch.hi):
        ref r = batch.records[k - batch.lo]
        if not r.present:
            continue
        if region and (r.pos > region.value()[1] or r.last() < region.value()[0]):
            continue
        var lv = 0
        var ps_at = len(batch.ps_text)
        var ps_length = 0
        var c = parent[refs[k].snarl]
        while c != NONE:
            var pi = Int(c)
            if emitted[pi]:
                lv += 1
                if ps_length == 0:
                    put_visit(batch.ps_text, snarls[pi].start)
                    put_visit(batch.ps_text, snarls[pi].end)
                    ps_length = len(batch.ps_text) - ps_at
            c = parent[pi]
        var length = r.line_length(chrom.byte_length(), lv, ps_length)
        batch.items.append(_Item(k - batch.lo, lv, ps_at, ps_length, length, 0, 0))
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
    for ref it in batch.items:
        ref r = batch.records[it.record]
        if it.length > STREAM_BYTES:
            it.unit = w.place_giant(r.pos, r.last(), r.min_allele(), r.max_walk(), it.length + 1)
            it.offset = -1
        else:
            var place = w.place(r.pos, r.last(), r.min_allele(), r.max_walk(), it.length)
            it.unit = place[0]
            it.offset = place[1]
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
    var region: Optional[Tuple[Int, Int]]
    var helper_builders: List[Builder]
    var counters: Allocation[Atomic[Int64]]
    """The helpers' next heavy snarl, and their stop."""
    var counters_ptr: MutPointer[Atomic[Int64], MutUntrackedOrigin]
    var bad: List[Int]
    """1 once a task failed."""
    var helpers: Background[_HelperTask]
    var timing: _Phases
    var t0: Int

    def __init__(
        out self,
        graph: Graph,
        var reference: ReferencePath,
        snarls_path: String,
        region: Optional[Tuple[Int, Int]],
        threads: Int,
    ) raises:
        self.t0 = perf_counter_ns()
        self.timing = _Phases()
        self.reference = OwnedPointer(reference^)
        self.region = region
        self.bad = List[Int](length=1, fill=0)
        self.helper_builders = List[Builder]()
        self.heavy = List[_Heavy]()
        self.snarls = read_snarls(snarls_path, threads)
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
        # need one record per crossing) and, for a region, when it can overlap.
        var found = List[Traversal](unsafe_uninit_length=n)
        var pass1 = _TraverseTask(
            untracked(Pointer(to=graph)),
            rp,
            imm_untracked(self.snarls.unsafe_ptr()),
            n,
            mut_untracked(found.unsafe_ptr()),
            Bool(region),
            region.value()[0] if region else 0,
            region.value()[1] if region else 0,
            mut_untracked(self.bad.unsafe_ptr()),
        )
        parallel_for(pass1, (n + TRAVERSE_SPAN - 1) // TRAVERSE_SPAN, threads)
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
        for k in range(len(self.refs)):
            var t = self.refs[k].t
            var sn = self.snarls[self.refs[k].snarl]
            var span = abs(r.pos[t.s0] - r.pos[t.e0])
            var nodes = abs(Int(sn.end >> 1) - Int(sn.start >> 1))
            if span >= HEAVY_SPAN or nodes >= HEAVY_NODES:
                var first = graph.parts.search_state(Int(sn.start >> 1), sn.start & 1 == 1)
                var haplotypes = first.value().count() if first else 0
                self.heavy_of[k] = len(self.heavy)
                self.heavy.append(_Heavy(k, _parts(span, nodes, haplotypes, threads)))
        var helpers_n = max(threads // 2, 1)
        for _ in range(helpers_n):
            self.helper_builders.append(Builder())
        self.counters = alloc(Layout[Atomic[Int64]](length=2))
        self.counters_ptr = mut_untracked(self.counters.unsafe_ptr())
        self.counters_ptr.unsafe_write(Atomic[Int64](0))
        self.counters_ptr.unsafe_offset(1).unsafe_write(Atomic[Int64](0))
        self.helpers = Background(
            _HelperTask(
                _Shared(
                    untracked(Pointer(to=graph)),
                    rp,
                    untracked(self.snarls.unsafe_ptr()),
                    untracked(self.refs.unsafe_ptr()),
                    mut_untracked(self.heavy.unsafe_ptr()),
                    mut_untracked(self.bad.unsafe_ptr()),
                ),
                len(self.heavy),
                self.counters_ptr,
                self.counters_ptr.unsafe_offset(1),
                mut_untracked(self.helper_builders.unsafe_ptr()),
            ),
            helpers_n,
            helpers_n if len(self.heavy) > 0 else 0,
        )
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
        self.counters_ptr.unsafe_offset(1)[].store(1)
        self.helpers.join()

    def __deinit__(deinit self):
        self.counters_ptr.unsafe_offset(1)[].store(1)
        self.helpers.join()
        dealloc(self.counters^)


struct ContigTail(Movable):
    """A catalog whose records are all placed: the giant records still compressing on threads of their own, and the
    blocks left to write."""

    var plan: ContigPlan
    var writer: Writer
    var giants: List[OwnedPointer[_GiantOut]]
    var giant_threads: List[Background[_GiantTask]]
    var records: Int

    def __init__(
        out self,
        var plan: ContigPlan,
        var writer: Writer,
        var giants: List[OwnedPointer[_GiantOut]],
        var giant_threads: List[Background[_GiantTask]],
        records: Int,
    ):
        self.plan = plan^
        self.writer = writer^
        self.giants = giants^
        self.giant_threads = giant_threads^
        self.records = records


def _giants_done(mut giants: List[OwnedPointer[_GiantOut]], mut w: Writer):
    """The frames of the giants made since, to their units."""
    for i in range(len(giants)):
        if giants[i][].state.load() == 1 and giants[i][].unit >= w.base:
            var u = w.unit(giants[i][].unit)
            if not u[].done:
                swap(u[].frame, giants[i][].frame)
                u[].done = True


def run_steps(graph: Graph, var plan: ContigPlan, out_prefix: String, level: Int, threads: Int) raises -> ContigTail:
    """A prepared catalog's records built, placed, written and compressed, in steps: a step's tasks build one batch's
    records and, at the same time, write the previous batch's lines into their blocks and compress the sealed ones;
    between steps the batch just built is placed (one thread) and the blocks done are written, in order. Heavy
    snarls' parts no helper has taken are searched in the step that needs them; giant records are compressed on
    threads of their own."""
    plan.timing.lap(8)
    var sh = plan.shared(graph)
    var batches = len(plan.bounds) - 1
    ref r = plan.reference[]
    var chrom = r.chrom
    var w = Writer(out_prefix + ".pgsites", vcf_header(chrom, r.length), level, threads)
    var emitted = List[Bool](length=len(plan.snarls), fill=False)
    var builders = List[Builder]()
    for _ in range(threads):
        builders.append(Builder())
    var giants = List[OwnedPointer[_GiantOut]]()
    var giant_threads = List[Background[_GiantTask]]()
    var records_written = 0
    var prev = _Batch(0, 0)
    for step in range(batches + 1):
        var cur = _Batch(plan.bounds[step], plan.bounds[step + 1]) if step < batches else _Batch(0, 0)
        if step == batches:
            w.seal()
        # The jobs: heavy snarls' parts no helper has taken first, then units, the other records, and last the heavy
        # snarls' records awaited (a helper may still be searching a part).
        var jobs = List[_Job]()
        var unit_lines = List[List[Int]](length=len(w.units), fill=List[Int]())
        for k in range(len(prev.items)):
            var it = prev.items[k]
            if it.offset >= 0:
                unit_lines[it.unit - w.base].append(k)
        for k in range(cur.lo, cur.hi):
            var h = plan.heavy_of[k]
            if h >= 0:
                var taken = Int(plan.heavy[h].state_ptr[unsafe_offset=_TAKEN].load())
                for _ in range(len(plan.heavy[h].ranges) - taken):
                    jobs.append(_Job(_PART, h))
        for u in range(len(w.units)):
            ref unit = w.units[u]
            if not unit.giant and (len(unit_lines[u]) > 0 or (unit.sealed and not unit.done)):
                jobs.append(_Job(_UNIT, w.base + u))
        for k in range(cur.lo, cur.hi):
            if plan.heavy_of[k] < 0:
                jobs.append(_Job(_BUILD, k - cur.lo))
        for k in range(cur.lo, cur.hi):
            if plan.heavy_of[k] >= 0:
                jobs.append(_Job(_HEAVY, plan.heavy_of[k]))
        var next = Atomic[Int64](0)
        var task = _StepTask(
            _Shared(sh.graph, sh.reference, sh.snarls, sh.refs, sh.heavy, sh.bad),
            imm_untracked(jobs.unsafe_ptr()),
            len(jobs),
            mut_untracked(Pointer(to=next)),
            mut_untracked(Pointer(to=cur)),
            imm_untracked(Pointer(to=prev)),
            imm_untracked(unit_lines.unsafe_ptr()),
            mut_untracked(Pointer(to=w)),
            mut_untracked(builders.unsafe_ptr()),
            imm_untracked(Pointer(to=chrom)),
            level,
            threads,
        )
        plan.timing.lap(4)
        parallel_for(task, min(threads, len(jobs)), threads)
        plan.timing.lap(5)
        # The tasks reach these through untracked pointers: reads after them keep them alive until they are done.
        _ = len(jobs) + len(unit_lines) + len(builders) + Int(next.load()) + len(prev.items) + chrom.byte_length()
        if plan.bad[0] != 0:
            plan.stop()
            raise Error(t"{chrom}: a record could not be built or written")
        _giants_done(giants, w)
        w.write_done()
        plan.timing.lap(6)
        if step < batches:
            for k in range(cur.lo, cur.hi):
                if plan.heavy_of[k] >= 0:
                    swap(cur.records[k - cur.lo], plan.heavy[plan.heavy_of[k]].record)
            records_written += _place(cur, emitted, plan.parent, plan.snarls, plan.refs, plan.region, chrom, w)
            # Giant records leave the batch for threads of their own.
            for ref it in cur.items:
                if it.offset >= 0:
                    continue
                giants.append(OwnedPointer(_GiantOut(it.unit)))
                var ps = List[UInt8](Span(cur.ps_text)[it.ps_at : it.ps_at + it.ps_length])
                var record = Record()
                swap(record, cur.records[it.record])
                giant_threads.append(
                    Background(
                        _GiantTask(
                            sh.graph,
                            record^,
                            chrom,
                            it.lv,
                            ps^,
                            it.length,
                            level,
                            threads,
                            mut_untracked(giants[len(giants) - 1].ptr()),
                        ),
                        1,
                        1,
                    )
                )
                if len(giant_threads[len(giant_threads) - 1].handles) == 0:
                    raise Error("cannot start a thread")
        prev = cur^
        plan.timing.lap(4)
    return ContigTail(plan^, w^, giants^, giant_threads^, records_written)


def finish_catalog(var tail: ContigTail) raises -> Summary:
    """The giants still compressing awaited, the last blocks and the index written, the file synced and renamed."""
    ref w = tail.writer
    for i in range(len(tail.giant_threads)):
        tail.giant_threads[i].join()
    for i in range(len(tail.giants)):
        if tail.giants[i][].state.load() != 1:
            raise Error(
                t"{tail.plan.reference[].chrom}: a giant record could not be compressed: {tail.giants[i][].error}"
            )
    _giants_done(tail.giants, w)
    w.write_done()
    tail.plan.stop()
    var blocks = w.finish()
    tail.plan.timing.lap(7)
    tail.plan.timing.report(tail.plan.reference[].chrom)
    return Summary(
        len(tail.plan.snarls),
        len(tail.plan.refs),
        tail.records,
        blocks,
        Float64(perf_counter_ns() - tail.plan.t0) * 1e-9,
    )


def build_catalog(
    graph: Graph,
    reference: ReferencePath,
    snarls_path: String,
    region: Optional[Tuple[Int, Int]],
    out_prefix: String,
    level: Int,
    threads: Int,
) raises -> Summary:
    """The catalog `<out_prefix>.pgsites` of the reference's contig, or of `region` of it (1-based, inclusive: only
    the records whose [POS, POS + len(REF) - 1] overlaps it, as `tabix contig:start-end` would return them): prepared,
    built and finished in turn (`dromos catalog` overlaps one selection's with the next's)."""
    var plan = ContigPlan(graph, reference.copy(), snarls_path, region, threads)
    return finish_catalog(run_steps(graph, plan^, out_prefix, level, threads))
