"""Tests: protobuf varints, GAM group framing, GAF lines against hand-worked ones, placement on the reference as
pgcat places alignments, tabix bins as htslib numbers them, surject's alignments as pgcatalog makes them, a GBZ read
through the gbz crate (tests/data/micb-kir3dl1.gbz, two reference samples), snarls streams, and site catalogs as
pgcatalog writes and reads them (tests/data/micb-kir3dl1.catalog: its snarls and pgcatalog's catalogs)."""

from std.memory import bitcast
from std.testing import TestSuite, assert_equal, assert_false, assert_true

from dromos.gam import Alignment, Framer, frame, gaf_line, slim_line, varint
from dromos.nodes import NodeTable, complement
from dromos.cache import (
    CURRENT,
    DAMAGED,
    MISSING,
    NONE,
    NOT_A_CACHE,
    OTHER_VERSION,
    CacheFile,
    CacheWriter,
    Header,
    examine,
)
from dromos.graph import (
    GbzGraph,
    GbzParts,
    GraphId,
    NodeSource,
    SearchState,
    fingerprint,
    graph_id,
    is_gbz,
    reference_paths,
    reference_samples,
    resolve_reference_sample,
)
from dromos.prepare import ReferenceSample, folder_of, make_cache, name_of, pack_letters, write_nodes, write_placements
from dromos.build import natural_less, parse_region
from dromos.catalog import Graph, ReferencePath, Traversal, build_catalog, traversals, vcf_header
from dromos.pgsites import Catalog, Writer, compress_sealed, records, same_contig
from dromos.snarls import NO_PARENT, Snarl, canonical, parse_snarls, read_snarls, snarl_name
from dromos.gaf_text import MALFORMED, PLACED, SKIP, UNMAPPED, append_slim, classify
from dromos.output import decimal_digits
from dromos.bgzf import load_libdeflate
from dromos.stream import BGZF, GZIP, PLAIN, looks_like_gaf
from dromos.refmap import AMBIGUOUS, INCONSISTENT, MISSING_NODE, NO_ANCHOR, Anchor, Placed, RefMap
from dromos.tabix import N_BINS, bin_level, bin_parent, reg2bin
from dromos.linearize import (
    CLIP_MATCH,
    CLIP_MISMATCH,
    D,
    I,
    M,
    S,
    Block,
    CigarOp,
    Window,
    align_end,
    align_global,
    chain,
    clip_cigar,
    clip_gap,
    clip_lengths,
    linearize,
)


def _varint(mut out: List[UInt8], value: Int):
    var v = value
    while v >= 0x80:
        out.append(UInt8(v & 0x7F | 0x80))
        v >>= 7
    out.append(UInt8(v))


def _bytes_field(mut out: List[UInt8], field: Int, data: List[UInt8]):
    _varint(out, field << 3 | 2)
    _varint(out, len(data))
    out.extend(Span(data))


def _int_field(mut out: List[UInt8], field: Int, value: Int):
    _varint(out, field << 3)
    _varint(out, value)


def _text(s: StaticString) -> List[UInt8]:
    return List[UInt8](s.as_bytes())


def _edit(from_length: Int, to_length: Int, sequence: StaticString) -> List[UInt8]:
    var e = List[UInt8]()
    _int_field(e, 1, from_length)
    _int_field(e, 2, to_length)
    if sequence.byte_length() > 0:
        _bytes_field(e, 3, _text(sequence))
    return e^


def _mapping(node: Int, offset: Int, var edits: List[List[UInt8]]) -> List[UInt8]:
    var position = List[UInt8]()
    _int_field(position, 1, node)
    if offset > 0:
        _int_field(position, 2, offset)
    var m = List[UInt8]()
    _bytes_field(m, 1, position)
    for e in edits:
        _bytes_field(m, 2, e)
    return m^


struct _MemoryNodes(NodeSource):
    """Nodes held in memory: lengths (NONE for an id the graph lacks) and their letters one after another."""

    var lengths: List[UInt32]
    var letters: List[UInt8]

    def __init__(out self, var lengths: List[UInt32], letters: StaticString):
        self.lengths = lengths^
        self.letters = _text(letters)

    def _length(self, node: Int) -> Int:
        if node >= len(self.lengths) or self.lengths[node] == UInt32(NONE):
            return 0
        return Int(self.lengths[node])

    def lengths_into(
        self, first: Int, count: Int, dst: MutPointer[UInt32, MutUntrackedOrigin]
    ) raises -> Tuple[Int, Int]:
        var total = 0
        for i in range(count):
            dst[unsafe_offset=i] = self.lengths[first + i]
            total += self._length(first + i)
        return (total, -1)

    def letters_into(
        self, first: Int, count: Int, dst: MutPointer[UInt8, MutUntrackedOrigin], capacity: Int
    ) raises -> Int:
        var at = 0
        for i in range(first):
            at += self._length(i)
        var n = 0
        for i in range(first, first + count):
            n += self._length(i)
        for i in range(n):
            dst[unsafe_offset=i] = self.letters[at + i]
        return n


def _cache(
    path: String, var lengths: List[UInt32], letters: StaticString, var samples: List[ReferenceSample], width: Int = 4
) raises -> String:
    """A graph cache of the nodes and reference samples given, as `prepare` writes one."""
    var count = len(lengths)
    var graph = _MemoryNodes(lengths^, letters)
    var w = CacheWriter(path)
    var nodes = write_nodes(w, graph, count, width, 2)
    write_placements(w, samples, nodes.lengths(), count, 2)
    w.finish(GraphId(1, 2, 3), "test.gbz")
    _ = nodes.count
    return path


def _nodes_file(width: Int = 4) raises -> String:
    """A graph cache of nodes 1 ACGTACGT, 2 GGG and 3 TTTTT (node 0 absent), no reference sample."""
    var lengths: List[UInt32] = [UInt32(NONE), 8, 3, 5]
    return _cache(
        String(t"/tmp/dromos_test_nodes_{width}.bin"), lengths^, "ACGTACGTGGGTTTTT", List[ReferenceSample](), width
    )


def test_varint() raises:
    var b: List[UInt8] = [0x96, 0x01, 0x05]
    var at = 0
    assert_equal(Int(varint(b.unsafe_ptr(), at, len(b))), 150)
    assert_equal(at, 2)
    assert_equal(Int(varint(b.unsafe_ptr(), at, len(b))), 5)


def test_frame_resumes_inside_a_group() raises:
    # Two groups: 2 messages, then 1. Cut inside the first group's second message, the first message is handed out
    # and the rest read again with more data, the framer still inside the group.
    var stream = List[UInt8]()
    for count in [2, 1]:
        _varint(stream, count + 1)
        _varint(stream, 3)
        stream.extend(Span(_text("GAM")))
        for _ in range(count):
            _varint(stream, 3)
            stream.extend(Span(_text("abc")))
    var messages = List[Int]()
    assert_equal(frame(Span(stream), messages), len(stream))
    assert_equal(len(messages), 6)
    var f = Framer()
    messages.clear()
    var cut = 2 + 3 + 4 + 2
    var used = f.frame(Span(stream)[:cut], messages)
    assert_equal(used, 2 + 3 + 4)
    assert_equal(len(messages), 2)
    var rest = List[UInt8](Span(stream)[used:])
    messages.clear()
    assert_equal(f.frame(Span(rest), messages), len(rest))
    assert_equal(len(messages), 4)
    assert_equal(messages[0], 1)
    assert_equal(f.remaining, 0)


def test_frame_skips_other_tags_and_reads_untagged_groups() raises:
    # A PARAMS group (skipped), then an untagged group of two messages (an older GAM).
    var stream = List[UInt8]()
    _varint(stream, 2)
    _varint(stream, 6)
    stream.extend(Span(_text("PARAMS")))
    _varint(stream, 2)
    stream.extend(Span(_text("{}")))
    _varint(stream, 2)
    for _ in range(2):
        _varint(stream, 40)
        for _ in range(40):
            stream.append(0x0A)
    var messages = List[Int]()
    assert_equal(frame(Span(stream), messages), len(stream))
    assert_equal(len(messages), 4)
    assert_equal(messages[1] - messages[0], 40)


def test_gaf_line_as_vg_writes_it() raises:
    var nodes = NodeTable(_nodes_file())
    # r1 over >1>2>3: 4 matches, A->T, 3 matches on node 1; node 2 whole; 2 of node 3's 5.
    var edits1: List[List[UInt8]] = [_edit(4, 4, ""), _edit(1, 1, "T"), _edit(3, 3, "")]
    var edits2: List[List[UInt8]] = [_edit(3, 3, "")]
    var edits3: List[List[UInt8]] = [_edit(2, 2, "")]
    var path = List[UInt8]()
    _bytes_field(path, 2, _mapping(1, 0, edits1^))
    _bytes_field(path, 2, _mapping(2, 0, edits2^))
    _bytes_field(path, 2, _mapping(3, 0, edits3^))
    var msg = List[UInt8]()
    _bytes_field(msg, 1, _text("ACGTTCGTGGGTT"))
    _bytes_field(msg, 2, path)
    _bytes_field(msg, 3, _text("r1"))
    var quality = List[UInt8](length=13, fill=30)
    _bytes_field(msg, 4, quality)
    _int_field(msg, 5, 60)
    _int_field(msg, 6, 9)
    _varint(msg, 16 << 3 | 1)
    var bits = bitcast[DType.uint64](Float64(0.9876))
    for k in range(8):
        msg.append(UInt8((bits >> UInt64(8 * k)) & 0xFF))
    var a = Alignment()
    a.decode(msg.unsafe_ptr(), 0, len(msg))
    var out = List[UInt8]()
    gaf_line(msg.unsafe_ptr(), a, nodes, out)
    var want = _text("r1\t13\t0\t13\t+\t>1>2>3\t16\t0\t13\t12\t13\t60\tAS:i:9\tbq:Z:?????????????")
    want.extend(Span(_text("\tcs:Z::4*AT:8\tdv:f:0.0124\n")))
    assert_true(out == want)
    # --slim: the 12 columns alone, from the full decode and from the slim one, with the walk and path span.
    var slim = List[UInt8]()
    gaf_line[False](msg.unsafe_ptr(), a, nodes, slim)
    assert_true(slim == _text("r1\t13\t0\t13\t+\t>1>2>3\t16\t0\t13\t12\t13\t60\n"))
    var b = Alignment()
    b.decode[slim=True](msg.unsafe_ptr(), 0, len(msg))
    var fast = List[UInt8]()
    var walk = List[Int]()
    var span = slim_line(msg.unsafe_ptr(), b, nodes, fast, walk)
    assert_true(fast == slim)
    assert_equal(len(walk), 3)
    assert_equal(walk[2], 3 << 1)
    assert_true(span[0] == 0 and span[1] == 13)


def test_unaligned_line() raises:
    var nodes = NodeTable(_nodes_file())
    var msg = List[UInt8]()
    _bytes_field(msg, 1, _text("ACG"))
    _bytes_field(msg, 3, _text("u"))
    var a = Alignment()
    a.decode(msg.unsafe_ptr(), 0, len(msg))
    var out = List[UInt8]()
    gaf_line(msg.unsafe_ptr(), a, nodes, out)
    assert_true(out == _text("u\t3\t0\t3\t*\t*\t*\t*\t*\t*\t*\t0\tcs:Z:+ACG\n"))


def _sample(
    name: String,
    var walks: List[List[UInt64]],
    var fragments: List[Int],
    var path_contig: List[Int],
    var contigs: List[String],
) -> ReferenceSample:
    var s = ReferenceSample(name)
    s.walks = walks^
    s.fragments = fragments^
    s.path_contig = path_contig^
    s.contigs = contigs^
    return s^


def _refmap_file(twice: Bool) raises -> String:
    """pgcat's test map: reference path >1 >2 >3 (lengths 10, 5, 10) at 100 on chrA; node 4 (5 bp) off the reference;
    node 0 absent. With `twice`, node 2 is also on a second path, at 500 of chrB."""
    var walks: List[List[UInt64]] = [[2, 4, 6]]
    var fragments: List[Int] = [100]
    var path_contig: List[Int] = [0]
    if twice:
        walks.append([4])
        fragments.append(500)
        path_contig.append(1)
    var samples = List[ReferenceSample]()
    samples.append(_sample("CHM13", walks^, fragments^, path_contig^, ["chrA", "chrB"]))
    var lengths: List[UInt32] = [UInt32(NONE), 10, 5, 10, 5]
    return _cache(
        String(t"/tmp/dromos_test_refmap_{Int(twice)}.bin"), lengths^, "AAAAAAAAAACCCCCGGGGGGGGGGTTTTT", samples^
    )


def _project(m: RefMap, walk: List[Int], ts: Int, te: Int) -> Tuple[Int, Int, Int, Int]:
    """Why (0 placed), contig, start, end."""
    var anchors = List[Anchor]()
    var p = Placed(-1, -1, -1)
    var why = m.project(walk, ts, te, anchors, p)
    return (why, p.contig, p.start, p.end)


def _step(node: Int, reverse: Bool) -> Int:
    return node << 1 | Int(reverse)


def test_projection_as_pgcat() raises:
    var m = RefMap(_refmap_file(False), "CHM13")
    assert_equal(m.sample, "CHM13")
    assert_equal(m.length(2), 5)
    assert_equal(m.length(0), 0)
    # Forward, all three nodes, aligned from base 3 to 22 of the walk.
    var r = _project(m, [_step(1, False), _step(2, False), _step(3, False)], 3, 22)
    assert_true(r[0] == 0 and r[1] == 0 and r[2] == 103 and r[3] == 122)
    # Reverse walk <3 <2 <1: the same bases, read backwards.
    r = _project(m, [_step(3, True), _step(2, True), _step(1, True)], 3, 22)
    assert_true(r[0] == 0 and r[1] == 0 and r[2] == 103 and r[3] == 122)
    # An off-reference node between anchors is skipped.
    r = _project(m, [_step(1, False), _step(4, False), _step(3, False)], 0, 25)
    assert_true(r[0] == 0 and r[2] == 100 and r[3] == 125)
    # Steps going backwards along a forward walk, or the orientation flipping part way: no place.
    assert_equal(_project(m, [_step(3, False), _step(1, False)], 0, 20)[0], INCONSISTENT)
    assert_equal(_project(m, [_step(1, False), _step(3, True)], 0, 20)[0], INCONSISTENT)
    # A node missing from the graph, or only off-reference nodes: no place.
    assert_equal(_project(m, [_step(1, False), _step(9, False)], 0, 12)[0], MISSING_NODE)
    assert_equal(_project(m, [_step(4, False)], 0, 5)[0], NO_ANCHOR)


def test_ambiguous_placements() raises:
    var m = RefMap(_refmap_file(True), "CHM13")
    assert_equal(m.placement(3).path, 1)
    assert_equal(m.placement(3).start, 500)
    # Node 2 alone has two places; with node 1 before it, only path 0 fits.
    assert_equal(_project(m, [_step(2, False)], 0, 5)[0], AMBIGUOUS)
    var r = _project(m, [_step(1, False), _step(2, False)], 0, 15)
    assert_true(r[0] == 0 and r[1] == 0 and r[2] == 100 and r[3] == 115)


def test_bins_as_htslib() raises:
    assert_equal(N_BINS, 37449)
    assert_equal(reg2bin(0, 1), 4681)
    assert_equal(reg2bin(0, 1 << 14), 4681)
    assert_equal(reg2bin(0, (1 << 14) + 1), 585)
    assert_equal(reg2bin(0, 1 << 29), 0)
    assert_equal(bin_level(0), 0)
    assert_equal(bin_level(4681), 5)
    assert_equal(bin_parent(4681), 585)


def test_decimal_digits() raises:
    var bound = 1
    for d in range(1, 19):
        assert_equal(decimal_digits(bound), d)
        assert_equal(decimal_digits(bound * 10 - 1), d)
        bound *= 10
    assert_equal(decimal_digits(0), 1)
    assert_equal(decimal_digits(9223372036854775807), 19)


def _classify(m: RefMap, line: StaticString) -> Tuple[Int, Int, Int, Int]:
    """Kind, contig, start, end (as pgcat's classify gives them)."""
    var walk = List[Int]()
    var anchors = List[Anchor]()
    var placed = Placed(-1, -1, -1)
    var why = 0
    var aligned = 0
    var b = _text(line)
    var kind = classify(m, b.unsafe_ptr(), len(b), walk, anchors, placed, why, aligned)
    return (kind, placed.contig, placed.start, placed.end)


def test_lines_as_pgcat_classifies_them() raises:
    # pgcat's `lines` test, on its test map.
    var m = RefMap(_refmap_file(False), "CHM13")
    var r = _classify(m, "r1\t20\t0\t20\t+\t>1>2>3\t25\t2\t22\t20\t20\t60")
    assert_true(r[0] == PLACED and r[1] == 0 and r[2] == 102 and r[3] == 122)
    assert_equal(_classify(m, "r2\t20\t*\t*\t*\t*\t*\t*\t*\t0\t0\t255")[0], UNMAPPED)
    assert_equal(_classify(m, "r3\t20\t0\t20\t+\t>1>2")[0], MALFORMED)
    assert_equal(_classify(m, "r4\t20\t0\t20\t+\tchr1:5-9\t25\t2\t22\t20\t20\t60")[0], MALFORMED)
    assert_equal(_classify(m, "@HD\tVN:1")[0], SKIP)
    assert_equal(_classify(m, "")[0], SKIP)


def test_slim_keeps_columns_and_pgphase_tags() raises:
    var out = List[UInt8]()
    var line = _text("r\t1\t2\t3\t+\t>1\t5\t6\t7\t8\t9\t60\tcs:Z::5\trq:i:7\tbq:Z:II\tpu:Z:no_anchor\trs:Z:x\trqx\t")
    append_slim(out, line.unsafe_ptr(), len(line))
    assert_true(out == _text("r\t1\t2\t3\t+\t>1\t5\t6\t7\t8\t9\t60\trq:i:7\tpu:Z:no_anchor\trs:Z:x"))
    # Fewer than 12 columns: those there are; an empty last column is kept as one.
    out.clear()
    var short = _text("a\tb\t")
    append_slim(out, short.unsafe_ptr(), len(short))
    assert_true(out == _text("a\tb\t"))


def test_gam_told_from_gaf() raises:
    var lib = load_libdeflate()
    var gaf = _text("r1\t20\t0\t20\t+\t>1\t25\t2\t22\t20\t20\t60\n")
    assert_true(looks_like_gaf(lib, Span(gaf), PLAIN, "t"))
    assert_true(looks_like_gaf(lib, Span(_text("@HD\tVN:Z:1.0\n")), PLAIN, "t"))
    assert_true(looks_like_gaf(lib, Span(List[UInt8]()), PLAIN, "t"))
    assert_true(looks_like_gaf(lib, Span(gaf), GZIP, "t"))
    # A GAM group: count, the tag's length, GAM, then a message.
    var gam = List[UInt8]()
    _varint(gam, 2)
    _varint(gam, 3)
    gam.extend(Span(_text("GAM")))
    _varint(gam, 3)
    gam.extend(Span(_text("abc")))
    assert_true(not looks_like_gaf(lib, Span(gam), PLAIN, "t"))
    # Text without a tab on its first line is no GAF either.
    assert_true(not looks_like_gaf(lib, Span(_text("just words\n")), PLAIN, "t"))


# pgcatalog's surject.rs tests, ported: the alignments, chaining, linearizing and clipping dromos surject makes.


def _cigar(c: List[CigarOp]) -> String:
    var out = String()
    for o in c:
        out += String(o.length)
        out += chr(Int(o.op))
    return out^


def _linear(
    query: List[UInt8], blocks: List[Block], reference: List[UInt8], graph: Bool, windows: List[Window] = []
) -> String:
    """linearize's result as "start CIGAR", or "none"."""
    var start = 0
    var cigar = List[CigarOp]()
    if not linearize(Span(query), blocks, Span(reference), graph, windows, start, cigar):
        return "none"
    return String(t"{start} {_cigar(cigar)}")


def test_global_alignment() raises:
    assert_equal(_cigar(align_global(Span(_text("ACGT")), Span(_text("ACGT")))), "4M")
    assert_equal(_cigar(align_global(Span(_text("ACTT")), Span(_text("ACGT")))), "4M")
    # A deletion in a run goes to its left end.
    assert_equal(_cigar(align_global(Span(_text("GAAC")), Span(_text("GAAAC")))), "1M1D3M")
    assert_equal(_cigar(align_global(Span(_text("GAAAAC")), Span(_text("GAAC")))), "1M2I3M")
    assert_equal(_cigar(align_global(Span(_text("")), Span(_text("ACG")))), "3D")
    assert_equal(_cigar(align_global(Span(_text("ACG")), Span(_text("")))), "3I")


def test_long_gap_is_aligned_in_a_band() raises:
    # A 3 kb stretch with two mismatches and a 2-letter deletion.
    var x: UInt64 = 7
    var r = List[UInt8]()
    var letters = _text("ACGT")
    for _ in range(3000):
        x = x * 6364136223846793005 + 1442695040888963407
        r.append(letters[Int(x >> 33) % 4])
    var q = r.copy()
    q[500] = UInt8(ord("C")) if q[500] == UInt8(ord("A")) else UInt8(ord("A"))
    q[2500] = UInt8(ord("C")) if q[2500] == UInt8(ord("A")) else UInt8(ord("A"))
    _ = q.pop(1500)
    _ = q.pop(1500)
    var got = align_global(Span(q), Span(r))
    var m = 0
    var i = 0
    var d = 0
    for o in got:
        if o.op == M:
            m += o.length
        elif o.op == I:
            i += o.length
        else:
            d += o.length
    assert_true(m == 2998 and i == 0 and d == 2)
    # A read end of 3 kb against the reference beyond it: aligned, not clipped.
    var rr = r.copy()
    rr.extend(Span(_text("TTTTTTTTTT")))
    var ops = List[CigarOp]()
    assert_equal(align_end(Span(q), Span(rr), ops), 3000)
    m = 0
    for o in ops:
        if o.op == M:
            m += o.length
    assert_equal(m, 2998)


def test_chaining() raises:
    # A block far off the others' diagonal (a walk touching a reference node 50 kb away) is left out.
    var got = chain([Block(0, 1000, 500), Block(500, 51_500, 20), Block(520, 1520, 500)])
    assert_true(got == [Block(0, 1000, 500), Block(520, 1520, 500)])
    # A real 3 kb deletion between long blocks is kept.
    assert_equal(len(chain([Block(0, 1000, 5000), Block(5000, 9000, 5000)])), 2)


def test_linearized_blocks() raises:
    var reference = _text("AAAACCCCGGGGTTTTAAAACCCC")
    # Two blocks with a substitution bubble between (equal gaps): one match run.
    assert_equal(_linear(_text("AAAACCCCGGGG"), [Block(0, 0, 4), Block(6, 6, 6)], reference, False), "0 12M")
    # An insertion, then a deletion, between blocks; read ends that match nothing beyond the blocks clipped.
    var got = _linear(_text("TTAAAAXXCCCCGGNN"), [Block(2, 0, 4), Block(8, 4, 4), Block(12, 10, 2)], reference, False)
    assert_equal(got, "0 2S4M2I4M2D2M2S")
    # A block running backwards on the reference is dropped into the gap.
    got = _linear(_text("AAAACCCCGGGG"), [Block(0, 0, 4), Block(4, 2, 2), Block(6, 6, 6)], reference, False)
    assert_equal(got, "0 12M")
    # An indel bubble: aligned against the reference between the blocks.
    assert_equal(_linear(_text("AAAACCCGGGGTTTT"), [Block(0, 0, 4), Block(11, 12, 4)], reference, False), "0 4M1D11M")
    # Read ends that match the reference beyond the blocks are aligned, not clipped.
    assert_equal(_linear(_text("CCCGGGGTTTTAAA"), [Block(3, 8, 8)], reference, False), "5 14M")
    assert_equal(_linear(_text("ACGT"), [], reference, False), "none")


def test_window_aligned_as_an_allele() raises:
    # Reference AAAA GAGAGAGA CCCC; the read's graph alignment put a 2-letter deletion in one place and the same 2
    # letters back as an insertion further on (blocks either side of a window 4..12): aligned as the window's allele,
    # it is matches and one mismatch, not an insertion and a deletion.
    var r = _text("AAAAGAGAGAGACCCC")
    var blocks: List[Block] = [Block(0, 0, 6), Block(6, 8, 2), Block(10, 10, 6)]
    assert_equal(_linear(_text("AAAAGAGAGTGACCCC"), blocks, r, False, [Window(4, 12)]), "0 16M")


def test_divergent_ends_clipped() raises:
    # 40 bases mismatched every third, 100 matching, 40 mismatched every third: both ends clipped.
    var score = List[Int](length=180, fill=CLIP_MATCH)
    for i in range(180):
        if (i < 40 or i >= 140) and i % 3 == 0:
            score[i] = CLIP_MISMATCH
    var clips = clip_lengths(score)
    assert_true(clips[0] == 40 and clips[1] == 39)  # base 140 matches
    # One mismatch near an end is HiFi's; a long insertion (a structural variant) is free.
    var one = List[Int](length=180, fill=CLIP_MATCH)
    one[3] = CLIP_MISMATCH
    clips = clip_lengths(one)
    assert_true(clips[0] == 0 and clips[1] == 0)
    var sv = List[Int](length=180, fill=CLIP_MATCH)
    sv[19] += clip_gap(80)
    clips = clip_lengths(sv)
    assert_true(clips[0] == 0 and clips[1] == 0)
    # The CIGAR: clipped bases soft-clipped, the start past the clipped reference, a gap at the cut joins the clip.
    var c = clip_cigar(100, [CigarOp(M, 10), CigarOp(I, 2), CigarOp(M, 10), CigarOp(D, 3), CigarOp(M, 20)], 10, 5)
    assert_equal(c[0], 110)
    assert_equal(_cigar(c[1]), "12S10M3D15M5S")
    c = clip_cigar(100, [CigarOp(S, 4), CigarOp(M, 30)], 6, 0)
    assert_equal(c[0], 102)
    assert_equal(_cigar(c[1]), "6S28M")


def test_graph_placement() raises:
    var r = List[UInt8]()
    var letters = _text("ACGT")
    for i in range(400):
        r.append(letters[(i * 7 + i // 3) % 4])
    # Read = reference 0..100, 80 letters of another copy in place of reference 100..180, then 180..300.
    var read = List[UInt8](Span(r)[:100])
    for _ in range(80):
        read.append(UInt8(ord("T")))
    read.extend(Span(r)[180:300])
    var blocks: List[Block] = [Block(0, 0, 100), Block(180, 180, 120)]
    assert_equal(_linear(read, blocks, r, True), "0 100M80I80D120M")
    # Without it the equal-length stretch is aligned base for base.
    assert_equal(_linear(read, blocks, r, False), "0 300M")
    # A read end of 60 letters beyond the last block is clipped, even where it matches the reference.
    var head = List[UInt8](Span(r)[:260])
    assert_equal(_linear(head, [Block(0, 0, 200)], r, True), "0 200M60S")
    assert_equal(_linear(head, [Block(0, 0, 200)], r, False), "0 260M")


comptime GBZ = "tests/data/micb-kir3dl1.gbz"


def test_gbz_told_by_its_bytes() raises:
    assert_true(is_gbz(GBZ))
    assert_false(is_gbz("tests/test_dromos.mojo"))
    assert_false(is_gbz("tests/no_such_file"))


def test_reference_samples() raises:
    var found = reference_samples(GBZ)
    assert_equal(len(found), 2)
    assert_equal(found[0], "CHM13")
    assert_equal(found[1], "GRCh38")
    assert_equal(resolve_reference_sample(GBZ, String("GRCh38")), "GRCh38")
    var message = String()
    try:
        _ = resolve_reference_sample(GBZ, None)
    except e:
        message = String(e)
    assert_equal(
        message,
        String(t"{GBZ} declares several reference samples (CHM13, GRCh38); choose one with --reference-sample"),
    )


def test_reference_paths() raises:
    var by_contig = reference_paths(GBZ, "CHM13")
    assert_true(len(by_contig) > 0)
    for i in range(len(by_contig)):
        if i > 0:
            assert_true(by_contig[i - 1][0] < by_contig[i][0])
        for name in by_contig[i][1]:
            assert_true(name.startswith(String(t"CHM13#0#{by_contig[i][0]}")))
    assert_equal(len(reference_paths(GBZ, "nobody")), 0)


def test_fingerprint() raises:
    # FNV-1a of the file's first and last 4 MB (here the whole file, twice).
    assert_equal(fingerprint(GBZ), "size=73920 fnv1a=9736178bc136e029")


def test_paths_and_edges() raises:
    var g = GbzGraph(GBZ)
    var md = g.metadata()
    assert_equal(md.paths(), g.paths())
    var next = List[UInt64]()
    var prev = List[UInt64]()
    for id in range(g.paths()):
        var walk = g.path(id)
        var back = g.path(id, True)
        assert_equal(len(walk), len(back))
        for i in range(len(walk)):
            # The reverse walk visits the same nodes backwards, each flipped.
            assert_equal(back[len(walk) - 1 - i], walk[i] ^ 1)
            assert_true(g.has_node(Int(walk[i] >> 1)))
            if i + 1 < len(walk):
                assert_true(g.successors(Int(walk[i] >> 1), walk[i] & 1 == 1, next))
                assert_true(walk[i + 1] in next)
                assert_true(g.predecessors(Int(walk[i + 1] >> 1), walk[i + 1] & 1 == 1, prev))
                assert_true(walk[i] in prev)
    assert_false(g.successors(g.max_node() + 1, False, next))


def test_search_counts_haplotypes() raises:
    """A search state holds every visit of its node in that orientation, over both strands of every path; following
    it forward, every visit followed by each next node."""
    var g = GbzGraph(GBZ)
    var walks = List[List[UInt64]]()
    for id in range(g.paths()):
        walks.append(g.path(id))
        walks.append(g.path(id, True))
    var v = walks[0][3]
    var count = 0
    var after = Dict[Int, Int]()
    for w in walks:
        for i in range(len(w)):
            if w[i] == v:
                count += 1
                if i + 1 < len(w):
                    after[Int(w[i + 1])] = after.get(Int(w[i + 1]), 0) + 1
    var state = g.search_state(Int(v >> 1), v & 1 == 1)
    assert_true(Bool(state))
    assert_equal(state.value().count(), count)
    var nexts = List[SearchState]()
    assert_true(g.follow_forward(state.value(), nexts))
    assert_equal(len(nexts), len(after))
    for s in nexts:
        assert_equal(s.count(), after[Int(s.forward_node)])
    assert_false(Bool(g.search_state(g.max_node() + 1, False)))


def test_parts_read_as_the_graph() raises:
    """The parts loader gives the crate's GBZ's nodes, lengths, letters and paths."""
    var g = GbzGraph(GBZ)
    var p = GbzParts(GBZ)
    var count = g.max_node() + 1
    assert_equal(p.max_node(), g.max_node())
    var a = List[UInt64](length=count + 3, fill=0)
    var b = List[UInt64](length=count + 3, fill=0)
    g.lengths(0, count + 3, a.unsafe_ptr())
    p.lengths(0, count + 3, b.unsafe_ptr())
    var letters = List[UInt8]()
    for node in range(count + 3):
        assert_equal(a[node], b[node])
        if g.has_node(node):
            letters.extend(g.sequence(node))
    var got = List[UInt8](length=len(letters), fill=0)
    assert_equal(p.letters(0, count + 3, got.unsafe_ptr(), len(got)), len(letters))
    assert_true(got == letters)
    assert_equal(p.paths(), g.paths())
    for id in range(g.paths()):
        assert_true(p.path(id) == g.path(id))


def test_cache_from_the_graph() raises:
    """A graph cache made from the graph reads back as its nodes and each reference sample's paths, its letters the
    same at 4 bits and 8."""
    var g = GbzGraph(GBZ)
    make_cache(GBZ, "/tmp/dromos_test.dromos", ["CHM13", "GRCh38"], None, 4, 3)
    make_cache(GBZ, "/tmp/dromos_test8.dromos", ["GRCh38"], None, 8, 3)
    var t = NodeTable("/tmp/dromos_test.dromos")
    var t8 = NodeTable("/tmp/dromos_test8.dromos")
    assert_equal(t.count, g.max_node() + 1)
    assert_false(t.wide)
    assert_true(t8.wide)
    for node in range(t.count):
        var s = g.sequence(node)
        assert_equal(t.length(node), len(s))
        for i in range(len(s)):
            assert_equal(t.letter(node, False, i), s[i])
            assert_equal(t8.letter(node, False, i), s[i])
            assert_equal(t.letter(node, True, i), complement(s[len(s) - 1 - i]))
    var md = g.metadata()
    var names = md.path_names()
    for sample in ["CHM13", "GRCh38"]:
        var m = RefMap("/tmp/dromos_test.dromos", sample)
        for node in range(t.count):
            assert_equal(m.length(node), len(g.sequence(node)))
        var id_of = md.sample_id(sample).value()
        var path = 0
        for id in range(len(names)):
            if names[id].sample != id_of:
                continue
            var walk = g.path(id)
            var start = names[id].fragment
            var k = m.first_placement_of(path)
            for i in range(len(walk)):
                var pl = m.placement(k + i)
                assert_true(pl.path == path and pl.step == i and pl.start == start and pl.reverse == (walk[i] & 1 == 1))
                start += len(g.sequence(Int(walk[i] >> 1)))
            path += 1
        assert_equal(len(m.path_contig), path)
    var only = CacheFile("/tmp/dromos_test8.dromos")
    assert_equal(len(only.samples), 1)
    var raised = False
    try:
        _ = RefMap("/tmp/dromos_test8.dromos", "CHM13")
    except:
        raised = True
    assert_true(raised)


def test_letters_packed() raises:
    var letters = _text("ACGTNRYKMSWBDHVUACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACG")
    var packed = List[UInt8](length=(len(letters) + 1) // 2, fill=0)
    assert_true(pack_letters(letters.unsafe_ptr(), len(letters), packed.unsafe_ptr()))
    assert_equal(packed[0], 0x10)
    assert_equal(packed[len(packed) - 1], 0x02)
    var t = NodeTable(_nodes_file(4))
    var t8 = NodeTable(_nodes_file(8))
    for node in range(4):
        for i in range(t.length(node)):
            assert_equal(t.letter(node, False, i), t8.letter(node, False, i))
            assert_equal(t.letter(node, True, i), t8.letter(node, True, i))
    letters[70] = UInt8(ord("x"))
    assert_false(pack_letters(letters.unsafe_ptr(), len(letters), packed.unsafe_ptr()))
    letters[70] = UInt8(ord("A"))
    letters[len(letters) - 1] = UInt8(ord("a"))
    assert_false(pack_letters(letters.unsafe_ptr(), len(letters), packed.unsafe_ptr()))


def test_cache_examined() raises:
    """A cache is told from a missing file, another file, another version's cache and a damaged one."""
    var h = Header()
    assert_equal(examine("/tmp/dromos_no_such_cache", h), MISSING)
    assert_equal(examine(GBZ, h), NOT_A_CACHE)
    var path = _nodes_file()
    assert_equal(examine(path, h), CURRENT)
    assert_true(h.graph == GraphId(1, 2, 3))
    assert_equal(h.name, "test.gbz")
    var bytes: List[UInt8]
    with open(path, "r") as f:
        bytes = f.read_bytes()
    var damaged = bytes.copy()
    damaged[100] ^= 1
    var cut = bytes.copy()
    cut.resize(len(bytes) - 10, 0)
    var other = bytes.copy()
    other[8] = 99
    var cases: List[Tuple[Int, List[UInt8]]] = [(DAMAGED, damaged^), (DAMAGED, cut^), (OTHER_VERSION, other^)]
    for c in cases:
        with open("/tmp/dromos_test_bad.dromos", "w") as f:
            f.write_bytes(Span(c[1]))
        assert_equal(examine("/tmp/dromos_test_bad.dromos", h), c[0])
        var raised = False
        try:
            _ = NodeTable("/tmp/dromos_test_bad.dromos")
        except:
            raised = True
        assert_true(raised)


def test_folders_and_names() raises:
    assert_equal(folder_of("a/b/g.gbz"), "a/b")
    assert_equal(folder_of("/g.gbz"), "/")
    assert_equal(folder_of("g.gbz"), ".")
    assert_equal(name_of("a/b/g.gbz"), "g.gbz")
    assert_equal(name_of("g.gbz"), "g.gbz")


def test_parts_search_as_the_graph() raises:
    """The parts' edges and bidirectional search, many states at a time, give what the crate's GBZ gives one at a
    time."""
    var g = GbzGraph(GBZ)
    var p = GbzParts(GBZ)
    p.find_nodes(2)
    var visits = List[UInt64]()
    for node in range(g.max_node() + 2):
        visits.append(UInt64(node << 1))
        visits.append(UInt64(node << 1 | 1))
    for predecessors in [False, True]:
        var out = List[UInt64]()
        var counts = List[Int64]()
        p.edges(visits.unsafe_ptr(), len(visits), predecessors, out, counts)
        var at = 0
        var one = List[UInt64]()
        for i in range(len(visits)):
            var v = visits[i]
            var has: Bool
            if predecessors:
                has = g.predecessors(Int(v >> 1), v & 1 == 1, one)
            else:
                has = g.successors(Int(v >> 1), v & 1 == 1, one)
            if not has:
                assert_equal(counts[i], -1)
                continue
            assert_equal(Int(counts[i]), len(one))
            for k in range(len(one)):
                assert_equal(out[at + k], one[k])
            at += len(one)
    var states = List[SearchState]()
    for v in visits:
        var s = g.search_state(Int(v >> 1), v & 1 == 1)
        var t = p.search_state(Int(v >> 1), v & 1 == 1)
        assert_equal(Bool(s), Bool(t))
        if s:
            assert_true(s.value() == t.value())
            states.append(s.value())
    var out = List[SearchState]()
    var counts = List[Int64]()
    p.follow_forward(states.unsafe_ptr(), len(states), out, counts)
    var at = 0
    var one = List[SearchState]()
    for i in range(len(states)):
        assert_true(g.follow_forward(states[i], one))
        assert_equal(Int(counts[i]), len(one))
        for k in range(len(one)):
            assert_true(out[at + k] == one[k])
        at += len(one)
    assert_equal(at, len(out))


def _snarl_message(start: UInt64, end: UInt64, parent: Optional[Tuple[UInt64, UInt64]]) -> List[UInt8]:
    def visit(v: UInt64) -> List[UInt8]:
        var m = List[UInt8]()
        _varint(m, 1 << 3)
        _varint(m, Int(v >> 1))
        if v & 1 == 1:
            _varint(m, 3 << 3)
            _varint(m, 1)
        return m^

    var m = List[UInt8]()
    if parent:
        var pm = List[UInt8]()
        _bytes_field(pm, 2, visit(parent.value()[0]))
        _bytes_field(pm, 3, visit(parent.value()[1]))
        _bytes_field(m, 4, pm)
    _bytes_field(m, 2, visit(start))
    _varint(m, 6 << 3)  # a varint field read past
    _varint(m, 300)
    _bytes_field(m, 3, visit(end))
    return m^


def _group(mut out: List[UInt8], tag: String, messages: List[List[UInt8]]):
    _varint(out, len(messages) + 1)
    _varint(out, tag.byte_length())
    out.extend(tag.as_bytes())
    for m in messages:
        _varint(out, len(m))
        out.extend(Span(m))


def test_snarls_stream() raises:
    """The groups of vg snarls: SNARL messages (in any case of the tag) read for their visits and parents, other groups
    skipped; a damaged stream fails."""
    var stream = List[UInt8]()
    _group(
        stream,
        "SNARL",
        [
            _snarl_message(10 << 1, 12 << 1, None),
            _snarl_message(11 << 1 | 1, 7 << 1, (UInt64(14 << 1 | 1), UInt64(9 << 1))),
        ],
    )
    _group(stream, "OTHER", [List[UInt8](length=5, fill=1)])
    _group(stream, "snarl", [_snarl_message(3 << 1, 4 << 1 | 1, (UInt64(1 << 1), UInt64(2 << 1)))])
    var s = parse_snarls(Span(stream), 2)
    assert_equal(len(s), 3)
    assert_true(s[0].start == 20 and s[0].end == 24 and not s[0].has_parent())
    assert_true(s[1].start == 23 and s[1].end == 14 and s[1].parent_start == 29 and s[1].parent_end == 18)
    assert_true(s[2].start == 6 and s[2].end == 9 and s[2].parent_start == 2 and s[2].parent_end == 4)
    var cut = List[UInt8](Span(stream)[: len(stream) - 3])
    var failed = False
    try:
        _ = parse_snarls(Span(cut), 2)
    except e:
        failed = String(e).startswith("malformed snarls file")
    assert_true(failed)
    var no_end = List[UInt8]()
    var m = List[UInt8]()
    _bytes_field(m, 2, List[UInt8]([UInt8(8), UInt8(5)]))
    _group(no_end, "SNARL", [m^])
    failed = False
    try:
        _ = parse_snarls(Span(no_end), 1)
    except e:
        failed = String(e) == "malformed snarls file: snarl without end"
    assert_true(failed)


def test_snarl_names_and_keys() raises:
    assert_equal(snarl_name(12 << 1, 34 << 1), ">12>34")
    assert_equal(snarl_name(12 << 1 | 1, 3 << 1 | 1), "<12<3")
    # A snarl read backwards has the same key: (start, end) and (end flipped, start flipped).
    var a = canonical(5 << 1, 9 << 1)
    var b = canonical(9 << 1 | 1, 5 << 1 | 1)
    assert_true(a[0] == b[0] and a[1] == b[1])
    assert_true(a[0] == 10 and a[1] == 18)


def test_regions() raises:
    var s = parse_region("chr20", "CHM13")
    assert_true(s.contig == "chr20" and not s.region)
    s = parse_region("CHM13#0#chr20:10,000,000-12,000,000", "CHM13")
    assert_true(s.contig == "chr20" and s.region.value()[0] == 10_000_000 and s.region.value()[1] == 12_000_000)
    for bad in ["chr20:0-5", "chr20:9-5", "chr20:5", "GRCh38#0#chr20", "chr20:a-5"]:
        var failed = False
        try:
            _ = parse_region(bad, "CHM13")
        except:
            failed = True
        assert_true(failed)


def test_contigs_in_natural_order() raises:
    var v: List[String] = ["chrX", "chr10", "chr2", "chr1", "chrM", "chr22", "chrY"]
    sort(Span(v), natural_less)
    var want: List[String] = ["chr1", "chr2", "chr10", "chr22", "chrM", "chrX", "chrY"]
    assert_true(v == want)


def test_same_contig() raises:
    assert_true(same_contig("chr20", "CHM13#0#chr20"))
    assert_true(same_contig("CHM13#0#chr20", "chr20"))
    assert_true(same_contig("CHM13#0#chr20", "CHM13#0#chr20"))
    assert_false(same_contig("GRCh38#0#chr20", "CHM13#0#chr20"))
    assert_false(same_contig("chr2", "CHM13#0#chr20"))


def test_region_queries_return_the_overlapping_records() raises:
    """Records across several windows, one with a REF reaching far past its window: every region query returns
    exactly the overlapping records (pgsites.rs's test)."""
    var path = String("/tmp/dromos_test_region.pgsites")
    var header = String("##fileformat=VCFv4.2\n##contig=<ID=CHM13#0#chrT,length=1000000>\n#CHROM\tPOS\tID\tREF\tALT\n")
    var positions = List[Int]()
    var lengths = List[Int]()
    var w = Writer(path, header, 3, 2)
    for k in range(400):
        var pos = 1 + k * 1_337
        var length = 250_000 if k == 7 else 2
        positions.append(pos)
        lengths.append(length)
        var line = String(t"CHM13#0#chrT\t{pos}\t>{k}>{k + 1}\t")
        line += "A" * length if k == 7 else "AC"
        line += "\tG"
        var place = w.place(pos, pos + length - 1, length, 3, line.byte_length())
        var at = w.at(place[0], place[1])
        for i in range(line.byte_length()):
            at[unsafe_offset=i] = line.as_bytes()[i]
        if k % 50 == 0:
            compress_sealed(w)
            w.write_done()
    w.seal()
    compress_sealed(w)
    w.write_done()
    assert_true(w.finish() > 1)
    var c = Catalog(path)
    assert_equal(c.chrom, "CHM13#0#chrT")
    assert_equal(c.header, header)
    var queries: List[Tuple[Int, Int]] = [
        (1, 1),
        (2, 2),
        (3, 3),
        (100_000, 100_001),
        (200_000, 260_000),
        (1, Int.MAX),
        (600_000, 700_000),
    ]
    for q in queries:
        var got = List[UInt8]()
        records(c, q[0], q[1], 2, -1, got)
        var got_pos = List[Int]()
        var text = String(unsafe_from_utf8=got)
        for line in text.split("\n"):
            if line.byte_length() > 0:
                got_pos.append(Int(line.split("\t")[1]))
        var want = List[Int]()
        for k in range(len(positions)):
            if positions[k] <= q[1] and positions[k] + lengths[k] > q[0]:
                want.append(positions[k])
        assert_true(got_pos == want)


def _read_all(path: String) raises -> List[UInt8]:
    with open(path, "r") as f:
        return f.read_bytes()


def test_catalogs_of_the_small_graph() raises:
    """The small graph's catalogs, built from its snarls, are pgcatalog's byte for byte, on one thread and on
    several; a region holds the records a query of the whole catalog returns."""
    var data = String("tests/data/micb-kir3dl1.catalog")
    var graph = Graph(GBZ, 2)
    var md = graph.parts.metadata()
    for threads in [1, 3]:
        for contig in ["chr6", "chr19"]:
            var r = ReferencePath(graph, md, "CHM13", contig, threads)
            var out = String(t"/tmp/dromos_test_{contig}")
            _ = build_catalog(graph, r, data + "/snarls/" + contig + ".snarls", None, out, 6, threads)
            assert_true(_read_all(out + ".pgsites") == _read_all(data + "/" + contig + ".pgsites"))
    var r = ReferencePath(graph, md, "CHM13", "chr6", 2)
    var whole = Catalog(data + "/chr6.pgsites")
    var beg = whole.blocks[0] + 2000
    var end = beg + 5000
    _ = build_catalog(graph, r, data + "/snarls/chr6.snarls", (beg, end), "/tmp/dromos_test_region", 6, 2)
    var want = List[UInt8]()
    records(whole, beg, end, 2, -1, want)
    var got = List[UInt8]()
    records(Catalog("/tmp/dromos_test_region.pgsites"), 1, Int.MAX, 2, -1, got)
    assert_true(len(want) > 0 and got == want)


def test_reference_traversals() raises:
    """Every snarl of the small graph whose boundary nodes the reference visits in order has a traversal from the
    start's step to the end's."""
    var graph = Graph(GBZ, 2)
    var md = graph.parts.metadata()
    var r = ReferencePath(graph, md, "CHM13", "chr6", 2)
    var snarls = read_snarls("tests/data/micb-kir3dl1.catalog/snarls/chr6.snarls", 2)
    var found = 0
    var block = List[UInt64]()
    for sn in snarls:
        var t = Traversal(-1, 0, False)
        if traversals(r, sn, block, t) == 1:
            found += 1
            var first = r.steps[t.s0] ^ UInt64(Int(t.backward))
            var last = r.steps[t.e0] ^ UInt64(Int(t.backward))
            assert_true(first >> 1 == sn.start >> 1 and last == sn.end)
    assert_true(found > 100)
    assert_true(
        vcf_header("X", 5).endswith("##contig=<ID=X,length=5>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    )


def main() raises:
    TestSuite.discover_tests[__functions_in_module()]().run()
