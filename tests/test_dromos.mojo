"""Tests: protobuf varints, GAM group framing, GAF lines against hand-worked ones, placement on the reference as
pgcat places alignments, tabix bins as htslib numbers them, surject's alignments as pgcatalog makes them, a GBZ read
through the gbz crate (tests/data/micb-kir3dl1.gbz, two reference samples), snarls streams, site catalogs (their
columns, node tables, index search and checks; tests/data/micb-kir3dl1.catalog: the small graph's vg snarls and
pgcatalog view's text of pgcatalog's catalogs of it), and snarls from a distance index (tests/data/micb-kir3dl1.dist,
vg 1.73's)."""

from std.memory import OwnedPointer, bitcast
from std.os import getenv
from std.random import random_ui64, seed
from std.testing import TestSuite, assert_equal, assert_false, assert_true

from dromos.gam import Alignment, Framer, frame, gaf_line, slim_line, varint
from dromos.nodes import NodeTable, comp, complement
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
from dromos.build import check_dist, natural_less, parse_region, reference_md5
from dromos.catalog import Encoder, Graph, Placed as SitePlaced, Record, ReferencePath, Traversal, build_catalog
from dromos.catalog import encode_block, traversals
from dromos.dist import DistIndex
from dromos.pool import Pool
from dromos.sites import (
    COLUMNS,
    EXC,
    NODEID,
    NODELEN,
    LONG_BYTES,
    SEQ,
    WALK,
    ContigEntry,
    Cursor,
    EncodedBlock,
    LetterPacker,
    Md5,
    NodeTable as SiteNodes,
    Sites,
    SitesWriter,
    hex_string,
    put_oriented,
    put_varint,
    records,
    same_contig,
    unzigzag,
    vcf_header,
    zigzag,
)
from dromos.snarls import NO_PARENT, Snarl, canonical, parse_snarls, read_snarls, snarl_name
from dromos.gaf_text import MALFORMED, PLACED, SKIP, UNMAPPED, append_slim, classify
from dromos.output import decimal_digits
from dromos.bgzf import inflate_blocks, load_libdeflate, scan_blocks
from dromos.stream import BGZF, GZIP, PLAIN, looks_like_gaf, mut_untracked
from dromos.refmap import AMBIGUOUS, INCONSISTENT, MISSING_NODE, NO_ANCHOR, Anchor, Placed, RefMap
from dromos.tabix import N_BINS, bin_level, bin_parent, reg2bin
from dromos.linearize import (
    ASM_COSTS,
    CLIP_MATCH,
    CLIP_MISMATCH,
    D,
    I,
    M,
    READ_COSTS,
    S,
    Banded,
    Block,
    CigarOp,
    Costs,
    Window,
    align_end,
    align_global,
    chain,
    clip_cigar,
    clip_gap,
    clip_lengths,
    half_band,
    linearize,
    upper,
)
from dromos.bam import BamIndex, RecordWriter, bam_header, csi_levels, needs_csi, nt16, record_end
from dromos.bam_sort import BamStaging, write_bam



def _tmp(name: String) -> String:
    """A scratch file's path: in $DROMOS_TEST_DIR, else /tmp."""
    var dir = getenv("DROMOS_TEST_DIR")
    return (dir if dir else String("/tmp")) + "/" + name

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
        String(_tmp(String(t"dromos_test_nodes_{width}.bin"))), lengths^, "ACGTACGTGGGTTTTT", List[ReferenceSample](), width
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
        String(_tmp(String(t"dromos_test_refmap_{Int(twice)}.bin"))), lengths^, "AAAAAAAAAACCCCCGGGGGGGGGGTTTTT", samples^
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
    make_cache(GBZ, _tmp("dromos_test.dromos"), ["CHM13", "GRCh38"], None, 4, 3)
    make_cache(GBZ, _tmp("dromos_test8.dromos"), ["GRCh38"], None, 8, 3)
    var t = NodeTable(_tmp("dromos_test.dromos"))
    var t8 = NodeTable(_tmp("dromos_test8.dromos"))
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
        # Stretches decoded two letters a lookup: every start and length parity, both strands.
        var start = t.start(node)
        var start8 = t8.start(node)
        for off in range(min(3, len(s))):
            for take in [len(s) - off, (len(s) - off) // 2, 1]:
                if take <= 0:
                    continue
                var got = List[UInt8](length=take, fill=0)
                var got8 = List[UInt8](length=take, fill=0)
                t.copy_forward(got.unsafe_ptr(), start + off, take)
                t8.copy_forward(got8.unsafe_ptr(), start8 + off, take)
                for k in range(take):
                    assert_equal(got[k], upper(t.at(start + off + k)))
                    assert_equal(got8[k], got[k])
                var last = start + len(s) - 1 - off
                t.copy_reverse(got.unsafe_ptr(), last, take)
                t8.copy_reverse(got8.unsafe_ptr(), start8 + len(s) - 1 - off, take)
                for k in range(take):
                    assert_equal(got[k], comp(t.at(last - k)))
                    assert_equal(got8[k], got[k])
    var md = g.metadata()
    var names = md.path_names()
    for sample in ["CHM13", "GRCh38"]:
        var m = RefMap(_tmp("dromos_test.dromos"), sample)
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
    var only = CacheFile(_tmp("dromos_test8.dromos"))
    assert_equal(len(only.samples), 1)
    var raised = False
    try:
        _ = RefMap(_tmp("dromos_test8.dromos"), "CHM13")
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
    assert_equal(examine(_tmp("dromos_no_such_cache"), h), MISSING)
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
        with open(_tmp("dromos_test_bad.dromos"), "w") as f:
            f.write_bytes(Span(c[1]))
        assert_equal(examine(_tmp("dromos_test_bad.dromos"), h), c[0])
        var raised = False
        try:
            _ = NodeTable(_tmp("dromos_test_bad.dromos"))
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


def _read_all(path: String) raises -> List[UInt8]:
    with open(path, "r") as f:
        return f.read_bytes()


def test_site_varints() raises:
    """Column integers: LEB128 varints and zigzag differences read back, and a column cut short is an error."""
    var values: List[UInt64] = [0, 1, 127, 128, 300, 16383, 16384, UInt64(1) << 42, UInt64.MAX]
    var col = List[UInt8]()
    for v in values:
        put_varint(col, v)
    var diffs: List[Int] = [0, 1, -1, 63, -64, 1 << 40, -(1 << 40), Int.MAX, Int.MIN]
    for d in diffs:
        put_varint(col, zigzag(d))
    var c = Cursor(col)
    for v in values:
        assert_equal(c.varint(), v)
    for d in diffs:
        assert_equal(unzigzag(c.varint()), d)
    assert_equal(c.at, c.end)
    assert_equal(zigzag(-1), 1)
    assert_equal(zigzag(1), 2)
    var cut = List[UInt8]()
    put_varint(cut, 1 << 20)
    cut.resize(len(cut) - 1, 0)
    var failed = False
    try:
        var k = Cursor(cut)
        _ = k.varint()
    except:
        failed = True
    assert_true(failed)


def _site_node_cols(ids: List[UInt64], letters: List[String]) -> List[List[UInt8]]:
    var cols = List[List[UInt8]](length=COLUMNS, fill=List[UInt8]())
    var packer = LetterPacker()
    var prev = UInt64(0)
    var seq = List[UInt8]()
    var exc = List[UInt8]()
    for k in range(len(ids)):
        put_varint(cols[NODEID], ids[k] - prev)
        prev = ids[k]
        put_varint(cols[NODELEN], UInt64(letters[k].byte_length()))
        packer.add(letters[k].as_bytes(), seq, exc)
    packer.finish(seq)
    cols[SEQ] = seq^
    cols[EXC] = exc^
    return cols^


def test_site_node_table() raises:
    """A block's node table: ids as differences, lengths, letters two bits each with the others (N, IUPAC codes,
    anything) as exceptions, lower case stored as capitals; read back exactly."""
    var ids: List[UInt64] = [3, 4, 9, 1000, 1 << 33]
    var letters: List[String] = ["ACGT", "", "acgtNNa", "RYKM*-ACG", "T" * 37]
    var cols = _site_node_cols(ids, letters)
    var t = SiteNodes()
    t.read(cols, True)
    assert_equal(t.count(), len(ids))
    for k in range(len(ids)):
        assert_equal(t.ids[k], ids[k])
        var got = String(unsafe_from_utf8=Span(t.letters)[t.starts[k] : t.starts[k + 1]])
        assert_equal(got, letters[k].upper())
    # Six exceptions in the third node and six in the fourth.
    var exc = Cursor(cols[EXC])
    var n = 0
    while exc.at < exc.end:
        _ = exc.varint()
        _ = exc.byte()
        n += 1
    assert_equal(n, 2 + 6)
    # A SEQ column of the wrong size is refused.
    cols[SEQ].append(0)
    var failed = False
    try:
        t.read(cols, True)
    except:
        failed = True
    assert_true(failed)


def test_md5() raises:
    var h = Md5()
    assert_equal(hex_string(Span(h.digest())), "d41d8cd98f00b204e9800998ecf8427e")
    var fox = String("The quick brown fox jumps over the lazy dog")
    var whole = Md5()
    whole.update(fox.as_bytes())
    assert_equal(hex_string(Span(whole.digest())), "9e107d9d372bb6826bd81d3542a419d6")
    # In pieces across 64-byte blocks.
    var long = String("ACGT" * 1000)
    var once = Md5()
    once.update(long.as_bytes())
    var pieces = Md5()
    var b = long.as_bytes()
    var at = 0
    var step = 1
    while at < len(b):
        var n = min(step, len(b) - at)
        pieces.update(b[at : at + n])
        at += n
        step = step * 3 + 1
    assert_true(once.digest() == pieces.digest())


def _synthetic_catalog(path: String, firsts: List[Int], ends: List[Int], split: Int) raises:
    """A catalog of empty blocks (no records), two contigs: the first `split` blocks the first contig's."""
    var props: List[Tuple[String, String]] = [("graph.fingerprint", "size=1 fnv1a=0"), ("x", "y")]
    var entries = List[ContigEntry]()
    entries.append(ContigEntry("chrA", "S#0#chrA", 0, List[UInt8](), List[Tuple[Int, Int]](), 0, 0, 0, 0))
    entries.append(ContigEntry("chrB", "S#0#chrB", 0, List[UInt8](), [(5, 9)], 0, 0, 0, 0))
    var w = SitesWriter(path, props^, entries^)
    var md5 = List[UInt8](length=16, fill=7)
    for k in [1, 0]:
        var blocks = List[EncodedBlock]()
        for b in range(len(firsts)):
            if (b < split) == (k == 0):
                var e = EncodedBlock()
                e.first_pos = firsts[b]
                e.last_pos = firsts[b]
                e.max_end = ends[b]
                blocks.append(e^)
        w.complete(k, blocks^, md5, 0, 1000 * (k + 1))
    _ = w.finish()


def test_site_index_search() raises:
    """A region's blocks found by binary search (first_pos, and the running largest max_end, `reach`) are exactly the
    blocks a scan would read, a block reaching far past its window (a giant record's) included; the header reads
    back."""
    var path = String(_tmp("dromos_test_index.sites"))
    var firsts = List[Int]()
    var ends = List[Int]()
    seed(7)
    var pos = 1
    for b in range(300):
        pos += Int(random_ui64(0, 3000))
        firsts.append(pos)
        ends.append(pos + (250_000 if b == 40 else Int(random_ui64(0, 5000))))
    _synthetic_catalog(path, firsts, ends, 120)
    var s = Sites(path)
    s.load_all()
    assert_equal(len(s.contigs), 2)
    assert_equal(s.contigs[0].name, "chrA")
    assert_equal(s.contigs[1].regions[0][1], 9)
    assert_equal(s.contigs[1].length, 2000)
    assert_equal(s.prop("x").value(), "y")
    assert_equal(s.contig("S#0#chrB").value(), 1)
    assert_false(Bool(s.contig("T#0#chrB")))
    for k in range(2):
        var lo = 0 if k == 0 else 120
        var hi = 120 if k == 0 else 300
        for _ in range(500):
            var a = Int(random_ui64(1, UInt64(pos)))
            var e = a + Int(random_ui64(0, 20_000))
            var got = List[Int]()
            var span = s.overlapping(k, a, e)
            for b in range(span[0], span[1]):
                if s.max_end(k, b) >= a:
                    got.append(lo + b)
            var want = List[Int]()
            for b in range(lo, hi):
                if firsts[b] <= e and ends[b] >= a:
                    want.append(b)
            assert_true(got == want)


def _expect_refused(path: String, words: String) raises:
    var failed = False
    try:
        var s = Sites(path)
        s.load_all()
    except e:
        failed = String(e).find(words) >= 0
        if not failed:
            print(e)
    assert_true(failed)


def test_site_catalog_refused_when_damaged() raises:
    """A truncated catalog, a damaged one, another version's, a .pgsites, an empty file: each refused, saying so."""
    var path = String(_tmp("dromos_test_bad.sites"))
    var firsts: List[Int] = [1, 10, 20]
    var ends: List[Int] = [5, 15, 25]
    _synthetic_catalog(path, firsts, ends, 2)
    var good = _read_all(path)
    _ = Sites(path)

    def put(path: String, data: List[UInt8]) raises:
        with open(path, "w") as f:
            f.write_bytes(data)

    var bad = path + ".bad"
    put(bad, List[UInt8](Span(good)[: len(good) - 10]))
    _expect_refused(bad, "truncated")
    var flipped = good.copy()
    flipped[len(flipped) - 3] ^= 0xFF
    put(bad, flipped)
    _expect_refused(bad, "damaged")
    var version = good.copy()
    version[8] = 2
    put(bad, version)
    _expect_refused(bad, "version 2")
    var pre = good.copy()
    pre[30] ^= 1
    put(bad, pre)
    _expect_refused(bad, "damaged")
    put(bad, List[UInt8]("PGSITES2xxxxxxxxxxxxxxxx".as_bytes()))
    _expect_refused(bad, ".pgsites")
    put(bad, List[UInt8]())
    _expect_refused(bad, "not a dromos site catalog")


def _view_of(path: String, contig: String, beg: Int = 1, end: Int = Int.MAX) raises -> String:
    var s = Sites(path)
    var k = s.contig(contig).value()
    s.load(k)
    var c: List[Tuple[String, Int]] = [(s.contigs[k].chrom, s.contigs[k].length)]
    var out = List[UInt8](vcf_header(c).as_bytes())
    records(s, k, beg, end, 2, -1, out)
    return String(unsafe_from_utf8=out)


def _small_catalog(path: String, contigs: List[String], regions: List[Tuple[Int, Int]], threads: Int) raises:
    """The small graph's catalog of `contigs` (each limited to `regions`), from its committed vg snarls."""
    var data = String("tests/data/micb-kir3dl1.catalog")
    var graph = Graph(GBZ, 2)
    var md = graph.parts.metadata()
    var entries = List[ContigEntry]()
    for c in contigs:
        entries.append(ContigEntry(c, "CHM13#0#" + c, 0, List[UInt8](), regions.copy(), 0, 0, 0, 0))
    var props: List[Tuple[String, String]] = [("graph.fingerprint", fingerprint(GBZ))]
    var w = SitesWriter(path, props^, entries^)
    var pool = OwnedPointer(Pool(threads))
    pool[].start()
    for k in range(len(contigs)):
        var r = ReferencePath(graph, md, "CHM13", contigs[k], threads)
        var length = r.length
        var md5 = reference_md5(graph, r)
        var snarls = read_snarls(data + "/snarls/" + contigs[k] + ".snarls", 2)
        var done = build_catalog(graph, r, snarls^, regions, 9, threads, mut_untracked(pool.ptr()))
        var blocks = List[EncodedBlock]()
        swap(blocks, done[1])
        w.complete(k, blocks^, md5, done[0].records, length)
    pool[].close()
    _ = w.finish()


def test_catalogs_of_the_small_graph() raises:
    """The small graph's catalog, built from its snarls, prints pgcatalog view's text of pgcatalog's catalogs (as
    committed) on one thread and on several; a region holds the records a query of the whole catalog returns."""
    var data = String("tests/data/micb-kir3dl1.catalog")
    var contigs: List[String] = ["chr6", "chr19"]
    var sizes = List[Int]()
    for threads in [1, 3]:
        var path = String(_tmp(String(t"dromos_test_small_{threads}.sites")))
        _small_catalog(path, contigs, List[Tuple[Int, Int]](), threads)
        sizes.append(len(_read_all(path)))
        for c in contigs:
            assert_equal(_view_of(path, c), String(unsafe_from_utf8=_read_all(data + "/" + c + ".vcf")))
    assert_true(_read_all(_tmp("dromos_test_small_1.sites")) == _read_all(_tmp("dromos_test_small_3.sites")))
    # A region's catalog: the records the whole one returns for it.
    var whole = Sites(_tmp("dromos_test_small_1.sites"))
    whole.load(0)
    var beg = whole.first_pos(0, 0) + 2000
    var end = beg + 5000
    var regions: List[Tuple[Int, Int]] = [(beg, end)]
    var chr6: List[String] = ["chr6"]
    _small_catalog(_tmp("dromos_test_region.sites"), chr6, regions, 2)
    var want = _view_of(_tmp("dromos_test_small_1.sites"), "chr6", beg, end)
    assert_true(want.count("\n") > 10)
    assert_equal(_view_of(_tmp("dromos_test_region.sites"), "chr6"), want)
    var r = Sites(_tmp("dromos_test_region.sites"))
    assert_equal(r.contigs[0].regions[0][0], beg)
    # Its MD5 is the reference sequence's.
    assert_equal(hex_string(Span(r.contigs[0].md5)), hex_string(Span(whole.contigs[0].md5)))


def test_giant_record() raises:
    """A record of a few hundred thousand walks (its WALK column long enough for long-distance matching and zstd's
    own threads) encoded, written and read back as the VCF line its walks spell."""
    var graph = Graph(GBZ, 2)
    # A walk over a few nodes of the graph, varied per allele.
    var nodes: List[UInt64] = [5, 6, 7, 8, 9, 10, 11, 12]
    var r = Record()
    r.present = True
    r.pos = 1000
    r.start = UInt64(2 << 1)
    r.end = UInt64(40 << 1)
    r.pad = -1
    var walks = 200_000
    for k in range(walks):
        r.starts.append(len(r.visits))
        r.visits.append(r.start)
        for j in range(30):
            var v = nodes[(k * 7 + j * (k % 5 + 1)) % len(nodes)] << 1 | UInt64((k + j) & 1)
            r.visits.append(v)
        r.visits.append(r.end)
        r.ends.append(len(r.visits))
    for k in range(walks):
        var n = 0
        for i in range(r.starts[k] + 1, r.ends[k] - 1):
            n += graph.length(Int(r.visits[i] >> 1))
        r.lens.append(n)
    var items = List[SitePlaced]()
    items.append(SitePlaced(r^, 3, UInt64(1 << 1), UInt64(99 << 1 | 1), True))
    var e = Encoder()
    var block = EncodedBlock()
    encode_block(graph, Span(items), e, 3, 2, block)
    assert_true(block.sizes[2 * WALK + 1] > LONG_BYTES)
    block.first_pos = 1000
    block.last_pos = 1000
    block.max_end = 1000
    block.records = 1
    var path = String(_tmp("dromos_test_giant.sites"))
    var entries: List[ContigEntry] = [
        ContigEntry("chrG", "S#0#chrG", 0, List[UInt8](), List[Tuple[Int, Int]](), 0, 0, 0, 0)
    ]
    var w = SitesWriter(path, List[Tuple[String, String]](), entries^)
    var blocks = List[EncodedBlock]()
    blocks.append(block^)
    w.complete(0, blocks^, List[UInt8](length=16, fill=0), 1, 5000)
    _ = w.finish()
    # The line, made here from the record's walks and the graph's letters.
    ref g = items[0].record
    var want = String("S#0#chrG\t1000\t>2>40\t")
    for k in range(walks):
        if k > 1:
            want += ","
        var allele = List[UInt8]()
        for i in range(g.starts[k] + 1, g.ends[k] - 1):
            var v = g.visits[i]
            var s = graph.letters(Int(v >> 1))
            var at = len(allele)
            allele.resize(at + len(s), 0)
            put_oriented(allele.unsafe_ptr().unsafe_offset(at), s, v & 1 == 1)
        want += String(unsafe_from_utf8=allele)
        if k == 0:
            want += "\t"
    want += "\t60\t.\tAT="
    for k in range(walks):
        if k > 0:
            want += ","
        for i in range(g.starts[k], g.ends[k]):
            var v = g.visits[i]
            want += String(t"{'<' if v & 1 == 1 else '>'}{v >> 1}")
    want += ";LV=3;PS=>1<99\n"
    var got = _view_of(path, "chrG")
    assert_equal(String(got[byte = got.find("\nS#0#chrG") + 1 :]), want)


def test_snarls_from_the_distance_index() raises:
    """The small graph's distance index (vg 1.73's, version 4) gives vg snarls' snarls: the same start and end visits,
    parents the same snarls."""
    var d = DistIndex("tests/data/micb-kir3dl1.dist")
    assert_equal(d.version, 4)
    var graph = Graph(GBZ, 2)
    check_dist(d, graph)
    var md = graph.parts.metadata()
    for c in ["chr6", "chr19"]:
        var r = ReferencePath(graph, md, "CHM13", c, 2)
        var got = List[Snarl]()
        d.snarls(d.top(Int(r.steps[0] >> 1)), got)
        var want = read_snarls("tests/data/micb-kir3dl1.catalog/snarls/" + c + ".snarls", 2)
        assert_equal(len(got), len(want))

        def key(s: Snarl) -> Tuple[UInt64, UInt64, UInt64, UInt64]:
            var p = canonical(s.parent_start, s.parent_end) if s.has_parent() else (NO_PARENT, NO_PARENT)
            return (s.start, s.end, p[0], p[1])

        def less(a: Snarl, b: Snarl) -> Bool:
            var x = key(a)
            var y = key(b)
            return x[0] < y[0] or (x[0] == y[0] and x[1] < y[1])

        sort(Span(got), less)
        sort(Span(want), less)
        for i in range(len(got)):
            var x = key(got[i])
            var y = key(want[i])
            assert_true(x[0] == y[0] and x[1] == y[1] and x[2] == y[2] and x[3] == y[3])


def test_distance_index_refused() raises:
    """Not a distance index, a truncated one, and one of another graph: each refused."""
    var good = _read_all("tests/data/micb-kir3dl1.dist")
    var bad = String(_tmp("dromos_test_bad.dist"))

    def put(path: String, data: List[UInt8]) raises:
        with open(path, "w") as f:
            f.write_bytes(data)

    put(bad, List[UInt8](Span(good)[: len(good) // 2]))
    var failed = False
    try:
        _ = DistIndex(bad)
    except e:
        failed = String(e).find("damaged") >= 0
    assert_true(failed)
    put(bad, List[UInt8]("not a distance index at all, just text, long enough to be read as one".as_bytes()))
    failed = False
    try:
        _ = DistIndex(bad)
    except e:
        failed = String(e).find("not a distance index") >= 0
    assert_true(failed)
    # The graph's own passes its checks.
    var d = DistIndex("tests/data/micb-kir3dl1.dist")
    var graph = Graph(GBZ, 2)
    check_dist(d, graph)
    assert_true(d.node_count() > 1000)


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
        vcf_header([(String("X"), 5)]).endswith(
            "##contig=<ID=X,length=5>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        )
    )


def _hex(data: Span[UInt8, _]) -> String:
    var digits = "0123456789abcdef".as_bytes()
    var out = List[UInt8]()
    for b in data:
        out.append(digits[Int(b >> 4)])
        out.append(digits[Int(b & 15)])
    return String(unsafe_from_utf8=out)


def _ops(text: StaticString) -> List[CigarOp]:
    """A CIGAR string's operations."""
    var out = List[CigarOp]()
    var n = 0
    for b in text.as_bytes():
        if b >= UInt8(ord("0")) and b <= UInt8(ord("9")):
            n = n * 10 + Int(b - UInt8(ord("0")))
        else:
            out.append(CigarOp(b, n))
            n = 0
    return out^


def _mutated(src: List[UInt8], rate: Int) -> List[UInt8]:
    """`src` with deletions, substitutions (N and lower case among them) and insertions, each at `rate` per mil."""
    var out = List[UInt8]()
    var letters = "ACGTacgtN".as_bytes()
    for b in src:
        var x = Int(random_ui64(0, 999))
        if x < rate:
            continue
        if x < 2 * rate:
            out.append(letters[Int(random_ui64(0, 8))])
            continue
        if x < 3 * rate:
            for _ in range(Int(random_ui64(1, 30))):
                out.append(letters[Int(random_ui64(0, 3))])
        out.append(b)
    return out^


def test_band_filled_in_vectors_as_cell_by_cell() raises:
    """`Banded.fill` (a row in vectors) gives `fill_reference`'s traceback bytes and last row, for any costs, bands
    and shapes: reads against a whole, a longer and a shorter reference, unrelated ones, empty ones."""
    seed(42)
    var letters = "ACGT".as_bytes()
    for trial in range(600):
        var n = Int(random_ui64(0, 400)) if trial % 10 != 0 else Int(random_ui64(0, 3000))
        var r = List[UInt8]()
        for _ in range(n):
            r.append(letters[Int(random_ui64(0, 3))])
        var q = _mutated(r, Int(random_ui64(0, 120)))
        if trial % 7 == 0:
            q = List[UInt8]()
            for _ in range(Int(random_ui64(0, 300))):
                q.append(letters[Int(random_ui64(0, 3))])
        var c = READ_COSTS if trial % 2 == 0 else ASM_COSTS
        if trial % 11 == 0:
            c = Costs(Int(random_ui64(0, 9)), Int(random_ui64(0, 9)), Int(random_ui64(0, 9)))
        if trial % 3 == 1:
            for _ in range(Int(random_ui64(0, 200))):
                r.append(letters[Int(random_ui64(0, 3))])
        if trial % 3 == 2 and len(r) > 10:
            r.resize(len(r) // 2, 0)
        var halves: List[Int] = [half_band(len(q), len(r)), Int(random_ui64(0, 40)), 0, 1, 17]
        for h in halves:
            var a = Banded()
            var b = Banded()
            var filled = a.fill(Span(q), Span(r), h, c)
            assert_equal(filled, b.fill_reference(Span(q), Span(r), h, c))
            if not filled:
                continue
            assert_equal(len(a.tb), len(b.tb))
            for k in range(len(a.tb)):
                if a.tb[k] != b.tb[k]:
                    assert_equal(Int(a.tb[k]), Int(b.tb[k]), String(t"trial {trial}, band {h}, cell {k}"))
            assert_equal(len(a.last), len(b.last))
            for k in range(len(a.last)):
                assert_equal(a.last[k], b.last[k])


def _clip_lengths_by_base(score: List[Int]) -> Tuple[Int, Int]:
    """clip_lengths base by base, as pgcatalog computes it."""
    var n = len(score)
    var total = 0
    var best = 0
    var left = 0
    for i in range(n):
        if score[i] == CLIP_MATCH and total < best:
            best = total
            left = i
        total += score[i]
    if -best <= 5:
        left = 0
    total = 0
    best = 0
    var right = n
    for i in range(n - 1, left - 1, -1):
        if score[i] == CLIP_MATCH and total < best:
            best = total
            right = i + 1
        total += score[i]
    if -best <= 5:
        right = n
    return (left, n - right)


def test_clip_lengths_in_vectors_as_base_by_base() raises:
    """Vectors of scores give the clips a base at a time gives, ties and all, whatever the length."""
    seed(7)
    for trial in range(4000):
        var n = Int(random_ui64(0, 70)) if trial % 4 != 0 else Int(random_ui64(0, 3000))
        var score = List[Int]()
        var bad = Int(random_ui64(0, 60))
        for _ in range(n):
            var x = Int(random_ui64(0, 99))
            if x < bad:
                score.append(CLIP_MISMATCH if x % 3 != 0 else -Int(random_ui64(5, 30)))
            elif x < bad + 5:
                score.append(0)
            else:
                score.append(CLIP_MATCH)
        var got = clip_lengths(score)
        var want = _clip_lengths_by_base(score)
        assert_equal(got[0], want[0], String(t"trial {trial}"))
        assert_equal(got[1], want[1], String(t"trial {trial}"))


def test_nt16_as_htslib() raises:
    assert_equal(nt16(UInt8(ord("A"))), 1)
    assert_equal(nt16(UInt8(ord("c"))), 2)
    assert_equal(nt16(UInt8(ord("G"))), 4)
    assert_equal(nt16(UInt8(ord("T"))), 8)
    assert_equal(nt16(UInt8(ord("U"))), 8)
    assert_equal(nt16(UInt8(ord("N"))), 15)
    assert_equal(nt16(UInt8(ord("R"))), 5)
    assert_equal(nt16(UInt8(ord("="))), 0)
    assert_equal(nt16(UInt8(ord("X"))), 15)
    assert_equal(nt16(UInt8(ord("2"))), 4)


def test_bam_records_as_samtools() raises:
    """Records as `samtools view -b` encodes the same SAM lines (htslib 1.24; the bytes below are its output)."""
    var w = RecordWriter()
    var out = List[UInt8]()
    var name = "read/1".as_bytes()
    var qual = '!"#$%&5?@ABCDEFGH'.as_bytes()
    w.put(
        out,
        name.unsafe_ptr(),
        len(name),
        16,
        0,
        99,
        60,
        _ops("3S5M2I4M1D3M"),
        "ACGTNACGTRACGTACG".as_bytes(),
        qual.unsafe_ptr(),
        0,
        False,
    )
    assert_equal(
        _hex(Span(out)),
        (
            "590000000000000063000000073c49120600100011000000ffffffffffffffff00000000726561642f31003400000050000000"
            "210000004000000012000000300000001248f1248512481240000102030405141e1f2021222324252627"
        ),
    )
    out.clear()
    # No qualities ('I' throughout), no reference bases (the bin of [pos, pos + 1)), U written as T.
    var r2 = "r2".as_bytes()
    w.put(out, r2.unsafe_ptr(), 2, 0, 0, 19999, 255, _ops("2S3I"), "ACGTU".as_bytes(), r2.unsafe_ptr(), -1, False)
    assert_equal(
        _hex(Span(out)),
        "33000000000000001f4e000003ff4a120200000005000000ffffffffffffffff0000000072320024000000310000001248802828282828",
    )
    assert_equal(record_end(out.unsafe_ptr()), 20000)
    out.clear()
    # Qualities reversed: the last character first.
    var q3 = "ABC".as_bytes()
    w.put(out, r2.unsafe_ptr(), 2, 16, 1, 5, 0, _ops("3M"), "ACG".as_bytes(), q3.unsafe_ptr(), 0, True)
    var n = len(out)
    assert_equal(Int(out[n - 3]), ord("C") - 33)
    assert_equal(Int(out[n - 1]), ord("A") - 33)
    assert_equal(record_end(out.unsafe_ptr()), 8)


def test_bam_long_cigar_goes_to_cg() raises:
    """More than 65,535 operations: <length>S<reference length>N in place and the CIGAR in CG:B,I, as bam_write1."""
    var cigar = List[CigarOp]()
    var seq = List[UInt8]()
    for _ in range(35000):
        cigar.append(CigarOp(UInt8(ord("M")), 1))
        cigar.append(CigarOp(UInt8(ord("I")), 1))
        seq.append(UInt8(ord("A")))
        seq.append(UInt8(ord("C")))
    var w = RecordWriter()
    var out = List[UInt8]()
    var name = "long".as_bytes()
    w.put(out, name.unsafe_ptr(), 4, 0, 0, 1000, 60, cigar, Span(seq), name.unsafe_ptr(), -1, False)
    var p = out.unsafe_ptr()
    assert_equal(Int(p[16]) | Int(p[17]) << 8, 2)
    assert_equal(record_end(p), 1000 + 35000)
    var tag = 36 + 5 + 8 + 35000 + 70000
    assert_equal(String(unsafe_from_utf8=Span(out)[tag : tag + 4]), "CGBI")
    assert_equal(len(out), tag + 8 + 4 * 70000)


def test_bam_index_kinds() raises:
    assert_equal(csi_levels(248387328, 14), 5)
    assert_equal(csi_levels(600_000_000, 14), 6)
    assert_false(needs_csi([1 << 29]))
    assert_true(needs_csi([10, (1 << 29) + 1]))
    assert_equal(BamIndex(["a"], [10], False).suffix(), ".bai")
    assert_equal(BamIndex(["a"], [10], True).suffix(), ".csi")


def _inflate_file(path: String) raises -> List[UInt8]:
    var data: List[UInt8]
    with open(path, "r") as f:
        data = f.read_bytes()
    var blocks = List[Int]()
    var got = scan_blocks(Span(data), 0, blocks, path)
    assert_equal(got[0], len(data))
    var out = List[UInt8](length=got[1], fill=0)
    assert_true(
        inflate_blocks(load_libdeflate(), data.unsafe_ptr(), blocks.unsafe_ptr(), out.unsafe_ptr(), 0, len(blocks) // 4)
    )
    return out^


def test_bam_sorted_through_spills() raises:
    """Records staged out of order, spilled at every step, come back sorted as samtools sorts: contig, position,
    reverse strand, then input order; the BAM's header first."""
    var w = RecordWriter()
    var names: List[String] = ["c1", "c2"]
    var lengths: List[Int] = [3_000_000, 1000]
    var header = bam_header("@HD\tVN:1.6\tSO:coordinate\n", names, lengths)
    # Every record a frame of its own, and the frames spilled after every step but the last.
    var staging = BamStaging(2, 3, 1, "build", frame=1)
    # (contig, position, flag) in input order; names number them.
    var input: List[Tuple[Int, Int, Int]] = [
        (1, 50, 0),
        (0, 2_500_000, 16),
        (0, 7, 16),
        (0, 7, 0),
        (1, 50, 16),
        (0, 2_500_000, 0),
        (0, 7, 0),
        (1, 10, 0),
        (0, 1_200_000, 0),
        (0, 7, 16),
    ]
    for i in range(len(input)):
        var recs = List[UInt8]()
        var name = String(t"r{i}")
        var q = "IIII".as_bytes()
        w.put(
            recs,
            name.unsafe_ptr(),
            name.byte_length(),
            input[i][2],
            input[i][0],
            input[i][1],
            60,
            _ops("4M"),
            "ACGT".as_bytes(),
            q.unsafe_ptr(),
            0,
            False,
        )
        staging.make_room(recs.unsafe_ptr(), len(recs))
        for share in range(3):
            staging.take_share(share, recs.unsafe_ptr(), len(recs))
        if i < len(input) - 1:
            staging.after_shares()
    assert_equal(staging.spills, len(input) - 1)
    var path = String("build/test_sorted.bam")
    write_bam(staging, header, BamIndex(names, lengths, False), path, 1, 2)
    var raw = _inflate_file(path)
    assert_equal(String(unsafe_from_utf8=Span(raw)[0:4]), "BAM\x01")
    var at = len(header)
    var order = List[Int]()
    while at < len(raw):
        var size = Int(raw[at]) | Int(raw[at + 1]) << 8
        var name_len = Int(raw[at + 12])
        order.append(Int(String(unsafe_from_utf8=Span(raw)[at + 37 : at + 36 + name_len - 1])))
        at += 4 + size
    var want: List[Int] = [3, 6, 2, 9, 8, 5, 1, 7, 0, 4]
    assert_equal(len(order), len(want))
    for k in range(len(want)):
        assert_equal(order[k], want[k])
    with open(path + ".bai", "r") as f:
        assert_equal(String(unsafe_from_utf8=f.read_bytes()[0:4]), "BAI\x01")


def main() raises:
    TestSuite.discover_tests[__functions_in_module()]().run()
