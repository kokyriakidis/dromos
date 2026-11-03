"""Tests: protobuf varints, GAM group framing, GAF lines against hand-worked ones, placement on the reference as
pgcat places alignments, and tabix bins as htslib numbers them."""

from std.memory import bitcast
from std.testing import TestSuite, assert_equal, assert_true

from dromos.gam import Alignment, Framer, frame, gaf_line, varint
from dromos.nodes import NodeTable
from dromos.refmap import AMBIGUOUS, INCONSISTENT, MISSING_NODE, NO_ANCHOR, Anchor, Placed, RefMap
from dromos.tabix import N_BINS, bin_level, bin_parent, reg2bin


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


def _nodes_file() raises -> String:
    """A node table of nodes 1 ACGTACGT, 2 GGG and 3 TTTTT (node 0 absent)."""
    var path = String("/tmp/dromos_test_nodes.bin")
    var t = _text("DROMOSN1")
    var lengths: List[Int] = [0, 8, 3, 5]
    for k in range(8):
        t.append(UInt8((4 >> (8 * k)) & 0xFF))
    for n in lengths:
        for k in range(4):
            t.append(UInt8((n >> (8 * k)) & 0xFF))
    t.append(1)
    for _ in range(7):
        t.append(0)
    var at = 0
    for n in [0, 8, 3, 5, 0]:
        for k in range(8):
            t.append(UInt8((at >> (8 * k)) & 0xFF))
        at += n
    t.extend(Span(_text("ACGTACGTGGGTTTTT")))
    with open(path, "w") as f:
        f.write_bytes(Span(t))
    return path


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
    # --slim: the 12 columns alone.
    var slim = List[UInt8]()
    gaf_line[False](msg.unsafe_ptr(), a, nodes, slim)
    assert_true(slim == _text("r1\t13\t0\t13\t+\t>1>2>3\t16\t0\t13\t12\t13\t60\n"))


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


def _le(mut out: List[UInt8], v: Int, bytes: Int):
    for k in range(bytes):
        out.append(UInt8((v >> (8 * k)) & 0xFF))


def _refmap_file(twice: Bool) raises -> String:
    """pgcat's test map: reference path >1 >2 >3 (lengths 10, 5, 10) at 100 on chrA; node 4 (5 bp) off the reference;
    node 0 absent. With `twice`, node 2 is also on a second path, at 500 of chrB."""
    var path = String(t"/tmp/dromos_test_refmap_{Int(twice)}.bin")
    var t = _text("PGCATRM1")
    _le(t, 5, 8)
    t.extend(Span(_text("CHM13")))
    _le(t, 5, 8)
    for n in [0xFFFFFFFF, 10, 5, 10, 5]:
        _le(t, n, 4)
    for f in [0xFFFFFFFF, 0, 1, 2, 0xFFFFFFFF]:
        _le(t, f, 4)
    var placements: List[List[Int]] = [[0, 0, 100, 0], [0, 1, 110, 0], [0, 2, 115, 0]]
    if twice:
        placements.append([1, 0, 500, 0])
    _le(t, len(placements), 8)
    for pl in placements:
        _le(t, pl[0], 4)
        _le(t, pl[1], 4)
        _le(t, pl[2], 8)
        _le(t, pl[3], 1)
    _le(t, 1 if twice else 0, 8)
    if twice:
        _le(t, 2, 8)
        _le(t, 1, 8)
        _le(t, 3, 4)
    _le(t, 2, 8)
    _le(t, 0, 4)
    _le(t, 1, 4)
    _le(t, 2, 8)
    _le(t, 4, 8)
    t.extend(Span(_text("chrA")))
    _le(t, 4, 8)
    t.extend(Span(_text("chrB")))
    with open(path, "w") as f:
        f.write_bytes(Span(t))
    return path


def _project(m: RefMap, walk: List[Int], ts: Int, te: Int) -> Tuple[Int, Int, Int, Int]:
    """Why (0 placed), contig, start, end."""
    var anchors = List[Anchor]()
    var p = Placed(-1, -1, -1)
    var why = m.project(walk, ts, te, anchors, p)
    return (why, p.contig, p.start, p.end)


def _step(node: Int, reverse: Bool) -> Int:
    return node << 1 | Int(reverse)


def test_projection_as_pgcat() raises:
    var m = RefMap(_refmap_file(False))
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
    var m = RefMap(_refmap_file(True))
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


def main() raises:
    TestSuite.discover_tests[__functions_in_module()]().run()
