"""Tests: protobuf varints, GAM group framing, and GAF lines against hand-worked ones."""

from std.memory import bitcast
from std.testing import TestSuite, assert_equal, assert_true

from dromos.gam import Alignment, frame, gaf_line, varint
from dromos.nodes import NodeTable


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


def test_frame_keeps_whole_groups() raises:
    # Two groups: 2 messages, then 1; the second cut short is left for the next read.
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
    var cut = stream.copy()
    _ = cut.pop()
    messages.clear()
    var first = 2 + 3 + 2 * 4
    assert_equal(frame(Span(cut), messages), first)
    assert_equal(len(messages), 4)


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


def main() raises:
    TestSuite.discover_tests[__functions_in_module()]().run()
