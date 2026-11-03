"""vg's GAM alignments decoded, and written as the GAF lines `vg convert -G` writes for them.

A GAM file is BGZF holding groups: a varint64 count (the messages plus one), a varint32-long type tag ("GAM"),
then each message as a varint32 length and its bytes (libvgio's message_emitter.hpp). Each message is a protobuf
`vg.Alignment` (libvgio's vg.proto). `gaf_line` follows libvgio's `alignment_to_gaf` (src/alignment_io.cpp) with
convert's settings -- cs and bq on, no segment translation -- so its lines match vg's byte for byte.
"""

from std.math import floor
from std.memory import bitcast

from dromos.nodes import Nodes
from dromos.output import InlineDigits, put_padded

comptime TAG = "GAM"
"""The type tag of a GAM group of alignments."""
comptime MAX_TAG_LENGTH = 32
"""The longest first item still read as a type tag (registry tags are a few capitals; a message is far longer)."""


@always_inline
def varint(p: ImmPointer[UInt8, _], mut at: Int, end: Int) raises -> UInt64:
    """The varint at `at`, which moves past it."""
    if at < end:
        var c = p[unsafe_offset=at]
        if c < 0x80:
            at += 1
            return UInt64(c)
    return _varint_long(p, at, end)


def _varint_long(p: ImmPointer[UInt8, _], mut at: Int, end: Int) raises -> UInt64:
    """`varint`, for one of more than a byte (or none)."""
    var x: UInt64 = 0
    var shift: UInt64 = 0
    while True:
        if at >= end:
            raise Error("GAM: truncated varint")
        var c = p[unsafe_offset=at]
        at += 1
        x |= UInt64(c & 0x7F) << shift
        if c < 0x80:
            return x
        shift += 7
        if shift > 63:
            raise Error("GAM: varint too long")


@always_inline
def skip(p: ImmPointer[UInt8, _], mut at: Int, end: Int, wire: Int) raises:
    """Moves `at` past a field's value of wire type `wire`."""
    if wire == 0:
        _ = varint(p, at, end)
    elif wire == 1:
        at += 8
    elif wire == 2:
        at += Int(varint(p, at, end))
    elif wire == 5:
        at += 4
    else:
        raise Error(t"GAM: wire type {wire}")
    if at > end:
        raise Error("GAM: truncated field")


@always_inline
def int32_of(x: UInt64) -> Int:
    """A protobuf int32: the low 32 bits, signed (negative values come as 10-byte varints)."""
    return Int(x.cast[DType.int32]())


@fieldwise_init
struct Edit(TrivialRegisterPassable):
    var from_length: Int
    var to_length: Int
    var sequence_at: Int
    var sequence_length: Int


@fieldwise_init
struct Mapping(TrivialRegisterPassable):
    var node: Int
    var offset: Int
    var reverse: Bool
    var first_edit: Int
    var edit_count: Int


struct Alignment(Movable):
    """One decoded `vg.Alignment`: its fields as spans of the message, and its path's mappings and edits."""

    var name_at: Int
    var name_length: Int
    var sequence_at: Int
    var sequence_length: Int
    var quality_at: Int
    var quality_length: Int
    var mapping_quality: Int
    var score: Int
    var identity: Float64
    var has_path: Bool
    var mappings: List[Mapping]
    var edits: List[Edit]

    def __init__(out self):
        self.name_at = 0
        self.name_length = 0
        self.sequence_at = 0
        self.sequence_length = 0
        self.quality_at = 0
        self.quality_length = 0
        self.mapping_quality = 0
        self.score = 0
        self.identity = 0.0
        self.has_path = False
        self.mappings = List[Mapping]()
        self.edits = List[Edit]()

    def reset(mut self):
        """Empties every field, keeping the lists' memory."""
        self.name_at = 0
        self.name_length = 0
        self.sequence_at = 0
        self.sequence_length = 0
        self.quality_at = 0
        self.quality_length = 0
        self.mapping_quality = 0
        self.score = 0
        self.identity = 0.0
        self.has_path = False
        self.mappings.clear()
        self.edits.clear()

    def decode(mut self, p: ImmPointer[UInt8, _], start: Int, end: Int) raises:
        """Reads the message at `p[start:end]`, reusing this alignment's lists."""
        self.reset()
        var at = start
        while at < end:
            var key = Int(varint(p, at, end))
            var field = key >> 3
            var wire = key & 7
            if wire == 2 and (field == 1 or field == 2 or field == 3 or field == 4):
                var n = Int(varint(p, at, end))
                if at + n > end:
                    raise Error("GAM: truncated field")
                if field == 1:
                    self.sequence_at = at
                    self.sequence_length = n
                elif field == 2:
                    self.has_path = True
                    self._path(p, at, at + n)
                elif field == 3:
                    self.name_at = at
                    self.name_length = n
                else:
                    self.quality_at = at
                    self.quality_length = n
                at += n
            elif wire == 0 and field == 5:
                self.mapping_quality = int32_of(varint(p, at, end))
            elif wire == 0 and field == 6:
                self.score = int32_of(varint(p, at, end))
            elif wire == 1 and field == 16:
                if at + 8 > end:
                    raise Error("GAM: truncated field")
                var bits: UInt64 = 0
                for k in range(8):
                    bits |= UInt64(p[unsafe_offset=at + k]) << UInt64(8 * k)
                self.identity = bitcast[DType.float64](bits)
                at += 8
            else:
                skip(p, at, end, wire)

    def _path(mut self, p: ImmPointer[UInt8, _], start: Int, end: Int) raises:
        var at = start
        while at < end:
            var key = Int(varint(p, at, end))
            if key == (2 << 3) | 2:
                var n = Int(varint(p, at, end))
                self._mapping(p, at, at + n)
                at += n
            else:
                skip(p, at, end, key & 7)

    def _mapping(mut self, p: ImmPointer[UInt8, _], start: Int, end: Int) raises:
        var m = Mapping(0, 0, False, len(self.edits), 0)
        var at = start
        while at < end:
            var key = Int(varint(p, at, end))
            if key == (1 << 3) | 2:
                var n = Int(varint(p, at, end))
                var stop = at + n
                while at < stop:
                    var k = Int(varint(p, at, stop))
                    if k == (1 << 3) | 0:
                        m.node = Int(varint(p, at, stop))
                    elif k == (2 << 3) | 0:
                        m.offset = Int(varint(p, at, stop))
                    elif k == (4 << 3) | 0:
                        m.reverse = varint(p, at, stop) != 0
                    else:
                        skip(p, at, stop, k & 7)
            elif key == (2 << 3) | 2:
                var n = Int(varint(p, at, end))
                var stop = at + n
                var e = Edit(0, 0, 0, 0)
                while at < stop:
                    var k = Int(varint(p, at, stop))
                    if k == (1 << 3) | 0:
                        e.from_length = int32_of(varint(p, at, stop))
                    elif k == (2 << 3) | 0:
                        e.to_length = int32_of(varint(p, at, stop))
                    elif k == (3 << 3) | 2:
                        var s = Int(varint(p, at, stop))
                        e.sequence_at = at
                        e.sequence_length = s
                        at += s
                    else:
                        skip(p, at, stop, k & 7)
                self.edits.append(e)
                m.edit_count += 1
            else:
                skip(p, at, end, key & 7)
        self.mappings.append(m)


@always_inline
def _is_match(e: Edit) -> Bool:
    return e.from_length == e.to_length and e.sequence_length == 0


@always_inline
def _is_sub(e: Edit) -> Bool:
    return e.from_length == e.to_length and e.sequence_length > 0


@always_inline
def _is_insertion(e: Edit) -> Bool:
    return e.from_length == 0 and e.to_length > 0 and e.sequence_length > 0


@always_inline
def _is_deletion(e: Edit) -> Bool:
    return e.from_length > 0 and e.to_length == 0


@always_inline
def put_int(mut out: List[UInt8], value: Int):
    """`value` in decimal."""
    if value >= 0:
        var at = len(out)
        if out.capacity() < at + 20:
            out.reserve(max(at + 20, 2 * out.capacity()))
        out.resize(unsafe_uninit_length=at + put_padded(out.unsafe_ptr().unsafe_offset(at), value, 0))
        return
    var digits = InlineDigits(value)
    for i in range(digits.count):
        out.append(digits.bytes[digits.count - 1 - i])


@always_inline
def put_text(mut out: List[UInt8], text: StaticString):
    out.extend(text.as_bytes())


@always_inline
def put_bytes(mut out: List[UInt8], p: ImmPointer[UInt8, _], at: Int, length: Int):
    out.extend(Span(unsafe_ptr=p.unsafe_offset(at), length=length))


def _put_divergence(mut out: List[UInt8], identity: Float64):
    """vg's dv: `floor((1 - identity) * 10000 + 0.5) / 10000` as an ostream prints a double (six significant
    digits, %g): a multiple of 0.0001 at most 1, so at most four decimals, trailing zeros dropped."""
    var k = Int(floor((1.0 - identity) * 10000.0 + 0.5))
    if k < 0:
        out.append(UInt8(ord("-")))
        k = -k
    put_int(out, k // 10000)
    var frac = k % 10000
    if frac == 0:
        return
    out.append(UInt8(ord(".")))
    var digits = SIMD[DType.uint8, 4](0)
    for i in range(4):
        digits[3 - i] = UInt8(ord("0") + frac % 10)
        frac //= 10
    var last = 3
    while digits[last] == UInt8(ord("0")):
        last -= 1
    for i in range(last + 1):
        out.append(digits[i])


struct _Cs:
    """The cs string as alignment_to_gaf builds it: runs of matches merged, deletions run together."""

    var text: List[UInt8]
    var matches: Int
    var deleting: Bool

    def __init__(out self):
        self.text = List[UInt8]()
        self.matches = 0
        self.deleting = False

    @always_inline
    def flush(mut self):
        if self.matches > 0:
            self.text.append(UInt8(ord(":")))
            put_int(self.text, self.matches)
            self.matches = 0

    @always_inline
    def delete_mark(mut self):
        """The "-" a deletion starts with, unless one runs on."""
        if not self.deleting:
            self.text.append(UInt8(ord("-")))
        self.deleting = True


def gaf_line[
    N: Nodes, //, tags: Bool = True
](p: ImmPointer[UInt8, _], a: Alignment, nodes: N, mut out: List[UInt8]) raises:
    """Appends the GAF line (with its newline) vg convert writes for alignment `a` of the message at `p`; with
    `tags` false, its 12 columns alone (no cs, bq, AS or dv, and no node letters needed)."""
    var walk = List[Int]()
    _ = _gaf_line[tags, False](p, a, nodes, out, walk)


def slim_line[
    N: Nodes, //
](p: ImmPointer[UInt8, _], a: Alignment, nodes: N, mut out: List[UInt8], mut walk: List[Int]) raises -> Tuple[Int, Int]:
    """`gaf_line` without tags, also giving the line's walk (its steps `node << 1 | reverse`, column 6) and its path
    start and end (columns 8 and 9; (0, 0) for an unaligned read, whose columns are `*`)."""
    return _gaf_line[False, True](p, a, nodes, out, walk)


def _gaf_line[
    N: Nodes, //, tags: Bool, steps: Bool
](p: ImmPointer[UInt8, _], a: Alignment, nodes: N, mut out: List[UInt8], mut walk: List[Int]) raises -> Tuple[Int, Int]:
    comptime if steps:
        walk.clear()
    put_bytes(out, p, a.name_at, a.name_length)
    if not a.has_path or len(a.mappings) == 0:
        # Unaligned: the read as one insertion, no place.
        out.append(UInt8(ord("\t")))
        put_int(out, a.sequence_length)
        put_text(out, "\t0\t")
        put_int(out, a.sequence_length)
        put_text(out, "\t*\t*\t*\t*\t*\t*\t*\t")
        put_int(out, a.mapping_quality)
        comptime if tags:
            put_text(out, "\tcs:Z:+")
            put_bytes(out, p, a.sequence_at, a.sequence_length)
        out.append(UInt8(ord("\n")))
        return (0, 0)
    if tags and not nodes.has_letters():
        raise Error("GAM: the node table has no letters, which cs needs")
    var n_maps = len(a.mappings)
    var final_mapping = n_maps - 1
    if final_mapping > 0:
        ref last = a.mappings[final_mapping]
        if last.edit_count == 1 and _is_insertion(a.edits[last.first_edit]):
            final_mapping -= 1
    var query_length = a.sequence_length
    if query_length == 0:
        for e in a.edits:
            query_length += e.to_length
    out.append(UInt8(ord("\t")))
    put_int(out, query_length)
    put_text(out, "\t0\t")
    put_int(out, query_length)
    put_text(out, "\t+\t")
    var cs = _Cs()
    var path_length = 0
    var path_start = 0
    var path_end = 0
    var matches = 0
    var prev_offset = 0
    for i in range(n_maps):
        ref m = a.mappings[i]
        var start = m.offset
        var offset = start
        var node_length = nodes.length(m.node)
        var skip_step = False
        if i > 0 and start > 0:
            ref pm = a.mappings[i - 1]
            if start == prev_offset and m.node == pm.node and m.reverse == pm.reverse:
                skip_step = True
            else:
                # A split mapping: the start of the node it arrives at is deleted.
                comptime if tags:
                    var deleted_from = prev_offset if m.node == pm.node else 0
                    if start > deleted_from:
                        cs.flush()
                        cs.delete_mark()
                        for k in range(deleted_from, start):
                            cs.text.append(nodes.letter(m.node, m.reverse, k))
        for j in range(m.first_edit, m.first_edit + m.edit_count):
            var e = a.edits[j]
            if _is_match(e):
                matches += e.from_length
                comptime if tags:
                    cs.matches += e.from_length
                    cs.deleting = False
            elif tags:
                cs.flush()
                if _is_sub(e):
                    for k in range(e.from_length):
                        cs.text.append(UInt8(ord("*")))
                        cs.text.append(nodes.letter(m.node, m.reverse, offset + k))
                        cs.text.append(p[unsafe_offset=e.sequence_at + k])
                    cs.deleting = False
                elif _is_deletion(e):
                    cs.delete_mark()
                    for k in range(e.from_length):
                        cs.text.append(nodes.letter(m.node, m.reverse, offset + k))
                elif _is_insertion(e):
                    cs.text.append(UInt8(ord("+")))
                    put_bytes(cs.text, p, e.sequence_at, e.sequence_length)
                    cs.deleting = False
            offset += e.from_length
        if i == 0:
            path_start = start
        elif i > final_mapping:
            # giraffe's empty mapping for a softclip at the end.
            skip_step = True
        if i < n_maps - 1 and offset != node_length:
            ref nm = a.mappings[i + 1]
            if m.node != nm.node or m.reverse != nm.reverse:
                # Hopping off the middle of a node: the rest of it is deleted.
                comptime if tags:
                    cs.flush()
                    cs.delete_mark()
                    for k in range(offset, node_length):
                        cs.text.append(nodes.letter(m.node, m.reverse, k))
            else:
                # A duplicate mapping of the same node: no step of its own.
                skip_step = True
        if not skip_step:
            path_length += node_length
            if i == 0:
                path_start = start
            out.append(UInt8(ord("<")) if m.reverse else UInt8(ord(">")))
            put_int(out, m.node)
            comptime if steps:
                walk.append(m.node << 1 | Int(m.reverse))
        if i == final_mapping:
            path_end = path_length - (node_length - offset)
        prev_offset = offset
    cs.flush()
    out.append(UInt8(ord("\t")))
    put_int(out, path_length)
    out.append(UInt8(ord("\t")))
    put_int(out, path_start)
    out.append(UInt8(ord("\t")))
    put_int(out, path_end)
    out.append(UInt8(ord("\t")))
    put_int(out, matches)
    out.append(UInt8(ord("\t")))
    put_int(out, max(path_end - path_start, query_length))
    out.append(UInt8(ord("\t")))
    put_int(out, a.mapping_quality)
    comptime if not tags:
        out.append(UInt8(ord("\n")))
        return (path_start, path_end)
    # The optional fields in key order (gafkluge writes its std::map): AS, bq, cs, dv.
    if a.score > 0:
        put_text(out, "\tAS:i:")
        put_int(out, a.score)
    if a.quality_length > 0:
        put_text(out, "\tbq:Z:")
        for k in range(a.quality_length):
            out.append(p[unsafe_offset=a.quality_at + k] + 33)
    put_text(out, "\tcs:Z:")
    out.extend(Span(cs.text))
    if a.identity > 0:
        put_text(out, "\tdv:f:")
        _put_divergence(out, a.identity)
    out.append(UInt8(ord("\n")))
    return (path_start, path_end)


@always_inline
def _varint_in(p: ImmPointer[UInt8, _], mut at: Int, n: Int, mut x: Int) raises -> Bool:
    """Reads the varint at `at` into `x` and moves past it; false, `at` unmoved, if `p[:n]` ends inside it."""
    var v: UInt64 = 0
    var shift: UInt64 = 0
    var k = at
    while k < n:
        var c = p[unsafe_offset=k]
        k += 1
        v |= UInt64(c & 0x7F) << shift
        if c < 0x80:
            x = Int(v)
            at = k
            return True
        shift += 7
        if shift > 63:
            raise Error("GAM: varint too long")
    return False


struct Framer(Movable):
    """Finds the messages of a GAM read a piece at a time, a group's messages handed out as soon as each is whole:
    what is left of a piece is at most an unfinished message or group header, read again with more data."""

    var remaining: Int
    """The items of the group being read still to come."""
    var is_gam: Bool
    """Whether they are alignments (a GAM group's, or an untagged group's) or another tag's, skipped."""

    def __init__(out self):
        self.remaining = 0
        self.is_gam = False

    def frame(mut self, data: Span[UInt8, _], mut messages: List[Int]) raises -> Int:
        """Appends (start, end) of each whole message at the front of `data`; returns the bytes they (and their
        groups' headers) take."""
        var p = data.unsafe_ptr()
        var n = len(data)
        var tag = TAG.as_bytes()
        var done = 0
        while done < n:
            var at = done
            if self.remaining > 0:
                var length = 0
                if not _varint_in(p, at, n, length) or at + length > n:
                    return done
                if self.is_gam:
                    messages.append(at)
                    messages.append(at + length)
                self.remaining -= 1
                done = at + length
                continue
            var count = 0
            var tag_length = 0
            if not _varint_in(p, at, n, count) or (count > 0 and not _varint_in(p, at, n, tag_length)):
                return done
            if count == 0:
                done = at
                continue
            if at + tag_length > n:
                return done
            # The first item is a type tag when it looks like one (libvgio's registry tags: short, capitals, digits,
            # _); a GAM group's messages are alignments, another tag's (giraffe's PARAMS) are skipped, and an untagged
            # group (an older GAM) is alignments from its first item.
            var looks_like_tag = tag_length > 0 and tag_length <= MAX_TAG_LENGTH
            for k in range(tag_length if looks_like_tag else 0):
                var c = p[unsafe_offset=at + k]
                looks_like_tag = looks_like_tag and (
                    (c >= UInt8(ord("A")) and c <= UInt8(ord("Z")))
                    or (c >= UInt8(ord("0")) and c <= UInt8(ord("9")))
                    or c == UInt8(ord("_"))
                )
            var is_gam = looks_like_tag and tag_length == len(tag)
            for k in range(tag_length if is_gam else 0):
                is_gam = is_gam and p[unsafe_offset=at + k] == tag[k]
            if not looks_like_tag:
                # No tag: the first item was a message's length, so the group holds `count` messages.
                at = done
                _ = _varint_in(p, at, n, count)
                self.remaining = count
                self.is_gam = True
            else:
                at += tag_length
                self.remaining = count - 1
                self.is_gam = is_gam
            done = at
        return done


def frame(data: Span[UInt8, _], mut messages: List[Int]) raises -> Int:
    """`Framer.frame` from the start of a GAM."""
    var f = Framer()
    return f.frame(data, messages)
