"""A GAF line as pgcat reads it (pgphase's tools/pgcat, gaf.rs): its columns found, classified (placed on the
reference, or why not), and kept slim (its 12 columns and the tags pgphase reads)."""

from std.ffi import external_call

from dromos.refmap import Anchor, Placed, RefMap

comptime SKIP = 0
"""A line that is no alignment: empty, a comment (#) or a header (@)."""
comptime PLACED = 1
comptime UNMAPPED = 2
comptime UNPLACED = 3
comptime MALFORMED = 4


def append_slim(mut out: List[UInt8], p: ImmPointer[UInt8, _], n: Int):
    """pgcat's append_line with slim: the GAF line `p[:n]`'s first 12 columns, then of the rest only its SLIM_TAGS
    (rq:i:, the rescue's MAPQ; pu:Z:, why a line had no place; rs:Z:, where the rescue put it), in their order. Tabs
    are found with memchr: a tag's bytes after its first five are never looked at."""
    var field = 0
    var s = 0
    while True:
        var hit = external_call["memchr", Int](p.unsafe_offset(s), Int32(ord("\t")), n - s) if s < n else 0
        var e = n if hit == 0 else hit - Int(p)
        var length = e - s
        if field < 12:
            if field > 0:
                out.append(UInt8(ord("\t")))
            out.extend(Span(unsafe_ptr=p.unsafe_offset(s), length=length))
        elif length >= 5 and p[unsafe_offset=s + 2] == UInt8(ord(":")) and p[unsafe_offset=s + 4] == UInt8(ord(":")):
            var a = p[unsafe_offset=s]
            var b = p[unsafe_offset=s + 1]
            var c = p[unsafe_offset=s + 3]
            var keep = (
                (a == UInt8(ord("r")) and b == UInt8(ord("q")) and c == UInt8(ord("i")))
                or (a == UInt8(ord("p")) and b == UInt8(ord("u")) and c == UInt8(ord("Z")))
                or (a == UInt8(ord("r")) and b == UInt8(ord("s")) and c == UInt8(ord("Z")))
            )
            if keep:
                out.append(UInt8(ord("\t")))
                out.extend(Span(unsafe_ptr=p.unsafe_offset(s), length=length))
        field += 1
        if e >= n:
            break
        s = e + 1


def parse_u64(p: ImmPointer[UInt8, _], at: Int, end: Int) -> Int:
    """The decimal number `p[at:end]`, or -1 (pgcat's parse_u64: 1 to 19 digits, nothing else)."""
    if end <= at or end - at > 19:
        return -1
    var v = 0
    for i in range(at, end):
        var c = Int(p[unsafe_offset=i]) - ord("0")
        if c < 0 or c > 9:
            return -1
        v = v * 10 + c
    return v if v >= 0 else -1


def parse_walk(p: ImmPointer[UInt8, _], at: Int, end: Int, mut walk: List[Int]) -> Bool:
    """A GAF walk (`>12<34`) as steps `node << 1 | reverse`; false if it is not one."""
    walk.clear()
    var i = at
    while i < end:
        var c = p[unsafe_offset=i]
        var reverse: Int
        if c == UInt8(ord(">")):
            reverse = 0
        elif c == UInt8(ord("<")):
            reverse = 1
        else:
            return False
        i += 1
        var s = i
        while i < end and p[unsafe_offset=i] != UInt8(ord(">")) and p[unsafe_offset=i] != UInt8(ord("<")):
            i += 1
        var node = parse_u64(p, s, i)
        if node < 0:
            return False
        walk.append(node << 1 | reverse)
    return True


def classify(
    map: RefMap,
    p: ImmPointer[UInt8, _],
    n: Int,
    mut walk: List[Int],
    mut anchors: List[Anchor],
    mut placed: Placed,
    mut why: Int,
    mut aligned: Int,
) -> Int:
    """pgcat's classify, for one GAF line `p[:n]` (no newline): SKIP, MALFORMED, UNMAPPED, UNPLACED (and `why`) or
    PLACED (and `placed`, and the graph bases `aligned`)."""
    if n == 0 or p[unsafe_offset=0] == UInt8(ord("#")) or p[unsafe_offset=0] == UInt8(ord("@")):
        return SKIP
    # The 12 columns' starts, and the end of the 12th.
    var starts = SIMD[DType.int64, 16](0)
    var ends = SIMD[DType.int64, 16](0)
    var field = 0
    var s = 0
    while field < 12:
        var hit = external_call["memchr", Int](p.unsafe_offset(s), Int32(ord("\t")), n - s) if s < n else 0
        var e = n if hit == 0 else hit - Int(p)
        starts[field] = Int64(s)
        ends[field] = Int64(e)
        field += 1
        if e >= n:
            break
        s = e + 1
    if field < 12:
        return MALFORMED
    var strand_at = Int(starts[4])
    var strand_len = Int(ends[4]) - strand_at
    if strand_len == 1 and p[unsafe_offset=strand_at] == UInt8(ord("*")):
        return UNMAPPED
    if strand_len != 1 or (
        p[unsafe_offset=strand_at] != UInt8(ord("+")) and p[unsafe_offset=strand_at] != UInt8(ord("-"))
    ):
        return MALFORMED
    var ts = parse_u64(p, Int(starts[7]), Int(ends[7]))
    var te = parse_u64(p, Int(starts[8]), Int(ends[8]))
    if ts < 0 or te < 0:
        return MALFORMED
    if not parse_walk(p, Int(starts[5]), Int(ends[5]), walk):
        return MALFORMED
    why = map.project(walk, ts, te, anchors, placed)
    if why != 0:
        return UNPLACED
    aligned = te - ts
    return PLACED


def parse_unsigned(p: ImmPointer[UInt8, _], at: Int, end: Int, limit: Int = Int.MAX) -> Int:
    """The number `p[at:end]` as Rust's `str::parse` reads an unsigned integer (an optional +, then one or more
    digits, at most `limit`), or -1."""
    var i = at
    if i < end and p[unsafe_offset=i] == UInt8(ord("+")):
        i += 1
    if i >= end:
        return -1
    var v = 0
    for k in range(i, end):
        var c = Int(p[unsafe_offset=k]) - ord("0")
        if c < 0 or c > 9 or v > (limit - c) // 10:
            return -1
        v = v * 10 + c
    return v
