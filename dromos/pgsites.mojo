"""A contig's site catalog (pgcatalog's `.pgsites`, its pgsites.rs) read for its snarl windows, as pgcatalog pileup's
`load_windows` makes them: each record whose ID names its snarl's two boundary nodes, both visited once by the
contig, gives the reference bases between them (at most MAX_SITE); windows inside another are dropped, the rest padded
by WINDOW_PAD and those that overlap merged up to MAX_WINDOW.

The catalog (integers little-endian): "PGSITES2"; u64 header length and the zstd VCF header; zstd blocks of VCF
lines; the index (u64 count, then per block nine u64s: first_pos, last_pos, max_end, offset, compressed length, raw
length, records, min_allele, max_walk); u64 the index's offset. Every record is read (no region)."""

from std.ffi import OwnedDLHandle, external_call

from dromos.bgzf import MappedFile, load_zstd, zstd_decompress
from dromos.contig import Contig
from dromos.gaf_text import parse_unsigned
from dromos.linearize import Window
from dromos.refmap import RefMap

comptime MAGIC = "PGSITES2"
comptime INDEX_FIELDS = 9
comptime MAX_SITE = 500
"""The most reference bases a snarl's window may span before padding."""
comptime WINDOW_PAD = 10
comptime MAX_WINDOW = 1000
"""The widest a merged window grows."""


@always_inline
def _u64(p: ImmPointer[UInt8, _], at: Int) -> Int:
    var v = 0
    for k in range(8):
        v |= Int(p[unsafe_offset=at + k]) << (8 * k)
    return v


def _inflate(lib: OwnedDLHandle, frame: Span[UInt8, _]) raises -> List[UInt8]:
    """A zstd frame whose content size it records, inflated."""
    var size = lib.get_function[UInt64]("ZSTD_getFrameContentSize")(frame.unsafe_ptr(), len(frame))
    if size >= UInt64(0xFFFF_FFFF_FFFF_FFFE):
        raise Error("zstd: a frame without its size")
    var out = List[UInt8](unsafe_uninit_length=Int(size))
    zstd_decompress(lib, frame, out.unsafe_ptr(), Int(size))
    return out^


def _snarl(contig: Contig, map: RefMap, p: ImmPointer[UInt8, _], at: Int, end: Int, mut all: List[Window]):
    """The window of the record `p[at:end]` (a VCF line), if its ID (third column) names two nodes the contig visits
    once and the bases between them are few enough."""
    var field = 0
    var s = at
    while field < 2:
        var hit = external_call["memchr", Int](p.unsafe_offset(s), Int32(ord("\t")), end - s) if s < end else 0
        if hit == 0:
            return
        s = hit - Int(p) + 1
        field += 1
    var e = s
    while e < end and p[unsafe_offset=e] != UInt8(ord("\t")):
        e += 1
    # The ID split at > and <, empty pieces and pieces that are not numbers dropped: exactly two must be left.
    var ends = List[Int]()
    var k = s
    while k <= e:
        var piece = k
        while k < e and p[unsafe_offset=k] != UInt8(ord(">")) and p[unsafe_offset=k] != UInt8(ord("<")):
            k += 1
        if k > piece:
            var v = parse_unsigned(p, piece, k)
            if v >= 0:
                ends.append(v)
        k += 1
    if len(ends) != 2:
        return
    var pa = contig.placed(ends[0])
    var pb = contig.placed(ends[1])
    if pa < 0 or pb < 0:
        return
    var sa = pa >> 1
    var sb = pb >> 1
    var lo: Int
    var hi: Int
    if sa < sb:
        lo = sa + map.length(ends[0])
        hi = sb
    else:
        lo = sb + map.length(ends[1])
        hi = sa
    if hi >= lo and hi - lo <= MAX_SITE:
        all.append(Window(lo, hi))


def load_windows(path: String, contig: Contig, map: RefMap) raises -> List[Window]:
    """The catalog's snarl windows on `contig` (sorted, disjoint)."""
    var file = MappedFile(path)
    var p = file.ptr
    var n = file.size
    var magic = MAGIC.as_bytes()
    var prefix = n >= 7
    var has_magic = n >= 8
    for i in range(8):
        if i < 7:
            prefix = prefix and p[unsafe_offset=i] == magic[i]
        has_magic = has_magic and p[unsafe_offset=i] == magic[i]
    if prefix and not has_magic:
        raise Error(t"{path}: an older .pgsites format; rebuild it with this pgcatalog")
    if not has_magic or n < 24:
        raise Error(t"{path}: not a .pgsites file")
    var lib = load_zstd()
    var header_length = _u64(p, 8)
    var header = _inflate(lib, Span(unsafe_ptr=p.unsafe_offset(16), length=header_length))
    if StringSlice(unsafe_from_utf8=Span(header)).find("##contig=<ID=") < 0:
        raise Error(t"{path}: no ##contig line")
    var index_at = _u64(p, n - 8)
    var count = _u64(p, index_at)
    var all = List[Window]()
    var raw = List[UInt8]()
    for b in range(count):
        var e = index_at + 8 + b * INDEX_FIELDS * 8
        var offset = _u64(p, e + 3 * 8)
        var compressed = _u64(p, e + 4 * 8)
        var raw_length = _u64(p, e + 5 * 8)
        raw.resize(unsafe_uninit_length=raw_length)
        zstd_decompress(lib, Span(unsafe_ptr=p.unsafe_offset(offset), length=compressed), raw.unsafe_ptr(), raw_length)
        var q = raw.unsafe_ptr().as_imm()
        var at = 0
        while at < raw_length:
            var hit = external_call["memchr", Int](q.unsafe_offset(at), Int32(ord("\n")), raw_length - at)
            var end = raw_length if hit == 0 else hit - Int(q)
            if end > at:
                _snarl(contig, map, q, at, end, all)
            at = end + 1
    # The file is read through `p`, untracked: it stays mapped until here.
    _ = file^

    # By start, the widest first; repeats dropped; then windows inside the one before (but for an empty one at its
    # end) dropped.
    def before(a: Window, b: Window) -> Bool:
        return a.lo < b.lo or (a.lo == b.lo and a.hi > b.hi)

    sort(Span(all), before)
    var outer = List[Window]()
    for k in range(len(all)):
        var w = all[k]
        if k > 0 and w == all[k - 1]:
            continue
        if len(outer) > 0:
            var last = outer[len(outer) - 1]
            if w.lo >= last.lo and w.hi <= last.hi and not (w.lo == w.hi and w.lo == last.hi):
                continue
        outer.append(w)
    var length = len(contig.seq)
    var out = List[Window]()
    for w in outer:
        var lo = max(w.lo - WINDOW_PAD, 1)
        var hi = min(w.hi + WINDOW_PAD, max(length - 1, 0))
        var k = len(out)
        if k > 0 and lo <= out[k - 1].hi and max(hi, out[k - 1].hi) - out[k - 1].lo <= MAX_WINDOW:
            out[k - 1].hi = max(out[k - 1].hi, hi)
        else:
            out.append(Window(lo, hi))
    return out^
