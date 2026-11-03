"""The site catalog: one file per graph and reference sample, every contig's records in position blocks of columns,
each column its own zstd frame (docs/catalog-format.md is the format's specification; keep the two in step).

Written by `dromos catalog` (`SitesWriter`, fed `EncodedBlock`s that catalog.mojo encodes), read by `dromos view`
(`view`: pgcatalog view's VCF text, rebuilt from the columns), `dromos catalog-info` (`info`) and `dromos surject
--catalog` (`load_all_windows`: each record whose ID names its snarl's two boundary nodes, both visited once by the
contig, gives the reference bases between them, at most MAX_SITE; windows inside another are dropped, the rest padded
by WINDOW_PAD and those that overlap merged up to MAX_WINDOW, as pgcatalog pileup's `load_windows` makes them).

A reader needs nothing but the file: the node letters the alleles spell are in each block's node table. The header
names the graph the catalog was built from (its fingerprint, tags and node ids), so a catalog is never used with
another graph's alignments.
"""

from std.ffi import OwnedDLHandle, external_call
from std.memory import unsafe_memcpy

from dromos.bgzf import MappedFile, load_libdeflate, load_zstd
from dromos.contig import Contig, Reference
from dromos.linearize import Window
from dromos.output import OutFile, decimal_digits, put_padded, write_all
from dromos.parallel import ParallelTask, parallel_for
from dromos.refmap import RefMap
from dromos.stream import imm_untracked, mut_untracked, untracked

comptime DROMOS_VERSION = "0.1.0"
comptime MAGIC = "DRMSITES"
comptime VERSION = 1
"""The format version this dromos writes and reads; a file of any other is refused."""
comptime PREAMBLE = 64
"""The fixed bytes at the file's start: magic, version, the file's size, where the header and index are, a hash."""
comptime COLUMNS = 13
comptime POS = 0
comptime FLAGS = 1
comptime PAD = 2
comptime BOUND = 3
comptime NWALK = 4
comptime WLEN = 5
comptime WALK = 6
comptime LV = 7
comptime PS = 8
comptime NODEID = 9
comptime NODELEN = 10
comptime SEQ = 11
comptime EXC = 12


def column_name(c: Int) -> String:
    var names: List[String] = [
        "POS", "FLAGS", "PAD", "BOUND", "NWALK", "WLEN", "WALK", "LV", "PS", "NODEID", "NODELEN", "SEQ", "EXC"
    ]
    return names[c] if c >= 0 and c < len(names) else String(t"column {c}")


comptime FLAG_REVERSED = 1
"""The ID is the first walk's ends flipped and swapped: the snarl was found from its end node."""
comptime FLAG_PAD = 2
"""Every allele starts with the record's padding base (PAD)."""
comptime FLAG_PS = 4
"""The record has a PS (a parent snarl with a record)."""
comptime INDEX_ARRAYS = 9
"""Per-block arrays of the index before the column sizes: first_pos, last_pos, max_end, reach, offset, records,
min_allele, max_walk, crc32."""
comptime BLOCK_BP = 100_000
"""A block holds the records whose POS falls in one BLOCK_BP window..."""
comptime BLOCK_BYTES = 4 << 20
"""...up to BLOCK_BYTES of their VCF text..."""
comptime OWN_BLOCK_BYTES = 1 << 20
"""...and a record whose VCF line is longer than OWN_BLOCK_BYTES has a block of its own."""
comptime LONG_BYTES = 1 << 20
"""A column longer than this is compressed with long-distance matching and a window as long as the column."""
comptime MT_BYTES = 16 << 20
"""A column longer than this is compressed with zstd's own threads (whose output does not depend on how many)."""
comptime DEFAULT_LEVEL = 9
comptime META_LEVEL = 9
"""The header's and the index's zstd level: 19 makes the HPRC index 2% smaller in 17 times as long (7.9 s)."""
comptime MAX_SITE = 500
"""The most reference bases a snarl's window may span before padding."""
comptime WINDOW_PAD = 10
comptime MAX_WINDOW = 1000
"""The widest a merged window grows."""
comptime VIEW_BYTES = 256 << 20
"""The most column bytes a step of `view` decodes at once (blocks in parallel)."""
comptime STREAM_TEXT = 64 << 20
"""A block whose columns hold more than this is printed on the calling thread, its text written as it is made."""


# Letters as alleles are written: upper case, or reverse-complemented.


@always_inline
def upper(c: UInt8) -> UInt8:
    return c - 32 if c >= UInt8(ord("a")) and c <= UInt8(ord("z")) else c


@always_inline
def complement(c: UInt8) -> UInt8:
    """The crate's COMPLEMENT: A, C, G, T in either case to their complement in capitals, any other letter to N."""
    var u = upper(c)
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
            dst[unsafe_offset=i] = complement(s[unsafe_offset=n - 1 - i])
            i += 1
    else:
        while i + _W <= n:
            dst.unsafe_store(i, _upper_lanes(s.unsafe_load[width=_W](i)))
            i += _W
        while i < n:
            dst[unsafe_offset=i] = upper(s[unsafe_offset=i])
            i += 1


def vcf_header(contigs: List[Tuple[String, Int]]) -> String:
    """pgcatalog's VCF header, with a ##contig line for each (CHROM, length)."""
    var out = String(
        "##fileformat=VCFv4.2\n",
        '##INFO=<ID=AT,Number=R,Type=String,Description="Allele Traversal as path in graph">\n',
        '##INFO=<ID=LV,Number=1,Type=Integer,Description="Level in the snarl tree (0=top level)">\n',
        '##INFO=<ID=PS,Number=1,Type=String,Description="ID of variant corresponding to parent snarl">\n',
    )
    for c in contigs:
        out += String(t"##contig=<ID={c[0]},length={c[1]}>\n")
    out += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    return out^


# Integers in columns.


@always_inline
def put_varint(mut out: List[UInt8], value: UInt64):
    """LEB128: seven bits a byte, low first, the high bit set on every byte but the last."""
    var v = value
    while v >= 0x80:
        out.append(UInt8(v & 0x7F) | 0x80)
        v >>= 7
    out.append(UInt8(v))


@always_inline
def zigzag(v: Int) -> UInt64:
    return UInt64(v << 1) ^ UInt64(v >> 63)


@always_inline
def unzigzag(v: UInt64) -> Int:
    return Int(v >> 1) ^ -Int(v & 1)


def _damaged[W: Writable](what: W) -> Error:
    return Error(t"a damaged catalog block: {what}")


struct Cursor[mut: Bool, //, origin: Origin[mut=mut]](Movable):
    """A column read from its start: varints and bytes, never past its end. It keeps the column alive."""

    var p: Pointer[UInt8, Self.origin]
    var at: Int
    var end: Int

    def __init__(out self, ref[Self.origin] column: List[UInt8]):
        self.p = column.unsafe_ptr()
        self.at = 0
        self.end = len(column)

    def __init__(out self, p: Pointer[UInt8, Self.origin], length: Int):
        self.p = p
        self.at = 0
        self.end = length

    @always_inline
    def varint(mut self) raises -> UInt64:
        var v = UInt64(0)
        var shift = UInt64(0)
        while True:
            if self.at >= self.end or shift > 63:
                raise _damaged("a column ends early")
            var c = self.p[unsafe_offset=self.at]
            self.at += 1
            v |= UInt64(c & 0x7F) << shift
            if c < 0x80:
                return v
            shift += 7

    @always_inline
    def byte(mut self) raises -> UInt8:
        if self.at >= self.end:
            raise _damaged("a column ends early")
        var c = self.p[unsafe_offset=self.at]
        self.at += 1
        return c


@always_inline
def put_u64(mut out: List[UInt8], value: Int):
    for i in range(8):
        out.append(UInt8((value >> (8 * i)) & 0xFF))


def put_str(mut out: List[UInt8], s: String):
    put_u64(out, s.byte_length())
    out.extend(s.as_bytes())


@always_inline
def u64_at(p: ImmPointer[UInt8, _], at: Int) -> Int:
    var v = 0
    for k in range(8):
        v |= Int(p[unsafe_offset=at + k]) << (8 * k)
    return v


struct Reader(Movable):
    """Fields of the header or the index, read in order, never past the end."""

    var data: List[UInt8]
    var at: Int
    var what: String

    def __init__(out self, var data: List[UInt8], what: String):
        self.data = data^
        self.at = 0
        self.what = what

    def u64(mut self) raises -> Int:
        if self.at + 8 > len(self.data):
            raise Error(t"{self.what}: damaged (it ends early)")
        var v = u64_at(self.data.unsafe_ptr(), self.at)
        self.at += 8
        if v < 0:
            raise Error(t"{self.what}: damaged (a number out of range)")
        return v

    def bytes(mut self, n: Int) raises -> List[UInt8]:
        if n < 0 or self.at + n > len(self.data):
            raise Error(t"{self.what}: damaged (it ends early)")
        var out = List[UInt8](Span(self.data)[self.at : self.at + n])
        self.at += n
        return out^

    def string(mut self) raises -> String:
        return String(unsafe_from_utf8=self.bytes(self.u64()))


# zstd, as the writer sets it up for each frame.

comptime _C_LEVEL: Int32 = 100
comptime _C_WINDOW_LOG: Int32 = 101
comptime _C_LDM: Int32 = 160
comptime _C_CONTENT_SIZE: Int32 = 200
comptime _C_CHECKSUM: Int32 = 201
comptime _C_WORKERS: Int32 = 400
comptime _RESET_ALL: Int32 = 3
comptime _D_WINDOW_LOG_MAX: Int32 = 100


def window_log(length: Int) -> Int:
    """The window a long column of `length` bytes is compressed with: as long as the column, from 2^20 to 2^31."""
    var x = max(length, 2) - 1
    var bits = 0
    while x > 0:
        bits += 1
        x >>= 1
    return min(max(bits, 20), 31)


struct Compressor(Movable):
    """A zstd compression context, set up afresh for each frame (a context made for every frame costs as much as
    compressing a small one)."""

    var lib: OwnedDLHandle
    var cctx: Int

    def __init__(out self) raises:
        self.lib = load_zstd()
        self.cctx = self.lib.get_function[Int]("ZSTD_createCCtx")()
        if self.cctx == 0:
            raise Error("zstd: no compression context")

    def __deinit__(deinit self):
        try:
            _ = self.lib.get_function[Int]("ZSTD_freeCCtx")(self.cctx)
        except:
            pass

    def _check(self, code: Int) raises:
        if self.lib.get_function[UInt32]("ZSTD_isError")(code) != 0:
            var name = self.lib.get_function[Int]("ZSTD_getErrorName")(code)
            var length = external_call["strlen", Int](name)
            var text = String(
                StringSlice(
                    unsafe_from_utf8=Span(
                        unsafe_ptr=ImmPointer[UInt8, ImmUntrackedOrigin](unsafe_from_address=name), length=length
                    )
                )
            )
            raise Error(t"zstd: {text}")

    def _set(self, param: Int32, value: Int) raises:
        self._check(self.lib.get_function[Int]("ZSTD_CCtx_setParameter")(self.cctx, param, Int32(value)))

    def append(
        self, data: Span[UInt8, _], level: Int, threads: Int, mut out: List[UInt8], checksum: Bool = False
    ) raises -> Int:
        """`data` as one frame appended to `out`: at `level`, its content size recorded, and for a column over
        LONG_BYTES long-distance matching with a window as long as it (zstd's own threads over MT_BYTES). Returns
        the frame's size."""
        var n = len(data)
        self._check(self.lib.get_function[Int]("ZSTD_CCtx_reset")(self.cctx, _RESET_ALL))
        self._set(_C_LEVEL, level)
        self._set(_C_CONTENT_SIZE, 1)
        if checksum:
            self._set(_C_CHECKSUM, 1)
        if n > LONG_BYTES:
            self._set(_C_LDM, 1)
            self._set(_C_WINDOW_LOG, window_log(n))
            if n > MT_BYTES:
                self._set(_C_WORKERS, max(threads, 1))
        var bound = self.lib.get_function[Int]("ZSTD_compressBound")(n)
        var at = len(out)
        out.resize(unsafe_uninit_length=at + bound)
        var size = self.lib.get_function[Int]("ZSTD_compress2")(
            self.cctx, out.unsafe_ptr().unsafe_offset(at), bound, data.unsafe_ptr(), n
        )
        if self.lib.get_function[UInt32]("ZSTD_isError")(size) != 0:
            out.resize(unsafe_uninit_length=at)
            self._check(size)
        out.resize(unsafe_uninit_length=at + size)
        return size

    def column(self, data: Span[UInt8, _], level: Int, threads: Int, mut out: List[UInt8]) raises -> Int:
        """A column appended to `out` as a zstd frame, or as it is when that is no smaller (a reader tells them
        apart by the sizes: compressed equal to raw is stored). Returns the bytes appended; none for an empty
        column."""
        var n = len(data)
        if n == 0:
            return 0
        var at = len(out)
        var size = self.append(data, level, threads, out)
        if size < n:
            return size
        out.resize(unsafe_uninit_length=at)
        out.extend(data)
        return n


struct Decompressor(Movable):
    """A zstd decompression context (windows up to 2^31 allowed), and libdeflate for blocks' CRCs."""

    var lib: OwnedDLHandle
    var deflate: OwnedDLHandle
    var dctx: Int

    def __init__(out self) raises:
        self.lib = load_zstd()
        self.deflate = load_libdeflate()
        self.dctx = self.lib.get_function[Int]("ZSTD_createDCtx")()
        if self.dctx == 0:
            raise Error("zstd: no decompression context")
        _ = self.lib.get_function[Int]("ZSTD_DCtx_setParameter")(self.dctx, _D_WINDOW_LOG_MAX, Int32(31))

    def __deinit__(deinit self):
        try:
            _ = self.lib.get_function[Int]("ZSTD_freeDCtx")(self.dctx)
        except:
            pass

    def into(self, frame: Span[UInt8, _], dst: MutPointer[UInt8, _], length: Int) -> Bool:
        """Inflates `frame` into `dst`, which it must fill exactly."""
        try:
            var n = self.lib.get_function[Int]("ZSTD_decompressDCtx")(
                self.dctx, dst, length, frame.unsafe_ptr(), len(frame)
            )
            return self.lib.get_function[UInt32]("ZSTD_isError")(n) == 0 and n == length
        except:
            return False

    def frame(self, frame: Span[UInt8, _], what: String) raises -> List[UInt8]:
        """A frame that records its content size, inflated (its checksum checked if it has one)."""
        var size = self.lib.get_function[UInt64]("ZSTD_getFrameContentSize")(frame.unsafe_ptr(), len(frame))
        if size >= UInt64(0xFFFF_FFFF_FFFF_FFFE) or size > UInt64(1 << 40):
            raise Error(t"{what}: damaged")
        var out = List[UInt8](unsafe_uninit_length=Int(size))
        if not self.into(frame, out.unsafe_ptr(), Int(size)):
            raise Error(t"{what}: damaged")
        return out^


def crc32(lib: OwnedDLHandle, data: Span[UInt8, _]) raises -> UInt32:
    """The CRC-32 zlib's crc32 computes (libdeflate's)."""
    return lib.get_function[UInt32]("libdeflate_crc32")(UInt32(0), data.unsafe_ptr(), len(data))


def fnv1a(data: Span[UInt8, _]) -> UInt64:
    var h = UInt64(0xCBF29CE484222325)
    for b in data:
        h = (h ^ UInt64(b)) * 0x100000001B3
    return h


# MD5 (RFC 1321), for each contig's reference sequence: SAM's M5.

comptime _MD5_K = SIMD[DType.uint32, 64](
    0xD76AA478,
    0xE8C7B756,
    0x242070DB,
    0xC1BDCEEE,
    0xF57C0FAF,
    0x4787C62A,
    0xA8304613,
    0xFD469501,
    0x698098D8,
    0x8B44F7AF,
    0xFFFF5BB1,
    0x895CD7BE,
    0x6B901122,
    0xFD987193,
    0xA679438E,
    0x49B40821,
    0xF61E2562,
    0xC040B340,
    0x265E5A51,
    0xE9B6C7AA,
    0xD62F105D,
    0x02441453,
    0xD8A1E681,
    0xE7D3FBC8,
    0x21E1CDE6,
    0xC33707D6,
    0xF4D50D87,
    0x455A14ED,
    0xA9E3E905,
    0xFCEFA3F8,
    0x676F02D9,
    0x8D2A4C8A,
    0xFFFA3942,
    0x8771F681,
    0x6D9D6122,
    0xFDE5380C,
    0xA4BEEA44,
    0x4BDECFA9,
    0xF6BB4B60,
    0xBEBFBC70,
    0x289B7EC6,
    0xEAA127FA,
    0xD4EF3085,
    0x04881D05,
    0xD9D4D039,
    0xE6DB99E5,
    0x1FA27CF8,
    0xC4AC5665,
    0xF4292244,
    0x432AFF97,
    0xAB9423A7,
    0xFC93A039,
    0x655B59C3,
    0x8F0CCC92,
    0xFFEFF47D,
    0x85845DD1,
    0x6FA87E4F,
    0xFE2CE6E0,
    0xA3014314,
    0x4E0811A1,
    0xF7537E82,
    0xBD3AF235,
    0x2AD7D2BB,
    0xEB86D391,
)
comptime _MD5_S = SIMD[DType.uint32, 16](7, 12, 17, 22, 5, 9, 14, 20, 4, 11, 16, 23, 6, 10, 15, 21)


struct Md5(Movable):
    """MD5, fed in pieces."""

    var state: SIMD[DType.uint32, 4]
    var buffer: List[UInt8]
    var length: UInt64

    def __init__(out self):
        self.state = SIMD[DType.uint32, 4](0x67452301, 0xEFCDAB89, 0x98BADCFE, 0x10325476)
        self.buffer = List[UInt8](capacity=64)
        self.length = 0

    @staticmethod
    @always_inline
    def _rotl(x: UInt32, c: UInt32) -> UInt32:
        return (x << c) | (x >> (32 - c))

    def _block(mut self, p: ImmPointer[UInt8, _]):
        var m = SIMD[DType.uint32, 16]()
        comptime for i in range(16):
            m[i] = (
                UInt32(p[unsafe_offset=4 * i])
                | (UInt32(p[unsafe_offset=4 * i + 1]) << 8)
                | (UInt32(p[unsafe_offset=4 * i + 2]) << 16)
                | (UInt32(p[unsafe_offset=4 * i + 3]) << 24)
            )
        var a = self.state[0]
        var b = self.state[1]
        var c = self.state[2]
        var d = self.state[3]
        comptime for i in range(64):
            var f: UInt32
            comptime if i < 16:
                f = (b & c) | (~b & d)
            elif i < 32:
                f = (d & b) | (~d & c)
            elif i < 48:
                f = b ^ c ^ d
            else:
                f = c ^ (b | ~d)
            comptime g = i if i < 16 else ((5 * i + 1) & 15 if i < 32 else ((3 * i + 5) & 15 if i < 48 else (7 * i) & 15))
            f = f + a + _MD5_K[i] + m[g]
            a = d
            d = c
            c = b
            b = b + Self._rotl(f, _MD5_S[(i >> 4) * 4 + (i & 3)])
        self.state += SIMD[DType.uint32, 4](a, b, c, d)

    def update(mut self, data: Span[UInt8, _]):
        var n = len(data)
        var p = data.unsafe_ptr()
        var i = 0
        self.length += UInt64(n)
        if len(self.buffer) > 0:
            while i < n and len(self.buffer) < 64:
                self.buffer.append(p[unsafe_offset=i])
                i += 1
            if len(self.buffer) == 64:
                var full = self.buffer^
                self.buffer = List[UInt8](capacity=64)
                self._block(full.unsafe_ptr())
        while i + 64 <= n:
            self._block(p.unsafe_offset(i))
            i += 64
        while i < n:
            self.buffer.append(p[unsafe_offset=i])
            i += 1

    def digest(mut self) -> List[UInt8]:
        var bits = self.length * 8
        var pad = List[UInt8]()
        pad.append(0x80)
        while (len(self.buffer) + len(pad)) % 64 != 56:
            pad.append(0)
        for k in range(8):
            pad.append(UInt8((bits >> UInt64(8 * k)) & 0xFF))
        var total = self.length
        self.update(Span(pad))
        self.length = total
        var out = List[UInt8]()
        for k in range(4):
            for j in range(4):
                out.append(UInt8((self.state[k] >> UInt32(8 * j)) & 0xFF))
        return out^


def hex_string(data: Span[UInt8, _]) -> String:
    var digits = "0123456789abcdef".as_bytes()
    var out = List[UInt8]()
    for b in data:
        out.append(digits[Int(b >> 4)])
        out.append(digits[Int(b & 15)])
    return String(unsafe_from_utf8=out)


struct LetterPacker(Movable):
    """Node letters into the SEQ and EXC columns: two bits each (A C G T as 0 to 3, four a byte from the low bits),
    and each letter that is none of those (after upper-casing) as an exception: its index as a difference from the
    exception before (the first from 0), then the letter, its two bits 0."""

    var letters: Int
    var last_exc: Int
    var acc: UInt8
    var bits: Int

    def __init__(out self):
        self.letters = 0
        self.last_exc = 0
        self.acc = 0
        self.bits = 0

    def add(mut self, s: Span[UInt8, _], mut seq: List[UInt8], mut exc: List[UInt8]):
        for c in s:
            var u = upper(c)
            var code: UInt8
            if u == UInt8(ord("A")):
                code = 0
            elif u == UInt8(ord("C")):
                code = 1
            elif u == UInt8(ord("G")):
                code = 2
            elif u == UInt8(ord("T")):
                code = 3
            else:
                code = 0
                put_varint(exc, UInt64(self.letters - self.last_exc))
                exc.append(u)
                self.last_exc = self.letters
            self.acc |= code << UInt8(self.bits)
            self.bits += 2
            if self.bits == 8:
                seq.append(self.acc)
                self.acc = 0
                self.bits = 0
            self.letters += 1

    def finish(mut self, mut seq: List[UInt8]):
        if self.bits > 0:
            seq.append(self.acc)
            self.acc = 0
            self.bits = 0


# Writing.


struct EncodedBlock(Movable):
    """A block as stored: its columns' frames one after another, their sizes, and what the index says of it."""

    var first_pos: Int
    var last_pos: Int
    var max_end: Int
    var records: Int
    var min_allele: Int
    var max_walk: Int
    var sizes: List[Int]
    """Per column, its compressed then its raw size (equal: stored as it is; 0: empty)."""
    var data: List[UInt8]

    def __init__(out self):
        self.first_pos = 0
        self.last_pos = 0
        self.max_end = 0
        self.records = 0
        self.min_allele = 0
        self.max_walk = 0
        self.sizes = List[Int](length=2 * COLUMNS, fill=0)
        self.data = List[UInt8]()


@fieldwise_init
struct ContigEntry(Copyable, Movable):
    """A contig as the header lists it."""

    var name: String
    """chr20."""
    var chrom: String
    """CHM13#0#chr20: the records' CHROM."""
    var length: Int
    var md5: List[UInt8]
    """MD5 of its reference sequence (16 bytes; SAM's M5)."""
    var regions: List[Tuple[Int, Int]]
    """The regions built (1-based, inclusive); none: the whole contig."""
    var first_block: Int
    var blocks: Int
    var records: Int


struct SitesWriter(Movable):
    """Writes `<path>.tmp` and renames it to `path` once complete, so a catalog is never seen half written. Contigs
    are built in any order (any thread, one `complete` at a time, which a lock keeps) and written in the header's
    order as each one and those before it are complete."""

    var path: String
    var file: OutFile
    var offset: Int
    var props: List[Tuple[String, String]]
    var contigs: List[ContigEntry]
    var staged: List[List[EncodedBlock]]
    var ready: List[Bool]
    var next: Int
    var arrays: List[List[Int]]
    """The index's per-block arrays (INDEX_ARRAYS), in file order."""
    var sizes: List[Int]
    """Per block, per column: compressed, raw."""
    var crc_lib: OwnedDLHandle

    def __init__(out self, path: String, var props: List[Tuple[String, String]], var contigs: List[ContigEntry]) raises:
        self.path = path
        self.file = OutFile(path + ".tmp")
        self.offset = PREAMBLE
        self.props = props^
        self.contigs = contigs^
        self.staged = List[List[EncodedBlock]]()
        for _ in range(len(self.contigs)):
            self.staged.append(List[EncodedBlock]())
        self.ready = List[Bool](length=len(self.contigs), fill=False)
        self.next = 0
        self.arrays = List[List[Int]]()
        for _ in range(INDEX_ARRAYS):
            self.arrays.append(List[Int]())
        self.sizes = List[Int]()
        self.crc_lib = load_libdeflate()
        # The preamble's room, filled in by `finish`.
        self.file.put(Span(List[UInt8](length=PREAMBLE, fill=0)), 0)

    def complete(
        mut self, k: Int, var blocks: List[EncodedBlock], md5: List[UInt8], records: Int, length: Int
    ) raises:
        """Contig `k`'s blocks, in position order, its reference sequence's MD5, records and length; written now if
        every contig before it is."""
        self.staged[k] = blocks^
        self.contigs[k].length = length
        self.contigs[k].md5 = md5.copy()
        self.contigs[k].records = records
        self.ready[k] = True
        while self.next < len(self.contigs) and self.ready[self.next]:
            self._write(self.next)
            self.next += 1

    def _write(mut self, k: Int) raises:
        ref c = self.contigs[k]
        c.first_block = len(self.arrays[0])
        c.blocks = len(self.staged[k])
        var reach = 0
        var out = List[UInt8]()
        for i in range(len(self.staged[k])):
            ref b = self.staged[k][i]
            reach = max(reach, b.max_end)
            var crc = crc32(self.crc_lib, Span(b.data))
            var fields: List[Int] = [
                b.first_pos,
                b.last_pos,
                b.max_end,
                reach,
                self.offset + len(out),
                b.records,
                b.min_allele,
                b.max_walk,
                Int(crc),
            ]
            for i in range(INDEX_ARRAYS):
                self.arrays[i].append(fields[i])
            self.sizes.extend(Span(b.sizes))
            out.extend(Span(b.data))
            if len(out) >= 64 << 20:
                self.file.put(Span(out), self.offset)
                self.offset += len(out)
                out.clear()
        self.file.put(Span(out), self.offset)
        self.offset += len(out)
        self.staged[k] = List[EncodedBlock]()

    def finish(mut self) raises -> Int:
        """Writes the header and the index after the blocks, then the preamble; syncs the file and renames it.
        Returns the file's size."""
        if self.next != len(self.contigs):
            raise Error(t"{self.path}: contigs left unwritten")
        var c = Compressor()
        var header = List[UInt8]()
        put_u64(header, VERSION)
        put_u64(header, len(self.props))
        for p in self.props:
            put_str(header, p[0])
            put_str(header, p[1])
        put_u64(header, len(self.contigs))
        for ref e in self.contigs:
            put_str(header, e.name)
            put_str(header, e.chrom)
            put_u64(header, e.length)
            header.extend(Span(e.md5))
            put_u64(header, e.first_block)
            put_u64(header, e.blocks)
            put_u64(header, e.records)
            put_u64(header, len(e.regions))
            for r in e.regions:
                put_u64(header, r[0])
                put_u64(header, r[1])
        var index = List[UInt8]()
        var blocks = len(self.arrays[0])
        put_u64(index, blocks)
        put_u64(index, COLUMNS)
        for a in self.arrays:
            for v in a:
                put_u64(index, v)
        for which in range(2):
            for b in range(blocks):
                for col in range(COLUMNS):
                    put_u64(index, self.sizes[b * 2 * COLUMNS + 2 * col + which])
        var tail = List[UInt8]()
        var header_size = c.append(Span(header), META_LEVEL, 1, tail, checksum=True)
        var index_size = c.append(Span(index), META_LEVEL, 1, tail, checksum=True)
        var header_at = self.offset
        self.file.put(Span(tail), self.offset)
        self.offset += len(tail)
        var pre = List[UInt8](MAGIC.as_bytes())
        var v = List[UInt8]()
        put_u64(v, VERSION)
        pre.extend(Span(v)[:4])
        pre.extend(Span(List[UInt8](length=4, fill=0)))
        put_u64(pre, self.offset)
        put_u64(pre, header_at)
        put_u64(pre, header_size)
        put_u64(pre, header_at + header_size)
        put_u64(pre, index_size)
        put_u64(pre, Int(fnv1a(Span(pre))))
        self.file.put(Span(pre), 0)
        if external_call["fsync", Int32](self.file.fd) != 0:
            raise Error(t"cannot write {self.path}.tmp")
        self.file.close()
        var a = List[UInt8]((self.path + ".tmp").as_bytes())
        a.append(0)
        var b = List[UInt8](self.path.as_bytes())
        b.append(0)
        if external_call["rename", Int32](a.unsafe_ptr(), b.unsafe_ptr()) != 0:
            raise Error(t"cannot rename {self.path}.tmp to {self.path}")
        return self.offset


# Reading.


def _not_a_catalog(path: String, p: ImmPointer[UInt8, _], n: Int) -> String:
    var pgsites = "PGSITES".as_bytes()
    var old = n >= 8
    for i in range(min(7, n)):
        old = old and p[unsafe_offset=i] == pgsites[i]
    if old:
        return String(
            t"{path}: a .pgsites file (pgcatalog's catalog of one contig); dromos reads its own site catalog: make one"
            t" with dromos catalog"
        )
    return String(t"{path}: not a dromos site catalog")


def _map(path: String) raises -> MappedFile:
    try:
        return MappedFile(path, ahead=False)
    except e:
        if String(e).endswith(": empty"):
            raise Error(t"{path}: not a dromos site catalog (an empty file)")
        raise e^


struct _Parsed(Movable):
    var props: List[Tuple[String, String]]
    var contigs: List[ContigEntry]
    var arrays: List[List[Int]]
    var sizes: List[Int]

    def __init__(out self):
        self.props = List[Tuple[String, String]]()
        self.contigs = List[ContigEntry]()
        self.arrays = List[List[Int]]()
        self.sizes = List[Int]()


def _read(path: String, file: MappedFile) raises -> _Parsed:
    """A catalog's preamble checked, its header and index read: `file` is borrowed for the whole call, so it stays
    mapped on every path, raising ones included."""
    var p = file.ptr
    var n = file.size
    var magic = MAGIC.as_bytes()
    var is_magic = n >= 8
    for i in range(min(8, n)):
        is_magic = is_magic and p[unsafe_offset=i] == magic[i]
    if not is_magic:
        raise Error(_not_a_catalog(path, p, n))
    if n < PREAMBLE:
        raise Error(t"{path}: truncated (a catalog's first {PREAMBLE} bytes are its preamble; it has {n})")
    var version = 0
    for k in range(4):
        version |= Int(p[unsafe_offset=8 + k]) << (8 * k)
    if version != VERSION:
        raise Error(
            t"{path}: a site catalog of format version {version}; this dromos reads version {VERSION}: rebuild it"
            t" with this dromos (dromos catalog) or read it with the dromos that wrote it"
        )
    if UInt64(u64_at(p, 56)) != fnv1a(Span(unsafe_ptr=p, length=56)):
        raise Error(t"{path}: damaged (its preamble's hash differs)")
    var size = u64_at(p, 16)
    if size != n:
        if n < size:
            raise Error(t"{path}: truncated ({n} of its {size} bytes)")
        raise Error(t"{path}: {n} bytes where its preamble says {size}: damaged")
    var header_at = u64_at(p, 24)
    var header_size = u64_at(p, 32)
    var index_at = u64_at(p, 40)
    var index_size = u64_at(p, 48)
    for r in [(header_at, header_size), (index_at, index_size)]:
        if r[0] < PREAMBLE or r[1] <= 0 or r[0] + r[1] > n:
            raise Error(t"{path}: damaged (its preamble points outside the file)")
    var parsed = _Parsed()
    var z = Decompressor()
    var h = Reader(
        z.frame(Span(unsafe_ptr=p.unsafe_offset(header_at), length=header_size), path + ": its header"),
        path + ": its header",
    )
    if h.u64() != VERSION:
        raise Error(t"{path}: damaged (its header's version differs from its preamble's)")
    parsed.props = List[Tuple[String, String]]()
    for _ in range(h.u64()):
        var key = h.string()
        parsed.props.append((key, h.string()))
    parsed.contigs = List[ContigEntry]()
    for _ in range(h.u64()):
        var name = h.string()
        var chrom = h.string()
        var length = h.u64()
        var md5 = h.bytes(16)
        var first = h.u64()
        var blocks = h.u64()
        var records = h.u64()
        var regions = List[Tuple[Int, Int]]()
        for _ in range(h.u64()):
            var a = h.u64()
            regions.append((a, h.u64()))
        parsed.contigs.append(ContigEntry(name, chrom, length, md5^, regions^, first, blocks, records))
    var x = Reader(
        z.frame(Span(unsafe_ptr=p.unsafe_offset(index_at), length=index_size), path + ": its index"),
        path + ": its index",
    )
    var count = x.u64()
    if x.u64() != COLUMNS or len(x.data) != 16 + 8 * count * (INDEX_ARRAYS + 2 * COLUMNS):
        raise Error(t"{path}: damaged (its index has the wrong size)")
    parsed.arrays = List[List[Int]]()
    for _ in range(INDEX_ARRAYS):
        var a = List[Int](capacity=count)
        for _ in range(count):
            a.append(x.u64())
        parsed.arrays.append(a^)
    parsed.sizes = List[Int](length=2 * COLUMNS * count, fill=0)
    for which in range(2):
        for b in range(count):
            for col in range(COLUMNS):
                parsed.sizes[b * 2 * COLUMNS + 2 * col + which] = x.u64()
    var blocks_seen = 0
    for ref c in parsed.contigs:
        if c.first_block != blocks_seen or c.first_block + c.blocks > count:
            raise Error(t"{path}: damaged (its contigs' blocks)")
        blocks_seen += c.blocks
    if blocks_seen != count:
        raise Error(t"{path}: damaged (its contigs' blocks)")
    for b in range(count):
        var total = 0
        for col in range(COLUMNS):
            var cs = parsed.sizes[b * 2 * COLUMNS + 2 * col]
            var rs = parsed.sizes[b * 2 * COLUMNS + 2 * col + 1]
            if (cs == 0) != (rs == 0) or cs > rs:
                raise Error(t"{path}: damaged (block {b}'s column sizes)")
            total += cs
        var at = parsed.arrays[4][b]
        if at < PREAMBLE or at + total > header_at:
            raise Error(t"{path}: damaged (block {b} lies outside the blocks)")
    return parsed^


struct Sites(Movable):
    """A catalog mapped into memory: its preamble checked, its header and index read."""

    var path: String
    var file: MappedFile
    var props: List[Tuple[String, String]]
    var contigs: List[ContigEntry]
    var arrays: List[List[Int]]
    var sizes: List[Int]
    """Per block, per column: compressed, raw."""

    def __init__(out self, path: String) raises:
        self.path = path
        var file = _map(path)
        var parsed = _read(path, file)
        self.file = file^
        self.props = List[Tuple[String, String]]()
        self.contigs = List[ContigEntry]()
        self.arrays = List[List[Int]]()
        self.sizes = List[Int]()
        swap(self.props, parsed.props)
        swap(self.contigs, parsed.contigs)
        swap(self.arrays, parsed.arrays)
        swap(self.sizes, parsed.sizes)

    @always_inline
    def blocks(self) -> Int:
        return len(self.arrays[0])

    @always_inline
    def first_pos(self, b: Int) -> Int:
        return self.arrays[0][b]

    @always_inline
    def max_end(self, b: Int) -> Int:
        return self.arrays[2][b]

    @always_inline
    def reach(self, b: Int) -> Int:
        return self.arrays[3][b]

    @always_inline
    def offset(self, b: Int) -> Int:
        return self.arrays[4][b]

    @always_inline
    def records(self, b: Int) -> Int:
        return self.arrays[5][b]

    @always_inline
    def crc(self, b: Int) -> Int:
        return self.arrays[8][b]

    def compressed(self, b: Int) -> Int:
        var n = 0
        for col in range(COLUMNS):
            n += self.sizes[b * 2 * COLUMNS + 2 * col]
        return n

    def raw(self, b: Int) -> Int:
        var n = 0
        for col in range(COLUMNS):
            n += self.sizes[b * 2 * COLUMNS + 2 * col + 1]
        return n

    def prop(self, key: String) -> Optional[String]:
        for p in self.props:
            if p[0] == key:
                return p[1]
        return None

    def contig(self, name: String) -> Optional[Int]:
        """The contig `name` names: chr20, or CHM13#0#chr20 (of the reference sample alone)."""
        for k in range(len(self.contigs)):
            if same_contig(name, self.contigs[k].chrom):
                return k
        return None

    def overlapping(self, k: Int, beg: Int, end: Int) -> Tuple[Int, Int]:
        """Contig `k`'s blocks that can hold records overlapping [beg, end]: those from the first whose reach (the
        largest max_end so far) is at least `beg` to the last whose first_pos is at most `end`, by binary search;
        each must still have max_end >= beg."""
        var lo = self.contigs[k].first_block
        var hi = lo + self.contigs[k].blocks
        var a = lo
        var b = hi
        while a < b:
            var m = (a + b) // 2
            if self.reach(m) < beg:
                a = m + 1
            else:
                b = m
        var first = a
        a = first
        b = hi
        while a < b:
            var m = (a + b) // 2
            if self.first_pos(m) <= end:
                a = m + 1
            else:
                b = m
        return (first, a)

    def columns(self, b: Int, wanted: Int, z: Decompressor, mut out: List[List[UInt8]], check: Bool = True) raises:
        """Block `b`'s columns in the `wanted` mask (bit c for column c) into `out` (COLUMNS lists); the others
        left empty. With `check`, the block's CRC is checked first."""
        var base = self.file.ptr.unsafe_offset(self.offset(b))
        if check:
            if Int(crc32(z.deflate, Span(unsafe_ptr=base, length=self.compressed(b)))) != self.crc(b):
                raise Error(t"{self.path}: damaged (block {b}'s CRC differs)")
        var at = 0
        for col in range(COLUMNS):
            var cs = self.sizes[b * 2 * COLUMNS + 2 * col]
            var rs = self.sizes[b * 2 * COLUMNS + 2 * col + 1]
            ref dst = out[col]
            dst.clear()
            if (wanted >> col) & 1 == 1 and rs > 0:
                dst.resize(unsafe_uninit_length=rs)
                if cs == rs:
                    unsafe_memcpy(dest=dst.unsafe_ptr(), src=base.unsafe_offset(at), count=rs)
                elif not z.into(Span(unsafe_ptr=base.unsafe_offset(at), length=cs), dst.unsafe_ptr(), rs):
                    raise Error(t"{self.path}: damaged (block {b}'s {column_name(col)} column)")
            at += cs


def same_contig(a: String, b: String) -> Bool:
    """`chr20` names the contig of any sample's path; `CHM13#0#chr20` only CHM13's."""

    def split(s: String) -> Tuple[String, String, Bool]:
        var i = s.find("#")
        if i < 0:
            return (String(), s, False)
        var j = s.rfind("#")
        return (String(s[byte=:i]), String(s[byte = j + 1 :]), True)

    var x = split(a)
    var y = split(b)
    return x[1] == y[1] and (not x[2] or not y[2] or x[0] == y[0])


# Blocks decoded: the node table, then each record.


struct NodeTable(Movable):
    """A block's nodes: ids in order, and each one's letters (capitals; other characters as stored)."""

    var ids: List[UInt64]
    var starts: List[Int]
    """Node i's letters are letters[starts[i] ..< starts[i + 1]]."""
    var letters: List[UInt8]

    def __init__(out self):
        self.ids = List[UInt64]()
        self.starts = List[Int]()
        self.letters = List[UInt8]()

    def read(mut self, cols: List[List[UInt8]], with_letters: Bool) raises:
        """From the NODEID and NODELEN columns, and with `with_letters` SEQ and EXC."""
        self.ids.clear()
        self.starts.clear()
        self.starts.append(0)
        var ids = Cursor(cols[NODEID])
        var lens = Cursor(cols[NODELEN])
        var id = UInt64(0)
        var total = 0
        while ids.at < ids.end:
            id += ids.varint()
            self.ids.append(id)
            total += Int(lens.varint())
            self.starts.append(total)
        if lens.at != lens.end:
            raise _damaged("node lengths left over")
        if not with_letters:
            return
        ref seq = cols[SEQ]
        if len(seq) != (total + 3) // 4:
            raise _damaged("the letters' column has the wrong size")
        self.letters.resize(unsafe_uninit_length=total + 4)
        var o = self.letters.unsafe_ptr()
        var s = seq.unsafe_ptr()
        comptime ACGT = SIMD[DType.uint8, 4](65, 67, 71, 84)
        var i = 0
        while i < total:
            var x = s[unsafe_offset=i >> 2]
            o[unsafe_offset=i] = ACGT[Int(x & 3)]
            o[unsafe_offset=i + 1] = ACGT[Int((x >> 2) & 3)]
            o[unsafe_offset=i + 2] = ACGT[Int((x >> 4) & 3)]
            o[unsafe_offset=i + 3] = ACGT[Int(x >> 6)]
            i += 4
        self.letters.resize(unsafe_uninit_length=total)
        var exc = Cursor(cols[EXC])
        var at = 0
        while exc.at < exc.end:
            at += Int(exc.varint())
            if at >= total:
                raise _damaged("an exception past the letters")
            self.letters[at] = exc.byte()

    @always_inline
    def count(self) -> Int:
        return len(self.ids)


def _column(cols: List[List[UInt8]], c: Int) -> Cursor[ImmUntrackedOrigin]:
    """Column `c`, read through a pointer whose owner (the Block holding `cols`) keeps it alive."""
    return Cursor(untracked(cols[c].unsafe_ptr()), len(cols[c]))


struct Block[mut: Bool, //, origin: Origin[mut=mut]](Movable):
    """One block's records, read in turn (`next`): each one's POS, ID, padding base, walks (node-table index << 1 |
    reverse), LV and PS. It keeps the columns alive."""

    var cols: Pointer[List[List[UInt8]], Self.origin]
    var pos: Cursor[ImmUntrackedOrigin]
    var flags: Cursor[ImmUntrackedOrigin]
    var pad: Cursor[ImmUntrackedOrigin]
    var bound: Cursor[ImmUntrackedOrigin]
    var nwalk: Cursor[ImmUntrackedOrigin]
    var wlen: Cursor[ImmUntrackedOrigin]
    var walk: Cursor[ImmUntrackedOrigin]
    var lv: Cursor[ImmUntrackedOrigin]
    var ps: Cursor[ImmUntrackedOrigin]
    var nodes: Int
    var at_pos: Int
    var prev_first: UInt64
    var widx: Int
    # The record read last.
    var f: UInt8
    var pad_base: UInt8
    var first: UInt64
    var last: UInt64
    var level: Int
    var ps_start: UInt64
    var ps_end: UInt64
    var visits: List[UInt64]
    var bounds: List[Int]

    def __init__(out self, ref[Self.origin] cols: List[List[UInt8]], first_pos: Int, nodes: Int):
        self.cols = Pointer(to=cols)
        self.pos = _column(cols, POS)
        self.flags = _column(cols, FLAGS)
        self.pad = _column(cols, PAD)
        self.bound = _column(cols, BOUND)
        self.nwalk = _column(cols, NWALK)
        self.wlen = _column(cols, WLEN)
        self.walk = _column(cols, WALK)
        self.lv = _column(cols, LV)
        self.ps = _column(cols, PS)
        self.nodes = nodes
        self.at_pos = first_pos
        self.prev_first = 0
        self.widx = 0
        self.f = 0
        self.pad_base = 0
        self.first = 0
        self.last = 0
        self.level = 0
        self.ps_start = 0
        self.ps_end = 0
        self.visits = List[UInt64]()
        self.bounds = List[Int]()

    def next(mut self, walks: Bool) raises:
        """The next record; its walks only with `walks` (the WALK, WLEN and NWALK columns must then be read)."""
        self.at_pos += Int(self.pos.varint())
        self.f = self.flags.byte()
        self.pad_base = self.pad.byte() if self.f & FLAG_PAD != 0 else 0
        self.first = UInt64(Int(self.prev_first) + unzigzag(self.bound.varint()))
        self.last = UInt64(Int(self.first) + unzigzag(self.bound.varint()))
        self.prev_first = self.first
        if walks:
            var n = Int(self.nwalk.varint())
            self.visits.clear()
            self.bounds.clear()
            self.bounds.append(0)
            for _ in range(n):
                var m = Int(self.wlen.varint())
                for _ in range(m):
                    var x = self.walk.varint()
                    self.widx += unzigzag(x >> 1)
                    if self.widx < 0 or self.widx >= self.nodes:
                        raise _damaged("a walk visits a node outside its table")
                    self.visits.append(UInt64(self.widx << 1) | (x & 1))
                self.bounds.append(len(self.visits))
            self.level = Int(self.lv.varint())
            if self.f & FLAG_PS != 0:
                var base = self.id_start()
                self.ps_start = UInt64(Int(base) + unzigzag(self.ps.varint()))
                self.ps_end = UInt64(Int(self.ps_start) + unzigzag(self.ps.varint()))

    @always_inline
    def id_start(self) -> UInt64:
        return self.last ^ 1 if self.f & FLAG_REVERSED != 0 else self.first

    @always_inline
    def id_end(self) -> UInt64:
        return self.first ^ 1 if self.f & FLAG_REVERSED != 0 else self.last

    @always_inline
    def walks(self) -> Int:
        return len(self.bounds) - 1


@always_inline
def _reserve_end(mut out: List[UInt8], n: Int) -> MutPointer[UInt8, MutUntrackedOrigin]:
    """Room for `n` more bytes at the end of `out` (its length grows by `n`): where they go."""
    var at = len(out)
    out.resize(unsafe_uninit_length=at + n)
    return mut_untracked(out.unsafe_ptr().unsafe_offset(at))


@always_inline
def _put_number(mut out: List[UInt8], v: Int):
    var at = len(out)
    out.resize(unsafe_uninit_length=at + 20)
    var n = put_padded(out.unsafe_ptr().unsafe_offset(at), v, 0)
    out.resize(unsafe_uninit_length=at + n)


@always_inline
def _put_visit(mut out: List[UInt8], v: UInt64):
    out.append(UInt8(ord("<")) if v & 1 == 1 else UInt8(ord(">")))
    _put_number(out, Int(v >> 1))


@always_inline
def _put_text(mut out: List[UInt8], s: StaticString):
    out.extend(s.as_bytes())


def block_text(
    s: Sites, b: Int, chrom: String, beg: Int, end: Int, z: Decompressor, mut cols: List[List[UInt8]],
    mut table: NodeTable, mut out: List[UInt8], fd: Int32,
) raises:
    """Block `b`'s records overlapping [beg, end] (POS ..= POS + len(REF) - 1, a REF of no bases as one), as
    pgcatalog view prints them, appended to `out`; with `fd` >= 0, written there whenever `out` grows past
    STREAM_TEXT."""
    s.columns(b, (1 << COLUMNS) - 1, z, cols)
    table.read(cols, True)
    var r = Block(cols, s.first_pos(b), table.count())
    var lt = table.letters.unsafe_ptr()
    for _ in range(s.records(b)):
        r.next(True)
        if r.walks() < 1:
            raise _damaged("a record without walks")
        var pad = Int(r.f & FLAG_PAD != 0)
        var ref_length = pad
        for j in range(r.bounds[0], r.bounds[1]):
            var k = Int(r.visits[j] >> 1)
            ref_length += table.starts[k + 1] - table.starts[k]
        if r.at_pos > end or r.at_pos + max(ref_length, 1) - 1 < beg:
            continue
        out.extend(chrom.as_bytes())
        out.append(UInt8(ord("\t")))
        _put_number(out, r.at_pos)
        out.append(UInt8(ord("\t")))
        _put_visit(out, r.id_start())
        _put_visit(out, r.id_end())
        out.append(UInt8(ord("\t")))
        for k in range(r.walks()):
            if k > 1:
                out.append(UInt8(ord(",")))
            if pad == 1:
                out.append(r.pad_base)
            for j in range(r.bounds[k], r.bounds[k + 1]):
                var v = r.visits[j]
                var node = Int(v >> 1)
                var a = table.starts[node]
                var n = table.starts[node + 1] - a
                var dst = _reserve_end(out, n)
                put_oriented(dst, Span(unsafe_ptr=lt.unsafe_offset(a), length=n), v & 1 == 1)
                if fd >= 0 and len(out) > STREAM_TEXT:
                    write_all(fd, out.unsafe_ptr(), len(out))
                    out.clear()
            if k == 0:
                out.append(UInt8(ord("\t")))
        _put_text(out, "\t60\t.\tAT=")
        for k in range(r.walks()):
            if k > 0:
                out.append(UInt8(ord(",")))
            _put_visit(out, r.first)
            for j in range(r.bounds[k], r.bounds[k + 1]):
                var v = r.visits[j]
                _put_visit(out, table.ids[Int(v >> 1)] << 1 | (v & 1))
            _put_visit(out, r.last)
            if fd >= 0 and len(out) > STREAM_TEXT:
                write_all(fd, out.unsafe_ptr(), len(out))
                out.clear()
        _put_text(out, ";LV=")
        _put_number(out, r.level)
        if r.f & FLAG_PS != 0:
            _put_text(out, ";PS=")
            _put_visit(out, r.ps_start)
            _put_visit(out, r.ps_end)
        out.append(UInt8(ord("\n")))
    if (
        r.pos.at != r.pos.end
        or r.flags.at != r.flags.end
        or r.pad.at != r.pad.end
        or r.bound.at != r.bound.end
        or r.nwalk.at != r.nwalk.end
        or r.wlen.at != r.wlen.end
        or r.walk.at != r.walk.end
        or r.lv.at != r.lv.end
        or r.ps.at != r.ps.end
    ):
        raise _damaged("columns left over")


@fieldwise_init
struct _ViewTask(ParallelTask):
    """Block `order[first + index]`'s lines overlapping [beg, end], into `texts[index]`."""

    var sites: ImmPointer[Sites, ImmUntrackedOrigin]
    var order: ImmPointer[Int, ImmUntrackedOrigin]
    var first: Int
    var chrom: ImmPointer[String, ImmUntrackedOrigin]
    var beg: Int
    var end: Int
    var texts: MutPointer[List[UInt8], MutUntrackedOrigin]
    var error: MutPointer[String, MutUntrackedOrigin]

    def run(self, index: Int):
        var b = self.order[unsafe_offset=self.first + index]
        try:
            var z = Decompressor()
            var cols = List[List[UInt8]](length=COLUMNS, fill=List[UInt8]())
            var table = NodeTable()
            block_text(
                self.sites[], b, self.chrom[], self.beg, self.end, z, cols, table, self.texts[unsafe_offset=index], -1
            )
        except e:
            self.error[] = String(e)


def records(s: Sites, k: Int, beg: Int, end: Int, threads: Int, fd: Int32, mut kept: List[UInt8]) raises:
    """Contig `k`'s records overlapping [beg, end] (1-based, inclusive), in order, as tabix would return them from the
    VCF: written to `fd`, or (`fd` -1) appended to `kept`. Blocks are decoded in parallel, up to VIEW_BYTES of
    columns at a time; one with more than STREAM_TEXT of columns alone, on this thread, its text written as made."""
    var chrom = s.contigs[k].chrom
    var span = s.overlapping(k, beg, end)
    var order = List[Int]()
    for b in range(span[0], span[1]):
        if s.max_end(b) >= beg:
            order.append(b)
    var first = 0
    while first < len(order):
        if s.raw(order[first]) > STREAM_TEXT and fd >= 0:
            var z = Decompressor()
            var cols = List[List[UInt8]](length=COLUMNS, fill=List[UInt8]())
            var table = NodeTable()
            var out = List[UInt8]()
            block_text(s, order[first], chrom, beg, end, z, cols, table, out, fd)
            write_all(fd, out.unsafe_ptr(), len(out))
            first += 1
            continue
        var last = first
        var bytes = 0
        while last < len(order) and (last == first or bytes + s.raw(order[last]) <= VIEW_BYTES):
            if s.raw(order[last]) > STREAM_TEXT and fd >= 0:
                break
            bytes += s.raw(order[last])
            last += 1
        var texts = List[List[UInt8]](length=last - first, fill=List[UInt8]())
        var error = String()
        var task = _ViewTask(
            untracked(Pointer(to=s)),
            imm_untracked(order.unsafe_ptr()),
            first,
            untracked(Pointer(to=chrom)),
            beg,
            end,
            mut_untracked(texts.unsafe_ptr()),
            mut_untracked(Pointer(to=error)),
        )
        parallel_for(task, last - first, threads)
        # The tasks reach these through untracked pointers: reads after them keep them alive until they are done.
        _ = len(order) + chrom.byte_length()
        if error:
            raise Error(error)
        for t in texts:
            if fd >= 0:
                write_all(fd, t.unsafe_ptr(), len(t))
            else:
                kept.extend(Span(t))
        first = last


def _number(text: String, what: String) raises -> Int:
    """A region bound as pgcatalog reads it: commas dropped, then an unsigned integer."""
    var v = 0
    var digits = 0
    for b in text.as_bytes():
        if b == UInt8(ord(",")):
            continue
        if b < UInt8(ord("0")) or b > UInt8(ord("9")) or v > (Int.MAX - 9) // 10:
            raise Error(what)
        v = v * 10 + Int(b - UInt8(ord("0")))
        digits += 1
    if digits == 0:
        raise Error(what)
    return v


def view(path: String, region: Optional[String], threads: Int) raises:
    """`dromos view`: with a region (CONTIG or CONTIG:START-END), what pgcatalog view prints of that contig's
    catalog: its VCF header, then its records overlapping the region. Without one, every contig: one VCF whose
    header has a ##contig line for each, then each contig's records in turn."""
    var s = Sites(path)
    if not region:
        var contigs = List[Tuple[String, Int]]()
        for ref c in s.contigs:
            contigs.append((c.chrom, c.length))
        var h = vcf_header(contigs)
        write_all(1, h.unsafe_ptr(), h.byte_length())
        var none = List[UInt8]()
        for k in range(len(s.contigs)):
            records(s, k, 1, Int.MAX, threads, 1, none)
        return
    var q = region.value()
    var beg = 1
    var end = Int.MAX
    var name = q
    var colon = q.rfind(":")
    if colon >= 0:
        name = String(q[byte=:colon])
        var span = String(q[byte = colon + 1 :])
        var dash = span.find("-")
        if dash < 0:
            raise Error("region: CONTIG:START-END")
        beg = _number(String(span[byte=:dash]), "region start")
        end = _number(String(span[byte = dash + 1 :]), "region end")
    var k = s.contig(name)
    if not k:
        var names = List[String]()
        for ref c in s.contigs:
            names.append(c.name)
        raise Error(t"{path} has no contig {name} (it has {', '.join(names)})")
    ref c = s.contigs[k.value()]
    var contigs: List[Tuple[String, Int]] = [(c.chrom, c.length)]
    var h = vcf_header(contigs)
    write_all(1, h.unsafe_ptr(), h.byte_length())
    var none = List[UInt8]()
    records(s, k.value(), beg, end, threads, 1, none)


def info(path: String) raises:
    """`dromos catalog-info`: the header (the graph, its tags, the reference sample, the snarls' source, how it was
    built), then a line per contig, on stdout."""
    var s = Sites(path)
    var out = String(t"format\tdromos site catalog, version {VERSION}\n")
    out += String(t"file\t{path}\t{s.file.size} bytes\n")
    for p in s.props:
        out += String(t"{p[0]}\t{p[1]}\n")
    var total = 0
    for ref c in s.contigs:
        total += c.records
    out += String(t"contigs\t{len(s.contigs)}\nblocks\t{s.blocks()}\nrecords\t{total}\n")
    out += "#contig\tCHROM\tlength\tM5\trecords\tblocks\tregions\n"
    for ref c in s.contigs:
        var regions = String("whole")
        if len(c.regions) > 0:
            var parts = List[String]()
            for r in c.regions:
                parts.append(String(t"{r[0]}-{r[1]}"))
            regions = ",".join(parts)
        out += String(t"{c.name}\t{c.chrom}\t{c.length}\t{hex_string(Span(c.md5))}\t{c.records}\t{c.blocks}\t{regions}\n")
    write_all(1, out.unsafe_ptr(), out.byte_length())


# Surject's windows.


def _window(contig: Contig, map: RefMap, a: Int, b: Int, mut all: List[Window]):
    """The window between a record's ID nodes `a` and `b`, if the contig visits both once and the bases between
    them are few enough."""
    var pa = contig.placed(a)
    var pb = contig.placed(b)
    if pa < 0 or pb < 0:
        return
    var sa = pa >> 1
    var sb = pb >> 1
    var lo: Int
    var hi: Int
    if sa < sb:
        lo = sa + map.length(a)
        hi = sb
    else:
        lo = sb + map.length(b)
        hi = sa
    if hi >= lo and hi - lo <= MAX_SITE:
        all.append(Window(lo, hi))


def load_windows(s: Sites, k: Int, contig: Contig, map: RefMap) raises -> List[Window]:
    """Contig `k`'s snarl windows (sorted, disjoint): of each record, its ID's nodes (the FLAGS and BOUND columns
    alone are read)."""
    var all = List[Window]()
    var z = Decompressor()
    var cols = List[List[UInt8]](length=COLUMNS, fill=List[UInt8]())
    ref c = s.contigs[k]
    for b in range(c.first_block, c.first_block + c.blocks):
        s.columns(b, (1 << FLAGS) | (1 << BOUND) | (1 << POS) | (1 << PAD), z, cols)
        var r = Block(cols, s.first_pos(b), 0)
        for _ in range(s.records(b)):
            r.next(False)
            _window(contig, map, Int(r.id_start() >> 1), Int(r.id_end() >> 1), all)

    # By start, the widest first; repeats dropped; then windows inside the one before (but for an empty one at its
    # end) dropped.
    def before(a: Window, b: Window) -> Bool:
        return a.lo < b.lo or (a.lo == b.lo and a.hi > b.hi)

    sort(Span(all), before)
    var outer = List[Window]()
    for i in range(len(all)):
        var w = all[i]
        if i > 0 and w == all[i - 1]:
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
        var n = len(out)
        if n > 0 and lo <= out[n - 1].hi and max(hi, out[n - 1].hi) - out[n - 1].lo <= MAX_WINDOW:
            out[n - 1].hi = max(out[n - 1].hi, hi)
        else:
            out.append(Window(lo, hi))
    return out^


@fieldwise_init
struct _LoadWindows[map_origin: ImmOrigin, ref_origin: ImmOrigin](ParallelTask):
    """The catalog's windows on contig `i` of the reference, a contig a task."""

    var sites: ImmPointer[Sites, ImmUntrackedOrigin]
    var map: ImmPointer[RefMap, Self.map_origin]
    var reference: ImmPointer[Reference, Self.ref_origin]
    var which: ImmPointer[Int, ImmUntrackedOrigin]
    var out: MutPointer[List[Window], MutUntrackedOrigin]
    var error: MutPointer[String, MutUntrackedOrigin]

    def run(self, index: Int):
        try:
            self.out[unsafe_offset=index] = load_windows(
                self.sites[], self.which[unsafe_offset=index], self.reference[].contigs[index], self.map[]
            )
        except e:
            self.error[] = String(e)


def load_all_windows(
    s: Sites, names: List[String], reference: Reference, map: RefMap, threads: Int
) raises -> List[List[Window]]:
    """Each contig's windows (`names[i]` for the reference's contig i), contigs read at once. Every contig must be
    in the catalog."""
    var which = List[Int]()
    for n in names:
        var k = s.contig(n)
        if not k:
            raise Error(t"the catalog {s.path} has no contig {n}")
        which.append(k.value())
    var out = List[List[Window]](length=len(names), fill=List[Window]())
    var error = String()
    var task = _LoadWindows(
        untracked(Pointer(to=s)),
        Pointer(to=map),
        Pointer(to=reference),
        imm_untracked(which.unsafe_ptr()),
        out.unsafe_ptr().unsafe_origin_cast[MutUntrackedOrigin](),
        mut_untracked(Pointer(to=error)),
    )
    parallel_for(task, len(names), threads)
    _ = len(which)
    if error:
        raise Error(error)
    return out^
