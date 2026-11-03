"""BGZF in pieces: a file read with pread and the kernel's readahead, the whole blocks of a chunk found, and a range
of blocks inflated by libdeflate with every block's CRC32 checked (after sporos's fasta.mojo, `_inflate_bgzf`)."""

from std.ffi import OwnedDLHandle, external_call
from std.memory import unsafe_memcpy
from std.os import getenv, stat
from std.sys.info import CompilationTarget


comptime CHUNK = 64 << 20
"""Compressed bytes read at a time."""


def load_libdeflate() raises -> OwnedDLHandle:
    """libdeflate, from the environment pixi sets up or the system's."""
    comptime suffix = ".dylib" if CompilationTarget.is_macos() else ".so"
    var names: List[String] = [
        "@rpath/libdeflate.0.dylib" if CompilationTarget.is_macos() else "libdeflate.so.0",
        String("libdeflate") + suffix,
        "/opt/homebrew/lib/libdeflate.dylib",
    ]
    var prefix = getenv("CONDA_PREFIX")
    if prefix:
        names.append(prefix + "/lib/libdeflate" + suffix)
    for name in names:
        try:
            return OwnedDLHandle(name)
        except:
            pass
    raise Error("libdeflate not found (run under pixi, or install it)")


comptime O_RDONLY: Int32 = 0
comptime O_DIRECT: Int32 = 0o40000
"""Linux's O_DIRECT (x86-64 and aarch64): reads straight from the device into the caller's memory."""
comptime DIRECT_ALIGN = 4096
"""What O_DIRECT reads' memory, offsets and lengths are aligned to."""
comptime POSIX_FADV_WILLNEED: Int32 = 3


struct RawFile(Movable):
    """A file read with pread. With `direct` on Linux it is opened O_DIRECT where the file system allows it: reads go
    from the device into the caller's memory, past the page cache (whose copying and reclaiming hold a buffered read
    of a file larger than memory to a fraction of an NVMe's speed). Otherwise the kernel reads ahead of `prefetch`."""

    var fd: Int32
    var path: String
    var size: Int
    var direct: Bool

    def __init__(out self, path: String, direct: Bool = False) raises:
        var c_path = List[UInt8](path.as_bytes())
        c_path.append(0)
        comptime open_name = "open64" if CompilationTarget.is_linux() else "open"
        self.direct = False
        self.fd = -1
        comptime if CompilationTarget.is_linux():
            if direct:
                self.fd = external_call[open_name, Int32](c_path.unsafe_ptr(), O_RDONLY | O_DIRECT, Int32(0))
                self.direct = self.fd >= 0
        if self.fd < 0:
            self.fd = external_call[open_name, Int32](c_path.unsafe_ptr(), O_RDONLY, Int32(0))
        if self.fd < 0:
            raise Error(t"cannot open {path}")
        self.path = path
        self.size = Int(stat(path).st_size)

    def read(self, dst: MutPointer[UInt8, _], count: Int, offset: Int) -> Int:
        """Up to `count` bytes from `offset` into `dst`: fewer only at the end of the file; -1 on an error. A direct
        file's `dst`, `count` and `offset` are multiples of DIRECT_ALIGN."""
        var got = 0
        while got < count and offset + got < self.size:
            var n = external_call["pread", Int](self.fd, dst.unsafe_offset(got), count - got, offset + got)
            if n < 0:
                return -1
            if n == 0:
                break
            got += n
        return got

    def prefetch(self, offset: Int, count: Int):
        if self.direct:
            return
        comptime if CompilationTarget.is_linux():
            _ = external_call["posix_fadvise", Int32](self.fd, offset, count, POSIX_FADV_WILLNEED)

    def __deinit__(deinit self):
        _ = external_call["close", Int32](self.fd)


def scan_blocks(data: Span[UInt8, _], base: Int, mut blocks: List[Int], path: String) raises -> Tuple[Int, Int]:
    """Appends, for each whole BGZF block at the front of `data`, its deflate data's offset and size and its output's
    offset (from `base`) and size; returns the bytes those blocks take and their output's total size."""
    var n = len(data)
    var total = 0
    var at = 0
    while at + 18 <= n:
        if data[at] != 0x1F or data[at + 1] != 0x8B or data[at + 3] & 4 == 0:
            raise Error(t"{path}: not BGZF")
        var xlen = Int(data[at + 10]) | (Int(data[at + 11]) << 8)
        if at + 12 + xlen > n:
            break
        var bsize = -1
        var x = at + 12
        while x + 4 <= at + 12 + xlen:
            var sub_len = Int(data[x + 2]) | (Int(data[x + 3]) << 8)
            if data[x] == 66 and data[x + 1] == 67 and sub_len == 2:
                bsize = Int(data[x + 4]) | (Int(data[x + 5]) << 8)
            x += 4 + sub_len
        if bsize < xlen + 19:
            raise Error(t"{path}: corrupt BGZF block")
        var end = at + bsize + 1
        if end > n:
            break
        var isize = (
            Int(data[end - 4]) | (Int(data[end - 3]) << 8) | (Int(data[end - 2]) << 16) | (Int(data[end - 1]) << 24)
        )
        blocks.append(at + 12 + xlen)
        blocks.append(end - 8 - (at + 12 + xlen))
        blocks.append(base + total)
        blocks.append(isize)
        total += isize
        at = end
    return (at, total)


def inflate_blocks(
    lib: OwnedDLHandle,
    data: ImmPointer[UInt8, _],
    blocks: ImmPointer[Int, _],
    dst: MutPointer[UInt8, _],
    first: Int,
    last: Int,
) -> Bool:
    """Inflates blocks `first ..< last` (as `scan_blocks` lists them) from `data` into `dst`, each checked against
    its CRC32; false if any fails."""
    try:
        var decompressor = lib.get_function[Int]("libdeflate_alloc_decompressor")()
        var inflate = lib.get_function[Int32]("libdeflate_deflate_decompress")
        var crc32 = lib.get_function[UInt32]("libdeflate_crc32")
        var ok = True
        for b in range(first, last):
            var block = blocks.unsafe_offset(4 * b)
            var out = dst.unsafe_offset(block[unsafe_offset=2])
            var status = inflate(
                decompressor,
                data.unsafe_offset(block[unsafe_offset=0]),
                block[unsafe_offset=1],
                out,
                block[unsafe_offset=3],
                Int(0),
            )
            if status != 0:
                ok = False
                continue
            var t = data.unsafe_offset(block[unsafe_offset=0] + block[unsafe_offset=1])
            var stored = (
                UInt32(t[unsafe_offset=0])
                | (UInt32(t[unsafe_offset=1]) << 8)
                | (UInt32(t[unsafe_offset=2]) << 16)
                | (UInt32(t[unsafe_offset=3]) << 24)
            )
            ok = ok and crc32(UInt32(0), out, block[unsafe_offset=3]) == stored
        lib.get_function[NoneType]("libdeflate_free_decompressor")(decompressor)
        return ok
    except:
        return False


comptime BLOCK = 0xFF00
"""Uncompressed bytes per block written, as bgzip and htslib cut them."""
comptime MAX_BLOCK_SIZE = 0x10000
comptime _HEADER = SIMD[DType.uint8, 16](0x1F, 0x8B, 0x08, 0x04, 0, 0, 0, 0, 0, 0xFF, 0x06, 0, 0x42, 0x43, 0x02, 0)
"""A BGZF block's gzip header with its BC extra field, the block size to follow."""


def eof_block(mut out: List[UInt8]):
    """Appends the empty block that ends every BGZF file."""
    for i in range(16):
        out.append(_HEADER[i])
    for b in [0x1B, 0, 0x03, 0, 0, 0, 0, 0, 0, 0, 0, 0]:
        out.append(UInt8(b))


struct Deflater(Movable):
    """libdeflate compressors for BGZF blocks: one at the level asked for, and a level-0 one for data that will not
    fit a block compressed (stored, as htslib does)."""

    var compressor: Int
    var stored: Int
    var level: Int

    def __init__(out self, lib: OwnedDLHandle, level: Int) raises:
        var alloc = lib.get_function[Int]("libdeflate_alloc_compressor")
        self.level = level
        self.compressor = alloc(Int32(level))
        self.stored = alloc(Int32(0))
        if self.compressor == 0 or self.stored == 0:
            raise Error("libdeflate: no compressor")

    def block(self, lib: OwnedDLHandle, data: ImmPointer[UInt8, _], length: Int, mut out: List[UInt8]) raises:
        """Appends `data[:length]` (at most BLOCK bytes) as one BGZF block: header, raw deflate, CRC32, ISIZE."""
        var compress = lib.get_function[Int]("libdeflate_deflate_compress")
        var crc32 = lib.get_function[UInt32]("libdeflate_crc32")
        var at = len(out)
        var room = MAX_BLOCK_SIZE - 26
        out.resize(unsafe_uninit_length=at + 18 + room + 8)
        var body = out.unsafe_ptr().unsafe_offset(at + 18)
        var n = compress(self.compressor, data, length, body, room)
        if n == 0:
            n = compress(self.stored, data, length, body, room)
        if n == 0:
            raise Error("BGZF block overflow")
        var size = n + 26
        var p = out.unsafe_ptr().unsafe_offset(at)
        for i in range(16):
            p[unsafe_offset=i] = _HEADER[i]
        p[unsafe_offset=16] = UInt8((size - 1) & 0xFF)
        p[unsafe_offset=17] = UInt8((size - 1) >> 8)
        var crc = Int(crc32(UInt32(0), data, length))
        for k in range(4):
            p[unsafe_offset=18 + n + k] = UInt8((crc >> (8 * k)) & 0xFF)
            p[unsafe_offset=22 + n + k] = UInt8((length >> (8 * k)) & 0xFF)
        out.resize(unsafe_uninit_length=at + size)

    def free(self, lib: OwnedDLHandle):
        """Releases the compressors (once, when done with them)."""
        try:
            var release = lib.get_function[NoneType]("libdeflate_free_compressor")
            release(self.compressor)
            release(self.stored)
        except:
            pass


def bgzf_block(lib: OwnedDLHandle, level: Int, data: ImmPointer[UInt8, _], length: Int, mut out: List[UInt8]) raises:
    """Appends one BGZF block at `level` (a compressor made for the one block)."""
    var d = Deflater(lib, level)
    d.block(lib, data, length, out)
    d.free(lib)


def load_zstd() raises -> OwnedDLHandle:
    """libzstd, from the environment pixi sets up or the system's."""
    comptime suffix = ".dylib" if CompilationTarget.is_macos() else ".so"
    var names: List[String] = [
        "@rpath/libzstd.1.dylib" if CompilationTarget.is_macos() else "libzstd.so.1",
        String("libzstd") + suffix,
        "/opt/homebrew/lib/libzstd.dylib",
    ]
    var prefix = getenv("CONDA_PREFIX")
    if prefix:
        names.append(prefix + "/lib/libzstd" + suffix)
    for name in names:
        try:
            return OwnedDLHandle(name)
        except:
            pass
    raise Error("libzstd not found (run under pixi, or install it)")


def zstd_compress(lib: OwnedDLHandle, data: ImmPointer[UInt8, _], length: Int, level: Int) raises -> List[UInt8]:
    """`data[:length]` as one zstd frame."""
    var bound = lib.get_function[Int]("ZSTD_compressBound")(length)
    var out = List[UInt8](unsafe_uninit_length=bound)
    var n = lib.get_function[Int]("ZSTD_compress")(out.unsafe_ptr(), bound, data, length, Int32(level))
    if lib.get_function[UInt32]("ZSTD_isError")(n) != 0:
        raise Error("zstd: compression failed")
    # A copy of the frame's size: the bound's capacity would stay with every frame kept.
    var frame = List[UInt8](unsafe_uninit_length=n)
    unsafe_memcpy(dest=frame.unsafe_ptr(), src=out.unsafe_ptr(), count=n)
    return frame^


def zstd_decompress(lib: OwnedDLHandle, frame: Span[UInt8, _], dst: MutPointer[UInt8, _], length: Int) raises:
    """Inflates `frame` into `dst`, which it must fill exactly (`length` bytes)."""
    var n = lib.get_function[Int]("ZSTD_decompress")(dst, length, frame.unsafe_ptr(), len(frame))
    if lib.get_function[UInt32]("ZSTD_isError")(n) != 0 or n != length:
        raise Error("zstd: a frame did not decompress")


comptime PROT_READ: Int32 = 1
comptime MAP_PRIVATE: Int32 = 2
comptime MADV_WILLNEED: Int32 = 3


struct MappedFile(Movable):
    """A file mapped read-only into memory: its pages are read when first touched, by whichever thread touches them,
    and the kernel is asked to read the whole file ahead."""

    var ptr: ImmPointer[UInt8, ImmUntrackedOrigin]
    var size: Int

    def __init__(out self, path: String) raises:
        var file = RawFile(path)
        self.size = Int(stat(path).st_size)
        if self.size == 0:
            raise Error(t"{path}: empty")
        var addr = external_call["mmap", Int](Int(0), self.size, PROT_READ, MAP_PRIVATE, file.fd, Int(0))
        if addr == -1:
            raise Error(t"cannot map {path}")
        self.ptr = ImmPointer[UInt8, ImmUntrackedOrigin](unsafe_from_address=addr)
        _ = external_call["madvise", Int32](addr, self.size, MADV_WILLNEED)

    def __deinit__(deinit self):
        _ = external_call["munmap", Int32](self.ptr, self.size)
