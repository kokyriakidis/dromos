"""BGZF in pieces: a file read with pread and the kernel's readahead, the whole blocks of a chunk found, and a range
of blocks inflated by libdeflate with every block's CRC32 checked (after sporos's fasta.mojo, `_inflate_bgzf`)."""

from std.ffi import OwnedDLHandle, external_call
from std.os import getenv
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
comptime POSIX_FADV_WILLNEED: Int32 = 3


struct RawFile(Movable):
    """A file read with pread, its next bytes fetched ahead by the kernel on Linux (posix_fadvise WILLNEED)."""

    var fd: Int32
    var path: String

    def __init__(out self, path: String) raises:
        var c_path = List[UInt8](path.as_bytes())
        c_path.append(0)
        self.fd = external_call["open64" if CompilationTarget.is_linux() else "open", Int32](
            c_path.unsafe_ptr(), O_RDONLY, Int32(0)
        )
        if self.fd < 0:
            raise Error(t"cannot open {path}")
        self.path = path

    def read(self, dst: MutPointer[UInt8, _], count: Int, offset: Int) -> Int:
        """Up to `count` bytes from `offset` into `dst`: fewer only at the end of the file; -1 on an error."""
        var got = 0
        while got < count:
            var n = external_call["pread", Int](self.fd, dst.unsafe_offset(got), count - got, offset + got)
            if n < 0:
                return -1
            if n == 0:
                break
            got += n
        return got

    def prefetch(self, offset: Int, count: Int):
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
