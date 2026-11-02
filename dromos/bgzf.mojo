"""BGZF read a chunk at a time, each chunk's blocks inflated in parallel by libdeflate, every block's CRC32 checked.

From sporos's fasta.mojo (`_inflate_bgzf`), streamed: a file too large to hold is read in chunks of whole blocks,
the tail of a block cut by a chunk's end carried to the next.
"""

from std.ffi import OwnedDLHandle
from std.math import ceildiv
from std.memory import unsafe_memmove
from std.os import getenv
from std.sys.info import CompilationTarget

from dromos.parallel import ParallelTask, parallel_for

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


struct _InflateTask[data_origin: ImmOrigin, blocks_origin: ImmOrigin, lib_origin: ImmOrigin](ParallelTask):
    """Inflates one run of BGZF blocks into their slots of the output."""

    var data: ImmPointer[UInt8, Self.data_origin]
    var blocks: ImmPointer[Int, Self.blocks_origin]
    var dst: MutPointer[UInt8, MutUntrackedOrigin]
    var lib: ImmPointer[OwnedDLHandle, Self.lib_origin]
    var per_task: Int
    var count: Int
    var failed: MutPointer[Int, MutUntrackedOrigin]

    def __init__(
        out self,
        data: ImmPointer[UInt8, Self.data_origin],
        blocks: ImmPointer[Int, Self.blocks_origin],
        dst: MutPointer[UInt8, MutUntrackedOrigin],
        lib: ImmPointer[OwnedDLHandle, Self.lib_origin],
        per_task: Int,
        count: Int,
        failed: MutPointer[Int, MutUntrackedOrigin],
    ):
        self.data = data
        self.blocks = blocks
        self.dst = dst
        self.lib = lib
        self.per_task = per_task
        self.count = count
        self.failed = failed

    def run(self, index: Int):
        try:
            ref lib = self.lib[]
            var decompressor = lib.get_function[Int]("libdeflate_alloc_decompressor")()
            var inflate = lib.get_function[Int32]("libdeflate_deflate_decompress")
            var crc32 = lib.get_function[UInt32]("libdeflate_crc32")
            for b in range(index * self.per_task, min(self.count, (index + 1) * self.per_task)):
                # Block b: its deflate data's offset and size, and its output's offset and size.
                var block = self.blocks.unsafe_offset(4 * b)
                var status = inflate(
                    decompressor,
                    self.data.unsafe_offset(block[unsafe_offset=0]),
                    block[unsafe_offset=1],
                    self.dst.unsafe_offset(block[unsafe_offset=2]),
                    block[unsafe_offset=3],
                    Int(0),
                )
                if status != 0:
                    self.failed[unsafe_offset=0] = 1
                    continue
                var t = self.data.unsafe_offset(block[unsafe_offset=0] + block[unsafe_offset=1])
                var stored = (
                    UInt32(t[unsafe_offset=0])
                    | (UInt32(t[unsafe_offset=1]) << 8)
                    | (UInt32(t[unsafe_offset=2]) << 16)
                    | (UInt32(t[unsafe_offset=3]) << 24)
                )
                if crc32(UInt32(0), self.dst.unsafe_offset(block[unsafe_offset=2]), block[unsafe_offset=3]) != stored:
                    self.failed[unsafe_offset=0] = 1
            lib.get_function[NoneType]("libdeflate_free_decompressor")(decompressor)
        except:
            self.failed[unsafe_offset=0] = 1


struct BgzfReader(Movable):
    """A BGZF file, inflated a chunk of whole blocks at a time."""

    var file: FileHandle
    var path: String
    var lib: OwnedDLHandle
    var threads: Int
    var raw: List[UInt8]
    """Compressed bytes read but not yet inflated: whole blocks then, at the end, part of one."""
    var at_end: Bool

    def __init__(out self, path: String, threads: Int) raises:
        self.file = open(path, "r")
        self.path = path
        self.lib = load_libdeflate()
        self.threads = threads
        self.raw = List[UInt8]()
        self.at_end = False

    def next(mut self, mut out: List[UInt8]) raises -> Bool:
        """Appends the letters of the next chunk's whole blocks to `out`; false once the file is done."""
        if self.at_end and len(self.raw) == 0:
            return False
        if not self.at_end:
            var have = len(self.raw)
            self.raw.resize(unsafe_uninit_length=have + CHUNK)
            var got = self.file.read(Span(self.raw)[have : have + CHUNK])
            if got < 0:
                raise Error(t"cannot read {self.path}")
            self.raw.resize(unsafe_uninit_length=have + Int(got))
            self.at_end = got == 0
        var data = Span(self.raw)
        var n = len(data)
        var blocks = List[Int]()
        var total = 0
        var at = 0
        while at + 18 <= n:
            if data[at] != 0x1F or data[at + 1] != 0x8B or data[at + 3] & 4 == 0:
                raise Error(t"{self.path}: not BGZF")
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
                raise Error(t"{self.path}: corrupt BGZF block")
            var end = at + bsize + 1
            if end > n:
                break
            var isize = (
                Int(data[end - 4]) | (Int(data[end - 3]) << 8) | (Int(data[end - 2]) << 16) | (Int(data[end - 1]) << 24)
            )
            blocks.append(at + 12 + xlen)
            blocks.append(end - 8 - (at + 12 + xlen))
            blocks.append(len(out) + total)
            blocks.append(isize)
            total += isize
            at = end
        if at < n and self.at_end:
            raise Error(t"{self.path}: ends inside a BGZF block")
        var start = len(out)
        out.resize(unsafe_uninit_length=start + total)
        var count = len(blocks) // 4
        var per_task = max(1, ceildiv(count, max(1, self.threads) * 8))
        var failed = List[Int](length=1, fill=0)
        if count > 0:
            var task = _InflateTask(
                self.raw.unsafe_ptr().as_imm(),
                blocks.unsafe_ptr().as_imm(),
                out.unsafe_ptr().unsafe_origin_cast[MutUntrackedOrigin](),
                Pointer(to=self.lib),
                per_task,
                count,
                failed.unsafe_ptr().unsafe_origin_cast[MutUntrackedOrigin](),
            )
            parallel_for(task, ceildiv(count, per_task), self.threads)
        if failed[0] != 0:
            raise Error(t"{self.path}: BGZF inflation failed")
        # The cut block's bytes move to the front for the next chunk.
        var rest = n - at
        if rest > 0 and at > 0:
            var p = self.raw.unsafe_ptr().unsafe_origin_cast[MutUntrackedOrigin]()  # memmove allows the overlap
            unsafe_memmove(dest=p, src=p.unsafe_offset(at).as_imm(), count=rest)
        self.raw.resize(unsafe_uninit_length=rest)
        return True
