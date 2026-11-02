"""Dromos: pangenome graph alignments, fast.

    dromos gam-to-gaf [-t THREADS] [--slim] NODES GAM    GAM as GAF on stdout, the lines vg convert -G writes

NODES is the graph's node table (`pgcat node-table`; see dromos/nodes.mojo).
"""

from std.math import ceildiv
from std.ffi import OwnedDLHandle
from std.memory import unsafe_memcpy
from std.os import getenv
from std.sys import argv, exit, stderr
from std.time import perf_counter_ns

from dromos.bgzf import CHUNK, RawFile, inflate_blocks, load_libdeflate, scan_blocks
from dromos.gam import Alignment, frame, gaf_line
from dromos.nodes import NodeTable
from dromos.output import OutputBuffer
from dromos.parallel import ParallelTask, default_threads, parallel_for

comptime USAGE = """dromos gam-to-gaf [-t THREADS] [--slim] NODES GAM
  GAM as GAF on stdout, the lines vg convert -G writes for it (NODES: the graph's node table, pgcat node-table)
  --slim   each line's 12 GAF columns alone (no cs, bq, AS, dv; NODES needs no letters)
"""


@fieldwise_init
struct _StepTask[
    tags: Bool, file_origin: ImmOrigin, lib_origin: ImmOrigin, blocks_origin: ImmOrigin, nodes_origin: ImmOrigin
](ParallelTask):
    """One step of the pipeline, its tasks run together: task 0 reads the chunk after next (so the disk starts first),
    the next `inflates` inflate the next chunk's blocks, and the last `formats` write the GAF lines of slices of the
    current chunk's messages. The file, libdeflate, the block list and the node table are reached through tracked
    pointers, which keep them alive while the step runs; the chunk buffers, used after every step, through untracked
    ones."""

    var file: ImmPointer[RawFile, Self.file_origin]
    var reads: Int
    var read_dst: MutPointer[UInt8, MutUntrackedOrigin]
    var read_count: Int
    var read_offset: Int
    var read_got: MutPointer[Int, MutUntrackedOrigin]
    var lib: ImmPointer[OwnedDLHandle, Self.lib_origin]
    var raw: ImmPointer[UInt8, ImmUntrackedOrigin]
    var blocks: ImmPointer[Int, Self.blocks_origin]
    var inflate_dst: MutPointer[UInt8, MutUntrackedOrigin]
    var block_count: Int
    var blocks_per: Int
    var inflates: Int
    var data: ImmPointer[UInt8, ImmUntrackedOrigin]
    var messages: ImmPointer[Int, ImmUntrackedOrigin]
    var message_count: Int
    var messages_per: Int
    var formats: Int
    var nodes: ImmPointer[NodeTable, Self.nodes_origin]
    var texts: MutPointer[List[UInt8], MutUntrackedOrigin]
    var failed: MutPointer[Int, MutUntrackedOrigin]

    def run(self, index: Int):
        var i = index
        if i < self.reads:
            var got = self.file[].read(self.read_dst, self.read_count, self.read_offset)
            self.read_got[] = got
            if got > 0:
                self.file[].prefetch(self.read_offset + got, CHUNK)
            return
        i -= self.reads
        if i < self.inflates:
            var first = i * self.blocks_per
            if not inflate_blocks(
                self.lib[],
                self.raw,
                self.blocks,
                self.inflate_dst,
                first,
                min(self.block_count, first + self.blocks_per),
            ):
                self.failed[] = 1
            return
        i -= self.inflates
        ref text = self.texts[unsafe_offset=i]
        text.clear()
        var a = Alignment()
        try:
            for m in range(i * self.messages_per, min(self.message_count, (i + 1) * self.messages_per)):
                a.decode(self.data, self.messages[unsafe_offset=2 * m], self.messages[unsafe_offset=2 * m + 1])
                gaf_line[Self.tags](self.data, a, self.nodes[], text)
        except:
            self.failed[] = 2


@always_inline
def _imm[T: AnyType](p: MutPointer[T, _]) -> ImmPointer[T, ImmUntrackedOrigin]:
    return p.unsafe_origin_cast[MutUntrackedOrigin]().as_imm().unsafe_origin_cast[ImmUntrackedOrigin]()


@always_inline
def _mut[T: AnyType](p: MutPointer[T, _]) -> MutPointer[T, MutUntrackedOrigin]:
    return p.unsafe_origin_cast[MutUntrackedOrigin]()


def gam_to_gaf[tags: Bool](nodes_path: String, gam_path: String, threads: Int) raises:
    """GAM to GAF, a pipeline of steps: while one chunk's messages are written as GAF, the next chunk's blocks inflate
    and the chunk after that is read, all on the same threads."""
    var nodes = NodeTable(nodes_path)
    var file = RawFile(gam_path)
    var lib = load_libdeflate()
    var out = OutputBuffer(1)
    out.write("@HD\tVN:Z:1.0\n")
    var timing = Bool(getenv("DROMOS_TIMING"))
    var step_ns = 0
    var serial_ns = 0
    # Two of each buffer, swapping roles each step: raws[r] the next chunk's compressed bytes, bufs[c] the current
    # chunk's letters (an unfinished group first).
    var raws = List[List[UInt8]](length=2, fill=List[UInt8]())
    var bufs = List[List[UInt8]](length=2, fill=List[UInt8]())
    var slices = 4 * threads
    var texts = List[List[UInt8]](length=slices, fill=List[UInt8]())
    var failed = List[Int](length=1, fill=0)
    var got = List[Int](length=1, fill=0)
    var r = 0
    var c = 0
    # The first chunk is read before the pipeline starts.
    raws[r].resize(unsafe_uninit_length=CHUNK)
    var n0 = file.read(_mut(raws[r].unsafe_ptr()), CHUNK, 0)
    if n0 < 0:
        raise Error(t"cannot read {gam_path}")
    raws[r].resize(unsafe_uninit_length=n0)
    var offset = n0
    var at_end = n0 < CHUNK
    var messages = List[Int]()
    var used = 0
    var first = True
    while True:
        var t0 = Int(perf_counter_ns())
        # The next chunk's whole blocks; the cut block's bytes carry to the chunk after it.
        var blocks = List[Int]()
        var leftover = len(bufs[c]) - used
        var scanned = scan_blocks(Span(raws[r]), leftover, blocks, gam_path)
        var consumed = scanned[0]
        var total = scanned[1]
        var nb = 1 - c
        bufs[nb].resize(unsafe_uninit_length=leftover + total)
        var nr = 1 - r
        var carry = len(raws[r]) - consumed
        var reads = 0 if at_end else 1
        raws[nr].resize(unsafe_uninit_length=carry + (CHUNK if reads == 1 else 0))
        if carry > 0:
            unsafe_memcpy(
                dest=_mut(raws[nr].unsafe_ptr()), src=_imm(raws[r].unsafe_ptr()).unsafe_offset(consumed), count=carry
            )
        var block_count = len(blocks) // 4
        var blocks_per = max(1, ceildiv(block_count, 8 * threads))
        var inflates = ceildiv(block_count, blocks_per)
        var message_count = 0 if first else len(messages) // 2
        var messages_per = max(1, ceildiv(message_count, slices))
        var formats = ceildiv(message_count, messages_per)
        got[0] = 0
        var t1 = Int(perf_counter_ns())
        var task = _StepTask[tags](
            Pointer(to=file).as_imm(),
            reads,
            _mut(raws[nr].unsafe_ptr().unsafe_offset(carry)),
            CHUNK,
            offset,
            _mut(got.unsafe_ptr()),
            Pointer(to=lib).as_imm(),
            _imm(raws[r].unsafe_ptr()),
            blocks.unsafe_ptr().as_imm(),
            _mut(bufs[nb].unsafe_ptr()),
            block_count,
            blocks_per,
            inflates,
            _imm(bufs[c].unsafe_ptr()),
            _imm(messages.unsafe_ptr()),
            message_count,
            messages_per,
            formats,
            Pointer(to=nodes).as_imm(),
            _mut(texts.unsafe_ptr()),
            _mut(failed.unsafe_ptr()),
        )
        parallel_for(task, reads + inflates + formats, threads)
        var t2 = Int(perf_counter_ns())
        if failed[0] == 1:
            raise Error(t"{gam_path}: BGZF inflation failed")
        if failed[0] == 2:
            raise Error(t"{gam_path}: a GAM message could not be read")
        if got[0] < 0:
            raise Error(t"cannot read {gam_path}")
        for t in range(formats):
            out.write_bytes(Span(texts[t]))
        # The current chunk's unfinished group goes in front of the next chunk's letters, then the next is framed.
        if leftover > 0:
            unsafe_memcpy(
                dest=_mut(bufs[nb].unsafe_ptr()), src=_imm(bufs[c].unsafe_ptr()).unsafe_offset(used), count=leftover
            )
        messages.clear()
        used = frame(Span(bufs[nb]), messages)
        c = nb
        r = nr
        if reads == 1:
            raws[r].resize(unsafe_uninit_length=carry + got[0])
            offset += got[0]
            at_end = got[0] < CHUNK
        var done = block_count == 0 and len(messages) == 0 and not first
        first = False
        var t3 = Int(perf_counter_ns())
        step_ns += t2 - t1
        serial_ns += (t1 - t0) + (t3 - t2)
        if done:
            break
    out.flush()
    if timing:
        print(
            t"dromos timing (s): parallel steps {Float64(step_ns) * 1e-9}, serial {Float64(serial_ns) * 1e-9}",
            file=stderr,
        )
    if len(raws[r]) > 0:
        raise Error(t"{gam_path}: ends inside a BGZF block")
    if len(bufs[c]) - used > 0:
        raise Error(t"{gam_path}: ends inside a GAM group")


def main():
    var args = argv()
    var threads = default_threads()
    var slim = False
    var positional = List[String]()
    var i = 1
    while i < len(args):
        var a = String(args[i])
        if a == "-t" and i + 1 < len(args):
            try:
                threads = max(1, Int(String(args[i + 1])))
            except:
                print("dromos: -t needs a number", file=stderr)
                exit(2)
            i += 2
            continue
        if a == "--slim":
            slim = True
            i += 1
            continue
        if a == "-h" or a == "--help":
            print(USAGE, end="")
            return
        positional.append(a)
        i += 1
    if len(positional) != 3 or positional[0] != "gam-to-gaf":
        print(USAGE, end="", file=stderr)
        exit(2)
    try:
        if slim:
            gam_to_gaf[False](positional[1], positional[2], threads)
        else:
            gam_to_gaf[True](positional[1], positional[2], threads)
    except e:
        print(t"dromos: {e}", file=stderr)
        exit(1)
