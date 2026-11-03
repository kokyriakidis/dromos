"""A GAM read as a pipeline of steps, each a set of tasks run together on every thread: one reads the chunk after next
(so the disk starts first), some inflate the next chunk's BGZF blocks, and some hand slices of the current chunk's
messages to a sink, beside any tasks of the sink's own. Between steps the sink takes what its slices made, on one
thread, while the next chunk's groups are framed."""

from std.math import ceildiv
from std.ffi import OwnedDLHandle
from std.memory import unsafe_memcpy
from std.os import getenv
from std.sys import stderr
from std.time import perf_counter_ns

from dromos.bgzf import CHUNK, DIRECT_ALIGN, RawFile, inflate_blocks, load_libdeflate, scan_blocks
from dromos.gam import Framer
from dromos.parallel import ParallelTask, parallel_for


trait MessageSink:
    """What a GAM's messages become."""

    def slices(self) -> Int:
        """How many slices of a chunk's messages the sink takes at once (each with state of its own)."""
        ...

    def take(
        mut self,
        slice: Int,
        data: ImmPointer[UInt8, ImmUntrackedOrigin],
        messages: ImmPointer[Int, ImmUntrackedOrigin],
        first: Int,
        last: Int,
    ) raises:
        """Messages `first ..< last` of the chunk (`messages` holds each one's start and end in `data`), as slice
        `slice`: slices run at once, each touching only its own state."""
        ...

    def side_tasks(self) -> Int:
        """How many tasks of the sink's own the next step runs beside its slices."""
        ...

    def side_task(mut self, index: Int) raises:
        """Runs side task `index`, at once with the slices and the other side tasks."""
        ...

    def after_step(mut self, slices: Int) raises:
        """Takes what slices `0 ..< slices` made in the step just run (on one thread)."""
        ...


@fieldwise_init
struct _StepTask[S: MessageSink, file_origin: ImmOrigin, lib_origin: ImmOrigin, blocks_origin: ImmOrigin](ParallelTask):
    """One step: `reads` reads, then `inflates` inflates, then `sides` side tasks, then `formats` slices. The file,
    libdeflate and the block list are reached through tracked pointers, which keep them alive while the step runs; the
    chunk buffers and the sink, used after every step, through untracked ones."""

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
    var sides: Int
    var sink: MutPointer[Self.S, MutUntrackedOrigin]
    var failed: MutPointer[Int, MutUntrackedOrigin]

    def run(self, index: Int):
        var i = index
        if i < self.reads:
            # Read k of the chunk's READS, a share of it each, for the device to have several at once.
            var share = self.read_count // self.reads
            var got = self.file[].read(self.read_dst.unsafe_offset(i * share), share, self.read_offset + i * share)
            self.read_got[unsafe_offset=i] = got
            if got > 0 and i == self.reads - 1:
                self.file[].prefetch(self.read_offset + self.read_count, CHUNK)
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
        if i < self.sides:
            try:
                self.sink[].side_task(i)
            except:
                self.failed[] = 3
            return
        i -= self.sides
        try:
            self.sink[].take(
                i,
                self.data,
                self.messages,
                i * self.messages_per,
                min(self.message_count, (i + 1) * self.messages_per),
            )
        except:
            self.failed[] = 2


comptime READS = 4
"""The reads a chunk is read in, at once."""
comptime CARRY_ROOM = 1 << 16
"""Room before a chunk for the last one's cut block (a BGZF block is at most 64 KiB)."""


struct _Raw(Movable):
    """A chunk of the GAM as read: room for the last chunk's cut block, then the chunk, at an address aligned for
    O_DIRECT. Its bytes are buf[lo:hi]."""

    var buf: List[UInt8]
    var start: Int
    var lo: Int
    var hi: Int

    def __init__(out self):
        self.buf = List[UInt8](unsafe_uninit_length=DIRECT_ALIGN + CARRY_ROOM + CHUNK)
        var addr = Int(self.buf.unsafe_ptr())
        self.start = CARRY_ROOM + (DIRECT_ALIGN - addr % DIRECT_ALIGN) % DIRECT_ALIGN
        self.lo = self.start
        self.hi = self.start


@always_inline
def imm_untracked[T: AnyType](p: MutPointer[T, _]) -> ImmPointer[T, ImmUntrackedOrigin]:
    return p.unsafe_origin_cast[MutUntrackedOrigin]().as_imm().unsafe_origin_cast[ImmUntrackedOrigin]()


@always_inline
def mut_untracked[T: AnyType](p: MutPointer[T, _]) -> MutPointer[T, MutUntrackedOrigin]:
    return p.unsafe_origin_cast[MutUntrackedOrigin]()


def stream_gam[S: MessageSink](gam_path: String, threads: Int, mut sink: S) raises:
    """Hands every message of the GAM at `gam_path`, in order, to `sink`."""
    var file = RawFile(gam_path, direct=True)
    var lib = load_libdeflate()
    var timing = Bool(getenv("DROMOS_TIMING"))
    var step_ns = 0
    var serial_ns = 0
    var sink_ns = 0
    # Two of each buffer, swapping roles each step: raws[r] the next chunk's compressed bytes (the cut block before
    # them), bufs[c] the current chunk's letters (an unfinished message first).
    var raws = List[_Raw]()
    raws.append(_Raw())
    raws.append(_Raw())
    var bufs = List[List[UInt8]](length=2, fill=List[UInt8]())
    var slices = sink.slices()
    var failed = List[Int](length=1, fill=0)
    var got = List[Int](length=READS, fill=0)
    var r = 0
    var c = 0
    # The first chunk is read before the pipeline starts.
    var n0 = file.read(mut_untracked(raws[r].buf.unsafe_ptr().unsafe_offset(raws[r].start)), CHUNK, 0)
    if n0 < 0:
        raise Error(t"cannot read {gam_path}")
    raws[r].hi = raws[r].start + n0
    var offset = n0
    var at_end = n0 < CHUNK
    var messages = List[Int]()
    var framer = Framer()
    var used = 0
    var first = True
    while True:
        var t0 = Int(perf_counter_ns())
        # The next chunk's whole blocks; the cut block's bytes carry to the chunk after it.
        var blocks = List[Int]()
        var leftover = len(bufs[c]) - used
        var scanned = scan_blocks(Span(raws[r].buf)[raws[r].lo : raws[r].hi], leftover, blocks, gam_path)
        var consumed = scanned[0]
        var total = scanned[1]
        var nb = 1 - c
        bufs[nb].resize(unsafe_uninit_length=leftover + total)
        var nr = 1 - r
        var carry = raws[r].hi - raws[r].lo - consumed
        var reads = 0 if at_end else READS
        raws[nr].lo = raws[nr].start - carry
        raws[nr].hi = raws[nr].start
        if carry > 0:
            unsafe_memcpy(
                dest=mut_untracked(raws[nr].buf.unsafe_ptr()).unsafe_offset(raws[nr].lo),
                src=imm_untracked(raws[r].buf.unsafe_ptr()).unsafe_offset(raws[r].lo + consumed),
                count=carry,
            )
        var block_count = len(blocks) // 4
        var blocks_per = max(1, ceildiv(block_count, 8 * threads))
        var inflates = ceildiv(block_count, blocks_per)
        var message_count = 0 if first else len(messages) // 2
        var messages_per = max(1, ceildiv(message_count, slices))
        var formats = ceildiv(message_count, messages_per)
        var sides = sink.side_tasks()
        for k in range(READS):
            got[k] = 0
        var t1 = Int(perf_counter_ns())
        var task = _StepTask(
            Pointer(to=file).as_imm(),
            reads,
            mut_untracked(raws[nr].buf.unsafe_ptr()).unsafe_offset(raws[nr].start),
            CHUNK,
            offset,
            mut_untracked(got.unsafe_ptr()),
            Pointer(to=lib).as_imm(),
            imm_untracked(raws[r].buf.unsafe_ptr()).unsafe_offset(raws[r].lo),
            blocks.unsafe_ptr().as_imm(),
            mut_untracked(bufs[nb].unsafe_ptr()),
            block_count,
            blocks_per,
            inflates,
            imm_untracked(bufs[c].unsafe_ptr()),
            imm_untracked(messages.unsafe_ptr()),
            message_count,
            messages_per,
            formats,
            sides,
            mut_untracked(Pointer(to=sink)),
            mut_untracked(failed.unsafe_ptr()),
        )
        parallel_for(task, reads + inflates + formats + sides, threads)
        var t2 = Int(perf_counter_ns())
        if failed[0] == 1:
            raise Error(t"{gam_path}: BGZF inflation failed")
        if failed[0] == 2:
            raise Error(t"{gam_path}: a GAM message could not be read")
        if failed[0] == 3:
            raise Error(t"{gam_path}: a side task failed")
        # The chunk's bytes: its reads' shares, up to the first short one (the end of the file).
        var read_total = 0
        for k in range(reads):
            if got[k] < 0:
                raise Error(t"cannot read {gam_path}")
            read_total += got[k]
            if got[k] < CHUNK // READS:
                break
        sink.after_step(formats)
        var t3 = Int(perf_counter_ns())
        # The current chunk's unfinished message goes in front of the next chunk's letters, then the next is framed.
        if leftover > 0:
            unsafe_memcpy(
                dest=mut_untracked(bufs[nb].unsafe_ptr()),
                src=imm_untracked(bufs[c].unsafe_ptr()).unsafe_offset(used),
                count=leftover,
            )
        messages.clear()
        used = framer.frame(Span(bufs[nb]), messages)
        c = nb
        r = nr
        if reads > 0:
            raws[r].hi = raws[r].start + read_total
            offset += read_total
            at_end = read_total < CHUNK
        var done = block_count == 0 and len(messages) == 0 and not first and sink.side_tasks() == 0
        first = False
        var t4 = Int(perf_counter_ns())
        step_ns += t2 - t1
        sink_ns += t3 - t2
        serial_ns += (t1 - t0) + (t4 - t3)
        if done:
            break
    if timing:
        print(
            (
                t"dromos timing (s): steps {Float64(step_ns) * 1e-9}, sink {Float64(sink_ns) * 1e-9}, framing"
                t" {Float64(serial_ns) * 1e-9}"
            ),
            file=stderr,
        )
    if raws[r].hi > raws[r].lo:
        raise Error(t"{gam_path}: ends inside a BGZF block")
    if len(bufs[c]) - used > 0 or framer.remaining > 0:
        raise Error(t"{gam_path}: ends inside a GAM group")
