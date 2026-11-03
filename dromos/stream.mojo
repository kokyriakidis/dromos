"""A GAM (or a GAF) read as a pipeline of steps, each a set of tasks run together on every thread: some read the chunk
after next (so the disk starts first), some inflate the next chunk's BGZF blocks, and some hand slices of the current
chunk's messages (a GAF's lines) to a sink, beside any tasks of the sink's own. Between steps the sink takes what its
slices made, on one thread, while the next chunk is framed.

A GAF may also be plain text, its chunks copied in parallel, or gzip that is not BGZF, inflated as one stream by one
task of each step (gzip cannot be split)."""

from std.math import ceildiv
from std.ffi import OwnedDLHandle
from std.ffi import external_call
from std.memory import unsafe_memcpy
from std.os import getenv
from std.sys import stderr
from std.time import perf_counter_ns

from dromos.bgzf import (
    CHUNK,
    DIRECT_ALIGN,
    GzipStream,
    RawFile,
    inflate_blocks,
    load_libdeflate,
    load_zlib,
    scan_blocks,
)
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


comptime BGZF = 0
comptime PLAIN = 1
comptime GZIP = 2
comptime PIECE = 1 << 20
"""Bytes of plain text one task copies."""


def encoding(head: Span[UInt8, _]) -> Int:
    """How a file starting with `head` is compressed (pgcat's rule): BGZF when its first gzip member has the BC
    extra field first, gzip otherwise when it starts with gzip's magic, else plain."""
    var gz = len(head) >= 2 and head[0] == 0x1F and head[1] == 0x8B
    if gz and len(head) >= 16 and head[3] & 4 != 0 and head[12] == UInt8(ord("B")) and head[13] == UInt8(ord("C")):
        return BGZF
    return GZIP if gz else PLAIN


comptime _DECODE_FAILED = 1
comptime _TAKE_FAILED = 2
comptime _SIDE_FAILED = 3


@fieldwise_init
struct _Reads[file_origin: ImmOrigin](Copyable):
    """The reads of the chunk after next: `tasks` of them, a share of `length` bytes each from `offset` into `dst`,
    for the device to have several at once (one, for a pipe). Each task's byte count goes to `got`."""

    var file: ImmPointer[RawFile, Self.file_origin]
    var tasks: Int
    var dst: MutPointer[UInt8, MutUntrackedOrigin]
    var length: Int
    var offset: Int
    var got: MutPointer[Int, MutUntrackedOrigin]

    def run(self, i: Int):
        var share = self.length // self.tasks
        var got = self.file[].read(self.dst.unsafe_offset(i * share), share, self.offset + i * share)
        self.got[unsafe_offset=i] = got
        if got > 0 and i == self.tasks - 1:
            self.file[].prefetch(self.offset + self.length, CHUNK)


@fieldwise_init
struct _Plan(Copyable):
    """What the next chunk's bytes take to become letters (`_plan`): `blocks` (four ints each: source offset and
    size, destination offset and size) in `tasks` tasks of `per`, the bytes `consumed` (the rest carries to the chunk
    after) and the letters' `total` (unknown for gzip)."""

    var blocks: List[Int]
    var per: Int
    var tasks: Int
    var consumed: Int
    var total: Int


def _plan(mode: Int, raw: Span[UInt8, _], leftover: Int, threads: Int, path: String) raises -> _Plan:
    """The next chunk's whole BGZF blocks, plain text's pieces, or for gzip one task (it is one stream). Letters go
    after `leftover` bytes (the last chunk's unfinished message)."""
    var blocks = List[Int]()
    var length = len(raw)
    var consumed = length
    var total = 0
    if mode == BGZF:
        var scanned = scan_blocks(raw, leftover, blocks, path)
        consumed = scanned[0]
        total = scanned[1]
    elif mode == PLAIN:
        var at = 0
        while at < length:
            var piece = min(PIECE, length - at)
            blocks.append(at)
            blocks.append(piece)
            blocks.append(leftover + at)
            blocks.append(piece)
            at += piece
        total = length
    if mode == GZIP:
        var tasks = 1 if length > 0 else 0
        return _Plan(blocks^, 1, tasks, consumed, total)
    var count = len(blocks) // 4
    var per = max(1, ceildiv(count, 8 * threads))
    return _Plan(blocks^, per, ceildiv(count, per), consumed, total)


@fieldwise_init
struct _Decode[lib_origin: ImmOrigin, blocks_origin: ImmOrigin](Copyable):
    """The next chunk's bytes `raw[:length]` made letters in `dst` (`out`, a list that may grow, for gzip): task i
    inflates blocks i * per onward (BGZF), copies pieces (plain), or inflates the whole chunk (gzip)."""

    var mode: Int
    var lib: ImmPointer[OwnedDLHandle, Self.lib_origin]
    var zlib: ImmPointer[OwnedDLHandle, ImmUntrackedOrigin]
    var gzip: MutPointer[GzipStream, MutUntrackedOrigin]
    var raw: ImmPointer[UInt8, ImmUntrackedOrigin]
    var length: Int
    var blocks: ImmPointer[Int, Self.blocks_origin]
    var count: Int
    var per: Int
    var tasks: Int
    var dst: MutPointer[UInt8, MutUntrackedOrigin]
    var out: MutPointer[List[UInt8], MutUntrackedOrigin]

    def run(self, i: Int) -> Bool:
        """Task `i`; false if its input does not inflate."""
        if self.mode == GZIP:
            try:
                self.gzip[].inflate(self.zlib[], self.raw, self.length, self.out[])
                return True
            except:
                return False
        var first = i * self.per
        var last = min(self.count, first + self.per)
        if self.mode == PLAIN:
            for b in range(first, last):
                var block = self.blocks.unsafe_offset(4 * b)
                unsafe_memcpy(
                    dest=self.dst.unsafe_offset(block[unsafe_offset=2]),
                    src=self.raw.unsafe_offset(block[unsafe_offset=0]),
                    count=block[unsafe_offset=1],
                )
            return True
        return inflate_blocks(self.lib[], self.raw, self.blocks, self.dst, first, last)


@fieldwise_init
struct _Takes[S: MessageSink](Copyable):
    """The current chunk's messages handed to the sink in `tasks` slices of `per` (messages[2k], messages[2k + 1]:
    message k's start and end in `data`), beside `sides` side tasks of the sink's own."""

    var sink: MutPointer[Self.S, MutUntrackedOrigin]
    var data: ImmPointer[UInt8, ImmUntrackedOrigin]
    var messages: ImmPointer[Int, ImmUntrackedOrigin]
    var count: Int
    var per: Int
    var tasks: Int
    var sides: Int


@fieldwise_init
struct _StepTask[S: MessageSink, file_origin: ImmOrigin, lib_origin: ImmOrigin, blocks_origin: ImmOrigin](ParallelTask):
    """One step: the reads, then the decoding tasks, then the sink's side tasks, then its slices. The file, the
    libraries and the block list are reached through tracked pointers, which keep them alive while the step runs; the
    chunk buffers and the sink, used after every step, through untracked ones."""

    var reads: _Reads[Self.file_origin]
    var decode: _Decode[Self.lib_origin, Self.blocks_origin]
    var takes: _Takes[Self.S]
    var failed: MutPointer[Int, MutUntrackedOrigin]

    def run(self, index: Int):
        var i = index
        if i < self.reads.tasks:
            self.reads.run(i)
            return
        i -= self.reads.tasks
        if i < self.decode.tasks:
            if not self.decode.run(i):
                self.failed[] = _DECODE_FAILED
            return
        i -= self.decode.tasks
        ref t = self.takes
        if i < t.sides:
            try:
                t.sink[].side_task(i)
            except:
                self.failed[] = _SIDE_FAILED
            return
        i -= t.sides
        try:
            t.sink[].take(i, t.data, t.messages, i * t.per, min(t.count, (i + 1) * t.per))
        except:
            self.failed[] = _TAKE_FAILED


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


def frame_lines(data: Span[UInt8, _], mut messages: List[Int], final: Bool, pieces: Int) -> Int:
    """Cuts the whole lines at the front of `data` into up to `pieces` pieces of whole lines, appending each piece's
    start and end (its last newline included): a slice finds the lines of its piece itself. Returns the bytes the
    pieces take. With `final` (no more data), what follows the last newline is a line too."""
    var p = data.unsafe_ptr()
    var n = len(data)
    var used = n
    if not final:
        var last = external_call["memrchr", Int](p, Int32(ord("\n")), n)
        used = 0 if last == 0 else last - Int(p) + 1
    var target = max(1, (used + pieces - 1) // pieces)
    var at = 0
    while at < used:
        var cut = min(used, at + target)
        if cut < used:
            var hit = external_call["memchr", Int](p.unsafe_offset(cut - 1), Int32(ord("\n")), used - cut + 1)
            cut = used if hit == 0 else hit - Int(p) + 1
        messages.append(at)
        messages.append(cut)
        at = cut
    return used


def stream_gam[S: MessageSink](gam_path: String, threads: Int, mut sink: S) raises:
    """Hands every message of the GAM at `gam_path`, in order, to `sink`."""
    stream_input[lines=False](gam_path, threads, sink)


def stream_gaf[S: MessageSink](gaf_path: String, threads: Int, mut sink: S) raises:
    """Hands every line of the GAF at `gaf_path` (plain, gzip or BGZF), in order, to `sink`: each message is a piece
    of whole lines (`frame_lines`), its lines the sink's to find."""
    stream_input[lines=True](gaf_path, threads, sink)


def stream_input[S: MessageSink, //, lines: Bool](path: String, threads: Int, mut sink: S) raises:
    """`stream_gam`, or with `lines` `stream_gaf`."""
    var file = RawFile(path, direct=True)
    var lib = load_libdeflate()
    var zlib = load_zlib()
    var gzip = GzipStream()
    var timing = Bool(getenv("DROMOS_TIMING"))
    var step_ns = 0
    var serial_ns = 0
    var sink_ns = 0
    # Two of each buffer, swapping roles each step: raws[r] the next chunk as read (the cut block before it), bufs[c]
    # the current chunk's letters (an unfinished message or line first).
    var raws = List[_Raw]()
    raws.append(_Raw())
    raws.append(_Raw())
    var bufs = List[List[UInt8]](length=2, fill=List[UInt8]())
    var slices = sink.slices()
    var failed = List[Int](length=1, fill=0)
    var got = List[Int](length=READS, fill=0)
    # A pipe is read in order, one read a step.
    var read_tasks = READS if file.seekable else 1
    var r = 0
    var c = 0
    # The first chunk is read before the pipeline starts.
    var n0 = file.read(mut_untracked(raws[r].buf.unsafe_ptr().unsafe_offset(raws[r].start)), CHUNK, 0)
    if n0 < 0:
        raise Error(t"cannot read {path}")
    raws[r].hi = raws[r].start + n0
    var mode = encoding(Span(raws[r].buf)[raws[r].lo : raws[r].hi])
    comptime what = "a GAF line" if lines else "a GAM message"
    if not lines and mode != BGZF:
        raise Error(t"{path}: not BGZF (a GAM is)")
    var offset = n0
    var at_end = n0 < CHUNK
    var messages = List[Int]()
    var framer = Framer()
    var used = 0
    var first = True
    while True:
        var t0 = Int(perf_counter_ns())
        # The next chunk's whole blocks (pieces, for plain text); the cut block's bytes carry to the chunk after it.
        var leftover = len(bufs[c]) - used
        var length = raws[r].hi - raws[r].lo
        var plan = _plan(mode, Span(raws[r].buf)[raws[r].lo : raws[r].hi], leftover, threads, path)
        var consumed = plan.consumed
        var total = plan.total
        var nb = 1 - c
        bufs[nb].resize(unsafe_uninit_length=leftover + total)
        var nr = 1 - r
        var carry = length - consumed
        var reads = 0 if at_end else read_tasks
        raws[nr].lo = raws[nr].start - carry
        raws[nr].hi = raws[nr].start
        if carry > 0:
            unsafe_memcpy(
                dest=mut_untracked(raws[nr].buf.unsafe_ptr()).unsafe_offset(raws[nr].lo),
                src=imm_untracked(raws[r].buf.unsafe_ptr()).unsafe_offset(raws[r].lo + consumed),
                count=carry,
            )
        var message_count = 0 if first else len(messages) // 2
        var messages_per = max(1, ceildiv(message_count, slices))
        var formats = ceildiv(message_count, messages_per)
        var sides = sink.side_tasks()
        for k in range(READS):
            got[k] = 0
        var t1 = Int(perf_counter_ns())
        var task = _StepTask(
            _Reads(
                Pointer(to=file).as_imm(),
                reads,
                mut_untracked(raws[nr].buf.unsafe_ptr()).unsafe_offset(raws[nr].start),
                CHUNK,
                offset,
                mut_untracked(got.unsafe_ptr()),
            ),
            _Decode(
                mode,
                Pointer(to=lib).as_imm(),
                imm_untracked(Pointer(to=zlib)),
                mut_untracked(Pointer(to=gzip)),
                imm_untracked(raws[r].buf.unsafe_ptr()).unsafe_offset(raws[r].lo),
                length,
                plan.blocks.unsafe_ptr().as_imm(),
                len(plan.blocks) // 4,
                plan.per,
                plan.tasks,
                mut_untracked(bufs[nb].unsafe_ptr()),
                mut_untracked(Pointer(to=bufs[nb])),
            ),
            _Takes(
                mut_untracked(Pointer(to=sink)),
                imm_untracked(bufs[c].unsafe_ptr()),
                imm_untracked(messages.unsafe_ptr()),
                message_count,
                messages_per,
                formats,
                sides,
            ),
            mut_untracked(failed.unsafe_ptr()),
        )
        parallel_for(task, reads + plan.tasks + sides + formats, threads)
        var t2 = Int(perf_counter_ns())
        if failed[0] == _DECODE_FAILED:
            raise Error(t"{path}: does not inflate (corrupt, or not what it is named)")
        if failed[0] == _TAKE_FAILED:
            raise Error(t"{path}: {what} could not be read")
        if failed[0] == _SIDE_FAILED:
            raise Error(t"{path}: a side task failed")
        # The chunk's bytes: its reads' shares, up to the first short one (the end of the file).
        var read_total = 0
        for k in range(reads):
            if got[k] < 0:
                raise Error(t"cannot read {path}")
            read_total += got[k]
            if got[k] < CHUNK // reads:
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
        comptime if lines:
            # The last of the input is in hand once nothing more was read and nothing is carried.
            used = frame_lines(Span(bufs[nb]), messages, reads == 0 and carry == 0, slices)
        else:
            used = framer.frame(Span(bufs[nb]), messages)
        c = nb
        r = nr
        if reads > 0:
            raws[r].hi = raws[r].start + read_total
            offset += read_total
            at_end = read_total < CHUNK
        var done = plan.tasks == 0 and len(messages) == 0 and not first and sink.side_tasks() == 0
        first = False
        var t4 = Int(perf_counter_ns())
        step_ns += t2 - t1
        sink_ns += t3 - t2
        serial_ns += (t1 - t0) + (t4 - t3)
        if done:
            break
    gzip.free(zlib)
    if timing:
        print(
            (
                t"dromos timing: reading {Float64(step_ns + sink_ns + serial_ns) * 1e-9} s (steps"
                t" {Float64(step_ns) * 1e-9}, between steps {Float64(sink_ns + serial_ns) * 1e-9})"
            ),
            file=stderr,
        )
    if raws[r].hi > raws[r].lo:
        raise Error(t"{path}: ends inside a BGZF block")
    if mode == GZIP and not gzip.finished():
        raise Error(t"{path}: ends inside a gzip member")
    if len(bufs[c]) - used > 0 or framer.remaining > 0:
        raise Error(t"{path}: ends inside a GAM group")
