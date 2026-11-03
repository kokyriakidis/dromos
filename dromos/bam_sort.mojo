"""BAM records sorted as `samtools sort` sorts them (by reference, position, then reverse strand, then input order:
its comparator, a stable merge sort, and a merge that breaks ties by input order), written as a BGZF file with its
.bai or .csi index.

Records are staged by bucket as they come (a contig's BUCKET_SHIFT bases), compressed by zstd a FRAME at a time.
Once the frames held reach half the memory allowed, all of them are written to a spill file (a bucket's frames
together) by a thread of its own while staging goes on; the file is unlinked as soon as it is made, so nothing is
left behind. Each bucket keeps its records in input order: its spilled frames, its frames in memory, then its
records not yet compressed. At the end, windows of buckets are read back, inflated and sorted, laid out, deflated and
written in a three-deep pipeline (sort_write.mojo's), the index built as each block's address becomes known."""

from std.ffi import OwnedDLHandle, external_call
from std.memory import unsafe_memcpy
from std.os import getenv
from std.sys import stderr
from std.time import perf_counter_ns

from dromos.bam import BamIndex, get32, get_signed32, record_end
from dromos.bgzf import BLOCK, Deflater, load_libdeflate, load_zstd, zstd_compress_with, zstd_decompress_with
from dromos.output import OutFile, errno
from dromos.parallel import Background, ParallelTask, parallel_for
from dromos.sort_write import BLOCKS_PER_TASK, IndexWriter, Tab
from dromos.stream import imm_untracked, mut_untracked

comptime BUCKET_SHIFT = 20
"""A bucket is 2^20 reference bases of a contig."""
comptime FRAME = 1 << 20
"""A bucket's records are compressed once this many bytes of them gather."""
comptime ZSTD_LEVEL = 1
comptime WINDOW_BYTES = 1 << 30
"""The records a window of buckets holds at most (but for one bucket larger alone)."""


@fieldwise_init
struct _Frame(TrivialRegisterPassable):
    """A spilled frame: where it is in the spill file, its size there, and its size inflated."""

    var at: Int
    var size: Int
    var raw: Int


struct BamBucket(Copyable, Movable):
    """A bucket's records while the input streams: frames spilled, frames in memory (and each one's size inflated),
    then those not yet compressed."""

    var spilled: List[_Frame]
    var frames: List[List[UInt8]]
    var sizes: List[Int]
    var pending: List[UInt8]
    var records: Int
    var bytes: Int
    """The records' bytes, inflated."""

    def __init__(out self):
        self.spilled = List[_Frame]()
        self.frames = List[List[UInt8]]()
        self.sizes = List[Int]()
        self.pending = List[UInt8]()
        self.records = 0
        self.bytes = 0


@fieldwise_init
struct BamRec(TrivialRegisterPassable):
    """A record of a bucket being sorted: its key (position << 1 | reverse), where it is in the bucket's bytes, its
    length there (block_size included), and the end of its reference span."""

    var key: Int
    var off: Int
    var len: Int
    var end: Int


struct BamSlot(Movable):
    """A bucket of the window being written: its records inflated, sorted, and where they go in the window."""

    var tid: Int
    var bucket: Int
    var raw: List[UInt8]
    var recs: List[BamRec]
    var text_len: Int
    var at: Int

    def __init__(out self):
        self.tid = 0
        self.bucket = 0
        self.raw = List[UInt8]()
        self.recs = List[BamRec]()
        self.text_len = 0
        self.at = 0


@fieldwise_init
struct _SpillWrite(Movable, ParallelTask):
    """Writes frames at their places in the spill file (on a thread of its own)."""

    var fd: Int32
    var frames: List[List[UInt8]]
    var offsets: List[Int]
    var failed: MutPointer[Int, MutUntrackedOrigin]

    def run(self, index: Int):
        for k in range(len(self.frames)):
            if not _pwrite(self.fd, self.frames[k], self.offsets[k]):
                self.failed[] = Int(errno())
                return


def _pwrite(fd: Int32, data: List[UInt8], offset: Int) -> Bool:
    var done = 0
    var n = len(data)
    while done < n:
        var w = external_call["pwrite", Int](fd, data.unsafe_ptr().unsafe_offset(done), n - done, offset + done)
        if w <= 0:
            return False
        done += w
    return True


def _pread(fd: Int32, dst: MutPointer[UInt8, _], count: Int, offset: Int) -> Bool:
    var done = 0
    while done < count:
        var r = external_call["pread", Int](fd, dst.unsafe_offset(done), count - done, offset + done)
        if r <= 0:
            return False
        done += r
    return True


struct BamStaging(Movable):
    """The buckets, by contig number then bucket, and the spill file once one is needed."""

    var buckets: List[List[BamBucket]]
    var zstd: OwnedDLHandle
    var shares: Int
    """How many tasks the buckets are shared among when records are handed over (`take_share`)."""
    var held: List[Int]
    """Bytes of frames each share has made since they were last counted."""
    var in_memory: Int
    """Bytes of frames held in memory."""
    var limit: Int
    """Frames are spilled once they hold half this."""
    var spill_dir: String
    var spill: OutFile
    var spill_size: Int
    var spiller: List[Background[_SpillWrite]]
    """The spill being written, if one is (none or one)."""
    var spill_failed: List[Int]
    var spills: Int
    var spilled_bytes: Int
    var frame: Int
    """Bytes of a bucket's records compressed at a time."""
    var cctx: List[Int]
    """A zstd compression context for each share."""

    def __init__(out self, contigs: Int, shares: Int, limit: Int, spill_dir: String, frame: Int = FRAME) raises:
        self.buckets = List[List[BamBucket]](length=contigs, fill=List[BamBucket]())
        self.zstd = load_zstd()
        self.shares = shares
        self.held = List[Int](length=shares, fill=0)
        self.in_memory = 0
        self.limit = limit
        self.spill_dir = spill_dir
        self.spill = OutFile()
        self.spill_size = 0
        self.spiller = List[Background[_SpillWrite]]()
        self.spill_failed = List[Int](length=1, fill=0)
        self.spills = 0
        self.spilled_bytes = 0
        self.frame = frame
        self.cctx = List[Int]()
        for _ in range(shares):
            var c = self.zstd.get_function[Int]("ZSTD_createCCtx")()
            if c == 0:
                raise Error("zstd: no compression context")
            self.cctx.append(c)

    def __deinit__(deinit self):
        try:
            var release = self.zstd.get_function[Int]("ZSTD_freeCCtx")
            for c in self.cctx:
                _ = release(c)
        except:
            pass

    def make_room(mut self, recs: ImmPointer[UInt8, _], total: Int) raises:
        """Makes the buckets the records `recs[:total]` need (on one thread, before they are handed over)."""
        var at = 0
        while at < total:
            var tid = get32(recs.unsafe_offset(at + 4))
            var b = get_signed32(recs.unsafe_offset(at + 8)) >> BUCKET_SHIFT
            if tid >= len(self.buckets) or b < 0:
                raise Error("a BAM record of no contig")
            if len(self.buckets[tid]) <= b:
                self.buckets[tid].resize(b + 1, BamBucket())
            at += 4 + get32(recs.unsafe_offset(at))

    def take_share(mut self, share: Int, recs: ImmPointer[UInt8, _], total: Int) raises:
        """Takes the records of `recs[:total]` in share `share`'s buckets (bucket b of contig t: (t << 12 | b) mod
        shares), compressing a bucket's records as a frame once `frame` bytes of them gather. Shares run at once."""
        var at = 0
        while at < total:
            var length = 4 + get32(recs.unsafe_offset(at))
            var tid = get32(recs.unsafe_offset(at + 4))
            var b = get_signed32(recs.unsafe_offset(at + 8)) >> BUCKET_SHIFT
            if (tid << 12 | b) % self.shares == share:
                var bucket = mut_untracked(Pointer(to=self.buckets[tid][b]))
                bucket[].pending.extend(Span(unsafe_ptr=recs.unsafe_offset(at), length=length))
                bucket[].records += 1
                bucket[].bytes += length
                if len(bucket[].pending) >= self.frame:
                    var frame = zstd_compress_with(
                        self.zstd, self.cctx[share], bucket[].pending.unsafe_ptr(), len(bucket[].pending), ZSTD_LEVEL
                    )
                    self.held[share] += len(frame)
                    bucket[].frames.append(frame^)
                    bucket[].sizes.append(len(bucket[].pending))
                    bucket[].pending = List[UInt8]()
            at += length

    def after_shares(mut self) raises:
        """Counts the frames the shares made; once they hold half the memory allowed, spills them (on one thread,
        after the shares have run)."""
        for k in range(self.shares):
            self.in_memory += self.held[k]
            self.held[k] = 0
        if self.in_memory * 2 >= self.limit:
            self._spill()

    def _spill(mut self) raises:
        """Every frame in memory to the spill file, a bucket's together, written by a thread of its own; the spill
        before it waited for first."""
        self._join()
        if self.spill.fd < 0:
            var path = self.spill_dir + "/.dromos-spill." + String(external_call["getpid", Int32]())
            self.spill = OutFile(path, mappable=True)
            var c = List[UInt8](path.as_bytes())
            c.append(0)
            _ = external_call["unlink", Int32](c.unsafe_ptr())
        var frames = List[List[UInt8]]()
        var offsets = List[Int]()
        for ref list in self.buckets:
            for ref bucket in list:
                var held = List[List[UInt8]]()
                swap(held, bucket.frames)
                for k in range(len(held)):
                    var size = len(held[k])
                    bucket.spilled.append(_Frame(self.spill_size, size, bucket.sizes[k]))
                    offsets.append(self.spill_size)
                    self.spill_size += size
                frames.extend(held^)
                bucket.sizes = List[Int]()
        self.spills += 1
        self.spilled_bytes += self.in_memory
        self.in_memory = 0
        var task = _SpillWrite(self.spill.fd, frames^, offsets^, mut_untracked(self.spill_failed.unsafe_ptr()))
        self.spiller.append(Background(task^, 1, 1))

    def _join(mut self) raises:
        """Waits for the spill being written, if one is."""
        if len(self.spiller) > 0:
            self.spiller[0].join()
            self.spiller.clear()
        if self.spill_failed[0] != 0:
            raise Error(t"cannot write the spill file in {self.spill_dir} (errno {self.spill_failed[0]})")

    def finish_staging(mut self) raises:
        """Waits for the last spill, before the buckets are read back."""
        self._join()

    def order(self) -> List[Tuple[Int, Int]]:
        """The buckets with records, in output order: by contig, then position."""
        var order = List[Tuple[Int, Int]]()
        for t in range(len(self.buckets)):
            for b in range(len(self.buckets[t])):
                if self.buckets[t][b].records > 0:
                    order.append((t, b))
        return order^

    def bucket_bytes(self, tid: Int, b: Int) -> Int:
        return self.buckets[tid][b].bytes

    def sort(mut self, mut slot: BamSlot) raises:
        """Slot `slot`'s bucket read back and inflated into one buffer (in input order) and its records sorted as
        samtools sorts them; the bucket's memory freed."""
        var bucket = mut_untracked(Pointer(to=self.buckets[slot.tid][slot.bucket]))
        var dctx = self.zstd.get_function[Int]("ZSTD_createDCtx")()
        if dctx == 0:
            raise Error("zstd: no decompression context")
        try:
            self._inflate(bucket[], slot, dctx)
        finally:
            _ = self.zstd.get_function[Int]("ZSTD_freeDCtx")(dctx)
        slot.recs = List[BamRec](capacity=bucket[].records)
        bucket[].spilled = List[_Frame]()
        bucket[].frames = List[List[UInt8]]()
        bucket[].sizes = List[Int]()
        bucket[].pending = List[UInt8]()
        var total = bucket[].bytes
        var p = imm_untracked(slot.raw.unsafe_ptr())
        var at = 0
        slot.text_len = total
        while at < total:
            var length = 4 + get32(p.unsafe_offset(at))
            var r = p.unsafe_offset(at)
            var pos = get_signed32(r.unsafe_offset(8))
            var rev = (get32(r.unsafe_offset(16)) >> 20) & 1
            slot.recs.append(BamRec(pos << 1 | rev, at, length, record_end(r)))
            at += length

        def less(a: BamRec, b: BamRec) -> Bool:
            if a.key != b.key:
                return a.key < b.key
            return a.off < b.off

        sort(Span(slot.recs), less)

    def _inflate(self, bucket: BamBucket, mut slot: BamSlot, dctx: Int) raises:
        """A bucket's records into `slot.raw`, in input order: spilled frames read back, frames in memory, then the
        records not yet compressed."""
        var total = bucket.bytes
        slot.raw.resize(unsafe_uninit_length=total)
        var dst = slot.raw.unsafe_ptr()
        var at = 0
        # Spilled frames: consecutive ones read at once.
        var packed = List[UInt8]()
        var f = 0
        while f < len(bucket.spilled):
            var g = f + 1
            while (
                g < len(bucket.spilled)
                and bucket.spilled[g].at == bucket.spilled[g - 1].at + bucket.spilled[g - 1].size
            ):
                g += 1
            var start = bucket.spilled[f].at
            var span = bucket.spilled[g - 1].at + bucket.spilled[g - 1].size - start
            packed.resize(unsafe_uninit_length=span)
            if not _pread(self.spill.fd, packed.unsafe_ptr(), span, start):
                raise Error(t"cannot read the spill file back (errno {errno()})")
            for k in range(f, g):
                ref fr = bucket.spilled[k]
                zstd_decompress_with(
                    self.zstd,
                    dctx,
                    Span(unsafe_ptr=packed.unsafe_ptr().unsafe_offset(fr.at - start), length=fr.size),
                    dst.unsafe_offset(at),
                    fr.raw,
                )
                at += fr.raw
            f = g
        for k in range(len(bucket.frames)):
            zstd_decompress_with(self.zstd, dctx, Span(bucket.frames[k]), dst.unsafe_offset(at), bucket.sizes[k])
            at += bucket.sizes[k]
        if len(bucket.pending) > 0:
            unsafe_memcpy(dest=dst.unsafe_offset(at), src=bucket.pending.unsafe_ptr(), count=len(bucket.pending))
            at += len(bucket.pending)
        if at != total:
            raise Error("a bucket's records did not read back whole")

    def lay_out(self, slot: BamSlot, dst: MutPointer[UInt8, _]):
        """Slot `slot`'s records, sorted, at `dst`."""
        var src = slot.raw.unsafe_ptr()
        var at = slot.at
        for rec in slot.recs:
            unsafe_memcpy(dest=dst.unsafe_offset(at), src=src.unsafe_offset(rec.off), count=rec.len)
            at += rec.len


def write_bam(
    mut staging: BamStaging,
    header: List[UInt8],
    var index: BamIndex,
    out_path: String,
    level: Int,
    threads: Int,
) raises:
    """Writes `out_path` and its index (each through a temporary name): the BAM header, then every staged record,
    sorted."""
    var t0 = Int(perf_counter_ns())
    staging.finish_staging()
    var order = staging.order()
    var lib = load_libdeflate()
    var suffix = index.suffix()
    var w = IndexWriter(out_path, level, index^, suffix)
    var end = w.header_bytes(lib, Span(header))
    # Windows of buckets, by bytes: each step inflates and sorts window j, lays out window j - 1, deflates the whole
    # blocks of window j - 2 and writes window j - 3, all at once.
    var bounds = List[Int]()
    bounds.append(0)
    var held = 0
    for k in range(len(order)):
        var b = staging.bucket_bytes(order[k][0], order[k][1])
        if k > bounds[len(bounds) - 1] and (held + b > WINDOW_BYTES or k - bounds[len(bounds) - 1] >= 4 * threads):
            bounds.append(k)
            held = 0
        held += b
    if len(order) > 0:
        bounds.append(len(order))
    var windows = len(bounds) - 1
    var most = 0
    for j in range(windows):
        most = max(most, bounds[j + 1] - bounds[j])
    var slots = List[List[BamSlot]]()
    for _ in range(2):
        var group = List[BamSlot]()
        for _ in range(most):
            group.append(BamSlot())
        slots.append(group^)
    var win_start = List[Int](length=windows + 1, fill=0)
    var win_count = List[Int](length=max(windows, 1), fill=0)
    var deflated = List[Int](length=2, fill=0)
    var bad = List[Int](length=1, fill=0)
    var step_ns = 0
    var serial_ns = 0
    for j in range(windows + 3):
        var s0 = Int(perf_counter_ns())
        var gathers = 0
        if j < windows:
            gathers = bounds[j + 1] - bounds[j]
            win_count[j] = gathers
            for k in range(gathers):
                slots[j % 2][k].tid = order[bounds[j] + k][0]
                slots[j % 2][k].bucket = order[bounds[j] + k][1]
        var lay = j - 1
        var layouts = 0
        if lay >= 0 and lay < windows:
            layouts = win_count[lay]
            var x = lay % 2
            var start = end
            var base = (start // BLOCK) * BLOCK
            var at = start - base
            for k in range(layouts):
                slots[x][k].at = at
                at += slots[x][k].text_len
            w.texts[x].resize(unsafe_uninit_length=at)
            if lay > 0 and start > base:
                var y = 1 - x
                unsafe_memcpy(
                    dest=mut_untracked(w.texts[x].unsafe_ptr()),
                    src=imm_untracked(w.texts[y].unsafe_ptr()).unsafe_offset(base - w.text_base[y]),
                    count=start - base,
                )
            w.text_base[x] = base
            win_start[lay] = start
            end = base + at
            win_start[lay + 1] = end
        var dw = j - 2
        var deflates = 0
        var deflate_len = 0
        if dw >= 0 and dw < windows:
            var x = dw % 2
            deflate_len = (win_start[dw + 1] // BLOCK) * BLOCK - w.text_base[x]
            var blocks = deflate_len // BLOCK
            deflates = (blocks + BLOCKS_PER_TASK - 1) // BLOCKS_PER_TASK
            while len(w.outs[x]) < deflates:
                w.outs[x].append(List[UInt8]())
            deflated[x] = deflates
        var wr = j - 3
        var writes = 1 if wr >= 0 and wr < windows else 0
        var task = _BamTask(
            mut_untracked(Pointer(to=staging)),
            Pointer(to=lib).as_imm(),
            level,
            writes,
            mut_untracked(Pointer(to=w.file)),
            mut_untracked(w.outs[(j + 1) % 2].unsafe_ptr()),
            deflated[(j + 1) % 2] if writes == 1 else 0,
            gathers,
            mut_untracked(slots[j % 2].unsafe_ptr()),
            layouts,
            mut_untracked(slots[(j + 1) % 2].unsafe_ptr()),
            mut_untracked(w.texts[(j + 1) % 2].unsafe_ptr()),
            deflates,
            imm_untracked(w.texts[j % 2].unsafe_ptr()),
            deflate_len,
            mut_untracked(w.outs[j % 2].unsafe_ptr()),
            mut_untracked(bad.unsafe_ptr()),
        )
        var s1 = Int(perf_counter_ns())
        parallel_for(task, writes + gathers + layouts + deflates, threads)
        var s2 = Int(perf_counter_ns())
        if bad[0] == 1:
            raise Error("a bucket's records could not be read back")
        if bad[0] == 2:
            raise Error("BGZF: a block did not deflate")
        if layouts > 0:
            var x = lay % 2
            for k in range(layouts):
                var tid = slots[x][k].tid
                var stream_end = w.text_base[x] + slots[x][k].at
                for rec in slots[x][k].recs:
                    stream_end += rec.len
                    w.pending.append(Tab(tid, rec.key >> 1, rec.end, stream_end))
                slots[x][k].raw = List[UInt8]()
                slots[x][k].recs = List[BamRec]()
        if deflates > 0:
            w.addressed(dw % 2, deflates, w.text_base[dw % 2] + deflate_len)
            w.index_written()
        step_ns += s2 - s1
        serial_ns += (s1 - s0) + (Int(perf_counter_ns()) - s2)
    var last = (windows - 1) % 2 if windows > 0 else 0
    w.finish_tail(lib, last, end)
    w.index_written()
    w.close(lib)
    var t1 = Int(perf_counter_ns())
    if Bool(getenv("DROMOS_TIMING")):
        print(
            (
                t"dromos timing: sorting and writing {Float64(t1 - t0) * 1e-9} s (steps {Float64(step_ns) * 1e-9},"
                t" between steps {Float64(serial_ns) * 1e-9}; {staging.spills} spills,"
                t" {Float64(staging.spilled_bytes) * 1e-9} GB)"
            ),
            file=stderr,
        )


@fieldwise_init
struct _BamTask[lib_origin: ImmOrigin](ParallelTask):
    """One step of the final phase: `writes` (0 or 1) writes a window's deflated blocks in order, `gathers` read
    back, inflate and sort a window's buckets, `layouts` lay out the window before it, and `deflates` deflate the
    whole blocks of the window before that, BLOCKS_PER_TASK blocks each."""

    var staging: MutPointer[BamStaging, MutUntrackedOrigin]
    var lib: ImmPointer[OwnedDLHandle, Self.lib_origin]
    var level: Int
    var writes: Int
    var file: MutPointer[FileHandle, MutUntrackedOrigin]
    var write_outs: MutPointer[List[UInt8], MutUntrackedOrigin]
    var write_count: Int
    var gathers: Int
    var gather_slots: MutPointer[BamSlot, MutUntrackedOrigin]
    var layouts: Int
    var layout_slots: MutPointer[BamSlot, MutUntrackedOrigin]
    var layout_dst: MutPointer[UInt8, MutUntrackedOrigin]
    var deflates: Int
    var deflate_src: ImmPointer[UInt8, ImmUntrackedOrigin]
    var deflate_len: Int
    var deflate_outs: MutPointer[List[UInt8], MutUntrackedOrigin]
    var bad: MutPointer[Int, MutUntrackedOrigin]

    def run(self, index: Int):
        var i = index
        if i < self.writes:
            for t in range(self.write_count):
                self.file[].write_bytes(Span(self.write_outs[unsafe_offset=t]))
            return
        i -= self.writes
        if i < self.gathers:
            try:
                self.staging[].sort(self.gather_slots[unsafe_offset=i])
            except:
                self.bad[] = 1
            return
        i -= self.gathers
        if i < self.layouts:
            self.staging[].lay_out(self.layout_slots[unsafe_offset=i], self.layout_dst)
            return
        i -= self.layouts
        ref out = self.deflate_outs[unsafe_offset=i]
        out.clear()
        try:
            var d = Deflater(self.lib[], self.level)
            var at = i * BLOCKS_PER_TASK * BLOCK
            var stop = min(self.deflate_len, at + BLOCKS_PER_TASK * BLOCK)
            while at < stop:
                var n = min(BLOCK, stop - at)
                d.block(self.lib[], self.deflate_src.unsafe_offset(at), n, out)
                at += n
            d.free(self.lib[])
        except:
            self.bad[] = 2
