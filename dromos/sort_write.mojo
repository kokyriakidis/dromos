"""The final phase of an index: the staged buckets sorted and written as BGZF cut into blocks of exactly BLOCK bytes
(pgcat's IndexedWriter), the tabix index built as they go.

A window of buckets at a time, three deep: each step inflates and sorts window j, a bucket a task, lays out window
j - 1, deflates the whole blocks of window j - 2 and writes window j - 3, all at once. Block addresses, and the
index, follow from the deflated sizes, not from the writes."""

from std.ffi import OwnedDLHandle, external_call
from std.memory import unsafe_memcpy
from std.os import getenv
from std.sys import stderr
from std.time import perf_counter_ns

from dromos.bgzf import BLOCK, Deflater, eof_block, load_libdeflate
from dromos.parallel import ParallelTask, parallel_for
from dromos.staging import Rec, Slot, Staging
from dromos.stream import imm_untracked, mut_untracked
from dromos.tabix import TbiBuilder

comptime BLOCKS_PER_TASK = 16


def write_sorted(mut staging: Staging, max_span: List[Int], out_path: String, level: Int, threads: Int) raises:
    """Writes `out_path` and `out_path`.tbi (each through a temporary name) from the staged buckets: header lines
    `#max_span <contig> <bases>` (`max_span` by contig rank), then every line, sorted."""
    var t0 = Int(perf_counter_ns())
    var order = staging.order()
    var lib = load_libdeflate()
    var w = IndexWriter(out_path, level)
    # The longest alignment on each contig, for pgphase: a site whose every allele spans more graph bases cannot
    # be traversed whole by any read. Every contig is listed, 0 for one without alignments.
    var header = String()
    for r in range(len(staging.names)):
        header += String(t"#max_span\t{staging.names[r]}\t{max_span[r]}\n")
    w.header(header)
    # Windows of buckets, three deep: each step inflates and sorts window j, lays out window j - 1, deflates the
    # whole blocks of window j - 2 and writes window j - 3, all at once. Two slot sets, text buffers and output
    # sets take turns.
    var window = 4 * threads
    var windows = (len(order) + window - 1) // window
    var slots = List[List[Slot]]()
    for _ in range(2):
        var group = List[Slot]()
        for _ in range(window):
            group.append(Slot())
        slots.append(group^)
    # Each window's first stream offset (its lines'), then the end of the last; each window's buckets; how many
    # outputs each output set holds.
    var win_start = List[Int](length=windows + 1, fill=0)
    var win_count = List[Int](length=windows, fill=0)
    var deflated = List[Int](length=2, fill=0)
    var bad = List[Int](length=1, fill=0)
    var end = len(w.texts[0])
    var step_ns = 0
    var serial_ns = 0
    for j in range(windows + 3):
        var s0 = Int(perf_counter_ns())
        # Window j's buckets, to inflate and sort.
        var gathers = 0
        if j < windows:
            gathers = min(window, len(order) - j * window)
            win_count[j] = gathers
            for k in range(gathers):
                slots[j % 2][k].rank = order[j * window + k][0]
                slots[j % 2][k].bucket = order[j * window + k][1]
        # Window j - 1, sorted, laid out after what window j - 2 leaves of its last block.
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
        # Window j - 2's whole blocks, to deflate.
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
        # Window j - 3's blocks, to write.
        var wr = j - 3
        var writes = 1 if wr >= 0 and wr < windows else 0
        var task = _FinishTask(
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
            raise Error("a bucket's lines did not inflate")
        if bad[0] == 2:
            raise Error("BGZF: a block did not deflate")
        # Window j - 1's lines wait for their blocks; window j - 2's blocks have their addresses, and the lines
        # ending in them are indexed.
        if layouts > 0:
            var x = lay % 2
            for k in range(layouts):
                var contig = slots[x][k].rank
                var stream_end = w.text_base[x] + slots[x][k].at
                for rec in slots[x][k].recs:
                    stream_end += rec.out_len
                    w.pending.append(_Tab(contig, rec.start, rec.end, stream_end))
                slots[x][k].raw = List[UInt8]()
                slots[x][k].recs = List[Rec]()
        if deflates > 0:
            w.addressed(dw % 2, deflates, w.text_base[dw % 2] + deflate_len)
            w.index_written(staging.names)
        step_ns += s2 - s1
        serial_ns += (s1 - s0) + (Int(perf_counter_ns()) - s2)
    # The last window's partial block (or the header alone), then the end.
    var last = (windows - 1) % 2 if windows > 0 else 0
    w.finish_tail(lib, last, end)
    w.index_written(staging.names)
    w.close(lib)
    var t1 = Int(perf_counter_ns())
    if Bool(getenv("DROMOS_TIMING")):
        print(
            (
                t"dromos timing: sorting and writing {Float64(t1 - t0) * 1e-9} s (steps {Float64(step_ns) * 1e-9},"
                t" between steps {Float64(serial_ns) * 1e-9})"
            ),
            file=stderr,
        )


@fieldwise_init
struct _Tab(TrivialRegisterPassable):
    """A line written, waiting for the block holding its end: its contig, interval, and its end in the stream."""

    var contig: Int
    var beg: Int
    var end: Int
    var stream_end: Int


struct IndexWriter(Movable):
    """The BGZF output, cut into blocks of exactly BLOCK bytes, so the virtual offset of any byte follows from its
    stream offset and the blocks' addresses (pgcat's IndexedWriter): lines are indexed once the block holding their
    end has its address."""

    var path: String
    var tmp: String
    var file: FileHandle
    var level: Int
    var texts: List[List[UInt8]]
    """Two windows' text, taking turns: the one being laid out, the one being deflated."""
    var text_base: List[Int]
    """The stream offset of each text's first byte, a block boundary."""
    var outs: List[List[List[UInt8]]]
    """Two windows' deflated blocks, taking turns: the ones being deflated, the ones being written."""
    var base: Int
    """The stream bytes in blocks with addresses."""
    var written: Int
    """The compressed bytes of those blocks: the address of the next block."""
    var addrs: List[Int]
    """Every block's address, by block number."""
    var pending: List[_Tab]
    var head: Int
    var index: TbiBuilder

    def __init__(out self, path: String, level: Int) raises:
        self.path = path
        self.tmp = path + ".tmp"
        self.file = open(self.tmp, "w")
        self.level = level
        self.texts = List[List[UInt8]](length=2, fill=List[UInt8]())
        self.text_base = List[Int](length=2, fill=0)
        self.outs = List[List[List[UInt8]]](length=2, fill=List[List[UInt8]]())
        self.base = 0
        self.written = 0
        self.addrs = List[Int]()
        self.pending = List[_Tab]()
        self.head = 0
        self.index = TbiBuilder()

    def header(mut self, text: String) raises:
        """Header lines (each starting with '#', with its newline), before any line: tabix skips them, and the index
        starts after them."""
        self.texts[0].extend(text.as_bytes())
        if len(self.texts[0]) >= BLOCK:
            raise Error("header longer than a BGZF block")
        # Block 0 starts at address 0: the first line starts at (0, length).
        self.index.set_first_offset(len(self.texts[0]))

    def addressed(mut self, x: Int, outputs: Int, through: Int):
        """The blocks in output set `x`'s first `outputs` outputs, deflated from the stream up to `through`, get
        their addresses."""
        for t in range(outputs):
            var p = self.outs[x][t].unsafe_ptr().as_imm()
            var at = 0
            while at < len(self.outs[x][t]):
                self.addrs.append(self.written)
                var size = (Int(p[unsafe_offset=at + 16]) | Int(p[unsafe_offset=at + 17]) << 8) + 1
                self.written += size
                at += size
        self.base = through

    def finish_tail(mut self, lib: OwnedDLHandle, x: Int, end: Int) raises:
        """Deflates and writes what text `x` holds past the last whole block, up to stream offset `end`."""
        var rest = end - self.base
        if rest <= 0:
            return
        self.outs[x].clear()
        self.outs[x].append(List[UInt8]())
        var d = Deflater(lib, self.level)
        d.block(lib, self.texts[x].unsafe_ptr().unsafe_offset(self.base - self.text_base[x]), rest, self.outs[x][0])
        d.free(lib)
        self.file.write_bytes(Span(self.outs[x][0]))
        self.addressed(x, 1, end)

    def _voffset(self, s: Int) -> Int:
        """The virtual offset of stream offset `s`, in a block with an address or at the start of the next: a
        block's end reads as the start of the next block, as htslib reports it after reading a line."""
        if s == self.base:
            return self.written << 16
        return self.addrs[s // BLOCK] << 16 | (s % BLOCK)

    def index_written(mut self, contigs: List[String]) raises:
        """Indexes the lines whose ends lie in blocks with addresses."""
        while self.head < len(self.pending):
            var t = self.pending[self.head]
            if t.stream_end > self.base:
                break
            self.index.push(t.contig, contigs[t.contig], t.beg, t.end, self._voffset(t.stream_end))
            self.head += 1
        if self.head == len(self.pending):
            self.pending.clear()
            self.head = 0

    def close(mut self, lib: OwnedDLHandle) raises:
        """The EOF block, the index, and both files under their names."""
        var final_offset = self.written << 16
        var eof = List[UInt8]()
        eof_block(eof)
        self.file.write_bytes(Span(eof))
        self.file.close()
        self.index.finish(final_offset)
        var tbi = self.path + ".tbi"
        self.index.write(lib, tbi + ".tmp")
        _rename(self.tmp, self.path)
        _rename(tbi + ".tmp", tbi)


def _rename(src: String, dst: String) raises:
    var a = List[UInt8](src.as_bytes())
    a.append(0)
    var b = List[UInt8](dst.as_bytes())
    b.append(0)
    if external_call["rename", Int32](a.unsafe_ptr(), b.unsafe_ptr()) != 0:
        raise Error(t"cannot rename {src} to {dst}")


@fieldwise_init
struct _FinishTask[lib_origin: ImmOrigin](ParallelTask):
    """One step of the final phase: `writes` (0 or 1) writes a window's deflated blocks in order, `gathers` inflate
    and sort a window's buckets, `layouts` lay out the window before it, and `deflates` deflate the whole blocks of
    the window before that, BLOCKS_PER_TASK blocks each."""

    var staging: MutPointer[Staging, MutUntrackedOrigin]
    var lib: ImmPointer[OwnedDLHandle, Self.lib_origin]
    var level: Int
    var writes: Int
    var file: MutPointer[FileHandle, MutUntrackedOrigin]
    var write_outs: MutPointer[List[UInt8], MutUntrackedOrigin]
    var write_count: Int
    var gathers: Int
    var gather_slots: MutPointer[Slot, MutUntrackedOrigin]
    var layouts: Int
    var layout_slots: MutPointer[Slot, MutUntrackedOrigin]
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
