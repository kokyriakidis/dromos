"""The `.pgsites` site catalog (pgcatalog's pgsites.rs): VCF lines in independently zstd-compressed position blocks,
so a region query decompresses only the blocks it needs. Written by `dromos catalog` (`Writer`, byte for byte as
pgcatalog writes it), printed by `dromos view` (`view`), and read for its snarl windows by `dromos surject --catalog`
(`load_windows`, as pgcatalog pileup's `load_windows` makes them: each record whose ID names its snarl's two boundary
nodes, both visited once by the contig, gives the reference bases between them, at most MAX_SITE; windows inside
another are dropped, the rest padded by WINDOW_PAD and those that overlap merged up to MAX_WINDOW).

Layout (integers little-endian): "PGSITES2"; u64 header length and the zstd VCF header; zstd blocks of VCF lines;
the index (u64 count, then per block nine u64s: first_pos, last_pos, max_end, offset, compressed length, raw length,
records, min_allele, max_walk); u64 the index's offset.

A block holds the records whose POS falls in one BLOCK_BP window, split every BLOCK_BYTES of text; a record over
OWN_BLOCK_BYTES has a block of its own. max_end is the largest POS + len(REF) - 1, the span tabix uses for a VCF
record, so a query for [beg, end] reads the blocks with first_pos <= end and max_end >= beg and keeps the records
that overlap. min_allele is the fewest bases any allele of the block's records spans (padding base not counted);
max_walk the most nodes any allele walk of its records visits. Blocks are compressed with long-distance matching and
a window as long as the block (log <= 31); one over MT_BYTES with zstd's own threads (whose output does not depend on
how many)."""

from std.ffi import OwnedDLHandle, external_call
from std.memory import unsafe_memcpy, unsafe_memmove

from dromos.bgzf import MappedFile, load_zstd, zstd_decompress
from dromos.contig import Contig, Reference
from dromos.gaf_text import parse_unsigned
from dromos.linearize import Window
from dromos.output import OutFile, write_all
from dromos.parallel import ParallelTask, parallel_for
from dromos.refmap import RefMap
from dromos.stream import imm_untracked, mut_untracked, untracked

comptime MAGIC = "PGSITES2"
comptime INDEX_FIELDS = 9
comptime WINDOW_BUFFER = 8 << 20
comptime _ZSTD_D_WINDOW_LOG_MAX = 100
"""The bytes of a block surject's windows are read from at a time."""
comptime BLOCK_BP = 100_000
comptime BLOCK_BYTES = 4 << 20
comptime OWN_BLOCK_BYTES = 1 << 20
comptime MT_BYTES = 16 << 20
"""Blocks longer than this are compressed with zstd's own threads (a giant record's hundreds of near-identical
alleles, which only a window as long as the block finds)."""
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


@fieldwise_init
struct _LoadWindows[map_origin: ImmOrigin, ref_origin: ImmOrigin](ParallelTask):
    """Catalog `i`'s windows on contig `i` of the reference, a contig a task."""

    var map: ImmPointer[RefMap, Self.map_origin]
    var reference: ImmPointer[Reference, Self.ref_origin]
    var paths: ImmPointer[String, ImmUntrackedOrigin]
    var out: MutPointer[List[Window], MutUntrackedOrigin]
    var failed: MutPointer[Int, MutUntrackedOrigin]

    def run(self, index: Int):
        try:
            self.out[unsafe_offset=index] = load_windows(
                self.paths[unsafe_offset=index], self.reference[].contigs[index], self.map[]
            )
        except:
            self.failed[] = index + 1


def load_all_windows(paths: List[String], reference: Reference, map: RefMap, threads: Int) raises -> List[List[Window]]:
    """Each contig's windows from its catalog (`paths[i]` for contig i), contigs read at once."""
    var out = List[List[Window]](length=len(paths), fill=List[Window]())
    var failed = List[Int](length=1, fill=0)
    var task = _LoadWindows(
        Pointer(to=map),
        Pointer(to=reference),
        paths.unsafe_ptr().as_imm().unsafe_origin_cast[ImmUntrackedOrigin](),
        out.unsafe_ptr().unsafe_origin_cast[MutUntrackedOrigin](),
        failed.unsafe_ptr().unsafe_origin_cast[MutUntrackedOrigin](),
    )
    parallel_for(task, len(paths), threads)
    if failed[0] > 0:
        # Read again on this thread, for its error.
        var k = failed[0] - 1
        _ = load_windows(paths[k], reference.contigs[k], map)
        raise Error(t"{paths[k]}: its windows could not be read")
    return out^


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
    # Each block streamed through a buffer of WINDOW_BUFFER bytes: its whole lines read there, and of a line longer
    # than the buffer (a giant record) the ID at its start, the rest passed over.
    var dctx = lib.get_function[Int]("ZSTD_createDCtx")()
    if dctx == 0:
        raise Error("zstd: no decompression context")
    # Catalog frames may use long-distance windows past the streaming default (2^27): allow zstd's largest.
    _ = lib.get_function[Int]("ZSTD_DCtx_setParameter")(dctx, Int32(_ZSTD_D_WINDOW_LOG_MAX), Int32(31))
    var buffer = List[UInt8](length=WINDOW_BUFFER, fill=0)
    var q = buffer.unsafe_ptr().as_imm()
    var input = List[Int](length=3, fill=0)
    var output = List[Int](length=3, fill=0)
    var stream = lib.get_function[Int]("ZSTD_decompressStream")
    var is_error = lib.get_function[UInt32]("ZSTD_isError")
    var reset = lib.get_function[Int]("ZSTD_DCtx_reset")
    try:
        for b in range(count):
            var e = index_at + 8 + b * INDEX_FIELDS * 8
            var offset = _u64(p, e + 3 * 8)
            var compressed = _u64(p, e + 4 * 8)
            var raw_length = _u64(p, e + 5 * 8)
            if is_error(reset(dctx, Int32(1))) != 0:
                raise Error("zstd: a context did not reset")
            input[0] = Int(p.unsafe_offset(offset))
            input[1] = compressed
            input[2] = 0
            var have = 0
            var skipping = False
            var produced = 0
            while True:
                output[0] = Int(buffer.unsafe_ptr())
                output[1] = WINDOW_BUFFER
                output[2] = have
                var left = stream(dctx, output.unsafe_ptr(), input.unsafe_ptr())
                if is_error(left) != 0:
                    raise Error(t"{path}: a block does not inflate (zstd error {-left})")
                var filled = output[2]
                produced += filled - have
                var at = 0
                while at < filled:
                    var hit = external_call["memchr", Int](q.unsafe_offset(at), Int32(ord("\n")), filled - at)
                    if hit == 0:
                        break
                    var end = hit - Int(q)
                    if skipping:
                        skipping = False
                    elif end > at:
                        _snarl(contig, map, q, at, end, all)
                    at = end + 1
                var tail = filled - at
                if skipping:
                    tail = 0
                elif tail == WINDOW_BUFFER:
                    # A line longer than the buffer: its ID is in the buffered start; the rest is passed over.
                    _snarl(contig, map, q, 0, tail, all)
                    skipping = True
                    tail = 0
                if tail > 0 and at > 0:
                    var dst = mut_untracked(buffer.unsafe_ptr())
                    var src = imm_untracked(buffer.unsafe_ptr()).unsafe_offset(at)
                    unsafe_memmove(dest=dst, src=src, count=tail)
                var grew = filled > have
                have = tail
                if left == 0 and input[2] == input[1]:
                    break
                if not grew and input[2] == input[1]:
                    raise Error(t"{path}: a block ends early")
            # The block's last line, without a newline.
            if have > 0 and not skipping:
                _snarl(contig, map, q, 0, have, all)
            if produced != raw_length:
                raise Error(t"{path}: a block's size differs from its index")
    finally:
        _ = lib.get_function[Int]("ZSTD_freeDCtx")(dctx)
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


# Writing (pgsites.rs's Writer).

comptime _C_LEVEL: Int32 = 100
comptime _C_WINDOW_LOG: Int32 = 101
comptime _C_LDM: Int32 = 160
comptime _C_CONTENT_SIZE: Int32 = 200
comptime _C_WORKERS: Int32 = 400
comptime _RESET_ALL: Int32 = 3
comptime _E_CONTINUE: Int32 = 0
comptime _E_END: Int32 = 2


def window_log(length: Int) -> Int:
    """The window a block of `length` bytes is compressed with: as long as the block, from 2^20 to 2^31."""
    var x = max(length, 2) - 1
    var bits = 0
    while x > 0:
        bits += 1
        x >>= 1
    return min(max(bits, 20), 31)


struct Compressor(Movable):
    """A zstd compression context, set up for each frame as pgcatalog's zstd crate sets up its own: the level, the
    content size recorded, and for a block long-distance matching, its window and, over MT_BYTES, `threads`
    workers."""

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

    def start(self, level: Int, block_length: Int, threads: Int, block: Bool) raises:
        """A new frame's parameters: for a block of `block_length` bytes when `block`, else the header's."""
        self._check(self.lib.get_function[Int]("ZSTD_CCtx_reset")(self.cctx, _RESET_ALL))
        self._set(_C_LEVEL, level)
        if block:
            self._set(_C_CONTENT_SIZE, 1)
            self._set(_C_LDM, 1)
            self._set(_C_WINDOW_LOG, window_log(block_length))
            if block_length > MT_BYTES:
                self._set(_C_WORKERS, threads)

    def frame(self, data: ImmPointer[UInt8, _], length: Int) raises -> List[UInt8]:
        """`data[:length]` as one frame, with the parameters `start` set."""
        var bound = self.lib.get_function[Int]("ZSTD_compressBound")(length)
        var out = List[UInt8](unsafe_uninit_length=bound)
        var n = self.lib.get_function[Int]("ZSTD_compress2")(self.cctx, out.unsafe_ptr(), bound, data, length)
        self._check(n)
        # A copy of the frame's size: the bound's capacity would stay with every frame kept.
        var frame = List[UInt8](unsafe_uninit_length=n)
        unsafe_memcpy(dest=frame.unsafe_ptr(), src=out.unsafe_ptr(), count=n)
        return frame^


@fieldwise_init
struct _ZBuffer(TrivialRegisterPassable):
    """ZSTD_inBuffer and ZSTD_outBuffer: memory, its size, how far zstd has got."""

    var ptr: Int
    var size: Int
    var pos: Int


comptime STREAM_PIECE = 1 << 20
"""The bytes a streamed record hands zstd at a time."""


struct FrameStream(Movable):
    """One record's block compressed as its line is written (pgsites.rs's compress_record): the line never exists in
    memory whole. `put` its bytes, then `finish`, which checks it got the `length` it was promised."""

    var compressor: Compressor
    var length: Int
    var written: Int
    var piece: List[UInt8]
    var out: List[UInt8]

    def __init__(out self, length: Int, level: Int, threads: Int) raises:
        """For a line of `length` bytes, newline included."""
        self.compressor = Compressor()
        self.compressor.start(level, length, threads, True)
        self.compressor._check(
            self.compressor.lib.get_function[Int]("ZSTD_CCtx_setPledgedSrcSize")(self.compressor.cctx, UInt64(length))
        )
        self.length = length
        self.written = 0
        self.piece = List[UInt8](capacity=STREAM_PIECE)
        self.out = List[UInt8]()

    def _compress(mut self, data: Int, length: Int, mode: Int32) raises:
        var src = _ZBuffer(data, length, 0)
        var call = self.compressor.lib.get_function[Int]("ZSTD_compressStream2")
        while True:
            var at = len(self.out)
            var room = max(STREAM_PIECE, at // 4)
            self.out.resize(unsafe_uninit_length=at + room)
            var dst = _ZBuffer(Int(self.out.unsafe_ptr()), at + room, at)
            var left = call(self.compressor.cctx, Pointer(to=dst), Pointer(to=src), mode)
            self.compressor._check(left)
            self.out.resize(unsafe_uninit_length=dst.pos)
            if mode == _E_END:
                if left == 0:
                    break
            elif src.pos == src.size:
                break
        self.written += length

    def _run(mut self, mode: Int32) raises:
        self._compress(Int(self.piece.unsafe_ptr()), len(self.piece), mode)
        self.piece.clear()

    def feed(mut self, data: ImmPointer[UInt8, _], length: Int) raises:
        """`length` bytes at `data`, handed to zstd as they are (after the piece so far)."""
        if len(self.piece) > 0:
            self._run(_E_CONTINUE)
        self._compress(Int(data), length, _E_CONTINUE)

    @always_inline
    def room(self) -> Int:
        return STREAM_PIECE - len(self.piece)

    def end_ptr(mut self) -> MutPointer[UInt8, MutUntrackedOrigin]:
        """Where the piece's next byte goes (`room` of them, then `grew`)."""
        return mut_untracked(self.piece.unsafe_ptr().unsafe_offset(len(self.piece)))

    def grew(mut self, n: Int) raises:
        self.piece.resize(unsafe_uninit_length=len(self.piece) + n)
        if len(self.piece) == STREAM_PIECE:
            self._run(_E_CONTINUE)

    def finish(mut self) raises -> List[UInt8]:
        self._run(_E_END)
        if self.written != self.length:
            raise Error(t"record wrote {self.written} bytes, expected {self.length}")
        var out = self.out^
        self.out = List[UInt8]()
        return out^


struct Unit(Movable):
    """A block of the catalog as it is written: lines placed in its text, then compressed into its frame. A giant
    unit is one record compressed from its walks elsewhere (FrameStream), its frame set from there."""

    var window: Int
    var raw: List[UInt8]
    var frame: List[UInt8]
    var first: Int
    var last: Int
    var max_end: Int
    var records: Int
    var min_allele: Int
    var max_walk: Int
    var raw_length: Int
    """A giant's line length, newline included."""
    var sealed: Bool
    """Takes no further records: an own block (a record over OWN_BLOCK_BYTES), a giant, or followed by another."""
    var giant: Bool
    var done: Bool
    """Its frame is made."""

    def __init__(out self, window: Int, pos: Int, sealed: Bool):
        self.window = window
        self.raw = List[UInt8]()
        self.frame = List[UInt8]()
        self.first = pos
        self.last = pos
        self.max_end = 0
        self.records = 0
        self.min_allele = Int.MAX
        self.max_walk = 0
        self.raw_length = 0
        self.sealed = sealed
        self.giant = False
        self.done = False

    def compress(mut self, c: Compressor, level: Int, threads: Int) raises:
        """The block's text into its frame (one zstd frame, `c` set up afresh for it), the text then freed."""
        c.start(level, len(self.raw), threads, True)
        self.frame = c.frame(self.raw.unsafe_ptr(), len(self.raw))
        self.raw_length = len(self.raw)
        self.raw = List[UInt8]()
        self.done = True


struct Writer(Movable):
    """Writes `<path>.tmp` and renames it to `path` once complete, so a catalog is never seen half written (an
    interrupted build leaves only the .tmp). Records arrive in POS order: `place` gives each line its room in a
    block (a unit), or `place_giant` a unit of its own; the caller writes the lines there and compresses the sealed
    units (any thread, each its own), and `write_done` writes the units done, in order. Units are numbered from 0
    in file order; `units` holds those not yet written, from `base`."""

    var path: String
    var file: OutFile
    var offset: Int
    var level: Int
    var threads: Int
    var units: List[Unit]
    var base: Int
    var index: List[UInt64]

    def __init__(out self, path: String, header: String, level: Int, threads: Int) raises:
        self.path = path
        self.file = OutFile(path + ".tmp")
        self.level = level
        self.threads = threads
        self.units = List[Unit]()
        self.base = 0
        self.index = List[UInt64]()
        var c = Compressor()
        c.start(level, 0, threads, False)
        var h = c.frame(header.unsafe_ptr(), header.byte_length())
        var head = List[UInt8](MAGIC.as_bytes())
        _put_u64(head, len(h))
        head.extend(Span(h))
        self.file.put(Span(head), 0)
        self.offset = len(head)

    def _seal_last(mut self):
        if len(self.units) > 0:
            self.units[len(self.units) - 1].sealed = True

    def place(mut self, pos: Int, max_end: Int, min_allele: Int, max_walk: Int, length: Int) -> Tuple[Int, Int]:
        """Room for a record's line of `length` bytes (newline excluded, written here): its unit and offset there.
        A record starts a new block when it has one of its own, or the block filling is sealed, of another
        BLOCK_BP window, or BLOCK_BYTES long."""
        var window = pos // BLOCK_BP
        var own = length > OWN_BLOCK_BYTES
        var n = len(self.units)
        var new_block = own or n == 0
        if not new_block:
            ref b = self.units[n - 1]
            new_block = b.sealed or b.window != window or len(b.raw) >= BLOCK_BYTES
        if new_block:
            self._seal_last()
            self.units.append(Unit(window, pos, own))
            # Room for the block's lines at once: it takes lines until it is BLOCK_BYTES long, the last at most
            # OWN_BLOCK_BYTES and its newline.
            self.units[len(self.units) - 1].raw.reserve(length + 1 if own else BLOCK_BYTES + OWN_BLOCK_BYTES + 1)
        var k = len(self.units) - 1
        ref b = self.units[k]
        var at = len(b.raw)
        b.raw.resize(unsafe_uninit_length=at + length + 1)
        b.raw[at + length] = UInt8(ord("\n"))
        b.last = pos
        b.max_end = max(b.max_end, max_end)
        b.min_allele = min(b.min_allele, min_allele)
        b.max_walk = max(b.max_walk, max_walk)
        b.records += 1
        return (self.base + k, at)

    def place_giant(mut self, pos: Int, max_end: Int, min_allele: Int, max_walk: Int, raw_length: Int) -> Int:
        """A unit for one record compressed elsewhere, after every unit so far (which it seals)."""
        self._seal_last()
        var u = Unit(pos // BLOCK_BP, pos, True)
        u.giant = True
        u.max_end = max_end
        u.min_allele = min_allele
        u.max_walk = max_walk
        u.records = 1
        u.raw_length = raw_length
        self.units.append(u^)
        return self.base + len(self.units) - 1

    def unit(mut self, id: Int) -> MutPointer[Unit, MutUntrackedOrigin]:
        """Unit `id`: valid until the next `place`, `place_giant` or `write_done`."""
        return mut_untracked(self.units.unsafe_ptr().unsafe_offset(id - self.base))

    def at(mut self, id: Int, offset: Int) -> MutPointer[UInt8, MutUntrackedOrigin]:
        """Where `place` put a line: valid until the next `place`."""
        return mut_untracked(self.units[id - self.base].raw.unsafe_ptr().unsafe_offset(offset))

    def seal(mut self):
        """No further records: the block filling is sealed too."""
        self._seal_last()

    def write_done(mut self) raises:
        """Writes the units done, up to the first not done, in order."""
        var n = 0
        while n < len(self.units) and self.units[n].done:
            n += 1
        if n == 0:
            return
        var out = List[UInt8]()
        for i in range(n):
            ref b = self.units[i]
            var entry: List[Int] = [
                b.first,
                b.last,
                b.max_end,
                self.offset + len(out),
                len(b.frame),
                b.raw_length,
                b.records,
                b.min_allele,
                b.max_walk,
            ]
            for e in entry:
                self.index.append(UInt64(e))
            out.extend(Span(b.frame))
        self.file.put(Span(out), self.offset)
        self.offset += len(out)
        var rest = List[Unit]()
        while len(self.units) > n:
            rest.append(self.units.pop(n))
        self.units = rest^
        self.base += n

    def blocks(self) -> Int:
        return len(self.index) // INDEX_FIELDS

    def finish(mut self) raises -> Int:
        """Writes the index (every unit must be written), syncs the file and renames it; returns the block count."""
        if len(self.units) > 0:
            raise Error(t"{self.path}: blocks left unwritten")
        var tail = List[UInt8]()
        _put_u64(tail, self.blocks())
        for e in self.index:
            _put_u64(tail, Int(e))
        _put_u64(tail, self.offset)
        self.file.put(Span(tail), self.offset)
        if external_call["fsync", Int32](self.file.fd) != 0:
            raise Error(t"cannot write {self.path}.tmp")
        self.file.close()
        var a = List[UInt8]((self.path + ".tmp").as_bytes())
        a.append(0)
        var b = List[UInt8](self.path.as_bytes())
        b.append(0)
        if external_call["rename", Int32](a.unsafe_ptr(), b.unsafe_ptr()) != 0:
            raise Error(t"cannot rename {self.path}.tmp to {self.path}")
        return self.blocks()


@fieldwise_init
struct _CompressTask(ParallelTask):
    """Unit `first + index` compressed, if sealed and not done."""

    var units: MutPointer[Unit, MutUntrackedOrigin]
    var level: Int
    var threads: Int
    var bad: MutPointer[Int, MutUntrackedOrigin]

    def run(self, index: Int):
        ref u = self.units[unsafe_offset=index]
        if u.sealed and not u.done and not u.giant:
            try:
                u.compress(Compressor(), self.level, self.threads)
            except:
                self.bad[] = 1


def compress_sealed(mut w: Writer) raises:
    """Compresses every sealed block not done, in parallel (for a writer filled on one thread)."""
    var bad = 0
    var task = _CompressTask(mut_untracked(w.units.unsafe_ptr()), w.level, w.threads, mut_untracked(Pointer(to=bad)))
    parallel_for(task, len(w.units), w.threads)
    if bad != 0:
        raise Error(t"{w.path}: zstd compression failed")


@always_inline
def _put_u64(mut out: List[UInt8], value: Int):
    for i in range(8):
        out.append(UInt8((value >> (8 * i)) & 0xFF))


# Reading (pgsites.rs's Reader) and `dromos view`.


struct Catalog(Movable):
    """A catalog mapped into memory: its header, its CHROM (from the header's ##contig line) and its blocks."""

    var file: MappedFile
    var header: String
    var chrom: String
    var blocks: List[Int]
    """Per block: first_pos, max_end, offset, compressed length, raw length."""

    def __init__(out self, path: String) raises:
        self.file = MappedFile(path)
        var p = self.file.ptr
        var n = self.file.size
        var magic = MAGIC.as_bytes()
        var prefix = n >= 7
        var has_magic = n >= 8
        for i in range(min(8, n)):
            if i < 7:
                prefix = prefix and p[unsafe_offset=i] == magic[i]
            has_magic = has_magic and p[unsafe_offset=i] == magic[i]
        if prefix and not has_magic:
            raise Error(t"{path}: an older .pgsites format; rebuild it with dromos catalog")
        if not has_magic or n < 24:
            raise Error(t"{path}: not a .pgsites file")
        var lib = load_zstd()
        var header_length = _u64(p, 8)
        if 16 + header_length > n:
            raise Error(t"{path}: not a .pgsites file")
        var text = _inflate(lib, Span(unsafe_ptr=p.unsafe_offset(16), length=header_length))
        self.header = String(unsafe_from_utf8=text)
        var at = self.header.find("##contig=<ID=")
        if at < 0:
            raise Error(t"{path}: no ##contig line")
        at += 13
        var rest = self.header.as_bytes()
        var stop = at
        while stop < len(rest) and rest[stop] != UInt8(ord(",")) and rest[stop] != UInt8(ord(">")):
            stop += 1
        # pgcatalog's find(...).unwrap_or(0): no , or > after the ID leaves it empty.
        self.chrom = String(self.header[byte=at:stop]) if stop < len(rest) else String()
        var index_at = _u64(p, n - 8)
        if index_at + 8 > n:
            raise Error(t"{path}: a damaged index")
        var count = _u64(p, index_at)
        if index_at + 8 + count * INDEX_FIELDS * 8 > n:
            raise Error(t"{path}: a damaged index")
        self.blocks = List[Int](capacity=5 * count)
        for k in range(count):
            var e = index_at + 8 + k * INDEX_FIELDS * 8
            self.blocks.append(_u64(p, e))
            self.blocks.append(_u64(p, e + 2 * 8))
            self.blocks.append(_u64(p, e + 3 * 8))
            self.blocks.append(_u64(p, e + 4 * 8))
            self.blocks.append(_u64(p, e + 5 * 8))
            if self.blocks[5 * k + 2] + self.blocks[5 * k + 3] > n:
                raise Error(t"{path}: a damaged index")


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


@fieldwise_init
struct _ViewTask(ParallelTask):
    """Block `order[first + index]` inflated and its lines overlapping [beg, end] kept, into `texts[index]`."""

    var catalog: ImmPointer[Catalog, ImmUntrackedOrigin]
    var order: ImmPointer[Int, ImmUntrackedOrigin]
    var first: Int
    var beg: Int
    var end: Int
    var texts: MutPointer[List[UInt8], MutUntrackedOrigin]
    var bad: MutPointer[Int, MutUntrackedOrigin]

    def run(self, index: Int):
        ref c = self.catalog[]
        var b = self.order[unsafe_offset=self.first + index]
        var raw_length = c.blocks[5 * b + 4]
        var raw = List[UInt8](unsafe_uninit_length=raw_length)
        try:
            var lib = load_zstd()
            var frame = Span(unsafe_ptr=c.file.ptr.unsafe_offset(c.blocks[5 * b + 2]), length=c.blocks[5 * b + 3])
            zstd_decompress(lib, frame, raw.unsafe_ptr(), raw_length)
        except:
            self.bad[] = 1
            return
        # Lines kept in place: the block is its own output when every line overlaps.
        var q = raw.unsafe_ptr()
        var kept = 0
        var at = 0
        while at < raw_length:
            var hit = external_call["memchr", Int](q.unsafe_offset(at), Int32(ord("\n")), raw_length - at)
            var stop = raw_length if hit == 0 else hit - Int(q)
            if stop > at and _overlaps(q.as_imm(), at, stop, self.beg, self.end):
                if kept != at:
                    _ = external_call["memmove", Int](Int(q) + kept, Int(q) + at, stop - at)
                kept += stop - at
                q[unsafe_offset=kept] = UInt8(ord("\n"))
                kept += 1
            at = stop + 1
        raw.resize(unsafe_uninit_length=kept)
        self.texts[unsafe_offset=index] = raw^


def _overlaps(p: ImmPointer[UInt8, _], at: Int, stop: Int, beg: Int, end: Int) -> Bool:
    """Whether the VCF line `p[at:stop]` overlaps [beg, end] as tabix sees it: POS ..= POS + len(REF) - 1 (a POS
    that is not a number read as 0, a missing REF as one base)."""
    var tabs = Array[Int, 4](fill=stop)
    var found = 0
    var k = at
    while k < stop and found < 4:
        if p[unsafe_offset=k] == UInt8(ord("\t")):
            tabs[found] = k
            found += 1
        k += 1
    var pos = 0
    if found >= 1:
        pos = max(parse_unsigned(p, tabs[0] + 1, tabs[1]), 0)
    var ref_length = 1
    if found >= 3:
        ref_length = max(tabs[3] - tabs[2] - 1, 1)
    return pos <= end and pos + ref_length > beg


comptime VIEW_BYTES = 1 << 30
"""The most raw bytes a step of `view` inflates at once."""


def view(path: String, region: Optional[String], threads: Int) raises:
    """`pgcatalog view`: the header, then every record (or those overlapping `region`, CONTIG[:START-END], as tabix
    would return them), on stdout. Blocks are inflated in parallel, a step's worth at a time, and written in order."""
    var c = Catalog(path)
    var beg = 1
    var end = Int.MAX
    var contig = Optional[String]()
    if region:
        var q = region.value()
        var colon = q.rfind(":")
        if colon < 0:
            contig = q
        else:
            contig = String(q[byte=:colon])
            var span = String(q[byte = colon + 1 :])
            var dash = span.find("-")
            if dash < 0:
                raise Error("region: CONTIG:START-END")
            beg = _number(String(span[byte=:dash]), "region start")
            end = _number(String(span[byte = dash + 1 :]), "region end")
    write_all(1, c.header.unsafe_ptr(), c.header.byte_length())
    if contig and not same_contig(contig.value(), c.chrom):
        return
    var none = List[UInt8]()
    records(c, beg, end, threads, 1, none)
    _ = c^


def records(c: Catalog, beg: Int, end: Int, threads: Int, fd: Int32, mut kept: List[UInt8]) raises:
    """The record lines overlapping [beg, end] (1-based, inclusive), in order, as tabix would return them from the VCF:
    written to `fd`, or (`fd` -1) appended to `kept`. The blocks that can hold them are inflated in parallel, up to
    VIEW_BYTES at a time."""
    var order = List[Int]()
    for b in range(len(c.blocks) // 5):
        if not (c.blocks[5 * b] > end or c.blocks[5 * b + 1] < beg):
            order.append(b)
    var first = 0
    while first < len(order):
        var last = first
        var bytes = 0
        while last < len(order) and (last == first or bytes + c.blocks[5 * order[last] + 4] <= VIEW_BYTES):
            bytes += c.blocks[5 * order[last] + 4]
            last += 1
        var texts = List[List[UInt8]](length=last - first, fill=List[UInt8]())
        var bad = 0
        var task = _ViewTask(
            untracked(Pointer(to=c)),
            imm_untracked(order.unsafe_ptr()),
            first,
            beg,
            end,
            mut_untracked(texts.unsafe_ptr()),
            mut_untracked(Pointer(to=bad)),
        )
        parallel_for(task, last - first, threads)
        _ = len(order)
        if bad != 0:
            raise Error("a catalog block did not decompress")
        for t in texts:
            if fd >= 0:
                write_all(fd, t.unsafe_ptr(), len(t))
            else:
                kept.extend(Span(t))
        first = last


def _number(text: String, what: String) raises -> Int:
    """A region bound as pgcatalog reads it: commas dropped, then an unsigned integer."""
    var digits = List[UInt8]()
    for b in text.as_bytes():
        if b != UInt8(ord(",")):
            digits.append(b)
    var v = parse_unsigned(digits.unsafe_ptr(), 0, len(digits))
    if v < 0:
        raise Error(what)
    return v
