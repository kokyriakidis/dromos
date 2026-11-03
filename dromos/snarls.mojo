"""Snarls files written by `vg snarls` (pgcatalog's snarls.rs): a libvgio stream of `vg.Snarl` messages, read for each
snarl's start and end visits and its parent's.

The stream is plain, gzip or BGZF (a contig's file is its components' files one after another) and is a sequence of
groups, each a varint count, a varint-long type tag ("SNARL") and that many minus one varint-long messages. A visit
is a u64 as the GBWT encodes it, `node << 1 | reverse`; compared as numbers, visits order as pgcatalog's (node,
orientation) handles do.
"""

from std.ffi import OwnedDLHandle

from dromos.bgzf import GzipStream, RawFile, inflate_blocks, load_libdeflate, load_zlib, scan_blocks
from dromos.parallel import ParallelTask, parallel_for
from dromos.stream import BGZF, GZIP, encoding, imm_untracked, mut_untracked, untracked

comptime NO_PARENT = UInt64.MAX
"""`Snarl.parent_start` of a snarl the file gives no parent."""
comptime SNARL_TAG = "SNARL"
comptime DECODE_SPAN = 1 << 16
"""Messages one task decodes."""


@fieldwise_init
struct Snarl(Copyable, TrivialRegisterPassable):
    """A snarl's start and end visits, and its parent's as the file records them (NO_PARENT for none)."""

    var start: UInt64
    var end: UInt64
    var parent_start: UInt64
    var parent_end: UInt64

    @always_inline
    def has_parent(self) -> Bool:
        return self.parent_start != NO_PARENT


@always_inline
def canonical(start: UInt64, end: UInt64) -> Tuple[UInt64, UInt64]:
    """A snarl read in either direction, as one key: files often record a parent reversed."""
    var rs = end ^ 1
    var re = start ^ 1
    if rs < start or (rs == start and re < end):
        return (rs, re)
    return (start, end)


@always_inline
def put_visit(mut out: List[UInt8], v: UInt64):
    """`>12` or `<12`, as vg writes a visit in a snarl's name."""
    out.append(UInt8(ord("<")) if v & 1 == 1 else UInt8(ord(">")))
    var node = Int(v >> 1)
    var digits = 1
    var x = node
    while x >= 10:
        x //= 10
        digits += 1
    var at = len(out)
    out.resize(at + digits, 0)
    for k in range(digits):
        out[at + digits - 1 - k] = UInt8(ord("0") + node % 10)
        node //= 10


def snarl_name(start: UInt64, end: UInt64) -> String:
    """`>12>34`: the snarl's id as vg writes it."""
    var out = List[UInt8]()
    put_visit(out, start)
    put_visit(out, end)
    return String(unsafe_from_utf8=out)


def _bad[W: Writable](what: W) -> Error:
    return Error(t"malformed snarls file: {what}")


@always_inline
def _varint(p: ImmPointer[UInt8, _], mut at: Int, end: Int) raises -> UInt64:
    var shift: UInt64 = 0
    var v: UInt64 = 0
    while True:
        if at >= end:
            raise _bad("truncated varint")
        var c = p[unsafe_offset=at]
        at += 1
        v |= UInt64(c & 0x7F) << shift
        if c < 0x80:
            return v
        shift += 7
        if shift > 63:
            raise _bad("varint too long")


@always_inline
def _take(mut at: Int, end: Int, length: Int) raises -> Int:
    """The start of the `length` bytes at `at`, which moves past them."""
    if length > end - at:
        raise _bad("truncated message")
    var start = at
    at += length
    return start


def _visit(p: ImmPointer[UInt8, _], start: Int, end: Int) raises -> UInt64:
    """A `vg.Visit`'s node id (field 1) and backward flag (field 3), as a visit."""
    var id: UInt64 = 0
    var back = False
    var at = start
    while at < end:
        var key = _varint(p, at, end)
        var wire = key & 7
        var field = key >> 3
        if wire == 0:
            var v = _varint(p, at, end)
            if field == 1:
                id = v
            elif field == 3:
                back = v != 0
        elif wire == 2:
            _ = _take(at, end, Int(_varint(p, at, end)))
        elif wire == 1:
            _ = _take(at, end, 8)
        elif wire == 5:
            _ = _take(at, end, 4)
        else:
            raise _bad(t"protobuf wire type {wire}")
    return id << 1 | UInt64(Int(back))


def _snarl(p: ImmPointer[UInt8, _], start: Int, end: Int, nested: Bool) raises -> Snarl:
    """A `vg.Snarl` message's start and end visits (fields 2 and 3) and, unless `nested`, its parent's (field 4,
    itself a Snarl holding only those)."""
    var s = NO_PARENT
    var e = NO_PARENT
    var has_start = False
    var has_end = False
    var parent_start = NO_PARENT
    var parent_end = NO_PARENT
    var at = start
    while at < end:
        var key = _varint(p, at, end)
        var wire = key & 7
        var field = key >> 3
        if wire == 0:
            _ = _varint(p, at, end)
        elif wire == 2:
            var length = Int(_varint(p, at, end))
            var b = _take(at, end, length)
            if field == 2:
                s = _visit(p, b, b + length)
                has_start = True
            elif field == 3:
                e = _visit(p, b, b + length)
                has_end = True
            elif field == 4:
                var parent = _snarl(p, b, b + length, True)
                parent_start = parent.start
                parent_end = parent.end
        elif wire == 1:
            _ = _take(at, end, 8)
        elif wire == 5:
            _ = _take(at, end, 4)
        else:
            raise _bad(t"protobuf wire type {wire}")
    if not has_start:
        raise _bad("snarl without start")
    if not has_end:
        raise _bad("snarl without end")
    return Snarl(s, e, parent_start, parent_end)


@fieldwise_init
struct _InflateTask(ParallelTask):
    """BGZF blocks `index * span ..` of the file inflated into place."""

    var lib: ImmPointer[OwnedDLHandle, ImmUntrackedOrigin]
    var data: ImmPointer[UInt8, ImmUntrackedOrigin]
    var blocks: ImmPointer[Int, ImmUntrackedOrigin]
    var count: Int
    var span: Int
    var out: MutPointer[UInt8, MutUntrackedOrigin]
    var bad: MutPointer[Int, MutUntrackedOrigin]

    def run(self, index: Int):
        var first = index * self.span
        if not inflate_blocks(self.lib[], self.data, self.blocks, self.out, first, min(first + self.span, self.count)):
            self.bad[] = 1


def _decompressed(var raw: List[UInt8], path: String, threads: Int) raises -> List[UInt8]:
    """The file's bytes, inflated if gzip (BGZF on every thread, its blocks' CRCs checked; other gzip as one
    stream, its members one after another)."""
    var kind = encoding(Span(raw)[: min(len(raw), 18)])
    if kind == BGZF:
        var blocks = List[Int]()
        var whole = True
        var got = (0, 0)
        try:
            got = scan_blocks(Span(raw), 0, blocks, path)
        except:
            whole = False
        if whole and got[0] == len(raw):
            var lib = load_libdeflate()
            var out = List[UInt8](unsafe_uninit_length=got[1])
            var count = len(blocks) // 4
            var span = 256
            var bad = 0
            var task = _InflateTask(
                imm_untracked(Pointer(to=lib)),
                imm_untracked(raw.unsafe_ptr()),
                imm_untracked(blocks.unsafe_ptr()),
                count,
                span,
                mut_untracked(out.unsafe_ptr()),
                mut_untracked(Pointer(to=bad)),
            )
            parallel_for(task, (count + span - 1) // span, threads)
            # The tasks reach these through untracked pointers: reads after them keep them alive until they are done.
            _ = len(raw) + len(blocks)
            _ = lib^
            if bad != 0:
                raise Error(t"{path}: corrupt BGZF block")
            return out^
        kind = GZIP
    if kind == GZIP:
        var zlib = load_zlib()
        var stream = GzipStream()
        var out = List[UInt8]()
        stream.inflate(zlib, raw.unsafe_ptr(), len(raw), out)
        var done = stream.finished()
        stream.free(zlib)
        if not done:
            raise Error(t"{path}: truncated gzip")
        return out^
    return raw^


@fieldwise_init
struct _DecodeTask(ParallelTask):
    """Messages `index * DECODE_SPAN ..` decoded into place; a malformed one marks `bad`."""

    var data: ImmPointer[UInt8, ImmUntrackedOrigin]
    var messages: ImmPointer[Int, ImmUntrackedOrigin]
    """Each message's start and end."""
    var count: Int
    var out: MutPointer[Snarl, MutUntrackedOrigin]
    var bad: MutPointer[Int, MutUntrackedOrigin]

    def run(self, index: Int):
        var first = index * DECODE_SPAN
        for i in range(first, min(first + DECODE_SPAN, self.count)):
            try:
                self.out[unsafe_offset=i] = _snarl(
                    self.data, self.messages[unsafe_offset=2 * i], self.messages[unsafe_offset=2 * i + 1], False
                )
            except:
                self.bad[] = 1
                return


def parse_snarls(data: Span[UInt8, _], threads: Int) raises -> List[Snarl]:
    """The snarls of a decompressed stream, in file order (messages decoded on `threads` threads)."""
    var p = data.unsafe_ptr()
    var n = len(data)
    var messages = List[Int]()
    var at = 0
    while at < n:
        var count = _varint(p, at, n)
        var tag_length = Int(_varint(p, at, n))
        var tag = _take(at, n, tag_length)
        var is_snarl = tag_length == SNARL_TAG.byte_length()
        if is_snarl:
            var name = SNARL_TAG.as_bytes()
            for k in range(tag_length):
                var c = p[unsafe_offset=tag + k]
                if c >= UInt8(ord("a")) and c <= UInt8(ord("z")):
                    c -= 32
                is_snarl = is_snarl and c == name[k]
        var m = UInt64(1)
        while m < count:
            var length = Int(_varint(p, at, n))
            var start = _take(at, n, length)
            if is_snarl:
                messages.append(start)
                messages.append(start + length)
            m += 1
    var count = len(messages) // 2
    var out = List[Snarl](unsafe_uninit_length=count)
    var bad = 0
    var task = _DecodeTask(
        untracked(p),
        imm_untracked(messages.unsafe_ptr()),
        count,
        mut_untracked(out.unsafe_ptr()),
        mut_untracked(Pointer(to=bad)),
    )
    parallel_for(task, (count + DECODE_SPAN - 1) // DECODE_SPAN, threads)
    _ = len(messages)
    if bad != 0:
        # Again in order, for the first malformed message's error.
        for i in range(count):
            _ = _snarl(p, messages[2 * i], messages[2 * i + 1], False)
    return out^


def read_snarls(path: String, threads: Int) raises -> List[Snarl]:
    """The snarls of a `vg snarls` file, in file order."""
    var f = RawFile(path)
    var raw = List[UInt8](unsafe_uninit_length=f.size)
    if f.read(raw.unsafe_ptr(), f.size, 0) != f.size:
        raise Error(t"{path}: cannot read")
    var data = _decompressed(raw^, path, threads)
    return parse_snarls(Span(data), threads)
