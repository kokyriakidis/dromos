"""Buffered byte output to a file descriptor with fast integer formatting."""

from std.ffi import external_call
from std.sys import exit, stderr
from std.sys.info import CompilationTarget

comptime _FLUSH_AT = 1 << 20


struct OutputBuffer(Movable):
    """Accumulates bytes and writes them to `fd` in large blocks; with
    `fd = -1` it only accumulates (read them from `buf`)."""

    var buf: List[UInt8]
    var fd: Int32

    def __init__(out self, fd: Int32 = 1):
        self.buf = List[UInt8](capacity=_FLUSH_AT + 4096)
        self.fd = fd

    def write(mut self, text: StringSlice):
        self.write_bytes(text.as_bytes())

    def write_bytes(mut self, bytes: Span[UInt8, _]):
        if self.fd >= 0 and len(bytes) >= _FLUSH_AT:
            # Large blocks go out directly, without a copy into the buffer.
            self.flush()
            _write_all(self.fd, bytes.unsafe_ptr(), len(bytes))
            return
        self.buf.extend(bytes)
        if len(self.buf) >= _FLUSH_AT and self.fd >= 0:
            self.flush()

    @inline(.always)
    def write_byte(mut self, b: Int):
        self.buf.append(UInt8(b))

    def write_int(mut self, value: Int, width: Int = 0, left: Bool = False):
        """Writes `value` in decimal, padded with spaces to `width`."""
        if value >= 0 and not left:
            # Two digits per step, straight into the buffer.
            var at = len(self.buf)
            var need = at + max(width, 20)
            if self.buf.capacity() < need:
                self.buf.reserve(max(need, 2 * self.buf.capacity()))
            self.buf.resize(unsafe_uninit_length=at + put_padded(self.buf.unsafe_ptr().unsafe_offset(at), value, width))
            if len(self.buf) >= _FLUSH_AT and self.fd >= 0:
                self.flush()
            return
        var digits = InlineDigits(value)
        var pad = max(0, width - digits.count)
        if not left:
            for _ in range(pad):
                self.buf.append(UInt8(ord(" ")))
        for i in range(digits.count):
            self.buf.append(digits.bytes[digits.count - 1 - i])
        if left:
            for _ in range(pad):
                self.buf.append(UInt8(ord(" ")))
        if len(self.buf) >= _FLUSH_AT and self.fd >= 0:
            self.flush()

    def take(mut self) -> List[UInt8]:
        """The accumulated bytes (for `fd = -1`), leaving the buffer empty."""
        var bytes = self.buf^
        self.buf = List[UInt8]()
        return bytes^

    def flush(mut self):
        if self.fd < 0:
            return
        _write_all(self.fd, self.buf.unsafe_ptr(), len(self.buf))
        self.buf.clear()


def _write_all(fd: Int32, data: ImmPointer[UInt8, _], length: Int):
    """Writes every byte, retrying interrupted writes; any other failure (a
    full disk, a closed pipe) ends the program with an error, never silently
    truncated output."""
    var done = 0
    while done < length:
        var n = external_call["write", Int](fd, data.unsafe_offset(done), length - done)
        if n > 0:
            done += n
            continue
        if n < 0 and _errno() == _EINTR:
            continue
        print("dromos: write failed (errno ", _errno(), ")", sep="", file=stderr)
        exit(1)


comptime _EINTR: Int32 = 4


def _errno() -> Int32:
    comptime if CompilationTarget.is_macos():
        return external_call["__error", MutPointer[Int32, MutAnyOrigin]]()[]
    else:
        return external_call["__errno_location", MutPointer[Int32, MutAnyOrigin]]()[]


@fieldwise_init
struct InlineDigits(TrivialRegisterPassable):
    """Decimal digits of an integer, least significant first."""

    var bytes: SIMD[DType.uint8, 32]
    var count: Int

    def __init__(out self, value: Int):
        self.bytes = SIMD[DType.uint8, 32](0)
        self.count = 0
        var v = value
        var negative = v < 0
        if negative:
            v = -v
        while True:
            self.bytes[self.count] = UInt8(ord("0") + v % 10)
            self.count += 1
            v //= 10
            if v == 0:
                break
        if negative:
            self.bytes[self.count] = UInt8(ord("-"))
            self.count += 1


comptime _PAIRS = "00010203040506070809101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899"
"""The decimal digit pairs 00 to 99, for writing two digits per step."""


@inline(.always)
def put_padded(dst: MutPointer[UInt8, _], value: Int, width: Int) -> Int:
    """Writes non-negative `value` right-aligned in `width` columns (wider
    when it needs more digits); returns the bytes written."""
    var digits = 1
    var bound = 10
    while value >= bound and digits < 19:
        digits += 1
        bound *= 10
    var n = max(width, digits)
    for i in range(n - digits):
        dst[unsafe_offset=i] = UInt8(ord(" "))
    var table = _PAIRS.as_bytes().unsafe_ptr()
    var v = value
    var at = n
    while v >= 100:
        var q = v // 100
        var r = 2 * (v - 100 * q)
        at -= 2
        dst[unsafe_offset=at] = table[unsafe_offset=r]
        dst[unsafe_offset=at + 1] = table[unsafe_offset=r + 1]
        v = q
    if v >= 10:
        at -= 2
        dst[unsafe_offset=at] = table[unsafe_offset=2 * v]
        dst[unsafe_offset=at + 1] = table[unsafe_offset=2 * v + 1]
    else:
        at -= 1
        dst[unsafe_offset=at] = UInt8(ord("0") + v)
    return n
