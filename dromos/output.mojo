"""Output: whole writes to a file descriptor, and fast integer formatting."""

from std.bit import count_leading_zeros
from std.ffi import external_call
from std.sys import exit, stderr
from std.sys.info import CompilationTarget


def write_all(fd: Int32, data: ImmPointer[UInt8, _], length: Int):
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


comptime _POW10 = SIMD[DType.int64, 32](
    1,
    10,
    100,
    1000,
    10000,
    100000,
    1000000,
    10000000,
    100000000,
    1000000000,
    10000000000,
    100000000000,
    1000000000000,
    10000000000000,
    100000000000000,
    1000000000000000,
    10000000000000000,
    100000000000000000,
    1000000000000000000,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
)


@always_inline
def decimal_digits(value: Int) -> Int:
    """How many decimal digits non-negative `value` has: from its bit length (times log10 2, as 1233 / 4096), less
    one when it falls short of that power of ten (Hacker's Delight)."""
    var bits = 64 - Int(count_leading_zeros(UInt64(value | 1)))
    var t = (bits * 1233) >> 12
    return max(1, t + Int(Int64(value) >= _POW10[t]))


@inline(.always)
def put_padded(dst: MutPointer[UInt8, _], value: Int, width: Int) -> Int:
    """Writes non-negative `value` right-aligned in `width` columns (wider
    when it needs more digits); returns the bytes written."""
    var digits = decimal_digits(value)
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
