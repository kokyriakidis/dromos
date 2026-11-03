"""BAM as htslib writes it (SAMv1 section 4.2): the header, records encoded the way `samtools view -b` encodes the
same SAM lines (sam_parse1, bam_write1: the bin of [pos, pos + max(1, reference length)) at 14 bits and 5 levels,
sequence letters by seq_nt16_table, a CIGAR of more than 65,535 operations moved to a CG tag), and the .bai or .csi
index `samtools index` makes of a sorted BAM."""

from std.ffi import OwnedDLHandle
from std.memory import unsafe_memcpy

from dromos.linearize import D, I, M, S, CigarOp
from dromos.tabix import MIN_SHIFT, N_LVLS, BinIndex, BlockIndex, csi_levels, reg2bin, write_bgzf_file

comptime BAI_MAX = 1 << (MIN_SHIFT + 3 * N_LVLS)
"""The farthest end a .bai can index (2^29); a longer contig takes a .csi."""
comptime MAX_NAME = 254
"""The longest read name BAM holds."""
comptime MAX_OPS = 0xFFFF
"""The most CIGAR operations a record holds in place; more go to its CG tag."""
comptime FIXED = 36
"""A record's bytes before its name: block_size, then the 32 bytes of its fixed fields."""

comptime _NT16 = "=ACMGRSVTWYHKDBN"


def nt16(c: UInt8) -> UInt8:
    """htslib's seq_nt16_table: a letter's 4-bit code (either case; U as T; 0-3 as ACGT; = as 0; anything else 15)."""
    var u = c & 0xDF if c >= UInt8(ord("a")) and c <= UInt8(ord("z")) else c
    if u >= UInt8(ord("0")) and u <= UInt8(ord("3")):
        return UInt8(1) << (u - UInt8(ord("0")))
    if u == UInt8(ord("U")):
        return 8
    if u == UInt8(ord("=")):
        return 0
    var t = _NT16.as_bytes()
    for k in range(1, 15):
        if t[k] == u:
            return UInt8(k)
    return 15


def _nt16_table() -> Array[UInt8, 256]:
    var t = Array[UInt8, 256](fill=15)
    for c in range(256):
        t[c] = nt16(UInt8(c))
    return t^


@always_inline
def cigar_code(op: UInt8) -> Int:
    """A CIGAR operation's BAM code (MIDNSHP=X)."""
    if op == M:
        return 0
    if op == I:
        return 1
    if op == D:
        return 2
    if op == S:
        return 4
    var k = "MIDNSHP=X".as_bytes()
    for i in range(len(k)):
        if k[i] == op:
            return i
    return 0


@always_inline
def ref_length(cigar: List[CigarOp]) -> Int:
    """Reference bases a CIGAR covers (M, D, N, = and X)."""
    var n = 0
    for o in cigar:
        var c = cigar_code(o.op)
        if c == 0 or c == 2 or c == 3 or c == 7 or c == 8:
            n += o.length
    return n


@always_inline
def put32(p: MutPointer[UInt8, _], v: Int):
    p[unsafe_offset=0] = UInt8(v & 0xFF)
    p[unsafe_offset=1] = UInt8((v >> 8) & 0xFF)
    p[unsafe_offset=2] = UInt8((v >> 16) & 0xFF)
    p[unsafe_offset=3] = UInt8((v >> 24) & 0xFF)


@always_inline
def get32(p: ImmPointer[UInt8, _]) -> Int:
    return (
        Int(p[unsafe_offset=0])
        | Int(p[unsafe_offset=1]) << 8
        | Int(p[unsafe_offset=2]) << 16
        | Int(p[unsafe_offset=3]) << 24
    )


@always_inline
def get_signed32(p: ImmPointer[UInt8, _]) -> Int:
    var v = get32(p)
    return v - (1 << 32) if v >= 1 << 31 else v


struct RecordWriter(Movable):
    """Writes BAM records: their letters coded through a table made once."""

    var table: Array[UInt8, 256]

    def __init__(out self):
        self.table = _nt16_table()

    def put(
        self,
        mut out: List[UInt8],
        name: ImmPointer[UInt8, _],
        name_len: Int,
        flag: Int,
        tid: Int,
        pos: Int,
        mapq: Int,
        cigar: List[CigarOp],
        seq: Span[UInt8, _],
        qual: ImmPointer[UInt8, _],
        qual_at: Int,
        reverse_qual: Bool,
    ) raises:
        """Appends a mapped record without mate or tags. Its qualities are the len(seq) SAM characters at
        `qual[qual_at:]` (reversed with `reverse_qual`), or all 'I' when `qual_at` is negative."""
        if name_len > MAX_NAME:
            raise Error("a read name longer than BAM's 254 characters")
        if mapq > 255:
            raise Error(t"a mapping quality of {mapq}, above BAM's 255")
        if pos > 0x7FFF_FFFF:
            raise Error("a position beyond BAM's 2^31")
        var n_ops = len(cigar)
        var rlen = ref_length(cigar)
        var long_cigar = n_ops > MAX_OPS
        var l_seq = len(seq)
        var stored_ops = 2 if long_cigar else n_ops
        var size = FIXED + name_len + 1 + 4 * stored_ops + (l_seq + 1) // 2 + l_seq
        if long_cigar:
            size += 8 + 4 * n_ops
        var at = len(out)
        out.resize(unsafe_uninit_length=at + size)
        var p = out.unsafe_ptr().unsafe_offset(at)
        put32(p, size - 4)
        put32(p.unsafe_offset(4), tid)
        put32(p.unsafe_offset(8), pos)
        var bin = reg2bin(pos, pos + max(rlen, 1)) & 0xFFFF
        put32(p.unsafe_offset(12), bin << 16 | mapq << 8 | (name_len + 1))
        put32(p.unsafe_offset(16), flag << 16 | stored_ops)
        put32(p.unsafe_offset(20), l_seq)
        put32(p.unsafe_offset(24), -1)
        put32(p.unsafe_offset(28), -1)
        put32(p.unsafe_offset(32), 0)
        var q = p.unsafe_offset(FIXED)
        unsafe_memcpy(dest=q, src=name, count=name_len)
        q[unsafe_offset=name_len] = 0
        q = q.unsafe_offset(name_len + 1)
        if long_cigar:
            # bam_write1: <read length>S<reference length>N in place, the real CIGAR in CG:B,I after the qualities.
            put32(q, l_seq << 4 | 4)
            put32(q.unsafe_offset(4), rlen << 4 | 3)
            q = q.unsafe_offset(8)
        else:
            for o in cigar:
                put32(q, o.length << 4 | cigar_code(o.op))
                q = q.unsafe_offset(4)
        var s = seq.unsafe_ptr()
        for k in range(l_seq // 2):
            q[unsafe_offset=k] = (
                self.table[Int(s[unsafe_offset=2 * k])] << 4 | self.table[Int(s[unsafe_offset=2 * k + 1])]
            )
        if l_seq % 2 == 1:
            q[unsafe_offset=l_seq // 2] = self.table[Int(s[unsafe_offset=l_seq - 1])] << 4
        var qs = q.unsafe_offset((l_seq + 1) // 2)
        if qual_at < 0:
            for k in range(l_seq):
                qs[unsafe_offset=k] = UInt8(ord("I") - 33)
        elif reverse_qual:
            var last = qual_at + l_seq - 1
            for k in range(l_seq):
                qs[unsafe_offset=k] = qual[unsafe_offset=last - k] - 33
        else:
            var first = qual.unsafe_offset(qual_at)
            for k in range(l_seq):
                qs[unsafe_offset=k] = first[unsafe_offset=k] - 33
        if long_cigar:
            var t = qs.unsafe_offset(l_seq)
            t[unsafe_offset=0] = UInt8(ord("C"))
            t[unsafe_offset=1] = UInt8(ord("G"))
            t[unsafe_offset=2] = UInt8(ord("B"))
            t[unsafe_offset=3] = UInt8(ord("I"))
            put32(t.unsafe_offset(4), n_ops)
            t = t.unsafe_offset(8)
            for o in cigar:
                put32(t, o.length << 4 | cigar_code(o.op))
                t = t.unsafe_offset(4)


def record_end(p: ImmPointer[UInt8, _]) -> Int:
    """The end of the reference span of the record at `p` (its block_size first), as htslib's bam_endpos gives it:
    pos + reference length, at least pos + 1. A CIGAR moved to CG keeps its length in the N of the one in place."""
    var pos = get_signed32(p.unsafe_offset(8))
    var name_len = get32(p.unsafe_offset(12)) & 0xFF
    var n_ops = get32(p.unsafe_offset(16)) & 0xFFFF
    var c = p.unsafe_offset(FIXED + name_len)
    var rlen = 0
    for k in range(n_ops):
        var v = get32(c.unsafe_offset(4 * k))
        var code = v & 15
        if code == 0 or code == 2 or code == 3 or code == 7 or code == 8:
            rlen += v >> 4
    return pos + max(rlen, 1)


def bam_header(text: String, names: List[String], lengths: List[Int]) -> List[UInt8]:
    """A BAM's header: magic, the SAM header text, then each reference's name and length."""
    var out = List[UInt8]()
    out.extend(Span("BAM\x01".as_bytes()))
    var t = text.as_bytes()
    _append32(out, len(t))
    out.extend(t)
    _append32(out, len(names))
    for i in range(len(names)):
        var n = names[i].as_bytes()
        _append32(out, len(n) + 1)
        out.extend(n)
        out.append(0)
        _append32(out, lengths[i])
    return out^


@always_inline
def _append32(mut out: List[UInt8], v: Int):
    for k in range(4):
        out.append(UInt8((v >> (8 * k)) & 0xFF))


@always_inline
def _append64(mut out: List[UInt8], v: Int):
    for k in range(8):
        out.append(UInt8((v >> (8 * k)) & 0xFF))


def needs_csi(lengths: List[Int]) -> Bool:
    """Whether a contig is too long for a .bai (samtools index would refuse its records)."""
    for n in lengths:
        if n > BAI_MAX:
            return True
    return False


struct BamIndex(BlockIndex):
    """A BAM's index as `samtools index` makes it: a .bai (14 bits, 5 levels, every reference listed), or, with
    `csi`, a .csi as `samtools index -c` makes it (14 bits, as deep as the longest reference needs)."""

    var index: BinIndex
    var csi: Bool

    def __init__(out self, names: List[String], lengths: List[Int], csi: Bool):
        var max_len = 0
        for n in lengths:
            max_len = max(max_len, n)
        self.csi = csi
        self.index = BinIndex(MIN_SHIFT, csi_levels(max_len, MIN_SHIFT) if csi else N_LVLS)
        for n in names:
            self.index.add_contig(n)

    def suffix(self) -> String:
        return ".csi" if self.csi else ".bai"

    def set_first_offset(mut self, offset: Int):
        self.index.set_first_offset(offset)

    def push(mut self, contig: Int, beg: Int, end: Int, offset: Int) raises:
        self.index.push(contig, beg, end, offset)

    def finish(mut self, final_offset: Int) raises:
        self.index.finish(final_offset)

    def bytes(self) -> List[UInt8]:
        """The index, uncompressed: a .bai's bytes, or a .csi's before BGZF."""
        var raw = List[UInt8]()
        if self.csi:
            raw.extend(Span("CSI\x01".as_bytes()))
            _append32(raw, self.index.min_shift)
            _append32(raw, self.index.n_lvls)
            _append32(raw, 0)
        else:
            raw.extend(Span("BAI\x01".as_bytes()))
        _append32(raw, len(self.index.contigs))
        for ci in range(len(self.index.contigs)):
            self.index.put_bins(ci, self.csi, raw)
            if not self.csi:
                self.index.put_linear(ci, raw)
        _append64(raw, 0)  # records without coordinates
        return raw^

    def write(self, lib: OwnedDLHandle, path: String) raises:
        """A .bai as is (htslib writes it uncompressed), a .csi as BGZF."""
        var raw = self.bytes()
        if self.csi:
            write_bgzf_file(lib, raw, path)
            return
        with open(path, "w") as f:
            f.write_bytes(Span(raw))
