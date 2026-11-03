"""A tabix (.tbi) index, built the way htslib's `tabix` builds one for `-s 1 -b 2 -e 3` (hts_idx_push, hts_idx_finish,
compress_binning): a port of pgcat's tabix.rs, so the two write the same index for the same file. Bins are written
in numeric order (htslib writes them in hash order; readers do not depend on it)."""

from std.ffi import OwnedDLHandle

from dromos.bgzf import BLOCK, bgzf_block, eof_block

comptime MIN_SHIFT = 14
comptime N_LVLS = 5
comptime N_BINS = ((1 << (3 * N_LVLS + 3)) - 1) // 7
comptime META_BIN = N_BINS + 1
"""The pseudo-bin holding each contig's offsets and record counts."""
comptime MIN_MARKER_DIST = 0x10000
comptime NO_OFFSET = -1
"""A linear-index window no record has reached (u64::MAX on disk)."""
comptime NO_BIN = -1


def reg2bin(beg: Int, end: Int) -> Int:
    var e = end - 1
    var s = MIN_SHIFT
    var t = ((1 << (3 * N_LVLS)) - 1) // 7
    for l in range(N_LVLS, 0, -1):
        if beg >> s == e >> s:
            return t + (beg >> s)
        # htslib's loop decrements the level before this step.
        s += 3
        t -= 1 << (3 * (l - 1))
    return 0


def bin_first(level: Int) -> Int:
    return ((1 << (3 * level)) - 1) // 7


def bin_parent(bin: Int) -> Int:
    return (bin - 1) >> 3


def bin_level(bin: Int) -> Int:
    var l = 0
    var b = bin
    while b != 0:
        b = bin_parent(b)
        l += 1
    return l


@fieldwise_init
struct Chunk(TrivialRegisterPassable):
    var beg: Int
    var end: Int


struct _Bins(Movable):
    """A contig's bins: each bin's chunks, by bin number."""

    var ids: List[Int]
    """The bins, in the order they were first made (NO_BIN once one has joined its parent)."""
    var chunks: List[List[Chunk]]
    var at: Dict[Int, Int]
    """Each bin's place in `ids`."""
    var linear: List[Int]

    def __init__(out self):
        self.ids = List[Int]()
        self.chunks = List[List[Chunk]]()
        self.at = Dict[Int, Int]()
        self.linear = List[Int]()

    def is_empty(self) -> Bool:
        return len(self.ids) == 0

    def insert(mut self, bin: Int, c: Chunk):
        var k = self.at.get(bin)
        if k:
            self.chunks[k.value()].append(c)
            return
        self.at[bin] = len(self.ids)
        self.ids.append(bin)
        self.chunks.append([c])

    def sorted_bins(self) -> List[Int]:
        """The bins present, in numeric order (a BTreeMap's)."""
        var out = List[Int](capacity=len(self.ids))
        for b in self.ids:
            if b != NO_BIN:
                out.append(b)
        sort(Span(out))
        return out^


def _sort_by_beg(mut list: List[Chunk]):
    """Sorts chunks by their start, keeping the order of equal ones (Rust's sort_by_key)."""

    def less(a: Chunk, b: Chunk) -> Bool:
        return a.beg < b.beg

    sort[stable=True](Span(list), less)


def _compress_binning(mut bins: _Bins) raises:
    """htslib's compress_binning: a bin whose chunks span under 64 kb of compressed file joins its parent, deepest
    level first; then chunks that start in the same BGZF block as the previous one ends are merged."""
    for level in range(N_LVLS, 0, -1):
        var start = bin_first(level)
        var at_level = List[Int]()
        for b in bins.sorted_bins():
            if b >= start and b < N_BINS and bin_level(b) == level:
                at_level.append(b)
        for b in at_level:
            var k = bins.at.get(b).value()
            _sort_by_beg(bins.chunks[k])
            ref list = bins.chunks[k]
            var span = (list[len(list) - 1].end >> 16) - (list[0].beg >> 16)
            var parent = bin_parent(b)
            if span < MIN_MARKER_DIST and parent in bins.at:
                var moved = bins.chunks[k].copy()
                bins.chunks[k].clear()
                bins.ids[k] = NO_BIN
                _ = bins.at.pop(b)
                bins.chunks[bins.at.get(parent).value()].extend(Span(moved))
    for b in bins.sorted_bins():
        if b >= N_BINS:
            continue
        var k = bins.at.get(b).value()
        _sort_by_beg(bins.chunks[k])
        ref list = bins.chunks[k]
        var m = 0
        for l in range(1, len(list)):
            if list[m].end >> 16 >= list[l].beg >> 16:
                if list[m].end < list[l].end:
                    list[m].end = list[l].end
            else:
                m += 1
                list[m] = list[l]
        list.resize(m + 1, Chunk(0, 0))


struct TbiBuilder(Movable):
    """A .tbi index built record by record (hts_idx_t)."""

    var names: List[String]
    var tids: Dict[Int, Int]
    """Each contig number's place in `names`."""
    var contigs: List[_Bins]
    var save_tid: Int
    var last_tid: Int
    var save_bin: Int
    var last_bin: Int
    var save_off: Int
    var last_off: Int
    var off_beg: Int
    var last_coor: Int
    var n_mapped: Int
    var n_unmapped: Int

    def __init__(out self):
        self.names = List[String]()
        self.tids = Dict[Int, Int]()
        self.contigs = List[_Bins]()
        self.save_tid = -1
        self.last_tid = -1
        self.save_bin = NO_BIN
        self.last_bin = NO_BIN
        self.save_off = 0
        self.last_off = 0
        self.off_beg = 0
        self.last_coor = -1
        self.n_mapped = 0
        self.n_unmapped = 0

    def set_first_offset(mut self, offset: Int):
        """The virtual offset where the first record starts (after header lines)."""
        self.save_off = offset
        self.last_off = offset
        self.off_beg = offset

    def push(mut self, contig: Int, name: String, beg_in: Int, end_in: Int, offset: Int) raises:
        """A record on contig `contig` (named `name`) covering [beg, end) (0-based), whose line ends at virtual offset
        `offset`. Records must come contig by contig, sorted."""
        var tid: Int
        var known = self.tids.get(contig)
        if known:
            tid = known.value()
        else:
            tid = len(self.names)
            self.tids[contig] = tid
            self.names.append(name)
            self.contigs.append(_Bins())
        if end_in > 1 << (MIN_SHIFT + 3 * N_LVLS):
            raise Error(t"{name}:{end_in} is beyond what a .tbi index can hold")
        if self.last_tid != tid:
            if not self.contigs[tid].is_empty():
                raise Error(t"records of {name} are not together")
            self.last_tid = tid
            self.last_bin = NO_BIN
        elif self.last_coor > beg_in:
            raise Error(t"records of {name} are not sorted")
        if end_in < beg_in:
            raise Error(t"record on {name} ends before it begins")
        var beg = max(beg_in, 0)
        var end = max(end_in, 1)
        # Linear index: the first record overlapping each 16 kb window.
        ref linear = self.contigs[tid].linear
        var wb = beg >> MIN_SHIFT
        var we = (end - 1) >> MIN_SHIFT
        if len(linear) < we + 1:
            linear.resize(we + 1, NO_OFFSET)
        for w in range(wb, we + 1):
            if linear[w] == NO_OFFSET:
                linear[w] = self.last_off
        var bin = reg2bin(beg, end)
        if self.last_bin != bin:
            if self.save_bin != NO_BIN:
                self.contigs[self.save_tid].insert(self.save_bin, Chunk(self.save_off, self.last_off))
            if self.last_bin == NO_BIN and self.save_bin != NO_BIN:
                # A new contig: close the previous one's metadata.
                var off_end = self.last_off
                self.contigs[self.save_tid].insert(META_BIN, Chunk(self.off_beg, off_end))
                self.contigs[self.save_tid].insert(META_BIN, Chunk(self.n_mapped, self.n_unmapped))
                self.n_mapped = 0
                self.n_unmapped = 0
                self.off_beg = off_end
            self.save_off = self.last_off
            self.save_bin = bin
            self.last_bin = bin
            self.save_tid = tid
        self.n_mapped += 1
        self.last_off = offset
        self.last_coor = beg

    def finish(mut self, final_offset: Int) raises:
        """`final_offset`: the virtual offset after the last record."""
        if self.save_tid >= 0:
            self.contigs[self.save_tid].insert(self.save_bin, Chunk(self.save_off, final_offset))
            self.contigs[self.save_tid].insert(META_BIN, Chunk(self.off_beg, final_offset))
            self.contigs[self.save_tid].insert(META_BIN, Chunk(self.n_mapped, self.n_unmapped))
        for ci in range(len(self.contigs)):
            ref c = self.contigs[ci]
            # Windows no record starts in point to the next one that does.
            for i in range(len(c.linear) - 2, -1, -1):
                if c.linear[i] == NO_OFFSET:
                    c.linear[i] = c.linear[i + 1]
            _compress_binning(c)

    def bytes(self) -> List[UInt8]:
        """The index, uncompressed."""
        var raw = List[UInt8]()
        raw.extend(Span("TBI\x01".as_bytes()))
        _put32(raw, len(self.names))
        # tbx_conf_t: TBX_GENERIC, sequence column 1, begin 2, end 3, meta '#', skip 0.
        for v in [0, 1, 2, 3, ord("#"), 0]:
            _put32(raw, v)
        var names_len = 0
        for n in self.names:
            names_len += n.byte_length() + 1
        _put32(raw, names_len)
        for n in self.names:
            raw.extend(n.as_bytes())
            raw.append(0)
        for ci in range(len(self.contigs)):
            ref c = self.contigs[ci]
            var bins = c.sorted_bins()
            _put32(raw, len(bins))
            for b in bins:
                ref chunks = c.chunks[c.at.get(b).value()]
                _put32(raw, b)
                _put32(raw, len(chunks))
                for ch in chunks:
                    _put64(raw, ch.beg)
                    _put64(raw, ch.end)
            _put32(raw, len(c.linear))
            for o in c.linear:
                _put64(raw, o)
        _put64(raw, 0)  # records without coordinates
        return raw^

    def write(self, lib: OwnedDLHandle, path: String) raises:
        """The index as a BGZF file (level 6, as pgcat writes it)."""
        var raw = self.bytes()
        var out = List[UInt8]()
        var at = 0
        while at < len(raw):
            var n = min(BLOCK, len(raw) - at)
            bgzf_block(lib, 6, raw.unsafe_ptr().unsafe_offset(at), n, out)
            at += n
        eof_block(out)
        with open(path, "w") as f:
            f.write_bytes(Span(out))


@always_inline
def _put32(mut out: List[UInt8], v: Int):
    for k in range(4):
        out.append(UInt8((v >> (8 * k)) & 0xFF))


@always_inline
def _put64(mut out: List[UInt8], v: Int):
    for k in range(8):
        out.append(UInt8((v >> (8 * k)) & 0xFF))
