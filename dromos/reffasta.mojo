"""`dromos reference-fasta`: the reference FASTA of a graph's reference sample, every contig its paths spell, letters
from the graph's cache (its reference map and node table), and its .fai, written directly.

Each contig is laid out as the site catalog lays its reference sequence (docs/catalog-format.md): the sample's phase-0
paths of the contig (its fragments) in fragment-offset order, then path order, each step at its placement, a node's
letters in capitals (reverse-complemented when visited backwards: A, C, G, T to their complements, any other letter to
N), N where no fragment lies, as long as its farthest fragment end. So each contig's MD5 is the one `catalog-info`
prints for it, and contigs come in the catalog's order (natural: chr1 < chr2 < ... < chr10 < ... < chrM < chrX).

The file is laid out first (headers, LINE letters a line, each contig's offset), then cut into chunks of whole lines
that threads fill and write at their offsets."""

from std.math import ceildiv
from std.memory import unsafe_memcpy
from std.os import remove
from std.os.path import exists
from std.sys import stderr
from std.time import perf_counter_ns

from dromos.build import graph_stem, natural_less
from dromos.graph import Metadata, read_metadata, resolve_reference_sample
from dromos.nodes import NodeTable
from dromos.output import OutFile
from dromos.parallel import ParallelTask, parallel_for
from dromos.prepare import graph_cache
from dromos.refmap import ContigPath, RefMap
from dromos.stream import imm_untracked, mut_untracked, untracked

comptime LINE = 80
"""Letters a line, as chm13v2.0.fa has them."""
comptime CHUNK_LINES = 1 << 16
"""Lines a task fills and writes (5.2 MB of letters)."""


struct ReferenceFastaOptions(Movable):
    var graph: String
    var cache: Optional[String]
    var reference_sample: Optional[String]
    var out: Optional[String]
    """The FASTA [GRAPH.SAMPLE.fa beside the graph]; its index is OUT.fai."""
    var full_names: Bool
    """Contigs named `<sample>#0#<contig>` (the catalog's CHROM), not `<contig>`."""
    var threads: Int

    def __init__(out self):
        self.graph = String()
        self.cache = None
        self.reference_sample = None
        self.out = None
        self.full_names = False
        self.threads = 1


struct FastaContig(Movable):
    """A contig as written: its name, its length, and its non-empty fragments in the order they are laid."""

    var name: String
    var length: Int
    var fragments: List[ContigPath]

    def __init__(out self, var name: String, length: Int, var fragments: List[ContigPath]):
        self.name = name^
        self.length = length
        self.fragments = fragments^


def reference_contigs(
    map: RefMap, nodes: NodeTable, md: Metadata, full_names: Bool, threads: Int
) raises -> List[FastaContig]:
    """Every contig of the map's sample with a phase-0 path, in natural order, each with its phase-0 fragments (the
    catalog's reference path). The cache keeps the sample's paths of every phase in path id order without their
    phases or offsets: those come from the graph's metadata (`md`), whose paths of the sample are the cache's, one for
    one."""
    var sample_id = md.sample_id(map.sample)
    if not sample_id:
        raise Error(t"no paths of sample {map.sample} in the graph")
    var phase = List[Int]()
    var offset = List[Int]()
    var contig_id = List[Int]()
    for p in md.path_names():
        if p.sample == sample_id.value():
            phase.append(p.phase)
            offset.append(p.fragment)
            contig_id.append(p.contig)
    var mismatch = String(t"the graph's cache does not hold the graph's paths of {map.sample}: make it again")
    if len(phase) != len(map.path_contig):
        raise Error(mismatch)
    var checked = List[Bool](length=len(map.contigs), fill=False)
    for k in range(len(phase)):
        var c = map.path_contig[k]
        if not checked[c]:
            if md.contig_name(contig_id[k]) != map.contigs[c]:
                raise Error(mismatch)
            checked[c] = True

    var by_contig = map.paths_by_contig(threads)
    var names = List[String]()
    var at = Dict[String, Int]()
    for k in range(len(phase)):
        var c = map.path_contig[k]
        if phase[k] == 0 and map.contigs[c] not in at:
            at[map.contigs[c]] = c
            names.append(map.contigs[c])
    sort(Span(names), natural_less)

    def by_offset(a: Tuple[Int, Int], b: Tuple[Int, Int]) -> Bool:
        return a[0] < b[0] or (a[0] == b[0] and a[1] < b[1])

    var out = List[FastaContig](capacity=len(names))
    for name in names:
        var c = at[name]
        var chrom = String(t"{map.sample}#0#{name}")
        # (offset, path) of each fragment, in that order.
        var order = List[Tuple[Int, Int]]()
        for k in range(len(phase)):
            if map.path_contig[k] == c and phase[k] == 0:
                order.append((offset[k], k))
        sort(Span(order), by_offset)
        var length = 0
        var fragments = List[ContigPath]()
        for f in order:
            length = max(length, f[0])
            var j = 0
            while j < len(by_contig[c]) and by_contig[c][j].path != f[1]:
                j += 1
            if j == len(by_contig[c]):
                raise Error(mismatch)
            var cp = ContigPath(f[1], 0, List[Int]())
            swap(cp, by_contig[c][j])
            if len(cp.nodes) == 0:
                continue
            for node in cp.nodes:
                if node < 0 or nodes.length(node) == 0:
                    raise Error(t"{chrom}: a step on a node the graph's cache lacks")
            var last = len(cp.nodes) - 1
            length = max(length, map.placement(cp.first + last).start + nodes.length(cp.nodes[last]))
            fragments.append(cp^)
        out.append(FastaContig(chrom if full_names else name, length, fragments^))
    return out^


def fill_letters(map: RefMap, nodes: NodeTable, contig: FastaContig, a: Int, b: Int, dst: MutPointer[UInt8, _]):
    """Letters `a ..< b` of the contig at `dst`, where a fragment lies (a later fragment over an earlier one); the
    rest of `dst` is left as it is."""
    for ref cp in contig.fragments:
        var n = len(cp.nodes)
        var last = map.placement(cp.first + n - 1)
        if map.placement(cp.first).start >= b or last.start + nodes.length(cp.nodes[n - 1]) <= a:
            continue
        # The last step starting at or before `a` (step 0 if none does).
        var lo = 0
        var hi = n - 1
        while lo < hi:
            var mid = (lo + hi + 1) // 2
            if map.placement(cp.first + mid).start <= a:
                lo = mid
            else:
                hi = mid - 1
        for s in range(lo, n):
            var p = map.placement(cp.first + s)
            if p.start >= b:
                break
            var node = cp.nodes[s]
            var length = nodes.length(node)
            if p.start + length <= a:
                continue
            var i0 = max(a - p.start, 0)
            var i1 = min(b - p.start, length)
            var to = dst.unsafe_offset(p.start + i0 - a)
            var letters = nodes.start(node)
            if p.reverse:
                nodes.copy_reverse(to, letters + length - 1 - i0, i1 - i0)
            else:
                nodes.copy_forward(to, letters + i0, i1 - i0)


def fasta_layout(contigs: List[FastaContig]) -> List[Int]:
    """Each contig's offset in the file (of its first letter: its header line before it), then the file's size."""
    var out = List[Int](capacity=len(contigs) + 1)
    var at = 0
    for k in range(len(contigs)):
        ref c = contigs[k]
        at += c.name.byte_length() + 2
        out.append(at)
        at += c.length + ceildiv(c.length, LINE)
    out.append(at)
    return out^


def fai_text(contigs: List[FastaContig], layout: List[Int]) -> String:
    """The .fai, as samtools faidx writes it: name, length, offset, letters and bytes a line."""
    var out = String()
    for k in range(len(contigs)):
        out += String(t"{contigs[k].name}\t{contigs[k].length}\t{layout[k]}\t{LINE}\t{LINE + 1}\n")
    return out^


@fieldwise_init
struct _ChunkTask(ParallelTask):
    """Chunk `index`: its letters filled, cut into lines and written at their offset."""

    var map: ImmPointer[RefMap, ImmUntrackedOrigin]
    var nodes: ImmPointer[NodeTable, ImmUntrackedOrigin]
    var contigs: ImmPointer[List[FastaContig], ImmUntrackedOrigin]
    var layout: ImmPointer[List[Int], ImmUntrackedOrigin]
    var chunks: ImmPointer[List[Tuple[Int, Int]], ImmUntrackedOrigin]
    """Each chunk's contig and first letter."""
    var file: ImmPointer[OutFile, ImmUntrackedOrigin]
    var chunk_letters: Int
    var failed: MutPointer[Int, MutUntrackedOrigin]

    def run(self, index: Int):
        var c = self.chunks[][index][0]
        var a = self.chunks[][index][1]
        ref contig = self.contigs[][c]
        var b = min(a + self.chunk_letters, contig.length)
        var letters = List[UInt8](length=b - a, fill=UInt8(ord("N")))
        fill_letters(self.map[], self.nodes[], contig, a, b, letters.unsafe_ptr())
        var lines = ceildiv(b - a, LINE)
        var text = List[UInt8](unsafe_uninit_length=b - a + lines)
        var src = letters.unsafe_ptr()
        var dst = text.unsafe_ptr()
        for i in range(lines):
            var n = min(LINE, b - a - i * LINE)
            unsafe_memcpy(dest=dst.unsafe_offset(i * (LINE + 1)), src=src.unsafe_offset(i * LINE), count=n)
            dst[unsafe_offset=i * (LINE + 1) + n] = UInt8(ord("\n"))
        # `a` starts a line: the letters before it in the contig take a newline every LINE.
        if not self.file[].write_at(text.unsafe_ptr(), len(text), self.layout[][c] + a + a // LINE):
            self.failed[] = 1


def write_reference_fasta(
    map: RefMap,
    nodes: NodeTable,
    contigs: List[FastaContig],
    path: String,
    threads: Int,
    chunk_lines: Int = CHUNK_LINES,
) raises -> Int:
    """The contigs as a FASTA at `path` and its index at `path`.fai (written once the FASTA is: an index is never left
    beside a FASTA it does not describe); returns the FASTA's size."""
    var fai = path + ".fai"
    if exists(fai):
        remove(fai)
    var layout = fasta_layout(contigs)
    var file = OutFile(path)
    var headers = List[UInt8]()
    for k in range(len(contigs)):
        headers.clear()
        headers.append(UInt8(ord(">")))
        headers.extend(contigs[k].name.as_bytes())
        headers.append(UInt8(ord("\n")))
        file.put(Span(headers), layout[k] - len(headers))
    var chunk_letters = chunk_lines * LINE
    var chunks = List[Tuple[Int, Int]]()
    for k in range(len(contigs)):
        for a in range(0, contigs[k].length, chunk_letters):
            chunks.append((k, a))
    var failed = 0
    var task = _ChunkTask(
        untracked(Pointer(to=map)),
        untracked(Pointer(to=nodes)),
        untracked(Pointer(to=contigs)),
        imm_untracked(Pointer(to=layout)),
        imm_untracked(Pointer(to=chunks)),
        imm_untracked(Pointer(to=file)),
        chunk_letters,
        mut_untracked(Pointer(to=failed)),
    )
    parallel_for(task, len(chunks), threads)
    # The tasks reach these through untracked pointers: kept alive until here.
    _ = chunks^
    if failed != 0:
        raise Error(t"cannot write {path}")
    file.close()
    var index = OutFile(fai)
    index.put(Span(fai_text(contigs, layout).as_bytes()), 0)
    index.close()
    return layout[len(contigs)]


def run_reference_fasta(o: ReferenceFastaOptions) raises -> String:
    """`dromos reference-fasta`: the reference sample's contigs as a FASTA and its .fai; returns the FASTA's path."""
    var t0 = perf_counter_ns()
    var sample = resolve_reference_sample(o.graph, o.reference_sample)
    # A sample the graph lacks is refused before the cache is looked at (a cache would be made again for it).
    var md = read_metadata(o.graph)
    if not md.sample_id(sample):
        raise Error(t"no paths of sample {sample} in the graph")
    var cache = graph_cache(o.graph, o.cache, sample, o.threads)
    var map = RefMap(cache, sample)
    var nodes = NodeTable(cache)
    var contigs = reference_contigs(map, nodes, md, o.full_names, o.threads)
    if len(contigs) == 0:
        raise Error(t"the graph has no {sample}#0#<contig> reference paths")
    var path = o.out.or_else(graph_stem(o.graph) + "." + sample + ".fa")
    var size = write_reference_fasta(map, nodes, contigs, path, o.threads)
    var letters = 0
    for k in range(len(contigs)):
        letters += contigs[k].length
    print(
        (
            t"[dromos] {path}: {len(contigs)} contigs of {sample}, {letters} letters, {size} bytes, and its .fai"
            t" ({Float64(perf_counter_ns() - t0) * 1e-9} s)"
        ),
        file=stderr,
    )
    return path
