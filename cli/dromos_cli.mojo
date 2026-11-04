"""Dromos: pangenome graph alignments, fast. `dromos COMMAND -h` for each command's usage."""

from std.ffi import external_call
from std.os.path import dirname
from std.sys import argv, exit, stderr
from std.time import perf_counter_ns

from dromos.build import BuildOptions, run_build
from dromos.bam import BamIndex, bam_header, needs_csi
from dromos.bam_sort import BamStaging, write_bam
from dromos.contig import Reference
from dromos.convert import GafSink
from dromos.graph import fingerprint, reference_samples, resolve_reference_sample
from dromos.index import IndexSink
from dromos.linearize import Window
from dromos.nodes import NodeTable
from dromos.output import OutFile
from dromos.parallel import default_threads
from dromos.prepare import graph_cache
from dromos.sites import Sites, info, load_all_windows, view
from dromos.snarlfind import find_all_snarls, write_records_tsv, write_snarls_tsv
from dromos.snarls import read_snarls
from dromos.refmap import RefMap
from dromos.stream import AUTO, GAF, GAM, stream_input
from dromos.surject import SurjectSink, Surjector

comptime OVERVIEW = """dromos: pangenome graph alignments, fast

  dromos index       a GAM or GAF placed on the reference, sorted and tabix-indexed, as pgphase reads it
  dromos gam-to-gaf  a GAM as GAF, the lines vg convert -G writes
  dromos surject     a GAM or GAF on the reference contigs, the records pgcatalog surject writes: SAM, or a
                     sorted, indexed BAM

  dromos catalog            a GBZ graph's site catalog: every snarl on the reference, every contig, one file
  dromos view               a catalog as VCF (pgcatalog view's text), whole or a region of it
  dromos catalog-info       a catalog's header: the graph it belongs to, its contigs, its snarls' source
  dromos snarls             a GBZ graph's snarls, the decomposition vg snarls makes, as a table
  dromos prepare            a GBZ graph's cache, made ahead of the commands that read it
  dromos reference-samples  the reference samples a GBZ graph declares

index, gam-to-gaf and surject take the GBZ graph (--graph) and read what they need of it from its cache, one file
made on first use and kept beside the graph (GRAPH.gbz.dromos; --cache FILE puts it elsewhere).

dromos COMMAND -h shows a command's usage. DROMOS_TIMING=1 reports where the time goes, on stderr.
"""

comptime INDEX_USAGE = """dromos index [options] --graph GRAPH.gbz IN OUT

IN's alignments placed on the reference, sorted, bgzipped and tabix-indexed as pgphase reads them: OUT and OUT.tbi,
the files pgphase's previous indexer (pgcat index-gaf --slim) wrote from the same alignments as GAF.

  --graph FILE              the pangenome graph (GBZ)
  --cache FILE              the graph's cache: made there on first use (and again when the graph changes), read
                            after [GRAPH.gbz.dromos, beside the graph]
  IN                        a GAM, or a GAF (plain, gzip or BGZF), told apart by their first bytes; - reads stdin
  --reference-sample NAME   the sample whose paths are the reference: CHM13 for CHM13#0#chr20 [the graph's declared
                            reference sample, if it declares exactly one; HPRC graphs declare CHM13 and GRCh38 and
                            need this]
  -t THREADS                threads [the machine's cores]
  -l LEVEL                  BGZF compression level, 0-12 [1]: 1 deflates about 3 times as fast as 6 into a file 40%
                            larger; 6 for a file to keep
  --gaf                     IN is a GAF (not told from its bytes)
  --gam                     IN is a GAM (not told from its bytes)
"""

comptime GAM_TO_GAF_USAGE = """dromos gam-to-gaf [options] --graph GRAPH.gbz GAM

GAM as GAF on stdout: the lines vg convert -G writes for it, in the GAM's order.

  --graph FILE              the pangenome graph (GBZ)
  --cache FILE              the graph's cache: made there on first use (and again when the graph changes), read
                            after [GRAPH.gbz.dromos, beside the graph]
  GAM                       the alignments; - reads stdin
  -t THREADS                threads [the machine's cores]
  --slim                    each line's 12 GAF columns alone (no cs, bq, AS or dv; the node letters are not read)
"""

comptime SURJECT_USAGE = """dromos surject [options] --graph GRAPH.gbz IN -o OUT.bam
dromos surject [options] --graph GRAPH.gbz --contig NAME IN > OUT.sam

IN's alignments as alignments to the reference contigs: the records pgcatalog surject writes for the same graph,
reads and options, one run for each contig, in one pass over IN. A read with bases on several contigs has a record on
each; reads with no base on any are left out.
  -o OUT.bam                a BAM sorted as samtools sort sorts it, and its index (OUT.bam.bai, or .csi when a
                            contig is longer than a .bai holds); without it, SAM on stdout in IN's order (a read's
                            records in contig order)
  --csi                     index the BAM as .csi (as samtools index -c) even when a .bai would hold it

  --graph FILE              the pangenome graph (GBZ)
  --cache FILE              the graph's cache: made there on first use (and again when the graph changes), read
                            after [GRAPH.gbz.dromos, beside the graph]
  IN                        a GAM, or a GAF with cs (and bq for qualities; plain, gzip or BGZF), told apart by their
                            first bytes; - reads stdin. dromos index's output keeps no cs and cannot be used
  --contig NAME             a reference contig (chr20); repeat it for several [every contig of the reference sample]
  --reference-sample NAME   the sample whose paths are the reference [the graph's declared reference sample, if it
                            declares exactly one]
  --sam-names full|short    the contigs' names in the SAM or BAM: CHM13#0#chr20, or chr20 [full]
  --sam-contig NAME         the one contig's name in the SAM or BAM (with a single --contig)
  --graph-placement         place read bases only where the graph puts them: a stretch of more than 50 bases off the
                            reference is not aligned to it, and such read ends are clipped
  --catalog FILE            the graph's site catalog (dromos catalog), holding every contig surjected to: a read
                            spanning one of a contig's snarl windows is aligned there as its allele, at minimap2's
                            asm5 costs
  -t THREADS                threads [the machine's cores]
  -l LEVEL                  the BAM's deflate level, 0 to 12 [4]
  -m SIZE                   memory for records waiting to be sorted, compressed (K, M or G); past half of it they
                            are spilled to --tmp-dir [8G]
  --tmp-dir DIR             where records are spilled [OUT.bam's folder]
  --gaf                     IN is a GAF (not told from its bytes)
  --gam                     IN is a GAM (not told from its bytes)

Each read is rebuilt from its walk and cs; a base on a node a contig visits once takes that node's position.
Between placed stretches, a substitution bubble is aligned base for base and an indel bubble globally; read ends off
the reference are aligned or soft-clipped, and ends diverging from the read's own graph alignment clipped.
"""

comptime PREPARE_USAGE = """dromos prepare [options] --graph GRAPH.gbz

The graph's cache, made now rather than by the first command to need it: for a graph in a folder those commands
cannot write to, made by someone who can. Node lengths and letters, and the placements of every reference sample the
graph declares, in one file read in place (memory-mapped). Nothing is done when the cache is current.

  --graph FILE              the pangenome graph (GBZ)
  --cache FILE              the graph's cache: made there on first use (and again when the graph changes), read
                            after [GRAPH.gbz.dromos, beside the graph]
  --reference-sample NAME   a sample to place beside those the graph declares (a graph declaring none needs one for
                            index and surject)
  -t THREADS                threads [the machine's cores]
"""


comptime CATALOG_USAGE = """dromos catalog --graph GRAPH.gbz [options]

The graph's site catalog, one file: every snarl the reference path passes through, at every nesting level, with all
its alleles, for every reference contig (or the regions asked for), as binary columns in position blocks
(docs/catalog-format.md). dromos view prints it as the VCF pgcatalog view prints; its header names the graph it was
built from, and dromos refuses it with another graph.

  --graph FILE              the pangenome graph (GBZ)
  --out PATH                the catalog's file, or a folder to put it in [GRAPH.SAMPLE.sites beside the graph;
                            GRAPH.SAMPLE.chr20_1-1000.sites for regions]
  --reference-sample NAME   the sample whose paths are the reference: CHM13 for CHM13#0#chr20 [the graph's declared
                            reference sample, if it declares exactly one; HPRC graphs declare CHM13 and GRCh38 and
                            need this]
  --region REGION           a contig (chr20) or part of one (chr20:10000000-12000000, 1-based, inclusive),
                            repeatable [every reference contig]; a region holds the records a query of the whole
                            contig's catalog returns
  -t, --threads N           threads [the machine's cores]
  --level N                 zstd level, 1-22 [9]; higher is smaller and slower

Snarls are computed from the graph (vg snarls' decomposition of each contig's component: no vg, no distance index),
unless one of these gives them:
  --dist FILE               a distance index of the graph (giraffe's .dist, read in place)
  --snarls-cache DIR        vg snarls of each contig (<contig>.snarls), kept for later runs on the same graph; vg chunk
                            and vg snarls make those missing
  --vg PATH                 the vg binary, for vg chunk and vg snarls, run for contigs not in the cache [vg]
  --snarls-jobs N           vg snarls runs at a time [6, at most --threads]; fewer when memory is short (each needs
                            about 20 times its component's size)

For a tabix-indexed VCF of one contig:
  dromos view GRAPH.CHM13.sites chr20 | bgzip > chr20.sites.vcf.gz && tabix -p vcf chr20.sites.vcf.gz
"""

comptime SNARLS_USAGE = """dromos snarls --graph GRAPH.gbz -o OUT.tsv [options]
dromos snarls --vg-file SNARLS -o OUT.tsv

The snarls vg snarls writes for a GBZ graph, computed (vg's IntegratedSnarlFinder and SnarlManager rules, every
connected component of the graph), as a table, a line per snarl: component, start, end, parent start, parent end
(visits as id << 1 | reverse; - for none), chain, rank in the chain, the chain's length (its trivial snarls too), the
chain's kind (0 a cactus cycle, 1 a bridge path, 2 a pinched cycle, 3 the root cycle) and whether its orientation
vote was tied (1). Each component's snarls come in vg snarls' order. With --vg-file, a vg snarls file's snarls as
the first five columns (component 0), to compare with (tests/integration/compare.py).

  --graph FILE              the pangenome graph (GBZ)
  -o FILE                   the table
  -t THREADS                threads [the machine's cores]
  --stats                   a line per component on stderr: sizes, each step's time, the ties met
"""

comptime VIEW_USAGE = """dromos view [-t THREADS] CATALOG [CONTIG[:START-END]]

A site catalog as VCF on stdout. With a contig (chr20 or CHM13#0#chr20), and a region of it (1-based, inclusive),
what pgcatalog view prints of that contig's catalog: its header, then the records overlapping the region, as tabix
would return them from the VCF. Without one, every contig: one VCF, a ##contig line for each.

  -t THREADS  threads, decoding blocks [the machine's cores]
"""

comptime CATALOG_INFO_USAGE = """dromos catalog-info CATALOG

A site catalog's header on stdout, a field a line (key, tab, value): the format, the graph it was built from (file,
fingerprint, node ids, GBZ and GBWT tags), the reference sample, where its snarls came from and how it was built;
then a line per contig: its name, CHROM, length, the MD5 of its reference sequence (SAM's M5), records, blocks and the
regions built.
"""

comptime REFERENCE_SAMPLES_USAGE = """dromos reference-samples GRAPH.gbz

The reference samples a GBZ graph declares, one per line, read from the start of the file.
"""


def _fail(message: String, usage_of: StaticString):
    print(t"dromos: {message}", file=stderr)
    print(t"see dromos {usage_of} -h", file=stderr)
    exit(2)


def _number(text: String, option: String, low: Int, high: Int, command: StaticString) -> Int:
    """`text` as a number from `low` to `high`, for `option`."""
    var n = 0
    try:
        n = Int(text)
    except:
        _fail(String(t"{option} needs a number, not {text}"), command)
    if n < low or n > high:
        _fail(String(t"{option} must be from {low} to {high}, not {n}"), command)
    return n


def _input(path: String) -> String:
    return String("/dev/stdin") if path == "-" else path


struct _Args(Movable):
    """A command's arguments: its options (each taken once named), then what is left, in order."""

    var items: List[String]
    var command: StaticString

    def __init__(out self, command: StaticString):
        self.items = List[String]()
        var args = argv()
        for i in range(2, len(args)):
            self.items.append(String(args[i]))
        self.command = command

    def help(self) -> Bool:
        for a in self.items:
            if a == "-h" or a == "--help":
                return True
        return False

    def flag(mut self, name: StaticString) -> Bool:
        """Whether `name` is given (taken out)."""
        for i in range(len(self.items)):
            if self.items[i] == name:
                _ = self.items.pop(i)
                return True
        return False

    def value(mut self, name: StaticString) -> Optional[String]:
        """The value after `name`, if given (both taken out)."""
        for i in range(len(self.items)):
            if self.items[i] == name:
                if i + 1 >= len(self.items):
                    _fail(String(t"{name} needs a value"), self.command)
                var v = self.items.pop(i + 1)
                _ = self.items.pop(i)
                return v^
        return None

    def values(mut self, name: StaticString) -> List[String]:
        """Every value given after `name` (all taken out), in order."""
        var out = List[String]()
        while True:
            var v = self.value(name)
            if not v:
                return out^
            out.append(v.value())

    def threads(mut self) -> Int:
        var t = self.value("-t")
        return _number(t.value(), "-t", 1, 1024, self.command) if t else default_threads()

    def positional_between(mut self, low: Int, high: Int) -> List[String]:
        """What is left: from `low` to `high` arguments, none of them an option."""
        for a in self.items:
            if a.startswith("-") and a != "-":
                _fail(String(t"unknown option {a}"), self.command)
        if len(self.items) < low or len(self.items) > high:
            _fail(String(t"{low} to {high} arguments wanted, {len(self.items)} given"), self.command)
        return self.items.copy()

    def positional(mut self, count: Int) -> List[String]:
        """What is left: exactly `count` arguments, none of them an option."""
        for a in self.items:
            if a.startswith("-") and a != "-":
                _fail(String(t"unknown option {a}"), self.command)
        if len(self.items) != count:
            _fail(String(t"{count} arguments wanted, {len(self.items)} given"), self.command)
        return self.items.copy()


def _graph(mut args: _Args, usage: StaticString) -> String:
    """The --graph value; without it, the usage on stderr."""
    var graph = args.value("--graph")
    if not graph:
        print(usage, end="", file=stderr)
        exit(1)
    return graph.value()


def index_main() raises:
    var args = _Args("index")
    if args.help():
        print(INDEX_USAGE, end="")
        return
    var threads = args.threads()
    var level_arg = args.value("-l")
    var level = _number(level_arg.value(), "-l", 0, 12, "index") if level_arg else 1
    var sample_arg = args.value("--reference-sample")
    var cache_arg = args.value("--cache")
    var gaf = args.flag("--gaf")
    var gam = args.flag("--gam")
    if gaf and gam:
        _fail("--gaf or --gam, not both", "index")
    var graph = _graph(args, INDEX_USAGE)
    var pos = args.positional(2)
    var sample = resolve_reference_sample(graph, sample_arg)
    var cache = graph_cache(graph, cache_arg, sample, threads)
    var t0 = perf_counter_ns()
    var sink = IndexSink(RefMap(cache, sample), threads)
    print(
        (
            t"[dromos] the graph's cache: {len(sink.map.contigs)} contigs of {sample}"
            t" ({Float64(perf_counter_ns() - t0) * 1e-9} s)"
        ),
        file=stderr,
    )
    stream_input(_input(pos[0]), threads, sink, GAF if gaf else (GAM if gam else AUTO))
    sink.finish(pos[1], level, threads)


def gam_to_gaf_main() raises:
    var args = _Args("gam-to-gaf")
    if args.help():
        print(GAM_TO_GAF_USAGE, end="")
        return
    var threads = args.threads()
    var slim = args.flag("--slim")
    var cache_arg = args.value("--cache")
    var graph = _graph(args, GAM_TO_GAF_USAGE)
    var pos = args.positional(1)
    var cache = graph_cache(graph, cache_arg, None, threads)
    if slim:
        var sink = GafSink[False](NodeTable(cache, letters=False), threads)
        try:
            stream_input(_input(pos[0]), threads, sink, AUTO)
        finally:
            sink.finish()
    else:
        var sink = GafSink[True](NodeTable(cache), threads)
        try:
            stream_input(_input(pos[0]), threads, sink, AUTO)
        finally:
            sink.finish()


def _size(text: String, option: String, command: StaticString) -> Int:
    """`text` as bytes: a number with an optional K, M or G (powers of 1024)."""
    var n = text.byte_length()
    var scale = 1
    var digits = text
    if n > 0:
        var last = String(text[byte = n - 1 :]).upper()
        if last == "K" or last == "M" or last == "G":
            scale = 1 << 10 if last == "K" else (1 << 20 if last == "M" else 1 << 30)
            digits = String(text[byte = : n - 1])
    return _number(digits, option, 1, Int.MAX // scale, command) * scale


def _writable(path: String) raises:
    """Makes and removes `path`: an error now if it cannot be written."""
    _ = OutFile(path)
    var c = List[UInt8](path.as_bytes())
    c.append(0)
    _ = external_call["unlink", Int32](c.unsafe_ptr())


def _same_graph(sites: Sites, graph: String) raises:
    """Refuses a catalog built from another graph than `graph` (their fingerprints differ)."""
    var want = sites.prop("graph.fingerprint")
    var have = fingerprint(graph)
    if not want or want.value() != have:
        var built = sites.prop("graph.file").or_else(String("?"))
        raise Error(
            t"{sites.path} was built from another graph ({built}, {want.or_else(String('no fingerprint'))}); {graph}"
            t" is {have}"
        )


def surject_main() raises:
    var args = _Args("surject")
    if args.help():
        print(SURJECT_USAGE, end="")
        return
    var threads = args.threads()
    var contig_args = args.values("--contig")
    var sample_arg = args.value("--reference-sample")
    var sam_arg = args.value("--sam-contig")
    var names_arg = args.value("--sam-names")
    var catalog_arg = args.value("--catalog")
    var cache_arg = args.value("--cache")
    var out_arg = args.value("-o")
    var level_arg = args.value("-l")
    var memory_arg = args.value("-m")
    var tmp_arg = args.value("--tmp-dir")
    var placement = args.flag("--graph-placement")
    var csi = args.flag("--csi")
    var gaf = args.flag("--gaf")
    var gam = args.flag("--gam")
    if gaf and gam:
        _fail("--gaf or --gam, not both", "surject")
    var short_names = False
    if names_arg:
        if names_arg.value() == "short":
            short_names = True
        elif names_arg.value() != "full":
            _fail(String(t"--sam-names is full or short, not {names_arg.value()}"), "surject")
    var level = _number(level_arg.value(), "-l", 0, 12, "surject") if level_arg else 4
    var memory = _size(memory_arg.value(), "-m", "surject") if memory_arg else 8 << 30
    if not out_arg and (level_arg or memory_arg or tmp_arg or csi):
        _fail("-l, -m, --tmp-dir and --csi are for a BAM (-o)", "surject")
    var graph = _graph(args, SURJECT_USAGE)
    var pos = args.positional(1)
    var sample = resolve_reference_sample(graph, sample_arg)
    var cache = graph_cache(graph, cache_arg, sample, threads)
    var t0 = perf_counter_ns()
    var map = RefMap(cache, sample)
    var nodes = NodeTable(cache)
    # The contigs, in reference order: those named, or all.
    var names = List[String]()
    for c in map.contigs:
        var wanted = len(contig_args) == 0
        for a in contig_args:
            wanted = wanted or a == c
        if wanted:
            names.append(c)
    for a in contig_args:
        var known = False
        for c in map.contigs:
            known = known or a == c
        if not known:
            raise Error(t"no reference path {sample}#0#{a} in the graph")
    if len(names) == 0:
        raise Error(t"no reference paths of {sample} in the graph")
    if sam_arg and len(names) != 1:
        _fail("--sam-contig names one contig: give one --contig (or --sam-names for several)", "surject")
    if sam_arg and names_arg:
        _fail("--sam-contig or --sam-names, not both", "surject")
    var reference = Reference(map, nodes, names, threads)
    var sam_names = List[String]()
    var lengths = List[Int]()
    for k in range(len(names)):
        if sam_arg:
            sam_names.append(sam_arg.value())
        elif short_names:
            sam_names.append(names[k])
        else:
            sam_names.append(reference.contigs[k].chrom)
        lengths.append(len(reference.contigs[k].seq))
    var windows = List[List[Window]](length=len(names), fill=List[Window]())
    if catalog_arg:
        var sites = Sites(catalog_arg.value())
        _same_graph(sites, graph)
        windows = load_all_windows(sites, names, reference, map, threads)
        var count = 0
        for w in windows:
            count += len(w)
        print(t"[dromos] surject: {count} snarl windows", file=stderr)
    if len(names) == 1:
        print(
            (
                t"[dromos] {reference.contigs[0].chrom}: {len(reference.contigs[0].seq)} bases"
                t" ({Float64(perf_counter_ns() - t0) * 1e-9} s)"
            ),
            file=stderr,
        )
    else:
        var bases = 0
        for n in lengths:
            bases += n
        print(
            t"[dromos] {len(names)} contigs of {sample}: {bases} bases ({Float64(perf_counter_ns() - t0) * 1e-9} s)",
            file=stderr,
        )
    var what = names[0] if len(names) == 1 else String(t"{len(names)} contigs")
    var surjector = Surjector(reference^, nodes^, sam_names.copy(), placement, windows^)
    var kind = GAF if gaf else (GAM if gam else AUTO)
    if not out_arg:
        var header = String("@HD\tVN:1.6\tSO:unsorted\n")
        for k in range(len(sam_names)):
            header += String(t"@SQ\tSN:{sam_names[k]}\tLN:{lengths[k]}\n")
        header += "@PG\tID:dromos-surject\tPN:dromos\tVN:0.1.0\n"
        var sink = SurjectSink[False](surjector^, threads, header, BamStaging(0, 1, memory, "."))
        try:
            stream_input(_input(pos[0]), threads, sink, kind)
        except e:
            sink.finish(threads, complete=False)
            raise e^
        sink.finish(threads, complete=True)
        sink.report(what)
        return
    var out = out_arg.value()
    var folder = dirname(out)
    var spill_dir = tmp_arg.value() if tmp_arg else (folder if folder else String("."))
    # Where the BAM and the spills go is tried now, not after the input has been read.
    _writable(out + ".tmp")
    _writable(spill_dir + "/.dromos-spill." + String(external_call["getpid", Int32]()))
    var text = String("@HD\tVN:1.6\tSO:coordinate\n")
    for k in range(len(sam_names)):
        text += String(t"@SQ\tSN:{sam_names[k]}\tLN:{lengths[k]}\n")
    text += "@PG\tID:dromos-surject\tPN:dromos\tVN:0.1.0\n"
    var header = bam_header(text, sam_names, lengths)
    var staging = BamStaging(len(names), threads, memory, spill_dir)
    var sink = SurjectSink[True](surjector^, threads, String(), staging^)
    try:
        stream_input(_input(pos[0]), threads, sink, kind)
    except e:
        sink.finish(threads, complete=False)
        raise e^
    sink.finish(threads, complete=True)
    sink.report(what)
    write_bam(sink.staging, header, BamIndex(sam_names, lengths, csi or needs_csi(lengths)), out, level, threads)


def prepare_main() raises:
    var args = _Args("prepare")
    if args.help():
        print(PREPARE_USAGE, end="")
        return
    var threads = args.threads()
    var sample = args.value("--reference-sample")
    var cache_arg = args.value("--cache")
    var graph = _graph(args, PREPARE_USAGE)
    _ = args.positional(0)
    var t0 = perf_counter_ns()
    var cache = graph_cache(graph, cache_arg, sample, threads)
    print(t"[dromos] the graph's cache: {cache} ({Float64(perf_counter_ns() - t0) * 1e-9} s)", file=stderr)


def catalog_main() raises:
    var args = _Args("catalog")
    if args.help():
        print(CATALOG_USAGE, end="")
        return
    var o = BuildOptions()
    var t = args.value("--threads")
    if not t:
        t = args.value("-t")
    o.threads = _number(t.value(), "--threads", 1, 1024, "catalog") if t else default_threads()
    o.graph = args.value("--graph").or_else(String())
    o.out = args.value("--out").or_else(String())
    o.dist = args.value("--dist")
    o.reference_sample = args.value("--reference-sample")
    while True:
        var r = args.value("--region")
        if not r:
            break
        o.regions.append(r.value())
    var jobs = args.value("--snarls-jobs")
    if jobs:
        o.snarls_jobs = _number(jobs.value(), "--snarls-jobs", 1, 1024, "catalog")
    o.snarls_cache = args.value("--snarls-cache")
    var vg = args.value("--vg")
    if vg:
        o.vg = vg.value()
    var level = args.value("--level")
    if level:
        o.level = _number(level.value(), "--level", 1, 22, "catalog")
    _ = args.positional(0)
    if not o.graph:
        print(CATALOG_USAGE, end="", file=stderr)
        exit(1)
    if o.dist and o.snarls_cache:
        _fail("--dist or --snarls-cache (the snarls' source), not both", "catalog")
    run_build(o)


def snarls_main() raises:
    var args = _Args("snarls")
    if args.help() or len(args.items) == 0:
        print(SNARLS_USAGE, end="")
        return
    var threads = args.threads()
    var stats = args.flag("--stats")
    var out = args.value("-o")
    var vg_file = args.value("--vg-file")
    var graph = args.value("--graph")
    _ = args.positional(0)
    if not out:
        _fail("-o OUT.tsv is needed", "snarls")
    if vg_file and graph:
        _fail("--graph or --vg-file, not both", "snarls")
    if vg_file:
        write_records_tsv(out.value(), read_snarls(vg_file.value(), threads))
        return
    if not graph:
        _fail("--graph or --vg-file is needed", "snarls")
    var t0 = perf_counter_ns()
    var trees = find_all_snarls(graph.value(), threads, stats)
    write_snarls_tsv(out.value(), trees, stats)
    var total = 0
    for k in range(len(trees)):
        total += len(trees[k])
    print(
        t"[dromos] {total} snarls in {len(trees)} component(s), {Float64(perf_counter_ns() - t0) * 1e-9} s", file=stderr
    )


def view_main() raises:
    var args = _Args("view")
    if args.help() or len(args.items) == 0:
        print(VIEW_USAGE, end="")
        return
    var threads = args.threads()
    var pos = args.positional_between(1, 2)
    view(pos[0], Optional(pos[1]) if len(pos) == 2 else None, threads)


def catalog_info_main() raises:
    var args = _Args("catalog-info")
    if args.help() or len(args.items) == 0:
        print(CATALOG_INFO_USAGE, end="")
        return
    info(args.positional(1)[0])


def reference_samples_main() raises:
    var args = _Args("reference-samples")
    if args.help():
        print(REFERENCE_SAMPLES_USAGE, end="")
        return
    if len(args.items) == 0:
        raise Error("usage: dromos reference-samples GRAPH.gbz")
    var graph = args.positional(1)[0]
    for s in reference_samples(graph):
        print(s)


def main():
    var args = argv()
    if len(args) < 2 or String(args[1]) == "-h" or String(args[1]) == "--help":
        print(OVERVIEW, end="")
        if len(args) < 2:
            exit(2)
        return
    var command = String(args[1])
    try:
        if command == "index":
            index_main()
        elif command == "gam-to-gaf":
            gam_to_gaf_main()
        elif command == "surject":
            surject_main()
        elif command == "catalog":
            catalog_main()
        elif command == "view":
            view_main()
        elif command == "catalog-info":
            catalog_info_main()
        elif command == "snarls":
            snarls_main()
        elif command == "prepare":
            prepare_main()
        elif command == "reference-samples":
            reference_samples_main()
        else:
            print(t"dromos: no command {command}", file=stderr)
            print(OVERVIEW, end="", file=stderr)
            exit(2)
    except e:
        print(t"dromos: {e}", file=stderr)
        exit(1)
