"""Dromos: pangenome graph alignments, fast. `dromos COMMAND -h` for each command's usage."""

from std.sys import argv, exit, stderr
from std.time import perf_counter_ns

from dromos.build import BuildOptions, run_build
from dromos.contig import Contig
from dromos.convert import GafSink
from dromos.graph import reference_samples, resolve_reference_sample
from dromos.index import IndexSink
from dromos.linearize import Window
from dromos.nodes import NodeTable
from dromos.parallel import default_threads
from dromos.pgsites import load_windows, view
from dromos.prepare import graph_cache
from dromos.refmap import RefMap
from dromos.stream import AUTO, GAF, GAM, stream_input
from dromos.surject import SurjectSink, Surjector

comptime OVERVIEW = """dromos: pangenome graph alignments, fast

  dromos index       a GAM or GAF placed on the reference, sorted and tabix-indexed, as pgphase reads it
  dromos gam-to-gaf  a GAM as GAF, the lines vg convert -G writes
  dromos surject     a GAM or GAF as SAM on one reference contig, the records pgcatalog surject writes

  dromos catalog            a GBZ graph's site catalogs, one per reference contig (<contig>.pgsites), as pgphase
                            reads them: the files pgcatalog build writes
  dromos view               a catalog as VCF, whole or a region of it
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

comptime SURJECT_USAGE = """dromos surject [options] --graph GRAPH.gbz --contig NAME IN > OUT.sam

IN's alignments as alignments to one reference contig, as SAM on stdout: the records pgcatalog surject writes for the
same graph, reads and options, in IN's order. Reads with no base on the contig are left out.
  | samtools sort -o OUT.bam && samtools index OUT.bam

  --graph FILE              the pangenome graph (GBZ)
  --cache FILE              the graph's cache: made there on first use (and again when the graph changes), read
                            after [GRAPH.gbz.dromos, beside the graph]
  IN                        a GAM, or a GAF with cs (and bq for qualities; plain, gzip or BGZF), told apart by their
                            first bytes; - reads stdin. dromos index's output keeps no cs and cannot be used
  --contig NAME             the reference contig (chr20)
  --reference-sample NAME   the sample whose path is the reference [the graph's declared reference sample, if it
                            declares exactly one]
  --sam-contig NAME         the contig's name in the SAM [<sample>#0#<contig>]
  --graph-placement         place read bases only where the graph puts them: a stretch of more than 50 bases off the
                            reference is not aligned to it, and such read ends are clipped
  --catalog FILE            the contig's catalog (.pgsites): a read spanning one of its snarl windows is aligned
                            there as its allele, at minimap2's asm5 costs
  -t THREADS                threads [the machine's cores]
  --gaf                     IN is a GAF (not told from its bytes)
  --gam                     IN is a GAM (not told from its bytes)

Each read is rebuilt from its walk and cs; a base on a node the contig visits once takes that node's position.
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

Site catalogs (DIR/<contig>.pgsites) from a GBZ graph, as pgphase reads them: the files pgcatalog build writes, byte
for byte. Every snarl the reference path passes through is kept, at every nesting level, with all its alleles.

  --graph FILE              the pangenome graph (GBZ)
  --out DIR                 where catalogs are written: <contig>.pgsites, or <contig>_<start>-<end>.pgsites for a
                            region [GRAPH.catalog, in the current folder]
  --reference-sample NAME   the sample whose paths are the reference: CHM13 for CHM13#0#chr20 [the graph's declared
                            reference sample, if it declares exactly one; HPRC graphs declare CHM13 and GRCh38 and
                            need this]
  --region REGION           a contig (chr20) or part of one (chr20:10000000-12000000, 1-based, inclusive),
                            repeatable [every reference contig]; a region holds the records a tabix query of the whole
                            contig's catalog returns
  -t, --threads N           threads [the machine's cores]
  --snarls-jobs N           vg snarls runs at a time [6, at most --threads]; fewer when memory is short (each needs
                            about 20 times its component's size)
  --snarls-cache DIR        where each contig's snarls are kept and reused by later runs on the same graph
                            [DIR/snarls]
  --vg PATH                 the vg binary, for vg chunk and vg snarls, run only for contigs whose snarls are not in
                            the cache [vg]
  --level N                 zstd level, 1-22 [6]; higher is smaller and slower

dromos view prints a catalog as VCF; for a tabix-indexed copy:
  dromos view chr20.pgsites | bgzip > chr20.sites.vcf.gz && tabix -p vcf chr20.sites.vcf.gz
"""

comptime VIEW_USAGE = """dromos view [-t THREADS] FILE.pgsites [CONTIG[:START-END]]

A catalog as VCF on stdout: its header, then every record, or the records overlapping a region (1-based, inclusive),
as tabix would return them from the VCF: what pgcatalog view prints. CONTIG is chr20 or CHM13#0#chr20; another
contig's region prints the header alone.

  -t THREADS  threads, inflating blocks [the machine's cores]
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


def surject_main() raises:
    var args = _Args("surject")
    if args.help():
        print(SURJECT_USAGE, end="")
        return
    var threads = args.threads()
    var contig_arg = args.value("--contig")
    var sample_arg = args.value("--reference-sample")
    var sam_arg = args.value("--sam-contig")
    var catalog_arg = args.value("--catalog")
    var cache_arg = args.value("--cache")
    var placement = args.flag("--graph-placement")
    var gaf = args.flag("--gaf")
    var gam = args.flag("--gam")
    if gaf and gam:
        _fail("--gaf or --gam, not both", "surject")
    if not contig_arg:
        _fail("--contig is required", "surject")
    var graph = _graph(args, SURJECT_USAGE)
    var pos = args.positional(1)
    var sample = resolve_reference_sample(graph, sample_arg)
    var cache = graph_cache(graph, cache_arg, sample, threads)
    var contig_name = contig_arg.value()
    var t0 = perf_counter_ns()
    var map = RefMap(cache, sample)
    var nodes = NodeTable(cache)
    var contig = Contig(map, nodes, contig_name)
    var windows = List[Window]()
    if catalog_arg:
        windows = load_windows(catalog_arg.value(), contig, map)
        print(t"[dromos] surject: {len(windows)} snarl windows", file=stderr)
    var sam_contig = sam_arg.value() if sam_arg else contig.chrom
    print(t"[dromos] {contig.chrom}: {len(contig.seq)} bases ({Float64(perf_counter_ns() - t0) * 1e-9} s)", file=stderr)
    var sink = SurjectSink(Surjector(contig^, nodes^, sam_contig, placement, windows^), threads)
    try:
        stream_input(_input(pos[0]), threads, sink, GAF if gaf else (GAM if gam else AUTO))
    except e:
        sink.finish(contig_name, threads, complete=False)
        raise e^
    sink.finish(contig_name, threads, complete=True)


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
    if not o.out:
        # The graph's file name without its last extension, then .catalog, in the current folder.
        var slash = o.graph.rfind("/")
        var name = String(o.graph[byte = slash + 1 :]) if slash >= 0 else o.graph
        var dot = name.rfind(".")
        o.out = (String(name[byte=:dot]) if dot > 0 else name) + ".catalog"
    run_build(o)


def view_main() raises:
    var args = _Args("view")
    if args.help() or len(args.items) == 0:
        print(VIEW_USAGE, end="")
        return
    var threads = args.threads()
    var pos = args.positional_between(1, 2)
    view(pos[0], Optional(pos[1]) if len(pos) == 2 else None, threads)


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
