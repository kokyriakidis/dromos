"""Dromos: pangenome graph alignments, fast. `dromos COMMAND -h` for each command's usage."""

from std.sys import argv, exit, stderr
from std.time import perf_counter_ns

from dromos.contig import Contig
from dromos.convert import GafSink
from dromos.graph import is_gbz, reference_samples, resolve_reference_sample
from dromos.graph_files import graph_node_table, graph_ref_map, write_node_table, write_ref_map
from dromos.index import IndexSink
from dromos.linearize import Window
from dromos.nodes import NodeTable
from dromos.parallel import default_threads
from dromos.pgsites import load_windows
from dromos.refmap import RefMap
from dromos.stream import AUTO, GAF, GAM, stream_input
from dromos.surject import SurjectSink, Surjector

comptime OVERVIEW = """dromos: pangenome graph alignments, fast

  dromos index       a GAM or GAF placed on the reference, sorted and tabix-indexed, as pgphase reads it
  dromos gam-to-gaf  a GAM as GAF, the lines vg convert -G writes
  dromos surject     a GAM or GAF as SAM on one reference contig, the records pgcatalog surject writes

  dromos ref-map            a GBZ graph's reference map, which index and surject place alignments by
  dromos node-table         a GBZ graph's node lengths and letters, which gam-to-gaf and surject read
  dromos reference-samples  the reference samples a GBZ graph declares

Where a command takes the graph's reference map or node table, the GBZ itself may be given: dromos makes the file
once and keeps it beside the graph (G.gbz.SAMPLE.refmap, G.gbz.nodes).

dromos COMMAND -h shows a command's usage. DROMOS_TIMING=1 reports where the time goes, on stderr.
"""

comptime INDEX_USAGE = """dromos index [-t THREADS] [-l LEVEL] [--reference-sample NAME] [--gaf | --gam] REFMAP IN OUT

IN's alignments placed on the reference, sorted, bgzipped and tabix-indexed as pgphase reads them: OUT and OUT.tbi,
the files pgphase's previous indexer (pgcat index-gaf --slim) wrote from the same alignments as GAF.

  REFMAP                   the graph's reference map (dromos ref-map), or the GBZ graph itself: its reference map
                           is made once and kept beside it as GRAPH.gbz.SAMPLE.refmap (beside OUT if the graph's
                           folder is not writable)
  IN                       a GAM, or a GAF (plain, gzip or BGZF), told apart by their first bytes; - reads stdin
  --reference-sample NAME  the sample whose paths are the reference: for a graph, the map to make [the graph's
                           declared reference sample, if it declares exactly one]; for a map, checked against its
  -t THREADS               threads [the machine's cores]
  -l LEVEL                 BGZF compression level, 0-12 [1]: 1 deflates about 3 times as fast as 6 into a file 40%
                           larger; 6 for a file to keep
  --gaf                    IN is a GAF (not told from its bytes)
  --gam                    IN is a GAM (not told from its bytes)
"""

comptime GAM_TO_GAF_USAGE = """dromos gam-to-gaf [-t THREADS] [--slim] NODES GAM

GAM as GAF on stdout: the lines vg convert -G writes for it, in the GAM's order.

  NODES       the graph's node table (dromos node-table; with --slim it needs no letters), or the GBZ graph itself:
              its node table is made once and kept beside it as GRAPH.gbz.nodes
  GAM         the alignments; - reads stdin
  -t THREADS  threads [the machine's cores]
  --slim      each line's 12 GAF columns alone (no cs, bq, AS or dv)
"""

comptime SURJECT_USAGE = """dromos surject [options] --contig NAME REFMAP NODES IN > OUT.sam
dromos surject [options] --contig NAME GRAPH.gbz IN > OUT.sam

IN's alignments as alignments to one reference contig, as SAM on stdout: the records pgcatalog surject writes for the
same reads and options (the graph's reference map and node table in place of the graph), in IN's order. Reads with no
base on the contig are left out.
  | samtools sort -o OUT.bam && samtools index OUT.bam

  REFMAP                    the graph's reference map (dromos ref-map)
  NODES                     the graph's node table, with letters (dromos node-table)
  GRAPH.gbz                 the GBZ graph, in place of either or both: its reference map and node table are made
                            once and kept beside it (GRAPH.gbz.SAMPLE.refmap, GRAPH.gbz.nodes)
  IN                        a GAM, or a GAF with cs (and bq for qualities; plain, gzip or BGZF), told apart by their
                            first bytes; - reads stdin. dromos index's output keeps no cs and cannot be used
  --contig NAME             the reference contig (chr20)
  --reference-sample NAME   the sample whose path is the reference; REFMAP's, which it must be [REFMAP's; for a
                            graph, its declared reference sample, if it declares exactly one]
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


comptime REF_MAP_USAGE = """dromos ref-map [-t THREADS] [--reference-sample NAME] --graph GRAPH.gbz --out FILE

Where each node of a GBZ graph lies on its reference sample's paths, saved once per graph for dromos index and
surject, which place alignments on the reference by it: the file pgcatalog ref-map writes.

  --graph FILE              the pangenome graph (GBZ)
  --out FILE                the reference map
  --reference-sample NAME   the sample whose paths are the reference: CHM13 for CHM13#0#chr20 [the graph's declared
                            reference sample, if it declares exactly one; HPRC graphs declare CHM13 and GRCh38 and
                            need this]
  -t THREADS                threads [the machine's cores]
"""

comptime NODE_TABLE_USAGE = """dromos node-table [-t THREADS] [--no-letters] --graph GRAPH.gbz --out FILE

A GBZ graph's node lengths and letters as one flat file, the node table dromos gam-to-gaf and surject read: the file
pgcatalog node-table writes.

  --graph FILE   the pangenome graph (GBZ)
  --out FILE     the node table
  --no-letters   lengths alone (enough for GAF's columns, not for its cs)
  -t THREADS     threads [the machine's cores]
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


def index_main() raises:
    var args = _Args("index")
    if args.help():
        print(INDEX_USAGE, end="")
        return
    var threads = args.threads()
    var level_arg = args.value("-l")
    var level = _number(level_arg.value(), "-l", 0, 12, "index") if level_arg else 1
    var sample = args.value("--reference-sample")
    var gaf = args.flag("--gaf")
    var gam = args.flag("--gam")
    if gaf and gam:
        _fail("--gaf or --gam, not both", "index")
    var pos = args.positional(3)
    var map_path = graph_ref_map(pos[0], sample, pos[2], threads) if is_gbz(pos[0]) else pos[0]
    var t0 = perf_counter_ns()
    var map = RefMap(map_path)
    if sample and sample.value() != map.sample:
        raise Error(t"{map_path} is the reference map of {map.sample}, not {sample.value()}")
    var sink = IndexSink(map^, threads)
    print(
        (
            t"[dromos] reference map loaded: {len(sink.map.contigs)} contigs of {sink.map.sample}"
            t" ({Float64(perf_counter_ns() - t0) * 1e-9} s)"
        ),
        file=stderr,
    )
    stream_input(_input(pos[1]), threads, sink, GAF if gaf else (GAM if gam else AUTO))
    sink.finish(pos[2], level, threads)


def gam_to_gaf_main() raises:
    var args = _Args("gam-to-gaf")
    if args.help():
        print(GAM_TO_GAF_USAGE, end="")
        return
    var threads = args.threads()
    var slim = args.flag("--slim")
    var pos = args.positional(2)
    var nodes = graph_node_table(pos[0], None, threads) if is_gbz(pos[0]) else pos[0]
    if slim:
        var sink = GafSink[False](NodeTable(nodes), threads)
        try:
            stream_input(_input(pos[1]), threads, sink, AUTO)
        finally:
            sink.finish()
    else:
        var sink = GafSink[True](NodeTable(nodes), threads)
        try:
            stream_input(_input(pos[1]), threads, sink, AUTO)
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
    var graph = args.flag("--graph-placement")
    var gaf = args.flag("--gaf")
    var gam = args.flag("--gam")
    if gaf and gam:
        _fail("--gaf or --gam, not both", "surject")
    if not contig_arg:
        _fail("--contig is required", "surject")
    var pos = args.positional_between(2, 3)
    if len(pos) == 2:
        if not is_gbz(pos[0]):
            _fail("3 arguments wanted (REFMAP NODES IN), or 2 with a GBZ graph first", "surject")
        pos.insert(1, pos[0])
    var map_path = graph_ref_map(pos[0], sample_arg, None, threads) if is_gbz(pos[0]) else pos[0]
    var nodes = graph_node_table(pos[1], None, threads) if is_gbz(pos[1]) else pos[1]
    var contig_name = contig_arg.value()
    var t0 = perf_counter_ns()
    var map = RefMap(map_path)
    if sample_arg and sample_arg.value() != map.sample:
        raise Error(t"{map_path} is the reference map of {map.sample}, not {sample_arg.value()}")
    var contig = Contig(map, NodeTable(nodes), contig_name)
    var windows = List[Window]()
    if catalog_arg:
        windows = load_windows(catalog_arg.value(), contig, map)
        print(t"[dromos] surject: {len(windows)} snarl windows", file=stderr)
    var sam_contig = sam_arg.value() if sam_arg else contig.chrom
    print(t"[dromos] {contig.chrom}: {len(contig.seq)} bases ({Float64(perf_counter_ns() - t0) * 1e-9} s)", file=stderr)
    var sink = SurjectSink(Surjector(contig^, NodeTable(nodes), sam_contig, graph, windows^), threads)
    try:
        stream_input(_input(pos[2]), threads, sink, GAF if gaf else (GAM if gam else AUTO))
    except e:
        sink.finish(contig_name, threads, complete=False)
        raise e^
    sink.finish(contig_name, threads, complete=True)


def _required(value: Optional[String], usage: StaticString) -> String:
    """A required option's value; without it, the usage on stderr (as pgcatalog does)."""
    if not value:
        print(usage, end="", file=stderr)
        exit(1)
    return value.value()


def ref_map_main() raises:
    var args = _Args("ref-map")
    if args.help():
        print(REF_MAP_USAGE, end="")
        return
    var threads = args.threads()
    var graph_arg = args.value("--graph")
    var out_arg = args.value("--out")
    var sample = args.value("--reference-sample")
    _ = args.positional(0)
    var graph = _required(graph_arg, REF_MAP_USAGE)
    var out = _required(out_arg, REF_MAP_USAGE)
    # The sample first: a graph declaring several fails at once, not after its load.
    var s = resolve_reference_sample(graph, sample)
    write_ref_map(graph, s, out, threads)


def node_table_main() raises:
    var args = _Args("node-table")
    if args.help():
        print(NODE_TABLE_USAGE, end="")
        return
    var threads = args.threads()
    var graph_arg = args.value("--graph")
    var out_arg = args.value("--out")
    var letters = not args.flag("--no-letters")
    _ = args.positional(0)
    var graph = _required(graph_arg, NODE_TABLE_USAGE)
    var out = _required(out_arg, NODE_TABLE_USAGE)
    write_node_table(graph, out, letters, threads)


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
        elif command == "ref-map":
            ref_map_main()
        elif command == "node-table":
            node_table_main()
        elif command == "reference-samples":
            reference_samples_main()
        else:
            print(t"dromos: no command {command}", file=stderr)
            print(OVERVIEW, end="", file=stderr)
            exit(2)
    except e:
        print(t"dromos: {e}", file=stderr)
        exit(1)
