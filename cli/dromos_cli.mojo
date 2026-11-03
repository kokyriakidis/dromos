"""Dromos: pangenome graph alignments, fast. `dromos COMMAND -h` for each command's usage."""

from std.sys import argv, exit, stderr
from std.time import perf_counter_ns

from dromos.convert import GafSink
from dromos.index import IndexSink
from dromos.nodes import NodeTable
from dromos.parallel import default_threads
from dromos.refmap import RefMap
from dromos.stream import stream_gaf, stream_gam

comptime OVERVIEW = """dromos: pangenome graph alignments, fast

  dromos index       a GAM or GAF placed on the reference, sorted and tabix-indexed, as pgphase reads it
  dromos gam-to-gaf  a GAM as GAF, the lines vg convert -G writes

dromos COMMAND -h shows a command's usage. DROMOS_TIMING=1 reports where the time goes, on stderr.
"""

comptime INDEX_USAGE = """dromos index [-t THREADS] [-l LEVEL] [--gaf] REFMAP IN OUT

IN's alignments placed on the reference, sorted, bgzipped and tabix-indexed as pgphase reads them: OUT and OUT.tbi,
the files pgcat index-gaf --slim writes from the same alignments as GAF.

  REFMAP      the graph's reference map (pgcat ref-map)
  IN          a GAM, or a GAF (plain, gzip or BGZF) when named *.gaf, *.gaf.gz, *.gaf.bgz or given --gaf;
              - reads stdin
  -t THREADS  threads [the machine's cores]
  -l LEVEL    BGZF compression level, 0-12 [1]: 1 deflates about 3 times as fast as 6 into a file 40% larger
              (as large as pgcat's, zlib at 6); 6 for a file to keep
  --gaf       IN is a GAF, whatever its name
"""

comptime GAM_TO_GAF_USAGE = """dromos gam-to-gaf [-t THREADS] [--slim] NODES GAM

GAM as GAF on stdout: the lines vg convert -G writes for it, in the GAM's order.

  NODES       the graph's node table (pgcat node-table; with --slim it needs no letters)
  GAM         the alignments; - reads stdin
  -t THREADS  threads [the machine's cores]
  --slim      each line's 12 GAF columns alone (no cs, bq, AS or dv)
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


def _gaf_name(path: String) -> Bool:
    """Whether `path` is named as a GAF."""
    for suffix in [".gaf", ".gaf.gz", ".gaf.bgz", ".gaf.bgzf"]:
        if path.endswith(suffix):
            return True
    return False


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
    var gaf = args.flag("--gaf")
    var pos = args.positional(3)
    var t0 = perf_counter_ns()
    var lines = gaf or _gaf_name(pos[1])
    var sink = IndexSink(RefMap(pos[0]), threads, lines)
    print(
        (
            t"[dromos] reference map loaded: {len(sink.map.contigs)} contigs of {sink.map.sample}"
            t" ({Float64(perf_counter_ns() - t0) * 1e-9} s)"
        ),
        file=stderr,
    )
    if lines:
        stream_gaf(_input(pos[1]), threads, sink)
    else:
        stream_gam(_input(pos[1]), threads, sink)
    sink.finish(pos[2], level, threads)


def gam_to_gaf_main() raises:
    var args = _Args("gam-to-gaf")
    if args.help():
        print(GAM_TO_GAF_USAGE, end="")
        return
    var threads = args.threads()
    var slim = args.flag("--slim")
    var pos = args.positional(2)
    if slim:
        var sink = GafSink[False](NodeTable(pos[0]), threads)
        try:
            stream_gam(_input(pos[1]), threads, sink)
        finally:
            sink.finish()
    else:
        var sink = GafSink[True](NodeTable(pos[0]), threads)
        try:
            stream_gam(_input(pos[1]), threads, sink)
        finally:
            sink.finish()


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
        else:
            print(t"dromos: no command {command}", file=stderr)
            print(OVERVIEW, end="", file=stderr)
            exit(2)
    except e:
        print(t"dromos: {e}", file=stderr)
        exit(1)
