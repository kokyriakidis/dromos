"""Dromos: pangenome graph alignments, fast.

    dromos gam-to-gaf [-t THREADS] [--slim] NODES GAM       GAM as GAF on stdout, the lines vg convert -G writes
    dromos index [-t THREADS] [-l LEVEL] REFMAP GAM OUT     GAM as pgphase's index: OUT and OUT.tbi

NODES is the graph's node table (`pgcat node-table`; see dromos/nodes.mojo); REFMAP its reference map
(`pgcat ref-map`; see dromos/refmap.mojo).
"""

from std.sys import argv, exit, stderr
from std.time import perf_counter_ns

from dromos.convert import GafSink
from dromos.index import IndexSink
from dromos.nodes import NodeTable
from dromos.parallel import default_threads
from dromos.refmap import RefMap
from dromos.stream import stream_gam

comptime USAGE = """dromos gam-to-gaf [-t THREADS] [--slim] NODES GAM
  GAM as GAF on stdout, the lines vg convert -G writes for it (NODES: the graph's node table, pgcat node-table)
  --slim   each line's 12 GAF columns alone (no cs, bq, AS, dv; NODES needs no letters)

dromos index [-t THREADS] [-l LEVEL] REFMAP GAM OUT
  GAM placed on the reference, sorted, bgzipped and tabix-indexed as pgphase reads it: OUT and OUT.tbi, the files
  pgcat index-gaf --slim writes from its GAF (REFMAP: the graph's reference map, pgcat ref-map)
  -l LEVEL   BGZF compression level, 0-12 [6]
"""


def _number(text: String, what: StaticString) -> Int:
    try:
        return Int(text)
    except:
        print(t"dromos: {what} needs a number", file=stderr)
        exit(2)
    return 0


def main():
    var args = argv()
    var threads = default_threads()
    var slim = False
    var level = 6
    var positional = List[String]()
    var i = 1
    while i < len(args):
        var a = String(args[i])
        if a == "-t" and i + 1 < len(args):
            threads = max(1, _number(String(args[i + 1]), "-t"))
            i += 2
            continue
        if a == "-l" and i + 1 < len(args):
            level = min(12, max(0, _number(String(args[i + 1]), "-l")))
            i += 2
            continue
        if a == "--slim":
            slim = True
            i += 1
            continue
        if a == "-h" or a == "--help":
            print(USAGE, end="")
            return
        positional.append(a)
        i += 1
    try:
        if len(positional) == 3 and positional[0] == "gam-to-gaf":
            if slim:
                var sink = GafSink[False](NodeTable(positional[1]), threads)
                try:
                    stream_gam(positional[2], threads, sink)
                finally:
                    sink.finish()
            else:
                var sink = GafSink[True](NodeTable(positional[1]), threads)
                try:
                    stream_gam(positional[2], threads, sink)
                finally:
                    sink.finish()
            return
        if len(positional) == 4 and positional[0] == "index":
            var t0 = perf_counter_ns()
            var sink = IndexSink(RefMap(positional[1]), threads)
            print(
                (
                    t"[dromos] reference map loaded: {len(sink.map.contigs)} contigs of {sink.map.sample}"
                    t" ({Float64(perf_counter_ns() - t0) * 1e-9} s)"
                ),
                file=stderr,
            )
            stream_gam(positional[2], threads, sink)
            sink.finish(positional[3], level, threads)
            return
    except e:
        print(t"dromos: {e}", file=stderr)
        exit(1)
    print(USAGE, end="", file=stderr)
    exit(2)
