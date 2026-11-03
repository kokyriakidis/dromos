"""`dromos catalog` (pgcatalog's build.rs, ported, into one file): the graph's site catalog of every reference contig,
of chosen contigs, or of regions of them.

0. The reference paths come from the GBZ's path metadata, read alone.
1. Snarls come from the graph's distance index (dist.mojo; GRAPH.dist beside the graph unless told otherwise), checked
   against the graph, or from vg snarls of each contig's connected component, cached: one `vg chunk -C` pass extracts
   the components of the contigs not yet in the snarls cache, and `vg snarls` runs on several at once (it uses about
   one core each). Whole-genome snarls would make every contig pay for the genome: vg enumerates a snarl's haplotypes
   before checking which contig it is on. `SnarlSource` is the one interface to both.
2. The graph's parts are read once (the GBWT, and the node letters' pieces, unpacked as they are needed), then the
   contigs are built CONTIGS_AT_ONCE at a time, the biggest first, their steps on one pool of workers (catalog.mojo),
   and written into one file in their order (sites.mojo).

The snarls cache records the graph's fingerprint and refuses another graph's snarls.
"""

from std.atomic import Atomic
from std.ffi import OwnedDLHandle, external_call
from std.memory import OwnedPointer
from std.os import getenv, listdir, makedirs, remove, rmdir
from std.os.path import exists, isdir, isfile, islink
from std.sys import stderr
from std.sys.info import CompilationTarget
from std.time import perf_counter_ns

from dromos.bgzf import RawFile
from dromos.catalog import ContigPlan, Graph, ReferencePath, finish_catalog, run_steps
from dromos.dist import CHAIN, ORIENTATION, DistIndex
from dromos.graph import Metadata, fingerprint, graph_tags, reference_paths, resolve_reference_sample
from dromos.sites import (
    BLOCK_BP,
    BLOCK_BYTES,
    DEFAULT_LEVEL,
    DROMOS_VERSION,
    OWN_BLOCK_BYTES,
    ContigEntry,
    EncodedBlock,
    Md5,
    SitesWriter,
    put_oriented,
)
from dromos.snarls import Snarl, read_snarls
from dromos.parallel import Background, ParallelTask
from dromos.pool import Pool
from dromos.stream import mut_untracked, untracked
from dromos.output import OutFile, write_all

comptime SNARLS_MEMORY_PER_BYTE = 22
"""vg snarls peaks at 16-21 times its component's file size (HPRC v2.1, all 25 CHM13 contigs)."""
comptime WORK_PREFIX = ".dromos-work-"


struct BuildOptions(Movable):
    var graph: String
    var out: String
    var dist: Optional[String]
    var reference_sample: Optional[String]
    var regions: List[String]
    var snarls_jobs: Optional[Int]
    var snarls_cache: Optional[String]
    var vg: String
    var threads: Int
    var level: Int

    def __init__(out self):
        self.graph = String()
        self.out = String()
        self.dist = None
        self.reference_sample = None
        self.regions = List[String]()
        self.snarls_jobs = None
        self.snarls_cache = None
        self.vg = String("vg")
        self.threads = 1
        self.level = DEFAULT_LEVEL


@fieldwise_init
struct Selection(Copyable, Movable):
    var contig: String
    var region: Optional[Tuple[Int, Int]]


struct _Log:
    var t0: Int

    def __init__(out self):
        self.t0 = perf_counter_ns()

    def __init__(out self, t0: Int):
        self.t0 = t0

    def __call__[W: Writable](self, message: W):
        var s = Float64(perf_counter_ns() - self.t0) * 1e-9
        var t = String(Int(s * 10 + 0.5))
        if t.byte_length() == 1:
            t = "0" + t
        var secs = String(t[byte = : t.byte_length() - 1]) + "." + String(t[byte = t.byte_length() - 1 :])
        while secs.byte_length() < 7:
            secs = " " + secs
        print("[dromos ", secs, "s] ", message, sep="", file=stderr)


def _number(text: StringSlice, region: String) raises -> Int:
    var v = 0
    var digits = 0
    for b in text.as_bytes():
        if b == UInt8(ord(",")):
            continue
        if b < UInt8(ord("0")) or b > UInt8(ord("9")) or v > (Int.MAX - 9) // 10:
            raise Error(t"bad region {region}")
        v = v * 10 + Int(b - UInt8(ord("0")))
        digits += 1
    if digits == 0:
        raise Error(t"bad region {region}")
    return v


def parse_region(text: String, sample: String) raises -> Selection:
    """`chr20`, `chr20:10,000,000-12,000,000` or `CHM13#0#chr20:...`."""
    var prefix = sample + "#0#"
    var r = String(text[byte = prefix.byte_length() :]) if text.startswith(prefix) else text
    if r.find("#") >= 0:
        raise Error(t"region {r} is not on the reference sample {sample}'s paths")
    var colon = r.rfind(":")
    if colon < 0:
        return Selection(r, None)
    var contig = String(r[byte=:colon])
    var span = String(r[byte = colon + 1 :])
    var dash = span.find("-")
    if dash < 0:
        raise Error(t"bad region {r}: use contig:start-end")
    var a = _number(span[byte=:dash], r)
    var b = _number(span[byte = dash + 1 :], r)
    if a == 0 or b < a:
        raise Error(t"bad region {r}: start-end are 1-based and inclusive")
    return Selection(contig, (a, b))


def natural_less(a: String, b: String) -> Bool:
    """chr1 < chr2 < ... < chr10 < chrX: digits compared as numbers (pgcatalog's natural_key: runs of text each with
    the number after it; the last text's number, if none, counts as the largest)."""
    var x = _natural_key(a)
    var y = _natural_key(b)
    for i in range(min(len(x), len(y))):
        if x[i][0] != y[i][0]:
            return x[i][0] < y[i][0]
        if x[i][1] != y[i][1]:
            return x[i][1] < y[i][1]
    return len(x) < len(y)


def _natural_key(s: String) -> List[Tuple[String, UInt64]]:
    var out = List[Tuple[String, UInt64]]()
    var text = String()
    var num = String()
    for c in s.codepoint_slices():
        var b = c.as_bytes()[0]
        if c.byte_length() == 1 and b >= UInt8(ord("0")) and b <= UInt8(ord("9")):
            num += c
        else:
            if num.byte_length() > 0:
                out.append((text, _parse_u64(num)))
                text = String()
                num = String()
            text += c
    out.append((text, _parse_u64(num) if num.byte_length() > 0 else UInt64.MAX))
    return out^


def _parse_u64(digits: String) -> UInt64:
    """Rust's u64 parse of digits, 0 when it overflows (pgcatalog's unwrap_or(0))."""
    var v = UInt64(0)
    for b in digits.as_bytes():
        var d = UInt64(b - UInt8(ord("0")))
        if v > (UInt64.MAX - d) // 10:
            return 0
        v = v * 10 + d
    return v


def _meminfo(key: String) -> Optional[Int]:
    """/proc/meminfo's value for `key`, in bytes."""
    try:
        var text = open("/proc/meminfo", "r").read()
        for line in text.split("\n"):
            if line.startswith(key):
                var fields = line.split()
                if len(fields) >= 2:
                    return Int(fields[1]) * 1024
    except:
        pass
    return None


def _c(s: String) -> List[UInt8]:
    var c = List[UInt8](s.as_bytes())
    c.append(0)
    return c^


def _remove_tree(path: String):
    """Removes a file or folder and everything in it, ignoring failures."""
    try:
        if isdir(path) and not _is_link(path):
            for name in listdir(path):
                _remove_tree(path + "/" + name)
            rmdir(path)
        else:
            remove(path)
    except:
        pass


def _is_link(path: String) -> Bool:
    return islink(path)


def remove_abandoned_work(dir: String):
    """Removes the `WORK_PREFIX<pid>` work folders in `dir` of earlier runs that were killed (their process is gone),
    which would otherwise stay indefinitely: a build's contig graphs. Needs /proc to tell; elsewhere they are left
    alone."""
    if not exists("/proc/self"):
        return
    try:
        for name in listdir(dir):
            if not name.startswith(WORK_PREFIX):
                continue
            var pid = String(name[byte = WORK_PREFIX.byte_length() :])
            var ok = pid.byte_length() > 0
            for b in pid.as_bytes():
                ok = ok and b >= UInt8(ord("0")) and b <= UInt8(ord("9"))
            if ok and not exists("/proc/" + pid):
                print(t"[dromos] removing {dir}/{name}, left by an interrupted run", file=stderr)
                _remove_tree(dir + "/" + name)
    except:
        pass


comptime _O_WRONLY: Int32 = 1
comptime _O_CREAT: Int32 = 0x200 if CompilationTarget.is_macos() else 0o100
comptime _O_TRUNC: Int32 = 0x400 if CompilationTarget.is_macos() else 0o1000


@fieldwise_init
struct _Child(Copyable, Movable):
    var pid: Int32
    var what: String
    var need: Int


def _spawn(argv: List[String], stdout_path: String) raises -> Int32:
    """Starts `argv` (found on PATH, in this process's environment) with its stdout to `stdout_path`. Everything the
    child needs is made before the fork: it only opens, duplicates and execs."""
    var path = _c(stdout_path)
    var strings = List[List[UInt8]]()
    for a in argv:
        strings.append(_c(a))
    var pointers = List[Int]()
    for s in strings:
        pointers.append(Int(s.unsafe_ptr()))
    pointers.append(0)
    var pid = external_call["fork", Int32]()
    if pid == 0:
        comptime open_name = "open64" if CompilationTarget.is_linux() else "open"
        var fd = external_call[open_name, Int32](path.unsafe_ptr(), _O_WRONLY | _O_CREAT | _O_TRUNC, Int32(0o644))
        if fd < 0 or external_call["dup2", Int32](fd, Int32(1)) < 0:
            external_call["_exit", NoneType](Int32(127))
        _ = external_call["execvp", Int32](strings[0].unsafe_ptr(), pointers.unsafe_ptr())
        external_call["_exit", NoneType](Int32(127))
    _ = len(strings) + len(pointers) + len(path)
    if pid < 0:
        raise Error(t"running {argv[0]} {argv[1]}: cannot fork")
    return pid


def _finished(pid: Int32) -> Optional[Bool]:
    """Whether the child has exited (and succeeded); none while it runs."""
    var status = Int32(0)
    var got = external_call["waitpid", Int32](pid, Pointer(to=status), Int32(1))
    if got == 0:
        return None
    return got == pid and status == 0


def _wait(pid: Int32) -> Bool:
    var status = Int32(0)
    var got = external_call["waitpid", Int32](pid, Pointer(to=status), Int32(0))
    return got == pid and status == 0


def _copy_into(src: String, mut out: OutFile, mut at: Int) raises:
    var f = RawFile(src)
    var buf = List[UInt8](unsafe_uninit_length=min(f.size, 64 << 20))
    var done = 0
    while done < f.size:
        var n = f.read(buf.unsafe_ptr(), min(len(buf), f.size - done), done)
        if n <= 0:
            raise Error(t"cannot read {src}")
        out.put(Span(buf)[:n], at)
        at += n
        done += n


def _rename(a: String, b: String) raises:
    var x = _c(a)
    var y = _c(b)
    if external_call["rename", Int32](x.unsafe_ptr(), y.unsafe_ptr()) != 0:
        raise Error(t"cannot rename {a} to {b}")


def _compute_snarls_in(
    log: _Log,
    o: BuildOptions,
    sample: String,
    paths: List[Tuple[String, List[String]]],
    contigs: List[String],
    cache: String,
    work: String,
) raises:
    var names = List[String]()
    var by_contig = Dict[String, List[String]]()
    for p in paths:
        by_contig[p[0]] = p[1].copy()
    var list_text = String()
    for c in contigs:
        for name in by_contig[c]:
            list_text += name + "\n"
            names.append(name)
    var path_list = work + "/paths.txt"
    var lf = OutFile(path_list)
    lf.put(list_text.as_bytes(), 0)
    lf.close()
    log(t"extracting {len(contigs)} contig component(s) with vg chunk -C")
    var prefix = work + "/comp"
    var chunk: List[String] = [
        o.vg,
        "chunk",
        "-x",
        o.graph,
        "-C",
        "-P",
        path_list,
        "-O",
        "pg",
        "-t",
        String(o.threads),
        "-b",
        prefix,
    ]
    if not _wait(_spawn(chunk, "/dev/null")):
        raise Error(t"{o.vg} chunk failed")
    # The vg the snarls come from, for the catalog's header.
    var version: List[String] = [o.vg, "version"]
    if _wait(_spawn(version, work + "/vg.version")):
        var text = _read_text(work + "/vg.version")
        if text:
            var lines = text.value().split("\n")
            if len(lines) > 0:
                _write_text(cache + "/vg.version", String(lines[0]) + "\n")
    # Then each component's snarls, up to --snarls-jobs at once, a job starting only if the running jobs' and its
    # own predicted peaks fit in the memory available (one job always may).
    var ready = List[Tuple[String, String, Int]]()
    for c in contigs:
        for n in by_contig[c]:
            var f = prefix + "_" + n + ".vg"
            if isfile(f):
                ready.append((c, f, RawFile(f).size))
    var max_jobs = max(min(o.snarls_jobs.or_else(6), o.threads), 1)
    var avail = _meminfo("MemAvailable:")
    var budget = max(avail.value() - (2 << 30), 0) if avail else Int.MAX
    log(t"computing snarls of {len(ready)} component(s), up to {min(max_jobs, len(ready))} at a time")
    var per_job = 2 if o.threads >= 2 * max_jobs else 1
    var running = List[_Child]()
    var jobs = List[Tuple[String, String]]()  # contig, snarls part
    while len(ready) > 0 or len(running) > 0:
        var k = 0
        while k < len(running):
            var done = _finished(running[k].pid)
            if done:
                var child = running.pop(k)
                if not done.value():
                    for r in running:
                        _ = external_call["kill", Int32](r.pid, Int32(9))
                        _ = _wait(r.pid)
                    raise Error(t"{child.what} failed")
            else:
                k += 1
        while len(running) < max_jobs and len(ready) > 0:
            var need = ready[0][2] * SNARLS_MEMORY_PER_BYTE
            var committed = 0
            for r in running:
                committed += r.need
            if len(running) > 0 and committed + need > budget:
                break
            var job = ready.pop(0)
            var parts = 0
            for j in jobs:
                if j[0] == job[0]:
                    parts += 1
            var part = work + "/" + job[0] + "." + String(parts) + ".snarls"
            var argv: List[String] = [o.vg, "snarls", "-t", String(per_job), job[1]]
            running.append(_Child(_spawn(argv, part), String(t"vg snarls on {job[0]}"), need))
            jobs.append((job[0], part))
        _ = external_call["usleep", Int32](UInt32(500_000))
    # Fragments of a contig that share a component come out in one file, so a path may have none; a contig must have
    # one.
    for c in contigs:
        var any = False
        for j in jobs:
            any = any or j[0] == c
        if not any:
            raise Error(t"vg chunk wrote no component for {sample}#0#{c}")
    log(t"snarls of {len(jobs)} component(s) computed")
    # A contig's snarls file is its components' files one after another (both the stream format and gzip allow
    # concatenation).
    for c in contigs:
        var tmp = cache + "/" + c + ".snarls.tmp"
        var out = OutFile(tmp)
        var at = 0
        for j in jobs:
            if j[0] == c:
                _copy_into(j[1], out, at)
        out.close()
        _rename(tmp, cache + "/" + c + ".snarls")


def _compute_snarls(
    log: _Log,
    o: BuildOptions,
    sample: String,
    paths: List[Tuple[String, List[String]]],
    contigs: List[String],
    cache: String,
    folder: String,
) raises:
    """Computes the snarls of `contigs` into the cache, in a work folder beside the catalog that is removed
    afterwards, also when a step fails (it holds the contigs' graphs: 24 GB for HPRC v2.1)."""
    remove_abandoned_work(folder)
    var work = folder + "/" + WORK_PREFIX + String(external_call["getpid", Int32]())
    makedirs(work, exist_ok=True)
    try:
        _compute_snarls_in(log, o, sample, paths, contigs, cache, work)
    except e:
        _remove_tree(work)
        raise e^
    _remove_tree(work)


def _read_text(path: String) -> Optional[String]:
    try:
        return open(path, "r").read()
    except:
        return None


def _write_text(path: String, text: String) raises:
    var f = OutFile(path)
    f.put(text.as_bytes(), 0)
    f.close()


def _base_name(path: String) -> String:
    var slash = path.rfind("/")
    return String(path[byte = slash + 1 :]) if slash >= 0 else path


def _folder(path: String) -> String:
    var slash = path.rfind("/")
    if slash < 0:
        return String(".")
    return String(path[byte=:slash]) if slash > 0 else String("/")


def graph_stem(graph: String) -> String:
    """The graph's path without its .gbz."""
    return String(graph[byte = : graph.byte_length() - 4]) if graph.endswith(".gbz") else graph


def _selection_tag(s: Selection) -> String:
    if s.region:
        return String(t"{s.contig}_{s.region.value()[0]}-{s.region.value()[1]}")
    return s.contig


def catalog_path(o: BuildOptions, sample: String) -> String:
    """Where the catalog goes: GRAPH.SAMPLE.sites beside the graph (GRAPH.SAMPLE.chr20_1-1000.sites for regions, a
    tag each, joined by +); under --out if it is a folder (or ends with /), else at --out."""
    var name = graph_stem(o.graph) + "." + sample
    if len(o.regions) > 0:
        var tags = List[String]()
        for r in o.regions:
            try:
                tags.append(_selection_tag(parse_region(r, sample)))
            except:
                tags.append(r)
        name += "." + "+".join(tags)
    name += ".sites"
    if not o.out:
        return name
    if o.out.endswith("/") or isdir(o.out):
        return o.out.rstrip("/") + "/" + _base_name(name) if o.out != "/" else "/" + _base_name(name)
    return o.out


comptime VG_SNARLS = 0
comptime DIST_SNARLS = 1


struct SnarlSource(Movable):
    """Where a contig's snarls come from: vg snarls files (a folder of <contig>.snarls, made by vg chunk and vg
    snarls when missing), or the graph's distance index (dist.mojo). A snarls source computed in dromos itself would
    be a third kind here."""

    var kind: Int
    var cache: String
    var dist: Optional[DistIndex]

    def __init__(out self, cache: String):
        self.kind = VG_SNARLS
        self.cache = cache
        self.dist = None

    def __init__(out self, var dist: DistIndex):
        self.kind = DIST_SNARLS
        self.cache = String()
        self.dist = dist^

    def read(self, graph: Graph, r: ReferencePath, threads: Int) raises -> List[Snarl]:
        """The snarls of `r`'s contig (its connected component)."""
        var contig = String(r.chrom[byte = r.chrom.rfind("#") + 1 :])
        if self.kind == VG_SNARLS:
            return read_snarls(self.cache + "/" + contig + ".snarls", threads)
        ref d = self.dist.value()
        if len(r.steps) == 0:
            return List[Snarl]()
        var top = d.top(Int(r.steps[0] >> 1))
        var last = d.top(Int(r.steps[len(r.steps) - 1] >> 1))
        if last != top:
            raise Error(t"{d.path}: {r.chrom} runs through two components of the distance index: another graph's?")
        check_chain_edges(d, graph, top, r.chrom)
        var out = List[Snarl]()
        d.snarls(top, out)
        return out^

    def props(self, mut out: List[Tuple[String, String]]):
        if self.kind == DIST_SNARLS:
            ref d = self.dist.value()
            out.append(("snarls.source", "distance index"))
            out.append(("snarls.file", _base_name(d.path)))
            out.append(("snarls.dist.version", String(d.version)))
            out.append(("snarls.dist.bytes", String(d.file.size)))
            out.append(("snarls.orientation", ORIENTATION))
        else:
            out.append(("snarls.source", "vg snarls"))
            var v = _read_text(self.cache + "/vg.version")
            out.append(("snarls.vg", String(v.value().strip()) if v else String("unknown (snarls found cached)")))
            out.append(("snarls.orientation", "vg snarls'"))


def check_dist(d: DistIndex, graph: Graph) raises:
    """Refuses a distance index whose node ids are not the graph's: the same smallest and largest id and the same
    number of nodes."""
    var g = graph.parts.node_range()
    var count = d.node_count()
    if d.min_id != g[1] or d.max_id() != g[2] or count != g[0]:
        raise Error(
            t"{d.path} is not this graph's distance index: it holds {count} nodes, ids {d.min_id} to {d.max_id()};"
            t" the graph {g[0]} nodes, ids {g[1]} to {g[2]}"
        )


def check_chain_edges(d: DistIndex, graph: Graph, top: Int, chrom: String) raises:
    """Refuses a distance index whose top-level chain (the one `chrom` runs along) steps between nodes the graph has
    no edge between: each pair of consecutive chain nodes with no snarl between them is one of the graph's edges."""
    var words = List[UInt64]()
    var snarl_at = List[Int]()
    if d.kind(top) < CHAIN:
        return
    d.chain_words(top, words, snarl_at)
    var skip = List[Bool](length=len(words) + 1, fill=False)
    for i in range(len(snarl_at) // 2):
        skip[snarl_at[2 * i + 1]] = True
    var froms = List[UInt64]()
    var tos = List[UInt64]()
    var succ = List[UInt64]()
    var counts = List[Int64]()

    for k in range(1, len(words)):
        if skip[k]:
            continue
        froms.append(words[k - 1])
        tos.append(words[k])
        if len(froms) == 1 << 16:
            _check_edges(d, graph, chrom, froms, tos, succ, counts)
    if len(froms) > 0:
        _check_edges(d, graph, chrom, froms, tos, succ, counts)


def _check_edges(
    d: DistIndex,
    graph: Graph,
    chrom: String,
    mut froms: List[UInt64],
    mut tos: List[UInt64],
    mut succ: List[UInt64],
    mut counts: List[Int64],
) raises:
    """Each `froms[j]` has an edge to `tos[j]`; the lists are emptied."""
    graph.parts.edges(froms.unsafe_ptr(), len(froms), False, succ, counts)
    var at = 0
    for j in range(len(froms)):
        var n = Int(counts[j])
        if n < 0:
            raise Error(t"{d.path} is not this graph's distance index: node {froms[j] >> 1} is not in the graph")
        var found = False
        for x in range(n):
            found = found or succ[at + x] == tos[j]
        at += n
        if not found:
            raise Error(
                t"{d.path} is not this graph's distance index: its chain along {chrom} steps from node"
                t" {froms[j] >> 1} to {tos[j] >> 1}, which the graph does not join"
            )
    froms.clear()
    tos.clear()


def run_build(o: BuildOptions) raises:
    var log = _Log()
    var sample = resolve_reference_sample(o.graph, o.reference_sample)
    if not o.reference_sample:
        log(t"reference sample {sample} (the only one the graph declares)")
    var path = catalog_path(o, sample)
    var folder = _folder(path)
    makedirs(folder, exist_ok=True)
    var id = fingerprint(o.graph)

    var selections = List[Selection]()
    for r in o.regions:
        selections.append(parse_region(r, sample))
    var paths = reference_paths(o.graph, sample)
    if len(paths) == 0:
        raise Error(t"the graph has no {sample}#0#<contig> reference paths")
    var contigs = List[String]()
    if len(selections) == 0:
        for p in paths:
            contigs.append(p[0])
    else:
        for s in selections:
            var seen = False
            for c in contigs:
                seen = seen or c == s.contig
            if not seen:
                contigs.append(s.contig)
    sort(Span(contigs), natural_less)
    for c in contigs:
        var known = False
        for p in paths:
            known = known or p[0] == c
        if not known:
            raise Error(t"no reference path {sample}#0#{c} in the graph")
    # Each contig's regions: none when it is selected whole.
    var regions = List[List[Tuple[Int, Int]]]()
    for c in contigs:
        var whole = len(selections) == 0
        var rs = List[Tuple[Int, Int]]()
        for s in selections:
            if s.contig == c:
                if s.region:
                    rs.append(s.region.value())
                else:
                    whole = True
        regions.append(List[Tuple[Int, Int]]() if whole else rs^)

    # The snarls: the graph's distance index, or vg snarls (cached, or made now).
    var dist_path = o.dist
    if not dist_path and not o.snarls_cache and isfile(graph_stem(o.graph) + ".dist"):
        dist_path = graph_stem(o.graph) + ".dist"
    var source: SnarlSource
    if dist_path:
        source = SnarlSource(DistIndex(dist_path.value()))
        log(t"snarls from the distance index {dist_path.value()} (version {source.dist.value().version})")
    else:
        var cache = o.snarls_cache.value() if o.snarls_cache else graph_stem(folder + "/" + _base_name(o.graph))
        if not o.snarls_cache:
            cache += ".snarls"
        makedirs(cache, exist_ok=True)
        var id_file = cache + "/graph.id"
        var cached_id = _read_text(id_file)
        if cached_id:
            var c = String(cached_id.value().strip())
            if c != id:
                raise Error(t"{cache} holds snarls of another graph ({c}); use another --snarls-cache")
        else:
            _write_text(id_file, id + "\n")
        var missing = List[String]()
        for c in contigs:
            if not isfile(cache + "/" + c + ".snarls"):
                missing.append(c)
        if len(missing) > 0:
            _compute_snarls(log, o, sample, paths, missing, cache, folder)
        source = SnarlSource(cache)
        log(t"snarls from vg snarls in {cache}")
    log(t"{len(contigs)} contig(s); writing {path}")

    var graph = Graph(o.graph, o.threads)
    log("graph loaded")
    if source.kind == DIST_SNARLS:
        check_dist(source.dist.value(), graph)
    var md = graph.parts.metadata()
    var props = List[Tuple[String, String]]()
    var nodes = graph.parts.node_range()
    props.append(("dromos.version", DROMOS_VERSION))
    props.append(("graph.file", _base_name(o.graph)))
    props.append(("graph.fingerprint", id))
    props.append(("graph.nodes", String(nodes[0])))
    props.append(("graph.min_node", String(nodes[1])))
    props.append(("graph.max_node", String(nodes[2])))
    for t in graph_tags(o.graph):
        props.append((t[0] + "." + t[1], t[2]))
    props.append(("reference.sample", sample))
    source.props(props)
    props.append(("build.block_bp", String(BLOCK_BP)))
    props.append(("build.block_bytes", String(BLOCK_BYTES)))
    props.append(("build.own_block_bytes", String(OWN_BLOCK_BYTES)))
    props.append(("build.level", String(o.level)))
    var entries = List[ContigEntry]()
    for k in range(len(contigs)):
        entries.append(
            ContigEntry(
                contigs[k], String(t"{sample}#0#{contigs[k]}"), 0, List[UInt8](), regions[k].copy(), 0, 0, 0, 0
            )
        )
    var writer = OwnedPointer(SitesWriter(path, props^, entries^))
    run_contigs(log, graph, md, sample, contigs, regions, source, writer[], o.level, o.threads)
    var size = writer[].finish()
    log(t"{path}: {size} bytes")
    _ = graph.max_node
    _ = md^
    log("done")


def reference_md5(graph: Graph, r: ReferencePath) raises -> List[UInt8]:
    """MD5 of the contig's reference sequence (SAM's M5): its path's letters in capitals at their positions, N
    where no fragment lies."""
    var seq = List[UInt8](length=r.length, fill=UInt8(ord("N")))
    var p = seq.unsafe_ptr()
    for k in range(len(r.steps)):
        var v = r.steps[k]
        var s = graph.letters(Int(v >> 1))
        put_oriented(p.unsafe_offset(r.pos[k]), s, v & 1 == 1)
    var h = Md5()
    h.update(Span(seq))
    return h.digest()


comptime CONTIGS_AT_ONCE = 3
"""Contigs built at a time ($DROMOS_CONTIGS_AT_ONCE to change): their steps share the pool's workers."""


def contig_weight(graph: Graph, md: Metadata, sample: String, contig: String, source: SnarlSource) raises -> Int:
    """How much work a contig is, roughly, to start the biggest first: its chain's nodes in the distance index, or
    its vg snarls file's size."""
    if source.kind == VG_SNARLS:
        return RawFile(source.cache + "/" + contig + ".snarls").size
    var sample_id = md.sample_id(sample)
    var contig_id = md.contig_id(contig)
    if not sample_id or not contig_id:
        return 0
    var names = md.path_names()
    var best = -1
    var best_fragment = Int.MAX
    for id in range(len(names)):
        var p = names[id]
        if p.sample == sample_id.value() and p.contig == contig_id.value() and p.phase == 0:
            if p.fragment < best_fragment:
                best = id
                best_fragment = p.fragment
    if best < 0:
        return 0
    var v = graph.parts.path_first(best)
    ref d = source.dist.value()
    if v == 0 or d.node_record(Int(v >> 1)) == 0:
        return 0
    var top = d.top(Int(v >> 1))
    return d.entry(top + 1) if d.kind(top) >= CHAIN else 0


struct _Shared(Movable):
    """What every contig's driver reads, and where it writes."""

    var graph: ImmPointer[Graph, ImmUntrackedOrigin]
    var md: ImmPointer[Metadata, ImmUntrackedOrigin]
    var source: ImmPointer[SnarlSource, ImmUntrackedOrigin]
    var sample: String
    var contigs: List[String]
    var regions: List[List[Tuple[Int, Int]]]
    var order: List[Int]
    """Contig indexes, the biggest first."""
    var writer: MutPointer[SitesWriter, MutUntrackedOrigin]
    var pool: MutPointer[Pool, MutUntrackedOrigin]
    var lock: Atomic[Int64]
    """Held while a driver hands its contig to the writer."""
    var failed: Atomic[Int64]
    var errors: List[String]
    var level: Int
    var threads: Int
    var t0: Int

    def __init__(
        out self,
        graph: ImmPointer[Graph, ImmUntrackedOrigin],
        md: ImmPointer[Metadata, ImmUntrackedOrigin],
        source: ImmPointer[SnarlSource, ImmUntrackedOrigin],
        sample: String,
        var contigs: List[String],
        var regions: List[List[Tuple[Int, Int]]],
        var order: List[Int],
        writer: MutPointer[SitesWriter, MutUntrackedOrigin],
        pool: MutPointer[Pool, MutUntrackedOrigin],
        level: Int,
        threads: Int,
        t0: Int,
    ):
        self.graph = graph
        self.md = md
        self.source = source
        self.sample = sample
        self.contigs = contigs^
        self.regions = regions^
        self.errors = List[String](length=len(order), fill=String())
        self.order = order^
        self.writer = writer
        self.pool = pool
        self.lock = Atomic[Int64](0)
        self.failed = Atomic[Int64](0)
        self.level = level
        self.threads = threads
        self.t0 = t0


@fieldwise_init
struct _ContigTask(Movable, ParallelTask):
    """A contig's driver, on a thread of its own: the contig at `order[index]` prepared (its reference path, the MD5
    of its sequence, its snarls, its plan), built in steps on the pool, finished and handed to the writer."""

    var sh: MutPointer[_Shared, MutUntrackedOrigin]

    def run(self, index: Int):
        ref sh = self.sh[]
        if sh.failed.load() != 0:
            return
        var k = sh.order[index]
        try:
            var t = perf_counter_ns()
            ref graph = sh.graph[]
            var r = ReferencePath(graph, sh.md[], sh.sample, sh.contigs[k], sh.threads)
            if getenv("DROMOS_TIMING"):
                print(
                    t"dromos timing: {r.chrom}: reference path {Float64(perf_counter_ns() - t) * 1e-9} s", file=stderr
                )
            var chrom = r.chrom
            var length = r.length
            var md5 = reference_md5(graph, r)
            var snarls = sh.source[].read(graph, r, sh.threads)
            var plan = ContigPlan(graph, r^, snarls^, sh.regions[k], sh.threads, sh.pool, index)
            var done = finish_catalog(run_steps(graph, plan^, sh.level, sh.threads))
            var blocks = List[EncodedBlock]()
            swap(blocks, done[1])
            ref sm = done[0]
            while True:
                var expected = Int64(0)
                if sh.lock.compare_exchange(expected, 1):
                    break
                _ = external_call["usleep", Int32](UInt32(100))
            try:
                sh.writer[].complete(k, blocks^, md5, sm.records, length)
            finally:
                sh.lock.store(0)
            var log = _Log(sh.t0)
            log(
                t"{chrom}: {sm.records} records from {sm.snarls} snarls ({sm.on_reference} on the reference),"
                t" {sm.blocks} blocks, {_tenths(sm.seconds)}s"
            )
        except e:
            sh.errors[index] = String(t"{sh.contigs[k]}: {e}")
            sh.failed.store(1)


def _tenths(seconds: Float64) -> String:
    var t = String(Int(seconds * 10 + 0.5))
    if t.byte_length() == 1:
        t = "0" + t
    return String(t[byte = : t.byte_length() - 1]) + "." + String(t[byte = t.byte_length() - 1 :])


def run_contigs(
    log: _Log,
    graph: Graph,
    md: Metadata,
    sample: String,
    contigs: List[String],
    regions: List[List[Tuple[Int, Int]]],
    source: SnarlSource,
    mut writer: SitesWriter,
    level: Int,
    threads: Int,
) raises:
    """Every contig built, CONTIGS_AT_ONCE at a time, the biggest first, each by a driver thread of its own; their
    parallel steps run on one pool of `threads` workers, which takes a job of a contig started earlier first. The
    writer puts the contigs in the file in their order, whatever order they finish in."""
    var weights = List[Int]()
    for c in contigs:
        weights.append(contig_weight(graph, md, sample, c, source))
    var order = List[Int]()
    for k in range(len(contigs)):
        order.append(k)

    def heavier(a: Int, b: Int) {weights} -> Bool:
        return weights[a] > weights[b] or (weights[a] == weights[b] and a < b)

    sort(Span(order), heavier)
    var at_once = CONTIGS_AT_ONCE
    var env = getenv("DROMOS_CONTIGS_AT_ONCE")
    if env:
        at_once = max(Int(env), 1)
    var pool = OwnedPointer(Pool(threads))
    pool[].start()
    var sh = OwnedPointer(
        _Shared(
            untracked(Pointer(to=graph)),
            untracked(Pointer(to=md)),
            untracked(Pointer(to=source)),
            sample,
            contigs.copy(),
            regions.copy(),
            order^,
            mut_untracked(Pointer(to=writer)),
            mut_untracked(pool.ptr()),
            level,
            threads,
            log.t0,
        )
    )
    var drivers = Background(_ContigTask(mut_untracked(sh.ptr())), len(contigs), min(at_once, len(contigs)))
    if len(drivers.handles) == 0:
        raise Error("cannot start a thread")
    drivers.join()
    pool[].close()
    for e in sh[].errors:
        if e:
            raise Error(e)
    # The drivers reach these through untracked pointers: uses after the join keep them alive until they are done.
    _ = graph.max_node
    _ = source.kind
    _ = writer.next
