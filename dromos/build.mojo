"""`dromos catalog` (pgcatalog's build.rs, ported): catalogs for every reference contig of a graph, for chosen contigs,
or for regions of them.

0. The reference paths come from the GBZ's path metadata, read alone.
1. Snarls come from each contig's connected component alone: one `vg chunk -C` pass extracts the components of the
   contigs not yet in the snarls cache, and `vg snarls` runs on several at once (it uses about one core each).
   Whole-genome snarls would make every contig pay for the genome: vg enumerates a snarl's haplotypes before
   checking which contig it is on.
2. The graph's parts are read once (the GBWT, and the node letters' pieces, unpacked as they are needed), then each
   selection is built (catalog.mojo).

The cache records the graph's fingerprint and refuses another graph's snarls.
"""

from std.ffi import OwnedDLHandle, external_call
from std.memory import OwnedPointer
from std.os import getenv, listdir, makedirs, remove, rmdir
from std.os.path import exists, isdir, isfile, islink
from std.sys import stderr
from std.sys.info import CompilationTarget
from std.time import perf_counter_ns

from dromos.bgzf import RawFile
from dromos.catalog import ContigPlan, ContigTail, Graph, ReferencePath, Summary, finish_catalog, run_steps
from dromos.graph import Metadata, fingerprint, reference_paths, resolve_reference_sample
from dromos.parallel import Background, ParallelTask
from dromos.stream import mut_untracked, untracked
from dromos.output import OutFile, write_all

comptime SNARLS_MEMORY_PER_BYTE = 22
"""vg snarls peaks at 16-21 times its component's file size (HPRC v2.1, all 25 CHM13 contigs)."""
comptime WORK_PREFIX = ".dromos-work-"


struct BuildOptions(Movable):
    var graph: String
    var out: String
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
        self.reference_sample = None
        self.regions = List[String]()
        self.snarls_jobs = None
        self.snarls_cache = None
        self.vg = String("vg")
        self.threads = 1
        self.level = 6


@fieldwise_init
struct Selection(Copyable, Movable):
    var contig: String
    var region: Optional[Tuple[Int, Int]]


struct _Log:
    var t0: Int

    def __init__(out self):
        self.t0 = perf_counter_ns()

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
) raises:
    """Computes the snarls of `contigs` into the cache, in a work folder under --out that is removed afterwards, also
    when a step fails (it holds the contigs' graphs: 24 GB for HPRC v2.1)."""
    remove_abandoned_work(o.out)
    var work = o.out + "/" + WORK_PREFIX + String(external_call["getpid", Int32]())
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


def run_build(o: BuildOptions) raises:
    var log = _Log()
    var sample = resolve_reference_sample(o.graph, o.reference_sample)
    if not o.reference_sample:
        log(t"reference sample {sample} (the only one the graph declares)")
    makedirs(o.out, exist_ok=True)
    var cache = o.snarls_cache.value() if o.snarls_cache else o.out + "/snarls"
    makedirs(cache, exist_ok=True)
    var id = fingerprint(o.graph)
    var id_file = cache + "/graph.id"
    var cached_id = _read_text(id_file)
    if cached_id:
        var c = String(cached_id.value().strip())
        if c != id:
            raise Error(t"{cache} holds snarls of another graph ({c}); use another --snarls-cache")
    else:
        var f = OutFile(id_file)
        f.put((id + "\n").as_bytes(), 0)
        f.close()

    var selections = List[Selection]()
    for r in o.regions:
        selections.append(parse_region(r, sample))
    var contigs = List[String]()
    for s in selections:
        var seen = False
        for c in contigs:
            seen = seen or c == s.contig
        if not seen:
            contigs.append(s.contig)
    sort(contigs)
    var need_paths = len(selections) == 0
    for c in contigs:
        need_paths = need_paths or not isfile(cache + "/" + c + ".snarls")
    if need_paths:
        var paths = reference_paths(o.graph, sample)
        if len(paths) == 0:
            raise Error(t"the graph has no {sample}#0#<contig> reference paths")
        if len(selections) == 0:
            var all = List[String]()
            for p in paths:
                all.append(p[0])
            sort(Span(all), natural_less)
            for c in all:
                selections.append(Selection(c, None))
            contigs = all^
        for c in contigs:
            var known = False
            for p in paths:
                known = known or p[0] == c
            if not known:
                raise Error(t"no reference path {sample}#0#{c} in the graph")
        var missing = List[String]()
        for c in contigs:
            if not isfile(cache + "/" + c + ".snarls"):
                missing.append(c)
        if len(missing) > 0:
            _compute_snarls(log, o, sample, paths, missing, cache)
    log(t"{len(selections)} selection(s) on {len(contigs)} contig(s); snarls in {cache}")

    var graph = Graph(o.graph, o.threads)
    log("graph loaded")
    var md = graph.parts.metadata()
    # Each selection is prepared (its reference path, snarls and traversals; its heavy snarls' helpers started)
    # while the one before is built, and finished (its giant records awaited, its last blocks written) while the one
    # after is: three at a time.
    var prepared = OwnedPointer(_Prepared())
    var preparing = Background(
        _PrepareTask(
            untracked(Pointer(to=graph)),
            untracked(Pointer(to=md)),
            sample,
            selections[0].copy(),
            cache,
            o.threads,
            mut_untracked(prepared.ptr()),
        ),
        1,
        1,
    )
    var finishing = OwnedPointer(_Finishing())
    var finisher = Optional[Background[_FinishTask]]()
    for i in range(len(selections)):
        preparing.join()
        if prepared[].error:
            raise Error(prepared[].error)
        var plan = prepared[].plan.take()
        if i + 1 < len(selections):
            prepared = OwnedPointer(_Prepared())
            preparing = Background(
                _PrepareTask(
                    untracked(Pointer(to=graph)),
                    untracked(Pointer(to=md)),
                    sample,
                    selections[i + 1].copy(),
                    cache,
                    o.threads,
                    mut_untracked(prepared.ptr()),
                ),
                1,
                1,
            )
        var tail = run_steps(graph, plan^, _prefix(o.out, selections[i]), o.level, o.threads)
        if finisher:
            finisher.value().join()
            _report(log, finishing[], o.out)
        finishing = OwnedPointer(_Finishing())
        finishing[].tail = tail^
        finishing[].selection = selections[i].copy()
        finisher = Background(_FinishTask(mut_untracked(finishing.ptr())), 1, 1)
    if finisher:
        finisher.value().join()
        _report(log, finishing[], o.out)
    # The tasks reach the graph and its metadata through untracked pointers: uses after the last join keep them alive
    # until every task is done.
    _ = graph.max_node
    _ = md^
    log("done")


def _prefix(folder: String, s: Selection) -> String:
    if s.region:
        return String(t"{folder}/{s.contig}_{s.region.value()[0]}-{s.region.value()[1]}")
    return String(t"{folder}/{s.contig}")


struct _Prepared(Movable):
    var plan: Optional[ContigPlan]
    var error: String

    def __init__(out self):
        self.plan = None
        self.error = String()


@fieldwise_init
struct _PrepareTask(Movable, ParallelTask):
    """A selection prepared on a thread of its own: its reference path walked, its plan made."""

    var graph: ImmPointer[Graph, ImmUntrackedOrigin]
    var md: ImmPointer[Metadata, ImmUntrackedOrigin]
    var sample: String
    var selection: Selection
    var cache: String
    var threads: Int
    var result: MutPointer[_Prepared, MutUntrackedOrigin]

    def run(self, index: Int):
        try:
            var t = perf_counter_ns()
            var r = ReferencePath(self.graph[], self.md[], self.sample, self.selection.contig, self.threads)
            if getenv("DROMOS_TIMING"):
                print(
                    t"dromos timing: {r.chrom}: reference path {Float64(perf_counter_ns() - t) * 1e-9} s", file=stderr
                )
            var snarls = self.cache + "/" + self.selection.contig + ".snarls"
            self.result[].plan = ContigPlan(self.graph[], r^, snarls, self.selection.region, self.threads)
        except e:
            self.result[].error = String(e)


struct _Finishing(Movable):
    var tail: Optional[ContigTail]
    var selection: Selection
    var summary: Optional[Summary]
    var error: String

    def __init__(out self):
        self.tail = None
        self.selection = Selection(String(), None)
        self.summary = None
        self.error = String()


@fieldwise_init
struct _FinishTask(Movable, ParallelTask):
    """A selection finished on a thread of its own."""

    var slot: MutPointer[_Finishing, MutUntrackedOrigin]

    def run(self, index: Int):
        try:
            self.slot[].summary = finish_catalog(self.slot[].tail.take())
        except e:
            self.slot[].error = String(e)


def _report(log: _Log, f: _Finishing, folder: String) raises:
    if f.error:
        raise Error(f.error)
    ref sm = f.summary.value()
    var also = " and the region" if f.selection.region else ""
    var secs = String(Int(sm.seconds * 10 + 0.5))
    if secs.byte_length() == 1:
        secs = "0" + secs
    secs = String(secs[byte = : secs.byte_length() - 1]) + "." + String(secs[byte = secs.byte_length() - 1 :])
    log(
        t"{_prefix(folder, f.selection)}.pgsites: {sm.records} records from {sm.snarls} snarls ({sm.on_reference} on"
        t" the reference{also}), {sm.blocks} blocks, {secs}s"
    )
