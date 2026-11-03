"""A pool of worker threads shared by several jobs at once: `dromos catalog` builds several contigs at a time, each
contig's steps handed to the same workers, so one contig's tail (its last heavy snarls, its giant record) overlaps
the next contigs' work instead of leaving cores idle.

A job is a task and a count: the workers run `task.run(index, worker)` for each index, taking indices one at a time
from the job that ranks first (the lowest priority, then the lowest rank: the contig started first) among those with
indices left. `run` waits for its job; `submit` returns at once (a background job, `cancel`led or `wait`ed for
later). The job table is kept under a spin lock, taken only to pick an index and to finish one.
"""

from std.atomic import Atomic
from std.ffi import external_call

from dromos.parallel import MAX_THREADS

comptime MAX_JOBS = 256
"""Jobs in the table at once."""
comptime _Runner = def(Int, Int, Int) thin -> None
"""Runs index `i` of the task at an address on worker `w`."""


trait PoolTask:
    def run(self, index: Int, worker: Int):
        ...


def _run_one[T: PoolTask](task: Int, index: Int, worker: Int):
    MutPointer[T, MutUntrackedOrigin](unsafe_from_address=task)[].run(index, worker)


@fieldwise_init
struct _Start(Copyable, TrivialRegisterPassable):
    var pool: Int
    var worker: Int


@fieldwise_init
struct _Job(Copyable, Movable):
    var runner: _Runner
    var task: Int
    var count: Int
    var next: Int
    var done: Int
    var priority: Int
    var rank: Int
    var live: Bool


def _nothing(task: Int, index: Int, worker: Int):
    pass


struct Pool(Movable):
    """`threads` workers and the job table. Make it on the heap (OwnedPointer) and keep it there: the workers hold
    its address. `close` stops them once every job is done."""

    var jobs: List[_Job]
    var lock: Atomic[Int64]
    var stopping: Atomic[Int64]
    var threads: Int
    var handles: List[UInt64]
    var starts: List[_Start]
    """What each worker is started with: the pool and its number (kept here, so it outlives the worker)."""

    def __init__(out self, threads: Int):
        self.jobs = List[_Job]()
        for _ in range(MAX_JOBS):
            self.jobs.append(_Job(_nothing, 0, 0, 0, 0, 0, 0, False))
        self.lock = Atomic[Int64](0)
        self.stopping = Atomic[Int64](0)
        self.threads = max(min(threads, MAX_THREADS), 1)
        self.handles = List[UInt64]()
        self.starts = List[_Start]()

    def start(mut self):
        """Starts the workers (the pool must be where it stays)."""
        var at = Int(Pointer(to=self))
        self.starts = List[_Start](capacity=self.threads)
        for w in range(self.threads):
            self.starts.append(_Start(at, w))
        for w in range(self.threads):
            var h = UInt64(0)
            var arg = MutOpaquePointer[MutUntrackedOrigin](unsafe_from_address=Int(self.starts.unsafe_ptr().unsafe_offset(w)))
            if external_call["pthread_create", Int32](Pointer(to=h), Int(0), _worker, arg) == 0:
                self.handles.append(h)

    @always_inline
    def _acquire(mut self):
        var spins = 0
        while True:
            var expected = Int64(0)
            if self.lock.compare_exchange(expected, 1):
                return
            spins += 1
            if spins > 64:
                _ = external_call["sched_yield", Int32]()

    @always_inline
    def _release(mut self):
        self.lock.store(0)

    def submit[T: PoolTask](mut self, mut task: T, count: Int, priority: Int, rank: Int) -> Int:
        """Adds a job of `count` indices; returns its slot. The task must outlive the job (`wait` or `cancel`)."""
        while True:
            self._acquire()
            for s in range(MAX_JOBS):
                if not self.jobs[s].live:
                    self.jobs[s] = _Job(_run_one[T], Int(Pointer(to=task)), count, 0, 0, priority, rank, True)
                    self._release()
                    return s
            self._release()
            _ = external_call["usleep", Int32](UInt32(1000))

    def wait(mut self, slot: Int):
        """Waits until every index of the job in `slot` is done, and frees the slot."""
        var naps = 0
        while True:
            self._acquire()
            ref j = self.jobs[slot]
            if j.done >= j.count:
                j.live = False
                self._release()
                return
            self._release()
            naps += 1
            _ = external_call["usleep", Int32](UInt32(50 if naps < 100 else 200))

    def cancel(mut self, slot: Int):
        """No index of the job in `slot` not yet taken will be; waits for those taken, and frees the slot."""
        self._acquire()
        ref j = self.jobs[slot]
        j.done += j.count - j.next
        j.count = j.next
        self._release()
        self.wait(slot)

    def run[T: PoolTask](mut self, mut task: T, count: Int, priority: Int, rank: Int):
        """`task.run(index, worker)` for every index below `count`, on the pool's workers; returns when all are
        done."""
        if count <= 0:
            return
        var slot = self.submit(task, count, priority, rank)
        self.wait(slot)

    def _take(mut self) -> Tuple[_Runner, Int, Int, Int]:
        """The next index to run: (runner, task, index, slot); slot -1 when there is none."""
        self._acquire()
        var best = -1
        for s in range(MAX_JOBS):
            ref j = self.jobs[s]
            if j.live and j.next < j.count:
                if best < 0:
                    best = s
                else:
                    ref b = self.jobs[best]
                    if j.priority < b.priority or (j.priority == b.priority and j.rank < b.rank):
                        best = s
        if best < 0:
            self._release()
            return (_nothing, 0, 0, -1)
        ref j = self.jobs[best]
        var index = j.next
        j.next += 1
        var out = (j.runner, j.task, index, best)
        self._release()
        return out

    def _done(mut self, slot: Int):
        self._acquire()
        self.jobs[slot].done += 1
        self._release()

    def close(mut self):
        """Stops the workers (after every job) and waits for them."""
        self.stopping.store(1)
        for h in self.handles:
            _ = external_call["pthread_join", Int32](h, Int(0))
        self.handles.clear()

    def __deinit__(deinit self):
        self.stopping.store(1)
        for h in self.handles:
            _ = external_call["pthread_join", Int32](h, Int(0))


def _worker(arg: MutOpaquePointer[MutUntrackedOrigin]) abi("C") -> MutOpaquePointer[MutUntrackedOrigin]:
    var start = arg.unsafe_bitcast[_Start]()[]
    var pool = MutPointer[Pool, MutUntrackedOrigin](unsafe_from_address=start.pool)
    var worker = start.worker
    var naps = 0
    while True:
        var t = pool[]._take()
        if t[3] >= 0:
            t[0](t[1], t[2], worker)
            pool[]._done(t[3])
            naps = 0
            continue
        if pool[].stopping.load() != 0:
            break
        naps += 1
        _ = external_call["usleep", Int32](UInt32(20 if naps < 50 else 200))
    return arg
