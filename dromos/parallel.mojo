"""Dependency-free thread pool primitive built on pthreads.

`parallel_for` runs `task.run(i)` for every `i` in `[0, count)` on up to
`num_threads` OS threads that pull indices from a shared atomic counter.
"""

from std.atomic import Atomic
from std.ffi import external_call
from std.os import abort
from std.sys import num_physical_cores


comptime MAX_THREADS = 1024
"""Thread counts are capped here: more only costs memory and setup."""


trait ParallelTask(Deinitable):
    def run(self, index: Int):
        ...


struct _Job[T: ParallelTask]:
    var task: Pointer[Self.T, MutUntrackedOrigin]
    var next: Atomic[Int64]
    var count: Int

    def __init__(out self, task: Pointer[Self.T, MutUntrackedOrigin], count: Int):
        self.task = task
        self.next = Atomic[Int64](0)
        self.count = count


def _worker[
    T: ParallelTask
](arg: MutOpaquePointer[MutUntrackedOrigin]) abi("C") -> MutOpaquePointer[MutUntrackedOrigin]:
    var job = arg.unsafe_bitcast[_Job[T]]()
    while True:
        var i = Int(job[].next.fetch_add(1))
        if i >= job[].count:
            break
        job[].task[].run(i)
    return arg


def default_threads() -> Int:
    return num_physical_cores()


def parallel_for[T: ParallelTask](mut task: T, count: Int, num_threads: Int):
    """Runs `task.run(i)` for `i` in `[0, count)` across `num_threads` threads."""
    var threads = min(num_threads, count, MAX_THREADS)
    if threads <= 1:
        for i in range(count):
            task.run(i)
        return
    var job = _Job[T](Pointer(to=task).unsafe_origin_cast[MutUntrackedOrigin](), count)
    var arg = Pointer(to=job).unsafe_origin_cast[MutUntrackedOrigin]().unsafe_bitcast[NoneType]()
    # pthread_t is an 8-byte handle on every supported platform.
    var handles = List[UInt64](length=threads, fill=0)
    var started = List[Bool](length=threads, fill=False)
    for t in range(threads):
        started[t] = external_call["pthread_create", Int32](Pointer(to=handles[t]), Int(0), _worker[T], arg) == 0
    # Threads that could not start leave their share to this one.
    var all_started = True
    for t in range(threads):
        all_started = all_started and started[t]
    if not all_started:
        _ = _worker[T](arg)
    for t in range(threads):
        if started[t]:
            _ = external_call["pthread_join", Int32](handles[t], Int(0))
    # A real read keeps `job` alive until every worker has joined; the workers
    # reach it through an untracked pointer.
    if Int(job.next.load()) < count:
        abort("parallel_for: workers stopped early")
