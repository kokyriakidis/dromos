"""Dromos: pangenome graph alignments, fast.

    dromos gam-to-gaf [-t THREADS] NODES GAM    GAM as GAF on stdout, the lines vg convert -G writes

NODES is the graph's node table (`pgcat node-table`; see dromos/nodes.mojo).
"""

from std.math import ceildiv
from std.memory import unsafe_memmove
from std.sys import argv, exit, stderr

from dromos.bgzf import BgzfReader
from dromos.gam import Alignment, frame, gaf_line
from dromos.nodes import NodeTable
from dromos.output import OutputBuffer
from dromos.parallel import ParallelTask, default_threads, parallel_for

comptime USAGE = """dromos gam-to-gaf [-t THREADS] NODES GAM
  GAM as GAF on stdout, the lines vg convert -G writes for it (NODES: the graph's node table, pgcat node-table)
"""


struct _GafTask[data_origin: ImmOrigin, messages_origin: ImmOrigin, nodes_origin: ImmOrigin, out_origin: MutOrigin](
    ParallelTask
):
    """Writes the GAF lines of messages `index * per ..` into text `index`."""

    var data: ImmPointer[UInt8, Self.data_origin]
    var messages: ImmPointer[Int, Self.messages_origin]
    var count: Int
    var per: Int
    var nodes: ImmPointer[NodeTable, Self.nodes_origin]
    var texts: MutPointer[List[UInt8], Self.out_origin]
    var failed: MutPointer[Int, MutUntrackedOrigin]

    def __init__(
        out self,
        data: ImmPointer[UInt8, Self.data_origin],
        messages: ImmPointer[Int, Self.messages_origin],
        count: Int,
        per: Int,
        nodes: ImmPointer[NodeTable, Self.nodes_origin],
        texts: MutPointer[List[UInt8], Self.out_origin],
        failed: MutPointer[Int, MutUntrackedOrigin],
    ):
        self.data = data
        self.messages = messages
        self.count = count
        self.per = per
        self.nodes = nodes
        self.texts = texts
        self.failed = failed

    def run(self, index: Int):
        ref text = self.texts[unsafe_offset=index]
        text.clear()
        var a = Alignment()
        try:
            for i in range(index * self.per, min(self.count, (index + 1) * self.per)):
                a.decode(self.data, self.messages[unsafe_offset=2 * i], self.messages[unsafe_offset=2 * i + 1])
                gaf_line(self.data, a, self.nodes[], text)
        except:
            self.failed[unsafe_offset=0] = 1


def gam_to_gaf(nodes_path: String, gam_path: String, threads: Int) raises:
    var nodes = NodeTable(nodes_path)
    var reader = BgzfReader(gam_path, threads)
    var out = OutputBuffer(1)
    out.write("@HD\tVN:Z:1.0\n")
    var buf = List[UInt8]()
    var slices = 4 * threads
    var texts = List[List[UInt8]](length=slices, fill=List[UInt8]())
    var failed = List[Int](length=1, fill=0)
    while reader.next(buf):
        var messages = List[Int]()
        var used = frame(Span(buf), messages)
        var count = len(messages) // 2
        if count > 0:
            var per = ceildiv(count, slices)
            var task = _GafTask(
                buf.unsafe_ptr().as_imm(),
                messages.unsafe_ptr().as_imm(),
                count,
                per,
                Pointer(to=nodes),
                texts.unsafe_ptr(),
                failed.unsafe_ptr().unsafe_origin_cast[MutUntrackedOrigin](),
            )
            var tasks = ceildiv(count, per)
            parallel_for(task, tasks, threads)
            if failed[0] != 0:
                raise Error(t"{gam_path}: a GAM message could not be read")
            for t in range(tasks):
                out.write_bytes(Span(texts[t]))
        # The unfinished group moves to the front, to be read with the next chunk.
        var rest = len(buf) - used
        if rest > 0 and used > 0:
            var p = buf.unsafe_ptr().unsafe_origin_cast[MutUntrackedOrigin]()  # memmove allows the overlap
            unsafe_memmove(dest=p, src=p.unsafe_offset(used).as_imm(), count=rest)
        buf.resize(unsafe_uninit_length=rest)
    if len(buf) > 0:
        raise Error(t"{gam_path}: ends inside a GAM group")
    out.flush()


def main():
    var args = argv()
    var threads = default_threads()
    var positional = List[String]()
    var i = 1
    while i < len(args):
        var a = String(args[i])
        if a == "-t" and i + 1 < len(args):
            try:
                threads = max(1, Int(String(args[i + 1])))
            except:
                print("dromos: -t needs a number", file=stderr)
                exit(2)
            i += 2
            continue
        if a == "-h" or a == "--help":
            print(USAGE, end="")
            return
        positional.append(a)
        i += 1
    if len(positional) != 3 or positional[0] != "gam-to-gaf":
        print(USAGE, end="", file=stderr)
        exit(2)
    try:
        gam_to_gaf(positional[1], positional[2], threads)
    except e:
        print(t"dromos: {e}", file=stderr)
        exit(1)
