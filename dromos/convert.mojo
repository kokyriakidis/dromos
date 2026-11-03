"""`dromos gam-to-gaf`: a GAM's alignments as the GAF lines `vg convert -G` writes for them, in the input's order."""

from dromos.gam import Alignment, gaf_line, slim_line
from dromos.nodes import NodeTable
from dromos.output import write_all
from dromos.stream import MessageSink, mut_untracked


struct GafSink[tags: Bool](MessageSink, Movable):
    """Each slice writes its messages' lines into a text of its own. The step after, a side task writes them out in
    order (on stdout) while the next chunk's lines are made: two texts a slice, taking turns."""

    var nodes: NodeTable
    var texts: List[List[UInt8]]
    var prev: List[List[UInt8]]
    """The last step's texts, being written."""
    var pending: Int
    """How many of `prev` are still to write."""
    var alignments: List[Alignment]
    var walks: List[List[Int]]

    def __init__(out self, var nodes: NodeTable, threads: Int):
        self.nodes = nodes^
        self.texts = List[List[UInt8]](length=4 * threads, fill=List[UInt8]())
        self.prev = List[List[UInt8]](length=4 * threads, fill=List[UInt8]())
        self.pending = 0
        self.alignments = List[Alignment]()
        for _ in range(4 * threads):
            self.alignments.append(Alignment())
        self.walks = List[List[Int]](length=4 * threads, fill=List[Int]())
        var header = "@HD\tVN:Z:1.0\n".as_bytes()
        write_all(1, header.unsafe_ptr(), len(header))

    def begin(mut self, lines: Bool) raises:
        if lines:
            raise Error("a GAF, not a GAM (gam-to-gaf reads a GAM)")

    def slices(self) -> Int:
        return len(self.texts)

    def take(
        mut self,
        slice: Int,
        data: ImmPointer[UInt8, ImmUntrackedOrigin],
        messages: ImmPointer[Int, ImmUntrackedOrigin],
        first: Int,
        last: Int,
    ) raises:
        var text = mut_untracked(Pointer(to=self.texts[slice]))
        var a = mut_untracked(Pointer(to=self.alignments[slice]))
        var walk = mut_untracked(Pointer(to=self.walks[slice]))
        text[].clear()
        for m in range(first, last):
            comptime if Self.tags:
                a[].decode(data, messages[unsafe_offset=2 * m], messages[unsafe_offset=2 * m + 1])
                gaf_line(data, a[], self.nodes, text[])
            else:
                a[].decode[slim=True](data, messages[unsafe_offset=2 * m], messages[unsafe_offset=2 * m + 1])
                _ = slim_line(data, a[], self.nodes, text[], walk[])

    def side_tasks(self) -> Int:
        return 1 if self.pending > 0 else 0

    def side_task(mut self, index: Int) raises:
        self._write_pending()

    def _write_pending(mut self):
        for t in range(self.pending):
            write_all(1, self.prev[t].unsafe_ptr(), len(self.prev[t]))
        self.pending = 0

    def after_step(mut self, slices: Int) raises:
        # `prev` was written in the step just run; the step's texts are the next to write.
        for t in range(slices):
            swap(self.texts[t], self.prev[t])
        self.pending = slices

    def finish(mut self):
        """Writes what is left (also after an error, so the lines made so far are out)."""
        self._write_pending()
