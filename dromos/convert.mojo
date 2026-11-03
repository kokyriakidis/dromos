"""`dromos gam-to-gaf`: a GAM's alignments as the GAF lines `vg convert -G` writes for them, in the input's order."""

from dromos.gam import Alignment, gaf_line
from dromos.nodes import NodeTable
from dromos.output import OutputBuffer
from dromos.stream import MessageSink, mut_untracked


struct GafSink[tags: Bool](MessageSink, Movable):
    """Each slice writes its messages' lines into a text of its own; between steps the texts go out in order."""

    var nodes: NodeTable
    var texts: List[List[UInt8]]
    var alignments: List[Alignment]
    var out: OutputBuffer

    def __init__(out self, var nodes: NodeTable, threads: Int):
        self.nodes = nodes^
        self.texts = List[List[UInt8]](length=4 * threads, fill=List[UInt8]())
        self.alignments = List[Alignment]()
        for _ in range(4 * threads):
            self.alignments.append(Alignment())
        self.out = OutputBuffer(1)
        self.out.write("@HD\tVN:Z:1.0\n")

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
        text[].clear()
        for m in range(first, last):
            a[].decode(data, messages[unsafe_offset=2 * m], messages[unsafe_offset=2 * m + 1])
            gaf_line[Self.tags](data, a[], self.nodes, text[])

    def side_tasks(self) -> Int:
        return 0

    def side_task(mut self, index: Int) raises:
        pass

    def after_step(mut self, slices: Int) raises:
        for t in range(slices):
            self.out.write_bytes(Span(self.texts[t]))

    def finish(mut self):
        self.out.flush()
