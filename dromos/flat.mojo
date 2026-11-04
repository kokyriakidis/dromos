"""`Flat`: a List read and written without bounds checks, for hot loops over arrays whose indices are known good
(Mojo's default build checks every List index, which doubles the snarl finder's time)."""


struct Flat[T: TrivialRegisterPassable](Movable, Sized):
    var data: List[Self.T]

    def __init__(out self):
        self.data = List[Self.T]()

    def __init__(out self, *, length: Int, fill: Self.T):
        self.data = List[Self.T](length=length, fill=fill)

    def __init__(out self, *, unsafe_uninit_length: Int):
        self.data = List[Self.T](unsafe_uninit_length=unsafe_uninit_length)

    def __init__(out self, *, capacity: Int):
        self.data = List[Self.T](capacity=capacity)

    def __init__(out self, var data: List[Self.T]):
        self.data = data^

    @always_inline
    def __getitem__(self, i: Int) -> Self.T:
        return self.data.unsafe_ptr()[unsafe_offset=i]

    @always_inline
    def __setitem__(mut self, i: Int, value: Self.T):
        self.data.unsafe_ptr()[unsafe_offset=i] = value

    @always_inline
    def __len__(self) -> Int:
        return len(self.data)

    @always_inline
    def append(mut self, value: Self.T):
        self.data.append(value)

    @always_inline
    def pop(mut self) -> Self.T:
        var n = len(self.data) - 1
        var v = self.data.unsafe_ptr()[unsafe_offset=n]
        self.data.shrink(n)
        return v

    def clear(mut self):
        self.data.clear()

    def resize(mut self, *, unsafe_uninit_length: Int):
        self.data.resize(unsafe_uninit_length=unsafe_uninit_length)

    def copy(self) -> Self:
        return Self(self.data.copy())

    def unsafe_ptr(self) -> ImmPointer[Self.T, origin_of(self.data)]:
        return self.data.unsafe_ptr()

    def unsafe_mut_ptr(mut self) -> MutPointer[Self.T, origin_of(self.data)]:
        return self.data.unsafe_ptr()
