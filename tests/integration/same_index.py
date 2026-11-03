"""Whether two BAM indexes (.bai, or .csi, BGZF) say the same: every reference's bins with their chunks (and a CSI's
linear offsets), its linear index, and the unplaced count. Bins may come in any order (htslib writes them in hash
order, dromos in numeric order); everything else must be equal.

    python3 -I same_index.py A.bai B.bai     exit 0 if the same, 1 (with the first difference) if not
"""

import gzip
import struct
import sys


def read(path):
    with open(path, "rb") as f:
        data = f.read()
    if data[:2] == b"\x1f\x8b":
        data = gzip.decompress(data)
    at = 0

    def take(fmt):
        nonlocal at
        v = struct.unpack_from("<" + fmt, data, at)
        at += struct.calcsize(fmt)
        return v

    magic = data[:4]
    at = 4
    csi = magic == b"CSI\x01"
    if not csi and magic != b"BAI\x01":
        raise SystemExit(f"{path}: not a .bai or .csi")
    head = ()
    if csi:
        min_shift, depth, l_aux = take("iiI")
        at += l_aux
        head = (min_shift, depth)
    (n_ref,) = take("i")
    refs = []
    for _ in range(n_ref):
        (n_bin,) = take("i")
        bins = {}
        for _ in range(n_bin):
            if csi:
                b, loff, n_chunk = take("IQi")
            else:
                b, n_chunk = take("Ii")
                loff = None
            chunks = [take("QQ") for _ in range(n_chunk)]
            bins[b] = (loff, chunks)
        linear = []
        if not csi:
            (n_intv,) = take("i")
            linear = [take("Q")[0] for _ in range(n_intv)]
        refs.append((bins, linear))
    no_coor = take("Q")[0] if at + 8 <= len(data) else None
    return magic, head, refs, no_coor


def main():
    a, b = read(sys.argv[1]), read(sys.argv[2])
    if a[0] != b[0] or a[1] != b[1]:
        print(f"different kinds or depths: {a[0]} {a[1]} / {b[0]} {b[1]}")
        return 1
    if len(a[2]) != len(b[2]):
        print(f"{len(a[2])} references against {len(b[2])}")
        return 1
    for i, (x, y) in enumerate(zip(a[2], b[2])):
        if x[0] != y[0]:
            only = sorted(set(x[0]) ^ set(y[0]))
            diff = [k for k in x[0] if k in y[0] and x[0][k] != y[0][k]]
            print(f"reference {i}: bins differ (only in one: {only[:5]}, different: {diff[:5]})")
            return 1
        if x[1] != y[1]:
            print(f"reference {i}: linear index differs")
            return 1
    if a[3] != b[3]:
        print(f"unplaced counts differ: {a[3]} / {b[3]}")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
