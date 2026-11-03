"""The MD5 of a BAM's records, header left out: an uncompressed BAM stream on stdin (bgzip -dc FILE.bam), so two
BAMs with different headers (@PG lines) compare by their record bytes, every field as encoded.

    bgzip -dc A.bam | python3 -I bam_records.py
"""

import hashlib
import struct
import sys


def main():
    f = sys.stdin.buffer
    if f.read(4) != b"BAM\x01":
        raise SystemExit("not a BAM")
    (l_text,) = struct.unpack("<i", f.read(4))
    f.read(l_text)
    (n_ref,) = struct.unpack("<i", f.read(4))
    for _ in range(n_ref):
        (l_name,) = struct.unpack("<i", f.read(4))
        f.read(l_name + 4)
    h = hashlib.md5()
    while True:
        block = f.read(1 << 22)
        if not block:
            break
        h.update(block)
    print(h.hexdigest())


if __name__ == "__main__":
    main()
