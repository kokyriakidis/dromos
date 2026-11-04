"""Checks that every snarl dromos computes (`dromos snarls -o MINE.tsv`) is one of BubbleFinder's (all snarls, on the
GFA segment graph). BubbleFinder (GPL-3.0) is used only as a program whose output is compared.

    python3 -I bfcheck.py TRANSLATION MINE.tsv BUBBLEFINDER.txt [--examples N]

TRANSLATION is `vg gbwt -Z G.gbz --translation`'s table (T segment node,node,...); empty when the GBZ has none
(BubbleFinder then names nodes by their ids). A snarl's boundary is mapped to
a segment incidence when it is a segment's end: the inward side of a start `>a` is a's right side (`a+` if a is the
segment's last node), of `<a` its left side (`a-` if a is the first); of an end `>b` b's left side, of `<b` its right
side. A BubbleFinder line is a set of incidences any two of which are a snarl. Snarls with a boundary inside a segment cannot be mapped (counted).
"""

import sys


def main():
    trans_path, mine_path, bf_path = sys.argv[1:4]
    examples = 10
    if "--examples" in sys.argv:
        examples = int(sys.argv[sys.argv.index("--examples") + 1])
    first_of, last_of = {}, {}
    for line in open(trans_path):
        p = line.rstrip("\n").split("\t")
        if p[0] != "T":
            continue
        nodes = [int(x) for x in p[2].split(",")]
        first_of[nodes[0]] = p[1]
        last_of[nodes[-1]] = p[1]
    bf = set()
    with open(bf_path) as f:
        next(f)
        for line in f:
            p = line.split()
            # A line is a set of incidences any two of which are a snarl (BubbleFinder's compact output).
            for i in range(len(p)):
                for j in range(i + 1, len(p)):
                    bf.add(frozenset((p[i], p[j])))

    identity = not first_of

    def side(node, right):
        if identity:
            return str(node) + ("+" if right else "-")
        if right:
            s = last_of.get(node)
            return None if s is None else s + "+"
        s = first_of.get(node)
        return None if s is None else s + "-"

    total = found = unmapped = 0
    missing = []
    for line in open(mine_path):
        p = line.split("\t")
        a, b = int(p[1]), int(p[2])
        total += 1
        x = side(a >> 1, (a & 1) == 0)
        y = side(b >> 1, (b & 1) == 1)
        if x is None or y is None:
            unmapped += 1
            continue
        if frozenset((x, y)) in bf:
            found += 1
        else:
            missing.append((a, b, x, y))
    print(f"mine {total}; in BubbleFinder's {found}; not in it {len(missing)}; boundary inside a segment {unmapped};"
          f" BubbleFinder snarls {len(bf)}")
    for a, b, x, y in missing[:examples]:
        print(f"  not in BubbleFinder: {'<' if a & 1 else '>'}{a >> 1}{'<' if b & 1 else '>'}{b >> 1} = {x} {y}")


main()
