"""Compares dromos's computed snarls with vg snarls' (`dromos snarls --vg-file SNARLS -o VG.tsv` and
`dromos snarls --graph G.gbz -o MINE.tsv`).

    python3 -I compare.py VG.tsv MINE.tsv [--examples N] [--all]

(--all: compare every component of MINE.tsv, not only those vg's snarls touch.)

Checks, for the snarls of the components vg's file touches:
- the same set of (start, end) visits, oriented exactly;
- the same parent, compared as a snarl in either orientation (vg writes the parent's unregularized orientation),
  and how often vg's parent orientation equals the parent's own written orientation;
- chains: each of mine (its non-trivial snarls, in order) is a run of consecutive siblings in vg's file, in the same
  order, and siblings with one parent split into exactly my chains.
"""

import sys
from collections import defaultdict


def canon(s, e):
    rs, re_ = e ^ 1, s ^ 1
    return (rs, re_) if (rs, re_) < (s, e) else (s, e)


def name(v):
    return ("<" if v & 1 else ">") + str(v >> 1)


def sname(k):
    return name(k[0]) + name(k[1])


def load(path, mine):
    rows = []
    with open(path) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            s, e = int(p[1]), int(p[2])
            parent = None if p[3] == "-" else (int(p[3]), int(p[4]))
            row = {"comp": int(p[0]), "key": (s, e), "parent": parent}
            if mine:
                row["chain"] = int(p[5])
                row["rank"] = int(p[6])
                row["kind"] = int(p[8]) if len(p) > 8 else 0
                row["tied"] = len(p) > 9 and p[9] == "1"
            rows.append(row)
    return rows


def main():
    vg_path, mine_path = sys.argv[1], sys.argv[2]
    examples = 10
    if "--examples" in sys.argv:
        examples = int(sys.argv[sys.argv.index("--examples") + 1])
    vg = load(vg_path, False)
    mine_all = load(mine_path, True)
    # The components vg's snarls touch: those of my snarls sharing a canonical key with one of vg's.
    by_canon_mine = {canon(*r["key"]): r for r in mine_all}
    comps = set()
    for r in vg:
        m = by_canon_mine.get(canon(*r["key"]))
        if m is not None:
            comps.add(m["comp"])
    if "--all" in sys.argv:
        comps = {r["comp"] for r in mine_all}
    mine = [r for r in mine_all if r["comp"] in comps]
    print(f"vg {len(vg)} snarls; mine {len(mine)} in {len(comps)} component(s) ({len(mine_all)} in all)")
    vg_keys = {r["key"]: r for r in vg}
    mine_keys = {r["key"]: r for r in mine}
    if len(vg_keys) != len(vg):
        print(f"  vg has duplicate snarls: {len(vg) - len(vg_keys)}")
    only_vg = [k for k in vg_keys if k not in mine_keys]
    only_mine = [k for k in mine_keys if k not in vg_keys]
    vg_canon = {canon(*k) for k in vg_keys}
    reversed_ = [k for k in only_mine if canon(*k) in vg_canon]
    tied = [k for k in reversed_ if mine_keys[k]["tied"]]
    print(f"  snarls: same {len(vg_keys) - len(only_vg)}, only vg {len(only_vg)}, only mine {len(only_mine)}"
          f" (of which reversed {len(reversed_)}, in chains with a tied vote {len(tied)})")
    for k in only_vg[:examples]:
        print(f"    only vg: {sname(k)} parent {sname(vg_keys[k]['parent']) if vg_keys[k]['parent'] else '-'}")
    for k in only_mine[:examples]:
        r = mine_keys[k]
        print(f"    only mine: {sname(k)} parent {sname(r['parent']) if r['parent'] else '-'} chain {r['chain']}"
              f" rank {r['rank']}")
    # Parents.
    same_parent = differ = oriented_same = 0
    bad = []
    for k, r in vg_keys.items():
        m = mine_keys.get(k)
        if m is None:
            continue
        a, b = r["parent"], m["parent"]
        if (a is None) != (b is None) or (a is not None and canon(*a) != canon(*b)):
            differ += 1
            if len(bad) < examples:
                bad.append((k, a, b))
            continue
        same_parent += 1
        if a is not None and a == b:
            oriented_same += 1
    with_parent = sum(1 for r in vg if r["parent"] is not None)
    print(f"  parents: same {same_parent}, differ {differ}; vg's parent as written (regularized) for"
          f" {oriented_same} of {with_parent} with one")
    for k, a, b in bad:
        print(f"    {sname(k)}: vg parent {sname(a) if a else '-'}, mine {sname(b) if b else '-'}")
    # Chains: vg's siblings in file order.
    siblings = defaultdict(list)
    for r in vg:
        p = r["parent"]
        siblings[canon(*p) if p else None].append(r["key"])
    position = {}
    for p, keys in siblings.items():
        for i, k in enumerate(keys):
            position[k] = (p, i)
    chains = defaultdict(list)
    circular = set()
    for r in mine:
        c = r["chain"] + 1000000000 * r["comp"]
        chains[c].append((r["rank"], r["key"]))
        if r["kind"] == 3:
            circular.add(c)
    good = broken = rotated = 0
    examples_b = []
    covered = defaultdict(int)
    for c, members in chains.items():
        members.sort()
        keys = [k for _, k in members]
        pos = [position.get(k) for k in keys]
        if c in circular and all(x is not None for x in pos) and len(keys) > 1:
            # A circular chain (the root cycle): vg starts it where its DFS closed the cycle; compare cyclically.
            first = min(range(len(pos)), key=lambda i: pos[i][1])
            if first != 0:
                rotated += 1
                keys = keys[first:] + keys[:first]
                pos = pos[first:] + pos[:first]
        ok = all(x is not None for x in pos)
        if ok:
            ok = all(x[0] == pos[0][0] for x in pos) and all(pos[i + 1][1] == pos[i][1] + 1 for i in range(len(pos) - 1))
        if ok:
            good += 1
            covered[pos[0][0]] += len(keys)
        else:
            broken += 1
            if len(examples_b) < examples:
                examples_b.append((c, keys, pos))
    print(f"  chains: {good} are runs of vg siblings in order, {broken} are not;"
          f" {rotated} circular chains start elsewhere in vg")
    for c, keys, pos in examples_b:
        print(f"    chain {c}: {[sname(k) for k in keys[:6]]} at {pos[:6]}")
    unsplit = sum(1 for p, keys in siblings.items() if covered[p] != len(keys))
    print(f"  sibling groups not covered exactly by my chains: {unsplit} of {len(siblings)}")
    exact = not only_vg and not only_mine and differ == 0 and broken == 0 and unsplit == 0
    # Reversed snarls of chains whose vote tied (and chains broken only by them) are vg's hash-order choice.
    ties_only = len(only_vg) == len(only_mine) == len(reversed_) == len(tied) and differ == 0
    print("EXACT" if exact else ("TIES ONLY" if ties_only else "DIFFERENT"))


main()
