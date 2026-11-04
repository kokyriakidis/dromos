"""Random graphs: dromos's computed snarls (`dromos snarls`) against vg snarls'.

    python3 -I fuzz.py WORKDIR SEED COUNT [--vg VG] [--dromos BIN] [--small]

Each graph is the union of random haplotype walks (GFA W lines) over a backbone: substitutions, insertions,
deletions, inversions (a stretch walked in reverse), duplications (a stretch walked twice: cycles), walks starting or
ending inside the graph (tips), self loops and reversing edges, sometimes a second component. vg gbwt makes the GBZ,
vg snarls (-t 1) the expected snarls. Node lengths are random so that ties in path and cycle lengths (which vg
breaks by hash order) are rare. Prints a line per graph that differs and a summary.
"""

import os
import random
import subprocess
import sys


def make_graph(rng, small):
    n_backbone = rng.randint(3, 12 if small else 60)
    nodes = {}
    next_id = [1]

    def new_node():
        i = next_id[0]
        next_id[0] += 1
        nodes[i] = "".join(rng.choice("ACGT") for _ in range(rng.randint(1, 40)))
        return i

    backbone = [new_node() for _ in range(n_backbone)]
    walks = [[(v, False) for v in backbone]]
    alt_pool = []
    for _ in range(rng.randint(1, 6 if small else 12)):
        w = []
        i = 0
        while i < len(backbone):
            r = rng.random()
            if r < 0.08:
                # substitution: an alternative node (shared among walks sometimes)
                if alt_pool and rng.random() < 0.5:
                    w.append((rng.choice(alt_pool), False))
                else:
                    a = new_node()
                    alt_pool.append(a)
                    w.append((a, False))
                i += 1
            elif r < 0.13:
                i += rng.randint(1, 3)  # deletion
            elif r < 0.17:
                a = new_node()
                w.append((a, False))  # insertion
                w.append((backbone[i], False))
                i += 1
            elif r < 0.21 and i + 2 < len(backbone):
                k = rng.randint(1, min(4, len(backbone) - i))
                for v in reversed(backbone[i:i + k]):
                    w.append((v, True))  # inversion
                i += k
            elif r < 0.24 and i + 2 < len(backbone):
                k = rng.randint(1, min(4, len(backbone) - i))
                seg = [(v, False) for v in backbone[i:i + k]]
                w.extend(seg)
                w.extend(seg)  # duplication: a cycle
                i += k
            elif r < 0.25:
                w.append((backbone[i], False))
                w.append((backbone[i], False))  # self loop
                i += 1
            elif r < 0.26:
                w.append((backbone[i], False))
                w.append((backbone[i], True))  # reversing edge
                i += 1
            else:
                w.append((backbone[i], False))
                i += 1
        if not w:
            continue
        # tips: start or end inside
        if rng.random() < 0.2:
            w = w[rng.randint(0, len(w) - 1):]
        if rng.random() < 0.2:
            w = w[:rng.randint(1, len(w))]
        if rng.random() < 0.1:
            t = new_node()
            w.append((t, False))  # a tip node
        walks.append(w)
    if rng.random() < 0.15:
        # a second component
        other = [new_node() for _ in range(rng.randint(2, 8))]
        walks.append([(v, False) for v in other])
        if rng.random() < 0.5:
            walks.append([(other[0], False), (new_node(), False), (other[-1], False)])
    if rng.random() < 0.1:
        # a circular walk: the backbone closing on itself
        walks.append([(v, False) for v in backbone] + [(backbone[0], False)])
    return nodes, walks


def write_gfa(path, nodes, walks):
    edges = set()
    for w in walks:
        for (a, ra), (b, rb) in zip(w, w[1:]):
            e = (a, ra, b, rb)
            f = (b, not rb, a, not ra)
            edges.add(min(e, f))
    used = {v for w in walks for v, _ in w}
    with open(path, "w") as f:
        f.write("H\tVN:Z:1.1\n")
        for i in sorted(used):
            f.write(f"S\t{i}\t{nodes[i]}\n")
        for a, ra, b, rb in sorted(edges):
            f.write(f"L\t{a}\t{'-' if ra else '+'}\t{b}\t{'-' if rb else '+'}\t0M\n")
        for k, w in enumerate(walks):
            text = "".join(("<" if r else ">") + str(v) for v, r in w)
            f.write(f"W\tsample{k}\t0\tchr\t0\t{len(w)}\t{text}\n")


def run(cmd, **kw):
    return subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, **kw)


def main():
    work, seed, count = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
    vg = "vg"
    dromos = "dromos"
    if "--vg" in sys.argv:
        vg = sys.argv[sys.argv.index("--vg") + 1]
    if "--dromos" in sys.argv:
        dromos = sys.argv[sys.argv.index("--dromos") + 1]
    small = "--small" in sys.argv
    here = os.path.dirname(os.path.abspath(__file__))
    os.makedirs(work, exist_ok=True)
    rng = random.Random(seed)
    exact = differ = failed = tied_votes = tied_root = 0
    for k in range(count):
        nodes, walks = make_graph(rng, small)
        base = os.path.join(work, f"g{seed}_{k}")
        write_gfa(base + ".gfa", nodes, walks)
        r = run([vg, "gbwt", "-G", base + ".gfa", "--gbz-format", "-g", base + ".gbz"])
        if r.returncode != 0:
            failed += 1
            print(f"{base}: vg gbwt failed: {r.stderr.decode()[-300:]}")
            continue
        with open(base + ".vg.snarls", "wb") as f:
            r = subprocess.run([vg, "snarls", "-t", "1", base + ".gbz"], stdout=f, stderr=subprocess.PIPE)
        if r.returncode != 0:
            failed += 1
            print(f"{base}: vg snarls failed: {r.stderr.decode()[-300:]}")
            continue
        r1 = run([dromos, "snarls", "--vg-file", base + ".vg.snarls", "-o", base + ".vg.tsv", "-t", "1"])
        r2 = run([dromos, "snarls", "--graph", base + ".gbz", "-o", base + ".mine.tsv", "--stats", "-t", "1"])
        ties = sum(int(x) for line in r2.stderr.decode().splitlines() if "; ties cycles " in line
                   for x in line.split("; ties cycles ")[1].replace("leaves", "").replace("paths", "").split())
        if r1.returncode != 0 or r2.returncode != 0:
            failed += 1
            print(f"{base}: dromos snarls failed: {(r1.stderr + r2.stderr).decode()[-500:]}")
            continue
        c = run([sys.executable, "-I", os.path.join(here, "compare.py"), base + ".vg.tsv", base + ".mine.tsv", "--all"])
        text = c.stdout.decode()
        if text.strip().endswith("TIES ONLY"):
            tied_votes += 1
            print(f"{base}: TIES ONLY")
        elif text.strip().endswith("EXACT"):
            exact += 1
            for ext in (".gfa", ".gbz", ".vg.snarls", ".vg.tsv", ".mine.tsv"):
                os.remove(base + ext)
        elif ties > 0:
            # Broken some other way (each tie met, first or last of the equals: DROMOS_SNARL_TIES=mask:N), do the ties
            # give vg's snarls?
            other = None
            for mask in range(1, min(1 << (ties + 2), 1024)):
                r3 = run([dromos, "snarls", "--graph", base + ".gbz", "-o", base + ".tie.tsv", "-t", "1"],
                         env=dict(os.environ, DROMOS_SNARL_TIES=f"mask:{mask}"))
                c3 = run([sys.executable, "-I", os.path.join(here, "compare.py"), base + ".vg.tsv", base + ".tie.tsv",
                          "--all"])
                if r3.returncode == 0 and c3.stdout.decode().strip().splitlines()[-1] in ("EXACT", "TIES ONLY"):
                    other = mask
                    break
            if other is not None:
                tied_root += 1
                print(f"{base}: DIFFERENT, with {ties} length ties; vg's with ties broken as mask {other}")
            else:
                differ += 1
                print(f"{base}: DIFFERENT, with {ties} length ties, however broken\n{text}")
        else:
            differ += 1
            print(f"{base}: DIFFERENT\n{text}")
    print(f"seed {seed}: {count} graphs, exact {exact}, orientation ties only {tied_votes}, different with length ties"
          f" {tied_root}, different {differ}, failed {failed}")


main()
