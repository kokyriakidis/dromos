# dromos

**Pangenome graph alignments, read as fast as the disk delivers them.**

dromos reads [vg](https://github.com/vgteam/vg)'s graph alignments, GAM or GAF, and turns them into what downstream
tools need, in one pass over the file, on every core. It is written in [Mojo](https://www.modular.com/mojo) and is
the indexing engine behind [pgphase](https://github.com/kokyriakidis/pgphase) (`pgphase index`).

| Command | What it writes |
|---|---|
| `dromos index` | the alignments placed on the reference, sorted, bgzipped and tabix-indexed: the index pgphase reads |
| `dromos gam-to-gaf` | a GAM as the GAF lines `vg convert -G` writes, in the GAM's order |

## Speed

HG002 PacBio HiFi on the HPRC v2.1 Minigraph-Cactus graph (a 155 GB GAM, 12.8 million alignments), on one 10-core
workstation (Intel i9-7900X, NVMe):

| From the GAM to pgphase's index | Wall time |
|---|---|
| `vg convert -G`, then `pgcat index-gaf` (before dromos) | about 2 h 50 min |
| `dromos gam-to-gaf --slim \| pgcat index-gaf` | 412 s |
| **`dromos index`** | **196 s**, in about 16 GB of memory |

Inflating vg's DEFLATE is half the work. It runs at about 870 MB/s of output per core with
[libdeflate](https://github.com/ebiggers/libdeflate), as fast as DEFLATE decoders go. A 10 GB GAF indexes in under
5 s, about the speed of the disk.

## Install

dromos builds with [pixi](https://pixi.sh), which fetches the pinned Mojo toolchain and libraries:

```sh
git clone https://github.com/kokyriakidis/dromos && cd dromos
pixi run build          # build/dromos
pixi run test           # unit tests
```

It is built and tested on Linux x86-64.

## Usage

Each command reads a file that pgcatalog (pgphase's graph tool, `tools/pgcatalog` in
[pgphase](https://github.com/kokyriakidis/pgphase)) writes once per graph, so dromos never loads the graph itself:

```sh
pgcatalog ref-map --graph graph.gbz --reference-sample CHM13 --out graph.refmap   # for index
pgcatalog node-table --graph graph.gbz --out graph.nodes                         # for gam-to-gaf
```

### `dromos index`

```sh
dromos index graph.refmap reads.gam reads.coord.gaf.gz       # writes reads.coord.gaf.gz and its .tbi
dromos index graph.refmap reads.gaf.gz reads.coord.gaf.gz    # a GAF: plain, gzip or bgzip
vg giraffe ... -o gaf | dromos index graph.refmap - reads.coord.gaf.gz    # straight from the aligner
```

- **Any input.** A GAM, or a GAF that is plain, gzip or bgzip, told apart by its first bytes; `--gam` or `--gaf`
  says which. `-` reads stdin. BGZF is read in parallel. Gzip that is not BGZF is a single stream and is read on
  one thread, so bgzip a large GAF.
- **What it writes.** Each alignment is placed on the reference by the nodes of its walk that lie on the reference
  sample's paths, and written as `contig  start  end  <its 12 GAF columns>` (1-based, inclusive; the tags `rq`, `pu`
  and `rs` are kept), sorted, bgzipped and tabix-indexed. Header lines `#max_span <contig> <bases>` give the longest
  alignment on each contig. Unmapped alignments, and those with no unique place on the reference, are left out and
  counted by reason.
- **Options.** `-t THREADS` (all cores); `-l LEVEL`, the BGZF level, 0 to 12. The default, 1, writes about three
  times as fast as 6, into a file 40% larger.

### `dromos gam-to-gaf`

```sh
dromos gam-to-gaf graph.nodes reads.gam > reads.gaf     # the lines vg convert -G writes
dromos gam-to-gaf --slim graph.nodes reads.gam          # their 12 columns alone
```

`dromos COMMAND -h` lists each command's options. `DROMOS_TIMING=1` reports where the time goes.

## How it works

- **Read once, straight from the disk.** The input is read with O_DIRECT, several requests at a time, 64 MiB per
  chunk. Through the page cache, a file larger than memory reads at a fifth of the NVMe's speed.
- **One pipeline, every core busy.** Each step reads the chunk after next, inflates the next chunk's BGZF blocks
  (checking every CRC), and hands slices of the current chunk to the command, all at once.
- **A GAM decoded only as far as needed.** A small protobuf reader frames GAM groups, resuming inside a group across
  chunks. For the index it sums each mapping's edits as it reads them. GAF lines follow libvgio's
  `alignment_to_gaf` line for line.
- **No graph in memory.** The reference map and node table are memory-mapped, so their pages load as they are used.
- **Sorted without spilling.** Placed alignments gather by megabase of the reference, compressed with zstd as they
  arrive. At the end, windows of buckets are sorted, laid out, deflated and written in a three-deep pipeline, with
  the tabix index built as each block's address becomes known.

## Correctness

dromos is checked byte for byte against the tools it replaces, on HG002 HiFi alignments to the HPRC graph:

- `gam-to-gaf` writes the same lines as `vg convert -G` (vg 1.70.0).
- `index` writes, at `-l 0`, the same file and `.tbi` as pgphase's previous indexer (`pgcat index-gaf --slim`),
  from the GAM and from its GAF in every encoding. At other levels it writes the same decompressed content, with the
  same answer to every tabix query. The whole-genome index matches the previous one record for record.
- `tests/integration/edges.sh` covers odd and damaged inputs: empty files, multi-member gzip, lines across chunk
  boundaries, any thread count, stdin, and truncated or corrupt files, which fail instead of leaving a partial
  index.

```sh
pixi run check ENV_FILE [--full]        # end-to-end checks; ENV_FILE names the data (see tests/integration/check.sh)
tests/integration/edges.sh ENV_FILE     # odd and damaged inputs
```

## License

MIT; see [`LICENSE`](LICENSE). At run time dromos loads [libdeflate](https://github.com/ebiggers/libdeflate) (MIT),
[zstd](https://github.com/facebook/zstd) (BSD-3-Clause) and [zlib](https://zlib.net) (zlib licence).
