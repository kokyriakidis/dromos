# dromos

**Pangenome graph alignments, read as fast as the disk delivers them.**

dromos reads [vg](https://github.com/vgteam/vg)'s graph alignments, GAM or GAF, and turns them into what downstream
tools need, in one pass over the file, on every core. It reads the GBZ pangenome graph itself, through the
[gbz](https://github.com/jltsiren/gbwt-rs) crate, makes the rest of pgphase's inputs from it (the site catalog),
and keeps what the other commands need of it in one cache file beside it. It is written in [Mojo](https://www.modular.com/mojo) and is the indexing engine behind
[pgphase](https://github.com/kokyriakidis/pgphase) (`pgphase index`).

| Command | What it writes |
|---|---|
| `dromos index` | the alignments placed on the reference, sorted, bgzipped and tabix-indexed: the index pgphase reads |
| `dromos gam-to-gaf` | a GAM as the GAF lines `vg convert -G` writes, in the GAM's order |
| `dromos surject` | the alignments on every reference contig in one pass, the records `pgcatalog surject` writes: a sorted, indexed BAM, or SAM |
| `dromos catalog` | a GBZ graph's site catalog: every snarl on the reference, every contig, in one file of binary columns |
| `dromos view` | a catalog as VCF, a contig or a region of it: the text `pgcatalog view` prints |
| `dromos catalog-info` | a catalog's header: the graph it belongs to, its contigs, where its snarls came from |
| `dromos prepare` | a GBZ graph's cache, made ahead of the commands that read it |
| `dromos reference-samples` | the reference samples a GBZ graph declares |

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

`dromos surject` against `pgcatalog surject` (which loads the 6 GB graph), on the chr10 alignments of a 5.8 GB
slice of that GAM (487,571 reads; SAM to /dev/null, 10 threads), and against dromos before its aligners ran on
threads of their own and its band was filled in vectors (both dromos columns with turbo off):

| chr10 of the slice | `pgcatalog surject` (from its GAF) | dromos before | **`dromos surject`** (from the GAM) |
|---|---|---|---|
| default | 200 s, 11.5 GB | 78.4 s, 561 CPU s, 4.7 GB | **22.5 s**, 205 CPU s, 3.6 GB |
| `--graph-placement` | 46 s, 10.5 GB | 21.4 s, 186 CPU s, 3.6 GB | **15.2 s**, 134 CPU s, 2.3 GB |
| `--catalog` (chr10's windows) | 329 s, 14.1 GB | 178.2 s, 1,561 CPU s, 5.9 GB | **88.3 s**, 818 CPU s, 5.9 GB |

Every contig in one pass, into a sorted, indexed BAM (`-o`, level 4; output and spills on a hard disk that
writes 181 MB/s; 10 threads, turbo off):

| Every contig, one pass | Wall | CPU | Peak memory |
|---|---|---|---|
| the slice (489,489 reads; a 2.6 GB BAM) | 46.9 s | 396 s | 17.1 GB |
| the slice, `--graph-placement` | 36.2 s | 327 s | 17.2 GB |
| the slice, `--catalog` (every contig's windows) | 124 s ¹ | 1,160 s ¹ | 21.5 GB |
| **HG002, the whole GAM** (155 GB, 12,800,930 reads; a 67.9 GB BAM and its .bai) | **25 min 26 s** | 11,358 s | 27.7 GB |

The whole GAM takes 13.8 minutes to read and align (every core busy) and 11.3 to sort and write: 43 GB of records
spilled in 10 rounds (`-m 8G`) are read back while the 68 GB BAM is written, both on the one hard disk, which bounds
that phase. Before, a contig was a pass over the whole GAM, then `samtools sort`: chr20 alone takes 204 s with the
previous dromos (140 s now). In the slice's pass the BAM's deflate takes 27% of the CPU, the band's fill 14%,
inflating the GAM 8%, the records' zstd staging 7% and rebuilding the reads 5%.

¹ Measured while each catalog was inflated whole (31.7 GB at the peak); read a buffer at a time, as now, they cost
about 4 s more on the small data, and the peak is 21.5 GB.

`dromos catalog` from the same cached snarls, against `pgcatalog build` (which loads the 6 GB graph whole first)
and the previous dromos (both write pgcatalog's 25 per-contig `.pgsites` files), 10 threads, back to back on a quiet
machine, writing to a hard disk:

| HPRC v2.1, CHM13, every contig | Wall | CPU | Peak memory | Catalog |
|---|---|---|---|---|
| `pgcatalog build` (measured earlier) | 499 s | 2,588 s | 17.8 GB | 25 files, 1.51 GB |
| dromos before | 260.6 s | 2,162 s | 24.8 GB | 25 files, 1.51 GB |
| **`dromos catalog`** | **184.7 s** | 1,732 s | 19.2 GB | **one file, 487 MB** |
| **`dromos catalog`**, snarls from the graph's `.dist` | **191.0 s** | 1,737 s | 27.5 GB ² | the same records |
| `dromos catalog` with no snarls (`vg chunk` and `vg snarls` first) | VGPATH | | | |

The contigs are built three at a time, the biggest first, on one pool of workers, busy 98% of the build. On a busy
machine (two of the ten cores taken by other work) the same build takes 239 s, the previous one 322 s. The catalog's
view of every contig is pgcatalog view's text of pgcatalog's catalogs, and the catalog from the distance index prints
the same text as the one from vg snarls.

² The distance index (8.5 GB) is mapped and read as needed; its pages count in the peak.

A region query (`dromos view CATALOG chr20:START-END`, a process each) opens the catalog and reads its contig's
index in 10 ms; 100 kb of chr6 takes 24 ms against `pgcatalog view`'s 41 ms. Regions whose records make tens of MB
of VCF take up to twice as long as pgcatalog, which stores that text: dromos rebuilds it from the columns. In process,
the records of a region with their walks and alleles are read 1.0 to 1.5 times as fast as from `.pgsites` blocks
inflated and parsed for the same (chr6, 1 Mb: 18.9 ms against 28.4 ms; chr20, 100 kb: 35.4 ms against 36.0 ms).

## Install

dromos builds with [pixi](https://pixi.sh), which fetches the pinned Mojo toolchain, the Rust toolchain and the
libraries:

```sh
git clone https://github.com/kokyriakidis/dromos && cd dromos
pixi run build          # build/dromos, and beside it build/libdromos_gbz.so (the gbz crate)
pixi run test           # unit tests
```

dromos finds `libdromos_gbz` beside its binary, then at `$DROMOS_GBZ_LIB`; keep the two together.

It is built and tested on Linux x86-64.

## Usage

### The graph and its cache

`index`, `gam-to-gaf` and `surject` take the GBZ graph (`--graph`) and read what they need of it from its cache:
one file, `graph.gbz.dromos`, made the first time a command needs it and kept beside the graph. It holds the node
lengths and letters (4 bits a letter) and, for every reference sample the graph declares (CHM13 and GRCh38 in HPRC
graphs), where each node lies on the sample's paths. Commands map it into memory and read only the parts they use.

```sh
dromos prepare --graph graph.gbz                       # made now, for a graph folder others will only read
dromos index --graph graph.gbz --reference-sample CHM13 reads.gam reads.coord.gaf.gz
dromos gam-to-gaf --graph graph.gbz --cache /scratch/graph.dromos reads.gam > reads.gaf
```

- **Kept current.** The cache records the graph it was made from (its size, time and a hash of its first and last
  4 MB). A graph that has changed, a cache cut short or damaged, or one of another dromos version is made again in
  its place. A file there that is not a dromos cache, or another graph's cache, is refused, not written over.
- **Where.** Beside the graph, or wherever `--cache FILE` says (on every command that reads it). It is written under
  a temporary name in the same folder and renamed, so two runs making it at once both finish with a whole one. When
  the folder cannot be written, the command stops and says so: give `--cache`, or have the cache made there once
  with `dromos prepare`.
- **Reference samples.** `--reference-sample` names the reference, needed when the graph declares more than one
  (`dromos reference-samples graph.gbz` lists them). For a graph that declares none, the sample given is added to
  the cache.

On the HPRC graph the cache takes 8 s and 6.0 GB (both reference samples), in 9 GB of memory; the files it replaces,
a reference map per sample and a node table, took 7 s each and 3.3 GB, 3.2 GB and 6.4 GB.

### `dromos catalog`, `dromos view` and `dromos catalog-info`

```sh
dromos catalog --graph graph.gbz --reference-sample CHM13              # graph.CHM13.sites, every reference contig
dromos catalog --graph graph.gbz --reference-sample CHM13 --region chr20 --region chr6:28,000,000-34,000,000
dromos view graph.CHM13.sites chr20:10,000,000-12,000,000 | bgzip > sites.vcf.gz
dromos catalog-info graph.CHM13.sites                                  # the graph it belongs to, its contigs
```

- **What it writes.** One file per graph and reference sample (`GRAPH.SAMPLE.sites` beside the graph, or `--out`
  a file or folder): every snarl the reference path passes through, at every nesting level, with all its alleles,
  as `vg deconstruct -a` reports them without what needs haplotype names, for every reference contig, or for the
  regions asked for (`--region chr20:START-END`, 1-based, inclusive: the records a query of the whole contig
  returns). Records are kept in blocks by position, each block's columns (positions, flags, boundary visits, walks
  as indexes into the block's node table, LV, PS, the node table's ids, lengths and 2-bit letters) a zstd frame
  each, so a reader decodes only what it needs; the letters are in the file, so it needs nothing else to be read.
  [docs/catalog-format.md](docs/catalog-format.md) specifies the format (version 1).
- **The graph it belongs to.** The header names the graph (its fingerprint, node ids, GBZ and GBWT tags), the
  reference sample, where the snarls came from and each contig's length and MD5 (SAM's M5 of its reference
  sequence). `dromos surject --catalog` refuses a catalog of another graph; any reader refuses another format
  version, a truncated file or a damaged one (every block has a CRC).
- **Snarls.** From the graph's distance index (giraffe's `.dist`, `GRAPH.dist` beside the graph, or `--dist`), read
  in place and checked against the graph (its node ids and count, and each contig's chain against the graph's
  edges); its snarls are vg snarls', orientation included. Without one, from vg (`vg chunk -C` for each contig's
  component, then `vg snarls`), kept in `--snarls-cache` (`GRAPH.snarls` beside the catalog), which records the graph
  it belongs to; later runs read them from there. A cache made by `pgcatalog build` serves as well.
- **Options.** `--graph`, `--out`, `--reference-sample`, `--region` (repeatable), `-t`/`--threads`, `--level` (zstd,
  9), `--dist`, `--snarls-cache`, `--vg`, `--snarls-jobs` (vg snarls runs at a time, 6).
- **`dromos view CATALOG [CONTIG[:START-END]]`** prints, for a contig and a region of it, the text `pgcatalog view`
  printed of that contig's `.pgsites` catalog: its VCF header, then the records overlapping the region as tabix
  would return them. Without a contig, every contig, in one VCF.
- **`dromos catalog-info CATALOG`** prints the header, a field a line, then a line per contig.

### `dromos index`

```sh
dromos index --graph graph.gbz --reference-sample CHM13 reads.gam reads.coord.gaf.gz    # and its .tbi
dromos index --graph graph.gbz --reference-sample CHM13 reads.gaf.gz reads.coord.gaf.gz # a GAF: plain, gzip or bgzip
vg giraffe ... -o gaf | dromos index --graph graph.gbz --reference-sample CHM13 - reads.coord.gaf.gz
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
  times as fast as 6, into a file 40% larger. `--reference-sample NAME`, the reference; `--cache FILE`.

### `dromos gam-to-gaf`

```sh
dromos gam-to-gaf --graph graph.gbz reads.gam > reads.gaf     # the lines vg convert -G writes
dromos gam-to-gaf --slim --graph graph.gbz reads.gam          # their 12 columns alone
```

### `dromos surject`

```sh
dromos surject --graph graph.gbz --reference-sample CHM13 reads.gam -o reads.bam   # every contig: BAM, .bai
dromos surject --graph graph.gbz --reference-sample CHM13 --contig chr20 --contig chr6 reads.gaf.gz -o two.bam
dromos surject --graph graph.gbz --reference-sample CHM13 --contig chr20 reads.gam > chr20.sam   # SAM, IN's order
```

- **What it writes.** The alignments of IN on the reference contigs, the records `pgcatalog surject` writes for the
  same reads and options (a port of its surject.rs), every contig of the reference sample in one pass over IN
  (`--contig`, repeated, chooses some). Each read is rebuilt from its walk and `cs`; a base on a node a contig visits
  once takes that node's position, and the placed stretches are joined by aligning what lies between them to the
  reference (affine gaps, in a band). Read ends off the reference are aligned or soft-clipped, and ends diverging
  from the read's own graph alignment clipped. A read with bases on several contigs gets a record on each, the one
  that contig's own run makes; reads with no base on any are left out.
- **BAM.** `-o OUT.bam` writes a BAM sorted as `samtools sort` sorts (by contig, position, then reverse strand,
  then input order) and its index, `OUT.bam.bai`, or `.csi` when a contig is longer than a .bai holds (2^29 bases)
  or with `--csi`. The header is `@HD SO:coordinate`, an `@SQ` for each contig in the order of the reference's
  paths, and dromos's `@PG`. Records wait for the sort compressed by zstd, and past half of `-m` (8G) they spill to
  a file in `--tmp-dir` (OUT.bam's folder), unlinked as soon as it is made.
- **Level.** `-l` sets the BAM's deflate level, 4 by default, chosen by measurement on the slice below (11.6 GB of
  BAM records): libdeflate levels 1, 4 and 6 compress them to 2.91, 2.63 and 2.59 GB at 187, 101 and 52 MB/s a
  core, and the final sort and write takes 10.1, 15.0 and 19.4 s on 10 cores. 4 is within 2% of 6's size at twice
  its speed; on a hard disk (180 MB/s) the whole genome's final write waits on the disk at 4 anyway.
- **SAM.** Without `-o`, SAM on stdout in IN's order, a read's records in contig order; one contig's SAM is what it
  always was.
- **Names.** `--sam-names full` (`CHM13#0#chr20`, the default) or `short` (`chr20`); `--sam-contig NAME` names a
  single contig.
- **Any input.** A GAM, or a GAF with `cs` (and `bq` for qualities) in any encoding, told apart by its first bytes.
  A GAM's alignments are first written as vg's GAF lines, so a GAM and its GAF give the same records. dromos index's
  output keeps no `cs` and cannot be used.
- **Options.** `--graph-placement`, read bases placed only where the graph puts them; `--catalog FILE`, the graph's
  site catalog (`dromos catalog`; it must hold every contig surjected, and be of the same graph): reads spanning one
  of a contig's snarl windows aligned there as alleles; `--reference-sample NAME`, the
  reference; `--cache FILE`; `-t THREADS`.

`dromos COMMAND -h` lists each command's options. `DROMOS_TIMING=1` reports where the time goes.

## How it works

- **Read once, straight from the disk.** The input is read with O_DIRECT, several requests at a time, 64 MiB per
  chunk. Through the page cache, a file larger than memory reads at a fifth of the NVMe's speed.
- **One pipeline, every core busy.** Each step reads the chunk after next, inflates the next chunk's BGZF blocks
  (checking every CRC), and hands slices of the current chunk to the command, all at once.
- **A GAM decoded only as far as needed.** A small protobuf reader frames GAM groups, resuming inside a group across
  chunks. For the index it sums each mapping's edits as it reads them. GAF lines follow libvgio's
  `alignment_to_gaf` line for line.
- **No graph in memory.** The graph's cache is memory-mapped, so its pages load as they are used, and a command
  touches only its own parts: `index` the lengths and one sample's placements, `gam-to-gaf` the lengths and
  letters.
- **The graph read in parts, once.** The gbz crate reads the GBZ; dromos does the rest in Mojo. For the cache it has
  the crate load the GBWT and the node letters' pieces at the same time; then the reference samples' paths are
  walked on threads of their own while the nodes are found, measured and their letters unpacked and packed by range
  on every core, and each path's placements are laid out by a task of its own. Every section is written through a
  mapping of the file, by the threads that make it (buffered writes to one file take its lock in turn).
- **Surjection without barriers.** A read's alignment to the reference costs from nothing to a second. Each step's
  slices gather the reads that touch a contig, as GAF lines, into pieces of a ring; aligner threads that live for
  the whole run take the pieces as they are published, and each step writes out the prefix of pieces done since the
  last, in order. No thread waits for the slowest read of a batch; reading waits once 512 MB of lines are in flight.
- **The band in vectors.** The affine-gap band is filled a row at a time, 16 cells a vector (AVX-512): the
  insertion and diagonal scores come from the row before, and the deletion scores along the row from a running
  minimum, which gives exactly the integers, and so the tracebacks, of the cell-by-cell fill.
- **Sorted with spills.** surject's BAM records gather by megabase of a contig in zstd frames as they are made;
  past half the memory allowed they are written to a spill file by a thread of their own, a bucket's frames
  together. At the end, windows of buckets (a gigabyte of records at most) are read back, sorted, deflated and
  written in the index's three-deep pipeline, the .bai built as each block's address becomes known.
- **Catalogs on every core.** The graph's GBWT is loaded by the crate and node letters unpacked a page of ids at a
  time, as needed. The contigs are built three at a time, the biggest first; every parallel step of each (snarls
  traversed, records built, blocks encoded) is a job of one pool of workers, which runs the job of the contig
  started first, so a contig's tail overlaps the next contigs' work. Snarls are built in pgcatalog's batches, each
  step building one batch while the blocks sealed so far are encoded; workers with nothing left build the next
  batches. Snarls spanning megabases are searched ahead, their haplotypes split into parts; a record of gigabytes of
  text is encoded as a job of its own, its long columns by zstd's own threads with long-distance matching (chr2's
  catalog is 38 MB with it, 54 MB without, for the same CPU). The finished contigs are written in
  their order whatever order they finish in.
- **Sorted without spilling.** Placed alignments gather by megabase of the reference, compressed with zstd as they
  arrive. At the end, windows of buckets are sorted, laid out, deflated and written in a three-deep pipeline, with
  the tabix index built as each block's address becomes known.

## Correctness

dromos is checked byte for byte against the tools it replaces, on HG002 HiFi alignments to the HPRC graph:

- `gam-to-gaf` writes the same lines as `vg convert -G` (vg 1.70.0).
- `index` writes, at `-l 0`, the same file and `.tbi` as pgphase's previous indexer (`pgcat index-gaf --slim`),
  from the GAM and from its GAF in every encoding. At other levels it writes the same decompressed content, with the
  same answer to every tabix query. The whole-genome index matches the previous one record for record.
- `surject` writes the same SAM records as `pgcatalog surject` (pgcatalog 0.2.0), from the GAM and from its GAF,
  with and without `--graph-placement` and `--catalog`: chr10 of the slice (487,571 records), and slices of the whole
  GAM on chr20 (48,706) and chr6 across the MHC (47,027).
- `surject -o` writes, for every contig, the records its single-contig run gives piped through `samtools sort`
  (samtools 1.24, `samtools view OUT.bam CONTIG`), in every mode, on the slice, the chr20 and chr6 slices and the
  small data; a single contig's BAM holds samtools sort's records byte for byte, and its .bai (or .csi) says what
  `samtools index` says of the same file; spilling, threads and stdin change nothing. The whole HG002 GAM's BAM
  (12,733,701 records) indexes as `samtools index` indexes it, and its chr20 is chr20's single run over the whole GAM
  sorted by samtools.
- `reference-samples` prints the lines `pgcatalog reference-samples` prints.
- `catalog` writes, for every contig, the records `pgcatalog build` writes from the same snarls: `dromos view` of
  the one file prints pgcatalog view's text of each `.pgsites` catalog, byte for byte, for all 25 CHM13 contigs of
  the HPRC graph (43,494,146 records, 152 GB of VCF text, the largest record a 6.6 GB line on chr2), regions of
  them, and the small graph's catalogs pgphase's tests commit; the file is the same whatever the threads. Snarls
  read from the graph's distance index give the same file's records as vg snarls (all 51,186,896 snarls, vg's
  orientation included).
- `tests/integration/edges.sh` covers odd and damaged inputs: empty files, multi-member gzip, lines across chunk
  boundaries, any thread count, stdin, and truncated or corrupt files, which fail instead of leaving a partial
  index; and the graph's cache made, reused, made again (a changed graph, a cache cut short or damaged, another
  version's), kept elsewhere, made by several runs at once, and refused (a read-only folder, a file that is not a
  cache, another graph's cache); site catalogs of odd snarls files, regions and thread counts (the same file on any),
  and catalogs and distance indexes refused (truncated, damaged, another format version, another graph's).

```sh
pixi run check ENV_FILE [--full]        # end-to-end checks; ENV_FILE names the data (see tests/integration/check.sh)
tests/integration/edges.sh ENV_FILE     # odd and damaged inputs
```

## License

MIT; see [`LICENSE`](LICENSE). At run time dromos loads [libdeflate](https://github.com/ebiggers/libdeflate) (MIT),
[zstd](https://github.com/facebook/zstd) (BSD-3-Clause) and [zlib](https://zlib.net) (zlib licence).
`libdromos_gbz` builds in the [gbz](https://github.com/jltsiren/gbwt-rs) and
[simple-sds](https://github.com/jltsiren/simple-sds) crates (both MIT) and, through them, zstd.
