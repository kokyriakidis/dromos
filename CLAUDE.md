# dromos

Pangenome graph alignments in Mojo: vg's GAM read as fast as the disk delivers it, written as pgphase's index (the
file `pgcat index-gaf --slim` writes), as the GAF lines `vg convert -G` writes, or surjected onto the reference
contigs (`pgcatalog surject`'s records, every contig in one pass, as SAM or a sorted, indexed BAM)
(github.com/kokyriakidis/pgphase). It reads the GBZ graph itself, through the gbz crate: it makes the graph's site
catalog, one file of every contig's snarls in binary columns (`catalog`; `view` prints pgcatalog view's text of it,
`catalog-info` its header), and keeps what the other commands read of the graph in one cache file beside it
(`G.gbz.dromos`, `prepare`). dromos alone turns a graph and its alignments into pgphase's inputs.

dromos is one repository holding a graph-algorithms library (`dromos/`) and its command line (`cli/`): a
dependency-light Mojo alternative to vg and libbdsg for what pgphase needs (the GBZ read through the gbz crate, the
snarl decomposition, haplotype walks, reference positions, linearization).

## Layering

- Library modules (`dromos/*.mojo`: graph access, formats, the snarl tree, haplotype walks, reference positions,
  linearization) never import `cli/`. They raise errors; they do not parse arguments or print usage.
- `cli/dromos_cli.mojo` is a thin layer: it parses each command's options, prints its usage and calls the command's
  entry point in the library (`run_build`, `run_index`, `run_gam_to_gaf`, `run_surject`, `graph_cache`, `view`,
  `info`, `find_all_snarls`, `run_dist_snarls`). No algorithm lives there.
- Each algorithm has an importable entry point and its own tests (tests/test_dromos.mojo, and the integration
  checks), and is accepted only once checked against its reference (vg, libbdsg, pgcatalog) on the HPRC graph.

## Mojo skills

Modular's Mojo skills live in `.claude/skills` and correct pretrained Mojo, which is out of date. Hooks in
`.claude/settings.json` inject them: the full `mojo-syntax` skill and `references/idiomatic-mojo.md` at session
start (also after `/clear` and compaction) and into every subagent, and a reminder on every prompt
(`DROMOS_SKILLS_EVERY_PROMPT=full` injects the full text on every prompt instead). Follow them over pretrained
habits:

- `mojo-syntax` before writing, editing or reviewing any Mojo, with `references/idiomatic-mojo.md`.
- Include the skills rule in every subagent prompt.

`pixi run sync-skills` refreshes them from upstream; `.claude/skills/UPSTREAM` records the commit.

## Layout

- `cli/dromos_cli.mojo`: the commands, `index`, `gam-to-gaf`, `surject`, `catalog`, `view`, `catalog-info`,
  `snarls`, `prepare` and `reference-samples`, each with its own options and `-h`. index, gam-to-gaf and surject take
  the graph as `--graph G.gbz` (and `--cache FILE`) and read its cache.
- `gbz-ffi/`: the gbz crate (=0.6.2, with simple-sds =0.4.2) behind a C interface, a Rust cdylib
  (`libdromos_gbz`, built by `pixi run build` beside build/dromos; dromos finds it there, then at $DROMOS_GBZ_LIB).
  A thin layer: the whole graph (`dgbz_load`), header-only reads (tags, metadata), and the parts the graph's cache and
  catalogs need (`dgbz_parts_*`: the GBWT and the node letters' pieces, read at once, letters unpacked by range; edges
  and the bidirectional search for many visits or states a call, as the crate's GBZ gives them; every node's edges
  and length by ranges of ids, `dgbz_parts_node_edges`, and the ids' bounds, `dgbz_parts_id_bounds`).
  `ABI` in lib.rs and graph.mojo must agree (4); raise both on any change to the interface.
- `dromos/graph.mojo`: the interface in Mojo (`GbzGraph`, `GbzParts`, `Metadata`, `SearchState`, `NodeSource`), and
  pgcatalog's graph.rs ported (reference samples and paths, `resolve_reference_sample`, `fingerprint`, `graph_id`).
- `dromos/cache.mojo`: the graph's cache, its format (documented there), `CacheWriter`, `CacheFile` (mapped, nothing
  read ahead but the sections a reader asks for) and `examine` (current, missing, not a cache, another version,
  damaged). `dromos/prepare.mojo`: the cache made from a GBZ (`make_cache`: one load of the graph's parts; paths
  walked in the background while the nodes are found and read; letters packed in the background while the
  placements are laid out; every section written through a mapping) and found or made again (`graph_cache`).
- `dromos/snarls.mojo`: vg snarls files (plain, gzip, BGZF), pgcatalog's snarls.rs. `dromos/dist.mojo`: snarls from a
  giraffe distance index (libbdsg's SnarlDistanceIndex, versions 3 and 4), read in place, oriented as vg snarls writes
  them; a dependency-light reader and a second oracle for the snarl finder, not a source of the catalog's snarls.
  `dromos/distgraph.mojo`: `dromos snarls --dist` (`run_dist_snarls`): every component's snarls of the index
  (`dist_snarls`), and the checks that refuse another graph's index (`check_graph`: node ids and count, each reference
  path in one component, top-level chains' steps the graph's edges). `dromos/bigraph.mojo`: the graph's nodes and edges
  read from the GBZ by ranges of ids on every thread and split into weakly connected components (`Components`: a CSR
  `Bigraph` each, and each node's component). `dromos/snarlfind.mojo`: vg snarls' snarls computed (vg's
  IntegratedSnarlFinder and SnarlManager rules on flat arrays; its notes say which step is which), the components found
  in parallel, the biggest first (`find_all_snarls`, `snarls_of_nodes`, `write_snarls_tsv`). `dromos/flat.mojo`: `Flat`,
  a List read and written without bounds checks, for the finder's hot loops. `dromos/catalog.mojo`: one contig's records
  (pgcatalog's catalog.rs, record.rs and reference.rs): node letters a page at a time, the reference path and its
  traversals, the haplotype walk search, records, their blocks and `encode_block` (a block's columns), and the steps
  that build and encode them on the pool. `dromos/build.mojo`: the `catalog` command (build.rs: regions, `SnarlSource`,
  the one interface to the snarls' two sources: computed from the graph's own load, the default, `compute_snarls`; or
  the snarls cache, `--snarls-cache`, vg chunk and vg snarls when a contig's snarls are not cached; contigs built three
  at a time, the biggest first by their computed snarls' count or their snarls file's size; `--dist` is refused).
  `dromos/sites.mojo`: the site catalog's format (docs/catalog-format.md is its specification): the writer, the reader,
  `view`, `info` and surject's snarl windows. `dromos/pool.mojo`: worker threads shared by several jobs at once (every
  contig's steps).
- `dromos/stream.mojo`: the pipeline every command runs on: reads (O_DIRECT), decoding (BGZF, plain, gzip), GAM
  told from GAF by the first bytes (`looks_like_gaf`), framing (GAM messages or GAF lines), and a `MessageSink`
  taking slices of each chunk.
- `dromos/gam.mojo`: GAM framing and protobuf decoding; `gaf_line` (vg's lines) and `slim_line` (the 12 columns).
- `dromos/convert.mojo`: the `gam-to-gaf` sink and command (`run_gam_to_gaf`). `dromos/index.mojo`: the `index` sink
  (placement, counts) and command (`run_index`).
- `dromos/refmap.mojo`: pgcat's reference map (a reference sample's sections of the cache) and placement rule.
  `dromos/gaf_text.mojo`: a GAF line as pgcat reads it. `dromos/staging.mojo`: placed lines by bucket.
  `dromos/sort_write.mojo`: the final sort, BGZF and tabix; `IndexWriter` cuts BGZF into fixed blocks and indexes
  records as their blocks get addresses, for any `BlockIndex`. `dromos/tabix.mojo`: htslib's binning index at any
  depth (`BinIndex`), as a .tbi (`TbiBuilder`).
- `dromos/surject.mojo`: the `surject` command (`SurjectOptions`, `run_surject`) and sink: reads rebuilt from GAF
  lines, a record on each contig a read places a base on, pieces published to a ring and aligned by aligner threads
  that run for the whole input, records written out in order (SAM on stdout, or BAM records to the staging buckets)
  as each prefix of pieces is done;
  `dromos/linearize.mojo`: its alignment to the reference (banded affine gaps, a row of the band in vectors,
  chaining, clipping); `dromos/contig.mojo`: a contig's letters and once-visited nodes from the cache's map and node
  table, and `Reference`, the contigs surjected to (made at once) and which of them a node lies on; with
  `load_windows` in pgsites.mojo they port pgcatalog's surject.rs, reference.rs and pileup.rs.
- `dromos/bam.mojo`: BAM records encoded as `samtools view -b` encodes the same SAM lines, the header, and the
  .bai/.csi `BamIndex`. `dromos/bam_sort.mojo`: `surject -o`'s sort: records staged by megabase bucket in zstd
  frames, spilled to an unlinked file past half the memory allowed, read back, sorted, deflated and indexed.
- `dromos/bgzf.mojo`: files, BGZF blocks, libdeflate, zstd and zlib. `dromos/nodes.mojo`: the node table (the
  cache's lengths and letters). `dromos/output.mojo`, `dromos/parallel.mojo`: writing, integers, threads.

## Working here

- Build: `pixi run build` (build/dromos and build/libdromos_gbz; Rust comes from pixi, pinned in pixi.toml, and
  gbz-ffi/Cargo.lock is committed); test: `pixi run test`; format: `pixi run format` (120 columns).
- The gbz crate reads the GBZ: headers and tags, path metadata, GBWT path walks, node sequences, edges, the
  bidirectional search. Everything above it is Mojo. gbz-ffi holds no pgcatalog logic, and nothing reimplements the
  GBZ, GBWT or simple-sds formats in Mojo.
- The graph's cache is dromos's own format, versioned: `VERSION` in cache.mojo is 1. Any change to its layout or
  meaning raises it (a cache of another version is made again, never misread), and the format notes at the top of
  cache.mojo say what it holds. The cache is judged by what is made from it: the checks compare index, gam-to-gaf,
  surject, catalog and view against pgcat, vg and pgcatalog, not the cache itself. `reference-samples` prints
  pgcatalog's lines (pgcatalog from pgphase's plegma-msa branch, 60a3d09, `cargo build --release` in
  tools/pgcatalog).
- A cache is written under a temporary name in its folder and renamed; a file that is not a dromos cache, or
  another graph's (its name and fingerprint differ), is never written over; a damaged one, another version's, or one
  of an older state of the graph is made again. edges.sh checks each case.
- End to end: `pixi run check ENV_FILE [--full] [--only GROUP,...]` (tests/integration/check.sh) checks dromos's
  outputs against vg's and pgcat's byte for byte; ENV_FILE names the data on the machine that has it. Groups: cache,
  catalog, snarls, convert, index, surject; each prints its time. tests/integration/edges.sh ENV_FILE runs odd and
  damaged inputs (empty files, gzip and BGZF shapes, the chunk boundary, threads, stdin, truncation, the graph's
  cache): every line must say ok.
- The test loop must be fast: run the groups a change touches (`--only`), and all of them, with `--full`, only for a
  change to shared code (stream.mojo, the gbz FFI, the cache) or a release. Reference outputs (pgcat's, pgcatalog's,
  vg's snarls as tables) are made once in WORK and kept; never empty WORK to check. Cases run in parallel where they
  are independent (contigs), never repeat what another case implies, and use a byte comparison before any slow
  explaining tool. A slow check is a bug: time it and fix it.
- Toolchain: one Mojo nightly, pinned exactly in `pixi.toml` (1.2.0.dev2026102005), the one karan, sporos and plegma
  build with. Move it on purpose: change the pin, `pixi lock`, build and test.
- Build and test where the data are (a GAM and its GBZ graph): `pixi run ...`.
- `gaf_line` follows libvgio's `alignment_to_gaf` (src/alignment_io.cpp) line for line; keep it so, and check any
  change against `vg convert -G` (the same lines; vg's order varies, dromos keeps the input's).
- `dromos index` must write what pgphase's previous indexer, `pgcat index-gaf --slim --ref-map`, writes (pgcat, now
  pgcatalog, dropped index-gaf for dromos; the checks run a pgcat built from pgphase's rescue-experiment branch,
  ff2b4b2; refmap.mojo and tabix.mojo port its gaf.rs and tabix.rs): at `-l 0` the same file byte for byte, and the
  same .tbi decompressed (the check does this). At other levels only the deflate differs: compare decompressed content
  and tabix queries.
- The site catalog is dromos's own format, versioned: `VERSION` in sites.mojo is 1, and docs/catalog-format.md says
  exactly what a reader needs (pgphase's C++ reader is written from it). Any change to the layout or meaning raises
  VERSION and the document together; a reader refuses other versions, and truncated or damaged files, saying which.
- `dromos view` of a `dromos catalog` file must print, for each contig and region, pgcatalog view's text of pgcatalog
  build's catalog of it from the same snarls, byte for byte; the file must be the same for any threads. Check any change
  with `pixi run check` (its catalog and view cases; `--full` for the whole genome, against pgcatalog's catalogs, and
  the catalog from vg snarls against the computed snarls') and tests/integration/edges.sh. `DROMOS_CONTIGS_AT_ONCE`
  changes how many contigs a build runs at once (3; the file is the same for any). Haplotype walks may be found in any
  order (they are sorted); a heavy snarl's haplotypes are searched in parts whose walks and counts add up.
- Snarls computed by dromos (the catalog's default) and snarls from a distance index (`dromos snarls --dist`) must be vg
  snarls' (the same start and end visits, the same parents read in either direction, the same chains): the catalog from
  either of its two sources is the same file but for its header's `snarls.*` keys, on the whole HPRC genome, and `pixi
  run check ENV --full` compares the distance index's table (DIST in the env file names the graph's .dist) with vg
  snarls' files and with the computed snarls (compare_genome.sh with MINE or ORACLE: every contig EXACT); the small
  graph's tests/data/micb-kir3dl1.dist is compared always. vg writes a parent in its orientation before `regularize`,
  which follows its hash orders; the catalog reads a parent in either direction (`canonical`). Where vg breaks a tie by
  hash order (equal-length cycles or leaf paths, a tied orientation vote) the finder takes its own first;
  `DROMOS_SNARL_TIES=last|mask:N` breaks them otherwise, for tests/integration/fuzz.py. Check any change to the finder
  with the unit tests, `pixi run check ENV --full` (compare_genome.sh: every contig of the cached vg snarls EXACT),
  fuzz.py (random graphs against vg snarls: differences only from ties), and bfcheck.py (every snarl is one of
  BubbleFinder's). BubbleFinder is GPL-3.0: run its binary as an oracle, never read or copy its source.
- `dromos surject` must write the SAM records `pgcatalog surject` writes (pgphase's tools/pgcatalog), byte for byte,
  from a GAM and from its GAF: check any change with `pixi run check` (its surject cases; `--full` for the slice).
  Float ties in `chain` use libm's log2, as Rust does. Without `--contig` it surjects every contig of the reference
  sample in one pass; a read gets a record on each contig it places a base on, exactly the record that contig's
  single run makes (one contig's SAM stays byte for byte what it was).
- `dromos surject -o OUT.bam`: for every contig, `samtools view OUT.bam CONTIG` must give that contig's single run
  piped through `samtools sort` (which sorts by reference, position, reverse strand, then input order: checked
  against samtools 1.24 itself); a single contig's records are samtools sort's byte for byte; the .bai or .csi
  must say what `samtools index` (`-c`) says of the same file (tests/integration/same_index.py: bins may come in
  any order). check.sh and edges.sh hold these cases (SAMTOOLS names samtools in the env file).
- `Banded.fill` (vectors) must give `fill_reference`'s (cell by cell, pgcatalog's) traceback bytes and last row
  exactly, and `clip_lengths` its base-by-base results: the unit tests compare them on random inputs; keep both in
  step with any change to the costs or the band.
- The pipeline (stream.mojo) is shared: a sink says what messages (a GAF's lines) become (convert.mojo, index.mojo).
  GAF input must index as pgcat indexes the same file, in every encoding (the check covers plain, gzip, BGZF, stdin
  and edge cases).
- Pointers handed to parallel tasks carry their owner's origin, or the owner is an explicit allocation that outlives
  the tasks (see sporos's CLAUDE.md for the patterns, copied here from sporos: parallel.mojo, output.mojo, bgzf.mojo).
- A raw (untracked) pointer does not keep its owner alive: Mojo destroys a value after its last named use, on raising
  paths too. Read a mapped file through a function that borrows the mapping for the whole call (sites.mojo's `_read`,
  dist.mojo's `_header`), and give readers of a buffer its origin (`Cursor`).
- Unit tests write their scratch files in `$DROMOS_TEST_DIR` (else /tmp): on xaros, a folder on /mnt/sdb.
- Finished work merges straight into `main` and is pushed; no pull requests.
