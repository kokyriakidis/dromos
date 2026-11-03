# dromos

Pangenome graph alignments in Mojo: vg's GAM read as fast as the disk delivers it, written as pgphase's index (the
file `pgcat index-gaf --slim` writes), as the GAF lines `vg convert -G` writes, or surjected onto the reference
contigs (`pgcatalog surject`'s records, every contig in one pass, as SAM or a sorted, indexed BAM)
(github.com/kokyriakidis/pgphase). It reads the GBZ graph itself, through the gbz crate: it makes the site catalogs
pgphase's pgcatalog made (`catalog`, `view`), and keeps what the other commands read of the graph in one cache file
beside it (`G.gbz.dromos`, `prepare`). dromos alone turns a graph and its alignments into pgphase's inputs.

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

- `cli/dromos_cli.mojo`: the commands, `index`, `gam-to-gaf`, `surject`, `catalog`, `view`, `prepare` and
  `reference-samples`, each with its own options and `-h`. index, gam-to-gaf and surject take the graph as
  `--graph G.gbz` (and `--cache FILE`) and read its cache.
- `gbz-ffi/`: the gbz crate (=0.6.2, with simple-sds =0.4.2) behind a C interface, a Rust cdylib
  (`libdromos_gbz`, built by `pixi run build` beside build/dromos; dromos finds it there, then at $DROMOS_GBZ_LIB).
  A thin layer: the whole graph (`dgbz_load`), header-only reads (tags, metadata), and the parts the graph's cache and
  catalogs need (`dgbz_parts_*`: the GBWT and the node letters' pieces, read at once, letters unpacked by range; edges
  and the bidirectional search for many visits or states a call, as the crate's GBZ gives them).
  `ABI` in lib.rs and graph.mojo must agree; raise both on any change to the interface.
- `dromos/graph.mojo`: the interface in Mojo (`GbzGraph`, `GbzParts`, `Metadata`, `SearchState`, `NodeSource`), and
  pgcatalog's graph.rs ported (reference samples and paths, `resolve_reference_sample`, `fingerprint`, `graph_id`).
- `dromos/cache.mojo`: the graph's cache, its format (documented there), `CacheWriter`, `CacheFile` (mapped, nothing
  read ahead but the sections a reader asks for) and `examine` (current, missing, not a cache, another version,
  damaged). `dromos/prepare.mojo`: the cache made from a GBZ (`make_cache`: one load of the graph's parts; paths
  walked in the background while the nodes are found and read; letters packed in the background while the
  placements are laid out; every section written through a mapping) and found or made again (`graph_cache`).
- `dromos/snarls.mojo`: vg snarls files (plain, gzip, BGZF), pgcatalog's snarls.rs. `dromos/catalog.mojo`: one
  contig's catalog (pgcatalog's catalog.rs, record.rs and reference.rs): node letters a page at a time, the reference
  path and its traversals, the haplotype walk search, records and their lines, and the steps that build, write and
  compress them. `dromos/build.mojo`: the `catalog` command (build.rs: regions, the snarls cache, vg chunk and vg
  snarls when a contig's snarls are not cached). `dromos/pgsites.mojo`: the `.pgsites` format: the writer, the reader,
  `view`, and surject's snarl windows.
- `dromos/stream.mojo`: the pipeline every command runs on: reads (O_DIRECT), decoding (BGZF, plain, gzip), GAM
  told from GAF by the first bytes (`looks_like_gaf`), framing (GAM messages or GAF lines), and a `MessageSink`
  taking slices of each chunk.
- `dromos/gam.mojo`: GAM framing and protobuf decoding; `gaf_line` (vg's lines) and `slim_line` (the 12 columns).
- `dromos/convert.mojo`: the `gam-to-gaf` sink. `dromos/index.mojo`: the `index` sink (placement, counts).
- `dromos/refmap.mojo`: pgcat's reference map (a reference sample's sections of the cache) and placement rule.
  `dromos/gaf_text.mojo`: a GAF line as pgcat reads it. `dromos/staging.mojo`: placed lines by bucket.
  `dromos/sort_write.mojo`: the final sort, BGZF and tabix; `IndexWriter` cuts BGZF into fixed blocks and indexes
  records as their blocks get addresses, for any `BlockIndex`. `dromos/tabix.mojo`: htslib's binning index at any
  depth (`BinIndex`), as a .tbi (`TbiBuilder`).
- `dromos/surject.mojo`: the `surject` sink: reads rebuilt from GAF lines, a record on each contig a read places a
  base on, pieces published to a ring and aligned by aligner threads that run for the whole input, records written
  out in order (SAM on stdout, or BAM records to the staging buckets) as each prefix of pieces is done;
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
- End to end: `pixi run check ENV_FILE [--full]` (tests/integration/check.sh) checks dromos's outputs against vg's and
  pgcat's byte for byte; ENV_FILE names the data on the machine that has it. Run it after any change that
  could touch output. tests/integration/edges.sh ENV_FILE runs odd and damaged inputs (empty files, gzip and BGZF
  shapes, the chunk boundary, threads, stdin, truncation, the graph's cache): every line must say ok.
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
- `dromos catalog` must write pgcatalog build's `.pgsites` byte for byte (the same records, blocks, zstd frames and
  index; zstd 1.5.7, pgcatalog's, with its parameters per frame), and `dromos view` print pgcatalog view's text, for
  any threads: check any change with `pixi run check` (its catalog and view cases; `--full` for the whole genome)
  and tests/integration/edges.sh. Haplotype walks may be found in any order (they are sorted); a heavy snarl's
  haplotypes are searched in parts whose walks and counts add up.
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
- Finished work merges straight into `main` and is pushed; no pull requests.
