# dromos

Pangenome graph alignments in Mojo: vg's GAM read as fast as the disk delivers it, written as pgphase's index (the
file `pgcat index-gaf --slim` writes) or as the GAF lines `vg convert -G` writes (github.com/kokyriakidis/pgphase).

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

- `cli/dromos_cli.mojo`: the commands, `index` and `gam-to-gaf`, each with its own options and `-h`.
- `dromos/stream.mojo`: the pipeline every command runs on: reads (O_DIRECT), decoding (BGZF, plain, gzip), GAM
  told from GAF by the first bytes (`looks_like_gaf`), framing (GAM messages or GAF lines), and a `MessageSink`
  taking slices of each chunk.
- `dromos/gam.mojo`: GAM framing and protobuf decoding; `gaf_line` (vg's lines) and `slim_line` (the 12 columns).
- `dromos/convert.mojo`: the `gam-to-gaf` sink. `dromos/index.mojo`: the `index` sink (placement, counts).
- `dromos/refmap.mojo`: pgcat's reference map and placement rule. `dromos/gaf_text.mojo`: a GAF line as pgcat reads
  it. `dromos/staging.mojo`: placed lines by bucket. `dromos/sort_write.mojo`: the final sort, BGZF and tabix
  (`dromos/tabix.mojo`).
- `dromos/bgzf.mojo`: files, BGZF blocks, libdeflate, zstd and zlib. `dromos/nodes.mojo`: the node table.
  `dromos/output.mojo`, `dromos/parallel.mojo`: writing, integers, threads.

## Working here

- Build: `pixi run build` (build/dromos); test: `pixi run test`; format: `pixi run format` (120 columns).
- End to end: `pixi run check ENV_FILE [--full]` (tests/integration/check.sh) checks dromos's outputs against vg's and
  pgcat's byte for byte; ENV_FILE names the data on the machine that has it. Run it after any change that
  could touch output. tests/integration/edges.sh ENV_FILE runs odd and damaged inputs (empty files, gzip and BGZF
  shapes, the chunk boundary, threads, stdin, truncation): every line must say ok.
- Toolchain: one Mojo nightly, pinned exactly in `pixi.toml` (1.2.0.dev2026102005), the one karan, sporos and plegma
  build with. Move it on purpose: change the pin, `pixi lock`, build and test.
- Build and test where the data are (a GAM, the graph's ref map and node table): `pixi run ...`.
- `gaf_line` follows libvgio's `alignment_to_gaf` (src/alignment_io.cpp) line for line; keep it so, and check any
  change against `vg convert -G` (the same lines; vg's order varies, dromos keeps the input's).
- `dromos index` must write what pgphase's previous indexer, `pgcat index-gaf --slim --ref-map`, writes (pgcat, now
  pgcatalog, dropped index-gaf for dromos; the checks run a pgcat built from pgphase's rescue-experiment branch,
  ff2b4b2; refmap.mojo and tabix.mojo port its gaf.rs and tabix.rs): at `-l 0` the same file byte for byte, and the
  same .tbi decompressed (the check does this). At other levels only the deflate differs: compare decompressed content
  and tabix queries.
- The pipeline (stream.mojo) is shared: a sink says what messages (a GAF's lines) become (convert.mojo, index.mojo).
  GAF input must index as pgcat indexes the same file, in every encoding (the check covers plain, gzip, BGZF, stdin
  and edge cases).
- Pointers handed to parallel tasks carry their owner's origin, or the owner is an explicit allocation that outlives
  the tasks (see sporos's CLAUDE.md for the patterns, copied here from sporos: parallel.mojo, output.mojo, bgzf.mojo).
- Finished work merges straight into `main` and is pushed; no pull requests.
