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

## Working here

- Build: `pixi run build` (build/dromos); test: `pixi run test`; format: `pixi run format` (120 columns).
- End to end: `pixi run check ENV_FILE [--full]` (tests/integration/check.sh) checks dromos's outputs against vg's and
  pgcat's byte for byte; ENV_FILE names the data on the machine that has it. Run it after any change that
  could touch output.
- Toolchain: one Mojo nightly, pinned exactly in `pixi.toml` (1.2.0.dev2026102005), the one karan, sporos and plegma
  build with. Move it on purpose: change the pin, `pixi lock`, build and test.
- Build and test where the data are (a GAM, the graph's ref map and node table): `pixi run ...`.
- `gaf_line` follows libvgio's `alignment_to_gaf` (src/alignment_io.cpp) line for line; keep it so, and check any
  change against `vg convert -G` (the same lines; vg's order varies, dromos keeps the input's).
- `dromos index` must write what `pgcat index-gaf --slim --ref-map` writes (pgphase's tools/pgcat: gaf.rs, bgzf.rs,
  tabix.rs; refmap.mojo and tabix.mojo port them): at `-l 0` the same file byte for byte, and the same .tbi
  decompressed (the check does this). At other levels only the deflate differs: compare decompressed content and
  tabix queries.
- The pipeline (stream.mojo) is shared: a sink says what messages (a GAF's lines) become (convert.mojo, index.mojo).
  GAF input must index as pgcat indexes the same file, in every encoding (the check covers plain, gzip, BGZF, stdin
  and edge cases).
- Pointers handed to parallel tasks carry their owner's origin, or the owner is an explicit allocation that outlives
  the tasks (see sporos's CLAUDE.md for the patterns, copied here from sporos: parallel.mojo, output.mojo, bgzf.mojo).
- Finished work merges straight into `main` and is pushed; no pull requests.
