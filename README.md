# dromos

Pangenome graph alignments in Mojo. `dromos gam-to-gaf` reads vg's GAM, block-parallel, and writes the GAF lines
`vg convert -G` writes for it, in the input's order.

```sh
pgcat node-table --graph graph.gbz --out graph.nodes     # once per graph: node lengths and letters
pixi run build
pixi run ./build/dromos gam-to-gaf -t 16 graph.nodes reads.gam > reads.gaf
```

## How

- **BGZF, block-parallel.** The GAM is read 64 MiB at a time; each chunk's blocks inflate in parallel with libdeflate,
  every block's CRC32 checked (after sporos).
- **GAM groups framed, messages decoded on every thread.** A group is a varint count, the `GAM` tag and its messages;
  each `vg.Alignment` is decoded by a small protobuf reader and written as GAF by a port of libvgio's
  `alignment_to_gaf`, each thread its slice of the chunk, the slices written in order.
- **No graph in memory.** GAF needs each node's length (and, for cs, its letters): a flat node table, written once per
  graph by `pgcat node-table`, stands in for the GBZ.

## Checked

The same lines as `vg convert -G` on HG002 HiFi alignments to the HPRC v2.1 Minigraph-Cactus graph (vg 1.70.0).
