# dromos

Pangenome graph alignments in Mojo, read as fast as the disk delivers them.

- `dromos index` turns vg's GAM, or a GAF, into pgphase's index in one pass: each alignment placed on the reference,
  sorted, bgzipped and tabix-indexed. It writes the same file `pgcat index-gaf --slim` writes from the same
  alignments as GAF, with no pipe and no second process. A GAF may be plain, gzip or BGZF, or come on stdin (`-`).
- `dromos gam-to-gaf` writes the GAF lines `vg convert -G` writes for a GAM, in the input's order.

```sh
pgcat ref-map --graph graph.gbz --out graph.refmap         # once per graph: where each node lies on the reference
pgcat node-table --graph graph.gbz --out graph.nodes       # once per graph, for gam-to-gaf: node lengths and letters
pixi run build
pixi run ./build/dromos index -t 16 graph.refmap reads.gam reads.coord.gaf.gz      # and reads.coord.gaf.gz.tbi
pixi run ./build/dromos index -t 16 graph.refmap reads.gaf.gz reads.coord.gaf.gz   # a GAF: *.gaf[.gz|.bgz] or --gaf
vg giraffe ... -o gaf | pixi run ./build/dromos index --gaf graph.refmap - reads.coord.gaf.gz
pixi run ./build/dromos gam-to-gaf -t 16 graph.nodes reads.gam > reads.gaf
```

## How

- **BGZF, block-parallel.** The GAM is read 64 MiB at a time; each chunk's blocks inflate in parallel with libdeflate,
  every block's CRC32 checked (after sporos).
- **GAM groups framed, messages decoded on every thread.** A group is a varint count, the `GAM` tag and its messages;
  each `vg.Alignment` is decoded by a small protobuf reader and written as GAF by a port of libvgio's
  `alignment_to_gaf`, each thread its slice of the chunk, the slices written in order.
- **No graph in memory.** GAF needs each node's length (and, for cs, its letters): a flat node table, written once per
  graph by `pgcat node-table`, stands in for the GBZ. The index needs where each node lies on the reference, and the
  node lengths, which pgcat's reference map holds.
- **GAF too.** A GAF's lines are framed as a GAM's messages are, then classified and kept slim as pgcat does (its
  12 columns, and its rq/pu/rs tags). BGZF inflates in parallel and plain text is copied in parallel; gzip that is
  not BGZF is one stream, inflated by zlib on one thread of each step.
- **One pass to the index.** Each alignment's 12 GAF columns are laid out and placed on the reference as pgcat places
  them (a port of its rule), on every thread. Placed lines gather by 1 Mb bucket of the reference, compressed by
  zstd a megabyte at a time beside the GAM's own tasks. Then a window of buckets at a time is inflated and sorted,
  a bucket per thread, and its blocks deflated in parallel (libdeflate) and written in order, cut at exactly
  0xff00 bytes as pgcat cuts them, the tabix index built as they go (a port of htslib's rules, via pgcat).

## Checked

On HG002 HiFi alignments to the HPRC v2.1 Minigraph-Cactus graph (vg 1.70.0):

- `gam-to-gaf`: the same lines as `vg convert -G`.
- `index`: at level 0 (stored blocks) the same file as `pgcat index-gaf --slim --ref-map --level 0` byte for byte, and
  a .tbi of the same content, from the GAM and from its GAF (plain, gzip, BGZF, stdin; vg's full GAF with its tags;
  comments, CRLF, short lines, a last line without a newline); at other levels the same decompressed content and the same tabix queries. The default,
  level 1, deflates about 3 times as fast as level 6 into a file the size of pgcat's (zlib at 6); at `-l 6`
  libdeflate makes it about 30% smaller than pgcat's.
