# The dromos site catalog, format version 1

`dromos catalog` writes one file per graph and reference sample: every snarl the reference path passes through, at
every nesting level, with all its alleles, for every reference contig (or the regions asked for). `dromos view`
prints it as the VCF pgcatalog view printed of pgcatalog's per-contig `.pgsites` files, byte for byte; this document
says how, so a reader can be written from it alone (`dromos/sites.mojo` is the reference implementation, and keeps
this document in step).

A reader must refuse a file whose version is not one it implements. dromos reads version 1 only.

## Conventions

- Integers in the preamble, header and index are unsigned 64-bit little-endian (u64) unless noted.
- A string is a u64 byte length, then that many bytes of UTF-8, no terminator.
- Columns use LEB128 varints: seven bits a byte, least significant group first, the high bit set on every byte but
  the last (at most 10 bytes for a u64).
- zigzag(d) = (d << 1) XOR (d >> 63) (arithmetic shift) maps a signed difference to an unsigned one; its inverse is
  (z >> 1) XOR -(z AND 1).
- A visit is a node id and an orientation, coded `id << 1 | reverse`. In text a visit is `>ID` (forward) or `<ID`
  (reverse), the id in decimal.
- Positions are 1-based, as in VCF.
- zstd frames are standard frames (RFC 8878) with their content size recorded. A frame may use a window up to 2^31
  bytes (long-distance matching for long columns): decompress into a buffer of the known size in one call
  (`ZSTD_decompressDCtx`), or raise `ZSTD_d_windowLogMax` to 31 for streaming.

## Layout

```
offset 0                    preamble (64 bytes)
offset 64                   blocks, one after another
index_offset                each contig's index (one zstd frame each), in contig order
header_offset               header (one zstd frame), ending the file
```

### Preamble

| offset | size | field |
|---|---|---|
| 0 | 8 | magic `DRMSITES` (ASCII) |
| 8 | 4 | format version, u32 (1) |
| 12 | 4 | zero |
| 16 | 8 | the file's size in bytes |
| 24 | 8 | header offset |
| 32 | 8 | header frame size |
| 40 | 8 | index offset: where the contigs' index frames start |
| 48 | 8 | the index frames' total size |
| 56 | 8 | FNV-1a 64 of bytes 0 to 55 (offset basis 0xCBF29CE484222325, prime 0x100000001B3, a byte at a time) |

A reader checks, in order: the magic (a file starting `PGSITES` is pgcatalog's per-contig format, not this one); the
file has 64 bytes; the version; the hash; the recorded size against the file's (smaller: truncated); that the index
frames lie after the preamble and before the header, and the header inside the file. The header and index frames
carry zstd's content checksum, which decompression checks.

### Header

The header frame decompresses to:

```
u64        format version (equal to the preamble's)
u64        P, then P properties: string key, string value
u64        C, then C contigs:
             string   name           chr20
             string   CHROM          CHM13#0#chr20: the records' CHROM, <sample>#0#<name>
             u64      length         the contig's length (one past its last reference base)
             16 bytes MD5            of the contig's reference sequence (below)
             u64      blocks         how many blocks it has (B below)
             u64      records
             u64      index offset   where its index frame is
             u64      index size     the frame's size (inside the index region the preamble gives)
             u64      R, then R regions: u64 start, u64 end (1-based, inclusive; R = 0: the whole contig)
```

Contigs are in natural order (`chr1 < chr2 < … < chr10 < … < chrX`: runs of digits compared as numbers); a reader
looks a contig up by name (`chr20`, or `CHM13#0#chr20` when the sample matches) and reads that contig's index alone.

A contig with regions holds only the records overlapping at least one of them (the records a query of the whole
contig's catalog returns for that region).

The MD5 is of the contig's reference sequence as the graph spells it: the reference path's fragments laid at their
offsets, each node's letters in capitals (reverse-complemented when the path visits it backwards), and `N` where no
fragment lies, `length` letters in all. Where the path spells the whole contig it is the reference FASTA's SAM M5
(HPRC v2.1's chr20 has CHM13's f1a31f30d3da50893470ff47bf0b5a38; its chrY path stops 4.8 Mb short of CHM13's chrY, so
its MD5 is that of the 57.7 Mb the graph has).

Properties (all values are text; a reader should ignore keys it does not know):

| key | value |
|---|---|
| `dromos.version` | the dromos that wrote it |
| `graph.file` | the GBZ's file name |
| `graph.fingerprint` | `size=<bytes> fnv1a=<16 hex digits>`: the GBZ's size and the FNV-1a 64 of its first and last 4 MiB (dromos's `fingerprint`); dromos refuses to use the catalog with a graph whose fingerprint differs |
| `graph.nodes`, `graph.min_node`, `graph.max_node` | the graph's node count and id range |
| `gbz.<key>`, `gbwt.<key>` | each of the GBZ's tags and its GBWT's (`gbz.pggname`: the graph's content name) |
| `reference.sample` | the sample whose paths are the reference (CHM13) |
| `snarls.source` | `distance index` or `vg snarls` |
| `snarls.file`, `snarls.dist.version`, `snarls.dist.bytes` | the distance index: its file name, its SnarlDistanceIndex version, its size |
| `snarls.vg` | the vg that made the snarls (`vg version`'s first line), when known |
| `snarls.orientation` | how snarls are oriented (vg snarls' orientation either way) |
| `build.block_bp`, `build.block_bytes`, `build.own_block_bytes`, `build.level` | the block rules and zstd level |

### A contig's index

Each contig's index frame decompresses to exactly `16 + 8 * B * (9 + 2K)` bytes:

```
u64        B, the number of its blocks (equal to the header's)
u64        K, the number of columns (13)
9 arrays of B u64, in this order:
  first_pos    POS of the block's first record
  last_pos     POS of its last record
  max_end      the largest POS + max(len(REF), 1) - 1 of its records (the span tabix gives a VCF record)
  reach        the largest max_end of the contig's blocks up to this one (non-decreasing within a contig)
  offset       the block's offset in the file (its bytes lie before the index region)
  records      how many records it holds
  min_allele   the fewest bases any allele of its records spans (the padding base not counted)
  max_walk     the most visits any walk of its records has (the boundary visits counted)
  crc32        the CRC-32 (zlib's) of the block's bytes
B * K u64      the compressed size of each block's columns (block-major: block b's column c at b * K + c)
B * K u64      their raw sizes
```

A contig's blocks are in position order: first_pos and reach never decrease. The records overlapping [beg, end] are
in the blocks from the first whose `reach` ≥ beg to the last whose `first_pos` ≤ end (two binary searches), among
those whose own `max_end` ≥ beg; a reader keeps each such block's records with `POS ≤ end` and
`POS + max(len(REF), 1) - 1 ≥ beg`.

pgphase's first pass may also skip a block whose `min_allele` exceeds the largest allele span it handles, and size
its buffers by `max_walk`.

### Blocks

A block's bytes are its 13 columns, one after another in column order, each stored as:

- nothing, when the column is empty (compressed and raw size 0);
- its raw bytes, when the compressed size equals the raw size;
- otherwise one zstd frame of the compressed size, decompressing to the raw size.

Column c starts at `offset + Σ compressed sizes of the columns before it`. A reader decompresses only the columns it
needs: positions and IDs need columns 0 to 3; walks 4 to 6 and the node ids (9); alleles also the lengths and letters
(10 to 12), which hold most of a block's bytes.

Records within a block are in file order (POS, then the ID's text in byte order). Each column below is read in
record order; "per record" values appear once per record, "per walk" values once per walk of each record, and so on.

| # | name | contents |
|---|---|---|
| 0 | POS | per record, varint: POS minus the previous record's POS (the block's first record: minus first_pos, so 0) |
| 1 | FLAGS | per record, one byte: bit 0 `REVERSED`, bit 1 `PAD`, bit 2 `PS`; other bits 0 |
| 2 | PAD | per record with `PAD`, one byte: the padding base (A, C, G, T or N) |
| 3 | BOUND | per record, two varints: zigzag(first − previous record's first) (the block's first record: minus 0), then zigzag(last − first), where first and last are the first and last visits of the record's walks |
| 4 | NWALK | per record, varint: its walks (alleles), the reference's first; at least 2 |
| 5 | WLEN | per walk, varint: its inner visits (the visits between first and last) |
| 6 | WALK | per inner visit, varint: zigzag(i − previous i) << 1 \| reverse, where i is the visit's node's index in the block's node table and the previous i runs across walks and records (the block's first: minus 0) |
| 7 | LV | per record, varint: the LV field |
| 8 | PS | per record with `PS`, two varints: zigzag(ps_start − id_start), then zigzag(ps_end − ps_start), the PS field's two visits |
| 9 | NODEID | per table node, varint: its id minus the previous node's (the first: minus 0); ids strictly increase |
| 10 | NODELEN | per table node, varint: its letters |
| 11 | SEQ | the table's letters, nodes one after another, 2 bits each, A C G T as 0 1 2 3, four a byte from the low bits; ceil(letters / 4) bytes |
| 12 | EXC | per letter that is not A, C, G or T: varint, its index among the table's letters minus the previous exception's (the first: minus 0), then one byte, the letter (its 2 bits in SEQ are 0) |

Every walk of a record starts with the same visit (first) and ends with the same visit (last). The table holds
every node an inner visit of the block visits, once, by id. Its letters are the node's forward letters in capitals
(a lower-case letter in the graph is stored upper-cased; any other byte as it is, through EXC).

A reader checks that each column is consumed exactly and that table indexes are below the table's size; the CRC
covers the whole block.

### A record as VCF text

With the record's first and last visits, its flags, and its walks (walk k: first, its inner visits, last):

- `id_start, id_end` = (last XOR 1, first XOR 1) when `REVERSED` (the snarl was found from its end), else
  (first, last).
- The allele of walk k is the padding base (when `PAD`), then for each inner visit in order its node's letters:
  forward, as stored; reverse, reverse-complemented (A↔T, C↔G, any other letter to N).
- REF is walk 0's allele; ALT the others' joined by `,` (empty when there is one walk).
- len(REF), for overlap, is walk 0's allele's length.

The line, fields separated by tabs and ended by a newline:

```
CHROM  POS  <id_start><id_end>  REF  ALT  60  .  AT=<walk 0>,<walk 1>,...;LV=<LV>[;PS=<ps_start><ps_end>]
```

where `<walk k>` is its visits in text (first, inner visits, last), `<id_start><id_end>` and `<ps_start><ps_end>` two
visits in text each, and `;PS=…` is present only with `PS`. CHROM is the contig's (header).

### `dromos view`'s text

With a contig (`chr20` or `CHM13#0#chr20`) and optionally a region `START-END` (1-based, inclusive; commas allowed),
`dromos view` prints pgcatalog view's text of that contig's catalog: the VCF header with that contig alone, then the
overlapping records in order. Without one, every contig: one header with a `##contig` line each, then each contig's
records in turn. The header:

```
##fileformat=VCFv4.2
##INFO=<ID=AT,Number=R,Type=String,Description="Allele Traversal as path in graph">
##INFO=<ID=LV,Number=1,Type=Integer,Description="Level in the snarl tree (0=top level)">
##INFO=<ID=PS,Number=1,Type=String,Description="ID of variant corresponding to parent snarl">
##contig=<ID=<CHROM>,length=<length>>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
```

## How the writer cuts blocks

(Not needed to read a catalog.) Records are placed in POS order. A record starts a new block when its VCF line is
longer than `build.own_block_bytes` (1 MiB: it has a block of its own, which no later record joins), or the block
filling has one of its own, or it falls in another `build.block_bp` (100,000) window of POS, or the block's lines
already add up to `build.block_bytes` (4 MiB, newlines counted). These are pgcatalog's PGSITES2 rules: the blocks are
its blocks. A column over 1 MiB is compressed with long-distance matching and a window as long as it (2^20 to 2^31),
one over 16 MiB with zstd's own threads (whose output does not depend on how many); a column is stored as it is when
zstd does not make it smaller. The header and index are compressed at level 9 with checksums.
