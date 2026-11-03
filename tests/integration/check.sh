#!/bin/bash
# dromos's end-to-end checks: its outputs against vg's and pgcat's, byte for byte.
#
#   tests/integration/check.sh ENV_FILE [--full]
#
# ENV_FILE sets (paths on the machine with the data):
#   DROMOS     the dromos binary to check (default: build/dromos)
#   PGCAT      pgcat from pgphase's rescue-experiment branch, ff2b4b2 (tools/pgcat; `cargo build --release --no-default-features`;
#              the last with index-gaf --ref-map --slim), the indexer dromos
#              replaced, which makes the reference indexes (pgcatalog, its successor, has no index-gaf)
#   REFMAP     the graph's reference map (pgcatalog ref-map); NODES its node table with letters, LENGTHS one without
#   SMALL_GAM  a few thousand alignments; SMALL_GAF  vg convert -G of it (any order)
#   SLICE_GAM  a GAM ending at a group boundary, some GB; SLICE_SLIM  dromos gam-to-gaf --slim of it, as checked;
#   SLICE_FULL dromos gam-to-gaf of it, as checked (only with --full)
#   PGCATALOG  pgcatalog (pgphase's tools/pgcatalog, plegma-msa branch, 60a3d09, `cargo build --release`), whose
#              surject dromos surject ports and whose ref-map and node-table dromos's must match; GBZ the graph the
#              alignments are on (REFMAP and NODES are its); GBZ_SMALL a small GBZ (tests/data/micb-kir3dl1.gbz)
#   CATALOG_DIR the graph's catalogs (<contig>.pgsites); SURJECT_CONTIG a contig the alignments lie on (chr10)
#   WORK       a directory for outputs and the pgcat references (made once, kept)
#
# Every index is made at -l 0 (stored blocks), where dromos and pgcat write the same bytes: the files must be equal,
# and their .tbi equal decompressed.
set -u
[ $# -ge 1 ] || { sed -n 2,20p "$0"; exit 2; }
source "$1"
FULL=${2:-}
DROMOS=${DROMOS:-build/dromos}
mkdir -p "$WORK"
fails=0
pass() { echo "ok    $1"; }
fail() { echo "FAIL  $1"; fails=$((fails + 1)); }
same_index() {  # name, dromos index, pgcat index
  if cmp -s "$2" "$3" && cmp -s <(bgzip -dc "$2.tbi") <(bgzip -dc "$3.tbi"); then pass "$1"; else fail "$1"; fi
}
# surject's SAM is large (a read's letters and qualities each): records are compared by checksum, the @PG lines
# (each tool's own) left out.
records() { grep -v '^@PG'; }
pgcatalog_surject() {  # checksum file, GAF, options...: pgcatalog surject's records, made once
  local out=$1 gaf=$2; shift 2
  [ -s "$out" ] || (set -o pipefail; "$PGCATALOG" surject -t 10 --graph "$GBZ" --gaf "$gaf" "$@" 2> "$out.err" \
    | records | md5sum > "$out") || { rm -f "$out"; echo "pgcatalog failed on $gaf"; exit 1; }
}
same_surject() {  # name, pgcatalog checksum file, input, [sort], dromos options...: the same records
  local name=$1 want=$2 in=$3 order=$4; shift 4
  local got
  if [ "$order" == "sort" ]; then
    got=$("$DROMOS" surject "$@" "$REFMAP" "$NODES" "$in" 2> /dev/null | records | LC_ALL=C sort | md5sum)
  else
    got=$("$DROMOS" surject "$@" "$REFMAP" "$NODES" "$in" 2> /dev/null | records | md5sum)
  fi
  [ "$got" == "$(cat "$want")" ] && pass "$name" || fail "$name"
}
pgcat_index() {  # GAF, output: pgcat's slim index at level 0, made once
  [ -s "$2" ] || "$PGCAT" index-gaf --ref-map "$REFMAP" --gaf "$1" --slim --level 0 -t 10 --out "$2" \
    --tmp-dir "$WORK" 2> "$2.err" || { echo "pgcat failed on $1"; exit 1; }
}

# The graph files: dromos ref-map, node-table and reference-samples write what pgcatalog writes, on both graphs, for
# each declared reference sample.
pgcatalog_file() {  # output, pgcatalog args...: pgcatalog's file, made once
  local out=$1; shift
  [ -s "$out" ] || "$PGCATALOG" "$@" --out "$out" 2> "$out.err" || { echo "pgcatalog failed: $*"; exit 1; }
}
for g in small:$GBZ_SMALL hprc:$GBZ; do
  n=${g%%:*}; f=${g#*:}
  "$PGCATALOG" reference-samples "$f" > "$WORK/$n.samples.want"
  "$DROMOS" reference-samples "$f" | cmp -s - "$WORK/$n.samples.want" && pass "reference-samples: $n" \
    || fail "reference-samples: $n"
  for s in $(cat "$WORK/$n.samples.want"); do
    pgcatalog_file "$WORK/pgc.$n.$s.refmap" ref-map --graph "$f" --reference-sample "$s"
    "$DROMOS" ref-map --graph "$f" --reference-sample "$s" --out "$WORK/d.$n.$s.refmap" 2> /dev/null
    cmp -s "$WORK/d.$n.$s.refmap" "$WORK/pgc.$n.$s.refmap" && pass "ref-map: $n, $s" || fail "ref-map: $n, $s"
  done
  rm -f "$WORK"/d.$n.*.refmap
  pgcatalog_file "$WORK/pgc.$n.nodes" node-table --graph "$f"
  "$DROMOS" node-table --graph "$f" --out "$WORK/d.$n.nodes" 2> /dev/null
  cmp -s "$WORK/d.$n.nodes" "$WORK/pgc.$n.nodes" && pass "node-table: $n" || fail "node-table: $n"
  pgcatalog_file "$WORK/pgc.$n.lengths" node-table --graph "$f" --no-letters
  "$DROMOS" node-table -t 3 --graph "$f" --no-letters --out "$WORK/d.$n.lengths" 2> /dev/null
  cmp -s "$WORK/d.$n.lengths" "$WORK/pgc.$n.lengths" && pass "node-table --no-letters: $n" \
    || fail "node-table --no-letters: $n"
  rm -f "$WORK/d.$n.nodes" "$WORK/d.$n.lengths"
done

# gam-to-gaf: vg's lines, and the slices as checked before.
"$DROMOS" gam-to-gaf "$NODES" "$SMALL_GAM" > "$WORK/small.gaf"
cmp -s <(sort "$WORK/small.gaf") <(sort "$SMALL_GAF") && pass "gam-to-gaf: vg convert's lines" || fail "gam-to-gaf: vg convert's lines"
"$DROMOS" gam-to-gaf --slim "$LENGTHS" "$SLICE_GAM" | cmp -s - "$SLICE_SLIM" && pass "gam-to-gaf --slim: the slice" || fail "gam-to-gaf --slim: the slice"
if [ "$FULL" == "--full" ]; then
  "$DROMOS" gam-to-gaf "$NODES" "$SLICE_GAM" | cmp -s - "$SLICE_FULL" && pass "gam-to-gaf: the slice, full" || fail "gam-to-gaf: the slice, full"
fi

# index from the GAM, and from its slim GAF in every encoding.
pgcat_index "$SLICE_SLIM" "$WORK/ref.slice.gz"
rm -f "$WORK"/d.*.gz "$WORK"/d.*.gz.tbi
"$DROMOS" index -l 0 "$REFMAP" "$SLICE_GAM" "$WORK/d.gam.gz" 2> /dev/null
same_index "index: the slice GAM" "$WORK/d.gam.gz" "$WORK/ref.slice.gz"
[ -s "$WORK/slice.gaf.gz" ] || gzip -1 -c "$SLICE_SLIM" > "$WORK/slice.gaf.gz"
[ -s "$WORK/slice.gaf.bgz" ] || bgzip -@ 10 -c "$SLICE_SLIM" > "$WORK/slice.gaf.bgz"
"$DROMOS" index -l 0 --gaf "$REFMAP" "$SLICE_SLIM" "$WORK/d.plain.gz" 2> /dev/null
same_index "index: the slice's GAF, plain" "$WORK/d.plain.gz" "$WORK/ref.slice.gz"
"$DROMOS" index -l 0 "$REFMAP" "$WORK/slice.gaf.gz" "$WORK/d.gzip.gz" 2> /dev/null
same_index "index: the slice's GAF, gzip" "$WORK/d.gzip.gz" "$WORK/ref.slice.gz"
"$DROMOS" index -l 0 "$REFMAP" "$WORK/slice.gaf.bgz" "$WORK/d.bgzf.gz" 2> /dev/null
same_index "index: the slice's GAF, BGZF" "$WORK/d.bgzf.gz" "$WORK/ref.slice.gz"
"$DROMOS" index -l 0 --gaf "$REFMAP" - "$WORK/d.stdin.gz" < "$SLICE_SLIM" 2> /dev/null
same_index "index: the slice's GAF, stdin" "$WORK/d.stdin.gz" "$WORK/ref.slice.gz"
# The kind told from the bytes, not the name or a flag: a GAF named .gam, gzip named .gam, and on stdin.
ln -sf "$SLICE_SLIM" "$WORK/slice-gaf.gam"
"$DROMOS" index -l 0 "$REFMAP" "$WORK/slice-gaf.gam" "$WORK/d.named.gz" 2> /dev/null
same_index "index: a GAF named .gam, told by its bytes" "$WORK/d.named.gz" "$WORK/ref.slice.gz"
ln -sf "$WORK/slice.gaf.gz" "$WORK/slice-gz.gam"
"$DROMOS" index -l 0 "$REFMAP" "$WORK/slice-gz.gam" "$WORK/d.namedgz.gz" 2> /dev/null
same_index "index: a gzip GAF named .gam, told by its bytes" "$WORK/d.namedgz.gz" "$WORK/ref.slice.gz"
"$DROMOS" index -l 0 "$REFMAP" - "$WORK/d.stdin2.gz" < "$SLICE_SLIM" 2> /dev/null
same_index "index: a GAF on stdin, told by its bytes" "$WORK/d.stdin2.gz" "$WORK/ref.slice.gz"
"$DROMOS" index -l 0 "$REFMAP" - "$WORK/d.gamstdin.gz" < "$SLICE_GAM" 2> /dev/null
same_index "index: a GAM on stdin, told by its bytes" "$WORK/d.gamstdin.gz" "$WORK/ref.slice.gz"

# index from vg's full GAF (its tags dropped), and from edge cases: tags kept and dropped, CRLF, empty and comment
# lines, names starting with @, a short line, no final newline.
pgcat_index "$SMALL_GAF" "$WORK/ref.small.gz"
"$DROMOS" index -l 0 --gaf "$REFMAP" "$SMALL_GAF" "$WORK/d.small.gz" 2> /dev/null
same_index "index: vg's full GAF" "$WORK/d.small.gz" "$WORK/ref.small.gz"
E="$WORK/edge.gaf"
{
  echo "# a comment"
  echo ""
  head -40 "$SMALL_GAF" | awk 'BEGIN{OFS="\t"} NR%4==1{print $0, "rq:i:7", "zz:i:1", "pu:Z:no_anchor", "rs:Z:chr1:5:+:10M"; next} NR%4==2{printf "%s\r\n", $0; next} NR%4==3{$1="@" $1; print; next} {print}'
  printf "short\t5\t0\t5\t+\n"
  head -45 "$SMALL_GAF" | tail -1 | tr -d '\n'
} > "$E"
pgcat_index "$E" "$WORK/ref.edge.gz"
"$DROMOS" index -l 0 --gaf "$REFMAP" "$E" "$WORK/d.edge.gz" 2> /dev/null
same_index "index: edge-case GAF" "$WORK/d.edge.gz" "$WORK/ref.edge.gz"

# surject: pgcatalog surject's records, from the GAF, the GAM (vg's GAF lines are in another order: both sorted) and
# the GAF compressed; with graph placement and with the contig's catalog.
SC=$SURJECT_CONTIG
CAT="$CATALOG_DIR/$SC.pgsites"
pgcatalog_surject "$WORK/surject.small.md5" "$SMALL_GAF" --contig "$SC"
same_surject "surject: vg's GAF" "$WORK/surject.small.md5" "$SMALL_GAF" keep --contig "$SC"
[ -s "$WORK/small.gaf.bgz" ] || bgzip -c "$SMALL_GAF" > "$WORK/small.gaf.bgz"
same_surject "surject: vg's GAF, BGZF, 3 threads" "$WORK/surject.small.md5" "$WORK/small.gaf.bgz" keep \
  --contig "$SC" -t 3
[ -s "$WORK/surject.small.sorted.md5" ] || (set -o pipefail; "$PGCATALOG" surject --graph "$GBZ" --gaf "$SMALL_GAF" \
  --contig "$SC" 2> /dev/null | records | LC_ALL=C sort | md5sum > "$WORK/surject.small.sorted.md5") || exit 1
same_surject "surject: the GAM" "$WORK/surject.small.sorted.md5" "$SMALL_GAM" sort --contig "$SC"
pgcatalog_surject "$WORK/surject.small.gp.md5" "$SMALL_GAF" --contig "$SC" --graph-placement --sam-contig "$SC"
same_surject "surject --graph-placement: vg's GAF" "$WORK/surject.small.gp.md5" "$SMALL_GAF" keep --contig "$SC" \
  --graph-placement --sam-contig "$SC"
pgcatalog_surject "$WORK/surject.small.cat.md5" "$SMALL_GAF" --contig "$SC" --catalog "$CAT"
same_surject "surject --catalog: vg's GAF" "$WORK/surject.small.cat.md5" "$SMALL_GAF" keep --contig "$SC" \
  --catalog "$CAT"
pgcatalog_surject "$WORK/surject.small.other.md5" "$SMALL_GAF" --contig chr20
same_surject "surject: a contig the reads are not on" "$WORK/surject.small.other.md5" "$SMALL_GAF" keep --contig chr20

if [ "$FULL" == "--full" ]; then
  pgcatalog_surject "$WORK/surject.full.md5" "$SLICE_FULL" --contig "$SC" --sam-contig "$SC"
  same_surject "surject: the slice GAM" "$WORK/surject.full.md5" "$SLICE_GAM" keep --contig "$SC" --sam-contig "$SC"
  same_surject "surject: the slice's full GAF" "$WORK/surject.full.md5" "$SLICE_FULL" keep --contig "$SC" \
    --sam-contig "$SC"
  pgcatalog_surject "$WORK/surject.full.gp.md5" "$SLICE_FULL" --contig "$SC" --sam-contig "$SC" --graph-placement
  same_surject "surject --graph-placement: the slice GAM" "$WORK/surject.full.gp.md5" "$SLICE_GAM" keep \
    --contig "$SC" --sam-contig "$SC" --graph-placement
  pgcatalog_surject "$WORK/surject.full.cat.md5" "$SLICE_FULL" --contig "$SC" --sam-contig "$SC" --catalog "$CAT"
  same_surject "surject --catalog: the slice GAM" "$WORK/surject.full.cat.md5" "$SLICE_GAM" keep --contig "$SC" \
    --sam-contig "$SC" --catalog "$CAT"
  pgcat_index "$SLICE_FULL" "$WORK/ref.full.gz"
  "$DROMOS" index -l 0 --gaf "$REFMAP" "$SLICE_FULL" "$WORK/d.full.gz" 2> /dev/null
  same_index "index: the slice's full GAF" "$WORK/d.full.gz" "$WORK/ref.full.gz"
  cmp -s "$WORK/d.full.gz" "$WORK/d.gam.gz" && pass "index: the full GAF gives the GAM's file" || fail "index: the full GAF gives the GAM's file"
fi
echo "$fails failed"
[ $fails -eq 0 ]
