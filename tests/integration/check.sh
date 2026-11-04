#!/bin/bash
# dromos's end-to-end checks: its outputs against vg's and pgcat's, byte for byte.
#
#   tests/integration/check.sh ENV_FILE [--full] [--only GROUP,...]
#
# Groups: cache (reference-samples, prepare, --cache), catalog (catalog, view, catalog-info), snarls (computed and
# --dist snarls against vg snarls), convert (gam-to-gaf), index, surject. --only runs those groups alone: run the ones
# a change touches. --full adds the whole-genome and slice-wide cases of the groups run. Each group's time is printed.
#
# ENV_FILE sets (paths on the machine with the data):
#   DROMOS     the dromos binary to check (default: build/dromos)
#   PGCAT      pgcat from pgphase's rescue-experiment branch, ff2b4b2 (tools/pgcat; `cargo build --release --no-default-features`;
#              the last with index-gaf --ref-map --slim), the indexer dromos
#              replaced, which makes the reference indexes (pgcatalog, its successor, has no index-gaf)
#   GBZ        the graph the alignments are on; GBZ_SMALL a small GBZ (tests/data/micb-kir3dl1.gbz)
#   SMALL_GAM  a few thousand alignments; SMALL_GAF  vg convert -G of it (any order)
#   SLICE_GAM  a GAM ending at a group boundary, some GB; SLICE_SLIM  dromos gam-to-gaf --slim of it, as checked;
#   SLICE_FULL dromos gam-to-gaf of it, as checked (only with --full)
#   PGCATALOG  pgcatalog (pgphase's tools/pgcatalog, plegma-msa branch, 60a3d09, `cargo build --release`), whose
#              surject, build and view dromos's surject, catalog and view port, and whose ref-map makes pgcat's
#              reference map
#   CATALOG_DIR pgcatalog's catalogs of the graph (<contig>.pgsites); SURJECT_CONTIG a contig the alignments lie on
#              (chr10)
#   SNARLS     the graph's snarls cache (<contig>.snarls and graph.id, as pgcatalog build keeps it), read, not changed:
#              vg snarls, the oracle dromos's computed snarls are checked against
#   DIST       the graph's distance index (giraffe's .dist; optional: with --full, dromos snarls --dist of it is checked
#              against vg snarls and the computed snarls, every contig)
#   WORK       a directory for outputs and the pgcat and pgcatalog references (made once, kept), and for the graph's
#              cache: dromos reads GBZ through a link in WORK/graph, its cache made beside it
#   SAMTOOLS   samtools (1.24; default: samtools), whose sort and index surject -o's BAM is checked against
#   CARVED     a folder of carved slices of the whole GAM, chr20s.gam and chr6s.gam with their .gaf (only with --full)
#
# Every index is made at -l 0 (stored blocks), where dromos and pgcat write the same bytes: the files must be equal,
# and their .tbi equal decompressed.
set -u
[ $# -ge 1 ] || { sed -n 2,28p "$0"; exit 2; }
source "$1"
shift
FULL="" ONLY=""
while [ $# -gt 0 ]; do
  case $1 in
    --full) FULL=--full ;;
    --only) ONLY=$2; shift ;;
    *) echo "unknown option $1"; exit 2 ;;
  esac
  shift
done
DROMOS=${DROMOS:-build/dromos}
SAMTOOLS=${SAMTOOLS:-samtools}
HERE=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$WORK"
fails=0
pass() { echo "ok    $1"; }
fail() { echo "FAIL  $1"; fails=$((fails + 1)); }
want() { [ -z "$ONLY" ] || [[ ",$ONLY," == *",$1,"* ]]; }
GROUP="" STARTED=0
group() {  # name: the group that follows (and the time of the one before)
  [ -n "$GROUP" ] && echo "==    $GROUP: $((SECONDS - STARTED)) s"
  GROUP=$1 STARTED=$SECONDS
}
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
    got=$("$DROMOS" surject --graph "$G" --reference-sample CHM13 "$@" "$in" 2> /dev/null | records | LC_ALL=C sort \
      | md5sum)
  else
    got=$("$DROMOS" surject --graph "$G" --reference-sample CHM13 "$@" "$in" 2> /dev/null | records | md5sum)
  fi
  [ "$got" == "$(cat "$want")" ] && pass "$name" || fail "$name"
}
REFMAP="$WORK/pgc.hprc.CHM13.refmap"
pgcat_index() {  # GAF, output: pgcat's slim index at level 0, made once (and pgcat's reference map, pgcatalog ref-map's
  # file, for it)
  [ -s "$REFMAP" ] || "$PGCATALOG" ref-map --graph "$GBZ" --reference-sample CHM13 --out "$REFMAP" 2> "$REFMAP.err" \
    || { echo "pgcatalog ref-map failed"; exit 1; }
  [ -s "$2" ] || "$PGCAT" index-gaf --ref-map "$REFMAP" --gaf "$1" --slim --level 0 -t 10 --out "$2" \
    --tmp-dir "$WORK" 2> "$2.err" || { echo "pgcat failed on $1"; exit 1; }
}
# The graph, through a link in WORK/graph, its cache beside it.
mkdir -p "$WORK/graph"
G="$WORK/graph/hprc.gbz"
ln -sf "$GBZ" "$G"

if want cache; then
group cache
# reference-samples prints what pgcatalog prints, on both graphs.
for g in small:$GBZ_SMALL hprc:$GBZ; do
  n=${g%%:*}; f=${g#*:}
  [ -s "$WORK/$n.samples.want" ] || "$PGCATALOG" reference-samples "$f" > "$WORK/$n.samples.want"
  "$DROMOS" reference-samples "$f" | cmp -s - "$WORK/$n.samples.want" && pass "reference-samples: $n" \
    || fail "reference-samples: $n"
done
# The graph's cache: made by dromos prepare beside the graph, current after, and the same file when a command makes it
# on first use, beside the graph or elsewhere (--cache).
rm -f "$G".* "$WORK/elsewhere.dromos"
"$DROMOS" prepare --graph "$G" 2> "$WORK/prepare.err" && [ -s "$G.dromos" ] && grep -q "making the graph's cache" \
  "$WORK/prepare.err" && pass "prepare: the graph's cache made beside it" || fail "prepare: the graph's cache made beside it"
"$DROMOS" prepare --graph "$G" 2> "$WORK/prepare.err" && ! grep -q "making" "$WORK/prepare.err" \
  && pass "prepare: a current cache left as it is" || fail "prepare: a current cache left as it is"
pgcat_index "$SLICE_SLIM" "$WORK/ref.slice.gz"
"$DROMOS" index -l 0 --graph "$G" --reference-sample CHM13 --cache "$WORK/elsewhere.dromos" "$SLICE_SLIM" \
  "$WORK/d.graph.gz" 2> /dev/null
same_index "index: the cache made on first use, elsewhere (--cache)" "$WORK/d.graph.gz" "$WORK/ref.slice.gz"
cmp -s "$WORK/elsewhere.dromos" "$G.dromos" && pass "the cache made on first use is prepare's, byte for byte" \
  || fail "the cache made on first use is prepare's, byte for byte"
rm -f "$WORK/elsewhere.dromos" "$WORK/d.graph.gz" "$WORK/d.graph.gz.tbi"
fi
# The other groups read the graph's cache: made now if missing (current otherwise: nothing is done).
"$DROMOS" prepare --graph "$G" 2> /dev/null || { echo "dromos prepare failed"; exit 1; }

# catalog and view: dromos view of dromos's one catalog prints, for each contig and region, what pgcatalog view prints
# of pgcatalog build's catalog of it from the same snarls (the small graph from its committed snarls and from snarls
# dromos computes; HPRC contigs, a region, and the whole genome); the same file whatever the threads. Computed snarls
# must give the catalog vg snarls give: the same bytes but for the header, whose fields are the same but for its
# snarls.* keys.
pgcatalog_catalog() {  # folder, region options...: pgcatalog's catalogs, made once
  local out=$1; shift
  [ -e "$out/done" ] || { rm -rf "$out"; mkdir -p "$out" && cp -r "$SNARLS_OF" "$out/snarls" \
    && "$PGCATALOG" build --graph "$GRAPH_OF" --reference-sample CHM13 --out "$out" "$@" 2> "$out/err" \
    && touch "$out/done"; } || { echo "pgcatalog build failed: $*"; exit 1; }
}
dromos_catalog() {  # file, options...: dromos's catalog, from SNARLS_OF (linked, never written: every contig is in it)
  local out=$1; shift
  local from
  from=$(cd "$SNARLS_OF" && pwd)
  rm -rf "$out" "$out.snarls"; mkdir -p "$out.snarls"
  for f in "$from"/*; do ln -s "$f" "$out.snarls/"; done
  "$DROMOS" catalog --graph "$GRAPH_OF" --reference-sample CHM13 --snarls-cache "$out.snarls" --out "$out" "$@" \
    2> "$out.err"
}
same_views() {  # name, dromos catalog, pgcatalog folder: each pgcatalog catalog's view, dromos's of its contig/region,
  # the contigs at once (JOBS, 8)
  local name=$1 cat=$2 n=0 bad f b r tmp
  tmp=$(mktemp -d "$WORK/views.XXXXXX")
  for f in $(ls -S "$3"/*.pgsites); do
    b=$(basename "$f" .pgsites)
    r=$b
    [[ "$b" == *_*-* ]] && r="${b%_*}:${b##*_}"
    n=$((n + 1))
    { cmp -s <("$DROMOS" view -t 2 "$cat" "$r" 2> /dev/null) <("$PGCATALOG" view "$f" 2> /dev/null) \
      || touch "$tmp/$b"; } &
    while [ "$(jobs -rp | wc -l)" -ge "${JOBS:-8}" ]; do wait -n || true; done
  done
  wait
  bad=$(ls "$tmp" | wc -l)
  rm -rf "$tmp"
  [ $n -gt 0 ] && [ "$bad" -eq 0 ] && pass "$name ($n contigs)" || fail "$name ($bad of $n differ)"
}
same_catalog() {  # name, dromos catalog, dromos catalog: byte for byte the same but for the header (and the preamble
  # that locates it), whose fields are the same but for the snarls.* keys
  local a b ha hb
  a=$("$DROMOS" catalog-info "$2" 2> /dev/null | grep -v '^snarls\.\|^file' | md5sum)
  b=$("$DROMOS" catalog-info "$3" 2> /dev/null | grep -v '^snarls\.\|^file' | md5sum)
  ha=$(od -An -t u8 -j 24 -N 8 "$2" | tr -d ' ')
  hb=$(od -An -t u8 -j 24 -N 8 "$3" | tr -d ' ')
  [ -n "$ha" ] && [ "$ha" == "$hb" ] && [ "$a" == "$b" ] \
    && cmp -s <(head -c "$ha" "$2" | tail -c +65) <(head -c "$hb" "$3" | tail -c +65) && pass "$1" || fail "$1"
}
same_view() {  # name, dromos catalog, pgcatalog catalog, region...: dromos view prints what pgcatalog view prints
  local name=$1 cat=$2 pcat=$3; shift 3
  cmp -s <("$DROMOS" view "$cat" "$@" 2> /dev/null) <("$PGCATALOG" view "$pcat" "$@" 2> /dev/null) && pass "$name" \
    || fail "$name"
}
SMALL_SNARLS=tests/data/micb-kir3dl1.catalog/snarls
ALL=$WORK/dcat.all.sites
if want catalog; then
group catalog
GRAPH_OF=$GBZ_SMALL SNARLS_OF=$SMALL_SNARLS
pgcatalog_catalog "$WORK/pcat.small" --region chr6 --region chr19 --threads 2
dromos_catalog "$WORK/dcat.small.sites" --region chr6 --region chr19 --threads 2
same_views "catalog: the small graph" "$WORK/dcat.small.sites" "$WORK/pcat.small"
bad=0
for c in chr6 chr19; do
  cmp -s <("$DROMOS" view "$WORK/dcat.small.sites" $c) tests/data/micb-kir3dl1.catalog/$c.vcf || bad=1
done
[ $bad -eq 0 ] && pass "catalog: the small graph's committed views" || fail "catalog: the small graph's committed views"
"$DROMOS" catalog --graph "$GBZ_SMALL" --reference-sample CHM13 --region chr6 --region chr19 -t 2 \
  --out "$WORK/dcat.small.computed.sites" 2> /dev/null
same_views "catalog: the small graph, snarls computed" "$WORK/dcat.small.computed.sites" "$WORK/pcat.small"
same_catalog "catalog: the small graph, computed snarls as vg snarls" "$WORK/dcat.small.computed.sites" \
  "$WORK/dcat.small.sites"
GRAPH_OF=$GBZ SNARLS_OF=$SNARLS
for c in chr20 chr6 chr1 chrY "chr20:10,000,000-12,000,000"; do
  d=$(echo "$c" | tr ':,' '__')
  pgcatalog_catalog "$WORK/pcat.$d" --region "$c" --threads 10
  dromos_catalog "$WORK/dcat.$d.sites" --region "$c" -t 10
  same_views "catalog: $c" "$WORK/dcat.$d.sites" "$WORK/pcat.$d"
done
"$DROMOS" catalog --graph "$GBZ" --reference-sample CHM13 --region chr20 -t 10 --out "$WORK/dcat.chr20.computed.sites" \
  2> /dev/null
same_catalog "catalog: chr20, computed snarls as vg snarls" "$WORK/dcat.chr20.computed.sites" "$WORK/dcat.chr20.sites"
dromos_catalog "$WORK/dcat.chr20_3.sites" --region chr20 -t 3
cmp -s "$WORK/dcat.chr20_3.sites" "$WORK/dcat.chr20.sites" && pass "catalog: chr20, 3 threads, the same file" \
  || fail "catalog: chr20, 3 threads, the same file"
dromos_catalog "$WORK/dcat.two.sites" --region chr20:10,000,000-12,000,000 --region chrY -t 10
same_views "catalog: two contigs, one a region" "$WORK/dcat.two.sites" "$WORK/pcat.chr20_10_000_000-12_000_000"
same_views "catalog: two contigs, the other whole" "$WORK/dcat.two.sites" "$WORK/pcat.chrY"
P20=$WORK/pcat.chr20/chr20.pgsites
same_view "view: chr20, whole" "$WORK/dcat.chr20.sites" "$P20"
same_view "view: chr20, a region" "$WORK/dcat.chr20.sites" "$P20" "chr20:10,000,000-12,000,000"
same_view "view: chr20, a region named by sample" "$WORK/dcat.chr20.sites" "$P20" "CHM13#0#chr20:30000000-30001000"
same_view "view: chr20, a contig named alone" "$WORK/dcat.chr20.sites" "$P20" chr20
same_view "view: chr6, the MHC" "$WORK/dcat.chr6.sites" "$WORK/pcat.chr6/chr6.pgsites" chr6:29,000,000-34,000,000
same_view "view: chr1, past its end" "$WORK/dcat.chr1.sites" "$WORK/pcat.chr1/chr1.pgsites" chr1:300000000-300000100
# The whole genome, once, from snarls dromos computes (the default); with --full, each contig's view against
# pgcatalog's catalogs (made from vg snarls: so the computed snarls give vg snarls' catalog; chr2's giant record is in
# its contig's view).
"$DROMOS" catalog --graph "$GBZ" --reference-sample CHM13 -t 10 --out "$ALL" 2> "$ALL.err"
"$DROMOS" catalog-info "$ALL" | grep -q "^graph.fingerprint" && pass "catalog-info: the whole genome's header" \
  || fail "catalog-info: the whole genome's header"
"$DROMOS" catalog-info "$ALL" | grep -q "^snarls.source.dromos" && pass "catalog-info: the snarls computed" \
  || fail "catalog-info: the snarls computed"
[ "$FULL" == "--full" ] && same_views "catalog: the whole genome, computed snarls" "$ALL" "$CATALOG_DIR"
find "$WORK" -maxdepth 1 -name 'dcat.*' ! -name 'dcat.all.sites' -exec rm -rf {} +
fi

if want snarls; then
group snarls
# snarls --dist: the small graph's distance index gives vg snarls' snarls and the computed snarls (compare.py: the same
# snarls, parents, and chains or children's order), and the same table checked against its graph or not.
exact() {  # name, table, table: compare.py says EXACT
  [ "$(python3 -I "$HERE/compare.py" "$2" "$3" | tail -1)" == EXACT ] && pass "$1" || fail "$1"
}
SN=$WORK/snarls.small
rm -rf "$SN"; mkdir -p "$SN"
"$DROMOS" snarls --dist tests/data/micb-kir3dl1.dist --graph "$GBZ_SMALL" -o "$SN/dist.tsv" 2> /dev/null
"$DROMOS" snarls --dist tests/data/micb-kir3dl1.dist -t 3 -o "$SN/alone.tsv" 2> /dev/null
"$DROMOS" snarls --vg-file tests/data/micb-kir3dl1.vg.snarls -o "$SN/vg.tsv" 2> /dev/null
"$DROMOS" snarls --graph "$GBZ_SMALL" -o "$SN/computed.tsv" 2> /dev/null
[ -s "$SN/dist.tsv" ] && cmp -s "$SN/dist.tsv" "$SN/alone.tsv" \
  && pass "snarls --dist: the small graph's, the same table checked against its graph or not" \
  || fail "snarls --dist: the small graph's, the same table checked against its graph or not"
exact "snarls --dist: the small graph's as vg snarls'" "$SN/vg.tsv" "$SN/dist.tsv"
exact "snarls --dist: the small graph's as the computed snarls" "$SN/dist.tsv" "$SN/computed.tsv"
# With --full, every contig (tests/integration/compare_genome.sh, contigs at once): the computed snarls against vg
# snarls' files, then the distance index's against the computed ones. The two say the index's are vg's: no third
# comparison. vg's files as tables are references, made once (WORK/snarls.vg).
if [ "$FULL" == "--full" ]; then
  n=$(ls "$SNARLS" | grep -c '\.snarls$')
  rm -rf "$WORK/snarls.cmp"
  VG_TABLES=$WORK/snarls.vg "$HERE/compare_genome.sh" "$DROMOS" "$GBZ" "$SNARLS" "$WORK/snarls.cmp" \
    > "$WORK/snarls.cmp.log" 2>&1
  tail -1 "$WORK/snarls.cmp.log" | grep -qx "EXACT $n of $n" \
    && pass "snarls: the whole genome's computed snarls as vg snarls' ($n contigs)" \
    || fail "snarls: the whole genome's computed snarls as vg snarls' ($(tail -1 "$WORK/snarls.cmp.log"))"
  if [ -n "${DIST:-}" ]; then
    "$DROMOS" snarls --dist "$DIST" --graph "$GBZ" -t 10 -o "$WORK/dist.tsv" 2> "$WORK/dist.err" \
      && pass "snarls --dist: the whole genome's distance index, checked against its graph" \
      || fail "snarls --dist: the whole genome's distance index, checked against its graph"
    VG_TABLES=$WORK/snarls.vg ORACLE=$WORK/dist.tsv "$HERE/compare_genome.sh" "$DROMOS" "$GBZ" "$SNARLS" \
      "$WORK/snarls.cmp" > "$WORK/dist.computed.log" 2>&1
    tail -1 "$WORK/dist.computed.log" | grep -qx "EXACT $n of $n" \
      && pass "snarls --dist: the whole genome's distance index as the computed snarls ($n contigs)" \
      || fail "snarls --dist: the whole genome's distance index as the computed snarls ($(tail -1 \
        "$WORK/dist.computed.log"))"
    rm -rf "$WORK/dist.tsv" "$WORK/snarls.cmp/oracle"
  fi
fi
fi

if want convert; then
group convert
# gam-to-gaf: vg's lines, and the slices as checked before.
"$DROMOS" gam-to-gaf --graph "$G" "$SMALL_GAM" > "$WORK/small.gaf"
cmp -s <(sort "$WORK/small.gaf") <(sort "$SMALL_GAF") && pass "gam-to-gaf: vg convert's lines" || fail "gam-to-gaf: vg convert's lines"
"$DROMOS" gam-to-gaf --slim --graph "$G" "$SLICE_GAM" | cmp -s - "$SLICE_SLIM" && pass "gam-to-gaf --slim: the slice" || fail "gam-to-gaf --slim: the slice"
if [ "$FULL" == "--full" ]; then
  "$DROMOS" gam-to-gaf --graph "$G" "$SLICE_GAM" | cmp -s - "$SLICE_FULL" && pass "gam-to-gaf: the slice, full" || fail "gam-to-gaf: the slice, full"
fi
fi

if want index; then
group index
# index from the GAM, and from its slim GAF in every encoding.
DG="--graph $G --reference-sample CHM13"
pgcat_index "$SLICE_SLIM" "$WORK/ref.slice.gz"
rm -f "$WORK"/d.*.gz "$WORK"/d.*.gz.tbi
"$DROMOS" index -l 0 $DG "$SLICE_GAM" "$WORK/d.gam.gz" 2> /dev/null
same_index "index: the slice GAM" "$WORK/d.gam.gz" "$WORK/ref.slice.gz"
[ -s "$WORK/slice.gaf.gz" ] || gzip -1 -c "$SLICE_SLIM" > "$WORK/slice.gaf.gz"
[ -s "$WORK/slice.gaf.bgz" ] || bgzip -@ 10 -c "$SLICE_SLIM" > "$WORK/slice.gaf.bgz"
"$DROMOS" index -l 0 --gaf $DG "$SLICE_SLIM" "$WORK/d.plain.gz" 2> /dev/null
same_index "index: the slice's GAF, plain" "$WORK/d.plain.gz" "$WORK/ref.slice.gz"
"$DROMOS" index -l 0 $DG "$WORK/slice.gaf.gz" "$WORK/d.gzip.gz" 2> /dev/null
same_index "index: the slice's GAF, gzip" "$WORK/d.gzip.gz" "$WORK/ref.slice.gz"
"$DROMOS" index -l 0 $DG "$WORK/slice.gaf.bgz" "$WORK/d.bgzf.gz" 2> /dev/null
same_index "index: the slice's GAF, BGZF" "$WORK/d.bgzf.gz" "$WORK/ref.slice.gz"
"$DROMOS" index -l 0 --gaf $DG - "$WORK/d.stdin.gz" < "$SLICE_SLIM" 2> /dev/null
same_index "index: the slice's GAF, stdin" "$WORK/d.stdin.gz" "$WORK/ref.slice.gz"
# The kind told from the bytes, not the name or a flag: a GAF named .gam, gzip named .gam, and on stdin.
ln -sf "$SLICE_SLIM" "$WORK/slice-gaf.gam"
"$DROMOS" index -l 0 $DG "$WORK/slice-gaf.gam" "$WORK/d.named.gz" 2> /dev/null
same_index "index: a GAF named .gam, told by its bytes" "$WORK/d.named.gz" "$WORK/ref.slice.gz"
ln -sf "$WORK/slice.gaf.gz" "$WORK/slice-gz.gam"
"$DROMOS" index -l 0 $DG "$WORK/slice-gz.gam" "$WORK/d.namedgz.gz" 2> /dev/null
same_index "index: a gzip GAF named .gam, told by its bytes" "$WORK/d.namedgz.gz" "$WORK/ref.slice.gz"
"$DROMOS" index -l 0 $DG - "$WORK/d.stdin2.gz" < "$SLICE_SLIM" 2> /dev/null
same_index "index: a GAF on stdin, told by its bytes" "$WORK/d.stdin2.gz" "$WORK/ref.slice.gz"
"$DROMOS" index -l 0 $DG - "$WORK/d.gamstdin.gz" < "$SLICE_GAM" 2> /dev/null
same_index "index: a GAM on stdin, told by its bytes" "$WORK/d.gamstdin.gz" "$WORK/ref.slice.gz"

# index from vg's full GAF (its tags dropped), and from edge cases: tags kept and dropped, CRLF, empty and comment
# lines, names starting with @, a short line, no final newline.
pgcat_index "$SMALL_GAF" "$WORK/ref.small.gz"
"$DROMOS" index -l 0 --gaf $DG "$SMALL_GAF" "$WORK/d.small.gz" 2> /dev/null
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
"$DROMOS" index -l 0 --gaf $DG "$E" "$WORK/d.edge.gz" 2> /dev/null
same_index "index: edge-case GAF" "$WORK/d.edge.gz" "$WORK/ref.edge.gz"
if [ "$FULL" == "--full" ]; then
  pgcat_index "$SLICE_FULL" "$WORK/ref.full.gz"
  "$DROMOS" index -l 0 --gaf $DG "$SLICE_FULL" "$WORK/d.full.gz" 2> /dev/null
  same_index "index: the slice's full GAF" "$WORK/d.full.gz" "$WORK/ref.full.gz"
  cmp -s "$WORK/d.full.gz" "$WORK/d.gam.gz" && pass "index: the full GAF gives the GAM's file" || fail "index: the full GAF gives the GAM's file"
fi
fi

if want surject; then
group surject
# surject --catalog reads the whole genome's catalog: the catalog group's, or made now.
[ -s "$ALL" ] || "$DROMOS" catalog --graph "$GBZ" --reference-sample CHM13 -t 10 --out "$ALL" 2> "$ALL.err"
# surject: pgcatalog surject's records, from the GAF, the GAM (vg's GAF lines are in another order: both sorted) and
# the GAF compressed; with graph placement and with the contig's catalog (pgcatalog's, against dromos's of the same
# sites).
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
  --catalog "$ALL"
pgcatalog_surject "$WORK/surject.small.other.md5" "$SMALL_GAF" --contig chr20
same_surject "surject: a contig the reads are not on" "$WORK/surject.small.other.md5" "$SMALL_GAF" keep --contig chr20

# surject -o: one pass over every contig into a BAM sorted as samtools sort sorts and indexed as samtools index
# indexes. Each contig's records (samtools view OUT CONTIG) are its single-contig run's, sorted by samtools; the
# index says what samtools index says of the same file; a single contig's BAM holds the very records samtools makes
# of its SAM.
single_sorted() {  # input, contig, options...: the checksum of the single-contig run's records, sorted by samtools
  local in=$1 c=$2; shift 2
  "$DROMOS" surject --graph "$G" --reference-sample CHM13 --contig "$c" "$@" "$in" 2> /dev/null \
    | "$SAMTOOLS" sort -@ 4 -m 1G -T "$WORK/sorting" -O SAM 2> /dev/null | "$SAMTOOLS" view - | md5sum
}
same_bam() {  # name, BAM, input, all|placed (every contig, or those with records and one without), options...
  local name=$1 bam=$2 in=$3 which=$4; shift 4
  local bad=0 n=0 c empty=""
  "$SAMTOOLS" quickcheck "$bam" || bad=$((bad + 1))
  local idx=$bam.bai flag=""
  [ -e "$bam.csi" ] && idx=$bam.csi flag=-c
  rm -f "$WORK/again.bam"*; cp "$bam" "$WORK/again.bam"
  "$SAMTOOLS" index $flag "$WORK/again.bam" && python3 -I "$HERE/same_index.py" "$idx" "$WORK/again.bam.${idx##*.}" \
    > /dev/null || bad=$((bad + 1))
  cmp -s <("$SAMTOOLS" idxstats "$bam") <("$SAMTOOLS" idxstats "$WORK/again.bam") || bad=$((bad + 1))
  while read -r c len mapped unmapped; do
    [ "$c" == "*" ] && continue
    if [ "$which" == "placed" ] && [ "$mapped" == 0 ]; then
      [ -n "$empty" ] && continue
      empty=$c
    fi
    n=$((n + 1))
    [ "$("$SAMTOOLS" view "$bam" "$c" | md5sum)" == "$(single_sorted "$in" "${c##*#}" "$@")" ] || bad=$((bad + 1))
  done < <("$SAMTOOLS" idxstats "$bam")
  [ $bad -eq 0 ] && pass "$name ($n contigs)" || fail "$name ($bad of $n checks differ)"
}
bam_records() { bgzip -dc "$1" | python3 -I "$HERE/bam_records.py"; }
SO="--graph $G --reference-sample CHM13"
"$DROMOS" surject $SO --contig "$SC" "$SMALL_GAF" -o "$WORK/one.bam" 2> /dev/null
"$DROMOS" surject $SO --contig "$SC" "$SMALL_GAF" 2> /dev/null | "$SAMTOOLS" sort -O BAM -T "$WORK/sorting" \
  -o "$WORK/one.samtools.bam" 2> /dev/null
[ "$(bam_records "$WORK/one.bam")" == "$(bam_records "$WORK/one.samtools.bam")" ] \
  && pass "surject -o, one contig: samtools sort's records, byte for byte" \
  || fail "surject -o, one contig: samtools sort's records, byte for byte"
same_bam "surject -o, one contig: index and records" "$WORK/one.bam" "$SMALL_GAF" all
"$DROMOS" surject $SO "$SMALL_GAF" -o "$WORK/all.bam" 2> /dev/null
same_bam "surject -o: every contig, vg's GAF" "$WORK/all.bam" "$SMALL_GAF" all
"$DROMOS" surject $SO "$SMALL_GAM" -o "$WORK/all.gam.bam" 2> /dev/null
same_bam "surject -o: every contig, the GAM" "$WORK/all.gam.bam" "$SMALL_GAM" all
"$DROMOS" surject $SO "$SMALL_GAF" -o "$WORK/all.spill.bam" -m 1M 2> /dev/null
cmp -s "$WORK/all.bam" "$WORK/all.spill.bam" && cmp -s "$WORK/all.bam.bai" "$WORK/all.spill.bam.bai" \
  && pass "surject -o: spilled to disk, the same BAM and index" || fail "surject -o: spilled to disk, the same BAM and index"
"$DROMOS" surject $SO "$SMALL_GAF" -o "$WORK/all.csi.bam" --csi 2> /dev/null
cmp -s "$WORK/all.bam" "$WORK/all.csi.bam" && ! test -e "$WORK/all.csi.bam.bai" \
  && same_bam "surject -o --csi: a .csi as samtools index -c makes it" "$WORK/all.csi.bam" "$SMALL_GAF" placed \
  || fail "surject -o --csi: a .csi as samtools index -c makes it"
"$DROMOS" surject $SO "$SMALL_GAF" -o "$WORK/all.gp.bam" --graph-placement 2> /dev/null
same_bam "surject -o --graph-placement: every contig" "$WORK/all.gp.bam" "$SMALL_GAF" all --graph-placement
"$DROMOS" surject $SO "$SMALL_GAF" -o "$WORK/all.cat.bam" --catalog "$ALL" 2> /dev/null
same_bam "surject -o --catalog: every contig" "$WORK/all.cat.bam" "$SMALL_GAF" placed --catalog "$ALL"
"$DROMOS" surject $SO "$SMALL_GAF" -o "$WORK/short.bam" --sam-names short 2> /dev/null
[ "$("$SAMTOOLS" view "$WORK/short.bam" "$SC" | cut -f 1,2,4- | md5sum)" \
  == "$("$SAMTOOLS" view "$WORK/all.bam" "CHM13#0#$SC" | cut -f 1,2,4- | md5sum)" ] \
  && "$SAMTOOLS" view -H "$WORK/short.bam" | grep -q "^@SQ.SN:$SC.LN" \
  && pass "surject -o --sam-names short" || fail "surject -o --sam-names short"
# Several contigs as SAM on stdout: each contig's records, in IN's order, are its single-contig run's.
"$DROMOS" surject $SO --contig chr20 --contig "$SC" "$SMALL_GAF" > "$WORK/two.sam" 2> /dev/null
bad=0
for c in "$SC" chr20; do
  cmp -s <(awk -v c="CHM13#0#$c" '$3 == c' "$WORK/two.sam") <("$DROMOS" surject $SO --contig "$c" "$SMALL_GAF" \
    2> /dev/null | grep -v '^@') || bad=1
done
[ "$(grep '^@SQ' "$WORK/two.sam" | cut -f 2 | tr '\n' ' ')" == "SN:CHM13#0#$SC SN:CHM13#0#chr20 " ] || bad=1
[ $bad -eq 0 ] && pass "surject: two contigs as SAM, each its single run's" || fail "surject: two contigs as SAM, each its single run's"

if [ "$FULL" == "--full" ]; then
  "$DROMOS" surject $SO "$SLICE_GAM" -o "$WORK/slice.bam" 2> /dev/null
  same_bam "surject -o: the slice GAM, every contig" "$WORK/slice.bam" "$SLICE_GAM" all
  "$DROMOS" surject $SO "$SLICE_FULL" -o "$WORK/slice.gaf.bam" 2> /dev/null
  cmp -s "$WORK/slice.bam" "$WORK/slice.gaf.bam" && pass "surject -o: the slice's full GAF gives the GAM's BAM" \
    || fail "surject -o: the slice's full GAF gives the GAM's BAM"
  "$DROMOS" surject $SO "$SLICE_GAM" -o "$WORK/slice.gp.bam" --graph-placement 2> /dev/null
  same_bam "surject -o --graph-placement: the slice GAM" "$WORK/slice.gp.bam" "$SLICE_GAM" placed --graph-placement
  "$DROMOS" surject $SO "$SLICE_GAM" -o "$WORK/slice.cat.bam" --catalog "$ALL" 2> /dev/null
  same_bam "surject -o --catalog: the slice GAM" "$WORK/slice.cat.bam" "$SLICE_GAM" placed --catalog "$ALL"
  "$DROMOS" surject $SO "$SLICE_GAM" -o "$WORK/slice.spill.bam" -m 256M 2> /dev/null
  cmp -s "$WORK/slice.bam" "$WORK/slice.spill.bam" && cmp -s "$WORK/slice.bam.bai" "$WORK/slice.spill.bam.bai" \
    && pass "surject -o: the slice spilled to disk, the same BAM" || fail "surject -o: the slice spilled to disk, the same BAM"
  rm -f "$WORK"/slice.*bam*
  for s in chr20s chr6s; do
    for k in gam gaf; do
      "$DROMOS" surject $SO "$CARVED/$s.$k" -o "$WORK/$s.$k.bam" 2> /dev/null
      same_bam "surject -o: the carved $s, its ${k^^}" "$WORK/$s.$k.bam" "$CARVED/$s.$k" placed
    done
    cmp -s <(bam_records "$WORK/$s.gam.bam") <(bam_records "$WORK/$s.gaf.bam") \
      && pass "surject -o: the carved $s, the GAM gives its GAF's BAM" \
      || fail "surject -o: the carved $s, the GAM gives its GAF's BAM"
    rm -f "$WORK/$s".*.bam*
  done
  pgcatalog_surject "$WORK/surject.full.md5" "$SLICE_FULL" --contig "$SC" --sam-contig "$SC"
  same_surject "surject: the slice GAM" "$WORK/surject.full.md5" "$SLICE_GAM" keep --contig "$SC" --sam-contig "$SC"
  same_surject "surject: the slice's full GAF" "$WORK/surject.full.md5" "$SLICE_FULL" keep --contig "$SC" \
    --sam-contig "$SC"
  pgcatalog_surject "$WORK/surject.full.gp.md5" "$SLICE_FULL" --contig "$SC" --sam-contig "$SC" --graph-placement
  same_surject "surject --graph-placement: the slice GAM" "$WORK/surject.full.gp.md5" "$SLICE_GAM" keep \
    --contig "$SC" --sam-contig "$SC" --graph-placement
  pgcatalog_surject "$WORK/surject.full.cat.md5" "$SLICE_FULL" --contig "$SC" --sam-contig "$SC" --catalog "$CAT"
  same_surject "surject --catalog: the slice GAM" "$WORK/surject.full.cat.md5" "$SLICE_GAM" keep --contig "$SC" \
    --sam-contig "$SC" --catalog "$ALL"
fi
fi
group ""
echo "$fails failed"
[ $fails -eq 0 ]
