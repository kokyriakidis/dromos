#!/bin/bash
# dromos's computed snarls of a whole graph against vg snarls' files of its contigs, with compare.py, the contigs in
# parallel.
#
#   tests/integration/compare_genome.sh DROMOS GRAPH.gbz SNARLS_DIR OUT_DIR [CONTIG...]
#
# SNARLS_DIR is a snarls cache (<contig>.snarls: vg chunk -C and vg snarls of the contig's component(s)). The graph's
# snarls are computed once (`dromos snarls`) and split by component; each contig is compared with the components
# whose node ids overlap the span of its snarls' start nodes (compare.py keeps those its snarls touch). Writes
# OUT_DIR/<contig>.cmp, prints a line per contig, then `EXACT n of m`: n of the m contigs exact.
#   ORACLE=TABLE  compare with TABLE's components overlapping each contig's span (a distance index's: `dromos snarls
#                 --dist`), in place of its vg snarls file (which still names the contigs and their spans);
#                 OUT_DIR/<contig>.oracle.cmp
#   VG_TABLES=DIR each contig's vg snarls file as a table (dromos snarls --vg-file), made once and kept there
#                 [OUT_DIR/vg]: a reference, like the snarls files it is made from
#   JOBS=N        contigs compared at once [8]; each compare.py holds about 2.5 GB
set -euo pipefail
[ $# -ge 4 ] || { sed -n 2,16p "$0"; exit 2; }
dromos=$1; graph=$2; dir=$3; out=$4; shift 4
here=$(cd "$(dirname "$0")" && pwd)
vgt=${VG_TABLES:-$out/vg}
jobs=${JOBS:-8}
split() {  # table, folder: the table's lines by component (<k>.tsv), each component's smallest and largest boundary
  # node (ranges)
  local table=$1 d=$2
  rm -rf "$d"; mkdir -p "$d"
  awk -F'\t' -v d="$d" '{ f = d "/" $1 ".tsv"; print > f }' "$table"
  awk -F'\t' '{a=int($2/2); b=int($3/2); if(!($1 in lo)){lo[$1]=a; hi[$1]=a} if(a<lo[$1])lo[$1]=a; if(b<lo[$1])lo[$1]=b;
      if(a>hi[$1])hi[$1]=a; if(b>hi[$1])hi[$1]=b} END {for (c in lo) print c"\t"lo[c]"\t"hi[c]}' "$table" > "$d/ranges"
  touch "$d/.split"
}
overlapping() {  # folder, vg's table of a contig: the folder's components whose span of node ids overlaps that of
  # vg's start nodes
  awk -F'\t' 'NR==FNR {lo[$1]=$2; hi[$1]=$3; next} {n=int($2/2); if (FNR==1 || n<a) a=n; if (FNR==1 || n>b) b=n}
      END {for (k in lo) if (lo[k]<=b && a<=hi[k]) print k}' "$1/ranges" "$2" | sort -n | paste -sd' '
}
mkdir -p "$out" "$vgt"
mine=$out/comp
if [ ! -e "$mine/.split" ]; then
  "$dromos" snarls --graph "$graph" -o "$out/all.tsv" 2> "$out/snarls.err"
  split "$out/all.tsv" "$mine"
  rm -f "$out/all.tsv"
fi
[ -z "${ORACLE:-}" ] || split "$ORACLE" "$out/oracle"
suffix=${ORACLE:+.oracle}
contigs=("$@")
[ ${#contigs[@]} -gt 0 ] || contigs=($(ls "$dir" | sed -n 's/\.snarls$//p'))
one() {  # contig: its comparison, OUT_DIR/<contig><suffix>.cmp and a line of it in .line
  local c=$1 comps want vg
  vg=$vgt/$c.tsv
  [ -s "$vg" ] || { "$dromos" snarls --vg-file "$dir/$c.snarls" -o "$vg.tmp" -t 2 && mv "$vg.tmp" "$vg"; }
  comps=$(overlapping "$mine" "$vg")
  : > "$out/$c.mine.tsv"
  for k in $comps; do cat "$mine/$k.tsv" >> "$out/$c.mine.tsv"; done
  if [ -n "${ORACLE:-}" ]; then
    want=$(overlapping "$out/oracle" "$vg")
    vg=$out/$c.oracle.tsv
    : > "$vg"
    for k in $want; do cat "$out/oracle/$k.tsv" >> "$vg"; done
  fi
  python3 -I "$here/compare.py" "$vg" "$out/$c.mine.tsv" > "$out/$c$suffix.cmp"
  echo "$c (component $comps): $(tail -1 "$out/$c$suffix.cmp"); $(grep 'snarls:' "$out/$c$suffix.cmp" | sed 's/^ *//');" \
    "$(grep 'parents:' "$out/$c$suffix.cmp" | sed 's/^ *//' | cut -d';' -f1)" > "$out/$c$suffix.line"
  rm -f "$out/$c.mine.tsv" "$out/$c.oracle.tsv"
}
# Biggest contigs first (their vg files), so the last to finish are small.
order=($(for c in "${contigs[@]}"; do echo "$(stat -c %s "$dir/$c.snarls") $c"; done | sort -rn | cut -d' ' -f2))
for c in "${order[@]}"; do
  one "$c" &
  while [ "$(jobs -rp | wc -l)" -ge "$jobs" ]; do wait -n || true; done
done
wait
exact=0
for c in "${contigs[@]}"; do
  cat "$out/$c$suffix.line"
  [ "$(tail -1 "$out/$c$suffix.cmp")" == "EXACT" ] && exact=$((exact + 1))
done
echo "EXACT $exact of ${#contigs[@]}"
