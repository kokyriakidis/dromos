#!/bin/bash
# dromos's computed snarls of a whole graph against vg snarls' files of its contigs, with compare.py.
#
#   tests/integration/compare_genome.sh DROMOS GRAPH.gbz SNARLS_DIR OUT_DIR [CONTIG...]
#
# SNARLS_DIR is a snarls cache (<contig>.snarls: vg chunk -C and vg snarls of the contig's component(s)). The graph's
# snarls are computed once (`dromos snarls`) and split by component; each contig is compared with the components
# whose node ids overlap the span of its snarls' start nodes (compare.py keeps those its snarls touch). Writes
# OUT_DIR/<contig>.cmp, prints a line per contig, then `EXACT n of m`: n of the m contigs exact.
set -euo pipefail
[ $# -ge 4 ] || { sed -n 2,9p "$0"; exit 2; }
dromos=$1; graph=$2; dir=$3; out=$4; shift 4
here=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$out/comp"
if [ ! -e "$out/comp/.split" ]; then
  "$dromos" snarls --graph "$graph" -o "$out/all.tsv" 2> "$out/snarls.err"
  rm -f "$out"/comp/*.tsv
  awk -F'\t' -v d="$out/comp" '{ f = d "/" $1 ".tsv"; print > f }' "$out/all.tsv"
  # Each component's smallest and largest boundary node.
  awk -F'\t' '{a=int($2/2); b=int($3/2); if(!($1 in lo)){lo[$1]=a; hi[$1]=a} if(a<lo[$1])lo[$1]=a; if(b<lo[$1])lo[$1]=b;
      if(a>hi[$1])hi[$1]=a; if(b>hi[$1])hi[$1]=b} END {for (c in lo) print c"\t"lo[c]"\t"hi[c]}' "$out/all.tsv" \
    > "$out/comp/ranges"
  rm -f "$out/all.tsv"
  touch "$out/comp/.split"
fi
contigs=("$@")
[ ${#contigs[@]} -gt 0 ] || contigs=($(ls "$dir" | sed -n 's/\.snarls$//p'))
exact=0
for c in "${contigs[@]}"; do
  "$dromos" snarls --vg-file "$dir/$c.snarls" -o "$out/$c.vg.tsv" -t 4
  # The components whose span of node ids overlaps that of vg's start nodes.
  comps=$(awk -F'\t' 'NR==FNR {lo[$1]=$2; hi[$1]=$3; next} {n=int($2/2); if (FNR==1 || n<a) a=n; if (FNR==1 || n>b) b=n}
      END {for (k in lo) if (lo[k]<=b && a<=hi[k]) print k}' "$out/comp/ranges" "$out/$c.vg.tsv" | sort -n | paste -sd' ')
  : > "$out/$c.mine.tsv"
  for k in $comps; do cat "$out/comp/$k.tsv" >> "$out/$c.mine.tsv"; done
  python3 -I "$here/compare.py" "$out/$c.vg.tsv" "$out/$c.mine.tsv" > "$out/$c.cmp"
  verdict=$(tail -1 "$out/$c.cmp")
  [ "$verdict" == "EXACT" ] && exact=$((exact + 1))
  echo "$c (component $comps): $verdict; $(grep 'snarls:' "$out/$c.cmp" | sed 's/^ *//'); $(grep 'parents:' \
    "$out/$c.cmp" | sed 's/^ *//' | cut -d';' -f1)"
  rm -f "$out/$c.vg.tsv" "$out/$c.mine.tsv"
done
echo "EXACT $exact of ${#contigs[@]}"
