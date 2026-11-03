#!/bin/bash
# Adversarial inputs: dromos index against pgcat index-gaf --slim at level 0 (empty and odd files, gzip and BGZF
# shapes, the chunk boundary), the same output whatever the threads or stdin, and damaged input failing loudly.
#
#   tests/integration/edges.sh ENV_FILE        (the ENV_FILE check.sh takes; WORK is used for scratch)
set -u
[ $# -ge 1 ] || { sed -n 2,5p "$0"; exit 2; }
source "$1"
DROMOS=${DROMOS:-build/dromos}
D=$DROMOS; P=$PGCAT; M=$REFMAP
W=$WORK/edges; rm -rf $W; mkdir -p $W; cd $W
same() {  # name, input, [dromos options]
  local name=$1 in=$2; shift 2
  rm -f d.gz* p.gz*
  $P index-gaf --ref-map $M --gaf $in --slim --level 0 -t 4 --out p.gz --tmp-dir $W > p.err 2>&1; local prc=$?
  $D index -l 0 "$@" $M $in d.gz > d.err 2>&1; local drc=$?
  if [ $prc -ne 0 ] || [ $drc -ne 0 ]; then
    echo "??    $name: pgcat rc=$prc, dromos rc=$drc"; tail -1 p.err; tail -1 d.err; return
  fi
  if cmp -s d.gz p.gz && cmp -s <(bgzip -dc d.gz.tbi) <(bgzip -dc p.gz.tbi); then echo "ok    $name"; else echo "FAIL  $name"; ls -l d.gz p.gz; fi
}
fails() {  # name, command...: must exit non-zero
  local name=$1; shift
  "$@" > f.out 2> f.err && echo "FAIL  $name: succeeded" || echo "ok    $name: $(tail -1 f.err)"
}
G=$SMALL_GAF
# 1. Inputs with nothing to index.
: > empty.gaf; same "empty GAF" empty.gaf
printf "@HD\tVN:Z:1.0\n# nothing\n" > header.gaf; same "header-only GAF" header.gaf
awk 'BEGIN{OFS="\t"} {$6=">999999999999"; print}' $G | head -50 > missing.gaf; same "only missing-node lines" missing.gaf
awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,"*","*","*","*","*","*","*","0"}' $G | head -50 > unmapped.gaf; same "only unmapped lines" unmapped.gaf
# 2. Compression shapes.
head -2000 $G > a.gaf; tail -n +2001 $G > b.gaf
gzip -c a.gaf > m.gaf.gz; gzip -c b.gaf >> m.gaf.gz; same "two-member gzip" m.gaf.gz
gzip -c empty.gaf > e.gaf.gz; same "gzip of nothing" e.gaf.gz
bgzip -c $G > g.gaf.bgz; same "BGZF" g.gaf.bgz
cat g.gaf.bgz g.gaf.bgz > gg.gaf.bgz; same "two BGZF files concatenated (EOF block between)" gg.gaf.bgz
python3 - <<PY
import zlib, struct
def block(data):
    c = zlib.compressobj(6, zlib.DEFLATED, -15); body = c.compress(data) + c.flush()
    size = len(body) + 26
    return bytes([31,139,8,4,0,0,0,0,0,255,6,0,66,67,2,0]) + struct.pack("<H", size - 1) + body + struct.pack("<II", zlib.crc32(data), len(data))
lines = open("$G", "rb").read().split(b"\n")[:300]
text = b"\n".join(lines) + b"\n"
out = block(b"")
for i in range(0, len(text), 7000):
    out += block(text[i:i+7000]) + block(b"")
out += block(b"")
open("eb.gaf.bgz", "wb").write(out)
PY
same "BGZF with empty blocks, lines cut across blocks" eb.gaf.bgz
# 3. Sizes against the 64 MiB chunk: a line across the boundary, and exactly one chunk.
python3 - <<PY
src = open("$SLICE_SLIM", "rb")
chunk = 64 << 20
data = src.read(chunk + (8 << 20))
cut = data.rfind(b"\n", 0, chunk + (1 << 20)) + 1
open("over.gaf", "wb").write(data[:cut])
# exactly one chunk: pad the last line's name so the file is 64 MiB to the byte
lines = data[:chunk].split(b"\n")[:-1]
body = b"\n".join(lines[:-1]) + b"\n"
last = lines[-1]
pad = chunk - len(body) - len(last) - 1
name, rest = last.split(b"\t", 1)
exact = body + name + b"x" * pad + b"\t" + rest + b"\n"
assert len(exact) == chunk, len(exact)
open("exact.gaf", "wb").write(exact)
PY
same "a line across the chunk boundary" over.gaf
same "exactly one chunk" exact.gaf
head -c -1 exact.gaf > exactnl.gaf; same "exactly one chunk less its last newline" exactnl.gaf
# 4. Threads: the same file however many.
for t in 1 3 7; do rm -f t$t.gz*; $D index -l 0 -t $t --gaf $M $G t$t.gz 2>/dev/null; done
cmp -s t1.gz t3.gz && cmp -s t1.gz t7.gz && echo "ok    the same index with 1, 3 and 7 threads" || echo "FAIL  thread counts differ"
for t in 1 3; do $D gam-to-gaf -t $t $NODES $SMALL_GAM > g$t.gaf 2>/dev/null; done
cmp -s g1.gaf g3.gaf && echo "ok    gam-to-gaf the same with 1 and 3 threads" || echo "FAIL  gam-to-gaf thread counts differ"
# 5. A GAM on stdin.
rm -f s.gz* f.gz*; $D index -l 0 $M - s.gz < $SMALL_GAM 2>/dev/null; $D index -l 0 $M $SMALL_GAM f.gz 2>/dev/null
cmp -s s.gz f.gz && echo "ok    a GAM on stdin indexes as from its file" || echo "FAIL  GAM on stdin"
$D gam-to-gaf $NODES - < $SMALL_GAM | cmp -s - g1.gaf && echo "ok    gam-to-gaf from stdin" || echo "FAIL  gam-to-gaf from stdin"
# 6. Damaged input must fail, not write a partial index.
head -c 100000 $SMALL_GAM > cut.gam; rm -f x.gz*; fails "a GAM cut short" $D index $M cut.gam x.gz; [ -e x.gz ] && echo "FAIL  an index was left for a GAM cut short"
head -c 3000 m.gaf.gz > cut.gaf.gz; rm -f x.gz*; fails "a gzip GAF cut short" $D index $M cut.gaf.gz x.gz; [ -e x.gz ] && echo "FAIL  an index was left for a gzip GAF cut short"
printf "garbage" | cat m.gaf.gz - > junk.gaf.gz; fails "junk after a gzip member" $D index $M junk.gaf.gz x.gz
head -c 30000 g.gaf.bgz > cut.gaf.bgz; fails "a BGZF GAF cut short" $D index $M cut.gaf.bgz x.gz
python3 -c "
d = bytearray(open('g.gaf.bgz','rb').read()); d[5000] ^= 0xFF; open('bad.gaf.bgz','wb').write(d)"
fails "a BGZF GAF with a flipped byte" $D index $M bad.gaf.bgz x.gz
cp a.gaf a.gam; rm -f n.gz* o.gz*; $D index -l 0 $M a.gam n.gz 2>/dev/null; $D index -l 0 --gaf $M a.gaf o.gz 2>/dev/null
cmp -s n.gz o.gz && echo "ok    a GAF named as a GAM is read as a GAF" || echo "FAIL  a GAF named as a GAM"
fails "a GAF given --gam" $D index --gam $M a.gaf x.gz
fails "gam-to-gaf given a GAF" $D gam-to-gaf $NODES a.gaf
fails "a missing input" $D index $M /nonexistent.gam x.gz
fails "an unwritable output" $D index $M $SMALL_GAM /proc/nope/x.gz
fails "a node table that is not one" $D gam-to-gaf $SMALL_GAF $SMALL_GAM
fails "a ref map that is not one" $D index $SMALL_GAF $SMALL_GAM x.gz
