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
# 7. A GBZ given in place of its files (GBZ_SMALL): made once beside it, remade when older than it, beside the output
# when the graph's folder is not writable, and a clear error otherwise.
ok() { "$@" && echo "ok    $label" || echo "FAIL  $label"; }
mkdir -p gw gr; cp $GBZ_SMALL gw/g.gbz; cp $GBZ_SMALL gr/g.gbz; chmod 555 gr
$D index --reference-sample CHM13 gw/g.gbz empty.gaf e1.gz 2> e1.err
label="a graph's reference map made beside it"; ok test -s gw/g.gbz.CHM13.refmap
$PGCATALOG ref-map --graph $GBZ_SMALL --reference-sample CHM13 --out p.refmap
label="its reference map is pgcatalog's"; ok cmp -s gw/g.gbz.CHM13.refmap p.refmap
$D index --reference-sample CHM13 gw/g.gbz empty.gaf e2.gz 2> e2.err
label="a current reference map reused"; ok bash -c "! grep -q 'kept for the next run' e2.err"
touch -d "2000-01-01" gw/g.gbz.CHM13.refmap
$D index --reference-sample CHM13 gw/g.gbz empty.gaf e3.gz 2> e3.err
label="a reference map older than the graph remade"; ok grep -q "kept for the next run" e3.err
label="no temporary files left"; ok bash -c "! ls gw | grep -q part"
$D index --reference-sample GRCh38 gw/g.gbz empty.gaf e4.gz 2> /dev/null
label="each reference sample its own map"; ok test -s gw/g.gbz.GRCh38.refmap
fails "a graph declaring two reference samples, none chosen" $D index gw/g.gbz empty.gaf x.gz
fails "a reference sample the graph lacks" $D index --reference-sample NOBODY gw/g.gbz empty.gaf x.gz
fails "a reference map of another sample" $D index --reference-sample GRCh38 gw/g.gbz.CHM13.refmap empty.gaf x.gz
$D index --reference-sample CHM13 gr/g.gbz empty.gaf e5.gz 2> /dev/null
label="a read-only graph folder: the reference map beside the output"; ok test -s g.gbz.CHM13.refmap
fails "a read-only graph folder: no node table without an output" $D gam-to-gaf gr/g.gbz $SMALL_GAM
printf "@x\n" | $D gam-to-gaf gw/g.gbz - > /dev/null 2>&1
$PGCATALOG node-table --graph $GBZ_SMALL --out p.nodes
label="a graph's node table made beside it, pgcatalog's"; ok cmp -s gw/g.gbz.nodes p.nodes
chmod 755 gr
# 8. catalog and view on the small graph (GBZ_SMALL, its committed snarls), against pgcatalog build and view: odd
# snarls files, regions and contigs, threads; a missing snarls cache made with vg (VG, the vg binary pgcatalog runs).
S=$OLDPWD/tests/data/micb-kir3dl1.catalog/snarls
[ -d "$S" ] || S=$HOME/dromos/tests/data/micb-kir3dl1.catalog/snarls
cat_pair() {  # name, snarls folder, region options...: dromos catalog and pgcatalog build, the same files or failures
  local name=$1 snarls=$2; shift 2
  rm -rf ca cb; mkdir -p ca cb; cp -r $snarls ca/snarls; cp -r $snarls cb/snarls
  $D catalog --graph $GBZ_SMALL --reference-sample CHM13 --out ca "$@" > ca.err 2>&1; local drc=$?
  $PGCATALOG build --graph $GBZ_SMALL --reference-sample CHM13 --out cb "$@" > cb.err 2>&1; local prc=$?
  if [ $drc -ne 0 ] || [ $prc -ne 0 ]; then
    [ $drc -ne 0 ] && [ $prc -ne 0 ] && echo "ok    $name: both fail ($(tail -1 ca.err))" \
      || echo "FAIL  $name: dromos rc=$drc, pgcatalog rc=$prc"; return
  fi
  local n=0 bad=0
  for f in cb/*.pgsites; do n=$((n + 1)); cmp -s $f ca/$(basename $f) || bad=$((bad + 1)); done
  [ $(ls ca/*.pgsites 2> /dev/null | wc -l) -eq $n ] && [ $bad -eq 0 ] && echo "ok    $name ($n catalogs)" || echo "FAIL  $name"
}
cat_pair "catalog: both contigs" $S --region chr6 --region chr19 --threads 2
for t in 1 3 7; do cat_pair "catalog: $t threads" $S --region chr6 -t $t; done
cat_pair "catalog: a region" $S --region chr6:31,500,000-31,520,000
cat_pair "catalog: a region named by sample" $S --region CHM13#0#chr19:54,800,000-54,900,000
cat_pair "catalog: a region past the contig's end" $S --region chr6:200000000-200000100
cat_pair "catalog: a region before the contig's fragment" $S --region chr6:1-1000
cat_pair "catalog: a one-base region" $S --region chr6:31500001-31500001
cat_pair "catalog: the same contig twice" $S --region chr6 --region chr6:31,500,000-31,600,000
cat_pair "catalog: a contig absent" $S --region chrZZ
cat_pair "catalog: another sample's region" $S --region GRCh38#0#chr6
cat_pair "catalog: a region start of 0" $S --region chr6:0-5
mkdir -p odd; cp $S/graph.id odd/
: > odd/chr6.snarls; cat_pair "catalog: an empty snarls file" odd --region chr6
gzip -dc $S/chr6.snarls > odd/chr6.snarls; cat_pair "catalog: a plain snarls file" odd --region chr6
gzip -dc $S/chr6.snarls | gzip -c > odd/chr6.snarls; cat_pair "catalog: a gzip (not BGZF) snarls file" odd --region chr6
cat $S/chr6.snarls $S/chr19.snarls > odd/chr6.snarls; cat_pair "catalog: two contigs' snarls in one file" odd --region chr6
head -c 600 $S/chr6.snarls > odd/chr6.snarls; cat_pair "catalog: a snarls file cut short" odd --region chr6
gzip -dc $S/chr6.snarls | head -c 500 > odd/chr6.snarls; cat_pair "catalog: a plain snarls file cut short" odd --region chr6
printf 'size=1 fnv1a=0\n' > odd/graph.id; cp $S/chr6.snarls odd/; cat_pair "catalog: another graph's snarls cache" odd --region chr6
if [ -n "${VG:-}" ]; then
  mkdir -p none; cat_pair "catalog: no snarls cached (vg chunk and snarls)" none --region chr6 --region chr19 --vg $VG -t 4
  cmp -s ca/snarls/chr6.snarls cb/snarls/chr6.snarls && cmp -s ca/snarls/chr19.snarls cb/snarls/chr19.snarls \
    && echo "ok    catalog: vg's snarls cached as pgcatalog caches them" || echo "FAIL  catalog: the snarls cache"
  label="catalog: vg's work folder removed"; ok bash -c "! ls -a ca | grep -q dromos-work"
  cat_pair "catalog: a vg that fails" none --region chr6 --vg /bin/false
fi
fails "catalog: no graph" $D catalog --out cx
fails "catalog: a graph that is not one" $D catalog --graph $SMALL_GAF --out cx --reference-sample CHM13
rm -rf cv; mkdir -p cv; cp -r $S cv/snarls
$D catalog --graph $GBZ_SMALL --reference-sample CHM13 --out cv --region chr6 --region chr6:200000000-200000100 -t 2 \
  2> /dev/null
V=cv/chr6.pgsites
view_pair() {  # name, catalog, region...: dromos view prints what pgcatalog view prints, or both fail
  local name=$1 cat=$2; shift 2
  $D view $cat "$@" > vd.out 2> vd.err; local drc=$?
  $PGCATALOG view $cat "$@" > vp.out 2> vp.err; local prc=$?
  if [ $drc -ne 0 ] || [ $prc -ne 0 ]; then
    [ $drc -ne 0 ] && [ $prc -ne 0 ] && echo "ok    $name: both fail ($(tail -1 vd.err))" \
      || echo "FAIL  $name: dromos rc=$drc, pgcatalog rc=$prc"; return
  fi
  cmp -s vd.out vp.out && echo "ok    $name" || echo "FAIL  $name"
}
view_pair "view: whole" $V
view_pair "view: a region" $V chr6:31,500,000-31,520,000
view_pair "view: a region by sample" $V CHM13#0#chr6:31500000-31600000
view_pair "view: the contig alone" $V chr6
view_pair "view: another contig" $V chr19:1-100000000
view_pair "view: another sample's contig" $V GRCh38#0#chr6
view_pair "view: past the end" $V chr6:900000000-900000001
view_pair "view: one base" $V chr6:31500001-31500001
view_pair "view: a region without an end" $V chr6:5
view_pair "view: a region with a bad start" $V chr6:x-5
view_pair "view: an empty catalog" cv/chr6_200000000-200000100.pgsites
printf 'PGSITES1xxxxxxxxxxxxxxxxxxxxxxxx' > old.pgsites; view_pair "view: an older catalog format" old.pgsites
fails "view: not a catalog" $D view $SMALL_GAF
fails "view: a missing catalog" $D view /nonexistent.pgsites
head -c 2000 $V > cut.pgsites; fails "view: a catalog cut short" $D view cut.pgsites
