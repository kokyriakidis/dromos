#!/bin/bash
# Adversarial inputs: dromos index against pgcat index-gaf --slim at level 0 (empty and odd files, gzip and BGZF
# shapes, the chunk boundary), the same output whatever the threads or stdin, and damaged input failing loudly; the
# graph's cache made, reused, made again and refused as it should be; catalogs and views of odd inputs.
#
#   tests/integration/edges.sh ENV_FILE        (the ENV_FILE check.sh takes, after check.sh has run: its pgcat
#                                               reference map and the graph's cache in WORK; WORK/edges for scratch)
set -u
[ $# -ge 1 ] || { sed -n 2,5p "$0"; exit 2; }
source "$1"
DROMOS=${DROMOS:-build/dromos}
D=$DROMOS; P=$PGCAT; M=$WORK/pgc.hprc.CHM13.refmap
[ -s "$M" ] || "$PGCATALOG" ref-map --graph "$GBZ" --reference-sample CHM13 --out "$M" 2> /dev/null
mkdir -p "$WORK/graph"; ln -sf "$GBZ" "$WORK/graph/hprc.gbz"
DG="--graph $WORK/graph/hprc.gbz --reference-sample CHM13"; DN="--graph $WORK/graph/hprc.gbz"
W=$WORK/edges; rm -rf $W; mkdir -p $W; cd $W
same() {  # name, input, [dromos options]
  local name=$1 in=$2; shift 2
  rm -f d.gz* p.gz*
  $P index-gaf --ref-map $M --gaf $in --slim --level 0 -t 4 --out p.gz --tmp-dir $W > p.err 2>&1; local prc=$?
  $D index -l 0 "$@" $DG $in d.gz > d.err 2>&1; local drc=$?
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
for t in 1 3 7; do rm -f t$t.gz*; $D index -l 0 -t $t --gaf $DG $G t$t.gz 2>/dev/null; done
cmp -s t1.gz t3.gz && cmp -s t1.gz t7.gz && echo "ok    the same index with 1, 3 and 7 threads" || echo "FAIL  thread counts differ"
for t in 1 3; do $D gam-to-gaf -t $t $DN $SMALL_GAM > g$t.gaf 2>/dev/null; done
cmp -s g1.gaf g3.gaf && echo "ok    gam-to-gaf the same with 1 and 3 threads" || echo "FAIL  gam-to-gaf thread counts differ"
# 5. A GAM on stdin.
rm -f s.gz* f.gz*; $D index -l 0 $DG - s.gz < $SMALL_GAM 2>/dev/null; $D index -l 0 $DG $SMALL_GAM f.gz 2>/dev/null
cmp -s s.gz f.gz && echo "ok    a GAM on stdin indexes as from its file" || echo "FAIL  GAM on stdin"
$D gam-to-gaf $DN - < $SMALL_GAM | cmp -s - g1.gaf && echo "ok    gam-to-gaf from stdin" || echo "FAIL  gam-to-gaf from stdin"
# 6. Damaged input must fail, not write a partial index.
head -c 100000 $SMALL_GAM > cut.gam; rm -f x.gz*; fails "a GAM cut short" $D index $DG cut.gam x.gz; [ -e x.gz ] && echo "FAIL  an index was left for a GAM cut short"
head -c 3000 m.gaf.gz > cut.gaf.gz; rm -f x.gz*; fails "a gzip GAF cut short" $D index $DG cut.gaf.gz x.gz; [ -e x.gz ] && echo "FAIL  an index was left for a gzip GAF cut short"
printf "garbage" | cat m.gaf.gz - > junk.gaf.gz; fails "junk after a gzip member" $D index $DG junk.gaf.gz x.gz
head -c 30000 g.gaf.bgz > cut.gaf.bgz; fails "a BGZF GAF cut short" $D index $DG cut.gaf.bgz x.gz
python3 -c "
d = bytearray(open('g.gaf.bgz','rb').read()); d[5000] ^= 0xFF; open('bad.gaf.bgz','wb').write(d)"
fails "a BGZF GAF with a flipped byte" $D index $DG bad.gaf.bgz x.gz
cp a.gaf a.gam; rm -f n.gz* o.gz*; $D index -l 0 $DG a.gam n.gz 2>/dev/null; $D index -l 0 --gaf $DG a.gaf o.gz 2>/dev/null
cmp -s n.gz o.gz && echo "ok    a GAF named as a GAM is read as a GAF" || echo "FAIL  a GAF named as a GAM"
fails "a GAF given --gam" $D index --gam $DG a.gaf x.gz
fails "gam-to-gaf given a GAF" $D gam-to-gaf $DN a.gaf
fails "a missing input" $D index $DG /nonexistent.gam x.gz
fails "an unwritable output" $D index $DG $SMALL_GAM /proc/nope/x.gz
fails "a graph that is not a GBZ" $D gam-to-gaf --graph $SMALL_GAF $SMALL_GAM
fails "no graph" $D index $SMALL_GAM x.gz
# 7. The graph's cache (GBZ_SMALL): made on first use beside the graph, reused, made again when the graph changes,
# when it is damaged or of another version, kept elsewhere with --cache, made by two first uses at once; a read-only
# folder, a file that is not a cache and another graph's cache refused, untouched.
ok() { "$@" && echo "ok    $label" || echo "FAIL  $label"; }
made() { grep -q "making the graph's cache" "$1"; }
mkdir -p gw gr; cp $GBZ_SMALL gw/g.gbz; cp $GBZ_SMALL gr/g.gbz; chmod 555 gr
GW="--graph gw/g.gbz --reference-sample CHM13"
$D index $GW empty.gaf e1.gz 2> e1.err
label="a cache made beside the graph on first use"; ok bash -c "test -s gw/g.gbz.dromos && grep -q 'first use' e1.err"
$D index $GW empty.gaf e2.gz 2> e2.err
label="a current cache reused"; ok bash -c "! grep -q making e2.err"
$D index --graph gw/g.gbz --reference-sample GRCh38 empty.gaf e2b.gz 2> e2b.err
label="one cache for both reference samples"; ok bash -c "! grep -q making e2b.err"
cp gw/g.gbz.dromos first.dromos
touch -d "2001-01-01" gw/g.gbz
$D index $GW empty.gaf e3.gz 2> e3.err
label="a graph changed (a new time): its cache made again"; ok grep -q "the graph has changed" e3.err
label="no temporary files left"; ok bash -c "! ls gw | grep -q part"
head -c 5000 first.dromos > gw/g.gbz.dromos
$D gam-to-gaf --graph gw/g.gbz $SMALL_GAM > /dev/null 2> e4.err || echo "FAIL  gam-to-gaf failed: $(tail -1 e4.err)"
label="a cache cut short: made again"; ok grep -q "(damaged)" e4.err
$D gam-to-gaf --graph gw/g.gbz $SMALL_GAM > /dev/null 2> e5.err || echo "FAIL  gam-to-gaf failed: $(tail -1 e5.err)"
label="then reused"; ok bash -c "! grep -q making e5.err"
python3 -c "
d = bytearray(open('gw/g.gbz.dromos','rb').read()); d[200] ^= 1; open('gw/g.gbz.dromos','wb').write(d)"
$D gam-to-gaf --graph gw/g.gbz $SMALL_GAM > /dev/null 2> e6.err || echo "FAIL  gam-to-gaf failed: $(tail -1 e6.err)"
label="a cache with a damaged header: made again"; ok grep -q "(damaged)" e6.err
python3 -c "
d = bytearray(open('gw/g.gbz.dromos','rb').read()); d[8] = 0; open('gw/g.gbz.dromos','wb').write(d)"
$D gam-to-gaf --graph gw/g.gbz $SMALL_GAM > /dev/null 2> e7.err || echo "FAIL  gam-to-gaf failed: $(tail -1 e7.err)"
label="another version's cache: made again"; ok grep -q "another dromos version" e7.err
echo "not a cache" > gw/g.gbz.dromos
fails "a file that is not a cache, at the cache's place" $D gam-to-gaf --graph gw/g.gbz $SMALL_GAM
label="that file left as it was"; ok grep -q "not a cache" gw/g.gbz.dromos
rm gw/g.gbz.dromos
fails "another graph's cache" $D gam-to-gaf --graph $GBZ --cache first.dromos $SMALL_GAM
label="another graph's cache left as it was"; ok bash -c "test \$(stat -c %s first.dromos) -gt 5000"
$D index --graph gr/g.gbz --reference-sample CHM13 empty.gaf x.gz 2> ro.err; rc=$?
label="a read-only folder: an error naming --cache ($(tail -1 ro.err | cut -c1-60)...)"
ok bash -c "[ $rc -ne 0 ] && grep -q -- --cache ro.err && ! test -e gr/g.gbz.dromos"
$D index --graph gr/g.gbz --reference-sample CHM13 --cache elsewhere.dromos empty.gaf e8.gz 2> /dev/null
label="a read-only folder: the cache elsewhere (--cache)"; ok bash -c "test -s elsewhere.dromos && ! test -e gr/g.gbz.dromos"
chmod 755 gr; $D prepare --graph gr/g.gbz 2> /dev/null; chmod 555 gr
$D gam-to-gaf --graph gr/g.gbz $SMALL_GAM > ro.gaf 2> e9.err || echo "FAIL  gam-to-gaf failed: $(tail -1 e9.err)"
label="a read-only folder with its cache prepared: read"; ok bash -c "test -s ro.gaf && ! grep -q making e9.err"
chmod 755 gr
mkdir -p gc; cp $GBZ_SMALL gc/g.gbz
for i in 1 2 3 4; do $D index --graph gc/g.gbz --reference-sample CHM13 $G c$i.gz 2> c$i.err & done; wait
label="four first uses at once: the same index, one cache, no temporary files"
ok bash -c "cmp -s c1.gz c2.gz && cmp -s c1.gz c3.gz && cmp -s c1.gz c4.gz && test -s gc/g.gbz.dromos \
  && ! ls gc | grep -q part"
$D prepare --graph gc/g.gbz 2> e10.err
label="and that cache current"; ok bash -c "! grep -q making e10.err"
fails "a graph declaring two reference samples, none chosen" $D index --graph gw/g.gbz empty.gaf x.gz
fails "a reference sample the graph lacks" $D index --graph gw/g.gbz --reference-sample NOBODY empty.gaf x.gz
label="no cache left by a sample the graph lacks"; ok bash -c "! ls gw | grep -q part"
fails "prepare: no graph" $D prepare
fails "surject: no contig" $D surject $DG $SMALL_GAF
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
