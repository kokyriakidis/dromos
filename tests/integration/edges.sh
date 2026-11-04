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
# reference-fasta refuses a sample the graph lacks (before its cache is looked at) and an output it cannot write, and
# leaves no index of a FASTA it did not write.
$D prepare --graph gw/g.gbz 2> /dev/null; cp gw/g.gbz.dromos rf.before
fails "reference-fasta: a reference sample the graph lacks" $D reference-fasta --graph gw/g.gbz \
  --reference-sample NOBODY -o rf.fa
label="reference-fasta: the graph's cache untouched"; ok bash -c "cmp -s gw/g.gbz.dromos rf.before && ! test -e rf.fa"
fails "reference-fasta: an unwritable output" $D reference-fasta $GW -o /proc/nope/x.fa
$D reference-fasta $GW -o rf.fa 2> /dev/null && chmod 444 rf.fa
fails "reference-fasta: a FASTA it cannot write over" $D reference-fasta $GW -o rf.fa
label="reference-fasta: no index left beside it"; ok bash -c "! test -e rf.fa.fai"
chmod 644 rf.fa
# 7b. surject: every contig without --contig; -o's BAM from odd inputs, whatever the threads; options that do not fit.
SAMTOOLS=${SAMTOOLS:-samtools}
HERE=$OLDPWD/tests/integration
[ -e "$HERE/same_index.py" ] || HERE=$HOME/dromos/tests/integration
$D surject $DG $SMALL_GAF 2> /dev/null | grep -c '^@SQ' > sq.n
label="surject: no contig, every contig in the SAM's header"; ok bash -c "[ \$(cat sq.n) -eq 25 ]"
bam_ok() {  # BAM: samtools reads it, and indexes it as dromos did
  rm -f again.bam*; cp $1 again.bam; $SAMTOOLS quickcheck $1 && $SAMTOOLS index again.bam \
    && python3 -I $HERE/same_index.py $1.bai again.bam.bai > /dev/null \
    && cmp -s <($SAMTOOLS idxstats $1) <($SAMTOOLS idxstats again.bam)
}
$D surject $DG empty.gaf -o s.empty.bam 2> /dev/null
label="surject -o: an empty GAF, a BAM of no records"
ok bash -c "$(declare -f bam_ok); SAMTOOLS=$SAMTOOLS HERE=$HERE bam_ok s.empty.bam && [ \$($SAMTOOLS view -c s.empty.bam) -eq 0 ] \
  && [ \$($SAMTOOLS view -H s.empty.bam | grep -c '^@SQ') -eq 25 ]"
$D surject $DG --contig chr20 $G -o s.none.bam 2> /dev/null
label="surject -o: a contig no read is on"
ok bash -c "$(declare -f bam_ok); SAMTOOLS=$SAMTOOLS HERE=$HERE bam_ok s.none.bam && [ \$($SAMTOOLS view -c s.none.bam) -eq 0 ]"
$D surject $DG --contig chr10 $G -o s.one.bam 2> /dev/null
label="surject -o: a single contig"
ok bash -c "$(declare -f bam_ok); SAMTOOLS=$SAMTOOLS HERE=$HERE bam_ok s.one.bam && [ \$($SAMTOOLS view -H s.one.bam | grep -c '^@SQ') -eq 1 ]"
$D surject $DG $G -o s.t1.bam -t 1 2> /dev/null
$D surject $DG $G -o s.t7.bam -t 7 2> /dev/null
$D surject $DG - -o s.stdin.bam < $G 2> /dev/null
label="surject -o: 1 thread and 7, and stdin, the same BAM and index"
ok bash -c "cmp -s s.t1.bam s.t7.bam && cmp -s s.t1.bam.bai s.t7.bam.bai && cmp -s s.t1.bam s.stdin.bam"
$D surject $DG $G -o s.tiny.bam -m 1K 2> /dev/null
label="surject -o: spilled at every step, the same BAM"; ok cmp -s s.t1.bam s.tiny.bam
$D surject $DG $G -o s.l0.bam -l 0 2> /dev/null
label="surject -o: level 0, the same records"
ok bash -c "cmp -s <(bgzip -dc s.l0.bam) <(bgzip -dc s.t1.bam)"
$D surject $DG $G -o s.cut.bam --gam 2> /dev/null
label="surject -o: a failed run leaves no BAM"; ok bash -c "! test -e s.cut.bam && ! ls | grep -q 's.cut.bam'"
fails "surject: --sam-contig with two contigs" $D surject $DG --contig chr10 --contig chr20 --sam-contig x $G
fails "surject: a contig the graph lacks" $D surject $DG --contig chrNone $G
fails "surject: --sam-names neither full nor short" $D surject $DG --sam-names long $G
fails "surject: -l without -o" $D surject $DG -l 3 $G
fails "surject -o: a folder that does not exist" $D surject $DG $G -o /nonexistent/x.bam
# 8. catalog and view on the small graph (GBZ_SMALL, its committed snarls), against pgcatalog build and view: odd
# snarls files, regions and contigs, threads (the same file whatever their number); a missing snarls cache made with
# vg (VG, the vg binary pgcatalog runs); --dist refused; snarls computed by dromos (the default: any threads, a contig
# with no snarls, a tiny graph); damaged catalogs refused. dromos snarls --dist --graph refuses distance indexes of
# another graph, cut short or of another version.
S=$OLDPWD/tests/data/micb-kir3dl1.catalog/snarls
[ -d "$S" ] || S=$HOME/dromos/tests/data/micb-kir3dl1.catalog/snarls
SD=$(dirname $(dirname $S))/micb-kir3dl1.dist
views_of() {  # dromos catalog, pgcatalog folder: each pgcatalog catalog's view, dromos's of its contig or region
  local n=0 bad=0 f b r
  for f in $2/*.pgsites; do
    b=$(basename $f .pgsites); r=$b
    [[ "$b" == *_*-* ]] && r="${b%_*}:${b##*_}"
    n=$((n + 1))
    cmp -s <($D view $1 "$r" 2> /dev/null) <($PGCATALOG view $f 2> /dev/null) || bad=$((bad + 1))
  done
  [ $n -gt 0 ] && [ $bad -eq 0 ] && echo "$n" || echo "-$bad"
}
cat_pair() {  # name, snarls folder, region options...: dromos catalog and pgcatalog build, the same views or failures
  local name=$1 snarls=$2; shift 2
  rm -rf ca cb; mkdir -p ca cb; cp -r $snarls ca/snarls; cp -r $snarls cb/snarls
  $D catalog --graph $GBZ_SMALL --reference-sample CHM13 --snarls-cache ca/snarls --out ca/c.sites "$@" > ca.err 2>&1
  local drc=$?
  $PGCATALOG build --graph $GBZ_SMALL --reference-sample CHM13 --out cb "$@" > cb.err 2>&1; local prc=$?
  if [ $drc -ne 0 ] || [ $prc -ne 0 ]; then
    [ $drc -ne 0 ] && [ $prc -ne 0 ] && echo "ok    $name: both fail ($(tail -1 ca.err))" \
      || echo "FAIL  $name: dromos rc=$drc, pgcatalog rc=$prc"; return
  fi
  local n=$(views_of ca/c.sites cb)
  [ "$n" -gt 0 ] && echo "ok    $name ($n catalogs)" || echo "FAIL  $name ($n)"
}
cat_pair "catalog: both contigs" $S --region chr6 --region chr19 --threads 2
for t in 1 3 7; do cat_pair "catalog: $t threads" $S --region chr6 --region chr19 -t $t; cp ca/c.sites t$t.sites; done
label="catalog: the same file on 1, 3 and 7 threads"; ok bash -c "cmp -s t1.sites t3.sites && cmp -s t1.sites t7.sites"
cat_pair "catalog: a region" $S --region chr6:31,500,000-31,520,000
cat_pair "catalog: a region named by sample" $S --region CHM13#0#chr19:54,800,000-54,900,000
cat_pair "catalog: a region past the contig's end" $S --region chr6:200000000-200000100
cat_pair "catalog: a region before the contig's fragment" $S --region chr6:1-1000
cat_pair "catalog: a one-base region" $S --region chr6:31500001-31500001
cat_pair "catalog: the same contig twice" $S --region chr6 --region chr6:31,500,000-31,600,000
cat_pair "catalog: two regions of a contig" $S --region chr19:54,800,000-54,850,000 --region chr19:54,840,000-54,900,000
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
  label="catalog: the vg named in the header"; ok bash -c "$D catalog-info ca/c.sites | grep -q '^snarls.vg.vg version'"
  cat_pair "catalog: a vg that fails" none --region chr6 --vg /bin/false
fi
fails "catalog: no graph" $D catalog --out cx
fails "catalog: a graph that is not one" $D catalog --graph $SMALL_GAF --out cx --reference-sample CHM13
# A distance index is no longer the catalog's: --dist refused, naming where it went; one beside the graph not read.
fails "catalog: --dist" $D catalog --graph $GBZ_SMALL --reference-sample CHM13 --dist $SD --out cx
label="catalog: --dist refused, the snarls computed from the graph and dromos snarls --dist named"
ok bash -c "grep -q 'snarls are computed from the graph' f.err && grep -q 'dromos snarls --dist' f.err && ! test -e cx"
rm -rf cb; mkdir -p cb; cp -r $S cb/snarls
$PGCATALOG build --graph $GBZ_SMALL --reference-sample CHM13 --out cb --region chr6 --region chr19 > cb.err 2>&1
mkdir -p dg; ln -sf $GBZ_SMALL dg/small.gbz; cp $SD dg/small.dist
$D catalog --graph dg/small.gbz --reference-sample CHM13 --region chr6 --region chr19 -t 3 > dg.err 2>&1
label="catalog: a .dist beside the graph left alone, the snarls computed"
ok bash -c "test -s dg/small.CHM13.chr6+chr19.sites && $D catalog-info dg/small.CHM13.chr6+chr19.sites \
  | grep -q '^snarls.source.dromos'"
# Snarls computed by dromos: pgcatalog's records from vg snarls, the same file on any threads.
for t in 1 3 7; do
  $D catalog --graph $GBZ_SMALL --reference-sample CHM13 --region chr6 --region chr19 -t $t --out computed$t.sites \
    > computed$t.err 2>&1
done
n=$(views_of computed3.sites cb)
[ "$n" -gt 0 ] && echo "ok    catalog: from computed snarls, vg snarls' records ($n contigs)" \
  || echo "FAIL  catalog: from computed snarls ($n)"
label="catalog: computed snarls, the same file on 1, 3 and 7 threads"
ok bash -c "cmp -s computed1.sites computed3.sites && cmp -s computed1.sites computed7.sites"
label="catalog: computed snarls named in the header"
ok bash -c "$D catalog-info computed1.sites | grep -q '^snarls.ties' && $D catalog-info computed1.sites \
  | grep -q '^snarls.orientation'"
if [ -n "${VG:-}" ]; then
  # Graphs made here (vg gbwt from GFA): computed snarls give the catalog vg chunk and vg snarls give.
  computed_pair() {  # name, GFA, contig...: the catalog from computed snarls and from vg's, the same views
    local name=$1 gfa=$2; shift 2
    local regions=() c
    for c in "$@"; do regions+=(--region "$c"); done
    rm -rf tg; mkdir -p tg/snarls
    $VG gbwt -G $gfa --set-reference CHM13 --gbz-format -g tg/g.gbz > tg/gbwt.err 2>&1 \
      || { echo "??    $name: vg gbwt failed"; return; }
    $D catalog --graph tg/g.gbz --reference-sample CHM13 "${regions[@]}" -t 3 --out tg/c.sites > tg/c.err 2>&1
    local crc=$?
    $D catalog --graph tg/g.gbz --reference-sample CHM13 "${regions[@]}" -t 3 --snarls-cache tg/snarls --vg $VG \
      --out tg/v.sites > tg/v.err 2>&1
    local vrc=$?
    if [ $crc -ne 0 ] || [ $vrc -ne 0 ]; then
      echo "FAIL  $name: computed rc=$crc, vg rc=$vrc ($(tail -1 tg/c.err))"; return
    fi
    local bad=0 records
    for c in "$@"; do cmp -s <($D view tg/c.sites $c) <($D view tg/v.sites $c) || bad=1; done
    records=$($D catalog-info tg/c.sites | awk -F'\t' '$1=="records" {print $2}')
    [ $bad -eq 0 ] && echo "ok    $name ($records records)" || echo "FAIL  $name"
  }
  # A contig with no snarls: one linear path, a second haplotype along it.
  { printf 'H\tVN:Z:1.1\nS\t1\tACGTACGTAC\nS\t2\tGGTTAACC\nS\t3\tTTGACA\n'
    printf 'L\t1\t+\t2\t+\t0M\nL\t2\t+\t3\t+\t0M\n'
    printf 'W\tCHM13\t0\tchrA\t0\t24\t>1>2>3\nW\tHG1\t1\tchrA\t0\t24\t>1>2>3\n'; } > nosnarls.gfa
  computed_pair "catalog: computed snarls, a contig with no snarls" nosnarls.gfa chrA
  label="catalog: computed snarls, a contig with no snarls has no records"
  ok bash -c "$D catalog-info tg/c.sites | grep -qx 'records.0'"
  # A tiny graph: a SNP and an insertion on chrA, a lone node on chrB (a second component).
  { printf 'H\tVN:Z:1.1\nS\t1\tACGTAC\nS\t2\tA\nS\t3\tG\nS\t4\tTTGA\nS\t5\tCCC\nS\t6\tGATTACA\nS\t7\tAC\n'
    printf 'L\t1\t+\t2\t+\t0M\nL\t1\t+\t3\t+\t0M\nL\t2\t+\t4\t+\t0M\nL\t3\t+\t4\t+\t0M\n'
    printf 'L\t4\t+\t5\t+\t0M\nL\t5\t+\t6\t+\t0M\nL\t4\t+\t6\t+\t0M\n'
    printf 'W\tCHM13\t0\tchrA\t0\t18\t>1>2>4>6\nW\tHG1\t1\tchrA\t0\t21\t>1>3>4>5>6\n'
    printf 'W\tCHM13\t0\tchrB\t0\t2\t>7\n'; } > tiny.gfa
  computed_pair "catalog: computed snarls, a tiny graph" tiny.gfa chrA chrB
fi
# dromos snarls --dist --graph: the graph's own index read, others refused, no table written.
$D snarls --dist $SD --graph $GBZ_SMALL -o sx.tsv 2> sx.err
label="snarls --dist --graph: the graph's own distance index"; ok test -s sx.tsv
rm -f sx.tsv
if [ -n "${DIST:-}" ]; then
  fails "snarls --dist --graph: another graph's distance index (HPRC's .dist, the small graph)" $D snarls --dist $DIST \
    --graph $GBZ_SMALL -o sx.tsv
  label="snarls --dist --graph: another graph's refused, named so"; ok grep -q "not this graph's distance index" f.err
fi
fails "snarls --dist --graph: another graph's distance index (the small graph's, HPRC)" $D snarls --dist $SD $DN \
  -t 4 -o sx.tsv
label="snarls --dist --graph: another graph's refused, named so (HPRC)"
ok grep -q "not this graph's distance index" f.err
label="snarls --dist --graph: no table when refused"; ok bash -c "! test -e sx.tsv"
head -c 50000 $SD > cut.dist; fails "snarls --dist --graph: a distance index cut short" $D snarls --dist cut.dist \
  --graph $GBZ_SMALL -o sx.tsv
python3 -I - $SD v5.dist << 'PY'
import struct, sys
d = bytearray(open(sys.argv[1], 'rb').read())
obj = 36 + struct.unpack_from('<q', d, 36)[0]
width = struct.unpack_from('<Q', d, obj + 8)[0]
data = obj + 40 + struct.unpack_from('<q', d, obj + 40)[0]
bit = width  # entry 1 holds 1023 ^ version: version 4 to 5 clears its lowest bit
d[data + bit // 8] &= ~(1 << (bit % 8)) & 0xFF
open(sys.argv[2], 'wb').write(d)
PY
fails "snarls --dist --graph: a distance index of version 5" $D snarls --dist v5.dist --graph $GBZ_SMALL -o sx.tsv
label="snarls --dist: version 5 refused, naming the way out"; ok grep -q "dromos snarls --graph computes" f.err
fails "snarls: --dist with --vg-file" $D snarls --dist $SD --vg-file $S/chr6.snarls -o sx.tsv
fails "snarls: --dist with --stats" $D snarls --dist $SD --stats -o sx.tsv
# The catalog's views, against pgcatalog's catalog of chr6 (and of a region past chr6's end: no records).
rm -rf cb; mkdir -p cb; cp -r $S cb/snarls
$PGCATALOG build --graph $GBZ_SMALL --reference-sample CHM13 --out cb --region chr6 --region chr6:200000000-200000100 \
  > cb.err 2>&1
$D catalog --graph $GBZ_SMALL --reference-sample CHM13 --snarls-cache $S --out cv.sites --region chr6 -t 2 2> /dev/null
$D catalog --graph $GBZ_SMALL --reference-sample CHM13 --snarls-cache $S --out ce.sites \
  --region chr6:200000000-200000100 -t 2 2> /dev/null
V=cv.sites; PV=cb/chr6.pgsites
view_pair() {  # name, region...: dromos view prints what pgcatalog view prints of chr6's catalog, or both fail
  local name=$1; shift
  $D view $V "$@" > vd.out 2> vd.err; local drc=$?
  $PGCATALOG view $PV "$@" > vp.out 2> vp.err; local prc=$?
  if [ $drc -ne 0 ] || [ $prc -ne 0 ]; then
    [ $drc -ne 0 ] && [ $prc -ne 0 ] && echo "ok    $name: both fail ($(tail -1 vd.err))" \
      || echo "FAIL  $name: dromos rc=$drc, pgcatalog rc=$prc"; return
  fi
  cmp -s vd.out vp.out && echo "ok    $name" || echo "FAIL  $name"
}
view_pair "view: whole"
view_pair "view: a region" chr6:31,500,000-31,520,000
view_pair "view: a region by sample" CHM13#0#chr6:31500000-31600000
view_pair "view: the contig alone" chr6
view_pair "view: past the end" chr6:900000000-900000001
view_pair "view: one base" chr6:31500001-31500001
view_pair "view: a region without an end" chr6:5
view_pair "view: a region with a bad start" chr6:x-5
label="view: a catalog of no records"
ok cmp -s <($D view ce.sites chr6) <($PGCATALOG view cb/chr6_200000000-200000100.pgsites)
fails "view: another contig" $D view $V chr19:1-100000000
fails "view: another sample's contig" $D view $V GRCh38#0#chr6
fails "view: not a catalog" $D view $SMALL_GAF
fails "view: a missing catalog" $D view /nonexistent.sites
fails "view: a .pgsites" $D view $PV
label="view: a .pgsites, named"; ok grep -q "pgsites" f.err
head -c $(( $(stat -c %s $V) / 2 )) $V > cut.sites; fails "view: a catalog cut short" $D view cut.sites
label="view: a catalog cut short, named"; ok grep -q "truncated" f.err
python3 -I -c "import sys; d=bytearray(open(sys.argv[1],'rb').read()); d[200]^=0xFF; open(sys.argv[2],'wb').write(d)" \
  $V flip.sites
fails "view: a damaged block" $D view flip.sites
python3 -I -c "import sys; d=bytearray(open(sys.argv[1],'rb').read()); d[8]=2; open(sys.argv[2],'wb').write(d)" \
  $V v2.sites
fails "view: another format version" $D view v2.sites
label="view: another format version, named"; ok grep -q "version 2" f.err
fails "catalog-info: a damaged catalog" $D catalog-info cut.sites
# surject refuses a catalog of another graph, or one lacking a contig.
fails "surject: another graph's catalog" $D surject $DG --contig chr6 --catalog $V $G
$D catalog $DG --region chrY:1-2000000 --out y.sites -t 4 2> /dev/null
fails "surject: a catalog lacking the contig" $D surject $DG --contig chr20 --catalog y.sites $G
label="surject: a catalog of a region of its contig"; ok bash -c "$D surject $DG --contig chrY --catalog y.sites $G \
  > /dev/null 2> y.err"
