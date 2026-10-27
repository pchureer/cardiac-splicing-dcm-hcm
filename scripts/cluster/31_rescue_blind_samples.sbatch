#!/bin/bash
#SBATCH -p cpu
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH --qos=ngs16
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#
# 31 -- rescue the samples the identity check could not see.
#
# THE SITUATION
# -------------
# Script 30 found 10 MAGNet samples with zero called genotypes. Their bqsr BAMs exist
# and are indexed but are EMPTY (0 MB), so bcftools mpileup read nothing from them and
# reported no error. They were never compared to any other sample: neither cleared of
# duplication nor implicated.
#
# WHY THIS IS PROBABLY NOT A SPLICING PROBLEM
# -------------------------------------------
# rMATS runs on the RAW STAR bam (*Aligned.sortedByCoord.out.bam), not on bqsr. bqsr
# output is SplitNCigar-processed and is only used for variant calling. So an empty
# bqsr bam breaks genotyping while leaving the splicing analysis untouched -- PROVIDED
# the raw STAR bam is healthy. Step 1 establishes which of the two is true, and that
# distinction decides whether this is a five-minute fix or a real hole in the discovery
# cohort. Do not skip reading step 1's table.
#
# WHAT IT DOES
#   1) report both bam types, with exact byte sizes, for every zero-call sample, and
#      say whether each sample is actually used by phase 2
#   2) genotype them from the best available bam at the SAME panel (geno/panel.tsv)
#   3) splice the result into geno/merged.vcf.gz, replacing the empty columns
#   4) tell you to re-run 28 and 29
#
# Nothing is overwritten: the old matrix is kept as geno/merged.preRescue.vcf.gz.
#
# Run:  bash 31_rescue_blind_samples.sbatch        (a few minutes, 10 samples)
date
export PATH=/data/apps/anaconda3/envs/mamba/envs/rnaseq_tools/bin:$PATH
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
ROOT=/data/project/chureeratp/RNA/HEART
cd "$ROOT"

MINBAM=${MINBAM:-1000000}     # a bam smaller than this is treated as unusable
command -v bcftools >/dev/null || { echo "bcftools not on PATH"; exit 1; }
for f in geno/panel.tsv geno/samples.tsv geno/sample_qc.tsv geno/merged.vcf.gz; do
  [ -s "$f" ] || { echo "*** missing $f -- run 27, 28 and 29 first"; exit 1; }
done
FA=$(ls /data/project/chureeratp/hg38/ref/*.fa /data/project/chureeratp/hg38/ref/*.fasta \
     2>/dev/null | head -1)
[ -n "$FA" ] || { echo "*** no reference FASTA"; exit 1; }
NS=$(awk 'NF' geno/panel.tsv | wc -l)
echo "panel: $NS sites   reference: $FA"

echo
echo "############################################################"
echo "#####  1) which bam is broken, and is the sample even used?"
echo "############################################################"
python3 - "$MINBAM" <<'PY' || exit 1
import glob, os, re, sys
ROOT = "/data/project/chureeratp/RNA/HEART"
MINBAM = int(sys.argv[1])
ACC = re.compile(r"(SRR\d+|ERR\d+|SAMEA\d+|SAMN\d+)")

qc = [l.split("\t") for l in open("geno/sample_qc.tsv").read().splitlines()]
hdr, rows = qc[0], qc[1:]
si, ni = hdr.index("sample"), hdr.index("n_called")
blind = [r[si] for r in rows if int(float(r[ni])) == 0]
if not blind:
    print("[info] no zero-call samples -- nothing to rescue")
    open("geno/rescue.tsv", "w").write("sample\tbam\n")
    raise SystemExit(0)

# is the sample referenced by any phase 2 bam list?
used = {}
for f in glob.glob(os.path.join(ROOT, "*_bams.txt")):
    if "_dedup" in f:
        continue
    raw = open(f).read()
    for p in re.split(r"[,\n]", raw):
        p = p.strip()
        if not p:
            continue
        m = ACC.search(p)
        sid = m.group(1) if m else os.path.basename(os.path.dirname(p))
        used.setdefault(sid, []).append(os.path.basename(f))

print(f"  {len(blind)} zero-call sample(s)\n")
print(f"    {'sample':<15}{'bqsr bam':>13}{'raw STAR bam':>15}  {'use':<7}{'phase2 lists'}")
pick, broken_both = [], []
for s in blind:
    d = os.path.join(ROOT, s)
    bq = sorted(glob.glob(os.path.join(d, "*_starAligned.bqsr.bam")))
    st = sorted(glob.glob(os.path.join(d, "*Aligned.sortedByCoord.out.bam")))
    sz = lambda L: os.path.getsize(L[0]) if L and os.path.exists(L[0]) else -1
    zb, zs = sz(bq), sz(st)
    fb = f"{zb:,}" if zb >= 0 else "absent"
    fs = f"{zs:,}" if zs >= 0 else "absent"
    if zs >= MINBAM:
        chosen, tag = st[0], "STAR"
    elif zb >= MINBAM:
        chosen, tag = bq[0], "bqsr"
    else:
        chosen, tag = None, "NONE"
        broken_both.append(s)
    lists = ",".join(sorted(set(used.get(s, [])))) or "-- not used --"
    print(f"    {s:<15}{fb:>13}{fs:>15}  {tag:<7}{lists}")
    if chosen:
        pick.append((s, chosen))

with open("geno/rescue.tsv", "w") as fh:
    fh.write("sample\tbam\n")
    for s, b in pick:
        fh.write(f"{s}\t{b}\n")
print(f"\n  rescuable from a healthy bam : {len(pick)}")
print(f"  both bams unusable           : {len(broken_both)}")
if broken_both:
    print("\n  *** These samples have NO usable alignment at all:")
    for s in broken_both:
        print(f"      {s}  ({','.join(sorted(set(used.get(s, [])))) or 'not in any bam list'})")
    print("  *** If any is listed in a phase 2 bam list, rMATS read nothing from it")
    print("  *** either, so that contrast's true n is lower than the .ids file says.")
    print("  *** Check the rMATS summary for those runs before trusting any n.")
print("\n  [write] geno/rescue.tsv")
PY

NR=$(( $(wc -l < geno/rescue.tsv) - 1 ))
if [ "$NR" -lt 1 ]; then
  echo
  echo "Nothing rescuable. Report the identity check as covering 660 of 670 samples,"
  echo "and resolve the alignment failure separately."
  exit 0
fi

echo
echo "############################################################"
echo "#####  2) genotype the rescued samples at the same panel"
echo "############################################################"
awk 'NR>1 {print $2}' geno/rescue.tsv > geno/rescue_bams.txt
awk 'NR>1 {print $1}' geno/rescue.tsv > geno/rescue_names.txt
echo "genotyping $NR sample(s) at $NS sites"
bcftools mpileup -f "$FA" -R geno/panel.tsv -b geno/rescue_bams.txt \
    -q 10 -Q 15 -I -a FORMAT/DP --ignore-RG 2>geno/rescue.err \
  | bcftools call -m -Oz -o geno/rescue.raw.vcf.gz 2>>geno/rescue.err
[ -s geno/rescue.raw.vcf.gz ] || { echo "*** mpileup failed:"; tail -10 geno/rescue.err; exit 1; }
bcftools reheader -s geno/rescue_names.txt -o geno/rescue.vcf.gz geno/rescue.raw.vcf.gz
bcftools index -f -t geno/rescue.vcf.gz
rm -f geno/rescue.raw.vcf.gz
echo "rescued VCF: $(bcftools view -H geno/rescue.vcf.gz | wc -l) sites"

echo
echo "############################################################"
echo "#####  3) splice into the merged matrix"
echo "############################################################"
# The empty columns are already in merged.vcf.gz, so they must be REMOVED before the
# rescued ones are merged in -- bcftools merge refuses duplicate sample names, and
# silently keeping both would be worse.
cp -n geno/merged.vcf.gz geno/merged.preRescue.vcf.gz 2>/dev/null || true
DROP=$(paste -sd, geno/rescue_names.txt)
bcftools view -s "^${DROP}" -Oz -o geno/merged.trimmed.vcf.gz geno/merged.vcf.gz \
  || { echo "*** could not drop the empty columns"; exit 1; }
bcftools index -f -t geno/merged.trimmed.vcf.gz
bcftools merge -m none -Oz -o geno/merged.rescued.vcf.gz \
    geno/merged.trimmed.vcf.gz geno/rescue.vcf.gz 2>geno/rescue_merge.err \
  || { echo "*** merge failed:"; tail -5 geno/rescue_merge.err; exit 1; }
bcftools index -f -t geno/merged.rescued.vcf.gz
NSAMP=$(bcftools view -h geno/merged.rescued.vcf.gz | grep -m1 '^#CHROM' | awk '{print NF-9}')
echo "rescued matrix: $NSAMP samples, $(bcftools view -H geno/merged.rescued.vcf.gz | wc -l) sites"

mv -f geno/merged.rescued.vcf.gz geno/merged.vcf.gz
mv -f geno/merged.rescued.vcf.gz.tbi geno/merged.vcf.gz.tbi
rm -f geno/merged.trimmed.vcf.gz geno/merged.trimmed.vcf.gz.tbi
rm -f geno/gt.txt geno/gt.samples        # force 28 to re-extract from the new matrix
echo "old matrix kept as geno/merged.preRescue.vcf.gz"

echo
echo "############################################################"
echo "#####  4) re-run the verdict"
echo "############################################################"
python3 28_identity_metrics.py > step5_rescued.txt 2>&1 \
  || python 28_identity_metrics.py > step5_rescued.txt 2>&1
python3 29_identity_triage.py > triage_rescued.txt 2>&1 \
  || python 29_identity_triage.py > triage_rescued.txt 2>&1
echo
grep -E "informative sites|call rate|CONFIRMED|CANDIDATE|ARTEFACT|zero called" \
  step5_rescued.txt triage_rescued.txt 2>/dev/null | head -20
echo
echo "Compare against the first run:"
echo "  * CONFIRMED should still be 26 -- the same 23 MAGNet/LBNL pairs plus the 3"
echo "    within-cohort replicates. A DIFFERENT number means one of the rescued"
echo "    samples is itself a duplicate, which is exactly what could not be seen before."
echo "  * section 1 of triage_rescued.txt should now report 0 zero-call samples."
echo "  * if it does, the methods sentence becomes 670 of 670, not 660."
date
