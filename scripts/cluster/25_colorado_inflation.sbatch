#!/bin/bash
#SBATCH -p cpu
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH --qos=ngs16
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#
# WHY COLORADO CALLS SO MANY EVENTS SIGNIFICANT -- and whether it is real.
#
# THE OBSERVATION (from master_numbers.tsv, exon-level)
#   cohort    contrast      n        SE      MXE
#   MAGNet    DCM vs NF   166/166  1.47%    0.37%
#   Colorado  DCM vs NF    37/14   2.23%   11.87%
#   Colorado  ICM vs NF    13/14   2.37%   21.38%
#
# A 27-sample contrast calling 21% of MXE events significant, where a 332-sample
# cohort calls 0.37%, is not a biological result.
#
# MIN-SAMPLES IS THE WRONG SUSPECT. Colorado's --min-samples 10 against groups of
# 13 and 14 demands coverage in 77% and 71% of samples; MAGNet's --min-samples 20
# against 166/166 demands 12%. Colorado is STRICTER in fractional terms. Raising it
# further would only shrink the event set, not fix a false-positive rate.
#
# THE LIKELY CAUSE: Colorado runs --no-covariates while MAGNet adjusts for RIN, age,
# sex and pool. If sequencing depth or RNA degradation differs systematically between
# the DCM, ICM and NF groups -- plausible for a cohort prepared in batches -- that
# difference lands entirely in the disease coefficient. MXE suffers most because it is
# the noisiest event class and the hardest to quantify from single-end 51 bp reads.
#
# HOW THIS SCRIPT SETTLES IT
#   A) builds a technical covariate from the rMATS table itself -- per-sample median
#      junction coverage, log10. No BAM access, no samtools, no new alignment. This is
#      a depth/degradation proxy specific to the splicing measurement.
#   B) re-runs each contrast with and without it, at two coverage stringencies.
#   C) PERMUTATION NULL: shuffles the group labels and counts "significant" events.
#      This is the decisive number. A well-calibrated analysis returns ~0. Whatever it
#      returns IS the false-positive count, and no argument can talk it down.
#   D) runs the same permutation on MAGNet SE as the calibration reference, so the
#      comparison is against this project's own well-behaved contrast rather than an
#      assumption.
#
# WHAT THIS DOES *NOT* THREATEN. The etiology-independence claim is a WITHIN-Colorado
# comparison: the same 37 DCM hearts, the same pipeline, the same settings, thousands of
# events against NF and zero against ICM. An inflated filter applies to both arms, so
# the contrast stays fair and the positive control holds. What is at risk is only the
# absolute counts -- which must not be printed in a table beside MAGNet's until this
# runs.
#
# Submit: sbatch 25_colorado_inflation.sbatch
# RUNNING WITHOUT SLURM
#   The Colorado sections (A-C) are small -- 25k SE and 14k MXE events over ~51
#   samples -- and are fine on a login node. Section D is NOT: MAGNet is 100k events
#   over 332 samples, twice. Skip it there and submit it separately when CPUs free up:
#
#     cd /data/project/chureeratp/RNA/HEART
#     SKIP_MAGNET=1 nohup bash 25_colorado_inflation.sbatch > inflate_run.log 2>&1 &
#     tail -f inflate_run.log
#
#   nohup + & means it survives losing the connection. Better still, if any node has a
#   spare core:  srun -p cpu -c 4 --mem=32G --pty bash 25_colorado_inflation.sbatch
date
export PATH=/data/apps/anaconda3/envs/mamba/envs/rnaseq_tools/bin:$PATH
# one BLAS thread per process: phase 2 is not parallel, and letting numpy grab every
# core is what gets a login-node job killed
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
SKIP_MAGNET=${SKIP_MAGNET:-0}
ROOT=/data/project/chureeratp/RNA/HEART
cd "$ROOT"
mkdir -p inflate

for f in rmats_Colorado/SE.MATS.JCEC.txt Colorado_dcm_bams.txt Colorado_nf_bams.txt \
         colorado_dcm_nf.cov; do
  [ -s "$f" ] || { echo "missing $f"; exit 1; }
done

echo "############################################################"
echo "#####  A) technical covariate from the rMATS table"
echo "############################################################"
python3 - <<'PY'
import numpy as np, pandas as pd, re, os

def ids(path):
    txt = open(path).read().strip()
    out = []
    for p in re.split(r"[,\n]", txt):
        p = p.strip()
        if not p:
            continue
        m = re.search(r"(SRR\d+|ERR\d+|SAMEA\d+|SAMN\d+)", p)
        out.append(m.group(1) if m else p)
    return out

def jcov(mats, order, nmax=40000):
    """per-sample median (IJC+SJC) across events -- depth/degradation proxy"""
    t = pd.read_csv(mats, sep="\t", nrows=nmax)
    I = t["IJC_SAMPLE_1"].astype(str) + "," + t["IJC_SAMPLE_2"].astype(str)
    S = t["SJC_SAMPLE_1"].astype(str) + "," + t["SJC_SAMPLE_2"].astype(str)
    def mat(col):
        return np.array([[float(x) if x not in ("", "NA") else np.nan
                          for x in row.split(",")] for row in col])
    tot = mat(I) + mat(S)
    if tot.shape[1] != len(order):
        raise SystemExit(f"[fatal] {mats}: {tot.shape[1]} sample columns but "
                         f"{len(order)} ids from the bam lists")
    med = np.nanmedian(tot, axis=0)
    return pd.DataFrame({"sample": order,
                         "med_jcov": med,
                         "log10_jcov": np.log10(np.maximum(med, 1.0))})

JOBS = [
    ("colorado_dcmnf", "rmats_Colorado/SE.MATS.JCEC.txt",
     "Colorado_dcm_bams.txt", "Colorado_nf_bams.txt"),
    ("colorado_icmnf", "rmats_Colorado_ICMvsNF/SE.MATS.JCEC.txt",
     "Colorado_icm_bams.txt", "Colorado_nf_bams.txt"),
]
for tag, mats, b1, b2 in JOBS:
    if not (os.path.exists(mats) and os.path.exists(b1) and os.path.exists(b2)):
        print(f"[skip] {tag}: missing input")
        continue
    order = ids(b1) + ids(b2)
    d = jcov(mats, order)
    d.to_csv(f"inflate/{tag}_tech.tsv", sep="\t", index=False, float_format="%.5f")
    n1 = len(ids(b1))
    g1, g2 = d.log10_jcov.values[:n1], d.log10_jcov.values[n1:]
    from scipy import stats
    t, p = stats.ttest_ind(g1, g2, equal_var=False)
    print(f"[write] inflate/{tag}_tech.tsv   n={len(d)}")
    print(f"  median junction coverage: group1 {10**np.mean(g1):.0f} "
          f"vs group2 {10**np.mean(g2):.0f}   Welch p = {p:.3g}")
    if p < 0.05:
        print("  *** DEPTH DIFFERS BY GROUP. Without this covariate the difference goes")
        print("  *** straight into the disease coefficient. This is very likely the cause.")
    else:
        print("  (no significant depth difference -- if the inflation persists after")
        print("   adjustment, look at degradation or at MXE quantification instead)")
PY

echo
echo "############################################################"
echo "#####  B+C) settings sweep and permutation null"
echo "############################################################"

# $1 tag  $2 mats  $3 b1  $4 b2  $5 cov  $6 g1  $7 g2  $8 min-samples
#   $9 extra ("" or "tech")  $10 min-reads  $11 min-psi-sd  $12 permseed (0 = none)
run () {
  local tag=$1 mats=$2 b1=$3 b2=$4 cov=$5 g1=$6 g2=$7 ms=$8 ex=$9
  local mr=${10} sd=${11} seed=${12} et=${13}
  local out=inflate/${tag}
  local COVFILE=$cov
  if [ "$seed" != "0" ]; then
    COVFILE=inflate/${tag}.permcov
    python3 - "$cov" "$COVFILE" "$seed" "$g1" "$g2" <<'PY'
import random, sys
src, dst, seed, g1, g2 = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4], sys.argv[5]
rows = [l.rstrip("\n").split("\t") for l in open(src) if l.strip()]
labs = [r[1] for r in rows]
random.Random(seed).shuffle(labs)          # break the disease-label association
for r, l in zip(rows, labs):
    r[1] = l
open(dst, "w").write("\n".join("\t".join(r) for r in rows) + "\n")
PY
  fi
  local XTRA=""
  if [ "$ex" = "tech" ]; then
    local TF=inflate/$(echo $tag | sed 's/_.*//')_tech.tsv
    # tag is like colorado_dcmnf_S1 -> tech file is colorado_dcmnf_tech.tsv
    TF=inflate/$(echo $tag | cut -d_ -f1-2)_tech.tsv
    [ -s "$TF" ] || { echo "  missing $TF -- skipping $tag"; return; }
    XTRA="--extra-cov-file $TF --extra-cov log10_jcov"
  fi
  python phase2_diffsplice_v2.py \
    --mats "$mats" --event-type $et \
    --b1 "$b1" --b2 "$b2" --cov "$COVFILE" --group1 $g1 --group2 $g2 \
    --min-samples $ms --min-reads $mr --min-psi-sd $sd \
    --no-covariates $XTRA \
    --out ${out} > ${out}.log 2>&1
  # VERIFY the covariate actually entered the model rather than being silently dropped
  local model=$(grep -m1 "^model:" ${out}.log)
  if [ "$ex" = "tech" ] && ! echo "$model" | grep -q "log10_jcov"; then
    echo "  *** $tag: log10_jcov NOT in the model -- '$model'"
    echo "  *** --no-covariates may be stripping --extra-cov. STOP and tell me."
  fi
  local n=$(grep -m1 "events tested:" ${out}.log | grep -o '[0-9]*')
  local s=$(python3 -c "
import pandas as pd,sys
try:
    t=pd.read_csv('${out}.all_events.tsv',sep='\t').drop_duplicates('key')
    print(int(((t.FDR<0.05)&(t.dPSI_adj.abs()>=0.10)).sum()))
except Exception: print('NA')
" 2>/dev/null)
  printf "  %-34s %-5s tested=%-8s sig=%-7s  %s\n" "$tag" "$et" "${n:-NA}" "${s:-NA}" "$model"
}

for ET in SE MXE; do
  echo
  echo "===== $ET : Colorado DCM vs NF ====="
  M=rmats_Colorado/${ET}.MATS.JCEC.txt
  [ -s "$M" ] || { echo "missing $M"; continue; }
  run colorado_dcmnf_S0_$ET  "$M" Colorado_dcm_bams.txt Colorado_nf_bams.txt colorado_dcm_nf.cov DCM NF 10 ""     10 0.02 0 $ET
  run colorado_dcmnf_S1_$ET  "$M" Colorado_dcm_bams.txt Colorado_nf_bams.txt colorado_dcm_nf.cov DCM NF 10 "tech" 10 0.02 0 $ET
  run colorado_dcmnf_S2_$ET  "$M" Colorado_dcm_bams.txt Colorado_nf_bams.txt colorado_dcm_nf.cov DCM NF 12 "tech" 20 0.05 0 $ET
  for s in 101 102 103; do
    run colorado_dcmnf_P0s${s}_$ET "$M" Colorado_dcm_bams.txt Colorado_nf_bams.txt colorado_dcm_nf.cov DCM NF 10 ""     10 0.02 $s $ET
    run colorado_dcmnf_P2s${s}_$ET "$M" Colorado_dcm_bams.txt Colorado_nf_bams.txt colorado_dcm_nf.cov DCM NF 12 "tech" 20 0.05 $s $ET
  done

  echo
  echo "===== $ET : Colorado ICM vs NF ====="
  M2=rmats_Colorado_ICMvsNF/${ET}.MATS.JCEC.txt
  if [ -s "$M2" ]; then
    run colorado_icmnf_S0_$ET "$M2" Colorado_icm_bams.txt Colorado_nf_bams.txt colorado_icm_nf.cov ICM NF 10 ""     10 0.02 0 $ET
    run colorado_icmnf_S2_$ET "$M2" Colorado_icm_bams.txt Colorado_nf_bams.txt colorado_icm_nf.cov ICM NF 9  "tech" 20 0.05 0 $ET
    for s in 101 102; do
      run colorado_icmnf_P0s${s}_$ET "$M2" Colorado_icm_bams.txt Colorado_nf_bams.txt colorado_icm_nf.cov ICM NF 10 ""     10 0.02 $s $ET
      run colorado_icmnf_P2s${s}_$ET "$M2" Colorado_icm_bams.txt Colorado_nf_bams.txt colorado_icm_nf.cov ICM NF 9  "tech" 20 0.05 $s $ET
    done
  fi
done

echo
echo "############################################################"
echo "#####  D) MAGNet SE permutation -- the calibration reference"
echo "############################################################"
if [ "$SKIP_MAGNET" = "1" ]; then
  echo "SKIPPED (SKIP_MAGNET=1). Sections A-C stand on their own: Colorado's own"
  echo "permutation null is the false-positive count, with or without the MAGNet"
  echo "reference. Run section D later with:"
  echo "    sbatch 25_colorado_inflation.sbatch     # full, when CPUs are free"
else
echo "A well-behaved contrast should return ~0 significant events under permutation."
for s in 101 102; do
  python3 - magnet_covariates.tsv inflate/magnet.P${s}.cov $s <<'PY'
import random, sys
src, dst, seed = sys.argv[1], sys.argv[2], int(sys.argv[3])
rows = [l.rstrip("\n").split("\t") for l in open(src) if l.strip()]
labs = [r[1] for r in rows]
random.Random(seed).shuffle(labs)
for r, l in zip(rows, labs):
    r[1] = l
open(dst, "w").write("\n".join("\t".join(r) for r in rows) + "\n")
PY
  python phase2_diffsplice_v2.py \
    --mats rmats_magnet/SE.MATS.JCEC.txt --event-type SE \
    --b1 dcm_bams.txt --b2 nf_bams.txt --cov inflate/magnet.P${s}.cov \
    --min-samples 20 --out inflate/magnet_P${s}_SE > inflate/magnet_P${s}_SE.log 2>&1
  n=$(grep -m1 "events tested:" inflate/magnet_P${s}_SE.log | grep -o '[0-9]*')
  sg=$(python3 -c "
import pandas as pd
t=pd.read_csv('inflate/magnet_P${s}_SE.all_events.tsv',sep='\t').drop_duplicates('key')
print(int(((t.FDR<0.05)&(t.dPSI_adj.abs()>=0.10)).sum()))" 2>/dev/null)
  printf "  MAGNet SE permutation seed %-5s tested=%-8s sig=%s\n" "$s" "${n:-NA}" "${sg:-NA}"
done
fi

echo
echo "############################################################"
echo "HOW TO READ IT"
echo "############################################################"
echo "S0 = current settings (reproduces master_numbers.tsv)"
echo "S1 = + log10_jcov technical covariate"
echo "S2 = + log10_jcov, min-samples up, min-reads 20, PSI sd 0.05"
echo "P0/P2 = the SAME settings with group labels SHUFFLED"
echo
echo "The permutation rows are the answer. sig(P) is the false-positive count at"
echo "FDR<0.05 and |dPSI|>=0.10 with no real signal present."
echo
echo "  sig(P0) near 0 and sig(S0) large   -> the signal is real and the counts stand."
echo "  sig(P0) large                      -> S0 is uncalibrated. Report the setting"
echo "                                        whose permutation returns ~0, and give"
echo "                                        the permutation counts in the supplement."
echo "  sig(P2) near 0 while sig(P0) large -> the covariate/stringency fix WORKED."
echo "                                        Re-run master_numbers with S2 for Colorado"
echo "                                        and note the change in the methods."
echo
echo "Compare against the MAGNet permutation. If MAGNet returns ~0 and Colorado S0"
echo "returns hundreds, that contrast is the cleanest possible statement of the problem"
echo "and belongs in the supplement exactly as measured."
echo
echo "REMEMBER: the etiology-independence claim is unaffected either way -- it is a"
echo "within-Colorado comparison at identical settings. Only the absolute counts move."
date
