#!/bin/bash
#SBATCH -p cpu
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH --qos=ngs16
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#
# THE QUESTION: how much of the ~6000-event DCM-vs-NF signature survives
# adjustment for cell composition?
#
# Why it matters: DCM hearts have less myocardium and more fibrosis than donor
# hearts. Your top genome-wide hits are DCN (fibroblast/ECM), FABP3, PLA2G5 --
# cell-type markers. If most events vanish on composition adjustment, the
# signature is largely tissue makeup and the paper's central claim changes.
# A reviewer WILL ask this. Better to know first.
#
# Runs three models per event type and compares them:
#   base      PSI ~ disease + sex + age + RIN                      (already have)
#   compadj   + CM_score + FB_score                                (the test)
#   fulladj   + CM_score + FB_score + log10_libsize + TIN          (everything)
#
# Needs magnet_expr.per_sample.tsv WITH the TIN column -- 03b added it.
# Submit:  sbatch 08_composition_adjusted_discovery.sbatch
# ~50 min.

date
export PATH=/data/apps/anaconda3/envs/mamba/envs/rnaseq_tools/bin:$PATH
ROOT=/data/project/chureeratp/RNA/HEART
cd "$ROOT"

EXPR=magnet_expr.per_sample.tsv
[ -s "$EXPR" ] || { echo "missing $EXPR -- run 02 first"; exit 1; }
echo "per-sample columns: $(head -1 $EXPR)"
head -1 "$EXPR" | grep -q "TIN" || echo "WARNING: no TIN column; fulladj will fail. Run 03b first."

for ET in SE MXE; do
  MATS=rmats_magnet_rl125/${ET}.MATS.JCEC.txt
  [ -s "$MATS" ] || { echo "missing $MATS -- skipping"; continue; }

  echo
  echo "################ $ET : + composition ################"
  python phase2_diffsplice_v2.py \
    --mats "$MATS" --event-type $ET \
    --b1 dcm_bams.txt --b2 nf_bams.txt \
    --cov magnet_covariates.tsv --min-samples 20 \
    --extra-cov-file "$EXPR" --extra-cov CM_score,FB_score \
    --out ${ET}_v2compadj_rl125

  echo
  echo "################ $ET : + composition + technical ################"
  python phase2_diffsplice_v2.py \
    --mats "$MATS" --event-type $ET \
    --b1 dcm_bams.txt --b2 nf_bams.txt \
    --cov magnet_covariates.tsv --min-samples 20 \
    --extra-cov-file "$EXPR" --extra-cov CM_score,FB_score,log10_libsize,TIN \
    --out ${ET}_v2fulladj_rl125
done

echo
echo "################ CONCORDANCE ################"
for ET in SE MXE; do
  BASE=${ET}_v2filtered_rl125.all_events.tsv
  [ -s "$BASE" ] || { echo "missing $BASE -- run 04 first"; continue; }
  echo
  echo "===== $ET ====="
  python 06_sensitivity_concordance.py \
    --primary "$BASE" \
    --compare composition=${ET}_v2compadj_rl125.all_events.tsv \
              full=${ET}_v2fulladj_rl125.all_events.tsv \
    --out ${ET}_compositionCheck_rl125
done

echo
echo "################ WHAT TO READ ################"
echo "In each concordance block, the line that matters is:"
echo "  'primary-significant events retained: X% by full key'"
echo "  X > 70%  -> the signature is robust; composition is a caveat, not the story"
echo "  X 40-70% -> composition explains a substantial minority; report both models"
echo "  X < 40%  -> the signature is largely tissue makeup; the paper changes again"
echo "Also check 'RBM20-program events lost' -- if the program genes drop out"
echo "preferentially, that is a finding in itself."
date
