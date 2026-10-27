#!/bin/bash
#SBATCH -p cpu
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH --qos=ngs16
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#
# RBM20 expression <-> target-exon PSI, SE and MXE.
# Needs magnet_expr.per_sample.tsv from 02 first.
# Submit:  sbatch 03_rbm20_psi.sbatch

date
export PATH=/data/apps/anaconda3/envs/mamba/envs/rnaseq_tools/bin:$PATH
ROOT=/data/project/chureeratp/RNA/HEART
cd "$ROOT"

if [ ! -s magnet_expr.per_sample.tsv ]; then
  echo "magnet_expr.per_sample.tsv missing -- run 02_build_expr_matrix.sbatch first"
  exit 1
fi

# --------------------------------------------------------------------------
# Read length. RLSUF drives the rMATS directory and the output prefix together,
# so a run is consistent at ONE read length or it does not run. Empty = the
# original 50 bp behaviour, reproducing the 7 Sep output.
#
# Called RLSUF, not SUF: 56_rbm20_targets45.sbatch already uses SUF for the
# model variant ("" | "_techadj"), and a collision there would silently
# mis-name every input and output.
RLSUF="${TBL_SUFFIX:-}"
RMATS_DIR="${RMATS_DIR:-rmats_magnet${RLSUF}}"
# --------------------------------------------------------------------------

# --------------------------------------------------------------------------
# A wrong rMATS directory must be an ERROR, not a skip: the per-class guard
# below only `continue`s, which turns a mistyped directory into a job that
# completes having done nothing and reports no failure.
if [ ! -d "$RMATS_DIR" ]; then
  echo "*** rMATS directory not found: $RMATS_DIR"
  echo "*** TBL_SUFFIX='$RLSUF'. Set RMATS_DIR explicitly if the directory for"
  echo "*** this read length has another name."
  exit 1
fi

# DON'T-CLOBBER GUARD (03, unadjusted).
# The 7 Sep tables are the 50 bp comparators, and 56's *_g45 / *_repro9 reports
# were built from them. Overwriting them orphans those too, and nothing errors.
for ET in SE MXE; do
  for EXT in events.tsv report.txt; do
    F="rbm20_${ET}${RLSUF}.${EXT}"
    if [ -e "$F" ] && [ -z "${OVERWRITE_OK:-}" ]; then
      echo "*** REFUSING — $F already exists"
      echo "***   modified: $(date -r "$F" '+%Y-%m-%d %H:%M' 2>/dev/null)"
      echo "***   size:     $(stat -c %s "$F" 2>/dev/null) bytes"
      echo "***"
      echo "*** Set OVERWRITE_OK=1 only if you mean to replace it."
      exit 1
    fi
  done
done
echo "[guard] 03, unadjusted: no existing output for TBL_SUFFIX='${RLSUF}' — safe to write"
# --------------------------------------------------------------------------
# JC (junction-only), NOT JCEC: JCEC folds exon-body reads into PSI and body
# coverage scales with expression, which is the axis we are correlating against.
for ET in SE MXE; do
  echo "=================== $ET ==================="
  MATS=${RMATS_DIR}/${ET}.MATS.JC.txt
  if [ ! -s "$MATS" ]; then echo "missing $MATS -- skipping"; continue; fi
  python 03_rbm20_psi_correlation.py \
    --rmats "$MATS" --event-type $ET \
    --b1 dcm_bams.txt --b2 nf_bams.txt \
    --expr magnet_expr.per_sample.tsv \
    --covariates magnet_covariates.tsv \
    --out-prefix rbm20_${ET}${RLSUF} \
    --n-perm 200
done

echo "--- reports written ---"
ls -lh rbm20_*.report.txt rbm20_*.events.tsv 2>/dev/null
date
