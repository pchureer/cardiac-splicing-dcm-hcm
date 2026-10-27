#!/bin/bash
#SBATCH -p cpu
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH --qos=ngs16
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#
# RBM20 <-> PSI, WITH technical covariates.
#
# Why: in the first run the matched-background median |t| was 1.0 (SE) and
# 2.0-3.1 (MXE). Under the null it should be ~0.67. That inflation means
# expression and PSI share technical variation -- they come from the same
# library -- so everything correlates with everything a little. It could
# both bury a real RBM20 signal and manufacture the GAPDH one.
#
# This run adds log10_libsize (already in magnet_expr.per_sample.tsv) and
# TIN.median (added below if the pheno file has it) to every model.
#
# READ THE OUTPUT LIKE THIS:
#   - background median |t| should drop toward 0.67. If it does, the
#     inflation was technical and this run is the one to trust.
#   - then compare RBM20 against GAPDH in section 4. If RBM20 still does
#     not beat GAPDH, the mechanistic null is final and can be written up.
#   - section 2 (direction test) cell-adjusted is the mechanism claim.
#
# Submit:  sbatch 03b_rbm20_psi_techadj.sbatch

date
export PATH=/data/apps/anaconda3/envs/mamba/envs/rnaseq_tools/bin:$PATH
ROOT=/data/project/chureeratp/RNA/HEART
PHENO=$ROOT/filereport_read_run_PRJNA595151_withpheno.tsv
cd "$ROOT"

[ -s magnet_expr.per_sample.tsv ] || { echo "run 02 first"; exit 1; }

# --- add TIN.median to the per-sample table if it is not already there ---
# TIN = transcript integrity number, computed from the alignments. A better
# degradation measure than RIN for splicing, and complete for all samples.
if ! head -1 magnet_expr.per_sample.tsv | grep -q "TIN"; then
  if [ -s "$PHENO" ]; then
    TINCOL=$(head -1 "$PHENO" | tr '\t' '\n' | grep -n "TIN" | head -1 | cut -d: -f1)
    if [ -n "$TINCOL" ]; then
      echo "adding TIN from column $TINCOL of the pheno file"
      awk -F'\t' -v c="$TINCOL" 'NR==FNR{t[$1]=$c; next}
           FNR==1{print $0"\tTIN"; next}
           {print $0"\t"(($1 in t && t[$1]!="")?t[$1]:"NA")}' OFS='\t' \
           "$PHENO" magnet_expr.per_sample.tsv > .tmp_expr && \
      mv .tmp_expr magnet_expr.per_sample.tsv
      echo "TIN missing for $(awk -F'\t' 'NR>1 && $NF=="NA"' magnet_expr.per_sample.tsv | wc -l) samples"
    else
      echo "no TIN column found in pheno -- proceeding with libsize only"
    fi
  else
    echo "pheno file not found -- proceeding with libsize only"
  fi
fi
echo "per-sample columns: $(head -1 magnet_expr.per_sample.tsv)"

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

# DON'T-CLOBBER GUARD (03b, techadj).
# The 7 Sep tables are the 50 bp comparators, and 56's *_g45 / *_repro9 reports
# were built from them. Overwriting them orphans those too, and nothing errors.
for ET in SE MXE; do
  for EXT in events.tsv report.txt; do
    F="rbm20_${ET}_techadj${RLSUF}.${EXT}"
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
echo "[guard] 03b, techadj: no existing output for TBL_SUFFIX='${RLSUF}' — safe to write"
# --------------------------------------------------------------------------
# TIN is silently ignored by the script if the column is absent.
for ET in SE MXE; do
  echo
  echo "=================== $ET (technical-covariate adjusted) ==================="
  MATS=${RMATS_DIR}/${ET}.MATS.JC.txt
  [ -s "$MATS" ] || { echo "missing $MATS -- skipping"; continue; }
  python 03_rbm20_psi_correlation.py \
    --rmats "$MATS" --event-type $ET \
    --b1 dcm_bams.txt --b2 nf_bams.txt \
    --expr magnet_expr.per_sample.tsv \
    --covariates magnet_covariates.tsv \
    --numeric-cov age,RIN,log10_libsize,TIN \
    --out-prefix rbm20_${ET}_techadj${RLSUF}
done

echo
echo "--- compare against the unadjusted run ---"
for ET in SE MXE; do
  echo "### $ET background |t| before vs after ###"
  grep -h "matched-bg median" rbm20_${ET}${RLSUF}.report.txt 2>/dev/null | head -2
  echo "  ---"
  grep -h "matched-bg median" rbm20_${ET}_techadj${RLSUF}.report.txt 2>/dev/null | head -2
done
date
