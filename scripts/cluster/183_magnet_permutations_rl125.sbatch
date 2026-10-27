#!/bin/bash
#SBATCH -J magnet_perm_rl125
#SBATCH -p cpu
#SBATCH -c 8
#SBATCH --mem=16G
#SBATCH --qos=ngs16
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#
# 183_magnet_permutations_rl125.sbatch -- the two MAGNet label permutations of
# Supplementary Table 8(a), repeated on the 125 bp quantification.
#
# 25_colorado_inflation.sbatch section D ran them on rmats_magnet (50 bp): seeds 101
# and 102, SE only, phase2_diffsplice_v2.py with the covariate file's disease labels
# shuffled, --min-samples 20. Both returned exactly zero significant events across
# 103,054 tested. Only MAGNet was re-quantified, so only MAGNet's two rows need
# repeating; the twenty Colorado analyses are unchanged.
#
# WHAT THIS JOB DOES, IN ORDER
#   0. refuses unless every input exists and 25's shuffle block can be extracted
#   1. REAL model (unshuffled labels) on the 125 bp table -- the calibration. It must
#      reproduce the discovery count (944 full-coordinate / 813 exon-level) or the job
#      stops before any shuffle runs: a permutation null is only interpretable beside
#      the positive control from the same invocation.
#   2. seeds 101 and 102, shuffled with 25's OWN python block, extracted from 25 at run
#      time so the shuffle is identical by construction rather than by transcription.
#   3. counts significant events under BOTH conventions (full-coordinate as 25 did,
#      exon-level as the manuscript does) and writes a summary in 25's printf format.
#
# Writes only under inflate_rl125/. The 7 Oct S8a_permutations.tsv is untouched; its
# two MAGNet rows are replaced afterwards by 184 (next step), not by this job.
#
#   sbatch 183_magnet_permutations_rl125.sbatch
#   # ~1-2 h: three phase2 runs on ~100k events x 332 samples. 25's whole run was 2h08
#   # for 23 analyses; the 50 bp phase2 needed 1.6 GB, so 16G is generous.

set -euo pipefail
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate rnaseq_tools
cd /data/project/chureeratp/RNA/HEART

MATS=rmats_magnet_rl125/SE.MATS.JCEC.txt
COV=magnet_covariates.tsv
OUTDIR=inflate_rl125
EXPECT_FULL=944      # NUMBERS_125bp_REDERIVATION: SE sig full-coordinate at 125 bp
EXPECT_EXON=813      # the manuscript's discovery count
SEEDS="101 102"

echo "##### 0) inputs"
for f in "$MATS" "$COV" dcm_bams.txt nf_bams.txt phase2_diffsplice_v2.py \
         25_colorado_inflation.sbatch; do
  [ -s "$f" ] || { echo "[REFUSE] missing: $f"; exit 2; }
  echo "  ok $f"
done
mkdir -p "$OUTDIR"

# 25's shuffle block, verbatim. It sits between the line
#     python3 - magnet_covariates.tsv inflate/magnet.P${s}.cov $s <<'PY'
# and the closing PY. Extracting it means the permutation is the same code 25 ran.
awk '/python3 - magnet_covariates.tsv/{grab=1; next} grab && /^PY$/{exit} grab' \
    25_colorado_inflation.sbatch > "$OUTDIR/shuffle_from_25.py"
if [ ! -s "$OUTDIR/shuffle_from_25.py" ] || ! grep -q "shuffle" "$OUTDIR/shuffle_from_25.py"; then
  echo "[REFUSE] could not extract 25's MAGNet shuffle block; check lines ~233-241 of 25"
  exit 2
fi
echo "  shuffle block extracted from 25 ($(wc -l < "$OUTDIR/shuffle_from_25.py") lines):"
sed 's/^/      | /' "$OUTDIR/shuffle_from_25.py"

count_sig () {   # $1 = all_events.tsv ; prints "tested_full sig_full tested_exon sig_exon"
python3 - "$1" <<'PY'
import sys, pandas as pd
t = pd.read_csv(sys.argv[1], sep="\t", low_memory=False)
d = next(c for c in ("dPSI_adj", "dPSI_raw", "dPSI") if c in t.columns)
t["dPSI"] = pd.to_numeric(t[d], errors="coerce"); t["FDR"] = pd.to_numeric(t["FDR"], errors="coerce")
full = t.drop_duplicates("key")                                   # 25's convention
sig_full = int(((full.FDR < 0.05) & (full.dPSI.abs() >= 0.10)).sum())
t["ekey"] = [":".join(str(k).split(":")[:3]) for k in t["key"]]   # exon level, lowest FDR kept
exon = t.sort_values("FDR", kind="mergesort").drop_duplicates("ekey")
sig_exon = int(((exon.FDR < 0.05) & (exon.dPSI.abs() >= 0.10)).sum())
print(len(full), sig_full, len(exon), sig_exon)
PY
}

run_phase2 () {    # $1 = covariate file ; $2 = output prefix
  python phase2_diffsplice_v2.py \
    --mats "$MATS" --event-type SE \
    --b1 dcm_bams.txt --b2 nf_bams.txt --cov "$1" \
    --min-samples 20 --out "$2" > "$2.log" 2>&1
}

echo
echo "##### 1) REAL model, unshuffled labels -- the calibration"
run_phase2 "$COV" "$OUTDIR/magnet_REAL_SE"
read tf sf te se <<< "$(count_sig "$OUTDIR/magnet_REAL_SE.all_events.tsv")"
printf "  MAGNet SE real model          tested=%-8s sig=%s   (exon level: %s / %s)\n" "$tf" "$sf" "$te" "$se"
if [ "$sf" != "$EXPECT_FULL" ] || [ "$se" != "$EXPECT_EXON" ]; then
  echo "[REFUSE] real model gives $sf / $se, expected $EXPECT_FULL / $EXPECT_EXON."
  echo "         The positive control does not reproduce; the shuffles would be uninterpretable."
  echo "         Check --cov and --min-samples against 08_composition_adjusted_discovery / 04."
  exit 3
fi
echo "  [ok] real model reproduces the discovery count at 125 bp"

echo
echo "##### 2) label permutations, seeds $SEEDS"
SUMMARY="$OUTDIR/magnet_permutations_rl125.summary.txt"
: > "$SUMMARY"
printf "  MAGNet SE real model          tested=%-8s sig=%s\n" "$tf" "$sf" >> "$SUMMARY"
for s in $SEEDS; do
  python3 "$OUTDIR/shuffle_from_25.py" "$COV" "$OUTDIR/magnet.P${s}.cov" "$s"
  run_phase2 "$OUTDIR/magnet.P${s}.cov" "$OUTDIR/magnet_P${s}_SE"
  n=$(grep -m1 "events tested:" "$OUTDIR/magnet_P${s}_SE.log" | grep -o '[0-9]*' | head -1 || true)
  read tf sf te se <<< "$(count_sig "$OUTDIR/magnet_P${s}_SE.all_events.tsv")"
  # 25's own line format, so the record reads the same way at both read lengths
  printf "  MAGNet SE permutation seed %-5s tested=%-8s sig=%s\n" "$s" "${n:-$tf}" "$sf" | tee -a "$SUMMARY"
  printf "      (exon level: tested=%s sig=%s)\n" "$te" "$se" | tee -a "$SUMMARY"
done

echo
echo "##### 3) verdict"
cat "$SUMMARY"
if grep "permutation" "$SUMMARY" | grep -qv "sig=0$"; then
  echo "  >>> a permutation returned a NON-ZERO count. Do not update the caption until this is understood."
  exit 4
fi
echo "  both permutations return exactly zero at 125 bp, beside a real model of $EXPECT_FULL / $EXPECT_EXON."
echo "  Next: python 184_update_s8a_magnet_rl125.py   (replaces the two MAGNet rows of S8a)"
