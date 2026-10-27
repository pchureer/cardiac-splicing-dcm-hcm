#!/bin/bash
#SBATCH --job-name=cmfrac
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --time=12:00:00
#SBATCH --qos=ngs16
#SBATCH --output=logs/cmfrac_%j.out
#SBATCH --error=logs/cmfrac_%j.err
#
# Re-run the composition-adjusted discovery contrast with a cardiomyocyte score
# restricted to markers that are NOT upregulated in disease.
#
# WHY: 70_marker_panel_audit.py showed the nine-marker CM score is part myocyte
# fraction and part per-cell sarcomeric upregulation -- four of nine markers are
# strongly up in DCM, and corr(full, fraction-only) = 0.699, so the two are not
# interchangeable as covariates. This produces the number that answers the
# circularity objection instead of arguing it.
#
# NOTHING IN phase2_diffsplice_v2.py CHANGES. It already takes --extra-cov-file
# and --extra-cov, which is how SE_v2compadj was built in the first place.
#
# PREREQUISITE:  python 71_cmfrac_score.py build
#   That writes magnet_expr.cmfrac.tsv and GATES on reconstructing the stored
#   nine-marker CM_score (r > 0.999). If the gate fails, do not run this.
#
# AFTER:         python 71_cmfrac_score.py compare --adj SE_v2cmfrac.all_events.tsv

# RESOURCES: phase2_diffsplice_v2.py is a SINGLE SERIAL LOOP --
#   `for _, row in mats.iterrows()` with no multiprocessing, joblib or thread pool
#   anywhere in the file. Extra CPUs do nothing except make the job queue.
#   2 CPUs / 24 GB leaves room for an 8-CPU job (e.g. the permutation re-run) to hold
#   the rest of a 16-CPU / 96-GB allowance, so this starts immediately instead of
#   pending behind it.

# ⚠ INPUT MUST BE JCEC, NOT JC.
#   04_phase2_v2_checks.sbatch:21 built SE_v2filtered from rmats_magnet/SE.MATS.JCEC.txt
#   ("same file your phase-2 used"), and scripts 05/08/11/12/15b/16/19 all use JCEC too.
#   The first submission of this script (job 406355) passed JC and was cancelled: it
#   would have produced a JC-derived table compared against a JCEC-derived base, and the
#   926 gate would still have passed because the base is read from the existing file.
#   A comparison across two quantifications is not a comparison.
#   NOTE: the Methods currently states JC was used in preference to JCEC. That is wrong
#   and is being corrected -- see MANUSCRIPT_TODO.

set -euo pipefail

cd /data/project/chureeratp/RNA/HEART
mkdir -p logs

SCORE=magnet_expr.cmfrac.tsv
OUT=SE_v2cmfrac_rl125

if [[ ! -s "$SCORE" ]]; then
    echo "[fatal] $SCORE missing. Run: python 71_cmfrac_score.py build" >&2
    exit 1
fi

# Confirm the column exists before burning the allocation on a run that will
# exit at the argument check.
if ! head -1 "$SCORE" | grep -qw CM_score_frac; then
    echo "[fatal] $SCORE has no CM_score_frac column" >&2
    head -1 "$SCORE" >&2
    exit 1
fi

echo "=== $(date) :: SE, DCM vs NF, fraction-only cardiomyocyte score ==="
echo "model: PSI ~ disease + sex + age + RIN + CM_score_frac + FB_score"
echo

python phase2_diffsplice_v2.py \
    --mats rmats_magnet_rl125/SE.MATS.JCEC.txt \
    --b1 dcm_bams.txt \
    --b2 nf_bams.txt \
    --cov magnet_covariates.tsv \
    --event-type SE \
    --min-samples 20 \
    --extra-cov-file "$SCORE" \
    --extra-cov CM_score_frac,FB_score \
    --out "$OUT"

echo
echo "=== $(date) :: done ==="
ls -la "${OUT}.all_events.tsv"
echo
echo "Now run:  python 71_cmfrac_score.py compare --adj ${OUT}.all_events.tsv"
echo
echo "The comparison GATES on the base table still returning 926. The number to"
echo "watch is retention against the recorded 590/926 = 63.7%:"
echo "  * unchanged (within ~3 points) -> the objection is answered, quote both"
echo "  * HIGHER                       -> the nine-marker score was over-adjusting;"
echo "                                    that difference IS the over-adjustment"
echo "  * LOWER                        -> not the predicted direction. Check whether"
echo "                                    five markers simply made the covariate"
echo "                                    noisier before concluding anything."
