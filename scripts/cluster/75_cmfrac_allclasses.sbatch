#!/bin/bash
#SBATCH --job-name=cmfrac_all
#SBATCH --array=0-3%2
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --time=12:00:00
#SBATCH --qos=ngs16
#SBATCH --output=logs/cmfrac_all_%A_%a.out
#SBATCH --error=logs/cmfrac_all_%A_%a.err
#
# The fraction-only cardiomyocyte score, for the FOUR event classes 72 did not run.
#
# 72_run_cmfrac.sbatch covered SE only (job 406360): 1,054 significant, 716 of 926
# retained = 77.3%, against the nine-marker model's 899 / 590 / 63.7%.
# → RESULT_composition_score_is_a_mediator.md
#
# S1(a) is a FIVE-class table. Without these four it would carry five classes under one
# composition model and one class under the other, which is not a table a reviewer can
# read. This completes it.
#
# RESOURCES: phase2_diffsplice_v2.py is a single serial `for _, row in mats.iterrows()`
#   loop -- no multiprocessing anywhere -- so 2 CPUs is 2 more than it uses. `%2` caps
#   concurrency at two tasks, i.e. 4 CPUs / 40 GB at once, well inside the 16 CPU / 96 GB
#   allowance and leaving room for interactive work.
#
# ⚠ INPUT IS JCEC. 04_phase2_v2_checks.sbatch built every *_v2filtered table from
#   rmats_magnet/${ET}.MATS.JCEC.txt. A JC-derived table compared against a JCEC-derived
#   base is not a comparison, and the 926-style gate would not catch it because the gate
#   reads the base from the existing file. → RESULT_jcec_not_jc.md
#
# ⚠ --event-type IS REQUIRED. A3SS and A5SS share rMATS column names, so detection is
#   ambiguous and phase 2 exits rather than guessing. Passed explicitly for all four.
#
# PREREQUISITE: magnet_expr.cmfrac.tsv (from `python 71_cmfrac_score.py build`).
#   It already exists -- 72 used it -- and its gate passed at r = 1.000000 against the
#   stored nine-marker CM_score.

set -euo pipefail

cd /data/project/chureeratp/RNA/HEART
mkdir -p logs

ETS=(MXE A3SS A5SS RI)
ET=${ETS[$SLURM_ARRAY_TASK_ID]}

SCORE=magnet_expr.cmfrac.tsv
MATS=rmats_magnet_rl125/${ET}.MATS.JCEC.txt
OUT=${ET}_v2cmfrac_rl125

if [[ ! -s "$SCORE" ]]; then
    echo "[fatal] $SCORE missing. Run: python 71_cmfrac_score.py build" >&2
    exit 1
fi
if ! head -1 "$SCORE" | grep -qw CM_score_frac; then
    echo "[fatal] $SCORE has no CM_score_frac column" >&2
    exit 1
fi
if [[ ! -s "$MATS" ]]; then
    echo "[fatal] $MATS missing" >&2
    exit 1
fi

echo "=== $(date) :: ${ET}, DCM vs NF, fraction-only cardiomyocyte score ==="
echo "input : $MATS"
echo "model : PSI ~ disease + sex + age + RIN + CM_score_frac + FB_score"
echo

python phase2_diffsplice_v2.py \
    --mats "$MATS" \
    --b1 dcm_bams.txt \
    --b2 nf_bams.txt \
    --cov magnet_covariates.tsv \
    --event-type "$ET" \
    --min-samples 20 \
    --extra-cov-file "$SCORE" \
    --extra-cov CM_score_frac,FB_score \
    --out "$OUT"

echo
echo "=== $(date) :: ${ET} done ==="
ls -la "${OUT}.all_events.tsv"
echo
echo "When all four array tasks are COMPLETED, rebuild S1 in one pass:"
echo "    python 76_build_s1_v2.py --list"
echo "    python 76_build_s1_v2.py"
echo
echo "Exon-level base counts the rebuild will gate on:"
echo "    SE 926 | MXE 231 | A3SS 90 | A5SS 84 | RI 114"
