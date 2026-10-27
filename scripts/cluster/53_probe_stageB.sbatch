#!/bin/bash
#SBATCH -J probeB
#SBATCH -p cpu
#SBATCH -c 2
#SBATCH --mem=24G
#SBATCH --qos=ngs16
#SBATCH -t 08:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#
# 53_probe_stageB.sbatch -- Stage B of the read-length probe.
#
# Runs phase2_diffsplice_v2.py on BOTH probe rMATS outputs, then compares them with
# 50_readlength_impact.py --probe. This is the covariate-adjusted half of the test:
# Stage A showed the read-length constant does not move raw rMATS dPSI; this asks
# whether the model that actually produces the manuscript's 926 reacts differently.
#
# WHY IT IS A JOB AND NOT A ONE-LINER
# -----------------------------------
# phase2 fits a linear model per event. The probe has ~197,000 SE events, so each run
# is ~197k model fits on 60 samples. That is tens of minutes of single-threaded Python
# per read length, which is too long for a login node and looks like a hang because
# Python buffers its output. PYTHONUNBUFFERED fixes the second problem.
#
# IDEMPOTENT: if an output table already exists it is reused rather than recomputed,
# so a resubmission after a timeout resumes instead of starting over. Nothing is
# deleted. Remove a table yourself if you want it rebuilt.
#
# Usage:
#   cd /data/project/chureeratp/RNA/HEART
#   sbatch 53_probe_stageB.sbatch
#   sbatch 53_probe_stageB.sbatch --event-type RI    # a second event class
#
set -uo pipefail
date
echo "host: $(hostname)   job: ${SLURM_JOB_ID:-NA}"
echo "cores: ${SLURM_CPUS_PER_TASK:-unset}   mem: ${SLURM_MEM_PER_NODE:-unset} MB"

export PATH=/data/apps/anaconda3/envs/mamba/envs/rnaseq_tools/bin:$PATH
export PYTHONUNBUFFERED=1          # phase2 prints progress; buffered, it looks hung

ROOT=/data/project/chureeratp/RNA/HEART
COV=magnet_covariates.tsv
B1=probe_dcm.txt
B2=probe_nf.txt
ET=SE

while [[ $# -gt 0 ]]; do
  case "$1" in
    --event-type) ET="$2"; shift 2 ;;
    --cov)        COV="$2"; shift 2 ;;
    *) echo "unknown option: $1"; exit 2 ;;
  esac
done

cd "$ROOT" || { echo "cannot cd to $ROOT"; exit 1; }

echo "=============================================================="
echo "STAGE B — probe, event type $ET"
echo "=============================================================="

fail=0
for f in "$B1" "$B2" "$COV" phase2_diffsplice_v2.py 50_readlength_impact.py; do
  [[ -s "$f" ]] || { echo "  *** missing: $f"; fail=1; }
done
for d in probe_rl50 probe_rl125; do
  [[ -s "$d/${ET}.MATS.JC.txt" ]] || { echo "  *** missing: $d/${ET}.MATS.JC.txt"; fail=1; }
done
# Both runs MUST have used the same samples, or the comparison is meaningless. The
# probe printed the list hashes; this re-derives them so the job can refuse on its own.
echo "  sample lists: $(md5sum "$B1" "$B2" | awk '{printf "%s ", $1}')"
n1=$(tr ',' '\n' < "$B1" | grep -c '\.bam$'); n2=$(tr ',' '\n' < "$B2" | grep -c '\.bam$')
echo "  group sizes : $n1 vs $n2"
[[ "$n1" -lt 20 || "$n2" -lt 20 ]] && {
  echo "  *** phase2 --min-samples defaults to 20; a group below that returns 0 events"
  fail=1; }
[[ $fail -eq 1 ]] && { echo; echo "PREFLIGHT FAILED — nothing was run."; exit 1; }
echo "  preflight ok"

for RL in 50 125; do
  OUT="probe_rl${RL}_${ET}"
  if [[ -s "${OUT}.all_events.tsv" ]]; then
    echo
    echo "  [skip] ${OUT}.all_events.tsv already exists — reusing it."
    echo "         Delete it yourself if you want it rebuilt."
    continue
  fi
  echo
  echo "--------------------------------------------------------------"
  echo "phase2 on probe_rl${RL}  ($(date +%H:%M:%S))"
  echo "--------------------------------------------------------------"
  python phase2_diffsplice_v2.py \
    --mats "probe_rl${RL}/${ET}.MATS.JC.txt" \
    --b1 "$B1" --b2 "$B2" \
    --cov "$COV" \
    --out "$OUT"
  rc=$?
  echo "phase2 rl${RL} exit status: $rc"
  [[ $rc -ne 0 ]] && { echo "*** phase2 failed — stopping"; date; exit $rc; }
done

echo
echo "=============================================================="
echo "THE COMPARISON"
echo "=============================================================="
python 50_readlength_impact.py \
  --old probe_rl50 --new probe_rl125 \
  --types "$ET" \
  --probe \
  --old-adj "probe_rl50_${ET}.all_events.tsv" \
  --new-adj "probe_rl125_${ET}.all_events.tsv"

echo
echo "=============================================================="
echo "HOW TO READ IT"
echo "=============================================================="
echo "  Stage B's two counts are the covariate-adjusted ones, on 30 vs 30. They are NOT"
echo "  the manuscript's 926 and must never be quoted as such."
echo
echo "  The only question is whether the ADJUSTED counts move more than the raw rMATS"
echo "  counts did (SE was -2.9% in Stage A). If they move comparably, the covariate"
echo "  model inherits the read-length constant exactly as expected and verdict [A]"
echo "  stands for the pipeline as a whole, not just for raw rMATS."
echo
echo "  If the adjusted counts move MUCH more than the raw ones, something in the model"
echo "  is amplifying the constant and the full cohort rerun becomes necessary after"
echo "  all. That would be a surprise — say so rather than explaining it away."
date
