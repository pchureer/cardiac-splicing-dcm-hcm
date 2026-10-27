#!/bin/bash
#SBATCH -J rlprobe
#SBATCH -p cpu
#SBATCH -c 6
#SBATCH --mem=80G
#SBATCH --qos=ngs16
#SBATCH -t 12:00:00
#SBATCH --array=1-2%1
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#
# 52_rl_probe.sbatch -- does --readLength 50 vs 125 change anything? Answer it on a
# subset, cheaply, before spending the cluster on the full cohort again.
#
# WHY A SUBSET IS A VALID TEST HERE
# ---------------------------------
# The question is not "what are MAGNet's true counts". It is "does the declared read
# length change PSI". That is a property of the normalisation constant, and it is
# isolated perfectly by running THE SAME SAMPLES at both read lengths. The subset is
# irrelevant to the comparison: every event is measured twice on identical reads, and
# the only thing that differs is the number rMATS was told.
#
# A subset is NOT valid for anything else. The counts it produces are not MAGNet's
# counts and must never be quoted. They exist to answer one yes/no question.
#
# WHY THIS EXISTS AT ALL
# ----------------------
# The full run (51) was SIGKILLed twice at --mem=64G with -c 12. The account caps at
# 96 GB and 16 cores, so there is little headroom to brute-force. rMATS prep holds
# reads per thread, so 6 threads on 60 bams fits easily where 12 threads on 332 did
# not. If the probe says the constant does not matter, the full rerun is unnecessary.
# If it says it does, the full rerun is unavoidable and worth queueing properly.
#
# Usage:
#   cd /data/project/chureeratp/RNA/HEART
#   sbatch 52_rl_probe.sbatch
#   sbatch 52_rl_probe.sbatch --n 20        # smaller still
#
# Then:
#   python 50_readlength_impact.py --old probe_rl50 --new probe_rl125

set -uo pipefail
date
echo "host: $(hostname)   job: ${SLURM_JOB_ID:-NA} task: ${SLURM_ARRAY_TASK_ID:-NA}"
echo "cores: ${SLURM_CPUS_PER_TASK:-unset}   mem: ${SLURM_MEM_PER_NODE:-unset} MB"

export PATH=/data/apps/anaconda3/envs/mamba/envs/rnaseq_tools/bin:$PATH
export PYTHONUNBUFFERED=1          # or a SIGKILL takes the progress log with it

ROOT=/data/project/chureeratp/RNA/HEART
GTF=/data/project/chureeratp/hg38/ref/gencode.v49.annotation.gtf
N=30                               # samples per group
NTHREAD=${SLURM_CPUS_PER_TASK:-6}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --n) N="$2"; shift 2 ;;
    *) echo "unknown option: $1"; exit 2 ;;
  esac
done

cd "$ROOT" || { echo "cannot cd to $ROOT"; exit 1; }

case "${SLURM_ARRAY_TASK_ID:-1}" in
  1) RL=50;  OD=probe_rl50 ;;
  2) RL=125; OD=probe_rl125 ;;
  *) echo "unknown array index"; exit 2 ;;
esac
TMP="${OD}_tmp"
B1=probe_dcm.txt
B2=probe_nf.txt

echo "=============================================================="
echo "PROBE  --readLength $RL   ($N vs $N samples)"
echo "=============================================================="

# --------------------------------------------------------------------------
# Build the subset ONCE, deterministically, and reuse it for both tasks. Both read
# lengths MUST see byte-identical lists or the comparison is not a comparison.
# --------------------------------------------------------------------------
if [[ ! -s "$B1" || ! -s "$B2" ]]; then
  echo "  building subset lists (first $N of each group, file order)"
  tr ',' '\n' < dcm_bams.txt | sed 's/[[:space:]]//g' | grep '\.bam$' | head -n "$N" \
    | paste -sd, - > "$B1"
  tr ',' '\n' < nf_bams.txt  | sed 's/[[:space:]]//g' | grep '\.bam$' | head -n "$N" \
    | paste -sd, - > "$B2"
else
  echo "  reusing existing $B1 / $B2 — both tasks must use the same samples"
fi

nbam () { tr ',' '\n' < "$1" | sed 's/[[:space:]]//g' | grep -c '\.bam$'; }
n1=$(nbam "$B1"); n2=$(nbam "$B2")
echo "  $B1: $n1 bams"
echo "  $B2: $n2 bams"
echo "  md5: $(md5sum "$B1" "$B2" | awk '{printf "%s ", $1}')"

fail=0
[[ "$n1" -ne "$N" || "$n2" -ne "$N" ]] && { echo "  *** subset is not $N vs $N"; fail=1; }
[[ -s "$GTF" ]] || { echo "  *** GTF missing: $GTF"; fail=1; }
while IFS= read -r b; do
  [[ -z "$b" ]] && continue
  [[ -f "$b" ]] || { echo "  *** bam not on disk: $b"; fail=1; }
done < <(cat "$B1" "$B2" | tr ',' '\n' | sed 's/[[:space:]]//g' | grep '\.bam$')
for d in "$OD" "$TMP"; do
  [[ -e "$d" ]] && { echo "  *** $d exists — refusing to overwrite. rm -rf it first."; fail=1; }
done
[[ $fail -eq 1 ]] && { echo; echo "PREFLIGHT FAILED — nothing was run."; exit 1; }
echo "  preflight ok"

echo
rmats.py \
  --b1 "$B1" \
  --b2 "$B2" \
  --gtf "$GTF" \
  -t paired \
  --readLength "$RL" \
  --variable-read-length \
  --nthread "$NTHREAD" \
  --od "$OD" \
  --tmp "$TMP"
rc=$?
echo
echo "rmats exit status: $rc"
if [[ $rc -ne 0 ]]; then
  echo "*** failed. If sacct shows ExitCode 0:9 this was SIGKILL = memory."
  echo "    sacct -j \$SLURM_JOB_ID -o JobID%14,State%12,ExitCode,MaxRSS,ReqMem"
  echo "    Drop -c to 4 and resubmit, or --n 20. Remove $OD $TMP first."
  date; exit $rc
fi

for et in SE MXE A3SS A5SS RI; do
  f="$OD/${et}.MATS.JC.txt"
  [[ -s "$f" ]] && printf '  %-5s %8s rows\n' "$et" "$(( $(wc -l < "$f") - 1 ))" \
                || echo "  *** $et missing"
done

echo
echo "=============================================================="
echo "NEXT — after BOTH tasks finish"
echo "=============================================================="
echo "  python 50_readlength_impact.py --old probe_rl50 --new probe_rl125"
echo
echo "  Read the Inc/Skip median ratio first: that is what the constant changes."
echo "  Then sign concordance, then the count change."
echo
echo "  Verdict [A]  -> the constant does not matter. The frozen 926 / 1,445 / 0 stand"
echo "                  and the full rerun is UNNECESSARY. Fix Table 1, rewrite the"
echo "                  attenuation sentence, disclose it in the Methods, move on."
echo "  Verdict [B]/[C] -> the full MAGNet rerun is unavoidable. Queue it with -c 4"
echo "                  and --mem 90G and expect it to be slow."
echo
echo "  The probe's own counts are NOT MAGNet's counts. Never quote them."
date
