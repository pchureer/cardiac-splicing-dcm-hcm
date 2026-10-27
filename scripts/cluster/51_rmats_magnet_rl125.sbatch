#!/bin/bash
#SBATCH -J rmats_rl125
#SBATCH -p cpu
#SBATCH -c 8
#SBATCH --mem=96G
#SBATCH --qos=ngs16
#SBATCH -t 8-00:00:00
#SBATCH --array=1-2%1
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#
# 51_rmats_magnet_rl125.sbatch -- re-run the two MAGNet contrasts with the read length
# the data actually has.
#
# WHY
# ---
# Both MAGNet contrasts were quantified with `--readLength 50`; the reads are 125 bp
# (48_provenance_versions.sh section 5, 49_readlength_audit.sh across samples). rMATS
# normalises PSI by effective form lengths computed from --readLength, so the declared
# value enters PSI itself. 50_readlength_impact.py measures the consequence once these
# finish.
#
# THE ONE RULE THIS SCRIPT EXISTS TO ENFORCE
# ------------------------------------------
# The new run must differ from the original in --readLength AND NOTHING ELSE. Same bam
# lists, same GTF, same -t, same --variable-read-length, same --nthread, same rMATS
# build, single-step (no --task prep/post). If anything else changes, the comparison
# stops being a test of the read length and becomes a test of two pipelines.
#
# That is also why this does NOT use prep/post to share the 166 NF samples between the
# two contrasts, even though it would save real time: the originals were produced
# single-step, and a second difference is not worth the hours.
#
# IT WILL NOT OVERWRITE ANYTHING
# ------------------------------
# This project has already lost an rMATS output to an `rm -rf` in a job script. There is
# no rm here. If an output directory exists the task refuses to start, and the original
# rmats_magnet / rmats_HCMvsNF directories are never touched or named.
#
# CORE BUDGET
# -----------
# Sized for a 16-core / 96 GB account limit. The first attempt at -c 12 --mem 64G was
# SIGKILLed twice. rMATS prep holds reads PER THREAD, so the thread count is the main
# driver of peak memory: -c 4 cuts it to roughly a third while --mem 90G takes almost
# all the headroom the account allows. It will be slow. Run 52_rl_probe.sbatch FIRST --
# if the probe says the read-length constant does not change the result, this full
# rerun is not needed at all.
#
# Running them sequentially at 12 is also FASTER here than running both at once on 6
# each, because the two tasks are unequal: task 1 has 332 samples and task 2 has 194,
# so a concurrent pair is bottlenecked by task 1 at half the threads while task 2 sits
# finished. Sequential at 12 beats parallel at 6 for any split this lopsided.
#
# If your limit changes, edit `#SBATCH -c` only -- --nthread reads it.
#
# Usage:
#   cd /data/project/chureeratp/RNA/HEART
#   sbatch 51_rmats_magnet_rl125.sbatch
#   sbatch --array=1 51_rmats_magnet_rl125.sbatch     # just DCM vs NF

set -uo pipefail
date
echo "host: $(hostname)   job: ${SLURM_JOB_ID:-NA} task: ${SLURM_ARRAY_TASK_ID:-NA}"
echo "cores: ${SLURM_CPUS_PER_TASK:-unset}   mem: ${SLURM_MEM_PER_NODE:-unset} MB"

export PATH=/data/apps/anaconda3/envs/mamba/envs/rnaseq_tools/bin:$PATH
# rmats.py is Python and buffers stdout. The first attempt was SIGKILLed and took its
# buffered output to the grave, so the .out file ended at the command echo and said
# nothing about WHERE it died. Unbuffered, the log survives the kill.
export PYTHONUNBUFFERED=1
ROOT=/data/project/chureeratp/RNA/HEART
GTF=/data/project/chureeratp/hg38/ref/gencode.v49.annotation.gtf
RL=125                 # the whole point. Everything else mirrors the original run.
# One source of truth for the core count: change `#SBATCH -c` and --nthread follows.
# Thread count is parallelism only -- it does not change what rMATS computes -- so
# running at 12 where the original ran at 16 is NOT a second analytical difference.
NTHREAD=${SLURM_CPUS_PER_TASK:-12}

cd "$ROOT" || { echo "cannot cd to $ROOT"; exit 1; }

case "${SLURM_ARRAY_TASK_ID:-1}" in
  1) B1=dcm_bams.txt; B2=nf_bams.txt; OD=rmats_magnet_rl125;  N1=166; N2=166
     LABEL="MAGNet DCM vs NF" ;;
  2) B1=hcm_bams.txt; B2=nf_bams.txt; OD=rmats_HCMvsNF_rl125; N1=28;  N2=166
     LABEL="MAGNet HCM vs NF" ;;
  3) B1=dcm_bams.txt; B2=hcm_bams.txt; OD=rmats_DCMvsHCM_rl125; N1=166; N2=28
     LABEL="MAGNet DCM vs HCM" ;;
  *) echo "unknown array index"; exit 2 ;;
esac
TMP="${OD}_tmp"

echo "=============================================================="
echo "$LABEL   --readLength $RL"
echo "=============================================================="

# --------------------------------------------------------------------------
# preflight. Every check here corresponds to a failure this project has already
# had, or to one rMATS performs silently.
# --------------------------------------------------------------------------
fail=0

nbam () { tr ',' '\n' < "$1" | sed 's/[[:space:]]//g' | grep -c '\.bam$'; }

for f in "$B1" "$B2" "$GTF"; do
  if [[ ! -s "$f" ]]; then echo "  *** missing or empty: $f"; fail=1; fi
done
[[ $fail -eq 1 ]] && { echo "preflight failed"; exit 1; }

n1=$(nbam "$B1"); n2=$(nbam "$B2")
echo "  $B1: $n1 bams (expected $N1)"
echo "  $B2: $n2 bams (expected $N2)"
# A short group is this project's documented silent failure: the list builder used
# `ls | head -1` per id and did not fail on a missing bam, so the group just shrank.
[[ "$n1" -ne "$N1" ]] && { echo "  *** group 1 is not $N1 samples"; fail=1; }
[[ "$n2" -ne "$N2" ]] && { echo "  *** group 2 is not $N2 samples"; fail=1; }

missing=0
while IFS= read -r b; do
  [[ -z "$b" ]] && continue
  [[ -f "$b" ]] || { echo "  *** bam not on disk: $b"; missing=$((missing+1)); }
done < <(cat "$B1" "$B2" | tr ',' '\n' | sed 's/[[:space:]]//g' | grep '\.bam$')
[[ $missing -gt 0 ]] && { echo "  *** $missing bam path(s) do not exist"; fail=1; }

# rMATS must read the RAW STAR alignments. The bqsr bams have been through
# SplitNCigarReads, which strips the N-CIGAR junction reads splicing depends on.
if grep -qiE 'bqsr|SplitNC|dedup\.bam' <(cat "$B1" "$B2"); then
  echo "  *** a list points at bqsr/SplitNCigar bams — rMATS needs the raw STAR bams"
  fail=1
fi

for d in "$OD" "$TMP"; do
  if [[ -e "$d" ]]; then
    echo "  *** $d already exists. REFUSING to start — it would be overwritten."
    echo "      Move it aside yourself if you meant to redo this. Nothing is deleted"
    echo "      by this script, deliberately."
    fail=1
  fi
done

[[ $fail -eq 1 ]] && { echo; echo "PREFLIGHT FAILED — nothing was run."; exit 1; }
echo "  preflight ok"

# --------------------------------------------------------------------------
echo
echo "command:"
echo "  rmats.py --b1 $B1 --b2 $B2 --gtf $GTF -t paired \\"
echo "    --readLength $RL --variable-read-length --nthread $NTHREAD \\"
echo "    --od $OD --tmp $TMP"
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
  echo "*** rMATS failed — do not use this output"
  echo
  echo "    Exit 137, or sacct showing ExitCode 0:9, means SIGKILL — on SLURM that is"
  echo "    almost always the memory limit rather than a bug. Check what it used:"
  echo "      sacct -j \$SLURM_JOB_ID -o JobID%14,State%12,ExitCode,MaxRSS,ReqMem"
  echo "    then raise --mem above that MaxRSS and resubmit. Lowering -c also lowers"
  echo "    peak memory, because rMATS prep holds reads per thread."
  echo
  echo "    The partial output dirs are NOT deleted. Move or remove them yourself"
  echo "    before resubmitting — the preflight guard will otherwise refuse to start:"
  echo "      rm -rf $OD $TMP"
  date
  exit $rc
fi

# --------------------------------------------------------------------------
echo "=============================================================="
echo "POST-RUN CHECKS"
echo "=============================================================="
for et in SE MXE A3SS A5SS RI; do
  f="$OD/${et}.MATS.JC.txt"
  if [[ -s "$f" ]]; then
    printf '  %-5s %8s rows\n' "$et" "$(( $(wc -l < "$f") - 1 ))"
  else
    echo "  *** $et: $f missing or empty"
  fi
done

# "USED: 0 / NOT_PAIRED: huge" is what a wrong -t looks like, and it still writes
# output files that look normal.
echo
echo "  read-pairing summary from the prep logs:"
grep -rhoE '(USED|NOT_PAIRED|NOT_NH_1|NOT_EXPECTED_CIGAR)[: ]+[0-9]+' "$TMP" 2>/dev/null \
  | awk -F'[: ]+' '{s[$1]+=$2} END {for (k in s) printf "    %-20s %s\n", k, s[k]}'
echo "  USED near zero with NOT_PAIRED large means -t is wrong. It is NOT wrong here"
echo "  unless the original run was also wrong, since -t paired mirrors it."

# group widths: the other silent failure mode
hdr=$(head -1 "$OD/SE.MATS.JC.txt" 2>/dev/null)
if [[ -n "$hdr" ]]; then
  i1=$(awk -F'\t' -v h="$hdr" 'BEGIN{n=split(h,a,"\t"); for(i=1;i<=n;i++) if(a[i]=="IncLevel1") print i}')
  i2=$(awk -F'\t' -v h="$hdr" 'BEGIN{n=split(h,a,"\t"); for(i=1;i<=n;i++) if(a[i]=="IncLevel2") print i}')
  if [[ -n "$i1" && -n "$i2" ]]; then
    w1=$(sed -n '2p' "$OD/SE.MATS.JC.txt" | cut -f"$i1" | awk -F',' '{print NF}')
    w2=$(sed -n '2p' "$OD/SE.MATS.JC.txt" | cut -f"$i2" | awk -F',' '{print NF}')
    echo
    echo "  IncLevel1 width: $w1 (expected $N1)"
    echo "  IncLevel2 width: $w2 (expected $N2)"
    [[ "$w1" != "$N1" || "$w2" != "$N2" ]] && \
      echo "  *** PSI column widths do not match the group sizes — do not use this run"
  fi
fi

echo
echo "=============================================================="
echo "NEXT"
echo "=============================================================="
echo "  When BOTH array tasks have finished:"
echo
echo "    python 50_readlength_impact.py --old <original DCMvsNF dir> \\"
echo "        --new rmats_magnet_rl125"
echo
echo "  Stage A of that script answers whether the read length mattered at all."
echo "  Stage B needs phase2_diffsplice_v2.py rerun on this output first; until then"
echo "  the manuscript's 926 / 1,445 are untested and must not be changed."
echo "  Change nothing in Table 1, Figure 1 or the text before Stage A is read."
date
