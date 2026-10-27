#!/bin/bash
#SBATCH -p cpu
#SBATCH -c 1
#SBATCH --mem=6G
#SBATCH --qos=ngs16
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#
# A3SS / A5SS / RI -- the same battery already run for SE and MXE.
#
# Job 04 produced <ET>_v2filtered.all_events.tsv for all five event types but
# nothing downstream was run on these three. This fills the gap so every claim
# in the paper can be stated for all five types instead of two.
#
# Per event type:
#   1. HCM-vs-NF phase 2                     -> <ET>_HCM_v2
#   2. composition-adjusted DCM-vs-NF        -> <ET>_v2compadj
#   3. HCM specificity comparison            -> <ET>_specificity
#   4. DIRECT DCM vs HCM, no NF reference    -> <ET>_DCMvsHCM
#   5. composition concordance               -> <ET>_compositionCheck
#   6. why events were lost                  -> <ET>_lostwhy
#
# NOTE ON A3SS vs A5SS: they share identical rMATS column names, so
# --event-type must be passed explicitly. The scripts refuse to guess.
#
# NOTE ON RI: retained introns are the event type most sensitive to library
# prep and RNA degradation -- intronic signal in poly-A data partly reflects
# incomplete splicing and degradation rather than regulated retention. The
# TIN/RIN adjustment matters more here than anywhere else, and RI will NOT be
# comparable to the total-RNA LBNL cohort. Treat RI results with extra caution.
#
# Submit: sbatch 11_other_event_types.sbatch      (~1 hour; these files are
# much smaller than SE)

date
export PATH=/data/apps/anaconda3/envs/mamba/envs/rnaseq_tools/bin:$PATH
ROOT=/data/project/chureeratp/RNA/HEART
cd "$ROOT"

EXPR=magnet_expr.per_sample.tsv
COVH=hcm_nf_covariates.tsv
[ -s "$EXPR" ] || { echo "missing $EXPR -- run 02"; exit 1; }
[ -s "$COVH" ] || { echo "missing $COVH -- run rebuild_hcm_covariates.sh"; exit 1; }

for ET in A3SS A5SS RI; do
  echo
  echo "############################################################"
  echo "#####                    $ET                           #####"
  echo "############################################################"

  BASE=${ET}_v2filtered_rl125.all_events.tsv
  if [ ! -s "$BASE" ]; then
    echo "missing $BASE -- job 04 should have made it. Skipping $ET."
    continue
  fi

  DMATS=rmats_magnet_rl125/${ET}.MATS.JCEC.txt
  HMATS=rmats_HCMvsNF_rl125/${ET}.MATS.JCEC.txt

  # ---- 1. HCM vs NF -------------------------------------------------
  if [ -s "$HMATS" ]; then
    echo; echo "----- $ET : HCM vs NF -----"
    python phase2_diffsplice_v2.py \
      --mats "$HMATS" --event-type $ET \
      --b1 hcm_bams.txt --b2 nf_bams.txt \
      --cov "$COVH" --group1 HCM --group2 NF --min-samples 15 \
      --out ${ET}_HCM_v2_rl125
  else
    echo "missing $HMATS -- skipping HCM arm for $ET"
  fi

  # ---- 2. composition-adjusted discovery ----------------------------
  echo; echo "----- $ET : DCM vs NF + composition -----"
  python phase2_diffsplice_v2.py \
    --mats "$DMATS" --event-type $ET \
    --b1 dcm_bams.txt --b2 nf_bams.txt \
    --cov magnet_covariates.tsv --min-samples 20 \
    --extra-cov-file "$EXPR" --extra-cov CM_score,FB_score \
    --out ${ET}_v2compadj_rl125

  echo; echo "----- $ET : DCM vs NF + composition + technical -----"
  python phase2_diffsplice_v2.py \
    --mats "$DMATS" --event-type $ET \
    --b1 dcm_bams.txt --b2 nf_bams.txt \
    --cov magnet_covariates.tsv --min-samples 20 \
    --extra-cov-file "$EXPR" --extra-cov CM_score,FB_score,log10_libsize,TIN \
    --out ${ET}_v2fulladj_rl125

  # ---- 3. HCM specificity -------------------------------------------
  if [ -s "${ET}_HCM_v2_rl125.all_events.tsv" ]; then
    echo; echo "----- $ET : HCM specificity -----"
    python 05_hcm_specificity.py \
      --dcm "$BASE" --hcm ${ET}_HCM_v2_rl125.all_events.tsv \
      --out ${ET}_specificity_rl125
  fi

  # ---- 4. direct DCM vs HCM -----------------------------------------
  if [ -s "$HMATS" ]; then
    echo; echo "----- $ET : DIRECT DCM vs HCM (no NF reference) -----"
    python 07_dcm_vs_hcm_direct.py \
      --mats-dcm "$DMATS" --mats-hcm "$HMATS" --event-type $ET \
      --b1-dcm dcm_bams.txt --b1-hcm hcm_bams.txt --b2-nf nf_bams.txt \
      --cov-dcm magnet_covariates.tsv --cov-hcm "$COVH" \
      --disc "$BASE" --out ${ET}_DCMvsHCM_rl125
  fi

  # ---- 5/6. composition concordance and loss diagnosis --------------
  echo; echo "----- $ET : composition concordance -----"
  python 06_sensitivity_concordance.py \
    --primary "$BASE" \
    --compare composition=${ET}_v2compadj_rl125.all_events.tsv \
              full=${ET}_v2fulladj_rl125.all_events.tsv \
    --out ${ET}_compositionCheck_rl125

  echo; echo "----- $ET : why events were lost -----"
  python 10_why_events_were_lost.py \
    --base "$BASE" --adj ${ET}_v2compadj_rl125.all_events.tsv \
    --label "$ET composition" --out ${ET}_lostwhy_rl125
done

echo
echo "############### SUMMARY ACROSS ALL FIVE EVENT TYPES ###############"
printf "%-6s %10s %10s %12s %14s\n" TYPE base_sig comp_sig retained DCMvsHCM_sig
for ET in SE MXE A3SS A5SS RI; do
  b=$(awk -F'\t' 'NR>1 && $(NF)<0.05' ${ET}_v2filtered_rl125.all_events.tsv 2>/dev/null | wc -l)
  r=$(grep -m1 "primary-significant events retained" ${ET}_compositionCheck_rl125.report.txt 2>/dev/null \
        | sed 's/.*retained: //; s/ by full key.*//')
  d=$(grep -m1 "significant, FDR<0.05 and |dPSI|>=0.10" ${ET}_DCMvsHCM_rl125.report.txt 2>/dev/null \
        | sed 's/.*: //')
  printf "%-6s %10s %10s %12s %14s\n" "$ET" "${b:-NA}" "-" "${r:-NA}" "${d:-NA}"
done
echo
echo "The DCMvsHCM column is the etiology-independence result. For SE it was 0"
echo "and for MXE 0. If A3SS/A5SS/RI also come back at or near 0, the claim"
echo "holds across all five event types, which is a much stronger statement."
date
