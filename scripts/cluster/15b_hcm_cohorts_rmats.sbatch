#!/bin/bash
#SBATCH -p cpu
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH --qos=ngs16
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#
# Three independent HCM cohorts, each against its OWN controls.
#
#   Brigham   PRJNA675090   18 HCM vs 5 ctrl    101bp paired   v49
#   Fuwai     PRJNA522931   28 HCM vs 9 donor   150bp paired   v46
#   Stanford  PRJNA747806   13 HCM vs 7 donor   150bp paired   v49
#
# Run 15a first and confirm COUNTS OK / 0 missing BAMs.
#
# WHY: the MAGNet DCM-vs-HCM comparison shares the same 166 NF hearts between
# both contrasts, so it carries the whole "explanted failing heart vs donor"
# axis. These cohorts have their own controls -- comparing MAGNet DCM-vs-NF
# against, say, Fuwai HCM-vs-donor shares NO samples at all.
#
# WHAT IS *NOT* DONE HERE, DELIBERATELY: the direct DCM-vs-HCM contrast
# (07_dcm_vs_hcm_direct.py). rMATS normalises PSI by IncFormLen, which depends
# on read length. MAGNet is 50bp; these are 101/150bp. Raw PSI is therefore NOT
# comparable across them and a direct per-sample contrast would be invalid.
# The comparison must happen at the dPSI level, where the read-length effect is
# a cohort-level constant that cancels. That is what 05 does.
#
# The GTF passed to rMATS is v49 for every cohort, so EVENT COORDINATES (and
# therefore the `key` column) are identical across runs even for Fuwai, whose
# STAR index was v46. The v46 issue affects alignment sensitivity only.
#
# Submit: sbatch 15b_hcm_cohorts_rmats.sbatch     (allow a few hours)

date
export PATH=/data/apps/anaconda3/envs/mamba/envs/rnaseq_tools/bin:$PATH
ROOT=/data/project/chureeratp/RNA/HEART
GTF=/data/project/chureeratp/hg38/ref/gencode.v49.annotation.gtf
cd "$ROOT"

# cohort  readlength  min-samples (limited by the SMALLER group)
declare -A RL=(   [Brigham]=101 [Fuwai]=150 [Stanford]=150 )
declare -A MINS=( [Brigham]=4   [Fuwai]=7   [Stanford]=5 )

for COH in Brigham Fuwai Stanford; do
  echo
  echo "############################################################"
  echo "#####  $COH : rMATS  (${RL[$COH]}bp paired)"
  echo "############################################################"
  for f in ${COH}_hcm_bams.txt ${COH}_ctrl_bams.txt; do
    [ -s "$f" ] || { echo "missing $f -- run 15a"; exit 1; }
  done
  echo "HCM: $(tr ',' '\n' < ${COH}_hcm_bams.txt | grep -c .)  ctrl: $(tr ',' '\n' < ${COH}_ctrl_bams.txt | grep -c .)"

  # index any BAM missing an index
  cat ${COH}_hcm_bams.txt ${COH}_ctrl_bams.txt | tr ',' '\n' | while read b; do
    [ -z "$b" ] && continue
    [ -f "$b.bai" ] || echo "$b"
  done | xargs -r -P 16 -I{} samtools index {}

  rm -rf ${ROOT}/rmats_${COH} ${ROOT}/rmats_${COH}_tmp
  rmats.py \
    --b1 ${COH}_hcm_bams.txt --b2 ${COH}_ctrl_bams.txt \
    --gtf "$GTF" \
    -t paired --readLength ${RL[$COH]} --variable-read-length --nthread 16 \
    --od ${ROOT}/rmats_${COH} \
    --tmp ${ROOT}/rmats_${COH}_tmp

  # covariate file: group labels only (no age/sex/RIN available)
  printf '%s\tHCM\n'  $(cat ${COH}_hcm.ids)  >  ${COH}.cov
  printf '%s\tCTRL\n' $(cat ${COH}_ctrl.ids) >> ${COH}.cov
done

echo
echo "############################################################"
echo "#####  phase 2 + comparison against MAGNet DCM"
echo "############################################################"

for COH in Brigham Fuwai Stanford; do
  for ET in SE MXE A3SS A5SS RI; do
    M=rmats_${COH}/${ET}.MATS.JCEC.txt
    [ -s "$M" ] || { echo "missing $M -- skipping"; continue; }
    echo
    echo "----- $COH $ET : HCM vs control -----"
    python phase2_diffsplice_v2.py \
      --mats "$M" --event-type $ET \
      --b1 ${COH}_hcm_bams.txt --b2 ${COH}_ctrl_bams.txt \
      --cov ${COH}.cov --group1 HCM --group2 CTRL \
      --min-samples ${MINS[$COH]} --no-covariates \
      --out ${COH}_${ET}_HCM

    DCM=${ET}_v2filtered.all_events.tsv
    [ -s "$DCM" ] || continue
    echo
    echo "----- $COH $ET : MAGNet DCM vs $COH HCM  (INDEPENDENT CONTROLS) -----"
    python 05_hcm_specificity.py \
      --dcm "$DCM" --hcm ${COH}_${ET}_HCM.all_events.tsv \
      --no-shared-controls \
      --out ${COH}_${ET}_vsMAGNetDCM
  done
done

echo
echo "############### SUMMARY ###############"
printf "%-9s %-5s %9s %11s %9s %9s\n" cohort type HCMsig measured Deming sameDir
for COH in Brigham Fuwai Stanford; do
  for ET in SE MXE A3SS A5SS RI; do
    R=${COH}_${ET}_vsMAGNetDCM.report.txt
    [ -s "$R" ] || continue
    s=$(grep -m1 "sig FDR<0.05" ${COH}_${ET}_HCM.all_events.tsv 2>/dev/null)
    n=$(grep -m1 "measured in the HCM run as well" "$R" | sed 's/.*well: //; s/ .*//')
    d=$(grep -m1 "Deming slope" "$R" | sed 's/.*Deming slope = //; s/ .*//')
    sd=$(grep -m1 "same direction:" "$R" | sed 's/.*same direction: //; s/ .*//')
    printf "%-9s %-5s %9s %11s %9s %9s\n" "$COH" "$ET" "-" "${n:-NA}" "${d:-NA}" "${sd:-NA}"
  done
done
echo
echo "READ THIS: with independent controls the shared-control null baseline is 0,"
echo "so the observed correlation is judged against 0 rather than +0.27. A Deming"
echo "slope near 1 across three independent cohorts would mean the failing-"
echo "myocardium signature is genuinely shared with HCM and is NOT an artefact of"
echo "MAGNet's shared donor hearts. A slope near 0 would mean the MAGNet result"
echo "WAS the shared-reference artefact -- which would change the paper."
echo
echo "ALSO: these are septal myectomy samples with preserved EF, not end-stage"
echo "explants. If the slope is high anyway, that argues the signature is not"
echo "merely end-stage remodelling. If it is much lower than MAGNet's ~1.0, that"
echo "difference is the composition/disease-stage effect showing itself."
date
