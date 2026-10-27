#!/bin/bash
#SBATCH -p cpu
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH --qos=ngs16
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#
# HCM-vs-NF phase 2 + DCM-specificity comparison.
# Run VERIFY_HCM_RUN.sh first -- if the HCM group is not 28 samples, or rMATS
# reported USED: 0, nothing below is meaningful.
# Submit:  sbatch 05_hcm_specificity.sbatch

date
export PATH=/data/apps/anaconda3/envs/mamba/envs/rnaseq_tools/bin:$PATH
ROOT=/data/project/chureeratp/RNA/HEART
cd "$ROOT"

# The HCM covariate file must contain the 28 HCM samples AND the 166 NF samples,
# with the disease column reading HCM / NF. magnet_covariates.tsv does NOT --
# it is exactly 166 DCM + 166 NF.
COV=hcm_nf_covariates.tsv
if [ ! -s "$COV" ]; then
  echo "MISSING $COV -- see the note in the chat about building it"; exit 1
fi
echo "HCM covariate rows: $(wc -l < $COV)   (expect 28 HCM + 166 NF = 194)"
echo "labels present:"; cut -f2 "$COV" | sort | uniq -c

for ET in SE MXE; do
  echo
  echo "############ $ET : HCM vs NF ############"
  MATS=rmats_HCMvsNF/${ET}.MATS.JCEC.txt
  [ -s "$MATS" ] || { echo "missing $MATS -- skipping"; continue; }
  python phase2_diffsplice_v2.py \
    --mats "$MATS" --b1 hcm_bams.txt --b2 nf_bams.txt \
    --cov "$COV" --group1 HCM --group2 NF --min-samples 15 \
    --out ${ET}_HCM_v2

  echo
  echo "############ $ET : DCM-specificity ############"
  DCMTAB=${ET}_v2filtered.all_events.tsv
  [ -s "$DCMTAB" ] || { echo "missing $DCMTAB -- run 04_phase2_v2_checks.sbatch first"; continue; }
  python 05_hcm_specificity.py \
    --dcm "$DCMTAB" \
    --hcm ${ET}_HCM_v2.all_events.tsv \
    --out ${ET}_specificity
done

echo
echo "--- outputs ---"
ls -lh *_HCM_v2.all_events.tsv *_specificity.report.txt 2>/dev/null
date
