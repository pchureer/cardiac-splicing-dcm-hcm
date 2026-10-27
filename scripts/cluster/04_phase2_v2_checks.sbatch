#!/bin/bash
#SBATCH -p cpu
#SBATCH -c 1
#SBATCH --mem=6G
#SBATCH --qos=ngs16
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#
# Three independent checks on the discovery analysis. Nothing here depends on
# the RBM20 expression work -- it can run at the same time.
#   A) does phase2_diffsplice_v2 --legacy reproduce your existing SE numbers?
#   B) do the results survive adding pool (batch) to the model?
#   C) do they survive dropping the young NF donors?
# Submit:  sbatch 04_phase2_v2_checks.sbatch

date
export PATH=/data/apps/anaconda3/envs/mamba/envs/rnaseq_tools/bin:$PATH
ROOT=/data/project/chureeratp/RNA/HEART
cd "$ROOT"

MATS=rmats_magnet_rl125/SE.MATS.JCEC.txt      # same file your phase-2 used
[ -s "$MATS" ] || { echo "missing $MATS"; exit 1; }

echo "############ A) LEGACY REPRODUCTION ############"
echo "compare 'sig FDR<0.05' below against your existing SE_adjusted table"
python phase2_diffsplice_v2.py \
  --mats "$MATS" --b1 dcm_bams.txt --b2 nf_bams.txt \
  --cov magnet_covariates.tsv --min-samples 20 \
  --legacy --out SE_v2legacy_rl125

echo
echo "############ A2) FILTERED (dedup + coverage + PSI-variance) ############"
echo "this is the locked-table version -- expect FEWER events, the drop is artifacts"
echo "run for every event type; these tables feed 05_hcm_specificity"
for ET in SE MXE A3SS A5SS RI; do
  M=rmats_magnet_rl125/${ET}.MATS.JCEC.txt
  [ -s "$M" ] || { echo "  missing $M -- skipping"; continue; }
  echo "----- $ET -----"
  python phase2_diffsplice_v2.py \
    --mats "$M" --event-type $ET --b1 dcm_bams.txt --b2 nf_bams.txt \
    --cov magnet_covariates.tsv --min-samples 20 \
    --out ${ET}_v2filtered_rl125
done

echo
echo "############ B) POOL-ADJUSTED SENSITIVITY ############"
echo "Magnet_03 is 21 DCM/7 NF and Magnet_10 is 9/19 -- does batch change anything?"
python phase2_diffsplice_v2.py \
  --mats "$MATS" --b1 dcm_bams.txt --b2 nf_bams.txt \
  --cov magnet_covariates.tsv --min-samples 20 \
  --adjust-pool --out SE_v2pooladj_rl125

echo
echo "############ C) AGE SENSITIVITY: drop donors under 25 ############"
# NF has 8 donors under 25 (4 under 20, youngest 15); DCM has 4 and none under 20.
# PSI is per-sample, so --exclude drops them from the existing rMATS output.
# NO rMATS RERUN NEEDED.
awk -F'\t' '$4+0 < 25 {print $1}' magnet_covariates.tsv > young_under25.ids
echo "donors under 25 excluded: $(wc -l < young_under25.ids)"
python phase2_diffsplice_v2.py \
  --mats "$MATS" --b1 dcm_bams.txt --b2 nf_bams.txt \
  --cov magnet_covariates.tsv --min-samples 20 \
  --exclude young_under25.ids --out SE_v2age25_rl125

echo
echo "--- outputs ---"
ls -lh SE_v2*.all_events.tsv 2>/dev/null
date
