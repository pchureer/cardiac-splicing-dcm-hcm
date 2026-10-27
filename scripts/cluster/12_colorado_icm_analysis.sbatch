#!/bin/bash
#SBATCH -p cpu
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH --qos=ngs16
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#
# Colorado: independent replication of etiology-independence, second cohort,
# second etiology pair (DCM vs ICM rather than DCM vs HCM).
#
# Three contrasts, all five event types, all --no-covariates (Colorado has no
# age/sex/RIN):
#   A) DCM vs NF   (from the existing rmats_Colorado run)
#   B) ICM vs NF   (from 09b)
#   C) DCM vs ICM  DIRECTLY (from 09b) -- no shared reference, the clean test
#
# Then:
#   - A vs B comparison (shares the 14 NF hearts, so the shared-control
#     baseline applies -- script 05 computes it: with 37/13/14 it will be
#     LARGER than the 0.27 seen in MAGNet, because 14 controls are noisier)
#   - equivalence: how big is the DCM-ICM difference relative to DCM-NF?
#
# POWER WARNING: ICM n=13, NF n=14. This is a far weaker test than MAGNet.
# A null here is much softer evidence than the DCM-vs-HCM null. Read the
# mde80 / powered-event numbers, not the p-values.
#
# Submit: sbatch 12_colorado_icm_analysis.sbatch

date
export PATH=/data/apps/anaconda3/envs/mamba/envs/rnaseq_tools/bin:$PATH
ROOT=/data/project/chureeratp/RNA/HEART
cd "$ROOT"

for f in Colorado_dcm_bams.txt Colorado_icm_bams.txt Colorado_nf_bams.txt; do
  [ -s "$f" ] || { echo "missing $f -- run 09b first"; exit 1; }
done

echo "=== covariate files (group labels only) ==="
printf '%s\tDCM\n' $(cat DCM_Colorado.ids) >  colorado_dcm_icm.cov
printf '%s\tICM\n' $(cat ICM_Colorado.ids) >> colorado_dcm_icm.cov
printf '%s\tICM\n' $(cat ICM_Colorado.ids) >  colorado_icm_nf.cov
printf '%s\tNF\n'  $(cat NF_Colorado.ids)  >> colorado_icm_nf.cov
printf '%s\tDCM\n' $(cat DCM_Colorado.ids) >  colorado_dcm_nf.cov
printf '%s\tNF\n'  $(cat NF_Colorado.ids)  >> colorado_dcm_nf.cov
wc -l colorado_*.cov

for ET in SE MXE A3SS A5SS RI; do
  echo
  echo "############################################################"
  echo "#####                 $ET                              #####"
  echo "############################################################"

  # ---- A) DCM vs NF (existing replication run) ----------------------
  M=rmats_Colorado/${ET}.MATS.JCEC.txt
  if [ -s "$M" ]; then
    echo; echo "----- $ET Colorado DCM vs NF -----"
    python phase2_diffsplice_v2.py --mats "$M" --event-type $ET \
      --b1 Colorado_dcm_bams.txt --b2 Colorado_nf_bams.txt \
      --cov colorado_dcm_nf.cov --group1 DCM --group2 NF \
      --min-samples 10 --no-covariates --out ${ET}_Colorado_DCMvsNF
  else
    echo "missing $M -- skipping DCM-vs-NF for $ET"
  fi

  # ---- B) ICM vs NF -------------------------------------------------
  M=rmats_Colorado_ICMvsNF/${ET}.MATS.JCEC.txt
  if [ -s "$M" ]; then
    echo; echo "----- $ET Colorado ICM vs NF -----"
    python phase2_diffsplice_v2.py --mats "$M" --event-type $ET \
      --b1 Colorado_icm_bams.txt --b2 Colorado_nf_bams.txt \
      --cov colorado_icm_nf.cov --group1 ICM --group2 NF \
      --min-samples 10 --no-covariates --out ${ET}_Colorado_ICMvsNF
  fi

  # ---- C) DCM vs ICM directly ---------------------------------------
  M=rmats_Colorado_DCMvsICM/${ET}.MATS.JCEC.txt
  if [ -s "$M" ]; then
    echo; echo "----- $ET Colorado DCM vs ICM (DIRECT) -----"
    python phase2_diffsplice_v2.py --mats "$M" --event-type $ET \
      --b1 Colorado_dcm_bams.txt --b2 Colorado_icm_bams.txt \
      --cov colorado_dcm_icm.cov --group1 DCM --group2 ICM \
      --min-samples 10 --no-covariates --out ${ET}_Colorado_DCMvsICM
  fi

  # ---- effect-size agreement, DCM-vs-NF against ICM-vs-NF ------------
  if [ -s "${ET}_Colorado_DCMvsNF.all_events.tsv" ] && \
     [ -s "${ET}_Colorado_ICMvsNF.all_events.tsv" ]; then
    echo; echo "----- $ET Colorado: does the DCM signature appear in ICM? -----"
    python 05_hcm_specificity.py \
      --dcm ${ET}_Colorado_DCMvsNF.all_events.tsv \
      --hcm ${ET}_Colorado_ICMvsNF.all_events.tsv \
      --out ${ET}_Colorado_ICMspecificity
  fi
done

echo
echo "############### EQUIVALENCE: DCM-vs-ICM against DCM-vs-NF ###############"
python - <<'PY'
import glob, os
import numpy as np, pandas as pd

print(f"{'type':<6}{'DCMvNF sig':>11}{'DCMvICM sig':>12}{'median ratio':>14}"
      f"{'powered':>9}{'<half':>8}")
for ET in ["SE", "MXE", "A3SS", "A5SS", "RI"]:
    fa, fc = f"{ET}_Colorado_DCMvsNF.all_events.tsv", f"{ET}_Colorado_DCMvsICM.all_events.tsv"
    if not (os.path.exists(fa) and os.path.exists(fc)):
        continue
    a = pd.read_csv(fa, sep="\t").drop_duplicates("key")
    c = pd.read_csv(fc, sep="\t").drop_duplicates("key")
    m = a.merge(c, on="key", suffixes=("_nf", "_icm"))
    sig = m[(m.FDR_nf < 0.05) & (m.dPSI_adj_nf.abs() >= 0.10)]
    nsig_c = int(((c.FDR < 0.05) & (c.dPSI_adj.abs() >= 0.10)).sum())
    if len(sig) < 5:
        print(f"{ET:<6}{len(sig):>11}{nsig_c:>12}{'--':>14}{'--':>9}{'--':>8}")
        continue
    ratio = (sig.dPSI_adj_icm.abs() / sig.dPSI_adj_nf.abs())
    powered = sig[sig.mde80_icm <= sig.dPSI_adj_nf.abs()]
    small = powered[powered.dPSI_adj_icm.abs() < 0.5 * powered.dPSI_adj_nf.abs()]
    print(f"{ET:<6}{len(sig):>11}{nsig_c:>12}{np.median(ratio):>14.3f}"
          f"{f'{len(powered)}/{len(sig)}':>9}"
          f"{f'{100*len(small)/max(len(powered),1):.0f}%':>8}")
print()
print("median ratio near 0 -> DCM and ICM are indistinguishable on the events")
print("that define the DCM signature: etiology-independence replicates.")
print("'powered' = events where this comparison COULD have detected a")
print("DCM-vs-NF-sized difference. With ICM n=13 that fraction will be low;")
print("if it is small, the null is weak evidence and must be reported as such.")
PY
date
