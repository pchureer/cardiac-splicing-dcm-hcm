#!/bin/bash
#SBATCH -p cpu
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH --qos=ngs16
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#
# SPLIT-CONTROL CHECK -- the missing control for the 15b result.
#
# THE PROBLEM THIS SOLVES
# -----------------------
# 15b found the MAGNet DCM signature transfers to three independent HCM cohorts
# at a Deming slope of 0.44-0.70, versus ~0.97 for MAGNet's own DCM-vs-HCM
# comparison. Two explanations are confounded:
#
#   (a) MAGNet's ~0.97 is INFLATED because both contrasts divide by the SAME
#       166 NF hearts, so every NF sampling fluctuation appears in both dPSIs.
#   (b) The external slope is genuinely LOWER because those cohorts are septal
#       myectomy tissue with preserved EF, not end-stage explants.
#
# You cannot tell these apart from 15b alone, and the difference matters: (a)
# means MAGNet's number was an artefact, (b) means it was real and the gap is
# disease stage. This script separates them WITHOUT any new rMATS run.
#
# THE DESIGN
# ----------
# Split the 166 MAGNet NF hearts at random into halves A and B, then build two
# versions of the same DCM-vs-HCM comparison:
#
#   SHARED   : DCM vs NF(A)   and   HCM vs NF(A)     -- same controls, n=83
#   DISJOINT : DCM vs NF(A)   and   HCM vs NF(B)     -- NO shared samples
#
# Both arms have identical control counts, identical cases, identical tissue,
# identical read length and annotation. The ONLY thing that differs is whether
# the two contrasts share control samples. So:
#
#   slope(SHARED) - slope(DISJOINT)  =  the shared-control inflation, measured.
#
# Then compare slope(DISJOINT) against the 15b external slopes (0.44-0.70):
#   * DISJOINT lands near 0.9-1.0  -> the gap to the external cohorts is REAL
#                                     biology (disease stage / myectomy vs
#                                     explant). Explanation (b).
#   * DISJOINT lands near 0.5-0.7  -> MAGNet's ~0.97 was mostly the shared-
#                                     control artefact, and the external
#                                     cohorts actually AGREE with a properly
#                                     controlled MAGNet. Explanation (a).
#   * DISJOINT lands in between    -> both contribute; report the decomposition.
#
# Repeated over NSPLIT random splits so the answer is not one lucky partition.
# The seed is fixed, so this is reproducible; record it in the methods.
#
# Cheap: reuses rmats_magnet/ and rmats_HCMvsNF/, only phase 2 is re-run.
#
# Submit: sbatch 16_split_control_check.sbatch
date
export PATH=/data/apps/anaconda3/envs/mamba/envs/rnaseq_tools/bin:$PATH
ROOT=/data/project/chureeratp/RNA/HEART
cd "$ROOT"

NSPLIT=5
TYPES="SE MXE"          # the two types the MAGNet HCM run was phase-2'd for

for f in magnet_covariates.tsv hcm_nf_covariates.tsv nf_bams.txt dcm_bams.txt hcm_bams.txt; do
  [ -s "$f" ] || { echo "missing $f"; exit 1; }
done

# NF sample ids, in the order they appear in nf_bams.txt
tr ',' '\n' < nf_bams.txt | grep . | sed 's#.*/\([^/]*\)/[^/]*$#\1#' > nf_all.ids
NNF=$(grep -c . nf_all.ids)
echo "NF samples: $NNF   (expect 166)"

mkdir -p splitctl

for S in $(seq 1 $NSPLIT); do
  echo
  echo "############################################################"
  echo "#####  SPLIT $S"
  echo "############################################################"

  python3 - "$S" <<'PY'
import random, sys
seed = 1000 + int(sys.argv[1])
ids = [l.strip() for l in open("nf_all.ids") if l.strip()]
random.Random(seed).shuffle(ids)
h = len(ids) // 2
open(f"splitctl/nfA_{sys.argv[1]}.ids", "w").write("\n".join(ids[:h]) + "\n")
open(f"splitctl/nfB_{sys.argv[1]}.ids", "w").write("\n".join(ids[h:]) + "\n")
print(f"  seed {seed}: A={h}  B={len(ids)-h}")
PY

  A=splitctl/nfA_${S}.ids
  B=splitctl/nfB_${S}.ids

  for ET in $TYPES; do
    MD=rmats_magnet/${ET}.MATS.JCEC.txt
    MH=rmats_HCMvsNF/${ET}.MATS.JCEC.txt
    [ -s "$MD" ] && [ -s "$MH" ] || { echo "  missing $ET MATS -- skip"; continue; }

    # DCM vs NF(A)  -- drop half B from the controls
    python phase2_diffsplice_v2.py --mats "$MD" --event-type $ET \
      --b1 dcm_bams.txt --b2 nf_bams.txt --cov magnet_covariates.tsv \
      --group1 DCM --group2 NF --min-samples 15 --exclude "$B" \
      --out splitctl/${ET}_DCMvsA_${S} > splitctl/${ET}_DCMvsA_${S}.log 2>&1

    # HCM vs NF(A)  -- SHARED controls
    python phase2_diffsplice_v2.py --mats "$MH" --event-type $ET \
      --b1 hcm_bams.txt --b2 nf_bams.txt --cov hcm_nf_covariates.tsv \
      --group1 HCM --group2 NF --min-samples 15 --exclude "$B" \
      --out splitctl/${ET}_HCMvsA_${S} > splitctl/${ET}_HCMvsA_${S}.log 2>&1

    # HCM vs NF(B)  -- DISJOINT controls
    python phase2_diffsplice_v2.py --mats "$MH" --event-type $ET \
      --b1 hcm_bams.txt --b2 nf_bams.txt --cov hcm_nf_covariates.tsv \
      --group1 HCM --group2 NF --min-samples 15 --exclude "$A" \
      --out splitctl/${ET}_HCMvsB_${S} > splitctl/${ET}_HCMvsB_${S}.log 2>&1

    D=splitctl/${ET}_DCMvsA_${S}.all_events.tsv
    [ -s "$D" ] || { echo "  $ET split $S: DCM table empty -- see the .log"; continue; }

    # SHARED arm: the standard shared-control baseline applies, so let 05
    # compute it (do NOT pass --no-shared-controls here).
    python 05_hcm_specificity.py --dcm "$D" \
      --hcm splitctl/${ET}_HCMvsA_${S}.all_events.tsv \
      --out splitctl/${ET}_SHARED_${S} > splitctl/${ET}_SHARED_${S}.log 2>&1

    # DISJOINT arm: no shared samples, so the null baseline is 0.
    python 05_hcm_specificity.py --dcm "$D" \
      --hcm splitctl/${ET}_HCMvsB_${S}.all_events.tsv --no-shared-controls \
      --out splitctl/${ET}_DISJOINT_${S} > splitctl/${ET}_DISJOINT_${S}.log 2>&1
  done
done

echo
echo "############### SUMMARY ###############"
echo "(n / r / Deming / sameDir are the RBM20-PROGRAM numbers, i.e. the same"
echo " fields the 15b summary reported, so they are directly comparable.)"
printf "%-5s %-9s %-6s %8s %9s %9s %10s\n" type arm split n r Deming sameDir
for ET in $TYPES; do
  for ARM in SHARED DISJOINT; do
    for S in $(seq 1 $NSPLIT); do
      R=splitctl/${ET}_${ARM}_${S}.report.txt
      [ -s "$R" ] || continue
      n=$(grep -m1 "measured in the HCM run as well" "$R" | sed 's/.*well: //; s/ .*//')
      r=$(grep -m1 "Pearson r = " "$R" | sed 's/.*Pearson r = //; s/ .*//')
      d=$(grep -m1 "Deming slope" "$R" | sed 's/.*Deming slope = //; s/ .*//')
      sd=$(grep -m1 "same direction:" "$R" | sed 's/.*same direction: //; s/ .*//')
      printf "%-5s %-9s %-6s %8s %9s %9s %10s\n" "$ET" "$ARM" "$S" "${n:-NA}" "${r:-NA}" "${d:-NA}" "${sd:-NA}"
    done
  done
done

echo
echo "HOW TO READ IT"
echo "  mean Deming(SHARED) - mean Deming(DISJOINT) = the shared-control inflation."
echo "  Compare mean Deming(DISJOINT) with the 15b external slopes:"
echo "     Brigham 0.587 | Fuwai 0.440 | Stanford 0.703   (SE, program events)"
echo "  If DISJOINT sits in that range, MAGNet and the external cohorts agree once"
echo "  the shared reference is removed, and the 'attenuation' in 15b is largely"
echo "  the artefact -- not disease stage. If DISJOINT stays near 1.0, the gap is"
echo "  real and is the myectomy/preserved-EF vs end-stage-explant difference."
echo
echo "CAVEAT to state either way: events entering the comparison are selected on"
echo "DCM significance, so regression to the mean pulls BOTH arms' slopes down by"
echo "the same amount. The SHARED-minus-DISJOINT difference is unaffected by it;"
echo "the absolute slopes are lower bounds."
date
