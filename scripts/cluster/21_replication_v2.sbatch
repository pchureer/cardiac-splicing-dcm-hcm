#!/bin/bash
#SBATCH -p cpu
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH --qos=ngs16
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#
# Colorado + LBNL replication, redone on v2 tables.
#
# WHY: the replication numbers in the handoff (Colorado SE 22/28 same-direction,
# LBNL SE 16/25) came from the v1 pipeline -- no coordinate dedup, no coverage mask,
# no PSI-variance filter, no `key` column. Multi-transcript rMATS rows duplicate the
# same event, so a v1 "28 events" may be fewer real events counted more than once, and
# constant-PSI events gave nan correlations. Everything else in this project is on v2
# now; these two have to match before they share a table.
#
# Colorado DCM-vs-NF v2 tables already exist from script 12
# (${ET}_Colorado_DCMvsNF.all_events.tsv). They are reused unless --force-colorado
# is set below. LBNL has never been run through v2 at all.
#
# SAFE TO RUN ALONGSIDE script 20 -- different rMATS directories, different output
# names, no shared files.
#
# COHORTS (both SINGLE-end, both ~50bp like MAGNet -- this matters, see the note at
# the end of the summary):
#   Colorado  PRJNA477855 / GSE116250   37 DCM vs 14 NF   SINGLE 51bp  intermediate prep
#   LBNL      PRJNA522417               15 DCM vs 18 NF   SINGLE 50bp  total-RNA
#
# Submit: sbatch 21_replication_v2.sbatch
date
export PATH=/data/apps/anaconda3/envs/mamba/envs/rnaseq_tools/bin:$PATH
ROOT=/data/project/chureeratp/RNA/HEART

# --------------------------------------------------------------------------
# Read length. Called RLSUF, NOT SUF: 21_replication_v2.py's --suffix means the
# COHORT (_Colorado_DCMvsNF), and a collision would mis-name every table.
#
# Only the MAGNet DISCOVERY table carries this suffix. rmats_LBNL/ and
# rmats_Colorado/ below are deliberately unsuffixed — neither cohort was
# re-quantified, and a 125 bp discovery set against 50 bp cohort tables is what
# cross-cohort replication means here.
RLSUF="${TBL_SUFFIX:-}"

# Mirrors FORCE_COLORADO. Without it the LBNL loop re-runs phase2 five times on
# every invocation and rewrites ${ET}_LBNL_DCMvsNF.all_events.tsv for nothing.
FORCE_LBNL="${FORCE_LBNL:-0}"
# --------------------------------------------------------------------------
cd "$ROOT"

FORCE_COLORADO=0          # set to 1 to rebuild Colorado's v2 tables from scratch

# ---------------------------------------------------------------- LBNL bam lists
# Build them if they do not exist. Accepts either naming convention.
find_bamlist () {   # $1 = cohort, $2 = group  -> echoes a filename or nothing
  for f in ${1}_${2}_bams.txt ${2}_bams_${1}.txt ${1}_${2}.bams.txt; do
    [ -s "$f" ] && { echo "$f"; return; }
  done
}

for GRP in dcm nf; do
  F=$(find_bamlist LBNL $GRP)
  if [ -z "$F" ]; then
    IDS=""
    case $GRP in
      dcm) IDS=DCM_LBNL.ids ;;
      nf)  IDS=NF_LBNL.ids  ;;
    esac
    if [ ! -s "$IDS" ]; then
      echo "MISSING $IDS -- cannot build the LBNL bam list. Stop and check the ids."
      exit 1
    fi
    echo "building LBNL_${GRP}_bams.txt from $IDS"
    for s in $(grep . "$IDS"); do
      ls "$ROOT/$s/"*Aligned.sortedByCoord.out.bam 2>/dev/null | head -1
    done | paste -sd, > LBNL_${GRP}_bams.txt
  fi
done
echo "LBNL DCM bams: $(tr ',' '\n' < LBNL_dcm_bams.txt | grep -c .)  (expect 15)"
echo "LBNL NF  bams: $(tr ',' '\n' < LBNL_nf_bams.txt  | grep -c .)  (expect 18)"

# !! LBNL also contains 5 FETAL hearts. They must not be in NF_LBNL.ids. Check:
echo
echo "CHECK THIS: if NF_LBNL.ids is 23 rather than 18, the 5 fetal hearts are in the"
echo "control group and every LBNL number below is wrong. NF_LBNL.ids has"
echo "$(grep -c . NF_LBNL.ids 2>/dev/null || echo '?') entries."

# ---------------------------------------------------------------- covariates
# Neither cohort has age/sex/RIN, so group labels only and --no-covariates.
printf '%s\tDCM\n' $(cat DCM_LBNL.ids) >  lbnl_dcm_nf.cov
printf '%s\tNF\n'  $(cat NF_LBNL.ids)  >> lbnl_dcm_nf.cov
wc -l lbnl_dcm_nf.cov

if [ ! -s colorado_dcm_nf.cov ] || [ "$FORCE_COLORADO" = "1" ]; then
  printf '%s\tDCM\n' $(cat DCM_Colorado.ids) >  colorado_dcm_nf.cov
  printf '%s\tNF\n'  $(cat NF_Colorado.ids)  >> colorado_dcm_nf.cov
fi

# ---------------------------------------------------------------- phase 2
echo
echo "############################################################"
echo "#####  LBNL : DCM vs NF, v2, all five event types"
echo "############################################################"
for ET in SE MXE A3SS A5SS RI; do
  OUT=${ET}_LBNL_DCMvsNF.all_events.tsv
  if [ -s "$OUT" ] && [ "$FORCE_LBNL" != "1" ]; then
    echo "----- $ET : reusing $OUT (FORCE_LBNL=1 to rebuild) -----"
    continue
  fi
  M=rmats_LBNL/${ET}.MATS.JCEC.txt
  [ -s "$M" ] || { echo "missing $M -- skipping $ET"; continue; }
  echo "----- $ET -----"
  # min-samples 8: limited by the smaller group (15 DCM), kept permissive because
  # LBNL is the smallest cohort in the project.
  python phase2_diffsplice_v2.py \
    --mats "$M" --event-type $ET \
    --b1 LBNL_dcm_bams.txt --b2 LBNL_nf_bams.txt \
    --cov lbnl_dcm_nf.cov --group1 DCM --group2 NF \
    --min-samples 8 --no-covariates \
    --out ${ET}_LBNL_DCMvsNF
done

echo
echo "############################################################"
echo "#####  Colorado : DCM vs NF, v2"
echo "############################################################"
for ET in SE MXE A3SS A5SS RI; do
  OUT=${ET}_Colorado_DCMvsNF.all_events.tsv
  if [ -s "$OUT" ] && [ "$FORCE_COLORADO" != "1" ]; then
    echo "----- $ET : reusing $OUT (from script 12) -----"
    continue
  fi
  M=rmats_Colorado/${ET}.MATS.JCEC.txt
  [ -s "$M" ] || { echo "missing $M -- skipping $ET"; continue; }
  echo "----- $ET -----"
  python phase2_diffsplice_v2.py \
    --mats "$M" --event-type $ET \
    --b1 Colorado_dcm_bams.txt --b2 Colorado_nf_bams.txt \
    --cov colorado_dcm_nf.cov --group1 DCM --group2 NF \
    --min-samples 10 --no-covariates \
    --out ${ET}_Colorado_DCMvsNF
done

# ---------------------------------------------------------------- concordance
echo
echo "############################################################"
echo "#####  CONCORDANCE AGAINST THE MAGNet DISCOVERY TABLES"
echo "############################################################"
for ET in SE MXE A3SS A5SS RI; do
  [ -s "${ET}_v2filtered${RLSUF}.all_events.tsv" ] || \
    echo "WARNING: ${ET}_v2filtered${RLSUF}.all_events.tsv missing -- run 04 first"
done

python 21_replication_v2.py --rep Colorado --suffix _Colorado_DCMvsNF \
    --disc-suffix _v2filtered${RLSUF} \
    --out replication_v2_Colorado${RLSUF}
echo
python 21_replication_v2.py --rep LBNL     --suffix _LBNL_DCMvsNF \
    --disc-suffix _v2filtered${RLSUF} \
    --out replication_v2_LBNL${RLSUF}

echo
echo "############################################################"
echo "NOTE FOR THE DISCUSSION -- this doubles as a read-length control."
echo "Colorado (51bp) and LBNL (50bp) match MAGNet's read length, while the three"
echo "HCM cohorts are 101/150bp. If script 20 shows read length IS driving the"
echo "0.44-0.70 slopes, then these two DCM cohorts should come in NEAR 0.9. If they"
echo "also come in low, read length is not the explanation and something about"
echo "cross-cohort replication generally is attenuating the slope -- which would"
echo "affect how the HCM result is described too."
echo "############################################################"
date
