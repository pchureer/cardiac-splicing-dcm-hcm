#!/bin/bash
#SBATCH -p cpu
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH --qos=ngs16
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#
# 36 -- PHASE 1.5: the matched technical-replicate panel.
#
# THE OPPORTUNITY
# ---------------
# The genotype check found 24 LBNL samples that are the SAME INDIVIDUALS as MAGNet
# samples: 14 NF and 10 DCM, every pair matching on disease group. That cost the paper
# a replication cohort, but it created something nothing else in this project provides
# -- the same hearts measured twice, under two library preps and two read layouts:
#
#     MAGNet arm   poly-A,    PAIRED 50bp
#     LBNL arm     total-RNA, SINGLE 50bp
#
# Biology is held EXACTLY constant. Any attenuation between the two arms is technical.
#
# THE DECOMPOSITION
# -----------------
# Three slopes, all at the same group sizes (10 DCM vs 14 NF), so power is matched and
# the comparison is fair:
#
#   A) MAGNet_matched -> LBNL_matched    SAME people, DIFFERENT platform  = TECHNICAL
#   B) MAGNet_matched -> MAGNet_other    DIFFERENT people, SAME platform  = SAMPLING
#   C) MAGNet_full    -> Colorado        different people AND platform    = 0.521 (known)
#
# If A is near 1.0, technical differences do not attenuate and the 0.52 cross-cohort
# slope is biological plus sampling. If A is near 0.5, technical differences account for
# essentially all of it, and RESULT_magnitude_not_estimable.md becomes a quantified
# technical attenuation term rather than a shrug.
#
# SELECTION IS INDEPENDENT OF BOTH ARMS
# -------------------------------------
# Events are chosen from a MAGNet discovery run that EXCLUDES every sample used in any
# comparison arm. Selecting on the full MAGNet would bias the MAGNet-matched arm (whose
# samples are inside the selection set) more than the LBNL arm, manufacturing exactly
# the asymmetry this script exists to measure.
#
# NO COVARIATES ANYWHERE, DELIBERATELY. MAGNet has age/sex/RIN and LBNL has none.
# Adjusting one arm and not the other would confound platform with covariate
# adjustment. Symmetry matters more than adjustment here; say so in the methods.
#
# LBNL USES THE ORIGINAL BAM LISTS, NOT THE _dedup ONES -- the 24 twins are exactly the
# samples the dedup removed, so this analysis needs them back.
#
# Run:  bash 36_matched_pair_panel.sbatch 2>&1 | tail -60
date
export PATH=/data/apps/anaconda3/envs/mamba/envs/rnaseq_tools/bin:$PATH
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
ROOT=/data/project/chureeratp/RNA/HEART
cd "$ROOT"

TYPES=${TYPES:-SE}
MINS=${MINS:-8}
for f in geno/duplicates_triaged.tsv dcm_bams.txt nf_bams.txt \
         LBNL_dcm_bams.txt LBNL_nf_bams.txt magnet_covariates.tsv lbnl_dcm_nf.cov; do
  [ -s "$f" ] || { echo "*** missing $f"; exit 1; }
done

echo "############################################################"
echo "#####  0) build the matched sets and the exclusion lists"
echo "############################################################"
python3 - <<'PY' || exit 1
import pandas as pd

T = pd.read_csv("geno/duplicates_triaged.tsv", sep="\t")
C = T[(T.triage == "confirmed") & (T.cohort1 != T.cohort2)].copy()
pairs = []
for _, r in C.iterrows():
    if {r.cohort1, r.cohort2} != {"MAGNet", "LBNL"}:
        continue
    m, l = (r.s1, r.s2) if r.cohort1 == "MAGNet" else (r.s2, r.s1)
    g = r.group1 if r.cohort1 == "MAGNet" else r.group2
    pairs.append({"magnet": m, "lbnl": l, "group": str(g).upper()})
P = pd.DataFrame(pairs)
if P.empty:
    raise SystemExit("*** no confirmed MAGNet<->LBNL pairs in duplicates_triaged.tsv")
P.to_csv("matched_pairs.tsv", sep="\t", index=False)
print(f"[info] matched MAGNet<->LBNL pairs: {len(P)}")
print(P.group.value_counts().to_string())

nD = int((P.group == "DCM").sum())
nN = int((P.group == "NF").sum())
if nD < 5 or nN < 5:
    raise SystemExit(f"*** only {nD} DCM and {nN} NF pairs -- too few for this analysis")

# arm 1: the matched MAGNet samples
open("mp_magnet_matched.keep", "w").write("\n".join(P.magnet) + "\n")
# arm 2: the matched LBNL twins
open("mp_lbnl_matched.keep", "w").write("\n".join(P.lbnl) + "\n")

# arm 3: a same-size MAGNet subset drawn from NON-matched samples.
# Deterministic (sorted, then evenly spaced) so the run is reproducible.
def ids(path):
    raw = open(path).read()
    import re, os
    out = []
    for p in re.split(r"[,\n]", raw):
        p = p.strip()
        if not p:
            continue
        m = re.search(r"(SRR\d+|ERR\d+|SAMEA\d+|SAMN\d+)", p)
        out.append(m.group(1) if m else os.path.basename(os.path.dirname(p)))
    return sorted(set(out))

mdcm, mnf = ids("dcm_bams.txt"), ids("nf_bams.txt")
used = set(P.magnet)
freeD = [x for x in mdcm if x not in used]
freeN = [x for x in mnf if x not in used]
print(f"[info] MAGNet non-matched pool: {len(freeD)} DCM, {len(freeN)} NF")
if len(freeD) < nD * 2 or len(freeN) < nN * 2:
    raise SystemExit("*** not enough non-matched MAGNet samples for arm 3 + selection")
pick = lambda L, k: [L[round(i * (len(L) - 1) / (k - 1))] for i in range(k)]
otherD, otherN = pick(freeD, nD), pick(freeN, nN)
open("mp_magnet_other.keep", "w").write("\n".join(otherD + otherN) + "\n")

# selection set: everything NOT used in any arm
sel_excl = sorted(used | set(otherD) | set(otherN))
open("mp_selection.exclude", "w").write("\n".join(sel_excl) + "\n")
print(f"[info] selection run excludes {len(sel_excl)} samples; "
      f"leaves {len(freeD) - nD} DCM and {len(freeN) - nN} NF")

# phase2 has --exclude but no --include, so an arm's exclude file is
# "everything in that cohort that is not in the keep list"
def write_excl(name, all_ids, keep):
    drop = sorted(set(all_ids) - set(keep))
    open(name, "w").write("\n".join(drop) + "\n")
    print(f"[write] {name}  ({len(drop)} excluded, {len(set(keep) & set(all_ids))} kept)")

write_excl("mp_magnet_matched.exclude", mdcm + mnf, P.magnet)
write_excl("mp_magnet_other.exclude", mdcm + mnf, otherD + otherN)
ldcm, lnf = ids("LBNL_dcm_bams.txt"), ids("LBNL_nf_bams.txt")
write_excl("mp_lbnl_matched.exclude", ldcm + lnf, P.lbnl)
PY

NPAIR=$(( $(wc -l < matched_pairs.tsv) - 1 ))
echo "matched pairs: $NPAIR"

echo
echo "############################################################"
echo "#####  1) four contrasts, identical group sizes, no covariates"
echo "############################################################"
for ET in $TYPES; do
  echo; echo "----- $ET -----"

  echo "  [selection]  MAGNet, all arms excluded"
  python phase2_diffsplice_v2.py --mats rmats_magnet_rl125/${ET}.MATS.JCEC.txt --event-type $ET \
    --b1 dcm_bams.txt --b2 nf_bams.txt --cov magnet_covariates.tsv \
    --group1 DCM --group2 NF --min-samples 20 --no-covariates \
    --exclude mp_selection.exclude --out mp_${ET}_selection_rl125 2>&1 | grep -E "tested|sig FDR" | head -3

  echo "  [arm A1]     MAGNet, matched samples only"
  python phase2_diffsplice_v2.py --mats rmats_magnet_rl125/${ET}.MATS.JCEC.txt --event-type $ET \
    --b1 dcm_bams.txt --b2 nf_bams.txt --cov magnet_covariates.tsv \
    --group1 DCM --group2 NF --min-samples $MINS --no-covariates \
    --exclude mp_magnet_matched.exclude --out mp_${ET}_magnet_matched_rl125 2>&1 | grep -E "tested|sig FDR" | head -3

  echo "  [arm A2]     LBNL, the same individuals"
  python phase2_diffsplice_v2.py --mats rmats_LBNL/${ET}.MATS.JCEC.txt --event-type $ET \
    --b1 LBNL_dcm_bams.txt --b2 LBNL_nf_bams.txt --cov lbnl_dcm_nf.cov \
    --group1 DCM --group2 NF --min-samples $MINS --no-covariates \
    --exclude mp_lbnl_matched.exclude --out mp_${ET}_lbnl_matched_rl125 2>&1 | grep -E "tested|sig FDR" | head -3

  echo "  [arm B]      MAGNet, different people, same size"
  python phase2_diffsplice_v2.py --mats rmats_magnet_rl125/${ET}.MATS.JCEC.txt --event-type $ET \
    --b1 dcm_bams.txt --b2 nf_bams.txt --cov magnet_covariates.tsv \
    --group1 DCM --group2 NF --min-samples $MINS --no-covariates \
    --exclude mp_magnet_other.exclude --out mp_${ET}_magnet_other_rl125 2>&1 | grep -E "tested|sig FDR" | head -3
done

echo
echo "############################################################"
echo "#####  2) the decomposition"
echo "############################################################"
TYPES="$TYPES" TBL_SUFFIX=_rl125 python3 36_matched_pair_analysis.py
date
