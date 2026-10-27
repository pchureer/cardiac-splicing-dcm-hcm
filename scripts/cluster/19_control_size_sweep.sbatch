#!/bin/bash
#SBATCH -p cpu
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH --qos=ngs16
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#
# CONTROL-GROUP SIZE SWEEP -- the third explanation for the external attenuation.
#
# THE PROBLEM
# -----------
# 17 established that MAGNet's HCM slope is ~0.87 with 83 disjoint controls, while
# the three external HCM cohorts give 0.44-0.70. Three explanations are live:
#
#   (1) disease stage      end-stage explant vs preserved-EF myectomy
#   (2) read length        50bp vs 101/150bp changes IncFormLen normalisation
#   (3) CONTROL-GROUP SIZE MAGNet has 83 controls. Brigham has 5, Stanford 7,
#                          Fuwai 9. Noisier dPSI lowers the slope -- and it does
#                          so SPECIFICALLY, because the Deming lambda is built
#                          from se_adj, which omits the control-group MEAN
#                          sampling error. That omission is proportionally far
#                          worse at n=5 than at n=83, so Deming under-corrects
#                          more in exactly the cohorts with the low slopes.
#
# (3) is the most boring explanation and the cheapest to test, so it goes first.
# If it explains the gap, the "stage gradient" reading in RESULT_split_control_check
# has to be withdrawn.
#
# THE DESIGN
# ----------
# Hold everything constant except the number of controls. The DCM reference arm is
# fixed (DCM vs NF half-A, from script 16). The HCM arm is rebuilt against
# m controls drawn from NF half-B -- always disjoint from the DCM arm -- for
#
#     m = 5, 7, 9, 14, 28, 56, 83     (5 independent draws each)
#
# 5 / 7 / 9 are exactly Brigham / Stanford / Fuwai. Then plot slope against m.
#
#   * slope stays ~0.87 down to m=5  -> control size is NOT the explanation, and
#                                       the stage reading survives.
#   * slope falls to 0.44-0.70 at m=5-9 -> the external "attenuation" is a
#                                       small-control-group artefact. WITHDRAW the
#                                       stage-gradient claim. This would be the
#                                       single most important correction left.
#   * partial fall                    -> report the decomposition, as with 16.
#
# A case-count arm (k=13, matching Stanford) is included at m=7 to confirm that
# case count is not doing the work -- MAGNet HCM n=28 is already close to the
# external 13-28, so this should barely move.
#
# SE only, to bound runtime: ~40 phase-2 runs on the existing rMATS tables.
# No new alignment. Requires 16 to have run (uses splitctl/nfA_1.ids, nfB_1.ids).
#
# Submit: sbatch 19_control_size_sweep.sbatch
date
export PATH=/data/apps/anaconda3/envs/mamba/envs/rnaseq_tools/bin:$PATH
ROOT=/data/project/chureeratp/RNA/HEART
cd "$ROOT"

REPS=5
ET=SE
MATS_H=rmats_HCMvsNF/${ET}.MATS.JCEC.txt
DCMREF=splitctl/${ET}_DCMvsA_1.all_events.tsv

for f in "$MATS_H" "$DCMREF" splitctl/nfA_1.ids splitctl/nfB_1.ids hcm_bams.txt \
         nf_bams.txt hcm_nf_covariates.tsv; do
  [ -s "$f" ] || { echo "missing $f -- run 16 first"; exit 1; }
done
mkdir -p sizesweep

echo "NF half-A: $(grep -c . splitctl/nfA_1.ids)   half-B: $(grep -c . splitctl/nfB_1.ids)"
echo "DCM reference arm: $DCMREF  (fixed for the whole sweep)"

# HCM sample ids, same rule phase2 uses
tr ',' '\n' < hcm_bams.txt | grep . | \
  sed -E 's#.*(SRR[0-9]+|ERR[0-9]+|SAMEA[0-9]+|SAMN[0-9]+).*#\1#; t; s#.*/([^/]*)/[^/]*$#\1#' \
  > sizesweep/hcm_all.ids
echo "HCM cases: $(grep -c . sizesweep/hcm_all.ids)"

run_one () {   # $1=m controls  $2=k cases  $3=rep
  local m=$1 k=$2 r=$3
  local tag="m${m}_k${k}_r${r}"
  python3 - "$m" "$k" "$r" <<'PY'
import random, sys
m, k, r = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
rng = random.Random(90000 + m * 1000 + k * 100 + r)
A  = [l.strip() for l in open("splitctl/nfA_1.ids") if l.strip()]
B  = [l.strip() for l in open("splitctl/nfB_1.ids") if l.strip()]
H  = [l.strip() for l in open("sizesweep/hcm_all.ids") if l.strip()]
keep_ctrl = set(rng.sample(B, min(m, len(B))))
keep_case = set(rng.sample(H, min(k, len(H))))
# exclude: ALL of half A (keeps the arm disjoint from the DCM reference),
# the unselected controls in half B, and the unselected HCM cases
drop = set(A) | (set(B) - keep_ctrl) | (set(H) - keep_case)
open(f"sizesweep/excl_m{m}_k{k}_r{r}.ids", "w").write("\n".join(sorted(drop)) + "\n")
PY
  # min-samples must not exceed the smaller group
  local ms=$(( m < k ? m : k ))
  python phase2_diffsplice_v2.py --mats "$MATS_H" --event-type $ET \
    --b1 hcm_bams.txt --b2 nf_bams.txt --cov hcm_nf_covariates.tsv \
    --group1 HCM --group2 NF --min-samples $ms \
    --exclude sizesweep/excl_m${m}_k${k}_r${r}.ids \
    --out sizesweep/${ET}_${tag} > sizesweep/${ET}_${tag}.log 2>&1
  echo "  $tag  $(grep -c . sizesweep/${ET}_${tag}.all_events.tsv 2>/dev/null || echo 0) rows"
}

echo
echo "=== control-size sweep (cases fixed at 28) ==="
for m in 5 7 9 14 28 56 83; do
  for r in $(seq 1 $REPS); do run_one $m 28 $r; done
done

echo
echo "=== case-count control arm (m=7) ==="
for k in 13 18; do
  for r in $(seq 1 $REPS); do run_one 7 $k $r; done
done

echo
echo "############### SUMMARY ###############"
python3 - <<'PY'
import glob, os, re
import numpy as np, pandas as pd
from scipy import stats

PROGRAM = ["TTN","CAMK2D","LDB3","LMO7","SORBS1","PDLIM3","DTNA","MLIP","TNNT2"]
EXTERNAL = {"Brigham": (5, 0.662, 0.587), "Stanford": (7, 0.579, 0.703),
            "Fuwai": (9, 0.328, 0.440)}   # n_ctrl, background slope, program slope

def deming(x, y, lam):
    if len(x) < 3: return np.nan
    sxx, syy = np.var(x, ddof=1), np.var(y, ddof=1)
    sxy = np.cov(x, y, ddof=1)[0, 1]
    if abs(sxy) < 1e-12: return np.nan
    a = syy - lam * sxx
    return (a + np.sqrt(a*a + 4*lam*sxy*sxy)) / (2*sxy)

d = pd.read_csv("splitctl/SE_DCMvsA_1.all_events.tsv", sep="\t").drop_duplicates("key")
dsig = d[(d.FDR < 0.05) & (d.dPSI_adj.abs() >= 0.10)][["key","gene","dPSI_adj","se_adj"]]
dsig = dsig.rename(columns={"dPSI_adj":"x","se_adj":"sx"})
print(f"DCM reference: {len(dsig)} significant events\n")

rows = []
for f in sorted(glob.glob("sizesweep/SE_m*_k*_r*.all_events.tsv")):
    mm = re.search(r"m(\d+)_k(\d+)_r(\d+)", f)
    m, k, r = (int(g) for g in mm.groups())
    h = pd.read_csv(f, sep="\t").drop_duplicates("key")
    j = dsig.merge(h[["key","dPSI_adj","se_adj"]], on="key").dropna(
        subset=["x","dPSI_adj"])
    for label, sub in (("background", j), ("program", j[j.gene.isin(PROGRAM)])):
        if len(sub) < 8: continue
        x, y = sub.x.values, sub.dPSI_adj.values
        lam = np.nanmean(sub.se_adj.values**2) / max(np.nanmean(sub.sx.values**2), 1e-12)
        rows.append({"m": m, "k": k, "rep": r, "set": label, "n": len(sub),
                     "deming": deming(x, y, lam),
                     "r": stats.pearsonr(x, y)[0],
                     "pct_same": 100*np.mean(np.sign(x) == np.sign(y))})
R = pd.DataFrame(rows)
R.to_csv("control_size_sweep.tsv", sep="\t", index=False, float_format="%.6g")

for st in ("background", "program"):
    S = R[(R.set == st) & (R.k == 28)]
    if not len(S): continue
    print(f"--- {st}, cases fixed at 28 ---")
    print(f"  {'n_ctrl':>7}{'events':>8}{'Deming':>18}{'r':>16}{'%same':>16}")
    for m, g in S.groupby("m"):
        print(f"  {m:>7}{int(g.n.mean()):>8}"
              f"{f'{g.deming.mean():.3f} ± {g.deming.std(ddof=1):.3f}':>18}"
              f"{f'{g.r.mean():.3f}':>16}{f'{g.pct_same.mean():.1f}%':>16}")
    print()

print("--- case-count control (m=7) ---")
for st in ("background", "program"):
    S = R[(R.set == st) & (R.m == 7)]
    for k, g in S.groupby("k"):
        print(f"  {st:<11} cases={k:<4} Deming = {g.deming.mean():.3f} "
              f"± {g.deming.std(ddof=1):.3f}")
print()

print("--- THE COMPARISON THAT DECIDES IT ---")
print(f"  {'cohort':<10}{'n_ctrl':>7}{'observed bg':>13}{'observed prog':>15}"
      f"{'MAGNet at same n_ctrl (bg)':>29}")
for c, (n, bg, pg) in sorted(EXTERNAL.items(), key=lambda t: t[1][0]):
    S = R[(R.set == "background") & (R.k == 28) & (R.m == n)]
    mag = f"{S.deming.mean():.3f} ± {S.deming.std(ddof=1):.3f}" if len(S) else "--"
    print(f"  {c:<10}{n:>7}{bg:>13.3f}{pg:>15.3f}{mag:>29}")
print()
print("READ IT LIKE THIS:")
print("  If MAGNet's slope at n_ctrl=5-9 is still ~0.87, small control groups do")
print("  NOT explain the external cohorts and the disease-stage reading survives.")
print("  If MAGNet's slope at n_ctrl=5-9 has fallen into 0.44-0.70, the external")
print("  'attenuation' is a control-group-size artefact. The stage-gradient claim")
print("  in RESULT_split_control_check.md must then be WITHDRAWN -- the cohorts")
print("  would be telling us about their own noise, not about disease stage.")
print("  Also check %same: it should track the slope. If concordance holds up while")
print("  the slope collapses, the effect is specifically the errors-in-variables")
print("  under-correction, not a loss of real signal.")
PY
date
