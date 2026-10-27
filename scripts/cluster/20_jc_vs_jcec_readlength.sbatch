#!/bin/bash
#SBATCH -p cpu
#SBATCH -c 4
#SBATCH --mem=8G          # was 32G; phase 2's measured MaxRSS is
                          # 1.57 GB (job 408099), a 5x margin here
#SBATCH --qos=ngs16
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#
# READ-LENGTH CHECK -- the last explanation standing.
#
# THE PROBLEM
# -----------
# The three external HCM cohorts recapitulate the MAGNet DCM signature at a Deming
# slope of 0.44-0.70, while MAGNet's own disjoint-control comparison gives 0.87-0.99.
# The shared-control artefact is measured (+0.11, script 16/17) and control-group size
# is excluded (script 19). What is left is that MAGNet is 50bp and the external
# cohorts are 101/150bp.
#
# WHY READ LENGTH COULD DO IT
# ---------------------------
# rMATS turns read counts into PSI by dividing by effective isoform lengths:
#
#     PSI = (I / IncFormLen) / (I / IncFormLen + S / SkipFormLen)
#
# Both lengths depend on read length. The correction is meant to make PSI comparable
# across runs, but it is imperfect -- and any systematic compression of dPSI in the
# longer-read cohorts would look exactly like biological attenuation.
#
# THE TEST
# --------
# rMATS already wrote TWO count tables per event type:
#
#     JCEC  junction reads PLUS reads on the exon body   <- what we have used
#     JC    junction-spanning reads ONLY
#
# In the JC table the exon body does not contribute, so the length normalisation is
# far simpler and the read-length sensitivity is much reduced. Re-running the whole
# comparison on JC therefore separates the two:
#
#   * slopes stay at 0.44-0.70  -> read length is NOT the cause. Every technical
#                                  explanation is now excluded and the disease-stage
#                                  reading can be stated.
#   * slopes jump toward 0.9    -> read length WAS the cause. The stage-gradient claim
#                                  must be withdrawn and panel B relabelled.
#   * slopes move partway       -> report the decomposition.
#
# No new rMATS run and no realignment: the JC tables are already on disk.
#
# WHAT IS COMPARED
# ----------------
# The full ladder, in JC and in JCEC, on matched parameters:
#     MAGNet HCM vs NF(B)   against MAGNet DCM vs NF(A)   [disjoint, from script 16]
#     Brigham / Fuwai / Stanford HCM vs own controls      against MAGNet DCM vs NF
#
# The summary reports each slope three ways: on each table's own event set, and on the
# INTERSECTION of events measurable in both JC and JCEC. The intersection is the
# comparison that matters -- JC has fewer reads per event, so it loses events to the
# coverage filter, and a slope change caused by a different event set would otherwise
# be mistaken for a read-length effect.
#
# Submit: sbatch 20_jc_vs_jcec_readlength.sbatch
date
export PATH=/data/apps/anaconda3/envs/mamba/envs/rnaseq_tools/bin:$PATH
ROOT=/data/project/chureeratp/RNA/HEART
cd "$ROOT"

ET=SE
declare -A MINS=( [Brigham]=4 [Fuwai]=7 [Stanford]=5 )

# ---------------------------------------------------------------- read length
# Only the two MAGNet rMATS directories move. The three cohort directories
# (rmats_Brigham/ etc.) do NOT -- no replication cohort was re-quantified --
# which is why part C below must come out byte-identical at both read lengths,
# and why the run ends by checking exactly that.
SUF=${TBL_SUFFIX:-}
MAGDIR=rmats_magnet${SUF}
HCMDIR=rmats_HCMvsNF${SUF}
JCDIR=jccheck${SUF}
FORCE=${FORCE:-0}
echo "read length suffix : '${SUF}'   (empty = 50 bp)"
echo "MAGNet DCM-vs-NF   : ${MAGDIR}"
echo "MAGNet HCM-vs-NF   : ${HCMDIR}"
echo "output directory   : ${JCDIR}"
echo "cohort rMATS dirs  : rmats_<COH>/   (unsuffixed, deliberately)"

for f in ${MAGDIR}/${ET}.MATS.JC.txt magnet_covariates.tsv dcm_bams.txt nf_bams.txt; do
  [ -s "$f" ] || { echo "missing $f"; exit 1; }
done

# Never overwrite an existing set. jccheck/ is the ONLY copy of the 50 bp JC
# quantification and this script is the only thing that makes it.
if [ -n "$(ls -A "$JCDIR" 2>/dev/null)" ] && [ "$FORCE" != "1" ]; then
  echo
  echo "REFUSING — ${JCDIR}/ already contains files:"
  ls -1 "$JCDIR" | sed 's/^/    /'
  echo "  Re-run with FORCE=1 only after checking nothing downstream has"
  echo "  consumed them. S5(a) reads this directory."
  exit 1
fi
mkdir -p "$JCDIR"

echo "=== A) MAGNet DCM vs NF, JC table (matches the min-samples 20 used in 04) ==="
python phase2_diffsplice_v2.py \
  --mats ${MAGDIR}/${ET}.MATS.JC.txt --event-type $ET \
  --b1 dcm_bams.txt --b2 nf_bams.txt --cov magnet_covariates.tsv \
  --min-samples 20 --out ${JCDIR}/${ET}_v2filtered_JC

echo
echo "=== B) MAGNet disjoint arm, JC  (needs script 16's NF halves) ==="
if [ -s splitctl/nfA_1.ids ] && [ -s splitctl/nfB_1.ids ] && \
   [ -s ${HCMDIR}/${ET}.MATS.JC.txt ]; then
  python phase2_diffsplice_v2.py \
    --mats ${MAGDIR}/${ET}.MATS.JC.txt --event-type $ET \
    --b1 dcm_bams.txt --b2 nf_bams.txt --cov magnet_covariates.tsv \
    --group1 DCM --group2 NF --min-samples 15 --exclude splitctl/nfB_1.ids \
    --out ${JCDIR}/${ET}_DCMvsA_JC
  python phase2_diffsplice_v2.py \
    --mats ${HCMDIR}/${ET}.MATS.JC.txt --event-type $ET \
    --b1 hcm_bams.txt --b2 nf_bams.txt --cov hcm_nf_covariates.tsv \
    --group1 HCM --group2 NF --min-samples 15 --exclude splitctl/nfA_1.ids \
    --out ${JCDIR}/${ET}_HCMvsB_JC
else
  echo "  skipping the MAGNet disjoint arm (run 16 first, or no JC table for HCMvsNF)"
fi

echo
echo "=== C) the three external cohorts, JC ==="
for COH in Brigham Fuwai Stanford; do
  M=rmats_${COH}/${ET}.MATS.JC.txt
  [ -s "$M" ] || { echo "  missing $M -- skipping $COH"; continue; }
  echo "----- $COH -----"
  python phase2_diffsplice_v2.py \
    --mats "$M" --event-type $ET \
    --b1 ${COH}_hcm_bams.txt --b2 ${COH}_ctrl_bams.txt \
    --cov ${COH}.cov --group1 HCM --group2 CTRL \
    --min-samples ${MINS[$COH]} --no-covariates \
    --out ${JCDIR}/${COH}_${ET}_HCM_JC
done

# ---------------------------------------------------------------- control
# Part C reads ONLY the cohorts' own rMATS, which were never re-quantified,
# so at any read length these three files must be byte-identical to the 50 bp
# ones. Free, and a real check on phase-2 reproducibility.
if [ -n "$SUF" ]; then
  echo
  echo "=== CONTROL: part C must be byte-identical to the 50 bp run ==="
  _cfail=0
  for COH in Brigham Fuwai Stanford; do
    _new=${JCDIR}/${COH}_${ET}_HCM_JC.all_events.tsv
    _old=jccheck/${COH}_${ET}_HCM_JC.all_events.tsv
    if [ ! -s "$_old" ]; then
      echo "  [skip] $COH -- no 50 bp counterpart at $_old"
    elif cmp -s "$_new" "$_old"; then
      echo "  [ok]   $COH byte-identical"
    else
      echo "  [FAIL] $COH DIFFERS from the 50 bp run -- $_new vs $_old"
      _cfail=1
    fi
  done
  if [ "$_cfail" = "1" ]; then
    echo
    echo "  ⚠ The cohort arms depend on NOTHING that changed between read"
    echo "    lengths. If they differ, phase 2 is not reproducible and the"
    echo "    MAGNet arms from this run should not be used either."
  fi
fi

echo
echo "############### SUMMARY : JC vs JCEC ###############"
# <<'PY' is QUOTED, so bash does NOT substitute into the block.
# SUF therefore has to arrive by environment, not interpolation.
TBL_SUFFIX="$SUF" python3 - <<'PY'
import os
import numpy as np, pandas as pd
from scipy import stats

ET = "SE"
SUF = os.environ.get("TBL_SUFFIX", "")   # passed in, not substituted
PROGRAM = ["TTN","CAMK2D","LDB3","LMO7","SORBS1","PDLIM3","DTNA","MLIP","TNNT2"]

def deming(x, y, lam):
    if len(x) < 5: return np.nan
    sxx, syy = np.var(x, ddof=1), np.var(y, ddof=1)
    sxy = np.cov(x, y, ddof=1)[0, 1]
    if abs(sxy) < 1e-12: return np.nan
    a = syy - lam * sxx
    return (a + np.sqrt(a*a + 4*lam*sxy*sxy)) / (2*sxy)

def load(p):
    if not (os.path.exists(p) and os.path.getsize(p)): return None
    t = pd.read_csv(p, sep="\t").drop_duplicates("key")
    return t if {"key","gene","dPSI_adj","se_adj","FDR"} <= set(t.columns) else None

def slope(d, h, keys=None, prog=False):
    """Deming of HCM dPSI on DCM dPSI over DCM-significant events."""
    s = d[(d.FDR < 0.05) & (d.dPSI_adj.abs() >= 0.10)]
    m = s[["key","gene","dPSI_adj","se_adj"]].rename(
            columns={"dPSI_adj":"x","se_adj":"sx"}).merge(
        h[["key","dPSI_adj","se_adj"]], on="key").dropna(subset=["x","dPSI_adj"])
    if keys is not None:
        m = m[m.key.isin(keys)]
    if prog:
        m = m[m.gene.isin(PROGRAM)]
    if len(m) < 5:
        return np.nan, np.nan, len(m), np.nan
    x, y = m.x.values, m.dPSI_adj.values
    lam = np.nanmean(m.se_adj.values**2) / max(np.nanmean(m.sx.values**2), 1e-12)
    return (deming(x, y, lam), stats.pearsonr(x, y)[0], len(m),
            100*np.mean(np.sign(x) == np.sign(y)))

# (label, DCM JCEC, DCM JC, HCM JCEC, HCM JC)
# MAGNet paths take SUF; the three cohorts' own tables never do.
ARMS = [
 ("MAGNet disjoint",
  f"splitctl/{ET}_DCMvsA_1{SUF}.all_events.tsv",  f"jccheck{SUF}/{ET}_DCMvsA_JC.all_events.tsv",
  f"splitctl/{ET}_HCMvsB_1{SUF}.all_events.tsv",  f"jccheck{SUF}/{ET}_HCMvsB_JC.all_events.tsv"),
]
for C in ("Brigham", "Fuwai", "Stanford"):
    ARMS.append((C,
        f"{ET}_v2filtered{SUF}.all_events.tsv",   f"jccheck{SUF}/{ET}_v2filtered_JC.all_events.tsv",
        f"{C}_{ET}_HCM.all_events.tsv",           f"jccheck{SUF}/{C}_{ET}_HCM_JC.all_events.tsv"))

rows = []
for lab, dce, djc, hce, hjc in ARMS:
    D_ce, D_jc, H_ce, H_jc = load(dce), load(djc), load(hce), load(hjc)
    miss = [p for p, t in ((dce,D_ce),(djc,D_jc),(hce,H_ce),(hjc,H_jc)) if t is None]
    if miss:
        print(f"[skip] {lab}: missing/unreadable {miss}")
        continue
    common = (set(D_ce.key) & set(D_jc.key) & set(H_ce.key) & set(H_jc.key))
    for prog in (False, True):
        a = slope(D_ce, H_ce, None, prog)
        b = slope(D_jc, H_jc, None, prog)
        c = slope(D_ce, H_ce, common, prog)
        e = slope(D_jc, H_jc, common, prog)
        rows.append({"arm": lab, "set": "program" if prog else "background",
                     "JCEC_own": a[0], "n_JCEC": a[2], "JC_own": b[0], "n_JC": b[2],
                     "JCEC_matched": c[0], "JC_matched": e[0], "n_matched": c[2],
                     "pct_JCEC": c[3], "pct_JC": e[3]})
R = pd.DataFrame(rows)
if R.empty:
    raise SystemExit("[fatal] nothing to compare -- see the [skip] lines above")
_out = f"jc_vs_jcec{SUF}.tsv"
R.to_csv(_out, sep="\t", index=False, float_format="%.6g")
print(f"[write] {_out}")

for st in ("background", "program"):
    S = R[R.set == st]
    if S.empty: continue
    print(f"\n--- {st} ---")
    print(f"  {'arm':<17}{'JCEC own':>10}{'n':>7}{'JC own':>9}{'n':>7}"
          f"{'JCEC match':>12}{'JC match':>10}{'n':>7}{'%same JC':>10}")
    for _, r in S.iterrows():
        print(f"  {r.arm:<17}{r.JCEC_own:>10.3f}{r.n_JCEC:>7}{r.JC_own:>9.3f}{r.n_JC:>7}"
              f"{r.JCEC_matched:>12.3f}{r.JC_matched:>10.3f}{r.n_matched:>7}"
              f"{r.pct_JC:>9.1f}%")

print("""
HOW TO READ IT -- use the MATCHED columns, they are scored on the identical events.

  The external cohorts' JC-matched slopes stay near their JCEC values (0.44-0.70)
      -> read length is NOT the cause. All three technical explanations are now
         excluded and the disease-stage reading can be stated in the paper.

  The external cohorts' JC-matched slopes rise toward the MAGNet arm (~0.9)
      -> read length WAS the cause. Withdraw the stage-gradient claim and relabel
         panel B of fig_split_control.

  Watch the MAGNet disjoint row as the internal control: MAGNet is 50bp in both
  tables, so ITS slope should barely move between JC and JCEC. If it moves as much
  as the external cohorts do, the change is something about JC itself rather than
  read length, and neither column settles the question.
""")
PY
date
