#!/bin/bash
#SBATCH -p cpu
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH --qos=ngs16
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#
# 32 -- rerun every contrast affected by the genotype-duplicate exclusions.
#
# WHY NO rMATS RERUN
# ------------------
# PSI is per-sample, so phase2_diffsplice_v2.py --exclude drops samples from the
# EXISTING rMATS output. This is the same mechanism 04_phase2_v2_checks.sbatch already
# used for the age sensitivity analysis (--exclude young_under25.ids against the
# unmodified magnet_covariates.tsv), so the .cov files are left alone here too: an
# extra row in a .cov file for a sample that is no longer in the PSI matrix is inert.
#
# WHAT IS RERUN, AND WHY EACH
#   Colorado DCM vs NF     36 vs 13   SRR7426833 (DCM) and SRR7426789 (NF) were
#   Colorado ICM vs NF     13 vs 13   technical replicates of SRR7426832 and
#   Colorado DCM vs ICM    36 vs 13   SRR7426790 -- they inflated effective n
#   Stanford HCM vs ctrl   13 vs  6   SRR15183558 replicated SRR15183553
#   LBNL DCM vs NF          5 vs  4   OPTIONAL, see below
#
# MAGNET DISCOVERY IS NOT RERUN. No MAGNet sample was excluded -- the keep rule always
# retains the discovery copy, and no MAGNet-MAGNet duplicate exists. Every MAGNet number,
# the event convention and the composition/GO results stand unchanged.
#
# LBNL IS NOT A REPLICATION ANALYSIS ANY MORE. 24 of its 33 samples are MAGNet hearts.
# At 5 vs 4 it has no power and must never be quoted as independent replication. It is
# rerun only so the library-prep robustness panel has a dedup'd table; set RUN_LBNL=0
# to skip it. --min-samples drops to 4 because the original 8 exceeds the group sizes.
#
# OUTPUTS ARE NEW FILES (_dedup suffix). Nothing existing is overwritten, so the
# before/after comparison at the end is meaningful and you can revert by deleting them.
#
# Run:  bash 32_rerun_dedup_contrasts.sbatch 2>&1 | tail -60
#   or: sbatch 32_rerun_dedup_contrasts.sbatch
date
export PATH=/data/apps/anaconda3/envs/mamba/envs/rnaseq_tools/bin:$PATH
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
ROOT=/data/project/chureeratp/RNA/HEART
cd "$ROOT"

EXCL=${EXCL:-exclude_dups.ids}
RUN_LBNL=${RUN_LBNL:-1}
TYPES="SE MXE A3SS A5SS RI"

[ -s "$EXCL" ] || { echo "*** $EXCL missing. Run:"
                    echo "***   python 30_build_dedup_lists.py --apply-to-bams"
                    echo "***   awk 'NR>1 {print \$1}' geno/EXCLUDE_samples.tsv | sort -u > $EXCL"
                    exit 1; }
NE=$(awk 'NF' "$EXCL" | wc -l)
echo "exclusion list: $EXCL ($NE samples)"
[ "$NE" -ge 20 ] || echo "[warn] only $NE ids -- expected 27; is the list current?"

# Stanford's --min-samples lives in an associative array in 15b. Read it rather than
# guess: a wrong value silently changes how many events are testable, which would
# contaminate the before/after comparison this script exists to produce.
SMIN=$(grep -oP 'MINS\[Stanford\]=\K[0-9]+' 15b_hcm_cohorts_rmats.sbatch 2>/dev/null | head -1)
if [ -z "$SMIN" ]; then
  SMIN=5
  echo "[warn] could not read MINS[Stanford] from 15b_hcm_cohorts_rmats.sbatch;"
  echo "       defaulting to --min-samples $SMIN. Check this before quoting the result."
else
  echo "Stanford --min-samples read from 15b: $SMIN"
fi
if [ "$SMIN" -gt 6 ]; then
  echo "*** MINS[Stanford]=$SMIN but the control arm is now n=6. Every event would be"
  echo "*** dropped. Lower it deliberately and say so in the methods, or skip Stanford."
  exit 1
fi

run1 () {   # run1 <mats> <et> <b1> <b2> <cov> <g1> <g2> <minsamp> <outprefix>
  local M="$1" ET="$2"
  if [ ! -s "$M" ]; then echo "  missing $M -- skipping"; return; fi
  for f in "$3" "$4" "$5"; do
    [ -s "$f" ] || { echo "  missing $f -- skipping $9"; return; }
  done
  echo "  -> $9"
  python phase2_diffsplice_v2.py --mats "$M" --event-type "$ET" \
    --b1 "$3" --b2 "$4" --cov "$5" --group1 "$6" --group2 "$7" \
    --min-samples "$8" --no-covariates --exclude "$EXCL" --out "$9" \
    2>&1 | grep -Ei "error|traceback|events (kept|tested)|significant" | head -4
}

for ET in $TYPES; do
  echo
  echo "############################################################"
  echo "#####  $ET"
  echo "############################################################"
  run1 rmats_Colorado/${ET}.MATS.JCEC.txt          "$ET" \
       Colorado_dcm_bams.txt Colorado_nf_bams.txt  colorado_dcm_nf.cov \
       DCM NF 10 ${ET}_Colorado_DCMvsNF_dedup
  run1 rmats_Colorado_ICMvsNF/${ET}.MATS.JCEC.txt  "$ET" \
       Colorado_icm_bams.txt Colorado_nf_bams.txt  colorado_icm_nf.cov \
       ICM NF 10 ${ET}_Colorado_ICMvsNF_dedup
  run1 rmats_Colorado_DCMvsICM/${ET}.MATS.JCEC.txt "$ET" \
       Colorado_dcm_bams.txt Colorado_icm_bams.txt colorado_dcm_icm.cov \
       DCM ICM 10 ${ET}_Colorado_DCMvsICM_dedup
  run1 rmats_Stanford/${ET}.MATS.JCEC.txt          "$ET" \
       Stanford_hcm_bams.txt Stanford_ctrl_bams.txt Stanford.cov \
       HCM CTRL "$SMIN" Stanford_${ET}_HCM_dedup
  if [ "$RUN_LBNL" = "1" ]; then
    run1 rmats_LBNL/${ET}.MATS.JCEC.txt            "$ET" \
         LBNL_dcm_bams.txt LBNL_nf_bams.txt        lbnl_dcm_nf.cov \
         DCM NF 4 ${ET}_LBNL_DCMvsNF_dedup
  fi
done

echo
echo "############################################################"
echo "#####  BEFORE vs AFTER"
echo "############################################################"
python3 - <<'PY'
import glob
import os

import numpy as np
import pandas as pd

TYPES = ["SE", "MXE", "A3SS", "A5SS", "RI"]
PAIRS = [("Colorado DCMvsNF",  "{ET}_Colorado_DCMvsNF"),
         ("Colorado ICMvsNF",  "{ET}_Colorado_ICMvsNF"),
         ("Colorado DCMvsICM", "{ET}_Colorado_DCMvsICM"),
         ("Stanford HCMvsCTRL", "Stanford_{ET}_HCM"),
         ("LBNL DCMvsNF",      "{ET}_LBNL_DCMvsNF")]


def stats(path):
    if not (os.path.exists(path) and os.path.getsize(path)):
        return None
    t = pd.read_csv(path, sep="\t")
    fdr = next((c for c in ("FDR", "fdr", "padj") if c in t.columns), None)
    dps = next((c for c in ("dPSI_adj", "dPSI", "IncLevelDifference")
                if c in t.columns), None)
    if fdr is None:
        return None
    sig = t[fdr] < 0.05
    out = {"tested": len(t), "sig": int(sig.sum())}
    if dps is not None:
        out["sig05"] = int((sig & (t[dps].abs() >= 0.05)).sum())
        out["med_abs_dpsi"] = float(t.loc[sig, dps].abs().median()) if sig.any() else np.nan
    return out


rows = []
print(f"  {'contrast':<20}{'type':<6}{'tested':>18}{'FDR<0.05':>16}"
      f"{'+|dPSI|>=.05':>16}{'median |dPSI|':>22}")
for label, pat in PAIRS:
    seen = False
    for ET in TYPES:
        base = pat.format(ET=ET)
        a = stats(f"{base}.all_events.tsv")
        b = stats(f"{base}_dedup.all_events.tsv")
        if a is None and b is None:
            continue
        seen = True
        f = lambda d, k: ("--" if d is None or k not in d
                          or (isinstance(d.get(k), float) and np.isnan(d[k]))
                          else (f"{d[k]:.4f}" if isinstance(d[k], float) else f"{d[k]}"))
        cell = lambda k: f"{f(a, k)} -> {f(b, k)}"
        print(f"  {label:<20}{ET:<6}{cell('tested'):>18}{cell('sig'):>16}"
              f"{cell('sig05'):>16}{cell('med_abs_dpsi'):>22}")
        if a and b:
            rows.append({"contrast": label, "event_type": ET,
                         "tested_before": a["tested"], "tested_after": b["tested"],
                         "sig_before": a["sig"], "sig_after": b["sig"],
                         "sig_pct_change": (100 * (b["sig"] - a["sig"]) / a["sig"]
                                            if a["sig"] else np.nan)})
    if not seen:
        print(f"  {label:<20}{'--':<6}{'no tables found':>16}")

if rows:
    R = pd.DataFrame(rows)
    R.to_csv("dedup_before_after.tsv", sep="\t", index=False, float_format="%.4g")
    print(f"\n  -> dedup_before_after.tsv")
    print()
    print("  HOW TO READ THIS")
    print("  Colorado and Stanford each lost ONE sample per affected arm, so the")
    print("  significant counts should move by a few percent at most. A large move")
    print("  means the result was resting on that one sample and the original number")
    print("  was not robust -- that is worth knowing, not worth hiding.")
    # A percentage on a count of 8 is not informative: losing two events is a 25%
    # "move" that means nothing. Only judge contrasts with enough events for the
    # percentage to carry information, and list the small ones separately as raw
    # counts so they are still visible.
    MINSIG = 20
    big = R[(R.sig_pct_change.abs() > 15) & (R.sig_before >= MINSIG)]
    small = R[R.sig_before < MINSIG]
    print()
    if len(big):
        print("  *** MOVED BY MORE THAN 15% -- inspect before using either number:")
        for _, r in big.iterrows():
            print(f"      {r.contrast} {r.event_type}: {int(r.sig_before)} -> "
                  f"{int(r.sig_after)}  ({r.sig_pct_change:+.1f}%)")
        print("      Removing one sample should not do this. Check whether the")
        print("      excluded sample was an outlier in that arm before quoting either.")
    else:
        print(f"  No contrast with >={MINSIG} significant events moved by more than 15%.")
        print("  The exclusions did not change any conclusion; quote the _dedup numbers,")
        print("  state that duplicates were removed, and keep the originals as a")
        print("  sensitivity check.")
    if len(small):
        print()
        print(f"  Contrasts with <{MINSIG} significant events -- percentages are not")
        print("  meaningful at this scale, so read the raw counts:")
        for _, r in small.iterrows():
            print(f"      {r.contrast} {r.event_type}: {int(r.sig_before)} -> "
                  f"{int(r.sig_after)}")
    print()
    print("  LBNL: do NOT read its row as replication. At 5 vs 4 it is a library-prep")
    print("  panel only. 24 of its 33 samples are MAGNet hearts.")
PY
echo
echo "Next: rebuild the downstream summaries that read these tables"
echo "  -- 12_colorado_icm_analysis.sbatch concordance section (Colorado)"
echo "  -- 15b_hcm_cohorts_rmats.sbatch specificity section (Stanford)"
echo "  -- 24_lock_event_convention.py to refresh master_numbers.tsv"
echo "  -- 26_flag_master_numbers.py to re-apply the reliability verdicts"
date
