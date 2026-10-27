#!/usr/bin/env python
"""
78_lvef_skew_origin.py -- where does the global LVEF sign skew come from?

THE OBSERVATION
---------------
77_lvef_s3.py found ZERO events associated with LVEF within DCM, and at the same time
a global tilt: 53.5% of the 38,084 covariate-adjusted LVEF coefficients are negative
(mean rho -0.0075, median -0.0070 -- the whole distribution shifted by a trivial
amount, not a tail of real events). Adding CM_score_frac + FB_score does NOT remove it;
the excess goes 3.5 -> 4.3 points. **Cell composition is not the explanation.**

THE QUESTION THIS SCRIPT ANSWERS, AND WHY IT MATTERS MORE THAN THE SKEW
-----------------------------------------------------------------------
A negative LVEF coefficient means PSI rises as ejection fraction falls -- i.e. as
disease worsens. So the skew has a direction, and that direction can be compared with
the DISEASE direction from the discovery contrast.

  * If events that shift in DCM-vs-NF shift FURTHER as LVEF falls, this is a
    severity DOSE-RESPONSE: a supportive finding, not a nuisance.
  * If the tilt is unrelated to the disease direction, it is a nuisance to disclose.

And one check matters regardless of the skew question:

  * (A) IS dPSI ITSELF GLOBALLY SKEWED in the discovery contrast? Nothing in the
    manuscript reports the aggregate sign balance of the 62,971 tested events. If
    that balance is off 50/50, it is worth knowing before submission, because it
    bears on the headline contrast and not just on a supplement.

PRIORITY. (A) is worth running whatever you decide about the skew. (B)-(E) are for the
skew question, which is currently claimed NOWHERE in the manuscript -- no section
reports aggregate sign balance -- so a reviewer will not raise it. Treat them as
optional, and do not let a 3.5-point tilt with a mean rho of -0.0075 become a
manuscript claim without a mechanism.

WHAT IT DOES
------------
  (A) sign balance of dPSI_adj: all 62,971 representatives, and the 38,084 analysed
  (B) 2x2 agreement between sign(dPSI_adj) and sign(-beta_rank), plus the same on the
      926 significant events, where the disease direction is actually confident
  (C) skew stratified by mean PSI -- a bounded variable has asymmetric room to move,
      so if the tilt lives at one end of the PSI range it is an estimator artefact
  (D) skew stratified by coverage (n) -- if it lives in the low-coverage events,
      library depth is implicated, and depth is NOT in the S3 model
  (E) skew stratified by |dPSI| -- if it tracks disease effect size, that is the
      dose-response reading again, quantified

GATES
  G1  the base table collapses to 62,971 exon-level events and carries 926 significant
  G2  the S3 file joins to the base table on ekey with no unmatched rows
  G3  the S3 file has the expected number of rows

USAGE
-----
  conda activate rnaseq_tools
  cd /data/project/chureeratp/RNA/HEART

  python 78_lvef_skew_origin.py                        # uses S3_lvef_compadj.tsv
  python 78_lvef_skew_origin.py --s3 S3_lvef.tsv       # the unadjusted run
  python 78_lvef_skew_origin.py --beta-col beta_rank_adj   # composition-adjusted coef

EVERY p-VALUE HERE IS DEPENDENCE-INFLATED. 38k exons share 160 samples and their
coverage, so a binomial p over events is not a valid test and is printed only as an
order of magnitude. Read the SIZE of each excess and whether it moves across strata.
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

EXPECT_TESTED = 62_971
EXPECT_SIG = 926
EXPECT_S3_ROWS = 38_084

NKEEP_SE = 3
C_KEY = "key"
FDR_MAX, DPSI_MIN = 0.05, 0.10
C_DPSI_CANDIDATES = ("dPSI_adj", "dPSI_raw", "dPSI")


def die(msg):
    print("\n[REFUSE] " + msg + "\n", file=sys.stderr)
    sys.exit(2)


def ok(msg):
    print("  [ok] " + msg)


def balance(v, label, indent="  "):
    """Sign balance of a vector, with an explicitly discounted p-value."""
    v = pd.Series(v).dropna()
    v = v[v != 0]
    if len(v) < 10:
        print(f"{indent}{label:<34}: too few events ({len(v)})")
        return np.nan
    neg = int((v < 0).sum())
    f = neg / len(v)
    p = stats.binomtest(neg, len(v), 0.5).pvalue
    print(f"{indent}{label:<34}: {neg:,}/{len(v):,} negative = {f:6.1%}"
          f"   excess {abs(f - 0.5) * 100:4.1f} pts   (p~{p:.1g})")
    return f


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--s3", default="S3_lvef_compadj.tsv")
    ap.add_argument("--events", default="SE_v2filtered.all_events.tsv")
    ap.add_argument("--beta-col", default="beta_rank",
                    help="beta_rank (base) or beta_rank_adj (composition-adjusted)")
    ap.add_argument("--strata", type=int, default=4)
    a = ap.parse_args()

    for p in (a.s3, a.events):
        if not os.path.exists(p):
            die(f"{p} not found")

    # ---------------------------------------------------------------- base table
    print("=" * 74)
    print("GATES")
    print("=" * 74)
    ev = pd.read_csv(a.events, sep="\t")
    dcol = next((c for c in C_DPSI_CANDIDATES if c in ev.columns), None)
    if dcol is None:
        die(f"{a.events} has none of {C_DPSI_CANDIDATES}")
    ev["dPSI"] = pd.to_numeric(ev[dcol], errors="coerce")
    ev["ekey"] = ev[C_KEY].astype(str).str.split(":").str[:NKEEP_SE].str.join(":")
    rep = (ev.sort_values("FDR", kind="mergesort")
             .drop_duplicates(subset=["ekey"], keep="first"))
    if len(rep) != EXPECT_TESTED:
        die(f"G1 FAILED: {len(rep):,} exon-level events, expected {EXPECT_TESTED:,}")
    nsig = int(((rep.FDR < FDR_MAX) & (rep.dPSI.abs() >= DPSI_MIN)).sum())
    if nsig != EXPECT_SIG:
        die(f"G1 FAILED: {nsig} significant, expected {EXPECT_SIG}")
    ok(f"G1 base table {len(rep):,} tested / {nsig} significant  (effect '{dcol}')")

    # ---------------------------------------------------------------- S3 table
    s3 = pd.read_csv(a.s3, sep="\t")
    if a.beta_col not in s3.columns:
        die(f"{a.s3} has no '{a.beta_col}'. Columns: {list(s3.columns)}")
    if len(s3) != EXPECT_S3_ROWS:
        print(f"  [warn] {a.s3} has {len(s3):,} rows, expected {EXPECT_S3_ROWS:,}")
    else:
        ok(f"G3 S3 rows = {len(s3):,}")

    j = s3.merge(rep[["ekey", "dPSI", "FDR", "mean_psi"]], on="ekey", how="left",
                 suffixes=("", "_base"))
    unmatched = int(j["dPSI"].isna().sum())
    if unmatched:
        die(f"G2 FAILED: {unmatched:,} S3 rows did not join to the base table")
    ok(f"G2 joined {len(j):,} S3 rows to the base table on ekey")

    print(f"\n  coefficient column: {a.beta_col}")
    print("  a NEGATIVE coefficient means PSI RISES as LVEF FALLS, i.e. with severity")

    # ---------------------------------------------------------------- (A)
    print("\n" + "=" * 74)
    print("(A) IS dPSI ITSELF GLOBALLY SKEWED?  <- worth knowing regardless")
    print("=" * 74)
    print("  Sign balance of the DISCOVERY contrast's effect. Nothing in the")
    print("  manuscript reports this, and it bears on the headline, not a supplement.")
    print("  (negative dPSI = exon LESS included in DCM than non-failing)\n")
    balance(rep["dPSI"], f"all tested ({len(rep):,})")
    balance(j["dPSI"], f"the {len(j):,} analysed in S3")
    sig = rep[(rep.FDR < FDR_MAX) & (rep.dPSI.abs() >= DPSI_MIN)]
    balance(sig["dPSI"], f"the {len(sig):,} significant")
    print("\n  If the first row is near 50% the discovery contrast is balanced and the")
    print("  LVEF tilt is NOT a shadow of it. If it is off 50% by a similar margin in")
    print("  the SAME direction, the two are likely one phenomenon -- and the 926 row")
    print("  says whether real events share it or only the undifferentiated bulk.")

    # ---------------------------------------------------------------- (B)
    print("\n" + "=" * 74)
    print("(B) DOES THE LVEF TILT FOLLOW THE DISEASE DIRECTION?")
    print("=" * 74)
    print("  Dose-response predicts AGREEMENT: an exon more included in DCM becomes")
    print("  more included still as LVEF falls -> dPSI > 0 with coefficient < 0.\n")

    b = j[a.beta_col]
    agree = ((j["dPSI"] > 0) & (b < 0)) | ((j["dPSI"] < 0) & (b > 0))
    usable = j["dPSI"].notna() & b.notna() & (j["dPSI"] != 0) & (b != 0)
    fa = float(agree[usable].mean())
    n_us = int(usable.sum())
    pa = stats.binomtest(int(agree[usable].sum()), n_us, 0.5).pvalue
    print(f"  all {n_us:,} analysed: agreement = {fa:.1%}   "
          f"excess {abs(fa - 0.5) * 100:.1f} pts   (p~{pa:.1g})")

    for lbl, sub in ((f"the {nsig} significant",
                      usable & (j.FDR < FDR_MAX) & (j["dPSI"].abs() >= DPSI_MIN)),
                     ("FDR<0.05, any effect size", usable & (j.FDR < FDR_MAX))):
        if int(sub.sum()) >= 10:
            f = float(agree[sub].mean())
            pp = stats.binomtest(int(agree[sub].sum()), int(sub.sum()), 0.5).pvalue
            print(f"  {lbl:<26} ({int(sub.sum()):,}): agreement = {f:.1%}   "
                  f"excess {abs(f - 0.5) * 100:.1f} pts   (p~{pp:.1g})")

    print("\n  READ IT LIKE THIS:")
    print("    agreement near 50%  -> the tilt is independent of disease direction.")
    print("                           A nuisance. Disclose it or ignore it; do not")
    print("                           dress it as biology.")
    print("    agreement well >50% AND holding on the 926 -> dose-response, and the")
    print("                           926 row is the one that counts, because there")
    print("                           the disease direction is actually established.")
    print("    agreement >50% on the bulk but NOT on the 926 -> a property of noise,")
    print("                           not of the signature. The likeliest trap here.")

    # ---------------------------------------------------------------- (C)(D)(E)
    print("\n" + "=" * 74)
    print("(C)(D)(E) WHERE DOES THE TILT LIVE?")
    print("=" * 74)
    print("  A tilt that is flat across strata is global. One concentrated in a")
    print("  stratum names its own cause.\n")

    for col, label, note in (
        ("mean_psi", "(C) mean PSI", "boundedness: PSI near 0 or 1 has asymmetric room"),
        ("n", "(D) coverage n", "library depth -- NOT a covariate in the S3 model"),
    ):
        if col not in j.columns:
            print(f"  {label}: column '{col}' absent, skipped")
            continue
        print(f"  {label} -- {note}")
        try:
            q = pd.qcut(j[col], a.strata, duplicates="drop")
            n_lev = int(q.cat.categories.size)
        except (ValueError, AttributeError):
            n_lev = 0
        if n_lev < 2:
            print(f"    '{col}' has too few distinct values to stratify "
                  f"({j[col].nunique()} unique) -- skipped\n")
            continue
        for lev, grp in j.groupby(q, observed=True):
            balance(grp[a.beta_col], f"  {lev}", indent="    ")
        print()

    print("  (E) |dPSI| -- does the tilt track disease effect size?")
    try:
        q = pd.qcut(j["dPSI"].abs(), a.strata, duplicates="drop")
        for lev, grp in j.groupby(q, observed=True):
            sub = grp[grp["dPSI"] != 0]
            ag = ((sub["dPSI"] > 0) & (sub[a.beta_col] < 0)) | \
                 ((sub["dPSI"] < 0) & (sub[a.beta_col] > 0))
            if len(sub) >= 10:
                print(f"    {str(lev):<34}: agreement = {ag.mean():6.1%}"
                      f"   (n = {len(sub):,})")
        print("\n    Rising agreement with |dPSI| is the dose-response signature.")
        print("    Flat agreement means the tilt has nothing to do with disease.")
    except ValueError:
        print("    not enough distinct values to stratify")

    # ---------------------------------------------------------------- verdict
    print("\n" + "=" * 74)
    print("WHAT TO DO WITH THIS")
    print("=" * 74)
    print("  Whatever the numbers, the per-event null stands: 0 of 38,084 at")
    print("  FDR < 0.05, PLA2G5 nearest at 0.0501 with an effect of only -0.083,")
    print("  and the one event reaching 0.1009 nowhere near significance.")
    print()
    print("  This tilt is reportable ONLY with a mechanism. mean rho = -0.0075 is")
    print("  not a finding on its own, and the binomial p-values above are inflated")
    print("  by dependence and must not be quoted. If (B) shows dose-response ON THE")
    print("  926, that is a real and supportive observation worth one sentence. If")
    print("  (C) or (D) concentrates the tilt, it is an artefact of boundedness or")
    print("  depth and belongs nowhere. If everything is flat and near 50%, the")
    print("  honest conclusion is that a 3.5-point tilt of unknown origin exists,")
    print("  is too small to matter, and goes unmentioned.")
    print()
    print("  Composition is already EXCLUDED as the explanation: adding")
    print("  CM_score_frac + FB_score moved the excess 3.5 -> 4.3 points.")
    print()


if __name__ == "__main__":
    main()
