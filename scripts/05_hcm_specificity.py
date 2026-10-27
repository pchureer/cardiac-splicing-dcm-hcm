#!/usr/bin/env python3
"""
05_hcm_specificity.py -- is the RBM20 splicing program DCM-specific or shared with HCM?

Takes two phase2_diffsplice_v2 output tables (DCM-vs-NF and HCM-vs-NF), merges them
on the full event coordinate `key`, and answers the question the way it has to be
answered when the two cohorts have wildly different power (166 vs 28).

THREE THINGS THIS DOES THAT A NAIVE COMPARISON GETS WRONG
---------------------------------------------------------
1. It does not compare significance. "Significant in DCM, not significant in HCM"
   is not evidence of specificity when HCM has a sixth of the sample size -- absence
   of significance is mostly absence of power. Every event is instead labelled using
   its HCM *minimum detectable effect* (mde80): an event is only called DCM-specific
   if HCM could have detected an effect the size of the DCM one and did not.

2. It corrects for SHARED CONTROLS. Your HCM run reuses nf_bams.txt, the same 166 NF
   samples as discovery. Both dPSI estimates therefore subtract the same NF mean, so
   they are positively correlated even when HCM has no real effect at all. With
   n_DCM=166, n_HCM=28, n_NF=166 the expected null correlation is about +0.27. The
   script computes this baseline from the actual per-event group sizes and reports
   the observed correlation against it. Ignoring this is how a null result gets
   published as "a shared program".

3. It reports the regression slope of HCM effect on DCM effect (how much of the DCM
   program HCM recapitulates), with the errors-in-variables (Deming) slope alongside
   the OLS one, because OLS is biased toward zero when the x-axis is itself estimated.

Usage:
  python 05_hcm_specificity.py \
     --dcm SE_v2filtered.all_events.tsv \
     --hcm SE_HCM_v2.all_events.tsv \
     --genes TTN,CAMK2D,LDB3,LMO7,SORBS1,PDLIM3,DTNA,MLIP,TNNT2 \
     --well-powered CAMK2D,LDB3,LMO7,SORBS1 \
     --out SE_specificity
"""
import argparse
import sys

import numpy as np
import pandas as pd
from scipy import stats


def deming_slope(x, y, lam):
    """Errors-in-variables slope; lam = var(err_y)/var(err_x)."""
    sxx = np.var(x, ddof=1)
    syy = np.var(y, ddof=1)
    sxy = np.cov(x, y, ddof=1)[0, 1]
    if abs(sxy) < 1e-12:
        return np.nan
    a = syy - lam * sxx
    return (a + np.sqrt(a * a + 4 * lam * sxy * sxy)) / (2 * sxy)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dcm", required=True, help="phase2_v2 table, DCM vs NF")
    ap.add_argument("--hcm", required=True, help="phase2_v2 table, HCM vs NF")
    ap.add_argument("--genes", default="TTN,CAMK2D,LDB3,LMO7,SORBS1,PDLIM3,DTNA,MLIP,TNNT2")
    ap.add_argument("--well-powered", default="CAMK2D,LDB3,LMO7,SORBS1",
                    help="genes to lead the specificity claim with (titin is noisy per-exon)")
    ap.add_argument("--fdr", type=float, default=0.05)
    ap.add_argument("--min-dpsi", type=float, default=0.10)
    ap.add_argument("--hcm-alpha", type=float, default=0.05)
    ap.add_argument("--shared-controls", action="store_true", default=True,
                    help="HCM run reuses nf_bams.txt (default true; --no-shared-controls "
                         "if you split the NF group)")
    ap.add_argument("--no-shared-controls", dest="shared_controls", action="store_false")
    ap.add_argument("--out", default="specificity")
    a = ap.parse_args()

    d = pd.read_csv(a.dcm, sep="\t")
    h = pd.read_csv(a.hcm, sep="\t")
    for name, t in (("--dcm", d), ("--hcm", h)):
        need = {"key", "gene", "dPSI_adj", "pvalue", "FDR", "n_g1", "n_g2"}
        miss = need - set(t.columns)
        if miss:
            sys.exit(f"[fatal] {name} is missing {sorted(miss)}. Both tables must come "
                     "from phase2_diffsplice_v2.py (v1 has no `key` column).")
    if "mde80" not in h.columns:
        if "se_adj" not in h.columns:
            sys.exit("[fatal] --hcm table has neither mde80 nor se_adj. Re-run "
                     "phase2_diffsplice_v2.py (current version emits both); without a "
                     "standard error there is no way to tell 'no effect in HCM' from "
                     "'no power in HCM', and the specificity claim is unsupportable.")
        h["mde80"] = 2.802 * h["se_adj"]

    lines = []
    def say(s=""):
        lines.append(s)
        print(s)

    genes = [g.strip() for g in a.genes.split(",") if g.strip()]
    wp = [g.strip() for g in a.well_powered.split(",") if g.strip()]

    # discovery set: DCM-significant, real-effect events in the RBM20 program
    dsig = d[(d.FDR < a.fdr) & (d.dPSI_adj.abs() >= a.min_dpsi)].copy()
    dprog = dsig[dsig.gene.isin(genes)].copy()
    say("=" * 68)
    say(f"DCM-vs-NF significant events (FDR<{a.fdr}, |dPSI|>={a.min_dpsi}): {len(dsig)}")
    say(f"  of those in the RBM20 program genes: {len(dprog)}")

    m = dprog.merge(h, on="key", suffixes=("_DCM", "_HCM"), how="left")
    measured = m[m.dPSI_adj_HCM.notna()].copy()
    say(f"  measured in the HCM run as well: {len(measured)}"
        f"   (not measured: {len(m) - len(measured)})")
    if len(measured) < 5:
        say("\n[stop] too few events measured in both runs to say anything. Check that "
            "both rMATS runs used the same GTF and that the HCM run produced output.")
        open(f"{a.out}.report.txt", "w").write("\n".join(lines) + "\n")
        return
    say(f"  HCM group size (median across events): n={int(measured.n_g1_HCM.median())}")
    say()

    # ---- the shared-control null baseline -------------------------------
    if a.shared_controls:
        n_d = measured.n_g1_DCM.values.astype(float)
        n_h = measured.n_g1_HCM.values.astype(float)
        n_nf = measured.n_g2_DCM.values.astype(float)
        rho0 = ((1 / n_nf) /
                np.sqrt((1 / n_d + 1 / n_nf) * (1 / n_h + 1 / n_nf)))
        rho_null = float(np.median(rho0))
        say("--- 0. Shared-control baseline (READ THIS FIRST) ---")
        say(f"  The HCM run reuses the same {int(np.median(n_nf))} NF samples as discovery,")
        say("  so both dPSI estimates subtract the same control mean and are correlated")
        say("  even under a true null.")
        say(f"  EXPECTED null correlation from shared controls alone: rho = {rho_null:+.3f}")
        say("  Any observed correlation must be judged against that number, not against 0.")
        say()
    else:
        rho_null = 0.0
        say("--- 0. Controls are independent between the two runs (rho_null = 0) ---")
        say()

    x = measured.dPSI_adj_DCM.values.astype(float)
    y = measured.dPSI_adj_HCM.values.astype(float)

    # ---- 1. effect-size agreement ---------------------------------------
    rho, prho = stats.spearmanr(x, y)
    r, pr = stats.pearsonr(x, y)
    same = int(np.sum(np.sign(x) == np.sign(y)))
    pb = stats.binomtest(same, len(x), 0.5, alternative="greater").pvalue
    say("--- 1. Do the DCM and HCM effect sizes agree? ---")
    say(f"  Spearman rho = {rho:+.3f} (p={prho:.3g})   Pearson r = {r:+.3f} (p={pr:.3g})")
    say(f"  same direction: {same}/{len(x)} ({100*same/len(x):.0f}%)  binomial p = {pb:.3g}")
    say(f"  vs shared-control null of rho = {rho_null:+.3f}")
    if r <= rho_null + 0.05:
        say("  => at or below the null baseline: NO evidence of a shared program.")
        say("     This is the DCM-SPECIFIC result.")
    elif r < rho_null + 0.25:
        say("  => only marginally above the null baseline. Weak/no shared signal;")
        say("     do not describe the program as shared on this evidence.")
    else:
        say("  => clearly above the null baseline: the two diseases share part of the")
        say("     program. Report the slope below as 'how much', not 'whether'.")
    say()

    # ---- 2. how much of the DCM effect does HCM show? --------------------
    slope, icpt, rv, pv, sterr = stats.linregress(x, y)
    ci = 1.96 * sterr
    lam = 1.0
    if "se_adj_HCM" in measured.columns and "se_adj_DCM" in measured.columns:
        vh = np.nanmean(measured.se_adj_HCM.values ** 2)
        vd = np.nanmean(measured.se_adj_DCM.values ** 2)
        lam = float(vh / vd) if vd > 0 else 1.0
    dem = deming_slope(x, y, lam)
    say("--- 2. How much of the DCM effect does HCM recapitulate? ---")
    say(f"  OLS slope    = {slope:+.3f}  (95% CI {slope-ci:+.3f} to {slope+ci:+.3f})")
    say(f"  Deming slope = {dem:+.3f}  (errors-in-variables, lambda={lam:.2f})")
    say("  OLS is biased toward 0 here because the x-axis is itself an estimate;")
    say("  quote the Deming slope, or quote both.")
    say("  slope ~1 = same program same size; ~0.3 = attenuated; ~0 = DCM-specific.")
    say()

    # ---- 3. per-event classification with power awareness ---------------
    def classify(row):
        dd, hh = row.dPSI_adj_DCM, row.dPSI_adj_HCM
        if pd.isna(hh):
            return "not_measured_in_HCM"
        if row.pvalue_HCM < a.hcm_alpha and np.sign(dd) == np.sign(hh):
            return "shared"
        if row.pvalue_HCM < a.hcm_alpha and np.sign(dd) != np.sign(hh):
            return "opposite"
        mde = row.get("mde80_HCM", row.get("mde80", np.nan))
        if pd.isna(mde) or mde > abs(dd):
            return "underpowered_in_HCM"
        if abs(hh) < 0.5 * abs(dd):
            return "DCM_specific"
        return "inconclusive"

    measured["class"] = measured.apply(classify, axis=1)
    counts = measured["class"].value_counts()
    say("--- 3. Per-event classification (power-aware) ---")
    for k in ["shared", "DCM_specific", "underpowered_in_HCM", "inconclusive", "opposite"]:
        say(f"  {k:<22} {int(counts.get(k, 0))}")
    npow = int(counts.get("underpowered_in_HCM", 0))
    say(f"\n  {npow}/{len(measured)} events could NOT have been detected in HCM even if")
    say("  the effect were as large as in DCM. Those are silent on specificity --")
    say("  never count them as evidence for it.")
    say()

    # ---- 4. per-gene ----------------------------------------------------
    say("--- 4. Per-gene ---")
    say(f"  {'gene':<10}{'n_ev':>5}{'med dPSI_DCM':>14}{'med dPSI_HCM':>14}"
        f"{'shared':>8}{'DCMspec':>9}{'underpwr':>10}")
    for g, gg in measured.groupby("gene_DCM"):
        cc = gg["class"].value_counts()
        say(f"  {g:<10}{len(gg):>5}{np.median(gg.dPSI_adj_DCM):>+14.3f}"
            f"{np.median(gg.dPSI_adj_HCM):>+14.3f}"
            f"{int(cc.get('shared',0)):>8}{int(cc.get('DCM_specific',0)):>9}"
            f"{int(cc.get('underpowered_in_HCM',0)):>10}")
    say()

    # ---- 5. the claim you can actually make -----------------------------
    sub = measured[measured.gene_DCM.isin(wp)]
    say(f"--- 5. Lead claim, restricted to well-powered genes {wp} ---")
    if len(sub) >= 5:
        rr, pp = stats.pearsonr(sub.dPSI_adj_DCM.values, sub.dPSI_adj_HCM.values)
        sl = stats.linregress(sub.dPSI_adj_DCM.values, sub.dPSI_adj_HCM.values).slope
        cc = sub["class"].value_counts()
        say(f"  events: {len(sub)}   Pearson r = {rr:+.3f} (p={pp:.3g}, null {rho_null:+.3f})"
            f"   slope = {sl:+.3f}")
        say(f"  shared {int(cc.get('shared',0))} / DCM-specific {int(cc.get('DCM_specific',0))}"
            f" / underpowered {int(cc.get('underpowered_in_HCM',0))}")
    else:
        say(f"  only {len(sub)} events -- too few. Widen --well-powered or report titin"
            "  with explicit caveats.")
    say()
    say("  Titin is deliberately excluded from the lead claim: per-exon titin PSI is")
    say("  noisy and n=28 cannot support exon-level specificity statements.")
    say()

    # ---- 6. THE CONTEXT NUMBER: is the program special at all? ----------
    say("--- 6. Is the RBM20 program different from the genome-wide background? ---")
    allm = dsig.merge(h, on="key", suffixes=("_DCM", "_HCM"), how="inner")
    allm = allm[allm.dPSI_adj_HCM.notna()]
    if len(allm) >= 20:
        xa = allm.dPSI_adj_DCM.values.astype(float)
        ya = allm.dPSI_adj_HCM.values.astype(float)
        ra, pa = stats.pearsonr(xa, ya)
        sa = stats.linregress(xa, ya)
        lam_a = 1.0
        if "se_adj_HCM" in allm.columns and "se_adj_DCM" in allm.columns:
            va = np.nanmean(allm.se_adj_HCM.values ** 2)
            vb = np.nanmean(allm.se_adj_DCM.values ** 2)
            lam_a = float(va / vb) if vb > 0 else 1.0
        da = deming_slope(xa, ya, lam_a)
        samea = int(np.sum(np.sign(xa) == np.sign(ya)))
        say(f"  ALL DCM-significant events measured in HCM: n={len(allm)}")
        say(f"    Pearson r = {ra:+.3f}   OLS slope = {sa.slope:+.3f}   "
            f"Deming = {da:+.3f}   same direction {samea}/{len(xa)} "
            f"({100*samea/len(xa):.0f}%)")
        say(f"  RBM20 program only: n={len(measured)}")
        say(f"    Pearson r = {r:+.3f}   OLS slope = {slope:+.3f}   Deming = {dem:+.3f}"
            f"   same direction {same}/{len(x)} ({100*same/len(x):.0f}%)")
        say()
        say("  THIS IS THE COMPARISON THAT MATTERS FOR THE PAPER'S CLAIM.")
        if abs(da) > 0.5 and abs(ra) > 0.5:
            say("  The genome-wide background ALSO transfers to HCM at a similar slope.")
            say("  That means the RBM20 program is not behaving differently from")
            say("  diseased myocardium in general -- a shared end-stage signature,")
            say("  which is equally consistent with cell-composition differences")
            say("  (less myocardium, more fibrosis in any failing heart) as with a")
            say("  shared splicing program. You cannot distinguish those two without")
            say("  adjusting both models for cardiomyocyte fraction.")
        else:
            say("  The background does NOT transfer at this slope, so the RBM20")
            say("  program is behaving differently from the genome-wide average.")
            say("  That is a real, specific observation about the program.")
    else:
        say(f"  only {len(allm)} background events measured in both -- skipped")
    say("=" * 68)

    measured.to_csv(f"{a.out}.merged_events.tsv", sep="\t", index=False)
    with open(f"{a.out}.report.txt", "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"\nwrote {a.out}.merged_events.tsv and {a.out}.report.txt")


if __name__ == "__main__":
    main()
