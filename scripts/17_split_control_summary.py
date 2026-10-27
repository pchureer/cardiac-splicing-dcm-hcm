#!/usr/bin/env python3
"""
17_split_control_summary.py -- read the output of 16_split_control_check.sbatch and
answer the one question it was built for.

THE QUESTION
------------
Three independent HCM cohorts recapitulate the MAGNet DCM signature at a Deming slope
of 0.44-0.70. MAGNet's OWN DCM-vs-HCM comparison gives ~0.97. Why the gap?

  (a) MAGNet's 0.97 is INFLATED, because both of its contrasts subtract the same 166
      NF hearts, so every NF sampling fluctuation lands in both dPSI vectors; or
  (b) the external slope is genuinely LOWER, because myectomy / preserved-EF tissue
      carries a milder version of the signature than end-stage explants.

Script 16 built two arms that are identical in cases, tissue, read length, annotation
and control COUNT, and differ ONLY in whether the two contrasts share control samples:

    SHARED    DCM vs NF(A)   against   HCM vs NF(A)
    DISJOINT  DCM vs NF(A)   against   HCM vs NF(B)

So  slope(SHARED) - slope(DISJOINT)  IS the shared-control inflation, measured.

WHY THIS SCRIPT RECOMPUTES INSTEAD OF GREPPING THE REPORTS
----------------------------------------------------------
Two reasons. First, the two arms must be evaluated on the IDENTICAL event set --
otherwise part of the difference is just a different set of events, and script 05
selects each arm's set independently. This script intersects all three tables first.
Second, that shared event set makes the comparison PAIRED, which lets the bootstrap
resample events jointly and put a confidence interval directly on the DIFFERENCE. A
CI on the difference is the number that decides the question; two separate CIs do not.

Usage:
  python 17_split_control_summary.py                       # SE and MXE, 5 splits
  python 17_split_control_summary.py --types SE --nsplit 5
  python 17_split_control_summary.py --boot 5000 --plot
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

PROGRAM = ["TTN", "CAMK2D", "LDB3", "LMO7", "SORBS1", "PDLIM3", "DTNA", "MLIP", "TNNT2"]

# 15b, RBM20-program SE slopes, for the verdict at the end
EXTERNAL = {"Brigham": 0.587, "Fuwai": 0.440, "Stanford": 0.703}
MAGNET_SHARED_REPORTED = 0.97


def deming(x, y, lam):
    """Errors-in-variables slope; lam = var(err_y)/var(err_x)."""
    if len(x) < 3:
        return np.nan
    sxx, syy = np.var(x, ddof=1), np.var(y, ddof=1)
    sxy = np.cov(x, y, ddof=1)[0, 1]
    if abs(sxy) < 1e-12:
        return np.nan
    a = syy - lam * sxx
    return (a + np.sqrt(a * a + 4 * lam * sxy * sxy)) / (2 * sxy)


def load(path):
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return None
    t = pd.read_csv(path, sep="\t")
    need = {"key", "gene", "dPSI_adj", "FDR"}
    if need - set(t.columns):
        return None
    return t.drop_duplicates("key")


def stats_for(x, y, sx, sy):
    """r, OLS slope, Deming slope, direction concordance."""
    out = {"n": len(x)}
    if len(x) < 5:
        return {**out, "r": np.nan, "ols": np.nan, "deming": np.nan,
                "same": np.nan, "pct": np.nan}
    r = stats.pearsonr(x, y)[0]
    ols = stats.linregress(x, y).slope
    vx, vy = np.nanmean(sx ** 2), np.nanmean(sy ** 2)
    lam = float(vy / vx) if vx > 0 else 1.0
    same = int(np.sum(np.sign(x) == np.sign(y)))
    return {**out, "r": r, "ols": ols, "deming": deming(x, y, lam),
            "same": same, "pct": 100.0 * same / len(x), "lam": lam}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default="splitctl")
    ap.add_argument("--types", default="SE,MXE")
    ap.add_argument("--nsplit", type=int, default=5)
    ap.add_argument("--fdr", type=float, default=0.05)
    ap.add_argument("--min-dpsi", type=float, default=0.10)
    ap.add_argument("--boot", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--plot", action="store_true",
                    help="also write split_control_check.png if matplotlib is present")
    ap.add_argument("--out", default="split_control_summary")
    a = ap.parse_args()

    types = [t.strip() for t in a.types.split(",") if t.strip()]
    rng = np.random.default_rng(a.seed)
    lines = []

    def say(s=""):
        lines.append(s)
        print(s)

    rows = []
    qc = []

    for et in types:
        for s in range(1, a.nsplit + 1):
            d = load(f"{a.dir}/{et}_DCMvsA_{s}.all_events.tsv")
            ha = load(f"{a.dir}/{et}_HCMvsA_{s}.all_events.tsv")
            hb = load(f"{a.dir}/{et}_HCMvsB_{s}.all_events.tsv")
            if d is None or ha is None or hb is None:
                say(f"[skip] {et} split {s}: a table is missing or malformed "
                    f"(check {a.dir}/{et}_*_{s}.log)")
                continue

            # ---- QC: did the exclusion actually bite? -----------------------
            n_ctrl_a = float(np.nanmedian(ha["n_g2"])) if "n_g2" in ha else np.nan
            n_ctrl_b = float(np.nanmedian(hb["n_g2"])) if "n_g2" in hb else np.nan
            n_ctrl_d = float(np.nanmedian(d["n_g2"])) if "n_g2" in d else np.nan
            qc.append({"type": et, "split": s, "NF_in_DCM_arm": n_ctrl_d,
                       "NF_in_HCM_A": n_ctrl_a, "NF_in_HCM_B": n_ctrl_b})

            # ---- the SHARED event set -------------------------------------
            dsig = d[(d.FDR < a.fdr) & (d.dPSI_adj.abs() >= a.min_dpsi)]
            # rename BEFORE merging -- pandas `suffixes` would collide on the
            # second merge (left already carries dPSI_adj_A) and silently
            # duplicate column names.
            L = dsig[["key", "gene", "dPSI_adj", "se_adj"]].rename(
                columns={"dPSI_adj": "x", "se_adj": "sx"})
            A = ha[["key", "dPSI_adj", "se_adj"]].rename(
                columns={"dPSI_adj": "dPSI_adj_A", "se_adj": "se_adj_A"})
            B = hb[["key", "dPSI_adj", "se_adj"]].rename(
                columns={"dPSI_adj": "dPSI_adj_B", "se_adj": "se_adj_B"})
            m = L.merge(A, on="key").merge(B, on="key")
            m = m.dropna(subset=["x", "dPSI_adj_A", "dPSI_adj_B"])
            if len(m) < 10:
                say(f"[skip] {et} split {s}: only {len(m)} events in all three tables")
                continue

            for label, sub in (("background", m),
                               ("program", m[m.gene.isin(PROGRAM)])):
                if len(sub) < 5:
                    continue
                x = sub.x.values.astype(float)
                ya = sub.dPSI_adj_A.values.astype(float)
                yb = sub.dPSI_adj_B.values.astype(float)
                sx = sub.sx.values.astype(float)
                sa = sub.se_adj_A.values.astype(float)
                sb = sub.se_adj_B.values.astype(float)

                sh = stats_for(x, ya, sx, sa)
                dj = stats_for(x, yb, sx, sb)

                # paired bootstrap on the DIFFERENCE (same events resampled jointly)
                idx = np.arange(len(x))
                diffs = np.empty(a.boot)
                diffs[:] = np.nan
                for b in range(a.boot):
                    j = rng.choice(idx, size=len(idx), replace=True)
                    ds = deming(x[j], ya[j], np.nanmean(sa[j] ** 2) / max(np.nanmean(sx[j] ** 2), 1e-12))
                    dd = deming(x[j], yb[j], np.nanmean(sb[j] ** 2) / max(np.nanmean(sx[j] ** 2), 1e-12))
                    diffs[b] = ds - dd
                good = diffs[np.isfinite(diffs)]
                lo, hi = (np.percentile(good, [2.5, 97.5]) if len(good) > 50
                          else (np.nan, np.nan))

                # how often does A/B disagree in sign? a direct look at NF noise
                rows.append({
                    "type": et, "split": s, "set": label, "n": len(x),
                    "r_shared": sh["r"], "r_disjoint": dj["r"],
                    "deming_shared": sh["deming"], "deming_disjoint": dj["deming"],
                    "delta": sh["deming"] - dj["deming"],
                    "delta_lo": lo, "delta_hi": hi,
                    "pct_same_shared": sh["pct"], "pct_same_disjoint": dj["pct"],
                    "r_AB": stats.pearsonr(ya, yb)[0],
                })

    if not rows:
        sys.exit("[fatal] nothing to summarise -- check that 16 finished and that "
                 f"{a.dir}/ contains *_DCMvsA_*.all_events.tsv etc.")

    R = pd.DataFrame(rows)
    Q = pd.DataFrame(qc)

    # ---- QC first: if the split did not bite, nothing below means anything ---
    say("=" * 78)
    say("0. QC -- did the NF split actually take effect?")
    say()
    say(f"  {'type':<5}{'split':>6}{'NF in DCM arm':>15}{'NF in HCM(A)':>14}{'NF in HCM(B)':>14}")
    for _, q in Q.iterrows():
        say(f"  {q['type']:<5}{q['split']:>6}{q['NF_in_DCM_arm']:>15.0f}"
            f"{q['NF_in_HCM_A']:>14.0f}{q['NF_in_HCM_B']:>14.0f}")
    say()
    bad = Q[(Q.NF_in_HCM_A > 120) | (Q.NF_in_HCM_B > 120)]
    if len(bad):
        say("  *** WARNING: an HCM arm still has >120 controls. The --exclude did NOT")
        say("  *** bite -- both arms are effectively the shared one and the whole")
        say("  *** comparison is void. Check the [info] excluding ... lines in the logs.")
    else:
        say("  All arms near 83 controls: the split took effect.")
    say()
    say("  Correlation between the two HCM arms' dPSI (same cases, different controls):")
    for et in R.type.unique():
        for st in R.set.unique():
            sub = R[(R.type == et) & (R.set == st)]
            if len(sub):
                say(f"    {et:<5} {st:<11} r(A,B) = {sub.r_AB.mean():+.3f}")
    say("  EXPECT THIS TO BE HIGH (0.9+), not low. Both arms share the same 28 HCM")
    say("  cases, and with n_HCM=28 against n_NF=83 the case-group sampling error is")
    say("  ~3x the control-group one -- so the two arms share their DOMINANT noise")
    say("  source and the true effect on top of it. A high r(A,B) is NOT evidence")
    say("  that the split failed. The real proof the split bit is the systematic")
    say("  shared > disjoint gap in section 1: if --exclude had done nothing, the")
    say("  two arms would be identical and every delta would be exactly 0.")
    say()

    # ---- the answer ---------------------------------------------------------
    say("=" * 78)
    say("1. THE MEASUREMENT -- shared-control inflation of the Deming slope")
    say()
    say(f"  {'type':<5}{'set':<12}{'split':>6}{'n':>6}{'shared':>9}{'disjoint':>10}"
        f"{'delta':>9}{'95% CI on delta':>22}")
    for _, r in R.sort_values(["type", "set", "split"]).iterrows():
        ci = (f"[{r.delta_lo:+.3f}, {r.delta_hi:+.3f}]"
              if np.isfinite(r.delta_lo) else "--")
        say(f"  {r.type:<5}{r['set']:<12}{r.split:>6}{r.n:>6}"
            f"{r.deming_shared:>9.3f}{r.deming_disjoint:>10.3f}{r.delta:>9.3f}{ci:>22}")
    say()

    say("  ACROSS SPLITS (mean +/- sd; the splits share most events, so treat this")
    say("  as a stability check, NOT as five independent replicates):")
    say()
    say(f"  {'type':<5}{'set':<12}{'shared':>16}{'disjoint':>16}{'delta':>16}")
    agg = {}
    for (et, st), g in R.groupby(["type", "set"]):
        agg[(et, st)] = g
        say(f"  {et:<5}{st:<12}"
            f"{f'{g.deming_shared.mean():.3f} ± {g.deming_shared.std(ddof=1):.3f}':>16}"
            f"{f'{g.deming_disjoint.mean():.3f} ± {g.deming_disjoint.std(ddof=1):.3f}':>16}"
            f"{f'{g.delta.mean():.3f} ± {g.delta.std(ddof=1):.3f}':>16}")
    say()

    # ---- verdict ------------------------------------------------------------
    say("=" * 78)
    say("2. VERDICT")
    say()
    # judge on the BACKGROUND set: ~700 events with tight CIs. The program set
    # has ~20 events and its slope is unstable (MXE/program can swing past 2).
    key = ("SE", "background")
    g = agg.get(key)
    if g is None:
        g = agg.get(("SE", "program"))
        key = ("SE", "program")
    if g is None:
        say("  No SE result to judge. Read the table above directly.")
    else:
        dis = g.deming_disjoint.mean()
        sha = g.deming_shared.mean()
        dlt = g.delta.mean()
        anyci = g.dropna(subset=["delta_lo"])
        ci_excludes_zero = bool(len(anyci)) and bool((anyci.delta_lo > 0).all())
        ext_lo, ext_hi = min(EXTERNAL.values()), max(EXTERNAL.values())

        say(f"  Using {key[0]} / {key[1]}:")
        say(f"    shared-control slope   : {sha:.3f}   (MAGNet as originally computed, ~{MAGNET_SHARED_REPORTED})")
        say(f"    disjoint-control slope : {dis:.3f}   <-- the honest MAGNet number")
        say(f"    inflation              : {dlt:+.3f}"
            + ("   (95% CI excludes 0 in every split)" if ci_excludes_zero
               else "   (CI includes 0 in at least one split)"))
        say(f"    external cohorts (15b) : {ext_lo:.3f} - {ext_hi:.3f}"
            f"  (Brigham {EXTERNAL['Brigham']}, Fuwai {EXTERNAL['Fuwai']},"
            f" Stanford {EXTERNAL['Stanford']})")
        say()
        if not np.isfinite(dis):
            say("  Slope not estimable -- read the per-split table.")
        elif ext_lo - 0.10 <= dis <= ext_hi + 0.10:
            say("  >>> INTERPRETATION (a): the gap was largely the ARTEFACT.")
            say("  The disjoint-control MAGNet slope lands inside the range the three")
            say("  independent cohorts gave. Once the shared reference is removed,")
            say("  MAGNet and the external cohorts AGREE. The signature is shared at")
            say("  a uniform magnitude and there is no disease-stage attenuation to")
            say("  claim. ACTION: stop quoting the ~0.97 anywhere; report the disjoint")
            say("  value as MAGNet's slope and say why in the methods.")
        elif dis >= 0.85:
            say("  >>> INTERPRETATION (b): the gap is REAL.")
            say("  Removing the shared controls barely moves MAGNet's slope, so the")
            say("  0.44-0.70 seen externally is genuine attenuation -- myectomy /")
            say("  preserved-EF tissue carries a milder version of the same signature.")
            say("  ACTION: this is a quantitative severity axis. It partly rescues the")
            say("  LVEF null: severity shows up ACROSS disease stages even though it")
            say("  does not grade WITHIN end-stage DCM. Worth its own figure panel.")
        else:
            say("  >>> INTERPRETATION (c): BOTH contribute.")
            say(f"  Of the {MAGNET_SHARED_REPORTED - ext_lo:.2f}-{MAGNET_SHARED_REPORTED - ext_hi:.2f}"
                f" total gap, about {dlt:.2f} is the shared-control")
            say("  artefact and the remainder is real attenuation in the external")
            say("  cohorts. ACTION: report the decomposition explicitly -- it is a")
            say("  more honest and more interesting paragraph than either extreme,")
            say("  and it pre-empts the reviewer who spots the shared reference.")
        say()
        say("  TWO CAVEATS, both of which push the ABSOLUTE slopes down and neither")
        say("  of which touches the DELTA:")
        say("   1. Events are selected on MAGNet DCM significance, so regression to")
        say("      the mean pulls both arms down by the same amount.")
        say("   2. The Deming lambda is built from se_adj, which is the WITHIN-group")
        say("      standard error. It does not include the control-group mean")
        say("      sampling error, so the x-error is understated and Deming")
        say("      under-corrects. On simulated data with a planted slope of 0.60")
        say("      this design reads the disjoint arm at ~0.51.")
        say("  So the disjoint slope printed above is a LOWER BOUND, and the test")
        say("  'does it land in the external range' is mildly biased toward")
        say("  interpretation (a). The external slopes carry the same downward bias,")
        say("  which partly cancels -- but do not treat a borderline (a) as settled.")
        say("  The DELTA is immune to both, because the selection and the lambda")
        say("  construction are identical on the two sides. That is the whole reason")
        say("  this design was used instead of eyeballing 0.97 against 0.44.")
        say()
        say("  WHICH SET TO QUOTE: lead the DELTA with the background set -- it has")
        say("  ~1500 events and its CI is tight. Quote the PROGRAM set for continuity")
        say("  with 15b, but its n is ~20 and its CI will be wide; do not read a")
        say("  program-vs-background difference into that width.")

    say("=" * 78)

    R.to_csv(f"{a.out}.tsv", sep="\t", index=False, float_format="%.6g")
    with open(f"{a.out}.report.txt", "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"\n[write] {a.out}.tsv  {a.out}.report.txt", file=sys.stderr)

    if a.plot:
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            sub = R[R.set == ("program" if (R.set == "program").any() else "background")]
            fig, ax = plt.subplots(figsize=(6, 4))
            for i, (_, r) in enumerate(sub.iterrows()):
                ax.plot([0, 1], [r.deming_shared, r.deming_disjoint],
                        "-o", color="0.4", alpha=.7, ms=4)
            for v, c, lab in [(EXTERNAL["Fuwai"], "tab:red", "Fuwai"),
                              (EXTERNAL["Brigham"], "tab:orange", "Brigham"),
                              (EXTERNAL["Stanford"], "tab:green", "Stanford")]:
                ax.axhline(v, ls="--", lw=1, color=c)
                ax.text(1.02, v, lab, color=c, va="center", fontsize=8)
            ax.set_xticks([0, 1])
            ax.set_xticklabels(["shared\ncontrols", "disjoint\ncontrols"])
            ax.set_ylabel("Deming slope (HCM on DCM)")
            ax.set_title("Shared-control inflation, 5 random NF splits")
            ax.set_xlim(-.3, 1.45)
            fig.tight_layout()
            fig.savefig("split_control_check.png", dpi=200)
            print("[write] split_control_check.png", file=sys.stderr)
        except Exception as e:
            print(f"[warn] plot skipped: {e}", file=sys.stderr)


if __name__ == "__main__":
    main()
