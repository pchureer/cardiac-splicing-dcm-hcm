#!/usr/bin/env python3
"""
22_symmetry_winners_curse.py -- why is the cross-cohort slope ~0.56?

THE SITUATION
-------------
Every independent cohort recapitulates the MAGNet signature at a Deming slope near
0.56 AT 50 bp (Colorado 0.571, LBNL 0.525, Brigham 0.757, Fuwai 0.359, Stanford
0.578), while MAGNet compared with itself gives 0.92.

⚠ THOSE SIX FIGURES ARE 50 bp VALUES AND THIS SCRIPT RECOMPUTES THEM. Do not
quote them from this docstring. They are nevertheless transcribed into
82_build_s4.py (GENOME_WIDE), 62_build_s5.py and 98_build_replication_table.py
line 393; at 125 bp every one of them moves and all four copies need the
run's own output.

Disease stage, COHORT read length, control-group size and shared controls are all
excluded (scripts 16/17/19/20/21).

⚠ "read length" there means the REPLICATION COHORTS' read lengths — Colorado
~50 bp against Fuwai 150 bp and so on. It does NOT mean MAGNet's own
quantification read length, which none of those scripts varied and which the
125 bp promotion has since shown moves every cohort slope by +15% to +29%. The
exclusion stands for the axis it tested and says nothing about the axis it
could not test.

The remaining candidate is that the comparison itself is biased:

  WINNER'S CURSE. Every event enters the analysis BECAUSE it was significant in
  MAGNet. Selection on x truncates the x distribution at the extremes of its own
  noise, so MAGNet's effect sizes are upward-biased and any independent re-measurement
  regresses toward the truth. Deming regression corrects for measurement error in x;
  it does NOT and cannot correct for selection ON x.

THE TEST -- and a warning about the obvious version of it
---------------------------------------------------------
The intuitive test is "selection bias is symmetric, a real effect-size difference is
not": compare FORWARD (select on MAGNet, regress the cohort on MAGNet) against
REVERSE (select on the cohort, regress MAGNet on it).

THAT TEST DOES NOT WORK HERE, and this script was calibrated on synthetic data to
find out. With a planted TRUE transfer of exactly 1.00 and perfectly symmetric truth,
differing only in precision (MAGNet se=0.030, replication se=0.080), it returned:

    forward 0.949    reverse 0.674    asymmetry 0.275

Selection truncates the SELECTED cohort's own noise, so the noisier the cohort you
select on, the more of its extremes are noise and the more its partner regresses.
MAGNet is ~2.7x more precise than any replication cohort, so unequal precision
GENERATES forward/reverse asymmetry with no real difference in effect size at all.
The forward/reverse columns are therefore reported as diagnostics only. DO NOT read
asymmetry as evidence that something inflates MAGNet.

WHAT DOES WORK is the UNSELECTED slope -- every event measured in both cohorts, no
significance filter, so no selection to be cursed by. In the same simulation it
returned 0.999 against a planted 1.000, as did the precision-filtered version.

The second calibration result is the one that matters most for reading the real data:
with a true transfer of 1.00, winner's curse cost only 1.00 -> 0.95. MAGNet is simply
too precise for selection on FDR<0.05 and |dPSI|>=0.10 to cost much. So a real forward
slope of 0.56 already implies a true transfer near 0.56/0.95 = 0.59, and the
UNSELECTED column is the direct measurement of it.

  * unselected ~= 1.0    -> the 0.56 is selection. The signature transfers at full
                            magnitude; report the unselected slope and explain the bias.
  * unselected ~= 0.6    -> the attenuation is REAL and is not selection, not stage,
                            not read length, not control-group size. Report it as an
                            unexplained cross-cohort attenuation.

A selection-strength sweep is also run (|dPSI| thresholds 0.00-0.20 at FDR<0.05).
Winner's curse predicts the slope falls as the selection gets more extreme; in the
simulation it fell only 0.96 -> 0.92, so a steep fall in the real data would itself
be informative.

MAGNet-internal (DCM vs HCM on disjoint NF halves) is included as the control: it
shares platform and batch with the discovery arm, so if the asymmetric branch is right
it should behave differently from the five independent cohorts.

CAVEAT ON THE UNSELECTED SLOPE: with no selection, most events have a true effect near
zero, so the regression is dominated by noise and the Deming lambda matters much more.
The script therefore also reports the unselected slope restricted to events with a
minimum PRECISION (small se in both cohorts) rather than a minimum effect -- precision
filtering is not selection on the effect, so it does not reintroduce the curse.

Runs on the login node. No new phase 2, no rMATS.

Usage:
  python 22_symmetry_winners_curse.py
  python 22_symmetry_winners_curse.py --event-type MXE
"""
import argparse
import os

import numpy as np
import pandas as pd
from scipy import stats

PROGRAM = ["TTN", "CAMK2D", "LDB3", "LMO7", "SORBS1", "PDLIM3", "DTNA", "MLIP", "TNNT2"]


def deming(x, y, lam):
    if len(x) < 10:
        return np.nan
    sxx, syy = np.var(x, ddof=1), np.var(y, ddof=1)
    sxy = np.cov(x, y, ddof=1)[0, 1]
    if abs(sxy) < 1e-12:
        return np.nan
    a = syy - lam * sxx
    return (a + np.sqrt(a * a + 4 * lam * sxy * sxy)) / (2 * sxy)


def load(p):
    if not (os.path.exists(p) and os.path.getsize(p)):
        return None
    t = pd.read_csv(p, sep="\t")
    if not {"key", "gene", "dPSI_adj", "FDR", "se_adj"} <= set(t.columns):
        return None
    return t.drop_duplicates("key")[["key", "gene", "dPSI_adj", "FDR", "se_adj"]]


def fit(m, xcol, ycol, sxcol, sycol):
    """Regress ycol on xcol. Returns Deming, OLS, r, %same, n."""
    x = m[xcol].values.astype(float)
    y = m[ycol].values.astype(float)
    if len(x) < 10:
        return dict(n=len(x), deming=np.nan, ols=np.nan, r=np.nan, pct=np.nan)
    lam = (np.nanmean(m[sycol].values ** 2) /
           max(np.nanmean(m[sxcol].values ** 2), 1e-12))
    return dict(n=len(x), deming=deming(x, y, lam),
                ols=stats.linregress(x, y).slope,
                r=stats.pearsonr(x, y)[0],
                pct=100 * np.mean(np.sign(x) == np.sign(y)))


# --------------------------------------------------------------------------
# READ LENGTH. Three different rules apply inside PAIRS below — see the
# comments there — and `out` carries SUF because all three of this script's
# outputs derive from it ({out}.tsv, {out}_sweep.tsv, {out}.report.txt) and an
# unsuffixed 125 bp run would overwrite the 50 bp ones in place.
SUF = os.environ.get("TBL_SUFFIX", "")


def _require_internal(pairs):
    """The MAGNet-internal row needs 136's output. Refuse rather than skip.

    `load()` returns None for a missing file and the row drops out silently,
    which would leave 22's headline contrast — cohort slopes against MAGNet
    compared with itself — with nothing on one side of it.
    """
    if not SUF:
        return
    label, x, y = pairs[0]
    missing = [p for p in (x, y) if not os.path.exists(p)]
    if missing:
        raise SystemExit(
            "REFUSING — TBL_SUFFIX=%s but the MAGNet-internal arm is absent:\n"
            "%s\n\n"
            "  Build it with 136_splitctl_rl125.sbatch, which re-runs phase 2\n"
            "  against rmats_magnet%s/ and rmats_HCMvsNF%s/ on the SAME 83/83\n"
            "  NF partition (splitctl/nfA_1.ids, nfB_1.ids — reused verbatim,\n"
            "  never re-drawn):\n\n"
            "    ET=%s sbatch 136_splitctl_rl125.sbatch\n\n"
            "  Without it this script's headline contrast has no MAGNet side,\n"
            "  and running it anyway would compare 125 bp cohort slopes\n"
            "  against a 50 bp MAGNet-internal slope."
            % (SUF, "".join("    missing: %s\n" % m for m in missing),
               SUF, SUF, os.environ.get("ET_HINT", "SE")))
# --------------------------------------------------------------------------


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--event-type", default="SE")
    ap.add_argument("--fdr", type=float, default=0.05)
    ap.add_argument("--min-dpsi", type=float, default=0.10)
    ap.add_argument("--prec-quantile", type=float, default=0.25,
                    help="for the precision-filtered unselected slope: keep events in "
                         "the best this-fraction of se_adj in BOTH cohorts")
    ap.add_argument("--program", action="store_true",
                    help="restrict to the RBM20 program genes")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    ET = a.event_type
    # SUF LAST, after the _program variant, so the stem reads
    #   symmetry_SE_program_rl125
    out = a.out or f"symmetry_{ET}" + ("_program" if a.program else "") + SUF

    # label: (x table = MAGNet side, y table = the other cohort)
    # THREE RULES, NOT ONE.
    #
    # Row 1, MAGNet internal: BOTH sides are MAGNet, so both carry SUF. Built
    #   by 136_splitctl_rl125.sbatch on the same 83/83 NF partition.
    # Rows 2-6, the cohorts: the x side is the MAGNet discovery table and
    #   carries SUF; the y side is the replication cohort's own quantification
    #   and must NOT — no replication cohort was re-quantified. Suffixing it
    #   would read a file that does not exist; not suffixing the x side would
    #   compare 50 bp discovery against a 125 bp claim.
    PAIRS = [
        ("MAGNet internal", f"splitctl/{ET}_DCMvsA_1{SUF}.all_events.tsv",
                            f"splitctl/{ET}_HCMvsB_1{SUF}.all_events.tsv"),
        ("Colorado",  f"{ET}_v2filtered{SUF}.all_events.tsv", f"{ET}_Colorado_DCMvsNF.all_events.tsv"),
        ("LBNL",      f"{ET}_v2filtered{SUF}.all_events.tsv", f"{ET}_LBNL_DCMvsNF.all_events.tsv"),
        ("Brigham",   f"{ET}_v2filtered{SUF}.all_events.tsv", f"Brigham_{ET}_HCM.all_events.tsv"),
        ("Fuwai",     f"{ET}_v2filtered{SUF}.all_events.tsv", f"Fuwai_{ET}_HCM.all_events.tsv"),
        ("Stanford",  f"{ET}_v2filtered{SUF}.all_events.tsv", f"Stanford_{ET}_HCM.all_events.tsv"),
    ]
    os.environ.setdefault("ET_HINT", ET)
    _require_internal(PAIRS)

    lines = []
    def say(s=""):
        lines.append(s); print(s)

    say("=" * 104)
    say(f"SYMMETRY TEST -- is the ~0.56 cross-cohort slope winner's curse?   "
        f"[{ET}{', RBM20 program' if a.program else ''}]")
    say(f"selection: FDR<{a.fdr} and |dPSI|>={a.min_dpsi}")
    say()

    rows, sweeps = [], []
    for lab, fx, fy in PAIRS:
        X, Y = load(fx), load(fy)
        if X is None or Y is None:
            say(f"[skip] {lab}: missing {[p for p, t in ((fx,X),(fy,Y)) if t is None]}")
            continue
        m = X.rename(columns={"dPSI_adj": "x", "FDR": "fx", "se_adj": "sx"}).merge(
            Y.rename(columns={"dPSI_adj": "y", "FDR": "fy", "se_adj": "sy"}),
            on="key", suffixes=("", "_y")).dropna(subset=["x", "y"])
        if a.program:
            m = m[m.gene.isin(PROGRAM)]
        if len(m) < 20:
            say(f"[skip] {lab}: only {len(m)} events measured in both")
            continue

        fwd = m[(m.fx < a.fdr) & (m.x.abs() >= a.min_dpsi)]
        rev = m[(m.fy < a.fdr) & (m.y.abs() >= a.min_dpsi)]
        # precision filter: best quartile of se in BOTH, no selection on effect size
        qx, qy = m.sx.quantile(a.prec_quantile), m.sy.quantile(a.prec_quantile)
        prec = m[(m.sx <= qx) & (m.sy <= qy)]

        r = {"cohort": lab, "n_common": len(m)}
        for tag, sub, xc, yc, sxc, syc in (
                ("fwd",  fwd,  "x", "y", "sx", "sy"),   # select on MAGNet, y on x
                ("rev",  rev,  "y", "x", "sy", "sx"),   # select on other,  x on y
                ("unsel", m,   "x", "y", "sx", "sy"),
                ("prec", prec, "x", "y", "sx", "sy")):
            f = fit(sub, xc, yc, sxc, syc)
            r.update({f"{tag}_{k}": v for k, v in f.items()})
        rows.append(r)

        for thr in (0.00, 0.05, 0.10, 0.15, 0.20):
            s = m[(m.fx < a.fdr) & (m.x.abs() >= thr)]
            f = fit(s, "x", "y", "sx", "sy")
            sweeps.append({"cohort": lab, "thr": thr, **f})

    if not rows:
        raise SystemExit("[fatal] no cohort pairs could be built -- check the filenames")
    R, S = pd.DataFrame(rows), pd.DataFrame(sweeps)
    R.to_csv(f"{out}.tsv", sep="\t", index=False, float_format="%.6g")
    S.to_csv(f"{out}_sweep.tsv", sep="\t", index=False, float_format="%.6g")

    say("1. SLOPES.  UNSEL and PREC are the trustworthy columns -- no significance")
    say("   selection, so no winner's curse. FORWARD/REVERSE are diagnostics only:")
    say("   unequal precision generates asymmetry by itself (planted 1.00 gave")
    say("   forward 0.95 / reverse 0.67 in simulation), so do NOT read asymmetry as")
    say("   evidence that something inflates MAGNet.")
    say()
    say(f"  {'cohort':<17}{'common':>8}{'fwd n':>8}{'FORWARD':>9}{'rev n':>8}"
        f"{'REVERSE':>9}{'fwd-rev':>9}{'unsel n':>9}{'UNSEL':>8}{'prec n':>8}{'PREC':>7}")
    for _, r in R.iterrows():
        say(f"  {r.cohort:<17}{r.n_common:>8}{r.fwd_n:>8}{r.fwd_deming:>9.3f}"
            f"{r.rev_n:>8}{r.rev_deming:>9.3f}{r.fwd_deming - r.rev_deming:>9.3f}"
            f"{r.unsel_n:>9}{r.unsel_deming:>8.3f}{r.prec_n:>8}{r.prec_deming:>7.3f}")
    say()

    say("2. SELECTION-STRENGTH SWEEP (forward direction, FDR<%.2f, rising |dPSI|)" % a.fdr)
    say("   winner's curse predicts the slope FALLS as the selection gets more extreme")
    say()
    thrs = sorted(S.thr.unique())
    say("  " + f"{'cohort':<17}" + "".join(f"{f'>={t:.2f}':>10}" for t in thrs))
    for coh in R.cohort:
        sub = S[S.cohort == coh].set_index("thr")
        say("  " + f"{coh:<17}" +
            "".join(f"{sub.loc[t,'deming']:>10.3f}" if t in sub.index else f"{'--':>10}"
                    for t in thrs))
    say()
    say("  (n at each threshold)")
    for coh in R.cohort:
        sub = S[S.cohort == coh].set_index("thr")
        say("  " + f"{coh:<17}" +
            "".join(f"{int(sub.loc[t,'n']):>10}" if t in sub.index else f"{'--':>10}"
                    for t in thrs))
    say()

    # ---- verdict ----------------------------------------------------------
    ind = R[R.cohort != "MAGNet internal"]
    say("=" * 104)
    say("3. VERDICT")
    say()
    if ind.empty:
        say("  no independent cohorts loaded -- read the table above directly")
    else:
        f_, r_, u_, p_ = (ind.fwd_deming.mean(), ind.rev_deming.mean(),
                          ind.unsel_deming.mean(), ind.prec_deming.mean())
        say(f"  independent cohorts, mean over {len(ind)}:")
        say(f"    forward (select on MAGNet)   {f_:.3f}")
        say(f"    reverse (select on cohort)   {r_:.3f}")
        say(f"    unselected                   {u_:.3f}")
        say(f"    precision-filtered           {p_:.3f}")
        say()
        say("  The verdict keys on UNSELECTED (calibrated: planted 1.000 -> read 0.999).")
        say("  FORWARD/REVERSE asymmetry is NOT used -- unequal precision creates it.")
        say()
        if not np.isfinite(u_):
            say("  >>> unselected slope not estimable; read the per-cohort table.")
        elif u_ > 0.80:
            say("  >>> WINNER'S CURSE. Without selection the slope is near 1, so the")
            say("  signature transfers at FULL magnitude and the 0.56 is created by")
            say("  selecting events on MAGNet significance.")
            say("  ACTION: report the unselected (or precision-filtered) slope as the")
            say("  effect-size transfer, state the selection bias explicitly in the")
            say("  methods, and drop all language about the signature being 'weaker'")
            say("  in the independent cohorts.")
        elif u_ < 0.75:
            say("  >>> REAL ATTENUATION. The unselected slope is well below 1, so this is")
            say("  NOT selection bias -- and scripts 19/20/21 already excluded disease")
            say("  stage, read length and control-group size.")
            say("  ACTION: report an unexplained cross-cohort attenuation of roughly")
            say(f"  {u_:.2f}. Attribute it to nothing. The replication claim rests on")
            say("  direction concordance (75-100%, all p<0.05), which is unaffected.")
            say("  A candidate worth naming in the limitations: cohort-specific")
            say("  technical variation (prep chemistry, RIN/degradation, batch) that")
            say("  no cross-cohort normalisation of PSI can remove.")
        else:
            say("  >>> INTERMEDIATE. Part selection, part real. Report the unselected")
            say(f"  slope ({u_:.2f}) as the effect-size transfer and note that selection")
            say("  on discovery significance accounts for some of the apparent gap.")
        say()
        say(f"  Sanity check on the precision-filtered column ({p_:.3f}): it should")
        say("  agree with UNSELECTED. If they disagree by more than ~0.1, the Deming")
        say("  lambda is misestimated and neither number should be quoted.")
        say()
        say("  Compare the MAGNet-internal row: it shares platform and batch with the")
        say("  discovery arm. If it stays high where the independent cohorts fall, that")
        say("  is the shared-platform component, separate from selection.")
    say("=" * 104)

    with open(f"{out}.report.txt", "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"\n[write] {out}.tsv  {out}_sweep.tsv  {out}.report.txt")


if __name__ == "__main__":
    main()
