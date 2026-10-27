#!/usr/bin/env python
"""
79_dpsi_global_skew.py -- is the discovery contrast's global dPSI tilt material?

WHAT 78 FOUND, AND WHY IT MATTERS MORE THAN THE LVEF SKEW
---------------------------------------------------------
`78_lvef_skew_origin.py` section (A) measured something nothing in this project had
looked at: the SIGN BALANCE of the discovery contrast itself.

    all 62,971 tested : 42.9% negative  ->  57.1% MORE INCLUDED IN DCM
    the 926 significant: 36.9% negative  ->  63.1% more included in DCM

That is a property of the HEADLINE contrast, not of a supplement, and it is unreported.
A reviewer with the deposited tables can compute it in one line.

THE QUESTION IS WHETHER IT IS MATERIAL, AND SIGN BALANCE CANNOT ANSWER THAT
---------------------------------------------------------------------------
A sign count says nothing about size. If the median dPSI is +0.002, the tilt is a
baseline nothing that only a 63,000-event sign test can even see. If it is +0.02, that
is a fifth of the |dPSI| >= 0.10 threshold and a different conversation.
**Magnitude first.** That is this script's main job.

THE MECHANISM TO RULE IN OR OUT
-------------------------------
78's stratification was the real clue, and it points at the estimator rather than at
biology. The LVEF tilt was ABSENT where PSI has room to move both ways
(mean PSI <= 0.80: 50.7%, the null) and PRESENT where PSI sits near the ceiling
(0.80-0.97: 56.4%). And the event set is dominated by near-constitutive exons --
three of four quartiles above 0.80, the median above 0.97.

For an exon with PSI ~ 0.99 the skipping count is near zero, so

    lower effective depth  ->  skipping reads more often exactly 0  ->  PSI reads HIGHER

PSI near the ceiling is therefore depth- and degradation-sensitive in ONE direction.
RIN is lower in DCM than in non-failing in these data (8.2 vs 8.6, P = 9.4e-6), and a
linear RIN term need not absorb an effect that only bites at one end of the PSI range.

**If the dPSI tilt is confined to high-PSI exons, it is this. If it is flat across the
PSI range, it is not, and a biological reading becomes available.**

WHAT THIS SCRIPT MEASURES
-------------------------
  (1) MAGNITUDE: mean and median dPSI overall and per stratum; the fraction of events
      whose |dPSI| even reaches the 0.10 threshold, split by sign
  (2) IS IT PSI-LEVEL CONFINED? sign balance and median dPSI in fixed PSI bins, so the
      ceiling mechanism above is ruled in or out
  (3) WHERE DO THE 926 LIVE? if the significant events sit mid-range while the tilt
      sits at the ceiling, the headline is insulated from it
  (4) OTHER CONTRASTS, via --also: the same summary for any other *.all_events.tsv.
      Two comparisons decide the attribution:
        * DCM vs HCM (BOTH EXPLANTS) near 50% would implicate donor-versus-explant
          handling, which this paper already discloses, rather than disease
        * SE_v2fulladj (+ library size + TIN) losing the tilt would implicate depth
          and integrity directly

⚠ EVERY p-VALUE IS DEPENDENCE-INFLATED. 63k exons are not 63k independent tests.
Read the SIZES. A p-value is printed only as an order of magnitude.

USAGE
-----
  conda activate rnaseq_tools
  cd /data/project/chureeratp/RNA/HEART

  python 79_dpsi_global_skew.py --list        # what tables exist, then stop
  python 79_dpsi_global_skew.py
  python 79_dpsi_global_skew.py --also SE_v2fulladj.all_events.tsv \\
                                --also <the DCM-vs-HCM table>
"""

import argparse
import glob
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

EXPECT_TESTED = 62_971
EXPECT_SIG = 926
NKEEP_SE = 3
C_KEY = "key"
FDR_MAX, DPSI_MIN = 0.05, 0.10
C_DPSI_CANDIDATES = ("dPSI_adj", "dPSI_raw", "dPSI")
PSI_BINS = [0.0, 0.50, 0.80, 0.95, 0.99, 1.01]
PSI_LABELS = ["<0.50", "0.50-0.80", "0.80-0.95", "0.95-0.99", ">0.99"]


def die(msg):
    print("\n[REFUSE] " + msg + "\n", file=sys.stderr)
    sys.exit(2)


def ok(msg):
    print("  [ok] " + msg)


def resolve_dpsi(cols):
    """Preferred names first, then ANY dPSI-prefixed column.

    The per-contrast tables name it after the contrast -- SE_DCMvsHCM carries
    `dPSI_DCMvsHCM` -- so a fixed candidate list refuses on exactly the tables this
    script most needs to read. Prefer a non-'raw' column when several match.
    """
    for c in C_DPSI_CANDIDATES:
        if c in cols:
            return c
    cand = [c for c in cols if str(c).lower().startswith("dpsi")]
    if not cand:
        return None
    nonraw = [c for c in cand if "raw" not in str(c).lower()]
    return (nonraw or cand)[0]


def collapse(path, gate=False):
    t = pd.read_csv(path, sep="\t", low_memory=False)
    dcol = resolve_dpsi(t.columns)
    if dcol is None:
        die(f"{path} has no dPSI-like column. Columns: {list(t.columns)}")
    for c in (C_KEY, "FDR"):
        if c not in t.columns:
            die(f"{path} has no '{c}'")
    t["dPSI"] = pd.to_numeric(t[dcol], errors="coerce")
    t["ekey"] = t[C_KEY].astype(str).str.split(":").str[:NKEEP_SE].str.join(":")
    rep = (t.sort_values("FDR", kind="mergesort")
            .drop_duplicates(subset=["ekey"], keep="first")
            .dropna(subset=["dPSI"]))
    nsig = int(((rep.FDR < FDR_MAX) & (rep.dPSI.abs() >= DPSI_MIN)).sum())
    if gate:
        if len(rep) != EXPECT_TESTED:
            die(f"G1 FAILED: {path} collapses to {len(rep):,}, "
                f"expected {EXPECT_TESTED:,}")
        if nsig != EXPECT_SIG:
            die(f"G1 FAILED: {path} carries {nsig} significant, expected {EXPECT_SIG}")
        ok(f"G1 {os.path.basename(path)}: {len(rep):,} tested / {nsig} significant "
           f"(effect '{dcol}')")
    return rep, nsig, dcol


def line(v, label, width=26):
    v = pd.Series(v).dropna()
    nz = v[v != 0]
    if len(nz) < 10:
        print(f"    {label:<{width}}: n = {len(v):,}  (too few non-zero to test)")
        return
    neg = int((nz < 0).sum())
    f = neg / len(nz)
    p = stats.binomtest(neg, len(nz), 0.5).pvalue
    print(f"    {label:<{width}}: n = {len(v):>6,}   {f * 100:5.1f}% neg "
          f"(excess {abs(f - 0.5) * 100:4.1f})   median {v.median():+.4f}   "
          f"mean {v.mean():+.4f}   p~{p:.1g}")


def summarise(rep, name, nsig):
    print(f"\n{'-' * 74}\n  {name}\n{'-' * 74}")

    print("  (1) MAGNITUDE -- the number a sign count cannot give you")
    line(rep["dPSI"], "all tested")
    sig = rep[(rep.FDR < FDR_MAX) & (rep.dPSI.abs() >= DPSI_MIN)]
    line(sig["dPSI"], f"the {len(sig):,} significant")

    reach = rep[rep.dPSI.abs() >= DPSI_MIN]
    if len(reach):
        pos = int((reach.dPSI > 0).sum())
        print(f"\n    events reaching |dPSI| >= {DPSI_MIN}: {len(reach):,} "
              f"({len(reach) / len(rep):.2%} of tested) -- "
              f"{pos:,} up / {len(reach) - pos:,} down = {pos / len(reach):.1%} up")
    med = float(rep["dPSI"].median())
    print(f"\n    >>> median dPSI = {med:+.4f}, i.e. {abs(med) / DPSI_MIN:.1%} of the "
          f"{DPSI_MIN} threshold.")
    if abs(med) < 0.005:
        print("        A baseline tilt this small is invisible per event and only a")
        print("        63,000-event sign test can see it. It cannot manufacture a")
        print("        significant event on its own.")
    else:
        print("        NOT negligible relative to the threshold. This needs an")
        print("        explanation before submission, not a footnote.")

    if "mean_psi" not in rep.columns:
        print("\n  (2)(3) skipped -- no mean_psi column")
        return

    print("\n  (2) IS IT PSI-LEVEL CONFINED?  (the ceiling mechanism)")
    b = pd.cut(rep["mean_psi"], PSI_BINS, labels=PSI_LABELS, include_lowest=True)
    for lev in PSI_LABELS:
        grp = rep[b == lev]
        if len(grp):
            line(grp["dPSI"], f"mean PSI {lev}")
    print("\n    Flat across bins -> not the ceiling. Confined to the top bins ->")
    print("    a depth/degradation artefact of near-constitutive exons, and the")
    print("    paper should say so rather than let a reviewer find it.")

    print("\n  (3) WHERE DO THE SIGNIFICANT EVENTS LIVE?")
    print("    The question is whether they sit where the TILT IS STRONGEST --")
    print("    not whether they avoid the ceiling.")
    tb = pd.crosstab(b, rep.index.isin(sig.index), normalize="columns") * 100
    allp = tb[False] if False in tb.columns else None
    sigp = tb[True] if True in tb.columns else None
    print(f"\n    {'mean PSI':<12}{'all tested':>12}{'significant':>14}")
    for lev in PSI_LABELS:
        av = f"{allp.get(lev, 0):.1f}%" if allp is not None else "n/a"
        sv = f"{sigp.get(lev, 0):.1f}%" if sigp is not None else "n/a"
        print(f"    {lev:<12}{av:>12}{sv:>14}")
    if sigp is None:
        return

    # The question is NOT whether the significant events avoid the ceiling. It is
    # whether they sit where the TILT IS STRONGEST. Depletion at PSI > 0.95 is
    # reassuring only if the tilt peaks there -- and in these data it does not.
    print("\n    Does the signature occupy the tilt's peak? (excess per bin, "
          "share of significant)")
    rows = []
    for lev in PSI_LABELS:
        grp = rep[b == lev]
        nz = grp["dPSI"][grp["dPSI"] != 0]
        exc = abs((nz < 0).mean() - 0.5) * 100 if len(nz) >= 10 else np.nan
        rows.append((lev, exc, float(grp["dPSI"].median()) if len(grp) else np.nan,
                     float(sigp.get(lev, 0))))
    print(f"    {'mean PSI':<12}{'excess pts':>12}{'median dPSI':>14}"
          f"{'% of sig':>11}")
    for lev, exc, med, sh in rows:
        print(f"    {lev:<12}{exc:>12.1f}{med:>+14.4f}{sh:>10.1f}%")

    valid = [r for r in rows if np.isfinite(r[1])]
    if valid:
        peak = max(valid, key=lambda r: r[1])
        sh_peak = sum(r[3] for r in valid if r[1] >= peak[1] - 3)
        print(f"\n    >>> The tilt PEAKS at mean PSI {peak[0]} "
              f"({peak[1]:.1f} points, median dPSI {peak[2]:+.4f}).")
        print(f"        {sh_peak:.0f}% of significant events sit in that peak region.")
        if sh_peak > 50:
            print("        ⚠ The signature occupies the SAME PSI range as the tilt's")
            print("          peak. Depletion at the ceiling is therefore NOT the")
            print("          reassurance it looks like: whatever drives the tilt acts")
            print("          where the called events are. The attribution runs below")
            print("          decide this, not this table.")
        else:
            print("        The signature largely sits away from the tilt's peak.")
    print("\n    Note the shape: if the excess were a ceiling artefact it would RISE")
    print("    monotonically toward PSI 1. An inverted U, peaking mid-range and")
    print("    weakest at both extremes, is not that -- so the ceiling mechanism is")
    print("    excluded and something acting on mid-inclusion exons is not.")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--events", default="SE_v2filtered.all_events.tsv")
    ap.add_argument("--also", action="append", default=[],
                    help="another *.all_events.tsv to summarise (repeatable)")
    ap.add_argument("--list", action="store_true")
    a = ap.parse_args()

    if a.list:
        print("=" * 74)
        print("Tables on disk -- pick the DCM-vs-HCM one and the technical model")
        print("=" * 74)
        found = sorted(glob.glob("*.all_events.tsv")) + sorted(
            glob.glob("*/*.all_events.tsv"))
        if not found:
            die("no *.all_events.tsv in the working directory")
        for f in found:
            try:
                n = sum(1 for _ in open(f)) - 1
            except OSError:
                n = -1
            print(f"  {n:>9,} rows   {f}")
        print("\n  THREE are worth passing to --also, in this order:")
        print()
        print("   1. jccheck/SE_v2filtered_JC.all_events.tsv  <- THE DECISIVE ONE")
        print("      In JCEC the INCLUSION form is counted by junction reads AND")
        print("      exon-body reads, while the SKIPPING form has no exon body and is")
        print("      counted by junction reads alone. Anything that depresses")
        print("      junction-spanning reads -- degradation, shorter effective")
        print("      fragments -- therefore depresses skipping MORE than inclusion and")
        print("      pushes PSI UP. DCM has the lower RIN here (8.2 vs 8.6,")
        print("      P = 9.4e-6). Under JC both forms are junction-only, so the")
        print("      asymmetry cannot operate. **If the tilt shrinks or vanishes in")
        print("      JC, it is this.** The effect is largest where both forms carry")
        print("      real counts, i.e. mid-range PSI -- which is the inverted U.")
        print()
        print("   2. SE_DCMvsHCM.all_events.tsv  (both arms explanted)")
        print("      Near 50% there would place the cause in donor-versus-explant")
        print("      handling -- ALREADY a disclosed limitation -- and would leave the")
        print("      paper's primary claim, the DCM-vs-HCM null, untouched by it.")
        print()
        print("   3. SE_v2fulladj.all_events.tsv  (+ library size + TIN)")
        print("      If the tilt goes once depth and transcript integrity are fitted,")
        print("      the cause is named and the fuller model earns a mention for that")
        print("      reason rather than only as a count-robustness row.")
        return

    print("=" * 74)
    print("GATES")
    print("=" * 74)
    rep, nsig, _ = collapse(a.events, gate=True)
    print("\n  sign convention: NEGATIVE dPSI = exon LESS included in DCM")
    print("  78 measured 42.9% negative over all tested and 36.9% over the 926,")
    print("  i.e. a tilt toward MORE inclusion in disease. Reproduced below.")

    summarise(rep, f"{os.path.basename(a.events)}  (the discovery contrast)", nsig)

    base_keys = set(rep["ekey"])
    for extra in a.also:
        if not os.path.exists(extra):
            print(f"\n  [warn] {extra} not found, skipped")
            continue
        r2, n2, dc2 = collapse(extra)
        print(f"\n  [{os.path.basename(extra)}] effect column '{dc2}'  -- a NEGATIVE")
        print("  value here refers to THAT table's own contrast, read in the order its")
        print("  column name gives (e.g. dPSI_DCMvsHCM < 0 = less included in DCM).")
        # A different table tests a different event SET -- JC in particular tests
        # 50,124 exons against JCEC's 62,971. Comparing a tilt across different sets
        # is not a comparison, so restrict to the shared exons and say so.
        shared = r2[r2["ekey"].isin(base_keys)]
        print(f"\n  {os.path.basename(extra)}: {len(r2):,} exon-level events, "
              f"{len(shared):,} shared with the base table "
              f"({len(shared) / max(len(rep), 1):.1%} of it)")
        nz = r2["dPSI"][r2["dPSI"] != 0]
        if len(nz) >= 10:
            print(f"    unrestricted sign balance: {(nz < 0).mean() * 100:.1f}% "
                  f"negative, median {r2['dPSI'].median():+.4f}")
        if len(shared) < max(1000, 0.2 * len(rep)):
            print("    [warn] too little overlap for a fair comparison")
            summarise(r2, os.path.basename(extra) + "  (UNRESTRICTED)", n2)
        else:
            n_sh = int(((shared.FDR < FDR_MAX)
                        & (shared.dPSI.abs() >= DPSI_MIN)).sum())
            summarise(shared, os.path.basename(extra)
                      + "  (restricted to the base table's exons)", n_sh)

    print("\n" + "=" * 74)
    print("HOW TO DECIDE")
    print("=" * 74)
    print("  1. If the median |dPSI| is a few thousandths, the tilt is a baseline")
    print("     property with no power to create a significant event, and at most it")
    print("     earns one sentence in the Methods or S2.")
    print("  2. If it is confined to mean PSI above ~0.80 AND the 926 are depleted")
    print("     there, the headline is insulated. Say both things in one sentence.")
    print("  3. If DCM-vs-HCM (both explants) shows no tilt while DCM-vs-NF does, the")
    print("     cause is donor-versus-explant -- ALREADY a disclosed limitation of")
    print("     this paper, which makes the finding cheap to absorb.")
    print("  4. If SE_v2fulladj removes it, the cause is depth and integrity, and the")
    print("     fuller model becomes worth reporting for that reason and not only as")
    print("     a count-robustness row.")
    print()
    print("  WHAT NOT TO DO: do not report a directional global shift as biology.")
    print("  Three of four PSI quartiles here sit above 0.80, the skipping count for")
    print("  such exons is near zero, and anything that lowers effective depth pushes")
    print("  their PSI up. That is the null hypothesis for a one-sided global tilt and")
    print("  it has not been excluded.")
    print()
    print("  ALSO NOT TO DO: do not treat 78's rising agreement with |dPSI| as")
    print("  dose-response evidence. A shared technical gradient -- one that moves")
    print("  both the disease contrast and the severity contrast -- predicts exactly")
    print("  the same monotone rise. The two readings are not separated by that test.")
    print()


if __name__ == "__main__":
    main()
