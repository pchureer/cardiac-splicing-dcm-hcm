#!/usr/bin/env python3
"""
21_replication_v2.py -- replication concordance of the MAGNet discovery signature in
Colorado and LBNL, recomputed on the v2 tables.

WHY THIS EXISTS
---------------
The replication numbers currently in the handoff (Colorado SE 22/28 same-direction,
LBNL SE 16/25) came from the v1 pipeline, which had no coordinate dedup, no coverage
mask, no PSI-variance filter and no `key` column. Multi-transcript rMATS rows duplicate
the same event, so a v1 "28 events" may be fewer real events counted more than once,
and constant-PSI events produced nan correlations. Every other result in this project
is now on v2 tables; these two must be brought onto the same footing before any of
them go in a table together.

WHAT IT REPORTS, AND WHY EACH COLUMN IS THERE
---------------------------------------------
  retention   what fraction of MAGNet-significant events are even measurable in the
              replication cohort. A low number is a power/coverage statement, not a
              replication failure, and must be reported before the concordance.
  same dir    the sign test. Reported with an exact binomial p against 0.5.
  rho / r     effect-size correlation. For small cohorts this is far better powered
              than the sign test and should be quoted alongside it.
  Deming      how much of the MAGNet effect size the cohort recapitulates. Comparable
              with the HCM-cohort slopes in RESULT_three_hcm_cohorts.md.
  powered     events where the replication cohort's mde80 is small enough that it
              COULD have detected an effect the size of MAGNet's. Concordance among
              powered events is the honest replication number; concordance over all
              events conflates disagreement with lack of power.

Both the genome-wide set and the RBM20-program subset are reported, because the
handoff's historical numbers (22/28, 16/25) are program-subset numbers and the
comparison has to be like for like.

Usage:
  python 21_replication_v2.py --rep Colorado --suffix _Colorado_DCMvsNF
  python 21_replication_v2.py --rep LBNL --suffix _LBNL_DCMvsNF
"""
import argparse
import os

import numpy as np
import pandas as pd
from scipy import stats

TYPES = ["SE", "MXE", "A3SS", "A5SS", "RI"]
PROGRAM = ["TTN", "CAMK2D", "LDB3", "LMO7", "SORBS1", "PDLIM3", "DTNA", "MLIP", "TNNT2"]


def deming(x, y, lam):
    if len(x) < 5:
        return np.nan
    sxx, syy = np.var(x, ddof=1), np.var(y, ddof=1)
    sxy = np.cov(x, y, ddof=1)[0, 1]
    if abs(sxy) < 1e-12:
        return np.nan
    a = syy - lam * sxx
    return (a + np.sqrt(a * a + 4 * lam * sxy * sxy)) / (2 * sxy)


def load(path):
    if not (os.path.exists(path) and os.path.getsize(path)):
        return None
    t = pd.read_csv(path, sep="\t")
    if not {"key", "gene", "dPSI_adj", "FDR"} <= set(t.columns):
        return None
    t = t.drop_duplicates("key")
    if "mde80" not in t.columns and "se_adj" in t.columns:
        t["mde80"] = 2.802 * t["se_adj"]
    return t


def block(disc, rep, prog=False):
    d = disc[(disc.FDR < 0.05) & (disc.dPSI_adj.abs() >= 0.10)]
    if prog:
        d = d[d.gene.isin(PROGRAM)]
    n_disc = len(d)
    m = d[["key", "gene", "dPSI_adj", "se_adj", "mde80"]].rename(
            columns={"dPSI_adj": "x", "se_adj": "sx", "mde80": "mx"}).merge(
        rep[["key", "dPSI_adj", "se_adj", "mde80"]], on="key").dropna(
        subset=["x", "dPSI_adj"])
    out = {"n_disc": n_disc, "n_meas": len(m),
           "retention": 100 * len(m) / n_disc if n_disc else np.nan}
    if len(m) < 5:
        return {**out, "same": np.nan, "pct": np.nan, "pbin": np.nan, "rho": np.nan,
                "r": np.nan, "pr": np.nan, "deming": np.nan, "n_pow": 0,
                "pct_pow": np.nan, "pbin_pow": np.nan}
    x, y = m.x.values, m.dPSI_adj.values
    same = int(np.sum(np.sign(x) == np.sign(y)))
    lam = np.nanmean(m.se_adj.values ** 2) / max(np.nanmean(m.sx.values ** 2), 1e-12)
    pw = m[m.mde80 <= m.x.abs()]
    same_pw = int(np.sum(np.sign(pw.x.values) == np.sign(pw.dPSI_adj.values))) if len(pw) else 0
    return {**out,
            "same": same, "pct": 100 * same / len(m),
            "pbin": stats.binomtest(same, len(m), 0.5).pvalue,
            "rho": stats.spearmanr(x, y)[0],
            "r": stats.pearsonr(x, y)[0], "pr": stats.pearsonr(x, y)[1],
            "deming": deming(x, y, lam),
            "n_pow": len(pw),
            "pct_pow": 100 * same_pw / len(pw) if len(pw) else np.nan,
            "pbin_pow": stats.binomtest(same_pw, len(pw), 0.5).pvalue if len(pw) >= 5 else np.nan}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rep", required=True, help="cohort label, e.g. Colorado")
    ap.add_argument("--suffix", required=True,
                    help="replication table suffix, e.g. _Colorado_DCMvsNF")
    ap.add_argument("--disc-suffix", default="_v2filtered",
                    help="MAGNet discovery table suffix")
    ap.add_argument("--types", default=",".join(TYPES))
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    out = a.out or f"replication_v2_{a.rep}"
    types = [t.strip() for t in a.types.split(",") if t.strip()]

    lines = []
    def say(s=""):
        lines.append(s); print(s)

    say("=" * 96)
    say(f"REPLICATION OF THE MAGNet SIGNATURE IN {a.rep.upper()}  (v2 tables)")
    say(f"discovery: *{a.disc_suffix}.all_events.tsv    replication: *{a.suffix}.all_events.tsv")
    say()

    rows = []
    for prog in (False, True):
        tag = "RBM20 program" if prog else "genome-wide"
        say(f"--- {tag} ---")
        say(f"  {'type':<6}{'MAGNet sig':>11}{'measured':>10}{'retain':>8}"
            f"{'same dir':>11}{'binom p':>10}{'rho':>8}{'Deming':>9}"
            f"{'powered':>10}{'same|pow':>10}")
        for et in types:
            D = load(f"{et}{a.disc_suffix}.all_events.tsv")
            R = load(f"{et}{a.suffix}.all_events.tsv")
            if D is None or R is None:
                say(f"  {et:<6}  (missing table -- skipped)")
                continue
            b = block(D, R, prog)
            b.update({"type": et, "set": tag})
            rows.append(b)
            sd = f"{b['same']}/{b['n_meas']}" if b["n_meas"] else "--"
            pw = f"{b['n_pow']}/{b['n_meas']}" if b["n_meas"] else "--"
            say(f"  {et:<6}{b['n_disc']:>11}{b['n_meas']:>10}"
                f"{b['retention']:>7.0f}%{sd:>11}"
                f"{b['pbin']:>10.2g}{b['rho']:>8.3f}{b['deming']:>9.3f}"
                f"{pw:>10}{b['pct_pow']:>9.0f}%")
        say()

    pd.DataFrame(rows).to_csv(f"{out}.tsv", sep="\t", index=False, float_format="%.6g")

    say("HOW TO READ IT")
    say("  RETENTION FIRST. If only 40% of MAGNet-significant events are measurable in")
    say("  this cohort, the concordance is computed on a coverage-selected subset and")
    say("  must be described that way.")
    say()
    say("  Quote rho AND the sign test. For 15-37 sample cohorts the effect-size")
    say("  correlation is the better-powered statistic; a non-significant sign test")
    say("  alongside a significant rho is concordance with low power, NOT a failure.")
    say()
    say("  'same|pow' is the honest headline: direction agreement restricted to events")
    say("  this cohort could actually have detected. If 'powered' is a small fraction")
    say("  of 'measured', say so -- that is the real limitation of these cohorts.")
    say()
    say("  Compare the Deming slope with the HCM cohorts (0.44-0.70) and MAGNet's own")
    say("  disjoint arm (0.87-0.99). Colorado and LBNL are ~50bp like MAGNet, so if")
    say("  the read-length hypothesis in script 20 is right, THESE slopes should sit")
    say("  near 0.9 while the 101/150bp HCM cohorts sit lower. That makes this script")
    say("  an independent read on the same question.")
    say("=" * 96)

    with open(f"{out}.report.txt", "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"\n[write] {out}.tsv  {out}.report.txt")


if __name__ == "__main__":
    main()
