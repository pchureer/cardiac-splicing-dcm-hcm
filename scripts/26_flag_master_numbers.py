#!/usr/bin/env python3
"""
26_flag_master_numbers.py -- add a reliability verdict to every row of
master_numbers.tsv, derived from the reverse-replication test rather than assumed.

WHY
---
Colorado calls 1,462 MXE events significant (exon-level) against MAGNet's 231, from a
fifth of the samples. Script 25 showed this is NOT a false-positive-rate problem: every
one of 13 permutations returned exactly 0 significant events, matching MAGNet, and
per-sample junction coverage does not differ by group (13 vs 13, p=0.89).

Script 22 settled it instead. Select events on COLORADO significance and regress MAGNet
-- five times the samples, far more power -- on them:

    Colorado SE    rev n =  453   slope = 0.447    reproduces
    Colorado MXE   rev n = 1477   slope = 0.061    DOES NOT REPRODUCE
    MAGNet internal MXE  n = 147  slope = 0.725    (control: shares platform)

Colorado's SE events reproduce; its MXE events do not. Same cohort, same pipeline, same
samples -- so this is an MXE-from-single-end-51bp problem, not a Colorado problem. MXE
requires assigning reads to one of two mutually exclusive exons, which single-end short
reads cannot do reliably.

The Deming columns agree independently: Colorado MXE returns UNSEL 11.0 / PREC 15.0. A
slope of 11 means the cohort's dPSI variance swamps MAGNet's with almost no shared
covariance -- noise with a number attached. It is also why script 22's automated verdict
is meaningless for MXE: its means are contaminated by those values.

THE RULE APPLIED HERE
---------------------
  exclude    rev_n >= 100 and rev_slope < 0.15      enough events to judge, and they
                                                    do not reproduce in the largest cohort
  caution    rev_n >=  30 and rev_slope < 0.35      reproduces weakly
  reliable   rev_slope >= 0.35
  not_tested fewer than 30 events, or no reverse run for that event type

Two propagation rules, because the reverse test only ever compares against MAGNet DCM:
  * a verdict earned by one contrast of a cohort+event_type applies to that cohort's
    OTHER contrasts of the same event type (Colorado ICM-vs-NF MXE inherits the
    Colorado MXE exclusion -- same data, same quantification problem);
  * the discovery cohort itself is never judged by this test, since it is the reference.

Usage (login node, seconds):
  python 26_flag_master_numbers.py
  python 26_flag_master_numbers.py --exclude-slope 0.20   # stricter
"""
import argparse
import os

import numpy as np
import pandas as pd

TYPES = ["SE", "MXE", "A3SS", "A5SS", "RI"]

# fallback if symmetry_*.tsv are absent: the values measured 18 Sep 2026
FALLBACK = {
    ("Colorado", "SE"):  (453, 0.447), ("Colorado", "MXE"):  (1477, 0.061),
    ("LBNL", "SE"):      (11, 0.784),  ("LBNL", "MXE"):      (1, np.nan),
    ("Brigham", "SE"):   (4, np.nan),  ("Brigham", "MXE"):   (1, np.nan),
    ("Fuwai", "SE"):     (156, 0.405), ("Fuwai", "MXE"):     (99, 0.281),
    ("Stanford", "SE"):  (9, np.nan),  ("Stanford", "MXE"):  (4, np.nan),
}


def load_symmetry(et):
    f = f"symmetry_{et}.tsv"
    if not (os.path.exists(f) and os.path.getsize(f)):
        return None
    t = pd.read_csv(f, sep="\t")
    need = {"cohort", "rev_n", "rev_deming"}
    if need - set(t.columns):
        return None
    t = t[["cohort", "rev_n", "rev_deming", "unsel_deming", "prec_deming"]].copy()
    t["event_type"] = et
    return t


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--master", default="master_numbers.tsv")
    ap.add_argument("--exclude-slope", type=float, default=0.15)
    ap.add_argument("--caution-slope", type=float, default=0.35)
    ap.add_argument("--min-n", type=int, default=30)
    ap.add_argument("--exclude-n", type=int, default=100)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    out = a.out or a.master.replace(".tsv", "_flagged.tsv")

    M = pd.read_csv(a.master, sep="\t")

    sym = [s for s in (load_symmetry(et) for et in TYPES) if s is not None]
    if sym:
        S = pd.concat(sym, ignore_index=True)
        print(f"[info] reverse-replication values read from symmetry_*.tsv "
              f"({len(S)} cohort x event-type rows)")
    else:
        S = pd.DataFrame([{"cohort": c, "event_type": e, "rev_n": n,
                           "rev_deming": s, "unsel_deming": np.nan,
                           "prec_deming": np.nan}
                          for (c, e), (n, s) in FALLBACK.items()])
        print("[warn] symmetry_*.tsv not found -- using the values measured 18 Sep 2026. "
              "Re-run 22_symmetry_winners_curse.py for SE and MXE to refresh.")

    M = M.merge(S, on=["cohort", "event_type"], how="left")

    def verdict(r):
        if r.cohort == "MAGNet":
            return "reference", ""
        n, s = r.rev_n, r.rev_deming
        if not np.isfinite(n) or n < a.min_n or not np.isfinite(s):
            return "not_tested", "too few events for a reverse-replication test"
        if n >= a.exclude_n and s < a.exclude_slope:
            return "exclude", (f"own significant events do not reproduce in MAGNet "
                               f"(reverse slope {s:.3f}, n={int(n)})")
        if s < a.caution_slope:
            return "caution", (f"reproduces weakly in MAGNet "
                               f"(reverse slope {s:.3f}, n={int(n)})")
        return "reliable", ""

    v = M.apply(verdict, axis=1, result_type="expand")
    M["reliability"], M["flag_reason"] = v[0], v[1]

    # propagation 1: a verdict for cohort+event_type applies to that cohort's other
    # contrasts of the same event type -- same data, same quantification behaviour
    worst = {"exclude": 0, "caution": 1, "not_tested": 2, "reliable": 3, "reference": 4}
    for (coh, et), g in M.groupby(["cohort", "event_type"]):
        if coh == "MAGNet":
            continue
        best = min(g.reliability, key=lambda x: worst[x])
        if best in ("exclude", "caution"):
            reason = g.loc[g.reliability == best, "flag_reason"].iloc[0]
            idx = M[(M.cohort == coh) & (M.event_type == et)].index
            M.loc[idx, "reliability"] = best
            M.loc[idx, "flag_reason"] = reason + " [applies to all contrasts of this "\
                                                 "cohort and event type]"

    # independent signal: a blown-out Deming means se_adj is not trustworthy either
    bad_var = (M.unsel_deming.abs() > 2) | (M.prec_deming.abs() > 2)
    M.loc[bad_var & (M.cohort != "MAGNet"), "flag_reason"] = (
        M.loc[bad_var & (M.cohort != "MAGNet"), "flag_reason"].fillna("")
        + " | Deming variance ratio blown out (>2): se_adj unreliable")

    # the n_g1/n_g2 trap: these are POST-MASK medians, not cohort sizes
    M = M.rename(columns={"n_g1": "median_samples_g1", "n_g2": "median_samples_g2"})

    cols = ["cohort", "contrast", "event_type", "convention", "reliability",
            "tested", "significant", "genes_sig",
            "median_samples_g1", "median_samples_g2",
            "rev_n", "rev_deming", "tested_full", "tested_exon",
            "sig_full", "sig_exon", "collapse_pct", "flag_reason", "file"]
    cols = [c for c in cols if c in M.columns]
    M[cols].to_csv(out, sep="\t", index=False, float_format="%.4g")

    print()
    print("=" * 92)
    print("RELIABILITY OF EVERY ROW")
    print("=" * 92)
    print(f"  {'cohort':<10}{'contrast':<34}{'type':<6}{'sig':>7}{'verdict':>12}"
          f"{'rev slope':>11}{'rev n':>7}")
    for _, r in M.iterrows():
        rs = f"{r.rev_deming:.3f}" if np.isfinite(r.get("rev_deming", np.nan)) else "--"
        rn = f"{int(r.rev_n)}" if np.isfinite(r.get("rev_n", np.nan)) else "--"
        mark = " <<<" if r.reliability == "exclude" else ""
        print(f"  {r.cohort:<10}{r.contrast:<34}{r.event_type:<6}"
              f"{int(r.significant):>7}{r.reliability:>12}{rs:>11}{rn:>7}{mark}")
    print()
    print("  counts:")
    for k, n in M.reliability.value_counts().items():
        print(f"    {k:<12}{n:>4}")
    print()
    ex = M[M.reliability == "exclude"]
    if len(ex):
        print("  EXCLUDED -- do not put these counts in any table or figure:")
        for _, r in ex.iterrows():
            print(f"    {r.cohort} / {r.contrast} / {r.event_type}: "
                  f"{int(r.significant)} events")
        print(f"    reason: {ex.flag_reason.iloc[0]}")
    print()
    print("  NOTE: n_g1/n_g2 have been renamed median_samples_g1/g2. They are per-event")
    print("  medians AFTER the coverage mask, NOT cohort sizes -- MAGNet reads 125/131,")
    print("  not 166/166. Build Table 1 from the true cohort counts, never from these.")
    print(f"\n[write] {out}")


if __name__ == "__main__":
    main()
