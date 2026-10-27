#!/usr/bin/env python3
"""
152_gapdh_rank_recount.py -- the GAPDH floor's rank, re-counted at 125 bp.

REVISED 22 Oct for the actual schema. The first version looked for a GAPDH ROW
and refused when it found none. The files instead carry:

    event_type  model  matching  gene  delta  floor  margin
    ci_lo  ci_hi  boot_p  holm_p  n_target  n_bg

160 rows = 2 classes x 4 models x 2 matchings x 10 candidate genes. **The floor
is a COLUMN, not a row**: `floor` is GAPDH's own delta, repeated within each
group, and `margin = delta - floor`. So the floor's rank among the eleven genes
scored is

    rank = 1 + (delta > floor).sum()          n = 10 candidates + 1

The first version returning "no rank" rather than a wrong one was the correct
failure, but it could not do the job.

WHY THIS IS NOT A NUMBER SWAP
-----------------------------
Three manuscript claims are **counts of a ranking**:

  [081] "In mutually exclusive exons the housekeeping control outscored nine of
         the ten candidates, RBM20's +0.009 the single exception ... In skipped
         exons the floor sits mid-table, fifth of the eleven genes scored."
  [096] "It outranks all but one of the regulators the field nominates."
  [007] "neither RBM20 nor nine further candidate regulators, QKI included, was
         distinguishable from a GAPDH negative control."

`RESULT_23_58_rbp_scan_rl125.md` records the largest 125 bp DCM/cell margins as
**PTBP1** (+0.0441 SE, +0.1027 MXE), where the manuscript names **QKI** (+0.033
SE) and **RBM20** (+0.009 MXE). The leader moves, so the rank must be re-counted.

WHICH CELLS ARE THE PAPER'S
---------------------------
`model = DCM/cell` and `matching = 3-D match (+|disease beta|)`.

⚠ The `matching` column labels the other scheme "2-D match (as published)", but
the manuscript's own Methods [052] describe **three** variables ("z-scored mean
PSI, PSI standard deviation and absolute disease effect size") and note that
"both schemes are reported and agree". `NUMBERS_125bp_REDERIVATION.md` section 5
is explicit: **"Keep the 3-D matching scheme ... Never mix with the 2-D
scheme"**, and gives the ranks to reproduce as MXE **2 of 11**, SE **5 of 11**.
This script calibrates on 3-D and prints 2-D beside it so the "both agree" claim
can be checked rather than assumed.

INTEGRITY CHECKS IT RUNS FIRST
------------------------------
1. `floor` is constant within each (class, model, matching) group -- if GAPDH's
   delta varied by candidate the column would not mean what its name says.
2. `margin == delta - floor` to 1e-6 everywhere.
3. exactly 10 candidate genes per group.

Any failure stops the run. A rank computed from a table that fails these is not
a rank.

USAGE
-----
  cd /data/project/chureeratp/RNA/HEART
  python3 152_gapdh_rank_recount.py
  python3 152_gapdh_rank_recount.py --schema-only
  python3 152_gapdh_rank_recount.py --self-test
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

NEED = ["event_type", "model", "matching", "gene", "delta", "floor", "margin",
        "ci_lo", "ci_hi", "holm_p"]

MODEL = os.environ.get("RBP_MODEL", "DCM/cell")
MATCH_3D = "3-D match (+|disease beta|)"
MATCH_2D = "2-D match (as published)"

# Published ranks, NUMBERS section 5, DCM/cell with 3-D matching.
CAL_RANK = {"MXE": int(os.environ.get("CAL_RANK_MXE", 2)),
            "SE": int(os.environ.get("CAL_RANK_SE", 5))}
CAL_N = int(os.environ.get("CAL_N_GENES", 11))

FILES = [("50 bp, published", "signature_rbp_margins.tsv"),
         ("50 bp, re-check", "signature_rbp_margins_50check.tsv"),
         ("125 bp", "signature_rbp_margins_rl125.tsv")]


# `delta`, `floor` and `margin` are stored rounded to five decimals, so
# margin - (delta - floor) can differ by up to 3 x 0.5e-5 = 1.5e-5 from
# rounding alone. My first version used 1e-6 and refused every real file at an
# observed deviation of exactly 1e-05 -- a tolerance error, not a data error.
# 2e-5 admits rounding and still catches a genuine inconsistency, and the
# observed maximum is always printed so drift cannot hide behind the threshold.
MARGIN_TOL = float(os.environ.get("MARGIN_TOL", 2e-5))


def integrity(T, label):
    """(problems, observed max |margin - (delta - floor)|)."""
    bad = []
    miss = [c for c in NEED if c not in T.columns]
    if miss:
        return [f"{label}: missing columns {miss}"], float("nan")
    worst = 0.0
    for (cl, mo, ma), g in T.groupby(["event_type", "model", "matching"]):
        f = g["floor"].dropna().unique()
        if len(f) != 1:
            bad.append(f"{label} {cl}/{mo}/{ma}: floor is not constant, "
                       f"{len(f)} distinct values {sorted(f)[:4]}")
        if len(g) != 10:
            bad.append(f"{label} {cl}/{mo}/{ma}: {len(g)} genes, expected 10")
        d = (g["margin"] - (g["delta"] - g["floor"])).abs().max()
        if pd.notna(d):
            worst = max(worst, float(d))
            if d > MARGIN_TOL:
                bad.append(f"{label} {cl}/{mo}/{ma}: margin != delta - floor "
                           f"(max deviation {d:.3g}, tolerance "
                           f"{MARGIN_TOL:.3g})")
    return bad, worst


def rank(g):
    """(floor rank, n, how many candidates beat it, who they are)."""
    f = float(g["floor"].iloc[0])
    above = g[g["delta"] > f].sort_values("delta", ascending=False)
    return 1 + len(above), len(g) + 1, len(above), list(above["gene"])


def cell(T, cl, match):
    g = T[(T.event_type == cl) & (T.model == MODEL) & (T.matching == match)]
    return g if len(g) else None


def show(label, cl, match, g):
    r, n, nabove, who = rank(g)
    f = float(g["floor"].iloc[0])
    top = g.sort_values("delta", ascending=False).iloc[0]
    print()
    print(f"  {label} — {cl}, {MODEL}, {match}")
    print(f"      floor (GAPDH) delta  {f:+.4f}")
    print(f"      FLOOR RANK           {r} of {n}"
          f"   ({nabove} candidate(s) above it"
          f"{': ' + ', '.join(who) if who else ''})")
    print(f"      leader               {top['gene']}"
          f"   delta {top['delta']:+.4f}"
          f"   margin {top['margin']:+.4f}"
          f"   [{top['ci_lo']:+.3f}, {top['ci_hi']:+.3f}]"
          f"   Holm {top['holm_p']:.3g}")
    cols = ["gene", "delta", "margin", "ci_lo", "ci_hi", "holm_p"]
    o = g.sort_values("delta", ascending=False)[cols].reset_index(drop=True)
    o.insert(0, "pos", [i + 1 if g["delta"].sort_values(
        ascending=False).iloc[i] > f else i + 2 for i in range(len(o))])
    with pd.option_context("display.width", 200):
        print(o.to_string(index=False))
    print(f"      (GAPDH floor would sit at position {r})")
    return r, n, nabove, who


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--schema-only", action="store_true")
    ap.add_argument("--self-test", action="store_true")
    a = ap.parse_args()
    if a.self_test:
        return self_test()

    self_test()
    print()
    print("=" * 78)
    print("152 — the GAPDH floor's rank, re-counted")
    print("=" * 78)
    print(f"  model    {MODEL}")
    print(f"  matching {MATCH_3D}   (2-D printed beside it for comparison)")
    print(f"  rank     1 + (delta > floor).sum(), out of 10 candidates + 1")
    print(f"  to reproduce at 50 bp: MXE {CAL_RANK['MXE']} of {CAL_N}, "
          f"SE {CAL_RANK['SE']} of {CAL_N}")

    loaded = []
    for label, path in FILES:
        if not (os.path.exists(path) and os.path.getsize(path)):
            print(f"\n  [skip] {label}: {path} not found")
            continue
        T = pd.read_csv(path, sep="\t", low_memory=False)
        loaded.append((label, path, T))
    if not loaded:
        print("\n  [fatal] no margins table found.", file=sys.stderr)
        return 2

    print()
    print("-" * 78)
    print("  INTEGRITY")
    print("-" * 78)
    allbad, worst = [], 0.0
    for label, path, T in loaded:
        bad, w = integrity(T, label)
        allbad += bad
        worst = max(worst, w if w == w else 0.0)
        print(f"    {label:<18} {path:<38} {len(T)} rows   "
              f"margin dev {w:.1e}   "
              f"{'OK' if not bad else str(len(bad)) + ' PROBLEM(S)'}")
    if allbad:
        print()
        print("  >>> REFUSING. A rank from a table that fails these is not a")
        print("      rank:")
        for b in allbad[:12]:
            print(f"        ✗ {b}")
        return 3
    print(f"    ✓ floor constant within every group, 10 candidates per")
    print(f"      group, and margin == delta - floor to {worst:.1e} "
          f"(tolerance {MARGIN_TOL:.1e};")
    print(f"      the columns are stored to 5 dp, so ~1.5e-05 is rounding)")

    results = {}
    for label, _p, T in loaded:
        print()
        print("-" * 78)
        print(f"  {label}")
        print("-" * 78)
        for cl in ("MXE", "SE"):
            for match in (MATCH_3D, MATCH_2D):
                g = cell(T, cl, match)
                if g is None:
                    print(f"\n  [warn] {label} {cl} {match}: no rows")
                    continue
                r, n, nabove, who = show(label, cl, match, g)
                results[(label, cl, match)] = (r, n, nabove, who)

    # ---- calibration on the 50 bp files, 3-D
    print()
    print("=" * 78)
    print("CALIBRATION")
    print("=" * 78)
    bad, checked = [], 0
    for (label, cl, match), (r, n, _na, _w) in results.items():
        if "50 bp" not in label or match != MATCH_3D:
            continue
        checked += 1
        if r != CAL_RANK[cl] or n != CAL_N:
            bad.append(f"{label} {cl}: floor rank {r} of {n}, published "
                       f"{CAL_RANK[cl]} of {CAL_N}")
    if not checked:
        print("  >>> NOT CALIBRATED — no 50 bp 3-D cell was computed.")
        return 3
    if bad:
        print("  >>> CALIBRATION FAILED:")
        for b in bad:
            print(f"      ✗ {b}")
        print()
        print("      Before assuming the files are wrong, check the matching")
        print("      scheme: NUMBERS section 5 records these ranks for 3-D,")
        print("      but the column labels the OTHER scheme 'as published'.")
        print("      Compare the 2-D rows printed above — if they reproduce")
        print("      2 and 5 instead, the published values are 2-D and")
        print("      NUMBERS section 5 needs correcting, not this script.")
        return 3
    print(f"  ✓ the published 50 bp 3-D ranks reproduce (MXE "
          f"{CAL_RANK['MXE']} of {CAL_N}, SE {CAL_RANK['SE']} of {CAL_N})")

    # ---- the three claims, at 125 bp
    print()
    print("=" * 78)
    print("THE THREE CLAIMS, AT 125 bp")
    print("=" * 78)
    for cl, claim in (("MXE", "[081]/[096] 'outscored nine of the ten "
                              "candidates, RBM20 the single exception'"),
                      ("SE", "[081] 'the floor sits mid-table, fifth of "
                             "the eleven'")):
        k = ("125 bp", cl, MATCH_3D)
        if k not in results:
            print(f"\n  {cl}: not computed")
            continue
        r, n, nabove, who = results[k]
        print()
        print(f"  {cl}  {claim}")
        print(f"      measured: floor rank {r} of {n}, "
              f"{nabove} above it{': ' + ', '.join(who) if who else ''}")
        if cl == "MXE":
            ok = (r == 2 and who == ["RBM20"])
            print(f"      the published clause is "
                  f"{'STILL TRUE' if ok else 'NO LONGER TRUE'}")
            if not ok:
                print(f"      -> rewrite as: the floor outscored "
                      f"{10 - nabove} of the ten candidates, with "
                      f"{', '.join(who) if who else 'none'} above it")
        else:
            print(f"      published position 5; measured {r}")
            if r != 5:
                print(f"      -> 'fifth of the eleven' becomes "
                      f"'{r}' of the eleven")
    # QKI, for the abstract
    k = ("125 bp", "SE", MATCH_3D)
    if k in results:
        g = cell(loaded[-1][2], "SE", MATCH_3D)
        o = g.sort_values("delta", ascending=False).reset_index(drop=True)
        qki = o.index[o.gene == "QKI"]
        if len(qki):
            i = int(qki[0])
            print()
            print(f"  [007] 'QKI included' — QKI is candidate "
                  f"{i + 1} of 10 by SE delta, margin "
                  f"{o.loc[i, 'margin']:+.4f}")
            print("        Naming it stays accurate (it is in the panel), but")
            print("        if it is no longer the leader do not imply it is")
            print("        the best case for a regulator.")
    print()
    print("  ⚠ Report the rank and the margin WITH its interval. The floor is")
    print("    a single-gene estimate, so a bare delta beside 58's bootstrap")
    print("    intervals is the inconsistency 58 itself warns about.")
    return 0


def self_test():
    import tempfile

    genes = ["RBM20", "PTBP1", "QKI", "CELF1", "MBNL1", "HNRNPF", "SRSF1",
             "RBM5", "ZRANB2", "HNRNPA0"]
    # MXE: floor 0.80, only RBM20 above -> rank 2 of 11
    mxe = dict(zip(genes, [0.90, 0.70, 0.60, 0.50, 0.40,
                           0.30, 0.20, 0.10, 0.05, 0.01]))
    # SE: floor 0.50, four above -> rank 5 of 11
    se = dict(zip(genes, [0.70, 0.80, 0.90, 0.60, 0.40,
                          0.30, 0.20, 0.10, 0.05, 0.01]))
    rows = []
    for cl, dmap, fl in (("MXE", mxe, 0.80), ("SE", se, 0.50)):
        for mo in ("DCM/cell", "ALL/nocell"):
            for ma in (MATCH_3D, MATCH_2D):
                for gname in genes:
                    d = dmap[gname]
                    # under the wrong model the floor is below everything
                    f = fl if mo == "DCM/cell" else -1.0
                    rows.append(dict(event_type=cl, model=mo, matching=ma,
                                     gene=gname, delta=d, floor=f,
                                     margin=d - f, ci_lo=d - f - 0.1,
                                     ci_hi=d - f + 0.1, boot_p=0.3,
                                     holm_p=0.9, n_target=100, n_bg=1000))
    T = pd.DataFrame(rows)
    assert len(T) == 2 * 2 * 2 * 10 == 80, len(T)

    assert not integrity(T, "fx")[0], integrity(T, "fx")

    g = cell(T, "MXE", MATCH_3D)
    r, n, nabove, who = rank(g)
    assert (r, n, nabove, who) == (2, 11, 1, ["RBM20"]), (r, n, nabove, who)

    g = cell(T, "SE", MATCH_3D)
    r, n, nabove, who = rank(g)
    assert (r, n) == (5, 11), (r, n)
    assert who == ["QKI", "PTBP1", "RBM20", "CELF1"], who

    # the wrong model must NOT reproduce the published rank
    gw = T[(T.event_type == "MXE") & (T.model == "ALL/nocell") &
           (T.matching == MATCH_3D)]
    rw, nw, _na, _w = rank(gw)
    assert rw == 11 and nw == 11, (rw, nw)   # floor last when it is lowest

    # --- integrity catches a non-constant floor
    T2 = T.copy()
    T2.loc[(T2.event_type == "MXE") & (T2.model == "DCM/cell") &
           (T2.matching == MATCH_3D) & (T2.gene == "QKI"), "floor"] = 0.99
    b2, _w2 = integrity(T2, "fx")
    assert any("floor is not constant" in x for x in b2), b2

    # --- and a margin that is not delta - floor
    T3 = T.copy()
    T3.loc[(T3.event_type == "SE") & (T3.gene == "MBNL1"), "margin"] = 42.0
    b3, _w3 = integrity(T3, "fx")
    assert any("margin != delta - floor" in x for x in b3), b3

    # --- rounding-sized noise must PASS. This is the case my first
    #     version refused on every real file.
    T_round = T.copy()
    T_round["margin"] = (T_round["margin"] + 1e-5).round(5)
    br, wr = integrity(T_round, "fx")
    assert not br, (br, wr)
    assert wr <= MARGIN_TOL, (wr, MARGIN_TOL)

    # --- and a group with the wrong number of genes
    T4 = T[T.gene != "SRSF1"]
    b4, _w4 = integrity(T4, "fx")
    assert any("9 genes, expected 10" in x for x in b4), b4

    # --- a missing column is reported, not crashed on
    b5, _w5 = integrity(T.drop(columns=["floor"]), "fx")
    assert b5 and "missing columns" in b5[0], b5

    # --- ties: a candidate exactly equal to the floor is NOT above it
    T6 = T.copy()
    m = ((T6.event_type == "MXE") & (T6.model == "DCM/cell") &
         (T6.matching == MATCH_3D) & (T6.gene == "PTBP1"))
    T6.loc[m, "delta"] = 0.80
    T6.loc[m, "margin"] = 0.0
    r6, n6, na6, _w = rank(cell(T6, "MXE", MATCH_3D))
    assert (r6, na6) == (2, 1), (r6, na6)

    d = tempfile.mkdtemp()
    T.to_csv(os.path.join(d, "signature_rbp_margins.tsv"),
             sep="\t", index=False)

    print("[ok] self-test passed: the floor is read from the `floor` COLUMN "
          "(not a GAPDH row, which is what the real files carry) and ranked "
          "as 1 + (delta > floor), giving 2 of 11 in MXE with RBM20 alone "
          "above it and 5 of 11 in SE; a candidate exactly equal to the floor "
          "counts as not above it; the wrong model does not reproduce the "
          "published rank; and the integrity checks catch a non-constant "
          "floor, a margin that is not delta minus floor (while 1e-05 "
          "rounding noise PASSES, the case the first version wrongly "
          "refused), a group with 9 "
          "genes, and a missing column")
    return 0


if __name__ == "__main__":
    sys.exit(main())
