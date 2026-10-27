#!/usr/bin/env python3
"""
121_selection_vs_precision.py -- why did every transfer slope rise at 125 bp?

THE QUESTION
------------
Re-quantifying MAGNet at 125 bp raised every Deming transfer slope:

    Colorado 0.521 -> 0.600   LBNL 0.445 -> 0.573   Brigham 0.476 -> 0.510
    Fuwai    0.243 -> 0.295   Stanford 0.386 -> 0.401

and A in the matched panel rose 0.584 -> 0.690. Naively a genuine read-length
difference should make the platforms MORE different and the slopes LOWER. They
went up, in all five cohorts.

Two candidate causes, and `RESULT_matched_panel_rl125.md` forbids asserting
either without measurement:

  SELECTION   the event set changed. MAGNet-significant SE fell 926 -> 813, and
              this project has a documented selection gradient of 0.867 -> 0.440,
              large enough to produce the whole effect on its own.

  PRECISION   MAGNet's PSI standard errors changed, moving the Deming error
              ratio lambda and hence the slope, with no change in biology.

That all five cohorts moved together does NOT discriminate: the selection set
shrank for all five simultaneously, so selection looks exactly like this.

THE DESIGN — a 2 x 2, five times over
-------------------------------------
The two factors are separable because the event set is just a list of ekeys:

                        |  event set = 50 bp (926)  |  event set = 125 bp (813)
    --------------------+---------------------------+--------------------------
    x from 50 bp table  |  cell 1 = the PUBLISHED   |  cell 4
                        |           50 bp slope     |
    x from 125 bp table |  cell 3                   |  cell 2 = the NEW
                        |                           |           125 bp slope

Cells 1 and 2 are the known answers and are used as CALIBRATION: if they do not
reproduce, the estimator or the selection rule has been replicated wrongly and
nothing is reported.

Cell 3 is the discriminating measurement. Holding the event set at the 50 bp
926 and changing only the quantification:

  * cell 3 close to cell 1  ->  changing the data alone does nothing.
                                The rise is SELECTION.
  * cell 3 close to cell 2  ->  changing the data alone reproduces the rise.
                                The rise is PRECISION.
  * in between              ->  both, and the script reports the shares.

Cell 4 gives the second path through the square, so the two decompositions can
be compared. If they disagree, the factors INTERACT and neither "share" is
meaningful on its own -- which is itself a reportable finding, not a failure.

ONE SOURCE OF TRUTH
-------------------
`load()`, `deming()`, `REPL` and `NKEEP` are imported from
`33_rebuild_summaries_dedup.py` rather than reimplemented, so the collapse keys,
the dedup preference and the errors-in-variables estimator are 33's own. Only
the two tiny helpers that live inside 33's `main()` (the dPSI column choice and
the significance filter) are replicated -- and the ten calibration values are
what prove the replication is faithful.

lambda is reported per cell as an independent diagnostic. It is computed here,
not taken from 33, so it corroborates rather than confirms.

USAGE
-----
  python 121_selection_vs_precision.py
  python 121_selection_vs_precision.py --v33 33_rebuild_summaries_dedup.py
  python 121_selection_vs_precision.py --self-test
"""
import argparse
import importlib.util
import os
import sys

import numpy as np
import pandas as pd

D50 = "SE_v2filtered"
D125 = "SE_v2filtered_rl125"

# calibration — section E of summaries_dedup.txt / summaries_dedup_rl125.txt
EXP_SLOPE_50 = {"Colorado": 0.521, "LBNL": 0.445, "Brigham": 0.476,
                "Fuwai": 0.243, "Stanford": 0.386}
EXP_SLOPE_125 = {"Colorado": 0.600, "LBNL": 0.573, "Brigham": 0.510,
                 "Fuwai": 0.295, "Stanford": 0.401}
EXP_N_50 = {"Colorado": 471, "LBNL": 199, "Brigham": 442,
            "Fuwai": 714, "Stanford": 536}
EXP_N_125 = {"Colorado": 385, "LBNL": 156, "Brigham": 364,
             "Fuwai": 596, "Stanford": 455}
EXP_SET_50, EXP_SET_125 = 926, 813

# Each cell carries its own data table, event set and label in ONE place, so
# the three cannot drift apart. The first version of this script unpacked the
# four cells positionally and printed them under headers in a different order,
# which mislabelled every column of the 2x2 while the decomposition beneath it
# was correct -- the worst kind of display bug, because the arithmetic looked
# right and the table did not.
CELLS = {
    "c1": dict(data="50",  eventset="50",  label="1: 50d/50s",
               note="(published)"),
    "c3": dict(data="125", eventset="50",  label="3: 125d/50s",
               note="DATA moved"),
    "c4": dict(data="50",  eventset="125", label="4: 50d/125s",
               note="SET moved"),
    "c2": dict(data="125", eventset="125", label="2: 125d/125s",
               note="(new)"),
}
ORDER = ["c1", "c3", "c4", "c2"]        # the path order across the square


def load_33(path):
    spec = importlib.util.spec_from_file_location("v33", path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)          # main() is behind __main__, safe
    for need in ("load", "deming", "REPL", "NKEEP"):
        if not hasattr(m, need):
            raise SystemExit(f"[fatal] {path} has no {need}()")
    return m


# --- the two helpers that live inside 33's main() and cannot be imported
def dcol(t):
    return next((c for c in ("dPSI_adj", "dPSI") if c in t.columns), None)


def fdrcol(t):
    return next((c for c in ("FDR", "FDR_adj", "fdr") if c in t.columns), None)


def significant(t, fdr_cut=0.05, dpsi_cut=0.10):
    """FDR<0.05 AND |dPSI|>=0.10, on an already-collapsed table."""
    d, q = dcol(t), fdrcol(t)
    if d is None or q is None:
        return None
    return ((pd.to_numeric(t[q], errors="coerce") < fdr_cut)
            & (pd.to_numeric(t[d], errors="coerce").abs() >= dpsi_cut))


def lam(sx, sy):
    """mean(sy^2)/mean(sx^2) -- the error-variance ratio, as a diagnostic."""
    a, b = np.nanmean(np.asarray(sy, float) ** 2), \
        np.nanmean(np.asarray(sx, float) ** 2)
    return float(a / b) if b and np.isfinite(a / b) else float("nan")


def cell(v33, disc, ekeys, coh_tbl):
    """Deming slope of coh_tbl's dPSI on disc's dPSI, over `ekeys`.

    -> (slope, n, lambda) or (None, 0, nan)
    """
    dd, d2 = dcol(disc), dcol(coh_tbl)
    if dd is None or d2 is None:
        return None, 0, float("nan")
    sx = "se_adj" if "se_adj" in disc.columns else None
    s2 = "se_adj" if "se_adj" in coh_tbl.columns else None
    if sx is None or s2 is None:
        # 33 falls back to OLS here; refuse rather than silently switch
        return None, 0, float("nan")
    X = disc.loc[disc.ekey.isin(ekeys), ["ekey", dd, sx]].rename(
        columns={dd: "x", sx: "sx"})
    Y = coh_tbl[["ekey", d2, s2]].rename(columns={d2: "y", s2: "sy"})
    m = X.merge(Y, on="ekey")
    m = m.dropna(subset=["x", "y", "sx", "sy"])
    if len(m) < 10:
        return None, len(m), float("nan")
    sl, n = v33.deming(m.x.values, m.y.values, m.sx.values, m.sy.values)
    return float(sl), int(n), lam(m.sx.values, m.sy.values)


def verdict(c1, c2, c3, c4):
    """-> (label, share_precision, share_selection, interaction)"""
    total = c2 - c1
    if abs(total) < 1e-9:
        return "no change to explain", 0.0, 0.0, 0.0
    # path 1: change the data first (set held at 50 bp), then the set
    data_first = c3 - c1
    set_second = c2 - c3
    # path 2: change the set first (data held at 50 bp), then the data
    set_first = c4 - c1
    data_second = c2 - c4
    # The 2x2 interaction term. NOT a difference of the path increments --
    # that expression cancels identically and would report zero interaction
    # for every input, which is how the first version of this function passed
    # a strong-interaction fixture as "BOTH".
    interaction = c2 - c3 - c4 + c1
    sp = 0.5 * (data_first + data_second) / total
    ss = 0.5 * (set_first + set_second) / total
    if abs(interaction) > 0.5 * abs(total):
        label = "INTERACTION dominates — shares are not separable"
    elif abs(sp) >= 2 * abs(ss):
        label = "PRECISION"
    elif abs(ss) >= 2 * abs(sp):
        label = "SELECTION"
    else:
        label = "BOTH"
    return label, sp, ss, interaction


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--v33", default="33_rebuild_summaries_dedup.py")
    ap.add_argument("--d50", default=D50)
    ap.add_argument("--d125", default=D125)
    ap.add_argument("--self-test", action="store_true")
    a = ap.parse_args()
    if a.self_test:
        return self_test()

    self_test()
    print()
    print("=" * 84)
    print("121 — selection or precision? The 2x2, five cohorts")
    print("=" * 84)
    if not os.path.exists(a.v33):
        print(f"[fatal] {a.v33} not found", file=sys.stderr)
        return 2
    v33 = load_33(a.v33)
    print(f"  estimator, collapse keys and cohort table imported from {a.v33}")
    print(f"  NKEEP = {v33.NKEEP}")
    print()

    d50, p50 = v33.load(a.d50, "SE", prefer_dedup=False)
    d125, p125 = v33.load(a.d125, "SE", prefer_dedup=False)
    for df, p in ((d50, p50), (d125, p125)):
        if df is None:
            print(f"[fatal] discovery table not found: {p}", file=sys.stderr)
            return 2
    print(f"  50 bp : {p50}   {len(d50):,} collapsed events")
    print(f"  125 bp: {p125}   {len(d125):,} collapsed events")

    s50, s125 = significant(d50, ), significant(d125)
    if s50 is None or s125 is None:
        print("[fatal] no FDR or dPSI column in a discovery table",
              file=sys.stderr)
        return 2
    set50 = set(d50.loc[s50, "ekey"].astype(str))
    set125 = set(d125.loc[s125, "ekey"].astype(str))
    print(f"  significant (FDR<0.05 & |dPSI|>=0.10):  "
          f"50 bp {len(set50):,}   125 bp {len(set125):,}")
    if len(set50) != EXP_SET_50 or len(set125) != EXP_SET_125:
        print(f"  ✗ CALIBRATION FAIL — expected {EXP_SET_50} and "
              f"{EXP_SET_125}.", file=sys.stderr)
        print("    The significance filter replicated from 33's main() does "
              "not match 33.\n    Nothing is reported.", file=sys.stderr)
        return 3
    print(f"  ✓ set sizes reproduce ({EXP_SET_50} / {EXP_SET_125})")
    shared = len(set50 & set125)
    print(f"    overlap {shared:,}  "
          f"({100 * shared / len(set125):.1f}% of the 125 bp set)")
    print()

    rows, fails = [], []
    src = {"50": (d50, set50), "125": (d125, set125)}
    for coh, pat, _flag, _why in v33.REPL:
        t, tp = v33.load(pat.format(ET="SE"), "SE")
        if t is None:
            print(f"  {coh:<10} cohort table not found: {tp}")
            continue
        # every cell computed FROM ITS SPEC -- no positional unpacking
        got = {}
        for key in ORDER:
            spec = CELLS[key]
            got[key] = cell(v33, src[spec["data"]][0],
                            src[spec["eventset"]][1], t)
        if any(got[k][0] is None for k in ORDER):
            print(f"  {coh:<10} a cell could not be fitted "
                  f"(se_adj missing, or fewer than 10 matched events)")
            continue
        c1, n1 = got["c1"][0], got["c1"][1]
        c2, n2 = got["c2"][0], got["c2"][1]
        # --- calibration on the two known cells
        if coh in EXP_SLOPE_50 and round(c1, 3) != EXP_SLOPE_50[coh]:
            fails.append(f"{coh} cell 1 = {c1:.3f}, recorded "
                         f"{EXP_SLOPE_50[coh]:.3f}")
        if coh in EXP_SLOPE_125 and round(c2, 3) != EXP_SLOPE_125[coh]:
            fails.append(f"{coh} cell 2 = {c2:.3f}, recorded "
                         f"{EXP_SLOPE_125[coh]:.3f}")
        if coh in EXP_N_50 and n1 != EXP_N_50[coh]:
            fails.append(f"{coh} cell 1 n = {n1}, recorded {EXP_N_50[coh]}")
        if coh in EXP_N_125 and n2 != EXP_N_125[coh]:
            fails.append(f"{coh} cell 2 n = {n2}, recorded {EXP_N_125[coh]}")
        rows.append((coh, got))

    if fails:
        print("=" * 84)
        print("✗ CALIBRATION FAILED — nothing is reported.")
        for f in fails:
            print(f"    {f}")
        print()
        print("  Cells 1 and 2 must reproduce section E of summaries_dedup.txt")
        print("  and summaries_dedup_rl125.txt exactly. They do not, so the")
        print("  significance filter or the merge replicated here differs from")
        print("  33's, and cells 3 and 4 cannot be trusted either.")
        return 3
    if not rows:
        print("[fatal] no cohort produced four fitted cells", file=sys.stderr)
        return 2

    print("=" * 84)
    print("THE 2x2 — Deming slope (n) [lambda]")
    print("=" * 84)
    print(f"  {'cohort':<10}"
          + "".join(f"{CELLS[k]['label']:>18}" for k in ORDER))
    print(f"  {'':<10}"
          + "".join(f"{CELLS[k]['note']:>18}" for k in ORDER))
    for coh, got in rows:
        print(f"  {coh:<10}" + "".join(
            f"{f'{got[k][0]:.3f} ({got[k][1]})':>18}" for k in ORDER))
    print()
    print(f"  {'lambda':<10}" + "".join(
        f"{np.nanmean([r[1][k][2] for r in rows]):>18.3f}" for k in ORDER)
        + "   (mean across cohorts)")

    print()
    print("=" * 84)
    print("DECOMPOSITION")
    print("=" * 84)
    print(f"  {'cohort':<10}{'total':>9}{'precision':>12}{'selection':>12}"
          f"{'inter':>9}   verdict")
    labels = []
    for coh, got in rows:
        c1, c2, c3, c4 = (got["c1"][0], got["c2"][0],
                          got["c3"][0], got["c4"][0])
        lab, sp, ss, inter = verdict(c1, c2, c3, c4)
        labels.append(lab)
        print(f"  {coh:<10}{c2 - c1:>+9.3f}{100 * sp:>11.0f}%"
              f"{100 * ss:>11.0f}%{inter:>+9.3f}   {lab}")
    print()
    print("  precision% and selection% are shares of the total change,")
    print("  averaged over the two paths through the square. They can exceed")
    print("  100% or go negative when the two factors push in opposite")
    print("  directions — that is information, not an error.")

    print()
    print("=" * 84)
    print("READING")
    print("=" * 84)
    from collections import Counter
    tally = Counter(labels)
    print(f"  verdicts: {dict(tally)}")
    top, cnt = tally.most_common(1)[0]
    if cnt == len(labels):
        print(f"  All {len(labels)} cohorts agree: **{top}**.")
        print("  Five independent cohorts reaching the same verdict is strong;")
        print("  they share the MAGNet x-axis, so they are not fully")
        print("  independent, but their y-axes are separate datasets.")
    else:
        print(f"  NOT unanimous — {cnt} of {len(labels)} say {top}.")
        print("  Report the mixture, not the majority. A cohort that disagrees")
        print("  is evidence that the two channels differ in size per cohort,")
        print("  which the manuscript cannot summarise in one word.")
    print()
    print("  ⚠ Whatever the verdict, it describes WHY THE SLOPE MOVED. It does")
    print("    not make either slope a magnitude estimate — both arms remain")
    print("    attenuated by regression to the mean, and the selection")
    print("    gradient (0.867 -> 0.440) is itself the reason a selection")
    print("    verdict would be unsurprising rather than reassuring.")
    print()
    print("  ⚠ If the verdict is SELECTION, the 125 bp slopes are not 'better")
    print("    measurements of the same thing' — they are measurements of a")
    print("    smaller, more stringently selected set. Say so where the")
    print("    slopes are quoted.")
    return 0


def self_test():
    # --- every cell's label must describe its own data/set combination.
    # This is what the first version got wrong: the labels were a separate
    # literal from the computation order.
    assert set(ORDER) == set(CELLS) == {"c1", "c2", "c3", "c4"}, ORDER
    for k, sp in CELLS.items():
        assert f"{sp['data']}d/{sp['eventset']}s" in sp["label"], (k, sp)
    assert CELLS["c1"]["data"] == CELLS["c1"]["eventset"] == "50"
    assert CELLS["c2"]["data"] == CELLS["c2"]["eventset"] == "125"
    # the two crosses must differ in exactly one factor from c1
    for k in ("c3", "c4"):
        diff = sum(CELLS[k][f] != CELLS["c1"][f]
                   for f in ("data", "eventset"))
        assert diff == 1, (k, CELLS[k])

    # --- decomposition arithmetic, with the answer known by construction
    # pure precision: moving the data does all the work
    lab, sp, ss, it = verdict(0.521, 0.600, 0.600, 0.521)
    assert lab == "PRECISION", (lab, sp, ss, it)
    assert abs(sp - 1.0) < 1e-9 and abs(ss) < 1e-9, (sp, ss)
    # pure selection: moving the set does all the work
    lab, sp, ss, it = verdict(0.521, 0.600, 0.521, 0.600)
    assert lab == "SELECTION", (lab, sp, ss, it)
    assert abs(ss - 1.0) < 1e-9 and abs(sp) < 1e-9, (sp, ss)
    # half and half
    mid = 0.5 * (0.521 + 0.600)
    lab, sp, ss, it = verdict(0.521, 0.600, mid, mid)
    assert lab == "BOTH", (lab, sp, ss)
    assert abs(sp - 0.5) < 1e-9 and abs(ss - 0.5) < 1e-9, (sp, ss)
    # strong interaction: both single-factor cells overshoot the endpoint
    lab, _sp, _ss, it = verdict(0.5, 0.6, 0.9, 0.9)
    assert "INTERACTION" in lab, (lab, it)
    # nothing to explain
    assert verdict(0.5, 0.5, 0.5, 0.5)[0] == "no change to explain"

    # --- the restriction and merge, on a frame where the slope is known
    class V:
        @staticmethod
        def deming(x, y, sx, sy):
            # OLS through the origin is enough to check plumbing
            x = np.asarray(x, float)
            y = np.asarray(y, float)
            return float((x * y).sum() / (x * x).sum()), len(x)

    disc = pd.DataFrame({"ekey": [f"e{i}" for i in range(20)],
                         "dPSI": np.linspace(0.1, 0.5, 20),
                         "se_adj": 0.01, "FDR": 0.001})
    coh = pd.DataFrame({"ekey": [f"e{i}" for i in range(20)],
                        "dPSI": np.linspace(0.05, 0.25, 20),
                        "se_adj": 0.02})
    sl, n, L = cell(V, disc, {f"e{i}" for i in range(20)}, coh)
    assert n == 20, n
    assert abs(sl - 0.5) < 1e-6, sl            # y = x/2 exactly
    assert abs(L - 4.0) < 1e-9, L              # (0.02^2)/(0.01^2)
    # restricting the ekeys must restrict n
    sl2, n2, _ = cell(V, disc, {f"e{i}" for i in range(12)}, coh)
    assert n2 == 12, n2
    # fewer than ten matched events -> no fit
    sl3, n3, _ = cell(V, disc, {"e0", "e1"}, coh)
    assert sl3 is None and n3 == 2, (sl3, n3)
    # se_adj missing -> refuse rather than fall back to OLS
    sl4, _n4, _ = cell(V, disc.drop(columns=["se_adj"]),
                       {f"e{i}" for i in range(20)}, coh)
    assert sl4 is None, "missing se_adj should not silently fit"

    # --- the significance filter
    t = pd.DataFrame({"ekey": list("abcd"),
                      "dPSI": [0.2, 0.05, 0.3, 0.4],
                      "FDR": [0.01, 0.01, 0.9, 0.04]})
    s = significant(t)
    assert list(t.loc[s, "ekey"]) == ["a", "d"], list(t.loc[s, "ekey"])

    print("[ok] self-test passed: pure-precision, pure-selection, 50/50 and "
          "strong-interaction cases each get the right verdict and shares; the "
          "ekey restriction, the 10-event floor, the lambda diagnostic and the "
          "FDR-and-dPSI filter all behave; a missing se_adj refuses instead of "
          "silently switching estimator")
    return 0


if __name__ == "__main__":
    sys.exit(main())
