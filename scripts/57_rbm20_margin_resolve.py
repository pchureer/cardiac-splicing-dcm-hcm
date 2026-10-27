#!/usr/bin/env python3
"""
57_rbm20_margin_resolve.py -- settle the RBM20-vs-GAPDH margin properly.

WHY THIS EXISTS. 56 compared RBM20's enrichment delta with GAPDH's by parsing the
printed report, which rounds medians to two decimals. The margins it returned were
+0.010 (SE_techadj), +0.000 (MXE_techadj), +0.070 (MXE) and -0.020 (SE) -- i.e. at or
below the rounding error of the numbers it was reading, and incoherent across runs
(the null held unadjusted and failed adjusted in SE, which is not a pattern biology
produces). A decision rule applied at finer resolution than the data supports is not
a decision rule. This script replaces it.

WHAT IT ADDS, in order of how much it matters:

  1. FULL PRECISION. Deltas recomputed from the events tables, not from the report.
  2. UNCERTAINTY. A paired bootstrap CI on (RBM20 delta - GAPDH delta). Paired because
     both genes are scored on the SAME target and background events, so their t-values
     are correlated across events; an unpaired bootstrap would overstate the spread.
  3. RANK. Where RBM20 sits among all fitted genes including GAPDH. This is the robust
     framing: "no candidate separates from the housekeeping control" survives any
     rounding argument if RBM20 is mid-pack, whatever the absolute deltas are.
  4. DECOMPOSITION. How much of the 9 -> 45 delta change came from the target median
     moving and how much from the matched background moving. If the target statistic
     is unchanged and only the comparator moved, the "enrichment" is a property of the
     background construction, not of the gene set.
  5. THE THREE-DIMENSIONAL MATCH. `RESULT_rbm20_mechanism_null.md` records that the
     background is matched on mean PSI and PSI sd but NOT on disease effect size, that
     this biases the test toward enrichment, and that the fix was judged unnecessary
     "because it can only strengthen a null that already holds". That reasoning expired
     the moment the null stopped clearly holding. This adds |disease beta| as a third
     matching dimension and reports the deltas under it.

READ IT LIKE THIS. The question is not "is RBM20's delta positive" -- at 45 genes
every gene's delta is positive, which is the inflation the project already documented.
The question is "is RBM20 distinguishable from a housekeeping gene". Two things answer
it: the margin CI covering zero, and RBM20's rank.

Needs: pandas, numpy, scipy. Runs in seconds on the existing tables; no sbatch needed.

  python 57_rbm20_margin_resolve.py
  python 57_rbm20_margin_resolve.py --n-boot 5000 --out margin_resolve.tsv
"""

import argparse
import os
import re
import sys

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

NINE = "TTN,CAMK2D,LDB3,LMO7,SORBS1,PDLIM3,DTNA,MLIP,TNNT2".split(",")
G45 = ("APTX,CACNA1C,CAMK2D,CAMK2G,DAB1,DNM3,DST,DTNA,ENAH,FHOD3,FNBP1,GIT2,IMMT,"
       "KALRN,KCNIP2,LDB3,LMO7,LRRFIP1,MECP2,MLIP,MTMR1,MYH7,MYOM1,NEXN,NFIA,NPRL3,"
       "NTRK3,OBSCN,PDLIM3,PDLIM5,PLEKHA5,RALGPS1,RTN4,RYR2,SEMA6D,SH3KBP1,SLC38A10,"
       "SORBS1,SPEN,TNNT2,TPM1,TRDN,TTN,UBE2F,ZNF451").split(",")

# the models script 03 fits; DCM/cell is the one the manuscript reports
MODELS = [("DCM", "cell"), ("DCM", "nocell"), ("ALL", "cell"), ("ALL", "nocell")]
REPORTED = ("DCM", "cell")
N_MATCH = 10           # same as script 03's default; changing it changes the background


def match_background(res, is_tgt, extra_feature=None):
    """Replicate script 03's matching exactly, optionally with a third dimension.

    Script 03: z-score [mean_psi, sd_psi], cKDTree over background rows, k=10 nearest
    per target, union. Reproduced here rather than imported so this script can be run
    against an events table without re-importing the fitting code.
    """
    feats = [res["mean_psi"].values, res["sd_psi"].values]
    if extra_feature is not None:
        feats.append(extra_feature)
    f = np.column_stack(feats).astype(float)
    sd = np.nanstd(f, axis=0)
    sd[sd == 0] = 1.0
    fz = (f - np.nanmean(f, axis=0)) / sd
    ok = np.isfinite(fz).all(axis=1)
    bg_idx = np.where(~is_tgt & ok)[0]
    tg_idx = np.where(is_tgt & ok)[0]
    if not len(bg_idx) or not len(tg_idx):
        return tg_idx, np.array([], int)
    tree = cKDTree(fz[bg_idx])
    k = min(N_MATCH, len(bg_idx))
    _, nb = tree.query(fz[tg_idx], k=k)
    return tg_idx, np.unique(bg_idx[np.atleast_2d(nb).ravel()])


def delta(tv, tg_idx, bg_idx):
    """median |t| on targets minus median |t| on matched background, full precision."""
    a = tv[tg_idx]; b = tv[bg_idx]
    a = a[np.isfinite(a)]; b = b[np.isfinite(b)]
    if len(a) < 3 or len(b) < 3:
        return np.nan, np.nan, np.nan, 0, 0
    mt, mb = float(np.median(a)), float(np.median(b))
    return mt - mb, mt, mb, len(a), len(b)


def boot_margin(tv_a, tv_b, tg_idx, bg_idx, n_boot, rng):
    """Paired bootstrap CI for delta(gene A) - delta(gene B).

    One resample of target rows and one of background rows per iteration, applied to
    BOTH genes. Pairing matters: RBM20's and GAPDH's t-values are correlated across
    events (they share the library-level technical variation that inflates the whole
    test), so resampling them independently would inflate the interval and make every
    comparison look inconclusive.
    """
    ok_t = np.isfinite(tv_a[tg_idx]) & np.isfinite(tv_b[tg_idx])
    ok_b = np.isfinite(tv_a[bg_idx]) & np.isfinite(tv_b[bg_idx])
    at, bt = tv_a[tg_idx][ok_t], tv_b[tg_idx][ok_t]
    ab, bb = tv_a[bg_idx][ok_b], tv_b[bg_idx][ok_b]
    if len(at) < 5 or len(ab) < 5:
        return np.nan, np.nan, np.nan
    nt, nb = len(at), len(ab)
    out = np.empty(n_boot)
    for i in range(n_boot):
        it = rng.integers(0, nt, nt)
        ib = rng.integers(0, nb, nb)
        out[i] = ((np.median(at[it]) - np.median(ab[ib])) -
                  (np.median(bt[it]) - np.median(bb[ib])))
    lo, hi = np.percentile(out, [2.5, 97.5])
    # two-sided bootstrap p for "margin != 0"
    p = 2 * min((out <= 0).mean(), (out >= 0).mean())
    return float(lo), float(hi), float(min(p, 1.0))


def genes_fitted(res, model, tag):
    """Which genes have a fitted t column in this table."""
    pat = re.compile(rf"^{model}_([A-Z0-9]+)_{tag}\.\1\.t$")
    return sorted({m.group(1) for c in res.columns for m in [pat.match(c)] if m})


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tables", nargs="*", default=None,
                    help="events tables; default rbm20_{SE,MXE}[_techadj].events.tsv")
    ap.add_argument("--n-boot", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", default="rbm20_margin_resolve.tsv")
    a = ap.parse_args()

    tabs = a.tables or [f"rbm20_{et}{s}.events.tsv"
                        for et in ("SE", "MXE") for s in ("", "_techadj")]
    tabs = [t for t in tabs if os.path.exists(t)]
    if not tabs:
        sys.exit("[fatal] no events tables found. Expected "
                 "rbm20_{SE,MXE}[_techadj].events.tsv in the working directory.")
    rng = np.random.default_rng(a.seed)
    rows = []

    for tab in tabs:
        res = pd.read_csv(tab, sep="\t", low_memory=False)
        need = {"gene", "mean_psi", "sd_psi"}
        if not need <= set(res.columns):
            print(f"[skip] {tab}: missing {sorted(need - set(res.columns))}")
            continue
        tag_name = os.path.basename(tab).replace("rbm20_", "").replace(".events.tsv", "")
        print("\n" + "=" * 78)
        print(f"{tag_name}   ({len(res)} events)")
        print("=" * 78)

        present = set(res["gene"].astype(str))
        g45 = [g for g in G45 if g in present]
        g9 = [g for g in NINE if g in present]
        print(f"  target genes with testable events: {len(g9)}/9 and {len(g45)}/45")

        for model, ctag in MODELS:
            fitted = genes_fitted(res, model, ctag)
            if "RBM20" not in fitted:
                continue
            star = "  <-- the model the manuscript reports" if (model, ctag) == REPORTED else ""
            print(f"\n  --- {model}/{ctag} ---{star}")
            if "GAPDH" not in fitted:
                print("      GAPDH not fitted in this table -- no calibration floor,"
                      " cannot decide")
            tv = {g: np.abs(res[f"{model}_{g}_{ctag}.{g}.t"].values) for g in fitted}

            # --- 2-D match (as published) at 9 and at 45 --------------------
            store = {}
            for lab, gl in (("9", g9), ("45", g45)):
                is_t = res["gene"].isin(gl).values
                tg, bg = match_background(res, is_t)
                per = {}
                for g in fitted:
                    d, mt, mb, nt, nb = delta(tv[g], tg, bg)
                    per[g] = (d, mt, mb, nt, nb)
                store[lab] = (tg, bg, per)
                print(f"      [{lab:>2} genes] {len(tg)} target events, "
                      f"{len(bg)} matched background")

            # --- the table, full precision, ranked -------------------------
            for lab in ("9", "45"):
                tg, bg, per = store[lab]
                order = sorted((g for g in fitted if np.isfinite(per[g][0])),
                               key=lambda g: -per[g][0])
                print(f"\n      delta = median|t| target - median|t| background, "
                      f"{lab}-gene target set")
                print(f"        {'rank':<5}{'gene':<10}{'target':>9}{'bg':>9}{'delta':>9}")
                for i, g in enumerate(order, 1):
                    d, mt, mb, _, _ = per[g]
                    mark = "  <-- control floor" if g == "GAPDH" else ""
                    print(f"        {i:<5}{g:<10}{mt:>9.4f}{mb:>9.4f}{d:>+9.4f}{mark}")
                if "GAPDH" in per and "RBM20" in per:
                    rk = [g for g in order].index("RBM20") + 1 if "RBM20" in order else None
                    gk = [g for g in order].index("GAPDH") + 1 if "GAPDH" in order else None
                    print(f"        RBM20 ranks {rk} of {len(order)}; "
                          f"GAPDH ranks {gk} of {len(order)}")

            # --- decomposition of the 9 -> 45 change for RBM20 -------------
            d9, t9, b9, _, _ = store["9"][2]["RBM20"]
            d45, t45, b45, _, _ = store["45"][2]["RBM20"]
            if np.isfinite(d9) and np.isfinite(d45):
                print(f"\n      WHERE THE 9 -> 45 CHANGE CAME FROM (RBM20):")
                print(f"        delta {d9:+.4f} -> {d45:+.4f}  (change {d45-d9:+.4f})")
                print(f"        target median  {t9:.4f} -> {t45:.4f}  "
                      f"({t45-t9:+.4f} of the change)")
                print(f"        background     {b9:.4f} -> {b45:.4f}  "
                      f"({-(b45-b9):+.4f} of the change)")
                if abs(t45 - t9) < 0.5 * abs(b45 - b9):
                    print("        >>> the target statistic barely moved; the delta changed")
                    print("            because the MATCHED BACKGROUND moved. That is a")
                    print("            property of the background construction, not")
                    print("            evidence about the gene set.")

            # --- margin with uncertainty, 2-D and 3-D matching -------------
            if "GAPDH" not in fitted:
                continue
            beta_col = f"dis_{ctag}.disease.beta"
            have_beta = beta_col in res.columns
            for lab, gl in (("45", g45), ("9", g9)):
                is_t = res["gene"].isin(gl).values
                variants = [("2-D match (as published)", None)]
                if have_beta:
                    variants.append(("3-D match (+|disease beta|)",
                                     np.abs(res[beta_col].values)))
                for vname, extra in variants:
                    tg, bg = match_background(res, is_t, extra)
                    dr = delta(tv["RBM20"], tg, bg)[0]
                    dg = delta(tv["GAPDH"], tg, bg)[0]
                    lo, hi, p = boot_margin(tv["RBM20"], tv["GAPDH"], tg, bg,
                                            a.n_boot, rng)
                    m = dr - dg
                    if not np.isfinite(m):
                        continue
                    crosses = (not np.isfinite(lo)) or (lo <= 0 <= hi)
                    verd = ("INDISTINGUISHABLE from the control"
                            if crosses else
                            ("RBM20 ABOVE the control" if m > 0
                             else "RBM20 BELOW the control"))
                    print(f"\n      margin, {lab}-gene set, {vname}")
                    print(f"        RBM20 {dr:+.4f}   GAPDH {dg:+.4f}   "
                          f"margin {m:+.4f}")
                    print(f"        95% CI [{lo:+.4f}, {hi:+.4f}]   "
                          f"bootstrap p = {p:.3g}")
                    print(f"        >>> {verd}")
                    rows.append(dict(table=tag_name, model=f"{model}/{ctag}",
                                     gene_set=lab, matching=vname,
                                     rbm20_delta=round(dr, 5),
                                     gapdh_delta=round(dg, 5),
                                     margin=round(m, 5),
                                     ci_lo=round(lo, 5) if np.isfinite(lo) else "",
                                     ci_hi=round(hi, 5) if np.isfinite(hi) else "",
                                     boot_p=round(p, 5) if np.isfinite(p) else "",
                                     n_target=len(tg), n_bg=len(bg),
                                     verdict=verd))
                if not have_beta:
                    print(f"        [note] {beta_col} absent -- 3-D match not run on"
                          f" this table")

    if rows:
        pd.DataFrame(rows).to_csv(a.out, sep="\t", index=False)
        print(f"\n[write] {a.out}")

    print("\n" + "=" * 78)
    print("HOW TO DECIDE")
    print("=" * 78)
    print("Look at DCM/cell on the two *_techadj tables -- those are the models the")
    print("manuscript reports -- at the 45-gene set, 3-D match.")
    print()
    print("  * margin CI covers 0  ->  the null stands. Write it as 'no candidate was")
    print("    distinguishable from a housekeeping negative control', drop the")
    print("    '+0.03' phrasing from the Results, and report the CI instead. A margin")
    print("    quoted without an interval was the weakness that produced 56's")
    print("    incoherent verdicts.")
    print()
    print("  * margin CI excludes 0 and RBM20 outranks GAPDH  ->  the null does not")
    print("    stand at 45. The honest write-up is that the RBM20 target set shows a")
    print("    small but resolvable excess over matched background beyond the")
    print("    housekeeping floor, with the inflation caveat attached -- and the paper's")
    print("    second claim narrows from 'no regulator explains it' to 'no regulator")
    print("    explains more than a few percent of it'. Do NOT report the nine-gene")
    print("    null in that case; it would be the subset that happened to look clean.")
    print()
    print("  * the two tables disagree  ->  report both event types separately. SE and")
    print("    MXE are different quantities and the project already reports them apart.")
    print()
    print("If the null does not stand, the next run is a refit of 03 with --n-perm 200")
    print("and a wider control-gene panel, so the floor is an empirical distribution of")
    print("housekeeping genes rather than GAPDH alone. One gene is not a null")
    print("distribution, and that is the real limitation of the published design.")


if __name__ == "__main__":
    main()
