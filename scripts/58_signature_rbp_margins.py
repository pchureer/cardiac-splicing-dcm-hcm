#!/usr/bin/env python3
"""
58_signature_rbp_margins.py -- put intervals on the inverted RBP scan.

WHY THIS EXISTS. The project runs two different tests and they are easy to confuse:

  TEST A  RBM20's own target genes vs matched background, with control genes scored
          through the same pipeline. Scripts 03 / 56 / 57. SETTLED 5 Oct: margins of
          -0.005 (SE) and +0.000 (MXE) against the GAPDH floor, 95% CIs covering zero.

  TEST B  the DISCOVERY SIGNATURE genes as the target set, eleven candidate regulators
          scored against it. Script 23. This is where "no candidate exceeded the GAPDH
          floor by more than 0.03" comes from -- the sentence in the Results, the
          abstract and the graphical abstract.

`57` fixed test A. Test B still reports bare deltas rounded to two decimals, and it
carries the identical two weaknesses: a single-gene floor, and a background matched on
mean PSI and PSI sd but not on disease effect size. A paper that gives a confidence
interval for test A and a bare delta for test B invites exactly one question, and the
answer should not be "we only did it for the test where it was asked".

WHAT THIS DOES
  1. REPRODUCTION CHECK against the published delta table, at full precision, before
     anything new is reported. If the recomputed deltas do not match the numbers in
     RESULT_rbm20_mechanism_null.md to within rounding, nothing downstream is used.
  2. Full-precision delta per candidate, ranked, with GAPDH marked.
  3. Paired bootstrap CI on (candidate delta - GAPDH delta) for EVERY candidate, with
     Holm adjustment across candidates. Multiplicity matters in the right direction
     here: the claim is a null, so the risk is one of eleven candidates crossing by
     chance, and Holm is the conservative guard against exactly that.
  4. The three-dimensional match, adding |disease beta|.
  5. A verdict that names any candidate whose interval clears the floor.

READ IT LIKE THIS. Every delta in this scan is positive in some model -- that is the
technical inflation the scan was built to expose, and GAPDH scoring is the evidence it
worked. The question is never "is the delta positive". It is "is any candidate
separable from a housekeeping gene", and only an interval answers that.

Read the DCM (within-disease) rows, cell-adjusted, as script 23's own header instructs:
target events were selected for large DCM-vs-NF dPSI, so an all-samples model lets any
disease-associated gene correlate with them for free.

  python 58_signature_rbp_margins.py
  python 58_signature_rbp_margins.py --n-boot 5000
"""

import argparse
import os
import re
import sys

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

CONTROL = "GAPDH"
MODELS = [("DCM", "cell"), ("DCM", "nocell"), ("ALL", "cell"), ("ALL", "nocell")]
REPORTED = ("DCM", "cell")
N_MATCH = 10          # script 03's default; changing it changes the background

# The published table from RESULT_rbm20_mechanism_null.md section 4, within-DCM,
# cell-adjusted. Used only as a reproduction check -- if these do not come back, the
# events table is not the one the manuscript was written from.
PUBLISHED = {
    "SE": {"PTBP1": 0.08, "QKI": 0.07, "HNRNPF": 0.07, "ZRANB2": 0.06,
           "HNRNPA0": 0.05, "MBNL1": 0.05, "GAPDH": 0.05, "SRSF1": 0.04,
           "RBM20": 0.03, "CELF1": 0.01, "RBM5": -0.01},
    "MXE": {"GAPDH": 0.10, "PTBP1": 0.10, "RBM20": 0.09, "HNRNPF": 0.00,
            "RBM5": -0.01, "CELF1": -0.01, "MBNL1": -0.02, "ZRANB2": -0.03,
            "HNRNPA0": -0.03, "QKI": -0.05, "SRSF1": -0.07},
}
REPRO_TOL = 0.015     # the published values are rounded to 2 dp; allow one unit plus a hair


def match_background(res, is_tgt, extra=None):
    """Script 03's matching, reproduced; optionally with a third dimension."""
    feats = [res["mean_psi"].values, res["sd_psi"].values]
    if extra is not None:
        feats.append(extra)
    f = np.column_stack(feats).astype(float)
    sd = np.nanstd(f, axis=0)
    sd[sd == 0] = 1.0
    fz = (f - np.nanmean(f, axis=0)) / sd
    ok = np.isfinite(fz).all(axis=1)
    bg = np.where(~is_tgt & ok)[0]
    tg = np.where(is_tgt & ok)[0]
    if not len(bg) or not len(tg):
        return tg, np.array([], int)
    _, nb = cKDTree(fz[bg]).query(fz[tg], k=min(N_MATCH, len(bg)))
    return tg, np.unique(bg[np.atleast_2d(nb).ravel()])


def delta(tv, tg, bg):
    a = tv[tg]; b = tv[bg]
    a = a[np.isfinite(a)]; b = b[np.isfinite(b)]
    if len(a) < 3 or len(b) < 3:
        return np.nan, np.nan, np.nan
    return float(np.median(a) - np.median(b)), float(np.median(a)), float(np.median(b))


def boot_margins(tv, genes, ctrl, tg, bg, n_boot, rng):
    """One resample of target rows and one of background rows per iteration, shared by
    every gene, so each candidate's margin against the control is paired. These genes'
    t-values are correlated across events through the library-level technical variation
    the scan exists to expose; resampling them independently would widen every interval
    and make the whole table look inconclusive for the wrong reason."""
    ok_t = np.ones(len(tg), bool)
    ok_b = np.ones(len(bg), bool)
    for g in genes:
        ok_t &= np.isfinite(tv[g][tg])
        ok_b &= np.isfinite(tv[g][bg])
    T = {g: tv[g][tg][ok_t] for g in genes}
    B = {g: tv[g][bg][ok_b] for g in genes}
    nt, nb = ok_t.sum(), ok_b.sum()
    if nt < 5 or nb < 5:
        return {g: (np.nan, np.nan, np.nan) for g in genes}, 0, 0
    draws = {g: np.empty(n_boot) for g in genes}
    for i in range(n_boot):
        it = rng.integers(0, nt, nt)
        ib = rng.integers(0, nb, nb)
        c = np.median(T[ctrl][it]) - np.median(B[ctrl][ib])
        for g in genes:
            draws[g][i] = (np.median(T[g][it]) - np.median(B[g][ib])) - c
    out = {}
    for g in genes:
        d = draws[g]
        lo, hi = np.percentile(d, [2.5, 97.5])
        p = 2 * min((d <= 0).mean(), (d >= 0).mean())
        out[g] = (float(lo), float(hi), float(min(p, 1.0)))
    return out, int(nt), int(nb)


def holm(pvals):
    """Holm-Bonferroni. Returns adjusted p in the input order."""
    items = sorted((p, i) for i, p in enumerate(pvals) if np.isfinite(p))
    m = len(items)
    adj = [np.nan] * len(pvals)
    run = 0.0
    for rank, (p, i) in enumerate(items):
        run = max(run, min(1.0, (m - rank) * p))
        adj[i] = run
    return adj


def fitted_genes(res, model, tag):
    pat = re.compile(rf"^{model}_([A-Z0-9]+)_{tag}\.\1\.t$")
    return sorted({m.group(1) for c in res.columns for m in [pat.match(c)] if m})


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tables", nargs="*",
                    default=["signature_rbp_SE.events.tsv",
                             "signature_rbp_MXE.events.tsv"])
    ap.add_argument("--target-file", default="compadj_multi_event_type.txt",
                    help="the signature gene list script 23 used as its target set")
    ap.add_argument("--n-boot", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--min-background", type=int, default=300,
                    help="script 23's own instruction: if the matched background is "
                         "not in the hundreds the enrichment test is void")
    ap.add_argument("--out", default="signature_rbp_margins.tsv")
    a = ap.parse_args()

    if not os.path.exists(a.target_file):
        sys.exit(f"[fatal] {a.target_file} not found. That file IS the target set; "
                 f"without it this script would silently test a different hypothesis. "
                 f"Rebuild it with 13_gene_lists_for_enrichment.py --tag compadj.")
    tgt_genes = sorted({t.strip() for t in open(a.target_file).read().split()
                        if t.strip()})
    print(f"[info] target set: {len(tgt_genes)} genes from {a.target_file}")

    tabs = [t for t in a.tables if os.path.exists(t)]
    missing = [t for t in a.tables if not os.path.exists(t)]
    for t in missing:
        print(f"[miss] {t}")
    if not tabs:
        sys.exit("[fatal] no signature_rbp_*.events.tsv found. Run "
                 "23_signature_rbp_scan.sbatch first.")

    rng = np.random.default_rng(a.seed)
    rows, verdicts = [], []

    for tab in tabs:
        res = pd.read_csv(tab, sep="\t", low_memory=False)
        need = {"gene", "mean_psi", "sd_psi"}
        if not need <= set(res.columns):
            print(f"[skip] {tab}: missing {sorted(need - set(res.columns))}")
            continue
        et = "MXE" if "MXE" in os.path.basename(tab) else "SE"
        print("\n" + "=" * 78)
        print(f"{os.path.basename(tab)}   ({len(res)} events, event type {et})")
        print("=" * 78)

        present = [g for g in tgt_genes if g in set(res["gene"].astype(str))]
        is_tgt = res["gene"].isin(tgt_genes).values
        print(f"  {len(present)}/{len(tgt_genes)} target genes contribute "
              f"{int(is_tgt.sum())} target events")
        if len(present) < len(tgt_genes):
            absent = [g for g in tgt_genes if g not in present]
            print(f"  absent from the tested set: {','.join(absent)}")
            print(f"  >>> the Methods must say {len(present)} genes, not {len(tgt_genes)}")

        for model, ctag in MODELS:
            genes = fitted_genes(res, model, ctag)
            if CONTROL not in genes:
                print(f"\n  --- {model}/{ctag}: {CONTROL} not fitted -- no floor, "
                      f"skipping ---")
                continue
            star = ("  <-- the model script 23 says to read"
                    if (model, ctag) == REPORTED else "")
            print(f"\n  --- {model}/{ctag} ---{star}")
            tv = {g: np.abs(res[f"{model}_{g}_{ctag}.{g}.t"].values) for g in genes}

            beta_col = f"dis_{ctag}.disease.beta"
            variants = [("2-D match (as published)", None)]
            if beta_col in res.columns:
                variants.append(("3-D match (+|disease beta|)",
                                 np.abs(res[beta_col].values)))
            else:
                print(f"      [note] {beta_col} absent -- 3-D match not available")

            for vname, extra in variants:
                tg, bg = match_background(res, is_tgt, extra)
                if len(bg) < a.min_background:
                    print(f"      {vname}: matched background only {len(bg)} events "
                          f"(< {a.min_background}). Script 23's own rule says the "
                          f"enrichment test is void here. NOT REPORTING.")
                    continue
                per = {g: delta(tv[g], tg, bg) for g in genes}
                ci, nt, nb = boot_margins(tv, genes, CONTROL, tg, bg, a.n_boot, rng)
                cand = [g for g in genes if g != CONTROL]
                adj = dict(zip(cand, holm([ci[g][2] for g in cand])))

                # reproduction check, published model only, published matching only
                if (model, ctag) == REPORTED and vname.startswith("2-D") and et in PUBLISHED:
                    pub = PUBLISHED[et]
                    bad = [(g, per[g][0], pub[g]) for g in genes
                           if g in pub and np.isfinite(per[g][0])
                           and abs(per[g][0] - pub[g]) > REPRO_TOL]
                    if bad:
                        print(f"      >>> REPRODUCTION CHECK FAILED against the "
                              f"published {et} table:")
                        for g, got, want in bad:
                            print(f"            {g:<9} recomputed {got:+.4f}, "
                                  f"published {want:+.2f}")
                        print("      >>> This events table is not the one the "
                              "manuscript was written from, or the matching has")
                        print("      >>> changed. Resolve before using anything below.")
                        verdicts.append((f"{et} {model}/{ctag}",
                                         "REPRODUCTION FAILED"))
                    else:
                        n_ck = sum(1 for g in genes if g in pub)
                        print(f"      reproduction: all {n_ck} published deltas "
                              f"recovered within {REPRO_TOL} -> MATCHES")

                order = sorted((g for g in genes if np.isfinite(per[g][0])),
                               key=lambda g: -per[g][0])
                print(f"\n      {vname}: {len(tg)} target / {len(bg)} background "
                      f"({nt}/{nb} usable in the bootstrap)")
                print(f"        {'rank':<5}{'gene':<9}{'target':>8}{'bg':>8}"
                      f"{'delta':>9}{'margin':>9}{'95% CI':>22}{'p':>8}{'Holm':>8}")
                for i, g in enumerate(order, 1):
                    d, mt, mb = per[g]
                    if g == CONTROL:
                        print(f"        {i:<5}{g:<9}{mt:>8.3f}{mb:>8.3f}{d:>+9.4f}"
                              f"{'floor':>9}{'':>22}{'':>8}{'':>8}")
                        continue
                    m = d - per[CONTROL][0]
                    lo, hi, p = ci[g]
                    ap_ = adj.get(g, np.nan)
                    print(f"        {i:<5}{g:<9}{mt:>8.3f}{mb:>8.3f}{d:>+9.4f}"
                          f"{m:>+9.4f}"
                          f"{f'[{lo:+.3f}, {hi:+.3f}]':>22}{p:>8.3g}"
                          f"{ap_:>8.3g}")
                    rows.append(dict(event_type=et, model=f"{model}/{ctag}",
                                     matching=vname, gene=g,
                                     delta=round(d, 5),
                                     floor=round(per[CONTROL][0], 5),
                                     margin=round(m, 5),
                                     ci_lo=round(lo, 5), ci_hi=round(hi, 5),
                                     boot_p=round(p, 5),
                                     holm_p=round(ap_, 5) if np.isfinite(ap_) else "",
                                     n_target=len(tg), n_bg=len(bg)))

                clear_raw = [g for g in cand
                             if np.isfinite(ci[g][0]) and ci[g][0] > 0]
                clear_holm = [g for g in clear_raw
                              if np.isfinite(adj.get(g, np.nan)) and adj[g] < 0.05]
                if not clear_raw:
                    v = "no candidate separates from the floor"
                elif not clear_holm:
                    v = (f"{','.join(clear_raw)} clear the floor raw but none survives "
                         f"Holm across {len(cand)} candidates")
                else:
                    v = f"CLEARS THE FLOOR after Holm: {','.join(clear_holm)}"
                print(f"        >>> {v}")
                if (model, ctag) == REPORTED:
                    verdicts.append((f"{et} {vname}", v))

    if rows:
        pd.DataFrame(rows).to_csv(a.out, sep="\t", index=False)
        print(f"\n[write] {a.out}")

    print("\n" + "=" * 78)
    print("VERDICT -- within-DCM, cell-adjusted, the model script 23 says to read")
    print("=" * 78)
    for k, v in verdicts:
        print(f"  {k:<34} {v}")
    if not verdicts:
        print("  nothing reportable -- check the notes above")

    print()
    print("WHAT TO WRITE")
    print("-" * 78)
    print("* no candidate separates, both event types  ->  the Results sentence becomes")
    print("  'no candidate regulator was distinguishable from a housekeeping negative")
    print("  control (largest margin X, 95% CI ...)'. Drop 'by more than 0.03' from the")
    print("  Results, the abstract and graphical-abstract block 3: a bare delta beside")
    print("  test A's interval is the inconsistency a reviewer will find.")
    print()
    print("* something clears Holm  ->  that gene is a positive and the paper's second")
    print("  claim narrows to name it. Do not bury it: follow script 23's own")
    print("  instruction and test that gene's OWN published targets as an independent")
    print("  target set, then check Colorado. A regulator that survives a housekeeping")
    print("  floor with an interval is a finding, not an inconvenience.")
    print()
    print("* raw-clear but Holm-fails  ->  report it as not surviving correction across")
    print("  eleven candidates, and say so explicitly rather than omitting it.")
    print()
    print("Either way the Methods gain one sentence: the floor is a single housekeeping")
    print("gene, so margins are reported with bootstrap intervals rather than as point")
    print("estimates. That is the limitation test A exposed and it applies here too.")


if __name__ == "__main__":
    main()
