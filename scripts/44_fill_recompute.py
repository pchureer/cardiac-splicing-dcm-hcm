#!/usr/bin/env python3
"""
44_fill_recompute.py -- resolve the ⚠ RECOMPUTE markers left in the project docs by
the Phase 0.5 sweep, from the post-exclusion tables rather than by hand.

Reads the per-cohort direction-concordance table (script 33's `_dedup` output) and,
if present, the per-cohort Deming slopes. Drops LBNL. Prints a paste-ready block
keyed to the marker numbers in PHASE05_edit_list.md.

WHAT IT RESOLVES
----------------
  3.4  pooled concordance, four cohorts          -- also the script's own self-check
  3.3  the DCM row (Colorado alone)
  3.5  "N of 18" cohort x event-type combinations individually significant
  2.2  the across-projects slope and concordance ranges
  2.3  the concordance range in the sentence for the paper
  1.2  mean slope over the four independent cohorts
  3.10 "MAGNet vs four other cohorts" row

THE SELF-CHECK IS THE POINT
---------------------------
The pooled four-cohort figure is already known independently: 1,636/2,163 = 75.6%
(RESULT_sample_identity_LBNL_overlap.md). If this script's parse of the table does not
reproduce that, its column detection is wrong and every other number it prints is
suspect. It says so loudly rather than printing a plausible wrong answer.

LBNL IS DROPPED, NOT ZEROED
---------------------------
LBNL rows are excluded from every statistic and reported separately, because they are a
technical replicate of the discovery cohort. Quietly folding them in at weight zero
would leave the denominators looking untouched.

Usage:
  export PATH=/data/apps/anaconda3/envs/mamba/envs/rnaseq_tools/bin:$PATH
  cd /data/project/chureeratp/RNA/HEART
  python 44_fill_recompute.py
  python 44_fill_recompute.py --conc direction_concordance_dedup.tsv --list
"""
import argparse
import glob
import os

import numpy as np
import pandas as pd

try:
    from scipy.stats import binomtest
    def binom_p(k, n):
        return float(binomtest(int(k), int(n), 0.5, alternative="greater").pvalue)
except Exception:                                   # normal approximation fallback
    from math import erfc, sqrt
    def binom_p(k, n):
        if n == 0:
            return np.nan
        z = (k - n / 2) / sqrt(n / 4)
        return float(0.5 * erfc(z / sqrt(2)))

EXCLUDE = {"LBNL"}
DCM_COHORTS = {"Colorado"}
HCM_COHORTS = {"Brigham", "Fuwai", "Stanford"}
KNOWN_POOLED = (1636, 2163)      # the independently computed four-cohort figure
# ...and it is SE-ONLY. A table carrying all five event types pools to a different
# number, which is not an error in either. Filter to REF_TYPE before the self-check.
REF_TYPE = "SE"
TYPES = ["SE", "MXE", "A3SS", "A5SS", "RI"]

# Colorado MXE does not reproduce and is excluded from every table and figure
DROP_COMBOS = {("Colorado", "MXE")}


def find(paths, patterns):
    out = [p for p in paths if p and os.path.exists(p) and os.path.getsize(p)]
    for pat in patterns:
        out += [f for f in sorted(glob.glob(pat))
                if os.path.isfile(f) and os.path.getsize(f) and f not in out]
    return out


def load_conc(path):
    """Per-cohort (and ideally per-event-type) direction concordance."""
    t = pd.read_csv(path, sep=None, engine="python")
    L = {str(c).lower(): c for c in t.columns}
    coh = next((L[k] for k in ("cohort", "dataset", "project") if k in L), None)
    same = next((L[k] for k in ("same_direction", "same_dir", "same", "k",
                                "n_same", "concordant") if k in L), None)
    tot = next((L[k] for k in ("n", "total", "n_total", "tested", "n_events",
                               "n_compared") if k in L), None)
    et = next((L[k] for k in ("event_type", "etype", "type") if k in L), None)
    if coh is None or same is None or tot is None:
        return None, (f"need cohort / same-direction / total columns; "
                      f"have {list(t.columns)[:10]}")
    d = pd.DataFrame({"cohort": t[coh].astype(str),
                      "same": pd.to_numeric(t[same], errors="coerce"),
                      "n": pd.to_numeric(t[tot], errors="coerce")})
    d["et"] = t[et].astype(str).str.upper() if et else "ALL"
    d = d.dropna(subset=["same", "n"])
    # a concordance count must satisfy 0 <= same <= n. A row that does not is a
    # parsing error (wrong column picked up) or a corrupt table, and silently
    # summing it would produce a plausible-looking wrong total.
    bad = d[(d.n <= 0) | (d["same"] < 0) | (d["same"] > d.n)]
    if len(bad):
        print(f"  *** {len(bad)} row(s) with an impossible count (same < 0, same > n "
              f"or n <= 0) — dropped, and the column detection is suspect:")
        for _, r in bad.head(6).iterrows():
            print(f"      {r.cohort:<12}{r.et if 'et' in bad.columns else '':<6}"
                  f"same={r['same']:g}  n={r.n:g}")
    d = d[(d.n > 0) & (d["same"] >= 0) & (d["same"] <= d.n)]
    return d, None


def load_slopes(path):
    t = pd.read_csv(path, sep=None, engine="python")
    L = {str(c).lower(): c for c in t.columns}
    coh = next((L[k] for k in ("cohort", "dataset") if k in L), None)
    sl = next((L[k] for k in ("deming", "slope", "deming_slope", "transfer_slope")
               if k in L), None)
    r = next((L[k] for k in ("pearson_r", "r", "pearson") if k in L), None)
    if coh is None or sl is None:
        return None
    d = pd.DataFrame({"cohort": t[coh].astype(str),
                      "slope": pd.to_numeric(t[sl], errors="coerce")})
    d["r"] = pd.to_numeric(t[r], errors="coerce") if r else np.nan
    return d.dropna(subset=["slope"])


def pct(k, n):
    return f"{k:,}/{n:,} ({100 * k / n:.1f}%)" if n else "—"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--conc", help="per-cohort direction concordance table")
    ap.add_argument("--slopes", help="per-cohort Deming slope table")
    ap.add_argument("--alpha", type=float, default=0.05)
    ap.add_argument("--list", dest="listonly", action="store_true")
    a = ap.parse_args()

    cands = find([a.conc], ["*direction_concordance*dedup*.tsv",
                            "*direction_concordance*.tsv", "*concordance*.tsv"])
    if not cands:
        print("[fatal] no direction-concordance table found. Point at it with --conc.")
        print("        Expected script 33's output, e.g. direction_concordance_dedup.tsv")
        return
    D = err = None
    used = None
    for f in cands:
        D, err = load_conc(f)
        if D is not None:
            used = f
            break
        print(f"  [skip] {os.path.basename(f)}: {err}")
    if D is None:
        print("[fatal] none of the candidates had the needed columns.")
        return
    print(f"[read] {used}   ({len(D)} rows, cohorts: "
          + ", ".join(sorted(D.cohort.unique())) + ")")

    S = None
    for f in find([a.slopes], ["*slope*.tsv", "master_numbers_flagged.tsv"]):
        S = load_slopes(f)
        if S is not None:
            print(f"[read] {f}   (slopes for {len(S)} cohorts)")
            break

    if a.listonly:
        print(D.to_string(index=False))
        return

    lb = D[D.cohort.isin(EXCLUDE)]
    D = D[~D.cohort.isin(EXCLUDE)]
    before = len(D)
    D = D[~D.apply(lambda r: (r.cohort, r.et) in DROP_COMBOS, axis=1)]
    if len(D) != before:
        print(f"  [info] dropped {before - len(D)} Colorado MXE row(s) — excluded "
              f"project-wide")

    # the headline replication figures are SE; per-event-type rows are used only for
    # the combination count further down
    ref = D[D.et == REF_TYPE] if (D.et != "ALL").any() else D
    if len(ref) and len(ref) != len(D):
        print(f"  [info] self-check and the headline rows use {REF_TYPE} only "
              f"({len(ref)} of {len(D)} rows); the combination count uses all types")
    dcm = ref[ref.cohort.isin(DCM_COHORTS)]
    hcm = ref[ref.cohort.isin(HCM_COHORTS)]
    k_all, n_all = int(ref["same"].sum()), int(ref.n.sum())

    print("\n" + "=" * 78)
    print("SELF-CHECK")
    print("=" * 78)
    print(f"  pooled, four cohorts:  {pct(k_all, n_all)}")
    print(f"  independently known:   {pct(*KNOWN_POOLED)}")
    if (k_all, n_all) == KNOWN_POOLED:
        print("  [ok ] exact match — the parse is right, trust what follows")
    else:
        d = abs(100 * k_all / n_all - 100 * KNOWN_POOLED[0] / KNOWN_POOLED[1])
        lvl = "close" if d < 0.5 else "*** DISAGREES"
        print(f"  [{lvl}] differs by {d:.2f} points.")
        print("  If this is not a rounding difference, the column detection or the")
        print("  event-type filtering is wrong and NOTHING BELOW SHOULD BE USED.")
    lbr = lb[lb.et == REF_TYPE] if (lb.et != "ALL").any() else lb
    if len(lbr):
        print(f"\n  LBNL, reported separately and in no total above: "
              f"{pct(int(lbr['same'].sum()), int(lbr.n.sum()))}")
    if (D.et != "ALL").any():
        ka, na = int(D["same"].sum()), int(D.n.sum())
        print(f"  all five event types pooled, four cohorts: {pct(ka, na)} "
              f"— a different statistic from the {REF_TYPE}-only headline")

    print("\n" + "=" * 78)
    print("PASTE-READY  —  keyed to PHASE05_edit_list.md")
    print("=" * 78)

    print("\n[3.4] PAPER_ONE_PAGE claim 4, pooled row")
    print(f"      | **all four pooled** | **{pct(k_all, n_all)}** | "
          f"p = {binom_p(k_all, n_all):.1e} |")

    print("\n[3.3] PAPER_ONE_PAGE claim 4, DCM row")
    if len(dcm):
        k, n = int(dcm["same"].sum()), int(dcm.n.sum())
        print(f"      | DCM cohort (Colorado) | {pct(k, n)} | "
              f"p = {binom_p(k, n):.1e} |")
    else:
        print("      *** no Colorado rows found")
    if len(hcm):
        k, n = int(hcm["same"].sum()), int(hcm.n.sum())
        print(f"      | HCM cohorts (Brigham, Fuwai, Stanford) | {pct(k, n)} | "
              f"p = {binom_p(k, n):.1e} |  (unchanged, for reference)")

    print("\n[3.5] \"N of M cohort x event-type combinations individually significant\"")
    if (D.et == "ALL").all():
        print("      *** the table has no event-type column, so this cannot be")
        print("      recomputed here. Rerun script 33 with per-event-type output.")
    else:
        rows = []
        for (c, e), g in D.groupby(["cohort", "et"]):
            k, n = int(g["same"].sum()), int(g.n.sum())
            rows.append((c, e, k, n, binom_p(k, n)))
        sig = [r for r in rows if r[4] < a.alpha]
        print(f"      {len(sig)} of {len(rows)} individually significant at "
              f"p < {a.alpha}")
        ns = [r for r in rows if r[4] >= a.alpha]
        if ns:
            print("      not significant: "
                  + ", ".join(f"{c} {e} ({k}/{n}, p={p:.2f})"
                              for c, e, k, n, p in ns))
        print(f"      (denominator was 23 with LBNL's five combinations included)")

    if S is not None:
        ind = S[~S.cohort.isin(EXCLUDE)]
        if len(ind):
            lo, hi = ind.slope.min(), ind.slope.max()
            print("\n[2.2 / 2.3 / 3.10] across-project ranges, LBNL removed")
            print(f"      Deming slope range: {lo:.2f}–{hi:.2f}   "
                  f"(was 0.32–0.87 with five cohorts)")
            print(f"      mean slope [1.2]:   {ind.slope.mean():.3f}   "
                  f"over {len(ind)} cohorts: "
                  + ", ".join(f"{r.cohort} {r.slope:.3f}"
                              for r in ind.itertuples()))
            if ind.r.notna().any():
                print(f"      Pearson r range:    {ind.r.min():.2f}–{ind.r.max():.2f}")
        per = ref.groupby("cohort")[["same", "n"]].sum()
        per = (100 * per["same"] / per["n"]).sort_values()
        print(f"      concordance range:  {per.min():.0f}–{per.max():.0f}%   "
              + ", ".join(f"{c} {v:.1f}%" for c, v in per.items()))
    else:
        print("\n[2.2 / 2.3 / 3.10 / 1.2] no slope table found — pass --slopes.")
        per = ref.groupby("cohort")[["same", "n"]].sum()
        per = (100 * per["same"] / per["n"]).sort_values()
        print(f"      concordance range, LBNL removed: {per.min():.0f}–{per.max():.0f}%")
        print("      " + ", ".join(f"{c} {v:.1f}%" for c, v in per.items()))

    print("\n" + "=" * 78)
    print("  Still NOT resolved by this script, because no file settles them:")
    print("   * [3.1]  what \"715 hearts\" counted")
    print("   * [3.6]  Colorado DCM-vs-ICM: take from script 33's _dedup output")
    print("   * [3.7]  MAGNet DCM-vs-HCM: 6 of 147,355 vs 0 of 45,380 — a convention")
    print("            decision, not a computation")
    print("=" * 78)


if __name__ == "__main__":
    main()
