#!/usr/bin/env python3
"""
45_cohort_slopes.py -- emit the two tables `44_fill_recompute.py` still needs:

  cohort_slopes_dedup.tsv            one row per cohort: Deming slope, Pearson r, n
  cohort_type_concordance_dedup.tsv  one row per cohort x event type: same/n/pct

The first closes markers 2.2 / 2.3 / 3.10 / 1.2. The second closes 3.5, which
`direction_concordance_dedup.tsv` cannot answer because it has no event-type column.

IT DOES NOT INVENT AN ESTIMATOR
-------------------------------
Everything here is lifted from `38_figure4_replication.py`, which already computes
these slopes for its panel B: the same `_dedup`-preferring loader, the same exon-level
collapse (NKEEP), the same MAGNet selection (FDR<0.05 and |dPSI|>=0.10), and the same
Deming estimator with lambda from se_adj. If this script and Figure 4 ever disagree,
one of them has been edited and the figure is what the manuscript shows.

THE SELF-CHECK
--------------
Per-cohort direction concordance is computed here from the same merged event sets the
slopes are computed on, and compared against the values already in the project
(Colorado 79.6%, Brigham 81.2%, Stanford 76.1%, Fuwai 69.2%, pooled 1,636/2,163).
A match means the slopes are computed on the right events. A mismatch means they are
not, and the slopes must not be used — which is the failure this script exists to
catch, because a Deming slope on the wrong event set looks perfectly reasonable.

Usage:
  export PATH=/data/apps/anaconda3/envs/mamba/envs/rnaseq_tools/bin:$PATH
  cd /data/project/chureeratp/RNA/HEART
  python 45_cohort_slopes.py
  python 44_fill_recompute.py --conc cohort_type_concordance_dedup.tsv \
                              --slopes cohort_slopes_dedup.tsv
"""
import argparse
import os

import numpy as np
import pandas as pd

NKEEP = {"SE": 3, "MXE": 5, "A3SS": 5, "A5SS": 5, "RI": 3}

# --------------------------------------------------------------------------
# Read length. Only the MAGNet DISCOVERY table carries this suffix; REPL below
# stays unsuffixed because no replication cohort was re-quantified.
SUF = os.environ.get("TBL_SUFFIX", "")
# --------------------------------------------------------------------------
TYPES = ["SE", "MXE", "A3SS", "A5SS", "RI"]
REPL = [("Colorado", "{ET}_Colorado_DCMvsNF"), ("Brigham", "Brigham_{ET}_HCM"),
        ("Fuwai", "Fuwai_{ET}_HCM"), ("Stanford", "Stanford_{ET}_HCM"),
        ("LBNL", "{ET}_LBNL_DCMvsNF")]
NOT_INDEPENDENT = {"LBNL"}
DROP_COMBOS = {("Colorado", "MXE")}        # does not reproduce; excluded project-wide

# Already in the project, used only to verify the event set is right.
# THESE ARE SE-ONLY. direction_concordance_dedup.tsv holds one row per cohort and that
# row is SE: the four independent cohorts' SE rows sum to exactly 1,636/2,163, the
# figure the manuscript quotes. The first version of this script compared all-five-
# event-type totals against them and reported three false mismatches -- the numbers
# were right and the reference was mislabelled. Check like against like.
KNOWN_PCT = {"Colorado": 79.6, "Brigham": 81.2, "Stanford": 76.1, "Fuwai": 69.2,
             "LBNL": 75.9}
KNOWN_POOLED = (1636, 2163)

# --------------------------------------------------------------------------
# 125 bp overrides, placed AFTER the dicts they rewrite.
#
# Measured (33 section D, confirmed by 61):
#   KNOWN_PCT=Colorado=80.8,Brigham=81.9,Stanford=75.8,Fuwai=70.8,LBNL=78.2
#   KNOWN_POOLED=1376,1800
#
# The override must name EXACTLY the keys already present, so a typo or a
# missing cohort fails loudly and the error names the expected list.
_kp = os.environ.get("KNOWN_PCT", "").strip()
if _kp:
    _d = dict(_p.split("=", 1) for _p in _kp.split(","))
    assert set(_d) == set(KNOWN_PCT), (
        f"KNOWN_PCT override must name exactly {sorted(KNOWN_PCT)}, "
        f"got {sorted(_d)}")
    KNOWN_PCT = {_k: float(_v) for _k, _v in _d.items()}

_kpo = os.environ.get("KNOWN_POOLED", "").strip()
if _kpo:
    _k, _n = (int(_x) for _x in _kpo.split(","))
    KNOWN_POOLED = (_k, _n)
# --------------------------------------------------------------------------
REF_TYPE = "SE"


def load(base, et):
    for p in (f"{base}_dedup.all_events.tsv", f"{base}.all_events.tsv"):
        if os.path.exists(p) and os.path.getsize(p):
            t = pd.read_csv(p, sep="\t")
            d = next((c for c in ("dPSI_adj", "dPSI") if c in t.columns), None)
            if d is None:
                cand = [c for c in t.columns if c.lower().startswith("dpsi")]
                d = next((c for c in cand if "adj" in c.lower()),
                         cand[0] if cand else None)
            if d is None or "key" not in t.columns:
                return None, p
            if "se_adj" not in t.columns and "se" in t.columns:
                t = t.rename(columns={"se": "se_adj"})
            n = NKEEP[et]
            t["ekey"] = [":".join(str(k).split(":")[:n]) for k in t["key"]]
            t = t.loc[t.groupby("ekey")["FDR"].idxmin()]
            return t.rename(columns={d: "dpsi"})[["ekey", "dpsi", "se_adj", "FDR"]], p
    return None, base


def deming(x, y, sx, sy):
    ok = np.isfinite(x) & np.isfinite(y) & np.isfinite(sx) & np.isfinite(sy)
    x, y, sx, sy = x[ok], y[ok], sx[ok], sy[ok]
    if len(x) < 10:
        return np.nan
    lam = np.mean(sy ** 2) / np.mean(sx ** 2) if np.mean(sx ** 2) > 0 else 1.0
    dx, dy = x - x.mean(), y - y.mean()
    Sxx, Syy, Sxy = (dx ** 2).sum(), (dy ** 2).sum(), (dx * dy).sum()
    if Sxy == 0:
        return np.nan
    return float((Syy - lam * Sxx + np.sqrt((Syy - lam * Sxx) ** 2
                                            + 4 * lam * Sxy ** 2)) / (2 * Sxy))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--slope-type", default="SE",
                    help="event type the per-cohort slope is computed on "
                         "(SE, to match Figure 4)")
    ap.add_argument("--fdr", type=float, default=0.05)
    ap.add_argument("--min-dpsi", type=float, default=0.10)
    ap.add_argument("--out-slopes", default="cohort_slopes_dedup" + SUF + ".tsv")
    ap.add_argument("--out-conc", default="cohort_type_concordance_dedup" + SUF + ".tsv")
    a = ap.parse_args()

    rows_conc, rows_slope = [], []
    missing = []
    for et in TYPES:
        disc, dp = load(f"{et}_v2filtered" + SUF, et)
        if disc is None:
            missing.append(f"{et} discovery ({dp})")
            continue
        sig = disc[(disc.FDR < a.fdr) & (disc.dpsi.abs() >= a.min_dpsi)]
        D = sig[["ekey", "dpsi", "se_adj"]].rename(
            columns={"dpsi": "x", "se_adj": "sx"})
        for coh, pat in REPL:
            if (coh, et) in DROP_COMBOS:
                continue
            T, tp = load(pat.format(ET=et), et)
            if T is None:
                missing.append(f"{coh} {et} ({tp})")
                continue
            m = D.merge(T[["ekey", "dpsi", "se_adj"]].rename(
                columns={"dpsi": "y", "se_adj": "sy"}), on="ekey")
            if not len(m):
                continue
            same = int((np.sign(m.x) == np.sign(m.y)).sum())
            rows_conc.append({"cohort": coh, "event_type": et,
                              "same_direction": same, "n": len(m),
                              "pct": round(100 * same / len(m), 2)})
            if et == a.slope_type:
                sl = deming(m.x.values, m.y.values, m.sx.values, m.sy.values)
                r = float(np.corrcoef(m.x, m.y)[0, 1]) if len(m) > 2 else np.nan
                rows_slope.append({"cohort": coh, "deming": round(sl, 4),
                                   "pearson_r": round(r, 4), "n": len(m),
                                   "independent": coh not in NOT_INDEPENDENT})

    if missing:
        print(f"  [warn] {len(missing)} table(s) not found: "
              + ", ".join(missing[:8]) + (" ..." if len(missing) > 8 else ""))
    if not rows_conc:
        print("[fatal] nothing loaded. Run this in the directory holding the "
              "*_dedup.all_events.tsv tables.")
        return

    C = pd.DataFrame(rows_conc)
    S = pd.DataFrame(rows_slope)
    C.to_csv(a.out_conc, sep="\t", index=False)
    S.to_csv(a.out_slopes, sep="\t", index=False)
    print(f"[write] {a.out_conc}   ({len(C)} cohort x event-type rows)")
    print(f"[write] {a.out_slopes}  ({len(S)} cohorts, {a.slope_type} slopes)")

    print("\n" + "=" * 78)
    print("SELF-CHECK — do these event sets reproduce the known concordances?")
    print("=" * 78)
    bad = 0
    percoh_bad = pooled_bad = 0
    ref = C[C.event_type == REF_TYPE]
    per = ref.groupby("cohort").apply(
        lambda g: pd.Series({"same": g.same_direction.sum(), "n": g.n.sum()}),
        include_groups=False)
    print(f"  comparing {REF_TYPE} rows only — the reference figures are {REF_TYPE}")
    for coh, r in per.iterrows():
        pct = 100 * r["same"] / r["n"]
        want = KNOWN_PCT.get(coh)
        if want is None:
            print(f"  [    ] {coh:<10}{pct:6.1f}%   (no reference value)")
            continue
        d = abs(pct - want)
        ok = d < 0.6
        bad += 0 if ok else 1
        percoh_bad = bad
        print(f"  [{'ok ' if ok else '***'}] {coh:<10}{pct:6.1f}%   "
              f"project says {want:.1f}%   diff {d:.2f}")
    ind = ref[~ref.cohort.isin(NOT_INDEPENDENT)]
    k, n = int(ind.same_direction.sum()), int(ind.n.sum())
    dp = abs(100 * k / n - 100 * KNOWN_POOLED[0] / KNOWN_POOLED[1])
    # exact is ideal but not required: the pooled denominator can differ slightly if
    # a cohort x type combination is filtered differently here than in script 33.
    # Demanding an exact match would raise a false alarm on a rounding-level gap and
    # get the whole check ignored, which is worse than the gap.
    exact = (k, n) == KNOWN_POOLED
    close = dp < 0.6
    print(f"  [{'ok ' if close else '***'}] pooled     "
          f"{k:,}/{n:,} = {100 * k / n:.1f}%   project says "
          f"{KNOWN_POOLED[0]:,}/{KNOWN_POOLED[1]:,}"
          + ("   exact" if exact else f"   diff {dp:.2f}"))
    pooled_bad = 0 if close else 1

    # Distinguish the two failures. A per-cohort mismatch means that cohort's merged
    # event set is wrong and its slope is meaningless. A pooled-only mismatch, with
    # every cohort passing, means the per-cohort sets are right and only the MIX of
    # events across cohorts differs from script 33's — the slopes are still valid.
    # Reporting both as "do not use" would have condemned a correct result.
    if bad:
        print(f"\n  *** {bad} cohort(s) mismatch. Those cohorts' merged event sets are")
        print("      NOT the ones the project's figures were computed on, so their")
        print(f"      slopes in {a.out_slopes} are on the wrong events — do not use")
        print("      them. Most likely a cohort table missing above, dropping its rows")
        print("      out of the denominator.")
    elif pooled_bad:
        print("\n  [ok ] every cohort matches individually, so each slope is computed")
        print("        on the right events and the per-cohort numbers are usable.")
        print("  [note] only the POOLED total differs, which means the mix of events")
        print("        across cohorts is not identical to script 33's. That affects the")
        print("        pooled percentage, not any single cohort's slope. Quote the")
        print("        pooled figure from script 33, not from here.")
    else:
        print("\n  [ok ] every cohort matches, pooled included — the slopes are on the")
        print("        right events.")

    allt = C[~C.cohort.isin(NOT_INDEPENDENT)]
    ka, na = int(allt.same_direction.sum()), int(allt.n.sum())
    print(f"\n  For reference, ALL FIVE EVENT TYPES pooled over the four independent")
    print(f"  cohorts: {ka:,}/{na:,} = {100 * ka / na:.1f}%. This is a DIFFERENT and")
    print(f"  additional statistic from the {REF_TYPE}-only {KNOWN_POOLED[0]:,}/"
          f"{KNOWN_POOLED[1]:,} the manuscript quotes — do not substitute one for the")
    print("  other, and say which the text means.")

    print("\n" + "=" * 78)
    print(f"PER-COHORT {a.slope_type} SLOPES")
    print("=" * 78)
    for r in S.sort_values("deming").itertuples():
        mark = "" if r.independent else "   * not independent — excluded from ranges"
        print(f"  {r.cohort:<10}{r.deming:6.3f}   r={r.pearson_r:5.2f}   "
              f"n={r.n:,}{mark}")
    I = S[S.independent]
    if len(I):
        print(f"\n  independent cohorts: slope {I.deming.min():.2f}–{I.deming.max():.2f}"
              f",  mean {I.deming.mean():.3f},  r {I.pearson_r.min():.2f}–"
              f"{I.pearson_r.max():.2f}")
        print(f"  (the docs currently say 0.32–0.87 and r 0.48–0.72, computed with "
              f"LBNL)")

    print(f"\n  Now: python 44_fill_recompute.py --conc {a.out_conc} "
          f"--slopes {a.out_slopes}")


if __name__ == "__main__":
    main()
