#!/usr/bin/env python3
"""
36_matched_pair_analysis.py -- the decomposition. Called by 36_matched_pair_panel.sbatch,
but re-runnable on its own once the mp_* tables exist.

    A) MAGNet_matched -> LBNL_matched    SAME hearts, DIFFERENT platform  = TECHNICAL
    B) MAGNet_matched -> MAGNet_other    DIFFERENT hearts, SAME platform  = SAMPLING
    C) MAGNet_full    -> Colorado        both                             = C_OBSERVED (ref; script 33 section E — 0.521 at 50 bp, 0.600 at 125 bp)

All three at matched group sizes, on events selected from a MAGNet run that shares no
sample with any arm. Deming (errors-in-variables) slopes, lambda estimated from se_adj,
with a paired bootstrap CI over events.

READ A AND B TOGETHER, NOT SEPARATELY
-------------------------------------
Neither slope means much alone -- both arms are n=10 vs 14 and heavily attenuated by
regression to the mean, so both will come in below 1. The informative quantity is the
RATIO and the gap: B is the ceiling this design can reach with different people on the
same platform, and A measures what changing platform costs on top of that. If A and B
are close, platform is not the problem. If A is well below B, it is.
"""
import os

import numpy as np
import pandas as pd

NKEEP = {"SE": 3, "MXE": 5, "A3SS": 5, "A5SS": 5, "RI": 3}


SUF = os.environ.get("TBL_SUFFIX", "")   # "" = 50 bp, "_rl125" = 125 bp

# ---------------------------------------------------------------------------
# C -- the MAGNet -> Colorado cross-cohort slope -- is NOT computed in this
# file. It comes from `33_rebuild_summaries_dedup.py` section E.
#
# Until this block existed, C was the literal 0.521 in three places here, and
# the Colorado direction concordance was the literal '79.6%'. After the
# $TBL_SUFFIX fix this script computed A and B from the 125 bp tables and
# printed them against the 50 bp C:
#
#     A x B = 0.668   against the observed cross-cohort C = 0.521   -> 28% gap
#
# That gap was a cross-read-length comparison, not a result. With C recomputed
# at 125 bp (0.600) the gap is 11% -- the same as at 50 bp -- and the
# decomposition closes.
#
# C cannot be verified from this script's own inputs, so it cannot be gated.
# The rule instead: IF THE TABLES ARE SUFFIXED, C MUST BE SUPPLIED.
C_OBS = os.environ.get("C_OBSERVED", "").strip()
C_CONCORD = os.environ.get("C_CONCORDANCE", "").strip()
if SUF:
    if not C_OBS or not C_CONCORD:
        raise SystemExit(
            "\nREFUSING — TBL_SUFFIX is set to '" + SUF + "', so A and B will "
            "be computed from\n  re-quantified tables. C cannot be computed "
            "here; it comes from script 33's\n  section E, and the 50 bp "
            "default (0.521) would make `A x B against C` a\n  "
            "cross-read-length comparison. That is what produced the spurious "
            "28% gap.\n\n"
            "  Re-run 33 against the same tables:\n\n"
            "    python3 33_rebuild_summaries_dedup.py \\\n"
            "        --discovery '{ET}_v2filtered" + SUF + "' \\\n"
            "        --out summaries_dedup" + SUF + ".txt\n\n"
            "  then pass its Colorado row from sections E and D:\n\n"
            "    C_OBSERVED=<slope> C_CONCORDANCE=<pct> TBL_SUFFIX=" + SUF
            + " \\\n    python3 36_matched_pair_analysis.py\n\n"
            "  At 125 bp those are 0.600 and 80.8%.\n")
else:
    C_OBS = C_OBS or "0.521"
    C_CONCORD = C_CONCORD or "79.6%"
C_OBS_F = float(C_OBS)
# ---------------------------------------------------------------------------


def load(base, et):
    p = f"{base}{SUF}.all_events.tsv"
    if not (os.path.exists(p) and os.path.getsize(p)):
        return None
    t = pd.read_csv(p, sep="\t")
    d = next((c for c in ("dPSI_adj", "dPSI") if c in t.columns), None)
    if d is None or "key" not in t.columns:
        return None
    n = NKEEP[et]
    t["ekey"] = [":".join(str(k).split(":")[:n]) for k in t["key"]]
    t = t.loc[t.groupby("ekey")["FDR"].idxmin()]
    return t[["ekey", "gene", d, "se_adj", "FDR"]].rename(
        columns={d: "y", "se_adj": "se"})


def deming(x, y, sx, sy):
    ok = np.isfinite(x) & np.isfinite(y) & np.isfinite(sx) & np.isfinite(sy)
    x, y, sx, sy = x[ok], y[ok], sx[ok], sy[ok]
    if len(x) < 8:
        return np.nan
    lam = np.mean(sy ** 2) / np.mean(sx ** 2) if np.mean(sx ** 2) > 0 else 1.0
    dx, dy = x - x.mean(), y - y.mean()
    Sxx, Syy, Sxy = (dx ** 2).sum(), (dy ** 2).sum(), (dx * dy).sum()
    if Sxy == 0:
        return np.nan
    return float((Syy - lam * Sxx + np.sqrt((Syy - lam * Sxx) ** 2
                                            + 4 * lam * Sxy ** 2)) / (2 * Sxy))


def boot(x, y, sx, sy, n=2000, seed=0):
    rng = np.random.default_rng(seed)
    out = []
    for _ in range(n):
        i = rng.integers(0, len(x), len(x))
        s = deming(x[i], y[i], sx[i], sy[i])
        if np.isfinite(s):
            out.append(s)
    if len(out) < 100:
        return np.nan, np.nan
    return float(np.percentile(out, 2.5)), float(np.percentile(out, 97.5))


def main():
    types = os.environ.get("TYPES", "SE").split()
    print()
    for et in types:
        S = load(f"mp_{et}_selection", et)
        if S is None:
            print(f"[{et}] selection table missing -- run the sbatch first")
            continue
        sig = S[(S.FDR < 0.05) & (S.y.abs() >= 0.10)]
        print("=" * 88)
        print(f"{et}   selection: {len(sig)} events (FDR<0.05, |dPSI|>=0.10) from a "
              f"MAGNet run\n     sharing no sample with any arm")
        print("=" * 88)
        if len(sig) < 20:
            print("  too few selected events -- the decomposition would be noise")
            continue

        ref = load(f"mp_{et}_magnet_matched", et)
        arms = [("A  technical   MAGNet_matched -> LBNL_matched",
                 f"mp_{et}_lbnl_matched"),
                ("B  sampling    MAGNet_matched -> MAGNet_other",
                 f"mp_{et}_magnet_other")]
        if ref is None:
            print("  MAGNet_matched table missing")
            continue
        base = sig[["ekey"]].merge(
            ref.rename(columns={"y": "x", "se": "sx"}), on="ekey")

        # A and B are otherwise computed on DIFFERENT event sets: LBNL tests far fewer
        # events than MAGNet, so arm A is restricted to the higher-expressed subset that
        # survives in LBNL. Comparing a slope measured on easy events with one measured
        # on all events is not a fair comparison, so everything is also reported on the
        # COMMON set -- events present in all three arms. Quote the common-set numbers.
        tabs = {lab[0]: load(t, et) for lab, t in arms}
        common = set(base.ekey)
        for T in tabs.values():
            if T is not None:
                common &= set(T.ekey)
        common = sorted(common)
        print(f"  {'arm':<52}{'n':>6}{'slope':>9}{'95% CI':>22}{'same dir':>16}")
        res, res_c = {}, {}
        for label, tab in arms:
            T = tabs.get(label[0])
            if T is None:
                print(f"  {label:<52}{'table missing':>6}")
                continue
            m = base.merge(T[["ekey", "y", "se"]], on="ekey")
            if len(m) < 20:
                print(f"  {label:<52}{len(m):>6}  too few matched events")
                continue
            sl = deming(m.x.values, m.y.values, m.sx.values, m.se.values)
            lo, hi = boot(m.x.values, m.y.values, m.sx.values, m.se.values)
            same = int((np.sign(m.x) == np.sign(m.y)).sum())
            res[label[0]] = sl
            print(f"  {label:<52}{len(m):>6}{sl:>9.3f}"
                  f"{f'{lo:.3f} to {hi:.3f}':>22}"
                  f"{f'{same}/{len(m)} ({100*same/len(m):.0f}%)':>16}")
            c = m[m.ekey.isin(common)]
            if len(c) >= 20:
                slc = deming(c.x.values, c.y.values, c.sx.values, c.se.values)
                loc, hic = boot(c.x.values, c.y.values, c.sx.values, c.se.values)
                samec = int((np.sign(c.x) == np.sign(c.y)).sum())
                res_c[label[0]] = slc
                print(f"  {'   on the COMMON event set':<52}{len(c):>6}{slc:>9.3f}"
                      f"{f'{loc:.3f} to {hic:.3f}':>22}"
                      f"{f'{samec}/{len(c)} ({100*samec/len(c):.0f}%)':>16}")
        print(f"  {'C  reference   MAGNet_full -> Colorado (script 33)':<52}"
              f"{'':>6}{C_OBS_F:>9.3f}{'':>22}{C_CONCORD:>16}")

        if len(res_c) == 2:
            print()
            print(f"  COMMON EVENT SET ({len(common)} events in all three arms) -- "
                  f"quote these: A={res_c['A']:.3f}  B={res_c['B']:.3f}  "
                  f"A/B={res_c['A'] / res_c['B']:.2f}")
        use = res_c if len(res_c) == 2 else res
        if "A" in use and "B" in use and np.isfinite(use["A"]) and use["B"] > 0:
            A_, B_ = use["A"], use["B"]
            r = A_ / B_
            print()
            print(f"  A / B = {r:.2f}   "
                  f"(platform costs {100*(1-r):.0f}%, changing people costs "
                  f"{100*(1-B_):.0f}%)")
            print(f"  A x B = {A_ * B_:.3f}   against the observed cross-cohort C = {C_OBS_F:.3f}")
            print("  If those two agree, the decomposition closes: cross-cohort")
            print("  attenuation is the product of a platform term and a sampling term.")
            if r >= 0.85 and B_ >= 0.85:
                print("  -> Changing platform costs almost nothing once sample size is")
                print("     held fixed. The cross-cohort attenuation is NOT technical;")
                print("     it is sampling variability and regression to the mean.")
                print("     RESULT_magnitude_not_estimable.md stands as written.")
            elif r < B_ - 0.15:
                print(f"  -> PLATFORM DOMINATES. Changing people costs {100*(1-B_):.0f}%;")
                print(f"     changing platform costs a further {100*(1-r):.0f}%. The same hearts")
                print("     measured two ways disagree nearly as much as different")
                print("     hearts do, so cross-cohort magnitude differences are")
                print("     largely a measurement artefact. This converts")
                print("     'magnitude is not estimable' into a quantified technical")
                print("     attenuation term -- the stronger methods result.")
            else:
                print("  -> Platform and sampling contribute comparably. Report both")
                print("     slopes and resist attributing the gap to either alone.")
        print()
        print("  CAVEATS, and they are not small:")
        print("   * both arms are 10 DCM vs 14 NF, so both are heavily attenuated by")
        print("     regression to the mean. Neither slope is a magnitude estimate;")
        print("     only their comparison is interpretable.")
        # Read length is a bundled component ONLY when MAGNet is 125 bp.
        # At 50 bp against LBNL's 51 bp it was effectively matched, so the
        # unsuffixed wording is not a simplification -- it was correct.
        _rl = "125bp" if SUF else "50bp"
        _bundle = ("prep, layout, READ LENGTH (125 vs 51 bp) and pipeline "
                   "batch") if SUF else "prep, layout and pipeline batch"
        print(f"   * the LBNL arm is single-end 51bp against MAGNet's paired {_rl}, so")
        print(f"     'platform' bundles {_bundle} together. This")
        print("     design cannot separate them.")
        print("   * arm B is one deterministic draw of non-matched samples. Re-run with")
        print("     a different draw before quoting B to two decimal places.")
        print()


if __name__ == "__main__":
    main()
