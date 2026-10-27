#!/usr/bin/env python3
"""
50_readlength_impact.py -- did declaring --readLength 50 on 125 bp reads change anything?

THE QUESTION
------------
Every MAGNet rMATS run was invoked with `--readLength 50`; the reads are 125 bp
(48_provenance_versions.sh section 5, confirmed across samples by 49_readlength_audit.sh).
rMATS normalises PSI by effective inclusion- and skipping-form lengths computed from
--readLength, so the declared value enters PSI itself. This script measures the
consequence instead of arguing about it.

WHAT IT COMPARES, AND WHY THE TWO STAGES ARE DIFFERENT
-----------------------------------------------------
Stage A, RAW rMATS, old vs new: same samples, same (no) covariate model, only the
normalisation constant differs. This is the clean isolation of the readLength effect and
is the test that matters. Run it as soon as the rMATS jobs finish.

Stage B, the PROJECT'S numbers: the manuscript's 926 SE / 1,445 all-types come from
`phase2_diffsplice_v2.py` with covariates and the locked exon-level convention, NOT from
raw rMATS. Stage A's counts are therefore NOT comparable to 926 -- they are a like-for-
like old-vs-new comparison at a simpler threshold. To compare against the frozen numbers
you must first rerun phase2 on the new rMATS output and pass both adjusted tables with
--old-adj/--new-adj. The script refuses to pretend otherwise.

EXPECTATION, RECORDED BEFORE THE RESULT
---------------------------------------
Within MAGNet both arms carry the same constant, so dPSI direction should be nearly
unchanged and the counts should shift modestly as events cross |dPSI| >= 0.10. If that
is what comes back, the frozen numbers stand and the Methods discloses the discrepancy.
If direction moves, every MAGNet number in the paper is provisional. The verdict at the
end applies fixed thresholds so the answer is not chosen after seeing the data.

Usage:
  export PATH=/data/apps/anaconda3/envs/mamba/envs/rnaseq_tools/bin:$PATH
  cd /data/project/chureeratp/RNA/HEART

  # stage A, straight after the rMATS jobs
  python 50_readlength_impact.py --old rmats_magnet --new rmats_magnet_rl125

  # stage B, after rerunning phase2 on the new output
  python 50_readlength_impact.py --old rmats_magnet --new rmats_magnet_rl125 \
      --old-adj magnet_SE_adjusted.tsv --new-adj magnet_rl125_SE_adjusted.tsv
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

try:
    from scipy.stats import spearmanr, binomtest
except ImportError:
    print("scipy is required (statsmodels is NOT needed and is not installed here)")
    raise

# Exon-level join keys, mirroring the project's collapse rule. These are built from the
# coordinate columns because raw rMATS output has no `key` column -- that is added later
# by phase2. Used for joining and counting here only; the manuscript convention is
# whatever 24_lock_event_convention.py defines.
KEYCOLS = {
    "SE":   ["chr", "strand", "exonStart_0base", "exonEnd"],
    "MXE":  ["chr", "strand", "1stExonStart_0base", "1stExonEnd",
             "2ndExonStart_0base", "2ndExonEnd"],
    "A3SS": ["chr", "strand", "longExonStart_0base", "longExonEnd",
             "shortES", "shortEE"],
    "A5SS": ["chr", "strand", "longExonStart_0base", "longExonEnd",
             "shortES", "shortEE"],
    "RI":   ["chr", "strand", "riExonStart_0base", "riExonEnd"],
}
TYPES = ["SE", "MXE", "A3SS", "A5SS", "RI"]

# fixed verdict thresholds -- set here, before any data is seen
SIGN_OK, RHO_OK, COUNT_OK = 0.99, 0.98, 0.05      # frozen numbers stand
SIGN_DIR, RHO_DIR = 0.95, 0.90                    # direction still safe


def load(d, et, which="JC"):
    """One rMATS result table, collapsed exon-level by lowest FDR.

    Says WHICH directory a file is missing from. The first version reported only
    "missing in one of the two runs", which cannot distinguish "the new run has not
    finished" from "--old points at a directory that was never the original output"."""
    for name in (f"{et}.MATS.{which}.txt", f"{et}.MATS.{which}.txt.gz"):
        p = os.path.join(d, name)
        if os.path.exists(p) and os.path.getsize(p):
            break
    else:
        print(f"    [none] {d}/{et}.MATS.{which}.txt — absent or empty")
        return None
    t = pd.read_csv(p, sep="\t")
    need = KEYCOLS[et] + ["FDR", "IncLevelDifference", "IncFormLen", "SkipFormLen"]
    missing = [c for c in need if c not in t.columns]
    if missing:
        print(f"  [skip] {p}: missing {missing}")
        return None
    t["ekey"] = t[KEYCOLS[et]].astype(str).agg(":".join, axis=1)
    t = t.loc[t.groupby("ekey")["FDR"].idxmin()].set_index("ekey")
    return t


def pct(x, n):
    return f"{x:,} ({100.0 * x / n:.1f}%)" if n else "0"


def compare(et, old, new, fdr, mind):
    """Old vs new on the shared event set. Returns the numbers the verdict uses."""
    shared = old.index.intersection(new.index)
    n = len(shared)
    print(f"  {et}: {len(old):,} old / {len(new):,} new events, "
          f"{n:,} shared exon-level")
    if n < 50:
        print("    too few shared events to judge")
        return None
    o, w = old.loc[shared], new.loc[shared]

    # 1. what the constant actually did to the normalisation
    for col in ("IncFormLen", "SkipFormLen"):
        mo, mn = float(o[col].median()), float(w[col].median())
        if mo:
            print(f"    {col:<12} median {mo:>8.0f} -> {mn:>8.0f}   "
                  f"ratio {mn / mo:.3f}")
        else:
            print(f"    {col:<12} old median is 0 — cannot form a ratio")
    # the ratio of the two forms is what PSI actually depends on
    ro = (o.IncFormLen / o.SkipFormLen).median()
    rn = (w.IncFormLen / w.SkipFormLen).median()
    print(f"    Inc/Skip     median {ro:>8.3f} -> {rn:>8.3f}   "
          f"ratio {rn / ro:.3f}   <- this is what moves PSI")

    # 2. agreement of the effect estimates
    do, dn = o.IncLevelDifference.astype(float), w.IncLevelDifference.astype(float)
    ok = do.notna() & dn.notna()
    do, dn = do[ok], dn[ok]
    rho = spearmanr(do, dn).correlation if len(do) > 2 else np.nan
    nz = (do != 0) & (dn != 0)
    same = int((np.sign(do[nz]) == np.sign(dn[nz])).sum())
    sign_frac = same / max(nz.sum(), 1)
    slope = float(np.polyfit(do, dn, 1)[0]) if len(do) > 2 else np.nan
    print(f"    dPSI         Spearman rho {rho:.4f}   same direction "
          f"{same:,}/{int(nz.sum()):,} ({100 * sign_frac:.2f}%)   "
          f"OLS new~old {slope:.3f}")
    print(f"    |dPSI|       median {do.abs().median():.4f} -> "
          f"{dn.abs().median():.4f}")

    # 3. counts at the locked rule, applied to raw rMATS on BOTH runs equally
    so = (o.FDR < fdr) & (o.IncLevelDifference.abs() >= mind)
    sn = (w.FDR < fdr) & (w.IncLevelDifference.abs() >= mind)
    co, cn = int(so.sum()), int(sn.sum())
    delta = (cn - co) / co if co else np.nan
    gained = int((~so & sn).sum())
    lost = int((so & ~sn).sum())
    print(f"    significant  {co:,} -> {cn:,}   change {100 * delta:+.1f}%")
    print(f"                 gained {pct(gained, n)}   lost {pct(lost, n)}   "
          f"kept {pct(int((so & sn).sum()), n)}")
    if co:
        b = binomtest(gained, gained + lost, 0.5) if (gained + lost) else None
        if b is not None:
            print(f"    churn        gained vs lost binomial P = {b.pvalue:.3g} "
                  f"(symmetric churn is threshold noise; one-sided is a real shift)")
    return dict(type=et, n=n, rho=rho, sign=sign_frac, slope=slope,
                old=co, new=cn, delta=delta, gained=gained, lost=lost)


# the project's locked convention, as in 44/45/46
NKEEP = {"SE": 3, "MXE": 5, "A3SS": 5, "A5SS": 5, "RI": 3}
FROZEN_SE = 926          # MAGNet DCM vs NF, exon-level, FDR<0.05 and |dPSI|>=0.10


def adj_count(path, et, fdr, mind):
    """Count a phase2-adjusted table the way the manuscript counts it."""
    t = pd.read_csv(path, sep="\t")
    d = next((c for c in ("dPSI_adj", "dPSI") if c in t.columns), None)
    if d is None:
        cand = [c for c in t.columns if c.lower().startswith("dpsi")]
        d = next((c for c in cand if "adj" in c.lower()), cand[0] if cand else None)
    if d is None or "key" not in t.columns or "FDR" not in t.columns:
        return None, None
    n = NKEEP.get(et, 3)
    t["ekey"] = [":".join(str(k).split(":")[:n]) for k in t["key"]]
    t = t.loc[t.groupby("ekey")["FDR"].idxmin()]
    return int(((t.FDR < fdr) & (t[d].abs() >= mind)).sum()), len(t)


def frozen_check(old_adj, new_adj, et, fdr, mind, probe=False):
    """Stage B: the project's own numbers, if the adjusted tables were supplied."""
    print("\n" + "=" * 78)
    print("STAGE B — against the FROZEN manuscript numbers")
    print("=" * 78)
    print("  Stage A is raw rMATS and is NOT comparable to 926 / 1,445, which are")
    print("  covariate-adjusted and exon-level from phase2_diffsplice_v2.py.")
    if not (old_adj and new_adj):
        print("\n  --old-adj/--new-adj not supplied, so this stage was skipped.")
        print("  To run it: rerun phase2_diffsplice_v2.py on the NEW rMATS output, then")
        print("  pass both adjusted tables. Until then the frozen numbers are untested")
        print("  and no manuscript count should be changed on Stage A alone.")
        return None
    out = {}
    for lbl, p in (("old (rl50)", old_adj), ("new (rl125)", new_adj)):
        if not (os.path.exists(p) and os.path.getsize(p)):
            print(f"  *** {p} not found")
            return None
        s, n = adj_count(p, et, fdr, mind)
        if s is None:
            print(f"  *** {p}: no key/FDR/dPSI columns — not a phase2 table?")
            return None
        out[lbl] = (s, n)
        print(f"  {lbl:<12} {s:>6,} significant of {n:,} tested   ({p})")
    (so, _), (sn, _) = out["old (rl50)"], out["new (rl125)"]
    print(f"\n  change {sn - so:+,} ({100 * (sn - so) / so:+.1f}%)" if so else "")
    if probe:
        print("\n  --probe: these are SUBSET tables. The frozen-number check is skipped,")
        print("  and these counts are NOT the manuscript's counts and must not be quoted.")
        print("  What is being tested here is whether the covariate-adjusted pipeline")
        print("  reacts to the read-length constant any differently than raw rMATS did.")
    elif et == "SE":
        flag = "ok " if so == FROZEN_SE else "***"
        print(f"  [{flag}] the old table should reproduce the frozen {FROZEN_SE:,}; "
              f"it gives {so:,}")
        if so != FROZEN_SE:
            print("        If it does not, --old-adj is not the table the manuscript's")
            print("        numbers came from, and this whole comparison is against the")
            print("        wrong baseline. Fix that before reading anything into it.")
    return (so, sn)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--old", default="rmats_magnet", help="rl50 rMATS output dir")
    ap.add_argument("--new", default="rmats_magnet_rl125", help="rl125 rMATS output dir")
    ap.add_argument("--which", default="JC", choices=["JC", "JCEC"])
    ap.add_argument("--types", nargs="+", default=TYPES)
    ap.add_argument("--fdr", type=float, default=0.05)
    ap.add_argument("--min-dpsi", type=float, default=0.10)
    ap.add_argument("--old-adj", default=None)
    ap.add_argument("--new-adj", default=None)
    ap.add_argument("--probe", action="store_true",
                    help="the adjusted tables are from a SUBSET run (52_rl_probe), not "
                         "the full cohort. Skips the check that --old-adj reproduces "
                         "the frozen 926, which a subset cannot and should not.")
    a = ap.parse_args()

    print("=" * 78)
    print("STAGE A — raw rMATS, --readLength 50 versus 125, same samples")
    print("=" * 78)
    print(f"  old: {a.old}\n  new: {a.new}\n  quantification: {a.which}, "
          f"rule FDR<{a.fdr} and |dPSI|>={a.min_dpsi}")
    for d in (a.old, a.new):
        if not os.path.isdir(d):
            print(f"\n  *** {d} is not a directory. Nothing to compare.")
            sys.exit(1)

    # what is actually in each directory, before anything is compared
    for d in (a.old, a.new):
        got = sorted(f for f in os.listdir(d) if ".MATS." in f)
        print(f"\n  {d}: {len(got)} MATS file(s)" + (f" — {', '.join(got[:6])}" if got
                                                     else " — EMPTY, still running?"))

    res = []
    for et in a.types:
        o, w = load(a.old, et, a.which), load(a.new, et, a.which)
        if o is None or w is None:
            print(f"  {et}: skipped (see the [none] line above for which run)")
            continue
        r = compare(et, o, w, a.fdr, a.min_dpsi)
        if r:
            res.append(r)

    if not res:
        print("\n  nothing comparable was found.")
        print("  If the NEW directory is empty, 51_rmats_magnet_rl125.sbatch has not")
        print("  finished — wait for 'rmats exit status: 0' in its .out file.")
        print("  If the OLD directory is empty, --old is not the original output. Find")
        print("  it with:  ls -d rmats_* */rmats_* 2>/dev/null")
        sys.exit(1)

    se = next((r for r in res if r["type"] == "SE"), res[0])
    stage_b = frozen_check(a.old_adj, a.new_adj, se["type"], a.fdr, a.min_dpsi,
                           probe=a.probe)

    # ------------------------------------------------------------------
    print("\n" + "=" * 78)
    print("VERDICT  (thresholds fixed in the script before any data was seen)")
    print("=" * 78)
    print(f"  judged on {se['type']}: rho {se['rho']:.4f}, direction "
          f"{100 * se['sign']:.2f}%, count change {100 * se['delta']:+.1f}%")
    print(f"  thresholds: rho>={RHO_OK}, direction>={100*SIGN_OK:.0f}%, "
          f"|count change|<={100*COUNT_OK:.0f}%")

    strong = (se["rho"] >= RHO_OK and se["sign"] >= SIGN_OK
              and abs(se["delta"]) <= COUNT_OK)
    direction_ok = se["rho"] >= RHO_DIR and se["sign"] >= SIGN_DIR

    if strong:
        print("\n  [A] THE DECLARED READ LENGTH DID NOT MATERIALLY CHANGE THE RESULT.")
        print("      Keep the frozen 926 / 1,445 / 0. Still required:")
        print("        * fix Table 1 (MAGNet is paired 125 bp) and rebuild Figure 1")
        print("        * rewrite the attenuation sentence — it tracks the DECLARED")
        print("          read length, not the sequenced one")
        print("        * disclose declared vs observed, and this check, in the Methods")
    elif direction_ok:
        print("\n  [B] DIRECTION HOLDS, COUNTS MOVE.")
        print("      Every claim resting on direction is safe: both etiology nulls, the")
        print("      double positive control, direction concordance, the regulator null,")
        print("      the enrichment. But the MAGNet COUNTS are provisional.")
        print("      Rerun phase2 on the new output, then 24/26 to rebuild the master")
        print("      file, then 05/07, 21/33, and every figure that quotes a count.")
        print("      Rerun this script with --old-adj/--new-adj before changing text.")
    else:
        print("\n  [C] DIRECTION MOVED. STOP AND RECOMPUTE EVERYTHING MAGNET.")
        print("      This is the outcome the expectation above did not predict, so do")
        print("      not rationalise it. Check first that the two runs really are the")
        print("      same samples in the same order (--b1/--b2 identical), then treat")
        print("      the rl125 run as the truth and rebuild the whole MAGNet chain.")

    if stage_b is None:
        print("\n  NOTE: this verdict rests on Stage A only. It says what the read-length")
        print("  constant did to raw rMATS, not what it did to the manuscript's counts.")
        print("  Rerun with --old-adj/--new-adj before changing any number in the paper.")

    print("\n  Whatever the verdict: the Methods states the declared --readLength, the")
    print("  observed read length, and what this comparison showed. A reviewer who")
    print("  opens PRJNA595151 sees 125 bp on the first screen.")

    out = pd.DataFrame(res)
    out.to_csv("readlength_impact.tsv", sep="\t", index=False)
    print("\n[write] readlength_impact.tsv")


if __name__ == "__main__":
    main()
