#!/usr/bin/env python3
"""
116_reconcile_analytic_null.py -- resolve 0.0538 (script 95) vs 0.0540 (script 77)
for the SAME quantity: the analytic null of |rho| in the within-DCM LVEF fit.

RESOLVED. Both values are explained, and the first version of this script was
testing the WRONG FUNCTIONAL.

THE ERROR IN v1
---------------
v1 tested E|rho| = sqrt(2 / (pi * df)) -- the null's MEAN. At the measured
n (median 157) that is 0.0637, which explains neither printed value, so v1
reported "neither reproduces". But the statistic being compared is the
**median** |rho| (0.0555), so the comparator must be the null's **median**:

    median |rho|  =  z(0.75) / sqrt(df)  =  0.6744898 / sqrt(df)

because rho ~ N(0, 1/df) under the null, and the median of a half-normal is
0.6745 sigma, not 0.7979 sigma. With that functional:

    df = n - 1 = 156   ->  0.054002   ->  prints 0.0540   = script 77
    df = n     = 157   ->  0.053830   ->  prints 0.0538   = script 95

So the two scripts use the same formula and differ ONLY in whether the variance
of Spearman's rho is 1/(n-1) or 1/n. The standard result is **1/(n-1)**, so
**77's 0.0540 is the defensible value** and 95's 0.0538 is off by one in df.

WHAT ELSE STEP 1 SHOWED
-----------------------
`95_s1c_column_and_s3_split.py` line 363:

    ("analytic null expectation for |rho|", "0.0538"),

It is a **quoted string literal**. 95 does not compute the null at all -- it
transcribes a constant into a deposited summary table. That is the defect, not
the df choice: the number cannot track the data, and at 50 bp it would have
printed the same 0.0538 whatever n was. Rule 5 territory (a deposited number
must derive from the gate constants, not from a transcription).

WHAT THIS SCRIPT NOW DOES
-------------------------
  1. prints each script's computing/printing lines, and flags any observed
     value that appears in the source as a QUOTED LITERAL -- transcribed, not
     computed;
  2. recomputes under BOTH functionals (null mean and null median) x every
     n-definition x df offset x per-event aggregation;
  3. reports which combination explains which printed value, and names the
     functional that matches the statistic actually being compared.

GATES
-----
G1  is n constant across events? (It is not: 144..160, median 157.)
G2  the median |rho| reproduces as 0.0555 -- refuse otherwise; the comparator
    is moot if the compared value is wrong.
G3  the median |rho| exceeds every candidate null, so the direction of the
    sentence does not depend on the choice.
G4  the functional matches the statistic: comparing a MEDIAN to a null MEAN is
    an apples-to-oranges error and is reported as such.

USAGE
-----
  python 116_reconcile_analytic_null.py \
      --data SupplementaryData_lvef_all_events_rl125.tsv \
      --s77 77_lvef_s3.py --s95 95_s1c_column_and_s3_split.py
  python 116_reconcile_analytic_null.py --n 157
  python 116_reconcile_analytic_null.py --self-test
"""
import argparse
import math
import os
import re
import sys

try:
    import pandas as pd
except ImportError:
    pd = None

try:
    from statistics import NormalDist
    Z75 = NormalDist().inv_cdf(0.75)
except Exception:                                        # py<3.8
    Z75 = 0.6744897501960817

OBSERVED = {"77_lvef_s3.py": 0.0540, "95_s1c_column_and_s3_split.py": 0.0538}
RECORDED_MEDIAN_ABS_RHO = 0.0555

N_COLS = ["n", "n_used", "n_obs", "nobs", "N", "n_samples", "n_pairs",
          "n_lvef", "n_fit", "n_complete"]
RHO_COLS = ["rho", "rho_lvef", "spearman_rho", "rho_spearman", "rho_rank"]

SRC_PAT = re.compile(
    r"analytic|null[_ ]?exp|exp[_ ]?null|exp_med_rho|E_?abs[_ ]?rho|"
    r"2\s*/\s*(?:math\.|np\.)?pi|(?:math\.|np\.)?pi\s*\*|sqrt\s*\(\s*2|"
    r"0\.674|ppf\s*\(\s*0?\.75|inv_cdf\s*\(\s*0?\.75", re.I)


# ---------------------------------------------------------------- functionals
def null_mean(df):
    """E|rho| when rho ~ N(0, 1/df).  0.7979 * sigma."""
    return math.sqrt(2.0 / (math.pi * df)) if df and df > 0 else None


def null_median(df):
    """median |rho| when rho ~ N(0, 1/df).  z(0.75) * sigma = 0.6745 * sigma."""
    return Z75 / math.sqrt(df) if df and df > 0 else None


FUNCTIONALS = [
    ("null MEDIAN  0.6745/sqrt(df)", null_median, True),
    ("null mean    sqrt(2/(pi*df))", null_mean, False),
]
# the third element: does this functional match the statistic being compared
# (the median |rho|)?


def show_source(path, label):
    """Print the computing lines; return observed values found as literals."""
    print(f"--- {label}: {path}")
    transcribed = []
    if not os.path.exists(path):
        print("      not found on disk -- cannot show its definition")
        return transcribed
    txt = open(path, errors="replace").read()
    keep = []
    for i, ln in enumerate(txt.splitlines(), 1):
        s = ln.strip()
        if not s or s.startswith("#"):
            continue
        if SRC_PAT.search(s):
            keep.append((i, s))
    if not keep:
        print("      no line matching the analytic-null pattern.")
    for i, s in keep[:14]:
        print(f"      {i:>5}  {s[:104]}")
    if len(keep) > 14:
        print(f"      ... and {len(keep) - 14} more matching line(s)")

    # a printed value that appears QUOTED in the source is transcribed
    for _src, obs in OBSERVED.items():
        lit = f"{obs:.4f}"
        for m in re.finditer(r"""["']\s*""" + re.escape(lit) + r"""\s*["']""",
                             txt):
            line = txt[:m.start()].count("\n") + 1
            transcribed.append((obs, line))
    for obs, line in transcribed:
        print(f"      ⚠ {obs:.4f} appears as a QUOTED LITERAL at line {line}.")
        print("        This script TRANSCRIBES the value; it does not compute")
        print("        it. It cannot track the data and would print the same")
        print("        number at any read length or any n.")
    return transcribed


def pick(df, names):
    for c in names:
        if c in df.columns:
            return c
    low = {c.lower(): c for c in df.columns}
    for c in names:
        if c.lower() in low:
            return low[c.lower()]
    return None


def build_candidates(nspecs, series=None):
    """nspecs: list of (label, scalar n). series: optional per-event n."""
    out = []
    for fname, fn, matches in FUNCTIONALS:
        for label, nval in nspecs:
            for k in (0, 1, 2, 3):
                v = fn(nval - k)
                if v is not None:
                    out.append((f"{fname} | df = {label}-{k} ({nval:g}-{k})",
                                v, matches))
        if series is not None:
            for k in (0, 1, 2, 3):
                vals = [fn(x - k) for x in series if x - k > 0]
                if not vals:
                    continue
                vals.sort()
                mean = sum(vals) / len(vals)
                mid = (vals[len(vals) // 2] if len(vals) % 2
                       else 0.5 * (vals[len(vals) // 2 - 1]
                                   + vals[len(vals) // 2]))
                out.append((f"{fname} | mean over events, df = n_i-{k}",
                            mean, matches))
                out.append((f"{fname} | median over events, df = n_i-{k}",
                            mid, matches))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data",
                    default="SupplementaryData_lvef_all_events_rl125.tsv")
    ap.add_argument("--s77", default="77_lvef_s3.py")
    ap.add_argument("--s95", default="95_s1c_column_and_s3_split.py")
    ap.add_argument("--n", type=int, default=None)
    ap.add_argument("--self-test", action="store_true")
    a = ap.parse_args()
    if a.self_test:
        return self_test()

    self_test()
    print()
    print("=" * 78)
    print("116 — one quantity, two printed values")
    print("=" * 78)
    for k, v in OBSERVED.items():
        print(f"  {k:<34} analytic null for |rho|  {v:.4f}")
    print()

    print("=" * 78)
    print("STEP 1 — what each script actually computes")
    print("=" * 78)
    tr77 = show_source(a.s77, "77")
    print()
    tr95 = show_source(a.s95, "95")
    print()

    print("=" * 78)
    print("STEP 2 — recompute from the data")
    print("=" * 78)
    cands, med_rho, n_varies, nmed = [], None, None, None

    if a.n is not None:
        print(f"  --n {a.n} given; formulas only, no file read.")
        nmed = a.n
        cands = build_candidates([("n", a.n)])
    elif pd is None:
        print("  [fatal] pandas unavailable and no --n given", file=sys.stderr)
        return 2
    elif not os.path.exists(a.data):
        print(f"  [fatal] {a.data} not found.", file=sys.stderr)
        return 2
    else:
        df = pd.read_csv(a.data, sep="\t")
        print(f"  {a.data}: {len(df):,} rows")
        rcol = pick(df, RHO_COLS)
        if rcol:
            r = pd.to_numeric(df[rcol], errors="coerce").abs().dropna()
            med_rho = float(r.median())
            print(f"  rho column '{rcol}': median |rho| = {med_rho:.4f}"
                  f"   (recorded {RECORDED_MEDIAN_ABS_RHO:.4f})")
            if round(med_rho, 4) != RECORDED_MEDIAN_ABS_RHO:
                print("  ✗ G2 FAIL — the median |rho| does not reproduce.")
                return 3
            print("  ✓ G2 — median |rho| reproduces")
        else:
            print(f"  [warn] no rho column; has {list(df.columns)[:10]}")

        ncol = pick(df, N_COLS)
        if ncol is None:
            print(f"  [warn] no sample-size column; has "
                  f"{list(df.columns)[:10]}. Re-run with --n.")
        else:
            nn = pd.to_numeric(df[ncol], errors="coerce").dropna()
            nmin, nmed, nmax = int(nn.min()), float(nn.median()), int(nn.max())
            nmean = float(nn.mean())
            n_varies = nmin != nmax
            print(f"  n column '{ncol}': min {nmin}  median {nmed:g}  "
                  f"mean {nmean:.2f}  max {nmax}")
            if n_varies:
                print("  ⚠ G1 — n VARIES across events (an analytic null "
                      "quoted as a bare scalar stands for a distribution)")
            else:
                print(f"  ✓ G1 — n constant at {nmin}")
            cands = build_candidates(
                [("median n", nmed), ("mean n", nmean),
                 ("max n", nmax), ("min n", nmin)],
                series=[float(x) for x in nn.tolist()])

    if not cands:
        print("  [fatal] nothing to test.", file=sys.stderr)
        return 2

    print()
    print("=" * 78)
    print("STEP 3 — which candidate explains which printed value")
    print("=" * 78)
    explained = {}
    hits = []
    for label, v, matches in cands:
        tag = ""
        for src, obs in OBSERVED.items():
            if round(v, 4) == obs:
                tag = f"<- {src.split('_')[0]}"
                explained.setdefault(obs, []).append((label, matches))
                hits.append((label, v, obs, matches))
        if tag:
            print(f"  {label:<58}{v:>10.6f} {tag}")
    if not hits:
        print("  (no candidate reproduces either printed value)")
    print()
    print("  full grid suppressed; only matches shown. Re-run with --n to see"
          " a small grid.")

    print()
    print("=" * 78)
    print("VERDICT")
    print("=" * 78)
    for src, obs in sorted(OBSERVED.items(), key=lambda x: -x[1]):
        who = explained.get(obs)
        short = src.split("_")[0]
        if who:
            print(f"  {obs:.4f}  ({short})  explained by:")
            for label, matches in who:
                flag = "" if matches else "   ✗ WRONG FUNCTIONAL"
                print(f"       {label}{flag}")
        else:
            print(f"  {obs:.4f}  ({short})  NOT reproduced by any candidate.")
        for o, line in (tr77 if short == "77" else tr95):
            if o == obs:
                print(f"       ⚠ and it is a QUOTED LITERAL at line {line} — "
                      f"transcribed, not computed")
    print()

    # --- G4
    matched = [h for h in hits if h[3]]
    mismatched = [h for h in hits if not h[3]]
    if matched:
        print("  ✓ G4 — the statistic compared is a MEDIAN |rho|, and the")
        print("        matching candidate(s) use the null's MEDIAN. Correct")
        print("        functional.")
    if mismatched:
        print("  ⚠ G4 — some match came from the null's MEAN. Comparing a")
        print("        median to a mean is an apples-to-oranges error; do not")
        print("        adopt that candidate even though it reproduces.")

    if med_rho is not None:
        worst = max(OBSERVED.values())
        if med_rho <= worst:
            print("  ✗ G3 FAIL — median |rho| does not exceed the larger null.")
            return 3
        print(f"  ✓ G3 — median |rho| {med_rho:.4f} exceeds both; direction of")
        print(f"        the sentence is invariant. Excess over {worst:.4f}: "
              f"{100 * (med_rho / worst - 1):.1f}%.")

    print()
    print("-" * 78)
    print("WHAT TO DO")
    print("-" * 78)
    print("  Var(rho) under the null is 1/(n-1), not 1/n, so df = n-1 is the")
    print("  defensible choice and the OTHER value is off by one in df.")
    print("  Make the off-by-one script DERIVE the number from the same")
    print("  expression, then re-run both so they print the same value.")
    if tr95 or tr77:
        print()
        print("  ⚠ A transcribed literal cannot be 'made to agree' — it agrees")
        print("    by hand until someone changes n. Replace it with the")
        print("    computation, or delete the row from the deposited summary.")
    print()
    print("  ⚠ Independently of which df is used: a ~2.8% excess of the median")
    print("    |rho| over an ASYMPTOTIC null, across tens of thousands of")
    print("    correlated events, is not evidence of signal. S3's result is the")
    print("    agreed two clauses (1 at FDR < 0.05, 0 at |dPSI| >= 0.10). The")
    print("    median-vs-null line is descriptive context. To carry weight it")
    print("    needs a PERMUTATION null. Dropping it is defensible.")
    return 0 if explained else 4


def self_test():
    assert abs(Z75 - 0.67448975) < 1e-7, Z75
    # the resolution this script exists to prove
    assert round(null_median(156), 4) == 0.0540, null_median(156)
    assert round(null_median(157), 4) == 0.0538, null_median(157)
    # and that the MEAN functional explains neither -- v1's error
    assert round(null_mean(156), 4) == 0.0639, null_mean(156)
    assert round(null_mean(157), 4) == 0.0637, null_mean(157)
    assert not any(round(null_mean(d), 4) in OBSERVED.values()
                   for d in range(120, 200))
    # median functional, df=n-1 at n=157, must be the one flagged as matching
    cs = build_candidates([("n", 157)])
    m = [(l, v, ok) for l, v, ok in cs if round(v, 4) == 0.0540]
    assert m and all(ok for _l, _v, ok in m), m
    assert "MEDIAN" in m[0][0], m[0][0]

    if pd is not None:
        import tempfile
        d = tempfile.mkdtemp()
        p = os.path.join(d, "s.py")
        open(p, "w").write('ROWS = [\n'
                           '    ("analytic null expectation for |rho|", "0.0538"),\n'
                           ']\n')
        tr = show_source(p, "fixture")
        assert tr == [(0.0538, 2)], tr
        os.remove(p)

    print("[ok] self-test passed: null MEDIAN 0.6745/sqrt(df) gives 0.0540 at "
          "df=n-1 and 0.0538 at df=n for n=157 — the two printed values differ "
          "by n vs n-1, not by rounding; the null MEAN functional (v1's test) "
          "explains neither at any plausible df; quoted-literal transcription "
          "detected")
    return 0


if __name__ == "__main__":
    sys.exit(main())
