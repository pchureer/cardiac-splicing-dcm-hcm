#!/usr/bin/env python3
"""
109_login_recounts_rl125.py -- the track-1 items that need NO job: they only
read `*_rl125.all_events.tsv`, which already exists.

Fills these rows of NUMBERS_125bp_REDERIVATION.md:

  §1  full-coordinate equivalents        (were 215,491 tested / 1,074 SE)
  §1  effective n                        (were 165 / 161)
  §3  duplication factor                 (was 1.588x for SE)
  §8  dPSI tilt, exon level              (was 57.14% inclusion, median +0.00214
                                          on 62,971)
  §8  dPSI tilt, full coordinate         (was 44.21% negative, median +0.00143
                                          on 100,014)

Runs on the login node in seconds. Nothing here needs SLURM.

CONVENTIONS, both reported side by side because the project has been bitten by
confusing them:
  exon level      = collapse-first, key chr:exonStart:exonEnd (SE, RI) or
                    key[:5] (MXE, A3SS, A5SS), lowest-FDR row, THEN the rule
  full coordinate = every row, no collapse

CALIBRATION: the script reproduces the nine asserted 50 bp SE values from
`SE_v2filtered.all_events.tsv` before it will print any 125 bp number. If any
one of them misses, it refuses.

USAGE
-----
  python 109_login_recounts_rl125.py
  python 109_login_recounts_rl125.py --frozen-pattern '{C}_v2filtered.all_events.tsv' \
                                     --new-pattern '{C}_rl125.all_events.tsv'
  python 109_login_recounts_rl125.py --self-test
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

CLASSES = ["SE", "MXE", "A3SS", "A5SS", "RI"]
KEY = {"SE": ("coords", None), "RI": ("coords", None),
       "MXE": ("keytrunc", 5), "A3SS": ("keytrunc", 5), "A5SS": ("keytrunc", 5)}

FDR_CUT, DPSI_CUT = 0.05, 0.10

# asserted 50 bp SE values, from NUMBERS_CANONICAL.md
ASSERT_SE = {
    "rows": 100014, "sig_full": 1074,
    "collapsed": 62971, "sig_coll": 926,
    "dup": 1.588,
    "incl_pct_exon": 57.14, "median_exon": 0.00214,
    "neg_pct_full": 44.21, "median_full": 0.00143,
}


def build_key(df, spec):
    kind, n = spec
    if kind == "coords":
        return (df["chr"].astype(str) + ":" + df["exonStart"].astype(str)
                + ":" + df["exonEnd"].astype(str))
    if "key" not in df.columns:
        return None
    return df["key"].astype(str).str.split(":").str[:n].str.join(":")


def dpsi_col(df):
    if "dPSI_adj" in df.columns:
        return "dPSI_adj"
    c = [x for x in df.columns if x.startswith("dPSI_") and x != "dPSI_raw"]
    return sorted(c)[0] if c else None


def stats(path, cls):
    df = pd.read_csv(path, sep="\t")
    dc = dpsi_col(df)
    if dc is None or "FDR" not in df.columns:
        return None, "not an events table"
    out = {}
    # --- full coordinate: every row
    out["rows"] = len(df)
    a_full = df[dc].abs()
    out["sig_full"] = int(((df.FDR < FDR_CUT) & (a_full >= DPSI_CUT)).sum())
    v = df[dc].dropna()
    out["neg_pct_full"] = round(100.0 * float((v < 0).mean()), 2)
    out["median_full"] = round(float(v.median()), 5)
    # --- exon level: collapse-first
    kk = build_key(df, KEY[cls])
    if kk is None:
        return None, "no key column for keytrunc"
    d = df.copy()
    d["_ck"] = kk
    d = d.sort_values("FDR", kind="stable").drop_duplicates("_ck", keep="first")
    out["collapsed"] = len(d)
    a = d[dc].abs()
    out["sig_coll"] = int(((d.FDR < FDR_CUT) & (a >= DPSI_CUT)).sum())
    out["dup"] = round(out["rows"] / out["collapsed"], 3) if out["collapsed"] else float("nan")
    w = d[dc].dropna()
    out["incl_pct_exon"] = round(100.0 * float((w > 0).mean()), 2)
    out["neg_pct_exon"] = round(100.0 * float((w < 0).mean()), 2)
    out["median_exon"] = round(float(w.median()), 5)
    # --- effective n
    for c in ("n_g1", "n_g2"):
        if c in df.columns:
            out[f"{c}_min"] = int(df[c].min())
            out[f"{c}_max"] = int(df[c].max())
    return out, None


def self_test():
    """A fixture where collapse-first and full-coordinate genuinely differ."""
    rows = [
        # chr, exonStart, exonEnd, key, FDR, dPSI_adj
        ("chr1", 100, 200, "chr1:100:200:1:2", 0.001,  0.30),  # kept, sig
        ("chr1", 100, 200, "chr1:100:200:3:4", 0.010,  0.40),  # dropped, sig
        ("chr2", 300, 400, "chr2:300:400:1:2", 0.001, -0.20),  # kept, sig
        ("chr3", 500, 600, "chr3:500:600:1:2", 0.900,  0.90),  # not sig
    ]
    df = pd.DataFrame(rows, columns=["chr", "exonStart", "exonEnd", "key",
                                     "FDR", "dPSI_adj"])
    p = ".__109_selftest.tsv"
    df.to_csv(p, sep="\t", index=False)
    try:
        s, err = stats(p, "SE")
        assert err is None, err
        assert s["rows"] == 4 and s["sig_full"] == 3, s
        assert s["collapsed"] == 3 and s["sig_coll"] == 2, s
        assert s["dup"] == round(4 / 3, 3), s
        # exon level: kept dPSI are +0.30, -0.20, +0.90 -> 66.67% positive
        assert s["incl_pct_exon"] == 66.67, s["incl_pct_exon"]
        assert s["median_exon"] == 0.3, s["median_exon"]
        # full coordinate: +0.30 +0.40 -0.20 +0.90 -> 25.00% negative
        assert s["neg_pct_full"] == 25.00, s["neg_pct_full"]
        assert s["median_full"] == 0.35, s["median_full"]
        # key[:5] keeps all four rows distinct
        s5, _ = stats(p, "MXE")
        assert s5["collapsed"] == 4 and s5["sig_coll"] == 3, s5
    finally:
        os.remove(p)
    print("[ok] self-test passed: full-coordinate 4/3 vs exon-level 3/2, "
          "duplication 1.333x, tilt percentages and medians differ by "
          "convention as they must, key[:5] keeps all four")
    return 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--frozen-pattern", default="{C}_v2filtered.all_events.tsv")
    ap.add_argument("--new-pattern", default="{C}_rl125.all_events.tsv")
    ap.add_argument("--self-test", action="store_true")
    a = ap.parse_args()

    if a.self_test:
        return self_test()
    self_test()
    print()

    # ---------- calibration on 50 bp SE ----------
    fz = a.frozen_pattern.format(C="SE")
    print("=" * 78)
    print("CALIBRATION -- nine asserted 50 bp SE values")
    print("=" * 78)
    if not os.path.exists(fz):
        print(f"[fatal] {fz} not found. Pass --frozen-pattern.", file=sys.stderr)
        return 2
    s, err = stats(fz, "SE")
    if err:
        print(f"[fatal] {fz}: {err}", file=sys.stderr)
        return 2
    bad = []
    for k, want in ASSERT_SE.items():
        got = s.get(k)
        hit = (got == want)
        if not hit:
            bad.append((k, want, got))
        print(f"  {k:<16} asserted {str(want):>10}   got {str(got):>10}"
              f"   {'ok' if hit else 'MISMATCH'}")
    if "n_g1_max" in s:
        print(f"  effective n      asserted    165/161   got "
              f"{s['n_g1_max']:>6}/{s['n_g2_max']}"
              f"   {'ok' if (s['n_g1_max'], s['n_g2_max']) == (165, 161) else 'MISMATCH'}")
    if bad:
        print()
        print("[fatal] calibration failed on: "
              + ", ".join(k for k, _, _ in bad), file=sys.stderr)
        print("        Do not use any 125 bp number from this script until the "
              "50 bp values reproduce. Check --frozen-pattern first.",
              file=sys.stderr)
        return 3
    print("  -> all asserted 50 bp SE values reproduced")

    # ---------- the 125 bp table ----------
    print()
    print("=" * 78)
    print("FULL-COORDINATE AND EXON-LEVEL, 50 bp -> 125 bp")
    print("=" * 78)
    print(f"{'cls':<6}{'key':<10}"
          f"{'rows50':>9}{'rows125':>9}{'sigF50':>8}{'sigF125':>8}"
          f"{'coll50':>8}{'coll125':>9}{'sigC50':>8}{'sigC125':>8}{'dup125':>8}")
    tot = {"rows50": 0, "rows125": 0, "sigF50": 0, "sigF125": 0,
           "coll50": 0, "coll125": 0, "sigC50": 0, "sigC125": 0}
    se_new = None
    for c in CLASSES:
        f, n = a.frozen_pattern.format(C=c), a.new_pattern.format(C=c)
        if not os.path.exists(f) or not os.path.exists(n):
            miss = [p for p in (f, n) if not os.path.exists(p)]
            print(f"{c:<6}{'--':<10} missing: {miss}")
            continue
        s0, e0 = stats(f, c)
        s1, e1 = stats(n, c)
        if e0 or e1:
            print(f"{c:<6}{'--':<10} {e0 or e1}")
            continue
        if c == "SE":
            se_new = s1
        kn = "coords" if KEY[c][0] == "coords" else "key[:5]"
        print(f"{c:<6}{kn:<10}"
              f"{s0['rows']:>9}{s1['rows']:>9}"
              f"{s0['sig_full']:>8}{s1['sig_full']:>8}"
              f"{s0['collapsed']:>8}{s1['collapsed']:>9}"
              f"{s0['sig_coll']:>8}{s1['sig_coll']:>8}{s1['dup']:>8.3f}")
        tot["rows50"] += s0["rows"];      tot["rows125"] += s1["rows"]
        tot["sigF50"] += s0["sig_full"];  tot["sigF125"] += s1["sig_full"]
        tot["coll50"] += s0["collapsed"]; tot["coll125"] += s1["collapsed"]
        tot["sigC50"] += s0["sig_coll"];  tot["sigC125"] += s1["sig_coll"]
    print(f"{'TOTAL':<16}"
          f"{tot['rows50']:>9}{tot['rows125']:>9}"
          f"{tot['sigF50']:>8}{tot['sigF125']:>8}"
          f"{tot['coll50']:>8}{tot['coll125']:>9}"
          f"{tot['sigC50']:>8}{tot['sigC125']:>8}")
    print()
    print(f"  §1 full-coordinate, 50 bp : {tot['rows50']:,} tested, "
          f"{tot['sigF50']:,} significant across the classes shown")
    print(f"     asserted in NUMBERS_CANONICAL.md: 215,491 tested and 1,074 SE.")
    print(f"     These match only when all five classes are present above.")
    print(f"  §1 full-coordinate, 125 bp: {tot['rows125']:,} tested, "
          f"{tot['sigF125']:,} significant")

    if se_new:
        print()
        print("-" * 78)
        print("SE DETAIL -- the rows §1, §3 and §8 of the sheet need")
        print("-" * 78)
        print(f"  effective n                 50 bp 165/161   "
              f"125 bp {se_new.get('n_g1_max')}/{se_new.get('n_g2_max')}")
        print(f"  duplication factor          50 bp {ASSERT_SE['dup']:.3f}x    "
              f"125 bp {se_new['dup']:.3f}x")
        print(f"  tilt, exon level  % toward inclusion   "
              f"50 bp {ASSERT_SE['incl_pct_exon']:.2f}%   "
              f"125 bp {se_new['incl_pct_exon']:.2f}%")
        print(f"  tilt, exon level  median dPSI          "
              f"50 bp {ASSERT_SE['median_exon']:+.5f}   "
              f"125 bp {se_new['median_exon']:+.5f}")
        print(f"  tilt, full coord  % negative           "
              f"50 bp {ASSERT_SE['neg_pct_full']:.2f}%   "
              f"125 bp {se_new['neg_pct_full']:.2f}%")
        print(f"  tilt, full coord  median dPSI          "
              f"50 bp {ASSERT_SE['median_full']:+.5f}   "
              f"125 bp {se_new['median_full']:+.5f}")
        print()
        print("  ⚠ The exon-level tilt is a PERCENTAGE OF A CHANGED DENOMINATOR.")
        print("    57.14% was 'toward inclusion' on 62,971 collapsed events. If the")
        print("    125 bp percentage moves, say which denominator it is on.")

    print()
    print("Paste these into NUMBERS_125bp_REDERIVATION.md §1, §3 and §8.")
    print("Nothing here needed SLURM. The remaining track-1 items do.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
