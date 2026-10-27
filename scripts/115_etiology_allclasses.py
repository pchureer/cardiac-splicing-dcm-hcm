#!/usr/bin/env python3
"""
115_etiology_allclasses.py -- the DCM vs HCM null across ALL FIVE event types,
at both read lengths, with the correct collapse key per class.

WHY THIS EXISTS
---------------
Job 407847's summary printed a `DCMvsHCM_sig` column:

    SE 0 · MXE 0 · A3SS 1 · A5SS 1 · RI 0

**The null is not zero in A3SS and A5SS at 125 bp.** One event each, at
FDR < 0.05 and |dPSI| >= 0.10. So a sentence of the form "no significant
differences between DCM and HCM in any event type" would be FALSE as written.

Before that can be reported either way, the 50 bp equivalent for the same five
classes under the same convention is needed — and it has never been measured
with the calibrated per-class keys. `107_etiology_table_rl125.py` is SE only.
`NUMBERS_CANONICAL.md` carries only the SE row (45,380 / 3 / 1 / 0) plus a
RETIRED full-coordinate figure ("6 of 147,355") that used no collapse at all.

This script measures all five classes at both read lengths, four columns each,
with `chr:exonStart:exonEnd` for SE and RI and `key[:5]` for MXE/A3SS/A5SS.
The SE 50 bp row is gated against 45,380 / 3 / 1 / 0; if it does not reproduce,
nothing is printed.

The effect column is auto-detected: `07_dcm_vs_hcm_direct.py` writes
`dPSI_DCMvsHCM`, not `dPSI_adj`.

USAGE
-----
  python 115_etiology_allclasses.py
  python 115_etiology_allclasses.py --find          # inventory *DCMvsHCM* files
  python 115_etiology_allclasses.py --self-test
"""
import argparse
import glob
import os
import re
import sys

import pandas as pd

FDR_CUT = 0.05
CLASSES = ["SE", "MXE", "A3SS", "A5SS", "RI"]
KEY = {"SE": ("coords", None), "RI": ("coords", None),
       "MXE": ("keytrunc", 5), "A3SS": ("keytrunc", 5), "A5SS": ("keytrunc", 5)}

# the ONLY asserted 50 bp row, from NUMBERS_CANONICAL.md
ASSERT_SE_50 = (45380, 3, 1, 0)


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
    c = [x for x in df.columns if re.match(r"^dPSI_", x) and x != "dPSI_raw"]
    return sorted(c)[0] if c else None


def counts(path, cls):
    """-> (tested, fdr_only, ge05, ge10, dpsi_column) or (None, reason)."""
    df = pd.read_csv(path, sep="\t")
    if "FDR" not in df.columns:
        return None, "no FDR column"
    dc = dpsi_col(df)
    if dc is None:
        return None, f"no dPSI_* column; has {list(df.columns)[:8]}"
    kk = build_key(df, KEY[cls])
    if kk is None:
        return None, "no key column for keytrunc"
    d = df.copy()
    d["_ck"] = kk
    d = d.sort_values("FDR", kind="stable").drop_duplicates("_ck", keep="first")
    s = d[d.FDR < FDR_CUT]
    a = s[dc].abs()
    return (len(d), len(s), int((a >= 0.05).sum()),
            int((a >= 0.10).sum()), dc), None


def self_test():
    rows = [
        # kept row lowest FDR in its group
        ("chr1", 100, 200, "chr1:100:200:1:2", 0.001, 0.30),  # sig >=0.10
        ("chr1", 100, 200, "chr1:100:200:3:4", 0.010, 0.90),  # dropped by coords
        ("chr2", 300, 400, "chr2:300:400:1:2", 0.001, 0.07),  # >=0.05 only
        ("chr3", 500, 600, "chr3:500:600:1:2", 0.001, 0.01),  # FDR only
        ("chr4", 700, 800, "chr4:700:800:1:2", 0.900, 0.80),  # not sig
    ]
    df = pd.DataFrame(rows, columns=["chr", "exonStart", "exonEnd", "key",
                                     "FDR", "dPSI_DCMvsHCM"])
    p = ".__115_selftest.tsv"
    df.to_csv(p, sep="\t", index=False)
    try:
        got, err = counts(p, "SE")
        assert err is None, err
        assert got[:4] == (4, 3, 2, 1), got
        assert got[4] == "dPSI_DCMvsHCM", got
        # key[:5] keeps both chr1 rows -> one more tested, one more sig
        got5, err5 = counts(p, "A3SS")
        assert err5 is None, err5
        assert got5[:4] == (5, 4, 3, 2), got5
    finally:
        os.remove(p)
    print("[ok] self-test passed: SE coords 4/3/2/1, A3SS key[:5] 5/4/3/2, "
          "dPSI_DCMvsHCM detected as the effect column")
    return 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--find", action="store_true")
    ap.add_argument("--self-test", action="store_true")
    ap.add_argument("--old-pattern", default="{C}_DCMvsHCM.all_events.tsv")
    ap.add_argument("--new-pattern", default="{C}_DCMvsHCM_rl125.all_events.tsv")
    a = ap.parse_args()
    if a.self_test:
        return self_test()

    if a.find:
        for f in sorted(glob.glob("*DCMvsHCM*.all_events.tsv")):
            cls = f.split("_")[0]
            if cls not in KEY:
                print(f"{f:<46}  class not recognised from the name")
                continue
            got, err = counts(f, cls)
            if err:
                print(f"{f:<46}  {err}")
                continue
            t, fo, g5, g10, dc = got
            print(f"{f:<46}{t:>8}{fo:>7}{g5:>6}{g10:>6}   {dc}")
        return 0

    self_test()
    print()

    # --- gate on the one asserted row
    se_old = a.old_pattern.format(C="SE")
    print("=" * 76)
    print("CALIBRATION — the SE 50 bp row must reproduce 45,380 / 3 / 1 / 0")
    print("=" * 76)
    if not os.path.exists(se_old):
        print(f"[fatal] {se_old} not found", file=sys.stderr)
        return 2
    got, err = counts(se_old, "SE")
    if err:
        print(f"[fatal] {se_old}: {err}", file=sys.stderr)
        return 2
    if got[:4] != ASSERT_SE_50:
        print(f"  got {got[:4]}, asserted {ASSERT_SE_50}   MISMATCH")
        print("[fatal] refusing to print any row until the asserted SE row "
              "reproduces.", file=sys.stderr)
        return 3
    print(f"  got {got[:4]}   MATCHES")
    print()

    print("=" * 76)
    print("DCM vs HCM — ALL FIVE EVENT TYPES, collapse-first")
    print("=" * 76)
    print(f"{'class':<6}{'key':<10}{'read':>7}{'tested':>9}{'FDRonly':>9}"
          f"{'+>=.05':>8}{'+>=.10':>8}")
    tot = {}
    for c in CLASSES:
        kn = "coords" if KEY[c][0] == "coords" else "key[:5]"
        for tag, pat in (("50 bp", a.old_pattern), ("125 bp", a.new_pattern)):
            f = pat.format(C=c)
            if not os.path.exists(f):
                print(f"{c:<6}{kn:<10}{tag:>7}   missing: {f}")
                continue
            got, err = counts(f, c)
            if err:
                print(f"{c:<6}{kn:<10}{tag:>7}   {err}")
                continue
            t, fo, g5, g10, _dc = got
            print(f"{c:<6}{kn:<10}{tag:>7}{t:>9}{fo:>9}{g5:>8}{g10:>8}")
            d = tot.setdefault(tag, [0, 0, 0, 0])
            for i, v in enumerate((t, fo, g5, g10)):
                d[i] += v
    print()
    for tag in ("50 bp", "125 bp"):
        if tag in tot:
            t, fo, g5, g10 = tot[tag]
            print(f"{'TOTAL':<16}{tag:>7}{t:>9}{fo:>9}{g5:>8}{g10:>8}")

    print()
    print("-" * 76)
    print("HOW TO WRITE THIS")
    print("-" * 76)
    if "125 bp" in tot and "50 bp" in tot:
        g10_50, g10_125 = tot["50 bp"][3], tot["125 bp"][3]
        t_125 = tot["125 bp"][0]
        if g10_125 == 0:
            print("  All five classes give ZERO at FDR<0.05 and |dPSI|>=0.10.")
            print('  "No event of any class differed between DCM and HCM" is')
            print("  then literally true and is the strongest honest form.")
        else:
            print(f"  ⚠ {g10_125} event(s) across the five classes pass both")
            print("    thresholds at 125 bp. A sentence saying NO events differ")
            print("    would be FALSE. Write the count:")
            print(f"      '{g10_125} of {t_125:,} tested events' ")
            print(f"    and give the 50 bp figure ({g10_50}) beside it.")
            print("  Name the classes. One event in a class with a few thousand")
            print("  tested is a different statement from one in SE.")
    print()
    print("  ⚠ Do NOT mix this with the RETIRED '6 of 147,355', which was")
    print("    full-coordinate with no collapse at all.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
