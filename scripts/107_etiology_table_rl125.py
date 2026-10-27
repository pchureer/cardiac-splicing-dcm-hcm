#!/usr/bin/env python3
"""
107_etiology_table_rl125.py -- fill section 2 of NUMBERS_125bp_REDERIVATION.md:
the etiology table, all three contrasts, all four columns, at 125 bp.

The table in NUMBERS_CANONICAL.md is SE ONLY and has four columns:

  contrast      tested   FDR only   +>=0.05   +>=0.10
  DCM vs HCM    45,380        3          1         0
  DCM vs NF     62,971    5,301      2,592       926
  HCM vs NF     48,024    1,355        909       468

`106_recount_hcm_rl125.py` gives tested and +>=0.10. It does not give the two
middle columns, and those are in the manuscript. This script gives all four,
for all three contrasts, and will not print a 125 bp row for a contrast whose
50 bp row it cannot reproduce from the frozen file first.

Convention: collapse-first, key = chr:exonStart:exonEnd (= key[:3]) for SE,
keep the lowest-FDR row, THEN threshold. SE only -- this table is SE only.

The dPSI column differs by provenance and is auto-detected:
  phase2_diffsplice_v2.py  -> dPSI_adj
  07_dcm_vs_hcm_direct.py  -> dPSI_DCMvsHCM   (so the DCM vs HCM files are NOT
                              "not a phase2 output"; they are a different
                              script's output with the same meaning)

USAGE
-----
  # if a filename below is wrong, find the right one first:
  python 107_etiology_table_rl125.py --find

  python 107_etiology_table_rl125.py
  python 107_etiology_table_rl125.py --new-dcmvsnf SE_rl125.all_events.tsv
  python 107_etiology_table_rl125.py --self-test
"""
import argparse
import glob
import os
import re
import sys

import pandas as pd

FDR_CUT = 0.05

# contrast -> (frozen 50 bp file, rl125 file, asserted 50 bp row)
# the asserted row is (tested, fdr_only, ge05, ge10) from NUMBERS_CANONICAL.md
CONTRASTS = [
    ("DCM vs HCM", "SE_DCMvsHCM.all_events.tsv",
     "SE_DCMvsHCM_rl125.all_events.tsv", (45380, 3, 1, 0)),
    ("DCM vs NF", "SE_v2filtered.all_events.tsv",
     "SE_rl125.all_events.tsv", (62971, 5301, 2592, 926)),
    ("HCM vs NF", "SE_HCM_v2.all_events.tsv",
     "SE_HCMvsNF_rl125.all_events.tsv", (48024, 1355, 909, 468)),
]


def dpsi_col(df):
    """Find the adjusted-effect column whatever script wrote the file."""
    if "dPSI_adj" in df.columns:
        return "dPSI_adj"
    cands = [c for c in df.columns
             if re.match(r"^dPSI_", c) and c not in ("dPSI_raw",)]
    if len(cands) == 1:
        return cands[0]
    if cands:
        return sorted(cands)[0]
    return None


def counts(path):
    """-> (tested, fdr_only, ge05, ge10, dpsi_column) or (None, reason)."""
    df = pd.read_csv(path, sep="\t")
    need = {"chr", "exonStart", "exonEnd", "FDR"}
    miss = need - set(df.columns)
    if miss:
        return None, f"missing {sorted(miss)}"
    dc = dpsi_col(df)
    if dc is None:
        return None, f"no dPSI_* column; has {list(df.columns)[:10]}"
    d = df.copy()
    d["_ck"] = (d["chr"].astype(str) + ":" + d["exonStart"].astype(str)
                + ":" + d["exonEnd"].astype(str))
    d = d.sort_values("FDR", kind="stable").drop_duplicates("_ck", keep="first")
    s = d[d.FDR < FDR_CUT]
    a = s[dc].abs()
    return (len(d), len(s), int((a >= 0.05).sum()), int((a >= 0.10).sum()), dc), None


def self_test():
    rows = [
        # kept row has lowest FDR in its group
        ("chr1", 100, 200, 0.001, 0.30),   # sig, >=0.10
        ("chr1", 100, 200, 0.010, 0.90),   # dropped by collapse
        ("chr2", 300, 400, 0.001, 0.07),   # FDR only, >=0.05, not >=0.10
        ("chr3", 500, 600, 0.001, 0.01),   # FDR only
        ("chr4", 700, 800, 0.900, 0.80),   # not FDR
    ]
    df = pd.DataFrame(rows, columns=["chr", "exonStart", "exonEnd",
                                     "FDR", "dPSI_adj"])
    p = ".__107_selftest.tsv"
    df.to_csv(p, sep="\t", index=False)
    try:
        got, err = counts(p)
        assert err is None, err
        tested, fo, g5, g10, dc = got
        assert (tested, fo, g5, g10) == (4, 3, 2, 1), (tested, fo, g5, g10)
        assert dc == "dPSI_adj"
        # same file, the DCM-vs-HCM column name
        df.rename(columns={"dPSI_adj": "dPSI_DCMvsHCM"}).to_csv(p, sep="\t", index=False)
        got2, err2 = counts(p)
        assert err2 is None, err2
        assert got2[:4] == (4, 3, 2, 1), got2
        assert got2[4] == "dPSI_DCMvsHCM"
    finally:
        os.remove(p)
    print("[ok] self-test passed: nested thresholds correct (4/3/2/1), "
          "and dPSI_DCMvsHCM is detected as the effect column")
    return 0


def do_find():
    files = sorted(glob.glob("SE*all_events.tsv"))
    if not files:
        print("[fatal] no SE*all_events.tsv here", file=sys.stderr)
        return 2
    print(f"{'file':<46}{'tested':>8}{'FDRonly':>9}{'>=.05':>7}{'>=.10':>7}  dPSI col")
    for f in files:
        got, err = counts(f)
        if err:
            print(f"{f:<46}  skipped: {err}")
            continue
        t, fo, g5, g10, dc = got
        hit = ""
        for name, fz, nw, rec in CONTRASTS:
            if (t, fo, g5, g10) == rec:
                hit = f"  <== the 50 bp {name} row"
        print(f"{f:<46}{t:>8}{fo:>9}{g5:>7}{g10:>7}  {dc}{hit}")
    return 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--self-test", action="store_true")
    ap.add_argument("--find", action="store_true",
                    help="count every SE*all_events.tsv and flag which files "
                         "reproduce the three asserted 50 bp rows")
    for name, fz, nw, rec in CONTRASTS:
        slug = name.replace(" vs ", "vs").replace(" ", "").lower()
        ap.add_argument(f"--frozen-{slug}", default=fz)
        ap.add_argument(f"--new-{slug}", default=nw)
    a = ap.parse_args()

    if a.self_test:
        return self_test()
    if a.find:
        return do_find()
    self_test()
    print()

    print("=" * 78)
    print("CALIBRATION -- reproduce each asserted 50 bp row before using its 125 bp row")
    print("=" * 78)
    ok, out = {}, {}
    for name, fz_d, nw_d, rec in CONTRASTS:
        slug = name.replace(" vs ", "vs").replace(" ", "").lower()
        fz = getattr(a, f"frozen_{slug}")
        if not os.path.exists(fz):
            print(f"{name:<12} frozen file missing: {fz}")
            ok[name] = False
            continue
        got, err = counts(fz)
        if err:
            print(f"{name:<12} {fz}: {err}")
            ok[name] = False
            continue
        t, fo, g5, g10, dc = got
        match = (t, fo, g5, g10) == rec
        ok[name] = match
        flag = "[MATCHES]" if match else "[DOES NOT MATCH]"
        print(f"{name:<12} {fz}")
        print(f"{'':<12}   got      {t:>7} {fo:>7} {g5:>7} {g10:>7}   ({dc})  {flag}")
        if not match:
            print(f"{'':<12}   asserted {rec[0]:>7} {rec[1]:>7} {rec[2]:>7} {rec[3]:>7}")

    print()
    print("=" * 78)
    print("ETIOLOGY TABLE -- SE only, collapse-first, key chr:exonStart:exonEnd")
    print("=" * 78)
    print(f"{'contrast':<12}{'read len':>9}{'tested':>9}{'FDR only':>10}"
          f"{'+>=0.05':>9}{'+>=0.10':>9}")
    for name, fz_d, nw_d, rec in CONTRASTS:
        slug = name.replace(" vs ", "vs").replace(" ", "").lower()
        print(f"{name:<12}{'50 bp':>9}{rec[0]:>9}{rec[1]:>10}{rec[2]:>9}{rec[3]:>9}")
        nw = getattr(a, f"new_{slug}")
        if not ok.get(name):
            print(f"{'':<12}{'125 bp':>9}{'WITHHELD -- 50 bp row not reproduced':>46}")
            continue
        if not os.path.exists(nw):
            cand = [f for f in sorted(glob.glob("SE*all_events.tsv"))
                    if "rl125" in f]
            print(f"{'':<12}{'125 bp':>9}   file not found: {nw}")
            if cand:
                print(f"{'':<12}            rl125 files present: {cand}")
                print(f"{'':<12}            pass the right one with --new-{slug}")
            continue
        got, err = counts(nw)
        if err:
            print(f"{'':<12}{'125 bp':>9}   {nw}: {err}")
            continue
        t, fo, g5, g10, dc = got
        print(f"{'':<12}{'125 bp':>9}{t:>9}{fo:>10}{g5:>9}{g10:>9}")
        out[name] = (rec, (t, fo, g5, g10))
        print(f"{'':<12}{'change':>9}"
              f"{t - rec[0]:>+9}{fo - rec[1]:>+10}{g5 - rec[2]:>+9}{g10 - rec[3]:>+9}")

    # --- the two ratios the manuscript quotes, recomputed
    print()
    print("-" * 78)
    print("RATIOS AT FDR ONLY (the manuscript quotes 452-fold and 1,767-fold)")
    print("-" * 78)
    if "DCM vs HCM" in out and "DCM vs NF" in out and "HCM vs NF" in out:
        dh = out["DCM vs HCM"][1][1]
        dn = out["DCM vs NF"][1][1]
        hn = out["HCM vs NF"][1][1]
        for lab, num, den in (("DCM vs NF : DCM vs HCM", dn, dh),
                              ("HCM vs NF : DCM vs HCM", hn, dh)):
            if den == 0:
                print(f"  {lab:<26} {num} vs 0 -- ratio UNDEFINED at 125 bp. "
                      f"Report the counts, not a fold-change.")
            else:
                print(f"  {lab:<26} {num}/{den} = {num/den:,.0f}-fold")
    else:
        print("  not all three 125 bp rows available")

    print()
    print("Paste the 125 bp rows into section 2 of NUMBERS_125bp_REDERIVATION.md. "
          "Do not move any gate until that page's column is full.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
