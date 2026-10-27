#!/usr/bin/env python3
"""
113_recount_ladder_rl125.py -- collapse-first recount of the covariate ladder
and the composition tables at 125 bp. Fills §3 of NUMBERS_125bp_REDERIVATION.md.

Each row CALIBRATES first: the frozen 50 bp file must reproduce its asserted
count from NUMBERS_CANONICAL.md before the 125 bp number for that row is
printed. A row whose calibration fails prints WITHHELD.

Convention: collapse-first, key `chr:exonStart:exonEnd` for SE and RI,
`key[:5]` for MXE/A3SS/A5SS, lowest-FDR row, THEN FDR < 0.05 and
|dPSI_adj| >= 0.10.

WHAT THIS DOES *NOT* GIVE YOU
-----------------------------
The "retained" figures -- 590/926 = 63.7% (nine-marker) and 716/1,054 = 77.3%
(restricted) -- are NOT recounts. They come from the concordance step:
`06_sensitivity_concordance.py` via `${ET}_compositionCheck_rl125`, and
`71_cmfrac_score.py compare`. Read them from those outputs, not from here.

USAGE
-----
  python 113_recount_ladder_rl125.py
  python 113_recount_ladder_rl125.py --base-new SE_rl125.all_events.tsv
  python 113_recount_ladder_rl125.py --self-test
"""
import argparse
import os
import sys

import pandas as pd

FDR_CUT, DPSI_CUT = 0.05, 0.10
CLASSES = ["SE", "MXE", "A3SS", "A5SS", "RI"]
KEY = {"SE": ("coords", None), "RI": ("coords", None),
       "MXE": ("keytrunc", 5), "A3SS": ("keytrunc", 5), "A5SS": ("keytrunc", 5)}

# (row label, class, frozen file, new file, asserted 50 bp tested, asserted sig)
# asserted values from NUMBERS_CANONICAL.md: ladder 926/899/835/866/888,
# row 5 on 62,482; restricted SE 1,054; nine-marker MXE 286; restricted MXE 337.
LADDER = [
    ("1 base (age+sex+RIN)", "SE",
     "SE_v2filtered.all_events.tsv", "SE_v2filtered_rl125.all_events.tsv",
     62971, 926),
    ("2 + CM/FB composition", "SE",
     "SE_v2compadj.all_events.tsv", "SE_v2compadj_rl125.all_events.tsv",
     62971, 899),
    ("3 + composition+libsize+TIN", "SE",
     "SE_v2fulladj.all_events.tsv", "SE_v2fulladj_rl125.all_events.tsv",
     62971, 835),
    ("4 + sequencing pool", "SE",
     "SE_v2pooladj.all_events.tsv", "SE_v2pooladj_rl125.all_events.tsv",
     62971, 866),
    ("5 excluding donors <25 y", "SE",
     "SE_v2age25.all_events.tsv", "SE_v2age25_rl125.all_events.tsv",
     62482, 888),
    ("restricted (cmfrac)", "SE",
     "SE_v2cmfrac.all_events.tsv", "SE_v2cmfrac_rl125.all_events.tsv",
     62971, 1054),
    ("2 + CM/FB composition", "MXE",
     "MXE_v2compadj.all_events.tsv", "MXE_v2compadj_rl125.all_events.tsv",
     None, 286),
    ("3 + composition+libsize+TIN", "MXE",
     "MXE_v2fulladj.all_events.tsv", "MXE_v2fulladj_rl125.all_events.tsv",
     None, None),
    ("restricted (cmfrac)", "MXE",
     "MXE_v2cmfrac.all_events.tsv", "MXE_v2cmfrac_rl125.all_events.tsv",
     None, 337),
    # added 15 Oct: 04 and 75 have now written base and cmfrac for every class,
    # so the "exceeds the base in three of five only" claim is testable.
    ("1 base (age+sex+RIN)", "MXE",
     "MXE_v2filtered.all_events.tsv", "MXE_v2filtered_rl125.all_events.tsv",
     62762, 231),
    ("1 base (age+sex+RIN)", "A3SS",
     "A3SS_v2filtered.all_events.tsv", "A3SS_v2filtered_rl125.all_events.tsv",
     15292, 90),
    ("restricted (cmfrac)", "A3SS",
     "A3SS_v2cmfrac.all_events.tsv", "A3SS_v2cmfrac_rl125.all_events.tsv",
     None, 108),
    ("1 base (age+sex+RIN)", "A5SS",
     "A5SS_v2filtered.all_events.tsv", "A5SS_v2filtered_rl125.all_events.tsv",
     8988, 84),
    ("restricted (cmfrac)", "A5SS",
     "A5SS_v2cmfrac.all_events.tsv", "A5SS_v2cmfrac_rl125.all_events.tsv",
     None, 67),
    ("1 base (age+sex+RIN)", "RI",
     "RI_v2filtered.all_events.tsv", "RI_v2filtered_rl125.all_events.tsv",
     6102, 114),
    ("restricted (cmfrac)", "RI",
     "RI_v2cmfrac.all_events.tsv", "RI_v2cmfrac_rl125.all_events.tsv",
     None, 112),
    # added 15 Oct after job 407847 (script 11) wrote compadj/fulladj for the
    # three remaining classes. 50 bp nine-marker values 105 / 63 / 111 from
    # NUMBERS_CANONICAL.md; no asserted 50 bp fulladj value exists for these.
    ("2 + CM/FB composition", "A3SS",
     "A3SS_v2compadj.all_events.tsv", "A3SS_v2compadj_rl125.all_events.tsv",
     None, 105),
    ("3 + composition+libsize+TIN", "A3SS",
     "A3SS_v2fulladj.all_events.tsv", "A3SS_v2fulladj_rl125.all_events.tsv",
     None, None),
    ("2 + CM/FB composition", "A5SS",
     "A5SS_v2compadj.all_events.tsv", "A5SS_v2compadj_rl125.all_events.tsv",
     None, 63),
    ("3 + composition+libsize+TIN", "A5SS",
     "A5SS_v2fulladj.all_events.tsv", "A5SS_v2fulladj_rl125.all_events.tsv",
     None, None),
    ("2 + CM/FB composition", "RI",
     "RI_v2compadj.all_events.tsv", "RI_v2compadj_rl125.all_events.tsv",
     None, 111),
    ("3 + composition+libsize+TIN", "RI",
     "RI_v2fulladj.all_events.tsv", "RI_v2fulladj_rl125.all_events.tsv",
     None, None),
]


def build_key(df, spec):
    kind, n = spec
    if kind == "coords":
        return (df["chr"].astype(str) + ":" + df["exonStart"].astype(str)
                + ":" + df["exonEnd"].astype(str))
    if "key" not in df.columns:
        return None
    return df["key"].astype(str).str.split(":").str[:n].str.join(":")


def dcol(df):
    if "dPSI_adj" in df.columns:
        return "dPSI_adj"
    c = [x for x in df.columns if x.startswith("dPSI_") and x != "dPSI_raw"]
    return sorted(c)[0] if c else None


def count(path, cls):
    df = pd.read_csv(path, sep="\t")
    dc = dcol(df)
    if dc is None or "FDR" not in df.columns:
        return None, "not an events table"
    kk = build_key(df, KEY[cls])
    if kk is None:
        return None, "no key column"
    d = df.copy()
    d["_ck"] = kk
    d = d.sort_values("FDR", kind="stable").drop_duplicates("_ck", keep="first")
    a = d[dc].abs()
    sig = int(((d.FDR < FDR_CUT) & (a >= DPSI_CUT)).sum())
    return (len(df), len(d), sig), None


def self_test():
    rows = [("chr1", 100, 200, "chr1:100:200:1:2", 0.001, 0.30),
            ("chr1", 100, 200, "chr1:100:200:3:4", 0.010, 0.90),
            ("chr2", 300, 400, "chr2:300:400:1:2", 0.001, 0.02),
            ("chr3", 500, 600, "chr3:500:600:1:2", 0.900, 0.80)]
    df = pd.DataFrame(rows, columns=["chr", "exonStart", "exonEnd", "key",
                                     "FDR", "dPSI_adj"])
    p = ".__113_selftest.tsv"
    df.to_csv(p, sep="\t", index=False)
    try:
        got, err = count(p, "SE")
        assert err is None and got == (4, 3, 1), got
        got5, _ = count(p, "MXE")
        assert got5 == (4, 4, 2), got5
    finally:
        os.remove(p)
    print("[ok] self-test passed: SE coords 4 rows -> 3 collapsed -> 1 sig; "
          "MXE key[:5] -> 4 collapsed -> 2 sig")
    return 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base-new", default=None,
                    help="override the 125 bp base table (e.g. SE_rl125.all_events.tsv)")
    ap.add_argument("--self-test", action="store_true")
    a = ap.parse_args()
    if a.self_test:
        return self_test()
    self_test()
    print()

    print("=" * 86)
    print("COVARIATE LADDER AND COMPOSITION, collapse-first, 50 bp -> 125 bp")
    print("=" * 86)
    print(f"{'row':<30}{'cls':<5}{'t50':>8}{'s50':>6}"
          f"{'t125':>9}{'s125':>7}{'Δsig':>8}  calib")

    se_ladder = {}
    per_class = {}
    for label, cls, fz, nw, a_tested, a_sig in LADDER:
        if a.base_new and label.startswith("1 base") and cls == "SE":
            nw = a.base_new

        f_res = n_res = None
        calib = "--"
        if os.path.exists(fz):
            got, err = count(fz, cls)
            if err:
                calib = f"frozen: {err}"
            else:
                rows, tested, sig = got
                f_res = (tested, sig)
                ok_t = (a_tested is None) or (tested == a_tested)
                ok_s = (a_sig is None) or (sig == a_sig)
                if a_sig is None and a_tested is None:
                    calib = "no asserted value"
                elif ok_t and ok_s:
                    calib = "MATCHES"
                else:
                    calib = (f"MISMATCH got {tested}/{sig}, "
                             f"asserted {a_tested}/{a_sig}")
        else:
            calib = "frozen file missing"

        if os.path.exists(nw):
            got, err = count(nw, cls)
            n_res = None if err else (got[1], got[2])
        else:
            n_res = None

        t50 = f"{f_res[0]}" if f_res else "--"
        s50 = f"{f_res[1]}" if f_res else "--"
        if calib.startswith("MISMATCH"):
            t125, s125 = "WITHHELD", "-"
            dsig = ""
        elif n_res:
            t125, s125 = str(n_res[0]), str(n_res[1])
            dsig = (f"{n_res[1] - f_res[1]:+d}" if f_res else "")
            if cls == "SE":
                se_ladder[label] = n_res[1]
            per_class.setdefault(cls, {})[label] = (n_res[1], f_res[1] if f_res else None)
        else:
            t125, s125 = "pending", "-"
            dsig = ""
        print(f"{label:<30}{cls:<5}{t50:>8}{s50:>6}{t125:>9}{s125:>7}"
              f"{dsig:>8}  {calib}")

    # --- the pattern claims, which can flip while the numbers look fine
    print()
    print("-" * 86)
    print("PATTERN CLAIMS — re-check these, do NOT carry them over")
    print("-" * 86)
    base = se_ladder.get("1 base (age+sex+RIN)")
    rows5 = [se_ladder.get(k) for k in
             ("1 base (age+sex+RIN)", "2 + CM/FB composition",
              "3 + composition+libsize+TIN", "4 + sequencing pool",
              "5 excluding donors <25 y")]
    if all(r is not None for r in rows5):
        spread = max(rows5) - min(rows5)
        print(f"  ladder at 125 bp: {' / '.join(str(r) for r in rows5)}")
        print(f"  spread {spread} = {100.0 * spread / max(rows5):.1f}%"
              f"   (50 bp was 926/899/835/866/888, spread 91 = 9.8%)")
    else:
        print("  ladder incomplete — run 04 and 08 first")
    restr = se_ladder.get("restricted (cmfrac)")
    if base is not None and restr is not None:
        print(f"  SE restricted {restr} vs base {base}: "
              f"restricted {'EXCEEDS' if restr > base else 'is BELOW'} the base")
        print("  50 bp: exceeded in THREE of five classes only "
              "(A5SS 67<84, RI 112<114). Re-check per class.")
    else:
        print("  restricted count not available — run 72 and 75")
    # --- the per-class claim, now testable
    print()
    print('  CLAIM: "the count EXCEEDS the base in three of five classes only"')
    print(f"    {'class':<7}{'base50':>8}{'restr50':>9}{'50bp':>9}"
          f"{'base125':>10}{'restr125':>10}{'125bp':>9}")
    n_exceed_50 = n_exceed_125 = n_testable = 0
    for c in CLASSES:
        d = per_class.get(c, {})
        b = d.get("1 base (age+sex+RIN)")
        r = d.get("restricted (cmfrac)")
        if not b or not r:
            print(f"    {c:<7}{'--':>8}{'--':>9}{'--':>7}{'pending':>9}")
            continue
        b125, b50 = b
        r125, r50 = r
        v50 = ("EXCEEDS" if (b50 and r50 and r50 > b50) else
               "below" if (b50 and r50) else "?")
        v125 = "EXCEEDS" if r125 > b125 else "below"
        print(f"    {c:<7}{str(b50):>8}{str(r50):>9}{v50:>9}"
              f"{b125:>10}{r125:>10}{v125:>9}")
        n_testable += 1
        if b50 and r50 and r50 > b50:
            n_exceed_50 += 1
        if r125 > b125:
            n_exceed_125 += 1
    if n_testable:
        print(f"    -> 50 bp: {n_exceed_50} of {n_testable} exceed; "
              f"125 bp: {n_exceed_125} of {n_testable} exceed")
        if n_exceed_50 != n_exceed_125:
            print("    ⚠ THE COUNT CHANGED. The sentence must be rewritten, "
                  "not renumbered.")
        else:
            print("    the pattern holds; only the numbers move")

    print()
    print('  ⚠ "retention higher under the restricted score in ALL FIVE classes"')
    print("    and the 590/926 = 63.7% and 716/1,054 = 77.3% figures are NOT")
    print("    recounts. Read them from ${ET}_compositionCheck_rl125 and from")
    print("    71_cmfrac_score.py compare. This script cannot give them.")
    print()
    print("Paste the 125 bp column into NUMBERS_125bp_REDERIVATION.md §3.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
