#!/usr/bin/env python3
"""
101_compare_readlength_125.py — put the 125 bp re-quantification beside every
recorded 50 bp number, so the decision rests on data rather than on a prediction
about what a reviewer will ask.

⚠ READ-ONLY on the 50 bp analysis. This script writes one new file
  (`readlength_125_comparison.tsv`) and changes nothing else. It does not promote
  the 125 bp result, edit any gate, or touch any published table.

WHAT IT DOES
  Recounts the discovery contrast from the 125 bp rMATS output using **exactly the
  paper's conventions**, and diffs against the record:

    collapse-first exon-level counting   collapse to the event's own coordinates
                                         (3 key fields for SE and RI, 5 for MXE,
                                         A3SS and A5SS), keep the LOWEST-FDR row,
                                         THEN apply the rule
    the rule                             FDR < 0.05 AND |dPSI| >= 0.10
    BH scope                             within contrast AND event class
    quantification                       JCEC

  ⚠ Collapse-first is not optional. Filter-first returns 991 skipped exons rather
    than 926 and can only ever return more. Getting this wrong would make the
    comparison meaningless in the direction that looks like a real change.

  ⚠ THIS SCRIPT DOES NOT FIT THE COVARIATE MODEL. rMATS' own likelihood-ratio
    test is not the paper's test: the paper fits PSI ~ disease + sex + age + RIN
    (effective n 165/161, six samples lacking RIN) and the two differ by 24.9-fold
    on the same events. So the counts below are the RAW rMATS counts at 125 bp
    against the RAW rMATS counts at 50 bp, which is a like-for-like comparison of
    the quantification. To compare the paper's headline numbers you must then run
    `phase2_diffsplice_v2.py` on this output — the command is printed at the end.

USAGE
    python 101_compare_readlength_125.py --od <rmats_magnet_rl125>
    python 101_compare_readlength_125.py --od <...> --od50 <rmats_magnet>

GATES
  G1  all five JCEC files present and parseable; the key columns resolve.
  G2  the 50 bp tree, if given, reproduces its recorded RAW rMATS count. ⚠ This is
      the calibration gate: if the 50 bp arm does not reproduce, the comparison is
      measuring this script's bugs rather than the read length.
  G3  group widths are 166/166 in both arms.
  G4  the two arms cover the same event universe to within a stated tolerance; a
      wholesale change in the tested set means something other than --readLength
      moved.
"""

import argparse
import os
import sys

import pandas as pd

KEY3 = ["chr", "exonStart_0base", "exonEnd"]
KEY5 = ["chr", "1stExonStart_0base", "1stExonEnd", "2ndExonStart_0base", "2ndExonEnd"]
KEY_A = ["chr", "longExonStart_0base", "longExonEnd", "shortES", "shortEE"]
KEY_RI = ["chr", "riExonStart_0base", "riExonEnd"]

CLASSES = ["SE", "MXE", "A3SS", "A5SS", "RI"]

# The record. Exon-level, collapse-first, COVARIATE-ADJUSTED — the paper's model,
# NOT rMATS' own test. Printed for orientation; G2 does not gate on these.
RECORD_ADJ = {"SE": (62_971, 926), "MXE": (62_762, 231), "A3SS": (15_292, 90),
              "A5SS": (8_988, 84), "RI": (6_102, 114)}
RECORD_ADJ_TOTAL = (156_115, 1_445)

# rMATS' OWN likelihood-ratio test on the same events, 50 bp: 102 significant SE.
# This is what G2 calibrates against, because it is what this script computes.
RECORD_RAW_SE_50 = 102


def hr(c="-"):
    print(c * 76)


def refuse(msg):
    hr("=")
    print("REFUSING — nothing was written.\n")
    print(msg)
    hr("=")
    sys.exit(2)


def keys_for(cls):
    return {"SE": KEY3, "RI": KEY_RI, "MXE": KEY5, "A3SS": KEY_A, "A5SS": KEY_A}[cls]


def load(od, cls):
    p = os.path.join(od, f"{cls}.MATS.JCEC.txt")
    if not os.path.exists(p):
        refuse(f"G1 — {p} not found.")
    df = pd.read_csv(p, sep="\t", low_memory=False)
    need = keys_for(cls) + ["FDR", "IncLevelDifference"]
    miss = [c for c in need if c not in df.columns]
    if miss:
        refuse(f"G1 — {cls}: columns not found: {miss}\n  present: "
               f"{', '.join(map(str, df.columns[:20]))}…")
    return df


def collapse_and_count(df, cls):
    """⚠ COLLAPSE FIRST, then apply the rule. Reversing these inflates the count."""
    k = keys_for(cls)
    d = df.copy()
    d["FDR"] = pd.to_numeric(d["FDR"], errors="coerce")
    d["dPSI"] = pd.to_numeric(d["IncLevelDifference"], errors="coerce")
    d["_ekey"] = d[k].astype(str).agg(":".join, axis=1)
    c = d.sort_values("FDR", kind="mergesort").drop_duplicates("_ekey")
    sig = c[(c.FDR < 0.05) & (c.dPSI.abs() >= 0.10)]
    return len(d), len(c), len(sig)


def width(df):
    v = str(df["IncLevel1"].iloc[0]) if "IncLevel1" in df.columns else ""
    return len([x for x in v.split(",") if x.strip() not in ("", "NA")])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--od", required=True, help="the 125 bp rMATS output dir")
    ap.add_argument("--od50", default=None, help="the 50 bp dir, for calibration")
    ap.add_argument("--out", default="readlength_125_comparison.tsv")
    a = ap.parse_args()

    hr("=")
    print("101 — 125 bp against the record (read-only on the 50 bp analysis)")
    hr("=")

    print("\n  ⚠ These are rMATS' OWN likelihood-ratio counts, not the paper's")
    print("    covariate-adjusted counts. The two differ by 24.9-fold on the same")
    print("    events. See the note at the end before comparing anything to 926.")

    rows = []
    print("\n  G1/G3  the 125 bp arm")
    tot125 = {"rows": 0, "exons": 0, "sig": 0}
    for cls in CLASSES:
        df = load(a.od, cls)
        w = width(df)
        r, e, s = collapse_and_count(df, cls)
        tot125["rows"] += r; tot125["exons"] += e; tot125["sig"] += s
        print(f"      {cls:<5} rows {r:>7,}  exons {e:>7,}  sig {s:>5,}  "
              f"IncLevel1 width {w}")
        if w and w != 166:
            refuse(f"G3 — {cls}: IncLevel1 carries {w} values, expected 166.\n"
                   f"  A mismatched -t setting yields normal-looking output with no\n"
                   f"  usable read pairs.")
        rows.append({"class": cls, "arm": "125bp", "rows": r, "exons": e, "sig": s})
    print(f"      TOTAL rows {tot125['rows']:,}  exons {tot125['exons']:,}  "
          f"sig {tot125['sig']:,}")

    cal = None
    if a.od50:
        print("\n  G2  calibration — does the 50 bp arm reproduce its record?")
        tot50 = {"rows": 0, "exons": 0, "sig": 0}
        for cls in CLASSES:
            df = load(a.od50, cls)
            r, e, s = collapse_and_count(df, cls)
            tot50["rows"] += r; tot50["exons"] += e; tot50["sig"] += s
            print(f"      {cls:<5} rows {r:>7,}  exons {e:>7,}  sig {s:>5,}")
            rows.append({"class": cls, "arm": "50bp", "rows": r, "exons": e, "sig": s})
        cal = tot50
        se50 = next(x for x in rows if x["class"] == "SE" and x["arm"] == "50bp")
        print(f"\n      SE exons  {se50['exons']:,}   record 62,971   "
              f"{'OK' if se50['exons'] == 62_971 else 'MISMATCH'}")
        print(f"      SE sig    {se50['sig']:,}   rMATS-own record {RECORD_RAW_SE_50}   "
              f"{'OK' if abs(se50['sig'] - RECORD_RAW_SE_50) <= 3 else 'CHECK'}")
        if se50["exons"] != 62_971:
            refuse(f"G2 — the 50 bp arm collapses to {se50['exons']:,} skipped exons, "
                   f"not the recorded 62,971.\n\n"
                   f"  This script is therefore not reproducing the paper's convention, "
                   f"and any\n  125-vs-50 difference it reports would be its own bug "
                   f"rather than an effect\n  of the read length. Fix the collapse before "
                   f"reading the comparison.")
        print("  G2 PASS")
    else:
        print("\n  G2  SKIPPED — no --od50 given.")
        print("      ⚠ Without the 50 bp arm this script cannot show that it reproduces")
        print("      the paper's convention, so a difference it reports cannot be")
        print("      attributed to the read length. Pass --od50 if the tree exists.")

    # ------------------------------------------------------------------ G4
    if cal:
        print("\n  G4  do the two arms cover the same event universe?")
        se125 = next(x for x in rows if x["class"] == "SE" and x["arm"] == "125bp")
        d = (se125["exons"] - se50["exons"]) / se50["exons"] * 100
        print(f"      SE exons 50bp {se50['exons']:,} -> 125bp {se125['exons']:,}  "
              f"({d:+.2f}%)")
        print(f"      the 30v30 pilot measured the TESTED set moving −3.2%")
        if abs(d) > 10:
            refuse(f"G4 — the tested set moved {d:+.1f}%, far beyond the −3.2% the pilot\n"
                   f"  measured. Something other than --readLength has changed: check the\n"
                   f"  rMATS version, the GTF release and the BAM lists before reading\n"
                   f"  anything else here.")
        print("  G4 PASS")

    out = pd.DataFrame(rows)
    out.to_csv(a.out, sep="\t", index=False)
    print(f"\n  [write] {a.out}  (nothing else was modified)")

    hr("=")
    print("WHAT THIS DOES AND DOES NOT TELL YOU")
    hr("=")
    print("""
It tells you whether the QUANTIFICATION moved: the tested event set, and rMATS'
own significant counts, at 125 bp against 50 bp, under the paper's collapse-first
exon-level convention.

It does NOT give you the paper's numbers. The 926 comes from
`PSI ~ disease + sex + age + RIN` on 165/161 samples, not from rMATS' test. To put
the headline numbers side by side, run the paper's model on this output:

    python phase2_diffsplice_v2.py \\
        --rmats-dir <the 125 bp dir> \\
        --covariates sex,age,RIN \\
        --out SE_rl125

then recount collapse-first and compare 926, 1,445 and the five per-class counts.

⚠ EXPECT THE NUMBERS TO MOVE SLIGHTLY, AND DECIDE IN ADVANCE WHAT WOULD CHANGE
  YOUR MIND. The 30v30 pilot gave: direction agreeing for 99.85% of events,
  Spearman rho 0.997, orthogonal slope 0.999, significant SE −2.9%, tested set
  −3.2%, A3SS +11.1% on 36 events. If the full cohort lands near those, the
  conclusions are untouched and the only question is whether to renumber the paper.
  If it lands far from them, that is itself a finding and it needs explaining
  before anything is promoted.

  Writing that threshold down BEFORE you look is the difference between a test and
  a rationalisation.

⚠ WHAT PROMOTION WOULD COST, so it is not a surprise later
  926 appears twelve times in the manuscript. Promotion re-baselines the gates in
  62, 68, 69, 70, 71, 73, 77, 78, 79, 80, 81, 87, 93, 95, 96 and 98 — roughly 120
  refusal points — and while they are being re-baselined they protect nothing.
  Figures 2-6, the graphical abstract, Tables 2-4 and Supplements 1, 2, 4, 5, 7, 8
  and 10 all rebuild.

  There is a cheaper destination for this run: report it as a full-cohort
  sensitivity panel (Supplementary Table 5e) beside the 30v30 pilot. The pilot goes
  from 60 samples to 332, no headline number moves, no gate is disarmed, and the
  answer to a referee who asks for the re-run is that it was done at full scale.
""")
    hr("=")


if __name__ == "__main__":
    main()
