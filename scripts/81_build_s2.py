#!/usr/bin/env python
"""
81_build_s2.py -- Supplementary Table S2: all five event classes, every contrast.

⚠ WHY THIS IS A REBUILD AND NOT A TRANSCRIPTION
-----------------------------------------------
`supplementary_plan.md` points S2 at `RESULT_all_five_event_types.md`. **That doc is from
8 September and every count in it is FULL-COORDINATE**, while the manuscript counts
exon-level everywhere:

    DCM vs HCM, "6 significant of 147,355 tested"   <- full-coordinate, RETIRED
    exon level                                       0 significant of  45,380

    composition base counts 1074 / 290 / 90 / 86 / 117   <- full-coordinate
    exon level               926 / 231 / 90 / 84 / 114

Transcribing that doc would put a sixth convention mismatch into the manuscript, after
four have already reached a headline number. So every number in S2 is recomputed from the
`*.all_events.tsv` tables under the locked rule and gated against the figures the
manuscript actually quotes. Nothing is typed in from prose.

(A3SS happens to be 90 in BOTH conventions, which is exactly the kind of coincidence that
makes a transcription look verified when it is not.)

WHAT S2 CONTAINS
----------------
  (a) per-class TESTED and SIGNIFICANT counts for every contrast in the paper, exon
      level, with the two MXE exclusions marked and the reason given
  (b) the equivalence ratios -- median |dPSI(DCM-HCM)| / |dPSI(DCM-NF)| over
      DCM-significant events, per class. The record's 0.19-0.31 is full-coordinate;
      this recomputes them exon-level and prints both so the supplement and the
      Results cannot disagree.

**S2 deliberately does NOT repeat the per-class composition table.** That is S1(a), built
by `76_build_s1_v2.py` across five classes and three models. Duplicating it in two
supplements is how two numbers for one quantity get into print.

GATES (refuses rather than guesses)
-----------------------------------
  G1  discovery per class = 926 / 231 / 90 / 84 / 114
  G2  discovery tested per class = 62,971 / 62,762 / 15,292 / 8,988 / 6,102, summing
      to exactly 156,115
  G3  MAGNet DCM vs HCM SE = 0 of 45,380;  HCM vs NF SE = 468 of 48,024
  G4  Colorado SE = 476 of 19,953 (DCM-NF), 407 of 17,019 (ICM-NF), 2 of 19,998 (DCM-ICM)
  G5  Colorado four usable classes total = 846 / 863 / 2

USAGE
-----
  conda activate rnaseq_tools
  cd /data/project/chureeratp/RNA/HEART

  python 81_build_s2.py --list        # which tables exist, then stop
  python 81_build_s2.py               # build, gate, deposit S2.tsv
  python 81_build_s2.py --allow-missing   # build a partial table anyway
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd

NKEEP = {"SE": 3, "MXE": 5, "A3SS": 5, "A5SS": 5, "RI": 3}
CLASSES = ["SE", "MXE", "A3SS", "A5SS", "RI"]
FDR_MAX, DPSI_MIN = 0.05, 0.10
C_KEY = "key"
C_DPSI_CANDIDATES = ("dPSI_adj", "dPSI_raw", "dPSI")
SUF = os.environ.get("TBL_SUFFIX", "")   # [178] read-length aware
OUT = f"S2{SUF}.tsv"

# contrast -> filename template. LBNL appears nowhere: 24 of its 33 samples are MAGNet
# individuals, so it is excluded from every table and figure in the paper.
CONTRASTS = [
    ("MAGNet",   "DCM vs NF",  "{ET}_v2filtered" + SUF + ".all_events.tsv",   True),
    ("MAGNet",   "HCM vs NF",  ("{ET}_HCMvsNF" + SUF if SUF else "{ET}_HCM_v2")
                               + ".all_events.tsv",                           False),
    ("MAGNet",   "DCM vs HCM", "{ET}_DCMvsHCM" + SUF + ".all_events.tsv",     False),
    ("Colorado", "DCM vs NF",  "{ET}_Colorado_DCMvsNF_dedup.all_events.tsv",   False),
    ("Colorado", "ICM vs NF",  "{ET}_Colorado_ICMvsNF_dedup.all_events.tsv",   False),
    ("Colorado", "DCM vs ICM", "{ET}_Colorado_DCMvsICM_dedup.all_events.tsv",  False),
    ("Brigham",  "HCM vs donor",  "Brigham_{ET}_HCM.all_events.tsv",           False),
    ("Fuwai",    "HCM vs donor",  "Fuwai_{ET}_HCM.all_events.tsv",             False),
    ("Stanford", "HCM vs donor",  "Stanford_{ET}_HCM_dedup.all_events.tsv",    False),
]

# MXE is excluded for these cohorts on reverse-replication grounds, project-wide.
MXE_EXCLUDED = {
    "Colorado": "reverse-replication slope 0.061 (n = 1,477) against SE's 0.447",
    "Fuwai": "lowest MXE correlation in the study, r = 0.088 (paired-end 150 bp)",
}

# ------------------------------------------------------------------ the record
G_SIG = {"SE": 926, "MXE": 231, "A3SS": 90, "A5SS": 84, "RI": 114}
G_TESTED = {"SE": 62_971, "MXE": 62_762, "A3SS": 15_292, "A5SS": 8_988, "RI": 6_102}
G_TESTED_TOTAL = 156_115
G_OTHER = {
    ("MAGNet", "DCM vs HCM", "SE"): (0, 45_380),
    ("MAGNet", "HCM vs NF", "SE"): (468, 48_024),
    ("Colorado", "DCM vs NF", "SE"): (476, 19_953),
    ("Colorado", "ICM vs NF", "SE"): (407, 17_019),
    ("Colorado", "DCM vs ICM", "SE"): (2, 19_998),
}
G_COLORADO_4 = {"DCM vs NF": 846, "ICM vs NF": 863, "DCM vs ICM": 2}
# the record's equivalence ratios -- FULL-COORDINATE, superseded by this build
REC_EQUIV_FULLCOORD = {"SE": 0.220, "MXE": 0.190, "A3SS": 0.238,
                       "A5SS": 0.204, "RI": 0.308}

# [178] read-length aware. The block above is the 50 bp record and stays the default.
# Under a suffix the MAGNet rows are re-based from the 125 bp record; the six
# non-requantified cohorts keep their values because their files did not change.
if SUF == "_rl125":
    G_SIG = {"SE": 813, "MXE": 229, "A3SS": 87, "A5SS": 74, "RI": 91}
    G_TESTED = {"SE": 53_456, "MXE": 62_537, "A3SS": 14_346, "A5SS": 8_360,
                "RI": 6_243}
    G_TESTED_TOTAL = 144_942
    G_OTHER[("MAGNet", "DCM vs HCM", "SE")] = (0, 38_273)
    G_OTHER[("MAGNet", "HCM vs NF", "SE")] = (431, 39_607)
    _fuwai_r = os.environ.get("S2_FUWAI_R", "").strip()
    if not _fuwai_r:
        sys.exit("[REFUSE] TBL_SUFFIX=_rl125 but S2_FUWAI_R is unset. The Fuwai MXE\n"
                 "exclusion text quotes a correlation that must be read, not assumed:\n"
                 "  awk -F'\\t' '$1==\"Fuwai\"{print $16, $21}' symmetry_MXE_rl125.tsv\n"
                 "(unsel_r and prec_r; the 50 bp values were 0.088 / 0.111), then\n"
                 "  TBL_SUFFIX=_rl125 S2_FUWAI_R=<unsel_r> python 81_build_s2.py")
    MXE_EXCLUDED = {
        "Colorado": "reverse-replication slope 0.060 (n = 1,478) against SE's 0.443",
        "Fuwai": f"lowest MXE correlation in the study, r = {_fuwai_r} (paired-end 150 bp)",
    }
elif SUF:
    sys.exit(f"[REFUSE] TBL_SUFFIX={SUF!r}: no record for this suffix. Add one to the\n"
             f"block marked [178] before building under it.")


def die(msg):
    print("\n[REFUSE] " + msg + "\n", file=sys.stderr)
    sys.exit(2)


def ok(msg):
    print("  [ok] " + msg)


def resolve_dpsi(cols):
    for c in C_DPSI_CANDIDATES:
        if c in cols:
            return c
    cand = [c for c in cols if str(c).lower().startswith("dpsi")]
    if not cand:
        return None
    nonraw = [c for c in cand if "raw" not in str(c).lower()]
    return (nonraw or cand)[0]


def load(path, et):
    """Collapse to exon level under the locked convention. None if absent."""
    if not (os.path.exists(path) and os.path.getsize(path)):
        return None
    t = pd.read_csv(path, sep="\t", low_memory=False)
    dcol = resolve_dpsi(t.columns)
    if dcol is None:
        die(f"{path} has no dPSI-like column. Columns: {list(t.columns)}")
    for c in (C_KEY, "FDR"):
        if c not in t.columns:
            die(f"{path} has no '{c}' column")
    t = t.copy()
    t["dPSI"] = pd.to_numeric(t[dcol], errors="coerce")
    t["FDR"] = pd.to_numeric(t["FDR"], errors="coerce")
    nk = NKEEP[et]
    t["ekey"] = [":".join(str(k).split(":")[:nk]) for k in t[C_KEY]]
    t = t.dropna(subset=["FDR", "dPSI"])
    return (t.sort_values("FDR", kind="mergesort")
             .drop_duplicates(subset=["ekey"], keep="first"))


def counts(rep):
    tested = len(rep)
    sig = int(((rep.FDR < FDR_MAX) & (rep.dPSI.abs() >= DPSI_MIN)).sum())
    sig_fdr = int((rep.FDR < FDR_MAX).sum())
    return tested, sig, sig_fdr


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=OUT)
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--allow-missing", action="store_true")
    a = ap.parse_args()

    if a.list:
        print("=" * 78)
        print("S2 inputs")
        print("=" * 78)
        missing = 0
        for cohort, contrast, tmpl, _ in CONTRASTS:
            row = []
            for et in CLASSES:
                p = tmpl.format(ET=et)
                if os.path.exists(p) and os.path.getsize(p):
                    row.append(f"{et}:ok")
                else:
                    row.append(f"{et}:MISSING")
                    missing += 1
            print(f"  {cohort:<9} {contrast:<12}  " + "  ".join(row))
        print(f"\n  {missing} missing. Template examples:")
        for _, _, tmpl, _ in CONTRASTS:
            print(f"    {tmpl.format(ET='SE')}")
        print("\n  MXE will be marked EXCLUDED for: " + ", ".join(MXE_EXCLUDED))
        print("  LBNL appears nowhere -- 24 of 33 samples are MAGNet individuals.")
        return

    print("=" * 78)
    print("SUPPLEMENTARY TABLE S2 -- all five event classes, every contrast")
    print("=" * 78)
    print("  exon level, NKEEP per class, lowest-FDR row kept")
    print(f"  significant = FDR < {FDR_MAX} and |dPSI| >= {DPSI_MIN}")

    rows, fail, warn = [], [], []
    tables = {}

    for cohort, contrast, tmpl, is_disc in CONTRASTS:
        for et in CLASSES:
            path = tmpl.format(ET=et)
            rep = load(path, et)
            if rep is None:
                warn.append(f"missing: {path}")
                continue
            tables[(cohort, contrast, et)] = rep
            tested, sig, sig_fdr = counts(rep)
            excl = (et == "MXE" and cohort in MXE_EXCLUDED)
            rows.append(dict(
                cohort=cohort, contrast=contrast, event_class=et,
                tested=tested, significant=sig, significant_fdr_only=sig_fdr,
                pct_significant=round(100 * sig / tested, 4) if tested else np.nan,
                excluded_from_paper="yes" if excl else "no",
                exclusion_reason=MXE_EXCLUDED[cohort] if excl else "",
                source_file=path, is_discovery="yes" if is_disc else "no",
            ))

            # ---- gates
            if is_disc:
                if sig != G_SIG[et]:
                    fail.append(f"G1 {et}: {sig} significant, expected {G_SIG[et]}")
                if tested != G_TESTED[et]:
                    fail.append(f"G2 {et}: {tested:,} tested, "
                                f"expected {G_TESTED[et]:,}")
            key = (cohort, contrast, et)
            if key in G_OTHER:
                es, et_ = G_OTHER[key]
                if (sig, tested) != (es, et_):
                    fail.append(f"G3/G4 {cohort} {contrast} {et}: "
                                f"{sig} of {tested:,}, expected {es} of {et_:,}")

    if not rows:
        die("no input tables found at all; run --list")

    df = pd.DataFrame(rows)

    # G2 total
    disc = df[df.is_discovery == "yes"]
    if len(disc) == 5:
        tot = int(disc.tested.sum())
        if tot != G_TESTED_TOTAL:
            fail.append(f"G2 total: {tot:,} tested, expected {G_TESTED_TOTAL:,}")
        else:
            ok(f"G2 discovery tested sums to exactly {tot:,}")
        if int(disc.significant.sum()) != sum(G_SIG.values()):
            fail.append(f"G1 total: {int(disc.significant.sum())} significant, "
                        f"expected {sum(G_SIG.values())}")
        else:
            ok(f"G1 discovery significant = {int(disc.significant.sum())}")
    else:
        warn.append(f"only {len(disc)} of 5 discovery classes present")

    # G5 Colorado four usable classes (MXE excluded)
    for contrast, expect in G_COLORADO_4.items():
        sub = df[(df.cohort == "Colorado") & (df.contrast == contrast)
                 & (df.event_class != "MXE")]
        if len(sub) == 4:
            got = int(sub.significant.sum())
            if got != expect:
                fail.append(f"G5 Colorado {contrast} four-class total: "
                            f"{got}, expected {expect}")
            else:
                ok(f"G5 Colorado {contrast} four usable classes = {got}")

    if fail:
        print("\n  GATE FAILURES:")
        for f in fail:
            print("    ✗ " + f)
        die(f"{len(fail)} gate(s) failed. S2 not written. A gate failure here means a "
            "table is not the one the manuscript's numbers came from -- find out which "
            "before changing the expected value.")
    ok("all gates passed")

    # ---------------------------------------------------------------- (a)
    print("\n" + "-" * 78)
    print("  (a) TESTED AND SIGNIFICANT, per class and contrast")
    print("-" * 78)
    print(f"  {'cohort':<9}{'contrast':<13}" +
          "".join(f"{et:>13}" for et in CLASSES))
    for (cohort, contrast), g in df.groupby(["cohort", "contrast"], sort=False):
        cells = []
        for et in CLASSES:
            r = g[g.event_class == et]
            if r.empty:
                cells.append(f"{'-':>13}")
            else:
                r = r.iloc[0]
                s = f"{int(r.significant)}/{int(r.tested):,}"
                if r.excluded_from_paper == "yes":
                    s = "[" + s + "]"
                cells.append(f"{s:>13}")
        print(f"  {cohort:<9}{contrast:<13}" + "".join(cells))
    print("\n  [brackets] = excluded from every table and figure in the paper")
    for c, why in MXE_EXCLUDED.items():
        print(f"    {c} MXE: {why}")

    # ---------------------------------------------------------------- (b)
    print("\n" + "-" * 78)
    print("  (b) EQUIVALENCE RATIOS -- recomputed at exon level")
    print("-" * 78)
    print("  median |dPSI(DCM-HCM)| / |dPSI(DCM-NF)| over DCM-significant events.")
    print("  The record's 0.19-0.31 is FULL-COORDINATE; both are shown so the")
    print("  supplement and the Results cannot disagree.\n")
    print(f"  {'class':<7}{'exon level':>12}{'full-coord (record)':>22}"
          f"{'n events':>10}")
    eq_rows = []
    for et in CLASSES:
        d = tables.get(("MAGNet", "DCM vs NF", et))
        h = tables.get(("MAGNet", "DCM vs HCM", et))
        if d is None or h is None:
            print(f"  {et:<7}{'-':>12}{REC_EQUIV_FULLCOORD[et]:>22.3f}{'-':>10}")
            continue
        dsig = d[(d.FDR < FDR_MAX) & (d.dPSI.abs() >= DPSI_MIN)]
        j = dsig.merge(h[["ekey", "dPSI"]].rename(columns={"dPSI": "dPSI_hcm"}),
                       on="ekey", how="inner")
        j = j[j.dPSI.abs() > 0]
        if len(j) < 5:
            print(f"  {et:<7}{'too few':>12}{REC_EQUIV_FULLCOORD[et]:>22.3f}"
                  f"{len(j):>10}")
            continue
        ratio = float((j.dPSI_hcm.abs() / j.dPSI.abs()).median())
        print(f"  {et:<7}{ratio:>12.3f}{REC_EQUIV_FULLCOORD[et]:>22.3f}{len(j):>10,}")
        eq_rows.append(dict(cohort="MAGNet", contrast="equivalence ratio",
                            event_class=et, tested=pd.NA,
                            equiv_n_events=len(j), significant=pd.NA,
                            significant_fdr_only=pd.NA,
                            pct_significant=np.nan,
                            excluded_from_paper="no", exclusion_reason="",
                            source_file="DCM vs NF x DCM vs HCM",
                            is_discovery="no",
                            equiv_ratio_exon=round(ratio, 4),
                            equiv_ratio_fullcoord_record=REC_EQUIV_FULLCOORD[et]))
    if eq_rows:
        vals = [r["equiv_ratio_exon"] for r in eq_rows]
        print(f"\n  >>> exon-level range {min(vals):.2f}-{max(vals):.2f} "
              f"(record, full-coordinate: 0.19-0.31)")
        print("      If this range differs materially, the Results sentence quoting")
        print("      0.19-0.31 must be updated -- it is currently a full-coordinate")
        print("      figure sitting in exon-level prose, which is the error that has")
        print("      reached a headline number four times in this project.")

    out = pd.concat([df, pd.DataFrame(eq_rows)], ignore_index=True) if eq_rows else df
    # Deposit counts as integers. "926.0" in a supplementary table is the kind of
    # detail a reviewer reads as carelessness, and the float came only from the NaN
    # rows the equivalence block adds.
    for c in ("tested", "significant", "significant_fdr_only", "equiv_n_events"):
        if c in out.columns:
            out[c] = out[c].astype("Int64")
    out.to_csv(a.out, sep="\t", index=False)

    print("\n" + "=" * 78)
    print(f"S2 BUILT -- {a.out}, {len(out)} rows")
    print("=" * 78)
    if warn:
        print(f"  {len(warn)} warning(s):")
        for w in warn[:12]:
            print("    - " + w)
        if len(warn) > 12:
            print(f"    ... and {len(warn) - 12} more")
    print("  CAPTION MUST STATE, because the table cannot say it for itself:")
    print("   * The two MXE exclusions rest on DIFFERENT evidence. Colorado MXE is")
    print("     inflated AND fails reverse replication; Fuwai MXE is NOT inflated")
    print("     (0.48%, in line with MAGNet's 0.37%) and is excluded on its")
    print("     reverse-replication correlation alone, r = 0.088. Do not let the")
    print("     caption imply Fuwai shows inflation it does not show.")
    print("   * The Colorado DCM-vs-ICM MXE zero appears in this table and is NOT")
    print("     evidence. It is excluded project-wide; a reader must not count it as")
    print("     a sixth null contrast.")
    print("   * The equivalence ratio is computed on discovery-significant events")
    print("     MEASURABLE IN BOTH contrasts, which is a subset -- 774 of the 926")
    print("     skipped exons, and less in the smaller classes. State the subset.")
    print("   * The five-class DCM-vs-HCM figure is 6 significant of 113,400 tested")
    print("     (11 at FDR alone). The numerator is unchanged from the retired")
    print("     full-coordinate '6 of 147,355'; only the denominator moved.")
    print()
    print("\n  S2 does NOT duplicate the per-class composition table -- that is S1(a),")
    print("  from 76_build_s1_v2.py. Two supplements carrying one quantity is how two")
    print("  numbers for it get into print.")
    print("  LBNL appears nowhere, by design.")
    print()


if __name__ == "__main__":
    main()
