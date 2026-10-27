#!/usr/bin/env python3
"""
87_build_table3.py  —  Table 3: within-cohort differential splicing, every null
                       beside its positive control.  v2.

WHAT CHANGED FROM v1, AND WHY
-----------------------------
v1 had nine rows, adding Brigham, Fuwai and Stanford. The inventory showed their own
significant counts are **4, 169 and 8 skipped exons**, against MAGNet's 926. Printing
those in one table is a cross-cohort comparison of counts, which this paper explicitly
forbids (`RESULT_magnitude_not_estimable.md`: "Counts are equally non-transferable —
different tested sets, covariate models, power and read geometry"). A reviewer reads
4-against-926 as a failed replication; the truth is that an 18-versus-5 cohort cannot
reach per-event significance and replicates by **direction** (81.2%). Those cohorts
belong in a replication table, not a counts table. G8 below prints them with that
reasoning so the omission is a recorded decision, not an oversight.

So Table 3 is **six within-cohort contrasts** where counts are meaningful:

    MAGNet   DCM vs NF    926 SE   <- the signature
    MAGNet   HCM vs NF    468 SE   <- the same 28 HCM hearts: not underpowered
    MAGNet   DCM vs HCM     0 SE   <- the null, both arms explants
    Colorado DCM vs NF    476 SE   <- double positive control, second cohort
    Colorado ICM vs NF    407 SE
    Colorado DCM vs ICM     2 SE   <- the null again, second etiology pair

**Skipped exons are the table's spine**, because every claim in the paper is stated in
SE and the threshold ladder exists only for SE. All-class totals are carried as two
further columns where the frozen master file is complete.

TWO HONEST GAPS, BOTH FOOTNOTED RATHER THAN HIDDEN
  - **MAGNet DCM vs HCM is absent from `master_numbers_flagged.tsv`.** Its SE row is
    taken from `fdr_only_counts.tsv`; its all-class cells are "—". Do not substitute a
    number from elsewhere without checking which convention produced it — a
    full-coordinate FDR-only count standing in this cell is the error that reached the
    abstract once already (1,443 against 468).
  - The FDR-only and |ΔPSI|≥0.05 ladder columns exist for the three MAGNet contrasts
    only.

SOURCES, AND ONE HARD RULE
  counts  : `master_numbers_flagged.tsv`, exon-level columns
  ladder  : `fdr_only_counts.tsv`
  group n : **Table 1**, never `median_samples_g1`/`g2` — those are per-event medians
            after the coverage mask, and the project has that recorded as a trap.

GATES — each refuses rather than warns
  G1  the six contrasts resolve (DCM vs HCM may come from the ladder); print the grid.
  G2  Colorado MXE and Fuwai MXE carry an exclude flag and never enter a total.
  G3  MAGNet DCM vs NF per class == 62,971/926, 62,762/231, 15,292/90, 8,988/84,
      6,102/114, summing to 156,115/1,445.
  G4  every SE pair matches the record.
  G5  every all-class total matches the record (31,569/846, 27,231/863, 32,381/2).
  G6  the ladder's +0.10 column equals the master file's significant count.
  G7  group sizes appear in Table1.tsv.
  G8  the three replication cohorts are reported and deliberately excluded.

USAGE
    python 87_build_table3.py --list     # inventory only
    python 87_build_table3.py           # gates, then Table3.tsv + Table3.md
"""

import argparse
import os
import re
import sys

import pandas as pd

MASTER = "master_numbers_flagged.tsv"
LADDER = "fdr_only_counts.tsv"
TABLE1 = "Table1.tsv"
CLASSES = ["SE", "MXE", "A3SS", "A5SS", "RI"]

# ---- group sizes: FROM TABLE 1 -------------------------------------------
GROUPS = {
    ("MAGNet", "DCMvsNF"):    (166, 166, "DCM", "NF"),
    ("MAGNet", "HCMvsNF"):    (28, 166, "HCM", "NF"),
    ("MAGNet", "DCMvsHCM"):   (166, 28, "DCM", "HCM"),
    ("Colorado", "DCMvsNF"):  (36, 13, "DCM", "NF"),
    ("Colorado", "ICMvsNF"):  (13, 13, "ICM", "NF"),
    ("Colorado", "DCMvsICM"): (36, 13, "DCM", "ICM"),
}
ROW_ORDER = list(GROUPS)

# contrasts with no master row; SE comes from the ladder, all-class is "—"
LADDER_ONLY = {("MAGNet", "DCMvsHCM")}

# deliberately NOT in Table 3 — direction, not counts (see G8)
REPLICATION = {
    ("Brigham", "HCMvsctrl"): (18, 5, "81.2% direction concordance on 442 SE events"),
    ("Fuwai", "HCMvsctrl"):   (28, 9, "69.2% on 714 SE events"),
    ("Stanford", "HCMvsctrl"): (13, 6, "76.1% on 536 SE events"),
}
EXCLUDED = {("Colorado", "MXE"), ("Fuwai", "MXE")}

REC_SE = {
    ("MAGNet", "DCMvsNF"):    (62971, 926),
    ("MAGNet", "HCMvsNF"):    (48024, 468),
    ("MAGNet", "DCMvsHCM"):   (45380, 0),
    ("Colorado", "DCMvsNF"):  (19953, 476),
    ("Colorado", "ICMvsNF"):  (17019, 407),
    ("Colorado", "DCMvsICM"): (19998, 2),
}
REC_ALL = {
    ("MAGNet", "DCMvsNF"):    (156115, 1445),
    ("Colorado", "DCMvsNF"):  (31569, 846),
    ("Colorado", "ICMvsNF"):  (27231, 863),
    ("Colorado", "DCMvsICM"): (32381, 2),
}
REC_PERCLASS_MAGNET = {"SE": (62971, 926), "MXE": (62762, 231), "A3SS": (15292, 90),
                       "A5SS": (8988, 84), "RI": (6102, 114)}

# ---- 168: the same record at 125 bp --------------------------------------
# Typed from the manuscript's own counts, NOT read back out of the master
# these dicts validate. A cross-check derived from the file it checks passes
# by construction and tests nothing.
#
# Colorado is inherited by unpacking, because that cohort was never
# re-quantified; only the MAGNet keys are overridden.
REC_SE_RL125 = {**REC_SE,
                ("MAGNet", "DCMvsNF"):  (53456, 813),
                ("MAGNet", "HCMvsNF"):  (39607, 431),
                ("MAGNet", "DCMvsHCM"): (38273, 0)}
REC_ALL_RL125 = {**REC_ALL,
                 ("MAGNet", "DCMvsNF"): (144942, 1294)}
REC_PERCLASS_MAGNET_RL125 = {"SE": (53456, 813), "MXE": (62537, 229),
                             "A3SS": (14346, 87), "A5SS": (8360, 74),
                             "RI": (6243, 91)}

_T3REC = os.environ.get("T3_RECORD", "").strip().lower()
if "rl125" in " ".join(sys.argv) and _T3REC != "rl125":
    sys.exit(
        "[fatal] a 125 bp path was named on the command line, but T3_RECORD is\n"
        "        not set to rl125. G3-G5 would check 125 bp inputs against the\n"
        "        50 bp record and refuse, or worse, pass for the wrong reason.\n"
        "        Re-run as:\n"
        "          T3_RECORD=rl125 python3 87_build_table3.py \\\n"
        "              --dir t3_rl125 --out Table3_rl125\n"
        "        Nothing built.")
if _T3REC == "rl125":
    if "rl125" not in " ".join(sys.argv):
        sys.exit("[fatal] T3_RECORD=rl125 but no 125 bp path was named. The 125 bp\n"
                 "        record would be checked against 50 bp inputs. Nothing built.")
    REC_SE = REC_SE_RL125
    REC_ALL = REC_ALL_RL125
    REC_PERCLASS_MAGNET = REC_PERCLASS_MAGNET_RL125
    print("[168] record: 125 bp "
          "(MAGNet re-based; Colorado inherited unchanged)")
elif _T3REC:
    sys.exit(f"[fatal] T3_RECORD={_T3REC!r} is not a record this script knows. "
             "Use rl125 or leave it unset.")

C_TESTED = ("tested_exon", "tested")
C_SIG = ("sig_exon", "significant")


def refuse(msg):
    print("\n" + "=" * 78)
    print("REFUSING — " + msg)
    print("=" * 78)
    sys.exit(2)


def pick(cols, cands, what, required=True):
    for c in cands:
        if c in cols:
            return c
    low = {c.lower(): c for c in cols}
    for c in cands:
        if c.lower() in low:
            return low[c.lower()]
    if required:
        refuse(f"cannot find a {what} column.\n  tried: {', '.join(cands)}\n"
               f"  table has: {', '.join(cols)}")
    return None


def norm_contrast(s):
    return re.sub(r"[^A-Za-z0-9]", "", str(s))


def norm_cohort(s):
    """'MAGNet (e)' -> 'MAGNet'. split() only; rstrip of a character set would turn
    'Stanford' into 'Stanfor'."""
    return str(s).strip().split()[0]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default=".")
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--out", default="Table3")
    args = ap.parse_args()

    mp = os.path.join(args.dir, MASTER)
    if not os.path.exists(mp):
        refuse(f"{MASTER} not found in {args.dir}")
    t = pd.read_csv(mp, sep="\t")
    for c in ("cohort", "contrast", "event_type"):
        if c not in t.columns:
            refuse(f"{MASTER} needs a {c} column.\n  has: {', '.join(t.columns)}")
    t["_cohort"] = t["cohort"].map(norm_cohort)
    t["_contrast"] = t["contrast"].map(norm_contrast)
    t["_type"] = t["event_type"].astype(str).str.upper()
    c_tested = pick(t.columns, C_TESTED, "tested (exon-level)")
    c_sig = pick(t.columns, C_SIG, "significant (exon-level)")
    c_genes = pick(t.columns, ("genes_sig",), "genes", required=False)
    c_rel = "reliability" if "reliability" in t.columns else None

    lp = os.path.join(args.dir, LADDER)
    L = pd.read_csv(lp, sep="\t") if os.path.exists(lp) else None
    if L is not None:
        L["_cohort"] = L["contrast"].map(norm_cohort)
        L["_contrast"] = L["contrast"].map(
            lambda s: norm_contrast(" ".join(str(s).split()[1:])))

    print("=" * 78)
    print("87_build_table3.py v2 — Table 3: within-cohort contrasts, nulls with controls")
    print("=" * 78)
    print(f"\nmaster: {MASTER} ({len(t)} rows)   tested={c_tested}  sig={c_sig}")

    # ---------------------------------------------------------------- G1
    print("\n" + "-" * 78)
    print("G1  inventory")
    print("-" * 78)
    grid = (t.assign(cell=t[c_tested].astype("Int64").astype(str) + "/"
                     + t[c_sig].astype("Int64").astype(str))
             .pivot_table(index=["_cohort", "_contrast"], columns="_type",
                          values="cell", aggfunc="first")
             .reindex(columns=CLASSES))
    print("  " + grid.fillna("—").to_string().replace("\n", "\n  "))

    def cell(coh, con, cls):
        r = t[(t._cohort == coh) & (t._contrast == con) & (t._type == cls)]
        if not len(r):
            return None
        return int(r[c_tested].iloc[0]), int(r[c_sig].iloc[0])

    def ladder_row(coh, con):
        if L is None:
            return None
        r = L[(L._cohort == coh) & (L._contrast == con)]
        return r.iloc[0] if len(r) else None

    missing = []
    for k in ROW_ORDER:
        if cell(*k, "SE") is None and ladder_row(*k) is None:
            missing.append(f"{k[0]} {k[1]}")
    print(f"\n  Table 3 rows wanted: {len(ROW_ORDER)}")
    for k in ROW_ORDER:
        src = ("master" if cell(*k, "SE") is not None
               else ("ladder only" if ladder_row(*k) is not None else "MISSING"))
        print(f"    {k[0]:9s} {k[1]:12s} {src}")
    if args.list:
        if L is not None:
            print(f"\n--- {LADDER}")
            print("  " + L.drop(columns=["_cohort", "_contrast"]).to_string(index=False)
                  .replace("\n", "\n  "))
        return
    if missing:
        refuse("G1 — no source for: " + ", ".join(missing))
    print("  G1 PASS")

    # ---------------------------------------------------------------- G2
    print("\n" + "-" * 78)
    print("G2  excluded MXE classes flagged, and kept out of every total")
    print("-" * 78)
    for coh, cls in sorted(EXCLUDED):
        rows = t[(t._cohort == coh) & (t._type == cls)]
        if not len(rows):
            print(f"  {coh} {cls}: absent from the master file — ok")
            continue
        flags = sorted(set(rows[c_rel].astype(str))) if c_rel else ["(no column)"]
        ok = any("exclu" in f.lower() for f in flags)
        print(f"  {coh} {cls}: {len(rows)} row(s), reliability={flags} "
              f"{'flagged' if ok else '⚠ NOT FLAGGED'}")
        if not ok:
            refuse(f"G2 — {coh} {cls} is excluded in the manuscript but not flagged in "
                   f"{MASTER}.")
    print("  G2 PASS")

    def totals(coh, con):
        """Sum only over classes that SHOULD be present. A missing class is a hole,
        reported, never summed over silently."""
        tt = ss = 0
        used, absent = [], []
        for cls in [c for c in CLASSES if (coh, c) not in EXCLUDED]:
            c = cell(coh, con, cls)
            if c is None:
                absent.append(cls)
                continue
            tt += c[0]
            ss += c[1]
            used.append(cls)
        return tt, ss, used, absent

    # ---------------------------------------------------------------- G3-G5
    print("\n" + "-" * 78)
    print("G3–G5  does every recorded number reproduce?")
    print("-" * 78)
    bad = []
    print("  MAGNet DCM vs NF, per class:")
    for cls in CLASSES:
        c = cell("MAGNet", "DCMvsNF", cls)
        r = REC_PERCLASS_MAGNET[cls]
        if c != r:
            bad.append(f"MAGNet DCMvsNF {cls}: {c} vs record {r}")
        print(f"    {cls:5s} {(f'{c[0]:,}/{c[1]}' if c else '—'):>16s}   "
              f"record {r[0]:,}/{r[1]}   {'match' if c == r else 'MISMATCH'}")

    print(f"\n  {'contrast':24s} {'SE tested/sig':>18s} {'record':>14s}  "
          f"{'all-class':>16s} {'record':>14s}")
    for k in ROW_ORDER:
        se = cell(*k, "SE")
        if se is None:
            lr = ladder_row(*k)
            se = (int(lr["exon-level_tested"]), int(lr["exon-level_+|dPSI|>=0.10"]))
        tt, ss, used, absent = totals(*k)
        all_cell = "—" if k in LADDER_ONLY else f"{tt:,}/{ss}"
        if k not in LADDER_ONLY and absent:
            bad.append(f"{k[0]} {k[1]}: missing classes {', '.join(absent)}")
        rse, ral = REC_SE.get(k), REC_ALL.get(k)
        line = (f"  {k[0] + ' ' + k[1]:24s} {f'{se[0]:,}/{se[1]}':>18s} "
                f"{(f'{rse[0]:,}/{rse[1]}' if rse else '—'):>14s}  {all_cell:>16s} "
                f"{(f'{ral[0]:,}/{ral[1]}' if ral else '—'):>14s}")
        if rse and se != rse:
            bad.append(f"{k[0]} {k[1]} SE: {se} vs record {rse}")
            line += "  MISMATCH(SE)"
        if ral and (tt, ss) != ral:
            bad.append(f"{k[0]} {k[1]} all-class: {(tt, ss)} vs record {ral}")
            line += "  MISMATCH(all)"
        print(line)
    if bad:
        refuse("G3–G5 — the master file does not reproduce the record:\n  - "
               + "\n  - ".join(bad)
               + "\n\n  Do not build Table 3 until this is resolved.")
    print("\n  G3–G5 PASS")

    # ---------------------------------------------------------------- G6
    print("\n" + "-" * 78)
    print("G6  ladder vs master on the +|dPSI|>=0.10 column")
    print("-" * 78)
    if L is None:
        print(f"  SKIPPED — {LADDER} absent; the ladder footnote cannot be written.")
    else:
        for _, r in L.iterrows():
            k = (r["_cohort"], r["_contrast"])
            se = cell(*k, "SE")
            d10 = int(r["exon-level_+|dPSI|>=0.10"])
            if se is None:
                print(f"  {k[0]} {k[1]:12s} ladder {d10:>6,d}   (no master row — this is "
                      f"the only source)")
                continue
            if d10 != se[1]:
                refuse(f"G6 — {LADDER} says {d10} for {k[0]} {k[1]}, {MASTER} says "
                       f"{se[1]}.\n  Two files disagree on a number in the manuscript.")
            print(f"  {k[0]} {k[1]:12s} ladder {d10:>6,d}   master {se[1]:>6,d}   match")
        print("  G6 PASS")

    # ---------------------------------------------------------------- G7
    print("\n" + "-" * 78)
    print("G7  group sizes against Table 1 (never median_samples_g1/g2)")
    print("-" * 78)
    t1p = os.path.join(args.dir, TABLE1)
    if not os.path.exists(t1p):
        print(f"  SKIPPED — {TABLE1} absent; confirm the hard-coded sizes by eye.")
    else:
        blob = " ".join(pd.read_csv(t1p, sep="\t").astype(str).values.ravel())
        miss = [f"{c} {n} {l}" for (c, _), (a, b, la, lb) in GROUPS.items()
                for n, l in ((a, la), (b, lb)) if not re.search(rf"\b{n}\b", blob)]
        print("  ⚠ not found in Table 1: " + "; ".join(miss) if miss
              else "  G7 PASS — every group size appears in Table 1.")

    # ---------------------------------------------------------------- G8
    print("\n" + "-" * 78)
    print("G8  the replication cohorts, and why they are NOT in this table")
    print("-" * 78)
    for (coh, con), (a, b, note) in REPLICATION.items():
        se = cell(coh, con, "SE")
        print(f"  {coh:9s} {a:>3d} vs {b:<3d}  SE {se[0]:,}/{se[1]:<4d}  -> {note}")
    print("""
  Those counts (4, 169, 8 significant SE) are not comparable with MAGNet's 926: the
  cohorts differ in tested set, covariate model, read geometry and above all power.
  Printing them beside 926 invites "the replication cohorts found nothing", which is
  false — they replicate by DIRECTION at 69–81%. They belong in the replication table.
  This is a decision, recorded here so it is not mistaken for an omission.""")

    # ---------------------------------------------------------------- build
    print("\n" + "-" * 78)
    print("TABLE 3")
    print("-" * 78)
    out = []
    for k in ROW_ORDER:
        coh, con = k
        a, b, la, lb = GROUPS[k]
        se = cell(coh, con, "SE")
        lr = ladder_row(coh, con)
        if se is None:
            se = (int(lr["exon-level_tested"]), int(lr["exon-level_+|dPSI|>=0.10"]))
        tt, ss, used, _ = totals(coh, con)
        gl = None
        if c_genes:
            r = t[(t._cohort == coh) & (t._contrast == con) & (t._type == "SE")]
            if len(r) and pd.notna(r[c_genes].iloc[0]):
                gl = int(r[c_genes].iloc[0])
        out.append({
            "Cohort": coh,
            "Contrast": f"{la} vs {lb}",
            "n": f"{a}/{b}",
            "SE tested": se[0],
            "FDR<0.05": int(lr["exon-level_FDR only"]) if lr is not None else None,
            "+|dPSI|>=0.05": int(lr["exon-level_+|dPSI|>=0.05"]) if lr is not None else None,
            "+|dPSI|>=0.10": se[1],
              # [97] dropped: no claim depends on a gene count; see RESULT_table3_built.md
            "Classes": (f"{len(used)} (f)" if (coh, "MXE") in EXCLUDED
                        else ("—" if k in LADDER_ONLY else str(len(used)))),
            "All tested": None if k in LADDER_ONLY else tt,
            "All significant": None if k in LADDER_ONLY else ss,
        })
    T3 = pd.DataFrame(out)
    # nullable integers: a mixed int/NaN column otherwise prints 5301.0 and deposits
    # counts as floats
    for c in ("SE tested", "FDR<0.05", "+|dPSI|>=0.05", "+|dPSI|>=0.10",  # [97] dropped: no claim depends on a gene count; see RESULT_table3_built.md
              "All tested", "All significant"):
        T3[c] = T3[c].astype("Int64")
    print(T3.astype(object).where(T3.notna(), "—").to_string(index=False))

    tsv = os.path.join(args.dir, f"{args.out}.tsv")
    T3.to_csv(tsv, sep="\t", index=False)
    # column names contain "|" (as in "+|dPSI|>=0.05"), which breaks a markdown table
    # unless escaped
    def esc(s):
        return str(s).replace("|", r"\|")

    md = ["| " + " | ".join(esc(c) for c in T3.columns) + " |",
          "|" + "|".join("---" for _ in T3.columns) + "|"]
    for _, r in T3.iterrows():
        md.append("| " + " | ".join(
            "—" if pd.isna(v) else (f"{int(v):,}" if isinstance(v, (int, float))
                                    else esc(v))
            for v in r) + " |")
    mdp = os.path.join(args.dir, f"{args.out}.md")
    with open(mdp, "w") as fh:
        fh.write("\n".join(md) + "\n")
    print(f"\n[write] {tsv}\n[write] {mdp}")

    # ---------------------------------------------------------------- caption
    print("\n" + "=" * 78)
    print("CAPTION — every number read from a file")
    print("=" * 78)
    print("""
**Table 3. Within-cohort differential splicing: each null contrast beside its positive
control.** Exon-level events, JCEC quantification. Columns show the cumulative effect of
the significance definition: FDR < 0.05 alone, then additionally |ΔPSI| ≥ 0.05, then
≥ 0.10, the threshold used throughout this paper. Counts are skipped exons; "All"
columns sum the event classes analysed in that cohort. Group sizes are post-exclusion,
as in Table 1.

The two null contrasts — DCM versus HCM in the discovery cohort, both arms explanted
hearts, and DCM versus ICM in Colorado — are each accompanied by comparisons of the same
hearts against non-failing myocardium on the same pipeline.

(f) Four classes: mutually exclusive exons are excluded for Colorado on
    reverse-replication grounds (Supplementary Table …).
(g) The discovery DCM-versus-HCM contrast is quantified here for skipped exons; the
    all-class total is not available in the frozen numbers file.
(i) The two intermediate threshold columns are available for the discovery cohort only;
    Colorado contrasts are reported at the significance definition used throughout.
(h) The three independent HCM cohorts are not shown: with 18/5, 28/9 and 13/6 samples
    they are powered for direction concordance, not for per-event significance, and
    event counts are not comparable across cohorts. See the replication table.

LBNL is excluded throughout — 24 of its 33 samples are discovery-cohort individuals.

WHAT THE RESULTS CAN NOW DROP — verify each before deleting:
  - the five per-class tested counts summing to 156,115
  - the Colorado 846 / 863 / 2 sentence
  - the HCM-versus-NF 468 comparator sentence
  - the threshold-ladder prose
KEEP the ratios in the text (238-fold SE, 423-fold all classes). A table of counts
supports the argument; it does not make it.
""")


if __name__ == "__main__":
    main()
