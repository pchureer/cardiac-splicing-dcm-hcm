#!/usr/bin/env python3
"""
167_build_table3_rl125_inputs.py -- Table 3 is still entirely 50 bp.

Table3.tsv, built 12 Oct:

    MAGNet  DCM vs NF   62971  5301  2592  926  5  156115  1445
    MAGNet  HCM vs NF   48024  1355   909  468  5  118234   695
    MAGNet  DCM vs HCM  45380     3     1    0  —

against a Results section that says 53,456 / 813 / 144,942 / 1,294 / 38,273 /
431. The table and the text contradict each other on the paper's headline
counts. The Colorado rows are correct -- that cohort was never re-quantified.

87_build_table3.py reads master_numbers_flagged.tsv and fdr_only_counts.tsv,
both frozen at 50 bp, which is why none of its eight gates noticed.

NO PATCH TO 87 IS NEEDED. It takes --dir, and resolves MASTER, LADDER and
TABLE1 against it. So this script writes a directory of correct 125 bp inputs
and 87 is run against that:

    python3 167_build_table3_rl125_inputs.py
    python3 87_build_table3.py --dir t3_rl125 --out Table3_rl125

WHAT IT MEASURES, AND WHY IT IS NOT TYPED IN
--------------------------------------------
Four of the needed values are on record (5,401 / 1,430 / 813 / 431). The
|dPSI|>=0.05 ladder column and the HCM all-class totals are recorded NOWHERE,
so every number here is recomputed from the event tables:

1. CALIBRATE. For each 50 bp master row, read THE FILE THAT ROW NAMES in its
   own `file` column and require (tested_exon, sig_exon) to reproduce exactly.
   Then reproduce all eight numbers of each ladder row. If the recipe cannot
   rebuild the frozen inputs it has no business writing new ones.
2. MEASURE the fifteen 125 bp tables.
3. GATE on the seven 125 bp values that ARE on record, plus the three
   FDR-only ladder values. A file that gives the wrong count is the wrong file.
4. WRITE t3_rl125/{master_numbers_flagged,fdr_only_counts,Table1}.tsv

CONVENTION: collapse to the event's own coordinates FIRST keeping the
lowest-FDR row, THEN threshold. NKEEP counts FIELDS of `key`: SE 3, MXE 5,
A3SS 5, A5SS 5, RI 3. full-coordinate = the uncollapsed rows.

⚠ THE TABLE'S SHAPE IS NOT CHANGED. DCM vs HCM stays SE-only with "—" for
Classes, because that is what footnote (g) explains and 87's docstring
requires. All five 125 bp DCM-vs-HCM tables now exist, so filling that row is
possible -- but it is a redesign, not a re-basing, and it would invalidate the
footnote. Raised, not taken.
"""
import argparse
import os
import shutil
import sys

import pandas as pd

TYPES = ["SE", "MXE", "A3SS", "A5SS", "RI"]
NKEEP = {"SE": 3, "MXE": 5, "A3SS": 5, "A5SS": 5, "RI": 3}
MASTER, LADDER, TABLE1 = ("master_numbers_flagged.tsv", "fdr_only_counts.tsv",
                          "Table1.tsv")
MAGNET_ROWS = ["DCM vs NF", "HCM vs NF"]          # the master's MAGNet contrasts
LADDER_ROWS = ["MAGNet DCM vs NF", "MAGNet HCM vs NF", "MAGNet DCM vs HCM"]

# 125 bp filenames, all confirmed present on disk.
F125 = {"DCM vs NF": "{ET}_rl125.all_events.tsv",
        "HCM vs NF": "{ET}_HCMvsNF_rl125.all_events.tsv",
        "DCM vs HCM": "{ET}_DCMvsHCM_rl125.all_events.tsv"}

# Values already on record elsewhere. Used ONLY as gates; what is written is
# what this script measures.
GATE = {("DCM vs NF", "SE"): (53456, 813), ("DCM vs NF", "MXE"): (62537, 229),
        ("DCM vs NF", "A3SS"): (14346, 87), ("DCM vs NF", "A5SS"): (8360, 74),
        ("DCM vs NF", "RI"): (6243, 91), ("HCM vs NF", "SE"): (39607, 431),
        ("DCM vs HCM", "SE"): (38273, 0)}
GATE_FDRONLY = {"DCM vs NF": 5401, "HCM vs NF": 1430, "DCM vs HCM": 1}


def dpsi_col(cols, path):
    """Resolve the effect column. NOT a guess -- two schemas exist in this
    project and they must not be conflated:

      * the standard 20-column phase-2 table carries BOTH dPSI_raw and
        dPSI_adj. dPSI_adj is the covariate-adjusted estimate and is the one
        every count in this paper uses, so it wins outright when present.
      * the DCM-vs-HCM tables are a 12-column product with a single effect
        column named for the contrast (dPSI_DCMvsHCM) and `se` rather than
        `se_adj`. There is no raw/adjusted pair to choose between.

    Anything else -- no candidate, or more than one after excluding dPSI_raw --
    is refused rather than resolved. Substituting an unadjusted effect for an
    adjusted one would manufacture or destroy a null silently.
    """
    cols = list(cols)
    if "dPSI_adj" in cols:
        return "dPSI_adj"
    cand = [c for c in cols if c.startswith("dPSI") and c != "dPSI_raw"]
    if len(cand) == 1:
        return cand[0]
    raise SystemExit(
        f"[fatal] {path}: cannot resolve the effect column.\n"
        f"        candidates after excluding dPSI_raw: {cand or 'none'}\n"
        f"        columns: {cols}\n"
        f"        Refusing rather than picking one.")


def measure(path, et, fdr=0.05):
    """Every count Table 3 and its ladder need, from one event table."""
    d = pd.read_csv(path, sep="\t")
    for c in ("key", "FDR"):
        if c not in d.columns:
            raise SystemExit(f"[fatal] {path}: no '{c}' column "
                             f"(has {list(d.columns)[:8]}...)")
    dc = dpsi_col(d.columns, path)
    if dc != "dPSI_adj":
        d = d.rename(columns={dc: "dPSI_adj"})
    sigf = d.FDR < fdr
    full = dict(tested_full=int(len(d)), full_fdr=int(sigf.sum()),
                full_05=int((sigf & (d.dPSI_adj.abs() >= 0.05)).sum()),
                full_10=int((sigf & (d.dPSI_adj.abs() >= 0.10)).sum()))
    nk = NKEEP[et]
    d["_ek"] = [":".join(str(k).split(":")[:nk]) for k in d["key"]]
    e = d.sort_values("FDR", kind="mergesort").drop_duplicates("_ek", keep="first")
    s = e.FDR < fdr
    s10 = s & (e.dPSI_adj.abs() >= 0.10)
    out = dict(tested=int(len(e)), fdr_only=int(s.sum()),
               ge05=int((s & (e.dPSI_adj.abs() >= 0.05)).sum()),
               ge10=int(s10.sum()),
               genes=int(e.loc[s10, "gene"].nunique()) if "gene" in e.columns else 0)
    out.update(full)
    out["collapse_pct"] = round(100 * (1 - out["tested"] / out["tested_full"]), 2)
    return out


def calibrate(M, L, verbose=True):
    bad = []
    if verbose:
        print("=" * 76)
        print("CALIBRATION — rebuild every frozen 50 bp number from its own source file")
        print("=" * 76)
    sub = M[(M.cohort.astype(str).str.strip() == "MAGNet")
            & (M.contrast.astype(str).str.strip().isin(MAGNET_ROWS))]
    for _, r in sub.iterrows():
        et, path = str(r.event_type).upper(), str(r["file"])
        if not os.path.exists(path):
            bad.append(f"{r.contrast} {et}: the master names {path}, which is absent")
            continue
        g = measure(path, et)
        want = (int(r["tested_exon"]), int(r["sig_exon"]))
        got = (g["tested"], g["ge10"])
        ok = got == want
        if not ok:
            bad.append(f"{r.contrast} {et}: got {got}, master says {want}  ({path})")
        if verbose:
            print(f"  {'ok ' if ok else '✗  '}{str(r.contrast):<11} {et:<5} "
                  f"{got[0]:>7,} / {got[1]:>5,}   master {want[0]:>7,} / {want[1]:>5,}")
    if L is not None:
        if verbose:
            print()
        for _, r in L.iterrows():
            et, path = str(r.event_type).upper(), str(r["file"])
            if not os.path.exists(path):
                bad.append(f"ladder {r.contrast}: {path} absent")
                continue
            g = measure(path, et)
            pairs = [("exon-level_tested", g["tested"]),
                     ("exon-level_FDR only", g["fdr_only"]),
                     ("exon-level_+|dPSI|>=0.05", g["ge05"]),
                     ("exon-level_+|dPSI|>=0.10", g["ge10"]),
                     ("full-coordinate_tested", g["tested_full"]),
                     ("full-coordinate_FDR only", g["full_fdr"]),
                     ("full-coordinate_+|dPSI|>=0.05", g["full_05"]),
                     ("full-coordinate_+|dPSI|>=0.10", g["full_10"])]
            miss = [(c, v, int(r[c])) for c, v in pairs
                    if c in L.columns and v != int(r[c])]
            for c, v, w in miss:
                bad.append(f"ladder {r.contrast} {c}: got {v}, file says {w}")
            if verbose:
                print(f"  {'ok ' if not miss else '✗  '}ladder {str(r.contrast):<20} "
                      + " ".join(f"{v:,}" for _, v in pairs[:4]))
    return bad


def measure125(verbose=True):
    out, bad, skipped = {}, [], []
    if verbose:
        print()
        print("=" * 76)
        print("MEASURING the 125 bp tables")
        print("=" * 76)
        print(f"  {'contrast':<12}{'cls':<6}{'tested':>9}{'FDR<.05':>9}"
              f"{'>=0.05':>8}{'>=0.10':>8}{'genes':>7}   gate")
    for con in ("DCM vs NF", "HCM vs NF", "DCM vs HCM"):
        for et in TYPES:
            p = F125[con].format(ET=et)
            # Only these are needed to WRITE the inputs: the master carries
            # DCM vs NF and HCM vs NF for all five classes, and the ladder
            # carries SE alone. DCM vs HCM's other four classes are printed
            # for information and must never block the deliverable.
            needed = not (con == "DCM vs HCM" and et != "SE")
            if not os.path.exists(p):
                (bad if needed else skipped).append(f"{con} {et}: {p} not found")
                continue
            try:
                g = measure(p, et)
            except SystemExit as e:
                (bad if needed else skipped).append(f"{con} {et}: {e}")
                continue
            out[(con, et)] = g
            mark = ""
            if (con, et) in GATE:
                want = GATE[(con, et)]
                got = (g["tested"], g["ge10"])
                if got == want:
                    mark = "MATCH"
                else:
                    mark = f"✗ want {want}"
                    bad.append(f"{con} {et}: got {got}, record says {want}  ({p})")
            if verbose:
                print(f"  {con:<12}{et:<6}{g['tested']:>9,}{g['fdr_only']:>9,}"
                      f"{g['ge05']:>8,}{g['ge10']:>8,}{g['genes']:>7,}   {mark}")
    for con, want in GATE_FDRONLY.items():
        g = out.get((con, "SE"))
        if g and g["fdr_only"] != want:
            bad.append(f"{con} SE FDR-only: got {g['fdr_only']}, record says {want}")
    if skipped and verbose:
        print("\n  informational only, skipped (nothing in Table 3 needs these):")
        for sk in skipped:
            print("    " + str(sk).replace("\n", " ")[:150])
    return out, bad


def write_inputs(M, L, G, outdir):
    os.makedirs(outdir, exist_ok=True)
    # ---- master: MAGNet rows re-based, every other cohort copied verbatim
    keep, dropped = [], []
    for _, r in M.iterrows():
        coh, con = str(r.cohort).strip(), str(r.contrast).strip()
        if coh != "MAGNet":
            keep.append(r)
            continue
        if con not in MAGNET_ROWS:
            dropped.append(f"{coh} {con} {r.event_type}")
            continue
        g = G[(con, str(r.event_type).upper())]
        r = r.copy()
        for col, val in (("tested", g["tested"]), ("significant", g["ge10"]),
                         ("tested_exon", g["tested"]), ("sig_exon", g["ge10"]),
                         ("tested_full", g["tested_full"]), ("sig_full", g["full_10"]),
                         ("genes_sig", g["genes"]),
                         ("collapse_pct", g["collapse_pct"]),
                         ("file", F125[con].format(ET=str(r.event_type).upper()))):
            if col in M.columns:
                r[col] = val
        keep.append(r)
    mo = pd.DataFrame(keep)[list(M.columns)]
    mo.to_csv(os.path.join(outdir, MASTER), sep="\t", index=False)

    # ---- ladder: the three SE rows
    rows = []
    for name in LADDER_ROWS:
        con = " ".join(name.split()[1:])
        g = G[(con, "SE")]
        rows.append({"contrast": name, "event_type": "SE",
                     "file": F125[con].format(ET="SE"),
                     "exon-level_FDR only": g["fdr_only"],
                     "exon-level_+|dPSI|>=0.05": g["ge05"],
                     "exon-level_+|dPSI|>=0.10": g["ge10"],
                     "full-coordinate_FDR only": g["full_fdr"],
                     "full-coordinate_+|dPSI|>=0.05": g["full_05"],
                     "full-coordinate_+|dPSI|>=0.10": g["full_10"],
                     "exon-level_tested": g["tested"],
                     "full-coordinate_tested": g["tested_full"]})
    cols = list(L.columns) if L is not None else list(rows[0])
    pd.DataFrame(rows)[cols].to_csv(os.path.join(outdir, LADDER), sep="\t", index=False)

    if os.path.exists(TABLE1):
        shutil.copy2(TABLE1, os.path.join(outdir, TABLE1))
    return mo, dropped


def self_test():
    f = []
    fx = pd.DataFrame({
        "key": ["chr1:1:2:3:4:5:6", "chr1:1:2:9:4:5:6", "chr2:1:2:3:4:5:6"],
        "FDR": [0.001, 0.004, 0.002],
        "dPSI_adj": [0.02, 0.40, 0.30],
        "gene": ["A", "A", "B"]})
    fx.to_csv("/tmp/_167.tsv", sep="\t", index=False)
    se = measure("/tmp/_167.tsv", "SE")
    # SE keeps 3 key fields -> rows 1 and 2 collapse; the LOWEST-FDR row wins and
    # it fails |dPSI|, so A is not significant. Filter-first would call it so.
    if (se["tested"], se["ge10"]) != (2, 1):
        f.append(f"SE collapse-first expected (2,1), got {(se['tested'], se['ge10'])}")
    if se["fdr_only"] != 2:
        f.append(f"SE FDR-only expected 2, got {se['fdr_only']}")
    if se["ge05"] != 1:
        f.append(f"SE >=0.05 expected 1, got {se['ge05']}")
    if se["tested_full"] != 3 or se["full_10"] != 2:
        f.append(f"full-coordinate wrong: {se['tested_full']}/{se['full_10']}")
    if se["genes"] != 1:
        f.append(f"genes expected 1, got {se['genes']}")
    mx = measure("/tmp/_167.tsv", "MXE")      # 5 key fields -> nothing collapses
    if (mx["tested"], mx["ge10"]) != (3, 2):
        f.append(f"MXE expected (3,2), got {(mx['tested'], mx['ge10'])}")
    if se["collapse_pct"] != round(100 * (1 - 2 / 3), 2):
        f.append("collapse_pct is wrong")
    # the ladder is a superset of the table: >=0.05 must never be below >=0.10
    if not (se["fdr_only"] >= se["ge05"] >= se["ge10"]):
        f.append("the ladder is not monotone")
    # the effect-column resolver: both real schemas, and both refusal cases
    std = ["key", "FDR", "dPSI_raw", "dPSI_adj", "se_adj"]
    dch = ["key", "FDR", "dPSI_DCMvsHCM", "se"]
    if dpsi_col(std, "x") != "dPSI_adj":
        f.append("resolver did not prefer dPSI_adj on the 20-column schema")
    if dpsi_col(dch, "x") != "dPSI_DCMvsHCM":
        f.append("resolver did not find dPSI_DCMvsHCM on the 12-column schema")
    for cols in (["key", "FDR"], ["key", "FDR", "dPSI_a", "dPSI_b"]):
        try:
            dpsi_col(cols, "x")
            f.append(f"resolver did not refuse {cols}")
        except SystemExit:
            pass
    # a table carrying ONLY dPSI_raw must be refused, never silently used
    try:
        dpsi_col(["key", "FDR", "dPSI_raw"], "x")
        f.append("resolver accepted a raw-only table")
    except SystemExit:
        pass
    # the recorded gates must be internally consistent
    if sum(v[0] for k, v in GATE.items() if k[0] == "DCM vs NF") != 144942:
        f.append("GATE's DCM vs NF tested do not sum to 144,942")
    if sum(v[1] for k, v in GATE.items() if k[0] == "DCM vs NF") != 1294:
        f.append("GATE's DCM vs NF significant do not sum to 1,294")
    if round(813 / 5401 * 100, 1) != 15.1 or round(431 / 1430 * 100, 1) != 30.1:
        f.append("813/5401 and 431/1430 do not give S5c's 15.1% and 30.1%")
    os.remove("/tmp/_167.tsv")
    print("167 self-test")
    for x in f:
        print("  ✗ " + x)
    if not f:
        print("  ✓ collapse-before-filter verified on a fixture that distinguishes")
        print("    the two orders; NKEEP counts key fields; the full-coordinate and")
        print("    exon-level ladders both compute; the ladder is monotone; the")
        print("    recorded 125 bp gates sum to 144,942 / 1,294 and reproduce S5c's")
        print("    15.1% and 30.1%; the effect-column resolver prefers dPSI_adj,")
        print("    finds dPSI_DCMvsHCM, and refuses raw-only and ambiguous tables")
    return 1 if f else 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--self-test", action="store_true")
    ap.add_argument("--outdir", default="t3_rl125")
    a = ap.parse_args()
    if a.self_test:
        return self_test()
    if self_test():
        print("\nREFUSING: 167's own self-test failed.")
        return 3
    print()
    for p in (MASTER, LADDER):
        if not os.path.exists(p):
            print(f"REFUSING: {p} not found")
            return 3
    M = pd.read_csv(MASTER, sep="\t")
    L = pd.read_csv(LADDER, sep="\t")

    bad = calibrate(M, L)
    if bad:
        print("\nREFUSING — calibration failed. The recipe does not rebuild the")
        print("frozen 50 bp inputs, so it must not write 125 bp ones:")
        for b in bad:
            print("    " + b)
        return 3
    print("\n✓ CALIBRATION PASSED — every frozen 50 bp number reproduced exactly")

    G, bad = measure125()
    if bad:
        print("\nREFUSING — the 125 bp measurement disagrees with the record:")
        for b in bad:
            print("    " + b)
        return 3
    print("\n✓ all ten recorded 125 bp values reproduced")

    mo, dropped = write_inputs(M, L, G, a.outdir)
    print()
    print("=" * 76)
    print(f"[write] {a.outdir}/{MASTER}  ({len(mo)} rows)")
    print(f"[write] {a.outdir}/{LADDER}  (3 rows)")
    print(f"[write] {a.outdir}/{TABLE1}  (copied unchanged)")
    print("=" * 76)
    if dropped:
        print(f"\n  {len(dropped)} MAGNet row(s) NOT carried over, because they are")
        print("  50 bp and this file must not contain a stale MAGNet row:")
        for d in dropped:
            print("    " + d)
        print("  87 does not read them (they are not Table 3 contrasts).")
    tot_t = sum(G[("DCM vs NF", e)]["tested"] for e in TYPES)
    tot_s = sum(G[("DCM vs NF", e)]["ge10"] for e in TYPES)
    ht = sum(G[("HCM vs NF", e)]["tested"] for e in TYPES)
    hs = sum(G[("HCM vs NF", e)]["ge10"] for e in TYPES)
    print(f"\n  DCM vs NF all classes : {tot_t:,} tested, {tot_s:,} significant")
    print("                          the manuscript says 144,942 and 1,294")
    print(f"  HCM vs NF all classes : {ht:,} tested, {hs:,} significant")
    print("                          (50 bp was 118,234 and 695 — NOT previously")
    print("                           recorded at 125 bp; this is the first")
    print("                           measurement, so check it before quoting it)")
    print()
    print("  NEXT:")
    print(f"    python3 87_build_table3.py --dir {a.outdir} --out Table3_rl125")
    print()
    print("  Then read the printed table: the MAGNet rows must show")
    print("    DCM vs NF   53,456  5,401  ...  813   5  144,942  1,294")
    print("    HCM vs NF   39,607  1,430  ...  431   5")
    print("    DCM vs HCM  38,273      1  ...    0   —")
    print("  and the Colorado rows must be UNCHANGED from Table3.tsv.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
