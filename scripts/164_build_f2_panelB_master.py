#!/usr/bin/env python3
"""
164_build_f2_panelB_master.py -- Figure 2 panel B reads the FROZEN master.

39_figure2_discovery.py is already read-length aware in panels A and C:
`--discovery` and `--mats` both flow through. Panel B does not. It reads
master_numbers_flagged.tsv and takes the MAGNet / "DCM vs NF" rows, which are
50 bp:

    MAGNet  DCM vs NF  SE    exon  reference  62971  926  ...
    MAGNet  DCM vs NF  MXE   exon  reference  62762  231  ...
    MAGNet  DCM vs NF  A3SS  exon  reference  15292   90  ...
    MAGNet  DCM vs NF  A5SS  exon  reference   8988   84  ...
    MAGNet  DCM vs NF  RI    exon  reference   6102  114  ...

39 TAKES `--master` AS A FLAG, so no patch to 39 is needed. What is needed is
a five-row master carrying the 125 bp counts -- and those counts must be
MEASURED FROM THE EVENT TABLES, not transcribed from a notes file. A constant
typed into a TSV by hand is the unsourced-constant failure mode wearing a
different hat.

WHAT THIS DOES
--------------
1. CALIBRATE. Collapse-then-filter the five 50 bp `{ET}_v2filtered` tables and
   require all ten numbers above to reproduce EXACTLY. If the recipe cannot
   rebuild the frozen master it has no business writing a new one.
2. FINGERPRINT. For each class, scan candidate 125 bp files and accept ONLY
   one that reproduces the recorded 125 bp (tested, significant) pair. Names
   are not trusted -- `SE_rl125` and `SE_v2filtered_rl125` are byte-identical
   twins and the other four classes' names are unknown to this script.
3. WRITE master_f2_panelB_rl125.tsv, five rows, same column names 39 reads.

CONVENTION (locked, and the same one 46/61/76/158 use independently):
  collapse to the event's own coordinates FIRST, keeping the lowest-FDR row,
  THEN apply FDR < 0.05 and |dPSI_adj| >= 0.10.
  NKEEP counts FIELDS OF `key`, not rows: SE 3, MXE 5, A3SS 5, A5SS 5, RI 3.

USAGE
-----
  python3 164_build_f2_panelB_master.py --self-test
  python3 164_build_f2_panelB_master.py --discover
  python3 164_build_f2_panelB_master.py            # calibrate, measure, write
"""
import argparse
import glob
import os
import sys

import numpy as np
import pandas as pd

TYPES = ["SE", "MXE", "A3SS", "A5SS", "RI"]
NKEEP = {"SE": 3, "MXE": 5, "A3SS": 5, "A5SS": 5, "RI": 3}

# The frozen master's own MAGNet / "DCM vs NF" rows, as pasted from the file.
# Calibration target: the recipe must reproduce every one of these ten.
RL50 = {"SE": (62971, 926), "MXE": (62762, 231), "A3SS": (15292, 90),
        "A5SS": (8988, 84), "RI": (6102, 114)}

# 125 bp, measured by 103_recount_rl125.py and reproduced independently by 61
# (SE 53,456 / 813) and 76. Used ONLY as a fingerprint to identify the right
# file -- the number written to the TSV is the one this script measures.
RL125 = {"SE": (53456, 813), "MXE": (62537, 229), "A3SS": (14346, 87),
         "A5SS": (8360, 74), "RI": (6243, 91)}

# A candidate must not name another cohort, contrast or model.
DENY = ("colorado", "brigham", "stanford", "fuwai", "lbnl", "magnet_",
        "hcm", "icm", "dcmvshcm", "compadj", "v2compadj", "splitctl",
        "perm", "shuffle", "null", "jc_", "_jc.", "matched", "mp_")


def collapse_filter(df, et, fdr=0.05, mind=0.10):
    """(tested, significant) under the locked rule. Collapse FIRST."""
    d = df.copy()
    nk = NKEEP[et]
    d["_ek"] = [":".join(str(k).split(":")[:nk]) for k in d["key"]]
    d = d.sort_values("FDR", kind="mergesort").drop_duplicates("_ek", keep="first")
    tested = int(len(d))
    sig = int(((d["FDR"] < fdr) & (d["dPSI_adj"].abs() >= mind)).sum())
    return tested, sig


def read(path):
    try:
        t = pd.read_csv(path, sep="\t")
    except Exception as e:
        return None, f"unreadable: {e}"
    for c in ("key", "FDR", "dPSI_adj"):
        if c not in t.columns:
            return None, f"no '{c}' column (has {list(t.columns)[:8]}...)"
    return t, None


def candidates(et):
    seen, out = set(), []
    for pat in (f"{et}_*rl125*.all_events.tsv", f"{et}_*125*.all_events.tsv",
                f"{et}rl125*.all_events.tsv"):
        for p in glob.glob(pat):
            b = os.path.basename(p).lower()
            if p in seen or any(x in b for x in DENY):
                continue
            seen.add(p)
            out.append(p)
    return sorted(out)


def calibrate(verbose=True):
    bad, got = [], {}
    for et in TYPES:
        p = f"{et}_v2filtered.all_events.tsv"
        if not os.path.exists(p):
            bad.append(f"{et}: {p} not found")
            continue
        t, err = read(p)
        if t is None:
            bad.append(f"{et}: {err}")
            continue
        got[et] = collapse_filter(t, et)
        if got[et] != RL50[et]:
            bad.append(f"{et}: got {got[et]}, master says {RL50[et]}  ({p})")
    if verbose:
        print("=" * 74)
        print("CALIBRATION — rebuild the frozen master's five 50 bp rows")
        print("=" * 74)
        for et in TYPES:
            if et in got:
                mark = "ok " if got[et] == RL50[et] else "✗  "
                print(f"  {mark}{et:<5} tested {got[et][0]:>7,}  sig {got[et][1]:>5,}"
                      f"      master {RL50[et][0]:>7,} / {RL50[et][1]:>5,}")
    return bad


def measure125(verbose=True):
    """Fingerprint the 125 bp table for each class. Returns {et: (path, t, s)}."""
    out, bad = {}, []
    if verbose:
        print()
        print("=" * 74)
        print("FINGERPRINTING the 125 bp tables (names are NOT trusted)")
        print("=" * 74)
    for et in TYPES:
        hits = []
        for p in candidates(et):
            t, err = read(p)
            if t is None:
                if verbose:
                    print(f"  {et:<5} {os.path.basename(p):<44} {err[:28]}")
                continue
            ts = collapse_filter(t, et)
            ok = ts == RL125[et]
            if verbose:
                print(f"  {et:<5} {os.path.basename(p):<44} "
                      f"{ts[0]:>7,} / {ts[1]:>5,}  {'<= MATCH' if ok else ''}")
            if ok:
                hits.append((p, ts))
        if not hits:
            bad.append(f"{et}: no candidate reproduces {RL125[et]}")
        else:
            out[et] = hits[0]
            if len(hits) > 1 and verbose:
                print(f"        {len(hits)} files give the identical pair "
                      f"(byte-identical twins); using {os.path.basename(hits[0][0])}")
    return out, bad


def write_master(found, path):
    rows = []
    for et in TYPES:
        p, (tested, sig) = found[et]
        rows.append({"cohort": "MAGNet", "contrast": "DCM vs NF", "event_type": et,
                     "convention": "exon", "reliability": "reference",
                     "tested": tested, "significant": sig,
                     "tested_exon": tested, "sig_exon": sig, "file": os.path.basename(p)})
    df = pd.DataFrame(rows)
    df.to_csv(path, sep="\t", index=False)
    return df


# ----------------------------------------------------------------- self-test
def _fix():
    """Two events. Event A has two rows; its LOWEST-FDR row fails |dPSI| but a
    sibling row passes. Collapse-first keeps the lowest-FDR row and A is NOT
    significant; filter-first would keep the sibling and call A significant.
    The two orders give different answers, so the test can tell them apart."""
    return pd.DataFrame({
        "key":      ["chr1:100:200:300:400:500:600", "chr1:100:200:999:400:500:600",
                     "chr2:100:200:300:400:500:600"],
        "FDR":      [0.001, 0.004, 0.002],
        "dPSI_adj": [0.02, 0.40, 0.30],
    })


def self_test():
    f = []
    # SE keeps 3 key fields -> rows 1 and 2 collapse together
    t, s = collapse_filter(_fix(), "SE")
    if (t, s) != (2, 1):
        f.append(f"SE collapse-first expected (2, 1), got ({t}, {s})")
    # MXE keeps 5 fields -> field 4 differs, so nothing collapses
    t, s = collapse_filter(_fix(), "MXE")
    if (t, s) != (3, 2):
        f.append(f"MXE (5 key fields) expected (3, 2), got ({t}, {s})")
    # the order matters: filter-then-collapse must give a DIFFERENT answer,
    # otherwise this fixture proves nothing
    d = _fix()
    d = d[(d.FDR < 0.05) & (d.dPSI_adj.abs() >= 0.10)]
    d["_ek"] = [":".join(str(k).split(":")[:3]) for k in d["key"]]
    if d["_ek"].nunique() != 2:
        f.append("fixture does not distinguish the two orders")
    # NKEEP counts FIELDS, not characters
    k = "chr1:100:200:300:400:500:600"
    if ":".join(k.split(":")[:3]) != "chr1:100:200":
        f.append("key truncation is wrong")
    # the denylist must reject a Colorado file and keep a plain one
    if not any(x in "se_colorado_dcmvsnf_rl125.all_events.tsv" for x in DENY):
        f.append("denylist would admit a Colorado table")
    if any(x in "se_rl125.all_events.tsv" for x in DENY):
        f.append("denylist rejects the plain SE table")
    # the two recorded columns must be arithmetically consistent
    if sum(v[0] for v in RL125.values()) != 144942:
        f.append("RL125 tested totals do not sum to 144,942")
    if sum(v[1] for v in RL125.values()) != 1294:
        f.append("RL125 significant totals do not sum to 1,294")
    if sum(v[0] for v in RL50.values()) != 156115:
        f.append("RL50 tested totals do not sum to 156,115")
    if sum(v[1] for v in RL50.values()) != 1445:
        f.append("RL50 significant totals do not sum to 1,445")
    print("164 self-test")
    for x in f:
        print("  ✗ " + x)
    if not f:
        print("  ✓ collapse-before-filter verified against a fixture that")
        print("    DISTINGUISHES the two orders, NKEEP counts key fields,")
        print("    the denylist holds, and both recorded columns sum to their")
        print("    published totals (144,942/1,294 and 156,115/1,445)")
    return 1 if f else 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--self-test", action="store_true")
    ap.add_argument("--discover", action="store_true")
    ap.add_argument("--out", default="master_f2_panelB_rl125.tsv")
    a = ap.parse_args()

    if a.self_test:
        return self_test()
    if self_test():
        print("\nREFUSING: 164's own self-test failed.")
        return 3
    print()

    if a.discover:
        for et in TYPES:
            print(f"{et}:")
            for p in candidates(et):
                print(f"    {p}")
            if not candidates(et):
                print("    (none)")
        return 0

    bad = calibrate()
    if bad:
        print("\nREFUSING — calibration failed. The recipe does not rebuild the")
        print("frozen master, so it must not write a new one:")
        for b in bad:
            print("    " + b)
        return 3
    print("\n✓ CALIBRATION PASSED — all ten 50 bp numbers reproduced exactly")

    found, bad = measure125()
    if bad:
        print("\nREFUSING — could not identify a 125 bp table for every class:")
        for b in bad:
            print("    " + b)
        print("\n  Run --discover to see what was considered. Do NOT relax the")
        print("  fingerprint: a file that gives the wrong count IS the wrong file.")
        return 3

    df = write_master(found, a.out)
    print()
    print("=" * 74)
    print(f"[write] {a.out}")
    print("=" * 74)
    print(df[["event_type", "tested_exon", "sig_exon", "file"]].to_string(index=False))
    print(f"\n  totals: {int(df.tested_exon.sum()):,} tested, "
          f"{int(df.sig_exon.sum()):,} significant")
    print("  the manuscript says 144,942 tested and 1,294 significant")
    print()
    print("  THEN REGENERATE FIGURE 2:")
    print()
    print("  python3 39_figure2_discovery.py --event-type SE \\")
    print(f"      --master {a.out} \\")
    print("      --discovery SE_rl125 \\")
    print("      --mats rmats_magnet_rl125/{ET}.MATS.JCEC.txt \\")
    print("      --out F2_discovery_rl125")
    print()
    print("  Watch for, in order:")
    print("    [info] volcano from SE_rl125...  (53,456 exon-level SE events)")
    print("    [info] event-type totals: SE 813, MXE 229, ... total 1,294")
    print("    [info] heatmap rows ranked by |Cohen's d| ... d ranges A to B")
    print("  The last line is the legend's (d 1.01-1.51) and WILL have moved.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
