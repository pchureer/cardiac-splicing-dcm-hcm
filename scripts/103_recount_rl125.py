#!/usr/bin/env python3
"""
103_recount_rl125.py (v2) — collapse-first recount of the 125 bp re-quantification,
with the exon-level key DETERMINED BY CALIBRATION rather than assumed.

WHY v2
  v1 assumed rMATS column names (exonStart_0base, longExonStart_0base, ...).
  phase2_diffsplice_v2.py writes its own schema:
      GeneID, gene, event_type, chr, strand, exonStart, exonEnd, key, n, n_g1,
      n_g2, mean_psi, sd_psi, dPSI_raw, dPSI_adj, se_adj, dof, pvalue, mde80, FDR
  v1's G2 caught that and refused, which is what it was for.

THE POINT OF THIS SCRIPT
  The collapse key is not guessed. Each candidate key is tried on the FROZEN 50 bp
  tables, and only a key that reproduces the recorded numbers is used on the 125 bp
  tables. If no candidate reproduces them, the script refuses — because then any
  125-vs-50 difference it reported would be its own bug rather than an effect of
  the read length. This is the same calibration discipline as 101's G2.

⚠ READ-ONLY. Writes `recount_rl125.tsv` and nothing else.

USAGE
    python 103_recount_rl125.py --list
        Discovery. Finds candidate adjusted tables, prints their schema and sample
        key values, and reports which candidate key reproduces the record on the
        50 bp arm. Runs no recount and writes nothing.

    python 103_recount_rl125.py --old <glob> --new-suffix _rl125
        The recount, once --list has shown which files are the frozen ones.
        e.g. --old 'SE_dedup.all_events.tsv' style patterns, one per class, given as
             --old-suffix _dedup
"""

import argparse
import glob
import os
import sys

import pandas as pd

CLASSES = ["SE", "MXE", "A3SS", "A5SS", "RI"]

# The frozen record: (tested exon-level, significant exon-level), collapse-first,
# covariate-adjusted. The manuscript's numbers.
REC = {"SE": (62_971, 926), "MXE": (62_762, 231), "A3SS": (15_292, 90),
       "A5SS": (8_988, 84), "RI": (6_102, 114)}
REC_TOTAL = (156_115, 1_445)

FDR_CUT, DPSI_CUT = 0.05, 0.10
MOVE_LIMIT, LEAD_GENE_LIMIT = 10.0, 2

# Candidate exon-level keys, in the order they will be tried.
# A spec is either a list of columns, or ("keytrunc", n) meaning the first n
# colon-separated fields of the 'key' column.
#
# Why keytrunc exists: phase2's 'key' is the FULL rMATS event identity, e.g. for A3SS
#   chr1:35442905:35444303:35442905:35443028:35448172:35448332
# which is chr plus SIX coordinates — long exon, short exon AND the flanking exon. The
# project's exon-level rule uses 3 key fields for SE and RI and 5 for MXE, A3SS and
# A5SS, so the collapse for those three must drop the flanking pair. Whether it does is
# settled by whether it reproduces the record, not by this comment.
KEYCANDS = [
    ("chr + exonStart + exonEnd",            ["chr", "exonStart", "exonEnd"]),
    ("chr + strand + exonStart + exonEnd",   ["chr", "strand", "exonStart", "exonEnd"]),
    ("GeneID + exonStart + exonEnd",         ["GeneID", "exonStart", "exonEnd"]),
    ("gene + exonStart + exonEnd",           ["gene", "exonStart", "exonEnd"]),
    ("key, first 3 fields",                  ("keytrunc", 3)),
    ("key, first 4 fields",                  ("keytrunc", 4)),
    ("key, first 5 fields",                  ("keytrunc", 5)),
    ("key, first 6 fields",                  ("keytrunc", 6)),
    ("key, first 7 fields",                  ("keytrunc", 7)),
    ("the 'key' column as written",          ["key"]),
]


def hr(c="-"):
    print(c * 78)


def refuse(msg):
    hr("=")
    print("REFUSING — nothing was written.\n")
    print(msg)
    hr("=")
    sys.exit(2)


def pct(new, old):
    return (new - old) / old * 100.0


def load(p, strict=True):
    """Returns (df, dPSI column) or, when strict=False, (None, reason)."""
    try:
        d = pd.read_csv(p, sep="\t", low_memory=False)
    except Exception as exc:
        if strict:
            refuse(f"{p}: unreadable — {exc}")
        return None, f"unreadable: {exc}"
    dcol = next((c for c in ("dPSI_adj", "dPSI_raw") if c in d.columns), None)
    if "FDR" not in d.columns or dcol is None:
        why = ("no FDR column" if "FDR" not in d.columns
               else "no dPSI_adj or dPSI_raw column")
        if strict:
            refuse(f"{p}: {why}.\n  columns: {', '.join(map(str, d.columns))}")
        return None, why
    d["_fdr"] = pd.to_numeric(d["FDR"], errors="coerce")
    d["_dpsi"] = pd.to_numeric(d[dcol], errors="coerce")
    return d, dcol


def describe(spec):
    return f"keytrunc:{spec[1]}" if isinstance(spec, tuple) else " + ".join(spec)


def keyseries(d, spec):
    """Build the collapse key, or None if the spec cannot apply to this table."""
    if isinstance(spec, tuple) and spec[0] == "keytrunc":
        if "key" not in d.columns:
            return None
        n = spec[1]
        s = d["key"].astype(str)
        nf = int(s.head(500).str.count(":").max()) + 1
        if nf <= n:                      # would not truncate: same as the whole key
            return None
        return s.str.split(":").str[:n].str.join(":")
    miss = [c for c in spec if c not in d.columns]
    if miss:
        return None
    return d[list(spec)].astype(str).agg(":".join, axis=1)


def count(d, spec):
    """collapse-first: one lowest-FDR row per key, THEN the rule."""
    k = keyseries(d, spec)
    if k is None:
        return None
    c = d.assign(_k=k).sort_values("_fdr", kind="mergesort").drop_duplicates("_k")
    sig = c[(c._fdr < FDR_CUT) & (c._dpsi.abs() >= DPSI_CUT)]
    return len(c), len(sig)


def discover():
    hr("=")
    print("103 --list — discovery and calibration. Nothing is recounted or written.")
    hr("=")
    allf = sorted(set(glob.glob("*all_events.tsv") + glob.glob("*/*all_events.tsv")))
    top = [f for f in allf if "/" not in f]
    if not top:
        refuse("No *all_events.tsv in this directory.")
    print(f"\n  {len(allf)} tables found, {len(top)} at the top level.")
    print("  Subdirectories (inflate/, sizesweep/, splitctl/, jccheck/) are permutation")
    print("  and sensitivity runs, not the frozen analysis, and are skipped.\n")

    probe = next((f for f in top if load(f, strict=False)[0] is not None), None)
    if probe:
        d, dcol = load(probe, strict=False)
        cols = [c for c in d.columns if not str(c).startswith("_")]
        print(f"  schema of {probe}:\n    {', '.join(map(str, cols))}")
        if "key" in d.columns:
            print("\n  sample 'key' values:")
            for v in d["key"].astype(str).head(3):
                print(f"    {v}")

    hr("=")
    print("CALIBRATION — searching for the file+key that reproduces the record")
    hr("=")
    print("  Only matches are printed. Everything else is counted, not listed.\n")

    found = {}
    for cls in CLASSES:
        rt, rs = REC[cls]
        tok = cls + "_"
        cands = [f for f in top
                 if (os.path.basename(f).startswith(tok) or f"_{cls}_" in f)
                 and "_rl125" not in f]
        tried = skipped = 0
        near = []
        for f in cands:
            d, info = load(f, strict=False)
            if d is None:
                skipped += 1
                continue
            if "event_type" in d.columns and d["event_type"].nunique() > 1:
                d = d[d["event_type"].astype(str).str.upper() == cls]
            for name, key in KEYCANDS:
                r = count(d, key)
                if r is None:
                    continue
                tried += 1
                tt, ss = r
                if tt == rt and ss == rs:
                    print(f"  ✅ {cls:<5} {os.path.basename(f):<34} {name:<34}"
                          f"  tested {tt:,}  sig {ss:,}")
                    found.setdefault(cls, []).append((f, name, tuple(key)))
                elif ss == rs or (rt and abs(tt - rt) / rt < 0.02):
                    near.append((cls, os.path.basename(f), name, tt, ss))
        print(f"     {cls:<5} {len(cands)} file(s), {tried} file x key combination(s) "
              f"tried, {skipped} skipped for schema, "
              f"{len(found.get(cls, []))} match(es)")
        for c, b, n, tt, ss in near[:4]:
            print(f"       near: {b:<34} {n:<30} tested {tt:,} sig {ss:,}")
        print()

    hr("=")
    missing = [c for c in CLASSES if c not in found]
    if missing:
        print(f"NO MATCH for: {', '.join(missing)}")
        print("""
  Do not recount. Without a key+file that reproduces the record, any 125-vs-50
  difference would be this script's bug rather than an effect of the read length.

  SE and RI calibrate on _v2filtered with chr + exonStart + exonEnd, so the
  comparator FILE is settled. What is unresolved is the KEY for the classes whose
  exon-level rule uses five fields rather than three.

  phase2's 'key' column is the full rMATS event identity, flanking exon included.
  Count its fields for an unmatched class:

    head -2 MXE_v2filtered.all_events.tsv  | cut -f8
    head -2 A3SS_v2filtered.all_events.tsv | cut -f8

  A key carrying chr plus six coordinates needs its first five fields to match the
  project's rule. This script tries keytrunc 3, 5, 6 and 7. If none reproduces the
  record, read the collapse out of whichever script produced 231 / 90 / 84 rather
  than guessing further — that script is the authority, not this candidate list.
""")
        sys.exit(2)
    print("CALIBRATED:")
    for cls, v in found.items():
        for f, name, key in v:
            print(f"  {cls:<5} {f:<36} key = {name}")
    amb = {c: v for c, v in found.items() if len(v) > 1}
    if amb:
        print(f"\n  ⚠ more than one match for: {', '.join(amb)}. Prefer the narrowest")
        print("    key (chr + exonStart + exonEnd); if two FILES match, they are the")
        print("    same analysis under two names and either will do — but say which.")
    sfx = sorted({os.path.basename(f).split("_", 1)[1].rsplit(".all_events", 1)[0]
                  for v in found.values() for f, _, _ in v})
    print(f"\n  suffix(es) seen: {', '.join('_' + s for s in sfx)}")
    print("  Next:  python 103_recount_rl125.py --old-suffix _<that> "
          "--key chr,exonStart,exonEnd")
    hr("=")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--list", action="store_true",
                    help="discovery and calibration only; recounts nothing")
    ap.add_argument("--old-suffix", default=None,
                    help="frozen tables are <CLASS><old-suffix>.all_events.tsv")
    ap.add_argument("--new-suffix", default="_rl125")
    ap.add_argument("--key", default=None,
                    help="comma-separated key columns, or keytrunc:N, as --list "
                         "reported them; applies to every class")
    ap.add_argument("--key-per-class", default=None,
                    help="different keys per class, e.g. 'SE=chr,exonStart,exonEnd;"
                         "MXE=keytrunc:5;A3SS=keytrunc:5;A5SS=keytrunc:5;"
                         "RI=chr,exonStart,exonEnd'")
    ap.add_argument("--table2", default="Table2.tsv")
    ap.add_argument("--out", default="recount_rl125.tsv")
    a = ap.parse_args()

    if a.list or (a.old_suffix is None and a.key is None):
        discover()
        if not a.list:
            print("\n  (ran --list because --old-suffix/--key were not given)")
        return

    def parse(spec):
        spec = spec.strip()
        if spec.startswith("keytrunc:"):
            return ("keytrunc", int(spec.split(":", 1)[1]))
        return [c.strip() for c in spec.split(",")]

    key = parse(a.key) if a.key else None
    perclass = {}
    if a.key_per_class:
        for part in a.key_per_class.split(";"):
            if part.strip():
                cls, spec = part.split("=", 1)
                perclass[cls.strip()] = parse(spec)
        gap = [c for c in CLASSES if c not in perclass]
        if gap:
            refuse(f"--key-per-class is missing: {', '.join(gap)}")
    if (key is None and not perclass) or a.old_suffix is None:
        refuse("Pass --old-suffix plus --key or --key-per-class, as --list reported\n"
               "  them. The key must be one that reproduced the record on the 50 bp\n"
               "  arm. SE and RI need 3 key fields and MXE/A3SS/A5SS need 5, so the\n"
               "  two groups may need different keys — that is what --key-per-class\n"
               "  is for.")

    hr("=")
    print("103 — collapse-first recount, calibrated")
    hr("=")
    if perclass:
        for c in CLASSES:
            print(f"  {c:<5} key: {describe(perclass[c])}")
    else:
        print(f"  key: {describe(key)}")

    rows = []
    tot_t = tot_s = 0
    for cls in CLASSES:
        po = f"{cls}{a.old_suffix}.all_events.tsv"
        pn = f"{cls}{a.new_suffix}.all_events.tsv"
        for p in (po, pn):
            if not os.path.exists(p):
                refuse(f"G1 — {p} not found.")
        do, _ = load(po)
        dn, _ = load(pn)
        k = perclass[cls] if perclass else key
        ro, rn = count(do, k), count(dn, k)
        if ro is None or rn is None:
            refuse(f"G2 — {cls}: key {describe(k)} does not resolve in both files.")
        rt, rs = REC[cls]
        if ro != (rt, rs):
            refuse(f"G3 CALIBRATION — {cls}: the 50 bp arm gives "
                   f"tested {ro[0]:,} / sig {ro[1]:,}, record is {rt:,} / {rs:,}.\n\n"
                   f"  The key does not reproduce the record for this class, so any\n"
                   f"  difference reported for it would be this script's bug. Rerun\n"
                   f"  --list and resolve {cls} before reading anything.")
        tot_t += rn[0]; tot_s += rn[1]
        rows.append({"class": cls, "tested_50": ro[0], "tested_125": rn[0],
                     "sig_50": ro[1], "sig_125": rn[1],
                     "pct_change": round(pct(rn[1], ro[1]), 1)})
    print(f"  G3 PASS  the 50 bp arm reproduces the record in all five classes")

    hr("=")
    print(f"{'class':<6}{'tested50':>10}{'tested125':>11}{'sig50':>8}{'sig125':>8}{'change':>9}")
    hr()
    for r in rows:
        print(f"{r['class']:<6}{r['tested_50']:>10,}{r['tested_125']:>11,}"
              f"{r['sig_50']:>8,}{r['sig_125']:>8,}{r['pct_change']:>8.1f}%")
    hr()
    print(f"{'TOTAL':<6}{REC_TOTAL[0]:>10,}{tot_t:>11,}{REC_TOTAL[1]:>8,}{tot_s:>8,}"
          f"{pct(tot_s, REC_TOTAL[1]):>8.1f}%")

    pd.DataFrame(rows).to_csv(a.out, sep="\t", index=False)
    print(f"\n  [write] {a.out}")

    hr("=")
    print("VERDICT against the threshold fixed 13 Oct, before the run finished")
    hr("=")
    se_move = abs(next(r["pct_change"] for r in rows if r["class"] == "SE"))
    tot_move = abs(pct(tot_s, REC_TOTAL[1]))
    breach = []
    if tot_move > MOVE_LIMIT:
        breach.append(f"TOTAL moved {tot_move:.1f}% (> {MOVE_LIMIT}%)")
    if se_move > MOVE_LIMIT:
        breach.append(f"SE moved {se_move:.1f}% (> {MOVE_LIMIT}%)")
    print(f"  TOTAL  {tot_move:>5.1f}%   limit {MOVE_LIMIT}%   "
          f"{'BREACH' if tot_move > MOVE_LIMIT else 'within'}")
    print(f"  SE     {se_move:>5.1f}%   limit {MOVE_LIMIT}%   "
          f"{'BREACH' if se_move > MOVE_LIMIT else 'within'}")
    print("  leads  not tested by this script — check Table 2 by hand, and across all")
    print("         event classes, not SE alone")
    print("\n  ⚠ DCM vs HCM: CANNOT BE TESTED. That contrast is a separate rMATS run and")
    print("    was NOT re-quantified at 125 bp, so promotion is unavailable regardless")
    print("    of the numbers above.")
    hr()
    if breach:
        print("  THRESHOLD BREACHED:")
        for b in breach:
            print("     " + b)
        print("\n  The rule says promote — but not until DCM vs HCM is re-quantified,")
        print("  and never partially: everything selecting on MAGNet significance must")
        print("  come from the new selection set.")
    else:
        print("  WITHIN THRESHOLD on both criteria this script can test.")
        print("\n  The rule says: report as a full-cohort sensitivity panel and keep")
        print("  50 bp primary. No headline number moves, no gate is disarmed, and the")
        print("  pilot goes from 60 samples to 332.")
    hr("=")


if __name__ == "__main__":
    main()
