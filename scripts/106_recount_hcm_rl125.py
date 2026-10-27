#!/usr/bin/env python3
"""
106_recount_hcm_rl125.py -- collapse-first recount of the HCM vs NF positive
control at 125 bp, calibrated against the frozen 50 bp comparator.

WHY THIS EXISTS
---------------
phase2_diffsplice_v2.py's printed "sig FDR<0.05: N (|adj dPSI|>=0.10: M)" is
NOT the manuscript number. It is computed on the UNCOLLAPSED set. The locked
rule is COLLAPSE-FIRST:

  1. collapse to the event's own coordinates, keeping the LOWEST-FDR row
  2. THEN apply FDR < 0.05 AND |dPSI_adj| >= 0.10

Filter-first inflates the count (991 SE instead of 926 in the discovery
contrast). The printed 464 for SE rl125 is a filter-first number.

THE CALIBRATED COLLAPSE KEY (discovered 14 Oct, do not assume it)
-----------------------------------------------------------------
  SE, RI               -> (chr, exonStart, exonEnd)
  MXE, A3SS, A5SS      -> key[:5]   (first 5 colon-separated fields of `key`)

phase2's `key` is the full rMATS event identity INCLUDING the flanking exon.
For MXE/A3SS/A5SS, exonStart/exonEnd alone are the WRONG key.

This script does not trust that table. It CALIBRATES: for each class it
searches the candidate keys against the frozen 50 bp file and uses only the
key that reproduces the frozen record. If no key reproduces it, the script
refuses to print an rl125 number for that class.

USAGE
-----
  # step 1 -- see what frozen comparators exist and what each key gives
  python 106_recount_hcm_rl125.py --list

  # step 2 -- the recount (runs calibration first, refuses on failure)
  python 106_recount_hcm_rl125.py

  # if your frozen filenames differ from the defaults:
  python 106_recount_hcm_rl125.py --frozen SE=SE_HCM_v2.all_events.tsv \
                                  --frozen MXE=MXE_HCM_v2.all_events.tsv

  # self-test of the collapse logic on synthetic data (no cluster files needed)
  python 106_recount_hcm_rl125.py --self-test

KNOWN FROZEN RECORD (the only one established in writing)
---------------------------------------------------------
  SE_HCM_v2.all_events.tsv : 70003 rows -> 48024 collapsed -> 468 significant
                             n_g1 min 15 max 28 ; n_g2 min 15 max 161
Anything else is discovered by this script, not supplied by me.
"""
import argparse
import glob
import os
import sys

import numpy as np
import pandas as pd

CLASSES = ["SE", "MXE", "A3SS", "A5SS", "RI"]

# the calibrated default, kept here only as the FIRST candidate to try
DEFAULT_KEY = {
    "SE": ("coords", None),
    "RI": ("coords", None),
    "MXE": ("keytrunc", 5),
    "A3SS": ("keytrunc", 5),
    "A5SS": ("keytrunc", 5),
}

# every key we are willing to consider, in the order we try them
CANDIDATES = [
    ("coords", None),      # chr, exonStart, exonEnd
    ("keytrunc", 5),
    ("keytrunc", 3),
    ("keytrunc", 4),
    ("keytrunc", 6),
    ("keyfull", None),     # the whole key column
]

# the frozen record we are allowed to assert. SE only.
FROZEN_RECORD = {"SE": {"rows": 70003, "collapsed": 48024, "sig": 468}}

REQUIRED = {"chr", "exonStart", "exonEnd", "FDR", "dPSI_adj"}

FDR_CUT = 0.05
DPSI_CUT = 0.10


def keyname(spec):
    kind, n = spec
    if kind == "coords":
        return "chr:exonStart:exonEnd"
    if kind == "keyfull":
        return "key (full)"
    return f"key[:{n}]"


def build_key(df, spec):
    kind, n = spec
    if kind == "coords":
        return (df["chr"].astype(str) + ":" + df["exonStart"].astype(str)
                + ":" + df["exonEnd"].astype(str))
    if "key" not in df.columns:
        return None
    k = df["key"].astype(str)
    if kind == "keyfull":
        return k
    return k.str.split(":").str[:n].str.join(":")


def collapse_then_filter(df, spec):
    """The locked rule. Returns (n_collapsed, n_sig, collapsed_frame)."""
    kk = build_key(df, spec)
    if kk is None:
        return None, None, None
    d = df.copy()
    d["_ck"] = kk
    # lowest FDR wins; stable sort so ties resolve by original file order,
    # which is what reproduced 926 and 468.
    d = d.sort_values("FDR", kind="stable").drop_duplicates("_ck", keep="first")
    sig = d[(d.FDR < FDR_CUT) & (d.dPSI_adj.abs() >= DPSI_CUT)]
    return len(d), len(sig), d


def filter_then_collapse(df, spec):
    """The WRONG order, reported only so the gap is visible."""
    kk = build_key(df, spec)
    if kk is None:
        return None
    d = df.copy()
    d["_ck"] = kk
    d = d[(d.FDR < FDR_CUT) & (d.dPSI_adj.abs() >= DPSI_CUT)]
    d = d.sort_values("FDR", kind="stable").drop_duplicates("_ck", keep="first")
    return len(d)


def load(path):
    df = pd.read_csv(path, sep="\t")
    miss = REQUIRED - set(df.columns)
    if miss:
        return None, f"not a phase2 output (missing {sorted(miss)})"
    return df, None


# ---------------------------------------------------------------- self-test
def self_test():
    """Prove collapse-first != filter-first, and that keytrunc works."""
    # two rows share chr:exonStart:exonEnd. One is significant, one is not,
    # and the NON-significant one has the lower FDR. Collapse-first must keep
    # the low-FDR row and therefore count 0. Filter-first keeps the
    # significant one and counts 1.
    rows = [
        # chr, exonStart, exonEnd, key,           FDR,   dPSI_adj
        ("chr1", 100, 200, "chr1:100:200:50:80",  0.001, 0.02),   # sig? no (dPSI)
        ("chr1", 100, 200, "chr1:100:200:60:90",  0.010, 0.50),   # sig? yes
        ("chr2", 300, 400, "chr2:300:400:10:20",  0.001, 0.30),   # sig
    ]
    df = pd.DataFrame(rows, columns=["chr", "exonStart", "exonEnd", "key",
                                     "FDR", "dPSI_adj"])
    nc, ns, _ = collapse_then_filter(df, ("coords", None))
    nf = filter_then_collapse(df, ("coords", None))
    assert nc == 2, nc
    assert ns == 1, f"collapse-first should give 1, got {ns}"
    assert nf == 2, f"filter-first should give 2, got {nf}"

    # keytrunc 5 keeps the two chr1 rows distinct -> both survive collapse
    nc5, ns5, _ = collapse_then_filter(df, ("keytrunc", 5))
    assert nc5 == 3, nc5
    assert ns5 == 2, ns5

    # keytrunc 3 == coords here
    nc3, ns3, _ = collapse_then_filter(df, ("keytrunc", 3))
    assert (nc3, ns3) == (2, 1), (nc3, ns3)

    print("[ok] self-test passed:")
    print("     collapse-first 1 sig vs filter-first 2 sig on the same 3 rows")
    print("     key[:5] 2 sig, key[:3] 1 sig -- the key changes the answer")
    return 0


# ---------------------------------------------------------------- listing
def do_list(args):
    pats = ["*HCM*all_events.tsv", "*HCMvsNF*all_events.tsv"]
    seen, files = set(), []
    for p in pats:
        for f in sorted(glob.glob(p)):
            if f not in seen:
                seen.add(f)
                files.append(f)
    if not files:
        print("[fatal] no *HCM*all_events.tsv in this directory. cd to the "
              "directory holding the phase2 outputs.", file=sys.stderr)
        return 2
    print(f"found {len(files)} candidate file(s)\n")
    for f in files:
        df, err = load(f)
        if err:
            print(f"{f}\n    skipped: {err}\n")
            continue
        ng1 = f"n_g1 {int(df.n_g1.min())}-{int(df.n_g1.max())}" if "n_g1" in df else ""
        ng2 = f"n_g2 {int(df.n_g2.min())}-{int(df.n_g2.max())}" if "n_g2" in df else ""
        print(f"{f}   {len(df)} rows   {ng1}   {ng2}")
        for spec in CANDIDATES:
            nc, ns, _ = collapse_then_filter(df, spec)
            if nc is None:
                continue
            nf = filter_then_collapse(df, spec)
            mark = ""
            for cls, rec in FROZEN_RECORD.items():
                if len(df) == rec["rows"] and nc == rec["collapsed"] and ns == rec["sig"]:
                    mark = f"   <== reproduces the frozen {cls} record"
            print(f"    {keyname(spec):<22} collapsed {nc:>7}  sig {ns:>5}"
                  f"   (filter-first would say {nf}){mark}")
        print()
    return 0


# ---------------------------------------------------------------- recount
def do_recount(args):
    frozen_map = {}
    for item in args.frozen or []:
        if "=" not in item:
            print(f"[fatal] --frozen wants CLASS=path, got {item!r}", file=sys.stderr)
            return 2
        c, p = item.split("=", 1)
        frozen_map[c.strip().upper()] = p.strip()
    for c in CLASSES:
        frozen_map.setdefault(c, f"{c}_HCM_v2.all_events.tsv")

    new_map = {c: args.new_pattern.replace("{C}", c) for c in CLASSES}

    print("=" * 74)
    print("CALIBRATION -- which collapse key reproduces the frozen 50 bp record?")
    print("=" * 74)

    chosen, cal_rows = {}, []
    for c in CLASSES:
        fp = frozen_map[c]
        if not os.path.exists(fp):
            print(f"{c:<6} frozen file not found: {fp}")
            cal_rows.append((c, fp, None, None, None, "missing"))
            continue
        df, err = load(fp)
        if err:
            print(f"{c:<6} {fp}: {err}")
            cal_rows.append((c, fp, None, None, None, err))
            continue

        rec = FROZEN_RECORD.get(c)
        picked, results = None, []
        for spec in CANDIDATES:
            nc, ns, _ = collapse_then_filter(df, spec)
            if nc is None:
                continue
            results.append((spec, nc, ns))
            if rec and len(df) == rec["rows"] and nc == rec["collapsed"] and ns == rec["sig"]:
                picked = spec
                break

        if picked is not None:
            chosen[c] = picked
            print(f"{c:<6} {fp}  {len(df)} rows  -> key {keyname(picked)}  "
                  f"collapsed {rec['collapsed']}  sig {rec['sig']}   [MATCHES FROZEN RECORD]")
            cal_rows.append((c, fp, len(df), rec["collapsed"], rec["sig"], "matched"))
        else:
            # no asserted record for this class -> fall back to the
            # documented default key, but say so loudly.
            spec = DEFAULT_KEY[c]
            nc, ns, _ = collapse_then_filter(df, spec)
            chosen[c] = spec
            why = ("no frozen record asserted for this class"
                   if not rec else "FROZEN RECORD NOT REPRODUCED")
            print(f"{c:<6} {fp}  {len(df)} rows  -> key {keyname(spec)} "
                  f"(documented default; {why})")
            print(f"       collapsed {nc}  sig {ns}")
            if rec:
                print(f"       expected rows {rec['rows']} collapsed {rec['collapsed']} "
                      f"sig {rec['sig']}")
                for spec2, nc2, ns2 in results:
                    print(f"         tried {keyname(spec2):<22} -> {nc2} / {ns2}")
                print(f"[fatal] {c}: the frozen record did not reproduce. Do not "
                      f"use any {c} number until this is settled.", file=sys.stderr)
                return 3
            cal_rows.append((c, fp, len(df), nc, ns, "default key, unasserted"))

    print()
    print("=" * 74)
    print("RECOUNT -- HCM vs NF at 125 bp, collapse-first")
    print("=" * 74)
    print(f"{'class':<6}{'key':<24}{'rows':>8}{'collapsed':>11}{'sig':>7}"
          f"{'(filter-1st)':>14}")

    out, tot_c, tot_s = [], 0, 0
    for c in CLASSES:
        np_ = new_map[c]
        if not os.path.exists(np_):
            print(f"{c:<6}{'--':<24}{'missing':>8}   {np_}")
            continue
        df, err = load(np_)
        if err:
            print(f"{c:<6}{'--':<24}{'bad':>8}   {err}")
            continue
        spec = chosen.get(c, DEFAULT_KEY[c])
        nc, ns, coll = collapse_then_filter(df, spec)
        nf = filter_then_collapse(df, spec)
        print(f"{c:<6}{keyname(spec):<24}{len(df):>8}{nc:>11}{ns:>7}{nf:>14}")
        out.append((c, len(df), nc, ns, nf, coll))
        tot_c += nc
        tot_s += ns
    print(f"{'TOTAL':<30}{'':>8}{tot_c:>11}{tot_s:>7}")

    # --- comparison against the one frozen number we can assert
    print()
    print("-" * 74)
    print("COMPARISON (only the SE row is against an asserted frozen number)")
    print("-" * 74)
    for c, nrow, nc, ns, nf, _ in out:
        rec = FROZEN_RECORD.get(c)
        if rec:
            d_c = 100.0 * (nc - rec["collapsed"]) / rec["collapsed"]
            d_s = 100.0 * (ns - rec["sig"]) / rec["sig"]
            print(f"{c:<6} tested {rec['collapsed']:>7} -> {nc:>7} ({d_c:+6.1f}%)   "
                  f"sig {rec['sig']:>5} -> {ns:>5} ({d_s:+6.1f}%)")
        else:
            print(f"{c:<6} tested      ? -> {nc:>7}            "
                  f"sig     ? -> {ns:>5}   (no asserted 50 bp number; run "
                  f"--list on the frozen file)")

    # --- where the significant events sit relative to the |dPSI| cut
    print()
    print("-" * 74)
    print("PROXIMITY TO THE |dPSI| >= 0.10 CUT (collapsed, FDR<0.05 only)")
    print("-" * 74)
    for c, nrow, nc, ns, nf, coll in out:
        sub = coll[coll.FDR < FDR_CUT]
        a = sub.dPSI_adj.abs()
        print(f"{c:<6} FDR<0.05 {len(sub):>5}   "
              f"|dPSI| in [0.095,0.105): {int(((a >= 0.095) & (a < 0.105)).sum()):>4}   "
              f"[0.08,0.12): {int(((a >= 0.08) & (a < 0.12)).sum()):>4}")

    if args.write:
        for c, nrow, nc, ns, nf, coll in out:
            p = f"{c}_HCMvsNF_rl125.collapsed.tsv"
            coll.drop(columns=["_ck"]).to_csv(p, sep="\t", index=False)
            sig = coll[(coll.FDR < FDR_CUT) & (coll.dPSI_adj.abs() >= DPSI_CUT)]
            ps = f"{c}_HCMvsNF_rl125.sig.tsv"
            sig.drop(columns=["_ck"]).to_csv(ps, sep="\t", index=False)
            print(f"[ok] wrote {p} ({nc}) and {ps} ({ns})")

    print()
    print("NOTE: the numbers phase2 printed to your terminal "
          "(SE 'sig FDR<0.05: 1528 (|adj dPSI|>=0.10: 464)') are FILTER-FIRST "
          "and are not the manuscript numbers. Use the 'sig' column above.")
    return 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--list", action="store_true",
                    help="inventory frozen comparators x candidate keys")
    ap.add_argument("--self-test", action="store_true",
                    help="prove the collapse logic on synthetic rows")
    ap.add_argument("--frozen", action="append", metavar="CLASS=PATH",
                    help="override a frozen 50 bp comparator path")
    ap.add_argument("--new-pattern", default="{C}_HCMvsNF_rl125.all_events.tsv")
    ap.add_argument("--write", action="store_true",
                    help="write per-class .collapsed.tsv and .sig.tsv")
    a = ap.parse_args()

    if a.self_test:
        return self_test()
    if a.list:
        return do_list(a)
    # always run the self-test before touching real data
    self_test()
    print()
    return do_recount(a)


if __name__ == "__main__":
    sys.exit(main())
