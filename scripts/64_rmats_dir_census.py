#!/usr/bin/env python3
"""
64_rmats_dir_census.py -- what contrast does each rMATS output directory actually encode?

WHY THIS EXISTS
---------------
`63_s5a_jc_vs_jcec.py` reported 102 significant skipped exons for `rmats_magnet` under
JC. `probe_rl50/` reports 172 -- and `probe_rl50` is recorded as THIRTY DCM against THIRTY
NF drawn from that same discovery cohort, quantified at the same declared read length.

A thirty-against-thirty subset cannot out-detect the cohort it was drawn from. One of
these is not what its name says, and until that is settled the "MAGNet internal control"
row of S5(a) is not usable -- it is the row the other four are read against.

The most likely explanation, and the one this script tests first: `rmats_magnet/` may not
be a DCM-versus-NF contrast at all. If that directory exists to hand phase 2 a PSI matrix
for the whole cohort, its own `b1`/`b2` grouping could be arbitrary, or lump disease
groups together, in which case its rMATS-internal FDR tests nothing the manuscript cares
about and its 102 is not a number to anchor a supplement on.

WHAT IT MEASURES, AND WHY EACH ONE ANSWERS THE QUESTION
-------------------------------------------------------
* **Samples per group**, counted from the comma-separated `IncLevel1` / `IncLevel2`
  fields. This is the direct test: if `rmats_magnet` is 30 against 30 it is the probe's
  contrast and the count discrepancy is a parameter difference; if it is 166 against 166
  it is the full cohort and 102 < 172 is impossible as stated; if the groups are lopsided
  or implausible, the grouping is not the manuscript's.
* **Unique `SkipFormLen`**, which for skipped exons is a function of the declared read
  length and nothing else. Two directories sharing it were quantified with the same
  `--readLength`. No formula is assumed -- the values are compared, not decoded.
* **Row counts and significant counts**, so the multiple-testing burden is visible. A
  directory testing 170,000 events carries a harsher FDR correction than one testing
  20,000, which can by itself invert which run looks more powered.
* **Event-set overlap with the probe**, and what becomes of the probe's 172 significant
  events inside `rmats_magnet`. If they are absent from its tested set, the two runs used
  different annotations. If they are present but not significant, it is a power or
  grouping difference.

This script reads. It computes nothing new and writes no table that could reach the
manuscript.

USAGE
-----
    conda activate rnaseq_tools
    cd /data/project/chureeratp/RNA/HEART
    python 64_rmats_dir_census.py
    python 64_rmats_dir_census.py --event MXE      # same census, another class
"""
import argparse
import glob
import os

import numpy as np
import pandas as pd

COORD = {"SE":   ["chr", "exonStart_0base", "exonEnd"],
         "MXE":  ["chr", "1stExonStart_0base", "1stExonEnd",
                  "2ndExonStart_0base", "2ndExonEnd"],
         "A3SS": ["chr", "longExonStart_0base", "longExonEnd", "shortES", "shortEE"],
         "A5SS": ["chr", "longExonStart_0base", "longExonEnd", "shortES", "shortEE"],
         "RI":   ["chr", "riExonStart_0base", "riExonEnd"]}
FDR_CUT, DPSI_CUT = 0.05, 0.10
PROBE = "probe_rl50"
PROBE_RECORD = 172          # significant SE in probe_rl50, confirmed by 62 --parts b
PROBE_DESIGN = "30 DCM vs 30 NF, per RESULT_magnet_readlength_discrepancy.md"


def ngroup(series):
    """Samples in a group, from rMATS' comma-separated per-sample PSI field."""
    vals = series.dropna().astype(str)
    if not len(vals):
        return None
    # a few rows can be truncated; the modal width is the sample count
    widths = vals.head(2000).str.count(",") + 1
    return int(widths.mode().iloc[0]) if len(widths) else None


def census(path, et):
    try:
        t = pd.read_csv(path, sep="\t", low_memory=False)
    except Exception as e:                                      # noqa: BLE001
        return None, f"unreadable ({e.__class__.__name__})"
    out = dict(rows=len(t))
    out["n1"] = ngroup(t["IncLevel1"]) if "IncLevel1" in t.columns else None
    out["n2"] = ngroup(t["IncLevel2"]) if "IncLevel2" in t.columns else None
    if "SkipFormLen" in t.columns:
        u = sorted(pd.to_numeric(t["SkipFormLen"], errors="coerce").dropna().unique())
        out["skiplen"] = u
    else:
        out["skiplen"] = []
    need = COORD[et] + ["FDR", "IncLevelDifference"]
    if any(c not in t.columns for c in need):
        out["sig"] = out["tested"] = None
        out["keys"] = out["sigkeys"] = set()
        return out, None
    t["ekey"] = t[COORD[et]].astype(str).agg(":".join, axis=1)
    t["FDR"] = pd.to_numeric(t["FDR"], errors="coerce")
    t["dpsi"] = pd.to_numeric(t["IncLevelDifference"], errors="coerce")
    t = t[np.isfinite(t["FDR"]) & np.isfinite(t["dpsi"])]
    t = t.loc[t.groupby("ekey")["FDR"].idxmin()]
    sig = t[(t["FDR"] < FDR_CUT) & (t["dpsi"].abs() >= DPSI_CUT)]
    out["tested"] = len(t)
    out["sig"] = len(sig)
    out["keys"] = set(t["ekey"])
    out["sigkeys"] = set(sig["ekey"])
    return out, None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--event", default="SE")
    a = ap.parse_args()
    et = a.event.strip().upper()
    if et not in COORD:
        raise SystemExit(f"unknown event class {et}; choose from {sorted(COORD)}")

    dirs = sorted({os.path.dirname(p) for p in
                   glob.glob(f"**/{et}.MATS.JC.txt", recursive=True)})
    print("=" * 94)
    print(f"rMATS directory census -- {et}")
    print("=" * 94)
    if not dirs:
        raise SystemExit(f"no */{et}.MATS.JC.txt found under {os.getcwd()}")

    res = {}
    print(f"  {'directory':<30}{'grp1':>6}{'grp2':>6}{'rows':>9}{'tested':>9}"
          f"{'signif':>8}   SkipFormLen")
    for d in dirs:
        c, err = census(os.path.join(d, f"{et}.MATS.JC.txt"), et)
        if c is None:
            print(f"  {d:<30} {err}")
            continue
        res[d] = c
        sl = (",".join(str(int(x)) for x in c["skiplen"][:3])
              + (f" (+{len(c['skiplen']) - 3} more)" if len(c["skiplen"]) > 3 else ""))
        print(f"  {d:<30}{str(c['n1'] or '?'):>6}{str(c['n2'] or '?'):>6}"
              f"{c['rows']:>9}{str(c['tested'] or '?'):>9}{str(c['sig'] or '?'):>8}"
              f"   {sl or '?'}")

    # ------------------------------------------------------------- the question
    probe, magnet = res.get(PROBE), res.get("rmats_magnet")
    print("\n" + "=" * 94)
    print("THE DISCREPANCY")
    print("=" * 94)
    if not (probe and magnet):
        print(f"  {PROBE} or rmats_magnet not present; cannot resolve.")
        return
    print(f"  {PROBE:<16} {probe['sig']:>6} significant   "
          f"{probe['n1']} vs {probe['n2']} samples   recorded design: {PROBE_DESIGN}")
    print(f"  {'rmats_magnet':<16} {magnet['sig']:>6} significant   "
          f"{magnet['n1']} vs {magnet['n2']} samples")
    if probe["sig"] != PROBE_RECORD:
        print(f"\n  [warn] {PROBE} gives {probe['sig']}, not the recorded "
              f"{PROBE_RECORD} -- resolve that first")

    # 1. is rmats_magnet the probe's contrast, the full cohort, or neither?
    pn = (probe["n1"] or 0) + (probe["n2"] or 0)
    mn = (magnet["n1"] or 0) + (magnet["n2"] or 0)
    print("\n  1. Group sizes")
    if magnet["n1"] == probe["n1"] and magnet["n2"] == probe["n2"]:
        print("     Same group sizes as the probe. So this is the probe's contrast and the")
        print("     count difference is a PARAMETER difference, not a sample difference --")
        print("     compare SkipFormLen below.")
    elif mn > pn:
        print(f"     rmats_magnet is the larger run ({mn} samples against the probe's"
              f" {pn}).")
        print(f"     A larger run finding FEWER significant events ({magnet['sig']} vs"
              f" {probe['sig']}) is not")
        print("     possible for the same contrast at the same parameters, so the contrast")
        print("     or the parameters differ. Check the grouping that produced b1/b2.")
    else:
        print(f"     rmats_magnet is the SMALLER run ({mn} samples against the probe's"
              f" {pn}).")
    if magnet["n1"] and magnet["n2"]:
        r = magnet["n1"] / magnet["n2"]
        if r > 3 or r < 1 / 3:
            print(f"     Groups are lopsided ({magnet['n1']}:{magnet['n2']}). A"
                  f" disease-vs-control contrast in this")
            print("     project is not that unbalanced -- this grouping is probably not the")
            print("     manuscript's.")

    # 2. same declared read length?
    print("\n  2. Declared read length (via SkipFormLen, which depends on it alone)")
    ps, ms = probe["skiplen"], magnet["skiplen"]
    if ps and ms and set(ps) == set(ms):
        print(f"     Identical ({', '.join(str(int(x)) for x in ps[:3])}). Same"
              f" --readLength, so read length does NOT")
        print("     explain the count difference.")
    elif ps and ms:
        print(f"     DIFFERENT: probe {[int(x) for x in ps[:3]]} vs rmats_magnet"
              f" {[int(x) for x in ms[:3]]}.")
        print("     The two runs used different --readLength. That alone can shift counts,")
        print("     and it means rmats_magnet is not the 50 bp arm the probe is.")
    else:
        print("     SkipFormLen unavailable in one or both tables.")

    # 3. multiple-testing burden
    print("\n  3. Multiple-testing burden")
    if probe["tested"] and magnet["tested"]:
        print(f"     tested: probe {probe['tested']:,} vs rmats_magnet"
              f" {magnet['tested']:,}"
              f"  ({magnet['tested'] / probe['tested']:.1f}x)")
        if magnet["tested"] > 2 * probe["tested"]:
            print("     rmats_magnet corrects across several times more events. Some of the")
            print("     gap is this, and it is a real effect -- but it is a reason the two")
            print("     numbers are not comparable, not a reason either is wrong.")

    # 4. do the probe's significant events even exist in rmats_magnet?
    print("\n  4. The probe's significant events, inside rmats_magnet")
    pk = probe["sigkeys"]
    if not pk:
        print("     none to trace")
    else:
        in_tested = pk & magnet["keys"]
        in_sig = pk & magnet["sigkeys"]
        print(f"     {len(in_tested)}/{len(pk)} are in rmats_magnet's tested set; "
              f"{len(in_sig)} are significant there")
        if len(in_tested) < 0.8 * len(pk):
            print("     >>> Most are ABSENT from its tested set. The two runs did not")
            print("         quantify the same events -- different annotation, or different")
            print("         filtering. They are not two views of one cohort.")
        elif len(in_sig) < 0.5 * len(in_tested):
            print("     >>> Present but mostly NOT significant there. So this is a power or")
            print("         grouping difference, not an annotation difference -- which points")
            print("         at the b1/b2 grouping, given the group sizes above.")
        else:
            print("     Largely preserved; the discrepancy is in the other direction and")
            print("     needs explaining from the counts above.")

    print("\n" + "=" * 94)
    print("WHAT TO DO WITH S5(a)")
    print("=" * 94)
    print("  The JC-versus-JCEC comparison is PAIRED inside each directory, so whatever")
    print("  rmats_magnet's contrast turns out to be, its JC/JCEC ratio is still a")
    print("  like-for-like comparison of two quantifications on one sample set. What is at")
    print("  risk is the INTERPRETATION of that row as 'the discovery cohort's internal")
    print("  control': if the grouping is not DCM versus NF, its significant events are")
    print("  mostly noise and the ratio describes how JCEC treats noise, not signal.")
    print()
    print(f"  {PROBE}/ has both a JC and a JCEC table and a documented design")
    print(f"  ({PROBE_DESIGN}).")
    print("  If rmats_magnet's grouping cannot be confirmed, use the probe as the internal")
    print("  control instead: it is a real DCM-versus-NF contrast with a known sample set,")
    print("  and 62_build_s5.py already anchors on it.")


if __name__ == "__main__":
    main()
