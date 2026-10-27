#!/usr/bin/env python3
"""
104_why_counts_moved.py — decompose the 125-vs-50 count change into a DENOMINATOR
part (events that stopped being testable) and a NUMERATOR part (events that stayed
testable but stopped being significant).

WHY THIS EXISTS
  103 found the threshold breached: SE -12.2%, total -10.4%. But SE's TESTED set fell
  15.1%, more than its significant count, and the significance RATE rose (1.470% ->
  1.521%). A count can fall because the biology moved or because the denominator moved,
  and those have opposite implications. This script separates them.

  It does NOT re-open the threshold. 103's verdict stands on the counts it was written
  for. This answers a different question: what the breach consists of.

WHAT IT REPORTS, per event class
  rate         significant / tested, before and after
  the 926      of the frozen significant events: still tested? still significant?
               lost significance while still testable? no longer testable at all?
  the dropout  sd_psi and n of events that left the tested set, against those that
               stayed — if the dropouts are low-variance, low-n events, the loss is
               the filter biting, not a change in effect
  direction    sign agreement on events significant in BOTH arms

⚠ READ-ONLY. Writes `why_counts_moved.tsv` and nothing else.

USAGE
  python 104_why_counts_moved.py --old-suffix _v2filtered \
      --key-per-class 'SE=chr,exonStart,exonEnd;MXE=keytrunc:5;A3SS=keytrunc:5;A5SS=keytrunc:5;RI=chr,exonStart,exonEnd'
"""

import argparse
import os
import sys

import pandas as pd

CLASSES = ["SE", "MXE", "A3SS", "A5SS", "RI"]
REC = {"SE": (62_971, 926), "MXE": (62_762, 231), "A3SS": (15_292, 90),
       "A5SS": (8_988, 84), "RI": (6_102, 114)}
FDR_CUT, DPSI_CUT = 0.05, 0.10


def hr(c="-"):
    print(c * 78)


def refuse(m):
    hr("="); print("REFUSING — nothing was written.\n"); print(m); hr("="); sys.exit(2)


def parse(spec):
    spec = spec.strip()
    if spec.startswith("keytrunc:"):
        return ("keytrunc", int(spec.split(":", 1)[1]))
    return [c.strip() for c in spec.split(",")]


def keyseries(d, spec):
    if isinstance(spec, tuple) and spec[0] == "keytrunc":
        if "key" not in d.columns:
            return None
        s = d["key"].astype(str)
        if int(s.head(500).str.count(":").max()) + 1 <= spec[1]:
            return None
        return s.str.split(":").str[:spec[1]].str.join(":")
    if [c for c in spec if c not in d.columns]:
        return None
    return d[list(spec)].astype(str).agg(":".join, axis=1)


def prep(p, spec):
    d = pd.read_csv(p, sep="\t", low_memory=False)
    dcol = next((c for c in ("dPSI_adj", "dPSI_raw") if c in d.columns), None)
    if dcol is None or "FDR" not in d.columns:
        refuse(f"{p}: needs FDR and dPSI_adj/dPSI_raw.\n"
               f"  columns: {', '.join(map(str, d.columns))}")
    d["_fdr"] = pd.to_numeric(d["FDR"], errors="coerce")
    d["_d"] = pd.to_numeric(d[dcol], errors="coerce")
    k = keyseries(d, spec)
    if k is None:
        refuse(f"{p}: key does not resolve.")
    d["_k"] = k
    c = d.sort_values("_fdr", kind="mergesort").drop_duplicates("_k").set_index("_k")
    c["_sig"] = (c._fdr < FDR_CUT) & (c._d.abs() >= DPSI_CUT)
    return c


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--old-suffix", required=True)
    ap.add_argument("--new-suffix", default="_rl125")
    ap.add_argument("--key-per-class", required=True)
    ap.add_argument("--out", default="why_counts_moved.tsv")
    a = ap.parse_args()

    pc = {}
    for part in a.key_per_class.split(";"):
        if part.strip():
            cls, spec = part.split("=", 1)
            pc[cls.strip()] = parse(spec)
    gap = [c for c in CLASSES if c not in pc]
    if gap:
        refuse(f"--key-per-class is missing: {', '.join(gap)}")

    hr("=")
    print("104 — what the count change is made of")
    hr("=")

    rows = []
    for cls in CLASSES:
        po, pn = f"{cls}{a.old_suffix}.all_events.tsv", f"{cls}{a.new_suffix}.all_events.tsv"
        for p in (po, pn):
            if not os.path.exists(p):
                refuse(f"{p} not found.")
        o, n = prep(po, pc[cls]), prep(pn, pc[cls])

        rt, rs = REC[cls]
        if len(o) != rt or int(o._sig.sum()) != rs:
            refuse(f"CALIBRATION — {cls}: the 50 bp arm gives tested {len(o):,} / "
                   f"sig {int(o._sig.sum()):,}, record {rt:,} / {rs:,}.\n"
                   f"  Wrong file or wrong key. Rerun 103 --list.")

        so = set(o.index[o._sig])
        tn = set(n.index)
        kept_tested = so & tn
        still_sig = {k for k in kept_tested if bool(n.at[k, "_sig"])}
        lost_while_tested = kept_tested - still_sig
        not_tested = so - tn

        rate_o = len(so) / len(o) * 100
        rate_n = int(n._sig.sum()) / len(n) * 100

        # dropout profile: events tested at 50 bp but not at 125 bp
        drop = o.loc[list(set(o.index) - tn)]
        stay = o.loc[list(set(o.index) & tn)]
        def med(df, col):
            return float(df[col].median()) if col in df.columns and len(df) else float("nan")

        # direction on events significant in both
        both = list(still_sig)
        agree = sum(1 for k in both
                    if (o.at[k, "_d"] > 0) == (n.at[k, "_d"] > 0)) if both else 0

        print(f"\n  {cls}")
        print(f"    rate          {rate_o:.3f}%  ->  {rate_n:.3f}%   "
              f"({(rate_n-rate_o)/rate_o*100:+.1f}%)")
        print(f"    tested        {len(o):>7,}  ->  {len(n):>7,}   "
              f"({(len(n)-len(o))/len(o)*100:+.1f}%)")
        print(f"    significant   {len(so):>7,}  ->  {int(n._sig.sum()):>7,}   "
              f"({(int(n._sig.sum())-len(so))/len(so)*100:+.1f}%)")
        print(f"    of the {len(so)} significant at 50 bp:")
        print(f"       still significant        {len(still_sig):>6,}")
        print(f"       tested, lost significance {len(lost_while_tested):>5,}")
        print(f"       NO LONGER TESTED          {len(not_tested):>5,}"
              f"   <- denominator loss")
        if both:
            print(f"    direction agreement on the {len(both)} significant in both: "
                  f"{agree}/{len(both)} = {agree/len(both)*100:.2f}%")
        print(f"    dropout profile  median sd_psi  left {med(drop,'sd_psi'):.4f}"
              f"   stayed {med(stay,'sd_psi'):.4f}")
        print(f"                     median n      left {med(drop,'n'):.0f}"
              f"      stayed {med(stay,'n'):.0f}")

        rows.append({"class": cls, "tested_50": len(o), "tested_125": len(n),
                     "sig_50": len(so), "sig_125": int(n._sig.sum()),
                     "rate_50_pct": round(rate_o, 4), "rate_125_pct": round(rate_n, 4),
                     "still_sig": len(still_sig),
                     "tested_lost_sig": len(lost_while_tested),
                     "no_longer_tested": len(not_tested),
                     "dir_agree": agree, "dir_n": len(both),
                     "med_sdpsi_dropped": round(med(drop, "sd_psi"), 5),
                     "med_sdpsi_stayed": round(med(stay, "sd_psi"), 5)})

    df = pd.DataFrame(rows)
    df.to_csv(a.out, sep="\t", index=False)

    hr("=")
    print("DECOMPOSITION — where the lost significant events went")
    hr("=")
    tl = df.sig_50.sum() - df.sig_125.sum()
    nlt = df.no_longer_tested.sum()
    lwt = df.tested_lost_sig.sum()
    print(f"  significant at 50 bp                     {df.sig_50.sum():>6,}")
    print(f"  significant at 125 bp                    {df.sig_125.sum():>6,}")
    print(f"  net change                               {-tl:>+6,}")
    print()
    print(f"  of the 50 bp significant events:")
    print(f"    no longer TESTABLE at 125 bp           {nlt:>6,}"
          f"   ({nlt/df.sig_50.sum()*100:.1f}% of {df.sig_50.sum():,})")
    print(f"    testable but no longer significant     {lwt:>6,}"
          f"   ({lwt/df.sig_50.sum()*100:.1f}%)")
    print(f"    still significant                      {df.still_sig.sum():>6,}")
    print(f"\n  [write] {a.out}")
    hr("=")
    print("""
HOW TO READ THIS

  If most of the loss sits in "no longer TESTABLE", the count fell because the
  denominator fell — the PSI-variance filter biting differently once PSI is
  computed with the correct form lengths. That is a statement about the filter,
  not about the splicing, and the significance RATE is then the comparable
  quantity.

  If most of the loss sits in "testable but no longer significant", effects or
  their standard errors moved, and the count change is real in the sense the
  threshold assumed.

  ⚠ Either way, 103's verdict stands as recorded. This does not reopen it. What
    it changes is the explanation that has to accompany whichever decision is
    taken — and if the loss is a denominator effect, that explanation belongs in
    the paper whether or not the 125 bp quantification is promoted.
""")


if __name__ == "__main__":
    main()
