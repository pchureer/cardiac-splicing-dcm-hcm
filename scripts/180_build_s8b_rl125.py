#!/usr/bin/env python3
"""
180_build_s8b_rl125.py -- Supplementary Table 8(b), reverse replication, at either
read length, from the symmetry tables.

67_build_s8.py builds S8b from symmetry_SE.tsv / symmetry_MXE.tsv with no read-length
parameter, gates it against a 50 bp record, and prints a 1.41x "synthetic null" that
is a hard-coded constant from a run that no longer exists on disk. Rather than patch
67 blind, this script re-implements 67's part (b) exactly (same columns, same
rounding, same row rule, same writer) so that it can be CALIBRATED:

  python 180_build_s8b_rl125.py                 # 50 bp: must reproduce the 6 Oct file
  cmp S8b_reverse_replication.tsv S8b_reverse_replication_50check.tsv
  TBL_SUFFIX=_rl125 python 180_build_s8b_rl125.py   # -> S8b_reverse_replication_rl125.tsv

At 50 bp it writes to *_50check.tsv so the 6 Oct deposit is never overwritten, and
refuses to proceed to a 125 bp build in the same invocation. Under a suffix it gates
the reverse slopes against the record established on 20-22 Oct (project docs
RESULT_22_symmetry_rl125, RESULT_colorado_mxe_reverse_rl125), which were written
from the same symmetry files but a week earlier and independently of this script.

The calibration line (forward/reverse against 1.41x) is NOT printed. See the caption
provenance note: it cannot be carried to 125 bp without re-running the synthetic.
S8a (permutations, from 74) is untouched.
"""

import os
import sys

import numpy as np
import pandas as pd

SUF = os.environ.get("TBL_SUFFIX", "")
SYM = "symmetry_{}" + SUF + ".tsv"
SYM_COLS = dict(cohort="cohort", n_common="n_common",
                fwd_n="fwd_n", fwd="fwd_deming", fwd_r="fwd_r",
                rev_n="rev_n", rev="rev_deming", rev_r="rev_r",
                uns_n="unsel_n", uns="unsel_deming", uns_r="unsel_r",
                pre_n="prec_n", pre="prec_deming", pre_r="prec_r")
CLASSES = ["SE", "MXE"]

# (cohort as it appears in the file, class, rev_n, rev slope) -- the record
RECORD = {
    "": [("MAGNet internal", "SE", 432, 0.695), ("Colorado", "SE", 453, 0.447),
         ("Fuwai", "SE", 156, 0.405), ("MAGNet internal", "MXE", 147, 0.725),
         ("Colorado", "MXE", 1477, 0.061), ("Fuwai", "MXE", 99, 0.281)],
    "_rl125": [("MAGNet internal", "SE", 417, 0.714), ("Colorado", "SE", 460, 0.443),
               ("Fuwai", "SE", 156, 0.413), ("MAGNet internal", "MXE", 146, 0.730),
               ("Colorado", "MXE", 1478, 0.060), ("Fuwai", "MXE", 99, 0.303)],
}
DEPOSIT_50 = "S8b_reverse_replication.tsv"


def die(msg):
    print("\n[REFUSE] " + msg + "\n", file=sys.stderr)
    sys.exit(2)


def main():
    if SUF not in RECORD:
        die(f"TBL_SUFFIX={SUF!r}: no record for this suffix.")
    out = (f"S8b_reverse_replication{SUF}.tsv" if SUF
           else "S8b_reverse_replication_50check.tsv")

    have = {}
    for cls in CLASSES:
        path = SYM.format(cls)
        if not (os.path.exists(path) and os.path.getsize(path)):
            die(f"{path} missing")
        t = pd.read_csv(path, sep="\t")
        miss = [c for c in SYM_COLS.values() if c not in t.columns]
        if miss:
            die(f"{path}: pinned columns absent: {miss}")
        have[cls] = t

    # ---- part (b) rows, exactly as 67 builds them
    rows = []
    print(f"\n  {'cohort':<18}{'class':<6}{'rev n':>7}{'rev slope':>11}{'rev r':>8}"
          f"{'fwd n':>7}{'fwd slope':>11}{'UNSEL slope':>13}{'UNSEL r':>9}"
          f"{'PREC slope':>12}{'PREC r':>8}")
    for cls in CLASSES:
        for _, r in have[cls].iterrows():
            c = str(r[SYM_COLS["cohort"]])
            g = lambda k: pd.to_numeric(r[SYM_COLS[k]], errors="coerce")
            if not np.isfinite(g("rev")):
                continue
            print(f"  {c:<18}{cls:<6}{int(g('rev_n')):>7}{g('rev'):>11.3f}"
                  f"{g('rev_r'):>8.3f}{int(g('fwd_n')):>7}{g('fwd'):>11.3f}"
                  f"{g('uns'):>13.3f}{g('uns_r'):>9.3f}{g('pre'):>12.3f}{g('pre_r'):>8.3f}")
            rows.append(dict(part="b. reverse replication", cohort=c, event_class=cls,
                             rev_n=int(g("rev_n")), rev_slope=round(float(g("rev")), 3),
                             rev_r=round(float(g("rev_r")), 3),
                             fwd_n=int(g("fwd_n")) if np.isfinite(g("fwd_n")) else "",
                             fwd_slope=round(float(g("fwd")), 3)
                                       if np.isfinite(g("fwd")) else "",
                             unsel_slope=round(float(g("uns")), 3),
                             unsel_r=round(float(g("uns_r")), 3),
                             prec_slope=round(float(g("pre")), 3),
                             prec_r=round(float(g("pre_r")), 3)))

    # ---- gate against the record (blocks the write; 67's version only noted)
    print()
    bad = []
    for coh, cls, n_, sl in RECORD[SUF]:
        m = have[cls][have[cls][SYM_COLS["cohort"]].astype(str) == coh]
        if not len(m):
            bad.append(f"{coh} {cls}: no row"); continue
        gn = int(pd.to_numeric(m.iloc[0][SYM_COLS["rev_n"]], errors="coerce"))
        gs = float(pd.to_numeric(m.iloc[0][SYM_COLS["rev"]], errors="coerce"))
        good = gn == n_ and abs(gs - sl) < 0.001
        print(f"  [{'ok' if good else 'MISMATCH':<8}] {coh} {cls}: n {gn} slope {gs:.3f}"
              f" (record n {n_} slope {sl})")
        if not good:
            bad.append(f"{coh} {cls}: n {gn}/{gs:.3f} vs record {n_}/{sl}")
    if bad:
        die("reverse slopes do not reproduce the record; nothing written:\n  "
            + "\n  ".join(bad))

    # ---- the two lines that need no calibration
    def pick(cls, coh, k):
        m = have[cls][have[cls][SYM_COLS["cohort"]].astype(str) == coh]
        return float(pd.to_numeric(m.iloc[0][SYM_COLS[k]], errors="coerce"))
    se_rev, mx_rev = pick("SE", "Colorado", "rev"), pick("MXE", "Colorado", "rev")
    print(f"\n  1. within cohort, within direction: Colorado SE {se_rev:.3f} against"
          f" Colorado MXE {mx_rev:.3f} ({se_rev / mx_rev:.1f}x); precision term cancels")
    for cls in CLASSES:
        rs = [pick(cls, "Colorado", k) for k in ("rev_r", "uns_r", "pre_r")]
        print(f"  2. correlation under every estimator: Colorado {cls:<4} r ="
              f" {min(rs):.3f} to {max(rs):.3f}")
    print(f"     Colorado MXE UNSEL/PREC Deming = {pick('MXE','Colorado','uns'):.2f} /"
          f" {pick('MXE','Colorado','pre'):.2f} -- Deming diverging, not attenuation")
    print("  (the forward/reverse asymmetry against a synthetic 1.41x null is NOT "
          "reported: that constant was calibrated at 50 bp and has no producer on disk)")

    T = pd.DataFrame(rows)
    T.to_csv(out, sep="\t", index=False)
    print(f"\n[write] {out}  ({len(T)} rows)")
    if not SUF:
        print(f"\nNOW:  cmp {DEPOSIT_50} {out} && echo IDENTICAL\n"
              f"Only if IDENTICAL:  TBL_SUFFIX=_rl125 python {os.path.basename(__file__)}")


if __name__ == "__main__":
    main()
