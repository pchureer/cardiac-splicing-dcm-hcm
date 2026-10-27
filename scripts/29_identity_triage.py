#!/usr/bin/env python3
"""
29_identity_triage.py -- separate the real duplicates from the artefacts in
geno/identity_metrics.tsv. Run after 28. Genotypes nothing; seconds.

WHY THIS EXISTS: 28's DUPLICATE LIST IS WRONG AND THIS FIXES IT
---------------------------------------------------------------
Script 28 reported 99 cross-project duplicate pairs. Only about a quarter of those
survive scrutiny. Two defects in 28, both mine:

1. THE IBS0 PROMOTION RULE WAS MISCALIBRATED. 28 promotes a pair to "duplicate" when
   its IBS0 falls below `ibs0_floor / 4`, where ibs0_floor is the MINIMUM IBS0over all
   216,237 pairs. A minimum over 216k draws of a discrete count is ~0 by construction,
   so the cut landed at 0.00045 and every pair with a literal zero opposite-homozygote
   count and kinship above 0.177 was promoted. 26 pairs exceed the nominal 0.354
   cut-point natively; 132 were reported. The extra ~106 are the rule firing, not data.

2. NO PER-SAMPLE HETEROZYGOSITY QC. A sample with an inflated het rate -- contamination,
   RNA editing, multi-mapping at low depth -- pairs as a relative with EVERYONE: it
   drives N(het,het) up with every partner while almost never being hom-alt, so IBS0
   stays at zero and the KING numerator stays positive. 28 had no way to see this.

The fingerprint is in 28's own output. SRR8575311 is listed as the "same individual" as
eleven different people across four cohorts. One heart cannot be eleven donors. Those
are sticky samples, and 28 reported 8,272 first-degree and 40,652 second-degree pairs --
22% of the cohort called second-degree relatives, which is not biology.

THE TEST THAT SEPARATES THEM: MUTUAL BEST PARTNER
-------------------------------------------------
An individual sequenced twice has exactly ONE twin. So for a true duplicate pair (i,j),
j must be i's highest-kinship partner AND i must be j's. A sticky sample has many
partners but is the best partner of at most one of them, so every spurious pair it
generates fails the test from one side. This costs nothing, needs no threshold tuning,
and cannot be satisfied by a sample that is merely noisy.

Applied together:
  CONFIRMED   mutual best partner, kinship >= --dup-kin, IBS0 <= --dup-ibs0,
              and neither sample flagged by het QC
  CANDIDATE   mutual best partner, kinship >= --cand-kin (dropout can deflate a true
              duplicate, so this band is reported, not discarded)
  ARTEFACT    everything else that 28 called a duplicate

WHAT THIS IS EXPECTED TO SHOW ON THE REAL DATA
----------------------------------------------
23 MAGNet-LBNL pairs at kinship 0.42-0.49 with IBS0 exactly 0 and concordance
0.94-0.99, every one of them matching on disease group (MAGNet/NF with LBNL/NF,
MAGNet/DCM with LBNL/DCM). That is 23 of LBNL's 33 samples. It is consistent with LBNL
having sourced tissue from the same bank, and it is exactly what this whole exercise was
built to detect. It is not an artefact and it will not go away.

Usage:
  python 29_identity_triage.py
  python 29_identity_triage.py --dup-kin 0.30     # looser, after seeing the dropout shift
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

CODE = {"0/0": 0, "0|0": 0, "0/1": 1, "0|1": 1, "1|0": 1, "1/1": 2, "1|1": 2}


def per_sample_qc(mindp):
    """Call rate, het rate and hom-alt rate per sample, from geno/gt.txt."""
    if not (os.path.exists("geno/gt.samples") and os.path.exists("geno/gt.txt")):
        return None
    samples = open("geno/gt.samples").readline().rstrip("\n").split("\t")
    n = len(samples)
    ncall = np.zeros(n, dtype=np.int64)
    nhet = np.zeros(n, dtype=np.int64)
    nalt = np.zeros(n, dtype=np.int64)
    for line in open("geno/gt.txt"):
        f = line.rstrip("\n").split("\t")
        fmt = f[0].split(":")
        gi = fmt.index("GT") if "GT" in fmt else 0
        di = fmt.index("DP") if "DP" in fmt else None
        for k, x in enumerate(f[1:]):
            p = x.split(":")
            g = CODE.get(p[gi] if gi < len(p) else "./.", -1)
            if g < 0:
                continue
            if di is not None:
                dp = p[di] if di < len(p) else "0"
                if not (dp.isdigit() and int(dp) >= mindp):
                    continue
            ncall[k] += 1
            if g == 1:
                nhet[k] += 1
            elif g == 2:
                nalt[k] += 1
    Q = pd.DataFrame({"sample": samples, "n_called": ncall,
                      "het_rate": nhet / np.maximum(ncall, 1),
                      "homalt_rate": nalt / np.maximum(ncall, 1)})
    # robust outlier test: MAD, because the contaminated samples are exactly the ones
    # that would inflate a standard deviation and hide themselves
    m = Q.het_rate.median()
    mad = np.median(np.abs(Q.het_rate - m)) * 1.4826
    Q["het_z"] = (Q.het_rate - m) / mad if mad > 0 else 0.0
    return Q


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pairs", default="geno/identity_metrics.tsv")
    ap.add_argument("--mindp", type=int, default=6)
    ap.add_argument("--dup-kin", type=float, default=0.354,
                    help="kinship for a CONFIRMED duplicate (KING nominal)")
    ap.add_argument("--cand-kin", type=float, default=0.25,
                    help="kinship floor for a reported CANDIDATE")
    ap.add_argument("--dup-ibs0", type=float, default=0.01)
    ap.add_argument("--min-conc", type=float, default=0.85,
                    help="raw concordance floor. Unlike kinship this is ancestry-"
                         "sensitive, so it is used ONLY as a necessary condition: a "
                         "genuine duplicate always scores high (0.94-0.99 in this "
                         "dataset), so anything below the floor cannot be one whatever "
                         "its kinship says.")
    ap.add_argument("--het-z", type=float, default=4.0,
                    help="|MAD-z| of het rate above which a sample is QC-flagged")
    ap.add_argument("--out", default="geno/duplicates_triaged.tsv")
    a = ap.parse_args()

    if not os.path.exists(a.pairs):
        sys.exit(f"[fatal] {a.pairs} not found -- run 28_identity_metrics.py first")
    R = pd.read_csv(a.pairs, sep="\t")
    print(f"[info] {len(R)} pairs from {a.pairs}")

    # ---- 1. per-sample QC ----------------------------------------------------------
    Q = per_sample_qc(a.mindp)
    flagged = set()
    print()
    print("=" * 86)
    print("1. PER-SAMPLE GENOTYPE QC")
    print("=" * 86)
    if Q is None:
        print("  geno/gt.txt not found -- skipping het QC (the mutual-best test below")
        print("  still applies, and is the load-bearing one)")
    else:
        Q.to_csv("geno/sample_qc.tsv", sep="\t", index=False, float_format="%.5g")
        flagged = set(Q.loc[Q.het_z.abs() >= a.het_z, "sample"])
        print(f"  het rate: median {Q.het_rate.median():.4f}  "
              f"IQR {Q.het_rate.quantile(.25):.4f}-{Q.het_rate.quantile(.75):.4f}")
        print(f"  call rate: median {Q.n_called.median():.0f} sites")
        print(f"  samples with |MAD-z| >= {a.het_z} on het rate: {len(flagged)}")
        if flagged:
            print(f"    {'sample':<16}{'het rate':>10}{'MAD-z':>8}{'n called':>10}")
            for _, r in (Q[Q.het_z.abs() >= a.het_z]
                         .sort_values("het_z", ascending=False).head(25).iterrows()):
                print(f"    {r['sample']:<16}{r.het_rate:>10.4f}{r.het_z:>8.1f}"
                      f"{int(r.n_called):>10}")
            print("  An inflated het rate makes a sample look related to everyone.")
            print("  -> geno/sample_qc.tsv")

    # ---- 2. sticky samples ---------------------------------------------------------
    print()
    print("=" * 86)
    print("2. STICKY SAMPLES -- how many high-kinship partners does each sample have?")
    print("=" * 86)
    hi = R[R.kinship >= a.cand_kin]
    deg = pd.concat([hi.s1, hi.s2]).value_counts()
    print(f"  pairs at kinship >= {a.cand_kin}: {len(hi)}")
    print(f"  samples involved: {len(deg)}")
    print(f"    exactly 1 partner : {int((deg == 1).sum())}   <- consistent with a twin")
    print(f"    2 or more         : {int((deg >= 2).sum())}   <- cannot all be the same person")
    if (deg >= 3).any():
        print(f"\n  worst offenders (a heart cannot be N donors):")
        print(f"    {'sample':<16}{'partners':>10}{'het MAD-z':>12}")
        for s, c in deg[deg >= 3].head(15).items():
            z = (Q.loc[Q["sample"] == s, "het_z"].iloc[0]
                 if Q is not None and (Q["sample"] == s).any() else np.nan)
            zs = f"{z:.1f}" if np.isfinite(z) else "--"
            print(f"    {s:<16}{int(c):>10}{zs:>12}")

    # ---- 3. mutual best partner ----------------------------------------------------
    best = {}
    for s, k, o in np.concatenate([
            R[["s1", "kinship", "s2"]].to_numpy(),
            R[["s2", "kinship", "s1"]].to_numpy()]):
        if s not in best or k > best[s][0]:
            best[s] = (k, o)
    R["mutual_best"] = [best.get(x, (0, None))[1] == y and best.get(y, (0, None))[1] == x
                        for x, y in zip(R.s1, R.s2)]
    R["qc_flagged"] = [(x in flagged) or (y in flagged) for x, y in zip(R.s1, R.s2)]

    def verdict(r):
        if not r.mutual_best:
            return "artefact"
        # concordance floor: two contaminated samples can pair with each other at high
        # kinship and zero IBS0 (both are mostly het, so neither is ever hom-alt), but
        # their genotypes still disagree everywhere. No duplicate scores this low.
        if r.concordance < a.min_conc:
            return "artefact"
        if r.kinship >= a.dup_kin and r.ibs0 <= a.dup_ibs0 and not r.qc_flagged:
            return "confirmed"
        if r.kinship >= a.cand_kin:
            return "candidate"
        return "artefact"

    D = R[R.kinship >= a.cand_kin].copy()
    D["triage"] = D.apply(verdict, axis=1)
    D = D.sort_values("kinship", ascending=False)
    D.to_csv(a.out, sep="\t", index=False, float_format="%.5g")

    conf = D[D.triage == "confirmed"]
    cand = D[D.triage == "candidate"]
    art = D[D.triage == "artefact"]

    print()
    print("=" * 86)
    print("3. TRIAGE")
    print("=" * 86)
    print(f"  pairs 28 would have called duplicates or near-duplicates : {len(D)}")
    print(f"    CONFIRMED  mutual best, kinship >= {a.dup_kin}, IBS0 <= {a.dup_ibs0} : "
          f"{len(conf)}")
    print(f"    CANDIDATE  mutual best, kinship {a.cand_kin}-{a.dup_kin}            : "
          f"{len(cand)}")
    print(f"    ARTEFACT   fails mutual-best                              : {len(art)}")
    print(f"  -> {a.out}")

    for name, T in [("CONFIRMED", conf), ("CANDIDATE", cand)]:
        if not len(T):
            continue
        print()
        print(f"  {name}")
        print(f"    {'sample 1':<16}{'sample 2':<16}{'cohorts':<22}{'groups':<14}"
              f"{'kinship':>9}{'IBS0':>8}{'conc':>7}{'sites':>7}  note")
        for _, r in T.iterrows():
            same = "==" if str(r.group1) == str(r.group2) else "!="
            note = "het-QC flagged" if r.qc_flagged else ""
            print(f"    {r.s1[:15]:<16}{r.s2[:15]:<16}"
                  f"{r.cohort1 + '/' + r.cohort2:<22}"
                  f"{str(r.group1) + same + str(r.group2):<14}"
                  f"{r.kinship:>9.4f}{r.ibs0:>8.4f}{r.concordance:>7.3f}"
                  f"{int(r.n_shared):>7}  {note}")

    # ---- 4. what it means per cohort ------------------------------------------------
    print()
    print("=" * 86)
    print("4. IMPACT BY COHORT PAIR (confirmed only)")
    print("=" * 86)
    cx = conf[conf.cross_project] if "cross_project" in conf.columns else \
        conf[conf.cohort1 != conf.cohort2]
    if not len(cx):
        print("  no confirmed cross-project duplicates -- independence holds")
    else:
        key = cx.apply(lambda r: "/".join(sorted([r.cohort1, r.cohort2])), axis=1)
        for k, c in key.value_counts().items():
            sub = cx[key == k]
            agree = int((sub.group1.astype(str) == sub.group2.astype(str)).sum())
            print(f"  {k:<24}{c:>4} pair(s), {agree} with matching disease group")
        print()
        print("  Matching disease group on both sides is the tell. Random relatedness")
        print("  does not respect case/control labels; the same heart deposited twice")
        print("  does. High agreement means shared tissue sourcing, not chance.")
        # how much of each cohort is affected
        print()
        print(f"  {'cohort':<12}{'samples in confirmed cross-project pairs':>44}")
        for coh in sorted(set(cx.cohort1) | set(cx.cohort2)):
            n = len(set(cx.loc[cx.cohort1 == coh, "s1"]) |
                    set(cx.loc[cx.cohort2 == coh, "s2"]))
            print(f"  {coh:<12}{n:>44}")

    cw = conf[conf.cohort1 == conf.cohort2]
    if len(cw):
        print()
        print(f"  WITHIN-cohort confirmed pairs: {len(cw)} -- technical replicates or")
        print("  re-sequenced samples. Not an independence problem, but they inflate")
        print("  effective n; drop one of each.")
        for _, r in cw.iterrows():
            print(f"    {r.s1} == {r.s2}  ({r.cohort1})  kinship {r.kinship:.3f}")

    print()
    print("=" * 86)
    print("5. WHAT TO DO")
    print("=" * 86)
    print("  * Treat CONFIRMED as real. Drop each duplicate from ONE cohort -- keep it")
    print("    in the discovery cohort, or the larger one -- rerun the affected")
    print("    contrast, and state the exclusion in the methods.")
    print("  * Treat CANDIDATE as unresolved. Het dropout deflates kinship, so a true")
    print("    duplicate can land here. Check each one's IBS0 and shared-site count by")
    print("    hand, and re-run 28 with --mindp 8 (stricter genotypes, less dropout) to")
    print("    see whether it moves up into the confirmed band.")
    print("  * Treat ARTEFACT as noise. Do not report it, and do NOT quote 28's")
    print("    first-degree or second-degree counts anywhere -- they are dominated by")
    print("    the same artefact and are not a relatedness estimate.")
    print("=" * 86)


if __name__ == "__main__":
    main()
