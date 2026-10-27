#!/usr/bin/env python3
"""
30_build_dedup_lists.py -- turn the 26 confirmed duplicate pairs into concrete
exclusion lists, filtered .ids / *_bams.txt files, and the exact list of contrasts
that must be rerun. Also audits the samples the identity check could not see.

Run after 29_identity_triage.py. Writes nothing destructive: every output is a NEW
file with a _dedup suffix, so the originals stay untouched and the change is reviewable
with diff before anything is rerun.

THE KEEP RULE
-------------
Cross-cohort pair: keep the sample in the higher-priority cohort, drop the other.
Priority is discovery first, then cohort size:
    MAGNet > Heinig > Colorado > Fuwai > LBNL > Brigham > Stanford > Amsterdam > Naples
This keeps the discovery cohort whole -- re-running discovery would invalidate every
downstream result -- and removes the sample from the cohort that was going to be used
as independent replication, which is the claim actually at risk.

Within-cohort pair: keep whichever sample has more called genotypes (from
geno/sample_qc.tsv), i.e. the better-covered library. Falls back to the first
accession alphabetically if the QC table is absent.

WHAT THE CHECK COULD NOT SEE -- report this as a limitation
-----------------------------------------------------------
Section 1 of the triage flagged 10 MAGNet samples with n_called = 0. They produced no
usable genotypes, so they were never compared to anything: they cannot be confirmed
free of duplication, and they cannot be confirmed as duplicates either. The methods
must say the check covered 660 of 670 samples, not all 670. This script re-checks their
BAMs so the cause can be stated rather than guessed.

Usage:
  python 30_build_dedup_lists.py
  python 30_build_dedup_lists.py --apply-to-bams   # also filter *_bams.txt
"""
import argparse
import glob
import os
import re
import subprocess
import sys

import pandas as pd

PRIORITY = ["MAGNet", "Heinig", "Colorado", "Fuwai", "LBNL", "Brigham", "Stanford",
            "Amsterdam", "Naples"]
RANK = {c: i for i, c in enumerate(PRIORITY)}
ACC = re.compile(r"(SRR\d+|ERR\d+|SAMEA\d+|SAMN\d+)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--triaged", default="geno/duplicates_triaged.tsv")
    ap.add_argument("--qc", default="geno/sample_qc.tsv")
    ap.add_argument("--samples", default="geno/samples.tsv")
    ap.add_argument("--root", default=".")
    ap.add_argument("--apply-to-bams", action="store_true")
    a = ap.parse_args()

    if not os.path.exists(a.triaged):
        sys.exit(f"[fatal] {a.triaged} not found -- run 29_identity_triage.py first")
    T = pd.read_csv(a.triaged, sep="\t")
    C = T[T.triage == "confirmed"].copy()
    print(f"[info] {len(C)} confirmed duplicate pairs")
    if C.empty:
        sys.exit("[info] nothing to exclude")

    ncall = {}
    if os.path.exists(a.qc):
        Q = pd.read_csv(a.qc, sep="\t")
        ncall = dict(zip(Q["sample"], Q.n_called))

    # ---- 1. decide which member of each pair to drop -------------------------------
    drops = []
    for _, r in C.iterrows():
        if r.cohort1 != r.cohort2:
            lo = r.s1 if RANK.get(r.cohort1, 99) > RANK.get(r.cohort2, 99) else r.s2
            keep = r.s2 if lo == r.s1 else r.s1
            why = "lower-priority cohort"
        else:
            n1, n2 = ncall.get(r.s1, -1), ncall.get(r.s2, -1)
            if n1 == n2:
                lo, keep = sorted([r.s1, r.s2])[1], sorted([r.s1, r.s2])[0]
                why = "tie on call count; kept alphabetically first"
            else:
                lo = r.s1 if n1 < n2 else r.s2
                keep = r.s2 if lo == r.s1 else r.s1
                why = f"fewer called sites ({min(n1, n2)} vs {max(n1, n2)})"
        drops.append({"drop": lo, "keep": keep,
                      "drop_cohort": r.cohort1 if lo == r.s1 else r.cohort2,
                      "keep_cohort": r.cohort2 if lo == r.s1 else r.cohort1,
                      "drop_group": r.group1 if lo == r.s1 else r.group2,
                      "kinship": r.kinship, "reason": why})
    D = pd.DataFrame(drops).sort_values(["drop_cohort", "drop_group", "drop"])
    D.to_csv("geno/EXCLUDE_samples.tsv", sep="\t", index=False, float_format="%.4f")

    print()
    print("=" * 88)
    print("1. SAMPLES TO EXCLUDE")
    print("=" * 88)
    for coh, g in D.groupby("drop_cohort"):
        print(f"\n  {coh}  ({len(g)} samples)")
        for _, r in g.iterrows():
            print(f"    drop {r['drop']:<16} ({r.drop_group:<5}) "
                  f"-- same individual as {r['keep']} ({r.keep_cohort}), "
                  f"kinship {r.kinship:.3f}")
    print(f"\n  -> geno/EXCLUDE_samples.tsv   ({len(D)} samples total)")
    dropped = set(D["drop"])

    # ---- 2. cohort sizes before and after ------------------------------------------
    print()
    print("=" * 88)
    print("2. COHORT SIZES BEFORE AND AFTER")
    print("=" * 88)
    if os.path.exists(a.samples):
        S = pd.read_csv(a.samples, sep="\t")
        S["kept"] = ~S["sample"].isin(dropped)
        print(f"  {'cohort':<12}{'group':<8}{'before':>8}{'after':>8}{'lost':>7}")
        hurt = []
        for (c, g), sub in S.groupby(["cohort", "group"]):
            b, k = len(sub), int(sub.kept.sum())
            if b != k:
                print(f"  {c:<12}{str(g):<8}{b:>8}{k:>8}{b - k:>7}")
                hurt.append((c, g, b, k))
        print()
        for c, g, b, k in hurt:
            if k < 8:
                print(f"  *** {c}/{g} drops to n={k}. Below any useful power; this arm")
                print(f"      cannot carry a replication claim any more.")

    # ---- 3. filtered .ids files ----------------------------------------------------
    print()
    print("=" * 88)
    print("3. FILTERED .ids FILES  (new files; originals untouched)")
    print("=" * 88)
    changed = []
    for f in sorted(glob.glob(os.path.join(a.root, "*.ids"))):
        if f.endswith("_dedup.ids"):
            continue
        ids = [l.strip() for l in open(f) if l.strip()]
        keep = [i for i in ids if i not in dropped]
        if len(keep) == len(ids):
            continue
        out = f[:-4] + "_dedup.ids"
        with open(out, "w") as fh:
            fh.write("\n".join(keep) + "\n")
        changed.append((os.path.basename(f), len(ids), len(keep), out))
        print(f"  {os.path.basename(f):<26}{len(ids):>5} -> {len(keep):<5} "
              f"{os.path.basename(out)}")
    if not changed:
        print("  no .ids file contained an excluded sample")

    # ---- 4. filtered bam lists -----------------------------------------------------
    if a.apply_to_bams:
        print()
        print("=" * 88)
        print("4. FILTERED *_bams.txt  (new files; originals untouched)")
        print("=" * 88)
        for f in sorted(glob.glob(os.path.join(a.root, "*_bams.txt"))):
            if f.endswith("_dedup.txt"):
                continue
            raw = open(f).read()
            sep = "," if "," in raw else "\n"
            items = [x.strip() for x in raw.replace("\n", sep).split(sep) if x.strip()]
            # Two naming conventions in these lists: an accession somewhere in the path,
            # or a per-sample directory whose basename IS the id. phase2 resolves both,
            # so this must too, or a dropped sample survives in the bam list while
            # disappearing from the .ids file -- the worst possible outcome, because
            # the two would then disagree silently.
            def is_dropped(p):
                m = ACC.search(p)
                if m and m.group(1) in dropped:
                    return True
                return os.path.basename(os.path.dirname(p)) in dropped

            keep = [p for p in items if not is_dropped(p)]
            if len(keep) == len(items):
                continue
            out = f.replace(".txt", "_dedup.txt")
            with open(out, "w") as fh:
                fh.write(sep.join(keep) + ("\n" if sep == "\n" else ""))
            print(f"  {os.path.basename(f):<26}{len(items):>5} -> {len(keep):<5} "
                  f"{os.path.basename(out)}")

    # ---- 5. what must be rerun -----------------------------------------------------
    print()
    print("=" * 88)
    print("5. CONTRASTS THAT MUST BE RERUN")
    print("=" * 88)
    for coh in sorted(set(D.drop_cohort)):
        n = int((D.drop_cohort == coh).sum())
        print(f"  {coh}: {n} sample(s) removed -> rerun every phase 2 contrast that "
              f"uses {coh}")
    print()
    print("  The discovery cohort is NOT in that list by construction, so the discovery")
    print("  results, the event convention and every MAGNet-internal number stand.")

    # ---- 6. the samples the check could not see ------------------------------------
    print()
    print("=" * 88)
    print("6. SAMPLES THE IDENTITY CHECK COULD NOT SEE")
    print("=" * 88)
    if not ncall:
        print("  geno/sample_qc.tsv absent -- cannot audit")
        return
    Q = pd.read_csv(a.qc, sep="\t")
    blind = Q[Q.n_called == 0]["sample"].tolist()
    if not blind:
        print("  none -- every sample produced genotypes")
        return
    print(f"  {len(blind)} sample(s) produced ZERO called genotypes. They were never")
    print("  compared to anything, so they are neither cleared nor implicated.")
    print("  State in the methods that the check covered "
          f"{len(Q) - len(blind)} of {len(Q)} samples.\n")
    S = pd.read_csv(a.samples, sep="\t").set_index("sample") \
        if os.path.exists(a.samples) else None
    print(f"    {'sample':<16}{'cohort':<11}{'bam':>6}{'index':>7}{'size (MB)':>11}")
    for s in blind:
        bam = S.loc[s, "bam"] if S is not None and s in S.index else ""
        coh = S.loc[s, "cohort"] if S is not None and s in S.index else "?"
        ex = os.path.exists(bam) if bam else False
        idx = (os.path.exists(bam + ".bai") or os.path.exists(bam[:-4] + ".bai")) \
            if bam else False
        sz = f"{os.path.getsize(bam) / 1e6:.0f}" if ex else "--"
        print(f"    {s:<16}{coh:<11}{'yes' if ex else 'NO':>6}"
              f"{'yes' if idx else 'NO':>7}{sz:>11}")
    print()
    print("  no index  -> bcftools mpileup skips the file silently. Run samtools index")
    print("               on each, then re-run 27 with RESUME=1 after deleting the")
    print("               chunk VCF that contains them, and re-run 28 and 29.")
    print("  no BAM    -> the sample is not analysable at all; check it is not in any")
    print("               contrast, because phase 2 would have skipped it too.")
    print("  BAM+index present but still zero calls -> different contig naming than the")
    print("               reference (chr1 vs 1), or a coordinate-unsorted file.")
    print("=" * 88)


if __name__ == "__main__":
    main()
