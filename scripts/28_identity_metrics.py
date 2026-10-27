#!/usr/bin/env python3
"""
28_identity_metrics.py -- ancestry-robust sample-identity check on the genotypes
script 27 already produced. Reads geno/merged.vcf.gz; genotypes nothing itself.

WHY THIS EXISTS
---------------
The 33-sample trial of script 27 returned:

    unrelated bulk   0.564 +/- 0.069
    cross-project max 0.829
    threshold        0.900   (bulk mean + 8 sd, capped)

and a verdict of "no duplicates". That verdict was saved by the cap. 0.829 against a
0.85 line is not a margin, and it gets worse, not better, with the full cohort: the
unrelated maximum is the maximum of a growing number of draws. 33 samples gave 523
pairs; 670 samples give 224,115. If the bulk keeps sd ~0.07, the expected maximum
moves from about mean + 3.8 sd to about mean + 4.5 sd -- roughly 0.88 -- and a real
duplicate would no longer be distinguishable from the top of the unrelated tail.

The instinct is to add sites. It does not work here. Decompose the observed sd:

    sampling noise at 332 shared sites, p ~= 0.56 : sqrt(.56*.44/332) = 0.027
    observed                                      : 0.069
    remainder (shared ancestry between pairs)     : sqrt(.069^2 - .027^2) = 0.064

Tripling the panel drops the noise term to 0.016 and the total only to 0.066. The
spread is almost entirely the fact that two people of similar ancestry share more
genotypes than two people of different ancestry, and no amount of sites removes it.
Raw concordance is simply the wrong statistic for a multi-ancestry sample set.

WHAT REPLACES IT
----------------
Two statistics built to cancel population structure, both from the same genotypes.

  IBS0 -- the fraction of shared sites where one sample is hom-ref and the other is
  hom-alt. The same individual CANNOT produce an opposite-homozygote pair except
  through genotyping error, so IBS0 is ~0 for a duplicate and substantially positive
  for anyone else, at every allele frequency.

  KING-robust kinship (Manichaikul et al. 2010, Bioinformatics 26:2867, eq. 11):

        phi = (N_het,het - 2 * N_homref,homalt) / (N_het_i + N_het_j)

  where the het counts are taken over the sites shared by the pair. For a duplicate,
  N_homref,homalt = 0 and the three het counts are equal, so phi = 0.5. For unrelated
  individuals the numerator has expectation zero at ANY allele frequency, which is the
  whole point of the estimator -- it does not need the allele frequencies, so it cannot
  be biased by getting them wrong for a structured sample.

  KING's standard cut-points:
        > 0.354   duplicate / monozygotic twin
        0.177-0.354  first degree   (parent-offspring, full sibs)
        0.0884-0.177 second degree
        0.0442-0.0884 third degree
        < 0.0442  unrelated

  For the paper's claim, only the first line is disqualifying: the same heart in two
  bioprojects. Related individuals across cohorts would be a curiosity, not a breach of
  independence, but they are reported because they are worth knowing about.

RNA-SEQ CAVEAT, STATED HONESTLY
-------------------------------
These are genotypes called from RNA. Allele-specific expression and low coverage make
a true heterozygote look homozygous more often than in DNA, which deflates KING kinship
for a true duplicate (fewer het,het pairs) and inflates IBS0 slightly. Both errors move
a duplicate TOWARD the unrelated range, so this test is conservative for identity:
it can in principle miss a duplicate, it will not invent one. The script therefore
reports the highest-scoring pairs whatever the threshold says, so a near-miss is
visible rather than silently swallowed.

Usage
-----
  python 28_identity_metrics.py                       # after 27 has run
  python 28_identity_metrics.py --mindp 4             # recover more calls, noisier
  python 28_identity_metrics.py --min-shared 100      # thin panels

Costs seconds. Re-run it as often as you like; it never touches the BAMs.
"""
import argparse
import gzip
import os
import subprocess
import sys

import numpy as np
import pandas as pd

CODE = {"0/0": 0, "0|0": 0, "0/1": 1, "0|1": 1, "1|0": 1, "1/1": 2, "1|1": 2}

DUP = 0.354          # KING duplicate / MZ
DEG1 = 0.177         # KING first degree
DEG2 = 0.0884        # KING second degree


def read_genotypes(vcf, mindp):
    """Return (samples, G int8 sites x samples with -1 = no call)."""
    if not os.path.exists(vcf):
        sys.exit(f"[fatal] {vcf} not found. Run 27_genotype_concordance.sbatch first.")

    # Prefer the flat dumps 27 leaves behind; fall back to bcftools on the VCF.
    if os.path.exists("geno/gt.samples") and os.path.getsize("geno/gt.txt" if
                                                             os.path.exists("geno/gt.txt")
                                                             else vcf):
        pass
    if not (os.path.exists("geno/gt.samples") and os.path.exists("geno/gt.txt")):
        print("[info] geno/gt.txt absent -- re-extracting with bcftools")
        with open("geno/gt.samples", "w") as fh:
            hdr = subprocess.run(["bcftools", "view", "-h", vcf], capture_output=True,
                                 text=True, check=True).stdout
            line = [l for l in hdr.splitlines() if l.startswith("#CHROM")][-1]
            fh.write("\t".join(line.split("\t")[9:]) + "\n")
        with open("geno/gt.txt", "w") as fh:
            p1 = subprocess.Popen(["bcftools", "view", "-H", vcf],
                                  stdout=subprocess.PIPE)
            subprocess.run(["cut", "-f9,10-"], stdin=p1.stdout, stdout=fh, check=True)
            p1.stdout.close()
            p1.wait()

    samples = open("geno/gt.samples").readline().rstrip("\n").split("\t")
    gts = []
    nodp = 0
    opener = gzip.open if "geno/gt.txt".endswith(".gz") else open
    for line in opener("geno/gt.txt", "rt"):
        f = line.rstrip("\n").split("\t")
        fmt = f[0].split(":")
        gi = fmt.index("GT") if "GT" in fmt else 0
        di = fmt.index("DP") if "DP" in fmt else None
        if di is None:
            nodp += 1
        row = np.empty(len(f) - 1, dtype=np.int8)
        for k, x in enumerate(f[1:]):
            p = x.split(":")
            g = CODE.get(p[gi] if gi < len(p) else "./.", -1)
            if g >= 0 and di is not None:
                dp = p[di] if di < len(p) else "0"
                if not (dp.isdigit() and int(dp) >= mindp):
                    g = -1
            row[k] = g
        gts.append(row)
    if not gts:
        sys.exit("[fatal] geno/gt.txt is empty")
    if nodp:
        print(f"[warn] {nodp} sites carried no FORMAT/DP -- depth filter not applied "
              f"there")
    G = np.vstack(gts)
    if G.shape[1] != len(samples):
        sys.exit(f"[fatal] {G.shape[1]} genotype columns but {len(samples)} sample names")
    return samples, G


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", default="geno/merged.vcf.gz")
    ap.add_argument("--samples-tsv", default="geno/samples.tsv")
    ap.add_argument("--mindp", type=int, default=6)
    ap.add_argument("--min-shared", type=int, default=200)
    ap.add_argument("--min-callrate", type=float, default=0.30,
                    help="site kept if called in this fraction of samples")
    ap.add_argument("--out", default="geno/identity_metrics.tsv")
    ap.add_argument("--top", type=int, default=20)
    a = ap.parse_args()

    samples, G = read_genotypes(a.vcf, a.mindp)
    n = len(samples)
    print(f"[info] {G.shape[0]} sites x {n} samples, DP>={a.mindp}")
    print(f"[info] call rate: {(G >= 0).mean():.3f}")

    meta_c = {s: "?" for s in samples}
    meta_g = {s: "?" for s in samples}
    if os.path.exists(a.samples_tsv):
        M = pd.read_csv(a.samples_tsv, sep="\t").set_index("sample")
        for s in samples:
            if s in M.index:
                meta_c[s] = M.loc[s, "cohort"]
                meta_g[s] = M.loc[s, "group"]
    else:
        print(f"[warn] {a.samples_tsv} missing -- cohorts unknown, cross-project "
              f"flagging disabled")

    called = G >= 0
    ncall = called.sum(axis=1)
    nonref = ((G == 1) | (G == 2)).sum(axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        frac = np.where(ncall > 0, nonref / np.maximum(ncall, 1), 0.0)
    keep = (ncall >= a.min_callrate * n) & (frac >= 0.10) & (frac <= 0.90)
    G, called = G[keep], called[keep]
    ns = G.shape[0]
    print(f"[info] informative sites: {ns}")
    if ns < 100:
        sys.exit("[fatal] <100 informative sites -- widen BAND_LO/BAND_HI in 27, "
                 "or lower --mindp.")

    # ---- all-pairs counts as matrix products -------------------------------------
    # every term below is restricted to sites CALLED IN BOTH members of the pair,
    # which is what makes KING's cancellation hold pairwise.
    C = called.astype(np.float32)
    HR = (G == 0).astype(np.float32)          # hom-ref
    HE = (G == 1).astype(np.float32)          # het
    HA = (G == 2).astype(np.float32)          # hom-alt

    n_shared = C.T @ C                        # sites called in both
    n_hethet = HE.T @ HE
    n_opphom = HR.T @ HA + HA.T @ HR          # IBS0 count
    n_het_i = HE.T @ C                        # het in i, called in j
    n_het_j = n_het_i.T
    n_match = HR.T @ HR + n_hethet + HA.T @ HA

    iu = np.triu_indices(n, k=1)
    shared = n_shared[iu]
    ok = shared >= a.min_shared
    print(f"[info] pairs with >={a.min_shared} shared sites: {int(ok.sum())} of "
          f"{len(shared)}  ({len(shared) - int(ok.sum())} skipped)")
    if not ok.any():
        sys.exit("[fatal] no pair shares enough sites -- lower --min-shared or "
                 "widen the panel")

    den = (n_het_i + n_het_j)[iu]
    with np.errstate(invalid="ignore", divide="ignore"):
        kin = np.where(den > 0, (n_hethet[iu] - 2 * n_opphom[iu]) / np.maximum(den, 1),
                       np.nan)
        ibs0 = n_opphom[iu] / np.maximum(shared, 1)
        conc = n_match[iu] / np.maximum(shared, 1)

    s1 = np.array(samples)[iu[0]]
    s2 = np.array(samples)[iu[1]]
    R = pd.DataFrame({
        "s1": s1, "s2": s2,
        "cohort1": [meta_c[x] for x in s1], "cohort2": [meta_c[x] for x in s2],
        "group1": [meta_g[x] for x in s1], "group2": [meta_g[x] for x in s2],
        "n_shared": shared, "n_het_sum": den,
        "kinship": kin, "ibs0": ibs0, "concordance": conc,
    })[ok].copy()
    R["cross_project"] = R.cohort1 != R.cohort2
    R = R.sort_values("kinship", ascending=False)
    R.to_csv(a.out, sep="\t", index=False, float_format="%.5g")
    print(f"[write] {a.out}  ({len(R)} pairs, {int(R.cross_project.sum())} cross-project)")
    print(f"[info] shared sites per pair: median {int(R.n_shared.median())}  "
          f"min {int(R.n_shared.min())}  max {int(R.n_shared.max())}")

    # Provisional classification on kinship alone, used only to define the unrelated
    # bulk. The final call below also uses IBS0, because RNA-seq het dropout deflates
    # kinship for everyone and can push a true duplicate under the nominal cut-point.
    def degree(k):
        if not np.isfinite(k):
            return "undetermined"
        if k > DUP:
            return "duplicate"
        if k > DEG1:
            return "1st_degree"
        if k > DEG2:
            return "2nd_degree"
        if k > 0.0442:
            return "3rd_degree"
        return "unrelated"

    R["relationship"] = [degree(k) for k in R.kinship]

    # --- IBS0 as an independent detector -------------------------------------------
    # Opposite homozygotes are impossible for one individual sequenced twice, so IBS0
    # is bounded by the genotype error rate for a duplicate and sits near 2*p^2*q^2
    # summed over sites for anyone else. Unlike kinship it is not deflated by het
    # dropout -- dropout raises it, which again errs toward missing a duplicate rather
    # than inventing one. A pair whose IBS0 is an order of magnitude below the unrelated
    # floor is a duplicate candidate whatever its kinship says.
    bulk_mask = R.kinship < DEG2
    ibs0_floor = R.loc[bulk_mask, "ibs0"].min() if bulk_mask.any() else np.nan
    ibs0_med = R.loc[bulk_mask, "ibs0"].median() if bulk_mask.any() else np.nan
    ibs0_cut = ibs0_floor / 4 if np.isfinite(ibs0_floor) else -1
    R["ibs0_flag"] = R.ibs0 < ibs0_cut
    # promote: very low IBS0 plus kinship well above the unrelated range
    promote = R.ibs0_flag & (R.kinship > DEG1) & (R.relationship != "duplicate")
    if promote.any():
        R.loc[promote, "relationship"] = "duplicate"

    # ---- the distribution ---------------------------------------------------------
    print()
    print("=" * 88)
    print("KING-ROBUST KINSHIP  (duplicate 0.5, 1st degree 0.25, unrelated 0)")
    print("=" * 88)
    kk = R.kinship.dropna()
    h, e = np.histogram(kk, bins=np.arange(-0.30, 0.55, 0.025))
    for c, lo in zip(h, e[:-1]):
        bar = "#" * min(60, int(60 * c / max(h.max(), 1)))
        print(f"  {lo:+.3f} to {lo + 0.025:+.3f} {c:>9}  {bar}")
    unrel = kk[kk < DEG2]
    print(f"\n  unrelated bulk : mean {unrel.mean():+.4f}  sd {unrel.std():.4f}  "
          f"max {unrel.max():+.4f}")
    print(f"  overall max    : {kk.max():+.4f}")
    gap = (kk.max() - unrel.mean()) / max(unrel.std(), 1e-9)
    print(f"  top pair sits {gap:.1f} sd above the unrelated mean")

    # Head-to-head against the statistic this replaces, on the same pairs. If raw
    # concordance separates just as well on the real data, say so -- the argument for
    # KING is empirical, not a matter of taste, and it should be checked here rather
    # than asserted.
    cbulk = R.loc[R.kinship < DEG2, "concordance"]
    cgap = (R.concordance.max() - cbulk.mean()) / max(cbulk.std(), 1e-9)
    print()
    print("  SEPARATION, HEAD TO HEAD (top pair, in sd above the unrelated bulk)")
    print(f"    KING kinship      {gap:>8.1f} sd   "
          f"(bulk {unrel.mean():+.4f} +/- {unrel.std():.4f}, top {kk.max():+.4f})")
    print(f"    raw concordance   {cgap:>8.1f} sd   "
          f"(bulk {cbulk.mean():.4f} +/- {cbulk.std():.4f}, top {R.concordance.max():.4f})")
    if cgap >= gap:
        print("    -> concordance separates at least as well here. Either statistic")
        print("       supports the verdict; quote kinship anyway, because its scale is")
        print("       interpretable (0.5 = same individual) and its threshold is not")
        print("       derived from this dataset.")
    else:
        print(f"    -> kinship separates {gap / max(cgap, 1e-9):.1f}x better. This is why")
        print("       the verdict should be quoted from kinship, not concordance.")
    print(f"\n  IBS0 (opposite homozygotes): unrelated median {ibs0_med:.4f}, "
          f"floor {ibs0_floor:.4f}; duplicate-candidate cut {ibs0_cut:.4f}")
    print("  A true duplicate should show IBS0 near 0 AND kinship near 0.5. One without")
    print("  the other is a red flag that the genotypes, not the samples, are the issue.")

    # --- het dropout: the systematic shift that matters ------------------------------
    # RNA-seq calls a true heterozygote as a homozygote whenever one allele is not
    # captured (allele-specific expression, or simply low depth). That removes het,het
    # pairs and creates opposite-homozygote pairs, so it drags EVERY kinship downward.
    # It is not population structure and it is not fixed by more sites.
    shift = unrel.mean()
    if shift < -0.02:
        print(f"\n  [dropout] the unrelated bulk is centred at {shift:+.4f}, not 0.")
        print("  That is heterozygote dropout, not population structure -- it subtracts")
        print("  roughly the same amount from every pair, duplicates included, so the")
        print(f"  nominal 0.354 cut-point is CONSERVATIVE here. A duplicate would land")
        print(f"  near {0.5 + shift:.2f} rather than 0.50. Any pair between 0.20 and 0.354")
        print("  with near-zero IBS0 must be treated as a duplicate candidate, not")
        print("  dismissed. Lowering --mindp makes this worse; raising it makes it better")
        print("  at the cost of shared sites.")

    print()
    print(f"  TOP {a.top} PAIRS BY KINSHIP")
    print(f"    {'sample 1':<16}{'sample 2':<16}{'cohorts':<22}{'kinship':>9}"
          f"{'IBS0':>8}{'conc':>7}{'sites':>7}  call")
    for _, r in R.head(a.top).iterrows():
        cc = f"{r.cohort1}/{r.cohort2}"
        mark = "  <<<" if r.cross_project and r.kinship > DEG1 else ""
        print(f"    {r.s1[:15]:<16}{r.s2[:15]:<16}{cc:<22}{r.kinship:>9.4f}"
              f"{r.ibs0:>8.4f}{r.concordance:>7.3f}{int(r.n_shared):>7}  "
              f"{r.relationship}{mark}")

    print()
    print("  RELATIONSHIP COUNTS")
    for k in ["duplicate", "1st_degree", "2nd_degree", "3rd_degree", "unrelated",
              "undetermined"]:
        sub = R[R.relationship == k]
        if len(sub):
            print(f"    {k:<14}{len(sub):>8}   cross-project {int(sub.cross_project.sum()):>6}")

    # ---- verdict -------------------------------------------------------------------
    dup = R[(R.relationship == "duplicate")]
    dup_x = dup[dup.cross_project]
    rel_x = R[R.cross_project & (R.relationship == "1st_degree")]

    print()
    print("=" * 88)
    if len(dup_x) == 0:
        print("VERDICT: NO cross-bioproject duplicates. Independence verified.")
        print(f"         highest cross-project kinship {R[R.cross_project].kinship.max():+.4f}"
              f"  (duplicate threshold {DUP})")
        print()
        print("Methods sentence -- check every number against the run above before use:")
        print(f"  'Sample independence across bioprojects was verified from genotypes.")
        print(f"   Genotypes were called with bcftools mpileup/call from the recalibrated")
        print(f"   alignments at {ns} common autosomal biallelic SNVs (minimum depth")
        print(f"   {a.mindp}), and KING-robust kinship coefficients and IBS0 computed for")
        print(f"   all {len(R)} sample pairs sharing at least {a.min_shared} called sites.")
        print(f"   No cross-cohort pair reached the duplicate threshold of 0.354 (maximum")
        print(f"   {R[R.cross_project].kinship.max():.3f}; unrelated pairs")
        print(f"   {unrel.mean():.3f} +/- {unrel.std():.3f}), excluding the presence of any")
        print(f"   individual in more than one cohort. KING kinship was used in preference")
        print(f"   to raw genotype concordance because it is invariant to the allele-")
        print(f"   frequency differences between the ancestrally distinct cohorts.'")
    else:
        print(f"*** VERDICT: {len(dup_x)} CROSS-PROJECT DUPLICATE PAIR(S).")
        for _, r in dup_x.iterrows():
            print(f"***   {r.s1} ({r.cohort1}/{r.group1})  ==  {r.s2} "
                  f"({r.cohort2}/{r.group2})   kinship {r.kinship:.3f}, IBS0 {r.ibs0:.4f}")
        print("*** The same individual is in two bioprojects. Drop each duplicate from")
        print("*** ONE cohort -- keep it where it contributes more, usually the larger")
        print("*** group -- rerun phase 2 for the affected contrast, and state the")
        print("*** exclusion in the methods. Until then the independent-replication")
        print("*** claim is overstated.")
        dup_x.to_csv("geno/DUPLICATES_cross_project.tsv", sep="\t", index=False,
                     float_format="%.5g")
        print("*** -> geno/DUPLICATES_cross_project.tsv")

    dup_w = dup[~dup.cross_project]
    if len(dup_w):
        print(f"\nnote: {len(dup_w)} WITHIN-project duplicate pair(s) -- technical")
        print("replicates or re-sequenced samples. Not an independence problem between")
        print("cohorts, but they inflate that cohort's effective n. Drop one of each.")
        dup_w.to_csv("geno/duplicates_within_project.tsv", sep="\t", index=False,
                     float_format="%.5g")
        for _, r in dup_w.head(10).iterrows():
            print(f"  {r.s1} == {r.s2}  ({r.cohort1})  kinship {r.kinship:.3f}")

    if len(rel_x):
        print(f"\nnote: {len(rel_x)} cross-project FIRST-DEGREE pair(s) (kinship "
              f"{DEG1}-{DUP}). Not the same person, so independence holds, but worth a")
        print("sentence if any of them fall in the same contrast.")
        for _, r in rel_x.head(10).iterrows():
            print(f"  {r.s1} ({r.cohort1}) ~ {r.s2} ({r.cohort2})  "
                  f"kinship {r.kinship:.3f}, IBS0 {r.ibs0:.4f}")

    # ---- sanity checks on the test itself -------------------------------------------
    print()
    print("SANITY CHECKS -- read these before trusting the verdict")

    # 1. does structure actually cancel? The honest test is not whether the global mean
    #    is zero -- dropout moves that -- but whether the mean is the SAME for pairs
    #    drawn from cohorts of different ancestry as for pairs from the same cohort.
    #    If it is, the estimator is doing its job and a duplicate cannot hide in a
    #    cohort-specific offset.
    B = R[bulk_mask].copy()
    B["pairkey"] = ["/".join(sorted([a, b])) for a, b in zip(B.cohort1, B.cohort2)]
    tab = (B.groupby("pairkey")["kinship"]
             .agg(n="size", mean="mean", sd="std", mx="max")
             .sort_values("mean"))
    tab = tab[tab.n >= 30]
    if len(tab) >= 2:
        spread = tab["mean"].max() - tab["mean"].min()
        print(f"  structure cancellation: unrelated kinship by cohort pair "
              f"(spread {spread:+.4f})")
        for k, r in tab.iterrows():
            print(f"    {k:<26}{int(r.n):>8}   mean {r['mean']:+.4f}  "
                  f"sd {r['sd']:.4f}  max {r['mx']:+.4f}")
        if spread < 0.03:
            print("    -> OK, structure cancels")
        else:
            lowest = tab.index[0]
            print(f"    -> cohort-dependent offset of {spread:.3f}. KING is robust to")
            print(f"       structure but not perfectly so at large Fst; ancestrally")
            print(f"       distant cohort pairs (lowest: {lowest}) sit lower.")
            print("       Note the DIRECTION: unrelated cross-ancestry pairs are pushed")
            print("       AWAY from the duplicate threshold, not toward it, so this makes")
            print("       the duplicate verdict conservative rather than optimistic. A")
            print("       duplicate is the same person in both samples, so it carries no")
            print("       cross-ancestry offset at all. What the offset does compromise is")
            print("       the intermediate degrees: a cross-cohort sib pair is understated")
            print("       by roughly this much, so quote first/second-degree calls with")
            print("       that caveat or not at all.")
    print(f"  unrelated kinship spread            : sd {unrel.std():.4f}   "
          f"{'OK' if unrel.std() < 0.05 else 'WIDE -- too few sites per pair'}")

    # 2. near misses. A verdict of "none" is only worth something if nothing came close.
    hi = R[(R.kinship > DEG2) & (R.relationship != "duplicate")].sort_values(
        "kinship", ascending=False)
    if len(hi):
        print(f"  {len(hi)} pair(s) above 2nd-degree but below the duplicate cut "
              f"-- top {min(5, len(hi))}:")
        for _, r in hi.head(5).iterrows():
            print(f"    {r.s1} / {r.s2}  ({r.cohort1}/{r.cohort2})  "
                  f"kinship {r.kinship:+.4f}  IBS0 {r.ibs0:.4f}  "
                  f"{'<-- LOW IBS0, inspect' if r.ibs0_flag else ''}")
    else:
        print("  near misses: none -- nothing sits between 2nd degree and the "
              "duplicate cut")

    # 3. positive control
    w = R[~R.cross_project]
    print("  positive control (within-cohort)    : ", end="")
    if len(w):
        print(f"max kinship {w.kinship.max():+.4f} in "
              f"{w.loc[w.kinship.idxmax(), 'cohort1']}")
    else:
        print("none -- no within-cohort pairs compared")
    if not len(dup):
        print("  NOTE: no duplicate pair anywhere, so this run contains no positive")
        print("  control. The test's power then rests entirely on the two diagnostics")
        print("  above: a tight unrelated distribution with no cohort-specific offset,")
        print("  and nothing sitting near the cut-point.")
    print("=" * 88)


if __name__ == "__main__":
    main()
