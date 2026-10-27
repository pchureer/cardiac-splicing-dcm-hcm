#!/usr/bin/env python3
"""
73_table2_recount.py -- Table 2's twenty lead genes under BOTH composition models.

WHY
---
The Results says: "18 of the 20 lead genes in Table 2 are among the 590 retained
(binomial P = 0.009 against that 63.7% background rate)."

`72_run_cmfrac.sbatch` showed the nine-marker cardiomyocyte score was over-adjusting:
retention rises from 590/926 (63.7%) to 716/926 (77.3%) when the score is restricted to
markers not upregulated in disease, and the adjusted model then returns 1,054 significant
events -- MORE than the unadjusted 926. → `RESULT_composition_score_is_a_mediator.md`

Promoting the fraction-only model to primary therefore moves the binomial background from
63.7% to 77.3%, and **the same 18 of 20 stops being significant**:

    18 of 20 against 63.7%  ->  P = 0.0090
    18 of 20 against 77.3%  ->  P = 0.1354     not significant
    19 of 20 against 77.3%  ->  P = 0.0401
    20 of 20 against 77.3%  ->  P = 0.0058

The count should rise too -- overall retention rose 13.6 points, and the lead genes were
selected for replication breadth, so they are the most robust events in the set. **But it
has to be counted.** If it comes back at 18, the honest move is to drop the binomial
sentence, not to keep quoting the weaker background in order to keep a p-value.

THE MATCHING PROBLEM, and why this is not a one-liner
-----------------------------------------------------
`key` is `chr:exonStart:exonEnd:upstreamES:upstreamEE:downstreamES:downstreamEE`, and the
exon-level convention collapses to the first three fields -- `chr:exonStart:exonEnd`.

`Table2.tsv`'s `locus` is only `chr:exonStart`. It does NOT uniquely identify an
exon-level event: `chr12:91179361` appears with exonEnd 91179686 AND 91179732. Matching on
locus alone would silently pick the wrong variant for any gene with several cassette
forms, and produce a plausible-looking wrong count.

So each Table 2 row is matched on (chr, exonStart) and then **disambiguated by its own
recorded dPSI and FDR**, which the table carries. The script REFUSES unless all twenty
match exactly one event each.

GATES -- five of them, because the output is a p-value in the manuscript
------------------------------------------------------------------------
1. `Table2.tsv` has exactly 20 rows
2. the base table returns exactly 926 significant exons
3. the nine-marker adjusted table retains exactly 590 of them
4. all 20 Table 2 rows match exactly one exon-level event
5. the retained count under the nine-marker model equals **18**, and equals the number of
   `comp_adj == "+"` rows in `Table2.tsv` -- i.e. the script reproduces the published
   figure by its own route AND agrees with the table's own flag

Only then does it report the fraction-only count. Anything else and it writes nothing.

USAGE
-----
    conda activate rnaseq_tools
    cd /data/project/chureeratp/RNA/HEART
    python 73_table2_recount.py
"""
import argparse
import os
import sys
from math import comb

import numpy as np
import pandas as pd

# --------------------------------------------------------------------------
# Read length. The chain is 61 -> Table2<SUF>.tsv -> 73, so the suffix must
# flow through both or this script validates the 125 bp lead genes against the
# 50 bp tables.
SUF = os.environ.get("TBL_SUFFIX", "")
# --------------------------------------------------------------------------
TABLE2 = "Table2" + SUF + ".tsv"
BASE = "SE_v2filtered" + SUF + ".all_events.tsv"
COMPADJ = "SE_v2compadj" + SUF + ".all_events.tsv"
CMFRAC = "SE_v2cmfrac" + SUF + ".all_events.tsv"

NKEEP_SE = 3
FDR_CUT, DPSI_CUT = 0.05, 0.10

# the record
N_LEAD = 20
REC_BASE_SIG = 926
REC_COMPADJ_KEPT = 590
REC_LEAD_KEPT = 18

# Measured at 125 bp (RESULT_s1_composition_rl125_complete.md):
#     REC_BASE_SIG=813 REC_COMPADJ_KEPT=498
# REC_LEAD_KEPT is the quantity being MEASURED, so 0 means "no expectation,
# report it". The Table2-flag-vs-recomputation check below still fires: that is
# internal consistency and needs no recorded value.
REC_BASE_SIG = int(os.environ.get("REC_BASE_SIG", str(REC_BASE_SIG)))
REC_COMPADJ_KEPT = int(os.environ.get("REC_COMPADJ_KEPT",
                                      str(REC_COMPADJ_KEPT)))
REC_LEAD_KEPT = int(os.environ.get("REC_LEAD_KEPT", str(REC_LEAD_KEPT)))
# tolerances for disambiguating a Table 2 row against the events table
DPSI_TOL = 0.0015
FDR_REL_TOL = 0.05


def read_table2(path):
    t = pd.read_csv(path, sep="\t", low_memory=False)
    if t.shape[1] == 1:                       # not tab-separated after all
        t = pd.read_csv(path, sep=r"\s+", engine="python")
    t.columns = [c.strip() for c in t.columns]
    return t


def representatives(path):
    """Exon-level representatives and the significant set. Project convention."""
    t = pd.read_csv(path, sep="\t", low_memory=False)
    for c in ("key", "FDR"):
        if c not in t.columns:
            sys.exit(f"[fatal] {path}: no '{c}' column")
    d = next((c for c in ("dPSI_adj", "dPSI") if c in t.columns), None)
    if d is None:
        sys.exit(f"[fatal] {path}: no dPSI column")
    t = t.copy()
    t["ekey"] = [":".join(str(k).split(":")[:NKEEP_SE]) for k in t["key"]]
    t["_f"] = pd.to_numeric(t["FDR"], errors="coerce")
    t["_d"] = pd.to_numeric(t[d], errors="coerce")
    t = t.dropna(subset=["_f"])
    rep = t.loc[t.groupby("ekey")["_f"].idxmin()].set_index("ekey")
    sig = set(rep[(rep._f < FDR_CUT) & (rep._d.abs() >= DPSI_CUT)].index)
    return rep, sig


def binom_ge(k, n, p):
    return sum(comb(n, i) * p**i * (1 - p)**(n - i) for i in range(k, n + 1))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="Table2_recount" + SUF + ".tsv")
    a = ap.parse_args()

    print("\n" + "=" * 76)
    print("TABLE 2's LEAD GENES UNDER BOTH COMPOSITION MODELS")
    print("=" * 76)

    for p in (TABLE2, BASE, COMPADJ, CMFRAC):
        if not os.path.exists(p):
            sys.exit(f"[fatal] {p} not found"
                     + ("  -- has 72_run_cmfrac.sbatch finished?" if p == CMFRAC else ""))

    t2 = read_table2(TABLE2)
    print(f"\n  {TABLE2}: {len(t2)} rows, columns: {', '.join(t2.columns)}")
    if len(t2) != N_LEAD:
        sys.exit(f"[fatal] expected {N_LEAD} lead genes, found {len(t2)}")
    print(f"  [ok] {N_LEAD} lead genes")

    base_rep, base_sig = representatives(BASE)
    comp_rep, comp_sig = representatives(COMPADJ)
    frac_rep, frac_sig = representatives(CMFRAC)

    print(f"\n  base          {len(base_sig):5d} significant of {len(base_rep):,} exons")
    print(f"  nine-marker   {len(comp_sig):5d} significant of {len(comp_rep):,} exons")
    print(f"  fraction-only {len(frac_sig):5d} significant of {len(frac_rep):,} exons")

    fail = []
    if len(base_sig) != REC_BASE_SIG:
        fail.append(f"base gives {len(base_sig)}, recorded {REC_BASE_SIG}")
    else:
        print(f"  [ok] base reproduces {REC_BASE_SIG}")

    kept_comp = len(base_sig & comp_sig)
    kept_frac = len(base_sig & frac_sig)
    if kept_comp != REC_COMPADJ_KEPT:
        fail.append(f"nine-marker retention {kept_comp}, recorded {REC_COMPADJ_KEPT}")
    else:
        print(f"  [ok] nine-marker retention reproduces {REC_COMPADJ_KEPT}")
    print(f"       fraction-only retention {kept_frac} "
          f"({100*kept_frac/len(base_sig):.1f}%)")

    # ------------------------------------------------ match the twenty, unambiguously
    print(f"\n  matching the {N_LEAD} lead events to exon-level keys")
    print(f"  (locus is chr:exonStart only, so each is disambiguated by its own "
          f"dPSI and FDR)")
    rows, unmatched, tied = [], [], []
    for _, r in t2.iterrows():
        locus = str(r["locus"]).strip()
        try:
            chrom, start = locus.split(":")
        except ValueError:
            unmatched.append(f"{r['gene']}: cannot parse locus '{locus}'")
            continue
        prefix = f"{chrom}:{start}:"
        cand = [k for k in base_rep.index if str(k).startswith(prefix)]
        if not cand:
            unmatched.append(f"{r['gene']} ({locus}): no exon with that chr:exonStart")
            continue
        want_d = float(r["dPSI"])
        want_f = float(r["FDR"])
        hits = [k for k in cand
                if abs(base_rep._d[k] - want_d) <= DPSI_TOL
                and abs(base_rep._f[k] - want_f) <= FDR_REL_TOL * max(want_f, 1e-300)]
        if len(hits) == 0:
            unmatched.append(
                f"{r['gene']} ({locus}): {len(cand)} candidate exon(s), none match "
                f"dPSI={want_d:+.3f} FDR={want_f:.3g}")
            continue
        if len(hits) > 1:
            # TIED: two or more exon-level keys with IDENTICAL dPSI and FDR. They are
            # the same measurement -- rMATS gives several annotated 3' boundaries of one
            # cassette exon the same inclusion-junction counts, so PSI is identical and
            # chr:exonStart:exonEnd splits one measurement into several keys.
            # The ambiguity is harmless IF every tied key has the same retention status,
            # because then it does not matter which one Table 2 meant. Refuse only if
            # they disagree.
            st_ = {(k in base_sig, k in comp_sig, k in frac_sig) for k in hits}
            if len(st_) > 1:
                unmatched.append(
                    f"{r['gene']} ({locus}): {len(hits)} tied keys with identical "
                    f"dPSI/FDR but DIFFERENT retention status -- genuinely ambiguous")
                continue
            tied.append(f"{r['gene']} ({locus}): {len(hits)} tied keys, identical "
                        f"dPSI/FDR and identical status -- resolved")
        k = hits[0]
        rows.append(dict(rank=int(r["rank"]), gene=r["gene"], ekey=k,
                         locus=locus, dPSI=want_d, FDR=want_f,
                         n_tied_keys=len(hits),
                         table2_comp_adj=str(r.get("comp_adj", "")).strip(),
                         in_base=k in base_sig,
                         kept_nine_marker=k in comp_sig,
                         kept_fraction_only=k in frac_sig))

    if unmatched:
        print("\n  REFUSING -- these lead events could not be matched to exactly one "
              "exon-level key:")
        for u in unmatched:
            print(f"    * {u}")
        print("\n  A wrong match would give a plausible-looking wrong count, which is")
        print("  the failure mode this project has been bitten by repeatedly. Nothing")
        print("  written. Widen DPSI_TOL/FDR_REL_TOL only after checking why.\n")
        return
    if tied:
        print(f"\n  {len(tied)} lead event(s) matched SEVERAL keys carrying the same "
              f"measurement:")
        for t_ in tied:
            print(f"    * {t_}")
        print("    Resolved: every tied key has identical retention status, so which one")
        print("    Table 2 meant does not affect the count.")
    print(f"  [ok] all {len(rows)} lead events resolved")

    # ---------------------------------------- residual duplication in the 926
    # The exon-level key is chr:exonStart:exonEnd. If one cassette exon has several
    # annotated 3' boundaries, rMATS can give them identical junction counts and hence
    # identical PSI -- so the "exon-level" set still contains repeated MEASUREMENTS.
    # FABP3 and ELOC showed this among the 20 lead genes alone, so it is worth counting.
    from collections import defaultdict
    grp = defaultdict(list)
    for k in base_sig:
        f = str(k).split(":")
        if len(f) >= 3:
            grp[(f[0], f[1])].append(k)
    excess, affected = 0, []
    for (c_, s_), ks in grp.items():
        if len(ks) < 2:
            continue
        byval = defaultdict(list)
        for k in ks:
            byval[(round(float(base_rep._d[k]), 6),
                   float(base_rep._f[k]))].append(k)
        for _, same in byval.items():
            if len(same) > 1:
                excess += len(same) - 1
                affected.append(f"{c_}:{s_} x{len(same)}")
    print(f"\n  RESIDUAL DUPLICATION in the {len(base_sig)} significant exons:")
    print(f"    {excess} are repeat MEASUREMENTS -- same chr:exonStart, different")
    print(f"    exonEnd, identical dPSI and FDR. Deduplicated count: "
          f"{len(base_sig) - excess}")
    if excess:
        print(f"    affected loci (first 8): {', '.join(affected[:8])}")
        print(f"    >>> The exon-level convention collapses transcript context but NOT")
        print(f"        alternative 3' boundaries of the same exon. {excess} of "
              f"{len(base_sig)} is "
              f"{100*excess/len(base_sig):.1f}% of the headline count.")

    df = pd.DataFrame(rows).sort_values("rank")

    # ------------------------------------------------ gate on the published 18
    n_flag = int((df.table2_comp_adj == "+").sum())
    n_nine = int(df.kept_nine_marker.sum())
    n_frac = int(df.kept_fraction_only.sum())
    print(f"\n  Table2.tsv's own comp_adj flag:        {n_flag} of {N_LEAD}")
    print(f"  recomputed against the nine-marker:   {n_nine} of {N_LEAD}")
    if REC_LEAD_KEPT and n_nine != REC_LEAD_KEPT:
        fail.append(f"recomputed lead retention {n_nine}, recorded {REC_LEAD_KEPT}")
    if n_flag != n_nine:
        fail.append(f"Table2.tsv's flag says {n_flag}, recomputation says {n_nine}")
    if not fail:
        print(f"  [ok] both reproduce the published {REC_LEAD_KEPT}")

    if fail:
        print("\n" + "=" * 76)
        print("NOT REPORTED -- a gate failed")
        print("=" * 76)
        for f in fail:
            print(f"  * {f}")
        print("\n  Nothing written.\n")
        return

    # ------------------------------------------------ the answer
    p_nine = kept_comp / len(base_sig)
    p_frac = kept_frac / len(base_sig)
    print("\n" + "=" * 76)
    print("THE RESULT")
    print("=" * 76)
    print(f"\n  {'model':<22}{'background':>12}{'lead kept':>12}{'binomial P':>14}")
    print("  " + "-" * 58)
    print(f"  {'nine-marker':<22}{p_nine*100:>11.1f}%{n_nine:>8} of {N_LEAD}"
          f"{binom_ge(n_nine, N_LEAD, p_nine):>14.4f}")
    print(f"  {'fraction-only':<22}{p_frac*100:>11.1f}%{n_frac:>8} of {N_LEAD}"
          f"{binom_ge(n_frac, N_LEAD, p_frac):>14.4f}")

    pv = binom_ge(n_frac, N_LEAD, p_frac)
    print()
    if pv < 0.05:
        print(f"  >>> THE CLAIM SURVIVES on the better model: {n_frac} of {N_LEAD} lead")
        print(f"      genes retained against a {p_frac*100:.1f}% background, P = {pv:.4f}.")
        print(f"      Quote this and drop the 63.7% version -- the fraction-only model is")
        print(f"      the one that adjusts for composition rather than for disease.")
    else:
        print(f"  >>> THE CLAIM DOES NOT SURVIVE: {n_frac} of {N_LEAD} against a")
        print(f"      {p_frac*100:.1f}% background gives P = {pv:.4f}.")
        print(f"      DROP THE BINOMIAL SENTENCE. Do not keep quoting the 63.7%")
        print(f"      background to preserve a p-value -- that background comes from a")
        print(f"      covariate now known to absorb part of the disease effect, and")
        print(f"      choosing it because it gives a better p-value is the exact practice")
        print(f"      this paper criticises in the regulator literature.")
        print(f"      The substantive point stands without the test: all {N_LEAD} lead")
        print(f"      genes are concordant in every cohort where measurable.")

    # events that changed status between the two models -- the interesting rows
    moved = df[df.kept_nine_marker != df.kept_fraction_only]
    if len(moved):
        print(f"\n  {len(moved)} lead gene(s) change status between the two models:")
        for _, r in moved.iterrows():
            a_, b_ = r.kept_nine_marker, r.kept_fraction_only
            print(f"    {r.gene:<10} rank {r['rank']:<3} nine-marker "
                  f"{'kept' if a_ else 'LOST':<5} -> fraction-only "
                  f"{'kept' if b_ else 'LOST'}")
        print("  These are the genes the nine-marker score was over-adjusting away.")

    df.to_csv(a.out, sep="\t", index=False)
    print(f"\n  wrote {a.out}  ({len(df)} rows)\n")


if __name__ == "__main__":
    main()
