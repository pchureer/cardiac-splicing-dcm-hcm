#!/usr/bin/env python3
"""
06_sensitivity_concordance.py

Turns the three 04 runs into the numbers a reviewer actually wants: not
"6037 vs 5619 significant events" (counts can agree while the underlying
events differ completely) but how much the *same events* survive, and how
well the effect sizes agree.

Also reports EXON-LEVEL DUPLICATION. rMATS emits one row per (cassette exon +
flanking exon pair), so the same cassette exon appears several times with
different flanking context. phase2_diffsplice_v2 dedups on the full coordinate
tuple, which keeps those as separate events -- but replication_concordance.py
and lvef_correlation.py key on (chr, exonStart, exonEnd) alone, i.e. they treat
them as one. That inconsistency inflates your FDR denominator and double-counts
genes in the RBM20 tables. This script quantifies it so you can decide.

Usage:
  python 06_sensitivity_concordance.py \
     --primary SE_v2filtered.all_events.tsv \
     --compare pool=SE_v2pooladj.all_events.tsv age25=SE_v2age25.all_events.tsv \
     --out SE_sensitivity
"""
import argparse
import sys

import numpy as np
import pandas as pd
from scipy import stats


def load(path):
    t = pd.read_csv(path, sep="\t")
    need = {"key", "gene", "chr", "exonStart", "exonEnd", "dPSI_adj", "FDR"}
    miss = need - set(t.columns)
    if miss:
        sys.exit(f"[fatal] {path} missing {sorted(miss)}")
    t["exon_key"] = (t["chr"].astype(str) + ":" + t["exonStart"].astype(str)
                     + ":" + t["exonEnd"].astype(str))
    ndup = len(t) - t["key"].nunique()
    if ndup:
        print(f"[warn] {path}: {ndup} duplicate full keys -- merges would fan out. "
              "Keeping the most significant row per key.", file=sys.stderr)
        t = t.sort_values("pvalue").drop_duplicates("key", keep="first")
    return t


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--primary", required=True)
    ap.add_argument("--compare", nargs="+", required=True,
                    help="label=path pairs")
    ap.add_argument("--fdr", type=float, default=0.05)
    ap.add_argument("--min-dpsi", type=float, default=0.10)
    ap.add_argument("--genes", default="TTN,CAMK2D,LDB3,LMO7,SORBS1,PDLIM3,DTNA,MLIP,TNNT2")
    ap.add_argument("--out", default="sensitivity")
    a = ap.parse_args()

    lines = []
    def say(s=""):
        lines.append(s)
        print(s)

    P = load(a.primary)
    genes = [g.strip() for g in a.genes.split(",") if g.strip()]

    say("=" * 70)
    say(f"PRIMARY: {a.primary}")
    say(f"  events (full coordinate key): {len(P)}")
    say(f"  unique cassette exons        : {P.exon_key.nunique()}")
    dup = len(P) / max(P.exon_key.nunique(), 1)
    say(f"  duplication factor           : {dup:.2f}x")
    say()

    sigP = P[(P.FDR < a.fdr) & (P.dPSI_adj.abs() >= a.min_dpsi)]
    say(f"--- 1. Significant at FDR<{a.fdr} and |dPSI|>={a.min_dpsi} ---")
    say(f"  by full coordinate key : {len(sigP)}")
    say(f"  by cassette exon       : {sigP.exon_key.nunique()}")
    say(f"  genes represented      : {sigP.gene.nunique()}")
    say("  The second number is what replication_concordance.py and")
    say("  lvef_correlation.py actually count, because they key on")
    say("  (chr, exonStart, exonEnd). Quote that one in the manuscript, or")
    say("  make every script use the same key. Do not mix them.")
    say()

    prog = sigP[sigP.gene.isin(genes)]
    say(f"--- 2. RBM20 program within the significant set ---")
    say(f"  events {len(prog)} across {prog.gene.nunique()} genes; "
        f"unique exons {prog.exon_key.nunique()}")
    if len(prog):
        say(f"  {'gene':<10}{'rows':>6}{'uniq exons':>12}{'median dPSI':>13}")
        for g, gg in prog.groupby("gene"):
            say(f"  {g:<10}{len(gg):>6}{gg.exon_key.nunique():>12}"
                f"{np.median(gg.dPSI_adj):>+13.3f}")
    say()

    say("--- 3. Sensitivity: do the SAME events survive? ---")
    setP_full = set(sigP.key)
    setP_exon = set(sigP.exon_key)
    for spec in a.compare:
        if "=" not in spec:
            sys.exit(f"[fatal] --compare needs label=path, got {spec}")
        lab, path = spec.split("=", 1)
        C = load(path)
        sigC = C[(C.FDR < a.fdr) & (C.dPSI_adj.abs() >= a.min_dpsi)]
        sf, se_ = set(sigC.key), set(sigC.exon_key)
        keep_f = len(setP_full & sf) / max(len(setP_full), 1)
        keep_e = len(setP_exon & se_) / max(len(setP_exon), 1)
        m = P[["key", "dPSI_adj"]].merge(C[["key", "dPSI_adj"]], on="key",
                                         suffixes=("_P", "_C"))
        rho, _ = stats.spearmanr(m.dPSI_adj_P, m.dPSI_adj_C)
        r, _ = stats.pearsonr(m.dPSI_adj_P, m.dPSI_adj_C)
        # restricted to primary-significant events, the number that matters
        ms = m[m.key.isin(setP_full)]
        rho_s, _ = stats.spearmanr(ms.dPSI_adj_P, ms.dPSI_adj_C)
        same = float(np.mean(np.sign(ms.dPSI_adj_P) == np.sign(ms.dPSI_adj_C)))
        say(f"  [{lab}]  {path}")
        say(f"    significant events: {len(sigC)} (primary {len(sigP)})")
        say(f"    primary-significant events retained: {keep_f:.1%} by full key, "
            f"{keep_e:.1%} by exon")
        say(f"    dPSI agreement, all common events      : rho={rho:+.3f} r={r:+.3f}"
            f"  (n={len(m)})")
        say(f"    dPSI agreement, primary-significant only: rho={rho_s:+.3f}"
            f"  same-sign {same:.1%}  (n={len(ms)})")
        prog_lost = sigP[sigP.gene.isin(genes) & ~sigP.key.isin(sf)]
        say(f"    RBM20-program events lost: {len(prog_lost)}"
            + (f" ({sorted(set(prog_lost.gene))})" if len(prog_lost) else ""))
        say()

    say("--- 4. How to report this ---")
    say("  'Adjusting additionally for sequencing pool, X% of significant events")
    say("   were retained with effect sizes correlated at rho=Y. Excluding donors")
    say("   under 25 retained Z% at rho=W.' Counts alone are not evidence of")
    say("   robustness -- two runs can report the same number of hits and share")
    say("   few of them.")
    say("=" * 70)

    with open(f"{a.out}.report.txt", "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"\nwrote {a.out}.report.txt")


if __name__ == "__main__":
    main()
