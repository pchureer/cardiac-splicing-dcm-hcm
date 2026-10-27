#!/usr/bin/env python3
"""
10_why_events_were_lost.py

~40% of significant events dropped out on composition adjustment. There are TWO
very different reasons that can happen, and they mean opposite things:

  (a) CONFOUNDING -- the event's effect was composition all along. Its adjusted
      dPSI collapses toward zero. This would mean the signature is partly an
      artefact of tissue makeup.

  (b) THRESHOLD ATTRITION -- the effect is unchanged but the standard error grew
      (extra covariates cost degrees of freedom and absorb variance correlated
      with disease), so borderline events slipped below FDR<0.05 or |dPSI|>=0.10.
      This means nothing was wrong; the test simply got more conservative.

The concordance report cannot tell these apart. This can: it compares the LOST
events against the KEPT events on how much their effect size actually moved, and
how much their standard error grew.

Read it like this:
  lost events, median dPSI ratio near 0    -> (a) confounding, real problem
  lost events, median dPSI ratio near 1    -> (b) attrition, not a problem
  lost events with a much larger SE increase than kept events -> (b)

Usage:
  python 10_why_events_were_lost.py --base SE_v2filtered.all_events.tsv \
      --adj SE_v2compadj.all_events.tsv --label composition --out SE_lostwhy
"""
import argparse

import numpy as np
import pandas as pd
from scipy import stats


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", required=True)
    ap.add_argument("--adj", required=True)
    ap.add_argument("--label", default="adjusted")
    ap.add_argument("--fdr", type=float, default=0.05)
    ap.add_argument("--min-dpsi", type=float, default=0.10)
    ap.add_argument("--genes", default="TTN,CAMK2D,LDB3,LMO7,SORBS1,PDLIM3,DTNA,MLIP,TNNT2")
    ap.add_argument("--out", default="lostwhy")
    a = ap.parse_args()

    b = pd.read_csv(a.base, sep="\t").drop_duplicates("key")
    c = pd.read_csv(a.adj, sep="\t").drop_duplicates("key")
    m = b.merge(c, on="key", suffixes=("_b", "_a"))

    sig = m[(m.FDR_b < a.fdr) & (m.dPSI_adj_b.abs() >= a.min_dpsi)].copy()
    keptmask = (sig.FDR_a < a.fdr) & (sig.dPSI_adj_a.abs() >= a.min_dpsi)
    kept, lost = sig[keptmask], sig[~keptmask]

    lines = []
    def say(s=""):
        lines.append(s)
        print(s)

    say("=" * 68)
    say(f"WHY EVENTS WERE LOST -- {a.label}")
    say(f"  base significant: {len(sig)}   kept: {len(kept)} "
        f"({100*len(kept)/max(len(sig),1):.1f}%)   lost: {len(lost)}")
    say()

    if len(lost) < 5:
        say("  too few lost events to characterise")
        open(f"{a.out}.report.txt", "w").write("\n".join(lines) + "\n")
        return

    say(f"  {'group':<8}{'n':>6}{'|dPSI| base':>13}{'|dPSI| adj':>12}"
        f"{'ratio adj/base':>16}{'SE base':>10}{'SE adj':>10}{'SE ratio':>10}")
    for lab, t in (("kept", kept), ("lost", lost)):
        rat = (t.dPSI_adj_a.abs() / t.dPSI_adj_b.abs()).replace([np.inf, -np.inf], np.nan)
        se_r = (t.se_adj_a / t.se_adj_b).replace([np.inf, -np.inf], np.nan)
        say(f"  {lab:<8}{len(t):>6}{np.nanmedian(t.dPSI_adj_b.abs()):>13.3f}"
            f"{np.nanmedian(t.dPSI_adj_a.abs()):>12.3f}{np.nanmedian(rat):>16.3f}"
            f"{np.nanmedian(t.se_adj_b):>10.4f}{np.nanmedian(t.se_adj_a):>10.4f}"
            f"{np.nanmedian(se_r):>10.3f}")
    say()

    lr = (lost.dPSI_adj_a.abs() / lost.dPSI_adj_b.abs()).replace([np.inf, -np.inf], np.nan)
    lr = lr.dropna()
    med = float(np.median(lr))
    frac_collapsed = float(np.mean(lr < 0.5))
    same_sign = float(np.mean(np.sign(lost.dPSI_adj_a) == np.sign(lost.dPSI_adj_b)))
    say("--- verdict on the LOST events ---")
    say(f"  median effect retained : {med:.3f}  (1.0 = unchanged, 0 = collapsed)")
    say(f"  fraction losing >half their effect: {frac_collapsed:.1%}")
    say(f"  still same direction   : {same_sign:.1%}")
    say(f"  median SE inflation    : "
        f"{float(np.nanmedian(lost.se_adj_a / lost.se_adj_b)):.3f}x")
    say()
    if med > 0.8 and frac_collapsed < 0.2:
        say("  => THRESHOLD ATTRITION. The lost events kept their effect sizes; they")
        say("     simply fell below the cutoffs as the model got more conservative.")
        say("     Composition is NOT explaining them away. Report the adjusted model")
        say("     as primary and note that effect estimates are unchanged.")
    elif med < 0.5:
        say("  => CONFOUNDING. The lost events genuinely lose their effect once")
        say("     composition is accounted for. That fraction of the signature is")
        say("     tissue makeup and must be described as such.")
    else:
        say("  => MIXED. Part attenuation, part attrition. Report both models and")
        say("     the effect-size correlation; do not claim the signature is")
        say("     composition-independent.")
    say()

    # do the borderline events explain it?
    say("--- were the lost events simply the marginal ones? ---")
    for lab, t in (("kept", kept), ("lost", lost)):
        say(f"  {lab:<6} median base |dPSI| {np.nanmedian(t.dPSI_adj_b.abs()):.3f}   "
            f"median base FDR {np.nanmedian(t.FDR_b):.2e}")
    u, pu = stats.mannwhitneyu(kept.dPSI_adj_b.abs(), lost.dPSI_adj_b.abs(),
                               alternative="greater")
    say(f"  kept events had larger baseline effects: MWU p = {pu:.3g}")
    say("  (a small p here means the losses are concentrated at the margin, which")
    say("   is the signature of attrition rather than confounding)")
    say()

    genes = [g.strip() for g in a.genes.split(",") if g.strip()]
    gb = sig[sig.gene_b.isin(genes)] if "gene_b" in sig.columns else sig[sig.gene.isin(genes)]
    if len(gb):
        gcol = "gene_b" if "gene_b" in sig.columns else "gene"
        klost = gb[~((gb.FDR_a < a.fdr) & (gb.dPSI_adj_a.abs() >= a.min_dpsi))]
        say(f"--- RBM20 program: {len(gb) - len(klost)}/{len(gb)} retained "
            f"({100*(len(gb)-len(klost))/len(gb):.0f}%) vs "
            f"{100*len(kept)/len(sig):.0f}% genome-wide ---")
        if len(klost):
            say(f"  lost: {sorted(klost[gcol].tolist())}")

    sig["kept"] = keptmask.values
    sig.to_csv(f"{a.out}.events.tsv", sep="\t", index=False)
    with open(f"{a.out}.report.txt", "w") as fh:
        fh.write("\n".join(lines) + "\n")
    say("=" * 68)
    print(f"\nwrote {a.out}.report.txt and {a.out}.events.tsv")


if __name__ == "__main__":
    main()
