#!/usr/bin/env python3
"""
13_gene_lists_for_enrichment.py

Writes the gene lists for GO / pathway enrichment from the phase2_v2 tables.

THE BACKGROUND IS THE WHOLE BALLGAME
------------------------------------
Do NOT run enrichment against "all human genes". The set of genes that can be
tested for alternative splicing is biased toward long, highly expressed,
multi-exon genes. In heart tissue those are sarcomeric, cytoskeletal and
mitochondrial. Against a whole-genome background you will get gorgeous
"muscle contraction" and "sarcomere organization" terms that are an artefact
of what is testable, not of what is dysregulated.

This script therefore writes `*_background.txt`: every gene with at least one
event that PASSED the coverage and variance filters and was actually tested.
That is the correct universe. Use it.

Outputs (per --tag):
  <tag>_background.txt        universe: all genes tested            <-- REQUIRED
  <tag>_foreground.txt        significant genes, any direction
  <tag>_foreground_up.txt     genes whose lead event gains inclusion in DCM
  <tag>_foreground_down.txt   genes whose lead event loses inclusion
  <tag>_multi_event_type.txt  genes significant in >=2 event types (high confidence)
  <tag>_gene_table.tsv        per gene: n events, best FDR, max |dPSI|, types, direction
  <tag>_summary.txt           counts and the exact settings to use in g:Profiler

Usage (composition-adjusted tables as primary):
  python 13_gene_lists_for_enrichment.py --suffix _v2compadj --tag compadj
  python 13_gene_lists_for_enrichment.py --suffix _v2filtered --tag base
"""
import argparse
import os

import numpy as np
import pandas as pd

TYPES = ["SE", "MXE", "A3SS", "A5SS", "RI"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--suffix", default="_v2compadj",
                    help="table suffix, e.g. _v2compadj / _v2filtered / _v2fulladj")
    ap.add_argument("--types", default=",".join(TYPES))
    ap.add_argument("--fdr", type=float, default=0.05)
    ap.add_argument("--min-dpsi", type=float, default=0.10)
    ap.add_argument("--tag", default="enrich")
    a = ap.parse_args()

    types = [t.strip() for t in a.types.split(",") if t.strip()]
    frames = []
    missing = []
    for et in types:
        f = f"{et}{a.suffix}.all_events.tsv"
        if not os.path.exists(f):
            missing.append(f)
            continue
        t = pd.read_csv(f, sep="\t", usecols=lambda c: c in
                        {"gene", "event_type", "key", "dPSI_adj", "FDR", "chr",
                         "exonStart", "exonEnd"})
        t["event_type"] = et
        frames.append(t)
    if missing:
        print(f"[warn] missing: {missing}")
    if not frames:
        raise SystemExit("[fatal] no tables found -- check --suffix")

    df = pd.concat(frames, ignore_index=True)
    df = df.dropna(subset=["gene"])
    df = df[~df.gene.astype(str).str.startswith("ENSG")]   # keep named genes only

    # ---- background: every gene with at least one TESTED event ----------
    background = sorted(df.gene.astype(str).unique())

    sig = df[(df.FDR < a.fdr) & (df.dPSI_adj.abs() >= a.min_dpsi)].copy()

    # ---- per-gene summary -----------------------------------------------
    rows = []
    for g, gg in sig.groupby("gene"):
        lead = gg.loc[gg.dPSI_adj.abs().idxmax()]
        rows.append({
            "gene": g,
            "n_events": len(gg),
            "n_event_types": gg.event_type.nunique(),
            "event_types": ",".join(sorted(gg.event_type.unique())),
            "best_FDR": gg.FDR.min(),
            "max_abs_dPSI": gg.dPSI_adj.abs().max(),
            "lead_dPSI": lead.dPSI_adj,
            "direction": "up" if lead.dPSI_adj > 0 else "down",
            "mixed_direction": (gg.dPSI_adj > 0).any() and (gg.dPSI_adj < 0).any(),
        })
    gt = pd.DataFrame(rows).sort_values(["n_event_types", "max_abs_dPSI"],
                                        ascending=[False, False])

    def write(fn, items):
        with open(fn, "w") as fh:
            fh.write("\n".join(items) + "\n")
        return fn

    out = []
    out.append(write(f"{a.tag}_background.txt", background))
    out.append(write(f"{a.tag}_foreground.txt", sorted(gt.gene)))
    out.append(write(f"{a.tag}_foreground_up.txt",
                     sorted(gt[gt.direction == "up"].gene)))
    out.append(write(f"{a.tag}_foreground_down.txt",
                     sorted(gt[gt.direction == "down"].gene)))
    out.append(write(f"{a.tag}_multi_event_type.txt",
                     sorted(gt[gt.n_event_types >= 2].gene)))
    gt.to_csv(f"{a.tag}_gene_table.tsv", sep="\t", index=False,
              float_format="%.6g")
    out.append(f"{a.tag}_gene_table.tsv")

    lines = []
    def say(s=""):
        lines.append(s)
        print(s)

    say("=" * 68)
    say(f"GENE LISTS FOR ENRICHMENT  (tables: *{a.suffix}, FDR<{a.fdr}, "
        f"|dPSI|>={a.min_dpsi})")
    say()
    say(f"  events tested (all types) : {len(df)}")
    say(f"  BACKGROUND genes tested   : {len(background)}   <-- use this as universe")
    say(f"  significant events        : {len(sig)}")
    say(f"  FOREGROUND genes          : {len(gt)}")
    say(f"    increased inclusion     : {(gt.direction == 'up').sum()}")
    say(f"    decreased inclusion     : {(gt.direction == 'down').sum()}")
    say(f"    mixed within gene       : {int(gt.mixed_direction.sum())}")
    say(f"  genes in >=2 event types  : {int((gt.n_event_types >= 2).sum())}")
    say()
    say("  per event type:")
    for et in types:
        s = sig[sig.event_type == et]
        b = df[df.event_type == et]
        say(f"    {et:<5} {s.gene.nunique():>5} significant genes "
            f"of {b.gene.nunique():>6} tested")
    say()
    say("--- TOP 25 GENES (by breadth then effect size) ---")
    say(f"  {'gene':<12}{'types':>6}{'events':>7}{'maxdPSI':>9}{'bestFDR':>11}  types")
    for _, r in gt.head(25).iterrows():
        say(f"  {r.gene:<12}{r.n_event_types:>6}{r.n_events:>7}"
            f"{r.lead_dPSI:>+9.3f}{r.best_FDR:>11.2e}  {r.event_types}")
    say()
    say("--- HOW TO RUN THE ENRICHMENT ---")
    say("  g:Profiler (https://biit.cs.ut.ee/gprofiler/gost) -- recommended,")
    say("  because it accepts a custom background properly:")
    say(f"    1. paste {a.tag}_foreground.txt into the query box")
    say("    2. open 'Advanced options'")
    say("    3. Statistical domain scope -> 'Custom over annotated genes'")
    say(f"    4. paste {a.tag}_background.txt into the background box")
    say("    5. Organism: Homo sapiens. Significance threshold: g:SCS.")
    say("    Run the up and down lists separately as well -- direction-specific")
    say("    enrichment is far more interpretable than the union.")
    say()
    say("  DAVID also accepts a background (upload list, then 'Use as background').")
    say("  ENRICHR DOES NOT support a custom background -- do not use it here,")
    say("  or your terms will be the testability artefact described above.")
    say()
    say("  In R: clusterProfiler::enrichGO(gene=fg, universe=bg, ...) -- the")
    say("  'universe' argument is the background and is NOT optional here.")
    say()
    say("--- SANITY CHECK BEFORE YOU BELIEVE ANY RESULT ---")
    say("  Run the SAME enrichment with the whole genome as background. If the")
    say("  top terms are muscle/sarcomere/mitochondrion in that run but change")
    say("  with the proper background, you have just demonstrated the artefact")
    say("  and you should say so in the methods.")
    say("=" * 68)

    write(f"{a.tag}_summary.txt", lines)
    print(f"\nwrote: {', '.join(out)}, {a.tag}_summary.txt")


if __name__ == "__main__":
    main()
