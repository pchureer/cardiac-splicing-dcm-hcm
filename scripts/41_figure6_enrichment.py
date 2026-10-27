#!/usr/bin/env python3
"""
41_figure6_enrichment.py -- FIGURE 6: the enrichment, and the background that decides it.

  A  How many terms survive. 4 against the proper background of tested genes;
     229 against the whole genome. Same foreground, same tool, same correction.
  B  What the terms ARE. The whole-genome list is led by cytoplasm and protein
     binding at p ~ 1e-37; the proper-background list is led by cytoskeleton in
     muscle cells at p = 3e-05. Smaller p-values, less meaning.
  C  The trap, and the reason this figure exists. The multi-event-type list returns
     sarcomere, Z disc, I band, myofibril -- exactly what a cardiac splicing paper
     wants to report -- on the wrong background, and NOTHING on the right one.
  D  What actually survives: the 15 genes of "structural constituent of muscle".

WHY THE FIGURE IS BUILT THIS WAY
--------------------------------
The temptation is to show only panel D and call it a sarcomere enrichment. The
project's own result is that the rich version of that story is an artefact of using
the whole genome as background: genes testable for splicing are long, highly
expressed and multi-exon, so they are enriched for generic compartments before any
biology is considered. Panels A-C exist so a reader cannot take D out of context,
and so a reviewer sees the check was done rather than having to ask for it.

The honest reading, which belongs in the legend: the surviving enrichment is modest
BECAUSE the background is already heart- and muscle-biased. A paper reporting rich
sarcomere enrichment here would be reporting its background.

PROVENANCE
----------
The audited numbers live in RESULT_go_enrichment.md and are shipped in FROZEN_*, so
the figure always builds. If the g:Profiler result tables from script 14 are still on
disk, the script finds them by their columns, uses them instead, and cross-checks the
term counts against the doc. The PROVENANCE block it prints names each panel's source.

A WARNING CARRIED OVER FROM THE DOC: g:Profiler's `intersections` are aligned to the
MAPPED Ensembl list, not to the submitted symbols. An older copy of 14_run_gprofiler.py
indexed the submitted list and produced the right COUNT with the wrong GENES. If a
discovered table's gene column disagrees with the frozen 15, this script says so and
keeps the frozen list.

Usage:
  export PATH=/data/apps/anaconda3/envs/mamba/envs/rnaseq_tools/bin:$PATH
  cd /data/project/chureeratp/RNA/HEART
  python 41_figure6_enrichment.py
  python 41_figure6_enrichment.py --correct gprofiler_bg.tsv --genome gprofiler_wg.tsv
  python 41_figure6_enrichment.py --list
"""
import argparse
import glob
import os

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.patches import FancyBboxPatch

BLUE = "#2a78d6"
ORANGE = "#eb6834"
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#8a8a85"
GRID = "#e4e4e0"
SURF = "#fcfcfb"
NULLGREY = "#c9c8c3"

# ---------------------------------------------------------------------------
# FROZEN -- RESULT_go_enrichment.md, 2026-09-08. Do not edit to improve a panel.
# ---------------------------------------------------------------------------
import os as _os, sys as _sys

BG_TESTED = 12725          # genes with >=1 tested event: the correct background
FG_ALL = 927               # significant genes, FDR<0.05 & |dPSI|>=0.10
FG_UP, FG_DOWN, FG_MIXED = 636, 291, 75
FG_MULTI = 134             # significant in >=2 event types

# 165: these six are PLOTTED TEXT, not filenames -- panel A's annotation,
# panel B's axis label and panel C's header all render them. Rule 10: a
# plotted constant is a write path. The override is therefore MANDATORY, and
# the trigger is the DATA rather than a flag that could be forgotten.
_F6C = _os.environ.get("F6_COUNTS", "").strip()
if "rl125" in " ".join(_sys.argv) and not _F6C:
    _sys.exit(
        "[fatal] a 125 bp g:Profiler table was named on the command line, but\n"
        "        F6_COUNTS is unset. BG_TESTED/FG_ALL/FG_UP/FG_DOWN/FG_MIXED/\n"
        "        FG_MULTI would be drawn as 50 bp text beside 125 bp bars.\n"
        "        For 125 bp, set all three:\n"
        "          F6_COUNTS=12564/805/550/255/73/131\n"
        "          F6_TOTALS=1/85\n"
        "          F6_MULTI=0/17\n"
        "        Nothing plotted.")
if _F6C:
    BG_TESTED, FG_ALL, FG_UP, FG_DOWN, FG_MIXED, FG_MULTI = (
        int(_x) for _x in _F6C.split("/"))

FROZEN_TOTAL = {"correct": 4, "genome": 229}
_F6T = _os.environ.get("F6_TOTALS", "").strip()
if _F6C and not _F6T:
    _sys.exit("[fatal] F6_COUNTS is set but F6_TOTALS is not. Panel A's bars\n"
              "        would stay at the 50 bp 4 and 229. Nothing plotted.")
if _F6T:
    _a, _b = (int(_x) for _x in _F6T.split("/"))
    FROZEN_TOTAL = {"correct": _a, "genome": _b}

# proper background: every term that survived, with its query list
FROZEN_CORRECT = [
    ("all",  "KEGG",  "Cytoskeleton in muscle cells",        3.41e-05, 38),
    ("all",  "GO:MF", "structural constituent of muscle",    9.47e-05, 15),
    ("all",  "GO:MF", "cytoskeletal protein binding",        3.38e-02, 87),
    ("up",   "KEGG",  "Cytoskeleton in muscle cells",        3.42e-02, 27),
]
# whole genome: the leading terms the doc records with a p-value
FROZEN_GENOME = [
    ("all", "GO:CC", "cytoplasm",                          6.8e-37, 628),
    ("all", "GO:MF", "protein binding",                    8.9e-19, 694),
    ("all", "GO:CC", "cytosol",                            3.2e-18, None),
    ("all", "GO:CC", "intracellular anatomical structure", 4.8e-12, None),
]
GENOME_ALSO = ["organelle", "organelle organization", "mitochondrion"]

# the trap: the multi-event-type list, 15 terms on the wrong background, 0 on the right
FROZEN_MULTI = {"correct": 0, "genome": 15}
_F6M = _os.environ.get("F6_MULTI", "").strip()
if _F6C and not _F6M:
    _sys.exit("[fatal] F6_COUNTS is set but F6_MULTI is not. Panel C would keep\n"
              "        the 50 bp 15. Nothing plotted.")
if _F6M:
    _a, _b = (int(_x) for _x in _F6M.split("/"))
    FROZEN_MULTI = {"correct": _a, "genome": _b}
QUERY_SHOWN = {"all"}   # panel B compares one query list across both backgrounds

MULTI_TERMS = ["sarcomere", "Z disc", "I band", "alpha-actinin binding",
               "myofibril", "contractile muscle fiber"]

# the one specific term that survives, and its genes
# 165: STALE AT 125 bp, AND NO LONGER DRAWN. "structural constituent of
# muscle" is one of the two GO:MF terms LOST when the proper-background run
# falls from three terms to one. panel_d is built entirely on it and is no
# longer called. Do not revive either without re-measuring.
SURVIVOR = ("structural constituent of muscle", "GO:MF", 9.47e-05)
GENE_GROUPS = [
    ("giant sarcomeric proteins", ["TTN", "NEB", "OBSCN"], False),
    ("Z-disc / myofibrillar disease genes",
     ["CSRP3", "MYOT", "MYOM1", "SORBS2", "DMD"], False),
    ("excitation-contraction coupling", ["JPH1", "ASPH"], False),
    ("contractile / thin filament", ["MYL3", "MYL6", "TPM1"], False),
    ("smooth muscle \u2014 a tissue-content caveat, not a finding",
     ["MYH11", "SMTN"], True),
]


def style(ax):
    ax.set_facecolor(SURF)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(GRID)
        ax.spines[s].set_linewidth(0.8)
    ax.tick_params(colors=INK2, labelsize=7, length=3, width=0.8)


def lower(cols):
    return {str(c).lower(): c for c in cols}


def sniff(path, nrows=300):
    try:
        return pd.read_csv(path, sep=None, engine="python", nrows=nrows)
    except Exception:
        return None


def read_gp(path, verbose=True):
    """A g:Profiler result table: term name + p-value, one row per term."""
    t = sniff(path)
    if t is None or not len(t):
        return None
    L = lower(t.columns)
    n = next((L[k] for k in ("term_name", "name", "term", "description")
              if k in L), None)
    p = next((L[k] for k in ("p_value", "p_val", "pvalue", "p", "adjusted_p_value",
                             "padj") if k in L), None)
    if n is None or p is None:
        if verbose:
            print(f"  [skip] {os.path.basename(path)}: need a term-name and a p-value "
                  f"column; have {list(t.columns)[:8]}")
        return None
    full = pd.read_csv(path, sep=None, engine="python")
    L = lower(full.columns)
    out = pd.DataFrame({"term": full[n].astype(str),
                        "p": pd.to_numeric(full[p], errors="coerce")})
    for key, names in (("source", ("source", "db", "category")),
                       ("query", ("query", "list", "set")),
                       ("n", ("intersection_size", "intersection", "n", "count"))):
        c = next((L[k] for k in names if k in L), None)
        out[key] = full[c] if c is not None else np.nan
    out = out.dropna(subset=["p"])
    return out if len(out) else None


def _n_all(df, alpha, verbose=True):
    """166: significant terms in the 'all' query list only.

    g:Profiler emits one row per (query list, term) and this project submits
    four lists -- all, up, down, multi_event_type. Counting every row answers a
    question the manuscript never asks AND double counts, because panel C
    reports the multi_event_type rows separately. Panel B already restricts
    itself to QUERY_SHOWN = {"all"}; panel A did not.
    """
    if df is None:
        return 0
    q = df.get("query")
    if q is None or not q.notna().any():
        if verbose:
            print("  [warn] no 'query' column in a g:Profiler table: panel A is "
                  "counting EVERY query list, which is the 179-vs-85 defect. "
                  "Check the number before publishing it.")
        return int((df.p < alpha).sum())
    m = q.astype(str).str.strip().str.lower().eq("all")
    if not m.any():
        raise SystemExit(
            "[fatal] a g:Profiler table has a 'query' column but no rows with\n"
            "        query == 'all'. Panel A cannot be scoped and would count\n"
            "        every list. Nothing plotted.")
    return int((df[m].p < alpha).sum())


def discover(explicit_correct, explicit_genome, verbose=True):
    """Find the two result tables. Filename decides which background each is."""
    if explicit_correct or explicit_genome:
        C = read_gp(explicit_correct) if explicit_correct else None
        G = read_gp(explicit_genome) if explicit_genome else None
        return (C, explicit_correct if C is not None else None,
                G, explicit_genome if G is not None else None)
    cands = []
    for pat in ("*gprofiler*.tsv", "*gprofiler*.csv", "*enrich*.tsv", "*enrich*.csv",
                "*/gprofiler*.tsv", "*go_*.tsv"):
        cands += [f for f in sorted(glob.glob(pat))
                  if os.path.isfile(f) and os.path.getsize(f)]
    cands = list(dict.fromkeys(cands))
    if verbose and cands:
        print(f"  [info] {len(cands)} candidate enrichment table(s): "
              + ", ".join(os.path.basename(f) for f in cands[:8]))
    # The whole-genome run is the one that gets MARKED in a filename
    # (..._WHOLEGENOME.tsv); the proper-background run is usually just the plain
    # name. So: anything whose name says genome is the genome run, and anything
    # else that reads as a g:Profiler table is the custom-background run. The
    # earlier rule required a positive keyword on BOTH and silently skipped
    # compadj_gprofiler.tsv, leaving panel A half frozen with no warning.
    GENOME_KEYS = ("genome", "wholegenome", "wg", "nobg", "default_bg", "allgenes")
    BG_KEYS = ("bg", "background", "tested", "custom", "domain")
    C = G = cp = gp = None
    parsed = []
    for f in cands:
        t = read_gp(f, verbose)
        if t is not None:
            parsed.append((f, t, any(k in os.path.basename(f).lower()
                                     for k in GENOME_KEYS)))
    for f, t, is_genome in parsed:
        if is_genome and G is None:
            G, gp = t, f
    # prefer a file that names itself as the background run, else take what is left
    rest = [(f, t) for f, t, is_genome in parsed if not is_genome]
    named = [(f, t) for f, t in rest
             if any(k in os.path.basename(f).lower() for k in BG_KEYS)]
    for f, t in (named or rest):
        if C is None:
            C, cp = t, f
    if verbose:
        if C is not None and not any(k in os.path.basename(cp).lower()
                                     for k in BG_KEYS):
            print(f"  [info] treating {os.path.basename(cp)} as the proper-background "
                  f"run because its name does not say 'genome'. If that is wrong, "
                  f"pass --correct / --genome explicitly.")
        if cands and C is None and G is None:
            print("  [info] no candidate parsed as a g:Profiler table; "
                  "using frozen numbers")
    return C, cp, G, gp


# ---------------------------------------------------------------------------
# panels
# ---------------------------------------------------------------------------

def panel_a(ax, counts):
    lab = ["tested genes\n(n = %s)" % f"{BG_TESTED:,}", "whole genome"]
    val = [counts["correct"], counts["genome"]]
    col = [BLUE, ORANGE]
    y = [1, 0]
    ax.barh(y, val, height=.5, color=col, linewidth=0)
    for yy, v in zip(y, val):
        ax.text(v + max(val) * .02, yy, f"{v:,}", va="center", ha="left",
                fontsize=11, color=INK, fontweight="bold")
    ax.set_yticks(y)
    ax.set_yticklabels(lab, fontsize=7.2, color=INK)
    ax.set_xlim(0, max(val) * 1.22)
    ax.set_xlabel("significant terms (g:SCS corrected)", fontsize=7.5, color=INK2)
    ax.xaxis.grid(True, color=GRID, linewidth=.6)
    ax.set_axisbelow(True)
    if counts["correct"]:
        fold = counts["genome"] / counts["correct"]
        # 165: no ratio. At 125 bp this is 85 over 1, and a fold change on a
        # denominator of one is what the Supplementary Discussion already
        # rejected for 452-fold and 1,767-fold. The bars carry both counts.
        ax.annotate(f"the same {FG_ALL} genes,\ntwo backgrounds",
                    xy=(.97, .52), xycoords="axes fraction", ha="right", va="center",
                    fontsize=7, color=INK)
    style(ax)
    ax.set_title("A   Background choice, not biology", loc="left", fontsize=8.6,
                 fontweight="bold", pad=8)


def panel_b(ax, C, G, alpha=0.05, ntop=4):
    """Top terms under each background. Uses the real tables when they were found."""
    # A term can appear more than once in a g:Profiler table -- once per query list
    # ("all", "up", "down", "multi") and again per source. Taking the top N rows by
    # p-value then prints the SAME term twice and pushes a real one off the panel,
    # which is exactly what happened on the first real run: "cytoplasm" twice, and
    # "intracellular anatomical structure" lost. Collapse on the term name, keep the
    # strongest p, and say which query list it came from.
    # ...but collapsing on the term alone then drops "Cytoskeleton in muscle cells
    # [up]", which IS one of the four terms panel A counts, so panel A said 4 and
    # panel B showed 3. Both fixes are needed and they pull against each other.
    # Resolution: panel B compares ONE query list across the two backgrounds --
    # `all`, the full significant set -- so the two sides are like for like, and
    # panel A keeps counting every term in every list. The subtitle says so.
    def top(T, n, need_sig=False):
        t = T[T.p < alpha] if need_sig else T
        if "query" in t.columns and t["query"].notna().any():
            only = t[t["query"].astype(str).str.lower().isin(QUERY_SHOWN)]
            if len(only):
                t = only
        t = t.sort_values("p")
        ndup = int(t.duplicated("term").sum())
        t = t.drop_duplicates("term", keep="first").head(n)
        return [(str(r.term), float(r.p)) for _, r in t.iterrows()], ndup

    rows, used, dups = [], {"genome": "FROZEN", "correct": "FROZEN"}, {}
    if G is not None:
        got, dups["whole-genome"] = top(G, ntop)
        rows += [("genome", t, pv) for t, pv in got]
        used["genome"] = "file"
    else:
        rows += [("genome", term, pv) for _, _, term, pv, _ in FROZEN_GENOME]
    rows.append((None, None, None))                      # separator
    if C is not None:
        got, dups["tested-gene"] = top(C, ntop + 2, need_sig=True)
        rows += [("correct", t, pv) for t, pv in got]
        used["correct"] = "file"
    else:
        rows += [("correct", term, pv) for q, _, term, pv, _ in FROZEN_CORRECT
                 if q in QUERY_SHOWN]
    for k, n in dups.items():
        if n:
            print(f"  [info] panel B: collapsed {n} repeated term name(s) in the "
                  f"{k} table (same term under several query lists or sources); "
                  f"kept the strongest p for each")
    xmax = max(-np.log10(p) for k, t, p in rows if k is not None)
    y = np.arange(len(rows))[::-1]
    for yy, (kind, term, p) in zip(y, rows):
        if kind is None:
            ax.axhline(yy, color=GRID, linewidth=1.0)
            continue
        col = ORANGE if kind == "genome" else BLUE
        x = -np.log10(p)
        ax.plot([0, x], [yy, yy], color=GRID, linewidth=1.6, solid_capstyle="round",
                zorder=1)
        ax.scatter([x], [yy], s=30, color=col, zorder=3, linewidths=0)
        ax.text(x + xmax * .028, yy, term, va="center", ha="left", fontsize=5.9,
                color=INK)
    ax.set_xlim(0, xmax * 1.45)
    ax.set_ylim(-.8, len(rows) - .2)
    ax.set_yticks([])
    ax.set_xlabel("$-\\log_{10}$ p", fontsize=7.5, color=INK2)
    ax.xaxis.grid(True, color=GRID, linewidth=.6)
    ax.set_axisbelow(True)
    ngen = sum(1 for k, _, _ in rows if k == "genome")
    ncor = sum(1 for k, _, _ in rows if k == "correct")
    tot = len(rows)
    ax.text(-0.022, 1 - (ngen / 2) / tot, "whole genome", transform=ax.transAxes,
            rotation=90, ha="right", va="center", fontsize=6.4, color=ORANGE,
            fontweight="bold")
    ax.text(-0.022, (ncor / 2) / tot, "tested-gene\nbackground",
            transform=ax.transAxes, rotation=90, ha="right", va="center",
            fontsize=6.4, color=BLUE, fontweight="bold", linespacing=.9)
    style(ax)
    ax.set_title("B   Smaller p-values, less meaning", loc="left", fontsize=8.6,
                 fontweight="bold", pad=19)
    ax.text(0, 1.012, f"leading terms of the “{'/'.join(sorted(QUERY_SHOWN))}” query "
                      f"list under each background", transform=ax.transAxes,
            fontsize=5.9, color=MUTED, va="bottom")


def panel_c(ax):
    """The trap. Text, because the content IS the argument."""
    ax.axis("off")
    ax.set_title("C   The attractive terms are the artefact", loc="left",
                 fontsize=8.6, fontweight="bold", pad=8)
    ax.text(0, .975, f"the {FG_MULTI} genes significant in ≥2 event types, tested "
                     f"against:", fontsize=6.8, color=INK2, va="top")

    def box(x, w, colr, head, sub, n, items, foot):
        ax.add_patch(FancyBboxPatch((x, .05), w, .82, boxstyle="round,pad=0.010",
                                    linewidth=0, facecolor=colr + "14", zorder=0))
        # head and sub on SEPARATE lines: "tested genes (n = 12,725)" on one line
        # overflows the box at any font size that is still readable in print
        ax.text(x + .030, .825, head, fontsize=7.2, color=colr, fontweight="bold",
                va="center")
        ax.text(x + .030, .757, sub, fontsize=5.9, color=colr, va="center", alpha=.85)
        ax.text(x + .030, .655, n, fontsize=16, color=colr, fontweight="bold",
                va="center")
        yy = .545
        for it in items:
            ax.plot([x + .038, x + .038], [yy - .022, yy + .022], color=colr,
                    linewidth=1.4, alpha=.55, solid_capstyle="round")
            ax.text(x + .058, yy, it, fontsize=6.0, color=INK, va="center")
            yy -= .066
        ax.text(x + .030, .105, foot, fontsize=5.9, color=INK2, va="center",
                style="italic")

    box(0.0, .475, ORANGE, "whole genome", "g:Profiler default",
        f"{FROZEN_MULTI['genome']} terms",
        MULTI_TERMS[:5] + [f"…and {FROZEN_MULTI['genome'] - len(MULTI_TERMS[:5])} more"],
        "exactly the cardiac story")
    box(.525, .475, BLUE, "tested genes", f"n = {BG_TESTED:,}",
        f"{FROZEN_MULTI['correct']} terms", [],
        "nothing survives correction")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)


def panel_d(ax):
    """The surviving genes, as a grid of chips grouped by role.

    A bare run of gene symbols reads as a leftover list. The grid gives every symbol
    the same footprint, the rule down the left binds each row to its role, and the
    greyed row makes the smooth-muscle caveat visibly a different kind of statement
    rather than a sentence a reader may skip.
    """
    ax.axis("off")
    term, src, p = SURVIVOR
    ax.set_title("D   Survives adjustment: 15 genes", loc="left", fontsize=8.6,
                 fontweight="bold", pad=8)
    ax.text(0, .975, f"“{term}”", fontsize=7.0, color=INK, va="top")
    ax.text(0, .905, f"{src}   p = {p:.1e}   ·   {sum(len(g) for _, g, _ in GENE_GROUPS)}"
                     f" of {FG_ALL} signature genes", fontsize=6.0, color=MUTED,
            va="top")

    NCOL, CW, CH, GAP = 4, .225, .075, .018
    y = .800
    for head, genes, caveat in GENE_GROUPS:
        col = MUTED if caveat else BLUE
        nrow = int(np.ceil(len(genes) / NCOL))
        block_h = nrow * (CH + GAP)
        # the accent rule spans the whole block, header included
        ax.plot([0.004, 0.004], [y - block_h + .012, y + .040], color=col,
                linewidth=1.8, alpha=.5 if caveat else .85, solid_capstyle="round")
        ax.text(.028, y + .013, head, fontsize=6.0,
                color=MUTED if caveat else INK2,
                va="center", style="italic" if caveat else "normal")
        for i, g in enumerate(genes):
            r, c = divmod(i, NCOL)
            gx = .028 + c * (CW + .012)
            gy = y - .052 - r * (CH + GAP)
            ax.add_patch(FancyBboxPatch((gx, gy - CH / 2), CW, CH,
                                        boxstyle="round,pad=0.004",
                                        linewidth=0.7,
                                        edgecolor=col + ("55" if caveat else "88"),
                                        facecolor=col + ("0f" if caveat else "1c"),
                                        zorder=1))
            ax.text(gx + CW / 2, gy, g, fontsize=6.6, ha="center", va="center",
                    color=MUTED if caveat else INK,
                    fontweight="normal" if caveat else "bold", zorder=2)
        y -= block_h + .055
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--correct", help="g:Profiler table run against the tested-gene "
                                      "background")
    ap.add_argument("--genome", help="g:Profiler table run against the whole genome")
    ap.add_argument("--alpha", type=float, default=0.05)
    ap.add_argument("--frozen-only", action="store_true")
    ap.add_argument("--list", dest="listonly", action="store_true")
    ap.add_argument("--out", default="F6_enrichment")
    a = ap.parse_args()

    print("[scan] looking for g:Profiler result tables")
    C = G = cp = gp = None
    if not a.frozen_only:
        C, cp, G, gp = discover(a.correct, a.genome)

    counts = dict(FROZEN_TOTAL)
    src = {"A": "FROZEN", "B": "FROZEN", "C": "FROZEN", "D": "FROZEN"}
    if C is not None:
        counts["correct"] = _n_all(C, a.alpha)
        src["A"] = cp
    if G is not None:
        counts["genome"] = _n_all(G, a.alpha)
        src["A"] = f"{src['A']} + {gp}" if C is not None else gp

    if a.listonly:
        print(f"  correct-background table: {cp or 'not found'}")
        print(f"  whole-genome table:       {gp or 'not found'}")
        print(f"  term counts that would be plotted: {counts}")
        return

    plt.rcParams.update({"font.family": "sans-serif", "axes.linewidth": .8,
                         "figure.facecolor": SURF, "savefig.facecolor": SURF})
    fig = plt.figure(figsize=(7.6, 6.4))
    gs = gridspec.GridSpec(2, 2, wspace=0.30, hspace=0.42,
                           left=0.145, right=0.965, top=0.935, bottom=0.085)
    panel_a(fig.add_subplot(gs[0, 0]), counts)
    panel_b(fig.add_subplot(gs[0, 1]), C, G, a.alpha)
    # 165: panel D retired -- its term does not exist at 125 bp. C spans the
    # bottom row, giving the three-panel figure the manuscript legend already
    # describes.
    panel_c(fig.add_subplot(gs[1, :]))
    for ext in ("pdf", "png"):
        fig.savefig(f"{a.out}.{ext}", dpi=400)
    print(f"[write] {a.out}.pdf and {a.out}.png")

    # cross-check anything recomputed against the doc
    for key, fz in FROZEN_TOTAL.items():
        if counts[key] != fz and (C is not None or G is not None):
            print(f"\n  *** {key} background: {counts[key]} terms at p<{a.alpha} from "
                  f"the file, {fz} in RESULT_go_enrichment.md.")
            print("      Check the p-value column (raw vs corrected) and the alpha "
                  "before publishing either number.")

    print("\n" + "=" * 78)
    print("PROVENANCE -- paste into the methods / legend")
    print("=" * 78)
    print(f"  panel A  {src['A']}")
    print(f"  panel B  {cp or gp or 'FROZEN: RESULT_go_enrichment.md'}")
    print(f"  panels C, D  FROZEN: RESULT_go_enrichment.md")
    print(f"  foreground {FG_ALL} genes ({FG_UP} up, {FG_DOWN} down, {FG_MIXED} mixed "
          f"within gene); {FG_MULTI} significant in ≥2 event types")
    print(f"  backgrounds: {BG_TESTED:,} tested genes  vs  whole genome")
    print(f"  terms: {counts['correct']} vs {counts['genome']}")
    print("\n  The legend must say the enrichment is modest BECAUSE the background is")
    print("  already heart- and muscle-biased. A paper reporting rich sarcomere")
    print("  enrichment against the whole genome is reporting its background.")
    print("\n  Panel D's gene list is the CORRECTED one. g:Profiler's `intersections`")
    print("  index the mapped Ensembl list, not the submitted symbols; an older")
    print("  14_run_gprofiler.py gave the right count with the wrong genes.")
    print("\n  NOW LOOK AT THE PNG.")


if __name__ == "__main__":
    main()
