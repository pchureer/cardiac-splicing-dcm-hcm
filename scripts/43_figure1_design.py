#!/usr/bin/env python3
"""
43_figure1_design.py -- FIGURE 1: study design and sample flow.

  A  Where the samples went. Seven BioProjects, the genotype-identity screen, and the
     510 analysed. The 24 MAGNet<->LBNL same-individual pairs are drawn as a link
     BETWEEN the two cohorts, not written in a footnote -- it is the one thing about
     this study's composition a reader must not miss.
  B  What was compared. Every contrast in the paper with its post-exclusion group
     sizes, grouped by the claim it supports.

WHY THE OVERLAP IS DRAWN AND NOT WRITTEN
----------------------------------------
Earlier drafts of this work described four independent replication cohorts. Twenty-four
of LBNL's 33 samples turned out to be the same individuals as MAGNet samples, so there
are three. A reader who takes "independent" at face value will mis-weigh every
replication number in the paper. A drawn link is checked in one second; a footnote is
skipped.

CONSISTENCY WITH TABLE 1
------------------------
Every count here also appears in Table 1. If Table1.tsv is present the script parses it
and verifies each number, printing a loud warning on any disagreement -- the two must
never drift apart, and they will if one is edited by hand.

Usage:
  export PATH=/data/apps/anaconda3/envs/mamba/envs/rnaseq_tools/bin:$PATH
  cd /data/project/chureeratp/RNA/HEART
  python 43_figure1_design.py
"""
import argparse
import os
import re

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

BLUE = "#2a78d6"
ORANGE = "#eb6834"
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#8a8a85"
GRID = "#e4e4e0"
SURF = "#fcfcfb"
PALE = "#c9c8c3"

# ---------------------------------------------------------------------------
# FROZEN, and identical to Table 1. Sources: PROJECT_STATE_handoff.md,
# RESULT_sample_identity_LBNL_overlap.md.
#   name, deposited n, analysed n, platform, kind
# kind: disc = discovery, dcm = DCM replication, hcm = HCM replication,
#       prep = prep-robustness only, none = not analysed
# ---------------------------------------------------------------------------
COHORTS = [
    ("MAGNet",   360, 360, "paired 125 bp", "disc"),
    ("Colorado",  64,  62, "single 51 bp",  "dcm"),
    ("LBNL",      33,   9, "single 50–51 bp","prep"),
    ("Brigham",   23,  23, "paired 101 bp", "hcm"),
    ("Fuwai",     37,  37, "paired 150 bp", "hcm"),
    ("Stanford",  20,  19, "paired 150 bp", "hcm"),
]
# Heinig (PRJEB8360, 128 DCM) is NOT drawn. It has no control group, so it contributes
# no dPSI contrast and was never analysed. It was also not the only dataset in that
# position: Amsterdam and Naples were genotyped too and are likewise unused. Drawing
# one unused dataset and not the others is less honest than drawing none, so the
# figure shows what was analysed and the Methods say in one sentence which public
# datasets were retrieved and set aside, and why.
GENOTYPED = 670
SNVS = 3820
N_OVERLAP = 24
N_REPLICATE = 3
N_EXCLUDED = 27
N_ANALYSED = 510

# ---- 169: the matched-pair slopes are PLOTTED LABEL TEXT, not filenames ----
# Rule 10. --out renames the file; it does not touch these two numbers, so a
# 125 bp filename can carry 50 bp slopes with every check passing. The override
# is therefore MANDATORY, and symmetric: a 125 bp name without the slopes is
# refused, and the slopes without a 125 bp name are refused too, because that
# would overwrite the published 50 bp figure in place.
#
# 0.690/0.968 are 36_matched_pair_analysis.py's estimates, which the Results
# quote. 82 gives the same points with different bootstrap intervals; the two
# must not be mixed (DEFECT_matched_panel_rl125_has_no_file.md).
import os as _os
import sys as _sys

_F1S = _os.environ.get("F1_SLOPES", "").strip()
_F1_RL125 = "rl125" in " ".join(_sys.argv)
if _F1_RL125 and not _F1S:
    _sys.exit(
        "[fatal] a 125 bp output was named, but F1_SLOPES is unset.\n"
        "        Figure 1 draws the matched-pair slopes as LABEL TEXT, so the\n"
        "        figure would carry 0.584 and 0.988 under a 125 bp filename.\n"
        "        Re-run as:\n"
        "          F1_SLOPES=0.690/0.968 python3 43_figure1_design.py \\\n"
        "              --out F1_design_rl125\n"
        "        Use 36's estimates, which the Results quote, not 82's.\n"
        "        Nothing plotted.")
if _F1S and not _F1_RL125:
    _sys.exit(
        "[fatal] F1_SLOPES is set but no 125 bp output was named. The 125 bp\n"
        "        slopes would be drawn into F1_design, replacing the published\n"
        "        50 bp figure in place. Pass --out F1_design_rl125.\n"
        "        Nothing plotted.")
if _F1S:
    _F1_PLATFORM, _F1_PEOPLE = (float(_x) for _x in _F1S.split("/"))
else:
    _F1_PLATFORM, _F1_PEOPLE = 0.584, 0.988
print(f"[169] matched-pair slopes drawn: platform {_F1_PLATFORM:.3f}, "
      f"people {_F1_PEOPLE:.3f}"
      f"{'  (125 bp)' if _F1S else '  (50 bp defaults)'}")

CONTRASTS = [
    ("Discovery", [
        ("MAGNet", "DCM 166  vs  NF 166", "disc"),
    ]),
    ("DCM replication", [
        ("Colorado", "DCM 36  vs  NF 13", "dcm"),
    ]),
    ("HCM replication", [
        ("Brigham", "HCM 18  vs  ctrl 5", "hcm"),
        ("Fuwai", "HCM 28  vs  donor 9", "hcm"),
        ("Stanford", "HCM 13  vs  donor 6", "hcm"),
    ]),
    ("Etiology specificity", [
        ("Colorado", "DCM 36  vs  ICM 13", "dcm"),
        ("MAGNet", "DCM 166  vs  HCM 28", "disc"),
    ]),
    # The LBNL 5 vs 4 contrast is NOT where the platform result comes from. The
    # reported platform slope is the matched-pair panel: the 24 shared
    # individuals, split
    # 10 DCM / 14 NF, run three ways at identical group sizes. A reader who sees only
    # "DCM 5 vs NF 4" cannot locate the number quoted in the Results.
    ("Platform vs sampling — the matched-pair panel", [
        ("MAGNet → LBNL",
         f"same 24 individuals, 10 DCM vs 14 NF   → slope {_F1_PLATFORM:.3f}",
         "prep"),
        ("MAGNet → MAGNet",
         f"different individuals, 10 DCM vs 14 NF   → slope {_F1_PEOPLE:.3f}",
         "prep"),
    ]),
    ("Prep robustness — not replication", [
        ("LBNL", "DCM 5  vs  NF 4   (the 9 samples left after exclusion)", "prep"),
    ]),
]

FILL = {"disc": BLUE, "dcm": BLUE, "hcm": BLUE, "prep": ORANGE, "none": PALE}


def box(ax, x, y, w, h, text, sub=None, color=BLUE, alpha=.16, fs=6.6,
        bold=True, edge=True, txtcol=None):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.4",
                                linewidth=.8 if edge else 0,
                                edgecolor=color if edge else "none",
                                facecolor=color, alpha=alpha, zorder=1))
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.4",
                                linewidth=.8, edgecolor=color, facecolor="none",
                                zorder=3))
    cy = y + h / 2 + (1.1 if sub else 0)
    ax.text(x + w / 2, cy, text, ha="center", va="center", fontsize=fs,
            color=txtcol or INK, fontweight="bold" if bold else "normal", zorder=4)
    if sub:
        ax.text(x + w / 2, y + h / 2 - 1.5, sub, ha="center", va="center",
                fontsize=fs - 1.3, color=INK2, zorder=4)


def arrow(ax, x1, y1, x2, y2, color=MUTED, lw=.9, style="-|>", rad=0.0):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle=style,
                                 mutation_scale=7, linewidth=lw, color=color,
                                 connectionstyle=f"arc3,rad={rad}", zorder=2))


def panel_a(ax):
    ax.axis("off")
    ax.set_xlim(0, 100)
    # the tiers occupy 30-101; without this the panel reserves a third of its height
    # for nothing and panel B ends up marooned at the bottom of the page
    ax.set_ylim(29, 107)
    # the count comes from COHORTS, so removing a dataset can never leave a stale
    # number in the title -- which is exactly what happened when Heinig came out
    ax.set_title(f"A   Samples: {len(COHORTS)} public datasets, one identity screen, "
                 f"{N_ANALYSED} analysed", loc="left", fontsize=8.6,
                 fontweight="bold", pad=6)

    n = len(COHORTS)
    W, GAP = 11.0, 1.9
    x0 = (100 - (n * W + (n - 1) * GAP)) / 2 + 2.0
    xs = [x0 + i * (W + GAP) for i in range(n)]

    def tier(label, y):
        ax.text(1.6, y, label, fontsize=5.8, color=MUTED, ha="center", va="center",
                rotation=90)

    # tier 1 -- deposited
    tier("retrieved", 88.5)   # NOT "deposited" -- see the note in main()
    for (name, dep, ana, plat, kind), x in zip(COHORTS, xs):
        box(ax, x, 84, W, 9, name, f"n = {dep}", color=FILL[kind],
            alpha=.10 if kind == "none" else .14, fs=6.8)
        ax.text(x + W / 2, 81.4, plat, ha="center", va="center", fontsize=5.1,
                color=MUTED)
        arrow(ax, x + W / 2, 80.2, x + W / 2, 72.5)

    # the overlap, drawn between MAGNet and LBNL rather than left to a footnote
    xm, xl = xs[0] + W / 2, xs[2] + W / 2
    ax.add_patch(FancyArrowPatch((xm, 93.2), (xl, 93.2), arrowstyle="<|-|>",
                                 mutation_scale=7, linewidth=1.5, color=ORANGE,
                                 connectionstyle="arc3,rad=-0.34", zorder=5))
    ax.text((xm + xl) / 2, 104.2, f"{N_OVERLAP} of LBNL's 33 samples are the SAME "
                                 f"INDIVIDUALS as MAGNet samples",
            ha="center", va="center", fontsize=6.4, color=ORANGE, fontweight="bold")

    # tier 2 -- the screen. The exclusions live INSIDE this band: as separate callout
    # boxes they sat on top of the columns' arrows and the figure read as a tangle.
    bw = n * W + (n - 1) * GAP
    box(ax, x0, 62.5, bw, 10, "", color=INK2, alpha=.055, fs=6.4, bold=False)
    shown = sum(c[1] for c in COHORTS)
    # this line runs the full width of the band; at 6.4 pt the last word crossed the
    # right border, which reads as a rendering fault rather than a tight fit
    ax.text(x0 + bw / 2, 69.4,
            f"Genotype-identity screen:  all {GENOTYPED} samples ({shown} shown),  "
            f"{SNVS:,} SNVs,  KING kinship + IBS0,  mutual-best-partner rule",
            ha="center", va="center", fontsize=6.0, color=INK)
    ax.text(x0 + bw / 2, 65.4,
            f"{N_EXCLUDED} samples excluded:   {N_OVERLAP} MAGNet↔LBNL same "
            f"individuals   +   {N_REPLICATE} within-cohort technical replicates",
            ha="center", va="center", fontsize=6.2, color=ORANGE, fontweight="bold")

    for x in xs:
        arrow(ax, x + W / 2, 62.1, x + W / 2, 53.5)

    # tier 3 -- analysed
    tier("analysed", 48.5)
    for (name, dep, ana, plat, kind), x in zip(COHORTS, xs):
        sub = f"n = {ana}" if ana else "not analysed"
        box(ax, x, 44, W, 9, name, sub, color=FILL[kind],
            alpha=.10 if kind == "none" else .20, fs=6.8,
            txtcol=MUTED if kind == "none" else INK)
        lost = dep - ana
        if lost and kind != "none":
            ax.text(x + W / 2, 41.4, f"−{lost}", ha="center", va="center",
                    fontsize=5.8, color=ORANGE, fontweight="bold")
        elif kind == "none":
            # NOT an exclusion -- Heinig's samples are real, they have no controls
            ax.text(x + W / 2, 41.4, "no controls", ha="center", va="center",
                    fontsize=5.2, color=MUTED)

    ax.text(50, 35.5, f"{N_ANALYSED} samples analysed across the six cohorts with "
                      f"controls", ha="center", va="center", fontsize=7.4,
            color=INK, fontweight="bold")
    ax.text(50, 31.8, "the screen also covered public datasets retrieved but not "
                      "analysed here; the Methods name them and say why",
            ha="center", va="center", fontsize=5.6, color=MUTED)


def panel_b(ax):
    ax.axis("off")
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    ax.set_title("B   Contrasts, with post-exclusion group sizes", loc="left",
                 fontsize=8.6, fontweight="bold", pad=6)
    # six groups and ten rows: the earlier spacing ran the last row off the axes.
    HEAD, ROW, GAP = 4.0, 5.5, 1.2
    y = 94
    for head, items in CONTRASTS:
        quiet = head.startswith("Prep")
        ax.plot([1.2, 1.2], [y - HEAD - ROW * len(items) + 2.6, y + 2.4],
                color=ORANGE if quiet else BLUE, linewidth=1.9,
                alpha=.55 if quiet else .85, solid_capstyle="round")
        ax.text(3.6, y + 1.0, head, fontsize=6.6, color=MUTED if quiet else INK2,
                va="center", style="italic" if quiet else "normal",
                fontweight="normal" if quiet else "bold")
        y -= HEAD
        for name, desc, kind in items:
            box(ax, 3.6, y - 2.3, 19.0, 4.9, name, color=FILL[kind],
                alpha=.10 if quiet else .18, fs=6.4,
                txtcol=MUTED if quiet else INK)
            ax.text(25.0, y + 0.15, desc, fontsize=6.5, va="center",
                    color=MUTED if quiet else INK,
                    style="italic" if quiet else "normal")
            y -= ROW
        y -= GAP


def verify(path="Table1.tsv"):
    """The figure and the table must not drift. Check the counts that appear in both."""
    if not (os.path.exists(path) and os.path.getsize(path)):
        print(f"  [skip] {path} not found — cannot cross-check against Table 1")
        return 0
    T = pd.read_csv(path, sep="\t")
    col = next((c for c in T.columns if c.lower().startswith("cohort")), None)
    if col is None:
        print(f"  [skip] {path} has no Cohort column")
        return 0
    bad = 0
    fig_names = {c[0] for c in COHORTS}
    extra = [n for n in T[col].astype(str)
             if not any(n.startswith(f) for f in fig_names)]
    if extra:
        print(f"  [note] in {path} but deliberately not drawn in Figure 1: "
              + ", ".join(e.split(" (")[0] for e in extra))
    for name, dep, ana, plat, kind in COHORTS:
        row = T[T[col].astype(str).str.startswith(name)]
        if not len(row):
            print(f"  *** {name} is in the figure but not in {path}")
            bad += 1
            continue
        r = row.iloc[0]
        got = sum(int(x) for x in re.findall(r"(\d+)", str(r.get("Cases (n)", ""))))
        got += sum(int(x) for x in re.findall(r"(\d+)",
                                              str(r.get("Controls (n)", ""))))
        # Table 1 lists a not-analysed cohort at its DEPOSITED size (Heinig's 128 DCM
        # are real samples, they just have no controls), so compare against that.
        want = dep if kind == "none" else ana
        if got != want:
            print(f"  *** {name}: figure says {want}, {path} adds to {got}")
            bad += 1

        # --- read layout and length, the field this check used to ignore -------
        # Counts agreed between the figure and the table for weeks while MAGNet's
        # read length was wrong in both: 50 bp declared against 125 bp sequenced.
        # A cross-check that compares only counts cannot catch that, so it now
        # compares the platform string too -- normalised, because the table writes
        # "paired, 125 bp" and the figure "paired 125 bp" and that difference is
        # punctuation, not disagreement.
        lcol = next((c for c in T.columns
                     if any(k in c.lower() for k in
                            ("layout", "read", "platform", "sequenc"))), None)
        if lcol is None:
            continue
        def norm(s):
            s = str(s).lower()
            lay = ("paired" if "pair" in s else "single" if "singl" in s else "?")
            return lay, set(re.findall(r"\d+", s))
        f_lay, f_len = norm(plat)
        t_lay, t_len = norm(r.get(lcol, ""))
        if not t_len:
            continue                      # nothing to compare against
        if (f_lay, f_len) != (t_lay, t_len):
            print(f"  *** {name}: figure says '{plat}', {path} says "
                  f"'{r.get(lcol)}' (layout {f_lay} vs {t_lay}, "
                  f"lengths {sorted(f_len)} vs {sorted(t_len)})")
            bad += 1
    print(f"  [{'ok ' if not bad else '***'}] cross-check against {path}: "
          f"{bad} disagreement(s)")
    return bad


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--table", default="Table1.tsv")
    ap.add_argument("--out", default="F1_design")
    a = ap.parse_args()

    tot = sum(c[2] for c in COHORTS)
    if tot != N_ANALYSED:
        print(f"  *** the cohort rows add to {tot}, not {N_ANALYSED}. Fix the "
              f"constants before using this figure.")

    plt.rcParams.update({"font.family": "sans-serif", "axes.linewidth": .8,
                         "figure.facecolor": SURF, "savefig.facecolor": SURF})
    fig = plt.figure(figsize=(7.8, 8.2))
    gs = gridspec.GridSpec(2, 1, height_ratios=[1.0, 0.95], hspace=0.13,
                           left=0.035, right=0.975, top=0.955, bottom=0.02)
    panel_a(fig.add_subplot(gs[0]))
    panel_b(fig.add_subplot(gs[1]))
    for ext in ("pdf", "png"):
        fig.savefig(f"{a.out}.{ext}", dpi=400)
    print(f"[write] {a.out}.pdf and {a.out}.png")

    print("\n" + "=" * 78)
    print("CONSISTENCY")
    print("=" * 78)
    verify(a.table)
    shown = sum(c[1] for c in COHORTS)
    print(f"\n  The drawn cohorts hold {shown} retrieved samples; the screen covered "
          f"{GENOTYPED}.\n  The {GENOTYPED - shown} difference is the retrieved-but-"
          f"unused datasets, and the Methods\n  must name them (Heinig PRJEB8360, "
          f"Amsterdam, Naples) and state why each was set aside.")
    print("\n  Figure 1 draws only the datasets that were analysed. Heinig and any "
          "other\n  retrieved-but-unused dataset belong in one Methods sentence, not "
          "in a column\n  that leads nowhere.")
    print("\n  CHECK: the top tier is labelled 'retrieved', not 'deposited'. MAGNet's")
    print("  360 is DCM + NF + HCM; PRJNA595151 also holds PPCM samples, and the other")
    print("  BioProjects may hold runs not used here. If you want 'deposited', put the")
    print("  true BioProject totals in COHORTS first.")
    print("\n  Figure 1 and Table 1 carry the same numbers by design. If you edit one,")
    print("  rerun the other and this check.")
    print("\n  NOW LOOK AT THE PNG — boxes and labels collide easily at this density.")


if __name__ == "__main__":
    main()
