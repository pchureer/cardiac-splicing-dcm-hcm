#!/usr/bin/env python3
"""
47_graphical_abstract.py -- GRAPHICAL ABSTRACT (required by Cardiovascular Research).

Three blocks, left to right, and one sentence underneath:

  1  INPUTS      510 hearts, six cohorts, the genotype screen that removed 27 samples
  2  THE RESULT  six contrasts in two cohorts -- each disease against donor myocardium
                 is large, each disease against another disease is ~0. This is the
                 double positive control, and it is the whole argument in one picture.
  3  MECHANISM   ten candidate regulators against the GAPDH floor, all inside it

NO NEW ANALYSIS. Every number is either frozen from the project docs or read back out
of master_numbers_flagged.tsv, and main() prints a comparison of the two.

DESIGNED FOR THE THUMBNAIL, NOT FOR THE PAGE
--------------------------------------------
A graphical abstract is read at table-of-contents width -- two inches, a quarter of
what you view it at. The first version of this script was a normal figure shrunk: six
small cohort chips, a dot plot, four lines of qualifying grey text. At thumbnail every
one of those turned to mush. This version is built the other way round:

  * the numbers ARE the graphic. 926 / 1,443 / 0 and 476 / 407 / 2 are set large
    enough to read small, and the two nulls are set largest of all, because they are
    the result.
  * one idea per block, and nothing in a block that is not that idea. The cohort
    chips, the aligner names and the FDR threshold came out -- they live in Table 1
    and the Methods, where a reader who wants them will look.
  * three type sizes, far apart. Anything that needed a fourth, smaller size was cut
    instead of shrunk.

The script writes a thumbnail alongside the full version. Judge it there. If a label
is unreadable in GA_thumbnail.png it does not belong in the panel.

WHAT IT MUST NOT DO
-------------------
Lead with a zero. The DCM-vs-HCM null is drawn next to the 1,443 events the same 28
hearts produce against non-failing myocardium, and the Colorado pair next to its 476
and 407. A bare zero reads as a failed experiment; a zero beside a positive control
from the same hearts reads as a result. Also: no Colorado MXE anywhere (it does not
reproduce -- RESULT_colorado_mxe_excluded.md), and no cross-cohort effect sizes
(not estimable -- RESULT_magnitude_not_estimable.md).

Usage:
  export PATH=/data/apps/anaconda3/envs/mamba/envs/rnaseq_tools/bin:$PATH
  cd /data/project/chureeratp/RNA/HEART
  python 47_graphical_abstract.py
  python 47_graphical_abstract.py --no-verify
"""
import argparse
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
PALE = "#c9c8c3"

# ---------------------------------------------------------------------------
# FROZEN. Identical to Table 1 and Figure 1 -- 43_figure1_design.py COHORTS.
# ---------------------------------------------------------------------------
COHORTS = [("MAGNet", 360), ("Colorado", 62), ("LBNL", 9),
           ("Brigham", 23), ("Fuwai", 37), ("Stanford", 19)]
N_ANALYSED = 510
GENOTYPED = 670
N_EXCLUDED = 27
N_OVERLAP = 24

# ---------------------------------------------------------------------------
# THE SIX CONTRASTS. All exon-level SE, FDR<0.05 and |dPSI|>=0.10 -- the locked rule.
#   label, significant, tested, cohort block, is_the_null
# Sources: PAPER_ONE_PAGE.md (30 Sep 2026), RESULT_hcm_specificity.md,
# RESULT_colorado_icm_replication.md (verified by 46_colorado_icm_lookup.py against
# both the frozen master file and a recount from the _dedup tables).
#
# THE HCM NUMBER WAS WRONG AND IS NOW FIXED (3 Oct 2026).
# It read 1,443, which is the FULL-COORDINATE, FDR-ONLY count. Every other number in
# this figure is exon-level at FDR<0.05 and |dPSI|>=0.10, so the panel was comparing a
# thresholded zero against an unthresholded full-coordinate count -- in the two rows
# that carry the paper's primary claim. The exon-level locked-rule value is 468
# (54_fdr_only_counts.py, SE_HCM_v2.all_events.tsv, 48,024 tested).
#
# Every count below is now ONE convention: exon-level, FDR<0.05, |dPSI|>=0.10.
# Do not mix in an FDR-only or full-coordinate figure without relabelling the panel.
# ---------------------------------------------------------------------------
CONTRASTS = [
    # MAGNet: both arms of the third row are explanted hearts, so that row is immune
    # to the donor-versus-explant confound. That is why it carries the claim.
    ("DCM  vs  non-failing",   926, 62971,  "MAGNet",   False),
    ("HCM  vs  non-failing",   468, 48024,  "MAGNet",   False),
    ("DCM  vs  HCM",             0, 45380,  "MAGNet",   True),
    ("DCM  vs  non-failing",   476, 19953,  "Colorado", False),
    ("ICM  vs  non-failing",   407, 17019,  "Colorado", False),
    ("DCM  vs  ICM",             2, 19998,  "Colorado", True),
]
# group composition, for the one line under each cohort's name
BLOCK_SUB = {"MAGNet": "166 DCM · 28 HCM · 166 non-failing",
             "Colorado": "36 DCM · 13 ICM · 13 non-failing"}

# Inverted scan, SE deltas -- 40_figure5_mechanism.py FROZEN_SCAN, which is
# RESULT_rbm20_mechanism_null.md. delta = median |t| target - median |t| background.
SCAN = {"PTBP1": +0.08, "QKI": +0.07, "HNRNPF": +0.07, "ZRANB2": +0.06,
        "HNRNPA0": +0.05, "MBNL1": +0.05, "GAPDH": +0.05, "SRSF1": +0.04,
        "RBM20": +0.03, "CELF1": +0.01, "RBM5": -0.01}
CONTROL = "GAPDH"
LABELLED = ("GAPDH", "RBM20", "QKI")   # the only three a thumbnail can carry
RBM20_DE = (+0.34, 0.0025)             # SD units, DCM vs NF, composition-adjusted

STRAP = ("Shared across cardiomyopathy etiologies  ·  not explained by RBM20 "
         "or nine other candidate regulators")

# three sizes, far apart. Nothing in this figure is set between them.
BIG, MID, SMALL = 19.0, 8.4, 6.8


# ---------------------------------------------------------------------------
# block 1 -- inputs. A hero number, not a diagram.
# ---------------------------------------------------------------------------

def block_inputs(ax):
    """The six cohort chips of the first version are gone. At thumbnail width they
    were six grey smudges, and the cohort composition belongs to Table 1 anyway. What
    survives is the one number a reader needs (510) and the one finding that makes
    this reanalysis different from the others (the identity screen)."""
    ax.axis("off")
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)

    ax.text(0, 93, "510", fontsize=BIG * 1.5, fontweight="bold", color=INK,
            va="top", ha="left")
    ax.text(0, 68, "human hearts", fontsize=MID * 1.15, color=INK, va="top",
            ha="left")
    ax.text(0, 57, "six public RNA-seq cohorts", fontsize=MID * 0.95, color=INK2,
            va="top", ha="left")
    ax.text(0, 47, "MAGNet · Colorado · Brigham\nFuwai · Stanford · LBNL",
            fontsize=SMALL, color=MUTED, va="top", ha="left", linespacing=1.5)

    for face in (ORANGE, "none"):
        ax.add_patch(FancyBboxPatch((0, 2), 100, 22, boxstyle="round,pad=0.8",
                                    linewidth=1.0, edgecolor=ORANGE,
                                    facecolor=face,
                                    alpha=0.11 if face != "none" else 1.0,
                                    zorder=1 if face != "none" else 2))
    boxed = [
        ax.text(50, 17.5, f"{N_EXCLUDED} samples excluded", ha="center",
                va="center", fontsize=MID * 1.08, fontweight="bold", color=ORANGE,
                zorder=3),
        ax.text(50, 7.5, "same individuals found\nin two bioprojects", ha="center",
                va="center", fontsize=SMALL, color=INK2, linespacing=1.5, zorder=3),
    ]
    # the callout text must sit INSIDE its tinted box. main() enforces it.
    return boxed, 0.0, 100.0


# ---------------------------------------------------------------------------
# block 2 -- the double positive control. The numbers are the graphic.
# ---------------------------------------------------------------------------

def block_result(ax, contrasts):
    """Two cohort blocks, three rows each, set as numerals with a bar behind them.

    EACH COHORT IS SCALED TO ITS OWN MAXIMUM. A shared axis would put MAGNet's 1,443
    beside Colorado's 476 and invite the reader to conclude that MAGNet carries the
    larger signature. It does not mean that: the two cohorts test different event sets
    at different depths under different covariate models, and
    RESULT_magnitude_not_estimable.md is the whole argument for why cross-cohort counts
    are not comparable. Each block gets its own tinted panel so that the separation is
    visible rather than only stated, and the exact counts are printed regardless."""
    ax.axis("off")
    rows = []
    for block in ("MAGNet", "Colorado"):
        rows += [c for c in contrasts if c[3] == block]

    X0, X1 = 0.0, 1.0            # bar track
    LAB = -0.06                  # right edge of the row labels
    PANEL_L, PANEL_R = -1.04, 1.46
    ax.set_xlim(PANEL_L - 0.02, PANEL_R + 0.02)

    # y layout, top down: header, three rows, gap, header, three rows
    ys, heads, y = {}, {}, 0.0
    for block in ("MAGNet", "Colorado"):
        heads[block] = y
        y -= 1.05
        ys[block] = []
        for _ in [c for c in rows if c[3] == block]:
            ys[block].append(y)
            y -= 1.0
        y -= 0.55
    bottom = y

    bmax = {b: max([c[1] for c in rows if c[3] == b] or [1]) or 1
            for b in ("MAGNet", "Colorado")}

    boxed = []
    for block in ("MAGNet", "Colorado"):
        bl = [c for c in rows if c[3] == block]
        top = heads[block] + 0.42
        low = ys[block][-1] - 0.60
        ax.add_patch(FancyBboxPatch((PANEL_L, low), PANEL_R - PANEL_L, top - low,
                                    boxstyle="round,pad=0.015", linewidth=0,
                                    facecolor=INK2, alpha=0.045, zorder=0))
        boxed.append(ax.text(PANEL_L + 0.07, heads[block], block, fontsize=MID,
                             fontweight="bold", color=INK, va="center", ha="left",
                             zorder=3))
        boxed.append(ax.text(PANEL_L + 0.07, heads[block] - 0.46, BLOCK_SUB[block],
                             fontsize=SMALL * 0.85, color=MUTED, va="center",
                             ha="left", zorder=3))

        for (label, sig, tested, _b, isnull), yy in zip(bl, ys[block]):
            if isnull:
                # no bar at all. A bar of length ~0 reads as a rendering fault; the
                # numeral, set larger than anything else in the figure, reads as the
                # result it is.
                boxed.append(ax.text(LAB, yy, label, fontsize=MID,
                                     fontweight="bold", color=INK, va="center",
                                     ha="right", zorder=3))
                boxed.append(ax.text(X0 + 0.04, yy, f"{sig:,}", fontsize=BIG,
                                     fontweight="bold", color=ORANGE, va="center",
                                     ha="left", zorder=3))
            else:
                w = X0 + (X1 - X0) * sig / bmax[_b]
                ax.barh(yy, w - X0, left=X0, height=0.52, color=BLUE, alpha=0.82,
                        edgecolor="none", zorder=2)
                boxed.append(ax.text(LAB, yy, label, fontsize=SMALL * 1.0,
                                     color=INK2, va="center", ha="right", zorder=3))
                boxed.append(ax.text(w + 0.03, yy, f"{sig:,}", fontsize=MID * 1.05,
                                     fontweight="bold", color=BLUE, va="center",
                                     ha="left", zorder=3))

    ax.set_ylim(bottom - 0.55, 1.35)
    ax.text(PANEL_L, 1.22, "significant splicing events", fontsize=MID * 1.12,
            fontweight="bold", color=INK, va="center", ha="left")
    ax.text(PANEL_L, 0.78, "each disease vs donor myocardium, then vs each other",
            fontsize=SMALL * 0.95, color=INK2, va="center", ha="left")
    ax.text(PANEL_L, bottom - 0.20,
            "bars scaled within each cohort · counts are exact",
            fontsize=SMALL * 0.82, color=MUTED, va="center", ha="left")
    return boxed, PANEL_L, PANEL_R


# ---------------------------------------------------------------------------
# block 3 -- mechanism. A verdict with a strip under it, not a dot plot.
# ---------------------------------------------------------------------------

def block_mechanism(ax):
    """Eleven points on one line, the control drawn as the line they fail to cross.

    The first version stacked the eleven vertically with an x-axis and tick labels.
    That is the right figure for the Results (it is Figure 5C) and the wrong one here:
    at thumbnail it was eleven pale specks. One row, larger marks, three labels."""
    ax.axis("off")
    vals = SCAN
    floor = vals[CONTROL]
    lo, hi = min(vals.values()), max(vals.values())
    pad = (hi - lo) * 0.22
    ax.set_xlim(lo - pad * 1.6, hi + pad * 1.6)
    ax.set_ylim(-4.0, 2.4)

    ax.text(lo - pad * 1.6, 2.25, "10 candidate regulators", fontsize=MID * 1.1,
            fontweight="bold", color=INK, va="center", ha="left")
    ax.text(lo - pad * 1.6, 1.84, "targets enriched in the signature?",
            fontsize=SMALL * 0.95, color=INK2, va="center", ha="left")

    # the control is a line, because the claim is that nothing crosses it
    ax.plot([floor, floor], [-0.55, 1.15], color=ORANGE, linewidth=2.0, zorder=2,
            solid_capstyle="round")
    ax.text(floor, 1.30, CONTROL, fontsize=SMALL * 1.05, fontweight="bold",
            color=ORANGE, ha="center", va="bottom", zorder=4)
    ax.plot([lo - pad, hi + pad], [0, 0], color=GRID, linewidth=1.0, zorder=1)
    # without this the strip has no scale and the reader cannot tell which way is
    # "more". One word at one end is enough; a full axis is not.
    ax.text(hi + pad, -0.95, "more enriched \u2192", fontsize=SMALL * 0.8,
            color=MUTED, ha="right", va="center")

    # ties get a small vertical offset so eleven points stay eleven points
    seen = {}
    for nm, v in sorted(vals.items(), key=lambda kv: kv[1]):
        k = round(v, 3)
        off = seen.get(k, 0)
        seen[k] = off + 1
        yy = 0 + off * 0.46
        hit = nm in LABELLED
        col = ORANGE if nm == CONTROL else (INK if nm == "RBM20" else BLUE)
        ax.scatter([v], [yy], s=92 if hit else 52, color=col,
                   alpha=1.0 if hit else 0.55, zorder=3, edgecolor=SURF,
                   linewidth=0.9 if hit else 0.6)
        if nm in ("RBM20", "QKI"):
            ax.text(v, yy - 0.50, nm, fontsize=SMALL * 1.0, fontweight="bold",
                    color=col, ha="center", va="top", zorder=4)

    ax.text((lo + hi) / 2, -1.95, "none beats the control",
            fontsize=MID * 1.15, fontweight="bold", color=INK,
            ha="center", va="center")
    ax.text((lo + hi) / 2, -3.20, f"RBM20 is HIGHER in DCM,\nnot lower  "
                                  f"({RBM20_DE[0]:+.2f} SD, P = {RBM20_DE[1]:.3f})",
            fontsize=SMALL * 1.0, fontweight="bold", color=ORANGE,
            ha="center", va="center", linespacing=1.6)


# ---------------------------------------------------------------------------
# keeping text inside the shapes behind it
# ---------------------------------------------------------------------------

def fit_inside(fig, ax, texts, lo_data, hi_data, pad_px=6.0, min_fs=4.6):
    """Shrink any label that sticks out of the panel drawn behind it.

    Text is measured in points and panels in data units, so whether a label fits is
    decided by the font the machine happens to have -- and DejaVu here is not
    necessarily what the cluster renders. Two labels overflowed their tinted boxes in
    an earlier version for exactly that reason. Rather than tune font sizes by eye
    against one machine's fonts, measure the drawn extent and shrink until it fits.
    Returns what it had to change, so a silent shrink cannot hide a layout that is
    actually too tight."""
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    lo = ax.transData.transform((lo_data, 0))[0]
    hi = ax.transData.transform((hi_data, 0))[0]
    changed, stuck = [], []
    for t in texts:
        # Only one side of a label can be fixed by shrinking it: the anchor does not
        # move. A left-aligned label's left edge stays where it was put however small
        # the type gets, so testing it is how an earlier version drove every heading
        # down to the floor trying to fix a margin that was never the text's fault.
        ha = t.get_horizontalalignment()
        free_lo = ha in ("right", "center")
        free_hi = ha in ("left", "center")
        anchor = ax.transData.transform((t.get_position()[0], 0))[0]
        name = t.get_text().split("\n")[0]
        if not (lo - 0.5 <= anchor <= hi + 0.5):
            stuck.append((name, "anchored outside the panel — fix the geometry"))
            continue
        fs0 = t.get_fontsize()
        for _ in range(120):
            bb = t.get_window_extent(renderer=r)
            if not ((free_lo and bb.x0 < lo + pad_px)
                    or (free_hi and bb.x1 > hi - pad_px)):
                break
            fs = t.get_fontsize()
            if fs <= min_fs:
                stuck.append((name, f"still overflows at the {min_fs:.1f} pt floor"))
                break
            t.set_fontsize(fs - 0.2)
        if t.get_fontsize() < fs0 - 1e-9:
            changed.append((name, fs0, t.get_fontsize()))
    return changed, stuck


# ---------------------------------------------------------------------------
# verification -- the frozen bars against the master file
# ---------------------------------------------------------------------------

NKEEP = {"SE": 3, "MXE": 5, "A3SS": 5, "A5SS": 5, "RI": 3}


def norm(s):
    return "".join(ch for ch in str(s).lower() if ch.isalnum())


def verify(path, contrasts):
    """Look every drawn number up in master_numbers_flagged.tsv. SE rows only -- the
    project's headline replication and specificity numbers are skipped-exon, and the
    all-type totals are a different statistic that must not be substituted."""
    print("\n" + "=" * 78)
    print("VERIFICATION — the drawn numbers against the frozen master file")
    print("=" * 78)
    if not (os.path.exists(path) and os.path.getsize(path)):
        print(f"  [skip] {path} not found — every number is unverified frozen text")
        return 1
    M = pd.read_csv(path, sep="\t")
    need = [c for c in ("cohort", "contrast") if c not in M.columns]
    if need:
        print(f"  [skip] {path} has no {' / '.join(need)} column — cannot match rows")
        return 1
    sc = "sig_exon" if "sig_exon" in M.columns else "significant"
    tc = "tested_exon" if "tested_exon" in M.columns else "tested"
    if sc not in M.columns:
        print(f"  [skip] {path} has neither sig_exon nor significant")
        return 1
    if "event_type" in M.columns:
        M = M[M.event_type.astype(str).str.upper() == "SE"]
        print(f"  matching SE rows only ({len(M)} in the file)")
    else:
        print("  *** no event_type column: cannot confirm these rows are SE. The "
              "headline\n      numbers are skipped-exon; if this file mixes types the "
              "comparison below\n      is meaningless. Check 26_flag_master_numbers.py.")

    bad = 0
    for label, sig, tested, block, isnull in contrasts:
        # 'DCM  vs  non-failing' -> dcmvsnonfailing / dcmvsnf, both spellings occur
        want = norm(label).replace("nonfailing", "nf")
        rows = M[M.cohort.astype(str).str.contains(block, case=False, na=False)]
        hit = None
        for _, r in rows.iterrows():
            if norm(r["contrast"]).replace("nonfailing", "nf") == want:
                hit = r
                break
        if hit is None:
            print(f"  [    ] {block:<9}{label:<24}no matching row — stays frozen")
            continue
        got_s = int(hit[sc])
        got_t = int(hit[tc]) if tc in M.columns and pd.notna(hit.get(tc)) else None
        ok = got_s == sig and (tested is None or got_t is None or got_t == tested)
        bad += 0 if ok else 1
        extra = ""
        if tested is not None and got_t is not None and got_t != tested:
            extra = f"   tested: drawn {tested:,}, file {got_t:,}"
        print(f"  [{'ok ' if ok else '***'}] {block:<9}{label:<24}"
              f"drawn {sig:>6,}   file {got_s:>6,}{extra}")

    if bad:
        print(f"\n  *** {bad} number(s) disagree with {path}. The master file is the "
              f"source of\n      truth — fix CONTRASTS here, not the file, and rerun "
              f"Figures 2 and 3 if they\n      share the number.")
    else:
        print("\n  [ok ] every number that could be matched agrees with the master "
              "file.")

    print("\n  CONVENTION: every count in this figure is exon-level, FDR<0.05 and")
    print("  |dPSI|>=0.10. The HCM vs non-failing cell was 1,443 until 3 Oct 2026,")
    print("  which was the full-coordinate FDR-only count — a different rule from")
    print("  every other number in the panel. It is 468. If any cell here disagrees")
    print("  with 54_fdr_only_counts.py's exon-level +|dPSI|>=0.10 column, the figure")
    print("  is wrong, not the script.")
    return bad


def check_constants():
    tot = sum(n for _, n in COHORTS)
    bad = 0
    if tot != N_ANALYSED:
        print(f"  *** the cohort rows add to {tot}, not {N_ANALYSED}. Fix COHORTS "
              f"against Table1.tsv before using this figure.")
        bad += 1
    if CONTROL not in SCAN:
        print(f"  *** {CONTROL} is not in SCAN — the control line has nothing to sit "
              f"on.")
        bad += 1
    over = {k: v - SCAN[CONTROL] for k, v in SCAN.items()
            if k != CONTROL and v - SCAN[CONTROL] > 0.03}
    if over:
        print("  *** SCAN has drifted: a candidate now sits more than 0.03 "
              "above the control: "
              + ", ".join(f"{k} {v:+.2f}" for k, v in over.items())
              + "\n      0.03 is a tripwire on these point estimates, NOT the "
                "paper's claim. The claim is that no candidate is "
                "\n      DISTINGUISHABLE from the control -- largest margin "
                "+0.033 (95% CI -0.025 to +0.082), see 58_signature_rbp_margins.py."
              + "\n      Block 3 draws 'none beats the control'. Re-check it "
                "against 58 before use.")
        bad += 1
    if not bad:
        print(f"  [ok ] constants self-consistent: {tot} samples, {len(SCAN)} "
              f"candidates, control = {CONTROL} at {SCAN[CONTROL]:+.2f}")
    return bad


# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--master", default="master_numbers_flagged.tsv")
    ap.add_argument("--out", default="GA_graphical_abstract")
    ap.add_argument("--thumb-width", type=float, default=2.0,
                    help="inches; the width a table of contents will show it at")
    ap.add_argument("--no-verify", dest="ver", action="store_false")
    a = ap.parse_args()

    print("=" * 78)
    print("CONSTANTS")
    print("=" * 78)
    check_constants()

    plt.rcParams.update({"font.family": "sans-serif", "axes.linewidth": 0.8,
                         "figure.facecolor": SURF, "savefig.facecolor": SURF})

    fig = plt.figure(figsize=(8.0, 4.0))
    gs = gridspec.GridSpec(1, 3, width_ratios=[0.80, 1.42, 0.98], wspace=0.30,
                           left=0.035, right=0.985, top=0.95, bottom=0.135)
    axes = [fig.add_subplot(gs[i]) for i in range(3)]
    for ax in axes:
        ax.set_facecolor(SURF)

    fits = [block_inputs(axes[0]), block_result(axes[1], CONTRASTS)]
    block_mechanism(axes[2])

    shrunk, stuck = [], []
    for ax, (texts, lo, hi) in zip(axes, fits):
        c, s = fit_inside(fig, ax, texts, lo, hi)
        shrunk += c
        stuck += s
    if shrunk:
        print("\n  [fit] these labels were wider than the panel behind them and have "
              "been\n        shrunk to fit. A drop of more than ~1 pt means the "
              "layout is too tight —\n        shorten the wording rather than letting "
              "it shrink again:")
        for txt, a0, b0 in shrunk:
            print(f"        {txt[:42]:<44}{a0:.1f} -> {b0:.1f} pt")
    if stuck:
        print("\n  *** these labels still do not fit their panel. Shrinking cannot "
              "fix them —\n      change the wording or the geometry:")
        for txt, why in stuck:
            print(f"      {txt[:42]:<44}{why}")
    if not shrunk and not stuck:
        print("\n  [fit] every label sits inside its panel at the intended size.")

    # flow between the blocks, so the panel reads left to right without a caption
    for x in (0.283, 0.722):
        fig.text(x, 0.545, "❯", ha="center", va="center", fontsize=13,
                 color=PALE)

    fig.text(0.5, 0.045, STRAP, ha="center", va="center",
             fontsize=MID * 1.0, color=INK, fontweight="bold")

    for ext in ("pdf", "png", "tif"):
        try:
            # LZW, or the 600-dpi TIFF lands at ~40 MB and most submission systems
            # cap uploads well below that
            kw = {"pil_kwargs": {"compression": "tiff_lzw"}} if ext == "tif" else {}
            fig.savefig(f"{a.out}.{ext}", dpi=600 if ext == "tif" else 400, **kw)
        except Exception as e:                       # TIFF needs Pillow
            print(f"  [skip] {ext}: {e}")
    print(f"\n[write] {a.out}.pdf, {a.out}.png (400 dpi), {a.out}.tif (600 dpi)")

    # the legibility check, because a graphical abstract is read small
    thumb = f"{a.out}_thumbnail.png"
    fig.savefig(thumb, dpi=a.thumb_width / fig.get_size_inches()[0] * 150)
    print(f"[write] {thumb} — this is roughly what a table of contents shows at "
          f"{a.thumb_width:g} in wide.")

    if a.ver:
        verify(a.master, CONTRASTS)

    print("\n" + "=" * 78)
    print("BEFORE YOU SUBMIT IT")
    print("=" * 78)
    print(f"  1. OPEN {thumb} FIRST, not the full-size PNG. Anything you cannot read "
          f"there\n     does not belong in the panel — cut it rather than shrink "
          f"something else.")
    print("  2. Confirm the MAGNet HCM vs non-failing event type (see the note "
          "above).")
    print("  3. Cardiovascular Research wants the graphical abstract as a SEPARATE "
          "file,\n     not embedded in the manuscript, and it does not count against "
          "the 7 figures.")
    print("  4. Check the strap line still matches the panel if any number moves.")
    print("  5. No Colorado MXE appears here, and no cross-cohort effect size. Keep "
          "it that\n     way — both are excluded results, not omissions.")


if __name__ == "__main__":
    main()
