#!/usr/bin/env python3
"""
38_figure4_replication.py -- FIGURE 4: the signature replicates in direction, and the
magnitude slopes are bounded by measurement rather than by biology.

  A  Direction concordance forest. Per cohort with Wilson intervals, the pooled
     estimate, and the 50% null. LBNL is shown BELOW a rule, greyed, excluded from
     the pool -- 24 of its 33 samples are MAGNet hearts.
  B  Transfer slopes against the Phase 1.5 ceilings. Each cohort's Deming slope is
     plotted with two reference bands measured on the 24 same-individual pairs:
     0.99 (changing people only) and 0.58 (changing platform only). The cohorts sit
     at or below the platform ceiling, which is the point -- the slopes are a
     measurement property, not weak biology.

WHY PANEL B IS BUILT THIS WAY
-----------------------------
Plotting per-cohort slopes side by side is exactly what script 22 showed you cannot
do as a magnitude comparison, and it is why fig_split_control panel B has to be
recast. It becomes legitimate here ONLY because the two horizontal reference lines
convert it from "cohort X attenuates more than cohort Y" (uninterpretable) into
"every cohort falls at or below the technical ceiling" (the claim). The caption must
carry that framing or the panel reverts to the anti-pattern.

Sources, all frozen:
  panel A  direction_concordance_dedup.tsv   (33_rebuild_summaries_dedup.py)
  panel B  recomputed from the same _dedup tables with the same Deming estimator
  ceilings RESULT_phase15_platform_vs_sampling.md  (0.584 and 0.988)

Usage:
  python 38_figure4_replication.py
  python 38_figure4_replication.py --event-type SE --out F4_replication
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import gridspec

BLUE = "#2a78d6"
ORANGE = "#eb6834"
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#8a8a85"
GRID = "#e4e4e0"
SURF = "#fcfcfb"

NKEEP = {"SE": 3, "MXE": 5, "A3SS": 5, "A5SS": 5, "RI": 3}
# measured on the 24 MAGNet<->LBNL same-individual pairs, common event set
# 162: PLOTTED CONSTANTS. Figure 5's panel D was drawn from literals exactly
# like these and carried 50 bp values into a 125 bp figure while every gate
# passed. At 125 bp the matched panel gives 0.690 and 0.968
# (RESULT_s4_rl125.md).  F4_CEILINGS="platform/people", e.g. "0.690/0.968".
_C4 = os.environ.get("F4_CEILINGS", "").strip()
if _C4:
    CEIL_PLATFORM, CEIL_PEOPLE = (float(_x) for _x in _C4.split("/"))
else:
    CEIL_PEOPLE = 0.988      # different people, same platform
    CEIL_PLATFORM = 0.584    # same people, different platform
# 162: only the MAGNet discovery arm re-bases. Every REPL entry below is a
# cohort that was NOT re-quantified and must stay unsuffixed.
F4SUF = os.environ.get("F4_SUFFIX", "")
REPL = [("Colorado", "{ET}_Colorado_DCMvsNF"), ("Brigham", "Brigham_{ET}_HCM"),
        ("Fuwai", "Fuwai_{ET}_HCM"), ("Stanford", "Stanford_{ET}_HCM"),
        ("LBNL", "{ET}_LBNL_DCMvsNF")]
NOT_INDEPENDENT = {"LBNL"}


def load(base, et):
    for p in (f"{base}_dedup.all_events.tsv", f"{base}.all_events.tsv"):
        if os.path.exists(p) and os.path.getsize(p):
            t = pd.read_csv(p, sep="\t")
            d = next((c for c in ("dPSI_adj", "dPSI") if c in t.columns), None)
            if d is None:
                cand = [c for c in t.columns if c.lower().startswith("dpsi")]
                d = next((c for c in cand if "adj" in c.lower()), cand[0] if cand else None)
            if d is None or "key" not in t.columns:
                return None
            if "se_adj" not in t.columns and "se" in t.columns:
                t = t.rename(columns={"se": "se_adj"})
            n = NKEEP[et]
            t["ekey"] = [":".join(str(k).split(":")[:n]) for k in t["key"]]
            t = t.loc[t.groupby("ekey")["FDR"].idxmin()]
            return t.rename(columns={d: "dpsi"})[["ekey", "dpsi", "se_adj", "FDR"]]
    return None


def wilson(k, n, z=1.96):
    if n == 0:
        return np.nan, np.nan
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (c - h) * 100, (c + h) * 100


def deming(x, y, sx, sy):
    ok = np.isfinite(x) & np.isfinite(y) & np.isfinite(sx) & np.isfinite(sy)
    x, y, sx, sy = x[ok], y[ok], sx[ok], sy[ok]
    if len(x) < 10:
        return np.nan
    lam = np.mean(sy ** 2) / np.mean(sx ** 2) if np.mean(sx ** 2) > 0 else 1.0
    dx, dy = x - x.mean(), y - y.mean()
    Sxx, Syy, Sxy = (dx ** 2).sum(), (dy ** 2).sum(), (dx * dy).sum()
    if Sxy == 0:
        return np.nan
    return float((Syy - lam * Sxx + np.sqrt((Syy - lam * Sxx) ** 2
                                            + 4 * lam * Sxy ** 2)) / (2 * Sxy))


def style(ax):
    ax.set_facecolor(SURF)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(GRID)
        ax.spines[s].set_linewidth(0.8)
    ax.tick_params(colors=INK2, labelsize=7, length=3, width=0.8)


def panel_a(ax, C):
    C = C.copy()
    C["indep"] = ~C.cohort.isin(NOT_INDEPENDENT)
    C = C.sort_values(["indep", "pct"], ascending=[True, True]).reset_index(drop=True)
    pooled = C[C.indep]
    tot_k, tot_n = int(pooled.same.sum()), int(pooled.n.sum())
    rows = list(C.itertuples())
    y = np.arange(len(rows) + 1)          # +1 for the pooled row at the top

    for i, r in enumerate(rows):
        lo, hi = wilson(r.same, r.n)
        col = BLUE if r.indep else MUTED
        ax.plot([lo, hi], [i, i], color=col, linewidth=1.6, solid_capstyle="round",
                zorder=2)
        ax.plot([r.pct], [i], "o", ms=5.5, color=col, zorder=3,
                markeredgecolor=SURF, markeredgewidth=1.1)
        ax.text(101, i, f"{r.same:,}/{r.n:,}", fontsize=6.4,
                color=INK if r.indep else MUTED, va="center", ha="left")

    yp = len(rows)
    plo, phi = wilson(tot_k, tot_n)
    ppct = 100 * tot_k / tot_n
    ax.plot([plo, phi], [yp, yp], color=ORANGE, linewidth=2.4,
            solid_capstyle="round", zorder=2)
    ax.plot([ppct], [yp], "D", ms=6, color=ORANGE, zorder=3,
            markeredgecolor=SURF, markeredgewidth=1.1)
    ax.text(101, yp, f"{tot_k:,}/{tot_n:,}", fontsize=6.6, color=INK,
            va="center", ha="left", fontweight="bold")

    ax.axvline(50, color=MUTED, linewidth=1, linestyle=(0, (4, 3)), zorder=1)
    ax.text(50, -0.85, "chance", fontsize=6.2, color=MUTED, ha="center", va="top")
    # separator: everything below it is not an independent cohort
    nsep = int((~C.indep).sum())
    if nsep:
        ax.axhline(nsep - 0.5, color=GRID, linewidth=0.9)

    labels = [(r.cohort if r.indep else f"{r.cohort} *") for r in rows] + ["pooled"]
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=7.5, color=INK)
    for t, r in zip(ax.get_yticklabels()[:len(rows)], rows):
        if not r.indep:
            t.set_color(MUTED)
    ax.get_yticklabels()[-1].set_fontweight("bold")
    ax.set_xlim(45, 100)
    ax.set_ylim(-1.4, len(rows) + 0.7)
    ax.set_xlabel("% of MAGNet-significant events with the same direction",
                  fontsize=7.5, color=INK2)
    ax.xaxis.grid(True, color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)
    style(ax)
    ax.set_title("A   Direction replicates", loc="left", fontsize=8.6,
                 fontweight="bold", pad=8)
    # Panel B reuses these exact y positions, so the two panels line up row for row
    # and a reader can scan across. Ordering each panel by its own values instead
    # would silently swap Brigham and Colorado between them.
    ypos = {r.cohort: i for i, r in enumerate(rows)}
    return ppct, tot_k, tot_n, ypos, yp, len(rows) + 0.7


def panel_b(ax, et, ypos, ypooled, ytop):
    rows = []
    disc = load((f"{et}_v2filtered" + F4SUF), et)
    if disc is None:
        ax.text(.5, .5, "discovery table missing", ha="center", va="center",
                transform=ax.transAxes, color=INK2, fontsize=8)
        style(ax)
        return
    sig = disc[(disc.FDR < 0.05) & (disc.dpsi.abs() >= 0.10)]
    D = sig[["ekey", "dpsi", "se_adj"]].rename(columns={"dpsi": "x", "se_adj": "sx"})
    for coh, pat in REPL:
        # LBNL is deliberately absent from panel B. Its remaining 5 DCM / 4 NF are the
        # samples LEFT OVER after the 24 matched ones were excluded, so its slope is a
        # different comparison from the 0.58 ceiling, which is measured ON those 24.
        # Two different sample sets adjacent on one axis invites exactly the wrong
        # reading. It stays in panel A, where it is informative.
        if coh in NOT_INDEPENDENT:
            continue
        T = load(pat.format(ET=et), et)
        if T is None:
            continue
        m = D.merge(T[["ekey", "dpsi", "se_adj"]].rename(
            columns={"dpsi": "y", "se_adj": "sy"}), on="ekey")
        sl = deming(m.x.values, m.y.values, m.sx.values, m.sy.values)
        if np.isfinite(sl):
            rows.append({"cohort": coh, "slope": sl, "n": len(m),
                         "indep": coh not in NOT_INDEPENDENT})
    if not rows:
        ax.text(.5, .5, "no slopes computed", ha="center", va="center",
                transform=ax.transAxes, color=INK2, fontsize=8)
        style(ax)
        return
    # Rows take panel A's y positions, not their own ranking. Ordering each panel by
    # its own values swaps Brigham and Colorado between them, and a reader scanning
    # across would silently compare the wrong cohorts.
    R = pd.DataFrame(rows)
    R = R[R.cohort.isin(ypos)].copy()
    R["y"] = R.cohort.map(ypos)
    R = R.sort_values("y").reset_index(drop=True)
    ax.axvspan(0, CEIL_PLATFORM, color=ORANGE, alpha=0.07, zorder=0, linewidth=0)
    # Ceiling labels sit on the row panel A uses for the pooled estimate, which is
    # empty here -- that keeps them clear of both the data and the panel title.
    # 162: the labels are derived from the constants. They were a SECOND
    # hard-coded copy of the same two numbers, so a re-based ceiling would have
    # been drawn with its old label.
    for xv, lab, col in ((CEIL_PLATFORM, f"platform {CEIL_PLATFORM:.2f}", ORANGE),
                         (CEIL_PEOPLE, f"people {CEIL_PEOPLE:.2f}", BLUE)):
        ax.axvline(xv, color=col, linewidth=1.3, linestyle=(0, (4, 3)), zorder=1)
        # 163: at 125 bp the ceilings are 0.69 and 0.97 and two centred labels
        # overlap. Anchor the lower one right and the upper one left so they
        # separate as the lines converge.
        _lo = min(CEIL_PLATFORM, CEIL_PEOPLE)
        _ha = "right" if xv == _lo else "left"
        _dx = -0.015 if _ha == "right" else 0.015
        ax.text(xv + _dx, ypooled, lab, fontsize=6.3, color=col, ha=_ha,
                va="center",
                bbox=dict(facecolor=SURF, edgecolor="none", pad=1.5))
    for r in R.itertuples():
        ax.plot([0, r.slope], [r.y, r.y], color=BLUE, linewidth=1.4, zorder=2)
        ax.plot([r.slope], [r.y], "o", ms=5.5, color=BLUE, zorder=3,
                markeredgecolor=SURF, markeredgewidth=1.1)
        ax.text(r.slope + 0.035, r.y, f"{r.slope:.2f}", fontsize=6.6,
                color=INK, va="center", ha="left")
    ax.set_yticks(list(R.y))
    ax.set_yticklabels(list(R.cohort), fontsize=7.5, color=INK)
    ax.set_xlim(0, 1.18)
    ax.set_ylim(-1.4, ytop)
    ax.set_xlabel("Deming slope against MAGNet", fontsize=7.5, color=INK2)
    ax.xaxis.grid(True, color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)
    style(ax)
    ax.set_title("B   Bounded by measurement", loc="left",
                 fontsize=8.6, fontweight="bold", pad=8)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--conc", default="direction_concordance_dedup.tsv")
    ap.add_argument("--event-type", default="SE")
    ap.add_argument("--out", default="F4_replication")
    a = ap.parse_args()
    if not os.path.exists(a.conc):
        sys.exit(f"[fatal] {a.conc} not found -- run 33_rebuild_summaries_dedup.py")
    C = pd.read_csv(a.conc, sep="\t")
    need = {"cohort", "n", "same", "pct"}
    if need - set(C.columns):
        sys.exit(f"[fatal] {a.conc} needs columns {need}, has {list(C.columns)}")
    print(f"[info] panel A from {a.conc}; panel B recomputed from the _dedup tables")

    plt.rcParams.update({"font.family": "sans-serif", "axes.linewidth": 0.8,
                         "figure.facecolor": SURF, "savefig.facecolor": SURF})
    fig = plt.figure(figsize=(7.4, 3.45))
    gs = gridspec.GridSpec(1, 2, width_ratios=[1.12, 1.0], wspace=0.42,
                           left=0.105, right=0.945, top=0.855, bottom=0.345)
    ppct, k, n, ypos, ypooled, ytop = panel_a(fig.add_subplot(gs[0]), C)
    panel_b(fig.add_subplot(gs[1]), a.event_type, ypos, ypooled, ytop)
    fig.text(0.105, 0.165,
             "* LBNL is not an independent cohort: 24 of its 33 samples are the same "
             "individuals as MAGNet samples. It is shown in panel A for\n"
             "   comparison, excluded from the pooled estimate, and omitted from panel "
             "B — its 9 remaining samples are a different comparison\n"
             "   from the reference lines, which are measured ON those 24 matched "
             "pairs. Slopes are interpretable only against the ceilings.\n"
             "   Rows are in the same order in both panels, so each cohort can be read "
             "across.",
             fontsize=6.1, color=MUTED, ha="left", va="top", linespacing=1.6)
    for ext in ("pdf", "png"):
        fig.savefig(f"{a.out}.{ext}", dpi=400)
    print(f"[write] {a.out}.pdf and {a.out}.png")
    print(f"\n  pooled (LBNL excluded): {k:,}/{n:,} = {ppct:.1f}%")
    print("\n  NOW LOOK AT THE PNG before using it.")


if __name__ == "__main__":
    main()
