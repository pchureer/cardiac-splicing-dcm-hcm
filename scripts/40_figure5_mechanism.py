#!/usr/bin/env python3
"""
40_figure5_mechanism.py -- FIGURE 5: the mechanism null, calibrated against a
housekeeping-gene floor.

  A  RBM20 expression is HIGHER in DCM, not lower. The positive finding that kills the
     loss-of-function framing before the null is argued.
  B  RBM20 targets against matched background: no enrichment, SE or MXE.
  C  The inverted 11-candidate scan as a delta plot, with GAPDH drawn as the
     calibration floor. THE panel -- it answers "if not RBM20, then what?".
  D  Why the p-values are not to be over-read: background median |t| before and after
     adjusting for library size and TIN, against the null expectation of ~0.67.

PROVENANCE, AND WHY THIS SCRIPT HAS A CONSTANTS BLOCK
-----------------------------------------------------
Scripts 03 and 23 produced their summaries months ago and the audited numbers live in
RESULT_rbm20_mechanism_null.md. Those numbers are frozen and are shipped below in
FROZEN_*, so the figure always builds and every plotted value is traceable to that doc.

Where the per-event tables still exist on disk the script PREFERS them, draws the real
distributions, and -- for panel C -- cross-checks the recomputed deltas against the
frozen ones, printing a loud warning on any disagreement above --tol. It never silently
mixes the two: the PROVENANCE block printed at the end names the source of every panel.

A figure built from constants is not a worse figure, but it IS a different claim about
where the numbers came from. Say in the legend which mode produced it; the printed
PROVENANCE block is written to be pasted into the methods.

Usage:
  export PATH=/data/apps/anaconda3/envs/mamba/envs/rnaseq_tools/bin:$PATH
  cd /data/project/chureeratp/RNA/HEART
  python 40_figure5_mechanism.py
  python 40_figure5_mechanism.py --scan rbp_scan_SE.tsv rbp_scan_MXE.tsv
  python 40_figure5_mechanism.py --frozen-only      # ignore anything on disk
  python 40_figure5_mechanism.py --list             # just report what it can find
"""
import argparse
import glob
import os
import re
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
NULLGREY = "#c9c8c3"

# ---------------------------------------------------------------------------
# FROZEN NUMBERS -- RESULT_rbm20_mechanism_null.md, revised 13 Sep 2026.
# Do not edit these to make a panel look better. If a number here is wrong, the doc is
# wrong, and the doc is what the manuscript cites.
# ---------------------------------------------------------------------------

# section 3: RBM20 differential expression, DCM vs NF, in SD units
FROZEN_DE = [("unadjusted", 0.243, 0.027),
             ("+ covariates", 0.382, 2.7e-04),
             ("+ composition", 0.340, 0.0025)]

# section 2: target vs matched background, median |t|, within-DCM, cell-adjusted
FROZEN_TARGET_BG = {"SE": (0.92, 0.96, 0.76), "MXE": (1.48, 1.65, 0.99)}

# section 4: inverted scan. delta = median |t| target - median |t| matched background
#            gene -> (SE delta, SE p, MXE delta, MXE p)
FROZEN_SCAN = {
    "PTBP1":    (+0.08, 8.4e-05, +0.10, 0.0019),
    "QKI":      (+0.07, 0.0019, -0.05, 0.69),
    "HNRNPF":   (+0.07, 0.0013, +0.00, 0.37),
    "ZRANB2":   (+0.06, 0.0016, -0.03, 0.71),
    "HNRNPA0":  (+0.05, 8.0e-04, -0.03, 0.94),
    "MBNL1":    (+0.05, 0.0021, -0.02, 0.72),
    "GAPDH":    (+0.05, 0.0046, +0.10, 0.0017),
    "SRSF1":    (+0.04, 0.013, -0.07, 0.87),
    "RBM20":    (+0.03, 0.015, +0.09, 0.0067),
    "CELF1":    (+0.01, 0.18, -0.01, 0.50),
    "RBM5":     (-0.01, 0.28, -0.01, 0.35),
}
FROZEN_SCAN_N = {"SE": (2504, 18685), "MXE": (2755, 20402)}  # target, background events

# section 1: background median |t| before and after adding log10_libsize + TIN.median.
# Two reported variants per event type (the doc gives them as "a / b").
FROZEN_TECH = {"SE": [(1.01, 0.93), (1.09, 0.96)],
               "MXE": [(2.01, 1.28), (2.65, 1.65)]}
NULL_EXPECT = 0.67  # median |t| expected under no association

CONTROL = "GAPDH"


def style(ax):
    ax.set_facecolor(SURF)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(GRID)
        ax.spines[s].set_linewidth(0.8)
    ax.tick_params(colors=INK2, labelsize=7, length=3, width=0.8)


def ids(path):
    """Sample IDs out of a bam list, the same way scripts 32 and 36 do it."""
    out = []
    for p in re.split(r"[,\n]", open(path).read()):
        p = p.strip()
        if not p:
            continue
        m = re.search(r"(SRR\d+|ERR\d+|SAMEA\d+|SAMN\d+)", p)
        out.append(m.group(1) if m else os.path.basename(os.path.dirname(p)))
    return sorted(set(out))


def lower(cols):
    return {str(c).lower(): c for c in cols}


# ---------------------------------------------------------------------------
# discovery. Column-driven, not name-driven: these files were written months ago and
# guessing their names is how a panel ends up silently empty.
# ---------------------------------------------------------------------------

def candidate_files(extra):
    seen, out = set(), []
    pats = list(extra) + [
        "*rbp*scan*.tsv", "*scan*rbp*.tsv", "*signature_rbp*.tsv",
        "*rbm20*.tsv", "*rbm20*.txt", "*target*background*.tsv",
        "*mechanism*.tsv", "rbp_*.tsv", "scan_*.tsv",
        "*/rbp*scan*.tsv", "*/rbm20*.tsv",
    ]
    for p in pats:
        for f in sorted(glob.glob(p)):
            if f not in seen and os.path.isfile(f) and os.path.getsize(f):
                seen.add(f)
                out.append(f)
    return out


def sniff(path, nrows=400):
    try:
        t = pd.read_csv(path, sep=None, engine="python", nrows=nrows)
    except Exception:
        return None
    return t


def find_scan(files, verbose=True):
    """A scan table: one row per candidate gene, with a delta (or the two medians)."""
    for f in files:
        t = sniff(f)
        if t is None or not len(t):
            continue
        L = lower(t.columns)
        gcol = next((L[k] for k in ("gene", "candidate", "rbp", "symbol") if k in L), None)
        if gcol is None:
            continue
        dcol = next((L[k] for k in ("delta", "enrichment", "diff", "delta_t",
                                    "median_diff") if k in L), None)
        tg = next((L[k] for k in ("median_t_target", "target_median_t", "target_t",
                                  "median_target") if k in L), None)
        bg = next((L[k] for k in ("median_t_background", "background_median_t",
                                  "background_t", "median_background") if k in L), None)
        if dcol is None and not (tg and bg):
            continue
        genes = set(str(x).upper() for x in t[gcol])
        if CONTROL not in genes:
            if verbose:
                print(f"  [skip] {f} looks like a scan table but has no {CONTROL} row; "
                      f"without the calibration floor it cannot be used")
            continue
        out = pd.DataFrame({"gene": [str(x).upper() for x in t[gcol]]})
        out["delta"] = (t[dcol].astype(float) if dcol
                        else t[tg].astype(float) - t[bg].astype(float))
        pcol = next((L[k] for k in ("p", "pval", "p_value", "pvalue", "mwu_p")
                     if k in L), None)
        out["p"] = t[pcol].astype(float) if pcol else np.nan
        etcol = next((L[k] for k in ("event_type", "etype", "type") if k in L), None)
        if etcol:
            out["et"] = [str(x).upper() for x in t[etcol]]
        else:  # one file per event type -- take it from the file name
            m = re.search(r"(SE|MXE|A3SS|A5SS|RI)", os.path.basename(f), re.I)
            if not m:
                if verbose:
                    print(f"  [skip] {f}: no event-type column and none in the filename")
                continue
            out["et"] = m.group(1).upper()
        return out, f
    return None, None


def find_events(files, verbose=True):
    """A per-event table with an |t| column and a target/background label."""
    for f in files:
        t = sniff(f)
        if t is None or not len(t):
            continue
        L = lower(t.columns)
        tcol = next((L[k] for k in ("abs_t", "abst", "t", "tstat", "t_stat",
                                    "t_value") if k in L), None)
        scol = next((L[k] for k in ("set", "group", "class", "is_target", "target")
                     if k in L), None)
        if tcol is None or scol is None:
            continue
        etcol = next((L[k] for k in ("event_type", "etype", "type") if k in L), None)
        vals = set(str(x).lower() for x in t[scol].head(200))
        if not (vals & {"target", "targets", "1", "true", "background", "bg", "0",
                        "false"}):
            continue
        full = pd.read_csv(f, sep=None, engine="python")
        L = lower(full.columns)
        d = pd.DataFrame({"abs_t": pd.to_numeric(full[tcol], errors="coerce").abs()})
        s = full[scol].astype(str).str.lower()
        d["is_target"] = s.isin(["target", "targets", "1", "true", "yes"])
        if etcol:
            d["et"] = full[etcol].astype(str).str.upper()
        else:
            m = re.search(r"(SE|MXE)", os.path.basename(f), re.I)
            if not m:
                continue
            d["et"] = m.group(1).upper()
        d = d.dropna(subset=["abs_t"])
        if len(d) < 200 or d.is_target.sum() < 20:
            if verbose:
                print(f"  [skip] {f}: only {len(d)} rows / {int(d.is_target.sum())} "
                      f"targets -- too few to draw a distribution")
            continue
        return d, f
    return None, None


def find_expr(expr_path, b1, b2, verbose=True):
    """Per-sample RBM20, plus DCM/NF membership from the bam lists."""
    if not (expr_path and os.path.exists(expr_path) and os.path.getsize(expr_path)):
        return None, None
    try:
        E = pd.read_csv(expr_path, sep="\t")
    except Exception as e:
        if verbose:
            print(f"  [skip] {expr_path}: {e}")
        return None, None
    # genes may be columns or the index
    if "RBM20" in E.columns:
        scol = next((c for c in E.columns
                     if str(c).lower() in ("sample", "sample_id", "id", "run")), None)
        if scol is None:
            scol = E.columns[0]
        V = pd.DataFrame({"sample": E[scol].astype(str),
                          "rbm20": pd.to_numeric(E["RBM20"], errors="coerce")})
    else:
        gcol = next((c for c in E.columns
                     if str(c).lower() in ("gene", "gene_name", "symbol")), E.columns[0])
        hit = E[E[gcol].astype(str).str.upper() == "RBM20"]
        if not len(hit):
            if verbose:
                print(f"  [skip] {expr_path}: no RBM20 row or column")
            return None, None
        r = hit.iloc[0].drop(labels=[gcol])
        V = pd.DataFrame({"sample": [str(i) for i in r.index],
                          "rbm20": pd.to_numeric(r.values, errors="coerce")})
    for p in (b1, b2):
        if not (p and os.path.exists(p)):
            if verbose:
                print(f"  [skip] panel A needs {b1} and {b2} for group membership")
            return None, None
    g1, g2 = set(ids(b1)), set(ids(b2))
    V["group"] = np.where(V["sample"].isin(g1), "DCM",
                          np.where(V["sample"].isin(g2), "NF", ""))
    V = V[(V.group != "") & np.isfinite(V.rbm20)]
    if (V.group == "DCM").sum() < 10 or (V.group == "NF").sum() < 10:
        if verbose:
            print(f"  [skip] {expr_path}: only {(V.group=='DCM').sum()} DCM / "
                  f"{(V.group=='NF').sum()} NF samples matched the bam lists")
        return None, None
    return V, expr_path


# ---------------------------------------------------------------------------
# panels
# ---------------------------------------------------------------------------

def panel_a(ax, V):
    """RBM20 is UP in DCM. Real per-sample values when available, else the estimates."""
    if V is not None:
        z = (V.rbm20 - V.rbm20.mean()) / V.rbm20.std(ddof=1)
        rng = np.random.default_rng(0)
        ns = {}
        for i, (g, col) in enumerate((("NF", NULLGREY), ("DCM", BLUE))):
            v = z[V.group.values == g].values
            x = i + rng.normal(0, 0.075, len(v))
            ax.scatter(x, v, s=6, color=col, alpha=.55, linewidths=0, zorder=2)
            bp = ax.boxplot([v], positions=[i], widths=.42, showfliers=False,
                            patch_artist=True, zorder=3)
            for b in bp["boxes"]:
                b.set(facecolor="none", edgecolor=INK, linewidth=.9)
            for k in ("whiskers", "caps", "medians"):
                for b in bp[k]:
                    b.set(color=INK, linewidth=.9)
            ns[g] = len(v)
        # A handful of very low-RBM20 libraries otherwise set the y range and squash
        # the comparison into the top third of the panel. Clip the AXIS, never the
        # data: the points are still in the boxplot statistics, and the count of
        # off-axis points is stated on the panel so nothing is hidden.
        # Tukey range, not a percentile: with 3 outliers in 330 samples a 0.4th
        # percentile still sits on top of them and nothing gets clipped.
        q1, q3 = np.percentile(z, [25, 75])
        iqr = q3 - q1
        lo = max(float(np.min(z)), q1 - 1.9 * iqr)
        hi = min(float(np.max(z)), q3 + 1.9 * iqr)
        pad = (hi - lo) * 0.10
        off = int(((z < lo) | (z > hi)).sum())
        ax.set_ylim(lo - pad, hi + pad)
        if off:
            ax.text(.01, .01, f"{off} sample(s) beyond axis", transform=ax.transAxes,
                    fontsize=5.6, color=MUTED, ha="left", va="bottom")
        # and name them, because an RBM20 value 3+ SD low is a QC question, not a dot
        odd = V.assign(z=z.values).loc[np.abs(z.values) > 3]
        if len(odd):
            print(f"  [check] {len(odd)} sample(s) with |z(RBM20)| > 3 -- confirm these "
                  f"are real and not failed libraries:")
            for _, r in odd.sort_values("z").iterrows():
                print(f"          {r['sample']:<14}{r.group:<5}z = {r.z:+.2f}")
        ax.set_xticks([0, 1])
        ax.set_xticklabels([f"non-failing\nn={ns.get('NF', 0)}",
                            f"DCM\nn={ns.get('DCM', 0)}"], fontsize=7.5, color=INK)
        ax.set_ylabel("RBM20 expression (z-scored log$_2$CPM)", fontsize=7.5, color=INK2)
        ax.set_xlim(-.55, 1.55)
        est, p = FROZEN_DE[1][1], FROZEN_DE[1][2]
        ax.annotate(f"+{est:.3f} SD in DCM\ncovariate-adjusted, p = {p:.1e}",
                    xy=(.99, .02), xycoords="axes fraction", ha="right", va="bottom",
                    fontsize=6.8, color=INK)
        sub = "per-sample values; adjusted estimate from script 03"
    else:
        y = np.arange(len(FROZEN_DE))[::-1]
        for i, (lab, est, p) in zip(y, FROZEN_DE):
            ax.plot([0, est], [i, i], color=GRID, linewidth=2.2, zorder=1,
                    solid_capstyle="round")
            ax.scatter([est], [i], s=46, color=BLUE, zorder=3, linewidths=0)
            ax.text(est + 0.012, i + .22, f"+{est:.3f}", fontsize=7, color=INK,
                    fontweight="bold", va="center")
            ax.text(est + 0.012, i - .22, f"p = {p:.2g}", fontsize=6, color=MUTED,
                    va="center")
        ax.axvline(0, color=INK, linewidth=.9)
        ax.set_yticks(y)
        ax.set_yticklabels([d[0] for d in FROZEN_DE], fontsize=7, color=INK)
        ax.set_xlim(-0.03, max(d[1] for d in FROZEN_DE) * 1.5)
        ax.set_ylim(-.7, len(FROZEN_DE) - .3)
        ax.set_xlabel("RBM20 in DCM vs non-failing (SD units)", fontsize=7.5, color=INK2)
        ax.xaxis.grid(True, color=GRID, linewidth=.6)
        ax.set_axisbelow(True)
        sub = "frozen estimates"
    style(ax)
    ax.set_title("A   RBM20 is UP in DCM, not down", loc="left", fontsize=8.6,
                 fontweight="bold", pad=8)
    return sub


def panel_b(ax, D):
    """RBM20 targets vs matched background. Null in both event types."""
    types = ["SE", "MXE"]
    if D is not None:
        pos, ticks, labs = [], [], []
        for i, et in enumerate(types):
            sub = D[D.et == et]
            if not len(sub):
                continue
            for j, (istg, col, lab) in enumerate(((False, NULLGREY, "background"),
                                                  (True, BLUE, "targets"))):
                v = sub.abs_t[sub.is_target == istg].values
                if not len(v):
                    continue
                x = i * 2.2 + j * 0.95
                bp = ax.boxplot([v], positions=[x], widths=.60, showfliers=False,
                                patch_artist=True)
                for b in bp["boxes"]:
                    b.set(facecolor=col, edgecolor=INK, linewidth=.8, alpha=.85)
                for k in ("whiskers", "caps", "medians"):
                    for b in bp[k]:
                        b.set(color=INK, linewidth=.8)
                pos.append(x)
                labs.append(f"{lab}\nn={len(v):,}")
            ax.text(i * 2.2 + .475, -0.155, et, transform=ax.get_xaxis_transform(),
                    ha="center", va="top", fontsize=7.8, color=INK,
                    fontweight="bold")
        ax.set_xticks(pos)
        ax.set_xticklabels(labs, fontsize=6.2, color=INK)
        ax.set_xlim(-.7, (len(types) - 1) * 2.2 + 1.65)
        ax.set_ylabel("|t|  PSI ~ RBM20, within DCM", fontsize=7.5, color=INK2)
        for i, et in enumerate(types):
            p = FROZEN_TARGET_BG[et][2]
            ax.text(i * 2.2 + .475, ax.get_ylim()[1] * .90, f"p = {p:.2f}",
                    ha="center", va="top", fontsize=6.8, color=INK)
        sub = "per-event |t| distributions"
    else:
        x = np.arange(len(types))
        w = .32
        tg = [FROZEN_TARGET_BG[e][0] for e in types]
        bg = [FROZEN_TARGET_BG[e][1] for e in types]
        ax.bar(x - w / 2, bg, w, color=NULLGREY, linewidth=0, label="matched background")
        ax.bar(x + w / 2, tg, w, color=BLUE, linewidth=0, label="RBM20 targets")
        for i, e in enumerate(types):
            ax.text(i - w / 2, bg[i] + .03, f"{bg[i]:.2f}", ha="center", fontsize=6.6,
                    color=INK2)
            ax.text(i + w / 2, tg[i] + .03, f"{tg[i]:.2f}", ha="center", fontsize=6.6,
                    color=INK, fontweight="bold")
            ax.text(i, max(tg[i], bg[i]) + .16, f"p = {FROZEN_TARGET_BG[e][2]:.2f}",
                    ha="center", fontsize=7, color=INK)
        ax.set_xticks(x)
        ax.set_xticklabels(types, fontsize=7.5, color=INK)
        ax.set_ylabel("median |t|, within DCM", fontsize=7.5, color=INK2)
        ax.set_ylim(0, max(tg + bg) * 1.35)
        ax.legend(fontsize=6.2, frameon=False, loc="upper left")
        ax.yaxis.grid(True, color=GRID, linewidth=.6)
        ax.set_axisbelow(True)
        sub = "frozen medians"
    style(ax)
    ax.set_title("B   RBM20 targets vs matched background", loc="left",
                 fontsize=8.6, fontweight="bold", pad=8)
    return sub


def panel_c(ax, S):
    """The inverted scan. GAPDH is the floor, drawn as a line, not a footnote."""
    src = "frozen deltas"
    tab = {g: {"SE": v[0], "MXE": v[2]} for g, v in FROZEN_SCAN.items()}
    if S is not None:
        got = {}
        for et in ("SE", "MXE"):
            sub = S[S.et == et]
            for _, r in sub.iterrows():
                got.setdefault(r.gene, {})[et] = float(r.delta)
        if CONTROL in got and len(got) >= 5:
            tab, src = got, "recomputed from the scan table"
    order = sorted(tab, key=lambda g: -tab[g].get("SE", -9))
    y = np.arange(len(order))[::-1]
    floor = {et: tab.get(CONTROL, {}).get(et, np.nan) for et in ("SE", "MXE")}
    for et, col, off in (("SE", BLUE, +.17), ("MXE", ORANGE, -.17)):
        f = floor[et]
        if np.isfinite(f):
            ax.axvline(f, color=col, linewidth=1.0, linestyle=(0, (2, 2)), zorder=1)
        for i, g in zip(y, order):
            v = tab[g].get(et, np.nan)
            if not np.isfinite(v):
                continue
            ax.plot([0, v], [i + off, i + off], color=GRID, linewidth=1.6, zorder=1,
                    solid_capstyle="round")
            ax.scatter([v], [i + off], s=26, color=col, zorder=3, linewidths=0)
    ax.axvline(0, color=INK, linewidth=.9, zorder=2)
    ax.set_yticks(y)
    ax.set_yticklabels(
        [f"$\\bf{{{g}}}$" if g == CONTROL else g for g in order],
        fontsize=6.4, color=INK)
    for i, g in zip(y, order):
        if g == CONTROL:
            ax.axhspan(i - .45, i + .45, color="#efeeea", zorder=0)
    ax.set_xlabel("enrichment $\\Delta$   (median |t| targets $-$ matched background)",
                  fontsize=7, color=INK2)
    ax.set_ylim(-.7, len(order) - .1)
    ax.xaxis.grid(True, color=GRID, linewidth=.6)
    ax.set_axisbelow(True)
    style(ax)
    h = [plt.Line2D([], [], marker="o", linestyle="", color=BLUE, markersize=4,
                    label="SE"),
         plt.Line2D([], [], marker="o", linestyle="", color=ORANGE, markersize=4,
                    label="MXE"),
         plt.Line2D([], [], linestyle=(0, (2, 2)), color=BLUE, linewidth=1,
                    label=f"{CONTROL} floor SE"),
         plt.Line2D([], [], linestyle=(0, (2, 2)), color=ORANGE, linewidth=1,
                    label=f"{CONTROL} floor MXE")]
    ax.legend(handles=h, fontsize=5.8, frameon=False, ncol=4, columnspacing=1.0,
              handletextpad=.4, loc="lower left", bbox_to_anchor=(0, 1.005, 1, .1),
              mode="expand", borderaxespad=0)
    # NOT "nothing beats the floor" -- four SE candidates sit above it, by 0.01-0.03.
    # The claim is that none SEPARATES from it, and an overstated title is the first
    # thing a reviewer checks against the panel.
    ax.set_title("C   Candidates against the housekeeping floor", loc="left",
                 fontsize=8.6, fontweight="bold", pad=18)
    return src, order, floor, tab


def panel_d(ax):
    """Why every p-value in C is inflated: shared technical variation."""
    xs = [0, 1]
    for et, col in (("SE", BLUE), ("MXE", ORANGE)):
        for k, (before, after) in enumerate(FROZEN_TECH[et]):
            ax.plot(xs, [before, after], color=col, linewidth=1.3, alpha=.9,
                    marker="o", markersize=4, zorder=3)
            if k == 0:
                ax.text(-0.06, before, et, ha="right", va="center", fontsize=7,
                        color=col, fontweight="bold")
    ax.axhline(NULL_EXPECT, color=INK, linewidth=.9, linestyle=(0, (3, 2)), zorder=2)
    ax.text(1.04, NULL_EXPECT, "null expectation\n≈0.67", fontsize=6, color=INK2,
            va="center", ha="left")
    ax.set_xticks(xs)
    ax.set_xticklabels(["before", "after\n+ libsize + TIN"], fontsize=7, color=INK)
    ax.set_xlim(-.35, 1.55)
    ax.set_ylabel("background median |t|", fontsize=7.5, color=INK2)
    ax.yaxis.grid(True, color=GRID, linewidth=.6)
    ax.set_axisbelow(True)
    style(ax)
    ax.set_title("D   Technical variation inflates everything", loc="left",
                 fontsize=8.6, fontweight="bold", pad=8)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--expr", default="magnet_expr.per_sample.tsv")
    ap.add_argument("--b1", default="dcm_bams.txt")
    ap.add_argument("--b2", default="nf_bams.txt")
    ap.add_argument("--scan", nargs="*", default=[],
                    help="explicit scan table(s) from script 23")
    ap.add_argument("--events", nargs="*", default=[],
                    help="explicit per-event table(s) from script 03")
    ap.add_argument("--frozen-only", action="store_true",
                    help="ignore files on disk; build entirely from the constants block")
    ap.add_argument("--list", dest="listonly", action="store_true",
                    help="report what could be found, draw nothing")
    ap.add_argument("--tol", type=float, default=0.02,
                    help="warn if a recomputed delta differs from the frozen one by more")
    ap.add_argument("--out", default="F5_mechanism")
    a = ap.parse_args()

    print("[scan] looking for per-event and scan tables")
    files = [] if a.frozen_only else candidate_files(a.scan + a.events)
    if files:
        print(f"  [info] {len(files)} candidate file(s): "
              + ", ".join(os.path.basename(f) for f in files[:8])
              + (" ..." if len(files) > 8 else ""))
    S, s_src = (None, None) if a.frozen_only else find_scan(files)
    D, d_src = (None, None) if a.frozen_only else find_events(files)
    V, v_src = (None, None) if a.frozen_only else find_expr(a.expr, a.b1, a.b2)

    if a.listonly:
        print(f"\n  panel A  {v_src or 'FROZEN (no usable per-sample RBM20)'}")
        print(f"  panel B  {d_src or 'FROZEN (no usable per-event table)'}")
        print(f"  panel C  {s_src or 'FROZEN (no usable scan table)'}")
        print(f"  panel D  FROZEN by construction (before/after summaries)")
        return

    plt.rcParams.update({"font.family": "sans-serif", "axes.linewidth": 0.8,
                         "figure.facecolor": SURF, "savefig.facecolor": SURF})
    fig = plt.figure(figsize=(7.6, 6.8))
    gs = gridspec.GridSpec(2, 2, wspace=0.40, hspace=0.52,
                           left=0.135, right=0.945, top=0.940, bottom=0.085)
    sa = panel_a(fig.add_subplot(gs[0, 0]), V)
    sb = panel_b(fig.add_subplot(gs[0, 1]), D)
    sc, order, floor, tab = panel_c(fig.add_subplot(gs[1, 0]), S)
    panel_d(fig.add_subplot(gs[1, 1]))
    for ext in ("pdf", "png"):
        fig.savefig(f"{a.out}.{ext}", dpi=400)
    print(f"[write] {a.out}.pdf and {a.out}.png")

    # -- the check that matters: did anything recomputed disagree with the doc? --
    if S is not None:
        bad = []
        for g, v in FROZEN_SCAN.items():
            for et, fz in (("SE", v[0]), ("MXE", v[2])):
                rc = tab.get(g, {}).get(et)
                if rc is not None and abs(rc - fz) > a.tol:
                    bad.append((g, et, fz, rc))
        if bad:
            print(f"\n  *** {len(bad)} delta(s) differ from RESULT_rbm20_mechanism_null.md "
                  f"by more than {a.tol}:")
            for g, et, fz, rc in bad:
                print(f"      {g:<10}{et:<5}doc {fz:+.3f}   file {rc:+.3f}")
            print("      DO NOT publish either number until you know which is current.")
        else:
            print(f"\n  [ok] every recomputed delta agrees with the doc to within {a.tol}")

    # -- what the reader needs to know, phrased for the methods --
    print("\n" + "=" * 78)
    print("PROVENANCE -- paste into the methods / legend")
    print("=" * 78)
    print(f"  panel A  {v_src or 'FROZEN: RESULT_rbm20_mechanism_null.md section 3'}"
          f"   ({sa})")
    print(f"  panel B  {d_src or 'FROZEN: RESULT_rbm20_mechanism_null.md section 2'}"
          f"   ({sb})")
    print(f"  panel C  {s_src or 'FROZEN: RESULT_rbm20_mechanism_null.md section 4'}"
          f"   ({sc})")
    print(f"  panel D  FROZEN: RESULT_rbm20_mechanism_null.md section 1")
    print("           NOTE: the doc reports two numbers per event type "
          "('1.01 / 1.09'), drawn\n           here as two lines. Confirm which model "
          "each is before writing the legend.")
    print(f"\n  panel C floor: {CONTROL} delta = "
          + ", ".join(f"{et} {floor[et]:+.2f}" for et in ("SE", "MXE")
                      if np.isfinite(floor[et])))
    for et in ("SE", "MXE"):
        ab = [(g, tab[g][et] - floor[et]) for g in order
              if g != CONTROL and np.isfinite(tab[g].get(et, np.nan))
              and tab[g][et] > floor[et] + 1e-9]
        ab.sort(key=lambda x: -x[1])
        print(f"  above the {et} floor: {len(ab)}/{len(order)-1}"
              + ("  " + ", ".join(f"{g} +{d:.2f}" for g, d in ab) if ab else "  none"))
    print("  Quote the MARGINS, not the count. Four SE candidates clear the floor by")
    print("  0.01-0.03 on a scale where the floor itself is +0.05; that is the point.")
    print(f"  target/background events: "
          + ", ".join(f"{et} {t:,}/{b:,}" for et, (t, b) in FROZEN_SCAN_N.items()))
    print("\n  The claim this figure supports is NOT 'no candidate is significant' --")
    print("  several are. It is that none is distinguishable from a housekeeping gene,")
    print("  and panel D shows why significance is the wrong reading.")
    print("\n  NOW LOOK AT THE PNG.")


if __name__ == "__main__":
    main()
