#!/usr/bin/env python3
"""
39_figure2_discovery.py -- FIGURE 2: the discovery signature in MAGNet.

  A  Volcano, SE. dPSI against -log10(FDR), with the two thresholds drawn and the
     strongest genes labelled. Up and down are a POLARITY, so they get a diverging
     pair (orange / blue) with non-significant events in neutral grey.
  B  Event-type breakdown. Significant counts per type, exon-level, read from the
     frozen master_numbers_flagged.tsv, with the tested counts alongside.
  C  PSI heatmap of the top events across all MAGNet samples, grouped DCM vs NF.
     Optional -- it needs the rMATS table, which is large. --no-heatmap skips it.

NOTHING HERE WAS AFFECTED BY THE GENOTYPE EXCLUSIONS. No MAGNet sample was removed
(the keep rule always retained the discovery copy, and no MAGNet-MAGNet duplicate
exists), so every number in this figure is the one already locked in
master_numbers_flagged.tsv: 926 SE / 231 MXE / 90 A3SS / 84 A5SS / 114 RI.

Colours are the validated palette: slot 1 blue #2a78d6, slot 2 orange #eb6834.
Panel C's diverging map is built from those two poles through a neutral grey
midpoint -- never a rainbow, and never a hue at the midpoint.

Usage:
  python 39_figure2_discovery.py
  python 39_figure2_discovery.py --no-heatmap        # skip the rMATS read
  python 39_figure2_discovery.py --top-heat 30 --label-genes 12
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
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

BLUE = "#2a78d6"
ORANGE = "#eb6834"
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#8a8a85"
GRID = "#e4e4e0"
SURF = "#fcfcfb"
NULLGREY = "#c9c8c3"

# minimum between-sample PSI standard deviation for a heatmap row. Below this the
# event is effectively constant across all samples and z-scoring turns it into either
# a blank band (sd exactly 0 -> NaN) or pure amplified noise. Overridable with
# --min-psi-sd; 0 restores the old behaviour.
MIN_PSI_SD = 0.02

# a heatmap row with missing samples reads as noise no matter how large its effect,
# so rows must have PSI in at least this fraction of samples. rMATS leaves PSI as NA
# wherever an event has too few reads, and lowly-expressed loci are full of them.
MIN_FINITE = 0.90

TYPES = ["SE", "MXE", "A3SS", "A5SS", "RI"]
NKEEP = {"SE": 3, "MXE": 5, "A3SS": 5, "A5SS": 5, "RI": 3}
DIVERGING = LinearSegmentedColormap.from_list(
    "psi", [BLUE, "#dcdcd8", ORANGE], N=256)
# samples with no PSI call are drawn in a flat off-white, distinct from both poles and
# from the neutral midpoint, so a gap reads as "not measured" rather than as a value
DIVERGING.set_bad("#f2f1ee")


def style(ax):
    ax.set_facecolor(SURF)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(GRID)
        ax.spines[s].set_linewidth(0.8)
    ax.tick_params(colors=INK2, labelsize=7, length=3, width=0.8)


def load(base, et):
    for p in (f"{base}_dedup.all_events.tsv", f"{base}.all_events.tsv"):
        if os.path.exists(p) and os.path.getsize(p):
            t = pd.read_csv(p, sep="\t")
            d = next((c for c in ("dPSI_adj", "dPSI") if c in t.columns), None)
            if d is None or "key" not in t.columns:
                return None, p
            n = NKEEP[et]
            t["ekey"] = [":".join(str(k).split(":")[:n]) for k in t["key"]]
            t = t.loc[t.groupby("ekey")["FDR"].idxmin()].copy()
            return t.rename(columns={d: "dpsi"}), p
    return None, base


def panel_a(ax, T, et, fdr, mind, nlab):
    nfloor = int((T.FDR <= 0).sum())
    if nfloor:
        print(f"  [note] {nfloor} events have FDR reported as exactly 0; they are "
              f"drawn at the\n         top of the axis (-log10 capped at 300). Say "
              f"'FDR < 1e-300' in the legend,\n         not 'FDR = 0'.")
    y = -np.log10(T.FDR.clip(lower=1e-300))
    sig = (T.FDR < fdr) & (T.dpsi.abs() >= mind)
    ns = ~sig
    ax.scatter(T.dpsi[ns], y[ns], s=3, color=NULLGREY, alpha=0.45, linewidth=0,
               zorder=1, rasterized=True)
    for m, col in ((sig & (T.dpsi > 0), ORANGE), (sig & (T.dpsi < 0), BLUE)):
        ax.scatter(T.dpsi[m], y[m], s=5, color=col, alpha=0.75, linewidth=0,
                   zorder=2, rasterized=True)
    ax.axhline(-np.log10(fdr), color=MUTED, linewidth=0.9, linestyle=(0, (4, 3)),
               zorder=1)
    for xv in (-mind, mind):
        ax.axvline(xv, color=MUTED, linewidth=0.9, linestyle=(0, (4, 3)), zorder=1)

    top = T[sig].assign(score=y[sig] * T.dpsi[sig].abs())
    top = (top.sort_values("score", ascending=False)
           .drop_duplicates("gene").head(nlab))
    # The x limit must contain the labelled points. A percentile alone clips exactly
    # the largest effects -- which are the ones being named -- and the labels vanish
    # silently with no error.
    lim = float(np.nanpercentile(T.dpsi.abs(), 99.8)) * 1.12
    if len(top):
        # extra margin so the leader-line labels sit inside the axes rather than
        # running off the left edge into the tick labels
        lim = max(lim, float(top.dpsi.abs().max()) * 1.42)
    ax.set_xlim(-lim, lim)
    ax.set_ylim(0, float(y.max()) * 1.08)
    # Greedy vertical de-overlap, per side. The strongest events cluster at the top of
    # the axis, so labels placed at their own y collide into an unreadable smear; each
    # side is walked top-down and any label too close to the previous one is pushed
    # below it, with a leader line back to the point.
    ylo, yhi = ax.get_ylim()
    dymin = (yhi - ylo) * 0.052
    top = top.assign(yv=-np.log10(top.FDR.clip(lower=1e-300)))
    for sign in (1, -1):
        sub = top[np.sign(top.dpsi) == sign].sort_values("yv", ascending=False)
        last = None
        for r in sub.itertuples():
            yy = r.yv if last is None else min(r.yv, last - dymin)
            last = yy
            ax.annotate(r.gene, xy=(r.dpsi, r.yv),
                        xytext=(r.dpsi + (0.055 * lim * sign), yy),
                        fontsize=5.8, color=INK,
                        ha="left" if sign > 0 else "right", va="center",
                        arrowprops=dict(arrowstyle="-", color=MUTED, linewidth=0.5,
                                        shrinkA=0, shrinkB=2))
    ax.set_xlabel(f"ΔPSI  DCM vs non-failing", fontsize=7.5, color=INK2)
    ax.set_ylabel("−log₁₀ FDR", fontsize=7.5, color=INK2)
    ax.grid(True, color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)
    style(ax)
    nup, ndn = int((sig & (T.dpsi > 0)).sum()), int((sig & (T.dpsi < 0)).sum())
    # The counts live in the title. Inside the axes they collide with the lobes at the
    # bottom and with the labelled points at the top -- there is no reliably empty
    # region in a volcano.
    ax.set_title(f"A   Discovery ({et}):  {ndn:,} down, {nup:,} up", loc="left",
                 fontsize=8.6, fontweight="bold", pad=8)
    return nup, ndn


def panel_b(ax, M, fdr):
    sigcol = "sig_exon" if "sig_exon" in M.columns else "significant"
    tescol = "tested_exon" if "tested_exon" in M.columns else "tested"
    rows = []
    for et in TYPES:
        q = M[(M.cohort == "MAGNet") & (M.contrast == "DCM vs NF")
              & (M.event_type == et)]
        if q.empty:
            continue
        r = q.iloc[0]
        if str(r.get("reliability", "")) == "exclude":
            print(f"  [skip] MAGNet {et} flagged exclude")
            continue
        rows.append({"et": et, "n": int(r[sigcol]), "tested": int(r[tescol])})
    if not rows:
        ax.text(.5, .5, "no rows", ha="center", va="center", transform=ax.transAxes,
                color=INK2, fontsize=8)
        style(ax)
        return
    R = pd.DataFrame(rows)[::-1].reset_index(drop=True)
    y = np.arange(len(R))
    ax.barh(y, R.n, height=0.6, color=BLUE, linewidth=0)
    for i, r in R.iterrows():
        ax.text(r.n + max(R.n) * 0.02, i, f"{r.n:,}", fontsize=7, color=INK,
                va="center", ha="left", fontweight="bold")
        ax.text(r.n + max(R.n) * 0.02, i - 0.29, f"of {r.tested:,}", fontsize=5.9,
                color=MUTED, va="center", ha="left")
    ax.set_yticks(y)
    ax.set_yticklabels(R.et, fontsize=7.5, color=INK)
    ax.set_xlim(0, max(R.n) * 1.33)
    ax.set_xlabel("significant events", fontsize=7.5, color=INK2)
    ax.xaxis.grid(True, color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)
    style(ax)
    ax.set_title(f"B   All event types", loc="left", fontsize=8.6,
                 fontweight="bold", pad=8)
    print(f"  [info] event-type totals: "
          + ", ".join(f"{r.et} {r.n:,}" for r in R[::-1].itertuples())
          + f"   total {int(R.n.sum()):,}")


def panel_c(ax, cax, T, et, mats, ntop, fdr, mind, minsd=MIN_PSI_SD,
            minfin=MIN_FINITE):
    """PSI heatmap. Needs the rMATS table for per-sample IncLevel values."""
    if not (os.path.exists(mats) and os.path.getsize(mats)):
        ax.text(.5, .5, f"{mats} not found\n(--no-heatmap to skip this panel)",
                ha="center", va="center", transform=ax.transAxes, color=INK2,
                fontsize=7.5)
        style(ax)
        cax.axis("off")
        return
    scol = {"SE": "exonStart_0base", "MXE": "1stExonStart_0base",
            "A3SS": "longExonStart_0base", "A5SS": "longExonStart_0base",
            "RI": "riExonStart_0base"}[et]
    ecol = scol.replace("Start_0base", "End").replace("Start", "End")
    try:
        head = pd.read_csv(mats, sep="\t", nrows=0)
        use = [c for c in ("chr", scol, ecol, "IncLevel1", "IncLevel2")
               if c in head.columns]
        if len(use) < 5:
            raise KeyError(f"need chr/{scol}/{ecol}/IncLevel1/IncLevel2; "
                           f"have {list(head.columns)[:12]}...")
        R = pd.read_csv(mats, sep="\t", usecols=use)
    except Exception as e:
        ax.text(.5, .5, f"could not read {mats}\n{e}", ha="center", va="center",
                transform=ax.transAxes, color=INK2, fontsize=6.5)
        style(ax)
        cax.axis("off")
        return
    sig = T[(T.FDR < fdr) & (T.dpsi.abs() >= mind)].copy()
    if not {"chr", "exonStart", "exonEnd"} <= set(sig.columns):
        ax.text(.5, .5, "phase-2 table has no chr/exonStart/exonEnd to match on",
                ha="center", va="center", transform=ax.transAxes, color=INK2,
                fontsize=7)
        style(ax)
        cax.axis("off")
        return
    # ROW SELECTION. FDR alone was wrong here. FDR is a function of effect AND of how
    # tightly the event is measured, so the lowest-FDR events include near-invariant
    # ones: an exon sitting at PSI 0.99 in everybody with tiny scatter clears FDR<0.05
    # on a dPSI of 0.10 and then renders as a flat or blank band once z-scored (sd = 0
    # -> all NaN -> white). The panel exists to SHOW the contrast, so rank by |dPSI|
    # among events that already passed FDR, and require a real amount of between-sample
    # PSI spread. Also collapse to one event per gene: two rows labelled DCN are
    # ambiguous to a reader and waste half the panel.
    sig = sig.reindex(sig.dpsi.abs().sort_values(ascending=False).index)
    R = R.rename(columns={scol: "exonStart", ecol: "exonEnd"})
    m = sig.merge(R.drop_duplicates(["chr", "exonStart", "exonEnd"]),
                  on=["chr", "exonStart", "exonEnd"], how="inner")
    if len(m) < 5:
        ax.text(.5, .5, f"only {len(m)} of {len(sig)} top events matched the rMATS "
                        f"table\n(coordinate convention mismatch)",
                ha="center", va="center", transform=ax.transAxes, color=INK2,
                fontsize=7)
        style(ax)
        cax.axis("off")
        return

    def parse(s):
        return np.array([np.nan if v in ("NA", "", "nan") else float(v)
                         for v in str(s).split(",")])

    # ROW SELECTION, take two. Ranking by |dPSI| replaced one failure mode with its
    # opposite: the largest-|dPSI| events are lncRNAs, antisense and unnamed ENSG loci
    # whose PSI is estimated from a handful of reads, so their between-sample scatter
    # swamps the group difference and the panel reads as salt-and-pepper.
    # What the eye actually judges in a heatmap is SEPARATION -- group difference
    # relative to within-group spread -- so rank on exactly that: Cohen's d computed
    # from the same per-sample PSI values the panel draws. FDR<0.05 and |dPSI|>=0.10
    # still gate entry, so every row is a project-significant event; d only decides
    # which of them are legible. Also require near-complete PSI, since a row riddled
    # with missing samples reads as noise whatever its effect size.
    cand, n1, drop_cov, drop_flat = [], None, 0, 0
    for _, r in m.iterrows():
        a, b = parse(r.IncLevel1), parse(r.IncLevel2)
        v = np.concatenate([a, b])
        fin = float(np.isfinite(v).mean())
        if fin < minfin:
            drop_cov += 1
            continue
        if np.nanstd(v) < minsd:
            drop_flat += 1
            continue
        n1 = len(a)
        s1, s2 = np.nanstd(a, ddof=1), np.nanstd(b, ddof=1)
        na, nb = int(np.isfinite(a).sum()), int(np.isfinite(b).sum())
        if na < 5 or nb < 5:
            drop_cov += 1
            continue
        sp = np.sqrt(((na - 1) * s1 ** 2 + (nb - 1) * s2 ** 2) / max(na + nb - 2, 1))
        if not np.isfinite(sp) or sp <= 0:
            drop_flat += 1
            continue
        d = (np.nanmean(a) - np.nanmean(b)) / sp
        if not np.isfinite(d):
            drop_flat += 1
            continue
        cand.append({"absd": abs(float(d)), "d": float(d), "dpsi": float(r.dpsi),
                     "gene": str(r.gene), "fin": fin, "v": v})
    if len(cand) < 5:
        ax.text(.5, .5, f"too few events with usable, non-flat PSI\n"
                        f"({drop_cov} dropped on coverage, {drop_flat} as invariant)",
                ha="center", va="center", transform=ax.transAxes, color=INK2,
                fontsize=7.5)
        style(ax)
        cax.axis("off")
        return
    cand.sort(key=lambda c: -c["absd"])
    # one row per gene: two rows labelled DCN are ambiguous and waste half the panel
    seen_g, dedup, ndup = set(), [], 0
    for c in cand:
        if c["gene"] in seen_g:
            ndup += 1
            continue
        seen_g.add(c["gene"])
        dedup.append(c)
    cand = dedup[:ntop]
    # rows ordered by signed dPSI: everything gained in DCM together, everything lost
    # together, so the block structure is visible instead of interleaved
    cand.sort(key=lambda c: -c["dpsi"])
    ok = [c["gene"] for c in cand]
    P = np.vstack([c["v"] for c in cand])
    if drop_cov or drop_flat:
        print(f"  [info] heatmap: dropped {drop_cov} event(s) with <{minfin:.0%} finite "
              f"PSI or <5 usable samples in a group, and {drop_flat} with PSI "
              f"sd < {minsd} (invariant -- these are the blank/flat bands)")
    if ndup:
        print(f"  [info] heatmap: collapsed {ndup} extra event(s) on genes already "
              f"represented (kept the best-separating event per gene)")
    print(f"  [info] heatmap rows ranked by |Cohen's d| on per-sample PSI; "
          f"d ranges {min(abs(c['d']) for c in cand):.2f} to "
          f"{max(abs(c['d']) for c in cand):.2f}")
    print(f"        {'gene':<18}{'dPSI':>8}{'Cohen d':>10}{'PSI finite':>12}")
    for c in cand:
        print(f"        {c['gene']:<18}{c['dpsi']:>8.3f}{c['d']:>10.2f}"
              f"{c['fin']:>11.0%}")
    # z-score per event so the map shows the DCM/NF contrast, not baseline PSI level
    mu = np.nanmean(P, axis=1, keepdims=True)
    sd = np.nanstd(P, axis=1, keepdims=True)
    Z = (P - mu) / np.where(sd > 0, sd, np.nan)
    im = ax.imshow(Z, aspect="auto", cmap=DIVERGING,
                   norm=TwoSlopeNorm(vcenter=0, vmin=-2.5, vmax=2.5),
                   interpolation="nearest")
    ax.axvline(n1 - 0.5, color=INK, linewidth=1.1)
    ax.set_yticks(np.arange(len(ok)))
    ax.set_yticklabels(ok, fontsize=5.2, color=INK)
    ax.set_xticks([n1 / 2, n1 + (P.shape[1] - n1) / 2])
    ax.set_xticklabels([f"DCM  (n={n1})", f"non-failing  (n={P.shape[1] - n1})"],
                       fontsize=7, color=INK)
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    short = "" if len(ok) >= ntop else f" of {ntop} requested"
    ax.set_title(f"C   {len(ok)} best-separating significant events{short}, one per "
                 f"gene, per-sample PSI (z-scored per event)",
                 loc="left", fontsize=8.6, fontweight="bold", pad=8)
    cb = plt.colorbar(im, cax=cax)
    cb.set_label("z-scored PSI", fontsize=6.5, color=INK2)
    cb.ax.tick_params(labelsize=6, colors=INK2, length=2)
    cb.outline.set_visible(False)
    print(f"  [info] heatmap: {len(ok)} events x {P.shape[1]} samples")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--master", default="master_numbers_flagged.tsv")
    ap.add_argument("--discovery", default="{ET}_v2filtered")
    ap.add_argument("--mats", default="rmats_magnet/{ET}.MATS.JCEC.txt")
    ap.add_argument("--event-type", default="SE")
    ap.add_argument("--fdr", type=float, default=0.05)
    ap.add_argument("--min-dpsi", type=float, default=0.10)
    ap.add_argument("--top-heat", type=int, default=25)
    ap.add_argument("--min-finite", type=float, default=MIN_FINITE,
                    help="heatmap rows need at least this fraction of samples with a "
                         "usable PSI (default 0.90)")
    ap.add_argument("--min-psi-sd", type=float, default=MIN_PSI_SD,
                    help="drop heatmap rows whose per-sample PSI sd is below this "
                         "(default 0.02; 0 keeps invariant events)")
    ap.add_argument("--label-genes", type=int, default=10)
    ap.add_argument("--no-heatmap", dest="heat", action="store_false")
    ap.add_argument("--out", default="F2_discovery")
    a = ap.parse_args()
    et = a.event_type

    if not os.path.exists(a.master):
        sys.exit(f"[fatal] {a.master} not found -- run 24 then 26")
    M = pd.read_csv(a.master, sep="\t")
    T, p = load(a.discovery.format(ET=et), et)
    if T is None:
        sys.exit(f"[fatal] discovery table not usable: {p}")
    print(f"[info] volcano from {p} ({len(T):,} exon-level {et} events); "
          f"panel B from {a.master}")

    plt.rcParams.update({"font.family": "sans-serif", "axes.linewidth": 0.8,
                         "figure.facecolor": SURF, "savefig.facecolor": SURF})
    if a.heat:
        fig = plt.figure(figsize=(7.4, 5.4))
        gs = gridspec.GridSpec(2, 3, width_ratios=[1.3, 1.0, 0.045],
                               height_ratios=[1.0, 1.25], wspace=0.42, hspace=0.34,
                               left=0.115, right=0.93, top=0.945, bottom=0.07)
        axa, axb = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])
        axc, cax = fig.add_subplot(gs[1, :2]), fig.add_subplot(gs[1, 2])
    else:
        fig = plt.figure(figsize=(7.4, 2.9))
        gs = gridspec.GridSpec(1, 2, width_ratios=[1.3, 1.0], wspace=0.42,
                               left=0.115, right=0.965, top=0.86, bottom=0.2)
        axa, axb = fig.add_subplot(gs[0]), fig.add_subplot(gs[1])
        axc = cax = None
    panel_a(axa, T, et, a.fdr, a.min_dpsi, a.label_genes)
    panel_b(axb, M, a.fdr)
    if axc is not None:
        panel_c(axc, cax, T, et, a.mats.format(ET=et), a.top_heat, a.fdr, a.min_dpsi,
                a.min_psi_sd, a.min_finite)
    for ext in ("pdf", "png"):
        fig.savefig(f"{a.out}.{ext}", dpi=400)
    print(f"[write] {a.out}.pdf and {a.out}.png")
    print("\n  NOW LOOK AT THE PNG -- gene labels in a volcano collide easily.")


if __name__ == "__main__":
    main()
