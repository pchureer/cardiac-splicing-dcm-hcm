#!/usr/bin/env python3
"""
37_figure3_etiology.py -- FIGURE 3: the signature is shared across etiologies.

THE FIGURE THAT CARRIES THE PAPER. Three panels:

  A  Significant-event counts. Disease-vs-control against disease-vs-other-etiology,
     at matched group sizes within each cohort. The claim is the RATIO.
  B  Power. ECDF of mde80 (the effect detectable with 80% power) for the etiology
     contrast, against the effects actually seen in the disease-vs-control contrast.
     This is what kills the "your null is just low power" objection.
  C  Effect sizes agree. MAGNet DCM-vs-NF against HCM-vs-NF on shared controls,
     with the identity line and the Deming slope.

DESIGN DECISIONS, so they are not re-litigated later
----------------------------------------------------
* Panel A uses a LINEAR axis, not log. The claim is "essentially zero", and a log
  axis would rescue the near-zero bars into visibility -- the opposite of the point.
  Every bar carries a direct label, so no value is lost to the reader.
* Counts are EXON-LEVEL at FDR<0.05 and |dPSI|>=0.10, read from
  master_numbers_flagged.tsv. Nothing is recomputed here: that file is the frozen
  source of truth, and a figure that recomputes its own numbers will eventually
  disagree with the text.
* Rows flagged `exclude` in that file are dropped automatically (Colorado MXE) and
  the script says so. It will not silently plot an excluded number.
* No dual axes anywhere. Grid and axes recessive; all text in ink, never series colour.
* Colours are slots 1 and 2 of the validated categorical palette (blue #2a78d6,
  orange #eb6834), the same pair as fig_split_control. Validated all-pairs, light
  mode: CVD dE 9.2, normal-vision dE 24.0, both above floor.

Usage:
  python 37_figure3_etiology.py
  python 37_figure3_etiology.py --event-type SE --out F3_etiology
  python 37_figure3_etiology.py --dcm-vs-hcm SE_DCMvsHCM   # if that table exists
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

BLUE = "#2a78d6"     # disease vs non-failing control
ORANGE = "#eb6834"   # disease vs other etiology
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#8a8a85"
GRID = "#e4e4e0"
SURF = "#fcfcfb"

NKEEP = {"SE": 3, "MXE": 5, "A3SS": 5, "A5SS": 5, "RI": 3}


def load_events(base, et, why=None):
    """why: optional list; on failure the exact reason is appended to it."""
    for p in (f"{base}_dedup.all_events.tsv", f"{base}.all_events.tsv"):
        if os.path.exists(p) and os.path.getsize(p):
            t = pd.read_csv(p, sep="\t")
            # Exact names first, then any dPSI* column: 07_dcm_vs_hcm_direct.py names
            # its effect column after the contrast (dPSI_DCMvsHCM), so a fixed list
            # will never cover every script in the project.
            d = next((c for c in ("dPSI_adj", "dPSI", "IncLevelDifference")
                      if c in t.columns), None)
            if d is None:
                cand = [c for c in t.columns if c.lower().startswith("dpsi")]
                d = cand[0] if len(cand) == 1 else None
                if len(cand) > 1:
                    d = next((c for c in cand if "adj" in c.lower()), cand[0])
            f = next((c for c in ("FDR", "fdr", "padj", "FDR_adj")
                      if c in t.columns), None)
            miss = []
            if d is None:
                miss.append("a dPSI column (dPSI_adj, dPSI, dPSI_*, "
                            "IncLevelDifference)")
            if f is None:
                miss.append("FDR / fdr / padj")
            if "key" not in t.columns:
                # No key means no exon-level collapse, and a full-coordinate count
                # next to exon-level ones would mix conventions inside one panel --
                # precisely what 24_lock_event_convention.py exists to prevent.
                miss.append("key (needed for the exon-level collapse)")
            if miss:
                if why is not None:
                    why.append(f"{p} is missing: {', '.join(miss)}\n"
                               f"         columns present: {', '.join(t.columns)}")
                return None, p
            if f != "FDR":
                t = t.rename(columns={f: "FDR"})
            if "se_adj" not in t.columns and "se" in t.columns:
                t = t.rename(columns={"se": "se_adj"})
            n = NKEEP[et]
            t["ekey"] = [":".join(str(k).split(":")[:n]) for k in t["key"]]
            t = t.loc[t.groupby("ekey")["FDR"].idxmin()].copy()
            return t.rename(columns={d: "dpsi"}), p
    if why is not None:
        why.append(f"no file at {base}_dedup.all_events.tsv or {base}.all_events.tsv")
    return None, f"{base}[_dedup].all_events.tsv"


def style(ax):
    ax.set_facecolor(SURF)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(GRID)
        ax.spines[s].set_linewidth(0.8)
    ax.tick_params(colors=INK2, labelsize=7, length=3, width=0.8)
    ax.title.set_color(INK)


# 162: the MAGNet arms re-base with the read length. Colorado does not: it was
# never re-quantified, which is why panel B must come out byte-identical.
F3SUF = os.environ.get("F3_SUFFIX", "")
# The 125 bp HCM-vs-NF table was RENAMED, not suffixed: SE_HCM_v2 became
# SE_HCMvsNF_rl125. Rule 8 -- a rename is only half an edit -- so it is named
# explicitly rather than derived.
F3_HCM = os.environ.get("F3_HCM_BASE", "")


def panel_a(ax, M, et, extra):
    """Counts: disease-vs-control against disease-vs-other-etiology."""
    sigcol = "sig_exon" if "sig_exon" in M.columns else "significant"
    tescol = "tested_exon" if "tested_exon" in M.columns else "tested"
    # 162: panel A reads the frozen master table, which has no 125 bp twin.
    # F3_PANEL_A overrides the plotted counts rather than inventing a frozen
    # file nothing else would maintain. Every named row must already exist, so
    # a typo refuses instead of silently plotting the 50 bp value.
    #   F3_PANEL_A="Colorado/DCMvsNF=476/19953,Colorado/DCMvsICM=2/19998,..."
    _ov = os.environ.get("F3_PANEL_A", "").strip()
    if _ov:
        M = M.copy()
        for _item in [x for x in _ov.split(",") if x.strip()]:
            _k, _v = _item.split("=")
            _coh, _con = _k.split("/")
            _s, _t = (int(x) for x in _v.split("/"))
            _m = ((M.cohort == _coh) & (M.contrast == _con)
                  & (M.event_type == et))
            if not _m.any():
                raise SystemExit(
                    f"REFUSING: F3_PANEL_A names {_k} for {et}, which is not a "
                    f"row of the master table. Nothing plotted.")
            M.loc[_m, sigcol] = _s
            M.loc[_m, tescol] = _t
        print(f"  [162] panel A overridden from F3_PANEL_A: {_ov}")
    rows = []

    def grab(cohort, contrast, colour, label):
        q = M[(M.cohort == cohort) & (M.contrast == contrast) & (M.event_type == et)]
        if q.empty:
            print(f"  [skip] {cohort} / {contrast} / {et} not in master_numbers")
            return
        r = q.iloc[0]
        if str(r.get("reliability", "")) == "exclude":
            print(f"  [skip] {cohort} / {contrast} / {et} is flagged EXCLUDE -- "
                  f"not plotted (this is correct)")
            return
        rows.append({"label": label, "n": int(r[sigcol]),
                     "tested": int(r[tescol]), "colour": colour})

    grab("Colorado", "DCM vs NF", BLUE, "Colorado\nDCM vs non-failing")
    grab("Colorado", "DCM vs ICM", ORANGE, "Colorado\nDCM vs ICM")
    grab("MAGNet", "DCM vs NF", BLUE, "MAGNet\nDCM vs non-failing")
    grab("MAGNet", "HCM vs NF", BLUE, "MAGNet\nHCM vs non-failing")
    if extra is not None:
        rows.append(extra)

    if not rows:
        ax.text(.5, .5, "no rows available", ha="center", va="center",
                transform=ax.transAxes, color=INK2, fontsize=8)
        return None
    R = pd.DataFrame(rows)[::-1].reset_index(drop=True)
    y = np.arange(len(R))
    ax.barh(y, R.n, height=0.62, color=R.colour, linewidth=0)
    for i, r in R.iterrows():
        off = max(R.n) * 0.015
        ax.text(r.n + off, i, f"{r.n:,}", va="center", ha="left",
                fontsize=7.5, color=INK, fontweight="bold")
        ax.text(r.n + off, i - 0.30, f"of {r.tested:,} tested", va="center",
                ha="left", fontsize=6, color=MUTED)
    ax.set_yticks(y)
    ax.set_yticklabels(R.label, fontsize=7, color=INK)
    ax.set_xlim(0, max(R.n) * 1.34)
    ax.set_xlabel(f"significant {et} events", fontsize=7.5, color=INK2)
    ax.xaxis.grid(True, color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)
    style(ax)
    ax.set_title("A   Etiology contrasts are null", loc="left",
                 fontsize=8.6, fontweight="bold", pad=8)
    return R


def panel_b(ax, et, a):
    """Power: could the etiology contrast have seen a disease-sized effect?"""
    E, pe = load_events(f"{et}_Colorado_DCMvsICM", et)
    D, pd_ = load_events(f"{et}_Colorado_DCMvsNF", et)
    if E is None or D is None or "mde80" not in (E.columns if E is not None else []):
        ax.text(.5, .5, "mde80 not available\n(needs the phase-2 tables)",
                ha="center", va="center", transform=ax.transAxes,
                color=INK2, fontsize=8)
        style(ax)
        return
    sig = D[(D.FDR < 0.05) & (D.dpsi.abs() >= 0.10)]
    m = E[E.ekey.isin(set(sig.ekey))]
    v = np.sort(m.mde80.dropna().values)
    if len(v) < 20:
        ax.text(.5, .5, f"only {len(v)} matched events", ha="center", va="center",
                transform=ax.transAxes, color=INK2, fontsize=8)
        style(ax)
        return
    ax.plot(v, np.arange(1, len(v) + 1) / len(v) * 100, color=ORANGE, linewidth=2,
            solid_capstyle="round")
    med = float(np.median(sig.dpsi.abs()))
    # vlines, not axvline, so the reference line stops at 100% and leaves a clear band
    # above the data for the annotation. Hunting for an empty quadrant fails as soon
    # as the data changes shape; reserved headroom does not.
    ax.vlines(med, 0, 100, color=BLUE, linewidth=1.4, linestyle=(0, (4, 3)))
    pct = 100 * (v <= med).mean()
    ax.plot([med], [pct], "o", ms=4.5, color=ORANGE, zorder=4,
            markeredgecolor=SURF, markeredgewidth=1.2)
    ax.text(0.02, 0.99,
            f"{pct:.0f}% could detect the median\nDCM-vs-NF effect (|ΔPSI| = {med:.2f})",
            transform=ax.transAxes, fontsize=6.4, color=INK, va="top",
            ha="left", linespacing=1.35)
    ax.set_xlabel("detectable |ΔPSI| in DCM vs ICM", fontsize=7.5, color=INK2)
    ax.set_ylabel("cumulative % of events", fontsize=7.5, color=INK2)
    ax.set_ylim(0, 132)          # 100-132 is reserved headroom for the annotation
    ax.set_yticks([0, 25, 50, 75, 100])
    ax.set_xlim(0, max(np.percentile(v, 99), med * 1.7))
    ax.grid(True, color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)
    style(ax)
    ax.set_title("B   The null is not low power", loc="left",
                 fontsize=8.6, fontweight="bold", pad=8)


def deming(x, y, sx, sy):
    lam = np.mean(sy ** 2) / np.mean(sx ** 2) if np.mean(sx ** 2) > 0 else 1.0
    dx, dy = x - x.mean(), y - y.mean()
    Sxx, Syy, Sxy = (dx ** 2).sum(), (dy ** 2).sum(), (dx * dy).sum()
    if Sxy == 0:
        return np.nan
    return float((Syy - lam * Sxx + np.sqrt((Syy - lam * Sxx) ** 2
                                            + 4 * lam * Sxy ** 2)) / (2 * Sxy))


def panel_c(ax, et, a):
    """MAGNet DCM-vs-NF against HCM-vs-NF: the two diseases share the program."""
    D, _ = load_events(f"{et}_v2filtered" + F3SUF, et)
    H, ph = load_events(F3_HCM or f"{et}_HCM_v2", et)
    if D is None or H is None:
        ax.text(.5, .5, f"missing {ph}", ha="center", va="center",
                transform=ax.transAxes, color=INK2, fontsize=7)
        style(ax)
        return
    sig = D[(D.FDR < 0.05) & (D.dpsi.abs() >= 0.10)][["ekey", "dpsi", "se_adj"]]
    m = sig.rename(columns={"dpsi": "x", "se_adj": "sx"}).merge(
        H[["ekey", "dpsi", "se_adj"]].rename(columns={"dpsi": "y", "se_adj": "sy"}),
        on="ekey")
    if len(m) < 20:
        ax.text(.5, .5, f"only {len(m)} shared events", ha="center", va="center",
                transform=ax.transAxes, color=INK2, fontsize=8)
        style(ax)
        return
    lim = float(np.nanpercentile(np.abs(np.r_[m.x, m.y]), 99.5)) * 1.1
    ax.axhline(0, color=GRID, linewidth=0.8)
    ax.axvline(0, color=GRID, linewidth=0.8)
    # Events enter only if |dPSI| >= 0.10 in the DCM contrast, which leaves an empty
    # band on the x-axis. Mark it, or a reader sees a hole in the data rather than a
    # selection boundary -- and that is the first thing a reviewer will ask about.
    for xv in (-0.10, 0.10):
        ax.axvline(xv, color=MUTED, linewidth=0.7, linestyle=(0, (1, 2)), zorder=1)
    ax.plot([-lim, lim], [-lim, lim], color=MUTED, linewidth=1,
            linestyle=(0, (4, 3)), zorder=1)
    ax.scatter(m.x, m.y, s=7, color=BLUE, alpha=0.45, linewidth=0, zorder=2)
    sl = deming(m.x.values, m.y.values, m.sx.values, m.sy.values)
    same = 100 * float((np.sign(m.x) == np.sign(m.y)).mean())
    if np.isfinite(sl):
        ax.plot([-lim, lim], [-lim * sl, lim * sl], color=ORANGE, linewidth=1.6,
                zorder=3)
    ax.text(0.04, 0.95,
            f"n = {len(m):,}\nDeming slope {sl:.2f}\n{same:.0f}% same direction",
            transform=ax.transAxes, fontsize=6.8, color=INK, va="top",
            linespacing=1.5, bbox=dict(facecolor="white", alpha=0.85, edgecolor="none", pad=2))
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    # equal aspect keeps the identity line at 45 degrees, which is the whole point
    # of the panel; anchor west so the shrunk box stays inside its gridspec cell
    # instead of drifting past the right margin and clipping the title.
    ax.set_aspect("equal", adjustable="box")
    ax.set_anchor("W")
    ax.set_xlabel("ΔPSI  DCM vs non-failing", fontsize=7.5, color=INK2)
    ax.set_ylabel("ΔPSI  HCM vs non-failing", fontsize=7.5, color=INK2)
    ax.grid(True, color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)
    style(ax)
    ax.set_title("C   Shared program", loc="left",
                 fontsize=8.6, fontweight="bold", pad=8)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--master", default="master_numbers_flagged.tsv")
    ap.add_argument("--event-type", default="SE")
    ap.add_argument("--dcm-vs-hcm", default=None,
                    help="prefix of the direct MAGNet DCM-vs-HCM table. Defaults to "
                         "'{ET}_DCMvsHCM', which is what 07_dcm_vs_hcm_direct.py "
                         "writes, so normally you do not need to pass this.")
    ap.add_argument("--no-dcm-vs-hcm", dest="want_hcm", action="store_false",
                    help="omit the direct DCM-vs-HCM bar even if the table exists")
    ap.add_argument("--out", default="F3_etiology")
    a = ap.parse_args()
    et = a.event_type

    if not os.path.exists(a.master):
        sys.exit(f"[fatal] {a.master} not found -- run 24 then 26 first")
    M = pd.read_csv(a.master, sep="\t")
    print(f"[info] panel A counts from {a.master} (frozen); "
          f"panels B and C from the event tables")

    extra = None
    if a.want_hcm:
        base = a.dcm_vs_hcm or f"{et}_DCMvsHCM"
        why = []
        T, p = load_events(base, et, why)
        if T is None:
            print(f"  [SKIP] the direct MAGNet DCM-vs-HCM bar is NOT in panel A.")
            print(f"         reason: {why[0] if why else 'unknown'}")
            print(f"         This is the paper's primary etiology comparison -- both "
                  f"arms are\n         explants, so it carries no donor-vs-explant "
                  f"confound. Fix the table or\n         pass --dcm-vs-hcm <prefix>.")
        else:
            n = int(((T.FDR < 0.05) & (T.dpsi.abs() >= 0.10)).sum())
            extra = {"label": "MAGNet\nDCM vs HCM", "n": n, "tested": len(T),
                     "colour": ORANGE}
            print(f"  [info] direct DCM vs HCM: {n} of {len(T):,} exon-level "
                  f"{et} events  ({p})")
            print(f"         NB this is {et} only. The 6-of-147,355 figure in the "
                  f"RESULT docs is\n         across all five event types -- do not "
                  f"conflate them in the legend.")

    plt.rcParams.update({"font.family": "sans-serif", "axes.linewidth": 0.8,
                         "figure.facecolor": SURF, "savefig.facecolor": SURF})
    fig = plt.figure(figsize=(7.4, 3.45))
    gs = gridspec.GridSpec(1, 3, width_ratios=[1.55, 1.0, 1.0], wspace=0.80,
                           left=0.155, right=0.965, top=0.845, bottom=0.275)
    R = panel_a(fig.add_subplot(gs[0]), M, et, extra)
    panel_b(fig.add_subplot(gs[1]), et, a)
    panel_c(fig.add_subplot(gs[2]), et, a)
    fig.text(0.155, 0.085,
             f"Significance: FDR < 0.05 and |ΔPSI| ≥ 0.10, exon-level, {et} events. "
             f"Group sizes are matched within each cohort.\nPanel C shows events "
             f"selected on the DCM contrast, so no point falls inside |ΔPSI| < 0.10 "
             f"on the x-axis (dotted guides).",
             fontsize=6.2, color=MUTED, ha="left", linespacing=1.6, va="top")

    for ext in ("pdf", "png"):
        fig.savefig(f"{a.out}.{ext}", dpi=400)
    print(f"\n[write] {a.out}.pdf and {a.out}.png")

    if R is not None and len(R) >= 2:
        pos = R[R.colour == BLUE].n.max()
        neg = R[R.colour == ORANGE].n
        if len(neg) and pos:
            print(f"\n  The sentence this figure makes: within one cohort and one "
                  f"pipeline,\n  at matched group sizes, the disease-vs-control "
                  f"contrast returns {pos:,}\n  significant {et} events and the "
                  f"etiology contrast returns {int(neg.min())}.")
    print("\n  NOW LOOK AT THE PNG. The validator checks colour, not layout --")
    print("  check for label collisions, clipped text and overlapping annotations")
    print("  before this goes anywhere near a manuscript.")


if __name__ == "__main__":
    main()
