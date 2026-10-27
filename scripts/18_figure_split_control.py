#!/usr/bin/env python3
"""
18_figure_split_control.py -- publication figure for the split-control result.

Replaces the throwaway plot in script 17, which was unreadable: a single unstable
MXE/program point at 2.357 (n=7) set the y-axis and flattened everything else into
a grey mat at the top.

Panel A  the measurement. Five random NF splits, shared vs disjoint controls.
         Every line drops, which IS the shared-control inflation.
Panel B  the ladder that makes the point for the paper:
              MAGNet HCM shared   ~0.98
              MAGNet HCM disjoint ~0.87   <- removing shared controls costs ~0.11
              Brigham / Fuwai / Stanford  0.44-0.70   <- a much bigger drop
         MAGNet HCM are explanted end-stage hearts; the other three are septal
         myectomy with preserved EF. Same disease, two stages.

Both panels plot the BACKGROUND and PROGRAM event sets side by side, because 15b
reported both and they must not be mixed: MAGNet program (0.99/0.96) is comparable
only with external program (0.587/0.440/0.703), and likewise for background.

SE only. MXE goes to a supplementary file; its program set (n=7-14) is too small
to plot and is dropped there with a note rather than shown as if it meant something.

Usage:
  python 18_figure_split_control.py                      # reads split_control_summary.tsv
  python 18_figure_split_control.py --tsv <path> --fmt pdf
"""
import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# categorical slots 1 and 2 of the validated default palette
C_BG, C_PROG = "#2a78d6", "#eb6834"
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#b8b7b2"

# 15b, HCM vs their OWN controls, against MAGNet DCM
EXTERNAL = {                    # cohort: (program slope, background slope, n_ctrl)
    "Brigham":  (0.587, 0.662, 5),
    "Fuwai":    (0.440, 0.328, 9),
    "Stanford": (0.703, 0.579, 7),
}


def style(ax):
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_color(MUTED)
    ax.tick_params(colors=INK2, labelsize=9, length=3, color=MUTED)
    ax.yaxis.grid(True, color=MUTED, lw=.6, alpha=.45)
    ax.set_axisbelow(True)


def panel_a(ax, R, et):
    sub = R[R.type == et]
    for st, col, lab in (("background", C_BG, "background"),
                         ("program", C_PROG, "RBM20 program")):
        s = sub[sub.set == st]
        if not len(s):
            continue
        for i, (_, r) in enumerate(s.iterrows()):
            ax.plot([0, 1], [r.deming_shared, r.deming_disjoint], "-",
                    color=col, lw=2, alpha=.55, solid_capstyle="round",
                    label=lab if i == 0 else None, zorder=2)
            ax.plot([0, 1], [r.deming_shared, r.deming_disjoint], "o",
                    color=col, ms=6, mec="white", mew=1.4, zorder=3)
        m0, m1 = s.deming_shared.mean(), s.deming_disjoint.mean()
        ax.plot([0, 1], [m0, m1], "-", color=col, lw=3.2, zorder=4,
                solid_capstyle="round")
        ax.annotate(f"Δ = −{m0 - m1:.3f}", xy=(1.04, m1),
                    color=col, fontsize=9, va="center", fontweight="bold")
    ax.set_xlim(-.22, 1.42)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["shared\ncontrols", "disjoint\ncontrols"], fontsize=9.5)
    ax.set_ylabel("Deming slope, HCM on DCM", fontsize=10, color=INK)
    ax.set_title(f"A   Shared-control inflation ({et}, 5 random NF splits)",
                 fontsize=11, color=INK, loc="left", pad=10)
    leg = ax.legend(frameon=False, fontsize=9, loc="lower left")
    for t in leg.get_texts():
        t.set_color(INK2)


def panel_b(ax, R, et):
    sub = R[R.type == et]
    labs = ["MAGNet\nshared", "MAGNet\ndisjoint",
            "Brigham", "Fuwai", "Stanford"]
    x = np.arange(5)
    for st, col, lab, k in (("background", C_BG, "background", 1),
                            ("program", C_PROG, "RBM20 program", 0)):
        s = sub[sub.set == st]
        if not len(s):
            continue
        vals = [s.deming_shared.mean(), s.deming_disjoint.mean()] + \
               [EXTERNAL[c][k] for c in ("Brigham", "Fuwai", "Stanford")]
        off = -.10 if st == "background" else .10
        ax.plot(x + off, vals, "o", color=col, ms=8, mec="white", mew=1.4,
                label=lab, zorder=3, ls="none")
        # spread across splits, as a vertical range on the two MAGNet points
        for xi, lo, hi in ((0, s.deming_shared.min(), s.deming_shared.max()),
                           (1, s.deming_disjoint.min(), s.deming_disjoint.max())):
            ax.plot([xi + off, xi + off], [lo, hi], "-", color=col, lw=2,
                    alpha=.55, solid_capstyle="round", zorder=2)

    ax.axvline(1.5, color=MUTED, lw=1, ls=":")
    ymin, ymax = ax.get_ylim()
    ax.set_ylim(ymin, ymax + .13 * (ymax - ymin))      # headroom for the labels
    ytop = ax.get_ylim()[1]
    ax.text(0.5, ytop, "end-stage\nexplant", ha="center", va="top",
            fontsize=9, color=INK2, style="italic", linespacing=1.3)
    ax.text(3.0, ytop, "septal myectomy,\npreserved EF", ha="center", va="top",
            fontsize=9, color=INK2, style="italic", linespacing=1.3)
    ax.set_xticks(x)
    ax.set_xticklabels(labs, fontsize=8.5)
    ax.set_xlim(-.6, 4.6)
    ax.set_ylabel("Deming slope, HCM on DCM", fontsize=10, color=INK)
    ax.set_title("B   Same disease, two stages", fontsize=11, color=INK,
                 loc="left", pad=10)
    leg = ax.legend(frameon=False, fontsize=9, loc="lower left")
    for t in leg.get_texts():
        t.set_color(INK2)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tsv", default="split_control_summary.tsv")
    ap.add_argument("--fmt", default="png", choices=["png", "pdf", "svg"])
    ap.add_argument("--out", default="fig_split_control")
    a = ap.parse_args()

    R = pd.read_csv(a.tsv, sep="\t")

    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.6))
    fig.patch.set_facecolor("white")
    panel_a(axes[0], R, "SE")
    panel_b(axes[1], R, "SE")
    for ax in axes:
        style(ax)
        ax.set_facecolor("white")
    fig.tight_layout(w_pad=3.0)
    fig.savefig(f"{a.out}.{a.fmt}", dpi=300, facecolor="white")
    print(f"[write] {a.out}.{a.fmt}")

    # ---- supplementary: MXE, background only -------------------------------
    mx = R[(R.type == "MXE") & (R.set == "background")]
    if len(mx):
        fig2, ax = plt.subplots(figsize=(4.6, 4.2))
        fig2.patch.set_facecolor("white")
        for i, (_, r) in enumerate(mx.iterrows()):
            ax.plot([0, 1], [r.deming_shared, r.deming_disjoint], "-o",
                    color=C_BG, lw=2, ms=6, alpha=.6, mec="white", mew=1.4)
        m0, m1 = mx.deming_shared.mean(), mx.deming_disjoint.mean()
        ax.plot([0, 1], [m0, m1], "-", color=C_BG, lw=3.2, zorder=4)
        ax.annotate(f"Δ = −{m0 - m1:.3f}", xy=(1.04, m1), color=C_BG,
                    fontsize=9, va="center", fontweight="bold")
        ax.set_xlim(-.22, 1.45)
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["shared\ncontrols", "disjoint\ncontrols"], fontsize=9.5)
        ax.set_ylabel("Deming slope, HCM on DCM", fontsize=10, color=INK)
        ax.set_title("MXE, background events", fontsize=11, color=INK, loc="left")
        style(ax)
        ax.set_facecolor("white")
        fig2.tight_layout()
        fig2.savefig(f"{a.out}_MXE_supp.{a.fmt}", dpi=300, facecolor="white")
        print(f"[write] {a.out}_MXE_supp.{a.fmt}")
        print("[note] the MXE RBM20-program set (n=7-14 events, one split returned a")
        print("       slope of 2.36) is deliberately NOT plotted. Report it as")
        print("       not estimable rather than showing it.")


if __name__ == "__main__":
    main()
