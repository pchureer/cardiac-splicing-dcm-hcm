#!/usr/bin/env python3
"""
179_supp_fig2_matched_panel.py -- draw Supplementary Figure 2, which the
manuscript cites and nobody has ever drawn.

The matched-pair decomposition (platform versus people, 24 individuals in both
MAGNet and LBNL) exists as numbers in 36's terminal output and in 82's
S4_rl125.tsv, and as four event tables:

    mp_SE_selection{SUF}.all_events.tsv       selection arm (no shared sample)
    mp_SE_magnet_matched{SUF}.all_events.tsv  the 24 people, MAGNet
    mp_SE_lbnl_matched{SUF}.all_events.tsv    the same 24 people, LBNL
    mp_SE_magnet_other{SUF}.all_events.tsv    different people, MAGNet

Panel A  same people, two platforms   : lbnl_matched  against magnet_matched
Panel B  different people, one platform: magnet_other against magnet_matched

This script does not re-derive the estimator. It imports 82_build_s4.py's own
loader, Deming and lambda conventions, and SEARCHES -- selection rule x lambda
convention x axis orientation -- for the combination that reproduces the
recorded A, B, common-event count and direction counts. If exactly that
combination is found, it draws the point clouds with those slopes. If not, it
refuses: a figure whose slope does not match the text is worse than no figure.

Environment (the same three 82 demands, so copy them from 82's own refusal):
    TBL_SUFFIX=_rl125
    S4_RECORD=common=161,A=0.690,B=0.968,A_CI=0.570:0.806,B_CI=0.875:1.069,dirA=130,dirB=140,product=0.668,C=0.600
    S4_PER_COHORT=Colorado=0.5999,Brigham=0.5099,Stanford=0.4008,Fuwai=0.2948

    python 179_supp_fig2_matched_panel.py --out SuppFig2_matched_panel_rl125

Writes <out>.png (300 dpi), <out>.pdf and <out>_points.tsv (the plotted events,
so the figure is depositable). Nothing else is written or modified.

The intervals printed on the panels are the ones passed in S4_RECORD -- 36's,
which are what the main text quotes. 82's bootstrap differs by ~0.005 (different
interval method, same point estimates); do not mix them.
"""

import argparse
import importlib.util
import os
import sys
from pathlib import Path

import numpy as np

# --- 82's import-time refusal fires if TBL_SUFFIX is set without its record.
#     Say so here, before the traceback.
_SUF = os.environ.get("TBL_SUFFIX", "")
_REC = os.environ.get("S4_RECORD", "").strip()
if _SUF and not _REC:
    sys.exit("[REFUSE] TBL_SUFFIX is set but S4_RECORD is not. 82 refuses to import\n"
             "without it, and this script needs the same record to gate the slopes.\n"
             "See the docstring for the three variables.")


def parse_record(s):
    rec = {}
    for kv in s.split(","):
        if "=" not in kv:
            continue
        k, v = kv.split("=", 1)
        k, v = k.strip(), v.strip()
        if ":" in v:
            lo, hi = v.split(":")
            rec[k] = (float(lo), float(hi))
        else:
            rec[k] = float(v)
    for k in ("common", "A", "B", "dirA", "dirB"):
        if k not in rec:
            sys.exit(f"[REFUSE] S4_RECORD lacks {k}=")
    return rec


def import_82(path="82_build_s4.py"):
    p = Path(path)
    if not p.is_file():
        sys.exit(f"[fatal] {p} not found; run in the HEART directory")
    src = p.read_text()
    if '__name__ == "__main__"' not in src and "__name__ == '__main__'" not in src:
        sys.exit("[REFUSE] 82 has no __main__ guard; importing it would run the build.")
    for name in ("def load(", "def deming(", "def lam_val(", "SEL_RULES", "CONVENTIONS"):
        if name not in src:
            sys.exit(f"[REFUSE] 82 lacks {name}; this script expects the 26 Oct listing.")
    spec = importlib.util.spec_from_file_location("s4", str(p))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=f"SuppFig2_matched_panel{_SUF}")
    ap.add_argument("--tol", type=float, default=0.003,
                    help="slope tolerance against the record (default 0.003)")
    ap.add_argument("--dpi", type=int, default=300)
    a = ap.parse_args()

    if not _SUF:
        print("⚠ TBL_SUFFIX unset: drawing the ORIGINAL 50 bp panel. Pass "
              "TBL_SUFFIX=_rl125 for the submitted figure.")
        rec = parse_record(_REC or
                           "common=203,A=0.584,B=0.988,A_CI=0.478:0.693,"
                           "B_CI=0.889:1.089,dirA=0,dirB=0")
    else:
        rec = parse_record(_REC)

    m = import_82()
    names = {k: f"mp_SE_{k}{_SUF}.all_events.tsv"
             for k in ("selection", "magnet_matched", "lbnl_matched", "magnet_other")}
    missing = [v for v in names.values() if not os.path.isfile(v)]
    if missing:
        sys.exit("[REFUSE] missing arm(s):\n  " + "\n  ".join(missing))

    tabs = {}
    for k, path in names.items():
        t, dcol, scol = m.load(path)
        tabs[k] = t.set_index("ekey")
        print(f"  {k:15s} {len(t):7d} exon-level events   dPSI={dcol}  SE={scol}")

    sel, mm, lm, mo = (tabs[k] for k in
                       ("selection", "magnet_matched", "lbnl_matched", "magnet_other"))

    # --- search: selection rule x convention x orientation
    want_n = int(rec["common"])
    hits = []
    for rname, rule in m.SEL_RULES.items():
        keep = sel.index[rule(sel.reset_index()).to_numpy()]
        common = keep.intersection(mm.index).intersection(lm.index).intersection(mo.index)
        common = [k for k in common
                  if np.isfinite(mm.at[k, "dPSI"]) and np.isfinite(lm.at[k, "dPSI"])
                  and np.isfinite(mo.at[k, "dPSI"])]
        n = len(common)
        print(f"\n  rule {rname!r}: {n} common events" + ("  <- matches record" if n == want_n else ""))
        if n != want_n:
            continue
        x0 = mm.loc[common]
        for conv in m.CONVENTIONS:
            for orient in ("matched on x", "matched on y"):
                def slope(other):
                    if orient == "matched on x":
                        x, y = x0, other.loc[common]
                    else:
                        x, y = other.loc[common], x0
                    lam = m.lam_val(conv, x["SE"].to_numpy(), y["SE"].to_numpy())
                    if not np.isfinite(lam):
                        return np.nan
                    return m.deming(x["dPSI"].to_numpy(), y["dPSI"].to_numpy(), lam)
                A, B = slope(lm), slope(mo)
                if np.isfinite(A) and np.isfinite(B) and \
                        abs(A - rec["A"]) <= a.tol and abs(B - rec["B"]) <= a.tol:
                    hits.append((rname, conv, orient, A, B, common))
                    print(f"    ✓ {conv:32s} {orient:13s} A={A:.4f} B={B:.4f}")

    if not hits:
        sys.exit(f"\n[REFUSE] no rule x convention x orientation reproduces "
                 f"A={rec['A']} B={rec['B']} on {want_n} events within ±{a.tol}. "
                 f"Not drawing. Run 82 with the same environment and compare its "
                 f"search output; the record and the tables disagree.")

    # Prefer the convention 83/82 identified as the project's estimator.
    pref = [h for h in hits if h[1] == "mean of squared SEs [inverted]"
            and h[2] == "matched on x"] or hits
    rname, conv, orient, A, B, common = pref[0]
    print(f"\nusing rule={rname!r}  convention={conv!r}  {orient}")

    # --- direction counts, orientation independent; a second gate
    xA, yA = mm.loc[common, "dPSI"].to_numpy(), lm.loc[common, "dPSI"].to_numpy()
    xB, yB = mm.loc[common, "dPSI"].to_numpy(), mo.loc[common, "dPSI"].to_numpy()
    dA = int((np.sign(xA) == np.sign(yA)).sum())
    dB = int((np.sign(xB) == np.sign(yB)).sum())
    print(f"direction agreement: A {dA}/{len(common)}  B {dB}/{len(common)}  "
          f"(record {int(rec['dirA'])}/{int(rec['dirB'])})")
    if rec["dirA"] and (dA != int(rec["dirA"]) or dB != int(rec["dirB"])):
        sys.exit("[REFUSE] direction counts do not match the record; the common set "
                 "is not the recorded one even though the slopes agree. Not drawing.")

    # --- draw
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    if orient == "matched on y":   # plot with the matched MAGNet arm on x regardless
        print("⚠ the estimator put the matched arm on y; the panel still plots it on "
              "x and the drawn line uses the recorded slope's reciprocal geometry.")

    plt.rcParams.update({"font.size": 8, "axes.linewidth": 0.6,
                         "xtick.major.width": 0.6, "ytick.major.width": 0.6})
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.5), sharex=True, sharey=True)
    panels = [
        ("A", xA, yA, A, rec.get("A_CI"), "tab:orange",
         "Same individuals, two platforms",
         "ΔPSI, LBNL (24 matched individuals)"),
        ("B", xB, yB, B, rec.get("B_CI"), "tab:blue",
         "Different individuals, one platform",
         "ΔPSI, MAGNet (other individuals)"),
    ]
    lim = float(np.nanmax(np.abs(np.concatenate([xA, yA, xB, yB])))) * 1.08
    for ax, (lab, x, y, s, ci, col, title, ylab) in zip(axes, panels):
        ax.axhline(0, color="0.85", lw=0.6, zorder=0)
        ax.axvline(0, color="0.85", lw=0.6, zorder=0)
        ax.plot([-lim, lim], [-lim, lim], ls="--", color="0.6", lw=0.8, zorder=1,
                label="identity")
        ax.scatter(x, y, s=9, color="0.25", alpha=0.55, lw=0, zorder=2)
        mx, my = float(np.mean(x)), float(np.mean(y))
        if orient == "matched on x":
            xs = np.array([-lim, lim]); ys = my + s * (xs - mx)
        else:
            ys = np.array([-lim, lim]); xs = mx + s * (ys - my)
        ax.plot(xs, ys, color=col, lw=1.6, zorder=3)
        txt = f"Deming slope {s:.3f}"
        if ci:
            txt += f"\n95% CI {ci[0]:.3f}–{ci[1]:.3f}"
        txt += f"\nn = {len(x)}"
        ax.text(0.04, 0.96, txt, transform=ax.transAxes, va="top", ha="left",
                fontsize=7.5, color=col,
                bbox=dict(facecolor="white", alpha=0.85, edgecolor="none", pad=2))
        ax.set_title(f"{lab}  {title}", loc="left", fontsize=8.5, fontweight="bold")
        ax.set_xlabel("ΔPSI, MAGNet (24 matched individuals)")
        ax.set_ylabel(ylab)
        ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
        ax.set_aspect("equal")
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
    fig.tight_layout(w_pad=2.0)
    fig.savefig(a.out + ".png", dpi=a.dpi)
    fig.savefig(a.out + ".pdf")

    import pandas as pd
    pd.DataFrame({"ekey": common,
                  "dPSI_magnet_matched": xA, "dPSI_lbnl_matched": yA,
                  "dPSI_magnet_other": yB}).to_csv(a.out + "_points.tsv",
                                                   sep="\t", index=False)
    print(f"\nwrote {a.out}.png  {a.out}.pdf  {a.out}_points.tsv")
    print(f"A = {A:.4f}  B = {B:.4f}  product = {A*B:.3f}"
          + (f"  (record product {rec['product']}, C {rec.get('C')})" if "product" in rec else ""))


if __name__ == "__main__":
    main()
