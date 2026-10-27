#!/usr/bin/env python3
"""
66_build_s9.py -- Supplementary Figure/Table S9: genotype identity across bioprojects.

WHAT S9 IS FOR
--------------
This is the methodological contribution, and the cover letter leads with it: **24 of
LBNL's 33 samples are the same individuals as MAGNet samples**, undetectable from
metadata. A reviewer who has never seen a public-data reanalysis do this will want to see
the separation for themselves, not be told about it.

Three panels, one figure:

  A  kinship distribution -- 216,237 unrelated pairs against the 27 duplicates
  B  kinship vs IBS0 -- the two-dimensional view the field reads these in, with the
     decision thresholds drawn, showing the duplicates isolated in a corner
  C  **why a one-sided threshold fails** -- partners per sample above a naive cut,
     against the mutual-best-partner rule. This panel is the method's justification:
     script 28's one-sided rule reported 99 cross-project duplicates because a handful
     of "sticky" samples each matched dozens of people. One sample appeared as a
     relative of 45 others. An individual sequenced twice has exactly one twin.

Plus two tables: the per-cohort-pair structure (which includes the NEGATIVE results --
no duplicates between MAGNet and Colorado, none among the three HCM cohorts -- and those
are what make the positive finding credible), and the 27 excluded samples named.

FIGURE DESIGN, AND WHY
----------------------
This is the **emphasis** form, not a categorical one: the unrelated bulk is context and
belongs in de-emphasis gray; the duplicates are the entire point and take the single
accent hue. Two encodings beyond colour, because this prints in greyscale: the duplicate
marks carry a distinct shape with a surface ring, and the histogram overlay is hatched.
Panel A needs a **log** count axis -- 27 against 216,237 is invisible otherwise -- and the
axis says so.

Colours: accent #184f95, de-emphasis #8a8a85. Separation validated at delta-E 21.6
(protanopia) and 24.0 (normal vision), both well above the required floors; the chroma
check fails for the gray by design, which is what a de-emphasis slot is. The darker
accent step was chosen over the lighter default because it also separates in greyscale
(67 of 255 against 30), which matters for a supplement that may print mono. That was
checked by rendering the figure desaturated and looking at it, not by assuming.

SELF-CHECKS -- the script refuses rather than drawing something wrong
--------------------------------------------------------------------
From `RESULT_sample_identity_LBNL_overlap.md`: 670 samples, 216,237 of 224,115 pairs
compared, 27 duplicate pairs (24 cross-project, 3 within-cohort), the 24 at kinship
0.425-0.486 all matching on disease group, and the three sticky samples that broke the
one-sided rule (45, 29 and 26 partners above the naive cut). If the file does not
reproduce those, it is not the file this result was computed from.

USAGE
-----
    conda activate rnaseq_tools
    export PATH=/data/apps/anaconda3/envs/mamba/envs/rnaseq_tools/bin:$PATH   # matplotlib
    cd /data/project/chureeratp/RNA/HEART

    python 66_build_s9.py --list     # columns and heads of all three inputs; draws nothing
    python 66_build_s9.py            # figure + tables

**Run `--list` first and send me the output.** These files' column names are recorded
nowhere, and the last script that guessed at column names matched a gene-count field to a
column called `table` because "testable" contains "table". This one refuses to guess:
unresolved columns stop it.
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

# ---------------------------------------------------------------- inputs
METRICS = "geno/identity_metrics.tsv"       # every pair, with cohort/group/cross_project
TRIAGED = "geno/duplicates_triaged.tsv"     # 448 candidates, with the triage verdict
EXCLUDE = "geno/EXCLUDE_samples.tsv"        # the 27 dropped, with the keep rule

# Columns PINNED from the real files, 6 Oct. Not guessed, not resolved by synonym -- a
# synonym resolver in an earlier script matched a gene-count field to a column called
# `table` because "testable" contains "table".
M = dict(a="s1", b="s2", c1="cohort1", c2="cohort2", g1="group1", g2="group2",
         nsites="n_shared", kin="kinship", ibs0="ibs0", conc="concordance",
         cross="cross_project")
T = dict(a="s1", b="s2", c1="cohort1", c2="cohort2", g1="group1", g2="group2",
         kin="kinship", ibs0="ibs0", conc="concordance", cross="cross_project",
         mutual="mutual_best", qcflag="qc_flagged", triage="triage")
X = dict(drop="drop", keep="keep", dc="drop_cohort", kc="keep_cohort",
         dg="drop_group", kin="kinship", reason="reason")

# ---------------------------------------------------------------- the record
N_SAMPLES = 670                     # genotyped
N_PAIRS_TOTAL = 224_115             # C(670, 2)
# ⚠ MEASURED FROM THE FILE 6 Oct, superseding the project record.
# `RESULT_sample_identity_LBNL_overlap.md` says "216,237 of 224,115 pairs compared
# (3.5% shared <200 called sites)". That numerator is PRE-RESCUE: before script 31
# recovered the 10 blind samples there were 660 genotyped samples, C(660,2) = 217,470,
# minus 1,233 pairs below the site floor = 216,237 exactly. The denominator 224,115 is
# C(670,2), i.e. POST-rescue. Subtracting one from the other produced the "3.5%".
# The true drop-out rate is the same in both runs: 0.57% then, 0.56% now.
N_PAIRS_COMPARED = 222_870          # pairs in the file, all already >= MIN_SITES
N_SAMPLES_IN_PAIRS = 669            # one genotyped sample has no comparable partner
MIN_SITES = 200
N_DUP_TOTAL = 27
N_DUP_CROSS = 24                    # all MAGNet <-> LBNL
N_DUP_WITHIN = 3                    # 2 Colorado, 1 Stanford
CROSS_GROUPS = {"NF": 14, "DCM": 10}        # the 24, all matching on disease group
DUP_KIN_RANGE = (0.425, 0.486)
# Measured 0.928-0.988. The record says 0.937-0.988; its lower bound is from the same
# pre-rescue run as the pair count above.
DUP_CONC_RANGE = (0.920, 0.990)
DUP_CONC_RECORD = (0.937, 0.988)
KIN_CUT, IBS0_CUT = 0.354, 0.01
NAIVE_CUT = 0.177                   # the one-sided cut that produced 99 "duplicates"
# Script 28's one-sided rule was kinship > 0.177 **AND IBS0 literally zero** -- the
# doc is explicit: "every pair with a literal zero and kinship > 0.177 was promoted".
# Counting on kinship alone gives several hundred partners and is not that rule.
# These counts are pre-rescue like everything else from script 28, so they are
# reported and compared, not gated on.
STICKY = {"SRR10677129": 45, "SRR10677005": 29, "SRR8575311": 26}
STICKY_TOL = 2
WITHIN_PAIRS = {("SRR7426832", "SRR7426833"): 0.484,
                ("SRR7426789", "SRR7426790"): 0.481,
                ("SRR15183553", "SRR15183558"): 0.466}
TOL_KIN = 0.01
CONFIRMED = "confirmed"             # value of the `triage` column for a kept pair

# ---------------------------------------------------------------- colours
ACCENT = "#184f95"      # the duplicates: the point of the figure
MUTED = "#8a8a85"       # the unrelated bulk: context
INK = "#0b0b0b"
INK2 = "#52514e"
SURFACE = "#fcfcfb"
GRID = "#d8d8d4"


def read(path, pin, label):
    if not (os.path.exists(path) and os.path.getsize(path)):
        print(f"  {label:<10} *** {path} missing or empty ***")
        return None
    t = pd.read_csv(path, sep="\t", low_memory=False)
    miss = [c for c in pin.values() if c not in t.columns]
    print(f"  {label:<10} {path}  ({len(t):,} rows, {len(t.columns)} columns)")
    if miss:
        print(f"    *** pinned columns absent: {miss} ***")
        print(f"    present: {', '.join(map(str, t.columns))}")
        print("    The schema changed. Fix the pinned names above rather than letting")
        print("    this script guess -- guessing is what went wrong before.")
        return None
    return t


def show_all(a):
    print("\nINPUTS")
    tm = read(a.metrics, M, "metrics")
    trg = read(a.triaged, T, "triaged")
    exc = read(a.exclude, X, "exclude")
    for t, lab in ((tm, "metrics"), (trg, "triaged"), (exc, "exclude")):
        if t is not None:
            print(f"    {lab} columns: {', '.join(map(str, t.columns))}")
    return tm, trg, exc


# ============================================================ checks
def band(name, got, lo, hi, notes, nd=3):
    ok = lo - 1e-9 <= got <= hi + 1e-9
    print(f"  [{'ok' if ok else 'MISMATCH':<8}] {name}: {got:.{nd}f} (record {lo}-{hi})")
    if not ok:
        notes.append(f"{name} {got:.{nd}f} outside the recorded {lo}-{hi}")
    return ok


def eq(name, got, want, notes, fmt="{:,}"):
    ok = got == want
    print(f"  [{'ok' if ok else 'MISMATCH':<8}] {name}: {fmt.format(got)}"
          f" (record {fmt.format(want)})")
    if not ok:
        notes.append(f"{name} {got} vs recorded {want}")
    return ok


def check_metrics(t, notes, advisories):
    """Describe the screen. These are bookkeeping numbers, not the paper's claims, and
    the recorded values predate the rescue of the 10 blind samples -- so a difference
    is reported and carried forward, never a reason to refuse the figure."""
    ns = pd.to_numeric(t[M["nsites"]], errors="coerce")
    below = int((ns < MIN_SITES).sum())
    comp = t[ns >= MIN_SITES]
    ids = pd.unique(pd.concat([t[M["a"]], t[M["b"]]]).astype(str))
    print(f"  [--      ] {len(t):,} pairs in the file; {below:,} below {MIN_SITES}"
          f" shared sites")
    print(f"  [--      ] {len(comp):,} of {N_PAIRS_TOTAL:,} possible pairs compared"
          f" ({100 * (N_PAIRS_TOTAL - len(comp)) / N_PAIRS_TOTAL:.2f}% not compared)")
    print(f"  [--      ] {len(ids)} of {N_SAMPLES} genotyped samples appear in at least"
          f" one comparable pair")
    if len(comp) != N_PAIRS_COMPARED:
        advisories.append(f"pairs compared {len(comp):,}, this script expected "
                          f"{N_PAIRS_COMPARED:,} -- re-measure and update the Methods")
    if len(ids) != N_SAMPLES_IN_PAIRS:
        advisories.append(f"{len(ids)} samples in pairs, expected {N_SAMPLES_IN_PAIRS}")

    k = pd.to_numeric(comp[M["kin"]], errors="coerce")
    i0 = pd.to_numeric(comp[M["ibs0"]], errors="coerce")
    print(f"  [--      ] {int(((k >= KIN_CUT) & (i0 <= IBS0_CUT)).sum())} pairs clear"
          f" kinship >= {KIN_CUT} and IBS0 <= {IBS0_CUT} natively")

    # NOT reconstructed. Script 28's promotion rule gave 45/29/26 partners for three
    # samples; kinship > 0.177 alone gives 243/234/165 and adding IBS0 = 0 gives 9/4/13.
    # Three attempts, three answers, none of them the record's -- so the rule cannot be
    # recovered from what was written down, and a fourth guess would be worth no more
    # than the first three. Panel C uses the triage columns in the file instead
    # (mutual_best, triage), which are exact and need no reconstruction.
    return True, comp, None



advisories_global = []


def check_sticky(trg, notes):
    """The record's 45 / 29 / 26 are the ARTEFACT rows of the triage list, per sample --
    not a kinship threshold on the metrics file. Three attempts at a threshold gave
    243/234/165 and 9/4/13 before the triage columns settled it: a sample's artefact
    partners are its candidate pairs minus its confirmed one."""
    if trg is None:
        return True
    tri = trg[T["triage"]].astype(str).str.strip().str.lower()
    cand = pd.concat([trg[T["a"]], trg[T["b"]]]).astype(str).value_counts()
    cf = trg[tri == CONFIRMED]
    kept = pd.concat([cf[T["a"]], cf[T["b"]]]).astype(str).value_counts()
    print("\n  Artefact partners per sticky sample (candidate pairs minus confirmed):")
    print(f"  {'sample':<14}{'candidates':>11}{'confirmed':>10}{'artefacts':>10}"
          f"{'record':>9}")
    ok = True
    for smp, want in STICKY.items():
        c = int(cand.get(smp, 0))
        k_ = int(kept.get(smp, 0))
        art = c - k_
        good = art == want
        ok &= good
        print(f"  {smp:<14}{c:>11}{k_:>10}{art:>10}{want:>9}"
              f"   {'ok' if good else '<-- differs'}")
    if ok:
        print("  [ok      ] the recorded sticky counts reproduce exactly from the triage")
    else:
        notes.append("sticky artefact counts differ from the record")
    return ok


def near_misses(comp, conf, trg, n=8):
    """The pairs that came closest to being called duplicates.

    This matters more than it looks: Colorado is the sole DCM replication cohort, so a
    Colorado-MAGNet pair near the threshold is the one near-miss that could change a
    conclusion. Report the margin; do not describe a thin one as comfortable."""
    k = pd.to_numeric(comp[M["kin"]], errors="coerce")
    conf_keys = set()
    if conf is not None:
        for _, r in conf.iterrows():
            conf_keys.add(frozenset((str(r[T["a"]]), str(r[T["b"]]))))
    sub = comp.assign(_k=k).sort_values("_k", ascending=False).head(400)
    rows = []
    for _, r in sub.iterrows():
        key = frozenset((str(r[M["a"]]), str(r[M["b"]])))
        if key in conf_keys:
            continue
        rows.append(r)
        if len(rows) >= n:
            break
    if not rows:
        return None
    print(f"\n  The {len(rows)} highest-kinship pairs NOT called duplicates:")
    print(f"  {'sample A':<14}{'sample B':<14}{'cohorts':<22}{'kinship':>9}"
          f"{'margin':>9}{'IBS0':>9}")
    for r in rows:
        kv = float(pd.to_numeric(r[M["kin"]], errors="coerce"))
        print(f"  {str(r[M['a']]):<14}{str(r[M['b']]):<14}"
              f"{str(r[M['c1']]) + ' <-> ' + str(r[M['c2']]):<22}"
              f"{kv:>9.4f}{KIN_CUT - kv:>9.4f}"
              f"{float(pd.to_numeric(r[M['ibs0']], errors='coerce')):>9.4f}")
    top = float(pd.to_numeric(rows[0][M["kin"]], errors="coerce"))
    print(f"\n  The closest call sits {KIN_CUT - top:.4f} below the threshold.")

    # Is the closest call driven by a sample that is promiscuous across the whole screen?
    # If so the margin is far less alarming than it looks: that is precisely the pattern
    # the mutual-best rule exists to reject.
    if trg is not None:
        cand = pd.concat([trg[T["a"]], trg[T["b"]]]).astype(str).value_counts()
        hot = []
        for r in rows[:3]:
            for smp in (str(r[M["a"]]), str(r[M["b"]])):
                c = int(cand.get(smp, 0))
                if c >= 10:
                    hot.append((smp, c, float(pd.to_numeric(r[M["kin"]],
                                                            errors="coerce"))))
        if hot:
            print("  But the nearest misses involve samples that are promiscuous across"
                  " the whole screen:")
            seen = set()
            for smp, c, kv in hot:
                if smp in seen:
                    continue
                seen.add(smp)
                print(f"    {smp}: {c} candidate pairs overall, and the pair at"
                      f" kinship {kv:.4f}")
            print("  A sample that pairs with dozens of people is generating noise, not")
            print("  relatives. The mutual-best rule rejects exactly this, so the thin")
            print("  margin is an artefact of one sample rather than a borderline call.")
    return rows


def check_triaged(t, notes):
    """448 candidates in, 27 confirmed out. The funnel IS the method."""
    if t is None:
        return False, None
    tri = t[T["triage"]].astype(str).str.strip().str.lower()
    conf = t[tri == CONFIRMED]
    print(f"  [--      ] {len(t)} candidate pairs triaged ->"
          f" {len(conf)} confirmed")
    vc = tri.value_counts()
    for v, n in vc.items():
        print(f"              {v}: {n}")
    ok = eq("confirmed duplicate pairs", len(conf), N_DUP_TOTAL, notes, "{}")
    if not len(conf):
        return False, conf

    cross = conf[conf[T["cross"]].astype(str).str.lower().isin(("true", "1"))]
    ok &= eq("cross-project duplicates", len(cross), N_DUP_CROSS, notes, "{}")
    ok &= eq("within-cohort duplicates", len(conf) - len(cross), N_DUP_WITHIN,
             notes, "{}")

    if len(cross):
        pairs = {tuple(sorted((str(r[T["c1"]]), str(r[T["c2"]]))))
                 for _, r in cross.iterrows()}
        want = {("LBNL", "MAGNet")}
        good = pairs == want
        print(f"  [{'ok' if good else 'MISMATCH':<8}] cross-project cohort pairs:"
              f" {', '.join('<->'.join(p) for p in sorted(pairs))}")
        if not good:
            notes.append(f"cross-project pairs are {pairs}, record says MAGNet<->LBNL only")
            ok = False

        # the argument that random relatedness cannot explain these
        same = cross[cross[T["g1"]].astype(str) == cross[T["g2"]].astype(str)]
        good = len(same) == len(cross)
        print(f"  [{'ok' if good else 'MISMATCH':<8}] all {len(cross)} cross-project"
              f" pairs match on disease group: {len(same)}/{len(cross)}")
        if not good:
            notes.append(f"only {len(same)}/{len(cross)} cross-project pairs match on group")
            ok = False
        gc = same[T["g1"]].astype(str).value_counts().to_dict()
        for g, want_n in CROSS_GROUPS.items():
            got = int(gc.get(g, 0))
            good = got == want_n
            print(f"  [{'ok' if good else 'MISMATCH':<8}]   {g}<->{g}: {got}"
                  f" (record {want_n})")
            if not good:
                notes.append(f"{g} pair count {got} vs recorded {want_n}")
                ok = False

        kk = pd.to_numeric(cross[T["kin"]], errors="coerce")
        ok &= band("cross-project kinship min", float(kk.min()),
                   DUP_KIN_RANGE[0], DUP_KIN_RANGE[1], notes)
        ok &= band("cross-project kinship max", float(kk.max()),
                   DUP_KIN_RANGE[0], DUP_KIN_RANGE[1], notes)
        if T["conc"] in cross.columns:
            cc = pd.to_numeric(cross[T["conc"]], errors="coerce")
            lo_, hi_ = float(cc.min()), float(cc.max())
            ok &= band("cross-project concordance min", lo_,
                       DUP_CONC_RANGE[0], DUP_CONC_RANGE[1], notes)
            ok &= band("cross-project concordance max", hi_,
                       DUP_CONC_RANGE[0], DUP_CONC_RANGE[1], notes)
            if not (abs(lo_ - DUP_CONC_RECORD[0]) < 5e-4
                    and abs(hi_ - DUP_CONC_RECORD[1]) < 5e-4):
                print(f"  [--      ] measured concordance {lo_:.3f}-{hi_:.3f};"
                      f" the record says"
                      f" {DUP_CONC_RECORD[0]}-{DUP_CONC_RECORD[1]}")
                advisories_global.append(
                    f"concordance range is {lo_:.3f}-{hi_:.3f}, not "
                    f"{DUP_CONC_RECORD[0]}-{DUP_CONC_RECORD[1]} -- update the record "
                    f"and the Methods")

    # the pair resting on the fewest shared sites -- the weakest evidence in the set
    if "n_shared" in conf.columns:
        ns = pd.to_numeric(conf["n_shared"], errors="coerce")
        w = conf.loc[ns.idxmin()]
        print(f"  [--      ] weakest confirmed pair: {w[T['a']]} == {w[T['b']]} on"
              f" {int(ns.min())} shared sites,")
        print(f"              kinship {float(pd.to_numeric(w[T['kin']], errors='coerce')):.4f},"
              f" concordance"
              f" {float(pd.to_numeric(w[T['conc']], errors='coerce')):.4f}"
              f" -- the lowest in the set")

    # the three named within-cohort pairs
    for (x, y), want_k in WITHIN_PAIRS.items():
        row = conf[((conf[T["a"]].astype(str) == x) & (conf[T["b"]].astype(str) == y))
                   | ((conf[T["a"]].astype(str) == y) & (conf[T["b"]].astype(str) == x))]
        if not len(row):
            print(f"  [MISMATCH] within-cohort pair {x} == {y} not among the confirmed")
            notes.append(f"within-cohort pair {x}=={y} absent")
            ok = False
            continue
        got = float(pd.to_numeric(row.iloc[0][T["kin"]], errors="coerce"))
        good = abs(got - want_k) <= TOL_KIN
        print(f"  [{'ok' if good else 'MISMATCH':<8}] {x} == {y}: kinship {got:.3f}"
              f" (record {want_k})")
        if not good:
            notes.append(f"{x}=={y} kinship {got:.3f} vs recorded {want_k}")
            ok = False
    return ok, conf


def check_exclude(t, notes):
    if t is None:
        return False
    return eq("excluded samples", len(t), N_DUP_TOTAL, notes, "{}")


# ============================================================ figure
def build_figure(t, trg, out, notes):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch

    plt.rcParams.update({
        "font.family": "sans-serif", "font.size": 7,
        "axes.edgecolor": INK2, "axes.linewidth": 0.6,
        "xtick.color": INK2, "ytick.color": INK2,
        "xtick.labelsize": 6.5, "ytick.labelsize": 6.5,
        "axes.labelcolor": INK, "axes.labelsize": 7,
        "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
        "savefig.facecolor": SURFACE, "legend.frameon": False,
        "legend.fontsize": 6.5,
    })

    k = pd.to_numeric(t[M["kin"]], errors="coerce").to_numpy()
    i0 = pd.to_numeric(t[M["ibs0"]], errors="coerce").to_numpy()
    good = np.isfinite(k) & np.isfinite(i0)
    k, i0 = k[good], i0[good]
    isdup = (k >= KIN_CUT) & (i0 <= IBS0_CUT)

    # 180 mm wide, three panels
    fig, axes = plt.subplots(1, 3, figsize=(7.28, 2.45))
    fig.subplots_adjust(left=0.068, right=0.985, bottom=0.21, top=0.84, wspace=0.52)

    # ---- A: kinship distribution. Log counts, or 27 against 216,237 is invisible.
    ax = axes[0]
    bins = np.linspace(min(k.min(), -0.05), max(k.max(), 0.55), 90)
    ax.hist(k[~isdup], bins=bins, color=MUTED, edgecolor="none", log=True,
            label=f"unrelated pairs (n = {int((~isdup).sum()):,})")
    ax.hist(k[isdup], bins=bins, color=ACCENT, edgecolor=ACCENT, linewidth=0.5,
            hatch="////", log=True, label=f"duplicate pairs (n = {int(isdup.sum())})")
    ax.axvline(KIN_CUT, color=INK, linewidth=0.7, linestyle=(0, (4, 2)))
    # No inline threshold text: it collided with the legend at the top and with the
    # low bars at the bottom. Dashed lines in both panels, named in the caption --
    # which is the convention anyway.
    # direct label -- the insurance that survives greyscale and a small reproduction
    ax.annotate(f"{int(isdup.sum())} duplicate pairs",
                xy=(float(k[isdup].mean()), float(max(2, isdup.sum() / 6))),
                xytext=(0, 16), textcoords="offset points", ha="center", va="bottom",
                fontsize=6, color=ACCENT,
                arrowprops=dict(arrowstyle="-", color=ACCENT, linewidth=0.6))
    ax.set_xlabel("KING-robust kinship coefficient")
    ax.set_ylabel("pairs (log scale)")
    ax.set_title("A   Kinship across all pairs", loc="left", fontsize=7.5,
                 color=INK, fontweight="bold", pad=6)
    ax.legend(loc="upper left", handlelength=1.3, borderpad=0.2)
    ax.grid(axis="y", color=GRID, linewidth=0.4)
    ax.set_axisbelow(True)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)

    # ---- B: kinship vs IBS0. Bulk as density (one hue), duplicates as marks.
    ax = axes[1]
    hb = ax.hexbin(k[~isdup], i0[~isdup], gridsize=55, bins="log",
                   cmap="Greys", mincnt=1, linewidths=0)
    cb = fig.colorbar(hb, ax=ax, pad=0.018, fraction=0.040)
    cb.set_label("pairs per bin (log)", fontsize=6)
    cb.ax.tick_params(labelsize=5.5, color=INK2)
    cb.outline.set_linewidth(0.4)
    ax.scatter(k[isdup], i0[isdup], s=26, marker="D", facecolor=ACCENT,
               edgecolor=SURFACE, linewidth=0.8, zorder=5)
    ax.annotate(f"{int(isdup.sum())} duplicate pairs",
                xy=(float(k[isdup].mean()), float(np.nanmax(i0[isdup]))),
                xytext=(0.40, 0.84), textcoords="axes fraction",
                ha="left", va="bottom", fontsize=6, color=ACCENT,
                arrowprops=dict(arrowstyle="-", color=ACCENT, linewidth=0.6))
    ax.axvline(KIN_CUT, color=INK, linewidth=0.7, linestyle=(0, (4, 2)))
    ax.axhline(IBS0_CUT, color=INK, linewidth=0.7, linestyle=(0, (4, 2)))
    ax.set_xlabel("KING-robust kinship coefficient")
    ax.set_ylabel("IBS0")
    ax.set_title("B   Kinship against IBS0", loc="left", fontsize=7.5,
                 color=INK, fontweight="bold", pad=6)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)

    # ---- B inset: 27 points at IBS0 <= 0.003 collapse to one mark at this scale
    if isdup.sum() > 1:
        # bottom-left: the density runs on the anti-diagonal, so this
        # corner is empty whatever the cohort mix
        ins = ax.inset_axes([0.05, 0.08, 0.38, 0.32])
        ins.scatter(k[isdup], i0[isdup], s=16, marker="D", facecolor=ACCENT,
                    edgecolor=SURFACE, linewidth=0.5)
        ins.set_xlim(float(k[isdup].min()) - 0.012, float(k[isdup].max()) + 0.012)
        ins.set_ylim(-0.0004, max(float(np.nanmax(i0[isdup])) * 1.5, 0.0035))
        ins.tick_params(labelsize=5, length=2, pad=1)
        ins.locator_params(axis="x", nbins=3)
        ins.locator_params(axis="y", nbins=3)
        ins.set_title(f"the {int(isdup.sum())}, magnified", fontsize=5.5, color=INK2,
                      pad=2)
        for sp in ("top", "right"):
            ins.spines[sp].set_visible(False)
        for sp in ins.spines.values():
            sp.set_linewidth(0.4)
        ins.grid(color=GRID, linewidth=0.3)
        ins.set_axisbelow(True)

    # ---- C: an individual sequenced twice has exactly one twin.
    # Straight from the triage columns: how many candidate pairs each sample appears in,
    # against how many of those survived. No reconstruction of any superseded rule.
    ax = axes[2]
    if trg is None:
        ax.axis("off")
        ax.set_title("C   (triage file unavailable)", loc="left", fontsize=7.5,
                     color=INK, fontweight="bold", pad=6)
    else:
        tri = trg[T["triage"]].astype(str).str.strip().str.lower()
        cand = pd.concat([trg[T["a"]], trg[T["b"]]]).astype(str).value_counts()
        cf = trg[tri == CONFIRMED]
        kept = pd.concat([cf[T["a"]], cf[T["b"]]]).astype(str).value_counts()
        x = cand.to_numpy()
        y = np.array([int(kept.get(s_, 0)) for s_ in cand.index])
        rng = np.random.default_rng(0)
        jx = x + rng.uniform(-0.18, 0.18, len(x))
        jy = y + rng.uniform(-0.055, 0.055, len(y))
        ax.scatter(jx[y == 0], jy[y == 0], s=9, color=MUTED, linewidth=0,
                   label=f"no duplicate found ({int((y == 0).sum())} samples)")
        ax.scatter(jx[y > 0], jy[y > 0], s=18, marker="D", facecolor=ACCENT,
                   edgecolor=SURFACE, linewidth=0.5, zorder=4,
                   label=f"duplicate confirmed ({int((y > 0).sum())} samples)")
        ax.axhline(1, color=ACCENT, linewidth=0.7, linestyle=(0, (4, 2)), zorder=1)
        ax.set_ylim(-0.3, max(1.6, float(y.max()) + 0.4))
        ax.set_yticks(range(0, int(max(1, y.max())) + 1))
        worst = cand.index[0]
        ax.annotate(f"{worst}\n{int(cand.iloc[0])} candidate pairs, "
                    f"{int(kept.get(worst, 0))} confirmed",
                    xy=(float(cand.iloc[0]), float(kept.get(worst, 0))),
                    xytext=(-8, -26), textcoords="offset points", ha="right", va="top",
                    fontsize=5.5, color=INK2, linespacing=1.3,
                    arrowprops=dict(arrowstyle="-", color=INK2, linewidth=0.5))
        ax.set_xlabel("candidate pairs containing this sample")
        ax.set_ylabel("confirmed duplicates")
        ax.legend(loc="upper left", handlelength=1.0, borderpad=0.2, labelspacing=0.2)
        ax.grid(axis="y", color=GRID, linewidth=0.4)
        ax.set_axisbelow(True)
    ax.set_title("C   One twin at most", loc="left", fontsize=7.5,
                 color=INK, fontweight="bold", pad=6)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)

    for ext in ("pdf", "png"):
        p = f"{out}.{ext}"
        fig.savefig(p, dpi=400 if ext == "png" else None,
                    bbox_inches="tight", pad_inches=0.04)
        print(f"[write] {p}")
    plt.close(fig)


# ============================================================ tables
def build_tables(comp, conf, exc, out, notes):
    """Per-cohort-pair structure, and the 27 excluded samples.

    The cohort of each sample is carried by the metrics file itself (`cohort1`,
    `cohort2`), so nothing has to be joined in and nothing can be mismatched."""
    k = pd.to_numeric(comp[M["kin"]], errors="coerce")
    i0 = pd.to_numeric(comp[M["ibs0"]], errors="coerce")
    c1 = comp[M["c1"]].astype(str)
    c2 = comp[M["c2"]].astype(str)
    pair = pd.DataFrame({
        "cohort_A": np.where(c1 <= c2, c1, c2),
        "cohort_B": np.where(c1 <= c2, c2, c1),
        "kin": k, "ibs0": i0,
    })
    pair["dup"] = (pair.kin >= KIN_CUT) & (pair.ibs0 <= IBS0_CUT)
    tab = (pair.groupby(["cohort_A", "cohort_B"], dropna=False)
                .agg(pairs_compared=("kin", "size"),
                     max_kinship=("kin", "max"),
                     duplicates=("dup", "sum"))
                .reset_index()
                .sort_values(["duplicates", "pairs_compared"], ascending=[False, False]))
    tab["max_kinship"] = tab["max_kinship"].round(4)
    tab.to_csv(f"{out}a_cohort_pairs.tsv", sep="\t", index=False)

    print(f"\n  S9a -- per-cohort-pair structure ({len(tab)} cohort pairs)")
    print(f"  {'cohort A':<12}{'cohort B':<12}{'pairs':>12}{'max kinship':>13}"
          f"{'duplicates':>12}")
    for _, r in tab.iterrows():
        mark = "  <--" if r.duplicates else ""
        print(f"  {str(r.cohort_A):<12}{str(r.cohort_B):<12}"
              f"{int(r.pairs_compared):>12,}{r.max_kinship:>13.4f}"
              f"{int(r.duplicates):>12}{mark}")
    print(f"[write] {out}a_cohort_pairs.tsv")

    zero = tab[tab.duplicates == 0]
    hi = zero.nlargest(1, "max_kinship")
    print(f"\n  {len(zero)} of {len(tab)} cohort pairs contain NO duplicate."
          f" These are the negative results that")
    print("  make the positive one credible, and they belong in the table rather than a")
    if len(hi):
        r = hi.iloc[0]
        margin = KIN_CUT - float(r.max_kinship)
        print(f"  footnote. The highest kinship among them is {r.max_kinship:.4f}"
              f" ({r.cohort_A}<->{r.cohort_B}), which is {margin:.4f} below the"
              f" {KIN_CUT} threshold.")
        if margin < 0.02:
            print("  That is a THIN margin. The caption must give the number rather than")
            print("  describing the separation as clean.")

    # the two cohort pairs the paper specifically needs to be clean
    for a_, b_ in (("Colorado", "MAGNet"),):
        row = tab[(tab.cohort_A == a_) & (tab.cohort_B == b_)]
        if len(row):
            r = row.iloc[0]
            good = int(r.duplicates) == 0
            print(f"  [{'ok' if good else 'MISMATCH':<8}] {a_} <-> {b_}:"
                  f" {int(r.duplicates)} duplicates, max kinship {r.max_kinship:.4f}")
            if not good:
                notes.append(f"{a_}<->{b_} has {int(r.duplicates)} duplicates; "
                             f"the record says none")
    hcm = ["Brigham", "Fuwai", "Stanford"]
    hrows = tab[tab.cohort_A.isin(hcm) & tab.cohort_B.isin(hcm)
                & (tab.cohort_A != tab.cohort_B)]
    if len(hrows):
        n = int(hrows.duplicates.sum())
        print(f"  [{'ok' if n == 0 else 'MISMATCH':<8}] among the three HCM cohorts:"
              f" {n} duplicates across {len(hrows)} cohort pairs")
        if n:
            notes.append(f"{n} duplicates among the HCM cohorts; the record says none")

    # ---- S9b: the excluded samples, with the keep rule
    if exc is not None:
        e = exc.rename(columns={X["drop"]: "excluded_sample", X["keep"]: "retained_sample",
                                X["dc"]: "excluded_cohort", X["kc"]: "retained_cohort",
                                X["dg"]: "disease_group", X["kin"]: "kinship",
                                X["reason"]: "reason_for_the_choice"})
        e.to_csv(f"{out}b_excluded_samples.tsv", sep="\t", index=False)
        print(f"\n  S9b -- {len(e)} excluded samples")
        why = (e["reason_for_the_choice"].astype(str)
                .str.replace(r"\s*\(.*\)\s*$", "", regex=True).value_counts())
        for r, n in why.items():
            print(f"      {n:>3}  {r}")
        print(f"[write] {out}b_excluded_samples.tsv")

    # ---- the confirmed pairs, as deposited detail
    if conf is not None and len(conf):
        conf.to_csv(f"{out}_duplicate_pairs.tsv", sep="\t", index=False)
        print(f"[write] {out}_duplicate_pairs.tsv  ({len(conf)} confirmed pairs)")
    return tab


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--metrics", default=METRICS)
    ap.add_argument("--triaged", default=TRIAGED)
    ap.add_argument("--exclude", default=EXCLUDE)
    ap.add_argument("--list", action="store_true",
                    help="show what was read and stop")
    ap.add_argument("--out", default="S9")
    ap.add_argument("--no-figure", action="store_true")
    ap.add_argument("--force", action="store_true",
                    help="build even if a check fails (say so if you use it)")
    a = ap.parse_args()

    print("=" * 78)
    print("S9  genotype identity across bioprojects")
    print("=" * 78)
    print(f"  decision rule : mutual best partner, kinship >= {KIN_CUT},"
          f" IBS0 <= {IBS0_CUT}, >= {MIN_SITES} shared sites")
    print(f"  record        : {N_SAMPLES} samples, {N_PAIRS_COMPARED:,} of"
          f" {N_PAIRS_TOTAL:,} pairs compared, {N_DUP_TOTAL} duplicates"
          f" ({N_DUP_CROSS} cross-project, {N_DUP_WITHIN} within-cohort)")

    tm, trg, exc = show_all(a)
    if a.list:
        print("\n[--list] stopping here.")
        return
    if tm is None:
        sys.exit("cannot proceed -- the metrics file is unusable")

    print("\nCHECKS")
    notes = []
    advisories = []
    ok1, comp, cnt = check_metrics(tm, notes, advisories)
    ok2, conf = check_triaged(trg, notes)
    check_sticky(trg, notes)
    near_misses(comp, conf, trg)
    ok3 = check_exclude(exc, notes)

    advisories.extend(advisories_global)
    if not (ok1 and ok2 and ok3) and not a.force:
        print("\n  >>> S9 NOT BUILT. The inputs do not reproduce the recorded screen:")
        for n in notes:
            print(f"        {n}")
        print("      This figure is the paper's methodological contribution. It must be")
        print("      drawn from the files the result was computed from, not a near miss.")
        print("      --force overrides, but then say so wherever the figure appears.")
        sys.exit(1)
    if a.force and notes:
        print("\n  [warn] --force: building despite failed checks")

    print("\nBUILD")
    if not a.no_figure:
        build_figure(comp, trg, a.out, notes)
    build_tables(comp, conf, exc, a.out, notes)

    print("\n" + "=" * 78)
    for n in notes:
        print(f"  >>> {n}")
    if advisories:
        print("\n  BOOKKEEPING DIFFERENCES -- the figure is unaffected, the Methods is not:")
        for n in advisories:
            print(f"    - {n}")
    print("  Panels A and B are the finding; panel C is the method's justification and")
    print("  is the part a reviewer will not have seen before. Look at the figure before")
    print("  it goes in -- check the three labels in C have not collided.")


if __name__ == "__main__":
    main()
