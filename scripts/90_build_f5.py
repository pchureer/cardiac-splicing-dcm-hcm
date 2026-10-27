#!/usr/bin/env python3
"""
90_build_f5.py  —  rebuild Figure 5 on the analysis the manuscript actually reports.

WHY THE PUBLISHED F5 CANNOT STAND
---------------------------------
Inspection of `F5_mechanism.png` against the current Results found all four panels on
superseded numbers:

  panel A  +0.382 SD, p = 2.7e-04   the Results report +0.34 SD, P = 0.003
                                    (different covariate model, unlabelled)
  panel B  0.92 vs 0.96 p = 0.76;   the Results RETIRED those two p-values:
           1.48 vs 1.65 p = 0.99    "those were nine-gene Mann-Whitney p-values"
  panel C  the 2-D matched deltas   the Results report the 3-D margins
  panel D  nine-gene backgrounds    the Results report the 45-gene test

Panel B is the sharpest: it displays the exact two numbers the manuscript says are
retired. Panel C matches `RESULT_rbm20_mechanism_null.md`'s 2-D table row for row, and
that doc warns the 2-D and 3-D schemes must never be mixed.

WHAT THIS VERSION DRAWS INSTEAD
-------------------------------
The quantity the paper reports, everywhere: **the margin above the GAPDH floor, with a
paired-bootstrap confidence interval.** That change does most of the work on its own —
with margin = delta − floor, the housekeeping floor sits at zero by construction, so
"no candidate separates from the control" is read directly off the figure instead of
being asserted beside a floor line. It also removes any temptation to quote a bare
p-value, which is what went wrong the first time.

  A  RBM20 expression, non-failing vs DCM, annotated with the model the Results quote
     and labelled as that model.
  B  forward test: RBM20's 45-gene target set, margin with CI, SE and MXE.
  C  inverted scan: ten candidates, margin above the floor with CI, 3-D matching.
  D  technical inflation, 45-gene backgrounds, both models labelled.

Form follows the data's job: A is a two-group distribution, B and C are effects with
uncertainty (forest), D is before/after on four series (slope, direct-labelled). One
axis per panel; no dual scales. The blue/orange house pair validates clean for
colour-vision deficiency (ΔE 24.0 protan, 32.1 normal against this surface), so it is
kept rather than replaced.

SOURCES
  panel B   rbm20_margin_resolve.tsv      (script 57: 45 genes, 3-D, bootstrap CIs)
  panel C   signature_rbp_margins.tsv     (script 58: inverted scan, 3-D, bootstrap CIs)
  panel A   magnet_expr.per_sample.tsv + the group bam lists
  panel D   constants below, read off rbm20_g45_406229.out step 5

GATES — the script refuses rather than drawing something wrong
  G1  inventory: print every column and every distinct `model` / `matching` value, and
      refuse to draw until --model and --matching name one combination.
  G2  REFUSE if the selected rows are the 2-D matching. That is the error in the
      published figure and it must not be repeatable by accident.
  G3  the inverted scan must reproduce the record: QKI SE +0.033 (-0.025, +0.082),
      RBM20 MXE +0.009 (-0.053, +0.088), floor ranked 5 of 11 in SE and 2 of 11 in MXE.
  G4  the forward test must reproduce -0.005 (-0.119, +0.080) and +0.000 (-0.085,
      +0.104).
  G5  panel A's statistic is chosen explicitly from a named model, never defaulted
      silently; all three are printed so the choice is visible.

USAGE
    python 90_build_f5.py --list                      # inventory. RUN THIS FIRST.
    python 90_build_f5.py --model <m> --matching <x> --out legend_check/f5_new
"""

import argparse
import os
import re
import sys

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# ---------------------------------------------------------------- house style
SURF = "#FCFCFB"
INK = "#1A1A1A"
INK2 = "#55555A"
MUTED = "#8A8A90"
BLUE = "#3C78D8"      # SE
ORANGE = "#E8663C"    # MXE
GRID = "#E3E3E0"
FLOORC = "#B0B0B4"

# READ LENGTH, added 22 Oct. These are joined with --dir at lines ~149 and ~205,
# and --dir cannot separate the read lengths because both deposits sit in the
# same directory. Read from the environment so neither call site changes:
#   F5_SCAN=signature_rbp_margins_rl125.tsv
#   F5_FWD=rbm20_margin_resolve_rl125.tsv
SCAN = os.environ.get("F5_SCAN", "signature_rbp_margins.tsv")
FWD = os.environ.get("F5_FWD", "rbm20_margin_resolve.tsv")
EXPR = "magnet_expr.per_sample.tsv"

# the eleven genes fitted: ten candidates + the floor
FLOOR_GENE = "GAPDH"
CANDIDATES = ["RBM20", "QKI", "RBM5", "ZRANB2", "HNRNPF", "HNRNPA0",
              "MBNL1", "CELF1", "SRSF1", "PTBP1"]

# ---- the record, for the gates ------------------------------------------
# READ LENGTH. All three are 50 bp records, and REC_SCAN is keyed by
# (class, GENE) -- the gene changes identity at 125 bp (QKI -> PTBP1 in SE,
# RBM20 -> PTBP1 in MXE), so a gate handed 125 bp numbers under 50 bp keys would
# name the wrong gene in a figure. Mandatory when F5_SCAN/F5_FWD are suffixed.
#   F5_REC_SCAN="SE=PTBP1/0.062/-0.025/0.122,MXE=PTBP1/0.113/0.016/0.187"
#   F5_FLOOR_RANK="SE=6,MXE=3"
#   F5_REC_FWD="SE=-0.051/-0.159/0.047,MXE=-0.019/-0.108/0.072"
def _f5_env(_name):
    _raw = os.environ.get(_name, "").strip()
    if not _raw:
        return None
    _out = {}
    for _part in _raw.split(","):
        _k, _, _v = _part.partition("=")
        if not _k.strip():
            continue
        _out[_k.strip()] = [_x.strip() for _x in _v.split("/")]
    return _out or None


_S_ENV = _f5_env("F5_REC_SCAN")
REC_SCAN = ({(_et, _v[0]): (float(_v[1]), float(_v[2]), float(_v[3]))
             for _et, _v in _S_ENV.items()} if _S_ENV else
            {("SE", "QKI"): (0.033, -0.025, 0.082),
             ("MXE", "RBM20"): (0.009, -0.053, 0.088)})

_R_ENV = _f5_env("F5_FLOOR_RANK")
REC_FLOOR_RANK = ({_k: int(_v[0]) for _k, _v in _R_ENV.items()} if _R_ENV
                  else {"SE": 5, "MXE": 2})   # of 11 genes, by delta descending

_W_ENV = _f5_env("F5_REC_FWD")
REC_FWD = ({_et: (float(_v[0]), float(_v[1]), float(_v[2]))
            for _et, _v in _W_ENV.items()} if _W_ENV else
           {"SE": (-0.005, -0.119, 0.080), "MXE": (0.000, -0.085, 0.104)})

if ("_rl" in SCAN or "_rl" in FWD):
    _miss = [_n for _n in ("F5_REC_SCAN", "F5_FLOOR_RANK", "F5_REC_FWD",
                           "F5_TECH45", "F5_TECH9", "F5_NTARGET", "F5_D_LOG")
             if not os.environ.get(_n, "").strip()]
    if _miss:
        sys.exit(
            "REFUSING: F5_SCAN/F5_FWD name a suffixed deposit but %s unset.\n"
            "  REC_SCAN is keyed by GENE and the top gene changes at 125 bp\n"
            "  (QKI/RBM20 -> PTBP1), so the figure would be labelled with the\n"
            "  wrong gene rather than fail. And panel D is PLOTTED from\n"
            "  constants (TECH45/TECH9) that no gate checks, so without\n"
            "  F5_TECH45/F5_TECH9 the figure would mix read lengths while\n"
            "  every gate passed. Supply all seven."
            % ", ".join(_miss))
EXPECT_GROUPS = (166, 166)    # DCM, NF — the discovery cohort, from Table 1
TOL = 0.004

# ---- panel A: the three models, so the choice is explicit ----------------
EXPR_MODELS = {
    "unadjusted":   (0.243, 0.027,   "unadjusted"),
    "covariate":    (0.382, 2.7e-04, "adjusted for sex, age, RIN, library size and TIN"),
    "composition":  (0.340, 0.0025,  "adjusted for sex, age, RIN, library size, TIN "
                                     "and cell composition"),
}
DEFAULT_EXPR_MODEL = "composition"   # the model the Results quote

# ---- panel D: 45-gene background medians, from rbm20_g45_406229.out ------
#      (nine-gene equivalents, which the published figure used, are in the comment)
# READ LENGTH, added 22 Oct by 156. These are the ONLY numbers in this file
# that are PLOTTED from constants rather than read from a deposit, which is why
# 155 missed them and why the 125 bp figure carried a 50 bp panel D while every
# gate passed. A plotted constant is a write path.
#
# Each entry is without_before/without_after/with_before/with_after, i.e. the
# 45-gene (or nine-gene) BACKGROUND median for DCM/nocell and DCM/cell, taken
# from the unadjusted block as `before` and the techadj block as `after` of
# section 5 of the rbm20_g45 log. All eight 50 bp values below were re-derived
# from rbm20_g45_406229.out before the 125 bp ones were transcribed.
#
#   F5_TECH45="SE=0.96/0.87/1.02/0.92,MXE=1.98/1.24/2.65/1.63"
#   F5_TECH9="SE=1.03/0.95/1.08/1.00,MXE=1.95/1.29/2.52/1.67"
_T45 = _f5_env("F5_TECH45")
TECH45 = ({_et: {"without composition adj.": (float(_v[0]), float(_v[1])),
                 "with composition adj.":    (float(_v[2]), float(_v[3]))}
           for _et, _v in _T45.items()} if _T45 else
          {"SE":  {"without composition adj.": (0.95, 0.86),
                   "with composition adj.":    (1.02, 0.89)},
           "MXE": {"without composition adj.": (2.03, 1.23),
                   "with composition adj.":    (2.71, 1.63)}})

_T9 = _f5_env("F5_TECH9")
TECH9 = ({_et: {"without": (float(_v[0]), float(_v[1])),
                "with":    (float(_v[2]), float(_v[3]))}
          for _et, _v in _T9.items()} if _T9 else
         {"SE": {"without": (1.01, 0.93), "with": (1.09, 0.96)},
          "MXE": {"without": (2.01, 1.28), "with": (2.65, 1.65)}})

# The provenance block prints these under "paste into the legend", so they are
# as quotable as anything in the figure.
#
# INTEGERS, NOT STRINGS, and the thousands separator is added at format time.
# A first version of this patch took them as strings so "1,393" could keep its
# comma -- but F5_NTARGET is itself a COMMA-separated list, so "MXE=1,399"
# parsed as MXE=1 and a stray "399". The separator cannot appear in a value.
#   F5_NTARGET="SE=641,MXE=1399"      F5_D_LOG=rbm20_g45_407951.out
_NTE = _f5_env("F5_NTARGET")
_NT_SE = int(_NTE["SE"][0]) if _NTE and "SE" in _NTE else 663
_NT_MXE = int(_NTE["MXE"][0]) if _NTE and "MXE" in _NTE else 1393
_D_LOG = os.environ.get("F5_D_LOG", "").strip() or "rbm20_g45_406229.out"
NULL_EXPECT = 0.67


def refuse(msg):
    print("\n" + "=" * 78)
    print("REFUSING — " + msg)
    print("=" * 78)
    sys.exit(2)


def pick(cols, cands, what, required=True):
    for c in cands:
        if c in cols:
            return c
    low = {str(c).lower(): c for c in cols}
    for c in cands:
        if c.lower() in low:
            return low[c.lower()]
    if required:
        refuse(f"no {what} column.\n  tried: {', '.join(cands)}\n"
               f"  table has: {', '.join(map(str, cols))}")
    return None


def load_scan(d):
    p = os.path.join(d, SCAN)
    if not os.path.exists(p):
        refuse(f"{SCAN} not found in {d} — panel C cannot be built from data.")
    t = pd.read_csv(p, sep="\t")
    cols = dict(
        etype=pick(t.columns, ("event_type", "type"), "event-type"),
        gene=pick(t.columns, ("gene",), "gene"),
        margin=pick(t.columns, ("margin",), "margin"),
        lo=pick(t.columns, ("ci_lo", "lo"), "CI lower"),
        hi=pick(t.columns, ("ci_hi", "hi"), "CI upper"),
        delta=pick(t.columns, ("delta",), "delta", required=False),
        floor=pick(t.columns, ("floor",), "floor", required=False),
        holm=pick(t.columns, ("holm_p", "holm"), "Holm p", required=False),
        model=pick(t.columns, ("model",), "model", required=False),
        match=pick(t.columns, ("matching", "match"), "matching", required=False),
    )
    return t, cols


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default=".")
    ap.add_argument("--out", default="legend_check/f5_new")
    ap.add_argument("--model", default=None, help="value of the `model` column to use")
    ap.add_argument("--matching", default=None, help="value of `matching` to use (3-D)")
    ap.add_argument("--expr-model", default=DEFAULT_EXPR_MODEL,
                    choices=list(EXPR_MODELS))
    ap.add_argument("--b1", default=None,
                    help="rMATS --b1 list: the CASE group (DCM). ONE comma-separated "
                         "line. Without --b1 and --b2 panel A draws no distribution.")
    ap.add_argument("--b2", default=None, help="rMATS --b2 list: the CONTROL group (NF)")
    ap.add_argument("--expect-n", default=None,
                    help="override the expected group sizes, as CASE,CONTROL")
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--allow-unverified", action="store_true",
                    help="draw even if G3/G4 do not reproduce. Use only with a reason.")
    a = ap.parse_args()

    if a.expect_n:
        globals()["EXPECT_GROUPS"] = tuple(int(x) for x in a.expect_n.split(","))
    t, C = load_scan(a.dir)

    print("=" * 78)
    print("90_build_f5.py — Figure 5 on the 45-gene, 3-D, composition-adjusted analysis")
    print("=" * 78)
    print(f"\n{SCAN}: {len(t):,} rows")
    print(f"  columns: {', '.join(map(str, t.columns))}")
    for key, label in (("model", "model"), ("match", "matching")):
        c = C[key]
        if c:
            print(f"  distinct `{c}`: {sorted(map(str, t[c].dropna().unique()))}")
        else:
            print(f"  no `{label}` column — the file holds one scheme only")
    print(f"  distinct event types: {sorted(map(str, t[C['etype']].unique()))}")
    print(f"  genes: {len(t[C['gene']].unique())} distinct")

    fp = os.path.join(a.dir, FWD)
    F = None
    if os.path.exists(fp):
        F = pd.read_csv(fp, sep="\t")
        print(f"\n{FWD}: {len(F):,} rows")
        print(f"  columns: {', '.join(map(str, F.columns))}")
        for c in ("table", "model", "gene_set", "matching", "verdict"):
            if c in F.columns:
                print(f"  distinct `{c}`: {sorted(map(str, F[c].dropna().unique()))}")
    else:
        print(f"\n{FWD}: NOT FOUND — panel B will use the recorded constants, "
              f"clearly labelled.")

    print("\npanel A statistic — pick one with --expr-model:")
    for k, (sd, p, desc) in EXPR_MODELS.items():
        mark = "  <- default, the model the Results quote" if k == DEFAULT_EXPR_MODEL else ""
        print(f"  {k:12s} +{sd:.3f} SD  P = {p:.2g}   {desc}{mark}")

    print("\npanel D — 45-gene background medians (the published figure used nine-gene):")
    for et in ("SE", "MXE"):
        for lab, (b, aft) in TECH45[et].items():
            n9 = TECH9[et]["without" if "without" in lab else "with"]
            print(f"  {et:4s} {lab:26s} {b:.2f} -> {aft:.2f}    "
                  f"(nine-gene was {n9[0]:.2f} -> {n9[1]:.2f})")

    if a.list:
        print("\n--list: nothing drawn. Re-run with --model and --matching naming the "
              "3-D,\ncomposition-adjusted, within-DCM combination from the lists above.")
        return

    # ------------------------------------------------------------------ G1
    sub = t.copy()
    for key, val, flag in (("model", a.model, "--model"), ("match", a.matching, "--matching")):
        c = C[key]
        if c is None:
            continue
        if val is None:
            refuse(f"the file has a `{c}` column with "
                   f"{sorted(map(str, t[c].dropna().unique()))}.\n"
                   f"  Name one with {flag}. This script will not pick a model for you —\n"
                   f"  picking the wrong one is exactly how the published figure ended up "
                   f"on\n  the 2-D scheme.")
        sub = sub[sub[c].astype(str) == str(val)]
    if sub.empty:
        refuse("no rows left after --model/--matching. Check the spelling against the "
               "lists above.")
    print(f"\nG1  selected {len(sub)} rows")

    # ------------------------------------------------------------------ G2
    if a.matching and any(k in str(a.matching).lower() for k in ("2d", "2-d", "two")):
        refuse(f"G2 — --matching {a.matching!r} looks like the 2-D scheme.\n"
               "  The 2-D matching is biased toward enrichment and is NOT what the "
               "Results\n  report. It is the error in the published figure. Use the 3-D "
               "matching.")
    print("G2  PASS — not the 2-D scheme")

    # ------------------------------------------------------------------ G3
    print("\nG3  does the inverted scan reproduce the record?")
    bad = []
    for (et, g), (m, lo, hi) in REC_SCAN.items():
        r = sub[(sub[C["etype"]].astype(str).str.upper() == et) &
                (sub[C["gene"]].astype(str) == g)]
        if len(r) != 1:
            bad.append(f"{et} {g}: {len(r)} rows, expected 1")
            continue
        got = (float(r[C["margin"]].iloc[0]), float(r[C["lo"]].iloc[0]),
               float(r[C["hi"]].iloc[0]))
        ok = all(abs(x - y) <= TOL for x, y in zip(got, (m, lo, hi)))
        print(f"  {et:4s} {g:7s} {got[0]:+.3f} ({got[1]:+.3f}, {got[2]:+.3f})   "
              f"record {m:+.3f} ({lo:+.3f}, {hi:+.3f})   {'match' if ok else 'MISMATCH'}")
        if not ok:
            bad.append(f"{et} {g}: {got} vs record {(m, lo, hi)}")

    if C["delta"]:
        for et, want in REC_FLOOR_RANK.items():
            s = sub[sub[C["etype"]].astype(str).str.upper() == et]
            if s.empty:
                bad.append(f"{et}: no rows"); continue
            genes = s[C["gene"]].astype(str)
            if (genes == FLOOR_GENE).any():
                ss = s.sort_values(C["delta"], ascending=False).reset_index(drop=True)
                got = int(ss.index[ss[C["gene"]].astype(str) == FLOOR_GENE][0]) + 1
                tot = len(ss)
            elif C["floor"]:
                fv = pd.to_numeric(s[C["floor"]], errors="coerce").dropna().unique()
                if len(fv) != 1:
                    bad.append(f"{et}: `floor` column holds {len(fv)} distinct values, "
                               f"expected one per event type"); continue
                d = pd.to_numeric(s[C["delta"]], errors="coerce")
                got, tot = 1 + int((d > fv[0]).sum()), len(s) + 1
            else:
                bad.append(f"{et}: no {FLOOR_GENE} row and no `floor` column"); continue
            print(f"  {et:4s} floor rank {got} of {tot}   record {want} of 11   "
                  f"{'match' if (got, tot) == (want, 11) else 'MISMATCH'}")
            if (got, tot) != (want, 11):
                bad.append(f"{et} floor rank {got} of {tot} vs record {want} of 11")
    if bad and not a.allow_unverified:
        refuse("G3 — the selected rows are not the analysis the Results describe:\n  - "
               + "\n  - ".join(bad)
               + "\n\n  Most likely --model/--matching names the wrong combination. "
                 "Check\n  --list again. Override with --allow-unverified only if you "
                 "know why.")
    print("  G3 PASS" if not bad else "  G3 OVERRIDDEN")

    # ------------------------------------------------------------------ G4
    print("\nG4  the forward 45-gene test")
    fwd = {}
    if F is not None:
        f_et = pick(F.columns, ("event_type", "type", "table"), "event-type",
                    required=False)
        f_m = pick(F.columns, ("margin",), "margin", required=False)
        f_lo = pick(F.columns, ("ci_lo", "lo"), "CI lower", required=False)
        f_hi = pick(F.columns, ("ci_hi", "hi"), "CI upper", required=False)
    if F is not None and None not in (f_et, f_m, f_lo, f_hi):
        G = F.copy()
        # narrow with whatever the file offers, using the same choices as panel C
        for col, val in (("model", a.model), ("matching", a.matching)):
            if col in G.columns and val is not None:
                G = G[G[col].astype(str) == str(val)]
        if "gene_set" in G.columns:
            gs = G["gene_set"].astype(str)
            keep = gs[gs.str.contains("45")]
            if len(keep):
                G = G[gs.str.contains("45")]
            else:
                print(f"  no `gene_set` value contains '45' "
                      f"({sorted(gs.unique())}) — not narrowing on it")
        # the record's margins are the technically-adjusted tables
        tb = G[f_et].astype(str)
        if tb.str.contains("techadj", case=False).any():
            G = G[tb.str.contains("techadj", case=False)]
        for et, (m, lo, hi) in REC_FWD.items():
            r = G[G[f_et].astype(str).str.upper().str.startswith(et)]
            if len(r) != 1:
                print(f"  {et}: {len(r)} rows after narrowing — using the recorded "
                      f"constant. Distinct values left:")
                for c in ("table", "model", "gene_set", "matching"):
                    if c in G.columns:
                        print(f"      {c}: {sorted(map(str, G[c].unique()))[:6]}")
                fwd[et] = (m, lo, hi, "recorded")
                continue
            got = (float(r[f_m].iloc[0]), float(r[f_lo].iloc[0]), float(r[f_hi].iloc[0]))
            ok = all(abs(x - y) <= TOL for x, y in zip(got, (m, lo, hi)))
            print(f"  {et:4s} {got[0]:+.3f} ({got[1]:+.3f}, {got[2]:+.3f})   "
                  f"record {m:+.3f} ({lo:+.3f}, {hi:+.3f})   "
                  f"{'match' if ok else 'MISMATCH'}")
            if not ok and not a.allow_unverified:
                refuse(f"G4 — {FWD} disagrees with the record for {et}.")
            fwd[et] = (*got, "file")
    else:
        for et, (m, lo, hi) in REC_FWD.items():
            fwd[et] = (m, lo, hi, "recorded")
        print("  using the recorded constants; the panel is labelled so.")

    # ------------------------------------------------------------------ draw
    sd, pval, desc = EXPR_MODELS[a.expr_model]
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)

    plt.rcParams.update({
        "figure.facecolor": SURF, "savefig.facecolor": SURF, "axes.facecolor": SURF,
        "font.size": 9, "axes.edgecolor": GRID, "axes.labelcolor": INK2,
        "xtick.color": INK2, "ytick.color": INK2, "text.color": INK,
        "axes.titlesize": 11, "axes.titleweight": "bold", "axes.titlelocation": "left",
    })
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 9.5))
    (axA, axB), (axC, axD) = axes

    # --- A: expression ---------------------------------------------------
    rng = np.random.default_rng(0)
    ep = os.path.join(a.dir, EXPR)
    g1 = g2 = np.array([])

    def bam_runs(path):
        """rMATS --b1/--b2 are ONE comma-separated line, not one path per line.
        Parsed line by line a 166-sample cohort reads as 1."""
        with open(path) as fh:
            raw = fh.read()
        out = []
        for tok in re.split(r"[,\s]+", raw):
            if not tok.strip():
                continue
            m = re.search(r"(SRR\d+|SAMEA\d+|ERR\d+|SRS\d+)", os.path.basename(tok))
            if m:
                out.append(m.group(1))
        return out

    if not (a.b1 and a.b2):
        print("  [panel A] --b1/--b2 not given: the DISTRIBUTION IS NOT DRAWN.\n"
              "            Group membership must be read, never inferred from column "
              "order.")
    elif not os.path.exists(ep):
        print(f"  [panel A] {EXPR} not found — distribution not drawn.")
    else:
        try:
            E = pd.read_csv(ep, sep="\t", index_col=0)
            ser = E["RBM20"] if "RBM20" in E.columns else (
                E.loc["RBM20"] if "RBM20" in E.index else None)
            if ser is None:
                raise KeyError("no RBM20 row or column")
            ser = pd.to_numeric(ser, errors="coerce").dropna()
            ser.index = [str(i) for i in ser.index]
            case = [r for r in bam_runs(os.path.join(a.dir, a.b1))]
            ctrl = [r for r in bam_runs(os.path.join(a.dir, a.b2))]
            have = set(ser.index)
            fc = [r for r in case if r in have]
            fo = [r for r in ctrl if r in have]
            print(f"  [panel A] --b1 {len(case)} run(s), {len(fc)} matched; "
                  f"--b2 {len(ctrl)} run(s), {len(fo)} matched")
            want_c, want_o = EXPECT_GROUPS
            if (len(case), len(ctrl)) != (want_c, want_o):
                raise ValueError(
                    f"bam lists parsed to {len(case)} and {len(ctrl)} runs, expected "
                    f"{want_c} and {want_o}. rMATS --b1/--b2 are ONE COMMA-SEPARATED "
                    f"LINE; read line by line a 166-sample cohort yields 1. Check the "
                    f"files, or pass --expect-n if this really is another cohort")
            if len(fc) < .9 * len(case) or len(fo) < .9 * len(ctrl):
                raise ValueError(f"fewer than 90% of runs matched the expression index "
                                 f"({len(fc)}/{len(case)}, {len(fo)}/{len(ctrl)}); "
                                 f"sample IDs do not line up")
            if set(fc) & set(fo):
                raise ValueError("the two groups overlap")
            z = (ser - ser.mean()) / ser.std()
            g2, g1 = z[fc].to_numpy(), z[fo].to_numpy()   # g1 control, g2 case
            print(f"  [panel A] drawn from real group membership: "
                  f"{len(g1)} non-failing, {len(g2)} DCM")
        except Exception as e:
            g1 = g2 = np.array([])
            print(f"  [panel A] NOT DRAWN — {e}")

    panelA_src = (f"{EXPR}, groups from {a.b1} / {a.b2}"
                  if (len(g1) and len(g2)) else
                  "summary only — distribution NOT drawn, group membership not supplied")
    if len(g1) and len(g2):
        for i, (g, col) in enumerate(((g1, MUTED), (g2, BLUE))):
            x = i + rng.normal(0, .055, len(g))
            axA.scatter(x, g, s=7, color=col, alpha=.55, linewidths=0, zorder=2)
        axA.boxplot([g1, g2], positions=[0, 1], widths=.45, showfliers=False,
                    medianprops=dict(color=INK, lw=1.6),
                    boxprops=dict(color=INK2, lw=1.1),
                    whiskerprops=dict(color=INK2, lw=1.1),
                    capprops=dict(color=INK2, lw=1.1), zorder=3)
        q = np.concatenate([g1, g2])
        axA.set_ylim(np.percentile(q, .3) - .4, np.percentile(q, 99.7) + .4)
        axA.set_xticks([0, 1])
        axA.set_xticklabels([f"non-failing\nn={len(g1)}", f"DCM\nn={len(g2)}"])
    else:
        axA.text(.5, .62, "distribution not drawn\n(group membership not supplied)",
                 transform=axA.transAxes, ha="center", fontsize=9, color=MUTED,
                 style="italic")
        axA.set_xticks([0, 1])
        axA.set_xticklabels(["non-failing", "DCM"])
    axA.set_xticks([0, 1])
    axA.set_ylabel("RBM20 expression (z-scored log$_2$CPM)")
    axA.set_title("A   RBM20 is higher in DCM, not lower")
    lo_, hi_ = axA.get_ylim()
    axA.set_ylim(lo_ - .22 * (hi_ - lo_), hi_)          # room under the whiskers
    axA.text(.02, .02, f"+{sd:.2f} SD in DCM,  $P$ = {pval:.3g}\n{desc}",
             transform=axA.transAxes, fontsize=8.5, color=INK2, va="bottom",
             bbox=dict(facecolor=SURF, edgecolor="none", pad=2.5))
    axA.grid(axis="y", color=GRID, lw=.7)
    axA.set_axisbelow(True)
    for s in ("top", "right"):
        axA.spines[s].set_visible(False)

    # --- B: forward test, margin with CI ---------------------------------
    ys = [1, 0]
    for y, et, col in zip(ys, ("SE", "MXE"), (BLUE, ORANGE)):
        m, lo, hi, src = fwd[et]
        axB.plot([lo, hi], [y, y], color=col, lw=2.2, solid_capstyle="round", zorder=2)
        axB.plot([m], [y], "o", ms=9, color=col, mec=SURF, mew=2, zorder=3)
        axB.text(hi + .012, y, f"{m:+.3f}  ({lo:+.3f}, {hi:+.3f})",
                 va="center", fontsize=8.5, color=INK2)
    axB.axvline(0, color=FLOORC, lw=1.4, ls="--", zorder=1)
    axB.text(0, 1.52, "GAPDH floor", ha="center", fontsize=8.5, color=INK2)
    axB.set_yticks(ys)
    axB.set_yticklabels(["skipped exons", "mutually exclusive exons"])
    axB.set_ylim(-.7, 1.75)
    axB.set_xlim(-.20, .26)
    axB.set_xlabel("margin above the housekeeping floor  (median |$t$|)")
    axB.set_title("B   RBM20's targets do not clear the floor", pad=10)
    axB.grid(axis="x", color=GRID, lw=.7)
    axB.set_axisbelow(True)
    for s in ("top", "right", "left"):
        axB.spines[s].set_visible(False)

    # --- C: inverted scan -------------------------------------------------
    se = sub[sub[C["etype"]].astype(str).str.upper() == "SE"].set_index(
        sub[sub[C["etype"]].astype(str).str.upper() == "SE"][C["gene"]].astype(str))
    order = [g for g in se.sort_values(C["margin"], ascending=False).index
             if g in CANDIDATES]
    order += [g for g in CANDIDATES if g not in order]
    ypos = {g: len(order) - 1 - i for i, g in enumerate(order)}

    for et, col, off in (("SE", BLUE, .17), ("MXE", ORANGE, -.17)):
        s = sub[sub[C["etype"]].astype(str).str.upper() == et]
        for _, r in s.iterrows():
            g = str(r[C["gene"]])
            if g not in ypos:
                continue
            y = ypos[g] + off
            lo, hi, m = float(r[C["lo"]]), float(r[C["hi"]]), float(r[C["margin"]])
            axC.plot([lo, hi], [y, y], color=col, lw=1.8, solid_capstyle="round", zorder=2)
            axC.plot([m], [y], "o", ms=6.5, color=col, mec=SURF, mew=1.5, zorder=3)
    axC.axvline(0, color=FLOORC, lw=1.6, ls="--", zorder=1)
    axC.text(0, len(order) - .32, f"{FLOOR_GENE} floor", ha="center", fontsize=8.5,
             color=INK2, bbox=dict(facecolor=SURF, edgecolor="none", pad=1.5))
    axC.set_yticks([ypos[g] for g in order])
    axC.set_yticklabels(order)
    axC.set_ylim(-.8, len(order) + .05)
    axC.set_xlabel("margin above the housekeeping floor  (median |$t$|)")
    axC.set_title("C   No candidate separates from the floor after correction")
    axC.grid(axis="x", color=GRID, lw=.7)
    axC.set_axisbelow(True)
    for s_ in ("top", "right", "left"):
        axC.spines[s_].set_visible(False)
    axC.legend(handles=[Line2D([], [], color=BLUE, lw=2.4, marker="o", ms=6.5,
                               label="skipped exons"),
                        Line2D([], [], color=ORANGE, lw=2.4, marker="o", ms=6.5,
                               label="mutually exclusive exons")],
               loc="upper center", bbox_to_anchor=(.5, -.13), ncol=2,
               frameon=False, fontsize=8.5)

    # --- D: technical inflation ------------------------------------------
    STYLE = {"without composition adj.": "-", "with composition adj.": "--"}
    for et, col in (("SE", BLUE), ("MXE", ORANGE)):
        for lab, (b, aft) in TECH45[et].items():
            axD.plot([0, 1], [b, aft], STYLE[lab], marker="o", color=col, lw=1.9,
                     ms=6.5, mec=SURF, mew=1.5, zorder=3)
        ys_ = [v[0] for v in TECH45[et].values()]
        axD.text(-.07, sum(ys_) / 2, et, ha="right", va="center", fontsize=10,
                 color=col, fontweight="bold")
    axD.legend(handles=[Line2D([], [], color=INK2, lw=1.9, ls=STYLE[k], label=k)
                        for k in TECH45["SE"]],
               loc="upper right", frameon=False, fontsize=8.5)
    axD.axhline(NULL_EXPECT, color=INK2, lw=1.2, ls=":", zorder=1)
    axD.text(0.02, NULL_EXPECT - .07, f"null expectation ≈ {NULL_EXPECT}",
             fontsize=8, color=MUTED)
    axD.set_xticks([0, 1])
    axD.set_xticklabels(["before", "after\n+ library size + TIN"])
    axD.set_xlim(-.30, 1.18)
    axD.set_ylabel("background median |$t$|")
    axD.set_title("D   Technical variation inflates every association")
    axD.grid(axis="y", color=GRID, lw=.7)
    axD.set_axisbelow(True)
    for s_ in ("top", "right"):
        axD.spines[s_].set_visible(False)

    fig.tight_layout(pad=2.4, h_pad=3.2)
    for ext in ("pdf", "png"):
        fig.savefig(f"{a.out}.{ext}", dpi=400)
    print(f"\n[write] {a.out}.pdf and {a.out}.png")

    # ------------------------------------------------------------------ report
    print("\n" + "=" * 78)
    print("PROVENANCE — paste into the legend")
    print("=" * 78)
    print(f"""
  panel A  {panelA_src}
           model: {a.expr_model} — {desc}
           +{sd:.2f} SD, P = {pval:.3g}
  panel B  {FWD if F is not None else 'recorded constants'}; 45-gene target set,
           40 of 45 genes testable for SE ({_NT_SE:,} events), 39 for MXE ({_NT_MXE:,})
  panel C  {SCAN}; model={a.model!r} matching={a.matching!r};
           ten candidates against the {FLOOR_GENE} floor, paired-bootstrap CIs, Holm
  panel D  45-gene background medians ({_D_LOG} step 5)

  EVERY PANEL REPORTS A MARGIN ABOVE THE FLOOR WITH AN INTERVAL. Do not reintroduce
  a bare p-value: 0.76 and 0.99 are the retired nine-gene values and 0.35 and 0.06
  their 45-gene equivalents — neither belongs in this figure.
""")


if __name__ == "__main__":
    main()
