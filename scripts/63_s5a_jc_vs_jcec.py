#!/usr/bin/env python3
"""
63_s5a_jc_vs_jcec.py -- build Supplementary Table S5(a) from files that already exist.

WHY THIS SCRIPT EXISTS, AND WHY IT DOES NOT REPORT SLOPES
---------------------------------------------------------
S5(a) was designed around Deming transfer slopes, for one bad reason: the project doc
that happened to record a JC-versus-JCEC comparison (`RESULT_magnitude_resolved.md`)
reported slopes, so the supplement copied that shape. Two problems followed.

1. Deming is an errors-in-variables fit. It needs a per-event standard error, which raw
   rMATS does not report, so slopes require phase 2 re-run over the JCEC tables -- a
   cluster job, for a comparison the main text never makes.
2. More importantly, a slope is the WRONG QUANTITY for this row. The paper states
   explicitly that effect-size magnitude is never compared across cohorts, because the
   read-length normalisation constant does not cancel across them. A supplement whose
   job is "does this pipeline choice change the answer" should therefore report the
   things the paper actually claims: how many events are called, and whether they move
   in the same direction. That is exactly what parts (b) and (c) report.

So this script answers (a) on the same footing as (b) and (c) -- significant counts,
tested-set size, direction agreement and Spearman rho -- directly from the
`{CLASS}.MATS.JC.txt` and `{CLASS}.MATS.JCEC.txt` files rMATS already wrote. No phase 2,
no standard errors, no new quantification. Minutes, not a job.

WHAT THIS ROW DOES NOT SUPPORT -- state these in the legend, do not let them be inferred
----------------------------------------------------------------------------------------
* It is UNADJUSTED. The manuscript's counts are covariate-adjusted; these are raw rMATS
  output. Part (b) has the precedent and the calibration: there the same choice moved the
  raw count -2.9% and the adjusted count -0.6%, i.e. adjustment DAMPS the difference. So
  an unadjusted comparison is the conservative one -- it overstates sensitivity, not
  understates it.
* It is computed on the rMATS run AS RUN, which for some arms includes samples later
  removed by the genotype-identity screen. This does not bias the comparison, because
  both arms are the same samples quantified two ways -- the design is paired and the
  sample set cancels. It does mean no count in this row is a manuscript number, and none
  should be quoted as one.
* It says nothing about magnitude. That is deliberate (see above), not an omission.

SELF-CHECK
----------
`probe_rl50/SE.MATS.JC.txt` is a file whose correct answer is already known: 172
significant skipped exons, confirmed twice by `62_build_s5.py --parts b`. If that
directory is present, this script must reproduce 172 from it or it refuses to report
anything. A reader that cannot hit a known value is not a reader to trust with an
unknown one.

USAGE
-----
    conda activate rnaseq_tools        # pandas + numpy; no scipy, no statsmodels
    cd /data/project/chureeratp/RNA/HEART
    python 63_s5a_jc_vs_jcec.py --list          # what it found, no computation
    python 63_s5a_jc_vs_jcec.py                 # SE, the five arms that carry the paper
    python 63_s5a_jc_vs_jcec.py --classes all   # all five event classes
    python 63_s5a_jc_vs_jcec.py --extra         # add the ICM, HCMvsNF and LBNL arms
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

# ---------------------------------------------------------------- locked conventions
# Identical to 62_build_s5.py and to every other script in this project. The exon-level
# collapse is the paper's convention; see Methods and S5(c).
COORD = {"SE":   ["chr", "exonStart_0base", "exonEnd"],
         "MXE":  ["chr", "1stExonStart_0base", "1stExonEnd",
                  "2ndExonStart_0base", "2ndExonEnd"],
         "A3SS": ["chr", "longExonStart_0base", "longExonEnd", "shortES", "shortEE"],
         "A5SS": ["chr", "longExonStart_0base", "longExonEnd", "shortES", "shortEE"],
         "RI":   ["chr", "riExonStart_0base", "riExonEnd"]}
CLASSES = ["SE", "MXE", "A3SS", "A5SS", "RI"]
FDR_CUT, DPSI_CUT = 0.05, 0.10

# The five arms that carry the paper. MAGNet first: it is the discovery cohort and the
# internal control -- it is quantified once, so anything that moves in it is the
# quantification choice and nothing else.
ARMS = [("MAGNet (discovery)", "rmats_magnet"),
        ("Colorado",           "rmats_Colorado"),
        ("Brigham",            "rmats_Brigham"),
        ("Fuwai",              "rmats_Fuwai"),
        ("Stanford",           "rmats_Stanford")]
# Available, deliberately not in the default table. LBNL shares 24 of 33 individuals with
# MAGNet and is excluded project-wide; the ICM and HCM-vs-NF arms are separate contrasts
# that no part of S5 is about.
EXTRA = [("MAGNet HCM vs NF",     "rmats_HCMvsNF"),
         ("Colorado DCM vs ICM",  "rmats_Colorado_DCMvsICM"),
         ("Colorado ICM vs NF",   "rmats_Colorado_ICMvsNF"),
         ("LBNL (not independent)", "rmats_LBNL")]
# Excluded project-wide on independent grounds -> RESULT_colorado_mxe_excluded.md
DROP_COMBOS = {("Colorado", "MXE")}
# A file whose answer is already known. Non-negotiable gate.
ANCHOR = ("probe_rl50", "SE", 172)

# Read length. Set TBL_SUFFIX=_rl125 to use each arm's re-quantified
# directory WHERE ONE EXISTS. Only MAGNet was mis-declared, so only MAGNet
# has one; the rest are identical at both read lengths and stay as published.
# ANCHOR above is deliberately NOT re-based: probe_rl50/SE = 172 is a fixed
# known value and the gate must pass at either setting.
SUF = os.environ.get("TBL_SUFFIX", "")


def load(path, et):
    """Collapse a raw rMATS table to one row per event under the paper's convention."""
    try:
        t = pd.read_csv(path, sep="\t", low_memory=False)
    except Exception as e:                                  # noqa: BLE001
        return None, f"unreadable ({e.__class__.__name__})"
    need = COORD[et] + ["FDR", "IncLevelDifference"]
    miss = [c for c in need if c not in t.columns]
    if miss:
        return None, f"missing columns {miss}"
    t["ekey"] = t[COORD[et]].astype(str).agg(":".join, axis=1)
    t["FDR"] = pd.to_numeric(t["FDR"], errors="coerce")
    t["dpsi"] = pd.to_numeric(t["IncLevelDifference"], errors="coerce")
    t = t[np.isfinite(t["FDR"]) & np.isfinite(t["dpsi"])]
    if not len(t):
        return None, "no rows with finite FDR and dPSI"
    t = t.loc[t.groupby("ekey")["FDR"].idxmin()]
    return t[["ekey", "FDR", "dpsi"]].reset_index(drop=True), None


def nsig(T):
    return int(((T["FDR"] < FDR_CUT) & (T["dpsi"].abs() >= DPSI_CUT)).sum())


def spearman(x, y):
    """Rank correlation without scipy -- Pearson on ranks. Ties get average ranks."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 3:
        return float("nan")
    rx = pd.Series(x[ok]).rank().to_numpy()
    ry = pd.Series(y[ok]).rank().to_numpy()
    if rx.std() == 0 or ry.std() == 0:
        return float("nan")
    return float(np.corrcoef(rx, ry)[0, 1])


def pct(new, old):
    return float("nan") if not old else 100.0 * (new - old) / old


def compare(d, et):
    """One arm, one event class: JC against JCEC on the same samples."""
    pj = os.path.join(d, f"{et}.MATS.JC.txt")
    pe = os.path.join(d, f"{et}.MATS.JCEC.txt")
    for p in (pj, pe):
        if not (os.path.exists(p) and os.path.getsize(p)):
            return None, f"{os.path.basename(p)} missing or empty"
    J, ej = load(pj, et)
    if J is None:
        return None, f"JC: {ej}"
    E, ee = load(pe, et)
    if E is None:
        return None, f"JCEC: {ee}"

    m = J.merge(E, on="ekey", suffixes=("_jc", "_ec"))
    # Direction agreement is only meaningful where there is a direction to agree about.
    # Events at dPSI == 0 exactly have no sign; drop them rather than scoring them.
    mm = m[(m["dpsi_jc"] != 0) & (m["dpsi_ec"] != 0)]
    same = int((np.sign(mm["dpsi_jc"]) == np.sign(mm["dpsi_ec"])).sum())
    # ... and the subset that actually carries the claim: events JC called significant.
    s = mm[(mm["FDR_jc"] < FDR_CUT) & (mm["dpsi_jc"].abs() >= DPSI_CUT)]
    same_s = int((np.sign(s["dpsi_jc"]) == np.sign(s["dpsi_ec"])).sum())

    return dict(
        tested_jc=len(J), tested_ec=len(E), tested_pct=pct(len(E), len(J)),
        sig_jc=nsig(J), sig_ec=nsig(E), sig_pct=pct(nsig(E), nsig(J)),
        shared=len(m), dir_n=len(mm), dir_same=same,
        dir_pct=100.0 * same / len(mm) if len(mm) else float("nan"),
        sig_dir_n=len(s), sig_dir_same=same_s,
        sig_dir_pct=100.0 * same_s / len(s) if len(s) else float("nan"),
        rho=spearman(m["dpsi_jc"], m["dpsi_ec"]),
        rho_sig=spearman(s["dpsi_jc"], s["dpsi_ec"]),
    ), None


def run_anchor():
    """Reproduce a known value or refuse to report. Returns (ok, message)."""
    d, et, want = ANCHOR
    p = os.path.join(d, f"{et}.MATS.JC.txt")
    if not os.path.exists(p):
        return None, (f"anchor directory {d}/ not present -- cannot verify the reader "
                      f"against a known value")
    T, err = load(p, et)
    if T is None:
        return False, f"anchor {p} unreadable: {err}"
    got = nsig(T)
    if got != want:
        return False, (f"anchor {p} gives {got} significant {et}, expected {want} "
                       f"(62_build_s5.py --parts b confirmed {want} twice)")
    return True, f"anchor ok: {p} -> {got} significant {et}, as recorded"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--classes", default="SE",
                    help="comma-separated, or 'all' (default: SE)")
    ap.add_argument("--extra", action="store_true",
                    help="also report the ICM, HCMvsNF and LBNL arms")
    ap.add_argument("--list", action="store_true", help="show what was found, then stop")
    ap.add_argument("--out", default="S5a" + SUF)
    ap.add_argument("--skip-anchor", action="store_true",
                    help="report without the known-value check (say so if you use it)")
    a = ap.parse_args()

    cls = CLASSES if a.classes.strip().lower() == "all" else \
        [c.strip().upper() for c in a.classes.split(",") if c.strip()]
    bad = [c for c in cls if c not in COORD]
    if bad:
        sys.exit(f"unknown event class(es): {bad}; choose from {CLASSES} or 'all'")
    arms = ARMS + (EXTRA if a.extra else [])

    # ---- READ-LENGTH RE-BASING, PER ARM, added 21 Oct ------------------
    # Only MAGNet was mis-declared (--readLength 50 on 125 bp reads), so
    # only MAGNet has an _rl125 twin. Colorado, Brigham, Fuwai and
    # Stanford were always quantified at their true lengths and their
    # files are IDENTICAL at both read lengths.
    #
    # Four of five rows of S5a_rl125.tsv are therefore byte-identical to
    # S5a.tsv. That is correct, and the table has to SAY it is correct --
    # a reader diffing the two files must be able to tell "nothing
    # happened" from "nothing should happen". So the basis goes into the
    # arm name, not just into this log.
    if SUF:
        _rebased, _kept, _arms2 = [], [], []
        for _n, _d in arms:
            if os.path.isdir(_d + SUF):
                _rebased.append(_n)
                _arms2.append((f"{_n} @{SUF.lstrip('_')}", _d + SUF))
            else:
                _kept.append(_n)
                _arms2.append((f"{_n} @as-published", _d))
        print(f"\n  TBL_SUFFIX={SUF}: re-based {len(_rebased)} of "
              f"{len(arms)} arms")
        for _n in _rebased:
            print(f"    re-based      {_n}")
        for _n in _kept:
            print(f"    as published  {_n}  (no {SUF} twin; this arm was "
                  f"never mis-declared)")
        if not _rebased:
            print(f"\n  >>> (a) NOT BUILT. TBL_SUFFIX={SUF} re-based NOTHING.")
            print("      No arm directory carries that suffix, so this would")
            print("      be the 50 bp table written under a 125 bp name.")
            return
        arms = _arms2

    print("=" * 78)
    print("S5(a)  junction counts (JC) versus junction + exon-body counts (JCEC)")
    print("=" * 78)
    print(f"  event classes : {', '.join(cls)}")
    print(f"  arms          : {len(arms)}"
          f"{'  (+extra)' if a.extra else '  (the five that carry the paper)'}")
    print(f"  rule          : FDR < {FDR_CUT} and |dPSI| >= {DPSI_CUT}, exon-level collapse")
    print("  NOTE          : unadjusted, and on the rMATS run as run. Paired design, so")
    print("                  the sample set cancels -- but no number here is a")
    print("                  manuscript number. See this script's header.")

    print("\n  arm / files found:")
    present = []
    for name, d in arms:
        if not os.path.isdir(d):
            print(f"    {name:<24} *** directory {d}/ not found ***")
            continue
        have = [c for c in cls
                if os.path.exists(os.path.join(d, f"{c}.MATS.JC.txt"))
                and os.path.exists(os.path.join(d, f"{c}.MATS.JCEC.txt"))]
        missing = [c for c in cls if c not in have]
        print(f"    {name:<24} {d}/  both tables: "
              f"{', '.join(have) or 'NONE'}"
              f"{'   missing: ' + ', '.join(missing) if missing else ''}")
        if have:
            present.append((name, d, have))
    if a.list:
        print("\n[--list] stopping here.")
        return

    if not present:
        print("\n  >>> (a) NOT BUILT. No arm has both a JC and a JCEC table for the")
        print("      requested classes. Nothing to compare.")
        return

    # ------------------------------------------------------------------ the gate
    ok, msg = run_anchor()
    if ok is False and not a.skip_anchor:
        print(f"\n  >>> (a) NOT BUILT. {msg}")
        print("      The reader is wrong, so every number it would print is suspect.")
        print("      Fix the reader before trusting it with an unknown value.")
        return
    print(f"\n  [{'ok' if ok else 'warn'}] {msg}")
    if a.skip_anchor:
        print("  [warn] --skip-anchor: no known-value check was performed")

    # ------------------------------------------------------------------ compute
    rows, notes = [], []
    for et in cls:
        res = {}
        for name, d, have in present:
            if et not in have:
                continue
            if (name, et) in DROP_COMBOS:
                print(f"\n  {name} {et}: excluded project-wide "
                      f"(RESULT_colorado_mxe_excluded.md) -- not computed")
                continue
            r, err = compare(d, et)
            if r is None:
                print(f"\n  [warn] {name} {et}: {err}")
                notes.append(f"(a) {name} {et}: {err}")
                continue
            res[name] = r
        if not res:
            continue

        print("\n" + "-" * 78)
        print(f"  {et}: significant events and tested-set size")
        print("-" * 78)
        print(f"  {'arm':<24}{'JC':>8}{'JCEC':>8}{'change':>9}"
              f"{'tested JC':>11}{'tested JCEC':>13}{'change':>9}")
        for name, r in res.items():
            print(f"  {name:<24}{r['sig_jc']:>8}{r['sig_ec']:>8}"
                  f"{r['sig_pct']:>+8.1f}%"
                  f"{r['tested_jc']:>11}{r['tested_ec']:>13}{r['tested_pct']:>+8.1f}%")

        print(f"\n  {et}: do the same events move the same way?")
        print("  The 'all shared' column includes every tested event, most of which sit")
        print("  near dPSI = 0 where the sign is noise in both quantifications -- expect")
        print("  it to be mediocre and do not read it as disagreement. The JC-significant")
        print("  column is the one that carries the claim.")
        print(f"  {'arm':<24}{'all shared':>22}{'rho':>7}"
              f"{'JC-significant':>22}{'rho':>7}")
        for name, r in res.items():
            allcol = f"{r['dir_same']}/{r['dir_n']} ({r['dir_pct']:.2f}%)"
            sigcol = f"{r['sig_dir_same']}/{r['sig_dir_n']} ({r['sig_dir_pct']:.2f}%)"
            print(f"  {name:<24}{allcol:>22}{r['rho']:>7.3f}"
                  f"{sigcol:>22}{r['rho_sig']:>7.3f}")

        for name, r in res.items():
            rows.append(dict(
                part="a. JC vs JCEC", event=et, arm=name,
                sig_JC=r["sig_jc"], sig_JCEC=r["sig_ec"],
                sig_change=f"{r['sig_pct']:+.1f}%",
                tested_JC=r["tested_jc"], tested_JCEC=r["tested_ec"],
                tested_change=f"{r['tested_pct']:+.1f}%",
                dir_agree=f"{r['dir_same']}/{r['dir_n']} ({r['dir_pct']:.2f}%)",
                rho=round(r["rho"], 3),
                dir_agree_sig=f"{r['sig_dir_same']}/{r['sig_dir_n']} "
                              f"({r['sig_dir_pct']:.2f}%)",
                rho_sig=round(r["rho_sig"], 3)))

        # ---------------------------------------------------- read it back honestly
        if "MAGNet (discovery)" in res:
            m = res["MAGNet (discovery)"]
            others = [r for n, r in res.items() if n != "MAGNet (discovery)"]
            print(f"\n  Internal control. MAGNet is quantified once, so its own JC-to-JCEC"
                  f"\n  change ({m['sig_pct']:+.1f}% significant,"
                  f" {m['sig_dir_pct']:.2f}% same direction among called"
                  f"\n  events) is the quantification choice and nothing else.")
            if others:
                sp = [abs(r["sig_pct"]) for r in others if np.isfinite(r["sig_pct"])]
                if sp:
                    print(f"  Replication arms move by {min(sp):.1f}% to {max(sp):.1f}%"
                          f" against MAGNet's {abs(m['sig_pct']):.1f}%.")
        # Judge the row on the JC-significant subset. Agreement over every tested event
        # is dominated by events with no real effect, whose sign is noise in both
        # quantifications -- a low value there is expected and means nothing.
        dmin = min((r["sig_dir_pct"] for r in res.values()
                    if np.isfinite(r["sig_dir_pct"])), default=float("nan"))
        if np.isfinite(dmin):
            if dmin >= 95:
                print(f"\n  [ok] among JC-significant events, direction agreement is"
                      f" >= {dmin:.1f}% in every\n       arm: the two quantifications"
                      f" disagree about which way a called event\n       moves in at most"
                      f" {100 - dmin:.1f}% of cases.")
            else:
                print(f"\n  >>> among JC-significant events, direction agreement falls to"
                      f" {dmin:.1f}% in the\n      worst arm. That is low enough to matter."
                      f" Do not summarise this row as\n      'the choice does not change"
                      f" the answer' -- report the number and let\n      the reader judge.")
                notes.append(f"(a) {et}: worst-arm agreement among JC-significant events "
                             f"{dmin:.1f}%")

    print("\n" + "=" * 78)
    print("RESULT")
    print("=" * 78)
    if not rows:
        print("  nothing built")
        for n in notes:
            print(f"  >>> {n}")
        return
    T = pd.DataFrame(rows)
    T.to_csv(f"{a.out}.tsv", sep="\t", index=False)
    with pd.option_context("display.width", 250, "display.max_columns", 50):
        print(T.to_string(index=False))
    print(f"\n[write] {a.out}.tsv")
    for n in notes:
        print(f"  >>> {n}")
    print("\n  Merge into S5 as part (a). The legend must say: unadjusted, paired on the")
    print("  same samples, and reporting counts and direction rather than magnitude --")
    print("  because this paper does not compare magnitude across cohorts.")


if __name__ == "__main__":
    main()
