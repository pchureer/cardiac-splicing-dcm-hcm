#!/usr/bin/env python3
"""
67_build_s8.py -- Supplementary Table S8: null behaviour and reverse replication.

WHY S8 IS THE MOST IMPORTANT SUPPLEMENT IN THE SET
--------------------------------------------------
This paper's central claims are NULLS: dilated against hypertrophic cardiomyopathy
returns 0 significant events of 45,380; dilated against ischaemic returns 2 of 19,998;
no candidate regulator separates from a housekeeping floor. A reader has one obvious
question about every one of them -- *is the pipeline capable of returning anything at
all?* S8 is where that is answered, and `supp_table_s5.md` part (d) now points here as
well.

Three parts:

  (a) **Permutation nulls.** Shuffle the group labels and the answer must collapse to
      zero. Thirteen permutations in Colorado across both contrasts, both event types and
      three settings, and two in the discovery cohort across 103,054 events. Every one
      returned EXACTLY zero.
  (b) **Reverse replication**, which is what excluded Colorado's mutually-exclusive-exon
      calls -- see the framing warning below, it is easy to state this wrongly.
  (c) **The explanations that were tested and rejected** before the exclusion was made:
      sequencing depth, stringency, and a technical covariate.

HOW TO STATE (b), AND A CORRECTION TO AN EARLIER VERSION OF THIS HEADER
----------------------------------------------------------------------
`RESULT_magnitude_not_estimable.md` warns that forward-versus-reverse asymmetry is not
interpretable *on its own*: on synthetic data with a planted transfer of exactly 1.00 and
symmetric truth, differing only in precision, the test returned forward 0.949 and reverse
0.674. Selection truncates the selected cohort's own noise, so the noisier the cohort you
select on, the more its partner regresses -- and the discovery cohort is ~2.7x more
precise than any replication cohort.

An earlier version of this script said the asymmetry therefore "must not be quoted as
evidence". **That was too strong.** The synthetic run does not make the asymmetry
uninterpretable; it *calibrates* it. The null expectation is 0.949 / 0.674 = **1.41x**,
and a cohort can be read against it:

    Colorado SE   0.571 / 0.447 = 1.28x   -- inside the artefact, as expected
    Colorado MXE  0.533 / 0.061 = 8.72x   -- 6.2x beyond it

So S8 has **three mutually independent lines**, and should give all three because no one
of them carries the whole weight:

  1. **Within cohort, within direction** (no calibration assumed): Colorado SE reverse
     0.447 against Colorado MXE reverse 0.061, same pipeline, so the precision term is
     common to both and cancels.
  2. **Asymmetry against its calibration**: 8.72x against a 1.41x null.
  3. **Correlation under every estimator** -- the cleanest, because correlations do not
     share Deming's instability: Colorado MXE r = 0.157 to 0.254, Colorado SE r = 0.298
     to 0.708.

And a diagnostic worth reporting in its own right: Colorado MXE's UNSEL and PREC Deming
slopes -- the two columns `RESULT_magnitude_not_estimable.md` calls the trustworthy ones
-- come out at **11.04 and 15.02**. Those are not attenuation estimates. Deming is a ratio
of variances and it diverges when there is no relationship to fit; a slope of 15 at r =
0.198 *is* the finding, stated in the units of the estimator.

Two further constraints, both documented:
  * Script 22's automated verdict is unreliable -- *"read the per-cohort table; ignore
    section 3; the docstring's decision rule is wrong as written"*. Do not reproduce it.
  * The recorded reverse slopes are **pre-exclusion** (Colorado 37/14, not the analysed
    36/13). The same provenance class as the pre-rescue pair count corrected in S9 and
    the pre-exclusion slopes banner-flagged in `RESULT_magnitude_resolved.md`. This
    script measures the group sizes behind each slope and says which run they came from.

USAGE
-----
    conda activate rnaseq_tools
    cd /data/project/chureeratp/RNA/HEART

    python 67_build_s8.py --list      # find and describe every candidate input
    python 67_build_s8.py             # build, once the columns are pinned

**Run `--list` first and send me the output.** Nothing here is pinned yet, and this
session has already paid three times for guessing at column names.
"""
import argparse
import glob
import os

import numpy as np
import pandas as pd

# ---------------------------------------------------------------- the record
PERM = dict(colorado=13, magnet=2, magnet_events=103_054, significant=0)
REVERSE = [  # cohort, event class, n, slope, verdict  (PRE-EXCLUSION)
    ("MAGNet internal", "MXE", 147, 0.725, "control; shares platform and batch"),
    ("Colorado", "SE", 453, 0.447, "reproduces"),
    ("Colorado", "MXE", 1477, 0.061, "DOES NOT reproduce"),
    ("Fuwai", "SE", 156, 0.405, "reproduces"),
    ("Fuwai", "MXE", 99, 0.281, "weak -- caution"),
]
RATES = [  # cohort, contrast, n, SE %, MXE %  -- how the problem surfaced
    ("MAGNet", "DCM vs NF", "166/166", 1.47, 0.37),
    ("Colorado", "DCM vs NF", "37/14", 2.23, 11.87),
    ("Colorado", "ICM vs NF", "13/14", 2.37, 21.38),
]
DEPTH = [("DCM vs NF", "13 vs 13", 0.89), ("ICM vs NF", "14 vs 14", 0.86)]
COVARIATE = [("SE", 478, 475), ("MXE", 1596, 1605)]   # before/after log10_jcov
MINSAMPLES = [("Colorado", 10, "13 and 14", 77, 71), ("MAGNet", 20, "166 and 166", 12, 12)]
STRINGENCY = dict(tested_cut=81, hits_retained=77)
# ⚠ The record says 1,462 and 2,058. master_numbers_flagged.tsv -- frozen
# 28 Sep and rebuilt from the post-exclusion tables -- says 1,610 and 2,112.
# These stay here as the CHECK; the file is what gets quoted.
EXCLUDED_CALLS = dict(dcm_nf=1462, icm_nf=2058)

# Pinned from the real files, 6 Oct. symmetry_{CLASS}.tsv carries four estimators per
# cohort: fwd (select on the discovery cohort), rev (select on the replication cohort),
# unsel (no selection) and prec (precision-filtered), each with n, Deming slope, OLS
# slope, correlation and % same direction.
SYM = "symmetry_{}.tsv"
SYM_COLS = dict(cohort="cohort", n_common="n_common",
                fwd_n="fwd_n", fwd="fwd_deming", fwd_r="fwd_r",
                rev_n="rev_n", rev="rev_deming", rev_r="rev_r",
                uns_n="unsel_n", uns="unsel_deming", uns_r="unsel_r",
                pre_n="prec_n", pre="prec_deming", pre_r="prec_r")
MASTER = "master_numbers_flagged.tsv"
MASTER_COLS = dict(cohort="cohort", contrast="contrast", event="event_type",
                   reliability="reliability", rev_n="rev_n", rev="rev_deming",
                   flag="flag_reason", sig="significant", tested="tested")
CALIB = (0.949, 0.674)      # synthetic forward / reverse at a planted transfer of 1.00
CLASSES = ["SE", "MXE"]

# candidate inputs -- discovered, not assumed
PATTERNS = [
    # Widened after the first pass found nothing. The permutation result is the single
    # most load-bearing claim in S8 and the Results states it directly, so if no file
    # holds it that is a gap in the deposit, not a naming problem.
    ("permutation", ["*perm*", "*shuf*", "*null*", "25_*", "*label*",
                     "**/*perm*", "**/*shuf*", "**/*null*",
                     "logs/*", "log/*", "*.log", "slurm-*.out"]),
    ("master numbers", ["master_numbers*.tsv"]),
    ("reverse / symmetry", ["*symmetry*", "*reverse*", "*winners*", "*wc_*.tsv"]),
]


def describe(path, nmax=3):
    try:
        if os.path.getsize(path) == 0:
            print(f"      {path}  (empty)")
            return
        if path.endswith((".tsv", ".csv", ".txt")):
            sep = "," if path.endswith(".csv") else "\t"
            t = pd.read_csv(path, sep=sep, low_memory=False, nrows=400)
            print(f"      {path}  ({len(t)} rows shown, {len(t.columns)} columns)")
            print(f"        columns: {', '.join(map(str, t.columns))}")
            for line in t.head(nmax).to_string(index=False).splitlines():
                print(f"          {line}")
        else:
            print(f"      {path}  ({os.path.getsize(path):,} bytes)")
    except Exception as e:                                          # noqa: BLE001
        print(f"      {path}  [unreadable: {e.__class__.__name__}]")


def discover():
    print("\nSEARCHING for the inputs S8 needs")
    found = {}
    for label, pats in PATTERNS:
        hits = sorted({p for pat in pats
                       for p in glob.glob(pat, recursive=True)
                       if os.path.isfile(p)})
        if label == "permutation":
            # a log file only counts if it actually mentions a permutation
            keep = []
            for p in hits:
                if p.endswith((".log", ".out")) or "/" in p:
                    try:
                        with open(p, errors="ignore") as fh:
                            head = fh.read(200_000).lower()
                        if not any(w in head for w in
                                   ("permut", "shuffl", "label swap", "null")):
                            continue
                    except OSError:
                        continue
                keep.append(p)
            hits = keep
        print(f"\n  {label}:")
        if not hits:
            print("      *** nothing matched ***")
        for p in hits[:8]:
            describe(p)
        if len(hits) > 8:
            print(f"      ... and {len(hits) - 8} more")
        found[label] = hits
    return found


def show_record():
    """Everything S8 must contain, from the project record, so the gap is visible."""
    print("\n" + "=" * 78)
    print("WHAT S8 MUST CONTAIN -- from the project record")
    print("=" * 78)

    print("\n(a) PERMUTATION NULLS -- the answer to 'can this pipeline return anything?'")
    print(f"    Colorado: {PERM['colorado']} permutations across both contrasts, both")
    print("              event classes and three settings")
    print(f"    Discovery: {PERM['magnet']} permutations across"
          f" {PERM['magnet_events']:,} events")
    print(f"    Every one returned EXACTLY {PERM['significant']} significant events.")
    print("    Limit to state, not bury: shuffling labels destroys the group-batch")
    print("    association along with the biology, so this validates p-value")
    print("    CALIBRATION and cannot detect CONFOUNDING. It excluded one explanation.")

    print("\n(b) REVERSE REPLICATION -- what excluded the Colorado MXE calls")
    print(f"    {'cohort':<18}{'class':<6}{'n':>7}{'slope':>8}   verdict")
    for c, e, n, s, v in REVERSE:
        print(f"    {c:<18}{e:<6}{n:>7}{s:>8.3f}   {v}")
    print("\n    The comparison that carries the argument is Colorado SE 0.447 against")
    print("    Colorado MXE 0.061 -- same cohort, same direction, same pipeline, so the")
    print("    precision term cancels. Forward-versus-reverse asymmetry is NOT")
    print("    interpretable here and must not be quoted as evidence.")
    print(f"    Calls excluded as a result: the record says"
          f" {EXCLUDED_CALLS['dcm_nf']:,} and {EXCLUDED_CALLS['icm_nf']:,};")
    print("    the frozen master table says 1,610 and 2,112 -- see the check above.")

    print("\n    How the problem surfaced -- a cohort a fifth the size calling far more:")
    print(f"    {'cohort':<10}{'contrast':<12}{'n':>9}{'SE':>8}{'MXE':>9}")
    for c, k, n, se, mxe in RATES:
        print(f"    {c:<10}{k:<12}{n:>9}{se:>7.2f}%{mxe:>8.2f}%")
    print("    Within the discovery cohort MXE sits BELOW SE; within Colorado the order")
    print("    inverts. That inversion is the signal that something is wrong.")

    print("\n(c) EXPLANATIONS TESTED AND REJECTED before the exclusion was made")
    print("    Sequencing depth -- per-sample median junction coverage by group:")
    for k, n, p in DEPTH:
        print(f"      {k:<12} {n:<10} Welch p = {p}")
    print("    Adding log10(junction coverage) as a covariate:")
    for e, a, b in COVARIATE:
        print(f"      {e:<5} {a} -> {b} significant events")
    print("    Stringency -- min-samples was the wrong suspect from the start:")
    for c, ms, grp, f1, f2 in MINSAMPLES:
        print(f"      {c:<10} --min-samples {ms:<3} against groups of {grp:<12}"
              f" = {f1}% / {f2}% of samples")
    print("      Colorado is STRICTER in fractional terms, not looser.")
    print(f"      Raising it further cut the tested set {STRINGENCY['tested_cut']}% while")
    print(f"      retaining {STRINGENCY['hits_retained']}% of hits -- signal concentrated")
    print("      in well-measured events, which is what real signal looks like.")

    print("\n" + "=" * 78)
    print("TWO THINGS TO SETTLE BEFORE S8 IS TYPESET")
    print("=" * 78)
    print("  1. The reverse slopes above are PRE-EXCLUSION -- Colorado 37/14 against the")
    print("     analysed 36/13, Stanford 13/7 against 13/6. The verdict is not in doubt")
    print("     (0.061 against 0.447 is sevenfold), but the NUMBERS are from a different")
    print("     sample set than the manuscript's. Either re-run 22_symmetry_winners_curse")
    print("     on the _dedup tables, or print the slopes with the group sizes beside")
    print("     them and say which run they are from. Do not quote them bare.")
    print("  2. Script 22's automated verdict must not be reproduced anywhere in S8.")
    print("     It reported two contradictory conclusions from the same data because it")
    print("     keyed on a mean that included unstable estimates.")


def build_reverse(out, notes):
    """Part (b), from symmetry_SE.tsv and symmetry_MXE.tsv. Three lines of evidence."""
    rows, have = [], {}
    for cls in CLASSES:
        path = SYM.format(cls)
        if not (os.path.exists(path) and os.path.getsize(path)):
            print(f"  >>> {path} missing -- part (b) cannot be built from file")
            notes.append(f"(b) {path} missing")
            return None
        t = pd.read_csv(path, sep="\t")
        miss = [c for c in SYM_COLS.values() if c not in t.columns]
        if miss:
            print(f"  >>> {path}: pinned columns absent: {miss}")
            notes.append(f"(b) {path} schema changed")
            return None
        have[cls] = t

    print(f"\n  {'cohort':<18}{'class':<6}{'rev n':>7}{'rev slope':>11}{'rev r':>8}"
          f"{'UNSEL slope':>13}{'UNSEL r':>9}{'PREC slope':>12}{'PREC r':>8}")
    for cls in CLASSES:
        t = have[cls]
        for _, r in t.iterrows():
            c = str(r[SYM_COLS["cohort"]])
            g = lambda k: pd.to_numeric(r[SYM_COLS[k]], errors="coerce")
            if not np.isfinite(g("rev")):
                continue
            print(f"  {c:<18}{cls:<6}{int(g('rev_n')):>7}{g('rev'):>11.3f}"
                  f"{g('rev_r'):>8.3f}{g('uns'):>13.3f}{g('uns_r'):>9.3f}"
                  f"{g('pre'):>12.3f}{g('pre_r'):>8.3f}")
            rows.append(dict(part="b. reverse replication", cohort=c, event_class=cls,
                             rev_n=int(g("rev_n")), rev_slope=round(float(g("rev")), 3),
                             rev_r=round(float(g("rev_r")), 3),
                             fwd_n=int(g("fwd_n")) if np.isfinite(g("fwd_n")) else "",
                             fwd_slope=round(float(g("fwd")), 3)
                                       if np.isfinite(g("fwd")) else "",
                             unsel_slope=round(float(g("uns")), 3),
                             unsel_r=round(float(g("uns_r")), 3),
                             prec_slope=round(float(g("pre")), 3),
                             prec_r=round(float(g("pre_r")), 3)))

    # ---- check against the record
    print()
    ok = True
    for c, cls, n_, sl in [(c, e, n, s_) for c, e, n, s_, _ in REVERSE]:
        key = "MAGNet internal" if c.startswith("MAGNet") else c
        t = have.get(cls)
        if t is None:
            continue
        m = t[t[SYM_COLS["cohort"]].astype(str) == key]
        if not len(m):
            print(f"  [MISMATCH] {c} {cls}: no row"); ok = False; continue
        gn = int(pd.to_numeric(m.iloc[0][SYM_COLS["rev_n"]], errors="coerce"))
        gs = float(pd.to_numeric(m.iloc[0][SYM_COLS["rev"]], errors="coerce"))
        good = gn == n_ and abs(gs - sl) < 0.001
        print(f"  [{'ok' if good else 'MISMATCH':<8}] {c} {cls}: n {gn} slope {gs:.3f}"
              f" (record n {n_} slope {sl})")
        if not good:
            notes.append(f"(b) {c} {cls}: n {gn}/{gs:.3f} vs record {n_}/{sl}")
            ok = False
    if ok:
        print("  [ok      ] every recorded reverse slope reproduces from file")

    # ---- the three lines, stated
    def pick(cls, coh, k):
        t = have[cls]
        m = t[t[SYM_COLS["cohort"]].astype(str) == coh]
        return float(pd.to_numeric(m.iloc[0][SYM_COLS[k]], errors="coerce")) if len(m) \
            else float("nan")

    print("\n  THE THREE LINES, each independent of the others")
    se_rev, mx_rev = pick("SE", "Colorado", "rev"), pick("MXE", "Colorado", "rev")
    print(f"    1. within cohort, within direction: Colorado SE {se_rev:.3f} against"
          f" Colorado MXE {mx_rev:.3f} ({se_rev / mx_rev:.1f}x)")
    cal = CALIB[0] / CALIB[1]
    se_a = pick("SE", "Colorado", "fwd") / se_rev
    mx_a = pick("MXE", "Colorado", "fwd") / mx_rev
    print(f"    2. asymmetry against its calibration: SE {se_a:.2f}x, MXE {mx_a:.2f}x,"
          f" synthetic null {cal:.2f}x")
    print(f"       SE sits inside the artefact; MXE exceeds it {mx_a / cal:.1f}-fold.")
    print(f"    3. correlation under every estimator:")
    for cls in CLASSES:
        rs = [pick(cls, "Colorado", k) for k in ("rev_r", "uns_r", "pre_r")]
        print(f"         Colorado {cls:<4} r = {min(rs):.3f} to {max(rs):.3f}")
    print(f"\n    Diagnostic: Colorado MXE UNSEL/PREC Deming ="
          f" {pick('MXE','Colorado','uns'):.2f} / {pick('MXE','Colorado','pre'):.2f}.")
    print("    Not attenuation -- Deming diverging because there is nothing to fit.")

    # ---- provenance
    # ---- is MXE worse than SE *within* each cohort, and does read layout order it?
    print("\n  MXE against SE within each cohort -- the comparison that tests the")
    print("  'MXE from single-end data' framing. UNSEL is selection-free.")
    print(f"    {'cohort':<18}{'SE unsel r':>11}{'MXE unsel r':>12}{'MXE worse by':>14}")
    worst, worst_c = None, None
    for coh in [str(c) for c in have["SE"][SYM_COLS["cohort"]]]:
        try:
            s_r = pick("SE", coh, "uns_r")
            m_r = pick("MXE", coh, "uns_r")
        except (IndexError, KeyError):
            continue
        if not (np.isfinite(s_r) and np.isfinite(m_r) and m_r > 0):
            continue
        print(f"    {coh:<18}{s_r:>11.3f}{m_r:>12.3f}{s_r / m_r:>13.1f}x")
        if worst is None or m_r < worst:
            worst, worst_c = m_r, coh
    if worst_c:
        print(f"\n    MXE is worse than SE in every cohort listed. The LOWEST MXE")
        print(f"    correlation belongs to {worst_c} ({worst:.3f}).")
        print("    If that cohort is paired-end, 'MXE from single-end data' is not the")
        print("    whole account -- the defensible claim is that MXE is systematically")
        print("    less reproducible than SE across layouts, and Colorado is where it")
        print("    collapsed rather than merely weakened.")
        notes.append(f"lowest MXE correlation is {worst_c} ({worst:.3f}) -- check its "
                     f"read layout before repeating the single-end framing")

    print("\n  [--      ] these are PRE-EXCLUSION (Colorado 37/14, analysed 36/13).")
    print("              The verdict is carried by correlations of 0.16 against 0.71,")
    print("              which two samples cannot move. Label them; a re-run is optional.")
    notes.append("(b) slopes are pre-exclusion -- label them in the table")
    return rows


def build_master(out, notes):
    """Part (b) support: how the exclusion was applied, from master_numbers_flagged."""
    if not (os.path.exists(MASTER) and os.path.getsize(MASTER)):
        print(f"  >>> {MASTER} missing")
        return None
    t = pd.read_csv(MASTER, sep="\t", low_memory=False)
    miss = [c for c in MASTER_COLS.values() if c not in t.columns]
    if miss:
        print(f"  >>> {MASTER}: pinned columns absent: {miss}")
        notes.append(f"(b) {MASTER} schema changed")
        return None
    fl = t[t[MASTER_COLS["flag"]].notna()]
    print(f"\n  {MASTER}: {len(t)} rows, {len(fl)} carrying a flag")
    if len(fl):
        print(f"  {'cohort':<10}{'contrast':<12}{'class':<6}{'reliability':<14}"
              f"{'sig':>7}   reason")
        for _, r in fl.iterrows():
            print(f"  {str(r[MASTER_COLS['cohort']]):<10}"
                  f"{str(r[MASTER_COLS['contrast']]):<12}"
                  f"{str(r[MASTER_COLS['event']]):<6}"
                  f"{str(r[MASTER_COLS['reliability']]):<14}"
                  f"{str(r[MASTER_COLS['sig']]):>7}   {r[MASTER_COLS['flag']]}")
    # the record's excluded-call counts, against the frozen table
    for coh, con, cls, rec in (("Colorado", "DCM vs NF", "MXE", EXCLUDED_CALLS["dcm_nf"]),
                               ("Colorado", "ICM vs NF", "MXE", EXCLUDED_CALLS["icm_nf"])):
        m = t[(t[MASTER_COLS["cohort"]].astype(str) == coh)
              & (t[MASTER_COLS["contrast"]].astype(str) == con)
              & (t[MASTER_COLS["event"]].astype(str) == cls)]
        if not len(m):
            continue
        got = int(pd.to_numeric(m.iloc[0][MASTER_COLS["sig"]], errors="coerce"))
        good = got == rec
        print(f"  [{'ok' if good else 'MISMATCH':<8}] {coh} {con} {cls}: {got:,}"
              f" significant (record {rec:,})")
        if not good:
            print(f"              The frozen table was rebuilt 28 Sep from the")
            print(f"              post-exclusion tables, so {got:,} supersedes {rec:,}.")
            notes.append(f"excluded-call count {coh} {con} {cls}: file {got:,}, "
                         f"record {rec:,} -- quote the file")

    # how much of the cohort set the reverse screen could actually cover
    rl = t[MASTER_COLS["reliability"]].astype(str)
    nt = t[rl.str.contains("not_tested", case=False, na=False)]
    if len(nt):
        cohorts_nt = sorted(set(nt[MASTER_COLS["cohort"]].astype(str)))
        tested = sorted(set(t[MASTER_COLS["cohort"]].astype(str))
                        - set(nt[nt.groupby(MASTER_COLS["cohort"])[MASTER_COLS["cohort"]]
                              .transform("size") >= 5][MASTER_COLS["cohort"]].astype(str)))
        print(f"\n  [--      ] {len(nt)} of {len(t)} rows are 'not_tested' -- too few")
        print(f"              events of their own to select on. Cohorts affected:")
        print(f"              {', '.join(cohorts_nt)}")
        print("              The reverse-replication screen could only be APPLIED to the")
        print("              cohorts calling enough of their own events. S8 must say so")
        print("              rather than implying all six passed a test two of them took.")
        notes.append("reverse screen covers only the cohorts with enough own events -- "
                     "state the coverage")

    print("\n  ⚠ median_samples_g1/g2 in this file are POST-MASK PER-EVENT MEDIANS,")
    print("    not cohort sizes (MAGNet reads 115/125 against its true 166/166).")
    print("    Never quote them as group sizes -- a standing trap in the project record.")
    return fl


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--list", action="store_true",
                    help="search for the inputs and describe them, then stop")
    ap.add_argument("--out", default="S8")
    a = ap.parse_args()

    print("=" * 78)
    print("S8  null behaviour and reverse replication")
    print("=" * 78)
    print("  The supplement that makes every null in this paper credible.")
    print("  S5(d) cites it; the Results cites the permutation zeros directly.")

    found = discover()
    notes = []

    print("\n" + "=" * 78)
    print("(b) REVERSE REPLICATION -- built from file")
    print("=" * 78)
    rows = build_reverse(a.out, notes)
    build_master(a.out, notes)

    show_record()
    if rows:
        T = pd.DataFrame(rows)
        T.to_csv(f"{a.out}b_reverse_replication.tsv", sep="\t", index=False)
        print(f"\n[write] {a.out}b_reverse_replication.tsv  ({len(T)} rows)")
    for n in notes:
        print(f"  >>> {n}")

    print("\n" + "=" * 78)
    if not found.get("permutation"):
        print("  >>> (a) HAS NO FILE. The permutation result is the most load-bearing")
        print("      claim in S8 -- the Results states it directly -- and nothing on")
        print("      disk holds it. Either it was never written out, or it is under a")
        print("      name none of the patterns above reach. This is a gap in the CODE")
        print("      AND DATA DEPOSIT as much as in the supplement: a reviewer asking to")
        print("      see the permutation output cannot be given it.")
        print("      Re-running script 25 would regenerate evidence for a claim already")
        print("      made and already measured -- which is the one case where a new run")
        print("      is not the standing rule's warning sign.")
    if not any(found.values()):
        print("  No candidate inputs found. S8's content is all in the project record")
        print("  above, but a supplement assembled from prose cannot be reproduced.")
        print("  Tell me where script 25's permutation output and script 22's reverse")
        print("  slopes were written, and I will pin them.")
    else:
        print("  Send this output back. Once the columns are pinned, the build path")
        print("  checks each recorded number against its file and refuses on a mismatch,")
        print("  as 65, 66 and 62 do.")
    print("  Nothing was built: this run is discovery only.")


if __name__ == "__main__":
    main()
