#!/usr/bin/env python3
"""
62_build_s5.py -- fill Supplementary Table S5 from the post-exclusion tables.

S5 has three parts and each asks the same question: does a pipeline choice change the
answer? They are computed here in one place so that all three are on the same event set
and the same exclusions, which is exactly what went wrong the first time.

  (a) JC versus JCEC          -- Deming transfer slopes under both quantifications
  (b) --readLength 50 vs 125  -- significant events per class, from the probe
  (c) exon-level vs full      -- replication concordance under both conventions

WHY THIS EXISTS RATHER THAN A RERUN OF SCRIPT 20
------------------------------------------------
`RESULT_magnitude_resolved.md` already holds a JC-versus-JCEC table, and it CONTRADICTS
the manuscript: it gives Brigham 0.757 where the Results and Figure 4 say 0.476.

⚠ CORRECTED 21 Oct. An earlier version of this paragraph said "the manuscript uses JC
on the post-exclusion `_dedup` tables". That is backwards. The manuscript is **JCEC**:
04, 05, 08 and 15b all run phase 2 on `<dir>/<ET>.MATS.JCEC.txt`, and 03/03b's use of
`.MATS.JC.txt` is a documented exception for the forward RBM20 target test alone. The
0.757 differs because that table is **pre-exclusion** and computed on the matched
JC-intersect-JCEC subset, not because of JC versus JCEC -- the same document records
that restricting to the matched subset alone raises Brigham from 0.662 to 0.757.
Putting it in the supplement would have published two different slopes for one cohort.

So part (a) is recomputed here, and the script REFUSES to emit it unless its own JC
column reproduces the manuscript's slopes. That check is the whole point: if the JC
column matches Figure 4, the JCEC column beside it was computed on the same footing and
the comparison is honest. If it does not match, nothing about the comparison is
trustworthy and the script says so instead of writing a table.

SELF-CHECKS, each of which aborts its own part rather than warning:
  (a) JCEC slopes must match 0.521 / 0.476 / 0.386 / 0.243 within 0.02
      (JCEC, not JC: the manuscript's tables are JCEC. The JC side comes
       from script 20's jccheck<SUF>/ output.)
      Each cohort emits TWO rows: own events, and the JC-intersect-JCEC
      matched subset. QUOTE THE MATCHED PAIR. The own-events pair is there
      to show how much of `change` is the event set rather than the
      estimator -- in script 20's Brigham background row that term is 0.09,
      larger than any JC-vs-JCEC difference in this panel.
      Each arm is a COMPLETE pipeline: the JC arm regresses the JC
      cohort table on the JC discovery (jccheck<SUF>/SE_v2filtered_JC),
      the JCEC arm on the JCEC discovery. The JCEC arm and the gate are
      unchanged, so 0.521 / 0.476 / 0.386 / 0.243 still has to hold.
  (b) the skipped-exon pair must be 172 (50 bp) and 167 (125 bp), i.e. -2.9%
  (c) exon-level concordance must match 79.6 / 81.2 / 69.2 / 76.1 within 0.5 points

Nothing here needs new rMATS. Every input already exists on disk.

  python 62_build_s5.py --list        # what it found for each part, then stop
  python 62_build_s5.py               # all three parts -> S5.tsv, S5.md
  python 62_build_s5.py --parts c     # just one part
"""

import argparse
import glob
import os
import re
import sys

import numpy as np
import pandas as pd

# -- conventions, copied from 45_cohort_slopes.py so the two cannot drift ------
NKEEP = {"SE": 3, "MXE": 5, "A3SS": 5, "A5SS": 5, "RI": 3}
TYPES = ["SE", "MXE", "A3SS", "A5SS", "RI"]
REPL = [("Colorado", "{ET}_Colorado_DCMvsNF"),
        ("Brigham", "Brigham_{ET}_HCM"),
        ("Fuwai", "Fuwai_{ET}_HCM"),
        ("Stanford", "Stanford_{ET}_HCM")]          # LBNL excluded: not independent
DROP_COMBOS = {("Colorado", "MXE")}
# Read length. DISC is a STEM: the readers below append
# "_dedup.all_events.tsv" or ".all_events.tsv" to it, and no _dedup variant
# exists at either read length, so this resolves to
# SE_v2filtered<SUF>.all_events.tsv. REPL below stays unsuffixed — no
# replication cohort was re-quantified.
SUF = os.environ.get("TBL_SUFFIX", "")
DISC = "{ET}_v2filtered" + SUF
# rMATS writes its own tables into each output directory. The probe directories hold
# raw rMATS output, not phase-2 `all_events` tables -- the first version of this script
# looked for the latter inside them and found nothing. Count from the raw files, which
# is what the recorded Stage A result (172 vs 167) was computed from.
COORD = {"SE": ["chr", "exonStart_0base", "exonEnd"],
         "MXE": ["chr", "1stExonStart_0base", "1stExonEnd",
                 "2ndExonStart_0base", "2ndExonEnd"],
         "A3SS": ["chr", "longExonStart_0base", "longExonEnd",
                  "shortES", "shortEE"],
         "A5SS": ["chr", "longExonStart_0base", "longExonEnd",
                  "shortES", "shortEE"],
         "RI": ["chr", "riExonStart_0base", "riExonEnd"]}


def count_rmats(path, et, fdr=0.05, dpsi=0.10):
    """Significant events in a raw rMATS MATS.JC.txt, collapsed to exon level."""
    t = pd.read_csv(path, sep="\t", low_memory=False)
    need = COORD[et] + ["FDR", "IncLevelDifference"]
    miss = [c for c in need if c not in t.columns]
    if miss:
        return None, f"missing columns {miss}"
    t["ekey"] = t[COORD[et]].astype(str).agg(":".join, axis=1)
    t = t.loc[t.groupby("ekey")["FDR"].idxmin()]
    d = pd.to_numeric(t["IncLevelDifference"], errors="coerce")
    return int(((t["FDR"] < fdr) & (d.abs() >= dpsi)).sum()), None

# the manuscript's own numbers, which these must reproduce
MANUSCRIPT_SLOPE = {"Colorado": 0.521, "Brigham": 0.476,
                    "Stanford": 0.386, "Fuwai": 0.243}
SLOPE_TOL = 0.02
MANUSCRIPT_CONC = {"Colorado": 79.6, "Brigham": 81.2, "Fuwai": 69.2, "Stanford": 76.1}
# the two conventions must give DIFFERENT discovery counts, or the switch did nothing
DISC_SE = {"exon": 926, "full": 1074}

# --------------------------------------------------------------------------
# 125 bp record, placed AFTER every literal it rewrites.
#
# Parts (a) and (c) ABORT rather than warn, so a 50 bp record against the
# 125 bp tables does not produce a wrong supplement — it refuses a correct
# one, reporting "JC slopes do not reproduce the manuscript", which reads as
# a pipeline fault. All three are therefore REQUIRED when TBL_SUFFIX is set.
#
# PROBE_SE, RECORDED_PCT and RECORDED_A3SS_BASE are deliberately NOT here:
# part (b) is already a 50-vs-125 comparison on the probe directories and has
# nothing to re-base.
_S5_REQ = ("S5_SLOPE", "S5_CONC", "S5_DISC_SE")
_s5env = {_k: os.environ.get(_k, "").strip() for _k in _S5_REQ}
if SUF:
    _miss = [_k for _k in _S5_REQ if not _s5env[_k]]
    if _miss:
        raise SystemExit(
            "REFUSING — TBL_SUFFIX=%s but these are unset: %s\n\n"
            "  Parts (a) and (c) abort rather than warn, so a 50 bp record\n"
            "  against the 125 bp tables refuses a CORRECT supplement.\n\n"
            "    %s\n" % (SUF, ", ".join(_miss), 'TBL_SUFFIX=_rl125 \\\n    S5_SLOPE=Colorado=0.5999,Brigham=0.5099,Stanford=0.4008,Fuwai=0.2948 \\\n    S5_CONC=Colorado=80.8,Brigham=81.9,Fuwai=70.8,Stanford=75.8 \\\n    S5_DISC_SE=exon=813,full=944 \\\n    python3 62_build_s5.py'))

if _s5env["S5_SLOPE"]:
    _d = dict(_p.split("=", 1) for _p in _s5env["S5_SLOPE"].split(","))
    assert set(_d) == set(MANUSCRIPT_SLOPE), (
        "S5_SLOPE must name exactly %s, got %s"
        % (sorted(MANUSCRIPT_SLOPE), sorted(_d)))
    MANUSCRIPT_SLOPE = {_k: float(_v) for _k, _v in _d.items()}

if _s5env["S5_CONC"]:
    _d = dict(_p.split("=", 1) for _p in _s5env["S5_CONC"].split(","))
    assert set(_d) == set(MANUSCRIPT_CONC), (
        "S5_CONC must name exactly %s, got %s"
        % (sorted(MANUSCRIPT_CONC), sorted(_d)))
    MANUSCRIPT_CONC = {_k: float(_v) for _k, _v in _d.items()}

if _s5env["S5_DISC_SE"]:
    _d = dict(_p.split("=", 1) for _p in _s5env["S5_DISC_SE"].split(","))
    assert set(_d) == set(DISC_SE), (
        "S5_DISC_SE must name exactly %s, got %s"
        % (sorted(DISC_SE), sorted(_d)))
    DISC_SE = {_k: int(_v) for _k, _v in _d.items()}
# --------------------------------------------------------------------------
CONC_TOL = 0.5
PROBE_SE = (172, 167)          # significant SE events at 50 bp and 125 bp
# The PERCENTAGES are the verified quantities in the project record; the raw counts are
# what S5(b) is missing. So check every class against its recorded percentage, not only
# the one class whose counts happen to be written down.
RECORDED_PCT = {"SE": -2.9, "MXE": 0.0, "RI": 0.0, "A5SS": 0.0, "A3SS": +11.1}
RECORDED_A3SS_BASE = 36        # "+11.1% on 36 events" -- the base is part of the claim
PCT_TOL = 0.15                 # recorded to one decimal place


def read_events(path, et, exon_level=True):
    """Load an all_events table and collapse by the chosen convention.

    exon_level=True  -> truncate `key` to NKEEP[et] fields (the paper's convention)
    exon_level=False -> keep the full coordinate tuple (the alternative)
    """
    if not (os.path.exists(path) and os.path.getsize(path)):
        return None, "missing"
    t = pd.read_csv(path, sep="\t", low_memory=False)
    if "key" not in t.columns or "FDR" not in t.columns:
        return None, f"no key/FDR column; have {list(t.columns)[:8]}"
    d = next((c for c in ("dPSI_adj", "dPSI") if c in t.columns), None)
    if d is None:
        cand = [c for c in t.columns if c.lower().startswith("dpsi")]
        d = next((c for c in cand if "adj" in c.lower()), cand[0] if cand else None)
    if d is None:
        return None, "no dPSI column"
    if "se_adj" not in t.columns and "se" in t.columns:
        t = t.rename(columns={"se": "se_adj"})
    t["ekey"] = ([":".join(str(k).split(":")[:NKEEP[et]]) for k in t["key"]]
                 if exon_level else t["key"].astype(str))
    t = t.loc[t.groupby("ekey")["FDR"].idxmin()]
    keep = ["ekey", "FDR"] + [c for c in ("se_adj",) if c in t.columns]
    return t.rename(columns={d: "dpsi"})[["ekey", "dpsi", "FDR"] +
                                         (["se_adj"] if "se_adj" in t.columns else [])], None


def first_existing(*paths):
    for p in paths:
        if os.path.exists(p) and os.path.getsize(p):
            return p
    return None


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


# ============================================================ part (a)
def part_a(a, rows, notes):
    print("\n" + "=" * 78)
    print("(a) JC versus JCEC")
    print("=" * 78)
    et = "SE"
    dpath = first_existing(f"{DISC.replace('{ET}', et)}_dedup.all_events.tsv",
                           f"{DISC.replace('{ET}', et)}.all_events.tsv")
    if dpath is None:
        notes.append("(a) skipped: no discovery table")
        print("  [fatal for (a)] no discovery table found"); return
    D, err = read_events(dpath, et)
    if D is None:
        notes.append(f"(a) skipped: discovery {err}"); print(f"  [fatal for (a)] {err}"); return
    sig = D[(D.FDR < 0.05) & (D.dpsi.abs() >= 0.10)]
    print(f"  discovery {dpath}: {len(sig)} significant of {len(D)}")

    # ---- THE JC DISCOVERY, added 21 Oct ------------------------------
    # Until now BOTH arms regressed against `sig`, the JCEC discovery, so
    # the JC column was a hybrid: JCEC on the x-axis, JC on the y-axis. It
    # measured the replication-side quantification only. "JC versus JCEC"
    # means two complete pipelines, so the JC arm takes the JC discovery
    # script 20 already writes -- and which 20 itself pairs with the JCEC
    # discovery at its own line 221.
    #
    # The significance rule is applied IDENTICALLY (FDR < 0.05 and
    # |dPSI| >= 0.10 on the collapsed table), so the two discovery sets
    # differ only because the quantification differs. That difference is
    # PART OF the JC-vs-JCEC effect and is deliberately not matched away
    # here; the matched arm below controls for it separately.
    jdpath = first_existing(
        f"jccheck{SUF}/{et}_v2filtered_JC_dedup.all_events.tsv",
        f"jccheck{SUF}/{et}_v2filtered_JC.all_events.tsv")
    SIG = {"JCEC": sig}
    if jdpath is None:
        print(f"  \u26a0 NO JC DISCOVERY in jccheck{SUF}/ -- the JC arm falls")
        print("    back to the JCEC discovery, making this a REPLICATION-SIDE-")
        print("    ONLY comparison. Either say so in the caption or run")
        print(f"      TBL_SUFFIX={SUF} sbatch 20_jc_vs_jcec_readlength.sbatch")
        notes.append("(a) JC arm fell back to the JCEC discovery -- "
                     "replication-side only, NOT a pipeline comparison")
        SIG["JC"] = sig
    else:
        JD, jerr = read_events(jdpath, et)
        if JD is None:
            print(f"  [fatal for (a)] JC discovery: {jerr}")
            notes.append(f"(a) skipped: JC discovery {jerr}"); return
        SIG["JC"] = JD[(JD.FDR < 0.05) & (JD.dpsi.abs() >= 0.10)]
        print(f"  discovery {jdpath}: {len(SIG['JC'])} significant of {len(JD)}")
        print(f"  note: DISC_SE['exon'] = {DISC_SE['exon']} is the JCEC")
        print(f"        discovery count; the JC discovery has "
              f"{len(SIG['JC'])}. The two are not expected to agree and no")
        print("        gate is applied to the JC one.")

    out = {}
    for coh, pat in REPL:
        base = pat.replace("{ET}", et)
        # THE LABELS WERE SWAPPED UNTIL 21 Oct. The plain/_dedup cohort
        # tables come from phase 2 on rmats_<COH>/<ET>.MATS.JCEC.txt -- they
        # ARE the JCEC side, and they are what 45 computed the manuscript's
        # slopes from. The JC side is script 20's output in jccheck<SUF>/,
        # whose names contain "_JC" and never "JCEC", which is why the old
        # glob could never match and (a) never built.
        #
        # Both sides prefer a _dedup variant where one exists, so they sit on
        # the same exclusion basis. Only Stanford has one: exclude_dups.ids
        # contains exactly one sample from any of these three cohorts.
        jcec = first_existing(f"{base}_dedup.all_events.tsv",
                              f"{base}.all_events.tsv")
        jc = first_existing(
            f"jccheck{SUF}/{coh}_{et}_HCM_JC_dedup.all_events.tsv",
            f"jccheck{SUF}/{coh}_{et}_HCM_JC.all_events.tsv")
        print(f"  {coh:<10} JCEC -> {jcec or '*** MISSING ***'}")
        print(f"  {'':<10} JC   -> {jc or '*** MISSING (script 20 writes these) ***'}")
        if jcec and jc and ("dedup" in jcec) != ("dedup" in jc):
            print(f"  {'':<10}      ⚠ ONE SIDE IS _dedup AND THE OTHER IS NOT --")
            print(f"  {'':<10}        the comparison would mix exclusion bases")
        out[coh] = (jc, jcec)

    if a.list:
        return
    have = {c: v for c, v in out.items() if v[0] and v[1]}
    gone = [c for c in out if c not in have]
    if not have:
        print("\n  >>> (a) NOT BUILT. No cohort has BOTH sides.")
        print("      The JC side comes from script 20, which writes")
        print(f"      jccheck{SUF}/<COH>_SE_HCM_JC.all_events.tsv. Run:")
        print(f"        TBL_SUFFIX={SUF} sbatch 20_jc_vs_jcec_readlength.sbatch")
        print("      The JCEC side is the plain/_dedup cohort table and should")
        print("      already exist. Do NOT re-quantify: raw SE.MATS.JCEC.txt is")
        print("      not what is needed -- Deming needs the per-event SE that")
        print("      phase 2 adds, and those tables are already on disk.")
        notes.append("(a) NOT BUILT -- no cohort has both sides")
        return
    if gone:
        print(f"\n  ⚠ SCOPE: no JC table for {', '.join(gone)}.")
        print("    Script 20's part C covers Brigham, Fuwai and Stanford only;")
        print("    Colorado is the DCM-vs-NF contrast and sits outside its HCM")
        print("    loop. (a) therefore covers "
              f"{len(have)} of {len(out)} cohorts.")
        print("    THE CAPTION MUST SAY SO. 62's own rule is that a supplement")
        print("    must not carry a blank row, so these cohorts are omitted")
        print("    rather than left empty.")
        notes.append(f"(a) covers {len(have)}/{len(out)} cohorts -- no JC table "
                     f"for {', '.join(gone)}")

    res = {}
    for coh, (jc, jcec) in have.items():
        r = {}
        M = {}
        for lab, p in (("JC", jc), ("JCEC", jcec)):
            R, err = read_events(p, et)
            if R is None:
                r[lab] = np.nan; continue
            # each basis against ITS OWN discovery. SIG["JCEC"] IS `sig`,
            # so the JCEC arm and the gate below are bit-for-bit unchanged.
            m = SIG[lab].merge(R, on="ekey", suffixes=("", "_r"))
            if "se_adj" not in m or "se_adj_r" not in m:
                r[lab] = np.nan; continue
            m = m[np.isfinite(m.dpsi.values) & np.isfinite(m.dpsi_r.values)
                  & np.isfinite(m.se_adj.values) & np.isfinite(m.se_adj_r.values)]
            M[lab] = m
            r["n_" + lab] = int(len(m))
            r[lab] = deming(m.dpsi.values, m.dpsi_r.values,
                            m.se_adj.values, m.se_adj_r.values)
        # ---- THE MATCHED ARM, added 21 Oct -----------------------------
        # Fitting each basis on its own events conflates the estimator
        # (JC vs JCEC) with the event set (the two tables differ in size by
        # 20-30%). Script 20 measures the second term directly: Brigham's
        # background slope is 0.662 on its own events and 0.757 on the
        # JC-intersect-JCEC subset -- a 0.09 gap, larger than any
        # JC-vs-JCEC difference in this panel. Without these columns a
        # reader cannot tell which term produced `change`.
        if len(M) == 2:
            common = set(M["JC"]["ekey"]) & set(M["JCEC"]["ekey"])
            for lab in ("JC", "JCEC"):
                mm = M[lab][M[lab]["ekey"].isin(common)]
                r[lab + "_m"] = deming(mm.dpsi.values, mm.dpsi_r.values,
                                       mm.se_adj.values, mm.se_adj_r.values)
            r["n_m"] = int(len(common))
        res[coh] = r

    # THE CHECK: the JCEC column must be the manuscript's slopes.
    # 45 computed 0.521 / 0.476 / 0.386 / 0.243 from the _dedup JCEC tables,
    # so JCEC is the column that has to reproduce them. Until 21 Oct this
    # gate read res[c]["JC"] -- and would have PASSED, because the column
    # then called "JC" was in fact the JCEC table. A gate that passes on a
    # mislabelled column certifies the mistake.
    bad = [(c, res[c]["JCEC"], MANUSCRIPT_SLOPE[c]) for c in res
           if np.isfinite(res[c]["JCEC"])
           and abs(res[c]["JCEC"] - MANUSCRIPT_SLOPE[c]) > SLOPE_TOL]
    print()
    print(f"  {'cohort':<10}{'JCEC':>9}{'manuscript':>12}{'JC':>9}{'change':>9}"
          f"{'|':>4}{'JCEC_m':>9}{'JC_m':>9}{'change_m':>10}{'n_m':>7}")
    for c in ("Colorado", "Brigham", "Fuwai", "Stanford"):
        r = res.get(c, {})
        f = lambda v: f"{v:.3f}" if np.isfinite(v) else "-"
        g = lambda v: f"{v:+.3f}" if np.isfinite(v) else "-"
        ch = (r.get("JC", np.nan) - r.get("JCEC", np.nan))
        chm = (r.get("JC_m", np.nan) - r.get("JCEC_m", np.nan))
        _n = r.get("n_m")
        print(f"  {c:<10}{f(r.get('JCEC', np.nan)):>9}{MANUSCRIPT_SLOPE[c]:>12.3f}"
              f"{f(r.get('JC', np.nan)):>9}{g(ch):>9}"
              f"{'|':>4}{f(r.get('JCEC_m', np.nan)):>9}"
              f"{f(r.get('JC_m', np.nan)):>9}{g(chm):>10}"
              f"{(str(_n) if _n else '-'):>7}")
    if bad:
        print("\n  >>> (a) NOT BUILT. The recomputed JCEC slopes do not reproduce")
        print("      the manuscript's:")
        for c, got, want in bad:
            print(f"        {c}: {got:.3f} here, {want:.3f} in the Results")
        print("      Until the JCEC column matches Figure 4, the JC column beside it")
        print("      is not computed on the same footing and the comparison means")
        print("      nothing.")
        notes.append("(a) NOT BUILT -- JCEC column does not reproduce the manuscript")
        return
    print("\n  [ok] the JCEC column reproduces the manuscript's slopes")
    print("  ⚠ EACH ARM IS A COMPLETE PIPELINE: JC against the JC discovery,")
    print("    JCEC against the JCEC discovery. The own-events pair therefore")
    print("    differs in BOTH the estimator and the event set -- the two")
    print("    discovery sets and the two cohort tables all differ in size.")
    print("    `change_m` is the like-for-like figure: same events, both ")
    print("    bases. QUOTE change_m; the own-events pair beside it shows")
    print("    how large the event-set term is.")
    print("  ⚠ DO NOT COMPARE THESE ROWS WITH script 20's jc_vs_jcec<SUF>.tsv")
    print("    ROW-FOR-ROW. Script 20 splits the DCM-significant events into")
    print("    `program` and `background` and reports each separately; this")
    print("    panel uses the manuscript's full significant set, which is why")
    print("    its Brigham JCEC is 0.476 and 20's background JCEC is 0.662.")
    print("    Script 20's MATCHED BACKGROUND arms are the source for the")
    print("    READ-LENGTH sentence: internal -0.006 vs externals ~+0.04 at")
    print("    50 bp, but internal -0.021 vs +0.026 / -0.007 / +0.009 at")
    print("    125 bp, which STRADDLE it -- so 'read length is largely")
    print("    excluded' does not survive and needs rewriting.")
    print("    THIS panel is the source for the JC-vs-JCEC sentence only.")
    for c in sorted(res):
        rows.append(dict(part="a. JC vs JCEC", row=c,
                         value_1=round(res[c]["JCEC"], 3), label_1="Deming, JCEC",
                         value_2=round(res[c]["JC"], 3), label_2="Deming, JC",
                         change=round(res[c]["JC"] - res[c]["JCEC"], 3)))
        # the matched pair as a SECOND ROW, not new columns: parts (b) and
        # (c) write the same seven-key schema into the same table.
        if (np.isfinite(res[c].get("JCEC_m", np.nan))
                and np.isfinite(res[c].get("JC_m", np.nan))):
            rows.append(dict(
                part="a. JC vs JCEC",
                row=f"{c}, matched subset (n={res[c].get('n_m', 0)})",
                value_1=round(res[c]["JCEC_m"], 3),
                label_1="Deming, JCEC (matched)",
                value_2=round(res[c]["JC_m"], 3),
                label_2="Deming, JC (matched)",
                change=round(res[c]["JC_m"] - res[c]["JCEC_m"], 3)))
        else:
            notes.append(f"(a) {c}: NO MATCHED ARM -- one basis unreadable, "
                         f"so its `change` mixes estimator and event set")

    # ---- THE COUNT ROWS, merged 21 Oct --------------------------------
    # 63_s5a_jc_vs_jcec.py measures what the slopes above do not: how MANY
    # events each quantification calls. The two halves must travel
    # together. RESULT_s5a_jc_vs_jcec.md's own warning is that S5(a)
    # "cannot be summarised as 'the choice does not change the answer'" --
    # and a slopes-only panel says precisely that, because the matched
    # slopes agree to within 0.024 while the counts differ by 9% to 58%.
    _s5a = f"S5a{SUF}.tsv"
    if not os.path.exists(_s5a):
        print(f"\n  ⚠ {_s5a} NOT FOUND. S5(a) will carry slopes only, which")
        print("    states only the reassuring half of the comparison. Build")
        print("    the counts with:")
        print(f"      TBL_SUFFIX={SUF} python3 63_s5a_jc_vs_jcec.py")
        notes.append(f"(a) COUNTS MISSING: {_s5a} not found, so S5(a) reports "
                     f"slopes without the count change")
    else:
        _C = pd.read_csv(_s5a, sep="\t")
        # BASIS GUARD. After 148 a re-based run of 63 tags every arm name
        # with "@". An S5a table without that tag predates 148 and is a
        # 50 bp measurement, so merging it into a 125 bp S5 would put right
        # numbers under a wrong name -- this project's recurring defect.
        if SUF and not _C["arm"].astype(str).str.contains("@").all():
            print(f"\n  >>> REFUSING to merge {_s5a}: its arm names carry no")
            print("      read-length basis, so it predates patch 148 and is a")
            print("      50 bp table. Re-run 63 with TBL_SUFFIX set.")
            notes.append(f"(a) COUNTS REFUSED: {_s5a} carries no per-arm "
                         f"read-length basis")
        else:
            for _, _r in _C.iterrows():
                rows.append(dict(
                    part="a. JC vs JCEC (counts)", row=str(_r["arm"]),
                    value_1=f"{_r['sig_JC']} significant", label_1="JC",
                    value_2=f"{_r['sig_JCEC']} significant", label_2="JCEC",
                    change=str(_r["sig_change"])))
                rows.append(dict(
                    part="a. JC vs JCEC (counts)",
                    row=f"{_r['arm']}, direction",
                    value_1=str(_r["dir_agree_sig"]),
                    label_1="same direction, JC-significant",
                    value_2=str(_r["rho_sig"]), label_2="Spearman rho",
                    change="-"))
            print(f"\n  [ok] merged {len(_C)} count rows from {_s5a}")
            print("    S5(a) now states BOTH halves: the number of events")
            print("    called is strongly quantification-dependent, and the")
            print("    direction and magnitude of those events are not.")


# ============================================================ part (b)
def part_b(a, rows, notes):
    print("\n" + "=" * 78)
    print("(b) declared read length, 50 bp versus 125 bp")
    print("=" * 78)
    dirs = sorted({d for d in glob.glob("*") if os.path.isdir(d)
                   and re.search(r"probe|rl\d{2,3}|readlen", d, re.I)})
    print(f"  candidate probe directories: {', '.join(dirs) or 'NONE FOUND'}")

    # ---- DIRECTORY SELECTION, rewritten 21 Oct ------------------------
    # The old picker took the FIRST sorted name containing "125". On
    # 21 Oct the 125 bp run of script 20 created jccheck_rl125/, which
    # sorts ahead of probe_rl125/ and holds all_events tables rather than
    # rMATS output. (b) had not changed; its directory had. A NEW name
    # silently satisfying an OLD glob is the mirror of Rule 8.
    #
    # Selection is now by CONTENT, and every rejection is printed:
    #   1. the directory must contain an SE probe file
    #   2. scratch and abandoned copies are dropped
    #   3. a name beginning "probe" is preferred
    #   4. anything still ambiguous REFUSES and asks for --probeNN
    def _probe_file(_d):
        return first_existing(os.path.join(_d, "SE.MATS.JC.txt"),
                              os.path.join(_d, "SE.MATS.JCEC.txt"),
                              os.path.join(_d, "SE.all_events.tsv"),
                              os.path.join(_d, "SE_probe.all_events.tsv"))

    def _pick(_want, _flag):
        _usable, _why = [], []
        for _d in [x for x in dirs if _want(x)]:
            if _probe_file(_d) is None:
                _why.append(f"{_d}: holds no SE probe file"); continue
            if re.search(r"tmp|abandon", _d, re.I):
                _why.append(f"{_d}: scratch or abandoned copy"); continue
            _usable.append(_d)
        for _w in _why:
            print(f"    rejected {_w}")
        _pref = [x for x in _usable if x.lower().startswith("probe")]
        _chosen = _pref or _usable
        if len(_chosen) > 1:
            print(f"    AMBIGUOUS, {len(_chosen)} usable: "
                  f"{', '.join(_chosen)}")
            print(f"    pass {_flag} DIR -- refusing to guess")
            return None
        return _chosen[0] if _chosen else None

    # an explicit directory always wins over the guess -- the help text promised this
    d50 = a.probe50 or _pick(
        lambda _d: re.search(r"(^|[^0-9])50([^0-9]|$)", _d), "--probe50")
    d125 = a.probe125 or _pick(lambda _d: "125" in _d, "--probe125")
    print(f"  50 bp  -> {d50 or '*** not identified ***'}")
    print(f"  125 bp -> {d125 or '*** not identified ***'}")
    if a.list:
        return
    for lab, d in (("--probe50", a.probe50), ("--probe125", a.probe125)):
        if d and not os.path.isdir(d):
            print(f"\n  >>> (b) NOT BUILT. {lab} {d} is not a directory.")
            notes.append(f"(b) NOT BUILT -- {lab} {d} is not a directory")
            return
    if not (d50 and d125):
        print("\n  >>> (b) NOT BUILT. Could not identify the two probe directories.")
        print("      Pass them explicitly: --probe50 DIR --probe125 DIR")
        notes.append("(b) NOT BUILT -- probe directories not identified")
        return

    counts = {}
    for et in TYPES:
        pair = []
        for d in (d50, d125):
            # raw rMATS first -- that is what these directories contain
            p = first_existing(os.path.join(d, f"{et}.MATS.JC.txt"),
                               os.path.join(d, f"{et}.MATS.JCEC.txt"),
                               os.path.join(d, f"{et}.all_events.tsv"),
                               os.path.join(d, f"{et}_probe.all_events.tsv"))
            if p is None:
                pair.append(None); continue
            if p.endswith(".MATS.JC.txt") or p.endswith(".MATS.JCEC.txt"):
                n, err = count_rmats(p, et)
                if err:
                    print(f"  [warn] {p}: {err}")
                pair.append(n)
            else:
                T, err = read_events(p, et)
                pair.append(None if T is None
                            else int(((T.FDR < 0.05) & (T.dpsi.abs() >= 0.10)).sum()))
        counts[et] = pair
        s = " / ".join("-" if v is None else str(v) for v in pair)
        print(f"  {et:<6} {s}")

    se = counts.get("SE", [None, None])
    if se[0] != PROBE_SE[0] or se[1] != PROBE_SE[1]:
        print(f"\n  >>> (b) NOT BUILT. Skipped exons give {se}, the recorded probe result")
        print(f"      is {list(PROBE_SE)} (-2.9%). Wrong directories, or a different")
        print("      filter. Resolve before using any row of (b).")
        notes.append("(b) NOT BUILT -- SE counts do not match the recorded probe")
        return
    print("\n  [ok] skipped exons reproduce the recorded 172 / 167")
    # provenance: record WHICH probe the built rows came from, so
    # it never has to be re-derived from a directory listing
    notes.append(f"(b) probe directories: 50 bp {d50}, 125 bp {d125}")
    # SCOPE. Added 21 Oct. (b) is a PROBE result and will be read
    # as a cohort result unless it says otherwise.
    print("\n  ⚠ SCOPE: these are PROBE counts, from probe_rl50/ and")
    print("    probe_rl125/ — roughly a fifth of the discovery set. The")
    print("    cohort-level change is FOUR TIMES larger:")
    print(f"      probe, SE            172 -> 167     -2.9%")
    print(f"      discovery, SE exon   {DISC_SE['exon']:>4} -> ???    see the 125 bp record")
    print(f"      discovery, SE full   {DISC_SE['full']:>4} -> ???    (script 21 agrees on both)")
    print("    They do not conflict — they are different event sets. But a")
    print("    sentence citing (b) for how much read length matters is")
    print("    generalising a probe to the cohort, and the full")
    print("    re-quantification contradicts that by a factor of four.")

    # Second check, and the stronger one. The counts are what (b) is missing; the
    # per-class PERCENTAGES are already on record and verified. If a recovered pair
    # gives a different percentage, the pair is wrong even when SE happens to match.
    pcts, drift = {}, []
    for et in TYPES:
        a50, a125 = counts[et]
        if a50 is None or a125 is None:
            drift.append((et, "one arm missing")); continue
        if not a50:
            drift.append((et, "zero events at 50 bp -- percentage undefined")); continue
        pcts[et] = 100 * (a125 - a50) / a50
        want = RECORDED_PCT.get(et)
        if want is not None and abs(pcts[et] - want) > PCT_TOL:
            drift.append((et, f"{pcts[et]:+.1f}% here, {want:+.1f}% on record"))
    if "A3SS" in counts and counts["A3SS"][0] not in (None, RECORDED_A3SS_BASE):
        drift.append(("A3SS", f"base {counts['A3SS'][0]} events, "
                              f"{RECORDED_A3SS_BASE} on record"))
    if drift:
        print("\n  >>> (b) NOT BUILT. The recovered counts disagree with the recorded")
        print("      per-class changes, which were verified from Stage A output:")
        for et, why in drift:
            print(f"        {et:<5} {why}")
        print("      Skipped exons matching is not sufficient -- these are the same")
        print("      measurement reported two ways, so they must agree. Resolve before")
        print("      using any row of (b).")
        notes.append("(b) NOT BUILT -- recovered counts disagree with recorded percentages")
        return
    print("  [ok] all five classes reproduce their recorded percentage change")

    for et in TYPES:
        if et not in pcts:
            continue
        a50, a125 = counts[et]
        rows.append(dict(part="b. read length", row=et,
                         value_1=a50, label_1="significant, 50 bp",
                         value_2=a125, label_2="significant, 125 bp",
                         change=f"{pcts[et]:+.1f}%"))


# ============================================================ part (c)
def part_c(a, rows, notes):
    print("\n" + "=" * 78)
    print("(c) exon-level versus full-coordinate counting")
    print("=" * 78)
    if a.list:
        print("  needs only the discovery and replication tables already used by 61")
        return
    et = "SE"
    out, disc_counts = {}, {}
    for exon_level, lab in ((True, "exon"), (False, "full")):
        dpath = first_existing(f"{DISC.replace('{ET}', et)}_dedup.all_events.tsv",
                               f"{DISC.replace('{ET}', et)}.all_events.tsv")
        D, err = read_events(dpath, et, exon_level)
        if D is None:
            notes.append(f"(c) skipped: discovery {err}"); print(f"  [fatal] {err}"); return
        sig = D[(D.FDR < 0.05) & (D.dpsi.abs() >= 0.10)]
        for coh, pat in REPL:
            base = pat.replace("{ET}", et)
            p = first_existing(f"{base}_dedup.all_events.tsv", f"{base}.all_events.tsv")
            if p is None:
                continue
            R, err = read_events(p, et, exon_level)
            if R is None:
                continue
            m = sig.merge(R, on="ekey", suffixes=("", "_r"))
            same = int((np.sign(m.dpsi) == np.sign(m.dpsi_r)).sum())
            out[(coh, lab)] = (same, len(m), 100 * same / len(m) if len(m) else np.nan)
        want = DISC_SE[lab]
        ok = abs(len(sig) - want) <= 2
        print(f"  {lab:<5} convention: discovery {len(sig)} significant "
              f"(expect {want}) {'ok' if ok else '*** MISMATCH ***'}")
        if not ok:
            notes.append(f"(c) {lab}-convention discovery count {len(sig)}, expect {want}")
        disc_counts[lab] = len(sig)

    if disc_counts.get("exon") == disc_counts.get("full"):
        print("\n  >>> (c) NOT BUILT. Both conventions returned the same discovery count")
        print(f"      ({disc_counts.get('exon')}). The truncation collapsed nothing, so")
        print("      the `key` column is not the multi-field coordinate tuple this")
        print("      convention assumes. Check the table before reading anything below.")
        notes.append("(c) NOT BUILT -- the two conventions gave identical counts")
        return

    bad = []
    print(f"\n  {'cohort':<10}{'exon-level':>22}{'full-coordinate':>22}{'change':>9}")
    for coh, _ in REPL:
        e = out.get((coh, "exon")); f_ = out.get((coh, "full"))
        if not e or not f_:
            continue
        print(f"  {coh:<10}{f'{e[0]}/{e[1]} ({e[2]:.1f}%)':>22}"
              f"{f'{f_[0]}/{f_[1]} ({f_[2]:.1f}%)':>22}{f_[2]-e[2]:>+9.1f}")
        if abs(e[2] - MANUSCRIPT_CONC[coh]) > CONC_TOL:
            bad.append((coh, e[2], MANUSCRIPT_CONC[coh]))
        rows.append(dict(part="c. counting convention", row=coh,
                         value_1=f"{e[0]}/{e[1]} ({e[2]:.1f}%)", label_1="exon-level",
                         value_2=f"{f_[0]}/{f_[1]} ({f_[2]:.1f}%)",
                         label_2="full-coordinate",
                         change=f"{f_[2]-e[2]:+.1f} pts"))
    if bad:
        print("\n  >>> (c) exon-level column does NOT reproduce the published concordance:")
        for c, got, want in bad:
            print(f"        {c}: {got:.1f}% here, {want}% published")
        print("      61_table2_lead_genes.py reproduced these exactly on 5 Oct, so a")
        print("      mismatch here means this script is reading different tables.")
        notes.append("(c) exon-level column does not reproduce the published figures")
        for r in [r for r in rows if r["part"].startswith("c.")]:
            rows.remove(r)
    else:
        print("\n  [ok] exon-level column reproduces the published concordance")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--parts", default="abc")
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--probe50"); ap.add_argument("--probe125")
    ap.add_argument("--out", default="S5" + SUF)
    a = ap.parse_args()
    rows, notes = [], []
    if "a" in a.parts: part_a(a, rows, notes)
    if "b" in a.parts: part_b(a, rows, notes)
    if "c" in a.parts: part_c(a, rows, notes)
    if a.list:
        print("\n[--list] stopping here."); return

    print("\n" + "=" * 78)
    print("RESULT")
    print("=" * 78)
    if rows:
        T = pd.DataFrame(rows)
        T.to_csv(f"{a.out}.tsv", sep="\t", index=False)
        print(T.to_string(index=False))
        print(f"\n[write] {a.out}.tsv")
    else:
        print("  nothing built")
    for n in notes:
        print(f"  >>> {n}")
    built = {r["part"][0] for r in rows}
    missing = sorted(set("abc") - built)
    if missing:
        print(f"\n  S5 is INCOMPLETE: part(s) {', '.join(missing)} not built.")
        print("  Do not submit a supplement with a blank row; either finish the part or")
        print("  drop it from S5 and say in the Methods that it was not assessed.")
    else:
        print("\n  All three parts built and all three self-checks passed.")


if __name__ == "__main__":
    main()
