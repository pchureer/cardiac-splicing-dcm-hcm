#!/usr/bin/env python3
"""
33_rebuild_summaries_dedup.py -- regenerate every downstream summary from the
post-exclusion tables, without re-running the sbatch scripts.

WHY NOT JUST RERUN 12 AND 15b
------------------------------
Both still point at the original bam lists, so they would recompute the PRE-exclusion
numbers and quietly overwrite the correct ones. Editing scripts you have already
validated is the other option and is worse. This reads the *_dedup.all_events.tsv
tables script 32 produced and derives everything from them.

For every table it prints WHICH FILE it used. If a _dedup table is missing it falls
back to the original and says so loudly, because a silent mix of pre- and
post-exclusion inputs is exactly the failure this script exists to prevent.

WHAT IT PRODUCES
  A. Colorado internal control at three effect thresholds -- the F3 panel A numbers,
     including the |dPSI|>=0.10 counts that script 32 did not report.
  B. Equivalence ratios, median |dPSI(DCM-ICM)| / |dPSI(DCM-NF)|.
  C. DCM-vs-NF against ICM-vs-NF: sign concordance and correlation, with the
     shared-control baseline stated.
  D. Direction concordance of MAGNet-significant SE events in each replication cohort,
     and the POOLED total WITHOUT LBNL. The published 1,960/2,574 = 76.1% included
     LBNL, whose samples are MAGNet hearts, so that figure double-counts the discovery
     cohort and is inflated. This recomputes it honestly.
  E. Transfer slopes (Deming, errors-in-variables) per cohort.

EVERYTHING IS EXON-LEVEL, matching the locked convention (script 24): the key is
collapsed to the first NKEEP fields and the lowest-FDR row per collapsed key is kept.

Usage:
  python 33_rebuild_summaries_dedup.py
  python 33_rebuild_summaries_dedup.py --fdr 0.05 --out summaries_dedup.txt
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

TYPES = ["SE", "MXE", "A3SS", "A5SS", "RI"]
NKEEP = {"SE": 3, "MXE": 5, "A3SS": 5, "A5SS": 5, "RI": 3}

# cohort -> (filename pattern, is a DCM-vs-NF contrast, excluded-from-pooling reason)
REPL = [
    ("Colorado", "{ET}_Colorado_DCMvsNF", True, None),
    ("LBNL", "{ET}_LBNL_DCMvsNF", True,
     "24 of 33 samples are MAGNet hearts -- pooling double-counts discovery"),
    ("Brigham", "Brigham_{ET}_HCM", False, None),
    ("Fuwai", "Fuwai_{ET}_HCM", False, None),
    ("Stanford", "Stanford_{ET}_HCM", False, None),
]
# event types whose quantification failed reverse replication (script 26)
BAD = {("Colorado", "MXE"), ("LBNL", "MXE")}


def load(base, et, prefer_dedup=True):
    """Return (DataFrame, path_used). Prefers the _dedup table; says so if it falls back."""
    cands = ([f"{base}_dedup.all_events.tsv", f"{base}.all_events.tsv"] if prefer_dedup
             else [f"{base}.all_events.tsv"])
    for p in cands:
        if os.path.exists(p) and os.path.getsize(p):
            t = pd.read_csv(p, sep="\t")
            if "key" not in t.columns or "FDR" not in t.columns:
                return None, f"{p} (UNUSABLE: no key/FDR column)"
            n = NKEEP[et]
            t["ekey"] = [":".join(str(k).split(":")[:n]) for k in t["key"]]
            t = t.loc[t.groupby("ekey")["FDR"].idxmin()].copy()
            return t, p
    return None, f"{base}[_dedup].all_events.tsv NOT FOUND"


def deming(x, y, sx, sy):
    """Errors-in-variables slope; lambda from the reported standard errors."""
    ok = np.isfinite(x) & np.isfinite(y) & np.isfinite(sx) & np.isfinite(sy)
    x, y, sx, sy = x[ok], y[ok], sx[ok], sy[ok]
    if len(x) < 5:
        return np.nan, len(x)
    lam = (np.mean(sy ** 2) / np.mean(sx ** 2)) if np.mean(sx ** 2) > 0 else 1.0
    dx, dy = x - x.mean(), y - y.mean()
    Sxx, Syy, Sxy = np.sum(dx ** 2), np.sum(dy ** 2), np.sum(dx * dy)
    if Sxy == 0:
        return np.nan, len(x)
    disc = (Syy - lam * Sxx) ** 2 + 4 * lam * Sxy ** 2
    return float((Syy - lam * Sxx + np.sqrt(disc)) / (2 * Sxy)), len(x)


# --------------------------------------------------------------------------
# READ LENGTH.  33 writes two files whose names were string literals, so its
# 125 bp run overwrote the 50 bp `direction_concordance_dedup.tsv` (read by 38,
# 44, 84 and 85) and `colorado_internal_control_dedup.tsv` in place, and exited
# 0.  `--out` was the only output that carried a suffix, because it was the
# only one that was a parameter.
#
# SUF is read here; _resolve_suffix() below refuses every combination in which
# --discovery, --out and SUF could disagree.
SUF = os.environ.get("TBL_SUFFIX", "")


def _looks_readlength(stem):
    """True if `stem` ends in a read-length suffix such as `_rl125`."""
    i = stem.rfind("_rl")
    return i >= 0 and stem[i + 3:].isdigit() and len(stem[i + 3:]) > 0


def _resolve_suffix(a):
    """Refuse any disagreement between --discovery, --out and TBL_SUFFIX.

    Returns `a` with a.out resolved.  Called immediately after parse_args, so
    nothing is opened before the check.
    """
    d = str(a.discovery)
    if SUF:
        if not d.endswith(SUF):
            raise SystemExit(
                "REFUSING — TBL_SUFFIX is %r but --discovery is %r, which does\n"
                "  not end with it. One of the two is wrong, and 33 writes two\n"
                "  files whose names follow TBL_SUFFIX. Set them consistently:\n"
                "    TBL_SUFFIX=%s python3 33_rebuild_summaries_dedup.py \\\n"
                "        --discovery '{ET}_v2filtered%s'" % (SUF, d, SUF, SUF))
    else:
        if _looks_readlength(d):
            raise SystemExit(
                "REFUSING — --discovery is %r, which carries a read-length\n"
                "  suffix, but TBL_SUFFIX is unset. That is exactly the run that\n"
                "  overwrote the 50 bp direction_concordance_dedup.tsv and\n"
                "  colorado_internal_control_dedup.tsv with 125 bp content under\n"
                "  50 bp names, and reported success.\n\n"
                "  Re-run as:\n"
                "    TBL_SUFFIX=%s python3 33_rebuild_summaries_dedup.py \\\n"
                "        --discovery '%s'" % (d, d[d.rfind("_rl"):], d))
    if a.out is None:
        a.out = "summaries_dedup%s.txt" % SUF
    elif SUF and SUF not in str(a.out):
        raise SystemExit(
            "REFUSING — --out %r does not contain the suffix %r, so the summary\n"
            "  and the two .tsv outputs would describe different read lengths.\n"
            "  Use --out summaries_dedup%s.txt, or omit --out and it is derived."
            % (a.out, SUF, SUF))
    return a
# --------------------------------------------------------------------------


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fdr", type=float, default=0.05)
    ap.add_argument("--min-dpsi", type=float, default=0.10,
                    help="EFFECT THRESHOLD. The project defines a significant event as "
                         "FDR<0.05 AND |dPSI|>=0.10 -- that is what produced the locked "
                         "926 SE / 1,445 all-type counts and the 1,960/2,574 pooled "
                         "concordance. Filtering on FDR alone gives 5,301 MAGNet SE "
                         "events and numbers that cannot be compared with anything "
                         "already written down. Set 0 only to reproduce an FDR-only "
                         "figure deliberately.")
    ap.add_argument("--discovery", default="{ET}_v2filtered")
    ap.add_argument("--out", default=None,
                    help="default: summaries_dedup<TBL_SUFFIX>.txt")
    a = ap.parse_args()

    a = _resolve_suffix(a)
    out = open(a.out, "w")

    def P(s=""):
        print(s)
        out.write(s + "\n")

    # Which cohorts actually lost samples? Only for those is the absence of a _dedup
    # table a problem. Brigham and Fuwai had no exclusions, so their original tables
    # ARE the post-exclusion tables and flagging them would be a false alarm that
    # trains you to ignore the warning.
    AFFECTED = set()
    if os.path.exists("geno/EXCLUDE_samples.tsv"):
        E = pd.read_csv("geno/EXCLUDE_samples.tsv", sep="\t")
        if "drop_cohort" in E.columns:
            AFFECTED = set(E.drop_cohort.astype(str))
    if not AFFECTED:
        AFFECTED = {"Colorado", "Stanford", "LBNL"}
        print("[warn] geno/EXCLUDE_samples.tsv not readable -- assuming the affected "
              "cohorts are Colorado, Stanford, LBNL")

    used, fell_back = [], []

    def dcol(t):
        return next((c for c in ("dPSI_adj", "dPSI") if c in t.columns), None)

    def sel(t):
        """The project's significance definition: FDR AND effect size."""
        d = dcol(t)
        m = t.FDR < a.fdr
        if a.min_dpsi > 0 and d:
            m = m & (t[d].abs() >= a.min_dpsi)
        return m

    def note(path, cohort):
        used.append(path)
        if "NOT FOUND" not in path and "_dedup" not in path and cohort in AFFECTED:
            fell_back.append(path)

    # ---------- A. Colorado internal control ----------------------------------------
    P("=" * 92)
    P("A. COLORADO INTERNAL CONTROL -- same 36 DCM, 13 NF vs 13 ICM  (F3 panel A)")
    P("=" * 92)
    P(f"  Sections B-E select events at FDR<{a.fdr}"
      + (f" AND |dPSI|>={a.min_dpsi} -- the project's definition." if a.min_dpsi > 0
         else " only (--min-dpsi 0) -- NOT the project's definition."))
    P(f"  {'type':<6}{'contrast':<12}{'tested':>9}{'FDR<.05':>10}"
      f"{'+|dPSI|>=.05':>14}{'+|dPSI|>=.10':>14}")
    ic = []
    for et in TYPES:
        for label, base in [("DCMvsNF", f"{et}_Colorado_DCMvsNF"),
                            ("DCMvsICM", f"{et}_Colorado_DCMvsICM")]:
            t, p = load(base, et)
            note(p, "Colorado")
            if t is None:
                P(f"  {et:<6}{label:<12}{'-- ' + p:>47}")
                continue
            d = next((c for c in ("dPSI_adj", "dPSI") if c in t.columns), None)
            sig = t.FDR < a.fdr
            n05 = int((sig & (t[d].abs() >= 0.05)).sum()) if d else -1
            n10 = int((sig & (t[d].abs() >= 0.10)).sum()) if d else -1
            mark = "  <<< Colorado MXE is EXCLUDED" if et == "MXE" else ""
            P(f"  {et:<6}{label:<12}{len(t):>9}{int(sig.sum()):>10}{n05:>14}{n10:>14}{mark}")
            ic.append({"event_type": et, "contrast": label, "tested": len(t),
                       "sig": int(sig.sum()), "sig05": n05, "sig10": n10})
    if ic:
        I = pd.DataFrame(ic)
        I.to_csv("colorado_internal_control_dedup" + SUF + ".tsv", sep="\t", index=False)
        ok = I[I.event_type != "MXE"]
        nf = ok[ok.contrast == "DCMvsNF"]
        icm = ok[ok.contrast == "DCMvsICM"]
        for col, name in [("sig", f"FDR<{a.fdr}"), ("sig10", "|dPSI|>=0.10")]:
            s1, s2 = int(nf[col].sum()), int(icm[col].sum())
            ratio = f"{s1 / s2:.0f}x" if s2 else "no ICM events"
            P(f"\n  EXCLUDING MXE, all four remaining types, {name}:")
            P(f"    DCM vs NF {s1}   vs   DCM vs ICM {s2}   ({ratio})")
        P(f"\n  -> colorado_internal_control_dedup{SUF}.tsv")
        P("  Lead with the RATIO, not an absolute zero: one cohort, one pipeline,")
        P("  identical group sizes, nearly identical numbers tested.")

    # ---------- B. equivalence ratios -------------------------------------------------
    P()
    P("=" * 92)
    P("B. EQUIVALENCE RATIOS  median |dPSI(DCM-ICM)| / |dPSI(DCM-NF)|")
    P("=" * 92)
    P("  computed on events significant in Colorado DCM-vs-NF, matched by exon-level key")
    for et in TYPES:
        A_, pa = load(f"{et}_Colorado_DCMvsNF", et)
        B_, pb = load(f"{et}_Colorado_DCMvsICM", et)
        if A_ is None or B_ is None:
            P(f"  {et:<6}-- missing table")
            continue
        d = next((c for c in ("dPSI_adj", "dPSI") if c in A_.columns), None)
        m = (A_.loc[sel(A_), ["ekey", d]].rename(columns={d: "nf"})
             .merge(B_[["ekey", d]].rename(columns={d: "icm"}), on="ekey"))
        if len(m) < 5:
            P(f"  {et:<6}n={len(m)} -- too few matched events")
            continue
        r = float(np.median(m.icm.abs() / m.nf.abs().replace(0, np.nan)))
        P(f"  {et:<6}n={len(m):<7}ratio {r:.3f}"
          + ("   (MXE excluded from the claim)" if et == "MXE" else ""))

    # ---------- C. DCM-vs-NF against ICM-vs-NF ---------------------------------------
    P()
    P("=" * 92)
    P("C. COLORADO DCM-vs-NF AGAINST ICM-vs-NF  (shared 13 NF hearts)")
    P("=" * 92)
    P("  These two contrasts SHARE their control group, so a positive correlation is")
    P("  expected under the null. Script 05 put that baseline at rho = +0.591 for 13")
    P("  controls. Report the DIRECT DCM-vs-ICM contrast as the primary evidence; this")
    P("  table is supporting only.")
    for et in TYPES:
        A_, _ = load(f"{et}_Colorado_DCMvsNF", et)
        B_, _ = load(f"{et}_Colorado_ICMvsNF", et)
        if A_ is None or B_ is None:
            continue
        d = next((c for c in ("dPSI_adj", "dPSI") if c in A_.columns), None)
        m = (A_.loc[sel(A_), ["ekey", d]].rename(columns={d: "x"})
             .merge(B_[["ekey", d]].rename(columns={d: "y"}), on="ekey"))
        if len(m) < 10:
            continue
        same = int((np.sign(m.x) == np.sign(m.y)).sum())
        r = float(np.corrcoef(m.x, m.y)[0, 1])
        P(f"  {et:<6}n={len(m):<7}same direction {same}/{len(m)} "
          f"({100 * same / len(m):.1f}%)   r={r:.3f}"
          + ("   <<< MXE EXCLUDED" if et == "MXE" else ""))

    # ---------- D. direction concordance, pooled WITHOUT LBNL -------------------------
    P()
    P("=" * 92)
    P("D. DIRECTION CONCORDANCE OF MAGNET-SIGNIFICANT SE EVENTS")
    P("=" * 92)
    disc, pd_ = load(a.discovery.format(ET="SE"), "SE", prefer_dedup=False)
    if disc is None:
        P(f"  discovery table not found: {pd_}")
    else:
        P(f"  discovery: {pd_}  ({len(disc)} exon-level events)")
        dd = next((c for c in ("dPSI_adj", "dPSI") if c in disc.columns), None)
        D = disc.loc[sel(disc), ["ekey", dd]].rename(columns={dd: "x"})
        P(f"  selection: FDR<{a.fdr}" + (f" AND |dPSI|>={a.min_dpsi}"
                                         if a.min_dpsi > 0 else " (FDR only)"))
        P(f"  MAGNet significant SE (exon-level): {len(D)}")
        if a.min_dpsi > 0 and not (700 <= len(D) <= 1200):
            P(f"  [warn] expected ~926 from the locked convention, got {len(D)}."
              f" Check that {pd_} is the table 24_lock_event_convention.py used.")
        P()
        P(f"  {'cohort':<11}{'table':<34}{'n':>7}{'same dir':>11}{'%':>8}{'binom p':>12}")
        tot_n = tot_s = 0
        rows = []
        for coh, pat, _, why in REPL:
            t, p = load(pat.format(ET="SE"), "SE")
            if t is None:
                P(f"  {coh:<11}{'-- not found':<34}")
                continue
            note(p, coh)
            d2 = next((c for c in ("dPSI_adj", "dPSI") if c in t.columns), None)
            m = D.merge(t[["ekey", d2]].rename(columns={d2: "y"}), on="ekey")
            if len(m) < 5:
                P(f"  {coh:<11}{os.path.basename(p):<34}{len(m):>7}  too few matched")
                continue
            same = int((np.sign(m.x) == np.sign(m.y)).sum())
            pv = stats.binomtest(same, len(m), 0.5, alternative="greater").pvalue
            flag = "  EXCLUDED FROM POOL" if why else ""
            P(f"  {coh:<11}{os.path.basename(p):<34}{len(m):>7}{same:>11}"
              f"{100 * same / len(m):>7.1f}%{pv:>12.3g}{flag}")
            rows.append({"cohort": coh, "n": len(m), "same": same,
                         "pct": 100 * same / len(m), "p": pv,
                         "pooled": why is None})
            if why is None:
                tot_n += len(m)
                tot_s += same
        if tot_n:
            pv = stats.binomtest(tot_s, tot_n, 0.5, alternative="greater").pvalue
            P()
            P(f"  POOLED WITHOUT LBNL : {tot_s}/{tot_n} = {100 * tot_s / tot_n:.1f}%"
              f"   binomial p = {pv:.3g}")
            P("  Replaces the previously reported 1,960/2,574 = 76.1%, which included")
            P("  LBNL. 24 of LBNL's 33 samples are MAGNet hearts, so that figure pooled")
            P("  the discovery cohort with itself and was inflated. Use this one.")
            pd.DataFrame(rows).to_csv("direction_concordance_dedup" + SUF + ".tsv", sep="\t",
                                      index=False, float_format="%.5g")
            P(f"  -> direction_concordance_dedup{SUF}.tsv")
        for coh, pat, _, why in REPL:
            if why:
                P(f"\n  [{coh} excluded] {why}")

    # ---------- E. transfer slopes ----------------------------------------------------
    P()
    P("=" * 92)
    P("E. TRANSFER SLOPES (Deming, errors-in-variables) on MAGNet-significant SE events")
    P("=" * 92)
    P("  These are LOWER BOUNDS: events are selected on MAGNet significance, so")
    P("  regression to the mean attenuates every slope. Do not read them as magnitude")
    P("  estimates, and do not compare them across cohorts (script 22).")
    if disc is not None:
        sx = "se_adj" if "se_adj" in disc.columns else None
        D2 = disc.loc[sel(disc), ["ekey", dd] + ([sx] if sx else [])]
        D2 = D2.rename(columns={dd: "x", sx: "sx"} if sx else {dd: "x"})
        P()
        P(f"  {'cohort':<11}{'n':>7}{'Deming slope':>15}")
        for coh, pat, _, why in REPL:
            t, p = load(pat.format(ET="SE"), "SE")
            if t is None:
                continue
            d2 = next((c for c in ("dPSI_adj", "dPSI") if c in t.columns), None)
            s2 = "se_adj" if "se_adj" in t.columns else None
            cols = ["ekey", d2] + ([s2] if s2 else [])
            m = D2.merge(t[cols].rename(columns={d2: "y", s2: "sy"} if s2
                                        else {d2: "y"}), on="ekey")
            if len(m) < 10:
                continue
            if "sx" in m.columns and "sy" in m.columns:
                sl, n = deming(m.x.values, m.y.values, m.sx.values, m.sy.values)
            else:
                sl, n = (float(np.polyfit(m.x, m.y, 1)[0]), len(m))
            note = "   (same hearts as MAGNet -- technical, not replication)" \
                if coh == "LBNL" else ""
            P(f"  {coh:<11}{n:>7}{sl:>15.3f}{note}")

    # ---------- provenance -------------------------------------------------------------
    P()
    P("=" * 92)
    P("PROVENANCE")
    P("=" * 92)
    fb = sorted(set(x for x in fell_back if "NOT FOUND" not in x))
    miss = sorted(set(x for x in used if "NOT FOUND" in x))
    P(f"  cohorts that lost samples: {', '.join(sorted(AFFECTED))}")
    P("  cohorts with no exclusions keep their original tables -- that is correct,")
    P("  not a fallback, and they are not flagged below.")
    if fb:
        P("\n  *** FELL BACK TO PRE-EXCLUSION TABLES for an AFFECTED cohort --")
        P("  *** these numbers are NOT post-dedup:")
        for x in fb:
            P(f"      {x}")
        P("  *** Rerun 32_rerun_dedup_contrasts.sbatch, or treat the affected rows as")
        P("  *** provisional. A mix of pre- and post-exclusion inputs is not quotable.")
    else:
        P("  All cohort tables read from _dedup versions. No pre-exclusion inputs used.")
    if miss:
        P("\n  missing tables (rows were skipped):")
        for x in miss:
            P(f"      {x}")
    P(f"\n[write] {a.out}")
    out.close()


if __name__ == "__main__":
    main()
