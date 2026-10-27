#!/usr/bin/env python3
"""
24_lock_event_convention.py -- settle how an "event" is counted, then freeze every
number the manuscript will quote into one file.

THE AMBIGUITY
-------------
phase2_diffsplice_v2.py builds `key` from the FULL coordinate tuple, e.g. for SE:
    chr : exonStart_0base : exonEnd : upstreamES : upstreamEE : downstreamES : downstreamEE
The same cassette exon therefore appears as SEVERAL "events" whenever different
transcripts give it different flanking exons. MXE carries eight coordinates and has the
most flanking context to vary, which is why its count moves ~31% between conventions.

  FULL       one row per (alternative exon + flanking context). What rMATS emits and what
             the pipeline has used throughout.
  EXON       one row per alternative splicing CHOICE, flanking context discarded:
               SE    chr + cassette exon
               MXE   chr + both alternative exons
               A3SS  chr + long and short forms   (long-vs-short IS the event)
               A5SS  chr + long and short forms
               RI    chr + retained intron

WHY EXON-LEVEL SHOULD BE PRIMARY
--------------------------------
Not tidiness. A full-coordinate count depends on how many TRANSCRIPTS happen to be
annotated for a gene, so it is not comparable across annotation versions -- and Fuwai's
STAR index was GENCODE v46 while everything else is v49. Exon-level counts are invariant
to that. It also matches what a reader means by "N splicing events changed".

Report full-coordinate counts in the supplement and state in the methods that FDR was
computed over full-coordinate rows -- the conservative direction, since more tests means
a stricter BH threshold -- and collapsed afterwards for reporting.

WHAT THIS SCRIPT DOES
---------------------
1. Counts every table under BOTH conventions, so the choice is documented rather than
   assumed.
2. Re-checks the replication conclusions under both. If direction concordance survives
   the switch, the convention is a reporting decision and not a scientific one -- which
   is exactly what you want to be able to say to a reviewer.
3. Writes `master_numbers.tsv` under the chosen convention. NO NUMBER SHOULD ENTER THE
   MANUSCRIPT THAT IS NOT IN THAT FILE. The v1/v2 drift already put stale figures into
   two project docs; a single frozen source is the fix.

Collapsing keeps the lowest-FDR row in each group, so significance is never invented by
the merge -- an exon called significant after collapsing was significant in at least one
of its transcript contexts.

Usage (login node, a few minutes):
  python 24_lock_event_convention.py
  python 24_lock_event_convention.py --convention full     # if you decide otherwise
"""
import argparse
import glob
import os
import re

import numpy as np
import pandas as pd
from scipy import stats

# how many ":"-separated fields of `key` define the splicing choice (chr counts as one)
NKEEP = {"SE": 3, "MXE": 5, "A3SS": 5, "A5SS": 5, "RI": 3}
TYPES = ["SE", "MXE", "A3SS", "A5SS", "RI"]

# label -> filename template. {ET} is substituted.
TABLES = [
    ("MAGNet",   "DCM vs NF",   "{ET}_v2filtered.all_events.tsv"),
    ("MAGNet",   "DCM vs NF (composition-adjusted)", "{ET}_v2compadj.all_events.tsv"),
    ("MAGNet",   "HCM vs NF",   "{ET}_HCM_v2.all_events.tsv"),
    ("Colorado", "DCM vs NF",   "{ET}_Colorado_DCMvsNF.all_events.tsv"),
    ("Colorado", "ICM vs NF",   "{ET}_Colorado_ICMvsNF.all_events.tsv"),
    ("Colorado", "DCM vs ICM",  "{ET}_Colorado_DCMvsICM.all_events.tsv"),
    ("LBNL",     "DCM vs NF",   "{ET}_LBNL_DCMvsNF.all_events.tsv"),
    ("Brigham",  "HCM vs ctrl", "Brigham_{ET}_HCM.all_events.tsv"),
    ("Fuwai",    "HCM vs ctrl", "Fuwai_{ET}_HCM.all_events.tsv"),
    ("Stanford", "HCM vs ctrl", "Stanford_{ET}_HCM.all_events.tsv"),
]

# replication checks: (label, discovery template, replication template)
REPS = [
    ("Colorado", "{ET}_v2filtered.all_events.tsv", "{ET}_Colorado_DCMvsNF.all_events.tsv"),
    ("LBNL",     "{ET}_v2filtered.all_events.tsv", "{ET}_LBNL_DCMvsNF.all_events.tsv"),
    ("Brigham",  "{ET}_v2filtered.all_events.tsv", "Brigham_{ET}_HCM.all_events.tsv"),
    ("Fuwai",    "{ET}_v2filtered.all_events.tsv", "Fuwai_{ET}_HCM.all_events.tsv"),
    ("Stanford", "{ET}_v2filtered.all_events.tsv", "Stanford_{ET}_HCM.all_events.tsv"),
]


def collapse(keys, et):
    """Strip flanking coordinates, keeping only the fields that define the choice."""
    n = NKEEP[et]
    out = []
    for k in keys.astype(str):
        p = k.split(":")
        out.append(":".join(p[:n]) if len(p) >= n else k)
    return np.array(out)


# Tables for cohorts that lost genotype-duplicate samples exist in two versions:
# the original and a _dedup one written by 32_rerun_dedup_contrasts.sbatch. The
# hardcoded names below are the originals, so without this resolver every number in
# master_numbers.tsv would silently be rebuilt from PRE-EXCLUSION data. Prefer the
# _dedup table whenever one exists, and say which was used.
_USED = []


def resolve(path):
    d = path.replace(".all_events.tsv", "_dedup.all_events.tsv")
    if d != path and os.path.exists(d) and os.path.getsize(d):
        _USED.append((path, d))
        return d
    return path


def load(path, et):
    path = resolve(path)
    if not (os.path.exists(path) and os.path.getsize(path)):
        return None
    t = pd.read_csv(path, sep="\t")
    if not {"key", "gene", "dPSI_adj", "FDR"} <= set(t.columns):
        return None
    t = t.drop_duplicates("key").copy()
    t["ekey"] = collapse(t["key"], et)
    return t


def to_exon(t):
    """One row per exon-level event: the lowest-FDR transcript context."""
    return t.loc[t.groupby("ekey")["FDR"].idxmin()].copy()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--convention", choices=["exon", "full"], default="exon")
    ap.add_argument("--fdr", type=float, default=0.05)
    ap.add_argument("--min-dpsi", type=float, default=0.10)
    ap.add_argument("--types", default=",".join(TYPES))
    ap.add_argument("--out", default="master_numbers")
    a = ap.parse_args()
    types = [t.strip() for t in a.types.split(",") if t.strip()]

    lines = []
    def say(s=""):
        lines.append(s); print(s)

    say("=" * 100)
    say(f"EVENT-COUNTING CONVENTION   (primary = {a.convention.upper()}-level, "
        f"FDR<{a.fdr}, |dPSI|>={a.min_dpsi})")
    say()

    # ---------- 1. counts under both conventions --------------------------
    rows = []
    for cohort, contrast, tmpl in TABLES:
        for et in types:
            f = tmpl.format(ET=et)
            t = load(f, et)
            if t is None:
                continue
            e = to_exon(t)
            sig = lambda d: int(((d.FDR < a.fdr) & (d.dPSI_adj.abs() >= a.min_dpsi)).sum())
            rows.append({
                "cohort": cohort, "contrast": contrast, "event_type": et, "file": f,
                "tested_full": len(t), "tested_exon": len(e),
                "collapse_pct": 100 * (1 - len(e) / len(t)) if len(t) else np.nan,
                "sig_full": sig(t), "sig_exon": sig(e),
                "sig_collapse_pct": (100 * (1 - sig(e) / sig(t)) if sig(t) else np.nan),
                "genes_sig": int(t[(t.FDR < a.fdr) &
                                   (t.dPSI_adj.abs() >= a.min_dpsi)].gene.nunique()),
                "n_g1": int(np.nanmedian(t["n_g1"])) if "n_g1" in t else np.nan,
                "n_g2": int(np.nanmedian(t["n_g2"])) if "n_g2" in t else np.nan,
            })
    if not rows:
        raise SystemExit("[fatal] no tables found -- run this in the HEART directory")
    C = pd.DataFrame(rows)

    say("1. HOW MUCH THE CONVENTION MOVES THE COUNTS")
    say()
    say(f"  {'cohort':<10}{'contrast':<34}{'type':<6}{'tested full':>12}{'exon':>9}"
        f"{'collapse':>10}{'sig full':>10}{'sig exon':>10}{'collapse':>10}")
    for _, r in C.iterrows():
        say(f"  {r.cohort:<10}{r.contrast:<34}{r.event_type:<6}{r.tested_full:>12}"
            f"{r.tested_exon:>9}{r.collapse_pct:>9.1f}%{r.sig_full:>10}{r.sig_exon:>10}"
            + (f"{r.sig_collapse_pct:>9.1f}%" if np.isfinite(r.sig_collapse_pct)
               else f"{'--':>10}"))
    say()
    say("  collapse % by event type, tested events (this is the number to quote when")
    say("  justifying the convention in the methods):")
    for et, g in C.groupby("event_type"):
        say(f"    {et:<6} {g.collapse_pct.mean():>5.1f}%   "
            f"(range {g.collapse_pct.min():.1f}-{g.collapse_pct.max():.1f}%)")
    say()

    # ---------- 2. do the conclusions survive the switch? ------------------
    say("2. DOES THE CHOICE CHANGE ANY CONCLUSION?")
    say("   Direction concordance of MAGNet-significant events, under both conventions.")
    say("   If these agree, the convention is a REPORTING decision, not a scientific one.")
    say()
    say(f"  {'cohort':<10}{'type':<6}{'n full':>8}{'same full':>12}{'n exon':>9}"
        f"{'same exon':>12}{'shift':>8}")
    rep_rows = []
    for lab, dt, rt in REPS:
        for et in types:
            D, R = load(dt.format(ET=et), et), load(rt.format(ET=et), et)
            if D is None or R is None:
                continue
            out = {}
            for conv, (dd, rr, on) in {
                    "full": (D, R, "key"),
                    "exon": (to_exon(D), to_exon(R), "ekey")}.items():
                s = dd[(dd.FDR < a.fdr) & (dd.dPSI_adj.abs() >= a.min_dpsi)]
                m = s[[on, "dPSI_adj"]].merge(rr[[on, "dPSI_adj"]], on=on,
                                              suffixes=("_d", "_r")).dropna()
                if len(m) < 5:
                    out[conv] = (len(m), np.nan)
                    continue
                same = float(np.mean(np.sign(m.dPSI_adj_d) == np.sign(m.dPSI_adj_r)))
                out[conv] = (len(m), 100 * same)
            nf, pf = out.get("full", (0, np.nan))
            ne, pe = out.get("exon", (0, np.nan))
            rep_rows.append({"cohort": lab, "event_type": et, "n_full": nf,
                             "same_full": pf, "n_exon": ne, "same_exon": pe,
                             "delta_pp": pe - pf})
            say(f"  {lab:<10}{et:<6}{nf:>8}"
                + (f"{pf:>11.1f}%" if np.isfinite(pf) else f"{'--':>12}")
                + f"{ne:>9}"
                + (f"{pe:>11.1f}%" if np.isfinite(pe) else f"{'--':>12}")
                + (f"{pe - pf:>+7.1f}" if np.isfinite(pf) and np.isfinite(pe)
                   else f"{'--':>8}"))
    RP = pd.DataFrame(rep_rows)
    say()
    worst = RP["delta_pp"].abs().max() if len(RP) else np.nan
    if len(RP) and worst < 5:
        say(f"  >>> Largest shift {worst:.1f} percentage points. The")
        say("  convention does NOT change any replication conclusion. Say so in the")
        say("  methods -- it turns an arbitrary-looking choice into a documented one.")
    elif len(RP):
        say(f"  *** Largest shift {worst:.1f} points -- inspect the rows")
        say("  above before committing. A conclusion that depends on the convention has")
        say("  to be reported under both.")
    say()

    # ---------- 3. freeze ---------------------------------------------------
    C["convention"] = a.convention
    C["tested"] = C["tested_exon"] if a.convention == "exon" else C["tested_full"]
    C["significant"] = C["sig_exon"] if a.convention == "exon" else C["sig_full"]
    cols = ["cohort", "contrast", "event_type", "convention", "tested", "significant",
            "genes_sig", "n_g1", "n_g2", "tested_full", "tested_exon", "sig_full",
            "sig_exon", "collapse_pct", "file"]
    C[cols].to_csv(f"{a.out}.tsv", sep="\t", index=False, float_format="%.4g")
    if len(RP):
        RP.to_csv(f"{a.out}_convention_sensitivity.tsv", sep="\t", index=False,
                  float_format="%.4g")

    say("3. FROZEN")
    say(f"  {a.out}.tsv   <-- every figure and table script reads from this file.")
    say("  NO NUMBER ENTERS THE MANUSCRIPT THAT IS NOT IN IT.")
    say()
    tot = C[C.contrast == "DCM vs NF"]
    mg = tot[tot.cohort == "MAGNet"]
    if len(mg):
        say(f"  MAGNet DCM vs NF, all event types, {a.convention}-level:")
        say(f"    tested      {int(mg.tested.sum()):>8}")
        say(f"    significant {int(mg.significant.sum()):>8}   "
            f"(FDR<{a.fdr} and |dPSI|>={a.min_dpsi})")
        say(f"    full-coordinate equivalent: {int(mg.sig_full.sum())}")
    say()
    say("METHODS SENTENCE")
    say("  'rMATS reports one row per alternative exon per flanking-transcript context.")
    say("   Because the number of contexts depends on annotation density, and one cohort")
    say("   was aligned against GENCODE v46 rather than v49, events were collapsed to the")
    say("   coordinates defining the splicing choice (the cassette exon for SE, both")
    say("   alternative exons for MXE, the long and short forms for A3SS/A5SS, the")
    az = "   retained intron for RI), retaining the lowest-FDR context. FDR was computed"
    say(az)
    say("   over uncollapsed rows, which is conservative. Counts under both conventions")
    say("   are given in Supplementary Table SX; replication concordance differs by less")
    say("   than N percentage points between them.'")
    say("=" * 100)

    with open(f"{a.out}.report.txt", "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"\n[write] {a.out}.tsv  {a.out}_convention_sensitivity.tsv  {a.out}.report.txt")


def _report_dedup():
    if _USED:
        print()
        print("=" * 78)
        print(f"POST-EXCLUSION TABLES USED ({len(_USED)} substitutions)")
        print("=" * 78)
        for a, b in sorted(set(_USED)):
            print(f"  {os.path.basename(a):<44}-> {os.path.basename(b)}")
        print("  Every other table had no _dedup version, which is correct for cohorts")
        print("  that lost no samples (MAGNet, Brigham, Fuwai).")
    else:
        print()
        print("*** NO _dedup TABLES FOUND. If 32_rerun_dedup_contrasts.sbatch has run,")
        print("*** these numbers are PRE-EXCLUSION and must not go in the manuscript.")


import atexit
atexit.register(_report_dedup)

if __name__ == "__main__":
    main()
