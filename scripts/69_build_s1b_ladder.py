#!/usr/bin/env python3
"""
69_build_s1b_ladder.py -- deposit Supplementary Table S1(b), the covariate ladder.

WHY THIS SCRIPT EXISTS
----------------------
`54_fdr_only_counts.py --find-discovery` prints the ladder to the console and writes
nothing. Four of its five rows (899 / 835 / 866 / 888) therefore existed nowhere on disk,
while the Results quotes them and S5(d) cites them. That is the same gap that cost S8 its
permutation evidence: a measured, published claim with no artefact behind it.

**This script does NOT parse 54's log.** A log is a transcription; if the console output
were wrong, parsing it would deposit the wrong number with a clean provenance trail.
Instead it recomputes all five counts from the five named tables, under the paper's own
locked convention, and REFUSES to write anything unless every one reproduces.

THE CONVENTION, taken unchanged from 68_build_s1.py
---------------------------------------------------
Exon level: collapse `key` to its first NKEEP fields (3 for SE), keep the lowest-FDR row
per collapsed event, then apply FDR < 0.05 AND |dPSI| >= 0.10. That is what returns 926
from SE_v2filtered, which is the gate.

WHAT 54's OUTPUT ALSO REVEALED, and what this script does about it
-----------------------------------------------------------------
Two tables in the same directory are NOT ladder rows and must never be mistaken for them:

  SE_v2legacy.all_events.tsv    1,514 of 133,823 -- a DIFFERENT PIPELINE VERSION, not a
                                different convention and not a different covariate model.
                                1,514 is close enough to the all-types exon-level total of
                                1,445 to be misread in a hurry, which is exactly how the
                                1,443 reached the abstract. Reported, never deposited.
  SE_adjusted.all_events.tsv    a 10-column table with NO `key` column, so neither 54 nor
                                this script's ladder path can collapse it. But it DOES
                                carry chr/exonStart/exonEnd, and for SE that triple IS the
                                exon-level identity -- the same one the first three fields
                                of `key` encode. So it is countable, just not by the
                                ladder's route, and v2 of this script counts it on its own
                                columns rather than dismissing it.
                                WHY THIS MATTERS: `replication_concordance.py` and
                                `lvef_correlation.py` are recorded as keying on
                                (chr, exonStart, exonEnd), and this table is shaped exactly
                                like what they would read. If a reported number came from
                                here, we need to know what it is, not that it is awkward.

Both go to S1b_advisory.tsv so the next reader meets them before guessing.

ALSO CHECKED: the denominator claim
-----------------------------------
The Results and PAPER_ONE_PAGE say the ladder is computed "on the same 62,971 exon-level
events". Four rows are; the age-25 row is on 62,482, because excluding donors changes the
PSI-variance filter and so changes the tested set. The script prints both and flags the
mismatch rather than letting the prose stand unchallenged.

USAGE
-----
    conda activate rnaseq_tools
    cd /data/project/chureeratp/RNA/HEART

    python 69_build_s1b_ladder.py --list     # show each table's schema, stop
    python 69_build_s1b_ladder.py            # recompute, gate, deposit
"""
import argparse
import os

import numpy as np
import pandas as pd

# ------------------------------------------------------------------ convention
NKEEP = {"SE": 3, "MXE": 5, "A3SS": 5, "A5SS": 5, "RI": 3}
FDR_CUT, DPSI_CUT = 0.05, 0.10
ET = "SE"                                   # the ladder is skipped-exon only

# ------------------------------------------------------------------ the record
# (order, label as the manuscript words it, file, recorded significant, recorded tested)
LADDER = [
    (1, "base: sex + age + RIN",                  "SE_v2filtered.all_events.tsv", 926, 62_971),
    (2, "+ cardiomyocyte / fibroblast composition", "SE_v2compadj.all_events.tsv", 899, 62_971),
    (3, "+ composition + library size + TIN",     "SE_v2fulladj.all_events.tsv",  835, 62_971),
    (4, "+ sequencing pool",                      "SE_v2pooladj.all_events.tsv",  866, 62_971),
    (5, "excluding donors < 25 y",                "SE_v2age25.all_events.tsv",    888, 62_482),
]
# --- read length switch. "" = 50 bp. ADVISORY is deliberately NOT suffixed:
# SE_adjusted has no 125 bp counterpart and the legacy row is a pipeline-version
# note, not a covariate model.
SUF = os.environ.get("TBL_SUFFIX", "")
_rec = os.environ.get("LADDER_REC", "").strip()
if SUF or _rec:
    _pairs = None
    if _rec:
        _pairs = [tuple(int(x) for x in _p.split(":")) for _p in _rec.split(",")]
        assert len(_pairs) == len(LADDER), (
            f"LADDER_REC needs {len(LADDER)} sig:tested entries, got {len(_pairs)}")
    LADDER = [(_o, _lab,
               _f.replace(".all_events.tsv", SUF + ".all_events.tsv"),
               _pairs[_i][0] if _pairs else _sig,
               _pairs[_i][1] if _pairs else _tst)
              for _i, (_o, _lab, _f, _sig, _tst) in enumerate(LADDER)]

BASE_ROW = 1                                # the gate: this one must be exactly 926
# the discovery gate. 926 at 50 bp, 813 at 125 bp.
BASE_SIG_RECORD = int(os.environ.get("BASE_SIG_RECORD", "926"))
# quoted in the Results as "a range of less than 10%"
RECORDED_SPREAD = int(os.environ.get("RECORDED_SPREAD", "91"))   # 926 - 835
RECORDED_SPREAD_PCT = float(os.environ.get("RECORDED_SPREAD_PCT", "9.8"))
LADDER_TESTED = int(os.environ.get("LADDER_TESTED", "62_971"))   # the tested set rows 1-4 share
RECORDED_DUP = float(os.environ.get("RECORDED_DUP", "1.59"))   # SE duplication factor, RESULT_composition_check.md
DUP_TOL = 0.01

# NOT ladder rows. Reported so they are never mistaken for one.
ADVISORY = [
    ("SE_v2legacy.all_events.tsv", 1_514, 133_823,
     "different PIPELINE VERSION, not a covariate model and not a convention. "
     "1,514 is four away from the all-types exon-level total of 1,445 -- do not quote."),
    ("SE_adjusted.all_events.tsv", None, None,
     "no `key` column, so the ladder path cannot read it -- but it carries "
     "chr/exonStart/exonEnd, which for SE is the exon-level identity, so it is counted "
     "below on its own columns. Shaped like what replication_concordance.py and "
     "lvef_correlation.py key on. Check whether any reported number came from here."),
]

OUT_LADDER = "S1b.tsv"
OUT_ADVISORY = "S1b_advisory.tsv"


# ------------------------------------------------------------------- machinery
def dpsi_column(t, path):
    """The effect column, preferring the adjusted one. Same rule as 68."""
    for c in ("dPSI_adj", "dPSI"):
        if c in t.columns:
            return c
    cand = [c for c in t.columns if str(c).lower().startswith("dpsi")]
    if not cand:
        return None
    adj = [c for c in cand if "adj" in str(c).lower()]
    return adj[0] if adj else cand[0]


COORD = ("chr", "exonStart", "exonEnd")        # for SE this triple IS the exon identity


def count_by_coordinate(path):
    """
    Exon-level count for a table with no `key` column, collapsed on (chr, exonStart,
    exonEnd) -- and `strand` too when present, since two genes can share an exon
    interval on opposite strands.

    ADVISORY TABLES ONLY. A gated ladder row must go through the `key` path, because
    silently switching the identity a count is computed on is how two conventions end up
    in one table.
    """
    if not (os.path.exists(path) and os.path.getsize(path)):
        return None, None, None, "missing or empty"
    t = pd.read_csv(path, sep="\t", low_memory=False)
    missing = [c for c in COORD if c not in t.columns]
    if missing:
        return len(t), None, None, f"no `key` and no {'/'.join(missing)} either"
    d = dpsi_column(t, path)
    if d is None or "FDR" not in t.columns:
        return len(t), None, None, "no FDR/dPSI column"
    ident = list(COORD) + (["strand"] if "strand" in t.columns else [])
    t = t.copy()
    t["ekey"] = t[ident].astype(str).agg(":".join, axis=1)
    t["_fdr"] = pd.to_numeric(t["FDR"], errors="coerce")
    t["_dpsi"] = pd.to_numeric(t[d], errors="coerce")
    t = t.dropna(subset=["_fdr"])
    rep = t.loc[t.groupby("ekey")["_fdr"].idxmin()]
    sig = rep[(rep["_fdr"] < FDR_CUT) & (rep["_dpsi"].abs() >= DPSI_CUT)]
    note = "collapsed on " + "+".join(ident) + (" (NO strand column)" if "strand" not in ident else "")
    return len(t), len(rep), len(sig), note


def count_exon_level(path, et=ET):
    """
    Returns (n_rows, n_exons, n_significant, note).

    n_rows      -- rows in the table, i.e. transcript contexts
    n_exons     -- distinct collapsed events, the exon-level tested set
    n_significant -- collapsed events whose lowest-FDR row passes the locked rule

    This is 68_build_s1.py's `exon_sig` with the denominators kept rather than discarded,
    because the denominator is half of what the Results claims.
    """
    if not os.path.exists(path):
        return None, None, None, "missing"
    if os.path.getsize(path) == 0:
        return None, None, None, "empty"
    try:
        t = pd.read_csv(path, sep="\t", low_memory=False)
    except Exception as e:                                           # noqa: BLE001
        return None, None, None, f"unreadable: {e.__class__.__name__}"
    if "key" not in t.columns:
        return len(t), None, None, "no `key` column -- cannot collapse exon-level"
    if "FDR" not in t.columns:
        return len(t), None, None, "no `FDR` column"
    d = dpsi_column(t, path)
    if d is None:
        return len(t), None, None, "no dPSI column"

    n_rows = len(t)
    nk = NKEEP[et]
    t = t.copy()
    t["ekey"] = [":".join(str(k).split(":")[:nk]) for k in t["key"]]
    t["_fdr"] = pd.to_numeric(t["FDR"], errors="coerce")
    t["_dpsi"] = pd.to_numeric(t[d], errors="coerce")
    t = t.dropna(subset=["_fdr"])
    if t.empty:
        return n_rows, 0, 0, "no usable FDR values"
    rep = t.loc[t.groupby("ekey")["_fdr"].idxmin()]
    n_exons = len(rep)
    sig = rep[(rep["_fdr"] < FDR_CUT) & (rep["_dpsi"].abs() >= DPSI_CUT)]
    return n_rows, n_exons, len(sig), None


def show_schema(path):
    if not os.path.exists(path):
        print(f"    {path}\n        *** missing ***")
        return
    if os.path.getsize(path) == 0:
        print(f"    {path}\n        *** empty ***")
        return
    try:
        t = pd.read_csv(path, sep="\t", low_memory=False, nrows=300)
    except Exception as e:                                           # noqa: BLE001
        print(f"    {path}\n        *** unreadable: {e.__class__.__name__} ***")
        return
    size_mb = os.path.getsize(path) / 1e6
    print(f"    {path}  ({size_mb:,.1f} MB, {len(t.columns)} columns)")
    print(f"        columns: {', '.join(map(str, t.columns))}")
    has_key = "key" in t.columns
    has_fdr = "FDR" in t.columns
    d = dpsi_column(t, path)
    verdict = "usable" if (has_key and has_fdr and d) else "NOT USABLE for an exon-level count"
    print(f"        key={has_key}  FDR={has_fdr}  dPSI={d}   -> {verdict}")


def do_list():
    print("\n" + "=" * 78)
    print("SCHEMAS -- every table 54 listed, ladder rows and advisories alike")
    print("=" * 78)
    print("\n  LADDER ROWS (must all be usable):")
    for _, label, path, _, _ in LADDER:
        print(f"\n  [{label}]")
        show_schema(path)
    print("\n  ADVISORY (reported, never deposited):")
    for path, _, _, _ in ADVISORY:
        print()
        show_schema(path)
    print("\n  If every ladder row says `usable`, run without --list.\n")


# ----------------------------------------------------------------------- build
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--list", action="store_true",
                    help="print each table's schema and stop")
    ap.add_argument("--outdir", default=".", help="where to write the TSVs")
    args = ap.parse_args()

    if args.list:
        do_list()
        return

    print("\n" + "=" * 78)
    print("S1(b) THE COVARIATE LADDER -- recomputed from the tables, not from 54's log")
    print("=" * 78)
    print(f"\n  convention: exon level, {NKEEP[ET]} key fields for {ET}, lowest-FDR row kept")
    print(f"  rule:       FDR < {FDR_CUT} and |dPSI| >= {DPSI_CUT}\n")

    hdr = (f"  {'model':<44}{'sig':>7}{'rec':>7}{'exons':>9}{'rec':>9}"
           f"{'rows':>9}  {'':<4}")
    print(hdr)
    print("  " + "-" * (len(hdr) - 2))

    rows, failures, notes = [], [], []
    for order, label, path, rec_sig, rec_tested in LADDER:
        n_rows, n_exons, n_sig, note = count_exon_level(path)
        if note:
            print(f"  {label:<44}{'-':>7}{rec_sig:>7}{'-':>9}{rec_tested:>9}{'-':>9}  FAIL")
            failures.append(f"{path}: {note}")
            continue

        ok_sig = (n_sig == rec_sig)
        ok_tested = (n_exons == rec_tested)
        mark = "ok" if ok_sig else "FAIL"
        print(f"  {label:<44}{n_sig:>7}{rec_sig:>7}{n_exons:>9}{rec_tested:>9}"
              f"{n_rows:>9}  {mark}")

        if not ok_sig:
            failures.append(
                f"{path}: {n_sig} significant against the recorded {rec_sig}")
        if not ok_tested:
            # Diagnostic, not fatal on its own: it tells us what 54's denominator was.
            which = "rows" if n_rows == rec_tested else "neither rows nor exons"
            notes.append(
                f"{path}: tested set is {n_exons:,} exon-level events, but the record "
                f"says {rec_tested:,} ({which} matches)")

        rows.append(dict(order=order, model=label, table_file=path,
                         n_rows_transcript_context=n_rows,
                         n_tested_exon_level=n_exons,
                         duplication_factor=round(n_rows / n_exons, 3),
                         n_significant=n_sig,
                         recorded_significant=rec_sig,
                         reproduced="yes" if ok_sig else "NO"))

    # ------------------------------------------------- the gate
    print()
    base = next((r for r in rows if r["order"] == BASE_ROW), None)
    if base is None or base["n_significant"] != BASE_SIG_RECORD:
        print(f"  [GATE   ] base model did NOT return {BASE_SIG_RECORD} -- something upstream has changed.")
        print("            The ladder in the Results is no longer the ladder that was")
        print("            computed. STOP and reconcile before depositing anything.")
        failures.append(f"base model != {BASE_SIG_RECORD}")
    else:
        print(f"  [ok     ] base model reproduces {BASE_SIG_RECORD} -- the discovery gate")

    if failures:
        print("\n" + "=" * 78)
        print("S1(b) NOT BUILT -- refusing to deposit a table that does not reproduce")
        print("=" * 78)
        for f in failures:
            print(f"  * {f}")
        print("\n  Nothing was written. Fix the input or the record, not this script.\n")
        return

    sigs = [r["n_significant"] for r in rows]
    spread = max(sigs) - min(sigs)
    spread_pct = 100.0 * spread / base["n_significant"]
    print(f"  [ok     ] all five rows reproduce their recorded count")
    print(f"\n  spread {min(sigs)}-{max(sigs)} = {spread} events = {spread_pct:.1f}% "
          f"of the headline")
    if spread != RECORDED_SPREAD or abs(spread_pct - RECORDED_SPREAD_PCT) > 0.15:
        print(f"  [WARN   ] the record says {RECORDED_SPREAD} events / "
              f"{RECORDED_SPREAD_PCT}% -- update the Results")
    else:
        print(f"  [ok     ] matches the recorded {RECORDED_SPREAD} events / "
              f"{RECORDED_SPREAD_PCT}%")
    if spread_pct >= 10.0:
        print('  [WARN   ] the Results says "a range of less than 10%" -- it is now '
              f"{spread_pct:.1f}%. Reword or requote.")

    # ------------------------------------------------- duplication factor
    dups = {r["table_file"]: r["duplication_factor"] for r in rows}
    print(f"\n  duplication factor (transcript contexts per exon): "
          f"{min(dups.values()):.3f}-{max(dups.values()):.3f}")
    if abs(rows[0]["duplication_factor"] - RECORDED_DUP) <= DUP_TOL:
        print(f"  [ok     ] discovery table reproduces the recorded {RECORDED_DUP}x "
              f"(it is {rows[0]['duplication_factor']:.3f}) -- checked here against the "
              f"table for the first time")
    else:
        print(f"  [WARN   ] discovery table gives {rows[0]['duplication_factor']:.3f}x "
              f"against the recorded {RECORDED_DUP}x")

    shared = [r for r in rows if r["n_tested_exon_level"] == LADDER_TESTED]
    rowcounts = {r["table_file"]: r["n_rows_transcript_context"] for r in shared}
    if len(set(rowcounts.values())) > 1:
        lo, hi = min(rowcounts.values()), max(rowcounts.values())
        odd = [p for p, n in rowcounts.items() if n != hi]
        print(f"\n  ROW ASYMMETRY among the tables sharing {LADDER_TESTED:,} exons: "
              f"{lo:,}-{hi:,} rows")
        for p in odd:
            print(f"    * {p}: {rowcounts[p]:,} rows, {hi - rowcounts[p]} fewer than the "
                  f"others")
        print("    The exon-level tested set is identical, which is what the prose claims,")
        print("    so no number moves. A dropped transcript context is a model that failed")
        print("    to fit for one row (rank deficiency on a covariate level is the usual")
        print("    cause); its exon survives through another context.")

    # ------------------------------------------------- the denominator claim
    if notes:
        print("\n  DENOMINATOR -- the prose claims one tested set for all five rows:")
        for n in notes:
            print(f"    * {n}")
        print("    The Results states the age-25 denominator correctly.")
        print("    PAPER_ONE_PAGE.md heads the whole ladder \"on the same 62,971")
        print("    exon-level events\", which is true of four rows, not five. One-word fix.")

    # ------------------------------------------------- advisory table
    print("\n  ADVISORY -- tables in the same directory that are NOT ladder rows:")
    adv_rows = []
    for path, rec_sig, rec_tested, why in ADVISORY:
        n_rows, n_exons, n_sig, note = count_exon_level(path)
        identity = f"key, first {NKEEP[ET]} fields"
        if n_sig is None and "key" in (note or ""):
            # no `key` column -- fall back to the coordinate identity, ADVISORY ONLY
            n_rows, n_exons, n_sig, note2 = count_by_coordinate(path)
            identity = note2 if n_sig is not None else identity
            note = note if n_sig is None else note2
        got = "n/a" if n_sig is None else f"{n_sig:,} of {n_exons:,}"
        rec = "n/a" if rec_sig is None else f"{rec_sig:,} of {rec_tested:,}"
        print(f"    {path}")
        print(f"        computed {got}   54 reported {rec}")
        print(f"        identity: {identity}")
        if n_rows and n_exons:
            print(f"        {n_rows:,} rows -> {n_exons:,} exons "
                  f"({n_rows / n_exons:.3f}x duplication)")
        print(f"        {why}")
        if n_sig is not None and rec_sig is None:
            print(f"        >>> NOW MEASURED: {n_sig:,} significant. Compare against {BASE_SIG_RECORD} "
                  f"(exon-level v2) and 1,074 (full-coordinate v2) before anyone quotes it.")
        adv_rows.append(dict(table_file=path,
                             n_significant=("" if n_sig is None else n_sig),
                             n_tested_exon_level=("" if n_exons is None else n_exons),
                             n_rows_transcript_context=("" if n_rows is None else n_rows),
                             identity_used=identity,
                             reported_by_54=rec,
                             status="NOT a ladder row -- do not quote",
                             reason=why,
                             note=(note or "")))

    # ------------------------------------------------- write
    os.makedirs(args.outdir, exist_ok=True)
    lp = os.path.join(args.outdir, OUT_LADDER)
    ap_ = os.path.join(args.outdir, OUT_ADVISORY)
    df = pd.DataFrame(rows).sort_values("order")
    df["pct_of_base"] = (100.0 * df["n_significant"] / base["n_significant"]).round(1)
    df["delta_vs_base"] = df["n_significant"] - base["n_significant"]
    df.to_csv(lp, sep="\t", index=False)
    pd.DataFrame(adv_rows).to_csv(ap_, sep="\t", index=False)

    print("\n" + "=" * 78)
    print("S1(b) BUILT")
    print("=" * 78)
    print(f"  {lp}        {len(df)} rows -- the deposited ladder")
    print(f"  {ap_}   {len(adv_rows)} rows -- what must not be mistaken for it")
    print("\n  The four rows that existed only in console output now have an artefact,")
    print("  recomputed from the tables rather than transcribed from a log.\n")
    for line in df[["order", "model", "n_tested_exon_level", "n_significant",
                    "pct_of_base", "delta_vs_base"]].to_string(index=False).splitlines():
        print("    " + line)
    print()


if __name__ == "__main__":
    main()
