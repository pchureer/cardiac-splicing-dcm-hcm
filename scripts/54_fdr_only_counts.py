#!/usr/bin/env python3
"""
54_fdr_only_counts.py -- the DCM-vs-HCM count at FDR alone, which the paper now needs.

WHY
---
Britto-Borges et al. (J Transl Med 2022;20:441) ran DCM versus HCM on essentially this
discovery cohort -- 165 DCM / 27 HCM / 162 non-failing against our 166 / 28 / 166 -- and
reported 201 differentially expressed genes and **149 differently used transcripts**, at
adjusted P <= 0.05 with NO effect-size filter.

This paper reports 0 of 45,380 at FDR<0.05 AND |dPSI|>=0.10. Those two numbers are not
comparable, and a Discussion paragraph that sets them side by side without saying so
would be misleading. The comparable number is this paper's count at FDR alone, exon
level, which has never been computed: `RESULT_hcm_specificity.md` carries a
full-coordinate "3" and the older docs a full-coordinate "6 of 147,355", both on a
different convention and a different significance rule.

WHAT IT DOES
------------
Counts every contrast it can find at a ladder of thresholds -- FDR only, FDR plus
|dPSI|>=0.05, FDR plus >=0.10 -- at BOTH the exon-level convention and the full
coordinate, so the old numbers can be reconciled rather than contradicted. It recomputes
nothing: these are filters over tables that already exist.

The positive controls are printed alongside deliberately. A DCM-vs-HCM count at FDR
alone means nothing without the same cohort's DCM-vs-NF and HCM-vs-NF counts at the same
threshold, because the question is never "is it zero" but "is it small relative to what
this pipeline finds when there IS a difference".

Usage:
  export PATH=/data/apps/anaconda3/envs/mamba/envs/rnaseq_tools/bin:$PATH
  cd /data/project/chureeratp/RNA/HEART
  python 54_fdr_only_counts.py
  python 54_fdr_only_counts.py --list          # just show which tables it found
  python 54_fdr_only_counts.py --table X.tsv --label "DCM vs HCM"
"""
import argparse
import glob
import os

import pandas as pd

NKEEP = {"SE": 3, "MXE": 5, "A3SS": 5, "A5SS": 5, "RI": 3}

# The three MAGNet contrasts, BY NAME. An earlier version matched filename fragments
# and would have picked A3SS_Colorado_DCMvsNF.all_events.tsv as "MAGNet DCM vs NF" --
# wrong cohort, wrong event class -- because 99 tables share these fragments and the
# glob is alphabetical. Fuzzy matching is fine for finding candidates and unacceptable
# for choosing between them. Override with --dcm-nf etc. if these names are wrong.
DEFAULTS = {
    "MAGNet DCM vs HCM": "SE_DCMvsHCM.all_events.tsv",
    "MAGNet DCM vs NF":  "SE_adjusted.all_events.tsv",   # the base covariate model
    "MAGNet HCM vs NF":  "SE_HCM_v2.all_events.tsv",
}
# tokens that mean a table belongs to another cohort or another analysis entirely
# A table from ANOTHER COHORT can never be a MAGNet contrast -- hard refusal.
FOREIGN = ("colorado", "lbnl", "brigham", "fuwai", "stanford", "probe", "mp_")
# A sensitivity VARIANT of the discovery model might legitimately be the one you want,
# so these only warn. The 926 self-check is the real guard.
VARIANT = ("compadj", "fulladj", "pooladj", "age25", "legacy")
# the discovery table must reproduce this, or it is not the discovery table
EXPECT_DCM_NF_EXON = 926

LADDER = [("FDR only", 0.0), ("+|dPSI|>=0.05", 0.05), ("+|dPSI|>=0.10", 0.10)]


def find_tables(extra):
    seen, out = set(), []
    pats = list(extra) + ["*all_events.tsv", "*adjusted*.tsv", "*.all_events.tsv",
                          "*/all_events.tsv", "*dcm*hcm*.tsv", "*DCM*HCM*.tsv"]
    for p in pats:
        for f in sorted(glob.glob(p)):
            if f not in seen and os.path.isfile(f) and os.path.getsize(f):
                seen.add(f)
                out.append(f)
    return out


def guess_type(t, path):
    """Event class, from a column fingerprint first and the filename only as a fallback."""
    cols = set(t.columns)
    if {"1stExonStart_0base"} & cols:
        return "MXE"
    if {"riExonStart_0base"} & cols:
        return "RI"
    if {"longExonStart_0base"} & cols:
        # A3SS and A5SS share their columns; only the name can separate them
        for et in ("A3SS", "A5SS"):
            if et.lower() in os.path.basename(path).lower():
                return et
        return "A3SS"
    if {"exonStart_0base"} & cols:
        return "SE"
    for et in ("SE", "MXE", "A3SS", "A5SS", "RI"):
        if et.lower() in os.path.basename(path).lower():
            return et
    return None


def load(path):
    try:
        t = pd.read_csv(path, sep="\t")
    except Exception as e:
        return None, None, f"unreadable ({type(e).__name__})"
    d = next((c for c in ("dPSI_adj", "dPSI") if c in t.columns), None)
    if d is None:
        cand = [c for c in t.columns if c.lower().startswith("dpsi")]
        d = next((c for c in cand if "adj" in c.lower()), cand[0] if cand else None)
    if d is None or "FDR" not in t.columns:
        return None, None, "no FDR/dPSI columns"
    if "key" not in t.columns:
        return None, None, "no key column — cannot collapse exon-level"
    return t, d, None


def counts(t, d, et, fdr):
    """Exon-level and full-coordinate counts at each rung of the ladder."""
    out = {}
    n = NKEEP.get(et, 3)
    full = t.drop_duplicates("key")
    ex = t.copy()
    ex["ekey"] = [":".join(str(k).split(":")[:n]) for k in ex["key"]]
    ex = ex.loc[ex.groupby("ekey")["FDR"].idxmin()]
    for name, df in (("exon-level", ex), ("full-coordinate", full)):
        row = {"tested": len(df)}
        for lbl, mind in LADDER:
            sig = (df.FDR < fdr) & (df[d].abs() >= mind) if mind else (df.FDR < fdr)
            row[lbl] = int(sig.sum())
        out[name] = row
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fdr", type=float, default=0.05)
    ap.add_argument("--glob", nargs="*", default=[])
    ap.add_argument("--table", default=None)
    ap.add_argument("--label", default="supplied table")
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--find-discovery", action="store_true",
                    help="try every SE_*.all_events.tsv and report which reproduces "
                         "the frozen 926 exon-level count")
    ap.add_argument("--dcm-hcm", default=DEFAULTS["MAGNet DCM vs HCM"])
    ap.add_argument("--dcm-nf", default=DEFAULTS["MAGNet DCM vs NF"])
    ap.add_argument("--hcm-nf", default=DEFAULTS["MAGNet HCM vs NF"])
    a = ap.parse_args()

    if a.find_discovery:
        print("=" * 78)
        print(f"WHICH TABLE IS THE DISCOVERY CONTRAST?  (looking for "
              f"{EXPECT_DCM_NF_EXON:,} exon-level SE at FDR<{a.fdr}, |dPSI|>=0.10)")
        print("=" * 78)
        hits = []
        for f in sorted(glob.glob("SE_*.all_events.tsv")):
            base = os.path.basename(f).lower()
            if any(t in base for t in FOREIGN):
                continue
            t, d, err = load(f)
            if err:
                print(f"  {os.path.basename(f):<40}{err}")
                continue
            c = counts(t, d, "SE", a.fdr)["exon-level"]
            got = c["+|dPSI|>=0.10"]
            mark = "  <== THIS ONE" if got == EXPECT_DCM_NF_EXON else ""
            if got == EXPECT_DCM_NF_EXON:
                hits.append(f)
            print(f"  {os.path.basename(f):<40}{got:>8,} of {c['tested']:>8,} "
                  f"tested{mark}")
        print()
        if len(hits) == 1:
            print(f"  Use:  --dcm-nf {hits[0]}")
        elif not hits:
            print("  *** NOTHING reproduces the frozen count. Either the convention has")
            print("      drifted or 926 is not what the locked rule gives on any of")
            print("      these tables. Resolve before writing any Results number.")
        else:
            print(f"  *** {len(hits)} tables reproduce it: {', '.join(hits)}")
            print("      Pick the base covariate model, not a sensitivity variant.")
        return

    if a.list:
        files = find_tables(a.glob)
        print(f"{len(files)} candidate table(s):")
        for f in files:
            print(f"  {f}")
        return

    if a.table:
        pairs = [(a.label, a.table)]
    else:
        chosen = {"MAGNet DCM vs HCM": a.dcm_hcm,
                  "MAGNet DCM vs NF": a.dcm_nf,
                  "MAGNet HCM vs NF": a.hcm_nf}
        pairs = []
        for label, path in chosen.items():
            if not (os.path.exists(path) and os.path.getsize(path)):
                print(f"  [none] {label}: {path} not found — pass the right path")
                continue
            base = os.path.basename(path).lower()
            bad = [t for t in FOREIGN if t in base]
            if bad:
                print(f"  *** {label}: {os.path.basename(path)} belongs to another "
                      f"cohort ({', '.join(bad)}). Refusing.")
                continue
            var = [t for t in VARIANT if t in base]
            if var:
                print(f"  [note] {label}: {os.path.basename(path)} is a sensitivity "
                      f"variant ({', '.join(var)}), not the base model. Allowed, but "
                      f"the 926 check below decides.")
            pairs.append((label, path))
        if not pairs:
            print("\n  nothing usable. Run --list, then pass --dcm-hcm/--dcm-nf/--hcm-nf")
            return

    print("=" * 78)
    print(f"COUNTS AT A THRESHOLD LADDER   (FDR < {a.fdr})")
    print("=" * 78)
    rows = []
    for label, path in pairs:
        t, d, err = load(path)
        if err:
            print(f"\n  {label:<20}{os.path.basename(path)}  — {err}")
            continue
        et = guess_type(t, path) or "SE"
        print(f"\n  {label}   [{et}]   {os.path.basename(path)}   (dPSI column: {d})")
        c = counts(t, d, et, a.fdr)
        hdr = f"    {'convention':<18}{'tested':>10}" + \
              "".join(f"{lbl:>16}" for lbl, _ in LADDER)
        print(hdr)
        for conv, r in c.items():
            line = f"    {conv:<18}{r['tested']:>10,}" + \
                   "".join(f"{r[lbl]:>16,}" for lbl, _ in LADDER)
            print(line)
        rows.append(dict(contrast=label, event_type=et, file=path,
                         **{f"{k}_{lbl}": v[lbl] for k, v in c.items()
                            for lbl, _ in LADDER},
                         **{f"{k}_tested": v["tested"] for k, v in c.items()}))

    # Self-check: the discovery table must reproduce the frozen 926 at the locked
    # rule. If it does not, --dcm-nf is pointing at the wrong table and every number
    # above is against the wrong baseline.
    disc = next((r for r in rows if r["contrast"] == "MAGNet DCM vs NF"), None)
    if disc is not None:
        got = disc.get("exon-level_+|dPSI|>=0.10")
        ok = got == EXPECT_DCM_NF_EXON
        print(f"\n  [{'ok ' if ok else '***'}] discovery self-check: "
              f"{os.path.basename(disc['file'])} gives {got:,} exon-level SE at "
              f"FDR<{a.fdr} and |dPSI|>=0.10; the frozen value is "
              f"{EXPECT_DCM_NF_EXON:,}")
        if not ok:
            print("        This is the WRONG discovery table, or the convention has")
            print("        drifted. Do not use any number above until it reconciles.")

    if rows:
        pd.DataFrame(rows).to_csv("fdr_only_counts.tsv", sep="\t", index=False)
        print("\n[write] fdr_only_counts.tsv")

    print("\n" + "=" * 78)
    print("FOR THE DISCUSSION — the Britto-Borges reconciliation")
    print("=" * 78)
    print("""
  Their numbers: 201 differentially expressed genes and 149 differently used
  transcripts, DCM vs HCM, adjusted P <= 0.05, NO effect-size filter, DRIMSeq on
  StringTie2 transcript abundances, 165 DCM / 27 HCM / 162 NFD.

  Use the DCM-vs-HCM "FDR only" exon-level cell above as the comparable figure, NOT
  the thresholded 0. Then make three points, in this order:

  1. Even at FDR alone, put the DCM-vs-HCM count next to the DCM-vs-NF and HCM-vs-NF
     counts in the same rows above. The claim is a RATIO, never an absolute zero.
  2. Differential transcript usage and exon-level PSI are different quantities.
     DRIMSeq asks whether a gene's transcript proportions shift — alternative
     promoters, 3' UTR and last-exon choice, NMD-targeted isoforms all qualify — while
     rMATS asks whether one exon's inclusion shifts. Neither test contains the other.
  3. Their 201 DEGs are expression, about which this paper claims nothing.

  DO NOT write that the two studies disagree, and DO NOT write that this paper is the
  first to contrast DCM with HCM in this cohort. Britto-Borges did it in 2022. The
  honest framing is corroboration with a different measurement, plus the ICM arm, the
  within-cohort positive controls and the genotype screen, which they did not have.
  See prior_art_britto_borges.md.
""")


if __name__ == "__main__":
    main()
