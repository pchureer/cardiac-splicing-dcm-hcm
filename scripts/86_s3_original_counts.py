#!/usr/bin/env python3
"""
86_s3_original_counts.py  —  source S3's "18 versus 7" from the ORIGINAL
                             `lvef_correlation.py` outputs, which survive.

Three files exist, written by `lvef_correlation.py` as `{--out}.lvef_events.tsv`:

    SE_lvef.lvef_events.tsv
    MXE_lvef.lvef_events.tsv
    se_lvef_v2.lvef_events.tsv

TWO THINGS TO ESTABLISH, AND THEY ARE DIFFERENT QUESTIONS
---------------------------------------------------------
1. **Where "18 versus 7" comes from.** It is quoted in the S3 material without a
   source. These files are the source. This script counts FDR<0.05 in each and says
   which pair of numbers reproduces it — or that none does.

2. **How many of those were constant-PSI artefacts.** `RESULT_s3_lvef_null.md` records
   the original hits as artefacts with beta = 0.0000 and rho = nan: a zero-variance PSI
   column makes the regression degenerate, the slope exactly zero and the correlation
   undefined, while the p-value can still come out small. So each "significant" row is
   classified:

     artefact   |beta| < 1e-9, or rho is NaN/non-finite
     genuine    everything else

   A count of significant events is meaningless until the artefacts are removed from it.

3. **THE DENOMINATOR, which is the real finding.** `lvef_correlation.py` prints
   "RBM20-program events tested within DCM", so the original analysis was restricted to
   RBM20-program events. `77_lvef_s3.py` tests 38,084 genome-wide skipped exons. These
   answer DIFFERENT QUESTIONS and the row counts below are what make that concrete.
   The rebuilt S3 broadens the test; it does not merely fix it. Any sentence comparing
   the two must name the event set for each.

GATES
  G1  inventory: every file readable, columns resolved by name, row counts printed.
  G2  does any pair of files reproduce "18 versus 7"? Reported, never assumed.
  G3  are the significant rows artefactual, as the record says? If a file's significant
      rows are NOT all artefacts, say so loudly — the record would then be wrong, and a
      genuine hit in the original analysis is a different story from a bug.

USAGE
    python 86_s3_original_counts.py --list     # columns and row counts only
    python 86_s3_original_counts.py
    python 86_s3_original_counts.py --files A.tsv B.tsv      # explicit
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd

DEFAULT_FILES = ["SE_lvef.lvef_events.tsv",
                 "MXE_lvef.lvef_events.tsv",
                 "se_lvef_v2.lvef_events.tsv"]

# what the record claims
REC_PAIR = (18, 7)
REC_ARTEFACT_N = 7          # "the seven significant events were constant-PSI artefacts"

FDR_CANDIDATES = ("FDR", "fdr", "qval", "q", "padj", "p_adj")
BETA_CANDIDATES = ("LVEF_beta", "lvef_beta", "beta", "slope", "coef")
RHO_CANDIDATES = ("rho", "spearman", "spearman_rho", "r", "corr")
P_CANDIDATES = ("p", "pvalue", "pval", "P", "p_value")
GENE_CANDIDATES = ("gene", "geneSymbol", "symbol")


def pick(cols, cands):
    for c in cands:
        if c in cols:
            return c
    low = {c.lower(): c for c in cols}
    for c in cands:
        if c.lower() in low:
            return low[c.lower()]
    return None


def classify(df, c_beta, c_rho):
    """artefact = zero slope or undefined correlation."""
    n = len(df)
    bad = np.zeros(n, dtype=bool)
    if c_beta is not None:
        b = pd.to_numeric(df[c_beta], errors="coerce").to_numpy(dtype=float)
        bad |= ~np.isfinite(b) | (np.abs(b) < 1e-9)
    if c_rho is not None:
        r = pd.to_numeric(df[c_rho], errors="coerce").to_numpy(dtype=float)
        bad |= ~np.isfinite(r)
    return bad


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default=".")
    ap.add_argument("--files", nargs="*", default=None)
    ap.add_argument("--fdr", type=float, default=0.05)
    ap.add_argument("--list", action="store_true")
    args = ap.parse_args()

    files = args.files or DEFAULT_FILES
    paths = [f if os.path.isabs(f) else os.path.join(args.dir, f) for f in files]

    print("=" * 74)
    print('86 — S3: sourcing "18 versus 7" from the original lvef_correlation outputs')
    print("=" * 74)

    rows = []
    for p in paths:
        name = os.path.basename(p)
        if not os.path.exists(p):
            print(f"\n--- {name}: MISSING")
            continue
        df = pd.read_csv(p, sep="\t")
        c_fdr = pick(df.columns, FDR_CANDIDATES)
        c_beta = pick(df.columns, BETA_CANDIDATES)
        c_rho = pick(df.columns, RHO_CANDIDATES)
        c_p = pick(df.columns, P_CANDIDATES)
        c_gene = pick(df.columns, GENE_CANDIDATES)

        print(f"\n--- {name}")
        print(f"    rows (events tested): {len(df):,d}")
        print(f"    columns: {', '.join(df.columns)}")
        print(f"    resolved: FDR={c_fdr}  beta={c_beta}  rho={c_rho}  p={c_p}")

        if args.list:
            continue
        if c_fdr is None:
            print("    SKIPPED — no FDR column; cannot count significance here.")
            continue

        q = pd.to_numeric(df[c_fdr], errors="coerce").to_numpy(dtype=float)
        sig = df.loc[np.isfinite(q) & (q < args.fdr)]
        bad = classify(sig, c_beta, c_rho)
        n_sig, n_bad = len(sig), int(bad.sum())
        n_ok = n_sig - n_bad
        print(f"    FDR < {args.fdr}: {n_sig}")
        print(f"      constant-PSI artefacts (|beta|<1e-9 or rho undefined): {n_bad}")
        print(f"      surviving as genuine: {n_ok}")
        if n_sig:
            show = sig.copy()
            show["_artefact"] = bad
            cshow = [c for c in (c_gene, c_beta, c_rho, c_p, c_fdr) if c] + ["_artefact"]
            print("      " + show[cshow].head(25).to_string(index=False).replace("\n", "\n      "))
        rows.append(dict(name=name, tested=len(df), sig=n_sig, bad=n_bad, ok=n_ok))

    if args.list or not rows:
        return

    # ------------------------------------------------------------------ G2
    print("\n" + "-" * 74)
    print('G2  does any pair reproduce "18 versus 7"?')
    print("-" * 74)
    counts = {r["name"]: r["sig"] for r in rows}
    print("  FDR<0.05 per file: " + ", ".join(f"{k} = {v}" for k, v in counts.items()))
    hit = [(a, b) for a in rows for b in rows
           if a is not b and (a["sig"], b["sig"]) == REC_PAIR]
    if hit:
        a, b = hit[0]
        print(f"  REPRODUCED — {REC_PAIR[0]} is {a['name']}, {REC_PAIR[1]} is {b['name']}.")
        print(f"  Quote it as: {a['sig']} ({a['name'].split('.')[0]}) versus "
              f"{b['sig']} ({b['name'].split('.')[0]}), on {a['tested']:,d} and "
              f"{b['tested']:,d} events tested respectively.")
    else:
        print(f"  NOT REPRODUCED — no pair of these files gives {REC_PAIR[0]} and "
              f"{REC_PAIR[1]}.")
        print("  The '18 versus 7' in the S3 material does not come from these outputs at")
        print("  this threshold. DO NOT quote it. Either it used a different threshold or")
        print("  a different pair of runs, or it was never measured. Drop the count and")
        print("  say the original hits were constant-PSI artefacts without quantifying.")

    # ------------------------------------------------------------------ G3
    print("\n" + "-" * 74)
    print("G3  were the original hits artefacts, as the record says?")
    print("-" * 74)
    for r in rows:
        if r["sig"] == 0:
            verdict = "no significant rows"
        elif r["ok"] == 0:
            verdict = "ALL artefactual — matches the record"
        elif r["bad"] == 0:
            verdict = "⚠ NONE artefactual — the record is WRONG for this file"
        else:
            verdict = f"⚠ MIXED — {r['bad']} artefact, {r['ok']} genuine"
        print(f"  {r['name']:34s} {r['sig']:>3d} sig  ->  {verdict}")
    anyok = sum(r["ok"] for r in rows)
    if anyok:
        print(f"""
  ⚠ {anyok} significant row(s) across these files are NOT constant-PSI artefacts.
  `RESULT_s3_lvef_null.md` records the original hits as artefacts; that is not true of
  all of them. A genuine PSI-LVEF association in the original restricted event set is a
  different claim from a bug, and the rebuilt genome-wide S3 finding 0 of 38,084 does
  not contradict it — the event sets differ. Resolve before writing the S3 caption.""")
    else:
        print("\n  Every significant row in every file is a constant-PSI artefact. The "
              "record\n  stands, and S3's null is safe to report as a correction.")

    # ------------------------------------------------------------------ scope
    print("\n" + "-" * 74)
    print("THE DENOMINATOR — this is the part that changes what S3 says")
    print("-" * 74)
    for r in rows:
        print(f"  {r['name']:34s} {r['tested']:>8,d} events tested")
    print(f"  {'77_lvef_s3.py (rebuilt, genome-wide SE)':34s} {38084:>8,d} events tested")
    print("""
  `lvef_correlation.py` prints "RBM20-program events tested within DCM": the original
  analysis was RESTRICTED to RBM20-program events. The rebuilt S3 tests genome-wide
  skipped exons. **These are different questions**, so the rebuilt S3 BROADENS the test
  rather than merely fixing it, and the two counts are not a before/after of one
  analysis. Name the event set for each wherever both are mentioned — the Methods, the
  S3 caption, and `RESULT_s3_lvef_null.md`.""")
    print("\n" + "=" * 74)


if __name__ == "__main__":
    main()
