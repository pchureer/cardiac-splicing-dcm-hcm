#!/usr/bin/env python3
"""
95_s1c_column_and_s3_split.py — two supplement jobs, independently gated.

  JOB s1c  Add the cell-type basis column to S1c_markers.tsv.
  JOB s3   Split S3_lvef.tsv (38,084 rows) into a Supplementary Data file and the
           typesettable summary that becomes Supplementary Table 10.

    python 95_s1c_column_and_s3_split.py --list            # inventory, change nothing
    python 95_s1c_column_and_s3_split.py --job s1c         # dry run
    python 95_s1c_column_and_s3_split.py --job s1c --apply
    python 95_s1c_column_and_s3_split.py --job s3 --apply
    python 95_s1c_column_and_s3_split.py --apply           # both

WHAT THE s1c COLUMN IS, AND WHAT IT IS NOT
------------------------------------------
`DECISION_m16_marker_panel.md` settled that **the panels were assembled for this study
and no published panel can honestly be cited.** So this column is NOT a citation column.
It states, for each marker, **what the gene product is and why that makes it
cell-type-restricted** — a checkable factual basis, not a borrowed authority.

The script carries a lookup of short descriptors for the markers it can state
defensibly, and **REFUSES if any marker in the file is not in that lookup**, printing
the misses. It will not invent a basis for a gene it does not know. Fill those by hand
and re-run.

GATES
  s1c G1  the file exists; the gene column resolves; the new column does not already
          exist (--force overrides, with a backup).
  s1c G2  every marker has a lookup entry. Any miss is a refusal, not a blank cell.
  s1c G3  after writing, reload: row count, column set and every pre-existing cell
          unchanged; the new column fully populated.

  s3  G1  the file exists and its columns resolve.
  s3  G2  exactly EXPECT_ROWS_S3 rows (38,084 at 50 bp, 33,130 at 125 bp) — the
          published tested count. A mismatch means this is
          not the table S3 reports, and the split would mislabel it.
  s3  G3  BOTH clauses of the published claim, checked separately: rows at
          FDR < 0.05 (S3_EXPECT_FDR) and rows that ALSO reach
          |dPSI per 10 pp| >= S3_DPSI_CUT (S3_EXPECT_THRESH). 50 bp: 0 and 0.
          125 bp: 1 and 0. The script refuses rather than writing a summary
          that contradicts the paper.
  s3  G4  median |rho| within RHO_TOL of EXPECT_MEDIAN_ABS_RHO
          (0.0559 at 50 bp, 0.0555 at 125 bp).
  s3  G5  duplicate gene symbols among the reported top rows are detected and annotated
          — the TRAPPC14 pair (two overlapping exons 9 bp apart) must not appear as one
          gene listed twice with no explanation.
"""

import argparse
import os
import shutil
import sys
import time

import pandas as pd

S1C = "S1c_markers.tsv"
S3 = os.environ.get("S3", "S3_lvef.tsv")
OUT_DATA = os.environ.get("OUT_DATA", "SupplementaryData_lvef_all_events.tsv")
OUT_SUMMARY = os.environ.get("OUT_SUMMARY", "SupplementaryTable10_lvef_summary.tsv")

NEW_COL = "cell_type_basis"
EXPECT_ROWS_S3 = int(os.environ.get("EXPECT_ROWS_S3", "38_084"))
# S3's claim has TWO clauses -> DECISION_s3_claim_two_clause.md
#   50 bp : 0 at FDR, 0 at the effect-size threshold
#   125 bp: 1 at FDR, 0 at the effect-size threshold
# Both are checked, because the claim asserts both.
S3_EXPECT_FDR = int(os.environ.get("S3_EXPECT_FDR", "0"))
S3_EXPECT_THRESH = int(os.environ.get("S3_EXPECT_THRESH", "0"))
S3_DPSI_CUT = float(os.environ.get("S3_DPSI_CUT", "0.10"))
EXPECT_MEDIAN_ABS_RHO = float(os.environ.get("EXPECT_MEDIAN_ABS_RHO", "0.0559"))
RHO_TOL = 0.0005

C_GENE = ["gene", "marker", "symbol", "Gene", "GeneSymbol", "gene_symbol"]
C_FDR = ["FDR_rank", "FDR_lvef", "FDR_rho", "FDR", "fdr", "padj", "qvalue"]
C_RHO = ["rho_lvef", "rho", "spearman_rho", "rho_LVEF"]
C_BETA = ["beta_lvef", "LVEF_beta", "beta"]

# --------------------------------------------------------------------------
# What each marker is. Textbook, checkable, and deliberately NOT a citation.
# --------------------------------------------------------------------------
BASIS = {
    # composite scores — rows of this table are not all genes
    "CM_SCORE": "composite, not a gene: mean z-score of the cardiomyocyte marker panel "
                "listed below",
    "FB_SCORE": "composite, not a gene: mean z-score of the fibroblast marker panel "
                "listed below",
    # cardiomyocyte
    "MYBPC3": "cardiac myosin-binding protein C; sarcomeric A-band, cardiac-restricted isoform",
    "TCAP":   "telethonin (titin-cap); Z-disc, striated-muscle restricted",
    "RYR2":   "cardiac ryanodine receptor; sarcoplasmic reticulum Ca2+ release channel, cardiac-restricted isoform",
    "PLN":    "phospholamban; SERCA2a regulator at the sarcoplasmic reticulum, cardiomyocyte-enriched",
    "MYL3":   "ventricular/slow essential myosin light chain, cardiomyocyte-restricted",
    "MYL2":   "ventricular regulatory myosin light chain, cardiomyocyte-restricted",
    "CSRP3":  "muscle LIM protein; Z-disc mechanosensor, striated-muscle restricted",
    "TNNI3":  "cardiac troponin I; thin-filament inhibitory subunit, cardiac-restricted isoform",
    "ACTN2":  "alpha-actinin-2; Z-disc crosslinker, striated-muscle restricted",
    "TNNT2":  "cardiac troponin T; thin filament, cardiac-restricted isoform",
    "NPPA":   "atrial natriuretic peptide; EXCLUDED a priori as a fetal-gene-programme marker",
    "NPPB":   "B-type natriuretic peptide; EXCLUDED a priori as a fetal-gene-programme marker",
    "MYH7":   "beta-myosin heavy chain; EXCLUDED a priori as a fetal-gene-programme marker",
    # fibroblast
    "DCN":    "decorin; small leucine-rich interstitial proteoglycan, fibroblast-derived",
    "LUM":    "lumican; interstitial proteoglycan, fibroblast-derived",
    "POSTN":  "periostin; matricellular protein of activated fibroblasts",
    "PDGFRA": "PDGF receptor alpha; canonical cardiac fibroblast surface marker",
    "COL1A1": "collagen I alpha-1; principal fibrillar collagen secreted by fibroblasts",
    "COL1A2": "collagen I alpha-2; principal fibrillar collagen secreted by fibroblasts",
    "COL3A1": "collagen III alpha-1; interstitial fibrillar collagen, fibroblast-derived",
    "COL6A3": "collagen VI alpha-3; pericellular matrix, fibroblast-derived",
    "FBLN1":  "fibulin-1; elastic-fibre-associated matrix glycoprotein, fibroblast-derived",
    "FBLN2":  "fibulin-2; elastic-fibre-associated matrix glycoprotein, fibroblast-derived",
    "VCAN":   "versican; interstitial chondroitin-sulfate proteoglycan, fibroblast-derived",
    "FN1":    "fibronectin 1; matrix glycoprotein, fibroblast-derived",
    "ELN":    "elastin; fibroblast-derived elastic-fibre core protein",
    "MMP2":   "matrix metalloproteinase 2; interstitial matrix turnover, fibroblast-derived",
    "TIMP1":  "TIMP metallopeptidase inhibitor 1; matrix turnover, fibroblast-associated",
    "SERPINF1": "pigment epithelium-derived factor; interstitial, fibroblast-enriched",
    "ABCA8":  "ABCA8 transporter; fibroblast-enriched in human myocardium",
    "C7":     "complement C7; interstitial fibroblast-enriched in human myocardium",
    # control
    "GAPDH":  "housekeeping gene; NOT a cell-type marker — carried as the negative control",
}


def hr(c="-"):
    print(c * 78)


def refuse(msg):
    print()
    hr("=")
    print("REFUSING — " + msg)
    hr("=")
    sys.exit(2)


def pick(df, cands, what, where):
    for c in cands:
        if c in df.columns:
            return c
    refuse(f"no {what} column in {where}\n  tried: {', '.join(cands)}\n"
           f"  has:   {', '.join(map(str, df.columns))}")


def backup(p):
    b = f"{p}.bak-{time.strftime('%Y%m%d-%H%M%S')}"
    shutil.copy2(p, b)
    print(f"  [backup] {b}")
    return b


# ============================================================ s1c ===========
def job_s1c(d, apply_, list_only, force):
    p = os.path.join(d, S1C)
    print()
    hr("=")
    print(f"s1c — cell-type basis column   ({S1C})")
    hr("=")
    if not os.path.exists(p):
        refuse(f"{S1C} not found in {d}")
    df = pd.read_csv(p, sep="\t", dtype=str, keep_default_na=False)
    print(f"  {len(df)} rows · columns: {', '.join(map(str, df.columns))}")

    g = pick(df, C_GENE, "gene", S1C)
    genes = [str(x).strip() for x in df[g]]
    print(f"  gene column: {g!r}")
    print(f"  markers: {', '.join(genes)}")

    if NEW_COL in df.columns and not force:
        print(f"\n  `{NEW_COL}` already exists. Nothing to do (use --force to rewrite).")
        return

    # ---------------------------------------------------------------- G2
    missing = sorted({x for x in genes if x.upper() not in BASIS})
    print(f"\n  G2  lookup coverage: {len(genes) - len(missing)}/{len(genes)}")
    if missing:
        print(f"\n  {len(missing)} row(s) have no entry in this script's lookup:")
        for m in missing:
            hint = ("  <- looks like a COMPOSITE SCORE, not a gene; it needs an entry "
                    "saying so" if m.lower().endswith(("_score", "score")) else "")
            print(f"    {m}{hint}")
        refuse("G2 — this script will not invent a cell-type basis for a gene it does "
               "not\n  know. Add an entry for each marker above to the BASIS dict (one "
               "line each,\n  stating what the gene product is and why it is "
               "cell-type-restricted), or add\n  the column by hand. Blank cells are "
               "not acceptable in a supplement whose\n  whole purpose is to show the "
               "panel was chosen on a principle.")

    vals = [BASIS[x.upper()] for x in genes]
    hr()
    print("  THE COLUMN")
    hr()
    pan = next((c for c in ("panel", "part") if c in df.columns), None)
    for i, (gene, v) in enumerate(zip(genes, vals)):
        tag = f"[{str(df[pan].iloc[i])[:14]:14s}] " if pan else ""
        print(f"    {gene:10s} {tag}{v}")
    print(f"""
  NOTE: this is a BASIS column, not a citation column. The panels were assembled
  for this study and no published panel is cited (DECISION_m16_marker_panel.md),
  so the honest claim is what each gene product is — which a reader can check —
  rather than a borrowed authority.""")

    if list_only or not apply_:
        print("\n  DRY RUN — nothing written.")
        return

    b = backup(p)
    orig = df.copy()
    df[NEW_COL] = vals
    df.to_csv(p, sep="\t", index=False)
    print(f"  [write]  {p}")

    # ---------------------------------------------------------------- G3
    print("\n  G3  reload and verify")
    re_df = pd.read_csv(p, sep="\t", dtype=str, keep_default_na=False)
    ok = len(re_df) == len(orig)
    print(f"    rows {len(re_df)} == {len(orig)}: {ok}")
    same_cols = list(re_df.columns) == list(orig.columns) + [NEW_COL]
    print(f"    columns unchanged plus `{NEW_COL}`: {same_cols}")
    untouched = all((re_df[c] == orig[c]).all() for c in orig.columns)
    print(f"    every pre-existing cell unchanged: {untouched}")
    filled = (re_df[NEW_COL].str.strip() != "").all()
    print(f"    new column fully populated: {filled}")
    if not (ok and same_cols and untouched and filled):
        refuse(f"G3 — the written file is not what was intended. Restore from {b}.")
    print("  G3 PASS")


# ============================================================ s3 ============
def job_s3(d, apply_, list_only, fdr_col):
    p = os.path.join(d, S3)
    print()
    hr("=")
    print(f"s3 — split the LVEF table   ({S3})")
    hr("=")
    if not os.path.exists(p):
        refuse(f"{S3} not found in {d}")
    df = pd.read_csv(p, sep="\t", low_memory=False)
    print(f"  {len(df):,} rows · columns: {', '.join(map(str, df.columns))}")

    g = pick(df, C_GENE, "gene", S3)
    r = pick(df, C_RHO, "rho", S3)
    bcol = next((c for c in C_BETA if c in df.columns), None)

    # ---- THREE FDR columns exist in this file. Do not guess which one the claim used.
    fdr_cols = [c for c in df.columns if c.upper().startswith("FDR")]
    print(f"\n  FDR columns present: {', '.join(fdr_cols) if fdr_cols else 'none'}")
    for c in fdr_cols:
        v = pd.to_numeric(df[c], errors="coerce")
        print(f"    {c:12s} significant at <0.05: {int((v < 0.05).sum()):>5}   "
              f"min {v.min():.4g}")
    if not fdr_cols:
        refuse(f"no FDR column in {S3}.\n  has: {', '.join(map(str, df.columns))}")
    if fdr_col:
        if fdr_col not in fdr_cols:
            refuse(f"--fdr-col {fdr_col!r} is not in the file.\n"
                   f"  choose one of: {', '.join(fdr_cols)}")
        q = fdr_col
    elif len(fdr_cols) == 1:
        q = fdr_cols[0]
    else:
        refuse(f"{len(fdr_cols)} FDR columns and no --fdr-col given.\n"
               f"  S3 reports a COVARIATE-ADJUSTED RANK REGRESSION, so `FDR_rank` is "
               f"almost\n  certainly the published column — but this script will not "
               f"assume it. The\n  counts above tell you which columns give the "
               f"published zero.\n\n  Re-run with, e.g.:  --fdr-col FDR_rank\n\n"
               f"  Whichever you choose, the caption must name the model it belongs to.")
    print(f"\n  resolved: gene={g}  FDR={q}  rho={r}  beta={bcol}")

    # ---------------------------------------------------------------- G2
    print(f"\n  G2  rows {len(df):,} vs published {EXPECT_ROWS_S3:,}")
    if len(df) != EXPECT_ROWS_S3:
        refuse(f"G2 — expected {EXPECT_ROWS_S3:,} rows, found {len(df):,}.\n"
               f"  S3 reports 'no event of {EXPECT_ROWS_S3:,} tested'. If this file has a "
               f"different\n  number of rows it is not the table that claim was made "
               f"from, and splitting\n  it would attach the published caption to the "
               f"wrong data.")
    print("  G2 PASS")

    # ---------------------------------------------------------------- G3
    fdr = pd.to_numeric(df[q], errors="coerce")
    nsig = int((fdr < 0.05).sum())
    _c10 = next((c for c in ("dPSI_per_10pp_lvef", "dPSI_per_10pp")
                 if c in df.columns), None)
    nthr = None
    if _c10:
        _eff = pd.to_numeric(df[_c10], errors="coerce").abs()
        nthr = int(((fdr < 0.05) & (_eff >= S3_DPSI_CUT)).sum())
    print(f"\n  G3  significant at FDR < 0.05: {nsig}   (claim {S3_EXPECT_FDR})")
    if nthr is None:
        print(f"      [warn] no dPSI-per-10pp column; the effect-size clause of the "
              f"claim CANNOT be checked")
    else:
        print(f"      of those, |{_c10}| >= {S3_DPSI_CUT}: {nthr}   "
              f"(claim {S3_EXPECT_THRESH})")
    if nsig != S3_EXPECT_FDR:
        refuse(f"G3 — {nsig} event(s) at FDR < 0.05, but the claim says "
               f"{S3_EXPECT_FDR}.\n"
               f"  Do not write a summary that contradicts the paper. Resolve which is "
               f"right\n  first — the table or the claim.")
    if nthr is not None and nthr != S3_EXPECT_THRESH:
        refuse(f"G3 — {nthr} event(s) at FDR < 0.05 AND |dPSI| >= {S3_DPSI_CUT}, but "
               f"the claim says {S3_EXPECT_THRESH}.\n"
               f"  This is the clause that carries the biology. Resolve it before "
               f"writing.")
    print(f"  G3 PASS — both clauses hold in the data, not only in the text")

    # ---------------------------------------------------------------- G4
    rho = pd.to_numeric(df[r], errors="coerce")
    med = rho.abs().median()
    print(f"\n  G4  median |rho| {med:.4f} vs recorded {EXPECT_MEDIAN_ABS_RHO}")
    if abs(med - EXPECT_MEDIAN_ABS_RHO) > RHO_TOL:
        refuse(f"G4 — median |rho| differs from the record by more than {RHO_TOL}.")
    print("  G4 PASS")

    # ---------------------------------------------------------------- G5
    # ⚠ BH FDR comes in long TIED RUNS (its monotonicity step assigns a run the same
    # value), so ordering by FDR alone makes "the ten lowest" an arbitrary ten out of
    # whichever tie block straddles the cut. Break ties on the underlying p-value.
    pcol = q.replace("FDR", "p").replace("fdr", "p")
    if pcol not in df.columns:
        pcol = next((c for c in df.columns
                     if c.lower().startswith("p_") or c.lower() == "p"), None)
    if pcol is None:
        refuse(f"cannot find the p-value column matching {q!r}; ordering the top rows "
               f"by a tied FDR alone would make the reported list arbitrary.")
    print(f"\n  ordering the reported rows by {pcol} within {q} (BH ties are long runs)")
    order = df.assign(_f=fdr, _p=pd.to_numeric(df[pcol], errors="coerce")) \
              .sort_values(["_f", "_p"], kind="mergesort")
    top = order.head(10).drop(columns=["_f", "_p"]).copy()
    cut = float(pd.to_numeric(order[q], errors="coerce").iloc[9])
    tied = int((fdr == cut).sum())
    if tied > 1:
        print(f"  ⚠ {tied:,} rows share the FDR at the tenth position ({cut:.6g}); "
              f"the p-value\n    ordering is what makes this list reproducible rather "
              f"than arbitrary.")
    dups = top[g].astype(str).duplicated(keep=False)
    print(f"\n  G5  duplicate gene symbols among the ten lowest-FDR rows: {int(dups.sum())}")
    if dups.any():
        for _, row in top[dups].iterrows():
            print(f"    {row[g]}  FDR={row[q]}")
        print("""    These are separate events in one gene, not a gene listed twice.
    The summary below carries an `event_note` column saying so — a top-N table
    that lists one gene twice with no explanation invites exactly one question.""")

    # ---------------------------------------------------------------- build
    note = []
    seen = {}
    for _, row in top.iterrows():
        gn = str(row[g])
        seen[gn] = seen.get(gn, 0) + 1
        note.append("" if list(top[g].astype(str)).count(gn) == 1
                    else f"event {seen[gn]} of {list(top[g].astype(str)).count(gn)} "
                         f"in this gene; separate overlapping exons")
    top["event_note"] = note

    # ------------------------------------------------------------------ G6
    # The null comparator must DERIVE from the data. It used to be the string
    # literal "0.0538", which nothing computed and nothing checked, so it
    # survived the 125 bp re-baselining unchanged.
    #
    # The statistic compared is a MEDIAN |rho|, so the comparator is the
    # null's MEDIAN, not its mean:
    #     rho ~ N(0, 1/(n-1))  under the null
    #     median |rho| = z(0.75) * sigma = 0.6745 / sqrt(n-1)
    #     mean   |rho| = sqrt(2/pi) * sigma = 0.7979 / sqrt(n-1)   <- NOT this
    # At n = 157 those are 0.0540 and 0.0637. Quoting the mean against a
    # measured median would overstate the null by 18%.
    #
    # Var(rho) = 1/(n-1), hence n-1 and not n. df = n gives 0.0538, which is
    # where this file's old literal came from.
    #
    # This expression is character-for-character the one in 77_lvef_s3.py
    # line 687. If one changes, both must.
    _ncol = next((c for c in ("n", "n_used", "n_obs", "nobs", "N")
                  if c in df.columns), None)
    if _ncol is None:
        refuse("no sample-size column in the S3 table, so the analytic null "
               "cannot be derived.\n  A transcribed constant is exactly what "
               "this gate exists to prevent; refusing to\n  write a summary "
               "carrying a number the data cannot support.")
    _nn = pd.to_numeric(df[_ncol], errors="coerce").dropna()
    _n_med = float(_nn.median())
    _n_min, _n_max = float(_nn.min()), float(_nn.max())
    _null_med = 0.6745 / (max(_n_med - 1.0, 2.0) ** 0.5)
    print(f"\n  G6  analytic null for the median |rho|: {_null_med:.4f}"
          f"  (df = n-1, n = {_n_med:g}, range {_n_min:g}-{_n_max:g})")
    if _n_min != _n_max:
        print(f"      n VARIES across events, so this scalar stands for a "
              f"distribution;\n      the deposited table names the n it was "
              f"evaluated at.")
    _exp_null = float(os.environ.get("S3_EXPECT_NULL", "0") or 0)
    if _exp_null > 0:
        # Equality at the DEPOSITED precision (4 dp), not a tolerance. 0.0538
        # and 0.0540 differ by 1.7e-4, so any tolerance at 5e-4 or looser would
        # pass the very disagreement this gate exists to catch.
        if round(_null_med, 4) != round(_exp_null, 4):
            refuse(f"G6 — derived null {_null_med:.4f} differs from the "
                   f"recorded {_exp_null:.4f}.\n  77_lvef_s3.py and this "
                   f"script compute THE SAME QUANTITY and must agree. Two "
                   f"values\n  for one quantity is how 0.0538 and 0.0540 both "
                   f"ended up in the record.")
        print(f"  G6 PASS — derived (not transcribed) and equal to the "
              f"recorded {_exp_null:.4f}")
    else:
        print("  G6 — derived, not transcribed. Set S3_EXPECT_NULL to gate it "
              "against 77's value.")

    summary = pd.DataFrame([
        ("events tested", f"{len(df):,}"),
        ("significant, FDR < 0.05", f"{nsig}"),
        ("median |rho|", f"{med:.4f}"),
        ("analytic null for the median |rho|", f"{_null_med:.4f}"),
        ("  evaluated at n (median across events)",
         f"{_n_med:g} (range {_n_min:g}-{_n_max:g})"),
        ("event set", "genome-wide skipped exons, exon-level, SE_v2filtered"),
        ("model", "rank regression of PSI on LVEF within DCM, adjusted for sex, age, RIN"),
        ("rho quartiles", f"{rho.abs().quantile(.25):.4f} / {med:.4f} / "
                          f"{rho.abs().quantile(.75):.4f}"),
        ("lowest FDR observed", f"{fdr.min():.4g}"),
        ("top rows ordered by", f"{q}, ties broken on {pcol}"),
        ("FDR column used", q),
        ("other FDR columns in the file",
         ", ".join(f"{c} ({int((pd.to_numeric(df[c], errors='coerce') < 0.05).sum())} "
                   f"at <0.05)" for c in fdr_cols if c != q) or "none"),
    ], columns=["quantity", "value"])

    hr()
    print("  SUMMARY TO BE WRITTEN")
    hr()
    for _, row in summary.iterrows():
        print(f"    {row['quantity']:38s} {row['value']}")
    print(f"\n  plus the ten lowest-FDR events with an `event_note` column "
          f"({len(top)} rows)")
    print(f"\n  files:")
    print(f"    {OUT_DATA:42s} full table, {len(df):,} rows -> Supplementary Data")
    print(f"    {OUT_SUMMARY:42s} summary -> Supplementary Table 10")

    if list_only or not apply_:
        print("\n  DRY RUN — nothing written.")
        return

    dp, sp = os.path.join(d, OUT_DATA), os.path.join(d, OUT_SUMMARY)
    for f in (dp, sp):
        if os.path.exists(f):
            backup(f)
    df.to_csv(dp, sep="\t", index=False)
    with open(sp, "w") as fh:
        summary.to_csv(fh, sep="\t", index=False)
        fh.write("\n")
        top.to_csv(fh, sep="\t", index=False)
    print(f"\n  [write]  {dp}")
    print(f"  [write]  {sp}")

    re_full = pd.read_csv(dp, sep="\t", low_memory=False)
    print(f"\n  verify: {OUT_DATA} reloads at {len(re_full):,} rows "
          f"{'ok' if len(re_full) == len(df) else '⚠ WRONG'}")
    if len(re_full) != len(df):
        refuse("the written data file does not round-trip.")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default=".")
    ap.add_argument("--job", choices=["s1c", "s3"], help="one job only (default: both)")
    ap.add_argument("--list", action="store_true", dest="list_only")
    ap.add_argument("--apply", action="store_true")
    ap.add_argument("--fdr-col", dest="fdr_col",
                    help="s3: which FDR column the published claim uses "
                         "(the file carries more than one)")
    ap.add_argument("--force", action="store_true",
                    help="s1c: rewrite the column if it already exists")
    a = ap.parse_args()
    if a.list_only and a.apply:
        refuse("--list and --apply are contradictory.")

    hr("=")
    print("95 — S1(c) cell-type basis column · S3 split")
    hr("=")
    print(f"  dir   {os.path.abspath(a.dir)}")
    print(f"  mode  {'--list' if a.list_only else ('--apply' if a.apply else 'dry run')}")

    if a.job in (None, "s1c"):
        job_s1c(a.dir, a.apply, a.list_only, a.force)
    if a.job in (None, "s3"):
        job_s3(a.dir, a.apply, a.list_only, a.fdr_col)

    print()
    hr("=")
    print("NEXT")
    hr("=")
    print(f"""
  1. S1(c): the basis column documents WHY each gene marks its cell type. It is not
     a citation — the panels were assembled for this study and no published panel is
     cited. Read the column before depositing; a wrong descriptor is worse than none.
  2. S3: the full {EXPECT_ROWS_S3:,}-row table goes to Supplementary Data, the summary
     becomes Supplementary Table 10. The caption is written and ready in
     RESULT_s3_original_counts_sourced.md — it quotes no count of the original
     restricted analysis, by decision.
  3. Neither job changes a published number. If anything in the manuscript moves after
     running this, stop and find out why.
""")


if __name__ == "__main__":
    main()
