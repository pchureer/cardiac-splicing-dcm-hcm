#!/usr/bin/env python3
"""
70_marker_panel_audit.py -- is the cell-composition score partly a DISEASE score?

THE QUESTION, AND WHY IT IS NOT RHETORICAL
------------------------------------------
The cardiomyocyte composition score is the mean z-score of nine markers: TNNI3, MYL2,
MYL3, ACTN2, TCAP, CSRP3, MYBPC3, RYR2, PLN. Eight of those nine appear in this
project's own `cardiomyopathy_genes.bed` as cardiomyopathy disease genes with ClinGen
Definitive or Moderate classifications.

NPPA, NPPB and MYH7 were deliberately EXCLUDED from the panel because they are "strongly
disease-regulated and would make the composition score partly a disease score" (Methods).
MYH7 sits in that BED with the same Definitive-HCM / Definitive-DCM status as TNNI3,
MYBPC3, MYL2 and MYL3 -- which did NOT get excluded. So the exclusion rule was applied by
hand and it is not obvious it went far enough.

**A disease GENE is not the same as a disease-REGULATED gene.** Mutations in MYBPC3 cause
HCM without its tissue expression necessarily shifting in failing myocardium, and these
structural genes are conventional cardiomyocyte identity markers. The overlap is a reason
to measure, not a verdict.

WHAT IS AT STAKE
----------------
S1's argument is that the ~36% of events losing significance under composition adjustment
are THRESHOLD ATTRITION, not confounding -- evidenced by uniform standard-error inflation
(1.214 kept against 1.215 lost) and lost events keeping 74-82% of their effect.

If the CM score is itself strongly differential between DCM and NF, then adjusting for it
partially adjusts for DISEASE, and some of that loss is OVER-adjustment -- removing real
signal rather than failing marginal events. The attrition argument would still stand on
the uniform inflation, but the interpretation changes and the Limitations would need a
sentence. If the score is NOT differential, the objection is closed with a number.

Either way the result goes in the Methods or the Limitations, because a reviewer who
notices the overlap will ask.

THE REAL FILE LAYOUTS, pinned 6 Oct from --list on the cluster
--------------------------------------------------------------
Both of v1's assumptions were wrong, and the correction splits the audit in two.

* `magnet_expr.per_sample.tsv` is **SAMPLES x FEATURES**, not genes x samples, and it
  holds only the regulator genes plus **CM_score, FB_score, libsize, log10_libsize,
  TIN**. The sixteen MARKER genes are not in it.
* `magnet_covariates.tsv` has **NO HEADER ROW** -- six positional columns:
  run, group, sex, age, RIN, pool. Read with header=None or the first sample becomes
  the column names.

So:

  PART A -- the SCORES.  Answerable from this file alone, and it is the half that
            matters: CM_score and FB_score are the actual regressors the model used,
            so testing them tests the adjustment itself rather than a reconstruction.
  PART B -- the per-MARKER audit.  Needs the full log2CPM gene matrix that script 02
            built (the one RBM5/ZRANB2/HNRNPF/HNRNPA0 were pulled from). Pass it with
            --gene-matrix. If it is absent, Part A still runs and Part B says so
            rather than guessing.

THREE GATES, not one
--------------------
`RESULT_rbm20_mechanism_null.md` records RBM20's disease effect under three nested
models, and every term is in this file:

    unadjusted                                  +0.243
    + age/sex/RIN/libsize/TIN   (technical)     +0.382
    + CM_score/FB_score         (composition)   +0.340

All three must reproduce. One gate can pass by luck; three nested ones passing means
the matrix, the covariate join, the disease coding and the model are all the ones the
recorded numbers came from.

WHAT THIS SCRIPT DOES
---------------------
Per marker gene, and for the two composition scores themselves, it fits

    log2CPM ~ disease + sex + age + RIN            (the discovery model)
    log2CPM ~ disease + sex + age + RIN + log10_libsize + TIN   (the technical model)

and reports the disease effect in SD units of expression, with t, p and a
Benjamini-Hochberg FDR across the panel.

SELF-CHECK: RBM20 must reproduce **+0.382 SD** under the technical model, which is the
recorded value in `RESULT_rbm20_mechanism_null.md` (and +0.243 unadjusted). If it does
not, the model or the matrix is not the one those numbers came from, and the script
refuses rather than reporting marker effects nobody can place.

THE PANEL IS READ FROM SCRIPT 02, NOT RE-TYPED
----------------------------------------------
Re-hard-coding the sixteen genes here would let the two copies drift, and a drifted panel
would silently audit the wrong thing. The script parses them out of
`02_build_expr_matrix.py` and refuses if it cannot. The fallback list exists only to be
compared against what it parsed.

USAGE
-----
    conda activate rnaseq_tools
    cd /data/project/chureeratp/RNA/HEART

    python 70_marker_panel_audit.py --list     # find inputs, show schemas and the panel
    python 70_marker_panel_audit.py            # fit, gate, deposit
"""
import argparse
import glob
import os
import re

import numpy as np
import pandas as pd

try:
    from scipy import stats as _st

    def _sf(t, df):
        return 2.0 * _st.t.sf(np.abs(t), df)
    _PMETHOD = "scipy t-distribution"
except Exception:                                                    # noqa: BLE001
    def _sf(t, df):
        # normal approximation; at n>300 the difference is in the 4th decimal
        from math import erfc, sqrt
        return erfc(abs(t) / sqrt(2.0))
    _PMETHOD = "normal approximation (scipy unavailable)"

# ------------------------------------------------------------------ the record
RECORDED_RBM20 = {"unadjusted": 0.243, "technical": 0.382, "composition": 0.340}
RECORDED_RBM20_TOL = 0.02
# magnet_covariates.tsv has NO header: six positional columns, pinned 6 Oct
COV_COLS = ["run", "group", "sex", "age", "RIN", "pool"]
CASE_LABEL, CTRL_LABEL = "DCM", "NF"
SCORE_COLS = ("CM_score", "FB_score")
SOURCE_SCRIPT = "02_build_expr_matrix.py"
BED = "cardiomyopathy_genes.bed"

# Comparison only -- the live panel is parsed from SOURCE_SCRIPT.
FALLBACK_CM = ["TNNI3", "MYL2", "MYL3", "ACTN2", "TCAP", "CSRP3", "MYBPC3", "RYR2", "PLN"]
FALLBACK_FB = ["COL1A1", "COL1A2", "COL3A1", "POSTN", "FN1", "DCN", "LUM"]
# excluded by hand, and the reason given in the Methods
EXCLUDED = ["NPPA", "NPPB", "MYH7"]

PATTERNS = [
    ("expression matrix", ["magnet_expr*.tsv", "*expr*matrix*.tsv", "*log2cpm*.tsv",
                           "*per_sample*.tsv"]),
    ("covariates", ["magnet_covariates.tsv", "*covariate*.tsv", "*pheno*.tsv",
                    "*.pheno.tsv"]),
    ("full gene matrix (for Part B)", ["*log2cpm*.tsv", "*gene*matrix*.tsv",
                                       "*expr*matrix*.tsv", "*ReadsPerGene*matrix*",
                                       "*counts*matrix*.tsv", "magnet_gene*.tsv"]),
    ("disease gene panel", [BED]),
    ("the panel's source script", [SOURCE_SCRIPT]),
]


# ------------------------------------------------------------- panel discovery
def parse_panel():
    """
    Pull the marker lists out of 02_build_expr_matrix.py rather than re-typing them.

    Looks for list literals assigned to a name containing CM/CARDIO/MYO or FB/FIBRO.
    Returns (cm, fb, note). Either list may be None, which is a refusal condition.
    """
    if not os.path.exists(SOURCE_SCRIPT):
        return None, None, f"{SOURCE_SCRIPT} not found"
    src = open(SOURCE_SCRIPT).read()
    cm = fb = None
    # NAME = [ "A", "B", ... ]   across one or more lines
    for m in re.finditer(r"(\w+)\s*=\s*\[([^\]]*)\]", src, re.S):
        name, body = m.group(1).upper(), m.group(2)
        genes = re.findall(r"[\"']([A-Z0-9\-]{2,15})[\"']", body)
        if len(genes) < 3:
            continue
        if cm is None and any(k in name for k in ("CM", "CARDIO", "MYO")):
            cm = genes
        elif fb is None and any(k in name for k in ("FB", "FIBRO")):
            fb = genes
    note = []
    if cm is None:
        note.append("could not parse a cardiomyocyte list")
    if fb is None:
        note.append("could not parse a fibroblast list")
    return cm, fb, "; ".join(note)


def read_bed_genes():
    if not os.path.exists(BED):
        return {}, f"{BED} not found"
    out = {}
    for line in open(BED):
        if line.startswith("#") or not line.strip():
            continue
        f = line.split()
        if len(f) >= 6:
            out[f[3]] = (f[4], f[5])          # gene -> (category, clingen_class)
    return out, None


# -------------------------------------------------------------------- the fits
def zscore(x):
    x = np.asarray(x, dtype=float)
    sd = np.nanstd(x, ddof=1)
    return (x - np.nanmean(x)) / sd if sd > 0 else x * 0.0


def fit(y, X, names):
    """OLS; returns the disease coefficient in SD units of y, its t, p and dof."""
    ok = np.isfinite(y) & np.isfinite(X).all(axis=1)
    y, X = y[ok], X[ok]
    n, p = X.shape
    if n <= p + 1:
        return None
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    dof = n - p
    s2 = resid @ resid / dof
    try:
        XtXi = np.linalg.inv(X.T @ X)
    except np.linalg.LinAlgError:
        return None
    se = np.sqrt(np.diag(XtXi) * s2)
    i = names.index("disease")
    sdy = np.std(y, ddof=1)
    if sdy == 0 or se[i] == 0:
        return None
    t = beta[i] / se[i]
    return dict(beta_sd=beta[i] / sdy, t=t, p=_sf(t, dof), n=n, dof=dof)


def bh(p):
    p = np.asarray(p, dtype=float)
    keep = np.isfinite(p)
    q = np.full(p.shape, np.nan)
    idx = np.where(keep)[0]
    if idx.size == 0:
        return q
    o = idx[np.argsort(p[idx])]
    m = o.size
    prev = 1.0
    for rank in range(m - 1, -1, -1):
        prev = min(prev, p[o[rank]] * m / (rank + 1))
        q[o[rank]] = prev
    return q


# -------------------------------------------------------------------- plumbing
def describe(path, nmax=3):
    try:
        t = pd.read_csv(path, sep="\t", low_memory=False, nrows=200)
    except Exception as e:                                           # noqa: BLE001
        print(f"      {path}  [unreadable: {e.__class__.__name__}]"); return
    print(f"      {path}  ({len(t.columns)} columns)")
    cols = list(map(str, t.columns))
    print(f"        columns: {', '.join(cols[:18])}"
          + (f" ... (+{len(cols) - 18} more)" if len(cols) > 18 else ""))
    for line in t.head(nmax).iloc[:, :8].to_string(index=False).splitlines():
        print(f"          {line}")


def do_list():
    print("\n" + "=" * 78)
    print("INPUTS -- what the marker audit needs")
    print("=" * 78)
    for label, pats in PATTERNS:
        hits = sorted({p for pat in pats for p in glob.glob(pat) if os.path.isfile(p)})
        print(f"\n  {label}:")
        if not hits:
            print("      *** nothing matched ***")
        for p in hits[:6]:
            if p.endswith((".py", ".bed")):
                print(f"      {p}  ({os.path.getsize(p):,} bytes)")
            else:
                describe(p)
        if len(hits) > 6:
            print(f"      ... and {len(hits) - 6} more")

    cm, fb, note = parse_panel()
    print("\n" + "=" * 78)
    print(f"THE PANEL, parsed from {SOURCE_SCRIPT}")
    print("=" * 78)
    print(f"  cardiomyocyte: {cm}")
    print(f"  fibroblast:    {fb}")
    if note:
        print(f"  *** {note} ***")
    if cm and sorted(cm) != sorted(FALLBACK_CM):
        print(f"  [WARN ] parsed CM panel differs from the recorded one {FALLBACK_CM}")
    if fb and sorted(fb) != sorted(FALLBACK_FB):
        print(f"  [WARN ] parsed FB panel differs from the recorded one {FALLBACK_FB}")

    bedg, berr = read_bed_genes()
    if berr:
        print(f"\n  {berr} -- the overlap table cannot be built")
    else:
        print(f"\n  OVERLAP with {BED} ({len(bedg)} disease genes):")
        for lab, panel in (("CM", cm or FALLBACK_CM), ("FB", fb or FALLBACK_FB)):
            hit = [g for g in panel if g in bedg]
            print(f"    {lab}: {len(hit)} of {len(panel)} are cardiomyopathy disease genes")
            for g in hit:
                print(f"       {g:<9} {bedg[g][0]:<14} {bedg[g][1]}")
        print(f"    excluded by hand: "
              + ", ".join(f"{g}{' [IN BED: ' + bedg[g][1] + ']' if g in bedg else ''}"
                          for g in EXCLUDED))
    print("\n  Pin --expr and --cov from the list above, then run without --list.\n")


def load_cov(path):
    """
    magnet_covariates.tsv has NO header row. Reading it with the default header=0
    silently turns the first sample into the column names -- which is exactly what
    happened on the first run, and the refusal messages then listed a sample accession
    as an available column.

    Detect it: if the first cell of the header looks like a run accession, re-read with
    header=None and apply the pinned positional names.
    """
    probe = pd.read_csv(path, sep="\t", nrows=1, dtype=str)
    first = str(probe.columns[0])
    headerless = bool(re.match(r"^(SRR|ERR|SAMEA|DRR)\d+$", first))
    if headerless:
        C = pd.read_csv(path, sep="\t", header=None, dtype=str)
        if C.shape[1] != len(COV_COLS):
            return None, (f"{path} is headerless with {C.shape[1]} columns; the pinned "
                          f"layout has {len(COV_COLS)}: {COV_COLS}")
        C.columns = COV_COLS
        return C, f"headerless, positional names applied: {COV_COLS}"
    C = pd.read_csv(path, sep="\t", dtype=str)
    return C, "has a header row"


def load_expr(path):
    """
    Returns (wide_dataframe_indexed_by_sample, note).

    magnet_expr.per_sample.tsv is SAMPLES x FEATURES with a `sample` column. A
    genes-x-samples matrix is transposed into the same shape so Part B can use either.
    """
    E = pd.read_csv(path, sep="\t", low_memory=False)
    c0 = str(E.columns[0])
    looks_wide = c0.lower() in ("sample", "sample_id", "run", "run_accession") or bool(
        re.match(r"^(SRR|ERR|SAMEA|DRR)\d+$", str(E.iloc[0, 0])))
    if looks_wide:
        E = E.set_index(E[c0].astype(str)).drop(columns=[c0])
        return E.apply(pd.to_numeric, errors="coerce"), f"wide (samples x features), index={c0}"
    E = E.drop_duplicates(subset=[c0]).set_index(E[c0].astype(str)).drop(columns=[c0])
    return (E.apply(pd.to_numeric, errors="coerce").T,
            f"long (features x samples), transposed; feature column={c0}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--expr", default="magnet_expr.per_sample.tsv")
    ap.add_argument("--cov", default="magnet_covariates.tsv")
    ap.add_argument("--gene-matrix", default=None,
                    help="full log2CPM gene matrix, for PART B (the per-marker audit)")
    ap.add_argument("--out", default="marker_panel_audit.tsv")
    args = ap.parse_args()

    if args.list:
        do_list()
        return

    print("\n" + "=" * 78)
    print("MARKER PANEL AUDIT -- is the composition score partly a disease score?")
    print("=" * 78)
    print(f"  p-values: {_PMETHOD}")

    cm, fb, note = parse_panel()
    if note:
        print(f"\n  REFUSING: {note}\n")
        return

    for pth in (args.expr, args.cov):
        if not os.path.exists(pth):
            print(f"\n  REFUSING: {pth} not found. Run --list.\n")
            return

    E, enote = load_expr(args.expr)
    C, cnote = load_cov(args.cov)
    if C is None:
        print(f"\n  REFUSING: {cnote}\n")
        return
    print(f"\n  expression: {args.expr}   {enote}")
    print(f"  covariates: {args.cov}   {cnote}")

    C = C.set_index(C["run"].astype(str))
    samples = [s for s in E.index if s in C.index]
    print(f"  samples matched: {len(samples)} of {len(E.index)} expression rows, "
          f"{len(C.index)} covariate rows")
    if len(samples) < 50:
        print("\n  REFUSING: too few samples matched. Check --expr / --cov.\n")
        return

    cv = C.loc[samples]
    grp = cv["group"].astype(str).str.upper()
    counts = grp.value_counts().to_dict()
    print(f"  groups: {counts}")
    keep = grp.isin([CASE_LABEL, CTRL_LABEL])
    if keep.sum() < 50:
        print(f"\n  REFUSING: fewer than 50 samples labelled {CASE_LABEL}/{CTRL_LABEL}.\n")
        return
    if (~keep).sum():
        print(f"  dropping {int((~keep).sum())} samples not labelled "
              f"{CASE_LABEL} or {CTRL_LABEL}")
    samples = [s for s, k in zip(samples, keep) if k]
    cv, grp = C.loc[samples], grp[keep]
    X_expr = E.loc[samples]

    dis = (grp == CASE_LABEL).astype(float).to_numpy()
    sex = cv["sex"].astype(str).str.upper().str.startswith("M").astype(float).to_numpy()
    age = pd.to_numeric(cv["age"], errors="coerce").to_numpy()
    rin = pd.to_numeric(cv["RIN"], errors="coerce").to_numpy()
    print(f"  coding: {CASE_LABEL}=1 (n={int(dis.sum())}), "
          f"{CTRL_LABEL}=0 (n={int((1-dis).sum())}); sex Male=1")
    nrin = int(np.isnan(rin).sum())
    if nrin:
        print(f"  {nrin} samples lack RIN and drop out of the adjusted models "
              f"(the recorded effective n is 165/161)")

    def col(name):
        if name not in X_expr.columns:
            return None
        return X_expr[name].to_numpy(float)

    lib, tin = col("log10_libsize"), col("TIN")
    missing = [n for n, v in (("log10_libsize", lib), ("TIN", tin)) if v is None]
    if missing:
        print(f"\n  REFUSING: {', '.join(missing)} absent, so the recorded "
              f"+{RECORDED_RBM20['technical']} gate cannot be evaluated.\n")
        return
    cmv, fbv = col(SCORE_COLS[0]), col(SCORE_COLS[1])
    if cmv is None or fbv is None:
        print(f"\n  REFUSING: {SCORE_COLS} absent -- these ARE the audit's subject.\n")
        return

    MODELS = {
        "unadjusted":  [("disease", dis)],
        "technical":   [("disease", dis), ("sex", sex), ("age", age), ("rin", rin),
                        ("lib", lib), ("tin", tin)],
        "composition": [("disease", dis), ("sex", sex), ("age", age), ("rin", rin),
                        ("lib", lib), ("tin", tin), ("cm", cmv), ("fb", fbv)],
    }

    def design(terms):
        names = ["const"] + [n for n, _ in terms]
        return np.column_stack([np.ones(len(samples))]
                               + [np.asarray(v, float) for _, v in terms]), names

    DES = {k: design(v) for k, v in MODELS.items()}

    # ---------------------------------------------------- THREE GATES
    print("\n  GATES -- RBM20's disease effect under three nested models")
    if "RBM20" not in X_expr.columns:
        print("  REFUSING: RBM20 not in the matrix.\n")
        return
    y = X_expr["RBM20"].to_numpy(float)
    failed = []
    for mname, rec in RECORDED_RBM20.items():
        Xm, nm = DES[mname]
        r = fit(y, Xm, nm)
        if r is None:
            failed.append(f"{mname}: fit failed"); continue
        d = abs(r["beta_sd"] - rec)
        ok = d <= RECORDED_RBM20_TOL
        print(f"    [{'ok     ' if ok else 'GATE   '}] {mname:<12} "
              f"{r['beta_sd']:+.3f} against the recorded {rec:+.3f}   n={r['n']}")
        if not ok:
            failed.append(f"{mname}: {r['beta_sd']:+.3f} against {rec:+.3f}")
    if failed:
        print("\n  REFUSING -- the model is not the one the record came from:")
        for f_ in failed:
            print(f"    * {f_}")
        print("  Three nested gates; any failure means the matrix, the covariate join,")
        print("  the disease coding or the model differs. Marker effects computed under")
        print("  a model that cannot reproduce a known value are not depositable.\n")
        return
    print("    all three reproduce -- matrix, join, coding and model are the recorded ones")

    bedg, _ = read_bed_genes()
    rows = []

    def add(panel, gene, y, extra=None):
        d = dict(part="", panel=panel, gene=gene, disease_gene="-", clingen="")
        for mname in ("unadjusted", "technical"):
            Xm, nm = DES[mname]
            r = fit(y, Xm, nm)
            d[f"beta_sd_{mname}"] = None if r is None else round(r["beta_sd"], 4)
            d[f"p_{mname}"] = None if r is None else r["p"]
            d["n"] = None if r is None else r["n"]
        d.update(extra or {})
        rows.append(d)

    # ---------------------------------------------------- PART A: the scores
    print("\n" + "=" * 78)
    print("PART A -- the composition SCORES themselves (the regressors the model used)")
    print("=" * 78)
    for nm_ in SCORE_COLS:
        add("score", nm_, X_expr[nm_].to_numpy(float), {"part": "A"})
    # GAPDH as an internal negative control on the machinery
    if "GAPDH" in X_expr.columns:
        add("control", "GAPDH", X_expr["GAPDH"].to_numpy(float), {"part": "A"})

    # ------------------------------------------- do the scores behave like fractions?
    cmz, fbz = X_expr[SCORE_COLS[0]].to_numpy(float), X_expr[SCORE_COLS[1]].to_numpy(float)
    ok = np.isfinite(cmz) & np.isfinite(fbz)
    r_cf = float(np.corrcoef(cmz[ok], fbz[ok])[0, 1]) if ok.sum() > 10 else float("nan")
    print(f"\n  corr(CM_score, FB_score) = {r_cf:+.3f}  over {int(ok.sum())} samples")
    if r_cf > 0.2:
        print("  >>> POSITIVE. Cardiomyocyte and fibroblast fractions should TRADE OFF:")
        print("      more fibrosis means proportionally fewer myocytes, so the two scores")
        print("      ought to be NEGATIVELY correlated. A positive correlation means they")
        print("      are tracking something global -- library composition, RNA integrity,")
        print("      or overall transcriptome concentration -- rather than cell fraction.")
        print("      If so, 'composition score' is the wrong name for the covariate, and")
        print("      the Methods must describe it as what it measures.")
    elif r_cf < -0.2:
        print("  >>> NEGATIVE, as cell fractions should be. The scores behave like")
        print("      competing compartments, which supports reading them as composition.")
    else:
        print("  >>> NEAR ZERO. The two compartments are not trading off, which is weak")
        print("      evidence either way; report it rather than assuming.")

    # ---------------------------------------------------- PART B: the markers
    print()
    gm = args.gene_matrix
    partB = False
    if gm and os.path.exists(gm):
        if os.path.abspath(gm) == os.path.abspath(args.out):
            print(f"  PART B REFUSED: --gene-matrix is this script's own OUTPUT file")
            print(f"  ({args.out}). Pass the full log2CPM gene matrix instead.")
            G = None
        else:
            try:
                G, gnote = load_expr(gm)
            except Exception as e:                                   # noqa: BLE001
                print(f"  PART B REFUSED: {gm} could not be read as an expression "
                      f"matrix ({e.__class__.__name__}).")
                print("  It needs one row per gene and one column per sample, or the")
                print("  transpose. Check it with: head -2 " + gm)
                G = None
        if G is not None and (G.shape[0] < 50 or G.shape[1] < 100):
            print(f"  PART B REFUSED: {gm} is {G.shape[0]} x {G.shape[1]} after "
                  f"orientation -- too small to be a gene matrix.")
            G = None
    if gm and os.path.exists(gm) and G is not None:
        print(f"  PART B gene matrix: {gm}   {gnote}")
        gs = [s for s in samples if s in G.index]
        if len(gs) < 50:
            print(f"  PART B SKIPPED: only {len(gs)} samples matched in {gm}")
        else:
            partB = True
            Gs = G.loc[gs]
            # rebuild the design on this sample subset
            idx = [samples.index(s) for s in gs]
            sub = {k: ([(n, np.asarray(v, float)[idx]) for n, v in v2])
                   for k, v2 in MODELS.items()}
            DESB = {k: (np.column_stack([np.ones(len(gs))]
                                        + [v for _, v in t]),
                        ["const"] + [n for n, _ in t]) for k, t in sub.items()}
            for panel, genes in (("cardiomyocyte", cm), ("fibroblast", fb),
                                 ("EXCLUDED by hand", EXCLUDED)):
                for g_ in genes:
                    if g_ not in Gs.columns:
                        rows.append(dict(part="B", panel=panel, gene=g_,
                                         note="not in gene matrix"))
                        continue
                    yy = Gs[g_].to_numpy(float)
                    d = dict(part="B", panel=panel, gene=g_,
                             disease_gene=("yes" if g_ in bedg else "no"),
                             clingen=(bedg[g_][1] if g_ in bedg else ""))
                    for mname in ("unadjusted", "technical"):
                        Xm, nm = DESB[mname]
                        r = fit(yy, Xm, nm)
                        d[f"beta_sd_{mname}"] = None if r is None else round(r["beta_sd"], 4)
                        d[f"p_{mname}"] = None if r is None else r["p"]
                        d["n"] = None if r is None else r["n"]
                    rows.append(d)
    else:
        print("  PART B NOT RUN -- the per-marker audit needs the full log2CPM gene")
        print("  matrix. `magnet_expr.per_sample.tsv` carries only the regulator genes")
        print("  plus CM_score/FB_score; the sixteen markers are not in it.")
        print("  Find the matrix script 02 built (the one RBM5/ZRANB2/HNRNPF/HNRNPA0")
        print("  were pulled from) and pass it:")
        print("      python 70_marker_panel_audit.py --gene-matrix <that file>")
        print("  Part A below is the half that matters and does not need it.")

    df = pd.DataFrame(rows)
    for c in ("p_unadjusted", "p_technical"):
        if c in df.columns:
            df[c.replace("p_", "fdr_")] = bh(df[c].to_numpy(dtype=float))

    # ---------------------------------------------------- report
    print("\n" + "-" * 78)
    print(f"  {'part':<6}{'panel':<18}{'gene':<12}{'dis.gene':<9}"
          f"{'unadj':>9}{'tech':>9}{'FDR tech':>11}")
    print("-" * 78)
    for _, r in df.iterrows():
        if r.get("note"):
            print(f"  {r['part']:<6}{r['panel']:<18}{r['gene']:<12}{r['note']}")
            continue
        bt, bu = r.get("beta_sd_technical"), r.get("beta_sd_unadjusted")
        q = r.get("fdr_technical")
        star = " *" if (pd.notna(q) and q < 0.05) else ""
        print(f"  {r['part']:<6}{r['panel']:<18}{r['gene']:<12}"
              f"{str(r.get('disease_gene', '-')):<9}"
              f"{(f'{bu:+.3f}' if pd.notna(bu) else '-'):>9}"
              f"{(f'{bt:+.3f}' if pd.notna(bt) else '-'):>9}"
              f"{(f'{q:.3g}' if pd.notna(q) else '-'):>11}{star}")
    print("\n  * FDR < 0.05 across everything audited here, technical model\n")

    # ---------------------------------------------------- verdict
    print("=" * 78)
    print("VERDICT")
    print("=" * 78)
    for nm_ in SCORE_COLS:
        r = df[df["gene"] == nm_].iloc[0]
        q = r["fdr_technical"]
        b = r["beta_sd_technical"]
        tag = ("IS differential" if (pd.notna(q) and q < 0.05) else "is NOT differential")
        print(f"\n  {nm_}: {b:+.3f} SD, FDR {q:.3g} -- {tag} between "
              f"{CASE_LABEL} and {CTRL_LABEL}")
    if "GAPDH" in set(df["gene"]):
        g = df[df["gene"] == "GAPDH"].iloc[0]
        print(f"\n  GAPDH (internal control on the machinery): "
              f"{g['beta_sd_technical']:+.3f} SD, FDR {g['fdr_technical']:.3g}")

    worst = df[df["gene"].isin(SCORE_COLS)]["beta_sd_technical"].abs().max()
    anysig = ((df[df["gene"].isin(SCORE_COLS)]["fdr_technical"] < 0.05).any())
    print()
    if not anysig:
        print("  >>> NEITHER SCORE is differentially expressed between groups. The")
        print("      composition adjustment is adjusting for cell content, not for")
        print("      disease, so S1's loss is not over-adjustment. This belongs in the")
        print(f"      METHODS as one sentence, with the largest effect ({worst:.3f} SD)")
        print("      quoted so a reviewer can see the margin.")
    else:
        print("  >>> AT LEAST ONE SCORE IS DIFFERENTIAL -- and this is NOT a flaw in")
        print("      the panel. Cell composition genuinely differs between failing and")
        print("      non-failing myocardium: fibrosis and cardiomyocyte loss ARE")
        print("      remodelling. A score that did not differ would be the suspicious")
        print("      result, because it would not be measuring composition.")
        print()
        print("      What it means is that composition is at least partly a MEDIATOR of")
        print("      disease, not an independent technical nuisance. Adjusting for a")
        print("      mediator does not correct confounding -- it asks a narrower")
        print("      question (does splicing differ BEYOND what the composition shift")
        print("      explains) and removes real disease signal along with it.")
        print()
        print("      So S1's ~36% loss is part threshold attrition and part CONSERVATIVE")
        print("      over-adjustment. The uniform SE inflation (1.214 kept / 1.215 lost)")
        print("      still rules out differential noise either way.")
        print()
        print("      STATED CORRECTLY THIS IS A STRONGER POSITION, not a weaker one:")
        print("      events surviving adjustment for part of the disease phenotype is a")
        print("      more demanding test than events surviving a nuisance correction.")
        print("      Required edits:")
        print("        * Methods/Results: reframe the composition analysis as a")
        print("          conservative sensitivity test against a mediator, not a")
        print("          confounding check.")
        print("        * Say explicitly that the headline 926 comes from the BASE model,")
        print("          which does not adjust for composition -- now a deliberate")
        print("          choice rather than an unexplained default.")
        print("        * Limitations: one sentence that composition change is part of")
        print("          the failing phenotype, so the adjusted counts are lower bounds.")
        print("        * Drop 'threshold attrition rather than confounding'. It is")
        print("          attrition plus conservative over-adjustment, not noise.")

    if partB:
        pb = df[(df["part"] == "B") & (df["panel"] != "EXCLUDED by hand")]
        ex = df[(df["part"] == "B") & (df["panel"] == "EXCLUDED by hand")]
        pb_b = pb["beta_sd_technical"].abs().dropna()
        ex_b = ex["beta_sd_technical"].abs().dropna()
        nsig = int((pb["fdr_technical"] < 0.05).sum())
        print(f"\n  PART B: {nsig} of {len(pb_b)} markers are themselves differential "
              f"at FDR < 0.05")
        if len(pb_b) and len(ex_b):
            if pb_b.max() < ex_b.min():
                print(f"  >>> The hand-drawn exclusion rule is COHERENT: every excluded")
                print(f"      gene ({ex_b.min():.3f} SD at minimum) exceeds every kept one")
                print(f"      (max {pb_b.max():.3f} SD). The line lands where a measured")
                print(f"      one would, and the Methods can say so.")
            else:
                print(f"  >>> The exclusion rule is NOT clean: a KEPT marker reaches")
                print(f"      {pb_b.max():.3f} SD against an excluded minimum of "
                      f"{ex_b.min():.3f}.")
                print(f"      The Methods cannot describe the exclusions as a clear line.")

    df.to_csv(args.out, sep="\t", index=False)
    print(f"\n  wrote {args.out}  ({len(df)} rows)\n")


if __name__ == "__main__":
    main()
