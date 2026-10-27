#!/usr/bin/env python3
"""
03_rbm20_psi_correlation.py

The mechanistic keystone: does RBM20 expression track PSI of its target exons
in *your* data, beyond disease label and beyond cell composition?

What it does, per event (rMATS, one event type per run):
  1. per-sample PSI from IncLevel1/IncLevel2, low-coverage masked, deduped by
     event coordinate, PSI-variance filtered (kills the constant-PSI nan class
     that produced the fake LVEF hits).
  2. base disease model      PSI ~ disease + covs            -> adjusted dPSI
  3. within-DCM model        PSI ~ RBM20 + covs              -> is PSI graded by RBM20?
  4. all-sample model        PSI ~ RBM20 + disease + covs    -> RBM20 beyond disease
  5. every model twice: without and with cell-composition adjustment (CM/FB score)
  6. target vs matched-background enrichment (matched on mean PSI and PSI sd)
  7. mechanism direction test: corr(adjusted dPSI, beta_RBM20) across targets
     -- if RBM20 loss drives the DCM shift this must be NEGATIVE
  8. disease-effect attenuation when RBM20 enters the model (mediation-flavoured,
     reported descriptively, NOT as a causal mediation test)
  9. negative-control RBPs run through the identical pipeline
 10. optional permutation null (shuffle RBM20 within disease group)

Honest-null design: every headline number has a matched-background comparison and
a control-gene comparison. If RBM20 targets are not more RBM20-correlated than
background exons, the script says so.

Example:
  python 03_rbm20_psi_correlation.py \
    --rmats /data/project/chureeratp/RNA/HEART/rmats_magnet/SE.MATS.JC.txt \
    --event-type SE \
    --b1 /data/project/chureeratp/RNA/HEART/dcm_bams.txt \
    --b2 /data/project/chureeratp/RNA/HEART/nf_bams.txt \
    --expr magnet_expr.per_sample.tsv \
    --covariates /data/project/chureeratp/RNA/HEART/magnet_covariates.tsv \
    --target-genes TTN,CAMK2D,LDB3,LMO7,SORBS1,PDLIM3,DTNA,MLIP,TNNT2 \
    --out-prefix rbm20_SE
"""
import argparse
import os
import re
import sys

import numpy as np
import pandas as pd
from scipy import stats
from scipy.spatial import cKDTree

COORD_COLS = {
    "SE":   ["chr", "strand", "exonStart_0base", "exonEnd",
             "upstreamES", "upstreamEE", "downstreamES", "downstreamEE"],
    "MXE":  ["chr", "strand", "1stExonStart_0base", "1stExonEnd",
             "2ndExonStart_0base", "2ndExonEnd",
             "upstreamES", "upstreamEE", "downstreamES", "downstreamEE"],
    "A3SS": ["chr", "strand", "longExonStart_0base", "longExonEnd",
             "shortES", "shortEE", "flankingES", "flankingEE"],
    "A5SS": ["chr", "strand", "longExonStart_0base", "longExonEnd",
             "shortES", "shortEE", "flankingES", "flankingEE"],
    "RI":   ["chr", "strand", "riExonStart_0base", "riExonEnd",
             "upstreamES", "upstreamEE", "downstreamES", "downstreamEE"],
}


# ---------------------------------------------------------------- utilities
def bh_fdr(p):
    p = np.asarray(p, float)
    ok = ~np.isnan(p)
    q = np.full(p.shape, np.nan)
    if ok.sum() == 0:
        return q
    pv = p[ok]
    n = pv.size
    order = np.argsort(pv)
    ranked = pv[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(ranked, 0, 1)
    q[ok] = out
    return q


def sample_id_from_bam(path, mode):
    path = path.strip()
    if not path:
        return None
    if mode == "parentdir":
        return os.path.basename(os.path.dirname(path))
    base = os.path.basename(path)
    for suf in ["Aligned.sortedByCoord.out.bam", ".bam"]:
        if base.endswith(suf):
            base = base[: -len(suf)]
            break
    return base.rstrip("_.")


def read_bam_list(path, mode):
    ids = []
    with open(path) as fh:
        for line in fh:
            for tok in line.strip().split(","):
                sid = sample_id_from_bam(tok, mode)
                if sid:
                    ids.append(sid)
    return ids


def split_num(cell, n_expected):
    if pd.isna(cell):
        return np.full(n_expected, np.nan)
    vals = str(cell).split(",")
    out = np.full(n_expected, np.nan)
    for i, v in enumerate(vals[:n_expected]):
        v = v.strip()
        if v not in ("", "NA", "nan", "None"):
            try:
                out[i] = float(v)
            except ValueError:
                pass
    return out


def ols(y, X, cols):
    """OLS with t-tests. Returns dict col->(beta, t, p) plus n, df. None if unusable."""
    n, p = X.shape
    df = n - p
    if df < 5:
        return None
    # A single NaN anywhere makes matrix_rank raise LinAlgError("SVD did not
    # converge") and kills the whole run. Never let one through.
    if not (np.all(np.isfinite(X)) and np.all(np.isfinite(y))):
        return None
    try:
        if np.linalg.matrix_rank(X) < p:
            return None
    except np.linalg.LinAlgError:
        return None
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    s2 = float(resid @ resid) / df
    if not np.isfinite(s2) or s2 <= 0:
        return None
    xtx_inv = np.linalg.pinv(X.T @ X)
    se = np.sqrt(np.maximum(np.diag(xtx_inv) * s2, 0))
    with np.errstate(divide="ignore", invalid="ignore"):
        t = np.where(se > 0, beta / se, np.nan)
    pv = 2 * stats.t.sf(np.abs(t), df)
    return {"n": n, "df": df,
            "res": {c: (beta[i], t[i], pv[i]) for i, c in enumerate(cols)}}


def build_design(meta, cols_numeric, cols_cat):
    """Full design frame (intercept added later); dummies drop-first."""
    parts = []
    names = []
    for c in cols_numeric:
        v = pd.to_numeric(meta[c], errors="coerce")
        parts.append(v.values.astype(float))
        names.append(c)
    for c in cols_cat:
        d = pd.get_dummies(meta[c].astype(str), prefix=c, drop_first=True, dtype=float)
        for cn in d.columns:
            parts.append(d[cn].values.astype(float))
            names.append(cn)
    if not parts:
        return np.zeros((len(meta), 0)), []
    return np.column_stack(parts), names


def fit_event(psi, Xfull, names, key_terms):
    """Subset to usable samples, drop constant columns, fit."""
    ok = np.isfinite(psi) & np.isfinite(Xfull).all(axis=1)
    if ok.sum() < 20:
        return None
    y = psi[ok]
    X = Xfull[ok]
    keep = [i for i in range(X.shape[1]) if np.nanstd(X[:, i]) > 1e-9]
    X = X[:, keep]
    nm = [names[i] for i in keep]
    if not all(k in nm for k in key_terms):
        return None
    X = np.column_stack([np.ones(X.shape[0]), X])
    nm = ["Intercept"] + nm
    return ols(y, X, nm)


# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rmats", default=None, help="e.g. SE.MATS.JC.txt")
    ap.add_argument("--from-events", default=None,
                    help="reuse a completed <prefix>.events.tsv and only "
                         "regenerate the summary report (skips ~30 min of "
                         "refitting; permutation test unavailable in this mode)")
    ap.add_argument("--event-type", required=True, choices=list(COORD_COLS))
    ap.add_argument("--b1", required=True, help="bam list for group1 (DCM)")
    ap.add_argument("--b2", required=True, help="bam list for group2 (NF)")
    ap.add_argument("--id-from", choices=["parentdir", "basename"], default="parentdir")
    ap.add_argument("--expr", required=True, help="magnet_expr.per_sample.tsv from step 02")
    ap.add_argument("--covariates", required=True, help="sample,disease,sex,age,RIN[,pool]")
    ap.add_argument("--cov-names", default="sample,disease,sex,age,RIN,pool",
                    help="column names used when the covariate file is headerless "
                         "(magnet_covariates.tsv is)")
    ap.add_argument("--numeric-cov", default="age,RIN")
    ap.add_argument("--cat-cov", default="sex")
    ap.add_argument("--cell-cov", default="CM_score,FB_score")
    ap.add_argument("--target-genes",
                    default="TTN,CAMK2D,LDB3,LMO7,SORBS1,PDLIM3,DTNA,MLIP,TNNT2")
    ap.add_argument("--control-genes", default="MBNL1,CELF1,SRSF1,PTBP1,QKI,GAPDH")
    ap.add_argument("--min-reads", type=int, default=10,
                    help="mask a sample's PSI if IJC+SJC < this")
    ap.add_argument("--min-frac-samples", type=float, default=0.7)
    ap.add_argument("--min-psi-sd", type=float, default=0.02,
                    help="THE variance filter; below this you get nan/artefact stats")
    ap.add_argument("--n-match", type=int, default=10, help="background exons per target")
    ap.add_argument("--n-perm", type=int, default=0,
                    help="permutations of RBM20 within disease group (200 is plenty)")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out-prefix", default="rbm20_psi")
    args = ap.parse_args()

    if not args.rmats and not args.from_events:
        ap.error('need --rmats (full run) or --from-events (report only)')
    rng = np.random.default_rng(args.seed)
    et = args.event_type
    log = lambda *a: print(*a, file=sys.stderr)

    # ---- samples -------------------------------------------------------
    g1 = read_bam_list(args.b1, args.id_from)
    g2 = read_bam_list(args.b2, args.id_from)
    samples = g1 + g2
    disease = np.array([1.0] * len(g1) + [0.0] * len(g2))
    log(f"[info] group1(DCM)={len(g1)} group2(NF)={len(g2)}")
    if len(set(samples)) != len(samples):
        log("[warn] duplicate sample ids across bam lists -- check --id-from")

    # magnet_covariates.tsv is HEADERLESS: sample disease sex age RIN pool.
    # Auto-detect so either flavour works.
    first = open(args.covariates).readline().rstrip("\n").split("\t")
    headerless = bool(first) and (
        re.match(r"^(SRR|ERR|SAMEA|SAMN)\d+", first[0].strip()) is not None)
    if headerless:
        names = args.cov_names.split(",")
        if len(first) != len(names):
            sys.exit(f"[fatal] covariates file has {len(first)} columns but "
                     f"--cov-names lists {len(names)}: {names}")
        cov = pd.read_csv(args.covariates, sep="\t", header=None, names=names)
        log(f"[info] covariates read as HEADERLESS with columns {names}")
    else:
        cov = pd.read_csv(args.covariates, sep=None, engine="python")
        cov.columns = [c.strip() for c in cov.columns]
        log(f"[info] covariates read WITH header: {list(cov.columns)}")
    scol = cov.columns[0]
    cov = cov.set_index(scol)
    for c in ["age", "RIN", "LVEF"]:
        if c in cov.columns:
            cov[c] = pd.to_numeric(cov[c], errors="coerce")
    expr = pd.read_csv(args.expr, sep="\t").set_index("sample")

    meta = pd.DataFrame(index=pd.Index(samples, name="sample"))
    meta["disease"] = disease

    # Guard: --b1/--b2 swapped flips every sign in this analysis. Check the
    # bam-list grouping against the covariate file's own disease labels.
    if "disease" in cov.columns:
        lab = cov["disease"].astype(str).reindex(meta.index)
        n_unknown = int(lab.isna().sum())
        if n_unknown:
            log(f"[warn] {n_unknown} bam-list samples absent from the covariate file")
        tab = pd.crosstab(lab.fillna("MISSING"), meta["disease"])
        log("[info] bam-list group (1=group1/--b1) vs covariate label:\n" + tab.to_string())
        if 1.0 in tab.columns and len(tab.index):
            top1 = tab[1.0].idxmax()
            purity = tab.loc[top1, 1.0] / max(tab[1.0].sum(), 1)
            log(f"[info] --b1 is predominantly '{top1}' ({purity:.0%} pure)")
            if purity < 0.9:
                log("[warn] --b1 is NOT a clean single group -- check your bam lists")
            if str(top1).upper() in ("NF", "CONTROL", "NORMAL", "DONOR"):
                log("[warn] --b1 looks like the CONTROL group. IncLevel1 is group1, so "
                    "every disease beta and every direction test below will be sign-"
                    "flipped relative to your phase-2 tables. Swap --b1/--b2.")

    for c in cov.columns:
        if c not in meta.columns:
            meta[c] = cov[c].reindex(meta.index)
    for c in expr.columns:
        meta[c] = expr[c].reindex(meta.index)

    n_missing_expr = meta["RBM20"].isna().sum()
    if n_missing_expr:
        log(f"[warn] {n_missing_expr} samples lack expression data; they drop out of "
            "every RBM20 model (disease models keep them).")

    # z-score RBM20 and controls so betas are 'PSI change per SD of expression'
    gene_cols = ["RBM20"] + [g for g in args.control_genes.split(",")
                             if g and g in meta.columns]
    for g in gene_cols:
        v = pd.to_numeric(meta[g], errors="coerce")
        meta[g] = (v - v.mean()) / v.std()

    num_cov = [c for c in args.numeric_cov.split(",") if c and c in meta.columns]
    cat_cov = [c for c in args.cat_cov.split(",") if c and c in meta.columns]
    cell_cov = [c for c in args.cell_cov.split(",") if c and c in meta.columns]
    log(f"[info] covariates numeric={num_cov} categorical={cat_cov} cell={cell_cov}")

    # ---- resume from an existing events table ---------------------------
    # The per-event fitting is the expensive part (~30 min). If it already
    # completed and only the summary failed, reuse it.
    if args.from_events:
        log(f"[info] RESUMING from {args.from_events} -- no refitting")
        res = pd.read_csv(args.from_events, sep="\t")
        log(f"[info] loaded {len(res)} events")
        psi = None
        designs = {}
        dcm_mask = None
        report(args, res, meta, num_cov, cat_cov, cell_cov, gene_cols,
               et, psi, designs, dcm_mask, rng, log)
        return

    # ---- rMATS ---------------------------------------------------------
    log(f"[info] reading {args.rmats}")
    df = pd.read_csv(args.rmats, sep="\t")
    need = COORD_COLS[et] + ["geneSymbol", "IncLevel1", "IncLevel2",
                             "IJC_SAMPLE_1", "SJC_SAMPLE_1", "IJC_SAMPLE_2", "SJC_SAMPLE_2"]
    miss = [c for c in need if c not in df.columns]
    if miss:
        sys.exit(f"[fatal] missing rMATS columns {miss} -- wrong --event-type?")
    df["key"] = df[COORD_COLS[et]].astype(str).agg(":".join, axis=1)

    # dedup multi-transcript duplicates: keep the row with most total reads
    def tot_reads(r):
        v = 0.0
        for c in ["IJC_SAMPLE_1", "SJC_SAMPLE_1", "IJC_SAMPLE_2", "SJC_SAMPLE_2"]:
            v += np.nansum(split_num(r[c], 10_000))
        return v
    n0 = len(df)
    df["_tot"] = df.apply(tot_reads, axis=1)
    df = df.sort_values("_tot", ascending=False).drop_duplicates("key", keep="first")
    log(f"[info] {n0} rows -> {len(df)} unique {et} events after coordinate dedup")

    n1, n2 = len(g1), len(g2)
    psi = np.full((len(df), n1 + n2), np.nan)
    cover = np.full((len(df), n1 + n2), np.nan)
    for i, (_, r) in enumerate(df.iterrows()):
        p1 = split_num(r["IncLevel1"], n1)
        p2 = split_num(r["IncLevel2"], n2)
        c1 = split_num(r["IJC_SAMPLE_1"], n1) + split_num(r["SJC_SAMPLE_1"], n1)
        c2 = split_num(r["IJC_SAMPLE_2"], n2) + split_num(r["SJC_SAMPLE_2"], n2)
        psi[i] = np.concatenate([p1, p2])
        cover[i] = np.concatenate([c1, c2])

    psi[~(cover >= args.min_reads)] = np.nan

    frac_ok = np.isfinite(psi).mean(axis=1)
    sd = np.nanstd(psi, axis=1)
    with np.errstate(invalid="ignore"):
        keep = (frac_ok >= args.min_frac_samples) & (sd >= args.min_psi_sd)
    log(f"[info] events passing coverage>={args.min_frac_samples:.0%} of samples "
        f"and PSI sd>={args.min_psi_sd}: {keep.sum()} / {len(df)}")
    log(f"[info]   (dropped for low coverage: {(frac_ok < args.min_frac_samples).sum()}, "
        f"for constant PSI: {(sd < args.min_psi_sd).sum()})")
    df = df[keep].reset_index(drop=True)
    psi = psi[keep]
    if len(df) == 0:
        sys.exit("[fatal] no events survive filtering")

    # ---- designs -------------------------------------------------------
    Xbase, nbase = build_design(meta, num_cov, cat_cov)
    Xcell, ncell = build_design(meta, num_cov + cell_cov, cat_cov)
    dis = meta["disease"].values.astype(float)[:, None]

    designs = {}
    for tag, (Xc, nc) in {"nocell": (Xbase, nbase), "cell": (Xcell, ncell)}.items():
        designs[f"dis_{tag}"] = (np.hstack([dis, Xc]), ["disease"] + nc, ["disease"])
        for g in gene_cols:
            gv = meta[g].values.astype(float)[:, None]
            designs[f"ALL_{g}_{tag}"] = (np.hstack([gv, dis, Xc]),
                                         [g, "disease"] + nc, [g, "disease"])
            designs[f"DCM_{g}_{tag}"] = (np.hstack([gv, Xc]), [g] + nc, [g])

    dcm_mask = dis[:, 0] == 1

    # ---- per-event fits -------------------------------------------------
    out = {"key": df["key"].values, "gene": df["geneSymbol"].values,
           "event_type": et, "n_ok": np.isfinite(psi).sum(axis=1),
           "mean_psi": np.nanmean(psi, axis=1), "sd_psi": np.nanstd(psi, axis=1)}
    stat = {}

    def run(name, rows_mask=None):
        X, nm, terms = designs[name]
        b = {t: np.full(len(df), np.nan) for t in terms}
        tt = {t: np.full(len(df), np.nan) for t in terms}
        pp = {t: np.full(len(df), np.nan) for t in terms}
        nn = np.full(len(df), np.nan)
        Xu = X if rows_mask is None else X[rows_mask]
        for i in range(len(df)):
            y = psi[i] if rows_mask is None else psi[i][rows_mask]
            fit = fit_event(y, Xu, nm, terms)
            if fit is None:
                continue
            nn[i] = fit["n"]
            for t in terms:
                b[t][i], tt[t][i], pp[t][i] = fit["res"][t]
        for t in terms:
            stat[f"{name}.{t}.beta"] = b[t]
            stat[f"{name}.{t}.t"] = tt[t]
            stat[f"{name}.{t}.p"] = pp[t]
            stat[f"{name}.{t}.fdr"] = bh_fdr(pp[t])
        stat[f"{name}.n"] = nn

    log("[info] fitting disease models")
    run("dis_nocell"); run("dis_cell")
    for tag in ["nocell", "cell"]:
        for g in gene_cols:
            log(f"[info] fitting {g} models ({tag})")
            run(f"ALL_{g}_{tag}")
            run(f"DCM_{g}_{tag}", rows_mask=dcm_mask)

    res = pd.DataFrame({**out, **stat})
    res_path = f"{args.out_prefix}.events.tsv"
    res.to_csv(res_path, sep="\t", index=False, float_format="%.6g")
    log(f"[write] {res_path}")

    report(args, res, meta, num_cov, cat_cov, cell_cov, gene_cols,
           et, psi, designs, dcm_mask, rng, log)


def report(args, res, meta, num_cov, cat_cov, cell_cov, gene_cols,
           et, psi, designs, dcm_mask, rng, log):
    """Everything after the per-event fitting. Split out so --from-events
    can reuse a completed events table without refitting."""
    # ---- target / background -------------------------------------------
    targets = set(g.strip() for g in args.target_genes.split(",") if g.strip())
    is_tgt = res["gene"].isin(targets).values
    log(f"[info] target events: {is_tgt.sum()} in genes {sorted(targets)}")
    if is_tgt.sum() < 5:
        log("[warn] very few target events -- enrichment tests will be weak")

    feats = np.column_stack([res["mean_psi"].values, res["sd_psi"].values])
    fz = (feats - np.nanmean(feats, axis=0)) / np.nanstd(feats, axis=0)
    bg_idx = np.where(~is_tgt & np.isfinite(fz).all(axis=1))[0]
    tg_idx = np.where(is_tgt & np.isfinite(fz).all(axis=1))[0]
    matched = np.array([], int)
    if len(bg_idx) and len(tg_idx):
        tree = cKDTree(fz[bg_idx])
        k = min(args.n_match, len(bg_idx))
        _, nb = tree.query(fz[tg_idx], k=k)
        matched = np.unique(bg_idx[np.atleast_2d(nb).ravel()])
    log(f"[info] matched background events: {len(matched)}")

    lines = []
    def say(s=""):
        lines.append(s)
        print(s)

    src = os.path.basename(args.from_events or args.rmats or "")
    say(f"=== RBM20 <-> PSI, {et}, {src} ===")
    say(f"events tested: {len(res)}   target: {is_tgt.sum()}   "
        f"matched background: {len(matched)}")
    say()

    # 0. Is RBM20 itself different between DCM and NF?
    say("--- 0. RBM20 expression, DCM vs NF (the premise) ---")
    for tag, cc in [("unadjusted", []), ("cov-adjusted", num_cov + cat_cov),
                    ("cov+cell-adjusted", num_cov + cat_cov + cell_cov)]:
        nc_ = [c for c in cc if c in num_cov + cell_cov]
        cc_ = [c for c in cc if c in cat_cov]
        # drop samples missing ANY covariate used here, not just RBM20 --
        # 6 MAGNet samples have RIN = NA
        sub = meta.dropna(subset=["RBM20", "disease"] + nc_ + cc_)
        if len(sub) < 20:
            say(f"  {tag:<20} skipped -- only {len(sub)} samples with complete data")
            continue
        Xc, nm = build_design(sub, nc_, cc_)
        X = np.column_stack([np.ones(len(sub)), sub["disease"].values.astype(float)] +
                            ([Xc] if Xc.shape[1] else []))
        nm = ["Intercept", "disease"] + nm
        f = ols(sub["RBM20"].values.astype(float), X, nm)
        if f:
            b, t, p = f["res"]["disease"]
            say(f"  {tag:<20} beta(DCM-NF) = {b:+.3f} SD   t={t:+.2f}  p={p:.3g}  n={f['n']}")
        else:
            say(f"  {tag:<20} model not estimable (n={len(sub)})")
    say("  (if RBM20 is NOT lower in DCM, say so plainly -- the program can still be")
    say("   RBM20-dependent via activity/phosphorylation, but you cannot claim")
    say("   'reduced RBM20 expression' as the driver.)")
    say()

    # 1. Target enrichment for RBM20 association
    for model in ["DCM", "ALL"]:
        for tag in ["nocell", "cell"]:
            col = f"{model}_RBM20_{tag}.RBM20.t"
            if col not in res:
                continue
            tv = np.abs(res[col].values)
            a = tv[tg_idx]; b_ = tv[matched] if len(matched) else tv[bg_idx]
            a = a[np.isfinite(a)]; b_ = b_[np.isfinite(b_)]
            if len(a) < 3 or len(b_) < 3:
                continue
            u, pu = stats.mannwhitneyu(a, b_, alternative="greater")
            nsig = int((res[f"{model}_RBM20_{tag}.RBM20.fdr"].values[tg_idx] < 0.05).sum())
            say(f"--- 1. {model} model, {tag}: |t| for RBM20 ---")
            say(f"  target median |t| = {np.median(a):.2f}   matched-bg median |t| = "
                f"{np.median(b_):.2f}   MWU p = {pu:.3g}")
            say(f"  target events with genome-wide FDR<0.05 for RBM20: {nsig}/{len(tg_idx)}")
    say()

    # 2. Mechanism direction
    say("--- 2. Direction test: does the disease shift go the way RBM20 predicts? ---")
    for tag in ["nocell", "cell"]:
        d = res[f"dis_{tag}.disease.beta"].values[tg_idx]
        for model in ["DCM", "ALL"]:
            b_ = res[f"{model}_RBM20_{tag}.RBM20.beta"].values[tg_idx]
            m = np.isfinite(d) & np.isfinite(b_)
            if m.sum() < 6:
                continue
            rho, pr = stats.spearmanr(d[m], b_[m])
            say(f"  {tag:<7} {model:<3}  rho(adjusted dPSI, beta_RBM20) = {rho:+.3f}  "
                f"p={pr:.3g}  n={m.sum()}")
    say("  NEGATIVE rho = exons that shift in DCM are exactly those RBM20 level")
    say("  moves the other way => consistent with loss of RBM20 function.")
    say("  Positive or null rho = the disease shift is NOT explained by RBM20 level.")
    say()

    # 3. Attenuation of disease effect
    say("--- 3. Disease effect with vs without RBM20 in the model (descriptive) ---")
    for tag in ["nocell", "cell"]:
        b0 = res[f"dis_{tag}.disease.beta"].values[tg_idx]
        b1 = res[f"ALL_RBM20_{tag}.disease.beta"].values[tg_idx]
        m = np.isfinite(b0) & np.isfinite(b1) & (np.abs(b0) > 0.02)
        if m.sum() < 5:
            continue
        att = 1 - (b1[m] / b0[m])
        say(f"  {tag:<7} median attenuation of |disease beta| = {np.median(att)*100:+.1f}% "
            f"over {m.sum()} target events")
    say("  This is NOT a mediation test (no causal identification, RBM20 and PSI come")
    say("  from the same library). Report as 'partially collinear', not 'mediated'.")
    say()

    # 4. Control genes
    say("--- 4. Specificity: control RBPs through the identical pipeline ---")
    say(f"  {'gene':<10}{'model':<6}{'adj':<8}{'median|t| tgt':>14}{'median|t| bg':>14}{'MWU p':>10}")
    for g in gene_cols:
        for model in ["DCM"]:
            for tag in ["nocell", "cell"]:
                col = f"{model}_{g}_{tag}.{g}.t"
                if col not in res:
                    continue
                tv = np.abs(res[col].values)
                a = tv[tg_idx]; b_ = tv[matched] if len(matched) else tv[bg_idx]
                a = a[np.isfinite(a)]; b_ = b_[np.isfinite(b_)]
                if len(a) < 3 or len(b_) < 3:
                    continue
                _, pu = stats.mannwhitneyu(a, b_, alternative="greater")
                say(f"  {g:<10}{model:<6}{tag:<8}{np.median(a):>14.2f}{np.median(b_):>14.2f}"
                    f"{pu:>10.3g}")
    say("  MBNL1/CELF1/QKI are real cardiac splicing regulators -- if they score as")
    say("  high as RBM20 that is a finding about co-regulation, not a failure, but it")
    say("  blocks any claim that the program is RBM20-SPECIFIC.")
    say()

    # 5. Permutation null
    if args.n_perm > 0 and psi is not None and designs:
        say(f"--- 5. Permutation null ({args.n_perm} shuffles of RBM20 within disease group) ---")
        sel = np.unique(np.concatenate([tg_idx, matched])) if len(matched) else tg_idx
        obs = np.nanmedian(np.abs(res["DCM_RBM20_cell.RBM20.t"].values[tg_idx]))
        X, nm, terms = designs["DCM_RBM20_cell"]
        gi = nm.index("RBM20")
        null = np.full(args.n_perm, np.nan)
        rows = np.where(dcm_mask)[0]
        for k in range(args.n_perm):
            Xp = X.copy()
            perm = rng.permutation(rows)
            Xp[rows, gi] = X[perm, gi]
            vals = []
            for i in tg_idx:
                f = fit_event(psi[i][dcm_mask], Xp[dcm_mask], nm, terms)
                if f:
                    vals.append(abs(f["res"]["RBM20"][1]))
            if vals:
                null[k] = np.median(vals)
        pe = (np.nansum(null >= obs) + 1) / (np.isfinite(null).sum() + 1)
        say(f"  observed median |t| = {obs:.2f}; null median = {np.nanmedian(null):.2f}; "
            f"empirical p = {pe:.4g}")
        say()

    # 6. Per-gene table
    say("--- 6. Per-gene summary (within-DCM, cell-adjusted) ---")
    col_t = "DCM_RBM20_cell.RBM20.t"
    col_q = "DCM_RBM20_cell.RBM20.fdr"
    if col_t in res:
        sub = res[is_tgt]
        say(f"  {'gene':<10}{'events':>7}{'median|t|':>11}{'FDR<0.05':>10}{'median beta':>13}")
        for g, gg in sub.groupby("gene"):
            tv = np.abs(gg[col_t].values)
            say(f"  {g:<10}{len(gg):>7}{np.nanmedian(tv):>11.2f}"
                f"{int((gg[col_q] < 0.05).sum()):>10}"
                f"{np.nanmedian(gg['DCM_RBM20_cell.RBM20.beta']):>13.3f}")

    with open(f"{args.out_prefix}.report.txt", "w") as fh:
        fh.write("\n".join(lines) + "\n")
    log(f"[write] {args.out_prefix}.report.txt")


if __name__ == "__main__":
    main()
