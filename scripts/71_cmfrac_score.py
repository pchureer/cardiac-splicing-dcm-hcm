#!/usr/bin/env python3
"""
71_cmfrac_score.py -- a cardiomyocyte score restricted to markers that are NOT
upregulated in disease, plus the retention comparison that uses it.

WHY
---
`70_marker_panel_audit.py` found that the cardiomyocyte composition score is a MIXTURE.
Four of its nine markers are strongly UP in DCM and five are flat or down:

    up in DCM     MYL2 +0.857   CSRP3 +0.770   TNNI3 +0.685   ACTN2 +0.668
    flat          MYBPC3 +0.135   TCAP -0.057   RYR2 -0.067   PLN -0.148
    down          MYL3 -0.186

Myocyte LOSS predicts DOWN. So the four up-markers are tracking per-cell sarcomeric
upregulation, not cell fraction, and the score is part fraction and part disease.
`corr(CM_score, fraction-only score) = +0.699` -- r^2 = 0.49, so the two are NOT
interchangeable as covariates and the question cannot be settled by argument.

MYL3 is KEPT in the fraction panel even though it is significant, because it is
significant in the direction myocyte loss requires. Selecting on |effect| rather than on
direction would have thrown away the one marker behaving as the biology demands.

**The fibroblast score is deliberately left alone.** Fibrosis IS the composition change,
so FB markers being up in DCM is the FB score working. Restricting it would remove the
signal it exists to carry.

WHAT THIS ANSWERS
-----------------
A reviewer's strongest version of the objection: "your composition-adjusted robustness
check is circular, because the covariate is partly a disease score." The answer is a
number -- how many of the 926 survive when the covariate is restricted to markers not
upregulated in disease -- not a paragraph.

THE GATE THAT MATTERS
---------------------
Before building anything, this script RECONSTRUCTS the existing nine-marker CM_score from
the gene matrix and checks it against the CM_score column that `02_build_expr_matrix.py`
already wrote. If the reconstruction does not match to r > 0.999, then this script's
z-scoring is not script 02's recipe, the new score is not comparable to the old one, and
it refuses. Without that check the comparison would be between two differently-built
scores and the difference would be uninterpretable.

USAGE
-----
    conda activate rnaseq_tools
    cd /data/project/chureeratp/RNA/HEART

    # 1. build the score file (seconds; gates on reproducing the stored CM_score)
    python 71_cmfrac_score.py build

    # 2. run phase 2 with it -- see 72_run_cmfrac.sbatch

    # 3. compare retention against the record
    python 71_cmfrac_score.py compare --adj SE_v2cmfrac.all_events.tsv
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

GENE_MATRIX = "magnet_expr.log2cpm.tsv.gz"
PER_SAMPLE = "magnet_expr.per_sample.tsv"
OUT_SCORE = "magnet_expr.cmfrac.tsv"

CM_FULL = ["TNNI3", "MYL2", "MYL3", "ACTN2", "TCAP", "CSRP3", "MYBPC3", "RYR2", "PLN"]
FB_FULL = ["COL1A1", "COL1A2", "COL3A1", "POSTN", "FN1", "DCN", "LUM"]
# flat or down in DCM -- i.e. not carrying the disease upregulation. MYL3 is included
# because its effect is NEGATIVE, the direction myocyte loss predicts.
CM_FRACTION = ["MYBPC3", "TCAP", "RYR2", "PLN", "MYL3"]
# measured 6 Oct by 70_marker_panel_audit.py, technical model, SD units -- printed so the
# panel choice is visible rather than asserted
MEASURED = {"MYL2": +0.857, "CSRP3": +0.770, "TNNI3": +0.685, "ACTN2": +0.668,
            "MYBPC3": +0.135, "TCAP": -0.057, "RYR2": -0.067, "PLN": -0.148,
            "MYL3": -0.186}

RECON_MIN_R = 0.999          # the gate: our z-score recipe must be script 02's
BASE_TABLE = os.environ.get("BASE_TABLE", "SE_v2filtered.all_events.tsv")
RECORD_BASE_SIG = int(os.environ.get("RECORD_BASE_SIG", "926"))
RECORD_COMPADJ_KEPT = int(os.environ.get("RECORD_COMPADJ_KEPT", "590"))    # the nine-marker composition result, exon level
NKEEP_SE = 3
FDR_CUT, DPSI_CUT = 0.05, 0.10


def zscore_panel(mat, genes, label):
    """Script 02's recipe: mean of per-gene z-scores across samples."""
    present = [g for g in genes if g in mat.index]
    absent = [g for g in genes if g not in mat.index]
    if absent:
        print(f"  [warn] {label}: absent from the matrix: {absent}")
    if not present:
        return None, f"{label}: no markers present"
    sub = mat.loc[present]
    z = sub.sub(sub.mean(axis=1), axis=0).div(sub.std(axis=1), axis=0)
    return z.mean(axis=0), None


def do_build(args):
    for p in (GENE_MATRIX, PER_SAMPLE):
        if not os.path.exists(p):
            sys.exit(f"[fatal] {p} not found")

    print(f"\n  reading {GENE_MATRIX}")
    G = pd.read_csv(GENE_MATRIX, sep="\t", index_col=0)
    print(f"  {G.shape[0]:,} genes x {G.shape[1]} samples")
    P = pd.read_csv(PER_SAMPLE, sep="\t")
    P = P.rename(columns={P.columns[0]: "sample"}).set_index("sample")

    # ---------------------------------------------------------------- the gate
    recon, err = zscore_panel(G, CM_FULL, "CM_full")
    if err:
        sys.exit(f"[fatal] {err}")
    if "CM_score" not in P.columns:
        sys.exit(f"[fatal] {PER_SAMPLE} has no CM_score column to check against")
    common = [s for s in recon.index if s in P.index]
    a = recon.loc[common].to_numpy(float)
    b = pd.to_numeric(P["CM_score"].loc[common], errors="coerce").to_numpy(float)
    ok = np.isfinite(a) & np.isfinite(b)
    r = float(np.corrcoef(a[ok], b[ok])[0, 1])
    mad = float(np.max(np.abs(a[ok] - b[ok])))
    print(f"\n  GATE -- reconstruct the stored nine-marker CM_score from the matrix")
    print(f"    r = {r:.6f}   max |difference| = {mad:.4f}   over {int(ok.sum())} samples")
    if r < RECON_MIN_R:
        print(f"\n  REFUSING: r < {RECON_MIN_R}. This script's z-scoring is not the")
        print(f"  recipe 02_build_expr_matrix.py used, so a score built here would not")
        print(f"  be comparable to the CM_score the recorded results were adjusted for.")
        print(f"  Read 02's zscore_panel() and match it before going further.\n")
        sys.exit(1)
    print(f"    [ok] the recipe matches -- the new score will be comparable\n")

    # ------------------------------------------------------- the restricted score
    print("  the cardiomyocyte panel, with each marker's measured disease effect:")
    for g in CM_FULL:
        tag = "KEPT in fraction panel" if g in CM_FRACTION else "dropped (up in DCM)"
        print(f"    {g:<8}{MEASURED.get(g, float('nan')):+.3f} SD   {tag}")
    print(f"\n  fibroblast panel left UNCHANGED -- fibrosis is the composition change,")
    print(f"  so FB markers being up in DCM is the score working, not failing.\n")

    frac, err = zscore_panel(G, CM_FRACTION, "CM_fraction")
    if err:
        sys.exit(f"[fatal] {err}")
    fbs, err = zscore_panel(G, FB_FULL, "FB_full")
    if err:
        sys.exit(f"[fatal] {err}")

    out = pd.DataFrame({"sample": frac.index,
                        "CM_score_frac": frac.values,
                        "CM_score_full": recon.reindex(frac.index).values,
                        "FB_score": fbs.reindex(frac.index).values})
    rr = float(np.corrcoef(out["CM_score_frac"], out["CM_score_full"])[0, 1])
    print(f"  corr(CM_score_frac, CM_score_full) = {rr:+.3f}   "
          f"(r^2 = {rr**2:.2f})")
    if rr > 0.95:
        print("  >>> above 0.95 the two are effectively the same covariate and the")
        print("      re-run cannot move anything. Report the correlation instead.")
    out.to_csv(OUT_SCORE, sep="\t", index=False)
    print(f"\n  wrote {OUT_SCORE}  ({len(out)} samples, "
          f"columns: {', '.join(out.columns[1:])})")
    print(f"\n  NEXT: sbatch 72_run_cmfrac.sbatch")
    print(f"        then: python 71_cmfrac_score.py compare "
          f"--adj SE_v2cmfrac.all_events.tsv\n")


def exon_sig(path):
    """Exon-level significant set, the paper's convention. Same as 68/69."""
    t = pd.read_csv(path, sep="\t", low_memory=False)
    for c in ("key", "FDR"):
        if c not in t.columns:
            sys.exit(f"[fatal] {path}: no '{c}' column")
    d = next((c for c in ("dPSI_adj", "dPSI") if c in t.columns), None)
    if d is None:
        sys.exit(f"[fatal] {path}: no dPSI column")
    t = t.copy()
    t["ekey"] = [":".join(str(k).split(":")[:NKEEP_SE]) for k in t["key"]]
    t["_f"] = pd.to_numeric(t["FDR"], errors="coerce")
    t["_d"] = pd.to_numeric(t[d], errors="coerce")
    t = t.dropna(subset=["_f"])
    rep = t.loc[t.groupby("ekey")["_f"].idxmin()]
    sig = rep[(rep["_f"] < FDR_CUT) & (rep["_d"].abs() >= DPSI_CUT)]
    return set(sig["ekey"]), len(rep)


def do_compare(args):
    if not os.path.exists(BASE_TABLE):
        sys.exit(f"[fatal] {BASE_TABLE} not found")
    if not os.path.exists(args.adj):
        sys.exit(f"[fatal] {args.adj} not found -- has 72_run_cmfrac.sbatch finished?")

    base, nbase = exon_sig(BASE_TABLE)
    adj, nadj = exon_sig(args.adj)
    print(f"\n  base   {BASE_TABLE}: {len(base)} significant of {nbase:,} exons")
    print(f"  new    {args.adj}: {len(adj)} significant of {nadj:,} exons")

    if len(base) != RECORD_BASE_SIG:
        print(f"\n  REFUSING: the base table gives {len(base)}, not the recorded "
              f"{RECORD_BASE_SIG}.")
        print("  Something upstream has changed; the comparison would not be the one")
        print("  the manuscript reports.\n")
        sys.exit(1)
    print(f"  [ok] base reproduces {RECORD_BASE_SIG}")

    kept = len(base & adj)
    pct = 100.0 * kept / len(base)
    rec_pct = 100.0 * RECORD_COMPADJ_KEPT / RECORD_BASE_SIG
    print("\n" + "=" * 70)
    print("  RETENTION under a fraction-only cardiomyocyte score")
    print("=" * 70)
    print(f"    nine-marker CM score (the record): {RECORD_COMPADJ_KEPT}/"
          f"{RECORD_BASE_SIG} = {rec_pct:.1f}%")
    print(f"    fraction-only CM score:            {kept}/{len(base)} = {pct:.1f}%")
    print(f"    difference: {kept - RECORD_COMPADJ_KEPT:+d} events, "
          f"{pct - rec_pct:+.1f} points")

    print()
    if abs(pct - rec_pct) < 3.0:
        print("  >>> RETENTION IS ESSENTIALLY UNCHANGED. The circularity objection is")
        print("      answered with a number: restricting the covariate to markers not")
        print("      upregulated in disease does not materially change how many events")
        print("      survive. Quote both figures in S1 and the Methods.")
    elif pct > rec_pct:
        print("  >>> RETENTION IS HIGHER with the fraction-only score, which is the")
        print("      predicted direction: the nine-marker score was removing real")
        print("      disease signal along with composition. The difference is the size")
        print("      of the over-adjustment. **Report the fraction-only figure as the")
        print("      primary composition sensitivity analysis** and the nine-marker one")
        print("      as the conservative bound.")
    else:
        print("  >>> RETENTION IS LOWER with the fraction-only score, which is NOT the")
        print("      predicted direction and needs explaining before it is published.")
        print("      Check that the restricted score has not simply become noisier")
        print("      (five markers instead of nine) -- a noisier covariate can absorb")
        print("      variance without being a better composition proxy. Compare the two")
        print("      scores' standard deviations before concluding anything.")

    pd.DataFrame([dict(
        base_table=BASE_TABLE, adj_table=args.adj,
        base_significant=len(base), adj_significant=len(adj),
        retained=kept, retained_pct=round(pct, 2),
        record_retained=RECORD_COMPADJ_KEPT, record_pct=round(rec_pct, 2),
        delta_events=kept - RECORD_COMPADJ_KEPT,
        delta_points=round(pct - rec_pct, 2))]).to_csv(
            "S1_cmfrac_comparison.tsv", sep="\t", index=False)
    print("\n  wrote S1_cmfrac_comparison.tsv\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="mode", required=True)
    sub.add_parser("build", help="build the fraction-only score file")
    c = sub.add_parser("compare", help="retention against the record")
    c.add_argument("--adj", default="SE_v2cmfrac.all_events.tsv")
    a = ap.parse_args()
    (do_build if a.mode == "build" else do_compare)(a)


if __name__ == "__main__":
    main()
