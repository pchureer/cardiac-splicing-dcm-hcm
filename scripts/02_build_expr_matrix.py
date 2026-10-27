#!/usr/bin/env python3
"""
02_build_expr_matrix.py

Merge per-sample gene counts (STAR ReadsPerGene.out.tab OR featureCounts) into a
log2-CPM matrix, and emit the per-sample covariate file the correlation analysis
needs: RBM20 level, control-gene levels, and cell-composition scores.

Why cell scores: DCM hearts are more fibrotic, so cardiomyocyte fraction drops.
RBM20 is cardiomyocyte-restricted and its target exons are cardiomyocyte exons,
so ANY RBM20-PSI correlation is confounded by composition unless you adjust for
it. We therefore compute a CM score and an FB score and carry both forward; the
analysis script reports every association with and without that adjustment.

Usage (STAR counts already exist -- preferred):
  python 02_build_expr_matrix.py \
    --heart-dir /data/project/chureeratp/RNA/HEART \
    --ids /data/project/chureeratp/RNA/HEART/magnet_all.ids \
    --source star \
    --gtf /data/project/chureeratp/hg38/ref/gencode.v49.annotation.gtf \
    --out-prefix magnet_expr

Usage (featureCounts fallback):
  ... --source featurecounts --fc-dir /data/project/chureeratp/RNA/HEART/fc_counts
"""
import argparse
import glob
import gzip
import os
import sys

import numpy as np
import pandas as pd

# Cardiomyocyte panel: structural/contractile genes that are CM-restricted.
# Deliberately excludes NPPA/NPPB/MYH7 (strongly disease-regulated, would make
# the score a disease proxy rather than a composition proxy).
CM_MARKERS = ["TNNI3", "MYL2", "MYL3", "ACTN2", "TCAP", "CSRP3", "MYBPC3", "RYR2", "PLN"]
FB_MARKERS = ["COL1A1", "COL1A2", "COL3A1", "POSTN", "FN1", "DCN", "LUM"]
DEFAULT_CONTROL_GENES = ["MBNL1", "MBNL2", "CELF1", "SRSF1", "PTBP1", "QKI", "HNRNPA1", "GAPDH"]


def open_maybe_gz(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def gene_symbol_map(gtf):
    """gene_id -> gene_name from the GTF gene lines. Keeps versioned and
    unversioned ids so it works with either counts flavour."""
    m = {}
    with open_maybe_gz(gtf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            attr = f[8]
            try:
                gid = attr.split('gene_id "')[1].split('"')[0]
                gname = attr.split('gene_name "')[1].split('"')[0]
            except IndexError:
                continue
            m[gid] = gname
            m[gid.split(".")[0]] = gname
    return m


def read_star(path):
    """STAR ReadsPerGene.out.tab -> DataFrame index=gene_id, cols unstr/fwd/rev."""
    df = pd.read_csv(path, sep="\t", header=None,
                     names=["gene", "unstr", "fwd", "rev"])
    df = df[~df["gene"].str.startswith("N_")]
    return df.set_index("gene")


def read_featurecounts(path):
    df = pd.read_csv(path, sep="\t", comment="#")
    counts = df.iloc[:, -1]
    counts.index = df["Geneid"]
    return counts


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--heart-dir", required=True)
    p.add_argument("--ids", required=True, help="one sample/run accession per line")
    p.add_argument("--source", choices=["star", "featurecounts"], default="star")
    p.add_argument("--fc-dir", default=None)
    p.add_argument("--gtf", required=True)
    p.add_argument("--out-prefix", default="magnet_expr")
    p.add_argument("--min-cpm", type=float, default=1.0,
                   help="gene kept if CPM>=this in >=--min-frac of samples")
    p.add_argument("--min-frac", type=float, default=0.2)
    p.add_argument("--control-genes", default=",".join(DEFAULT_CONTROL_GENES))
    args = p.parse_args()

    samples = [s.strip() for s in open(args.ids) if s.strip()]
    print(f"[info] {len(samples)} samples requested", file=sys.stderr)

    cols, kept, missing = {}, [], []
    strand_tot = np.zeros(3)

    for s in samples:
        if args.source == "star":
            hits = sorted(glob.glob(os.path.join(args.heart_dir, s, "*ReadsPerGene.out.tab")))
            if not hits:
                missing.append(s)
                continue
            df = read_star(hits[0])
            strand_tot += df[["unstr", "fwd", "rev"]].sum().values
            cols[s] = df  # decide the column after seeing all samples
        else:
            hits = sorted(glob.glob(os.path.join(args.fc_dir, f"{s}.counts.txt")))
            if not hits:
                missing.append(s)
                continue
            cols[s] = read_featurecounts(hits[0])
        kept.append(s)

    if missing:
        print(f"[warn] {len(missing)} samples had no counts file, e.g. {missing[:5]}",
              file=sys.stderr)
    if not kept:
        sys.exit("[fatal] no counts files found -- check --source/--heart-dir/--fc-dir")

    if args.source == "star":
        # Strandedness: in an UNstranded library fwd ~= rev ~= unstr/2.
        # In a stranded library one of fwd/rev carries nearly everything.
        u, f_, r_ = strand_tot
        ratio = max(f_, r_) / max(min(f_, r_), 1.0)
        if ratio > 2.0:
            name = "fwd" if f_ > r_ else "rev"
            print(f"[info] library looks STRANDED (fwd/rev ratio {ratio:.1f}) -> using '{name}'",
                  file=sys.stderr)
        else:
            name = "unstr"
            print(f"[info] library looks UNSTRANDED (fwd/rev ratio {ratio:.2f}) -> using 'unstr'",
                  file=sys.stderr)
        print(f"[info] totals unstr/fwd/rev = {(strand_tot/max(u,1)).round(3)}", file=sys.stderr)
        cols = {s: df[name] for s, df in cols.items()}

    counts = pd.DataFrame(cols)
    counts = counts.fillna(0).astype(float)
    print(f"[info] count matrix {counts.shape[0]} genes x {counts.shape[1]} samples",
          file=sys.stderr)

    libsize = counts.sum(axis=0)
    print(f"[info] library size median {libsize.median():,.0f}, "
          f"min {libsize.min():,.0f}, max {libsize.max():,.0f}", file=sys.stderr)
    if libsize.min() < 1e6:
        low = libsize[libsize < 1e6]
        print(f"[warn] {len(low)} sample(s) under 1M assigned reads: {list(low.index)[:5]}",
              file=sys.stderr)

    cpm = counts.divide(libsize, axis=1) * 1e6
    expressed = (cpm >= args.min_cpm).mean(axis=1) >= args.min_frac
    cpm = cpm[expressed]
    print(f"[info] {expressed.sum()} genes pass expression filter", file=sys.stderr)
    log2cpm = np.log2(cpm + 1.0)

    sym = gene_symbol_map(args.gtf)
    symbols = pd.Index([sym.get(g, g) for g in log2cpm.index], name="gene")
    log2cpm.index = symbols
    # collapse duplicate symbols by max mean expression
    log2cpm = log2cpm.groupby(level=0).max()

    out_mat = f"{args.out_prefix}.log2cpm.tsv.gz"
    log2cpm.to_csv(out_mat, sep="\t", compression="gzip", float_format="%.4f")
    print(f"[write] {out_mat}", file=sys.stderr)

    def zscore_panel(genes, label):
        present = [g for g in genes if g in log2cpm.index]
        absent = [g for g in genes if g not in log2cpm.index]
        if absent:
            print(f"[warn] {label}: markers absent/filtered: {absent}", file=sys.stderr)
        if not present:
            return pd.Series(np.nan, index=log2cpm.columns)
        sub = log2cpm.loc[present]
        z = sub.sub(sub.mean(axis=1), axis=0).div(sub.std(axis=1).replace(0, np.nan), axis=0)
        return z.mean(axis=0)

    per = pd.DataFrame(index=log2cpm.columns)
    per.index.name = "sample"
    if "RBM20" not in log2cpm.index:
        sys.exit("[fatal] RBM20 absent from the filtered matrix -- lower --min-cpm "
                 "or check the symbol map. RBM20 is lowly expressed; this matters.")
    per["RBM20"] = log2cpm.loc["RBM20"]
    for g in [x for x in args.control_genes.split(",") if x]:
        if g in log2cpm.index:
            per[g] = log2cpm.loc[g]
        else:
            print(f"[warn] control gene {g} not in matrix", file=sys.stderr)
    per["CM_score"] = zscore_panel(CM_MARKERS, "CM")
    per["FB_score"] = zscore_panel(FB_MARKERS, "FB")
    per["libsize"] = libsize.reindex(per.index)
    per["log10_libsize"] = np.log10(per["libsize"])

    out_per = f"{args.out_prefix}.per_sample.tsv"
    per.reset_index().to_csv(out_per, sep="\t", index=False, float_format="%.5f")
    print(f"[write] {out_per}", file=sys.stderr)

    print("\n=== quick diagnostics ===", file=sys.stderr)
    print(f"RBM20 log2CPM: mean {per['RBM20'].mean():.2f}, sd {per['RBM20'].std():.2f}, "
          f"range {per['RBM20'].min():.2f}-{per['RBM20'].max():.2f}", file=sys.stderr)
    r = np.corrcoef(per["RBM20"], per["CM_score"])[0, 1]
    print(f"cor(RBM20, CM_score) = {r:.3f}  <-- if this is high (>0.6) the RBM20 "
          "signal is largely cell composition; the cell-adjusted models are the "
          "ones you report.", file=sys.stderr)
    r2 = np.corrcoef(per["RBM20"], per["FB_score"])[0, 1]
    print(f"cor(RBM20, FB_score) = {r2:.3f}", file=sys.stderr)


if __name__ == "__main__":
    main()
