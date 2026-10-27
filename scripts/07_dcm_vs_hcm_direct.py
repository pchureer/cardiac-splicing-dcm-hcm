#!/usr/bin/env python3
"""
07_dcm_vs_hcm_direct.py -- compare DCM against HCM DIRECTLY, with no NF reference.

WHY THIS IS NEEDED
------------------
DCM-vs-NF and HCM-vs-NF are anchored on the SAME 166 donor hearts. They therefore
share every systematic difference between explanted failing myocardium and donor
myocardium: cardiomyocyte fraction, fibrosis, agonal state, ischaemic time,
medication, procurement and tissue handling. None of that is disease-specific and
none of it is captured by the shared-control sampling-error baseline in script 05
(which only models the shared NOISE in the NF mean, not shared SIGNAL).

That is why 98-99% of ALL significant events transfer between the two contrasts at
a slope of ~1. Two contrasts that share a reference will agree on everything the
reference differs from -- whether or not the two diseases differ from each other.

The only way to ask "do DCM and HCM differ?" is to contrast them directly and drop
NF entirely. That is what this does.

PSI is a per-sample quantity computed from that sample's own reads, so DCM PSI from
the discovery run and HCM PSI from the HCM run are directly comparable for the same
event, provided both runs used the same GTF and read length (they did). The script
VERIFIES this by checking the 166 NF samples, which appear in BOTH runs and must
give identical PSI. If that check fails, nothing else here is trustworthy.

Usage:
  python 07_dcm_vs_hcm_direct.py \
     --mats-dcm rmats_magnet/SE.MATS.JCEC.txt \
     --mats-hcm rmats_HCMvsNF/SE.MATS.JCEC.txt \
     --event-type SE \
     --b1-dcm dcm_bams.txt --b1-hcm hcm_bams.txt --b2-nf nf_bams.txt \
     --cov-dcm magnet_covariates.tsv --cov-hcm hcm_nf_covariates.tsv \
     --disc SE_v2filtered.all_events.tsv \
     --out SE_DCMvsHCM
"""
import argparse
import re
import sys

import numpy as np
import pandas as pd
from scipy import stats

FULL = {
    "SE":   ["exonStart_0base", "exonEnd", "upstreamES", "upstreamEE",
             "downstreamES", "downstreamEE"],
    "MXE":  ["1stExonStart_0base", "1stExonEnd", "2ndExonStart_0base", "2ndExonEnd",
             "upstreamES", "upstreamEE", "downstreamES", "downstreamEE"],
    "A3SS": ["longExonStart_0base", "longExonEnd", "shortES", "shortEE",
             "flankingES", "flankingEE"],
    "A5SS": ["longExonStart_0base", "longExonEnd", "shortES", "shortEE",
             "flankingES", "flankingEE"],
    "RI":   ["riExonStart_0base", "riExonEnd", "upstreamES", "upstreamEE",
             "downstreamES", "downstreamEE"],
}
PRIMARY = {"SE": ("exonStart_0base", "exonEnd"),
           "MXE": ("1stExonStart_0base", "1stExonEnd"),
           "A3SS": ("longExonStart_0base", "longExonEnd"),
           "A5SS": ("longExonStart_0base", "longExonEnd"),
           "RI": ("riExonStart_0base", "riExonEnd")}


def ids_from_bamlist(path):
    out = []
    for p in re.split(r"[,\n]", open(path).read().strip()):
        p = p.strip()
        if not p:
            continue
        m = re.search(r"(SRR\d+|ERR\d+|SAMEA\d+|SAMN\d+)", p)
        out.append(m.group(1) if m else p)
    return out


def parse_list(s):
    return [np.nan if x in ("NA", "", "nan") else float(x) for x in str(s).split(",")]


def read_cov(path):
    first = open(path).readline().rstrip("\n").split("\t")
    if re.match(r"^(SRR|ERR|SAMEA|SAMN)\d+", first[0].strip()):
        c = pd.read_csv(path, sep="\t", header=None,
                        names=["sample", "disease", "sex", "age", "RIN", "pool"][:len(first)])
    else:
        c = pd.read_csv(path, sep="\t")
        c = c.rename(columns={c.columns[0]: "sample"})
    for k in ("age", "RIN"):
        if k in c:
            c[k] = pd.to_numeric(c[k], errors="coerce")
    return c.set_index("sample")


def load_run(path, et, n1, n2, min_reads):
    """Return dict key -> (psi_g1 array, psi_g2 array, gene)."""
    df = pd.read_csv(path, sep="\t")
    df["key"] = df[["chr"] + FULL[et]].astype(str).agg(":".join, axis=1)
    cnt = ["IJC_SAMPLE_1", "SJC_SAMPLE_1", "IJC_SAMPLE_2", "SJC_SAMPLE_2"]
    have = all(c in df.columns for c in cnt)
    if have:
        df["_tot"] = sum(df[c].map(lambda s: np.nansum(parse_list(s))) for c in cnt)
    else:
        df["_tot"] = 0.0
    df = df.sort_values("_tot", ascending=False).drop_duplicates("key", keep="first")
    out = {}
    ps, pe = PRIMARY[et]
    for _, r in df.iterrows():
        p1 = np.array(parse_list(r["IncLevel1"]), float)
        p2 = np.array(parse_list(r["IncLevel2"]), float)
        if len(p1) != n1 or len(p2) != n2:
            continue
        if have and min_reads > 0:
            t1 = (np.array(parse_list(r["IJC_SAMPLE_1"]), float)
                  + np.array(parse_list(r["SJC_SAMPLE_1"]), float))
            t2 = (np.array(parse_list(r["IJC_SAMPLE_2"]), float)
                  + np.array(parse_list(r["SJC_SAMPLE_2"]), float))
            p1 = np.where(t1 >= min_reads, p1, np.nan)
            p2 = np.where(t2 >= min_reads, p2, np.nan)
        out[r["key"]] = (p1, p2, r["geneSymbol"], int(r[ps]), int(r[pe]), r["chr"])
    return out


def bh(p):
    p = np.asarray(p, float)
    ok = ~np.isnan(p)
    q = np.full(p.shape, np.nan)
    if ok.sum() == 0:
        return q
    pv = p[ok]
    n = pv.size
    o = np.argsort(pv)
    rk = pv[o] * n / (np.arange(n) + 1)
    rk = np.minimum.accumulate(rk[::-1])[::-1]
    out = np.empty(n)
    out[o] = np.clip(rk, 0, 1)
    q[ok] = out
    return q


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mats-dcm", required=True)
    ap.add_argument("--mats-hcm", required=True)
    ap.add_argument("--event-type", required=True, choices=list(FULL))
    ap.add_argument("--b1-dcm", required=True)
    ap.add_argument("--b1-hcm", required=True)
    ap.add_argument("--b2-nf", required=True)
    ap.add_argument("--cov-dcm", required=True)
    ap.add_argument("--cov-hcm", required=True)
    ap.add_argument("--disc", default=None,
                    help="phase2_v2 DCM-vs-NF table, for the equivalence comparison")
    ap.add_argument("--genes", default="TTN,CAMK2D,LDB3,LMO7,SORBS1,PDLIM3,DTNA,MLIP,TNNT2")
    ap.add_argument("--min-reads", type=int, default=10)
    ap.add_argument("--min-samples-hcm", type=int, default=15)
    ap.add_argument("--min-psi-sd", type=float, default=0.02)
    ap.add_argument("--out", default="dcm_vs_hcm")
    a = ap.parse_args()

    et = a.event_type
    log = lambda *s: print(*s, file=sys.stderr)
    lines = []
    def say(s=""):
        lines.append(s)
        print(s)

    dcm_ids = ids_from_bamlist(a.b1_dcm)
    hcm_ids = ids_from_bamlist(a.b1_hcm)
    nf_ids = ids_from_bamlist(a.b2_nf)
    log(f"[info] DCM {len(dcm_ids)}  HCM {len(hcm_ids)}  NF {len(nf_ids)}")
    ov = set(dcm_ids) & set(hcm_ids)
    if ov:
        sys.exit(f"[fatal] {len(ov)} samples appear in BOTH the DCM and HCM lists: "
                 f"{sorted(ov)[:5]}")

    log("[info] reading discovery run")
    runD = load_run(a.mats_dcm, et, len(dcm_ids), len(nf_ids), a.min_reads)
    log(f"[info]   {len(runD)} events")
    log("[info] reading HCM run")
    runH = load_run(a.mats_hcm, et, len(hcm_ids), len(nf_ids), a.min_reads)
    log(f"[info]   {len(runH)} events")

    shared = sorted(set(runD) & set(runH))
    say("=" * 70)
    say(f"DIRECT DCM vs HCM, {et}   (NF is NOT used as a reference here)")
    say(f"  events in both runs: {len(shared)}")
    say()

    # ---- QC: NF PSI must agree between the two runs ---------------------
    diffs, corrs = [], []
    for k in shared[:4000]:
        n1 = runD[k][1]
        n2 = runH[k][1]
        m = np.isfinite(n1) & np.isfinite(n2)
        if m.sum() > 20:
            diffs.append(np.nanmean(np.abs(n1[m] - n2[m])))
            if np.std(n1[m]) > 1e-6 and np.std(n2[m]) > 1e-6:
                corrs.append(np.corrcoef(n1[m], n2[m])[0, 1])
    say("--- QC: the 166 NF samples appear in BOTH runs; their PSI must match ---")
    if diffs:
        say(f"  mean |PSI difference| across {len(diffs)} events: {np.mean(diffs):.5f}")
        say(f"  median per-event correlation: {np.median(corrs):.5f}")
        if np.mean(diffs) > 0.01:
            say("  *** PSI does not agree between runs. The two rMATS runs are not")
            say("      comparable and everything below is invalid. Stop here.")
            open(f"{a.out}.report.txt", "w").write("\n".join(lines) + "\n")
            return
        say("  OK -- cross-run PSI is comparable, so a direct contrast is valid.")
    say()

    # ---- covariates -----------------------------------------------------
    covD = read_cov(a.cov_dcm)
    covH = read_cov(a.cov_hcm)
    meta = {}
    for s in dcm_ids:
        if s in covD.index:
            r = covD.loc[s]
            meta[s] = (1.0, str(r.get("sex", "")), r.get("age", np.nan), r.get("RIN", np.nan))
    for s in hcm_ids:
        r = covH.loc[s] if s in covH.index else (covD.loc[s] if s in covD.index else None)
        if r is not None:
            meta[s] = (0.0, str(r.get("sex", "")), r.get("age", np.nan), r.get("RIN", np.nan))
    n_missing = len([s for s in dcm_ids + hcm_ids if s not in meta])
    if n_missing:
        say(f"[warn] {n_missing} samples have no covariates and are dropped")

    order = dcm_ids + hcm_ids
    grp = np.array([meta[s][0] if s in meta else np.nan for s in order])
    sexm = np.array([1.0 if s in meta and meta[s][1].upper().startswith("M") else
                     (0.0 if s in meta else np.nan) for s in order])
    age = np.array([meta[s][2] if s in meta else np.nan for s in order], float)
    rin = np.array([meta[s][3] if s in meta else np.nan for s in order], float)

    # ---- per-event model -------------------------------------------------
    rows = []
    skips = {"low HCM n": 0, "low PSI sd": 0, "singular/df": 0}
    for k in shared:
        psi = np.concatenate([runD[k][0], runH[k][0]])
        ok = np.isfinite(psi) & np.isfinite(grp) & np.isfinite(sexm) & \
             np.isfinite(age) & np.isfinite(rin)
        if ok.sum() < 30:
            skips["low HCM n"] += 1
            continue
        nh = int(np.sum(ok & (grp == 0)))
        nd = int(np.sum(ok & (grp == 1)))
        if nh < a.min_samples_hcm or nd < 20:
            skips["low HCM n"] += 1
            continue
        y = psi[ok]
        if np.std(y) < a.min_psi_sd:
            skips["low PSI sd"] += 1
            continue
        X = np.column_stack([np.ones(ok.sum()), grp[ok], sexm[ok], age[ok], rin[ok]])
        keep = [0] + [i for i in range(1, X.shape[1]) if np.std(X[:, i]) > 1e-9]
        X = X[:, keep]
        if 1 not in keep or np.linalg.matrix_rank(X) < X.shape[1]:
            skips["singular/df"] += 1
            continue
        dof = len(y) - X.shape[1]
        if dof <= 5:
            skips["singular/df"] += 1
            continue
        inv = np.linalg.pinv(X.T @ X)
        b = inv @ X.T @ y
        res = y - X @ b
        s2 = float(res @ res) / dof
        if not np.isfinite(s2) or s2 <= 0:
            skips["singular/df"] += 1
            continue
        se = np.sqrt(np.maximum(np.diag(s2 * inv), 0))
        j = keep.index(1)
        if se[j] <= 0:
            skips["singular/df"] += 1
            continue
        t = b[j] / se[j]
        p = 2 * stats.t.sf(abs(t), dof)
        g, s_, e_, ch = runD[k][2], runD[k][3], runD[k][4], runD[k][5]
        rows.append((k, g, ch, s_, e_, nd, nh, b[j], se[j], p))

    res = pd.DataFrame(rows, columns=["key", "gene", "chr", "exonStart", "exonEnd",
                                      "n_DCM", "n_HCM", "dPSI_DCMvsHCM", "se", "pvalue"])
    if len(res) == 0:
        sys.exit("[fatal] no events tested")
    res["FDR"] = bh(res.pvalue.values)
    res["mde80"] = 2.802 * res["se"]
    res.to_csv(f"{a.out}.all_events.tsv", sep="\t", index=False)

    say(f"--- Direct contrast: PSI ~ DCM_vs_HCM + sex + age + RIN ---")
    say(f"  events tested: {len(res)}   (skipped: {skips})")
    sig = res[(res.FDR < 0.05) & (res.dPSI_DCMvsHCM.abs() >= 0.10)]
    say(f"  significant, FDR<0.05 and |dPSI|>=0.10: {len(sig)}")
    say(f"  significant at FDR<0.05, any effect size: {(res.FDR < 0.05).sum()}")
    say()

    genes = [g.strip() for g in a.genes.split(",") if g.strip()]
    prog = res[res.gene.isin(genes)]
    psig = prog[(prog.FDR < 0.05) & (prog.dPSI_DCMvsHCM.abs() >= 0.10)]
    say(f"--- RBM20 program, DCM vs HCM directly ---")
    say(f"  events tested: {len(prog)}   significant: {len(psig)}")
    if len(psig):
        say(f"  {'gene':<10}{'chr':<7}{'exonStart':>11}{'dPSI':>9}{'FDR':>11}")
        for _, r in psig.sort_values("pvalue").head(20).iterrows():
            say(f"  {r.gene:<10}{r.chr:<7}{r.exonStart:>11}"
                f"{r.dPSI_DCMvsHCM:>+9.3f}{r.FDR:>11.2e}")
    say()

    # ---- equivalence: is DCM-HCM SMALLER than DCM-NF? -------------------
    if a.disc:
        d = pd.read_csv(a.disc, sep="\t")
        m = res.merge(d[["key", "dPSI_adj", "FDR"]], on="key",
                      suffixes=("", "_disc"))
        m = m[(m.FDR_disc < 0.05) & (m.dPSI_adj.abs() >= 0.10)]
        say("--- Equivalence: DCM-vs-HCM effect against DCM-vs-NF effect ---")
        say(f"  on the {len(m)} DCM-significant events also testable here:")
        if len(m):
            ratio = np.abs(m.dPSI_DCMvsHCM.values) / np.abs(m.dPSI_adj.values)
            say(f"    median |dPSI(DCM-HCM)| / |dPSI(DCM-NF)| = {np.median(ratio):.3f}")
            say(f"    median |dPSI(DCM-NF)|  = {np.median(np.abs(m.dPSI_adj)):.3f}")
            say(f"    median |dPSI(DCM-HCM)| = {np.median(np.abs(m.dPSI_DCMvsHCM)):.3f}")
            powered = m[m.mde80 <= np.abs(m.dPSI_adj)]
            small = powered[np.abs(powered.dPSI_DCMvsHCM) < 0.5 * np.abs(powered.dPSI_adj)]
            say(f"    events where this test COULD have seen a DCM-NF-sized "
                f"difference: {len(powered)}/{len(m)}")
            say(f"    of those, DCM-HCM difference is <half the DCM-NF effect: "
                f"{len(small)} ({100*len(small)/max(len(powered),1):.0f}%)")
            say()
            say("  A ratio near 0 means DCM and HCM are indistinguishable on the very")
            say("  events that define the disease signature -- i.e. the signature is")
            say("  'failing myocardium vs donor', not 'DCM'. A ratio near 1 would mean")
            say("  DCM and HCM differ as much as DCM differs from NF.")
        mp = m[m.gene.isin(genes)]
        if len(mp):
            ratio_p = np.abs(mp.dPSI_DCMvsHCM.values) / np.abs(mp.dPSI_adj.values)
            say(f"  RBM20 program only (n={len(mp)}): median ratio = "
                f"{np.median(ratio_p):.3f}")
    say("=" * 70)

    with open(f"{a.out}.report.txt", "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"\nwrote {a.out}.all_events.tsv and {a.out}.report.txt")


if __name__ == "__main__":
    main()
