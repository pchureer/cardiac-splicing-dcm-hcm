#!/usr/bin/env python3
"""
34_meta_analysis.py -- the two meta-analyses. Phase 1.

TWO SEPARATE POOLS, NEVER COMBINED
-----------------------------------
  META A  the DCM signature   MAGNet + Colorado          (DCM vs NF)
  META B  the HCM transfer    Brigham + Fuwai + Stanford (HCM vs control)

A DCM-vs-NF contrast and an HCM-vs-control contrast answer different questions.
Pooling them would be a real methodological error, so this script cannot do it.

WHY META A HAS ONLY TWO COHORTS
--------------------------------
LBNL was the third. The genotype check found 24 of its 33 samples are the same
individuals as MAGNet samples, so pooling it would combine the discovery cohort with
itself and manufacture significance out of nothing. It is excluded here and has to
stay excluded. See RESULT_sample_identity_LBNL_overlap.md.

k=2 IS A REAL LIMITATION, STATED NOT HIDDEN
--------------------------------------------
With two studies the between-study variance tau^2 is estimated from a single degree of
freedom, so I^2 is close to uninformative -- it will read 0% or something large with
little in between, and neither should be interpreted as evidence about heterogeneity.
The script reports I^2 for META A but labels it unreliable. If you would rather not
defend a two-study meta-analysis, the honest alternative is to report Colorado as a
direct replication and drop META A; the per-event pooled estimates below are still
useful as a supplementary table either way.

MAGNET IS THE DISCOVERY COHORT. Pooling discovery with replication and then reporting
"significant in the meta-analysis" is circular for events already significant in MAGNet.
The non-circular result is the LAST column: events that reach significance POOLED but
were NOT significant in MAGNet alone. Those are genuinely new, and they are what META A
is for.

EXCLUSIONS APPLIED AUTOMATICALLY
  Colorado MXE   reverse slope 0.061 (n=1477) against MAGNet's SE 0.447 -- does not
                 reproduce; RESULT_colorado_mxe_excluded.md
  Fuwai MXE      flagged caution (reverse slope 0.281); included but marked, and the
                 SE-only pooled figure is reported alongside

Usage:
  python 34_meta_analysis.py
  python 34_meta_analysis.py --types SE            # SE only, the conservative report
  python 34_meta_analysis.py --model fixed         # fixed-effect sensitivity
"""
import argparse
import os

import numpy as np
import pandas as pd
from scipy import stats

NKEEP = {"SE": 3, "MXE": 5, "A3SS": 5, "A5SS": 5, "RI": 3}

# 160: only the MAGNet arm re-bases with the read length. Colorado was never
# re-quantified, and META_B is Brigham/Fuwai/Stanford — all 50/101/150 bp HCM
# cohorts, untouched by the MAGNet re-quantification. Same asymmetry that
# 129_patch_21_replication.py applied to script 21: discovery suffixed,
# replication cohorts not.
MSUF = os.environ.get("META_SUFFIX", "")
META_A = [("MAGNet", "{ET}_v2filtered" + MSUF),
          ("Colorado", "{ET}_Colorado_DCMvsNF")]
META_B = [("Brigham", "Brigham_{ET}_HCM"), ("Fuwai", "Fuwai_{ET}_HCM"),
          ("Stanford", "Stanford_{ET}_HCM")]
EXCLUDE = {("Colorado", "MXE")}
CAUTION = {("Fuwai", "MXE")}


def load(base, et):
    for p in (f"{base}_dedup.all_events.tsv", f"{base}.all_events.tsv"):
        if os.path.exists(p) and os.path.getsize(p):
            t = pd.read_csv(p, sep="\t")
            d = next((c for c in ("dPSI_adj", "dPSI") if c in t.columns), None)
            if d is None or "se_adj" not in t.columns or "key" not in t.columns:
                return None, f"{p} (missing dPSI/se_adj/key)"
            n = NKEEP[et]
            t["ekey"] = [":".join(str(k).split(":")[:n]) for k in t["key"]]
            t = t.loc[t.groupby("ekey")["FDR"].idxmin()].copy()
            t = t[["ekey", "gene", d, "se_adj", "FDR"]].rename(
                columns={d: "y", "se_adj": "se"})
            t = t[np.isfinite(t.y) & np.isfinite(t.se) & (t.se > 0)]
            return t, p
    return None, f"{base}[_dedup].all_events.tsv not found"


def dersimonian_laird(Y, S, random=True):
    """Y, S: n_events x k arrays (nan where a cohort lacks the event)."""
    W = 1.0 / S ** 2
    ok = np.isfinite(Y) & np.isfinite(W)
    W = np.where(ok, W, 0.0)
    Yz = np.where(ok, Y, 0.0)
    k = ok.sum(axis=1)
    sw = W.sum(axis=1)
    fixed = np.divide((W * Yz).sum(axis=1), sw, out=np.full(len(Y), np.nan),
                      where=sw > 0)
    Q = (W * (Yz - fixed[:, None]) ** 2 * ok).sum(axis=1)
    if not random:
        return fixed, np.sqrt(1 / np.where(sw > 0, sw, np.nan)), np.zeros(len(Y)), Q, k
    sw2 = (W ** 2).sum(axis=1)
    denom = sw - np.divide(sw2, sw, out=np.zeros_like(sw), where=sw > 0)
    tau2 = np.maximum(0.0, np.divide(Q - (k - 1), denom,
                                     out=np.zeros(len(Y)), where=denom > 0))
    Wr = np.where(ok, 1.0 / (S ** 2 + tau2[:, None]), 0.0)
    swr = Wr.sum(axis=1)
    pooled = np.divide((Wr * Yz).sum(axis=1), swr, out=np.full(len(Y), np.nan),
                       where=swr > 0)
    se = np.sqrt(np.divide(1.0, swr, out=np.full(len(Y), np.nan), where=swr > 0))
    return pooled, se, tau2, Q, k


def bh(p):
    p = np.asarray(p, float)
    ok = np.isfinite(p)
    q = np.full(p.shape, np.nan)
    v = p[ok]
    o = np.argsort(v)
    n = len(v)
    adj = np.minimum.accumulate((v[o] * n / (np.arange(n) + 1))[::-1])[::-1]
    out = np.empty(n)
    out[o] = np.minimum(adj, 1.0)
    q[ok] = out
    return q


def run(name, cohorts, types, a, P):
    P()
    P("=" * 94)
    P(f"{name}")
    P("=" * 94)
    allrows = []
    for et in types:
        frames, usedc = {}, []
        for coh, pat in cohorts:
            if (coh, et) in EXCLUDE:
                P(f"  [{et}] {coh} EXCLUDED (does not reproduce in reverse replication)")
                continue
            t, p = load(pat.format(ET=et), et)
            if t is None:
                P(f"  [{et}] {coh}: {p}")
                continue
            frames[coh] = t
            usedc.append(coh)
            if (coh, et) in CAUTION:
                P(f"  [{et}] {coh} included but flagged CAUTION")
        if len(frames) < 2:
            P(f"  [{et}] fewer than 2 cohorts -- cannot pool")
            continue
        keys = set.intersection(*[set(f.ekey) for f in frames.values()])
        if not keys:
            P(f"  [{et}] no shared events")
            continue
        keys = sorted(keys)
        idx = pd.Index(keys, name="ekey")
        Y = np.column_stack([frames[c].set_index("ekey").reindex(idx).y.values
                             for c in usedc])
        S = np.column_stack([frames[c].set_index("ekey").reindex(idx).se.values
                             for c in usedc])
        pooled, se, tau2, Q, k = dersimonian_laird(Y, S, random=(a.model == "random"))
        z = pooled / se
        pv = 2 * stats.norm.sf(np.abs(z))
        q = bh(pv)
        with np.errstate(invalid="ignore", divide="ignore"):
            I2 = np.where(Q > 0, np.maximum(0.0, (Q - (k - 1)) / Q) * 100, 0.0)
        sig = (q < a.fdr) & (np.abs(pooled) >= a.min_dpsi)

        R = pd.DataFrame({"ekey": keys, "event_type": et,
                          "gene": frames[usedc[0]].set_index("ekey")
                          .reindex(idx).gene.values,
                          "pooled_dPSI": pooled, "se": se, "ci_lo": pooled - 1.96 * se,
                          "ci_hi": pooled + 1.96 * se, "p": pv, "FDR": q,
                          "I2": I2, "tau2": tau2, "k": k, "significant": sig})
        for c in usedc:
            sub = frames[c].set_index("ekey").reindex(idx)
            R[f"dPSI_{c}"] = sub.y.values
            R[f"FDR_{c}"] = sub.FDR.values
        allrows.append(R)

        line = (f"  {et:<5} cohorts {'+'.join(usedc):<28} events {len(R):>6}"
                f"   pooled sig {int(sig.sum()):>5}")
        if "MAGNet" in usedc:
            disc = R[f"FDR_MAGNet"] < a.fdr
            discs = disc & (R.dPSI_MAGNet.abs() >= a.min_dpsi)
            new = int((sig & ~discs).sum())
            lost = int((~sig & discs).sum())
            line += f"   NEW (not sig in MAGNet) {new:>4}   lost {lost:>4}"
        P(line)
        P(f"        I2: median {np.nanmedian(I2):.0f}%  "
          f">50% in {int((I2 > 50).sum())} events ({100 * (I2 > 50).mean():.0f}%)"
          + ("   [k=2: I2 UNRELIABLE, do not interpret]" if k.max() <= 2 else ""))

    if not allrows:
        P("  nothing pooled")
        return None
    A = pd.concat(allrows, ignore_index=True)
    # 160: the NEW TOTAL is what [047] quotes as "51", and this script never
    # computed it — it printed the four per-class counts and stopped, so no
    # grep of this file or its output could find the published number.
    # Recomputed from the concatenated frame, which is identical to summing
    # the per-class counts because both predicates are per-row.
    if "FDR_MAGNet" in A.columns:
        _dsig = ((A.FDR_MAGNet < a.fdr)
                 & (A.dPSI_MAGNet.abs() >= a.min_dpsi))
        P(f"  NEW TOTAL across classes: {int((A.significant & ~_dsig).sum())}"
          f"   <- the count [047] quotes")
    # 160: the suffix MUST reach the output name. Without it a 125 bp run
    # overwrites meta_a.tsv in place, which 84_pool_b.py reads.
    out = f"meta_{name.split()[1].lower()}{MSUF}.tsv"
    A.to_csv(out, sep="\t", index=False, float_format="%.5g")
    P(f"  -> {out}   ({len(A)} events, {int(A.significant.sum())} significant)")
    return A


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--types", default="SE,MXE,A3SS,A5SS,RI")
    ap.add_argument("--fdr", type=float, default=0.05)
    ap.add_argument("--min-dpsi", type=float, default=0.10)
    ap.add_argument("--model", choices=["random", "fixed"], default="random")
    ap.add_argument("--out", default="meta_analysis.txt")
    a = ap.parse_args()
    types = [t for t in a.types.split(",") if t in NKEEP]
    if MSUF:
        # 148: a suffix that re-bases nothing exits 0 and publishes a 50 bp
        # answer under a 125 bp name.
        if not [p for et in types
                for p in (f"{et}_v2filtered{MSUF}_dedup.all_events.tsv",
                          f"{et}_v2filtered{MSUF}.all_events.tsv")
                if os.path.exists(p)]:
            raise SystemExit(
                f"REFUSING: META_SUFFIX={MSUF!r} but no "
                f"{{ET}}_v2filtered{MSUF}[_dedup].all_events.tsv exists for any "
                f"of {types}. The suffix would re-base nothing and the run "
                f"would report 50 bp numbers under a 125 bp name.")
        if a.out == ap.get_default("out"):
            raise SystemExit(
                f"REFUSING: META_SUFFIX is set but --out is still the default "
                f"{a.out!r}. Pass --out meta_analysis{MSUF}.txt, or the 50 bp "
                f"report is overwritten in place.")
    fh = open(a.out, "w")

    def P(s=""):
        print(s)
        fh.write(s + "\n")

    P(f"model: {a.model}-effects (DerSimonian-Laird)   "
      f"significance: FDR<{a.fdr} and |pooled dPSI|>={a.min_dpsi}")
    P("exon-level keys throughout, matching the locked convention")
    P(f"read length: MAGNet arm = {{ET}}_v2filtered{MSUF} "
      f"(META_SUFFIX={MSUF!r}); Colorado and all of META B are NOT "
      f"re-quantified")

    A = run("META A -- the DCM signature (MAGNet + Colorado, DCM vs NF)",
            META_A, types, a, P)
    B = run("META B -- the HCM transfer (Brigham + Fuwai + Stanford, HCM vs control)",
            META_B, types, a, P)

    P()
    P("=" * 94)
    P("HOW TO REPORT THIS")
    P("=" * 94)
    P("  META A is TWO cohorts, one of which is the discovery cohort. Report the NEW")
    P("  column -- events significant pooled but not in MAGNet alone. Those are not")
    P("  circular. Do NOT report 'N events significant in the meta-analysis' as if it")
    P("  were independent of discovery, and do not interpret I^2 at k=2.")
    P()
    P("  META B is a genuine three-cohort meta-analysis. Its I^2 is interpretable.")
    P("  High I^2 concentrated in Fuwai points at its GENCODE v46 index rather than")
    P("  biology -- check by rerunning with --types SE and comparing.")
    P()
    P("  LBNL is absent from both by design: 24 of its 33 samples are MAGNet hearts.")
    P("  Adding it would pool the discovery cohort with itself.")
    P()
    if A is not None and B is not None:
        P("  NEVER pool META A with META B. DCM-vs-NF and HCM-vs-control answer")
        P("  different questions; a combined figure would be uninterpretable.")
    P(f"\n[write] {a.out}")
    fh.close()


if __name__ == "__main__":
    main()
