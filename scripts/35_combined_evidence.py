#!/usr/bin/env python3
"""
35_combined_evidence.py -- combine cohorts by EVIDENCE, not by effect size.
Replaces 34_meta_analysis.py as the primary Phase 1 result; keep 34's output as a
supplementary sensitivity analysis.

WHY NOT INVERSE-VARIANCE META-ANALYSIS
---------------------------------------
An inverse-variance meta-analysis assumes every cohort estimates the SAME underlying
quantity, and pools the dPSI values directly. This project has already shown that
assumption is false here: cross-cohort dPSI magnitude does not transfer. Colorado's
transfer slope against MAGNet is 0.521, the three HCM cohorts run 0.24-0.48, and
RESULT_magnitude_not_estimable.md concludes magnitude is not reliably estimable
between independently generated splicing datasets.

Script 34 shows the damage this does. Of 471 MAGNet-significant SE events measured in
Colorado, only 191 survived pooling; 280 were lost. Those are mostly events near the
0.10 effect threshold in MAGNet whose Colorado estimate is systematically about half as
large, dragging the pooled estimate under the threshold. That is an artefact of pooling
non-comparable magnitudes, not disagreement about whether the event is real -- and a
reviewer who reads the magnitude section will say so.

WHAT THIS DOES INSTEAD
-----------------------
Stouffer's weighted Z combines TEST STATISTICS and DIRECTIONS. Each cohort contributes

    z_i = dPSI_i / se_i          (a standardised statistic: scale-free by construction)
    w_i = sqrt(n_i)              (sample size, NOT 1/se -- weighting by precision would
                                  quietly reintroduce the magnitude assumption)
    Z   = sum(w_i z_i) / sqrt(sum w_i^2)

No common effect size is assumed anywhere. The question becomes "do the cohorts agree
that this event changes, and in which direction", which is exactly the claim the paper
makes and the one the data support.

SIGNIFICANCE IS DEFINED WITHOUT ANY MAGNITUDE THRESHOLD
  FDR < 0.05 on the combined Z, AND every contributing cohort agrees in direction.
Unanimous direction is a strong requirement and it is the paper's actual claim. Counts
without it are printed alongside so the cost of the requirement is visible.

MAGNITUDE IS STILL REPORTED, BUT ONLY WITHIN A COHORT. The median |dPSI| in the
DISCOVERY cohort is a legitimate descriptive statistic; a pooled cross-cohort dPSI is
not. That distinction is the whole point.

EXCLUSIONS (same as 34): Colorado MXE excluded, Fuwai MXE flagged, LBNL absent.

Usage:
  python 35_combined_evidence.py
  python 35_combined_evidence.py --no-unanimous     # drop the direction requirement
"""
import argparse
import os

import numpy as np
import pandas as pd
from scipy import stats

NKEEP = {"SE": 3, "MXE": 5, "A3SS": 5, "A5SS": 5, "RI": 3}
# 161: only the MAGNet arm re-bases with the read length. Colorado was never
# re-quantified and POOL_B is the three HCM cohorts.
#
# RULE 9 — the suffix parameter is NOT uniform. This script's variable is
# COMB_SUFFIX, deliberately NOT the META_SUFFIX that 160 gave script 34.
# They are different producers writing different tables; one shared variable
# would re-base both from a single export and make it impossible to run one
# read length through 34 and another through 35 while diagnosing a
# disagreement between them — which is exactly what [047]'s claim requires.
CSUF = os.environ.get("COMB_SUFFIX", "")
POOL_A = [("MAGNet", "{ET}_v2filtered" + CSUF),
          ("Colorado", "{ET}_Colorado_DCMvsNF")]
# 161: the per-class NEW counts are printed but never summed, and [047]'s
# "48" IS that sum. Keyed by pool name so POOL A and POOL B cannot mix.
_NEW_TOTAL = {}
POOL_B = [("Brigham", "Brigham_{ET}_HCM"), ("Fuwai", "Fuwai_{ET}_HCM"),
          ("Stanford", "Stanford_{ET}_HCM")]
EXCLUDE = {("Colorado", "MXE")}
CAUTION = {("Fuwai", "MXE")}


def load(base, et):
    for p in (f"{base}_dedup.all_events.tsv", f"{base}.all_events.tsv"):
        if os.path.exists(p) and os.path.getsize(p):
            t = pd.read_csv(p, sep="\t")
            d = next((c for c in ("dPSI_adj", "dPSI") if c in t.columns), None)
            need = {"key", "se_adj", "FDR"}
            if d is None or need - set(t.columns):
                return None, f"{p} (missing columns)"
            n = NKEEP[et]
            t["ekey"] = [":".join(str(k).split(":")[:n]) for k in t["key"]]
            t = t.loc[t.groupby("ekey")["FDR"].idxmin()].copy()
            ncol = "n" if "n" in t.columns else None
            keep = ["ekey", "gene", d, "se_adj", "FDR"] + ([ncol] if ncol else [])
            t = t[keep].rename(columns={d: "y", "se_adj": "se", ncol: "n"}
                               if ncol else {d: "y", "se_adj": "se"})
            if ncol is None:
                t["n"] = np.nan
            t = t[np.isfinite(t.y) & np.isfinite(t.se) & (t.se > 0)]
            return t, p
    return None, f"{base}[_dedup].all_events.tsv not found"


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
    P("=" * 98)
    P(name)
    P("=" * 98)
    P(f"  {'type':<6}{'cohorts':<24}{'events':>8}{'FDR':>7}{'+effect':>9}"
      f"{'+unanim':>9}{'NEW':>6}{'med |dPSI| disc':>17}")
    allr = []
    for et in types:
        frames, used = {}, []
        for coh, pat in cohorts:
            if (coh, et) in EXCLUDE:
                P(f"  {et:<6}{coh} EXCLUDED (fails reverse replication)")
                continue
            t, p = load(pat.format(ET=et), et)
            if t is None:
                continue
            frames[coh] = t.set_index("ekey")
            used.append(coh)
        if len(frames) < 2:
            P(f"  {et:<6}fewer than 2 cohorts -- cannot combine")
            continue
        keys = sorted(set.intersection(*[set(f.index) for f in frames.values()]))
        if not keys:
            continue
        idx = pd.Index(keys, name="ekey")
        Y = np.column_stack([frames[c].reindex(idx).y.values for c in used])
        S = np.column_stack([frames[c].reindex(idx).se.values for c in used])
        N = np.column_stack([frames[c].reindex(idx).n.values for c in used])
        N = np.where(np.isfinite(N) & (N > 0), N, 1.0)

        z = Y / S
        w = np.sqrt(N)
        Z = (w * z).sum(axis=1) / np.sqrt((w ** 2).sum(axis=1))
        pv = 2 * stats.norm.sf(np.abs(Z))
        q = bh(pv)
        sgn = np.sign(Y)
        unan = (np.abs(sgn.sum(axis=1)) == sgn.shape[1]) & (sgn != 0).all(axis=1)
        sig_any = q < a.fdr
        absY = np.abs(Y)
        if a.min_dpsi <= 0:
            eff = np.ones(len(Y), bool)
        elif a.dpsi_rule == "all":
            eff = (absY >= a.min_dpsi).all(axis=1)
        elif a.dpsi_rule == "discovery":
            eff = absY[:, 0] >= a.min_dpsi
        else:
            eff = (absY >= a.min_dpsi).any(axis=1)
        sig_eff = sig_any & eff
        sig = sig_eff & unan if a.unanimous else sig_eff

        disc = used[0]
        dY = frames[disc].reindex(idx)
        R = pd.DataFrame({
            "ekey": keys, "event_type": et, "effect_pass": eff,
            "gene": frames[used[0]].reindex(idx).gene.values,
            "Z": Z, "p": pv, "FDR": q,
            "unanimous_direction": unan,
            "direction": np.where(Z > 0, "+", "-"),
            "significant": sig,
            f"dPSI_{disc}": dY.y.values})
        for c in used:
            s = frames[c].reindex(idx)
            R[f"dPSI_{c}"] = s.y.values
            R[f"FDR_{c}"] = s.FDR.values
        allr.append(R)

        new = ""
        if disc == "MAGNet":
            dsig = (dY.FDR.values < a.fdr) & (np.abs(dY.y.values) >= 0.10)
            new = f"{int((sig & ~dsig).sum())}"
            _NEW_TOTAL[name] = (_NEW_TOTAL.get(name, 0)
                                + int((sig & ~dsig).sum()))
        med = (np.nanmedian(np.abs(dY.y.values[sig])) if sig.any() else np.nan)
        flag = "  [CAUTION]" if any((c, et) in CAUTION for c in used) else ""
        P(f"  {et:<6}{'+'.join(used):<24}{len(R):>8}{int(sig_any.sum()):>7}"
          f"{int(sig_eff.sum()):>9}{int(sig.sum()):>9}{new:>6}{med:>17.3f}{flag}")

    if not allr:
        P("  nothing combined")
        return None
    A = pd.concat(allr, ignore_index=True)
    # 161: the sum [047] quotes. This script printed the per-class NEW counts
    # and stopped, so no grep of it or its output could find "48".
    if name in _NEW_TOTAL:
        P(f"  NEW TOTAL across classes: {_NEW_TOTAL[name]}"
          f"   <- the count [047] quotes for this pool")
    # 161: the suffix MUST reach the output name. Without it a 125 bp run
    # overwrites combined_a.tsv in place — the only copy of the Stouffer
    # pool, read by 84_pool_b.py. That is Rule 10, instance 1.
    out = f"combined_{name.split()[1].lower()}{CSUF}.tsv"
    A.to_csv(out, sep="\t", index=False, float_format="%.5g")
    P(f"  -> {out}   {len(A)} events, {int(A.significant.sum())} significant")
    P(f"     FDR = FDR<{a.fdr} alone; +effect additionally requires |dPSI|>="
      f"{a.min_dpsi} ({a.dpsi_rule}); +unanim additionally requires every cohort to")
    P(f"     agree in direction. Each column shows what the next requirement costs.")
    return A


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--types", default="SE,MXE,A3SS,A5SS,RI")
    ap.add_argument("--fdr", type=float, default=0.05)
    ap.add_argument("--min-dpsi", type=float, default=0.10,
                    help="Effect threshold applied WITHIN a cohort. Pooling magnitudes "
                         "across cohorts is invalid here; requiring a real effect "
                         "inside a cohort is not, and it is what the project's locked "
                         "926/1,445 counts use. Dropping it entirely admits thousands "
                         "of events with median |dPSI| ~0.03.")
    ap.add_argument("--dpsi-rule", choices=["any", "discovery", "all"], default="any",
                    help="any: |dPSI|>=threshold in at least one cohort (default). "
                         "discovery: in the first-listed cohort only. "
                         "all: in every cohort -- too strict here, because attenuation "
                         "means a 0.10 MAGNet effect is ~0.05 in Colorado.")
    ap.add_argument("--no-unanimous", dest="unanimous", action="store_false")
    ap.add_argument("--out", default="combined_evidence.txt")
    a = ap.parse_args()
    types = [t for t in a.types.split(",") if t in NKEEP]
    if CSUF:
        # 148: a suffix that re-bases nothing exits 0 and reports 50 bp
        # numbers under a 125 bp name.
        if not [p for et in types
                for p in (f"{et}_v2filtered{CSUF}_dedup.all_events.tsv",
                          f"{et}_v2filtered{CSUF}.all_events.tsv")
                if os.path.exists(p)]:
            raise SystemExit(
                f"REFUSING: COMB_SUFFIX={CSUF!r} but no "
                f"{{ET}}_v2filtered{CSUF}[_dedup].all_events.tsv exists for any "
                f"of {types}. The suffix would re-base nothing.")
        if a.out == ap.get_default("out"):
            raise SystemExit(
                f"REFUSING: COMB_SUFFIX is set but --out is still the default "
                f"{a.out!r}. Pass --out combined_evidence{CSUF}.txt, or the "
                f"28 Sep report is overwritten in place.")
    fh = open(a.out, "w")

    def P(s=""):
        print(s)
        fh.write(s + "\n")

    P("Stouffer weighted-Z on standardised statistics; weights sqrt(n).")
    P(f"read length: MAGNet arm = {{ET}}_v2filtered{CSUF} "
      f"(COMB_SUFFIX={CSUF!r}); Colorado and all of POOL B are NOT "
      f"re-quantified")
    P(f"dpsi-rule: {a.dpsi_rule}   (84_pool_b.py's REC_A record is the "
      f"`discovery` rule; this script defaults to `any`, and the 28 Sep "
      f"combined_evidence.txt and combined_a.tsv are five minutes apart, so "
      f"they may be from different invocations)")
    P("NO cross-cohort effect size is pooled anywhere -- see the docstring.")
    P(f"significant = FDR<{a.fdr} AND |dPSI|>={a.min_dpsi} ({a.dpsi_rule} cohort)"
      + (" AND unanimous direction" if a.unanimous else ""))
    P("The effect threshold is applied WITHIN a cohort, never to a pooled effect size.")
    P("'med |dPSI| disc' is the median absolute effect IN THE FIRST-LISTED COHORT,")
    P("a within-cohort descriptive. There is deliberately no pooled dPSI column.")

    run("POOL A -- the DCM signature (MAGNet + Colorado, DCM vs NF)", POOL_A, types, a, P)
    run("POOL B -- the HCM transfer (Brigham + Fuwai + Stanford, HCM vs control)",
        POOL_B, types, a, P)

    P()
    P("=" * 98)
    P("HOW TO REPORT")
    P("=" * 98)
    P("  This is a combined-evidence analysis, not a meta-analysis of effect sizes, and")
    P("  it should be described that way in the methods. The reason is in this paper's")
    P("  own results: dPSI magnitude does not transfer across independently generated")
    P("  datasets (transfer slopes 0.24-0.52), so the equal-effect assumption behind")
    P("  inverse-variance pooling does not hold. Stouffer's method assumes only that")
    P("  each cohort tests the same null.")
    P()
    P("  POOL A's NEW column -- significant combined but not significant in MAGNet")
    P("  alone -- is the non-circular result. Everything else in POOL A partly restates")
    P("  discovery, because MAGNet is in the pool and carries the largest weight.")
    P()
    P("  POOL B contains no discovery cohort, so every event in it is a result.")
    P()
    P("  34_meta_analysis.py's inverse-variance output is retained as a supplementary")
    P("  sensitivity analysis. If the two broadly agree, say so; where they differ,")
    P("  the difference is the magnitude assumption and this one is the defensible one.")
    P(f"\n[write] {a.out}")
    fh.close()


if __name__ == "__main__":
    main()
