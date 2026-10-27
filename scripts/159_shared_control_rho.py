#!/usr/bin/env python3
"""
159_shared_control_rho.py -- A8a: [042]'s shared-control rho = 0.27.

THE SENTENCE, AND WHY IT DOES NOT MATCH ITS PRODUCER
----------------------------------------------------
[042], verbatim:

  "Where two contrasts share a control group, their dPSI estimates subtract
   the same control mean and are positively correlated under a true null.
   That baseline was computed from the per-event group sizes and used in
   place of zero (Supplementary Methods). For MAGNet's 166 shared controls
   it is rho = 0.27."

The producer is `05_hcm_specificity.py` (lines 118-128):

    n_d  = measured.n_g1_DCM     # DCM group size, per event
    n_h  = measured.n_g1_HCM     # HCM group size, per event
    n_nf = measured.n_g2_DCM     # shared NF group size, per event
    rho0 = (1 / n_nf) / np.sqrt((1 / n_d + 1 / n_nf) * (1 / n_h + 1 / n_nf))
    rho_null = float(np.median(rho0))

Run on the published 50 bp tables, 05 prints **+0.265**, and says "the same
**152** NF samples", not 166. Two mismatches with the sentence:

  * **152, not 166.** 166 is the nominal group size; 152 is the per-event
    median after coverage dropout. The sentence names the nominal n but
    describes the per-event computation.
  * **25 events, not the signature.** 05 restricts to the nine RBM20
    programme genes (`--genes`) BEFORE computing the baseline:
    1,074 significant -> 30 in those genes -> 25 measured in the HCM run.
    A Methods statement about the shared-control baseline reads as a
    property of the design, not of a nine-gene subset.

0.27 IS the analytic value at the NOMINAL sizes: with n_DCM=166, n_HCM=28,
n_NF=166 the formula gives **0.2686**, which 05's own docstring states as
"about +0.27". So the published number is most likely the docstring
expectation, not the script's output — and if so it is read-length invariant,
because nominal group sizes do not change.

This script measures all three readings side by side at both read lengths and
refuses to choose for you. It does not search for a variant that reproduces
0.27: it reports what each defensible reading gives, having first proved it
can reproduce 05's own number.

05'S CONVENTION, WHICH THIS SCRIPT COPIES EXACTLY
-------------------------------------------------
05 does **no collapse**. `dsig` is raw rows passing FDR and |dPSI|, which is
why its count is 1,074 at 50 bp and 944 at 125 bp — not 926/813
(collapse-first) and not 991/864 (filter-first). Those same 1,074 and 944
appear as the "MAGNet sig" column of `RESULT_21_replication_rl125.md`. The
merge is on the **full** `key`, not the NKEEP-collapsed one. Reproducing 05's
0.265 is the gate; a script using the project's usual collapse would silently
compute a different quantity and still look plausible.

USAGE
-----
  python3 159_shared_control_rho.py \
      --dcm50  SE_v2filtered.all_events.tsv  --hcm50  SE_HCM_v2.all_events.tsv \
      --dcm125 SE_v2filtered_rl125.all_events.tsv \
      --hcm125 SE_HCMvsNF_rl125.all_events.tsv \
      --rec50 0.265 --rec125 0.261

  python3 159_shared_control_rho.py --self-test
"""
import argparse
import os
import sys

try:
    import numpy as np
    import pandas as pd
except ImportError:
    sys.exit("needs numpy and pandas (conda activate rnaseq_tools)")

FDR_MAX, DPSI_MIN = 0.05, 0.10
# 05's default --genes, the RBM20 programme
GENES = ["TTN", "CAMK2D", "LDB3", "LMO7", "SORBS1", "PDLIM3", "DTNA", "MLIP",
         "TNNT2"]
# the nominal sizes [042] names and 05's docstring uses
NOMINAL = (166, 28, 166)
REC_TOL = 5e-4
NEED = {"key", "gene", "dPSI_adj", "FDR", "n_g1", "n_g2"}


def rho0(n_d, n_h, n_nf):
    """05's formula, unaltered. Arrays or scalars."""
    n_d = np.asarray(n_d, float)
    n_h = np.asarray(n_h, float)
    n_nf = np.asarray(n_nf, float)
    return ((1 / n_nf) /
            np.sqrt((1 / n_d + 1 / n_nf) * (1 / n_h + 1 / n_nf)))


def load(path, what):
    if not os.path.exists(path):
        sys.exit(f"REFUSING: {what} table not found: {path}")
    t = pd.read_csv(path, sep="\t", low_memory=False)
    miss = NEED - set(t.columns)
    if miss:
        sys.exit(f"REFUSING: {path} is missing {sorted(miss)} — both tables "
                 f"must come from phase2_diffsplice_v2.py")
    return t


def arms(dcm_path, hcm_path):
    """05's pipeline, no collapse, merge on the FULL key."""
    d = load(dcm_path, "--dcm")
    h = load(hcm_path, "--hcm")
    dsig = d[(d.FDR < FDR_MAX) & (d.dPSI_adj.abs() >= DPSI_MIN)].copy()

    def merged(sub):
        m = sub.merge(h, on="key", suffixes=("_DCM", "_HCM"), how="left")
        return m[m.dPSI_adj_HCM.notna()].copy()

    out = {}
    # 1. 05's own set: the nine programme genes
    out["nine-gene (05's own)"] = merged(dsig[dsig.gene.isin(GENES)])
    # 2. the whole significant set — the natural reading of [042]
    out["all significant"] = merged(dsig)
    # 3. every tested event — the baseline as a property of the design
    out["all tested"] = merged(d)
    return out, len(dsig), len(d)


def describe(m):
    if not len(m):
        return None
    r = rho0(m.n_g1_DCM.values, m.n_g1_HCM.values, m.n_g2_DCM.values)
    return dict(n=len(m), med=float(np.median(r)), mean=float(np.mean(r)),
                n_d=float(np.median(m.n_g1_DCM)),
                n_h=float(np.median(m.n_g1_HCM)),
                n_nf=float(np.median(m.n_g2_DCM)))


def report(label, dcm, hcm, rec):
    print(f"\n{'=' * 76}\n{label}\n{'=' * 76}")
    sets, nsig, ntested = arms(dcm, hcm)
    print(f"  --dcm {dcm}")
    print(f"  --hcm {hcm}")
    print(f"  significant rows (NO collapse, 05's convention): {nsig:,}"
          f"   of {ntested:,} tested")
    print(f"\n  {'set':<24}{'n':>7}{'median rho':>12}{'mean':>9}"
          f"{'n_DCM':>8}{'n_HCM':>7}{'n_NF':>7}")
    print("  " + "-" * 74)
    got = {}
    for name, m in sets.items():
        s = describe(m)
        if s is None:
            print(f"  {name:<24}{'0':>7}   (no events)")
            continue
        got[name] = s
        print(f"  {name:<24}{s['n']:>7}{s['med']:>12.4f}{s['mean']:>9.4f}"
              f"{s['n_d']:>8.0f}{s['n_h']:>7.0f}{s['n_nf']:>7.0f}")
    an = float(rho0(*NOMINAL))
    print(f"  {'analytic, nominal n':<24}{'—':>7}{an:>12.4f}"
          f"{'':>9}{NOMINAL[0]:>8}{NOMINAL[1]:>7}{NOMINAL[2]:>7}")

    nine = got.get("nine-gene (05's own)")
    if rec is not None and nine is not None:
        dev = abs(nine["med"] - rec)
        if dev > REC_TOL:
            print(f"\n  ✗ CALIBRATION FAILED: the nine-gene median is "
                  f"{nine['med']:.4f}, but 05 reported {rec:.3f} "
                  f"(deviation {dev:.2g} > {REC_TOL:g}).")
            print("    159 is not computing 05's quantity. Most likely cause:")
            print("    a collapse was applied, or the merge used the collapsed")
            print("    key instead of the full one. No number below is usable.")
            return None
        print(f"\n  ✓ calibration: nine-gene median {nine['med']:.4f} "
              f"reproduces 05's {rec:.3f} (deviation {dev:.2g})")
    return got


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dcm50", default="SE_v2filtered.all_events.tsv")
    ap.add_argument("--hcm50", default="SE_HCM_v2.all_events.tsv")
    ap.add_argument("--dcm125", default="SE_v2filtered_rl125.all_events.tsv")
    ap.add_argument("--hcm125", default="SE_HCMvsNF_rl125.all_events.tsv")
    ap.add_argument("--rec50", type=float, default=0.265,
                    help="rho_null 05 printed at 50 bp")
    ap.add_argument("--rec125", type=float, default=0.261)
    ap.add_argument("--self-test", action="store_true")
    a = ap.parse_args()
    if a.self_test:
        return self_test()
    if self_test():
        print("REFUSING: 159's own self-test failed.")
        return 3

    g50 = report("50 bp", a.dcm50, a.hcm50, a.rec50)
    if g50 is None:
        return 1
    g125 = report("125 bp", a.dcm125, a.hcm125, a.rec125)
    if g125 is None:
        return 1

    print(f"\n{'=' * 76}\nA8a — the three readings of \"rho = 0.27\"\n{'=' * 76}")
    print(f"  {'reading':<26}{'50 bp':>10}{'125 bp':>10}   rounds to")
    for name in ("nine-gene (05's own)", "all significant", "all tested"):
        if name in g50 and name in g125:
            a_, b_ = g50[name]["med"], g125[name]["med"]
            print(f"  {name:<26}{a_:>10.4f}{b_:>10.4f}   "
                  f"{a_:.2f} / {b_:.2f}")
    an = float(rho0(*NOMINAL))
    print(f"  {'analytic, nominal n':<26}{an:>10.4f}{an:>10.4f}   "
          f"{an:.2f} / {an:.2f}   (invariant)")
    print("\n  WHICH ONE THE SENTENCE MEANS IS NOT A MEASUREMENT.")
    print("  * 'For MAGNet's 166 shared controls' names the NOMINAL n, and")
    print("    only the analytic reading uses 166 — the measured median NF")
    print("    is 152 at 50 bp and 148 at 125 bp. On that reading 0.27 is")
    print("    correct and READ-LENGTH INVARIANT: no edit.")
    print("  * 'computed from the per-event group sizes' describes 05's")
    print("    median, which is 0.265 at 50 bp — i.e. the published 0.27")
    print("    was never 05's output. On that reading the sentence needs")
    print("    both a new value AND the 166 corrected to the per-event n.")
    print("\n  These are different claims, not different roundings. Decide")
    print("  which the Methods means before editing; do not pick the one")
    print("  that needs no edit because it needs no edit.")
    return 0


def self_test():
    fails = []

    # the formula, against a value computed by hand from [042]'s nominal n
    want = (1 / 166) / ((1 / 166 + 1 / 166) * (1 / 28 + 1 / 166)) ** 0.5
    got = float(rho0(166, 28, 166))
    if abs(got - want) > 1e-12:
        fails.append(f"rho0(166,28,166) = {got}, hand value {want}")
    if abs(got - 0.2686) > 5e-4:
        fails.append(f"rho0 at nominal n is {got:.4f}; 05's docstring says "
                     f"'about +0.27', so ~0.2686 is expected")

    # independent controls must give exactly 0
    if abs(float(rho0(166, 28, 1e12))) > 1e-9:
        fails.append("an effectively infinite control group did not give 0")
    # The baseline FALLS as the shared control group grows, and RISES as the
    # disease arms grow. (An earlier version of this test asserted the
    # opposite, reasoning that "more shared controls means more sharing".
    # The algebra says otherwise: cov = sigma^2/n_C while
    # var_i = sigma^2(1/n_Di + 1/n_C), so growing n_C shrinks the covariance
    # faster than the variances and rho -> 0. The practical consequence is
    # that the baseline is LARGEST when controls are few relative to cases,
    # which is exactly why the 28-heart HCM arm drives it to 0.27.)
    if not float(rho0(166, 28, 300)) < float(rho0(166, 28, 100)):
        fails.append("rho0 does not fall as the shared control group grows")
    if not float(rho0(1000, 28, 166)) > float(rho0(100, 28, 166)):
        fails.append("rho0 does not rise as the disease arm grows")
    # the small arm dominates: halving n_HCM moves rho more than halving n_DCM
    base = float(rho0(166, 28, 166))
    if not abs(float(rho0(166, 14, 166)) - base) > \
            abs(float(rho0(83, 28, 166)) - base):
        fails.append("the smaller arm is not the dominant term")
    # and must be symmetric in the two disease arms
    if abs(float(rho0(30, 200, 166)) - float(rho0(200, 30, 166))) > 1e-12:
        fails.append("rho0 is not symmetric in n_DCM and n_HCM")

    # median is not mean: a skewed per-event distribution must show the gap,
    # because 05 reports the MEDIAN and quoting a mean would be a different
    # number wearing the same name
    r = rho0([166] * 9 + [20], [28] * 9 + [5], [166] * 10)
    if abs(float(np.median(r)) - float(np.mean(r))) < 1e-6:
        fails.append("the skew fixture does not separate median from mean")

    # 05's convention: NO collapse, merge on the FULL key. A fixture with two
    # rows sharing a collapsed key but differing in flanking coordinates must
    # yield TWO events, not one.
    d = pd.DataFrame({
        "key": ["c:1:2:3:4:5:6", "c:1:2:9:9:9:9"], "gene": ["TTN", "TTN"],
        "dPSI_adj": [0.3, 0.3], "FDR": [0.01, 0.01],
        "n_g1": [166, 166], "n_g2": [166, 166]})
    h = d.copy()
    h["dPSI_adj"] = [0.1, 0.1]
    h["n_g1"] = [28, 28]
    import tempfile
    with tempfile.TemporaryDirectory() as td:
        pd_, ph = os.path.join(td, "d.tsv"), os.path.join(td, "h.tsv")
        d.to_csv(pd_, sep="\t", index=False)
        h.to_csv(ph, sep="\t", index=False)
        sets, nsig, ntested = arms(pd_, ph)
        if nsig != 2:
            fails.append(f"05's uncollapsed convention gave {nsig} significant "
                         f"rows, want 2 — a collapse crept in")
        if len(sets["nine-gene (05's own)"]) != 2:
            fails.append("the full-key merge did not keep both flanking "
                         "variants")
        s = describe(sets["nine-gene (05's own)"])
        if abs(s["med"] - 0.2686) > 5e-4:
            fails.append(f"end-to-end median {s['med']:.4f}, want ~0.2686")
        # a gene outside the nine must drop from 05's set but stay in the others
        d2 = pd.concat([d, d.assign(key=["x:1:1:1:1:1:1", "y:1:1:1:1:1:1"],
                                    gene=["NOTAGENE", "NOTAGENE"])])
        d2.to_csv(pd_, sep="\t", index=False)
        h2 = pd.concat([h, h.assign(key=["x:1:1:1:1:1:1", "y:1:1:1:1:1:1"],
                                    gene=["NOTAGENE", "NOTAGENE"])])
        h2.to_csv(ph, sep="\t", index=False)
        sets2, _, _ = arms(pd_, ph)
        if len(sets2["nine-gene (05's own)"]) != 2:
            fails.append("a non-programme gene leaked into 05's nine-gene set")
        if len(sets2["all significant"]) != 4:
            fails.append("the all-significant set did not widen to 4")

    # the calibration gate must fire on a wrong record and pass on the right
    with tempfile.TemporaryDirectory() as td:
        pd_, ph = os.path.join(td, "d.tsv"), os.path.join(td, "h.tsv")
        d.to_csv(pd_, sep="\t", index=False)
        h.to_csv(ph, sep="\t", index=False)
        import io
        import contextlib
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            bad = report("t", pd_, ph, 0.900)
        if bad is not None or "CALIBRATION FAILED" not in buf.getvalue():
            fails.append("the calibration gate did not fire on a wrong record")
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            ok = report("t", pd_, ph, 0.2686)
        if ok is None:
            fails.append("the calibration gate fired on the correct record")

    print("159 self-test")
    for f in fails:
        print("  ✗ " + f)
    if not fails:
        print("  ✓ 05's formula reproduced by hand at nominal n (0.2686),")
        print("    zero for independent controls, FALLING in the shared n and")
        print("    rising in the disease arms, smaller arm dominant,")
        print("    symmetric in the two arms, median separated from mean,")
        print("    NO collapse and a full-key merge end to end, gene subset")
        print("    honoured, and the calibration gate fires both ways")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
