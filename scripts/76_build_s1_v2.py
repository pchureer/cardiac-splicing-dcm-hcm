#!/usr/bin/env python3
"""
76_build_s1_v2.py -- Supplementary Table S1, rebuilt. SUPERSEDES 68_build_s1.py.

WHY A REBUILD
-------------
68 built S1 against one composition model. `70_marker_panel_audit.py` then showed that
model's cardiomyocyte score is part cell fraction and part per-cell sarcomeric
upregulation -- four of its nine markers are strongly UP in DCM, and myocyte LOSS predicts
DOWN -- so it absorbs disease signal as well as composition. The fraction-only re-run
(72, 75) measured the cost:

    base, no composition term        926 significant
    nine-marker CM score             899 significant   590 of 926 retained = 63.7%
    fraction-only CM score         1,054 significant   716 of 926 retained = 77.3%

**Adjusting properly for cell fraction finds MORE events than not adjusting at all.** The
two composition models have identical dimension -- two covariates either way, same FB
score -- so nothing differs but which CM score. +13.6 points is the over-adjustment.
→ `RESULT_composition_score_is_a_mediator.md`

WHAT S1 NOW CONTAINS
--------------------
  (a) five event classes x three models: retention, DPSI concordance, and the attrition
      statistics -- recomputed for EACH composition model, not just one
  (b) the covariate ladder, READ from S1b.tsv (built and gated by 69) rather than recomputed
  (c) the per-marker audit: each marker's own disease effect and its keep/drop status
  (d) the both-scores-up observation, which DECISION_title_s3_composition.md places here
      rather than in the Results or the Limitations

THE TRAP THIS SCRIPT INHERITS FROM 68, AND KEEPS FIXED
------------------------------------------------------
**Retention is a SET question; every paired statistic is a ROW question.**

Retention asks "is this exon still called significant?", so each table is collapsed on its
own terms -- that is how each analysis would count. The paired statistics (rho, same-sign,
effect ratio, SE inflation) compare one event measured two ways, so BOTH values must come
from the SAME underlying row.

68's first pass took the lowest-FDR row per exon *independently in each table*. For an
exon in several transcript contexts, the base and adjusted representatives could be
DIFFERENT rows, so every paired statistic compared one transcript context against
another. Retention was unaffected, because it is a set intersection -- which is exactly
why the 590/926 gate passed while every paired number was wrong, in the same direction.
**A gate that passes can sit next to a broken number.**

Here: the base table picks the representative, and the adjusted values are looked up at
**that same full key**.

USAGE
-----
    conda activate rnaseq_tools
    cd /data/project/chureeratp/RNA/HEART

    python 76_build_s1_v2.py --list     # what exists, what is missing, then stop
    python 76_build_s1_v2.py            # build, gate, deposit
"""
import argparse
import os
import sys

# "" = 50 bp tables, "_rl125" = 125 bp. One switch for all three models.
SUF = os.environ.get("TBL_SUFFIX", "")

import numpy as np
import pandas as pd

# ------------------------------------------------------------------ convention
NKEEP = {"SE": 3, "MXE": 5, "A3SS": 5, "A5SS": 5, "RI": 3}
CLASSES = ["SE", "MXE", "A3SS", "A5SS", "RI"]
FDR_CUT, DPSI_CUT = 0.05, 0.10
C_KEY, C_FDR, C_DPSI, C_SE = "key", "FDR", "dPSI_adj", "se_adj"

MODELS = [
    ("base",          "{}_v2filtered" + SUF + ".all_events.tsv", "PSI ~ disease + sex + age + RIN"),
    ("nine_marker",   "{}_v2compadj" + SUF + ".all_events.tsv",
     "+ CM_score (9 markers) + FB_score"),
    ("fraction_only", "{}_v2cmfrac" + SUF + ".all_events.tsv",
     "+ CM_score_frac (5 markers, none upregulated in disease) + FB_score"),
]

# ------------------------------------------------------------------ the record
BASE_SIG = {"SE": 926, "MXE": 231, "A3SS": 90, "A5SS": 84, "RI": 114}
# override for another read length, e.g.
#   BASE_SIG=SE=813,MXE=229,A3SS=87,A5SS=74,RI=91
_ov = os.environ.get("BASE_SIG", "").strip()
if _ov:
    _d = dict(_p.split("=", 1) for _p in _ov.split(","))
    assert set(_d) == set(BASE_SIG), (
        f"BASE_SIG override must name exactly {sorted(BASE_SIG)}, "
        f"got {sorted(_d)}")
    BASE_SIG = {_k: int(_v) for _k, _v in _d.items()}
REC_NINE_SE_RETAINED = int(os.environ.get("REC_NINE_SE_RETAINED", "590"))
REC_FRAC_SE_RETAINED = int(os.environ.get("REC_FRAC_SE_RETAINED", "716"))
REC_FRAC_SE_TOTAL = int(os.environ.get("REC_FRAC_SE_TOTAL", "1054"))
LADDER_FILE = "S1b.tsv"
MARKER_FILE = "marker_panel_audit.tsv"
OUT_A, OUT_C = "S1a.tsv", "S1c_markers.tsv"


def spearman(x, y):
    """Pearson on ranks -- no scipy dependency, matching 68."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 3:
        return np.nan
    rx = pd.Series(x[ok]).rank().to_numpy()
    ry = pd.Series(y[ok]).rank().to_numpy()
    if rx.std() == 0 or ry.std() == 0:
        return np.nan
    return float(np.corrcoef(rx, ry)[0, 1])


def load(path, et):
    """Rows, with ekey attached. None if the file is absent."""
    if not (os.path.exists(path) and os.path.getsize(path)):
        return None
    t = pd.read_csv(path, sep="\t", low_memory=False)
    for c in (C_KEY, C_FDR):
        if c not in t.columns:
            sys.exit(f"[fatal] {path}: no '{c}' column")
    d = C_DPSI if C_DPSI in t.columns else ("dPSI" if "dPSI" in t.columns else None)
    if d is None:
        sys.exit(f"[fatal] {path}: no dPSI column")
    t = t.copy()
    t["_d"] = pd.to_numeric(t[d], errors="coerce")
    t["_f"] = pd.to_numeric(t[C_FDR], errors="coerce")
    t["_se"] = (pd.to_numeric(t[C_SE], errors="coerce")
                if C_SE in t.columns else np.nan)
    nk = NKEEP[et]
    t["ekey"] = [":".join(str(k).split(":")[:nk]) for k in t[C_KEY]]
    return t.dropna(subset=["_f"])


def representatives(t):
    """Lowest-FDR row per collapsed exon -- the project's convention."""
    return t.loc[t.groupby("ekey")["_f"].idxmin()]


def significant(rep):
    return set(rep[(rep._f < FDR_CUT) & (rep._d.abs() >= DPSI_CUT)]["ekey"])


def paired(base_t, adj_t, base_sig_keys):
    """
    Paired statistics, done correctly.

    The base table picks the representative row for each significant exon; the adjusted
    values are looked up at THAT SAME FULL KEY. Never at the adjusted table's own
    independently chosen representative -- see the module docstring.
    """
    brep = representatives(base_t)
    sig = brep[brep.ekey.isin(base_sig_keys)]
    a = adj_t.drop_duplicates(subset=[C_KEY]).set_index(C_KEY)
    j = sig.join(a[["_d", "_f", "_se"]].add_suffix("_a"), on=C_KEY)
    j = j.dropna(subset=["_d_a"])
    return j


def attrition(j, retained_keys):
    """Effect retained, SE inflation, >half lost, direction -- for kept vs lost."""
    out = {}
    kept = j[j.ekey.isin(retained_keys)]
    lost = j[~j.ekey.isin(retained_keys)]
    for lab, sub in (("kept", kept), ("lost", lost)):
        if not len(sub):
            out[f"median_effect_ratio_{lab}"] = np.nan
            out[f"se_inflation_{lab}"] = np.nan
            continue
        r = (sub._d_a.abs() / sub._d.abs()).replace([np.inf, -np.inf], np.nan)
        out[f"median_effect_ratio_{lab}"] = round(float(r.median()), 3)
        s = (sub._se_a / sub._se).replace([np.inf, -np.inf], np.nan)
        out[f"se_inflation_{lab}"] = (round(float(s.median()), 3)
                                      if s.notna().any() else np.nan)
    if len(lost):
        r = (lost._d_a.abs() / lost._d.abs()).replace([np.inf, -np.inf], np.nan)
        out["pct_lost_half"] = round(100.0 * float((r < 0.5).mean()), 1)
        out["pct_lost_same_direction"] = round(
            100.0 * float(((lost._d * lost._d_a) > 0).mean()), 1)
    else:
        out["pct_lost_half"] = np.nan
        out["pct_lost_same_direction"] = np.nan
    return out


# ----------------------------------------------------------------------- --list
def do_list():
    print("\n" + "=" * 78)
    print("S1 INPUTS")
    print("=" * 78)
    print(f"\n  {'class':<7}" + "".join(f"{m:<16}" for m, _, _ in MODELS))
    print("  " + "-" * 56)
    missing = []
    for et in CLASSES:
        row = f"  {et:<7}"
        for m, pat, _ in MODELS:
            p = pat.format(et)
            if os.path.exists(p) and os.path.getsize(p):
                row += f"{'ok':<16}"
            else:
                row += f"{'MISSING':<16}"
                missing.append(p)
        print(row)
    for f, what in ((LADDER_FILE, "covariate ladder, from 69"),
                    (MARKER_FILE, "per-marker audit, from 70")):
        print(f"\n  {f:<28}{'ok' if os.path.exists(f) else 'MISSING'}   ({what})")
    if missing:
        print(f"\n  {len(missing)} table(s) missing:")
        for p in missing[:8]:
            print(f"    * {p}")
        if any("cmfrac" in p for p in missing):
            print("\n  The *_v2cmfrac tables come from 75_cmfrac_allclasses.sbatch")
            print("  (SE came from 72). Wait for the array to finish, then rebuild.")
    else:
        print("\n  All inputs present. Run without --list.\n")


# ------------------------------------------------------------------------ build
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--allow-partial", action="store_true",
                    help="build with whatever fraction-only tables exist, and say so")
    a = ap.parse_args()
    if a.list:
        do_list()
        return

    print("\n" + "=" * 78)
    print("SUPPLEMENTARY TABLE S1 -- composition, the covariate ladder, the panel audit")
    print("=" * 78)
    print(f"\n  convention: exon level, NKEEP per class, lowest-FDR row kept")
    print(f"  rule:       FDR < {FDR_CUT} and |dPSI| >= {DPSI_CUT}")

    rows, fail, warn = [], [], []

    for et in CLASSES:
        base_t = load(MODELS[0][1].format(et), et)
        if base_t is None:
            fail.append(f"{et}: base table missing")
            continue
        brep = representatives(base_t)
        bsig = significant(brep)
        n_tested = len(brep)

        if len(bsig) != BASE_SIG[et]:
            fail.append(f"{et}: base {len(bsig)} significant, recorded {BASE_SIG[et]}")

        rows.append(dict(event_class=et, model="base",
                         model_formula=MODELS[0][2],
                         n_tested_exon=n_tested, n_significant=len(bsig),
                         n_base_significant=len(bsig), n_retained="",
                         pct_retained="", rho_dpsi="", pct_same_sign=""))

        for mname, pat, formula in MODELS[1:]:
            p = pat.format(et)
            adj_t = load(p, et)
            if adj_t is None:
                warn.append(f"{et}: {mname} table absent ({p})")
                continue
            arep = representatives(adj_t)
            asig = significant(arep)
            retained = bsig & asig
            j = paired(base_t, adj_t, bsig)
            r = dict(event_class=et, model=mname, model_formula=formula,
                     n_tested_exon=len(arep), n_significant=len(asig),
                     n_base_significant=len(bsig), n_retained=len(retained),
                     pct_retained=round(100.0 * len(retained) / len(bsig), 1),
                     rho_dpsi=round(spearman(j._d, j._d_a), 3) if len(j) else np.nan,
                     pct_same_sign=(round(100.0 * float(((j._d * j._d_a) > 0).mean()), 3)
                                    if len(j) else np.nan),
                     n_paired=len(j))
            r.update(attrition(j, retained))
            rows.append(r)

            if et == "SE" and mname == "nine_marker" and len(retained) != REC_NINE_SE_RETAINED:
                fail.append(f"SE nine-marker retention {len(retained)}, recorded "
                            f"{REC_NINE_SE_RETAINED}")
            if et == "SE" and mname == "fraction_only":
                if len(retained) != REC_FRAC_SE_RETAINED:
                    fail.append(f"SE fraction-only retention {len(retained)}, recorded "
                                f"{REC_FRAC_SE_RETAINED}")
                if len(asig) != REC_FRAC_SE_TOTAL:
                    fail.append(f"SE fraction-only total {len(asig)}, recorded "
                                f"{REC_FRAC_SE_TOTAL}")

    df = pd.DataFrame(rows)

    # ------------------------------------------------------------ report (a)
    print("\n" + "-" * 78)
    print("  (a) RETENTION AND CONCORDANCE, five classes x three models")
    print("-" * 78)
    print(f"  {'class':<6}{'model':<15}{'tested':>8}{'sig':>7}{'retained':>10}"
          f"{'%':>7}{'rho':>7}{'same sign':>11}")
    print("  " + "-" * 70)
    for _, r in df.iterrows():
        ret = r.n_retained if r.n_retained != "" else "-"
        pct = f"{r.pct_retained}" if r.pct_retained != "" else "-"
        rho = f"{r.rho_dpsi}" if r.rho_dpsi not in ("", None) and pd.notna(r.rho_dpsi) else "-"
        ss = (f"{r.pct_same_sign:.1f}%" if r.pct_same_sign not in ("", None)
              and pd.notna(r.pct_same_sign) else "-")
        print(f"  {r.event_class:<6}{r.model:<15}{r.n_tested_exon:>8}"
              f"{r.n_significant:>7}{str(ret):>10}{pct:>7}{rho:>7}{ss:>11}")

    adj = df[df.model != "base"]
    if len(adj):
        print("\n" + "-" * 78)
        print("  WHY THE LOSS IS NOT CONFOUNDING -- attrition, per model")
        print("-" * 78)
        print(f"  {'class':<6}{'model':<15}{'eff kept':>9}{'eff LOST':>9}"
              f"{'SE kept':>9}{'SE lost':>9}{'>half':>8}{'same dir':>10}")
        print("  " + "-" * 70)
        for _, r in adj.iterrows():
            g = lambda k: (f"{r[k]:.3f}" if pd.notna(r.get(k)) else "-")   # noqa: E731
            h = lambda k: (f"{r[k]:.1f}%" if pd.notna(r.get(k)) else "-")  # noqa: E731
            print(f"  {r.event_class:<6}{r.model:<15}"
                  f"{g('median_effect_ratio_kept'):>9}"
                  f"{g('median_effect_ratio_lost'):>9}"
                  f"{g('se_inflation_kept'):>9}{g('se_inflation_lost'):>9}"
                  f"{h('pct_lost_half'):>8}{h('pct_lost_same_direction'):>10}")
        print("\n  The SE-inflation columns are the argument: flat across kept and lost")
        print("  rules out differential noise. Confounding would shrink the LOST events'")
        print("  effects specifically.")

    # ------------------------------------------------------------ (b) ladder
    print("\n" + "-" * 78)
    print(f"  (b) THE COVARIATE LADDER -- read from {LADDER_FILE}, not recomputed")
    print("-" * 78)
    if not os.path.exists(LADDER_FILE):
        warn.append(f"{LADDER_FILE} absent -- run 69_build_s1b_ladder.py")
        print(f"  MISSING. Run 69_build_s1b_ladder.py.")
    else:
        lad = pd.read_csv(LADDER_FILE, sep="\t")
        print(f"  {'model':<42}{'tested':>9}{'sig':>7}{'% of base':>11}")
        print("  " + "-" * 68)
        for _, r in lad.sort_values("order").iterrows():
            print(f"  {str(r.model)[:41]:<42}{r.n_tested_exon_level:>9}"
                  f"{r.n_significant:>7}{r.pct_of_base:>10.1f}%")
        base_row = lad[lad.order == 1]
        if len(base_row) and int(base_row.iloc[0].n_significant) != BASE_SIG["SE"]:
            fail.append(f"{LADDER_FILE} base row is "
                        f"{int(base_row.iloc[0].n_significant)}, not {BASE_SIG['SE']}")
        sp = int(lad.n_significant.max() - lad.n_significant.min())
        print(f"\n  spread {sp} events = "
              f"{100.0*sp/BASE_SIG['SE']:.1f}% of the headline; every model within 10%")

    # ------------------------------------------------------------ (c) markers
    print("\n" + "-" * 78)
    print(f"  (c) THE MARKER PANEL AUDIT -- read from {MARKER_FILE}")
    print("-" * 78)
    mk = None
    if not os.path.exists(MARKER_FILE):
        warn.append(f"{MARKER_FILE} absent -- run 70_marker_panel_audit.py "
                    f"--gene-matrix magnet_expr.log2cpm.tsv.gz")
        print(f"  MISSING. Run 70_marker_panel_audit.py --gene-matrix "
              f"magnet_expr.log2cpm.tsv.gz")
    else:
        mk = pd.read_csv(MARKER_FILE, sep="\t")
        bcol = "beta_sd_technical"
        if bcol not in mk.columns:
            warn.append(f"{MARKER_FILE} has no {bcol} column")
        else:
            print(f"  {'panel':<18}{'gene':<12}{'disease gene':<14}{'beta SD':>9}"
                  f"{'FDR':>11}")
            print("  " + "-" * 66)
            for _, r in mk.iterrows():
                b = r.get(bcol)
                q = r.get("fdr_technical")
                print(f"  {str(r.get('panel',''))[:17]:<18}{str(r.get('gene',''))[:11]:<12}"
                      f"{str(r.get('disease_gene','-')):<14}"
                      f"{(f'{b:+.3f}' if pd.notna(b) else '-'):>9}"
                      f"{(f'{q:.3g}' if pd.notna(q) else '-'):>11}")
            # (d) the both-scores-up observation
            cm = mk[mk.gene == "CM_score"]
            fb = mk[mk.gene == "FB_score"]
            if len(cm) and len(fb):
                c_, f_ = float(cm.iloc[0][bcol]), float(fb.iloc[0][bcol])
                print(f"\n  (d) BOTH SCORES RISE IN DCM -- cardiomyocyte {c_:+.3f} SD, "
                      f"fibroblast {f_:+.3f} SD,")
                print(f"      while the two anti-correlate (r = -0.222). A simple")
                print(f"      two-compartment shift -- myocyte loss, fibroblast gain --")
                print(f"      predicts cardiomyocyte DOWN. Candidate account, UNTESTED:")
                print(f"      failing explant tissue may carry a greater total")
                print(f"      structural-gene fraction than donor left ventricle, with")
                print(f"      proportionally less endothelium, immune and adipose content.")
                print(f"      Report as an observation; do not explain it away.")
                if c_ <= 0:
                    warn.append("CM_score is no longer positive -- the (d) text assumes "
                                "both scores rise; recheck before quoting it")

    # ------------------------------------------------------------ gates
    print("\n" + "=" * 78)
    if warn:
        print("  ADVISORIES")
        for w in warn:
            print(f"    * {w}")
        print()
    if fail:
        print("S1 NOT BUILT -- a gate failed")
        print("=" * 78)
        for f_ in fail:
            print(f"  * {f_}")
        print("\n  Nothing written. Fix the input or the record, not this script.\n")
        return

    partial = [et for et in CLASSES
               if not os.path.exists(MODELS[2][1].format(et))]
    if partial and not a.allow_partial:
        print("S1 NOT BUILT -- the fraction-only model is incomplete")
        print("=" * 78)
        print(f"  missing for: {', '.join(partial)}")
        print("  S1(a) is a five-class table. Five classes under one composition model")
        print("  and one class under the other is not a table a reviewer can read.")
        print("  Finish 75_cmfrac_allclasses.sbatch, or pass --allow-partial to build")
        print("  anyway and have the gap stated in the output.\n")
        return
    if partial:
        print(f"  [PARTIAL] fraction-only model absent for: {', '.join(partial)}")
        print(f"            Built on request. S1 must say so in its caption.\n")

    df.to_csv(OUT_A, sep="\t", index=False)
    if mk is not None:
        mk.to_csv(OUT_C, sep="\t", index=False)

    print("S1 BUILT")
    print("=" * 78)
    print(f"  {OUT_A:<22}{len(df)} rows   (a) retention, concordance, attrition")
    print(f"  {LADDER_FILE:<22}{'(from 69)':<9}  (b) the covariate ladder")
    if mk is not None:
        print(f"  {OUT_C:<22}{len(mk)} rows   (c) the marker audit, (d) both-scores-up")

    se = df[(df.event_class == 'SE')]
    b = se[se.model == 'base'].iloc[0]
    n9 = se[se.model == 'nine_marker']
    nf = se[se.model == 'fraction_only']
    print(f"\n  THE HEADLINE, skipped exons:")
    print(f"    base, no composition term          {b.n_significant:>5} significant")
    if len(n9):
        print(f"    nine-marker CM score               "
              f"{int(n9.iloc[0].n_significant):>5} significant   "
              f"{int(n9.iloc[0].n_retained)}/{b.n_significant} retained = "
              f"{n9.iloc[0].pct_retained}%")
    if len(nf):
        print(f"    fraction-only CM score             "
              f"{int(nf.iloc[0].n_significant):>5} significant   "
              f"{int(nf.iloc[0].n_retained)}/{b.n_significant} retained = "
              f"{nf.iloc[0].pct_retained}%")
        if int(nf.iloc[0].n_significant) > int(b.n_significant):
            print(f"\n  >>> Adjusting properly for cell fraction finds MORE events than")
            print(f"      not adjusting at all. Both composition models have identical")
            print(f"      dimension, so nothing differs but which CM score. Report the")
            print(f"      fraction-only figures as primary and the nine-marker ones as")
            print(f"      the conservative bound.")
    print()


if __name__ == "__main__":
    main()
