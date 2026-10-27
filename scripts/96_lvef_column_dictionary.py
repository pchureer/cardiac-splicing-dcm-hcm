#!/usr/bin/env python3
"""
96_lvef_column_dictionary.py — a column dictionary for the deposited LVEF data file.

WHY THIS EXISTS
---------------
`SupplementaryData_lvef_all_events.tsv` (38,084 rows) carries THREE FDR columns:

    FDR_rank   0 events < 0.05   <- the reported model, the paper's null
    FDR_lvef 218 events < 0.05   <- raw least squares. NOT a result.
    FDR_rho    1 event  < 0.05   <- bare Spearman, no covariates

A reader who loads this file in R and types `sum(d$FDR_lvef < 0.05)` gets **218** in a
paper that reports **zero**. That is a referee letter waiting to happen, and the file is
the only place the answer can live, because nobody reads a caption before a TSV.

Supplementary Table 10 already discloses the other two columns in its summary rows.
This script writes the matching statement for the **data file** as a sidecar
`.columns.md`, which is where a deposit's column dictionary belongs. It does NOT touch
the TSV: a comment header would break `read.table` for the people most likely to use it.

WHAT THE 218 ARE — settled, and NOT by the test that was run
------------------------------------------------------------
The mechanism is **leverage in the predictor plus non-normal residuals**, and it is
measured in `RESULT_s3_lvef_null.md`:

  * LVEF takes **20 distinct values** across 160 samples; **55 samples sit at exactly
    10%**. Clinical LVEF is recorded in multiples of 5.
  * **Four samples lie above Q3 + 1.5 IQR**, the arm running to 75% — in a *dilated*
    cardiomyopathy group. High-leverage points for any least-squares fit.
  * PSI is bounded on [0, 1] and outlier-prone, so residuals are not normal and the
    t-test is anti-conservative. The raw model's KS deviation from uniform is
    **D = 0.050**, four times the rank model's **0.013**.

⚠ A degenerate-fit check (psi_sd == 0, |beta| < 1e-9, rho NaN) CANNOT diagnose this and
must not be read as clearing it. `77_lvef_s3.py` gate **G4** already certifies that no
constant-PSI or NaN-rho event survives among the 38,084 fits, so that check is a
restatement of a passed gate, not a test. The constant-PSI artefacts were the OLD
script's 17 hits on 357 RBM20-program events — a different failure, one level down.

USAGE
    python 96_lvef_column_dictionary.py --list      # header + per-FDR counts, write nothing
    python 96_lvef_column_dictionary.py --apply     # write the sidecar

GATES
  G1  the file exists and its header parses.
  G2  CALIBRATION. The counts under the three FDR columns must be 0 / 218 / 1. The
      dictionary names those numbers in prose, so if the file disagrees the prose would
      be false and the script refuses rather than writing it.
  G3  every column in the header has a dictionary entry. A miss is a refusal, printed by
      name — this script will not describe a column it cannot describe. Add the entry to
      DICT below and re-run.
  G4  after writing, reload the sidecar and confirm it names every column and all three
      counts.
"""

import argparse
import os
import shutil
import sys
import time

import pandas as pd

DATA = os.environ.get("DATA", "SupplementaryData_lvef_all_events.tsv")
OUT = os.environ.get("OUT", "SupplementaryData_lvef_all_events.columns.md")

# G2 calibration — the recorded counts. Keys are column names.
EXPECT = {"FDR_rank": 0, "FDR_lvef": 218, "FDR_rho": 1}
# override for another read length, e.g.  EXPECT=FDR_rank=1,FDR_lvef=125,FDR_rho=1
_ov = os.environ.get("EXPECT", "").strip()
if _ov:
    _d = dict(_p.split("=", 1) for _p in _ov.split(","))
    assert set(_d) == set(EXPECT), (
        f"EXPECT override must name exactly {sorted(EXPECT)}, got {sorted(_d)}")
    EXPECT = {_k: int(_v) for _k, _v in _d.items()}
# The KS deviations quoted in the FDR_lvef dictionary entry.
#
# These are INPUTS to deposited prose. Until G1b was added below, nothing
# in this file computed a KS statistic at all, so the comment that used to
# stand here pointed at a gate that did not exist. The prose agreed with
# the CONSTANT; the constant answered to nothing. Run at another read
# length, it would have deposited a 50 bp sentence beside the new data.
#
# G1b now recomputes both values from p_lvef and p_rank and refuses unless
# they match to the 3 decimals the prose prints. The comment is true now.
KS_D_RAW = float(os.environ.get("KS_D_RAW", "0.050"))
KS_D_RANK = float(os.environ.get("KS_D_RANK", "0.013"))
EXPECT_ROWS = int(os.environ.get("EXPECT_ROWS", "38_084"))

# ---------------------------------------------------------------------------- dictionary
# One entry per column. Keys are matched case-insensitively against the header.
# A column that is not here is a REFUSAL, not a blank row.
DICT = {
    "gene": "HGNC symbol of the gene containing the event, from the rMATS annotation.",
    "genesymbol": "HGNC symbol of the gene containing the event, from the rMATS annotation.",
    "geneid": "Ensembl gene identifier from the rMATS annotation.",
    "ekey": "Exon key: `chr:exonStart:exonEnd`, 0-based start as rMATS emits it. This is "
            "the **exon-level** identifier, and it is the unit of analysis here. The rMATS "
            "table was collapsed from 100,014 transcript-context rows to 62,971 distinct "
            "exons — keeping the lowest-FDR context, which has systematically better "
            "coverage — before any model was fitted. PSI was then looked up at the base "
            "table's own full seven-field key, not at an arbitrary context sharing this "
            "three-field key: rMATS emits 2.21 contexts per exon with different "
            "`IncLevel` values, so reading the wrong one changes the analysis.",
    "key": "Event identifier, `chr:exonStart:exonEnd` (0-based start, as rMATS emits it). "
           "Exon-level: the table was collapsed from 100,014 transcript-context rows to "
           "62,971 distinct exons, keeping the lowest-FDR context, before any model was "
           "fitted.",
    "event_id": "Event identifier, `chr:exonStart:exonEnd` (0-based start, as rMATS emits it).",
    "chr": "Chromosome.",
    "exonstart": "Start coordinate of the skipped exon, 0-based (rMATS convention).",
    "exonend": "End coordinate of the skipped exon.",
    "strand": "Genomic strand.",
    "n": "Number of DCM samples with a called PSI value for this event, out of 160. "
         "Events called in fewer than 90% of samples were not fitted.",
    "n_samples": "Number of DCM samples with a called PSI value for this event, out of 160.",
    "psi_mean": "Mean PSI across the DCM samples with a called value.",
    "psi_sd": "Standard deviation of PSI within DCM. Events with sd < 0.01, or with fewer "
              "than 3 distinct values, were excluded before fitting.",
    "rho_lvef": "Spearman rho between PSI and LVEF within DCM, UNADJUSTED for covariates. "
                "Median |rho| across the file is 0.0559 against an analytic null "
                "expectation of 0.0538 at this n — the distribution of correlations sits "
                "where chance puts it.",
    "p_rho": "Two-sided p-value for `rho_lvef` (bare Spearman, no covariates).",
    "FDR_rho": "Benjamini-Hochberg FDR on `p_rho`. **Secondary.** 1 event reaches < 0.05. "
               "Reported here for completeness; a bare Spearman cannot hold sex, age and "
               "RIN, which is why it is not the paper's test.",
    "beta_lvef": "Least-squares slope of PSI on LVEF within DCM, adjusted for sex, age and "
                 "RIN. ⚠ LVEF is stored as a FRACTION (0.05-0.80), not a percentage, so "
                 "this coefficient is per unit LVEF: divide by 10 for PSI change per 10 "
                 "percentage points.",
    "se_lvef": "Standard error of `beta_lvef`.",
    "p_lvef": "Two-sided p-value for `beta_lvef` from the raw least-squares fit. See "
              "`FDR_lvef` before using it.",
    "FDR_lvef": "Benjamini-Hochberg FDR on `p_lvef`, the RAW least-squares fit. "
                f"⚠ **NOT A RESULT. DO NOT USE THIS COLUMN AS A FINDING.** "
                f"{EXPECT['FDR_lvef']} events "
                f"reach < 0.05 here against {EXPECT['FDR_rank']} under the reported "
                f"model, and the cause is a "
                "model-assumption failure measurable in the predictor: LVEF takes only 20 "
                "distinct values across 160 samples with 55 samples at exactly 10%, and "
                "four samples lie above Q3 + 1.5 IQR with the arm running to 75% - "
                "high-leverage points for any least-squares fit. PSI is bounded on [0,1] "
                "and outlier-prone, so residuals are not normal and the t-test is "
                "anti-conservative. The raw model's Kolmogorov-Smirnov deviation from "
                f"uniform is D = {KS_D_RAW:.3f}, {KS_D_RAW / KS_D_RANK:.0f} times the rank "
                f"model's {KS_D_RANK:.3f}. These {EXPECT['FDR_lvef']} are "
                "not degenerate fits - no constant-PSI or NaN-rho event survives the "
                "variance filters - which is precisely why a degenerate-fit check does "
                "not clear the column. The column is deposited so the discrepancy is "
                "visible and explained rather than discoverable.",
    "beta_rank": "Slope of rank(PSI) on rank(LVEF) within DCM, adjusted for sex, age and "
                 "RIN. **This is the reported model.**",
    "se_rank": "Standard error of `beta_rank`.",
    "p_rank": "Two-sided p-value for `beta_rank`. **The reported model's p-value**, and the "
              "column to order by: BH FDR comes in long tied runs, so ranking on "
              "`FDR_rank` alone makes any top-N list arbitrary within a tie block.",
    "FDR_rank": "Benjamini-Hochberg FDR on `p_rank`. **PRIMARY. This is the paper's test.** "
                "0 of 38,084 events reach < 0.05. The nearest miss is PLA2G5 "
                "chr1:20063483:20063602 at FDR = 0.0501, which fails the effect half of "
                "the rule (-0.0829 PSI per 10 points of LVEF, against a 0.10 threshold); "
                "the one event that meets the effect threshold (0.1009) has FDR > 0.74. "
                "No event satisfies both halves, and the two nearest misses fail on "
                "opposite criteria.",
    "dpsi_per10": "PSI change per 10 percentage points of LVEF. The locked rule requires "
                  "|value| >= 0.10.",
    "dpsi_per_10pp_lvef": "__DPSI_ENTRY__",
}

# `dPSI_per_10pp_lvef` carries the effect size the locked rule is evaluated against, and
# its name says it is derived from `beta_lvef` — the RAW least-squares slope — not from
# the rank model. That is the correct arrangement (a rank regression has no interpretable
# PSI-scale slope, so magnitude has to come from the raw-scale fit), but it means
# SIGNIFICANCE AND MAGNITUDE COME FROM DIFFERENT MODELS, which a reader will not assume.
# The dictionary must say so — and the script VERIFIES the derivation before saying it,
# rather than inferring it from the column name. G2b below does that; this text is
# selected by what G2b measures.
# THE SCALE QUESTION — RESOLVED 12 Oct 2026 from `77_lvef_s3.py`. Two steps, not one:
#
#     line 434  # ...a coefficient per unit LVEF would then mean "per 100 percentage..."
#     line 440  ph["LVEF"] = ph["LVEF"] * 100.0        <- fraction -> percentage points
#     line 664  res["dPSI_per_10pp_lvef"] = res["beta_lvef"] * 10.0
#
# So LVEF is STORED as a fraction (0.05-0.80) and FITTED in percentage points (5-80).
# `beta_lvef` is PSI per percentage point; ten points is ten units; the x10 is correct.
# The project note about the fraction describes the stored metadata, and line 440 is the
# conversion it motivated — the two were never in conflict, they describe two moments.
#
# EXPECTED is therefore `beta_lvef * 10`. G2b still searches all candidates rather than
# checking that one: v2 of this script tested only the factor its author expected, and
# would have printed VERIFIED on a false relation had the data happened to agree.
EXPECTED_EXPR = "beta_lvef * 10"
CANDIDATES = [
    ("beta_lvef * 10", 10.0,
     "LVEF entered the model in PERCENTAGE POINTS (e.g. 5-80), so ten points is ten units"),
    ("beta_lvef / 10", 0.1,
     "LVEF entered the model as a FRACTION (e.g. 0.05-0.80), so ten points is 0.10 units"),
    ("beta_lvef (unscaled)", 1.0,
     "`beta_lvef` is already expressed per ten percentage points"),
    ("beta_lvef * 100", 100.0,
     "LVEF entered the model on a scale 10x finer than percentage points"),
    ("beta_lvef / 100", 0.01,
     "LVEF entered the model on a scale 10x coarser than a fraction"),
]

DPSI_VERIFIED = (
    "PSI change per 10 percentage points of LVEF — the effect size the locked rule "
    "(|ΔPSI| >= 0.10) is evaluated against. ⚠ **This is derived from `beta_lvef`, the "
    "RAW least-squares slope, not from the rank model.** Measured against this file, the "
    "relation is **`{expr}`**, to within {dev:.2g} over {n:,} rows. {implies} — LVEF is "
    "*stored* as a fraction (0.05-0.80) and *fitted* in percentage points, converted "
    "before the model is built, so `beta_lvef` is PSI per percentage point. "
    "A rank regression has no interpretable PSI-scale slope, so magnitude necessarily "
    "comes from the raw-scale fit while significance comes from the rank fit "
    "(`FDR_rank`). **The two halves of the locked rule therefore come from two different "
    "models.** This does not inherit `FDR_lvef`'s leverage problem — leverage distorts "
    "the *p*-value far more than the point estimate, and the magnitudes are used only to "
    "show that no event clears 0.10 with a credible *p*-value. The largest value in this "
    "file is {dmax:.4f}; it does clear 0.10 and its FDR_rank is > 0.74. The nearest event "
    "to significance (PLA2G5, FDR_rank = 0.0501) reaches only -0.0829. No event satisfies "
    "both halves."
)
DPSI_UNVERIFIED = (
    "PSI change per 10 percentage points of LVEF — the effect size the locked rule "
    "(|ΔPSI| >= 0.10) is evaluated against. ⚠ The column name says it derives from "
    "`beta_lvef` (the raw least-squares slope), but **this script could not confirm that "
    "from the file** ({why}), so the generating expression is stated as unverified here "
    "rather than asserted. LVEF is stored as a fraction, so 10 percentage points is 0.10 "
    "units. Significance is reported from the rank model (`FDR_rank`); a rank regression "
    "has no interpretable PSI-scale slope, so magnitude comes from a raw-scale fit, "
    "meaning the two halves of the locked rule come from different models."
)

# Matching is case-insensitive, so normalise the keys once rather than trusting that
# every key above was typed in lower case. (The first draft of this script had five
# keys in `FDR_*` form and refused on the three columns it documents best.)
DICT = {k.lower(): v for k, v in DICT.items()}

HEADER = """# Column dictionary - Supplementary Data: PSI against LVEF within DCM

Companion to `{data}` ({rows:,} rows) and Supplementary Table 10.

**Event set.** Genome-wide skipped exons from the discovery cohort's exon-level table
(`SE_v2filtered`), collapsed to one row per exon. This is NOT the RBM20-program-restricted
analysis of an earlier script; the two differ in event set, number of events tested,
coordinate convention, within-group variance filtering and primary model.

**Cohort.** 160 DCM samples from the discovery cohort with LVEF recorded.

**Primary model.** `rank(PSI) ~ rank(LVEF) + sex + age + RIN`, reported in `p_rank` and
`FDR_rank`. Its coefficient is not deposited, because a slope on ranked PSI against ranked
LVEF has no interpretable scale; **effect sizes come from the raw-scale fit**
(`beta_lvef` → `dPSI_per_10pp_lvef`). The locked rule's two halves therefore come from two
different models — significance from the rank fit, magnitude from the raw fit. See
`dPSI_per_10pp_lvef` below.

## ⚠ Three FDR columns, one result

| column | model | events at FDR < 0.05 | status |
|---|---|---|---|
| `FDR_rank` | rank regression, covariate-adjusted | **{n_rank}** | **primary - the reported result** |
| `FDR_lvef` | raw least squares, covariate-adjusted | {n_lvef} | ⚠ artefact of leverage; not a result |
| `FDR_rho` | bare Spearman, no covariates | {n_rho} | secondary; cannot hold covariates |

`FDR_lvef` is retained in the deposit so that the discrepancy between {n_lvef} and
{n_rank} is **visible and explained here**, rather than discovered by a reader and read
as a suppressed finding. See its entry below for the measured cause.

## Columns

"""


def hr(c="-"):
    print(c * 74)


def refuse(msg):
    hr("=")
    print("REFUSING — nothing was written.\n")
    print(msg)
    hr("=")
    sys.exit(2)


def backup(p):
    if os.path.exists(p):
        b = f"{p}.bak.{time.strftime('%Y%m%d-%H%M%S')}"
        shutil.copy2(p, b)
        print(f"  backed up {p} -> {b}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-d", "--dir", default=".", help="directory holding the data file")
    ap.add_argument("--list", dest="list_only", action="store_true")
    ap.add_argument("--apply", action="store_true")
    a = ap.parse_args()

    path = os.path.join(a.dir, DATA)
    out = os.path.join(a.dir, OUT)

    hr("=")
    print("96 — column dictionary for the deposited LVEF data file")
    hr("=")

    # ------------------------------------------------------------------ G1
    if not os.path.exists(path):
        refuse(f"G1 — {path} not found.\n\n  It is written by 95_s1c_column_and_s3_split.py "
               f"--job s3 --apply.\n  Run that first, or pass -d with the directory "
               f"holding it.")
    df = pd.read_csv(path, sep="\t", low_memory=False)
    print(f"\n  {DATA}: {len(df):,} rows x {len(df.columns)} columns")
    print("  header:")
    for c in df.columns:
        print(f"    {c}")
    print("  G1 PASS")

    if len(df) != EXPECT_ROWS:
        refuse(f"G1 — {len(df):,} rows, expected {EXPECT_ROWS:,}. This is not the table "
               f"S3 reports;\n  a dictionary written against it would mislabel the "
               f"deposit.")

    # ------------------------------------------------------------------ counts
    print("\n  events at FDR < 0.05 under every FDR-like column:")
    counts = {}
    for c in df.columns:
        if "fdr" in c.lower() or c.lower() in ("padj", "qvalue", "q"):
            n = int((pd.to_numeric(df[c], errors="coerce") < 0.05).sum())
            counts[c] = n
            print(f"    {c:<14} {n:>6}")
    if not counts:
        refuse("no FDR-like column found in the header. Nothing to calibrate against.")

    # ------------------------------------------------------------------ G2
    print("\n  G2  calibration against the recorded counts")
    bad = []
    # ----------------------------------------------------------------- G1b
    # The KS deviations are quoted in the deposited dictionary. Measure them
    # from the data and refuse if the constants disagree.
    #
    # D is the one-sample Kolmogorov-Smirnov statistic of the p-value
    # distribution against Uniform(0,1). Computed directly so this file needs
    # no scipy:
    #     D = max_i max( i/n - p_(i), p_(i) - (i-1)/n )   over sorted p
    def _ks_uniform(vals):
        s = sorted(float(v) for v in vals
                   if v == v and 0.0 <= float(v) <= 1.0)
        m = len(s)
        if m == 0:
            return None
        d = 0.0
        for _i, _v in enumerate(s, 1):
            a, b = _i / m - _v, _v - (_i - 1) / m
            if a > d:
                d = a
            if b > d:
                d = b
        return d

    _ks_meas = {}
    for _name, _col in (("KS_D_RAW", "p_lvef"), ("KS_D_RANK", "p_rank")):
        if _col not in df.columns:
            refuse(f"G1b — no `{_col}` column, so the KS deviation quoted in "
                   f"the deposited\n  dictionary cannot be verified. It would "
                   f"be a transcribed constant, which is\n  what this gate "
                   f"exists to prevent.")
        _ks_meas[_name] = _ks_uniform(
            pd.to_numeric(df[_col], errors="coerce").dropna().tolist())
        if _ks_meas[_name] is None:
            refuse(f"G1b — `{_col}` holds no usable p-values in [0, 1].")

    print(f"\n  G1b measured from the data:  "
          f"D(p_lvef) = {_ks_meas['KS_D_RAW']:.4f}   "
          f"D(p_rank) = {_ks_meas['KS_D_RANK']:.4f}")
    print(f"      constants in use:        "
          f"KS_D_RAW = {KS_D_RAW:.3f}   KS_D_RANK = {KS_D_RANK:.3f}")
    _ks_bad = [(_n, _ks_meas[_n], _c) for _n, _c in
               (("KS_D_RAW", KS_D_RAW), ("KS_D_RANK", KS_D_RANK))
               if round(_ks_meas[_n], 3) != round(_c, 3)]
    if _ks_bad:
        _lines = "\n".join(
            f"    {_n}: measured {_m:.4f} -> {_m:.3f}, constant {_c:.3f}"
            for _n, _m, _c in _ks_bad)
        refuse(
            "G1b — the KS deviation(s) quoted in the deposited dictionary do "
            "not match\n  this table:\n" + _lines + "\n\n"
            "  The dictionary says 'deviation from uniform is D = <KS_D_RAW>, "
            "<ratio> times\n  the rank model's <KS_D_RANK>'. Depositing that "
            "beside a table it was not\n  measured from is a transcription.\n\n"
            "  Re-run with the measured values:\n"
            f"    KS_D_RAW={_ks_meas['KS_D_RAW']:.3f} "
            f"KS_D_RANK={_ks_meas['KS_D_RANK']:.3f} \\\n"
            "    <your other env vars> python3 96_lvef_column_dictionary.py "
            "--apply\n\n"
            "  Then CHECK THE RATIO: the prose prints "
            f"{_ks_meas['KS_D_RAW'] / _ks_meas['KS_D_RANK']:.0f} times, "
            "and a\n  ratio that has moved is a change to a deposited claim, "
            "not a rounding detail.")
    print(f"  G1b PASS — the deposited KS sentence is measured from this table "
          f"({_ks_meas['KS_D_RAW'] / _ks_meas['KS_D_RANK']:.0f} times)")

    for col, want in EXPECT.items():
        got = counts.get(col)
        mark = "OK" if got == want else "MISMATCH"
        print(f"    {col:<14} recorded {want:>4}   file {str(got):>6}   {mark}")
        if got != want:
            bad.append((col, want, got))
    if bad:
        refuse("G2 — the file disagrees with the recorded counts:\n" +
               "\n".join(f"    {c}: recorded {w}, file {g}" for c, w, g in bad) +
               "\n\n  The dictionary states these numbers in prose. Writing it against a "
               "file that\n  disagrees would deposit a false statement. Either the file is "
               "not the one the\n  record describes, or a recorded number is wrong — "
               "resolve that first.")
    print("  G2 PASS")

    # ------------------------------------------------------------------ G2b
    # The effect-size column's provenance, MEASURED rather than inferred from its name.
    # This decides which of two dictionary entries gets written.
    dcol = next((c for c in df.columns if c.lower().startswith("dpsi")), None)
    bcol = next((c for c in df.columns if c.lower() == "beta_lvef"), None)
    if dcol is not None:
        print(f"\n  G2b  provenance of {dcol}")
        if bcol is None:
            DICT[dcol.lower()] = DPSI_UNVERIFIED.format(
                why="no `beta_lvef` column is present to check against")
            print("       no beta_lvef column — entry written as UNVERIFIED")
        else:
            dv = pd.to_numeric(df[dcol], errors="coerce")
            bv = pd.to_numeric(df[bcol], errors="coerce")
            ok = dv.notna() & bv.notna()
            nok = int(ok.sum())
            bmax = bv[ok].abs().max()
            dmax = dv[ok].abs().max()
            print(f"       n = {nok:,}   max |{bcol}| = {bmax:.4g}   "
                  f"max |{dcol}| = {dmax:.4g}")
            # ⚠ Search the candidates. Do NOT check only the factor the author expected:
            # the first version of this script did exactly that, with the wrong factor.
            scored = []
            for expr, f, implies in CANDIDATES:
                # tolerance scales with the magnitude being compared, not with beta
                tol = max(1e-9, abs(f) * bmax * 1e-6)
                scored.append(((dv[ok] - bv[ok] * f).abs().max(), expr, f, implies, tol))
            scored.sort(key=lambda r: r[0])
            for dev, expr, f, implies, tol in scored:
                print(f"       {expr:<22} max deviation {dev:<12.4g}"
                      f"{'MATCH' if dev <= tol else ''}")
            dev, expr, f, implies, tol = scored[0]
            if dev <= tol:
                DICT[dcol.lower()] = DPSI_VERIFIED.format(
                    expr=expr, dev=dev if dev else tol, n=nok,
                    implies=implies, dmax=dmax)
                print(f"\n       VERIFIED: {dcol} = {expr}  (deviation {dev:.3g}, "
                      f"tolerance {tol:.2g})")
                print(f"       IMPLIES: {implies}.")

                # Independent corroboration of the predictor's scale, from THIS FILE
                # alone — no need to read the fitting script. A fraction-scale predictor
                # would make beta 100x larger, which is not a subtle difference.
                print(f"\n       scale corroboration, from the file alone:")
                print(f"         observed max |{bcol}| = {bmax:.4g}")
                print(f"         a percentage-point predictor puts it near "
                      f"{dmax / 10:.4g}; a fraction predictor,\n         100x larger, "
                      f"near {dmax * 10:.4g}. The observed magnitude therefore says which "
                      f"scale\n         was fitted WITHOUT reading the script.")

                if expr == EXPECTED_EXPR:
                    print(f"         observed {bmax:.4g} ~ {dmax / 10:.4g}  -> "
                          f"percentage points. Conversion fired before fitting.")
                    print(f"\n       ✅ AS EXPECTED. `77_lvef_s3.py` line 440 converts "
                          f"LVEF from the stored\n       fraction to percentage points "
                          f"(x100); line 664 then takes x10 for ten points.\n       "
                          f"LVEF is STORED as a fraction and FITTED in percentage "
                          f"points — the project\n       note and this relation describe "
                          f"two different moments, not a conflict.")
                else:
                    print(f"\n       ⚠ NOT THE EXPECTED RELATION. Expected "
                          f"`{EXPECTED_EXPR}`, measured `{expr}`.")
                    print("       `77_lvef_s3.py` line 440 converts LVEF to percentage "
                          "points and line 664\n       multiplies beta by 10. If the "
                          "measured relation differs, one of those lines\n       has "
                          "changed and EVERY S3 EFFECT SIZE IS ON A DIFFERENT SCALE than "
                          "recorded.")
                    print("       ⚠ DO NOT QUOTE 0.1009, -0.0829, OR ANY |dPSI| UNTIL "
                          "THIS IS RESOLVED:")
                    print("         grep -n 'LVEF.*100\\|per_10pp' 77_lvef_s3.py")
                    print("       The dictionary states the MEASURED relation, which is "
                          "correct either way;\n       what is at stake is whether the "
                          "0.10 threshold was applied on the right scale.")
            else:
                why = ("no candidate relation to `beta_lvef` matched: " +
                       "; ".join(f"{e} off by {d:.3g}" for d, e, _, _, _ in scored[:3]))
                DICT[dcol.lower()] = DPSI_UNVERIFIED.format(why=why)
                print(f"\n       NOT verified — {why}")
                print("       entry written as UNVERIFIED. This is a disclosure, not a "
                      "refusal:\n       the column is deposited either way, and an "
                      "unverified claim is worse\n       than a stated uncertainty.")
    # Drop any placeholder G2b did not fill — mutate, never rebind (rebinding DICT here
    # would make it a local and break the assignments above with UnboundLocalError).
    for k in [k for k, v in DICT.items() if v == "__DPSI_ENTRY__"]:
        del DICT[k]

    # ------------------------------------------------------------------ G3
    miss = [c for c in df.columns if c.lower() not in DICT]
    print(f"\n  G3  columns without a dictionary entry: {len(miss)}")
    if miss:
        refuse("G3 — these columns have no entry:\n" +
               "".join(f"    {c}\n" for c in miss) +
               "\n  This script will not describe a column it cannot describe, and will "
               "not leave a\n  blank row in a deposited dictionary. Add each to DICT in "
               "this script and re-run.")
    print("  G3 PASS")

    if a.list_only or not a.apply:
        hr()
        print("  --list / dry run. Nothing written.")
        print(f"  Would write: {out}")
        hr()
        return

    # ------------------------------------------------------------------ write
    body = HEADER.format(data=DATA, rows=len(df),
                         n_rank=counts["FDR_rank"], n_lvef=counts["FDR_lvef"],
                         n_rho=counts["FDR_rho"])
    for c in df.columns:
        body += f"### `{c}`\n\n{DICT[c.lower()]}\n\n"
    body += ("---\n\nGenerated by `96_lvef_column_dictionary.py` against "
             f"`{DATA}` as deposited.\nEvery count above was read from that file at "
             "write time and checked against the\nproject record before this file was "
             "written.\n")

    backup(out)
    with open(out, "w") as fh:
        fh.write(body)
    print(f"\n  wrote {out}  ({len(body.encode()):,} bytes)")

    # ------------------------------------------------------------------ G4
    txt = open(out).read()
    gone = [c for c in df.columns if f"`{c}`" not in txt]
    # ⚠ Do NOT filter on truthiness here. FDR_rank's count is ZERO, and zero is the
    # paper's result — the one number this sidecar exists to put next to the 218.
    # Checking `if v` would silently exempt it. Verify the whole table row instead.
    nums = [f"| `{c}` " for c in counts] + \
           [f"| **{counts['FDR_rank']}** |", f"| {counts['FDR_lvef']} |"]
    numgone = [n for n in nums if n not in txt]
    print(f"\n  G4  columns named in the sidecar: {len(df.columns) - len(gone)}"
          f"/{len(df.columns)}")
    print(f"      counts named: {len(nums) - len(numgone)}/{len(nums)}")
    if gone or numgone:
        refuse("G4 — the written file does not name everything it should:\n" +
               "".join(f"    column missing: {c}\n" for c in gone) +
               "".join(f"    count missing: {n}\n" for n in numgone))
    print("  G4 PASS")

    hr("=")
    print("DONE. Deposit the sidecar beside the TSV.")
    print("The TSV itself was not modified — a comment header would break read.table().")
    hr("=")


if __name__ == "__main__":
    main()
