#!/usr/bin/env python
"""
77_lvef_s3.py -- Supplementary S3: does PSI track LVEF within DCM?

WHY THIS REPLACES lvef_correlation.py RATHER THAN PATCHING IT
-------------------------------------------------------------
The old script's seven "significant" events are artefacts, not results:
every one has beta = 0.0000 and rho = nan, which is what a regression on a
CONSTANT PSI vector returns. It also ran on the v1 table, and its own
docstring advertises a --cov argument that does not exist. Three defects, one
of which silently produced the headline, so the honest move is a rebuild with
the filter that should have been there and a gate on a number we already know.

WHAT THIS DOES
--------------
  PSI_sample ~ LVEF + sex + age + RIN,  within one disease group (default DCM)

per skipped-exon event, with a PSI-VARIANCE FILTER applied BEFORE any model is
fitted. The expected result is a null. A null is reportable; an artefact is not.

THE FILTER, and why each clause exists
--------------------------------------
  --min-sd 0.01      sd(PSI) over the analysed samples. Kills the constant-PSI
                     events outright -- that is the actual bug.
  --min-called 0.90  fraction of analysed samples with a non-missing PSI.
                     rMATS writes "NA" wherever coverage failed; an event called
                     in 12 of 160 hearts is not a measurement of anything.
  --min-distinct 3   number of DISTINCT PSI values. A vector that is 0.0 in 150
                     samples and 1.0 in two has sd > 0.01 and is still junk.

GATES (it refuses rather than guesses)
-------------------------------------
  G0  the key's first three fields ARE chr:exonStart:exonEnd, checked against the
      table's own columns -- if not, the exon-level collapse is collapsing on
      something else and the rMATS join is meaningless
  G1  SE_v2filtered collapses to exactly 62,971 exon-level events
  G2  SE_v2filtered carries exactly 926 significant events under the locked rule
  G3  the rMATS table joins to >= 95% of those exon keys (coordinate convention)
  G4  no retained event has sd == 0, and no retained event returns rho = nan
  G5  every analysed sample has LVEF, sex, age and RIN

The effect column is RESOLVED, not assumed: the phase-2 tables name it `dPSI_adj`
(with `dPSI_raw` beside it), not `dPSI`. v1 of this script hard-coded `dPSI` and was
refused on the real table -- which is the refusal working, but it should not have
needed to. 76_build_s1_v2.py resolves the same way.

Run --list FIRST. It prints the real columns of all four inputs and exits
without computing anything. Every file layout in this project has differed from
some docstring that described it (magnet_covariates.tsv has no header row;
magnet_expr.per_sample.tsv is samples x features), so nothing here is assumed.

USAGE
-----
  cd /data/project/chureeratp/RNA/HEART
  conda activate rnaseq_tools

  python 77_lvef_s3.py --list
  python 77_lvef_s3.py
  python 77_lvef_s3.py --group HCM --out S3_lvef_hcm.tsv

NOTE ON QUANTIFICATION: --mats defaults to rmats_magnet/SE.MATS.JCEC.txt,
because 04_phase2_v2_checks.sbatch:21 built SE_v2filtered from JCEC. Do not
point this at MATS.JC.txt -- it would be a different tested set from the base
table it gates against. -> RESULT_jcec_not_jc.md
"""

import argparse
import os
import re
import sys

import numpy as np
import pandas as pd
from scipy import stats

# ---------------------------------------------------------------- constants

EXPECT_TESTED = int(os.environ.get("EXPECT_TESTED", "62_971"))   # exon-level, SE_v2filtered           (G1)
EXPECT_SIG = int(os.environ.get("EXPECT_SIG", "926"))   # FDR<0.05 and |dPSI|>=0.10           (G2)
MIN_JOIN_FRAC = 0.95            # rMATS -> events-file key overlap    (G3)

NKEEP_SE = 3                    # chr:exonStart:exonEnd
C_KEY = "key"

# The phase-2 tables name the adjusted effect `dPSI_adj`, NOT `dPSI`. 76_build_s1_v2.py
# resolves this; v1 of this script hard-coded `dPSI` and was refused on the real table.
# Order matters: dPSI_adj is the quantity the locked rule and the 926 are defined on.
C_DPSI_CANDIDATES = ("dPSI_adj", "dPSI_raw", "dPSI")

FDR_MAX = 0.05
DPSI_MIN = 0.10

COV_COLS = ["run", "group", "sex", "age", "RIN", "pool"]   # NO header row
RUN_RE = re.compile(r"(SRR\d+|ERR\d+|DRR\d+|SAMEA\d+)")

MATS_COORD = ["chr", "exonStart_0base", "exonEnd",
              "upstreamES", "upstreamEE", "downstreamES", "downstreamEE"]


def die(msg):
    print("\n[REFUSE] " + msg + "\n", file=sys.stderr)
    sys.exit(2)


def ok(msg):
    print("  [ok] " + msg)


# ---------------------------------------------------------------- loaders

def load_covariates(path):
    """magnet_covariates.tsv has NO header row: six positional columns."""
    first = pd.read_csv(path, sep="\t", nrows=1, header=None).iloc[0, 0]
    has_header = not bool(RUN_RE.fullmatch(str(first).strip()))
    if has_header:
        cov = pd.read_csv(path, sep="\t")
    else:
        cov = pd.read_csv(path, sep="\t", header=None, names=COV_COLS)
    cov["run"] = cov["run"].astype(str).str.strip()
    return cov, has_header


def load_pheno(path, lvef_col=None):
    ph = pd.read_csv(path, sep="\t", dtype=str)
    ph.columns = [c.strip() for c in ph.columns]

    runcol = None
    for c in ph.columns:
        vals = ph[c].dropna().astype(str).head(20)
        if len(vals) and all(RUN_RE.fullmatch(v.strip()) for v in vals):
            runcol = c
            break
    if runcol is None:
        for c in ph.columns:
            if c.lower() in ("run_accession", "run", "sra_run"):
                runcol = c
                break

    if lvef_col is None:
        cands = [c for c in ph.columns if re.search(r"lvef|ejection|\bef\b", c, re.I)]
        lvef_col = cands[0] if len(cands) == 1 else (cands[0] if cands else None)
    return ph, runcol, lvef_col, [c for c in ph.columns
                                  if re.search(r"lvef|ejection|\bef\b", c, re.I)]


def bam_runs(path):
    """rMATS --b1/--b2 files are ONE line of COMMA-separated BAM paths.

    v1 of this script read them line by line and found 1 sample in a 166-sample
    cohort, which G5 caught. Split on commas AND whitespace, in order -- the order IS
    the IncLevel1/IncLevel2 column order and must not be disturbed.
    """
    with open(path) as fh:
        raw = fh.read()
    runs = []
    for tok in (t for t in re.split(r"[,\s]+", raw) if t.strip()):
        m = RUN_RE.search(os.path.basename(tok))
        runs.append(m.group(1) if m else None)
    return runs


def inclevel_widths(path):
    """(n1, n2) from one row of the rMATS table -- cheap, and it pins the bam lists."""
    h = pd.read_csv(path, sep="\t", nrows=1)
    out = []
    for c in ("IncLevel1", "IncLevel2"):
        out.append(len(str(h[c].iloc[0]).split(",")) if c in h.columns else None)
    return tuple(out)


def ekey_from_key(s):
    return s.astype(str).str.split(":").str[:NKEEP_SE].str.join(":")


def resolve_dpsi(cols, path):
    for c in C_DPSI_CANDIDATES:
        if c in cols:
            return c
    die(f"{path} has none of {C_DPSI_CANDIDATES}. Columns: {list(cols)}")


def ekey_from_mats(mats):
    return (mats["chr"].astype(str) + ":" +
            mats["exonStart_0base"].astype(str) + ":" +
            mats["exonEnd"].astype(str))


def fullkey_from_mats(mats):
    """The base table's `key`: all seven coordinate fields, no strand, in rMATS order."""
    return (mats["chr"].astype(str) + ":" +
            mats["exonStart_0base"].astype(str) + ":" +
            mats["exonEnd"].astype(str) + ":" +
            mats["upstreamES"].astype(str) + ":" +
            mats["upstreamEE"].astype(str) + ":" +
            mats["downstreamES"].astype(str) + ":" +
            mats["downstreamEE"].astype(str))


# ---------------------------------------------------------------- listing

def do_list(a):
    print("=" * 74)
    print("77_lvef_s3.py --list :: inputs as they actually are")
    print("=" * 74)

    for p in (a.events, a.mats, a.cov, a.pheno, a.b1, a.b2):
        print(f"  {'OK ' if os.path.exists(p) else 'MISSING'}  {p}")
    if not all(os.path.exists(p) for p in (a.events, a.mats, a.cov, a.pheno, a.b1, a.b2)):
        die("one or more inputs missing; nothing else can be reported")

    print("\n-- SE_v2filtered (events) --")
    ev = pd.read_csv(a.events, sep="\t", nrows=3)
    print("  columns:", list(ev.columns))
    print("  first key:", ev[C_KEY].iloc[0] if C_KEY in ev.columns else "<no 'key' column>")
    have = [c for c in C_DPSI_CANDIDATES if c in ev.columns]
    print(f"  effect columns present: {have or 'NONE -- this will refuse'}"
          f"{'  -> using ' + have[0] if have else ''}")
    print("  gene names available:",
          "yes" if "gene" in ev.columns else "no (output keyed by ekey only)")
    print("  own sd_psi column:",
          "yes (artefact census will run)" if "sd_psi" in ev.columns else "no")

    print("\n-- rMATS table --")
    mt = pd.read_csv(a.mats, sep="\t", nrows=3)
    print("  columns:", list(mt.columns))
    missing = [c for c in MATS_COORD if c not in mt.columns]
    print("  coordinate columns missing:", missing if missing else "none")
    for c in ("IncLevel1", "IncLevel2"):
        if c in mt.columns:
            n = len(str(mt[c].iloc[0]).split(","))
            print(f"  {c}: {n} per-sample values")

    print("\n-- covariates --")
    cov, has_header = load_covariates(a.cov)
    print(f"  header row present: {has_header}")
    print("  columns:", list(cov.columns))
    print("  groups:", dict(cov["group"].value_counts()) if "group" in cov else "n/a")

    print("\n-- pheno --")
    ph, runcol, lvef_col, lvef_cands = load_pheno(a.pheno, a.lvef_col)
    print(f"  rows: {len(ph)}")
    print(f"  run column detected: {runcol}")
    print(f"  LVEF-like columns:  {lvef_cands}")
    print(f"  LVEF column chosen: {lvef_col}")
    if lvef_col and lvef_col in ph.columns:
        v = pd.to_numeric(ph[lvef_col], errors="coerce")
        print(f"  LVEF non-missing: {v.notna().sum()} of {len(v)}; "
              f"range {v.min()}-{v.max()}")

    print("\n-- bam lists --")
    r1, r2 = bam_runs(a.b1), bam_runs(a.b2)
    print(f"  {a.b1}: {len(r1)} entries, {sum(x is None for x in r1)} unparsed")
    print(f"  {a.b2}: {len(r2)} entries, {sum(x is None for x in r2)} unparsed")
    w1, w2 = inclevel_widths(a.mats)
    print(f"  IncLevel widths: {w1} / {w2}")
    if (w1, w2) == (len(r1), len(r2)):
        print("  MATCH -- bam order will align to the PSI columns")
    else:
        print("  ⚠ MISMATCH -- the run will refuse (G6). The bam list is not being "
              "parsed the way rMATS read it.")

    print("\n-- can the join be made? --")
    if runcol is None:
        print("  NO -- no run-accession column found in the pheno file")
    elif lvef_col is None:
        print("  NO -- no LVEF column found; pass --lvef-col")
    else:
        ids = set(ph[runcol].astype(str).str.strip())
        print(f"  pheno runs matching {a.b1}: {sum(r in ids for r in r1 if r)}/{len(r1)}")
        print(f"  pheno runs matching covariates: "
              f"{len(ids & set(cov['run']))}/{len(cov)}")
        print("  YES -- rerun without --list")
    print()


# ---------------------------------------------------------------- stats

def ols_lvef(y, X):
    """Return (beta_lvef, se, p) for column 0 of X (LVEF), intercept appended."""
    A = np.column_stack([X, np.ones(len(y))])
    beta, *_ = np.linalg.lstsq(A, y, rcond=None)
    resid = y - A @ beta
    dof = len(y) - A.shape[1]
    if dof <= 0:
        return np.nan, np.nan, np.nan
    s2 = float(resid @ resid) / dof
    try:
        XtXi = np.linalg.inv(A.T @ A)
    except np.linalg.LinAlgError:
        return np.nan, np.nan, np.nan
    se = float(np.sqrt(s2 * XtXi[0, 0]))
    if not np.isfinite(se) or se == 0:
        return float(beta[0]), se, np.nan
    t = float(beta[0]) / se
    return float(beta[0]), se, float(2 * stats.t.sf(abs(t), dof))


def spearman(x, y):
    if len(np.unique(x)) < 3 or len(np.unique(y)) < 3:
        return np.nan, np.nan
    rx = stats.rankdata(x)
    ry = stats.rankdata(y)
    r = float(np.corrcoef(rx, ry)[0, 1])
    if not np.isfinite(r) or abs(r) >= 1:
        return r, np.nan
    n = len(x)
    t = r * np.sqrt((n - 2) / (1 - r * r))
    return r, float(2 * stats.t.sf(abs(t), n - 2))


def bh(p):
    p = np.asarray(p, float)
    out = np.full(p.shape, np.nan)
    m = np.isfinite(p)
    if not m.any():
        return out
    pv = p[m]
    n = len(pv)
    order = np.argsort(pv)
    q = pv[order] * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    adj = np.empty(n)
    adj[order] = np.clip(q, 0, 1)
    out[m] = adj
    return out


# ---------------------------------------------------------------- main

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--events", default="SE_v2filtered.all_events.tsv")
    ap.add_argument("--mats", default="rmats_magnet/SE.MATS.JCEC.txt")
    ap.add_argument("--cov", default="magnet_covariates.tsv")
    ap.add_argument("--pheno",
                    default="filereport_read_run_PRJNA595151_withpheno.tsv")
    ap.add_argument("--b1", default="dcm_bams.txt")
    ap.add_argument("--b2", default="nf_bams.txt")
    ap.add_argument("--lvef-col", default=None)
    ap.add_argument("--group", default="DCM")
    # Same interface as phase2_diffsplice_v2.py, deliberately: the file and column
    # names the user already has in hand for the composition models.
    ap.add_argument("--extra-cov-file", default=None,
                    help="samples x features TSV, e.g. magnet_expr.cmfrac.tsv")
    ap.add_argument("--extra-cov", default=None,
                    help="comma-separated columns, e.g. CM_score_frac,FB_score")
    ap.add_argument("--min-sd", type=float, default=0.01)
    ap.add_argument("--min-called", type=float, default=0.90)
    ap.add_argument("--min-distinct", type=int, default=3)
    ap.add_argument("--out", default="S3_lvef.tsv")
    ap.add_argument("--list", action="store_true")
    a = ap.parse_args()

    if a.list:
        do_list(a)
        return

    # ---- G1 / G2 : the base table is the one the paper reports -----------
    print("=" * 74)
    print("GATES on the base table")
    print("=" * 74)
    ev = pd.read_csv(a.events, sep="\t")
    for c in (C_KEY, "FDR"):
        if c not in ev.columns:
            die(f"{a.events} has no '{c}' column. Columns: {list(ev.columns)}")
    dcol = resolve_dpsi(ev.columns, a.events)
    ev["dPSI"] = pd.to_numeric(ev[dcol], errors="coerce")
    ok(f"effect column = '{dcol}'")

    ev["ekey"] = ekey_from_key(ev[C_KEY])

    # G0: the key's first three fields must BE chr:exonStart:exonEnd. If they are not,
    # the exon-level collapse is collapsing on something else and the rMATS join below
    # is meaningless. This is cheap and the table carries the columns to check it with.
    if all(c in ev.columns for c in ("chr", "exonStart", "exonEnd")):
        built = (ev["chr"].astype(str) + ":" + ev["exonStart"].astype(str) + ":" +
                 ev["exonEnd"].astype(str))
        bad = int((built != ev["ekey"]).sum())
        if bad:
            ex = ev.loc[built != ev["ekey"], [C_KEY, "chr", "exonStart", "exonEnd"]].head(3)
            die(f"G0 FAILED: {bad:,} of {len(ev):,} rows where key[:3] != "
                f"chr:exonStart:exonEnd.\n{ex.to_string(index=False)}")
        ok(f"G0 key[:3] == chr:exonStart:exonEnd for all {len(ev):,} rows")
    else:
        print("  [warn] no chr/exonStart/exonEnd columns; G0 skipped")

    rep = (ev.sort_values("FDR", kind="mergesort")
             .drop_duplicates(subset=["ekey"], keep="first"))
    if len(rep) != EXPECT_TESTED:
        die(f"G1 FAILED: {a.events} collapses to {len(rep):,} exon-level events, "
            f"expected {EXPECT_TESTED:,}. Wrong table or wrong convention.")
    ok(f"G1 exon-level tested = {len(rep):,}")

    nsig = int(((rep.FDR < FDR_MAX) & (rep.dPSI.abs() >= DPSI_MIN)).sum())
    if nsig != EXPECT_SIG:
        die(f"G2 FAILED: {nsig} significant exons, expected {EXPECT_SIG}.")
    ok(f"G2 significant = {nsig}")

    # The table's own sd_psi (over both groups) is an INDEPENDENT census of the flat
    # events the old script regressed on. Not a gate -- a different quantity from the
    # within-group sd computed below -- but it says how big the problem was.
    if "sd_psi" in rep.columns:
        s = pd.to_numeric(rep["sd_psi"], errors="coerce")
        n0 = int((s == 0).sum())
        print(f"  [census] base table sd_psi == 0 : {n0:,} exons")
        print(f"  [census] base table sd_psi < {a.min_sd}: "
              f"{int((s < a.min_sd).sum()):,} exons")
        if n0 == 0:
            print("  [note] ZERO flat events in the base table -- it is already")
            print("         variance-filtered, which independently confirms the")
            print("         Methods. The old script's constant-PSI artefacts therefore")
            print("         came from the v1 table, not from this one. The within-group")
            print("         filter below is STILL required: sd over both groups can be")
            print("         non-zero while sd within one group is exactly zero.")

    gene_map = None
    if "gene" in rep.columns:
        cols = ["ekey", "gene"] + (["GeneID"] if "GeneID" in rep.columns else [])
        gene_map = rep[cols].drop_duplicates(subset=["ekey"])

    # ---- samples ---------------------------------------------------------
    cov, _ = load_covariates(a.cov)
    ph, runcol, lvef_col, _ = load_pheno(a.pheno, a.lvef_col)
    if runcol is None:
        die("no run-accession column in the pheno file; run --list")
    if lvef_col is None:
        die("no LVEF column found; run --list and pass --lvef-col")

    ph = ph.rename(columns={runcol: "run", lvef_col: "LVEF"})
    ph["run"] = ph["run"].astype(str).str.strip()
    ph["LVEF"] = pd.to_numeric(ph["LVEF"], errors="coerce")

    # UNITS. This cohort records LVEF as a FRACTION (observed range 0.05-0.80), not a
    # percentage. A coefficient per unit LVEF would then mean "per 100 percentage
    # points", which is not a quantity anyone reports. Rescale and say so.
    _v = ph["LVEF"].dropna()
    if len(_v) and _v.max() <= 1.0:
        print(f"  LVEF read as a FRACTION (max {_v.max():.3g}) -> rescaled to "
              "percentage points")
        ph["LVEF"] = ph["LVEF"] * 100.0
    elif len(_v):
        print(f"  LVEF read as PERCENTAGE POINTS (max {_v.max():.3g})")

    r1, r2 = bam_runs(a.b1), bam_runs(a.b2)
    order = pd.DataFrame({"run": r1 + r2,
                          "arm": ["g1"] * len(r1) + ["g2"] * len(r2)})
    order["pos"] = list(range(len(r1))) + list(range(len(r2)))
    if order["run"].isna().any():
        die("some bam paths carry no run accession; cannot align IncLevel columns")

    # G6: the bam lists must be exactly as wide as the PSI columns. Checked BEFORE the
    # covariate join, because a mis-parsed bam list looks like missing metadata.
    w1, w2 = inclevel_widths(a.mats)
    if (w1, w2) != (len(r1), len(r2)):
        die(f"G6 FAILED: bam lists give {len(r1)}/{len(r2)} samples but IncLevel1/2 "
            f"carry {w1}/{w2}. The bam list is not parsed the way rMATS read it "
            "(these files are one comma-separated line, not one path per line).")
    ok(f"G6 bam lists {len(r1)}/{len(r2)} == IncLevel widths {w1}/{w2}")

    s = (order.merge(cov[["run", "group", "sex", "age", "RIN"]], on="run", how="left")
              .merge(ph[["run", "LVEF"]].drop_duplicates("run"), on="run", how="left"))

    extra = []
    if a.extra_cov_file:
        if not a.extra_cov:
            die("--extra-cov-file given without --extra-cov")
        extra = [c.strip() for c in a.extra_cov.split(",") if c.strip()]
        ex = pd.read_csv(a.extra_cov_file, sep="\t")
        rcol = ex.columns[0]
        ex = ex.rename(columns={rcol: "run"})
        ex["run"] = ex["run"].astype(str).str.strip()
        miss = [c for c in extra if c not in ex.columns]
        if miss:
            die(f"{a.extra_cov_file} lacks {miss}. Columns: {list(ex.columns)}")
        s = s.merge(ex[["run"] + extra].drop_duplicates("run"), on="run", how="left")
        print(f"  extra covariates: {extra} from {a.extra_cov_file}")

    want = s[s["group"].astype(str).str.upper() == a.group.upper()].copy()
    if want.empty:
        die(f"no samples in group '{a.group}'. Groups present: "
            f"{dict(s['group'].value_counts(dropna=False))}")

    before = len(want)
    want = want.dropna(subset=["LVEF", "sex", "age", "RIN"] + extra)
    print(f"\n{a.group}: {before} samples in the bam lists, "
          f"{len(want)} with LVEF + sex + age + RIN")
    if len(want):
        lv = want["LVEF"]
        print(f"  LVEF (percentage points): median {lv.median():.1f}, "
              f"IQR {lv.quantile(.25):.1f}-{lv.quantile(.75):.1f}, "
              f"range {lv.min():.1f}-{lv.max():.1f}")
        if lv.max() - lv.min() < 5:
            die("LVEF spans under 5 percentage points in this group; there is no "
                "gradient to regress against")
        # The predictor's shape decides which test is trustworthy. Clinical LVEF is
        # recorded in multiples of 5 and a transplant cohort has a long right tail, so
        # a least-squares fit can be carried by two or three high-leverage hearts.
        q1, q3 = lv.quantile(.25), lv.quantile(.75)
        hi = int((lv > q3 + 1.5 * (q3 - q1)).sum())
        print(f"  distinct LVEF values: {lv.nunique()} over {len(lv)} samples"
              f"  (most common: {lv.mode().iloc[0]:.0f}, "
              f"n = {int((lv == lv.mode().iloc[0]).sum())})")
        print(f"  above Q3 + 1.5 IQR  : {hi} samples"
              f"{'  <- high-leverage points for any least-squares fit' if hi else ''}")
    if len(want) < 20:
        die(f"G5 FAILED: only {len(want)} complete-covariate {a.group} samples")
    ok(f"G5 analysed n = {len(want)}")

    # ---- per-sample PSI from rMATS --------------------------------------
    print("\n" + "=" * 74)
    print("Per-sample PSI")
    print("=" * 74)
    usecols = MATS_COORD + ["IncLevel1", "IncLevel2"]
    mats = pd.read_csv(a.mats, sep="\t")
    miss = [c for c in usecols if c not in mats.columns]
    if miss:
        die(f"{a.mats} lacks {miss}. Columns: {list(mats.columns)}")
    mats = mats[usecols].copy()
    mats["ekey"] = ekey_from_mats(mats)
    mats["fkey"] = fullkey_from_mats(mats)

    # ------------------------------------------------------------------ THE ROW FIX
    # rMATS emits one row per TRANSCRIPT CONTEXT, so an exon appears several times with
    # the same chr:exonStart:exonEnd and different flanking coordinates -- and those rows
    # carry DIFFERENT IncLevel values. Taking "the first row per exon" would pick an
    # arbitrary context's PSI while the base table's dPSI and FDR came from the
    # LOWEST-FDR context. That is the project's standing trap: retention is a SET
    # question, every paired statistic is a ROW question.
    #
    # So the base table's representative row chooses, and PSI is looked up at that same
    # FULL key. v1 of this script joined on the 3-field exon key and kept the first
    # match, which was wrong in exactly the documented way.
    keep_e = set(rep["ekey"])
    matched_e = mats.loc[mats["ekey"].isin(keep_e), "ekey"].nunique()
    frac_e = matched_e / max(len(rep), 1)
    print(f"  rMATS rows: {len(mats):,}  ->  {mats['ekey'].nunique():,} distinct exons")
    print(f"  base exons found in rMATS   : {matched_e:,} of {len(rep):,} "
          f"({frac_e:.1%})")
    print(f"  mean transcript contexts/exon: "
          f"{len(mats) / max(mats['ekey'].nunique(), 1):.2f}")

    dup_f = int(mats["fkey"].duplicated().sum())
    if dup_f:
        print(f"  [warn] {dup_f:,} duplicate FULL keys in the rMATS table; "
              "keeping the first of each")
    midx = mats.drop_duplicates(subset=["fkey"], keep="first").set_index("fkey")

    want_full = rep[[C_KEY, "ekey"]].rename(columns={C_KEY: "fkey"})
    have = want_full["fkey"].isin(midx.index)
    frac_f = float(have.mean())
    print(f"  base REPRESENTATIVE rows matched at the full key: "
          f"{int(have.sum()):,} of {len(have):,} ({frac_f:.1%})")

    if frac_f < MIN_JOIN_FRAC:
        for shift in (+1, -1):
            alt = (mats["chr"].astype(str) + ":" +
                   (mats["exonStart_0base"] + shift).astype(str) + ":" +
                   mats["exonEnd"].astype(str))
            print(f"  diagnostic, exonStart{shift:+d}: "
                  f"{alt.isin(keep_e).sum():,} rows match on 3 fields")
        nochr = mats["ekey"].str.replace("^chr", "", regex=True)
        print(f"  diagnostic, no 'chr' prefix: {nochr.isin(keep_e).sum():,} match")
        die(f"G3 FAILED: only {frac_f:.1%} of base representative rows found in "
            f"{a.mats} at the full key. Do not proceed.")
    ok(f"G3 full-key join = {frac_f:.1%}  (exon-level coverage {frac_e:.1%})")

    sel = want_full[have]
    sub = midx.loc[sel["fkey"].to_numpy()]

    def psi_block(s, n_expected, name):
        M = (s.astype(str).str.split(",", expand=True)
             .apply(pd.to_numeric, errors="coerce"))
        if M.shape[1] != n_expected:
            die(f"{name} has {M.shape[1]} values per row, expected {n_expected}")
        return M.to_numpy(float)

    P1 = psi_block(sub["IncLevel1"], len(r1), "IncLevel1")
    P2 = psi_block(sub["IncLevel2"], len(r2), "IncLevel2")
    PSI = np.concatenate([P1, P2], axis=1)          # events x all samples

    colpos = {("g1", i): i for i in range(len(r1))}
    colpos.update({("g2", i): len(r1) + i for i in range(len(r2))})
    idx = [colpos[(row.arm, row.pos)] for row in want.itertuples()]
    PSI = PSI[:, idx]
    print(f"  PSI matrix: {PSI.shape[0]:,} events x {PSI.shape[1]} {a.group} samples")

    # ---- the filter that the old script did not have --------------------
    print("\n" + "=" * 74)
    print("PSI-variance filter")
    print("=" * 74)
    n_s = PSI.shape[1]
    called = np.isfinite(PSI).sum(axis=1) / n_s
    with np.errstate(invalid="ignore"):
        sd = np.nanstd(PSI, axis=1)
    distinct = np.array([len(np.unique(row[np.isfinite(row)])) for row in PSI])

    f_called = called >= a.min_called
    f_sd = np.nan_to_num(sd, nan=0.0) >= a.min_sd
    f_dist = distinct >= a.min_distinct
    keep = f_called & f_sd & f_dist

    print(f"  events joined                : {len(keep):,}")
    print(f"  fail called >= {a.min_called:<5}        : {(~f_called).sum():,}")
    print(f"  fail sd >= {a.min_sd:<9}        : {(~f_sd).sum():,}")
    print(f"  fail distinct >= {a.min_distinct:<3}          : {(~f_dist).sum():,}")
    print(f"  CONSTANT PSI (sd == 0)       : {int((np.nan_to_num(sd) == 0).sum()):,}"
          "   <- the old script's artefacts")
    print(f"  retained                     : {int(keep.sum()):,}")
    if keep.sum() < 100:
        die("fewer than 100 events survive the filter; the inputs are wrong")

    PSI = PSI[keep]
    kept_keys = sel["ekey"].to_numpy()[keep]

    # ---- regression ------------------------------------------------------
    print("\n" + "=" * 74)
    print(f"PSI ~ LVEF + sex + age + RIN, within {a.group}")
    print("=" * 74)
    sex = pd.Categorical(want["sex"].astype(str).str.upper()).codes.astype(float)
    base_cols = [want["LVEF"].to_numpy(float), sex,
                 want["age"].to_numpy(float), want["RIN"].to_numpy(float)]
    Xfull = np.column_stack(base_cols)
    Xext = (np.column_stack(base_cols + [want[c].to_numpy(float) for c in extra])
            if extra else None)
    if extra:
        print(f"  extended model adds: {', '.join(extra)}")

    rows = []
    for i in range(PSI.shape[0]):
        y = PSI[i]
        m = np.isfinite(y)
        yy, XX = y[m], Xfull[m]
        if len(np.unique(yy)) < a.min_distinct or yy.std() < a.min_sd:
            continue                                  # G4: never fit a flat vector
        b, se, p = ols_lvef(yy, XX)
        # PRIMARY: the same covariate-adjusted model on RANKS of PSI and of LVEF.
        # PSI is bounded and outlier-prone and LVEF is heavily tied, so a raw
        # least-squares t-test is anti-conservative here; ranking both fixes that
        # while keeping sex, age and RIN in the model (which a bare Spearman cannot).
        ry = stats.rankdata(yy)
        XR = XX.copy()
        XR[:, 0] = stats.rankdata(XX[:, 0])
        b_rank, _ser, p_rank = ols_lvef(ry, XR)
        rho, prho = spearman(XX[:, 0], yy)
        row = [kept_keys[i], int(m.sum()), float(yy.mean()), float(yy.std()),
               b, se, p, b_rank, p_rank, rho, prho]
        if extra:
            XE = Xext[m].copy()
            XE[:, 0] = stats.rankdata(XE[:, 0])
            b_e, _se_e, p_e = ols_lvef(ry, XE)
            row += [b_e, p_e]
        rows.append(tuple(row))

    cols = ["ekey", "n", "psi_mean", "psi_sd", "beta_lvef", "se_lvef", "p_lvef",
            "beta_rank", "p_rank", "rho_lvef", "p_rho"]
    if extra:
        cols += ["beta_rank_adj", "p_rank_adj"]
    res = pd.DataFrame(rows, columns=cols)
    if res.empty:
        die("no event produced a fit; nothing to report")

    # beta is PSI per percentage point of LVEF, which is tiny by construction.
    # Per 10 pp is the readable unit and the one to quote.
    res["dPSI_per_10pp_lvef"] = res["beta_lvef"] * 10.0

    # G4
    bad_sd = int((res.psi_sd == 0).sum())
    bad_rho = int(res.rho_lvef.isna().sum())
    if bad_sd or bad_rho:
        die(f"G4 FAILED: {bad_sd} events with sd == 0 and {bad_rho} with rho = nan "
            "survived the filter. This is the old script's bug; stop.")
    ok(f"G4 no constant-PSI or nan-rho events among {len(res):,} fits")

    res["FDR_lvef"] = bh(res.p_lvef.to_numpy())
    res["FDR_rank"] = bh(res.p_rank.to_numpy())
    res["FDR_rho"] = bh(res.p_rho.to_numpy())
    if gene_map is not None:
        res = res.merge(gene_map, on="ekey", how="left")
        front = ["ekey", "gene"] + (["GeneID"] if "GeneID" in res.columns else [])
        res = res[front + [c for c in res.columns if c not in front]]
    res = res.sort_values("p_rank", kind="mergesort")

    sig_ols = int((res.FDR_lvef < FDR_MAX).sum())
    sig_rank = int((res.FDR_rank < FDR_MAX).sum())
    sig_rho = int((res.FDR_rho < FDR_MAX).sum())
    n_med = int(res.n.median())
    exp_med_rho = 0.6745 / np.sqrt(max(n_med - 1, 2))  # half-normal median, null rho

    print(f"\n  events fitted                      : {len(res):,}")
    print(f"  FDR < 0.05, rank-adjusted (PRIMARY): {sig_rank}")
    print(f"  FDR < 0.05, raw least squares      : {sig_ols}")
    print(f"  FDR < 0.05, bare Spearman          : {sig_rho}")
    print(f"  smallest p, rank-adjusted          : {res.p_rank.min():.3g}")
    print(f"  median |rho|                       : "
          f"{res.rho_lvef.abs().median():.4f}   (null expectation "
          f"{exp_med_rho:.4f} at n = {n_med})")
    print(f"  |rho| max                          : {res.rho_lvef.abs().max():.3f}")
    print(f"  max |dPSI| per 10 pp of LVEF       : "
          f"{res.dPSI_per_10pp_lvef.abs().max():.4f}   (locked threshold 0.10)")

    for lbl, col in (("rank-adjusted", "p_rank"), ("raw least squares", "p_lvef"),
                     ("bare Spearman", "p_rho")):
        ks = stats.kstest(res[col].dropna(), "uniform")
        print(f"  KS uniformity, {lbl:<18}: D = {ks.statistic:.4f}, "
              f"p = {ks.pvalue:.3g}")
    print("  (KS assumes INDEPENDENT tests. These are not independent -- neighbouring")
    print("   exons share samples and coverage -- so a small p here is expected even")
    print("   under a true null. Read D, which is the size of the deviation; BH on the")
    print("   same p-values handles the dependence far better and is what decides.)")

    # How close is the nearest miss? BH's critical value at rank 1 is 0.05/m, so the
    # smallest p-value alone decides whether anything survives. Report the margin
    # rather than only "zero", so a borderline event is visible and not buried.
    m_tests = int(res.p_rank.notna().sum())
    bh1 = FDR_MAX / m_tests
    pmin = float(res.p_rank.min())
    print(f"\n  BH boundary at rank 1   : {bh1:.3g}   (= {FDR_MAX} / {m_tests:,})")
    print(f"  smallest p_rank         : {pmin:.3g}"
          f"   -> {'SURVIVES' if pmin <= bh1 else f'misses by {pmin / bh1:.3f}x'}")

    show = [c for c in ("gene", "ekey", "n", "psi_sd", "dPSI_per_10pp_lvef",
                        "rho_lvef", "p_rank", "FDR_rank") if c in res.columns]
    print("\n  Top 5 by the primary model -- inspect these before claiming a null:")
    print(res[show].head(5).to_string(index=False))

    # ---- GLOBAL SIGN SKEW -----------------------------------------------------
    # A null per event is not the same as a null overall. If PSI shifts genome-wide
    # with severity, every coefficient leans one way by a trivial amount and no single
    # event is significant. That is a real observation and a nuisance-shaped one, so
    # it is measured here and, when extra covariates are supplied, measured again with
    # them in the model -- same event set, one run, no cross-run bookkeeping.
    print("\n" + "-" * 74)
    print("  GLOBAL SIGN SKEW of the LVEF coefficient")
    print("-" * 74)

    def skew(col, label):
        v = res[col].dropna()
        neg = int((v < 0).sum())
        bt = stats.binomtest(neg, len(v), 0.5).pvalue if len(v) else float("nan")
        print(f"  {label:<34}: {neg:,}/{len(v):,} negative = {neg / len(v):.1%}"
              f"   binomial p = {bt:.3g}")
        return neg / len(v)

    f_base = skew("beta_rank", "base (LVEF + sex + age + RIN)")
    print(f"  {'mean rho / median rho':<34}: {res.rho_lvef.mean():+.5f} / "
          f"{res.rho_lvef.median():+.5f}")
    if extra:
        f_ext = skew("beta_rank_adj", f"+ {', '.join(extra)}")
        d_base, d_ext = abs(f_base - 0.5), abs(f_ext - 0.5)
        drop = 1 - d_ext / d_base if d_base > 0.001 else float("nan")
        print(f"\n  excess over 50%: {d_base * 100:.1f} points -> "
              f"{d_ext * 100:.1f} points"
              + (f"  ({drop:+.0%} change)" if np.isfinite(drop) else ""))
        if d_base < 0.01:
            print("  >>> There is no skew to explain: the base model's excess is under")
            print("      1 point. Whatever the adjusted model shows, it is not a")
            print("      correction to anything. Report neither.")
        elif d_ext < 0.01:
            print("  >>> The skew COLLAPSES once those covariates are in the model.")
            print("      Read it as mediated by them, not as splicing tracking")
            print("      severity. One Methods clause; no new claim.")
        elif drop > 0.4:
            print("  >>> The skew SHRINKS substantially but does not vanish. Report")
            print("      both numbers; attribute the explained part and leave the")
            print("      remainder as unexplained rather than naming a mechanism.")
        else:
            print("  >>> The skew SURVIVES. It is then NOT explained by these")
            print("      covariates. Do not report it as a finding on that basis --")
            print("      library depth and transcript integrity are the next")
            print("      candidates, and PSI is a count ratio sensitive to both.")
    else:
        print("\n  ⚠ The binomial p above assumes INDEPENDENT events and these are not:")
        print("    38k exons share 160 samples and their coverage, so the p-value is")
        print("    badly inflated and must not be quoted. What matters is the SIZE of")
        print("    the excess over 50% and whether it survives adjustment.")
        print("    Re-run with, e.g.:")
        print("      --extra-cov-file magnet_expr.cmfrac.tsv "
              "--extra-cov CM_score_frac,FB_score")
        print("    to measure both in one pass on the same event set.")

    res.to_csv(a.out, sep="\t", index=False)
    print(f"\n  wrote {a.out}  ({len(res):,} rows)")

    print("\n" + "=" * 74)
    print("WHAT TO WRITE IN S3")
    print("=" * 74)

    if sig_ols > 5 * max(sig_rank, 1):
        print(f"  ⚠ READ THIS FIRST. Raw least squares gives {sig_ols} hits and the")
        print(f"    rank-adjusted model gives {sig_rank}. That gap is not biology: PSI")
        print("    is bounded and outlier-prone, LVEF is recorded in multiples of 5")
        print("    with a long right tail, so the t-test on raw values is")
        print("    anti-conservative and a handful of hearts can carry a fit.")
        print("    **Quote the rank-adjusted count.** Do not report the raw count as")
        print("    a finding; report it, if at all, as the reason the rank test is")
        print("    primary. This is the same class of error as the constant-PSI")
        print("    artefacts -- a model assumption the data do not meet -- one level up.")
        print()

    if sig_rank == 0:
        mx = float(res.dPSI_per_10pp_lvef.abs().max())
        n_over = int((res.dPSI_per_10pp_lvef.abs() >= DPSI_MIN).sum())
        print(f"  A null, and reportable. Within {a.group} (n = {len(want)}), no")
        print(f"  skipped-exon PSI tracked LVEF at FDR < 0.05 across {len(res):,}")
        print("  events passing the variance filter, in a covariate-adjusted rank")
        print(f"  model. Median |rho| was {res.rho_lvef.abs().median():.3f} against a")
        print(f"  null expectation of {exp_med_rho:.3f}.")
        # Do NOT assert the effect is under the threshold without checking it.
        if mx < DPSI_MIN:
            print(f"  The largest effect of any kind was {mx:.3f} PSI per 10 points of")
            print(f"  LVEF -- below the {DPSI_MIN:.2f} threshold this paper requires of")
            print("  a real event, so nothing here is even of reportable size.")
        else:
            print(f"  ⚠ The largest effect reaches {mx:.4f} PSI per 10 points of LVEF,")
            print(f"    which is NOT below the {DPSI_MIN:.2f} threshold — "
                  f"{n_over} event(s)")
            print("    reach it. Do not write 'no event of reportable size'. Write that")
            print("    no event survived multiple testing, and that the largest")
            print("    slope was of borderline size and did not.")
        if pmin > bh1 and pmin / bh1 < 2:
            print(f"  ⚠ The nearest miss is close: p = {pmin:.3g} against a boundary of")
            print(f"    {bh1:.3g}. Name that event in the supplement rather than only")
            print("    reporting a zero — a reviewer who recomputes will find it.")
    else:
        print(f"  {sig_rank} events survive the rank-adjusted model. For each, check")
        print("  psi_sd, n and dPSI_per_10pp_lvef in the output before writing it up:")
        print("  a hit at psi_sd near the 0.01 floor, or an effect well under 0.10 per")
        print("  10 points, is statistically detectable and biologically nothing.")
    print()
    print(f"  State the filter: {int(len(keep) - keep.sum()):,} of {len(keep):,} joined")
    print("  events were dropped before any model was fitted "
          f"({int((np.nan_to_num(sd) == 0).sum()):,} with constant PSI within "
          f"{a.group}).")
    print()


if __name__ == "__main__":
    main()
