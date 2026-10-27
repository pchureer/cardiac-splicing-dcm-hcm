#!/usr/bin/env python
"""
82_build_s4.py -- Supplementary Table S4: cross-cohort magnitude, and the matched-pair
panel that explains it.

WHAT S4 CARRIES
---------------
The Discussion's strongest methods claim:

    arm A  MAGNet_matched -> LBNL_matched    platform only, same hearts   0.584 (0.478-0.693)
    arm B  MAGNet_matched -> MAGNet_other    people only, same platform   0.988 (0.889-1.089)
    A x B = 0.577, against an observed cross-cohort slope of 0.521

⚠ THOSE ARE 50 bp VALUES. At 125 bp, measured: arm A 0.690 (0.570-0.806),
arm B 0.968 (0.875-1.069), A x B = 0.668 against an observed C of 0.600 — a
closure gap of 11.3%, against 10.7% at 50 bp. The platform cost restates from
41% to 29%. Every EXPECT_ constant below is overridable and REQUIRED when
TBL_SUFFIX is set; do not read a read length off this docstring.

plus the per-cohort transfer slopes (Colorado 0.521, Brigham 0.476, Stanford 0.386,
Fuwai 0.243) and the reason no effect-size magnitude is compared across cohorts anywhere
in the paper.

WHY THIS RECOMPUTES RATHER THAN TRANSCRIBES
-------------------------------------------
The record gives slopes but not the convention they were computed under, and `19`/`36`
predate `24_lock_event_convention.py`. Transcribing a slope whose convention is unknown
is how `RESULT_magnitude_resolved.md` ended up contradicting Figure 4 by 0.281. So the
slopes are recomputed from `mp_SE_*.all_events.tsv` and the build **refuses unless it
reproduces EXPECT_A and EXPECT_B (0.584 and 0.988 at 50 bp)**.

⚠ THE LAMBDA CONVENTION IS NOT RECORDED, SO THE SCRIPT SEARCHES FOR IT
----------------------------------------------------------------------
Deming regression needs the error-variance ratio lambda = var(x_err)/var(y_err), and
"taken from the per-event standard errors" (Methods) does not pin down whether that means
the mean of squared SEs, the median, per-event weighting, or lambda = 1 (orthogonal).
Each gives a different slope.

Rather than guess, this script computes the slope under **every** convention and reports
which one reproduces the recorded EXPECT_A and EXPECT_B **in both arms at once**. That is the
same approach `65` took for S7, and it establishes the convention by numbers instead of
by label. The winning convention's name is what the Methods should state.

**If no convention reproduces both, the script refuses** -- that would mean the recorded
slopes came from an implementation detail not captured here, and the honest response is
to find it rather than publish a slope this script happened to produce.

GATES
-----
  G1  the three arms share exactly EXPECT_COMMON events
      (203 at 50 bp, 161 at 125 bp)
  G2  direction agreement EXPECT_DIR_A and EXPECT_DIR_B
      (162/203 and 174/203 at 50 bp; 130/161 and 140/161 at 125 bp)
  G3  one lambda convention reproduces BOTH recorded arm slopes within tolerance
  G4  arm A's CI upper bound < arm B's CI lower bound (the intervals must not overlap --
      that non-overlap IS the claim)

USAGE
-----
  conda activate rnaseq_tools
  cd /data/project/chureeratp/RNA/HEART

  python 82_build_s4.py --list
  python 82_build_s4.py
  python 82_build_s4.py --tol 0.02      # loosen the reproduction tolerance
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd

NKEEP_SE = 3
C_KEY = "key"
C_DPSI_CANDIDATES = ("dPSI_adj", "dPSI_raw", "dPSI")
SE_CANDIDATES = ("se_adj", "se", "se_dpsi")
# Read length. Drives the output name, all four input arms, and the record.
SUF = os.environ.get("TBL_SUFFIX", "")
OUT = f"S4{SUF}.tsv"

ARMS = {
    "selection": f"mp_SE_selection{SUF}.all_events.tsv",
    "magnet_matched": f"mp_SE_magnet_matched{SUF}.all_events.tsv",
    "lbnl_matched": f"mp_SE_lbnl_matched{SUF}.all_events.tsv",
    "magnet_other": f"mp_SE_magnet_other{SUF}.all_events.tsv",
}

# ------------------------------------------------------------------ the record
EXPECT_COMMON = 203
EXPECT_A = 0.584
EXPECT_B = 0.988
EXPECT_A_CI = (0.478, 0.693)
EXPECT_B_CI = (0.889, 1.089)
EXPECT_DIR_A = 162
EXPECT_DIR_B = 174
EXPECT_PRODUCT = 0.577
EXPECT_OBSERVED_C = 0.521
# per-cohort transfer slopes, on MAGNet-significant SE events (Figure 4's numbers)
PER_COHORT = {"Colorado": 0.521, "Brigham": 0.476, "Stanford": 0.386, "Fuwai": 0.243}
# --------------------------------------------------------------- RETIRED
# GENOME_WIDE = {"Colorado": 0.571, "Brigham": 0.757,
#                "Fuwai": 0.359, "Stanford": 0.578}
#
# DROPPED FROM S4 on 20 Oct. Deliberate, and recorded here rather than lost.
#
# These are script 21's PRE-EXCLUSION, JCEC-based genome-wide slopes on the
# matched JC-intersect-JCEC subset, from RESULT_magnitude_resolved.md — a
# document whose own banner reads "Quote nothing numeric from this document
# in the manuscript or the supplement". They were being written into a
# supplementary table in defiance of that. 62_build_s5.py exists for the same
# reason and guards against it by aborting unless its JC column reproduces
# the manuscript's slopes; 82 had no equivalent guard.
#
# They are neither wrong nor unsourced — they are superseded. Three distinct
# Brigham quantities exist and must never be interchanged:
#     0.476  manuscript — post-exclusion, JC, MAGNet-significant (45)
#     0.662  post-exclusion, JC, unrestricted forward (22)
#     0.757  pre-exclusion, JCEC, matched subset (21/20)   <- this one
# and the source document explains the third: restricting to the
# JC-intersect-JCEC set raises Brigham from 0.662 to 0.757, an event-set
# effect rather than a read-length one.
#
# No 125 bp counterpart exists, because nothing at 125 bp is pre-exclusion
# JCEC. Patch 138 briefly filled this with 22's 125 bp FORWARD slopes — a
# different quantity under the same label — which is the defect removing the
# column fixes.
# --------------------------------------------------------------------------

# --------------------------------------------------------------------------
# 125 bp RECORD, placed AFTER every literal it rewrites.
#
# 82 is reproduce-or-refuse: G1 dies if no selection-set definition gives
# EXPECT_COMMON events, G2 compares EXPECT_DIR_*, and G3 dies if no lambda
# convention reproduces BOTH EXPECT_A and EXPECT_B. A 50 bp record against the
# 125 bp arms therefore REFUSES A CORRECT TABLE, with a message that reads
# like a data problem. So all three groups are REQUIRED when TBL_SUFFIX is set.
#
# Thirteen constants grouped into three variables; each asserts SET EQUALITY
# against the keys already present, so a typo names what was expected.
_S4_KEYS = ("common", "A", "B", "A_CI", "B_CI", "dirA", "dirB",
            "product", "C")
_S4_REQ = ("S4_RECORD", "S4_PER_COHORT")
_s4env = {_k: os.environ.get(_k, "").strip() for _k in _S4_REQ}
if SUF:
    _miss = [_k for _k in _S4_REQ if not _s4env[_k]]
    if _miss:
        raise SystemExit(
            "REFUSING — TBL_SUFFIX=%s but these are unset: %s\n\n"
            "  82's gates DIE rather than warn, so a 50 bp record against the\n"
            "  125 bp arms refuses a CORRECT table — G1 would report 'no\n"
            "  selection-set definition gives 203 events' and G3 'no\n"
            "  convention reproduces both 0.584 and 0.988', neither of which\n"
            "  is a data problem.\n\n    %s\n"
            % (SUF, ", ".join(_miss), 'TBL_SUFFIX=_rl125 \\\n    S4_RECORD=common=161,A=0.690,B=0.968,A_CI=0.570:0.806,B_CI=0.875:1.069,dirA=130,dirB=140,product=0.668,C=0.600 \\\n    S4_PER_COHORT=Colorado=0.5999,Brigham=0.5099,Stanford=0.4008,Fuwai=0.2948 \\\n    python3 82_build_s4.py'))

if _s4env["S4_RECORD"]:
    _d = dict(_p.split("=", 1) for _p in _s4env["S4_RECORD"].split(","))
    assert set(_d) == set(_S4_KEYS), (
        "S4_RECORD must name exactly %s, got %s"
        % (sorted(_S4_KEYS), sorted(_d)))
    EXPECT_COMMON = int(_d["common"])
    EXPECT_A = float(_d["A"])
    EXPECT_B = float(_d["B"])
    # CIs use ':' internally because ',' separates the entries
    EXPECT_A_CI = tuple(float(_x) for _x in _d["A_CI"].split(":"))
    EXPECT_B_CI = tuple(float(_x) for _x in _d["B_CI"].split(":"))
    EXPECT_DIR_A = int(_d["dirA"])
    EXPECT_DIR_B = int(_d["dirB"])
    EXPECT_PRODUCT = float(_d["product"])
    EXPECT_OBSERVED_C = float(_d["C"])
    assert len(EXPECT_A_CI) == 2 and len(EXPECT_B_CI) == 2, \
        "A_CI and B_CI must each be lo:hi"

if _s4env["S4_PER_COHORT"]:
    _d = dict(_p.split("=", 1) for _p in _s4env["S4_PER_COHORT"].split(","))
    assert set(_d) == set(PER_COHORT), (
        "S4_PER_COHORT must name exactly %s, got %s"
        % (sorted(PER_COHORT), sorted(_d)))
    PER_COHORT = {_k: float(_v) for _k, _v in _d.items()}

# GENOME_WIDE is retired — see the comment above its former definition.
# There is deliberately no S4_GENOME_WIDE override: the quantity has no 125 bp
# counterpart, so there is nothing to override it with.


def _require_arms():
    """Refuse early if a suffixed arm is absent, rather than at G1."""
    if not SUF:
        return
    gone = [v for v in ARMS.values() if not os.path.exists(v)]
    if gone:
        raise SystemExit(
            "REFUSING — TBL_SUFFIX=%s but these arms are absent:\n%s\n"
            "  They come from the matched-pair build at 125 bp. Without them\n"
            "  G1's selection-set search would fail on missing data and read\n"
            "  as a stale-constant problem."
            % (SUF, "".join("    missing: %s\n" % g for g in gone)))
# --------------------------------------------------------------------------


def die(msg):
    print("\n[REFUSE] " + msg + "\n", file=sys.stderr)
    sys.exit(2)


def ok(msg):
    print("  [ok] " + msg)


def pick(cols, cands):
    for c in cands:
        if c in cols:
            return c
    return None


def load(path):
    if not (os.path.exists(path) and os.path.getsize(path)):
        return None
    t = pd.read_csv(path, sep="\t", low_memory=False)
    d = pick(t.columns, C_DPSI_CANDIDATES)
    if d is None:
        d = next((c for c in t.columns if str(c).lower().startswith("dpsi")), None)
    if d is None:
        die(f"{path}: no dPSI column. Columns: {list(t.columns)}")
    s = pick(t.columns, SE_CANDIDATES)
    if C_KEY not in t.columns:
        die(f"{path}: no '{C_KEY}'")
    t = t.copy()
    t["dPSI"] = pd.to_numeric(t[d], errors="coerce")
    t["SE"] = pd.to_numeric(t[s], errors="coerce") if s else np.nan
    t["FDR"] = pd.to_numeric(t["FDR"], errors="coerce") if "FDR" in t.columns else np.nan
    t["ekey"] = [":".join(str(k).split(":")[:NKEEP_SE]) for k in t[C_KEY]]
    srt = "FDR" if "FDR" in t.columns else d
    t = (t.sort_values(srt, kind="mergesort")
          .drop_duplicates(subset=["ekey"], keep="first")
          .dropna(subset=["dPSI"]))
    return t[["ekey", "dPSI", "SE", "FDR"]], d, s


# How the 203 is defined was NOT recorded. The record says the events were "selected
# from a MAGNet run sharing no sample with any arm", which is a SIGNIFICANCE filter on
# the selection set, not mere membership in it -- the raw three-arm intersection is
# 5,373. So the definition is searched for, the way 83 searched for S6's reference.
SEL_RULES = {
    "significant in selection (FDR<0.05 & |dPSI|>=0.10)":
        lambda t: (t.FDR < 0.05) & (t.dPSI.abs() >= 0.10),
    "significant in selection (FDR<0.05 & |dPSI|>=0.05)":
        lambda t: (t.FDR < 0.05) & (t.dPSI.abs() >= 0.05),
    "FDR<0.05 in selection only":
        lambda t: t.FDR < 0.05,
    "present in selection only (no significance filter)":
        lambda t: pd.Series(True, index=t.index),
}


def deming(x, y, lam):
    """Deming slope with error-variance ratio lam = var(x_err)/var(y_err)."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    n = len(x)
    if n < 3:
        return np.nan
    mx, my = x.mean(), y.mean()
    sxx = ((x - mx) ** 2).sum() / (n - 1)
    syy = ((y - my) ** 2).sum() / (n - 1)
    sxy = ((x - mx) * (y - my)).sum() / (n - 1)
    if sxy == 0:
        return np.nan
    disc = (syy - lam * sxx) ** 2 + 4 * lam * sxy ** 2
    return float((syy - lam * sxx + np.sqrt(disc)) / (2 * sxy))


_BASE = ["mean of squared SEs", "median of squared SEs", "squared mean SE",
         "mean of per-event ratios"]
CONVENTIONS = _BASE + [b + " [inverted]" for b in _BASE] + ["orthogonal (lambda = 1)"]


def lam_val(name, sx, sy):
    """lambda = var(x_err)/var(y_err), under one named convention.

    ⚠ THE INVERTED VARIANTS ARE NOT OPTIONAL. 83_build_s6.py identified this
    project's Deming estimator exactly, and it is 'mean of squared SEs [INVERTED]'
    with the reference on the x-axis -- 1 of 14 combinations, total error 0.001
    against 0.24-1.15 for every alternative. An earlier version of this script
    searched only the four non-inverted ratios plus orthogonal, which would have
    refused on G3 without ever trying the convention the project actually uses.
    → RESULT_s6_reconstructed.md
    """
    if name == "orthogonal (lambda = 1)":
        return 1.0
    inv = name.endswith(" [inverted]")
    base = name[: -len(" [inverted]")] if inv else name
    sx, sy = np.asarray(sx, float), np.asarray(sy, float)
    if not (np.isfinite(sx).all() and np.isfinite(sy).all()):
        return np.nan
    if inv:
        sx, sy = sy, sx
    if base == "mean of squared SEs":
        return float((sx ** 2).mean() / (sy ** 2).mean())
    if base == "median of squared SEs":
        return float(np.median(sx ** 2) / np.median(sy ** 2))
    if base == "squared mean SE":
        return float(sx.mean() ** 2 / sy.mean() ** 2)
    if base == "mean of per-event ratios":
        return float(np.mean((sx ** 2) / (sy ** 2)))
    return np.nan


def slope_for(x, y, sx, sy, conv, swap):
    """Deming slope. `swap` puts the OTHER arm on the x-axis -- a Deming slope is
    not symmetric, and which arm script 36 regressed on which was never recorded."""
    if swap:
        x, y, sx, sy = y, x, sy, sx
    lam = lam_val(conv, sx, sy)
    if not np.isfinite(lam):
        return np.nan, np.nan
    return deming(x, y, lam), lam


def boot_ci(x, y, sx, sy, conv, swap, n=2000, seed=11):
    rng = np.random.default_rng(seed)
    k = len(x)
    vals = []
    for _ in range(n):
        i = rng.integers(0, k, k)
        v, _l = slope_for(x[i], y[i], sx[i], sy[i], conv, swap)
        if np.isfinite(v):
            vals.append(v)
    if len(vals) < 100:
        return (np.nan, np.nan)
    return (float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=OUT)
    ap.add_argument("--tol", type=float, default=0.01)
    ap.add_argument("--boot", type=int, default=2000)
    ap.add_argument("--list", action="store_true")
    a = ap.parse_args()

    if a.list:
        print("=" * 74)
        print("S4 inputs -- the matched-pair panel")
        print("=" * 74)
        for name, path in ARMS.items():
            if os.path.exists(path):
                t, d, s = load(path)
                print(f"  OK       {path}")
                print(f"           {len(t):,} exon-level events, effect '{d}', "
                      f"SE '{s or 'ABSENT'}'")
            else:
                print(f"  MISSING  {path}")
        print(f"\n  The record: arm A {EXPECT_A} ({EXPECT_A_CI[0]}-{EXPECT_A_CI[1]}) platform only;")
        print(f"              arm B {EXPECT_B} ({EXPECT_B_CI[0]}-{EXPECT_B_CI[1]}) people only;")
        print(f"              {EXPECT_COMMON} common events, 10 DCM vs 14 NF in both arms.")
        print("\n  ⚠ If 'SE' is ABSENT above, only the orthogonal (lambda = 1)")
        print("    convention can be computed and the search cannot run.")
        return

    print("=" * 74)
    print("SUPPLEMENTARY TABLE S4 -- cross-cohort magnitude, matched-pair panel")
    print("=" * 74)

    tabs = {}
    for name, path in ARMS.items():
        got = load(path)
        if got is None:
            die(f"{path} missing. Run --list.")
        tabs[name] = got[0]
        print(f"  {name:<16} {len(got[0]):>7,} events   "
              f"effect '{got[1]}', SE '{got[2] or 'ABSENT'}'")

    # ---- G1: the common event set --------------------------------------
    inall = (set(tabs["magnet_matched"].ekey)
             & set(tabs["lbnl_matched"].ekey)
             & set(tabs["magnet_other"].ekey))
    sel = tabs["selection"]
    print(f"\n  events present in all three arms: {len(inall):,}")
    print(f"  SELECTION-SET DEFINITION SEARCH -- which gives {EXPECT_COMMON}?\n")
    common, selname = None, None
    for rname, rule in SEL_RULES.items():
        keys = set(sel.loc[rule(sel), "ekey"])
        n = len(inall & keys)
        hit = (n == EXPECT_COMMON)
        print(f"    {rname:<52}{n:>7,}" + ("   <<<" if hit else ""))
        if hit and common is None:
            common, selname = inall & keys, rname
    if common is None:
        die(f"G1 FAILED: no selection-set definition gives {EXPECT_COMMON} events. "
            "The arms may not be the ones the recorded slopes came from, or the "
            "selection applied a coverage filter this script does not model. Do NOT "
            "proceed on the closest -- the event set determines both slopes, and "
            "regression to the mean makes a different set give different numbers.")
    ok(f"G1 common event set = {len(common)}  via: {selname}")
    print(f"       >>> The Methods should state this selection rule; it is currently")
    print(f"           described only as 'selected from a MAGNet run sharing no")
    print(f"           sample with any arm', which does not pin it.")

    idx = sorted(common)
    def arm(name):
        t = tabs[name].set_index("ekey").loc[idx]
        return t.dPSI.to_numpy(float), t.SE.to_numpy(float)

    x, sx = arm("magnet_matched")
    ya, sya = arm("lbnl_matched")      # arm A: platform changes
    yb, syb = arm("magnet_other")      # arm B: people change

    # ---- G2: direction agreement ---------------------------------------
    dir_a = int(((x > 0) == (ya > 0)).sum())
    dir_b = int(((x > 0) == (yb > 0)).sum())
    print(f"\n  direction agreement: arm A {dir_a}/{len(idx)}, "
          f"arm B {dir_b}/{len(idx)}")
    bad = []
    if dir_a != EXPECT_DIR_A:
        bad.append(f"arm A {dir_a}, expected {EXPECT_DIR_A}")
    if dir_b != EXPECT_DIR_B:
        bad.append(f"arm B {dir_b}, expected {EXPECT_DIR_B}")
    if bad:
        die("G2 FAILED: " + "; ".join(bad))
    ok(f"G2 direction agreement {dir_a}/{EXPECT_DIR_A} and {dir_b}/{EXPECT_DIR_B}")

    # ---- G3: find the lambda convention --------------------------------
    print("\n" + "-" * 74)
    print("  LAMBDA CONVENTION SEARCH -- which one reproduces BOTH recorded slopes?")
    print("-" * 74)
    print("  S6's answer is the prior: 'mean of squared SEs [inverted]', reference")
    print("  on the x-axis. If S4 agrees, one estimator spec covers both supplements.")
    print(f"\n  {'x axis':<12}{'convention':<32}{'slope A':>9}{'slope B':>9}"
          f"{'err':>8}")
    winner, best, best_err = None, None, np.inf
    for swap in (False, True):
        for n in CONVENTIONS:
            sa, la_ = slope_for(x, ya, sx, sya, n, swap)
            sb, lb_ = slope_for(x, yb, sx, syb, n, swap)
            if not (np.isfinite(sa) and np.isfinite(sb)):
                continue
            err = abs(sa - EXPECT_A) + abs(sb - EXPECT_B)
            hit = abs(sa - EXPECT_A) <= a.tol and abs(sb - EXPECT_B) <= a.tol
            print(f"  {'LBNL/other' if swap else 'MAGNet':<12}{n:<32}"
                  f"{sa:>9.3f}{sb:>9.3f}{err:>8.3f}" + ("   <<<" if hit else ""))
            if err < best_err:
                best, best_err = (n, swap, sa, sb, la_, lb_), err
            # Among combinations that pass tolerance, take the BEST-fitting one, not
            # the first encountered. Several can pass at a loose tolerance and the
            # order of CONVENTIONS is not evidence about which script 36 used.
            if hit and (winner is None or err < winner[6]):
                winner = (n, swap, sa, sb, la_, lb_, err)

    if winner is None:
        die(f"G3 FAILED: no convention x direction reproduces both {EXPECT_A} and "
            f"{EXPECT_B} within +/-{a.tol}. Closest: "
            f"'{best[0]}'{' [swapped]' if best[1] else ''} at "
            f"A {best[2]:.3f}, B {best[3]:.3f} (err {best_err:.3f}). "
            "Do NOT publish whichever looks closest -- unlike S6, where the claim is "
            "flatness and an offset is harmless, HERE THE ABSOLUTE VALUES ARE THE "
            "CLAIM. A slope this script invented would contradict the Discussion.")
    cname, swap, slope_a, slope_b, lam_a, lam_b, _werr = winner
    ok(f"G3 convention = '{cname}'" + (" with the arms SWAPPED" if swap else "")
       + f"  ->  A {slope_a:.3f}, B {slope_b:.3f}")
    print(f"\n  >>> The Methods should name this convention: {cname}"
          + (", regressing the comparison arm on MAGNet_matched" if swap else ""))
    if cname == "mean of squared SEs [inverted]" and not swap:
        print("      Same as S6. ONE estimator spec covers both supplements, which is")
        print("      what the deposit README should say.")
    else:
        print("      ⚠ DIFFERENT from S6's estimator. Two supplements, two Deming")
        print("        conventions -- state each explicitly; do not let a reader")
        print("        assume one applies to both.")

    # ---- CIs and G4 -----------------------------------------------------
    print("\n  bootstrapping over events ...")
    ci_a = boot_ci(x, ya, sx, sya, cname, swap, a.boot)
    ci_b = boot_ci(x, yb, sx, syb, cname, swap, a.boot)
    print(f"    arm A  {slope_a:.3f}  95% CI {ci_a[0]:.3f}-{ci_a[1]:.3f}"
          f"   (record {EXPECT_A} , {EXPECT_A_CI[0]}-{EXPECT_A_CI[1]})")
    print(f"    arm B  {slope_b:.3f}  95% CI {ci_b[0]:.3f}-{ci_b[1]:.3f}"
          f"   (record {EXPECT_B} , {EXPECT_B_CI[0]}-{EXPECT_B_CI[1]})")

    if not (np.isfinite(ci_a[1]) and np.isfinite(ci_b[0])):
        die("G4 FAILED: bootstrap did not return usable intervals")
    if ci_a[1] >= ci_b[0]:
        die(f"G4 FAILED: the intervals OVERLAP ({ci_a[1]:.3f} >= {ci_b[0]:.3f}). "
            "Non-overlap is the claim; without it the Discussion's sentence is "
            "not supported and must be rewritten, not rounded.")
    ok(f"G4 intervals do not overlap: A <= {ci_a[1]:.3f} < {ci_b[0]:.3f} <= B")

    product = slope_a * slope_b
    print(f"\n  A x B = {product:.3f}   (record {EXPECT_PRODUCT}; "
          f"observed cross-cohort C = {EXPECT_OBSERVED_C})")
    print("  The multiplicative model is APPROXIMATE -- C involves a third platform,")
    print("  a different read length and a full-size MAGNet arm. Report it as")
    print("  'broadly consistent', never as a prediction.")

    # ---- deposit --------------------------------------------------------
    rows = [
        dict(part="matched-pair panel", arm="A", what="platform only, same hearts",
             comparison="MAGNet_matched -> LBNL_matched", n_events=len(idx),
             slope=round(slope_a, 4), ci_lo=round(ci_a[0], 4),
             ci_hi=round(ci_a[1], 4), same_direction=f"{dir_a}/{len(idx)}",
             lambda_convention=cname + (" [arms swapped]" if swap else ""),
             lambda_value=round(lam_a, 4),
             recorded=EXPECT_A, source="mp_SE_* tables, recomputed"),
        dict(part="matched-pair panel", arm="B", what="people only, same platform",
             comparison="MAGNet_matched -> MAGNet_other", n_events=len(idx),
             slope=round(slope_b, 4), ci_lo=round(ci_b[0], 4),
             ci_hi=round(ci_b[1], 4), same_direction=f"{dir_b}/{len(idx)}",
             lambda_convention=cname + (" [arms swapped]" if swap else ""),
             lambda_value=round(lam_b, 4),
             recorded=EXPECT_B, source="mp_SE_* tables, recomputed"),
        dict(part="matched-pair panel", arm="A x B", what="product",
             comparison="platform x sampling", n_events=len(idx),
             slope=round(product, 4), ci_lo=np.nan, ci_hi=np.nan,
             same_direction="",
             lambda_convention=cname + (" [arms swapped]" if swap else ""),
             lambda_value=np.nan,
             recorded=EXPECT_PRODUCT, source="derived"),
    ]
    for coh, sl in PER_COHORT.items():
        rows.append(dict(part="transfer slope, MAGNet-significant SE events",
                         arm=coh, what="cross-cohort, Figure 4",
                         comparison=f"MAGNet -> {coh}", n_events=np.nan,
                         slope=sl, ci_lo=np.nan, ci_hi=np.nan, same_direction="",
                         lambda_convention="", lambda_value=np.nan, recorded=sl,
                         source="45_cohort_slopes.py (NOT recomputed here)"))
    # The four GENOME-WIDE rows were removed on 20 Oct: they carried
    # pre-exclusion JCEC values from a document that forbids quoting them.
    # S4 now presents ONE slope quantity, which is why the "two different
    # quantities" caveat below is also gone. 11 rows -> 7.
    pd.DataFrame(rows).to_csv(a.out, sep="\t", index=False)

    print("\n" + "=" * 74)
    print(f"S4 BUILT -- {a.out}, {len(rows)} rows")
    print("=" * 74)
    print("  CAVEATS THAT MUST TRAVEL WITH THE TABLE:")
    print(f"   * {EXPECT_COMMON} events, both arms at 10 DCM vs 14 NF. Both slopes are attenuated")
    print("     by regression to the mean; ONLY THEIR COMPARISON is interpretable,")
    print("     never either alone as a magnitude estimate.")
    print("   * 'Platform' bundles library prep, read layout and pipeline batch. The")
    print("     design cannot separate them -- the claim is about changing measurement")
    print("     setup, not poly-A versus total-RNA specifically.")
    print("   * Arm B is ONE deterministic draw of non-matched MAGNet samples. The")
    print("     bootstrap covers event sampling, not draw-to-draw variability. Repeat")
    print("     with 3-5 draws before quoting B to three decimals.")
    print("   * No covariates in any arm, deliberately: MAGNet has age/sex/RIN and")
    print("     LBNL has none, so adjusting one and not the other would confound")
    print("     platform with adjustment.")
    print("   * Selection came from a MAGNet run excluding all 48 samples used in any")
    print("     arm. Without that, the MAGNet-matched arm would be biased relative to")
    print("     the LBNL arm and would MANUFACTURE the asymmetry being measured.")
    print()
    print("  ⚠ THIS TABLE CARRIES ONE SLOPE QUANTITY ONLY: the Figure 4")
    print("    transfer slopes, on MAGNet-significant SE events, post-exclusion,")
    print("    JC, from 45_cohort_slopes.py. The pre-exclusion JCEC genome-wide")
    print("    slopes that once sat beside them were removed on 20 Oct — their")
    print("    source document forbids quoting them. If a reviewer asks why")
    print("    Figure 4 differs from any earlier circulated table, that is the")
    print("    answer, and it belongs in a reply rather than in the supplement.")
    print()


if __name__ == "__main__":
    main()
