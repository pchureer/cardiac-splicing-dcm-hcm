#!/usr/bin/env python
"""
83_build_s6.py -- Supplementary Table S6: the control-group size sweep.

THE RESULT STANDS. THE VERDICT IN THE SOURCE DOC DOES NOT.
----------------------------------------------------------
`RESULT_control_size_sweep.md` measured the thing it set out to measure and the answer
has not changed: **control-group size does not depress the Deming slope.** Flat from 5 to
83 controls, across a 17-fold range.

But that doc is from 12 September and it closes with:

    "Disease stage is now the leading explanation; read length is the only
     alternative still untested."
    | read length / annotation | untested |
    | disease stage            | leading explanation, not yet established |

**Both have since been settled the other way.** Disease stage was REFUTED
(`RESULT_magnitude_resolved.md`: end-stage DCM sits at the same slope as preserved-EF
myectomy cohorts), read length was TESTED (S5(b)), and the actual explanation was
MEASURED -- the platform term, 0.584 against 0.988 at matched group sizes on the same 24
hearts (`RESULT_phase15_platform_vs_sampling.md`).

A supplement asserting "disease stage is leading" alongside a Discussion that refutes it
is the paper disagreeing with itself in print. **This script rebuilds the table and
replaces the verdict.** The sweep's own finding is untouched; only the paragraph that
pointed at what to try next is rewritten, because what was tried next has answers.

⚠ SCRIPT 19 PREDATES 24_lock_event_convention.py
------------------------------------------------
So its slopes may be full-coordinate. Every slope here is computed under BOTH
conventions and the build reports which reproduces the record -- the same approach `62`
took for S5(c), where the answer turned out to matter.

⚠ THE LAMBDA CONVENTION IS NOT RECORDED EITHER
----------------------------------------------
As in `82`, the script searches conventions rather than guessing, and names the winner.
Unlike `82`, a near-miss is survivable here: **S6's claim is FLATNESS across n, not any
absolute slope.** A systematic offset shifts every row equally and leaves the claim
intact. So an unmatched convention gives a loud warning and `--allow-unmatched`, not a
refusal -- but the warning must reach the caption.

GATES
-----
  G1  the file inventory matches the record's 7 failed runs exactly
      (m7_k28 x2, m14 x1, m28 x3, m7_k18 x1 -- 43 of 50 present)
  G2  the reference arm carries 752 significant events
  G3  the m = 5 / 7 / 9 mean slopes reproduce 0.948 / 0.891 / 0.986
  G4  the slope is FLAT -- max minus min across n_ctrl < 0.20, which is what the
      supplement actually claims

USAGE
-----
  conda activate rnaseq_tools
  cd /data/project/chureeratp/RNA/HEART

  python 83_build_s6.py --list
  python 83_build_s6.py
  python 83_build_s6.py --allow-unmatched    # proceed on the best lambda convention
"""

import argparse
import glob
import os
import sys

import numpy as np
import pandas as pd

NKEEP_SE = 3
C_KEY = "key"
C_DPSI_CANDIDATES = ("dPSI_adj", "dPSI_raw", "dPSI")
SE_CANDIDATES = ("se_adj", "se", "se_dpsi")
FDR_MAX, DPSI_MIN = 0.05, 0.10
OUT = "S6.tsv"

M_VALUES = [5, 7, 9, 14, 28, 56, 83]
REPS = [1, 2, 3, 4, 5]
K_ARM = [13, 18]                      # the case-count arm, all at m = 7
SWEEP = "sizesweep/SE_m{m}_k{k}_r{r}.all_events.tsv"
REF_GLOB = "splitctl/SE_DCMvsA_*.all_events.tsv"
# ⚠ THE SWEEP'S OWN OUTPUT TABLE EXISTS. Found 7 Oct, written 11 Sep 22:04 -- the same
# minute the last sweep run finished. It carries per-replicate deming/r/pct_same/n for
# BOTH arms, so S6 is built FROM IT and the recomputation above becomes an independent
# verification that the deposited script reproduces the deposited table. That is a
# stronger position than either alone, and it hands over the program arm for free.
SWEEP_TSV = "control_size_sweep.tsv"

# ------------------------------------------------------------------ the record
EXPECT_REF_SIG = 752
EXPECT_MISSING = {(7, 28, 2), (7, 28, 3), (14, 28, 2),
                  (28, 28, 2), (28, 28, 4), (28, 28, 5), (7, 18, 2)}
EXPECT_SLOPE = {5: 0.948, 7: 0.891, 9: 0.986, 14: 0.947,
                28: 1.011, 56: 1.006, 83: 0.992}
EXPECT_SD = {5: 0.115, 7: 0.115, 9: 0.097, 14: 0.104,
             28: 0.025, 56: 0.028, 83: 0.000}
EXPECT_EVENTS = {5: 357, 7: 378, 9: 353, 14: 212, 28: 83, 56: 114, 83: 114}
EXPECT_DIR = {5: 89.5, 7: 90.1, 9: 93.9, 14: 95.7, 28: 98.8, 56: 98.6, 83: 99.1}
EXPECT_REPS = {5: 5, 7: 3, 9: 5, 14: 4, 28: 2, 56: 5, 83: 5}
# the three cohorts whose own control groups are the matched comparison
COHORT_MATCH = {
    "Brigham": dict(n_ctrl=5, bg=0.662, prog=0.587),
    "Stanford": dict(n_ctrl=7, bg=0.579, prog=0.703),
    "Fuwai": dict(n_ctrl=9, bg=0.328, prog=0.440),
}
FLAT_TOL = 0.20


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


USE_EXON = True          # set by the reference search; the sweep must match it


def load(path, exon=None):
    exon = USE_EXON if exon is None else exon
    if not (os.path.exists(path) and os.path.getsize(path)):
        return None
    t = pd.read_csv(path, sep="\t", low_memory=False)
    d = pick(t.columns, C_DPSI_CANDIDATES)
    if d is None:
        d = next((c for c in t.columns if str(c).lower().startswith("dpsi")), None)
    if d is None or C_KEY not in t.columns:
        return None
    s = pick(t.columns, SE_CANDIDATES)
    t = t.copy()
    t["dPSI"] = pd.to_numeric(t[d], errors="coerce")
    t["SE"] = pd.to_numeric(t[s], errors="coerce") if s else np.nan
    t["FDR"] = pd.to_numeric(t["FDR"], errors="coerce") if "FDR" in t.columns else np.nan
    t["gkey"] = (t[C_KEY].astype(str).str.split(":").str[:NKEEP_SE].str.join(":")
                 if exon else t[C_KEY].astype(str))
    t = t.dropna(subset=["dPSI"])
    srt = "FDR" if t["FDR"].notna().any() else "dPSI"
    return (t.sort_values(srt, kind="mergesort")
             .drop_duplicates(subset=["gkey"], keep="first"))


RULES = {
    "FDR<0.05 & |dPSI|>=0.10 (locked)":
        lambda t: (t.FDR < 0.05) & (t.dPSI.abs() >= 0.10),
    "FDR<0.05 only":
        lambda t: t.FDR < 0.05,
    "FDR<0.05 & |dPSI|>=0.05":
        lambda t: (t.FDR < 0.05) & (t.dPSI.abs() >= 0.05),
}


def ref_grid(path):
    """Every (convention, rule) count for one candidate reference.

    The reference's own identity was NOT recorded -- only "752 significant events" --
    and script 19 predates 24_lock_event_convention.py, so neither the counting
    convention nor the significance rule can be assumed. Searching the grid finds the
    reference by its number, the way 65 found S7's model and 82 found the lambda.
    """
    out = {}
    for exon in (True, False):
        t = load(path, exon=exon)
        if t is None:
            continue
        for rname, rule in RULES.items():
            out[("exon" if exon else "full-coord", rname)] = int(rule(t).sum())
    return out


def deming(x, y, lam):
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


def lam_of(name, sx, sy):
    """lambda for Deming. The INVERTED variants exist because the direction of the
    ratio is a free choice nobody recorded, and it moves the slope a long way.

    Here the sweep runs have FEWER controls and therefore LARGER standard errors than
    the reference, so sx^2/sy^2 < 1 while sy^2/sx^2 > 1 -- and raising lambda lowers
    the fitted slope. A systematic +0.2 offset under every non-inverted convention is
    exactly what an inverted ratio looks like, so both directions are searched rather
    than assumed."""
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
    return np.nan


_BASE = ["mean of squared SEs", "median of squared SEs", "squared mean SE"]
CONVENTIONS = _BASE + [b + " [inverted]" for b in _BASE] + ["orthogonal (lambda = 1)"]


def arm_stats(ref, run, conv, swap=False):
    """Slope / r / direction on the events shared with the reference's significant set."""
    sig = ref[(ref.FDR < FDR_MAX) & (ref.dPSI.abs() >= DPSI_MIN)]
    j = sig.merge(run[["gkey", "dPSI", "SE"]].rename(
        columns={"dPSI": "y", "SE": "sy"}), on="gkey", how="inner")
    if len(j) < 5:
        return None
    xs, ys = (j.y.to_numpy(), j.dPSI.to_numpy()) if swap else \
             (j.dPSI.to_numpy(), j.y.to_numpy())
    sxs, sys = (j.sy.to_numpy(), j.SE.to_numpy()) if swap else \
               (j.SE.to_numpy(), j.sy.to_numpy())
    lam = lam_of(conv, sxs, sys)
    if not np.isfinite(lam):
        return None
    sl = deming(xs, ys, lam)
    r = float(np.corrcoef(j.dPSI, j.y)[0, 1]) if len(j) > 2 else np.nan
    same = float(((j.dPSI > 0) == (j.y > 0)).mean() * 100)
    return dict(n=len(j), slope=sl, r=r, same_dir=same, lam=lam)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=OUT)
    ap.add_argument("--tol", type=float, default=0.03)
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--allow-unmatched", action="store_true")
    ap.add_argument("--ref", default=None,
                    help="reference table, if it is not in splitctl/")
    a = ap.parse_args()

    # ---- inventory ------------------------------------------------------
    present, missing = {}, set()
    for m in M_VALUES:
        for r in REPS:
            p = SWEEP.format(m=m, k=28, r=r)
            (present.setdefault((m, 28), []).append(r) if os.path.exists(p)
             else missing.add((m, 28, r)))
    for k in K_ARM:
        for r in REPS:
            p = SWEEP.format(m=7, k=k, r=r)
            (present.setdefault((7, k), []).append(r) if os.path.exists(p)
             else missing.add((7, k, r)))

    if a.list:
        print("=" * 76)
        print("S6 inputs")
        print("=" * 76)
        for (m, k), rs in sorted(present.items()):
            print(f"  m={m:<3} k={k:<3}  replicates present: {sorted(rs)}")
        print(f"\n  {len(missing)} missing: {sorted(missing)}")
        print(f"  record expects exactly these 7: {sorted(EXPECT_MISSING)}")
        refs = sorted(glob.glob(REF_GLOB))
        print(f"\n  reference candidates ({len(refs)}), searching for "
              f"{EXPECT_REF_SIG} significant across conventions AND rules:")
        hits = []
        for f in refs:
            g = ref_grid(f)
            print(f"\n    {f}")
            for (convn, rname), n in g.items():
                mark = f"   <<< {EXPECT_REF_SIG}, THE REFERENCE" \
                    if n == EXPECT_REF_SIG else ""
                print(f"      {convn:<11} {rname:<34} {n:>6}{mark}")
                if n == EXPECT_REF_SIG:
                    hits.append((f, convn, rname))
        if hits:
            print(f"\n  {len(hits)} match(es). The sweep will be computed under the")
            print("  same convention, because a slope against a reference collapsed")
            print("  one way and runs collapsed the other compares different rows.")
        else:
            print(f"\n  ⚠ NO combination gives {EXPECT_REF_SIG}. The reference is not")
            print("    in splitctl/ under any convention or rule tried. Pass it with")
            print("    --ref, or check whether script 19 wrote its own reference")
            print("    table somewhere else. Do NOT pick the closest.")
        return

    print("=" * 76)
    print("SUPPLEMENTARY TABLE S6 -- control-group size sweep")
    print("=" * 76)

    # ---- G1 -------------------------------------------------------------
    if missing != EXPECT_MISSING:
        extra = missing - EXPECT_MISSING
        absent = EXPECT_MISSING - missing
        msg = []
        if extra:
            msg.append(f"unexpectedly missing: {sorted(extra)}")
        if absent:
            msg.append(f"expected missing but present: {sorted(absent)}")
        die("G1 FAILED: the run inventory does not match the record's 7 failed runs. "
            + "; ".join(msg) + ". The record says the failures were caused by "
            "--min-samples = min(m, k) demanding coverage in ALL controls; if the "
            "inventory differs, this is not the same sweep.")
    ok(f"G1 inventory matches: {50 - len(missing)} of 50 runs, the 7 known failures")

    # ---- G2: find the reference ------------------------------------------
    global USE_EXON
    cands = [a.ref] if a.ref else sorted(glob.glob(REF_GLOB))
    ref, refpath, refconv, refrule = None, None, None, None
    best = []
    for f in cands:
        for (convn, rname), n in ref_grid(f).items():
            best.append((abs(n - EXPECT_REF_SIG), f, convn, rname, n))
            # `refpath is None`, not `ref is None` -- ref is assigned after the
            # loop, so the old test was always true and this kept the LAST match
            # instead of the first. With RULES ordered locked-first and exon-first,
            # taking the first match prefers the paper's own convention and rule
            # whenever several reproduce the count.
            if n == EXPECT_REF_SIG and refpath is None:
                refpath, refconv, refrule = f, convn, rname
    if refpath is None:
        best.sort()
        print("\n  closest combinations found:")
        for d, f, convn, rname, n in best[:5]:
            print(f"    {n:>6}  ({d:+d})  {os.path.basename(f)}  "
                  f"{convn}, {rname}")
        die(f"G2 FAILED: no candidate carries {EXPECT_REF_SIG} significant events "
            "under any convention or rule. Run --list for the full grid. Do NOT "
            "substitute the closest -- every slope below would then be against a "
            "baseline the record never used, which is exactly how "
            "RESULT_magnitude_resolved.md came to contradict Figure 4 by 0.281. "
            "Pass the right table with --ref if it lives outside splitctl/.")
    USE_EXON = (refconv == "exon")
    ref = load(refpath, exon=USE_EXON)
    ok(f"G2 reference = {refpath}")
    print(f"       {EXPECT_REF_SIG} significant under: {refconv}, {refrule}")
    if not USE_EXON:
        print("       ⚠ FULL-COORDINATE. Script 19 predates the convention lock, so")
        print("         the whole sweep is recomputed full-coordinate to match the")
        print("         reference. The CAPTION must say S6 is full-coordinate while")
        print("         every count in the main text is exon-level -- two conventions")
        print("         in one paper is survivable only if each is labelled.")
    if refrule != "FDR<0.05 & |dPSI|>=0.10 (locked)":
        print(f"       ⚠ The reference was selected under '{refrule}', NOT the locked")
        print("         rule. Say so in the caption.")

    # ---- G3: the lambda convention ---------------------------------------
    print("\n" + "-" * 76)
    print("  LAMBDA CONVENTION SEARCH -- which reproduces 0.948 / 0.891 / 0.986?")
    print("-" * 76)
    runs = {m: [load(SWEEP.format(m=m, k=28, r=r)) for r in present.get((m, 28), [])]
            for m in (5, 7, 9)}
    best, best_err, best_swap = None, np.inf, False
    print(f"  {'x axis':<10}{'convention':<30}" +
          "".join(f"{m:>9}" for m in (5, 7, 9)) + "      err")
    for swap in (False, True):
        for conv in CONVENTIONS:
            means = {}
            for m in (5, 7, 9):
                sl = [arm_stats(ref, t, conv, swap) for t in runs[m] if t is not None]
                sl = [v["slope"] for v in sl if v and np.isfinite(v["slope"])]
                means[m] = float(np.mean(sl)) if sl else np.nan
            err = np.nansum([abs(means[m] - EXPECT_SLOPE[m]) for m in (5, 7, 9)])
            hit = all(np.isfinite(means[m])
                      and abs(means[m] - EXPECT_SLOPE[m]) <= a.tol for m in (5, 7, 9))
            print(f"  {'run' if swap else 'reference':<10}{conv:<30}"
                  + "".join(f"{means[m]:>9.3f}" for m in (5, 7, 9))
                  + f"{err:>9.3f}" + ("   <<<" if hit else ""))
            if err < best_err:
                best, best_err, best_swap = conv, err, swap
    matched = best_err <= a.tol * 3
    SWAP = best_swap
    if SWAP:
        print("\n  >>> The reference is the Y axis, not the X axis. Script 19")
        print("      regressed the reference ON the sweep run. State that in the")
        print("      caption -- a Deming slope is not symmetric and the reciprocal")
        print("      of 0.948 is 1.055, not 0.948.")
    if not matched:
        print(f"\n  ⚠ NO convention reproduces all three within +/-{a.tol}. "
              f"Best: '{best}' (total error {best_err:.3f}).")
        print("    S6's claim is FLATNESS, not an absolute slope, and a systematic")
        print("    offset shifts every row equally -- so this is survivable. But the")
        print("    CAPTION MUST SAY the slopes are recomputed and differ from script")
        print("    19's by up to "
              f"{best_err / 3:.3f} on average.")
        if not a.allow_unmatched:
            die("G3 FAILED. Re-run with --allow-unmatched to proceed on "
                f"'{best}', with that caption note mandatory.")
    else:
        ok(f"G3 convention = '{best}'")
    conv = best

    # ---- build the sweep table -------------------------------------------
    print("\n" + "-" * 76)
    print("  THE SWEEP -- does control-group size depress the slope?")
    print("-" * 76)
    print("  Deming here is the BACKGROUND arm (all events), recomputed. The PROGRAM")
    print("  arm is read from control_size_sweep.tsv below, so both are now in S6 --")
    print("  but the column still has to say WHICH, because the source doc quotes")
    print("  both and a bare 'Deming' header conflates them.")
    print(f"  {'n_ctrl':>7}{'reps':>6}{'events':>9}{'Deming (bg)':>18}"
          f"{'r':>8}{'same dir':>10}{'record':>9}")
    rows, slopes = [], {}
    for m in M_VALUES:
        st = [arm_stats(ref, load(SWEEP.format(m=m, k=28, r=r)), conv, SWAP)
              for r in present.get((m, 28), [])]
        st = [s for s in st if s]
        if not st:
            print(f"  {m:>7}{0:>6}{'-':>9}{'-':>18}{'-':>8}{'-':>10}"
                  f"{EXPECT_SLOPE[m]:>9.3f}")
            continue
        sl = np.array([s["slope"] for s in st], float)
        rows.append(dict(part="control-size sweep (background arm)", n_ctrl=m, n_cases=28,
                         replicates=len(st), events=int(np.mean([s["n"] for s in st])),
                         deming_mean=round(float(np.nanmean(sl)), 4),
                         deming_sd=round(float(np.nanstd(sl, ddof=1)), 4)
                         if len(sl) > 1 else 0.0,
                         pearson_r=round(float(np.mean([s["r"] for s in st])), 4),
                         pct_same_direction=round(
                             float(np.mean([s["same_dir"] for s in st])), 2),
                         recorded_deming=EXPECT_SLOPE[m],
                         recorded_sd=EXPECT_SD[m],
                         recorded_events=EXPECT_EVENTS[m],
                         recorded_pct_same_direction=EXPECT_DIR[m],
                         lambda_convention=conv,
                         x_axis="run" if SWAP else "reference"))
        slopes[m] = float(np.nanmean(sl))
        print(f"  {m:>7}{len(st):>6}{int(np.mean([s['n'] for s in st])):>9,}"
              f"{np.nanmean(sl):>12.3f} +/-"
              f"{(np.nanstd(sl, ddof=1) if len(sl) > 1 else 0.0):<5.3f}"
              f"{np.mean([s['r'] for s in st]):>8.3f}"
              f"{np.mean([s['same_dir'] for s in st]):>9.1f}%"
              f"{EXPECT_SLOPE[m]:>9.3f}")
        if len(st) != EXPECT_REPS[m]:
            print(f"          [warn] {len(st)} usable replicates, record says "
                  f"{EXPECT_REPS[m]}")

    # case-count arm
    print(f"\n  case-count arm at n_ctrl = 7:")
    for k in K_ARM + [28]:
        st = [arm_stats(ref, load(SWEEP.format(m=7, k=k, r=r)), conv, SWAP)
              for r in present.get((7, k), [])]
        st = [s for s in st if s]
        if not st:
            continue
        sl = np.nanmean([s["slope"] for s in st])
        print(f"    k = {k:<4} {len(st)} reps   Deming {sl:.3f}")
        rows.append(dict(part="case-count arm (background)", n_ctrl=7, n_cases=k,
                         replicates=len(st),
                         events=int(np.mean([s["n"] for s in st])),
                         deming_mean=round(float(sl), 4),
                         deming_sd=round(float(np.nanstd(
                             [s["slope"] for s in st], ddof=1)), 4)
                         if len(st) > 1 else 0.0,
                         pearson_r=np.nan, pct_same_direction=np.nan,
                         recorded_deming=np.nan, recorded_sd=np.nan,
                         recorded_events=np.nan,
                         recorded_pct_same_direction=np.nan,
                         lambda_convention=conv))

    # ---- G4: flatness, which is the actual claim -------------------------
    if len(slopes) >= 4:
        spread = max(slopes.values()) - min(slopes.values())
        print(f"\n  slope spread across n_ctrl {min(slopes)}-{max(slopes)}: "
              f"{spread:.3f}")
        if spread >= FLAT_TOL:
            die(f"G4 FAILED: the slope is NOT flat (spread {spread:.3f} >= "
                f"{FLAT_TOL}). S6's entire claim is flatness. If this fires, the "
                "supplement's conclusion is wrong, not its formatting.")
        ok(f"G4 slope is flat across a {max(slopes) / min(slopes):.0f}-fold range "
           f"of n (spread {spread:.3f})")

    # ---- the comparison that decides it ----------------------------------
    print("\n" + "-" * 76)
    print("  THE COMPARISON THAT DECIDES IT -- matched control size")
    print("-" * 76)
    print(f"  {'cohort':<10}{'n_ctrl':>7}{'observed bg':>13}{'observed prog':>15}"
          f"{'MAGNet here':>14}")
    for coh, d in COHORT_MATCH.items():
        mine = slopes.get(d["n_ctrl"], np.nan)
        print(f"  {coh:<10}{d['n_ctrl']:>7}{d['bg']:>13.3f}{d['prog']:>15.3f}"
              f"{mine:>14.3f}")
        rows.append(dict(part="matched control size", n_ctrl=d["n_ctrl"],
                         n_cases=np.nan, replicates=np.nan, events=np.nan,
                         deming_mean=round(mine, 4) if np.isfinite(mine) else np.nan,
                         deming_sd=np.nan, pearson_r=np.nan,
                         pct_same_direction=np.nan,
                         recorded_deming=d["bg"], recorded_sd=np.nan,
                         recorded_events=np.nan,
                         recorded_pct_same_direction=np.nan,
                         lambda_convention=f"{coh} observed bg/prog "
                                           f"{d['bg']}/{d['prog']}"))
    print("\n  At matched control size the gap is ~0.94 against ~0.52.")

    # ---- the sweep's OWN output table: authoritative source + cross-check ------
    file_rows = []
    if os.path.exists(SWEEP_TSV):
        print("\n" + "-" * 76)
        print(f"  {SWEEP_TSV} -- the sweep's own output, read as the SOURCE")
        print("-" * 76)
        st = pd.read_csv(SWEEP_TSV, sep=r"\s+")
        need = {"m", "k", "rep", "set", "n", "deming", "r", "pct_same"}
        if not need.issubset(st.columns):
            print(f"  [warn] unexpected columns {list(st.columns)}; skipped")
        else:
            print(f"  {len(st)} per-replicate rows, arms: "
                  f"{sorted(st['set'].unique())}")
            print(f"\n  {'n_ctrl':>7}{'k':>5}{'arm':>12}{'reps':>6}"
                  f"{'deming':>10}{'sd':>8}{'events':>8}{'record':>9}{'recomp':>9}")
            bad = []
            for (m, k, arm), g in st.groupby(["m", "k", "set"], sort=True):
                mu = float(g.deming.mean())
                sd = float(g.deming.std(ddof=1)) if len(g) > 1 else 0.0
                rec = EXPECT_SLOPE.get(m) if (k == 28 and arm == "background") else None
                mine = slopes.get(m) if (k == 28 and arm == "background") else None
                print(f"  {m:>7}{k:>5}{arm:>12}{len(g):>6}{mu:>10.3f}{sd:>8.3f}"
                      f"{g.n.mean():>8.0f}"
                      + (f"{rec:>9.3f}" if rec is not None else f"{'-':>9}")
                      + (f"{mine:>9.3f}" if mine is not None else f"{'-':>9}"))
                if rec is not None and abs(mu - rec) > 0.001:
                    bad.append(f"n_ctrl {m}: file {mu:.4f} vs record {rec}")
                if mine is not None and abs(mu - mine) > 0.01:
                    bad.append(f"n_ctrl {m}: file {mu:.4f} vs recomputed {mine:.4f}")
                file_rows.append(dict(
                    part=f"control-size sweep ({arm} arm, from {SWEEP_TSV})",
                    n_ctrl=m, n_cases=k, replicates=len(g),
                    events=int(round(g.n.mean())),
                    deming_mean=round(mu, 4), deming_sd=round(sd, 4),
                    pearson_r=round(float(g.r.mean()), 4),
                    pct_same_direction=round(float(g.pct_same.mean()), 2),
                    recorded_deming=rec if rec is not None else np.nan,
                    recorded_sd=EXPECT_SD.get(m, np.nan) if rec is not None else np.nan,
                    recorded_events=EXPECT_EVENTS.get(m, np.nan)
                    if rec is not None else np.nan,
                    recorded_pct_same_direction=EXPECT_DIR.get(m, np.nan)
                    if rec is not None else np.nan,
                    lambda_convention=f"as computed by 19 (see {SWEEP_TSV})",
                    x_axis="reference"))
            if bad:
                print("\n  ⚠ DISAGREEMENTS:")
                for b in bad:
                    print("     " + b)
                die("G5 FAILED: the sweep's own table disagrees with the record or "
                    "with the recomputation. The file is the source of truth -- find "
                    "out which of the other two is wrong before depositing anything.")
            ok(f"G5 {SWEEP_TSV} agrees with the record AND with the recomputation")
            print("\n  >>> S6 is now sourced from the deposited table and verified by")
            print("      an independent recomputation that reproduces it. The program")
            print("      arm comes from the file -- it no longer needs the 'not")
            print("      recomputed' caveat.")
    else:
        print(f"\n  [warn] {SWEEP_TSV} not found -- S6 rests on the recomputation "
              "alone, which reproduces the record but is not the original output.")

    pd.DataFrame(rows + file_rows).to_csv(a.out, sep="\t", index=False)

    # ---- the verdict, rewritten -----------------------------------------
    print("\n" + "=" * 76)
    print(f"S6 BUILT -- {a.out}, {len(rows)} rows")
    print("=" * 76)
    print("  THE VERDICT FOR THE CAPTION -- replacing the source doc's, which is stale")
    print("  " + "-" * 74)
    print("  Control-group size does not depress the transfer slope: it is flat across")
    print("  a 17-fold range of n, and at matched control size MAGNet returns ~0.94")
    print("  where the external cohorts return ~0.52. The errors-in-variables")
    print("  mechanism proposed in RESULT_split_control_check.md is negligible here")
    print("  and is withdrawn.")
    print()
    print("  Concordance DOES degrade with small n, as theory predicts -- 99.1% at")
    print("  n = 83 down to 89.5% at n = 5 -- which partly closes the gap with")
    print("  Brigham (86.3%) and Stanford (83.5%) but not Fuwai (73.7%). So noise")
    print("  accounts for some of the direction-concordance drop and none of the")
    print("  slope drop.")
    print()
    print("  ⚠ DO NOT CARRY OVER the source doc's closing paragraph. It says disease")
    print("     stage is the leading explanation and read length is untested. Since")
    print("     then: disease stage REFUTED (end-stage DCM sits at the same slope as")
    print("     the preserved-EF myectomy cohorts); read length TESTED (S5b); and the")
    print("     explanation MEASURED -- the platform term, 0.584 against 0.988 at")
    print("     matched group sizes on the same 24 hearts (S4). S6 excludes one")
    print("     candidate; S4 names the one that survives. Point the caption at S4.")
    print()
    print("  CAVEATS THAT MUST TRAVEL WITH THE TABLE:")
    print("   * 7 of 50 runs returned 0 rows. The decision rows at n = 5 / 7 / 9")
    print("     retain 5 / 3 / 5 replicates.")
    print("   * DESIGN FLAW in script 19: --min-samples was set to min(m, k), so the")
    print("     filter demands coverage in ALL controls at small m -- which is what")
    print("     killed those 7 runs AND means the event set changes across the sweep")
    print("     (357 events at m = 5 against 114 at m = 83). The slope is flat")
    print("     DESPITE that, not because the comparison was held fixed.")
    print("   * MAGNet's disjoint slope is THRESHOLD-DEPENDENT: 0.872 at")
    print("     min-samples 15 (637 events) against 0.992 here at min-samples 28")
    print("     (114 events). **Quote it as a range, 0.87-0.99, never a point.**")
    print("     Neither end is near the external cohorts' 0.44-0.70.")
    print("   * Requiring 5/5 controls selects the best-covered events, whereas")
    print("     Brigham's run required 4 of 5. Unlikely to move a slope by 0.3, but")
    print("     a per-cohort threshold-matched pass would close it.")
    if not matched:
        print(f"   * ⚠ Slopes recomputed under '{conv}'; script 19's convention is not")
        print("     recorded and these differ from its published values. The flatness")
        print("     claim is unaffected; the absolute values should not be quoted")
        print("     against script 19's.")
    print()


if __name__ == "__main__":
    main()
