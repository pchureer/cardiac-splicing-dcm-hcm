#!/usr/bin/env python3
"""
158_fig3_panels_bc.py -- A7: Figure 3's panel B and panel C at 125 bp.

v2, rewritten against the real schema. The first version refused on every
input, and the refusal was correct four times over:

  1. The event key column is `key`, seven colon-separated fields
     (chr:exonStart:exonEnd:upstreamStart:upstreamEnd:downstreamStart:
     downstreamEnd for SE). v1 collapsed on the WHOLE key, which is not a
     collapse at all — 100,014 rows stayed 100,014 events instead of becoming
     62,971. **NKEEP counts FIELDS OF THE KEY, not rows.** SE keeps 3
     (chr + the skipped exon's own two coordinates), MXE/A3SS/A5SS 5, RI 3.
     That is what "collapse to the event's own coordinates" means: the
     flanking exons vary, the event does not.
  2. dPSI is `dPSI_adj` (with `dPSI_raw` beside it), not `dpsi`.
  3. `mde80` is ALREADY A COLUMN. 2.802 x 0.006799149323365336 =
     0.019051216404069673 = column 19 exactly, so panel B uses the file's own
     published quantity and only cross-checks the constant.
  4. Colorado was never re-quantified, so there is NO `_rl125` Colorado table.
     The same file serves both read lengths. v1's suffix filter excluded it
     and then reported it missing.

THE FOUR NUMBERS
----------------
[108], verbatim:

  (B) 95% of the 19,998 DCM-versus-ICM events could have detected the median
      DCM-versus-non-failing effect (|dPSI| = 0.16).
      ^ the producer is 37_figure3_etiology.py:panel_b. Both of panel B's
        inputs are COLORADO (DCM-vs-NF for the median, DCM-vs-ICM for the
        ECDF), so panel B is READ-LENGTH INVARIANT. And the ECDF runs over
        the ~476 MATCHED events, not over 19,998 — the legend's denominator
        is wrong at both read lengths, which is a defect to fix, not a
        number to re-base.
  (C) The two disease contrasts against the shared non-failing group agree
      closely (n = 770 discovery-significant events measurable in both;
      Deming 0.82, 91% same direction).

An earlier draft of this file asserted that panel B's threshold was a MAGNet
quantity and therefore moved. It is not: the median comes from Colorado's own
DCM-vs-NF arm (MAGNet's 926 give 0.1318, which rounds to 0.13, not 0.16).
Panel C is the part that moves — all three of its numbers are MAGNet.

IDENTIFICATION
--------------
By content, never by name — but now with the tested count as well as the
significant count, which is a far stronger fingerprint:

    MAGNet DCM vs NF    62,971 tested /   926 sig (50 bp)
                        53,456 tested /   813 sig (125 bp)
    MAGNet HCM vs NF         ? tested /   468 sig (50 bp); 125 bp unrecorded
    Colorado DCM vs ICM 19,998 tested /     2 sig (both; not re-quantified)

Candidates are TOP-LEVEL tables only. v1 globbed subdirectories and swept in
`sizesweep/`, `splitctl/`, `inflate/`, `frozen_50bp/` and
`jccheck_50bp_BACKUP/` — 250 files of scratch, sweeps and backups, any of
which could have matched a count by coincidence.

USAGE
-----
  python3 158_fig3_panels_bc.py --schema      # columns + collapse arithmetic
  python3 158_fig3_panels_bc.py               # calibrate 50 bp, then 125 bp
  python3 158_fig3_panels_bc.py --self-test
"""
import argparse
import glob
import os
import re
import sys

try:
    import numpy as np
    import pandas as pd
except ImportError:
    sys.exit("needs numpy and pandas (conda activate rnaseq_tools)")

# ----------------------------------------------------------------- the record
NKEEP = {"SE": 3, "MXE": 5, "A3SS": 5, "A5SS": 5, "RI": 3}
FDR_MAX, DPSI_MIN = 0.05, 0.10
MDE80_K = 2.802
MDE80_TOL = 1e-9

# (tested, significant); None = unrecorded at that read length
FINGERPRINT = {
    "MAGNet DCM vs NF": {"": (62971, 926), "_rl125": (53456, 813)},
    "MAGNet HCM vs NF": {"": (None, 468), "_rl125": (None, None)},
    "Colorado DCM vs ICM": {"": (19998, 2), "_rl125": (19998, 2)},
    "Colorado DCM vs NF": {"": (19953, 476), "_rl125": (19953, 476)},
}
# Colorado is read-length-invariant: not re-quantified, so one file, both runs
INVARIANT = {"Colorado DCM vs ICM", "Colorado DCM vs NF"}

PUB = {"median_dpsi": 0.16, "power_frac": 0.95, "n_both": 770,
       "deming": 0.82, "same_dir": 0.91}

DPSI_NAMES = ["dPSI_adj", "dPSI_raw"]          # adjusted first: the paper's model
SE_NAMES = ["se_adj"]
FDR_NAMES = ["FDR", "fdr"]
KEY_NAMES = ["key", "ekey", "event_key"]
TYPE_NAMES = ["event_type", "etype", "type"]

DENY = ("sizesweep", "splitctl", "inflate", "frozen_50bp", "BACKUP",
        "tmp", "abandon", "probe", "mp_", "sweep")


def _col(df, names, what, required=True):
    for n in names:
        if n in df.columns:
            return n
    if not required:
        return None
    raise KeyError(f"no {what} column among {names}; columns are "
                   f"{list(df.columns)}")


def load(path, dpsi_col=None):
    df = pd.read_csv(path, sep="\t", low_memory=False)
    ren = {_col(df, KEY_NAMES, "event-key"): "key",
           _col(df, [dpsi_col] if dpsi_col else DPSI_NAMES, "dPSI"): "dpsi",
           _col(df, FDR_NAMES, "FDR"): "fdr",
           _col(df, SE_NAMES, "standard-error"): "se_adj"}
    t = _col(df, TYPE_NAMES, "event-type", required=False)
    if t:
        ren[t] = "etype"
    m = _col(df, ["mde80"], "mde80", required=False)
    if m:
        ren[m] = "mde80"
    out = df.rename(columns=ren)
    if "etype" not in out.columns:
        out["etype"] = _etype_from_name(path)
    return out


def _etype_from_name(path):
    b = os.path.basename(path)
    for et in ("A3SS", "A5SS", "MXE", "SE", "RI"):
        if re.search(rf"(^|[^A-Z]){et}([^A-Z]|$)", b):
            return et
    raise KeyError(f"no event_type column and none inferable from {b}")


def event_key(df):
    """Collapse key = the first NKEEP fields of `key`.

    NKEEP counts FIELDS, not rows. For SE that is chr + the skipped exon's own
    start and end; the four flanking coordinates are dropped, which is what
    makes one event out of its several flanking variants. Getting this wrong
    is not a small error: it is the difference between 100,014 "events" and
    62,971.
    """
    bad = set(df.etype.unique()) - set(NKEEP)
    if bad:
        raise KeyError(f"no NKEEP for event types {sorted(bad)}")
    nk = df.etype.map(NKEEP)
    return [":".join(str(k).split(":")[:n]) for k, n in zip(df.key, nk)]


def collapse_then_filter(df):
    """COLLAPSE-FIRST: lowest-FDR row per event, then FDR<0.05 & |dPSI|>=0.10."""
    d = df.assign(_ek=event_key(df))
    c = d.sort_values("fdr", kind="mergesort").drop_duplicates("_ek",
                                                               keep="first")
    sig = c[(c.fdr < FDR_MAX) & (c.dpsi.abs() >= DPSI_MIN)]
    return c, sig


def filter_then_collapse(df):
    """Reported only as the WRONG answer, so the 991-vs-926 gap stays visible."""
    d = df.assign(_ek=event_key(df))
    f = d[(d.fdr < FDR_MAX) & (d.dpsi.abs() >= DPSI_MIN)]
    return f.sort_values("fdr", kind="mergesort").drop_duplicates("_ek")


def check_mde80(df, path):
    """The file's mde80 must equal 2.802 x se_adj, or the constant has moved."""
    if "mde80" not in df.columns:
        return f"{os.path.basename(path)}: no mde80 column; using {MDE80_K} x se_adj"
    dev = float((df.mde80 - MDE80_K * df.se_adj).abs().max())
    if dev > MDE80_TOL:
        return (f"{os.path.basename(path)}: mde80 != {MDE80_K} x se_adj "
                f"(max deviation {dev:.3g}) — the power constant in the file "
                f"is not [042]'s")
    return None


# -------------------------------------------------------- 45's deming, reused
def load_deming(src="45_cohort_slopes.py"):
    """Extract and exec ONLY deming() from 45 — never reimplement it.

    45's and 33's agree to four decimals on five cohorts. A third
    implementation would be an unvalidated estimator wearing a published
    number's name, so if the extraction fails this refuses.
    """
    if not os.path.exists(src):
        return None, f"{src} not found — run from the HEART directory"
    text = open(src, errors="replace").read()
    # `|\Z` matters: without it a deming() that is the LAST block in the file
    # is not matched at all, and the loader reports "no top-level def deming"
    # for a file that plainly has one.
    m = re.search(r"^def deming\(.*?(?=^\S|\Z)", text, re.M | re.S)
    if not m:
        return None, f"no top-level `def deming(` in {src}"
    ns = {"np": np, "numpy": np, "math": __import__("math")}
    try:
        exec(compile(m.group(0), f"{src}:deming", "exec"), ns)
    except Exception as e:                                   # noqa: BLE001
        return None, f"could not exec 45's deming(): {e}"
    fn = ns.get("deming")
    if not callable(fn):
        return None, "extraction produced no callable"
    return fn, f"{src} line {text[:m.start()].count(chr(10)) + 1}"


def call_deming(fn, x, y, sx, sy):
    for args in ((x, y, sx, sy), (x, y)):
        try:
            r = fn(*args)
        except TypeError:
            continue
        return float(r[0]) if isinstance(r, (tuple, list, np.ndarray)) else float(r)
    raise TypeError("45's deming() accepted neither (x,y,sx,sy) nor (x,y)")


# ------------------------------------------------------------------ discovery
def candidates(dpsi_col):
    """Top-level SE event tables only. No subdirectories: v1 swept 250 files
    of sweeps, splits and backups, any of which could match a count."""
    rows, seen = [], set()
    for p in sorted(glob.glob("*all_events*.tsv")):
        if any(d in p for d in DENY) or not re.match(r"^(SE_|.*_SE_)", p):
            continue
        rp = os.path.realpath(p)
        if rp in seen:
            continue
        seen.add(rp)
        try:
            df = load(p, dpsi_col)
            col, sig = collapse_then_filter(df)
            rows.append([p, len(df), len(col), len(sig), ""])
        except Exception as e:                               # noqa: BLE001
            rows.append([p, None, None, None, str(e)[:90]])
    return rows


def identify(rows, suffix, ref_files=None):
    out, problems, notes = {}, [], []
    for lab, fps in FINGERPRINT.items():
        want_t, want_s = fps[suffix]
        if want_s is None:
            continue
        hits = [r for r in rows if r[3] == want_s
                and (want_t is None or r[2] == want_t)]
        if lab not in INVARIANT:
            is125 = [("_rl125" in r[0]) for r in hits]
            hits = [h for h, i in zip(hits, is125) if i == (suffix == "_rl125")]
        if len(hits) == 1:
            out[lab] = hits[0][0]
        elif not hits:
            problems.append(
                f"{lab}: no top-level table has {want_s} significant"
                + (f" out of {want_t:,} collapsed" if want_t else "")
                + (" (read-length invariant: expected the SAME file as 50 bp)"
                   if lab in INVARIANT else ""))
        else:
            # A TIE. Break it only by consistency with the 50 bp file that
            # was identified BY FINGERPRINT — never by which name looks
            # right. SE_rl125 and SE_v2filtered_rl125 both carry
            # 53,456/813, and only one of them is SE_v2filtered's twin.
            pick = None
            if ref_files and lab in ref_files:
                twins = set(suffixed_names(ref_files[lab], suffix))
                same = [h[0] for h in hits if h[0] in twins]
                if len(same) == 1:
                    pick = same[0]
            if pick:
                out[lab] = pick
                notes.append(
                    f"{lab}: {len(hits)} tables carry the same fingerprint "
                    f"({', '.join(h[0] for h in hits)}); took {pick} as the "
                    f"suffixed twin of {ref_files[lab]}, which was itself "
                    f"identified by fingerprint")
            else:
                problems.append(f"{lab}: {len(hits)} tables match — "
                                + ", ".join(h[0] for h in hits[:4])
                                + ("; none is the suffixed twin of "
                                   f"{ref_files[lab]}" if ref_files
                                   and lab in ref_files else ""))
    return out, problems, notes


def suffixed_names(path, suffix):
    """The candidate names of `path`'s suffixed twin. Name only, no I/O."""
    d, b = os.path.dirname(path), os.path.basename(path)
    return [os.path.join(d, b.replace(".all_events.tsv",
                                      f"{suffix}.all_events.tsv")),
            os.path.join(d, re.sub(r"\.tsv$", f"{suffix}.tsv", b))]


def find_renamed(rows, suffix, want_name, deny_name):
    """The 125 bp HCM-vs-NF table was RENAMED, not suffixed.

    At 50 bp it is `SE_HCM_v2`; at 125 bp `SE_HCMvsNF_rl125`. No suffix
    derivation finds that, so this matches on the CONTRAST in the name and
    refuses on ambiguity. It is the weakest identification in the script and
    is labelled as such in the output — Rule 8 in its original form: a rename
    is only half an edit.
    """
    hits = [r for r in rows if suffix in r[0] and want_name in r[0]
            and deny_name not in r[0] and r[2] is not None]
    if not hits:
        return None, (f"no {suffix} table whose name contains '{want_name}' "
                      f"(and not '{deny_name}')")
    if len(hits) > 1:
        return None, (f"{len(hits)} candidates named '{want_name}' at "
                      f"{suffix}: " + ", ".join(h[0] for h in hits))
    return hits[0][0], None


def derive_suffixed(path, suffix):
    """Carry a 50 bp identification across for the ONE arm with no recorded
    125 bp count. A name-derived input is Rule 8's mirror, so it is checked
    for existence, uniqueness and a plausible size — and labelled in output."""
    d, b = os.path.dirname(path), os.path.basename(path)
    cands = []
    for nb in (b.replace(".all_events.tsv", f"{suffix}.all_events.tsv"),
               re.sub(r"\.tsv$", f"{suffix}.tsv", b)):
        p = os.path.join(d, nb)
        if p != path and os.path.exists(p) and p not in cands:
            cands.append(p)
    if not cands:
        return None, f"no {suffix} twin of {path} exists"
    if len(cands) > 1:
        return None, f"{len(cands)} possible {suffix} twins: " + ", ".join(cands)
    p = cands[0]
    try:
        ref = len(collapse_then_filter(load(path))[0])
        got = len(collapse_then_filter(load(p))[0])
    except Exception as e:                                   # noqa: BLE001
        return None, f"{p} unreadable: {str(e)[:70]}"
    if ref and not 0.5 * ref <= got <= 1.5 * ref:
        return None, (f"{p} has {got:,} collapsed events against {ref:,} at "
                      f"50 bp — too far apart to be the same contrast")
    return p, None


# -------------------------------------------------------------- the two panels
def panel_b(cnf_sig, cicm_col, path):
    """EXACTLY 37_figure3_etiology.py:panel_b, lines 167-192.

        sig = D[(D.FDR < 0.05) & (D.dpsi.abs() >= 0.10)]   # Colorado DCM vs NF
        m   = E[E.ekey.isin(set(sig.ekey))]                # Colorado DCM vs ICM
        v   = np.sort(m.mde80.dropna().values)
        med = float(np.median(sig.dpsi.abs()))
        pct = 100 * (v <= med).mean()

    Three things this script had wrong before reading the producer:

    * the median is **Colorado's own** DCM-vs-NF effect, not MAGNet's
      (MAGNet's 926 give 0.1318, which rounds to 0.13, not 0.16);
    * the ECDF is over the DCM-vs-ICM events **matched to that significant
      set** — about 476 — and NOT over all 19,998;
    * `mde80.dropna()` drops events with no standard error, so the
      denominator is the matched-and-measurable count, smaller again.

    The legend says "95% of the 19,998 DCM-versus-ICM events". The figure
    computes 95% of the matched subset. That clause names the wrong
    denominator at BOTH read lengths — a defect the re-quantification did
    not cause and does not fix.
    """
    med = float(cnf_sig.dpsi.abs().median())
    m = cicm_col[cicm_col._ek.isin(set(cnf_sig._ek))]
    mde = (m.mde80 if "mde80" in m.columns else MDE80_K * m.se_adj).astype(float)
    v = mde.dropna()
    frac = float((v <= med).mean()) if len(v) else float("nan")
    # the denominator the legend claims, reported alongside so the gap is visible
    all_mde = (cicm_col.mde80 if "mde80" in cicm_col.columns
               else MDE80_K * cicm_col.se_adj).astype(float).dropna()
    frac_legend = float((all_mde <= med).mean()) if len(all_mde) else float("nan")
    return (med, frac, len(v), len(cicm_col), frac_legend,
            check_mde80(cicm_col, path))


def panel_c(disc_sig, other_col, deming_fn):
    both = disc_sig.merge(other_col[["_ek", "dpsi", "se_adj"]], on="_ek",
                          suffixes=("", "_o"))
    n = len(both)
    if not n:
        return 0, float("nan"), float("nan")
    same = float((np.sign(both.dpsi) == np.sign(both.dpsi_o)).mean())
    slope = call_deming(deming_fn,
                        both.dpsi.values.astype(float),
                        both.dpsi_o.values.astype(float),
                        both.se_adj.values.astype(float),
                        both.se_adj_o.values.astype(float))
    return n, slope, same


def run(suffix, deming_fn, dpsi_col, calibrate, rows, ref_files=None):
    tag = "125 bp" if suffix == "_rl125" else "50 bp"
    print(f"\n{'=' * 74}\n{tag}   (suffix {suffix or '<none>'}, dPSI = "
          f"{dpsi_col or DPSI_NAMES[0]})\n{'=' * 74}")
    files, problems, notes = identify(rows, suffix, ref_files)
    derived = {}
    if ref_files:
        for lab, fps in FINGERPRINT.items():
            if fps[suffix][1] is None and lab in ref_files and lab not in files:
                p, err = derive_suffixed(ref_files[lab], suffix)
                if err and lab == "MAGNet HCM vs NF":
                    # renamed at 125 bp: SE_HCM_v2 -> SE_HCMvsNF_rl125
                    p, err2 = find_renamed(rows, suffix, "HCMvsNF", "DCMvsHCM")
                    err = None if p else f"{err}; and {err2}"
                if err:
                    problems.append(f"{lab}: {err}")
                else:
                    files[lab] = p
                    derived[lab] = p
                    ref_n = next((r[2] for r in rows
                                  if r[0] == ref_files[lab]), None)
                    got_n = next((r[2] for r in rows if r[0] == p), None)
                    if ref_n and got_n and not 0.6 * ref_n <= got_n <= ref_n:
                        problems.append(
                            f"{lab}: {p} collapses to {got_n:,} against "
                            f"{ref_n:,} at 50 bp — outside the 0.85 shrinkage "
                            f"every other table shows; wrong file?")
    missing = [k for k in FINGERPRINT if k not in files]
    if missing and not problems:
        problems.append("not identified and not derivable: " + ", ".join(missing))
    for lab in sorted(files):
        mark = "  [name-derived, NOT fingerprinted]" if lab in derived else ""
        inv = "  [read-length invariant]" if lab in INVARIANT else ""
        print(f"  {lab:<22} {files[lab]}{inv}{mark}")
    for nt in notes:
        print(f"  ⚠ tie broken — {nt}")
    if problems:
        print("\nREFUSING — inputs not identified by fingerprint:")
        for q in problems:
            print("  * " + q)
        print("\n  Candidates and what they collapse to:")
        for p, nraw, ncol, nsig, err in rows:
            print(f"    {p:<46} rows {str(nraw):>7} collapsed {str(ncol):>7}"
                  f" sig {str(nsig):>5}  {err}")
        print("\n  Fingerprints are content, not names. Do not rename a file")
        print("  to make this pass — a count that cannot be reproduced means")
        print("  the collapse convention or the dPSI column differs.")
        return None

    disc = load(files["MAGNet DCM vs NF"], dpsi_col)
    disc_col, disc_sig = collapse_then_filter(disc)
    other_col, other_sig = collapse_then_filter(
        load(files["MAGNet HCM vs NF"], dpsi_col))
    col_col, _ = collapse_then_filter(load(files["Colorado DCM vs ICM"], dpsi_col))
    _, cnf_sig = collapse_then_filter(load(files["Colorado DCM vs NF"], dpsi_col))

    (med, frac, nmatch, ntested, frac_legend,
     mde_note) = panel_b(cnf_sig, col_col, files["Colorado DCM vs ICM"])
    n, slope, same = panel_c(disc_sig, other_col, deming_fn)
    ff = len(filter_then_collapse(disc))

    print(f"\n  discovery          rows {len(disc):>7,}  collapsed "
          f"{len(disc_col):>7,}  sig {len(disc_sig):>5}")
    print(f"                     (filter-first would give {ff} — not used)")
    if mde_note:
        print(f"  ⚠ {mde_note}")
    print(f"  panel B  median |dPSI|      {med:>8.4f}   (Colorado DCM vs NF, "
          f"{len(cnf_sig)} sig)")
    print(f"           matched+measurable {nmatch:>8}   of {ntested:,} "
          f"DCM-vs-ICM events")
    print(f"           mde80 <= median    {frac:>8.4f}   ({100 * frac:.1f}%)"
          f"   <- the figure")
    print(f"           same over all {ntested:,}  {frac_legend:>8.4f}   "
          f"({100 * frac_legend:.1f}%)   <- what the legend's wording claims")
    print(f"  panel C  n measurable both  {n:>8}")
    print(f"           Deming slope       {slope:>8.4f}")
    print(f"           same direction     {same:>8.4f}   ({100 * same:.1f}%)")
    print(f"  [A] HCM vs NF sig          {len(other_sig):>8}"
          f"   — panel A's \"the same 28 hypertrophic hearts return 468\"")

    got = {"median_dpsi": med, "power_frac": frac, "n_both": n,
           "deming": slope, "same_dir": same, "hcm_sig": len(other_sig),
           "disc_sig": len(disc_sig)}
    if calibrate:
        # every gate reads PUB. An earlier draft of this loop hard-coded the
        # five values with PUB sitting unused beside it — a check and its
        # target in two places, which is how 156's verify() came to refuse
        # the patch it had just written.
        bad = []
        for k in ("median_dpsi", "power_frac", "deming", "same_dir"):
            if round(got[k], 2) != PUB[k]:
                bad.append(f"{k} {got[k]:.4f} -> {round(got[k], 2)}, "
                           f"published {PUB[k]}")
        if n != PUB["n_both"]:
            bad.append(f"n measurable in both {n}, published {PUB['n_both']}")
        if bad:
            print("\n✗ CALIBRATION FAILED — the 50 bp record is not reproduced:")
            for b in bad:
                print("    " + b)
            alt = [c for c in DPSI_NAMES if c != (dpsi_col or DPSI_NAMES[0])]
            print(f"\n  No 125 bp number is believable until this passes.")
            print(f"  Try --dpsi-col {alt[0]} if the dPSI column is the issue;")
            print("  otherwise panel C's 'measurable in both' is a different")
            print("  set than significant-and-present, or the HCM-vs-NF table")
            print("  identified by 468 is not the one the figure used.")
            return None
        print("\n✓ CALIBRATION PASSED — all five published 50 bp values "
              "reproduced")
    return got


# ------------------------------------------------------------------ self-test
def _mk(rows):
    """rows: (key, etype, dpsi, fdr, se)"""
    d = pd.DataFrame(rows, columns=["key", "etype", "dpsi", "fdr", "se_adj"])
    d["mde80"] = MDE80_K * d.se_adj
    return d


def self_test():
    fails = []

    # THE v1 BUG: NKEEP counts fields of the key, not rows. Three rows that
    # are one SE event (same chr:start:end, different flanking coordinates)
    # must collapse to ONE.
    df = _mk([("chr1:100:200:1:2:3:4", "SE", 0.30, 0.02, 0.05),
              ("chr1:100:200:5:6:7:8", "SE", 0.02, 0.01, 0.05),
              ("chr1:100:200:9:9:9:9", "SE", 0.40, 0.03, 0.05),
              ("chr1:500:600:1:2:3:4", "SE", 0.40, 0.01, 0.05)])
    col, sig = collapse_then_filter(df)
    if len(col) != 2:
        fails.append(f"collapse on 3 key fields gave {len(col)} events, want 2 "
                     f"— NKEEP is being applied to rows, not key fields")
    # and the lowest-FDR row of the collapsed event fails the effect filter,
    # so that event drops: collapse-first 1, filter-first 2
    if len(sig) != 1 or len(filter_then_collapse(df)) != 2:
        fails.append(f"collapse/filter order not distinguished: "
                     f"{len(sig)} vs {len(filter_then_collapse(df))} (want 1, 2)")

    # MXE keeps 5 fields, so two rows differing in field 6 are ONE MXE event
    # but TWO under SE's rule — the per-type key length must be honoured
    mx = _mk([("chr1:1:2:3:4:9:9", "MXE", 0.3, 0.01, 0.05),
              ("chr1:1:2:3:4:8:8", "MXE", 0.3, 0.02, 0.05)])
    if len(collapse_then_filter(mx)[0]) != 1:
        fails.append("MXE did not collapse on 5 fields")
    se_same = _mk([("chr1:1:2:3:4:9:9", "SE", 0.3, 0.01, 0.05),
                   ("chr1:1:2:8:8:8:8", "SE", 0.3, 0.02, 0.05)])
    if len(collapse_then_filter(se_same)[0]) != 1:
        fails.append("SE did not collapse on 3 fields")
    if len(collapse_then_filter(_mk([
            ("chr1:1:2:3:4:9:9", "MXE", 0.3, 0.01, 0.05),
            ("chr1:1:2:8:8:8:8", "MXE", 0.3, 0.02, 0.05)]))[0]) != 2:
        fails.append("MXE collapsed two distinct events into one")

    # an unknown event type must raise, not silently pick a key length
    try:
        collapse_then_filter(_mk([("c:1:2:3:4:5:6", "XX", 0.3, 0.01, 0.05)]))
        fails.append("an unknown event type did not raise")
    except KeyError:
        pass

    # mde80 integrity: the real constant passes, a wrong one is caught
    ok = _mk([("c:1:2:3:4:5:6", "SE", 0.2, 0.01, 0.006799149323365336)])
    if check_mde80(ok, "x") is not None:
        fails.append("a correct mde80 column was rejected")
    bad = ok.copy()
    bad["mde80"] = 2.0 * bad.se_adj
    if "not [042]'s" not in (check_mde80(bad, "x") or ""):
        fails.append("an mde80 built from the wrong constant was accepted")

    # panel B must restrict the ECDF to the MATCHED events (37's line 173),
    # not take it over all of them. The fixture is built so the two
    # denominators give DIFFERENT answers — 50% matched against 95% overall —
    # so a regression to the old reading cannot pass.
    cicm = _mk([("c:1:%d:3:4:5:6" % i, "SE", 0.0, 1.0, 0.01)
                for i in range(95)]
               + [("c:9:%d:3:4:5:6" % i, "SE", 0.0, 1.0, 0.50)
                  for i in range(5)])
    cicm_col, _ = collapse_then_filter(cicm)
    # two significant DCM-vs-NF events: one cheap to detect, one expensive
    cnf = _mk([("c:1:0:9:9:9:9", "SE", 0.16, 0.01, 0.01),
               ("c:9:0:9:9:9:9", "SE", -0.16, 0.01, 0.50)])
    _, cnf_sig = collapse_then_filter(cnf)
    med, frac, nmatch, nt, frac_leg, note = panel_b(cnf_sig, cicm_col, "x")
    if abs(med - 0.16) > 1e-12:
        fails.append(f"panel B median {med}, want 0.16")
    if nmatch != 2 or abs(frac - 0.5) > 1e-12:
        fails.append(f"panel B matched {nmatch} events, frac {frac} "
                     f"(want 2 and 0.5) — the ECDF is not restricted to the "
                     f"DCM-vs-NF significant set")
    if nt != 100 or abs(frac_leg - 0.95) > 1e-12:
        fails.append(f"the legend-denominator figure is {frac_leg} over {nt} "
                     f"(want 0.95 over 100)")
    if abs(frac - frac_leg) < 1e-9:
        fails.append("the matched and legend denominators gave the same "
                     "answer; the fixture no longer distinguishes them")
    # <= is inclusive: the claim is "could have detected"
    edge_c = _mk([("c:1:2:3:4:5:6", "SE", 0.0, 1.0, 0.16 / MDE80_K)])
    edge_col, _ = collapse_then_filter(edge_c)
    edge_sig = collapse_then_filter(
        _mk([("c:1:2:9:9:9:9", "SE", 0.16, 0.01, 0.05)]))[1]
    if panel_b(edge_sig, edge_col, "x")[1] != 1.0:
        fails.append("an event whose mde80 exactly equals the median was "
                     "excluded; the claim is 'could have detected'")
    # mde80 NaN must drop from the denominator, as 37's .dropna() does
    nan_c = _mk([("c:1:2:3:4:5:6", "SE", 0.0, 1.0, 0.01),
                 ("c:3:4:3:4:5:6", "SE", 0.0, 1.0, 0.01)])
    nan_c.loc[1, "mde80"] = np.nan
    nan_col, _ = collapse_then_filter(nan_c)
    nan_sig = collapse_then_filter(
        _mk([("c:1:2:9:9:9:9", "SE", 0.16, 0.01, 0.05),
             ("c:3:4:9:9:9:9", "SE", 0.16, 0.01, 0.05)]))[1]
    if panel_b(nan_sig, nan_col, "x")[2] != 1:
        fails.append("an event with NaN mde80 stayed in the denominator")

    # panel C joins on the collapsed key, not the raw key
    a = _mk([("c:1:2:3:4:5:6", "SE", 0.2, 0.01, 0.05),
             ("c:3:4:3:4:5:6", "SE", -0.2, 0.01, 0.05),
             ("c:5:6:3:4:5:6", "SE", 0.4, 0.01, 0.05)])
    b = _mk([("c:1:2:9:9:9:9", "SE", 0.1, 0.01, 0.05),
             ("c:3:4:9:9:9:9", "SE", 0.1, 0.01, 0.05)])
    _, asig = collapse_then_filter(a)
    bcol, _ = collapse_then_filter(b)
    n, _, sd = panel_c(asig, bcol,
                       lambda x, y, sx, sy: (float(np.sum(x * y) / np.sum(x * x)),))
    if n != 2:
        fails.append(f"panel C joined {n} events, want 2 (different flanking "
                     f"coordinates must still match)")
    if abs(sd - 0.5) > 1e-12:
        fails.append(f"same-direction {sd}, want 0.5")

    # fingerprints: tested AND significant, invariance, ties
    rows = [["SE_v2filtered.all_events.tsv", 100014, 62971, 926, ""],
            ["SE_HCM_v2.all_events.tsv", 70003, 48024, 468, ""],
            ["SE_Colorado_DCMvsICM_dedup.all_events.tsv", 24288, 19998, 2, ""],
            ["SE_Colorado_DCMvsNF_dedup.all_events.tsv", 24314, 19953, 476, ""],
            # the plain Colorado tables must NOT be picked: 20,027/0 and
            # 20,590/460 are the un-deduplicated counts, and 37's
            # load_events() prefers _dedup
            ["SE_Colorado_DCMvsICM.all_events.tsv", 24336, 20027, 0, ""],
            ["SE_Colorado_DCMvsNF.all_events.tsv", 25135, 20590, 460, ""]]
    f, p, _ = identify(rows, "")
    if f.get("Colorado DCM vs NF") != "SE_Colorado_DCMvsNF_dedup.all_events.tsv":
        fails.append(f"the plain Colorado DCM-vs-NF table was picked: "
                     f"{f.get('Colorado DCM vs NF')}")
    if p or len(f) != 4:
        fails.append(f"clean 50 bp identification failed: {p}")
    # Colorado must be found at 125 bp from the SAME unsuffixed file
    f, p, _ = identify(rows + [["SE_v2filtered_rl125.all_events.tsv",
                             90000, 53456, 813, ""]], "_rl125")
    if "Colorado DCM vs ICM" not in f:
        fails.append("the read-length-invariant Colorado file was not reused "
                     "at 125 bp")
    if f.get("MAGNet DCM vs NF") != "SE_v2filtered_rl125.all_events.tsv":
        fails.append("the 125 bp discovery table was not the suffixed one")
    # a table with the right sig count but the wrong tested count is rejected
    f, p, _ = identify([["decoy.all_events.tsv", 1, 999, 926, ""]] + rows[1:], "")
    if not any("MAGNet DCM vs NF" in x for x in p):
        fails.append("926 significant was accepted without 62,971 collapsed")
    # a tie with NO 50 bp anchor must refuse
    tie = rows + [["SE_other.all_events.tsv", 100014, 62971, 926, ""]]
    f, p, _ = identify(tie, "")
    if not any("2 tables match" in x for x in p):
        fails.append("a tie was not refused")

    # a tie at 125 bp IS broken, but only by the suffixed twin of the file
    # fingerprinted at 50 bp — the real case: SE_rl125 and
    # SE_v2filtered_rl125 both carry 53,456/813
    r125 = rows + [["SE_rl125.all_events.tsv", 80973, 53456, 813, ""],
                   ["SE_v2filtered_rl125.all_events.tsv", 80973, 53456, 813, ""]]
    ref = {"MAGNet DCM vs NF": "SE_v2filtered.all_events.tsv"}
    f, p, nt = identify(r125, "_rl125", ref_files=ref)
    if f.get("MAGNet DCM vs NF") != "SE_v2filtered_rl125.all_events.tsv":
        fails.append(f"the 125 bp tie was not broken to the suffixed twin: "
                     f"{f.get('MAGNet DCM vs NF')}")
    if not any("tables carry the same fingerprint" in x for x in nt):
        fails.append("the tie-break was not reported as a note")
    # and with NO anchor, the same tie must refuse rather than pick
    f, p, nt = identify(r125, "_rl125", ref_files=None)
    if "MAGNet DCM vs NF" in f:
        fails.append("the 125 bp tie was broken without a 50 bp anchor")
    # an anchor whose twin is not among the tied files must also refuse
    f, p, nt = identify(r125, "_rl125",
                        ref_files={"MAGNet DCM vs NF":
                                   "SE_v2cmfrac.all_events.tsv"})
    if "MAGNet DCM vs NF" in f:
        fails.append("the tie was broken by a file that is not the anchor's "
                     "twin")
    if not any("none is the suffixed twin" in x for x in p):
        fails.append("the refusal did not say why the tie could not be broken")

    # the deming loader refuses rather than substituting
    fn, note = load_deming("/nonexistent_45.py")
    if fn is not None or "not found" not in note:
        fails.append("a missing 45_cohort_slopes.py did not refuse")

    # it must also find a deming() that is the LAST block in the file, with
    # nothing after it, and the extracted function must actually run
    import tempfile
    with tempfile.TemporaryDirectory() as td:
        for tail in ("", "\ndef after():\n    pass\n"):
            p = os.path.join(td, "f.py")
            open(p, "w").write(
                "import numpy as np\n"
                "def deming(x, y, sx, sy):\n"
                "    return float(np.sum(x * y) / np.sum(x * x)), 0.0\n" + tail)
            g, nt = load_deming(p)
            if g is None:
                fails.append(f"deming at end-of-file not found (tail={tail!r})"
                             f": {nt}")
                continue
            v = call_deming(g, np.array([1.0, 2.0]), np.array([2.0, 4.0]),
                            np.array([1.0, 1.0]), np.array([1.0, 1.0]))
            if abs(v - 2.0) > 1e-12:
                fails.append(f"extracted deming returned {v}, want 2.0")

    # the denylist keeps sweeps, splits and backups out
    for d in ("sizesweep/SE_x.all_events.tsv",
              "jccheck_50bp_BACKUP/SE_v2filtered_JC.all_events.tsv",
              "frozen_50bp/SE_v2filtered.all_events.tsv"):
        if not any(k in d for k in DENY):
            fails.append(f"{d} would not be excluded")

    # the calibration gate must read PUB, not a second copy of the numbers
    src = open(__file__, errors="replace").read() if "__file__" in dir() else ""
    if src:
        gate = src.split("# every gate reads PUB")[-1].split("if bad:")[0]
        for v in ("0.16", "0.95", "0.82", "0.91", "770"):
            if v in gate:
                fails.append(f"the calibration gate hard-codes {v} instead of "
                             f"reading PUB — one number, two places")

    if FINGERPRINT["MAGNet HCM vs NF"]["_rl125"][1] is not None:
        fails.append("the 125 bp HCM-vs-NF count is now recorded — "
                     "fingerprint it instead of deriving the name")

    print("158 self-test")
    for f_ in fails:
        print("  ✗ " + f_)
    if not fails:
        print("  ✓ NKEEP applied to KEY FIELDS (SE 3, MXE 5, per type),")
        print("    collapse-before-filter distinguished, unknown type raises,")
        print("    mde80 constant verified both ways, panel B uses the file's")
        print("    own mde80 with an inclusive <=, panel C joins on the")
        print("    collapsed key, fingerprints need tested AND significant,")
        print("    Colorado reused at 125 bp, denylist holds, and a tie is")
        print("    broken ONLY by the suffixed twin of a fingerprinted file")
    return 1 if fails else 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--schema", action="store_true")
    ap.add_argument("--self-test", action="store_true")
    ap.add_argument("--dpsi-col", default=None,
                    help=f"default {DPSI_NAMES[0]}; try {DPSI_NAMES[1]}")
    ap.add_argument("--deming-src", default="45_cohort_slopes.py")
    a = ap.parse_args()
    if a.self_test:
        return self_test()
    if self_test():
        print("REFUSING: 158's own self-test failed.")
        return 3
    print()

    rows = candidates(a.dpsi_col)
    if not rows:
        print("no top-level SE event tables found — is this the HEART "
              "directory?")
        return 3
    if a.schema:
        print("top-level SE event tables, with the collapse arithmetic")
        print(f"  NKEEP = {NKEEP}  (fields of `key`, not rows)")
        print(f"  dPSI column = {a.dpsi_col or DPSI_NAMES[0]}\n")
        for p, nraw, ncol, nsig, err in rows:
            print(f"  {p:<48} rows {str(nraw):>7} collapsed {str(ncol):>7}"
                  f" sig {str(nsig):>5}  {err}")
        return 0

    fn, note = load_deming(a.deming_src)
    if fn is None:
        print(f"REFUSING: {note}")
        print("  158 will not substitute its own Deming for 45's.")
        return 3
    print(f"[deming] reusing 45's implementation ({note})")

    ref_files, _, _ = identify(rows, "")
    got50 = run("", fn, a.dpsi_col, True, rows)
    if got50 is None:
        return 1
    got125 = run("_rl125", fn, a.dpsi_col, False, rows, ref_files=ref_files)
    if got125 is None:
        return 1

    print(f"\n{'=' * 74}\nA7 — the edits for [108]\n{'=' * 74}")
    print(f"  {'':<24}{'50 bp':>12}{'125 bp':>12}")
    for k, label, fmt in (("median_dpsi", "median |dPSI|", "{:.2f}"),
                          ("power_frac", "could have detected", "{:.0%}"),
                          ("n_both", "n measurable in both", "{:.0f}"),
                          ("deming", "Deming slope", "{:.2f}"),
                          ("same_dir", "same direction", "{:.0%}"),
                          ("disc_sig", "[A] DCM vs NF sig", "{:.0f}"),
                          ("hcm_sig", "[A] HCM vs NF sig", "{:.0f}")):
        print(f"  {label:<24}{fmt.format(got50[k]):>12}"
              f"{fmt.format(got125[k]):>12}")
    print("\n  PANEL B IS READ-LENGTH INVARIANT: both its inputs are")
    print("  Colorado, which was never re-quantified, so the 50 bp and")
    print("  125 bp columns above MUST be identical. If they are not, the")
    print("  wrong file was identified for one of the two Colorado arms.")
    print("\n  But panel B's legend is still wrong, at both read lengths:")
    print("  37_figure3_etiology.py:173 restricts the ECDF to the events")
    print("  matched to the DCM-vs-NF significant set, so the percentage is")
    print("  of ~476 events, not of 19,998. Compare the two rows printed")
    print("  above. Fixing that clause is an edit the re-quantification did")
    print("  not cause and does not fix.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
