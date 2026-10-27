#!/usr/bin/env python3
"""
98_build_replication_table.py — Table 4: replication by DIRECTION, not by count.

WHY THIS TABLE EXISTS
---------------------
Table 3's gate G8 prints, every run, the reason the three hypertrophic cohorts are not
in it: their own significant counts (Brigham 4, Fuwai 169, Stanford 8 skipped exons) are
not comparable with the discovery cohort's 926, and printing them side by side reads as
a failed replication. It is not one -- they replicate by DIRECTION at 69-81%.

But the figures carry those numbers and **nothing in the paper tabulates them.** A
reviewer who wants the replication evidence has to read it off Figure 4. This builds the
table that should sit beside Table 3.

    ⚠ THE RULE THIS TABLE MUST NOT BREAK. This paper forbids comparing effect-size
      magnitude or event counts across cohorts, on its own evidence: PSI is normalised
      by an inclusion-form length that depends on read length, and cohorts differ in
      tested set, depth, covariate model and power. So this table reports DIRECTION
      CONCORDANCE and a transfer SLOPE, and it does NOT carry a column of per-cohort
      significant-event counts. Adding one would contradict the Methods.

WHAT IS IN IT
  One row per independent replication cohort: group sizes (post-exclusion), the number of
  discovery-significant skipped exons measurable in that cohort, direction concordance
  with a Wilson 95% interval, a one-sided exact binomial P against 50%, and the Deming
  transfer slope. Then a pooled row.

USAGE
    python 98_build_replication_table.py --list            # find the inputs, change nothing
    python 98_build_replication_table.py --apply

GATES -- all against numbers recorded independently of this script
  G1  the input files exist and their columns resolve.
  G2  per-cohort concordance reproduces Colorado 375/471, Brigham 359/442,
      Stanford 408/536, Fuwai 494/714. ⚠ Three of those four NUMERATORS were derived
      from recorded percentages rather than read from a file (see DERIVED below), so
      this gate is what turns them into measurements. A mismatch is a refusal.
  G3  the pooled figure reproduces 1,636/2,163 = 75.6% AND the recomputed one-sided
      binomial P reproduces the recorded 8.4e-132 to one significant figure.
  G4  the Deming slopes reproduce 0.521 / 0.476 / 0.386 / 0.243, mean 0.406.
  G5  LBNL is ABSENT. It is a library-preparation replicate of discovery -- 24 of its 33
      samples are the same individuals -- and its 0.445 / 75.9% must never appear in a
      replication table.
  G6  no column of per-cohort significant-event counts is present in the output.
  G7  group sizes match Table 1 (post-exclusion). ⚠ Stanford is 13/6 here, not the
      13/7 recorded in RESULT_three_hcm_cohorts.md -- see STANFORD below.

DERIVED, AND WHY G2 MATTERS
  `PHASE05_edit_list.md` records Colorado as 375/471 explicitly. For the other three it
  records only percentages (81.2, 76.1, 69.2) and event counts (442, 536, 714). The
  numerators 359, 408 and 494 were obtained as round(n x pct) -- and they are almost
  certainly right, because four independently rounded percentages reconstruct BOTH
  recorded totals exactly (2,163 events, 1,636 concordant) and the recomputed pooled P
  lands on the recorded 8.4e-132. That is strong, but it is still arithmetic on rounded
  inputs. **G2 reads the real numerators from the file and refuses if they differ.**

STANFORD'S GROUP SIZE -- a conflict in the record, resolved
  `RESULT_three_hcm_cohorts.md` gives Stanford as 13 HCM / 7 donor. Table 3's G8 and
  Table 1 give 13 / 6. Both are right at different stages: `NUMBERS_CANONICAL.md`
  records three within-cohort replicate pairs (two Colorado, one Stanford) with one
  member of each retained, and the Stanford one was a control -- so 7 pre-exclusion and
  6 post. **This table is post-exclusion, like Tables 1 and 3.** G7 enforces it.
"""

import argparse
import glob
import os
import sys
from math import sqrt

import pandas as pd

try:
    from scipy.stats import binomtest
except ImportError:
    from scipy.stats import binom_test as _bt

    def binomtest(k, n, p, alternative):
        class R:
            pvalue = _bt(k, n, p, alternative=alternative)
        return R()

# Read length. Drives the output names, the input candidates, and the record.
SUF = os.environ.get("TBL_SUFFIX", "")
OUT_TSV = f"Table4{SUF}.tsv"
OUT_MD = f"Table4{SUF}.md"

# ---------------------------------------------------------------------- the record
# Independent of this script. Every one of these is a gate, not an input.
RECORD = {
    # cohort:    (contrast,      accession,     n,       k,   n_ev, deming)
    "Colorado": ("DCM vs NF",   "PRJNA477855", "36/13", 375,  471, 0.521),
    "Brigham":  ("HCM vs ctrl", "PRJNA675090", "18/5",  359,  442, 0.476),
    "Stanford": ("HCM vs ctrl", "PRJNA747806", "13/6",  408,  536, 0.386),
    "Fuwai":    ("HCM vs ctrl", "PRJNA522931", "28/9",  494,  714, 0.243),
}
POOLED_K, POOLED_N = 1636, 2163
POOLED_P_1SF = 8e-132          # 50 bp; display text and the gate value
                               # both come from REC_POOLED_P below
MEAN_SLOPE = 0.406

# --------------------------------------------------------------------------
# 125 bp RECORD, placed AFTER every literal it rewrites.
#
# These are gates, not inputs (98's own header says so). A 50 bp record against
# a 125 bp file does not make a wrong table — it makes G2/G3/G4 refuse a
# CORRECT one. So when TBL_SUFFIX is set, every value is required explicitly.
_REQ = ("REC_KN", "REC_SLOPE", "REC_POOLED", "REC_POOLED_P",
        "REC_MEAN_SLOPE", "REC_LBNL")
_env = {k: os.environ.get(k, "").strip() for k in _REQ}
if SUF:
    _miss = [k for k in _REQ if not _env[k]]
    if _miss:
        raise SystemExit(
            "REFUSING — TBL_SUFFIX=%s but these are unset: %s\n\n"
            "  RECORD is an independent gate. A 50 bp record read against a\n"
            "  125 bp file refuses a CORRECT table. Supply all six:\n\n"
            "    %s\n" % (SUF, ", ".join(_miss), 'TBL_SUFFIX=_rl125 \\\n    REC_KN=Colorado=311/385,Brigham=298/364,Stanford=345/455,Fuwai=422/596 \\\n    REC_SLOPE=Colorado=0.5999,Brigham=0.5099,Stanford=0.4008,Fuwai=0.2948 \\\n    REC_POOLED=1376,1800 \\\n    REC_POOLED_P=2.51e-117 \\\n    REC_MEAN_SLOPE=0.451 \\\n    REC_LBNL=0.5726/78.2 \\\n    python3 98_build_replication_table.py'))

if _env["REC_KN"] or _env["REC_SLOPE"]:
    def _parse(s):
        return dict(p.split("=", 1) for p in s.split(",")) if s else {}
    _kn, _sl = _parse(_env["REC_KN"]), _parse(_env["REC_SLOPE"])
    for _nm, _d in (("REC_KN", _kn), ("REC_SLOPE", _sl)):
        if _d and set(_d) != set(RECORD):
            raise SystemExit(
                "REFUSING — %s must name exactly %s, got %s"
                % (_nm, sorted(RECORD), sorted(_d)))
    for _c in list(RECORD):
        _contrast, _acc, _n, _k, _ne, _s = RECORD[_c]
        if _c in _kn:
            _a, _b = _kn[_c].split("/")
            _k, _ne = int(_a), int(_b)
        if _c in _sl:
            _s = float(_sl[_c])
        RECORD[_c] = (_contrast, _acc, _n, _k, _ne, _s)

if _env["REC_POOLED"]:
    POOLED_K, POOLED_N = (int(_x) for _x in _env["REC_POOLED"].split(","))

# ONE variable drives both the gate value and the two display sites. The
# recorded text is more precise than the 1-s.f. constant the factor-of-two
# gate needs, so deriving the text from the float would silently change the
# 50 bp output from "8.4e-132" to "8.0e-132".
POOLED_P_TEXT = _env["REC_POOLED_P"] or "8.4e-132"
if _env["REC_POOLED_P"]:
    POOLED_P_1SF = float(_env["REC_POOLED_P"])

if _env["REC_MEAN_SLOPE"]:
    MEAN_SLOPE = float(_env["REC_MEAN_SLOPE"])

# LBNL's recorded slope / concordance, quoted in the G5 note. Not a gate.
LBNL_SLOPE, LBNL_PCT = (_env["REC_LBNL"] or "0.445/75.9").split("/")

# At 50 bp only Colorado's numerator was recorded; Brigham, Stanford and Fuwai
# were back-derived from a rounded percentage. At 125 bp all four are COUNTED
# (33 section D, recounted independently by 45 and by 61), so the caveat does
# not apply and the label should not carry it.
SRC_LABEL = {_c: ("measured" if SUF else
                  ("recorded" if _c == "Colorado" else "DERIVED"))
             for _c in RECORD}
# --------------------------------------------------------------------------
SLOPE_TOL = 0.0015
FORBIDDEN = ("sig", "significant", "n_sig", "count")   # G6

CAND_CONC = [f"cohort_type_concordance_dedup{SUF}.tsv",
             f"cohort_type_concordance{SUF}.tsv"]
# _dedup FIRST IN BOTH LISTS. Before this patch CAND_SLOPE preferred the
# NON-dedup `cohort_slopes.tsv`, so had that file existed, 98 would have read
# slopes filter-first and concordances collapse-first — the 926-vs-991 error,
# silently, inside one table. Neither non-dedup file exists (checked 16 Oct),
# so the defect was latent; the order now makes it unreachable.
CAND_SLOPE = [f"cohort_slopes_dedup{SUF}.tsv",
              f"cohort_slopes{SUF}.tsv",
              f"Figure4_slopes{SUF}.tsv"]


def hr(c="-"):
    print(c * 76)


def refuse(msg):
    hr("=")
    print("REFUSING — nothing was written.\n")
    print(msg)
    hr("=")
    sys.exit(2)


def wilson(k, n, z=1.959963985):
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return 100 * (c - h), 100 * (c + h)


def one_sided(k, n):
    return binomtest(k, n, 0.5, alternative="greater").pvalue


def find(d, cands, what):
    for c in cands:
        p = os.path.join(d, c)
        if os.path.exists(p):
            return p
    hits = sorted(glob.glob(os.path.join(d, "*concord*")) +
                  glob.glob(os.path.join(d, "*slope*")))
    refuse(f"cannot find the {what} file.\n\n  Looked for: " +
           ", ".join(cands) +
           "\n\n  Files here that might be it:\n" +
           ("".join(f"    {os.path.basename(h)}\n" for h in hits) or "    (none)\n") +
           "\n  Re-run with the right name added to the candidate list, or pass -d.")


def pick(df, words, what, where, optional=False):
    """⚠ EXACT names first. The real concordance file's event-count column is named `n`,
    and "n" is a substring of "event_type" — a substring-only search picked the wrong
    column. Any single-letter or very short candidate must be matched exactly."""
    low = {str(c).lower(): c for c in df.columns}
    for w in words:
        if w in low:
            return low[w]
    for c in df.columns:
        if any(len(w) > 2 and w in str(c).lower() for w in words):
            return c
    if optional:
        return None
    refuse(f"cannot resolve the {what} column in {where}.\n  columns: " +
           ", ".join(map(str, df.columns)) +
           f"\n\n  candidates tried: {', '.join(words)}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-d", "--dir", default=".")
    ap.add_argument("--list", dest="list_only", action="store_true")
    ap.add_argument("--apply", action="store_true")
    a = ap.parse_args()

    hr("=")
    print("98 — Table 4: replication by direction")
    hr("=")

    # ------------------------------------------------------------------ G1
    cpath = find(a.dir, CAND_CONC, "concordance")
    spath = find(a.dir, CAND_SLOPE, "slope")
    conc = pd.read_csv(cpath, sep="\t")
    slop = pd.read_csv(spath, sep="\t")
    print(f"\n  concordance: {os.path.basename(cpath)}  "
          f"{len(conc)} rows x {len(conc.columns)} cols")
    print(f"    columns: {', '.join(map(str, conc.columns))}")
    print(f"  slopes:      {os.path.basename(spath)}  "
          f"{len(slop)} rows x {len(slop.columns)} cols")
    print(f"    columns: {', '.join(map(str, slop.columns))}")

    c_coh = pick(conc, ["cohort"], "cohort", os.path.basename(cpath))
    c_typ = pick(conc, ["type", "class", "event"], "event type", os.path.basename(cpath))
    c_k = pick(conc, ["same", "concord", "agree", "n_same"], "concordant count",
               os.path.basename(cpath))
    c_n = pick(conc, ["n", "n_ev", "n_tot", "n_match", "total", "measur"],
               "event count", os.path.basename(cpath))
    c_pct = pick(conc, ["pct", "percent", "concordance"], "percentage",
                 os.path.basename(cpath), optional=True)
    s_coh = pick(slop, ["cohort"], "cohort", os.path.basename(spath))
    s_sl = pick(slop, ["deming", "deming_slope", "slope"], "slope",
                os.path.basename(spath))
    s_r = pick(slop, ["pearson_r", "pearson", "r"], "Pearson r",
               os.path.basename(spath), optional=True)
    s_n = pick(slop, ["n"], "event count", os.path.basename(spath), optional=True)
    s_ind = pick(slop, ["independent"], "independence flag",
                 os.path.basename(spath), optional=True)
    print(f"\n  resolved: cohort={c_coh} type={c_typ} concordant={c_k} events={c_n}"
          f" pct={c_pct}")
    print(f"            slope={s_sl} pearson={s_r} slope_n={s_n} independent={s_ind}")
    if s_r:
        print("  ✅ `pearson_r` is present — per-cohort r goes in the table, closing an")
        print("     open item that the record held only as the range 0.36-0.59 at 50 bp / 0.41-0.60 at 125 bp.")
    if s_ind:
        print("  ✅ `independent` is present — G5 excludes on the FLAG, not on the name.")
    print("  G1 PASS")

    # SE rows only — this table is skipped exons, like Table 3
    se = conc[conc[c_typ].astype(str).str.upper().str.replace(" ", "") == "SE"]
    if se.empty:
        refuse(f"no rows with event type 'SE' in {os.path.basename(cpath)}.\n"
               f"  values present: {sorted(conc[c_typ].astype(str).unique())}")

    # ------------------------------------------------------------------ G5
    print("\n  G5  only INDEPENDENT cohorts belong in this table")
    if s_ind:
        flag = slop.set_index(slop[s_coh].astype(str).str.lower())[s_ind]
        def indep(name):
            v = flag.get(str(name).lower())
            if v is None:
                return None
            return str(v).strip().lower() in ("true", "yes", "1", "y", "t")
        nonind = [str(c) for c in slop[s_coh] if indep(c) is False]
        print(f"      flagged NOT independent in {os.path.basename(spath)}: "
              f"{nonind or 'none'}")
        if not any("lbnl" in x.lower() for x in nonind):
            refuse("G5 — LBNL is not flagged as non-independent in the slope file.\n\n"
                   "  24 of its 33 samples are the same individuals as discovery-cohort\n"
                   "  samples, verified from genotypes. If the file says it is "
                   "independent, the\n  file is wrong and every pooled number built "
                   "from it is suspect. Fix the flag\n  before building a replication "
                   "table.")
        for coh in list(RECORD):
            if indep(coh) is False:
                refuse(f"G5 — {coh} is flagged NOT independent but is in this script's "
                       f"record as a\n  replication cohort. Resolve before building.")
    else:
        nonind = [str(x) for x in conc[c_coh].unique() if "lbnl" in str(x).lower()]
        print(f"      no `independent` column; excluding by name: {nonind or 'none'}")
    print(f"      (LBNL is a library-prep replicate of discovery. Its {LBNL_SLOPE} / {LBNL_PCT}% was")
    print("       once read as independent corroboration — it is a technical replicate.)")
    print("  G5 PASS")

    # ------------------------------------------------------------------ G2
    print("\n  G2  per-cohort concordance against the record")
    rows, bad = [], []
    for coh, (contrast, acc, n, k_rec, n_rec, sl_rec) in RECORD.items():
        m = se[se[c_coh].astype(str).str.lower() == coh.lower()]
        if m.empty:
            bad.append(f"{coh}: no SE row found in the input")
            continue
        if len(m) > 1:
            bad.append(f"{coh}: {len(m)} SE rows, expected 1")
            continue
        k = int(pd.to_numeric(m[c_k].iloc[0]))
        nn = int(pd.to_numeric(m[c_n].iloc[0]))
        ok = (k == k_rec) and (nn == n_rec)
        if ok and c_pct is not None:
            # the file states the percentage too; confirm it agrees with k/n
            pf = float(pd.to_numeric(m[c_pct].iloc[0]))
            pf = pf * 100 if pf <= 1.0 else pf
            if abs(pf - 100 * k / nn) > 0.06:
                bad.append(f"{coh}: file pct {pf:.2f} but {k}/{nn} = "
                           f"{100 * k / nn:.2f} — the file is inconsistent with itself")
        tag = "OK" if ok else "MISMATCH"
        src = SRC_LABEL.get(coh, "recorded")
        print(f"    {coh:<9} file {k:>4}/{nn:<5}  record {k_rec:>4}/{n_rec:<5}"
              f" ({src:<8}) {tag}")
        if not ok:
            bad.append(f"{coh}: file {k}/{nn}, record {k_rec}/{n_rec}")
        rows.append((coh, contrast, acc, n, k, nn, sl_rec))
    if bad:
        refuse("G2 — the file disagrees with the record:\n" +
               "".join(f"    {b}\n" for b in bad) +
               "\n  For Brigham, Stanford and Fuwai the 50 bp recorded NUMERATOR was derived "
               "from a\n  rounded percentage, so a mismatch here may mean the record is "
               "wrong rather\n  than the file. Resolve which before writing a table "
               "either way.")
    print("  G2 PASS")

    # ------------------------------------------------------------------ G3
    K = sum(r[4] for r in rows)
    N = sum(r[5] for r in rows)
    pv = one_sided(K, N)
    print(f"\n  G3  pooled {K}/{N} = {100 * K / N:.1f}%   record {POOLED_K}/{POOLED_N}"
          f" = {100 * POOLED_K / POOLED_N:.1f}%")
    print(f"      one-sided exact binomial P = {pv:.1e}   record {POOLED_P_TEXT}")
    if (K, N) != (POOLED_K, POOLED_N):
        refuse(f"G3 — pooled {K}/{N} against the recorded {POOLED_K}/{POOLED_N}.")
    if not (POOLED_P_1SF / 2 < pv < POOLED_P_1SF * 2):
        refuse(f"G3 — the recomputed P ({pv:.2e}) is not within a factor of two of the "
               f"recorded {POOLED_P_TEXT}.\n  The concordance matches but the test does not, "
               f"which means the recorded P came\n  from a different test — check "
               f"sidedness before quoting it.")
    print("  G3 PASS")

    # ------------------------------------------------------------------ G4
    print("\n  G4  Deming slopes against the record")
    sl_bad = []
    got = {}
    for coh, _, _, _, _, _, sl_rec in rows:
        m = slop[slop[s_coh].astype(str).str.lower() == coh.lower()]
        if m.empty:
            sl_bad.append(f"{coh}: no row in {os.path.basename(spath)}")
            continue
        v = float(pd.to_numeric(m[s_sl]).iloc[0])
        got[coh] = v
        ok = abs(v - sl_rec) <= SLOPE_TOL
        print(f"    {coh:<9} file {v:.3f}   record {sl_rec:.3f}   "
              f"{'OK' if ok else 'MISMATCH'}")
        if not ok:
            sl_bad.append(f"{coh}: file {v:.4f}, record {sl_rec:.3f}")
    if sl_bad:
        refuse("G4 — slopes disagree with the record:\n" +
               "".join(f"    {b}\n" for b in sl_bad) +
               "\n  ⚠ These are the MAGNet-significant-event slopes, not the genome-wide "
               "ones\n  (Colorado 0.571, Brigham 0.757, Stanford 0.578, Fuwai 0.359). If "
               "the file\n  carries those instead, it is the wrong column or the wrong "
               "file.")
    # ⚠ NEW GATE. The slope file carries its own event count. If it differs from the
    # concordance file's, the slopes were computed on a DIFFERENT event set and the two
    # halves of each row do not describe the same events — which no value check catches.
    if s_n:
        print("\n  G4b the slopes sit on the SAME events as the concordance")
        mism = []
        for coh, _, _, _, _, nn, _ in rows:
            m = slop[slop[s_coh].astype(str).str.lower() == coh.lower()]
            if m.empty:
                continue
            v = int(pd.to_numeric(m[s_n]).iloc[0])
            print(f"    {coh:<9} slope file n={v:<6} concordance n={nn:<6}"
                  f" {'OK' if v == nn else 'MISMATCH'}")
            if v != nn:
                mism.append(f"{coh}: slope n={v}, concordance n={nn}")
        if mism:
            refuse("G4b — the slope and concordance files disagree on the event count:\n" +
                   "".join(f"    {x}\n" for x in mism) +
                   "\n  Each row would then pair a concordance and a slope computed on "
                   "different\n  event sets. Both numbers could be individually correct "
                   "and the row still\n  wrong — this is the project's 'two right "
                   "numbers, one comparison' error.")
        print("  G4b PASS")

    mean_sl = sum(got.values()) / len(got)
    print(f"    mean {mean_sl:.3f}   record {MEAN_SLOPE}")
    if abs(mean_sl - MEAN_SLOPE) > SLOPE_TOL:
        refuse(f"G4 — mean slope {mean_sl:.4f} against the recorded {MEAN_SLOPE}.")
    print("  G4 PASS")

    # ------------------------------------------------------------------ build
    recs = []
    for coh, contrast, acc, n, k, nn, _ in sorted(rows, key=lambda r: -r[4] / r[5]):
        lo, hi = wilson(k, nn)
        rec_r = ""
        if s_r:
            mm = slop[slop[s_coh].astype(str).str.lower() == coh.lower()]
            if not mm.empty:
                rec_r = f"{float(pd.to_numeric(mm[s_r]).iloc[0]):.2f}"
        recs.append({
            "Cohort": coh,
            "Accession": acc,
            "Contrast": contrast,
            "n (disease/control)": n,
            "SE events measurable": f"{nn:,}",
            "Same direction": f"{k:,}",
            "Concordance (%)": f"{100 * k / nn:.1f}",
            "95% CI (%)": f"{lo:.1f}–{hi:.1f}",
            "P (one-sided)": f"{one_sided(k, nn):.1e}".replace("e-0", "e−").replace("e-", "e−"),
            "Transfer slope": f"{got[coh]:.3f}",
            "Pearson r": rec_r,
        })
    lo, hi = wilson(K, N)
    recs.append({
        "Cohort": "Pooled",
        "Accession": "—",
        "Contrast": "—",
        "n (disease/control)": "—",
        "SE events measurable": f"{N:,}",
        "Same direction": f"{K:,}",
        "Concordance (%)": f"{100 * K / N:.1f}",
        "95% CI (%)": f"{lo:.1f}–{hi:.1f}",
        "P (one-sided)": f"{pv:.1e}".replace("e-0", "e−").replace("e-", "e−"),
        "Transfer slope": f"mean {mean_sl:.3f}",
        "Pearson r": "—",
    })
    tbl = pd.DataFrame(recs)
    if not s_r or (tbl["Pearson r"].astype(str).str.strip() == "").all():
        tbl = tbl.drop(columns=["Pearson r"])

    # ------------------------------------------------------------------ G6
    print("\n  G6  no per-cohort significant-event count column")
    offend = [c for c in tbl.columns
              if any(w in c.lower() for w in FORBIDDEN) and "measurable" not in c.lower()]
    if offend:
        refuse("G6 — these columns look like per-cohort significance counts: " +
               ", ".join(offend) +
               "\n  This paper forbids comparing counts across cohorts on its own "
               "evidence.\n  A counts column here would contradict the Methods.")
    print("  G6 PASS")

    hr()
    print(tbl.to_string(index=False))
    hr()

    if a.list_only or not a.apply:
        print("  --list / dry run. Nothing written.")
        print(f"  Would write: {OUT_TSV}, {OUT_MD}")
        hr()
        return

    tbl.to_csv(os.path.join(a.dir, OUT_TSV), sep="\t", index=False)
    with open(os.path.join(a.dir, OUT_MD), "w") as fh:
        fh.write(tbl.to_markdown(index=False))
    print(f"[write] {OUT_TSV}")
    print(f"[write] {OUT_MD}")

    # ------------------------------------------------------------------ caption
    # 170: pointer corrected to Supplementary Figure 2 (Mapping A), slope
    # point-estimate caveat added, pooled sentence dropped as duplicated.
    hr("=")
    print("CAPTION — every number read from a file or recomputed from it")
    hr("=")
    print(f"""
**Table 4. The discovery signature replicates by direction in four independent
cohorts.** Skipped exons significant in the discovery contrast
(FDR < 0.05, |ΔPSI| ≥ 0.10) and measurable in each replication cohort, exon-level, JCEC
quantification. Concordance is the proportion changing in the same direction as discovery,
tested against 50% with a one-sided exact binomial test; intervals are Wilson 95%
intervals. The transfer slope is a Deming (errors-in-variables) regression of the
replication cohort's ΔPSI on the discovery cohort's, with λ = mean(se_y²)/mean(se_x²) and
the discovery cohort on the x axis. Slopes are point estimates without confidence
intervals. Group sizes are post-exclusion, as in Table 1.

**Neither effect-size magnitude nor event counts are compared across cohorts**
(Methods): the cohorts differ in read length, sequencing depth, covariate model and
power. Every slope is attenuated by selection on discovery significance and is a lower
bound; the matched-replicate panel separates the technical component
(Supplementary Figure 2). Concordance and slope do not rank together.

ᵃ LBNL is excluded from this table and from every pooled analysis: 24 of its 33 samples
  are the same individuals as discovery-cohort samples, verified from genotypes, so it is
  a library-preparation replicate rather than an independent cohort.
ᵇ Mutually exclusive exons are excluded for Colorado and Fuwai on reverse-replication
  grounds, on different evidence in each case (Supplementary Table 8).
""")
    hr("=")
    print("""BEFORE THIS GOES IN THE PAPER — three things this script cannot check:

  1. The number of individually significant cohort x event-type combinations is
     deliberately NOT stated here. It moves with read length (15 of 18 at 50 bp,
     14 of 18 at 125 bp) and so do the identity and the event class of the cells
     that fail. Run `85_recount_combinations.py`, quote its scenario C, and say
     what it counts.
     Do NOT write that the failing cells are the lowest concordance rates. At
     125 bp two SIGNIFICANT cells fall below a non-significant one — Fuwai RI
     (51/80, 63.8%) and Fuwai A5SS (36/55, 65.5%) against Brigham RI (21/32,
     65.6%) — because significance there follows the number of testable events,
     not the rate. At 50 bp the two coincided, which is why the old claim held.
     Counted from `cohort_type_concordance_dedup[_rl125].tsv` by `85_recount_combinations.py`.
  2. Per-cohort Pearson r is now IN the table, read from the slope file's `pearson_r`
     column. The record held it only as the range 0.36-0.59 at 50 bp / 0.41-0.60 at 125 bp; check the per-cohort values
     fall inside that range before quoting either.
  3. Per-cohort slope confidence intervals are not in the record at all. The table gives
     point estimates only; say so, or compute them.

⚠ And the one to be careful about: the three HCM cohorts' OWN significant counts
  (Brigham 4, Fuwai 169, Stanford 8) are deliberately absent. Table 3's G8 exists to
  keep that a recorded decision. Do not add them here either.""")
    hr("=")


if __name__ == "__main__":
    main()
