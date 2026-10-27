#!/usr/bin/env python3
"""
85_recount_combinations.py  —  recount "sixteen of nineteen cohort x event-type
                               combinations", and the all-five-types pooled
                               concordance, now that Fuwai MXE is excluded.

Source: `cohort_type_concordance_dedup.tsv` (cohort, event_type, same_direction,
n, pct) — 24 combinations — cross-checked against `direction_concordance_dedup.tsv`
(one SE row per cohort).

TWO NUMBERS ARE IN QUESTION, AND A THIRD PROBLEM TURNED UP
----------------------------------------------------------
1. **"16 of 19"** (PAPER_ONE_PAGE claim 4). The record says the denominator "was 23
   with LBNL's five combinations included" and "drops to 18 now that Fuwai MXE is
   excluded". Those two statements cannot both be right: 24 rows minus LBNL's 5 is
   19, not 23, and 19 minus Fuwai MXE is 18 only if Fuwai MXE is one of the 19.
   This script counts, rather than patching the arithmetic.

2. **2,227/2,991 = 74.5%**, the all-five-types pooled concordance, which the record
   flags as including Fuwai MXE and therefore moving.

3. **THE P-VALUES ARE ONE-SIDED.** `direction_concordance_dedup.tsv`'s stored p is
   exactly half the two-sided binomial in all five cohort rows, and the record's
   three borderline figures (0.07, 0.08, 0.14) are one-sided too. That is defensible
   — the alternative of interest is better-than-chance agreement — but it is not
   stated in the Methods, AND the "individually significant at p<0.05" count depends
   on it. A combination at one-sided p = 0.03 is two-sided p = 0.06 and would stop
   being significant. **So the recount is reported under BOTH conventions**, and any
   combination whose verdict flips is named. If the count differs, the Methods must
   declare the sidedness before the sentence can be written.

GATES
  G1  structure: 24 rows, 5 cohorts, which cohort x type combinations are absent.
  G2  the file reproduces the three recorded borderline combinations exactly
      (Brigham A5SS 19/29, Fuwai A3SS 42/71, Fuwai A5SS 39/68). If not, this is
      not the table the record describes -> REFUSE.
  G3  cross-file: the SE row of every cohort must equal `direction_concordance_dedup.tsv`.
      Two files, two scripts, one answer -> REFUSE on mismatch.
  G4  sidedness: which convention reproduces the recorded 0.07 / 0.08 / 0.14.
  G5  the recount, four scenarios x two sidedness conventions.
  G6  the all-five-types pooled figure: reproduce 2,227/2,991, then recompute with
      Fuwai MXE (and Colorado MXE, already absent) out.

USAGE
    python 85_recount_combinations.py --list      # print the table and the structure
    python 85_recount_combinations.py             # the gates and the recount
"""

import argparse
import os
import sys

import pandas as pd
from scipy.stats import binomtest

TYPES = ["SE", "MXE", "A3SS", "A5SS", "RI"]
COHORTS = ["Colorado", "LBNL", "Brigham", "Fuwai", "Stanford"]
INDEPENDENT = ["Colorado", "Brigham", "Fuwai", "Stanford"]   # LBNL shares individuals
EXCLUDED_MXE = ["Colorado", "Fuwai"]                         # reverse-replication grounds

# the three combinations the record names as NOT individually significant
REC_BORDERLINE = {
    ("Brigham", "A5SS"): (19, 29, 0.07),
    ("Fuwai", "A3SS"):   (42, 71, 0.08),
    ("Fuwai", "A5SS"):   (39, 68, 0.14),
}
# direction_concordance_dedup.tsv, SE only
REC_SE = {"Colorado": (375, 471), "LBNL": (151, 199), "Brigham": (359, 442),
          "Fuwai": (494, 714), "Stanford": (408, 536)}
REC_POOLED_ALLTYPES = (2227, 2991)      # record's all-five-types pooled figure
REC_CLAIM = (16, 19)                    # "16 of 19"

# NEW: 1,636/2,163 occurred as an inline tuple (line ~328) and as prose
# (line ~331), neither of them the named REC_SE — Rule 8. One constant now.
REC_SE_POOLED = (1636, 2163)

# --------------------------------------------------------------------------
# READ LENGTH and the 125 bp record, placed AFTER every literal it rewrites.
#
# When TBL_SUFFIX is set all five records are REQUIRED. REC_BORDERLINE and
# REC_CLAIM may be the literal `none`: at a new read length they are unknown
# (G5 is what discovers them), and an EMPTY REC_BORDERLINE is not a neutral
# default — see the G2 and G4 guards further down, which this patch also adds.
SUF = os.environ.get("TBL_SUFFIX", "")
_REQ = ("REC_SE", "REC_SE_POOLED", "REC_ALLTYPES", "REC_BORDERLINE",
        "REC_CLAIM")
_env = {_k: os.environ.get(_k, "").strip() for _k in _REQ}
if SUF:
    _miss = [_k for _k in _REQ if not _env[_k]]
    if _miss:
        raise SystemExit(
            "REFUSING — TBL_SUFFIX=%s but these are unset: %s\n\n"
            "  REC_BORDERLINE and REC_CLAIM accept the literal `none` when the\n"
            "  value is not yet known at this read length. Say so explicitly —\n"
            "  an empty record is NOT a neutral default here.\n\n"
            "    %s\n" % (SUF, ", ".join(_miss), 'TBL_SUFFIX=_rl125 \\\n    REC_SE=Colorado=311/385,LBNL=122/156,Brigham=298/364,Fuwai=422/596,Stanford=345/455 \\\n    REC_SE_POOLED=1376,1800 \\\n    REC_ALLTYPES=1883,2505 \\\n    REC_BORDERLINE=none \\\n    REC_CLAIM=none \\\n    python3 85_recount_combinations.py'))

if _env["REC_SE"]:
    _d = dict(_p.split("=", 1) for _p in _env["REC_SE"].split(","))
    assert set(_d) == set(REC_SE), (
        "REC_SE must name exactly %s, got %s" % (sorted(REC_SE), sorted(_d)))
    REC_SE = {_c: tuple(int(_x) for _x in _v.split("/"))
              for _c, _v in _d.items()}

if _env["REC_SE_POOLED"]:
    REC_SE_POOLED = tuple(int(_x) for _x in _env["REC_SE_POOLED"].split(","))
if _env["REC_ALLTYPES"]:
    REC_POOLED_ALLTYPES = tuple(int(_x)
                                for _x in _env["REC_ALLTYPES"].split(","))

# REC_CLAIM: `none` -> (0, 0). The scenario() call already WARNS rather than
# refuses on a mismatch, and (0, 0) reads as "no published claim to compare".
if _env["REC_CLAIM"]:
    REC_CLAIM = ((0, 0) if _env["REC_CLAIM"].lower() == "none"
                 else tuple(int(_x) for _x in _env["REC_CLAIM"].split(",")))

# REC_BORDERLINE: `none` -> {}. Format otherwise
#   cohort:TYPE=k/n@p,cohort:TYPE=k/n@p
if _env["REC_BORDERLINE"]:
    if _env["REC_BORDERLINE"].lower() == "none":
        REC_BORDERLINE = {}
    else:
        _nb = {}
        for _item in _env["REC_BORDERLINE"].split(","):
            _key, _val = _item.split("=", 1)
            _c, _ty = _key.split(":", 1)
            # `@p` is OPTIONAL. With it, G2 and G4 both run. Without it, only
            # G2 runs — see the _RB_P note below. A recorded p at a NEW read
            # length can only come from this script's own one-sided column,
            # which would make G4 infer the convention from its own output.
            if "@" in _val:
                _kn, _p = _val.split("@", 1)
                _pv = float(_p)
            else:
                _kn, _pv = _val, None
            _a, _b = _kn.split("/")
            _nb[(_c, _ty)] = (int(_a), int(_b), _pv)
        REC_BORDERLINE = _nb

# G4 may use ONLY the entries that carry a PUBLISHED p-value. G2 is a fixity
# check on k/n and uses all of them; G4 is a provenance inference and cannot
# be fed a number this script produced. An entry written without `@p` is a
# deliberate statement that no published p exists for it.
_RB_P = {_k: _v for _k, _v in REC_BORDERLINE.items() if _v[2] is not None}
# --------------------------------------------------------------------------
ALPHA = 0.05

COMB = f"cohort_type_concordance_dedup{SUF}.tsv"
SE1 = f"direction_concordance_dedup{SUF}.tsv"


def refuse(msg):
    print("\n" + "=" * 74)
    print("REFUSING — " + msg)
    print("=" * 74)
    sys.exit(2)


def p_one(s, n):
    return binomtest(s, n, 0.5, alternative="greater").pvalue


def p_two(s, n):
    return binomtest(s, n, 0.5, alternative="two-sided").pvalue


def load(d):
    p = os.path.join(d, COMB)
    if not os.path.exists(p):
        refuse(f"{COMB} not found in {d}")
    t = pd.read_csv(p, sep="\t")
    need = {"cohort", "event_type", "same_direction", "n"}
    if not need.issubset(t.columns):
        refuse(f"{COMB} is missing columns.\n  need: {sorted(need)}\n"
               f"  has : {list(t.columns)}")
    t["cohort"] = t.cohort.astype(str)
    t["event_type"] = t.event_type.astype(str).str.upper()
    t["same_direction"] = t.same_direction.astype(int)
    t["n"] = t.n.astype(int)
    return t


def table(t, title, show_p=True):
    print("\n" + title)
    print(f"  {'cohort':10s} {'type':5s} {'same':>6s} {'n':>6s} {'pct':>7s} "
          f"{'p 1-sided':>11s} {'p 2-sided':>11s}  verdict(1s/2s)")
    for _, r in t.sort_values(["cohort", "event_type"]).iterrows():
        p1, p2 = p_one(r.same_direction, r.n), p_two(r.same_direction, r.n)
        v = ("sig " if p1 < ALPHA else "ns  ") + ("/ sig" if p2 < ALPHA else "/ ns ")
        flag = "   <-- FLIPS" if (p1 < ALPHA) != (p2 < ALPHA) else ""
        print(f"  {r.cohort:10s} {r.event_type:5s} {r.same_direction:6d} {r.n:6d} "
              f"{100*r.same_direction/r.n:6.2f}% {p1:11.4g} {p2:11.4g}  {v}{flag}")


def scenario(t, label, recn=None, recs=None):
    n_comb = len(t)
    s1 = int(sum(p_one(r.same_direction, r.n) < ALPHA for _, r in t.iterrows()))
    s2 = int(sum(p_two(r.same_direction, r.n) < ALPHA for _, r in t.iterrows()))
    ns1 = [f"{r.cohort} {r.event_type}" for _, r in t.iterrows()
           if p_one(r.same_direction, r.n) >= ALPHA]
    ns2 = [f"{r.cohort} {r.event_type}" for _, r in t.iterrows()
           if p_two(r.same_direction, r.n) >= ALPHA]
    tag = ""
    if recn is not None:
        tag = ("   [record " + f"{recs} of {recn}" +
               ("  MATCH" if (n_comb == recn and s1 == recs) else "  MISMATCH") + "]")
    print(f"\n  {label}")
    print(f"    combinations           : {n_comb}{tag}")
    print(f"    significant, 1-sided   : {s1} of {n_comb}")
    print(f"    significant, 2-sided   : {s2} of {n_comb}")
    print(f"    not significant 1-sided: {', '.join(ns1) if ns1 else '(none)'}")
    if ns2 != ns1:
        print(f"    not significant 2-sided: {', '.join(ns2)}")
    return dict(n=n_comb, s1=s1, s2=s2, ns1=ns1, ns2=ns2)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default=".")
    ap.add_argument("--list", action="store_true")
    args = ap.parse_args()

    t = load(args.dir)

    if args.list:
        print("=" * 74)
        print(f"INVENTORY — {COMB}")
        print("=" * 74)
        table(t, "all rows:")
        print("\nstructure (cohort x event type; '.' = absent):")
        print(f"  {'':10s} " + " ".join(f"{ty:>5s}" for ty in TYPES))
        for c in COHORTS:
            cells = []
            for ty in TYPES:
                m = t[(t.cohort == c) & (t.event_type == ty)]
                cells.append(f"{int(m.n.iloc[0]):5d}" if len(m) else "    .")
            print(f"  {c:10s} " + " ".join(cells))
        extra = sorted(set(t.cohort) - set(COHORTS))
        if extra:
            print(f"  unexpected cohorts: {extra}")
        print(f"\n  rows: {len(t)}   cohorts: {t.cohort.nunique()}   "
              f"types: {sorted(t.event_type.unique())}")
        return

    print("=" * 74)
    print('85 — recount "16 of 19", and the all-five-types pooled concordance')
    print("=" * 74)

    # ------------------------------------------------------------------ G1
    print("\n" + "-" * 74)
    print("G1  structure")
    print("-" * 74)
    print(f"  {'':10s} " + " ".join(f"{ty:>7s}" for ty in TYPES))
    absent = []
    for c in COHORTS:
        cells = []
        for ty in TYPES:
            m = t[(t.cohort == c) & (t.event_type == ty)]
            if len(m):
                cells.append(f"{int(m.same_direction.iloc[0])}/{int(m.n.iloc[0])}".rjust(7))
            else:
                cells.append("      .")
                absent.append((c, ty))
        print(f"  {c:10s} " + " ".join(cells))
    print(f"\n  rows: {len(t)}")
    print(f"  absent combinations: {', '.join(f'{c} {ty}' for c, ty in absent) or '(none)'}")
    if ("Colorado", "MXE") not in absent:
        print("  ⚠ Colorado MXE is PRESENT in this table but is excluded project-wide.")
    if ("Fuwai", "MXE") in absent:
        print("  note: Fuwai MXE is already absent — this table postdates the exclusion.")

    # ------------------------------------------------------------------ G2
    print("\n" + "-" * 74)
    print("G2  does the file reproduce the three recorded borderline combinations?")
    print("-" * 74)
    ok2 = True
    for (c, ty), (rs, rn, rp) in REC_BORDERLINE.items():
        m = t[(t.cohort == c) & (t.event_type == ty)]
        if not len(m):
            print(f"  {c} {ty}: ABSENT from the table, record says {rs}/{rn}")
            ok2 = False
            continue
        s, n = int(m.same_direction.iloc[0]), int(m.n.iloc[0])
        hit = (s, n) == (rs, rn)
        ok2 &= hit
        print(f"  {c:9s} {ty:5s} file {s}/{n}   record {rs}/{rn}   "
              f"{'match' if hit else 'MISMATCH'}")
    if not ok2:
        refuse("G2 — this is not the table the record describes. Do not recount from it.\n"
               "  Paste `--list` output.")
    if not REC_BORDERLINE:
        # An empty dict makes the loop above a no-op and leaves ok2 True, so
        # this gate used to print "G2 PASS" having checked nothing at all.
        print("  G2 NOT RUN — REC_BORDERLINE is empty, so nothing was")
        print("  reproduced. This is NOT a pass. At a new read length the")
        print("  borderline combinations are unknown until G5 below names")
        print("  them; read them off G5's 'Not individually significant'")
        print("  line and re-run with")
        print("    REC_BORDERLINE=Cohort:TYPE=k/n@p,Cohort:TYPE=k/n@p")
    else:
        print("  G2 PASS")

    # ------------------------------------------------------------------ G3
    print("\n" + "-" * 74)
    print(f"G3  cross-file: SE rows vs {SE1}")
    print("-" * 74)
    p1f = os.path.join(args.dir, SE1)
    if not os.path.exists(p1f):
        print(f"  SKIPPED — {SE1} not found. G5 still stands, but without this there is\n"
              f"  no second file agreeing with the first.")
    else:
        s1 = pd.read_csv(p1f, sep=r"\s+")
        ok3 = True
        for _, r in s1.iterrows():
            c = str(r.cohort)
            m = t[(t.cohort == c) & (t.event_type == "SE")]
            if not len(m):
                print(f"  {c:10s} SE absent from {COMB}")
                ok3 = False
                continue
            a = (int(m.same_direction.iloc[0]), int(m.n.iloc[0]))
            b = (int(r["same"]), int(r["n"]))
            rec = REC_SE.get(c)
            _rec_txt = f"{rec[0]}/{rec[1]}" if rec else "(none)"
            hit = a == b and (rec is None or a == rec)
            ok3 &= hit
            print(f"  {c:10s} {COMB.split('_')[0]+'...':>6s} {a[0]}/{a[1]}   "
                  f"{SE1.split('_')[0]+'...':>9s} {b[0]}/{b[1]}   "
                  f"record {_rec_txt}   {'match' if hit else 'MISMATCH'}")
        if not ok3:
            refuse("G3 — the two concordance files disagree on SE, or disagree with the\n"
                   "  record. One of them is stale. Resolve that before recounting.")
        print("  G3 PASS — two files, two scripts, identical SE rows.")

    # ------------------------------------------------------------------ G4
    print("\n" + "-" * 74)
    _rb_ps = " / ".join(f"{_rp:.2f}"
                       for (_, _, _rp) in _RB_P.values())
    print("G4  sidedness: which convention produced the recorded "
          + (_rb_ps + "?" if _rb_ps else "p-values? (none recorded)"))
    print("-" * 74)
    print(f"  {'combination':16s} {'same/n':>9s} {'1-sided':>10s} {'2-sided':>10s} {'record':>8s}")
    for (c, ty), (rs, rn, rp) in _RB_P.items():
        print(f"  {c+' '+ty:16s} {f'{rs}/{rn}':>9s} {p_one(rs,rn):10.4f} "
              f"{p_two(rs,rn):10.4f} {rp:8.2f}")
    d1 = sum(abs(p_one(rs, rn) - rp) for (rs, rn, rp) in _RB_P.values())
    d2 = sum(abs(p_two(rs, rn) - rp) for (rs, rn, rp) in _RB_P.values())
    if not _RB_P:
        # d1 == d2 == 0.0 for an empty record, and `d1 < d2` is then False, so
        # this gate used to ANNOUNCE "TWO-SIDED" from no data whatsoever.
        # Sidedness is settled independently: a one-sided exact binomial at
        # p = 0.5 reproduces every p in direction_concordance_dedup*.tsv
        # exactly — all five cohorts, both read lengths.
        sided = "ONE-SIDED"
        print("\n  G4 NOT RUN — no PUBLISHED p-value is recorded for any")
        print("  borderline combination, so there is no convention to infer.")
        print("  (A p-value taken from this script's own one-sided column")
        print("   would make this gate confirm itself. Entries written as")
        print("   k/n with no @p say so deliberately.)")
        print("  Taking ONE-SIDED as established independently")
        print("  elsewhere: a one-sided")
        print("  exact binomial at p=0.5 reproduces every p in")
        print(f"  {SE1} exactly, at both read lengths.")
    else:
        sided = "ONE-SIDED" if d1 < d2 else "TWO-SIDED"
        print(f"\n  total absolute deviation: one-sided {d1:.4f}, "
              f"two-sided {d2:.4f}")
    print(f"  G4 — the record's p-values are {sided}.")
    if sided == "ONE-SIDED":
        print("  ⚠ The Methods must say so. A one-sided binomial against 0.5 is\n"
              "    defensible here (the alternative of interest is better-than-chance\n"
              "    agreement), but it is not currently stated, and the counts below\n"
              "    show whether it changes any claim.")

    # ------------------------------------------------------------------ G5
    print("\n" + "-" * 74)
    print("G5  THE RECOUNT")
    print("-" * 74)
    table(t, "every combination in the file:")

    as_pub = scenario(t, "A. as the file stands (all cohorts, all present types)")
    indep = t[t.cohort.isin(INDEPENDENT)]
    # `or None` so an unrecorded claim (0, 0) prints as absent rather
    # than as a MISMATCH against zero — the other three scenario()
    # calls already pass None by omission. This was the THIRD site of
    # REC_CLAIM; 134 fixed the verdict line and the warning only.
    r_ind = scenario(indep, "B. independent cohorts only (LBNL removed)",
                     recn=REC_CLAIM[1] or None,
                     recs=REC_CLAIM[0] or None)
    final = indep[~((indep.cohort.isin(EXCLUDED_MXE)) & (indep.event_type == "MXE"))]
    r_fin = scenario(final, "C. independent cohorts, Colorado AND Fuwai MXE excluded"
                            "  <-- THE NUMBER TO QUOTE")
    lbnl = t[t.cohort == "LBNL"]
    if len(lbnl):
        scenario(lbnl, "D. LBNL alone (for the record — not quoted)")

    _claim_txt = ("not recorded at this read length"
                  if tuple(REC_CLAIM) == (0, 0)
                  else f'"{REC_CLAIM[0]} of {REC_CLAIM[1]}"')
    print(f"""
  VERDICT on the claim sentence
  -----------------------------
  As published : {_claim_txt}
  Recount (C)  : **{r_fin['s1']} of {r_fin['n']}** one-sided"""
          + (f", {r_fin['s2']} of {r_fin['n']} two-sided"
             if r_fin["s2"] != r_fin["s1"] else " (two-sided identical)"))
    if r_fin["s2"] != r_fin["s1"]:
        flip = sorted(set(r_fin["ns2"]) - set(r_fin["ns1"]))
        print(f"  ⚠ SIDEDNESS CHANGES THE COUNT. Combinations significant one-sided but\n"
              f"    not two-sided: {', '.join(flip)}.\n"
              f"    The Methods MUST declare the test before this sentence is written.")
    else:
        print("  Sidedness does not change the count, so the disclosure is a Methods\n"
              "  clause only and the sentence is safe either way.")
    print(f"  Not individually significant: {', '.join(r_fin['ns1'])}")
    if tuple(REC_CLAIM) != (0, 0) and (
            r_ind["n"] != REC_CLAIM[1] or r_ind["s1"] != REC_CLAIM[0]):
        print(f"\n  ⚠ Scenario B does NOT reproduce \"{REC_CLAIM[0]} of {REC_CLAIM[1]}\" "
              f"— it gives {r_ind['s1']} of {r_ind['n']}.\n"
              f"    So the published figure was not 'independent cohorts, all present\n"
              f"    types'. Whatever it was, quote the recount above and say what it\n"
              f"    counts; do not reuse the old sentence's denominator.")

    # ------------------------------------------------------------------ G6
    print("\n" + "-" * 74)
    print("G6  the all-five-types POOLED concordance")
    print("-" * 74)

    def pooled(sub, label, rec=None):
        s, n = int(sub.same_direction.sum()), int(sub.n.sum())
        tag = ""
        if rec:
            tag = ("   [record " + f"{rec[0]}/{rec[1]}" +
                   ("  MATCH" if (s, n) == rec else "  MISMATCH") + "]")
        print(f"  {label:52s} {s:>6d}/{n:<6d} = {100*s/n:6.3f}%   "
              f"p1={p_one(s,n):.3g}{tag}")
        return s, n

    pooled(indep, "independent cohorts, all present types (as published)",
           rec=REC_POOLED_ALLTYPES)
    new = pooled(final, "independent cohorts, Colorado + Fuwai MXE excluded")
    pooled(indep[indep.event_type == "SE"], "independent cohorts, SE only (the headline)",
           rec=REC_SE_POOLED)
    print(f"\n  The all-five-types pooled figure becomes **{new[0]}/{new[1]} = "
          f"{100*new[0]/new[1]:.1f}%**.")
    print(f"  The SE headline ({REC_SE_POOLED[0]:,}/{REC_SE_POOLED[1]:,} = "
          f"{100 * REC_SE_POOLED[0] / REC_SE_POOLED[1]:.1f}%) is unaffected — it never included MXE.")
    print("  Never substitute one for the other; say which the sentence means.")

    print("\n" + "=" * 74)
    print("Every number above is counted from the file. Nothing is inferred by")
    print("subtracting one from a published figure.")
    print("=" * 74)


if __name__ == "__main__":
    main()
