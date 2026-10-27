#!/usr/bin/env python3
"""
84_pool_b.py  —  close out POOL B = 526, the last number in the Results with
                 no file behind it.

WHAT IS ACTUALLY IN QUESTION
----------------------------
`RESULT_phase1_combined_evidence.md` records Pool B (Brigham + Fuwai +
Stanford, `any` effect rule) over five classes:

    type   events    FDR<0.05   +effect   +unanimous   retention
    SE     11,464       853       413        363          88%
    MXE     6,048       249        78         70          90%   <- now excluded
    A3SS    2,753        60        31         29          94%
    A5SS    1,757        45        24         21          88%
    RI      1,727       243       125        113          90%
    total  23,749                            596

`DECISION_fuwai_mxe_and_permutations.md` excludes Fuwai MXE, which removes the
MXE row (Pool B needs unanimity across all three cohorts), and the Results
quote 363 + 29 + 21 + 113 = 526.

That subtraction is only valid IF the Benjamini-Hochberg correction was applied
WITHIN each event type.  If BH was applied ACROSS all five types pooled, then
deleting 6,048 MXE rows changes the ranks, changes every q-value, and changes
all four surviving rows.  526 would then be wrong, not merely untraced.

So the question is not "re-run the pooling" — it is "what was the BH scope",
and that is answerable from the per-event Pool B table that already exists on
disk.  This script:

  G1  parses `combined_b.tsv` and reproduces the recorded five-class ladder
      (11,464 / 853 / 413 / 363 for SE, and the other four rows).  If it
      cannot, the parse is wrong and everything below it is meaningless -> REFUSE.

  G2  identifies the BH scope by recomputing q from the combined p under both
      scopes and seeing which one reproduces the recorded FDR<0.05 counts.
      If neither reproduces them -> REFUSE.  If both do -> says so and falls
      through to the conservative branch (treat as pooled).

  G3  applies the identified scope to the table with MXE dropped and reports
      the true four-class ladder and total, then compares it to 526.

  G4  cross-checks Pool A (`combined_a.tsv`, `discovery` rule) against its own
      record: 27,817 events, 409 SE unanimous, 550 total.  This ties the files
      on disk to the recorded run; if Pool A does not reproduce, the Pool B
      file is from a different vintage and nothing here can be trusted.

Nothing is guessed.  Every column is resolved by name against a candidate list
and the script refuses with the full column list when it cannot find one.

USAGE
    python 84_pool_b.py --list        # inventory only: files, headers, the
                                      # --types metavar from script 35's source,
                                      # and the re-run invocations.  Run this first.

    python 84_pool_b.py              # the gates

    python 84_pool_b.py --pool-b combined_b.tsv --pool-a combined_a.tsv
                                      # explicit paths if the defaults are wrong
"""

import argparse
import os
import re
import sys

import numpy as np
import pandas as pd

# ---------------------------------------------------------------- the record

TYPES5 = ["SE", "MXE", "A3SS", "A5SS", "RI"]
TYPES4 = ["SE", "A3SS", "A5SS", "RI"]

# RESULT_phase1_combined_evidence.md, Pool B, `any` rule, five classes
REC_B = {
    "SE":   dict(events=11464, fdr=853, effect=413, unan=363),
    "MXE":  dict(events=6048,  fdr=249, effect=78,  unan=70),
    "A3SS": dict(events=2753,  fdr=60,  effect=31,  unan=29),
    "A5SS": dict(events=1757,  fdr=45,  effect=24,  unan=21),
    "RI":   dict(events=1727,  fdr=243, effect=125, unan=113),
}
REC_B_EVENTS = 23749
REC_B_UNAN   = 596
QUOTED_526   = 526          # what the Results currently say

# RESULT_phase1_combined_evidence.md, Pool A, `discovery` rule
REC_A = {
    "SE":   dict(events=18111, fdr=4797, effect=491, unan=409, new=34),
    "A3SS": dict(events=4072,  fdr=531,  effect=51,  unan=45,  new=6),
    "A5SS": dict(events=2635,  fdr=459,  effect=45,  unan=37,  new=2),
    "RI":   dict(events=2999,  fdr=712,  effect=74,  unan=59,  new=6),
}
REC_A_EVENTS = 27817
REC_A_UNAN   = 550

FDR_CUT  = 0.05
DPSI_CUT = 0.10

# ------------------------------------------------------------- column naming

TYPE_CANDIDATES = ("type", "event_type", "etype", "class", "event_class", "EventType")

Q_CANDIDATES = ("FDR", "fdr", "FDR_combined", "FDR_stouffer", "qval", "q", "padj",
                "p_adj", "FDR_z", "FDR_meta")

P_CANDIDATES = ("p_stouffer", "P_stouffer", "p_combined", "P_combined", "pval",
                "p_value", "PValue", "p", "P", "pval_z", "p_z", "p_meta")

# per-cohort dPSI columns: anything that looks like dPSI_<cohort>
DPSI_RE = re.compile(r"^(?:dPSI|delta_?psi|dpsi)[._-]?(?P<coh>[A-Za-z0-9]+)$", re.I)
# columns that are NOT per-cohort (pooled / adjusted / reference variants)
DPSI_EXCLUDE = re.compile(r"(pooled|combined|mean|meta|adj|raw|abs|max|min)$", re.I)


def pick(cols, candidates, what):
    """First candidate present, exact then case-insensitive."""
    for c in candidates:
        if c in cols:
            return c
    low = {c.lower(): c for c in cols}
    for c in candidates:
        if c.lower() in low:
            return low[c.lower()]
    refuse(f"cannot find a {what} column.\n"
           f"  tried: {', '.join(candidates)}\n"
           f"  table has: {', '.join(cols)}")


def cohort_dpsi_cols(cols):
    out = {}
    for c in cols:
        m = DPSI_RE.match(c)
        if not m:
            continue
        coh = m.group("coh")
        if DPSI_EXCLUDE.search(coh):
            continue
        out[coh] = c
    return out


def refuse(msg):
    print("\n" + "=" * 72)
    print("REFUSING — " + msg)
    print("=" * 72)
    sys.exit(2)


# ------------------------------------------------------------------- BH

def bh(p):
    """Benjamini-Hochberg q-values. p is a 1-D float array, NaN allowed."""
    p = np.asarray(p, dtype=float)
    q = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    n = int(ok.sum())
    if n == 0:
        return q
    idx = np.where(ok)[0]
    order = idx[np.argsort(p[idx], kind="mergesort")]
    ranked = p[order] * n / np.arange(1, n + 1)
    q[order] = np.minimum.accumulate(ranked[::-1])[::-1].clip(max=1.0)
    return q


# --------------------------------------------------------------- the ladder

def ladder(df, types, c_type, c_q, dcols, q=None):
    """
    events -> FDR<0.05 -> +effect (|dPSI| >= 0.10 in ANY cohort) -> +unanimous.
    `q` overrides the q column (for recomputed scopes).
    Returns {type: dict(events, fdr, effect, unan)}.
    """
    qv = df[c_q].to_numpy(dtype=float) if q is None else np.asarray(q, dtype=float)
    D = df[list(dcols.values())].to_numpy(dtype=float)
    with np.errstate(invalid="ignore"):
        any_effect = (np.abs(D) >= DPSI_CUT).any(axis=1)
        pos = (D > 0).all(axis=1)
        neg = (D < 0).all(axis=1)
    unanimous = pos | neg

    sig = np.isfinite(qv) & (qv < FDR_CUT)
    t = df[c_type].astype(str).str.upper().to_numpy()

    out = {}
    for ty in types:
        m = (t == ty.upper())
        out[ty] = dict(
            events=int(m.sum()),
            fdr=int((m & sig).sum()),
            effect=int((m & sig & any_effect).sum()),
            unan=int((m & sig & any_effect & unanimous).sum()),
        )
    return out


def show(tab, types, rec=None, title=""):
    if title:
        print("\n" + title)
    hdr = f"  {'type':6s} {'events':>8s} {'FDR<.05':>8s} {'+effect':>8s} {'+unanim':>8s} {'reten':>6s}"
    if rec:
        hdr += "   vs record"
    print(hdr)
    tot = dict(events=0, fdr=0, effect=0, unan=0)
    allok = True
    for ty in types:
        r = tab[ty]
        for k in tot:
            tot[k] += r[k]
        ret = f"{100.0 * r['unan'] / r['effect']:.0f}%" if r["effect"] else "  -"
        line = (f"  {ty:6s} {r['events']:8,d} {r['fdr']:8,d} {r['effect']:8,d} "
                f"{r['unan']:8,d} {ret:>6s}")
        if rec:
            e = rec.get(ty)
            if e is None:
                line += "   (not in record)"
            else:
                diffs = [f"{k} {r[k]} vs {e[k]}" for k in ("events", "fdr", "effect", "unan")
                         if r[k] != e[k]]
                if diffs:
                    allok = False
                    line += "   MISMATCH: " + "; ".join(diffs)
                else:
                    line += "   match"
        print(line)
    print(f"  {'total':6s} {tot['events']:8,d} {'':8s} {'':8s} {tot['unan']:8,d}")
    return tot, allok


# ------------------------------------------------------------------- --list

def do_list(args):
    print("=" * 72)
    print("INVENTORY — Pool B / Pool A files and script 35's interface")
    print("=" * 72)

    wanted = ["combined_b.tsv", "combined_a.tsv", "combined_evidence.txt",
              "combined_evidence_disc.txt", "meta_a.tsv", "meta_b.tsv",
              "meta_analysis.txt", "35_combined_evidence.py",
              "master_numbers_flagged.tsv", "direction_concordance_dedup.tsv"]
    print("\nfiles:")
    for f in wanted:
        p = os.path.join(args.dir, f)
        if os.path.exists(p):
            print(f"  [x] {f:34s} {os.path.getsize(p):>12,d} bytes")
        else:
            print(f"  [ ] {f:34s} MISSING")

    for f in ("combined_b.tsv", "combined_a.tsv"):
        p = os.path.join(args.dir, f)
        if not os.path.exists(p):
            continue
        print(f"\n--- {f}: columns")
        df = pd.read_csv(p, sep="\t", nrows=5)
        for c in df.columns:
            print(f"      {c}")
        print(f"    rows (full): {sum(1 for _ in open(p)) - 1:,d}")
        d = cohort_dpsi_cols(df.columns)
        print(f"    per-cohort dPSI columns detected: "
              f"{', '.join(f'{k}->{v}' for k, v in d.items()) or 'NONE'}")

    for f in ("combined_evidence.txt", "combined_evidence_disc.txt"):
        p = os.path.join(args.dir, f)
        if not os.path.exists(p):
            continue
        print(f"\n--- {f}: first 40 lines")
        with open(p) as fh:
            for i, line in enumerate(fh):
                if i >= 40:
                    print("      ...")
                    break
                print("      " + line.rstrip())

    # the --types metavar, read off script 35's own source
    p35 = os.path.join(args.dir, "35_combined_evidence.py")
    if os.path.exists(p35):
        print("\n--- 35_combined_evidence.py: add_argument lines")
        src = open(p35).read()
        for m in re.finditer(r"add_argument\(\s*[\"'](--[a-z0-9-]+)[\"'][^)]*\)", src, re.S):
            print("      " + " ".join(m.group(0).split()))
        print("\n--- 35_combined_evidence.py: BH / multipletests calls "
              "(this is the scope question)")
        for i, line in enumerate(src.splitlines(), 1):
            if re.search(r"(fdrcorrection|multipletests|benjamini|\bbh\b|groupby.*type|"
                         r"for\s+ty|per[-_]type)", line, re.I):
                print(f"      {i:5d}: {line.strip()}")

    print("\n" + "=" * 72)
    print("RE-RUN INVOCATIONS (only needed if G2 says BH was POOLED across types)")
    print("=" * 72)
    print("""
  # Pool B, five classes, `any` rule — reproduces the recorded 596 baseline
  python 35_combined_evidence.py \\
      --types SE,MXE,A3SS,A5SS,RI \\
      --fdr 0.05 --min-dpsi 0.10 --dpsi-rule any \\
      --out combined_evidence_B5.txt

  # Pool B, four classes (Fuwai MXE excluded) — the 526 question
  python 35_combined_evidence.py \\
      --types SE,A3SS,A5SS,RI \\
      --fdr 0.05 --min-dpsi 0.10 --dpsi-rule any \\
      --out combined_evidence_B4.txt

  # Pool A — `discovery`, never `any`.  Ties the run to the record (409 SE, 550).
  python 35_combined_evidence.py \\
      --types SE,A3SS,A5SS,RI \\
      --fdr 0.05 --min-dpsi 0.10 --dpsi-rule discovery \\
      --out combined_evidence_A.txt

  NOTE on --types: the separator is printed above from script 35's own
  add_argument line.  If it uses nargs="+", drop the commas and use spaces.
  Do not guess — the line above is authoritative.
""")


# --------------------------------------------------------------------- main

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default=".", help="directory holding the combined_* files")
    ap.add_argument("--pool-b", default=None, help="per-event Pool B table (default combined_b.tsv)")
    ap.add_argument("--pool-a", default=None, help="per-event Pool A table (default combined_a.tsv)")
    ap.add_argument("--list", action="store_true", help="inventory only")
    args = ap.parse_args()

    if args.list:
        do_list(args)
        return

    pb = args.pool_b or os.path.join(args.dir, "combined_b.tsv")
    if not os.path.exists(pb):
        refuse(f"Pool B per-event table not found: {pb}\n"
               f"  run `python 84_pool_b.py --list --dir {args.dir}` first.")

    B = pd.read_csv(pb, sep="\t")
    c_type = pick(B.columns, TYPE_CANDIDATES, "event-type")
    c_q    = pick(B.columns, Q_CANDIDATES, "combined-FDR")
    c_p    = pick(B.columns, P_CANDIDATES, "combined-p")
    dcols  = cohort_dpsi_cols(B.columns)
    if len(dcols) < 3:
        refuse(f"found {len(dcols)} per-cohort dPSI columns, need 3 "
               f"(Brigham, Fuwai, Stanford).\n"
               f"  detected: {dcols}\n"
               f"  table has: {', '.join(B.columns)}")

    print("=" * 72)
    print("84_pool_b.py — is 526 arithmetic, or does dropping MXE move every row?")
    print("=" * 72)
    print(f"\nPool B table : {pb}  ({len(B):,d} rows)")
    print(f"  type column: {c_type}")
    print(f"  FDR column : {c_q}")
    print(f"  p column   : {c_p}")
    print(f"  cohorts    : {', '.join(f'{k} -> {v}' for k, v in sorted(dcols.items()))}")
    print(f"  rule       : FDR < {FDR_CUT}, |dPSI| >= {DPSI_CUT} in ANY cohort, "
          f"unanimous direction")

    # ---------------------------------------------------------------- G1
    print("\n" + "-" * 72)
    print("G1  does the parse reproduce the recorded five-class Pool B ladder?")
    print("-" * 72)
    tab5 = ladder(B, TYPES5, c_type, c_q, dcols)
    tot5, ok5 = show(tab5, TYPES5, rec=REC_B, title="recomputed from the file:")
    print(f"\n  recorded totals: {REC_B_EVENTS:,d} events, {REC_B_UNAN} unanimous")
    print(f"  recomputed     : {tot5['events']:,d} events, {tot5['unan']} unanimous")
    if not ok5 or tot5["events"] != REC_B_EVENTS or tot5["unan"] != REC_B_UNAN:
        refuse("G1 — the recomputation does not reproduce "
               "RESULT_phase1_combined_evidence.md.\n"
               "  Either the column resolution above is wrong, or this file is from a\n"
               "  different run than the one the record describes. Paste the table\n"
               "  above and the output of `--list`; do not proceed.")
    print("\n  G1 PASS — the file reproduces the record exactly. The parse is correct.")

    # ---------------------------------------------------------------- G2
    print("\n" + "-" * 72)
    print("G2  what was the BH scope — within each event type, or pooled across types?")
    print("-" * 72)
    q_per = np.full(len(B), np.nan)
    t = B[c_type].astype(str).str.upper().to_numpy()
    for ty in TYPES5:
        m = (t == ty.upper())
        q_per[m] = bh(B.loc[m, c_p].to_numpy(dtype=float))
    q_pool = bh(B[c_p].to_numpy(dtype=float))

    hits = {}
    for name, qq in (("per-type", q_per), ("pooled", q_pool)):
        tabx = ladder(B, TYPES5, c_type, c_q, dcols, q=qq)
        agree = all(tabx[ty]["fdr"] == REC_B[ty]["fdr"] for ty in TYPES5)
        hits[name] = agree
        line = "  ".join(f"{ty} {tabx[ty]['fdr']:,d}" for ty in TYPES5)
        print(f"  BH {name:9s}: {line}    "
              f"{'reproduces the record' if agree else 'does NOT reproduce the record'}")
    print("  record     : " + "  ".join(f"{ty} {REC_B[ty]['fdr']:,d}" for ty in TYPES5))

    if not hits["per-type"] and not hits["pooled"]:
        refuse("G2 — neither BH scope reproduces the recorded FDR<0.05 counts, even "
               "though\n  G1 passed on the stored FDR column. The stored q-values were "
               "produced by\n  something other than plain BH over this p column "
               "(filtering before correction,\n  a different p, or a different method). "
               "Paste this output — the scope has to be\n  read off script 35's source "
               "instead, and `--list` prints the relevant lines.")

    if hits["per-type"] and not hits["pooled"]:
        scope = "per-type"
    elif hits["pooled"] and not hits["per-type"]:
        scope = "pooled"
    else:
        scope = "pooled"
        print("\n  Both scopes reproduce the counts — the p-values are separated enough "
              "that the\n  rank change does not cross the threshold for any event. "
              "Taking the conservative\n  branch (pooled) below; if it agrees with 526, "
              "the scope question is moot.")
    print(f"\n  G2 PASS — BH scope identified as: {scope.upper()}")

    # ---------------------------------------------------------------- G3
    print("\n" + "-" * 72)
    print("G3  the four-class Pool B, MXE dropped, with BH applied at that scope")
    print("-" * 72)
    keep = (t != "MXE")
    B4 = B.loc[keep].reset_index(drop=True)
    t4 = B4[c_type].astype(str).str.upper().to_numpy()
    if scope == "per-type":
        q4 = np.full(len(B4), np.nan)
        for ty in TYPES4:
            m = (t4 == ty.upper())
            q4[m] = bh(B4.loc[m, c_p].to_numpy(dtype=float))
    else:
        q4 = bh(B4[c_p].to_numpy(dtype=float))

    tab4 = ladder(B4, TYPES4, c_type, c_q, dcols, q=q4)
    tot4, _ = show(tab4, TYPES4, title=f"four classes, BH {scope}, MXE removed:")

    moved = [ty for ty in TYPES4
             if (tab4[ty]["fdr"], tab4[ty]["effect"], tab4[ty]["unan"]) !=
                (REC_B[ty]["fdr"], REC_B[ty]["effect"], REC_B[ty]["unan"])]
    print(f"\n  rows that moved when MXE was removed: "
          f"{', '.join(moved) if moved else 'NONE'}")
    print(f"  four-class unanimous total: {tot4['unan']}")
    print(f"  currently quoted in Results: {QUOTED_526}")
    print(f"  arithmetic check 363+29+21+113 = "
          f"{sum(REC_B[ty]['unan'] for ty in TYPES4)}")

    if tot4["unan"] == QUOTED_526 and not moved:
        print(f"""
  VERDICT — 526 IS CORRECT AND IS NOW READ OFF A FILE.
  Removing MXE left all four surviving rows untouched, so the subtraction
  596 - 70 = 526 was valid. Pool B's four-class total is {tot4['unan']}
  ({' + '.join(str(tab4[ty]['unan']) for ty in TYPES4)}),
  recomputed from {os.path.basename(pb)}. No cluster re-run is needed.
  -> mark the 526 item in DRAFT_results.md and PAPER_ONE_PAGE.md as sourced.""")
    elif tot4["unan"] == QUOTED_526 and moved:
        print(f"""
  VERDICT — 526 holds, but not for the reason the record gives.
  Rows {', '.join(moved)} changed internally when MXE was removed and the
  total still lands on {tot4['unan']}. Report the per-class numbers from the
  table above, NOT the recorded five-class rows minus MXE.""")
    else:
        cause = (f"BH was applied {scope}, so deleting {int((~keep).sum()):,d} MXE rows "
                 f"re-ranked the\n  remaining p-values and moved {', '.join(moved)}"
                 if moved else
                 "no individual class row moved, so the discrepancy is in the recorded\n"
                 "  total itself, not in the pooling — check the record's arithmetic")
        print(f"""
  VERDICT — {QUOTED_526} IS WRONG. The four-class Pool B total is {tot4['unan']}.
  {cause}.
  Replace {QUOTED_526} with {tot4['unan']} and the per-class numbers above everywhere:
    - DRAFT_results.md       (the "526 across the four classes" sentence)
    - PAPER_ONE_PAGE.md      (section 4, combined evidence)
    - Figure 4 and its legend (the record calls 596 "the F4 number")
    - S2, if it carries the pooled totals
  Then confirm against the real script:
    python 35_combined_evidence.py --types SE,A3SS,A5SS,RI \\
        --fdr 0.05 --min-dpsi 0.10 --dpsi-rule any --out combined_evidence_B4.txt""")

    # ---------------------------------------------------------------- G4
    print("\n" + "-" * 72)
    print("G4  Pool A cross-check — ties these files to the recorded run")
    print("-" * 72)
    pa = args.pool_a or os.path.join(args.dir, "combined_a.tsv")
    if not os.path.exists(pa):
        print(f"  SKIPPED — {pa} not found. G1-G3 stand on their own, but without this\n"
              f"  there is no independent confirmation that these files are the vintage\n"
              f"  the record describes.")
    else:
        A = pd.read_csv(pa, sep="\t")
        try:
            a_type = pick(A.columns, TYPE_CANDIDATES, "event-type")
            a_q    = pick(A.columns, Q_CANDIDATES, "combined-FDR")
        except SystemExit:
            raise
        adc = cohort_dpsi_cols(A.columns)
        if len(adc) < 2:
            print(f"  SKIPPED — Pool A has {len(adc)} per-cohort dPSI columns, need 2 "
                  f"(MAGNet, Colorado).\n  columns: {', '.join(A.columns)}")
        else:
            # Pool A uses the `discovery` rule: |dPSI| >= 0.10 in MAGNet only.
            disc = None
            for k, v in adc.items():
                if re.search(r"magnet|magnt|mag", k, re.I):
                    disc = v
            if disc is None:
                print(f"  cannot identify the MAGNet (discovery) column among "
                      f"{list(adc.values())};\n  falling back to the ANY rule, which will "
                      f"over-count — treat a mismatch here as uninformative.")
                tabA = ladder(A, TYPES4, a_type, a_q, adc)
            else:
                # the `discovery` rule: effect threshold in MAGNet only,
                # unanimous direction across both cohorts
                qv = A[a_q].to_numpy(dtype=float)
                D = A[list(adc.values())].to_numpy(dtype=float)
                with np.errstate(invalid="ignore"):
                    eff = np.abs(A[disc].to_numpy(dtype=float)) >= DPSI_CUT
                    unan = (D > 0).all(axis=1) | (D < 0).all(axis=1)
                sig = np.isfinite(qv) & (qv < FDR_CUT)
                ta = A[a_type].astype(str).str.upper().to_numpy()
                tabA = {}
                for ty in TYPES4:
                    m = (ta == ty.upper())
                    tabA[ty] = dict(events=int(m.sum()),
                                    fdr=int((m & sig).sum()),
                                    effect=int((m & sig & eff).sum()),
                                    unan=int((m & sig & eff & unan).sum()))
                print(f"  discovery (effect) column: {disc}")
            totA, okA = show(tabA, TYPES4, rec=REC_A, title="Pool A, `discovery` rule:")
            print(f"\n  recorded: {REC_A_EVENTS:,d} events, {REC_A_UNAN} unanimous")
            print(f"  recomputed: {totA['events']:,d} events, {totA['unan']} unanimous")
            if okA and totA["events"] == REC_A_EVENTS and totA["unan"] == REC_A_UNAN:
                print("  G4 PASS — Pool A reproduces exactly. Both files are the recorded "
                      "vintage.")
            else:
                print("  G4 FAIL — Pool A does not reproduce. The Pool B conclusion above "
                      "may be\n  built on a file from a different run. Paste this output "
                      "before using it.")

    print("\n" + "=" * 72)
    print("REMINDER: Pool A is the `discovery` run and Pool B the `any` run.")
    print("Never quote one pool's number from the other's file.")
    print("=" * 72)


if __name__ == "__main__":
    main()
