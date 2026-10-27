#!/usr/bin/env python3
"""
88_fix_fuwai_mxe_flag.py  —  apply ONE decision to the frozen master numbers file:
                             Fuwai MXE reliability `caution` -> `exclude`.

WHY THIS EXISTS AS A SCRIPT AND NOT A HAND EDIT
-----------------------------------------------
`master_numbers_flagged.tsv` is the project's single source of truth — "nothing enters
the manuscript that is not in that file" — and it ships with the code deposit. Editing
it by hand leaves no record of what changed or why, and a one-character slip in a frozen
file is the worst kind of error to chase later. So: backup, gate, show the diff, and
write only with `--apply`.

THE DECISION BEING APPLIED
`DECISION_fuwai_mxe_and_permutations.md`, taken 6 Oct 2026 by CP: Fuwai MXE is excluded
on the same grounds as Colorado MXE. It is worse on both selection-free estimators
(r = 0.088 unselected and 0.111 precision-filtered, against Colorado's 0.157 and 0.198)
and carries the same QC flag. That doc's "what must change" list includes this file
edit; it was never applied, and `87_build_table3.py`'s G2 found it.

WHAT CHANGES, AND WHAT DOES NOT
- **No count changes.** This edits a `reliability` label, not a number.
- It does change what any script that respects the label will sum: S2's Fuwai MXE row,
  the all-five-types pooled concordance (2,227/2,991 -> 2,167/2,900), and Pool B's MXE
  row (already dropped in the manuscript).
- It makes the file agree with the manuscript, which is the point.

USAGE
    python 88_fix_fuwai_mxe_flag.py                 # dry run: show the diff, write nothing
    python 88_fix_fuwai_mxe_flag.py --apply         # back up, then write

GATES
  G1  the file has the needed columns, and exactly one Fuwai MXE row exists.
  G2  that row's reliability is currently `caution` — if it is already `exclude` the job
      is done and the script says so and stops; if it is anything else, REFUSE, because
      an unexpected starting value means this is not the file the decision describes.
  G3  report EVERY row still carrying `caution`, so that fixing one does not hide others.
  G4  after writing, reload from disk and verify the new value and that nothing else
      moved: row count, column set, and a checksum of every other cell.
"""

import argparse
import hashlib
import os
import shutil
import sys
import time

import pandas as pd

MASTER = "master_numbers_flagged.tsv"
TARGET = ("Fuwai", "MXE")
FROM, TO = "caution", "exclude"
REASON = ("excluded 6 Oct 2026: lowest MXE r in the study (0.088 unselected, 0.111 "
          "precision-filtered); see DECISION_fuwai_mxe_and_permutations.md")


def refuse(msg):
    print("\n" + "=" * 78)
    print("REFUSING — " + msg)
    print("=" * 78)
    sys.exit(2)


def fingerprint(df, skip=()):
    """Checksum every cell EXCEPT the intentionally changed ones, to prove nothing else
    moved. `skip` is a set of (index, column) pairs — it must list every cell this
    script writes, or the check raises a false alarm on its own correct edit."""
    skip = set(skip)
    h = hashlib.sha256()
    for i, row in df.iterrows():
        for c in df.columns:
            if (i, c) in skip:
                continue
            h.update(f"{i}|{c}|{row[c]!r}".encode())
    return h.hexdigest()[:16]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default=".")
    ap.add_argument("--apply", action="store_true", help="write the change")
    args = ap.parse_args()

    p = os.path.join(args.dir, MASTER)
    if not os.path.exists(p):
        refuse(f"{MASTER} not found in {args.dir}")
    df = pd.read_csv(p, sep="\t", dtype=str, keep_default_na=False)

    print("=" * 78)
    print("88 — Fuwai MXE: reliability `caution` -> `exclude`")
    print("=" * 78)
    print(f"\nfile: {p}   {len(df)} rows, {len(df.columns)} columns")

    # ---------------------------------------------------------------- G1
    for c in ("cohort", "event_type", "reliability"):
        if c not in df.columns:
            refuse(f"{MASTER} has no `{c}` column.\n  has: {', '.join(df.columns)}")
    has_reason = "flag_reason" in df.columns

    m = (df.cohort.str.split().str[0] == TARGET[0]) & \
        (df.event_type.str.upper() == TARGET[1])
    idx = list(df.index[m])
    print(f"\nG1  {TARGET[0]} {TARGET[1]} rows found: {len(idx)}")
    if len(idx) != 1:
        refuse(f"G1 — expected exactly one {TARGET[0]} {TARGET[1]} row, found "
               f"{len(idx)}.\n  Rows: {idx}\n  Resolve by hand; this script will not "
               f"guess which to change.")
    i = idx[0]
    show = ["cohort", "contrast", "event_type", "reliability"] + \
           (["flag_reason"] if has_reason else [])
    print("  current row:")
    for c in show:
        print(f"    {c:14s} {df.at[i, c]!r}")

    # ---------------------------------------------------------------- G2
    cur = str(df.at[i, "reliability"]).strip()
    print(f"\nG2  reliability is {cur!r}")
    if cur == TO:
        print(f"  Already `{TO}`. Nothing to do — the decision is applied.")
        return
    if cur != FROM:
        refuse(f"G2 — expected {FROM!r} before the change, found {cur!r}.\n"
               f"  The decision describes a `caution` row. An unexpected value means "
               f"this is\n  not that file, or something else already edited it. Do not "
               f"overwrite blind.")
    print(f"  G2 PASS — will become {TO!r}")

    # ---------------------------------------------------------------- G3
    still = df[df.reliability.str.strip() == FROM]
    print(f"\nG3  every row still carrying {FROM!r}: {len(still)}")
    for j, r in still.iterrows():
        mark = "  <- this one" if j == i else ""
        print(f"    {r['cohort']:10s} {r.get('contrast', ''):28s} "
              f"{r['event_type']:5s}{mark}")
    if len(still) > 1:
        print(f"  ⚠ {len(still) - 1} other `caution` row(s) remain after this change.\n"
              f"    Each is a judgement the manuscript must state or retract. Check "
              f"them;\n    this script deliberately changes only the one the decision "
              f"names.")

    # ---------------------------------------------------------------- diff
    # every cell this script writes must be listed, or G4 flags its own edit
    CHANGED = {(i, "reliability")} | ({(i, "flag_reason")} if has_reason else set())
    before_fp = fingerprint(df, skip=CHANGED)
    new_reason = df.at[i, "flag_reason"] if has_reason else None
    if has_reason:
        old = str(new_reason).strip()
        new_reason = f"{old} | {REASON}" if old and old.lower() != "nan" else REASON

    print("\n" + "-" * 78)
    print("THE CHANGE")
    print("-" * 78)
    print(f"  row {i}: reliability  {cur!r}  ->  {TO!r}")
    if has_reason:
        print(f"  row {i}: flag_reason  {str(df.at[i, 'flag_reason'])!r}")
        print(f"                   ->  {new_reason!r}")
    print(f"\n  every cell this script does NOT write, checksum: {before_fp}")

    if not args.apply:
        print("\n  DRY RUN — nothing written. Re-run with --apply to make the change.")
        return

    # ---------------------------------------------------------------- write
    stamp = time.strftime("%Y%m%d-%H%M%S")
    bak = f"{p}.bak-{stamp}"
    shutil.copy2(p, bak)
    print(f"\n[backup] {bak}")

    df.at[i, "reliability"] = TO
    if has_reason:
        df.at[i, "flag_reason"] = new_reason
    df.to_csv(p, sep="\t", index=False)
    print(f"[write]  {p}")

    # ---------------------------------------------------------------- G4
    print("\n" + "-" * 78)
    print("G4  reload from disk and verify")
    print("-" * 78)
    re_df = pd.read_csv(p, sep="\t", dtype=str, keep_default_na=False)
    ok = True
    if len(re_df) != len(df):
        print(f"  ⚠ row count {len(re_df)} != {len(df)}")
        ok = False
    if list(re_df.columns) != list(df.columns):
        print("  ⚠ column set changed")
        ok = False
    got = str(re_df.at[i, "reliability"]).strip()
    print(f"  reliability on row {i}: {got!r}   {'ok' if got == TO else '⚠ WRONG'}")
    ok &= got == TO
    after_fp = fingerprint(re_df, skip=CHANGED)
    same = after_fp == before_fp
    print(f"  every cell this script does NOT write, checksum: {after_fp}   "
          f"{'unchanged' if same else '⚠ CHANGED'}")
    ok &= same
    if not ok:
        refuse(f"G4 — the written file is not what was intended.\n"
               f"  Restore from {bak} and investigate before doing anything else.")
    print("\n  G4 PASS — one label changed, nothing else touched.")
    print("\n" + "=" * 78)
    print("NEXT")
    print("=" * 78)
    print(f"""
  1. python 87_build_table3.py            # G2 should now pass
  2. Anything that sums Fuwai's classes changes with this label:
       - S2's Fuwai MXE row must read excluded, with the reason
       - the all-five-types pooled concordance is 2,167/2,900 = 74.7%
         (-> RESULT_combination_recount_and_sidedness.md)
       - Pool B reports no MXE row (already applied in the manuscript)
  3. The deposit ships this file. Note in the README that `reliability` carries the
     exclusion decisions, and that `caution` is NOT the same as `exclude`.
  4. Keep {os.path.basename(bak)} until the submission is final.
""")


if __name__ == "__main__":
    main()
