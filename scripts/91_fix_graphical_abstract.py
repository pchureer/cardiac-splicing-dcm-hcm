#!/usr/bin/env python3
"""
91_fix_graphical_abstract.py  —  the two string edits the graphical abstract needs,
                                 applied safely to `47_graphical_abstract.py`.

WHAT IS WRONG, CONFIRMED FROM THE RENDERED PNG
----------------------------------------------
The graphical abstract counts the negative control as one of the things it calibrates:

  heading     "11 candidate regulators"             -> TEN candidates; GAPDH is the floor
  strap line  "...or ten other candidate regulators" -> NINE other, besides RBM20

Eleven genes are fitted. Ten are candidates. One of the ten is RBM20, so there are
**nine others**. This is the same miscount that `RESULT_s7_built.md` found in the
supplement, and the graphical abstract is the last place it survives.

Everything else in the abstract checks out against Table 3 and the record: 510 hearts,
27 excluded, 926 / 468 / 0, 476 / 407 / 2, and RBM20 at +0.34 SD, P = 0.003 — which is
the composition-adjusted model the Results quote, and therefore correct.

WHY A SCRIPT
A two-string edit is easy to get wrong in the direction that matters: changing one and
not the other leaves the figure internally inconsistent (ten candidates in the heading,
ten others in the strap line, which cannot both be true). This refuses unless it can
make both changes, backs up first, and prints the diff.

    python 91_fix_graphical_abstract.py              # dry run
    python 91_fix_graphical_abstract.py --apply      # back up, edit, verify

GATES
  G1  the target file exists and contains the expected strings.
  G2  each pattern matches EXACTLY the expected number of times. Zero means the
      wording differs from the PNG and must be found by hand; more than expected
      means an edit would hit something unintended. Either way: refuse.
  G3  after writing, reload and confirm every new string is present and no old one
      survives — including a scan for any other "eleven"/"11 candidate"/"ten other"
      wording the two patterns did not cover.
"""

import argparse
import os
import re
import shutil
import sys
import time

TARGET = "47_graphical_abstract.py"

# (description, regex, replacement, expected number of matches)
EDITS = [
    ("heading: eleven -> ten candidates",
     re.compile(r"\b11\s+candidate\s+regulators\b"),
     "10 candidate regulators", 1),
    ("strap line: ten other -> nine other",
     re.compile(r"\bten\s+other\s+candidate\s+regulators\b"),
     "nine other candidate regulators", 1),
]

# wording that must not survive anywhere in the file
RESIDUAL = re.compile(
    r"\b(eleven\s+candidate|11\s+candidate|ten\s+other\s+candidate"
    r"|eleven\s+regulators|11\s+regulators)\b", re.I)


def refuse(msg):
    print("\n" + "=" * 78)
    print("REFUSING — " + msg)
    print("=" * 78)
    sys.exit(2)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default=".")
    ap.add_argument("--file", default=TARGET)
    ap.add_argument("--apply", action="store_true")
    a = ap.parse_args()

    p = os.path.join(a.dir, a.file)
    if not os.path.exists(p):
        refuse(f"{a.file} not found in {a.dir}")
    src = open(p, encoding="utf-8").read()

    print("=" * 78)
    print("91 — graphical abstract: ten candidates, nine other than RBM20")
    print("=" * 78)
    print(f"\nfile: {p}  ({len(src):,} characters)")

    # ---------------------------------------------------------------- G1/G2
    print("\nG1/G2  pattern matches")
    plan = []
    bad = []
    for desc, rx, rep, want in EDITS:
        hits = list(rx.finditer(src))
        print(f"  {desc}")
        print(f"    pattern {rx.pattern}")
        print(f"    matches {len(hits)} (expected {want})")
        for h in hits:
            line = src.count("\n", 0, h.start()) + 1
            print(f"      line {line}: …{src[max(0,h.start()-45):h.end()+25]!r}…")
        if len(hits) != want:
            bad.append(f"{desc}: {len(hits)} matches, expected {want}")
        else:
            plan.append((desc, rx, rep))
    if bad:
        refuse("G2 — the file does not read as expected:\n  - " + "\n  - ".join(bad)
               + "\n\n  Zero matches means the wording differs from the rendered PNG — "
                 "find it by\n  hand with:  grep -n 'candidate\\|regulator' "
               + a.file +
               "\n  More than expected means an edit would hit something unintended.\n"
                 "  Both changes must be made together: ten candidates in the heading "
                 "and\n  ten others in the strap line cannot both be true.")

    print("\n" + "-" * 78)
    print("THE CHANGE")
    print("-" * 78)
    new = src
    for desc, rx, rep in plan:
        before = new
        new = rx.sub(rep, new)
        m = rx.search(before)
        print(f"  {desc}")
        print(f"    {before[m.start():m.end()]!r}  ->  {rep!r}")

    others = [m.group(0) for m in RESIDUAL.finditer(new)]
    if others:
        print(f"\n  ⚠ other wording that would survive these two edits: {others}")
        print("    Check each by hand — the count appears in more places than the PNG "
              "shows.")

    if not a.apply:
        print("\n  DRY RUN — nothing written. Re-run with --apply.")
        return

    stamp = time.strftime("%Y%m%d-%H%M%S")
    bak = f"{p}.bak-{stamp}"
    shutil.copy2(p, bak)
    with open(p, "w", encoding="utf-8") as fh:
        fh.write(new)
    print(f"\n[backup] {bak}\n[write]  {p}")

    # ---------------------------------------------------------------- G3
    print("\n" + "-" * 78)
    print("G3  reload and verify")
    print("-" * 78)
    chk = open(p, encoding="utf-8").read()
    ok = True
    for desc, rx, rep, _ in EDITS:
        present = rep in chk
        gone = not rx.search(chk)
        print(f"  {rep!r:38s} present={present}  old gone={gone}")
        ok &= present and gone
    left = [m.group(0) for m in RESIDUAL.finditer(chk)]
    print(f"  residual miscount wording: {left or 'none'}")
    ok &= not left
    if not ok:
        refuse(f"G3 — the written file is not what was intended. Restore from {bak}.")
    print("\n  G3 PASS")

    print(f"""
NEXT
  conda activate rnaseq_tools
  python {a.file} --out legend_check/ga_new      # then LOOK at the PNG
  # when it is right, copy over GA_graphical_abstract.png / .pdf

  The abstract's other numbers were checked against Table 3 and need no change:
  510 hearts · 27 excluded · 926 / 468 / 0 · 476 / 407 / 2 · RBM20 +0.34 SD, P = 0.003
  (that last is the composition-adjusted model, which is the one the Results quote —
  Figure 5A showed a different model and is being rebuilt by 90_build_f5.py).
""")


if __name__ == "__main__":
    main()
