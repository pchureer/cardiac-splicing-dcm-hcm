#!/usr/bin/env python3
"""
92_fix_f6_title_and_ga_constant.py  —  the last two figure edits.

TWO SEPARATE JOBS, ONE SCRIPT, EITHER CAN BE RUN ALONE
------------------------------------------------------
F6  Panel D's title reads "What survives: 15 genes" while panel C's reads
    "...15 terms". Two different fifteens, one glance apart. A reader who skims
    C then D will read them as the same fifteen things. The number in each is
    right; the collision is the problem, so only D's wording changes.

GA  `47_graphical_abstract.py` still holds the eleven-count in CODE. Script 91
    fixed the two DRAWN strings (heading and strap line) and warned that other
    instances survived. Those survivors matter for two reasons: one of them is
    a self-check that prints `[ok]` while asserting the miscount the figure no
    longer has, and all of them ship with the deposit.

    ⚠ 11 IS NOT ALWAYS WRONG IN THAT FILE. Eleven genes are FITTED; ten are
    CANDIDATES; nine are candidates OTHER THAN RBM20. A constant that counts
    fitted genes should stay 11. This script therefore changes `11` only where
    the adjacent word is candidate/regulator, and it REFUSES to guess anywhere
    else — it lists those lines for you to read.

USAGE
    python 92_fix_f6_title_and_ga_constant.py --list          # look first
    python 92_fix_f6_title_and_ga_constant.py                 # dry run: diffs
    python 92_fix_f6_title_and_ga_constant.py --apply --only f6
    python 92_fix_f6_title_and_ga_constant.py --apply --only ga
    python 92_fix_f6_title_and_ga_constant.py --apply         # both

    --new-title 'your wording'   override the default panel-D title

GATES
  G1  each target file exists.
  G2  the string to be replaced matches EXACTLY once. Zero means the wording
      differs from the rendered PNG (find it by hand); more than once means an
      edit would hit something unintended. Either way: refuse, write nothing.
  G3  F6 only: panel C's title is found and still contains its own number, so
      the collision being fixed is the one that actually exists.
  G4  after writing, reload from disk and confirm the new string is present and
      the old one is gone. Residual eleven-counts are REPORTED WITH LINE
      NUMBERS, never treated as failure — the file is correct as far as this
      script wrote it, and restoring a backup over good work is worse than a
      warning. (That mistake was made twice in this project.)
"""

import argparse
import os
import re
import shutil
import sys
import time

F6 = "41_figure6_enrichment.py"
GA = "47_graphical_abstract.py"

# ---------------------------------------------------------------- F6 ----------
F6_OLD = "What survives: 15 genes"
F6_NEW_DEFAULT = "Survives adjustment: 15 genes"
# panel C's title, used only to prove the collision is real
F6_SIBLING = re.compile(r"15\s+terms")

# ---------------------------------------------------------------- GA ----------
# `11` immediately beside the word candidate/regulator — safe to change to 10.
GA_SAFE = [
    (re.compile(r"\b11\s+candidate\b"),                 "10 candidate"),
    (re.compile(r"\beleven\s+candidate\b", re.I),       "ten candidate"),
    (re.compile(r"\b11\s+(?=regulators?\b)"),           "10 "),
    (re.compile(r"\beleven\s+(?=regulators?\b)", re.I), "ten "),
    (re.compile(r"\bten\s+other\s+candidate\b"),        "nine other candidate"),
    (re.compile(r"\bten\s+other\s+regulator"),          "nine other regulator"),
]
# every other place the number 11 / eleven appears — reported, never touched
GA_SURVEY = re.compile(r"\b(11|eleven)\b", re.I)
# a self-check that may be asserting the miscount
GA_CHECK = re.compile(r"(assert|\[ok\]|check|expect|n_cand|N_CAND|n_genes|N_GENES)",
                      re.I)

# A line that BINDS the number — an assert, a comparison, a constant assignment.
# On such a line the label and the number must change TOGETHER or not at all:
# rewriting `assert N_CAND == 11, "expected 11 candidate regulators"` into
# `... == 11, "expected 10 ..."` leaves the wrong number in force and makes it
# HARDER to spot. The test run of this script did exactly that. So: defer.
GA_BINDING = re.compile(r"(\bassert\b|[=!<>]=|^\s*[A-Za-z_][A-Za-z_0-9]*\s*=\s*(11|eleven)\b)",
                        re.M)


def hr(c="-"):
    print(c * 78)


def refuse(msg):
    print()
    hr("=")
    print("REFUSING — " + msg)
    hr("=")
    sys.exit(2)


def lineno(src, pos):
    return src.count("\n", 0, pos) + 1


def ctx(src, m, back=60, fwd=45):
    s = src[max(0, m.start() - back):m.end() + fwd].replace("\n", " ⏎ ")
    return re.sub(r"\s+", " ", s).strip()


def whole_line(src, pos):
    a = src.rfind("\n", 0, pos) + 1
    b = src.find("\n", pos)
    return src[a:(b if b != -1 else len(src))]


def read(path):
    if not os.path.exists(path):
        return None
    return open(path, encoding="utf-8").read()


# ============================================================ F6 =============
def do_f6(d, new_title, apply_, list_only):
    p = os.path.join(d, F6)
    print()
    hr("=")
    print(f"F6 — panel D title   ({F6})")
    hr("=")
    src = read(p)
    if src is None:
        print(f"  G1 FAIL — {F6} not found in {d}. Skipping F6.")
        return None
    print(f"  file: {p}  ({len(src):,} characters)")

    hits = list(re.finditer(re.escape(F6_OLD), src))
    print(f"\n  G2  {F6_OLD!r}: {len(hits)} match(es), expected 1")
    for h in hits:
        print(f"        line {lineno(src, h.start())}: …{ctx(src, h)}…")

    sib = list(F6_SIBLING.finditer(src))
    print(f"\n  G3  panel C's 15-terms wording: {len(sib)} match(es)")
    for h in sib:
        print(f"        line {lineno(src, h.start())}: …{ctx(src, h)}…")

    # show every set_title in the file, so you can see the whole panel strip
    titles = list(re.finditer(r"set_title\(\s*(?:f?[\"'])(.{0,90}?)[\"']", src))
    if titles:
        print(f"\n  every set_title in the file ({len(titles)}):")
        for t in titles:
            print(f"        line {lineno(src, t.start()):5d}  {t.group(1)!r}")

    if list_only:
        return None

    if len(hits) != 1:
        print(f"\n  SKIPPING F6 — expected exactly 1 match for the old title, "
              f"found {len(hits)}.")
        print(f"  Find it by hand:  grep -n 'survives\\|set_title' {F6}")
        print("  Then re-run with --new-title and the exact old wording edited "
              "into this\n  script's F6_OLD, or make the one-line change yourself.")
        return False
    if not sib:
        print("\n  ⚠ G3 — panel C's '15 terms' wording was not found. The collision "
              "this\n    edit exists to fix may already be gone, or panel C is worded "
              "differently.\n    The edit is still safe (it only reworks panel D), but "
              "LOOK at the PNG\n    afterwards and decide whether it was needed.")

    new = src.replace(F6_OLD, new_title)
    hr()
    print("  THE CHANGE")
    hr()
    print(f"    {F6_OLD!r}")
    print(f"      ->  {new_title!r}")
    print("\n    The number does not change. 15 genes survive adjustment; that is "
          "right.\n    What changes is that panel D no longer opens with the same two "
          "words as\n    panel C, so the two fifteens cannot be read as one.")

    if not apply_:
        print("\n  DRY RUN — nothing written.")
        return None

    bak = f"{p}.bak-{time.strftime('%Y%m%d-%H%M%S')}"
    shutil.copy2(p, bak)
    open(p, "w", encoding="utf-8").write(new)
    print(f"\n  [backup] {bak}\n  [write]  {p}")

    chk = read(p)
    present, gone = new_title in chk, F6_OLD not in chk
    print(f"\n  G4  new title present={present}   old title gone={gone}")
    if not (present and gone):
        refuse(f"G4 — the F6 edit did not apply cleanly. Restore from {bak}.")
    print("  G4 PASS")
    return True


# ============================================================ GA =============
def do_ga(d, apply_, list_only):
    p = os.path.join(d, GA)
    print()
    hr("=")
    print(f"GA — the eleven-count still in code   ({GA})")
    hr("=")
    src = read(p)
    if src is None:
        print(f"  G1 FAIL — {GA} not found in {d}. Skipping GA.")
        return None
    print(f"  file: {p}  ({len(src):,} characters)")

    # ---- what this script will change, decided MATCH BY MATCH
    spans, deferred, total = [], [], 0
    print("\n  G2  patterns this script is willing to change "
          "(11/eleven beside candidate|regulator)")
    for rx, rep in GA_SAFE:
        hits = list(rx.finditer(src))
        total += len(hits)
        print(f"    {rx.pattern:44s} {len(hits)} match(es)")
        for h in hits:
            ln = lineno(src, h.start())
            line = whole_line(src, h.start())
            if GA_BINDING.search(line):
                deferred.append((ln, line.strip()))
                print(f"        line {ln}: DEFERRED — the number is bound on this line")
            else:
                spans.append((h.start(), h.end(), rep, rx.pattern))
                print(f"        line {ln}: …{ctx(src, h)}…")
    if total == 0:
        print("\n    none — script 91's two edits are the only ones that were needed,")
        print("    or they have already been applied to the code as well as the text.")

    if deferred:
        hr()
        print(f"  {len(deferred)} LINE(S) NEED A HAND EDIT — this script will not touch them")
        hr()
        for ln, line in deferred:
            print(f"    line {ln}:  {line}")
        print("""
    Each of these BINDS the number: an assert, a comparison, or a constant set
    to 11. Rewriting only the words leaves the wrong number in force and makes
    it harder to spot, not easier — a check reading

        assert N_CAND == 11, "expected 10 candidate regulators"

    passes while contradicting itself. The number and its label have to move
    together. Decide what each constant actually counts:

        fitted genes      11   (10 candidates + GAPDH)      -> leave at 11
        candidates        10                                -> change to 10
        candidates ≠ RBM20 9                                -> change to 9

    then edit the constant AND its message in the same stroke.""")

    # ---- what it will NOT change
    seen = set()
    print("\n  SURVEY  every line holding '11' or 'eleven'")
    print("          (reported only — this script changes none of these directly)")
    for m in GA_SURVEY.finditer(src):
        ln = lineno(src, m.start())
        if ln in seen:
            continue
        seen.add(ln)
        line = whole_line(src, m.start()).strip()
        flag = "   <- self-check" if GA_CHECK.search(line) else ""
        print(f"        line {ln:5d}: {line[:110]}{flag}")
    print(f"          {len(seen)} line(s)")

    print("""
  HOW TO READ THAT LIST
    11 genes are FITTED      — a constant counting fitted genes is CORRECT at 11
    10 are CANDIDATES        — GAPDH is the floor, not a candidate
     9 are candidates OTHER THAN RBM20
    So `n_fit = 11` stays. `n_candidates = 11` is wrong. A check that prints
    [ok] against 11 candidates is the one that has been lying, and it will keep
    printing [ok] until someone reads it.""")

    if list_only:
        return None

    if not spans:
        print("\n  Nothing for this script to write in the GA — every match either "
              "does not\n  exist or is on a line that needs a hand edit (above).")
        return None

    hr()
    print("  THE CHANGE")
    hr()
    new = src
    for start, end, rep, pat in sorted(spans, reverse=True):   # back to front
        print(f"    line {lineno(src, start):5d}  {src[start:end]!r}  ->  {rep!r}")
        new = new[:start] + rep + new[end:]

    if not apply_:
        print("\n  DRY RUN — nothing written.")
        return None

    bak = f"{p}.bak-{time.strftime('%Y%m%d-%H%M%S')}"
    shutil.copy2(p, bak)
    open(p, "w", encoding="utf-8").write(new)
    print(f"\n  [backup] {bak}\n  [write]  {p}")

    chk = read(p)
    print("\n  G4  reload and verify")
    if chk != new:
        refuse(f"G4 — the file on disk is not what was written. Restore from {bak}.")
    print(f"    {len(spans)} replacement(s) written, file reloads identical to the "
          f"intended text")
    print("  G4 PASS")

    rest = sorted({lineno(chk, m.start()) for m in GA_SURVEY.finditer(chk)})
    if rest:
        print(f"\n  ⚠ 11/eleven still appears on {len(rest)} line(s) — and the file is "
              f"CORRECT\n    as far as this script wrote it. DO NOT RESTORE THE BACKUP. "
              f"Lines: "
              f"{', '.join(map(str, rest))}")
        if deferred:
            n = len(deferred)
            print(f"    {n} of them {'is' if n == 1 else 'are'} a bound line listed "
                  f"above, still needing the number\n    and its label changed "
                  f"together.")
    return True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default=".")
    ap.add_argument("--only", choices=["f6", "ga"], help="run one job only")
    ap.add_argument("--new-title", default=F6_NEW_DEFAULT,
                    help=f"panel D title (default: {F6_NEW_DEFAULT!r})")
    ap.add_argument("--list", action="store_true", dest="list_only",
                    help="report what is there and change nothing")
    ap.add_argument("--apply", action="store_true")
    a = ap.parse_args()

    if a.list_only and a.apply:
        refuse("--list and --apply are contradictory. Pick one.")

    hr("=")
    print("92 — F6 panel D title · the GA's eleven-count in code")
    hr("=")
    print(f"  dir   {os.path.abspath(a.dir)}")
    print(f"  mode  {'--list (read only)' if a.list_only else ('--apply' if a.apply else 'dry run')}")
    if a.only:
        print(f"  only  {a.only}")

    r6 = rg = None
    if a.only in (None, "f6"):
        r6 = do_f6(a.dir, a.new_title, a.apply, a.list_only)
    if a.only in (None, "ga"):
        rg = do_ga(a.dir, a.apply, a.list_only)

    print()
    hr("=")
    print("NEXT")
    hr("=")
    if a.list_only:
        print("""
  Read the lines above, then:
      python 92_fix_f6_title_and_ga_constant.py            # dry run, see the diffs
      python 92_fix_f6_title_and_ga_constant.py --apply
""")
        return
    if not a.apply:
        print("""
  DRY RUN — nothing on disk changed. When the diffs read right:
      python 92_fix_f6_title_and_ga_constant.py --apply
""")
        return

    print(f"""
  conda activate rnaseq_tools

  1. rebuild F6 and LOOK at it — the only check that has ever caught a wrong
     panel in this project was looking at the picture:
         python {F6} --out legend_check/f6_new
         # panels C and D should no longer open with the same two words

  2. when F6 is right, copy it over the published file:
         cp legend_check/f6_new.png F6_enrichment.png
         cp legend_check/f6_new.pdf F6_enrichment.pdf   # if the script writes one

  3. prove the whole set still reproduces:
         bash 89_make_all_figures.sh
     f6 is EXPECTED to report DIFFERS until step 2 is done, and IDENTICAL after.
     Every other figure must stay IDENTICAL — if one changes, stop.

  4. the GA's rendered text did not change in this edit{" (code only)" if rg else ""}, so
     GA_graphical_abstract.png should still be byte-identical. If 89 reports the
     GA as DIFFERS, a drawn string was touched — rebuild it and look before
     publishing anything.

  5. update FIGURE_LEGENDS.md: panel D's title in the F6 legend.
""")


if __name__ == "__main__":
    main()
