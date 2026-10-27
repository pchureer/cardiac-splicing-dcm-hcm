#!/usr/bin/env python3
"""
169_patch_43_figure1_slopes.py -- Figure 1 draws the matched-pair slopes as
LABEL TEXT, and they are 50 bp.

43_figure1_design.py, lines 101-105:

    ("MAGNet → LBNL",   "same 24 individuals, 10 DCM vs 14 NF   → slope 0.584",
    ("MAGNet → MAGNet", "different individuals, 10 DCM vs 14 NF   → slope 0.988",

At 125 bp those are 0.690 and 0.968. They are literals inside CONTRASTS, so
--out can rename the file while the figure keeps the old numbers.

⚠ THIS HAS ALREADY HAPPENED ONCE. Before this patch existed, a run of

    F1_SLOPES=0.690/0.968 python3 43_figure1_design.py --out F1_design_rl125

wrote F1_design_rl125.png containing 0.584 and 0.988 -- because nothing read
F1_SLOPES. The script printed [write] and PASSED its own consistency check,
which compares the figure against Table1.tsv and so cannot see a slope.

> An environment variable that no code reads is not a weak control. It is the
> absence of one, wearing the appearance of a strong one.

WHICH ESTIMATES
---------------
0.690 and 0.968 are 36_matched_pair_analysis.py's, the ones the Results quote.
82 gives the same point estimates with slightly different bootstrap intervals.
DEFECT_matched_panel_rl125_has_no_file.md: quote 36's in the manuscript and
Figure 1, 82's in S4, and never mix them.

THE GUARD IS SYMMETRIC
----------------------
* "rl125" on the command line without F1_SLOPES -> refuse.
* F1_SLOPES without "rl125" -> refuse, because that writes 125 bp slopes into
  F1_design.png and silently replaces the published 50 bp figure.

USAGE
-----
  python3 169_patch_43_figure1_slopes.py [--apply] [--self-test]
"""
import argparse
import difflib
import os
import shutil
import sys
import time

TARGET = "43_figure1_design.py"

A1 = "CONTRASTS = ["
B1 = '''# ---- 169: the matched-pair slopes are PLOTTED LABEL TEXT, not filenames ----
# Rule 10. --out renames the file; it does not touch these two numbers, so a
# 125 bp filename can carry 50 bp slopes with every check passing. The override
# is therefore MANDATORY, and symmetric: a 125 bp name without the slopes is
# refused, and the slopes without a 125 bp name are refused too, because that
# would overwrite the published 50 bp figure in place.
#
# 0.690/0.968 are 36_matched_pair_analysis.py's estimates, which the Results
# quote. 82 gives the same points with different bootstrap intervals; the two
# must not be mixed (DEFECT_matched_panel_rl125_has_no_file.md).
import os as _os
import sys as _sys

_F1S = _os.environ.get("F1_SLOPES", "").strip()
_F1_RL125 = "rl125" in " ".join(_sys.argv)
if _F1_RL125 and not _F1S:
    _sys.exit(
        "[fatal] a 125 bp output was named, but F1_SLOPES is unset.\\n"
        "        Figure 1 draws the matched-pair slopes as LABEL TEXT, so the\\n"
        "        figure would carry 0.584 and 0.988 under a 125 bp filename.\\n"
        "        Re-run as:\\n"
        "          F1_SLOPES=0.690/0.968 python3 43_figure1_design.py \\\\\\n"
        "              --out F1_design_rl125\\n"
        "        Use 36's estimates, which the Results quote, not 82's.\\n"
        "        Nothing plotted.")
if _F1S and not _F1_RL125:
    _sys.exit(
        "[fatal] F1_SLOPES is set but no 125 bp output was named. The 125 bp\\n"
        "        slopes would be drawn into F1_design, replacing the published\\n"
        "        50 bp figure in place. Pass --out F1_design_rl125.\\n"
        "        Nothing plotted.")
if _F1S:
    _F1_PLATFORM, _F1_PEOPLE = (float(_x) for _x in _F1S.split("/"))
else:
    _F1_PLATFORM, _F1_PEOPLE = 0.584, 0.988
print(f"[169] matched-pair slopes drawn: platform {_F1_PLATFORM:.3f}, "
      f"people {_F1_PEOPLE:.3f}"
      f"{'  (125 bp)' if _F1S else '  (50 bp defaults)'}")

CONTRASTS = ['''

A2 = '''        ("MAGNet → LBNL", "same 24 individuals, 10 DCM vs 14 NF   → slope 0.584",
         "prep"),
        ("MAGNet → MAGNet", "different individuals, 10 DCM vs 14 NF   → slope 0.988",
         "prep"),'''
B2 = '''        ("MAGNet → LBNL",
         f"same 24 individuals, 10 DCM vs 14 NF   → slope {_F1_PLATFORM:.3f}",
         "prep"),
        ("MAGNet → MAGNet",
         f"different individuals, 10 DCM vs 14 NF   → slope {_F1_PEOPLE:.3f}",
         "prep"),'''

A3 = "    # reported 0.584 is the matched-pair panel: the 24 shared individuals, split"
B3 = ("    # reported platform slope is the matched-pair panel: the 24 shared\n"
      "    # individuals, split")

EDITS = [("the override block", A1, B1),
         ("the two slope labels", A2, B2),
         ("the stale 0.584 in the comment", A3, B3)]


def apply_edits(text):
    if "169:" in text:
        return text, ["43 already carries 169 (restore the .bak to re-apply)"]
    probs = [f"{n}: anchor occurs {text.count(a)}x, expected 1x"
             for n, a, b in EDITS if text.count(a) != 1]
    if probs:
        return text, probs
    for n, a, b in EDITS:
        text = text.replace(a, b, 1)
    return text, []


def verify(text):
    bad = []
    for w in ('_F1S = _os.environ.get("F1_SLOPES", "").strip()',
              '_F1_RL125 = "rl125" in " ".join(_sys.argv)',
              'F1_SLOPES is unset',
              'F1_SLOPES is set but no 125 bp output',
              '_F1_PLATFORM, _F1_PEOPLE = 0.584, 0.988',
              'slope {_F1_PLATFORM:.3f}',
              'slope {_F1_PEOPLE:.3f}'):
        if w not in text:
            bad.append(f"missing: {w[:52]}")
    for w in ("→ slope 0.584", "→ slope 0.988"):
        if w in text:
            bad.append(f"a hard-coded slope label survived: {w}")
    if "CONTRASTS = [" not in text:
        bad.append("CONTRASTS was damaged")
    return bad


FIX = '''COHORTS = []

    # The LBNL 5 vs 4 contrast is NOT where the platform result comes from. The
    # reported 0.584 is the matched-pair panel: the 24 shared individuals, split
    # 10 DCM / 14 NF, run three ways.

CONTRASTS = [
    ("Discovery", [
        ("MAGNet", "DCM 166  vs  NF 166", "disc"),
    ]),
    ("Platform vs sampling — the matched-pair panel", [
        ("MAGNet → LBNL", "same 24 individuals, 10 DCM vs 14 NF   → slope 0.584",
         "prep"),
        ("MAGNet → MAGNet", "different individuals, 10 DCM vs 14 NF   → slope 0.988",
         "prep"),
    ]),
]
'''


def self_test():
    f = []
    new, probs = apply_edits(FIX)
    if probs:
        f.append(f"fixture did not patch: {probs}")
    else:
        bad = verify(new)
        if bad:
            f.append(f"verify rejected its own output: {bad}")
        try:
            compile(new, "x", "exec")
        except SyntaxError as e:
            f.append(f"does not compile: {e}")
        else:
            def run(env, argv):
                os.environ.pop("F1_SLOPES", None)
                if env:
                    os.environ["F1_SLOPES"] = env
                old = sys.argv
                sys.argv = argv
                ns = {}
                try:
                    exec(compile(new, "x", "exec"), ns)
                    return ns, None
                except SystemExit as e:
                    return None, str(e)
                finally:
                    sys.argv = old
                    os.environ.pop("F1_SLOPES", None)

            def labels(ns):
                out = []
                for _t, rows in ns["CONTRASTS"]:
                    for r in rows:
                        if "slope" in r[1]:
                            out.append(r[1].split("slope ")[1])
                return out

            # 1. plain 50 bp run: labels must still read 0.584 / 0.988
            ns, err = run(None, ["43_figure1_design.py", "--out", "F1_design"])
            if err:
                f.append(f"50 bp control refused: {err[:60]}")
            elif labels(ns) != ["0.584", "0.988"]:
                f.append(f"50 bp labels wrong: {labels(ns)}")

            # 2. THE REAL BUG: a 125 bp name with no F1_SLOPES must refuse
            ns, err = run(None, ["43_figure1_design.py", "--out", "F1_design_rl125"])
            if ns is not None or not err or "F1_SLOPES is unset" not in err:
                f.append("a 125 bp name without F1_SLOPES did not refuse "
                         "(this is the run that produced the wrong figure)")

            # 3. the reverse: slopes with no 125 bp name would overwrite the
            #    published figure in place
            ns, err = run("0.690/0.968", ["43_figure1_design.py", "--out", "F1_design"])
            if ns is not None or not err or "no 125 bp output" not in err:
                f.append("F1_SLOPES without a 125 bp name did not refuse")

            # 4. both: the labels must carry the 125 bp values
            ns, err = run("0.690/0.968",
                          ["43_figure1_design.py", "--out", "F1_design_rl125"])
            if err:
                f.append(f"the correct invocation refused: {err[:70]}")
            elif labels(ns) != ["0.690", "0.968"]:
                f.append(f"125 bp labels wrong: {labels(ns)}")

            # 5. three decimals, so 0.690 does not render as 0.69
            if not err and ns and not all(len(x) == 5 for x in labels(ns)):
                f.append("slopes are not rendered to three decimals")

        again, p2 = apply_edits(new)
        if not p2 or "already carries 169" not in p2[0]:
            f.append("re-application did not refuse")

    print("169 self-test")
    for x in f:
        print("  ✗ " + x)
    if not f:
        print("  ✓ the 50 bp labels still read 0.584 / 0.988 with no env set;")
        print("    a 125 bp --out WITHOUT F1_SLOPES refuses (the run that wrote")
        print("    the wrong figure); F1_SLOPES WITHOUT a 125 bp --out refuses")
        print("    too, so the published figure cannot be overwritten in place;")
        print("    with both, the labels read 0.690 / 0.968 to three decimals;")
        print("    re-application refuses")
    return 1 if f else 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("target", nargs="?", default=TARGET)
    ap.add_argument("--apply", action="store_true")
    ap.add_argument("--self-test", action="store_true")
    a = ap.parse_args()
    if a.self_test:
        return self_test()
    if self_test():
        print("\nREFUSING: 169's own self-test failed.")
        return 3
    print()
    if not os.path.exists(a.target):
        print(f"REFUSING: {a.target} not found")
        return 3
    src = open(a.target, errors="replace").read()
    new, probs = apply_edits(src)
    print("=" * 74)
    print("169 — Figure 1's matched-pair slopes")
    print("=" * 74)
    if probs:
        print("REFUSING:")
        for p in probs:
            print("    " + p)
        return 3
    print("".join(difflib.unified_diff(src.splitlines(True), new.splitlines(True),
                                       fromfile=a.target,
                                       tofile=a.target + " (patched)", n=2)))
    if not a.apply:
        print("DRY RUN. Re-run with --apply, then:")
        print()
        print("  # 0. the wrong file must go first")
        print("  rm -f F1_design_rl125.png F1_design_rl125.pdf")
        print()
        print("  # 1. 50 bp control — must reproduce the published figure")
        print("  python3 43_figure1_design.py --out /tmp/f1_cal")
        print()
        print("  # 2. the guard — must REFUSE, nothing written")
        print("  python3 43_figure1_design.py --out F1_design_rl125")
        print()
        print("  # 3. the reverse guard — must also REFUSE")
        print("  F1_SLOPES=0.690/0.968 python3 43_figure1_design.py --out F1_design")
        print()
        print("  # 4. the real one")
        print("  F1_SLOPES=0.690/0.968 python3 43_figure1_design.py \\")
        print("      --out F1_design_rl125")
        print()
        print("  Watch for  [169] matched-pair slopes drawn: platform 0.690,")
        print("  people 0.968  (125 bp)  — then read the two labels on the PNG.")
        return 0
    bak = f"{a.target}.bak.{time.strftime('%Y%m%d_%H%M%S')}"
    shutil.copy2(a.target, bak)
    open(a.target, "w").write(new)
    print(f"[write] {a.target}   (backup: {bak})")
    disk = open(a.target, errors="replace").read()
    bad = verify(disk)
    try:
        compile(disk, a.target, "exec")
    except SyntaxError as e:
        bad.append(f"does not compile: {e}")
    if bad or disk != new:
        print("✗ ON-DISK VERIFICATION FAILED — restoring:")
        for b in bad:
            print(f"    {b}")
        shutil.copy2(bak, a.target)
        return 3
    print("✓ verified: both slopes are now drawn from variables, the 50 bp")
    print("  values remain the defaults, and the guard refuses in both directions.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
