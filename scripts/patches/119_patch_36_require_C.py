#!/usr/bin/env python3
"""
119_patch_36_require_C.py -- stop 36 from printing a 125 bp A x B against a
hard-coded 50 bp C.

WHAT IS WRONG
-------------
`36_matched_pair_analysis.py` carries C as a literal in three places:

    line   8  docstring:  C) MAGNet_full -> Colorado  both  = 0.521 (ref)
    line 142  f"{'':>6}{0.521:>9.3f}{'':>22}{'79.6%':>16}")
    line 157  print(f"  A x B = {A_ * B_:.3f}   against the observed
              cross-cohort C = 0.521")

C is **not computed here.** It is the MAGNet -> Colorado Deming slope from
`33_rebuild_summaries_dedup.py` section E. Line 142 additionally hard-codes
`'79.6%'`, the 50 bp Colorado direction concordance.

After the `$TBL_SUFFIX` fix, 36 computed A and B from the 125 bp tables and
printed them against the 50 bp C, with no error:

    A x B = 0.668   against the observed cross-cohort C = 0.521   ->  28% gap

That 28% is an artifact of a cross-read-length comparison -- the exact error
class the 125 bp promotion exists to remove. With C recomputed at 125 bp
(0.600) the gap is 11%, the same as at 50 bp, and the decomposition closes.
See `RESULT_33_transfer_slopes_rl125.md`.

WHAT THIS PATCH DOES
--------------------
C cannot be derived inside 36, so it cannot be gated the way 117 gated 95's
null or 118 gated 96's KS values. Instead the patch makes the wrong comparison
**impossible to make by accident**:

    if TBL_SUFFIX is set, C_OBSERVED and C_CONCORDANCE MUST be supplied.
    There is no default. The script refuses and prints the command to run.

Unsuffixed runs keep the 50 bp values (0.521 / 79.6%) as defaults, so
reproducing the original output needs no environment at all -- the
`envify()` convention used throughout this project.

    TBL_SUFFIX unset                 -> C = 0.521, concordance 79.6%
    TBL_SUFFIX=_rl125, no C_OBSERVED -> REFUSES, non-zero exit
    TBL_SUFFIX=_rl125 C_OBSERVED=0.600 C_CONCORDANCE=80.8%  -> uses those

NOT DONE HERE
-------------
The printed caveat is still factually wrong:

    "the LBNL arm is single-end 51bp against MAGNet's paired 50bp, so
     'platform' bundles prep, layout and pipeline batch together"

MAGNet is 125 bp now, and 'platform' bundles a fourth component: read length.
That text is left alone deliberately -- it may sit inside a plain (non-f)
string, where inserting a `{...}` expression would print braces literally.
Locate and edit it by hand:

    grep -n "paired 50bp\\|bundles prep" 36_matched_pair_analysis.py

CONVENTIONS
-----------
  * dry run by default; --apply writes after a timestamped .bak
  * anchors asserted to occur EXACTLY once
  * re-read FROM DISK and re-verified; restore on failure
  * non-zero exit on any refusal

USAGE
-----
  python 119_patch_36_require_C.py
  python 119_patch_36_require_C.py --apply
  python 119_patch_36_require_C.py --self-test
"""
import argparse
import difflib
import os
import shutil
import sys
import time

TARGET = "36_matched_pair_analysis.py"

ANCHOR_SUF = ('SUF = os.environ.get("TBL_SUFFIX", "")'
              '   # "" = 50 bp, "_rl125" = 125 bp\n')

C_BLOCK = '''
# ---------------------------------------------------------------------------
# C -- the MAGNet -> Colorado cross-cohort slope -- is NOT computed in this
# file. It comes from `33_rebuild_summaries_dedup.py` section E.
#
# Until this block existed, C was the literal 0.521 in three places here, and
# the Colorado direction concordance was the literal '79.6%'. After the
# $TBL_SUFFIX fix this script computed A and B from the 125 bp tables and
# printed them against the 50 bp C:
#
#     A x B = 0.668   against the observed cross-cohort C = 0.521   -> 28% gap
#
# That gap was a cross-read-length comparison, not a result. With C recomputed
# at 125 bp (0.600) the gap is 11% -- the same as at 50 bp -- and the
# decomposition closes.
#
# C cannot be verified from this script's own inputs, so it cannot be gated.
# The rule instead: IF THE TABLES ARE SUFFIXED, C MUST BE SUPPLIED.
C_OBS = os.environ.get("C_OBSERVED", "").strip()
C_CONCORD = os.environ.get("C_CONCORDANCE", "").strip()
if SUF:
    if not C_OBS or not C_CONCORD:
        raise SystemExit(
            "\\nREFUSING — TBL_SUFFIX is set to '" + SUF + "', so A and B will "
            "be computed from\\n  re-quantified tables. C cannot be computed "
            "here; it comes from script 33's\\n  section E, and the 50 bp "
            "default (0.521) would make `A x B against C` a\\n  "
            "cross-read-length comparison. That is what produced the spurious "
            "28% gap.\\n\\n"
            "  Re-run 33 against the same tables:\\n\\n"
            "    python3 33_rebuild_summaries_dedup.py \\\\\\n"
            "        --discovery '{ET}_v2filtered" + SUF + "' \\\\\\n"
            "        --out summaries_dedup" + SUF + ".txt\\n\\n"
            "  then pass its Colorado row from sections E and D:\\n\\n"
            "    C_OBSERVED=<slope> C_CONCORDANCE=<pct> TBL_SUFFIX=" + SUF
            + " \\\\\\n    python3 36_matched_pair_analysis.py\\n\\n"
            "  At 125 bp those are 0.600 and 80.8%.\\n")
else:
    C_OBS = C_OBS or "0.521"
    C_CONCORD = C_CONCORD or "79.6%"
C_OBS_F = float(C_OBS)
# ---------------------------------------------------------------------------
'''

OLD_ROW = "              f\"{'':>6}{0.521:>9.3f}{'':>22}{'79.6%':>16}\")\n"
NEW_ROW = "              f\"{'':>6}{C_OBS_F:>9.3f}{'':>22}{C_CONCORD:>16}\")\n"

OLD_AXB = ('            print(f"  A x B = {A_ * B_:.3f}   against the observed '
           'cross-cohort C = 0.521")\n')
NEW_AXB = ('            print(f"  A x B = {A_ * B_:.3f}   against the observed '
           'cross-cohort C = {C_OBS_F:.3f}")\n')

# docstring: matched as a substring so no whitespace needs reproducing
OLD_DOC = "= 0.521 (ref)"
NEW_DOC = "= C_OBSERVED (ref; script 33 section E — 0.521 at 50 bp, 0.600 at 125 bp)"

EDITS = [
    ("C_OBSERVED / C_CONCORDANCE, required when the tables are suffixed",
     ANCHOR_SUF, ANCHOR_SUF + C_BLOCK, 1),
    ("the C reference row (0.521 and '79.6%' were both literals)",
     OLD_ROW, NEW_ROW, 1),
    ("the A x B vs C closure line", OLD_AXB, NEW_AXB, 1),
    ("the docstring's C reference", OLD_DOC, NEW_DOC, 1),
]


def apply_edits(text):
    problems, out = [], text
    for name, anchor, repl, want in EDITS:
        got = out.count(anchor)
        if got != want:
            problems.append(
                f"{name}: anchor occurs {got}x, expected {want}x\n"
                f"      anchor: {anchor.strip().splitlines()[0][:86]}")
            continue
        if "C_OBS = os.environ.get" in out and name.startswith("C_OBSERVED"):
            problems.append("C_OBS is already defined — already patched")
            continue
        out = out.replace(anchor, repl, want)
    return out, problems


def verify(text):
    bad = []
    code = "\n".join(l for l in text.splitlines()
                     if not l.lstrip().startswith("#"))
    if "C_OBS = os.environ.get" not in code:
        bad.append("C_OBS is not defined")
    if "C_OBS_F = float(C_OBS)" not in code:
        bad.append("C_OBS_F is not derived")
    if "{0.521:>9.3f}" in code:
        bad.append("the C reference row still prints the literal 0.521")
    if "'79.6%'" in code:
        bad.append("the 50 bp concordance literal '79.6%' survives in code")
    if "cross-cohort C = 0.521" in code:
        bad.append("the closure line still prints the literal 0.521")
    if "raise SystemExit" not in code:
        bad.append("the refusal is missing")
    # C_OBS_F must be defined above both consumers
    i_def = text.find("C_OBS_F = float(C_OBS)")
    for marker in ("{C_OBS_F:>9.3f}", "C = {C_OBS_F:.3f}"):
        j = text.find(marker)
        if j < 0:
            bad.append(f"consumer missing: {marker}")
        elif i_def < 0 or i_def > j:
            bad.append(f"C_OBS_F is defined below its consumer {marker}")
    return bad


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("target", nargs="?", default=TARGET)
    ap.add_argument("--apply", action="store_true")
    ap.add_argument("--self-test", action="store_true")
    a = ap.parse_args()
    if a.self_test:
        return self_test()

    self_test()
    print()
    print("=" * 78)
    print("119 — require C when the tables are suffixed")
    print("=" * 78)
    if not os.path.exists(a.target):
        print(f"[fatal] {a.target} not found.", file=sys.stderr)
        return 2
    src = open(a.target, errors="replace").read()
    new, problems = apply_edits(src)
    if problems:
        print("REFUSING — anchors do not match:")
        for p in problems:
            print(f"  ✗ {p}")
        print("\n  Nothing written. Send me the lines around 'SUF =', 0.521 and")
        print("  '79.6%' and I will re-anchor.")
        return 3
    bad = verify(new)
    if bad:
        print("REFUSING — patched text fails its own verification:")
        for b in bad:
            print(f"  ✗ {b}")
        return 3

    print("".join(difflib.unified_diff(
        src.splitlines(True), new.splitlines(True),
        fromfile=a.target, tofile=a.target + " (patched)", n=2)))
    print("=" * 78)

    if not a.apply:
        print("DRY RUN — nothing written. Re-run with --apply.")
        print()
        print("Then the 125 bp run becomes:")
        print()
        print("  C_OBSERVED=0.600 C_CONCORDANCE=80.8% \\")
        print("  TYPES=SE TBL_SUFFIX=_rl125 python3 36_matched_pair_analysis.py")
        print()
        print("  (confirm TYPES against the sbatch — line 180 passes it)")
        print()
        print("A run with TBL_SUFFIX set and no C_OBSERVED now REFUSES instead")
        print("of printing A x B against the 50 bp C. Check that first:")
        print()
        print("  TBL_SUFFIX=_rl125 python3 36_matched_pair_analysis.py")
        print()
        print("⚠ Still to fix BY HAND in the same file — the printed caveat")
        print("  says MAGNet is paired 50bp and that 'platform' bundles three")
        print("  components. It is 125 bp, and read length is a fourth:")
        print()
        print('  grep -n "paired 50bp\\|bundles prep" 36_matched_pair_analysis.py')
        return 0

    bak = f"{a.target}.bak.{time.strftime('%Y%m%d_%H%M%S')}"
    shutil.copy2(a.target, bak)
    open(a.target, "w").write(new)
    print(f"[write] {a.target}   (backup: {bak})")
    disk = open(a.target, errors="replace").read()
    bad = verify(disk)
    if bad or disk != new:
        print("✗ ON-DISK VERIFICATION FAILED — restoring:")
        for b in bad:
            print(f"    {b}")
        shutil.copy2(bak, a.target)
        print(f"  restored from {bak}")
        return 3
    print("✓ verified on disk: no 0.521 or '79.6%' literal remains in code,")
    print("  C_OBS_F is defined above both consumers, and a suffixed run")
    print("  without C_OBSERVED now refuses.")
    return 0


def self_test():
    import re
    import subprocess
    import tempfile

    # The fixture carries each anchor VERBATIM, at the real file's indentation:
    # the reference row is a continuation line of a print( at 14 spaces, and
    # the closure line sits at 12 spaces inside a nested block. Reproducing the
    # indentation is the point -- an anchor is whitespace-exact, and the first
    # version of this self-test reshaped the lines and so tested nothing.
    runnable = (
        '"""\n'
        '    C) MAGNet_full    -> Colorado    both    = 0.521 (ref)\n'
        '"""\n'
        "import os\n"
        'SUF = os.environ.get("TBL_SUFFIX", "")   # "" = 50 bp, "_rl125" = 125 bp\n'
        "\n"
        "def show(A_, B_):\n"
        "    if True:\n"
        "        print(f\"  {'C  reference   MAGNet_full -> Colorado "
        "(script 33)':<52}\"\n"
        "              f\"{'':>6}{0.521:>9.3f}{'':>22}{'79.6%':>16}\")\n"
        "        if True:\n"
        '            print(f"  A x B = {A_ * B_:.3f}   against the observed '
        'cross-cohort C = 0.521")\n')

    new, probs = apply_edits(runnable)
    assert not probs, probs
    assert not verify(new), verify(new)

    d = tempfile.mkdtemp()
    p = os.path.join(d, "fx.py")
    open(p, "w").write(new + "\nif __name__ == '__main__':\n"
                             "    show(0.690, 0.968)\n")

    # 1. unsuffixed -> the 50 bp defaults, output unchanged
    r = subprocess.run([sys.executable, "-I", p], capture_output=True,
                       text=True, cwd=d)
    assert r.returncode == 0, r.stdout + r.stderr
    assert "0.521" in r.stdout and "79.6%" in r.stdout, r.stdout

    # 2. suffixed with no C -> refuse, non-zero
    env = dict(os.environ, TBL_SUFFIX="_rl125")
    r2 = subprocess.run([sys.executable, "-I", p], capture_output=True,
                        text=True, cwd=d, env=env)
    assert r2.returncode != 0, r2.stdout + r2.stderr
    both = r2.stdout + r2.stderr
    assert "REFUSING" in both and "script 33" in both, both
    assert "28% gap" in both, both

    # 3. suffixed with C -> uses it, and 0.521 is gone
    env3 = dict(os.environ, TBL_SUFFIX="_rl125", C_OBSERVED="0.600",
                C_CONCORDANCE="80.8%")
    r3 = subprocess.run([sys.executable, "-I", p], capture_output=True,
                        text=True, cwd=d, env=env3)
    assert r3.returncode == 0, r3.stdout + r3.stderr
    assert "0.600" in r3.stdout and "80.8%" in r3.stdout, r3.stdout
    assert "0.521" not in r3.stdout and "79.6%" not in r3.stdout, r3.stdout
    # and A x B prints against the supplied C
    m = re.search(r"A x B = ([0-9.]+)\s+against the observed cross-cohort "
                  r"C = ([0-9.]+)", r3.stdout)
    assert m, r3.stdout
    assert m.group(1) == "0.668" and m.group(2) == "0.600", m.groups()

    # 4. double application refused
    _n4, p4 = apply_edits(new)
    assert p4, "second application was not refused"

    print("[ok] self-test passed: unsuffixed keeps 0.521/79.6%; suffixed "
          "without C_OBSERVED refuses non-zero and names script 33; suffixed "
          "with C prints A x B = 0.668 against C = 0.600 and no 50 bp literal "
          "survives; double application refused")
    return 0


if __name__ == "__main__":
    sys.exit(main())
