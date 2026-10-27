#!/usr/bin/env python3
"""
120_patch_36_caveat_readlength.py -- the printed caveat in 36 is factually
wrong at 125 bp, and it is wrong in the direction that understates a confound.

WHAT IS WRONG
-------------
`36_matched_pair_analysis.py` lines 223-224 print:

    * the LBNL arm is single-end 51bp against MAGNet's paired 50bp, so
      'platform' bundles prep, layout and pipeline batch together. This
      design cannot separate them.

At 50 bp that was right, and it was right in an important way: **51 bp against
50 bp means read length was effectively MATCHED**, so 'platform' bundled three
components and read length was not one of them.

At 125 bp MAGNet is paired 125 bp against LBNL's single-end 51 bp. Read length
is now a **fourth bundled component**, and a large one. The A term (0.690) is
the slope across that whole bundle.

This is printed output that gets quoted into the Methods and the Discussion, so
the sentence has to move with the tables. Leaving it says MAGNet is 50 bp in a
paper whose every other number is 125 bp.

It also matters for the interpretation. The caveat is what stops A from being
read as a clean 'platform' effect. Under-describing the bundle makes the
technical attenuation term look better identified than it is — so the honest
version is the stronger one, as the project note already observed.

WHY A PATCH AND NOT A SED
-------------------------
Both lines are plain `print("...")`, not f-strings, so a `{...}` expression
inserted into them would print braces literally. The two lines are replaced
with four: a conditional read length and a conditional bundle list, both keyed
on `SUF`, so the 50 bp wording still reproduces exactly when the suffix is
unset.

    TBL_SUFFIX unset   -> "MAGNet's paired 50bp"  / "prep, layout and
                           pipeline batch"
    TBL_SUFFIX=_rl125  -> "MAGNet's paired 125bp" / "prep, layout, READ LENGTH
                           (125 vs 51 bp) and pipeline batch"

Line 225 ("design cannot separate them") is left untouched.

USAGE
-----
  python 120_patch_36_caveat_readlength.py
  python 120_patch_36_caveat_readlength.py --apply
  python 120_patch_36_caveat_readlength.py --self-test
"""
import argparse
import difflib
import os
import shutil
import sys
import time

TARGET = "36_matched_pair_analysis.py"

OLD = (
    '        print("   * the LBNL arm is single-end 51bp against MAGNet\'s '
    'paired 50bp, so")\n'
    '        print("     \'platform\' bundles prep, layout and pipeline batch '
    'together. This")\n'
)

NEW = (
    '        # Read length is a bundled component ONLY when MAGNet is 125 bp.\n'
    '        # At 50 bp against LBNL\'s 51 bp it was effectively matched, so the\n'
    '        # unsuffixed wording is not a simplification -- it was correct.\n'
    '        _rl = "125bp" if SUF else "50bp"\n'
    '        _bundle = ("prep, layout, READ LENGTH (125 vs 51 bp) and pipeline "\n'
    '                   "batch") if SUF else "prep, layout and pipeline batch"\n'
    '        print(f"   * the LBNL arm is single-end 51bp against MAGNet\'s '
    'paired {_rl}, so")\n'
    '        print(f"     \'platform\' bundles {_bundle} together. This")\n'
)

EDITS = [("the read-length caveat, made conditional on the suffix", OLD, NEW, 1)]


def apply_edits(text):
    problems, out = [], text
    for name, anchor, repl, want in EDITS:
        got = out.count(anchor)
        if got != want:
            problems.append(f"{name}: anchor occurs {got}x, expected {want}x")
            continue
        if "_bundle = (" in out:
            problems.append("already patched (_bundle is defined)")
            continue
        out = out.replace(anchor, repl, want)
    return out, problems


def verify(text):
    bad = []
    code = "\n".join(l for l in text.splitlines()
                     if not l.lstrip().startswith("#"))
    if "_rl = \"125bp\" if SUF else \"50bp\"" not in code:
        bad.append("the conditional read length is missing")
    if "READ LENGTH (125 vs 51 bp)" not in code:
        bad.append("the fourth bundled component is not named")
    if "MAGNet's paired 50bp, so\")" in code:
        bad.append("the unconditional 50bp wording survives")
    if "bundles prep, layout and pipeline batch together. This\")" in code:
        bad.append("the unconditional three-component bundle survives")
    i_rl = text.find('_rl = "125bp"')
    i_use = text.find("MAGNet's paired {_rl}")
    if i_rl < 0 or i_use < 0 or i_rl > i_use:
        bad.append("_rl is not defined above its consumer")
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
    print("120 — the read-length caveat in 36")
    print("=" * 78)
    if not os.path.exists(a.target):
        print(f"[fatal] {a.target} not found.", file=sys.stderr)
        return 2
    src = open(a.target, errors="replace").read()
    if "SUF = os.environ.get" not in src:
        print("[fatal] no SUF in this file; the conditional depends on it. "
              "Apply 119 first.", file=sys.stderr)
        return 2
    new, problems = apply_edits(src)
    if problems:
        print("REFUSING — anchors do not match:")
        for p in problems:
            print(f"  ✗ {p}")
        print("\n  Nothing written. Send me lines 220-228 and I will re-anchor.")
        return 3
    bad = verify(new)
    if bad:
        print("REFUSING — patched text fails verification:")
        for b in bad:
            print(f"  ✗ {b}")
        return 3

    print("".join(difflib.unified_diff(
        src.splitlines(True), new.splitlines(True),
        fromfile=a.target, tofile=a.target + " (patched)", n=2)))
    print("=" * 78)
    if not a.apply:
        print("DRY RUN — nothing written. Re-run with --apply, then re-run 36:")
        print()
        print("  C_OBSERVED=0.600 C_CONCORDANCE=80.8% \\")
        print("  TYPES=SE TBL_SUFFIX=_rl125 python3 36_matched_pair_analysis.py")
        print()
        print("  The caveat should then read 'paired 125bp' and name READ")
        print("  LENGTH as a fourth bundled component. Nothing else moves —")
        print("  this edit changes printed text only, no computation.")
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
        return 3
    print("✓ verified on disk.")
    return 0


def self_test():
    import subprocess
    import tempfile

    runnable = (
        "import os\n"
        'SUF = os.environ.get("TBL_SUFFIX", "")\n'
        "def caveats():\n"
        "    if True:\n"
        + OLD +
        '        print("     design cannot separate them.")\n')
    new, probs = apply_edits(runnable)
    assert not probs, probs
    assert not verify(new), verify(new)

    d = tempfile.mkdtemp()
    p = os.path.join(d, "fx.py")
    open(p, "w").write(new + "\nif __name__ == '__main__':\n    caveats()\n")

    r = subprocess.run([sys.executable, "-I", p], capture_output=True,
                       text=True, cwd=d)
    assert r.returncode == 0, r.stdout + r.stderr
    assert "paired 50bp" in r.stdout, r.stdout
    assert "prep, layout and pipeline batch" in r.stdout, r.stdout
    assert "READ LENGTH" not in r.stdout, r.stdout

    env = dict(os.environ, TBL_SUFFIX="_rl125")
    r2 = subprocess.run([sys.executable, "-I", p], capture_output=True,
                        text=True, cwd=d, env=env)
    assert r2.returncode == 0, r2.stdout + r2.stderr
    assert "paired 125bp" in r2.stdout, r2.stdout
    assert "READ LENGTH (125 vs 51 bp)" in r2.stdout, r2.stdout
    assert "paired 50bp" not in r2.stdout, r2.stdout
    # the untouched third line must still print in both
    assert "design cannot separate them." in r.stdout
    assert "design cannot separate them." in r2.stdout

    _n, p3 = apply_edits(new)
    assert p3, "second application was not refused"

    print("[ok] self-test passed: unsuffixed prints 'paired 50bp' and the "
          "three-component bundle with no READ LENGTH mention; suffixed prints "
          "'paired 125bp' and names READ LENGTH (125 vs 51 bp) as a fourth; the "
          "third caveat line is untouched in both; double application refused")
    return 0


if __name__ == "__main__":
    sys.exit(main())
