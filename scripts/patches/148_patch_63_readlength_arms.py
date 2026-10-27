#!/usr/bin/env python3
"""
148_patch_63_readlength_arms.py -- 63 gains a read-length basis, PER ARM.

WHY PER ARM
-----------
Only MAGNet was mis-declared: rMATS was given `--readLength 50` for reads that
are 125 bp. Colorado, Brigham, Fuwai and Stanford were always quantified at
their true lengths, so they have no `_rl125` twin and their files are
IDENTICAL at both read lengths. `ls -d rmats_*` confirms it — the only
suffixed directories are rmats_magnet_rl125, rmats_HCMvsNF_rl125 and
rmats_DCMvsHCM_rl125.

A uniform suffix would therefore fail for four arms out of five. A uniform
suffix with a silent fallback would be worse: it would produce a table that
looks re-based and mostly is not, which is the hybrid trap 146 had to undo in
62's part (a) earlier today.

So re-basing is per arm, and it is LOUD in two places:

  1. the run prints which arms moved and which did not, with the reason;
  2. the arm NAME carries its basis into the table — "MAGNet (discovery)
     @rl125" against "Colorado @as-published".

(2) matters more than it looks. Four of the five rows of `S5a_rl125.tsv` will
be byte-identical to `S5a.tsv`, because nothing about those cohorts changed.
A reader diffing the two files has to be able to tell "nothing happened" from
"nothing should happen", and only the table itself can tell them.

If TBL_SUFFIX is set and NOTHING re-bases, the script refuses rather than
writing the 50 bp table under a 125 bp name.

WHAT IS NOT CHANGED
-------------------
ANCHOR stays ("probe_rl50", "SE", 172). It is the probe, not an arm, and 172
is a fixed known value at either read length — the gate still has to pass.

USAGE
-----
  python 148_patch_63_readlength_arms.py [--apply] [--self-test]
"""
import argparse
import difflib
import os
import shutil
import subprocess
import sys
import time

TARGET = "63_s5a_jc_vs_jcec.py"

A_SUF = 'ANCHOR = ("probe_rl50", "SE", 172)\n'

B_SUF = A_SUF + '''
# Read length. Set TBL_SUFFIX=_rl125 to use each arm's re-quantified
# directory WHERE ONE EXISTS. Only MAGNet was mis-declared, so only MAGNet
# has one; the rest are identical at both read lengths and stay as published.
# ANCHOR above is deliberately NOT re-based: probe_rl50/SE = 172 is a fixed
# known value and the gate must pass at either setting.
SUF = os.environ.get("TBL_SUFFIX", "")
'''

A_OUT = '    ap.add_argument("--out", default="S5a")\n'
B_OUT = '    ap.add_argument("--out", default="S5a" + SUF)\n'

A_ARMS = '    arms = ARMS + (EXTRA if a.extra else [])\n'

B_ARMS = '''    arms = ARMS + (EXTRA if a.extra else [])

    # ---- READ-LENGTH RE-BASING, PER ARM, added 21 Oct ------------------
    # Only MAGNet was mis-declared (--readLength 50 on 125 bp reads), so
    # only MAGNet has an _rl125 twin. Colorado, Brigham, Fuwai and
    # Stanford were always quantified at their true lengths and their
    # files are IDENTICAL at both read lengths.
    #
    # Four of five rows of S5a_rl125.tsv are therefore byte-identical to
    # S5a.tsv. That is correct, and the table has to SAY it is correct --
    # a reader diffing the two files must be able to tell "nothing
    # happened" from "nothing should happen". So the basis goes into the
    # arm name, not just into this log.
    if SUF:
        _rebased, _kept, _arms2 = [], [], []
        for _n, _d in arms:
            if os.path.isdir(_d + SUF):
                _rebased.append(_n)
                _arms2.append((f"{_n} @{SUF.lstrip('_')}", _d + SUF))
            else:
                _kept.append(_n)
                _arms2.append((f"{_n} @as-published", _d))
        print(f"\\n  TBL_SUFFIX={SUF}: re-based {len(_rebased)} of "
              f"{len(arms)} arms")
        for _n in _rebased:
            print(f"    re-based      {_n}")
        for _n in _kept:
            print(f"    as published  {_n}  (no {SUF} twin; this arm was "
                  f"never mis-declared)")
        if not _rebased:
            print(f"\\n  >>> (a) NOT BUILT. TBL_SUFFIX={SUF} re-based NOTHING.")
            print("      No arm directory carries that suffix, so this would")
            print("      be the 50 bp table written under a 125 bp name.")
            return
        arms = _arms2
'''

EDITS = [
    ("the SUF definition", A_SUF, B_SUF, 1),
    ("the --out default", A_OUT, B_OUT, 1),
    ("the arm list", A_ARMS, B_ARMS, 1),
]


def apply_edits(text):
    if 'SUF = os.environ.get("TBL_SUFFIX", "")' in text:
        return text, ["63 already carries 148 (restore the .bak to re-apply)"]
    problems, out = [], text
    for name, anchor, repl, want in EDITS:
        got = out.count(anchor)
        if got != want:
            problems.append(
                f"{name}: anchor occurs {got}x, expected {want}x\n"
                f"      first line: {anchor.splitlines()[0][:70]}")
            continue
        out = out.replace(anchor, repl, want)
    return out, problems


def verify(text):
    bad = []
    for w in ('SUF = os.environ.get("TBL_SUFFIX", "")',
              'default="S5a" + SUF',
              'if os.path.isdir(_d + SUF):',
              '@as-published',
              're-based NOTHING'):
        if w not in text:
            bad.append(f"missing: {w[:58]}")
    if 'ANCHOR = ("probe_rl50", "SE", 172)' not in text:
        bad.append("ANCHOR was changed — it is the probe, not an arm, and 172 "
                   "must hold at either read length")
    if 'ANCHOR = ("probe_rl50" + SUF' in text or 'ANCHOR[0] + SUF' in text:
        bad.append("ANCHOR was re-based — the gate would then test nothing")
    if 'default="S5a")' in text:
        bad.append("an unsuffixed --out default survives: a 125 bp run would "
                   "overwrite the 50 bp S5a.tsv")
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
    print("148 — 63 gains a per-arm read-length basis")
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
        print("DRY RUN — nothing written. Re-run with --apply, then:")
        print()
        print("  # 50 bp — must reproduce S5a.tsv EXACTLY (no SUF set)")
        print("  cp S5a.tsv S5a_50bp_BEFORE_148.tsv")
        print("  python3 63_s5a_jc_vs_jcec.py 2>&1 | tail -30")
        print("  diff S5a.tsv S5a_50bp_BEFORE_148.tsv && echo IDENTICAL")
        print()
        print("  # 125 bp")
        print("  TBL_SUFFIX=_rl125 python3 63_s5a_jc_vs_jcec.py 2>&1 | tail -35")
        print()
        print("⚠ EXPECTED at 125 bp: `re-based 1 of 5 arms`, MAGNet moving to")
        print("  rmats_magnet_rl125 and the four cohorts reported `as")
        print("  published`. Only the MAGNet row's numbers should differ from")
        print("  S5a.tsv; the other four are the same files and must match.")
        print()
        print("⚠ The anchor (probe_rl50 SE = 172) must pass at BOTH settings.")
        print("  It is not re-based on purpose.")
        print()
        print("⬜ THEN: the JCEC gain was previously said to track declared")
        print("  read length (rho = -0.975). That ordering used MAGNet's")
        print("  declared 50. With MAGNet at 125 the distances become")
        print("  Colorado 74, Brigham 24, Stanford/Fuwai 25 — so recompute")
        print("  before repeating the claim. It has already been decided to")
        print("  retire it; this run is the measurement that records why.")
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
    print("✓ verified on disk: each arm re-bases only if its suffixed")
    print("  directory exists, the basis is carried in the arm name into the")
    print("  table, a suffix that re-bases nothing refuses, --out follows the")
    print("  suffix, and ANCHOR is untouched.")
    return 0


def self_test():
    import tempfile

    fx = ("import argparse, os, sys\n"
          'ARMS = [("MAGNet (discovery)", "rmats_magnet"),\n'
          '        ("Colorado", "rmats_Colorado")]\n'
          'EXTRA = [("MAGNet HCM vs NF", "rmats_HCMvsNF")]\n'
          + A_SUF +
          "\n"
          "def main():\n"
          "    ap = argparse.ArgumentParser()\n"
          '    ap.add_argument("--extra", action="store_true")\n'
          + A_OUT +
          "    a = ap.parse_args()\n"
          + A_ARMS +
          '    print("ARMS", arms)\n'
          '    print("OUT", a.out)\n'
          '    print("ANCHOR", ANCHOR)\n'
          "main()\n")

    new, probs = apply_edits(fx)
    assert not probs, probs
    assert not verify(new), verify(new)
    compile(new, "fx", "exec")

    d = tempfile.mkdtemp()
    for n in ("rmats_magnet", "rmats_Colorado", "rmats_HCMvsNF",
              "rmats_magnet_rl125", "rmats_HCMvsNF_rl125"):
        os.makedirs(os.path.join(d, n), exist_ok=True)

    def run(src, name, env=None, argv=()):
        p = os.path.join(d, name)
        open(p, "w").write(src)
        e = dict(os.environ)
        e.pop("TBL_SUFFIX", None)
        e.update(env or {})
        r = subprocess.run([sys.executable, "-I", p, *argv],
                           capture_output=True, text=True, cwd=d, env=e)
        assert r.returncode == 0, r.stdout + r.stderr
        return r.stdout

    # --- unsuffixed: byte-identical behaviour to the pre-patch script
    o_pre = run(fx, "fx_pre.py")
    o_0 = run(new, "fx_post.py")
    assert "OUT S5a\n" in o_0, o_0
    for line in o_pre.splitlines():
        if line.startswith("ARMS ") or line.startswith("OUT ") \
                or line.startswith("ANCHOR "):
            assert line in o_0, (line, o_0)

    # --- suffixed: MAGNet re-bases, Colorado does not, and BOTH say so
    o_1 = run(new, "fx_post.py", {"TBL_SUFFIX": "_rl125"})
    assert "re-based 1 of 2 arms" in o_1, o_1
    assert "re-based      MAGNet (discovery)" in o_1, o_1
    assert "as published  Colorado" in o_1, o_1
    assert "never mis-declared" in o_1, o_1
    assert "('MAGNet (discovery) @rl125', 'rmats_magnet_rl125')" in o_1, o_1
    assert "('Colorado @as-published', 'rmats_Colorado')" in o_1, o_1
    assert "OUT S5a_rl125\n" in o_1, o_1
    # the anchor is NOT re-based
    assert "ANCHOR ('probe_rl50', 'SE', 172)" in o_1, o_1

    # --- --extra brings in a second re-basable arm
    o_2 = run(new, "fx_post.py", {"TBL_SUFFIX": "_rl125"}, ("--extra",))
    assert "re-based 2 of 3 arms" in o_2, o_2
    assert "('MAGNet HCM vs NF @rl125', 'rmats_HCMvsNF_rl125')" in o_2, o_2

    # --- A SUFFIX THAT RE-BASES NOTHING MUST REFUSE. Without this the
    #     script would write the 50 bp numbers into S5a_rl999.tsv.
    o_3 = run(new, "fx_post.py", {"TBL_SUFFIX": "_rl999"})
    assert "re-based NOTHING" in o_3, o_3
    assert "50 bp table written under a 125 bp name" in o_3, o_3
    assert "ARMS " not in o_3, o_3          # it returned before computing
    assert "OUT " not in o_3, o_3

    # --- double application refused without altering the text
    n4, p4 = apply_edits(new)
    assert p4 and "already carries 148" in p4[0], p4
    assert n4 == new, "a refused re-application still altered the text"

    print("[ok] self-test passed: unsuffixed behaviour is identical to the "
          "pre-patch script (arms, out and anchor lines all reproduced); "
          "TBL_SUFFIX re-bases ONLY the arms that have a suffixed directory "
          "and writes the basis into each arm name (@rl125 / @as-published) "
          "so the four unchanged rows are self-explaining; --extra re-bases "
          "its arm too; a suffix that matches nothing REFUSES instead of "
          "writing 50 bp numbers under a 125 bp name; ANCHOR stays "
          "probe_rl50/172 at both settings; double application refused")
    return 0


if __name__ == "__main__":
    sys.exit(main())
