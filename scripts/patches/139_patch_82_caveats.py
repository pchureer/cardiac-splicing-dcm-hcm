#!/usr/bin/env python3
"""
139_patch_82_caveats.py -- the literals patch 138 missed, in the block 82
itself labels "CAVEATS THAT MUST TRAVEL WITH THE TABLE".

WHY 138 MISSED THEM
-------------------
138's survey grepped for `EXPECT_`, slope values and filenames. None of these
four sites is any of those: three say "203" and one says "0.281", in printed
prose. Rule 8 — a rename is only half an edit — for the fourth time this week,
and this time inside the caveats the script insists must accompany the table.

  line  50  docstring  "G1  the three arms share exactly 203 events"
  line  51  docstring  "G2  direction agreement 162/203 ... 174/203"
  line 340  --list     "203 common events, 10 DCM vs 14 NF in both arms."
  line 525  caveat     "* 203 events, both arms at 10 DCM vs 14 NF"
  line 543  caveat     "They disagree by up to 0.281 for Brigham"

At 125 bp the common set is **161**, not 203.

THE 0.281 IS WORSE THAN STALE
-----------------------------
It is the largest disagreement between the Figure 4 transfer slopes and the
genome-wide slopes, and it was computed as Brigham 0.757 − 0.476 = 0.281.
**0.757 is the unsourced value** — it matches nothing in any surviving
symmetry output (see DEFECT_genome_wide_brigham_fuwai_unsourced.md). So 0.281
is a figure derived from an unreproducible number, printed in a caveat that
travels with a supplementary table.

Rather than substitute another literal, this patch DERIVES it:

    _gap  = {c: abs(GENOME_WIDE[c] - PER_COHORT[c]) for c in GENOME_WIDE}
    _worst = max(_gap, key=_gap.get)

which cannot go stale again, and names the cohort rather than assuming
Brigham. On the unsuffixed path it still prints exactly "0.281 for Brigham",
so the 50 bp calibration stays byte-identical. At 125 bp it prints **0.212**;
with the reproducible 50 bp Brigham value (0.662) it would print 0.186.

NOT TOUCHED, DELIBERATELY
-------------------------
  line  28  "ended up contradicting Figure 4 by 0.281" — a historical
            statement about what RESULT_magnitude_resolved.md did, not a
            current claim about the data.
  line 225  "How the 203 is defined was NOT recorded" — describes the
            recorded 50 bp fact and the reason G1 searches for it.

USAGE
-----
  python 139_patch_82_caveats.py [--apply] [--self-test]
"""
import argparse
import difflib
import os
import shutil
import sys
import time

TARGET = "82_build_s4.py"

A_DOC = """  G1  the three arms share exactly 203 events
  G2  direction agreement 162/203 (arm A) and 174/203 (arm B)
"""
B_DOC = """  G1  the three arms share exactly EXPECT_COMMON events
      (203 at 50 bp, 161 at 125 bp)
  G2  direction agreement EXPECT_DIR_A and EXPECT_DIR_B
      (162/203 and 174/203 at 50 bp; 130/161 and 140/161 at 125 bp)
"""

A_LIST = ('        print("              203 common events, '
          '10 DCM vs 14 NF in both arms.")\n')
B_LIST = ('        print(f"              {EXPECT_COMMON} common events, '
          '10 DCM vs 14 NF in both arms.")\n')

A_CAV = ('    print("   * 203 events, both arms at 10 DCM vs 14 NF. '
         'Both slopes are attenuated")\n')
B_CAV = ('    print(f"   * {EXPECT_COMMON} events, both arms at '
         '10 DCM vs 14 NF. Both slopes are attenuated")\n')

A_GAP = ('    print("    script 21\'s genome-wide slopes. '
         'They disagree by up to 0.281 for")\n'
         '    print("    Brigham. Both are correct for what they are; '
         'a sentence must say which.")\n')
B_GAP = ('    # DERIVED, not transcribed. The old literal 0.281 was Brigham\n'
         '    # 0.757 - 0.476, and 0.757 matches nothing in any surviving\n'
         '    # symmetry output — see\n'
         '    # DEFECT_genome_wide_brigham_fuwai_unsourced.md. Computing it\n'
         '    # also stops it going stale and stops it assuming Brigham.\n'
         '    _gap = {_c: abs(GENOME_WIDE[_c] - PER_COHORT[_c])\n'
         '            for _c in GENOME_WIDE if _c in PER_COHORT}\n'
         '    _worst = max(_gap, key=_gap.get)\n'
         '    print(f"    script 21\'s genome-wide slopes. '
         'They disagree by up to {_gap[_worst]:.3f} for")\n'
         '    print(f"    {_worst}. Both are correct for what they are; '
         'a sentence must say which.")\n')

EDITS = [
    ("the docstring's gate list", A_DOC, B_DOC, 1),
    ("the --list common-event count", A_LIST, B_LIST, 1),
    ("the caveat's event count", A_CAV, B_CAV, 1),
    ("the slope-disagreement caveat, derived", A_GAP, B_GAP, 1),
]


def apply_edits(text):
    if "_gap = {_c: abs(GENOME_WIDE[_c]" in text:
        return text, ["the derived gap is already present — 82 already "
                      "carries 139 (restore the .bak to re-apply)"]
    if 'SUF = os.environ.get("TBL_SUFFIX", "")' not in text:
        return text, ["138 has not been applied — apply 138_patch_82_s4.py "
                      "first; these edits reference EXPECT_COMMON as a "
                      "runtime value"]
    problems, out = [], text
    for name, anchor, repl, want in EDITS:
        got = out.count(anchor)
        if got != want:
            problems.append(
                f"{name}: anchor occurs {got}x, expected {want}x\n"
                f"      anchor: {anchor.strip().splitlines()[0][:72]}")
            continue
        out = out.replace(anchor, repl, want)
    return out, problems


def verify(text):
    bad = []
    for w in ("the three arms share exactly EXPECT_COMMON events",
              "(203 at 50 bp, 161 at 125 bp)",
              "direction agreement EXPECT_DIR_A and EXPECT_DIR_B",
              'print(f"              {EXPECT_COMMON} common events',
              'print(f"   * {EXPECT_COMMON} events, both arms',
              "_gap = {_c: abs(GENOME_WIDE[_c] - PER_COHORT[_c])",
              "_worst = max(_gap, key=_gap.get)",
              "{_gap[_worst]:.3f} for",
              'print(f"    {_worst}. Both are correct'):
        if w not in text:
            bad.append(f"missing: {w[:58]}")
    for gone in ("  G1  the three arms share exactly 203 events",
                 "direction agreement 162/203 (arm A)",
                 'print("              203 common events',
                 'print("   * 203 events, both arms',
                 "They disagree by up to 0.281 for"):
        if gone in text:
            bad.append(f"stale site survives: {gone[:50]}")
    # the two sites deliberately left alone must STILL be there
    for keep in ("ended up contradicting Figure 4 by 0.281",
                 "How the 203 is defined was NOT recorded"):
        if keep not in text:
            bad.append(f"a deliberately-kept site was altered: {keep[:46]}")
    # _gap must be computed before it is printed
    i_g = text.find("_gap = {_c: abs(")
    i_u = text.find("{_gap[_worst]:.3f}")
    if i_g < 0 or i_u < 0 or i_g > i_u:
        bad.append("_gap is used before it is computed")
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
    print("139 — 82's printed caveats: three '203's and one unsourced 0.281")
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
        print("  # 1. CALIBRATION — the 50 bp caveats must be UNCHANGED,")
        print("  #    including '203 events' and '0.281 for Brigham'")
        print("  cp S4.tsv S4.PRESERVED2")
        print("  python3 82_build_s4.py 2>&1 | tail -12")
        print("  cmp S4.tsv S4.PRESERVED2 && echo '50bp BYTE-IDENTICAL'")
        print()
        print("  # 2. 125 bp — expect '161 events' and '0.212 for Brigham'")
        print("  TBL_SUFFIX=_rl125 \\")
        print("  S4_RECORD=common=161,A=0.690,B=0.968,A_CI=0.570:0.806,"
              "B_CI=0.875:1.069,dirA=130,dirB=140,product=0.668,C=0.600 \\")
        print("  S4_PER_COHORT=Colorado=0.5999,Brigham=0.5099,"
              "Stanford=0.4008,Fuwai=0.2948 \\")
        print("  S4_GENOME_WIDE=Colorado=0.674,Brigham=0.722,"
              "Fuwai=0.365,Stanford=0.605 \\")
        print("  python3 82_build_s4.py 2>&1 | tail -12")
        print()
        print("⚠ The derived gap reproduces 0.281 at 50 bp ONLY because")
        print("  GENOME_WIDE still holds the unsourced Brigham 0.757 on the")
        print("  unsuffixed path. That keeps calibration honest, but the 50 bp")
        print("  caveat is still quoting a figure built from a number nothing")
        print("  reproduces. With 22's actual 50 bp Brigham (0.662) the gap is")
        print("  0.186. Decide whether the 50 bp table is being re-issued at")
        print("  all before quoting 0.281 anywhere.")
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
    print("✓ verified on disk: every printed '203' follows EXPECT_COMMON, the")
    print("  slope disagreement is computed and names its own cohort, and the")
    print("  two historical mentions are untouched.")
    return 0


def self_test():
    import subprocess
    import tempfile

    fx = ('"""\n' + A_DOC
          + "\nis how `RESULT_magnitude_resolved.md` ended up contradicting "
            "Figure 4 by 0.281. So the\n"
          + '"""\n'
            "import os\n"
            'SUF = os.environ.get("TBL_SUFFIX", "")\n'
            "EXPECT_COMMON = int(os.environ.get('COMMON', '203'))\n"
            "EXPECT_DIR_A = 162\n"
            "EXPECT_DIR_B = 174\n"
            'PER_COHORT = {"Colorado": 0.521, "Brigham": 0.476, '
            '"Stanford": 0.386, "Fuwai": 0.243}\n'
            'GENOME_WIDE = {"Colorado": 0.571, "Brigham": 0.757, '
            '"Fuwai": 0.359, "Stanford": 0.578}\n'
            "# How the 203 is defined was NOT recorded. The record says\n"
            'if os.environ.get("GW125"):\n'
            '    PER_COHORT = {"Colorado": 0.5999, "Brigham": 0.5099,\n'
            '                  "Stanford": 0.4008, "Fuwai": 0.2948}\n'
            '    GENOME_WIDE = {"Colorado": 0.674, "Brigham": 0.722,\n'
            '                   "Fuwai": 0.365, "Stanford": 0.605}\n'
            "def main():\n"
            "    if True:\n"
          + A_LIST
          + A_CAV
          + A_GAP
          + "if __name__ == '__main__':\n"
            "    main()\n")

    new, probs = apply_edits(fx)
    assert not probs, probs
    assert not verify(new), verify(new)
    compile(new, "fx", "exec")

    d = tempfile.mkdtemp()
    p = os.path.join(d, "fx.py")
    open(p, "w").write(new)

    def run(env=None):
        e = dict(os.environ)
        for k in ("TBL_SUFFIX", "COMMON", "GW125"):
            e.pop(k, None)
        e.update(env or {})
        return subprocess.run([sys.executable, "-I", p], capture_output=True,
                              text=True, cwd=d, env=e)

    # 1. unsuffixed -> the 50 bp caveats must be CHARACTER-IDENTICAL,
    #    including "0.281 for" and "Brigham."
    r = run()
    assert r.returncode == 0, r.stdout + r.stderr
    assert "              203 common events, 10 DCM vs 14 NF in both arms." \
        in r.stdout, r.stdout
    assert "   * 203 events, both arms at 10 DCM vs 14 NF. Both slopes are " \
        "attenuated" in r.stdout, r.stdout
    assert "They disagree by up to 0.281 for" in r.stdout, r.stdout
    assert "    Brigham. Both are correct for what they are; a sentence must " \
        "say which." in r.stdout, r.stdout

    # 2. 125 bp values -> 161 everywhere, and the gap recomputed to 0.212
    r2 = run({"COMMON": "161", "GW125": "1"})
    assert r2.returncode == 0, r2.stdout + r2.stderr
    assert "              161 common events" in r2.stdout, r2.stdout
    assert "   * 161 events, both arms" in r2.stdout, r2.stdout
    assert "203" not in r2.stdout, r2.stdout
    assert "They disagree by up to 0.212 for" in r2.stdout, r2.stdout
    assert "    Brigham. Both are correct" in r2.stdout, r2.stdout
    assert "0.281" not in r2.stdout, r2.stdout

    # 3. the worst cohort is NOT hard-coded: make Stanford the worst
    r3 = run({"COMMON": "161", "GW125": "1",
              "TBL_SUFFIX": ""})
    assert "Brigham" in r3.stdout, r3.stdout
    fx2 = new.replace('"Stanford": 0.605}', '"Stanford": 0.905}')
    p2 = os.path.join(d, "fx2.py")
    open(p2, "w").write(fx2)
    e = dict(os.environ, COMMON="161", GW125="1")
    e.pop("TBL_SUFFIX", None)
    r4 = subprocess.run([sys.executable, "-I", p2], capture_output=True,
                        text=True, cwd=d, env=e)
    assert r4.returncode == 0, r4.stdout + r4.stderr
    assert "They disagree by up to 0.504 for" in r4.stdout, r4.stdout
    assert "    Stanford. Both are correct" in r4.stdout, r4.stdout

    # 4. refuses if 138 was never applied
    _n, p5 = apply_edits("print('bare')\n")
    assert p5 and "138 has not been applied" in p5[0], p5

    # 5. double application refused, text unchanged
    n6, p6 = apply_edits(new)
    assert p6, "second application was not refused"
    assert n6 == new, "a refused re-application still altered the text"

    print("[ok] self-test passed: the 50 bp caveats print character-identically "
          "including '203 events' and '0.281 for Brigham'; with 125 bp values "
          "every count becomes 161 with no '203' left and the gap recomputes "
          "to 0.212; the worst cohort is derived, not assumed (a larger "
          "Stanford gap makes it say Stanford 0.504); the two historical "
          "mentions are preserved; refuses without 138; double application "
          "refused")
    return 0


if __name__ == "__main__":
    sys.exit(main())
