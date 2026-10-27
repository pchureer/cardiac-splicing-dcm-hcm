#!/usr/bin/env python3
"""
137_patch_22_symmetry.py -- 22 at 125 bp. THREE suffix behaviours in one list.

WHY 22 CANNOT BE PATCHED BY PATTERN
-----------------------------------
`22_symmetry_winners_curse.py`'s PAIRS list mixes three different rules, and
applying one rule to all of them would be wrong in two directions at once:

    row                x side                        y side
    MAGNet internal    splitctl/{ET}_DCMvsA_1   ->   splitctl/{ET}_HCMvsB_1
                       BOTH SIDES MOVE — both are MAGNet
    Colorado           {ET}_v2filtered          ->   {ET}_Colorado_DCMvsNF
    LBNL               {ET}_v2filtered          ->   {ET}_LBNL_DCMvsNF
    Brigham            {ET}_v2filtered          ->   Brigham_{ET}_HCM
    Fuwai              {ET}_v2filtered          ->   Fuwai_{ET}_HCM
    Stanford           {ET}_v2filtered          ->   Stanford_{ET}_HCM
                       X SIDE ONLY — no replication cohort was re-quantified

Suffix everything and the five cohort rows read files that do not exist.
Suffix nothing but the discovery stem and the MAGNet-internal row silently
compares 125 bp against 50 bp — which is the error that produced the false 28%
gap in the matched-pair closure.

THE OUTPUT DEFECT, SAME AS 33's
-------------------------------
    out = a.out or f"symmetry_{ET}" + ("_program" if a.program else "")

has no read-length component, and all three outputs derive from it:
`{out}.tsv`, `{out}_sweep.tsv`, `{out}.report.txt`. An unsuffixed 125 bp run
would overwrite all three in place and exit 0 — exactly how
`direction_concordance_dedup.tsv` was lost on 16 Oct. `out` now carries SUF.

PREREQUISITE — 136 MUST HAVE RUN
--------------------------------
The MAGNet-internal row needs `splitctl/{ET}_DCMvsA_1{SUF}.all_events.tsv` and
`splitctl/{ET}_HCMvsB_1{SUF}.all_events.tsv`, which `136_splitctl_rl125.sbatch`
produces. For SE those exist (job 408099, 17 Oct: both 50 bp arms reproduced
byte-identically first, and both read lengths report 385,826 / 332 samples
n=166 and 320,990 / 194 samples n=28 — same partition, same groups). For any
other event type they do not, so this patch adds a guard that refuses and
names 136 rather than letting `load()` skip the row and leave the headline
contrast quietly one-sided.

THE DOCSTRING CLAIM THAT MUST NOT SURVIVE UNQUALIFIED
------------------------------------------------------
Lines 9-10 say disease stage, "read length", control-group size and shared
controls are "all excluded (scripts 16/17/19/20/21)". That exclusion was about
the REPLICATION COHORTS' read lengths — Colorado ~50 bp against Fuwai 150 bp.
It was never about MAGNet's OWN quantification read length, which none of
those scripts varied, and which the promotion has shown moves every cohort
slope by +15% to +29%. Left as it stands, 22 reads as having ruled out the
effect this entire re-baselining demonstrates. The patch names which
hypothesis was tested instead of deleting the sentence.

The five slopes in line 8 are likewise 50 bp values that THIS SCRIPT
RECOMPUTES, yet they are transcribed into `82_build_s4.py` (GENOME_WIDE),
`62_build_s5.py` (S5's warning) and `98_build_replication_table.py:393`. The
docstring now says not to quote them from there.

USAGE
-----
  python 137_patch_22_symmetry.py [--apply] [--self-test]
"""
import argparse
import difflib
import os
import shutil
import sys
import time

TARGET = "22_symmetry_winners_curse.py"

# ------------------------------------------------------------------- edit 1
A_SUF = "def main():\n"
B_SUF = '''# --------------------------------------------------------------------------
# READ LENGTH. Three different rules apply inside PAIRS below — see the
# comments there — and `out` carries SUF because all three of this script's
# outputs derive from it ({out}.tsv, {out}_sweep.tsv, {out}.report.txt) and an
# unsuffixed 125 bp run would overwrite the 50 bp ones in place.
SUF = os.environ.get("TBL_SUFFIX", "")


def _require_internal(pairs):
    """The MAGNet-internal row needs 136's output. Refuse rather than skip.

    `load()` returns None for a missing file and the row drops out silently,
    which would leave 22's headline contrast — cohort slopes against MAGNet
    compared with itself — with nothing on one side of it.
    """
    if not SUF:
        return
    label, x, y = pairs[0]
    missing = [p for p in (x, y) if not os.path.exists(p)]
    if missing:
        raise SystemExit(
            "REFUSING — TBL_SUFFIX=%s but the MAGNet-internal arm is absent:\\n"
            "%s\\n\\n"
            "  Build it with 136_splitctl_rl125.sbatch, which re-runs phase 2\\n"
            "  against rmats_magnet%s/ and rmats_HCMvsNF%s/ on the SAME 83/83\\n"
            "  NF partition (splitctl/nfA_1.ids, nfB_1.ids — reused verbatim,\\n"
            "  never re-drawn):\\n\\n"
            "    ET=%s sbatch 136_splitctl_rl125.sbatch\\n\\n"
            "  Without it this script's headline contrast has no MAGNet side,\\n"
            "  and running it anyway would compare 125 bp cohort slopes\\n"
            "  against a 50 bp MAGNet-internal slope."
            % (SUF, "".join("    missing: %s\\n" % m for m in missing),
               SUF, SUF, os.environ.get("ET_HINT", "SE")))
# --------------------------------------------------------------------------


''' + A_SUF

# ------------------------------------------------------------------- edit 2
A_OUT = ('    out = a.out or f"symmetry_{ET}" '
         '+ ("_program" if a.program else "")\n')
B_OUT = ('    # SUF LAST, after the _program variant, so the stem reads\n'
         '    #   symmetry_SE_program_rl125\n'
         '    out = a.out or f"symmetry_{ET}" '
         '+ ("_program" if a.program else "") + SUF\n')

# ------------------------------------------------------------------- edit 3
A_PAIRS = '''    PAIRS = [
        ("MAGNet internal", f"splitctl/{ET}_DCMvsA_1.all_events.tsv",
                            f"splitctl/{ET}_HCMvsB_1.all_events.tsv"),
        ("Colorado",  f"{ET}_v2filtered.all_events.tsv", f"{ET}_Colorado_DCMvsNF.all_events.tsv"),
        ("LBNL",      f"{ET}_v2filtered.all_events.tsv", f"{ET}_LBNL_DCMvsNF.all_events.tsv"),
        ("Brigham",   f"{ET}_v2filtered.all_events.tsv", f"Brigham_{ET}_HCM.all_events.tsv"),
        ("Fuwai",     f"{ET}_v2filtered.all_events.tsv", f"Fuwai_{ET}_HCM.all_events.tsv"),
        ("Stanford",  f"{ET}_v2filtered.all_events.tsv", f"Stanford_{ET}_HCM.all_events.tsv"),
    ]
'''
B_PAIRS = '''    # THREE RULES, NOT ONE.
    #
    # Row 1, MAGNet internal: BOTH sides are MAGNet, so both carry SUF. Built
    #   by 136_splitctl_rl125.sbatch on the same 83/83 NF partition.
    # Rows 2-6, the cohorts: the x side is the MAGNet discovery table and
    #   carries SUF; the y side is the replication cohort's own quantification
    #   and must NOT — no replication cohort was re-quantified. Suffixing it
    #   would read a file that does not exist; not suffixing the x side would
    #   compare 50 bp discovery against a 125 bp claim.
    PAIRS = [
        ("MAGNet internal", f"splitctl/{ET}_DCMvsA_1{SUF}.all_events.tsv",
                            f"splitctl/{ET}_HCMvsB_1{SUF}.all_events.tsv"),
        ("Colorado",  f"{ET}_v2filtered{SUF}.all_events.tsv", f"{ET}_Colorado_DCMvsNF.all_events.tsv"),
        ("LBNL",      f"{ET}_v2filtered{SUF}.all_events.tsv", f"{ET}_LBNL_DCMvsNF.all_events.tsv"),
        ("Brigham",   f"{ET}_v2filtered{SUF}.all_events.tsv", f"Brigham_{ET}_HCM.all_events.tsv"),
        ("Fuwai",     f"{ET}_v2filtered{SUF}.all_events.tsv", f"Fuwai_{ET}_HCM.all_events.tsv"),
        ("Stanford",  f"{ET}_v2filtered{SUF}.all_events.tsv", f"Stanford_{ET}_HCM.all_events.tsv"),
    ]
    os.environ.setdefault("ET_HINT", ET)
    _require_internal(PAIRS)
'''

# ------------------------------------------------------------------- edit 4
A_DOC = """Every independent cohort recapitulates the MAGNet signature at a Deming slope near
0.56 (Colorado 0.571, LBNL 0.525, Brigham 0.757, Fuwai 0.359, Stanford 0.578), while
MAGNet compared with itself gives 0.92. Disease stage, read length, control-group size
and shared controls are all excluded (scripts 16/17/19/20/21). The remaining candidate
is that the comparison itself is biased:
"""
B_DOC = """Every independent cohort recapitulates the MAGNet signature at a Deming slope near
0.56 AT 50 bp (Colorado 0.571, LBNL 0.525, Brigham 0.757, Fuwai 0.359, Stanford
0.578), while MAGNet compared with itself gives 0.92.

⚠ THOSE SIX FIGURES ARE 50 bp VALUES AND THIS SCRIPT RECOMPUTES THEM. Do not
quote them from this docstring. They are nevertheless transcribed into
82_build_s4.py (GENOME_WIDE), 62_build_s5.py and 98_build_replication_table.py
line 393; at 125 bp every one of them moves and all four copies need the
run's own output.

Disease stage, COHORT read length, control-group size and shared controls are all
excluded (scripts 16/17/19/20/21).

⚠ "read length" there means the REPLICATION COHORTS' read lengths — Colorado
~50 bp against Fuwai 150 bp and so on. It does NOT mean MAGNet's own
quantification read length, which none of those scripts varied and which the
125 bp promotion has since shown moves every cohort slope by +15% to +29%. The
exclusion stands for the axis it tested and says nothing about the axis it
could not test.

The remaining candidate is that the comparison itself is biased:
"""

EDITS = [
    ("SUF and the MAGNet-internal guard", A_SUF, B_SUF, 1),
    ("the output stem, suffixed", A_OUT, B_OUT, 1),
    ("PAIRS — three rules", A_PAIRS, B_PAIRS, 1),
    ("the docstring's 50 bp values and the read-length exclusion",
     A_DOC, B_DOC, 1),
]


def apply_edits(text):
    if 'SUF = os.environ.get("TBL_SUFFIX", "")' in text:
        return text, ["SUF is already defined — 22 is already patched "
                      "(restore the .bak to re-apply)"]
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
    for w in ('SUF = os.environ.get("TBL_SUFFIX", "")',
              "def _require_internal(pairs):",
              "_require_internal(PAIRS)",
              'f"splitctl/{ET}_DCMvsA_1{SUF}.all_events.tsv"',
              'f"splitctl/{ET}_HCMvsB_1{SUF}.all_events.tsv"',
              'f"{ET}_v2filtered{SUF}.all_events.tsv"',
              '("_program" if a.program else "") + SUF',
              "COHORT read length",
              "THOSE SIX FIGURES ARE 50 bp VALUES"):
        if w not in text:
            bad.append(f"missing: {w[:58]}")
    # the replication stems must NOT have moved
    for y in ('f"{ET}_Colorado_DCMvsNF.all_events.tsv"',
              'f"{ET}_LBNL_DCMvsNF.all_events.tsv"',
              'f"Brigham_{ET}_HCM.all_events.tsv"',
              'f"Fuwai_{ET}_HCM.all_events.tsv"',
              'f"Stanford_{ET}_HCM.all_events.tsv"'):
        if y not in text:
            bad.append(f"a replication stem was altered: {y[:46]}")
    for suffixed in ('_Colorado_DCMvsNF{SUF}', '_LBNL_DCMvsNF{SUF}',
                     'Brigham_{ET}_HCM{SUF}', 'Fuwai_{ET}_HCM{SUF}',
                     'Stanford_{ET}_HCM{SUF}'):
        if suffixed in text:
            bad.append(f"a replication cohort was WRONGLY suffixed: {suffixed}")
    # the unsuffixed discovery / internal stems must be gone
    for gone in ('f"{ET}_v2filtered.all_events.tsv"',
                 'f"splitctl/{ET}_DCMvsA_1.all_events.tsv"',
                 'f"splitctl/{ET}_HCMvsB_1.all_events.tsv"',
                 '("_program" if a.program else "")\n'):
        if gone in text:
            bad.append(f"unsuffixed site survives: {gone[:50]}")
    # SUF module-level and above every use
    i_s = text.find('SUF = os.environ.get("TBL_SUFFIX", "")')
    i_m = text.find("def main():")
    if i_s < 0 or i_m < 0 or i_s > i_m:
        bad.append("SUF is not module-level / not above def main()")
    if "\n    SUF = os.environ.get" in text:
        bad.append("SUF is indented — it would be a local, not a global")
    # the guard must be called AFTER PAIRS is built
    i_p = text.find("PAIRS = [")
    i_g = text.find("_require_internal(PAIRS)")
    if i_p < 0 or i_g < 0 or i_p > i_g:
        bad.append("the guard runs before PAIRS is defined")
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
    print("137 — 22_symmetry_winners_curse.py: three suffix rules in one list")
    print("=" * 78)
    if not os.path.exists(a.target):
        print(f"[fatal] {a.target} not found.", file=sys.stderr)
        return 2
    src = open(a.target, errors="replace").read()
    if "import os" not in src:
        print("[fatal] no `import os`", file=sys.stderr)
        return 2
    new, problems = apply_edits(src)
    if problems:
        print("REFUSING — anchors do not match:")
        for p in problems:
            print(f"  ✗ {p}")
        print("\n  Nothing written. Send me lines 1-15 and 128-146.")
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
        print("  # 1. CALIBRATION — unsuffixed must reproduce symmetry_SE.tsv")
        print("  cp symmetry_SE.tsv symmetry_SE.PRESERVED")
        print("  python3 22_symmetry_winners_curse.py --event-type SE")
        print("  cmp symmetry_SE.tsv symmetry_SE.PRESERVED \\")
        print("    && echo '50bp BYTE-IDENTICAL' || diff symmetry_SE.PRESERVED "
              "symmetry_SE.tsv | head")
        print()
        print("  # 2. the guard must fire for a class 136 has not built")
        print("  TBL_SUFFIX=_rl125 python3 22_symmetry_winners_curse.py \\")
        print("      --event-type MXE 2>&1 | head -14")
        print()
        print("  # 3. 125 bp, SE — 136 built this one (job 408099)")
        print("  TBL_SUFFIX=_rl125 python3 22_symmetry_winners_curse.py \\")
        print("      --event-type SE")
        print("  cat symmetry_SE_rl125.tsv")
        print()
        print("⚠ THE UNSELECTED SLOPES ARE THE POINT. 22's ~0.56 headline is")
        print("  the mean of the UNSELECTED slopes (0.571/0.525/0.757/0.359/")
        print("  0.578 -> 0.558), NOT the selected ones (mean 0.414 at 50 bp,")
        print("  0.451 at 125 bp). Two different quantities; 22's is the one")
        print("  nobody has re-measured.")
        print()
        print("⬜ THEN the four transcriptions can be re-based: 82's")
        print("  GENOME_WIDE, 62's S5 warning, 98 line 393, and 22's own")
        print("  docstring. Until then they are 50 bp values.")
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
    print("✓ verified on disk: the MAGNet-internal arm moves on BOTH sides,")
    print("  the discovery stem on the x side only, the five replication")
    print("  stems not at all; all three outputs follow SUF; the guard")
    print("  refuses rather than letting the internal row drop out.")
    return 0


def self_test():
    import subprocess
    import tempfile

    fx = ('"""\n' + A_DOC + '"""\n'
          "import argparse, os\n"
          + A_SUF
          + "    ap = argparse.ArgumentParser()\n"
            '    ap.add_argument("--event-type", default="SE")\n'
            '    ap.add_argument("--program", action="store_true")\n'
            '    ap.add_argument("--out", default=None)\n'
            "    a = ap.parse_args()\n"
            "    ET = a.event_type\n"
          + A_OUT
          + A_PAIRS
          + '    print("OUT", out)\n'
            "    for lab, x, y in PAIRS:\n"
            '        print("PAIR", lab, "|", x, "|", y)\n'
            "if __name__ == '__main__':\n"
            "    main()\n")

    new, probs = apply_edits(fx)
    assert not probs, probs
    assert not verify(new), verify(new)
    compile(new, "fx", "exec")

    d = tempfile.mkdtemp()
    p = os.path.join(d, "fx.py")
    open(p, "w").write(new)
    os.makedirs(os.path.join(d, "splitctl"), exist_ok=True)

    def run(args=(), env=None):
        e = dict(os.environ)
        for k in ("TBL_SUFFIX", "ET_HINT"):
            e.pop(k, None)
        e.update(env or {})
        return subprocess.run([sys.executable, "-I", p, *args],
                              capture_output=True, text=True, cwd=d, env=e)

    # 1. unsuffixed -> every path and the output stem unchanged
    r = run()
    assert r.returncode == 0, r.stdout + r.stderr
    assert "OUT symmetry_SE\n" in r.stdout, r.stdout
    assert ("PAIR MAGNet internal | splitctl/SE_DCMvsA_1.all_events.tsv | "
            "splitctl/SE_HCMvsB_1.all_events.tsv") in r.stdout, r.stdout
    assert ("PAIR Colorado | SE_v2filtered.all_events.tsv | "
            "SE_Colorado_DCMvsNF.all_events.tsv") in r.stdout, r.stdout
    assert "Brigham_SE_HCM.all_events.tsv" in r.stdout, r.stdout

    # 2. suffixed WITHOUT 136's output -> REFUSE, naming 136 and the files.
    #    NB the guard runs after PAIRS and therefore BEFORE anything prints
    #    `out`, so the --program stem cannot be checked here — it is checked
    #    at step 5, once the arm exists. (An earlier version of this test
    #    asserted the stem against the refusal's output and failed.)
    r3 = run([], {"TBL_SUFFIX": "_rl125"})
    assert r3.returncode != 0, r3.stdout + r3.stderr
    b3 = r3.stdout + r3.stderr
    assert "REFUSING" in b3, b3
    assert "136_splitctl_rl125.sbatch" in b3, b3
    assert "splitctl/SE_DCMvsA_1_rl125.all_events.tsv" in b3, b3
    assert "splitctl/SE_HCMvsB_1_rl125.all_events.tsv" in b3, b3
    assert "nfA_1.ids" in b3, b3

    # 4. with 136's output present -> runs, and ONLY the right paths moved
    for f in ("splitctl/SE_DCMvsA_1_rl125.all_events.tsv",
              "splitctl/SE_HCMvsB_1_rl125.all_events.tsv"):
        open(os.path.join(d, f), "w").write("x\n")
    r4 = run([], {"TBL_SUFFIX": "_rl125"})
    assert r4.returncode == 0, r4.stdout + r4.stderr
    assert "OUT symmetry_SE_rl125" in r4.stdout, r4.stdout
    assert ("PAIR MAGNet internal | splitctl/SE_DCMvsA_1_rl125.all_events.tsv "
            "| splitctl/SE_HCMvsB_1_rl125.all_events.tsv") in r4.stdout, \
        r4.stdout
    assert ("PAIR Colorado | SE_v2filtered_rl125.all_events.tsv | "
            "SE_Colorado_DCMvsNF.all_events.tsv") in r4.stdout, r4.stdout
    # the five replication stems must be untouched
    for y in ("SE_Colorado_DCMvsNF.all_events.tsv",
              "SE_LBNL_DCMvsNF.all_events.tsv",
              "Brigham_SE_HCM.all_events.tsv",
              "Fuwai_SE_HCM.all_events.tsv",
              "Stanford_SE_HCM.all_events.tsv"):
        assert y in r4.stdout, (y, r4.stdout)
    for wrong in ("Colorado_DCMvsNF_rl125", "Brigham_SE_HCM_rl125",
                  "Fuwai_SE_HCM_rl125", "Stanford_SE_HCM_rl125",
                  "LBNL_DCMvsNF_rl125"):
        assert wrong not in r4.stdout, (wrong, r4.stdout)

    # 5. the --program variant keeps SUF LAST: symmetry_SE_program_rl125,
    #    not symmetry_SE_rl125_program
    r5 = run(["--program"], {"TBL_SUFFIX": "_rl125"})
    assert r5.returncode == 0, r5.stdout + r5.stderr
    assert "OUT symmetry_SE_program_rl125" in r5.stdout, r5.stdout

    # 6. an explicit --out still wins, suffix or not
    r6 = run(["--out", "mine"], {"TBL_SUFFIX": "_rl125"})
    assert "OUT mine" in r6.stdout, r6.stdout

    # 6. double application refused, text unchanged
    n6, p6 = apply_edits(new)
    assert p6, "second application was not refused"
    assert n6 == new, "a refused re-application still altered the text"

    print("[ok] self-test passed: unsuffixed paths and the output stem are "
          "unchanged; with TBL_SUFFIX the MAGNet-internal arm moves on BOTH "
          "sides and the discovery stem on the x side only, while all five "
          "replication stems stay put and none is wrongly suffixed; the "
          "output stem takes SUF last (symmetry_SE_program_rl125); a missing "
          "136 arm REFUSES and names both files, the sbatch and the .ids "
          "partition instead of silently dropping the row; an explicit --out "
          "still wins; double application refused")
    return 0


if __name__ == "__main__":
    sys.exit(main())
