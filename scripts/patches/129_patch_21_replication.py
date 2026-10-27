#!/usr/bin/env python3
"""
129_patch_21_replication.py -- 21 warns about the 125 bp table and then passes
the 50 bp one.

THE DEFECTS
-----------
**1. The warning was patched; the argument was not.** Line 130:

    [ -s "${ET}_v2filtered_rl125.all_events.tsv" ] || \\
      echo "WARNING: ${ET}_v2filtered_rl125.all_events.tsv missing -- run 04 first"

but lines 134/136 call

    python 21_replication_v2.py --rep Colorado --suffix _Colorado_DCMvsNF

with **no `--disc-suffix`**, so it falls back to its own default `_v2filtered` —
50 bp. 110 renamed the check and left the reader. Same shape as 23, where the
rMATS directory was repointed and the gene list and output prefix were not.

`21_replication_v2.py` needs **no patch**: `--disc-suffix` is already an
argument (line 109, default `_v2filtered`).

**2. The output name carries no read length.** Line 114 of the `.py`:

    out = a.out or f"replication_v2_{a.rep}"

so a 125 bp run overwrites `replication_v2_Colorado.{tsv,report.txt}` and
`replication_v2_LBNL.*`. Those have four real consumers —
`06_sensitivity_concordance.py`, `69_build_s1b_ladder.py`,
`74_build_s8a_permutations.py` and `98_build_replication_table.py` (Table 4).
The patch passes `--out` explicitly so the two read lengths coexist.

⬜ Those four will need the suffixed name when they run. That is track-2 work
and 98 is already on the list for the retired 0.521 slopes.

**3. The LBNL loop has no guard.** Colorado's has

    if [ -s "$OUT" ] && [ "$FORCE_COLORADO" != "1" ]; then continue

and LBNL's does not, so every run re-executes phase2 five times and rewrites
`${ET}_LBNL_DCMvsNF.all_events.tsv`. The inputs are unchanged so it should
reproduce, but it is five needless phase2 runs and a regeneration of a 50 bp
comparator. `FORCE_LBNL` mirrors the Colorado switch.

WHAT IS CORRECT ALREADY AND MUST NOT BE "FIXED"
-----------------------------------------------
`rmats_LBNL/` (line 90) and `rmats_Colorado/` (line 113) are **deliberately
unsuffixed**: neither cohort was re-quantified, and only MAGNet was. A run that
mixes a 125 bp discovery table against 50 bp cohort tables is the intended
design — that is what cross-cohort replication means here.

THE VARIABLE IS RLSUF, NOT SUF
------------------------------
`21_replication_v2.py --suffix` means the **cohort** (`_Colorado_DCMvsNF`), not
the read length. Same trap as 56, where `SUF` meant the techadj variant.

USAGE
-----
  python 129_patch_21_replication.py
  python 129_patch_21_replication.py --apply
  python 129_patch_21_replication.py --self-test

Then:
  TBL_SUFFIX=_rl125 sbatch 21_replication_v2.sbatch
"""
import argparse
import difflib
import os
import shutil
import sys
import time

TARGET = "21_replication_v2.sbatch"

A_ROOT = "ROOT=/data/project/chureeratp/RNA/HEART\n"
B_ROOT = A_ROOT + '''
# --------------------------------------------------------------------------
# Read length. Called RLSUF, NOT SUF: 21_replication_v2.py's --suffix means the
# COHORT (_Colorado_DCMvsNF), and a collision would mis-name every table.
#
# Only the MAGNet DISCOVERY table carries this suffix. rmats_LBNL/ and
# rmats_Colorado/ below are deliberately unsuffixed — neither cohort was
# re-quantified, and a 125 bp discovery set against 50 bp cohort tables is what
# cross-cohort replication means here.
RLSUF="${TBL_SUFFIX:-}"

# Mirrors FORCE_COLORADO. Without it the LBNL loop re-runs phase2 five times on
# every invocation and rewrites ${ET}_LBNL_DCMvsNF.all_events.tsv for nothing.
FORCE_LBNL="${FORCE_LBNL:-0}"
# --------------------------------------------------------------------------
'''

A_LBNL = ("for ET in SE MXE A3SS A5SS RI; do\n"
          "  M=rmats_LBNL/${ET}.MATS.JCEC.txt\n")
B_LBNL = ("for ET in SE MXE A3SS A5SS RI; do\n"
          "  OUT=${ET}_LBNL_DCMvsNF.all_events.tsv\n"
          '  if [ -s "$OUT" ] && [ "$FORCE_LBNL" != "1" ]; then\n'
          '    echo "----- $ET : reusing $OUT (FORCE_LBNL=1 to rebuild) -----"\n'
          "    continue\n"
          "  fi\n"
          "  M=rmats_LBNL/${ET}.MATS.JCEC.txt\n")

A_WARN = ('  [ -s "${ET}_v2filtered_rl125.all_events.tsv" ] || \\\n'
          '    echo "WARNING: ${ET}_v2filtered_rl125.all_events.tsv missing '
          '-- run 04 first"\n')
B_WARN = ('  [ -s "${ET}_v2filtered${RLSUF}.all_events.tsv" ] || \\\n'
          '    echo "WARNING: ${ET}_v2filtered${RLSUF}.all_events.tsv missing '
          '-- run 04 first"\n')

A_COL = ("python 21_replication_v2.py --rep Colorado "
         "--suffix _Colorado_DCMvsNF\n")
B_COL = ("python 21_replication_v2.py --rep Colorado "
         "--suffix _Colorado_DCMvsNF \\\n"
         "    --disc-suffix _v2filtered${RLSUF} \\\n"
         "    --out replication_v2_Colorado${RLSUF}\n")

A_LB = ("python 21_replication_v2.py --rep LBNL     "
        "--suffix _LBNL_DCMvsNF\n")
B_LB = ("python 21_replication_v2.py --rep LBNL     "
        "--suffix _LBNL_DCMvsNF \\\n"
        "    --disc-suffix _v2filtered${RLSUF} \\\n"
        "    --out replication_v2_LBNL${RLSUF}\n")

EDITS = [
    ("RLSUF and FORCE_LBNL", A_ROOT, B_ROOT, 1),
    ("the LBNL don't-rebuild guard", A_LBNL, B_LBNL, 1),
    ("the discovery-table warning, derived from RLSUF", A_WARN, B_WARN, 1),
    ("the Colorado call: --disc-suffix and --out", A_COL, B_COL, 1),
    ("the LBNL call: --disc-suffix and --out", A_LB, B_LB, 1),
]


def apply_edits(text):
    problems, out = [], text
    for name, anchor, repl, want in EDITS:
        got = out.count(anchor)
        if got != want:
            problems.append(f"{name}: anchor occurs {got}x, expected {want}x\n"
                            f"      anchor: "
                            f"{anchor.strip().splitlines()[0][:76]}")
            continue
        if 'RLSUF="${TBL_SUFFIX:-}"' in out and name.startswith("RLSUF"):
            problems.append("RLSUF is already defined — already patched")
            continue
        out = out.replace(anchor, repl, want)
    return out, problems


def verify(text):
    bad = []
    for w in ('RLSUF="${TBL_SUFFIX:-}"',
              'FORCE_LBNL="${FORCE_LBNL:-0}"',
              '[ "$FORCE_LBNL" != "1" ]',
              '${ET}_v2filtered${RLSUF}.all_events.tsv',
              "--disc-suffix _v2filtered${RLSUF}",
              "--out replication_v2_Colorado${RLSUF}",
              "--out replication_v2_LBNL${RLSUF}"):
        if w not in text:
            bad.append(f"missing: {w[:56]}")
    if "_v2filtered_rl125" in text:
        bad.append("a hard-coded _v2filtered_rl125 survives")
    # the cohort rMATS dirs must stay unsuffixed
    for d in ("rmats_LBNL/", "rmats_Colorado/"):
        if d not in text:
            bad.append(f"{d} was lost")
        if d.rstrip("/") + "${RLSUF}" in text or d.rstrip("/") + "_rl125" in text:
            bad.append(f"{d} was wrongly suffixed — that cohort was not "
                       f"re-quantified")
    # both calls must carry both new arguments
    if text.count("--disc-suffix _v2filtered${RLSUF}") != 2:
        bad.append("--disc-suffix is not on both calls")
    # RLSUF above every consumer
    i = text.find('RLSUF="${TBL_SUFFIX:-}"')
    for m in ("${ET}_v2filtered${RLSUF}", "--disc-suffix _v2filtered${RLSUF}",
              '[ "$FORCE_LBNL" != "1" ]'):
        j = text.find(m)
        if j < 0 or i > j:
            bad.append(f"RLSUF/FORCE_LBNL not defined above: {m[:40]}")
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
    print("129 — 21: pass the 125 bp discovery table it already checks for")
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
        print("\n  Nothing written. Send me lines 33-40, 88-92 and 128-137.")
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
        print("  bash -n 21_replication_v2.sbatch")
        print()
        print("  TBL_SUFFIX=_rl125 sbatch 21_replication_v2.sbatch")
        print()
        print("  It writes replication_v2_{Colorado,LBNL}_rl125.{tsv,report.txt}")
        print("  and reuses both cohorts' existing phase2 tables.")
        print()
        print("⚠ CROSS-CHECK, and this is the point of running it: 33's")
        print("  section D already gives the pooled concordance at 125 bp as")
        print("  1,376/1,800 = 76.4% (50 bp: 1,636/2,163 = 75.6%). 21 reaches")
        print("  the same quantity by a different route. If Colorado's row")
        print("  does not match 33's Colorado row (311/385, 80.8%), one of the")
        print("  two is wrong and that matters more than either number.")
        print()
        print("⚠ FOUR DOWNSTREAM SCRIPTS read replication_v2_* and will need")
        print("  the suffixed name — track 2:")
        print("    06_sensitivity_concordance.py")
        print("    69_build_s1b_ladder.py")
        print("    74_build_s8a_permutations.py")
        print("    98_build_replication_table.py   (Table 4; also carries the")
        print("                                     retired 0.521 slopes)")
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
    print("✓ verified on disk: the discovery suffix reaches both calls, both")
    print("  outputs are suffixed, LBNL is guarded, and the cohort rMATS")
    print("  directories are untouched.")
    return 0


def self_test():
    import subprocess
    import tempfile

    fx = ("#!/bin/bash\n"
          + A_ROOT
          + "FORCE_COLORADO=0\n"
          + 'echo "--- LBNL"\n'
          + A_LBNL
          + '  [ -s "$M" ] || { echo "missing $M -- skipping $ET"; continue; }\n'
          + '  echo "LBNL phase2 $ET -> ${ET}_LBNL_DCMvsNF"\n'
          + "done\n"
          + 'echo "--- Colorado"\n'
          + "for ET in SE; do\n"
          + "  M=rmats_Colorado/${ET}.MATS.JCEC.txt\n"
          + '  echo "COL $M"\n'
          + "done\n"
          + 'echo "--- concordance"\n'
          + "for ET in SE; do\n"
          + A_WARN
          + "done\n"
          + A_COL
          + A_LB)
    # Patch FIRST, then make the calls echoable. Substituting before the patch
    # destroys the very anchors being tested — which is what the first version
    # of this self-test did.
    new, probs = apply_edits(fx)
    assert not probs, probs
    assert not verify(new), verify(new)

    d = tempfile.mkdtemp()
    p = os.path.join(d, "fx.sh")
    open(p, "w").write(new.replace("python 21_replication_v2.py", "echo CALL"))
    os.makedirs(os.path.join(d, "rmats_LBNL"))
    for et in ("SE", "MXE", "A3SS", "A5SS", "RI"):
        with open(os.path.join(d, "rmats_LBNL", f"{et}.MATS.JCEC.txt"), "w") as f:
            f.write("stub\n")

    def run(env=None):
        return subprocess.run(["bash", p], capture_output=True, text=True,
                              cwd=d, env=dict(os.environ, **(env or {})))

    # 1. unsuffixed: original names, LBNL builds (no outputs exist yet)
    r = run()
    assert r.returncode == 0, r.stdout + r.stderr
    assert "CALL --rep Colorado --suffix _Colorado_DCMvsNF "\
           "--disc-suffix _v2filtered --out replication_v2_Colorado" \
           in r.stdout, r.stdout
    assert "--out replication_v2_LBNL\n" in r.stdout, r.stdout
    assert "SE_v2filtered.all_events.tsv missing" in r.stdout, r.stdout
    assert "LBNL phase2 SE" in r.stdout, r.stdout
    assert "rmats_Colorado/SE" in r.stdout, "cohort dir must stay unsuffixed"

    # 2. suffixed: discovery and outputs move, cohort dirs do not
    r2 = run({"TBL_SUFFIX": "_rl125"})
    assert r2.returncode == 0, r2.stdout + r2.stderr
    assert "--disc-suffix _v2filtered_rl125" in r2.stdout, r2.stdout
    assert "--out replication_v2_Colorado_rl125" in r2.stdout, r2.stdout
    assert "--out replication_v2_LBNL_rl125" in r2.stdout, r2.stdout
    assert "SE_v2filtered_rl125.all_events.tsv missing" in r2.stdout, r2.stdout
    assert "rmats_LBNL/SE" not in r2.stdout or "_rl125/SE" not in r2.stdout
    assert "rmats_Colorado/SE" in r2.stdout, r2.stdout

    # 3. the LBNL guard: once the output exists it is reused
    open(os.path.join(d, "SE_LBNL_DCMvsNF.all_events.tsv"), "w").write("x")
    r3 = run()
    assert r3.returncode == 0, r3.stdout + r3.stderr
    assert "reusing SE_LBNL_DCMvsNF.all_events.tsv" in r3.stdout, r3.stdout
    assert "LBNL phase2 SE\n" not in r3.stdout, r3.stdout
    # the other four classes still build
    assert "LBNL phase2 MXE" in r3.stdout, r3.stdout
    # FORCE_LBNL=1 rebuilds
    r4 = run({"FORCE_LBNL": "1"})
    assert "LBNL phase2 SE" in r4.stdout, r4.stdout

    _n, p5 = apply_edits(new)
    assert p5, "second application was not refused"

    print("[ok] self-test passed: --disc-suffix and --out both reach both "
          "calls and move with TBL_SUFFIX; the discovery warning derives from "
          "RLSUF; rmats_LBNL/ and rmats_Colorado/ stay unsuffixed; the LBNL "
          "guard reuses an existing table per class and FORCE_LBNL=1 rebuilds; "
          "double application refused")
    return 0


if __name__ == "__main__":
    sys.exit(main())
