#!/usr/bin/env python3
"""
126_patch_03_03b_readlength.py -- give 03 and 03b a read-length parameter, and
stop them overwriting ~320 MB of 50 bp comparators.

WHY
---
`56_rbm20_targets45.sbatch` reads `rbm20_${ET}${SUF}.events.tsv` via
`--from-events`. Those are 03's and 03b's outputs, so **56 is blocked on 03 and
03b**, neither of which was in 110's or 114's patched sets.

Both read `rmats_magnet/${ET}.MATS.JC.txt` and write unsuffixed:

    rbm20_SE.events.tsv            66.7 MB   7 Sep
    rbm20_MXE.events.tsv           93.8 MB   7 Sep
    rbm20_SE_techadj.events.tsv    66.6 MB   7 Sep
    rbm20_MXE_techadj.events.tsv   93.4 MB   7 Sep

and neither has a FORCE guard, so a 125 bp run overwrites all four plus their
reports. The eight `rbm20_*_g45.report.txt` / `rbm20_*_repro9.report.txt` files
from 5 Oct are 56's outputs built from those tables, so overwriting the events
tables orphans those as well. Third instance of the frozen-comparator pattern
after 23 and 11.

JC, NOT JCEC — DO NOT "FIX" THIS
--------------------------------
03's own comment states why:

    # JC (junction-only), NOT JCEC: JCEC folds exon-body reads into PSI and
    # body coverage scales with expression, which is the axis we are
    # correlating against.

This is the project's one documented JCEC exception. The patch repoints the
directory and leaves `.MATS.JC.txt` alone. `rmats_magnet_rl125/` carries the JC
files (SE 1.1 GB, MXE 381 MB, 13 Oct), so no new rMATS run is needed.

THE VARIABLE IS CALLED RLSUF, NOT SUF
-------------------------------------
`56` already uses `SUF` for the model variant (`"" | "_techadj"`). A read-length
variable called `SUF` would collide there and silently mis-name every input and
output. `RLSUF` is used here and must be used in 56 too.

NAME ORDER: the read-length suffix goes LAST
--------------------------------------------
    03   -> rbm20_${ET}${RLSUF}              e.g. rbm20_SE_rl125
    03b  -> rbm20_${ET}_techadj${RLSUF}      e.g. rbm20_SE_techadj_rl125

so 56's existing `rbm20_${ET}${SUF}` pattern extends to
`rbm20_${ET}${SUF}${RLSUF}` and resolves correctly for both variants.

THREE READERS IN 03b, NOT ONE
-----------------------------
03b's comparison section (line ~77) greps `rbm20_${ET}.report.txt` and
`rbm20_${ET}_techadj.report.txt`. Unsuffixed, a 125 bp run would print its
"before vs after" comparison from the September reports. Both are patched.

USAGE
-----
  python 126_patch_03_03b_readlength.py
  python 126_patch_03_03b_readlength.py --apply
  python 126_patch_03_03b_readlength.py --self-test
"""
import argparse
import difflib
import os
import shutil
import sys
import time


def guard(prefix_expr, label):
    """A don't-clobber guard for one output-prefix expression."""
    return f'''
# --------------------------------------------------------------------------
# A wrong rMATS directory must be an ERROR, not a skip: the per-class guard
# below only `continue`s, which turns a mistyped directory into a job that
# completes having done nothing and reports no failure.
if [ ! -d "$RMATS_DIR" ]; then
  echo "*** rMATS directory not found: $RMATS_DIR"
  echo "*** TBL_SUFFIX='$RLSUF'. Set RMATS_DIR explicitly if the directory for"
  echo "*** this read length has another name."
  exit 1
fi

# DON'T-CLOBBER GUARD ({label}).
# The 7 Sep tables are the 50 bp comparators, and 56's *_g45 / *_repro9 reports
# were built from them. Overwriting them orphans those too, and nothing errors.
for ET in SE MXE; do
  for EXT in events.tsv report.txt; do
    F="{prefix_expr}.${{EXT}}"
    if [ -e "$F" ] && [ -z "${{OVERWRITE_OK:-}}" ]; then
      echo "*** REFUSING — $F already exists"
      echo "***   modified: $(date -r "$F" '+%Y-%m-%d %H:%M' 2>/dev/null)"
      echo "***   size:     $(stat -c %s "$F" 2>/dev/null) bytes"
      echo "***"
      echo "*** Set OVERWRITE_OK=1 only if you mean to replace it."
      exit 1
    fi
  done
done
echo "[guard] {label}: no existing output for TBL_SUFFIX='${{RLSUF}}' — safe to write"
# --------------------------------------------------------------------------
'''


HEAD = '''# --------------------------------------------------------------------------
# Read length. RLSUF drives the rMATS directory and the output prefix together,
# so a run is consistent at ONE read length or it does not run. Empty = the
# original 50 bp behaviour, reproducing the 7 Sep output.
#
# Called RLSUF, not SUF: 56_rbm20_targets45.sbatch already uses SUF for the
# model variant ("" | "_techadj"), and a collision there would silently
# mis-name every input and output.
RLSUF="${TBL_SUFFIX:-}"
RMATS_DIR="${RMATS_DIR:-rmats_magnet${RLSUF}}"
# --------------------------------------------------------------------------
'''

A03_LOOP = (
    "# JC (junction-only), NOT JCEC: JCEC folds exon-body reads into PSI and "
    "body\n"
    "# coverage scales with expression, which is the axis we are correlating "
    "against.\n"
    "for ET in SE MXE; do\n")

A03B_LOOP = ("# TIN is silently ignored by the script if the column is "
             "absent.\n"
             "for ET in SE MXE; do\n")

MATS_OLD = "  MATS=rmats_magnet/${ET}.MATS.JC.txt\n"
MATS_NEW = "  MATS=${RMATS_DIR}/${ET}.MATS.JC.txt\n"

FILES = {
    "03_rbm20_psi.sbatch": [
        ("RLSUF, RMATS_DIR and the guard",
         A03_LOOP,
         HEAD + guard("rbm20_${ET}${RLSUF}", "03, unadjusted") + A03_LOOP, 1),
        ("the rMATS input path", MATS_OLD, MATS_NEW, 1),
        ("the output prefix",
         "    --out-prefix rbm20_${ET} \\\n",
         "    --out-prefix rbm20_${ET}${RLSUF} \\\n", 1),
    ],
    "03b_rbm20_psi_techadj.sbatch": [
        ("RLSUF, RMATS_DIR and the guard",
         A03B_LOOP,
         HEAD + guard("rbm20_${ET}_techadj${RLSUF}", "03b, techadj")
         + A03B_LOOP, 1),
        ("the rMATS input path", MATS_OLD, MATS_NEW, 1),
        ("the output prefix",
         "    --out-prefix rbm20_${ET}_techadj\n",
         "    --out-prefix rbm20_${ET}_techadj${RLSUF}\n", 1),
        ("the comparison reader, unadjusted report",
         '  grep -h "matched-bg median" rbm20_${ET}.report.txt 2>/dev/null '
         '| head -2\n',
         '  grep -h "matched-bg median" rbm20_${ET}${RLSUF}.report.txt '
         '2>/dev/null | head -2\n', 1),
        ("the comparison reader, techadj report",
         '  grep -h "matched-bg median" rbm20_${ET}_techadj.report.txt '
         '2>/dev/null | head -2\n',
         '  grep -h "matched-bg median" rbm20_${ET}_techadj${RLSUF}.report.txt '
         '2>/dev/null | head -2\n', 1),
    ],
}


def apply_edits(text, edits):
    problems, out = [], text
    for name, anchor, repl, want in edits:
        got = out.count(anchor)
        if got != want:
            problems.append(f"{name}: anchor occurs {got}x, expected {want}x\n"
                            f"      anchor: "
                            f"{anchor.strip().splitlines()[-1][:76]}")
            continue
        if 'RLSUF="${TBL_SUFFIX:-}"' in out and name.startswith("RLSUF"):
            problems.append("RLSUF is already defined — already patched")
            continue
        out = out.replace(anchor, repl, want)
    return out, problems


def verify(text, fname):
    bad = []
    if 'RLSUF="${TBL_SUFFIX:-}"' not in text:
        bad.append("RLSUF is not defined")
    if 'RMATS_DIR="${RMATS_DIR:-rmats_magnet${RLSUF}}"' not in text:
        bad.append("RMATS_DIR is not derived from RLSUF")
    if "rmats_magnet/" in text:
        bad.append("a hard-coded rmats_magnet/ path survives")
    if ".MATS.JC.txt" not in text:
        bad.append("the JC input was lost — it must stay JC, not JCEC")
    if ".MATS.JCEC.txt" in text:
        bad.append("JCEC appeared; this test is deliberately JC-only")
    if "OVERWRITE_OK" not in text:
        bad.append("the don't-clobber guard is missing")
    if 'if [ ! -d "$RMATS_DIR" ]; then' not in text:
        bad.append("the rMATS directory check is missing")
    if "03b" in fname:
        for m in ("--out-prefix rbm20_${ET}_techadj${RLSUF}",
                  "rbm20_${ET}${RLSUF}.report.txt",
                  "rbm20_${ET}_techadj${RLSUF}.report.txt"):
            if m not in text:
                bad.append(f"missing in 03b: {m}")
        if 'rbm20_${ET}.report.txt' in text:
            bad.append("an unsuffixed report reader survives in 03b")
    else:
        if "--out-prefix rbm20_${ET}${RLSUF}" not in text:
            bad.append("the 03 output prefix is not suffixed")
    i = text.find('RLSUF="${TBL_SUFFIX:-}"')
    j = text.find("--out-prefix rbm20_")
    if i < 0 or j < 0 or i > j:
        bad.append("RLSUF is not defined above the output prefix")
    g = text.find("OVERWRITE_OK")
    if g < 0 or g > j:
        bad.append("the guard does not precede the writing loop")
    return bad


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--apply", action="store_true")
    ap.add_argument("--self-test", action="store_true")
    ap.add_argument("--only", choices=list(FILES))
    a = ap.parse_args()
    if a.self_test:
        return self_test()

    self_test()
    print()
    print("=" * 78)
    print("126 — read-length parameter and output guard for 03 and 03b")
    print("=" * 78)

    targets = [a.only] if a.only else list(FILES)
    staged = {}
    for fname in targets:
        print(f"\n═══ {fname}")
        if not os.path.exists(fname):
            print(f"  [fatal] not found", file=sys.stderr)
            return 2
        src = open(fname, errors="replace").read()
        new, problems = apply_edits(src, FILES[fname])
        if problems:
            print("  REFUSING — anchors do not match:")
            for p in problems:
                print(f"    ✗ {p}")
            return 3
        bad = verify(new, fname)
        if bad:
            print("  REFUSING — verification failed:")
            for b in bad:
                print(f"    ✗ {b}")
            return 3
        print("".join(difflib.unified_diff(
            src.splitlines(True), new.splitlines(True),
            fromfile=fname, tofile=fname + " (patched)", n=2)))
        staged[fname] = new

    print("=" * 78)
    if not a.apply:
        print("DRY RUN — nothing written. Re-run with --apply, then:")
        print()
        print("  bash -n 03_rbm20_psi.sbatch 03b_rbm20_psi_techadj.sbatch")
        print()
        print("  # prove the guard fires on the 7 Sep files")
        print("  bash 03_rbm20_psi.sbatch 2>&1 | head -15")
        print()
        print("  # then the real runs")
        print("  TBL_SUFFIX=_rl125 sbatch 03_rbm20_psi.sbatch")
        print("  TBL_SUFFIX=_rl125 sbatch 03b_rbm20_psi_techadj.sbatch")
        print()
        print("  They write rbm20_{SE,MXE}_rl125.* and")
        print("  rbm20_{SE,MXE}_techadj_rl125.*, so the 7 Sep tables and 56's")
        print("  *_g45 / *_repro9 reports are untouched.")
        print()
        print("⚠ 03 passes --n-perm 200. Check the wall time of the 7 Sep run")
        print("  before assuming these are quick; SE.MATS.JC.txt is 1.1 GB at")
        print("  125 bp against whatever the 50 bp file was.")
        print()
        print("⚠ 56 still needs its own patch, with RLSUF (NOT SUF — it uses")
        print("  SUF for the techadj variant) and the suffix LAST:")
        print("    EV=\"rbm20_${ET}${SUF}${RLSUF}.events.tsv\"")
        return 0

    for fname, new in staged.items():
        bak = f"{fname}.bak.{time.strftime('%Y%m%d_%H%M%S')}"
        shutil.copy2(fname, bak)
        open(fname, "w").write(new)
        print(f"[write] {fname}   (backup: {bak})")
        disk = open(fname, errors="replace").read()
        bad = verify(disk, fname)
        if bad or disk != new:
            print("✗ ON-DISK VERIFICATION FAILED — restoring:")
            for b in bad:
                print(f"    {b}")
            shutil.copy2(bak, fname)
            return 3
    print("✓ both verified on disk.")
    return 0


def self_test():
    import subprocess
    import tempfile

    fx03 = ("#!/bin/bash\n"
            + A03_LOOP
            + '  echo "=================== $ET ==================="\n'
            + MATS_OLD
            + '  echo "MATS=$MATS"\n'
            + '  if [ ! -s "$MATS" ]; then echo "missing $MATS -- skipping"; '
              'continue; fi\n'
            + '  echo stub \\\n'
            + "    --out-prefix rbm20_${ET} \\\n"
            + "    --n-perm 200\n"
            + "done\n")

    fx03b = ("#!/bin/bash\n"
             + A03B_LOOP
             + '  MATS_UNUSED=1\n'
             + MATS_OLD
             + '  echo "MATS=$MATS"\n'
             + '  [ -s "$MATS" ] || { echo "missing $MATS -- skipping"; '
               'continue; }\n'
             + '  echo stub \\\n'
             + "    --out-prefix rbm20_${ET}_techadj\n"
             + "done\n"
             + 'echo "--- compare against the unadjusted run ---"\n'
             + "for ET in SE MXE; do\n"
             + '  grep -h "matched-bg median" rbm20_${ET}.report.txt '
               '2>/dev/null | head -2\n'
             + '  grep -h "matched-bg median" rbm20_${ET}_techadj.report.txt '
               '2>/dev/null | head -2\n'
             + "done\n")

    d = tempfile.mkdtemp()
    # The per-class guard is `[ -s "$MATS" ]`, so the stubs must be NON-EMPTY
    # or every class is skipped and the fixture tests nothing. The first
    # version of this self-test created empty directories and asserted on
    # output that never happened.
    for sub in ("rmats_magnet", "rmats_magnet_rl125"):
        os.makedirs(os.path.join(d, sub))
        for et in ("SE", "MXE"):
            with open(os.path.join(d, sub, f"{et}.MATS.JC.txt"), "w") as fh:
                fh.write("stub\n")

    for fname, fx in (("03_rbm20_psi.sbatch", fx03),
                      ("03b_rbm20_psi_techadj.sbatch", fx03b)):
        new, probs = apply_edits(fx, FILES[fname])
        assert not probs, (fname, probs)
        assert not verify(new, fname), (fname, verify(new, fname))
        p = os.path.join(d, fname)
        open(p, "w").write(new)

        def run(env=None):
            return subprocess.run(["bash", p], capture_output=True, text=True,
                                  cwd=d, env=dict(os.environ, **(env or {})))

        # unsuffixed -> 50 bp paths
        r = run()
        assert r.returncode == 0, r.stdout + r.stderr
        assert "rmats_magnet/SE.MATS.JC.txt" in r.stdout, (fname, r.stdout)
        assert ".MATS.JCEC" not in r.stdout, "must stay JC"
        assert "[guard]" in r.stdout, r.stdout

        # suffixed -> 125 bp paths, including the output
        r2 = run({"TBL_SUFFIX": "_rl125"})
        assert r2.returncode == 0, r2.stdout + r2.stderr
        assert "rmats_magnet_rl125/SE.MATS.JC.txt" in r2.stdout, r2.stdout
        want = ("rbm20_${ET}_techadj_rl125" if "03b" in fname
                else "rbm20_SE_rl125")
        if "03b" in fname:
            assert "--out-prefix rbm20_SE_techadj_rl125" in r2.stdout, r2.stdout
        else:
            assert "--out-prefix rbm20_SE_rl125" in r2.stdout, r2.stdout

        # an existing output refuses, before the loop
        out_name = ("rbm20_SE_techadj.events.tsv" if "03b" in fname
                    else "rbm20_SE.events.tsv")
        op = os.path.join(d, out_name)
        open(op, "w").write("x")
        r3 = run()
        assert r3.returncode == 1, (fname, r3.stdout)
        assert "REFUSING" in r3.stdout and out_name in r3.stdout, r3.stdout
        assert "--out-prefix" not in r3.stdout, "reached the loop anyway"
        # the suffixed run is unaffected
        assert run({"TBL_SUFFIX": "_rl125"}).returncode == 0
        # OVERWRITE_OK releases it
        assert run({"OVERWRITE_OK": "1"}).returncode == 0
        os.remove(op)

        # a missing rMATS dir errors rather than skipping
        r5 = run({"TBL_SUFFIX": "_nope"})
        assert r5.returncode == 1 and "not found" in r5.stdout, r5.stdout

        # double application refused
        _n, p6 = apply_edits(new, FILES[fname])
        assert p6, f"{fname}: second application not refused"
        _ = want

    # 03b's comparison readers must be suffixed too
    new03b, _ = apply_edits(fx03b, FILES["03b_rbm20_psi_techadj.sbatch"])
    assert 'rbm20_${ET}${RLSUF}.report.txt' in new03b
    assert 'rbm20_${ET}_techadj${RLSUF}.report.txt' in new03b
    assert 'rbm20_${ET}.report.txt' not in new03b

    print("[ok] self-test passed: both files take RLSUF (not SUF), keep "
          "MATS.JC, derive input and output together, refuse an existing "
          "output before the loop, error on a missing rMATS directory, release "
          "on OVERWRITE_OK, and 03b's two comparison readers are suffixed as "
          "well; double application refused for both")
    return 0


if __name__ == "__main__":
    sys.exit(main())
