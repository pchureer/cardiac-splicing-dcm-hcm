#!/usr/bin/env python3
"""
127_patch_56_readlength.py -- read-length parameter for 56, at all seven sites.

WHY ALL SEVEN
-------------
`grep -n 'rbm20_' 56_rbm20_targets45.sbatch` enumerates every reader and writer.
Reading windows of the file found three of them; the grep found seven. That is
the method rule 8 actually requires — enumerate, do not sample.

  105  EV="rbm20_${ET}${SUF}.events.tsv"              discovery read
  135  local ev="rbm20_${tag}.events.tsv"             read
  136  local pref="rbm20_${tag}${suf}"                WRITE
  183  ev = f"rbm20_{tag}.events.tsv"                 read, inside <<PY
  258  parse(f"rbm20_{tag}_repro9.report.txt")        read, inside <<'PY'
  259  parse(f"rbm20_{tag}_g45.report.txt")           read, inside <<'PY'
  262  printed "read rbm20_{tag}_{repro9,g45}..."     instruction to the user
  294  printed "read rbm20_{tag}_g45.report.txt..."   instruction to the user

Lines 118, 328 and 330 are globs or prose about `rbm20_*` and are left alone:
a glob matches both read lengths, which is what a diagnostic listing wants.

TWO HEREDOCS, TWO MECHANISMS
----------------------------
Line 178 is `python3 - $RUNS <<PY` — **unquoted**, so bash substitutes into it
(it already relies on that for `$G45` and `$NINE`). `$RLSUF` works directly.

Line 206 is `python3 - $RUNS <<'PY'` — **quoted**, so nothing is substituted.
`RLSUF` must be passed as an environment variable and read with `os.environ`.
Writing `$RLSUF` inside that block would emit a literal dollar sign into the
f-string and the script would look for a file called `rbm20_SE$RLSUF...`.

THE VARIABLE IS RLSUF, NOT SUF
------------------------------
Line 104 is `for SUF in "" "_techadj"` — `SUF` is the MODEL variant in this
file. A read-length variable called `SUF` would be overwritten by that loop on
its first iteration and every path would be silently wrong. `RLSUF` matches
what 126 put into 03 and 03b.

THE SUFFIX GOES LAST
--------------------
    03   writes  rbm20_SE_rl125.events.tsv
    03b  writes  rbm20_SE_techadj_rl125.events.tsv

so `rbm20_${ET}${SUF}${RLSUF}` resolves to both. `RUNS` keeps carrying the bare
tag (`SE`, `SE_techadj`) because `run_report` derives `et="${tag%%_*}"` from it —
RLSUF is applied where filenames are built, never folded into the tag.

56's own FORCE guard (`[ -s "${pref}.report.txt" ] && [ "$FORCE" != "1" ]`)
starts protecting the new outputs automatically, since `pref` now carries RLSUF.

USAGE
-----
  python 127_patch_56_readlength.py
  python 127_patch_56_readlength.py --apply
  python 127_patch_56_readlength.py --self-test

Then, once 03 and 03b have produced the 125 bp events tables:
  TBL_SUFFIX=_rl125 sbatch 56_rbm20_targets45.sbatch
"""
import argparse
import difflib
import os
import shutil
import sys
import time

TARGET = "56_rbm20_targets45.sbatch"

A_SCRIPT = "SCRIPT=03_rbm20_psi_correlation.py\n"
B_SCRIPT = '''# --------------------------------------------------------------------------
# Read length. Called RLSUF, NOT SUF: line ~104 uses SUF for the model variant
# ("" | "_techadj") and would overwrite a read-length SUF on its first
# iteration, silently mis-naming every path.
#
# The suffix goes LAST, matching what 126 put into 03 and 03b:
#     03   -> rbm20_SE_rl125.events.tsv
#     03b  -> rbm20_SE_techadj_rl125.events.tsv
# so rbm20_${ET}${SUF}${RLSUF} resolves to both. RUNS keeps the bare tag,
# because run_report derives et="${tag%%_*}" from it.
RLSUF="${TBL_SUFFIX:-}"
# --------------------------------------------------------------------------
''' + A_SCRIPT

EDITS = [
    ("RLSUF", A_SCRIPT, B_SCRIPT, 1),

    ("105 discovery read",
     '    EV="rbm20_${ET}${SUF}.events.tsv"\n',
     '    EV="rbm20_${ET}${SUF}${RLSUF}.events.tsv"\n', 1),

    ("135 run_report read",
     '  local ev="rbm20_${tag}.events.tsv"\n',
     '  local ev="rbm20_${tag}${RLSUF}.events.tsv"\n', 1),

    ("136 run_report WRITE prefix",
     '  local pref="rbm20_${tag}${suf}"\n',
     '  local pref="rbm20_${tag}${suf}${RLSUF}"\n', 1),

    # unquoted heredoc: bash substitutes $RLSUF into the python source
    ("183 section-4 read (unquoted heredoc, $RLSUF substituted by bash)",
     '    ev = f"rbm20_{tag}.events.tsv"\n',
     '    ev = f"rbm20_{tag}$RLSUF.events.tsv"\n', 1),

    # quoted heredoc: pass RLSUF through the environment instead
    ("206 pass RLSUF into the quoted heredoc",
     "python3 - $RUNS <<'PY'\n",
     'RLSUF="$RLSUF" python3 - $RUNS <<\'PY\'\n', 1),

    ("207 read RLSUF from the environment inside the quoted block",
     "import os, re, sys\n",
     "import os, re, sys\n"
     "\n"
     "# The heredoc above is QUOTED, so bash substitutes nothing into this\n"
     "# block. RLSUF arrives through the environment; a literal $RLSUF here\n"
     "# would be looked for as part of the filename.\n"
     "RLSUF = os.environ.get(\"RLSUF\", \"\")\n", 1),

    ("258 section-5 read, repro9",
     '    s9, c9   = parse(f"rbm20_{tag}_repro9.report.txt")\n',
     '    s9, c9   = parse(f"rbm20_{tag}_repro9{RLSUF}.report.txt")\n', 1),

    ("259 section-5 read, g45",
     '    s45, c45 = parse(f"rbm20_{tag}_g45.report.txt")\n',
     '    s45, c45 = parse(f"rbm20_{tag}_g45{RLSUF}.report.txt")\n', 1),

    ("262 printed instruction",
     '              f"read rbm20_{tag}_{{repro9,g45}}.report.txt by hand")\n',
     '              f"read rbm20_{tag}_{{repro9,g45}}{RLSUF}.report.txt '
     'by hand")\n', 1),

    ("294 printed instruction",
     '              f" section 4 -- read rbm20_{tag}_g45.report.txt section 4 '
     'and compare by hand")\n',
     '              f" section 4 -- read rbm20_{tag}_g45{RLSUF}.report.txt '
     'section 4 and compare by hand")\n', 1),
]


def apply_edits(text):
    problems, out = [], text
    for name, anchor, repl, want in EDITS:
        got = out.count(anchor)
        if got != want:
            problems.append(f"{name}: anchor occurs {got}x, expected {want}x\n"
                            f"      anchor: "
                            f"{anchor.strip().splitlines()[-1][:76]}")
            continue
        if 'RLSUF="${TBL_SUFFIX:-}"' in out and name == "RLSUF":
            problems.append("RLSUF is already defined — already patched")
            continue
        out = out.replace(anchor, repl, want)
    return out, problems


def verify(text):
    bad = []
    if 'RLSUF="${TBL_SUFFIX:-}"' not in text:
        bad.append("RLSUF is not defined")
    # every site must now carry the suffix
    for want in ('EV="rbm20_${ET}${SUF}${RLSUF}.events.tsv"',
                 'local ev="rbm20_${tag}${RLSUF}.events.tsv"',
                 'local pref="rbm20_${tag}${suf}${RLSUF}"',
                 'ev = f"rbm20_{tag}$RLSUF.events.tsv"',
                 'RLSUF="$RLSUF" python3 - $RUNS',
                 'RLSUF = os.environ.get("RLSUF", "")',
                 'parse(f"rbm20_{tag}_repro9{RLSUF}.report.txt")',
                 'parse(f"rbm20_{tag}_g45{RLSUF}.report.txt")'):
        if want not in text:
            bad.append(f"missing: {want[:58]}")
    # and no unsuffixed reader may survive
    for gone in ('EV="rbm20_${ET}${SUF}.events.tsv"',
                 'local ev="rbm20_${tag}.events.tsv"',
                 'local pref="rbm20_${tag}${suf}"',
                 'ev = f"rbm20_{tag}.events.tsv"',
                 'parse(f"rbm20_{tag}_repro9.report.txt")',
                 'parse(f"rbm20_{tag}_g45.report.txt")'):
        if gone in text:
            bad.append(f"unsuffixed site survives: {gone[:52]}")
    # The QUOTED heredoc body must not contain a literal $RLSUF: bash
    # substitutes nothing there, so it would be looked for as part of the
    # filename. Scan the body by delimiter rather than by a character offset —
    # the first version of this check used text[q+30:] and misfired.
    inq = False
    for ln in text.splitlines():
        if ln.rstrip().endswith("<<'PY'"):
            inq = True
            continue
        if inq and ln.strip() == "PY":
            inq = False
            continue
        if inq and ln.lstrip().startswith("#"):
            continue          # a python comment cannot reach a filename
        if inq and "$RLSUF" in ln:
            bad.append("a literal $RLSUF appears inside the QUOTED heredoc, "
                       "where bash substitutes nothing: " + ln.strip()[:90])
            break
    # RLSUF defined above every consumer
    i = text.find('RLSUF="${TBL_SUFFIX:-}"')
    for m in ('EV="rbm20_${ET}${SUF}${RLSUF}', 'local pref="rbm20_${tag}${suf}',
              'RLSUF="$RLSUF" python3'):
        j = text.find(m)
        if j < 0 or i > j:
            bad.append(f"RLSUF not defined above: {m[:40]}")
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
    print("127 — read-length parameter for 56, all seven sites")
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
        print("\n  Nothing written. Send me the exact lines for the failures")
        print("  above (grep -n 'rbm20_' names them all) and I will re-anchor.")
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
        print("  bash -n 56_rbm20_targets45.sbatch")
        print()
        print("  # 56 is a LOGIN-NODE-safe check: unsuffixed it should find the")
        print("  # 7 Sep tables and, because its own FORCE guard sees the 5 Oct")
        print("  # reports, [skip] every run_report rather than redo them.")
        print("  bash 56_rbm20_targets45.sbatch 2>&1 | head -30")
        print()
        print("  # the real run — ONLY after 03/03b have written")
        print("  #   rbm20_{SE,MXE}_rl125.events.tsv")
        print("  #   rbm20_{SE,MXE}_techadj_rl125.events.tsv")
        print("  TBL_SUFFIX=_rl125 sbatch 56_rbm20_targets45.sbatch")
        print()
        print("⚠ Section 4 prints a number the MANUSCRIPT quotes:")
        print('    ">>> the test is really a {N}-gene test, not a 45-gene test.')
        print('     The Methods must say {N}, not 45."')
        print("  N is however many of Gregorich's 45 have testable events. The")
        print("  tested gene set shrank at 125 bp, so N will move and the")
        print("  Methods sentence moves with it. Record the new N.")
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
    print("✓ verified on disk: all seven sites suffixed, no literal $RLSUF in")
    print("  the quoted heredoc, RLSUF above every consumer.")
    return 0


def self_test():
    import subprocess
    import tempfile

    # A fixture carrying every anchor verbatim, runnable with plain python3.
    fx = (
        "#!/bin/bash\n"
        'G45="A,B"\n'
        'NINE="A"\n'
        + A_SCRIPT
        + 'RUNS=""\n'
        + 'for ET in SE; do\n'
        + '  for SUF in "" "_techadj"; do\n'
        + '    EV="rbm20_${ET}${SUF}.events.tsv"\n'
        + '    echo "EV=$EV"\n'
        + '    RUNS="$RUNS ${ET}${SUF}"\n'
        + "  done\n"
        + "done\n"
        + 'run_report () {\n'
        + '  local tag="$1" genes="$2" suf="$3"\n'
        + '  local et="${tag%%_*}"\n'
        + '  local ev="rbm20_${tag}.events.tsv"\n'
        + '  local pref="rbm20_${tag}${suf}"\n'
        + '  echo "et=$et ev=$ev pref=$pref"\n'
        + "}\n"
        + 'for tag in $RUNS; do run_report "$tag" "$NINE" "_repro9"; done\n'
        + "python3 - $RUNS <<PY\n"
        + "import sys\n"
        + "for tag in sys.argv[1:]:\n"
        + '    ev = f"rbm20_{tag}.events.tsv"\n'
        + '    print("SEC4", ev)\n'
        + "PY\n"
        + "python3 - $RUNS <<'PY'\n"
        + "import os, re, sys\n"
        + "def parse(p):\n"
        + "    return p, None\n"
        + "for tag in sys.argv[1:]:\n"
        + '    s9, c9   = parse(f"rbm20_{tag}_repro9.report.txt")\n'
        + '    s45, c45 = parse(f"rbm20_{tag}_g45.report.txt")\n'
        + '    print("SEC5", s9, s45)\n'
        + '    print("MSG",\n'
        + '              f"read rbm20_{tag}_{{repro9,g45}}.report.txt by hand")\n'
        + '    print("MSG2",\n'
        + '              f" section 4 -- read rbm20_{tag}_g45.report.txt '
          'section 4 and compare by hand")\n'
        + "PY\n")

    new, probs = apply_edits(fx)
    assert not probs, probs
    assert not verify(new), verify(new)

    d = tempfile.mkdtemp()
    p = os.path.join(d, "fx.sh")
    open(p, "w").write(new)

    def run(env=None):
        return subprocess.run(["bash", p], capture_output=True, text=True,
                              cwd=d, env=dict(os.environ, **(env or {})))

    # 1. unsuffixed -> original names everywhere
    r = run()
    assert r.returncode == 0, r.stdout + r.stderr
    for want in ("EV=rbm20_SE.events.tsv",
                 "EV=rbm20_SE_techadj.events.tsv",
                 "et=SE ev=rbm20_SE.events.tsv pref=rbm20_SE_repro9",
                 "SEC4 rbm20_SE.events.tsv",
                 "SEC5 rbm20_SE_repro9.report.txt rbm20_SE_g45.report.txt",
                 "read rbm20_SE_{repro9,g45}.report.txt by hand"):
        assert want in r.stdout, (want, r.stdout)
    assert "RLSUF" not in r.stdout, "a literal RLSUF leaked into output"

    # 2. suffixed -> every one of the seven sites moves, suffix LAST
    r2 = run({"TBL_SUFFIX": "_rl125"})
    assert r2.returncode == 0, r2.stdout + r2.stderr
    for want in ("EV=rbm20_SE_rl125.events.tsv",
                 "EV=rbm20_SE_techadj_rl125.events.tsv",
                 "et=SE ev=rbm20_SE_rl125.events.tsv "
                 "pref=rbm20_SE_repro9_rl125",
                 "et=SE ev=rbm20_SE_techadj_rl125.events.tsv "
                 "pref=rbm20_SE_techadj_repro9_rl125",
                 "SEC4 rbm20_SE_rl125.events.tsv",
                 "SEC5 rbm20_SE_repro9_rl125.report.txt "
                 "rbm20_SE_g45_rl125.report.txt",
                 "read rbm20_SE_{repro9,g45}_rl125.report.txt by hand",
                 "read rbm20_SE_g45_rl125.report.txt section 4"):
        assert want in r2.stdout, (want, r2.stdout)
    # the tag must stay bare so et extraction still works
    assert "et=SE " in r2.stdout and "et=SE_rl125" not in r2.stdout, r2.stdout
    # no stray dollar sign survived into a filename
    assert "$RLSUF" not in r2.stdout, r2.stdout

    # 3. double application refused
    _n, p3 = apply_edits(new)
    assert p3, "second application was not refused"

    print("[ok] self-test passed: all seven sites move together with the "
          "suffix LAST; the unquoted heredoc takes $RLSUF by substitution and "
          "the quoted one by environment; the run tag stays bare so "
          "et=${tag%%_*} still resolves; unsuffixed output is byte-identical to "
          "the original names; no literal $RLSUF reaches a filename; double "
          "application refused")
    return 0


if __name__ == "__main__":
    sys.exit(main())
