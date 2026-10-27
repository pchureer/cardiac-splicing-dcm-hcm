#!/usr/bin/env python3
"""
105_patch_phase2_tally.py -- make phase2_diffsplice_v2.py print its skip tally
before it exits on "[fatal] no events tested".

WHY
---
phase2_diffsplice_v2.py line ~312:

    if len(res) == 0:
        sys.exit("[fatal] no events tested -- see skip tally above")

but the tally is printed at line ~326, i.e. AFTER that exit. So the one message
that tells you why nothing was tested is never printed, and the error instructs
you to look at output that does not exist. That cost a debugging round on
HCM vs NF rl125.

This patch is DIAGNOSTIC ONLY. It adds print statements on the zero-events
path. It cannot change any number in any .all_events.tsv, because that path
ends in sys.exit either way. Run it once; it is idempotent.

USAGE
-----
  python 105_patch_phase2_tally.py phase2_diffsplice_v2.py
  python 105_patch_phase2_tally.py phase2_diffsplice_v2.py --dry-run

A timestamped backup is written next to the file unless --no-backup.
"""
import argparse
import datetime
import shutil
import sys

MARKER = "# [105-patch] tally-before-exit"

OLD = '''    if len(res) == 0:
        sys.exit("[fatal] no events tested -- see skip tally above")
'''

NEW = '''    if len(res) == 0:
        ''' + MARKER + '''
        print("\\nevents NOT tested:", file=sys.stderr)
        for _k, _v in skips.most_common():
            print(f"  {_v:>7}  {_k}", file=sys.stderr)
        try:
            _present = sorted(str(x) for x in cov["disease"].dropna().unique())
        except Exception:
            _present = ["<no 'disease' column in covariate file>"]
        print(f"[hint] --group1={a.group1!r} --group2={a.group2!r}", file=sys.stderr)
        print(f"[hint] disease labels present in {a.cov}: {_present[:10]}",
              file=sys.stderr)
        if a.group1 not in _present or a.group2 not in _present:
            print("[hint] a requested group label is ABSENT from the covariate "
                  "file -- every event will fail --min-samples. Pass the right "
                  "--group1/--group2.", file=sys.stderr)
        sys.exit("[fatal] no events tested -- tally above")
'''


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("path")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--no-backup", action="store_true")
    a = ap.parse_args()

    src = open(a.path).read()

    if MARKER in src:
        print(f"[ok] {a.path} already carries the patch; nothing to do.")
        return 0

    n = src.count(OLD)
    if n == 0:
        print("[fatal] the expected block was not found verbatim. Do not patch "
              "blind -- show me:\n"
              f"    grep -n 'no events tested' {a.path}", file=sys.stderr)
        return 2
    if n > 1:
        print(f"[fatal] expected block appears {n} times; refusing.", file=sys.stderr)
        return 2

    out = src.replace(OLD, NEW)

    # --- gate: the patched file must still compile, and must contain the new
    #     symbols. Printing "OK" is not evidence; this is.
    try:
        compile(out, a.path, "exec")
    except SyntaxError as e:
        print(f"[fatal] patched text does not compile: {e}", file=sys.stderr)
        return 2
    for sym in (MARKER, "events NOT tested", "labels present in",
                "no events tested -- tally above"):
        if sym not in out:
            print(f"[fatal] post-patch text is missing {sym!r}", file=sys.stderr)
            return 2
    if "see skip tally above" in out:
        print("[fatal] the old message survived the replace", file=sys.stderr)
        return 2

    if a.dry_run:
        print("[dry-run] would patch. New block:\n")
        print(NEW)
        return 0

    if not a.no_backup:
        stamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        bak = f"{a.path}.bak_{stamp}"
        shutil.copy2(a.path, bak)
        print(f"[ok] backup: {bak}")

    open(a.path, "w").write(out)

    # --- verify on disk, not in memory
    chk = open(a.path).read()
    assert MARKER in chk and "see skip tally above" not in chk
    compile(chk, a.path, "exec")
    print(f"[ok] patched {a.path}; verified on disk and compiles.")
    print("[ok] diagnostic-only: the zero-events path still exits, so no "
          "result file can change.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
