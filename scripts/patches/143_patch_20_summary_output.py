#!/usr/bin/env python3
"""
143_patch_20_summary_output.py -- the write patch 142 missed.

THE MISS
--------
142 suffixed script 20's two MAGNet inputs, its five `--out` paths and the
ARMS comparison block. It did NOT suffix the summary table the heredoc writes
at the very end:

    line 244   R.to_csv("jc_vs_jcec.tsv", sep="\\t", index=False,
                        float_format="%.6g")

so a 125 bp run overwrites the 50 bp `jc_vs_jcec.tsv` in place.

142's verify checked for `--out jccheck/` and for the ARMS paths. It did not
enumerate every write inside the heredoc — which is precisely the failure
written up for script 33 on 16 Oct ("patching a producer for a suffix is not
the same as enumerating every path it writes"), committed inside the patch
meant to prevent it. Job 409724 was already running when it was found; the
50 bp file was backed up to `jc_vs_jcec_50bp_BACKUP.tsv` first.

WHAT THIS ADDS BEYOND THE ONE-LINE FIX
--------------------------------------
A verify that walks the WHOLE heredoc body and refuses if any `to_csv(` or
`open(..., "w")` writes to a quoted literal filename with no `{SUF}` in it.
A one-line fix would leave the next such write just as invisible; this makes
the class detectable instead of the instance.

USAGE
-----
  python 143_patch_20_summary_output.py [--apply] [--self-test]
"""
import argparse
import difflib
import os
import re
import shutil
import subprocess
import sys
import time

TARGET = "20_jc_vs_jcec_readlength.sbatch"

A_CSV = ('R.to_csv("jc_vs_jcec.tsv", sep="\\t", index=False, '
         'float_format="%.6g")')
B_CSV = ('_out = f"jc_vs_jcec{SUF}.tsv"\n'
         'R.to_csv(_out, sep="\\t", index=False, float_format="%.6g")\n'
         'print(f"[write] {_out}")')

EDITS = [("the summary table's filename", A_CSV, B_CSV, 1)]


def heredoc_body(text):
    """The quoted heredoc the summary runs in, or None."""
    i = text.find("<<'PY'\n")
    if i < 0:
        return None
    body = text[i + len("<<'PY'\n"):]
    j = body.find("\nPY\n")
    return body if j < 0 else body[:j]


def unsuffixed_writes(body):
    """Every write in `body` to a quoted literal filename with no {SUF}.

    This is the check 142 should have had. It looks at the heredoc as a
    whole rather than at the paths that happened to be on my mind.
    """
    out = []
    pat = re.compile(
        r'(?:\.to_csv|open)\(\s*(?:f?)(["\'])([^"\']+\.(?:tsv|txt|csv|md))\1')
    for n, line in enumerate(body.splitlines(), 1):
        s = line.strip()
        if s.startswith("#"):
            continue
        for m in pat.finditer(line):
            name = m.group(2)
            if "{SUF}" not in name:
                out.append((n, name, s[:68]))
    return out


def apply_edits(text):
    if 'f"jc_vs_jcec{SUF}.tsv"' in text:
        return text, ["the summary filename is already suffixed — 20 already "
                      "carries 143 (restore the .bak to re-apply)"]
    if "SUF=${TBL_SUFFIX:-}" not in text:
        return text, ["142 has not been applied — apply "
                      "142_patch_20_jc_readlength.py first; this fixes a "
                      "write it missed"]
    body = heredoc_body(text)
    if body is None:
        return text, ["the quoted heredoc could not be located"]
    if 'SUF = os.environ.get("TBL_SUFFIX", "")' not in body:
        return text, ["the heredoc does not define SUF — 142 should have "
                      "added it"]
    problems, out = [], text
    for name, anchor, repl, want in EDITS:
        got = out.count(anchor)
        if got != want:
            problems.append(
                f"{name}: anchor occurs {got}x, expected {want}x\n"
                f"      anchor: {anchor[:72]}")
            continue
        out = out.replace(anchor, repl, want)
    return out, problems


def verify(text):
    bad = []
    for w in ('_out = f"jc_vs_jcec{SUF}.tsv"',
              "R.to_csv(_out, sep=",
              'print(f"[write] {_out}")'):
        if w not in text:
            bad.append(f"missing: {w[:58]}")
    if 'R.to_csv("jc_vs_jcec.tsv"' in text:
        bad.append("the unsuffixed summary write survives")
    # THE GENERAL CHECK — every write in the heredoc, not just the one I knew
    body = heredoc_body(text)
    if body is None:
        bad.append("the quoted heredoc could not be located")
    else:
        for n, name, line in unsuffixed_writes(body):
            bad.append(f"heredoc line {n} still writes an unsuffixed literal "
                       f"{name!r}: {line}")
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
    print("143 — script 20's summary table, the write 142 missed")
    print("=" * 78)
    if not os.path.exists(a.target):
        print(f"[fatal] {a.target} not found.", file=sys.stderr)
        return 2
    src = open(a.target, errors="replace").read()

    # report what the general check finds BEFORE patching, so the scale of
    # the miss is visible rather than only its one known instance
    pre = heredoc_body(src)
    if pre is not None:
        found = unsuffixed_writes(pre)
        print(f"  unsuffixed literal writes inside the heredoc: {len(found)}")
        for n, name, line in found:
            print(f"    line {n:>4}  {name}")
        print()

    new, problems = apply_edits(src)
    if problems:
        print("REFUSING — " + problems[0] if len(problems) == 1
              else "REFUSING — anchors do not match:")
        for p in problems:
            print(f"  ✗ {p}")
        return 3
    bad = verify(new)
    if bad:
        print("REFUSING — patched text fails verification:")
        for b in bad:
            print(f"  ✗ {b}")
        return 3
    import tempfile
    p = os.path.join(tempfile.mkdtemp(), "chk.sbatch")
    open(p, "w").write(new)
    r = subprocess.run(["bash", "-n", p], capture_output=True, text=True)
    if r.returncode != 0:
        print("REFUSING — the patched script is not valid bash:")
        print(r.stderr)
        return 3

    print("".join(difflib.unified_diff(
        src.splitlines(True), new.splitlines(True),
        fromfile=a.target, tofile=a.target + " (patched)", n=3)))
    print("=" * 78)

    if not a.apply:
        print("DRY RUN — nothing written. Re-run with --apply.")
        print()
        print("⚠ DO THIS FIRST, once job 409724 has finished. It was running")
        print("  with the unsuffixed write, so it will have overwritten the")
        print("  50 bp table:")
        print()
        print("    mv jc_vs_jcec.tsv jc_vs_jcec_rl125.tsv")
        print("    cp jc_vs_jcec_50bp_BACKUP.tsv jc_vs_jcec.tsv")
        print("    ls -la jc_vs_jcec*")
        print()
        print("  Then both read lengths are correctly named and nothing is")
        print("  lost. After this patch, a re-run writes jc_vs_jcec_rl125.tsv")
        print("  directly and the rename is never needed again.")
        print()
        print("  bash -n 20_jc_vs_jcec_readlength.sbatch && echo BASH_OK")
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
    print("✓ verified on disk: the summary table follows SUF, and NO write")
    print("  anywhere in the heredoc targets an unsuffixed literal.")
    return 0


def self_test():
    import tempfile

    head = ("#!/bin/bash\n"
            "SUF=${TBL_SUFFIX:-}\n"
            'TBL_SUFFIX="$SUF" python3 - <<\'PY\'\n'
            "import os\n"
            "import pandas as pd\n"
            'SUF = os.environ.get("TBL_SUFFIX", "")\n'
            "R = pd.DataFrame({'a': [1]})\n")
    tail = "\nPY\n"
    fx = head + A_CSV + tail

    # the general check must SEE the problem before the patch
    found = unsuffixed_writes(heredoc_body(fx))
    assert len(found) == 1 and found[0][1] == "jc_vs_jcec.tsv", found

    new, probs = apply_edits(fx)
    assert not probs, probs
    assert not verify(new), verify(new)
    assert not unsuffixed_writes(heredoc_body(new)), \
        unsuffixed_writes(heredoc_body(new))

    d = tempfile.mkdtemp()
    p = os.path.join(d, "chk.sbatch")
    open(p, "w").write(new)
    r = subprocess.run(["bash", "-n", p], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr

    # the heredoc still compiles, and resolves the right name either way
    body = heredoc_body(new)
    compile(body, "hd", "exec")
    for suf, want in (("", "jc_vs_jcec.tsv"), ("_rl125", "jc_vs_jcec_rl125.tsv")):
        g = {"os": os}
        old = os.environ.get("TBL_SUFFIX")
        if suf:
            os.environ["TBL_SUFFIX"] = suf
        else:
            os.environ.pop("TBL_SUFFIX", None)
        try:
            src = body.replace("import pandas as pd\n", "")
            src = src.replace("R = pd.DataFrame({'a': [1]})\n", "")
            src = src.replace("R.to_csv(_out, sep=\"\\t\", index=False, "
                              "float_format=\"%.6g\")\n", "")
            exec(compile(src, "hd", "exec"), g)
        finally:
            if old is None:
                os.environ.pop("TBL_SUFFIX", None)
            else:
                os.environ["TBL_SUFFIX"] = old
        assert g["_out"] == want, (suf, g["_out"], want)

    # the general check must CATCH a second, different unsuffixed write —
    # the whole point of checking the class rather than the instance
    fx2 = (head + A_CSV + '\nopen("extra_summary.txt", "w").write("x")\n'
           + tail)
    n2, p2 = apply_edits(fx2)
    assert not p2, p2
    v2 = verify(n2)
    assert any("extra_summary.txt" in b for b in v2), v2

    # a suffixed write must NOT be flagged
    fx3 = (head + A_CSV + '\nopen(f"extra{SUF}.txt", "w").write("x")\n' + tail)
    n3, p3 = apply_edits(fx3)
    assert not p3, p3
    assert not verify(n3), verify(n3)

    # refuses without 142
    _n, p4 = apply_edits("#!/bin/bash\necho hi\n")
    assert p4 and "142 has not been applied" in p4[0], p4

    # double application refused
    n5, p5 = apply_edits(new)
    assert p5, "second application was not refused"
    assert n5 == new, "a refused re-application still altered the text"

    print("[ok] self-test passed: the unsuffixed summary write is found by a "
          "general scan of the heredoc before patching and gone after; the "
          "patched filename resolves to jc_vs_jcec.tsv and "
          "jc_vs_jcec_rl125.tsv under the two settings; the scan CATCHES a "
          "second unrelated unsuffixed write and does NOT flag a suffixed "
          "one; the script stays valid bash and the heredoc valid python; "
          "refuses without 142; double application refused")
    return 0


if __name__ == "__main__":
    sys.exit(main())
