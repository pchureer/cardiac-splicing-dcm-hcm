#!/usr/bin/env python3
"""
124_patch_14_provenance_fix.py -- fix two defects I introduced in 123.

Applies to `14_run_gprofiler.py` AFTER 123 has been applied. Anchors are lines
123 itself inserted, so this is a correction, not a reapplication.

DEFECT 1 — the provenance file was overwritten by the artefact run
------------------------------------------------------------------
14 names its results `<tag>_gprofiler.tsv` normally and
`<tag>_gprofiler_WHOLEGENOME.tsv` under `--whole-genome`. 123 named the
provenance file `<tag>_gprofiler_PROVENANCE.txt` in BOTH cases, so:

    python3 14_run_gprofiler.py --tag compadj_rl125                 # real run
    python3 14_run_gprofiler.py --tag compadj_rl125 --whole-genome  # artefact

left one provenance file on disk describing the SECOND query. The deposited
provenance then belongs to the run nobody reports, and nothing says so. This is
the defect class 123 existed to remove, reintroduced by 123.

DEFECT 2 — a false Methods sentence for the whole-genome run
------------------------------------------------------------
In whole-genome mode the payload carries no `background` key, so

    _nbg = len(payload.get("background") or [])

gave 0 and the draft read:

    "... a custom background of 0 genes carrying at least one event that
     passed the coverage and variance filters."

False twice: there is no custom background in that mode, and 0 is an artefact
of the absent key. A script must not hand over prose it cannot support — which
is exactly why 123 withholds the sentence when no version is captured. The same
rule was not applied to the whole-genome path.

THE FIX
-------
  * the provenance filename mirrors the results filename, so the two runs
    deposit two files
  * the file records `whole_genome` and the background mode explicitly
  * under `--whole-genome` NO Methods sentence is offered; a block explains that
    the run is the background-sensitivity check and must be reported as such

USAGE
-----
  python 124_patch_14_provenance_fix.py
  python 124_patch_14_provenance_fix.py --apply
  python 124_patch_14_provenance_fix.py --self-test

AFTER APPLYING, re-run the proper-background query to regenerate the provenance
the artefact run destroyed:

  python3 14_run_gprofiler.py --tag compadj_rl125
  python3 14_run_gprofiler.py --tag compadj_50check
"""
import argparse
import difflib
import os
import shutil
import sys
import time

TARGET = "14_run_gprofiler.py"

A1 = ('    import time as _time\n'
      '    _prov = f"{a.tag}_gprofiler_PROVENANCE.txt"\n')
B1 = ('    import time as _time\n'
      '    # The results file is <tag>_gprofiler_WHOLEGENOME.tsv in whole-genome\n'
      '    # mode, so the provenance file must carry the same distinction.\n'
      '    # Without it the artefact run OVERWRITES the real run\'s provenance\n'
      '    # and the deposited file describes a query nobody will report.\n'
      '    _wg = bool(getattr(a, "whole_genome", False))\n'
      '    _prov = (f"{a.tag}_gprofiler"\n'
      '             f"{\'_WHOLEGENOME\' if _wg else \'\'}_PROVENANCE.txt")\n')

A2 = '        _pf.write(f"sources requested  {a.sources}\\n")\n'
B2 = ('        _pf.write(f"sources requested  {a.sources}\\n")\n'
      '        _pf.write(f"whole_genome       {_wg}\\n")\n'
      '        _pf.write("background mode    "\n'
      '                  + ("WHOLE GENOME — no custom background; this run is "\n'
      '                     "the artefact\\n                   demonstration, "\n'
      '                     "NOT the reported result\\n"\n'
      '                     if _wg else\n'
      '                     "custom (genes with >=1 tested event)\\n"))\n')

A3 = ('    if _ver:\n'
      '        _gv = _ver.get("gprofiler_version", "?")\n')
B3 = ('    if _ver and _wg:\n'
      '        print("", file=sys.stderr)\n'
      '        print("  ** WHOLE-GENOME RUN — the ARTEFACT DEMONSTRATION, not "\n'
      '              "your result.", file=sys.stderr)\n'
      '        print("     No Methods sentence is offered: there is NO custom "\n'
      '              "background here,", file=sys.stderr)\n'
      '        print("     so describing one would be false. Report this run as "\n'
      '              "the", file=sys.stderr)\n'
      '        print("     background-sensitivity check, and cite the version "\n'
      '              "from the", file=sys.stderr)\n'
      '        print("     proper-background run.", file=sys.stderr)\n'
      '        print("", file=sys.stderr)\n'
      '    elif _ver:\n'
      '        _gv = _ver.get("gprofiler_version", "?")\n')

EDITS = [
    ("the provenance filename, made whole-genome aware", A1, B1, 1),
    ("the background mode, recorded in the file", A2, B2, 1),
    ("the draft sentence, withheld for the whole-genome run", A3, B3, 1),
]


def apply_edits(text):
    problems, out = [], text
    for name, anchor, repl, want in EDITS:
        got = out.count(anchor)
        if got != want:
            problems.append(f"{name}: anchor occurs {got}x, expected {want}x")
            continue
        if "_wg = bool(getattr" in out and name.startswith("the provenance"):
            problems.append("already fixed (_wg is defined)")
            continue
        out = out.replace(anchor, repl, want)
    return out, problems


def verify(text):
    bad = []
    code = "\n".join(l for l in text.splitlines()
                     if not l.lstrip().startswith("#"))
    if '_wg = bool(getattr(a, "whole_genome", False))' not in code:
        bad.append("_wg is not defined")
    if "_WHOLEGENOME' if _wg else" not in code:
        bad.append("the provenance filename is not whole-genome aware")
    if 'f"{a.tag}_gprofiler_PROVENANCE.txt"' in code:
        bad.append("the single-name provenance filename survives")
    if "if _ver and _wg:" not in code:
        bad.append("the whole-genome branch is missing")
    if "elif _ver:" not in code:
        bad.append("the normal draft branch was not preserved")
    if "whole_genome       {_wg}" not in code:
        bad.append("whole_genome is not recorded in the file")
    i_def = text.find('_wg = bool(getattr')
    for marker in ("_WHOLEGENOME' if _wg else", "if _ver and _wg:",
                   "whole_genome       {_wg}"):
        j = text.find(marker)
        if j < 0 or i_def > j:
            bad.append(f"_wg is not defined above its use: {marker[:28]}")
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
    print("124 — fix the provenance clobbering and the false draft sentence")
    print("=" * 78)
    if not os.path.exists(a.target):
        print(f"[fatal] {a.target} not found.", file=sys.stderr)
        return 2
    src = open(a.target, errors="replace").read()
    if "_gprofiler_PROVENANCE" not in src:
        print("[fatal] 123 has not been applied to this file; there is nothing",
              "to fix.", file=sys.stderr)
        return 2

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
        fromfile=a.target, tofile=a.target + " (fixed)", n=2)))
    print("=" * 78)
    if not a.apply:
        print("DRY RUN — nothing written. Re-run with --apply, then:")
        print()
        print("  # regenerate the provenance the artefact run destroyed")
        print("  python3 14_run_gprofiler.py --tag compadj_rl125")
        print("  python3 14_run_gprofiler.py --tag compadj_50check")
        print()
        print("  # the artefact run now writes its OWN provenance file")
        print("  python3 14_run_gprofiler.py --tag compadj_rl125 --whole-genome")
        print()
        print("  ls -l *_PROVENANCE.txt   # expect three distinct files")
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
    import json as _json
    import subprocess
    import tempfile

    # a fixture carrying 123's inserted block verbatim
    block = (A1
             + '    _vurl = "https://x/api/util/data_versions/"\n'
             + '    _ver = CANNED\n'
             + '    with open(_prov, "w") as _pf:\n'
             + '        _pf.write("VERSION CAPTURED\\n")\n'
             + A2
             + '        _pf.write(json.dumps(_ver, sort_keys=True) + "\\n")\n'
             + '    print(f"[prov] wrote {_prov}", file=sys.stderr)\n'
             + A3
             + '        _nbg = len(payload.get("background") or [])\n'
             + '        print(f"DRAFT bg={_nbg}", file=sys.stderr)\n'
             + '    else:\n'
             + '        print("NO VERSION", file=sys.stderr)\n')

    head = ("import json\n"
            "import sys\n"
            "class A:\n"
            '    tag = "t"\n'
            '    organism = "hsapiens"\n'
            '    sources = "GO:BP"\n'
            "    whole_genome = len(sys.argv) > 1\n"
            "a = A()\n"
            'payload = {"query": ["A", "B"]}\n'
            "if not a.whole_genome:\n"
            '    payload["background"] = ["A", "B", "C"]\n'
            "CANNED = " + _json.dumps({"gprofiler_version": "e114"}) + "\n"
            "def run():\n")
    tail = "\nif __name__ == '__main__':\n    run()\n"

    before = head + block + tail
    new, probs = apply_edits(before)
    assert not probs, probs
    assert not verify(new), verify(new)

    d = tempfile.mkdtemp()
    p = os.path.join(d, "fx.py")
    open(p, "w").write(new)

    # proper-background run
    r = subprocess.run([sys.executable, "-I", p], capture_output=True,
                       text=True, cwd=d)
    assert r.returncode == 0, r.stdout + r.stderr
    assert os.path.exists(os.path.join(d, "t_gprofiler_PROVENANCE.txt")), \
        os.listdir(d)
    assert "DRAFT bg=3" in r.stderr, r.stderr
    t1 = open(os.path.join(d, "t_gprofiler_PROVENANCE.txt")).read()
    assert "whole_genome       False" in t1, t1
    assert "custom (genes with >=1 tested event)" in t1, t1

    # whole-genome run must write a DIFFERENT file and offer NO draft
    r2 = subprocess.run([sys.executable, "-I", p, "wg"], capture_output=True,
                        text=True, cwd=d)
    assert r2.returncode == 0, r2.stdout + r2.stderr
    wgp = os.path.join(d, "t_gprofiler_WHOLEGENOME_PROVENANCE.txt")
    assert os.path.exists(wgp), os.listdir(d)
    assert "DRAFT" not in r2.stderr, r2.stderr
    assert "ARTEFACT DEMONSTRATION" in r2.stderr, r2.stderr
    t2 = open(wgp).read()
    assert "whole_genome       True" in t2, t2
    assert "WHOLE GENOME" in t2, t2
    # and the first file must be INTACT
    assert open(os.path.join(d, "t_gprofiler_PROVENANCE.txt")).read() == t1, \
        "the artefact run overwrote the real run's provenance again"

    _n, p3 = apply_edits(new)
    assert p3, "second application was not refused"

    print("[ok] self-test passed: the two runs write two distinct provenance "
          "files and the real one survives the artefact run; whole_genome and "
          "the background mode are recorded; the whole-genome run offers NO "
          "draft sentence and says why; the normal draft still fires with the "
          "real background count; double application refused")
    return 0


if __name__ == "__main__":
    sys.exit(main())
