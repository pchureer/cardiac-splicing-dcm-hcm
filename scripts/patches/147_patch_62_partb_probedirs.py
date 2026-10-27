#!/usr/bin/env python3
"""
147_patch_62_partb_probedirs.py -- (b) picks its probe directories by
CONTENT, not by name order.

WHAT BROKE, AND WHEN
--------------------
Part (b) selected directories like this:

    d125 = a.probe125 or next((d for d in dirs if "125" in d), None)

the first sorted directory name containing "125". That worked until 21 Oct,
when the 125 bp run of script 20 created `jccheck_rl125/`. It sorts ahead of
`probe_rl125/`, contains `all_events.tsv` tables rather than rMATS output,
and so (b) found no SE probe file and reported:

    125 bp -> jccheck_rl125
    SE     172 / -
    >>> (b) NOT BUILT. Skipped exons give [172, None] ...

Nothing in (b) changed. Its environment did. A NEW name quietly satisfying an
OLD glob is the mirror image of Rule 8 ("a rename is only half an edit"): here
nothing was renamed, and the reader still broke. The 50 bp side survived only
because `probe_rl50` happens to sort before `probe_rl50_tmp`.

THE FIX
-------
Selection by content, with every rejection printed:

  1. the directory must actually contain an SE probe file (the same four
     candidate names the counting loop below already tries)
  2. scratch and abandoned copies (`*tmp*`, `*abandon*`) are dropped
  3. a name beginning "probe" is preferred over one that merely matches
  4. anything still ambiguous REFUSES and asks for --probe50/--probe125,
     rather than taking the first one

An explicit --probe50/--probe125 still wins outright, as the help text says.

The chosen directories are also appended to the S5 notes, so the built table
records which probe it came from instead of leaving it to be re-derived.

WHY NOT JUST PASS --probe125
----------------------------
That fixes today's run and leaves the next new `*125*` directory to break it
again, silently, in a way that looks like a data problem rather than a
selection problem. The refusal in rule 4 is the part worth having.

USAGE
-----
  python 147_patch_62_partb_probedirs.py [--apply] [--self-test]
"""
import argparse
import difflib
import os
import shutil
import subprocess
import sys
import time

TARGET = "62_build_s5.py"

A_PICK = '''    dirs = sorted({d for d in glob.glob("*") if os.path.isdir(d)
                   and re.search(r"probe|rl\\d{2,3}|readlen", d, re.I)})
    print(f"  candidate probe directories: {', '.join(dirs) or 'NONE FOUND'}")
    # an explicit directory always wins over the guess -- the help text promised this
    d50 = a.probe50 or next((d for d in dirs
                             if re.search(r"(^|[^0-9])50([^0-9]|$)", d)), None)
    d125 = a.probe125 or next((d for d in dirs if "125" in d), None)
    print(f"  50 bp  -> {d50 or '*** not identified ***'}")
    print(f"  125 bp -> {d125 or '*** not identified ***'}")
'''

B_PICK = '''    dirs = sorted({d for d in glob.glob("*") if os.path.isdir(d)
                   and re.search(r"probe|rl\\d{2,3}|readlen", d, re.I)})
    print(f"  candidate probe directories: {', '.join(dirs) or 'NONE FOUND'}")

    # ---- DIRECTORY SELECTION, rewritten 21 Oct ------------------------
    # The old picker took the FIRST sorted name containing "125". On
    # 21 Oct the 125 bp run of script 20 created jccheck_rl125/, which
    # sorts ahead of probe_rl125/ and holds all_events tables rather than
    # rMATS output. (b) had not changed; its directory had. A NEW name
    # silently satisfying an OLD glob is the mirror of Rule 8.
    #
    # Selection is now by CONTENT, and every rejection is printed:
    #   1. the directory must contain an SE probe file
    #   2. scratch and abandoned copies are dropped
    #   3. a name beginning "probe" is preferred
    #   4. anything still ambiguous REFUSES and asks for --probeNN
    def _probe_file(_d):
        return first_existing(os.path.join(_d, "SE.MATS.JC.txt"),
                              os.path.join(_d, "SE.MATS.JCEC.txt"),
                              os.path.join(_d, "SE.all_events.tsv"),
                              os.path.join(_d, "SE_probe.all_events.tsv"))

    def _pick(_want, _flag):
        _usable, _why = [], []
        for _d in [x for x in dirs if _want(x)]:
            if _probe_file(_d) is None:
                _why.append(f"{_d}: holds no SE probe file"); continue
            if re.search(r"tmp|abandon", _d, re.I):
                _why.append(f"{_d}: scratch or abandoned copy"); continue
            _usable.append(_d)
        for _w in _why:
            print(f"    rejected {_w}")
        _pref = [x for x in _usable if x.lower().startswith("probe")]
        _chosen = _pref or _usable
        if len(_chosen) > 1:
            print(f"    AMBIGUOUS, {len(_chosen)} usable: "
                  f"{', '.join(_chosen)}")
            print(f"    pass {_flag} DIR -- refusing to guess")
            return None
        return _chosen[0] if _chosen else None

    # an explicit directory always wins over the guess -- the help text promised this
    d50 = a.probe50 or _pick(
        lambda _d: re.search(r"(^|[^0-9])50([^0-9]|$)", _d), "--probe50")
    d125 = a.probe125 or _pick(lambda _d: "125" in _d, "--probe125")
    print(f"  50 bp  -> {d50 or '*** not identified ***'}")
    print(f"  125 bp -> {d125 or '*** not identified ***'}")
'''

A_OK = '    print("\\n  [ok] skipped exons reproduce the recorded 172 / 167")\n'

B_OK = (A_OK +
        '    # provenance: record WHICH probe the built rows came from, so\n'
        '    # it never has to be re-derived from a directory listing\n'
        '    notes.append(f"(b) probe directories: 50 bp {d50}, '
        '125 bp {d125}")\n')

EDITS = [
    ("part (b)'s directory picker", A_PICK, B_PICK, 1),
    ("part (b)'s success path", A_OK, B_OK, 1),
]


def apply_edits(text):
    if "def _probe_file(_d):" in text:
        return text, ["62 already carries 147 — (b) selects by content "
                      "(restore the .bak to re-apply)"]
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
    for w in ('def _probe_file(_d):',
              'def _pick(_want, _flag):',
              'holds no SE probe file',
              'scratch or abandoned copy',
              'refusing to guess',
              'd125 = a.probe125 or _pick(lambda _d: "125" in _d, '
              '"--probe125")',
              '(b) probe directories: 50 bp {d50}, 125 bp {d125}'):
        if w not in text:
            bad.append(f"missing: {w[:60]}")
    # the first-match picker must be gone
    if 'next((d for d in dirs if "125" in d), None)' in text:
        bad.append("the first-match 125 bp picker survives")
    if 'next((d for d in dirs\n' in text and 'probe50 or next(' in text:
        bad.append("the first-match 50 bp picker survives")
    # --probeNN must still override
    if "a.probe50 or _pick(" not in text or "a.probe125 or _pick(" not in text:
        bad.append("an explicit --probeNN no longer takes precedence")
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
    print("147 — S5(b): probe directories chosen by content")
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
        print("  python3 62_build_s5.py --parts b 2>&1 | tail -30")
        print()
        print("⚠ EXPECTED: `125 bp -> probe_rl125`, a `rejected")
        print("  jccheck_rl125: holds no SE probe file` line above it, and")
        print("  SE 172 / 167. (b) then builds for the first time since")
        print("  script 20 ran this morning.")
        print()
        print("  If instead it prints AMBIGUOUS, read the names it lists and")
        print("  pass the right one with --probe125 — that refusal is the")
        print("  point of the patch, not a failure of it.")
        print()
        print("  Note (b) is a PROBE result at BOTH read lengths; it is not")
        print("  re-based by TBL_SUFFIX and 62 already prints that scope")
        print("  warning. 172 / 167 is the same pair either way.")
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
    print("✓ verified on disk: (b) requires an SE probe file, drops scratch")
    print("  and abandoned copies, prefers a probe* name, refuses when still")
    print("  ambiguous, prints every rejection, and records the chosen")
    print("  directories in the S5 notes.")
    return 0


def self_test():
    import tempfile

    head = (
        "import glob, os, re, sys\n"
        "def first_existing(*paths):\n"
        "    for p in paths:\n"
        "        if os.path.exists(p) and os.path.getsize(p):\n"
        "            return p\n"
        "    return None\n"
        "class A:\n"
        "    probe50 = os.environ.get('P50') or None\n"
        "    probe125 = os.environ.get('P125') or None\n"
        "notes = []\n"
        "a = A()\n"
        "def part_b():\n")
    fx = (head + A_PICK
          + '    print("D50", d50)\n    print("D125", d125)\n'
          + A_OK
          + '    print("NOTES", notes)\n'
          + "part_b()\n")

    new, probs = apply_edits(fx)
    assert not probs, probs
    assert not verify(new), verify(new)
    compile(new, "fx", "exec")

    d = tempfile.mkdtemp()

    def mk(name, files):
        os.makedirs(os.path.join(d, name), exist_ok=True)
        for f in files:
            open(os.path.join(d, name, f), "w").write("x\n")

    # the real layout, as `ls` showed it on the cluster
    mk("jccheck_rl125", ["Brigham_SE_HCM_JC.all_events.tsv",
                         "SE_v2filtered_JC.all_events.tsv"])
    mk("probe_rl125", ["SE.MATS.JC.txt", "SE.MATS.JCEC.txt"])
    mk("probe_rl125_tmp", ["SE.MATS.JC.txt"])
    mk("probe_rl50", ["SE.MATS.JC.txt", "SE.MATS.JCEC.txt"])
    mk("probe_rl50_tmp", ["SE.MATS.JC.txt"])
    mk("rmats_magnet_rl125", ["SE.MATS.JC.txt"])
    mk("rmats_DCMvsHCM_rl125", ["SE.MATS.JC.txt"])
    mk("rmats_magnet_rl125.abandoned_1207", ["SE.MATS.JC.txt"])

    def run(src, name, env=None):
        p = os.path.join(d, name)
        open(p, "w").write(src)
        e = dict(os.environ)
        e.pop("P50", None)
        e.pop("P125", None)
        e.update(env or {})
        r = subprocess.run([sys.executable, "-I", p], capture_output=True,
                           text=True, cwd=d, env=e)
        assert r.returncode == 0, r.stdout + r.stderr
        out = {}
        for line in r.stdout.splitlines():
            if line.startswith("D50 ") or line.startswith("D125 "):
                k, v = line.split(" ", 1)
                out[k] = v.strip()
        return r.stdout, out

    # --- the BUG must be reproduced by the unpatched fixture, or this
    #     test is not testing the thing that broke
    o0, p0 = run(fx, "fx_pre.py")
    assert p0["D125"] == "jccheck_rl125", p0
    assert p0["D50"] == "probe_rl50", p0

    # --- and fixed by the patched one
    o1, p1 = run(new, "fx_post.py")
    assert p1["D125"] == "probe_rl125", (p1, o1)
    assert p1["D50"] == "probe_rl50", (p1, o1)
    assert "rejected jccheck_rl125: holds no SE probe file" in o1, o1
    assert "rejected probe_rl125_tmp: scratch or abandoned copy" in o1, o1
    assert "rejected rmats_magnet_rl125.abandoned_1207: scratch or " \
           "abandoned copy" in o1, o1
    # rmats_magnet_rl125 is usable but NOT preferred, so it is not rejected
    # and not chosen -- the probe* preference resolves it silently
    assert "rejected rmats_magnet_rl125:" not in o1, o1
    assert "AMBIGUOUS" not in o1, o1

    # --- AMBIGUITY MUST REFUSE, not pick. Remove probe_rl125 and the two
    #     surviving rmats_*_rl125 dirs tie.
    shutil.rmtree(os.path.join(d, "probe_rl125"))
    o2, p2 = run(new, "fx_post.py")
    assert p2["D125"] == "None", (p2, o2)
    assert "AMBIGUOUS, 2 usable" in o2, o2
    assert "rmats_DCMvsHCM_rl125" in o2 and "rmats_magnet_rl125" in o2, o2
    assert "pass --probe125 DIR" in o2, o2
    assert p2["D50"] == "probe_rl50", p2       # 50 bp unaffected

    # --- an explicit --probe125 still wins, even over the refusal
    o3, p3 = run(new, "fx_post.py", {"P125": "rmats_magnet_rl125"})
    assert p3["D125"] == "rmats_magnet_rl125", (p3, o3)
    assert "AMBIGUOUS" not in o3, o3

    # --- and over a directory the content rule would have rejected, since
    #     the help text promises the explicit flag is obeyed
    o4, p4 = run(new, "fx_post.py", {"P125": "jccheck_rl125"})
    assert p4["D125"] == "jccheck_rl125", (p4, o4)

    # --- the 50 bp side is now content-checked too: gut probe_rl50 and it
    #     must fall through rather than be chosen on its name
    os.remove(os.path.join(d, "probe_rl50/SE.MATS.JC.txt"))
    os.remove(os.path.join(d, "probe_rl50/SE.MATS.JCEC.txt"))
    o5, p5 = run(new, "fx_post.py")
    assert p5["D50"] != "probe_rl50", (p5, o5)
    assert "rejected probe_rl50: holds no SE probe file" in o5, o5

    # --- double application refused without altering the text
    n6, p6 = apply_edits(new)
    assert p6 and "already carries 147" in p6[0], p6
    assert n6 == new, "a refused re-application still altered the text"

    print("[ok] self-test passed: the UNPATCHED fixture reproduces the bug "
          "(125 bp -> jccheck_rl125) and the patched one picks probe_rl125; "
          "jccheck_rl125 is rejected for holding no SE probe file, *_tmp and "
          "*.abandoned* for being scratch, and a usable-but-unpreferred "
          "rmats_* dir is resolved silently by the probe* preference; a "
          "genuine tie prints AMBIGUOUS and returns None instead of "
          "guessing; --probe125 overrides both the tie and the content rule; "
          "the 50 bp side is content-checked too; double application refused")
    return 0


if __name__ == "__main__":
    sys.exit(main())
