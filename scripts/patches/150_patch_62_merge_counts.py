#!/usr/bin/env python3
"""
150_patch_62_merge_counts.py -- S5(a) carries BOTH halves.

WHY
---
62's part (a) reports Deming slopes: JC and JCEC agree to within 0.024 on
matched events. 63_s5a_jc_vs_jcec.py measures the other half: how MANY events
each quantification calls, which differs by +9% to +58% across cohorts and
+33% within MAGNet, one-directional, never fewer under JCEC.

`RESULT_s5a_jc_vs_jcec.md` states the risk in its own words:

    "S5(a) CANNOT be summarised as 'the choice does not change the answer',
     the way (b) and (c) can. It changes the count substantially and leaves
     the direction alone, and the row has to say both."

A slopes-only panel says exactly the thing that document forbids. This patch
merges 63's rows into S5 part (a) so the two halves cannot be separated.

WHAT IT ADDS
------------
Two rows per arm under `part = "a. JC vs JCEC (counts)"`, in the SAME
seven-key schema parts (b) and (c) use, so nothing about the table's shape
changes:

    <arm>              93 significant | JC    124 significant | JCEC   +33.3%
    <arm>, direction   92/93 (98.92%) | ...   0.964 | Spearman rho      -

The slope rows keep `part = "a. JC vs JCEC"` and are untouched.

THE BASIS GUARD
---------------
After 148, a per-arm re-based run of 63 tags every arm name with `@` --
`MAGNet (discovery) @rl125`, `Colorado @as-published`. So when TBL_SUFFIX is
set, this patch REFUSES to merge an S5a table whose arm names carry no basis:
that table predates 148 and is a 50 bp measurement. Without the guard a stale
`S5a.tsv` copied to `S5a_rl125.tsv` would merge silently, which is the
defect class this project keeps hitting -- right numbers under a wrong name.

If the S5a table is simply absent, (a) still builds but records a note saying
it carries slopes without the count change, rather than quietly reporting the
reassuring half alone.

USAGE
-----
  python 150_patch_62_merge_counts.py [--apply] [--self-test]
"""
import argparse
import difflib
import os
import shutil
import subprocess
import sys
import time

TARGET = "62_build_s5.py"

A_TAIL = '''        else:
            notes.append(f"(a) {c}: NO MATCHED ARM -- one basis unreadable, "
                         f"so its `change` mixes estimator and event set")
'''

B_TAIL = A_TAIL + '''
    # ---- THE COUNT ROWS, merged 21 Oct --------------------------------
    # 63_s5a_jc_vs_jcec.py measures what the slopes above do not: how MANY
    # events each quantification calls. The two halves must travel
    # together. RESULT_s5a_jc_vs_jcec.md's own warning is that S5(a)
    # "cannot be summarised as 'the choice does not change the answer'" --
    # and a slopes-only panel says precisely that, because the matched
    # slopes agree to within 0.024 while the counts differ by 9% to 58%.
    _s5a = f"S5a{SUF}.tsv"
    if not os.path.exists(_s5a):
        print(f"\\n  ⚠ {_s5a} NOT FOUND. S5(a) will carry slopes only, which")
        print("    states only the reassuring half of the comparison. Build")
        print("    the counts with:")
        print(f"      TBL_SUFFIX={SUF} python3 63_s5a_jc_vs_jcec.py")
        notes.append(f"(a) COUNTS MISSING: {_s5a} not found, so S5(a) reports "
                     f"slopes without the count change")
    else:
        _C = pd.read_csv(_s5a, sep="\\t")
        # BASIS GUARD. After 148 a re-based run of 63 tags every arm name
        # with "@". An S5a table without that tag predates 148 and is a
        # 50 bp measurement, so merging it into a 125 bp S5 would put right
        # numbers under a wrong name -- this project's recurring defect.
        if SUF and not _C["arm"].astype(str).str.contains("@").all():
            print(f"\\n  >>> REFUSING to merge {_s5a}: its arm names carry no")
            print("      read-length basis, so it predates patch 148 and is a")
            print("      50 bp table. Re-run 63 with TBL_SUFFIX set.")
            notes.append(f"(a) COUNTS REFUSED: {_s5a} carries no per-arm "
                         f"read-length basis")
        else:
            for _, _r in _C.iterrows():
                rows.append(dict(
                    part="a. JC vs JCEC (counts)", row=str(_r["arm"]),
                    value_1=f"{_r['sig_JC']} significant", label_1="JC",
                    value_2=f"{_r['sig_JCEC']} significant", label_2="JCEC",
                    change=str(_r["sig_change"])))
                rows.append(dict(
                    part="a. JC vs JCEC (counts)",
                    row=f"{_r['arm']}, direction",
                    value_1=str(_r["dir_agree_sig"]),
                    label_1="same direction, JC-significant",
                    value_2=str(_r["rho_sig"]), label_2="Spearman rho",
                    change="-"))
            print(f"\\n  [ok] merged {len(_C)} count rows from {_s5a}")
            print("    S5(a) now states BOTH halves: the number of events")
            print("    called is strongly quantification-dependent, and the")
            print("    direction and magnitude of those events are not.")
'''

EDITS = [("part (a)'s row emit", A_TAIL, B_TAIL, 1)]


def apply_edits(text):
    if '_s5a = f"S5a{SUF}.tsv"' in text:
        return text, ["62 already carries 150 (restore the .bak to re-apply)"]
    if 'r[lab + "_m"] = deming(' not in text:
        return text, ["145 has not been applied — the matched arm must exist "
                      "before the counts are merged"]
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
    for w in ('_s5a = f"S5a{SUF}.tsv"',
              'part="a. JC vs JCEC (counts)"',
              'REFUSING to merge',
              'COUNTS MISSING',
              'COUNTS REFUSED',
              'str.contains("@").all()',
              'same direction, JC-significant'):
        if w not in text:
            bad.append(f"missing: {w[:58]}")
    # the slope rows must keep their own part label
    if 'part="a. JC vs JCEC", row=c,' not in text:
        bad.append("the slope rows' part label changed — they and the count "
                   "rows must stay distinguishable in S5.tsv")
    # the schema must not have grown
    for k in ("value_3", "label_3"):
        if k in text:
            bad.append(f"the row schema grew {k} — parts (b) and (c) write "
                       f"the same seven keys")
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
    print("150 — S5(a) carries the count rows as well as the slopes")
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
        print("  python3 62_build_s5.py 2>&1 | tail -45")
        print()
        print("  TBL_SUFFIX=_rl125 \\")
        print("  S5_SLOPE=Colorado=0.5999,Brigham=0.5099,"
              "Stanford=0.4008,Fuwai=0.2948 \\")
        print("  S5_CONC=Colorado=80.8,Brigham=81.9,Fuwai=70.8,Stanford=75.8 \\")
        print("  S5_DISC_SE=exon=813,full=944 \\")
        print("  python3 62_build_s5.py 2>&1 | tail -45")
        print()
        print("⚠ EXPECTED: `merged 5 count rows` at each read length, and ten")
        print("  new rows under `a. JC vs JCEC (counts)` beside the six slope")
        print("  rows. The slope rows must not move.")
        print()
        print("⚠ At 125 bp the arms must read `@rl125` / `@as-published`. If")
        print("  62 REFUSES the merge, S5a_rl125.tsv predates 148 — re-run")
        print("  `TBL_SUFFIX=_rl125 python3 63_s5a_jc_vs_jcec.py` first. That")
        print("  refusal is the guard working, not a failure.")
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
    print("✓ verified on disk: S5(a) merges 63's count rows under their own")
    print("  part label, refuses an S5a table with no per-arm read-length")
    print("  basis, records a note if the counts are missing, and leaves the")
    print("  slope rows and the seven-key schema alone.")
    return 0


def self_test():
    import tempfile

    fx = ("import os\n"
          "import numpy as np\n"
          "import pandas as pd\n"
          'SUF = os.environ.get("TBL_SUFFIX", "")\n'
          "# the 145 marker apply_edits requires\n"
          "def _unused(r, lab, mm):\n"
          '    r[lab + "_m"] = deming(mm.a, mm.b, mm.c, mm.d)\n'
          "def emit(res, rows, notes):\n"
          "    for c in sorted(res):\n"
          '        rows.append(dict(part="a. JC vs JCEC", row=c,\n'
          '                         value_1=res[c]["JCEC"], label_1="Deming, JCEC",\n'
          '                         value_2=res[c]["JC"], label_2="Deming, JC",\n'
          '                         change=res[c]["JC"] - res[c]["JCEC"]))\n'
          "        if True:\n"
          "            pass\n"
          + A_TAIL +
          "\n"
          'res = {"Brigham": {"JCEC": 0.476, "JC": 0.581}}\n'
          "rows, notes = [], []\n"
          "emit(res, rows, notes)\n"
          'print("NROWS", len(rows))\n'
          "for _r in rows:\n"
          '    print("ROW|", _r["part"], "|", _r["row"], "|", _r["value_1"],\n'
          '          "|", _r["label_1"], "|", _r["value_2"], "|",\n'
          '          _r["label_2"], "|", _r["change"])\n'
          'print("NOTES", notes)\n')

    new, probs = apply_edits(fx)
    assert not probs, probs
    assert not verify(new), verify(new)
    compile(new, "fx", "exec")

    d = tempfile.mkdtemp()

    import pandas as pd

    def s5a(path, arms):
        pd.DataFrame([{
            "part": "a. JC vs JCEC", "event": "SE", "arm": a_,
            "sig_JC": 93, "sig_JCEC": 124, "sig_change": "+33.3%",
            "tested_JC": 172662, "tested_JCEC": 174445,
            "tested_change": "+1.0%", "dir_agree": "x (96.17%)", "rho": 0.892,
            "dir_agree_sig": "92/93 (98.92%)", "rho_sig": 0.964,
        } for a_ in arms]).to_csv(os.path.join(d, path), sep="\t", index=False)

    def run(src, name, env=None):
        p = os.path.join(d, name)
        open(p, "w", encoding="utf-8").write(src)
        e = dict(os.environ)
        e.pop("TBL_SUFFIX", None)
        e.update(env or {})
        r = subprocess.run([sys.executable, "-I", p], capture_output=True,
                           text=True, cwd=d, env=e)
        assert r.returncode == 0, r.stdout + r.stderr
        return r.stdout

    # --- counts MISSING: (a) still builds, but records the shortfall
    o0 = run(new, "fx.py")
    assert "S5a.tsv NOT FOUND" in o0, o0
    assert "COUNTS MISSING" in o0, o0
    assert "NROWS 1" in o0, o0

    # --- counts present at 50 bp: merged, two rows per arm
    s5a("S5a.tsv", ["MAGNet (discovery)", "Colorado"])
    o1 = run(new, "fx.py")
    assert "merged 2 count rows from S5a.tsv" in o1, o1
    assert "NROWS 5" in o1, o1                       # 1 slope + 2 arms x 2
    assert "a. JC vs JCEC (counts) | MAGNet (discovery) | 93 significant" in o1, o1
    assert "| JC | 124 significant | JCEC | +33.3%" in o1, o1
    assert "MAGNet (discovery), direction | 92/93 (98.92%)" in o1, o1
    assert "Spearman rho" in o1, o1
    # the slope row is untouched and still distinguishable
    assert "ROW| a. JC vs JCEC | Brigham | 0.476 | Deming, JCEC" in o1, o1

    # --- 125 bp with an UNTAGGED S5a: must REFUSE, not merge
    s5a("S5a_rl125.tsv", ["MAGNet (discovery)", "Colorado"])
    o2 = run(new, "fx.py", {"TBL_SUFFIX": "_rl125"})
    assert "REFUSING to merge S5a_rl125.tsv" in o2, o2
    assert "predates patch 148" in o2, o2
    assert "COUNTS REFUSED" in o2, o2
    assert "NROWS 1" in o2, o2                       # nothing merged

    # --- 125 bp with a properly tagged S5a: merged
    s5a("S5a_rl125.tsv", ["MAGNet (discovery) @rl125", "Colorado @as-published"])
    o3 = run(new, "fx.py", {"TBL_SUFFIX": "_rl125"})
    assert "merged 2 count rows from S5a_rl125.tsv" in o3, o3
    assert "REFUSING" not in o3, o3
    assert "NROWS 5" in o3, o3
    assert "@as-published" in o3, o3

    # --- a PARTIALLY tagged table is also refused: one untagged arm is
    #     enough to make the basis of the whole table unknown
    s5a("S5a_rl125.tsv", ["MAGNet (discovery) @rl125", "Colorado"])
    o4 = run(new, "fx.py", {"TBL_SUFFIX": "_rl125"})
    assert "REFUSING to merge" in o4, o4

    # --- double application refused without altering the text
    n5, p5 = apply_edits(new)
    assert p5 and "already carries 150" in p5[0], p5
    assert n5 == new, "a refused re-application still altered the text"

    # --- refuses without 145
    _n, p6 = apply_edits("print('bare')\n")
    assert p6 and "145 has not been applied" in p6[0], p6

    print("[ok] self-test passed: the count rows merge two per arm under "
          "their own part label in the unchanged seven-key schema, the slope "
          "rows stay put and stay distinguishable; a missing S5a table "
          "records COUNTS MISSING instead of silently reporting slopes "
          "alone; with TBL_SUFFIX set an S5a whose arm names carry no "
          "read-length basis is REFUSED, including when only one arm is "
          "untagged, and a properly tagged one merges; refuses without 145; "
          "double application refused")
    return 0


if __name__ == "__main__":
    sys.exit(main())
