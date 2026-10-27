#!/usr/bin/env python3
"""
156_patch_90_panelD_and_provenance.py -- what 155 missed in Figure 5.

THE MISS
--------
155 parameterised 90's GATES (`REC_SCAN`, `REC_FLOOR_RANK`, `REC_FWD`) and its
input paths. The gates then matched at both read lengths — G4 returned
`SE -0.051 (-0.159, +0.047) record -0.051 ... match`. So the figure built, and
`legend_check/f5_rl125.png` looks finished.

It is not. Two things in the same file are neither gates nor inputs:

1. **Panel D is drawn from HARD-CODED CONSTANTS.** Line 46 of 90 says so
   outright: "panel D   constants below, read off rbm20_g45_406229.out step 5"
   — a **5 Oct, 50 bp** log. `TECH45` and `TECH9` are literals, so panel D
   printed byte-identical numbers in both runs and **the 125 bp figure carries
   a 50 bp panel**.
2. **The provenance block prints two 50 bp values**, under a heading that reads
   "PROVENANCE — paste into the legend":
       40 of 45 genes testable for SE (663 events), 39 for MXE (1,393)
       panel D  45-gene background medians (rbm20_g45_406229.out step 5)
   At 125 bp the event counts are **641** and **1,399**, and the log is
   `rbm20_g45_407951.out`.

`RULE_enumerate_every_write_path` says a missing output suffix "destroys the
comparator and exits 0". Here it exits 0 having **drawn a figure that mixes read
lengths** and printed a legend naming the wrong numbers. Same failure, different
coat: **a `print` is a write path, and so is a plotted constant.** This is the
second time in a week — 138 enumerated its gates and left three printed "203"s
for 139 to clean up.

THE VALUES, AND HOW THEY WERE VERIFIED
--------------------------------------
Panel D plots the **45-gene background median** (the `bg` of the `45: tgt/bg`
column) for DCM/nocell and DCM/cell, taken from the unadjusted block as
"before" and the techadj block as "after". That mapping was confirmed by
deriving **all eight** existing 50 bp constants from `rbm20_g45_406229.out`
section 5 before transcribing anything from the 125 bp log:

| | 50 bp (verified) | 125 bp |
|---|---|---|
| SE without comp. | 0.95 -> 0.86 | **0.96 -> 0.87** |
| SE with comp. | 1.02 -> 0.89 | **1.02 -> 0.92** |
| MXE without comp. | 2.03 -> 1.23 | **1.98 -> 1.24** |
| MXE with comp. | 2.71 -> 1.63 | **2.65 -> 1.63** |
| SE nine, without | 1.01 -> 0.93 | **1.03 -> 0.95** |
| SE nine, with | 1.09 -> 0.96 | **1.08 -> 1.00** |
| MXE nine, without | 2.01 -> 1.28 | **1.95 -> 1.29** |
| MXE nine, with | 2.65 -> 1.65 | **2.52 -> 1.67** |

Source: `rbm20_g45_407951.out` section 5 (job 407951, 16 Oct).

`NULL_EXPECT = 0.67` is **not** parameterised — it is a theoretical
expectation, not a measurement.

USAGE
-----
  python 156_patch_90_panelD_and_provenance.py [--apply] [--self-test]

  # 125 bp, all seven variables (155's three plus 156's four)
  F5_SCAN=signature_rbp_margins_rl125.tsv \\
  F5_FWD=rbm20_margin_resolve_rl125.tsv \\
  F5_REC_SCAN="SE=PTBP1/0.062/-0.025/0.122,MXE=PTBP1/0.113/0.016/0.187" \\
  F5_FLOOR_RANK="SE=6,MXE=3" \\
  F5_REC_FWD="SE=-0.051/-0.159/0.047,MXE=-0.019/-0.108/0.072" \\
  F5_TECH45="SE=0.96/0.87/1.02/0.92,MXE=1.98/1.24/2.65/1.63" \\
  F5_TECH9="SE=1.03/0.95/1.08/1.00,MXE=1.95/1.29/2.52/1.67" \\
  F5_NTARGET="SE=641,MXE=1399" \\
  F5_D_LOG=rbm20_g45_407951.out \\
  python3 90_build_f5.py --model DCM/cell \\
      --matching "3-D match (+|disease beta|)" --out legend_check/f5_rl125
"""
import argparse
import difflib
import os
import shutil
import subprocess
import sys
import time

TARGET = "90_build_f5.py"

A_TECH = '''TECH45 = {   # event type -> {model label: (before, after)}
    "SE":  {"without composition adj.": (0.95, 0.86),
            "with composition adj.":    (1.02, 0.89)},
    "MXE": {"without composition adj.": (2.03, 1.23),
            "with composition adj.":    (2.71, 1.63)},
}
TECH9 = {"SE": {"without": (1.01, 0.93), "with": (1.09, 0.96)},
         "MXE": {"without": (2.01, 1.28), "with": (2.65, 1.65)}}
'''

B_TECH = '''# READ LENGTH, added 22 Oct by 156. These are the ONLY numbers in this file
# that are PLOTTED from constants rather than read from a deposit, which is why
# 155 missed them and why the 125 bp figure carried a 50 bp panel D while every
# gate passed. A plotted constant is a write path.
#
# Each entry is without_before/without_after/with_before/with_after, i.e. the
# 45-gene (or nine-gene) BACKGROUND median for DCM/nocell and DCM/cell, taken
# from the unadjusted block as `before` and the techadj block as `after` of
# section 5 of the rbm20_g45 log. All eight 50 bp values below were re-derived
# from rbm20_g45_406229.out before the 125 bp ones were transcribed.
#
#   F5_TECH45="SE=0.96/0.87/1.02/0.92,MXE=1.98/1.24/2.65/1.63"
#   F5_TECH9="SE=1.03/0.95/1.08/1.00,MXE=1.95/1.29/2.52/1.67"
_T45 = _f5_env("F5_TECH45")
TECH45 = ({_et: {"without composition adj.": (float(_v[0]), float(_v[1])),
                 "with composition adj.":    (float(_v[2]), float(_v[3]))}
           for _et, _v in _T45.items()} if _T45 else
          {"SE":  {"without composition adj.": (0.95, 0.86),
                   "with composition adj.":    (1.02, 0.89)},
           "MXE": {"without composition adj.": (2.03, 1.23),
                   "with composition adj.":    (2.71, 1.63)}})

_T9 = _f5_env("F5_TECH9")
TECH9 = ({_et: {"without": (float(_v[0]), float(_v[1])),
                "with":    (float(_v[2]), float(_v[3]))}
          for _et, _v in _T9.items()} if _T9 else
         {"SE": {"without": (1.01, 0.93), "with": (1.09, 0.96)},
          "MXE": {"without": (2.01, 1.28), "with": (2.65, 1.65)}})

# The provenance block prints these under "paste into the legend", so they are
# as quotable as anything in the figure.
#
# INTEGERS, NOT STRINGS, and the thousands separator is added at format time.
# A first version of this patch took them as strings so "1,393" could keep its
# comma -- but F5_NTARGET is itself a COMMA-separated list, so "MXE=1,399"
# parsed as MXE=1 and a stray "399". The separator cannot appear in a value.
#   F5_NTARGET="SE=641,MXE=1399"      F5_D_LOG=rbm20_g45_407951.out
_NTE = _f5_env("F5_NTARGET")
_NT_SE = int(_NTE["SE"][0]) if _NTE and "SE" in _NTE else 663
_NT_MXE = int(_NTE["MXE"][0]) if _NTE and "MXE" in _NTE else 1393
_D_LOG = os.environ.get("F5_D_LOG", "").strip() or "rbm20_g45_406229.out"
'''

# 155 wrote this guard; 156 extends the list it checks.
A_GUARD = '''    _miss = [_n for _n in ("F5_REC_SCAN", "F5_FLOOR_RANK", "F5_REC_FWD")
             if not os.environ.get(_n, "").strip()]'''
B_GUARD = '''    _miss = [_n for _n in ("F5_REC_SCAN", "F5_FLOOR_RANK", "F5_REC_FWD",
                           "F5_TECH45", "F5_TECH9", "F5_NTARGET", "F5_D_LOG")
             if not os.environ.get(_n, "").strip()]'''

A_MSG = '''            "  wrong gene rather than fail. Supply all three."'''
B_MSG = '''            "  wrong gene rather than fail. And panel D is PLOTTED from\\n"
            "  constants (TECH45/TECH9) that no gate checks, so without\\n"
            "  F5_TECH45/F5_TECH9 the figure would mix read lengths while\\n"
            "  every gate passed. Supply all seven."'''

A_PROV_N = 'SE (663 events), 39 for MXE (1,393)'
B_PROV_N = 'SE ({_NT_SE:,} events), 39 for MXE ({_NT_MXE:,})'

A_PROV_L = '(rbm20_g45_406229.out step 5)'
B_PROV_L = '({_D_LOG} step 5)'

EDITS = [
    ("panel D's plotted constants", A_TECH, B_TECH, 1),
    ("155's guard list", A_GUARD, B_GUARD, 1),
    ("155's refusal message", A_MSG, B_MSG, 1),
    ("the provenance event counts", A_PROV_N, B_PROV_N, 1),
    ("the provenance log reference", A_PROV_L, B_PROV_L, 1),
]


def apply_edits(text):
    if '_T45 = _f5_env("F5_TECH45")' in text:
        return text, ["90 already carries 156 (restore the .bak to re-apply)"]
    if "def _f5_env(" not in text:
        return text, ["155 has not been applied — apply "
                      "155_patch_65_90_readlength_records.py first; 156 "
                      "reuses its _f5_env() and extends its guard"]
    problems, out = [], text
    for name, anchor, repl, want in EDITS:
        got = out.count(anchor)
        if got != want:
            problems.append(f"{name}: anchor occurs {got}x, expected {want}x\n"
                            f"      first line: {anchor.splitlines()[0][:66]}")
            continue
        out = out.replace(anchor, repl, want)
    return out, problems


def verify(text):
    bad = []
    for w in ('_T45 = _f5_env("F5_TECH45")', '_T9 = _f5_env("F5_TECH9")',
              '_NT_SE = int(_NTE["SE"][0])',
              '_D_LOG = os.environ.get("F5_D_LOG"',
              '"F5_TECH45", "F5_TECH9", "F5_NTARGET", "F5_D_LOG"',
              'A plotted constant is a write path',
              'SE ({_NT_SE:,} events), 39 for MXE ({_NT_MXE:,})',
              '({_D_LOG} step 5)',
              'Supply all seven'):
        if w not in text:
            bad.append(f"missing: {w[:56]}")
    # every 50 bp default must survive verbatim
    for w in ('"without composition adj.": (0.95, 0.86)',
              '"with composition adj.":    (1.02, 0.89)',
              '"without composition adj.": (2.03, 1.23)',
              '"with composition adj.":    (2.71, 1.63)',
              '{"without": (1.01, 0.93), "with": (1.09, 0.96)}',
              '{"without": (2.01, 1.28), "with": (2.65, 1.65)}',
              'else 663', 'else 1393', '"rbm20_g45_406229.out"'):
        if w not in text:
            bad.append(f"a 50 bp default was lost: {w[:56]}")
    # the stale literals must be GONE from the provenance
    if 'SE (663 events)' in text:
        bad.append("the provenance still hard-codes 663")
    if '(rbm20_g45_406229.out step 5)' in text:
        bad.append("the provenance still hard-codes the 50 bp log")
    # NULL_EXPECT is a theoretical value, not a measurement
    if "NULL_EXPECT = 0.67" not in text:
        bad.append("NULL_EXPECT was altered; it is a theoretical expectation "
                   "and must stay a literal")
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
    print("156 — Figure 5's panel D and provenance block")
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
    import tempfile
    p = os.path.join(tempfile.mkdtemp(), "chk.py")
    open(p, "w").write(new)
    r = subprocess.run([sys.executable, "-c",
                        f"import ast;ast.parse(open({p!r}).read())"],
                       capture_output=True, text=True)
    if r.returncode != 0:
        print("REFUSING — not valid python:\n" + r.stderr)
        return 3

    print("".join(difflib.unified_diff(
        src.splitlines(True), new.splitlines(True),
        fromfile=a.target, tofile=a.target + " (patched)", n=2)))
    print("=" * 78)

    if not a.apply:
        print("DRY RUN — nothing written. Re-run with --apply, then:")
        print()
        print('  MATCH="3-D match (+|disease beta|)"')
        print()
        print("  # 1. 50 bp must be unchanged — panel D still 0.95/0.86 etc.")
        print("  python3 90_build_f5.py --model DCM/cell --matching \"$MATCH\" \\")
        print("      --out legend_check/f5_56check 2>&1 | tail -22")
        print()
        print("  # 2. the refusal: suffixed input, panel D not supplied")
        print("  F5_SCAN=signature_rbp_margins_rl125.tsv \\")
        print("  F5_FWD=rbm20_margin_resolve_rl125.tsv \\")
        print('  F5_REC_SCAN="SE=PTBP1/0.062/-0.025/0.122,'
              'MXE=PTBP1/0.113/0.016/0.187" \\')
        print('  F5_FLOOR_RANK="SE=6,MXE=3" \\')
        print('  F5_REC_FWD="SE=-0.051/-0.159/0.047,MXE=-0.019/-0.108/0.072" \\')
        print("  python3 90_build_f5.py --model DCM/cell --matching \"$MATCH\" \\")
        print("      --out /tmp/x    # must REFUSE, naming F5_TECH45 etc.")
        print()
        print("  # 3. the real 125 bp figure, all seven variables")
        print("  F5_SCAN=signature_rbp_margins_rl125.tsv \\")
        print("  F5_FWD=rbm20_margin_resolve_rl125.tsv \\")
        print('  F5_REC_SCAN="SE=PTBP1/0.062/-0.025/0.122,'
              'MXE=PTBP1/0.113/0.016/0.187" \\')
        print('  F5_FLOOR_RANK="SE=6,MXE=3" \\')
        print('  F5_REC_FWD="SE=-0.051/-0.159/0.047,MXE=-0.019/-0.108/0.072" \\')
        print('  F5_TECH45="SE=0.96/0.87/1.02/0.92,MXE=1.98/1.24/2.65/1.63" \\')
        print('  F5_TECH9="SE=1.03/0.95/1.08/1.00,MXE=1.95/1.29/2.52/1.67" \\')
        print('  F5_NTARGET="SE=641,MXE=1399" \\')
        print("  F5_D_LOG=rbm20_g45_407951.out \\")
        print("  python3 90_build_f5.py --model DCM/cell --matching \"$MATCH\" \\")
        print("      --out legend_check/f5_rl125 2>&1 | tail -22")
        print()
        print("⚠ DELETE legend_check/f5_rl125.png AND .pdf FIRST. The copies on")
        print("  disk have a 50 bp panel D and a legend naming 663/1,393 and")
        print("  the 5 Oct log. Every gate passed on them, so nothing in the")
        print("  output says they are wrong.")
        print()
        print("⚠ Check panel D's printed numbers CHANGE between runs 1 and 3.")
        print("  If they do not, TECH45 is still being read from the literals")
        print("  and the figure is still mixed.")
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
    print("✓ verified on disk: panel D's constants and the provenance block")
    print("  follow the read length, the 50 bp values survive as defaults,")
    print("  NULL_EXPECT is untouched, and 155's guard now refuses a suffixed")
    print("  run that has not supplied panel D.")
    return 0


def self_test():
    import tempfile

    # a fixture in the POST-155 state: it must already carry _f5_env and the
    # guard, because 156 extends both.
    fx = ("import os\nimport sys\n"
          'SCAN = os.environ.get("F5_SCAN", "signature_rbp_margins.tsv")\n'
          'FWD = os.environ.get("F5_FWD", "rbm20_margin_resolve.tsv")\n'
          "def _f5_env(_name):\n"
          '    _raw = os.environ.get(_name, "").strip()\n'
          "    if not _raw:\n"
          "        return None\n"
          "    _out = {}\n"
          '    for _part in _raw.split(","):\n'
          '        _k, _, _v = _part.partition("=")\n'
          "        if not _k.strip():\n"
          "            continue\n"
          '        _out[_k.strip()] = [_x.strip() for _x in _v.split("/")]\n'
          "    return _out or None\n"
          "\n"
          '_W_ENV = _f5_env("F5_REC_FWD")\n'
          "REC_FWD = ({} if _W_ENV else {})\n"
          'if ("_rl" in SCAN or "_rl" in FWD):\n'
          + A_GUARD + "\n"
          "    if _miss:\n"
          "        sys.exit(\n"
          '            "REFUSING: F5_SCAN/F5_FWD name a suffixed deposit but '
          '%s unset.\\n"\n'
          + A_MSG + "\n"
          '            % ", ".join(_miss))\n'
          'DEFAULT_EXPR_MODEL = "composition"\n'
          + A_TECH +
          "NULL_EXPECT = 0.67\n"
          'print("T45", TECH45["SE"]["without composition adj."],\n'
          '      TECH45["MXE"]["with composition adj."])\n'
          'print("T9", TECH9["SE"]["without"], TECH9["MXE"]["with"])\n'
          'print("PROV", f"SE (663 events), 39 for MXE (1,393)")\n'
          'print("LOG", f"(rbm20_g45_406229.out step 5)")\n')

    new, probs = apply_edits(fx)
    assert not probs, probs
    assert not verify(new), verify(new)
    compile(new, "fx", "exec")

    d = tempfile.mkdtemp()
    KEYS = ("F5_SCAN", "F5_FWD", "F5_REC_SCAN", "F5_FLOOR_RANK", "F5_REC_FWD",
            "F5_TECH45", "F5_TECH9", "F5_NTARGET", "F5_D_LOG")

    def run(src, name, env=None):
        p = os.path.join(d, name)
        open(p, "w").write(src)
        e = dict(os.environ)
        for k in KEYS:
            e.pop(k, None)
        e.update(env or {})
        return subprocess.run([sys.executable, "-I", p], capture_output=True,
                              text=True, cwd=d, env=e)

    # --- unsuffixed: byte-identical to pre-patch
    r0 = run(fx, "pre.py")
    r1 = run(new, "post.py")
    assert r0.returncode == 0 and r1.returncode == 0, (r0.stderr, r1.stderr)
    assert r0.stdout == r1.stdout, (r0.stdout, r1.stdout)
    assert "T45 (0.95, 0.86) (2.71, 1.63)" in r1.stdout, r1.stdout
    assert "T9 (1.01, 0.93) (2.65, 1.65)" in r1.stdout, r1.stdout
    assert "SE (663 events), 39 for MXE (1,393)" in r1.stdout, r1.stdout
    assert "(rbm20_g45_406229.out step 5)" in r1.stdout, r1.stdout

    FULL = {"F5_SCAN": "signature_rbp_margins_rl125.tsv",
            "F5_FWD": "rbm20_margin_resolve_rl125.tsv",
            "F5_REC_SCAN": "x", "F5_FLOOR_RANK": "SE=6,MXE=3",
            "F5_REC_FWD": "SE=-0.051/-0.159/0.047",
            "F5_TECH45": "SE=0.96/0.87/1.02/0.92,MXE=1.98/1.24/2.65/1.63",
            "F5_TECH9": "SE=1.03/0.95/1.08/1.00,MXE=1.95/1.29/2.52/1.67",
            "F5_NTARGET": "SE=641,MXE=1399",
            "F5_D_LOG": "rbm20_g45_407951.out"}

    # --- all seven supplied: every plotted and printed value moves
    r2 = run(new, "post.py", FULL)
    assert r2.returncode == 0, r2.stderr
    assert "T45 (0.96, 0.87) (2.65, 1.63)" in r2.stdout, r2.stdout
    assert "T9 (1.03, 0.95) (2.52, 1.67)" in r2.stdout, r2.stdout
    assert "SE (641 events), 39 for MXE (1,399)" in r2.stdout, r2.stdout
    assert "(rbm20_g45_407951.out step 5)" in r2.stdout, r2.stdout

    # --- THE POINT OF 156: suffixed, gates supplied, panel D NOT -> refuse.
    #     This is the state 155 left, in which the figure built and every
    #     gate passed while panel D was 50 bp.
    part = {k: v for k, v in FULL.items()
            if k not in ("F5_TECH45", "F5_TECH9", "F5_NTARGET", "F5_D_LOG")}
    r3 = run(new, "post.py", part)
    assert r3.returncode != 0, r3.stdout
    out3 = r3.stdout + r3.stderr
    assert "F5_TECH45" in out3 and "F5_TECH9" in out3, out3
    assert "mix read lengths while" in out3, out3

    # --- missing only ONE of the four is still refused, and named
    for k in ("F5_TECH45", "F5_TECH9", "F5_NTARGET", "F5_D_LOG"):
        one = {kk: vv for kk, vv in FULL.items() if kk != k}
        rr = run(new, "post.py", one)
        assert rr.returncode != 0, (k, rr.stdout)
        assert k in rr.stdout + rr.stderr, (k, rr.stderr)

    # --- THE COMMA TRAP: F5_NTARGET is comma-separated, so a value
    #     containing a comma ("1,399") parses as two entries. Integers plus
    #     a format-time separator is why this now works; assert the output
    #     carries the comma even though the input does not.
    assert "MXE=1399" in FULL["F5_NTARGET"], FULL["F5_NTARGET"]
    assert "1,399" in r2.stdout, r2.stdout

    # --- an UNSUFFIXED run must not be forced to supply them
    r4 = run(new, "post.py", {"F5_TECH45": "SE=1/2/3/4,MXE=5/6/7/8"})
    assert r4.returncode == 0, r4.stderr
    assert "T45 (1.0, 2.0) (7.0, 8.0)" in r4.stdout, r4.stdout

    # --- refuses without 155
    _n, p5 = apply_edits("print('bare')\n")
    assert p5 and "155 has not been applied" in p5[0], p5

    # --- double application refused without altering the text
    n6, p6 = apply_edits(new)
    assert p6 and "already carries 156" in p6[0], p6
    assert n6 == new, "a refused re-application still altered the text"

    print("[ok] self-test passed: unsuffixed output is byte-identical to "
          "pre-patch (panel D still 0.95/0.86, provenance still 663/1,393 and "
          "the 5 Oct log); supplying all four moves every plotted AND printed "
          "value; a suffixed run with the GATES supplied but panel D NOT is "
          "REFUSED — the exact state 155 left, where the figure built and "
          "every gate passed while panel D was 50 bp; omitting any ONE of the "
          "four is refused and names it; an unsuffixed run is never forced to "
          "supply them; NULL_EXPECT untouched; refuses without 155; double "
          "application refused")
    return 0


if __name__ == "__main__":
    sys.exit(main())
