#!/usr/bin/env python3
"""
122_patch_93_parameterise.py -- give 93 a read-length parameter. It has none.

WHAT IS WRONG
-------------
`93_source_f6_foreground.py` exists to establish where Figure 6's **927** came
from — a number with no file behind it, printed twice on a published figure, one
digit from the Results headline 926. It is an audit script, and a good one.

But it has no suffix parameter at all. Its arguments are `--dir`, `--model`,
`--list`, and every table name is a literal:

    MODELS = {
        "compadj": "{cls}_v2compadj.all_events.tsv",
        "base":    "{cls}_v2filtered.all_events.tsv",
        "cmfrac":  "{cls}_v2cmfrac.all_events.tsv",
    }

So at 125 bp it silently audits the 50 bp tables. Unlike a renamed read that
fails loudly, this one succeeds and answers the wrong question.

It also carries two blocks of hard-coded counts:

    HARD = dict(ALL=927, UP=636, DOWN=291, MIXED=75, MULTI=134)
    KNOWN_SE = {"base": (926, ...), "cmfrac": (1054, ...),
                "compadj": (899, ...)}

and prints the KNOWN_SE values as a literal string at line 183:

    print(f"    ⚠ NOT CALIBRATED — reproduces no recorded count
            (926 / 1,054 / 899).")

`KNOWN_SE` is the calibration gate — the script's own rule is that its verdict
on 927 means nothing until it reproduces a known count. At 125 bp all three
knowns change, so the gate must be re-pointed or it declares every correct run
"NOT CALIBRATED".

`HARD` is different and must stay a supplied value: those are the numbers on the
figure. They have no file behind them, so they cannot be derived — only compared.

WHAT THIS PATCH DOES
--------------------
  * `SUF = os.environ.get("TBL_SUFFIX", "")`, appended to all three MODELS
    patterns
  * `HARD` overridable:     HARD=ALL=...,UP=...,DOWN=...,MIXED=...,MULTI=...
  * `KNOWN_SE` overridable: KNOWN_SE=base=813,cmfrac=939,compadj=777
  * the printed "(926 / 1,054 / 899)" derived from KNOWN_SE, so the message
    cannot disagree with the gate it describes

Both overrides assert the key set matches exactly, so a typo is a failure rather
than a silently ignored variable — the convention `96_lvef_column_dictionary.py`
already uses for its `EXPECT` override.

Defaults are the 50 bp values throughout. An unparameterised run reproduces the
original output exactly.

THE 125 bp KNOWN COUNTS
-----------------------
From `RESULT_s1_composition_rl125_complete.md` (76_build_s1_v2.py, confirmed
independently by 71 and 113):

    model                              50 bp   125 bp
    base         (v2filtered)            926      813
    cmfrac       (fraction_only)        1,054      939
    compadj      (nine_marker)            899      777

USAGE
-----
  python 122_patch_93_parameterise.py
  python 122_patch_93_parameterise.py --apply
  python 122_patch_93_parameterise.py --self-test
"""
import argparse
import difflib
import os
import shutil
import sys
import time

TARGET = "93_source_f6_foreground.py"

A_HARD = "HARD = dict(ALL=927, UP=636, DOWN=291, MIXED=75, MULTI=134)\n"

B_HARD = A_HARD + '''
# ---------------------------------------------------------------------------
# 125 bp parameterisation. Every default below is the 50 bp value, so an
# unparameterised run reproduces the original output exactly.
SUF = os.environ.get("TBL_SUFFIX", "")   # "" = 50 bp, "_rl125" = 125 bp

# HARD are the numbers PRINTED ON FIGURE 6 — the ones this script exists to
# source. They have no file behind them, so they cannot be derived, only
# supplied and compared. At a new read length there is no recorded equivalent;
# leave this unset and read the comparison column as "no record to check".
#     HARD=ALL=...,UP=...,DOWN=...,MIXED=...,MULTI=...
_hv = os.environ.get("HARD", "").strip()
if _hv:
    _hd = dict(_p.split("=", 1) for _p in _hv.split(","))
    assert set(_hd) == set(HARD), (
        f"HARD override must name exactly {sorted(HARD)}, got {sorted(_hd)}")
    HARD = {_k: int(_v) for _k, _v in _hd.items()}
# ---------------------------------------------------------------------------
'''

A_MODELS = [
    ('    "compadj": "{cls}_v2compadj.all_events.tsv",\n',
     '    "compadj": "{cls}_v2compadj" + SUF + ".all_events.tsv",\n'),
    ('    "base": "{cls}_v2filtered.all_events.tsv",\n',
     '    "base": "{cls}_v2filtered" + SUF + ".all_events.tsv",\n'),
    ('    "cmfrac": "{cls}_v2cmfrac.all_events.tsv",\n',
     '    "cmfrac": "{cls}_v2cmfrac" + SUF + ".all_events.tsv",\n'),
]

A_KNOWN = '''KNOWN_SE = {
    "base": (926, "base model, no composition term — the Results headline"),
    "cmfrac": (1054, "restricted cardiomyocyte-fraction score — Results primary"),
    "compadj": (899, "nine-marker score — Results conservative bound"),
}
'''

B_KNOWN = A_KNOWN + '''
# KNOWN_SE is the CALIBRATION GATE. This script's own rule is that its verdict
# on an unknown number (927) means nothing until it reproduces a known one, so
# at a new read length the gate must be re-pointed or every correct run is
# declared "NOT CALIBRATED".
#
# 125 bp (RESULT_s1_composition_rl125_complete.md; 76, confirmed by 71 and 113):
#     KNOWN_SE=base=813,cmfrac=939,compadj=777
#
#     model                      50 bp   125 bp
#     base     (v2filtered)        926      813
#     cmfrac   (fraction_only)   1,054      939
#     compadj  (nine_marker)       899      777
_kv = os.environ.get("KNOWN_SE", "").strip()
if _kv:
    _kd = dict(_p.split("=", 1) for _p in _kv.split(","))
    assert set(_kd) == set(KNOWN_SE), (
        f"KNOWN_SE override must name exactly {sorted(KNOWN_SE)}, "
        f"got {sorted(_kd)}")
    KNOWN_SE = {_k: (int(_v), KNOWN_SE[_k][1]) for _k, _v in _kd.items()}

# The printed list derives from the gate, so the message cannot quote a read
# length the gate is not checking. It used to be the literal "926 / 1,054 / 899".
KNOWN_LIST = " / ".join(f"{_v[0]:,}" for _v in KNOWN_SE.values())
'''

A_MSG = ('        print(f"    ⚠ NOT CALIBRATED — reproduces no recorded count '
         '(926 / 1,054 / 899).")\n')
B_MSG = ('        print(f"    ⚠ NOT CALIBRATED — reproduces no recorded count '
         '({KNOWN_LIST}).")\n')

EDITS = ([("SUF and the HARD override", A_HARD, B_HARD, 1)]
         + [(f"MODELS: {o.strip()[:28]}", o, n, 1) for o, n in A_MODELS]
         + [("the KNOWN_SE override and KNOWN_LIST", A_KNOWN, B_KNOWN, 1),
            ("the printed known-count list", A_MSG, B_MSG, 1)])


def apply_edits(text):
    problems, out = [], text
    for name, anchor, repl, want in EDITS:
        got = out.count(anchor)
        if got != want:
            problems.append(
                f"{name}: anchor occurs {got}x, expected {want}x\n"
                f"      anchor: {anchor.strip().splitlines()[0][:84]}")
            continue
        if "SUF = os.environ.get" in out and name.startswith("SUF"):
            problems.append("SUF is already defined — already patched")
            continue
        out = out.replace(anchor, repl, want)
    return out, problems


def verify(text):
    bad = []
    code = "\n".join(l for l in text.splitlines()
                     if not l.lstrip().startswith("#"))
    if 'SUF = os.environ.get("TBL_SUFFIX", "")' not in code:
        bad.append("SUF is not defined")
    for stem in ("v2compadj", "v2filtered", "v2cmfrac"):
        if f'"{{cls}}_{stem}" + SUF + ".all_events.tsv"' not in code:
            bad.append(f"MODELS entry for {stem} was not parameterised")
        if f'"{{cls}}_{stem}.all_events.tsv"' in code:
            bad.append(f"the unparameterised {stem} pattern survives")
    if 'os.environ.get("HARD"' not in code:
        bad.append("the HARD override is missing")
    if 'os.environ.get("KNOWN_SE"' not in code:
        bad.append("the KNOWN_SE override is missing")
    if "KNOWN_LIST" not in code:
        bad.append("KNOWN_LIST is not derived")
    if "926 / 1,054 / 899" in code:
        bad.append("the printed 50 bp known-count list survives on a code line")
    # SUF must be defined above MODELS, KNOWN_LIST above its consumer
    i_suf, i_mod = text.find("SUF = os.environ"), text.find("MODELS = {")
    if i_suf < 0 or i_mod < 0 or i_suf > i_mod:
        bad.append("SUF is not defined above MODELS")
    i_kl, i_use = text.find("KNOWN_LIST = "), text.find("({KNOWN_LIST})")
    if i_kl < 0 or i_use < 0 or i_kl > i_use:
        bad.append("KNOWN_LIST is not defined above its consumer")
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
    print("122 — give 93 a read-length parameter")
    print("=" * 78)
    if not os.path.exists(a.target):
        print(f"[fatal] {a.target} not found.", file=sys.stderr)
        return 2
    src = open(a.target, errors="replace").read()
    if "import os" not in src:
        print("[fatal] 93 does not import os; the overrides need it.",
              file=sys.stderr)
        return 2

    new, problems = apply_edits(src)
    if problems:
        print("REFUSING — anchors do not match:")
        for p in problems:
            print(f"  ✗ {p}")
        print("\n  Nothing written. Send me lines 55-85 and line 183 and I")
        print("  will re-anchor.")
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
        print("DRY RUN — nothing written. Re-run with --apply.")
        print()
        print("Then, in this order:")
        print()
        print("  # 1. CALIBRATION — unparameterised must still give the 50 bp")
        print("  #    answer. If this changes, the patch broke something.")
        print("  python3 93_source_f6_foreground.py --model all")
        print()
        print("  # 2. 125 bp. HARD left UNSET on purpose: Figure 6's numbers")
        print("  #    have no 125 bp equivalent yet, so the comparison column")
        print("  #    is there to be read, not matched.")
        print("  TBL_SUFFIX=_rl125 \\")
        print("  KNOWN_SE=base=813,cmfrac=939,compadj=777 \\")
        print("  python3 93_source_f6_foreground.py --model all")
        print()
        print("⚠ ALSO CHECK BY HAND — line ~236 carries the same")
        print("  '926 / 1,054 / 899' literal in a second message. This patch")
        print("  changed only line 183, because the other occurrence may sit in")
        print("  a plain (non-f) string where '{KNOWN_LIST}' would print braces:")
        print()
        print("  grep -n '926 / 1,054 / 899' 93_source_f6_foreground.py")
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
    print("✓ verified on disk: SUF above MODELS, all three patterns")
    print("  parameterised, both overrides present, KNOWN_LIST derived.")
    return 0


def self_test():
    import subprocess
    import tempfile

    runnable = (
        "import os\n"
        + A_HARD
        + 'CLASSES = ["SE"]\n'
        + "MODELS = {\n"
        + A_MODELS[0][0] + A_MODELS[1][0] + A_MODELS[2][0]
        + "}\n"
        + A_KNOWN
        + "def report(se_got):\n"
        "    if True:\n"
        + A_MSG
        + "    print('MODELS', MODELS['base'].format(cls='SE'))\n"
        "    print('HARD_ALL', HARD['ALL'])\n"
        "    print('KNOWN_BASE', KNOWN_SE['base'][0])\n")

    new, probs = apply_edits(runnable)
    assert not probs, probs
    assert not verify(new), verify(new)

    d = tempfile.mkdtemp()
    p = os.path.join(d, "fx.py")
    open(p, "w").write(new + "\nif __name__ == '__main__':\n    report(0)\n")

    # 1. unparameterised -> the 50 bp behaviour, unchanged
    r = subprocess.run([sys.executable, "-I", p], capture_output=True,
                       text=True, cwd=d)
    assert r.returncode == 0, r.stdout + r.stderr
    assert "SE_v2filtered.all_events.tsv" in r.stdout, r.stdout
    assert "926 / 1,054 / 899" in r.stdout, r.stdout
    assert "HARD_ALL 927" in r.stdout, r.stdout

    # 2. suffixed + re-pointed gate
    env = dict(os.environ, TBL_SUFFIX="_rl125",
               KNOWN_SE="base=813,cmfrac=939,compadj=777")
    r2 = subprocess.run([sys.executable, "-I", p], capture_output=True,
                        text=True, cwd=d, env=env)
    assert r2.returncode == 0, r2.stdout + r2.stderr
    assert "SE_v2filtered_rl125.all_events.tsv" in r2.stdout, r2.stdout
    assert "813 / 939 / 777" in r2.stdout, r2.stdout
    assert "926" not in r2.stdout, r2.stdout
    assert "HARD_ALL 927" in r2.stdout, "HARD should keep its default"

    # 3. HARD override works
    env3 = dict(os.environ, HARD="ALL=800,UP=500,DOWN=300,MIXED=60,MULTI=100")
    r3 = subprocess.run([sys.executable, "-I", p], capture_output=True,
                        text=True, cwd=d, env=env3)
    assert r3.returncode == 0, r3.stdout + r3.stderr
    assert "HARD_ALL 800" in r3.stdout, r3.stdout

    # 4. a typo in either override must FAIL, not be ignored
    for bad_env in ({"HARD": "ALL=1,UP=1,DOWN=1,MIXED=1"},
                    {"KNOWN_SE": "base=813,cmfrac=939"},
                    {"KNOWN_SE": "base=813,cmfrac=939,compad=777"}):
        rb = subprocess.run([sys.executable, "-I", p], capture_output=True,
                            text=True, cwd=d, env=dict(os.environ, **bad_env))
        assert rb.returncode != 0, (bad_env, rb.stdout + rb.stderr)
        assert "must name exactly" in (rb.stdout + rb.stderr), rb.stderr

    # 5. double application refused
    _n, p5 = apply_edits(new)
    assert p5, "second application was not refused"

    print("[ok] self-test passed: unparameterised keeps the 50 bp tables, the "
          "927 figure values and the '926 / 1,054 / 899' message; TBL_SUFFIX "
          "repoints all three MODELS patterns and KNOWN_SE repoints the gate "
          "and the printed list together; HARD overrides independently; a "
          "missing or misspelt key in either override fails loudly; double "
          "application refused")
    return 0


if __name__ == "__main__":
    sys.exit(main())
