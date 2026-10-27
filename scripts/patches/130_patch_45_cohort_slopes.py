#!/usr/bin/env python3
"""
130_patch_45_cohort_slopes.py -- 45 is the gate for 82, 85 and 98.

WHY THIS SCRIPT MATTERS
-----------------------
`45_cohort_slopes.py` writes `cohort_slopes_dedup.tsv` and
`cohort_type_concordance_dedup.tsv`, which are read by

    82_build_s4.py            Supplementary Table 4
    85_recount_combinations.py
    98_build_replication_table.py   Table 4

so **98 is not blocked on 21 — it is blocked on 45**, and so are the other two.
45 appeared on no task list; it is the fourth producer found by grepping the
scripts rather than reading the docs, after 11, 33 and 03/03b.

(61 and 62 also match a grep for these filenames, but only in a comment —
"conventions copied from 45_cohort_slopes.py so the two cannot drift" — and are
not consumers.)

WHAT MOVES
----------
  line 109  load(f"{et}_v2filtered", et)                 the discovery table
  line  55  KNOWN_PCT   79.6 / 81.2 / 76.1 / 69.2 / ...  per-cohort concordance
  line  57  KNOWN_POOLED = (1636, 2163)                  the pooled figure
  line 102  --out-slopes  cohort_slopes_dedup.tsv        output
  line 103  --out-conc    cohort_type_concordance_dedup.tsv

`REPL` (line 43) stays unsuffixed: Colorado, Brigham, Stanford and Fuwai were
never re-quantified. Only the MAGNet discovery table carries the suffix.

Measured at 125 bp, from 33's section D and confirmed independently by 61:

    Colorado 311/385 = 80.8%     Brigham  298/364 = 81.9%
    Stanford 345/455 = 75.8%     Fuwai    422/596 = 70.8%
    LBNL     122/156 = 78.2%  (not independent; excluded from the pool)
    POOLED without LBNL          1,376/1,800 = 76.4%   (was 1,636/2,163)

THE KNOWN_PCT OVERRIDE NAMES ITS OWN KEYS
-----------------------------------------
`KNOWN_PCT` spans two source lines and this patch does not read the second, so
the override asserts **set equality** against whatever keys are actually there.
A missing or misspelt cohort fails loudly and the error names the expected list
— the convention 96's `EXPECT` and 61's `KNOWN_PCT` already use. If LBNL is not
a key, the assert will say so.

THE CROSS-CHECK THIS ENABLES
----------------------------
45 has its own `deming()` (line 81), independent of 33's. Its slopes are the
same quantity 33's section E computes:

    Colorado 0.521 -> 0.600   Brigham 0.476 -> 0.510
    Stanford 0.386 -> 0.401   Fuwai   0.243 -> 0.295

45's own header says it: *"If this script and Figure 4 ever disagree…"*. Two
independent implementations of a published slope agreeing is worth more than
either value; disagreeing matters more still.

USAGE
-----
  python 130_patch_45_cohort_slopes.py
  python 130_patch_45_cohort_slopes.py --apply
  python 130_patch_45_cohort_slopes.py --self-test
"""
import argparse
import difflib
import os
import shutil
import sys
import time

TARGET = "45_cohort_slopes.py"

A_NKEEP = 'NKEEP = {"SE": 3, "MXE": 5, "A3SS": 5, "A5SS": 5, "RI": 3}\n'
B_NKEEP = A_NKEEP + '''
# --------------------------------------------------------------------------
# Read length. Only the MAGNet DISCOVERY table carries this suffix; REPL below
# stays unsuffixed because no replication cohort was re-quantified.
SUF = os.environ.get("TBL_SUFFIX", "")
# --------------------------------------------------------------------------
'''

A_POOLED = "KNOWN_POOLED = (1636, 2163)\n"
B_POOLED = A_POOLED + '''
# --------------------------------------------------------------------------
# 125 bp overrides, placed AFTER the dicts they rewrite.
#
# Measured (33 section D, confirmed by 61):
#   KNOWN_PCT=Colorado=80.8,Brigham=81.9,Stanford=75.8,Fuwai=70.8,LBNL=78.2
#   KNOWN_POOLED=1376,1800
#
# The override must name EXACTLY the keys already present, so a typo or a
# missing cohort fails loudly and the error names the expected list.
_kp = os.environ.get("KNOWN_PCT", "").strip()
if _kp:
    _d = dict(_p.split("=", 1) for _p in _kp.split(","))
    assert set(_d) == set(KNOWN_PCT), (
        f"KNOWN_PCT override must name exactly {sorted(KNOWN_PCT)}, "
        f"got {sorted(_d)}")
    KNOWN_PCT = {_k: float(_v) for _k, _v in _d.items()}

_kpo = os.environ.get("KNOWN_POOLED", "").strip()
if _kpo:
    _k, _n = (int(_x) for _x in _kpo.split(","))
    KNOWN_POOLED = (_k, _n)
# --------------------------------------------------------------------------
'''

EDITS = [
    ("SUF", A_NKEEP, B_NKEEP, 1),
    ("the KNOWN_PCT and KNOWN_POOLED overrides", A_POOLED, B_POOLED, 1),
    ("the discovery table, suffixed",
     'load(f"{et}_v2filtered", et)',
     'load(f"{et}_v2filtered" + SUF, et)', 1),
    # The two output anchors are the DEFAULT STRINGS ONLY, not the whole
    # add_argument call: indentation and the space after the comma are things
    # I would be guessing at, and an anchor that guesses is an anchor that
    # refuses. The default literal is unique in the file either way.
    ("the slopes output, suffixed",
     'default="cohort_slopes_dedup.tsv"',
     'default="cohort_slopes_dedup" + SUF + ".tsv"', 1),
    ("the concordance output, suffixed",
     'default="cohort_type_concordance_dedup.tsv"',
     'default="cohort_type_concordance_dedup" + SUF + ".tsv"', 1),
]


def apply_edits(text):
    # Checked FIRST, before any replacement: A_NKEEP and A_POOLED are both
    # PREFIXES of their replacements, so on an already-patched file their
    # counts are still 1 and the override block would be inserted twice.
    if 'SUF = os.environ.get("TBL_SUFFIX", "")' in text:
        return text, ["SUF is already defined — 45 is already patched; "
                      "nothing to do (restore the .bak to re-apply)"]
    problems, out = [], text
    for name, anchor, repl, want in EDITS:
        got = out.count(anchor)
        if got != want:
            problems.append(f"{name}: anchor occurs {got}x, expected {want}x\n"
                            f"      anchor: "
                            f"{anchor.strip().splitlines()[0][:76]}")
            continue
        out = out.replace(anchor, repl, want)
    return out, problems


def verify(text):
    bad = []
    for w in ('SUF = os.environ.get("TBL_SUFFIX", "")',
              'load(f"{et}_v2filtered" + SUF, et)',
              'os.environ.get("KNOWN_PCT"',
              'os.environ.get("KNOWN_POOLED"',
              'default="cohort_slopes_dedup" + SUF + ".tsv"',
              'default="cohort_type_concordance_dedup" + SUF + ".tsv"'):
        if w not in text:
            bad.append(f"missing: {w[:58]}")
    for gone in ('load(f"{et}_v2filtered", et)',
                 'default="cohort_slopes_dedup.tsv"',
                 'default="cohort_type_concordance_dedup.tsv"'):
        if gone in text:
            bad.append(f"unsuffixed site survives: {gone[:50]}")
    # REPL must stay unsuffixed — no cohort was re-quantified
    if "_Colorado_DCMvsNF\" + SUF" in text or "Brigham_{ET}_HCM\" + SUF" in text:
        bad.append("a replication cohort was wrongly suffixed")
    # the overrides must come AFTER the literals they rewrite
    i_pct = text.find("KNOWN_PCT = {")
    i_ov = text.find('os.environ.get("KNOWN_PCT"')
    if i_pct < 0 or i_ov < 0 or i_pct > i_ov:
        bad.append("the KNOWN_PCT override is not after the literal dict")
    i_pool = text.find("KNOWN_POOLED = (1636, 2163)")
    i_po = text.find('os.environ.get("KNOWN_POOLED"')
    if i_pool < 0 or i_po < 0 or i_pool > i_po:
        bad.append("the KNOWN_POOLED override is not after the literal")
    # SUF above every consumer
    i_s = text.find('SUF = os.environ.get("TBL_SUFFIX"')
    for m in ('load(f"{et}_v2filtered" + SUF, et)',
              'default="cohort_slopes_dedup" + SUF'):
        j = text.find(m)
        if j < 0 or i_s > j:
            bad.append(f"SUF not defined above: {m[:44]}")
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
    print("130 — 45_cohort_slopes.py, the gate for 82, 85 and 98")
    print("=" * 78)
    if not os.path.exists(a.target):
        print(f"[fatal] {a.target} not found.", file=sys.stderr)
        return 2
    src = open(a.target, errors="replace").read()
    if "import os" not in src:
        print("[fatal] no `import os`", file=sys.stderr)
        return 2
    new, problems = apply_edits(src)
    if problems:
        print("REFUSING — anchors do not match:")
        for p in problems:
            print(f"  ✗ {p}")
        print("\n  Nothing written. Send me lines 40-60 and 100-112.")
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
        print("  # 1. CALIBRATION — unsuffixed must still reproduce")
        print("  python3 45_cohort_slopes.py 2>&1 | tail -25")
        print()
        print("  # 2. 125 bp. If the KNOWN_PCT keys are wrong the assert names")
        print("  #    the expected list — paste it and I will correct the line.")
        print("  TBL_SUFFIX=_rl125 \\")
        print("  KNOWN_PCT=Colorado=80.8,Brigham=81.9,Stanford=75.8,"
              "Fuwai=70.8,LBNL=78.2 \\")
        print("  KNOWN_POOLED=1376,1800 \\")
        print("  python3 45_cohort_slopes.py")
        print()
        print("  # 3. the convergence test against 33's section E")
        print("  cat cohort_slopes_dedup_rl125.tsv")
        print()
        print("⚠ EXPECTED SLOPES, from 33: Colorado 0.600, Brigham 0.510,")
        print("  Stanford 0.401, Fuwai 0.295 (LBNL 0.573, not independent).")
        print("  45 has its own deming() — two implementations agreeing is")
        print("  worth more than either number. Disagreeing matters more.")
        print()
        print("⬜ THEN 82, 85 and 98 can be re-based. 98 carries SEVEN stale")
        print("  constants (RECORD × 4 cohorts, POOLED_K/N, POOLED_P_1SF,")
        print("  MEAN_SLOPE) and all of them reach Table 4.")
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
    print("✓ verified on disk: the discovery table and both outputs follow")
    print("  SUF, both recorded sets are overridable, REPL is untouched.")
    return 0


def self_test():
    import subprocess
    import tempfile

    fx = (
        "import os\n"
        + A_NKEEP
        + 'REPL = [("Colorado", "{ET}_Colorado_DCMvsNF"), '
          '("Brigham", "Brigham_{ET}_HCM")]\n'
        + 'KNOWN_PCT = {"Colorado": 79.6, "Brigham": 81.2, "Stanford": 76.1,\n'
        + '             "Fuwai": 69.2, "LBNL": 75.9}\n'
        + A_POOLED
        + "def load(base, et):\n"
        + '    print("LOAD", base)\n'
        + "    return None, None\n"
        + "class _AP:\n"
        + "    def add_argument(self, *a, **k):\n"
        + "        print('ARG', a[0], k.get('default'))\n"
        + "ap = _AP()\n"
        + "def main():\n"
        + '    for et in ["SE"]:\n'
        + '        disc, dp = load(f"{et}_v2filtered", et)\n'
        + '    ap.add_argument("--out-slopes", default="cohort_slopes_dedup.tsv")\n'
        + '    ap.add_argument("--out-conc", '
          'default="cohort_type_concordance_dedup.tsv")\n'
        + "if __name__ == '__main__':\n"
        + '    print("PCT", KNOWN_PCT["Colorado"], "POOLED", KNOWN_POOLED)\n'
        + "    main()\n")

    new, probs = apply_edits(fx)
    assert not probs, probs
    assert not verify(new), verify(new)

    d = tempfile.mkdtemp()
    p = os.path.join(d, "fx.py")
    open(p, "w").write(new)

    def run(env=None):
        return subprocess.run([sys.executable, "-I", p], capture_output=True,
                              text=True, cwd=d,
                              env=dict(os.environ, **(env or {})))

    # 1. unsuffixed -> original everything
    r = run()
    assert r.returncode == 0, r.stdout + r.stderr
    assert "LOAD SE_v2filtered\n" in r.stdout, r.stdout
    assert "PCT 79.6 POOLED (1636, 2163)" in r.stdout, r.stdout
    assert "cohort_slopes_dedup.tsv" in r.stdout, r.stdout
    assert "cohort_type_concordance_dedup.tsv" in r.stdout, r.stdout

    # 2. suffixed + both overrides
    env = {"TBL_SUFFIX": "_rl125",
           "KNOWN_PCT": "Colorado=80.8,Brigham=81.9,Stanford=75.8,"
                        "Fuwai=70.8,LBNL=78.2",
           "KNOWN_POOLED": "1376,1800"}
    r2 = run(env)
    assert r2.returncode == 0, r2.stdout + r2.stderr
    assert "LOAD SE_v2filtered_rl125" in r2.stdout, r2.stdout
    assert "PCT 80.8 POOLED (1376, 1800)" in r2.stdout, r2.stdout
    assert "cohort_slopes_dedup_rl125.tsv" in r2.stdout, r2.stdout
    assert "cohort_type_concordance_dedup_rl125.tsv" in r2.stdout, r2.stdout
    # the replication cohort stems must be untouched
    assert "_Colorado_DCMvsNF_rl125" not in r2.stdout, r2.stdout

    # 3. a wrong KNOWN_PCT key set fails loudly and names the expected keys
    rb = run({"KNOWN_PCT": "Colorado=80.8,Brigham=81.9"})
    assert rb.returncode != 0, rb.stdout + rb.stderr
    both = rb.stdout + rb.stderr
    assert "must name exactly" in both, both
    assert "LBNL" in both and "Fuwai" in both, both

    # 4. double application refused, and refused WITHOUT having inserted the
    #    override block a second time
    _n, p4 = apply_edits(new)
    assert p4, "second application was not refused"
    assert _n == new, "a refused re-application still altered the text"
    assert _n.count('os.environ.get("KNOWN_PCT"') == 1, _n.count("KNOWN_PCT")

    # 5. the two output anchors must survive a DIFFERENT add_argument layout —
    #    the real file's indentation and comma spacing are not mine to assume
    fx2 = fx.replace(
        '    ap.add_argument("--out-slopes", default="cohort_slopes_dedup.tsv")',
        '    ap.add_argument(\n'
        '        "--out-slopes",\n'
        '        default="cohort_slopes_dedup.tsv",\n'
        '        help="x")')
    assert fx2 != fx, "fixture-2 rewrite did not apply"
    n2, pr2 = apply_edits(fx2)
    assert not pr2, pr2
    assert not verify(n2), verify(n2)
    p2 = os.path.join(d, "fx2.py")
    open(p2, "w").write(n2)
    r5 = subprocess.run([sys.executable, "-I", p2], capture_output=True,
                        text=True, cwd=d,
                        env=dict(os.environ, TBL_SUFFIX="_rl125"))
    assert r5.returncode == 0, r5.stdout + r5.stderr
    assert "cohort_slopes_dedup_rl125.tsv" in r5.stdout, r5.stdout

    print("[ok] self-test passed: the discovery table and both output names "
          "follow SUF while the replication cohort stems do not; KNOWN_PCT and "
          "KNOWN_POOLED override after the literals they rewrite; a wrong key "
          "set fails and the error names the expected cohorts; unsuffixed "
          "behaviour is unchanged; double application refused")
    return 0


if __name__ == "__main__":
    sys.exit(main())
