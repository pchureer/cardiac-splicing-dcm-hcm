#!/usr/bin/env python3
"""
140_patch_82_drop_genome_wide.py -- remove GENOME_WIDE from S4.

CP's decision, 20 Oct: drop it entirely.

WHY
---
`GENOME_WIDE` holds script 21's **pre-exclusion, JCEC-based** genome-wide
slopes on the matched JC-intersect-JCEC subset, recorded in
`RESULT_magnitude_resolved.md` — a document that carries the banner

    ⚠⚠ THE DEMING SLOPES IN THIS DOCUMENT ARE SUPERSEDED …
    Quote nothing numeric from this document in the manuscript or the
    supplement. `45_cohort_slopes.py`'s output is the source for slopes.

They were being written into a supplementary table in defiance of that.
`62_build_s5.py` exists for exactly this reason — its header says putting the
recorded JC-versus-JCEC table in the supplement "would have published two
different slopes for the same cohort" — and 62 guards against it by aborting
unless its own JC column reproduces the manuscript's slopes. 82 had no
equivalent guard.

THE VALUES ARE NOT UNSOURCED — AN EARLIER CLAIM OF MINE WAS WRONG
-----------------------------------------------------------------
On 20 Oct I wrote that Brigham 0.757 and Fuwai 0.359 "match nothing in any
surviving output". They match `RESULT_magnitude_resolved.md` exactly. My
search covered only `symmetry_*.tsv`, which read `_dedup` tables — all
post-exclusion JC — so a pre-exclusion JCEC value could never have appeared
in it. `62_build_s5.py:16` named that document in plain text and I wrote a
defect note before reading it.

Three distinct Brigham quantities exist and must never be interchanged:

    0.476   manuscript — post-exclusion, JC, MAGNet-significant events (45)
    0.662   post-exclusion, JC, unrestricted forward slope (22)
    0.757   pre-exclusion, JCEC, matched JC-intersect-JCEC subset (21/20)

The same document explains the third: "restricting to the JC-intersect-JCEC
event set alone raises Brigham from 0.662 to 0.757 … That is an event-set
effect, not read length."

THE DEFECT THIS ALSO FIXES
--------------------------
`138_patch_82_s4.py` populated `S4_GENOME_WIDE` with 22's **125 bp FORWARD**
slopes (0.674/0.722/0.365/0.605) — post-exclusion, JC, unrestricted. The
50 bp column held pre-exclusion JCEC matched-subset values. So
`S4_rl125.tsv` as built on 20 Oct carried a column measuring a different
quantity from its 50 bp counterpart under the same label. That is the silent
redefinition this whole re-baselining exists to prevent, introduced by my own
patch, and removing the column removes it.

⚠ THIS CHANGES THE 50 bp TABLE TOO
----------------------------------
S4 goes from **11 rows to 7** at both read lengths (3 summary + 4 PER_COHORT).
The 50 bp calibration can therefore no longer be a plain `cmp`; it becomes
"the seven surviving rows are byte-identical and the four GENOME-WIDE rows
are gone". That is the correct consequence of the source document's
prohibition, but it means S4 is reissued at both read lengths, not only at
125 bp.

Nothing downstream reads S4 (`grep -ln 'S4\\.tsv|S4_rl125|S4\\.md'` matches
only patches 138 and 139), so no figure or tracker breaks.

WHAT MOVES
----------
  the dict                 -> retired to a comment, provenance written out
  _S4_REQ                  -> two groups, not three
  the refusal's command    -> S4_GENOME_WIDE segment removed
  the override block       -> removed
  the four table rows      -> removed
  the "two quantities" caveat and 139's derived gap -> removed

USAGE
-----
  python 140_patch_82_drop_genome_wide.py [--apply] [--self-test]
"""
import argparse
import difflib
import os
import shutil
import sys
import time

TARGET = "82_build_s4.py"

# ------------------------------------------------------------------- edit 1
A_DICT = ('# script 21\'s GENOME-WIDE slopes -- a DIFFERENT quantity. '
          'Never mix the two.\n'
          'GENOME_WIDE = {"Colorado": 0.571, "Brigham": 0.757, '
          '"Fuwai": 0.359, "Stanford": 0.578}\n')
B_DICT = '''# --------------------------------------------------------------- RETIRED
# GENOME_WIDE = {"Colorado": 0.571, "Brigham": 0.757,
#                "Fuwai": 0.359, "Stanford": 0.578}
#
# DROPPED FROM S4 on 20 Oct. Deliberate, and recorded here rather than lost.
#
# These are script 21's PRE-EXCLUSION, JCEC-based genome-wide slopes on the
# matched JC-intersect-JCEC subset, from RESULT_magnitude_resolved.md — a
# document whose own banner reads "Quote nothing numeric from this document
# in the manuscript or the supplement". They were being written into a
# supplementary table in defiance of that. 62_build_s5.py exists for the same
# reason and guards against it by aborting unless its JC column reproduces
# the manuscript's slopes; 82 had no equivalent guard.
#
# They are neither wrong nor unsourced — they are superseded. Three distinct
# Brigham quantities exist and must never be interchanged:
#     0.476  manuscript — post-exclusion, JC, MAGNet-significant (45)
#     0.662  post-exclusion, JC, unrestricted forward (22)
#     0.757  pre-exclusion, JCEC, matched subset (21/20)   <- this one
# and the source document explains the third: restricting to the
# JC-intersect-JCEC set raises Brigham from 0.662 to 0.757, an event-set
# effect rather than a read-length one.
#
# No 125 bp counterpart exists, because nothing at 125 bp is pre-exclusion
# JCEC. Patch 138 briefly filled this with 22's 125 bp FORWARD slopes — a
# different quantity under the same label — which is the defect removing the
# column fixes.
# --------------------------------------------------------------------------
'''

# ------------------------------------------------------------------- edit 2
A_REQ = '_S4_REQ = ("S4_RECORD", "S4_PER_COHORT", "S4_GENOME_WIDE")\n'
B_REQ = '_S4_REQ = ("S4_RECORD", "S4_PER_COHORT")\n'

# ------------------------------------------------------------------- edit 3
# 138 wrote this literal as repr(RUN125.replace("\n", "\n    ")), so the file
# holds THREE backslashes then 'n' per continuation. A first attempt here
# hand-escaped two and matched nothing. Rather than count backslashes off a
# terminal dump, rebuild both sides by the SAME computation 138 used: exact by
# construction, and immune to how any shell or pager renders them.
_RUN_OLD_SRC = (
    "TBL_SUFFIX=_rl125 \\\n"
    "S4_RECORD=common=161,A=0.690,B=0.968,A_CI=0.570:0.806,"
    "B_CI=0.875:1.069,dirA=130,dirB=140,product=0.668,C=0.600 \\\n"
    "S4_PER_COHORT=Colorado=0.5999,Brigham=0.5099,"
    "Stanford=0.4008,Fuwai=0.2948 \\\n"
    "S4_GENOME_WIDE=Colorado=0.674,Brigham=0.722,"
    "Fuwai=0.365,Stanford=0.605 \\\n"
    "python3 82_build_s4.py")
_RUN_NEW_SRC = (
    "TBL_SUFFIX=_rl125 \\\n"
    "S4_RECORD=common=161,A=0.690,B=0.968,A_CI=0.570:0.806,"
    "B_CI=0.875:1.069,dirA=130,dirB=140,product=0.668,C=0.600 \\\n"
    "S4_PER_COHORT=Colorado=0.5999,Brigham=0.5099,"
    "Stanford=0.4008,Fuwai=0.2948 \\\n"
    "python3 82_build_s4.py")
A_CMD = repr(_RUN_OLD_SRC.replace("\n", "\n    "))
B_CMD = repr(_RUN_NEW_SRC.replace("\n", "\n    "))
assert "S4_GENOME_WIDE" in A_CMD and "S4_GENOME_WIDE" not in B_CMD
assert A_CMD.count("S4_PER_COHORT") == B_CMD.count("S4_PER_COHORT") == 1

# ------------------------------------------------------------------- edit 4
A_OVR = '''# ⚠ Brigham and Fuwai are REPLACED, not updated. The recorded 0.757 and 0.359
# match nothing in any surviving symmetry output, at any read length, in any
# column, for any event type; 22's byte-identical 50 bp run gives 0.662 and
# 0.328. See DEFECT_genome_wide_brigham_fuwai_unsourced.md. And note that
# symmetry_SE_rl125.tsv's Brigham unsel_deming 0.757203 is NOT a confirmation:
# GENOME_WIDE holds FORWARD slopes, so that is the wrong column at the wrong
# read length.
if _s4env["S4_GENOME_WIDE"]:
    _d = dict(_p.split("=", 1) for _p in _s4env["S4_GENOME_WIDE"].split(","))
    assert set(_d) == set(GENOME_WIDE), (
        "S4_GENOME_WIDE must name exactly %s, got %s"
        % (sorted(GENOME_WIDE), sorted(_d)))
    GENOME_WIDE = {_k: float(_v) for _k, _v in _d.items()}
'''
B_OVR = '''# GENOME_WIDE is retired — see the comment above its former definition.
# There is deliberately no S4_GENOME_WIDE override: the quantity has no 125 bp
# counterpart, so there is nothing to override it with.
'''

# ------------------------------------------------------------------- edit 5
A_ROWS = '''    for coh, sl in GENOME_WIDE.items():
        rows.append(dict(part="GENOME-WIDE slope -- a DIFFERENT quantity",
                         arm=coh, what="do not mix with the row above",
                         comparison=f"MAGNet -> {coh}, all tested events",
                         n_events=np.nan, slope=sl, ci_lo=np.nan, ci_hi=np.nan,
                         same_direction="", lambda_convention="",
                         lambda_value=np.nan, recorded=sl,
                         source="script 21 / RESULT_magnitude_resolved.md "
                                "(pre-exclusion)"))
'''
B_ROWS = '''    # The four GENOME-WIDE rows were removed on 20 Oct: they carried
    # pre-exclusion JCEC values from a document that forbids quoting them.
    # S4 now presents ONE slope quantity, which is why the "two different
    # quantities" caveat below is also gone. 11 rows -> 7.
'''

# ------------------------------------------------------------------- edit 6
A_CAV = '''    print("  ⚠ THE TABLE HOLDS TWO DIFFERENT SLOPE QUANTITIES, LABELLED AS SUCH:")
    print("    the Figure 4 transfer slopes (on MAGNet-significant events) and")
    # DERIVED, not transcribed. The old literal 0.281 was Brigham
    # 0.757 - 0.476, and 0.757 matches nothing in any surviving
    # symmetry output — see
    # DEFECT_genome_wide_brigham_fuwai_unsourced.md. Computing it
    # also stops it going stale and stops it assuming Brigham.
    _gap = {_c: abs(GENOME_WIDE[_c] - PER_COHORT[_c])
            for _c in GENOME_WIDE if _c in PER_COHORT}
    _worst = max(_gap, key=_gap.get)
    print(f"    script 21's genome-wide slopes. They disagree by up to {_gap[_worst]:.3f} for")
    print(f"    {_worst}. Both are correct for what they are; a sentence must say which.")
'''
B_CAV = '''    print("  ⚠ THIS TABLE CARRIES ONE SLOPE QUANTITY ONLY: the Figure 4")
    print("    transfer slopes, on MAGNet-significant SE events, post-exclusion,")
    print("    JC, from 45_cohort_slopes.py. The pre-exclusion JCEC genome-wide")
    print("    slopes that once sat beside them were removed on 20 Oct — their")
    print("    source document forbids quoting them. If a reviewer asks why")
    print("    Figure 4 differs from any earlier circulated table, that is the")
    print("    answer, and it belongs in a reply rather than in the supplement.")
'''

EDITS = [
    ("the GENOME_WIDE dict, retired to a comment", A_DICT, B_DICT, 1),
    ("_S4_REQ, two groups", A_REQ, B_REQ, 1),
    ("the refusal's command line", A_CMD, B_CMD, 1),
    ("the S4_GENOME_WIDE override block", A_OVR, B_OVR, 1),
    ("the four GENOME-WIDE table rows", A_ROWS, B_ROWS, 1),
    ("the two-quantities caveat and 139's derived gap", A_CAV, B_CAV, 1),
]


def apply_edits(text):
    if "# --------------------------------------------------------------- RETIRED" in text:
        return text, ["GENOME_WIDE is already retired — 82 already carries "
                      "140 (restore the .bak to re-apply)"]
    for need, why in (
            ('SUF = os.environ.get("TBL_SUFFIX", "")',
             "138 has not been applied"),
            ("_gap = {_c: abs(GENOME_WIDE[_c]",
             "139 has not been applied")):
        if need not in text:
            return text, [f"{why} — apply it first; 140 removes what it added"]
    problems, out = [], text
    for name, anchor, repl, want in EDITS:
        got = out.count(anchor)
        if got != want:
            problems.append(
                f"{name}: anchor occurs {got}x, expected {want}x\n"
                f"      anchor: {anchor.strip().splitlines()[0][:72]}")
            continue
        out = out.replace(anchor, repl, want)
    return out, problems


def verify(text):
    bad = []
    # NOT ONE executable reference to GENOME_WIDE may survive. Comment lines
    # are allowed, and are the point.
    for i, line in enumerate(text.splitlines(), 1):
        if "GENOME_WIDE" not in line:
            continue
        if line.lstrip().startswith("#"):
            continue
        bad.append(f"line {i}: an executable GENOME_WIDE reference survives: "
                   f"{line.strip()[:56]}")
    for w in ("# --------------------------------------------------------------- RETIRED",
              "#                \"Fuwai\": 0.359, \"Stanford\": 0.578}",
              "0.757  pre-exclusion, JCEC, matched subset (21/20)",
              '_S4_REQ = ("S4_RECORD", "S4_PER_COHORT")',
              "There is deliberately no S4_GENOME_WIDE override",
              "11 rows -> 7",
              "ONE SLOPE QUANTITY ONLY"):
        if w not in text:
            bad.append(f"missing: {w[:58]}")
    for gone in ('_S4_REQ = ("S4_RECORD", "S4_PER_COHORT", "S4_GENOME_WIDE")',
                 "S4_GENOME_WIDE=Colorado=0.674",
                 "TWO DIFFERENT SLOPE QUANTITIES",
                 "part=\"GENOME-WIDE slope -- a DIFFERENT quantity\"",
                 "_worst = max(_gap"):
        if gone in text:
            bad.append(f"stale site survives: {gone[:50]}")
    # the PER_COHORT rows must be untouched
    if 'source="45_cohort_slopes.py (NOT recomputed here)"' not in text:
        bad.append("the PER_COHORT row block was altered — it must not be")
    if "for coh, sl in PER_COHORT.items():" not in text:
        bad.append("the PER_COHORT loop was removed")
    # the refusal must still name the two surviving groups
    if "S4_RECORD=common=161" not in text or "S4_PER_COHORT=Colorado=" not in text:
        bad.append("the refusal's command line lost a surviving group")
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
    print("140 — dropping GENOME_WIDE from S4 (CP's decision, 20 Oct)")
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
        print("  # 0. keep the 11-row 50 bp table as the comparator")
        print("  cp S4.tsv S4.ELEVENROW")
        print()
        print("  # 1. CALIBRATION, adjusted: 11 rows -> 7, and the seven")
        print("  #    survivors must be byte-identical")
        print("  python3 82_build_s4.py 2>&1 | tail -14")
        print("  wc -l S4.ELEVENROW S4.tsv       # expect 12 and 8")
        print("  grep -v 'GENOME-WIDE slope' S4.ELEVENROW > S4.expect")
        print("  cmp S4.tsv S4.expect \\")
        print("    && echo 'SURVIVING ROWS BYTE-IDENTICAL' \\")
        print("    || diff S4.expect S4.tsv")
        print()
        print("  # 2. the refusal must now name TWO groups, not three")
        print("  TBL_SUFFIX=_rl125 python3 82_build_s4.py 2>&1 | head -12")
        print()
        print("  # 3. rebuild 125 bp WITHOUT S4_GENOME_WIDE")
        print("  TBL_SUFFIX=_rl125 \\")
        print("  S4_RECORD=common=161,A=0.690,B=0.968,A_CI=0.570:0.806,"
              "B_CI=0.875:1.069,dirA=130,dirB=140,product=0.668,C=0.600 \\")
        print("  S4_PER_COHORT=Colorado=0.5999,Brigham=0.5099,"
              "Stanford=0.4008,Fuwai=0.2948 \\")
        print("  python3 82_build_s4.py 2>&1 | tail -14")
        print("  cut -f1 S4_rl125.tsv | sort -u")
        print()
        print("⚠ If passing S4_GENOME_WIDE now has no effect, that is correct —")
        print("  it is no longer read. The old S4_rl125.tsv was built WITH it")
        print("  and must be overwritten by step 3.")
        print()
        print("⬜ STILL TO DO on this thread: 62's line 16 keeps its statement")
        print("  (it was right all along), but 62's own MANUSCRIPT_SLOPE gate")
        print("  and its parts (b) 172/167 and (c) 79.6/81.2/69.2/76.1 are all")
        print("  50 bp and will refuse at 125 bp.")
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
    print("✓ verified on disk: no executable GENOME_WIDE reference survives,")
    print("  the provenance is written into the retiring comment, the")
    print("  PER_COHORT rows are untouched, and the refusal names two groups.")
    return 0


def self_test():
    import subprocess
    import tempfile

    fx = ("import os\n"
          "import numpy as np\n"
          "import pandas as pd\n"
          'SUF = os.environ.get("TBL_SUFFIX", "")\n'
          'PER_COHORT = {"Colorado": 0.521, "Brigham": 0.476, '
          '"Stanford": 0.386, "Fuwai": 0.243}\n'
          + A_DICT
          + A_REQ
          + '_s4env = {_k: os.environ.get(_k, "").strip() for _k in _S4_REQ}\n'
            "if SUF:\n"
            "    _miss = [_k for _k in _S4_REQ if not _s4env[_k]]\n"
            "    if _miss:\n"
            "        raise SystemExit(\n"
            '            "REFUSING — TBL_SUFFIX=%s but these are unset: '
            '%s\\n\\n    %s\\n"\n'
            # the embedded command is A_CMD itself, exactly as 138 wrote it,
            # so the fixture exercises the real byte sequence rather than a
            # hand-escaped lookalike
            '            % (SUF, ", ".join(_miss), ' + A_CMD + '))\n'
          + A_OVR
          + "def main():\n"
            "    rows = []\n"
            "    a = type('A', (), {'out': 'S4test.tsv'})()\n"
          + A_ROWS.replace("    for coh, sl in GENOME_WIDE",
                           "    for coh, sl in GENOME_WIDE")
          + "    for coh, sl in PER_COHORT.items():\n"
            "        rows.append(dict(part='transfer slope', arm=coh,\n"
            "                         slope=sl, recorded=sl,\n"
            '                         source="45_cohort_slopes.py '
            '(NOT recomputed here)"))\n'
            '    print("ROWS", len(rows))\n'
            '    print("PARTS", sorted({r["part"] for r in rows}))\n'
          + A_CAV
          + "if __name__ == '__main__':\n"
            "    main()\n")

    new, probs = apply_edits(fx)
    assert not probs, probs
    assert not verify(new), verify(new)
    compile(new, "fx", "exec")

    d = tempfile.mkdtemp()
    p = os.path.join(d, "fx.py")
    open(p, "w").write(new)

    def run(env=None):
        e = dict(os.environ)
        for k in ("TBL_SUFFIX", "S4_RECORD", "S4_PER_COHORT",
                  "S4_GENOME_WIDE"):
            e.pop(k, None)
        e.update(env or {})
        return subprocess.run([sys.executable, "-I", p], capture_output=True,
                              text=True, cwd=d, env=e)

    # 1. the four GENOME-WIDE rows are gone; the PER_COHORT rows remain
    r = run()
    assert r.returncode == 0, r.stdout + r.stderr
    assert "ROWS 4" in r.stdout, r.stdout          # fixture has only the 4 PC
    assert "GENOME-WIDE" not in r.stdout, r.stdout
    assert "transfer slope" in r.stdout, r.stdout
    assert "ONE SLOPE QUANTITY ONLY" in r.stdout, r.stdout
    assert "TWO DIFFERENT" not in r.stdout, r.stdout
    assert "disagree by up to" not in r.stdout, r.stdout

    # 2. the refusal names TWO groups and no longer mentions S4_GENOME_WIDE
    r2 = run({"TBL_SUFFIX": "_rl125"})
    assert r2.returncode != 0, r2.stdout + r2.stderr
    b2 = r2.stdout + r2.stderr
    assert "S4_RECORD" in b2 and "S4_PER_COHORT" in b2, b2
    assert "S4_GENOME_WIDE" not in b2, b2
    assert "S4_RECORD=common=161" in b2, b2
    assert "S4_PER_COHORT=Colorado=" in b2, b2

    # 3. passing S4_GENOME_WIDE now has NO effect and does not error
    r3 = run({"S4_GENOME_WIDE": "Colorado=9.9,Brigham=9.9,"
                                "Fuwai=9.9,Stanford=9.9"})
    assert r3.returncode == 0, r3.stdout + r3.stderr
    assert "9.9" not in r3.stdout, r3.stdout
    assert "ROWS 4" in r3.stdout, r3.stdout

    # 4. the provenance survives in the source, as comments only
    assert "0.757  pre-exclusion, JCEC, matched subset (21/20)" in new
    assert "0.662  post-exclusion, JC, unrestricted forward (22)" in new
    assert "Quote nothing numeric from this document" in new
    for line in new.splitlines():
        if "GENOME_WIDE" in line:
            assert line.lstrip().startswith("#"), line

    # 5. refuses if 138 or 139 is missing
    _n, p5 = apply_edits("print('bare')\n")
    assert p5 and "138 has not been applied" in p5[0], p5
    _n, p6 = apply_edits('SUF = os.environ.get("TBL_SUFFIX", "")\n')
    assert p6 and "139 has not been applied" in p6[0], p6

    # 6. double application refused, text unchanged
    n7, p7 = apply_edits(new)
    assert p7, "second application was not refused"
    assert n7 == new, "a refused re-application still altered the text"

    print("[ok] self-test passed: the four GENOME-WIDE rows are gone and the "
          "PER_COHORT rows are untouched; no executable GENOME_WIDE reference "
          "survives anywhere while the full provenance remains as comments, "
          "including the three-quantity table and the source document's "
          "banner; the refusal names two groups and never mentions "
          "S4_GENOME_WIDE; passing S4_GENOME_WIDE is now inert rather than an "
          "error; the two-quantities caveat and 139's derived gap are both "
          "removed; refuses without 138 or 139; double application refused")
    return 0


if __name__ == "__main__":
    sys.exit(main())
