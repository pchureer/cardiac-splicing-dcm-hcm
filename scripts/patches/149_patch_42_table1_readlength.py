#!/usr/bin/env python3
"""
149_patch_42_table1_readlength.py -- Table 1's read-length prose.

THREE PROBLEMS, ALL IN THE TEXT. The data cells are already right: MAGNet's
row reads `layout="paired, 125 bp"`.

1. DOCSTRING ITEM 3 CONTRADICTS THE FILE'S OWN SUBMISSION CHECK.
   Item 3 justifies the layout column by "Figure 4's slopes order
   monotonically with read-length distance from MAGNet's paired 50bp". The
   check at the end of the same file says MAGNet's reads are 125 bp and that
   50 was only the value passed to rMATS. So the justification's reference
   point is not a read length at all.

   Measured against MAGNet's true 125 bp the distances become Colorado 74,
   Brigham 24, Stanford 25, Fuwai 25 and the ordering collapses: Spearman
   rho +0.316, p = 0.68, with the FURTHEST cohort transferring best. At the
   old declared 50 it was rho -0.949, p = 0.051.

   The slope ordering itself is real and invariant -- but against each
   cohort's OWN declared length, which never changed. Only the
   distance-from-MAGNet framing needed MAGNet's 50. Item 3 is rewritten to
   state the retired claim, the reframed one, and that p = 0.051 on four
   cohorts is suggestive rather than significant.

2. FOOTNOTE (f) QUOTES THE PROBE AS IF IT WERE THE COHORT.
   It reports "99.9% of events (Spearman rho = 0.997) and the
   significant-event count by -2.9%" from "30 case and 30 control samples".
   Those are the pilot figures -- the same 172 -> 167 that Supplementary S5(b)
   reports, and which 62 itself warns is "roughly a fifth of the discovery
   set". The full re-quantification of all 332 samples gives rho = 0.996,
   direction 99.89%, and 926 -> 813, i.e. -12.2%: four times the effect the
   footnote states. The footnote is also framed as though 50 bp were primary
   and 125 bp a check, which the promotion reverses.

   The replacement is the text already drafted in
   NUMBERS_125bp_REDERIVATION.md section 9 ("the read-length row is REPLACED,
   not updated"), rendered as a footnote.

3. FOOTNOTE (d) CARRIES AN UN-RE-MEASURED 50 bp NUMBER.
   "reverse Deming slope 0.061, n = 1,477" justifies excluding Colorado MXE
   THROUGHOUT, including throughout the 125 bp analysis, and has never been
   re-measured. This patch does NOT invent a value: it states the basis in
   the footnote and books the measurement in the submission check. A
   correctly-reported 50 bp number is not wrong; silently reusing it as if
   read length were irrelevant is.

A GUARD IS ADDED TOO
--------------------
verify() refuses if MAGNet's layout cell is not "paired, 125 bp". The file's
own check says "Do not quietly revert this cell"; now something enforces it.

USAGE
-----
  python 149_patch_42_table1_readlength.py [--apply] [--self-test]
"""
import argparse
import difflib
import os
import shutil
import subprocess
import sys
import time

TARGET = "42_table1.py"

# ------------------------------------------------------------- 1. item 3
A_ITEM3 = """3. LAYOUT AND READ LENGTH GET A COLUMN. Figure 4's slopes order monotonically with
   read-length distance from MAGNet's paired 50bp (Colorado 0.52 at single 51bp,
   Brigham 0.48 at paired 101bp, Stanford 0.39 and Fuwai 0.24 at paired 150bp). Without
   this column that ordering is an assertion; with it the reader can check it.
"""

B_ITEM3 = """3. LAYOUT AND READ LENGTH GET A COLUMN. Two orderings have been tested against this
   column and they do NOT fare the same way, so the justification is stated precisely
   rather than left as "the slopes order with read length".

   RETIRED 21 Oct -- "the slopes order with read-length DISTANCE from MAGNet's paired
   50bp". MAGNet's reads are 125 bp; 50 was the value passed to rMATS and was never a
   read length, as the submission check at the end of this file already said. Against
   MAGNet's true 125 the distances become Colorado 74, Brigham 24, Stanford 25,
   Fuwai 25 and the ordering collapses -- Spearman rho +0.316, p = 0.68, with the
   FURTHEST cohort transferring best. At the old declared 50 it was rho -0.949,
   p = 0.051. DO NOT REINSTATE IT.

   KEPT, reframed -- the slopes order with each REPLICATION COHORT'S OWN declared
   length, which never changed: Colorado 0.600 at single 51bp, Brigham 0.510 at paired
   101bp, Stanford 0.401 and Fuwai 0.295 at paired 150bp (125 bp primary; 0.521 /
   0.476 / 0.386 / 0.243 at 50 bp, the same order). rho -0.949 at BOTH read lengths,
   because neither the cohorts' lengths nor the slope ordering moved.

   BUT p = 0.051 on four cohorts with two tied at 150 bp. SUGGESTIVE, NOT SIGNIFICANT.
   The column lets a reader check the ordering; the text must not state it as
   established, and must not imply a mechanism it has not tested.

   (The related finding that DOES carry a p-value and a control -- JCEC's power
   advantage over JC falling with declared read length, rho -0.975, p = 0.0048, with
   MAGNet's own re-declaration moving it from +66.7% to +33.3% on the same samples --
   belongs to Supplementary S5(a), not to this table.)
"""

# ------------------------------------------------------- 2. footnote (f)
A_FN_F = (
    '          "compared. MAGNet was quantified with --readLength 50 although its '
    'reads are "\n'
    '          "125 bp; a matched re-quantification at 125 bp on 30 case and 30 '
    'control "\n'
    '          "samples changed skipped-exon ΔPSI direction in 99.9% of events "\n'
    '          "(Spearman ρ = 0.997) and the significant-event count by '
    '−2.9%, because PSI "\n'
    '          "depends on the ratio of the two form lengths and that ratio is '
    'nearly "\n'
    '          "invariant to read length (Methods). LBNL read length is '
    'heterogeneous; "\n'
)

B_FN_F = (
    '          "compared. MAGNet was originally quantified with --readLength 50 '
    'although "\n'
    '          "its reads are 125 bp. The analyses reported here use a full '
    're-quantification "\n'
    '          "of all 332 samples at the sequenced 125 bp, differing from the '
    'original in "\n'
    '          "--readLength and nothing else. Skipped-exon ΔPSI agreed in '
    'direction for "\n'
    '          "99.89% of events (Spearman ρ = 0.996, orthogonal slope 0.999) '
    'and for "\n'
    '          "100.00% of events significant in both arms, because PSI depends on '
    'the "\n'
    '          "ratio of the two form lengths and that ratio moved 2% while the '
    'lengths "\n'
    '          "themselves more than doubled. Counts are more sensitive than '
    'direction: "\n'
    '          "skipped exons reaching FDR < 0.05 rose in both disease-versus-control '
    '"\n'
    '          "contrasts and fell in the null contrast, while the subset also '
    'clearing "\n'
    '          "|ΔPSI| ≥ 0.10 fell from 926 to 813, because the proportion '
    'of "\n'
    '          "FDR-significant events clearing that threshold dropped from 17.5% to '
    '"\n'
    '          "15.1% — a third of the significant set lies within 0.02 of the '
    'cut. "\n'
    '          "Every replication claim here therefore rests on direction rather than '
    'on "\n'
    '          "counts (Methods). LBNL read length is heterogeneous; "\n'
)

# ------------------------------------------------------- 3. footnote (d)
A_FN_D = (
    '          "Mutually exclusive exon (MXE) calls from Colorado are excluded '
    'throughout "\n'
    '          "(they do not reproduce: reverse Deming slope 0.061, n = 1,477), and '
    'Fuwai "\n'
    '          "MXE is flagged for caution."),\n'
)

B_FN_D = (
    '          "Mutually exclusive exon (MXE) calls from Colorado are excluded '
    'throughout "\n'
    '          "(they do not reproduce: reverse Deming slope 0.061, n = 1,477, '
    'measured "\n'
    '          "before the 125 bp re-quantification and not re-measured since), and '
    '"\n'
    '          "Fuwai MXE is flagged for caution."),\n'
)

# --------------------------------------------------- 4. submission check
A_CHECK = (
    '    print("   * MAGNet reads are 125 bp (measured, 49_readlength_audit.sh), '
    'NOT the")\n'
    '    print("     50 bp earlier drafts claimed. rMATS was given --readLength 50; '
    'the")\n'
    '    print("     consequence was measured and is in footnote (f). Do not quietly '
    'revert")\n'
    '    print("     this cell.")\n'
)

B_CHECK = (
    '    print("   * MAGNet reads are 125 bp (measured, 49_readlength_audit.sh), '
    'NOT the")\n'
    '    print("     50 bp earlier drafts claimed. rMATS was given --readLength 50; '
    'all")\n'
    '    print("     332 samples have since been re-quantified at 125 bp and that '
    'is")\n'
    '    print("     the primary analysis. Footnote (f) carries the comparison. Do '
    'not")\n'
    '    print("     quietly revert this cell.")\n'
    '    print("   * BOOKED MEASUREMENT, footnote (d): the Colorado-MXE exclusion '
    'rests")\n'
    '    print("     on a reverse Deming slope of 0.061 (n = 1,477) measured at the")\n'
    '    print("     ORIGINAL declared read length. The exclusion is applied '
    'throughout")\n'
    '    print("     the 125 bp analysis. Re-measure it, or say in the Methods that '
    'it")\n'
    '    print("     was not re-measured. Do not restate 0.061 as a 125 bp '
    'number.")\n'
    '    print("   * The slope-versus-read-length-DISTANCE justification was retired '
    'on")\n'
    '    print("     21 Oct (rho +0.316, p = 0.68 against MAGNet\'s true 125 bp). '
    'Item 3")\n'
    '    print("     of this file\'s header records why. Do not reinstate it, and do '
    'not")\n'
    '    print("     state the reframed own-length ordering (p = 0.051, n = 4) as")\n'
    '    print("     established.")\n'
)

EDITS = [
    ("item 3's read-length justification", A_ITEM3, B_ITEM3, 1),
    ("footnote (f)", A_FN_F, B_FN_F, 1),
    ("footnote (d)", A_FN_D, B_FN_D, 1),
    ("the submission check", A_CHECK, B_CHECK, 1),
]


def apply_edits(text):
    if "RETIRED 21 Oct" in text:
        return text, ["42 already carries 149 (restore the .bak to re-apply)"]
    problems, out = [], text
    for name, anchor, repl, want in EDITS:
        got = out.count(anchor)
        if got != want:
            problems.append(
                f"{name}: anchor occurs {got}x, expected {want}x\n"
                f"      first line: {anchor.splitlines()[0][:72]}")
            continue
        out = out.replace(anchor, repl, want)
    return out, problems


def verify(text):
    bad = []
    for w in ("RETIRED 21 Oct",
              "rho +0.316, p = 0.68",
              "KEPT, reframed",
              "SUGGESTIVE, NOT SIGNIFICANT",
              "full re-quantification",
              "99.89% of events",
              "fell from 926 to 813",
              "17.5% to ",
              "not re-measured since",
              "BOOKED MEASUREMENT, footnote (d)"):
        if w not in text:
            bad.append(f"missing: {w[:58]}")
    # the superseded probe figures must be GONE from footnote (f)
    for w in ("30 case and 30 control",
              "99.9% of events",
              "ρ = 0.997",
              "the significant-event count by −2.9%"):
        if w in text:
            bad.append(f"the probe figure survives in footnote (f): {w!r}")
    # the retired justification must be gone
    if "read-length distance from MAGNet's paired 50bp (Colorado 0.52" in text:
        bad.append("the retired distance justification survives verbatim")
    # THE DATA CELL GUARD. The file's own check says do not revert it; this
    # is what enforces it.
    if 'layout="paired, 125 bp"' not in text:
        bad.append("MAGNet's layout cell is not 'paired, 125 bp' — the whole "
                   "footnote is written against a 125 bp primary analysis")
    # 0.061 must still be present IN FOOTNOTE (d), and still attributed.
    # Checking for a bare "0.061" anywhere is not enough: the submission
    # check this patch adds also names it, so the footnote could be gutted
    # while the loose check still passed.
    if "0.061, n = 1,477, measured" not in text:
        bad.append("footnote (d) no longer reads '0.061, n = 1,477, "
                   "measured ...' — it is a correctly reported 50 bp "
                   "measurement, to be attributed and not deleted")
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
    print("149 — Table 1's read-length prose")
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
        print()
        print("  The footnotes contain U+2212 MINUS, U+0394 DELTA, U+03C1 RHO")
        print("  and U+2265 GREATER-OR-EQUAL. If an anchor missed, print the")
        print("  exact bytes before editing by hand:")
        print("    python3 -c \"import sys;"
              "[print(i+1, repr(l)) for i,l in "
              "enumerate(open('42_table1.py').read().splitlines()) "
              "if '0.997' in l or '1,477' in l]\"")
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
        print("  python3 -c \"import ast; ast.parse(open('42_table1.py')"
              ".read())\" && echo SYNTAX_OK")
        print("  python3 42_table1.py 2>&1 | tail -40")
        print()
        print("⚠ THE DATA CELLS MUST NOT MOVE. This patch changes the")
        print("  docstring, two footnotes and the submission check only.")
        print("  Diff Table1.tsv against a copy taken first:")
        print()
        print("    cp Table1.tsv Table1_BEFORE_149.tsv")
        print("    python3 42_table1.py >/dev/null 2>&1")
        print("    diff Table1.tsv Table1_BEFORE_149.tsv && echo IDENTICAL")
        print()
        print("  Figure 1 is built from Table1.tsv, so IDENTICAL also means")
        print("  Figure 1 does not need rebuilding for this change.")
        print()
        print("⬜ STILL OWED, and this patch only books it: footnote (d)'s")
        print("  reverse Deming slope 0.061 (n = 1,477) justifies excluding")
        print("  Colorado MXE throughout, was measured before the 125 bp")
        print("  re-quantification, and has never been re-measured. It is now")
        print("  attributed in the footnote and listed in the submission")
        print("  check, which is disclosure, not a measurement.")
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
    print("✓ verified on disk: the distance justification is retired with its")
    print("  numbers, the own-length reframing is marked suggestive, footnote")
    print("  (f) carries the 332-sample figures instead of the 60-sample")
    print("  probe, footnote (d) states its basis, and MAGNet's 125 bp cell")
    print("  is now guarded rather than only asked about.")
    return 0


def self_test():
    import tempfile

    fx = ('"""\n'
          "A fixture standing in for 42_table1.py.\n\n"
          + A_ITEM3 +
          '"""\n'
          "COHORTS = [\n"
          '    dict(cohort="MAGNet", layout="paired, 125 bp"),\n'
          "]\n"
          "FOOTNOTES = [\n"
          '    ("d", "Fuwai was quantified against GENCODE v46; all other cohorts '
          'against v49. "\n'
          + A_FN_D +
          '    ("f", "Read layout and length were measured from the aligned reads '
          'rather than "\n'
          '          "taken from the repository metadata. PSI is normalised by '
          'inclusive- and "\n'
          + A_FN_F +
          '          "--variable-read-length was set for every cohort."),\n'
          "]\n"
          "\n"
          "def main():\n"
          '    print("   * rows:", COHORTS)\n'
          + A_CHECK +
          "    for k, v in FOOTNOTES:\n"
          '        print(f"FN {k}: {v}")\n'
          "main()\n")

    new, probs = apply_edits(fx)
    assert not probs, probs
    assert not verify(new), verify(new)
    compile(new, "fx", "exec")

    d = tempfile.mkdtemp()

    def run(src, name):
        p = os.path.join(d, name)
        open(p, "w", encoding="utf-8").write(src)
        r = subprocess.run([sys.executable, "-I", p], capture_output=True,
                           text=True, cwd=d)
        assert r.returncode == 0, r.stdout + r.stderr
        return r.stdout

    o0 = run(fx, "fx_pre.py")
    o1 = run(new, "fx_post.py")

    # the probe figures are gone from the rendered footnote
    assert "30 case and 30 control" in o0, o0
    assert "30 case and 30 control" not in o1, o1
    assert "ρ = 0.997" in o0 and "ρ = 0.997" not in o1, o1
    # and the full-cohort figures are there
    assert "99.89% of events" in o1, o1
    assert "ρ = 0.996" in o1, o1
    assert "fell from 926 to 813" in o1, o1
    # footnote (d) keeps 0.061 and gains its basis
    assert "0.061, n = 1,477, measured" in o1, o1
    assert "not re-measured since" in o1, o1
    # the submission check gains both new items
    assert "BOOKED MEASUREMENT, footnote (d)" in o1, o1
    assert "retired on" in o1, o1
    assert "rho +0.316, p = 0.68" in o1 or "+0.316" in o1, o1
    # the data row is untouched
    assert "'layout': 'paired, 125 bp'" in o1, o1

    # --- THE GUARD: reverting MAGNet's cell must be caught, because every
    #     footnote this patch writes is phrased for a 125 bp primary
    reverted = new.replace('layout="paired, 125 bp"', 'layout="paired, 50 bp"')
    v = verify(reverted)
    assert any("layout cell is not" in b for b in v), v

    # --- deleting 0.061 instead of attributing it must also be caught
    # Deleting 0.061 from the FOOTNOTE must be caught even though the
    # submission check this patch adds also mentions 0.061 -- an earlier
    # version of verify() looked for a bare "0.061" anywhere and so passed
    # a gutted footnote.
    gutted = new.replace("0.061, n = 1,477, measured", "n = 1,477, measured")
    assert "0.061" in gutted, "the fixture no longer exercises the loose case"
    v2 = verify(gutted)
    assert any("footnote (d) no longer reads" in b for b in v2), v2

    # --- a file still carrying the retired justification verbatim fails
    stale = new.replace("RETIRED 21 Oct", "retired") + A_ITEM3
    v3 = verify(stale)
    assert any("retired distance justification survives" in b for b in v3), v3

    # --- double application refused without altering the text
    n4, p4 = apply_edits(new)
    assert p4 and "already carries 149" in p4[0], p4
    assert n4 == new, "a refused re-application still altered the text"

    print("[ok] self-test passed: the 60-sample probe figures (30+30, 99.9%, "
          "rho 0.997, -2.9%) are gone from footnote (f) and the 332-sample "
          "ones (99.89%, rho 0.996, 926 -> 813) are in; footnote (d) keeps "
          "0.061 and states it predates the re-quantification; item 3 "
          "records the retired distance claim with its numbers and marks the "
          "reframed own-length ordering suggestive; the submission check "
          "books the MXE re-measurement; the MAGNet 125 bp cell is guarded, "
          "and reverting it, deleting 0.061, or leaving the retired text in "
          "place are each caught; double application refused")
    return 0


if __name__ == "__main__":
    sys.exit(main())
