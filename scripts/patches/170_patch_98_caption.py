#!/usr/bin/env python3
"""
170_patch_98_caption.py -- Table 4's caption still prints a pointer that was
recorded as wrong on 12 October.

98_build_replication_table.py line ~493 prints the caption. It is correct in
almost every respect and better than the replacement drafted for it -- it
carries the Deming lambda definition, the attenuation caveat, footnote b and
"verified from genotypes". One line is wrong:

    bound; the matched-replicate panel separates the technical component
    (Supplementary Table 4).

DECISION_supplement_numbering_and_s5d.md records this exact string: under
Mapping A, "Supplementary Table 4" is the all-event-class table, so the pointer
is WRONG TYPE and WRONG CONTENT. The correct target is Supplementary Figure 2.

> A defect recorded is not a defect fixed. The correct replacement was written
> down eleven days ago. Nothing re-ran the script to check, because the string
> lives in a print() and no gate reads a caption.

WHY THE WHOLE TEMPLATE, NOT JUST THE POINTER
--------------------------------------------
98's caption is 236 words and the manuscript needs ~208 of them. Fixing only
the pointer leaves the script printing one version while the paper carries
another -- the two-live-versions problem this project has hit repeatedly. The
script should print what goes in the paper.

WHAT CHANGES
------------
  1. (Supplementary Table 4) -> (Supplementary Figure 2)          [the defect]
  2. + "Concordance and slope do not rank together."              [Brigham has
     the highest concordance, Colorado the highest slope, at BOTH read lengths]
  3. + "Slopes are point estimates without confidence intervals." [98's own
     closing item 3: "say so, or compute them". Nothing had said so.]
  4. title: "three independent cohorts and a second etiology pair" -> "four
     independent cohorts", matching the manuscript's own table title. The 3 + 1
     framing reads as a contradiction of it.
  5. the pooled sentence is removed -- it is the table's own Pooled row, and
     the Results states it.
  6. the read-length normalisation clause is shortened; it is in the Methods'
     Read length subsection and in Supplementary Methods M6.

The lambda definition and both footnotes are KEPT. Every anchor below is pure
ASCII, deliberately: the caption also contains lambda, superscripts, dPSI and
superscript footnote letters, and an anchor spanning those is a transcription
risk for no benefit.

USAGE
-----
  python3 170_patch_98_caption.py [--apply] [--self-test]
"""
import argparse
import difflib
import os
import shutil
import sys
import time

TARGET = "98_build_replication_table.py"

EDITS = [
    ("the stale supplementary pointer",
     """bound; the matched-replicate panel separates the technical component
(Supplementary Table 4).""",
     """bound; the matched-replicate panel separates the technical component
(Supplementary Figure 2). Concordance and slope do not rank together."""),

    ("the title's cohort count",
     """**Table 4. The discovery signature replicates by direction in three independent cohorts
and a second etiology pair.** Skipped exons significant in the discovery contrast""",
     """**Table 4. The discovery signature replicates by direction in four independent
cohorts.** Skipped exons significant in the discovery contrast"""),

    ("the slope confidence-interval statement",
     "the discovery cohort on the x axis. Group sizes are post-exclusion, as in Table 1.",
     "the discovery cohort on the x axis. Slopes are point estimates without confidence\n"
     "intervals. Group sizes are post-exclusion, as in Table 1."),

    ("the duplicated pooled sentence",
     "\nPooled concordance is {K:,} of {N:,} ({100 * K / N:.1f}%, P = {pv:.1e}).\n",
     ""),

    ("the read-length clause",
     """**Neither effect-size magnitude nor event counts are compared across cohorts**, here or
anywhere in this paper: PSI is normalised by an inclusion-form length that depends on
read length, and the cohorts differ in tested event set, depth, covariate model and""",
     """**Neither effect-size magnitude nor event counts are compared across cohorts**
(Methods): the cohorts differ in read length, sequencing depth, covariate model and"""),
]


def apply_edits(text):
    if "170:" in text:
        return text, ["98 already carries 170 (restore the .bak to re-apply)"]
    probs = [f"{n}: anchor occurs {text.count(a)}x, expected 1x"
             for n, a, b in EDITS if text.count(a) != 1]
    if probs:
        return text, probs
    for n, a, b in EDITS:
        text = text.replace(a, b, 1)
    # a marker so re-application refuses, placed in the caption's own comment
    text = text.replace(
        "    # ------------------------------------------------------------------ caption",
        "    # ------------------------------------------------------------------ caption\n"
        "    # 170: pointer corrected to Supplementary Figure 2 (Mapping A), slope\n"
        "    # point-estimate caveat added, pooled sentence dropped as duplicated.",
        1)
    return text, []


def verify(text):
    bad = []
    for w in ("(Supplementary Figure 2). Concordance and slope do not rank together.",
              "Slopes are point estimates without confidence",
              "in four independent\ncohorts.**",
              "(Methods): the cohorts differ in read length, sequencing depth",
              "170: pointer corrected to Supplementary Figure 2"):
        if w not in text:
            bad.append(f"missing: {w[:52]}")
    for w in ("(Supplementary Table 4)",
              "and a second etiology pair",
              "Pooled concordance is {K:,}",
              "here or\nanywhere in this paper"):
        if w in text:
            bad.append(f"survived, should be gone: {w[:46]}")
    # the parts that must NOT be touched
    for w in ("Supplementary Table 8", "verified from genotypes",
              "errors-in-variables", "attenuated by selection on discovery"):
        if w not in text:
            bad.append(f"content that should be KEPT was lost: {w[:44]}")
    return bad


FIX = '''    # ------------------------------------------------------------------ caption
    hr("=")
    print(f"""
**Table 4. The discovery signature replicates by direction in three independent cohorts
and a second etiology pair.** Skipped exons significant in the discovery contrast
(FDR < 0.05) and measurable in each replication cohort, exon-level, JCEC
quantification. Concordance is the proportion changing in the same direction as discovery,
tested against 50% with a one-sided exact binomial test; intervals are Wilson 95%
intervals. The transfer slope is a Deming (errors-in-variables) regression of the
replication cohort on the discovery cohort, with lambda as recorded and
the discovery cohort on the x axis. Group sizes are post-exclusion, as in Table 1.

Pooled concordance is {K:,} of {N:,} ({100 * K / N:.1f}%, P = {pv:.1e}).

**Neither effect-size magnitude nor event counts are compared across cohorts**, here or
anywhere in this paper: PSI is normalised by an inclusion-form length that depends on
read length, and the cohorts differ in tested event set, depth, covariate model and
power. Every slope is attenuated by selection on discovery significance and is a lower
bound; the matched-replicate panel separates the technical component
(Supplementary Table 4).

a LBNL is excluded from this table and from every pooled analysis: 24 of its 33 samples
  are the same individuals as discovery-cohort samples, verified from genotypes, so it is
  a library-preparation replicate rather than an independent cohort.
b Mutually exclusive exons are excluded for Colorado and Fuwai on reverse-replication
  grounds, on different evidence in each case (Supplementary Table 8).
""")
'''


def self_test():
    f = []
    new, probs = apply_edits(FIX)
    if probs:
        f.append(f"fixture did not patch: {probs}")
    else:
        bad = verify(new)
        if bad:
            f.append(f"verify rejected its own output: {bad}")
        # the f-string must still be a valid f-string after the pooled line goes
        body = new[new.index('print(f"""'):]
        body = body[len('print(f"""'):body.index('""")')]
        try:
            eval('f"""' + body.replace("{K:,}", "{0}") + '"""', {}, {})
        except Exception as e:
            f.append(f"the caption is no longer a valid f-string: {e}")
        # K/N/pv are gone, so the remaining braces must be none
        if "{" in body:
            f.append(f"a format placeholder survived the pooled-sentence removal: "
                     f"{body[body.index('{'):body.index('{') + 30]!r}")
        # the whole file must still compile
        try:
            compile("def _f():\n" + "\n".join("    " + l for l in new.splitlines()),
                    "x", "exec")
        except SyntaxError as e:
            f.append(f"does not compile: {e}")
        again, p2 = apply_edits(new)
        if not p2 or "already carries 170" not in p2[0]:
            f.append("re-application did not refuse")

    print("170 self-test")
    for x in f:
        print("  ✗ " + x)
    if not f:
        print("  ✓ the pointer now reads Supplementary Figure 2; the slope")
        print("    point-estimate caveat and the ranking sentence are in; the title")
        print("    says four independent cohorts; the duplicated pooled sentence is")
        print("    gone and NO format placeholder survives it; lambda, both")
        print("    footnotes, 'verified from genotypes' and the attenuation caveat")
        print("    are untouched; the caption is still a valid f-string; re-apply refuses")
    return 1 if f else 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("target", nargs="?", default=TARGET)
    ap.add_argument("--apply", action="store_true")
    ap.add_argument("--self-test", action="store_true")
    a = ap.parse_args()
    if a.self_test:
        return self_test()
    if self_test():
        print("\nREFUSING: 170's own self-test failed.")
        return 3
    print()
    if not os.path.exists(a.target):
        print(f"REFUSING: {a.target} not found")
        return 3
    src = open(a.target, errors="replace").read()
    new, probs = apply_edits(src)
    print("=" * 74)
    print("170 — Table 4's caption")
    print("=" * 74)
    if probs:
        print("REFUSING:")
        for p in probs:
            print("    " + p)
        print("\n  Every anchor is pure ASCII and copied from the script's own")
        print("  output. A miss means the file differs from what was read; paste")
        print("  the region rather than relaxing the anchor.")
        return 3
    print("".join(difflib.unified_diff(src.splitlines(True), new.splitlines(True),
                                       fromfile=a.target,
                                       tofile=a.target + " (patched)", n=3)))
    if not a.apply:
        print("DRY RUN. Re-run with --apply, then print the caption:")
        print()
        print("  TBL_SUFFIX=_rl125 \\")
        print("  REC_KN=Colorado=311/385,Brigham=298/364,Stanford=345/455,Fuwai=422/596 \\")
        print("  REC_SLOPE=Colorado=0.5999,Brigham=0.5099,Stanford=0.4008,Fuwai=0.2948 \\")
        print("  REC_POOLED=1376,1800 REC_POOLED_P=2.51e-117 \\")
        print("  REC_MEAN_SLOPE=0.451 REC_LBNL=0.5726/78.2 \\")
        print("  python3 98_build_replication_table.py --apply")
        print()
        print("  ⚠ --apply is needed: the caption prints after the dry-run return,")
        print("    which is why --list showed none. It rewrites Table4_rl125.tsv")
        print("    and .md with the same content as 17 Oct.")
        print()
        print("  Then check the caption says Supplementary FIGURE 2, names FOUR")
        print("  cohorts, and has no 'Pooled concordance is' line.")
        return 0
    bak = f"{a.target}.bak.{time.strftime('%Y%m%d_%H%M%S')}"
    shutil.copy2(a.target, bak)
    open(a.target, "w").write(new)
    print(f"[write] {a.target}   (backup: {bak})")
    disk = open(a.target, errors="replace").read()
    bad = verify(disk)
    try:
        compile(disk, a.target, "exec")
    except SyntaxError as e:
        bad.append(f"does not compile: {e}")
    if bad or disk != new:
        print("✗ ON-DISK VERIFICATION FAILED — restoring:")
        for b in bad:
            print(f"    {b}")
        shutil.copy2(bak, a.target)
        return 3
    print("✓ verified: the caption points at Supplementary Figure 2, states the")
    print("  slopes are point estimates, and keeps lambda and both footnotes.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
