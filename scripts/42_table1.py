#!/usr/bin/env python3
"""
42_table1.py -- TABLE 1: cohorts, with the true post-exclusion counts.

Writes Table1.tsv (machine-readable), Table1.md (paste into the manuscript) and prints
the footnotes. Then VERIFIES every count against the sample lists on disk and says
loudly where they disagree.

TWO DECISIONS BAKED IN, BOTH ARGUABLE -- change them here rather than in the output
-----------------------------------------------------------------------------------
1. LBNL IS IN THE MAIN TABLE, not a supplement, showing 5 DCM / 4 NF with the
   pre-exclusion 15/18 in parentheses and a footnote naming the overlap. Hiding a
   cohort that was in earlier drafts invites the question of what else moved; showing
   it shrink is the honest version and costs one footnote.
2. HEINIG IS NOT LISTED, and neither are Amsterdam or Naples. All three were retrieved
   and genotyped; none was analysed, Heinig because it has no control group and so
   supports no ΔPSI contrast. Listing one unused dataset and not the others implies the
   others do not exist. The Methods names all three in one sentence with the reason,
   which is the honest place for them. This also keeps Table 1 and Figure 1 identical:
   43_figure1_design.py draws only analysed cohorts.
3. LAYOUT AND READ LENGTH GET A COLUMN. Two orderings have been tested against this
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

WHAT THIS TABLE MUST NOT IMPLY
------------------------------
That the cohorts are comparable in effect size. They are not -- see
RESULT_magnitude_not_estimable.md and the Phase 1.5 decomposition. The table reports
composition and platform, never a per-cohort event count; those belong in Figure 4 with
their caveats attached.

Usage:
  export PATH=/data/apps/anaconda3/envs/mamba/envs/rnaseq_tools/bin:$PATH
  cd /data/project/chureeratp/RNA/HEART
  python 42_table1.py
  python 42_table1.py --no-verify
"""
import argparse
import os
import re

import pandas as pd

# ---------------------------------------------------------------------------
# FROZEN. Sources: PROJECT_STATE_handoff.md (accessions, platform, annotation),
# RESULT_sample_identity_LBNL_overlap.md (post-exclusion counts), PAPER_ONE_PAGE.md.
# `pre` is the count BEFORE the genotype exclusions; None where nothing was removed.
# ---------------------------------------------------------------------------
COHORTS = [
    dict(cohort="MAGNet", acc="PRJNA595151", tissue="LV, explant / donor",
         cases="166 DCM, 28 HCM", controls="166 NF", pre=None,
         layout="paired, 125 bp", prep="poly-A", annot="v49",
         role="Discovery. Only cohort with LVEF, RIN, age, sex and pool",
         lists=dict(dcm="dcm_bams.txt", nf="nf_bams.txt")),
    dict(cohort="Colorado", acc="PRJNA477855", tissue="LV, end-stage explant",
         cases="36 DCM, 13 ICM", controls="13 NF", pre="37 DCM, 13 ICM / 14 NF",
         layout="single, 51 bp", prep="intermediate", annot="v49",
         role="Independent DCM replication; second etiology pair. MXE excluded",
         lists=dict()),
    dict(cohort="LBNL", acc="PRJNA522417", tissue="LV, end-stage explant",
         cases="5 DCM", controls="4 NF", pre="15 DCM / 18 NF",
         layout="single, 50–51 bp", prep="total-RNA", annot="v49",
         role="Library-prep robustness check — NOT an independent cohort",
         lists=dict()),
    dict(cohort="Brigham", acc="PRJNA675090", tissue="Anterior septum, myectomy",
         cases="18 HCM", controls="5", pre=None,
         layout="paired, 101 bp", prep="—", annot="v49",
         role="Independent HCM", lists=dict()),
    dict(cohort="Fuwai", acc="PRJNA522931", tissue="Myocardium, myectomy",
         cases="28 HCM", controls="9 donor", pre=None,
         layout="paired, 150 bp", prep="—", annot="v46",
         role="Independent HCM, best powered", lists=dict()),
    dict(cohort="Stanford", acc="PRJNA747806", tissue="Septum, myectomy",
         cases="13 HCM", controls="6 donor", pre="13 HCM / 7 donor",
         layout="paired, 150 bp", prep="—", annot="v49",
         role="Independent HCM", lists=dict()),
]

FOOTNOTES = [
    ("a", "Counts are those analysed, after the genotype-identity exclusions described "
          "in the Methods. The pre-exclusion column gives the numbers as deposited, "
          "and is shown so that the difference is visible rather than silent."),
    ("b", "Twenty-four of LBNL's 33 samples are the same individuals as MAGNet samples "
          "(KING kinship 0.425–0.486, IBS0 ≈ 0, all 24 matching on disease group). "
          "LBNL is therefore not independent of the discovery cohort and is used only "
          "as a library-preparation robustness check; it is excluded from all pooled "
          "and combined-evidence analyses."),
    ("c", "Three within-cohort technical replicates were identified by the same "
          "analysis and one of each pair removed: two in Colorado (one DCM, one NF) "
          "and one in Stanford (control)."),
    ("d", "Fuwai was quantified against GENCODE v46; all other cohorts against v49. "
          "Mutually exclusive exon (MXE) calls from Colorado are excluded throughout "
          "(they do not reproduce: reverse Deming slope 0.061, n = 1,477, measured "
          "before the 125 bp re-quantification and not re-measured since), and "
          "Fuwai MXE is flagged for caution."),
    ("e", "MAGNet DCM and HCM contrasts share the same 166 non-failing hearts. Where "
          "that matters, the shared-control inflation is quantified rather than "
          "assumed (Methods)."),
    ("f", "Read layout and length were measured from the aligned reads rather than "
          "taken from the repository metadata. PSI is normalised by inclusion- and "
          "skipping-form lengths computed from the read length supplied to rMATS, so "
          "raw PSI is not comparable across cohorts and only ΔPSI direction is "
          "compared. MAGNet was originally quantified with --readLength 50 although "
          "its reads are 125 bp. The analyses reported here use a full re-quantification "
          "of all 332 samples at the sequenced 125 bp, differing from the original in "
          "--readLength and nothing else. Skipped-exon ΔPSI agreed in direction for "
          "99.89% of events (Spearman ρ = 0.996, orthogonal slope 0.999) and for "
          "100.00% of events significant in both arms, because PSI depends on the "
          "ratio of the two form lengths and that ratio moved 2% while the lengths "
          "themselves more than doubled. Counts are more sensitive than direction: "
          "skipped exons reaching FDR < 0.05 rose in both disease-versus-control "
          "contrasts and fell in the null contrast, while the subset also clearing "
          "|ΔPSI| ≥ 0.10 fell from 926 to 813, because the proportion of "
          "FDR-significant events clearing that threshold dropped from 17.5% to "
          "15.1% — a third of the significant set lies within 0.02 of the cut. "
          "Every replication claim here therefore rests on direction rather than on "
          "counts (Methods). LBNL read length is heterogeneous; "
          "--variable-read-length was set for every cohort."),
]


def ids(path):
    out = []
    for p in re.split(r"[,\n]", open(path).read()):
        p = p.strip()
        if not p:
            continue
        m = re.search(r"(SRR\d+|ERR\d+|SAMEA\d+|SAMN\d+)", p)
        out.append(m.group(1) if m else os.path.basename(os.path.dirname(p)))
    return sorted(set(out))


def nfrom(txt):
    """'166 DCM, 28 HCM' -> 194 ; '—'/'none' -> 0"""
    n = re.findall(r"(\d+)", txt or "")
    return sum(int(x) for x in n)


def verify(rows):
    """Check what can be checked against the files on disk."""
    print("\n" + "=" * 78)
    print("VERIFICATION against the sample lists on disk")
    print("=" * 78)
    problems = 0

    # MAGNet discovery counts come straight off the bam lists
    for r in rows:
        for key, path in r.get("lists", {}).items():
            if not os.path.exists(path):
                print(f"  [skip] {r['cohort']}: {path} not found")
                continue
            n = len(ids(path))
            want = {"dcm": 166, "nf": 166}.get(key)
            ok = "ok " if (want is None or n == want) else "***"
            if want is not None and n != want:
                problems += 1
            print(f"  [{ok}] {r['cohort']:<10}{path:<20}{n:>5} samples"
                  + (f"   (table says {want})" if want is not None else ""))

    # the exclusion list is the single source of truth for what was removed
    exc = "geno/EXCLUDE_samples.tsv"
    if os.path.exists(exc) and os.path.getsize(exc):
        E = pd.read_csv(exc, sep="\t")
        print(f"\n  [ok ] {exc}: {len(E)} samples excluded (table's footnotes "
              f"describe 24 + 3 = 27)")
        if len(E) != 27:
            print(f"  *** the doc says 27 exclusions, this file has {len(E)}. "
                  f"One of them is out of date.")
            problems += 1
        col = next((c for c in E.columns if c.lower() in ("cohort", "project")), None)
        if col:
            print("       by cohort: "
                  + ", ".join(f"{k} {v}" for k, v in E[col].value_counts().items()))
    else:
        print(f"\n  [skip] {exc} not found — cannot verify the exclusion count")

    # dedup id lists, where they exist, should match the table's post-exclusion n
    for f in sorted(os.listdir(".")):
        if f.endswith("_dedup.ids") and os.path.getsize(f):
            print(f"  [info] {f:<34}{len(ids(f)):>5} samples")

    print(f"\n  {problems} disagreement(s) between the table and the files on disk.")
    if problems:
        print("  Resolve them before submission: the table is what a reader checks "
              "first.")
    return problems


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="Table1")
    ap.add_argument("--no-verify", dest="ver", action="store_false")
    a = ap.parse_args()

    rows = []
    for c in COHORTS:
        # ASCII footnote markers, NOT Unicode superscript letters: those glyphs are
        # missing from Arial and Calibri, so Excel and Word substitute a fallback font
        # and the cell reads as a box. Real superscripts belong in the .docx, where
        # they are formatting rather than characters.
        cases = c["cases"]
        controls = c["controls"]
        if c["pre"]:
            cases += " (a)"          # superscript a, counts are post-exclusion
        mark = ""
        if c["cohort"] == "LBNL":
            mark = " (b)"
        if c["cohort"] in ("Colorado", "Stanford"):
            mark = " (c)"
        if c["cohort"] == "Fuwai":
            mark = " (d)"
        if c["cohort"] == "MAGNet":
            mark = " (e)"
        rows.append({
            "Cohort": c["cohort"] + mark,
            "BioProject": c["acc"],
            "Tissue": c["tissue"],
            "Cases (n)": cases,
            "Controls (n)": controls,
            "Pre-exclusion": c["pre"] or "—",
            "Layout, read length (f)": c["layout"],
            "Library prep": c["prep"],
            "GENCODE": c["annot"],
            "Role in this study": c["role"],
        })
    T = pd.DataFrame(rows)
    T.to_csv(f"{a.out}.tsv", sep="\t", index=False)

    with open(f"{a.out}.md", "w") as fh:
        fh.write("**Table 1. Cohorts.** Sample numbers are those analysed, after the "
                 "genotype-identity exclusions.\n\n")
        fh.write("| " + " | ".join(T.columns) + " |\n")
        fh.write("|" + "|".join(["---"] * len(T.columns)) + "|\n")
        for _, r in T.iterrows():
            fh.write("| " + " | ".join(str(x) for x in r.values) + " |\n")
        fh.write("\n")
        for k, v in FOOTNOTES:
            fh.write(f"^{k}^ {v}\n\n")
    print(f"[write] {a.out}.tsv and {a.out}.md")

    analysed = sum(nfrom(c["cases"]) + nfrom(c["controls"]) for c in COHORTS)
    print(f"\n  samples analysed across the {len(COHORTS)} cohorts with controls: "
          f"{analysed}")
    print("  Heinig (PRJEB8360, 128 DCM), Amsterdam and Naples were retrieved and")
    print("  genotyped but are NOT in this table: none was analysed, Heinig because it")
    print("  has no control group. The Methods names all three. Do not add one back")
    print("  without adding the others.")
    print("  NOTE: the genotype analysis covered 670 samples, which also included "
          "cohorts\n        not analysed here (Amsterdam, Naples). Do not present 670 "
          "as the study n.")

    print("\n" + "=" * 78)
    print("FOOTNOTES")
    print("=" * 78)
    for k, v in FOOTNOTES:
        print(f"  {k}) {v}")

    if a.ver:
        verify(COHORTS)

    print("\n  Check before submission:")
    print("   * Heinig is NOT listed. If it ever contributes to an analysis that")
    print("     survives into the paper, it needs a row — and so do Amsterdam and")
    print("     Naples if the same becomes true of them.")
    print("   * Brigham's 5 controls and the '—' library preps are the only cells with")
    print("     no source in the project docs. Fill them from the BioProject metadata")
    print("     or say 'not reported' rather than leaving an em dash.")
    print("   * Figure 1 is built from this table; rebuild it if any count changes.")
    print("   * MAGNet reads are 125 bp (measured, 49_readlength_audit.sh), NOT the")
    print("     50 bp earlier drafts claimed. rMATS was given --readLength 50; all")
    print("     332 samples have since been re-quantified at 125 bp and that is")
    print("     the primary analysis. Footnote (f) carries the comparison. Do not")
    print("     quietly revert this cell.")
    print("   * BOOKED MEASUREMENT, footnote (d): the Colorado-MXE exclusion rests")
    print("     on a reverse Deming slope of 0.061 (n = 1,477) measured at the")
    print("     ORIGINAL declared read length. The exclusion is applied throughout")
    print("     the 125 bp analysis. Re-measure it, or say in the Methods that it")
    print("     was not re-measured. Do not restate 0.061 as a 125 bp number.")
    print("   * The slope-versus-read-length-DISTANCE justification was retired on")
    print("     21 Oct (rho +0.316, p = 0.68 against MAGNet's true 125 bp). Item 3")
    print("     of this file's header records why. Do not reinstate it, and do not")
    print("     state the reframed own-length ordering (p = 0.051, n = 4) as")
    print("     established.")


if __name__ == "__main__":
    main()
