#!/usr/bin/env python3
"""
74_build_s8a_permutations.py -- Supplementary S8(a): the permutation nulls, deposited.

WHY THIS ONE PARSES A LOG
-------------------------
Every other builder in this project refuses to parse console output, because a log is a
transcription and recomputing from the tables is always available. **This one is the
exception, and deliberately so.** The permutation result cannot be recomputed without
re-running `25_colorado_inflation.sbatch` -- 2h08 on the cluster -- so the log IS the
primary record. That is precisely why the claim went missing in the first place: it was
published, measured, and never written to a file.

Because the log is the only source, the gates are stricter than usual: the script checks
the SHAPE of the run (how many shuffles, how many analyses, which seeds), the VALUE of
every permutation (all must be exactly zero), and one recorded figure it can cross-check
(MAGNet seed 101 at 103,054 events, which the Results quotes verbatim). It refuses rather
than depositing a partial or altered run.

WHAT THE RUN CONTAINS -- and the correction it forced
-----------------------------------------------------
    Colorado  5 label shuffles   3 x (DCM vs NF, seeds 101/102/103)
                                 2 x (ICM vs NF, seeds 101/102)
              x 2 models (P0 = no covariate, P2 = + log10_jcov)
              x 2 event classes (SE, MXE)
              = 20 analyses
    MAGNet    2 label shuffles   seeds 101/102, SE only
              = 2 analyses

**The Results says "thirteen label permutations in Colorado". It is FIVE shuffles / TWENTY
analyses.** MAGNet's "two permutations across 103,054 events" is exactly right. See
`RESULT_permutations_rerun.md`.

TWO CAVEATS THE TSV CARRIES PER ROW, BECAUSE THEY MATTER
--------------------------------------------------------
1. **Counting convention is FULL-COORDINATE, not exon-level.** Script 25's real-model row
   reads `S0 tested=25135 sig=478` against the manuscript's exon-level 476 of 19,953.
   Irrelevant to the permutations -- zero is zero under either convention, since
   collapsing can only hold or reduce a count -- but if the convention is not recorded,
   the next reader compares 478 with 476 and concludes something moved.
2. **Tested counts drift between shuffles**, because permuting the labels changes which
   events pass the coverage and PSI-variance filters. So "the same events were tested" is
   not true, and a single tested figure for all rows would be wrong. Hence one per row.

The real (non-permuted) models are deposited alongside with `run_kind = real`, because a
permutation null is only interpretable next to the positive control from the same script.

USAGE
-----
    conda activate rnaseq_tools
    cd /data/project/chureeratp/RNA/HEART

    python 74_build_s8a_permutations.py                  # default log, builds the TSV
    python 74_build_s8a_permutations.py --caption         # also print the S8(a) caption
"""
import argparse
import os
import re
import sys

import pandas as pd

LOG = "logs/permutation_rerun_406341.out"
ERR = "25_colorado_inflation.sbatch_406341.err"
JOB = "406341"
SCRIPT = "25_colorado_inflation.sbatch"
OUT = "S8a_permutations.tsv"

CONVENTION = "full_coordinate"   # a TOKEN, not prose -- the note lives in the caption
RULE = "FDR<0.05 and |dPSI|>=0.10"

# ------------------------------------------------------------------ the record
EXPECT_COLORADO_ANALYSES = 20
EXPECT_COLORADO_SHUFFLES = 5          # (contrast, seed) pairs
EXPECT_MAGNET_ANALYSES = 2
EXPECT_MAGNET_SEED101_TESTED = 103_054   # quoted verbatim in the Results
# script 25's ORIGINAL output, recorded in RESULT_colorado_mxe_excluded.md:
#   "adding log10_jcov moved SE from 478 to 475 and MXE from 1,596 to 1,605"
# If the re-run reproduces these, it is the same run the paper's claims came from.
EXPECT_REAL = {("colorado_dcmnf", "SE", 0): 478, ("colorado_dcmnf", "SE", 1): 475,
               ("colorado_dcmnf", "MXE", 0): 1596, ("colorado_dcmnf", "MXE", 1): 1605}

# `colorado_dcmnf_P0s101_SE` / `colorado_dcmnf_S0_SE`
LABEL = re.compile(
    r"^(?P<contrast>[a-z][a-z0-9_]*?)_(?P<kind>[SP])(?P<midx>\d+)"
    r"(?:s(?P<seed>\d+))?_(?P<et>SE|MXE|A3SS|A5SS|RI)$")
ROW = re.compile(
    r"^\s*(?P<label>\S+)\s+(?P<et>SE|MXE|A3SS|A5SS|RI)\s+"
    r"tested=(?P<tested>\d+)\s+sig=(?P<sig>\d+)\s*(?:model:\s*(?P<model>.*?))?\s*$")
MAGNET = re.compile(
    r"^\s*MAGNet\s+(?P<et>SE|MXE|A3SS|A5SS|RI)\s+permutation\s+seed\s+(?P<seed>\d+)\s+"
    r"tested=(?P<tested>\S+)\s+sig=(?P<sig>\S+)\s*$")


def parse(log):
    rows, skipped = [], []
    for ln, line in enumerate(open(log), 1):
        m = MAGNET.match(line)
        if m:
            t, s = m.group("tested"), m.group("sig")
            if not (t.isdigit() and s.isdigit()):
                skipped.append(f"line {ln}: MAGNet row with non-numeric "
                               f"tested={t} sig={s}")
                continue
            rows.append(dict(cohort="MAGNet", contrast="MAGNet_DCMvsNF_labels_shuffled",
                             event_type=m.group("et"), run_kind="permutation",
                             model_index="", model_formula="(script default)",
                             seed=int(m.group("seed")),
                             n_tested=int(t), n_significant=int(s), source_line=ln))
            continue
        m = ROW.match(line)
        if not m:
            continue
        lab = LABEL.match(m.group("label"))
        if not lab:
            skipped.append(f"line {ln}: unparsable label '{m.group('label')}'")
            continue
        kind = "permutation" if lab.group("kind") == "P" else "real"
        seed = lab.group("seed")
        if kind == "permutation" and seed is None:
            skipped.append(f"line {ln}: permutation row with no seed "
                           f"('{m.group('label')}')")
            continue
        rows.append(dict(cohort="Colorado", contrast=lab.group("contrast"),
                         event_type=m.group("et"), run_kind=kind,
                         model_index=int(lab.group("midx")),
                         model_formula=(m.group("model") or "").strip(),
                         seed=(int(seed) if seed else ""),
                         n_tested=int(m.group("tested")),
                         n_significant=int(m.group("sig")), source_line=ln))
    return pd.DataFrame(rows), skipped


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--log", default=LOG)
    ap.add_argument("--err", default=ERR)
    ap.add_argument("--out", default=OUT)
    ap.add_argument("--caption", action="store_true",
                    help="print the S8(a) caption text after building")
    a = ap.parse_args()

    print("\n" + "=" * 78)
    print("S8(a) -- PERMUTATION NULLS, deposited from the run log")
    print("=" * 78)
    if not os.path.exists(a.log):
        sys.exit(f"[fatal] {a.log} not found")
    print(f"\n  log: {a.log}")

    df, skipped = parse(a.log)
    if skipped:
        print(f"\n  {len(skipped)} line(s) matched a row pattern but could not be parsed:")
        for s in skipped[:10]:
            print(f"    * {s}")
        print("  REFUSING -- an unparsed row may be a permutation this deposit would")
        print("  silently omit. Fix the pattern, do not ignore them.\n")
        return

    perm = df[df.run_kind == "permutation"]
    real = df[df.run_kind == "real"]
    col_p = perm[perm.cohort == "Colorado"]
    mag_p = perm[perm.cohort == "MAGNet"]
    shuffles = col_p.groupby(["contrast", "seed"]).ngroups

    print(f"\n  parsed {len(df)} rows: {len(perm)} permutation, {len(real)} real")
    print(f"    Colorado permutations : {len(col_p)} analyses from {shuffles} shuffles")
    print(f"    MAGNet permutations   : {len(mag_p)} analyses")

    # --------------------------------------------------------------- the gates
    fail = []
    if len(col_p) != EXPECT_COLORADO_ANALYSES:
        fail.append(f"Colorado analyses {len(col_p)}, expected {EXPECT_COLORADO_ANALYSES}")
    if shuffles != EXPECT_COLORADO_SHUFFLES:
        fail.append(f"Colorado shuffles {shuffles}, expected {EXPECT_COLORADO_SHUFFLES}")
    if len(mag_p) != EXPECT_MAGNET_ANALYSES:
        fail.append(f"MAGNet analyses {len(mag_p)}, expected {EXPECT_MAGNET_ANALYSES}")

    nonzero = perm[perm.n_significant != 0]
    if len(nonzero):
        fail.append(f"{len(nonzero)} permutation(s) returned NON-ZERO")
        for _, r in nonzero.iterrows():
            fail.append(f"    {r.cohort} {r.contrast} seed {r.seed} {r.event_type}: "
                        f"sig={r.n_significant}")
    else:
        print(f"\n  [ok] all {len(perm)} permutations returned exactly 0")

    # does the re-run reproduce script 25's original real-model output?
    repro, missing = [], []
    for (ctr, et, mi), want in EXPECT_REAL.items():
        hit = real[(real.contrast == ctr) & (real.event_type == et) &
                   (real.model_index == mi)]
        if len(hit) != 1:
            missing.append(f"{ctr} {et} model {mi} not found once")
            continue
        got = int(hit.iloc[0].n_significant)
        repro.append((f"{ctr} {et} S{mi}", got, want, got == want))
    if repro:
        okn = sum(1 for _, _, _, o in repro if o)
        print(f"\n  REPRODUCTION against script 25's ORIGINAL output "
              f"(RESULT_colorado_mxe_excluded.md):")
        for lab, got, want, o in repro:
            print(f"    [{'ok     ' if o else 'MISMATCH'}] {lab:<24}{got:>6} "
                  f"against {want:>6}")
        if okn == len(repro):
            print(f"    All {okn} reproduce -- this is the same run the paper's claims")
            print(f"    came from, three weeks apart. The deposit rests on that.")
        else:
            fail.append(f"only {okn} of {len(repro)} real models reproduce script 25's "
                        f"original output -- the re-run is not the same analysis")
    if missing:
        print(f"    [note] {len(missing)} expected real model(s) absent: "
              f"{'; '.join(missing)}")

    m101 = mag_p[mag_p.seed == 101]
    if len(m101) == 1:
        got = int(m101.iloc[0].n_tested)
        if got != EXPECT_MAGNET_SEED101_TESTED:
            fail.append(f"MAGNet seed 101 tested={got:,}, the Results quotes "
                        f"{EXPECT_MAGNET_SEED101_TESTED:,}")
        else:
            print(f"  [ok] MAGNet seed 101 reproduces the quoted "
                  f"{EXPECT_MAGNET_SEED101_TESTED:,} tested events")
    else:
        fail.append("MAGNet seed 101 row not found exactly once")

    if len(real) and (real.n_significant <= 0).all():
        fail.append("no real (non-permuted) model returned any significant event -- "
                    "the positive control is missing, so the nulls are uninterpretable")
    elif len(real):
        print(f"  [ok] positive control present: real models return "
              f"{int(real.n_significant.min())}-{int(real.n_significant.max())} "
              f"significant")

    if os.path.exists(a.err):
        bad = sum(1 for l in open(a.err, errors="replace")
                  if re.search(r"error|traceback|exception", l, re.I))
        print(f"  [{'ok' if not bad else 'WARN'}] {a.err}: "
              f"{bad} error/traceback line(s)")
    else:
        print(f"  [note] {a.err} not found -- stderr not checked")

    if fail:
        print("\n" + "=" * 78)
        print("NOT DEPOSITED -- a gate failed")
        print("=" * 78)
        for f in fail:
            print(f"  * {f}")
        print("\n  Nothing written. A permutation deposit that is wrong is worse than")
        print("  none, because it looks like evidence.\n")
        return

    # --------------------------------------------------------------- deposit
    df = df.assign(job_id=JOB, script=SCRIPT, log_file=os.path.basename(a.log),
                   counting_convention=CONVENTION, significance_rule=RULE)
    # Colorado MXE is excluded from every table and figure in the manuscript
    # (reverse-replication slope 0.061, n=1,477). Its rows STAY in the deposit --
    # removing data from an artefact is worse than labelling it -- but they are
    # flagged, because an unflagged sig=2226 invites someone to quote it.
    df["excluded_from_manuscript"] = [
        "yes" if (c == "Colorado" and e == "MXE") else "no"
        for c, e in zip(df.cohort, df.event_type)]
    df["exclusion_reason"] = [
        "colorado_mxe_reverse_replication_slope_0.061" if x == "yes" else ""
        for x in df.excluded_from_manuscript]
    df["pct_significant"] = (100.0 * df.n_significant / df.n_tested).round(2)
    cols = ["job_id", "script", "log_file", "cohort", "contrast", "event_type",
            "run_kind", "model_index", "model_formula", "seed",
            "n_tested", "n_significant", "pct_significant",
            "excluded_from_manuscript", "exclusion_reason",
            "counting_convention", "significance_rule", "source_line"]
    df = df[cols].sort_values(["run_kind", "cohort", "contrast", "event_type",
                               "seed", "model_index"],
                              na_position="first")
    df.to_csv(a.out, sep="\t", index=False)

    print("\n" + "=" * 78)
    print("S8(a) DEPOSITED")
    print("=" * 78)
    print(f"  {a.out}  ({len(df)} rows)")
    print(f"\n  {len(perm)} permutations, every one zero. Seeds recorded, so each shuffle")
    print(f"  is individually reproducible -- better provenance than the original claim.")
    print(f"\n  {'cohort':<9}{'contrast':<16}{'ET':<5}{'model':<7}{'seed':>6}"
          f"{'tested':>9}{'sig':>5}")
    print("  " + "-" * 62)
    for _, r in df[df.run_kind == "permutation"].iterrows():
        print(f"  {r.cohort:<9}{str(r.contrast)[:15]:<16}{r.event_type:<5}"
              f"{str(r.model_index):<7}{str(r.seed):>6}{r.n_tested:>9}"
              f"{r.n_significant:>5}")

    rr = df[df.run_kind == "real"].sort_values("pct_significant", ascending=False)
    if len(rr):
        print(f"\n  REAL MODELS -- the positive control, by significant fraction:")
        print(f"  {'contrast':<16}{'ET':<5}{'model':<7}{'tested':>8}{'sig':>7}"
              f"{'%':>8}  excluded")
        print("  " + "-" * 60)
        for _, r in rr.iterrows():
            print(f"  {str(r.contrast)[:15]:<16}{r.event_type:<5}"
                  f"{str(r.model_index):<7}{r.n_tested:>8}{r.n_significant:>7}"
                  f"{r.pct_significant:>8.1f}  {r.excluded_from_manuscript}")
        hot = rr[(rr.pct_significant > 10) & (rr.excluded_from_manuscript == "no")]
        if len(hot):
            print(f"\n  [WARN   ] {len(hot)} NON-excluded real model(s) call >10% of")
            print(f"            tested events significant. That is implausible and is NOT")
            print(f"            explained by the permutation nulls, which are zero.")
            print(f"            Investigate before depositing.")
        mx = rr[rr.event_type == "MXE"]
        se = rr[rr.event_type == "SE"]
        if len(mx) and len(se) and se.pct_significant.mean() > 0:
            print(f"\n  >>> MXE calls {mx.pct_significant.mean():.1f}% of tested events "
                  f"significant against SE's {se.pct_significant.mean():.1f}% "
                  f"({mx.pct_significant.mean()/se.pct_significant.mean():.1f}x).")
            print(f"      This is NOT a new argument -- it is how the exclusion was")
            print(f"      first surfaced in Sep 2026 (RESULT_colorado_mxe_excluded.md,")
            print(f"      'How it surfaced': 11.87% and 21.38% against MAGNet's 0.37%).")
            print(f"      The 0.061 reverse slope is what SETTLED it. Reproduced here as")
            print(f"      confirmation, not as fresh evidence.")

    print("\n  ⚠ TWO THINGS TO CARRY WHEREVER THIS IS QUOTED")
    p0 = df[(df.run_kind == "permutation") & (df.model_index == 0) &
            (df.event_type == "SE")]
    p2 = df[(df.run_kind == "permutation") & (df.model_index == 2) &
            (df.event_type == "SE")]
    if len(p0) and len(p2):
        print(f"    1. The covariate-adjusted model tests far fewer events: SE P0 "
              f"{p0.n_tested.min():,}-{p0.n_tested.max():,} against P2 "
              f"{p2.n_tested.min():,}-{p2.n_tested.max():,}")
        print(f"       ({p0.n_tested.mean()/p2.n_tested.mean():.1f}x). A zero over "
              f"{p0.n_tested.min():,} events is strong; a zero over "
              f"{p2.n_tested.min():,} is weaker.")
        print(f"       Quote P0's counts, or quote both and say which is which.")
    print(f"    2. Counts are FULL-COORDINATE. Script 25's real row reads 478 of 25,135")
    print(f"       where the manuscript's exon-level figure is 476 of 19,953. Both are")
    print(f"       correct for what they are. The TSV records the convention per row.")

    if a.caption:
        print("\n" + "=" * 78)
        print("DRAFT S8(a) CAPTION")
        print("=" * 78)
        print(f"""
Supplementary Table S8(a). Label-permutation nulls. Disease labels were shuffled in the
covariate file (random.Random(seed).shuffle) and the full differential-splicing pipeline
re-run. Five independent shuffles in the Colorado cohort -- three of dilated versus
non-failing (seeds 101-103) and two of ischaemic versus non-failing (seeds 101-102) --
were each evaluated under two models and two event classes, giving twenty analyses; two
shuffles in the discovery cohort (seeds 101-102) were evaluated for skipped exons across
103,054 events. **Every one of the twenty-two analyses returned exactly zero significant
events** at FDR < 0.05 and |dPSI| >= 0.10. Tested-event counts differ between shuffles
because permutation changes which events pass the coverage and PSI-variance filters, and
are reported per analysis. Counts are on the full coordinate tuple rather than the
exon-level convention used elsewhere; a count of zero is unaffected by the choice, since
collapsing can only reduce it. The non-permuted models from the same run are included for
comparison as the positive control; these are this script's own full-coordinate counts
under its internal model series and are NOT the event counts reported in the main text,
which are exon-level (Table 1, Supplementary Table S5). Rows from Colorado mutually
exclusive exons are retained for completeness but flagged `excluded_from_manuscript`:
that class is excluded from every table and figure in the paper on reverse-replication
grounds, and its real models call 11-21% of tested events significant against 2% for
skipped exons in the same cohort. Job {JOB}, {SCRIPT}.
""".strip())
    print()


if __name__ == "__main__":
    main()
