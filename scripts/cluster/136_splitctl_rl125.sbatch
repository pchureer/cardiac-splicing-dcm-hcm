#!/bin/bash
#SBATCH --job-name=splitctl_rl125
#SBATCH -p cpu
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH --qos=ngs16
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#
# 136_splitctl_rl125.sbatch -- the 125 bp split-control arms that script 22's
# MAGNet-internal row needs.
#
# HEADER: copied from 16_split_control_check.sbatch, which is the job this one
# repeats against a different MATS file. `--qos=ngs16` is REQUIRED on this
# cluster's cpu partition (omitting it fails with "Invalid qos specification"),
# and 16 sets no --time, so none is set here either.
#
# ⚠ --mem=32G is 16's PROVEN value, not a measured one. An earlier draft of
#   this script asked for 8G, citing MaxRSS 3.45 G — but that measurement came
#   from script 56 and the RBM20 target tests, NOT from
#   phase2_diffsplice_v2.py, which reads 385,826 rows x 332 samples and builds
#   a PSI matrix from them (~1 GB as float64 before pandas overhead and before
#   the raw count strings). Carrying a memory figure across unrelated scripts
#   is the same mistake as carrying a slope across read lengths.
#   AFTER THIS RUNS:  sacct -j <jobid> --format=JobID,MaxRSS,Elapsed,State
#   and size the next run on that. If it queues on QOSMaxMemoryPerUser, drop
#   to 16G and re-submit rather than waiting.
#
# WHY
# ---
# 22_symmetry_winners_curse.py's first pair is
#     x  splitctl/{ET}_DCMvsA_1.all_events.tsv     from rmats_magnet/
#     y  splitctl/{ET}_HCMvsB_1.all_events.tsv     from rmats_HCMvsNF/
# and gives the 0.92 "MAGNet compared with itself" slope that 22's headline
# contrasts against the ~0.56 cohort slopes. The cohort rows re-base to 125 bp
# by suffixing their discovery stem; this pair cannot, because BOTH its sides
# are MAGNet. Leaving it at 50 bp while the cohort rows move would make 22's
# central contrast a mixed-read-length comparison — the same error that
# produced the false 28% gap in the matched-pair closure, where a 125 bp A x B
# was measured against a hard-coded 50 bp C.
#
# NO NEW ALIGNMENT AND NO NEW rMATS RUN.
# 16's splits are group re-assignments over an existing rMATS output, driven by
# `--exclude <ids>`; 19's header says so explicitly ("No new alignment").
# Everything here changes is which MATS file phase 2 reads.
#
# THE SPLIT IS REUSED VERBATIM, NOT RE-DRAWN.
# 16 wrote the partition to disk (`splitctl/nfA_1.ids`, `nfB_1.ids`, 83 + 83 of
# 166 NF hearts) rather than only seeding an RNG, so the 125 bp arms use the
# SAME hearts in the same halves. Nothing here re-shuffles anything; if either
# .ids file is missing the script refuses rather than regenerating it, because
# a regenerated split would differ from the published one.
#
# WHAT IT DOES, IN ORDER
#   0  verifies every input, and that the 125 bp MATS files carry the same
#      per-group sample arity as the 50 bp ones (a mismatch would mean phase 2
#      silently mis-assigns samples to groups, which a calibration run on the
#      50 bp MATS cannot detect)
#   1  CALIBRATION: re-runs the 50 bp DCMvsA_1 and HCMvsB_1 to .CALCHECK
#      prefixes and `cmp`s them against the existing tables. If the command
#      line here is not byte-identical in effect to 16's, we find out BEFORE
#      any 125 bp output is believed. Refuses on mismatch.
#   2  the 125 bp arms, written to *_1_rl125 — never over the 50 bp files
#   3  prints the two logs' group sizes and row counts side by side
#
# SCOPE: SE, split 1 only (what 22 reads). Test one before scaling.
#   ET=MXE sbatch 136_splitctl_rl125.sbatch     for another class
#   S=2    sbatch 136_splitctl_rl125.sbatch     for another split
# The five-split and SHARED/DISJOINT arms that 17 and 83 use are NOT built
# here. 83's disjoint range (0.87-0.99, threshold-dependent by its own note) is
# a separate claim and a separate decision.
#
set -uo pipefail

ET=${ET:-SE}
S=${S:-1}
SUF=${SUF:-_rl125}
FORCE=${FORCE:-0}

# live directories, pinned. NOTE: rmats_magnet_rl125.abandoned_1207 and
# rmats_magnet_rl125_tmp.abandoned_1207 also exist; they are named here only to
# make clear they are deliberately NOT used.
MD50=rmats_magnet/${ET}.MATS.JCEC.txt
MH50=rmats_HCMvsNF/${ET}.MATS.JCEC.txt
MD125=rmats_magnet_rl125/${ET}.MATS.JCEC.txt
MH125=rmats_HCMvsNF_rl125/${ET}.MATS.JCEC.txt

A=splitctl/nfA_${S}.ids
B=splitctl/nfB_${S}.ids

OUT_D50=splitctl/${ET}_DCMvsA_${S}
OUT_H50=splitctl/${ET}_HCMvsB_${S}
OUT_D125=splitctl/${ET}_DCMvsA_${S}${SUF}
OUT_H125=splitctl/${ET}_HCMvsB_${S}${SUF}

hr() { printf '%s\n' "----------------------------------------------------------------------"; }
die() { echo; echo "REFUSING — $*"; echo "Nothing was written."; exit 1; }

echo "=============================================================="
echo "136 — split-control arms at 125 bp.  ET=$ET  split=$S  suffix=$SUF"
echo "=============================================================="
date
echo "host $(hostname)   job ${SLURM_JOB_ID:-none}"
echo

# ---------------------------------------------------------------- stage 0
hr; echo "STAGE 0 — inputs"; hr
for f in "$MD50" "$MH50" "$MD125" "$MH125" "$A" "$B" \
         dcm_bams.txt nf_bams.txt hcm_bams.txt \
         magnet_covariates.tsv hcm_nf_covariates.tsv \
         phase2_diffsplice_v2.py; do
  [ -s "$f" ] || die "missing or empty: $f"
  printf '  ok  %-56s %s\n' "$f" "$(du -h "$f" | cut -f1)"
done

echo
echo "  the split, reused verbatim (NOT re-drawn):"
printf '    %-28s %s hearts\n' "$A" "$(grep -c . "$A")"
printf '    %-28s %s hearts\n' "$B" "$(grep -c . "$B")"
if [ "$(sort "$A" "$B" | sort -u | wc -l)" -ne \
     "$(( $(grep -c . "$A") + $(grep -c . "$B") ))" ]; then
  die "the two halves overlap — $A and $B share at least one sample"
fi
echo "    disjoint: yes"

# arity: a 125 bp MATS built from a different bam list or order would make
# phase 2 mis-assign samples, and the stage-1 calibration CANNOT see it
# because calibration reads the 50 bp MATS.
echo
echo "  per-group sample arity in the raw input (must match across read lengths):"
arity() {
  local d="$1" raw
  raw="${d}/JCEC.raw.input.${ET}.txt"
  [ -s "$raw" ] || { echo "MISSING"; return; }
  awk 'NR==2{n1=gsub(/,/,",",$2)+1; n2=gsub(/,/,",",$4)+1; print n1"/"n2; exit}' "$raw"
}
AR_MD50=$(arity rmats_magnet)
AR_MD125=$(arity rmats_magnet_rl125)
AR_MH50=$(arity rmats_HCMvsNF)
AR_MH125=$(arity rmats_HCMvsNF_rl125)
printf '    %-26s %-12s %-26s %s\n' rmats_magnet "$AR_MD50" rmats_magnet_rl125 "$AR_MD125"
printf '    %-26s %-12s %-26s %s\n' rmats_HCMvsNF "$AR_MH50" rmats_HCMvsNF_rl125 "$AR_MH125"
[ "$AR_MD50" = "$AR_MD125" ] || die "magnet arity differs: 50bp $AR_MD50 vs 125bp $AR_MD125.
  The 125 bp quantification is not column-compatible with the bam lists, so
  phase 2 would assign samples to the wrong groups. This is bigger than
  script 22 — every phase-2 script shares this machinery."
[ "$AR_MH50" = "$AR_MH125" ] || die "HCMvsNF arity differs: 50bp $AR_MH50 vs 125bp $AR_MH125.
  Same conclusion as above."
echo "    arity matches both ways."

# don't clobber
for o in "$OUT_D125" "$OUT_H125"; do
  if [ -s "${o}.all_events.tsv" ] && [ "$FORCE" != "1" ]; then
    die "${o}.all_events.tsv already exists.
  Re-run with FORCE=1 to overwrite, after checking nothing downstream has
  already consumed it."
  fi
done
echo "  outputs are free to write."

# ---------------------------------------------------------------- stage 1
echo
hr; echo "STAGE 1 — CALIBRATION on the 50 bp MATS"; hr
echo "  If this command line is not equivalent to 16's, the 125 bp numbers"
echo "  would be wrong in a way no later check would catch. So reproduce the"
echo "  EXISTING 50 bp tables first, byte for byte."
echo

cal_fail=0
for arm in DCM HCM; do
  if [ "$arm" = "DCM" ]; then
    mats="$MD50"; b1=dcm_bams.txt; cov=magnet_covariates.tsv
    g1=DCM; excl="$B"; ref="$OUT_D50"; tag=DCMvsA
  else
    mats="$MH50"; b1=hcm_bams.txt; cov=hcm_nf_covariates.tsv
    g1=HCM; excl="$A"; ref="$OUT_H50"; tag=HCMvsB
  fi
  if [ ! -s "${ref}.all_events.tsv" ]; then
    echo "  [skip] ${ref}.all_events.tsv absent — nothing to calibrate against"
    continue
  fi
  echo "  ${tag}: re-running 50 bp -> ${ref}.CALCHECK"
  python phase2_diffsplice_v2.py --mats "$mats" --event-type "$ET" \
    --b1 "$b1" --b2 nf_bams.txt --cov "$cov" \
    --group1 "$g1" --group2 NF --min-samples 15 --exclude "$excl" \
    --out "${ref}.CALCHECK" > "${ref}.CALCHECK.log" 2>&1
  if [ ! -s "${ref}.CALCHECK.all_events.tsv" ]; then
    echo "  ✗ ${tag}: produced no table — see ${ref}.CALCHECK.log"
    tail -5 "${ref}.CALCHECK.log"
    cal_fail=1
    continue
  fi
  if cmp -s "${ref}.CALCHECK.all_events.tsv" "${ref}.all_events.tsv"; then
    echo "  ✓ ${tag}: byte-identical to the existing 50 bp table"
    rm -f "${ref}.CALCHECK.all_events.tsv"
  else
    echo "  ✗ ${tag}: DIFFERS from ${ref}.all_events.tsv"
    echo "      rows: $(wc -l < "${ref}.CALCHECK.all_events.tsv") new vs $(wc -l < "${ref}.all_events.tsv") existing"
    diff <(head -3 "${ref}.all_events.tsv") \
         <(head -3 "${ref}.CALCHECK.all_events.tsv") | head -10
    cal_fail=1
  fi
done
[ "$cal_fail" = "0" ] || die "calibration failed.
  The 50 bp re-run does not reproduce 16's output, so this command line is not
  the one that produced the published tables. Do NOT trust a 125 bp run from
  it. Send me the .CALCHECK.log and the diff above.
  (The .CALCHECK tables were left in place for inspection.)"

# ---------------------------------------------------------------- stage 2
echo
hr; echo "STAGE 2 — the 125 bp arms"; hr
echo "  DCMvsA: $MD125   (controls = half A, excluding $B)"
python phase2_diffsplice_v2.py --mats "$MD125" --event-type "$ET" \
  --b1 dcm_bams.txt --b2 nf_bams.txt --cov magnet_covariates.tsv \
  --group1 DCM --group2 NF --min-samples 15 --exclude "$B" \
  --out "$OUT_D125" > "${OUT_D125}.log" 2>&1
[ -s "${OUT_D125}.all_events.tsv" ] || { tail -20 "${OUT_D125}.log"; \
  die "DCMvsA at 125 bp produced no table — see ${OUT_D125}.log"; }
echo "  [write] ${OUT_D125}.all_events.tsv"

echo "  HCMvsB: $MH125   (controls = half B, excluding $A)"
python phase2_diffsplice_v2.py --mats "$MH125" --event-type "$ET" \
  --b1 hcm_bams.txt --b2 nf_bams.txt --cov hcm_nf_covariates.tsv \
  --group1 HCM --group2 NF --min-samples 15 --exclude "$A" \
  --out "$OUT_H125" > "${OUT_H125}.log" 2>&1
[ -s "${OUT_H125}.all_events.tsv" ] || { tail -20 "${OUT_H125}.log"; \
  die "HCMvsB at 125 bp produced no table — see ${OUT_H125}.log"; }
echo "  [write] ${OUT_H125}.all_events.tsv"

# ---------------------------------------------------------------- stage 3
echo
hr; echo "STAGE 3 — 50 bp against 125 bp, from the logs"; hr
for pair in "$OUT_D50:$OUT_D125" "$OUT_H50:$OUT_H125"; do
  o50="${pair%%:*}"; o125="${pair##*:}"
  echo
  echo "  $(basename "$o125")"
  for o in "$o50" "$o125"; do
    [ -s "${o}.log" ] || { printf '    %-46s (no log)\n' "$(basename "$o")"; continue; }
    printf '    %-46s %s\n' "$(basename "$o")" \
      "$(grep -m1 -E 'event type|samples' "${o}.log" | sed 's/^\[info\] *//')"
    printf '    %-46s %s\n' "" \
      "$(grep -m1 'excluding' "${o}.log" | sed 's/^\[info\] *//')"
  done
  [ -s "${o50}.all_events.tsv" ] && [ -s "${o125}.all_events.tsv" ] && \
    printf '    rows: %s (50 bp)   %s (125 bp)\n' \
      "$(($(wc -l < "${o50}.all_events.tsv") - 1))" \
      "$(($(wc -l < "${o125}.all_events.tsv") - 1))"
done

echo
echo "=============================================================="
echo "DONE. Both 125 bp arms exist and the 50 bp tables are untouched."
echo "=============================================================="
cat <<'NOTE'

⚠ READ THE GROUP SIZES ABOVE BEFORE GOING FURTHER. The partition is the same
  83/83 of 166 NF hearts, but if the 125 bp quantification kept a different
  number of samples the effective control group differs and the two read
  lengths are then not the same comparison. "excluding 83 samples (83 of them
  present in the bam lists)" is what the 50 bp run reported; the 125 bp run
  must report the same.

NEXT, and NOT done here:
  * 22_symmetry_winners_curse.py needs patching before it can read these.
    Its PAIRS list (lines 135-143) hard-codes:
        splitctl/{ET}_DCMvsA_1.all_events.tsv      <- this pair moves
        splitctl/{ET}_HCMvsB_1.all_events.tsv      <- with BOTH sides
        {ET}_v2filtered.all_events.tsv             <- discovery side only
        {ET}_Colorado_DCMvsNF.all_events.tsv  etc. <- MUST stay unsuffixed
    and `--out` defaults to f"symmetry_{ET}" with no read-length component,
    so an unsuffixed 125 bp run would overwrite symmetry_SE.tsv,
    symmetry_SE_sweep.tsv and symmetry_SE.report.txt — the same defect 33 had.
    Do not run 22 at 125 bp until it is patched.
  * 22's docstring line 8 carries all five 50 bp slopes (Colorado 0.571,
    LBNL 0.525, Brigham 0.757, Fuwai 0.359, Stanford 0.578) and the "~0.56"
    and "0.92" headline figures. They are read by 82 (GENOME_WIDE), 62 (S5)
    and 98 (line 393, inside a refusal message).
  * 22's lines 9-10 say read length was "excluded (scripts 16/17/19/20/21)".
    That exclusion was about COHORT read length, not MAGNet's own
    quantification read length, which is the axis this whole re-baselining
    moves. The sentence needs to name which hypothesis was tested.
NOTE
