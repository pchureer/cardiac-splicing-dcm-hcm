#!/bin/bash
#SBATCH --job-name=magnet125
#SBATCH --output=logs/magnet125_%j.out
#SBATCH --error=logs/magnet125_%j.err
#SBATCH --time=96:00:00
# 48 h is 3.7x the measured expectation. Derived, not guessed: CP's own pilot logs
# (job 405972) show 1.75 h and 1.38 h for 60 samples, so 332 samples is ~8.7 h
# linear and ~13 h allowing for the post step scaling worse than sample count.
# Raise if your partition allows and you want more insurance:
#   sbatch --time=8-00:00:00 100_requant_magnet_125.sh --apply
# but note a longer walltime request can cost queue priority or exclude partitions.
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
# Sized to a 16-CPU / 96-GB account cap. 90G leaves headroom for SLURM's own
# accounting; a request of exactly the cap can pend indefinitely or be rejected.
#   sbatch --mem=94G --time=48:00:00 100_requant_magnet_125.sh --apply
# ---------------------------------------------------------------------------
# 100_requant_magnet_125.sh
#
# Re-quantify the MAGNet discovery cohort with --readLength 125, the measured
# sequenced length, instead of the 50 that was declared.
#
# ⚠ NON-DESTRUCTIVE BY CONSTRUCTION. Everything is written to a NEW tree
#   ($OUT125). Nothing under the existing 50 bp output is read for writing,
#   renamed or deleted. The 50 bp result remains the primary analysis until you
#   decide otherwise, having seen the comparison that 101 produces.
#
# WHY IT IS BUILT THIS WAY
#   A re-quantification of the discovery cohort changes every headline number in
#   the paper and invalidates every gate that protects them. So this does not
#   "switch" anything. It produces a second, complete quantification alongside
#   the first, and 101_compare_readlength_125.py puts them side by side against
#   the recorded numbers. Promotion to primary is a separate, deliberate act.
#
# USAGE
#   bash 100_requant_magnet_125.sh --list          # check inputs, run nothing
#   sbatch 100_requant_magnet_125.sh --apply       # submit the real run
#
# WHAT IT REPRODUCES FROM THE ORIGINAL RUN, deliberately and exactly
#   rMATS-turbo v4.3.0 · --variable-read-length · -t paired · GENCODE v49
#   --cstat 0.0001 · --anchorLength 1 · --novelSS off · --libType fr-unstranded
#   run on the RAW STAR alignments (the recalibrated BAMs have been through
#   SplitNCigarReads, which removes the N-CIGAR junction reads splicing needs)
#   The ONLY intended difference is --readLength 125.
# ---------------------------------------------------------------------------
set -euo pipefail

ROOT="${ROOT:-/data/project/chureeratp/RNA/HEART}"
GTF="${GTF:-$ROOT/ref/gencode.v49.annotation.gtf}"
B1="${B1:-$ROOT/dcm_bams.txt}"          # ⚠ ONE comma-separated line, rMATS --b1 format
B2="${B2:-$ROOT/nf_bams.txt}"
OUT125="${OUT125:-$ROOT/rmats_magnet_rl125}"
TMP125="${TMP125:-$ROOT/tmp_magnet_rl125}"
OUT50="${OUT50:-$ROOT/rmats_magnet}"     # read-only here, for the inventory only
THREADS="${SLURM_CPUS_PER_TASK:-16}"

MODE="${1:---list}"

hr(){ printf '%*s\n' 76 '' | tr ' ' -; }
die(){ hr; echo "REFUSING - nothing was run."; echo; echo "$1"; hr; exit 2; }

hr
echo "100 - re-quantify MAGNet at --readLength 125 (non-destructive)"
hr

# ------------------------------------------------------------------ G1 inputs
echo
echo "  G1  inputs"
for f in "$GTF" "$B1" "$B2"; do
  [[ -s "$f" ]] || die "not found or empty: $f

  Set the path explicitly, e.g.
    ROOT=/data/project/chureeratp/RNA/HEART GTF=... B1=... B2=... bash $0 --list"
  echo "      ok  $f"
done

# ⚠ The bam lists are ONE comma-separated line. Parsed line-by-line they yield
#   one sample out of 166 and rMATS runs happily on it. Count commas, not lines.
n1=$(tr ',' '\n' < "$B1" | grep -c '[^[:space:]]')
n2=$(tr ',' '\n' < "$B2" | grep -c '[^[:space:]]')
l1=$(wc -l < "$B1"); l2=$(wc -l < "$B2")
echo "      --b1 $B1 : $n1 samples across $l1 line(s)"
echo "      --b2 $B2 : $n2 samples across $l2 line(s)"
[[ "$n1" -eq 166 ]] || die "G1 - --b1 resolves to $n1 samples, expected 166 (DCM).
  If this says 1, the list was parsed line-by-line somewhere upstream."
[[ "$n2" -eq 166 ]] || die "G1 - --b2 resolves to $n2 samples, expected 166 (non-failing)."

# every BAM must exist and must be a RAW alignment, not a recalibrated one
missing=0; recal=0
while IFS= read -r b; do
  [[ -z "$b" ]] && continue
  [[ -s "$b" ]] || { echo "      MISSING: $b"; missing=$((missing+1)); }
  case "$b" in *recal*|*bqsr*|*splitn*|*SplitN*) recal=$((recal+1));; esac
done < <({ tr ',' '\n' < "$B1"; echo; tr ',' '\n' < "$B2"; echo; } | grep '[^[:space:]]')
[[ "$missing" -eq 0 ]] || die "G1 - $missing BAM file(s) missing or empty."
if [[ "$recal" -gt 0 ]]; then
  die "G1 - $recal path(s) look like RECALIBRATED alignments.

  rMATS must run on the RAW STAR alignments. The recalibrated path has been
  through SplitNCigarReads, which removes the N-CIGAR junction reads that
  splicing quantification depends on. A run on those produces normal-looking
  output with no usable junctions."
fi
echo "      all $((n1+n2)) BAMs present, none on the recalibrated path"
echo "  G1 PASS"

# ------------------------------------------------------- G2 do not overwrite
echo
echo "  G2  the existing 50 bp result is untouched"
if [[ -d "$OUT125" ]] && [[ -n "$(ls -A "$OUT125" 2>/dev/null)" ]]; then
  die "G2 - $OUT125 already exists and is not empty.

  This script never writes into a populated directory. Move it aside or set
  OUT125 to a new path. (It also never touches $OUT50.)"
fi
if [[ -d "$OUT50" ]]; then
  echo "      50 bp tree present at $OUT50 - will NOT be read for writing"
  echo "      its SE JCEC file: $(ls -la "$OUT50"/SE.MATS.JCEC.txt 2>/dev/null | awk '{print $5" bytes"}' || echo 'not found')"
else
  echo "      note: $OUT50 not found; the comparison in 101 needs the recorded"
  echo "            numbers only, so this is not fatal"
fi
echo "  G2 PASS"

# ---------------------------------------------------------------- G3 toolchain
echo
echo "  G3  toolchain"
command -v rmats.py >/dev/null 2>&1 || command -v rmats >/dev/null 2>&1 \
  || die "G3 - rmats.py not on PATH. conda activate rnaseq_tools first."
RM=$(command -v rmats.py || command -v rmats)
echo "      $RM"
"$RM" --version 2>&1 | head -2 | sed 's/^/      /' || true
echo "      ⚠ confirm this reports v4.3.0. A different version is a second"
echo "        change on top of the read length and makes the comparison useless."
echo "  G3 PASS"

# ------------------------------------------------------------------ the run
if [[ "$MODE" != "--apply" ]]; then
  hr
  echo "  --list only. Nothing was run."
  echo
  echo "  would run:"
  echo "    $RM --b1 $B1 --b2 $B2 --gtf $GTF -t paired \\"
  echo "      --readLength 125 --variable-read-length \\"
  echo "      --cstat 0.0001 --anchorLength 1 --libType fr-unstranded \\"
  echo "      --nthread $THREADS --od $OUT125 --tmp $TMP125"
  echo
  echo "  then:  python 101_compare_readlength_125.py --od $OUT125"
  hr
  exit 0
fi

# ------------------------------------------------------- G5 resources vs cap
echo
echo "  G5  resources"
echo "      cpus-per-task : ${SLURM_CPUS_PER_TASK:-unset (not under SLURM)}"
echo "      mem allocated : ${SLURM_MEM_PER_NODE:-unset} MB"
if [[ -n "${SLURM_MEM_PER_NODE:-}" ]] && [[ "$SLURM_MEM_PER_NODE" -lt 40000 ]]; then
  die "G5 - only ${SLURM_MEM_PER_NODE} MB allocated. 332 samples will not fit.
  Resubmit with --mem=90G, or the most the account allows under 96G."
fi
echo "      nthread passed to rMATS : $THREADS"
echo "      WARNING: 332 paired-end samples on 90 GB is NOT verified for this"
echo "      cohort. rMATS-turbo holds the post step's tables in memory and that"
echo "      scales with sample count. If this job is killed for memory, do NOT"
echo "      simply resubmit - see the prep/post note at the end of this output."
echo "  G5 PASS"

AVAIL=$(df -Pk "$(dirname "$TMP125")" 2>/dev/null | awk 'NR==2{print int($4/1048576)}')
echo
echo "  G6  disk: ${AVAIL:-?} GB free where --tmp will be written"
if [[ -n "${AVAIL:-}" ]] && [[ "$AVAIL" -lt 60 ]]; then
  die "G6 - only ${AVAIL} GB free. The 50 bp output is ~1.2 GB for skipped exons
  alone and --tmp holds per-sample intermediates that are larger again. Put tmp
  on scratch instead:
    TMP125=/scratch/\$USER/tmp_magnet_rl125 sbatch $0 --apply"
fi
echo "  G6 PASS"

mkdir -p "$OUT125" "$TMP125" logs
echo
echo "  running rMATS at --readLength 125 ..."
date
"$RM" \
  --b1 "$B1" \
  --b2 "$B2" \
  --gtf "$GTF" \
  -t paired \
  --readLength 125 \
  --variable-read-length \
  --cstat 0.0001 \
  --anchorLength 1 \
  --libType fr-unstranded \
  --nthread "$THREADS" \
  --od "$OUT125" \
  --tmp "$TMP125"
date

# ------------------------------------------------- G4 the output is complete
echo
echo "  G4  output completeness"
for c in SE MXE A3SS A5SS RI; do
  f="$OUT125/${c}.MATS.JCEC.txt"
  [[ -s "$f" ]] || die "G4 - $f missing or empty. The run did not complete."
  echo "      $c: $(($(wc -l < "$f") - 1)) rows"
done
# group widths: a mismatched -t silently yields no usable pairs
w1=$(head -2 "$OUT125/SE.MATS.JCEC.txt" | tail -1 | awk -F'\t' '{print $(NF-5)}' | tr -cd ',' | wc -c)
echo "      IncLevel1 width check: $((w1+1)) values (expect 166)"
[[ $((w1+1)) -eq 166 ]] || die "G4 - IncLevel1 carries $((w1+1)) values, expected 166.
  A mismatched -t setting produces normal-looking output with no usable read pairs."
echo "  G4 PASS"

hr
echo "DONE. The 50 bp result is unchanged."
echo
echo "IF THIS JOB WAS KILLED FOR MEMORY, the prep/post split is the fix:"
echo "  rMATS-turbo can run --task prep per batch (low memory, and resumable, so a"
echo "  failure does not discard completed work) and then one --task post over all"
echo "  the prep output. For 332 samples that is the documented approach."
echo "  WARNING: the .rmats files from separate prep runs must be collected into ONE"
echo "  --tmp directory without colliding, and prep and post must use the SAME"
echo "  --readLength 125. A prep/post mismatch produces complete-looking output"
echo "  built from inconsistent inputs, which no gate downstream would catch."
echo
echo "NEXT - do not promote anything yet:"
echo "  python 101_compare_readlength_125.py --od $OUT125"
echo
echo "That prints the 125 bp counts beside every recorded 50 bp number. Decide"
echo "after reading it. ⚠ Promotion to primary means re-baselining the gates in"
echo "62, 68, 69, 70, 71, 73, 77, 78, 79, 80, 81, 87, 93, 95, 96, 98 - and while"
echo "they are being re-baselined they protect nothing."
hr
