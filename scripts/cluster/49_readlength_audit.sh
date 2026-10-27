#!/bin/bash
# 49_readlength_audit.sh -- is MAGNet 50 bp or 125 bp, and what did rMATS think?
#
# WHY THIS EXISTS
# ---------------
# 48_provenance_versions.sh sampled ONE bam per cohort and found MAGNet reads at 125 bp.
# Table 1, Figure 1, PAPER_ONE_PAGE.md and the Methods all say "paired, 50 bp", and every
# MAGNet rMATS invocation passed `--readLength 50`. One bam is not a cohort, so this
# script checks EVERY sample in every list before anything in the manuscript is changed.
#
# WHY IT MATTERS BEYOND A TABLE CELL
# ----------------------------------
# rMATS normalises PSI by effective form lengths that are computed from --readLength:
# inclusion and skipping forms have different lengths, and the declared read length sets
# both. A declared 50 against true 125 therefore changes PSI itself, not just a label.
#
#   * WITHIN a cohort both arms carry the same constant, so dPSI direction is safe and
#     the discovery result is very unlikely to move qualitatively.
#   * ACROSS cohorts nothing cancels, and the paper's attenuation story is literally
#     "slopes order with read-length distance from MAGNet's paired 50 bp". If MAGNet is
#     125 bp, the slopes order with distance from the DECLARED 50, not the true length —
#     a different claim that must be stated differently.
#
# READ-ONLY. Reads a capped number of records per bam; writes only its own report.
#
# Usage:
#   cd /data/project/chureeratp/RNA/HEART
#   bash 49_readlength_audit.sh 2>&1 | tee readlength_audit.txt
#   bash 49_readlength_audit.sh --max-samples 0      # every sample, slow
#   sbatch 49_readlength_audit.sbatch                # the right way on a cluster
#   bash 49_readlength_audit.sh --nreads 50000

set -uo pipefail

ROOT="/data/project/chureeratp/RNA/HEART"
ENVBIN="/data/apps/anaconda3/envs/mamba/envs/rnaseq_tools/bin"
NREADS=5000           # records sampled per bam (plenty for a modal length)
MAXS=40               # bams per list; 0 = all. 40 settles a cohort mode

while [[ $# -gt 0 ]]; do
  case "$1" in
    --root)        ROOT="$2"; shift 2 ;;
    --env)         ENVBIN="$2"; shift 2 ;;
    --nreads)      NREADS="$2"; shift 2 ;;
    --max-samples) MAXS="$2"; shift 2 ;;
    -h|--help)     sed -n '2,30p' "$0"; exit 0 ;;
    *) echo "unknown option: $1"; exit 2 ;;
  esac
done

export PATH="$ENVBIN:$PATH"
cd "$ROOT" 2>/dev/null || { echo "cannot cd to $ROOT"; exit 1; }
command -v samtools >/dev/null 2>&1 || { echo "samtools not on PATH"; exit 1; }

rule () { printf '\n%s\n%s\n%s\n' "$(printf '=%.0s' {1..78})" "$1" "$(printf '=%.0s' {1..78})"; }

echo "read-length audit"
echo "date   : $(date -Is)"
echo "root   : $ROOT"
echo "sample : $NREADS records per bam$([[ $MAXS -gt 0 ]] && echo ", first $MAXS bams per list")"

# ---------------------------------------------------------------------------
rule "1. DECLARED --readLength, straight out of the job scripts"
# ---------------------------------------------------------------------------
# The job scripts live in THIS directory. An earlier version used `grep -r .`, which
# recursed into 716 sample directories and the rMATS *_tmp trees -- hundreds of
# thousands of small files, on a shared filesystem, while another job was writing more
# of them. It burned a two-hour wall limit without reading a single bam. Depth 1, by
# name, is the whole job.
mapfile -t SCRIPTS < <(find . -maxdepth 1 \( -name '*.sbatch' -o -name '*.sh' \) \
                            ! -name '49_readlength_audit.sh' 2>/dev/null)
echo "  scanning ${#SCRIPTS[@]} job script(s) in $ROOT (depth 1, no recursion)"

if [[ ${#SCRIPTS[@]} -eq 0 ]]; then
  echo "  *** no .sbatch/.sh found here — pass --root to where the job scripts live"
else
  echo
  echo "  literal --readLength values on uncommented rmats.py lines:"
  for f in "${SCRIPTS[@]}"; do
    sed -e :a -e '/\\$/N; s/\\\n//; ta' "$f" 2>/dev/null \
      | grep -E 'rmats\.py' | grep -vE '^[[:space:]]*#'
  done | grep -oE -- '--readLength[= ]+[^ ]+' | sort | uniq -c | sed 's/^/    /'

  echo
  echo "  RL array definitions (what \${RL[\$COH]} expands to):"
  grep -hnE 'RL\[|declare -A RL|RL=\(' "${SCRIPTS[@]}" 2>/dev/null \
    | grep -vE '^\s*#' | sort -u | head -20 | sed 's/^/    /'
  [[ $(grep -hlE 'RL\[' "${SCRIPTS[@]}" 2>/dev/null | wc -l) -eq 0 ]] \
    && echo "    (no RL array — every --readLength is a literal above)"
fi

# ---------------------------------------------------------------------------
rule "2. OBSERVED read length, every sample in every list"
# ---------------------------------------------------------------------------
printf '  %-30s %6s %10s %10s %8s  %s\n' \
       "list" "bams" "modal bp" "range" "layout" "samples off the mode"

shopt -s nullglob
for lst in dcm_bams.txt nf_bams.txt hcm_bams.txt *_bams_dedup.txt *_bams.txt; do
  [[ -f "$lst" ]] || continue
  mapfile -t BAMS < <(tr ',' '\n' < "$lst" | sed 's/[[:space:]]//g' | grep '\.bam$')
  [[ ${#BAMS[@]} -eq 0 ]] && { printf '  %-30s no bam paths parsed\n' "$lst"; continue; }
  n=${#BAMS[@]}
  lim=$n; [[ $MAXS -gt 0 && $MAXS -lt $n ]] && lim=$MAXS

  modes=(); layout="?"
  for ((i=0; i<lim; i++)); do
    b="${BAMS[$i]}"; [[ -f "$b" ]] || b="$ROOT/${BAMS[$i]}"
    [[ -f "$b" ]] || continue
    # ONE pass per bam for both answers. The first version called
    # `samtools view -c -f 1` to get the layout, which reads the WHOLE bam --
    # minutes per file on a 3-6 GB alignment, and the reason this script crawled.
    # Flag bit 1 is the lowest bit, so flag % 2 == 1 is "paired", and a thousand
    # records settle it as well as sixty million do.
    read -r m lay < <(samtools view "$b" 2>/dev/null | head -n "$NREADS" | awk '
      { c[length($10)]++; if ($2 % 2 == 1) p++ }
      END { best = 0; bl = 0
            for (k in c) if (c[k] > best) { best = c[k]; bl = k }
            print bl, (p > 0 ? "PAIRED" : "SINGLE") }')
    [[ -n "$m" && "$m" != "0" ]] && modes+=("$m")
    [[ "$layout" == "?" && -n "$lay" ]] && layout="$lay"
  done
  [[ ${#modes[@]} -eq 0 ]] && { printf '  %-30s no readable bam\n' "$lst"; continue; }

  overall=$(printf '%s\n' "${modes[@]}" | sort -n | uniq -c | sort -rn | head -1 | awk '{print $2}')
  lo=$(printf '%s\n' "${modes[@]}" | sort -n | head -1)
  hi=$(printf '%s\n' "${modes[@]}" | sort -n | tail -1)
  off=$(printf '%s\n' "${modes[@]}" | grep -vc "^${overall}$")
  printf '  %-30s %6s %10s %10s %8s  %s\n' \
         "$lst" "${#modes[@]}" "$overall" "${lo}-${hi}" "$layout" "$off"
done
shopt -u nullglob

# ---------------------------------------------------------------------------
rule "3. THE COMPARISON THAT MATTERS"
# ---------------------------------------------------------------------------
cat <<'TXT'
  Put section 1 next to section 2 and answer one question per cohort:
  does the DECLARED --readLength equal the OBSERVED modal read length?

  Known from the single-bam pass in script 48:
      Colorado   declared 51    observed 51     ok
      LBNL       declared ?     observed 50/51  ok (mixed, handled)
      Brigham    declared ?     observed 101
      Stanford   declared ?     observed 150
      Fuwai      declared ?     observed 150
      MAGNet     declared 50    observed 125    *** MISMATCH ***

  If section 2 confirms 125 bp across all 332 MAGNet samples, then:

  1. TABLE 1 AND FIGURE 1 ARE WRONG. MAGNet is "paired, 125 bp". Figure 1 is built
     from Table1.tsv, so rebuild it after fixing the table (43_figure1_design.py
     cross-checks and will complain if you do not).

  2. THE ATTENUATION SENTENCE CHANGES. It currently reads that slopes order with
     read-length distance from MAGNet's paired 50 bp. Against the TRUE lengths that
     ordering fails — Colorado is the most distant (51 vs 125) and has the HIGHEST
     slope, 0.521, while Brigham at 101 bp is closest and sits lower at 0.476. Against
     the DECLARED lengths it holds, because every other cohort was declared at its true
     length and only MAGNet was not. The honest statement is that transfer tracks the
     read length rMATS was TOLD to use, which is a statement about the normalisation
     constant rather than about sequencing chemistry. That is a more interesting
     methods result, not a weaker one — but it is a different sentence.

  3. RUN THE DECISIVE TEST BEFORE REWRITING ANYTHING. Re-run MAGNet discovery with
     --readLength 125 and compare against the frozen 926 SE / 1,445 all-types:

       rmats.py --b1 dcm_bams.txt --b2 nf_bams.txt \
         --gtf /data/project/chureeratp/hg38/ref/gencode.v49.annotation.gtf \
         -t paired --readLength 125 --variable-read-length --nthread 16 \
         --od rmats_magnet_rl125 --tmp rmats_magnet_rl125_tmp

     Then compare, at the locked convention (exon-level, FDR<0.05, |dPSI|>=0.10):
       * significant SE count against 926
       * Spearman rho and sign concordance of dPSI on the shared event set
       * whether the DCM-vs-HCM null still returns 0 of ~45,380

     EXPECTATION, stated in advance so the result is a test and not a rationalisation:
     within-cohort dPSI should be nearly unchanged in direction (both arms carry the
     same constant) with counts shifting modestly as events cross the 0.10 threshold.
     If the counts move a lot, every MAGNet number in the paper is provisional and the
     cross-cohort slopes must be recomputed from the corrected run. If they barely
     move, say so explicitly in the Methods and keep the frozen numbers.

  4. DO NOT quietly re-run and swap the numbers. Whatever happens, the Methods should
     state the declared --readLength, the observed read length, and what the check
     showed. A reviewer who downloads PRJNA595151 sees 125 bp immediately.
TXT

echo
echo "done. Save this: bash 49_readlength_audit.sh 2>&1 | tee readlength_audit.txt"
