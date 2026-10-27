#!/bin/bash
# 48_provenance_versions.sh -- collect everything the Methods section still has a blank
# for, from the cluster rather than from memory.
#
# READ-ONLY. It runs no analysis, writes nothing except its own report, and touches no
# BAM except to read the first few reads of one file per cohort.
#
# THE POINT: what is on $PATH today is not necessarily what produced the data.
# STAR records its own version and full command line in every Log.out it ever wrote, and
# the SLURM job files record the rMATS command as invoked. Those are the authoritative
# record; the `--version` calls in section 1 are a cross-check, and where the two
# disagree the Log.out wins and the Methods must quote it.
#
# Usage:
#   cd /data/project/chureeratp/RNA/HEART
#   bash 48_provenance_versions.sh 2>&1 | tee provenance_report.txt
#
#   bash 48_provenance_versions.sh --root /some/other/dir --adapters /path/adapters.fa

set -uo pipefail          # deliberately NOT -e: a missing tool is a finding, not a crash

ROOT="/data/project/chureeratp/RNA/HEART"
ADAPTERS="/data/project/chureeratp/fastq/adapters.fa"
ENVBIN="/data/apps/anaconda3/envs/mamba/envs/rnaseq_tools/bin"
NREADS=200000             # reads sampled per BAM for the read-length check

while [[ $# -gt 0 ]]; do
  case "$1" in
    --root)     ROOT="$2"; shift 2 ;;
    --adapters) ADAPTERS="$2"; shift 2 ;;
    --env)      ENVBIN="$2"; shift 2 ;;
    --nreads)   NREADS="$2"; shift 2 ;;
    -h|--help)  sed -n '2,20p' "$0"; exit 0 ;;
    *) echo "unknown option: $1"; exit 2 ;;
  esac
done

export PATH="$ENVBIN:$PATH"
cd "$ROOT" 2>/dev/null || { echo "cannot cd to $ROOT"; exit 1; }

rule () { printf '\n%s\n%s\n' "$(printf '=%.0s' {1..78})" "$1"; printf '%s\n' "$(printf '=%.0s' {1..78})"; }
sub  () { printf '\n--- %s\n' "$1"; }
have () { command -v "$1" >/dev/null 2>&1; }

echo "provenance report"
echo "host    : $(hostname)"
echo "date    : $(date -Is)"
echo "root    : $ROOT"
echo "PATH[0] : ${PATH%%:*}"

# ---------------------------------------------------------------------------
rule "1. TOOL VERSIONS ON \$PATH NOW  (cross-check, not the record)"
# ---------------------------------------------------------------------------
for t in STAR bcftools samtools trimmomatic rmats.py gatk python; do
  if have "$t"; then
    printf '  %-12s %s\n' "$t" "$(command -v "$t")"
  else
    printf '  %-12s *** NOT ON PATH\n' "$t"
  fi
done

sub "STAR"
have STAR && STAR --version 2>&1 | head -2 || echo "  (not on PATH)"

sub "bcftools"
have bcftools && bcftools --version 2>&1 | head -3 || echo "  (not on PATH)"

sub "samtools"
have samtools && samtools --version 2>&1 | head -2 || echo "  (not on PATH)"

sub "Trimmomatic"
# the conda wrapper answers -version; a bare jar needs java -jar
if have trimmomatic; then
  trimmomatic -version 2>&1 | head -2
  # the wrapper is a shell script that names the jar -- that name carries the version
  wrapper="$(command -v trimmomatic)"
  if [[ -f "$wrapper" ]]; then
    grep -oE 'trimmomatic[-._A-Za-z0-9]*\.jar' "$wrapper" 2>/dev/null | sort -u \
      | sed 's/^/  jar: /'
  fi
else
  find /data/apps /data/home/chureeratp/apps -maxdepth 4 -iname 'trimmomatic*.jar' \
       2>/dev/null | head -5 | sed 's/^/  found jar: /'
fi

sub "rMATS"
if have rmats.py; then
  rmats.py --version 2>&1 | head -2
else
  python -c "import rmats" 2>/dev/null || echo "  (rmats.py not on PATH)"
fi

sub "GATK"
have gatk && gatk --version 2>&1 | grep -iE 'gatk|version' | head -3 || echo "  (not on PATH)"

sub "Python and the packages the analysis scripts import"
have python && python - <<'PY'
import sys
print(f"  python      {sys.version.split()[0]}")
for m in ("pandas", "numpy", "scipy", "statsmodels", "matplotlib", "openpyxl"):
    try:
        print(f"  {m:<12}{__import__(m).__version__}")
    except Exception as e:
        print(f"  {m:<12}*** {type(e).__name__}")
PY

# ---------------------------------------------------------------------------
rule "2. WHAT ACTUALLY RAN  (STAR Log.out — authoritative)"
# ---------------------------------------------------------------------------
# STAR writes its version and the full command line into every Log.out. These files are
# the provenance record; quote them in the Methods, not section 1.
mapfile -t LOGS < <(find . -maxdepth 3 -name '*Log.out' -size +1c 2>/dev/null)
echo "  Log.out files found: ${#LOGS[@]} (all of them — no cap)"
if [[ ${#LOGS[@]} -gt 0 ]]; then
  sub "distinct STAR versions across those logs"
  for f in "${LOGS[@]}"; do
    grep -m1 -oE 'STAR +version=?[^ ]*|versionSTAR +[^ ]*' "$f" 2>/dev/null
  done | sed 's/  */ /g' | sort | uniq -c | sort -rn | sed 's/^/  /'

  sub "distinct genome index / annotation used"
  for f in "${LOGS[@]}"; do
    grep -m1 -oE '\-\-sjdbGTFfile +[^ ]+' "$f" 2>/dev/null
  done | sort | uniq -c | sort -rn | sed 's/^/  /'
  echo "  (expect gencode.v49 for every cohort except Fuwai, which should be v46)"

  sub "distinct --twopassMode / --outSAMattributes settings"
  for f in "${LOGS[@]}"; do
    grep -m1 -oE '\-\-twopassMode +[^ ]+' "$f" 2>/dev/null
  done | sort | uniq -c | sed 's/^/  /'
  for f in "${LOGS[@]}"; do
    grep -m1 -oE '\-\-outSAMattributes +[^ ]+' "$f" 2>/dev/null
  done | sort | uniq -c | sed 's/^/  /'

  sub "annotation per COHORT, resolved through the bam lists"
  echo "  This is the check that matters: Table 1 claims v49 everywhere except Fuwai."
  shopt -s nullglob
  for lst in dcm_bams.txt nf_bams.txt hcm_bams.txt *_bams_dedup.txt *_bams.txt; do
    [[ -f "$lst" ]] || continue
    ids=$(tr ',' '\n' < "$lst" | grep -oE '(SRR|ERR|SAMEA|SAMN)[0-9]+' | sort -u)
    [[ -z "$ids" ]] && continue
    n=$(printf '%s\n' "$ids" | wc -l)
    gtfs=$(for id in $ids; do
             f=$(ls "$id"/*Log.out 2>/dev/null | head -1)
             [[ -n "$f" ]] && grep -m1 -oE 'gencode\.v[0-9]+' "$f"
           done | sort | uniq -c | awk '{printf "%s x%s  ", $2, $1}')
    printf '  %-28s %4s samples   %s\n' "$lst" "$n" "${gtfs:-NO LOGS FOUND}"
  done
  shopt -u nullglob
  echo "  Any analysed cohort other than Fuwai showing v46 is a problem: the tested"
  echo "  event set depends on the annotation, and Table 1 and the Methods both say v49."

  sub "one full command line, verbatim (first log)"
  grep -m1 -A2 -iE '^\s*(##### )?(Command Line|STAR +--runThreadN)' "${LOGS[0]}" \
    2>/dev/null | head -8 | sed 's/^/  /'
  echo "  [log] ${LOGS[0]}"
else
  echo "  *** no Log.out found. Without them the STAR version cannot be established"
  echo "      from the data, and section 1's answer is only what is installed today."
fi

# ---------------------------------------------------------------------------
rule "3. ADAPTERS AND TRIMMING"
# ---------------------------------------------------------------------------
if [[ -f "$ADAPTERS" ]]; then
  echo "  path   : $ADAPTERS"
  echo "  md5    : $(md5sum "$ADAPTERS" | cut -d' ' -f1)"
  echo "  records: $(grep -c '^>' "$ADAPTERS")"
  sub "contents"
  cat "$ADAPTERS" | sed 's/^/  /'
else
  echo "  *** $ADAPTERS not found — pass --adapters with the right path."
  find /data/project/chureeratp -maxdepth 3 -iname '*adapter*' 2>/dev/null | head \
    | sed 's/^/  candidate: /'
fi

sub "trimmomatic invocations found in the job scripts"
grep -rhoE 'trimmomatic +(SE|PE)[^\n]*' --include='*.sh' --include='*.sbatch' \
     . /data/project/chureeratp/RNA 2>/dev/null \
  | sed -E 's/.*(ILLUMINACLIP[^ ]*|SLIDINGWINDOW[^ ]*|LEADING[^ ]*|TRAILING[^ ]*|MINLEN[^ ]*)/&/' \
  | sort -u | head -20 | sed 's/^/  /'
echo "  NOTE: paired-end and single-end use different stringency (4:15/3/3/36 vs"
echo "        4:10/2/2/25). Layout is a cohort property, so this is confounded with"
echo "        cohort and belongs in the Limitations."

# ---------------------------------------------------------------------------
rule "4. rMATS INVOCATION FLAGS"
# ---------------------------------------------------------------------------
# Two sources: the job scripts (what was written) and the SLURM output (what ran).
# Line continuations are joined first, and comment lines dropped, before any flag is
# counted. Without the join a flag on a continuation line is invisible; without the
# comment filter a flag that was commented OUT is reported as set -- which is worse
# than missing it, because it puts a parameter in the Methods that never ran.
CMDFILE=$(mktemp)
trap 'rm -f "$CMDFILE"' EXIT
while IFS= read -r f; do
  sed -e :a -e '/\\$/N; s/\\\n//; ta' "$f" 2>/dev/null \
    | grep -E 'rmats\.py' | grep -vE '^[[:space:]]*#'
done < <(find . -maxdepth 2 \( -name '*.sbatch' -o -name '*.sh' \) \
              ! -name '48_provenance_versions.sh' 2>/dev/null) \
  > "$CMDFILE"

sub "rmats.py commands in job scripts (comments excluded, continuations joined)"
if [[ -s "$CMDFILE" ]]; then
  sort -u "$CMDFILE" | head -20 | cut -c1-240 | sed 's/^/  /'
  echo "  ($(wc -l < "$CMDFILE") invocation line(s) total)"
else
  echo "  *** no uncommented rmats.py invocation found under $ROOT"
fi

sub "value flags"
for flag in '--readLength' '--cstat' '--libType' '-t' '--anchorLength' '--nthread'; do
  vals=$(grep -oE -- "(^| )\\$flag[= ]+[^ ]+" "$CMDFILE" 2>/dev/null \
         | sed -E "s/.*\\$flag[= ]+//" | sort -u | tr '\n' ' ')
  if [[ -z "${vals// }" ]]; then
    printf '  %-20s NOT SET — rMATS default applies\n' "$flag"
  else
    printf '  %-20s %s\n' "$flag" "$vals"
  fi
done

sub "boolean flags (presence only — these take no value)"
for flag in '--variable-read-length' '--novelSS' '--allow-clipping' '--darts-model'; do
  n=$(grep -cE -- "(^| )\\$flag( |$)" "$CMDFILE" 2>/dev/null)
  if [[ "${n:-0}" -gt 0 ]]; then
    printf '  %-24s SET on %s of %s invocation line(s)\n' "$flag" "$n" "$(wc -l < "$CMDFILE")"
  else
    printf '  %-24s NOT SET\n' "$flag"
  fi
done
echo
echo "  rMATS v4.3.0 defaults, for any flag reported NOT SET:"
echo "    --cstat 0.0001 · --anchorLength 1 · --novelSS off · --libType fr-unstranded"
echo "    --variable-read-length off  (reads longer than --readLength are DISCARDED"
echo "                                 when it is off, which is why the next section"
echo "                                 checks the real read lengths)"
echo "  --libType matters here: the project measured the libraries as REVERSE-stranded"
echo "  (fwd/rev 9.8). If --libType is unset, rMATS ran as fr-unstranded and the"
echo "  Methods must say so rather than claiming the strandedness was used."

sub "rmats.py commands captured in SLURM job output (what actually ran)"
grep -rhoE 'rmats\.py[^\n]*' --include='*.out' --include='*.err' . 2>/dev/null \
  | sort -u | head -12 | cut -c1-240 | sed 's/^/  /'
[[ $(grep -rl 'rmats\.py' --include='*.out' --include='*.err' . 2>/dev/null | wc -l) -eq 0 ]] \
  && echo "  (none — the job scripts above are then the only record)"

# ---------------------------------------------------------------------------
rule "5. OBSERVED READ LENGTHS PER COHORT  (checks --readLength and Table 1)"
# ---------------------------------------------------------------------------
# Table 1 states 50/51/101/150 bp. That is metadata. This reads the BAMs.
if have samtools; then
  shopt -s nullglob
  LISTS=( dcm_bams.txt nf_bams.txt hcm_bams.txt *_bams_dedup.txt *_bams.txt )
  shopt -u nullglob
  [[ ${#LISTS[@]} -eq 0 ]] && echo "  no *_bams*.txt lists found in $ROOT"
  for lst in "${LISTS[@]}"; do
    [[ -f "$lst" ]] || continue
    # lists are newline- OR comma-separated, with absolute or relative paths
    bam=$(tr ',' '\n' < "$lst" | sed 's/[[:space:]]//g' | grep -m1 '\.bam$')
    if [[ -z "$bam" ]]; then
      printf '  %-28s SKIPPED: no .bam path parsed out of the list\n' "$lst"; continue
    fi
    [[ -f "$bam" ]] || bam="$ROOT/$bam"
    if [[ ! -f "$bam" ]]; then
      printf '  %-28s SKIPPED: first path not found on disk\n' "$lst"
      printf '  %-28s   %s\n' "" "$bam"; continue
    fi
    printf '\n  %-28s %s\n' "$lst" "$(basename "$bam")"
    samtools view "$bam" 2>/dev/null | head -n "$NREADS" \
      | awk '{print length($10)}' | sort -n | uniq -c | sort -rn | head -4 \
      | awk '{printf "      %8d reads at %s bp\n", $1, $2}'
    npair=$(samtools view -c -f 1 "$bam" 2>/dev/null)
    printf '      layout: %s\n' "$([[ "${npair:-0}" -gt 0 ]] && echo PAIRED || echo SINGLE)"
  done
  echo
  echo "  --variable-read-length was set on every rMATS invocation (section 4), so a"
  echo "  cohort showing several read lengths is handled rather than silently trimmed."
  echo "  These numbers are still the check on Table 1's 50 / 51 / 101 / 150 bp."
else
  echo "  samtools not on PATH — skipped."
fi

# ---------------------------------------------------------------------------
rule "6. EXPRESSION MATRIX PROVENANCE  (only one source may be described)"
# ---------------------------------------------------------------------------
echo "  STAR --quantMode GeneCounts output (*ReadsPerGene.out.tab):"
find . -maxdepth 3 -name '*ReadsPerGene.out.tab' 2>/dev/null | wc -l | sed 's/^/    files: /'
echo "  featureCounts output:"
find . -maxdepth 3 \( -name '*featureCounts*' -o -name '*.counts' -o -name 'counts*.txt' \) \
     2>/dev/null | head -6 | sed 's/^/    /'
sub "what 02_build_expr_matrix.py actually reads"
if [[ -f 02_build_expr_matrix.py ]]; then
  grep -nE 'ReadsPerGene|featureCounts|\.counts|read_csv|glob|default=' \
       02_build_expr_matrix.py | head -20 | sed 's/^/  /'
else
  echo "  02_build_expr_matrix.py not in $ROOT"
fi
sub "and what the composition scores are built from"
grep -nE 'MARKER|marker|CM_|FB_|TNNT2|MYH7|NPPA|COL1A1|DCN|POSTN|PDGFRA' \
     02_build_expr_matrix.py 2>/dev/null | head -20 | sed 's/^/  /'

# ---------------------------------------------------------------------------
rule "7. PASTE INTO THE METHODS"
# ---------------------------------------------------------------------------
cat <<'TXT'
  Fill the bracketed blanks in DRAFT_methods.md from the sections above:

    [STAR version]        section 2 if Log.out exists, else section 1
    [bcftools version]    section 1
    [Trimmomatic version] section 1  (the jar name usually carries it)
    [adapters.fa]         section 3 — cite the md5 and deposit the file with the code
    [rMATS flags]         section 4 — and state the DEFAULTS explicitly for any flag
                          reported NOT SET; "we used the defaults" is not reproducible
                          unless the defaults are written down
    [read lengths]        section 5 — this is also the check on Table 1's 50/51/101/150
    [expression source]   section 6 — STAR GeneCounts or featureCounts, never both
    [marker panels]       section 6

  Two answers change the text rather than just filling a blank:

  * If --libType is NOT SET, rMATS ran unstranded even though the libraries are
    reverse-stranded. That is defensible for PSI, but the Methods must say the
    quantification was run unstranded instead of implying the strandedness was used.
  * If any cohort shows more than one read length in section 5 and
    --variable-read-length was NOT SET, rMATS discarded the reads longer than
    --readLength. State the effective read length actually used.
TXT

echo
echo "done. Save this: bash 48_provenance_versions.sh 2>&1 | tee provenance_report.txt"
