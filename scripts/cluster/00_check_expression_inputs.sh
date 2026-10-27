#!/usr/bin/env bash
# 00_check_expression_inputs.sh
# Run this FIRST. It decides whether you need featureCounts at all.
# Nothing here is expensive; it only inspects files.

set -euo pipefail
HEART=/data/project/chureeratp/RNA/HEART
GTF=/data/project/chureeratp/hg38/ref/gencode.v49.annotation.gtf

echo "=== 1. Do STAR gene counts already exist? ==="
# If your STAR runs used --quantMode GeneCounts you already have expression for free.
n_rpg=$(find "$HEART" -maxdepth 2 -name "*ReadsPerGene.out.tab" 2>/dev/null | wc -l)
echo "ReadsPerGene.out.tab files found: $n_rpg"
if [ "$n_rpg" -gt 0 ]; then
  echo "--- example file, first 8 lines ---"
  find "$HEART" -maxdepth 2 -name "*ReadsPerGene.out.tab" | head -1 | xargs head -8
  echo "NOTE: cols are gene / unstranded / fwd-stranded / rev-stranded."
  echo "      02_build_expr_matrix.py auto-detects which column to use."
fi

echo
echo "=== 2. How many MAGNet samples do we need? ==="
for f in DCM_MAGNet.ids NF_MAGNet.ids HCM_MAGNet.ids; do
  if [ -f "$HEART/$f" ]; then printf "%-18s %s\n" "$f" "$(wc -l < "$HEART/$f")"; fi
done

echo
echo "=== 3. Coverage check: how many of those samples have a counts file? ==="
if [ "$n_rpg" -gt 0 ]; then
  miss=0
  while read -r s; do
    [ -z "$s" ] && continue
    if ! ls "$HEART/$s/"*ReadsPerGene.out.tab >/dev/null 2>&1; then miss=$((miss+1)); fi
  done < <(cat "$HEART/DCM_MAGNet.ids" "$HEART/NF_MAGNet.ids")
  echo "MAGNet samples MISSING ReadsPerGene.out.tab: $miss"
  echo "(0 => skip featureCounts entirely, go straight to 02_build_expr_matrix.py)"
fi

echo
echo "=== 4. Is featureCounts available (fallback path)? ==="
which featureCounts 2>/dev/null || echo "featureCounts NOT in PATH -> mamba install -n rnaseq_tools -c bioconda subread"

echo
echo "=== 5. Sanity: RAW STAR bam present (NOT bqsr) for one sample? ==="
s=$(head -1 "$HEART/DCM_MAGNet.ids")
ls -lh "$HEART/$s/"*Aligned.sortedByCoord.out.bam 2>/dev/null || echo "no raw STAR bam for $s -- check naming"

echo
echo "=== 6. RBM20 is present in the GTF? ==="
grep -m1 -P '\tgene\t' "$GTF" | head -0 || true
awk -F'\t' '$3=="gene" && $9 ~ /gene_name "RBM20"/ {print $1,$4,$5,$7; exit}' "$GTF"
