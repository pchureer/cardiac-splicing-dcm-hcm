#!/bin/bash
#SBATCH -p cpu
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH --qos=ngs16
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --array=1-332%20
#
# YOU PROBABLY DO NOT NEED THIS JOB.
# 00_check_expression_inputs.sh found 716 ReadsPerGene.out.tab files, so STAR
# already produced gene counts. Run 02 directly. Keep this only as a fallback
# if 02 reports samples with no counts file.
#
# TEST ONE FIRST:
#   mkdir -p fc_counts
#   cat DCM_MAGNet.ids NF_MAGNet.ids > magnet_all.ids
#   sbatch --array=1-1 01_featurecounts_array.sbatch
#   # check fc_counts/<sample>.counts.txt.summary -> Assigned >50% of total
# THEN scale:
#   sbatch --array=1-$(wc -l < magnet_all.ids)%20 01_featurecounts_array.sbatch

date
export PATH=/data/apps/anaconda3/envs/mamba/envs/rnaseq_tools/bin:$PATH
ROOT=/data/project/chureeratp/RNA/HEART
GTF=/data/project/chureeratp/hg38/ref/gencode.v49.annotation.gtf
IDS=$ROOT/magnet_all.ids
OUT=$ROOT/fc_counts
cd "$ROOT"
mkdir -p "$OUT"

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$IDS")
[ -z "$SAMPLE" ] && { echo "no sample at line $SLURM_ARRAY_TASK_ID"; exit 1; }

# RAW STAR bam -- NOT the bqsr bam
BAM=$(ls "$ROOT/$SAMPLE/"*Aligned.sortedByCoord.out.bam 2>/dev/null | head -1)
[ -z "$BAM" ] && { echo "no raw STAR bam for $SAMPLE"; exit 1; }

if [ -s "$OUT/$SAMPLE.counts.txt" ]; then
  echo "$SAMPLE already done, skipping"; exit 0
fi

# MAGNet is PAIRED 50bp. -s 0 = unstranded, confirmed from the ReadsPerGene
# columns (N_noFeature 3.5M unstranded vs 29.2M/29.1M stranded).
featureCounts \
  -T 4 \
  -p --countReadPairs \
  -t exon -g gene_id \
  -s 0 \
  -a "$GTF" \
  -o "$OUT/$SAMPLE.counts.txt" \
  "$BAM"

echo "done $SAMPLE"
date
