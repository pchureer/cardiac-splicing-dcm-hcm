#!/bin/bash
#SBATCH -p cpu
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH --qos=ngs16
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#
# Build the MAGNet expression matrix + cell-composition scores.
# Submit:  sbatch 02_build_expr_matrix.sbatch
# Watch :  tail -f 02_build_expr_matrix_<jobid>.err

date
export PATH=/data/apps/anaconda3/envs/mamba/envs/rnaseq_tools/bin:$PATH
ROOT=/data/project/chureeratp/RNA/HEART
GTF=/data/project/chureeratp/hg38/ref/gencode.v49.annotation.gtf
cd "$ROOT"

cat DCM_MAGNet.ids NF_MAGNet.ids > magnet_all.ids
echo "MAGNet samples: $(wc -l < magnet_all.ids)   (expect 332)"

python 02_build_expr_matrix.py \
  --heart-dir "$ROOT" \
  --ids magnet_all.ids \
  --source star \
  --gtf "$GTF" \
  --out-prefix magnet_expr

echo "--- outputs ---"
ls -lh magnet_expr.log2cpm.tsv.gz magnet_expr.per_sample.tsv
date
