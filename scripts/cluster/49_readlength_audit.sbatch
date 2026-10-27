#!/bin/bash
#SBATCH -J rlaudit
#SBATCH -p cpu
#SBATCH -c 2
#SBATCH --mem=4G
#SBATCH --qos=ngs16
#SBATCH -t 010:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#
# Read-length audit as a batch job. It is I/O bound, not CPU bound, so two cores
# and 4 GB are plenty -- the cost is decompressing the head of each bam.
#
#   sbatch 49_readlength_audit.sbatch              # 40 bams per list (default)
#   sbatch 49_readlength_audit.sbatch --max-samples 0   # every sample
#
date
export PATH=/data/apps/anaconda3/envs/mamba/envs/rnaseq_tools/bin:$PATH
cd /data/project/chureeratp/RNA/HEART
bash 49_readlength_audit.sh "$@" 2>&1 | tee readlength_audit.txt
date
