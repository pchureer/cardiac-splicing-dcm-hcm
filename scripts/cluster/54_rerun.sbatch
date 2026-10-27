#!/bin/bash
#SBATCH --job-name=fdr_ladder
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --qos=ngs16
#SBATCH --output=logs/54_find_discovery_%j.out
set -euo pipefail
source activate rnaseq_tools
python 54_fdr_only_counts.py --find-discovery
