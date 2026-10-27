# Provenance

Every statement below was measured from the code and the files at
packaging time by `171_build_deposit.py`. Nothing here is typed in by
hand. A line reading `NOT MEASURED` means the script could not find
the evidence -- treat it as an open question, not as a null result.

## The `--dpsi-rule` that wrote the boolean significance columns

Declared as an argparse **default**:

- `any` — in 161_patch_35_combined_readlength.py, 35_combined_evidence.py

Seen **passed explicitly** on a command line inside a script:

- `any` — in 84_pool_b.py
- `discovery` — in 84_pool_b.py
- `that` — in 171_build_deposit.py

⚠ **3 different values appear: ['any', 'discovery', 'that'].** The boolean columns in the deposited tables were written by exactly one of them. Resolve which, and state it in README section 1 item 2 — two live rules is how a reader reproduces 864 instead of 813.

Mentioned in 4 script(s): 161_patch_35_combined_readlength.py, 171_build_deposit.py, 35_combined_evidence.py, 84_pool_b.py

## Reference genome and annotation

```
GTF   /data/project/chureeratp/hg38/ref/gencode.v49.annotation.gtf      (seen in 02_build_expr_matrix.py)
GTF   args.gtf      (seen in 02_build_expr_matrix.py)
GTF   ROOT/ref/gencode.v49.annotation.gtf      (seen in 100_requant_magnet_125.sh)
FASTA .bim/.fa      (seen in 171_build_deposit.py)
FASTA font.fa      (seen in 37_figure3_etiology.py)
FASTA figure.fa      (seen in 37_figure3_etiology.py)
FASTA savefig.fa      (seen in 37_figure3_etiology.py)
FASTA axes.fa      (seen in 66_build_s9.py)
FASTA /path/adapters.fa      (seen in 48_provenance_versions.sh)
FASTA /data/project/chureeratp/fastq/adapters.fa      (seen in 48_provenance_versions.sh)
FASTA adapters.fa      (seen in 48_provenance_versions.sh)
```

⚠ More than one annotation path appears in the scripts. Confirm which one the primary run used; a second annotation in the tree is how a cohort ends up quantified against different exon definitions.

## Read length passed to rMATS

- `125` — 3 script(s): 100_requant_magnet_125.sh, 125_patch_23_suffix_and_guard.py, 49_readlength_audit.sh
- `50` — 9 script(s): 148_patch_63_readlength_arms.py, 149_patch_42_table1_readlength.py, 42_table1.py, 49_readlength_audit.sh, 50_readlength_impact.py, 51_rmats_magnet_rl125.sbatch, 52_rl_probe.sbatch, 62_build_s5.py ...
- `51` — 1 script(s): 09b_colorado_icm_rmats.sbatch
- `variable-read-length` — 10 script(s): 09b_colorado_icm_rmats.sbatch, 100_requant_magnet_125.sh, 149_patch_42_table1_readlength.py, 15b_hcm_cohorts_rmats.sbatch, 171_build_deposit.py, 42_table1.py, 48_provenance_versions.sh, 49_readlength_audit.sh ...

⚠ Two read lengths are live in the deposited scripts (125, 50, 51 bp). That is expected — the 50 bp run is deposited so the read-length comparison can be reproduced — but confirm that every script reading a `_rl125` table is the 125 bp one.

## The two live names for the discovery table

`SE_rl125.all_events.tsv` and `SE_v2filtered_rl125.all_events.tsv` are byte-identical.

    SHA-256  symlink:SE_v2filtered_rl125.all_events.tsv

## Environment

Packaged with Python 3.11.14 on cuaim-fat-01.cuaim.hpc.
⚠ `conda env export -n rnaseq_tools` did not succeed, so `environment.yml` was NOT written. Run the packaging from inside the environment, or export it by hand -- a deposit without a pinned environment is not reproducible.
