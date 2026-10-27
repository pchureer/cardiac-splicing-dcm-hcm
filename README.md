# Alternative splicing in dilated and hypertrophic cardiomyopathy — analysis code and event-level results

This repository holds the complete analysis code and the event-level result tables for
the manuscript. It is written for a referee or a reader who wants to check a specific
number, so it is organised around **where each published number comes from** rather than
around the order the code was written in.

Nothing here is a pipeline you run end to end. The alignment and rMATS steps need a
SLURM cluster and ~550 GB of primary BAM files that are not redistributable; everything
downstream of rMATS runs from the deposited tables in minutes.

---

## 1. Read this first: the four things most likely to trip you up

1. **The exon-level rule is collapse-first.** Collapsing rMATS rows to event coordinates
   and *then* applying the significance filter gives 813 significant skipped exons in
   the discovery contrast; filtering first gives 864. At 50 bp the same pair is 926
   against 991. Every number in the paper is collapse-first.
2. **`dPSI_adj`, not `dPSI_raw`, carries every published effect size.** The boolean
   significance columns in the phase-2 tables were written against `dPSI_adj` — the
   exact flag value is recorded in `PROVENANCE.md`, measured from the producing script
   rather than remembered.
3. **`caution` is not `exclude`.** The `reliability` column has three values and the
   paper uses the `caution` rows. Code that keeps only `reliability == "reference"`
   drops data the paper reports.
4. **125 bp is primary.** The analysis was first run at a 50 bp read length and later
   re-quantified at 125 bp, which is now the reported baseline. Both sets of tables are
   deposited. Files and directories carrying `_rl125` are primary; the unsuffixed
   equivalents are the superseded 50 bp run, kept so the read-length comparison in the
   paper can be reproduced. **If a file exists in both forms, the `_rl125` one is the
   published one.**

`CONVENTIONS.md` states all of this formally, with the code that implements it. Read it
before reproducing any count.

---

## 2. What is in here

```
scripts/            numbered analysis scripts, 00-170
scripts/patches/    scripts whose only job is to edit another script in place
scripts/cluster/    SBATCH submission scripts and shell drivers
tables/             derived result tables (small; the ones the paper's numbers are read from)
figures/            the figure PNGs as submitted
events/             the full rMATS event tables (Zenodo only - see section 6)
environment.yml     the conda environment the analysis ran in
MANIFEST.tsv        every deposited file with its size and SHA-256
PROVENANCE.md       facts measured from the code at packaging time, not asserted by hand
CONVENTIONS.md      the analysis conventions - read this before reproducing a count
DATA_AVAILABILITY.md  the submission statements, and what is and is not redistributed
```

Two of the deposited scripts are about the deposit itself.
`scripts/172_verify_regeneration.py` re-runs every figure and table producer and
compares the result to the submitted artefact — by decoded pixels for a figure, since
matplotlib writes a timestamp into PNG metadata and byte-identity would report every
figure as changed. It backs each artefact up before running anything and restores it if
the producer fails.

The tree was assembled by `scripts/171_build_deposit.py`, which is itself deposited. It
refuses to package if a genotype-derived file would be copied, if the two names in
section 3 are not byte-identical, if the script numbering has an undeclared gap, or if
any documentation file still carries a placeholder. `PROVENANCE.md` is its output, not a
hand-written file.

### The script numbering has three gaps

There is no script **94**, **153** or **154**. These are not missing files. 94 was
written against a premise that turned out to be void and was deleted rather than left to
mislead; 153 and 154 were superseded before they produced anything used. A hole in a
numbering reads as an omission, so it is stated here.

Letter suffixes (`03b`, `09a`, `09b`, `15a`, `15b`) are second passes over the same step,
not alternatives — the lettered script runs after the bare number.

### Many of the scripts are patches

Roughly half of `scripts/` are `*_patch_*` scripts that edit an earlier script in place.
This is how the 50 bp → 125 bp re-baselining was done: rather than forking each producer,
each was patched, with the patch recorded as its own numbered, reviewable file. They are
in `scripts/patches/` and are **not re-runnable** — the scripts in `scripts/` are already
in their patched, final state. The patches are deposited because they are the audit trail
for how a plotted or printed constant changed.

`.bak` backup files written by those patches are **not** deposited (62 of them). They are
intermediate states, not results.

---

## 3. ⚠ Two naming traps in the rMATS output directories

These cost time if you assume the names mean what they look like.

| directory | what it actually is |
|---|---|
| `rmats_magnet_rl125` | MAGNet, **DCM vs NF** — the discovery contrast |
| `rmats_HCMvsNF_rl125` | MAGNet, **HCM vs NF** |
| `rmats_DCMvsHCM_rl125` | MAGNet, **DCM vs HCM** — the negative contrast |

The first is named for the cohort and the other two for the contrast. All three are the
same cohort. The asymmetry is historical and was not renamed, because renaming a
directory after the tables were written would have broken the provenance recorded in
the scripts.

**The `--group1` trap.** In the HCM vs NF contrast, `--group1` is *HCM*, so a positive
`dPSI` means higher inclusion in HCM. Several scripts take `--group1` on the command
line and the sign of every effect depends on it. Scripts **54, 62, 65, 69 and 77**
therefore verify the contrast from the file's own contents rather than trusting the
filename or the flag — if you adapt them, keep that check.

### Two live names for the same discovery table — one is a symlink

On the analysis machine `SE_rl125.all_events.tsv` is a **symbolic link** to
`SE_v2filtered_rl125.all_events.tsv`. They are not two copies that happen to agree;
they are one file reached by two names. Both names are live in the scripts because the
second was introduced mid-analysis and the first was never retired.

**The deposit carries the file once**, under the real name
`SE_v2filtered_rl125.all_events.tsv`. A script that asks for `SE_rl125.all_events.tsv`
is asking for the same bytes — recreate the link, or point the script's `--discovery`
flag at the real name:

```bash
ln -s SE_v2filtered_rl125.all_events.tsv SE_rl125.all_events.tsv
```

⚠ This also means **a hash comparison between the two names proves nothing** — it reads
the same bytes twice. `171_build_deposit.py` reports the link rather than asserting
byte-identity, because the gate that asserted it passed by construction.

### Superseded sample lists

`*_bams_dedup.txt` files, if you encounter them in an older copy of this analysis, are
**not deposited** and were never used. They came from an earlier exclusion strategy that
dropped both members of each technical-replicate pair, and so show 35/12/5 where the
paper reports 36/13/6. rMATS was run on the full per-cohort lists and exclusions were
applied at the modelling stage, recorded in `geno/EXCLUDE_samples.tsv`. The
technical-replicate tie-break was: **keep the member with more called genotype sites.**

---

## 4. Where each published number comes from

| artefact (Scientific Reports numbering) | producer | reads |
|---|---|---|
| Table 1 (lead genes) | `61_*`, `73_*` (patched by `128_`) | `tables/Table2_rl125.tsv` (file keeps its build name) |
| Table 2 (replication) | `98_*` (patched by `132_`, caption by `170_`) | `tables/Table4_rl125.tsv` |
| Supplementary Table 11 (cohorts; was main-text Table 1) | `42_*` (patched by `149_`) | `tables/Table1.tsv` |
| Supplementary Table 12 (etiology; was main-text Table 3) | `87_*` (patched by `168_`) | `tables/t3_rl125/` |
| Figure 1 | `43_*` (patched by `169_`) | slopes passed in, see `PROVENANCE.md` |
| Figure 2 | `39_*` | `tables/master_f2_panelB_rl125.tsv` (built by `164_`) |
| Figure 3 | `37_*` (patched by `162_`), panels B/C by `158_` | event tables |
| Figure 4 | `38_*` (patched by `162_`, `163_`) | cohort slopes |
| Figure 5 | `90_*` (patched by `155_`, `156_`, cosmetics `176_`) | event tables |
| Figure 6 | `41_*` (patched by `165_`, `166_`, cosmetics `176_`) | g:Profiler output in `tables/` |
| Supplementary Table 1 and 3 | `76_*` then `95_* --job s1c` (95 adds the basis column 76 strips) | `tables/S1a.tsv`, `tables/S1c_markers.tsv` |
| Supplementary Table 2 | `69_*` | `tables/S1b.tsv` |
| Supplementary Table 4 | `81_*` (patched by `178_`) | `tables/S2_rl125.tsv` |
| Supplementary Table 5 | `62_*` (run with all three parts; `--parts X` overwrites the full file) | `tables/S5_rl125.tsv` |
| Supplementary Table 6 | `83_*` — **50 bp by decision**, see its caption | `tables/S6.tsv` |
| Supplementary Table 7 | `65_*` with `--fwd/--inv` suffixed and the three `S7_*` checks | `tables/S7_rl125.tsv` |
| Supplementary Table 8 | (a) `183_*.sbatch` → `184_*`; (b) `180_*` | `tables/S8a_permutations_rl125.tsv`, `tables/S8b_reverse_replication_rl125.tsv` |
| Supplementary Table 9, Figure 1 | scripts `27_`–`31_` (genotype-level; summary only is deposited) | `tables/S9*.tsv`, `figures/S9.png` |
| Supplementary Table 10, Data 1 | `77_*` then `95_* --job s3 --fdr-col FDR_rank` | `tables/SupplementaryTable10_lvef_summary_rl125.tsv`, `events/` |
| Supplementary Figure 2 | `179_*` (imports `82_*`'s estimator; refuses unless the recorded slopes reproduce) | `events/mp_SE_*_rl125` |
| Supplementary Information file | `181_*` from `Supplementary_Captions_rl125.md`, `Supplementary_Methods_rl125.md`, `Supplementary_Discussion_rl125.md` | the tables above |

Every row above was regenerated from its command on 26 Oct 2026 and compared byte-for-byte
with the deposited file (`scripts/172_verify_regeneration.py` with `verify_spec.tsv`, which is
deposited in `tables/`; PDFs compared with their creation date stripped). The environment
variables each command needs are in that spec, not in this README, so they cannot drift from
what was actually run. `182_stale_number_sweep.py` lists every 50 bp value the record pairs with
a 125 bp replacement and searches a document for them; the manuscript and this deposit's text
were swept with it.

**Five of the six figure producers passed all of their own internal gates while drawing
at least one stale constant.** The gates checked inputs and computations but never a
plotted or printed literal. The patch scripts listed above are what fixed that, and
several of them replaced a hard-coded number with a mandatory command-line or
environment override that *refuses to plot* if the data named is 125 bp and the override
is absent. If you adapt a figure producer, keep those guards — they are the only thing
standing between a re-baselining and a figure that silently shows the old numbers.

---

## 5. What is deliberately not deposited, and why

**Genotype data.** The sample-identity screen called genotypes at 3,820 SNVs across 670
samples and produced a per-pair kinship matrix carrying sample identifiers. **Neither the
genotype VCF nor the per-pair kinship table is deposited.** A genotype matrix is
individually identifying in principle, and the source datasets' terms may restrict
redistribution of individual-level data. What *is* deposited is the complete code that
produces them from the primary BAMs, and the summary the paper reports: the kinship
ranges and the 24 identified cross-dataset duplicate pairs. Anyone with authorised access
to the primary data can regenerate the full output from the deposited scripts.

**Primary sequencing data.** Not ours to redistribute. All of it is public — see
`DATA_AVAILABILITY.md` for the accessions.

**Intermediate and scratch files.** 54 `.bak` patch backups, 2 `.pid` job-tracking
files, 155 run logs, a 61-file probe/scratch directory, and 5 superseded
`*_bams_dedup.txt` sample lists (section 3).

**Note on what the genotype exclusion does *not* cover.** The exclusion applies to
genotype *output*, never to code. `27_genotype_concordance.sbatch` is the script that
runs the screen and it is deposited, as section 1 promises. A gene-panel `.bed` is an
analysis input and is deposited; a `.bed` sitting beside a `.bim`/`.fam` is a PLINK
genotype triple and is not. The reported 24 duplicate pairs are deposited as the
published summary; the full pairwise kinship matrix they were drawn from is not.

---

## 6. The split between GitHub and Zenodo

The two targets hold the same tree with one difference:

- **GitHub** carries the code, the small derived tables, the figures and all the
  documentation. It is the readable, diffable record. `events/` is in `.gitignore`.
- **Zenodo** carries all of that *plus* `events/` — approximately 550 MB of full rMATS
  event tables. These are too large for either journal's supplementary-material limit
  (Scientific Reports caps supplementary files at 50 MB), so the Zenodo record, not the
  supplement, is where the event-level data lives. The DOI is what the manuscript cites.

`MANIFEST.tsv` covers both; the `in_git` column says which files GitHub carries.

---

## 7. Reproducing a number

```bash
conda env create -f environment.yml && conda activate rnaseq_tools

# the discovery count: 53,456 tested, 813 significant
python scripts/103_recount_rl125.py --dir events/

# Table 3
python scripts/87_*.py --dir tables/t3_rl125/
```

Most scripts print a calibration pass against a frozen record before they report a new
number, and **exit non-zero rather than report** if the frozen number does not reproduce.
A clean run is one that prints its calibration block and then its result. If a script
exits with `[fatal]`, read the message — it is telling you which assumption failed, and
it was put there deliberately after that assumption failed once.
