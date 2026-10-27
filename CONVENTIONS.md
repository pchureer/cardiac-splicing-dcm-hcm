# Analysis conventions

Every count in the paper follows the rules below. They are stated here once
because **the most likely way to fail to reproduce this analysis is to apply
them in the wrong order.**

## 1. The significance rule

An event is significant when

```
FDR < 0.05   AND   |dPSI_adj| >= 0.10
```

applied **within a cohort**, never to a pooled value.

## 2. ⚠ COLLAPSE FIRST, THEN FILTER

rMATS emits one row per transcript context, so the same cassette exon appears
several times. The analysis collapses to the event's own coordinates **before**
applying the significance rule, keeping the lowest-FDR row:

```python
NKEEP = {"SE": 3, "MXE": 5, "A3SS": 5, "A5SS": 5, "RI": 3}

ekey = ":".join(str(key).split(":")[:NKEEP[event_type]])
d = df.sort_values("FDR").drop_duplicates(ekey, keep="first")   # collapse
sig = d[(d.FDR < 0.05) & (d.dPSI_adj.abs() >= 0.10)]            # then filter
```

**Reversing the order changes the answer.** For the discovery contrast at
125 bp, collapse-first gives **813** significant skipped exons and filter-first
gives **864**. At 50 bp the same pair is 926 against 991. Every number in the
paper is collapse-first.

### NKEEP counts FIELDS of `key`, not rows

The `key` column is seven colon-separated fields (chromosome plus six
coordinates). `NKEEP` is how many of those fields define the event's identity,
which differs by class because the classes have different coordinate
structures. Truncating to a character count, or to `exonStart`/`exonEnd` alone,
is wrong for MXE, A3SS and A5SS.

## 3. Exon level, not full coordinate

"Exon-level" means after the collapse above. "Full-coordinate" means the
uncollapsed rows. Both appear in the deposited tables and they are not
interchangeable:

| discovery contrast, 125 bp | tested | significant |
|---|---|---|
| exon-level (reported) | 53,456 | 813 |
| full-coordinate | 80,973 | 1,074-equivalent |

Supplementary Table 6 (the control-group size sweep) is the one analysis
computed full-coordinate, because it predates the convention being locked. Its
reference arm holds 752 events full-coordinate against 665 exon-level.
**Re-reporting it at exon level requires recomputation, not rescaling.**

## 4. The `reliability` column

`master_numbers_flagged.tsv` carries a `reliability` column with three values,
and **`caution` is not `exclude`**:

| value | meaning |
|---|---|
| `reference` | used in the paper |
| `caution` | flagged, still used |
| `exclude` | removed from every table and figure |

Colorado MXE is `exclude` (reverse Deming slope 0.060 against an internal
control of 0.730). Fuwai MXE is `caution`. Code that filters on
`reliability != "reference"` will silently drop the `caution` rows, which the
paper does use.

## 5. dPSI columns

The standard 20-column phase-2 table carries **both** `dPSI_raw` and
`dPSI_adj`. **`dPSI_adj` is the covariate-adjusted estimate and is the one
every count in the paper uses.**

The DCM-versus-HCM tables are a different 12-column product with a single
effect column named for the contrast (`dPSI_DCMvsHCM`) and `se` rather than
`se_adj`. There is no raw/adjusted pair to choose between in those files.
Substituting an unadjusted effect for an adjusted one would change a null.

## 6. Three conventions that were not recorded at the time

Recovered by re-reading the producers, and stated here because nothing in the
file names carries them:

1. **Supplementary Table 6's reference arm is draw 1 of 5** from the
   split-control sweep, not a mean over the five draws.
2. **Replicate spread in that sweep is the sample standard deviation**
   (`ddof = 1`), not the standard error.
3. **`median_samples_g1` / `median_samples_g2` in the master table are
   per-event medians after the coverage mask**, not cohort group sizes. For
   MAGNet they read 125/131 against true group sizes of 166/166. Use Table 1
   for group sizes.

## 7. mde80

`mde80 = 2.802 x se_adj` is the minimum effect detectable at 80% power and is
present as a column in the phase-2 tables. It is not recomputed downstream.
