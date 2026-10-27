#!/usr/bin/env bash
# 09a_get_colorado_icm_ids.sh -- run on the LOGIN NODE (needs network).
#
# Your Colorado cohort is PRJNA477855 = GSE116250 (Sweet et al. 2018):
# 64 LV samples = 14 non-failing + 37 DCM + 13 ISCHEMIC cardiomyopathy.
# You are using the 37 DCM and 14 NF. The 13 ICM are sitting unused.
#
# They give you a second etiology pair, in a second cohort, same batch and
# same prep -- an independent test of whether the splicing signature is
# etiology-independent. No public cohort has DCM+HCM in one batch, so this
# is the only replication available.
#
# METHOD: rather than parsing free-text disease labels (fragile), we take
# every run in the project and subtract the DCM and NF lists you already
# have. The remainder must be the ICM samples. Then we sanity-check the
# count and, if titles are available, confirm they say ischemic.

set -uo pipefail
ROOT=/data/project/chureeratp/RNA/HEART
cd "$ROOT"

PRJ=PRJNA477855
META=colorado_filereport.tsv

echo "=== 1. Fetch the run list for $PRJ from ENA ==="
if [ ! -s "$META" ]; then
  curl -sS -o "$META" \
    "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${PRJ}&result=read_run&fields=run_accession,sample_accession,sample_title,sample_alias,library_layout,read_count&format=tsv&limit=0" \
    || { echo "download failed -- fetch it in a browser and save as $META"; exit 1; }
fi
echo "rows: $(( $(wc -l < "$META") - 1 ))   (expect 64)"
head -1 "$META"

echo
echo "=== 2. What you already have ==="
for f in DCM_Colorado.ids NF_Colorado.ids; do
  [ -s "$f" ] && printf "  %-20s %s\n" "$f" "$(grep -c . $f)" || echo "  $f MISSING"
done

echo
echo "=== 3. ICM = everything else ==="
tail -n +2 "$META" | cut -f1 | grep . | sort -u > .all_colorado.ids
cat DCM_Colorado.ids NF_Colorado.ids | tr -d '\r' | grep . | sort -u > .known_colorado.ids
comm -23 .all_colorado.ids .known_colorado.ids > ICM_Colorado.ids
echo "  total runs in project : $(grep -c . .all_colorado.ids)"
echo "  DCM + NF already used : $(grep -c . .known_colorado.ids)"
echo "  remaining (= ICM)     : $(grep -c . ICM_Colorado.ids)   <-- expect 13"

echo
echo "  any of your DCM/NF ids NOT in the project? (should be none):"
comm -13 .all_colorado.ids .known_colorado.ids | sed 's/^/    /' | head

echo
echo "=== 4. Confirm the remainder really is ischemic ==="
TITLECOL=$(head -1 "$META" | tr '\t' '\n' | grep -n "sample_title" | cut -d: -f1)
if [ -n "$TITLECOL" ]; then
  echo "  sample_title for the inferred ICM runs:"
  awk -F'\t' -v c="$TITLECOL" 'NR==FNR{k[$1];next} $1 in k {print "    "$1"\t"$c}' \
      ICM_Colorado.ids "$META" | head -20
  echo
  echo "  titles of your DCM list, for contrast:"
  awk -F'\t' -v c="$TITLECOL" 'NR==FNR{k[$1];next} $1 in k {print "    "$1"\t"$c}' \
      DCM_Colorado.ids "$META" | head -5
fi
rm -f .all_colorado.ids .known_colorado.ids

echo
echo "=== 5. Do the BAMs exist? ==="
miss=0
while read -r s; do
  [ -z "$s" ] && continue
  ls "$ROOT/$s/"*Aligned.sortedByCoord.out.bam >/dev/null 2>&1 || { echo "  no BAM: $s"; miss=$((miss+1)); }
done < ICM_Colorado.ids
echo "  ICM samples without a BAM: $miss"
if [ "$miss" -gt 0 ]; then
  echo "  >>> These were never downloaded/aligned. You will need to pull them"
  echo "      from ENA and run them through the same STAR pipeline first."
fi

echo
echo "If the count is 13, the titles say ischemic, and 0 BAMs are missing,"
echo "you are ready for: sbatch 09b_colorado_icm_rmats.sbatch"
