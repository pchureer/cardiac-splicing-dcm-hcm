#!/bin/bash
#SBATCH -p cpu
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH --qos=ngs16
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#
# Colorado DCM vs ICM -- independent replication of etiology-independence.
#
# Run 09a_get_colorado_icm_ids.sh FIRST and confirm 13 ICM ids with BAMs.
#
# TWO rMATS runs:
#   A) DCM vs ICM directly (37 vs 13) -- the clean test, no shared reference,
#      exactly the design that gave the DCM-vs-HCM null in MAGNet.
#   B) ICM vs NF (13 vs 14) -- lets you also ask whether the DCM-vs-NF
#      signature transfers to ICM-vs-NF, the parallel of the HCM analysis.
#
# Colorado is SINGLE-END, readLength 51 (per run_replication_rmats.sh).
# Getting -t wrong gives "USED: 0" -- check the log.

date
export PATH=/data/apps/anaconda3/envs/mamba/envs/rnaseq_tools/bin:$PATH
ROOT=/data/project/chureeratp/RNA/HEART
GTF=/data/project/chureeratp/hg38/ref/gencode.v49.annotation.gtf
cd "$ROOT"

[ -s ICM_Colorado.ids ] || { echo "run 09a first"; exit 1; }

build_list () {   # $1=idlist  $2=outfile
  > "$2"
  for s in $(cat "$1"); do
    ls "$ROOT/$s/"*Aligned.sortedByCoord.out.bam 2>/dev/null | head -1
  done | paste -sd, > "$2"
  echo "  $2: $(tr ',' '\n' < "$2" | grep -c .) bams (from $(grep -c . "$1") ids)"
}

echo "=== building bam lists ==="
build_list ICM_Colorado.ids Colorado_icm_bams.txt
[ -s Colorado_dcm_bams.txt ] || build_list DCM_Colorado.ids Colorado_dcm_bams.txt
[ -s Colorado_nf_bams.txt ]  || build_list NF_Colorado.ids  Colorado_nf_bams.txt

echo "  index any missing BAM indexes"
cat Colorado_icm_bams.txt | tr ',' '\n' | while read b; do
  [ -z "$b" ] && continue
  [ -f "$b.bai" ] || echo "$b"
done | xargs -r -P 16 -I{} samtools index {}

echo
echo "=== A) DCM vs ICM (the direct test) ==="
rm -rf "$ROOT/rmats_Colorado_DCMvsICM" "$ROOT/rmats_Colorado_DCMvsICM_tmp"
rmats.py \
  --b1 Colorado_dcm_bams.txt --b2 Colorado_icm_bams.txt \
  --gtf "$GTF" \
  -t single --readLength 51 --variable-read-length --nthread 16 \
  --od "$ROOT/rmats_Colorado_DCMvsICM" \
  --tmp "$ROOT/rmats_Colorado_DCMvsICM_tmp"

echo
echo "=== B) ICM vs NF ==="
rm -rf "$ROOT/rmats_Colorado_ICMvsNF" "$ROOT/rmats_Colorado_ICMvsNF_tmp"
rmats.py \
  --b1 Colorado_icm_bams.txt --b2 Colorado_nf_bams.txt \
  --gtf "$GTF" \
  -t single --readLength 51 --variable-read-length --nthread 16 \
  --od "$ROOT/rmats_Colorado_ICMvsNF" \
  --tmp "$ROOT/rmats_Colorado_ICMvsNF_tmp"

echo
echo "=== rMATS read outcomes (USED must be large; USED: 0 = wrong -t) ==="
grep -h -E "USED|NOT_" 09b_colorado_icm_rmats.sbatch_*.out 2>/dev/null | tail -6

echo
echo "=== output sizes ==="
for d in rmats_Colorado_DCMvsICM rmats_Colorado_ICMvsNF; do
  for ET in SE MXE; do
    f=$d/${ET}.MATS.JCEC.txt
    [ -s "$f" ] && printf "  %-45s %s rows\n" "$f" "$(( $(wc -l < $f) - 1 ))" \
                || printf "  %-45s MISSING\n" "$f"
  done
done

echo
echo "NEXT (Colorado has no age/sex/RIN, so use --no-covariates):"
echo "  printf '%s\\tDCM\\n' \$(cat DCM_Colorado.ids) >  colorado_dcm_icm.cov"
echo "  printf '%s\\tICM\\n' \$(cat ICM_Colorado.ids) >> colorado_dcm_icm.cov"
echo "  for ET in SE MXE; do"
echo "    python phase2_diffsplice_v2.py \\"
echo "      --mats rmats_Colorado_DCMvsICM/\${ET}.MATS.JCEC.txt --event-type \$ET \\"
echo "      --b1 Colorado_dcm_bams.txt --b2 Colorado_icm_bams.txt \\"
echo "      --cov colorado_dcm_icm.cov --group1 DCM --group2 ICM \\"
echo "      --min-samples 10 --no-covariates --out \${ET}_Colorado_DCMvsICM"
echo "  done"
date
