#!/bin/bash
#SBATCH -p cpu
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH --qos=ngs16
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#
# INVERTED RBP SCAN -- does ANY candidate regulator track the DISCOVERY SIGNATURE?
#
# WHAT WAS ALREADY DONE, AND WHY IT IS NOT THIS
# ---------------------------------------------
# Script 03 tested RBM20 against the nine RBM20-target genes, and ran MBNL1, CELF1,
# SRSF1, PTBP1, QKI and GAPDH through the identical pipeline as SPECIFICITY CONTROLS --
# i.e. "does QKI correlate with RBM20's targets more than with background?" All were
# null (p>0.2 cell-adjusted), which is the correct control for the RBM20 claim.
#
# It is NOT a test of whether QKI regulates OUR signature. The signature is ~1,074
# significant SE events across hundreds of genes, not the nine-gene RBM20 program.
# Nobody has asked whether any candidate RBP tracks THOSE events.
#
# This script inverts the design: the composition-adjusted signature genes become the
# target set, and RBM20 is demoted to just another candidate.
#
# THE SILENT-NULL TRAP THIS SCRIPT EXISTS TO AVOID
# ------------------------------------------------
# magnet_expr.per_sample.tsv holds only RBM20 plus the control genes step 02 was
# ORIGINALLY run with. Script 03 selects candidates with `if g in meta.columns`, so a
# gene missing from that file is dropped WITHOUT ERROR and reported as no result --
# indistinguishable from a real null. RBM5, ZRANB2, HNRNPF and HNRNPA0 are not in the
# file. Step A below pulls them from the full log2CPM matrix and FAILS LOUDLY if any is
# absent, so "not expressed above threshold" can never be mistaken for "tested, null".
#
# CANDIDATES AND WHY EACH IS THERE
#   QKI                     proposed for exactly this signature (Chen 2021 Nat Commun;
#                           Furtado 2026 Cardiovasc Res)
#   RBM5 ZRANB2 HNRNPF      the four factors Li 2024 (Heliyon e35408) selected by
#   HNRNPA0                 LASSO + SVM-RFE -- from a study that used our Colorado cohort
#   MBNL1 CELF1 SRSF1 PTBP1 established cardiac/developmental splicing regulators
#   RBM20                   auto-included by the script; now a candidate, not the hypothesis
#   GAPDH                   negative control. It is what exposed the technical inflation
#                           the first time round -- if GAPDH scores, the model is broken
#
# READ THE "DCM" (WITHIN-DISEASE) MODEL, NOT THE ALL-SAMPLES ONE. Target events were
# selected for large DCM-vs-NF dPSI, so in an all-samples model ANY gene whose expression
# differs by disease correlates with them simply because both track disease. Section 4
# already restricts to the within-DCM model; do not quote the ALL model here.
#
# Submit: sbatch 23_signature_rbp_scan.sbatch
date
export PATH=/data/apps/anaconda3/envs/mamba/envs/rnaseq_tools/bin:$PATH
ROOT=/data/project/chureeratp/RNA/HEART
cd "$ROOT"

# --------------------------------------------------------------------------
# Read length. Every path below derives from SUF, so a run is consistent at ONE
# read length or it does not run. Empty = the original 50 bp behaviour, which
# reproduces the 12 Sep output; _rl125 = the re-quantified tables.
#
# 110 renamed the rMATS directory on the M= line and left the gene list and the
# --out-prefix alone, which made this file a three-way mix whose 125 bp run
# would have overwritten the 50 bp comparators.
SUF="${TBL_SUFFIX:-}"
RMATS_DIR="${RMATS_DIR:-rmats_magnet${SUF}}"
# --------------------------------------------------------------------------
CANDIDATES="QKI,RBM5,ZRANB2,HNRNPF,HNRNPA0,MBNL1,CELF1,SRSF1,PTBP1,GAPDH"
MAT=magnet_expr.log2cpm.tsv.gz
PER=magnet_expr.per_sample.tsv
TGT="${TGT:-compadj${SUF}_multi_event_type.txt}"      # the 134 multi-event-type signature genes (count verified 5 Oct 2026; an earlier comment said 38)

for f in "$MAT" "$PER" magnet_covariates.tsv dcm_bams.txt nf_bams.txt; do
  [ -s "$f" ] || { echo "missing $f"; exit 1; }
done

echo "############################################################"
echo "#####  A) augment the per-sample expression table"
echo "############################################################"
python3 - "$MAT" "$PER" "$CANDIDATES" <<'PY'
import sys
import pandas as pd

mat, per_f, cands = sys.argv[1], sys.argv[2], sys.argv[3].split(",")
per = pd.read_csv(per_f, sep="\t")
have = set(per.columns)
need = [g for g in cands if g and g not in have]
print(f"[info] per-sample table already has: "
      f"{sorted(g for g in cands if g in have)}")
if not need:
    print("[info] nothing to add")
    raise SystemExit(0)
print(f"[info] pulling from {mat}: {need}")

m = pd.read_csv(mat, sep="\t", index_col=0)
absent = [g for g in need if g not in m.index]
if absent:
    # A gene can be missing for two very different reasons and they must not be
    # conflated: filtered out by --min-cpm (genuinely low in myocardium, so a null
    # would be uninterpretable), or a symbol mismatch. Either way, stop.
    print(f"\n*** FATAL: {absent} are not in the filtered log2CPM matrix.")
    print("*** Either they fall below step 02's --min-cpm/--min-frac in heart tissue,")
    print("*** or the symbol differs (check HNRNPF vs HNRNPH1, HNRNPA0 vs HNRNPA1).")
    print("*** Do NOT proceed and report these as null -- an untested gene is not a")
    print("*** negative result. Rebuild step 02 with a lower --min-cpm, or drop them")
    print("*** from CANDIDATES and say in the methods that they were not testable.")
    raise SystemExit(1)

sub = m.loc[need].T
sub.index.name = "sample"
out = per.set_index("sample").join(sub, how="left")
miss = out[need].isna().sum()
if miss.any():
    print(f"[warn] samples with no value after join:\n{miss[miss > 0]}")
out.reset_index().to_csv(per_f, sep="\t", index=False, float_format="%.5f")
print(f"[write] {per_f}  now carries {len(out.columns)} columns")
print(out[["RBM20"] + need].describe().T[["count", "mean", "std", "min", "max"]]
      .to_string(float_format=lambda v: f"{v:.3f}"))
PY
[ $? -eq 0 ] || { echo "step A failed -- see the FATAL note above"; exit 1; }

echo
echo "############################################################"
echo "#####  B) the inverted scan"
echo "############################################################"
if [ ! -s "$TGT" ]; then
  echo "*** $TGT missing. Run 13_gene_lists_for_enrichment.py --tag compadj first."
  echo "*** Falling back to compadj_foreground.txt is NOT advised: a target set of"
  echo "*** 500+ genes leaves too little clean matched background and the"
  echo "*** enrichment comparison becomes meaningless."
  exit 1
fi
NT=$(grep -c . "$TGT")
echo "target set: $NT genes from $TGT"
echo "candidates: RBM20 (auto) + $CANDIDATES"
GENES=$(paste -sd, "$TGT")

# --------------------------------------------------------------------------
# A wrong rMATS directory must be an ERROR, not a skip. The per-class guard
# below only `continue`s, which turns a mistyped directory into a job that
# completes having done nothing and reports no failure.
if [ ! -d "$RMATS_DIR" ]; then
  echo "*** rMATS directory not found: $RMATS_DIR"
  echo "*** TBL_SUFFIX='$SUF'. Set RMATS_DIR explicitly if the directory for"
  echo "*** this read length has another name."
  exit 1
fi

# DON'T-CLOBBER GUARD.
# signature_rbp_{SE,MXE}.{events.tsv,report.txt} from 12 Sep are the 50 bp
# comparators that 58_signature_rbp_margins.py consumed to build
# signature_rbp_margins.tsv. Overwriting them makes that table unreproducible
# from its own inputs, and nothing errors.
for ET in SE MXE; do
  for EXT in events.tsv report.txt; do
    F="signature_rbp_${ET}${SUF}.${EXT}"
    if [ -e "$F" ] && [ -z "${OVERWRITE_OK:-}" ]; then
      echo "*** REFUSING — $F already exists"
      echo "***   modified: $(date -r "$F" '+%Y-%m-%d %H:%M' 2>/dev/null)"
      echo "***   size:     $(stat -c %s "$F" 2>/dev/null) bytes"
      echo "***"
      echo "*** Writing over it would destroy a comparator and nothing would"
      echo "*** error. If you mean to replace it, move it aside or re-run with"
      echo "***   OVERWRITE_OK=1 TBL_SUFFIX='$SUF' sbatch $0"
      exit 1
    fi
  done
done
echo "[guard] no output of signature_rbp_*${SUF} exists — safe to write"
# --------------------------------------------------------------------------

for ET in SE MXE; do
  M=${RMATS_DIR}/${ET}.MATS.JCEC.txt
  [ -s "$M" ] || { echo "missing $M -- skipping $ET"; continue; }
  echo
  echo "----- $ET -----"
  python 03_rbm20_psi_correlation.py \
    --rmats "$M" --event-type $ET \
    --b1 dcm_bams.txt --b2 nf_bams.txt \
    --expr "$PER" \
    --covariates magnet_covariates.tsv \
    --numeric-cov age,RIN,log10_libsize,TIN \
    --target-genes "$GENES" \
    --control-genes "$CANDIDATES" \
    --out-prefix signature_rbp_${ET}${SUF}
done

echo
echo "############################################################"
echo "HOW TO READ IT"
echo "############################################################"
echo "1. FIRST check the '[info] matched background events:' line in each run. The"
echo "   target set is larger than the RBM20 nine, so confirm the matched background"
echo "   is still in the hundreds. If it is tiny, the enrichment test is void."
echo
echo "2. FIRST check GAPDH in section 4. If GAPDH scores as high as any real"
echo "   regulator, the model is still technically inflated and NOTHING in the table"
echo "   means anything. That control is why the original RBM20 analysis survived"
echo "   scrutiny -- honour it again."
echo
echo "3. Read section 4, DCM (within-disease) rows, 'cell' adjustment only."
echo
echo "   A CANDIDATE SCORES  -> the mechanism section becomes positive. Follow up with"
echo "     that gene's own known targets as a second, independent target set, and check"
echo "     it replicates in Colorado."
echo "   NOTHING SCORES      -> write: 'the signature is not explained by RBM20, nor by"
echo "     five further regulators nominated for cardiac splicing, including the two"
echo "     proposed for this signature in the current literature.' That closes the"
echo "     reviewer's question instead of leaving it open, which the RBM20 null alone"
echo "     does not."
echo
echo "4. Effect sizes matter more than p-values here. Median |beta| is PSI change per SD"
echo "   of expression; the RBM20 run gave 0.001-0.018, i.e. negligible. A p-value with"
echo "   a beta of 0.002 is not a mechanism."
date
