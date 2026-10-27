#!/usr/bin/env bash
# 89_make_all_figures.sh — regenerate every manuscript figure and verify it against
#                          the published file. For the code deposit and for any
#                          "can this be reproduced?" question at review.
#
# The figure scripts themselves already exist (43, 39, 37, 38, 40, 41, 47). What has
# never existed is one command that runs them all and PROVES the published PNGs come
# out again. This is that command.
#
# It writes to a scratch directory and NEVER touches a published file.
#
#   bash 89_make_all_figures.sh              # rebuild into figrepro/ and compare
#   bash 89_make_all_figures.sh --publish    # ALSO copy verified rebuilds over the
#                                            # published files (refuses on mismatch)
#
# Verified 9 Oct 2026: F2, F3, F4 and F6 are byte-identical on md5 after rebuild.
#
# ⚠ FIGURE 5 IS BUILT BY 90_build_f5.py, NOT BY 40_figure5_mechanism.py.
#   Script 40 drew all four panels on superseded numbers — the retired nine-gene
#   p-values (0.76 / 0.99) in panel B, the 2-D matching in panel C, the nine-gene
#   backgrounds in panel D, and a different covariate model in panel A. It is kept
#   for provenance and must NOT be used to regenerate the published figure. This
#   driver therefore calls 90 for F5, with the model and matching the Results
#   report; change them here and the figure changes, so they live in one place.

set -uo pipefail

OUT=figrepro
PUBLISH=0
[[ "${1:-}" == "--publish" ]] && PUBLISH=1

command -v python >/dev/null || { echo "no python on PATH"; exit 1; }
python -c "import matplotlib" 2>/dev/null || {
  echo "matplotlib not importable — activate the rnaseq_tools environment first:"
  echo "    conda activate rnaseq_tools"; exit 1; }

mkdir -p "$OUT"

# script : scratch stem : published stem : md5-comparable?
FIGS=(
  "43_figure1_design.py:f1:F1_design:yes"
  "39_figure2_discovery.py:f2:F2_discovery:yes"
  "37_figure3_etiology.py:f3:F3_etiology:yes"
  "38_figure4_replication.py:f4:F4_replication:yes"
  "90_build_f5.py:f5:F5_mechanism:yes"
  "41_figure6_enrichment.py:f6:F6_enrichment:yes"
  "47_graphical_abstract.py:ga:GA_graphical_abstract:yes"
)

echo "============================================================"
echo "REBUILD  -> $OUT/   (nothing published is touched)"
echo "============================================================"
declare -A STATUS
for row in "${FIGS[@]}"; do
  IFS=: read -r script stem pub cmp <<< "$row"
  if [[ ! -f "$script" ]]; then
    printf '  %-28s SCRIPT MISSING\n' "$script"; STATUS[$stem]="missing-script"; continue
  fi
  # per-figure arguments, kept here so the published figure has one definition
  EXTRA=()
  case "$stem" in
    f5) EXTRA=(--model 'DCM/cell' --matching '3-D match (+|disease beta|)'
               --b1 dcm_bams.txt --b2 nf_bams.txt) ;;
  esac
  if python "$script" ${EXTRA[@]+"${EXTRA[@]}"} --out "$OUT/$stem" > "$OUT/$stem.log" 2>&1; then
    printf '  %-28s built\n' "$script"; STATUS[$stem]="built"
  else
    printf '  %-28s FAILED (see %s)\n' "$script" "$OUT/$stem.log"; STATUS[$stem]="failed"
  fi
done

echo
echo "============================================================"
echo "COMPARE against the published PNGs"
echo "============================================================"
printf '  %-22s %-34s %s\n' "figure" "published" "result"
ALLOK=1
for row in "${FIGS[@]}"; do
  IFS=: read -r script stem pub cmp <<< "$row"
  [[ "${STATUS[$stem]:-}" == "built" ]] || {
    printf '  %-22s %-34s %s\n' "$stem" "$pub.png" "skipped (${STATUS[$stem]:-?})"
    ALLOK=0; continue; }
  if [[ ! -f "$pub.png" ]]; then
    printf '  %-22s %-34s %s\n' "$stem" "$pub.png" "no published file"; continue
  fi
  a=$(md5sum "$pub.png"       | cut -d' ' -f1)
  b=$(md5sum "$OUT/$stem.png" | cut -d' ' -f1)
  if [[ "$a" == "$b" ]]; then
    printf '  %-22s %-34s IDENTICAL\n' "$stem" "$pub.png"
    STATUS[$stem]="identical"
  elif [[ "$cmp" == "NO" ]]; then
    printf '  %-22s %-34s differs — EXPECTED (see header)\n' "$stem" "$pub.png"
    STATUS[$stem]="differs-expected"
  else
    printf '  %-22s %-34s *** DIFFERS ***\n' "$stem" "$pub.png"
    STATUS[$stem]="differs"; ALLOK=0
  fi
done

echo
echo "  A byte-identical rebuild proves the figure is REPRODUCIBLE from current"
echo "  inputs. It does NOT prove the numbers drawn on it match the current text —"
echo "  read $OUT/*.log for that. Both checks are needed."

if (( PUBLISH )); then
  echo
  echo "============================================================"
  echo "--publish"
  echo "============================================================"
  n=0
  for row in "${FIGS[@]}"; do
    IFS=: read -r script stem pub cmp <<< "$row"
    case "${STATUS[$stem]:-}" in
      identical) ;;                      # nothing to copy
      differs)
        echo "  REFUSING to publish $stem: it differs and is md5-comparable."
        echo "  Read $OUT/$stem.log, decide whether the change is wanted, and copy by hand."
        ;;
      differs-expected)
        echo "  SKIPPING $stem: expected to differ for a known reason; publish by hand"
        echo "  only after checking its provenance block."
        ;;
      *) echo "  SKIPPING $stem (${STATUS[$stem]:-?})" ;;
    esac
  done
  echo "  nothing copied automatically ($n files) — by design."
fi

echo
(( ALLOK )) && echo "ALL COMPARABLE FIGURES REPRODUCE." \
            || echo "SOME FIGURES DID NOT REPRODUCE — see above."
echo "logs: $OUT/*.log   (each ends with the numbers that figure drew)"
