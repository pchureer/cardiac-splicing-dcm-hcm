#!/usr/bin/env bash
# 59_text_consistency_fixes.sh -- close the three loose ends `58` left, and sweep for
# any other place the superseded phrasing survives.
#
# DRY RUN BY DEFAULT. Nothing is written without --apply, and --apply backs up every
# file it touches to <file>.bak.<timestamp> first. No rm anywhere.
#
#   bash 59_text_consistency_fixes.sh            # show what would change
#   bash 59_text_consistency_fixes.sh --apply    # make the changes
#
# ---------------------------------------------------------------------------
# A CORRECTION FIRST, because it changes what item 1 actually is.
#
# I said the graphical abstract "carries the old +0.03 phrasing" and needed
# regenerating. Reading 47_graphical_abstract.py properly: the only visible text in
# block 3 is "11 candidate regulators", the GAPDH label, "none beats the control" and
# the RBM20 expression line. **The figure never prints 0.03.** The string lives in the
# self-check guard and the file docstring, not on the canvas.
#
# So the figure needs NO visual change and NO regeneration. "none beats the control" is
# in fact better supported now than when it was written, because it has intervals behind
# it instead of a point estimate. What needs fixing is the guard's message, which states
# the paper's claim as "by more than 0.03" -- a future reader would take that as the
# claim rather than as a drift tripwire. Step B rewords it and leaves the 0.03 threshold
# in place doing its real job: catching a silent edit to SCAN.
# ---------------------------------------------------------------------------
#
# WHAT IT DOES
#   A  23_signature_rbp_scan.sbatch -- the stale "38 high-confidence genes" comment,
#      replaced with the count actually in the file
#   B  47_graphical_abstract.py -- reword the guard message and the docstring; verify
#      SCAN against 58's recomputed deltas rather than trusting it
#   C  verify the 130 / 124 gene counts from the events tables and print the Methods
#      sentence to paste
#   D  sweep every script and text file for superseded phrasing and report file:line

set -uo pipefail
APPLY=0
[[ "${1:-}" == "--apply" ]] && APPLY=1
ROOT="${ROOT:-/data/project/chureeratp/RNA/HEART}"
cd "$ROOT" || { echo "[fatal] cannot cd $ROOT"; exit 1; }
STAMP=$(date +%Y%m%d-%H%M%S)
(( APPLY )) && echo "[mode] APPLY -- files will be changed, backups at *.bak.$STAMP" \
            || echo "[mode] DRY RUN -- nothing will be written (use --apply)"
echo "[info] $(pwd)"
echo

backup () {   # $1 = file
  if (( APPLY )); then
    cp -p "$1" "$1.bak.$STAMP" || return 1
    echo "      backed up -> $1.bak.$STAMP"
  fi
}

# ===========================================================================
echo "=============================================================="
echo "A. 23_signature_rbp_scan.sbatch -- the stale gene-count comment"
echo "=============================================================="
F=23_signature_rbp_scan.sbatch
TGTFILE=compadj_multi_event_type.txt
if [[ ! -s $F ]]; then
  echo "  [miss] $F -- skipping"
elif [[ ! -s $TGTFILE ]]; then
  echo "  [miss] $TGTFILE -- cannot establish the true count, skipping"
else
  NTRUE=$(tr -s '[:space:]' '\n' < "$TGTFILE" | grep -c .)
  echo "  $TGTFILE holds $NTRUE gene symbols"
  HIT=$(grep -n "high-confidence multi-event-type genes" "$F" || true)
  if [[ -z $HIT ]]; then
    echo "  [ok] no stale count comment found (already fixed?)"
  else
    echo "  current:  $HIT"
    backup "$F"
    if (( APPLY )); then
      python3 - "$F" "$NTRUE" <<'PY'
import re, sys
f, n = sys.argv[1], sys.argv[2]
t = open(f).read()
t2 = re.sub(r"the \d+ high-confidence multi-event-type genes",
            f"the {n} multi-event-type signature genes "
            f"(count verified 5 Oct 2026; an earlier comment said 38)", t)
open(f, "w").write(t2)
print(f"      rewritten to name {n} genes")
PY
    else
      echo "      would become: the $NTRUE multi-event-type signature genes"
      echo "                    (count verified 5 Oct 2026; an earlier comment said 38)"
    fi
  fi
fi

# ===========================================================================
echo
echo "=============================================================="
echo "B. 47_graphical_abstract.py -- guard wording, and SCAN verified"
echo "=============================================================="
G=47_graphical_abstract.py
if [[ ! -s $G ]]; then
  echo "  [miss] $G -- skipping"
else
  echo "  what the figure actually draws in block 3:"
  grep -nE 'ax\.text\(.*(none beats|candidate regulators|HIGHER in DCM)' "$G" \
    | sed 's/^/    /' || echo "    (none matched -- read the file)"
  if grep -q 'exceeds the control by more than 0.03' "$G"; then
    echo
    echo "  the guard states the paper's claim as a 0.03 threshold. Rewording so it"
    echo "  reads as a drift tripwire, which is what it is."
    backup "$G"
    if (( APPLY )); then
      python3 - "$G" <<'PY'
import sys
f = sys.argv[1]
t = open(f).read()
old = ('        print("  *** a candidate now exceeds the control by more than 0.03: "\n'
       '              + ", ".join(f"{k} {v:+.2f}" for k, v in over.items())\n'
       '              + "\\n      Block 3 says \'none beats the control\' and the strap line "\n'
       '                "agrees with it.\\n      Resolve before use.")')
new = ('        print("  *** SCAN has drifted: a candidate now sits more than 0.03 "\n'
       '              "above the control: "\n'
       '              + ", ".join(f"{k} {v:+.2f}" for k, v in over.items())\n'
       '              + "\\n      0.03 is a tripwire on these point estimates, NOT the "\n'
       '                "paper\'s claim. The claim is that no candidate is "\n'
       '                "\\n      DISTINGUISHABLE from the control -- largest margin "\n'
       '                "+0.033 (95% CI -0.025 to +0.082), see 58_signature_rbp_margins.py."\n'
       '              + "\\n      Block 3 draws \'none beats the control\'. Re-check it "\n'
       '                "against 58 before use.")')
if old not in t:
    print("      [warn] the guard block does not match byte-for-byte -- not touching it.")
    print("      [warn] edit the message by hand; the threshold itself is fine.")
    raise SystemExit(0)
open(f, "w").write(t.replace(old, new))
print("      guard message reworded")
PY
    else
      echo "      would reword the message to call 0.03 a tripwire and cite 58's CI"
    fi
  else
    echo "  [ok] guard message already reworded, or absent"
  fi

  # --- SCAN against 58's recomputed deltas --------------------------------
  echo
  echo "  verifying the figure's SCAN constants against 58's output:"
  if [[ ! -s signature_rbp_margins.tsv ]]; then
    echo "    [miss] signature_rbp_margins.tsv -- run 58_signature_rbp_margins.py first"
  else
    python3 - "$G" <<'PY'
import ast, re, sys
import pandas as pd

src = open(sys.argv[1]).read()
m = re.search(r"^SCAN\s*=\s*(\{.*?\})\s*$", src, re.S | re.M)
if not m:
    print("    [warn] could not find a SCAN = {...} literal; check by hand")
    raise SystemExit(0)
try:
    scan = ast.literal_eval(m.group(1))
except Exception as e:
    print(f"    [warn] SCAN is not a plain literal ({e}); check by hand")
    raise SystemExit(0)

df = pd.read_csv("signature_rbp_margins.tsv", sep="\t")
sub = df[(df.event_type == "SE") & (df.model == "DCM/cell") &
         (df.matching.str.startswith("2-D"))]
if sub.empty:
    print("    [warn] no SE DCM/cell 2-D rows in 58's output; check by hand")
    raise SystemExit(0)
got = dict(zip(sub.gene, sub.delta))
floor = float(sub.floor.iloc[0])
got["GAPDH"] = floor

bad, ok = [], 0
for g, v in scan.items():
    if g not in got:
        bad.append((g, v, None)); continue
    if abs(float(v) - got[g]) > 0.015:
        bad.append((g, v, got[g]))
    else:
        ok += 1
print(f"    {ok} of {len(scan)} SCAN values agree with 58 within 0.015")
for g, v, r in bad:
    if r is None:
        print(f"    >>> {g}: in the figure at {v:+.2f} but NOT scored by 58")
    else:
        print(f"    >>> {g}: figure {v:+.2f}, recomputed {r:+.4f}  -- MISMATCH")
if bad:
    print("    >>> the figure's constants no longer match the analysis. Fix SCAN in")
    print("        47_graphical_abstract.py and regenerate the figure.")
else:
    print("    [ok] figure constants match the analysis; NO regeneration needed")
PY
  fi
fi

# ===========================================================================
echo
echo "=============================================================="
echo "C. the gene counts the Methods must quote, computed not asserted"
echo "=============================================================="
if [[ ! -s $TGTFILE ]]; then
  echo "  [miss] $TGTFILE"
else
  python3 - "$TGTFILE" <<'PY'
import os, sys
import pandas as pd

tgt = sorted({t.strip() for t in open(sys.argv[1]).read().split() if t.strip()})
print(f"  target set: {len(tgt)} genes")
out = {}
for et in ("SE", "MXE"):
    f = f"signature_rbp_{et}.events.tsv"
    if not os.path.exists(f):
        print(f"  [miss] {f}")
        continue
    d = pd.read_csv(f, sep="\t", usecols=lambda c: c in ("gene",))
    present = set(d["gene"].astype(str))
    hit = [g for g in tgt if g in present]
    nev = int(d["gene"].isin(tgt).sum())
    out[et] = (len(hit), nev, [g for g in tgt if g not in present])
    print(f"  {et}: {len(hit)}/{len(tgt)} genes contribute {nev} events; "
          f"absent: {','.join(out[et][2]) or 'none'}")
if len(out) == 2:
    se, mx = out["SE"], out["MXE"]
    print()
    print("  PASTE INTO THE METHODS:")
    print(f'    "The signature gene set comprised {len(tgt)} genes significant across')
    print(f'     multiple event classes after composition adjustment, of which')
    print(f'     {se[0]} contributed testable skipped-exon events ({se[1]:,}) and {mx[0]}')
    print(f'     contributed testable mutually exclusive exons ({mx[1]:,}). Genes with no')
    print(f'     event passing the PSI-variance filter contribute nothing to the test and')
    print(f'     are listed in Supplementary Table S7."')
PY
fi

# ===========================================================================
echo
echo "=============================================================="
echo "D. sweep for superseded phrasing anywhere else"
echo "=============================================================="
echo "  Every hit is a place the old framing may still be stated. Judge each one:"
echo "  a NUMBER that is still correct can stay; a CLAIM phrased as a bare margin"
echo "  should become the margin plus its interval."
echo
# maxdepth 1 -- never recurse into the sample tree from here
for pat in \
  'by more than 0\.03' \
  'more than 0\.03' \
  '38 high-confidence' \
  'nine RBM20' \
  'nine published RBM20' \
  'p *= *0\.76' \
  'p *= *0\.99' \
  '1,443'
do
  echo "  --- /$pat/"
  found=0
  while IFS= read -r line; do
    echo "      $line"; found=1
  done < <(find . -maxdepth 1 -type f \
             \( -name '*.py' -o -name '*.sh' -o -name '*.sbatch' -o -name '*.md' \
                -o -name '*.txt' -o -name '*.tex' \) \
             ! -name '59_text_consistency_fixes.sh' -print0 \
           | xargs -0 grep -nE "$pat" 2>/dev/null | head -20)
  (( found )) || echo "      none"
done

echo
echo "=============================================================="
echo "SUMMARY"
echo "=============================================================="
(( APPLY )) && echo "  Changes applied. Backups: *.bak.$STAMP" \
            || echo "  Nothing written. Re-run with --apply once the diffs above look right."
cat <<'EOT'

  Not handled by this script, because they are prose rather than code:
    * DRAFT_methods.md -- the gene counts from section C, and one sentence saying the
      GAPDH floor is a single-gene estimate so margins are reported with bootstrap
      intervals. That sentence is the whole reason 57 and 58 exist.
    * Supplementary S7 -- the RBP delta tables should now carry the CI and Holm columns
      from rbm20_margin_resolve.tsv and signature_rbp_margins.tsv, not bare deltas.

  The figure: if section B reports the constants match, the graphical abstract is
  correct as it stands and does not need regenerating. Only its guard message was
  misleading, and only to whoever reads the source.
EOT
date
