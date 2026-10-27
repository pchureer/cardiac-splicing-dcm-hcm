#!/bin/bash
#SBATCH -p cpu
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -t 02:00:00
#SBATCH --qos=ngs16
#SBATCH -J rbm20_g45
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#
# Re-report the RBM20 target-enrichment test against the FULL published reference set
# of 45 genes (Gregorich et al., Circ Genom Precis Med 2024;17:e004355, Table 2) instead
# of the nine-gene subset the frozen result used.
#
# WHY THIS IS CHEAP: --target-genes only affects the REPORT. It decides which fitted
# events are called "target" and which form the matched background; it does not touch
# the per-event fitting, which is the ~30 min part. So this runs in resume mode
# (--from-events) over the events tables that already exist, and the per-event
# statistics are byte-identical to the published run. The ONLY thing that changes is
# the target/background split. That makes this a cleaner comparison than a refit would
# be, because a refit would also carry fitting noise.
#
# WHAT IT DOES, in order:
#   1. preflight -- the 45-gene file matches the list hard-coded here, and every
#      events table exists
#   2. REPRODUCTION CHECK -- re-report at the NINE and confirm the frozen numbers come
#      back. If the resume path cannot reproduce p = 0.76 / 0.99, nothing downstream is
#      trustworthy and the job stops.
#   3. re-report at the 45
#   4. gene coverage -- how many of the 45 actually contributed testable events
#   5. side-by-side comparison, with the decision rule applied
#
# THE EXPECTATION, RECORDED BEFORE THE RESULT:
#   Still null. Two reasons. The nine-gene test was already biased TOWARD enrichment
#   (GAPDH scored on it), and the 45 adds non-sarcomeric genes -- MECP2, ZNF451, APTX,
#   NFIA, SLC38A10 -- whose events carry less disease effect, which should dilute that
#   bias rather than add signal. A null at 45 is therefore a stronger null than a null
#   at 9, and that is the point of running it.
#
#   THE DECISION RULE, also fixed before the result: the null HOLDS if RBM20's
#   enrichment delta (target median |t| minus matched-background median |t|) does not
#   exceed the GAPDH control's delta in the same model. It FAILS if RBM20 beats GAPDH.
#   An MWU p-value on its own does not decide it -- the test is inflated and the paper's
#   own standard is the housekeeping floor.
#
# NO rm ANYWHERE. Writes only new files under new prefixes; the published
# rbm20_*.report.txt and rbm20_*.events.tsv are never touched.
#
# Submit:  sbatch 56_rbm20_targets45.sbatch
#          FORCE=1 sbatch 56_rbm20_targets45.sbatch     # redo existing outputs

set -uo pipefail
date
export PATH=/data/apps/anaconda3/envs/mamba/envs/rnaseq_tools/bin:$PATH
export PYTHONUNBUFFERED=1
ROOT=/data/project/chureeratp/RNA/HEART
cd "$ROOT" || { echo "[fatal] cannot cd $ROOT"; exit 1; }

FORCE="${FORCE:-0}"
G45FILE=gregorich2024_table2.txt
# --------------------------------------------------------------------------
# Read length. Called RLSUF, NOT SUF: line ~104 uses SUF for the model variant
# ("" | "_techadj") and would overwrite a read-length SUF on its first
# iteration, silently mis-naming every path.
#
# The suffix goes LAST, matching what 126 put into 03 and 03b:
#     03   -> rbm20_SE_rl125.events.tsv
#     03b  -> rbm20_SE_techadj_rl125.events.tsv
# so rbm20_${ET}${SUF}${RLSUF} resolves to both. RUNS keeps the bare tag,
# because run_report derives et="${tag%%_*}" from it.
RLSUF="${TBL_SUFFIX:-}"
# --------------------------------------------------------------------------
SCRIPT=03_rbm20_psi_correlation.py

# The 45, hard-coded so a truncated or edited file is caught rather than silently used.
G45_EXPECT="APTX,CACNA1C,CAMK2D,CAMK2G,DAB1,DNM3,DST,DTNA,ENAH,FHOD3,FNBP1,GIT2,IMMT,KALRN,KCNIP2,LDB3,LMO7,LRRFIP1,MECP2,MLIP,MTMR1,MYH7,MYOM1,NEXN,NFIA,NPRL3,NTRK3,OBSCN,PDLIM3,PDLIM5,PLEKHA5,RALGPS1,RTN4,RYR2,SEMA6D,SH3KBP1,SLC38A10,SORBS1,SPEN,TNNT2,TPM1,TRDN,TTN,UBE2F,ZNF451"
NINE="TTN,CAMK2D,LDB3,LMO7,SORBS1,PDLIM3,DTNA,MLIP,TNNT2"

echo "=============================================================="
echo "1. PREFLIGHT"
echo "=============================================================="
[ -s "$SCRIPT" ] || { echo "[fatal] $SCRIPT not found in $ROOT"; exit 1; }

if [ -s "$G45FILE" ]; then
  G45=$(tr -s '[:space:],' '\n' < "$G45FILE" | grep -E '^[A-Z][A-Z0-9]{1,9}$' | sort -u | paste -sd,)
  N45=$(echo "$G45" | tr ',' '\n' | grep -c .)
  EXP_SORTED=$(echo "$G45_EXPECT" | tr ',' '\n' | sort | paste -sd,)
  if [ "$G45" != "$EXP_SORTED" ]; then
    echo "[fatal] $G45FILE does not match the 45 hard-coded in this script."
    echo "        file has $N45 symbols. Differences:"
    diff <(echo "$G45" | tr ',' '\n') <(echo "$EXP_SORTED" | tr ',' '\n') | sed 's/^/        /'
    echo "        Fix the file, or edit G45_EXPECT here if the reference list changed"
    echo "        on purpose -- but then say so in the Methods."
    exit 1
  fi
  echo "[ok] $G45FILE matches the hard-coded list: $N45 genes"
else
  echo "[warn] $G45FILE not found -- using the list hard-coded in this script"
  G45=$(echo "$G45_EXPECT" | tr ',' '\n' | sort | paste -sd,)
  N45=45
fi

# the nine must be a subset of the 45, or the comparison is not nested
MISSING=$(comm -23 <(echo "$NINE" | tr ',' '\n' | sort) <(echo "$G45" | tr ',' '\n') | paste -sd,)
if [ -n "$MISSING" ]; then
  echo "[fatal] the nine are not a subset of the 45; missing: $MISSING"
  exit 1
fi
echo "[ok] the nine are a clean subset of the 45"

# events tables: base (03) and technical-covariate-adjusted (03b). The published
# numbers come from the techadj run; both are re-reported so neither can be confused
# for the other later.
RUNS=""
for ET in SE MXE; do
  for SUF in "" "_techadj"; do
    EV="rbm20_${ET}${SUF}${RLSUF}.events.tsv"
    if [ -s "$EV" ]; then
      n=$(($(wc -l < "$EV") - 1))
      echo "[ok] $EV  ($n events)"
      RUNS="$RUNS ${ET}${SUF}"
    else
      echo "[miss] $EV -- will be skipped"
    fi
  done
done
if [ -z "$RUNS" ]; then
  echo "[fatal] no events tables found. Expected rbm20_{SE,MXE}[_techadj].events.tsv"
  echo "        in $ROOT. Run 03_rbm20_psi.sbatch / 03b first."
  ls -1 rbm20_*.events.tsv 2>/dev/null | sed 's/^/        found: /'
  exit 1
fi
echo "[info] runs in scope:$RUNS"

# the covariate flags must match how each events table was fitted, or section 0's
# RBM20-expression models are refitted under the wrong design
numcov_for () {
  case "$1" in
    *_techadj) echo "age,RIN,log10_libsize,TIN" ;;
    *)         echo "age,RIN" ;;
  esac
}

run_report () {   # $1=run tag  $2=gene list  $3=out suffix
  local tag="$1" genes="$2" suf="$3"
  local et="${tag%%_*}"
  local ev="rbm20_${tag}${RLSUF}.events.tsv"
  local pref="rbm20_${tag}${suf}${RLSUF}"
  if [ -s "${pref}.report.txt" ] && [ "$FORCE" != "1" ]; then
    echo "[skip] ${pref}.report.txt exists (FORCE=1 to redo)"
    return 0
  fi
  python "$SCRIPT" \
    --from-events "$ev" --event-type "$et" \
    --b1 dcm_bams.txt --b2 nf_bams.txt \
    --expr magnet_expr.per_sample.tsv \
    --covariates magnet_covariates.tsv \
    --numeric-cov "$(numcov_for "$tag")" \
    --target-genes "$genes" \
    --out-prefix "$pref" \
    > "${pref}.stdout" 2>&1
  local rc=$?
  if [ $rc -ne 0 ] || [ ! -s "${pref}.report.txt" ]; then
    echo "[fail] $pref  (exit $rc) -- tail of its log:"
    tail -20 "${pref}.stdout" | sed 's/^/        /'
    return 1
  fi
  echo "[ok] ${pref}.report.txt"
  return 0
}

echo
echo "=============================================================="
echo "2. REPRODUCTION CHECK -- re-report at the NINE from the same tables"
echo "=============================================================="
echo "If the resume path reproduces the frozen numbers, the 45-gene numbers that"
echo "follow are comparable to them. If it does not, stop and find out why."
for tag in $RUNS; do run_report "$tag" "$NINE" "_repro9" || exit 1; done

echo
echo "=============================================================="
echo "3. RE-REPORT AT THE 45"
echo "=============================================================="
for tag in $RUNS; do run_report "$tag" "$G45" "_g45" || exit 1; done

echo
echo "=============================================================="
echo "4. GENE COVERAGE -- how many of the 45 have testable events"
echo "=============================================================="
python3 - $RUNS <<PY
import os, sys, pandas as pd
G45 = "$G45".split(",")
NINE = "$NINE".split(",")
for tag in sys.argv[1:]:
    ev = f"rbm20_{tag}$RLSUF.events.tsv"
    if not os.path.exists(ev):
        continue
    df = pd.read_csv(ev, sep="\t", usecols=lambda c: c in ("gene",))
    present = set(df["gene"].astype(str))
    h45 = [g for g in G45 if g in present]
    h9  = [g for g in NINE if g in present]
    n45 = int(df["gene"].isin(G45).sum())
    n9  = int(df["gene"].isin(NINE).sum())
    print(f"\n{tag}:  {len(df)} events tested")
    print(f"  the 45: {len(h45)}/45 genes contribute {n45} target events")
    print(f"  the  9: {len(h9)}/9 genes contribute {n9} target events")
    miss = [g for g in G45 if g not in present]
    if miss:
        print(f"  absent from the tested set ({len(miss)}): {','.join(miss)}")
    print(f"  >>> the test is really a {len(h45)}-gene test, not a 45-gene test."
          f" The Methods must say {len(h45)}, not 45.")
PY

echo
echo "=============================================================="
echo "5. SIDE BY SIDE, AND THE DECISION"
echo "=============================================================="
RLSUF="$RLSUF" python3 - $RUNS <<'PY'
import os, re, sys

# The heredoc above is QUOTED, so bash substitutes nothing into this
# block. RLSUF arrives through the environment; a literal $RLSUF here
# would be looked for as part of the filename.
RLSUF = os.environ.get("RLSUF", "")

# section 1 lines look like:
#   --- 1. DCM model, cell: |t| for RBM20 ---
#     target median |t| = 0.92   matched-bg median |t| = 0.96   MWU p = 0.76
HDR = re.compile(r"---\s*1\.\s*(\w+)\s*model,\s*(\w+):")
VAL = re.compile(r"target median \|t\|\s*=\s*([-\d.eE+]+).*?matched-bg median \|t\|\s*=\s*"
                 r"([-\d.eE+]+).*?MWU p\s*=\s*([-\d.eEnan+]+)")
# section 4 is a fixed-width table, not prose:
#   gene      model adj      median|t| tgt  median|t| bg     MWU p
#   GAPDH     DCM   cell              0.92          0.96      0.76
CTRL = re.compile(r"^\s+([A-Z][A-Z0-9]{1,9})\s+(DCM|ALL)\s+(nocell|cell)\s+"
                  r"(-?[\d.]+|nan)\s+(-?[\d.]+|nan)\s+(-?[\d.eE+]+|nan)\s*$")

def parse(path):
    """-> {(model, celltag): (target, bg, p)}, {(gene, model, celltag): (target, bg, p)}"""
    if not os.path.exists(path):
        return None, None
    sec1, ctrl, cur, in4 = {}, {}, None, False
    for ln in open(path, errors="replace"):
        h = HDR.search(ln)
        if h:
            cur = (h.group(1), h.group(2)); in4 = False; continue
        if cur:
            v = VAL.search(ln)
            if v:
                sec1[cur] = tuple(float(x) if x != "nan" else float("nan")
                                  for x in v.groups())
                cur = None
                continue
        if re.search(r"---\s*4\.", ln):
            in4 = True; continue
        if in4:
            if ln.strip().startswith("---"):
                in4 = False; continue
            c = CTRL.match(ln.rstrip("\n"))
            if c:
                g, model, ctag = c.group(1).upper(), c.group(2), c.group(3)
                num = lambda s: float("nan") if s == "nan" else float(s)
                try:
                    ctrl[(g, model, ctag)] = (num(c.group(4)), num(c.group(5)),
                                              num(c.group(6)))
                except ValueError:
                    pass
    return sec1, ctrl

FROZEN = {("SE_techadj", ("DCM", "cell")): 0.76,
          ("MXE_techadj", ("DCM", "cell")): 0.99}

verdicts = []
for tag in sys.argv[1:]:
    s9, c9   = parse(f"rbm20_{tag}_repro9{RLSUF}.report.txt")
    s45, c45 = parse(f"rbm20_{tag}_g45{RLSUF}.report.txt")
    if not s9 or not s45:
        print(f"\n{tag}: one of the two reports could not be parsed -- "
              f"read rbm20_{tag}_{{repro9,g45}}{RLSUF}.report.txt by hand")
        continue
    print(f"\n### {tag} " + "#" * (58 - len(tag)))
    print(f"{'model':<16}{'nine: tgt/bg (p)':<30}{'45: tgt/bg (p)':<30}{'delta nine':>11}{'delta 45':>10}")
    for key in sorted(set(s9) | set(s45)):
        a = s9.get(key); b = s45.get(key)
        f = lambda v: "-" if v is None else f"{v[0]:.2f}/{v[1]:.2f} (p={v[2]:.3g})"
        d = lambda v: "-" if v is None else f"{v[0]-v[1]:+.3f}"
        print(f"{key[0]+'/'+key[1]:<16}{f(a):<30}{f(b):<30}{d(a):>11}{d(b):>10}")

    # reproduction check against the frozen p-values
    fz = FROZEN.get((tag, ("DCM", "cell")))
    if fz is not None and ("DCM", "cell") in s9:
        got = s9[("DCM", "cell")][2]
        ok = abs(got - fz) <= 0.02
        print(f"\n  reproduction: DCM/cell MWU p at nine = {got:.3g}, frozen = {fz}  "
              f"-> {'MATCHES' if ok else '*** DOES NOT MATCH ***'}")
        if not ok:
            print("  *** The resume path does not reproduce the published number.")
            print("  *** Do not use the 45-gene result until this is explained.")
            verdicts.append((tag, "REPRODUCTION FAILED"))
            continue

    # the decision rule: RBM20 vs the GAPDH floor, same model, 45-gene target set
    key = ("DCM", "cell")
    if key not in s45:
        print("  no DCM/cell model in the 45-gene report -- cannot apply the rule")
        continue
    rbm_delta = s45[key][0] - s45[key][1]
    gap = (c45 or {}).get(("GAPDH", "DCM", "cell"))
    if gap is None:
        print(f"\n  RBM20 delta at 45 = {rbm_delta:+.3f}; GAPDH floor not parsed from"
              f" section 4 -- read rbm20_{tag}_g45{RLSUF}.report.txt section 4 and compare by hand")
        verdicts.append((tag, f"delta {rbm_delta:+.3f}, floor unread"))
        continue
    gap_delta = gap[0] - gap[1]
    print(f"\n  RBM20 delta = {rbm_delta:+.3f}   GAPDH floor = {gap_delta:+.3f}   "
          f"margin = {rbm_delta - gap_delta:+.3f}")
    if rbm_delta <= gap_delta:
        v = "NULL HOLDS at 45 -- RBM20 does not beat the housekeeping floor"
    else:
        v = (f"NULL FAILS at 45 -- RBM20 exceeds GAPDH by "
             f"{rbm_delta - gap_delta:+.3f}; the mechanism claim has to change")
    print(f"  >>> {v}")
    verdicts.append((tag, v))

print("\n" + "=" * 62)
print("VERDICT")
print("=" * 62)
for tag, v in verdicts:
    print(f"  {tag:<16} {v}")
if not verdicts:
    print("  nothing to report -- check the per-run logs above")
print()
print("If the null holds at 45: the Methods and Discussion say the test used the")
print("published reference set of 45 -- but quote the gene count from STEP 4 of this")
print("job, not 45, since genes with no testable events contribute nothing. The")
print("nine-of-45 selection question disappears, and the comparison with Radke et al.'s")
print("17/45 and 25/45 is on the same gene set.")
print("If it fails: the 45-gene result is the one to report, the nine-gene null was a")
print("subset artefact, and this needs a full refit with --n-perm 200, since the")
print("permutation null is NOT available in resume mode.")
PY

echo
echo "--- files written ---"
ls -lh rbm20_*_repro9.report.txt rbm20_*_g45.report.txt 2>/dev/null
echo
echo "[note] the published rbm20_*.report.txt and rbm20_*.events.tsv were not touched."
date
