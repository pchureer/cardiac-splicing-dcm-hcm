#!/usr/bin/env bash
# 15a_setup_hcm_cohorts.sh -- run on the LOGIN NODE (needs network).
#
# Prepares the three independent HCM cohorts. INSPECTION AND LIST-BUILDING ONLY,
# nothing expensive. Verify the output before running 15b.
#
#   PRJNA675090  Brigham   18 HCM + 5 control   (+18 WES runs that MUST be excluded)
#   PRJNA522931  Fuwai     28 HCM + 9 donor
#   PRJNA747806  Stanford  13 HCM + 7 donor
#
# WHY THESE MATTER: the MAGNet DCM-vs-HCM result uses the same 166 NF hearts for
# both contrasts, so they share the whole "explanted failing heart vs donor" axis.
# These three cohorts bring their OWN controls, making the comparison fully
# independent. They are also mostly septal myectomy tissue from patients with
# preserved ejection fraction rather than end-stage explants -- a natural
# experiment on whether the signature is end-stage remodelling.
#
# The script:
#   1. fetches each ENA filereport (if not already present)
#   2. classifies HCM vs control, EXCLUDING WES runs by library_strategy
#   3. checks a BAM exists for every sample
#   4. auto-detects read length and paired/single from the BAMs themselves
#   5. reports the GENCODE annotation version each cohort was aligned with
#   6. writes <cohort>_hcm.ids / <cohort>_ctrl.ids and the bam lists

set -uo pipefail
ROOT=/data/project/chureeratp/RNA/HEART
cd "$ROOT"

FIELDS="run_accession,sample_title,experiment_title,library_strategy,library_layout,read_count"

declare -A PRJ=( [Brigham]=PRJNA675090 [Fuwai]=PRJNA522931 [Stanford]=PRJNA747806 )
declare -A EXP_HCM=( [Brigham]=18 [Fuwai]=28 [Stanford]=13 )
declare -A EXP_CTL=( [Brigham]=5  [Fuwai]=9  [Stanford]=7 )

for COH in Brigham Fuwai Stanford; do
  P=${PRJ[$COH]}
  META=${COH}_filereport.tsv
  echo "############################################################"
  echo "#####  $COH  ($P)"
  echo "############################################################"

  if [ ! -s "$META" ]; then
    echo "fetching $P from ENA ..."
    curl -sS -o "$META" \
      "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${P}&result=read_run&fields=${FIELDS}&format=tsv&limit=0" \
      || { echo "  download failed -- fetch in a browser and save as $META"; continue; }
  fi
  echo "rows: $(( $(wc -l < "$META") - 1 ))"

  python3 - "$COH" "$META" <<'PY'
import csv, sys, re
coh, meta = sys.argv[1], sys.argv[2]
rows = list(csv.DictReader(open(meta), delimiter="\t"))
hdr = rows[0].keys() if rows else []

def field(r, *names):
    for n in names:
        if n in r and r[n]:
            return r[n]
    return ""

hcm, ctl, dropped = [], [], []
for r in rows:
    acc = r["run_accession"]
    strat = field(r, "library_strategy").upper()
    txt = " ".join([field(r, "sample_title"), field(r, "experiment_title")])
    low = txt.lower()

    # exclude anything that is not RNA-seq (Brigham mixes in 18 WES runs)
    if strat and strat not in ("RNA-SEQ",):
        dropped.append((acc, f"library_strategy={strat}"))
        continue
    if "exome" in low:
        dropped.append((acc, "title says exome"))
        continue

    if coh == "Brigham":
        if re.search(r"\bctrl\d+", low):            ctl.append(acc)
        elif re.search(r"\bhcm\d+", low):           hcm.append(acc)
        else: dropped.append((acc, "unclassified"))
    elif coh == "Fuwai":
        if "healthy donor" in low:                  ctl.append(acc)
        elif "hcm" in low:                          hcm.append(acc)
        else: dropped.append((acc, "unclassified"))
    else:  # Stanford: names like M#### (HCM) and D#### (donor)
        m = re.search(r"GSM\d+:\s*([A-Za-z]+)\d+", txt)
        pref = (m.group(1) if m else "").upper()
        if pref.startswith("D"):                    ctl.append(acc)
        elif pref.startswith("M"):                  hcm.append(acc)
        else: dropped.append((acc, f"unclassified (prefix '{pref}')"))

open(f"{coh}_hcm.ids", "w").write("\n".join(hcm) + "\n")
open(f"{coh}_ctrl.ids", "w").write("\n".join(ctl) + "\n")
print(f"  HCM      : {len(hcm)}")
print(f"  control  : {len(ctl)}")
print(f"  excluded : {len(dropped)}")
for a, why in dropped[:6]:
    print(f"      {a}  {why}")
if len(dropped) > 6:
    print(f"      ... and {len(dropped)-6} more")
print("  example titles:")
for r in rows[:2]:
    print(f"      {r['run_accession']}  {field(r,'sample_title')[:60]} | "
          f"{field(r,'experiment_title')[:60]}")
PY

  nh=$(grep -c . ${COH}_hcm.ids 2>/dev/null || echo 0)
  nc=$(grep -c . ${COH}_ctrl.ids 2>/dev/null || echo 0)
  echo "  expected: HCM ${EXP_HCM[$COH]}, control ${EXP_CTL[$COH]}"
  [ "$nh" -eq "${EXP_HCM[$COH]}" ] && [ "$nc" -eq "${EXP_CTL[$COH]}" ] \
    && echo "  COUNTS OK" || echo "  *** COUNT MISMATCH -- check the classification above"

  echo
  echo "  --- BAMs, read length, layout, annotation ---"
  miss=0; first=""
  while read -r s; do
    [ -z "$s" ] && continue
    b=$(ls "$ROOT/$s/"*Aligned.sortedByCoord.out.bam 2>/dev/null | head -1)
    if [ -z "$b" ]; then echo "    no BAM: $s"; miss=$((miss+1));
    elif [ -z "$first" ]; then first="$b"; fi
  done < <(cat ${COH}_hcm.ids ${COH}_ctrl.ids)
  echo "    samples without a BAM: $miss"

  if [ -n "$first" ]; then
    echo "    annotation: $(basename "$first" | grep -o 'gencode\.v[0-9]*')"
    RL=$(samtools view "$first" 2>/dev/null | head -2000 | awk '{print length($10)}' \
         | sort -n | uniq -c | sort -rn | head -1 | awk '{print $2}')
    NP=$(samtools view -c -f 1 "$first" 2>/dev/null | head -1)
    echo "    modal read length: ${RL:-unknown}"
    if [ "${NP:-0}" -gt 0 ]; then echo "    layout: PAIRED  (use -t paired)";
    else echo "    layout: SINGLE  (use -t single)"; fi
    echo "$RL" > ${COH}_readlength.txt
  fi

  # bam lists
  for grp in hcm ctrl; do
    > ${COH}_${grp}_bams.txt
    for s in $(cat ${COH}_${grp}.ids); do
      ls "$ROOT/$s/"*Aligned.sortedByCoord.out.bam 2>/dev/null | head -1
    done | paste -sd, > ${COH}_${grp}_bams.txt
    echo "    ${COH}_${grp}_bams.txt: $(tr ',' '\n' < ${COH}_${grp}_bams.txt | grep -c .)"
  done
  echo
done

echo "############################################################"
echo "If all three show COUNTS OK and 0 missing BAMs, run 15b."
echo
echo "NOTE ON ANNOTATION: Fuwai was aligned with GENCODE v46 while MAGNet used"
echo "v49. That is acceptable because each cohort is internally consistent and"
echo "the cross-cohort comparison is at the dPSI level, where a cohort-level"
echo "constant cancels. It does mean junctions annotated only in v49 are"
echo "recovered with more noise in Fuwai, which ATTENUATES its effect sizes --"
echo "conservative for a 'shared program' conclusion. State it in the methods."
