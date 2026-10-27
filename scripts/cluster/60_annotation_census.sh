#!/usr/bin/env bash
# 60_annotation_census.sh -- settle the GENCODE v46 / v49 split, per cohort, per sample.
#
# THE OPEN QUESTION, from DRAFT_methods.md's notes:
#   "Of the first 400 Log.out files, 306 used GENCODE v49 and 94 used v46 -- far more
#    than Fuwai's 37 samples, and the example log shown was a SAMEA accession, which is
#    Heinig. If those 94 are Heinig (never analysed) plus Fuwai, Table 1 is correct as
#    written. If any analysed cohort other than Fuwai is on v46, Table 1 and the Methods
#    are both wrong, and the tested event set for that cohort is not what the paper
#    claims."
#
# This is the only outstanding item that could put a wrong fact in a table, so it gets
# its own script rather than a rerun of all seven sections of 48.
#
# WHAT THIS DOES THAT SECTION 2 OF 48 DOES NOT
#   1. ACCOUNTS FOR EVERY Log.out, not just the ones reachable from a bam list. The
#      actual question is whether the v46 logs belong to analysed cohorts or to the
#      three retrieved-but-never-analysed datasets, and that needs the orphans counted,
#      not ignored.
#   2. Deduplicates the bam lists. 48 globs `dcm_bams.txt ... *_bams.txt`, so the
#      explicitly named lists match twice and print twice.
#   3. Reports the full GTF PATH, not just the vNN token -- two different v49 files in
#      two directories would read as agreement otherwise.
#   4. Cross-checks Table1.tsv's annotation column, if that file is present.
#   5. Writes annotation_census.tsv (one row per sample) so the claim is auditable.
#
# Reads only Log.out files. Touches no BAM, writes one TSV, deletes nothing.
#
#   bash 60_annotation_census.sh
#   bash 60_annotation_census.sh --reuse       # reuse a previous harvest, skip the scan

set -uo pipefail
ROOT="${ROOT:-/data/project/chureeratp/RNA/HEART}"
cd "$ROOT" || { echo "[fatal] cannot cd $ROOT"; exit 1; }
REUSE=0
[[ "${1:-}" == "--reuse" ]] && REUSE=1
HARVEST=.annotation_harvest.tsv
OUT=annotation_census.tsv
date
echo "[info] $(pwd)"

# --- harvest: sample id -> GTF, from every Log.out ---------------------------
# -maxdepth 3 only. A project script once burned a 2 h wall limit recursing into the
# sample tree; do not widen this.
if (( REUSE )) && [[ -s $HARVEST ]]; then
  echo "[info] reusing $HARVEST ($(( $(wc -l < $HARVEST) - 1 )) rows)"
else
  echo "[info] scanning for Log.out files..."
  mapfile -t LOGS < <(find . -maxdepth 3 -name '*Log.out' -size +1c 2>/dev/null | sort)
  echo "[info] ${#LOGS[@]} Log.out files found (no cap)"
  if (( ${#LOGS[@]} == 0 )); then
    echo "[fatal] no Log.out files. Without them the annotation cannot be established"
    echo "        from the data at all, and Table 1's annotation column is unsupported."
    exit 1
  fi
  printf 'sample\tgtf_path\tgtf_tag\tlog\n' > "$HARVEST"
  for f in "${LOGS[@]}"; do
    # sample id from the path; fall back to the log's own directory name
    id=$(printf '%s\n' "$f" | grep -oE '(SRR|ERR|SAMEA|SAMN)[0-9]+' | head -1)
    [[ -z $id ]] && id=$(basename "$(dirname "$f")")
    p=$(grep -m1 -oE '\-\-sjdbGTFfile[[:space:]]+[^[:space:]]+' "$f" 2>/dev/null \
        | awk '{print $2}')
    [[ -z $p ]] && p="(no --sjdbGTFfile in log)"
    tag=$(printf '%s\n' "$p" | grep -oE 'gencode\.v[0-9]+' | head -1)
    [[ -z $tag ]] && tag="UNKNOWN"
    printf '%s\t%s\t%s\t%s\n' "$id" "$p" "$tag" "$f" >> "$HARVEST"
  done
  echo "[write] $HARVEST"
fi

# --- the census --------------------------------------------------------------
python3 - "$HARVEST" "$OUT" <<'PY'
import glob, os, re, sys
from collections import Counter, defaultdict

harvest, out_path = sys.argv[1], sys.argv[2]

rows = []
with open(harvest) as fh:
    next(fh)
    for ln in fh:
        p = ln.rstrip("\n").split("\t")
        if len(p) >= 4:
            rows.append(dict(sample=p[0], gtf=p[1], tag=p[2], log=p[3]))
# one row per sample; if a sample has two logs with different GTFs that is itself news
by_sample = defaultdict(set)
log_of = {}
for r in rows:
    by_sample[r["sample"]].add((r["tag"], r["gtf"]))
    log_of.setdefault(r["sample"], r["log"])
print(f"\n[info] {len(rows)} logs covering {len(by_sample)} distinct samples")
multi = {s: v for s, v in by_sample.items() if len(v) > 1}
if multi:
    print(f"  >>> {len(multi)} sample(s) have logs with MORE THAN ONE annotation:")
    for s, v in list(multi.items())[:10]:
        print(f"      {s}: {sorted(t for t, _ in v)}")
    print("      These were re-aligned at some point. Decide which run fed rMATS.")

# --- cohort membership from the bam lists -----------------------------------
# dedupe by resolved path: 48's glob matched the named lists twice
lists = sorted({os.path.realpath(p) for p in glob.glob("*_bams*.txt")})
COHORT = [
    (re.compile(r"^(dcm|nf|hcm)_bams", re.I), "MAGNet"),
    (re.compile(r"^colorado", re.I), "Colorado"),
    (re.compile(r"^lbnl", re.I), "LBNL"),
    (re.compile(r"^brigham", re.I), "Brigham"),
    (re.compile(r"^fuwai", re.I), "Fuwai"),
    (re.compile(r"^stanford", re.I), "Stanford"),
]
ID = re.compile(r"(SRR|ERR|SAMEA|SAMN)[0-9]+")

cohort_ids = defaultdict(set)
per_list = []
unmapped_lists = []
for p in lists:
    base = os.path.basename(p)
    coh = next((c for rx, c in COHORT if rx.match(base)), None)
    ids = {m.group(0) for m in ID.finditer(open(p, errors="replace").read())}
    per_list.append((base, coh, len(ids)))
    if coh is None:
        unmapped_lists.append(base)
    else:
        cohort_ids[coh] |= ids

print(f"\n[info] {len(lists)} bam list(s), deduplicated by path")
for base, coh, n in sorted(per_list):
    print(f"    {base:<34} {coh or '*** UNRECOGNISED ***':<10} {n:>4} ids")
if unmapped_lists:
    print("  >>> Lists above marked UNRECOGNISED were not assigned to a cohort, so")
    print("      their samples are NOT checked. Extend COHORT in this script.")

EXPECTED = {"MAGNet": "gencode.v49", "Colorado": "gencode.v49", "LBNL": "gencode.v49",
            "Brigham": "gencode.v49", "Fuwai": "gencode.v46", "Stanford": "gencode.v49"}

print("\n" + "=" * 78)
print("ANNOTATION PER ANALYSED COHORT")
print("=" * 78)
print(f"  {'cohort':<10}{'samples':>8}{'with logs':>11}  annotation(s)")
problems, missing_logs = [], []
for coh in sorted(cohort_ids):
    ids = cohort_ids[coh]
    tags = Counter()
    nolog = 0
    for s in ids:
        if s in by_sample:
            for t, _g in by_sample[s]:
                tags[t] += 1
        else:
            nolog += 1
    shown = ", ".join(f"{t} x{n}" for t, n in tags.most_common()) or "NO LOGS"
    print(f"  {coh:<10}{len(ids):>8}{len(ids)-nolog:>11}  {shown}")
    exp = EXPECTED.get(coh)
    if nolog:
        missing_logs.append((coh, nolog, len(ids)))
    if not tags:
        problems.append((coh, "no logs at all -- annotation unverifiable"))
    elif len(tags) > 1:
        problems.append((coh, f"MIXED annotation: {dict(tags)}"))
    elif exp and next(iter(tags)) != exp:
        problems.append((coh, f"is {next(iter(tags))}, Table 1 says {exp}"))

# --- what the orphan logs are: the actual open question ---------------------
assigned = set().union(*cohort_ids.values()) if cohort_ids else set()
orphans = {s: v for s, v in by_sample.items() if s not in assigned}
print("\n" + "=" * 78)
print("LOGS NOT BELONGING TO ANY ANALYSED COHORT")
print("=" * 78)
print("  These should be Heinig, Amsterdam and Naples -- retrieved, genotyped, never")
print("  analysed. This is the accounting that answers the open question.")
otags = Counter()
pref = Counter()
for s, v in orphans.items():
    for t, _g in v:
        otags[t] += 1
    pref[re.match(r"[A-Z]+", s).group(0) if re.match(r"[A-Z]+", s) else "?"] += 1
print(f"  {len(orphans)} orphan samples; annotations: "
      f"{', '.join(f'{t} x{n}' for t, n in otags.most_common()) or 'none'}")
print(f"  accession prefixes: {', '.join(f'{k} x{v}' for k, v in pref.most_common())}")

v46_orphan = sum(n for t, n in otags.items() if t.endswith("v46"))
v46_cohort = 0
for coh in cohort_ids:
    for s in cohort_ids[coh]:
        for t, _g in by_sample.get(s, ()):
            if t.endswith("v46"):
                v46_cohort += 1
print(f"\n  v46 total: {v46_orphan} in orphans + {v46_cohort} in analysed cohorts")

# the open question, answered rather than described
v46_by_cohort = Counter()
for coh in cohort_ids:
    for s in cohort_ids[coh]:
        for t, _g in by_sample.get(s, ()):
            if t.endswith("v46"):
                v46_by_cohort[coh] += 1
if not v46_by_cohort:
    print("  >>> No analysed cohort uses v46 at all -- including Fuwai, which Table 1")
    print("      says is v46. Check Fuwai's logs before trusting either statement.")
elif set(v46_by_cohort) == {"Fuwai"}:
    print(f"  >>> ANSWERED: v46 inside the analysed data is Fuwai ONLY "
          f"({v46_by_cohort['Fuwai']} samples). The rest of the v46 logs "
          f"({v46_orphan}) belong to datasets that were never analysed.")
    print("      The 94 figure in DRAFT_methods.md's notes was Heinig (plus the other")
    print("      retrieved-but-unanalysed sets) plus Fuwai. Table 1 stands.")
else:
    print(f"  >>> ANSWERED, and badly: v46 appears in {dict(v46_by_cohort)}.")
    print("      Any cohort here other than Fuwai contradicts Table 1 and the Methods.")

# --- distinct GTF paths, in full --------------------------------------------
paths = Counter(g for v in by_sample.values() for _t, g in v)
print("\n" + "=" * 78)
print("DISTINCT --sjdbGTFfile PATHS (full, not just the version token)")
print("=" * 78)
for g, n in paths.most_common():
    print(f"  x{n:<6} {g}")
if len([g for g in paths if "v49" in g]) > 1:
    print("  >>> more than one v49 file path. Same version, different files: confirm")
    print("      they are the same annotation before calling the cohorts comparable.")

# --- cross-check Table1.tsv -------------------------------------------------
print("\n" + "=" * 78)
print("CROSS-CHECK AGAINST Table1.tsv")
print("=" * 78)
if not os.path.exists("Table1.tsv"):
    print("  Table1.tsv not found -- skipping. Check the annotation column by hand.")
else:
    txt = open("Table1.tsv", errors="replace").read()
    for coh in sorted(cohort_ids):
        line = next((l for l in txt.splitlines() if coh.lower() in l.lower()), None)
        if line is None:
            print(f"  {coh:<10} no row found in Table1.tsv")
            continue
        claims = set(re.findall(r"v(?:ersion)?\s*(\d{2})", line)) | \
                 set(re.findall(r"gencode[ .]?v?(\d{2})", line, re.I))
        tags = {t for s in cohort_ids[coh] for t, _g in by_sample.get(s, ())}
        obs = {t.split(".v")[-1] for t in tags}
        verdict = "agrees" if claims and claims == obs else (
            "NO ANNOTATION IN THE ROW" if not claims else
            f">>> DISAGREES: table says v{sorted(claims)}, logs say v{sorted(obs)}")
        print(f"  {coh:<10} table v{sorted(claims) or '?'}  logs v{sorted(obs) or '?'}  {verdict}")

# --- write the audit trail --------------------------------------------------
owner = {}
for coh, ids in cohort_ids.items():
    for s in ids:
        owner[s] = coh
with open(out_path, "w") as fh:
    fh.write("sample\tcohort\tgtf_tag\tgtf_path\tlog\n")
    for s in sorted(by_sample):
        for t, g in sorted(by_sample[s]):
            fh.write(f"{s}\t{owner.get(s, 'NOT_ANALYSED')}\t{t}\t{g}\t"
                     f"{log_of.get(s, '')}\n")
print(f"\n[write] {out_path}")

# --- verdict ----------------------------------------------------------------
print("\n" + "=" * 78)
print("VERDICT")
print("=" * 78)
if missing_logs:
    for coh, n, tot in missing_logs:
        print(f"  [warn] {coh}: {n} of {tot} samples have no Log.out -- their annotation")
        print(f"         is unverified. If that is most of the cohort, say so in the")
        print(f"         Methods rather than claiming an annotation for it.")
if not problems:
    print("  [ok] every analysed cohort is homogeneous and matches Table 1.")
    print("       Table 1 and the Methods stand as written; close the open question in")
    print("       DRAFT_methods.md's notes.")
else:
    for coh, msg in problems:
        print(f"  >>> {coh}: {msg}")
    print()
    print("  If an analysed cohort is on an annotation Table 1 does not claim, three")
    print("  things are wrong at once and all three must be fixed:")
    print("    1. Table 1's annotation column and the Methods sentence.")
    print("    2. That cohort's TESTED EVENT SET -- a different GENCODE release means a")
    print("       different event universe, so its denominators are not what the paper")
    print("       says and its event counts are not directly comparable.")
    print("    3. The event-counting convention paragraph, which already notes that the")
    print("       duplication factor varies with annotation density. If more than one")
    print("       analysed cohort is affected, that sentence understates the problem.")
    print("  A MIXED cohort is worse than a wrong one: it cannot be described by any")
    print("  single annotation, and the honest fix is to report the split or re-run.")
PY
rc=$?
date
exit $rc
