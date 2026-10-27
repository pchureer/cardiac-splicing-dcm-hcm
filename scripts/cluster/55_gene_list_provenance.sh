#!/usr/bin/env bash
# 55_gene_list_provenance.sh -- audit every hard-coded gene list in the pipeline and
# say which ones have a citable source on record. This is the check behind reference
# budget slots M16 (marker panels) and D10-D14 (regulator target lists).
#
# It answers one question per list: if a reviewer asks "where did these genes come
# from", is the answer in the repository, or only in someone's memory?
#
# Usage:
#   bash 55_gene_list_provenance.sh [SCRIPTDIR]
#   bash 55_gene_list_provenance.sh --check-against published_list.txt [SCRIPTDIR]
#
# --check-against takes a plain text file of gene symbols (one per line, or
# whitespace/comma separated -- it does not care) transcribed from a published table,
# and tests each audited list for subset/extras against it. That converts "these look
# like Maatz's targets" into a verified statement or a disproved one.
#
# Writes gene_list_provenance.tsv. Reads only; creates nothing else; deletes nothing.

set -uo pipefail

REF=""
if [[ "${1:-}" == "--check-against" ]]; then
  REF="${2:?--check-against needs a file}"; shift 2
fi
DIR="${1:-.}"
SELF="$(basename "$0")"
OUT="gene_list_provenance.tsv"

cd "$DIR" || { echo "[fatal] no such directory: $DIR" >&2; exit 1; }
echo "[info] auditing $(pwd)"
[[ -n "$REF" ]] && echo "[info] reference list: $REF"

# maxdepth 1 on purpose. An earlier script in this project recursed into 716 sample
# directories and burned a 2-hour wall limit; never scan the data tree from here.
mapfile -t SCRIPTS < <(find . -maxdepth 1 -type f \
  \( -name '*.py' -o -name '*.sh' -o -name '*.sbatch' -o -name '*.R' \) \
  ! -name "$SELF" | sort)
echo "[info] ${#SCRIPTS[@]} scripts in scope (maxdepth 1, excluding $SELF)"
if (( ${#SCRIPTS[@]} == 0 )); then
  echo "[fatal] no scripts found -- wrong directory?" >&2; exit 1
fi

# ---------------------------------------------------------------------------
# The lists the manuscript has to cite, named explicitly rather than discovered:
# a list that has gone missing must be reported as missing, and a pattern scan
# cannot report an absence.
# ---------------------------------------------------------------------------
EXPECT="$(mktemp)"
cat > "$EXPECT" <<'EOF'
RBM20_targets	published RBM20 target genes	D10-D14	TTN,CAMK2D,LDB3,LMO7,SORBS1,PDLIM3,DTNA,MLIP,TNNT2
candidate_RBPs	candidate regulators + GAPDH control	D10-D14	MBNL1,CELF1,SRSF1,PTBP1,QKI,GAPDH,RBM5,ZRANB2,HNRNPF,HNRNPA0,RBM20
CM_markers	cardiomyocyte composition score panel	M16	TNNI3,MYL2,MYL3,ACTN2,TCAP,CSRP3,MYBPC3,RYR2,PLN
FB_markers	fibroblast composition score panel	M16	COL1A1,COL1A2,COL3A1,POSTN,FN1,DCN,LUM
EOF

python3 - "$REF" "$OUT" "$EXPECT" "${SCRIPTS[@]}" <<'PY'
import re, sys, os
from collections import defaultdict

ref_path, out_path, expect_path = sys.argv[1], sys.argv[2], sys.argv[3]
scripts = sys.argv[4:]

# A gene symbol as it appears in these scripts: 2-10 chars, starts with a letter,
# upper case, may carry digits. Strict on purpose -- better to miss an exotic symbol
# than to report every SLURM variable as a gene.
SYM = re.compile(r"^[A-Z][A-Z0-9]{1,9}$")

# Tokens that pass SYM but are not genes. Two groups, both needed: pipeline jargon,
# and the common English words that appear in upper-case banner strings ("HOW TO READ
# IT"). Without the second group every print() banner in the project is a false
# positive. Nothing here is a gene symbol used in this pipeline.
STOP = set("""
TRUE FALSE NONE NULL NA NAN INF AND OR NOT IF FI DO DONE THEN ELSE
SE MXE A3SS A5SS RI JC JCEC PSI DPSI FDR TSV CSV TXT PDF PNG SVG
BAM SAM VCF GTF GFF BED FASTQ PATH HOME USER DIR OUT IN ID IDS
DCM HCM ICM NF NFD LV RV RNA DNA CPM TPM TIN RIN SD CM FB GB MB KB TB
CPU RAM QOS JOB ERR LOG MD5 OK YES NO R1 R2 P1 P2 B1 B2 N1 N2
X Y Z T P Q N M A B C D E F G H I J K L O R S U V W
KEGG REAC WP GO HP CORUM TF MIRNA HPA
THE THIS THAT THESE THOSE IS ARE WAS WERE BE BEEN IT ITS OF TO FOR FROM
WITH WITHOUT AT ON BY AS AN A HOW WHY WHAT WHEN WHERE WHICH WHO
READ WRITE RUN RUNS RAN USE USED USING MAKE MADE SET GET PUT ADD ALL ANY
EVERY EACH ROW ROWS COL COLS LINE LINES FILE FILES NOTE NOTES STEP STEPS
POOL POOLED EXCLUDED EXCLUDE INCLUDE ONLY BOTH SAME DIFF SPLIT JOIN
CONTROL CONTROLS NORMAL DONOR CASE CASES GROUP GROUPS TOTAL COUNT COUNTS
NUMBER NUMBERS VALUE VALUES RESULT RESULTS CHECK CHECKS TEST TESTS
PAPER TABLE FIGURE PANEL LEGEND TEXT MAIN SUPP
MATTERS COMPARISON COMPARE REPORT REPORTED FROZEN FILL BLUE ORANGE PALE
GREY GRAY INK MUTED GRID SURF BOLD SCAN SCANS STILL NEXT PREV FIRST LAST
WARN INFO FATAL ERROR DEBUG PASS FAIL SKIP KEEP DROP OLD NEW RAW ADJ
REPL ET MEAN MEDIAN MIN MAX SUM AVG STD VAR CI SEM
""".split())

# A source comment is one that actually names a provenance: an author-year, a DOI, a
# PMID, a URL, a journal, a table reference, or an explicit statement that the list
# was assembled by hand.
SRC = re.compile(
    r"(doi[:\s]|10\.\d{4}/|pmid|pmc\d|https?://|et\s+al\.?|"
    r"\b(19|20)\d{2}\b[^\n]{0,60}\b(journal|nature|circ|cell|jci|science|heliyon|"
    r"genome|lancet|eur|am j|j clin|j mol|j transl|bioinform)|"
    r"table\s*\d|tbl\s*\d|supplement|suppl\.|"
    r"hand[- ]?(picked|chosen|curated|selected|assembled)|chosen by (hand|eye)|"
    r"\bno source\b|\bsource unknown\b|provenance)",
    re.I,
)

def genes_in(text):
    toks = [t.strip().strip("'\"`") for t in re.split(r"[,\s]+", text)]
    return [t for t in toks if SYM.match(t) and t not in STOP]

def line_of(text, pos):
    return text.count("\n", 0, pos) + 1

def comment_near(lines, lineno, width=10):
    """Contiguous comment lines immediately above, plus a trailing comment on the line
    itself. A source note is written next to the list or it does not count."""
    out, j, seen = [], lineno - 2, 0
    while j >= 0 and seen < width:
        s = lines[j].strip()
        if s.startswith(("#", "//", '"""', "'''", "*")):
            out.append(s.lstrip("#/*\"' ")); j -= 1; seen += 1
        elif s == "":
            j -= 1; seen += 1                 # blanks do not break the block
        else:
            break
    tail = lines[lineno - 1] if lineno - 1 < len(lines) else ""
    for mark in ("#", "//"):
        if mark in tail:
            out.insert(0, tail.split(mark, 1)[1].strip()); break
    return " | ".join(reversed(out))

# ---- find candidate lists over the WHOLE file, not line by line -------------
# Line-by-line scanning splits a list written across two lines into two half-lists,
# which then look like two different gene sets. Both patterns below are multi-line.
Q_STR  = re.compile(r"""(['"])([A-Z0-9_,\s]{8,600}?)\1""")              # "A,B,C"
BRACKET = re.compile(r"""[\[\(\{]([^\[\]\(\){}]{8,1200})[\]\)\}]""", re.S)  # [ "A", "B" ]

found = []   # dict per find
for path in scripts:
    try:
        text = open(path, encoding="utf-8", errors="replace").read()
    except OSError as e:
        print(f"[warn] cannot read {path}: {e}"); continue
    lines = text.splitlines()
    head  = " ".join(lines[:60])
    raw = []
    for m in Q_STR.finditer(text):
        raw.append((m.start(), genes_in(m.group(2))))
    for m in BRACKET.finditer(text):
        raw.append((m.start(), genes_in(m.group(1))))
    # keep sets of >=3; then drop any find whose genes are a subset of another find
    # in the same file (the bracket match usually subsumes its quoted members)
    cands = [(p, g) for p, g in raw if len(set(g)) >= 3]
    keep = []
    for p, g in cands:
        gs = set(g)
        if any(gs < set(g2) for _p2, g2 in cands):
            continue
        keep.append((p, g))
    seen_sets = set()
    for p, g in sorted(keep):
        key = frozenset(g)
        if key in seen_sets:
            continue
        seen_sets.add(key)
        ln = line_of(text, p)
        # SORTED, not source order. The same set written in a different order is the
        # same list, and keying on source order reports it as a second membership --
        # which sent a reader chasing a difference that did not exist.
        found.append(dict(script=os.path.basename(path), line=ln,
                          genes=sorted(set(g)),
                          comment=comment_near(lines, ln), head=head))

# ---- label each find against the lists the manuscript must cite -------------
expect = []
for ln in open(expect_path):
    if not ln.strip(): continue
    name, desc, slot, genes = ln.rstrip("\n").split("\t")
    expect.append((name, desc, slot, set(genes.split(","))))

VARNAME = re.compile(r"([A-Za-z_][A-Za-z0-9_]*)\s*=\s*$|([A-Za-z_][A-Za-z0-9_]*)\s*=")

def label(genes, script, lineno):
    gs = set(genes)
    best, score = None, 0.0
    for name, _d, _s, known in expect:
        j = len(gs & known) / max(len(gs | known), 1)
        if j > score:
            best, score = name, j
    if score >= 0.45:
        return best, score
    return "unclassified", score

rows = []
for f in found:
    name, score = label(f["genes"], f["script"], f["line"])
    rows.append(dict(
        list=name, script=f["script"], line=f["line"], n=len(f["genes"]),
        overlap=f"{score:.2f}", genes=",".join(f["genes"]),
        source=("YES" if SRC.search(f["comment"] or "") else
                ("file-docstring" if SRC.search(f["head"]) else "NONE")),
        comment=re.sub(r"\s+", " ", (f["comment"] or ""))[:300],
    ))

with open(out_path, "w") as fh:
    cols = ["list","script","line","n","overlap","source","genes","comment"]
    fh.write("\t".join(cols) + "\n")
    for r in sorted(rows, key=lambda r: (r["list"], r["script"], r["line"])):
        fh.write("\t".join(str(r[c]) for c in cols) + "\n")

# ---- report ----------------------------------------------------------------
named = [r for r in rows if r["list"] != "unclassified"]
other = [r for r in rows if r["list"] == "unclassified"]

print()
print("=" * 78)
print("SECTION 1 -- the lists the manuscript must cite")
print("=" * 78)
for name, desc, slot, known in expect:
    hits = [r for r in rows if r["list"] == name]
    print(f"\n  {name}  (budget slot {slot})  -- {desc}")
    if not hits:
        print("      >>> NOT FOUND IN ANY SCRIPT.")
        print("          Either it lives in a data file rather than a script, or it was")
        print("          applied by hand. Either way the Methods must say which.")
        continue
    memberships = defaultdict(list)
    for r in hits:
        memberships[r["genes"]].append(f"{r['script']}:{r['line']}")
    srcs = {r["source"] for r in hits}
    print(f"      found in {len(hits)} place(s), {len(memberships)} distinct membership(s)")
    for mv, where in sorted(memberships.items(), key=lambda kv: -len(kv[1])):
        gl = mv.split(",")
        print(f"        n={len(gl):<3} {mv}")
        print(f"             in: {', '.join(where)}")
    if len(memberships) > 1:
        # Nested subsets of one parent are a different situation from genuinely
        # divergent lists: the provenance question collapses to the parent, and the
        # subsets need a selection rule rather than sources of their own.
        sets = {mv: set(mv.split(",")) for mv in memberships}
        parents = [mv for mv, s in sets.items()
                   if all(s2 <= s for s2 in sets.values())]
        if parents:
            par = parents[0]
            kids = [mv for mv in sets if mv != par]
            print(f"      >>> NESTED, not divergent: every membership is a subset of the")
            print(f"          n={len(sets[par].copy())} list in "
                  f"{', '.join(memberships[par])}.")
            for mv in kids:
                extra = sorted(sets[par] - sets[mv])
                print(f"            n={len(sets[mv])} = parent minus {','.join(extra)}")
            print("          So there is ONE list to source, and the subsets need a")
            print("          stated selection rule, not citations of their own.")
        else:
            print("      >>> INCONSISTENT MEMBERSHIP, and not nested. The paper can")
            print("          describe only one of these. Decide which produced the")
            print("          reported result, and check whether any reported number")
            print("          came from another.")
            allg = set().union(*sets.values())
            common = set.intersection(*sets.values())
            print(f"          union {len(allg)}, shared by all {len(common)}: "
                  f"{','.join(sorted(common))}")
    if srcs == {"YES"}:
        print("      verdict: SOURCE ON RECORD")
    elif "YES" in srcs:
        print("      verdict: SOURCE ON RECORD IN SOME PLACES ONLY -- copies disagree")
    elif srcs == {"file-docstring"}:
        print("      verdict: NO SOURCE AT THE LIST; the file docstring cites something,")
        print("               which may or may not be the list's source. Read it.")
    else:
        print("      verdict: >>> NO SOURCE ANYWHERE. This is the hole.")
    for r in hits:
        if r["comment"]:
            print(f"      comment @ {r['script']}:{r['line']}: {r['comment'][:180]}")

print()
print("=" * 78)
print(f"SECTION 2 -- {len(other)} other gene-like lists, for completeness")
print("=" * 78)
print("  Scored against the four expected lists and matched none. Most will be")
print("  unrelated; scan for anything that is actually a fifth gene list.")
for r in sorted(other, key=lambda r: (r["script"], r["line"])):
    print(f"    {r['script']}:{r['line']}  n={r['n']}  {r['genes'][:160]}")

# ---- subset test against a published list ----------------------------------
if ref_path:
    print()
    print("=" * 78)
    print(f"SECTION 3 -- subset test against {os.path.basename(ref_path)}")
    print("=" * 78)
    if not os.path.exists(ref_path):
        print(f"  [fatal] {ref_path} does not exist"); sys.exit(2)
    ref = set(genes_in(open(ref_path, encoding="utf-8", errors="replace").read()))
    print(f"  reference list: {len(ref)} symbols parsed")
    if len(ref) < 3:
        print("  [fatal] fewer than 3 symbols parsed -- is that file gene symbols?")
        sys.exit(2)
    any_tested = False
    for name, _d, _s, _k in expect:
        for mv in sorted({r["genes"] for r in rows if r["list"] == name}):
            any_tested = True
            g = set(mv.split(","))
            inside, outside = sorted(g & ref), sorted(g - ref)
            if not inside:
                # a reference table for one list says nothing about the others
                print(f"\n  {name} (n={len(g)}): 0 overlap -- not this reference, skipped")
                continue
            print(f"\n  {name} (n={len(g)}): {len(inside)}/{len(g)} present in the reference")
            print(f"      present: {','.join(inside) or '-'}")
            print(f"      ABSENT : {','.join(outside) or '-'}")
            if not outside:
                print("      -> clean subset. Cite the reference, and state in the Methods")
                print("         how the subset was chosen -- a subset still needs a rule.")
            else:
                print("      -> NOT a clean subset. Each absent gene needs its own source,")
                print("         or is dropped, or the list is declared hand-assembled.")
    if not any_tested:
        print("  nothing to test -- no expected list was found in the scripts")

print(f"\n[write] {out_path}")
print("[next] for each list with no source: find the publication its membership")
print("       matches, re-run with --check-against that table, and record the outcome")
print("       in the Methods. A list that matches nothing is hand-assembled and must")
print("       say so.")
PY
rc=$?
rm -f "$EXPECT"
exit $rc
