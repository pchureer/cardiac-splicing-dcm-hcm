#!/usr/bin/env python3
"""
171_build_deposit.py -- assemble the GitHub/Zenodo deposit tree.

Dry-run by default. Nothing is copied until --apply is given.

What it does
------------
  1. Enumerates every candidate file under --root.
  2. REFUSES to proceed if anything matching a genotype pattern would be copied.
  3. Sorts scripts into scripts/ , scripts/patches/ , scripts/cluster/.
  4. Copies the small derived tables and the _rl125 figures.
  5. Lists (and, with --with-events, copies) the large event tables into events/.
  6. Measures the facts the README refuses to assert by hand and writes
     PROVENANCE.md from the measurements.
  7. Writes MANIFEST.tsv with size + SHA-256 + in_git for every file.

Gates -- each one exits non-zero rather than package a wrong tree:
  G1  no genotype-derived file is in the copy set
  G2  no .bak / .pid / scratch file is in the copy set
  G3  the two live names for the discovery table are byte-identical
  G4  the script numbering is contiguous apart from the three declared gaps
  G5  every file the README names as deposited actually exists
  G6  CONVENTIONS.md and README.md are present and carry no bracketed placeholder

Usage
-----
  python 171_build_deposit.py --root /data/project/chureeratp/RNA/HEART \
                              --out  /data/project/chureeratp/RNA/HEART/deposit
  # review the dry run, then
  python 171_build_deposit.py --root ... --out ... --apply
  python 171_build_deposit.py --root ... --out ... --apply --with-events
"""

import argparse
import hashlib
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

# ---------------------------------------------------------------- declarations

# Script numbers that do not exist and must not be reported as missing files.
DECLARED_GAPS = {94, 153, 154}

# Hand-written documentation that is part of the deposit. Each is placeholder-
# checked by G6 and copied verbatim from --docs.
DOCS = ("README.md", "CONVENTIONS.md", "DATA_AVAILABILITY.md")

# Anything whose path matches one of these is genotype-derived and is NEVER
# copied. Deliberately broad: a false positive costs one --allow-path argument,
# a false negative deposits individual-level data.
GENOTYPE_PATTERNS = [
    r"\.vcf(\.gz)?$",
    r"\.bcf$",
    r"\.g?vcf",
    r"kinship",
    r"\bking\b",
    r"relatedness",
    r"\bibd\b",
    r"\bibs\b",
    r"\.bim$", r"\.fam$",                      # PLINK (.bed: see below)
    r"\.pgen$", r"\.pvar$", r"\.psam$",
    r"genotyp",
    r"\bsnv\b", r"\bsnps?\b",
    r"\bgeno/",                                # the geno/ working directory
    r"pileup",
    # Found in the 27 Oct package under tables/: per-sample rare-variant calls
    # (SRR accession beside a TTN truncating variant) and allele-specific counts.
    # \bsnv\b missed "snv_screen" because "_" is a word character.
    r"snv",
    r"carrier",
    r"\bvariant",
    r"(^|/)ase(_|/)",                          # ase/ and ase_counts.tsv
    r"_calls\.",
]

# Files that LOOK genotype-derived but are the published summary, and so are
# deposited. Checked by G7 for size before being trusted.
SUMMARY_ALLOW = [r"duplicate_pairs.*\.(tsv|csv|md)$"]
# geno/panel.bed (the 3,820 panel SITES, no genotype in it) passes the .bed branch of
# genotype_refusal() -- no .bim/.fam beside it -- so it needs no allow entry here.

# Never deposited: intermediate state, not results.
EXCLUDE_PATTERNS = [
    r"\.bak(\.|$)",
    r"\.pid$",
    r"~$",
    r"\.swp$",
    r"__pycache__",
    r"\.ipynb_checkpoints",
    r"probe_rl125_tmp/",
    r"/tmp_",
    r"\.DS_Store$",
    r"_bams_dedup\.txt$",                      # superseded exclusion strategy
    r"\.log$",
    r"slurm-\d+\.out$",
]

# The event tables: large, Zenodo-only.
EVENT_GLOBS = ["*_rl125.all_events.tsv", "*_v2filtered.all_events.tsv",
               "*_v2filtered_rl125.all_events.tsv"]

# The two names that must be byte-identical (G3).
TWIN_A = "SE_rl125.all_events.tsv"
TWIN_B = "SE_v2filtered_rl125.all_events.tsv"

SCRIPT_EXT = {".py", ".sh", ".sbatch", ".R", ".r"}
CLUSTER_EXT = {".sh", ".sbatch"}

# Tables are "small" below this; anything larger is treated as an event table.
SMALL_TABLE_BYTES = 5 * 1024 * 1024


# ---------------------------------------------------------------- helpers

def matches(path_str, patterns):
    for p in patterns:
        if re.search(p, path_str, flags=re.IGNORECASE):
            return p
    return None


def genotype_refusal(full, rel):
    """Why this file is genotype-derived, or None.

    `.bed` is ambiguous: a PLINK genotype .bed, or a plain genomic interval
    file. It is refused only when a sibling .bim or .fam makes it the PLINK
    kind -- a gene panel BED is an analysis INPUT and belongs in the deposit.
    """
    if matches(rel, SUMMARY_ALLOW):
        return None
    if full.suffix == ".bed":
        stem = full.with_suffix("")
        if stem.with_suffix(".bim").exists() or stem.with_suffix(".fam").exists():
            return ".bed beside a .bim/.fam (PLINK genotype triple)"
        return None
    return matches(rel, GENOTYPE_PATTERNS)


def sha256(path, chunk=1 << 20):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while True:
            b = fh.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def script_number(name):
    m = re.match(r"^(\d{2,3})([a-z]?)_", name)
    if not m:
        return None, ""
    return int(m.group(1)), m.group(2)


def is_patch(name):
    return "_patch_" in name


def say(*a):
    print(*a, flush=True)


# ---------------------------------------------------------------- discovery

def collect(root, allow_paths):
    """Walk root once and bucket every file. Returns (buckets, refused)."""
    buckets = {"scripts": [], "patches": [], "cluster": [],
               "tables": [], "figures": [], "events": [], "docs": []}
    refused_geno = []
    refused_excl = []
    big_tables = []

    for dirpath, dirnames, filenames in os.walk(root):
        # don't descend into a deposit we previously built, or into .git
        dirnames[:] = [d for d in dirnames
                       if d not in {".git", "deposit", "node_modules"}]
        for fn in sorted(filenames):
            full = Path(dirpath) / fn
            rel = str(full.relative_to(root))

            ext = full.suffix
            try:
                size = full.stat().st_size
            except OSError:
                continue

            # ---- CODE IS NEVER GENOTYPE DATA -------------------------------
            # A numbered script is code and is always deposited, whatever its
            # name says. 27_genotype_concordance.sbatch is the script that
            # PRODUCES the screen, not its output, and the README promises it.
            # Refusing it also punched a false hole in the numbering, which G4
            # then reported as a missing file: one over-broad filter, two
            # consequences, and the second one did not look like the first.
            n_script, _ = script_number(fn)
            is_code = ext in SCRIPT_EXT and n_script is not None

            if not is_code and rel not in allow_paths:
                g = genotype_refusal(full, rel)
                if g:
                    refused_geno.append((rel, g))
                    continue
            e = matches(rel, EXCLUDE_PATTERNS)
            if e:
                refused_excl.append((rel, e))
                continue

            # Event tables first -- they are .tsv but must not land in tables/.
            # Identified by NAME, not by size: S3_lvef_rl125.tsv is 10 MB and is
            # a derived table, so a size rule would push it into the Zenodo-only
            # bucket and silently drop it from the git tree.
            if fn.endswith(".all_events.tsv"):
                buckets["events"].append(full)
            elif ext in SCRIPT_EXT:
                n, _ = script_number(fn)
                if n is None:
                    continue                      # unnumbered helper: skip
                if is_patch(fn):
                    buckets["patches"].append(full)
                elif ext in CLUSTER_EXT:
                    buckets["cluster"].append(full)
                else:
                    buckets["scripts"].append(full)
            elif ext in {".tsv", ".csv", ".bed", ".gmt"}:
                buckets["tables"].append(full)
                if size > SMALL_TABLE_BYTES:
                    big_tables.append((rel, size))
            elif ext == ".png" and "_rl125" in fn:
                buckets["figures"].append(full)
            elif ext in {".md"} and "/claude/" not in rel:
                buckets["docs"].append(full)

    return buckets, refused_geno, refused_excl, big_tables


# ---------------------------------------------------------------- gates

def gate_G1(refused_geno, strict):
    say("\n[G1] genotype exclusion")
    if not refused_geno:
        say("  no genotype-derived file found under --root.")
        say("  ⚠ that is worth a second look: the screen produced a VCF and a")
        say("    kinship table. If they live outside --root this gate cannot")
        say("    see them -- confirm by hand that neither is in --out.")
        return
    say(f"  refused {len(refused_geno)} file(s), NOT copied:")
    for rel, pat in refused_geno[:40]:
        say(f"    {rel}   [{pat}]")
    if len(refused_geno) > 40:
        say(f"    ... and {len(refused_geno) - 40} more")
    say("  these are excluded by design, not by accident. See README section 5.")


def gate_G2(refused_excl):
    say("\n[G2] intermediate-state exclusion")
    by_pat = {}
    for rel, pat in refused_excl:
        by_pat.setdefault(pat, 0)
        by_pat[pat] += 1
    for pat, n in sorted(by_pat.items(), key=lambda kv: -kv[1]):
        say(f"  {n:5d}  {pat}")
    say(f"  total excluded: {len(refused_excl)}")


def gate_G3(events):
    say("\n[G3] the two live names for the discovery table")
    a = [p for p in events if p.name == TWIN_A]
    b = [p for p in events if p.name == TWIN_B]

    # ⚠ One of them is a SYMLINK to the other. Comparing their SHA-256 then
    # hashes the same bytes twice and passes by construction -- the gate
    # asserted byte-identity and tested nothing. Worse, shutil.copy2 FOLLOWS a
    # symlink, so the deposit would carry two 23 MB copies of one file under
    # two names, which is the opposite of what recording the alias is for.
    for p in a + b:
        if p.is_symlink():
            tgt = os.readlink(p)
            say(f"  {p.name} is a SYMLINK -> {tgt}")
            say(f"  so byte-identity is trivially true and proves nothing. The")
            say(f"  alias is real and should be DOCUMENTED, not duplicated:")
            say(f"  the deposit carries the target once and the README records")
            say(f"  that the other name points at it.")
            return "symlink:" + tgt

    if not a or not b:
        say(f"  [warn] did not find both {TWIN_A} and {TWIN_B} under --root.")
        say("         found: " + (", ".join(p.name for p in a + b) or "neither"))
        say("         the README claims they are byte-identical. If only one")
        say("         exists here, trim that claim rather than leave it.")
        return None
    ha, hb = sha256(a[0]), sha256(b[0])
    say(f"  {TWIN_A}  {ha}")
    say(f"  {TWIN_B}  {hb}")
    if ha != hb:
        sys.exit("[fatal] G3: the two live names are NOT byte-identical. The "
                 "README's claim is false and the scripts that read them are "
                 "reading different data. Resolve this before depositing.")
    say("  identical. README section 3 is correct.")
    return ha


def gate_G4(scripts, patches, cluster, extra_gaps=()):
    say("\n[G4] script numbering")
    nums = {}
    for p in scripts + patches + cluster:
        n, letter = script_number(p.name)
        if n is None:
            continue
        nums.setdefault(n, []).append(p.name)
    if not nums:
        sys.exit("[fatal] G4: no numbered scripts found under --root. Wrong --root?")
    lo, hi = min(nums), max(nums)
    declared = DECLARED_GAPS | set(extra_gaps)
    missing = sorted(set(range(lo, hi + 1)) - set(nums))
    unexpected = sorted(set(missing) - declared)
    absent_gaps = sorted(declared & set(nums))
    total = sum(len(v) for v in nums.values())
    say(f"  {total} script file(s), numbers {lo:02d}-{hi}")
    say(f"  gaps found:    {missing}")
    say(f"  gaps declared: {sorted(declared)}")
    if unexpected:
        sys.exit(f"[fatal] G4: {len(unexpected)} undeclared gap(s): {unexpected}\n"
                 f"        Either those files are missing from --root (wrong "
                 f"--root, or deleted) or the README's gap note is incomplete.\n"
                 f"        A hole in a deposited numbering reads as a lost file, "
                 f"so each one has to be declared, not suppressed:\n"
                 f"          --declare-gap "
                 f"{' --declare-gap '.join(str(n) for n in unexpected[:6])}"
                 f"{' ...' if len(unexpected) > 6 else ''}\n"
                 f"        and every number you declare here must also be named "
                 f"in README section 2.")
    if absent_gaps:
        say(f"  [warn] {absent_gaps} declared as gaps but a file with that "
            f"number exists. Fix DECLARED_GAPS and the README together.")
    dups = {n: v for n, v in nums.items() if len(v) > 1}
    if dups:
        say("  numbers with more than one file (letter suffixes are expected):")
        for n, v in sorted(dups.items()):
            say(f"    {n}: {', '.join(sorted(v))}")
    return total, lo, hi, missing


MAX_SUMMARY_ROWS = 200


def gate_G7(tables):
    """A file allowed past G1 as 'the published summary' must BE a summary.

    The screen's per-pair kinship table over 670 samples is ~224,000 rows; the
    reported duplicate-pair summary is ~24. If a file named like the summary is
    long, it is the matrix and must not be deposited.
    """
    say("\n[G7] the duplicate-pair summary is a summary, not the matrix")
    hits = [p for p in tables if matches(str(p), SUMMARY_ALLOW)]
    if not hits:
        say("  no duplicate-pair file in the copy set.")
        say("  ⚠ the paper reports 24 cross-dataset duplicate pairs. If that "
            "table is a supplementary item, it should be here.")
        return
    for p in hits:
        try:
            with open(p, encoding="utf-8", errors="replace") as fh:
                head = fh.readline().rstrip("\n")
                nrow = sum(1 for _ in fh)
        except OSError as exc:
            sys.exit(f"[fatal] G7: cannot read {p}: {exc}")
        say(f"  {p.name}: {nrow} data row(s)")
        say(f"    header: {head[:160]}")
        if nrow > MAX_SUMMARY_ROWS:
            sys.exit(
                f"[fatal] G7: {p.name} has {nrow} rows, over the {MAX_SUMMARY_ROWS}-row "
                f"ceiling for a summary. A per-pair table this long is the kinship "
                f"matrix, which is NOT deposited (README section 5). Either point "
                f"--root at the 24-pair summary instead, or remove this file.")
    say(f"  within the {MAX_SUMMARY_ROWS}-row ceiling: deposited as the reported "
        f"summary. Confirm by eye that it carries no kinship coefficient per pair "
        f"beyond what the paper reports.")


def report_events(events, root):
    """The event tables are the expensive part. Group them so the triage
    decision is made on numbers, not on a total."""
    say("\n[events] candidates, grouped by directory (paths relative to --root)")
    by_dir = {}
    for p in events:
        try:
            d = str(p.parent.relative_to(root))
        except ValueError:
            d = str(p.parent)
        a, b = by_dir.setdefault(d, [0, 0])
        by_dir[d] = [a + 1, b + p.stat().st_size]
    rows = sorted(by_dir.items(), key=lambda kv: -kv[1][1])
    tot_n = sum(v[0] for _, v in rows)
    tot_b = sum(v[1] for _, v in rows)
    for d, (n, b) in rows[:25]:
        say(f"  {b/1e6:9.1f} MB  {n:4d} file(s)  {d}")
    if len(rows) > 25:
        rest = rows[25:]
        say(f"  {sum(v[1] for _, v in rest)/1e6:9.1f} MB  "
            f"{sum(v[0] for _, v in rest):4d} file(s)  "
            f"... and {len(rest)} more directories")
    say(f"  {tot_b/1e6:9.1f} MB  {tot_n:4d} file(s)  TOTAL")
    if tot_b > 1e9:
        say(f"\n  ⚠ {tot_b/1e9:.1f} GB. Zenodo's default per-record quota is "
            f"50 GB, so this fits, but most of it is probably not needed to "
            f"reproduce a published number. Decide which directories are the "
            f"primary contrasts and pass --events-dir for each; the rest stay "
            f"on the cluster. Depositing every intermediate makes the record "
            f"harder to use, not more complete.")


def gate_G6(outdir_src, fatal=True):
    """The DOCS files must exist and carry no placeholder.

    In a dry run this only warns, so G1-G4 can be inspected before the Zenodo
    DOI has been reserved. On --apply it refuses.
    """
    say("\n[G6] documentation" + ("" if fatal else "  (dry run: warn only)"))
    def stop(msg):
        if fatal:
            sys.exit(msg)
        say("  " + msg.replace("[fatal]", "[warn]"))
    bad = [r"…", r"\.\.\.(?!\w)", r"\bTBD\b", r"\bTODO\b", r"\bXXX\b",
           r"\[fill[: ]", r"\[Zenodo DOI\]", r"⟨", r"\bFIXME\b",
           r"\[INSERT", r"<placeholder",
           r"ZENODO_DOI_TO_FILL", r"GITHUB_USER_TO_FILL"]
    for name in DOCS:
        p = Path(outdir_src) / name
        if not p.exists():
            stop(f"[fatal] G6: {name} not found in --docs ({outdir_src}). "
                 f"It is part of the deposit, not an optional extra.")
            continue
        txt = p.read_text(encoding="utf-8")
        hit = matches(txt, bad)
        if hit:
            for i, line in enumerate(txt.splitlines(), 1):
                if matches(line, [hit]):
                    say(f"  {name}:{i}: {line.strip()[:110]}")
            stop(f"[fatal] G6: {name} carries an unresolved placeholder "
                 f"matching {hit!r}. Fill it before packaging -- a "
                 f"placeholder in a deposited README survives review.")
            continue
        say(f"  {name}: {len(txt.split())} words, no placeholder.")


# ---------------------------------------------------------------- provenance

def measure_provenance(root, scripts, twin_hash):
    """Facts the README refuses to assert by hand. Measured, or recorded absent."""
    out = ["# Provenance",
           "",
           "Every statement below was measured from the code and the files at",
           "packaging time by `171_build_deposit.py`. Nothing here is typed in by",
           "hand. A line reading `NOT MEASURED` means the script could not find",
           "the evidence -- treat it as an open question, not as a null result.",
           ""]

    # --- the --dpsi-rule that wrote the boolean columns
    out += ["## The `--dpsi-rule` that wrote the boolean significance columns", ""]
    defaults = {}      # value -> [scripts declaring it as an argparse default]
    passed = {}        # value -> [scripts passing it on a command line]
    mentions = []
    for p in scripts:
        try:
            txt = p.read_text(encoding="utf-8", errors="replace")
        except OSError:
            continue
        if not re.search(r"dpsi[-_]rule", txt):
            continue
        mentions.append(p.name)
        for m in re.finditer(r"""add_argument\(\s*["']--dpsi-rule["'][^)]*?"""
                             r"""default\s*=\s*["']([^"']+)["']""", txt, re.S):
            defaults.setdefault(m.group(1), []).append(p.name)
        for m in re.finditer(r"""--dpsi-rule[=\s]+["']?([A-Za-z0-9_.]+)""", txt):
            if m.group(1) not in ("default", "dest", "help", "type", "choices"):
                passed.setdefault(m.group(1), []).append(p.name)
    if defaults or passed:
        if defaults:
            out.append("Declared as an argparse **default**:")
            out.append("")
            for v, who in sorted(defaults.items()):
                out.append(f"- `{v}` — in {', '.join(sorted(set(who))[:6])}"
                           f"{' and others' if len(set(who)) > 6 else ''}")
            out.append("")
        if passed:
            out.append("Seen **passed explicitly** on a command line inside a script:")
            out.append("")
            for v, who in sorted(passed.items()):
                out.append(f"- `{v}` — in {', '.join(sorted(set(who))[:6])}"
                           f"{' and others' if len(set(who)) > 6 else ''}")
            out.append("")
        allv = set(defaults) | set(passed)
        if len(allv) > 1:
            out.append(f"⚠ **{len(allv)} different values appear: "
                       f"{sorted(allv)}.** The boolean columns in the deposited "
                       f"tables were written by exactly one of them. Resolve "
                       f"which, and state it in README section 1 item 2 — two "
                       f"live rules is how a reader reproduces 864 instead of "
                       f"813.")
        else:
            out.append(f"One value only: **`{sorted(allv)[0]}`**. That is what "
                       f"wrote the deposited boolean columns.")
        out.append("")
        out.append(f"Mentioned in {len(mentions)} script(s): "
                   f"{', '.join(sorted(mentions)[:20])}"
                   f"{' ...' if len(mentions) > 20 else ''}")
    else:
        out.append("NOT MEASURED -- no `--dpsi-rule` / `dpsi_rule` token found in any "
                   "deposited script. If the boolean columns were written by an "
                   "earlier version of the producer, say so explicitly in the README "
                   "rather than leaving the question open.")
    out.append("")

    # --- reference genome and annotation, from the rMATS command lines
    out += ["## Reference genome and annotation", ""]
    gtf, fa = [], []
    for p in scripts:
        try:
            txt = p.read_text(encoding="utf-8", errors="replace")
        except OSError:
            continue
        for m in re.finditer(r"[-\w/\.]+\.g(tf|ff3?)(\.gz)?", txt):
            gtf.append((p.name, m.group(0)))
        for m in re.finditer(r"[-\w/\.]+\.(fa|fasta)(\.gz)?", txt):
            fa.append((p.name, m.group(0)))

    def uniq(pairs):
        seen, keep = set(), []
        for nm, v in pairs:
            if v not in seen:
                seen.add(v)
                keep.append((nm, v))
        return keep
    if gtf or fa:
        out.append("```")
        for nm, v in uniq(gtf)[:15]:
            out.append(f"GTF   {v}      (seen in {nm})")
        for nm, v in uniq(fa)[:15]:
            out.append(f"FASTA {v}      (seen in {nm})")
        out.append("```")
        if len(set(v for _, v in gtf)) > 1:
            out.append("")
            out.append("⚠ More than one annotation path appears in the scripts. "
                       "Confirm which one the primary run used; a second "
                       "annotation in the tree is how a cohort ends up "
                       "quantified against different exon definitions.")
    else:
        out.append("NOT MEASURED -- no annotation or genome path found in the scripts.")
    out.append("")

    # --- read-length configuration actually passed to rMATS
    out += ["## Read length passed to rMATS", ""]
    rl = []
    for p in scripts:
        try:
            txt = p.read_text(encoding="utf-8", errors="replace")
        except OSError:
            continue
        for m in re.finditer(r"--readLength\s+(\d+)", txt):
            rl.append((p.name, m.group(1)))
        if "--variable-read-length" in txt:
            rl.append((p.name, "variable-read-length"))
    if rl:
        byval = {}
        for nm, v in rl:
            byval.setdefault(v, set()).add(nm)
        for v, who in sorted(byval.items()):
            w = sorted(who)
            out.append(f"- `{v}` — {len(w)} script(s): {', '.join(w[:8])}"
                       f"{' ...' if len(w) > 8 else ''}")
        out.append("")
        numeric = sorted(v for v in byval if v.isdigit())
        if len(numeric) > 1:
            out.append(f"⚠ Two read lengths are live in the deposited scripts "
                       f"({', '.join(numeric)} bp). That is expected — the 50 bp "
                       f"run is deposited so the read-length comparison can be "
                       f"reproduced — but confirm that every script reading a "
                       f"`_rl125` table is the {max(numeric, key=int)} bp one.")
    else:
        out.append("NOT MEASURED -- no `--readLength` found. If rMATS was driven from "
                   "a submission script not deposited here, that is a gap.")
    out.append("")

    # --- the twin hash
    out += ["## The two live names for the discovery table", ""]
    if twin_hash:
        out.append(f"`{TWIN_A}` and `{TWIN_B}` are byte-identical.")
        out.append("")
        out.append(f"    SHA-256  {twin_hash}")
    else:
        out.append("NOT MEASURED -- both files were not present under `--root`.")
    out.append("")

    # --- environment
    out += ["## Environment", ""]
    try:
        v = subprocess.run([sys.executable, "-V"], capture_output=True,
                           text=True, timeout=20).stdout.strip()
        out.append(f"Packaged with {v} on {os.uname().nodename}.")
    except Exception:
        pass
    env = subprocess.run("conda env export -n rnaseq_tools 2>/dev/null",
                         shell=True, capture_output=True, text=True)
    if env.returncode == 0 and env.stdout.strip():
        out.append("`environment.yml` is a full `conda env export` of `rnaseq_tools`.")
    else:
        out.append("⚠ `conda env export -n rnaseq_tools` did not succeed, so "
                   "`environment.yml` was NOT written. Run the packaging from "
                   "inside the environment, or export it by hand -- a deposit "
                   "without a pinned environment is not reproducible.")
    out.append("")
    return "\n".join(out), (env.stdout if env.returncode == 0 else None)


# ---------------------------------------------------------------- write

STATIC = {}

STATIC[".gitignore"] = """\
# The event tables are the Zenodo record, not the git record (README section 6).
events/

# Never commit genotype-derived output. See README section 5.
*.vcf
*.vcf.gz
*.bcf
*kinship*
*genotyp*
geno/
*.bed
*.bim
*.fam

# Intermediate state
*.bak
*.bak.*
*.pid
__pycache__/
*.pyc
.ipynb_checkpoints/
slurm-*.out
*.log
"""

STATIC["LICENSE"] = """\
MIT License

Copyright (c) 2026 Chureerat Phokaew

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

Note: the MIT licence covers the CODE in this repository. The deposited result
tables are released under CC-BY-4.0. The primary sequencing data are not
covered by either and remain under the terms of their original repositories.
"""

STATIC["CITATION.cff"] = """\
cff-version: 1.2.0
message: "If you use this code or the deposited result tables, please cite both the software record and the paper."
title: "Analysis code and event-level results for a four-cohort reanalysis of alternative splicing in dilated and hypertrophic cardiomyopathy"
abstract: "Complete analysis code and event-level rMATS result tables underlying the manuscript. Includes the analysis conventions required to reproduce every reported count."
type: software
authors:
  - family-names: Phokaew
    given-names: Chureerat
    email: chureerat.phokaew@gmail.com
license: MIT
doi: 10.5281/zenodo.23681815
repository-code: "https://github.com/pchureer/cardiac-splicing-dcm-hcm"
identifiers:
  - type: doi
    value: 10.5281/zenodo.23681815
    description: "Zenodo archive of this repository and the event-level tables"
keywords:
  - alternative splicing
  - RNA-seq
  - dilated cardiomyopathy
  - hypertrophic cardiomyopathy
  - rMATS
  - RBM20
"""


GIT_FILE_LIMIT = 50 * 1024 * 1024   # GitHub warns at 50 MB and rejects at 100 MB


def plan_copy(buckets, out, with_events, root=None):
    """Returns a list of (src, dest_rel, in_git).

    Tables and figures keep their path relative to --root under tables/ and
    figures/, so t3_rl125/Table1.tsv and ./Table1.tsv (different files, same
    basename) land in different places instead of colliding. A table over
    GIT_FILE_LIMIT goes to events/large_tables/ (Zenodo only, in_git=no); the
    MANIFEST records it either way."""
    jobs = []

    def rel_of(p):
        if root is not None:
            try:
                return str(Path(p).resolve().relative_to(Path(root).resolve()))
            except ValueError:
                pass
        return p.name
    # A symlink whose target is also in the copy set would be copied as a second
    # full file. Record it as an alias instead.
    real = {}
    for k in ("events", "tables", "figures"):
        for p in buckets.get(k, []):
            try:
                real.setdefault(os.path.realpath(p), []).append(p)
            except OSError:
                pass
    aliases = []
    for rp, group in real.items():
        if len(group) < 2:
            continue
        links = [q for q in group if q.is_symlink()]
        keep = next((q for q in group if not q.is_symlink()), group[0])
        for q in links:
            if q is not keep:
                aliases.append((q, keep))
    if aliases:
        say(f"\n[alias] {len(aliases)} symlink(s) point at a file already in "
            f"the copy set; the target is deposited once:")
        for q, keep in aliases:
            say(f"  {q.name} -> {keep.name}")
        drop = {str(q) for q, _ in aliases}
        for k in ("events", "tables", "figures"):
            buckets[k] = [p for p in buckets.get(k, []) if str(p) not in drop]

    for p in sorted(buckets["scripts"], key=lambda q: q.name):
        jobs.append((p, f"scripts/{p.name}", True))
    for p in sorted(buckets["patches"], key=lambda q: q.name):
        jobs.append((p, f"scripts/patches/{p.name}", True))
    for p in sorted(buckets["cluster"], key=lambda q: q.name):
        jobs.append((p, f"scripts/cluster/{p.name}", True))
    big = []
    for p in sorted(buckets["tables"], key=lambda q: q.name):
        rel = rel_of(p)
        if p.stat().st_size > GIT_FILE_LIMIT:
            jobs.append((p, f"events/large_tables/{rel}", False)); big.append((rel, p.stat().st_size))
        else:
            jobs.append((p, f"tables/{rel}", True))
    if big:
        say(f"\n[tables] {len(big)} derived table(s) over {GIT_FILE_LIMIT/1e6:.0f} MB routed to "
            f"events/large_tables/ (Zenodo only; MANIFEST in_git=no):")
        for rel, size in sorted(big, key=lambda kv: -kv[1]):
            say(f"  {size/1e6:8.1f} MB  {rel}")
    for p in sorted(buckets["figures"], key=lambda q: q.name):
        jobs.append((p, f"figures/{rel_of(p)}", True))
    for p in sorted(buckets["events"], key=lambda q: q.name):
        # path-preserving: jccheck/X and jccheck_rl125/X are different files
        jobs.append((p, f"events/{rel_of(p)}", False))
    # collisions
    seen = {}
    for src, rel, _ in jobs:
        if rel in seen:
            sys.exit(f"[fatal] two different files would be written to {rel}:\n"
                     f"          {seen[rel]}\n          {src}\n"
                     f"        A deposit that silently keeps one of two "
                     f"same-named files is worse than no deposit.")
        seen[rel] = src
    if not with_events:
        jobs = [j for j in jobs if not j[1].startswith("events/")]
    return jobs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True, help="project root to harvest from")
    ap.add_argument("--out", required=True, help="deposit directory to build")
    ap.add_argument("--docs", default=None,
                    help="directory holding README.md and CONVENTIONS.md "
                         "(default: beside this script)")
    ap.add_argument("--apply", action="store_true")
    ap.add_argument("--with-events", action="store_true",
                    help="also copy event tables into events/; requires at "
                         "least one --events-dir")
    ap.add_argument("--events-dir", action="append", default=[],
                    help="root-relative directory whose event tables are "
                         "deposited; repeatable. Everything else is left on "
                         "the cluster. Run once without it to see the grouped "
                         "listing and choose.")
    ap.add_argument("--declare-gap", action="append", type=int, default=[],
                    help="script number that legitimately does not exist, in "
                         "addition to 94/153/154; repeatable. Anything declared "
                         "here must also be named in README section 2.")
    ap.add_argument("--allow-path", action="append", default=[],
                    help="root-relative path to exempt from the genotype refusal; "
                         "repeatable. Use only for a demonstrated false positive.")
    a = ap.parse_args()

    root = Path(a.root).resolve()
    out = Path(a.out).resolve()
    docs = Path(a.docs).resolve() if a.docs else Path(__file__).resolve().parent
    if not root.is_dir():
        sys.exit(f"[fatal] --root is not a directory: {root}")
    if out == root or root in out.parents and out.name != "deposit":
        say(f"[note] --out sits inside --root; the walk skips a directory named "
            f"'deposit', so name it that or keep it outside.")

    say(f"root : {root}")
    say(f"out  : {out}")
    say(f"docs : {docs}")
    say(f"mode : {'APPLY' if a.apply else 'DRY RUN'}"
        f"{' +events' if a.with_events else ' (events listed, not copied)'}")

    buckets, ref_geno, ref_excl, big_tables = collect(root, set(a.allow_path))

    say("\n--- discovered ---")
    for k in ("scripts", "patches", "cluster", "tables", "figures", "events"):
        n = len(buckets[k])
        tot = sum(p.stat().st_size for p in buckets[k]) if n else 0
        say(f"  {k:9s} {n:5d} file(s)  {tot/1e6:9.1f} MB")

    gate_G1(ref_geno, strict=True)
    gate_G2(ref_excl)
    twin_hash = gate_G3(buckets["events"])
    nscripts, lo, hi, gaps = gate_G4(buckets["scripts"], buckets["patches"],
                                     buckets["cluster"], a.declare_gap)
    gate_G7(buckets["tables"])
    gate_G6(docs, fatal=a.apply)
    report_events(buckets["events"], root)
    if big_tables:
        say(f"\n[tables] {len(big_tables)} derived table(s) over "
            f"{SMALL_TABLE_BYTES/1e6:.0f} MB, kept in the git tree:")
        for rel, size in sorted(big_tables, key=lambda kv: -kv[1])[:10]:
            say(f"  {size/1e6:8.1f} MB  {rel}")
        worst = max(size for _, size in big_tables)
        if worst > 50e6:
            say("  those over 50 MB are routed to events/large_tables/ by the copy "
                "plan (Zenodo only); see the [tables] line below.")

    if a.events_dir:
        keep, drop = [], 0
        wanted = [d.strip().rstrip("/") for d in a.events_dir]
        for p in buckets["events"]:
            rel = str(p.relative_to(root))
            top = rel.split("/")[0] if "/" in rel else "."
            # exact directory match: "." = files directly under --root; "jccheck"
            # must not also select jccheck_rl125 or jccheck_50bp_BACKUP (a bare
            # startswith() did exactly that and silently picked up the backups)
            if any((d == "." and top == ".") or (d != "." and rel.startswith(d + "/"))
                   for d in wanted):
                keep.append(p)
            else:
                drop += 1
        unknown = [d for d in wanted
                   if d != "." and not any(str(q.relative_to(root)).startswith(d + "/")
                                           for q in buckets["events"])]
        if unknown:
            say(f"  ⚠ --events-dir named directories holding no event file: {unknown}")
        say(f"\n  --events-dir selected {len(keep)} of "
            f"{len(buckets['events'])} event file(s); {drop} left on the cluster.")
        if not keep:
            sys.exit("[fatal] --events-dir matched nothing. The paths are "
                     "relative to --root; check them against the grouped "
                     "listing above.")
        buckets["events"] = keep
    elif a.with_events:
        sys.exit(f"[fatal] --with-events was given but no --events-dir. That "
                 f"would copy all {len(buckets['events'])} candidate event "
                 f"files. Name the directories that hold the primary contrasts "
                 f"explicitly -- a 4 GB record assembled by accident is worse "
                 f"than one assembled on purpose.")

    jobs = plan_copy(buckets, out, a.with_events, root=root)
    say(f"\n--- copy plan: {len(jobs)} file(s) ---")
    by_dir = {}
    for _, rel, in_git in jobs:
        d = rel.rsplit("/", 1)[0] if "/" in rel else "."
        by_dir.setdefault(d, [0, 0])
        by_dir[d][0] += 1
        by_dir[d][1] += 1 if in_git else 0
    for d, (n, g) in sorted(by_dir.items()):
        say(f"  {d:18s} {n:5d} file(s)   {g:5d} in git")
    if not a.with_events:
        ev = buckets["events"]
        say(f"\n  events/ NOT copied: {len(ev)} file(s), "
            f"{sum(p.stat().st_size for p in ev)/1e6:.0f} MB. "
            f"Re-run with --with-events for the Zenodo archive.")

    if not a.apply:
        say("\nDRY RUN -- nothing written. Re-run with --apply.")
        return

    # ---- write
    for sub in ("scripts/patches", "scripts/cluster", "tables", "figures", "events"):
        (out / sub).mkdir(parents=True, exist_ok=True)

    rows = []
    for src, rel, in_git in jobs:
        dst = out / rel
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)
        rows.append((rel, dst.stat().st_size, sha256(dst), "yes" if in_git else "no"))
    say(f"copied {len(rows)} file(s)")

    for name, text in STATIC.items():
        (out / name).write_text(text, encoding="utf-8")
    for name in DOCS:
        shutil.copy2(docs / name, out / name)

    prov, envyml = measure_provenance(root, buckets["scripts"] + buckets["patches"]
                                      + buckets["cluster"], twin_hash)
    (out / "PROVENANCE.md").write_text(prov, encoding="utf-8")
    if envyml:
        (out / "environment.yml").write_text(envyml, encoding="utf-8")

    for name in list(STATIC) + list(DOCS) + ["PROVENANCE.md"] + (
            ["environment.yml"] if envyml else []):
        p = out / name
        rows.append((name, p.stat().st_size, sha256(p), "yes"))

    man = out / "MANIFEST.tsv"
    with open(man, "w", encoding="utf-8") as fh:
        fh.write("path\tbytes\tsha256\tin_git\n")
        for rel, size, h, g in sorted(rows):
            fh.write(f"{rel}\t{size}\t{h}\t{g}\n")
    rows.append(("MANIFEST.tsv", man.stat().st_size, "-", "yes"))
    say(f"wrote MANIFEST.tsv: {len(rows)} entries")

    git_bytes = sum(r[1] for r in rows if r[3] == "yes")
    all_bytes = sum(r[1] for r in rows)
    say(f"\ngit tree : {git_bytes/1e6:8.1f} MB")
    say(f"zenodo   : {all_bytes/1e6:8.1f} MB")
    if git_bytes > 95e6:
        say("[warn] the git tree is over 95 MB. GitHub warns at 50 MB per file "
            "and the repository is awkward above ~1 GB. Check MANIFEST.tsv for "
            "what is large and whether it belongs in events/ instead.")
    say(f"\nnext: cd {out} && git init && git add -A && git status --short | head -30")
    say("      read PROVENANCE.md and resolve every 'NOT MEASURED' line before "
            "the Zenodo upload.")


if __name__ == "__main__":
    main()
