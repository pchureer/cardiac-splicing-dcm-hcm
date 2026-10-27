#!/usr/bin/env python3
"""
108_rl125_rerun_audit.py -- read the nine SBATCH scripts that must be re-run at
125 bp and print, per script, exactly what has to change. EDITS NOTHING.

WHY THIS AND NOT A REWRITE
--------------------------
Each of those scripts points at `rmats_magnet/` (the 50 bp output) and several
carry hard-coded expected counts as gates. Rewriting them blind would be the
same mistake as inferring a third rMATS run from the fact that the paper
reports a third contrast. This script reads what is actually there, checks the
proposed replacement path against the directories that actually exist, and
hands back a patch plan for a human to approve.

WHAT IT REPORTS, PER SCRIPT
---------------------------
  1. every rMATS input directory referenced, and whether an rl125 counterpart
     exists on disk
  2. every `--out` prefix, and the proposed rl125 name
  3. every `--cov` / `--b1` / `--b2` / `--extra-cov-file` argument (these do
     NOT change -- they are sample metadata, not quantification)
  4. hard-coded numbers that match a known 50 bp gate value, with the 125 bp
     replacement where it is already measured
  5. whether the rMATS directory is held in a shell variable (one edit) or
     written out repeatedly (N edits)
  6. missing `--event-type` on an A3SS/A5SS call -- a silent failure

USAGE
-----
  python 108_rl125_rerun_audit.py                 # audit the default nine
  python 108_rl125_rerun_audit.py --plan plan.txt # also write the plan to a file
  python 108_rl125_rerun_audit.py 21_replication_v2.sbatch
  python 108_rl125_rerun_audit.py --self-test

Nothing is submitted and nothing is modified. Approve the plan, make the edits,
then submit by hand.
"""
import argparse
import glob
import os
import re
import sys

# the nine that re-run, from REMAINING_TASKS_14oct.md track 1
DEFAULT_SCRIPTS = [
    "04_phase2_v2_checks.sbatch",
    "08_composition_adjusted_discovery.sbatch",
    "72_run_cmfrac.sbatch",
    "75_cmfrac_allclasses.sbatch",
    "56_rbm20_targets45.sbatch",
    "23_signature_rbp_scan.sbatch",
    "21_replication_v2.sbatch",
    "36_matched_pair_panel.sbatch",
    "20_jc_vs_jcec_readlength.sbatch",
]

# 50 bp gate values -> 125 bp replacement, where measured.
# None means "measured value not yet available; this gate blocks".
GATE_VALUES = {
    "926": "813",      "231": "229",      "90": "87",
    "84": "74",        "114": "91",       "1445": "1294",
    "62971": "53456",  "62762": "62537",  "15292": "14346",
    "8988": "8360",    "6102": "6243",    "156115": "144942",
    "48024": "39607",  "468": "431",      "45380": "38273",
    "1355": "1430",    "909": "913",      "5301": "5401",
    "2592": "2373",    "62482": None,     "899": None,
    "835": None,       "866": None,       "888": None,
    "1054": None,      "716": None,       "590": None,
    "1074": None,      "215491": None,    "2163": None,
    "1636": None,      "526": None,       "927": None,
    "12725": None,     "38084": None,     "134": None,
}

# arguments that are sample metadata and must NOT be renamed
METADATA_ARGS = ("--cov", "--b1", "--b2", "--extra-cov-file", "--exclude")

# require a real suffix after rmats_ so bare `rmats` and the shell
# interpolation `rmats_${COH}` are not mistaken for directories
RMATS_DIR = re.compile(r"\b(rmats_[A-Za-z0-9][A-Za-z0-9_.]*)\b")
OUT_ARG = re.compile(r"--out\s+(\S+)")
META_ARG = re.compile(r"(--cov|--b1|--b2|--extra-cov-file|--exclude)\s+(\S+)")
NUMBER = re.compile(r"(?<![\w.])(\d[\d_]{1,8})(?![\w.])")
VAR_ASSIGN = re.compile(r"^\s*([A-Z][A-Z0-9_]*)=(.+)$")


def norm_num(s):
    return s.replace("_", "")


def rl125_dir_for(d, present):
    """Propose the rl125 counterpart of an rMATS directory, checked on disk."""
    if "rl125" in d:
        return d, "already rl125"
    cand = f"{d}_rl125"
    if cand in present:
        return cand, "exists"
    # a near-match must still be a real run: never a scratch or abandoned dir
    near = [p for p in present
            if "rl125" in p and d.split("_", 1)[-1] in p
            and not p.endswith("_tmp") and "abandoned" not in p]
    if near:
        return near[0], "exists (name differs -- CONFIRM)"
    return cand, "DOES NOT EXIST ON DISK"


def audit(path, present, out):
    def w(s=""):
        out.append(s)

    if not os.path.exists(path):
        w(f"═══ {path}")
        w("    NOT FOUND in this directory")
        w()
        return

    lines = open(path, errors="replace").read().splitlines()
    w(f"═══ {path}   ({len(lines)} lines)")

    # --- 1. rMATS directories, and whether they sit in a variable
    dirs, dir_lines = {}, {}
    varholds = {}
    for i, ln in enumerate(lines, 1):
        if ln.lstrip().startswith("#"):
            continue
        for d in RMATS_DIR.findall(ln):
            dirs.setdefault(d, 0)
            dirs[d] += 1
            dir_lines.setdefault(d, []).append(i)
        m = VAR_ASSIGN.match(ln)
        if m and RMATS_DIR.search(m.group(2)):
            varholds[m.group(1)] = (i, m.group(2).strip())

    if dirs:
        w("  rMATS input directories")
        for d, n in sorted(dirs.items()):
            new, status = rl125_dir_for(d, present)
            flag = "" if status in ("exists", "already rl125") else "   <-- CHECK"
            w(f"    {d:<28} {n:>2}x  lines {dir_lines[d][:8]}")
            w(f"      -> {new:<26} [{status}]{flag}")
        if varholds:
            w("    held in a shell variable (ONE edit each):")
            for v, (i, val) in sorted(varholds.items()):
                w(f"      line {i}: {v}={val}")
        else:
            hard = sum(n for d, n in dirs.items() if "rl125" not in d)
            if hard > 1:
                w(f"    ⚠ no variable holds it -- {hard} separate edits, "
                  f"and missing one gives a SILENT 50/125 bp mix")
    else:
        w("  no rMATS directory referenced")

    # --- 2. --out prefixes
    outs = [(i, m.group(1)) for i, ln in enumerate(lines, 1)
            if not ln.lstrip().startswith("#")
            for m in [OUT_ARG.search(ln)] if m]
    if outs:
        w("  --out prefixes")
        for i, o in outs:
            if "rl125" in o:
                w(f"    line {i:>4}: {o}   [already rl125]")
            else:
                w(f"    line {i:>4}: {o}  ->  {o}_rl125")
        w("    ⚠ renaming the output is NOT optional. Writing 125 bp results "
          "over the 50 bp filenames destroys the comparator the gates "
          "calibrate against.")

    # --- 3. metadata arguments: must not change
    metas = {}
    for i, ln in enumerate(lines, 1):
        if ln.lstrip().startswith("#"):
            continue
        for m in META_ARG.finditer(ln):
            metas.setdefault((m.group(1), m.group(2)), []).append(i)
    if metas:
        w("  sample metadata -- DO NOT RENAME (same samples at both read lengths)")
        for (a, v), ls in sorted(metas.items()):
            w(f"    {a} {v}   lines {ls[:6]}")

    # --- 4. hard-coded gate values
    hits = {}
    for i, ln in enumerate(lines, 1):
        for m in NUMBER.finditer(ln):
            n = norm_num(m.group(1))
            if n in GATE_VALUES:
                hits.setdefault(n, []).append((i, ln.strip()[:88]))
    if hits:
        w("  hard-coded 50 bp values (gate constants)")
        blocked = []
        for n, occ in sorted(hits.items(), key=lambda kv: -len(kv[1])):
            rep = GATE_VALUES[n]
            tag = f"-> {rep}" if rep else "-> NOT YET MEASURED (blocks)"
            if rep is None:
                blocked.append(n)
            w(f"    {n:<8} {tag:<30} {len(occ)}x")
            for i, txt in occ[:3]:
                w(f"        line {i:>4}: {txt}")
        if blocked:
            w(f"    ⚠ {len(blocked)} value(s) here have no 125 bp replacement yet: "
              f"{sorted(blocked)}")
            w("      This script CANNOT be re-baselined until track 1 produces them.")
            w("      Re-run it to GENERATE numbers; do not edit its gates yet.")

    # --- 5. --event-type on A3SS / A5SS
    for i, ln in enumerate(lines, 1):
        if ln.lstrip().startswith("#"):
            continue
        if ("A3SS" in ln or "A5SS" in ln) and "phase2_diffsplice" in ln \
                and "--event-type" not in ln:
            w(f"  ⚠ line {i}: phase2 call mentions A3SS/A5SS without "
              f"--event-type -- silent failure of that class")

    # --- 6. resource header, for planning concurrency
    res = [ln.strip() for ln in lines
           if ln.startswith("#SBATCH") and
           any(k in ln for k in ("-c", "--cpus", "--mem", "-t ", "--time",
                                 "-p ", "--partition", "--qos", "--array"))]
    if res:
        w("  resources as written")
        for r in res:
            w(f"    {r}")
        if any("--cpus" in r or re.search(r"-c\s*\d", r) for r in res):
            w("    note: phase2_diffsplice_v2.py is a SINGLE SERIAL LOOP. "
              "CPUs above 1 buy nothing; memory is the binding constraint.")
    w()


def self_test():
    import tempfile
    sb = """#!/bin/bash
#SBATCH -p cpu
#SBATCH -c 8
#SBATCH --mem=48G
# comment mentioning rmats_magnet should be ignored
MATSDIR=rmats_magnet
python phase2_diffsplice_v2.py --mats $MATSDIR/SE.MATS.JCEC.txt \\
    --b1 dcm_bams.txt --b2 nf_bams.txt --cov magnet_covariates.tsv \\
    --out SE_v2compadj
python phase2_diffsplice_v2.py --mats rmats_magnet/A3SS.MATS.JCEC.txt --out A3
EXPECT=926
"""
    d = tempfile.mkdtemp()
    p = os.path.join(d, "t.sbatch")
    open(p, "w").write(sb)
    out = []
    audit(p, {"rmats_magnet", "rmats_magnet_rl125"}, out)
    txt = "\n".join(out)
    assert "rmats_magnet" in txt
    assert "rmats_magnet_rl125" in txt and "[exists]" in txt
    assert "SE_v2compadj  ->  SE_v2compadj_rl125" in txt
    assert "--cov magnet_covariates.tsv" in txt
    assert "DO NOT RENAME" in txt
    assert "926" in txt and "-> 813" in txt
    assert "without --event-type" in txt, "A3SS check did not fire"
    assert "MATSDIR=rmats_magnet" in txt, "variable not detected"
    assert "SINGLE SERIAL LOOP" in txt
    os.remove(p)
    os.rmdir(d)
    print("[ok] self-test passed: paths, --out rename, metadata protection, "
          "gate values, A3SS --event-type, variable detection, resources")
    return 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("scripts", nargs="*", default=None)
    ap.add_argument("--plan", help="also write the report to this file")
    ap.add_argument("--self-test", action="store_true")
    a = ap.parse_args()

    if a.self_test:
        return self_test()

    present = {os.path.basename(p.rstrip("/"))
               for p in glob.glob("rmats*") if os.path.isdir(p)}
    # scratch and abandoned directories are not candidates for anything
    present = {p for p in present
               if not p.endswith("_tmp") and "abandoned" not in p}
    out = []
    out.append("=" * 78)
    out.append("rl125 RE-RUN AUDIT -- nothing is edited or submitted by this script")
    out.append("=" * 78)
    out.append(f"rMATS directories present: {sorted(present) or 'NONE FOUND'}")
    if not present:
        out.append("  ⚠ no rmats* directory here. cd to the analysis directory.")
    out.append("")

    for s in (a.scripts or DEFAULT_SCRIPTS):
        audit(s, present, out)

    out.append("=" * 78)
    out.append("ORDER OF WORK")
    out.append("=" * 78)
    out.append("1. For each script: repoint the rMATS directory AND rename --out.")
    out.append("   Leave --cov/--b1/--b2/--extra-cov-file/--exclude alone.")
    out.append("2. Submit. Do NOT touch any hard-coded gate value yet -- these runs")
    out.append("   GENERATE the replacements.")
    out.append("3. Recount collapse-first (chr:exonStart:exonEnd for SE and RI,")
    out.append("   key[:5] for MXE/A3SS/A5SS) and fill NUMBERS_125bp_REDERIVATION.md.")
    out.append("4. Only then re-baseline gates, one at a time, old value in a comment.")
    out.append("")
    out.append("phase2 is single-threaded: request 1 CPU and the memory one run")
    out.append("actually used (sacct -j <id> -o JobID,Elapsed,MaxRSS; MaxRSS is a")
    out.append("SAMPLED maximum, so leave headroom). Within 96 GB you can run")
    out.append("several concurrently.")

    txt = "\n".join(out)
    print(txt)
    if a.plan:
        open(a.plan, "w").write(txt + "\n")
        print(f"\n[ok] plan written to {a.plan}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
