#!/usr/bin/env python3
"""
142_patch_20_jc_readlength.py -- script 20 at 125 bp, so S5(a) can be built.

WHY THIS MUST HAPPEN BEFORE 20 IS RE-RUN
----------------------------------------
20 hard-codes `rmats_magnet/` and `rmats_HCMvsNF/` and writes every output to
an unsuffixed `jccheck/`. Running it with 125 bp inputs would overwrite the
50 bp JC tables in place — the `direction_concordance_dedup.tsv` defect of
16 Oct, except `jccheck/` is the ONLY copy of the 50 bp JC quantification and
script 20 is the only thing that produces it.

THREE SUFFIX RULES, AGAIN
-------------------------
  rmats_magnet/            MAGNet      -> rmats_magnet<SUF>/
  rmats_HCMvsNF/           MAGNet      -> rmats_HCMvsNF<SUF>/
  rmats_Brigham/ etc.      cohort      -> UNCHANGED, never re-quantified
  splitctl/{ET}_DCMvsA_1   MAGNet      -> ..._1<SUF>   (built by 136)
  splitctl/{ET}_HCMvsB_1   MAGNet      -> ..._1<SUF>
  {ET}_v2filtered          MAGNet      -> {ET}_v2filtered<SUF>
  {C}_{ET}_HCM             cohort      -> UNCHANGED
  jccheck/                 output      -> jccheck<SUF>/

PART C IS A FREE CALIBRATION
----------------------------
The three cohort arms read ONLY `rmats_${COH}/${ET}.MATS.JC.txt` — the
replication cohorts' own quantifications, which were never re-quantified.
So part C at 125 bp must produce files **byte-identical** to the 50 bp ones.
The patch adds that `cmp` at the end of the run. If they differ, phase 2 is
not reproducible and nothing else in this job should be believed. It costs
three extra phase-2 calls and buys a real control.

THE HEREDOC IS QUOTED
---------------------
The summary block opens with `python3 - <<'PY'`, so bash does NOT substitute
into it. SUF is passed by environment — the mechanism 127 had to use for
script 56's second heredoc — not by interpolation.

SCOPE: THREE COHORTS, NOT FOUR
------------------------------
Line 101 is `for COH in Brigham Fuwai Stanford`. **Colorado has no JC table
and 20 does not make one** — it is the DCM-vs-NF contrast and sits outside
20's HCM loop. 62's `REPL` and `MANUSCRIPT_SLOPE` both include Colorado, so
S5(a) will cover three of four cohorts and must say so. Adding Colorado means
a fourth arm with its own rMATS source, bam lists and --min-samples; it is
not done here.

ET=SE is hard-coded in 20, so jccheck is SE-only at both read lengths and
S5(a) is an SE-only comparison regardless.

MEMORY
------
20 asks for 32G. phase 2's measured MaxRSS is 1.57 GB (job 408099), so this
drops to 8G — still a 5x margin. The 32G ask is where the
QOSMaxMemoryPerUser queue waits came from.

USAGE
-----
  python 142_patch_20_jc_readlength.py [--apply] [--self-test]
"""
import argparse
import difflib
import os
import shutil
import subprocess
import sys
import time

TARGET = "20_jc_vs_jcec_readlength.sbatch"

# ------------------------------------------------------------------- edit 1
A_HDR = "#SBATCH --mem=32G\n"
B_HDR = ("#SBATCH --mem=8G          # was 32G; phase 2's measured MaxRSS is\n"
         "                          # 1.57 GB (job 408099), a 5x margin here\n")

# ------------------------------------------------------------------- edit 2
A_SETUP = '''ET=SE
declare -A MINS=( [Brigham]=4 [Fuwai]=7 [Stanford]=5 )

for f in rmats_magnet/${ET}.MATS.JC.txt magnet_covariates.tsv dcm_bams.txt nf_bams.txt; do
  [ -s "$f" ] || { echo "missing $f"; exit 1; }
done
mkdir -p jccheck
'''
B_SETUP = '''ET=SE
declare -A MINS=( [Brigham]=4 [Fuwai]=7 [Stanford]=5 )

# ---------------------------------------------------------------- read length
# Only the two MAGNet rMATS directories move. The three cohort directories
# (rmats_Brigham/ etc.) do NOT -- no replication cohort was re-quantified --
# which is why part C below must come out byte-identical at both read lengths,
# and why the run ends by checking exactly that.
SUF=${TBL_SUFFIX:-}
MAGDIR=rmats_magnet${SUF}
HCMDIR=rmats_HCMvsNF${SUF}
JCDIR=jccheck${SUF}
FORCE=${FORCE:-0}
echo "read length suffix : '${SUF}'   (empty = 50 bp)"
echo "MAGNet DCM-vs-NF   : ${MAGDIR}"
echo "MAGNet HCM-vs-NF   : ${HCMDIR}"
echo "output directory   : ${JCDIR}"
echo "cohort rMATS dirs  : rmats_<COH>/   (unsuffixed, deliberately)"

for f in ${MAGDIR}/${ET}.MATS.JC.txt magnet_covariates.tsv dcm_bams.txt nf_bams.txt; do
  [ -s "$f" ] || { echo "missing $f"; exit 1; }
done

# Never overwrite an existing set. jccheck/ is the ONLY copy of the 50 bp JC
# quantification and this script is the only thing that makes it.
if [ -n "$(ls -A "$JCDIR" 2>/dev/null)" ] && [ "$FORCE" != "1" ]; then
  echo
  echo "REFUSING — ${JCDIR}/ already contains files:"
  ls -1 "$JCDIR" | sed 's/^/    /'
  echo "  Re-run with FORCE=1 only after checking nothing downstream has"
  echo "  consumed them. S5(a) reads this directory."
  exit 1
fi
mkdir -p "$JCDIR"
'''

# ------------------------------------------------------------------- edit 3
A_ARM_A = '''python phase2_diffsplice_v2.py \\
  --mats rmats_magnet/${ET}.MATS.JC.txt --event-type $ET \\
  --b1 dcm_bams.txt --b2 nf_bams.txt --cov magnet_covariates.tsv \\
  --min-samples 20 --out jccheck/${ET}_v2filtered_JC
'''
B_ARM_A = '''python phase2_diffsplice_v2.py \\
  --mats ${MAGDIR}/${ET}.MATS.JC.txt --event-type $ET \\
  --b1 dcm_bams.txt --b2 nf_bams.txt --cov magnet_covariates.tsv \\
  --min-samples 20 --out ${JCDIR}/${ET}_v2filtered_JC
'''

# ------------------------------------------------------------------- edit 4
A_ARM_B = '''if [ -s splitctl/nfA_1.ids ] && [ -s splitctl/nfB_1.ids ] && \\
   [ -s rmats_HCMvsNF/${ET}.MATS.JC.txt ]; then
  python phase2_diffsplice_v2.py \\
    --mats rmats_magnet/${ET}.MATS.JC.txt --event-type $ET \\
    --b1 dcm_bams.txt --b2 nf_bams.txt --cov magnet_covariates.tsv \\
    --group1 DCM --group2 NF --min-samples 15 --exclude splitctl/nfB_1.ids \\
    --out jccheck/${ET}_DCMvsA_JC
  python phase2_diffsplice_v2.py \\
    --mats rmats_HCMvsNF/${ET}.MATS.JC.txt --event-type $ET \\
    --b1 hcm_bams.txt --b2 nf_bams.txt --cov hcm_nf_covariates.tsv \\
    --group1 HCM --group2 NF --min-samples 15 --exclude splitctl/nfA_1.ids \\
    --out jccheck/${ET}_HCMvsB_JC
'''
B_ARM_B = '''if [ -s splitctl/nfA_1.ids ] && [ -s splitctl/nfB_1.ids ] && \\
   [ -s ${HCMDIR}/${ET}.MATS.JC.txt ]; then
  python phase2_diffsplice_v2.py \\
    --mats ${MAGDIR}/${ET}.MATS.JC.txt --event-type $ET \\
    --b1 dcm_bams.txt --b2 nf_bams.txt --cov magnet_covariates.tsv \\
    --group1 DCM --group2 NF --min-samples 15 --exclude splitctl/nfB_1.ids \\
    --out ${JCDIR}/${ET}_DCMvsA_JC
  python phase2_diffsplice_v2.py \\
    --mats ${HCMDIR}/${ET}.MATS.JC.txt --event-type $ET \\
    --b1 hcm_bams.txt --b2 nf_bams.txt --cov hcm_nf_covariates.tsv \\
    --group1 HCM --group2 NF --min-samples 15 --exclude splitctl/nfA_1.ids \\
    --out ${JCDIR}/${ET}_HCMvsB_JC
'''

# ------------------------------------------------------------------- edit 5
A_ARM_C = "    --out jccheck/${COH}_${ET}_HCM_JC\ndone\n"
B_ARM_C = '''    --out ${JCDIR}/${COH}_${ET}_HCM_JC
done

# ---------------------------------------------------------------- control
# Part C reads ONLY the cohorts' own rMATS, which were never re-quantified,
# so at any read length these three files must be byte-identical to the 50 bp
# ones. Free, and a real check on phase-2 reproducibility.
if [ -n "$SUF" ]; then
  echo
  echo "=== CONTROL: part C must be byte-identical to the 50 bp run ==="
  _cfail=0
  for COH in Brigham Fuwai Stanford; do
    _new=${JCDIR}/${COH}_${ET}_HCM_JC.all_events.tsv
    _old=jccheck/${COH}_${ET}_HCM_JC.all_events.tsv
    if [ ! -s "$_old" ]; then
      echo "  [skip] $COH -- no 50 bp counterpart at $_old"
    elif cmp -s "$_new" "$_old"; then
      echo "  [ok]   $COH byte-identical"
    else
      echo "  [FAIL] $COH DIFFERS from the 50 bp run -- $_new vs $_old"
      _cfail=1
    fi
  done
  if [ "$_cfail" = "1" ]; then
    echo
    echo "  ⚠ The cohort arms depend on NOTHING that changed between read"
    echo "    lengths. If they differ, phase 2 is not reproducible and the"
    echo "    MAGNet arms from this run should not be used either."
  fi
fi
'''

# ------------------------------------------------------------------- edit 6
A_HD = ('echo "############### SUMMARY : JC vs JCEC ###############"\n'
        "python3 - <<'PY'\n")
B_HD = ('echo "############### SUMMARY : JC vs JCEC ###############"\n'
        '# <<\'PY\' is QUOTED, so bash does NOT substitute into the block.\n'
        '# SUF therefore has to arrive by environment, not interpolation.\n'
        'TBL_SUFFIX="$SUF" python3 - <<\'PY\'\n')

# ------------------------------------------------------------------- edit 7
A_PYSUF = 'ET = "SE"\n'
B_PYSUF = ('ET = "SE"\n'
           'SUF = os.environ.get("TBL_SUFFIX", "")   # passed in, not '
           'substituted\n')

# ------------------------------------------------------------------- edit 8
A_ARMS = '''ARMS = [
 ("MAGNet disjoint",
  f"splitctl/{ET}_DCMvsA_1.all_events.tsv",  f"jccheck/{ET}_DCMvsA_JC.all_events.tsv",
  f"splitctl/{ET}_HCMvsB_1.all_events.tsv",  f"jccheck/{ET}_HCMvsB_JC.all_events.tsv"),
]
for C in ("Brigham", "Fuwai", "Stanford"):
    ARMS.append((C,
        f"{ET}_v2filtered.all_events.tsv",        f"jccheck/{ET}_v2filtered_JC.all_events.tsv",
        f"{C}_{ET}_HCM.all_events.tsv",           f"jccheck/{C}_{ET}_HCM_JC.all_events.tsv"))
'''
B_ARMS = '''# MAGNet paths take SUF; the three cohorts' own tables never do.
ARMS = [
 ("MAGNet disjoint",
  f"splitctl/{ET}_DCMvsA_1{SUF}.all_events.tsv",  f"jccheck{SUF}/{ET}_DCMvsA_JC.all_events.tsv",
  f"splitctl/{ET}_HCMvsB_1{SUF}.all_events.tsv",  f"jccheck{SUF}/{ET}_HCMvsB_JC.all_events.tsv"),
]
for C in ("Brigham", "Fuwai", "Stanford"):
    ARMS.append((C,
        f"{ET}_v2filtered{SUF}.all_events.tsv",   f"jccheck{SUF}/{ET}_v2filtered_JC.all_events.tsv",
        f"{C}_{ET}_HCM.all_events.tsv",           f"jccheck{SUF}/{C}_{ET}_HCM_JC.all_events.tsv"))
'''

EDITS = [
    ("the memory request", A_HDR, B_HDR, 1),
    ("SUF, the suffixed dirs and the clobber guard", A_SETUP, B_SETUP, 1),
    ("arm A", A_ARM_A, B_ARM_A, 1),
    ("arm B", A_ARM_B, B_ARM_B, 1),
    ("arm C and the byte-identity control", A_ARM_C, B_ARM_C, 1),
    ("the quoted heredoc gets SUF by environment", A_HD, B_HD, 1),
    ("SUF inside the heredoc", A_PYSUF, B_PYSUF, 1),
    ("the comparison pairs", A_ARMS, B_ARMS, 1),
]


def apply_edits(text):
    if "SUF=${TBL_SUFFIX:-}" in text:
        return text, ["SUF is already defined — 20 is already patched "
                      "(restore the .bak to re-apply)"]
    problems, out = [], text
    for name, anchor, repl, want in EDITS:
        got = out.count(anchor)
        if got != want:
            problems.append(
                f"{name}: anchor occurs {got}x, expected {want}x\n"
                f"      anchor: {anchor.strip().splitlines()[0][:72]}")
            continue
        out = out.replace(anchor, repl, want)
    return out, problems


def verify(text):
    bad = []
    for w in ("SUF=${TBL_SUFFIX:-}",
              "MAGDIR=rmats_magnet${SUF}",
              "HCMDIR=rmats_HCMvsNF${SUF}",
              "JCDIR=jccheck${SUF}",
              "REFUSING — ${JCDIR}/ already contains files",
              "--mats ${MAGDIR}/${ET}.MATS.JC.txt",
              "--mats ${HCMDIR}/${ET}.MATS.JC.txt",
              "--out ${JCDIR}/${ET}_v2filtered_JC",
              "--out ${JCDIR}/${ET}_DCMvsA_JC",
              "--out ${JCDIR}/${ET}_HCMvsB_JC",
              "--out ${JCDIR}/${COH}_${ET}_HCM_JC",
              "CONTROL: part C must be byte-identical",
              'TBL_SUFFIX="$SUF" python3 - ',
              'SUF = os.environ.get("TBL_SUFFIX", "")',
              'f"splitctl/{ET}_DCMvsA_1{SUF}.all_events.tsv"',
              'f"{ET}_v2filtered{SUF}.all_events.tsv"',
              "#SBATCH --mem=8G"):
        if w not in text:
            bad.append(f"missing: {w[:58]}")
    for gone in ("--mats rmats_magnet/${ET}.MATS.JC.txt",
                 "--mats rmats_HCMvsNF/${ET}.MATS.JC.txt",
                 "--out jccheck/",
                 "mkdir -p jccheck\n",
                 "#SBATCH --mem=32G",
                 'f"jccheck/{ET}_DCMvsA_JC.all_events.tsv"',
                 'f"jccheck/{ET}_v2filtered_JC.all_events.tsv"'):
        if gone in text:
            bad.append(f"unsuffixed site survives: {gone[:50]}")
    # the COHORT inputs and tables must NOT have moved
    if "M=rmats_${COH}/${ET}.MATS.JC.txt" not in text:
        bad.append("the cohort rMATS path was altered — it must not be")
    if "rmats_${COH}${SUF}" in text:
        bad.append("a cohort rMATS dir was wrongly suffixed")
    if 'f"{C}_{ET}_HCM.all_events.tsv"' not in text:
        bad.append("a cohort table was altered — it must not be")
    if "{C}_{ET}_HCM{SUF}" in text:
        bad.append("a cohort table was wrongly suffixed")
    # SUF must be set before every use
    i_s = text.find("SUF=${TBL_SUFFIX:-}")
    for m in ("MAGDIR=rmats_magnet${SUF}", "--out ${JCDIR}/${ET}_v2filtered_JC",
              'TBL_SUFFIX="$SUF" python3 - '):
        j = text.find(m)
        if j < 0 or i_s > j:
            bad.append(f"SUF is not set above: {m[:44]}")
    return bad


def _bash_ok(text, tmpdir):
    p = os.path.join(tmpdir, "chk.sbatch")
    open(p, "w").write(text)
    r = subprocess.run(["bash", "-n", p], capture_output=True, text=True)
    return r.returncode == 0, r.stderr


def _heredoc_py(text):
    """Extract the quoted heredoc body so it can be compiled."""
    i = text.find("<<'PY'\n")
    if i < 0:
        return None
    body = text[i + len("<<'PY'\n"):]
    j = body.find("\nPY\n")
    return body if j < 0 else body[:j]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("target", nargs="?", default=TARGET)
    ap.add_argument("--apply", action="store_true")
    ap.add_argument("--self-test", action="store_true")
    a = ap.parse_args()
    if a.self_test:
        return self_test()

    self_test()
    print()
    print("=" * 78)
    print("142 — script 20 at 125 bp, so S5(a) can be built")
    print("=" * 78)
    if not os.path.exists(a.target):
        print(f"[fatal] {a.target} not found.", file=sys.stderr)
        return 2
    src = open(a.target, errors="replace").read()
    new, problems = apply_edits(src)
    if problems:
        print("REFUSING — anchors do not match:")
        for p in problems:
            print(f"  ✗ {p}")
        return 3
    bad = verify(new)
    if bad:
        print("REFUSING — patched text fails verification:")
        for b in bad:
            print(f"  ✗ {b}")
        return 3
    import tempfile
    ok, err = _bash_ok(new, tempfile.mkdtemp())
    if not ok:
        print("REFUSING — the patched script is not valid bash:")
        print(err)
        return 3

    print("".join(difflib.unified_diff(
        src.splitlines(True), new.splitlines(True),
        fromfile=a.target, tofile=a.target + " (patched)", n=2)))
    print("=" * 78)

    if not a.apply:
        print("DRY RUN — nothing written. Re-run with --apply, then:")
        print()
        print("  bash -n 20_jc_vs_jcec_readlength.sbatch && echo BASH_OK")
        print()
        print("  # 1. the guard must refuse to touch the existing 50 bp set")
        print("  sbatch 20_jc_vs_jcec_readlength.sbatch")
        print("  #   -> REFUSING, jccheck/ already contains files")
        print()
        print("  # 2. 125 bp, into a NEW directory")
        print("  TBL_SUFFIX=_rl125 sbatch 20_jc_vs_jcec_readlength.sbatch")
        print("  #   log: jc_vs_jcec_<jobid>.out, or %x_%j with 20's job name")
        print()
        print("  # 3. read the control block at the end of the log")
        print("  grep -A8 'CONTROL: part C' *_<jobid>.out")
        print("  sacct -j <jobid>.batch --format=JobID,MaxRSS,Elapsed")
        print()
        print("⚠ THE CONTROL IS THE POINT. Part C reads only the cohorts' own")
        print("  rMATS, which never moved, so its three files MUST come out")
        print("  byte-identical to jccheck/'s. If they do not, phase 2 is not")
        print("  reproducible and the MAGNet arms from the same run are not")
        print("  usable either.")
        print()
        print("⬜ THEN 62's part (a) still needs its JC/JCEC labels corrected")
        print("  and pointing at jccheck<SUF>/ — that is a separate patch, and")
        print("  (a) will cover THREE cohorts: Colorado has no JC table and")
        print("  script 20 does not make one.")
        return 0

    bak = f"{a.target}.bak.{time.strftime('%Y%m%d_%H%M%S')}"
    shutil.copy2(a.target, bak)
    open(a.target, "w").write(new)
    print(f"[write] {a.target}   (backup: {bak})")
    disk = open(a.target, errors="replace").read()
    bad = verify(disk)
    if bad or disk != new:
        print("✗ ON-DISK VERIFICATION FAILED — restoring:")
        for b in bad:
            print(f"    {b}")
        shutil.copy2(bak, a.target)
        return 3
    print("✓ verified on disk: both MAGNet rMATS dirs and every output follow")
    print("  SUF, the three cohort dirs and tables do not, the heredoc gets")
    print("  SUF by environment, an existing output directory is refused, and")
    print("  part C now checks itself against the 50 bp run.")
    return 0


def self_test():
    import tempfile

    fx = ("#!/bin/bash\n"
          "#SBATCH -p cpu\n"
          + A_HDR
          + "#SBATCH --qos=ngs16\n"
            "date\n"
          + A_SETUP
          + '\necho "=== A ==="\n'
          + A_ARM_A
          + '\necho "=== B ==="\n'
          + A_ARM_B
          + 'else\n  echo "  skipping"\nfi\n'
            '\necho "=== C ==="\n'
            "for COH in Brigham Fuwai Stanford; do\n"
            "  M=rmats_${COH}/${ET}.MATS.JC.txt\n"
            '  [ -s "$M" ] || { echo "  missing $M"; continue; }\n'
            "  python phase2_diffsplice_v2.py \\\n"
            '    --mats "$M" --event-type $ET \\\n'
            "    --min-samples ${MINS[$COH]} --no-covariates \\\n"
          + A_ARM_C
          + "\n"
          + A_HD
          + "import os\n"
            "import numpy as np, pandas as pd\n"
            "from scipy import stats\n"
            "\n"
          + A_PYSUF
          + 'PROGRAM = ["TTN"]\n'
            "\n"
          + A_ARMS
          + 'print("ARMS", ARMS)\n'
            "PY\n")

    new, probs = apply_edits(fx)
    assert not probs, probs
    assert not verify(new), verify(new)

    d = tempfile.mkdtemp()
    ok, err = _bash_ok(new, d)
    assert ok, err

    # the heredoc body must still be valid python
    body = _heredoc_py(new)
    assert body is not None, "heredoc not found"
    compile(body, "heredoc", "exec")

    # and it must resolve the right paths for each read length
    for suf, want, forbid in (
            ("", ["splitctl/SE_DCMvsA_1.all_events.tsv",
                  "jccheck/SE_v2filtered_JC.all_events.tsv",
                  "SE_v2filtered.all_events.tsv",
                  "Brigham_SE_HCM.all_events.tsv"],
             ["_rl125"]),
            ("_rl125", ["splitctl/SE_DCMvsA_1_rl125.all_events.tsv",
                        "jccheck_rl125/SE_v2filtered_JC.all_events.tsv",
                        "SE_v2filtered_rl125.all_events.tsv",
                        "Brigham_SE_HCM.all_events.tsv"],
             ["Brigham_SE_HCM_rl125.all_events.tsv",
              "jccheck/SE_v2filtered_JC"])):
        g = {"os": os}
        env = dict(os.environ)
        env.pop("TBL_SUFFIX", None)
        if suf:
            env["TBL_SUFFIX"] = suf
        src = body.replace("import numpy as np, pandas as pd\n", "")
        src = src.replace("from scipy import stats\n", "")
        src = src.replace('print("ARMS", ARMS)', "")
        old = os.environ.get("TBL_SUFFIX")
        if suf:
            os.environ["TBL_SUFFIX"] = suf
        else:
            os.environ.pop("TBL_SUFFIX", None)
        try:
            exec(compile(src, "hd", "exec"), g)
        finally:
            if old is None:
                os.environ.pop("TBL_SUFFIX", None)
            else:
                os.environ["TBL_SUFFIX"] = old
        flat = " ".join(str(x) for arm in g["ARMS"] for x in arm)
        for w in want:
            assert w in flat, (suf, w, flat)
        for f in forbid:
            assert f not in flat, (suf, f, flat)

    # the cohort rMATS dir and the clobber guard survive in the bash part
    assert "M=rmats_${COH}/${ET}.MATS.JC.txt" in new
    assert "rmats_${COH}${SUF}" not in new
    assert 'REFUSING — ${JCDIR}/ already contains files' in new
    assert "cmp -s \"$_new\" \"$_old\"" in new

    # double application refused, text unchanged
    n2, p2 = apply_edits(new)
    assert p2, "second application was not refused"
    assert n2 == new, "a refused re-application still altered the text"

    print("[ok] self-test passed: the patched script is valid bash and its "
          "heredoc is valid python; both MAGNet rMATS dirs and every output "
          "follow SUF while the three cohort dirs and tables do not; the "
          "heredoc resolves splitctl/_rl125, jccheck_rl125/ and "
          "SE_v2filtered_rl125 under TBL_SUFFIX while Brigham_SE_HCM stays "
          "unsuffixed, and resolves all of them unsuffixed without it; an "
          "existing output directory is refused; part C's byte-identity "
          "control is present; double application refused")
    return 0


if __name__ == "__main__":
    sys.exit(main())
