#!/usr/bin/env python3
"""
114_scan_called_scripts.py -- find the 50/125 bp mixes that live in the PYTHON
scripts the patched SBATCH files call.

THE HOLE THIS CLOSES
--------------------
`111_rl125_consistency_check.py` scanned only the .sbatch files. But an sbatch
can hand the real work to a python script that hard-codes its own filenames:

    36_matched_pair_panel.sbatch line 180:
        TYPES="$TYPES" python3 36_matched_pair_analysis.py

`36_matched_pair_analysis.py` takes no file arguments, so every table it reads
is hard-coded. The sbatch was patched to WRITE `mp_${ET}_*_rl125`; the analysis
script still READ `mp_${ET}_*`. Job 407779 therefore reported the 50 bp answer
(platform 0.584, 95% CI 0.478-0.693; people 0.988, 0.889-1.089; n=203) --
identical to the recorded 50 bp values -- while printing "sig ... 952" from the
new run beside "selection: 956 events" from the old tables. It did not error.

Two classes, and they are not equally urgent:

  CRITICAL   a python script INVOKED BY a patched sbatch, reading a prefix that
             the same sbatch now writes as _rl125. Silent mix inside one job.
             This is what happened to 36.

  TRACK 2    any other script reading a renamed prefix. It will FAIL LOUDLY on
             a missing file when run later. Safe direction; handle with the
             rest of the gates.

A script that receives its filenames as ARGUMENTS from the sbatch (e.g.
`06_sensitivity_concordance.py --primary "$BASE"`) is fine -- the sbatch
already carries the right name. Only hard-coded names matter, so this script
reports whether each hit sits in a script that also takes file arguments.

USAGE
-----
  python 114_scan_called_scripts.py
  python 114_scan_called_scripts.py --sbatch 36_matched_pair_panel.sbatch
  python 114_scan_called_scripts.py --self-test
"""
import argparse
import glob
import os
import re
import sys

SUFFIX = "_rl125"

# sbatch files that have been repointed to 125 bp (from 110's PLAN), plus 11,
# which was found later to be the A3SS/A5SS/RI producer.
PATCHED = [
    "04_phase2_v2_checks.sbatch",
    "08_composition_adjusted_discovery.sbatch",
    "72_run_cmfrac.sbatch",
    "75_cmfrac_allclasses.sbatch",
    "23_signature_rbp_scan.sbatch",
    "36_matched_pair_panel.sbatch",
]
# known to need the same treatment but NOT yet patched
PENDING = ["11_other_event_types.sbatch", "20_jc_vs_jcec_readlength.sbatch"]

# table-name stems the 125 bp work renames. A read of one of these without
# _rl125 is suspect.
PREFIXES = [
    "v2filtered", "v2compadj", "v2fulladj", "v2cmfrac", "v2pooladj",
    "v2age25", "v2legacy", "compositionCheck", "specificity", "lostwhy",
    "mp_", "jccheck/", "HCM_v2", "DCMvsHCM",
]

# these are legitimately 50 bp wherever they appear: other cohorts, never
# re-quantified.
COHORT_SAFE = ["Colorado", "Brigham", "Fuwai", "Stanford", "LBNL"]

PYCALL = re.compile(r"python[23]?\s+(?:-\S+\s+)*([A-Za-z0-9_./-]+\.py)")
FILEARG = re.compile(r"--(?:mats|primary|compare|base|adj|disc|out|cov|tab)\b")


def called_scripts(sbatch):
    """-> {py_path: True if the sbatch passes it file arguments}"""
    if not os.path.exists(sbatch):
        return {}
    out = {}
    lines = open(sbatch, errors="replace").read().splitlines()
    for i, ln in enumerate(lines):
        if ln.lstrip().startswith("#"):
            continue
        # `echo "    python 76_build_s1_v2.py"` prints INSTRUCTIONS for a human.
        # It is not an invocation. Treat echo/printf lines as text.
        if re.match(r"^\s*(echo|printf)\b", ln):
            continue
        for m in PYCALL.finditer(ln):
            py = m.group(1)
            if "phase2_diffsplice" in py:
                continue  # takes everything as arguments
            # look at this line and the next few (backslash continuations)
            blob = ln
            j = i
            while blob.rstrip().endswith("\\") and j + 1 < len(lines):
                j += 1
                blob += " " + lines[j]
            out[py] = bool(FILEARG.search(blob))
    return out


def hits(py, prefixes):
    res = []
    if not os.path.exists(py):
        return None
    for i, ln in enumerate(open(py, errors="replace").read().splitlines(), 1):
        s = ln.strip()
        if s.startswith("#"):
            continue
        for p in prefixes:
            for _m in re.finditer(re.escape(p) + r"(?!" + SUFFIX + ")", ln):
                safe = any(c in ln for c in COHORT_SAFE)
                res.append((i, p, s[:92], safe))
                break
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sbatch", nargs="*")
    ap.add_argument("--self-test", action="store_true")
    a = ap.parse_args()
    if a.self_test:
        return self_test()

    targets = a.sbatch or PATCHED
    print("=" * 80)
    print("114 — python scripts called by the patched SBATCH files")
    print("=" * 80)
    print()

    critical = []
    for sb in targets:
        calls = called_scripts(sb)
        print(f"═══ {sb}")
        if not os.path.exists(sb):
            print("    not found")
            print()
            continue
        if not calls:
            print("    calls no python script other than phase2_diffsplice_v2.py")
            print()
            continue
        for py, takes_args in sorted(calls.items()):
            h = hits(py, PREFIXES)
            if h is None:
                print(f"    {py}: called but NOT FOUND on disk")
                continue
            live = [x for x in h if not x[3]]
            tag = "takes file arguments" if takes_args else \
                  "NO file arguments — every name is hard-coded"
            print(f"    {py}   ({tag})")
            if not live:
                print("        no unrenamed table name — clean")
                continue
            for i, p, txt, _safe in live[:14]:
                mark = "✗ CRITICAL" if not takes_args else "? check"
                print(f"        {mark} line {i:>4}  {txt}")
                if not takes_args:
                    critical.append((sb, py, i, txt))
            if len(live) > 14:
                print(f"        ... and {len(live) - 14} more")
            safe = [x for x in h if x[3]]
            if safe:
                print(f"        ({len(safe)} further hit(s) name another cohort "
                      f"— Colorado/Brigham/Fuwai/Stanford/LBNL — and are "
                      f"correctly 50 bp)")
        print()

    print("=" * 80)
    if critical:
        print(f"✗ {len(critical)} CRITICAL hit(s): a hard-coded 50 bp table name")
        print("  inside a script that a patched sbatch runs in the same job.")
        print("  The job writes 125 bp and reads 50 bp, and does NOT error.")
        print()
        print("  For each: rename the read to _rl125, OR give the script the")
        print("  filename as an argument from the sbatch so there is one source")
        print("  of truth. The second is better but is a code change.")
        print()
        print("  ⚠ ANY RESULT ALREADY PRODUCED BY THESE JOBS IS THE 50 bp ANSWER.")
        print("    Re-run after fixing. Job 407779's decomposition is the known case.")
    else:
        print("✅ no hard-coded unrenamed table name in any called python script")
    print()
    print("NOT scanned, and they still need the same treatment:")
    for p in PENDING:
        print(f"  {p}")
    print()
    print("⚠ 11_other_event_types.sbatch WRITES THE FROZEN COMPARATORS:")
    print("    line  67  --out ${ET}_HCM_v2      -> SE_HCM_v2.all_events.tsv")
    print("    line 104  --out ${ET}_DCMvsHCM    -> SE_DCMvsHCM.all_events.tsv")
    print("  Those two files are the calibration targets for 106, 107 and 113")
    print("  (70,003 / 48,024 / 468 and 45,380 / 3 / 1 / 0). Re-running 11")
    print("  without renaming EVERY output destroys them silently.")
    return 1 if critical else 0


def self_test():
    import tempfile
    d = tempfile.mkdtemp()
    cwd = os.getcwd()
    os.chdir(d)
    try:
        open("x_panel.sbatch", "w").write(
            "#!/bin/bash\n"
            "# python3 commented_out.py should be ignored\n"
            'python phase2_diffsplice_v2.py --mats a --out mp_${ET}_sel_rl125\n'
            'TYPES="$T" python3 x_analysis.py\n'
            'python 06_sensitivity_concordance.py --primary "$BASE" \\\n'
            '    --compare composition=${ET}_v2compadj_rl125.all_events.tsv\n')
        open("x_analysis.py", "w").write(
            "# mp_SE_selection.all_events.tsv in a comment is ignored\n"
            'A = "mp_{}_selection.all_events.tsv".format(t)\n'
            'B = f"mp_{t}_magnet_matched.all_events.tsv"\n'
            'C = "Colorado_SE_v2filtered.all_events.tsv"\n')
        open("06_sensitivity_concordance.py", "w").write(
            "import sys\n# reads whatever --primary gives it\n")

        calls = called_scripts("x_panel.sbatch")
        assert "x_analysis.py" in calls, calls
        assert calls["x_analysis.py"] is False, "should have no file args"
        assert calls["06_sensitivity_concordance.py"] is True, calls
        assert "commented_out.py" not in calls, "comment was parsed"
        assert "phase2_diffsplice_v2.py" not in calls

        h = hits("x_analysis.py", PREFIXES)
        live = [x for x in h if not x[3]]
        assert any(x[0] == 2 for x in live), h
        assert any(x[0] == 3 for x in live), h
        assert not any(x[0] == 1 for x in live), "comment flagged"
        safe = [x for x in h if x[3]]
        assert any("Colorado" in x[2] for x in safe), h
        print("[ok] self-test passed: python call parsed from the sbatch, "
              "comments ignored, phase2 excluded, hard-coded vs argument-taking "
              "distinguished, .format() and f-string reads both found, "
              "another cohort's table marked correctly 50 bp")
        return 0
    finally:
        os.chdir(cwd)


if __name__ == "__main__":
    sys.exit(main())
