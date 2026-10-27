#!/usr/bin/env python3
"""
111_rl125_consistency_check.py -- find the 50/125 bp mixes that
110_apply_rl125_patches.py left behind.

THE DEFECT
----------
110 renamed every `--out PREFIX` to `PREFIX_rl125`. It did NOT look at whether
the same script later READS `PREFIX.all_events.tsv`. Where it does, the script
now writes 125 bp tables and reads 50 bp ones.

Confirmed in 08_composition_adjusted_discovery.sbatch, around line 70:

    --out ${ET}_v2compadj_rl125        <- 110 renamed this
    --out ${ET}_v2fulladj_rl125        <- and this
    ...
    full=${ET}_v2fulladj.all_events.tsv    <- but NOT this read

The compositionCheck step would compare the 125 bp run against the 50 bp
tables. It would not error. It would print a plausible number.

This script finds every such read. It does NOT auto-fix, because a reference to
a 50 bp table is sometimes DELIBERATE -- a script that compares the new result
against the frozen baseline is correct to name the old file. Only a human can
tell those apart, so each hit is printed with its line and a classification.

USAGE
-----
  python 111_rl125_consistency_check.py
  python 111_rl125_consistency_check.py --fix 08_composition_adjusted_discovery.sbatch
  python 111_rl125_consistency_check.py --self-test

--fix renames reads of prefixes THIS SCRIPT WROTE in the same file, which is
the unambiguous case. Everything else stays manual.
"""
import argparse
import datetime
import difflib
import os
import re
import shutil
import sys

SUFFIX = "_rl125"

# the prefixes 110 renamed, per script -- must match 110's PLAN
RENAMED = {
    "04_phase2_v2_checks.sbatch":
        ["SE_v2legacy", "${ET}_v2filtered", "SE_v2pooladj", "SE_v2age25"],
    "08_composition_adjusted_discovery.sbatch":
        ["${ET}_v2compadj", "${ET}_v2fulladj", "${ET}_compositionCheck"],
    "36_matched_pair_panel.sbatch":
        ["mp_${ET}_selection", "mp_${ET}_magnet_matched",
         "mp_${ET}_lbnl_matched", "mp_${ET}_magnet_other"],
    "20_jc_vs_jcec_readlength.sbatch":
        ["jccheck/${ET}_v2filtered_JC", "jccheck/${ET}_DCMvsA_JC",
         "jccheck/${ET}_HCMvsB_JC", "jccheck/${COH}_${ET}_HCM_JC"],
    # these two take --out "$OUT"; the OUT= assignment is the prefix
    "72_run_cmfrac.sbatch": ["SE_v2cmfrac"],
    "75_cmfrac_allclasses.sbatch": ["${ET}_v2cmfrac"],
    "23_signature_rbp_scan.sbatch": [],
    "56_rbm20_targets45.sbatch": [],
    "21_replication_v2.sbatch": [],
}

# prefixes that are legitimately 50 bp comparators ANYWHERE: a read of these
# is a baseline, not a bug. Keep this list short and justified.
FROZEN_BASELINES = {
    "SE_v2filtered.all_events.tsv": "the frozen 50 bp discovery comparator",
    "SE_HCM_v2.all_events.tsv": "the frozen 50 bp HCM comparator",
    "SE_DCMvsHCM.all_events.tsv": "the frozen 50 bp null comparator",
}

# every table prefix in the project that the 125 bp work renames. A read of one
# of these in ANY patched script is suspect even if that script did not write it.
ALL_RENAMED_PREFIXES = [
    "_v2compadj", "_v2fulladj", "_v2cmfrac", "_v2pooladj", "_v2age25",
    "_v2legacy", "_v2filtered", "_compositionCheck",
    "mp_", "jccheck/",
]


def hits_in(path, prefixes):
    """Occurrences of PREFIX not already followed by _rl125, on code lines."""
    out = []
    for i, ln in enumerate(open(path, errors="replace").read().splitlines(), 1):
        if ln.lstrip().startswith("#"):
            continue
        for p in prefixes:
            for m in re.finditer(re.escape(p) + r"(?!" + SUFFIX + ")", ln):
                out.append((i, p, ln.rstrip()))
    return out


def classify(prefix, line, wrote_here):
    for base, why in FROZEN_BASELINES.items():
        if base in line:
            return "BASELINE", f"deliberate 50 bp comparator ({why})"
    if "--out" in line and prefix in line.split("--out", 1)[1]:
        return "OUTPUT", "an --out that 110 should have renamed; re-run 110"
    if wrote_here:
        return "MIX", ("this script WRITES this prefix at 125 bp and READS it "
                       "at 50 bp -- silent mix, must be renamed")
    return "REVIEW", ("read of a table that the 125 bp work renames, written "
                      "elsewhere -- confirm which read length is intended")


def do_fix(path, prefixes, report):
    text = open(path, errors="replace").read()
    lines = text.splitlines(keepends=True)
    new = []
    n = 0
    for ln in lines:
        if ln.lstrip().startswith("#"):
            new.append(ln)
            continue
        cur = ln
        skip = any(b in cur for b in FROZEN_BASELINES)
        if not skip:
            for p in prefixes:
                pat = re.escape(p) + r"(?!" + SUFFIX + ")"
                cur, k = re.subn(pat, p + SUFFIX, cur)
                n += k
        new.append(cur)
    out = "".join(new)
    if out == text:
        report.append("    nothing to fix")
        return
    diff = difflib.unified_diff(text.splitlines(), out.splitlines(),
                                fromfile=path + " (mixed)",
                                tofile=path + " (consistent)",
                                lineterm="", n=1)
    report.append("    --- diff ---")
    for d in diff:
        report.append(f"    {d}")
    stamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    bak = f"{path}.bak111_{stamp}"
    shutil.copy2(path, bak)
    open(path, "w").write(out)
    chk = open(path, errors="replace").read()
    left = hits_in(path, prefixes)
    left = [h for h in left
            if not any(b in h[2] for b in FROZEN_BASELINES)]
    if left:
        report.append(f"    ✗ POST-CHECK: {len(left)} occurrence(s) remain")
        for i, p, ln in left[:6]:
            report.append(f"        line {i}: {ln}")
        report.append(f"    restore with: cp {bak} {path}")
    else:
        report.append(f"    [ok] {n} rename(s) applied; backup {bak}")
    if SUFFIX + SUFFIX in chk:
        report.append("    ✗ POST-CHECK: double suffix _rl125_rl125 present")
        report.append(f"    restore with: cp {bak} {path}")


def self_test():
    import tempfile
    d = tempfile.mkdtemp()
    p = os.path.join(d, "08_composition_adjusted_discovery.sbatch")
    open(p, "w").write(
        "# a comment naming ${ET}_v2fulladj.all_events.tsv must be ignored\n"
        "python phase2.py --out ${ET}_v2compadj_rl125\n"
        "python phase2.py --out ${ET}_v2fulladj_rl125\n"
        "python cmp.py nine=${ET}_v2compadj.all_events.tsv \\\n"
        "              full=${ET}_v2fulladj.all_events.tsv \\\n"
        "    --base SE_v2filtered.all_events.tsv \\\n"
        "    --out ${ET}_compositionCheck_rl125\n")
    pre = RENAMED["08_composition_adjusted_discovery.sbatch"]
    h = hits_in(p, pre)
    # the two reads must be found; the comment must not
    found = {(i, pp) for i, pp, _ in h}
    assert any(pp == "${ET}_v2compadj" for _, pp in found), h
    assert any(pp == "${ET}_v2fulladj" for _, pp in found), h
    assert not any(i == 1 for i, _ in found), "comment line was flagged"
    # the frozen baseline line must classify as BASELINE
    kinds = {classify(pp, ln, True)[0] for i, pp, ln in h}
    assert "MIX" in kinds, kinds
    base_line = [ln for ln in open(p) if "SE_v2filtered" in ln][0]
    assert classify("_v2filtered", base_line, False)[0] == "BASELINE"

    rep = []
    do_fix(p, pre, rep)
    txt = open(p).read()
    assert "nine=${ET}_v2compadj_rl125.all_events.tsv" in txt, txt
    assert "full=${ET}_v2fulladj_rl125.all_events.tsv" in txt, txt
    assert "--base SE_v2filtered.all_events.tsv" in txt, \
        "the frozen baseline was renamed!"
    assert "_rl125_rl125" not in txt, "double suffix"
    assert "# a comment naming ${ET}_v2fulladj.all_events.tsv" in txt, \
        "comment was edited"
    # idempotent
    rep2 = []
    do_fix(p, pre, rep2)
    assert "_rl125_rl125" not in open(p).read()
    print("[ok] self-test passed: the two 50 bp reads found and fixed, the "
          "frozen SE_v2filtered baseline left alone, comments untouched, "
          "idempotent, no double suffix")
    return 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fix", nargs="*", metavar="SCRIPT",
                    help="rename same-file reads of prefixes that script writes")
    ap.add_argument("--self-test", action="store_true")
    a = ap.parse_args()
    if a.self_test:
        return self_test()

    rep = ["=" * 78, "rl125 CONSISTENCY CHECK -- 50/125 bp mixes left by 110",
           "=" * 78, ""]
    any_mix = False
    for s, wrote in RENAMED.items():
        if not os.path.exists(s):
            continue
        rep.append(f"═══ {s}")
        # the prefixes this script wrote, plus every renamed prefix anywhere
        h_own = hits_in(s, wrote)
        h_any = hits_in(s, ALL_RENAMED_PREFIXES)
        seen, rows = set(), []
        for i, p, ln in h_own + h_any:
            if (i, ln) in seen:
                continue
            seen.add((i, ln))
            kind, why = classify(p, ln, p in wrote or
                                 any(w.endswith(p) for w in wrote))
            rows.append((kind, i, p, ln, why))
        if not rows:
            rep.append("    clean")
        order = {"MIX": 0, "OUTPUT": 1, "REVIEW": 2, "BASELINE": 3}
        for kind, i, p, ln, why in sorted(rows, key=lambda r: (order[r[0]], r[1])):
            mark = {"MIX": "✗ MIX     ", "OUTPUT": "✗ OUTPUT  ",
                    "REVIEW": "? REVIEW  ", "BASELINE": "✓ BASELINE"}[kind]
            rep.append(f"    {mark} line {i:>4}  {ln.strip()[:84]}")
            rep.append(f"               {why}")
            if kind == "MIX":
                any_mix = True
        if [r for r in rows if r[0] == "MIX"]:
            rep.append(f"    -> fix with: python {os.path.basename(__file__)} "
                       f"--fix {s}")
        rep.append("")

        if a.fix is not None and (not a.fix or s in a.fix):
            if [r for r in rows if r[0] == "MIX"]:
                rep.append(f"    FIXING {s}")
                do_fix(s, wrote, rep)
                rep.append("")

    rep += ["=" * 78,
            "✗ MIX      the script writes this prefix at 125 bp and reads it at",
            "           50 bp. Produces a wrong number WITHOUT erroring. Must fix.",
            "✗ OUTPUT   an --out 110 missed. Re-run 110.",
            "? REVIEW   a read of a renamed table written by another script.",
            "           Decide which read length is intended and say so in a",
            "           comment on that line.",
            "✓ BASELINE a deliberate 50 bp comparator. Correct as written.",
            "",
            "DO NOT SUBMIT ANY SCRIPT WITH A MIX.",
            ]
    print("\n".join(rep))
    return 1 if any_mix and a.fix is None else 0


if __name__ == "__main__":
    sys.exit(main())
