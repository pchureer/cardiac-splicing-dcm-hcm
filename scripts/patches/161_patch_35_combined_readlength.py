#!/usr/bin/env python3
"""
161_patch_35_combined_readlength.py -- the Stouffer half of [047]'s "51 vs 48".

WHY 35 AND NOT ONLY 34
----------------------
[047]: "An inverse-variance analysis was run as a sensitivity check and agrees
on the count of newly significant events (51 versus 48)."

  51 = `34_meta_analysis.py`  META A, NEW summed: 44+3+2+2   (inverse-variance)
  48 = `35_combined_evidence.py` POOL A, NEW summed: 34+6+2+6 (Stouffer)

160 measured the 51 side: **51 → 64** at 125 bp. This patch does the 48 side.

**The clause is an AGREEMENT claim.** At 50 bp the two methods differ by 3
events (6%) and "agrees" is fair. If the Stouffer side stays near 48 while the
inverse-variance side has moved to 64, the gap becomes 33% and the word
"agrees" is false — with both numbers correctly re-based. So this cannot be
edited by swapping digits; the claim has to be re-tested.

35 WAS INVISIBLE TO EVERY OBVIOUS GREP
--------------------------------------
`grep -ln "combined_a" *.py` returns only the **consumer**, `84_pool_b.py`.
`grep -rn "combined_a" *.sh *.sbatch` returns nothing. The producer builds the
name at runtime:

    out = f"combined_{name.split()[1].lower()}.tsv"      # line 182

That is Rule 10's corollary A — a computed output name is invisible to a
filename grep. The `unsuffixed_writes()` helper from 143 does not catch it
either: it looks for quoted literal filenames, and an f-string is not one.

THE FIVE EDITS
--------------
1. **`POOL_A`'s MAGNet arm only.** Colorado was never re-quantified; `POOL_B`
   is Brigham/Fuwai/Stanford. Same asymmetry as 129 (script 21) and 160
   (script 34).

2. **The suffix variable is `COMB_SUFFIX`, not `META_SUFFIX`.** Rule 9: the
   suffix parameter is not uniform. 34 and 35 are different producers;
   sharing one variable would re-base both from a single export, which is
   precisely the uniformity Rule 9 says does not hold here.

3. **The output name takes the suffix.** Without it a 125 bp run overwrites
   `combined_a.tsv` **in place** — the only copy of the Stouffer pool, dated
   28 Sep, read by `84_pool_b.py`. Rule 10 instance 1 all over again.

4. **Two guards**: a suffix that matches no table refuses (148), and a suffix
   set while `--out` is still `combined_evidence.txt` refuses.

5. **The NEW total is printed**, and the provenance line records
   `--dpsi-rule`. On disk `combined_evidence.txt` is 21:21 and
   `combined_a.tsv` is 21:26 — **five minutes apart, so from different
   invocations** — and 84 says Pool A uses the `discovery` rule while 35
   defaults to `any`. Every future run now states which rule produced it.

CALIBRATION — which rule reproduces REC_A?
------------------------------------------
`84_pool_b.py`'s record, Pool A, `discovery` rule:

    SE    events 18,111  fdr 4,797  effect 491  unan 409  new 34
    A3SS  events  4,072  fdr   531  effect  51  unan  45  new  6
    A5SS  events  2,635  fdr   459  effect  45  unan  37  new  2
    RI    events  2,999  fdr   712  effect  74  unan  59  new  6
    totals: 27,817 events, 550 unanimous, **NEW TOTAL 48**

⚠ `effect 51` for A3SS is **not** [047]'s 51, and at 125 bp 34's SE-class NEW
is **48**, which is not this 48 either. Three live collisions in one clause.

USAGE
-----
  python3 161_patch_35_combined_readlength.py            # dry run
  python3 161_patch_35_combined_readlength.py --apply
  python3 161_patch_35_combined_readlength.py --self-test
"""
import argparse
import difflib
import os
import shutil
import sys
import time

TARGET = "35_combined_evidence.py"

A_POOL = '''POOL_A = [("MAGNet", "{ET}_v2filtered"), ("Colorado", "{ET}_Colorado_DCMvsNF")]'''
B_POOL = '''# 161: only the MAGNet arm re-bases with the read length. Colorado was never
# re-quantified and POOL_B is the three HCM cohorts.
#
# RULE 9 — the suffix parameter is NOT uniform. This script's variable is
# COMB_SUFFIX, deliberately NOT the META_SUFFIX that 160 gave script 34.
# They are different producers writing different tables; one shared variable
# would re-base both from a single export and make it impossible to run one
# read length through 34 and another through 35 while diagnosing a
# disagreement between them — which is exactly what [047]'s claim requires.
CSUF = os.environ.get("COMB_SUFFIX", "")
POOL_A = [("MAGNet", "{ET}_v2filtered" + CSUF),
          ("Colorado", "{ET}_Colorado_DCMvsNF")]
# 161: the per-class NEW counts are printed but never summed, and [047]'s
# "48" IS that sum. Keyed by pool name so POOL A and POOL B cannot mix.
_NEW_TOTAL = {}'''

A_NEW = '''            new = f"{int((sig & ~dsig).sum())}"'''
B_NEW = '''            new = f"{int((sig & ~dsig).sum())}"
            _NEW_TOTAL[name] = (_NEW_TOTAL.get(name, 0)
                                + int((sig & ~dsig).sum()))'''

A_OUT = '''    A = pd.concat(allr, ignore_index=True)
    out = f"combined_{name.split()[1].lower()}.tsv"
    A.to_csv(out, sep="\\t", index=False, float_format="%.5g")'''
B_OUT = '''    A = pd.concat(allr, ignore_index=True)
    # 161: the sum [047] quotes. This script printed the per-class NEW counts
    # and stopped, so no grep of it or its output could find "48".
    if name in _NEW_TOTAL:
        P(f"  NEW TOTAL across classes: {_NEW_TOTAL[name]}"
          f"   <- the count [047] quotes for this pool")
    # 161: the suffix MUST reach the output name. Without it a 125 bp run
    # overwrites combined_a.tsv in place — the only copy of the Stouffer
    # pool, read by 84_pool_b.py. That is Rule 10, instance 1.
    out = f"combined_{name.split()[1].lower()}{CSUF}.tsv"
    A.to_csv(out, sep="\\t", index=False, float_format="%.5g")'''

A_GUARD = '''    ap.add_argument("--out", default="combined_evidence.txt")
    a = ap.parse_args()
    types = [t for t in a.types.split(",") if t in NKEEP]
    fh = open(a.out, "w")'''
B_GUARD = '''    ap.add_argument("--out", default="combined_evidence.txt")
    a = ap.parse_args()
    types = [t for t in a.types.split(",") if t in NKEEP]
    if CSUF:
        # 148: a suffix that re-bases nothing exits 0 and reports 50 bp
        # numbers under a 125 bp name.
        if not [p for et in types
                for p in (f"{et}_v2filtered{CSUF}_dedup.all_events.tsv",
                          f"{et}_v2filtered{CSUF}.all_events.tsv")
                if os.path.exists(p)]:
            raise SystemExit(
                f"REFUSING: COMB_SUFFIX={CSUF!r} but no "
                f"{{ET}}_v2filtered{CSUF}[_dedup].all_events.tsv exists for any "
                f"of {types}. The suffix would re-base nothing.")
        if a.out == ap.get_default("out"):
            raise SystemExit(
                f"REFUSING: COMB_SUFFIX is set but --out is still the default "
                f"{a.out!r}. Pass --out combined_evidence{CSUF}.txt, or the "
                f"28 Sep report is overwritten in place.")
    fh = open(a.out, "w")'''

A_PROV = '''    P("Stouffer weighted-Z on standardised statistics; weights sqrt(n).")'''
B_PROV = '''    P("Stouffer weighted-Z on standardised statistics; weights sqrt(n).")
    P(f"read length: MAGNet arm = {{ET}}_v2filtered{CSUF} "
      f"(COMB_SUFFIX={CSUF!r}); Colorado and all of POOL B are NOT "
      f"re-quantified")
    P(f"dpsi-rule: {a.dpsi_rule}   (84_pool_b.py's REC_A record is the "
      f"`discovery` rule; this script defaults to `any`, and the 28 Sep "
      f"combined_evidence.txt and combined_a.tsv are five minutes apart, so "
      f"they may be from different invocations)")'''

EDITS = [
    ("POOL_A's MAGNet arm", A_POOL, B_POOL, 1),
    ("the per-class NEW accumulator", A_NEW, B_NEW, 1),
    ("the output name and the NEW total", A_OUT, B_OUT, 1),
    ("the suffix guards", A_GUARD, B_GUARD, 1),
    ("the provenance lines", A_PROV, B_PROV, 1),
]


def apply_edits(text):
    if "CSUF = os.environ.get" in text:
        return text, ["35 already carries 161 (restore the .bak to re-apply)"]
    if "import os" not in text:
        return text, ["35 does not import os at module level; 161 needs it"]
    problems, out = [], text
    for name, anchor, repl, want in EDITS:
        got = out.count(anchor)
        if got != want:
            problems.append(f"{name}: anchor occurs {got}x, expected {want}x\n"
                            f"      first line: {anchor.splitlines()[0][:70]}")
            continue
        out = out.replace(anchor, repl, want)
    return out, problems


def verify(text):
    bad = []
    for w in ('CSUF = os.environ.get("COMB_SUFFIX", "")',
              '("MAGNet", "{ET}_v2filtered" + CSUF)',
              'f"combined_{name.split()[1].lower()}{CSUF}.tsv"',
              "_NEW_TOTAL = {}",
              "NEW TOTAL across classes",
              "REFUSING: COMB_SUFFIX",
              'ap.get_default("out")',
              "read length: MAGNet arm",
              "dpsi-rule:"):
        if w not in text:
            bad.append(f"missing: {w[:60]}")
    if 'out = f"combined_{name.split()[1].lower()}.tsv"' in text:
        bad.append("the output name is still un-suffixed; a 125 bp run would "
                   "overwrite combined_a.tsv in place")
    # Rule 9: this must NOT reach for 34's variable. The test is for a READ
    # of it, not for the word: 161's own comment explains why it does not use
    # META_SUFFIX, and a bare substring check rejected the patch it had just
    # written — the fourth time this week a check and its target drifted.
    for pat in ('os.environ.get("META_SUFFIX"', "os.environ.get('META_SUFFIX'",
                'os.environ["META_SUFFIX"]', "os.environ['META_SUFFIX']"):
        if pat in text:
            bad.append("35 reads META_SUFFIX; it must have its own "
                       "COMB_SUFFIX (Rule 9 — the suffix parameter is not "
                       "uniform)")
    if '"{ET}_Colorado_DCMvsNF" + CSUF' in text:
        bad.append("Colorado was suffixed; it was never re-quantified")
    for w in ('("Brigham", "Brigham_{ET}_HCM")',
              '("Fuwai", "Fuwai_{ET}_HCM")'):
        if w not in text:
            bad.append(f"POOL_B was altered: {w} no longer verbatim")
    return bad


FIXTURE = '''#!/usr/bin/env python3
"""35_combined_evidence.py -- Stouffer combined evidence. Phase 1."""
import argparse
import os

import numpy as np
import pandas as pd

NKEEP = {"SE": 3, "MXE": 5, "A3SS": 5, "A5SS": 5, "RI": 3}
POOL_A = [("MAGNet", "{ET}_v2filtered"), ("Colorado", "{ET}_Colorado_DCMvsNF")]
POOL_B = [("Brigham", "Brigham_{ET}_HCM"), ("Fuwai", "Fuwai_{ET}_HCM"),
          ("Stanford", "Stanford_{ET}_HCM")]
CAUTION = {("Fuwai", "MXE")}


def run(name, POOL, types, a, P):
    allr = []
    for et in types:
        sig = np.array([True, False])
        dsig = np.array([False, False])
        disc = POOL[0][0]      # as the real 35: the first-listed cohort
        new = ""
        if disc == "MAGNet":
            new = f"{int((sig & ~dsig).sum())}"
        allr.append(pd.DataFrame({"significant": sig}))
    if not allr:
        P("  nothing combined")
        return None
    A = pd.concat(allr, ignore_index=True)
    out = f"combined_{name.split()[1].lower()}.tsv"
    A.to_csv(out, sep="\\t", index=False, float_format="%.5g")
    P(f"  -> {out}   {len(A)} events, {int(A.significant.sum())} significant")
    return A


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--types", default="SE,MXE,A3SS,A5SS,RI")
    ap.add_argument("--fdr", type=float, default=0.05)
    ap.add_argument("--min-dpsi", type=float, default=0.10)
    ap.add_argument("--dpsi-rule", choices=["any", "discovery", "all"],
                    default="any")
    ap.add_argument("--no-unanimous", dest="unanimous", action="store_false")
    ap.add_argument("--out", default="combined_evidence.txt")
    a = ap.parse_args()
    types = [t for t in a.types.split(",") if t in NKEEP]
    fh = open(a.out, "w")

    def P(s=""):
        print(s)
        fh.write(s + "\\n")

    P("Stouffer weighted-Z on standardised statistics; weights sqrt(n).")
    run("POOL A -- the DCM signature", POOL_A, types, a, P)
    run("POOL B -- the HCM transfer", POOL_B, types, a, P)
'''


def self_test():
    import subprocess
    import tempfile
    fails = []

    new, probs = apply_edits(FIXTURE)
    if probs:
        fails.append(f"the fixture did not patch cleanly: {probs}")
    else:
        bad = verify(new)
        if bad:
            fails.append(f"verify() rejected its own output: {bad}")
        try:
            compile(new, "patched", "exec")
        except SyntaxError as e:
            fails.append(f"patched fixture does not compile: {e}")
        if 'lower()}.tsv"' in new:
            fails.append("the un-suffixed output name survived")
        if '("Colorado", "{ET}_Colorado_DCMvsNF")' not in new:
            fails.append("Colorado's entry was altered")
        if "Brigham_{ET}_HCM" not in new:
            fails.append("POOL_B was altered")
        again, probs2 = apply_edits(new)
        if not probs2 or "already carries 161" not in probs2[0]:
            fails.append("re-applying 161 did not refuse")

    # Rule 9: borrowing 34's variable must be caught
    borrowed = new.replace('os.environ.get("COMB_SUFFIX", "")',
                           'os.environ.get("META_SUFFIX", "")')
    if not any("Rule 9" in x for x in verify(borrowed)):
        fails.append("verify() did not catch 35 reading META_SUFFIX")
    # and the explanatory comment naming META_SUFFIX must NOT trip it
    if any("Rule 9" in x for x in verify(new)):
        fails.append("verify() tripped on its own comment mentioning "
                     "META_SUFFIX")
    if not any("overwrite combined_a.tsv" in x for x in
               verify(new.replace('lower()}{CSUF}.tsv"', 'lower()}.tsv"'))):
        fails.append("verify() missed an un-suffixed output name")

    noos = FIXTURE.replace("import os\n", "")
    _, p = apply_edits(noos)
    if not p or "does not import os" not in p[0]:
        fails.append("a file lacking `import os` was patched anyway")

    # run the patched main(): the NEW TOTAL must print, both pools must stay
    # separate, and the guards must fire
    with tempfile.TemporaryDirectory() as td:
        f = os.path.join(td, "fake35.py")
        open(f, "w").write(new + "\n\nif __name__ == '__main__':\n    main()\n")

        r = subprocess.run([sys.executable, f, "--out", "o.txt"], cwd=td,
                           env={k: v for k, v in os.environ.items()
                                if k != "COMB_SUFFIX"},
                           capture_output=True, text=True)
        o = r.stdout + r.stderr
        if "REFUSING" in o:
            fails.append(f"the 50 bp default run was refused: {o[:160]}")
        if o.count("NEW TOTAL across classes") != 1:
            fails.append(f"the NEW TOTAL printed "
                         f"{o.count('NEW TOTAL across classes')} times; POOL A "
                         f"should print it and POOL B should not (its first "
                         f"cohort is not MAGNet)")
        if "dpsi-rule: any" not in o:
            fails.append("the dpsi-rule provenance line did not print")

        env = dict(os.environ, COMB_SUFFIX="_rl125")
        r = subprocess.run([sys.executable, f, "--out", "o2.txt"], cwd=td,
                           env=env, capture_output=True, text=True)
        if "REFUSING: COMB_SUFFIX=" not in (r.stdout + r.stderr):
            fails.append("the no-such-table guard did not fire")
        open(os.path.join(td, "SE_v2filtered_rl125.all_events.tsv"), "w").write("x")
        r = subprocess.run([sys.executable, f], cwd=td, env=env,
                           capture_output=True, text=True)
        if "--out is still the default" not in (r.stdout + r.stderr):
            fails.append("the default --out guard did not fire")
        # clear the unsuffixed table the 50 bp run above left behind, or the
        # overwrite assertion below proves nothing
        for leftover in ("combined_a.tsv", "combined_b.tsv"):
            lp = os.path.join(td, leftover)
            if os.path.exists(lp):
                os.remove(lp)
        r = subprocess.run([sys.executable, f, "--out", "c_rl125.txt"],
                           cwd=td, env=env, capture_output=True, text=True)
        o = r.stdout + r.stderr
        if "REFUSING" in o:
            fails.append(f"the guards fired on a correct invocation: {o[:160]}")
        if not os.path.exists(os.path.join(td, "combined_a_rl125.tsv")):
            fails.append("the suffixed output table was not written")
        if os.path.exists(os.path.join(td, "combined_a.tsv")):
            fails.append("the SUFFIXED run wrote combined_a.tsv — the exact "
                         "overwrite this patch exists to prevent")

    print("161 self-test")
    for f in fails:
        print("  ✗ " + f)
    if not fails:
        print("  ✓ patches cleanly and compiles, output name takes the suffix,")
        print("    Colorado and POOL_B untouched, re-application refuses,")
        print("    Rule 9 enforced (COMB_SUFFIX, never META_SUFFIX), both")
        print("    guards fire, NEW TOTAL prints once for POOL A only, and a")
        print("    suffixed run leaves combined_a.tsv untouched")
    return 1 if fails else 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("target", nargs="?", default=TARGET)
    ap.add_argument("--apply", action="store_true")
    ap.add_argument("--self-test", action="store_true")
    a = ap.parse_args()
    if a.self_test:
        return self_test()
    if self_test():
        print("REFUSING: 161's own self-test failed.")
        return 3
    print()

    if not os.path.exists(a.target):
        print(f"REFUSING: {a.target} not found — run from the HEART directory")
        return 3
    src = open(a.target, errors="replace").read()
    new, probs = apply_edits(src)
    print("=" * 78)
    print("161 — script 35's read length, output names, and the missing total")
    print("=" * 78)
    if probs:
        print("REFUSING — anchors did not match:")
        for p in probs:
            print("    " + p)
        return 3
    print("".join(difflib.unified_diff(
        src.splitlines(True), new.splitlines(True),
        fromfile=a.target, tofile=a.target + " (patched)", n=2)))
    print("=" * 78)

    if not a.apply:
        print("DRY RUN — nothing written. Re-run with --apply, then:")
        print()
        print("  # 0. the 28 Sep record is already backed up in record_28sep/")
        print("  #    (combined_a.tsv md5 ddae24234a56868fcd588648210b627b)")
        print()
        print("  # 1. WHICH RULE made the record? 84 says `discovery`, 35")
        print("  #    defaults to `any`, and the .txt and .tsv on disk are")
        print("  #    five minutes apart. Try both, into scratch names.")
        print("  for R in any discovery; do")
        print("    python3 35_combined_evidence.py --dpsi-rule $R \\")
        print("        --out cal_$R.txt >/dev/null")
        print("    echo \"--- $R\"; grep -E 'POOL A|SE |A3SS|A5SS|RI |NEW TOTAL' \\")
        print("        cal_$R.txt | head -8")
        print("  done")
        print()
        print("  #    REC_A (discovery rule) expects:")
        print("  #      SE 18,111 / 4,797 / 491 / 409 / new 34")
        print("  #      A3SS 4,072/531/51/45/6   A5SS 2,635/459/45/37/2")
        print("  #      RI 2,999/712/74/59/6     NEW TOTAL 48")
        print()
        print("  # 2. the two refusals")
        print("  COMB_SUFFIX=_nope python3 35_combined_evidence.py --out /tmp/x")
        print("  COMB_SUFFIX=_rl125 python3 35_combined_evidence.py")
        print()
        print("  # 3. the 125 bp run, with the rule that calibrated")
        print("  COMB_SUFFIX=_rl125 python3 35_combined_evidence.py \\")
        print("      --dpsi-rule <the one that reproduced REC_A> \\")
        print("      --out combined_evidence_rl125.txt")
        print()
        print("  # 4. the record must have survived")
        print("  md5sum combined_a.tsv record_28sep/combined_a.tsv")
        print("  ls -l combined_a_rl125.tsv combined_b_rl125.tsv")
        print()
        print("⚠ Step 1 rewrites combined_a.tsv (the unsuffixed output name IS")
        print("  that file, so no guard can protect it) — which is why step 0")
        print("  exists. Compare md5s against record_28sep/ afterwards.")
        print()
        print("⚠ THE CLAIM, NOT THE DIGITS. 34 moved 51 -> 64. If Stouffer's")
        print("  48 does not move comparably, \"agrees\" becomes false and the")
        print("  clause needs rewriting, not renumbering.")
        return 0

    bak = f"{a.target}.bak.{time.strftime('%Y%m%d_%H%M%S')}"
    shutil.copy2(a.target, bak)
    open(a.target, "w").write(new)
    print(f"[write] {a.target}   (backup: {bak})")
    disk = open(a.target, errors="replace").read()
    bad = verify(disk)
    try:
        compile(disk, a.target, "exec")
    except SyntaxError as e:
        bad.append(f"does not compile: {e}")
    if bad or disk != new:
        print("✗ ON-DISK VERIFICATION FAILED — restoring:")
        for b in bad:
            print(f"    {b}")
        shutil.copy2(bak, a.target)
        return 3
    print("✓ verified on disk: POOL_A's MAGNet arm follows COMB_SUFFIX, the")
    print("  output names carry it, Colorado and POOL B are untouched, the")
    print("  guards are in place, the NEW total prints, and the provenance")
    print("  records the dpsi-rule.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
