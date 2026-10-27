#!/usr/bin/env python3
"""
160_patch_34_meta_readlength.py -- A8b: [047]'s "51 versus 48" at 125 bp.

WHAT 51 IS, AND WHY NO GREP FOUND IT
------------------------------------
[047]: "An inverse-variance analysis was run as a sensitivity check and agrees
on the count of newly significant events (51 versus 48)."

  51 = `34_meta_analysis.py`, META A, the NEW column summed over classes:
       SE 44 + A3SS 3 + A5SS 2 + RI 2 = 51        (inverse-variance, DL random)
  48 = Pool A, Stouffer, `84_pool_b.py`'s REC_A `new` fields:
       SE 34 + A3SS 6 + A5SS 2 + RI 6 = 48

`meta_analysis.txt` (28 Sep) IS the published run — `--model random`, which is
inverse-variance weighted with a between-study term. Grepping it for `51`
returned nothing because **51 is a sum the script never computes.** It prints
the four per-class NEW counts and stops. A published number can be an
arithmetic step that exists only in the manuscript; grep finds the inputs, not
the total. This patch makes the total a printed quantity so it is auditable
next time.

THE FOUR EDITS
--------------
1. **`META_A`'s MAGNet arm re-bases; nothing else does.**
   `META_A = [("MAGNet", "{ET}_v2filtered"), ("Colorado", ...)]` — Colorado was
   never re-quantified and `META_B` is all HCM cohorts, so only the MAGNet
   entry takes the suffix. The same asymmetry 129 applied to script 21.

2. **THE OUTPUT NAMES. This is the dangerous one.** Line 189 is
   `out = f"meta_{name.split()[1].lower()}.tsv"` — hard-coded. A 125 bp run
   would overwrite `meta_a.tsv` and `meta_b.tsv` **in place**, and
   `84_pool_b.py` reads `meta_a.tsv`. That is how `symmetry_MXE.tsv` was lost:
   the computation was re-runnable, the original bytes were not. The suffix now
   goes into the output name.

3. **A refusal when the suffix re-bases nothing**, and a refusal when
   `--out` is left at its default while a suffix is set (which would overwrite
   the 50 bp report). 148's lesson: a suffix that re-bases nothing exits 0 and
   publishes a 50 bp answer under a 125 bp name.

4. **The NEW TOTAL is printed**, recomputed from the concatenated frame — so
   the figure the manuscript quotes comes out of the script that owns it.

CALIBRATION — the 50 bp run must reproduce ALL of this
------------------------------------------------------
META A: SE 18,111 / sig 235 / NEW 44 / lost 280
        A3SS 4,072 / 26 / 3 / 24      A5SS 2,635 / 20 / 2 / 25
        RI 2,999 / 29 / 2 / 44        meta_a 27,817 events, 310 significant
        **NEW TOTAL 51**
META B: SE 11,464 / 93   MXE 6,048 / 14   A3SS 2,753 / 9
        A5SS 1,757 / 1   RI 1,727 / 9     meta_b 23,749 events, 126 significant

USAGE
-----
  python3 160_patch_34_meta_readlength.py            # dry run
  python3 160_patch_34_meta_readlength.py --apply
  python3 160_patch_34_meta_readlength.py --self-test
"""
import argparse
import difflib
import os
import shutil
import sys
import time

TARGET = "34_meta_analysis.py"

A_META = '''META_A = [("MAGNet", "{ET}_v2filtered"), ("Colorado", "{ET}_Colorado_DCMvsNF")]'''
B_META = '''# 160: only the MAGNet arm re-bases with the read length. Colorado was never
# re-quantified, and META_B is Brigham/Fuwai/Stanford — all 50/101/150 bp HCM
# cohorts, untouched by the MAGNet re-quantification. Same asymmetry that
# 129_patch_21_replication.py applied to script 21: discovery suffixed,
# replication cohorts not.
MSUF = os.environ.get("META_SUFFIX", "")
META_A = [("MAGNet", "{ET}_v2filtered" + MSUF),
          ("Colorado", "{ET}_Colorado_DCMvsNF")]'''

A_OUT = '''    A = pd.concat(allrows, ignore_index=True)
    out = f"meta_{name.split()[1].lower()}.tsv"
    A.to_csv(out, sep="\\t", index=False, float_format="%.5g")'''
B_OUT = '''    A = pd.concat(allrows, ignore_index=True)
    # 160: the NEW TOTAL is what [047] quotes as "51", and this script never
    # computed it — it printed the four per-class counts and stopped, so no
    # grep of this file or its output could find the published number.
    # Recomputed from the concatenated frame, which is identical to summing
    # the per-class counts because both predicates are per-row.
    if "FDR_MAGNet" in A.columns:
        _dsig = ((A.FDR_MAGNet < a.fdr)
                 & (A.dPSI_MAGNet.abs() >= a.min_dpsi))
        P(f"  NEW TOTAL across classes: {int((A.significant & ~_dsig).sum())}"
          f"   <- the count [047] quotes")
    # 160: the suffix MUST reach the output name. Without it a 125 bp run
    # overwrites meta_a.tsv in place, which 84_pool_b.py reads.
    out = f"meta_{name.split()[1].lower()}{MSUF}.tsv"
    A.to_csv(out, sep="\\t", index=False, float_format="%.5g")'''

A_GUARD = '''    types = [t for t in a.types.split(",") if t in NKEEP]
    fh = open(a.out, "w")'''
B_GUARD = '''    types = [t for t in a.types.split(",") if t in NKEEP]
    if MSUF:
        # 148: a suffix that re-bases nothing exits 0 and publishes a 50 bp
        # answer under a 125 bp name.
        if not [p for et in types
                for p in (f"{et}_v2filtered{MSUF}_dedup.all_events.tsv",
                          f"{et}_v2filtered{MSUF}.all_events.tsv")
                if os.path.exists(p)]:
            raise SystemExit(
                f"REFUSING: META_SUFFIX={MSUF!r} but no "
                f"{{ET}}_v2filtered{MSUF}[_dedup].all_events.tsv exists for any "
                f"of {types}. The suffix would re-base nothing and the run "
                f"would report 50 bp numbers under a 125 bp name.")
        if a.out == ap.get_default("out"):
            raise SystemExit(
                f"REFUSING: META_SUFFIX is set but --out is still the default "
                f"{a.out!r}. Pass --out meta_analysis{MSUF}.txt, or the 50 bp "
                f"report is overwritten in place.")
    fh = open(a.out, "w")'''

A_PROV = '''    P("exon-level keys throughout, matching the locked convention")'''
B_PROV = '''    P("exon-level keys throughout, matching the locked convention")
    P(f"read length: MAGNet arm = {{ET}}_v2filtered{MSUF} "
      f"(META_SUFFIX={MSUF!r}); Colorado and all of META B are NOT "
      f"re-quantified")'''

EDITS = [
    ("META_A's MAGNet arm", A_META, B_META, 1),
    ("the output name and the NEW total", A_OUT, B_OUT, 1),
    ("the suffix guards", A_GUARD, B_GUARD, 1),
    ("the read-length provenance line", A_PROV, B_PROV, 1),
]


def apply_edits(text):
    if "MSUF = os.environ.get" in text:
        return text, ["34 already carries 160 (restore the .bak to re-apply)"]
    if "import os" not in text:
        return text, ["34 does not import os at module level; 160 needs it "
                      "for environ and path checks — add the import first"]
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
    for w in ('MSUF = os.environ.get("META_SUFFIX", "")',
              '("MAGNet", "{ET}_v2filtered" + MSUF)',
              'f"meta_{name.split()[1].lower()}{MSUF}.tsv"',
              "NEW TOTAL across classes",
              "REFUSING: META_SUFFIX",
              'ap.get_default("out")',
              "read length: MAGNet arm"):
        if w not in text:
            bad.append(f"missing: {w[:60]}")
    # the un-suffixed output name must be GONE — this is edit 2's whole point
    if 'out = f"meta_{name.split()[1].lower()}.tsv"' in text:
        bad.append("the output name is still un-suffixed; a 125 bp run would "
                   "overwrite meta_a.tsv in place")
    # Colorado and META_B must NOT have been suffixed
    if '"{ET}_Colorado_DCMvsNF" + MSUF' in text or \
            '"{ET}_Colorado_DCMvsNF" +MSUF' in text:
        bad.append("Colorado was suffixed; it was never re-quantified")
    for w in ('("Brigham", "Brigham_{ET}_HCM")',
              '("Fuwai", "Fuwai_{ET}_HCM")'):
        if w not in text:
            bad.append(f"META_B was altered: {w} no longer present verbatim")
    # the 50 bp default must survive
    if 'os.environ.get("META_SUFFIX", "")' not in text:
        bad.append("META_SUFFIX has no empty default; 50 bp would break")
    return bad


FIXTURE = '''#!/usr/bin/env python3
"""34_meta_analysis.py -- the two meta-analyses. Phase 1."""
import argparse
import os

import numpy as np
import pandas as pd
from scipy import stats

NKEEP = {"SE": 3, "MXE": 5, "A3SS": 5, "A5SS": 5, "RI": 3}
META_A = [("MAGNet", "{ET}_v2filtered"), ("Colorado", "{ET}_Colorado_DCMvsNF")]
META_B = [("Brigham", "Brigham_{ET}_HCM"), ("Fuwai", "Fuwai_{ET}_HCM"),
          ("Stanford", "Stanford_{ET}_HCM")]
EXCLUDE = {("Colorado", "MXE")}


def run(name, META, types, a, P):
    allrows = []
    for et in types:
        pass
    if not allrows:
        P("  nothing pooled")
        return None
    A = pd.concat(allrows, ignore_index=True)
    out = f"meta_{name.split()[1].lower()}.tsv"
    A.to_csv(out, sep="\\t", index=False, float_format="%.5g")
    P(f"  -> {out}   ({len(A)} events, {int(A.significant.sum())} significant)")
    return A


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--types", default="SE,MXE,A3SS,A5SS,RI")
    ap.add_argument("--fdr", type=float, default=0.05)
    ap.add_argument("--min-dpsi", type=float, default=0.10)
    ap.add_argument("--model", choices=["random", "fixed"], default="random")
    ap.add_argument("--out", default="meta_analysis.txt")
    a = ap.parse_args()
    types = [t for t in a.types.split(",") if t in NKEEP]
    fh = open(a.out, "w")

    def P(s=""):
        print(s)
        fh.write(s + "\\n")

    P(f"model: {a.model}-effects (DerSimonian-Laird)   "
      f"significance: FDR<{a.fdr} and |pooled dPSI|>={a.min_dpsi}")
    P("exon-level keys throughout, matching the locked convention")
    A = run("META A -- the DCM signature", META_A, types, a, P)
    B = run("META B -- the HCM transfer", META_B, types, a, P)
'''


def self_test():
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

        # the output name must now carry the suffix, and the bare form be gone
        if 'lower()}{MSUF}.tsv' not in new:
            fails.append("the output name did not take the suffix")
        if 'lower()}.tsv"' in new:
            fails.append("the un-suffixed output name survived")

        # Colorado and META_B untouched
        if '("Colorado", "{ET}_Colorado_DCMvsNF")' not in new:
            fails.append("Colorado's entry was altered")
        if 'Brigham_{ET}_HCM' not in new:
            fails.append("META_B was altered")

        # idempotence: a second application must refuse, not double-patch
        again, probs2 = apply_edits(new)
        if not probs2 or "already carries 160" not in probs2[0]:
            fails.append("re-applying 160 did not refuse")

    # a file without `import os` must refuse rather than emit broken code
    noos = FIXTURE.replace("import os\n", "")
    _, p = apply_edits(noos)
    if not p or "does not import os" not in p[0]:
        fails.append("a file lacking `import os` was patched anyway")

    # verify() must catch the two regressions it exists for
    v = verify(new.replace('f"meta_{name.split()[1].lower()}{MSUF}.tsv"',
                           'f"meta_{name.split()[1].lower()}.tsv"'))
    if not any("overwrite meta_a.tsv" in x for x in v):
        fails.append("verify() missed an un-suffixed output name")
    v = verify(new.replace('("Colorado", "{ET}_Colorado_DCMvsNF")',
                           '("Colorado", "{ET}_Colorado_DCMvsNF" + MSUF)'))
    if not any("Colorado was suffixed" in x for x in v):
        fails.append("verify() missed a suffixed Colorado arm")

    # the guards must actually fire: exercise the patched main() in a temp dir
    # with META_SUFFIX set and no suffixed table present
    with tempfile.TemporaryDirectory() as td:
        p = os.path.join(td, "fake34.py")
        open(p, "w").write(new + "\n\nif __name__ == '__main__':\n    main()\n")
        import subprocess
        env = dict(os.environ, META_SUFFIX="_rl125")
        r = subprocess.run([sys.executable, p], cwd=td, env=env,
                           capture_output=True, text=True)
        if "REFUSING: META_SUFFIX=" not in (r.stderr + r.stdout):
            fails.append(f"the no-such-table guard did not fire: "
                         f"{(r.stderr or r.stdout)[:160]}")
        # now create a suffixed table, so only the --out guard should fire
        open(os.path.join(td, "SE_v2filtered_rl125.all_events.tsv"), "w").write("x")
        r = subprocess.run([sys.executable, p], cwd=td, env=env,
                           capture_output=True, text=True)
        if "--out is still the default" not in (r.stderr + r.stdout):
            fails.append(f"the default --out guard did not fire: "
                         f"{(r.stderr or r.stdout)[:160]}")
        # and with both supplied it must get past the guards
        r = subprocess.run([sys.executable, p, "--out", "meta_rl125.txt"],
                           cwd=td, env=env, capture_output=True, text=True)
        if "REFUSING" in (r.stderr + r.stdout):
            fails.append(f"the guards fired on a correct invocation: "
                         f"{(r.stderr or r.stdout)[:160]}")
        # 50 bp (no suffix) must need no arguments at all
        r = subprocess.run([sys.executable, p], cwd=td,
                           env={k: v for k, v in os.environ.items()
                                if k != "META_SUFFIX"},
                           capture_output=True, text=True)
        if "REFUSING" in (r.stderr + r.stdout):
            fails.append(f"the 50 bp default run was refused: "
                         f"{(r.stderr or r.stdout)[:160]}")

    print("160 self-test")
    for f in fails:
        print("  ✗ " + f)
    if not fails:
        print("  ✓ patches cleanly and compiles, output name takes the suffix")
        print("    and the bare form is gone, Colorado and META_B untouched,")
        print("    re-application refuses, a file without `import os` refuses,")
        print("    verify() catches both regressions, and all three guards")
        print("    were exercised by running the patched main()")
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
        print("REFUSING: 160's own self-test failed.")
        return 3
    print()

    if not os.path.exists(a.target):
        print(f"REFUSING: {a.target} not found — run from the HEART directory")
        return 3
    src = open(a.target, errors="replace").read()
    new, probs = apply_edits(src)
    print("=" * 78)
    print("160 — script 34's read length, output names, and the missing total")
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
        print("  # 1. 50 bp must reproduce the 28 Sep record EXACTLY.")
        print("  #    Note the distinct --out: meta_analysis.txt is the")
        print("  #    published file and must not be touched.")
        print("  python3 34_meta_analysis.py --out meta_analysis_cal.txt")
        print("  diff <(grep -E 'META|events|->' meta_analysis.txt) \\")
        print("       <(grep -E 'META|events|->' meta_analysis_cal.txt)")
        print()
        print("  #    expected, and NEW TOTAL must print 51:")
        print("  #      SE 18,111 / 235 / NEW 44 / lost 280")
        print("  #      A3SS 4,072/26/3/24   A5SS 2,635/20/2/25")
        print("  #      RI 2,999/29/2/44     meta_a 27,817 / 310")
        print()
        print("  # 2. the two refusals")
        print("  META_SUFFIX=_nosuchthing python3 34_meta_analysis.py \\")
        print("      --out /tmp/x.txt          # must REFUSE: re-bases nothing")
        print("  META_SUFFIX=_rl125 python3 34_meta_analysis.py")
        print("      # must REFUSE: --out still default")
        print()
        print("  # 3. the 125 bp run")
        print("  META_SUFFIX=_rl125 python3 34_meta_analysis.py \\")
        print("      --out meta_analysis_rl125.txt")
        print()
        print("⚠ Check that meta_a.tsv and meta_b.tsv are UNCHANGED after")
        print("  step 3 and that meta_a_rl125.tsv / meta_b_rl125.tsv appeared.")
        print("  84_pool_b.py reads meta_a.tsv; overwriting it would corrupt")
        print("  the Pool A cross-check with 125 bp content under a 50 bp name.")
        print()
        print("⚠ META B is unchanged by design — Brigham/Fuwai/Stanford were")
        print("  never re-quantified. If its numbers move, the patch is wrong.")
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
    print("✓ verified on disk: the MAGNet arm follows META_SUFFIX, the output")
    print("  names carry it, Colorado and META B are untouched, the suffix")
    print("  guards are in place, and the NEW TOTAL is printed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
