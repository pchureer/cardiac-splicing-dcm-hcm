#!/usr/bin/env python3
"""
168_patch_87_record_rl125.py -- 87 carries a THIRD frozen copy of the 50 bp
numbers, and refuses correctly when the inputs disagree with it.

After 167 wrote correct 125 bp inputs, 87 refused:

    REFUSING — G3-G5 — the master file does not reproduce the record:
      - MAGNet DCMvsNF SE: (53456, 813) vs record (62971, 926)
      ...

That is the gate working. 87 holds REC_SE, REC_ALL and REC_PERCLASS_MAGNET --
an independent statement of what Table 3's numbers must be -- and will not
build a table whose inputs contradict it.

⚠ THE RECORD MUST NOT BE DERIVED FROM THE FILES IT CHECKS.
Every other patch in this series removed a hard-coded constant. This one adds
some, deliberately. A cross-check computed from the same master it validates
passes by construction and tests nothing. The 125 bp values below are typed
from the project's own independent sources -- the manuscript's own counts and
the ten gates 167 measured against -- not read back out of t3_rl125/.

COLORADO IS INHERITED, NOT RETYPED. The 125 bp record is built by unpacking
the 50 bp one and overriding only the MAGNet keys, so "Colorado was never
re-quantified" is stated in the code rather than relied on in transcription.

THE TRIGGER IS THE DATA, NOT A FLAG. If "rl125" appears anywhere on the command
line and T3_RECORD is not set to rl125, 87 exits before G1. A gate that can be
bypassed by changing a path is not a gate -- pointing 87 at a 125 bp directory
must REQUIRE a 125 bp record, never skip the check.

NOT ADDED: a ("MAGNet", "HCMvsNF") entry in REC_ALL. There is none at 50 bp,
and 108,478/670 was first measured by 167 today. Putting it in the record would
only echo that measurement back at itself.

USAGE
-----
  python3 168_patch_87_record_rl125.py [--apply] [--self-test]

  T3_RECORD=rl125 python3 87_build_table3.py --dir t3_rl125 --out Table3_rl125
"""
import argparse
import difflib
import os
import shutil
import sys
import time

TARGET = "87_build_table3.py"

A = '''REC_PERCLASS_MAGNET = {"SE": (62971, 926), "MXE": (62762, 231), "A3SS": (15292, 90),
                       "A5SS": (8988, 84), "RI": (6102, 114)}'''

B = '''REC_PERCLASS_MAGNET = {"SE": (62971, 926), "MXE": (62762, 231), "A3SS": (15292, 90),
                       "A5SS": (8988, 84), "RI": (6102, 114)}

# ---- 168: the same record at 125 bp --------------------------------------
# Typed from the manuscript's own counts, NOT read back out of the master
# these dicts validate. A cross-check derived from the file it checks passes
# by construction and tests nothing.
#
# Colorado is inherited by unpacking, because that cohort was never
# re-quantified; only the MAGNet keys are overridden.
REC_SE_RL125 = {**REC_SE,
                ("MAGNet", "DCMvsNF"):  (53456, 813),
                ("MAGNet", "HCMvsNF"):  (39607, 431),
                ("MAGNet", "DCMvsHCM"): (38273, 0)}
REC_ALL_RL125 = {**REC_ALL,
                 ("MAGNet", "DCMvsNF"): (144942, 1294)}
REC_PERCLASS_MAGNET_RL125 = {"SE": (53456, 813), "MXE": (62537, 229),
                             "A3SS": (14346, 87), "A5SS": (8360, 74),
                             "RI": (6243, 91)}

_T3REC = os.environ.get("T3_RECORD", "").strip().lower()
if "rl125" in " ".join(sys.argv) and _T3REC != "rl125":
    sys.exit(
        "[fatal] a 125 bp path was named on the command line, but T3_RECORD is\\n"
        "        not set to rl125. G3-G5 would check 125 bp inputs against the\\n"
        "        50 bp record and refuse, or worse, pass for the wrong reason.\\n"
        "        Re-run as:\\n"
        "          T3_RECORD=rl125 python3 87_build_table3.py \\\\\\n"
        "              --dir t3_rl125 --out Table3_rl125\\n"
        "        Nothing built.")
if _T3REC == "rl125":
    if "rl125" not in " ".join(sys.argv):
        sys.exit("[fatal] T3_RECORD=rl125 but no 125 bp path was named. The 125 bp\\n"
                 "        record would be checked against 50 bp inputs. Nothing built.")
    REC_SE = REC_SE_RL125
    REC_ALL = REC_ALL_RL125
    REC_PERCLASS_MAGNET = REC_PERCLASS_MAGNET_RL125
    print("[168] record: 125 bp "
          "(MAGNet re-based; Colorado inherited unchanged)")
elif _T3REC:
    sys.exit(f"[fatal] T3_RECORD={_T3REC!r} is not a record this script knows. "
             "Use rl125 or leave it unset.")'''


def apply_edits(text):
    if "168:" in text:
        return text, ["87 already carries 168 (restore the .bak to re-apply)"]
    n = text.count(A)
    if n != 1:
        return text, [f"REC_PERCLASS_MAGNET anchor occurs {n}x, expected 1x"]
    return text.replace(A, B, 1), []


def verify(text):
    bad = []
    for w in ('REC_SE_RL125 = {**REC_SE,', '("MAGNet", "DCMvsNF"):  (53456, 813)',
              '("MAGNet", "HCMvsNF"):  (39607, 431)',
              '("MAGNet", "DCMvsHCM"): (38273, 0)',
              '("MAGNet", "DCMvsNF"): (144942, 1294)',
              '_T3REC = os.environ.get("T3_RECORD", "").strip().lower()',
              'T3_RECORD is\\n'):
        if w not in text:
            bad.append(f"missing: {w[:50]}")
    # the 50 bp record must survive as the default
    for w in ('("MAGNet", "DCMvsNF"):    (62971, 926)',
              '("MAGNet", "DCMvsNF"):    (156115, 1445)',
              '"MXE": (62762, 231)'):
        if w not in text:
            bad.append(f"the 50 bp record was destroyed: {w[:46]}")
    # the gate itself must still be there
    if "REC_PERCLASS_MAGNET" not in text or "REC_ALL" not in text:
        bad.append("a record dict was lost")
    return bad


FIX = '''import os
import sys

REC_SE = {
    ("MAGNet", "DCMvsNF"):    (62971, 926),
    ("MAGNet", "HCMvsNF"):    (48024, 468),
    ("MAGNet", "DCMvsHCM"):   (45380, 0),
    ("Colorado", "DCMvsNF"):  (19953, 476),
    ("Colorado", "ICMvsNF"):  (17019, 407),
    ("Colorado", "DCMvsICM"): (19998, 2),
}
REC_ALL = {
    ("MAGNet", "DCMvsNF"):    (156115, 1445),
    ("Colorado", "DCMvsNF"):  (31569, 846),
    ("Colorado", "ICMvsNF"):  (27231, 863),
    ("Colorado", "DCMvsICM"): (32381, 2),
}
REC_PERCLASS_MAGNET = {"SE": (62971, 926), "MXE": (62762, 231), "A3SS": (15292, 90),
                       "A5SS": (8988, 84), "RI": (6102, 114)}
'''


def self_test():
    f = []
    new, probs = apply_edits(FIX)
    if probs:
        f.append(f"fixture did not patch: {probs}")
    else:
        bad = verify(new)
        if bad:
            f.append(f"verify rejected its own output: {bad}")
        try:
            compile(new, "x", "exec")
        except SyntaxError as e:
            f.append(f"does not compile: {e}")
        else:
            def run(env, argv):
                os.environ.pop("T3_RECORD", None)
                if env:
                    os.environ["T3_RECORD"] = env
                old = sys.argv
                sys.argv = argv
                ns = {}
                try:
                    exec(compile(new, "x", "exec"), ns)
                    return ns, None
                except SystemExit as e:
                    return None, str(e)
                finally:
                    sys.argv = old
                    os.environ.pop("T3_RECORD", None)

            # 1. the 50 bp default is untouched
            ns, err = run(None, ["87_build_table3.py"])
            if err:
                f.append(f"50 bp default refused: {err[:60]}")
            elif ns["REC_SE"][("MAGNet", "DCMvsNF")] != (62971, 926):
                f.append("50 bp default did not keep (62971, 926)")
            elif ns["REC_PERCLASS_MAGNET"]["RI"] != (6102, 114):
                f.append("50 bp per-class default was changed")

            # 2. a 125 bp path WITHOUT the record: must refuse
            ns, err = run(None, ["87_build_table3.py", "--dir", "t3_rl125"])
            if ns is not None or not err or "T3_RECORD is" not in err:
                f.append("125 bp path without T3_RECORD did not refuse")

            # 3. the record WITHOUT a 125 bp path: must also refuse
            ns, err = run("rl125", ["87_build_table3.py"])
            if ns is not None or not err or "no 125 bp path" not in err:
                f.append("T3_RECORD without a 125 bp path did not refuse")

            # 4. an unknown record value: must refuse
            ns, err = run("rl50", ["87_build_table3.py", "--dir", "t3_rl125"])
            if ns is not None or not err:
                f.append("an unknown T3_RECORD value did not refuse")

            # 5. both together: MAGNet re-based, Colorado inherited unchanged
            ns, err = run("rl125", ["87_build_table3.py", "--dir", "t3_rl125",
                                    "--out", "Table3_rl125"])
            if err:
                f.append(f"the correct invocation refused: {err[:70]}")
            else:
                S, AL, P = (ns["REC_SE"], ns["REC_ALL"],
                            ns["REC_PERCLASS_MAGNET"])
                for k, want in ((("MAGNet", "DCMvsNF"), (53456, 813)),
                                (("MAGNet", "HCMvsNF"), (39607, 431)),
                                (("MAGNet", "DCMvsHCM"), (38273, 0))):
                    if S[k] != want:
                        f.append(f"REC_SE{k} = {S[k]}, want {want}")
                if AL[("MAGNet", "DCMvsNF")] != (144942, 1294):
                    f.append("REC_ALL MAGNet DCMvsNF wrong")
                for k in (("Colorado", "DCMvsNF"), ("Colorado", "ICMvsNF"),
                          ("Colorado", "DCMvsICM")):
                    if S[k] != {("Colorado", "DCMvsNF"): (19953, 476),
                                ("Colorado", "ICMvsNF"): (17019, 407),
                                ("Colorado", "DCMvsICM"): (19998, 2)}[k]:
                        f.append(f"Colorado {k} was altered: {S[k]}")
                    if k in AL and AL[k] != {("Colorado", "DCMvsNF"): (31569, 846),
                                             ("Colorado", "ICMvsNF"): (27231, 863),
                                             ("Colorado", "DCMvsICM"): (32381, 2)}[k]:
                        f.append(f"REC_ALL Colorado {k} was altered")
                if sum(v[0] for v in P.values()) != 144942:
                    f.append("the 125 bp per-class tested do not sum to 144,942")
                if sum(v[1] for v in P.values()) != 1294:
                    f.append("the 125 bp per-class significant do not sum to 1,294")
                if P["SE"] != S[("MAGNet", "DCMvsNF")]:
                    f.append("per-class SE disagrees with REC_SE")
                if ("MAGNet", "HCMvsNF") in AL:
                    f.append("REC_ALL gained an HCM entry it should not have")

        again, p2 = apply_edits(new)
        if not p2 or "already carries 168" not in p2[0]:
            f.append("re-application did not refuse")

    print("168 self-test")
    for x in f:
        print("  ✗ " + x)
    if not f:
        print("  ✓ the 50 bp record is untouched by default; a 125 bp path without")
        print("    T3_RECORD refuses, as does T3_RECORD without a 125 bp path and")
        print("    an unknown record value; with both, MAGNet is re-based and")
        print("    Colorado is byte-for-byte inherited; the 125 bp per-class record")
        print("    sums to 144,942 / 1,294 and agrees with REC_SE; REC_ALL gains no")
        print("    HCM entry; re-application refuses")
    return 1 if f else 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("target", nargs="?", default=TARGET)
    ap.add_argument("--apply", action="store_true")
    ap.add_argument("--self-test", action="store_true")
    a = ap.parse_args()
    if a.self_test:
        return self_test()
    if self_test():
        print("\nREFUSING: 168's own self-test failed.")
        return 3
    print()
    if not os.path.exists(a.target):
        print(f"REFUSING: {a.target} not found")
        return 3
    src = open(a.target, errors="replace").read()
    new, probs = apply_edits(src)
    print("=" * 74)
    print("168 — 87's frozen record, at 125 bp")
    print("=" * 74)
    if probs:
        print("REFUSING:")
        for p in probs:
            print("    " + p)
        return 3
    print("".join(difflib.unified_diff(src.splitlines(True), new.splitlines(True),
                                       fromfile=a.target,
                                       tofile=a.target + " (patched)", n=2)))
    if not a.apply:
        print("DRY RUN. Re-run with --apply, then:")
        print()
        print("  # 50 bp control — must still build the published table")
        print("  python3 87_build_table3.py --out /tmp/T3cal")
        print()
        print("  # the guard — must REFUSE, nothing built")
        print("  python3 87_build_table3.py --dir t3_rl125 --out Table3_rl125")
        print()
        print("  # 125 bp")
        print("  T3_RECORD=rl125 python3 87_build_table3.py \\")
        print("      --dir t3_rl125 --out Table3_rl125")
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
    print("✓ verified: the 50 bp record is the default, the 125 bp record is")
    print("  selected only with T3_RECORD=rl125, and G3-G5 still runs either way.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
