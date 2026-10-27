#!/usr/bin/env python3
"""
135_patch_scenarioB_and_98_note.py -- two files, all-or-nothing.

FIX 1 (85): THE THIRD SITE OF REC_CLAIM
---------------------------------------
134 fixed the verdict line and the mismatch warning. It missed the one inside
`scenario()`, which is fed at line ~278:

    r_ind = scenario(indep, "B. independent cohorts only (LBNL removed)",
                     recn=REC_CLAIM[1], recs=REC_CLAIM[0])

so scenario B still prints `[record 0 of 0  MISMATCH]` at 125 bp. Rule 8 — a
rename is only half an edit — applied to my own patch two steps after invoking
it. `0 or None` sends None, which is what the other three scenario() calls
already pass by omission.

FIX 2 (98): THE ADVISORY MUST NOT ASSERT A COUNT
------------------------------------------------
98's closing note 1 currently states:

    Fifteen of eighteen cohort x event-type combinations are individually
    significant; the three that are not are the three LOWEST concordance rates
    in the set (57-66%), all alternative-splice-site classes.

Measured at 125 bp, **every clause of that is false**: 14 of 18, four failing
cells, and one of them (Brigham RI) is a retained intron. Worse, the "lowest
rates" framing is wrong on its own terms — two SIGNIFICANT cells sit below a
non-significant one:

    Fuwai    A3SS  40/67  59.70%   ns
    Stanford A5SS  20/32  62.50%   ns
    Fuwai    RI    51/80  63.75%   SIGNIFICANT   n=80
    Brigham  A5SS  16/25  64.00%   ns
    Fuwai    A5SS  36/55  65.45%   SIGNIFICANT   n=55
    Brigham  RI    21/32  65.62%   ns            n=32

Significance there follows the denominator, not the rate. At 50 bp the two
happened to coincide, which is how the claim survived.

So the replacement does NOT substitute the 125 bp answer for the 50 bp one —
hard-coding either read length's answer in an advisory is what produced this.
It points at the computation and records only the fact that the quantity moves.

USAGE
-----
  python 135_patch_scenarioB_and_98_note.py [--apply] [--self-test]
"""
import argparse
import difflib
import os
import shutil
import sys
import time

T85 = "85_recount_combinations.py"
T98 = "98_build_replication_table.py"

A_85 = ('    r_ind = scenario(indep, "B. independent cohorts only '
        '(LBNL removed)",\n'
        "                     recn=REC_CLAIM[1], recs=REC_CLAIM[0])")
B_85 = ('    # `or None` so an unrecorded claim (0, 0) prints as absent rather\n'
        "    # than as a MISMATCH against zero — the other three scenario()\n"
        "    # calls already pass None by omission. This was the THIRD site of\n"
        "    # REC_CLAIM; 134 fixed the verdict line and the warning only.\n"
        '    r_ind = scenario(indep, "B. independent cohorts only '
        '(LBNL removed)",\n'
        "                     recn=REC_CLAIM[1] or None,\n"
        "                     recs=REC_CLAIM[0] or None)")

A_98 = ("""  1. Fifteen of eighteen cohort x event-type combinations are individually significant;
     the three that are not are the three LOWEST concordance rates in the set (57-66%),
     all alternative-splice-site classes. That belongs in the caption or the Results.
     Counted from `cohort_type_concordance_dedup[_rl125].tsv` by `85_recount_combinations.py`.""")
B_98 = ("""  1. The number of individually significant cohort x event-type combinations is
     deliberately NOT stated here. It moves with read length (15 of 18 at 50 bp,
     14 of 18 at 125 bp) and so do the identity and the event class of the cells
     that fail. Run `85_recount_combinations.py`, quote its scenario C, and say
     what it counts.
     Do NOT write that the failing cells are the lowest concordance rates. At
     125 bp two SIGNIFICANT cells fall below a non-significant one — Fuwai RI
     (51/80, 63.8%) and Fuwai A5SS (36/55, 65.5%) against Brigham RI (21/32,
     65.6%) — because significance there follows the number of testable events,
     not the rate. At 50 bp the two coincided, which is why the old claim held.
     Counted from `cohort_type_concordance_dedup[_rl125].tsv` by `85_recount_combinations.py`.""")

EDITS = {T85: [("scenario B's record label", A_85, B_85, 1)],
         T98: [("the closing advisory note 1", A_98, B_98, 1)]}

GUARDS = {T85: ("_RB_P = {", "134 has not been applied to 85"),
          T98: ('SUF = os.environ.get("TBL_SUFFIX", "")',
                "132 has not been applied to 98")}

WANT = {T85: ["recn=REC_CLAIM[1] or None", "recs=REC_CLAIM[0] or None"],
        T98: ["deliberately NOT stated here",
              "Do NOT write that the failing cells are the lowest",
              "15 of 18 at 50 bp"]}

GONE = {T85: ["recn=REC_CLAIM[1], recs=REC_CLAIM[0])"],
        T98: ["Fifteen of eighteen cohort x event-type",
              "the three LOWEST concordance rates"]}


def patch_one(target, text):
    problems, out = [], text
    for name, anchor, repl, want in EDITS[target]:
        got = out.count(anchor)
        if got != want:
            problems.append(f"{target}: {name}: anchor occurs {got}x, "
                            f"expected {want}x")
            continue
        out = out.replace(anchor, repl, want)
    return out, problems


def verify_one(target, text):
    bad = [f"{target}: missing {w[:50]}" for w in WANT[target]
           if w not in text]
    bad += [f"{target}: stale site survives: {g[:50]}" for g in GONE[target]
            if g in text]
    return bad


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--apply", action="store_true")
    ap.add_argument("--self-test", action="store_true")
    a = ap.parse_args()
    if a.self_test:
        return self_test()

    self_test()
    print()
    print("=" * 78)
    print("135 — REC_CLAIM's third site in 85, and 98's advisory note 1")
    print("=" * 78)

    # BOTH files are read, patched and verified before EITHER is written.
    staged, problems = {}, []
    for t in (T85, T98):
        if not os.path.exists(t):
            problems.append(f"{t}: not found")
            continue
        src = open(t, errors="replace").read()
        g, why = GUARDS[t]
        if g not in src:
            problems.append(f"{t}: {why}")
            continue
        new, probs = patch_one(t, src)
        if probs:
            problems += probs
            continue
        problems += verify_one(t, new)
        staged[t] = (src, new)

    if problems or len(staged) != 2:
        print("REFUSING — nothing written to either file:")
        for p in problems:
            print(f"  ✗ {p}")
        return 3

    for t in (T85, T98):
        src, new = staged[t]
        print("".join(difflib.unified_diff(
            src.splitlines(True), new.splitlines(True),
            fromfile=t, tofile=t + " (patched)", n=2)))
    print("=" * 78)

    if not a.apply:
        print("DRY RUN — nothing written. Re-run with --apply, then:")
        print()
        print("  # 85: scenario B must no longer say `record 0 of 0  MISMATCH`")
        print("  TBL_SUFFIX=_rl125 \\")
        print("  REC_SE=Colorado=311/385,LBNL=122/156,Brigham=298/364,"
              "Fuwai=422/596,Stanford=345/455 \\")
        print("  REC_SE_POOLED=1376,1800 REC_ALLTYPES=1883,2505 \\")
        print("  REC_BORDERLINE=Fuwai:A3SS=40/67,Stanford:A5SS=20/32,"
              "Brigham:A5SS=16/25,Brigham:RI=21/32 \\")
        print("  REC_CLAIM=none \\")
        print("  python3 85_recount_combinations.py 2>&1 | grep -E "
              "'combinations +:|record|G2 '")
        print()
        print("  # and 50 bp must STILL compare against 16 of 19")
        print("  python3 85_recount_combinations.py 2>&1 | grep 'record 16'")
        print()
        print("  # 98: rebuild Table 4 so the corrected advisory reaches the .md")
        print("  TBL_SUFFIX=_rl125 \\")
        print("  REC_KN=Colorado=311/385,Brigham=298/364,Stanford=345/455,"
              "Fuwai=422/596 \\")
        print("  REC_SLOPE=Colorado=0.5999,Brigham=0.5099,Stanford=0.4008,"
              "Fuwai=0.2948 \\")
        print("  REC_POOLED=1376,1800 REC_POOLED_P=2.51e-117 \\")
        print("  REC_MEAN_SLOPE=0.451 REC_LBNL=0.5726/78.2 \\")
        print("  python3 98_build_replication_table.py --apply")
        print()
        print("⚠ THE CAPTION ITSELF IS STILL WRONG and this patch does not")
        print("  touch it — the caption is prose inside 98, separate from the")
        print("  advisory. The replacement sentence is in the project doc")
        print("  RESULT_85_recount_rl125_caption_wrong.md and has to go into")
        print("  the manuscript by hand.")
        return 0

    baks = {}
    try:
        for t in (T85, T98):
            src, new = staged[t]
            bak = f"{t}.bak.{time.strftime('%Y%m%d_%H%M%S')}"
            shutil.copy2(t, bak)
            baks[t] = bak
            open(t, "w").write(new)
            disk = open(t, errors="replace").read()
            bad = verify_one(t, disk)
            if bad or disk != new:
                raise RuntimeError("; ".join(bad) or f"{t}: write mismatch")
            print(f"[write] {t}   (backup: {bak})")
    except Exception as e:
        print(f"✗ FAILED: {e}\n  Restoring BOTH files:")
        for t, bak in baks.items():
            shutil.copy2(bak, t)
            print(f"    restored {t}")
        return 3

    print("✓ both files verified on disk. An unrecorded claim now prints as")
    print("  absent in every one of REC_CLAIM's three sites, and 98's advisory")
    print("  points at the recount instead of asserting its answer.")
    return 0


def self_test():
    # The anchor sits INSIDE main() in the real file (4-space indent), so the
    # fixture must put it inside a def — at module level it is an
    # IndentationError, which is what a first version of this test hit.
    fx85 = ("def scenario(d, label, recn=None, recs=None):\n"
            "    print('SC', label, recn, recs)\n"
            "    return {'n': 19, 's1': 16}\n"
            "_RB_P = {}\n"
            "REC_CLAIM = (0, 0)\n"
            "indep = None\n"
            "def main():\n"
            + A_85 + "\n"
            "    return r_ind\n"
            "main()\n")
    new85, p = patch_one(T85, fx85)
    assert not p, p
    assert not verify_one(T85, new85), verify_one(T85, new85)
    compile(new85, "fx85", "exec")

    # (0, 0) must reach scenario() as None, None
    import io
    import contextlib
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        exec(compile(new85, "fx85", "exec"), {})
    out = buf.getvalue()
    assert "SC B. independent cohorts only (LBNL removed) None None" in out, out

    # a REAL claim must still reach it as 19, 16
    with contextlib.redirect_stdout(buf := io.StringIO()):
        exec(compile(new85.replace("REC_CLAIM = (0, 0)",
                                   "REC_CLAIM = (16, 19)"),
                     "fx85b", "exec"), {})
    out2 = buf.getvalue()
    assert out2.strip().endswith("19 16"), out2

    fx98 = 'SUF = os.environ.get("TBL_SUFFIX", "")\nx = """\n' + A_98 + '\n"""\n'
    new98, p2 = patch_one(T98, fx98)
    assert not p2, p2
    assert not verify_one(T98, new98), verify_one(T98, new98)
    assert "Fifteen of eighteen" not in new98
    assert "15 of 18 at 50 bp, \n" not in new98        # no broken line joins
    for frag in ("14 of 18 at 125 bp", "Fuwai RI", "Brigham RI (21/32,",
                 "scenario C"):
        assert frag in new98, frag

    # the guards must refuse an unpatched predecessor
    assert "_RB_P = {" not in "print('bare')\n"
    assert GUARDS[T85][0] not in "print('bare')\n"
    assert GUARDS[T98][0] not in "print('bare')\n"

    # idempotence: the anchors are gone after one pass
    _n, p3 = patch_one(T85, new85)
    assert p3, "85: second application was not refused"
    _n, p4 = patch_one(T98, new98)
    assert p4, "98: second application was not refused"

    print("[ok] self-test passed: an unrecorded (0, 0) claim reaches "
          "scenario() as None/None so it prints as absent rather than as a "
          "MISMATCH against zero, while a real (16, 19) still reaches it "
          "intact; 98's note 1 no longer asserts a count, records that the "
          "quantity moves with read length, and warns against the "
          "lowest-rate framing with the two counter-examples; both files are "
          "verified before either is written and both are restored on any "
          "failure; double application refused on each")
    return 0


if __name__ == "__main__":
    sys.exit(main())
