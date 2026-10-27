#!/usr/bin/env python3
"""
134_patch_85_separate_g2_g4.py -- fixes two flaws in my own patch 133.

FLAW 1: SUPPLYING REC_BORDERLINE WOULD MAKE G4 CIRCULAR
-------------------------------------------------------
133 made `REC_BORDERLINE` overridable and guarded G2 and G4 against an empty
one. But the two gates want different things, and 133 tied them together:

  * **G2 is a FIXITY check.** "Does the file still contain the combinations the
    record describes?" It needs k and n only. Running it at 125 bp is
    worthwhile — it catches the input table changing under us later.
  * **G4 is a PROVENANCE inference.** "Which convention produced the
    *published* p-values?" It compares a recorded p against the one- and
    two-sided values it computes. At 125 bp there is no published p, so
    supplying one means supplying a number 85 itself just printed in its
    one-sided column — and G4 would then "discover" ONE-SIDED from its own
    output. A guaranteed pass that tests nothing.

So `REC_BORDERLINE` entries may now omit `@p`:

    REC_BORDERLINE=Fuwai:A3SS=40/67,Stanford:A5SS=20/32      # G2 runs, G4 does not
    REC_BORDERLINE=Brigham:A5SS=19/29@0.07                   # both run

G4 iterates only the entries that carry a recorded p (`_RB_P`), and reports
NOT RUN when there are none — keeping 133's independent justification, which
stands on its own: a one-sided exact binomial at p = 0.5 reproduces every p in
`direction_concordance_dedup*.tsv` exactly, five cohorts, both read lengths.

This also fixes a crash 133 introduced: G4's table prints `{rp:8.2f}`, which
raises TypeError on a None.

FLAW 2: REC_CLAIM=none PRINTS "0 of 0" AND WARNS ABOUT IT
---------------------------------------------------------
With `REC_CLAIM=(0, 0)` the verdict block prints

    As published : "0 of 0"
    ⚠ Scenario B does NOT reproduce "0 of 0" — it gives 15 of 19.

Working as designed — scenario() warns rather than refuses — but it is noise,
and noise in a mismatch warning is how a real mismatch gets ignored later.
`(0, 0)` now prints "not recorded at this read length" and the comparison is
skipped.

THE MEASURED 125 bp BORDERLINE SET (k/n only, deliberately no p)
----------------------------------------------------------------
    Fuwai    A3SS  40/67  59.70%   p1 0.0710
    Stanford A5SS  20/32  62.50%   p1 0.1077
    Brigham  A5SS  16/25  64.00%   p1 0.1148
    Brigham  RI    21/32  65.62%   p1 0.0551

The p-values are recorded in the project doc, NOT passed to G4. They came from
85; they cannot also be 85's check on itself.

USAGE
-----
  python 134_patch_85_separate_g2_g4.py [--apply] [--self-test]
"""
import argparse
import difflib
import os
import shutil
import sys
import time

TARGET = "85_recount_combinations.py"

RUN125 = (
    "TBL_SUFFIX=_rl125 \\\n"
    "REC_SE=Colorado=311/385,LBNL=122/156,Brigham=298/364,"
    "Fuwai=422/596,Stanford=345/455 \\\n"
    "REC_SE_POOLED=1376,1800 \\\n"
    "REC_ALLTYPES=1883,2505 \\\n"
    "REC_BORDERLINE=Fuwai:A3SS=40/67,Stanford:A5SS=20/32,"
    "Brigham:A5SS=16/25,Brigham:RI=21/32 \\\n"
    "REC_CLAIM=none \\\n"
    "python3 85_recount_combinations.py")

# ------------------------------------------------------------------- edit 1
# the parser gains an optional @p, and _RB_P is derived
A_PARSE = '''        _nb = {}
        for _item in _env["REC_BORDERLINE"].split(","):
            _key, _val = _item.split("=", 1)
            _c, _ty = _key.split(":", 1)
            _kn, _p = _val.split("@", 1)
            _a, _b = _kn.split("/")
            _nb[(_c, _ty)] = (int(_a), int(_b), float(_p))
        REC_BORDERLINE = _nb
'''
B_PARSE = '''        _nb = {}
        for _item in _env["REC_BORDERLINE"].split(","):
            _key, _val = _item.split("=", 1)
            _c, _ty = _key.split(":", 1)
            # `@p` is OPTIONAL. With it, G2 and G4 both run. Without it, only
            # G2 runs — see the _RB_P note below. A recorded p at a NEW read
            # length can only come from this script's own one-sided column,
            # which would make G4 infer the convention from its own output.
            if "@" in _val:
                _kn, _p = _val.split("@", 1)
                _pv = float(_p)
            else:
                _kn, _pv = _val, None
            _a, _b = _kn.split("/")
            _nb[(_c, _ty)] = (int(_a), int(_b), _pv)
        REC_BORDERLINE = _nb
'''

# ------------------------------------------------------------------- edit 2
# _RB_P: the subset G4 may use. Appended to the record block's closing rule.
A_CLOSE = ("        REC_BORDERLINE = _nb\n"
           "# -----------------------------------------------------------"
           "---------------\n")
B_CLOSE = ('''        REC_BORDERLINE = _nb

# G4 may use ONLY the entries that carry a PUBLISHED p-value. G2 is a fixity
# check on k/n and uses all of them; G4 is a provenance inference and cannot
# be fed a number this script produced. An entry written without `@p` is a
# deliberate statement that no published p exists for it.
_RB_P = {_k: _v for _k, _v in REC_BORDERLINE.items() if _v[2] is not None}
'''
           "# -----------------------------------------------------------"
           "---------------\n")

# ------------------------------------------------------------------- edit 3
A_G4HDR = ('''    _rb_ps = " / ".join(f"{_rp:.2f}"
                       for (_, _, _rp) in REC_BORDERLINE.values())''')
B_G4HDR = ('''    _rb_ps = " / ".join(f"{_rp:.2f}"
                       for (_, _, _rp) in _RB_P.values())''')

# ------------------------------------------------------------------- edit 4
A_G4TBL = ('    for (c, ty), (rs, rn, rp) in REC_BORDERLINE.items():\n'
           "        print(f\"  {c+' '+ty:16s} {f'{rs}/{rn}':>9s} "
           '{p_one(rs,rn):10.4f} "\n')
B_G4TBL = ('    for (c, ty), (rs, rn, rp) in _RB_P.items():\n'
           "        print(f\"  {c+' '+ty:16s} {f'{rs}/{rn}':>9s} "
           '{p_one(rs,rn):10.4f} "\n')

# ------------------------------------------------------------------- edit 5
A_D12 = ('    d1 = sum(abs(p_one(rs, rn) - rp) for (rs, rn, rp) '
         'in REC_BORDERLINE.values())\n'
         '    d2 = sum(abs(p_two(rs, rn) - rp) for (rs, rn, rp) '
         'in REC_BORDERLINE.values())\n'
         '    if not REC_BORDERLINE:\n')
B_D12 = ('    d1 = sum(abs(p_one(rs, rn) - rp) for (rs, rn, rp) '
         'in _RB_P.values())\n'
         '    d2 = sum(abs(p_two(rs, rn) - rp) for (rs, rn, rp) '
         'in _RB_P.values())\n'
         '    if not _RB_P:\n')

# ------------------------------------------------------------------- edit 6
A_G4MSG = ('        print("\\n  G4 NOT RUN — no recorded p-values to infer '
           'a convention")\n')
B_G4MSG = ('        print("\\n  G4 NOT RUN — no PUBLISHED p-value is recorded '
           'for any")\n'
           '        print("  borderline combination, so there is no convention '
           "to infer.\")\n"
           '        print("  (A p-value taken from this script\'s own one-sided '
           'column")\n'
           '        print("   would make this gate confirm itself. Entries '
           'written as")\n'
           '        print("   k/n with no @p say so deliberately.)")\n'
           '        print("  Taking ONE-SIDED as established independently")\n')

# ------------------------------------------------------------------- edit 7
A_G4TAIL = ('        print("  from. Taking ONE-SIDED as already established: '
            'a one-sided")\n')
B_G4TAIL = ('        print("  elsewhere: a one-sided")\n')

# ------------------------------------------------------------------- edit 8
A_CLAIM = ('  As published : "{REC_CLAIM[0]} of {REC_CLAIM[1]}"')
B_CLAIM = ('  As published : {_claim_txt}')

EDITS = [
    ("the optional @p in the REC_BORDERLINE parser", A_PARSE, B_PARSE, 1),
    ("_RB_P, the subset G4 may use", A_CLOSE, B_CLOSE, 1),
    ("the G4 header, over _RB_P", A_G4HDR, B_G4HDR, 1),
    ("the G4 table, over _RB_P (and no TypeError on a None p)",
     A_G4TBL, B_G4TBL, 1),
    ("the sidedness sums, over _RB_P", A_D12, B_D12, 1),
    ("the G4 NOT RUN message, saying why", A_G4MSG, B_G4MSG, 1),
    ("the G4 NOT RUN tail", A_G4TAIL, B_G4TAIL, 1),
    ("the verdict's published-claim line", A_CLAIM, B_CLAIM, 1),
]

# edit 9 is inserted separately because it needs an anchor inside main()
A_VERDICT = '    print(f"""\n  VERDICT on the claim sentence'
B_VERDICT = ('    _claim_txt = ("not recorded at this read length"\n'
             '                  if tuple(REC_CLAIM) == (0, 0)\n'
             '                  else f\'"{REC_CLAIM[0]} of {REC_CLAIM[1]}"\')\n'
             '    print(f"""\n  VERDICT on the claim sentence')
EDITS.append(("_claim_txt, before the verdict block", A_VERDICT, B_VERDICT, 1))

# edit 10 -- the mismatch warning must not fire against (0, 0)
A_WARN = '    if r_ind["n"] != REC_CLAIM[1] or r_ind["s1"] != REC_CLAIM[0]:'
B_WARN = ('    if tuple(REC_CLAIM) != (0, 0) and (\n'
          '            r_ind["n"] != REC_CLAIM[1] or '
          'r_ind["s1"] != REC_CLAIM[0]):')
EDITS.append(("the claim mismatch warning, skipped when unrecorded",
              A_WARN, B_WARN, 1))


def apply_edits(text):
    if "_RB_P = {" in text:
        return text, ["_RB_P is already defined — 85 already carries 134 "
                      "(restore the .bak to re-apply)"]
    if 'SUF = os.environ.get("TBL_SUFFIX", "")' not in text:
        return text, ["133 has not been applied — apply "
                      "133_patch_85_recount.py first"]
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
    for w in ("_RB_P = {_k: _v for _k, _v in REC_BORDERLINE.items()",
              'if "@" in _val:',
              "_kn, _pv = _val, None",
              "for (c, ty), (rs, rn, rp) in _RB_P.items():",
              "in _RB_P.values())",
              "    if not _RB_P:",
              "no PUBLISHED p-value is recorded",
              "would make this gate confirm itself",
              "_claim_txt = (",
              "not recorded at this read length",
              "  As published : {_claim_txt}",
              'if tuple(REC_CLAIM) != (0, 0) and ('):
        if w not in text:
            bad.append(f"missing: {w[:58]}")
    for gone in ("in REC_BORDERLINE.values())",
                 "for (c, ty), (rs, rn, rp) in REC_BORDERLINE.items():\n"
                 "        print(f\"  {c+' '+ty:16s}",
                 '  As published : "{REC_CLAIM[0]} of {REC_CLAIM[1]}"',
                 "    if not REC_BORDERLINE:\n        # d1 == d2"):
        if gone in text:
            bad.append(f"stale site survives: {gone[:52]}")
    # G2's loop must STILL use the full dict -- it is a fixity check on k/n
    if "for (c, ty), (rs, rn, rp) in REC_BORDERLINE.items():" not in text:
        bad.append("G2's loop no longer iterates the full REC_BORDERLINE — "
                   "it is a fixity check and must use every entry")
    # G2's own guard must still key off the full dict, not _RB_P
    i_g2 = text.find("G2 NOT RUN")
    if i_g2 < 0 or "if not REC_BORDERLINE:" not in text[:i_g2]:
        bad.append("G2's empty-guard no longer keys off the full dict")
    # _RB_P must be defined above every use
    i_p = text.find("_RB_P = {_k: _v")
    for m in ("in _RB_P.values())", "if not _RB_P:",
              "for (c, ty), (rs, rn, rp) in _RB_P.items():"):
        j = text.find(m)
        if j < 0 or i_p > j:
            bad.append(f"_RB_P not defined above: {m[:40]}")
    # _claim_txt must be assigned before it is formatted
    i_a = text.find("_claim_txt = (")
    i_u = text.find("As published : {_claim_txt}")
    if i_a < 0 or i_u < 0 or i_a > i_u:
        bad.append("_claim_txt is used before it is assigned")
    return bad


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
    print("134 — separating G2's fixity check from G4's provenance inference")
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
        print("\n  Nothing written.")
        return 3
    bad = verify(new)
    if bad:
        print("REFUSING — patched text fails verification:")
        for b in bad:
            print(f"  ✗ {b}")
        return 3

    print("".join(difflib.unified_diff(
        src.splitlines(True), new.splitlines(True),
        fromfile=a.target, tofile=a.target + " (patched)", n=2)))
    print("=" * 78)

    if not a.apply:
        print("DRY RUN — nothing written. Re-run with --apply, then:")
        print()
        print("  # 1. CALIBRATION — 50 bp has published p-values, so G4 must")
        print("  #    STILL RUN and still conclude ONE-SIDED")
        print("  python3 85_recount_combinations.py 2>&1 | "
              "grep -A8 'G4  sidedness'")
        print("  python3 85_recount_combinations.py 2>&1 | grep 'As published'")
        print()
        print("  # 2. 125 bp — G2 RUNS (k/n fixity), G4 does NOT (no published p)")
        print("  " + RUN125.replace("\n", "\n  "))
        print()
        print("⚠ Expect, at 125 bp:")
        print("    G2 PASS                     <- the four cells reproduce")
        print("    G4 NOT RUN                  <- and it says why")
        print("    As published : not recorded at this read length")
        print("    no '0 of 0' and no spurious mismatch warning")
        print()
        print("⚠ The four p-values (0.0710 / 0.1077 / 0.1148 / 0.0551) are")
        print("  recorded in the project doc and deliberately NOT passed to")
        print("  G4. They came out of 85; they cannot also be 85's check on")
        print("  itself.")
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
    print("✓ verified on disk: G2 checks every recorded combination, G4 only")
    print("  those with a published p, and an unrecorded claim prints as")
    print("  unrecorded instead of as 0 of 0.")
    return 0


def self_test():
    import subprocess
    import tempfile

    # The fixture reproduces 85 AS 133 LEFT IT, so 134's anchors are tested
    # against 133's output rather than against text I invent here.
    fx = (
        "import os\n"
        "def p_one(k, n):\n"
        "    return {(19, 29): 0.0712, (42, 71): 0.0804, "
        "(39, 68): 0.1399}.get((k, n), 0.5)\n"
        "def p_two(k, n):\n"
        "    return 2 * p_one(k, n)\n"
        'REC_BORDERLINE = {\n'
        '    ("Brigham", "A5SS"): (19, 29, 0.07),\n'
        '    ("Fuwai", "A3SS"):   (42, 71, 0.08),\n'
        '    ("Fuwai", "A5SS"):   (39, 68, 0.14),\n'
        "}\n"
        'REC_CLAIM = (16, 19)\n'
        'SUF = os.environ.get("TBL_SUFFIX", "")\n'
        '_env = {"REC_BORDERLINE": os.environ.get("REC_BORDERLINE", "").strip(),\n'
        '        "REC_CLAIM": os.environ.get("REC_CLAIM", "").strip()}\n'
        'if _env["REC_CLAIM"]:\n'
        '    REC_CLAIM = ((0, 0) if _env["REC_CLAIM"].lower() == "none"\n'
        '                 else tuple(int(_x) for _x in _env["REC_CLAIM"].split(",")))\n'
        'if _env["REC_BORDERLINE"]:\n'
        '    if _env["REC_BORDERLINE"].lower() == "none":\n'
        "        REC_BORDERLINE = {}\n"
        "    else:\n"
        + A_PARSE
        + A_CLOSE
        + "def main():\n"
          "    r_ind = {'n': 19, 's1': 16}\n"
          "    for (c, ty), (rs, rn, rp) in REC_BORDERLINE.items():\n"
          '        print("G2ROW", c, ty, rs, rn)\n'
          "    if not REC_BORDERLINE:\n"
          '        print("  G2 NOT RUN — REC_BORDERLINE is empty, so nothing was")\n'
          "    else:\n"
          '        print("  G2 PASS")\n'
        + "    " + A_G4HDR.strip() + "\n"
          '    print("G4HDR", _rb_ps or "(empty)")\n'
        + A_G4TBL
        + '              f"{p_two(rs,rn):10.4f} {rp:8.2f}")\n'
        + A_D12
        + '        sided = "ONE-SIDED"\n'
        + A_G4MSG
        + A_G4TAIL
        + '        print("  exact binomial ...")\n'
          "    else:\n"
          '        sided = "ONE-SIDED" if d1 < d2 else "TWO-SIDED"\n'
          '        print(f"\\n  total absolute deviation: one-sided {d1:.4f}, "\n'
          '              f"two-sided {d2:.4f}")\n'
          '    print(f"  G4 — the record\'s p-values are {sided}.")\n'
          "    if r_ind['n'] != REC_CLAIM[1] or r_ind['s1'] != REC_CLAIM[0]:\n"
          '        print("WARNFIRED")\n'
        + A_VERDICT + '\n'
          '  As published : "{REC_CLAIM[0]} of {REC_CLAIM[1]}"\n'
          '  end""")\n'
          "if __name__ == '__main__':\n"
          "    main()\n")

    # the fixture's warn line uses r_ind[...] with single quotes; normalise so
    # edit 10's anchor (double quotes, as in the real file) is present
    fx = fx.replace(
        "    if r_ind['n'] != REC_CLAIM[1] or r_ind['s1'] != REC_CLAIM[0]:",
        '    if r_ind["n"] != REC_CLAIM[1] or r_ind["s1"] != REC_CLAIM[0]:')

    new, probs = apply_edits(fx)
    assert not probs, probs
    assert not verify(new), verify(new)
    compile(new, "fx", "exec")

    d = tempfile.mkdtemp()
    p = os.path.join(d, "fx.py")
    open(p, "w").write(new)

    def run(env=None):
        e = dict(os.environ)
        for k in ("TBL_SUFFIX", "REC_BORDERLINE", "REC_CLAIM"):
            e.pop(k, None)
        e.update(env or {})
        return subprocess.run([sys.executable, "-I", p], capture_output=True,
                              text=True, cwd=d, env=e)

    # 1. unsuffixed: three published p-values -> G2 and G4 BOTH run
    r = run()
    assert r.returncode == 0, r.stdout + r.stderr
    assert "G2 PASS" in r.stdout, r.stdout
    assert "G4HDR 0.07 / 0.08 / 0.14" in r.stdout, r.stdout
    assert "total absolute deviation" in r.stdout, r.stdout
    assert "G4 NOT RUN" not in r.stdout, r.stdout
    assert "p-values are ONE-SIDED" in r.stdout, r.stdout
    assert 'As published : "16 of 19"' in r.stdout, r.stdout
    assert "WARNFIRED" not in r.stdout, r.stdout

    # 2. k/n WITHOUT @p -> G2 runs, G4 does NOT, and no TypeError on the None
    r2 = run({"REC_BORDERLINE": "Fuwai:A3SS=40/67,Stanford:A5SS=20/32,"
                                "Brigham:A5SS=16/25,Brigham:RI=21/32",
              "REC_CLAIM": "none"})
    assert r2.returncode == 0, r2.stdout + r2.stderr
    assert "G2ROW Fuwai A3SS 40 67" in r2.stdout, r2.stdout
    assert "G2ROW Brigham RI 21 32" in r2.stdout, r2.stdout
    assert "G2 PASS" in r2.stdout, r2.stdout
    assert "G2 NOT RUN" not in r2.stdout, r2.stdout
    assert "G4 NOT RUN" in r2.stdout, r2.stdout
    assert "no PUBLISHED p-value is recorded" in r2.stdout, r2.stdout
    assert "confirm itself" in r2.stdout, r2.stdout
    assert "total absolute deviation" not in r2.stdout, r2.stdout
    assert "TWO-SIDED" not in r2.stdout, r2.stdout
    assert "G4HDR (empty)" in r2.stdout, r2.stdout
    # the unrecorded claim must print as unrecorded, with no warning
    assert "As published : not recorded at this read length" in r2.stdout, \
        r2.stdout
    assert "0 of 0" not in r2.stdout, r2.stdout
    assert "WARNFIRED" not in r2.stdout, r2.stdout

    # 3. a MIXED set: only the @p entry reaches G4, both reach G2
    r3 = run({"REC_BORDERLINE": "Fuwai:A3SS=42/71@0.08,Brigham:RI=21/32"})
    assert r3.returncode == 0, r3.stdout + r3.stderr
    assert "G2ROW Fuwai A3SS 42 71" in r3.stdout, r3.stdout
    assert "G2ROW Brigham RI 21 32" in r3.stdout, r3.stdout
    assert "G4HDR 0.08" in r3.stdout, r3.stdout
    assert "G4 NOT RUN" not in r3.stdout, r3.stdout
    assert "total absolute deviation" in r3.stdout, r3.stdout
    assert "p-values are ONE-SIDED" in r3.stdout, r3.stdout

    # 4. REC_BORDERLINE=none still gives G2 NOT RUN (133's behaviour kept)
    r4 = run({"REC_BORDERLINE": "none", "REC_CLAIM": "none"})
    assert r4.returncode == 0, r4.stdout + r4.stderr
    assert "G2 NOT RUN" in r4.stdout and "G2 PASS" not in r4.stdout, r4.stdout
    assert "G4 NOT RUN" in r4.stdout, r4.stdout

    # 5. a real REC_CLAIM still compares, and a genuine mismatch still warns
    r5 = run({"REC_CLAIM": "15,18"})
    assert r5.returncode == 0, r5.stdout + r5.stderr
    assert 'As published : "15 of 18"' in r5.stdout, r5.stdout
    assert "WARNFIRED" in r5.stdout, r5.stdout

    # 6. double application refused, text unchanged
    n6, p6 = apply_edits(new)
    assert p6, "second application was not refused"
    assert n6 == new, "a refused re-application still altered the text"

    # 7. refuses outright if 133 was never applied
    _n7, p7 = apply_edits("print('no SUF here')\n")
    assert p7 and "133 has not been applied" in p7[0], p7

    print("[ok] self-test passed: with published p-values both gates run and "
          "the 50 bp output is unchanged; with k/n and no @p, G2 runs on every "
          "entry while G4 reports NOT RUN and says why; a mixed set sends only "
          "the @p entry to G4; REC_BORDERLINE=none still leaves G2 NOT RUN; an "
          "unrecorded REC_CLAIM prints as unrecorded with no '0 of 0' and no "
          "spurious warning, while a real one still compares and still warns; "
          "no TypeError on a None p; double application refused; refuses if "
          "133 was never applied")
    return 0


if __name__ == "__main__":
    sys.exit(main())
