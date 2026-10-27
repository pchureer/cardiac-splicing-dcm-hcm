#!/usr/bin/env python3
"""
133_patch_85_recount.py -- 85 at 125 bp. Two gates would report a conclusion
from no data, and one would crash.

85 is the script 98's caption cites for "fifteen of eighteen cohort x
event-type combinations are individually significant". Everything it checks is
recorded 50 bp values, and three of its gates fail unsafely rather than loudly.

THREE DEFECTS, NOT JUST STALE CONSTANTS
---------------------------------------
1. **G4 announces a false conclusion from an empty record.** It infers whether
   the recorded p-values are one- or two-sided by

       d1 = sum(abs(p_one(rs, rn) - rp) for ... in REC_BORDERLINE.values())
       d2 = sum(abs(p_two(rs, rn) - rp) for ... in REC_BORDERLINE.values())
       sided = "ONE-SIDED" if d1 < d2 else "TWO-SIDED"

   With `REC_BORDERLINE` empty both sums are 0.0, `d1 < d2` is False, and the
   script PRINTS "the record's p-values are TWO-SIDED" having tested nothing.
   At a new read length the three borderline combinations are unknown — that
   is what G5 exists to discover — so an empty dict is the expected state.

   Sidedness is settled independently anyway: a ONE-SIDED exact binomial at
   p = 0.5 reproduces every p in `direction_concordance_dedup*.tsv` exactly,
   all five cohorts, both read lengths. G4 now says so instead of guessing.

2. **G2 passes vacuously.** Empty `REC_BORDERLINE` means the loop body never
   runs, `ok2` stays True, and it prints "G2 PASS". It now prints
   "G2 NOT RUN" and says what to do to turn it back on.

3. **G3 crashes when `REC_SE` lacks a cohort.** `rec = REC_SE.get(c)` may be
   None; `hit` handles that, the f-string `record {rec[0]}/{rec[1]}` does not.

WHAT MOVES
----------
  62-65  REC_BORDERLINE   3 combinations     -> REC_BORDERLINE (may be `none`)
  68-69  REC_SE           5 cohorts          -> REC_SE
  70     REC_POOLED_ALLTYPES (2227, 2991)    -> REC_ALLTYPES
  71     REC_CLAIM        (16, 19)           -> REC_CLAIM (may be `none`)
  74/75  COMB / SE1       input filenames    -> suffixed
  ~257   the G4 header    "0.07 / 0.08 / 0.14" -> derived from REC_BORDERLINE
  ~328   inline rec=(1636, 2163)             -> REC_SE_POOLED
  ~331   prose "(1,636/2,163 = 75.6%)"       -> derived from REC_SE_POOLED

`REC_SE_POOLED` is new: 1,636/2,163 appeared as an **inline tuple** at line 328
and again as **prose** at line 331, neither of them the named `REC_SE` — Rule 8.
One constant now drives both.

`INDEPENDENT` and `EXCLUDED_MXE` are not touched: which cohorts are
independent, and which lose MXE on reverse-replication grounds, do not depend
on read length.

WHY `none` IS REQUIRED RATHER THAN DEFAULTED
--------------------------------------------
When `TBL_SUFFIX` is set, all five records must be given. Two of them —
`REC_BORDERLINE` and `REC_CLAIM` — may be the literal string `none`, because
at a new read length they are genuinely unknown. Requiring the operator to
type `none` forces an explicit decision instead of letting a stale 50 bp
expectation refuse a correct recount, or an empty one pass silently.

MEASURED AT 125 bp
------------------
    Colorado 311/385   LBNL 122/156   Brigham 298/364
    Fuwai    422/596   Stanford 345/455
    SE pooled, 4 independent      1,376/1,800 = 76.4%
    all five types, 4 independent 1,883/2,505 = 75.2%   (was 2,227/2,991)

USAGE
-----
  python 133_patch_85_recount.py [--apply] [--self-test]
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
    "REC_BORDERLINE=none \\\n"
    "REC_CLAIM=none \\\n"
    "python3 85_recount_combinations.py")

# ------------------------------------------------------------------- edit 1
A_REC = ('REC_POOLED_ALLTYPES = (2227, 2991)      '
         "# record's all-five-types pooled figure\n"
         'REC_CLAIM = (16, 19)                    # "16 of 19"\n')
B_REC = A_REC + '''
# NEW: 1,636/2,163 occurred as an inline tuple (line ~328) and as prose
# (line ~331), neither of them the named REC_SE — Rule 8. One constant now.
REC_SE_POOLED = (1636, 2163)

# --------------------------------------------------------------------------
# READ LENGTH and the 125 bp record, placed AFTER every literal it rewrites.
#
# When TBL_SUFFIX is set all five records are REQUIRED. REC_BORDERLINE and
# REC_CLAIM may be the literal `none`: at a new read length they are unknown
# (G5 is what discovers them), and an EMPTY REC_BORDERLINE is not a neutral
# default — see the G2 and G4 guards further down, which this patch also adds.
SUF = os.environ.get("TBL_SUFFIX", "")
_REQ = ("REC_SE", "REC_SE_POOLED", "REC_ALLTYPES", "REC_BORDERLINE",
        "REC_CLAIM")
_env = {_k: os.environ.get(_k, "").strip() for _k in _REQ}
if SUF:
    _miss = [_k for _k in _REQ if not _env[_k]]
    if _miss:
        raise SystemExit(
            "REFUSING — TBL_SUFFIX=%s but these are unset: %s\\n\\n"
            "  REC_BORDERLINE and REC_CLAIM accept the literal `none` when the\\n"
            "  value is not yet known at this read length. Say so explicitly —\\n"
            "  an empty record is NOT a neutral default here.\\n\\n"
            "    %s\\n" % (SUF, ", ".join(_miss), @@RUN125@@))

if _env["REC_SE"]:
    _d = dict(_p.split("=", 1) for _p in _env["REC_SE"].split(","))
    assert set(_d) == set(REC_SE), (
        "REC_SE must name exactly %s, got %s" % (sorted(REC_SE), sorted(_d)))
    REC_SE = {_c: tuple(int(_x) for _x in _v.split("/"))
              for _c, _v in _d.items()}

if _env["REC_SE_POOLED"]:
    REC_SE_POOLED = tuple(int(_x) for _x in _env["REC_SE_POOLED"].split(","))
if _env["REC_ALLTYPES"]:
    REC_POOLED_ALLTYPES = tuple(int(_x)
                                for _x in _env["REC_ALLTYPES"].split(","))

# REC_CLAIM: `none` -> (0, 0). The scenario() call already WARNS rather than
# refuses on a mismatch, and (0, 0) reads as "no published claim to compare".
if _env["REC_CLAIM"]:
    REC_CLAIM = ((0, 0) if _env["REC_CLAIM"].lower() == "none"
                 else tuple(int(_x) for _x in _env["REC_CLAIM"].split(",")))

# REC_BORDERLINE: `none` -> {}. Format otherwise
#   cohort:TYPE=k/n@p,cohort:TYPE=k/n@p
if _env["REC_BORDERLINE"]:
    if _env["REC_BORDERLINE"].lower() == "none":
        REC_BORDERLINE = {}
    else:
        _nb = {}
        for _item in _env["REC_BORDERLINE"].split(","):
            _key, _val = _item.split("=", 1)
            _c, _ty = _key.split(":", 1)
            _kn, _p = _val.split("@", 1)
            _a, _b = _kn.split("/")
            _nb[(_c, _ty)] = (int(_a), int(_b), float(_p))
        REC_BORDERLINE = _nb
# --------------------------------------------------------------------------
'''

# ------------------------------------------------------------------- edit 2
A_FILES = ('COMB = "cohort_type_concordance_dedup.tsv"\n'
           'SE1 = "direction_concordance_dedup.tsv"\n')
B_FILES = ('COMB = f"cohort_type_concordance_dedup{SUF}.tsv"\n'
           'SE1 = f"direction_concordance_dedup{SUF}.tsv"\n')

# ------------------------------------------------------------------- edit 3
A_G2PASS = '    print("  G2 PASS")\n'
B_G2PASS = '''    if not REC_BORDERLINE:
        # An empty dict makes the loop above a no-op and leaves ok2 True, so
        # this gate used to print "G2 PASS" having checked nothing at all.
        print("  G2 NOT RUN — REC_BORDERLINE is empty, so nothing was")
        print("  reproduced. This is NOT a pass. At a new read length the")
        print("  borderline combinations are unknown until G5 below names")
        print("  them; read them off G5's 'Not individually significant'")
        print("  line and re-run with")
        print("    REC_BORDERLINE=Cohort:TYPE=k/n@p,Cohort:TYPE=k/n@p")
    else:
        print("  G2 PASS")
'''

# ------------------------------------------------------------------- edit 4
A_G4HDR = ('print("G4  sidedness: which convention produced the recorded '
           '0.07 / 0.08 / 0.14?")')
B_G4HDR = ('_rb_ps = " / ".join(f"{_rp:.2f}"\n'
           '                       for (_, _, _rp) in REC_BORDERLINE.values())\n'
           '    print("G4  sidedness: which convention produced the recorded "\n'
           '          + (_rb_ps + "?" if _rb_ps else "p-values? (none recorded)"))')

# ------------------------------------------------------------------- edit 5
A_SIDED = ('    sided = "ONE-SIDED" if d1 < d2 else "TWO-SIDED"\n'
           '    print(f"\\n  total absolute deviation: one-sided {d1:.4f}, '
           'two-sided {d2:.4f}")\n')
B_SIDED = '''    if not REC_BORDERLINE:
        # d1 == d2 == 0.0 for an empty record, and `d1 < d2` is then False, so
        # this gate used to ANNOUNCE "TWO-SIDED" from no data whatsoever.
        # Sidedness is settled independently: a one-sided exact binomial at
        # p = 0.5 reproduces every p in direction_concordance_dedup*.tsv
        # exactly — all five cohorts, both read lengths.
        sided = "ONE-SIDED"
        print("\\n  G4 NOT RUN — no recorded p-values to infer a convention")
        print("  from. Taking ONE-SIDED as already established: a one-sided")
        print("  exact binomial at p=0.5 reproduces every p in")
        print(f"  {SE1} exactly, at both read lengths.")
    else:
        sided = "ONE-SIDED" if d1 < d2 else "TWO-SIDED"
        print(f"\\n  total absolute deviation: one-sided {d1:.4f}, "
              f"two-sided {d2:.4f}")
'''

EDITS = [
    ("SUF, REC_SE_POOLED and the 125 bp record", A_REC, B_REC, 1),
    ("the two input filenames, suffixed", A_FILES, B_FILES, 1),
    ("G2 must not pass vacuously", A_G2PASS, B_G2PASS, 1),
    ("the G4 header, derived", A_G4HDR, B_G4HDR, 1),
    ("G4 must not conclude from an empty record", A_SIDED, B_SIDED, 1),
    ("G3's crash when REC_SE lacks a cohort",
     "            rec = REC_SE.get(c)\n",
     "            rec = REC_SE.get(c)\n"
     '            _rec_txt = f"{rec[0]}/{rec[1]}" if rec else "(none)"\n', 1),
    ("G3's record column, None-tolerant",
     'f"record {rec[0]}/{rec[1]}   ',
     'f"record {_rec_txt}   ', 1),
    ("the inline SE pooled tuple",
     "rec=(1636, 2163))",
     "rec=REC_SE_POOLED)", 1),
    ("the SE headline prose, derived",
     'print("  The SE headline (1,636/2,163 = 75.6%) is unaffected '
     '— it never included MXE.")',
     'print(f"  The SE headline ({REC_SE_POOLED[0]:,}/{REC_SE_POOLED[1]:,} = "\n'
     '          f"{100 * REC_SE_POOLED[0] / REC_SE_POOLED[1]:.1f}%) is '
     'unaffected — it never included MXE.")', 1),
]

B_REC = B_REC.replace("@@RUN125@@", repr(RUN125.replace("\n", "\n    ")))
assert "@@RUN125@@" not in B_REC and "TBL_SUFFIX=_rl125" in B_REC
EDITS[0] = ("SUF, REC_SE_POOLED and the 125 bp record", A_REC, B_REC, 1)


def apply_edits(text):
    if 'SUF = os.environ.get("TBL_SUFFIX", "")' in text:
        return text, ["SUF is already defined — 85 is already patched "
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
    for w in ('SUF = os.environ.get("TBL_SUFFIX", "")',
              "REC_SE_POOLED = (1636, 2163)",
              '_REQ = ("REC_SE", "REC_SE_POOLED", "REC_ALLTYPES", '
              '"REC_BORDERLINE",',
              'COMB = f"cohort_type_concordance_dedup{SUF}.tsv"',
              'SE1 = f"direction_concordance_dedup{SUF}.tsv"',
              "G2 NOT RUN",
              "G4 NOT RUN",
              "rec=REC_SE_POOLED)",
              "_rec_txt = ",
              "f\"record {_rec_txt}",
              "REC_SE_POOLED[0]:,"):
        if w not in text:
            bad.append(f"missing: {w[:58]}")
    for gone in ('COMB = "cohort_type_concordance_dedup.tsv"',
                 'SE1 = "direction_concordance_dedup.tsv"',
                 "rec=(1636, 2163))",
                 "recorded 0.07 / 0.08 / 0.14?",
                 "The SE headline (1,636/2,163 = 75.6%)",
                 'f"record {rec[0]}/{rec[1]}'):
        if gone in text:
            bad.append(f"stale site survives: {gone[:50]}")
    # the empty-record guard must come BEFORE `sided` is used
    i_g = text.find("G4 NOT RUN")
    i_u = text.find("the record's p-values are {sided}")
    if i_g < 0 or i_u < 0 or i_g > i_u:
        bad.append("the G4 guard does not precede the use of `sided`")
    # neither guard may have been dropped into the wrong gate
    if text.find("G2 NOT RUN") > text.find("G4 NOT RUN"):
        bad.append("the G2 and G4 guards are in the wrong order")
    # SUF module-level, above def main()
    i_s = text.find('SUF = os.environ.get("TBL_SUFFIX", "")')
    i_m = text.find("def main(")
    if i_s < 0 or i_m < 0 or i_s > i_m:
        bad.append("SUF is not module-level / not above def main()")
    if "\n    SUF = os.environ.get" in text:
        bad.append("SUF is indented — it would be a local, not a global")
    # the record block must follow every literal it rewrites
    i_o = text.find("_env = {_k: os.environ.get(_k")
    for lit, nm in (("REC_SE = {", "REC_SE"),
                    ("REC_POOLED_ALLTYPES = (2227, 2991)", "REC_POOLED_ALLTYPES"),
                    ("REC_CLAIM = (16, 19)", "REC_CLAIM"),
                    ("REC_BORDERLINE = {", "REC_BORDERLINE")):
        i_l = text.find(lit)
        if i_l < 0 or i_o < 0 or i_l > i_o:
            bad.append(f"the record block is not after the {nm} literal")
    # COMB/SE1 must be defined after SUF
    for m in ('COMB = f"cohort_type_concordance_dedup{SUF}.tsv"',
              'SE1 = f"direction_concordance_dedup{SUF}.tsv"'):
        j = text.find(m)
        if j < 0 or i_s > j:
            bad.append(f"SUF not defined above: {m[:40]}")
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
    print("133 — 85_recount_combinations.py: two unsafe gates and a crash")
    print("=" * 78)
    if not os.path.exists(a.target):
        print(f"[fatal] {a.target} not found.", file=sys.stderr)
        return 2
    src = open(a.target, errors="replace").read()
    if "import os" not in src:
        print("[fatal] no `import os`", file=sys.stderr)
        return 2
    new, problems = apply_edits(src)
    if problems:
        print("REFUSING — anchors do not match:")
        for p in problems:
            print(f"  ✗ {p}")
        print("\n  Nothing written. Send me the lines the failing anchors name.")
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
        print("  # 1. CALIBRATION — unsuffixed must reproduce, G2 and G4 both RUN")
        print("  python3 85_recount_combinations.py 2>&1 | tail -60")
        print()
        print("  # 2. the refusal must name both `none`-able records")
        print("  TBL_SUFFIX=_rl125 python3 85_recount_combinations.py "
              "2>&1 | head -12")
        print()
        print("  # 3. 125 bp — G2 and G4 must say NOT RUN, not PASS")
        print("  " + RUN125.replace("\n", "\n  "))
        print()
        print("⚠ READ G5's 'Not individually significant' LINE. That is the")
        print("  125 bp borderline set, and 98's caption claims 'fifteen of")
        print("  eighteen ... the three that are not are the three LOWEST")
        print("  concordance rates (57-66%), all alternative-splice-site")
        print("  classes'. Every clause of that sentence has to be rechecked")
        print("  against the recount, including the count, the rates and the")
        print("  claim that all three are A3SS/A5SS.")
        print()
        print("⚠ Then re-run with REC_BORDERLINE set from that line, so G2 and")
        print("  G4 are real gates again rather than NOT RUN.")
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
    print("✓ verified on disk: G2 can no longer pass vacuously, G4 can no")
    print("  longer announce a sidedness from an empty record, G3 cannot crash")
    print("  on a missing cohort, and 1,636/2,163 has one source instead of")
    print("  three.")
    return 0


def self_test():
    import subprocess
    import tempfile

    fx = (
        "import os\n"
        "def p_one(k, n):\n"
        "    return 0.07 if (k, n) == (19, 29) else 0.50\n"
        "def p_two(k, n):\n"
        "    return 0.99\n"
        "def refuse(m):\n"
        "    raise SystemExit('REFUSING ' + m)\n"
        'REC_BORDERLINE = {\n'
        '    ("Brigham", "A5SS"): (19, 29, 0.07),\n'
        '    ("Fuwai", "A3SS"):   (42, 71, 0.08),\n'
        '    ("Fuwai", "A5SS"):   (39, 68, 0.14),\n'
        "}\n"
        'REC_SE = {"Colorado": (375, 471), "LBNL": (151, 199), '
        '"Brigham": (359, 442),\n'
        '          "Fuwai": (494, 714), "Stanford": (408, 536)}\n'
        + A_REC
        + A_FILES
        + "def main():\n"
          '    print("COMB", COMB, "SE1", SE1)\n'
          '    print("ALLTYPES", REC_POOLED_ALLTYPES, "CLAIM", REC_CLAIM)\n'
          '    print("SEPOOL", REC_SE_POOLED, "SE_COL", REC_SE["Colorado"])\n'
          "    ok2 = True\n"
          "    for (c, ty), (rs, rn, rp) in REC_BORDERLINE.items():\n"
          '        print("BORD", c, ty, rs, rn, rp)\n'
          "    if not ok2:\n"
          '        refuse("G2 — x")\n'
        + A_G2PASS
        + "    " + A_G4HDR + "\n"
          "    d1 = sum(abs(p_one(rs, rn) - rp) for (rs, rn, rp) "
          "in REC_BORDERLINE.values())\n"
          "    d2 = sum(abs(p_two(rs, rn) - rp) for (rs, rn, rp) "
          "in REC_BORDERLINE.values())\n"
        + A_SIDED
        + '    print(f"  G4 — the record\'s p-values are {sided}.")\n'
          "    for c in list(REC_SE) + ['Ghost']:\n"
          "            rec = REC_SE.get(c)\n"
          "            a = (1, 2)\n"
          "            b = (1, 2)\n"
          "            hit = a == b and (rec is None or a == rec)\n"
          '            print(f"  {c:10s} x {a[0]}/{a[1]}   '
          'f {b[0]}/{b[1]}   "\n'
          '                  f"record {rec[0]}/{rec[1]}   '
          "{'match' if hit else 'MISMATCH'}\")\n"
          "    def pooled(sub, label, rec=None):\n"
          '        print("POOLED", label, rec)\n'
          "        return (1, 2)\n"
          '    pooled(None, "SE only (the headline)",\n'
          "           rec=(1636, 2163))\n"
          '    print("  The SE headline (1,636/2,163 = 75.6%) is unaffected '
          '— it never included MXE.")\n'
          "if __name__ == '__main__':\n"
          "    main()\n")

    new, probs = apply_edits(fx)
    assert not probs, probs
    assert not verify(new), verify(new)
    compile(new, "fx", "exec")

    d = tempfile.mkdtemp()
    p = os.path.join(d, "fx.py")
    open(p, "w").write(new)

    def run(env=None):
        e = dict(os.environ)
        for k in ("TBL_SUFFIX",) + tuple(
                ("REC_SE", "REC_SE_POOLED", "REC_ALLTYPES",
                 "REC_BORDERLINE", "REC_CLAIM")):
            e.pop(k, None)
        e.update(env or {})
        return subprocess.run([sys.executable, "-I", p], capture_output=True,
                              text=True, cwd=d, env=e)

    # 1. unsuffixed -> original values, BOTH gates actually run
    r = run()
    assert r.returncode == 0, r.stdout + r.stderr
    assert "COMB cohort_type_concordance_dedup.tsv" in r.stdout, r.stdout
    assert "SE1 direction_concordance_dedup.tsv" in r.stdout, r.stdout
    assert "ALLTYPES (2227, 2991) CLAIM (16, 19)" in r.stdout, r.stdout
    assert "SEPOOL (1636, 2163) SE_COL (375, 471)" in r.stdout, r.stdout
    assert "G2 PASS" in r.stdout and "G2 NOT RUN" not in r.stdout, r.stdout
    assert "recorded 0.07 / 0.08 / 0.14?" in r.stdout, r.stdout
    assert "total absolute deviation" in r.stdout, r.stdout
    assert "G4 NOT RUN" not in r.stdout, r.stdout
    assert "p-values are ONE-SIDED" in r.stdout, r.stdout
    assert "POOLED SE only (the headline) (1636, 2163)" in r.stdout, r.stdout
    assert "The SE headline (1,636/2,163 = 75.6%)" in r.stdout, r.stdout
    # G3's None-tolerance: a cohort absent from REC_SE must print, not crash
    assert "record (none)" in r.stdout, r.stdout

    # 2. TBL_SUFFIX alone must refuse and name all five
    r2 = run({"TBL_SUFFIX": "_rl125"})
    assert r2.returncode != 0, r2.stdout + r2.stderr
    b2 = r2.stdout + r2.stderr
    for k in ("REC_SE", "REC_SE_POOLED", "REC_ALLTYPES", "REC_BORDERLINE",
              "REC_CLAIM"):
        assert k in b2, (k, b2)
    assert "`none`" in b2 and "NOT a neutral default" in b2, b2

    full = {"TBL_SUFFIX": "_rl125",
            "REC_SE": "Colorado=311/385,LBNL=122/156,Brigham=298/364,"
                      "Fuwai=422/596,Stanford=345/455",
            "REC_SE_POOLED": "1376,1800",
            "REC_ALLTYPES": "1883,2505",
            "REC_BORDERLINE": "none",
            "REC_CLAIM": "none"}

    # 3. the 125 bp run -> BOTH gates must say NOT RUN, and G4 must NOT
    #    announce TWO-SIDED from the empty record
    r3 = run(full)
    assert r3.returncode == 0, r3.stdout + r3.stderr
    assert "COMB cohort_type_concordance_dedup_rl125.tsv" in r3.stdout, r3.stdout
    assert "SE1 direction_concordance_dedup_rl125.tsv" in r3.stdout, r3.stdout
    assert "ALLTYPES (1883, 2505) CLAIM (0, 0)" in r3.stdout, r3.stdout
    assert "SEPOOL (1376, 1800) SE_COL (311, 385)" in r3.stdout, r3.stdout
    assert "G2 NOT RUN" in r3.stdout and "G2 PASS" not in r3.stdout, r3.stdout
    assert "G4 NOT RUN" in r3.stdout, r3.stdout
    assert "TWO-SIDED" not in r3.stdout, r3.stdout        # <-- the whole point
    assert "p-values are ONE-SIDED" in r3.stdout, r3.stdout
    assert "total absolute deviation" not in r3.stdout, r3.stdout
    assert "p-values? (none recorded)" in r3.stdout, r3.stdout
    assert "The SE headline (1,376/1,800 = 76.4%)" in r3.stdout, r3.stdout
    assert "POOLED SE only (the headline) (1376, 1800)" in r3.stdout, r3.stdout
    # NB not a bare "BORD" — that substring also occurs in the guard's own
    # "REC_BORDERLINE is empty" message, which is output, not a loop iteration.
    assert "BORD Brigham" not in r3.stdout, r3.stdout
    assert "BORD Fuwai" not in r3.stdout, r3.stdout

    # 4. a REAL 125 bp borderline set turns both gates back on
    r4 = run(dict(full, REC_BORDERLINE="Brigham:A5SS=19/29@0.07"))
    assert r4.returncode == 0, r4.stdout + r4.stderr
    assert "G2 PASS" in r4.stdout and "G2 NOT RUN" not in r4.stdout, r4.stdout
    assert "G4 NOT RUN" not in r4.stdout, r4.stdout
    assert "recorded 0.07?" in r4.stdout, r4.stdout
    assert "BORD Brigham A5SS 19 29 0.07" in r4.stdout, r4.stdout
    assert "p-values are ONE-SIDED" in r4.stdout, r4.stdout

    # 5. a wrong REC_SE key set fails loudly and names the expected cohorts
    r5 = run(dict(full, REC_SE="Colorado=311/385,Brigham=298/364"))
    assert r5.returncode != 0, r5.stdout
    b5 = r5.stdout + r5.stderr
    assert "must name exactly" in b5 and "Stanford" in b5 and "LBNL" in b5, b5

    # 6. REC_CLAIM may still be given a real pair
    r6 = run(dict(full, REC_CLAIM="15,18"))
    assert r6.returncode == 0, r6.stdout + r6.stderr
    assert "CLAIM (15, 18)" in r6.stdout, r6.stdout

    # 7. double application refused, text unchanged
    n7, p7 = apply_edits(new)
    assert p7, "second application was not refused"
    assert n7 == new, "a refused re-application still altered the text"

    print("[ok] self-test passed: unsuffixed output is unchanged including "
          '"recorded 0.07 / 0.08 / 0.14?" and "(1,636/2,163 = 75.6%)"; '
          "TBL_SUFFIX alone refuses and names all five records and the `none` "
          "option; with REC_BORDERLINE=none, G2 prints NOT RUN instead of PASS "
          "and G4 prints NOT RUN instead of announcing TWO-SIDED from an empty "
          "record; a real borderline set turns both gates back on; a cohort "
          "missing from REC_SE prints 'record (none)' instead of crashing; a "
          "wrong REC_SE key set names the expected cohorts; 1,636/2,163 now "
          "has one source; double application refused")
    return 0


if __name__ == "__main__":
    sys.exit(main())
