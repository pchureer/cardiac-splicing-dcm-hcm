#!/usr/bin/env python3
"""
132_patch_98_table4.py -- Table 4 at 125 bp. Seven stale constants, two
duplicated literals, and one latent collapse-order defect.

WHY THE OVERRIDES ARE MANDATORY, NOT DEFAULTED
----------------------------------------------
`RECORD`, `POOLED_K/N`, `POOLED_P_1SF` and `MEAN_SLOPE` are not inputs — 98's
own header calls each one "a gate, not an input". They exist to be checked
against the files. So leaving a 50 bp record in place while reading a 125 bp
file does not produce a wrong table: it makes G2/G3/G4 **refuse a correct
one**. That is a safe failure, but a confusing and expensive one, so when
`TBL_SUFFIX` is set this patch requires every record value explicitly and
refuses with the full command line if any is missing. Silence is not a default.

WHAT MOVES
----------
  84/85  OUT_TSV / OUT_MD          Table4.tsv / .md          -> suffixed
  89-95  RECORD                    k, n_ev, deming x4        -> REC_KN, REC_SLOPE
  96     POOLED_K, POOLED_N        1636, 2163                -> REC_POOLED
  97     POOLED_P_1SF              8e-132                    -> REC_POOLED_P
  98     MEAN_SLOPE                0.406                     -> REC_MEAN_SLOPE
  102/3  CAND_CONC / CAND_SLOPE    input candidates          -> suffixed + REORDERED
  ~243   the G5 note               "Its 0.445 / 75.9%"       -> REC_LBNL
  ~270   the G2 provenance label   "recorded"/"DERIVED"      -> SRC_LABEL
  ~289   the G3 print              "record 8.4e-132"         -> POOLED_P_TEXT
  ~294   the G3 refusal            "recorded 8.4e-132"       -> POOLED_P_TEXT
  ~455   the closing advisory      "0.36-0.59", filename     -> both named

`8.4e-132` appears THREE times: as `POOLED_P_1SF` (rounded to 8e-132, which is
all the factor-of-two gate needs) and twice as display text at full precision.
Changing the constant alone would leave "record 8.4e-132" printed beside a
different number — Rule 8. `REC_POOLED_P` now drives both: its string is the
display text and its float is the gate value, so the 50 bp display stays
byte-identical at "8.4e-132" while the constant stays 8e-132.

THE LATENT DEFECT IN CAND_SLOPE
-------------------------------
    CAND_CONC  = [..._dedup.tsv, ....tsv]        dedup FIRST
    CAND_SLOPE = [cohort_slopes.tsv, ..._dedup.tsv, ...]   dedup SECOND

If `cohort_slopes.tsv` ever existed, 98 would read **slopes filter-first and
concordances collapse-first** — the 926-vs-991 error, silently, inside one
table. Neither non-dedup file exists today (checked 16 Oct), so this is latent,
not active. Both lists now put `_dedup` first.

THE PROVENANCE LABEL ACTUALLY CHANGES STATUS
--------------------------------------------
98 prints `src = "recorded" if coh == "Colorado" else "DERIVED"`, and its G2
refusal warns that Brigham's, Stanford's and Fuwai's numerators "was derived
from a rounded percentage". **At 125 bp all four numerators are counted** — 33
section D measured them and 45 and 61 recounted them independently. So the
label becomes "measured" for all four when suffixed, rather than carrying a
caveat that no longer applies.

MEASURED VALUES (33 section D; 45 and 61 independently; slopes from 45's own
deming(), which agrees with 33's to 4 dp)

    Colorado 311/385  0.5999      Brigham  298/364  0.5099
    Stanford 345/455  0.4008      Fuwai    422/596  0.2948
    pooled   1376/1800 = 76.4%    P = 2.51e-117     mean slope 0.451
    LBNL     122/156 = 78.2%      0.5726   (not independent)

`POOLED_P` recomputed from k and n (one-sided exact binomial, p = 0.5), method
confirmed by reproducing the recorded 8.4e-132 from 1636/2163 exactly.
`MEAN_SLOPE` is 0.451, NOT the 0.452 I first predicted:
(0.5999+0.5099+0.4008+0.2948)/4 = 0.45135.

USAGE
-----
  python 132_patch_98_table4.py
  python 132_patch_98_table4.py --apply
  python 132_patch_98_table4.py --self-test
"""
import argparse
import difflib
import os
import shutil
import sys
import time

TARGET = "98_build_replication_table.py"

RUN125 = (
    "TBL_SUFFIX=_rl125 \\\n"
    "REC_KN=Colorado=311/385,Brigham=298/364,Stanford=345/455,Fuwai=422/596 \\\n"
    "REC_SLOPE=Colorado=0.5999,Brigham=0.5099,Stanford=0.4008,Fuwai=0.2948 \\\n"
    "REC_POOLED=1376,1800 \\\n"
    "REC_POOLED_P=2.51e-117 \\\n"
    "REC_MEAN_SLOPE=0.451 \\\n"
    "REC_LBNL=0.5726/78.2 \\\n"
    "python3 98_build_replication_table.py")

# ------------------------------------------------------------------- edit 1
A_OUT = 'OUT_TSV = "Table4.tsv"\nOUT_MD = "Table4.md"\n'
B_OUT = '''# Read length. Drives the output names, the input candidates, and the record.
SUF = os.environ.get("TBL_SUFFIX", "")
OUT_TSV = f"Table4{SUF}.tsv"
OUT_MD = f"Table4{SUF}.md"
'''

# ------------------------------------------------------------------- edit 2
A_MEAN = "MEAN_SLOPE = 0.406\n"
B_MEAN = A_MEAN + '''
# --------------------------------------------------------------------------
# 125 bp RECORD, placed AFTER every literal it rewrites.
#
# These are gates, not inputs (98's own header says so). A 50 bp record against
# a 125 bp file does not make a wrong table — it makes G2/G3/G4 refuse a
# CORRECT one. So when TBL_SUFFIX is set, every value is required explicitly.
_REQ = ("REC_KN", "REC_SLOPE", "REC_POOLED", "REC_POOLED_P",
        "REC_MEAN_SLOPE", "REC_LBNL")
_env = {k: os.environ.get(k, "").strip() for k in _REQ}
if SUF:
    _miss = [k for k in _REQ if not _env[k]]
    if _miss:
        raise SystemExit(
            "REFUSING — TBL_SUFFIX=%s but these are unset: %s\\n\\n"
            "  RECORD is an independent gate. A 50 bp record read against a\\n"
            "  125 bp file refuses a CORRECT table. Supply all six:\\n\\n"
            "    %s\\n" % (SUF, ", ".join(_miss), @@RUN125@@))

if _env["REC_KN"] or _env["REC_SLOPE"]:
    def _parse(s):
        return dict(p.split("=", 1) for p in s.split(",")) if s else {}
    _kn, _sl = _parse(_env["REC_KN"]), _parse(_env["REC_SLOPE"])
    for _nm, _d in (("REC_KN", _kn), ("REC_SLOPE", _sl)):
        if _d and set(_d) != set(RECORD):
            raise SystemExit(
                "REFUSING — %s must name exactly %s, got %s"
                % (_nm, sorted(RECORD), sorted(_d)))
    for _c in list(RECORD):
        _contrast, _acc, _n, _k, _ne, _s = RECORD[_c]
        if _c in _kn:
            _a, _b = _kn[_c].split("/")
            _k, _ne = int(_a), int(_b)
        if _c in _sl:
            _s = float(_sl[_c])
        RECORD[_c] = (_contrast, _acc, _n, _k, _ne, _s)

if _env["REC_POOLED"]:
    POOLED_K, POOLED_N = (int(_x) for _x in _env["REC_POOLED"].split(","))

# ONE variable drives both the gate value and the two display sites. The
# recorded text is more precise than the 1-s.f. constant the factor-of-two
# gate needs, so deriving the text from the float would silently change the
# 50 bp output from "8.4e-132" to "8.0e-132".
POOLED_P_TEXT = _env["REC_POOLED_P"] or "8.4e-132"
if _env["REC_POOLED_P"]:
    POOLED_P_1SF = float(_env["REC_POOLED_P"])

if _env["REC_MEAN_SLOPE"]:
    MEAN_SLOPE = float(_env["REC_MEAN_SLOPE"])

# LBNL's recorded slope / concordance, quoted in the G5 note. Not a gate.
LBNL_SLOPE, LBNL_PCT = (_env["REC_LBNL"] or "0.445/75.9").split("/")

# At 50 bp only Colorado's numerator was recorded; Brigham, Stanford and Fuwai
# were back-derived from a rounded percentage. At 125 bp all four are COUNTED
# (33 section D, recounted independently by 45 and by 61), so the caveat does
# not apply and the label should not carry it.
SRC_LABEL = {_c: ("measured" if SUF else
                  ("recorded" if _c == "Colorado" else "DERIVED"))
             for _c in RECORD}
# --------------------------------------------------------------------------
'''

# The run command is embedded in the refusal as a Python literal, so the
# message 98 prints and the command in this patch's dry run can never drift.
# B_MEAN is a plain ''' string, so this must be a real substitution — an
# inline `""" + repr(...) + """` would have been literal text (it was, and the
# self-test printed it verbatim).
B_MEAN = B_MEAN.replace("@@RUN125@@", repr(RUN125.replace("\n", "\n    ")))
assert "@@RUN125@@" not in B_MEAN and "TBL_SUFFIX=_rl125" in B_MEAN

# ------------------------------------------------------------------- edit 3
A_CAND = ('CAND_CONC = ["cohort_type_concordance_dedup.tsv", '
          '"cohort_type_concordance.tsv"]\n'
          'CAND_SLOPE = ["cohort_slopes.tsv", "cohort_slopes_dedup.tsv", '
          '"Figure4_slopes.tsv"]\n')
B_CAND = '''CAND_CONC = [f"cohort_type_concordance_dedup{SUF}.tsv",
             f"cohort_type_concordance{SUF}.tsv"]
# _dedup FIRST IN BOTH LISTS. Before this patch CAND_SLOPE preferred the
# NON-dedup `cohort_slopes.tsv`, so had that file existed, 98 would have read
# slopes filter-first and concordances collapse-first — the 926-vs-991 error,
# silently, inside one table. Neither non-dedup file exists (checked 16 Oct),
# so the defect was latent; the order now makes it unreachable.
CAND_SLOPE = [f"cohort_slopes_dedup{SUF}.tsv",
              f"cohort_slopes{SUF}.tsv",
              f"Figure4_slopes{SUF}.tsv"]
'''

EDITS = [
    ("SUF and the suffixed output names", A_OUT, B_OUT, 1),
    ("the 125 bp record block", A_MEAN, B_MEAN, 1),
    ("the input candidates, suffixed and reordered", A_CAND, B_CAND, 1),
    ("the G3 print, from POOLED_P_TEXT",
     'print(f"      one-sided exact binomial P = {pv:.1e}   record 8.4e-132")',
     'print(f"      one-sided exact binomial P = {pv:.1e}   '
     'record {POOLED_P_TEXT}")', 1),
    ("the G3 refusal, from POOLED_P_TEXT",
     'f"recorded 8.4e-132.',
     'f"recorded {POOLED_P_TEXT}.', 1),
    ("the G5 LBNL note, from REC_LBNL",
     'print("      (LBNL is a library-prep replicate of discovery. '
     'Its 0.445 / 75.9% was")',
     'print(f"      (LBNL is a library-prep replicate of discovery. '
     'Its {LBNL_SLOPE} / {LBNL_PCT}% was")', 1),
    ("the G2 provenance label",
     'src = "recorded" if coh == "Colorado" else "DERIVED"',
     'src = SRC_LABEL.get(coh, "recorded")', 1),
    # prose inside a NON-f triple-quoted block: name both read lengths rather
    # than convert the block to an f-string, which would break on any brace.
    # TWO occurrences, found by the anchor guard refusing on a count of 2:
    # line 206 (a plain print in G1) and line 457 (the closing advisory). Both
    # make the same claim, so both move, and both get the literal text rather
    # than an f-string — line 457 sits inside a non-f triple-quoted block, and
    # naming both read lengths explicitly cannot be wrong about which is which.
    ("the recorded r range, both occurrences",
     "the range 0.36-0.59",
     "the range 0.36-0.59 at 50 bp / 0.41-0.60 at 125 bp", 2),
    ("the closing advisory filename",
     "Counted from `cohort_type_concordance_dedup.tsv` by",
     "Counted from `cohort_type_concordance_dedup[_rl125].tsv` by", 1),
    ("the G2 refusal caveat, dated",
     "the recorded NUMERATOR was derived ",
     "the 50 bp recorded NUMERATOR was derived ", 1),
    # line 97's own comment documents the constant and must point at the new
    # mechanism — and note the "stale site" checks below deliberately match
    # CODE forms only (`record 8.4e-132")`, `f"recorded 8.4e-132.`) so that
    # this comment cannot make a correct patch look like a failed one.
    ("the POOLED_P_1SF comment",
     "POOLED_P_1SF = 8e-132          # recorded 8.4e-132",
     "POOLED_P_1SF = 8e-132          # 50 bp; display text and the gate value\n"
     "                               # both come from REC_POOLED_P below", 1),
]


def apply_edits(text):
    if 'SUF = os.environ.get("TBL_SUFFIX", "")' in text:
        return text, ["SUF is already defined — 98 is already patched "
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
              'OUT_TSV = f"Table4{SUF}.tsv"',
              'OUT_MD = f"Table4{SUF}.md"',
              '_REQ = ("REC_KN", "REC_SLOPE", "REC_POOLED", "REC_POOLED_P",',
              "POOLED_P_TEXT = _env",
              "LBNL_SLOPE, LBNL_PCT = ",
              "SRC_LABEL = {_c: ",
              'f"cohort_slopes_dedup{SUF}.tsv"',
              'f"cohort_type_concordance_dedup{SUF}.tsv"',
              "record {POOLED_P_TEXT}",
              "recorded {POOLED_P_TEXT}",
              "Its {LBNL_SLOPE} / {LBNL_PCT}%",
              'src = SRC_LABEL.get(coh, "recorded")'):
        if w not in text:
            bad.append(f"missing: {w[:58]}")
    # CODE forms only. A bare "recorded 8.4e-132" also occurs in line 97's
    # comment, which this patch rewrites but which must not be able to make a
    # correct patch report a stale site.
    for gone in ('OUT_TSV = "Table4.tsv"',
                 '"cohort_slopes.tsv", "cohort_slopes_dedup.tsv"',
                 'record 8.4e-132")',
                 'f"recorded 8.4e-132.',
                 "Its 0.445 / 75.9%",
                 'if coh == "Colorado" else "DERIVED"'):
        if gone in text:
            bad.append(f"stale site survives: {gone[:50]}")
    # _dedup must precede the non-dedup candidate in BOTH lists
    i_cs = text.find("CAND_SLOPE = [")
    if i_cs >= 0:
        seg = text[i_cs:i_cs + 260]
        if seg.find("cohort_slopes_dedup") > seg.find('f"cohort_slopes{SUF}'):
            bad.append("CAND_SLOPE still prefers the non-dedup file")
    # SUF module-level, above def main()
    i_s = text.find('SUF = os.environ.get("TBL_SUFFIX", "")')
    i_m = text.find("def main(")
    if i_s < 0 or i_m < 0 or i_s > i_m:
        bad.append("SUF is not module-level / not above def main()")
    if "\n    SUF = os.environ.get" in text:
        bad.append("SUF is indented — it would be a local, not a global")
    # the record block must come after every literal it rewrites
    for lit, nm in (("MEAN_SLOPE = 0.406", "MEAN_SLOPE"),
                    ("POOLED_K, POOLED_N = 1636, 2163", "POOLED_K/N"),
                    ("POOLED_P_1SF = 8e-132", "POOLED_P_1SF")):
        i_l = text.find(lit)
        i_o = text.find("_env = {k: os.environ.get(k")
        if i_l < 0 or i_o < 0 or i_l > i_o:
            bad.append(f"the override block is not after the {nm} literal")
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
    print("132 — 98_build_replication_table.py, Table 4 at 125 bp")
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
        print("  # 1. CALIBRATION — unsuffixed must still reproduce Table 4")
        print("  python3 98_build_replication_table.py --list 2>&1 | tail -40")
        print("  python3 98_build_replication_table.py 2>&1 | tail -40")
        print()
        print("  # 2. the refusal must fire when the record is not supplied")
        print("  TBL_SUFFIX=_rl125 python3 98_build_replication_table.py "
              "2>&1 | head -12")
        print()
        print("  # 3. 125 bp")
        print("  " + RUN125.replace("\n", "\n  "))
        print()
        print("  # 4. then --apply to write Table4_rl125.tsv / .md")
        print()
        print("⚠ G2/G3/G4 are gates on values I supplied, so a PASS means the")
        print("  files and the record agree — it does not independently confirm")
        print("  either. The independent confirmations already happened: 33, 45")
        print("  and 61 each counted 311/385 etc., and 45's own deming() agrees")
        print("  with 33's to 4 dp. What G2-G4 add is that Table 4 cannot be")
        print("  built from a file that disagrees with them.")
        print()
        print("⚠ POOLED_P 2.51e-117 is recomputed from k and n, one-sided exact")
        print("  binomial at p=0.5 — the same code reproduces the recorded")
        print("  8.4e-132 from 1636/2163 exactly, so the method is verified.")
        print()
        print("⬜ STILL STALE AFTER THIS: 85 (REC_SE, REC_POOLED_ALLTYPES ->")
        print("  1883/2505, REC_BORDERLINE, REC_CLAIM, an inline rec=(1636,")
        print("  2163) at line 328 and the prose at 331), 82 (EXPECT_OBSERVED_C")
        print("  0.521 -> 0.600, PER_COHORT), 44 (KNOWN_POOLED line 61), and")
        print("  the slope literals in 26, 35, 42 (Table 1!), 67 and 83.")
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
    print("✓ verified on disk: the record is overridable and REQUIRED when")
    print("  suffixed; both display copies of the recorded P follow one")
    print("  variable; both candidate lists prefer _dedup; the provenance")
    print("  label reflects that 125 bp numerators are counted, not derived.")
    return 0


def self_test():
    import subprocess
    import tempfile

    fx = (
        "import argparse, os\n"
        + A_OUT
        + 'RECORD = {\n'
          '    "Colorado": ("DCM vs NF",   "PRJNA477855", "36/13", 375,  471, 0.521),\n'
          '    "Brigham":  ("HCM vs ctrl", "PRJNA675090", "18/5",  359,  442, 0.476),\n'
          '    "Stanford": ("HCM vs ctrl", "PRJNA747806", "13/6",  408,  536, 0.386),\n'
          '    "Fuwai":    ("HCM vs ctrl", "PRJNA522931", "28/9",  494,  714, 0.243),\n'
          '}\n'
          "POOLED_K, POOLED_N = 1636, 2163\n"
          "POOLED_P_1SF = 8e-132          # recorded 8.4e-132\n"
        + A_MEAN
        + "SLOPE_TOL = 0.0015\n"
        + A_CAND
        + "def main():\n"
          "    pv = 1e-120\n"
          '    print("OUT", OUT_TSV, OUT_MD)\n'
          '    print("CONC", CAND_CONC[0], "SLOPE", CAND_SLOPE[0])\n'
          '    print("POOLED", POOLED_K, POOLED_N, "MEAN", MEAN_SLOPE)\n'
          '    print("KN", RECORD["Colorado"][3], RECORD["Colorado"][4],\n'
          '          "SL", RECORD["Fuwai"][5])\n'
          '    print(f"      one-sided exact binomial P = {pv:.1e}   '
          'record 8.4e-132")\n'
          '    print("      (LBNL is a library-prep replicate of discovery. '
          'Its 0.445 / 75.9% was")\n'
          "    for coh in RECORD:\n"
          '        src = "recorded" if coh == "Colorado" else "DERIVED"\n'
          '        print("SRC", coh, src)\n'
          "    if pv < 0:\n"
          '        raise SystemExit(f"recorded 8.4e-132. tail")\n'
          '    print("     open item that the record held only as the range '
          '0.36-0.59.")\n'
          '    print("""x\n'
          "  2. the range 0.36-0.59; check\n"
          "     Counted from `cohort_type_concordance_dedup.tsv` by x\n"
          '  the recorded NUMERATOR was derived from a rounded pct\n'
          '""")\n'
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
        for k in ("TBL_SUFFIX", "REC_KN", "REC_SLOPE", "REC_POOLED",
                  "REC_POOLED_P", "REC_MEAN_SLOPE", "REC_LBNL"):
            e.pop(k, None)
        e.update(env or {})
        return subprocess.run([sys.executable, "-I", p], capture_output=True,
                              text=True, cwd=d, env=e)

    # 1. unsuffixed -> every original value and every original display string
    r = run()
    assert r.returncode == 0, r.stdout + r.stderr
    assert "OUT Table4.tsv Table4.md" in r.stdout, r.stdout
    assert "POOLED 1636 2163 MEAN 0.406" in r.stdout, r.stdout
    assert "KN 375 471 SL 0.243" in r.stdout, r.stdout
    assert "record 8.4e-132" in r.stdout, r.stdout          # NOT 8.0e-132
    assert "Its 0.445 / 75.9%" in r.stdout, r.stdout
    assert "SRC Colorado recorded" in r.stdout, r.stdout
    assert "SRC Fuwai DERIVED" in r.stdout, r.stdout
    assert "CONC cohort_type_concordance_dedup.tsv" in r.stdout, r.stdout
    assert "SLOPE cohort_slopes_dedup.tsv" in r.stdout, r.stdout   # REORDERED
    assert r.stdout.count("0.36-0.59 at 50 bp / 0.41-0.60 at 125 bp") == 2, \
        r.stdout
    assert "the range 0.36-0.59." not in r.stdout, r.stdout
    assert "the 50 bp recorded NUMERATOR" in r.stdout, r.stdout

    # 2. TBL_SUFFIX alone must REFUSE and name every missing variable
    r2 = run({"TBL_SUFFIX": "_rl125"})
    assert r2.returncode != 0, r2.stdout + r2.stderr
    b2 = r2.stdout + r2.stderr
    assert "REFUSING" in b2, b2
    for k in ("REC_KN", "REC_SLOPE", "REC_POOLED", "REC_POOLED_P",
              "REC_MEAN_SLOPE", "REC_LBNL"):
        assert k in b2, (k, b2)
    assert "311/385" in b2 and "2.51e-117" in b2 and "0.451" in b2, b2

    # 3. a partial record must still refuse, naming only what is missing
    r3 = run({"TBL_SUFFIX": "_rl125", "REC_POOLED": "1376,1800"})
    assert r3.returncode != 0, r3.stdout
    assert "REC_KN" in (r3.stdout + r3.stderr), r3.stdout + r3.stderr

    full = {"TBL_SUFFIX": "_rl125",
            "REC_KN": "Colorado=311/385,Brigham=298/364,"
                      "Stanford=345/455,Fuwai=422/596",
            "REC_SLOPE": "Colorado=0.5999,Brigham=0.5099,"
                         "Stanford=0.4008,Fuwai=0.2948",
            "REC_POOLED": "1376,1800",
            "REC_POOLED_P": "2.51e-117",
            "REC_MEAN_SLOPE": "0.451",
            "REC_LBNL": "0.5726/78.2"}

    # 4. the full record -> every value and every display string moves
    r4 = run(full)
    assert r4.returncode == 0, r4.stdout + r4.stderr
    assert "OUT Table4_rl125.tsv Table4_rl125.md" in r4.stdout, r4.stdout
    assert "POOLED 1376 1800 MEAN 0.451" in r4.stdout, r4.stdout
    assert "KN 311 385 SL 0.2948" in r4.stdout, r4.stdout
    assert "record 2.51e-117" in r4.stdout, r4.stdout
    assert "Its 0.5726 / 78.2%" in r4.stdout, r4.stdout
    assert "CONC cohort_type_concordance_dedup_rl125.tsv" in r4.stdout, r4.stdout
    assert "SLOPE cohort_slopes_dedup_rl125.tsv" in r4.stdout, r4.stdout
    # the provenance caveat must be gone for ALL FOUR cohorts
    assert "DERIVED" not in r4.stdout, r4.stdout
    for c in ("Colorado", "Brigham", "Stanford", "Fuwai"):
        assert f"SRC {c} measured" in r4.stdout, r4.stdout

    # 5. a misspelt / short cohort set fails loudly and names the expected keys
    bad = dict(full, REC_KN="Colorado=311/385,Brigham=298/364")
    r5 = run(bad)
    assert r5.returncode != 0, r5.stdout
    b5 = r5.stdout + r5.stderr
    assert "must name exactly" in b5 and "Stanford" in b5 and "Fuwai" in b5, b5
    bad2 = dict(full, REC_SLOPE="Colorado=0.6,Brigham=0.5,"
                                "Stanford=0.4,Fuwaii=0.3")
    r6 = run(bad2)
    assert r6.returncode != 0, r6.stdout
    assert "must name exactly" in (r6.stdout + r6.stderr), r6.stdout

    # 6. the override may be used WITHOUT a suffix (a 50 bp re-record), and
    #    then only what is given moves
    r7 = run({"REC_MEAN_SLOPE": "0.451"})
    assert r7.returncode == 0, r7.stdout + r7.stderr
    assert "POOLED 1636 2163 MEAN 0.451" in r7.stdout, r7.stdout
    assert "OUT Table4.tsv" in r7.stdout, r7.stdout

    # 7. double application refused, text unchanged
    n8, p8 = apply_edits(new)
    assert p8, "second application was not refused"
    assert n8 == new, "a refused re-application still altered the text"

    print("[ok] self-test passed: unsuffixed output is byte-unchanged down to "
          '"record 8.4e-132" and "Its 0.445 / 75.9%"; TBL_SUFFIX alone REFUSES '
          "and names all six missing variables plus the measured values; a "
          "partial or misspelt record refuses and names the expected cohorts; "
          "the full record moves every constant, both output names, both "
          "candidate lists and both display strings; CAND_SLOPE now prefers "
          "_dedup; the provenance label becomes 'measured' for all four "
          "cohorts at 125 bp; overrides work without a suffix; double "
          "application refused")
    return 0


if __name__ == "__main__":
    sys.exit(main())
