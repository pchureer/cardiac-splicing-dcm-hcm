#!/usr/bin/env python3
"""
138_patch_82_s4.py -- Supplementary Table 4 at 125 bp. THIRTEEN constants.

82 is the last blocked deliverable in the 45-consumer chain
(45 -> 98 -> 85 -> 136 -> 22 -> **82**).

WHY THE OVERRIDES ARE MANDATORY
-------------------------------
82 is a reproduce-or-refuse build, by its own design. Three of its gates DIE
rather than warn:

  G1  searches selection-set definitions for the one giving EXPECT_COMMON
      events, and dies if none does
  G2  compares direction agreement against EXPECT_DIR_A / EXPECT_DIR_B
  G3  searches lambda conventions x arm directions for the one reproducing
      BOTH EXPECT_A and EXPECT_B within --tol (default 0.01), and dies if none
      does

So a 50 bp record read against the 125 bp `mp_SE_*_rl125` tables does not
produce a wrong table — it **refuses a correct one**, and the refusal arrives
as "no selection-set definition gives 203 events" or "no convention
reproduces both 0.584 and 0.988", which reads like a data problem rather than
a stale constant. Hence: when TBL_SUFFIX is set, every record value is
required and the refusal prints the whole command line.

THIRTEEN CONSTANTS, GROUPED INTO THREE VARIABLES
------------------------------------------------
Thirteen separate env vars would be unusable, so they are grouped. Each group
asserts SET EQUALITY against the keys already in the file, so a typo or an
omission fails loudly and the error names what was expected.

  line  71  OUT = "S4.tsv"                                  -> suffixed
  lines 73-78  ARMS, four mp_SE_* tables                     -> suffixed
  line  81  EXPECT_COMMON   203   -> 161      S4_RECORD common=
  line  82  EXPECT_A      0.584   -> 0.690              ... A=
  line  83  EXPECT_B      0.988   -> 0.968              ... B=
  line  84  EXPECT_A_CI  (0.478, 0.693) -> (0.570, 0.806)  ... A_CI=
  line  85  EXPECT_B_CI  (0.889, 1.089) -> (0.875, 1.069)  ... B_CI=
  line  86  EXPECT_DIR_A    162   -> 130                ... dirA=
  line  87  EXPECT_DIR_B    174   -> 140                ... dirB=
  line  88  EXPECT_PRODUCT 0.577  -> 0.668              ... product=
  line  89  EXPECT_OBSERVED_C 0.521 -> 0.600            ... C=
  line  91  PER_COHORT                                  S4_PER_COHORT
  line  93  GENOME_WIDE                                 S4_GENOME_WIDE
  lines 249-250  the --list record print                 derived from the above

`EXPECT_COMMON` was nearly missed: it does not contain a slope and does not
look like a read-length quantity, but G1 dies without it.

MEASURED VALUES AND THEIR SOURCES
---------------------------------
  A, B, CIs, DIR_A, DIR_B   36_matched_pair_analysis.py at 125 bp
                            (130/161 and 140/161 on the common set)
  COMMON 161                same run
  PRODUCT 0.668  C 0.600    A x B against 45/33's Colorado slope; the closure
                            gap is 11.3% against 10.7% at 50 bp
  PER_COHORT                45_cohort_slopes.py, confirmed against 33's
                            independent deming() to 4 dp
  GENOME_WIDE               22_symmetry_winners_curse.py FORWARD column at
                            125 bp, after 136 built the MAGNet-internal arm

⚠ GENOME_WIDE's Brigham AND Fuwai entries are REPLACED, NOT UPDATED.
The recorded 0.757 and 0.359 match **nothing** in any surviving symmetry
output, at any read length, in any column, for any event type. 22's 50 bp
output was verified byte-identical on 20 Oct and gives Brigham 0.662 and
Fuwai 0.328; Colorado 0.571, LBNL 0.525 and Stanford 0.578 all reproduce. The
run that produced the other two was overwritten before 22's outputs carried a
version. `62_build_s5.py:16` already flags Brigham 0.757 as contradicting
Figure 4 — it is not a genome-wide-vs-selected distinction, it is an
unreproducible number. See DEFECT_genome_wide_brigham_fuwai_unsourced.md.

⚠ Do NOT be reassured that `symmetry_SE_rl125.tsv` contains Brigham
unsel_deming 0.757203. GENOME_WIDE is a dict of FORWARD slopes; that is the
UNSELECTED column at the other read length. Two differences at once.

WATCH THE LAMBDA CONVENTION
---------------------------
82's line 26 warns the convention "is not recorded, so the script searches for
it". 121 measured lambda moving 5.37 -> 7.19 (+34%) between read lengths, so
G3 may select a DIFFERENT convention at 125 bp than at 50 bp. That is not a
bug — but if it happens, the two tables' slopes were computed under different
conventions and S4 must say so. The patch leaves the search alone and prints
the winner either way.

USAGE
-----
  python 138_patch_82_s4.py [--apply] [--self-test]
"""
import argparse
import difflib
import os
import shutil
import sys
import time

TARGET = "82_build_s4.py"

RUN125 = (
    "TBL_SUFFIX=_rl125 \\\n"
    "S4_RECORD=common=161,A=0.690,B=0.968,A_CI=0.570:0.806,"
    "B_CI=0.875:1.069,dirA=130,dirB=140,product=0.668,C=0.600 \\\n"
    "S4_PER_COHORT=Colorado=0.5999,Brigham=0.5099,"
    "Stanford=0.4008,Fuwai=0.2948 \\\n"
    "S4_GENOME_WIDE=Colorado=0.674,Brigham=0.722,"
    "Fuwai=0.365,Stanford=0.605 \\\n"
    "python3 82_build_s4.py")

# ------------------------------------------------------------------- edit 1
A_ARMS = '''OUT = "S4.tsv"

ARMS = {
    "selection": "mp_SE_selection.all_events.tsv",
    "magnet_matched": "mp_SE_magnet_matched.all_events.tsv",
    "lbnl_matched": "mp_SE_lbnl_matched.all_events.tsv",
    "magnet_other": "mp_SE_magnet_other.all_events.tsv",
}
'''
B_ARMS = '''# Read length. Drives the output name, all four input arms, and the record.
SUF = os.environ.get("TBL_SUFFIX", "")
OUT = f"S4{SUF}.tsv"

ARMS = {
    "selection": f"mp_SE_selection{SUF}.all_events.tsv",
    "magnet_matched": f"mp_SE_magnet_matched{SUF}.all_events.tsv",
    "lbnl_matched": f"mp_SE_lbnl_matched{SUF}.all_events.tsv",
    "magnet_other": f"mp_SE_magnet_other{SUF}.all_events.tsv",
}
'''

# ------------------------------------------------------------------- edit 2
A_REC = ('GENOME_WIDE = {"Colorado": 0.571, "Brigham": 0.757, '
         '"Fuwai": 0.359, "Stanford": 0.578}\n')
B_REC = A_REC + '''
# --------------------------------------------------------------------------
# 125 bp RECORD, placed AFTER every literal it rewrites.
#
# 82 is reproduce-or-refuse: G1 dies if no selection-set definition gives
# EXPECT_COMMON events, G2 compares EXPECT_DIR_*, and G3 dies if no lambda
# convention reproduces BOTH EXPECT_A and EXPECT_B. A 50 bp record against the
# 125 bp arms therefore REFUSES A CORRECT TABLE, with a message that reads
# like a data problem. So all three groups are REQUIRED when TBL_SUFFIX is set.
#
# Thirteen constants grouped into three variables; each asserts SET EQUALITY
# against the keys already present, so a typo names what was expected.
_S4_KEYS = ("common", "A", "B", "A_CI", "B_CI", "dirA", "dirB",
            "product", "C")
_S4_REQ = ("S4_RECORD", "S4_PER_COHORT", "S4_GENOME_WIDE")
_s4env = {_k: os.environ.get(_k, "").strip() for _k in _S4_REQ}
if SUF:
    _miss = [_k for _k in _S4_REQ if not _s4env[_k]]
    if _miss:
        raise SystemExit(
            "REFUSING — TBL_SUFFIX=%s but these are unset: %s\\n\\n"
            "  82's gates DIE rather than warn, so a 50 bp record against the\\n"
            "  125 bp arms refuses a CORRECT table — G1 would report 'no\\n"
            "  selection-set definition gives 203 events' and G3 'no\\n"
            "  convention reproduces both 0.584 and 0.988', neither of which\\n"
            "  is a data problem.\\n\\n    %s\\n"
            % (SUF, ", ".join(_miss), @@RUN125@@))

if _s4env["S4_RECORD"]:
    _d = dict(_p.split("=", 1) for _p in _s4env["S4_RECORD"].split(","))
    assert set(_d) == set(_S4_KEYS), (
        "S4_RECORD must name exactly %s, got %s"
        % (sorted(_S4_KEYS), sorted(_d)))
    EXPECT_COMMON = int(_d["common"])
    EXPECT_A = float(_d["A"])
    EXPECT_B = float(_d["B"])
    # CIs use ':' internally because ',' separates the entries
    EXPECT_A_CI = tuple(float(_x) for _x in _d["A_CI"].split(":"))
    EXPECT_B_CI = tuple(float(_x) for _x in _d["B_CI"].split(":"))
    EXPECT_DIR_A = int(_d["dirA"])
    EXPECT_DIR_B = int(_d["dirB"])
    EXPECT_PRODUCT = float(_d["product"])
    EXPECT_OBSERVED_C = float(_d["C"])
    assert len(EXPECT_A_CI) == 2 and len(EXPECT_B_CI) == 2, \\
        "A_CI and B_CI must each be lo:hi"

if _s4env["S4_PER_COHORT"]:
    _d = dict(_p.split("=", 1) for _p in _s4env["S4_PER_COHORT"].split(","))
    assert set(_d) == set(PER_COHORT), (
        "S4_PER_COHORT must name exactly %s, got %s"
        % (sorted(PER_COHORT), sorted(_d)))
    PER_COHORT = {_k: float(_v) for _k, _v in _d.items()}

# ⚠ Brigham and Fuwai are REPLACED, not updated. The recorded 0.757 and 0.359
# match nothing in any surviving symmetry output, at any read length, in any
# column, for any event type; 22's byte-identical 50 bp run gives 0.662 and
# 0.328. See DEFECT_genome_wide_brigham_fuwai_unsourced.md. And note that
# symmetry_SE_rl125.tsv's Brigham unsel_deming 0.757203 is NOT a confirmation:
# GENOME_WIDE holds FORWARD slopes, so that is the wrong column at the wrong
# read length.
if _s4env["S4_GENOME_WIDE"]:
    _d = dict(_p.split("=", 1) for _p in _s4env["S4_GENOME_WIDE"].split(","))
    assert set(_d) == set(GENOME_WIDE), (
        "S4_GENOME_WIDE must name exactly %s, got %s"
        % (sorted(GENOME_WIDE), sorted(_d)))
    GENOME_WIDE = {_k: float(_v) for _k, _v in _d.items()}


def _require_arms():
    """Refuse early if a suffixed arm is absent, rather than at G1."""
    if not SUF:
        return
    gone = [v for v in ARMS.values() if not os.path.exists(v)]
    if gone:
        raise SystemExit(
            "REFUSING — TBL_SUFFIX=%s but these arms are absent:\\n%s\\n"
            "  They come from the matched-pair build at 125 bp. Without them\\n"
            "  G1's selection-set search would fail on missing data and read\\n"
            "  as a stale-constant problem."
            % (SUF, "".join("    missing: %s\\n" % g for g in gone)))
# --------------------------------------------------------------------------
'''

# ------------------------------------------------------------------- edit 3
A_LIST = ('        print("\\n  The record: arm A 0.584 (0.478-0.693) '
          'platform only;")\n'
          '        print("              arm B 0.988 (0.889-1.089) '
          'people only;")\n')
B_LIST = ('        print(f"\\n  The record: arm A {EXPECT_A} '
          '({EXPECT_A_CI[0]}-{EXPECT_A_CI[1]}) platform only;")\n'
          '        print(f"              arm B {EXPECT_B} '
          '({EXPECT_B_CI[0]}-{EXPECT_B_CI[1]}) people only;")\n')

# ------------------------------------------------------------------- edit 4
A_DOC = """    arm A  MAGNet_matched -> LBNL_matched    platform only, same hearts   0.584 (0.478-0.693)
    arm B  MAGNet_matched -> MAGNet_other    people only, same platform   0.988 (0.889-1.089)
    A x B = 0.577, against an observed cross-cohort slope of 0.521
"""
B_DOC = """    arm A  MAGNet_matched -> LBNL_matched    platform only, same hearts   0.584 (0.478-0.693)
    arm B  MAGNet_matched -> MAGNet_other    people only, same platform   0.988 (0.889-1.089)
    A x B = 0.577, against an observed cross-cohort slope of 0.521

⚠ THOSE ARE 50 bp VALUES. At 125 bp, measured: arm A 0.690 (0.570-0.806),
arm B 0.968 (0.875-1.069), A x B = 0.668 against an observed C of 0.600 — a
closure gap of 11.3%, against 10.7% at 50 bp. The platform cost restates from
41% to 29%. Every EXPECT_ constant below is overridable and REQUIRED when
TBL_SUFFIX is set; do not read a read length off this docstring.
"""

EDITS = [
    ("SUF, the output name and the four arms", A_ARMS, B_ARMS, 1),
    ("the 125 bp record block", A_REC, B_REC, 1),
    ("the --list record print, derived", A_LIST, B_LIST, 1),
    ("the docstring's 50 bp values", A_DOC, B_DOC, 1),
    ("the G3 refusal note in the header",
     "reproduces 0.584 and 0.988**.",
     "reproduces EXPECT_A and EXPECT_B (0.584 and 0.988 at 50 bp)**.", 1),
    ("the convention-search explanation",
     "which one reproduces the recorded 0.584 and 0.988 "
     "**in both arms at once**.",
     "which one reproduces the recorded EXPECT_A and EXPECT_B "
     "**in both arms at once**.", 1),
    ("the G3 gate list entry",
     "G3  one lambda convention reproduces BOTH 0.584 and 0.988 "
     "within tolerance",
     "G3  one lambda convention reproduces BOTH recorded arm slopes "
     "within tolerance", 1),
]

B_REC = B_REC.replace("@@RUN125@@", repr(RUN125.replace("\n", "\n    ")))
assert "@@RUN125@@" not in B_REC and "TBL_SUFFIX=_rl125" in B_REC
EDITS[1] = ("the 125 bp record block", A_REC, B_REC, 1)


def apply_edits(text):
    if 'SUF = os.environ.get("TBL_SUFFIX", "")' in text:
        return text, ["SUF is already defined — 82 is already patched "
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
              'OUT = f"S4{SUF}.tsv"',
              'f"mp_SE_selection{SUF}.all_events.tsv"',
              'f"mp_SE_magnet_matched{SUF}.all_events.tsv"',
              'f"mp_SE_lbnl_matched{SUF}.all_events.tsv"',
              'f"mp_SE_magnet_other{SUF}.all_events.tsv"',
              '_S4_KEYS = ("common", "A", "B", "A_CI", "B_CI", "dirA", "dirB",',
              '_S4_REQ = ("S4_RECORD", "S4_PER_COHORT", "S4_GENOME_WIDE")',
              "EXPECT_COMMON = int(_d[\"common\"])",
              "EXPECT_OBSERVED_C = float(_d[\"C\"])",
              "def _require_arms():",
              "The record: arm A {EXPECT_A}",
              "THOSE ARE 50 bp VALUES",
              "REPLACED, not updated"):
        if w not in text:
            bad.append(f"missing: {w[:58]}")
    for gone in ('OUT = "S4.tsv"',
                 '"mp_SE_selection.all_events.tsv"',
                 '"mp_SE_magnet_other.all_events.tsv"',
                 'arm A 0.584 (0.478-0.693) platform only;")',
                 "reproduces 0.584 and 0.988**.",
                 "BOTH 0.584 and 0.988 within tolerance"):
        if gone in text:
            bad.append(f"stale site survives: {gone[:50]}")
    # SUF module-level, above def main() and above every literal it feeds
    i_s = text.find('SUF = os.environ.get("TBL_SUFFIX", "")')
    i_m = text.find("def main():")
    if i_s < 0 or i_m < 0 or i_s > i_m:
        bad.append("SUF is not module-level / not above def main()")
    if "\n    SUF = os.environ.get" in text:
        bad.append("SUF is indented — it would be a local, not a global")
    for m in ('OUT = f"S4{SUF}.tsv"', 'f"mp_SE_selection{SUF}'):
        j = text.find(m)
        if j < 0 or i_s > j:
            bad.append(f"SUF not defined above: {m[:40]}")
    # the record block must follow every literal it rewrites
    i_o = text.find("_s4env = {_k: os.environ.get(_k")
    for lit, nm in (("EXPECT_COMMON = 203", "EXPECT_COMMON"),
                    ("EXPECT_A = 0.584", "EXPECT_A"),
                    ("EXPECT_B = 0.988", "EXPECT_B"),
                    ("EXPECT_A_CI = (0.478, 0.693)", "EXPECT_A_CI"),
                    ("EXPECT_DIR_A = 162", "EXPECT_DIR_A"),
                    ("EXPECT_PRODUCT = 0.577", "EXPECT_PRODUCT"),
                    ("EXPECT_OBSERVED_C = 0.521", "EXPECT_OBSERVED_C"),
                    ("PER_COHORT = {", "PER_COHORT"),
                    ("GENOME_WIDE = {", "GENOME_WIDE")):
        i_l = text.find(lit)
        if i_l < 0 or i_o < 0 or i_l > i_o:
            bad.append(f"the record block is not after the {nm} literal")
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
    print("138 — 82_build_s4.py: Supplementary Table 4, thirteen constants")
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
        print("\n  Nothing written. Send me lines 8-50, 68-95 and 245-255.")
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
        print("  # 1. CALIBRATION — unsuffixed must still reproduce S4.tsv")
        print("  cp S4.tsv S4.PRESERVED 2>/dev/null")
        print("  python3 82_build_s4.py --list 2>&1 | tail -20")
        print("  python3 82_build_s4.py 2>&1 | tail -40")
        print("  cmp S4.tsv S4.PRESERVED && echo '50bp BYTE-IDENTICAL'")
        print()
        print("  # 2. the refusal must name all three groups")
        print("  TBL_SUFFIX=_rl125 python3 82_build_s4.py 2>&1 | head -14")
        print()
        print("  # 3. 125 bp")
        print("  " + RUN125.replace("\n", "\n  "))
        print()
        print("⚠ WATCH THE LAMBDA CONVENTION G3 PICKS. 121 measured lambda")
        print("  moving 5.37 -> 7.19 (+34%) between read lengths, so G3 may")
        print("  select a DIFFERENT convention at 125 bp. Not a bug — but if")
        print("  it differs from the 50 bp winner, the two tables' slopes were")
        print("  computed under different conventions and S4 must say so.")
        print()
        print("⚠ G1/G2/G3 are gates on values supplied here, so a PASS means")
        print("  the arms and the record agree. The independent confirmations")
        print("  are separate: 45 reproduced 33's slopes to 4 dp from its own")
        print("  deming(), and 136's split-control arms reproduced the 50 bp")
        print("  tables byte-identically before the 125 bp ones were believed.")
        print()
        print("⬜ GENOME_WIDE's Brigham and Fuwai are REPLACED, not updated —")
        print("  0.757 and 0.359 match nothing in any surviving output. Record")
        print("  the substitution, and rewrite 62's line-16 warning, which")
        print("  explains 0.757 as a genome-wide/selected distinction when it")
        print("  is simply unreproducible.")
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
    print("✓ verified on disk: all four arms and the output follow SUF; all")
    print("  thirteen record constants are overridable and required when")
    print("  suffixed; the --list print derives from them; the docstring no")
    print("  longer states a read length as if it were the only one.")
    return 0


def self_test():
    import subprocess
    import tempfile

    fx = ('"""\n' + A_DOC
          + "\nreproduces 0.584 and 0.988**.\n"
          + "which one reproduces the recorded 0.584 and 0.988 "
            "**in both arms at once**.\n"
          + "  G3  one lambda convention reproduces BOTH 0.584 and 0.988 "
            "within tolerance\n"
          + '"""\n'
            "import os\n"
          + A_ARMS
          + "EXPECT_COMMON = 203\n"
            "EXPECT_A = 0.584\n"
            "EXPECT_B = 0.988\n"
            "EXPECT_A_CI = (0.478, 0.693)\n"
            "EXPECT_B_CI = (0.889, 1.089)\n"
            "EXPECT_DIR_A = 162\n"
            "EXPECT_DIR_B = 174\n"
            "EXPECT_PRODUCT = 0.577\n"
            "EXPECT_OBSERVED_C = 0.521\n"
            'PER_COHORT = {"Colorado": 0.521, "Brigham": 0.476, '
            '"Stanford": 0.386, "Fuwai": 0.243}\n'
          + A_REC
          + "def main():\n"
            "    if True:\n"
          + A_LIST
          + '    print("OUT", OUT)\n'
            '    print("ARM", ARMS["selection"], ARMS["magnet_other"])\n'
            '    print("COMMON", EXPECT_COMMON, "A", EXPECT_A, "B", EXPECT_B)\n'
            '    print("ACI", EXPECT_A_CI, "BCI", EXPECT_B_CI)\n'
            '    print("DIR", EXPECT_DIR_A, EXPECT_DIR_B)\n'
            '    print("PROD", EXPECT_PRODUCT, "C", EXPECT_OBSERVED_C)\n'
            '    print("PC", PER_COHORT["Colorado"], PER_COHORT["Fuwai"])\n'
            '    print("GW", GENOME_WIDE["Brigham"], GENOME_WIDE["Fuwai"])\n'
            "    _require_arms()\n"
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
        for k in ("TBL_SUFFIX", "S4_RECORD", "S4_PER_COHORT",
                  "S4_GENOME_WIDE"):
            e.pop(k, None)
        e.update(env or {})
        return subprocess.run([sys.executable, "-I", p], capture_output=True,
                              text=True, cwd=d, env=e)

    # 1. unsuffixed -> every original value, every original string
    r = run()
    assert r.returncode == 0, r.stdout + r.stderr
    assert "OUT S4.tsv" in r.stdout, r.stdout
    assert ("ARM mp_SE_selection.all_events.tsv "
            "mp_SE_magnet_other.all_events.tsv") in r.stdout, r.stdout
    assert "COMMON 203 A 0.584 B 0.988" in r.stdout, r.stdout
    assert "ACI (0.478, 0.693) BCI (0.889, 1.089)" in r.stdout, r.stdout
    assert "DIR 162 174" in r.stdout, r.stdout
    assert "PROD 0.577 C 0.521" in r.stdout, r.stdout
    assert "PC 0.521 0.243" in r.stdout, r.stdout
    assert "GW 0.757 0.359" in r.stdout, r.stdout
    assert "The record: arm A 0.584 (0.478-0.693) platform only;" in r.stdout, \
        r.stdout

    # 2. TBL_SUFFIX alone must refuse and name all three groups + the reason
    r2 = run({"TBL_SUFFIX": "_rl125"})
    assert r2.returncode != 0, r2.stdout + r2.stderr
    b2 = r2.stdout + r2.stderr
    assert "REFUSING" in b2, b2
    for k in ("S4_RECORD", "S4_PER_COHORT", "S4_GENOME_WIDE"):
        assert k in b2, (k, b2)
    assert "refuses a CORRECT table" in b2, b2
    assert "common=161" in b2 and "C=0.600" in b2, b2

    full = {"TBL_SUFFIX": "_rl125",
            "S4_RECORD": "common=161,A=0.690,B=0.968,A_CI=0.570:0.806,"
                         "B_CI=0.875:1.069,dirA=130,dirB=140,"
                         "product=0.668,C=0.600",
            "S4_PER_COHORT": "Colorado=0.5999,Brigham=0.5099,"
                             "Stanford=0.4008,Fuwai=0.2948",
            "S4_GENOME_WIDE": "Colorado=0.674,Brigham=0.722,"
                              "Fuwai=0.365,Stanford=0.605"}

    # 3. full record, arms ABSENT -> _require_arms refuses, naming them
    r3 = run(full)
    assert r3.returncode != 0, r3.stdout + r3.stderr
    b3 = r3.stdout + r3.stderr
    assert "these arms are absent" in b3, b3
    assert "mp_SE_selection_rl125.all_events.tsv" in b3, b3
    # but every constant moved before the refusal
    assert "COMMON 161 A 0.69 B 0.968" in r3.stdout, r3.stdout
    assert "ACI (0.57, 0.806) BCI (0.875, 1.069)" in r3.stdout, r3.stdout
    assert "DIR 130 140" in r3.stdout, r3.stdout
    assert "PROD 0.668 C 0.6" in r3.stdout, r3.stdout
    assert "PC 0.5999 0.2948" in r3.stdout, r3.stdout
    assert "GW 0.722 0.365" in r3.stdout, r3.stdout
    assert "OUT S4_rl125.tsv" in r3.stdout, r3.stdout
    assert "ARM mp_SE_selection_rl125.all_events.tsv" in r3.stdout, r3.stdout
    assert "The record: arm A 0.69 (0.57-0.806) platform only;" in r3.stdout, \
        r3.stdout

    # 4. with the arms present it runs clean
    for v in ("mp_SE_selection_rl125.all_events.tsv",
              "mp_SE_magnet_matched_rl125.all_events.tsv",
              "mp_SE_lbnl_matched_rl125.all_events.tsv",
              "mp_SE_magnet_other_rl125.all_events.tsv"):
        open(os.path.join(d, v), "w").write("x\n")
    r4 = run(full)
    assert r4.returncode == 0, r4.stdout + r4.stderr

    # 5. a short / misspelt key set fails loudly and names the expected keys
    for var, bad_val, want in (
            ("S4_RECORD", "common=161,A=0.690", "A_CI"),
            ("S4_PER_COHORT", "Colorado=0.6,Brigham=0.5,Stanford=0.4",
             "Fuwai"),
            ("S4_GENOME_WIDE", "Colorado=0.674,Brigham=0.722,Fuwai=0.365,"
                               "Stanfrod=0.605", "Stanford")):
        rb = run(dict(full, **{var: bad_val}))
        assert rb.returncode != 0, (var, rb.stdout)
        bb = rb.stdout + rb.stderr
        assert "must name exactly" in bb, (var, bb)
        assert want in bb, (var, want, bb)

    # 6. the groups work WITHOUT a suffix too (a 50 bp re-record)
    r6 = run({"S4_RECORD": "common=203,A=0.584,B=0.988,A_CI=0.478:0.693,"
                           "B_CI=0.889:1.089,dirA=162,dirB=174,"
                           "product=0.577,C=0.521"})
    assert r6.returncode == 0, r6.stdout + r6.stderr
    assert "OUT S4.tsv" in r6.stdout, r6.stdout
    assert "COMMON 203 A 0.584" in r6.stdout, r6.stdout

    # 7. double application refused, text unchanged
    n7, p7 = apply_edits(new)
    assert p7, "second application was not refused"
    assert n7 == new, "a refused re-application still altered the text"

    print("[ok] self-test passed: unsuffixed output is unchanged across all "
          "thirteen constants, both arms strings and the --list record line; "
          "TBL_SUFFIX alone REFUSES, names all three groups, says it would "
          "refuse a CORRECT table and prints the full command; the full "
          "record moves every constant including EXPECT_COMMON 203->161 and "
          "the CIs; absent 125 bp arms refuse early and by name rather than "
          "at G1; a short or misspelt key set in any of the three groups "
          "names the expected keys; the groups work without a suffix; double "
          "application refused")
    return 0


if __name__ == "__main__":
    sys.exit(main())
