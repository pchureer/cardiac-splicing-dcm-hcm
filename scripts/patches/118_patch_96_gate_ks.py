#!/usr/bin/env python3
"""
118_patch_96_gate_ks.py -- put a GATE under 96's KS constants, which currently
answer to nothing.

WHAT IS WRONG
-------------
`96_lvef_column_dictionary.py` lines 78-81:

    # the KS deviations quoted in the FDR_lvef dictionary entry. The PROSE reads
    # these, so the deposited text can never disagree with the gate.
    KS_D_RAW  = float(os.environ.get("KS_D_RAW",  "0.050"))
    KS_D_RANK = float(os.environ.get("KS_D_RANK", "0.013"))

and lines 139-140 interpolate them into the deposited dictionary:

    f"uniform is D = {KS_D_RAW:.3f}, {KS_D_RAW / KS_D_RANK:.0f} times the rank "
    f"model's {KS_D_RANK:.3f}. These {EXPECT['FDR_lvef']} are "

The comment is true in the harmless direction and false in the dangerous one.
The prose cannot disagree with the CONSTANT -- but **there is no gate.** Nothing
in the file computes a KS statistic (`grep -n 'ks_2samp\\|kstest'` is empty).
`EXPECT` is gated at line 331 and `EXPECT_ROWS` at line 312; the KS values are
an ungated input that goes straight into a deposited file.

So running 96 at 125 bp with the defaults deposits "D = 0.050, 4 times the rank
model's 0.013" -- a 50 bp sentence beside 125 bp data. That is exactly the
defect 117 removed from 95, one layer down and in the same supplement.

WHY IT CAN BE GATED
-------------------
The table 96 reads carries both underlying p-value columns:

    ekey gene GeneID n psi_mean psi_sd beta_lvef se_lvef p_lvef
    beta_rank p_rank rho_lvef p_rho dPSI_per_10pp_lvef
    FDR_lvef FDR_rank FDR_rho

D is the one-sample Kolmogorov-Smirnov statistic of each p-value distribution
against Uniform(0,1) -- computable from the file, with no scipy:

    D = max_i  max( i/n - p_(i),  p_(i) - (i-1)/n )     over sorted p

WHAT THIS PATCH DOES
--------------------
Inserts **G1b** between the row gate (G1) and the counts gate (G2): recompute
both D values from `p_lvef` and `p_rank`, and refuse unless they agree with
`KS_D_RAW` / `KS_D_RANK` at the 3 decimals the prose prints.

The refusal prints the MEASURED values and the exact re-run line, so the first
run at a new read length tells you what to pass. The deposited sentence then
cannot be a transcription: either the constants match the data or nothing is
written.

It also corrects the comment, which claimed a gate that did not exist.

Why a gate rather than deriving in place (as 117 did for 95): `DICT` is a
module-level dict literal built at import time, before the data is read. Making
the prose derive would mean restructuring `DICT` into templates filled at
runtime -- a larger change to a file this patch cannot see in full. A gate that
refuses is sufficient: an unverified number can no longer reach the deposit.

CONVENTIONS
-----------
  * dry run by default; --apply writes after a timestamped .bak
  * anchors asserted to occur EXACTLY once
  * the patched file is re-read FROM DISK and re-verified; restore on failure
  * non-zero exit on any refusal
  * no new imports -- pure python, so it cannot fail on a missing dependency

USAGE
-----
  python 118_patch_96_gate_ks.py
  python 118_patch_96_gate_ks.py --apply
  python 118_patch_96_gate_ks.py --self-test
"""
import argparse
import difflib
import os
import shutil
import sys
import time

TARGET = "96_lvef_column_dictionary.py"

OLD_COMMENT = (
    "# the KS deviations quoted in the FDR_lvef dictionary entry. The PROSE reads\n"
    "# these, so the deposited text can never disagree with the gate.\n"
)
NEW_COMMENT = (
    "# The KS deviations quoted in the FDR_lvef dictionary entry.\n"
    "#\n"
    "# These are INPUTS to deposited prose. Until G1b was added below, nothing\n"
    "# in this file computed a KS statistic at all, so the comment that used to\n"
    "# stand here pointed at a gate that did not exist. The prose agreed with\n"
    "# the CONSTANT; the constant answered to nothing. Run at another read\n"
    "# length, it would have deposited a 50 bp sentence beside the new data.\n"
    "#\n"
    "# G1b now recomputes both values from p_lvef and p_rank and refuses unless\n"
    "# they match to the 3 decimals the prose prints. The comment is true now.\n"
)

ANCHOR_G2 = "    for col, want in EXPECT.items():\n"

G1B = '''    # ----------------------------------------------------------------- G1b
    # The KS deviations are quoted in the deposited dictionary. Measure them
    # from the data and refuse if the constants disagree.
    #
    # D is the one-sample Kolmogorov-Smirnov statistic of the p-value
    # distribution against Uniform(0,1). Computed directly so this file needs
    # no scipy:
    #     D = max_i max( i/n - p_(i), p_(i) - (i-1)/n )   over sorted p
    def _ks_uniform(vals):
        s = sorted(float(v) for v in vals
                   if v == v and 0.0 <= float(v) <= 1.0)
        m = len(s)
        if m == 0:
            return None
        d = 0.0
        for _i, _v in enumerate(s, 1):
            a, b = _i / m - _v, _v - (_i - 1) / m
            if a > d:
                d = a
            if b > d:
                d = b
        return d

    _ks_meas = {}
    for _name, _col in (("KS_D_RAW", "p_lvef"), ("KS_D_RANK", "p_rank")):
        if _col not in df.columns:
            refuse(f"G1b — no `{_col}` column, so the KS deviation quoted in "
                   f"the deposited\\n  dictionary cannot be verified. It would "
                   f"be a transcribed constant, which is\\n  what this gate "
                   f"exists to prevent.")
        _ks_meas[_name] = _ks_uniform(
            pd.to_numeric(df[_col], errors="coerce").dropna().tolist())
        if _ks_meas[_name] is None:
            refuse(f"G1b — `{_col}` holds no usable p-values in [0, 1].")

    print(f"\\n  G1b measured from the data:  "
          f"D(p_lvef) = {_ks_meas['KS_D_RAW']:.4f}   "
          f"D(p_rank) = {_ks_meas['KS_D_RANK']:.4f}")
    print(f"      constants in use:        "
          f"KS_D_RAW = {KS_D_RAW:.3f}   KS_D_RANK = {KS_D_RANK:.3f}")
    _ks_bad = [(_n, _ks_meas[_n], _c) for _n, _c in
               (("KS_D_RAW", KS_D_RAW), ("KS_D_RANK", KS_D_RANK))
               if round(_ks_meas[_n], 3) != round(_c, 3)]
    if _ks_bad:
        _lines = "\\n".join(
            f"    {_n}: measured {_m:.4f} -> {_m:.3f}, constant {_c:.3f}"
            for _n, _m, _c in _ks_bad)
        refuse(
            "G1b — the KS deviation(s) quoted in the deposited dictionary do "
            "not match\\n  this table:\\n" + _lines + "\\n\\n"
            "  The dictionary says 'deviation from uniform is D = <KS_D_RAW>, "
            "<ratio> times\\n  the rank model's <KS_D_RANK>'. Depositing that "
            "beside a table it was not\\n  measured from is a transcription.\\n\\n"
            "  Re-run with the measured values:\\n"
            f"    KS_D_RAW={_ks_meas['KS_D_RAW']:.3f} "
            f"KS_D_RANK={_ks_meas['KS_D_RANK']:.3f} \\\\\\n"
            "    <your other env vars> python3 96_lvef_column_dictionary.py "
            "--apply\\n\\n"
            "  Then CHECK THE RATIO: the prose prints "
            f"{_ks_meas['KS_D_RAW'] / _ks_meas['KS_D_RANK']:.0f} times, "
            "and a\\n  ratio that has moved is a change to a deposited claim, "
            "not a rounding detail.")
    print(f"  G1b PASS — the deposited KS sentence is measured from this table "
          f"({_ks_meas['KS_D_RAW'] / _ks_meas['KS_D_RANK']:.0f} times)")

'''

EDITS = [
    ("the comment that claimed a gate which did not exist",
     OLD_COMMENT, NEW_COMMENT, 1),
    ("G1b — measure the KS deviations and gate the constants",
     ANCHOR_G2, G1B + ANCHOR_G2, 1),
]


def apply_edits(text):
    problems, out = [], text
    for name, anchor, repl, want in EDITS:
        got = out.count(anchor)
        if got != want:
            problems.append(
                f"{name}: anchor occurs {got}x, expected {want}x\n"
                f"      anchor: {anchor.strip().splitlines()[0][:88]}")
            continue
        if anchor not in repl and repl in out:
            problems.append(f"{name}: replacement already present")
            continue
        # Check for the CODE marker, not the string "G1b": the corrected
        # comment inserted by the first edit mentions G1b, so testing for the
        # name would make this misfire on its own earlier edit.
        if "def _ks_uniform(" in out and name.startswith("G1b"):
            problems.append("the G1b estimator is already present — "
                            "already patched")
            continue
        out = out.replace(anchor, repl, want)
    return out, problems


def verify(text):
    bad = []
    code = "\n".join(l for l in text.splitlines()
                     if not l.lstrip().startswith("#"))
    if "def _ks_uniform(vals):" not in code:
        bad.append("the KS estimator was not inserted")
    if "_ks_bad" not in code:
        bad.append("the comparison against the constants is missing")
    if "the deposited text can never disagree with the gate" in text:
        bad.append("the false comment survives")
    # the gate must sit BEFORE the counts gate it precedes
    i_ks = text.find("def _ks_uniform(vals):")
    i_g2 = text.find("for col, want in EXPECT.items():")
    if i_ks < 0 or i_g2 < 0 or i_ks > i_g2:
        bad.append("G1b is not above the EXPECT counts gate")
    # and AFTER the row gate, so a wrong table refuses on rows first
    i_g1 = text.find("if len(df) != EXPECT_ROWS:")
    if i_g1 < 0:
        bad.append("the EXPECT_ROWS gate is missing from the file")
    elif i_g1 > i_ks:
        bad.append("G1b was inserted above the row gate; a wrong table would "
                   "report a KS mismatch instead of a row-count mismatch")
    if "round(_ks_meas[_n], 3) != round(_c, 3)" not in code:
        bad.append("the comparison is not at the 3 decimals the prose prints")
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
    print("118 — gate 96's KS constants against the data")
    print("=" * 78)
    if not os.path.exists(a.target):
        print(f"[fatal] {a.target} not found.", file=sys.stderr)
        return 2
    src = open(a.target, errors="replace").read()
    if "refuse(" not in src:
        print("[fatal] no refuse() in this file; the gate calls it.",
              file=sys.stderr)
        return 2

    new, problems = apply_edits(src)
    if problems:
        print("REFUSING — anchors do not match:")
        for p in problems:
            print(f"  ✗ {p}")
        print("\n  Nothing written. Send me the lines around 'KS_D_RAW' and")
        print("  'for col, want in EXPECT' and I will re-anchor.")
        return 3
    bad = verify(new)
    if bad:
        print("REFUSING — patched text fails its own verification:")
        for b in bad:
            print(f"  ✗ {b}")
        return 3

    print("".join(difflib.unified_diff(
        src.splitlines(True), new.splitlines(True),
        fromfile=a.target, tofile=a.target + " (patched)", n=2)))
    print("=" * 78)

    if not a.apply:
        print("DRY RUN — nothing written. Re-run with --apply.")
        print()
        print("Then run 96 at 125 bp. The KS constants are still the 50 bp")
        print("defaults, so G1b will REFUSE and print the measured values:")
        print()
        print("  DATA=SupplementaryData_lvef_all_events_rl125.tsv \\")
        print("  OUT=SupplementaryData_lvef_all_events_rl125.columns.md \\")
        print("  EXPECT_ROWS=33130 \\")
        print("  EXPECT=FDR_rank=1,FDR_lvef=125,FDR_rho=1 \\")
        print("  python3 96_lvef_column_dictionary.py --list")
        print()
        print("Take the two numbers from that refusal, add them to the prefix,")
        print("and run again with --apply. CHECK THE RATIO before depositing:")
        print("the prose says 'N times', and N moving is a changed claim.")
        print()
        print("⚠ First, run the SAME patched script against the 50 bp table")
        print("  with the defaults. It must PASS — that is the calibration:")
        print()
        print("  DATA=SupplementaryData_lvef_all_events.tsv \\")
        print("  OUT=/dev/null EXPECT_ROWS=38084 \\")
        print("  python3 96_lvef_column_dictionary.py --list")
        print()
        print("  If 0.050 / 0.013 do NOT reproduce from the 50 bp table, the")
        print("  published sentence was never measured from it either, and that")
        print("  is a finding about the 50 bp record, not a bug in this gate.")
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
        print(f"  restored from {bak}")
        return 3
    print("✓ verified on disk: G1b sits between the row gate and the counts")
    print("  gate, compares at 3 decimals, and the false comment is gone.")
    return 0


def self_test():
    # the estimator itself, against values computable by hand
    src_fix = (
        "import os\n"
        "import pandas as pd\n"
        "def refuse(m):\n    raise SystemExit('REFUSED: ' + m)\n"
        "KS_D_RAW = float(os.environ.get('KS_D_RAW', '0.050'))\n"
        "KS_D_RANK = float(os.environ.get('KS_D_RANK', '0.013'))\n"
        "EXPECT_ROWS = int(os.environ.get('EXPECT_ROWS', '4'))\n"
        "EXPECT = {'FDR_rank': 0}\n"
        "# the KS deviations quoted in the FDR_lvef dictionary entry. "
        "The PROSE reads\n"
        "# these, so the deposited text can never disagree with the gate.\n"
        "def run(df):\n"
        "    if len(df) != EXPECT_ROWS:\n"
        "        refuse('G1 rows')\n"
        "    for col, want in EXPECT.items():\n"
        "        pass\n"
        "    return 'ok'\n")
    new, probs = apply_edits(src_fix)
    assert not probs, probs
    assert not verify(new), verify(new)
    assert "never disagree with the gate" not in new

    import subprocess
    import tempfile
    d = tempfile.mkdtemp()
    p = os.path.join(d, "fx.py")
    # p_lvef uniform -> D small ; p_rank heavily skewed -> D large.
    runner = (
        "\nif __name__ == '__main__':\n"
        "    import pandas as pd\n"
        "    n = 400\n"
        "    pl = [(i + 0.5) / n for i in range(n)]\n"
        "    pr = [((i + 0.5) / n) ** 3 for i in range(n)]\n"
        "    df = pd.DataFrame({'p_lvef': pl, 'p_rank': pr})\n"
        "    print(run(df))\n")
    open(p, "w").write(new.replace("EXPECT_ROWS', '4'", "EXPECT_ROWS', '400'")
                       + runner)

    # with the 50 bp defaults the fixture must REFUSE (its data is not that)
    r = subprocess.run([sys.executable, "-I", p], capture_output=True,
                       text=True, cwd=d)
    out = r.stdout + r.stderr
    assert r.returncode != 0, out
    assert "G1b" in out and "Re-run with the measured values" in out, out

    # the refusal must print values that then PASS
    import re as _re
    m = _re.search(r"KS_D_RAW=([0-9.]+) KS_D_RANK=([0-9.]+)", out)
    assert m, out
    env = dict(os.environ, KS_D_RAW=m.group(1), KS_D_RANK=m.group(2))
    r2 = subprocess.run([sys.executable, "-I", p], capture_output=True,
                        text=True, cwd=d, env=env)
    assert r2.returncode == 0, r2.stdout + r2.stderr
    assert "G1b PASS" in r2.stdout, r2.stdout
    # a uniform column must give a small D, a cubed one a large D
    m2 = _re.search(r"D\(p_lvef\) = ([0-9.]+)\s+D\(p_rank\) = ([0-9.]+)",
                    r2.stdout)
    assert m2, r2.stdout
    assert float(m2.group(1)) < 0.01, m2.group(1)
    assert float(m2.group(2)) > 0.20, m2.group(2)

    # a missing p column must refuse, not silently skip
    p3 = os.path.join(d, "fx3.py")
    open(p3, "w").write(new.replace("EXPECT_ROWS', '4'", "EXPECT_ROWS', '2'")
                        + "\nif __name__ == '__main__':\n"
                          "    import pandas as pd\n"
                          "    print(run(pd.DataFrame({'p_lvef': [.1, .2]})))\n")
    r3 = subprocess.run([sys.executable, "-I", p3], capture_output=True,
                        text=True, cwd=d)
    assert r3.returncode != 0 and "p_rank" in (r3.stdout + r3.stderr), \
        r3.stdout + r3.stderr

    # double application refused
    _n4, p4 = apply_edits(new)
    assert p4, "second application was not refused"

    print("[ok] self-test passed: estimator separates a uniform column (D<0.01) "
          "from a skewed one (D>0.20); the gate refuses the 50 bp defaults on "
          "other data, prints measured values that then PASS, refuses a missing "
          "p column, and double application is refused")
    return 0


if __name__ == "__main__":
    sys.exit(main())
