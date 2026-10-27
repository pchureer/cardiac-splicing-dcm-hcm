#!/usr/bin/env python3
"""
117_patch_95_null_derive.py -- make 95's analytic null DERIVE from the data
instead of being a transcribed string literal.

WHAT IS WRONG
-------------
`95_s1c_column_and_s3_split.py` builds Supplementary Table 10 and includes:

    ("analytic null expectation for |rho|", "0.0538"),

A QUOTED STRING. Nothing computes it, nothing checks it, and nothing reads it,
so it survived the entire 125 bp re-baselining untouched. It would print 0.0538
at any read length, at any n, forever.

Meanwhile `77_lvef_s3.py` line 687 does compute it:

    exp_med_rho = 0.6745 / np.sqrt(max(n_med - 1, 2))   # half-normal median

which at n_med = 157 gives 0.054002 -> 0.0540. The two scripts therefore print
0.0540 and 0.0538 for one quantity. 0.0538 is the same formula with df = n
instead of df = n-1. Under the null Var(rho) = 1/(n-1), so **77 is right**.

WHAT THIS PATCH DOES
--------------------
Replaces the literal with a derivation that reproduces 77's expression exactly,
adds a gate (G6) that refuses when the derived value disagrees with a recorded
expectation, and splits the deposited row in two so the table states the n the
null was evaluated at and the range of n across events -- because n is NOT
constant (144..160), and a bare scalar null stands for a distribution.

Deposited row, before:

    analytic null expectation for |rho|        0.0538

after:

    analytic null for the median |rho|         0.0540
      evaluated at n (median across events)    157 (range 144-160)

The label changes too: "expectation for |rho|" invites the null's MEAN
(0.0637 here). The quantity is the null's MEDIAN. The imprecise label is what
let a mean/median mix-up go unnoticed for this long.

CONVENTIONS THIS FOLLOWS
------------------------
  * dry run by default; --apply writes, after a timestamped .bak
  * every anchor asserted to occur EXACTLY the expected number of times
  * the patched file is re-read FROM DISK and the edits re-verified
  * NON-ZERO EXIT ON ANY REFUSAL  (lesson from 110: an exit-0 refusal led to
    job 407846 being submitted against an unpatched script)
  * no numpy/math import added -- uses ** 0.5, so the patch cannot fail on a
    missing import in a file we cannot see

G6's expectation comes from the environment, so an unparameterised run still
derives and prints the value but does not gate:

    S3_EXPECT_NULL=0.0540 python3 95_s1c_column_and_s3_split.py --apply ...

USAGE
-----
  python 117_patch_95_null_derive.py                      # dry run, shows diff
  python 117_patch_95_null_derive.py --apply
  python 117_patch_95_null_derive.py --self-test
"""
import argparse
import difflib
import os
import shutil
import sys
import time

TARGET = "95_s1c_column_and_s3_split.py"

# ---------------------------------------------------------------- the edits
# (name, anchor, replacement, expected occurrences)

OLD_ROW = '        ("analytic null expectation for |rho|", "0.0538"),\n'
NEW_ROW = (
    '        ("analytic null for the median |rho|", f"{_null_med:.4f}"),\n'
    '        ("  evaluated at n (median across events)",\n'
    '         f"{_n_med:g} (range {_n_min:g}-{_n_max:g})"),\n'
)

ANCHOR_SUMMARY = "    summary = pd.DataFrame([\n"

DERIVE_BLOCK = '''    # ------------------------------------------------------------------ G6
    # The null comparator must DERIVE from the data. It used to be the string
    # literal "0.0538", which nothing computed and nothing checked, so it
    # survived the 125 bp re-baselining unchanged.
    #
    # The statistic compared is a MEDIAN |rho|, so the comparator is the
    # null's MEDIAN, not its mean:
    #     rho ~ N(0, 1/(n-1))  under the null
    #     median |rho| = z(0.75) * sigma = 0.6745 / sqrt(n-1)
    #     mean   |rho| = sqrt(2/pi) * sigma = 0.7979 / sqrt(n-1)   <- NOT this
    # At n = 157 those are 0.0540 and 0.0637. Quoting the mean against a
    # measured median would overstate the null by 18%.
    #
    # Var(rho) = 1/(n-1), hence n-1 and not n. df = n gives 0.0538, which is
    # where this file's old literal came from.
    #
    # This expression is character-for-character the one in 77_lvef_s3.py
    # line 687. If one changes, both must.
    _ncol = next((c for c in ("n", "n_used", "n_obs", "nobs", "N")
                  if c in df.columns), None)
    if _ncol is None:
        refuse("no sample-size column in the S3 table, so the analytic null "
               "cannot be derived.\\n  A transcribed constant is exactly what "
               "this gate exists to prevent; refusing to\\n  write a summary "
               "carrying a number the data cannot support.")
    _nn = pd.to_numeric(df[_ncol], errors="coerce").dropna()
    _n_med = float(_nn.median())
    _n_min, _n_max = float(_nn.min()), float(_nn.max())
    _null_med = 0.6745 / (max(_n_med - 1.0, 2.0) ** 0.5)
    print(f"\\n  G6  analytic null for the median |rho|: {_null_med:.4f}"
          f"  (df = n-1, n = {_n_med:g}, range {_n_min:g}-{_n_max:g})")
    if _n_min != _n_max:
        print(f"      n VARIES across events, so this scalar stands for a "
              f"distribution;\\n      the deposited table names the n it was "
              f"evaluated at.")
    _exp_null = float(os.environ.get("S3_EXPECT_NULL", "0") or 0)
    if _exp_null > 0:
        # Equality at the DEPOSITED precision (4 dp), not a tolerance. 0.0538
        # and 0.0540 differ by 1.7e-4, so any tolerance at 5e-4 or looser would
        # pass the very disagreement this gate exists to catch.
        if round(_null_med, 4) != round(_exp_null, 4):
            refuse(f"G6 — derived null {_null_med:.4f} differs from the "
                   f"recorded {_exp_null:.4f}.\\n  77_lvef_s3.py and this "
                   f"script compute THE SAME QUANTITY and must agree. Two "
                   f"values\\n  for one quantity is how 0.0538 and 0.0540 both "
                   f"ended up in the record.")
        print(f"  G6 PASS — derived (not transcribed) and equal to the "
              f"recorded {_exp_null:.4f}")
    else:
        print("  G6 — derived, not transcribed. Set S3_EXPECT_NULL to gate it "
              "against 77's value.")

'''

EDITS = [
    ("the transcribed literal -> a derived value, label corrected",
     OLD_ROW, NEW_ROW, 1),
    ("the derivation + G6, inserted before the summary is built",
     ANCHOR_SUMMARY, DERIVE_BLOCK + ANCHOR_SUMMARY, 1),
]


def apply_edits(text):
    """-> (new_text, [problems]). Does not touch disk."""
    problems = []
    out = text
    for name, anchor, repl, want in EDITS:
        got = out.count(anchor)
        if got != want:
            problems.append(
                f"{name}: anchor occurs {got}x, expected {want}x\n"
                f"      anchor: {anchor.strip()[:88]}")
            continue
        if repl in out and anchor not in repl:
            problems.append(f"{name}: replacement already present")
            continue
        out = out.replace(anchor, repl, want)
    return out, problems


def verify(text):
    """Checks that must hold in the FINAL text. -> [problems]

    The quoted-literal check looks at CODE lines only: the patch's own
    explanatory comment quotes the old literal to say why it was removed,
    and a comment cannot put a number into a deposited table.
    """
    bad = []
    code = "\n".join(l for l in text.splitlines()
                     if not l.lstrip().startswith("#"))
    if '"analytic null expectation for |rho|", "0.0538"' in text:
        bad.append("the transcribed literal 0.0538 is still present")
    if '"0.0538"' in code or "'0.0538'" in code:
        bad.append("a quoted 0.0538 survives on a code line")
    if "_null_med = 0.6745 / (max(_n_med - 1.0, 2.0) ** 0.5)" not in text:
        bad.append("the derivation line is missing")
    if "S3_EXPECT_NULL" not in text:
        bad.append("the G6 expectation constant is missing")
    if text.count("_null_med") < 3:
        bad.append("_null_med is defined but not used in the summary row")
    # the derivation must come BEFORE the row that uses it
    i_def = text.find("_null_med = 0.6745")
    i_use = text.find('("analytic null for the median |rho|"')
    if i_def < 0 or i_use < 0 or i_def > i_use:
        bad.append("the derivation is not above the row that consumes it")
    # 0.6745 / sqrt(n) -- the off-by-one -- must not appear
    if "0.6745 / (max(_n_med, 2.0)" in text:
        bad.append("df = n (the off-by-one) was written instead of df = n-1")
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
    print(f"117 — derive 95's analytic null instead of transcribing it")
    print("=" * 78)
    if not os.path.exists(a.target):
        print(f"[fatal] {a.target} not found. Run from the directory holding it.",
              file=sys.stderr)
        return 2

    src = open(a.target, errors="replace").read()
    if "def refuse" not in src and "refuse(" not in src:
        print("[fatal] no refuse() in this file; the inserted gate calls it.",
              file=sys.stderr)
        return 2

    new, problems = apply_edits(src)
    if problems:
        print("REFUSING — anchors do not match:")
        for p in problems:
            print(f"  ✗ {p}")
        print()
        print("  Nothing was written. The file on the cluster differs from what")
        print("  this patch expects. Send me the lines around")
        print("  'analytic null expectation' and 'summary = pd.DataFrame' and")
        print("  I will re-anchor it.")
        return 3

    bad = verify(new)
    if bad:
        print("REFUSING — the patched text fails its own verification:")
        for b in bad:
            print(f"  ✗ {b}")
        return 3

    diff = list(difflib.unified_diff(
        src.splitlines(True), new.splitlines(True),
        fromfile=a.target, tofile=a.target + " (patched)", n=2))
    print("".join(diff))
    print("=" * 78)

    if not a.apply:
        print("DRY RUN — nothing written. Re-run with --apply.")
        print()
        print("Then, to gate it against 77's value:")
        print("  S3_EXPECT_NULL=0.0540 \\")
        print("  S3_EXPECT_FDR=1 S3_EXPECT_THRESH=0 \\")
        print("  EXPECT_ROWS_S3=33130 EXPECT_MEDIAN_ABS_RHO=0.0555 \\")
        print("  TBL_SUFFIX=_rl125 python3 95_s1c_column_and_s3_split.py --apply")
        print()
        print("  ⚠ the other env vars above are from the run that already")
        print("    passed; confirm them against your own command before reusing.")
        return 0

    bak = f"{a.target}.bak.{time.strftime('%Y%m%d_%H%M%S')}"
    shutil.copy2(a.target, bak)
    open(a.target, "w").write(new)
    print(f"[write] {a.target}   (backup: {bak})")

    # re-read FROM DISK and verify again
    disk = open(a.target, errors="replace").read()
    bad = verify(disk)
    if bad:
        print("✗ ON-DISK VERIFICATION FAILED — restoring the backup:")
        for b in bad:
            print(f"    {b}")
        shutil.copy2(bak, a.target)
        print(f"  restored from {bak}")
        return 3
    if disk != new:
        print("✗ the file on disk differs from what was written; restoring.")
        shutil.copy2(bak, a.target)
        return 3
    print("✓ verified on disk: the literal is gone, the derivation is above the")
    print("  row that uses it, df = n-1, and G6 is wired to S3_EXPECT_NULL.")
    return 0


def self_test():
    import tempfile
    fixture = (
        "import os\n"
        "import pandas as pd\n"
        "def refuse(m):\n    raise SystemExit(m)\n"
        "def build(df, rho, med, fdr, q, pcol, g, top):\n"
        "    summary = pd.DataFrame([\n"
        '        ("events tested", f"{len(df):,}"),\n'
        '        ("analytic null expectation for |rho|", "0.0538"),\n'
        '        ("event set", "SE_v2filtered"),\n'
        '    ], columns=["quantity", "value"])\n'
        "    return summary\n")
    new, probs = apply_edits(fixture)
    assert not probs, probs
    assert not verify(new), verify(new)
    _code = "\n".join(l for l in new.splitlines()
                      if not l.lstrip().startswith("#"))
    assert '"0.0538"' not in _code, "a quoted literal survives on a code line"
    assert new.index("_null_med = 0.6745") < new.index(
        '("analytic null for the median |rho|"')

    # the patched fixture must actually RUN and produce 0.0540 at n=157
    d = tempfile.mkdtemp()
    p = os.path.join(d, "fx.py")
    open(p, "w").write(new + (
        "\nif __name__ == '__main__':\n"
        "    import pandas as pd\n"
        "    df = pd.DataFrame({'n': [144, 157, 157, 160],\n"
        "                       'rho_lvef': [.01, .02, .03, .04]})\n"
        "    s = build(df, None, None, None, None, None, None, None)\n"
        "    v = dict(zip(s.quantity, s.value))\n"
        "    assert v['analytic null for the median |rho|'] == '0.0540', v\n"
        "    assert v['  evaluated at n (median across events)'] == "
        "'157 (range 144-160)', v\n"
        "    print('[fixture] 0.0540 at n=157, range reported')\n"))
    import subprocess
    r = subprocess.run([sys.executable, "-I", p], capture_output=True, text=True,
                       cwd=d)
    assert r.returncode == 0, r.stdout + r.stderr
    assert "0.0540 at n=157" in r.stdout, r.stdout

    # and G6 must REFUSE when the recorded value disagrees
    env = dict(os.environ, S3_EXPECT_NULL="0.0538")
    r2 = subprocess.run([sys.executable, "-I", p], capture_output=True,
                        text=True, cwd=d, env=env)
    assert r2.returncode != 0, "G6 did not refuse on a mismatched expectation"
    assert "G6" in (r2.stdout + r2.stderr), r2.stdout + r2.stderr

    # a file missing the anchor must be refused, not half-patched
    _n2, p2 = apply_edits("x = 1\n")
    assert p2 and len(p2) == 2, p2

    # applying twice must refuse rather than duplicate
    _n3, p3 = apply_edits(new)
    assert p3, "second application was not refused"

    os.remove(p)
    print("[ok] self-test passed: literal replaced, derivation placed above its "
          "consumer, patched fixture runs and yields 0.0540 at n=157 with the n "
          "range reported, G6 refuses a mismatched expectation, missing anchor "
          "and double application both refused")
    return 0


if __name__ == "__main__":
    sys.exit(main())
