#!/usr/bin/env python3
"""
131_patch_33_output_names.py -- 33 writes two HARD-CODED files. It already
overwrote both.

WHAT HAPPENED
-------------
`33_rebuild_summaries_dedup.py` takes `--discovery` and `--out`, and has no
other read-length parameter. The 125 bp run was therefore

    python3 33_rebuild_summaries_dedup.py \\
        --discovery "{ET}_v2filtered_rl125" --out summaries_dedup_rl125.txt

and `--out` protected the summary — `summaries_dedup.txt` is still the 28 Sep
50 bp file. But sections C and D write their names as string literals:

    line 172   I.to_csv("colorado_internal_control_dedup.tsv", ...)
    line 284   pd.DataFrame(rows).to_csv("direction_concordance_dedup.tsv", ...)

so BOTH were overwritten in place with 125 bp content under 50 bp names, and
the run printed `->  direction_concordance_dedup.tsv` and exited 0.

`direction_concordance_dedup.tsv` is read by 38 (Figure 4 panel A), 44, 84 and
85 (gate G3). `colorado_internal_control_dedup.tsv` has no reader at all.

A missing INPUT suffix fails loudly — file not found. A missing OUTPUT suffix
destroys the comparator and reports success. That asymmetry is why 33 passed a
byte-identical check on sections A-C while quietly replacing a file in D.

WHAT THIS PATCH DOES
--------------------
It does NOT add a fourth suffix convention. It makes `--discovery`, `--out` and
the two hard-coded names **unable to disagree**, because their disagreement is
the whole defect:

  * `SUF` comes from `TBL_SUFFIX`, as in 36/71/76/69/45 (Rule 9).
  * if `TBL_SUFFIX` is set, `--discovery` must END with it, else REFUSE.
  * if `TBL_SUFFIX` is unset but `--discovery` ends with a read-length suffix
    (`_rl` + digits), REFUSE — that is precisely the 11:53 run, and it is the
    one combination that silently destroys data.
  * `--out` defaults to None and is DERIVED as `summaries_dedup{SUF}.txt`; if
    given explicitly it must contain `SUF`, else REFUSE.
  * both hard-coded names, and the two `->` messages that echo them, take SUF.

`--discovery` is only required to *end* with the suffix, not to start with a
known stem, so a sensitivity run on `{ET}_v2compadj` or `{ET}_v2compadj_rl125`
still works.

THE REPAIR, AFTER APPLYING THIS
-------------------------------
  0. cp -n direction_concordance_dedup.tsv direction_concordance_dedup_rl125.tsv
     cp -n colorado_internal_control_dedup.tsv \\
           colorado_internal_control_dedup_rl125.tsv
  1. re-run 33 UNSUFFIXED  -> regenerates both 50 bp files from the 50 bp
     tables. Do NOT rebuild them from 85's REC_SE literals: 85's gate G3
     checks the file against REC_SE, so a file built from REC_SE makes that
     gate check the record against itself.
  2. re-run 33 SUFFIXED    -> regenerates both _rl125 files, which must come
     out byte-identical to the copies from step 0.

Step 1 must reproduce, from nothing this pipeline supplied (one-sided exact
binomial at p = 0.5, verified by reproducing all five 125 bp p-values exactly):

    Colorado 375/471 79.618  2.7177e-40      Brigham  359/442 81.222  2.8730e-42
    Stanford 408/536 76.119  2.3790e-35      Fuwai    494/714 69.188  2.1577e-25
    LBNL     151/199 75.879  6.6605e-14

USAGE
-----
  python 131_patch_33_output_names.py
  python 131_patch_33_output_names.py --apply
  python 131_patch_33_output_names.py --self-test
"""
import argparse
import difflib
import os
import shutil
import sys
import time

TARGET = "33_rebuild_summaries_dedup.py"

# ---------------------------------------------------------------- the insert
A_MAIN = "def main():\n"
B_MAIN = '''# --------------------------------------------------------------------------
# READ LENGTH.  33 writes two files whose names were string literals, so its
# 125 bp run overwrote the 50 bp `direction_concordance_dedup.tsv` (read by 38,
# 44, 84 and 85) and `colorado_internal_control_dedup.tsv` in place, and exited
# 0.  `--out` was the only output that carried a suffix, because it was the
# only one that was a parameter.
#
# SUF is read here; _resolve_suffix() below refuses every combination in which
# --discovery, --out and SUF could disagree.
SUF = os.environ.get("TBL_SUFFIX", "")


def _looks_readlength(stem):
    """True if `stem` ends in a read-length suffix such as `_rl125`."""
    i = stem.rfind("_rl")
    return i >= 0 and stem[i + 3:].isdigit() and len(stem[i + 3:]) > 0


def _resolve_suffix(a):
    """Refuse any disagreement between --discovery, --out and TBL_SUFFIX.

    Returns `a` with a.out resolved.  Called immediately after parse_args, so
    nothing is opened before the check.
    """
    d = str(a.discovery)
    if SUF:
        if not d.endswith(SUF):
            raise SystemExit(
                "REFUSING — TBL_SUFFIX is %r but --discovery is %r, which does\\n"
                "  not end with it. One of the two is wrong, and 33 writes two\\n"
                "  files whose names follow TBL_SUFFIX. Set them consistently:\\n"
                "    TBL_SUFFIX=%s python3 33_rebuild_summaries_dedup.py \\\\\\n"
                "        --discovery '{ET}_v2filtered%s'" % (SUF, d, SUF, SUF))
    else:
        if _looks_readlength(d):
            raise SystemExit(
                "REFUSING — --discovery is %r, which carries a read-length\\n"
                "  suffix, but TBL_SUFFIX is unset. That is exactly the run that\\n"
                "  overwrote the 50 bp direction_concordance_dedup.tsv and\\n"
                "  colorado_internal_control_dedup.tsv with 125 bp content under\\n"
                "  50 bp names, and reported success.\\n\\n"
                "  Re-run as:\\n"
                "    TBL_SUFFIX=%s python3 33_rebuild_summaries_dedup.py \\\\\\n"
                "        --discovery '%s'" % (d, d[d.rfind("_rl"):], d))
    if a.out is None:
        a.out = "summaries_dedup%s.txt" % SUF
    elif SUF and SUF not in str(a.out):
        raise SystemExit(
            "REFUSING — --out %r does not contain the suffix %r, so the summary\\n"
            "  and the two .tsv outputs would describe different read lengths.\\n"
            "  Use --out summaries_dedup%s.txt, or omit --out and it is derived."
            % (a.out, SUF, SUF))
    return a
# --------------------------------------------------------------------------


''' + A_MAIN

EDITS = [
    ("SUF and the refusal gate", A_MAIN, B_MAIN, 1),
    # --out becomes derivable; None means "derive it from SUF"
    ("--out default becomes derived",
     '    ap.add_argument("--out", default="summaries_dedup.txt")',
     '    ap.add_argument("--out", default=None,\n'
     '                    help="default: summaries_dedup<TBL_SUFFIX>.txt")', 1),
    ("the suffix check, before anything is opened",
     '    out = open(a.out, "w")',
     '    a = _resolve_suffix(a)\n'
     '    out = open(a.out, "w")', 1),
    # section C -- the quoted literal only; line 181's string is `  -> name`,
    # preceded by a space rather than a quote, so these two anchors are each
    # unique and cannot match the other line.
    ("section C output name",
     '"colorado_internal_control_dedup.tsv"',
     '"colorado_internal_control_dedup" + SUF + ".tsv"', 1),
    ("section C echo",
     'P("\\n  -> colorado_internal_control_dedup.tsv")',
     'P(f"\\n  -> colorado_internal_control_dedup{SUF}.tsv")', 1),
    ("section D output name",
     '"direction_concordance_dedup.tsv"',
     '"direction_concordance_dedup" + SUF + ".tsv"', 1),
    ("section D echo",
     'P("  -> direction_concordance_dedup.tsv")',
     'P(f"  -> direction_concordance_dedup{SUF}.tsv")', 1),
]


def apply_edits(text):
    if 'SUF = os.environ.get("TBL_SUFFIX", "")' in text:
        return text, ["SUF is already defined — 33 is already patched "
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
              "def _resolve_suffix(a):",
              "def _looks_readlength(stem):",
              "a = _resolve_suffix(a)",
              '"colorado_internal_control_dedup" + SUF + ".tsv"',
              '"direction_concordance_dedup" + SUF + ".tsv"',
              'P(f"\\n  -> colorado_internal_control_dedup{SUF}.tsv")',
              'P(f"  -> direction_concordance_dedup{SUF}.tsv")',
              'ap.add_argument("--out", default=None,'):
        if w not in text:
            bad.append(f"missing: {w[:60]}")
    for gone in ('"colorado_internal_control_dedup.tsv"',
                 '"direction_concordance_dedup.tsv"',
                 'default="summaries_dedup.txt"'):
        if gone in text:
            bad.append(f"unsuffixed site survives: {gone[:52]}")
    # the gate must run BEFORE the output file is opened
    i_g = text.find("a = _resolve_suffix(a)")
    i_o = text.find('out = open(a.out, "w")')
    if i_g < 0 or i_o < 0 or i_g > i_o:
        bad.append("the suffix gate does not precede open(a.out)")
    # SUF must be MODULE-LEVEL — i.e. inserted above `def main()`, not inside
    # it — because the two hard-coded names may be written from other
    # functions. This is the only ordering constraint that matters: a
    # module-level global is visible to every function whatever the textual
    # order, since it is resolved at call time, so NOT requiring SUF to appear
    # above each individual use (an earlier version of this check did, and
    # would have refused a correct patch on any layout with helpers above
    # main). Runtime behaviour is proven by the self-test, which executes it.
    i_s = text.find('SUF = os.environ.get("TBL_SUFFIX", "")')
    i_m = text.find("def main():")
    if i_s < 0 or i_m < 0 or i_s > i_m:
        bad.append("SUF is not module-level / not above def main()")
    if "\n    SUF = os.environ.get" in text:
        bad.append("SUF is indented — it would be a local, not a global")
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
    print("131 — 33's two hard-coded outputs. Both were already overwritten.")
    print("=" * 78)
    if not os.path.exists(a.target):
        print(f"[fatal] {a.target} not found.", file=sys.stderr)
        return 2
    src = open(a.target, errors="replace").read()
    if "import os" not in src:
        print("[fatal] no `import os` in 33 — add it first.", file=sys.stderr)
        return 2
    new, problems = apply_edits(src)
    if problems:
        print("REFUSING — anchors do not match:")
        for p in problems:
            print(f"  ✗ {p}")
        print("\n  Nothing written. Send me lines 100-110, 168-184, 280-290.")
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
        print("  # 0. preserve the 125 bp data under a name step 2 will NOT")
        print("  #    overwrite. `_rl125.tsv` is where step 2 writes, so the")
        print("  #    comparator has to live somewhere else.")
        print("  cp -n colorado_internal_control_dedup.tsv \\")
        print("        colorado_internal_control_dedup_rl125.tsv")
        print("  cp direction_concordance_dedup_rl125.tsv \\")
        print("     direction_concordance_dedup_rl125.PRESERVED")
        print("  cp colorado_internal_control_dedup_rl125.tsv \\")
        print("     colorado_internal_control_dedup_rl125.PRESERVED")
        print()
        print("  # the gate must refuse the run that caused this")
        print("  python3 33_rebuild_summaries_dedup.py \\")
        print("      --discovery '{ET}_v2filtered_rl125' 2>&1 | head -12")
        print("  #  -> expect REFUSING, TBL_SUFFIX unset")
        print()
        print("  # 1. regenerate the 50 bp files FROM THE 50 bp TABLES.")
        print("  #    --out must NOT be summaries_dedup.txt: that 28 Sep file is")
        print("  #    the only surviving record of the 50 bp Colorado internal")
        print("  #    control, so it is a comparator, not a target.")
        print("  python3 33_rebuild_summaries_dedup.py \\")
        print("      --out summaries_dedup_RESTORE.txt 2>&1 | tail -30")
        print("  cat direction_concordance_dedup.tsv")
        print("  diff summaries_dedup_RESTORE.txt summaries_dedup.txt")
        print()
        print("  # 2. regenerate the 125 bp files and compare against step 0")
        print("  TBL_SUFFIX=_rl125 python3 33_rebuild_summaries_dedup.py \\")
        print("      --discovery '{ET}_v2filtered_rl125' >/dev/null")
        print("  for f in direction_concordance_dedup \\")
        print("           colorado_internal_control_dedup; do")
        print("    cmp \"${f}_rl125.tsv\" \"${f}_rl125.PRESERVED\" \\")
        print("      && echo \"$f  byte-identical\" || echo \"$f  DIFFERS\"")
        print("  done")
        print()
        print("⚠ STEP 1 MUST REPRODUCE (one-sided exact binomial, p=0.5 —")
        print("  independently recomputed, nothing in the pipeline supplied it):")
        print("    Colorado 375/471 79.618  2.7177e-40")
        print("    LBNL     151/199 75.879  6.6605e-14")
        print("    Brigham  359/442 81.222  2.8730e-42")
        print("    Fuwai    494/714 69.188  2.1577e-25")
        print("    Stanford 408/536 76.119  2.3790e-35")
        print()
        print("  Do NOT rebuild these from 85's REC_SE literals — 85's gate G3")
        print("  checks the file against REC_SE, so that would make the gate")
        print("  check the record against itself.")
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
    print("✓ verified on disk: both hard-coded names and both echoes follow")
    print("  SUF; --out is derived; the gate runs before any file is opened.")
    return 0


def self_test():
    import subprocess
    import tempfile

    fx = (
        "import argparse, os\n"
        "def P(s):\n"
        "    print(s)\n"
        "def section_c(I):\n"
        '    I.to_csv("colorado_internal_control_dedup.tsv", sep="\\t", index=False)\n'
        '    P("\\n  -> colorado_internal_control_dedup.tsv")\n'
        "def section_d(rows):\n"
        '    pd.DataFrame(rows).to_csv("direction_concordance_dedup.tsv", sep="\\t",\n'
        "                             index=False)\n"
        '    P("  -> direction_concordance_dedup.tsv")\n'
        + A_MAIN
        + "    ap = argparse.ArgumentParser()\n"
          '    ap.add_argument("--fdr", type=float, default=0.05)\n'
          '    ap.add_argument("--discovery", default="{ET}_v2filtered")\n'
          '    ap.add_argument("--out", default="summaries_dedup.txt")\n'
          "    a = ap.parse_args()\n"
          '    out = open(a.out, "w")\n'
          "    out.close()\n"
          '    print("OUT", a.out)\n'
          '    print("DISC", a.discovery)\n'
          "    class _I:\n"
          "        def to_csv(self, p, **k):\n"
          '            print("WRITE", p)\n'
          "    section_c(_I())\n"
          '    print("DNAME", "direction_concordance_dedup" + SUF + ".tsv")\n'
          "if __name__ == '__main__':\n"
          "    main()\n")

    new, probs = apply_edits(fx)
    assert not probs, probs
    assert not verify(new), verify(new)
    compile(new, "fx", "exec")          # the inserted block must parse

    d = tempfile.mkdtemp()
    p = os.path.join(d, "fx.py")
    open(p, "w").write(new)

    def run(args=(), env=None):
        e = dict(os.environ)
        e.pop("TBL_SUFFIX", None)
        e.update(env or {})
        return subprocess.run([sys.executable, "-I", p, *args],
                              capture_output=True, text=True, cwd=d, env=e)

    # 1. unsuffixed, --out omitted -> the ORIGINAL names, unchanged behaviour
    r = run()
    assert r.returncode == 0, r.stdout + r.stderr
    assert "OUT summaries_dedup.txt" in r.stdout, r.stdout
    assert "WRITE colorado_internal_control_dedup.tsv\n" in r.stdout, r.stdout
    assert "DNAME direction_concordance_dedup.tsv" in r.stdout, r.stdout
    assert "-> colorado_internal_control_dedup.tsv\n" in r.stdout, r.stdout

    # 2. suffixed consistently -> every name moves, including the echo
    r2 = run(["--discovery", "{ET}_v2filtered_rl125"],
             {"TBL_SUFFIX": "_rl125"})
    assert r2.returncode == 0, r2.stdout + r2.stderr
    assert "OUT summaries_dedup_rl125.txt" in r2.stdout, r2.stdout
    assert "WRITE colorado_internal_control_dedup_rl125.tsv" in r2.stdout, r2.stdout
    assert "DNAME direction_concordance_dedup_rl125.tsv" in r2.stdout, r2.stdout
    assert "-> colorado_internal_control_dedup_rl125.tsv" in r2.stdout, r2.stdout

    # 3. THE RUN THAT CAUSED THE DEFECT must now refuse:
    #    a read-length --discovery with TBL_SUFFIX unset
    r3 = run(["--discovery", "{ET}_v2filtered_rl125",
              "--out", "summaries_dedup_rl125.txt"])
    assert r3.returncode != 0, r3.stdout + r3.stderr
    both3 = r3.stdout + r3.stderr
    assert "REFUSING" in both3 and "TBL_SUFFIX is unset" in both3, both3
    assert "TBL_SUFFIX=_rl125" in both3, both3

    # 4. TBL_SUFFIX set but --discovery not suffixed -> refuse
    r4 = run([], {"TBL_SUFFIX": "_rl125"})
    assert r4.returncode != 0, r4.stdout + r4.stderr
    both4 = r4.stdout + r4.stderr
    # NB the message wraps, so match a fragment that cannot straddle a newline
    assert "REFUSING" in both4 and "not end with it" in both4, both4

    # 5. --out explicitly given WITHOUT the suffix -> refuse
    r5 = run(["--discovery", "{ET}_v2filtered_rl125", "--out", "summaries.txt"],
             {"TBL_SUFFIX": "_rl125"})
    assert r5.returncode != 0, r5.stdout + r5.stderr
    assert "does not contain the suffix" in (r5.stdout + r5.stderr), r5.stderr

    # 6. a non-default stem without a read-length suffix is still allowed —
    #    sensitivity runs on {ET}_v2compadj must not be blocked
    r6 = run(["--discovery", "{ET}_v2compadj"])
    assert r6.returncode == 0, r6.stdout + r6.stderr
    assert "DISC {ET}_v2compadj" in r6.stdout, r6.stdout
    r7 = run(["--discovery", "{ET}_v2compadj_rl125"], {"TBL_SUFFIX": "_rl125"})
    assert r7.returncode == 0, r7.stdout + r7.stderr
    assert "DNAME direction_concordance_dedup_rl125.tsv" in r7.stdout, r7.stdout

    # 7. double application refused, text unchanged
    n8, p8 = apply_edits(new)
    assert p8, "second application was not refused"
    assert n8 == new, "a refused re-application still altered the text"

    print("[ok] self-test passed: both hard-coded names and both echo lines "
          "follow SUF; --out is derived when omitted; the 11:53 run "
          "(read-length --discovery, TBL_SUFFIX unset) now REFUSES; so do a "
          "suffix/--discovery mismatch and a suffix/--out mismatch; a "
          "non-read-length stem such as {ET}_v2compadj is still allowed; "
          "unsuffixed behaviour is byte-unchanged; double application refused")
    return 0


if __name__ == "__main__":
    sys.exit(main())
