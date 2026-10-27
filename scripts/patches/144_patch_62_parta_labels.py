#!/usr/bin/env python3
"""
144_patch_62_parta_labels.py -- S5(a): the JC and JCEC columns are swapped.

THE DEFECT
----------
`62_build_s5.py` part (a) does:

    jc = first_existing(f"{base}_dedup.all_events.tsv",
                        f"{base}.all_events.tsv")          # labelled JC
    cands = glob(f"*{coh}*JCEC*.tsv") ...                  # sought as JCEC

The plain/`_dedup` cohort tables come from phase 2 on
`rmats_<COH>/<ET>.MATS.JCEC.txt` — **they are the JCEC side**, and they are
what 45 computed the manuscript's 0.521 / 0.476 / 0.386 / 0.243 from. The real
JC side is script 20's output in `jccheck<SUF>/`, whose filenames contain
`_JC` and never `JCEC`, so the glob could never match and part (a) has never
built — at either read length.

**The gate is the dangerous part.** It reads

    bad = [... for c in res if abs(res[c]["JC"] - MANUSCRIPT_SLOPE[c]) > TOL]

and would have PASSED, because the column called "JC" really is the
manuscript's quantification. A gate that passes on a mislabelled column is
worse than one that fails: it certifies the mistake. Nothing wrong has been
published only because the JCEC side was never found.

62's own error text shows the same confusion — "Script 20 produced them"
refers to tables 20 does not produce — and its fallback points at raw
`SE.MATS.JCEC.txt`, suggesting a quantification job that is not needed,
because the JCEC side is already on disk as those `all_events` tables.

WHAT THIS CHANGES
-----------------
  * `jcec` <- the plain/_dedup cohort tables (what was called `jc`)
  * `jc`   <- `jccheck<SUF>/{coh}_SE_HCM_JC[_dedup].all_events.tsv`
  * both sides prefer a `_dedup` variant where one exists, so the two are on
    the same exclusion basis
  * the gate moves to the JCEC column, where the manuscript's slopes live
  * the printed table and the emitted rows swap to match
  * missing cohorts are reported as SCOPE, not as a failure

THE DEDUP BASIS — WHY STANFORD NEEDED AN EXTRA TABLE
-----------------------------------------------------
`exclude_dups.ids` holds 27 ids, of which exactly ONE is in any of the three
cohorts' bam lists: a Stanford control. So `Stanford_SE_HCM_dedup.all_events.tsv`
exists and Brigham's and Fuwai's do not. Without a matching dedup JC table,
Stanford alone would have compared a dedup JCEC side against a non-dedup JC
side. `jccheck/Stanford_SE_HCM_JC_dedup.all_events.tsv` and its `_rl125` twin
were built on 21 Oct for exactly this, and verified byte-identical to each
other (they must be: `rmats_Stanford/` never moved).

SCOPE: THREE COHORTS
--------------------
Script 20's part C is `for COH in Brigham Fuwai Stanford`. **Colorado has no
JC table** — it is the DCM-vs-NF contrast and sits outside 20's HCM loop. Part
(a) therefore covers three of four cohorts. Rather than leave a blank row,
which 62 forbids, it now prints the scope, records it in `notes`, and emits
rows only for the cohorts it has.

WHAT (a) WILL REPORT — the conclusion does not survive
-------------------------------------------------------
Script 20's own 125 bp run (job 409724) already gives the answer:

    arm               JCEC->JC at 50 bp      at 125 bp
    MAGNet disjoint     -0.006                 -0.021
    Brigham             +0.061                 +0.026
    Fuwai               +0.036                 -0.007
    Stanford            +0.035                 +0.009

The 50 bp claim was "the internal control barely moves, so the external change
is a genuine read-length effect of ~0.05". At 125 bp the internal arm moves
3.5x MORE and the external shifts collapse and change sign. **S5(a)'s sentence
needs rewriting, not re-basing.**

USAGE
-----
  python 144_patch_62_parta_labels.py [--apply] [--self-test]
"""
import argparse
import difflib
import os
import shutil
import sys
import time

TARGET = "62_build_s5.py"

# ------------------------------------------------------------------- edit 1
A_SRC = '''        jc = first_existing(f"{base}_dedup.all_events.tsv", f"{base}.all_events.tsv")
        # JCEC tables: script 20 wrote them, naming unknown -- discover rather than guess
        cands = sorted(set(glob.glob(f"*{coh}*JCEC*.tsv") + glob.glob(f"*JCEC*{coh}*.tsv")
                           + glob.glob(f"{base}*JCEC*.tsv")))
        cands = [c for c in cands if et in c or coh in c]
        jcec = next((c for c in cands if "dedup" in c), cands[0] if cands else None)
        print(f"  {coh:<10} JC   -> {jc or '*** MISSING ***'}")
        print(f"  {'':<10} JCEC -> {jcec or '*** MISSING -- see note below ***'}")
        if cands and len(cands) > 1:
            print(f"  {'':<10}      ({len(cands)} JCEC candidates: {', '.join(cands[:4])})")
        out[coh] = (jc, jcec)
'''
B_SRC = '''        # THE LABELS WERE SWAPPED UNTIL 21 Oct. The plain/_dedup cohort
        # tables come from phase 2 on rmats_<COH>/<ET>.MATS.JCEC.txt -- they
        # ARE the JCEC side, and they are what 45 computed the manuscript's
        # slopes from. The JC side is script 20's output in jccheck<SUF>/,
        # whose names contain "_JC" and never "JCEC", which is why the old
        # glob could never match and (a) never built.
        #
        # Both sides prefer a _dedup variant where one exists, so they sit on
        # the same exclusion basis. Only Stanford has one: exclude_dups.ids
        # contains exactly one sample from any of these three cohorts.
        jcec = first_existing(f"{base}_dedup.all_events.tsv",
                              f"{base}.all_events.tsv")
        jc = first_existing(
            f"jccheck{SUF}/{coh}_{et}_HCM_JC_dedup.all_events.tsv",
            f"jccheck{SUF}/{coh}_{et}_HCM_JC.all_events.tsv")
        print(f"  {coh:<10} JCEC -> {jcec or '*** MISSING ***'}")
        print(f"  {'':<10} JC   -> {jc or '*** MISSING (script 20 writes these) ***'}")
        if jcec and jc and ("dedup" in jcec) != ("dedup" in jc):
            print(f"  {'':<10}      ⚠ ONE SIDE IS _dedup AND THE OTHER IS NOT --")
            print(f"  {'':<10}        the comparison would mix exclusion bases")
        out[coh] = (jc, jcec)
'''

# ------------------------------------------------------------------- edit 2
A_NB = '''    if any(v[1] is None for v in out.values()):
        print("\\n  >>> (a) NOT BUILT. No JCEC tables found for at least one cohort.")
        print("      Script 20 produced them; they are not in this directory under any")
        print("      name containing 'JCEC'. Either re-run 20_jc_vs_jcec_readlength.sbatch")
        print("      against the _dedup tables, or pass --jcec-glob with the real pattern.")
        raw = sorted(glob.glob("**/SE.MATS.JCEC.txt", recursive=True))[:12]
        if raw:
            print("\\n      RAW rMATS JCEC output DOES exist, here:")
            for r in raw:
                print(f"        {r}")
            print("      So (a) needs phase 2 run over these to produce the")
            print("      `all_events` tables with adjusted dPSI and standard errors --")
            print("      a job, not a lookup. Deming needs a per-event SE, which raw")
            print("      rMATS does not provide.")
        else:
            print("\\n      No SE.MATS.JCEC.txt anywhere either. rMATS was not run with")
            print("      JCEC output for these cohorts, so (a) needs the quantification")
            print("      step too, not just phase 2.")
        notes.append("(a) NOT BUILT -- JCEC tables not found")
        return

    res = {}
    for coh, (jc, jcec) in out.items():
'''
B_NB = '''    have = {c: v for c, v in out.items() if v[0] and v[1]}
    gone = [c for c in out if c not in have]
    if not have:
        print("\\n  >>> (a) NOT BUILT. No cohort has BOTH sides.")
        print("      The JC side comes from script 20, which writes")
        print(f"      jccheck{SUF}/<COH>_SE_HCM_JC.all_events.tsv. Run:")
        print(f"        TBL_SUFFIX={SUF} sbatch 20_jc_vs_jcec_readlength.sbatch")
        print("      The JCEC side is the plain/_dedup cohort table and should")
        print("      already exist. Do NOT re-quantify: raw SE.MATS.JCEC.txt is")
        print("      not what is needed -- Deming needs the per-event SE that")
        print("      phase 2 adds, and those tables are already on disk.")
        notes.append("(a) NOT BUILT -- no cohort has both sides")
        return
    if gone:
        print(f"\\n  ⚠ SCOPE: no JC table for {', '.join(gone)}.")
        print("    Script 20's part C covers Brigham, Fuwai and Stanford only;")
        print("    Colorado is the DCM-vs-NF contrast and sits outside its HCM")
        print("    loop. (a) therefore covers "
              f"{len(have)} of {len(out)} cohorts.")
        print("    THE CAPTION MUST SAY SO. 62's own rule is that a supplement")
        print("    must not carry a blank row, so these cohorts are omitted")
        print("    rather than left empty.")
        notes.append(f"(a) covers {len(have)}/{len(out)} cohorts -- no JC table "
                     f"for {', '.join(gone)}")

    res = {}
    for coh, (jc, jcec) in have.items():
'''

# ------------------------------------------------------------------- edit 3
A_GATE = '''    # THE CHECK: the JC column must be the manuscript's slopes
    bad = [(c, res[c]["JC"], MANUSCRIPT_SLOPE[c]) for c in res
           if np.isfinite(res[c]["JC"])
           and abs(res[c]["JC"] - MANUSCRIPT_SLOPE[c]) > SLOPE_TOL]
    print()
    print(f"  {'cohort':<10}{'JC':>9}{'manuscript':>12}{'JCEC':>9}{'change':>9}")
    for c in ("Colorado", "Brigham", "Fuwai", "Stanford"):
        r = res.get(c, {})
        f = lambda v: f"{v:.3f}" if np.isfinite(v) else "-"
        ch = (r.get("JCEC", np.nan) - r.get("JC", np.nan))
        print(f"  {c:<10}{f(r.get('JC', np.nan)):>9}{MANUSCRIPT_SLOPE[c]:>12.3f}"
              f"{f(r.get('JCEC', np.nan)):>9}{(f'{ch:+.3f}' if np.isfinite(ch) else '-'):>9}")
    if bad:
        print("\\n  >>> (a) NOT BUILT. The recomputed JC slopes do not reproduce the")
        print("      manuscript's:")
        for c, got, want in bad:
            print(f"        {c}: {got:.3f} here, {want:.3f} in the Results")
        print("      Until the JC column matches Figure 4, the JCEC column beside it is")
        print("      not computed on the same footing and the comparison means nothing.")
        notes.append("(a) NOT BUILT -- JC column does not reproduce the manuscript")
        return
    print("\\n  [ok] the JC column reproduces the manuscript's slopes; the JCEC column")
    print("       is therefore computed on the same event set and exclusions")
    for c in ("Colorado", "Brigham", "Fuwai", "Stanford"):
        rows.append(dict(part="a. JC vs JCEC", row=c,
                         value_1=round(res[c]["JC"], 3), label_1="Deming, JC",
                         value_2=round(res[c]["JCEC"], 3), label_2="Deming, JCEC",
                         change=round(res[c]["JCEC"] - res[c]["JC"], 3)))
'''
B_GATE = '''    # THE CHECK: the JCEC column must be the manuscript's slopes.
    # 45 computed 0.521 / 0.476 / 0.386 / 0.243 from the _dedup JCEC tables,
    # so JCEC is the column that has to reproduce them. Until 21 Oct this
    # gate read res[c]["JC"] -- and would have PASSED, because the column
    # then called "JC" was in fact the JCEC table. A gate that passes on a
    # mislabelled column certifies the mistake.
    bad = [(c, res[c]["JCEC"], MANUSCRIPT_SLOPE[c]) for c in res
           if np.isfinite(res[c]["JCEC"])
           and abs(res[c]["JCEC"] - MANUSCRIPT_SLOPE[c]) > SLOPE_TOL]
    print()
    print(f"  {'cohort':<10}{'JCEC':>9}{'manuscript':>12}{'JC':>9}{'change':>9}")
    for c in ("Colorado", "Brigham", "Fuwai", "Stanford"):
        r = res.get(c, {})
        f = lambda v: f"{v:.3f}" if np.isfinite(v) else "-"
        ch = (r.get("JC", np.nan) - r.get("JCEC", np.nan))
        print(f"  {c:<10}{f(r.get('JCEC', np.nan)):>9}{MANUSCRIPT_SLOPE[c]:>12.3f}"
              f"{f(r.get('JC', np.nan)):>9}{(f'{ch:+.3f}' if np.isfinite(ch) else '-'):>9}")
    if bad:
        print("\\n  >>> (a) NOT BUILT. The recomputed JCEC slopes do not reproduce")
        print("      the manuscript's:")
        for c, got, want in bad:
            print(f"        {c}: {got:.3f} here, {want:.3f} in the Results")
        print("      Until the JCEC column matches Figure 4, the JC column beside it")
        print("      is not computed on the same footing and the comparison means")
        print("      nothing.")
        notes.append("(a) NOT BUILT -- JCEC column does not reproduce the manuscript")
        return
    print("\\n  [ok] the JCEC column reproduces the manuscript's slopes; the JC")
    print("       column is therefore computed on the same event set and exclusions")
    print("  ⚠ READ THE SIGN AND SIZE OF `change` AGAINST THE MAGNet-INTERNAL ARM")
    print("    in script 20's own jc_vs_jcec<SUF>.tsv before writing a sentence.")
    print("    At 50 bp the internal arm moved -0.006 while the cohorts moved")
    print("    ~+0.04, which is what 'read length is largely excluded' rested on.")
    print("    At 125 bp the internal arm moves -0.021 and the cohorts +0.026 /")
    print("    -0.007 / +0.009 -- the premise does not hold and the sentence")
    print("    needs rewriting, not re-basing.")
    for c in sorted(res):
        rows.append(dict(part="a. JC vs JCEC", row=c,
                         value_1=round(res[c]["JCEC"], 3), label_1="Deming, JCEC",
                         value_2=round(res[c]["JC"], 3), label_2="Deming, JC",
                         change=round(res[c]["JC"] - res[c]["JCEC"], 3)))
'''

# ------------------------------------------------------------------- edit 4
A_D1 = """`RESULT_magnitude_resolved.md` already holds a JC-versus-JCEC table, and it CONTRADICTS
the manuscript: it gives Brigham 0.757 where the Results and Figure 4 say 0.476. That
table is pre-exclusion and JCEC-based; the manuscript uses JC on the post-exclusion
`_dedup` tables. Putting it in the supplement would have published two different slopes
for the same cohort.
"""
B_D1 = """`RESULT_magnitude_resolved.md` already holds a JC-versus-JCEC table, and it CONTRADICTS
the manuscript: it gives Brigham 0.757 where the Results and Figure 4 say 0.476.

⚠ CORRECTED 21 Oct. An earlier version of this paragraph said "the manuscript uses JC
on the post-exclusion `_dedup` tables". That is backwards. The manuscript is **JCEC**:
04, 05, 08 and 15b all run phase 2 on `<dir>/<ET>.MATS.JCEC.txt`, and 03/03b's use of
`.MATS.JC.txt` is a documented exception for the forward RBM20 target test alone. The
0.757 differs because that table is **pre-exclusion** and computed on the matched
JC-intersect-JCEC subset, not because of JC versus JCEC -- the same document records
that restricting to the matched subset alone raises Brigham from 0.662 to 0.757.
Putting it in the supplement would have published two different slopes for one cohort.
"""

A_D2 = ("  (a) JC slopes must match 0.521 / 0.476 / 0.386 / 0.243 within 0.02\n")
B_D2 = ("  (a) JCEC slopes must match 0.521 / 0.476 / 0.386 / 0.243 within 0.02\n"
        "      (JCEC, not JC: the manuscript's tables are JCEC. The JC side comes\n"
        "       from script 20's jccheck<SUF>/ output.)\n")

EDITS = [
    ("the JC/JCEC sources, swapped", A_SRC, B_SRC, 1),
    ("the NOT-BUILT block and the scope note", A_NB, B_NB, 1),
    ("the gate, moved to the JCEC column", A_GATE, B_GATE, 1),
    ("the header's provenance paragraph", A_D1, B_D1, 1),
    ("the self-check list", A_D2, B_D2, 1),
]


def apply_edits(text):
    if "THE LABELS WERE SWAPPED UNTIL 21 Oct" in text:
        return text, ["62 is already patched (restore the .bak to re-apply)"]
    if 'SUF = os.environ.get("TBL_SUFFIX", "")' not in text:
        return text, ["141 has not been applied — apply 141_patch_62_s5.py "
                      "first; this uses SUF to find jccheck<SUF>/"]
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
    for w in ("THE LABELS WERE SWAPPED UNTIL 21 Oct",
              'jcec = first_existing(f"{base}_dedup.all_events.tsv"',
              'f"jccheck{SUF}/{coh}_{et}_HCM_JC_dedup.all_events.tsv"',
              'f"jccheck{SUF}/{coh}_{et}_HCM_JC.all_events.tsv"',
              "ONE SIDE IS _dedup AND THE OTHER IS NOT",
              "have = {c: v for c, v in out.items() if v[0] and v[1]}",
              "for coh, (jc, jcec) in have.items():",
              'bad = [(c, res[c]["JCEC"], MANUSCRIPT_SLOPE[c]) for c in res',
              'if np.isfinite(res[c]["JCEC"])',
              "⚠ SCOPE: no JC table for",
              "READ THE SIGN AND SIZE OF `change`",
              "for c in sorted(res):",
              "(a) JCEC slopes must match",
              "CORRECTED 21 Oct"):
        if w not in text:
            bad.append(f"missing: {w[:58]}")
    for gone in ('jc = first_existing(f"{base}_dedup.all_events.tsv", '
                 'f"{base}.all_events.tsv")',
                 'glob.glob(f"*{coh}*JCEC*.tsv")',
                 'bad = [(c, res[c]["JC"], MANUSCRIPT_SLOPE[c]) for c in res',
                 "the manuscript uses JC on the post-exclusion",
                 "(a) JC slopes must match",
                 "for coh, (jc, jcec) in out.items():"):
        if gone in text:
            bad.append(f"stale site survives: {gone[:52]}")
    # the final emit must iterate res, never the hard-coded four, or a
    # three-cohort run raises KeyError on Colorado
    i_emit = text.find('rows.append(dict(part="a. JC vs JCEC"')
    if i_emit < 0:
        bad.append("the (a) row emit was removed")
    else:
        head = text[:i_emit].rsplit("\n", 2)[-2:]
        if not any("for c in sorted(res):" in h for h in head):
            bad.append("the (a) row emit does not iterate `res` — a missing "
                       "cohort would raise KeyError")
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
    print("144 — S5(a): the JC and JCEC columns were swapped")
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
        print("  # 1. 50 bp — (a) should BUILD for the first time")
        print("  python3 62_build_s5.py --parts a 2>&1 | tail -30")
        print()
        print("  # 2. 125 bp")
        print("  TBL_SUFFIX=_rl125 \\")
        print("  S5_SLOPE=Colorado=0.5999,Brigham=0.5099,"
              "Stanford=0.4008,Fuwai=0.2948 \\")
        print("  S5_CONC=Colorado=80.8,Brigham=81.9,Fuwai=70.8,Stanford=75.8 \\")
        print("  S5_DISC_SE=exon=813,full=944 \\")
        print("  python3 62_build_s5.py --parts a 2>&1 | tail -30")
        print()
        print("⚠ AT 50 bp THE GATE SHOULD PASS: the JCEC column is the")
        print("  manuscript's own quantification, so it must reproduce")
        print("  0.476 / 0.243 / 0.386 for Brigham / Fuwai / Stanford. If it")
        print("  does NOT, the labels are still wrong somewhere and nothing")
        print("  in (a) is usable.")
        print()
        print("⚠ AT 125 bp the gate compares against S5_SLOPE's 125 bp values.")
        print("  Colorado will be absent from (a) — three cohorts, stated.")
        print()
        print("⬜ THEN THE SENTENCE. Script 20's own run already shows the")
        print("  50 bp claim fails at 125 bp: internal arm -0.006 -> -0.021,")
        print("  cohorts ~+0.04 -> +0.026 / -0.007 / +0.009. S5(a) is a")
        print("  rewrite, not a re-basing.")
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
    print("✓ verified on disk: JCEC is the plain/_dedup table, JC comes from")
    print("  jccheck<SUF>/, both prefer _dedup, the gate is on JCEC, a")
    print("  mixed-basis pair warns, and a missing cohort is scope not failure.")
    return 0


def self_test():
    import subprocess
    import tempfile

    fx = ('"""\n' + A_D1 + "\n" + A_D2 + '"""\n'
          "import glob, os\n"
          "import numpy as np\n"
          'SUF = os.environ.get("TBL_SUFFIX", "")\n'
          'MANUSCRIPT_SLOPE = {"Colorado": 0.521, "Brigham": 0.476,\n'
          '                    "Stanford": 0.386, "Fuwai": 0.243}\n'
          "SLOPE_TOL = 0.02\n"
          'REPL = [("Colorado", "{ET}_Colorado_DCMvsNF"),\n'
          '        ("Brigham", "Brigham_{ET}_HCM"),\n'
          '        ("Fuwai", "Fuwai_{ET}_HCM"),\n'
          '        ("Stanford", "Stanford_{ET}_HCM")]\n'
          "def first_existing(*paths):\n"
          "    for p in paths:\n"
          "        if os.path.exists(p) and os.path.getsize(p):\n"
          "            return p\n"
          "    return None\n"
          "def read_events(p, et):\n"
          "    return ({'x': p}, None) if p else (None, 'missing')\n"
          "def part_a(a, rows, notes):\n"
          '    et = "SE"\n'
          "    sig = None\n"
          "    out = {}\n"
          "    for coh, pat in REPL:\n"
          '        base = pat.replace("{ET}", et)\n'
          + A_SRC
          + A_NB
          + "        r = {}\n"
            '        for lab, p in (("JC", jc), ("JCEC", jcec)):\n'
            "            r[lab] = {'Brigham': 0.476, 'Fuwai': 0.243,\n"
            "                      'Stanford': 0.386, 'Colorado': 0.521}[coh] \\\n"
            "                     + (0.03 if lab == 'JC' else 0.0)\n"
            "        res[coh] = r\n"
          + A_GATE
          + "class _A:\n"
            "    list = False\n"
            "rows, notes = [], []\n"
            "part_a(_A(), rows, notes)\n"
            'print("ROWS", [(r["row"], r["label_1"], r["value_1"], '
            'r["label_2"], r["value_2"]) for r in rows])\n'
            'print("NOTES", notes)\n')

    new, probs = apply_edits(fx)
    assert not probs, probs
    assert not verify(new), verify(new)
    compile(new, "fx", "exec")

    d = tempfile.mkdtemp()
    p = os.path.join(d, "fx.py")
    open(p, "w").write(new)
    os.makedirs(os.path.join(d, "jccheck"), exist_ok=True)
    os.makedirs(os.path.join(d, "jccheck_rl125"), exist_ok=True)

    def touch(rel):
        fp = os.path.join(d, rel)
        os.makedirs(os.path.dirname(fp), exist_ok=True)
        open(fp, "w").write("x\n")

    # the JCEC side for all four, the JC side for three (no Colorado)
    for c in ("Brigham", "Fuwai", "Stanford"):
        touch(f"{c}_SE_HCM.all_events.tsv")
        touch(f"jccheck/{c}_SE_HCM_JC.all_events.tsv")
        touch(f"jccheck_rl125/{c}_SE_HCM_JC.all_events.tsv")
    touch("SE_Colorado_DCMvsNF.all_events.tsv")
    # Stanford alone has a _dedup pair, as on the real cluster
    touch("Stanford_SE_HCM_dedup.all_events.tsv")
    touch("jccheck/Stanford_SE_HCM_JC_dedup.all_events.tsv")
    touch("jccheck_rl125/Stanford_SE_HCM_JC_dedup.all_events.tsv")

    def run(env=None):
        e = dict(os.environ)
        e.pop("TBL_SUFFIX", None)
        e.update(env or {})
        return subprocess.run([sys.executable, "-I", p], capture_output=True,
                              text=True, cwd=d, env=e)

    r = run()
    assert r.returncode == 0, r.stdout + r.stderr
    o = r.stdout
    # JCEC is the plain/_dedup table; JC comes from jccheck/
    assert "Brigham    JCEC -> Brigham_SE_HCM.all_events.tsv" in o, o
    assert "JC   -> jccheck/Brigham_SE_HCM_JC.all_events.tsv" in o, o
    # Stanford picks the _dedup pair on BOTH sides, so no mixed-basis warning
    assert "Stanford_SE_HCM_dedup.all_events.tsv" in o, o
    assert "jccheck/Stanford_SE_HCM_JC_dedup.all_events.tsv" in o, o
    assert "ONE SIDE IS _dedup" not in o, o
    # Colorado has no JC table -> scope, not failure
    assert "SCOPE: no JC table for Colorado" in o, o
    assert "covers 3 of 4 cohorts" in o, o
    assert "NOT BUILT" not in o, o
    # the gate is on JCEC and passes; the emitted rows exclude Colorado
    assert "[ok] the JCEC column reproduces" in o, o
    assert "'Deming, JCEC'" in o or "Deming, JCEC" in o, o
    # the ROWS line ALONE -- the NOTES line legitimately names Colorado
    # (that is the scope note), so splitting on "ROWS" and taking everything
    # after it tests the wrong text.
    _rowline = [ln for ln in o.splitlines() if ln.startswith("ROWS ")]
    assert len(_rowline) == 1, o
    _rowline = _rowline[0]
    assert "Colorado" not in _rowline, _rowline
    for c in ("Brigham", "Fuwai", "Stanford"):
        assert c in _rowline, (c, _rowline)
    # and the scope note IS emitted, on the notes line
    _noteline = [ln for ln in o.splitlines() if ln.startswith("NOTES ")]
    assert len(_noteline) == 1 and "Colorado" in _noteline[0], o
    # the change column is JC - JCEC, i.e. +0.03 in the fixture
    assert "READ THE SIGN AND SIZE" in o, o

    # suffixed: the JC side must come from jccheck_rl125/
    r2 = run({"TBL_SUFFIX": "_rl125"})
    assert r2.returncode == 0, r2.stdout + r2.stderr
    assert "jccheck_rl125/Brigham_SE_HCM_JC.all_events.tsv" in r2.stdout, r2.stdout
    assert "jccheck/Brigham" not in r2.stdout, r2.stdout

    # a mixed basis must WARN: give Fuwai a _dedup JCEC and no _dedup JC
    touch("Fuwai_SE_HCM_dedup.all_events.tsv")
    r3 = run()
    assert "ONE SIDE IS _dedup AND THE OTHER IS NOT" in r3.stdout, r3.stdout
    os.remove(os.path.join(d, "Fuwai_SE_HCM_dedup.all_events.tsv"))

    # no JC tables at all -> NOT BUILT, naming script 20
    for c in ("Brigham", "Fuwai", "Stanford"):
        os.remove(os.path.join(d, f"jccheck/{c}_SE_HCM_JC.all_events.tsv"))
    os.remove(os.path.join(d, "jccheck/Stanford_SE_HCM_JC_dedup.all_events.tsv"))
    r4 = run()
    assert "NOT BUILT. No cohort has BOTH sides" in r4.stdout, r4.stdout
    assert "20_jc_vs_jcec_readlength.sbatch" in r4.stdout, r4.stdout
    assert "Do NOT re-quantify" in r4.stdout, r4.stdout

    # refuses without 141
    _n, p5 = apply_edits("print('bare')\n")
    assert p5 and "141 has not been applied" in p5[0], p5

    # double application refused
    n6, p6 = apply_edits(new)
    assert p6, "second application was not refused"
    assert n6 == new, "a refused re-application still altered the text"

    print("[ok] self-test passed: JCEC resolves to the plain/_dedup cohort "
          "table and JC to jccheck<SUF>/, with Stanford taking the _dedup pair "
          "on both sides and no mixed-basis warning; a deliberately mixed pair "
          "DOES warn; Colorado's missing JC table is reported as scope with "
          "'covers 3 of 4' and omitted from the emitted rows rather than "
          "raising KeyError; the gate is on JCEC and passes; TBL_SUFFIX moves "
          "the JC side to jccheck_rl125/; no JC tables at all gives NOT BUILT "
          "naming script 20 and warning against re-quantifying; refuses "
          "without 141; double application refused")
    return 0


if __name__ == "__main__":
    sys.exit(main())
