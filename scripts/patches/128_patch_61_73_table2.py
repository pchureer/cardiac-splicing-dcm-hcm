#!/usr/bin/env python3
"""
128_patch_61_73_table2.py -- rebuild Table 2 at 125 bp.

WHY REBUILD RATHER THAN EDIT
----------------------------
Table 2 is **derived, not curated**. `61_table2_lead_genes.py` line 283 is
`.head(a.top)` with `--top` defaulting to 20, and its docstring fixes the rule:

  1. events passing FDR < 0.05 AND |dPSI| >= 0.10, exon-level
  2. one row per GENE, that gene's largest-|dPSI| significant event
  3. rank by number of independent replication cohorts agreeing in direction,
     then by |dPSI| in discovery

So CSRP3, STK39 and ECHDC1 do not "drop out of Table 2" — they drop out of the
**eligible pool**, and the same rule promotes whoever is next. The table still
has twenty rows. `DECISION_rl125_verdict_confirmed.md`'s objection — "the
replacements would come from the new selection set while the counts came from
the old" — was an argument against fixing Table 2 while 50 bp stayed primary.
With 125 bp primary, rebuilding is the consistent move.

THE CHAIN
---------
    61  ->  Table2.tsv  ->  73  ->  Table2_recount.tsv

73 reads 61's output, so the suffix must flow through both or 73 validates the
125 bp lead genes against the 50 bp tables.

WHAT MOVES, AND WHAT DOES NOT
-----------------------------
61's cohort self-check **passes unchanged**. `KNOWN_PCT` is
Colorado 79.6 / Brigham 81.2 / Fuwai 69.2 / Stanford 76.1 with `PCT_TOL = 4.0`,
and 33's 125 bp rerun gave 80.8 / 81.9 / 70.8 / 75.8 — every one inside the
tolerance. It is parameterised anyway, because the measured values now exist and
a gate should check what it claims to check.

`FROZEN_SE = 926` is the only check in 61 that aborts: 813 at 125 bp.

73's constants, measured (`RESULT_s1_composition_rl125_complete.md`):

    REC_BASE_SIG      926  ->  813
    REC_COMPADJ_KEPT  590  ->  498      (nine-marker retention, 61.3% of 813)
    REC_LEAD_KEPT      18  ->  UNKNOWN — this is the number being measured

`REC_LEAD_KEPT=0` therefore means "no expectation; report it". The companion
check — Table2's own `comp_adj` flag against the recomputation — still fires,
because that is an internal-consistency test and does not depend on a recorded
value.

THE BINOMIAL IS ALREADY DROPPED AND 125 bp DOES NOT REVIVE IT
-------------------------------------------------------------
73's docstring records 18 of 20 against a 77.3% background giving *P* = 0.1354.
At 125 bp the backgrounds barely move — nine-marker 63.7% -> **61.3%**
(498/813), fraction-only 77.3% -> **77.1%** (627/813) — so that decision stands
on its own terms whatever the recount returns.

⚠ NOT GUARDED, BY CHOICE
------------------------
61 writes `<prefix>.{tsv,md,xlsx}`, and the suffix is what keeps the 125 bp run
off the published `Table2.tsv`. There is no don't-clobber guard here: adding one
blind would mean patching write sites this patch has not read. **Do not run 61
unsuffixed** unless you mean to regenerate the published table.

USAGE
-----
  python 128_patch_61_73_table2.py
  python 128_patch_61_73_table2.py --apply
  python 128_patch_61_73_table2.py --self-test
"""
import argparse
import difflib
import os
import shutil
import sys
import time

HEAD_61 = '''
# --------------------------------------------------------------------------
# Read length. Empty = the original 50 bp behaviour, reproducing the published
# Table2.{tsv,md,xlsx}; _rl125 = the re-quantified tables, written to
# Table2_rl125.*.
#
# Table 2 is DERIVED (line ~283 is .head(--top)), so a read-length change is a
# rebuild, not an edit: three genes leave the eligible pool and the same rule
# promotes whoever is next.
SUF = os.environ.get("TBL_SUFFIX", "")

# FROZEN_SE is the only self-check here that aborts. 926 at 50 bp, 813 at 125.
FROZEN_SE = int(os.environ.get("FROZEN_SE", "926"))

# KNOWN_PCT passes unchanged at 125 bp (measured 80.8 / 81.9 / 70.8 / 75.8
# against PCT_TOL 4.0), but is overridable so the gate checks what it claims:
#     KNOWN_PCT=Colorado=80.8,Brigham=81.9,Fuwai=70.8,Stanford=75.8
_kp = os.environ.get("KNOWN_PCT", "").strip()
# --------------------------------------------------------------------------
'''

TAIL_61 = '''
if _kp:
    _d = dict(_p.split("=", 1) for _p in _kp.split(","))
    assert set(_d) == set(KNOWN_PCT), (
        f"KNOWN_PCT override must name exactly {sorted(KNOWN_PCT)}, "
        f"got {sorted(_d)}")
    KNOWN_PCT = {_k: float(_v) for _k, _v in _d.items()}
'''

HEAD_73 = '''# --------------------------------------------------------------------------
# Read length. The chain is 61 -> Table2<SUF>.tsv -> 73, so the suffix must
# flow through both or this script validates the 125 bp lead genes against the
# 50 bp tables.
SUF = os.environ.get("TBL_SUFFIX", "")
# --------------------------------------------------------------------------
'''

TAIL_73 = '''
# Measured at 125 bp (RESULT_s1_composition_rl125_complete.md):
#     REC_BASE_SIG=813 REC_COMPADJ_KEPT=498
# REC_LEAD_KEPT is the quantity being MEASURED, so 0 means "no expectation,
# report it". The Table2-flag-vs-recomputation check below still fires: that is
# internal consistency and needs no recorded value.
REC_BASE_SIG = int(os.environ.get("REC_BASE_SIG", str(REC_BASE_SIG)))
REC_COMPADJ_KEPT = int(os.environ.get("REC_COMPADJ_KEPT",
                                      str(REC_COMPADJ_KEPT)))
REC_LEAD_KEPT = int(os.environ.get("REC_LEAD_KEPT", str(REC_LEAD_KEPT)))
'''

FILES = {
    "61_table2_lead_genes.py": [
        ("SUF, FROZEN_SE and KNOWN_PCT",
         'NKEEP = {"SE": 3, "MXE": 5, "A3SS": 5, "A5SS": 5, "RI": 3}\n',
         'NKEEP = {"SE": 3, "MXE": 5, "A3SS": 5, "A5SS": 5, "RI": 3}\n'
         + HEAD_61, 1),
        ("the discovery table, suffixed",
         'DISC = "{ET}_v2filtered"',
         'DISC = "{ET}_v2filtered" + SUF', 1),
        ("the old FROZEN_SE literal",
         "FROZEN_SE = 926\n", "", 1),
        ("the KNOWN_PCT override, applied after the dict exists",
         'KNOWN_PCT = {"Colorado": 79.6, "Brigham": 81.2, "Fuwai": 69.2, '
         '"Stanford": 76.1}\n',
         'KNOWN_PCT = {"Colorado": 79.6, "Brigham": 81.2, "Fuwai": 69.2, '
         '"Stanford": 76.1}\n' + TAIL_61, 1),
        ("the composition-adjusted table, suffixed",
         '    ap.add_argument("--compadj", default="{ET}_v2compadj",',
         '    ap.add_argument("--compadj", default="{ET}_v2compadj" + SUF,', 1),
        ("the output prefix, suffixed",
         '    ap.add_argument("--prefix", default="Table2")',
         '    ap.add_argument("--prefix", default="Table2" + SUF)', 1),
    ],
    "73_table2_recount.py": [
        ("SUF", 'TABLE2 = "Table2.tsv"\n',
         HEAD_73 + 'TABLE2 = "Table2" + SUF + ".tsv"\n', 1),
        ("the base table", 'BASE = "SE_v2filtered.all_events.tsv"\n',
         'BASE = "SE_v2filtered" + SUF + ".all_events.tsv"\n', 1),
        ("the nine-marker table", 'COMPADJ = "SE_v2compadj.all_events.tsv"\n',
         'COMPADJ = "SE_v2compadj" + SUF + ".all_events.tsv"\n', 1),
        ("the fraction-only table", 'CMFRAC = "SE_v2cmfrac.all_events.tsv"\n',
         'CMFRAC = "SE_v2cmfrac" + SUF + ".all_events.tsv"\n', 1),
        ("the three recorded counts, overridable",
         "REC_LEAD_KEPT = 18\n", "REC_LEAD_KEPT = 18\n" + TAIL_73, 1),
        ("REC_LEAD_KEPT=0 means no expectation",
         "    if n_nine != REC_LEAD_KEPT:\n",
         "    if REC_LEAD_KEPT and n_nine != REC_LEAD_KEPT:\n", 1),
        ("the output name, suffixed",
         '    ap.add_argument("--out", default="Table2_recount.tsv")',
         '    ap.add_argument("--out", default="Table2_recount" + SUF '
         '+ ".tsv")', 1),
    ],
}


def apply_edits(text, edits):
    problems, out = [], text
    for name, anchor, repl, want in edits:
        got = out.count(anchor)
        if got != want:
            problems.append(f"{name}: anchor occurs {got}x, expected {want}x\n"
                            f"      anchor: "
                            f"{anchor.strip().splitlines()[0][:76] if anchor.strip() else '(empty)'}")
            continue
        if 'SUF = os.environ.get("TBL_SUFFIX", "")' in out and \
                name.startswith("SUF"):
            problems.append("SUF is already defined — already patched")
            continue
        out = out.replace(anchor, repl, want)
    return out, problems


def verify(text, fname):
    bad = []
    if 'SUF = os.environ.get("TBL_SUFFIX", "")' not in text:
        bad.append("SUF is not defined")
    if "61_" in fname:
        for w in ('DISC = "{ET}_v2filtered" + SUF',
                  'int(os.environ.get("FROZEN_SE", "926"))',
                  'default="{ET}_v2compadj" + SUF',
                  'default="Table2" + SUF',
                  'KNOWN_PCT = {_k: float(_v)'):
            if w not in text:
                bad.append(f"missing in 61: {w[:52]}")
        if text.count("FROZEN_SE = 926\n") != 0:
            bad.append("the bare FROZEN_SE literal survives")
        i_s = text.find('SUF = os.environ.get("TBL_SUFFIX"')
        for m in ('DISC = "{ET}_v2filtered" + SUF', 'default="Table2" + SUF'):
            j = text.find(m)
            if j < 0 or i_s > j:
                bad.append(f"SUF not above its use: {m[:40]}")
        # the override must come AFTER the dict it rewrites
        i_d = text.find('KNOWN_PCT = {"Colorado"')
        i_o = text.find("KNOWN_PCT = {_k: float(_v)")
        if i_d < 0 or i_o < 0 or i_d > i_o:
            bad.append("the KNOWN_PCT override is not after the literal dict")
    else:
        for w in ('TABLE2 = "Table2" + SUF + ".tsv"',
                  'BASE = "SE_v2filtered" + SUF',
                  'COMPADJ = "SE_v2compadj" + SUF',
                  'CMFRAC = "SE_v2cmfrac" + SUF',
                  'os.environ.get("REC_BASE_SIG"',
                  'os.environ.get("REC_COMPADJ_KEPT"',
                  'os.environ.get("REC_LEAD_KEPT"',
                  "if REC_LEAD_KEPT and n_nine != REC_LEAD_KEPT:",
                  'default="Table2_recount" + SUF'):
            if w not in text:
                bad.append(f"missing in 73: {w[:52]}")
        # the internal-consistency check must survive
        if "recomputation says" not in text:
            bad.append("the flag-vs-recomputation check was lost")
        for stem in ("SE_v2filtered.all_events.tsv",
                     "SE_v2compadj.all_events.tsv",
                     "SE_v2cmfrac.all_events.tsv"):
            if f'"{stem}"' in text:
                bad.append(f"an unsuffixed literal survives: {stem}")
    return bad


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--apply", action="store_true")
    ap.add_argument("--self-test", action="store_true")
    ap.add_argument("--only", choices=list(FILES))
    a = ap.parse_args()
    if a.self_test:
        return self_test()

    self_test()
    print()
    print("=" * 78)
    print("128 — rebuild Table 2 at 125 bp (61 -> Table2<SUF>.tsv -> 73)")
    print("=" * 78)

    targets = [a.only] if a.only else list(FILES)
    staged = {}
    for fname in targets:
        print(f"\n═══ {fname}")
        if not os.path.exists(fname):
            print("  [fatal] not found", file=sys.stderr)
            return 2
        src = open(fname, errors="replace").read()
        if "import os" not in src:
            print("  [fatal] no `import os`", file=sys.stderr)
            return 2
        new, problems = apply_edits(src, FILES[fname])
        if problems:
            print("  REFUSING — anchors do not match:")
            for p in problems:
                print(f"    ✗ {p}")
            return 3
        bad = verify(new, fname)
        if bad:
            print("  REFUSING — verification failed:")
            for b in bad:
                print(f"    ✗ {b}")
            return 3
        print("".join(difflib.unified_diff(
            src.splitlines(True), new.splitlines(True),
            fromfile=fname, tofile=fname + " (patched)", n=2)))
        staged[fname] = new

    print("=" * 78)
    if not a.apply:
        print("DRY RUN — nothing written. Re-run with --apply, then:")
        print()
        print("  # 1. CALIBRATION — unsuffixed must still reproduce 926 and the")
        print("  #    published twenty. --list touches nothing.")
        print("  python3 61_table2_lead_genes.py --list")
        print()
        print("  # 2. rebuild at 125 bp -> Table2_rl125.{tsv,md,xlsx}")
        print("  TBL_SUFFIX=_rl125 FROZEN_SE=813 \\")
        print("  KNOWN_PCT=Colorado=80.8,Brigham=81.9,Fuwai=70.8,Stanford=75.8 \\")
        print("  python3 61_table2_lead_genes.py")
        print()
        print("  # 3. recount -> REC_LEAD_KEPT=0 reports rather than refuses")
        print("  TBL_SUFFIX=_rl125 REC_BASE_SIG=813 REC_COMPADJ_KEPT=498 \\")
        print("  REC_LEAD_KEPT=0 python3 73_table2_recount.py")
        print()
        print("  # 4. WHICH GENES CHANGED — the thing this is all for")
        print("  cut -f1 Table2.tsv | tail -n +2 | sort > /tmp/t2_50.txt")
        print("  cut -f1 Table2_rl125.tsv | tail -n +2 | sort > /tmp/t2_125.txt")
        print("  echo '--- left the table:'; comm -23 /tmp/t2_50.txt /tmp/t2_125.txt")
        print("  echo '--- entered the table:'; comm -13 /tmp/t2_50.txt /tmp/t2_125.txt")
        print()
        print("  (check the gene column is field 1 first: head -1 Table2.tsv)")
        print()
        print("⚠ DO NOT run 61 unsuffixed without --list. The suffix is what")
        print("  keeps the 125 bp run off the published Table2.{tsv,md,xlsx};")
        print("  there is no don't-clobber guard in this patch.")
        print()
        print("⚠ Expect CSRP3, STK39, ECHDC1 to leave. Three others enter by")
        print("  the same rule. The count stays 20.")
        return 0

    for fname, new in staged.items():
        bak = f"{fname}.bak.{time.strftime('%Y%m%d_%H%M%S')}"
        shutil.copy2(fname, bak)
        open(fname, "w").write(new)
        print(f"[write] {fname}   (backup: {bak})")
        disk = open(fname, errors="replace").read()
        bad = verify(disk, fname)
        if bad or disk != new:
            print("✗ ON-DISK VERIFICATION FAILED — restoring:")
            for b in bad:
                print(f"    {b}")
            shutil.copy2(bak, fname)
            return 3
    print("✓ both verified on disk.")
    return 0


def self_test():
    import subprocess
    import tempfile

    fx61 = (
        "import os\n"
        'NKEEP = {"SE": 3, "MXE": 5, "A3SS": 5, "A5SS": 5, "RI": 3}\n'
        'DISC = "{ET}_v2filtered"                 # a trailing comment\n'
        "FROZEN_SE = 926\n"
        'KNOWN_PCT = {"Colorado": 79.6, "Brigham": 81.2, "Fuwai": 69.2, '
        '"Stanford": 76.1}\n'
        "def main():\n"
        '    ap_compadj = None\n'
        '    ap.add_argument("--compadj", default="{ET}_v2compadj",\n'
        '                    help="x")\n'
        '    ap.add_argument("--prefix", default="Table2")\n'
        "if __name__ == '__main__':\n"
        '    print("DISC", DISC.format(ET="SE"))\n'
        '    print("FROZEN_SE", FROZEN_SE)\n'
        '    print("KNOWN_PCT", KNOWN_PCT["Colorado"])\n'
        # main() must actually RUN or the add_argument stubs never print and
        # the assertions below test output that cannot exist.
        "    main()\n")

    fx73 = (
        "import os\n"
        'TABLE2 = "Table2.tsv"\n'
        'BASE = "SE_v2filtered.all_events.tsv"\n'
        'COMPADJ = "SE_v2compadj.all_events.tsv"\n'
        'CMFRAC = "SE_v2cmfrac.all_events.tsv"\n'
        "N_LEAD = 20\n"
        "REC_BASE_SIG = 926\n"
        "REC_COMPADJ_KEPT = 590\n"
        "REC_LEAD_KEPT = 18\n"
        "def check(n_nine, n_flag):\n"
        "    fail = []\n"
        "    if n_nine != REC_LEAD_KEPT:\n"
        '        fail.append("lead retention")\n'
        "    if n_flag != n_nine:\n"
        '        fail.append("recomputation says mismatch")\n'
        "    return fail\n"
        "class _AP:\n"
        "    def add_argument(self, *a, **k):\n"
        "        print('ARG', a, k.get('default'))\n"
        "ap = _AP()\n"
        "def main():\n"
        '    ap.add_argument("--out", default="Table2_recount.tsv")\n'
        "if __name__ == '__main__':\n"
        "    main()\n"
        '    print("TABLE2", TABLE2)\n'
        '    print("BASE", BASE)\n'
        '    print("RECS", REC_BASE_SIG, REC_COMPADJ_KEPT, REC_LEAD_KEPT)\n'
        '    print("FAIL18", check(18, 18))\n'
        '    print("FAIL17", check(17, 17))\n'
        '    print("FLAGMISMATCH", check(17, 16))\n')

    d = tempfile.mkdtemp()
    for fname, fx in (("61_table2_lead_genes.py", fx61),
                      ("73_table2_recount.py", fx73)):
        # the fixture's add_argument lines need a stub `ap`
        body = fx.replace("def main():\n    ap_compadj = None\n",
                          "class _AP:\n"
                          "    def add_argument(self, *a, **k):\n"
                          "        print('ARG', a, k.get('default'))\n"
                          "ap = _AP()\n"
                          "def main():\n")
        new, probs = apply_edits(body, FILES[fname])
        assert not probs, (fname, probs)
        assert not verify(new, fname), (fname, verify(new, fname))
        p = os.path.join(d, fname)
        open(p, "w").write(new)

        def run(env=None):
            return subprocess.run([sys.executable, "-I", p],
                                  capture_output=True, text=True, cwd=d,
                                  env=dict(os.environ, **(env or {})))

        r = run()
        assert r.returncode == 0, r.stdout + r.stderr
        if "61_" in fname:
            assert "DISC SE_v2filtered" in r.stdout, r.stdout
            assert "FROZEN_SE 926" in r.stdout, r.stdout
            assert "KNOWN_PCT 79.6" in r.stdout, r.stdout
            r2 = run({"TBL_SUFFIX": "_rl125", "FROZEN_SE": "813",
                      "KNOWN_PCT": "Colorado=80.8,Brigham=81.9,"
                                   "Fuwai=70.8,Stanford=75.8"})
            assert r2.returncode == 0, r2.stdout + r2.stderr
            assert "DISC SE_v2filtered_rl125" in r2.stdout, r2.stdout
            assert "FROZEN_SE 813" in r2.stdout, r2.stdout
            assert "KNOWN_PCT 80.8" in r2.stdout, r2.stdout
            # the two argparse defaults, printed by the stub inside main()
            assert "{ET}_v2compadj_rl125" in r2.stdout, r2.stdout
            assert "Table2_rl125" in r2.stdout, r2.stdout
            assert "{ET}_v2compadj\n" in r.stdout or \
                "{ET}_v2compadj " in r.stdout, r.stdout
            # a bad KNOWN_PCT key must fail loudly
            rb = run({"KNOWN_PCT": "Colorado=1,Brigham=2,Fuwai=3"})
            assert rb.returncode != 0, rb.stdout + rb.stderr
            assert "must name exactly" in (rb.stdout + rb.stderr), rb.stderr
        else:
            assert "TABLE2 Table2.tsv" in r.stdout, r.stdout
            assert "RECS 926 590 18" in r.stdout, r.stdout
            assert "FAIL18 []" in r.stdout, r.stdout
            assert "FAIL17 ['lead retention']" in r.stdout, r.stdout
            r2 = run({"TBL_SUFFIX": "_rl125", "REC_BASE_SIG": "813",
                      "REC_COMPADJ_KEPT": "498", "REC_LEAD_KEPT": "0"})
            assert r2.returncode == 0, r2.stdout + r2.stderr
            assert "TABLE2 Table2_rl125.tsv" in r2.stdout, r2.stdout
            assert "BASE SE_v2filtered_rl125" in r2.stdout, r2.stdout
            assert "RECS 813 498 0" in r2.stdout, r2.stdout
            # REC_LEAD_KEPT=0 -> no expectation, but the flag check still fires
            assert "FAIL17 []" in r2.stdout, r2.stdout
            assert "FLAGMISMATCH ['recomputation says mismatch']" \
                in r2.stdout, r2.stdout

        _n, p3 = apply_edits(new, FILES[fname])
        assert p3, f"{fname}: second application not refused"

    print("[ok] self-test passed: 61 takes TBL_SUFFIX/FROZEN_SE/KNOWN_PCT and "
          "writes Table2<SUF>; a misspelt KNOWN_PCT key fails loudly; 73's four "
          "table names and three recorded counts follow the suffix; "
          "REC_LEAD_KEPT=0 suppresses the expectation while the "
          "flag-vs-recomputation check still fires; unsuffixed behaviour is "
          "unchanged for both; double application refused")
    return 0


if __name__ == "__main__":
    sys.exit(main())
