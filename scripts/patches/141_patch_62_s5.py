#!/usr/bin/env python3
"""
141_patch_62_s5.py -- Supplementary Table S5 at 125 bp.

Three of 62's self-checks ABORT their own part rather than warn, and all
three hold 50 bp values:

  (a) JC slopes must match 0.521 / 0.476 / 0.386 / 0.243 within 0.02
  (b) the skipped-exon pair must be 172 (50 bp) and 167 (125 bp)
  (c) exon-level concordance must match 79.6 / 81.2 / 69.2 / 76.1 within 0.5

(a) and (c) are re-based here. **(b) is NOT, and deliberately so.**

WHY (b) IS LEFT ALONE
---------------------
S5(b) is already a 50-versus-125 comparison — it is the read-length panel. Its
`PROBE_SE = (172, 167)` is a PROBE result, from `probe_rl50/` and
`probe_rl125/`, not the discovery cohort. The same script knows the
difference: `DISC_SE = {"exon": 926, "full": 1074}` is the real discovery set,
roughly five times larger.

So there is nothing in (b) to re-base, and re-basing it would be wrong. What
(b) does need is a SCOPE statement, because its headline −2.9% and the
promotion's cohort-level change are different populations and will be read as
the same number:

    probe, SE          172 -> 167    -2.9%     (S5(b))
    discovery, SE      926 -> 813   -12.2%     exon-level, collapse-first
    discovery, SE full 1074 -> 944  -12.1%     (62's DISC_SE; 21 agrees exactly)

**They do not conflict — they measure different sets.** The risk is entirely
in the text: if the manuscript cites S5(b) to say read length barely changes
the event count, it generalises a probe to the cohort, and the full
re-quantification contradicts that generalisation by a factor of four. This
patch makes (b) print its own scope so the table cannot be read that way.

A useful convergence found while checking: `DISC_SE["full"] = 1074` is exactly
what script 21 reports as MAGNet-significant SE at 50 bp, and 21 gives 944 at
125 bp — two independent scripts agreeing on both ends.

WHAT MOVES
----------
  line  56  DISC = "{ET}_v2filtered"            -> suffixed
  line  84  MANUSCRIPT_SLOPE                    -> S5_SLOPE
  line  87  MANUSCRIPT_CONC                     -> S5_CONC
  line  89  DISC_SE {"exon": 926, "full": 1074} -> S5_DISC_SE
  line 423  --out default "S5"                  -> S5<SUF>
  part (b)                                      -> scope note, no constants

`REPL` stays unsuffixed — no replication cohort was re-quantified — and
`PROBE_SE`, `RECORDED_PCT` and `RECORDED_A3SS_BASE` stay exactly as they are.

MEASURED 125 bp VALUES
----------------------
  S5_SLOPE   Colorado 0.5999  Brigham 0.5099  Stanford 0.4008  Fuwai 0.2948
             (45, agreeing with 33's independent deming() to 4 dp)
  S5_CONC    Colorado 80.8    Brigham 81.9    Fuwai 70.8       Stanford 75.8
             (33 section D; recounted by 45 and 61)
  S5_DISC_SE exon 813         full 944        (813 from the discovery
             re-quantification; 944 from script 21)

`DISC` resolves to `SE_v2filtered_rl125.all_events.tsv` — checked: no
`_dedup` variant exists at either read length, so `first_existing` falls
through to the plain name at both.

USAGE
-----
  python 141_patch_62_s5.py [--apply] [--self-test]
"""
import argparse
import difflib
import os
import shutil
import sys
import time

TARGET = "62_build_s5.py"

RUN125 = (
    "TBL_SUFFIX=_rl125 \\\n"
    "S5_SLOPE=Colorado=0.5999,Brigham=0.5099,"
    "Stanford=0.4008,Fuwai=0.2948 \\\n"
    "S5_CONC=Colorado=80.8,Brigham=81.9,Fuwai=70.8,Stanford=75.8 \\\n"
    "S5_DISC_SE=exon=813,full=944 \\\n"
    "python3 62_build_s5.py")

# ------------------------------------------------------------------- edit 1
A_DISC = 'DISC = "{ET}_v2filtered"\n'
B_DISC = '''# Read length. DISC is a STEM: the readers below append
# "_dedup.all_events.tsv" or ".all_events.tsv" to it, and no _dedup variant
# exists at either read length, so this resolves to
# SE_v2filtered<SUF>.all_events.tsv. REPL below stays unsuffixed — no
# replication cohort was re-quantified.
SUF = os.environ.get("TBL_SUFFIX", "")
DISC = "{ET}_v2filtered" + SUF
'''

# ------------------------------------------------------------------- edit 2
A_REC = '''# the manuscript's own numbers, which these must reproduce
MANUSCRIPT_SLOPE = {"Colorado": 0.521, "Brigham": 0.476,
                    "Stanford": 0.386, "Fuwai": 0.243}
SLOPE_TOL = 0.02
MANUSCRIPT_CONC = {"Colorado": 79.6, "Brigham": 81.2, "Fuwai": 69.2, "Stanford": 76.1}
# the two conventions must give DIFFERENT discovery counts, or the switch did nothing
DISC_SE = {"exon": 926, "full": 1074}
'''
B_REC = A_REC + '''
# --------------------------------------------------------------------------
# 125 bp record, placed AFTER every literal it rewrites.
#
# Parts (a) and (c) ABORT rather than warn, so a 50 bp record against the
# 125 bp tables does not produce a wrong supplement — it refuses a correct
# one, reporting "JC slopes do not reproduce the manuscript", which reads as
# a pipeline fault. All three are therefore REQUIRED when TBL_SUFFIX is set.
#
# PROBE_SE, RECORDED_PCT and RECORDED_A3SS_BASE are deliberately NOT here:
# part (b) is already a 50-vs-125 comparison on the probe directories and has
# nothing to re-base.
_S5_REQ = ("S5_SLOPE", "S5_CONC", "S5_DISC_SE")
_s5env = {_k: os.environ.get(_k, "").strip() for _k in _S5_REQ}
if SUF:
    _miss = [_k for _k in _S5_REQ if not _s5env[_k]]
    if _miss:
        raise SystemExit(
            "REFUSING — TBL_SUFFIX=%s but these are unset: %s\\n\\n"
            "  Parts (a) and (c) abort rather than warn, so a 50 bp record\\n"
            "  against the 125 bp tables refuses a CORRECT supplement.\\n\\n"
            "    %s\\n" % (SUF, ", ".join(_miss), @@RUN125@@))

if _s5env["S5_SLOPE"]:
    _d = dict(_p.split("=", 1) for _p in _s5env["S5_SLOPE"].split(","))
    assert set(_d) == set(MANUSCRIPT_SLOPE), (
        "S5_SLOPE must name exactly %s, got %s"
        % (sorted(MANUSCRIPT_SLOPE), sorted(_d)))
    MANUSCRIPT_SLOPE = {_k: float(_v) for _k, _v in _d.items()}

if _s5env["S5_CONC"]:
    _d = dict(_p.split("=", 1) for _p in _s5env["S5_CONC"].split(","))
    assert set(_d) == set(MANUSCRIPT_CONC), (
        "S5_CONC must name exactly %s, got %s"
        % (sorted(MANUSCRIPT_CONC), sorted(_d)))
    MANUSCRIPT_CONC = {_k: float(_v) for _k, _v in _d.items()}

if _s5env["S5_DISC_SE"]:
    _d = dict(_p.split("=", 1) for _p in _s5env["S5_DISC_SE"].split(","))
    assert set(_d) == set(DISC_SE), (
        "S5_DISC_SE must name exactly %s, got %s"
        % (sorted(DISC_SE), sorted(_d)))
    DISC_SE = {_k: int(_v) for _k, _v in _d.items()}
# --------------------------------------------------------------------------
'''

# ------------------------------------------------------------------- edit 3
A_OUT = '    ap.add_argument("--out", default="S5")\n'
B_OUT = '    ap.add_argument("--out", default="S5" + SUF)\n'

# ------------------------------------------------------------------- edit 4
A_PROBE = ('    print("\\n  [ok] skipped exons reproduce the recorded '
           '172 / 167")\n')
B_PROBE = ('    print("\\n  [ok] skipped exons reproduce the recorded '
           '172 / 167")\n'
           '    # SCOPE. Added 21 Oct. (b) is a PROBE result and will be read\n'
           '    # as a cohort result unless it says otherwise.\n'
           '    print("\\n  ⚠ SCOPE: these are PROBE counts, from '
           'probe_rl50/ and")\n'
           '    print("    probe_rl125/ — roughly a fifth of the discovery '
           'set. The")\n'
           '    print("    cohort-level change is FOUR TIMES larger:")\n'
           '    print(f"      probe, SE            172 -> 167     -2.9%")\n'
           '    print(f"      discovery, SE exon   {DISC_SE[\'exon\']:>4} -> '
           '???    see the 125 bp record")\n'
           '    print(f"      discovery, SE full   {DISC_SE[\'full\']:>4} -> '
           '???    (script 21 agrees on both)")\n'
           '    print("    They do not conflict — they are different event '
           'sets. But a")\n'
           '    print("    sentence citing (b) for how much read length '
           'matters is")\n'
           '    print("    generalising a probe to the cohort, and the full")\n'
           '    print("    re-quantification contradicts that by a factor of '
           'four.")\n')

EDITS = [
    ("SUF and the suffixed discovery stem", A_DISC, B_DISC, 1),
    ("the 125 bp record block", A_REC, B_REC, 1),
    ("the output stem, suffixed", A_OUT, B_OUT, 1),
    ("part (b)'s scope note", A_PROBE, B_PROBE, 1),
]

B_REC = B_REC.replace("@@RUN125@@", repr(RUN125.replace("\n", "\n    ")))
assert "@@RUN125@@" not in B_REC and "TBL_SUFFIX=_rl125" in B_REC
EDITS[1] = ("the 125 bp record block", A_REC, B_REC, 1)


def apply_edits(text):
    if 'SUF = os.environ.get("TBL_SUFFIX", "")' in text:
        return text, ["SUF is already defined — 62 is already patched "
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
              'DISC = "{ET}_v2filtered" + SUF',
              '_S5_REQ = ("S5_SLOPE", "S5_CONC", "S5_DISC_SE")',
              "MANUSCRIPT_SLOPE = {_k: float(_v)",
              "MANUSCRIPT_CONC = {_k: float(_v)",
              "DISC_SE = {_k: int(_v)",
              'ap.add_argument("--out", default="S5" + SUF)',
              "⚠ SCOPE: these are PROBE counts",
              "generalising a probe to the cohort"):
        if w not in text:
            bad.append(f"missing: {w[:58]}")
    for gone in ('DISC = "{ET}_v2filtered"\n',
                 'ap.add_argument("--out", default="S5")'):
        if gone in text:
            bad.append(f"unsuffixed site survives: {gone[:50]}")
    # (b)'s own constants must be UNTOUCHED
    for keep in ("PROBE_SE = (172, 167)",
                 'RECORDED_PCT = {"SE": -2.9',
                 "RECORDED_A3SS_BASE = 36"):
        if keep not in text:
            bad.append(f"part (b)'s record was altered: {keep[:44]}")
    if "S5_PROBE" in text:
        bad.append("part (b) was given an override it must not have")
    # REPL must stay unsuffixed
    for y in ('("Colorado", "{ET}_Colorado_DCMvsNF")',
              '("Brigham", "Brigham_{ET}_HCM")'):
        if y not in text:
            bad.append(f"a replication stem was altered: {y[:44]}")
    if "_Colorado_DCMvsNF\" + SUF" in text:
        bad.append("a replication cohort was wrongly suffixed")
    # SUF module-level, above DISC and above main()
    i_s = text.find('SUF = os.environ.get("TBL_SUFFIX", "")')
    i_d = text.find('DISC = "{ET}_v2filtered" + SUF')
    i_m = text.find("def main(")
    if i_s < 0 or i_d < 0 or i_s > i_d:
        bad.append("SUF is not defined above DISC")
    if i_m < 0 or i_s > i_m:
        bad.append("SUF is not above def main()")
    if "\n    SUF = os.environ.get" in text:
        bad.append("SUF is indented — it would be a local, not a global")
    # the record block must follow every literal it rewrites
    i_o = text.find("_s5env = {_k: os.environ.get(_k")
    for lit, nm in (("MANUSCRIPT_SLOPE = {\"Colorado\": 0.521", "SLOPE"),
                    ("MANUSCRIPT_CONC = {\"Colorado\": 79.6", "CONC"),
                    ('DISC_SE = {"exon": 926, "full": 1074}', "DISC_SE")):
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
    print("141 — 62_build_s5.py: parts (a) and (c) re-based, (b) scoped")
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
        print("  # 1. CALIBRATION — all three parts must still build at 50 bp")
        print("  cp S5.tsv S5.PRESERVED 2>/dev/null")
        print("  python3 62_build_s5.py 2>&1 | tail -40")
        print("  cmp S5.tsv S5.PRESERVED && echo '50bp BYTE-IDENTICAL'")
        print()
        print("  # 2. the refusal must name all three groups")
        print("  TBL_SUFFIX=_rl125 python3 62_build_s5.py 2>&1 | head -10")
        print()
        print("  # 3. 125 bp")
        print("  " + RUN125.replace("\n", "\n  "))
        print()
        print("⚠ PART (b) WILL PRODUCE THE SAME NUMBERS AT BOTH READ LENGTHS,")
        print("  and that is correct — it reads the probe directories, which")
        print("  do not move. Only its scope note is new. If (b) ever DID")
        print("  change with TBL_SUFFIX, something would be wrong.")
        print()
        print("⚠ WATCH (a)'s JC/JCEC COMPARISON. Its claim is that the")
        print("  internal control barely moves while the external cohorts")
        print("  shift by ~0.05. 22 has now shown the same control is flat")
        print("  across 50/125 bp too (1.072 -> 1.071). If (a) reproduces at")
        print("  125 bp, the two results agree by the same logic on two axes.")
        print()
        print("⬜ The Methods sentence that cites S5(b) for how much read")
        print("  length matters needs rewriting, not renumbering — the probe")
        print("  says -2.9% and the cohort says -12.2%.")
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
    print("✓ verified on disk: (a) and (c)'s records are overridable and")
    print("  required when suffixed; the discovery stem and output follow")
    print("  SUF; REPL and all of part (b)'s constants are untouched; (b)")
    print("  now prints its own scope.")
    return 0


def self_test():
    import subprocess
    import tempfile

    fx = ("import os\n"
          'REPL = [("Colorado", "{ET}_Colorado_DCMvsNF"),\n'
          '        ("Brigham", "Brigham_{ET}_HCM")]\n'
          + A_DISC
          + A_REC
          + "CONC_TOL = 0.5\n"
            "PROBE_SE = (172, 167)\n"
            'RECORDED_PCT = {"SE": -2.9, "MXE": 0.0}\n'
            "RECORDED_A3SS_BASE = 36\n"
            "def part_b():\n"
          + A_PROBE
          + "def main():\n"
            "    class _AP:\n"
            "        def add_argument(self, *a, **k):\n"
            "            print('ARG', a[0], k.get('default'))\n"
            "    ap = _AP()\n"
          + A_OUT
          + '    print("DISC", DISC.replace("{ET}", "SE"))\n'
            '    print("SLOPE", MANUSCRIPT_SLOPE["Colorado"], '
            'MANUSCRIPT_SLOPE["Fuwai"])\n'
            '    print("CONC", MANUSCRIPT_CONC["Colorado"], '
            'MANUSCRIPT_CONC["Stanford"])\n'
            '    print("DISCSE", DISC_SE["exon"], DISC_SE["full"])\n'
            '    print("PROBE", PROBE_SE)\n'
            "    part_b()\n"
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
        for k in ("TBL_SUFFIX", "S5_SLOPE", "S5_CONC", "S5_DISC_SE"):
            e.pop(k, None)
        e.update(env or {})
        return subprocess.run([sys.executable, "-I", p], capture_output=True,
                              text=True, cwd=d, env=e)

    # 1. unsuffixed -> every original value
    r = run()
    assert r.returncode == 0, r.stdout + r.stderr
    assert "ARG --out S5\n" in r.stdout, r.stdout
    assert "DISC SE_v2filtered" in r.stdout, r.stdout
    assert "DISC SE_v2filtered_" not in r.stdout, r.stdout
    assert "SLOPE 0.521 0.243" in r.stdout, r.stdout
    assert "CONC 79.6 76.1" in r.stdout, r.stdout
    assert "DISCSE 926 1074" in r.stdout, r.stdout
    assert "PROBE (172, 167)" in r.stdout, r.stdout
    # the scope note appears at BOTH read lengths
    assert "SCOPE: these are PROBE counts" in r.stdout, r.stdout
    assert "172 -> 167     -2.9%" in r.stdout, r.stdout
    assert "926 -> ???" in r.stdout, r.stdout

    # 2. TBL_SUFFIX alone refuses, naming all three
    r2 = run({"TBL_SUFFIX": "_rl125"})
    assert r2.returncode != 0, r2.stdout + r2.stderr
    b2 = r2.stdout + r2.stderr
    for k in ("S5_SLOPE", "S5_CONC", "S5_DISC_SE"):
        assert k in b2, (k, b2)
    assert "refuses a CORRECT supplement" in b2, b2
    assert "exon=813,full=944" in b2, b2

    full = {"TBL_SUFFIX": "_rl125",
            "S5_SLOPE": "Colorado=0.5999,Brigham=0.5099,"
                        "Stanford=0.4008,Fuwai=0.2948",
            "S5_CONC": "Colorado=80.8,Brigham=81.9,Fuwai=70.8,Stanford=75.8",
            "S5_DISC_SE": "exon=813,full=944"}

    # 3. the full record moves (a) and (c) but NOT (b)
    r3 = run(full)
    assert r3.returncode == 0, r3.stdout + r3.stderr
    assert "ARG --out S5_rl125" in r3.stdout, r3.stdout
    assert "DISC SE_v2filtered_rl125" in r3.stdout, r3.stdout
    assert "SLOPE 0.5999 0.2948" in r3.stdout, r3.stdout
    assert "CONC 80.8 75.8" in r3.stdout, r3.stdout
    assert "DISCSE 813 944" in r3.stdout, r3.stdout
    assert "PROBE (172, 167)" in r3.stdout, r3.stdout      # UNCHANGED
    assert "172 -> 167     -2.9%" in r3.stdout, r3.stdout  # UNCHANGED
    assert "813 -> ???" in r3.stdout, r3.stdout            # the note follows

    # 4. a wrong key set in any group names the expected keys
    for var, bad_val, want in (
            ("S5_SLOPE", "Colorado=0.6,Brigham=0.5", "Stanford"),
            ("S5_CONC", "Colorado=80.8,Brigham=81.9,Fuwai=70.8", "Stanford"),
            ("S5_DISC_SE", "exon=813", "full")):
        rb = run(dict(full, **{var: bad_val}))
        assert rb.returncode != 0, (var, rb.stdout)
        bb = rb.stdout + rb.stderr
        assert "must name exactly" in bb and want in bb, (var, bb)

    # 5. the replication stems never move
    for out in (r.stdout, r3.stdout):
        assert "Colorado_DCMvsNF_rl125" not in out, out
    assert '("Brigham", "Brigham_{ET}_HCM")' in new

    # 6. double application refused, text unchanged
    n6, p6 = apply_edits(new)
    assert p6, "second application was not refused"
    assert n6 == new, "a refused re-application still altered the text"

    print("[ok] self-test passed: unsuffixed values and the output stem are "
          "unchanged; TBL_SUFFIX alone refuses naming all three groups; the "
          "full record moves the slopes, concordances, discovery counts, the "
          "discovery stem and the output name, while PROBE_SE and part (b)'s "
          "counts stay EXACTLY as they are at both read lengths; (b) prints "
          "its probe-versus-cohort scope either way; a wrong key set names "
          "the expected keys; the replication stems never move; double "
          "application refused")
    return 0


if __name__ == "__main__":
    sys.exit(main())
