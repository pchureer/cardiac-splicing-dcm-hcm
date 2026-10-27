#!/usr/bin/env python3
"""
146_patch_62_parta_bothaxes.py -- S5(a): JC on BOTH axes.

THE DEFECT
----------
As built by 144/145, part (a) regresses BOTH cohort tables against the same
discovery set -- `sig`, read from `SE_v2filtered<SUF>`, which is JCEC. The JC
column was therefore a hybrid: JCEC on the x-axis, JC on the y-axis. It
measures whether the REPLICATION-side quantification matters, which is not
what "JC versus JCEC" claims and not what a reader will assume.

It is also not what script 20 does: 20 pairs the two discoveries at its own
line 221, which is part of why its cohort arms move the other way.

THE FIX
-------
The JC discovery already exists at both read lengths, written by script 20:

    jccheck/SE_v2filtered_JC.all_events.tsv
    jccheck_rl125/SE_v2filtered_JC.all_events.tsv

so no new quantification is needed. The JC arm takes it, under the IDENTICAL
significance rule (FDR < 0.05, |dPSI| >= 0.10 on the collapsed table). The
two discovery sets will differ in size; that difference IS part of the
JC-vs-JCEC effect and is deliberately left in the own-events columns. The
matched arm from 145 continues to control for it.

WHAT MOVES AND WHAT MUST NOT
----------------------------
  MOVES:     the JC column, change, JC_m, change_m, n_m
  MUST NOT:  the JCEC column -- SIG["JCEC"] IS `sig`, the same object, so
             the gate against 0.521 / 0.476 / 0.386 / 0.243 is untouched.
             If JCEC moves, revert.

If the JC discovery is absent the patch does NOT fail: the arm falls back to
the JCEC discovery, prints that it has done so, and records a note saying the
panel is replication-side only. A silent fallback would be worse than either.

USAGE
-----
  python 146_patch_62_parta_bothaxes.py [--apply] [--self-test]
"""
import argparse
import difflib
import os
import shutil
import subprocess
import sys
import time

TARGET = "62_build_s5.py"

A_DISC = '    D, err = read_events(dpath, et)\n    if D is None:\n        notes.append(f"(a) skipped: discovery {err}"); print(f"  [fatal for (a)] {err}"); return\n    sig = D[(D.FDR < 0.05) & (D.dpsi.abs() >= 0.10)]\n    print(f"  discovery {dpath}: {len(sig)} significant of {len(D)}")\n'
B_DISC = '    D, err = read_events(dpath, et)\n    if D is None:\n        notes.append(f"(a) skipped: discovery {err}"); print(f"  [fatal for (a)] {err}"); return\n    sig = D[(D.FDR < 0.05) & (D.dpsi.abs() >= 0.10)]\n    print(f"  discovery {dpath}: {len(sig)} significant of {len(D)}")\n\n    # ---- THE JC DISCOVERY, added 21 Oct ------------------------------\n    # Until now BOTH arms regressed against `sig`, the JCEC discovery, so\n    # the JC column was a hybrid: JCEC on the x-axis, JC on the y-axis. It\n    # measured the replication-side quantification only. "JC versus JCEC"\n    # means two complete pipelines, so the JC arm takes the JC discovery\n    # script 20 already writes -- and which 20 itself pairs with the JCEC\n    # discovery at its own line 221.\n    #\n    # The significance rule is applied IDENTICALLY (FDR < 0.05 and\n    # |dPSI| >= 0.10 on the collapsed table), so the two discovery sets\n    # differ only because the quantification differs. That difference is\n    # PART OF the JC-vs-JCEC effect and is deliberately not matched away\n    # here; the matched arm below controls for it separately.\n    jdpath = first_existing(\n        f"jccheck{SUF}/{et}_v2filtered_JC_dedup.all_events.tsv",\n        f"jccheck{SUF}/{et}_v2filtered_JC.all_events.tsv")\n    SIG = {"JCEC": sig}\n    if jdpath is None:\n        print(f"  \\u26a0 NO JC DISCOVERY in jccheck{SUF}/ -- the JC arm falls")\n        print("    back to the JCEC discovery, making this a REPLICATION-SIDE-")\n        print("    ONLY comparison. Either say so in the caption or run")\n        print(f"      TBL_SUFFIX={SUF} sbatch 20_jc_vs_jcec_readlength.sbatch")\n        notes.append("(a) JC arm fell back to the JCEC discovery -- "\n                     "replication-side only, NOT a pipeline comparison")\n        SIG["JC"] = sig\n    else:\n        JD, jerr = read_events(jdpath, et)\n        if JD is None:\n            print(f"  [fatal for (a)] JC discovery: {jerr}")\n            notes.append(f"(a) skipped: JC discovery {jerr}"); return\n        SIG["JC"] = JD[(JD.FDR < 0.05) & (JD.dpsi.abs() >= 0.10)]\n        print(f"  discovery {jdpath}: {len(SIG[\'JC\'])} significant of {len(JD)}")\n        print(f"  note: DISC_SE[\'exon\'] = {DISC_SE[\'exon\']} is the JCEC")\n        print(f"        discovery count; the JC discovery has "\n              f"{len(SIG[\'JC\'])}. The two are not expected to agree and no")\n        print("        gate is applied to the JC one.")\n'
A_MERGE = '            if R is None:\n                r[lab] = np.nan; continue\n            m = sig.merge(R, on="ekey", suffixes=("", "_r"))\n'
B_MERGE = '            if R is None:\n                r[lab] = np.nan; continue\n            # each basis against ITS OWN discovery. SIG["JCEC"] IS `sig`,\n            # so the JCEC arm and the gate below are bit-for-bit unchanged.\n            m = SIG[lab].merge(R, on="ekey", suffixes=("", "_r"))\n'
A_OKHEAD = '    print("  ⚠ `change` IS NOT PURELY A JC-vs-JCEC EFFECT. The own-events pair")\n    print("    is fit on two different event sets -- the JC and JCEC tables")\n    print("    differ in size by 20-30%. `change_m` is the like-for-like figure:")\n    print("    same events, both bases. QUOTE change_m. The own-events pair is")\n    print("    reported beside it to show how large the event-set term is.")\n'
B_OKHEAD = '    print("  ⚠ EACH ARM IS A COMPLETE PIPELINE: JC against the JC discovery,")\n    print("    JCEC against the JCEC discovery. The own-events pair therefore")\n    print("    differs in BOTH the estimator and the event set -- the two")\n    print("    discovery sets and the two cohort tables all differ in size.")\n    print("    `change_m` is the like-for-like figure: same events, both ")\n    print("    bases. QUOTE change_m; the own-events pair beside it shows")\n    print("    how large the event-set term is.")\n'
A_DOCTAIL = "      estimator -- in script 20's Brigham background row that term is 0.09,\n      larger than any JC-vs-JCEC difference in this panel.\n"
B_DOCTAIL = "      estimator -- in script 20's Brigham background row that term is 0.09,\n      larger than any JC-vs-JCEC difference in this panel.\n      Each arm is a COMPLETE pipeline: the JC arm regresses the JC\n      cohort table on the JC discovery (jccheck<SUF>/SE_v2filtered_JC),\n      the JCEC arm on the JCEC discovery. The JCEC arm and the gate are\n      unchanged, so 0.521 / 0.476 / 0.386 / 0.243 still has to hold.\n"

# the post-145 text of the rest of part (a), needed to build the fixture
P145_LOOP = '        M = {}\n        for lab, p in (("JC", jc), ("JCEC", jcec)):\n            R, err = read_events(p, et)\n            if R is None:\n                r[lab] = np.nan; continue\n            m = sig.merge(R, on="ekey", suffixes=("", "_r"))\n            if "se_adj" not in m or "se_adj_r" not in m:\n                r[lab] = np.nan; continue\n            m = m[np.isfinite(m.dpsi.values) & np.isfinite(m.dpsi_r.values)\n                  & np.isfinite(m.se_adj.values) & np.isfinite(m.se_adj_r.values)]\n            M[lab] = m\n            r["n_" + lab] = int(len(m))\n            r[lab] = deming(m.dpsi.values, m.dpsi_r.values,\n                            m.se_adj.values, m.se_adj_r.values)\n        # ---- THE MATCHED ARM, added 21 Oct -----------------------------\n        # Fitting each basis on its own events conflates the estimator\n        # (JC vs JCEC) with the event set (the two tables differ in size by\n        # 20-30%). Script 20 measures the second term directly: Brigham\'s\n        # background slope is 0.662 on its own events and 0.757 on the\n        # JC-intersect-JCEC subset -- a 0.09 gap, larger than any\n        # JC-vs-JCEC difference in this panel. Without these columns a\n        # reader cannot tell which term produced `change`.\n        if len(M) == 2:\n            common = set(M["JC"]["ekey"]) & set(M["JCEC"]["ekey"])\n            for lab in ("JC", "JCEC"):\n                mm = M[lab][M[lab]["ekey"].isin(common)]\n                r[lab + "_m"] = deming(mm.dpsi.values, mm.dpsi_r.values,\n                                       mm.se_adj.values, mm.se_adj_r.values)\n            r["n_m"] = int(len(common))\n        res[coh] = r\n'
P145_PRINT = '    print(f"  {\'cohort\':<10}{\'JCEC\':>9}{\'manuscript\':>12}{\'JC\':>9}{\'change\':>9}"\n          f"{\'|\':>4}{\'JCEC_m\':>9}{\'JC_m\':>9}{\'change_m\':>10}{\'n_m\':>7}")\n    for c in ("Colorado", "Brigham", "Fuwai", "Stanford"):\n        r = res.get(c, {})\n        f = lambda v: f"{v:.3f}" if np.isfinite(v) else "-"\n        g = lambda v: f"{v:+.3f}" if np.isfinite(v) else "-"\n        ch = (r.get("JC", np.nan) - r.get("JCEC", np.nan))\n        chm = (r.get("JC_m", np.nan) - r.get("JCEC_m", np.nan))\n        _n = r.get("n_m")\n        print(f"  {c:<10}{f(r.get(\'JCEC\', np.nan)):>9}{MANUSCRIPT_SLOPE[c]:>12.3f}"\n              f"{f(r.get(\'JC\', np.nan)):>9}{g(ch):>9}"\n              f"{\'|\':>4}{f(r.get(\'JCEC_m\', np.nan)):>9}"\n              f"{f(r.get(\'JC_m\', np.nan)):>9}{g(chm):>10}"\n              f"{(str(_n) if _n else \'-\'):>7}")\n'
P145_OK = '    print("\\n  [ok] the JCEC column reproduces the manuscript\'s slopes")\n    print("  ⚠ `change` IS NOT PURELY A JC-vs-JCEC EFFECT. The own-events pair")\n    print("    is fit on two different event sets -- the JC and JCEC tables")\n    print("    differ in size by 20-30%. `change_m` is the like-for-like figure:")\n    print("    same events, both bases. QUOTE change_m. The own-events pair is")\n    print("    reported beside it to show how large the event-set term is.")\n    print("  ⚠ DO NOT COMPARE THESE ROWS WITH script 20\'s jc_vs_jcec<SUF>.tsv")\n    print("    ROW-FOR-ROW. Script 20 splits the DCM-significant events into")\n    print("    `program` and `background` and reports each separately; this")\n    print("    panel uses the manuscript\'s full significant set, which is why")\n    print("    its Brigham JCEC is 0.476 and 20\'s background JCEC is 0.662.")\n    print("    Script 20\'s MATCHED BACKGROUND arms are the source for the")\n    print("    READ-LENGTH sentence: internal -0.006 vs externals ~+0.04 at")\n    print("    50 bp, but internal -0.021 vs +0.026 / -0.007 / +0.009 at")\n    print("    125 bp, which STRADDLE it -- so \'read length is largely")\n    print("    excluded\' does not survive and needs rewriting.")\n    print("    THIS panel is the source for the JC-vs-JCEC sentence only.")\n'
P145_ROWS = '    for c in sorted(res):\n        rows.append(dict(part="a. JC vs JCEC", row=c,\n                         value_1=round(res[c]["JCEC"], 3), label_1="Deming, JCEC",\n                         value_2=round(res[c]["JC"], 3), label_2="Deming, JC",\n                         change=round(res[c]["JC"] - res[c]["JCEC"], 3)))\n        # the matched pair as a SECOND ROW, not new columns: parts (b) and\n        # (c) write the same seven-key schema into the same table.\n        if (np.isfinite(res[c].get("JCEC_m", np.nan))\n                and np.isfinite(res[c].get("JC_m", np.nan))):\n            rows.append(dict(\n                part="a. JC vs JCEC",\n                row=f"{c}, matched subset (n={res[c].get(\'n_m\', 0)})",\n                value_1=round(res[c]["JCEC_m"], 3),\n                label_1="Deming, JCEC (matched)",\n                value_2=round(res[c]["JC_m"], 3),\n                label_2="Deming, JC (matched)",\n                change=round(res[c]["JC_m"] - res[c]["JCEC_m"], 3)))\n        else:\n            notes.append(f"(a) {c}: NO MATCHED ARM -- one basis unreadable, "\n                         f"so its `change` mixes estimator and event set")\n'
P145_DOC = "  (a) JCEC slopes must match 0.521 / 0.476 / 0.386 / 0.243 within 0.02\n      (JCEC, not JC: the manuscript's tables are JCEC. The JC side comes\n       from script 20's jccheck<SUF>/ output.)\n      Each cohort emits TWO rows: own events, and the JC-intersect-JCEC\n      matched subset. QUOTE THE MATCHED PAIR. The own-events pair is there\n      to show how much of `change` is the event set rather than the\n      estimator -- in script 20's Brigham background row that term is 0.09,\n      larger than any JC-vs-JCEC difference in this panel.\n"

EDITS = [
    ("the discovery block", A_DISC, B_DISC, 1),
    ("the per-basis merge", A_MERGE, B_MERGE, 1),
    ("the advisory's first paragraph", A_OKHEAD, B_OKHEAD, 1),
    ("the docstring's (a) self-check", A_DOCTAIL, B_DOCTAIL, 1),
]


def apply_edits(text):
    if 'SIG = {"JCEC": sig}' in text:
        return text, ["62 already carries 146 — the JC arm has its own "
                      "discovery (restore the .bak to re-apply)"]
    if 'r[lab + "_m"] = deming(' not in text:
        return text, ["145 has not been applied — apply "
                      "145_patch_62_parta_matched.py first; this builds on it"]
    problems, out = [], text
    for name, anchor, repl, want in EDITS:
        got = out.count(anchor)
        if got != want:
            problems.append(
                f"{name}: anchor occurs {got}x, expected {want}x\n"
                f"      first line: {anchor.splitlines()[0][:70]}")
            continue
        out = out.replace(anchor, repl, want)
    return out, problems


def verify(text):
    bad = []
    for w in ('SIG = {"JCEC": sig}',
              '_v2filtered_JC.all_events.tsv"',
              'SIG["JC"] = JD[(JD.FDR < 0.05) & (JD.dpsi.abs() >= 0.10)]',
              'm = SIG[lab].merge(R, on="ekey", suffixes=("", "_r"))',
              'EACH ARM IS A COMPLETE PIPELINE',
              'NO JC DISCOVERY in jccheck',
              'replication-side only, NOT a pipeline comparison'):
        if w not in text:
            bad.append(f"missing: {w[:60]}")
    # the JCEC arm must still run off `sig` itself, not a copy or a re-read
    if 'SIG = {"JCEC": sig}' in text and 'SIG = {"JCEC": sig.copy()' in text:
        bad.append("SIG['JCEC'] is a copy — it must BE `sig` so the gate is "
                   "provably unchanged")
    # Part (a)'s shared-discovery merge must be gone. The check is on the
    # THREE-LINE anchor, not on the bare `sig.merge(...)` line: parts (b)
    # and (c) use that same line against the same discovery and are right
    # to. An earlier version of this check tested the bare string and
    # refused a correctly patched file.
    if A_MERGE in text:
        bad.append("part (a)'s shared-discovery merge survives — the JC arm "
                   "would still use the JCEC discovery")
    if text.count('m = SIG[lab].merge(R, on="ekey", suffixes=("", "_r"))') != 1:
        bad.append("part (a) does not merge against SIG[lab] exactly once")
    # the significance rule must be written identically for both
    if text.count('.FDR < 0.05) & (') < 2:
        bad.append("the JC discovery does not apply the same FDR/|dPSI| rule")
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
    print("146 — S5(a): the JC arm gets the JC discovery")
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
        print("  python3 62_build_s5.py --parts a 2>&1 | tail -45")
        print()
        print("  TBL_SUFFIX=_rl125 \\")
        print("  S5_SLOPE=Colorado=0.5999,Brigham=0.5099,"
              "Stanford=0.4008,Fuwai=0.2948 \\")
        print("  S5_CONC=Colorado=80.8,Brigham=81.9,Fuwai=70.8,Stanford=75.8 \\")
        print("  S5_DISC_SE=exon=813,full=944 \\")
        print("  python3 62_build_s5.py --parts a 2>&1 | tail -45")
        print()
        print("⚠ THE JCEC COLUMN MUST NOT MOVE: 0.476 / 0.243 / 0.386 at")
        print("  50 bp and 0.510 / 0.295 / 0.401 at 125 bp. SIG['JCEC'] is")
        print("  the same object `sig` was, so a shift means something else")
        print("  changed and the patch should be reverted from the .bak.")
        print()
        print("⚠ THE JC COLUMN IS EXPECTED TO MOVE off 0.511 / 0.256 /")
        print("  0.357 (50 bp) and 0.562 / 0.306 / 0.383 (125 bp). Those")
        print("  were the hybrid. If it does NOT move at all, the JC")
        print("  discovery was not found and the fallback fired — check for")
        print("  'NO JC DISCOVERY' in the output.")
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
    print("✓ verified on disk: the JC arm reads its own discovery under the")
    print("  same FDR/|dPSI| rule, the JCEC arm still runs off `sig` itself,")
    print("  and a missing JC discovery degrades loudly instead of silently.")
    return 0


# ------------------------------------------------------------- self-test

_REAL_HELPERS = '\nNKEEP = {"SE": 3, "MXE": 5, "A3SS": 5, "A5SS": 5, "RI": 3}\n\n\ndef read_events(path, et, exon_level=True):\n    if not (os.path.exists(path) and os.path.getsize(path)):\n        return None, "missing"\n    t = pd.read_csv(path, sep="\\t", low_memory=False)\n    if "key" not in t.columns or "FDR" not in t.columns:\n        return None, "no key/FDR column"\n    d = next((c for c in ("dPSI_adj", "dPSI") if c in t.columns), None)\n    if d is None:\n        return None, "no dPSI column"\n    t["ekey"] = ([":".join(str(k).split(":")[:NKEEP[et]]) for k in t["key"]]\n                 if exon_level else t["key"].astype(str))\n    t = t.loc[t.groupby("ekey")["FDR"].idxmin()]\n    return t.rename(columns={d: "dpsi"})[["ekey", "dpsi", "FDR"] +\n                                         (["se_adj"] if "se_adj" in t.columns\n                                          else [])], None\n\n\ndef first_existing(*paths):\n    for p in paths:\n        if os.path.exists(p) and os.path.getsize(p):\n            return p\n    return None\n\n\ndef deming(x, y, sx, sy):\n    ok = np.isfinite(x) & np.isfinite(y) & np.isfinite(sx) & np.isfinite(sy)\n    x, y, sx, sy = x[ok], y[ok], sx[ok], sy[ok]\n    if len(x) < 10:\n        return np.nan\n    lam = np.mean(sy ** 2) / np.mean(sx ** 2) if np.mean(sx ** 2) > 0 else 1.0\n    dx, dy = x - x.mean(), y - y.mean()\n    Sxx, Syy, Sxy = (dx ** 2).sum(), (dy ** 2).sum(), (dx * dy).sum()\n    if Sxy == 0:\n        return np.nan\n    return float((Syy - lam * Sxx + np.sqrt((Syy - lam * Sxx) ** 2\n                                            + 4 * lam * Sxy ** 2)) / (2 * Sxy))\n'


def _fixture():
    """A post-145 62, with every anchor 146 touches in its real context."""
    return (
        '"""\n' + P145_DOC + '"""\n'
        "import os\n"
        "import numpy as np\n"
        "import pandas as pd\n"
        + _REAL_HELPERS +
        'SUF = os.environ.get("TBL_SUFFIX", "")\n'
        'DISC = "{ET}_v2filtered" + SUF\n'
        'DISC_SE = {"exon": 0, "full": 0}\n'
        'MANUSCRIPT_SLOPE = {"Colorado": 0.5, "Brigham": 0.5,\n'
        '                    "Stanford": 0.5, "Fuwai": 0.5}\n'
        "SLOPE_TOL = 99.0\n"
        'REPL = [("Colorado", "{ET}_Colorado_DCMvsNF"),\n'
        '        ("Brigham", "Brigham_{ET}_HCM"),\n'
        '        ("Fuwai", "Fuwai_{ET}_HCM"),\n'
        '        ("Stanford", "Stanford_{ET}_HCM")]\n'
        "\n"
        "def part_a(rows, notes):\n"
        '    et = "SE"\n'
        "    dpath = first_existing("
        'f"{DISC.replace(\'{ET}\', et)}_dedup.all_events.tsv",\n'
        '                           f"{DISC.replace(\'{ET}\', et)}'
        '.all_events.tsv")\n'
        + A_DISC +
        "    out = {}\n"
        "    for coh, pat in REPL:\n"
        '        base = pat.replace("{ET}", et)\n'
        '        jcec = first_existing(f"{base}_dedup.all_events.tsv",\n'
        '                              f"{base}.all_events.tsv")\n'
        "        jc = first_existing(\n"
        '            f"jccheck{SUF}/{coh}_{et}_HCM_JC_dedup.all_events.tsv",\n'
        '            f"jccheck{SUF}/{coh}_{et}_HCM_JC.all_events.tsv")\n'
        '        print(f"  {coh:<10} JCEC -> {jcec or \'*** MISSING ***\'}")\n'
        '        print(f"  {\'\':<10} JC   -> '
        '{jc or \'*** MISSING (script 20 writes these) ***\'}")\n'
        "        out[coh] = (jc, jcec)\n"
        "    have = {c: v for c, v in out.items() if v[0] and v[1]}\n"
        "    res = {}\n"
        "    for coh, (jc, jcec) in have.items():\n"
        "        r = {}\n"
        + P145_LOOP +
        "\n"
        "    print()\n"
        + P145_PRINT
        + P145_OK
        + P145_ROWS +
        "\n"
        "rows, notes = [], []\n"
        "part_a(rows, notes)\n"
        'print("NROWS", len(rows))\n'
        "for _r in rows:\n"
        '    print("ROW|", _r["row"], "|", _r["label_1"], _r["value_1"],\n'
        '          "|", _r["label_2"], _r["value_2"], "|", _r["change"])\n'
        'print("NOTES", notes)\n')


def _rows(stdout):
    own, mat = {}, {}
    for l in stdout.splitlines():
        if not l.startswith("ROW| "):
            continue
        _, row, c1, c2, ch = [s.strip() for s in l.split("|")]
        name = row.split(",")[0]
        (mat if "matched" in row else own)[name] = (
            float(c1.split()[-1]), float(c2.split()[-1]), float(ch))
    return own, mat


def self_test():
    import tempfile
    import numpy as np
    import pandas as pd

    fx = _fixture()
    new, probs = apply_edits(fx)
    assert not probs, probs
    assert not verify(new), verify(new)
    compile(new, "fx", "exec")

    d = tempfile.mkdtemp()
    os.makedirs(os.path.join(d, "jccheck"), exist_ok=True)

    rng = np.random.default_rng(7)
    n = 200
    keys = [f"chr1:{i}:{i + 50}:a:b" for i in range(n)]
    # every discovery event must clear |dPSI| >= 0.10, or the significance
    # filter in the anchor silently empties the fixture
    dx = rng.choice([-1.0, 1.0], n) * rng.uniform(0.12, 0.60, n)

    def tsv(path, idx, vals):
        t = pd.DataFrame({"key": [keys[i] for i in idx], "FDR": 1e-4,
                          "dPSI_adj": vals, "se_adj": 0.05})
        fp = os.path.join(d, path)
        t.to_csv(fp, sep="\t", index=False)

    # JCEC discovery: all 200, dpsi = dx
    tsv("SE_v2filtered.all_events.tsv", range(n), dx)
    # JC discovery: 180 events, dpsi = 0.8 * dx. Regressing the SAME cohort
    # y on this x must give a slope 1/0.8 = 1.25x larger, which is how the
    # test detects that the JC arm switched axes.
    jd = list(range(0, 180))
    tsv("jccheck/SE_v2filtered_JC.all_events.tsv", jd, [0.8 * dx[i] for i in jd])

    jcec_idx = list(range(0, 160))
    jc_idx = list(range(40, 200))
    for c in ("Brigham", "Fuwai", "Stanford"):
        tsv(f"{c}_SE_HCM.all_events.tsv", jcec_idx,
            [0.50 * dx[i] for i in jcec_idx])
        t = pd.DataFrame({"key": [keys[i] for i in jc_idx], "FDR": 1e-4,
                          "dPSI_adj": [0.50 * dx[i] for i in jc_idx],
                          "se_adj": 0.05})
        t.to_csv(os.path.join(d, f"jccheck/{c}_SE_HCM_JC.all_events.tsv"),
                 sep="\t", index=False)
    tsv("SE_Colorado_DCMvsNF.all_events.tsv", jcec_idx,
        [0.50 * dx[i] for i in jcec_idx])

    e = dict(os.environ)
    e.pop("TBL_SUFFIX", None)

    def run(src, name):
        p = os.path.join(d, name)
        open(p, "w").write(src)
        r = subprocess.run([sys.executable, "-I", p], capture_output=True,
                           text=True, cwd=d, env=e)
        assert r.returncode == 0, r.stdout + r.stderr
        return r.stdout

    o0 = run(fx, "fx_pre.py")
    o1 = run(new, "fx_post.py")
    own0, mat0 = _rows(o0)
    own1, mat1 = _rows(o1)

    assert set(own0) == set(own1) == {"Brigham", "Fuwai", "Stanford"}, own1

    for c in own1:
        # the JCEC arm must be BIT-IDENTICAL — this is the promise the dry
        # run makes to CP, so it is tested, not asserted
        assert own0[c][0] == own1[c][0], (c, own0[c][0], own1[c][0])
        assert mat0[c][0] == mat1[c][0], (c, mat0[c][0], mat1[c][0])
        # the JC arm must MOVE, by the 1/0.8 the fixture builds in
        assert own1[c][1] > own0[c][1] + 0.05, (c, own0[c], own1[c])
        assert abs(own1[c][1] / own0[c][1] - 1.25) < 0.02, (c, own1[c])

    # The JC discovery is read, reported, and SIGNIFICANCE-FILTERED ON ITS
    # OWN. The fixture's JC discovery holds 180 events at 0.8x the JCEC
    # dPSI, so at least one falls below the |dPSI| >= 0.10 cut that it
    # cleared in JCEC. A count equal to 180 would mean the filter was
    # inherited rather than re-applied, which is the bug this arm exists to
    # avoid -- so the test requires it to be strictly smaller.
    import re as _re
    _mm = _re.search(
        r"jccheck/SE_v2filtered_JC\.all_events\.tsv: (\d+) significant "
        r"of (\d+)", o1)
    assert _mm, o1
    _nsig, _ntot = int(_mm.group(1)), int(_mm.group(2))
    assert _ntot == 180, (_ntot, o1)
    assert 160 <= _nsig < 180, (_nsig, o1)
    assert "SE_v2filtered.all_events.tsv: 200 significant of 200" in o1, o1
    assert "is the JCEC" in o1 and "discovery count" in o1, o1
    assert "NO JC DISCOVERY" not in o1, o1
    assert "EACH ARM IS A COMPLETE PIPELINE" in o1, o1
    # the superseded wording is gone
    assert "is fit on two different event sets" not in o1, o1

    # --- THE FALLBACK: remove the JC discovery, it must degrade LOUDLY
    os.remove(os.path.join(d, "jccheck/SE_v2filtered_JC.all_events.tsv"))
    o2 = run(new, "fx_post.py")
    assert "NO JC DISCOVERY" in o2, o2
    assert "replication-side only, NOT a pipeline comparison" in o2, o2
    own2, _m2 = _rows(o2)
    # and with the fallback the numbers return to the pre-patch hybrid
    for c in own2:
        assert own2[c][:2] == own0[c][:2], (c, own2[c], own0[c])

    # --- THE REGRESSION THAT BIT US: parts (b) and (c) use the SAME
    #     `m = sig.merge(R, ...)` line against the same discovery, and are
    #     right to. An earlier verify() tested for that bare string and
    #     refused a correctly patched 62. The fixture above has only part
    #     (a), so it could not catch it. This one can.
    fx_bc = fx + (
        "\ndef part_c_lookalike(sig, R, r, lab):\n"
        '    m = sig.merge(R, on="ekey", suffixes=("", "_r"))\n'
        "    return m\n")
    n_bc, p_bc = apply_edits(fx_bc)
    assert not p_bc, p_bc
    v_bc = verify(n_bc)
    assert not v_bc, ("a legitimate sig.merge elsewhere in 62 was treated "
                      "as part (a)'s", v_bc)
    assert 'm = sig.merge(R, on="ekey", suffixes=("", "_r"))' in n_bc, \
        "the other part's merge was rewritten — only part (a) may change"
    compile(n_bc, "fx_bc", "exec")

    # --- refuses without 145
    _n, p4 = apply_edits("print('bare')\n")
    assert p4 and "145 has not been applied" in p4[0], p4

    # --- double application refused without altering the text
    n6, p6 = apply_edits(new)
    assert p6 and "already carries 146" in p6[0], p6
    assert n6 == new, "a refused re-application still altered the text"

    print("[ok] self-test passed: the JC arm reads jccheck<SUF>/"
          "SE_v2filtered_JC under the same FDR/|dPSI| rule and its slope "
          "moves by exactly the 1.25x the fixture builds in, while THE JCEC "
          "ARM AND ITS MATCHED SLOPE ARE BIT-IDENTICAL to the pre-patch run; "
          "both discovery counts are printed and differ (180 vs 200); "
          "removing the JC discovery degrades LOUDLY, records the "
          "replication-side-only note, and returns exactly the pre-patch "
          "hybrid numbers; the superseded advisory wording is gone; refuses "
          "without 145; double application refused; and a legitimate "
          "sig.merge() elsewhere in 62 is left alone instead of "
          "refusing the patch")
    return 0


if __name__ == "__main__":
    sys.exit(main())
