#!/usr/bin/env python3
"""
145_patch_62_parta_matched.py -- give S5(a) a matched arm.

WHY
---
As built by 144, part (a) fits each basis on its own event set: the JC table
and the JCEC table are merged separately onto the discovery-significant set,
and the two merges have different sizes (20-30% apart in script 20's own
counts). The reported `change` therefore mixes two effects:

    change = (JC estimator - JCEC estimator) + (JC event set - JCEC event set)

and script 20 measures the second term directly. Brigham's background slope
is 0.662489 on its own events and 0.756999 on the JC-intersect-JCEC subset --
a 0.09 gap, LARGER than any JC-vs-JCEC difference anywhere in this panel.
(Those two cells are also the long-standing 0.662 / 0.757 puzzle, now closed.)

A supplement panel that cannot distinguish the two is not interpretable.

WHAT IT ADDS
------------
For each cohort, on the JC-intersect-JCEC subset of the discovery-significant
events, a second Deming fit per basis:

    JCEC_m, JC_m, change_m, n_m

`change_m` is the like-for-like JC-vs-JCEC figure. The own-events pair stays,
because the difference between the two IS the event-set cost and is worth
showing.

The printed table gains four columns. S5.tsv gains a SECOND ROW per cohort
("<cohort>, matched subset (n=...)") rather than new columns, so the
part/row/value_1/label_1/value_2/label_2/change schema that parts (b) and (c)
also write is unchanged.

IT ALSO CORRECTS TWO FALSE LINES 144 LEFT BEHIND
------------------------------------------------
1. "the JC column is therefore computed on the same event set and exclusions"
   -- it is not, which is the whole reason for this patch.
2. The advisory printed script 20's matched BACKGROUND arms directly beneath
   this panel's unmatched ALL-EVENT table, inviting a row-for-row comparison
   that does not hold: 20 splits DCM-significant events into `program` and
   `background` and reports each separately, where this panel uses the
   manuscript's full significant set. The replacement says which table is the
   source for which sentence.

USAGE
-----
  python 145_patch_62_parta_matched.py [--apply] [--self-test]
"""
import argparse
import difflib
import os
import shutil
import subprocess
import sys
import time

TARGET = "62_build_s5.py"

# ---------------------------------------------------------------- anchors
# All five are the text as 144 left it. 144 must be applied first.

A_LOOP = '''        for lab, p in (("JC", jc), ("JCEC", jcec)):
            R, err = read_events(p, et)
            if R is None:
                r[lab] = np.nan; continue
            m = sig.merge(R, on="ekey", suffixes=("", "_r"))
            if "se_adj" not in m or "se_adj_r" not in m:
                r[lab] = np.nan; continue
            r[lab] = deming(m.dpsi.values, m.dpsi_r.values,
                            m.se_adj.values, m.se_adj_r.values)
        res[coh] = r
'''

B_LOOP = '''        M = {}
        for lab, p in (("JC", jc), ("JCEC", jcec)):
            R, err = read_events(p, et)
            if R is None:
                r[lab] = np.nan; continue
            m = sig.merge(R, on="ekey", suffixes=("", "_r"))
            if "se_adj" not in m or "se_adj_r" not in m:
                r[lab] = np.nan; continue
            m = m[np.isfinite(m.dpsi.values) & np.isfinite(m.dpsi_r.values)
                  & np.isfinite(m.se_adj.values) & np.isfinite(m.se_adj_r.values)]
            M[lab] = m
            r["n_" + lab] = int(len(m))
            r[lab] = deming(m.dpsi.values, m.dpsi_r.values,
                            m.se_adj.values, m.se_adj_r.values)
        # ---- THE MATCHED ARM, added 21 Oct -----------------------------
        # Fitting each basis on its own events conflates the estimator
        # (JC vs JCEC) with the event set (the two tables differ in size by
        # 20-30%). Script 20 measures the second term directly: Brigham's
        # background slope is 0.662 on its own events and 0.757 on the
        # JC-intersect-JCEC subset -- a 0.09 gap, larger than any
        # JC-vs-JCEC difference in this panel. Without these columns a
        # reader cannot tell which term produced `change`.
        if len(M) == 2:
            common = set(M["JC"]["ekey"]) & set(M["JCEC"]["ekey"])
            for lab in ("JC", "JCEC"):
                mm = M[lab][M[lab]["ekey"].isin(common)]
                r[lab + "_m"] = deming(mm.dpsi.values, mm.dpsi_r.values,
                                       mm.se_adj.values, mm.se_adj_r.values)
            r["n_m"] = int(len(common))
        res[coh] = r
'''

A_PRINT = '''    print(f"  {'cohort':<10}{'JCEC':>9}{'manuscript':>12}{'JC':>9}{'change':>9}")
    for c in ("Colorado", "Brigham", "Fuwai", "Stanford"):
        r = res.get(c, {})
        f = lambda v: f"{v:.3f}" if np.isfinite(v) else "-"
        ch = (r.get("JC", np.nan) - r.get("JCEC", np.nan))
        print(f"  {c:<10}{f(r.get('JCEC', np.nan)):>9}{MANUSCRIPT_SLOPE[c]:>12.3f}"
              f"{f(r.get('JC', np.nan)):>9}{(f'{ch:+.3f}' if np.isfinite(ch) else '-'):>9}")
'''

B_PRINT = '''    print(f"  {'cohort':<10}{'JCEC':>9}{'manuscript':>12}{'JC':>9}{'change':>9}"
          f"{'|':>4}{'JCEC_m':>9}{'JC_m':>9}{'change_m':>10}{'n_m':>7}")
    for c in ("Colorado", "Brigham", "Fuwai", "Stanford"):
        r = res.get(c, {})
        f = lambda v: f"{v:.3f}" if np.isfinite(v) else "-"
        g = lambda v: f"{v:+.3f}" if np.isfinite(v) else "-"
        ch = (r.get("JC", np.nan) - r.get("JCEC", np.nan))
        chm = (r.get("JC_m", np.nan) - r.get("JCEC_m", np.nan))
        _n = r.get("n_m")
        print(f"  {c:<10}{f(r.get('JCEC', np.nan)):>9}{MANUSCRIPT_SLOPE[c]:>12.3f}"
              f"{f(r.get('JC', np.nan)):>9}{g(ch):>9}"
              f"{'|':>4}{f(r.get('JCEC_m', np.nan)):>9}"
              f"{f(r.get('JC_m', np.nan)):>9}{g(chm):>10}"
              f"{(str(_n) if _n else '-'):>7}")
'''

A_OK = (
    '    print("\\n  [ok] the JCEC column reproduces the manuscript\'s slopes;'
    ' the JC")\n'
    '    print("       column is therefore computed on the same event set and'
    ' exclusions")\n'
    '    print("  ⚠ READ THE SIGN AND SIZE OF `change` AGAINST THE'
    ' MAGNet-INTERNAL ARM")\n'
    '    print("    in script 20\'s own jc_vs_jcec<SUF>.tsv before writing a'
    ' sentence.")\n'
    '    print("    At 50 bp the internal arm moved -0.006 while the cohorts'
    ' moved")\n'
    '    print("    ~+0.04, which is what \'read length is largely excluded\''
    ' rested on.")\n'
    '    print("    At 125 bp the internal arm moves -0.021 and the cohorts'
    ' +0.026 /")\n'
    '    print("    -0.007 / +0.009 -- the premise does not hold and the'
    ' sentence")\n'
    '    print("    needs rewriting, not re-basing.")\n'
)

B_OK = (
    '    print("\\n  [ok] the JCEC column reproduces the manuscript\'s'
    ' slopes")\n'
    '    print("  ⚠ `change` IS NOT PURELY A JC-vs-JCEC EFFECT. The'
    ' own-events pair")\n'
    '    print("    is fit on two different event sets -- the JC and JCEC'
    ' tables")\n'
    '    print("    differ in size by 20-30%. `change_m` is the like-for-like'
    ' figure:")\n'
    '    print("    same events, both bases. QUOTE change_m. The own-events'
    ' pair is")\n'
    '    print("    reported beside it to show how large the event-set term'
    ' is.")\n'
    '    print("  ⚠ DO NOT COMPARE THESE ROWS WITH script 20\'s'
    ' jc_vs_jcec<SUF>.tsv")\n'
    '    print("    ROW-FOR-ROW. Script 20 splits the DCM-significant events'
    ' into")\n'
    '    print("    `program` and `background` and reports each separately;'
    ' this")\n'
    '    print("    panel uses the manuscript\'s full significant set, which is'
    ' why")\n'
    '    print("    its Brigham JCEC is 0.476 and 20\'s background JCEC is'
    ' 0.662.")\n'
    '    print("    Script 20\'s MATCHED BACKGROUND arms are the source for the'
    '")\n'
    '    print("    READ-LENGTH sentence: internal -0.006 vs externals ~+0.04'
    ' at")\n'
    '    print("    50 bp, but internal -0.021 vs +0.026 / -0.007 / +0.009 at'
    '")\n'
    '    print("    125 bp, which STRADDLE it -- so \'read length is largely'
    '")\n'
    '    print("    excluded\' does not survive and needs rewriting.")\n'
    '    print("    THIS panel is the source for the JC-vs-JCEC sentence'
    ' only.")\n'
)

A_ROWS = '''    for c in sorted(res):
        rows.append(dict(part="a. JC vs JCEC", row=c,
                         value_1=round(res[c]["JCEC"], 3), label_1="Deming, JCEC",
                         value_2=round(res[c]["JC"], 3), label_2="Deming, JC",
                         change=round(res[c]["JC"] - res[c]["JCEC"], 3)))
'''

B_ROWS = '''    for c in sorted(res):
        rows.append(dict(part="a. JC vs JCEC", row=c,
                         value_1=round(res[c]["JCEC"], 3), label_1="Deming, JCEC",
                         value_2=round(res[c]["JC"], 3), label_2="Deming, JC",
                         change=round(res[c]["JC"] - res[c]["JCEC"], 3)))
        # the matched pair as a SECOND ROW, not new columns: parts (b) and
        # (c) write the same seven-key schema into the same table.
        if (np.isfinite(res[c].get("JCEC_m", np.nan))
                and np.isfinite(res[c].get("JC_m", np.nan))):
            rows.append(dict(
                part="a. JC vs JCEC",
                row=f"{c}, matched subset (n={res[c].get('n_m', 0)})",
                value_1=round(res[c]["JCEC_m"], 3),
                label_1="Deming, JCEC (matched)",
                value_2=round(res[c]["JC_m"], 3),
                label_2="Deming, JC (matched)",
                change=round(res[c]["JC_m"] - res[c]["JCEC_m"], 3)))
        else:
            notes.append(f"(a) {c}: NO MATCHED ARM -- one basis unreadable, "
                         f"so its `change` mixes estimator and event set")
'''

A_DOC = '''  (a) JCEC slopes must match 0.521 / 0.476 / 0.386 / 0.243 within 0.02
      (JCEC, not JC: the manuscript's tables are JCEC. The JC side comes
       from script 20's jccheck<SUF>/ output.)
'''

B_DOC = '''  (a) JCEC slopes must match 0.521 / 0.476 / 0.386 / 0.243 within 0.02
      (JCEC, not JC: the manuscript's tables are JCEC. The JC side comes
       from script 20's jccheck<SUF>/ output.)
      Each cohort emits TWO rows: own events, and the JC-intersect-JCEC
      matched subset. QUOTE THE MATCHED PAIR. The own-events pair is there
      to show how much of `change` is the event set rather than the
      estimator -- in script 20's Brigham background row that term is 0.09,
      larger than any JC-vs-JCEC difference in this panel.
'''

EDITS = [
    ("part (a)'s per-basis slope loop", A_LOOP, B_LOOP, 1),
    ("the printed table", A_PRINT, B_PRINT, 1),
    ("the post-gate advisory", A_OK, B_OK, 1),
    ("the S5 row emit", A_ROWS, B_ROWS, 1),
    ("the docstring's (a) self-check", A_DOC, B_DOC, 1),
]


def apply_edits(text):
    if 'r[lab + "_m"] = deming(' in text:
        return text, ["62 already carries 145 — the matched arm is present "
                      "(restore the .bak to re-apply)"]
    if "'*** MISSING (script 20 writes these) ***'" not in text:
        return text, ["144 has not been applied — apply "
                      "144_patch_62_parta_labels.py first; this builds on it"]
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
    for w in ('M = {}',
              'common = set(M["JC"]["ekey"]) & set(M["JCEC"]["ekey"])',
              'r[lab + "_m"] = deming(',
              'r["n_m"] = int(len(common))',
              "{'JCEC_m':>9}",
              "{'change_m':>10}",
              'matched subset (n=',
              'Deming, JCEC (matched)',
              'QUOTE change_m',
              'DO NOT COMPARE THESE ROWS'):
        if w not in text:
            bad.append(f"missing: {w[:58]}")
    # the two false lines 144 left must be GONE
    if "computed on the same event set and exclusions" in text:
        bad.append("the false 'same event set and exclusions' line survives")
    if "READ THE SIGN AND SIZE OF `change` AGAINST" in text:
        bad.append("the old advisory survives — it invites the row-for-row "
                   "comparison with script 20 that does not hold")
    # the matched fit must not reuse the own-events frame
    if 'mm = M[lab][M[lab]["ekey"].isin(common)]' not in text:
        bad.append("the matched fit does not restrict to the common ekeys")
    # the schema must not have grown new keys
    for k in ("value_3", "label_3", "n_matched="):
        if k in text:
            bad.append(f"the row schema grew a new key ({k}) — parts (b) and "
                       f"(c) write the same table")
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
    print("145 — S5(a): the matched arm")
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
        print("  # 1. 50 bp — the JCEC column must STILL be 0.476/0.243/0.386")
        print("  python3 62_build_s5.py --parts a 2>&1 | tail -40")
        print()
        print("  # 2. 125 bp")
        print("  TBL_SUFFIX=_rl125 \\")
        print("  S5_SLOPE=Colorado=0.5999,Brigham=0.5099,"
              "Stanford=0.4008,Fuwai=0.2948 \\")
        print("  S5_CONC=Colorado=80.8,Brigham=81.9,Fuwai=70.8,Stanford=75.8 \\")
        print("  S5_DISC_SE=exon=813,full=944 \\")
        print("  python3 62_build_s5.py --parts a 2>&1 | tail -40")
        print()
        print("⚠ THE OWN-EVENTS COLUMNS MUST NOT MOVE. This patch adds four")
        print("  columns and changes nothing about the existing fit. If JCEC")
        print("  or JC shifts even in the third decimal, the finite-value")
        print("  mask introduced here is dropping rows the old code kept and")
        print("  the patch must be reverted — the .bak is written beside it.")
        print()
        print("  Expected unchanged at 50 bp:  0.476/0.511, 0.243/0.256,")
        print("                                0.386/0.357")
        print("  Expected unchanged at 125 bp: 0.510/0.562, 0.295/0.306,")
        print("                                0.401/0.383")
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
    print("✓ verified on disk: (a) now fits each basis twice — on its own")
    print("  events and on the JC∩JCEC subset — the table carries JCEC_m,")
    print("  JC_m, change_m and n_m, S5.tsv gains a matched ROW per cohort")
    print("  with the schema unchanged, and both false lines are gone.")
    return 0


# ------------------------------------------------------------- self-test
# The fixture carries the REAL read_events and deming, copied verbatim from
# 62, so the matched arm is tested against the functions it will actually
# run against rather than against stand-ins.

_REAL_HELPERS = '''
NKEEP = {"SE": 3, "MXE": 5, "A3SS": 5, "A5SS": 5, "RI": 3}


def read_events(path, et, exon_level=True):
    if not (os.path.exists(path) and os.path.getsize(path)):
        return None, "missing"
    t = pd.read_csv(path, sep="\\t", low_memory=False)
    if "key" not in t.columns or "FDR" not in t.columns:
        return None, f"no key/FDR column; have {list(t.columns)[:8]}"
    d = next((c for c in ("dPSI_adj", "dPSI") if c in t.columns), None)
    if d is None:
        return None, "no dPSI column"
    if "se_adj" not in t.columns and "se" in t.columns:
        t = t.rename(columns={"se": "se_adj"})
    t["ekey"] = ([":".join(str(k).split(":")[:NKEEP[et]]) for k in t["key"]]
                 if exon_level else t["key"].astype(str))
    t = t.loc[t.groupby("ekey")["FDR"].idxmin()]
    return t.rename(columns={d: "dpsi"})[["ekey", "dpsi", "FDR"] +
                                         (["se_adj"] if "se_adj" in t.columns
                                          else [])], None


def first_existing(*paths):
    for p in paths:
        if os.path.exists(p) and os.path.getsize(p):
            return p
    return None


def deming(x, y, sx, sy):
    ok = np.isfinite(x) & np.isfinite(y) & np.isfinite(sx) & np.isfinite(sy)
    x, y, sx, sy = x[ok], y[ok], sx[ok], sy[ok]
    if len(x) < 10:
        return np.nan
    lam = np.mean(sy ** 2) / np.mean(sx ** 2) if np.mean(sx ** 2) > 0 else 1.0
    dx, dy = x - x.mean(), y - y.mean()
    Sxx, Syy, Sxy = (dx ** 2).sum(), (dy ** 2).sum(), (dx * dy).sum()
    if Sxy == 0:
        return np.nan
    return float((Syy - lam * Sxx + np.sqrt((Syy - lam * Sxx) ** 2
                                            + 4 * lam * Sxy ** 2)) / (2 * Sxy))
'''


def _fixture():
    """A module that contains all five anchors in their real context."""
    return (
        '"""\n' + A_DOC + '"""\n'
        "import os\n"
        "import numpy as np\n"
        "import pandas as pd\n"
        + _REAL_HELPERS +
        'SUF = os.environ.get("TBL_SUFFIX", "")\n'
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
        '    sig, _e = read_events("discovery.all_events.tsv", et)\n'
        "    out = {}\n"
        "    for coh, pat in REPL:\n"
        '        base = pat.replace("{ET}", et)\n'
        '        jcec = first_existing(f"{base}_dedup.all_events.tsv",\n'
        '                              f"{base}.all_events.tsv")\n'
        "        jc = first_existing(\n"
        '            f"jccheck{SUF}/{coh}_{et}_HCM_JC_dedup.all_events.tsv",\n'
        '            f"jccheck{SUF}/{coh}_{et}_HCM_JC.all_events.tsv")\n'
        # the 144 marker apply_edits requires, so the fixture is a
        # post-144 file like the one on the cluster
        '        print(f"  {coh:<10} JCEC -> {jcec or \'*** MISSING ***\'}")\n'
        '        print(f"  {\'\':<10} JC   -> '
        '{jc or \'*** MISSING (script 20 writes these) ***\'}")\n'
        "        out[coh] = (jc, jcec)\n"
        "    have = {c: v for c, v in out.items() if v[0] and v[1]}\n"
        "    gone = [c for c in out if c not in have]\n"
        "    if gone:\n"
        '        notes.append(f"(a) covers {len(have)}/{len(out)} cohorts")\n'
        "    res = {}\n"
        "    for coh, (jc, jcec) in have.items():\n"
        "        r = {}\n"
        + A_LOOP +
        "\n"
        "    print()\n"
        + A_PRINT
    )


def _fixture_tail():
    return (A_OK + A_ROWS)


def self_test():
    import tempfile

    fx = _fixture() + _fixture_tail() + (
        "\n"
        "rows, notes = [], []\n"
        "part_a(rows, notes)\n"
        'print("NROWS", len(rows))\n'
        "for _r in rows:\n"
        '    print("ROW|", _r["row"], "|", _r["label_1"], _r["value_1"],\n'
        '          "|", _r["label_2"], _r["value_2"], "|", _r["change"])\n'
        'print("NOTES", notes)\n')

    new, probs = apply_edits(fx)
    assert not probs, probs
    assert not verify(new), verify(new)
    compile(new, "fx", "exec")

    d = tempfile.mkdtemp()
    os.makedirs(os.path.join(d, "jccheck"), exist_ok=True)

    import numpy as np
    import pandas as pd

    rng = np.random.default_rng(0)
    n = 200
    keys = [f"chr1:{i}:{i + 50}:a:b" for i in range(n)]
    dx = rng.normal(0, 0.25, n)

    def write(path, idx, slope, noise):
        t = pd.DataFrame({
            "key": [keys[i] for i in idx],
            "FDR": 1e-4,
            "dPSI_adj": [slope * dx[i] + rng.normal(0, noise) for i in idx],
            "se_adj": 0.05,
        })
        fp = os.path.join(d, path)
        os.makedirs(os.path.dirname(fp) or d, exist_ok=True)
        t.to_csv(fp, sep="\t", index=False)

    # discovery: all 200
    pd.DataFrame({"key": keys, "FDR": 1e-4, "dPSI_adj": dx,
                  "se_adj": 0.05}).to_csv(
        os.path.join(d, "discovery.all_events.tsv"), sep="\t", index=False)

    # JCEC side: events 0..159, slope 0.50 throughout
    # JC side:   events 40..199, slope 0.50 on the overlap but 1.60 on its
    #            exclusive tail, so own != matched by a wide margin
    jcec_idx = list(range(0, 160))
    jc_overlap = list(range(40, 160))
    jc_tail = list(range(160, 200))
    for c in ("Brigham", "Fuwai", "Stanford"):
        write(f"{c}_SE_HCM.all_events.tsv", jcec_idx, 0.50, 0.01)
        t1 = pd.DataFrame({
            "key": [keys[i] for i in jc_overlap + jc_tail],
            "FDR": 1e-4,
            "dPSI_adj": ([0.50 * dx[i] + rng.normal(0, 0.01)
                          for i in jc_overlap]
                         + [1.60 * dx[i] + rng.normal(0, 0.01)
                            for i in jc_tail]),
            "se_adj": 0.05})
        t1.to_csv(os.path.join(d, f"jccheck/{c}_SE_HCM_JC.all_events.tsv"),
                  sep="\t", index=False)
    # Colorado: JCEC only, so it is `gone` and emits no row at all
    write("SE_Colorado_DCMvsNF.all_events.tsv", jcec_idx, 0.50, 0.01)

    p = os.path.join(d, "fx.py")
    open(p, "w").write(new)
    e = dict(os.environ)
    e.pop("TBL_SUFFIX", None)
    r = subprocess.run([sys.executable, "-I", p], capture_output=True,
                       text=True, cwd=d, env=e)
    assert r.returncode == 0, r.stdout + r.stderr
    o = r.stdout

    # --- the printed table gained the four columns
    for h in ("JCEC_m", "JC_m", "change_m", "n_m"):
        assert h in o, (h, o)

    # --- the matched subset is exactly the 120-event overlap
    assert o.count("    120") >= 3 or "120" in o, o

    # --- own and matched MUST differ: JC's own fit is dragged up by its
    #     exclusive steep tail, the matched fit is not
    rowlines = [l for l in o.splitlines() if l.startswith("ROW| ")]
    assert len(rowlines) == 6, rowlines          # 3 cohorts x (own, matched)
    own = {}
    mat = {}
    for l in rowlines:
        _, row, c1, c2, ch = [s.strip() for s in l.split("|")]
        name = row.split(",")[0]
        tgt = mat if "matched" in row else own
        tgt[name] = (float(c1.split()[-1]), float(c2.split()[-1]), float(ch))
    for c in ("Brigham", "Fuwai", "Stanford"):
        assert c in own and c in mat, (c, own, mat)
        # JCEC is unaffected by JC's tail, so JCEC_own ~= JCEC_m
        assert abs(own[c][0] - mat[c][0]) < 0.03, (c, own[c], mat[c])
        # JC own is pulled up; matched is not
        assert own[c][1] - mat[c][1] > 0.05, (c, own[c], mat[c])
        # and therefore change_m is much smaller than change
        assert abs(mat[c][2]) < abs(own[c][2]) - 0.05, (c, own[c], mat[c])
        # the matched pair is like-for-like, so it should be ~0 here
        assert abs(mat[c][2]) < 0.03, (c, mat[c])

    # --- Colorado has no JC table: no own row and no matched row
    assert not any("Colorado" in l for l in rowlines), rowlines

    # --- the matched row label carries n
    assert "matched subset (n=120)" in o, o

    # --- the two false lines are gone, the corrected advice is present
    assert "same event set and exclusions" not in o, o
    assert "QUOTE change_m" in o, o
    assert "DO NOT COMPARE THESE ROWS" in o, o
    assert "0.476 and 20's background JCEC is 0.662" in o, o

    # --- THE REGRESSION GUARD: the own-events slopes must be byte-identical
    #     to what the pre-patch code computes. This is the claim the dry run
    #     makes to CP, so it has to be tested, not asserted.
    pre = os.path.join(d, "fx_pre.py")
    open(pre, "w").write(fx)
    r0 = subprocess.run([sys.executable, "-I", pre], capture_output=True,
                        text=True, cwd=d, env=e)
    assert r0.returncode == 0, r0.stdout + r0.stderr
    pre_own = {}
    for l in r0.stdout.splitlines():
        if not l.startswith("ROW| "):
            continue
        _, row, c1, c2, ch = [s.strip() for s in l.split("|")]
        pre_own[row.split(",")[0]] = (float(c1.split()[-1]),
                                      float(c2.split()[-1]))
    assert set(pre_own) == set(own), (pre_own, own)
    for c in pre_own:
        assert pre_own[c] == own[c][:2], (c, pre_own[c], own[c])

    # --- refuses without 144
    _n, p4 = apply_edits("print('bare')\n")
    assert p4 and "144 has not been applied" in p4[0], p4

    # --- double application refused, and refuses without altering the text
    n6, p6 = apply_edits(new)
    assert p6 and "already carries 145" in p6[0], p6
    assert n6 == new, "a refused re-application still altered the text"

    print("[ok] self-test passed: the matched arm restricts both bases to the "
          "JC∩JCEC subset (n=120 of 160/160) and is fit with 62's own "
          "read_events and deming, copied verbatim; a steep JC-exclusive tail "
          "moves the own-events JC slope but NOT the matched one, so "
          "change_m is ~0 where change is not — the decomposition works; "
          "THE OWN-EVENTS SLOPES ARE BIT-IDENTICAL TO THE PRE-PATCH RUN; "
          "Colorado emits neither row; the matched row carries n in its "
          "label and the schema is unchanged at seven keys; both false lines "
          "are gone; refuses without 144; double application refused")
    return 0


if __name__ == "__main__":
    sys.exit(main())
