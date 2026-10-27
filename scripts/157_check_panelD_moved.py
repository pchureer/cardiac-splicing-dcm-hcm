#!/usr/bin/env python3
"""
157_check_panelD_moved.py -- the one check 156's own dry run asked for.

WHY THIS EXISTS
---------------
156's dry run ended with:

    ⚠ Check panel D's printed numbers CHANGE between runs 1 and 3.
      If they do not, TECH45 is still being read from the literals
      and the figure is still mixed.

That check was never made. The run-3 output we looked at was a `tail -22`
window, which showed the PROVENANCE block (`SE (641 events), 39 for MXE
(1,399)`, `rbm20_g45_407951.out step 5`) but not panel D's own summary rows.
The provenance block and panel D are **separate write paths** — 156 edited
both, but they are read from different variables (`_NT_SE`/`_D_LOG` vs
`TECH45`/`TECH9`). A provenance block that moved proves `_f5_env` works; it
does **not** prove `TECH45` was read from the environment rather than from
the literals. That is exactly the failure 156 was written to catch, so
confirming it by eye on the one output 156 named is not optional.

WHAT IT CHECKS
--------------
Against the eight-value table 156 derived and verified (50 bp re-derived from
rbm20_g45_406229.out before the 125 bp values were transcribed from
rbm20_g45_407951.out):

                      50 bp            125 bp
  SE  45, without     0.95 -> 0.86     0.96 -> 0.87
  SE  45, with        1.02 -> 0.89     1.02 -> 0.92
  MXE 45, without     2.03 -> 1.23     1.98 -> 1.24
  MXE 45, with        2.71 -> 1.63     2.65 -> 1.63
  SE  nine, without   1.01 -> 0.93     1.03 -> 0.95
  SE  nine, with      1.09 -> 0.96     1.08 -> 1.00
  MXE nine, without   2.01 -> 1.28     1.95 -> 1.29
  MXE nine, with      2.65 -> 1.65     2.52 -> 1.67

The two sets of eight PAIRS are disjoint, which is what makes this a usable
test. Individual numbers are not: `1.02` is the SE-with "before" at both read
lengths, and `1.63` is an "after" in both tables under different rows. So the
unit of comparison is the ordered pair, never a single float. (`0.281` and
`0.690` and `864` have each already collided once in this project.)

It refuses rather than reports when the inputs cannot support a verdict:
two identical logs, two logs from the same run, or logs passed in the wrong
order (a 50 bp log containing 125 bp pairs). A negative result from a
comparison of a file against itself is not evidence.

  --pdf50/--pdf125 additionally read the PLOTTED text out of the figures with
  pdftotext, which is the actual deliverable; the printed rows are a proxy
  for it. Absent pdftotext the check degrades to the logs alone and says so.

USAGE
-----
  # on the cluster, after re-running steps 1 and 3 with tee (see below)
  python3 157_check_panelD_moved.py --log50 f5_rl50.log --log125 f5_rl125.log \\
      [--pdf50 legend_check/f5_56check.pdf --pdf125 legend_check/f5_rl125.pdf]

  python3 157_check_panelD_moved.py --self-test
"""
import argparse
import os
import re
import shutil
import subprocess
import sys

# ---------------------------------------------------------------- the table
# key -> (before, after), as "%.2f" strings so 1.00/1.0/1 all compare equal.
TECH = {
    50: {
        ("SE", 45, "without"): ("0.95", "0.86"),
        ("SE", 45, "with"): ("1.02", "0.89"),
        ("MXE", 45, "without"): ("2.03", "1.23"),
        ("MXE", 45, "with"): ("2.71", "1.63"),
        ("SE", 9, "without"): ("1.01", "0.93"),
        ("SE", 9, "with"): ("1.09", "0.96"),
        ("MXE", 9, "without"): ("2.01", "1.28"),
        ("MXE", 9, "with"): ("2.65", "1.65"),
    },
    125: {
        ("SE", 45, "without"): ("0.96", "0.87"),
        ("SE", 45, "with"): ("1.02", "0.92"),
        ("MXE", 45, "without"): ("1.98", "1.24"),
        ("MXE", 45, "with"): ("2.65", "1.63"),
        ("SE", 9, "without"): ("1.03", "0.95"),
        ("SE", 9, "with"): ("1.08", "1.00"),
        ("MXE", 9, "without"): ("1.95", "1.29"),
        ("MXE", 9, "with"): ("2.52", "1.67"),
    },
}
PAIRS = {rl: set(TECH[rl].values()) for rl in TECH}

# the provenance block is a SECOND write path; 156 edited it from different
# variables, so it is checked separately and never used to stand in for panel D
PROV = {
    50: {"must": ["663 event", "1,393", "406229"],
         "must_not": ["641 event", "1,399", "407951"]},
    125: {"must": ["641 event", "1,399", "407951"],
          "must_not": ["663 event", "1,393", "406229"]},
}

ARROW = re.compile(r"(-?\d+(?:\.\d+)?)\s*(?:->|-->|→|\bto\b)\s*(-?\d+(?:\.\d+)?)")
FLOAT = re.compile(r"-?\d+\.\d+")


def _norm(x):
    try:
        return "%.2f" % float(x)
    except ValueError:
        return None


def extract_pairs(text):
    """Every ordered float pair the text states as 'X -> Y', normalised to 2dp.

    Format-agnostic on purpose: 90's panel D summary rows are matched by the
    arrow, not by a header string. A header can be reworded; `X -> Y` is the
    number itself. Returns (pairs, lines_that_produced_them).
    """
    pairs, where = set(), []
    for ln in text.splitlines():
        for m in ARROW.finditer(ln):
            a, b = _norm(m.group(1)), _norm(m.group(2))
            if a is None or b is None:
                continue
            pairs.add((a, b))
            where.append((a, b, ln.strip()[:96]))
    return pairs, where


def extract_pairs_fallback(text):
    """If no arrows are found, pair floats positionally on 'composition' rows.

    Used only to produce a usable diagnostic, never a PASS: a positional pair
    is a guess at the format, and a guess cannot license a verdict.
    """
    out = []
    for ln in text.splitlines():
        if "composition" not in ln.lower() and "median" not in ln.lower():
            continue
        f = [_norm(x) for x in FLOAT.findall(ln)]
        if len(f) >= 2:
            out.append((f[0], f[1], ln.strip()[:96]))
    return out


def pdftotext(path):
    exe = shutil.which("pdftotext")
    if not exe:
        return None, "pdftotext not on PATH"
    if not os.path.exists(path):
        return None, f"no such file: {path}"
    try:
        r = subprocess.run([exe, "-layout", path, "-"],
                           capture_output=True, text=True, timeout=120)
    except Exception as e:                                   # noqa: BLE001
        return None, f"pdftotext failed: {e}"
    if r.returncode != 0:
        return None, f"pdftotext exit {r.returncode}: {r.stderr.strip()[:120]}"
    return r.stdout, None


def _table(p50, p125):
    out = ["panel D's PLOTTED constants",
           "-" * 62,
           f"{'row':<22}{'50 bp':<18}{'125 bp':<18}moved"]
    for k in TECH[50]:
        a, b = TECH[50][k]
        c, d = TECH[125][k]
        tag = f"{k[0]} {k[1]}, {k[2]}"
        ok50, ok125 = (a, b) in p50, (c, d) in p125
        mark = "yes" if (ok50 and ok125) else "NO"
        note = ""
        if not ok50:
            note += "   <-- 50 bp value absent"
        if not ok125:
            note += "   <-- 125 bp value absent"
        out.append(f"{tag:<22}{a + ' -> ' + b:<18}{c + ' -> ' + d:<18}"
                   f"{mark}{note}")
    out.append("")
    return out


def _prov(t50, t125, name50, name125):
    """The provenance block is its own write path and is never allowed to
    stand in for panel D: 156 edited both, from different variables."""
    bad = []
    for rl, (nm, txt) in ((50, (name50, t50)), (125, (name125, t125))):
        for s in PROV[rl]["must"]:
            if s not in txt:
                bad.append(f"{nm}: provenance is missing '{s}'")
        for s in PROV[rl]["must_not"]:
            if s in txt:
                bad.append(f"{nm}: provenance still carries '{s}'")
    out = ["provenance block (a SEPARATE write path from panel D)", "-" * 62]
    out += ["  ✗ " + b for b in bad] if bad else \
           ["  ✓ 50 bp names 663/1,393/406229; "
            "125 bp names 641/1,399/407951"]
    out.append("")
    return bad, out


# ---------------------------------------------------------------- the check
def check(log50_text, log125_text, name50="log50", name125="log125"):
    """Returns (code, lines). code 0 pass, 1 fail, 2 inconclusive, 3 refuse."""
    out, refuse = [], []

    if log50_text == log125_text:
        refuse.append(
            "the two logs are byte-identical. Either the same file was passed\n"
            "    twice or the same run was tee'd twice. A file compared against\n"
            "    itself cannot show a change.")
    for nm, txt, other_rl in ((name50, log50_text, 125), (name125, log125_text, 50)):
        if not txt.strip():
            refuse.append(f"{nm} is empty")
    if refuse:
        return 3, ["REFUSING — the inputs cannot support a verdict:"] + \
                  ["  * " + r for r in refuse]

    p50, w50 = extract_pairs(log50_text)
    p125, w125 = extract_pairs(log125_text)

    # order check BEFORE anything else: a swap would read as a clean pass on
    # the 125 bp side and a failure on the 50 bp side, which is the most
    # misleading output this script could produce.
    swapped = (PAIRS[125] & p50) and (PAIRS[50] & p125) and not (PAIRS[50] & p50)
    if swapped:
        return 3, ["REFUSING — the logs look swapped:",
                   f"  * {name50} carries 125 bp panel-D pairs "
                   f"({sorted(PAIRS[125] & p50)[:3]})",
                   f"  * {name125} carries 50 bp pairs "
                   f"({sorted(PAIRS[50] & p125)[:3]})",
                   "  Re-run with --log50 and --log125 the other way round."]

    hit50 = PAIRS[50] & p50
    hit125 = PAIRS[125] & p125

    # cross-contamination is decided BEFORE completeness. A 125 bp log that
    # prints the 50 bp panel has 0/8 of its own pairs, so a completeness test
    # run first calls the one failure this script exists to catch
    # "inconclusive" and sends it back for a better capture.
    stale125 = PAIRS[50] & p125
    stale50 = PAIRS[125] & p50
    if stale125 or stale50:
        out.append("✗ FAIL — a log prints the OTHER read length's panel D.")
        out.append("")
        out.extend(_table(p50, p125))
        if stale125:
            out.append(f"  {name125} still prints 50 bp panel-D pairs:")
            for pr in sorted(stale125):
                out.append(f"      {pr[0]} -> {pr[1]}")
            out.append("  TECH45/TECH9 are still being read from the literals.")
            out.append("  The figure on disk mixes read lengths. Check that")
            out.append("  F5_TECH45 and F5_TECH9 were actually exported into")
            out.append("  the run-3 environment — they sit on backslash-")
            out.append("  continued lines, and one missing backslash drops")
            out.append("  both while every gate still passes.")
        if stale50:
            out.append(f"  {name50} prints 125 bp panel-D pairs: "
                       f"{sorted(stale50)}")
            out.append("  The 50 bp control is no longer a control.")
        out.append("")
        out.extend(_prov(log50_text, log125_text, name50, name125)[1])
        return 1, out

    if len(hit50) < 8 or len(hit125) < 8:
        out.append("INCONCLUSIVE — panel D's rows were not found as 'X -> Y'.")
        out.append(f"  {name50}:  {len(hit50)}/8 expected 50 bp pairs found")
        out.append(f"  {name125}: {len(hit125)}/8 expected 125 bp pairs found")
        for nm, txt in ((name50, log50_text), (name125, log125_text)):
            fb = extract_pairs_fallback(txt)
            out.append(f"  candidate rows in {nm} (format guess, NOT a verdict):")
            for a, b, ln in fb[:10]:
                out.append(f"      {a} / {b}   | {ln}")
            if not fb:
                out.append("      (none — did the log capture the full run? a")
                out.append("       `tail -22` window cuts panel D off)")
        out.append("")
        out.append("  Fix the capture (tee the WHOLE run, not a tail) and re-run.")
        out.append("  If the rows are present but worded differently, send one")
        out.append("  of them and 157's extractor gets the real format.")
        return 2, out

    # the actual question
    out.extend(_table(p50, p125))
    prov_bad, prov_lines = _prov(log50_text, log125_text, name50, name125)
    out.extend(prov_lines)

    if prov_bad:
        return 1, out
    if hit50 == hit125:
        out.append("✗ FAIL — panel D did not move between the runs.")
        return 1, out
    out.append("✓ PASS — all eight panel-D pairs moved, each log carries only")
    out.append("  its own read length's values, and the provenance block")
    out.append("  moved independently. 156 is confirmed on the figure, not")
    out.append("  only on the gates.")
    return 0, out


def check_pdfs(pdf50, pdf125):
    out, verdict = [], None
    out.append("the figure itself (pdftotext — the plotted text, not the print)")
    out.append("-" * 60)
    txts = {}
    for rl, p in ((50, pdf50), (125, pdf125)):
        if not p:
            continue
        t, err = pdftotext(p)
        if err:
            out.append(f"  skipped {rl} bp: {err}")
            continue
        txts[rl] = t
    if len(txts) < 2:
        out.append("  degraded to the logs alone — the printed rows are a PROXY")
        out.append("  for the figure, so this leaves the plotted panel unverified")
        return None, out
    for rl in (50, 125):
        own = PAIRS[rl] & extract_pairs(txts[rl])[0]
        other = PAIRS[125 if rl == 50 else 50] & extract_pairs(txts[rl])[0]
        out.append(f"  {rl} bp pdf: {len(own)}/8 own pairs, "
                   f"{len(other)} foreign pairs")
        if other:
            verdict = 1
            out.append(f"      ✗ foreign: {sorted(other)}")
    if verdict is None:
        n = sum(len(PAIRS[rl] & extract_pairs(txts[rl])[0]) for rl in (50, 125))
        if n == 0:
            out.append("  no panel-D pairs found as text in either pdf (the")
            out.append("  values may be drawn as marks, not labels) — the logs")
            out.append("  remain the evidence")
            return None, out
        verdict = 0
        out.append("  ✓ each figure carries only its own read length's values")
    return verdict, out


# ---------------------------------------------------------------- self-test
def _fixture(rl, prov=True, extra=""):
    t = TECH[rl]
    ls = ["building figure 5 ...",
          "G4 SE -0.051 (-0.159, +0.047) record -0.051 (-0.159, +0.047)  match",
          "",
          "panel D - 45-gene background medians (the published figure used nine-gene):"]
    for (et, n, w), (a, b) in t.items():
        ls.append(f"    {et:<4} {'nine-gene' if n == 9 else '45-gene'} "
                  f"{w} composition adj.   {a} -> {b}")
    ls.append("")
    ls.append("PROVENANCE - paste into the legend")
    if prov:
        p = PROV[rl]["must"]
        ls.append(f"  40 of 45 genes testable for SE ({p[0].split()[0]} events),"
                  f" 39 for MXE ({p[1]})")
        ls.append(f"  panel D  45-gene background medians ({p[2]} step 5)")
    ls.append(extra)
    return "\n".join(ls) + "\n"


def self_test():
    fails = []

    def expect(name, code, got, lines):
        if got != code:
            fails.append(f"{name}: expected code {code}, got {got}\n      "
                         + "\n      ".join(lines[:6]))

    f50, f125 = _fixture(50), _fixture(125)

    c, l = check(f50, f125)
    expect("moved -> pass", 0, c, l)

    # the failure 156 exists to catch: TECH read from the literals, so the
    # 125 bp run prints the 50 bp panel while its provenance moved
    mixed = _fixture(50, prov=False) + "\n".join(
        [f"  40 of 45 genes testable for SE (641 events), 39 for MXE (1,399)",
         "  panel D  45-gene background medians (rbm20_g45_407951.out step 5)"])
    c, l = check(f50, mixed)
    expect("125 bp prints the 50 bp panel -> fail", 1, c, l)
    if not any("still prints 50 bp panel-D pairs" in x for x in l):
        fails.append("the mixed case did not name the stale pairs")

    # identical inputs must refuse, not pass and not fail
    c, l = check(f50, f50)
    expect("identical logs -> refuse", 3, c, l)

    # swapped order must refuse rather than report a half-pass
    c, l = check(f125, f50)
    expect("swapped logs -> refuse", 3, c, l)

    # a tail-truncated capture must be inconclusive, with the reason named
    c, l = check("panel D ... (cut off)\nPROVENANCE\n  SE (663 events) 1,393 406229\n",
                 "panel D ... (cut off)\nPROVENANCE\n  SE (641 events) 1,399 407951\n")
    expect("truncated capture -> inconclusive", 2, c, l)
    if not any("tail -22" in x for x in l):
        fails.append("the truncated case did not point at the capture")

    # panel D right, provenance not re-based -> fail, and the two must be
    # reported separately (they are different write paths)
    c, l = check(f50, _fixture(125, prov=False)
                 + "  40 of 45 genes testable for SE (663 events), 39 for MXE (1,393)\n"
                   "  panel D  45-gene background medians (rbm20_g45_406229.out step 5)\n")
    expect("provenance stale -> fail", 1, c, l)
    if not any("provenance still carries" in x for x in l):
        fails.append("the stale-provenance case did not name it")

    # a pair must be the unit: 1.02 and 1.63 and 2.65 each appear in BOTH
    # tables as single values, so a single-float test would mis-call this
    for v in ("1.02", "1.63", "2.65"):
        in50 = any(v in p for p in PAIRS[50])
        in125 = any(v in p for p in PAIRS[125])
        if not (in50 and in125):
            fails.append(f"fixture assumption broken: {v} not in both tables")
    if PAIRS[50] & PAIRS[125]:
        fails.append("the two pair sets overlap — the test is not discriminating")

    print("157 self-test")
    for f in fails:
        print("  ✗ " + f)
    if not fails:
        print("  ✓ 8 cases: pass, mixed-panel fail, identical refuse, swapped")
        print("    refuse, truncated inconclusive, stale provenance fail, and")
        print("    the pair-not-float and disjointness assumptions hold")
    return 1 if fails else 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--log50")
    ap.add_argument("--log125")
    ap.add_argument("--pdf50")
    ap.add_argument("--pdf125")
    ap.add_argument("--self-test", action="store_true")
    a = ap.parse_args()

    if a.self_test:
        return self_test()
    if self_test():
        print("REFUSING: 157's own self-test failed.")
        return 3
    print()
    if not (a.log50 and a.log125):
        print("REFUSING: --log50 and --log125 are both required.")
        print("  The question is whether panel D CHANGED, which no single")
        print("  log can answer. See the header for the two tee'd re-runs.")
        return 3
    for p in (a.log50, a.log125):
        if not os.path.exists(p):
            print(f"REFUSING: no such file: {p}")
            return 3
    if os.path.realpath(a.log50) == os.path.realpath(a.log125):
        print("REFUSING: --log50 and --log125 are the same file.")
        return 3

    t50 = open(a.log50, errors="replace").read()
    t125 = open(a.log125, errors="replace").read()
    print("=" * 70)
    print("157 — did panel D actually move? (156's own open check)")
    print("=" * 70)
    print(f"  50 bp log : {a.log50}  ({len(t50.splitlines())} lines)")
    print(f"  125 bp log: {a.log125}  ({len(t125.splitlines())} lines)")
    print()
    code, lines = check(t50, t125, a.log50, a.log125)
    for ln in lines:
        print(ln)

    if a.pdf50 or a.pdf125:
        print()
        pcode, plines = check_pdfs(a.pdf50, a.pdf125)
        for ln in plines:
            print(ln)
        if pcode == 1:
            code = 1
    print()
    print("=" * 70)
    print({0: "PASS", 1: "FAIL", 2: "INCONCLUSIVE", 3: "REFUSED"}[code])
    return code


if __name__ == "__main__":
    sys.exit(main())
