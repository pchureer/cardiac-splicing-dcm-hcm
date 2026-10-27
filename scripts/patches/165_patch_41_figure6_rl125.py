#!/usr/bin/env python3
r"""
165_patch_41_figure6_rl125.py -- Figure 6 at 125 bp, and panel D is retired.

41_figure6_enrichment.py takes --correct and --genome as flags, so its DATA is
already parameterised. Its NUMBERS are not. Six of them are module-level
literals that get drawn as TEXT on three panels:

    BG_TESTED = 12725          FG_ALL = 927
    FG_UP, FG_DOWN, FG_MIXED = 636, 291, 75
    FG_MULTI = 134
    FROZEN_TOTAL = {"correct": 4, "genome": 229}
    FROZEN_MULTI = {"correct": 0, "genome": 15}

At 125 bp they are 12,564 / 805 / 550 / 255 / 73 / 131, totals 1 and 85, multi
0 and 17. Point 41 at the _rl125 tables without changing these and it draws 50
bp constants beside 125 bp bars and says nothing -- the failure 37 and 38 both
had. Rule 10: a plotted constant is a write path.

So the override is MANDATORY, and the trigger is the DATA, not a flag anyone
could forget: if "rl125" appears anywhere on the command line and F6_COUNTS is
unset, 41 exits before plotting.

PANEL D IS RETIRED
------------------
    SURVIVOR = ("structural constituent of muscle", "GO:MF", 9.47e-05)
    ax.set_title("D   Survives adjustment: 15 genes", ...)

That term DOES NOT EXIST at 125 bp. The proper-background run goes from three
terms to one and both GO:MF terms are lost, this one among them. Panel D cannot
be re-based because its subject is gone, so it is no longer called and panel C
spans the bottom row instead.

This matches the manuscript, whose Figure 6 legend already describes exactly
three panels (A, B, C). Today the figure has four and the legend has three.

PANEL A'S RATIO
---------------
Panel A annotates "{fold}x more terms from the same {FG_ALL} genes". At 125 bp
that is 85 over 1 -- a fold change on a denominator of one, which the
Supplementary Discussion already rejected for the 452-fold and 1,767-fold
splicing ratios on identical grounds. The annotation now names the gene set and
the two backgrounds; the bars already carry both counts.

USAGE
-----
  python3 165_patch_41_figure6_rl125.py [--apply] [--self-test]
"""
import argparse
import difflib
import os
import shutil
import sys
import time

TARGET = "41_figure6_enrichment.py"

EDITS = []

# ---------------------------------------------------------------- 1. counts
EDITS.append((
    'the six plotted counts',
r'''BG_TESTED = 12725          # genes with >=1 tested event: the correct background
FG_ALL = 927               # significant genes, FDR<0.05 & |dPSI|>=0.10
FG_UP, FG_DOWN, FG_MIXED = 636, 291, 75
FG_MULTI = 134             # significant in >=2 event types''',
r'''import os as _os, sys as _sys

BG_TESTED = 12725          # genes with >=1 tested event: the correct background
FG_ALL = 927               # significant genes, FDR<0.05 & |dPSI|>=0.10
FG_UP, FG_DOWN, FG_MIXED = 636, 291, 75
FG_MULTI = 134             # significant in >=2 event types

# 165: these six are PLOTTED TEXT, not filenames -- panel A's annotation,
# panel B's axis label and panel C's header all render them. Rule 10: a
# plotted constant is a write path. The override is therefore MANDATORY, and
# the trigger is the DATA rather than a flag that could be forgotten.
_F6C = _os.environ.get("F6_COUNTS", "").strip()
if "rl125" in " ".join(_sys.argv) and not _F6C:
    _sys.exit(
        "[fatal] a 125 bp g:Profiler table was named on the command line, but\n"
        "        F6_COUNTS is unset. BG_TESTED/FG_ALL/FG_UP/FG_DOWN/FG_MIXED/\n"
        "        FG_MULTI would be drawn as 50 bp text beside 125 bp bars.\n"
        "        For 125 bp, set all three:\n"
        "          F6_COUNTS=12564/805/550/255/73/131\n"
        "          F6_TOTALS=1/85\n"
        "          F6_MULTI=0/17\n"
        "        Nothing plotted.")
if _F6C:
    BG_TESTED, FG_ALL, FG_UP, FG_DOWN, FG_MIXED, FG_MULTI = (
        int(_x) for _x in _F6C.split("/"))'''))

# ---------------------------------------------------------------- 2. totals
EDITS.append((
    "panel A's two bar totals",
r'''FROZEN_TOTAL = {"correct": 4, "genome": 229}''',
r'''FROZEN_TOTAL = {"correct": 4, "genome": 229}
_F6T = _os.environ.get("F6_TOTALS", "").strip()
if _F6C and not _F6T:
    _sys.exit("[fatal] F6_COUNTS is set but F6_TOTALS is not. Panel A's bars\n"
              "        would stay at the 50 bp 4 and 229. Nothing plotted.")
if _F6T:
    _a, _b = (int(_x) for _x in _F6T.split("/"))
    FROZEN_TOTAL = {"correct": _a, "genome": _b}'''))

# ----------------------------------------------------------------- 3. multi
EDITS.append((
    "panel C's two totals",
r'''FROZEN_MULTI = {"correct": 0, "genome": 15}''',
r'''FROZEN_MULTI = {"correct": 0, "genome": 15}
_F6M = _os.environ.get("F6_MULTI", "").strip()
if _F6C and not _F6M:
    _sys.exit("[fatal] F6_COUNTS is set but F6_MULTI is not. Panel C would keep\n"
              "        the 50 bp 15. Nothing plotted.")
if _F6M:
    _a, _b = (int(_x) for _x in _F6M.split("/"))
    FROZEN_MULTI = {"correct": _a, "genome": _b}'''))

# -------------------------------------------------------------- 4. SURVIVOR
EDITS.append((
    "a stale-term warning on SURVIVOR",
r'''SURVIVOR = ("structural constituent of muscle", "GO:MF", 9.47e-05)''',
r'''# 165: STALE AT 125 bp, AND NO LONGER DRAWN. "structural constituent of
# muscle" is one of the two GO:MF terms LOST when the proper-background run
# falls from three terms to one. panel_d is built entirely on it and is no
# longer called. Do not revive either without re-measuring.
SURVIVOR = ("structural constituent of muscle", "GO:MF", 9.47e-05)'''))

# ---------------------------------------------------------- 5. retire panel D
EDITS.append((
    'panel D retired; panel C spans the bottom row',
r'''    panel_c(fig.add_subplot(gs[1, 0]))
    panel_d(fig.add_subplot(gs[1, 1]))''',
r'''    # 165: panel D retired -- its term does not exist at 125 bp. C spans the
    # bottom row, giving the three-panel figure the manuscript legend already
    # describes.
    panel_c(fig.add_subplot(gs[1, :]))'''))

# ------------------------------------------------------- 6. panel A's ratio
EDITS.append((
    "panel A's fold-change annotation",
r'''        ax.annotate(f"{fold:.0f}× more terms\nfrom the same {FG_ALL} genes",''',
r'''        # 165: no ratio. At 125 bp this is 85 over 1, and a fold change on a
        # denominator of one is what the Supplementary Discussion already
        # rejected for 452-fold and 1,767-fold. The bars carry both counts.
        ax.annotate(f"the same {FG_ALL} genes,\ntwo backgrounds",'''))


def apply_edits(text):
    if "165:" in text:
        return text, ["41 already carries 165 (restore the .bak to re-apply)"]
    probs = []
    for name, a, b in EDITS:
        n = text.count(a)
        if n != 1:
            probs.append(f"{name}: anchor occurs {n}x, expected 1x")
    if probs:
        return text, probs
    for name, a, b in EDITS:
        text = text.replace(a, b, 1)
    return text, []


def verify(text):
    bad = []
    for w in ('_F6C = _os.environ.get("F6_COUNTS", "").strip()',
              'F6_COUNTS is unset',
              '_F6T = _os.environ.get("F6_TOTALS", "").strip()',
              '_F6M = _os.environ.get("F6_MULTI", "").strip()',
              'panel_c(fig.add_subplot(gs[1, :]))',
              'two backgrounds'):
        if w not in text:
            bad.append(f"missing: {w[:52]}")
    for w in ('panel_d(fig.add_subplot(gs[1, 1]))',
              'panel_c(fig.add_subplot(gs[1, 0]))',
              '× more terms'):
        if w in text:
            bad.append(f"survived, should be gone: {w[:46]}")
    # the literals themselves must still be there as the 50 bp defaults
    for w in ("BG_TESTED = 12725", "FG_ALL = 927", "FG_MULTI = 134",
              'FROZEN_TOTAL = {"correct": 4, "genome": 229}',
              'FROZEN_MULTI = {"correct": 0, "genome": 15}'):
        if w not in text:
            bad.append(f"the 50 bp default was destroyed: {w[:46]}")
    return bad


FIX = r'''import sys
BG_TESTED = 12725          # genes with >=1 tested event: the correct background
FG_ALL = 927               # significant genes, FDR<0.05 & |dPSI|>=0.10
FG_UP, FG_DOWN, FG_MIXED = 636, 291, 75
FG_MULTI = 134             # significant in >=2 event types

FROZEN_TOTAL = {"correct": 4, "genome": 229}
FROZEN_MULTI = {"correct": 0, "genome": 15}
SURVIVOR = ("structural constituent of muscle", "GO:MF", 9.47e-05)


def panel_a(ax, counts):
    fold = counts["genome"] / max(counts["correct"], 1)
    if True:
        ax.annotate(f"{fold:.0f}× more terms\nfrom the same {FG_ALL} genes",
                    xy=(0, 0))


def panel_c(ax): pass
def panel_d(ax): pass


def build(fig, gs, counts, C=None, G=None):
    panel_a(fig.add_subplot(gs[0, 0]), counts)
    panel_b(fig.add_subplot(gs[0, 1]), C, G, 0.05)
    panel_c(fig.add_subplot(gs[1, 0]))
    panel_d(fig.add_subplot(gs[1, 1]))
'''


def self_test():
    f = []
    new, probs = apply_edits(FIX)
    if probs:
        f.append(f"fixture did not patch: {probs}")
    else:
        bad = verify(new)
        if bad:
            f.append(f"verify rejected its own output: {bad}")
        try:
            compile(new, "x", "exec")
        except SyntaxError as e:
            f.append(f"does not compile: {e}")

        def run(env, argv):
            for k in ("F6_COUNTS", "F6_TOTALS", "F6_MULTI"):
                os.environ.pop(k, None)
            os.environ.update(env)
            old_argv = sys.argv
            sys.argv = argv
            ns = {}
            try:
                exec(compile(new, "x", "exec"), ns)
                return ns, None
            except SystemExit as e:
                return None, str(e)
            finally:
                sys.argv = old_argv
                for k in ("F6_COUNTS", "F6_TOTALS", "F6_MULTI"):
                    os.environ.pop(k, None)

        # 1. 50 bp, no env: the published constants survive untouched
        ns, err = run({}, ["41_figure6_enrichment.py", "--correct", "compadj_gprofiler.tsv"])
        if err:
            f.append(f"50 bp control refused when it should not: {err[:60]}")
        elif (ns["BG_TESTED"], ns["FG_ALL"], ns["FG_MULTI"]) != (12725, 927, 134):
            f.append("50 bp control did not keep 12725/927/134")
        elif ns["FROZEN_TOTAL"] != {"correct": 4, "genome": 229}:
            f.append("50 bp control did not keep FROZEN_TOTAL 4/229")

        # 2. 125 bp data, no env: MUST refuse
        ns, err = run({}, ["41_figure6_enrichment.py", "--correct",
                           "compadj_rl125_gprofiler.tsv"])
        if ns is not None or not err or "F6_COUNTS is unset" not in err:
            f.append("125 bp data WITHOUT F6_COUNTS did not refuse")

        # 3. counts set but totals missing: MUST refuse
        ns, err = run({"F6_COUNTS": "12564/805/550/255/73/131"},
                      ["41_figure6_enrichment.py"])
        if ns is not None or not err or "F6_TOTALS is not" not in err:
            f.append("F6_COUNTS without F6_TOTALS did not refuse")

        # 4. totals set but multi missing: MUST refuse
        ns, err = run({"F6_COUNTS": "12564/805/550/255/73/131", "F6_TOTALS": "1/85"},
                      ["41_figure6_enrichment.py"])
        if ns is not None or not err or "F6_MULTI is not" not in err:
            f.append("F6_TOTALS without F6_MULTI did not refuse")

        # 5. all three: every constant takes the 125 bp value
        ns, err = run({"F6_COUNTS": "12564/805/550/255/73/131",
                       "F6_TOTALS": "1/85", "F6_MULTI": "0/17"},
                      ["41_figure6_enrichment.py", "--correct",
                       "compadj_rl125_gprofiler.tsv"])
        if err:
            f.append(f"fully overridden run refused: {err[:70]}")
        else:
            got = (ns["BG_TESTED"], ns["FG_ALL"], ns["FG_UP"], ns["FG_DOWN"],
                   ns["FG_MIXED"], ns["FG_MULTI"])
            if got != (12564, 805, 550, 255, 73, 131):
                f.append(f"125 bp counts wrong: {got}")
            if ns["FROZEN_TOTAL"] != {"correct": 1, "genome": 85}:
                f.append(f"125 bp totals wrong: {ns['FROZEN_TOTAL']}")
            if ns["FROZEN_MULTI"] != {"correct": 0, "genome": 17}:
                f.append(f"125 bp multi wrong: {ns['FROZEN_MULTI']}")

        again, p2 = apply_edits(new)
        if not p2 or "already carries 165" not in p2[0]:
            f.append("re-application did not refuse")

    print("165 self-test")
    for x in f:
        print("  ✗ " + x)
    if not f:
        print("  ✓ patches and compiles; the 50 bp control keeps 12725/927/134")
        print("    and 4/229 untouched; 125 bp DATA without F6_COUNTS refuses;")
        print("    a half-set override refuses at both steps; a full override")
        print("    gives 12564/805/550/255/73/131, 1/85, 0/17; panel D is no")
        print("    longer called and C spans; the ratio is gone; re-apply refuses")
    return 1 if f else 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("target", nargs="?", default=TARGET)
    ap.add_argument("--apply", action="store_true")
    ap.add_argument("--self-test", action="store_true")
    a = ap.parse_args()
    if a.self_test:
        return self_test()
    if self_test():
        print("\nREFUSING: 165's own self-test failed.")
        return 3
    print()
    if not os.path.exists(a.target):
        print(f"REFUSING: {a.target} not found")
        return 3
    src = open(a.target, errors="replace").read()
    new, probs = apply_edits(src)
    print("=" * 74)
    print("165 — Figure 6 at 125 bp; panel D retired")
    print("=" * 74)
    if probs:
        print("REFUSING:")
        for p in probs:
            print("    " + p)
        print("\n  An anchor that does not appear exactly once means the file is")
        print("  not the one this patch was written against. Paste the region.")
        return 3
    print("".join(difflib.unified_diff(src.splitlines(True), new.splitlines(True),
                                       fromfile=a.target,
                                       tofile=a.target + " (patched)", n=2)))
    if not a.apply:
        print("DRY RUN. Re-run with --apply, then:")
        print()
        print("  # 50 bp control — must reproduce the published figure")
        print("  python3 41_figure6_enrichment.py \\")
        print("      --correct compadj_gprofiler.tsv \\")
        print("      --genome  compadj_gprofiler_WHOLEGENOME.tsv \\")
        print("      --out /tmp/f6_cal")
        print()
        print("  # the guard — this MUST refuse, nothing written")
        print("  python3 41_figure6_enrichment.py \\")
        print("      --correct compadj_rl125_gprofiler.tsv \\")
        print("      --genome  compadj_rl125_gprofiler_WHOLEGENOME.tsv \\")
        print("      --out /tmp/f6_should_not_exist")
        print()
        print("  # 125 bp")
        print("  F6_COUNTS=12564/805/550/255/73/131 F6_TOTALS=1/85 F6_MULTI=0/17 \\")
        print("  python3 41_figure6_enrichment.py \\")
        print("      --correct compadj_rl125_gprofiler.tsv \\")
        print("      --genome  compadj_rl125_gprofiler_WHOLEGENOME.tsv \\")
        print("      --out F6_enrichment_rl125")
        print()
        print("  Pass BOTH files explicitly. 41's own line 186 warns that")
        print("  auto-discovery can leave panel A half frozen with no warning,")
        print("  and there are now 50 bp and 125 bp tables side by side.")
        return 0
    bak = f"{a.target}.bak.{time.strftime('%Y%m%d_%H%M%S')}"
    shutil.copy2(a.target, bak)
    open(a.target, "w").write(new)
    print(f"[write] {a.target}   (backup: {bak})")
    disk = open(a.target, errors="replace").read()
    bad = verify(disk)
    try:
        compile(disk, a.target, "exec")
    except SyntaxError as e:
        bad.append(f"does not compile: {e}")
    if bad or disk != new:
        print("✗ ON-DISK VERIFICATION FAILED — restoring:")
        for b in bad:
            print(f"    {b}")
        shutil.copy2(bak, a.target)
        return 3
    print("✓ verified: overrides are mandatory on 125 bp data, panel D is no")
    print("  longer called, C spans the bottom row, and the ratio is gone.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
