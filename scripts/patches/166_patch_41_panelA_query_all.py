#!/usr/bin/env python3
r"""
166_patch_41_figure6_panelA.py -- panel A counts every query list. It should
count one.

THE DEFECT
----------
g:Profiler writes one row per (query list, term), and this project submits FOUR
lists. The 125 bp whole-genome table holds:

    85 all · 69 up · 8 down · 17 multi_event_type  =  179 rows

41 lines 474/477 count rows with no filter:

    counts["correct"] = int((C.p < a.alpha).sum())
    counts["genome"]  = int((G.p < a.alpha).sum())

so panel A draws 179 where the manuscript -- abstract, Results, Discussion and
legend -- says 85. Two consequences, the second worse than the first:

1. Panel A plots a quantity the paper never discusses.
2. It DOUBLE COUNTS. Panel C reports the multi_event_type list as 17 terms, and
   those same 17 rows are inside panel A's 179. One number, two panels,
   presented as different findings.

read_gp already parses a `query` column (it is in `out`), and panel B already
scopes itself to QUERY_SHOWN = {"all"}. Panel A simply never applied the same
restriction. This is a defect at BOTH read lengths -- 50 bp's frozen 4 is
likewise 3 all + 1 up -- so it was not caused by the re-quantification and is
not fixed by it.

WHAT THIS DOES *NOT* DO, DELIBERATELY
-------------------------------------
It does not make F6_TOTALS authoritative over the file. That was the first
instinct and it is wrong: line 503 compares the file-derived count against
FROZEN_TOTAL and prints a warning on disagreement, and that comparison is what
caught a bad F6_TOTALS in the first place. An override that always wins makes
the two agree by construction and silences the check permanently.

The honest fix is to make the MEASURED value and the DECLARED value describe
the same quantity, so that when they disagree it means something.

⚠ THE 50 bp CONTROL WILL CHANGE, AND THAT IS THE POINT
------------------------------------------------------
FROZEN_TOTAL = {"correct": 4, "genome": 229} are four-list totals. After this
patch the 50 bp run reports the `all`-list counts instead and the existing
warning WILL fire, naming both numbers. That is a pre-existing defect becoming
visible, not a regression. Do not silence it by editing FROZEN_TOTAL to match;
record what it prints.

USAGE
-----
  python3 166_patch_41_panelA_query_all.py [--apply] [--self-test]
"""
import argparse
import difflib
import os
import shutil
import sys
import time

TARGET = "41_figure6_enrichment.py"

HELPER = r'''def _n_all(df, alpha, verbose=True):
    """166: significant terms in the 'all' query list only.

    g:Profiler emits one row per (query list, term) and this project submits
    four lists -- all, up, down, multi_event_type. Counting every row answers a
    question the manuscript never asks AND double counts, because panel C
    reports the multi_event_type rows separately. Panel B already restricts
    itself to QUERY_SHOWN = {"all"}; panel A did not.
    """
    if df is None:
        return 0
    q = df.get("query")
    if q is None or not q.notna().any():
        if verbose:
            print("  [warn] no 'query' column in a g:Profiler table: panel A is "
                  "counting EVERY query list, which is the 179-vs-85 defect. "
                  "Check the number before publishing it.")
        return int((df.p < alpha).sum())
    m = q.astype(str).str.strip().str.lower().eq("all")
    if not m.any():
        raise SystemExit(
            "[fatal] a g:Profiler table has a 'query' column but no rows with\n"
            "        query == 'all'. Panel A cannot be scoped and would count\n"
            "        every list. Nothing plotted.")
    return int((df[m].p < alpha).sum())


'''

EDITS = [
    ("the _n_all helper",
     "def discover(explicit_correct, explicit_genome, verbose=True):",
     HELPER + "def discover(explicit_correct, explicit_genome, verbose=True):"),
    ("panel A's proper-background count",
     '        counts["correct"] = int((C.p < a.alpha).sum())',
     '        counts["correct"] = _n_all(C, a.alpha)'),
    ("panel A's whole-genome count",
     '        counts["genome"] = int((G.p < a.alpha).sum())',
     '        counts["genome"] = _n_all(G, a.alpha)'),
]


def apply_edits(text):
    if "166:" in text:
        return text, ["41 already carries 166 (restore the .bak to re-apply)"]
    probs = [f"{n}: anchor occurs {text.count(a)}x, expected 1x"
             for n, a, b in EDITS if text.count(a) != 1]
    if probs:
        return text, probs
    for n, a, b in EDITS:
        text = text.replace(a, b, 1)
    return text, []


def verify(text):
    bad = []
    for w in ("def _n_all(df, alpha, verbose=True):",
              'counts["correct"] = _n_all(C, a.alpha)',
              'counts["genome"] = _n_all(G, a.alpha)',
              'str.strip().str.lower().eq("all")'):
        if w not in text:
            bad.append(f"missing: {w[:52]}")
    for w in ('counts["correct"] = int((C.p < a.alpha).sum())',
              'counts["genome"] = int((G.p < a.alpha).sum())'):
        if w in text:
            bad.append(f"unscoped count survived: {w[:48]}")
    # the cross-check must NOT have been disabled
    if "terms at p<" not in text:
        bad.append("the file-vs-FROZEN_TOTAL warning was removed")
    for w in ("FROZEN_TOTAL", "QUERY_SHOWN"):
        if w not in text:
            bad.append(f"{w} was damaged")
    return bad


FIX = r'''import pandas as pd
import numpy as np

FROZEN_TOTAL = {"correct": 4, "genome": 229}
QUERY_SHOWN = {"all"}


def discover(explicit_correct, explicit_genome, verbose=True):
    return None, None, None, None


def count(C, G, a_alpha):
    counts = {}
    if C is not None:
        counts["correct"] = int((C.p < a_alpha).sum())
    if G is not None:
        counts["genome"] = int((G.p < a_alpha).sum())
    return counts


def warn(counts, key, fz, alpha):
    if counts[key] != fz:
        print(f"  *** {key} background: {counts[key]} terms at p<{alpha}")
'''


def _frame(rows):
    import pandas as pd
    return pd.DataFrame(rows, columns=["term", "p", "query"])


def self_test():
    # the no-query-column case prints its warning here; that is the fixture
    # exercising the fallback, not a problem with the patch.
    f = []
    # the fixture's count() body uses `a_alpha`, so patch against the real names
    fix = FIX.replace("a_alpha", "a.alpha").replace(
        "def count(C, G, a.alpha):", "def count(C, G, a):")
    # rebuild cleanly: the anchors must appear verbatim
    fix = FIX
    fix = fix.replace('counts["correct"] = int((C.p < a_alpha).sum())',
                      '        counts["correct"] = int((C.p < a.alpha).sum())'.strip())
    # simplest: construct the fixture with the exact anchor text
    fix = r'''import pandas as pd
import numpy as np

FROZEN_TOTAL = {"correct": 4, "genome": 229}
QUERY_SHOWN = {"all"}


def discover(explicit_correct, explicit_genome, verbose=True):
    return None, None, None, None


class A:
    alpha = 0.05


def count(C, G):
    a = A()
    counts = {}
    if True:
        counts["correct"] = int((C.p < a.alpha).sum())
    if True:
        counts["genome"] = int((G.p < a.alpha).sum())
    return counts


def warn(counts, key, fz, alpha):
    if counts[key] != fz:
        print(f"  *** {key} background: {counts[key]} terms at p<{alpha}")
'''
    new, probs = apply_edits(fix)
    if probs:
        f.append(f"fixture did not patch: {probs}")
        print("166 self-test")
        for x in f:
            print("  ✗ " + x)
        return 1
    bad = verify(new)
    if bad:
        f.append(f"verify rejected its own output: {bad}")
    try:
        compile(new, "x", "exec")
    except SyntaxError as e:
        f.append(f"does not compile: {e}")
    else:
        ns = {}
        exec(compile(new, "x", "exec"), ns)
        # the real 125 bp shape: 85 all, 69 up, 8 down, 17 multi = 179
        G = _frame([("t%d" % i, 0.001, q) for q, k in
                    (("all", 85), ("up", 69), ("down", 8), ("multi_event_type", 17))
                    for i in range(k)])
        C = _frame([("Cytoskeleton in muscle cells", 0.0112, "all")])
        got = ns["count"](C, G)
        if got != {"correct": 1, "genome": 85}:
            f.append(f"expected correct 1 / genome 85, got {got}")
        if len(G) != 179:
            f.append("fixture is not 179 rows")
        # alpha still applies WITHIN the all list
        G2 = _frame([("a", 0.001, "all"), ("b", 0.9, "all"), ("c", 0.001, "up")])
        if ns["count"](C, G2)["genome"] != 1:
            f.append("alpha is not applied within the 'all' list")
        # a table with NO query column falls back and warns rather than crashing
        import pandas as pd
        G3 = pd.DataFrame({"term": ["a", "b"], "p": [0.001, 0.002]})
        G3["query"] = float("nan")
        if ns["count"](C, G3)["genome"] != 2:
            f.append("a table without a usable query column did not fall back")
        # a query column with no 'all' rows must REFUSE, not silently miscount
        G4 = _frame([("a", 0.001, "up"), ("b", 0.001, "down")])
        try:
            ns["count"](C, G4)
            f.append("a table with no 'all' rows did not refuse")
        except SystemExit:
            pass
    again, p2 = apply_edits(new)
    if not p2 or "already carries 166" not in p2[0]:
        f.append("re-application did not refuse")
    print("166 self-test")
    for x in f:
        print("  ✗ " + x)
    if not f:
        print("  ✓ on the real 125 bp shape (85 all + 69 up + 8 down + 17 multi")
        print("    = 179 rows) panel A now counts 85, and the proper background 1;")
        print("    alpha still applies within the 'all' list; a table with no query")
        print("    column falls back and warns; a query column with no 'all' rows")
        print("    refuses; the file-vs-FROZEN_TOTAL cross-check is intact;")
        print("    re-application refuses")
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
        print("\nREFUSING: 166's own self-test failed.")
        return 3
    print()
    if not os.path.exists(a.target):
        print(f"REFUSING: {a.target} not found")
        return 3
    src = open(a.target, errors="replace").read()
    new, probs = apply_edits(src)
    print("=" * 74)
    print("166 — panel A counts the 'all' query list only")
    print("=" * 74)
    if probs:
        print("REFUSING:")
        for p in probs:
            print("    " + p)
        return 3
    print("".join(difflib.unified_diff(src.splitlines(True), new.splitlines(True),
                                       fromfile=a.target,
                                       tofile=a.target + " (patched)", n=2)))
    if not a.apply:
        print("DRY RUN. Re-run with --apply, then:")
        print()
        print("  # 125 bp — panel A must now read 1 and 85")
        print("  F6_COUNTS=12564/805/550/255/73/131 F6_TOTALS=1/85 F6_MULTI=0/17 \\")
        print("  python3 41_figure6_enrichment.py \\")
        print("      --correct compadj_rl125_gprofiler.tsv \\")
        print("      --genome  compadj_rl125_gprofiler_WHOLEGENOME.tsv \\")
        print("      --out F6_enrichment_rl125")
        print()
        print("  The provenance block must end  'terms: 1 vs 85'  and NO ***")
        print("  warning should appear: the measured value and F6_TOTALS now")
        print("  describe the same quantity, so agreement is informative.")
        print()
        print("  # 50 bp control — it WILL now differ, and the *** warning WILL")
        print("  #   fire. FROZEN_TOTAL's 4 and 229 are four-list totals. Record")
        print("  #   what it prints; do NOT edit FROZEN_TOTAL to match.")
        print("  python3 41_figure6_enrichment.py \\")
        print("      --correct compadj_gprofiler.tsv \\")
        print("      --genome  compadj_gprofiler_WHOLEGENOME.tsv \\")
        print("      --out /tmp/f6_cal2")
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
    print("✓ verified: panel A is scoped to the 'all' query list and the")
    print("  file-vs-FROZEN_TOTAL cross-check is still in place.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
