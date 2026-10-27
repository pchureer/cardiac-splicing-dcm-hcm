#!/usr/bin/env python3
"""
162_patch_figures_37_38_readlength.py -- Figures 3 and 4 at 125 bp.

WHY THIS EXISTS
---------------
On 23 Oct both figures were regenerated with `--out F3_etiology_rl125` and
`--out F4_replication_rl125`. Both came out carrying **50 bp data under a
125 bp filename**: Figure 3 showed 926 / 468 / 62,971 / 48,024 / 45,380 and
Figure 4 showed 1,636/2,163 with slopes 0.48/0.52/0.39/0.24.

`--out` renames the OUTPUT. Neither script takes a suffix for its INPUTS.

That is Rule 10 inverted and worse than an overwrite: an overwrite destroys a
comparator, but this *asserts* a read length the content does not have, and
nothing in either script's output says so. It was caught only by looking at the
PNG. Both files have been deleted.

WHAT EACH FIGURE NEEDS
----------------------
**Figure 3** (`37_figure3_etiology.py`)
  panel A  reads the FROZEN `master_numbers_flagged.tsv`. There is no
           `master_numbers_flagged_rl125.tsv`, so panel A has no 125 bp source
           at all. `F3_PANEL_A` supplies the five counts directly rather than
           inventing a frozen file that nothing else would maintain.
  panel B  Colorado DCM-vs-ICM and Colorado DCM-vs-NF. Colorado was never
           re-quantified, so panel B is READ-LENGTH INVARIANT and must come out
           byte-identical. That is also the check that the patch did not
           over-reach.
  panel C  `{et}_v2filtered` and `{et}_HCM_v2`, both MAGNet. The first takes the
           suffix. The second was RENAMED at 125 bp to `SE_HCMvsNF_rl125`, so a
           suffix alone will not find it and `F3_HCM_BASE` names it.

**Figure 4** (`38_figure4_replication.py`)
  panel A  `--conc` is already a flag and `direction_concordance_dedup_rl125.tsv`
           exists. No patch needed: pass the flag.
  panel B  recomputes Deming slopes against MAGNet. Every entry of `REPL`
           (Colorado, Brigham, Fuwai, Stanford, LBNL) is a cohort that was NOT
           re-quantified and must stay unsuffixed; only the MAGNet discovery
           table moves.
  ceilings `CEIL_PEOPLE = 0.988` and `CEIL_PLATFORM = 0.584` are PLOTTED
           CONSTANTS, and their labels "people 0.99" / "platform 0.58" are
           hard-coded a second time at the draw call. This is exactly Figure 5's
           panel D, which carried 50 bp literals into a 125 bp figure while every
           gate passed. At 125 bp the matched panel gives **0.690 and 0.968**
           (RESULT_s4_rl125.md), so both the constants and the labels move.

USAGE
-----
  python3 162_patch_figures_37_38_readlength.py            # dry run
  python3 162_patch_figures_37_38_readlength.py --apply
  python3 162_patch_figures_37_38_readlength.py --self-test
"""
import argparse
import difflib
import os
import re
import shutil
import sys
import time

T37 = "37_figure3_etiology.py"
T38 = "38_figure4_replication.py"

# ----------------------------------------------------------------- figure 3
A37_SUF = '''def panel_a(ax, M, et, extra):'''
B37_SUF = '''# 162: the MAGNet arms re-base with the read length. Colorado does not: it was
# never re-quantified, which is why panel B must come out byte-identical.
F3SUF = os.environ.get("F3_SUFFIX", "")
# The 125 bp HCM-vs-NF table was RENAMED, not suffixed: SE_HCM_v2 became
# SE_HCMvsNF_rl125. Rule 8 -- a rename is only half an edit -- so it is named
# explicitly rather than derived.
F3_HCM = os.environ.get("F3_HCM_BASE", "")


def panel_a(ax, M, et, extra):'''

A37_PANELA = '''    sigcol = "sig_exon" if "sig_exon" in M.columns else "significant"
    tescol = "tested_exon" if "tested_exon" in M.columns else "tested"
    rows = []'''
B37_PANELA = '''    sigcol = "sig_exon" if "sig_exon" in M.columns else "significant"
    tescol = "tested_exon" if "tested_exon" in M.columns else "tested"
    # 162: panel A reads the frozen master table, which has no 125 bp twin.
    # F3_PANEL_A overrides the plotted counts rather than inventing a frozen
    # file nothing else would maintain. Every named row must already exist, so
    # a typo refuses instead of silently plotting the 50 bp value.
    #   F3_PANEL_A="Colorado/DCMvsNF=476/19953,Colorado/DCMvsICM=2/19998,..."
    _ov = os.environ.get("F3_PANEL_A", "").strip()
    if _ov:
        M = M.copy()
        for _item in [x for x in _ov.split(",") if x.strip()]:
            _k, _v = _item.split("=")
            _coh, _con = _k.split("/")
            _s, _t = (int(x) for x in _v.split("/"))
            _m = ((M.cohort == _coh) & (M.contrast == _con)
                  & (M.event_type == et))
            if not _m.any():
                raise SystemExit(
                    f"REFUSING: F3_PANEL_A names {_k} for {et}, which is not a "
                    f"row of the master table. Nothing plotted.")
            M.loc[_m, sigcol] = _s
            M.loc[_m, tescol] = _t
        print(f"  [162] panel A overridden from F3_PANEL_A: {_ov}")
    rows = []'''

A37_DISC = '''    D, _ = load_events(f"{et}_v2filtered", et)'''
B37_DISC = '''    D, _ = load_events(f"{et}_v2filtered" + F3SUF, et)'''

A37_HCM = '''    H, ph = load_events(f"{et}_HCM_v2", et)'''
B37_HCM = '''    H, ph = load_events(F3_HCM or f"{et}_HCM_v2", et)'''

# ----------------------------------------------------------------- figure 4
A38_CEIL = '''CEIL_PEOPLE = 0.988      # different people, same platform
CEIL_PLATFORM = 0.584    # same people, different platform'''
B38_CEIL = '''# 162: PLOTTED CONSTANTS. Figure 5's panel D was drawn from literals exactly
# like these and carried 50 bp values into a 125 bp figure while every gate
# passed. At 125 bp the matched panel gives 0.690 and 0.968
# (RESULT_s4_rl125.md).  F4_CEILINGS="platform/people", e.g. "0.690/0.968".
_C4 = os.environ.get("F4_CEILINGS", "").strip()
if _C4:
    CEIL_PLATFORM, CEIL_PEOPLE = (float(_x) for _x in _C4.split("/"))
else:
    CEIL_PEOPLE = 0.988      # different people, same platform
    CEIL_PLATFORM = 0.584    # same people, different platform
# 162: only the MAGNet discovery arm re-bases. Every REPL entry below is a
# cohort that was NOT re-quantified and must stay unsuffixed.
F4SUF = os.environ.get("F4_SUFFIX", "")'''

A38_LAB = '''    for xv, lab, col in ((CEIL_PLATFORM, "platform 0.58", ORANGE),
                         (CEIL_PEOPLE, "people 0.99", BLUE)):'''
B38_LAB = '''    # 162: the labels are derived from the constants. They were a SECOND
    # hard-coded copy of the same two numbers, so a re-based ceiling would have
    # been drawn with its old label.
    for xv, lab, col in ((CEIL_PLATFORM, f"platform {CEIL_PLATFORM:.2f}", ORANGE),
                         (CEIL_PEOPLE, f"people {CEIL_PEOPLE:.2f}", BLUE)):'''

EDITS = {
    T37: [("the suffix variables", A37_SUF, B37_SUF, 1),
          ("panel A's override", A37_PANELA, B37_PANELA, 1),
          ("the discovery table", A37_DISC, B37_DISC, 1),
          ("the HCM table (renamed at 125 bp)", A37_HCM, B37_HCM, 1)],
    T38: [("the plotted ceilings", A38_CEIL, B38_CEIL, 1),
          ("the ceiling labels", A38_LAB, B38_LAB, 1)],
}

# 38's MAGNet discovery base is not visible in the greps run so far, so it is
# located by pattern rather than by a fixed anchor, and the patch REFUSES
# rather than guessing if the pattern does not appear exactly once.
RE38_DISC = re.compile(r'f"\{et\}_v2filtered"')


def apply_edits(path, text):
    if "162:" in text:
        return text, [f"{path} already carries 162 (restore the .bak to re-apply)"]
    if "import os" not in text:
        return text, [f"{path} does not import os at module level"]
    problems, out = [], text
    for name, anchor, repl, want in EDITS[path]:
        got = out.count(anchor)
        if got != want:
            problems.append(f"{path}: {name}: anchor occurs {got}x, expected "
                            f"{want}x\n      first line: "
                            f"{anchor.splitlines()[0][:66]}")
            continue
        out = out.replace(anchor, repl, want)
    if path == T38:
        n = len(RE38_DISC.findall(out))
        if n != 1:
            problems.append(
                f"{T38}: found {n} occurrences of f\"{{et}}_v2filtered\", "
                f"expected exactly 1. Panel B's MAGNet discovery table could "
                f"not be located, so the slopes would stay at 50 bp while the "
                f"ceilings moved. Send `grep -n v2filtered {T38}` and 162 will "
                f"be given the real anchor.")
        else:
            out = RE38_DISC.sub('(f"{et}_v2filtered" + F4SUF)', out)
    return out, problems


def verify(path, text):
    bad = []
    need = {T37: ['F3SUF = os.environ.get("F3_SUFFIX", "")',
                  'F3_HCM = os.environ.get("F3_HCM_BASE", "")',
                  'f"{et}_v2filtered" + F3SUF',
                  'F3_HCM or f"{et}_HCM_v2"',
                  "F3_PANEL_A", "row of the master table"],
            T38: ['_C4 = os.environ.get("F4_CEILINGS", "").strip()',
                  'F4SUF = os.environ.get("F4_SUFFIX", "")',
                  'f"platform {CEIL_PLATFORM:.2f}"',
                  'f"people {CEIL_PEOPLE:.2f}"',
                  '(f"{et}_v2filtered" + F4SUF)']}[path]
    for w in need:
        if w not in text:
            bad.append(f"missing: {w[:58]}")
    # the 50 bp defaults must survive
    if path == T38:
        for w in ("CEIL_PEOPLE = 0.988", "CEIL_PLATFORM = 0.584"):
            if w not in text:
                bad.append(f"a 50 bp default was lost: {w}")
        if '"platform 0.58"' in text or '"people 0.99"' in text:
            bad.append("a hard-coded ceiling LABEL survived; a re-based ceiling "
                       "would be drawn with its old number")
        # Colorado and the HCM cohorts must NOT be suffixed
        for w in ('("Colorado", "{ET}_Colorado_DCMvsNF")',
                  '("Brigham", "Brigham_{ET}_HCM")',
                  '("Fuwai", "Fuwai_{ET}_HCM")',
                  '("Stanford", "Stanford_{ET}_HCM")'):
            if w not in text:
                bad.append(f"a REPL entry was altered: {w[:46]} -- those "
                           f"cohorts were never re-quantified")
    if path == T37:
        # Colorado's two panel-B arms must stay untouched
        for w in ('load_events(f"{et}_Colorado_DCMvsICM", et)',
                  'load_events(f"{et}_Colorado_DCMvsNF", et)'):
            if w not in text:
                bad.append(f"Colorado was altered: {w[:52]} -- panel B is "
                           f"read-length invariant and is the control")
    return bad


# ------------------------------------------------------------------ self-test
F37 = '''#!/usr/bin/env python3
import argparse
import os
import pandas as pd


def load_events(base, et, why=None):
    return None, base


def panel_a(ax, M, et, extra):
    """Counts: disease-vs-control against disease-vs-other-etiology."""
    sigcol = "sig_exon" if "sig_exon" in M.columns else "significant"
    tescol = "tested_exon" if "tested_exon" in M.columns else "tested"
    rows = []
    return rows


def panel_b(ax, et, a):
    E, pe = load_events(f"{et}_Colorado_DCMvsICM", et)
    D, pd_ = load_events(f"{et}_Colorado_DCMvsNF", et)


def panel_c(ax, et, a):
    D, _ = load_events(f"{et}_v2filtered", et)
    H, ph = load_events(f"{et}_HCM_v2", et)
'''

F38 = '''#!/usr/bin/env python3
import argparse
import os
import pandas as pd

ORANGE = "#eb6834"
BLUE = "#2a78d6"
NKEEP = {"SE": 3}
# measured on the 24 MAGNet<->LBNL same-individual pairs, common event set
CEIL_PEOPLE = 0.988      # different people, same platform
CEIL_PLATFORM = 0.584    # same people, different platform
REPL = [("Colorado", "{ET}_Colorado_DCMvsNF"), ("Brigham", "Brigham_{ET}_HCM"),
        ("Fuwai", "Fuwai_{ET}_HCM"), ("Stanford", "Stanford_{ET}_HCM"),
        ("LBNL", "{ET}_LBNL_DCMvsNF")]


def panel_b(ax, et, a):
    M, _ = load(f"{et}_v2filtered", et)
    for xv, lab, col in ((CEIL_PLATFORM, "platform 0.58", ORANGE),
                         (CEIL_PEOPLE, "people 0.99", BLUE)):
        ax.axvline(xv, color=col)
'''


def self_test():
    fails = []
    for path, fx in ((T37, F37), (T38, F38)):
        new, probs = apply_edits(path, fx)
        if probs:
            fails.append(f"{path} fixture did not patch: {probs}")
            continue
        bad = verify(path, new)
        if bad:
            fails.append(f"{path} verify rejected its own output: {bad}")
        try:
            compile(new, path, "exec")
        except SyntaxError as e:
            fails.append(f"{path} patched fixture does not compile: {e}")
        again, p2 = apply_edits(path, new)
        if not p2 or "already carries 162" not in p2[0]:
            fails.append(f"{path}: re-application did not refuse")

    # 38: the ceiling LABELS must follow the constants
    n38, _ = apply_edits(T38, F38)
    ns = {}
    exec(compile(n38.replace("def panel_b", "def _unused_panel_b"), "x", "exec"), ns)
    if (ns["CEIL_PLATFORM"], ns["CEIL_PEOPLE"]) != (0.584, 0.988):
        fails.append(f"50 bp defaults wrong: {ns['CEIL_PLATFORM']}, {ns['CEIL_PEOPLE']}")
    os.environ["F4_CEILINGS"] = "0.690/0.968"
    ns2 = {}
    exec(compile(n38.replace("def panel_b", "def _unused_panel_b"), "x", "exec"), ns2)
    if (ns2["CEIL_PLATFORM"], ns2["CEIL_PEOPLE"]) != (0.690, 0.968):
        fails.append(f"F4_CEILINGS not honoured: {ns2.get('CEIL_PLATFORM')}, "
                     f"{ns2.get('CEIL_PEOPLE')}")
    del os.environ["F4_CEILINGS"]
    if 'f"platform {CEIL_PLATFORM:.2f}"' not in n38:
        fails.append("the ceiling label is not derived from the constant")

    # 38: REPL must be untouched and the discovery base suffixed exactly once
    if n38.count('+ F4SUF') != 1:
        fails.append(f"F4SUF applied {n38.count('+ F4SUF')}x, want 1")
    for coh in ("Colorado", "Brigham", "Fuwai", "Stanford", "LBNL"):
        if f'("{coh}", ' not in n38:
            fails.append(f"REPL lost {coh}")
    if "F4SUF" in n38.split("REPL = [")[1].split("]")[0]:
        fails.append("a REPL entry was suffixed; none of those cohorts moved")

    # 38: refuse when the discovery base cannot be located
    _, p = apply_edits(T38, F38.replace('f"{et}_v2filtered"', 'SOMETHING_ELSE'))
    if not any("could not be located" in x for x in p):
        fails.append("a missing discovery base did not refuse")

    # 37: panel A override applies and a bad key refuses
    n37, _ = apply_edits(T37, F37)
    import pandas as pd
    ns3 = {"os": os, "pd": pd}
    exec(compile(n37, "x", "exec"), ns3)
    M = pd.DataFrame({"cohort": ["MAGNet"], "contrast": ["DCMvsNF"],
                      "event_type": ["SE"], "significant": [926],
                      "tested": [62971]})
    os.environ["F3_PANEL_A"] = "MAGNet/DCMvsNF=813/53456"
    ns3["panel_a"](None, M, "SE", None)
    os.environ["F3_PANEL_A"] = "MAGNet/NOSUCH=1/2"
    try:
        ns3["panel_a"](None, M, "SE", None)
        fails.append("an unknown F3_PANEL_A row did not refuse")
    except SystemExit as e:
        if "not a row of the master table" not in str(e):
            fails.append(f"wrong refusal message: {e}")
    del os.environ["F3_PANEL_A"]
    if 'load_events(f"{et}_Colorado_DCMvsICM", et)' not in n37:
        fails.append("37: Colorado's panel B arm was altered")

    print("162 self-test")
    for f in fails:
        print("  ✗ " + f)
    if not fails:
        print("  ✓ both files patch and compile, 50 bp defaults survive,")
        print("    F4_CEILINGS moves BOTH the constants and their labels,")
        print("    REPL and Colorado untouched, the discovery base suffixed")
        print("    exactly once, a missing base refuses, and F3_PANEL_A")
        print("    overrides a real row and refuses an unknown one")
    return 1 if fails else 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--apply", action="store_true")
    ap.add_argument("--self-test", action="store_true")
    a = ap.parse_args()
    if a.self_test:
        return self_test()
    if self_test():
        print("REFUSING: 162's own self-test failed.")
        return 3
    print()

    rc = 0
    staged = {}
    for path in (T37, T38):
        if not os.path.exists(path):
            print(f"REFUSING: {path} not found -- run from the HEART directory")
            return 3
        src = open(path, errors="replace").read()
        new, probs = apply_edits(path, src)
        print("=" * 78)
        print(f"162 -- {path}")
        print("=" * 78)
        if probs:
            print("REFUSING -- anchors did not match:")
            for p in probs:
                print("    " + p)
            rc = 3
            continue
        print("".join(difflib.unified_diff(
            src.splitlines(True), new.splitlines(True),
            fromfile=path, tofile=path + " (patched)", n=2)))
        staged[path] = (src, new)
    if rc:
        return rc

    if not a.apply:
        print("=" * 78)
        print("DRY RUN -- nothing written. Re-run with --apply, then:")
        print()
        print("  # 1. 50 bp controls. BOTH must be byte-identical to the")
        print("  #    published figures; 89_make_all_figures.sh checks md5.")
        print("  python3 37_figure3_etiology.py --event-type SE --out /tmp/f3_cal")
        print("  python3 38_figure4_replication.py --out /tmp/f4_cal")
        print()
        print("  # 2. Figure 3 at 125 bp")
        print('  F3_SUFFIX=_rl125 F3_HCM_BASE=SE_HCMvsNF_rl125 \\')
        print('  F3_PANEL_A="Colorado/DCMvsNF=476/19953,'
              'Colorado/DCMvsICM=2/19998,\\')
        print('MAGNet/DCMvsNF=813/53456,MAGNet/HCMvsNF=431/39607,'
              'MAGNet/DCMvsHCM=0/38273" \\')
        print("  python3 37_figure3_etiology.py --event-type SE \\")
        print("      --out F3_etiology_rl125")
        print()
        print("  # 3. Figure 4 at 125 bp")
        print("  F4_SUFFIX=_rl125 F4_CEILINGS=0.690/0.968 \\")
        print("  python3 38_figure4_replication.py \\")
        print("      --conc direction_concordance_dedup_rl125.tsv \\")
        print("      --out F4_replication_rl125")
        print()
        print("⚠ THEN LOOK AT THE PNGs. That is what caught this today.")
        print("  Figure 3 must show 813 / 431 / 0 of 38,273 / 53,456 / 39,607,")
        print("  panel C n = 658, Deming 0.82, 92% -- and panel B MUST be")
        print("  UNCHANGED (95%, 0.16), because both its arms are Colorado.")
        print("  Figure 4 must show 1,376/1,800 = 76.4%, slopes 0.60 / 0.51 /")
        print("  0.40 / 0.29, LBNL 122/156, ceilings 0.69 and 0.97.")
        print()
        print("⚠ If panel B of Figure 3 moves, the patch over-reached and the")
        print("  run must be discarded: Colorado was never re-quantified.")
        return 0

    for path, (src, new) in staged.items():
        bak = f"{path}.bak.{time.strftime('%Y%m%d_%H%M%S')}"
        shutil.copy2(path, bak)
        open(path, "w").write(new)
        print(f"[write] {path}   (backup: {bak})")
        disk = open(path, errors="replace").read()
        bad = verify(path, disk)
        try:
            compile(disk, path, "exec")
        except SyntaxError as e:
            bad.append(f"does not compile: {e}")
        if bad or disk != new:
            print("✗ ON-DISK VERIFICATION FAILED -- restoring:")
            for b in bad:
                print(f"    {b}")
            shutil.copy2(bak, path)
            return 3
    print("✓ verified on disk: the MAGNet arms follow their suffix, Colorado")
    print("  and every replication cohort are untouched, Figure 4's ceilings")
    print("  and their labels move together, and panel A refuses an unknown row.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
