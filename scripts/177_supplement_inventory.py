#!/usr/bin/env python3
"""
177_supplement_inventory.py -- is the SUPPLEMENT finished at 125 bp?

The main figures (F1-F6) and tables (1-4) were verified on 26 Oct: every file
regenerates from its documented command. The supplement was never put through
that check. This script answers the narrower question first: for each of the
twelve supplementary items (10 tables, 2 figures) plus the Supplementary
Methods and Discussion, DOES A 125 bp FILE EXIST, when was it written, and
does the caption text still carry 50 bp numbers?

It reads. It writes nothing. Run it in the HEART directory:

  python 177_supplement_inventory.py
  python 177_supplement_inventory.py --captions supplementary_captions.md \
         --captions Supplementary_Methods.md --captions Supplementary_Discussion.md

Verdicts per item
  OK           a 125 bp file exists (or the item is read-length independent)
  50bp-ONLY    only a file without the _rl125 suffix exists -> needs rebuild
  MISSING      no file matches any pattern -> producer never wrote one
  CHECK        item exists but the record says something specific to confirm

Exit 1 if anything is 50bp-ONLY or MISSING, so it can gate the packaging step.
"""

import argparse
import glob
import os
import re
import signal
import sys
import time

try:
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)   # `| head` must not traceback
except (AttributeError, ValueError):
    pass

# Each item: (label, description, independent?, rl125_globs, any_globs, note)
# independent=True means the content does not depend on the rMATS read-length
# declaration (expression-level or genotype-level), so a 50 bp-dated file is
# the right file.
ITEMS = [
    ("Supp Table 1", "composition and covariate-model robustness (76_build_s1_v2.py)",
     False, ["S1a*rl125*", "*S1a_rl125*", "S1_rl125*"], ["S1a*", "S1_*"],
     "16 Oct record: measured at 125 bp but 76's fourth hard-coded gate "
     "(926/231/90/84/114) refused before writing. Confirm the file exists."),
    ("Supp Table 2", "the covariate ladder (69_build_s1b_ladder.py)",
     False, ["S1b*rl125*"], ["S1b*"],
     "16 Oct record: 69 printed 50 bp values (926/899/835/866/888). The 125 bp "
     "ladder is 813/777/732/744/770; check the file carries THOSE."),
    ("Supp Table 3", "marker panels with cell-type source (70_marker_panel_audit.py)",
     True, [], ["S1c_markers*", "marker_panel_audit*"],
     "Expression-level; invariant (confirmed 23 Oct). The 50 bp file is correct."),
    ("Supp Table 4", "all five event classes, all contrasts (81_build_s2.py)",
     False, ["S2*rl125*"], ["S2.tsv", "S2_*"],
     "Counts measured 15 Oct (RESULT_etiology_allclasses_measured); confirm 81 "
     "was re-run with the suffix."),
    ("Supp Table 5", "quantification choices (a)-(d) (62/63)",
     False, ["S5_rl125*"], ["S5.tsv", "S5_*"],
     "Complete at both read lengths 21 Oct, all self-checks passing."),
    ("Supp Table 6", "control-group size sweep (83_build_s6.py)",
     False, ["S6*rl125*"], ["S6*"],
     "No 125 bp record exists in the project. Likely still 50 bp."),
    ("Supp Table 7", "candidate regulators (65 / 23 / 58)",
     False, ["signature_rbp_margins_rl125*", "rbm20_margin_resolve_rl125*"],
     ["signature_rbp_margins*", "rbm20_margin_resolve*"],
     "Files exist (F5 reads them). CAPTION: leader is PTBP1 +0.062 SE / +0.113 "
     "MXE at 125 bp, not QKI +0.033 -- see --captions scan."),
    ("Supp Table 8", "null behaviour and reverse replication (74 / 22 / Colorado MXE)",
     False, ["S8a*rl125*", "symmetry_*rl125*", "*colorado*mxe*rl125*",
             "*reverse*rl125*"],
     ["S8a*", "symmetry_SE*", "symmetry_MXE*"],
     "(b) symmetry re-run 20 Oct; (c) Colorado MXE reverse slope re-measured "
     "22 Oct; (a) permutations: confirm 74 ran at 125 bp (DEFECT_03_n_perm)."),
    ("Supp Table 9", "genotype identity: pair structure, the 27 excluded (27-31)",
     True, [], ["S9*", "*duplicate_pairs*", "*excluded*"],
     "Genotype-level; read-length independent."),
    ("Supp Table 10", "LVEF null, summary (77_lvef_s3.py)",
     False, ["S3*rl125*", "*lvef*rl125*"], ["S3_lvef*", "*lvef*"],
     "RESULT_s3_lvef_rl125 (16 Oct) measured it; confirm the summary file."),
    ("Supp Figure 1", "genotype identity figure (kinship, IBS0)",
     True, [], ["S9*.png", "S9*.pdf", "*kinship*.png", "*kinship*.pdf"],
     "Genotype-level; read-length independent."),
    ("Supp Figure 2", "platform versus sampling, 24 matched individuals (82 / 36)",
     False, ["*matched*rl125*.png", "*matched*rl125*.pdf", "S4*rl125*.png",
             "S4*rl125*.pdf", "mp_panel_rl125*"],
     ["S4*.png", "S4*.pdf", "*matched*.png", "*matched*.pdf", "mp_panel*"],
     "22 Oct: 36 writes NO figure file; the 125 bp panel exists only as "
     "mp_panel_rl125_20261022.out. 82's figure must be regenerated from the "
     "mp_SE_*_rl125 tables -- confirm a 125 bp image exists."),
    ("Supp Methods", "M1-M10, rebuilt 23 Oct at 125 bp (4,121 words)",
     False, ["Supplementary_Methods*rl125*", "SUPP_METHODS*rl125*"],
     ["Supplementary_Methods*", "supp_methods*"],
     "The rebuilt text lives in the project doc SUPP_METHODS_rebuilt_rl125.md. "
     "A 1,030-word file holding only S1-S4 is the OLD one."),
    ("Supp Discussion", "Britto-Borges reconciliation, rebuilt 23 Oct",
     False, ["Supplementary_Discussion*rl125*", "SUPP_DISCUSSION*rl125*"],
     ["Supplementary_Discussion*"],
     "Rebuilt text is SUPP_DISCUSSION_rebuilt_rl125.md; 1,355 -> 1,430."),
]

# 50 bp values that must not appear in a 125 bp caption unless explicitly
# labelled as the 50 bp comparator. Reported with context; you judge.
STALE = [
    (r"\b926\b", "50 bp SE count (now 813)"),
    (r"\b1,?054\b", "50 bp restricted count (now 939)"),
    (r"\b62,?971\b", "50 bp SE universe (now 53,456)"),
    (r"\b0\.584\b", "50 bp platform slope (now 0.690)"),
    (r"\b0\.988\b", "50 bp people slope (now 0.968)"),
    (r"\+?0\.033\b", "50 bp QKI margin (leader now PTBP1 +0.062/+0.113)"),
    (r"\bQKI\b", "QKI led the scan at 50 bp; PTBP1 leads at 125 bp"),
    (r"three of five", "pattern claim; now FOUR of five (RI crossed over)"),
    (r"at least 99", "false at 125 bp (RI 97.8%); write 'at least 98%'"),
    (r"\b0\.997\b", "pilot rho; retired for 0.996 on 332 samples"),
    (r"\b1,?355\b", "Supp Discussion count; now 1,430"),
    (r"\b75\.6%", "50 bp replication (now 76.4%)"),
    (r"\b2,?163\b", "50 bp measurable events (now 1,800)"),
    (r"\b1,?636\b", "50 bp concordant (now 1,376)"),
    (r"\b86%", "50 bp median retention (now 84%)"),
    (r"Every interval covers", "false at 125 bp: PTBP1 MXE excludes the floor"),
    (r"\b899\b", "50 bp nine-marker count (now 777)"),
    (r"\b1\.124\b|\b1\.130\b|\b1\.214\b|\b1\.215\b", "50 bp SE inflation values"),
]


def rows_in(path):
    try:
        with open(path, "rb") as fh:
            return sum(1 for _ in fh)
    except OSError:
        return None


def find(patterns, root):
    out = set()
    for p in patterns:
        out.update(glob.glob(os.path.join(root, p)))
        out.update(glob.glob(os.path.join(root, "**", p), recursive=True))
    # drop backups and verify dirs
    return sorted(f for f in out
                  if "verify_backup" not in f and ".bak" not in f
                  and os.path.isfile(f))


def fmt(path):
    st = os.stat(path)
    when = time.strftime("%Y-%m-%d %H:%M", time.localtime(st.st_mtime))
    n = rows_in(path) if path.endswith((".tsv", ".csv", ".txt", ".md")) else None
    size = f"{st.st_size/1024:.0f} KB"
    return f"{when}  {size:>8}  {('%d rows' % n) if n else '':>10}  {os.path.relpath(path)}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=".")
    ap.add_argument("--captions", action="append", default=[],
                    help="caption / supplement text file to scan for 50 bp "
                         "numbers (repeatable)")
    ap.add_argument("--max-files", type=int, default=6)
    a = ap.parse_args()

    bad = 0
    print(f"supplement inventory under {os.path.abspath(a.root)}\n")
    for label, desc, indep, rl, anyg, note in ITEMS:
        rl_files = find(rl, a.root) if rl else []
        any_files = find(anyg, a.root)
        if indep:
            verdict = "OK (read-length independent)" if any_files else "MISSING"
        elif rl_files:
            verdict = "OK"
        elif any_files:
            verdict = "50bp-ONLY"
        else:
            verdict = "MISSING"
        if verdict in ("50bp-ONLY", "MISSING"):
            bad += 1
        flag = {"OK": "✅", "OK (read-length independent)": "✅",
                "50bp-ONLY": "⚠", "MISSING": "❌"}[verdict]
        print(f"{flag} {label:16s} {verdict:14s} {desc}")
        shown = rl_files or any_files
        for f in shown[:a.max_files]:
            print(f"      {fmt(f)}")
        if len(shown) > a.max_files:
            print(f"      … {len(shown) - a.max_files} more")
        if note:
            print(f"      note: {note}")
        print()

    if a.captions:
        print("=" * 70)
        print("50 bp numbers in caption / supplement text (judge each in context)")
        for cf in a.captions:
            if not os.path.isfile(cf):
                print(f"\n  {cf}: not found")
                continue
            print(f"\n  {cf}")
            hits = 0
            with open(cf, encoding="utf-8", errors="replace") as fh:
                for i, line in enumerate(fh, 1):
                    if line.lstrip().startswith("<!--"):
                        continue
                    for pat, why in STALE:
                        if re.search(pat, line):
                            hits += 1
                            print(f"    L{i:<5d} [{why}]")
                            print(f"           {line.strip()[:140]}")
                            break
            if not hits:
                print("    no stale tokens found")
    else:
        print("(pass --captions FILE to scan caption text for 50 bp numbers)")

    print("\n" + "=" * 70)
    if bad:
        print(f"{bad} item(s) are 50bp-ONLY or MISSING. The supplement is NOT "
              f"finished at 125 bp.")
        sys.exit(1)
    print("Every supplementary item has a 125 bp file (or is read-length "
          "independent). Next: caption scan, then assemble the single PDF.")


if __name__ == "__main__":
    main()
