#!/usr/bin/env python3
"""
178_patch_81_s2_rl125.py -- make 81_build_s2.py (Supplementary Table 4, the
all-classes/all-contrasts table) read-length aware.

81 hard-codes the 50 bp MAGNet table names (lines 74-76) and the 50 bp record it
gates against (G_SIG, G_TESTED, G_TESTED_TOTAL, G_OTHER). Run on 125 bp data it
would either not find the tables or refuse on every gate. The six
non-requantified cohorts (Colorado x3, Brigham, Fuwai, Stanford) keep their
names and their record: those files did not change.

After the patch 81 behaves exactly as before with TBL_SUFFIX unset, and with
TBL_SUFFIX=_rl125 it:
  * reads {ET}_v2filtered_rl125, {ET}_HCMvsNF_rl125, {ET}_DCMvsHCM_rl125
  * writes S2_rl125.tsv
  * gates against the 125 bp record (813/229/87/74/91 on 53,456/62,537/14,346/
    8,360/6,243 = 144,942; HCM vs NF SE 431/39,607; DCM vs HCM SE 0/38,273)
  * REFUSES unless S2_FUWAI_R is set, because the Fuwai MXE exclusion text
    quotes a correlation that has to be read from symmetry_MXE_rl125.tsv
  * refuses any other suffix (no record for it)

  python 178_patch_81_s2_rl125.py            # dry run: shows what would change
  python 178_patch_81_s2_rl125.py --apply    # backs up, patches, self-checks

Then calibrate before believing anything:
  python 81_build_s2.py --out S2_50check.tsv && cmp S2.tsv S2_50check.tsv
and only if IDENTICAL:
  TBL_SUFFIX=_rl125 S2_FUWAI_R=0.0xx python 81_build_s2.py
"""

import argparse
import re
import shutil
import sys
import time
from pathlib import Path

TARGET = Path("81_build_s2.py")
MARK = "# [178] read-length aware"

EDITS = [
    # 1. suffix + output name
    ('OUT = "S2.tsv"\n',
     'SUF = os.environ.get("TBL_SUFFIX", "")   ' + MARK + '\n'
     'OUT = f"S2{SUF}.tsv"\n'),
    # 2. the three MAGNet rows
    ('    ("MAGNet",   "DCM vs NF",  "{ET}_v2filtered.all_events.tsv",               True),\n'
     '    ("MAGNet",   "HCM vs NF",  "{ET}_HCM_v2.all_events.tsv",                   False),\n'
     '    ("MAGNet",   "DCM vs HCM", "{ET}_DCMvsHCM.all_events.tsv",                 False),\n',
     '    ("MAGNet",   "DCM vs NF",  "{ET}_v2filtered" + SUF + ".all_events.tsv",   True),\n'
     '    ("MAGNet",   "HCM vs NF",  ("{ET}_HCMvsNF" + SUF if SUF else "{ET}_HCM_v2")\n'
     '                               + ".all_events.tsv",                           False),\n'
     '    ("MAGNet",   "DCM vs HCM", "{ET}_DCMvsHCM" + SUF + ".all_events.tsv",     False),\n'),
    # 3. the record, re-based under the suffix
    ('REC_EQUIV_FULLCOORD = {"SE": 0.220, "MXE": 0.190, "A3SS": 0.238,\n'
     '                       "A5SS": 0.204, "RI": 0.308}\n',
     'REC_EQUIV_FULLCOORD = {"SE": 0.220, "MXE": 0.190, "A3SS": 0.238,\n'
     '                       "A5SS": 0.204, "RI": 0.308}\n'
     '\n'
     + MARK + '. The block above is the 50 bp record and stays the default.\n'
     '# Under a suffix the MAGNet rows are re-based from the 125 bp record; the six\n'
     '# non-requantified cohorts keep their values because their files did not change.\n'
     'if SUF == "_rl125":\n'
     '    G_SIG = {"SE": 813, "MXE": 229, "A3SS": 87, "A5SS": 74, "RI": 91}\n'
     '    G_TESTED = {"SE": 53_456, "MXE": 62_537, "A3SS": 14_346, "A5SS": 8_360,\n'
     '                "RI": 6_243}\n'
     '    G_TESTED_TOTAL = 144_942\n'
     '    G_OTHER[("MAGNet", "DCM vs HCM", "SE")] = (0, 38_273)\n'
     '    G_OTHER[("MAGNet", "HCM vs NF", "SE")] = (431, 39_607)\n'
     '    _fuwai_r = os.environ.get("S2_FUWAI_R", "").strip()\n'
     '    if not _fuwai_r:\n'
     '        sys.exit("[REFUSE] TBL_SUFFIX=_rl125 but S2_FUWAI_R is unset. The Fuwai MXE\\n"\n'
     '                 "exclusion text quotes a correlation that must be read, not assumed:\\n"\n'
     '                 "  awk -F\'\\\\t\' \'$1==\\"Fuwai\\"{print $16, $21}\' symmetry_MXE_rl125.tsv\\n"\n'
     '                 "(unsel_r and prec_r; the 50 bp values were 0.088 / 0.111), then\\n"\n'
     '                 "  TBL_SUFFIX=_rl125 S2_FUWAI_R=<unsel_r> python 81_build_s2.py")\n'
     '    MXE_EXCLUDED = {\n'
     '        "Colorado": "reverse-replication slope 0.060 (n = 1,478) against SE\'s 0.443",\n'
     '        "Fuwai": f"lowest MXE correlation in the study, r = {_fuwai_r} (paired-end 150 bp)",\n'
     '    }\n'
     'elif SUF:\n'
     '    sys.exit(f"[REFUSE] TBL_SUFFIX={SUF!r}: no record for this suffix. Add one to the\\n"\n'
     '             f"block marked [178] before building under it.")\n'),
]

# 50 bp literals that must appear ONLY inside the record block after the patch.
LITERALS = [r"\b926\b", r"62_?971", r"156_?115", r"45_?380", r"48_?024", r"\b468\b",
            r"0\.061", r"1,477", r"0\.447", r"0\.088"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--apply", action="store_true")
    ap.add_argument("--target", default=str(TARGET))
    a = ap.parse_args()
    p = Path(a.target)
    if not p.is_file():
        sys.exit(f"[fatal] {p} not found; run in the HEART directory")
    src = p.read_text()

    if MARK in src:
        print(f"{p}: already patched ({MARK!r} present). Nothing to do.")
        return
    if "import os" not in src:
        sys.exit("[fatal] 81 does not import os (it should: load() uses os.path). "
                 "Check the file before patching.")

    new = src
    for i, (old, rep) in enumerate(EDITS, 1):
        n = new.count(old)
        if n != 1:
            sys.exit(f"[fatal] edit {i}: anchor found {n} times, expected 1. The file "
                     f"differs from the 26 Oct listing; nothing written.\n--- anchor ---\n{old}")
        new = new.replace(old, rep, 1)
        print(f"  edit {i}: ok")

    # Where do 50 bp literals still live? Only the record block is allowed.
    rec_start = new.index("# MXE is excluded for these cohorts")
    _tail = 'block marked [178] before building under it.")\n'
    rec_end = new.index(_tail) + len(_tail)
    body = new[:rec_start] + new[rec_end:]
    stray = []
    for pat in LITERALS:
        for m in re.finditer(pat, body):
            line = body.count("\n", 0, m.start()) + 1
            ctx = body[max(0, m.start() - 40):m.end() + 40].replace("\n", " ")
            stray.append((pat, ctx))
    print()
    if stray:
        print("⚠ 50 bp literals OUTSIDE the record block (docstring mentions are fine; a\n"
              "  literal inside gate code is not -- read each):")
        for pat, ctx in stray:
            print(f"    {pat:10s} …{ctx}…")
    else:
        print("  no 50 bp literal outside the record block")

    if not a.apply:
        print("\nDRY RUN. Re-run with --apply to write.")
        return
    bak = p.with_suffix(p.suffix + ".bak." + time.strftime("%Y%m%d_%H%M%S"))
    shutil.copy2(p, bak)
    p.write_text(new)
    try:
        compile(new, str(p), "exec")
    except SyntaxError as e:
        shutil.copy2(bak, p)
        sys.exit(f"[fatal] patched file does not compile ({e}); original restored from {bak}")
    print(f"\npatched {p}; backup {bak}")
    print("\nNEXT, in this order:\n"
          "  python 81_build_s2.py --out S2_50check.tsv && cmp S2.tsv S2_50check.tsv\n"
          "      (must print nothing: the patch changes nothing at 50 bp)\n"
          "  awk -F'\\t' '$1==\"Fuwai\"{print $16, $21}' symmetry_MXE_rl125.tsv\n"
          "  TBL_SUFFIX=_rl125 S2_FUWAI_R=<unsel_r> python 81_build_s2.py\n"
          "      -> S2_rl125.tsv")


if __name__ == "__main__":
    main()
