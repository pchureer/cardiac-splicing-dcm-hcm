#!/usr/bin/env python3
"""
97_drop_se_genes_column.py — actually remove Table 3's "SE genes" column.

WHY THIS EXISTS
---------------
`RESULT_table3_built.md` states, under "A COLUMN DELIBERATELY REMOVED":

    v2 carried "SE genes" -- 671, 341, 320, 274 and 1 ... **Dropped** ...
    The column is commented out in the script rather than deleted.

**It is not.** The 12 Oct run of `87_build_table3.py v2` printed:

    Cohort   Contrast  ... +|dPSI|>=0.10 SE genes Classes ...
    MAGNet  DCM vs NF  ...           926      671       5 ...

So the doc recorded an intention as an accomplishment, and re-running 87 could never
drop the column because nothing in 87 had changed. Re-running a script does not apply
an edit that was never made.

    ⚠ A decision written in a doc is not a change made to a file. Record the edit and
      the verification separately, or the record will claim work that does not exist.

WHY THE COLUMN GOES (the reasoning is sound; only the application was missing)
  1. No claim in the paper depends on a gene count. The Results quote gene SETS --
     134 signature genes, 20 lead genes, 40 of 45 RBM20 targets, 130 of 134 -- and
     never a count of genes per contrast.
  2. 671 would sit beside the Methods' 134-gene signature set and invite a question the
     table cannot answer. "Genes with >=1 significant skipped exon" and "genes
     significant across multiple event classes" are different quantities and a reader
     cannot tell that from the table.

USAGE
    python 97_drop_se_genes_column.py --list     # show what it would change
    python 97_drop_se_genes_column.py --apply    # comment the column out, with a backup
    python 87_build_table3.py                    # then re-run 87 to rebuild Table3.tsv

It COMMENTS OUT rather than deletes, so the column can be restored if a sentence later
needs it -- which is what the record says was done.

GATES
  G1  87_build_table3.py exists and is readable.
  G2  the column is still present -- if it is already commented out, this script says so
      and changes nothing.
  G3  every line that mentions the column is shown, and the edit touches ONLY lines that
      bind it to the table. A mention inside a comment or a docstring is NOT edited.
  G4  the file still parses as Python after the edit (ast.parse). A syntax error here
      would break the only script that builds Table 3.
  G5  after writing, reload and confirm: the column's binding lines are commented, the
      line count is unchanged, and no other line differs.
"""

import argparse
import ast
import os
import re
import shutil
import sys
import time

SCRIPT = "87_build_table3.py"

# The header label as 87 prints it. Matching is on the quoted string, not on the
# rendered table, so a change in column order cannot hide it.
LABEL_PATTERNS = [
    re.compile(r'["\']SE[ _]genes["\']', re.I),
    re.compile(r'\bse_genes\b', re.I),
    re.compile(r'\bn_genes_se\b', re.I),
    re.compile(r'\bgenes_se\b', re.I),
]

# ⚠ TWO SHAPES, TWO DIFFERENT EDITS. v1 of this script tried to comment out whichever
# line mentioned the column, and got both cases wrong on the first test:
#
#   * `COLS = [..., "SE genes", "Classes"]` spread over two lines -- commenting the
#     inner line breaks the list literal. The right edit is to delete the ELEMENT.
#   * `rec["SE genes"] = r.get("n_genes_se", "-")` IS a binding, but v1's assignment
#     pattern excluded it because the subscript contains a space. It was left in place
#     while the unsafe line was chosen. Exactly inverted.
#
# So: an occurrence inside a list/tuple of strings has its element removed; a whole
# single-line statement that assigns to or from it is commented out; anything else is
# reported and left alone.
# ⚠ A THIRD SHAPE, found on the real file: a DICT ENTRY.
#
#     374            "SE genes": gl,
#
# Treating that as an element removes only the key and leaves `: gl,` — a syntax error,
# which G4 caught before anything was written. A dict entry has to be removed WHOLE,
# key and value together. Order matters below: test for the entry before the element,
# or the element pattern claims the key first.
DICT_ENTRY = re.compile(r'\s*["\']SE[ _]genes["\']\s*:\s*[^,}\n]+,?')
ELEMENT = re.compile(r'\s*["\']SE[ _]genes["\']\s*,?')
STATEMENT = re.compile(
    r"""^\s*(
          [A-Za-z_]\w*\s*\[\s*["']SE[ _]genes["']\s*\]\s*=      # rec["SE genes"] = ...
        | [A-Za-z_][\w\.]*\s*=.*["']SE[ _]genes["']             # x = ... "SE genes" ...
        | [A-Za-z_][\w\.]*\s*=.*\b(se_genes|n_genes_se|genes_se)\b
        | .*\.append\s*\(.*["']SE[ _]genes["']
    )""", re.X | re.I)


def hr(c="-"):
    print(c * 74)


def refuse(msg):
    hr("=")
    print("REFUSING — nothing was written.\n")
    print(msg)
    hr("=")
    sys.exit(2)


def is_comment(line):
    return line.lstrip().startswith("#")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-d", "--dir", default=".")
    ap.add_argument("--list", dest="list_only", action="store_true")
    ap.add_argument("--apply", action="store_true")
    a = ap.parse_args()

    path = os.path.join(a.dir, SCRIPT)

    hr("=")
    print('97 — remove Table 3\'s "SE genes" column (for real this time)')
    hr("=")

    # ---------------------------------------------------------------- G1
    if not os.path.exists(path):
        refuse(f"G1 — {path} not found. Pass -d with the directory holding it.")
    src = open(path, encoding="utf-8").read()
    lines = src.splitlines(keepends=True)
    print(f"\n  {SCRIPT}: {len(lines):,} lines")
    print("  G1 PASS")

    # ---------------------------------------------------------------- G2/G3
    hits = []
    for i, line in enumerate(lines, 1):
        if any(p.search(line) for p in LABEL_PATTERNS):
            hits.append((i, line))

    if "[97]" in src and not any(
            (not is_comment(l)) and any(p.search(l) for p in LABEL_PATTERNS)
            for l in lines):
        hr()
        print("  ALREADY APPLIED — this script's marker is present and no live mention")
        print("  of the column remains. Nothing to do.")
        print("  Re-run 87 if Table3.tsv predates the edit.")
        hr()
        return

    if not hits:
        refuse("G2 — no mention of the column found at all.\n\n"
               "  Searched for: \"SE genes\", se_genes, n_genes_se, genes_se.\n"
               "  87 printed an 'SE genes' header on the 12 Oct run, so the column is\n"
               "  built under a name this script does not know. Run:\n\n"
               "    grep -n -i 'gene' 87_build_table3.py\n\n"
               "  and send the output; the pattern list above needs the real name.")

    print(f"\n  G3  lines mentioning the column: {len(hits)}")
    to_edit, left = [], []
    for i, line in enumerate(lines, 1):
        if not any(p.search(line) for p in LABEL_PATTERNS):
            continue
        if is_comment(line):
            left.append((i, line, "already a comment"))
        elif DICT_ENTRY.search(line):
            to_edit.append((i, line, "drop dict entry"))
        elif STATEMENT.match(line):
            to_edit.append((i, line, "comment"))
        elif ELEMENT.search(line):
            to_edit.append((i, line, "drop element"))
        else:
            left.append((i, line, "mentions but does not bind — not edited"))

    for i, line, how in to_edit:
        print(f"    {i:>5}  {how:<12} {line.rstrip()[:80]}")
    for i, line, why in left:
        print(f"    {i:>5}  {'leave':<12} {line.rstrip()[:62]}   ({why})")

    if not to_edit and any("[97]" in l for l in lines):
        hr()
        print("  ALREADY APPLIED — this script's marker is present and no binding line")
        print("  remains. Nothing to do. Re-run 87 if Table3.tsv predates the edit.")
        hr()
        return

    if not to_edit:
        print("\n  Every mention is already a comment or does not bind the column.")
        print("  Nothing to do — but 87 still printed the column on 12 Oct, so the")
        print("  column is built elsewhere under another name. Run:")
        print("\n    grep -n -i 'gene' 87_build_table3.py\n")
        refuse("G2 — the column is present in 87's output but no binding line was found.\n"
               "  Refusing rather than reporting success on a file that still carries it.")
    print("  G3 PASS")

    if a.list_only or not a.apply:
        hr()
        print("  --list / dry run. Nothing written.")
        print(f"  Would comment out {len(to_edit)} line(s), then re-run 87.")
        hr()
        return

    # ---------------------------------------------------------------- edit
    out = list(lines)
    marker = "  # [97] dropped: no claim depends on a gene count; see RESULT_table3_built.md"
    for i, line, how in to_edit:
        nl = "\n" if line.endswith("\n") else ""
        if how == "comment":
            indent = line[:len(line) - len(line.lstrip())]
            out[i - 1] = f"{indent}# {line.rstrip()}{marker}{nl}".replace(
                f"# {indent}", "# ", 1)
        else:  # drop element / dict entry — keep the literal valid
            pat = DICT_ENTRY if how == "drop dict entry" else ELEMENT
            stripped = pat.sub("", line.rstrip("\n"), count=1)
            # a removal at the end of a line can leave ",]" or ", ]"
            stripped = re.sub(r",(\s*[\]\)])", r"\1", stripped)
            if not stripped.strip():
                stripped = line[:len(line) - len(line.lstrip())].rstrip("\n")
            out[i - 1] = f"{stripped}{marker}{nl}"

    new = "".join(out)

    # ---------------------------------------------------------------- G4
    try:
        ast.parse(new)
    except SyntaxError as e:
        refuse(f"G4 — the edited file would not parse: {e}\n\n"
               f"  Commenting out line {e.lineno} breaks the syntax, which means the "
               f"column is\n  bound inside a multi-line expression that cannot be "
               f"commented line by line.\n  Send these lines and they can be edited by "
               f"hand:\n\n" +
               "".join(f"    {i}: {l.rstrip()}\n" for i, l, _ in to_edit))
    print("\n  G4 PASS — the edited file parses.")

    b = f"{path}.bak.{time.strftime('%Y%m%d-%H%M%S')}"
    shutil.copy2(path, b)
    print(f"  backed up {path} -> {b}")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(new)
    print(f"  wrote {path}")

    # ---------------------------------------------------------------- G5
    chk = open(path, encoding="utf-8").read().splitlines(keepends=True)
    bad = []
    if len(chk) != len(lines):
        bad.append(f"line count changed: {len(lines)} -> {len(chk)}")
    edited = {i: how for i, _, how in to_edit}
    for i, (old, cur) in enumerate(zip(lines, chk), 1):
        how = edited.get(i)
        if how == "comment":
            if not is_comment(cur):
                bad.append(f"line {i} was not commented")
        elif how in ("drop element", "drop dict entry"):
            if any(pp.search(cur) for pp in LABEL_PATTERNS) and "[97]" not in cur:
                bad.append(f"line {i} still carries the column label")
        elif old != cur:
            bad.append(f"line {i} changed but should not have")
    print(f"\n  G5  lines edited: {len(edited)}   unintended changes: "
          f"{len([x for x in bad if 'should not' in x])}")
    if bad:
        shutil.copy2(b, path)
        refuse("G5 — the written file is not what was intended; the backup has been "
               "restored:\n" + "".join(f"    {x}\n" for x in bad))
    print("  G5 PASS")

    hr("=")
    print("DONE. Now rebuild the table:\n\n    python 87_build_table3.py\n")
    print('The "SE genes" column should be absent from the printed table and from')
    print("Table3.tsv. ⚠ If it is still there, the column is built in a second place —")
    print("send `grep -n -i gene 87_build_table3.py` rather than assuming it worked.")
    hr("=")


if __name__ == "__main__":
    main()
