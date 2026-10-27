#!/usr/bin/env python3
"""
181_assemble_supplement.py -- assemble the single Supplementary Information file that
Scientific Reports requires (one PDF, <= 50 MB) from the deposited pieces.

Order: Supplementary Methods (M1-M10) -> Supplementary Discussion -> Supplementary
Figures 1-2 -> Supplementary Tables 1-10 -> (optionally) the main-text Tables 1 and 3
moved to the supplement as Supplementary Tables 11-12, for the 8-display-item limit.

Nothing is computed here. Every caption comes from Supplementary_Captions_rl125.md,
every table from its deposited TSV, every figure from its deposited image. The script
REFUSES if a file is missing, if a caption still carries a [[...]] marker, if a table is
too long to typeset (that is Supplementary Data, not a table), or if a caption names a
file that is not the one being placed under it.

  python 181_assemble_supplement.py --list            # manifest, dates, row counts
  python 181_assemble_supplement.py                   # -> Supplementary_Information_rl125.md + .docx
  python 181_assemble_supplement.py --with-main-tables  # + Tables 1 and 3 as Supp Tables 11-12
  python 181_assemble_supplement.py --pdf             # also try pandoc -> PDF

The .docx is produced with pandoc. For the PDF, open the .docx in Word and export, or
  soffice --headless --convert-to pdf Supplementary_Information_rl125.docx
Then check the size (SR: <= 50 MB) and that every figure rendered.
"""

import argparse
import csv
import os
import re
import shutil
import subprocess
import sys
import time

CAPTIONS = "Supplementary_Captions_rl125.md"
METHODS = "Supplementary_Methods_rl125.md"
DISCUSSION = "Supplementary_Discussion_rl125.md"
TITLE = ("Alternative splicing in failing human myocardium is shared across "
         "cardiomyopathy etiologies and not explained by RBM20 or other candidate regulators")
MAX_ROWS = 120          # longer than this is Supplementary Data, not a typeset table

# kind, number, [(sub-label or None, file)], note
ITEMS = [
    ("Figure", 1, [(None, "S9.png")], "genotype identity; read-length independent"),
    ("Figure", 2, [(None, "SuppFig2_matched_panel_rl125.png")], "drawn by 179"),
    ("Table", 1, [(None, "S1a_rl125.tsv")], ""),
    ("Table", 2, [(None, "S1b_rl125.tsv")], ""),
    ("Table", 3, [(None, "S1c_markers.tsv")], "expression-level; invariant"),
    ("Table", 4, [(None, "S2_rl125.tsv")], "built by patched 81"),
    ("Table", 5, [(None, "S5_rl125.tsv")], ""),
    ("Table", 6, [(None, "S6.tsv")], "KEPT AT 50 bp BY DECISION; caption says so"),
    ("Table", 7, [(None, "S7_rl125.tsv")], ""),
    ("Table", 8, [("(a) label permutations", "S8a_permutations_rl125.tsv"),
                  ("(b) reverse replication", "S8b_reverse_replication_rl125.tsv")],
     "(a) MAGNet rows from job 411394 (183/184); Colorado rows as 7 Oct"),
    ("Table", 9, [("(a) cohort pairs", "S9a_cohort_pairs.tsv"),
                  ("(b) excluded samples", "S9b_excluded_samples.tsv")],
     "genotype-level; invariant"),
    ("Table", 10, [(None, "SupplementaryTable10_lvef_summary_rl125.tsv")],
     "per-event results are Supplementary Data, not typeset"),
]
MAIN_TABLES = [  # (new number, markdown file from the main text, old label)
    (11, "Table1.md", "Table 1"),
    (12, os.path.join("t3_rl125", "Table3_rl125.md"), "Table 3"),
]
SUPP_DATA = [("Supplementary Data 1",
              "SupplementaryData_lvef_all_events_rl125.tsv",
              "per-event ejection-fraction results, 33,130 skipped exons; columns "
              "described in SupplementaryData_lvef_all_events_rl125.columns.md")]


def die(msg):
    print("\n[REFUSE] " + msg + "\n", file=sys.stderr)
    sys.exit(2)


def strip_comments(text):
    return re.sub(r"<!--.*?-->", "", text, flags=re.S)


def parse_captions(path):
    """{('Table', 4): caption text} from the ## headings and their blockquotes."""
    if not os.path.isfile(path):
        die(f"{path} not found")
    text = strip_comments(open(path, encoding="utf-8").read())
    caps = {}
    blocks = re.split(r"^## Supplementary (Table|Figure) (\d+)[^\n]*\n", text, flags=re.M)
    # blocks = [preamble, kind, num, body, kind, num, body, ...]
    for i in range(1, len(blocks), 3):
        kind, num, body = blocks[i], int(blocks[i + 1]), blocks[i + 2]
        lines = []
        for ln in body.splitlines():
            if ln.startswith(">"):
                lines.append(ln[1:].lstrip())
            elif lines and ln.strip() == "":
                lines.append("")
        cap = "\n".join(lines).strip()
        cap = re.sub(r"\n{3,}", "\n\n", cap)
        cap = cap.replace("⚠ ", "").replace("⚠", "")   # Word renders the glyph as an emoji image
        if cap:
            caps[(kind, num)] = cap
    return caps


def read_tsv(path):
    with open(path, encoding="utf-8", newline="") as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    rows = [r for r in rows if any(c.strip() for c in r)]
    return rows


PAGE_BREAK = ('```{=openxml}\n<w:p><w:r><w:br w:type="page"/></w:r></w:p>\n```\n\n'
              '```{=latex}\n\\newpage\n```\n')
MAX_COLS = 7            # wider tables are split into column blocks that share the id columns
ID_COLS = 2


def _one_table(rows):
    def esc(c):
        c = (c or "").replace("|", "\\|").replace("\n", " ").strip()
        return c if c.lower() not in ("nan", "none") else ""
    width = max(len(r) for r in rows)
    rows = [r + [""] * (width - len(r)) for r in rows]
    # header tokens like n_tested_exon wrap mid-token in Word; spaces wrap cleanly
    head = "| " + " | ".join(esc(c).replace("_", " ") for c in rows[0]) + " |"
    # pandoc sizes columns by the relative length of the dash runs when any line is
    # wider than the page, so make each run proportional to the column's widest cell
    # (header counted, floor of 3) -- otherwise "contrast" wraps one word per line.
    widths = [max(12, min(24, max(len(esc(r[j])) for r in rows))) for j in range(width)]   # floor 12 / cap 24: short id columns must not wrap letter by letter
    # pandoc applies the relative widths only when the separator line is wider than
    # its --columns (72); scale the runs up so it always is.
    scale = max(1.0, 96.0 / sum(widths))
    widths = [int(round(w * scale)) for w in widths]
    sep = "|" + "|".join("-" * w for w in widths) + "|"
    body = ["| " + " | ".join(esc(c) for c in r) + " |" for r in rows[1:]]
    return "\n".join([head, sep] + body)


def pipe_table(rows):
    """Pipe table; a table wider than MAX_COLS is split into blocks of columns, each
    block repeating the first ID_COLS identifier columns, because a portrait page
    cannot carry 15 columns legibly and pandoc cannot rotate one."""
    width = max(len(r) for r in rows)
    rows = [r + [""] * (width - len(r)) for r in rows]
    if width <= MAX_COLS:
        return _one_table(rows)
    out, start, k = [], ID_COLS, 1
    per = MAX_COLS - ID_COLS
    while start < width:
        cols = list(range(ID_COLS)) + list(range(start, min(width, start + per)))
        block = [[r[c] for c in cols] for r in rows]
        out.append(f"*columns {start + 1}–{min(width, start + per)} of {width}*\n\n"
                   + _one_table(block))
        start += per
        k += 1
    return "\n\n".join(out)


def when(path):
    return time.strftime("%Y-%m-%d", time.localtime(os.stat(path).st_mtime))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--with-main-tables", action="store_true",
                    help="append main-text Tables 1 and 3 as Supplementary Tables 11-12")
    ap.add_argument("--allow-markers", action="store_true",
                    help="do not refuse on [[...]] markers in captions (drafts only)")
    ap.add_argument("--pdf", action="store_true", help="also try pandoc -> PDF")
    ap.add_argument("--out", default="Supplementary_Information_rl125")
    a = ap.parse_args()

    # ---- manifest and gates
    problems = []
    print(f"{'item':<26}{'file':<48}{'date':<12}rows   note")
    for kind, num, files, note in ITEMS:
        for sub, f in files:
            if not os.path.isfile(f):
                problems.append(f"Supplementary {kind} {num}: {f} missing")
                print(f"{'Supp ' + kind + ' ' + str(num):<26}{f:<48}{'MISSING':<12}")
                continue
            n = len(read_tsv(f)) - 1 if f.endswith(".tsv") else ""
            if f.endswith(".tsv") and n > MAX_ROWS:
                problems.append(f"Supplementary {kind} {num}: {f} has {n} rows (> {MAX_ROWS}); "
                                f"that is Supplementary Data, not a typeset table")
            print(f"{'Supp ' + kind + ' ' + str(num):<26}{f:<48}{when(f):<12}{str(n):<7}{note}")
    for src in (METHODS, DISCUSSION, CAPTIONS):
        if not os.path.isfile(src):
            problems.append(f"{src} missing")
        else:
            print(f"{'text':<26}{src:<48}{when(src)}")
    if a.with_main_tables:
        for num, f, old in MAIN_TABLES:
            if not os.path.isfile(f):
                problems.append(f"Supplementary Table {num}: {f} missing")
            else:
                print(f"{'Supp Table ' + str(num):<26}{f:<48}{when(f):<12}{'':<7}was {old}")
    caps = parse_captions(CAPTIONS) if os.path.isfile(CAPTIONS) else {}
    wanted = [(k, n) for k, n, _, _ in ITEMS]
    if a.with_main_tables:
        wanted += [("Table", n) for n, _, _ in MAIN_TABLES]
    for kind, num in wanted:
        c = caps.get((kind, num))
        if not c:
            problems.append(f"no caption for Supplementary {kind} {num} in {CAPTIONS}")
        elif "[[" in c and not a.allow_markers:
            m = re.search(r"\[\[[^\]]{0,60}", c)
            problems.append(f"Supplementary {kind} {num} caption still carries a marker: "
                            f"{m.group(0)}…")
    if problems:
        print()
        for p in problems:
            print("  ✗ " + p)
        die(f"{len(problems)} problem(s); nothing assembled.")
    print("\n  every file present, every caption clean")
    if a.list:
        return

    # ---- assemble markdown
    md = []
    md.append(f"# Supplementary Information\n\n**{TITLE}**\n")
    md.append("## Contents\n")
    md.append("Supplementary Methods (M1–M10) · Supplementary Discussion · "
              "Supplementary Figures 1–2 · Supplementary Tables 1–"
              + ("12" if a.with_main_tables else "10") + " · Supplementary Data\n")
    md.append(PAGE_BREAK)

    meth = strip_comments(open(METHODS, encoding="utf-8").read()).strip()
    meth = re.sub(r"^# Supplementary Methods\s*$", "# Supplementary Methods", meth, flags=re.M)
    md.append(meth + "\n\n" + PAGE_BREAK)
    disc = strip_comments(open(DISCUSSION, encoding="utf-8").read()).strip()
    md.append(disc + "\n\n" + PAGE_BREAK)

    md.append("# Supplementary Figures\n")
    for kind, num, files, note in ITEMS:
        if kind != "Figure":
            continue
        f = files[0][1]
        md.append(f"![]({f}){{width=100%}}\n")
        md.append(caps[(kind, num)] + "\n\n" + PAGE_BREAK)

    md.append("# Supplementary Tables\n")
    for kind, num, files, note in ITEMS:
        if kind != "Table":
            continue
        md.append(caps[(kind, num)] + "\n")
        for sub, f in files:
            if sub:
                md.append(f"**{sub}**\n")
            md.append(pipe_table(read_tsv(f)) + "\n")
        md.append(PAGE_BREAK)

    if a.with_main_tables:
        for num, f, old in MAIN_TABLES:
            body = strip_comments(open(f, encoding="utf-8").read()).strip()
            # The legend comes from the captions file (the manuscript's Tables section,
            # re-based); the md file's own one-line label, if any, is dropped so the
            # table is not captioned twice.
            lines = body.splitlines()
            if lines and lines[0].startswith("**" + old):
                lines = lines[1:]
            body = "\n".join(lines).strip()
            md.append(caps[("Table", num)] + "\n")
            md.append(body + "\n\n" + PAGE_BREAK)

    md.append("# Supplementary Data\n")
    for label, f, desc in SUPP_DATA:
        state = "deposited" if os.path.isfile(f) else "⚠ NOT FOUND"
        md.append(f"**{label}.** `{f}` — {desc}. ({state})\n")

    text = "\n".join(md)
    with open(a.out + ".md", "w", encoding="utf-8") as fh:
        fh.write(text)
    print(f"\n[write] {a.out}.md  ({len(text.split()):,} words)")

    if not shutil.which("pandoc"):
        print("pandoc not found: convert the .md yourself, or `conda install -c conda-forge pandoc`")
        return
    cmd = ["pandoc", a.out + ".md", "-o", a.out + ".docx", "--from", "markdown+pipe_tables"]
    if os.path.isfile("reference_supp.docx"):
        cmd += ["--reference-doc", "reference_supp.docx"]   # 8-pt table text, 2 cm margins
        print("using reference_supp.docx (table text 8 pt, margins 2 cm)")
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode:
        die("pandoc -> docx failed:\n" + r.stderr)
    print(f"[write] {a.out}.docx  ({os.path.getsize(a.out + '.docx')/1e6:.1f} MB)")
    if a.pdf:
        r = subprocess.run(["pandoc", a.out + ".md", "-o", a.out + ".pdf", "-V", "geometry:margin=2cm"],
                           capture_output=True, text=True)
        if r.returncode:
            print("pandoc -> pdf failed (no LaTeX engine?). Export the .docx from Word, or:\n"
                  f"  soffice --headless --convert-to pdf {a.out}.docx")
        else:
            print(f"[write] {a.out}.pdf  ({os.path.getsize(a.out + '.pdf')/1e6:.1f} MB; SR limit 50 MB)")
    print("\nOPEN IT. Check: both figures render; no table lost columns; Supplementary Table 6 "
          "caption carries its 50 bp label; the title matches the manuscript.")


if __name__ == "__main__":
    main()
