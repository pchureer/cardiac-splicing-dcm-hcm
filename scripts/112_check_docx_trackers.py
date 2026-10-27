#!/usr/bin/env python3
"""
112_check_docx_trackers.py -- close the loop between the MANUSCRIPT and
`102_verify_pmids.py`.

THE GAP THIS FILLS
------------------
`102` verifies its own hard-coded REFS table against PubMed and Crossref. That
table is a MIRROR of the manuscript's reference list, transcribed by hand. If
the mirror is edited and the document is not, `102` prints "every tracker
matches the record it is attached to" while the .docx still carries the wrong
identifier. Both files then check out against themselves and disagree with each
other -- this project's signature failure mode, applied to references.

This script reads the trackers out of the .docx and diffs them against `102`'s
table. It imports the table from `102` directly, so there is one source of
truth and no second transcription.

CHECKS
------
  D1  every reference number 1..N present in the document's reference list
  D2  the tracker in the document equals the tracker in 102's table
  D3  a reference 102 declares NO_IDENTIFIER carries no tracker in the document
  D4  no FORBIDDEN identifier appears anywhere in the document
  D5  reference-list length matches 102's table

USAGE
-----
  python 112_check_docx_trackers.py Manuscript_with_TP.docx
  python 112_check_docx_trackers.py Manuscript_with_TP.docx --102 102_verify_pmids.py
  python 112_check_docx_trackers.py --self-test

Exit 0 only when the document and the table agree.
"""
import argparse
import importlib.util
import os
import re
import sys
import zipfile

PARA = re.compile(r"<w:p[ >].*?</w:p>", re.S)
TAG = re.compile(r"<[^>]+>")
# a reference-list entry starts with "12." or "12)" possibly after whitespace
REFNUM = re.compile(r"^\s*(\d{1,3})[.)]\s")
PMID = re.compile(r"PMID\s*:?\s*(\d{6,9})", re.I)
DOI = re.compile(r"(?:doi\s*:?\s*|doi\.org/)(10\.\d{4,9}/[^\s,;)\]]+)", re.I)


def load_102(path):
    spec = importlib.util.spec_from_file_location("v102", path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def paragraphs(docx_path):
    with zipfile.ZipFile(docx_path) as z:
        names = [n for n in z.namelist()
                 if n == "word/document.xml" or
                 re.match(r"word/(header|footer)\d*\.xml$", n)]
        if "word/document.xml" not in names:
            raise SystemExit("[fatal] no word/document.xml -- not a .docx?")
        xml = z.read("word/document.xml").decode("utf8", "replace")
    out = []
    for p in PARA.findall(xml):
        # drop deleted (tracked-change) runs so a deletion is not read as text
        p = re.sub(r"<w:del[ >].*?</w:del>", "", p, flags=re.S)
        txt = TAG.sub("", p)
        txt = (txt.replace("&amp;", "&").replace("&lt;", "<")
                  .replace("&gt;", ">").replace("&quot;", '"')
                  .replace("&#39;", "'").replace("&apos;", "'"))
        out.append(txt.strip())
    return out


def reference_entries(paras):
    """Numbered paragraphs, keeping the LAST run of consecutive numbering.

    A numbered paragraph can also be a numbered heading ('1. Introduction'),
    so we take the longest ascending run that reaches the highest number --
    the reference list.
    """
    cand = [(i, int(m.group(1)), t)
            for i, t in enumerate(paras)
            for m in [REFNUM.match(t)] if m]
    if not cand:
        return []
    runs, cur = [], []
    for item in cand:
        if cur and item[1] == cur[-1][1] + 1:
            cur.append(item)
        else:
            if cur:
                runs.append(cur)
            cur = [item]
    if cur:
        runs.append(cur)
    runs.sort(key=lambda r: (len(r), r[-1][1]))
    return runs[-1]


def tracker_of(text):
    m = PMID.search(text)
    if m:
        return m.group(1)
    m = DOI.search(text)
    if m:
        return m.group(1).rstrip(".")
    return None


def check(docx_path, v102, report):
    def w(s=""):
        report.append(s)

    refs = {r[0]: r for r in v102.REFS}
    no_id = getattr(v102, "NO_IDENTIFIER", {})
    forbidden = getattr(v102, "FORBIDDEN", [])

    paras = paragraphs(docx_path)
    entries = reference_entries(paras)
    w(f"document: {docx_path}")
    w(f"reference list: {len(entries)} numbered entries"
      + (f", {entries[0][1]}..{entries[-1][1]}" if entries else ""))
    w(f"102 table:      {len(refs)} references")
    w()

    fails = []

    # --- D5
    if len(entries) != len(refs):
        fails.append(f"D5 length: document has {len(entries)} entries, "
                     f"102 has {len(refs)}")

    doc_nums = {n for _, n, _ in entries}
    # --- D1
    missing = sorted(set(refs) - doc_nums)
    if missing:
        fails.append(f"D1 missing from the document's reference list: {missing}")
    extra = sorted(doc_nums - set(refs))
    if extra:
        fails.append(f"D1 in the document but not in 102's table: {extra}")

    # --- D2, D3
    for _, n, txt in entries:
        if n not in refs:
            continue
        doc_tr = tracker_of(txt)
        tbl_tr = refs[n][1]
        if n in no_id:
            if doc_tr:
                fails.append(
                    f"D3 ref {n}: 102 declares NO identifier, but the document "
                    f"still carries '{doc_tr}'. The manuscript was not edited.")
            continue
        if doc_tr is None:
            fails.append(f"D2 ref {n}: no PMID or DOI found in the document "
                         f"entry; 102 expects '{tbl_tr}'")
        elif doc_tr.lower() != tbl_tr.lower():
            fails.append(f"D2 ref {n}: document '{doc_tr}' vs 102 '{tbl_tr}'")

    # --- D4
    whole = "\n".join(paras)
    for f in forbidden:
        if f in whole:
            fails.append(f"D4 FORBIDDEN identifier '{f}' still appears in the "
                         f"document. It resolves to a different work.")

    if fails:
        w("FAILURES")
        for f in fails:
            w(f"  ✗ {f}")
        w()
        w("The manuscript and 102's table disagree. 102 passing means only that")
        w("its TABLE matches PubMed -- not that the document does. Fix the")
        w("document, or correct the table if the document is right.")
    else:
        w("✅ D1-D5 pass: every tracker in the document matches 102's table,")
        w("   declared no-identifier references carry none, and no retired")
        w("   identifier survives anywhere in the file.")
    return fails


def self_test():
    import tempfile

    def make_docx(path, paras):
        body = "".join(
            "<w:p><w:r><w:t>" + p.replace("&", "&amp;") + "</w:t></w:r></w:p>"
            for p in paras)
        xml = ('<?xml version="1.0"?><w:document '
               'xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main">'
               "<w:body>" + body + "</w:body></w:document>")
        with zipfile.ZipFile(path, "w") as z:
            z.writestr("word/document.xml", xml)

    class V:
        REFS = [(1, "20124440", "Kong SW", "J", "2010", "3", "138"),
                (2, "10.1007/x", "Akhigbe", "J", "2026", "28", "62"),
                (3, "NO-ID-BOOK", "Deming WE", "BOOK", "1943", "-", "-")]
        NO_IDENTIFIER = {3: ("Deming WE. ...", "a book")}
        FORBIDDEN = ["11682119", "10.2307/2965635"]

    d = tempfile.mkdtemp()
    good = os.path.join(d, "good.docx")
    make_docx(good, [
        "1. Introduction",                      # a numbered HEADING, not a ref
        "2. Methods",
        "References",
        "1. Kong SW, et al. J 2010;3:138-46. PMID: 20124440",
        "2. Akhigbe A, et al. J 2026;28:62. doi: 10.1007/x",
        "3. Deming WE. Statistical Adjustment of Data. New York: Wiley; 1943.",
    ])
    rep = []
    f = check(good, V, rep)
    assert not f, "\n".join(rep)

    bad = os.path.join(d, "bad.docx")
    make_docx(bad, [
        "References",
        "1. Kong SW, et al. PMID: 20124440",
        "2. Akhigbe A, et al. doi: 10.1007/WRONG",
        "3. Deming WE. Statistical Adjustment of Data. DOI: 10.2307/2965635",
    ])
    rep2 = []
    f2 = check(bad, V, rep2)
    joined = "\n".join(f2)
    assert "D2 ref 2" in joined, joined
    assert "D3 ref 3" in joined, joined
    assert "D4 FORBIDDEN" in joined, joined
    assert "D5" not in joined, "length should match"

    short = os.path.join(d, "short.docx")
    make_docx(short, ["References", "1. Kong SW. PMID: 20124440"])
    rep3 = []
    f3 = check(short, V, rep3)
    j3 = "\n".join(f3)
    assert "D5 length" in j3 and "D1 missing" in j3, j3

    print("[ok] self-test passed: numbered headings not mistaken for the "
          "reference list; wrong DOI, un-removed no-identifier tracker, "
          "forbidden identifier and a truncated list all caught")
    return 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("docx", nargs="?")
    ap.add_argument("--102", dest="v102", default="102_verify_pmids.py")
    ap.add_argument("--self-test", action="store_true")
    a = ap.parse_args()
    if a.self_test:
        return self_test()
    if not a.docx:
        print("[fatal] give the manuscript .docx", file=sys.stderr)
        return 2
    self_test()
    print()
    if not os.path.exists(a.v102):
        print(f"[fatal] {a.v102} not found", file=sys.stderr)
        return 2
    v = load_102(a.v102)
    rep = ["=" * 74, "112 — document vs 102's table", "=" * 74, ""]
    fails = check(a.docx, v, rep)
    print("\n".join(rep))
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
