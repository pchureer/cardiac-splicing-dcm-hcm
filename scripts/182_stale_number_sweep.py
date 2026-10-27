#!/usr/bin/env python3
"""
182_stale_number_sweep.py -- find 50 bp numbers that survived the 125 bp re-basing.

Seven sentences in the manuscript quoted 50 bp values after the re-quantification
had been adopted (QKI +0.033, "every interval covers the floor", +0.28/+0.19,
0.061 "not re-measured since", 38,084, "QKI included", the Table 1 booking).
Each was found by reading. This script is the mechanical version: every 50 bp
value the project record pairs with a 125 bp replacement, searched for in a
.docx or .md, printed with context. It cannot judge -- a labelled comparator
("from 926 to 813") is a legitimate hit -- but it guarantees nothing is missed
for want of looking.

  python 182_stale_number_sweep.py Manuscript_125bp_SR.docx
  python 182_stale_number_sweep.py Supplementary_Captions_rl125.md --skip-comments

Exit 1 if any hit; read each and decide.
"""
import argparse, re, sys

STALE = [
 (r"\b926\b","SE discovery count","813"), (r"\b1,?054\b","restricted count","939"), (r"\b899\b","nine-marker count","777"),
 (r"\b835\b","ladder row 3","732"), (r"\b866\b","ladder row 4","744"), (r"\b888\b","ladder row 5","770"),
 (r"62,?971","SE universe","53,456"), (r"156,?115","all-class universe","144,942"), (r"\b1,?445\b","all-class sig","1,294"),
 (r"\b468\b","HCM vs NF SE","431"), (r"48,?024","HCM universe","39,607"), (r"45,?380","DCMvsHCM universe","38,273"),
 (r"\b5,?301\b","FDR-only DCM","5,401"), (r"\b1,?355\b","FDR-only HCM","1,430"),
 (r"\b2,?163\b","measurable events","1,800"), (r"\b1,?636\b","concordant events","1,376"), (r"75\.6\s?%","pooled concordance","76.4%"), (r"10⁻¹³²|10-132","pooled P","2.5 × 10⁻¹¹⁷"),
 (r"\b0\.584\b","platform slope","0.690"), (r"\b0\.988\b","people slope (also a genotype-concordance bound: judge)","0.968"), (r"\b203\b","common events","161"), (r"\b0\.577\b","product","0.668"), (r"\b0\.521\b","Colorado slope / observed C","0.5999 / 0.600"),
 (r"\b0\.476\b","Brigham slope","0.510"), (r"\b0\.386\b","Stanford slope","0.401"), (r"\b0\.243\b","Fuwai slope","0.295"), (r"\b0\.406\b","mean slope","0.451"), (r"\b0\.445\b","LBNL slope","0.573"), (r"75\.9\s?%","LBNL concordance","78.2%"),
 (r"79\.6\s?%|81\.2\s?%|69\.2\s?%|76\.1\s?%","per-cohort concordance","80.8/81.9/70.8/75.8"), (r"375/471|359/442|494/714|408/536","concordance fractions","311/385, 298/364, 422/596, 345/455"),
 (r"\+?0\.033\b","QKI margin","+0.062 PTBP1"), (r"\bQKI\b","QKI led at 50 bp","PTBP1"), (r"[−-]0\.005\b","RBM20 SE margin","−0.051"), (r"[−-]0\.119","RBM20 SE CI","−0.159"), (r"\+0\.080\b","RBM20 SE CI","+0.047"), (r"[−-]0\.085|\+0\.104\b","RBM20 MXE CI","−0.108/+0.072"),
 (r"\b663\b","RBM20 SE targets","641"), (r"\b1,?393\b","RBM20 MXE targets","1,399"), (r"\b134\b","signature genes","131"), (r"\b130\b(?! and)","genes with events","123"), (r"\b2,?504\b","scan targets","1,738"), (r"\b5,?848\b","scan background","5,640/12,518"), (r"\b2,?755\b","MXE targets","2,717"), (r"\b11,?023\b","MXE background","11,103/18,718"), (r"53,?288","scan universe","40,272"),
 (r"\b927\b","enrichment foreground","805"), (r"12,?725","enrichment universe","12,564"),
 (r"\b590\b","retained nine-marker","498"), (r"63\.7\s?%","retention nine","61.3%"), (r"\b716\b","retained restricted","627"), (r"77\.3\s?%","retention restricted (also MXE 125 bp: judge)","77.1%"), (r"1\.588","duplication","1.515"), (r"100,?014","rows","80,973"), (r"\b991\b","filter-first","864"), (r"\b1,?074\b","full-coordinate","944"),
 (r"38,?084","LVEF universe","33,130"), (r"53\.6\s?%","LVEF skew","53.2%"), (r"0\.0559|0\.0538","LVEF rho/null","0.0555/0.0540"), (r"\b155\b","LVEF n","160"),
 (r"\b0\.061\b","Colorado MXE reverse","0.060"), (r"\b1,?477\b","its n","1,478"), (r"\b0\.447\b","Colorado SE reverse","0.443"), (r"\b453\b","its n","460"), (r"\b0\.725\b","MXE internal","0.730"), (r"\b0\.157\b","Colorado MXE r","0.158"),
 (r"\b0\.639\b","unselected mean","0.722"), (r"\b0\.533\b","forward mean","0.604"), (r"24\.9","model effect","28.2"),
 (r"thirteen (label )?permutations","wrong count","five shuffles, twenty analyses"), (r"\b715\b|50,?124|97\.11|\b864\b","JC figures","97.8% of 762"),
 (r"452-fold|1,?767-fold","retired ratios","orders of magnitude"), (r"at least 99\s?%","five-class direction (SE-only 100%: judge)","at least 98%"), (r"three of five","flipped","four of five"),
 (r"\b0\.997\b","pilot rho (ok only if labelled pilot)","0.996 on 332"), (r"by more than 0\.03","bare delta","—"),
 (r"\+0\.28\b|\+0\.19\b","sensitivity margins","+0.227/+0.314"), (r"1\.41×|1\.41x|1\.28×","asymmetry calibration","dropped"),
 (r"216,?2[13]\d","pre-rescue pair count","222,870 / 222,843"), (r"\b0\.937\b","pre-rescue concordance bound","0.928"),
 (r"not re-measured since","booking","measured"),
 (r"103,?054|103,?073","50 bp MAGNet permutation tested count","83,830/83,870 full; 54,879/54,874 exon"),
]

def read_docx(path):
    """Paragraph texts straight from word/document.xml -- no python-docx needed, so
    it runs in any environment (the base env on the cluster lacks docx)."""
    import zipfile, html
    xml = zipfile.ZipFile(path).read("word/document.xml").decode("utf-8")
    xml = re.sub(r"<w:del\b.*?</w:del>", "", xml, flags=re.S)        # tracked deletions
    out = []
    for i, para in enumerate(re.findall(r"<w:p[ >].*?</w:p>", xml, flags=re.S)):
        text = html.unescape("".join(re.findall(r"<w:t[^>]*>(.*?)</w:t>", para, flags=re.S)))
        if text.strip():
            out.append((i, text))
    return out


def read(path, skip_comments):
    if path.endswith(".docx"):
        return read_docx(path)
    t = open(path, encoding="utf-8").read()
    if skip_comments:
        t = re.sub(r"<!--.*?-->", "", t, flags=re.S)
    return [(i + 1, l) for i, l in enumerate(t.splitlines()) if l.strip()]

def main():
    ap = argparse.ArgumentParser(); ap.add_argument("file"); ap.add_argument("--skip-comments", action="store_true")
    ap.add_argument("--stop-at", default="References", help="paragraph text at which to stop (docx)")
    a = ap.parse_args()
    units = read(a.file, a.skip_comments)
    hits = 0
    for i, t in units:
        if a.file.endswith(".docx") and t.strip() == a.stop_at: break
        for pat, what, new in STALE:
            for m in re.finditer(pat, t):
                hits += 1; s = max(0, m.start() - 70); e = min(len(t), m.end() + 50)
                print(f"{i:5d} [{what} -> {new}]\n      …{t[s:e]}…")
    print(f"\n{hits} hit(s) in {a.file}; read each: a labelled 50 bp comparator is allowed, a bare one is not.")
    sys.exit(1 if hits else 0)

if __name__ == "__main__":
    main()
