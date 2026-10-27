#!/usr/bin/env python3
"""
102_verify_pmids.py — check every PMID and DOI tracker in the manuscript against
the record it is attached to.

WHY THIS IS A SCRIPT AND NOT AN ANSWER
  PubMed refused the session that drafted this (reCAPTCHA interstitial, then HTTP
  429). One record got through before the block and it was wrong, which is exactly
  why the other 37 cannot be taken on trust. This runs from the cluster, which can
  reach NCBI.

WHAT IT CHECKS, per reference
  surname of the first author   HARD — a mismatch means the wrong paper
  publication year             HARD — same
  volume                       HARD
  first page / article number   HARD, compared against BOTH `pages` and `elocationid`
  journal                      SOFT — PubMed uses NLM abbreviations, which differ
                               legitimately from house style

  For DOI-only trackers it asks PubMed whether a PMID exists at all, which is the
  question that matters: "no PMID" must be a measured fact, not an assumption.

⚠ READ-ONLY. Writes one file (`pmid_verification.tsv`) and touches nothing else.

USAGE
    python 102_verify_pmids.py --list          # the table, no network
    python 102_verify_pmids.py                 # verify
    python 102_verify_pmids.py --only 27 28 29
"""

import argparse
import json
import re
import sys
import time
import urllib.parse
import urllib.request

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
CROSSREF = "https://api.crossref.org/works/"
TOOL = "manuscript_ref_check"

# n, tracker, first-author surname+initials, journal, year, volume, first page
REFS = [
    (1,  "20124440", "Kong SW", "Circ Cardiovasc Genet", "2010", "3", "138"),
    (2,  "30051286", "Beqqali A", "Biophys Rev", "2018", "10", "1061"),
    (3,  "23775418", "Lara-Pezzi E", "J Cardiovasc Transl Res", "2013", "6", "945"),
    (4,  "19712804", "Brauch KM", "J Am Coll Cardiol", "2009", "54", "930"),
    (5,  "38288598", "Gregorich ZR", "Circ Genom Precis Med", "2024", "17", "e004355"),
    (6,  "42177204", "Radke MH", "Nat Commun", "2026", "17", "4607"),
    (7,  "41877661", "Furtado M", "Cardiovasc Res", "2026", "122", "906"),
    (8,  "41902389", "Wang Y", "Cardiovasc Res", "2026", "122", "808"),
    (9,  "10.1007/s11886-026-02386-0", "Akhigbe F", "Curr Cardiol Rep", "2026", "28", "62"),
    (10, "33517678", "Liu CF", "Circ Genom Precis Med", "2021", "14", "e003142"),
    (11, "35732739", "Chaffin M", "Nature", "2022", "608", "174"),
    (12, "39170450", "Li J", "Heliyon", "2024", "10", "e35408"),
    (13, "41017471", "Pan T", "Circulation", "2025", "152", "1501"),
    (14, "36345035", "Britto-Borges T", "J Transl Med", "2022", "20", "513"),
    (15, "25528681", "Liu Y", "Genomics", "2015", "105", "83"),
    (16, "30419824", "Sweet ME", "BMC Genomics", "2018", "19", "812"),
    (17, "36130500", "Spurrell CH", "Cell Rep", "2022", "40", "111400"),
    (18, "33558530", "Maron BA", "Nat Commun", "2021", "12", "873"),
    (19, "31197155", "Liu X", "Sci Data", "2019", "6", "90"),
    (20, "33590861", "Danecek P", "Gigascience", "2021", "10", "giab008"),
    (21, "20644199", "McKenna A", "Genome Res", "2010", "20", "1297"),
    (22, "20926424", "Manichaikul A", "Bioinformatics", "2010", "26", "2867"),
    (23, "24695404", "Bolger AM", "Bioinformatics", "2014", "30", "2114"),
    (24, "23104886", "Dobin A", "Bioinformatics", "2013", "29", "15"),
    (25, "39565199", "Mudge JM", "Nucleic Acids Res", "2025", "53", "D966"),
    (26, "38396040", "Wang Y", "Nat Protoc", "2024", "19", "1083"),
    (27, "10.1111/j.2517-6161.1995.tb02031.x", "Benjamini Y",
     "J R Stat Soc Series B", "1995", "57", "289"),
    (28, "10.2307/2276774", "Wilson EB", "J Am Stat Assoc", "1927", "22", "209"),
    (29, "NO-ID-BOOK", "Deming WE", "BOOK", "1943", "-", "-"),
    (30, "21605215", "Zaykin DV", "J Evol Biol", "2011", "24", "1836"),
    (31, "24960161", "Maatz H", "J Clin Invest", "2014", "124", "3419"),
    (32, "33805770", "Maimaiti R", "Int J Mol Sci", "2021", "22", "2928"),
    (33, "33397958", "Chen X", "Nat Commun", "2021", "12", "89"),
    (34, "32276354", "Fochi S", "Genes", "2020", "11", "402"),
    (35, "19075228", "Kalsotra A", "Proc Natl Acad Sci U S A", "2008", "105", "20333"),
    (36, "15988035", "Ladd AN", "Mol Cell Biol", "2005", "25", "6267"),
    (37, "15652482", "Xu X", "Cell", "2005", "120", "59"),
    (38, "34403477", "Zhao W", "Nucleic Acids Res", "2022", "50", "D287"),
    (39, "37144459", "Kolberg L", "Nucleic Acids Res", "2023", "51", "W207"),
    (40, "42761222", "Li J", "Front Bioinform", "2026", "6", "1934275"),
    (41, "40034311", "Conning-Rowland M", "Heliyon", "2025", "11", "e42499"),
]

# A correction already applied to the manuscript, kept here so the reason survives and
# the wrong value cannot drift back in. G0 reprints it.
CORRECTED = {
    27: ("PMID 11682119 -> DOI 10.1111/j.2517-6161.1995.tb02031.x",
         "11682119 is Behav Brain Res 2001;125:279-84, 'Controlling the false "
         "discovery rate in behavior genetics research' - a different paper, six "
         "years later, by an overlapping author.",
         "Benjamini & Hochberg 1995 (J R Stat Soc B 57:289-300) is NOT indexed in "
         "PubMed; the DOI is the only usable identifier. Applied 13 Oct 2026."),
}

# A page/article-number that PubMed does not carry, verified against another
# registry instead. The waiver applies ONLY while PubMed records no page at all:
# if PubMed later supplies a page that disagrees, this still fails. The evidence
# is printed on every run so the waiver cannot become blind.
PAGE_VERIFIED_ELSEWHERE = {
    6: ("Crossref 10.1038/s41467-026-73230-w",
        "article-number 4607, volume 17, 2026, Radke, Nature Communications; "
        "Crossref records no page field either. Nature Communications uses "
        "article numbers, and PubMed's elocationid carries the DOI rather than "
        "the article number for this record. Verified 14 Oct 2026."),
}

# References with NO usable identifier, with the evidence that none exists.
# G2 exempts these; G4 reprints the reason every run so the gap stays deliberate.
NO_IDENTIFIER = {
    29: ("Deming WE. Statistical Adjustment of Data. New York: John Wiley & "
         "Sons; 1943.",
         "A BOOK. Not in PubMed and it has no DOI of its own. The DOI "
         "10.2307/2965635 previously used here resolves to Smith, J Am Stat "
         "Assoc 1945;40:380 -- a ONE-PAGE item carrying the book's exact title, "
         "i.e. a REVIEW of the book, by a different author two years later. "
         "The title matches perfectly, so a title check passes it; only the "
         "author and year disagree. Removed 14 Oct 2026."),
}

# If any of these strings reappears in a tracker, the correction has been undone.
FORBIDDEN = ["11682119", "10.2307/2965635"]


def hr(c="-"):
    print(c * 78)


def get(url):
    req = urllib.request.Request(url, headers={"User-Agent": TOOL})
    with urllib.request.urlopen(req, timeout=60) as fh:
        return json.loads(fh.read().decode("utf-8", "replace"))


def norm(s):
    return re.sub(r"[^a-z0-9]", "", str(s).lower())


def surname(author_field):
    """'Kong SW' -> 'kong'. Handles compound surnames such as 'Lara-Pezzi E'."""
    parts = str(author_field).strip().split()
    if len(parts) >= 2 and re.fullmatch(r"[A-Za-z]{1,4}", parts[-1]):
        return norm(" ".join(parts[:-1]))
    return norm(author_field)


def esummary(pmids):
    out = {}
    for i in range(0, len(pmids), 20):
        chunk = pmids[i:i + 20]
        q = {"db": "pubmed", "retmode": "json", "id": ",".join(chunk), "tool": TOOL}
        d = get(EUTILS + "esummary.fcgi?" + urllib.parse.urlencode(q))
        res = d.get("result", {})
        for p in chunk:
            out[p] = res.get(p)
        time.sleep(0.4)
    return out


def esearch_doi(doi):
    q = {"db": "pubmed", "retmode": "json", "tool": TOOL,
         "term": '"' + doi + '"[AID] OR "' + doi + '"[DOI]'}
    d = get(EUTILS + "esearch.fcgi?" + urllib.parse.urlencode(q))
    time.sleep(0.4)
    return d.get("esearchresult", {}).get("idlist", [])


def crossref(doi):
    try:
        d = get(CROSSREF + urllib.parse.quote(doi))
    except Exception as exc:
        return {"error": str(exc)}
    m = d.get("message", {})
    auth = m.get("author") or [{}]
    yr = ""
    for k in ("published-print", "published-online", "issued"):
        if m.get(k, {}).get("date-parts"):
            yr = str(m[k]["date-parts"][0][0])
            break
    return {"title": (m.get("title") or [""])[0],
            "author": auth[0].get("family", ""),
            "journal": (m.get("container-title") or [""])[0],
            "year": yr, "volume": m.get("volume", ""),
            "page": m.get("page", "") or m.get("article-number", ""),
            "type": m.get("type", "")}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--list", action="store_true", help="print the table, no network")
    ap.add_argument("--only", nargs="*", type=int, default=None)
    ap.add_argument("--out", default="pmid_verification.tsv")
    a = ap.parse_args()

    refs = REFS
    if a.only:
        refs = [r for r in REFS if r[0] in a.only]

    hr("=")
    print("102 — verify every PMID / DOI tracker against the cited record")
    hr("=")

    # ------------------------------------------------------------------ G1
    ns = [r[0] for r in REFS]
    if ns != list(range(1, 42)):
        print("G1 REFUSING — reference numbers are not 1..41 contiguous:", ns)
        sys.exit(2)
    print("  G1 PASS  41 references, numbered 1..41 with no gaps")

    # ------------------------------------------------------------------ G2
    bad = [r[0] for r in REFS
           if r[0] not in NO_IDENTIFIER
           and not (r[1].isdigit() or r[1].startswith("10."))]
    if bad:
        print("G2 REFUSING — trackers that are neither a PMID nor a DOI:", bad)
        sys.exit(2)
    print("  G2 PASS  every tracker is a PMID or a DOI, or a declared no-identifier reference")

    # ------------------------------------------------------------------ G0
    print("\n  G0  corrections already applied, reprinted so they cannot drift back:")
    for n, (change, why, note) in CORRECTED.items():
        print("      ref " + str(n) + "  " + change)
        print("          " + why)
        print("          " + note)
    back = [r[0] for r in REFS for f in FORBIDDEN if f in r[1]]
    if back:
        print("\n  G0 REFUSING — a corrected tracker has reappeared in refs: "
              + str(back))
        print("     " + str(FORBIDDEN) + " must not appear in any tracker.")
        sys.exit(2)
    print("  G0 PASS  no corrected tracker has reappeared")

    # ------------------------------------------------------------------ G3
    if PAGE_VERIFIED_ELSEWHERE:
        print("\n  G3  page/article-numbers PubMed does not carry, verified "
              "against another registry:")
        for _n in sorted(PAGE_VERIFIED_ELSEWHERE):
            _src, _why = PAGE_VERIFIED_ELSEWHERE[_n]
            print("      ref " + str(_n) + "  via " + _src)
            print("          " + _why)
        print("      The waiver holds ONLY while PubMed records no page at all. "
              "A disagreeing PubMed page still fails.")

    # ------------------------------------------------------------------ G4
    if NO_IDENTIFIER:
        print("\n  G4  references with NO identifier, deliberately:")
        for _n in sorted(NO_IDENTIFIER):
            _cite, _why = NO_IDENTIFIER[_n]
            print("      ref " + str(_n) + "  " + _cite)
            print("          " + _why)

    if a.list:
        print("\n  --list: the table, nothing fetched\n")
        hr()
        for n, tr, au, jr, yr, vol, pg in refs:
            kind = "PMID" if tr.isdigit() else "DOI "
            print(kind + " " + tr.ljust(30) + " ref " + str(n).rjust(2) + "  "
                  + au.ljust(20) + " " + jr.ljust(26) + " " + yr + ";" + vol + ":" + pg)
            if len(tr) > 30:
                pass
        hr()
        print("\n  " + str(sum(1 for r in refs if r[1].isdigit())) + " PMIDs to check, "
              + str(sum(1 for r in refs
                          if not r[1].isdigit() and r[0] not in NO_IDENTIFIER))
              + " DOIs to resolve, "
              + str(len(NO_IDENTIFIER)) + " with no identifier by design")
        return

    pmids = [r[1] for r in refs if r[1].isdigit()]
    print("\n  fetching " + str(len(pmids)) + " PubMed records ...")
    try:
        recs = esummary(pmids)
    except Exception as exc:
        print("\nREFUSING — PubMed could not be reached: " + str(exc))
        print("  Nothing was written. This is a network failure, not a result.")
        sys.exit(2)

    rows = []
    hard_fail, soft = [], []
    hr("=")
    print("PMID TRACKERS")
    hr("=")
    for n, tr, au, jr, yr, vol, pg in refs:
        if not tr.isdigit():
            continue
        r = recs.get(tr)
        if not r or "error" in (r or {}):
            hard_fail.append((n, tr, "no PubMed record returned"))
            print("  ref " + str(n).rjust(2) + "  PMID " + tr + "  *** NO RECORD")
            continue
        got_au = (r.get("sortfirstauthor") or
                  (r.get("authors") or [{}])[0].get("name", ""))
        got_jr = r.get("source", "")
        got_yr = (r.get("pubdate", "") or "")[:4]
        got_vol = r.get("volume", "")
        got_pg = r.get("pages", "") or ""
        got_el = r.get("elocationid", "") or ""
        got_ti = r.get("title", "")

        probs = []
        if surname(au) != surname(got_au):
            probs.append("AUTHOR: want " + au + ", PubMed " + got_au)
        if yr != got_yr:
            probs.append("YEAR: want " + yr + ", PubMed " + got_yr)
        if norm(vol) != norm(got_vol) and vol != "-":
            probs.append("VOLUME: want " + vol + ", PubMed " + got_vol)
        pgs = norm(got_pg) + " " + norm(got_el)
        if norm(pg) not in pgs and pg != "-":
            if n in PAGE_VERIFIED_ELSEWHERE and not got_pg.strip():
                # PubMed carries no page for this record; verified elsewhere.
                src, _why = PAGE_VERIFIED_ELSEWHERE[n]
                soft.append((n, "page/article-number " + pg
                             + " not in PubMed (no page field); verified against "
                             + src))
            else:
                probs.append("PAGE: want " + pg + ", PubMed pages='" + got_pg
                             + "' eloc='" + got_el + "'")
        jr_ok = norm(jr) == norm(got_jr)

        flag = "OK  " if not probs else "FAIL"
        print("  ref " + str(n).rjust(2) + "  PMID " + tr + "  " + flag + "  "
              + got_au.ljust(20) + " " + got_jr.ljust(26) + " "
              + got_yr + ";" + got_vol + ":" + (got_pg or got_el))
        if probs:
            for p in probs:
                print("        *** " + p)
            print("        PubMed title: " + got_ti[:90])
            hard_fail.append((n, tr, "; ".join(probs)))
        if not jr_ok:
            soft.append((n, "journal: manuscript '" + jr + "' vs NLM '" + got_jr + "'"))
        rows.append([n, tr, flag, got_au, got_jr, got_yr, got_vol,
                     got_pg or got_el, got_ti])

    hr("=")
    print("DOI TRACKERS — does a PMID exist for these?")
    hr("=")
    for n, tr, au, jr, yr, vol, pg in refs:
        if tr.isdigit():
            continue
        if n in NO_IDENTIFIER:
            # declared no-identifier: nothing to resolve. Printed by G4.
            print("  ref " + str(n).rjust(2) + "  NO IDENTIFIER (book) — "
                  "nothing to resolve; see G4 above")
            continue
        ids = esearch_doi(tr)
        cr = crossref(tr)
        print("  ref " + str(n).rjust(2) + "  DOI " + tr)
        if "error" in cr:
            print("        Crossref: could not resolve (" + cr["error"] + ")")
            soft.append((n, "Crossref did not resolve " + tr))
        else:
            print("        Crossref: " + cr["author"] + " | " + cr["journal"]
                  + " | " + cr["year"] + ";" + str(cr["volume"]) + ":"
                  + str(cr["page"]) + " | " + cr["type"])
            print("        title: " + cr["title"][:90])
            if surname(au) != norm(cr["author"]):
                print("        *** AUTHOR: want " + au + ", Crossref "
                      + cr["author"])
                hard_fail.append((n, tr, "Crossref author mismatch"))
            if yr != cr["year"] and yr != "-":
                print("        *** YEAR: want " + yr + ", Crossref " + cr["year"])
                hard_fail.append((n, tr, "Crossref year mismatch"))
        if ids:
            print("        >>> a PMID DOES exist: " + ", ".join(ids)
                  + "  (add it to your working copy)")
        else:
            print("        no PMID in PubMed — correct to leave this one DOI-only")
        rows.append([n, tr, "DOI", cr.get("author", ""), cr.get("journal", ""),
                     cr.get("year", ""), cr.get("volume", ""), cr.get("page", ""),
                     cr.get("title", "")])

    with open(a.out, "w") as fh:
        fh.write("ref\ttracker\tstatus\tauthor\tjournal\tyear\tvolume\tpages\ttitle\n")
        for r in rows:
            fh.write("\t".join(str(x).replace("\t", " ") for x in r) + "\n")

    hr("=")
    print("SUMMARY")
    hr("=")
    print("  references checked : " + str(len(refs)))
    print("  hard failures      : " + str(len(hard_fail)))
    for n, tr, why in hard_fail:
        print("      ref " + str(n) + " (" + tr + "): " + why)
    print("  journal-name notes : " + str(len(soft))
          + "  (NLM abbreviations differ from house style; not errors)")
    for n, why in soft:
        print("      ref " + str(n) + ": " + why)
    print("\n  [write] " + a.out)
    hr("=")
    if hard_fail:
        print("\n⚠ Fix the hard failures before the EndNote import. A wrong PMID")
        print("  imports a different paper and the error survives every later check,")
        print("  because the reference list will still read correctly.")
        sys.exit(1)
    print("\n✅ every tracker matches the record it is attached to")


if __name__ == "__main__":
    main()
