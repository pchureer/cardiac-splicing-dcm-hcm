#!/usr/bin/env python3
"""
99_pubmed_initials_pass.py — settle every author initial in one command.

WHY THIS EXISTS
---------------
The reference list was built with full author lists read from **Crossref publisher
deposits** and **publisher pages**, because PubMed, PMC and Europe PMC were all
captcha- or policy-blocked from the machine that drafted it. NCBI E-utilities is
blocked too (`CONNECT tunnel failed, response 403`).

Crossref and OpenAlex disagree on *middle* initials for two references:

    GATK / McKenna 2010   Hanna M vs MG · Kernytsky A vs AM · Garimella K vs KV
                          Gabriel S vs SL · Daly M vs MJ
    KING / Manichaikul    Manichaikul A vs AW · Sale M vs MM

and GENCODE 2025 has two compound Spanish/Catalan surnames whose Vancouver form is a
house-style question (Gonzalez Martinez J, Martinez Gomez L).

This script runs one E-utilities pass over every reference and prints **PubMed's own
rendering** of each author list, diffed against the drafted list. Run it from any machine
with ordinary internet — the cluster frontend, or a laptop. It needs only Python 3 and
the standard library; no pip install, no API key.

    python 99_pubmed_initials_pass.py                  # all references
    python 99_pubmed_initials_pass.py --only M6 M7 M3  # just the disputed ones
    python 99_pubmed_initials_pass.py --email you@example.com   # polite to NCBI

⚠ WHAT PUBMED IS AND IS NOT THE AUTHORITY FOR
---------------------------------------------
**PubMed is the conventional Vancouver source, not the ultimate one.** NLM derives
initials from the article's own title page, so where they disagree, **the paper's title
page wins** — that is what the journal printed and what a copy-editor will check against.
This script's output is therefore a *strong default*, and for any reference it flags,
the final answer is the PDF's byline.

That distinction is not academic here. The project has already seen an institutional
repository record list an author who is not on the paper at all (the MDC Berlin record
for Guo 2012 gives "Lincoln, S.E." where the paper has Shijun Li). **Aggregators are
convenient and occasionally wrong; title pages are not.**

HOW TO READ THE OUTPUT
  AGREES    PubMed renders the list exactly as drafted. Nothing to do.
  DIFFERS   the differing positions are printed side by side. Decide from the title page.
  NO PMID   not indexed, or the PMID is unknown. Falls back to a DOI lookup and says so.
"""

import argparse
import json
import sys
import time
import urllib.parse
import urllib.request

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

# slot: (pmid or None, doi, drafted author string)
REFS = {
    "M1":  ("24695404", "10.1093/bioinformatics/btu170",
            "Bolger AM, Lohse M, Usadel B"),
    "M2":  ("23104886", "10.1093/bioinformatics/bts635",
            "Dobin A, Davis CA, Schlesinger F, Drenkow J, Zaleski C, Jha S, Batut P, "
            "Chaisson M, Gingeras TR"),
    "M3":  ("39565199", "10.1093/nar/gkae1078", "[56 authors — see REFERENCES_draft.md]"),
    "M4":  ("38396040", "10.1038/s41596-023-00944-2",
            "Wang Y, Xie Z, Kutschera E, Adams JI, Kadash-Edmondson KE, Xing Y"),
    "M5":  ("33590861", "10.1093/gigascience/giab008",
            "Danecek P, Bonfield JK, Liddle J, Marshall J, Ohan V, Pollard MO, "
            "Whitwham A, Keane T, McCarthy SA, Davies RM, Li H"),
    "M6":  ("20644199", "10.1101/gr.107524.110",
            "McKenna A, Hanna M, Banks E, Sivachenko A, Cibulskis K, Kernytsky A, "
            "Garimella K, Altshuler D, Gabriel S, Daly M, DePristo MA"),
    "M7":  ("20926424", "10.1093/bioinformatics/btq559",
            "Manichaikul A, Mychaleckyj JC, Rich SS, Daly K, Sale M, Chen WM"),
    "M10": ("21605215", "10.1111/j.1420-9101.2011.02297.x", "Zaykin DV"),
    "M11": ("37144459", "10.1093/nar/gkad347",
            "Kolberg L, Raudvere U, Kuzmin I, Adler P, Vilo J, Peterson H"),
    "I1":  (None, "10.1161/CIRCGENETICS.109.904698",
            "Kong SW, Hu YW, Ho JWK, Ikeda S, Polster S, John R, Hall JL, Bisping E, "
            "Pieske B, dos Remedios CG, Pu WT"),
    # ⚠ The PMID here was 29995277 — WHICH IS A DIFFERENT PAPER (esummary returns
    # Vimalanathan AK, Ehler E, Gehmlich K). It was never reported by any source: it
    # came from an illustrative example invented while writing the brief for the
    # verification pass, and was then copied in here as if it were data. Set to None
    # so the DOI resolves the real record. See RESULT_pubmed_initials_pass.md.
    "I2":  (None, "10.1007/s12551-018-0439-y", "Beqqali A"),
    "I3":  (None, "10.1007/s12265-013-9482-z",
            "Lara-Pezzi E, Gomez-Salinero J, Gatto A, Garcia-Pavia P"),
    "I4":  (None, "10.1016/j.jacc.2009.05.038",
            "Brauch KM, Karst ML, Herron KJ, de Andrade M, Pellikka PA, Rodeheffer RJ, "
            "Michels VV, Olson TM"),
    "I6":  ("42177204", "10.1038/s41467-026-73230-w",
            "Radke MH, Badillo Lisakowski V, Meinke S, Britto-Borges T, "
            "Schneider-Lunitz V, Hummel O, van Heesch S, Ruiz Orera J, Hubner N, "
            "Granzier H, Dieterich C, Gotthardt M"),
    "I7":  (None, "10.1038/s41467-020-20327-5", "[23 authors — see REFERENCES_draft.md]"),
    "I8":  ("41877661", "10.1093/cvr/cvag068",
            "Furtado M, Barbosa P, Wemans A, Zhang L, Lumley A, Carvalho T, Napoleao P, "
            "Leszek P, Carmo-Fonseca M, Devaux Y, Martins S"),
    "I9":  ("41902389", "10.1093/cvr/cvag070", "Wang Y, Tay WT, Gao C"),
    "I10": (None, "10.1007/s11886-026-02386-0",
            "Akhigbe F, Song N, Alagarsamy J, Li H, Gao C"),
    "I11": ("41017471", "10.1161/CIRCULATIONAHA.125.074959",
            "[18 authors — see REFERENCES_draft.md]"),
    "I12": (None, "10.1016/j.heliyon.2024.e35408",
            "Li J, Tu D, Li S, Guo Z, Song X"),
    "I13": (None, "10.1186/s12967-022-03694-z",
            "Britto-Borges T, Ludt A, Boileau E, Gjerga E, Marini F, Dieterich C"),
    "D1":  ("33517678", "10.1161/CIRCGEN.120.003142",
            "Liu CF, Ni Y, Moravec CS, Morley M, Ashley EA, Cappola TP, Margulies KB, "
            "Tang WHW"),
    "D1b": (None, "10.1016/j.ygeno.2014.12.002",
            "Liu Y, Morley M, Brandimarto J, Hannenhalli S, Hu Y, Ashley EA, Tang WHW, "
            "Moravec CS, Margulies KB, Cappola TP, Li M"),
    "D2":  ("30419824", "10.1186/s12864-018-5213-9",
            "Sweet ME, Cocciolo A, Slavov D, Jones KL, Sweet JR, Graw SL, Reece TB, "
            "Ambardekar AV, Bristow MR, Mestroni L, Taylor MRG"),
    "D3":  ("33558530", "10.1038/s41467-021-21146-y",
            "[15 authors — see REFERENCES_draft.md]"),
    "D4":  ("31197155", "10.1038/s41597-019-0094-6",
            "[12 authors — see REFERENCES_draft.md]"),
    "D5b": ("36130500", "10.1016/j.celrep.2022.111400",
            "[26 authors — see REFERENCES_draft.md]"),
    "D7":  ("35732739", "10.1038/s41586-022-04817-8",
            "[19 authors — see REFERENCES_draft.md]"),
    "D9":  (None, "10.3389/fbinf.2026.1934275",
            "Li J, He Y, Zhao X, Guo Z, Song J, Song X"),
    "D10": ("38288598", "10.1161/CIRCGEN.123.004355",
            "Gregorich ZR, Zhang Y, Kamp TJ, Granzier HL, Guo W"),
    "D11": ("24960161", "10.1172/JCI74523",
            "[17 authors — see REFERENCES_draft.md]"),
    "D12": ("33805770", "10.3390/ijms22062928",
            "Maimaiti R, Zhu C, Zhang Y, Ding Q, Guo W"),
    "D13": ("32276354", "10.3390/genes11040402",
            "Fochi S, Lorenzi P, Galasso M, Stefani C, Trabetti E, Zipeto D, "
            "Romanelli MG"),
    "D14": ("19075228", "10.1073/pnas.0809045105",
            "Kalsotra A, Xiao X, Ward AJ, Castle JC, Johnson JM, Burge CB, Cooper TA"),
    "D15": ("15988035", "10.1128/MCB.25.14.6267-6278.2005",
            "Ladd AN, Taffet G, Hartley C, Kearney DL, Cooper TA"),
    "D16": ("15652482", "10.1016/j.cell.2004.11.036",
            "[17 authors — see REFERENCES_draft.md]"),
    "D17": (None, "10.1093/nar/gkab702",
            "Zhao W, Zhang S, Zhu Y, Xi X, Bao P, Ma Z, Kapral TH, Chen S, Zagrovic B, "
            "Yang YT, Lu ZJ"),
    "L3":  ("40034311", "10.1016/j.heliyon.2025.e42499",
            "Conning-Rowland M, Cheng CW, Brown O, Giannoudi M, Levelt E, Roberts LD, "
            "Griffin KJ, Cubbon RM"),
}

# The two references this pass exists to settle, plus the house-style one.
DISPUTED = {
    "M6": "Crossref gives bare initials; OpenAlex gives Hanna MG, Kernytsky AM, "
          "Garimella KV, Gabriel SL, Daly MJ",
    "M7": "Crossref gives Manichaikul A / Sale M; OpenAlex gives Manichaikul AW / Sale MM",
    "M3": "two compound Spanish/Catalan surnames. ⚠ RESOLVED AGAINST THE EARLIER "
          "GUIDANCE: PubMed renders them 'Martinez JG' and 'Gomez LM', not "
          "'Gonzalez Martinez J' and 'Martinez Gomez L'. Decide from the title page.",
}


def get(url):
    req = urllib.request.Request(url, headers={"User-Agent": "refs-check/1.0"})
    with urllib.request.urlopen(req, timeout=30) as r:
        return r.read().decode("utf-8", "replace")


def pmid_from_doi(doi, email):
    q = urllib.parse.urlencode({"db": "pubmed", "term": f"{doi}[DOI]",
                                "retmode": "json", **({"email": email} if email else {})})
    try:
        d = json.loads(get(f"{EUTILS}/esearch.fcgi?{q}"))
        ids = d.get("esearchresult", {}).get("idlist", [])
        return ids[0] if ids else None
    except Exception:
        return None


def authors_from_pmid(pmid, email):
    q = urllib.parse.urlencode({"db": "pubmed", "id": pmid, "retmode": "json",
                                **({"email": email} if email else {})})
    d = json.loads(get(f"{EUTILS}/esummary.fcgi?{q}"))
    rec = d["result"][pmid]
    names = [a["name"] for a in rec.get("authors", [])
             if a.get("authtype") == "Author"]
    return names, rec.get("sorttitle") or rec.get("title", "")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", nargs="*", help="slots to check (default: all)")
    ap.add_argument("--email", default=None, help="your email, passed to NCBI")
    ap.add_argument("--sleep", type=float, default=0.4,
                    help="seconds between calls (NCBI asks for <=3/sec unkeyed)")
    a = ap.parse_args()

    slots = a.only or list(REFS)
    bad = [s for s in slots if s not in REFS]
    if bad:
        print(f"unknown slot(s): {', '.join(bad)}\nknown: {', '.join(REFS)}")
        sys.exit(2)

    print("=" * 78)
    print("99 — PubMed initials pass")
    print("=" * 78)
    print("\n⚠ PubMed is the conventional Vancouver source, NOT the final authority.")
    print("  NLM derives initials from the article's title page, so where they differ")
    print("  the TITLE PAGE wins. Treat DIFFERS below as 'open the PDF', not 'PubMed")
    print("  is right'.\n")

    agree = differ = missing = failed = 0
    notes = []

    for slot in slots:
        pmid, doi, drafted = REFS[slot]
        print("-" * 78)
        print(f"{slot}   doi:{doi}")
        if slot in DISPUTED:
            print(f"  ⚠ DISPUTED: {DISPUTED[slot]}")

        if not pmid:
            pmid = pmid_from_doi(doi, a.email)
            time.sleep(a.sleep)
            if pmid:
                print(f"  PMID resolved from DOI: {pmid}")
            else:
                print("  NO PMID — not indexed, or the DOI lookup found nothing.")
                print("  → settle from the article's title page.")
                missing += 1
                notes.append((slot, "no PMID"))
                continue

        try:
            names, title = authors_from_pmid(pmid, a.email)
            time.sleep(a.sleep)
        except Exception as e:
            print(f"  LOOKUP FAILED ({e.__class__.__name__}: {e})")
            print("  → retry, or settle from the title page.")
            failed += 1
            notes.append((slot, f"lookup failed: {e}"))
            continue

        pm = ", ".join(names)
        print(f"  PubMed ({len(names)}): {pm}")

        if drafted.startswith("["):
            print(f"  drafted: {drafted}")
            print("  → compare by hand against REFERENCES_draft.md (long list).")
            notes.append((slot, "long list — compare by hand"))
            continue

        want = [x.strip() for x in drafted.split(",")]
        if want == names:
            print("  AGREES with the drafted list.")
            agree += 1
            continue

        differ += 1
        print(f"  DIFFERS  (drafted {len(want)} authors, PubMed {len(names)})")
        for i in range(max(len(want), len(names))):
            w = want[i] if i < len(want) else "—"
            p = names[i] if i < len(names) else "—"
            if w != p:
                print(f"    pos {i + 1:>2}:  drafted {w:<28} PubMed {p}")
        notes.append((slot, "differs"))

    print("=" * 78)
    print(f"checked {len(slots)}   agrees {agree}   differs {differ}   "
          f"no PMID {missing}   failed {failed}")
    if notes:
        print("\nneeds attention:")
        for s, why in notes:
            print(f"  {s:<5} {why}")
    print("""
WHAT TO DO WITH 'DIFFERS'
  Open that paper's PDF and read the byline. Vancouver = surname + the initials of ALL
  given names as the title page prints them, so "Christine S. Moravec" is Moravec CS and
  "Michael Morley, MS" is Morley M (MS is a degree, not a name). Where a middle initial
  appears on the title page it belongs in the reference; where it does not, it does not.

  ⚠ COMPOUND SPANISH/CATALAN SURNAMES — THIS SCRIPT'S EARLIER ADVICE WAS WRONG.
  It said Vancouver treats both elements as the surname, so "Jose Gonzalez Martinez"
  would be "Gonzalez Martinez J". **PubMed renders it "Martinez JG"** — treating the
  first surname as a middle name. Both forms are in use; PubMed's is a known NLM
  handling of Spanish double surnames, and the journal's own byline is what a
  copy-editor checks against. **Decide from the title page, and say which you used.**
  Do not take either form as settled by style alone.
""")
    print("=" * 78)


if __name__ == "__main__":
    main()
