#!/usr/bin/env python3
"""
14_run_gprofiler.py -- run g:Profiler with a PROPER custom background, from the
command line, and write a TSV you can put straight in the supplement.

Avoids pasting 12,725 gene names into a browser textarea, and makes the
enrichment reproducible (the exact query is recorded in the output).

Runs all your lists as a multi-query in ONE call, so the terms are directly
comparable across them:
    all        <tag>_foreground.txt
    up         <tag>_foreground_up.txt
    down       <tag>_foreground_down.txt
    multi      <tag>_multi_event_type.txt

Background: <tag>_background.txt, with domain_scope=custom_annotated -- the API
equivalent of the web UI's "Custom over annotated genes".

Needs network. Run it on the LOGIN NODE, not a compute node.

Usage:
  python 14_run_gprofiler.py --tag compadj
  python 14_run_gprofiler.py --tag compadj --whole-genome   # the artefact check
"""
import argparse
import json
import os
import sys
import urllib.request

URL = "https://biit.cs.ut.ee/gprofiler/api/gost/profile/"
SOURCES = ["GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC", "WP"]


def read_list(path):
    with open(path) as fh:
        return [l.strip() for l in fh if l.strip()]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="compadj")
    ap.add_argument("--organism", default="hsapiens")
    ap.add_argument("--sources", default=",".join(SOURCES))
    ap.add_argument("--with-genes", action="store_true",
                    help="also return WHICH genes fall in each term, as an extra "
                         "'intersections' column -- needed to name the genes in "
                         "the paper")
    ap.add_argument("--whole-genome", action="store_true",
                    help="DELIBERATELY use the wrong background, to demonstrate "
                         "the testability artefact for the methods section")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()

    queries = {}
    for name, fn in [("all", f"{a.tag}_foreground.txt"),
                     ("up", f"{a.tag}_foreground_up.txt"),
                     ("down", f"{a.tag}_foreground_down.txt"),
                     ("multi_event_type", f"{a.tag}_multi_event_type.txt")]:
        if os.path.exists(fn):
            g = read_list(fn)
            if g:
                queries[name] = g
                print(f"[info] {name:<18} {len(g):>5} genes  ({fn})", file=sys.stderr)
    if not queries:
        sys.exit(f"[fatal] no {a.tag}_foreground*.txt files found")

    payload = {
        "organism": a.organism,
        "query": queries,
        "sources": [s.strip() for s in a.sources.split(",") if s.strip()],
        "user_threshold": 0.05,
        "significance_threshold_method": "g_SCS",
        "all_results": False,
        "ordered": False,
        "no_evidences": not a.with_genes,
    }

    if a.whole_genome:
        payload["domain_scope"] = "annotated"
        tagout = a.out or f"{a.tag}_gprofiler_WHOLEGENOME.tsv"
        print("[warn] using the WHOLE-GENOME background on purpose. These results "
              "are the artefact demonstration, NOT your result.", file=sys.stderr)
    else:
        bg_file = f"{a.tag}_background.txt"
        if not os.path.exists(bg_file):
            sys.exit(f"[fatal] {bg_file} missing -- run 13_gene_lists_for_enrichment.py")
        bg = read_list(bg_file)
        print(f"[info] {'background':<18} {len(bg):>5} genes  ({bg_file})", file=sys.stderr)
        payload["domain_scope"] = "custom_annotated"
        payload["background"] = bg
        tagout = a.out or f"{a.tag}_gprofiler.tsv"

    req = urllib.request.Request(
        URL, data=json.dumps(payload).encode(),
        headers={"Content-Type": "application/json",
                 "User-Agent": "cardiac-splicing-analysis"})
    print("[info] querying g:Profiler ...", file=sys.stderr)
    try:
        with urllib.request.urlopen(req, timeout=300) as r:
            res = json.load(r)
    except Exception as e:
        sys.exit(f"[fatal] request failed: {e}\n"
                 "        Run this on the LOGIN NODE (compute nodes have no network).\n"
                 "        If it still fails, paste the lists into the web UI instead:\n"
                 "        https://biit.cs.ut.ee/gprofiler/gost")

    # ---------------------------------------------------------- provenance
    # g:Profiler's annotation sources update several times a year, so a rerun
    # does NOT reproduce an earlier result even on an identical gene list. The
    # gost/profile response does not carry a version, so fetch it from the
    # util endpoint and deposit it beside the results.
    #
    # Verified live 16 Oct 2026:
    #   https://biit.cs.ut.ee/gprofiler/api/util/data_versions/
    #     gprofiler_version  e114_eg62_p19_9dc92ac8
    #     biomart            Ensembl 114 · genebuild GRCh38.p14
    #     GO classes         releases/2026-01-23
    #
    # The version URL is DERIVED from URL so the host is named once.
    import time as _time
    # The results file is <tag>_gprofiler_WHOLEGENOME.tsv in whole-genome
    # mode, so the provenance file must carry the same distinction.
    # Without it the artefact run OVERWRITES the real run's provenance
    # and the deposited file describes a query nobody will report.
    _wg = bool(getattr(a, "whole_genome", False))
    _prov = (f"{a.tag}_gprofiler"
             f"{'_WHOLEGENOME' if _wg else ''}_PROVENANCE.txt")
    _vurl = URL.split("/api/")[0] + "/api/util/data_versions/"
    _ver, _verr = None, None
    try:
        _vreq = urllib.request.Request(
            _vurl + "?organism=" + str(a.organism),
            headers={"User-Agent": "14_run_gprofiler"})
        with urllib.request.urlopen(_vreq, timeout=60) as _vr:
            _ver = json.load(_vr)
    except Exception as _e:
        _verr = repr(_e)

    _meta_ver = {_k: _v for _k, _v in (res.get("meta") or {}).items()
                 if "version" in str(_k).lower()}

    with open(_prov, "w") as _pf:
        if _ver:
            _pf.write("VERSION CAPTURED\n")
        else:
            _pf.write("** VERSION NOT CAPTURED — the Methods CANNOT cite a "
                      "g:Profiler version from this run **\n")
            _pf.write(f"   reason: {_verr}\n")
        _pf.write(f"accessed (UTC)     "
                  f"{_time.strftime('%Y-%m-%dT%H:%M:%SZ', _time.gmtime())}\n")
        _pf.write(f"query endpoint     {URL}\n")
        _pf.write(f"version endpoint   {_vurl}\n")
        _pf.write(f"organism           {a.organism}\n")
        _pf.write(f"sources requested  {a.sources}\n")
        _pf.write(f"whole_genome       {_wg}\n")
        _pf.write("background mode    "
                  + ("WHOLE GENOME — no custom background; this run is "
                     "the artefact\n                   demonstration, "
                     "NOT the reported result\n"
                     if _wg else
                     "custom (genes with >=1 tested event)\n"))
        for _k in sorted(payload):
            _v = payload[_k]
            if isinstance(_v, (list, tuple)):
                _pf.write(f"payload[{_k}]  {len(_v)} items\n")
            else:
                _pf.write(f"payload[{_k}]  {str(_v)[:120]}\n")
        if _meta_ver:
            _pf.write(f"meta version keys  {json.dumps(_meta_ver)[:400]}\n")
        if _ver:
            _pf.write("\n--- data_versions ---\n")
            _pf.write(json.dumps(_ver, indent=2, sort_keys=True) + "\n")
    print(f"[prov] wrote {_prov}", file=sys.stderr)

    if _ver and _wg:
        print("", file=sys.stderr)
        print("  ** WHOLE-GENOME RUN — the ARTEFACT DEMONSTRATION, not "
              "your result.", file=sys.stderr)
        print("     No Methods sentence is offered: there is NO custom "
              "background here,", file=sys.stderr)
        print("     so describing one would be false. Report this run as "
              "the", file=sys.stderr)
        print("     background-sensitivity check, and cite the version "
              "from the", file=sys.stderr)
        print("     proper-background run.", file=sys.stderr)
        print("", file=sys.stderr)
    elif _ver:
        _gv = _ver.get("gprofiler_version", "?")
        _bm = _ver.get("biomart", "?")
        _bv = _ver.get("biomart_version", "?")
        _gb = _ver.get("genebuild", "?")
        _nbg = len(payload.get("background") or [])
        _date = _time.strftime("%d %B %Y", _time.gmtime())
        print("", file=sys.stderr)
        print("  DRAFT METHODS SENTENCE — check the background wording against",
              file=sys.stderr)
        print("  how 13 defines it (genes with at least one TESTED event):",
              file=sys.stderr)
        print(f"    Functional enrichment was performed with g:Profiler "
              f"(version {_gv};", file=sys.stderr)
        print(f"    {_bm} {_bv}, {_gb}; accessed {_date}), with "
              f"g:SCS", file=sys.stderr)
        print(f"    multiple-testing correction and a custom background of "
              f"{_nbg:,} genes", file=sys.stderr)
        print(f"    carrying at least one event that passed the coverage and "
              f"variance filters.", file=sys.stderr)
        print("", file=sys.stderr)
    else:
        print("", file=sys.stderr)
        print("  ** NO VERSION CAPTURED. The results below are usable but "
              "UNPROVENANCED.", file=sys.stderr)
        print("     Do not cite a g:Profiler version in the Methods on the "
              "basis of this run.", file=sys.stderr)
        print(f"     Retry the endpoint by hand:  curl -s {_vurl}",
              file=sys.stderr)
        print("", file=sys.stderr)

    rows = res.get("result", [])
    if not rows:
        print("[warn] no significant terms returned", file=sys.stderr)

    # ---- gene mapping for --with-genes --------------------------------
    # CRITICAL: `intersections` is aligned to the MAPPED Ensembl gene list at
    # meta.genes_metadata.query[<query>].ensgs -- NOT to the symbols we
    # submitted. Indexing the submitted list produces silently wrong genes
    # (correct COUNT, wrong IDENTITY), which is undetectable without this
    # mapping. Build ensg -> input symbol from the mapping g:Profiler returns.
    ensg_order, rev_map = {}, {}
    if a.with_genes:
        gm = (res.get("meta", {}).get("genes_metadata", {}) or {}).get("query", {})
        for qname, qmeta in gm.items():
            ensg_order[qname] = qmeta.get("ensgs", []) or []
            rmap = {}
            for sym, ensgs in (qmeta.get("mapping", {}) or {}).items():
                for e in (ensgs if isinstance(ensgs, list) else [ensgs]):
                    rmap[e] = sym
            rev_map[qname] = rmap
        if not ensg_order:
            print("[warn] --with-genes: meta.genes_metadata.query missing from the "
                  "response; gene column will be left empty rather than guessed.",
                  file=sys.stderr)

    def genes_for(r):
        qn = r.get("query")
        inter = r.get("intersections") or []
        ensgs = ensg_order.get(qn, [])
        if not ensgs or len(ensgs) != len(inter):
            return "", False
        rmap = rev_map.get(qn, {})
        hits = [rmap.get(ensgs[i], ensgs[i]) for i, ev in enumerate(inter) if ev]
        return ",".join(hits), (len(hits) == r.get("intersection_size"))

    cols = ["query", "source", "native", "name", "p_value", "term_size",
            "query_size", "intersection_size", "precision", "recall"]
    if a.with_genes:
        cols.append("genes")
    n_bad = 0
    with open(tagout, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            gtxt, ok = ("", True)
            if a.with_genes:
                gtxt, ok = genes_for(r)
                if not ok:
                    n_bad += 1
            vals = [gtxt if c == "genes" else str(r.get(c, "")) for c in cols]
            fh.write("\t".join(vals) + "\n")
    if a.with_genes and n_bad:
        print(f"[warn] {n_bad} terms: extracted gene count != intersection_size. "
              "Do NOT trust the gene column for those rows.", file=sys.stderr)
    elif a.with_genes and rows:
        print("[info] gene column validated: count matches intersection_size for "
              "every term", file=sys.stderr)
    print(f"[write] {tagout}  ({len(rows)} significant terms)", file=sys.stderr)

    # readable summary
    print()
    for q in queries:
        sub = [r for r in rows if r.get("query") == q]
        sub.sort(key=lambda r: r["p_value"])
        print(f"=== {q}  ({len(sub)} significant terms) ===")
        print(f"  {'source':<7}{'p':>11}  {'n':>5}  term")
        for r in sub[:15]:
            print(f"  {r['source']:<7}{r['p_value']:>11.2e}  "
                  f"{r['intersection_size']:>5}  {r['name'][:70]}")
        if not sub:
            print("  (none)")
        print()

    print("REMEMBER: if you also ran --whole-genome, compare the two. Terms that "
          "appear ONLY with the whole-genome background are the testability "
          "artefact, not biology.")


if __name__ == "__main__":
    main()
