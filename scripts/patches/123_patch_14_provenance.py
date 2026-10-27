#!/usr/bin/env python3
"""
123_patch_14_provenance.py -- make 14 record WHICH g:Profiler it queried.

WHAT IS WRONG
-------------
`14_run_gprofiler.py` posts to

    https://biit.cs.ut.ee/gprofiler/api/gost/profile/

and writes the terms. It records **nothing** about the database it hit.

g:Profiler's annotation sources update several times a year, so a rerun does not
reproduce an earlier result even on a byte-identical gene list. An enrichment
figure with no version is not reproducible, and the gene list alone does not
make it so. Verified live from the cluster, 16 Oct 2026:

    curl -s https://biit.cs.ut.ee/gprofiler/api/util/data_versions/
    {"biomart":"Ensembl","biomart_version":"114","genebuild":"GRCh38.p14",
     "gprofiler_version":"e114_eg62_p19_9dc92ac8","organism":"hsapiens",
     "sources":{"CORUM":{"version":"28.11.2022 Corum 4.1"},
                "GO:BP":{"version":"annotations: BioMart
                                    classes: releases/2026-01-23"}, ...}}

**The GO classes are a January 2026 release.** If the original Figure 6
enrichment predates it, the terms differ for reasons that have nothing to do
with read length. That is a fact to record, not a defect to chase.

WHAT THIS PATCH DOES
--------------------
After a successful query, and before the results are written, it fetches the
`util/data_versions/` endpoint for the same organism and deposits

    <tag>_gprofiler_PROVENANCE.txt

carrying the access timestamp, both endpoints, the organism, the sources
requested, the domain scope, the foreground and background sizes, any
version-bearing key found in the response's own `meta`, and the full
`data_versions` JSON.

It also prints a **draft Methods sentence** with the captured values filled in,
because the point of the file is that a sentence can cite it.

Design choices, and why:

  * The version URL is derived from `URL` (`URL.split("/api/")[0]`), not written
    as a second literal. One source of truth for the host.
  * If the version fetch fails, the enrichment is still written — a transient
    network failure should not destroy a successful query — but the provenance
    file's **first line** says VERSION NOT CAPTURED and a warning block is
    printed. A reader of the file cannot miss it.
  * The payload is reported by iterating its keys, so this patch does not need
    to know g:Profiler's field names. List values are reported as counts; other
    values are truncated.
  * `import time` is local to the block, so no module-level import is added to
    a file this patch cannot see in full.

USAGE
-----
  python 123_patch_14_provenance.py
  python 123_patch_14_provenance.py --apply
  python 123_patch_14_provenance.py --self-test
"""
import argparse
import difflib
import os
import shutil
import sys
import time

TARGET = "14_run_gprofiler.py"

ANCHOR = '    rows = res.get("result", [])\n'

BLOCK = '''    # ---------------------------------------------------------- provenance
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
    _prov = f"{a.tag}_gprofiler_PROVENANCE.txt"
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
            _pf.write("VERSION CAPTURED\\n")
        else:
            _pf.write("** VERSION NOT CAPTURED — the Methods CANNOT cite a "
                      "g:Profiler version from this run **\\n")
            _pf.write(f"   reason: {_verr}\\n")
        _pf.write(f"accessed (UTC)     "
                  f"{_time.strftime('%Y-%m-%dT%H:%M:%SZ', _time.gmtime())}\\n")
        _pf.write(f"query endpoint     {URL}\\n")
        _pf.write(f"version endpoint   {_vurl}\\n")
        _pf.write(f"organism           {a.organism}\\n")
        _pf.write(f"sources requested  {a.sources}\\n")
        for _k in sorted(payload):
            _v = payload[_k]
            if isinstance(_v, (list, tuple)):
                _pf.write(f"payload[{_k}]  {len(_v)} items\\n")
            else:
                _pf.write(f"payload[{_k}]  {str(_v)[:120]}\\n")
        if _meta_ver:
            _pf.write(f"meta version keys  {json.dumps(_meta_ver)[:400]}\\n")
        if _ver:
            _pf.write("\\n--- data_versions ---\\n")
            _pf.write(json.dumps(_ver, indent=2, sort_keys=True) + "\\n")
    print(f"[prov] wrote {_prov}", file=sys.stderr)

    if _ver:
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

'''

EDITS = [("the provenance capture, before the results are read",
          ANCHOR, BLOCK + ANCHOR, 1)]


def apply_edits(text):
    problems, out = [], text
    for name, anchor, repl, want in EDITS:
        got = out.count(anchor)
        if got != want:
            problems.append(f"{name}: anchor occurs {got}x, expected {want}x\n"
                            f"      anchor: {anchor.strip()[:80]}")
            continue
        if "_gprofiler_PROVENANCE" in out:
            problems.append("already patched (_gprofiler_PROVENANCE present)")
            continue
        out = out.replace(anchor, repl, want)
    return out, problems


def verify(text):
    bad = []
    code = "\n".join(l for l in text.splitlines()
                     if not l.lstrip().startswith("#"))
    if "_gprofiler_PROVENANCE.txt" not in code:
        bad.append("the provenance filename is missing")
    if 'URL.split("/api/")[0]' not in code:
        bad.append("the version URL is not derived from URL")
    if "VERSION NOT CAPTURED" not in code:
        bad.append("the failure path does not mark the file")
    if "DRAFT METHODS SENTENCE" not in code:
        bad.append("the draft sentence is missing")
    # the block must sit ABOVE the results read
    i_b = text.find("_prov = f\"{a.tag}_gprofiler_PROVENANCE.txt\"")
    i_r = text.find('rows = res.get("result", [])')
    if i_b < 0 or i_r < 0 or i_b > i_r:
        bad.append("the provenance block is not above the results read")
    # a second hard-coded host would defeat the single-source-of-truth point
    if code.count("biit.cs.ut.ee/gprofiler/api/util") > 0:
        bad.append("the util endpoint was written as a literal URL in code")
    return bad


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("target", nargs="?", default=TARGET)
    ap.add_argument("--apply", action="store_true")
    ap.add_argument("--self-test", action="store_true")
    a = ap.parse_args()
    if a.self_test:
        return self_test()

    self_test()
    print()
    print("=" * 78)
    print("123 — record which g:Profiler was queried")
    print("=" * 78)
    if not os.path.exists(a.target):
        print(f"[fatal] {a.target} not found.", file=sys.stderr)
        return 2
    src = open(a.target, errors="replace").read()
    for need in ("import json", "import urllib.request", "URL ="):
        if need not in src:
            print(f"[fatal] {a.target} has no `{need}`; the block needs it.",
                  file=sys.stderr)
            return 2

    new, problems = apply_edits(src)
    if problems:
        print("REFUSING — anchors do not match:")
        for p in problems:
            print(f"  ✗ {p}")
        print("\n  Nothing written. Send me lines 95-115 and I will re-anchor.")
        return 3
    bad = verify(new)
    if bad:
        print("REFUSING — patched text fails verification:")
        for b in bad:
            print(f"  ✗ {b}")
        return 3

    print("".join(difflib.unified_diff(
        src.splitlines(True), new.splitlines(True),
        fromfile=a.target, tofile=a.target + " (patched)", n=2)))
    print("=" * 78)

    if not a.apply:
        print("DRY RUN — nothing written. Re-run with --apply.")
        print()
        print("Then, ON THE LOGIN NODE (compute nodes have no network):")
        print()
        print("  python3 14_run_gprofiler.py --tag compadj_rl125")
        print("  python3 14_run_gprofiler.py --tag compadj_rl125 --whole-genome")
        print()
        print("  # and the 50 bp comparator, so the two are provenanced alike")
        print("  python3 14_run_gprofiler.py --tag compadj_50check")
        print()
        print("Each writes <tag>_gprofiler_PROVENANCE.txt and prints a draft")
        print("Methods sentence. Expected version today: e114_eg62_p19_9dc92ac8,")
        print("Ensembl 114, GRCh38.p14, GO classes releases/2026-01-23.")
        print()
        print("⚠ The 50 bp and 125 bp enrichments MUST be run on the same day,")
        print("  against the same version. Comparing a 2026 enrichment of the")
        print("  125 bp list against whatever the original 927 run hit is a")
        print("  cross-version comparison and is not interpretable.")
        return 0

    bak = f"{a.target}.bak.{time.strftime('%Y%m%d_%H%M%S')}"
    shutil.copy2(a.target, bak)
    open(a.target, "w").write(new)
    print(f"[write] {a.target}   (backup: {bak})")
    disk = open(a.target, errors="replace").read()
    bad = verify(disk)
    if bad or disk != new:
        print("✗ ON-DISK VERIFICATION FAILED — restoring:")
        for b in bad:
            print(f"    {b}")
        shutil.copy2(bak, a.target)
        return 3
    print("✓ verified on disk.")
    return 0


def self_test():
    import json as _json
    import subprocess
    import tempfile

    head = (
        "import io\n"
        "import json\n"
        "import sys\n"
        "import urllib.request\n"
        'URL = "https://biit.cs.ut.ee/gprofiler/api/gost/profile/"\n'
        "class A:\n"
        '    tag = "t"\n'
        '    organism = "hsapiens"\n'
        '    sources = "GO:BP,GO:MF"\n'
        "a = A()\n"
        'payload = {"query": ["A", "B", "C"],\n'
        '           "background": ["A", "B", "C", "D", "E"],\n'
        '           "organism": "hsapiens",\n'
        '           "domain_scope": "custom_annotated",\n'
        '           "significance_threshold_method": "g_SCS"}\n'
        'res = {"result": [{"x": 1}], "meta": {"genes_metadata": {}}}\n'
        "CANNED = " + _json.dumps({
            "biomart": "Ensembl", "biomart_version": "114",
            "genebuild": "GRCh38.p14",
            "gprofiler_version": "e114_eg62_p19_9dc92ac8",
            "organism": "hsapiens",
            "sources": {"GO:BP": {"version": "releases/2026-01-23"}}}) + "\n"
        "FAIL = len(sys.argv) > 1\n"
        "def _fake(req, timeout=0):\n"
        "    if FAIL:\n"
        "        raise OSError('network down')\n"
        "    return io.BytesIO(json.dumps(CANNED).encode())\n"
        "urllib.request.urlopen = _fake\n"
        "def run():\n")
    tail = "    return rows\n\nif __name__ == '__main__':\n    run()\n"

    new, probs = apply_edits(head + ANCHOR + tail)
    assert not probs, probs
    assert not verify(new), verify(new)

    d = tempfile.mkdtemp()
    p = os.path.join(d, "fx.py")
    open(p, "w").write(new)

    # --- success path
    r = subprocess.run([sys.executable, "-I", p], capture_output=True,
                       text=True, cwd=d)
    assert r.returncode == 0, r.stdout + r.stderr
    prov = os.path.join(d, "t_gprofiler_PROVENANCE.txt")
    assert os.path.exists(prov), os.listdir(d)
    txt = open(prov).read()
    assert txt.startswith("VERSION CAPTURED"), txt[:200]
    assert "e114_eg62_p19_9dc92ac8" in txt, txt
    assert "GRCh38.p14" in txt, txt
    assert "releases/2026-01-23" in txt, txt
    assert "payload[background]  5 items" in txt, txt
    assert "payload[query]  3 items" in txt, txt
    assert "/api/util/data_versions/" in txt, txt
    assert "DRAFT METHODS SENTENCE" in r.stderr, r.stderr
    assert "5 genes" in r.stderr, r.stderr      # background count in the draft
    os.remove(prov)

    # --- failure path: results still produced, file marked, warning printed
    r2 = subprocess.run([sys.executable, "-I", p, "fail"],
                        capture_output=True, text=True, cwd=d)
    assert r2.returncode == 0, r2.stdout + r2.stderr
    txt2 = open(prov).read()
    assert txt2.startswith("** VERSION NOT CAPTURED"), txt2[:200]
    assert "network down" in txt2, txt2
    assert "UNPROVENANCED" in r2.stderr, r2.stderr
    assert "DRAFT METHODS SENTENCE" not in r2.stderr, \
        "a draft sentence must NOT be offered without a version"

    # --- double application refused
    _n, p3 = apply_edits(new)
    assert p3, "second application was not refused"

    print("[ok] self-test passed: the version is captured and deposited with "
          "the access time, both endpoints and the payload sizes; the draft "
          "Methods sentence carries the real version and background count; on "
          "a network failure the results still run but the file is marked "
          "NOT CAPTURED and no draft sentence is offered; the util URL is "
          "derived from URL, not a second literal; double application refused")
    return 0


if __name__ == "__main__":
    sys.exit(main())
