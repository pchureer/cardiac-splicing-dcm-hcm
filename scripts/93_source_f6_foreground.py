#!/usr/bin/env python3
"""
93_source_f6_foreground.py  (v2)  —  source the 927 on Figure 6.

WHAT v1 GOT WRONG — both corrections matter
-------------------------------------------
1. **v1's recompute was broken.** It globbed for event tables and kept the FIRST match
   per event class, so "the A3SS table" became `A3SS_Colorado_DCMvsICM`. The roll-up it
   printed (4 SE genes) was meaningless. v2 names the five tables explicitly and reads
   nothing else.

2. **v1 claimed the hard-coded numbers contradict each other. They do not.**
   `FG_ALL=927, FG_UP=636, FG_DOWN=291, FG_MIXED=75` is internally consistent if `up`
   and `down` PARTITION the 927 — each gene assigned to one arm — with `MIXED` merely
   recording how many of them carry events in both directions. 636 + 291 = 927 is then
   correct, not a double-count. v1 asserted a contradiction that is not there.

3. **`query_size` is not the foreground.** g:Profiler reports it per data source, as the
   number of query genes annotated in THAT source's domain, which is why one query shows
   several values (`all`: 357, 788, 799). It is always smaller than the submitted list.
   It is still useful as a consistency check, and it checks out: all four queries map at
   83–89% of the hard-coded sizes.

WHAT REMAINS, AND WHAT THIS SCRIPT DOES
---------------------------------------
927 is hard-coded with no file behind it, and it is printed twice on a published figure
one digit from the Results headline 926. Those are different quantities:

  926  skipped exons significant under the BASE model (no composition term)
  927  GENES, almost certainly from the COMPOSITION-ADJUSTED model — the enrichment
       input file is named `compadj_gprofiler.tsv`, and the Methods' inverted scan
       already describes "134 genes significant across multiple event classes AFTER
       composition adjustment", which is panel C's 134.

So this reads the five `*_v2compadj` tables, counts distinct significant genes with the
up/down split computed BOTH ways, and tests 927 / 636 / 291 / 75 against each. The base
model is recomputed alongside, so the two are never confused again.

    python 93_source_f6_foreground.py --list        # show which tables exist
    python 93_source_f6_foreground.py               # compadj (default) + base
    python 93_source_f6_foreground.py --model base  # base only

GATES
  G1  every named table exists and resolves its columns. Missing one is a refusal, not
      a silent skip — a four-class roll-up is not a five-class answer.
  G2  per class: exon-level collapse, FDR < 0.05 and |dPSI| >= 0.10, distinct genes.
  G3  the up/down split under both readings — partition (each gene to one arm) and
      union (a gene may be in both) — reported side by side.
  G4  test against the hard-coded 927/636/291/75 and name what reproduces.
"""

import argparse
import os
import sys

import pandas as pd

HARD = dict(ALL=927, UP=636, DOWN=291, MIXED=75, MULTI=134)

CLASSES = ["SE", "MXE", "A3SS", "A5SS", "RI"]
KEYFIELDS = {"SE": 3, "RI": 3, "MXE": 5, "A3SS": 5, "A5SS": 5}

MODELS = {
    # the composition-adjusted model — what `compadj_gprofiler.tsv` was built from
    "compadj": "{cls}_v2compadj.all_events.tsv",
    # the base model — the 926 / 1,445 the Results headline reports
    "base": "{cls}_v2filtered.all_events.tsv",
    # the restricted cardiomyocyte-fraction score, for completeness
    "cmfrac": "{cls}_v2cmfrac.all_events.tsv",
}

# CALIBRATION. Before this script's verdict on an UNKNOWN number (927) means anything,
# it has to reproduce a KNOWN one. These are the recorded skipped-exon counts for each
# model; if the script's SE count does not land on one of them, the script is reading
# the tables differently from however S1 was built, and nothing it says about the
# foreground can be trusted until that is resolved.
KNOWN_SE = {
    "base": (926, "base model, no composition term — the Results headline"),
    "cmfrac": (1054, "restricted cardiomyocyte-fraction score — Results primary"),
    "compadj": (899, "nine-marker score — Results conservative bound"),
}

C_GENE = ["gene", "GeneSymbol", "geneSymbol", "gene_symbol", "Gene", "symbol"]
C_FDR = ["FDR", "fdr", "padj", "qvalue"]
C_DPSI = ["dPSI_adj", "dPSI", "dPSI_raw", "dpsi"]
C_KEY = ["key", "event_key", "ID"]


def hr(c="-"):
    print(c * 78)


def refuse(msg):
    print()
    hr("=")
    print("REFUSING — " + msg)
    hr("=")
    sys.exit(2)


def pick(df, cands, what, where):
    for c in cands:
        if c in df.columns:
            return c
    refuse(f"no {what} column in {where}\n  tried: {', '.join(cands)}\n"
           f"  has:   {', '.join(map(str, df.columns))}")


def one_model(d, model, quiet=False):
    """Return dict(genes, up_any, down_any, both, up_part, down_part, events, per_class)."""
    pat = MODELS[model]
    print()
    hr("=")
    print(f"MODEL `{model}`   ({pat.format(cls='<CLASS>')})")
    hr("=")

    missing = [c for c in CLASSES if not os.path.exists(os.path.join(d, pat.format(cls=c)))]
    if missing:
        refuse(f"model `{model}` is missing {len(missing)} of 5 tables: "
               f"{', '.join(pat.format(cls=c) for c in missing)}\n"
               f"  A four-class roll-up is not a five-class answer. Either the model "
               f"name is\n  wrong or those runs do not exist; resolve before trusting "
               f"any total.")

    per_class, up_ev, down_ev = {}, {}, {}
    for cls in CLASSES:
        f = os.path.join(d, pat.format(cls=cls))
        df = pd.read_csv(f, sep="\t", low_memory=False)
        g = pick(df, C_GENE, "gene", os.path.basename(f))
        q = pick(df, C_FDR, "FDR", os.path.basename(f))
        k = pick(df, C_KEY, "key", os.path.basename(f))
        dp = pick(df, C_DPSI, "dPSI", os.path.basename(f))
        n = KEYFIELDS[cls]

        fdr = pd.to_numeric(df[q], errors="coerce")
        dps = pd.to_numeric(df[dp], errors="coerce")
        sub = df[(fdr < 0.05) & (dps.abs() >= 0.10)].copy()
        sub["_fdr"] = fdr[sub.index]
        sub["_dpsi"] = dps[sub.index]
        sub["ekey"] = sub[k].astype(str).str.split(":").str[:n].str.join(":")
        sub = sub.sort_values("_fdr").drop_duplicates("ekey")

        genes = set(sub[g].dropna().astype(str))
        per_class[cls] = dict(events=len(sub), genes=genes, tested=len(df),
                              col_dpsi=dp, col_gene=g)
        for _, r in sub.iterrows():
            tgt = up_ev if r["_dpsi"] > 0 else down_ev
            tgt[str(r[g])] = tgt.get(str(r[g]), 0) + 1
        print(f"  {cls:5s} {pat.format(cls=cls):34s} "
              f"{len(df):>7,} rows -> {len(sub):>5,} sig, {len(genes):>5,} genes "
              f"[dPSI={dp}]")

    allg = set().union(*[v["genes"] for v in per_class.values()])
    up_any = set(up_ev)
    down_any = set(down_ev)
    both = up_any & down_any
    # partition reading: a mixed gene goes to whichever direction has more events
    up_part = {g for g in allg
               if up_ev.get(g, 0) >= down_ev.get(g, 0) and up_ev.get(g, 0) > 0}
    down_part = allg - up_part
    multi = {g for g in allg
             if sum(1 for c in CLASSES if g in per_class[c]["genes"]) >= 2}

    # ------------------------------------------------------------- calibration
    se_got = per_class["SE"]["events"]
    print()
    hr()
    print("  CALIBRATION — does this reading reproduce a number already in the record?")
    want, why = KNOWN_SE.get(model, (None, ""))
    anymatch = [(m, v[0], v[1]) for m, v in KNOWN_SE.items() if v[0] == se_got]
    print(f"    SE significant events, this reading   {se_got:,}")
    if want is not None:
        print(f"    recorded for model `{model}`           {want:,}   ({why})")
    if want is not None and se_got == want:
        print("    CALIBRATED — this reading reproduces the recorded count. Totals below "
              "can be trusted.")
    elif anymatch:
        m2, v2, w2 = anymatch[0]
        print(f"    ⚠ matches the count recorded for `{m2}` instead ({v2:,} — {w2}).")
        print("      The model-to-filename mapping in this script is wrong. Fix MODELS "
              "before\n      reading anything below.")
    else:
        print(f"    ⚠ NOT CALIBRATED — reproduces no recorded count (926 / 1,054 / 899).")
        print("""      This script is reading the tables differently from however S1 was
      built, and NOTHING below is evidence about the foreground until that is
      resolved. The likely causes, cheapest first:
        - the exon-level key: this script collapses on key.split(":")[:3] for SE.
          Check one row -- `head -2 SE_v2compadj.all_events.tsv | cut -f8` -- and
          confirm the first three fields really are the event's own coordinates.
        - an extra screen S1 applied and this does not (a coverage mask, a
          reliability flag, a PSI-variance filter already baked into one table
          and not the other).
        - a different dPSI column: this used the one named in the per-class line
          above.""")

    ev = sum(v["events"] for v in per_class.values())
    print()
    hr()
    print(f"  total significant events, five classes   {ev:,}")
    print(f"  distinct genes                           {len(allg):,}")
    print(f"  genes significant in >=2 event types     {len(multi):,}")
    print()
    print("  up/down, UNION reading (a gene may appear in both)")
    print(f"    up {len(up_any):,} · down {len(down_any):,} · both {len(both):,} "
          f"-> union {len(up_any) + len(down_any) - len(both):,}")
    print("  up/down, PARTITION reading (each gene to its majority direction)")
    print(f"    up {len(up_part):,} · down {len(down_part):,} "
          f"-> sum {len(up_part) + len(down_part):,}")
    return dict(model=model, genes=allg, up_any=up_any, down_any=down_any, both=both,
                up_part=up_part, down_part=down_part, multi=multi, events=ev,
                per_class=per_class, se_events=se_got,
                calibrated=(want is not None and se_got == want))


def verdict(results):
    print()
    hr("=")
    print("G4  does anything reproduce the hard-coded numbers?")
    hr("=")
    print(f"  41_figure6_enrichment.py lines 74-75:")
    print(f"    FG_ALL={HARD['ALL']}  FG_UP={HARD['UP']}  FG_DOWN={HARD['DOWN']}  "
          f"FG_MIXED={HARD['MIXED']}")
    print(f"  and panel C draws {HARD['MULTI']} genes significant in >=2 event types.")
    print(f"  Note {HARD['UP']} + {HARD['DOWN']} = {HARD['UP'] + HARD['DOWN']} = "
          f"FG_ALL, so up/down is a PARTITION and FG_MIXED is a flag, not a third "
          f"class.\n  That is self-consistent; the question is only which model and "
          f"which file.")

    if not any(r["calibrated"] for r in results):
        print()
        hr("=")
        print("  ⚠ VERDICT WITHHELD — no model reproduced its own recorded SE count")
        hr("=")
        print("""
  Every reading below disagrees with the record on a number that is already settled
  (926 / 1,054 / 899). A reading that cannot reproduce a known number is not evidence
  about an unknown one, so the comparison that follows is printed for diagnosis only
  and must NOT be used to conclude anything about 927.

  Resolve the calibration first. The fastest route to the real answer is not this
  script at all -- it is the list that was actually submitted:

      grep -n "query\\|gene_list\\|foreground\\|to_csv\\|write\\|gprofiler" \\
           14_run_gprofiler.py
      ls -la *gene*list* *query*.txt *foreground* *_genes.txt 2>/dev/null

  That script built the 927. Whatever filter it applied is the definition, and the
  figure's constants should be read from its output rather than recomputed.""")

    hits = []
    for r in results:
        m = r["model"]
        print()
        hr()
        print(f"  model `{m}`" + ("" if r["calibrated"]
                                  else "   [NOT CALIBRATED — diagnostic only]"))
        checks = [
            ("distinct genes            ", len(r["genes"]), HARD["ALL"]),
            ("up, partition             ", len(r["up_part"]), HARD["UP"]),
            ("down, partition           ", len(r["down_part"]), HARD["DOWN"]),
            ("genes with both directions", len(r["both"]), HARD["MIXED"]),
            ("significant in >=2 classes", len(r["multi"]), HARD["MULTI"]),
            ("up, union                 ", len(r["up_any"]), HARD["UP"]),
            ("down, union               ", len(r["down_any"]), HARD["DOWN"]),
        ]
        for label, got, want in checks:
            ok = got == want
            print(f"    {label}  {got:>6,}   vs hard-coded {want:>6,}   "
                  f"{'MATCH' if ok else ''}")
            if ok and r["calibrated"]:
                hits.append(f"{m}: {label.strip()}")

    print()
    hr("=")
    if len(hits) >= 3:
        print(f"  SOURCED. {len(hits)} of the hard-coded numbers reproduce:")
        for h in hits:
            print(f"    {h}")
        print("""
  Then the figure's numbers are right and the problem was only that nobody could
  prove it. Two things to do:
    1. Replace the hard-coded constants with a read of these tables, or at minimum
       cite the five filenames in a comment beside them, so the next person can check.
    2. **The legend must name the quantity.** "927 signature genes" beside a Results
       headline of "926 skipped exons" reads as a typo. Say what it is: genes carrying
       a significant event in any of the five classes under the composition-adjusted
       model. The two numbers are close by coincidence, and the legend has to say so.""")
    elif hits:
        print(f"  PARTIAL — only {len(hits)} match:")
        for h in hits:
            print(f"    {h}")
        print("""
  Some of the figure's constants reproduce and some do not. Do not publish the ones
  that do not. Work out what the non-matching ones counted — a different model, a
  pre-exclusion run, or an earlier filter — before the figure goes anywhere.""")
    else:
        print("""  ⚠ NOTHING REPRODUCES.

  The foreground was not built from any of the models tried here. Find the gene list
  that was actually submitted to g:Profiler:

      grep -n "query\\|gene_list\\|foreground\\|to_csv\\|gprofiler" 14_run_gprofiler.py
      ls -la *gene*list* *query*.txt *foreground* 2>/dev/null

  Until it is found, 927 / 636 / 291 cannot go on a published figure.""")
    hr("=")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default=".")
    ap.add_argument("--model", choices=list(MODELS) + ["all"], default="all",
                    help="which model's tables to read (default: compadj and base)")
    ap.add_argument("--list", action="store_true", dest="list_only")
    a = ap.parse_args()

    hr("=")
    print("93 v2 — sourcing the 927 on Figure 6")
    hr("=")
    print(f"  dir  {os.path.abspath(a.dir)}")

    if a.list_only:
        print("\n  tables this script will read, by model:\n")
        for m, pat in MODELS.items():
            print(f"  {m}")
            for c in CLASSES:
                f = pat.format(cls=c)
                print(f"    {'OK  ' if os.path.exists(os.path.join(a.dir, f)) else 'MISS'}"
                      f"  {f}")
        print("\n  --list: nothing read. Re-run without --list.")
        return

    want = list(MODELS) if a.model == "all" else [a.model]
    results = [one_model(a.dir, m) for m in want]
    verdict(results)


if __name__ == "__main__":
    main()
