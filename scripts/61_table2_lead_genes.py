#!/usr/bin/env python3
"""
61_table2_lead_genes.py -- Table 2: lead genes with per-cohort replication status.

THE SELECTION RULE, fixed here and printed in the table footnote. A lead-gene table is
a ranking, so the rule is the whole credibility of it:

  1. Start from MAGNet DCM-vs-NF events passing the project's locked definition,
     FDR < 0.05 AND |dPSI| >= 0.10, exon-level.
  2. One row per GENE: keep that gene's largest-|dPSI| significant event, and report how
     many significant events the gene carries, because a gene with eleven is a different
     claim from a gene with one.
  3. Rank by the number of INDEPENDENT replication cohorts agreeing in direction, then
     by |dPSI| in discovery.

Ranking on replication breadth first is deliberate. Ranking on discovery effect size
alone would order the table by the same quantity it was selected on, which is how
Figure 2C ended up able to illustrate the signature but not validate it. Replication
breadth is measured in cohorts that took no part in the selection, so it is the one
axis in this table that is not circular.

WHAT THIS TABLE MUST NOT DO, enforced in code rather than remembered:
  * NO cross-cohort effect sizes. The project's position is that PSI magnitude is not
    comparable across cohorts (different read lengths, depths, covariate models), so
    replication appears as DIRECTION ONLY -- "+", "-", "." -- and never as a number.
  * NO LBNL. 24 of its 33 samples are MAGNet individuals; it is a library-preparation
    replicate, not a replication cohort.
  * NO Colorado MXE. Reverse-replication slope 0.061; excluded project-wide.

SELF-CHECKS, which abort rather than warn:
  * the discovery SE significant count must be 926
  * each cohort's concordance across the whole significant set must land within a few
    points of the published per-cohort figures (79.6 / 81.2 / 69.2 / 76.1). If the event
    matching is wrong, these move first and by a lot.

  python 61_table2_lead_genes.py --list          # what it would read, nothing else
  python 61_table2_lead_genes.py                 # Table2.{tsv,md,xlsx}
  python 61_table2_lead_genes.py --top 25 --types SE,MXE
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd

# -- conventions copied from 45_cohort_slopes.py so the two cannot drift ------
NKEEP = {"SE": 3, "MXE": 5, "A3SS": 5, "A5SS": 5, "RI": 3}

# --------------------------------------------------------------------------
# Read length. Empty = the original 50 bp behaviour, reproducing the published
# Table2.{tsv,md,xlsx}; _rl125 = the re-quantified tables, written to
# Table2_rl125.*.
#
# Table 2 is DERIVED (line ~283 is .head(--top)), so a read-length change is a
# rebuild, not an edit: three genes leave the eligible pool and the same rule
# promotes whoever is next.
SUF = os.environ.get("TBL_SUFFIX", "")

# FROZEN_SE is the only self-check here that aborts. 926 at 50 bp, 813 at 125.
FROZEN_SE = int(os.environ.get("FROZEN_SE", "926"))

# KNOWN_PCT passes unchanged at 125 bp (measured 80.8 / 81.9 / 70.8 / 75.8
# against PCT_TOL 4.0), but is overridable so the gate checks what it claims:
#     KNOWN_PCT=Colorado=80.8,Brigham=81.9,Fuwai=70.8,Stanford=75.8
_kp = os.environ.get("KNOWN_PCT", "").strip()
# --------------------------------------------------------------------------
REPL = [("Colorado", "{ET}_Colorado_DCMvsNF"),
        ("Brigham", "Brigham_{ET}_HCM"),
        ("Fuwai", "Fuwai_{ET}_HCM"),
        ("Stanford", "Stanford_{ET}_HCM")]          # LBNL deliberately absent
DROP_COMBOS = {("Colorado", "MXE")}
DISC = "{ET}_v2filtered" + SUF                            # 54_fdr_only_counts --find-discovery
KNOWN_PCT = {"Colorado": 79.6, "Brigham": 81.2, "Fuwai": 69.2, "Stanford": 76.1}

if _kp:
    _d = dict(_p.split("=", 1) for _p in _kp.split(","))
    assert set(_d) == set(KNOWN_PCT), (
        f"KNOWN_PCT override must name exactly {sorted(KNOWN_PCT)}, "
        f"got {sorted(_d)}")
    KNOWN_PCT = {_k: float(_v) for _k, _v in _d.items()}
PCT_TOL = 4.0                                       # percentage points

MARK = {1: "+", -1: "−", 0: "."}               # same / opposite / not measurable


def load(base, et, need_gene=False):
    """Load an all_events table, preferring the post-exclusion _dedup version.

    Collapses to exon level exactly as the rest of the project does: truncate `key` to
    NKEEP[et] fields and keep the lowest-FDR row per collapsed key.
    """
    for p in (f"{base}_dedup.all_events.tsv", f"{base}.all_events.tsv"):
        if not (os.path.exists(p) and os.path.getsize(p)):
            continue
        t = pd.read_csv(p, sep="\t", low_memory=False)
        if "key" not in t.columns:
            return None, p, "no `key` column"
        d = next((c for c in ("dPSI_adj", "dPSI") if c in t.columns), None)
        if d is None:
            cand = [c for c in t.columns if c.lower().startswith("dpsi")]
            d = next((c for c in cand if "adj" in c.lower()),
                     cand[0] if cand else None)
        if d is None:
            return None, p, f"no dPSI column; have {list(t.columns)[:8]}"
        if "FDR" not in t.columns:
            return None, p, "no FDR column"
        if need_gene and "gene" not in t.columns:
            return None, p, "no `gene` column (needed for the discovery table)"
        t["ekey"] = [":".join(str(k).split(":")[:NKEEP[et]]) for k in t["key"]]
        t = t.loc[t.groupby("ekey")["FDR"].idxmin()].copy()
        t = t.rename(columns={d: "dpsi"})
        keep = ["ekey", "dpsi", "FDR"]
        for extra in ("gene", "chr", "exonStart", "se_adj", "mde80"):
            if extra in t.columns:
                keep.append(extra)
        return t[keep], p, None
    return None, base, "file not found"


def fmt(v, col):
    """One place that decides how a value is printed, so .tsv, .md and .xlsx agree."""
    if v is None or (isinstance(v, float) and not np.isfinite(v)):
        return ""
    if col in ("ΔPSI", "dPSI"):
        return f"{float(v):+.3f}"
    if col == "FDR":
        return f"{float(v):.3g}"
    return str(v)


def to_md(df):
    """Pipe table without pandas.to_markdown, which needs `tabulate`.

    `tabulate` is not installed in this project's environment -- the same trap as
    statsmodels. Nothing optional gets imported here.
    """
    cols = list(df.columns)
    cells = [[fmt(v, c) for v, c in zip(row, cols)]
             for row in df.itertuples(index=False, name=None)]
    w = [max(len(c), *(len(r[i]) for r in cells)) if cells else len(c)
         for i, c in enumerate(cols)]
    num = {"ΔPSI", "FDR", "#", "Sig. events", "Conc.", "Meas."}
    def line(vals):
        return "| " + " | ".join(
            v.rjust(w[i]) if cols[i] in num else v.ljust(w[i])
            for i, v in enumerate(vals)) + " |"
    out = [line(cols),
           "|" + "|".join(("-" * (w[i] + 2)) if cols[i] not in num
                          else ("-" * (w[i] + 1)) + ":" for i in range(len(cols))) + "|"]
    out += [line(r) for r in cells]
    return "\n".join(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--types", default="SE",
                    help="event types to include (default SE, which is what the "
                         "published per-cohort concordance figures are computed on)")
    ap.add_argument("--top", type=int, default=20)
    ap.add_argument("--fdr", type=float, default=0.05)
    ap.add_argument("--dpsi", type=float, default=0.10)
    ap.add_argument("--disc", default=DISC)
    ap.add_argument("--compadj", default="{ET}_v2compadj" + SUF,
                    help="composition-adjusted discovery table stem. Adds a column "
                         "saying whether each lead event is still significant with "
                         "cardiomyocyte and fibroblast scores in the model. Pass '' "
                         "to omit the column")
    ap.add_argument("--list", action="store_true",
                    help="print which files would be read and stop")
    ap.add_argument("--no-check", action="store_true",
                    help="skip the 926 and concordance self-checks. Use ONLY when "
                         "deliberately running a non-standard event set, and say so "
                         "wherever the output is used")
    ap.add_argument("--prefix", default="Table2" + SUF)
    a = ap.parse_args()
    types = [t.strip().upper() for t in a.types.split(",") if t.strip()]
    bad = [t for t in types if t not in NKEEP]
    if bad:
        sys.exit(f"[fatal] unknown event type(s) {bad}; known: {list(NKEEP)}")

    # ---- what will be read -------------------------------------------------
    print("=" * 78)
    print("INPUTS")
    print("=" * 78)
    plan = []
    for et in types:
        plan.append(("discovery", et, a.disc.replace("{ET}", et)))
        for coh, pat in REPL:
            if (coh, et) in DROP_COMBOS:
                print(f"  [skip] {coh} {et} -- excluded project-wide "
                      f"(reverse-replication slope 0.061)")
                continue
            plan.append((coh, et, pat.replace("{ET}", et)))
    for role, et, base in plan:
        hit = next((p for p in (f"{base}_dedup.all_events.tsv",
                                f"{base}.all_events.tsv")
                    if os.path.exists(p) and os.path.getsize(p)), None)
        print(f"  {role:<10} {et:<5} {hit or '*** MISSING: ' + base + ' ***'}")
    if a.list:
        print("\n[--list] stopping here.")
        return

    # ---- build -------------------------------------------------------------
    rows, checks = [], []
    for et in types:
        D, dpath, err = load(a.disc.replace("{ET}", et), et, need_gene=True)
        if D is None:
            sys.exit(f"[fatal] discovery {et}: {err} ({dpath})\n"
                     f"        Run with --list, or pass --disc with the right stem.")
        sig = D[(D.FDR < a.fdr) & (D.dpsi.abs() >= a.dpsi)].copy()
        print(f"\n[{et}] discovery {dpath}: {len(D):,} exon-level events, "
              f"{len(sig):,} significant")
        if et == "SE" and not a.no_check and len(sig) != FROZEN_SE:
            sys.exit(f"[fatal] SE significant count is {len(sig)}, frozen value is "
                     f"{FROZEN_SE}.\n        Either the wrong discovery table was read "
                     f"({dpath}) or the\n        convention has drifted. Resolve before "
                     f"building a table from it.")

        # composition-adjusted significance of the same events
        comp = None
        if a.compadj:
            C, cpath, cerr = load(a.compadj.replace("{ET}", et), et)
            if C is None:
                print(f"  [warn] composition-adjusted table: {cerr} ({cpath}) -- the "
                      f"comp column will read '?' throughout")
            else:
                cs = C[(C.FDR < a.fdr) & (C.dpsi.abs() >= a.dpsi)]
                comp = set(cs.ekey)
                print(f"  comp-adj   {len(cs):>5} significant in {cpath} "
                      f"({len(comp & set(sig.ekey))} of the {len(sig)} base events)")

        # per-cohort direction on the significant set
        dirs, conc = {}, {}
        for coh, pat in REPL:
            if (coh, et) in DROP_COMBOS:
                continue
            R, rpath, rerr = load(pat.replace("{ET}", et), et)
            if R is None:
                print(f"  [warn] {coh}: {rerr} ({rpath}) -- scored as not measurable")
                dirs[coh] = {}
                continue
            m = sig.merge(R[["ekey", "dpsi"]], on="ekey", how="inner",
                          suffixes=("", "_r"))
            same = int((np.sign(m.dpsi) == np.sign(m.dpsi_r)).sum())
            n = len(m)
            conc[coh] = (same, n)
            dirs[coh] = dict(zip(m.ekey,
                                 np.where(np.sign(m.dpsi) == np.sign(m.dpsi_r), 1, -1)))
            pct = 100 * same / n if n else float("nan")
            ref = KNOWN_PCT.get(coh)
            flag = ""
            if et == "SE" and ref is not None and n:
                off = abs(pct - ref)
                flag = "  <-- off by %.1f pt" % off if off > PCT_TOL else "  ok"
                checks.append((coh, pct, ref, off))
            print(f"  {coh:<10} {same:>5}/{n:<6} {pct:5.1f}%   (published "
                  f"{ref if ref is not None else '-'}){flag}")

        # one row per gene: largest |dpsi| significant event
        sig["absd"] = sig.dpsi.abs()
        sig = sig.sort_values("absd", ascending=False)
        per_gene = sig.groupby("gene", sort=False)
        for gene, g in per_gene:
            top = g.iloc[0]
            cells, nsame, nmeas = {}, 0, 0
            for coh, _pat in REPL:
                if (coh, et) in DROP_COMBOS:
                    cells[coh] = "n/a"
                    continue
                v = dirs.get(coh, {}).get(top.ekey, 0)
                cells[coh] = MARK[int(v)]
                if v != 0:
                    nmeas += 1
                    nsame += int(v == 1)
            rows.append(dict(
                gene=gene, event=et,
                comp_adj=("?" if comp is None
                          else ("+" if top.ekey in comp else ".")),
                locus=(f"{top.get('chr', '?')}:{int(top['exonStart'])}"
                       if "exonStart" in top and pd.notna(top.get("exonStart"))
                       else "?"),
                dPSI=round(float(top.dpsi), 3), FDR=float(top.FDR),
                n_sig_events=int(len(g)),
                **{f"rep_{c}": cells[c] for c, _ in REPL},
                n_concordant=nsame, n_measurable=nmeas))

    if not rows:
        sys.exit("[fatal] no rows built")

    # ---- self-check verdict ------------------------------------------------
    if checks and not a.no_check:
        worst = max(checks, key=lambda r: r[3])
        if worst[3] > PCT_TOL:
            print(f"\n[fatal] {worst[0]} concordance {worst[1]:.1f}% against a published "
                  f"{worst[2]}% -- off by {worst[3]:.1f} points.")
            print("        The per-cohort figures move first when event matching is")
            print("        wrong, so this is a matching failure, not a rounding one.")
            print("        Do not build Table 2 on it. Re-run with --list and check that")
            print("        each cohort's table is the post-exclusion (_dedup) one.")
            sys.exit(2)
        print("\n[ok] every cohort's concordance reproduces its published figure "
              f"within {PCT_TOL} points")

    T = pd.DataFrame(rows).sort_values(
        ["n_concordant", "n_measurable", "FDR"],
        ascending=[False, False, True]).head(a.top).reset_index(drop=True)
    T.insert(0, "rank", range(1, len(T) + 1))

    # ---- write -------------------------------------------------------------
    T.to_csv(f"{a.prefix}.tsv", sep="\t", index=False, float_format="%.4g")
    print(f"\n[write] {a.prefix}.tsv")

    hdr = {"rank": "#", "gene": "Gene", "event": "Event", "locus": "Locus",
           "dPSI": "ΔPSI", "FDR": "FDR", "n_sig_events": "Sig. events",
           "rep_Colorado": "Col", "rep_Brigham": "Bri", "rep_Fuwai": "Fuw",
           "rep_Stanford": "Sta", "n_concordant": "Conc.", "n_measurable": "Meas.",
           "comp_adj": "CompAdj"}
    M = T.rename(columns=hdr)
    LEGEND = (
        "\n\nEvents significant in the discovery cohort (166 dilated versus 166 "
        "non-failing hearts) at FDR < 0.05 and |\u0394PSI| \u2265 0.10, exon-level, one "
        "row per gene (the gene's largest-|\u0394PSI| event; *Sig. events* gives how "
        "many the gene carries). Ranked by the number of independent replication "
        "cohorts agreeing in direction, then by FDR. Replication is reported as "
        "**direction only** \u2014 `+` concordant, `\u2212` discordant, `.` not "
        "measurable in that cohort, `n/a` excluded \u2014 because \u0394PSI magnitude "
        "is not comparable across cohorts with different read lengths, depths and "
        "covariate models. Col, Colorado (dilated versus non-failing); Bri, Brigham; "
        "Fuw, Fuwai; Sta, Stanford (each hypertrophic versus donor). LBNL is not "
        "shown: 24 of its 33 samples are discovery-cohort individuals, so it is a "
        "library-preparation replicate rather than an independent cohort. Colorado "
        "mutually exclusive exons are excluded project-wide (reverse-replication "
        "slope 0.061). *CompAdj* `+` marks events still significant under the same "
        "rule with cardiomyocyte and fibroblast composition scores in the model. Where "
        "a lead gene is also a member of those marker panels, the event listed is a "
        "splicing change in that gene and not a change in its expression.\n")
    try:
        with open(f"{a.prefix}.md", "w") as fh:
            fh.write("**Table 2. Lead genes of the failing-myocardium splicing "
                     "signature, with replication status.**\n\n")
            fh.write(to_md(M))
            fh.write(LEGEND)
        print(f"[write] {a.prefix}.md")
    except Exception as e:                      # noqa: BLE001 - never lose the .tsv
        print(f"[warn] .md not written ({type(e).__name__}: {e}); "
              f"{a.prefix}.tsv is complete and is the file that matters")
    try:
        from openpyxl import Workbook
        from openpyxl.styles import Alignment, Font
        wb = Workbook()
        ws = wb.active
        ws.title = "Table 2"
        ws.append(list(M.columns))
        for c in ws[1]:
            c.font = Font(bold=True)
            c.alignment = Alignment(horizontal="center")
        for r in M.itertuples(index=False):
            ws.append(list(r))
        for i, col in enumerate(M.columns, 1):
            w = max([len(str(col))] + [len(str(v)) for v in M[col]]) + 2
            ws.column_dimensions[ws.cell(row=1, column=i).column_letter].width = min(w, 28)
        ws.freeze_panes = "A2"
        wb.save(f"{a.prefix}.xlsx")
        print(f"[write] {a.prefix}.xlsx")
    except ImportError:
        print("[note] openpyxl not available -- .xlsx skipped, .tsv and .md written")

    print("\n" + "=" * 78)
    print(to_md(M))
    print("=" * 78)
    print("\nBEFORE THIS GOES IN THE PAPER")
    print("  * A gene whose top event is discordant in most cohorts is still shown. That")
    print("    is intentional -- a lead-gene table that hides its failures is a ranking,")
    print("    not a result. Check none of the top rows is mostly '-'.")
    print("  * `.` means the event was not in that cohort's table at all, which is")
    print("    usually coverage, not disagreement. The counts separate the two.")
    nfull = int((T.n_concordant == T.n_measurable).sum())
    if nfull == len(T):
        print(f"  * ALL {len(T)} rows replicate in every cohort they are measurable in,")
        print("    so the rank column is really 'concordant everywhere, then lowest")
        print("    FDR'. That is a strong definition but a saturated one -- say it that")
        print("    way in the legend rather than implying the ranking discriminates")
        print("    across all 20 rows.")
    print("  * Check whether any lead gene is also in the CM/FB marker panels (DCN,")
    print("    CSRP3, RYR2, TNNI3, MYL2, MYL3, ACTN2, TCAP, MYBPC3, PLN, COL1A1,")
    print("    COL1A2, COL3A1, POSTN, FN1, LUM). It is not circular -- these are")
    print("    splicing events, not expression -- but a reviewer will read a fibroblast")
    print("    marker at rank 1 as composition, so the CompAdj column has to answer it.")
    print("  * The .docx for submission is built from Table2.tsv the same way Table 1's")
    print("    was; the footnote text above is the legend to carry into it.")


if __name__ == "__main__":
    main()
