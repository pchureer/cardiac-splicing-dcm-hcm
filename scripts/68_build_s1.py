#!/usr/bin/env python3
"""
68_build_s1.py -- Supplementary Table S1: composition and covariate robustness.

WHAT S1 IS FOR
--------------
The Results summarises it in two sentences and S5(d) cites the covariate ladder, so S1 is
where a reviewer checks the obvious objection: *the signature is cell composition, and you
adjusted it away or failed to adjust for it.* Two parts:

  (a) **Composition adjustment** -- how many significant events survive adding
      cardiomyocyte and fibroblast scores, and what happens to the ones that do not.
      The answer is threshold attrition, not confounding, and the evidence for that
      distinction is the effect retained by the LOST events plus a uniform standard-error
      inflation.
  (b) **The covariate ladder** -- the discovery count under five covariate models, all on
      the same exon-level event set. 926 / 899 / 835 / 866 / 888.

⚠ THE TRAP IN THIS TABLE: TWO CONVENTIONS
-----------------------------------------
S1's two halves are recorded in different counting conventions, and putting them in one
table without noticing would be the project's signature failure mode.

    composition retention, RESULT_composition_check.md
        SE  672 / 1074 = 62.6%        1074 is FULL-COORDINATE
        MXE 172 /  290 = 59.3%         290 is FULL-COORDINATE

    the Results, already written, quotes the EXON-LEVEL version
        "590 of the 926 events remain significant" = 63.7%
        and Table 2's binomial test uses that 63.7% as its background rate

    the covariate ladder is exon-level throughout
        926 / 899 / 835 / 866 / 888 on the same 62,971 exon-level events

**The manuscript is exon-level everywhere. S1 must be too.** So:

* **SE is already available** at exon level -- 590 / 926 -- and that is what S1 reports.
* **MXE is NOT.** The record holds only the full-coordinate 172 / 290. The exon-level MXE
  base is 231 and its retained count is nowhere in the project. This script computes it
  from the tables if they are present, and refuses to print the full-coordinate number
  beside exon-level ones if they are not.

This script computes nothing else. Everything else is read, checked and formatted.

USAGE
-----
    conda activate rnaseq_tools
    cd /data/project/chureeratp/RNA/HEART

    python 68_build_s1.py --list     # find the inputs, describe them, stop
    python 68_build_s1.py            # build once the columns are pinned
"""
import argparse
import glob
import os

import numpy as np
import pandas as pd

# ---------------------------------------------------------------- the record
# Exon-level -- what the manuscript uses and what S1 must report.
EXON = dict(se_base=926, se_kept=590, se_pct=63.7,
            mxe_base=231, mxe_kept=None, mxe_pct=None)      # MXE kept is NOT on record
# Full-coordinate -- what RESULT_composition_check.md holds. For cross-checking only.
FULL = dict(se_base=1074, se_kept=672, se_pct=62.6,
            mxe_base=290, mxe_kept=172, mxe_pct=59.3)
LADDER = [("base: sex + age + RIN", 926),
          ("+ cardiomyocyte / fibroblast composition", 899),
          ("+ composition + library size + TIN", 835),
          ("+ sequencing pool", 866),
          ("excluding donors < 25 y (62,482 tested)", 888)]
LADDER_TESTED = 62_971
# the distinction that makes the ~37% loss attrition rather than confounding
LOST = [("median effect retained, kept events", 1.067, 0.996),
        ("median effect retained, LOST events", 0.802, 0.735),
        ("fraction of lost events losing >half their effect", 0.082, 0.136),
        ("lost events still in the same direction", 0.998, 0.992),
        ("standard-error inflation, kept", 1.213, 1.198),
        ("standard-error inflation, lost", 1.215, 1.187)]
CONCORD = dict(se_rho=0.915, se_sign=99.9, mxe_rho=0.928, mxe_sign=99.7)
PROGRAM = dict(se="23/30 (77%) against 63% genome-wide",
               mxe="7/13 (54%) against 59% -- indistinguishable at n=13")
NKEEP = {"SE": 3, "MXE": 5, "A3SS": 5, "A5SS": 5, "RI": 3}
FDR_CUT, DPSI_CUT = 0.05, 0.10

# Pinned from the real files, 6 Oct. The composition-adjusted tables carry the SAME
# 20-column schema as the filtered ones, so the whole of S1 can be computed at exon
# level rather than quoted from a full-coordinate record.
BASE = "{}_v2filtered.all_events.tsv"
ADJ = "{}_v2compadj.all_events.tsv"
LADDER_FILE = "fdr_only_counts.tsv"          # carries BOTH conventions side by side
COLS = dict(key="key", fdr="FDR", dpsi="dPSI_adj", se="se_adj", gene="gene")
ALL_CLASSES = ["SE", "MXE", "A3SS", "A5SS", "RI"]
# full-coordinate retention from RESULT_composition_check.md and the lostwhy reports --
# CHECKS on the full-coordinate pass, never what S1 prints
FULL_RETENTION = {"SE": (672, 1074), "MXE": (172, 290),
                  "A3SS": (62, 90), "A5SS": (44, 86)}

PATTERNS = [
    ("covariate ladder", ["*ladder*", "*covariate*", "*find_discovery*",
                          "*fdr_only*", "54_*"]),
    ("composition-adjusted tables", ["*compadj*", "*cell*adj*", "*_cmfb*",
                                     "SE_*compadj*", "MXE_*compadj*"]),
    ("why events were lost", ["*why_*lost*", "10_*", "*lost*"]),
    ("discovery tables", ["SE_v2filtered*", "MXE_v2filtered*"]),
]


def describe(path, nmax=3):
    try:
        if os.path.getsize(path) == 0:
            print(f"      {path}  (empty)"); return
        if path.endswith((".tsv", ".csv", ".txt")):
            sep = "," if path.endswith(".csv") else "\t"
            t = pd.read_csv(path, sep=sep, low_memory=False, nrows=400)
            print(f"      {path}  ({len(t)} rows shown, {len(t.columns)} columns)")
            print(f"        columns: {', '.join(map(str, t.columns))}")
            for line in t.head(nmax).to_string(index=False).splitlines():
                print(f"          {line}")
        else:
            print(f"      {path}  ({os.path.getsize(path):,} bytes)")
    except Exception as e:                                          # noqa: BLE001
        print(f"      {path}  [unreadable: {e.__class__.__name__}]")


def discover():
    print("\nSEARCHING for the inputs S1 needs")
    found = {}
    for label, pats in PATTERNS:
        hits = sorted({p for pat in pats for p in glob.glob(pat) if os.path.isfile(p)})
        print(f"\n  {label}:")
        if not hits:
            print("      *** nothing matched ***")
        for p in hits[:6]:
            describe(p)
        if len(hits) > 6:
            print(f"      ... and {len(hits) - 6} more")
        found[label] = hits
    return found


def exon_sig(path, et):
    """Significant events at exon level, the paper's convention."""
    if not (os.path.exists(path) and os.path.getsize(path)):
        return None, f"{path} missing"
    t = pd.read_csv(path, sep="\t", low_memory=False)
    if "key" not in t.columns or "FDR" not in t.columns:
        return None, f"{path}: no key/FDR column"
    d = next((c for c in ("dPSI_adj", "dPSI") if c in t.columns), None)
    if d is None:
        cand = [c for c in t.columns if c.lower().startswith("dpsi")]
        d = next((c for c in cand if "adj" in c.lower()), cand[0] if cand else None)
    if d is None:
        return None, f"{path}: no dPSI column"
    t["ekey"] = [":".join(str(k).split(":")[:NKEEP[et]]) for k in t["key"]]
    t = t.loc[t.groupby("ekey")["FDR"].idxmin()]
    dd = pd.to_numeric(t[d], errors="coerce")
    sig = t[(pd.to_numeric(t["FDR"], errors="coerce") < FDR_CUT)
            & (dd.abs() >= DPSI_CUT)]
    return set(sig["ekey"]), None


def show_record(notes):
    print("\n" + "=" * 78)
    print("WHAT S1 MUST CONTAIN -- and the convention it must be in")
    print("=" * 78)

    print("\n(a) COMPOSITION ADJUSTMENT")
    se_exon = f"{EXON['se_kept']}/{EXON['se_base']} ({EXON['se_pct']}%)"
    se_full = f"{FULL['se_kept']}/{FULL['se_base']} ({FULL['se_pct']}%)"
    mxe_exon = f"?/{EXON['mxe_base']}  NOT ON RECORD"
    mxe_full = f"{FULL['mxe_kept']}/{FULL['mxe_base']} ({FULL['mxe_pct']}%)"
    print(f"    {'':<34}{'exon-level':>22}{'full-coordinate':>20}")
    print(f"    {'skipped exons, significant':<34}{se_exon:>22}{se_full:>20}")
    print(f"    {'mutually exclusive exons':<34}{mxe_exon:>22}{mxe_full:>20}")
    print("\n    ⚠ The manuscript is EXON-LEVEL everywhere, and the Results already")
    print(f"      quotes {EXON['se_kept']}/{EXON['se_base']} = {EXON['se_pct']}% --")
    print("      which is also the background rate in Table 2's binomial test.")
    print("      The MXE exon-level retention is NOT on record. Do not fill it with")
    print("      the full-coordinate 172/290; that would put two conventions in one row.")
    notes.append("MXE exon-level composition retention is not on record -- compute it "
                 "or omit the MXE row")

    print(f"\n    Concordance among base-significant events:")
    print(f"      SE   rho = {CONCORD['se_rho']}, {CONCORD['se_sign']}% same sign")
    print(f"      MXE  rho = {CONCORD['mxe_rho']}, {CONCORD['mxe_sign']}% same sign")

    print("\n    Why the loss is threshold attrition and not confounding --")
    print("    this distinction IS the supplement's argument:")
    print(f"    {'':<50}{'SE':>8}{'MXE':>8}")
    for lab, se, mxe in LOST:
        print(f"      {lab:<48}{se:>8.3f}{mxe:>8.3f}")
    print("      Standard-error inflation is UNIFORM across kept and lost (~1.2x), so")
    print("      the losses are marginal events crossing a threshold, not differential")
    print("      noise. Strong events get STRONGER (1.067) -- composition was diluting")
    print("      them, not confounding them.")
    print(f"\n    RBM20 program: SE {PROGRAM['se']}")
    print(f"                   MXE {PROGRAM['mxe']}")

    print("\n(b) THE COVARIATE LADDER -- all on the same "
          f"{LADDER_TESTED:,} exon-level events")
    print(f"    {'model':<44}{'significant SE':>16}")
    for lab, n in LADDER:
        print(f"      {lab:<42}{n:>16,}")
    lo, hi = min(n for _, n in LADDER), max(n for _, n in LADDER)
    print(f"      spread {hi - lo} events = {100 * (hi - lo) / hi:.1f}% of the headline;")
    print("      every model within 10%. That is what answers 'is this your covariates?'")


def load_raw(path, et):
    """Every row, with both the full key and the exon-level key. No collapse."""
    if not (os.path.exists(path) and os.path.getsize(path)):
        return None, f"{path} missing"
    t = pd.read_csv(path, sep="\t", low_memory=False)
    miss = [c for c in (COLS["key"], COLS["fdr"], COLS["dpsi"]) if c not in t.columns]
    if miss:
        return None, f"{path}: columns absent {miss}"
    t["ekey"] = [":".join(str(k).split(":")[:NKEEP[et]]) for k in t[COLS["key"]]]
    for c in (COLS["fdr"], COLS["dpsi"], COLS["se"]):
        if c in t.columns:
            t[c] = pd.to_numeric(t[c], errors="coerce")
    return t[np.isfinite(t[COLS["fdr"]]) & np.isfinite(t[COLS["dpsi"]])], None


def load(path, et, exon_level=True):
    """One row per event under the chosen convention, lowest FDR kept."""
    if not (os.path.exists(path) and os.path.getsize(path)):
        return None, f"{path} missing"
    t = pd.read_csv(path, sep="\t", low_memory=False)
    miss = [c for c in (COLS["key"], COLS["fdr"], COLS["dpsi"]) if c not in t.columns]
    if miss:
        return None, f"{path}: columns absent {miss}"
    t["ekey"] = ([":".join(str(k).split(":")[:NKEEP[et]]) for k in t[COLS["key"]]]
                 if exon_level else t[COLS["key"]].astype(str))
    t[COLS["fdr"]] = pd.to_numeric(t[COLS["fdr"]], errors="coerce")
    t[COLS["dpsi"]] = pd.to_numeric(t[COLS["dpsi"]], errors="coerce")
    if COLS["se"] in t.columns:
        t[COLS["se"]] = pd.to_numeric(t[COLS["se"]], errors="coerce")
    t = t[np.isfinite(t[COLS["fdr"]]) & np.isfinite(t[COLS["dpsi"]])]
    t = t.loc[t.groupby("ekey")[COLS["fdr"]].idxmin()]
    return t.set_index("ekey"), None


def spearman(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 3:
        return float("nan")
    rx, ry = pd.Series(x[ok]).rank(), pd.Series(y[ok]).rank()
    return float(np.corrcoef(rx, ry)[0, 1])


def composition(et, exon_level, notes, verbose=True):
    """Everything RESULT_composition_check.md reports, in one convention.

    Two different operations, and conflating them is a real error:

    * RETENTION is a set question -- is this exon still called significant under the
      adjusted model? Each table is collapsed on its own terms, because that is how each
      analysis would count. This is what reproduces the Results' 590/926.
    * The PAIRED statistics -- rho, same-sign, effect ratio, SE ratio -- compare one
      event measured two ways, so both values must come from the SAME underlying row.
      An earlier version collapsed each table independently and then compared the two
      representatives, which for a multi-transcript exon compares different rows. Every
      paired figure came out worse than the record by roughly the amount that mistake
      predicts. The base table now picks the representative and the adjusted values are
      looked up at that same full key.
    """
    br, e1 = load_raw(BASE.format(et), et)
    ar, e2 = load_raw(ADJ.format(et), et)
    if br is None or ar is None:
        if verbose:
            print(f"  {et}: {e1 or e2}")
        return None

    if exon_level:
        base_rep = br.loc[br.groupby("ekey")[COLS["fdr"]].idxmin()]
        adj_rep = ar.loc[ar.groupby("ekey")[COLS["fdr"]].idxmin()]
        unit = "ekey"
    else:
        base_rep, adj_rep, unit = br, ar, COLS["key"]
    base_rep = base_rep.set_index(unit)
    adj_rep = adj_rep.set_index(unit)

    sig = base_rep[(base_rep[COLS["fdr"]] < FDR_CUT)
                   & (base_rep[COLS["dpsi"]].abs() >= DPSI_CUT)]
    csig = adj_rep[(adj_rep[COLS["fdr"]] < FDR_CUT)
                   & (adj_rep[COLS["dpsi"]].abs() >= DPSI_CUT)]
    kept_ix = sig.index.intersection(csig.index)
    lost_ix = sig.index.difference(csig.index)

    # --- paired values, matched on the FULL key of the base representative
    a_by_key = ar.drop_duplicates(subset=[COLS["key"]]).set_index(COLS["key"])
    pair = sig.join(a_by_key[[COLS["dpsi"], COLS["se"], COLS["fdr"]]].add_suffix("_a"),
                    on=COLS["key"]) if exon_level else sig.join(
        a_by_key[[COLS["dpsi"], COLS["se"], COLS["fdr"]]].add_suffix("_a"))
    pair = pair[np.isfinite(pair[COLS["dpsi"] + "_a"])]
    n_unpaired = len(sig) - len(pair)
    if n_unpaired and verbose:
        print(f"  [--      ] {et}: {n_unpaired} significant events have no matching row"
              f" in the adjusted table")

    def stats(ix):
        sub = pair[pair.index.isin(ix)]
        if not len(sub):
            return dict(n=0, ratio=np.nan, se_ratio=np.nan, half=np.nan, same=np.nan)
        r = (sub[COLS["dpsi"] + "_a"].abs() / sub[COLS["dpsi"]].abs()).replace(
            [np.inf, -np.inf], np.nan)
        se = (sub[COLS["se"] + "_a"] / sub[COLS["se"]]).replace(
            [np.inf, -np.inf], np.nan)
        same = (np.sign(sub[COLS["dpsi"]]) == np.sign(sub[COLS["dpsi"] + "_a"])).mean()
        return dict(n=len(sub), ratio=float(r.median()),
                    se_ratio=float(se.median()), half=float((r < 0.5).mean()),
                    same=float(same))

    k, l = stats(kept_ix), stats(lost_ix)
    rho = spearman(pair[COLS["dpsi"]], pair[COLS["dpsi"] + "_a"])
    sign = float((np.sign(pair[COLS["dpsi"]])
                  == np.sign(pair[COLS["dpsi"] + "_a"])).mean())
    return dict(et=et, base=len(sig), kept=len(kept_ix), lost=len(lost_ix),
                pct=100 * len(kept_ix) / len(sig) if len(sig) else np.nan,
                rho=rho, sign=100 * sign, k=k, l=l, unpaired=n_unpaired)


def part_a(notes):
    print("\n" + "=" * 78)
    print("(a) COMPOSITION ADJUSTMENT -- computed at EXON level, the paper's convention")
    print("=" * 78)
    rows = []
    print(f"  {'class':<6}{'base':>7}{'kept':>7}{'lost':>7}{'retained':>10}"
          f"{'rho':>8}{'same sign':>11}")
    res = {}
    for et in ALL_CLASSES:
        r = composition(et, True, notes)
        if r is None:
            continue
        res[et] = r
        print(f"  {et:<6}{r['base']:>7}{r['kept']:>7}{r['lost']:>7}"
              f"{r['pct']:>9.1f}%{r['rho']:>8.3f}{r['sign']:>10.1f}%")
        rows.append(dict(part="a. composition", event_class=et, base_significant=r["base"],
                         retained=r["kept"], lost=r["lost"],
                         retained_pct=round(r["pct"], 1), dpsi_rho=round(r["rho"], 3),
                         same_direction_pct=round(r["sign"], 1)))

    # the gate: SE must reproduce what the Results already published
    if "SE" in res:
        got = res["SE"]
        ok_b = got["base"] == EXON["se_base"]
        ok_k = abs(got["kept"] - EXON["se_kept"]) <= 2
        print(f"\n  [{'ok' if ok_b and ok_k else 'MISMATCH':<8}] SE: {got['kept']}/"
              f"{got['base']} against the Results' {EXON['se_kept']}/{EXON['se_base']}")
        if not (ok_b and ok_k):
            notes.append(f"SE composition retention {got['kept']}/{got['base']} does not "
                         f"reproduce the Results' {EXON['se_kept']}/{EXON['se_base']}")
    if "MXE" in res:
        m = res["MXE"]
        print(f"  [NEW     ] MXE exon-level retention: {m['kept']}/{m['base']}"
              f" = {m['pct']:.1f}% -- NOT previously on record")
        print(f"              (the record holds only the full-coordinate"
              f" {FULL_RETENTION['MXE'][0]}/{FULL_RETENTION['MXE'][1]}"
              f" = {100 * FULL_RETENTION['MXE'][0] / FULL_RETENTION['MXE'][1]:.1f}%)")
        notes.append(f"MXE exon-level retention computed: {m['kept']}/{m['base']} "
                     f"({m['pct']:.1f}%) -- add it to the record")

    # the argument: attrition, not confounding
    print("\n  Why the loss is threshold attrition and not confounding")
    print(f"  {'class':<6}{'effect kept':>13}{'effect LOST':>13}{'SE infl kept':>14}"
          f"{'SE infl lost':>14}{'lost >half':>12}{'lost same dir':>15}")
    for et, r in res.items():
        print(f"  {et:<6}{r['k']['ratio']:>13.3f}{r['l']['ratio']:>13.3f}"
              f"{r['k']['se_ratio']:>14.3f}{r['l']['se_ratio']:>14.3f}"
              f"{100 * r['l']['half']:>11.1f}%{100 * r['l']['same']:>14.1f}%")
        rows.append(dict(part="a. attrition", event_class=et,
                         effect_retained_kept=round(r["k"]["ratio"], 3),
                         effect_retained_lost=round(r["l"]["ratio"], 3),
                         se_inflation_kept=round(r["k"]["se_ratio"], 3),
                         se_inflation_lost=round(r["l"]["se_ratio"], 3),
                         lost_losing_over_half_pct=round(100 * r["l"]["half"], 1),
                         lost_same_direction_pct=round(100 * r["l"]["same"], 1)))
    print("\n  Standard-error inflation is near-identical for kept and lost events, so")
    print("  the losses are marginal events crossing a threshold rather than events")
    print("  whose noise grew. Where the kept ratio exceeds 1, composition was DILUTING")
    print("  those events, not confounding them.")

    # and the full-coordinate pass, purely as a cross-check on the record
    print("\n  Cross-check -- the same computation at FULL coordinate, against the record")
    print("  (the record's figures ARE full-coordinate, so this is the exact test)")
    print(f"  {'class':<6}{'computed':>14}{'record':>14}{'rho':>9}{'same sign':>12}")
    for et in ALL_CLASSES:
        r = composition(et, False, notes, verbose=False)
        if r is None or et not in FULL_RETENTION:
            continue
        rk, rb = FULL_RETENTION[et]
        good = abs(r["kept"] - rk) <= 2 and abs(r["base"] - rb) <= 2
        got_s = f"{r['kept']}/{r['base']}"
        want_s = f"{rk}/{rb}"
        print(f"  {et:<6}{got_s:>14}{want_s:>14}{r['rho']:>9.3f}{r['sign']:>11.1f}%"
              f"   {'ok' if good else '<-- differs'}")
        if et == "SE":
            for lab, got, rec in (("rho", r["rho"], CONCORD["se_rho"]),
                                  ("same sign", r["sign"], CONCORD["se_sign"])):
                flag = "ok" if abs(got - rec) < 0.01 * max(1, rec) else "<-- differs"
                print(f"         SE {lab:<10} {got:.3f} against the record's"
                      f" {rec}   {flag}")
    return rows


def part_b(notes):
    print("\n" + "=" * 78)
    print("(b) THE COVARIATE LADDER")
    print("=" * 78)
    if os.path.exists(LADDER_FILE):
        t = pd.read_csv(LADDER_FILE, sep="\t")
        col = "exon-level_+|dPSI|>=0.10"
        m = t[(t["contrast"].astype(str).str.contains("DCM vs NF"))
              & (t["event_type"].astype(str) == "SE")]
        if len(m) and col in t.columns:
            got = int(pd.to_numeric(m.iloc[0][col], errors="coerce"))
            good = got == LADDER[0][1]
            print(f"  [{'ok' if good else 'MISMATCH':<8}] {LADDER_FILE} gives the base"
                  f" model {got:,} exon-level (ladder says {LADDER[0][1]:,})")
            if not good:
                notes.append(f"base model {got} vs ladder {LADDER[0][1]}")
            fc = int(pd.to_numeric(m.iloc[0]["full-coordinate_+|dPSI|>=0.10"],
                                   errors="coerce"))
            print(f"  [--      ] the same row full-coordinate: {fc:,}"
                  f" -- the two conventions, side by side in one file")
    rows = []
    print(f"\n  {'model':<44}{'significant SE':>16}")
    for lab, n in LADDER:
        print(f"    {lab:<42}{n:>16,}")
        rows.append(dict(part="b. covariate ladder", model=lab, significant_SE=n,
                         tested=LADDER_TESTED))
    lo, hi = min(n for _, n in LADDER), max(n for _, n in LADDER)
    print(f"\n  Spread {hi - lo} events, {100 * (hi - lo) / hi:.1f}% of the headline."
          f" Every model within 10%.")
    print("  The ladder's counts are from 54_fdr_only_counts.py --find-discovery and are")
    print("  not in a TSV; they are checked above only at the base model.")
    notes.append("ladder rows beyond the base model are not in a file -- "
                 "re-run 54 --find-discovery to deposit them")
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--base", default="SE_v2filtered.all_events.tsv")
    ap.add_argument("--compadj", help="composition-adjusted SE table, if named oddly")
    ap.add_argument("--base-mxe", default="MXE_v2filtered.all_events.tsv")
    ap.add_argument("--compadj-mxe", help="composition-adjusted MXE table")
    ap.add_argument("--out", default="S1")
    a = ap.parse_args()

    print("=" * 78)
    print("S1  composition and covariate robustness")
    print("=" * 78)
    print("  Answers the obvious objection: is the signature cell composition?")
    print("  The Results summarises it; S5(d) cites the covariate ladder.")

    found = discover()
    notes = []
    if a.list:
        show_record(notes)
        print("\n[--list] stopping here. Send this back and I will pin the columns.")
        print("  What I most need: which file holds the COMPOSITION-ADJUSTED event")
        print("  tables, so the exon-level MXE retention can be computed rather than")
        print("  taken from the full-coordinate record.")
        return

    rows = part_a(notes) + part_b(notes)
    if rows:
        T = pd.DataFrame(rows)
        T.to_csv(f"{a.out}.tsv", sep="\t", index=False)
        print(f"\n[write] {a.out}.tsv  ({len(T)} rows)")
    print("\n" + "=" * 78)
    for n in notes:
        print(f"  >>> {n}")


if __name__ == "__main__":
    main()
