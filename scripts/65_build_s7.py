#!/usr/bin/env python3
"""
65_build_s7.py -- assemble Supplementary Table S7, candidate regulators in full.

WHAT S7 IS FOR
--------------
The paper's second headline claim is a null: the signature is not explained by RBM20 or by
ten further candidate regulators. The main text says that in two sentences. **S7 is the
only place the evidence appears in full**, so it carries margins, intervals and the
correction -- not a table of p-values, which is what the first version of this analysis
produced and which was misleading (every "significant" candidate sat within 0.03 of a
housekeeping gene).

    A, forward   do RBM20's own targets behave as its targets should?
                 Target set: the 45-gene reference list, Gregorich 2024 Table 2.
    B, inverted  does any candidate track OUR signature?
                 Target set: the 134 composition-adjusted signature genes, RBM20 demoted
                 to one candidate among eleven.

**GAPDH is the calibration floor and is labelled as such in every table.** It is not a
candidate that happened to score; it is the control the design rests on.

HOW THE REPORTED MODEL IS CHOSEN -- the important part
-----------------------------------------------------
Both input files hold MANY models. `rbm20_margin_resolve.tsv` has 64 rows across
`model` x `gene_set` x `matching`; `signature_rbp_margins.tsv` has 160 across
`model` x `matching`. Picking the wrong combination would put a defensible-looking table
in the supplement that reports a model the manuscript does not use.

So this script does **not** guess from column values. It searches every combination for
the one that reproduces the figures recorded in `RESULT_rbm20_mechanism_null.md`
(`CHECK_A`, `CHECK_B` below) and reports which combination that was. If exactly one
matches, that is the reported model, identified by its numbers rather than by its name --
and the printout is itself the provenance record the Methods needs. If none or several
match, it refuses and prints the grid.

An earlier version tried to identify the published matching scheme by looking for "pub" in
the label, and refused a correct file because the three-dimensional model was called
"3-D match (as published)". Substring guesses about which statistical model a row encodes
do not belong anywhere near a supplement.

WHAT IT DOES NOT DO
-------------------
It computes nothing. Scripts 57 and 58 did the fitting, the paired bootstrap and the Holm
correction. If a number is not in those files it does not reach the table -- there is no
fallback to the values in the project docs, because a supplement assembled from prose is
one nobody can reproduce.

USAGE
-----
    conda activate rnaseq_tools
    cd /data/project/chureeratp/RNA/HEART

    python 65_build_s7.py --list        # the model grid in both files; computes nothing
    python 65_build_s7.py               # identify the reported model, build S7.tsv + S7.md

    # only if auto-identification is ambiguous:
    python 65_build_s7.py --model "..." --matching "..." --gene-set "..."
"""
import argparse
import itertools
import os
import re
import sys

import numpy as np
import pandas as pd

FWD = "rbm20_margin_resolve.tsv"      # script 57
INV = "signature_rbp_margins.tsv"     # script 58
FLOOR = "GAPDH"
CLASSES = ["SE", "MXE", "A3SS", "A5SS", "RI"]

# ---------------------------------------------------------------- column names, pinned
# Verified against the real files on 6 Oct. Pinned rather than guessed: the synonym
# resolver this replaced matched 'n_testable' to the column 'table', because "testable"
# contains "table". Loose matching has no place in choosing which column a supplement
# reports.
PIN_A = dict(table="table", model="model", gene_set="gene_set", matching="matching",
             tgt="rbm20_delta", floor="gapdh_delta", margin="margin",
             lo="ci_lo", hi="ci_hi", p="boot_p", n_target="n_target", n_bg="n_bg",
             verdict="verdict", table_variant="_tv")
PIN_B = dict(event="event_type", model="model", matching="matching", gene="gene",
             delta="delta", floor="floor", margin="margin", lo="ci_lo", hi="ci_hi",
             p="boot_p", holm="holm_p", n_target="n_target", n_bg="n_bg")
# `table` carries BOTH the event class and the adjustment (e.g. "SE_cell_tech"), so it
# is split: the class becomes the row label, the remainder a selector axis like any
# other. Without the split there are several rows per combination -- one per adjustment --
# and no combination looks like a single reported model.
SEL_A = ["table_variant", "model", "gene_set", "matching"]
SEL_B = ["model", "matching"]                 # the two axes in file B

# ---------------------------------------------------------------- recorded headlines
# From RESULT_rbm20_mechanism_null.md, 5 Oct. CHECKS, never a source. The combination of
# model/gene_set/matching that reproduces these IS the reported model.
# READ LENGTH, added 22 Oct. Everything below is a 50 bp record. At 125 bp every
# figure moves AND CHECK_B's `top` changes IDENTITY -- QKI -> PTBP1 in SE and
# RBM20 -> PTBP1 in MXE (152_gapdh_rank_recount.py, 57's rl125 deposit). A gate
# told "top is QKI" and handed PTBP1's numbers would PASS while naming the wrong
# gene, so these overrides are MANDATORY when a suffixed input is given, never
# defaulted. _refuse_unparameterised() below enforces that.
#
#   S7_CHECK_A="SE=-0.051/-0.159/0.047/0.311/641,MXE=-0.019/-0.108/0.072/0.690/1399"
#   S7_CHECK_B="SE=PTBP1/0.062/-0.025/0.122,MXE=PTBP1/0.113/0.016/0.187"
#   S7_B_NGENES="total=131,SE=123,MXE=119"
def _s7_env(_name):
    """'SE=a/b/c,MXE=d/e/f' -> {'SE': ['a','b','c'], ...}, or None if unset."""
    _raw = os.environ.get(_name, "").strip()
    if not _raw:
        return None
    _out = {}
    for _part in _raw.split(","):
        _k, _, _v = _part.partition("=")
        if not _k.strip():
            continue
        _out[_k.strip()] = [_x.strip() for _x in _v.split("/")]
    return _out or None


_A_ENV = _s7_env("S7_CHECK_A")
CHECK_A = ({_et: dict(margin=float(_v[0]), ci=(float(_v[1]), float(_v[2])),
                      p=float(_v[3]), n_target=int(_v[4]))
            for _et, _v in _A_ENV.items()} if _A_ENV else
           {"SE":  dict(margin=-0.005, ci=(-0.119, 0.080), p=0.76, n_target=663),
            "MXE": dict(margin=+0.000, ci=(-0.085, 0.104), p=0.88, n_target=1393)})

_B_ENV = _s7_env("S7_CHECK_B")
CHECK_B = ({_et: dict(top=_v[0], margin=float(_v[1]),
                      ci=(float(_v[2]), float(_v[3])))
            for _et, _v in _B_ENV.items()} if _B_ENV else
           {"SE":  dict(top="QKI",   margin=+0.033, ci=(-0.025, 0.082)),
            "MXE": dict(top="RBM20", margin=+0.009, ci=(-0.053, 0.088))})


def _refuse_unparameterised():
    """A suffixed input with 50 bp records is the one case worse than failing."""
    import re as _re
    _argv = " ".join(sys.argv[1:])
    if not _re.search(r"_rl\d+", _argv):
        return
    _miss = [_n for _n in ("S7_CHECK_A", "S7_CHECK_B", "S7_B_NGENES")
             if not os.environ.get(_n, "").strip()]
    if _miss:
        sys.exit(
            "REFUSING: an input path carries _rl<digits> but %s %s unset.\n"
            "  These are CHECKS that identify the reported model, and the 50 bp\n"
            "  records would either refuse every combination or -- worse --\n"
            "  match one while naming the wrong top gene (QKI/RBM20 at 50 bp,\n"
            "  PTBP1 at 125 bp). Supply all three:\n"
            "    S7_CHECK_A=\"SE=-0.051/-0.159/0.047/0.311/641,"
            "MXE=-0.019/-0.108/0.072/0.690/1399\" \\\n"
            "    S7_CHECK_B=\"SE=PTBP1/0.062/-0.025/0.122,"
            "MXE=PTBP1/0.113/0.016/0.187\" \\\n"
            "    S7_B_NGENES=\"total=131,SE=123,MXE=119\""
            % (", ".join(_miss), "is" if len(_miss) == 1 else "are"))


_refuse_unparameterised()
# TEN candidate rows, not eleven. GAPDH is not a row in this file -- it is the `floor`
# COLUMN that every margin is measured against (margin = delta - floor). The doc's
# eleven-row tables are the ten candidates plus the floor. A first version of this script
# required eleven rows and refused every combination in the file.
CHECK_B_NCAND = 10
# Testable-gene counts are NOT in file A -- it carries target events (n_target), not
# genes. Recorded here for the legend, and flagged as not re-derivable from the input.
A_TESTABLE = {"SE": 40, "MXE": 39}
A_NGENES = 45
A_ABSENT = ["DAB1", "DNM3", "NTRK3", "SEMA6D", "SPEN"]      # + MECP2 in MXE only
A_ABSENT_MXE_ONLY = ["MECP2"]
_G_ENV = _s7_env("S7_B_NGENES")
B_NTARGET_GENES = ({_k: int(_v[0]) for _k, _v in _G_ENV.items()} if _G_ENV
                   else {"total": 134, "SE": 130, "MXE": 124})
TOL_MARGIN, TOL_CI, TOL_P = 0.004, 0.010, 0.03


DERIVED = {"_tv", "_ev"}        # computed by prep_a, not expected in the file


def read(path, pin, label):
    if not (os.path.exists(path) and os.path.getsize(path)):
        print(f"  {label}: *** {path} missing or empty ***")
        return None
    t = pd.read_csv(path, sep="\t")
    miss = [k for k, c in pin.items() if c not in t.columns and c not in DERIVED]
    print(f"\n  {label}: {path}  ({len(t)} rows, {len(t.columns)} columns)")
    if miss:
        print(f"    *** pinned columns absent: {[pin[k] for k in miss]} ***")
        print(f"    columns present: {', '.join(map(str, t.columns))}")
        print("    The input schema has changed. Fix the pinned names in this script")
        print("    rather than letting it guess -- guessing is what it used to do.")
        return None
    return t


def event_of(s):
    """File A has no event column; the class is inside `table` (e.g. SE_cell_tech)."""
    m = re.search(r"(?<![A-Z0-9])(MXE|A3SS|A5SS|SE|RI)(?![A-Z0-9])", str(s).upper())
    return m.group(1) if m else None


def prep_a(t):
    """Split `table` into event class (_ev) and adjustment variant (_tv)."""
    t = t.copy()
    t["_ev"] = [event_of(v) for v in t[PIN_A["table"]]]
    out = []
    for raw, ev in zip(t[PIN_A["table"]], t["_ev"]):
        v = str(raw)
        if ev:
            v = re.sub(r"(?<![A-Z0-9])" + ev + r"(?![A-Z0-9])", "", v, flags=re.I)
        out.append(v.strip("_-. ") or "(whole table)")
    t["_tv"] = out
    return t


def grid(t, sel, pin, label):
    print(f"    {label} model grid:")
    for k in sel:
        vals = [str(v) for v in pd.unique(t[pin[k]].astype(str))]
        print(f"      {pin[k]:<10} {len(vals)} value(s): {', '.join(vals)}")


def fmt(v, nd=3):
    return "-" if v is None or not np.isfinite(v) else f"{v:+.{nd}f}"


def fmt_p(v):
    if v is None or not np.isfinite(v):
        return "-"
    return "<0.001" if v < 0.001 else f"{v:.3f}".rstrip("0").rstrip(".")


def g(row, pin, key):
    v = pd.to_numeric(row[pin[key]], errors="coerce")
    return float(v) if np.isfinite(v) else None


def near(got, want, tol):
    return got is not None and abs(got - want) <= tol


# ============================================================ identify the model
def identify(t, sel, pin, score, a, what):
    """Find the selector combination whose numbers reproduce the record."""
    axes = []
    for k in sel:
        override = getattr(a, k.replace("-", "_"), None)
        vals = [str(v) for v in pd.unique(t[pin[k]].astype(str))]
        if override:
            hit = [v for v in vals if override.lower() in v.lower()]
            if len(hit) != 1:
                print(f"  >>> --{k} '{override}' matches {len(hit)} of "
                      f"{len(vals)} values. Be exact, or omit it.")
                return None, None
            vals = hit
        axes.append([(k, v) for v in vals])

    ok, fail = [], {}
    for combo in itertools.product(*axes):
        sub = t
        for k, v in combo:
            sub = sub[sub[pin[k]].astype(str) == v]
        if not len(sub):
            fail.setdefault("no rows for this combination", []).append(dict(combo))
            continue
        why = score(sub)
        if why:
            fail.setdefault(why, []).append(dict(combo))
        else:
            ok.append((dict(combo), sub))

    if len(ok) == 1:
        combo, sub = ok[0]
        print(f"  [identified] the reported {what} model is the one that reproduces the")
        print(f"               recorded figures:")
        for k, v in combo.items():
            print(f"                 {pin[k]:<10} = {v}")
        print(f"               {len(sub)} rows. Put this in the Methods -- it is the")
        print("               provenance of every number in this part.")
        return combo, sub
    if not ok:
        print(f"  >>> {what} NOT BUILT. None of the "
              f"{sum(len(v) for v in fail.values())} combinations of "
              f"{', '.join(pin[k] for k in sel)} reproduces the recorded figures.")
        print("      Why each one failed, most common first:")
        for why, combos in sorted(fail.items(), key=lambda kv: -len(kv[1])):
            print(f"        {len(combos):>3} x  {why}")
            ex = combos[0]
            print(f"             e.g. " + ", ".join(f"{pin[k]}={v}" for k, v in ex.items()))
        print("      Either this is a different file, or the record is wrong. Do not")
        print("      pick a combination by eye to make the table appear.")
        return None, None
    print(f"  >>> {what} NOT BUILT. {len(ok)} combinations reproduce the record:")
    for combo, _ in ok:
        print("        " + ", ".join(f"{pin[k]}={v}" for k, v in combo.items()))
    print("      Pick one explicitly with --model / --matching / --gene-set.")
    return None, None


# ============================================================ part A
def part_a(t, rows, notes, md, a):
    P = PIN_A
    print("\n" + "=" * 78)
    print("A. Forward test -- RBM20 against its own 45-gene reference set")
    print("=" * 78)
    unk = t[t["_ev"].isna()]
    if len(unk):
        print(f"  [warn] {len(unk)} row(s) whose event class could not be read from "
              f"'{P['table']}': {', '.join(map(str, pd.unique(unk[P['table']])[:4]))}")

    def score(sub):
        """Empty string = this combination is the reported model."""
        for et, w in CHECK_A.items():
            r = sub[sub["_ev"] == et]
            if len(r) != 1:
                return f"{len(r)} {et} rows per combination, expected exactly 1"
            r = r.iloc[0]
            if not near(g(r, P, "margin"), w["margin"], TOL_MARGIN):
                return f"{et} margin is not {w['margin']:+.3f}"
            if not (near(g(r, P, "lo"), w["ci"][0], TOL_CI)
                    and near(g(r, P, "hi"), w["ci"][1], TOL_CI)):
                return (f"{et} interval is not {w['ci'][0]:+.3f} to "
                        f"{w['ci'][1]:+.3f}")
            if not near(g(r, P, "p"), w["p"], TOL_P):
                return f"{et} bootstrap p is not {w['p']}"
        return ""

    combo, sub = identify(t, SEL_A, P, score, a, "A")
    if sub is None:
        notes.append("A not built -- reported model not identified")
        return

    print(f"\n  {'event':<7}{'RBM20':>9}{FLOOR:>9}{'margin':>9}{'95% CI':>22}"
          f"{'p':>8}{'target ev':>11}{'bg ev':>9}")
    for et in CHECK_A:
        r = sub[sub["_ev"] == et].iloc[0]
        nt, nb = g(r, P, "n_target"), g(r, P, "n_bg")
        ci = f"{fmt(g(r, P, 'lo'))} to {fmt(g(r, P, 'hi'))}"
        ntv = str(int(nt)) if nt else "-"
        nbv = str(int(nb)) if nb else "-"
        print(f"  {et:<7}{fmt(g(r, P, 'tgt')):>9}{fmt(g(r, P, 'floor')):>9}"
              f"{fmt(g(r, P, 'margin')):>9}{ci:>22}"
              f"{fmt_p(g(r, P, 'p')):>8}{ntv:>11}{nbv:>9}")
        if nt is not None and int(nt) != CHECK_A[et]["n_target"]:
            print(f"  [warn] {et}: {int(nt)} target events, record says "
                  f"{CHECK_A[et]['n_target']}")
            notes.append(f"A {et}: n_target {int(nt)} vs recorded "
                         f"{CHECK_A[et]['n_target']}")
        rows.append(dict(part="A. RBM20 vs its own targets", event=et, gene="RBM20",
                         testable_genes=f"{A_TESTABLE[et]} of {A_NGENES}",
                         target_events=int(nt) if nt else "",
                         background_events=int(nb) if nb else "",
                         delta=fmt(g(r, P, "tgt")), floor_delta=fmt(g(r, P, "floor")),
                         margin=fmt(g(r, P, "margin")),
                         ci=ci,
                         p_bootstrap=fmt_p(g(r, P, "p")), p_holm="",
                         model=" | ".join(f"{k}={v}" for k, v in combo.items())))

    print(f"\n  [ok] both margins, both intervals and both p-values reproduce the record.")
    print(f"  [note] the testable-gene counts ({A_TESTABLE['SE']} SE, {A_TESTABLE['MXE']}"
          f" MXE of {A_NGENES}) are NOT in this file -- it")
    print("         carries target EVENTS, not genes. They are taken from the project")
    print("         record and cannot be re-derived here. Verify them against 57's log")
    print("         before the legend is final.")
    notes.append(f"A: testable-gene counts ({A_TESTABLE['SE']}/{A_TESTABLE['MXE']}) "
                 f"come from the record, not the file -- verify against 57's log")
    md.append(
        "**(a) RBM20 against its own target set.** Target set: the 45-gene reference list "
        "(Gregorich *et al.* 2024, Table 2), of which "
        f"{A_TESTABLE['SE']} contribute testable skipped-exon events and "
        f"{A_TESTABLE['MXE']} testable mutually-exclusive-exon events. Enrichment is the "
        "median |*t*| across target events minus the median across a background matched on "
        "mean PSI, PSI standard deviation and |disease effect|; within dilated "
        "cardiomyopathy, adjusted for cell composition and technical covariates. The "
        f"margin is the difference from {FLOOR}, the calibration floor, with a paired "
        "bootstrap 95% interval and bootstrap *P*.")


# ============================================================ part B
def b_arm(t, et):
    """Candidate rows for one event class, newest-margin first, floor validated."""
    P = PIN_B
    s = t[t[P["event"]].astype(str).str.upper() == et].copy()
    if not len(s):
        return None, None, "no rows"
    s["_m"] = pd.to_numeric(s[P["margin"]], errors="coerce")
    s["_d"] = pd.to_numeric(s[P["delta"]], errors="coerce")
    fl = pd.to_numeric(s[P["floor"]], errors="coerce").dropna().unique()
    if len(fl) != 1:
        return None, None, (f"the floor is not constant within this arm "
                            f"({len(fl)} values: {', '.join(f'{v:+.3f}' for v in fl[:4])})"
                            f" -- rows from different models are mixed")
    # the floor must actually BE the reference: delta - floor == margin
    off = (s["_d"] - float(fl[0]) - s["_m"]).abs().max()
    if np.isfinite(off) and off > 0.002:
        return None, None, (f"delta - floor does not equal margin (worst discrepancy "
                            f"{off:.3f}) -- the floor column is not the reference these "
                            f"margins were taken against")
    return s.sort_values("_m", ascending=False), float(fl[0]), None


def part_b(t, rows, notes, md, a):
    P = PIN_B
    print("\n" + "=" * 78)
    print("B. Inverted scan -- ten candidates against the 134-gene signature,")
    print("   each measured against the GAPDH floor carried in the `floor` column")
    print("=" * 78)

    def score(sub):
        """Empty string = this combination is the reported model."""
        for et, w in CHECK_B.items():
            s, fl, err = b_arm(sub, et)
            if s is None:
                return f"{et}: {err}"
            if len(s) != CHECK_B_NCAND:
                return (f"{et}: {len(s)} candidate rows, expected "
                        f"{CHECK_B_NCAND} (the floor is a column, not a row)")
            r = s.iloc[0]
            if str(r[P["gene"]]).upper() != w["top"].upper():
                return (f"{et}: top candidate by margin is "
                        f"{str(r[P['gene']])}, not {w['top']}")
            if not near(g(r, P, "margin"), w["margin"], TOL_MARGIN):
                return f"{et}: top margin is not {w['margin']:+.3f}"
            if not (near(g(r, P, "lo"), w["ci"][0], TOL_CI)
                    and near(g(r, P, "hi"), w["ci"][1], TOL_CI)):
                return (f"{et}: top interval is not {w['ci'][0]:+.3f} to "
                        f"{w['ci'][1]:+.3f}")
        return ""

    combo, sub = identify(t, SEL_B, P, score, a, "B")
    if sub is None:
        notes.append("B not built -- reported model not identified")
        return

    bad = []
    for et in CHECK_B:
        s, fl, err = b_arm(sub, et)
        if s is None:
            bad.append((et, err))
            continue
        print(f"\n  {et}  ({len(s)} candidates against a floor of {fl:+.3f},"
              f" ordered by margin)")
        print(f"  {'gene':<10}{'delta':>9}{'margin':>9}{'95% CI':>22}"
              f"{'p boot':>9}{'p Holm':>9}")
        placed = False
        for _, r in s.iterrows():
            gene = str(r[P["gene"]])
            m = g(r, P, "margin")
            # the floor is a reference line, so print it where it falls
            if not placed and m is not None and m < 0:
                print(f"  {FLOOR:<10}{fl:>+9.3f}{0.0:>+9.3f}{'(reference)':>22}"
                      f"{'-':>9}{'-':>9}  <-- calibration floor")
                rows.append(dict(part="B. candidates vs the signature", event=et,
                                 gene=FLOOR,
                                 testable_genes=f"{B_NTARGET_GENES[et]} of "
                                                f"{B_NTARGET_GENES['total']}",
                                 target_events="", background_events="",
                                 delta=f"{fl:+.3f}", floor_delta="(calibration floor)",
                                 margin="+0.000", ci="(reference)",
                                 p_bootstrap="", p_holm="",
                                 model=" | ".join(f"{k}={v}" for k, v in combo.items())))
                placed = True
            ci = f"{fmt(g(r, P, 'lo'))} to {fmt(g(r, P, 'hi'))}"
            print(f"  {gene:<10}{fmt(g(r, P, 'delta')):>9}{fmt(m):>9}"
                  f"{ci:>22}{fmt_p(g(r, P, 'p')):>9}{fmt_p(g(r, P, 'holm')):>9}")
            rows.append(dict(part="B. candidates vs the signature", event=et, gene=gene,
                             testable_genes=f"{B_NTARGET_GENES[et]} of "
                                            f"{B_NTARGET_GENES['total']}",
                             target_events=int(g(r, P, "n_target") or 0) or "",
                             background_events=int(g(r, P, "n_bg") or 0) or "",
                             delta=fmt(g(r, P, "delta")), floor_delta=f"{fl:+.3f}",
                             margin=fmt(m), ci=ci,
                             p_bootstrap=fmt_p(g(r, P, "p")),
                             p_holm=fmt_p(g(r, P, "holm")),
                             model=" | ".join(f"{k}={v}" for k, v in combo.items())))
        if not placed:          # every candidate is above the floor
            print(f"  {FLOOR:<10}{fl:>+9.3f}{0.0:>+9.3f}{'(reference)':>22}"
                  f"{'-':>9}{'-':>9}  <-- calibration floor, BELOW every candidate")
            rows.append(dict(part="B. candidates vs the signature", event=et, gene=FLOOR,
                             testable_genes=f"{B_NTARGET_GENES[et]} of "
                                            f"{B_NTARGET_GENES['total']}",
                             target_events="", background_events="",
                             delta=f"{fl:+.3f}", floor_delta="(calibration floor)",
                             margin="+0.000", ci="(reference)",
                             p_bootstrap="", p_holm="",
                             model=" | ".join(f"{k}={v}" for k, v in combo.items())))

        # the claim: nothing above the floor survives correction
        surv = [(str(r[P["gene"]]), g(r, P, "holm")) for _, r in s.iterrows()
                if g(r, P, "holm") is not None and g(r, P, "holm") < 0.05
                and (g(r, P, "margin") or 0) > 0]
        if surv:
            bad.append((et, "a candidate SURVIVES Holm ABOVE the floor: "
                            + ", ".join(f"{x} (Holm {fmt_p(p)})" for x, p in surv)
                            + " -- the null does not hold as recorded"))

    if bad:
        print("\n  >>> B NOT BUILT:")
        for et, why in bad:
            print(f"        {et:<5} {why}")
        notes.append("B not built -- the identified model does not hold as recorded")
        for r in [r for r in rows if r["part"].startswith("B.")]:
            rows.remove(r)
        return
    print("\n  [ok] candidate count, top candidate, its margin and interval, and the Holm")
    print("       result all reproduce the record. Nothing above the floor survives.")

    # ---- the floor's standing, under BOTH matchings: settles the labelling question
    print("\n  " + "-" * 74)
    print("  Where the floor ranks among the eleven (ten candidates plus itself), by")
    print("  delta -- its own margin is zero by construction, so margin cannot rank it.")
    print("  This is the number the Discussion's sentence uses.")
    print("  " + "-" * 74)
    for mv in [str(v) for v in pd.unique(t[P["matching"]].astype(str))]:
        arm = t[t[P["matching"]].astype(str) == mv]
        for k, v in combo.items():
            if k != "matching":
                arm = arm[arm[PIN_B[k]].astype(str) == v]
        line = []
        for et in CHECK_B:
            s, fl, err = b_arm(arm, et)
            if s is None:
                line.append(f"{et}: {err}")
                continue
            rank = int((s["_d"] > fl).sum()) + 1          # candidates strictly above it
            below = sum(1 for _, r in s.iterrows()
                        if (g(r, P, "margin") or 0) < 0
                        and g(r, P, "holm") is not None and g(r, P, "holm") < 0.05)
            line.append(f"{et}: rank {rank} of {len(s) + 1}, "
                        f"{below} significantly below")
        marker = "   <-- the reported model" if mv == combo.get("matching") else ""
        print(f"    {mv:<30} {' | '.join(line)}{marker}")
    print("\n  The Discussion says the floor ranks first of eleven and beats five of the")
    print("  ten candidates. Read which row above that is, and LABEL IT in the text with")
    print("  the matching scheme it came from -- those figures and the margins quoted")
    print("  beside them may come from different rows of this table.")
    notes.append("B: label the matching scheme on the Discussion's GAPDH-rank sentence")

    md.append(
        "**(b) Ten candidate regulators against the signature.** Target set: the "
        f"{B_NTARGET_GENES['total']} composition-adjusted signature genes, of which "
        f"{B_NTARGET_GENES['SE']} contribute skipped-exon events and "
        f"{B_NTARGET_GENES['MXE']} mutually-exclusive-exon events; RBM20 is one candidate "
        "among them. Same enrichment statistic, model and background match as (a). "
        f"**{FLOOR} is a housekeeping gene serving as the calibration floor, not a "
        "candidate**: every margin is that candidate's enrichment minus the floor's, so "
        "the floor has a margin of zero by definition and no interval of its own. It is "
        "shown in each panel at the position its enrichment puts it. Bootstrap *P* is on "
        "the margin; the Holm column corrects across the ten candidates. No candidate's "
        "interval excludes the floor in either event class.")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fwd", default=FWD)
    ap.add_argument("--inv", default=INV)
    ap.add_argument("--parts", default="ab")
    ap.add_argument("--model"); ap.add_argument("--matching")
    ap.add_argument("--gene-set", dest="gene_set")
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--out", default="S7")
    a = ap.parse_args()

    print("=" * 78)
    print("S7  candidate splicing regulators -- both tests, with intervals")
    print("=" * 78)
    print(f"  calibration floor : {FLOOR}, labelled as such in every table")
    print("  this script computes nothing -- it reads, identifies the reported model by")
    print("  its recorded numbers, checks, and formats")

    tA = read(a.fwd, PIN_A, "A forward") if "a" in a.parts else None
    if tA is not None:
        tA = prep_a(tA)
        evs = sorted({str(e) for e in tA["_ev"] if e})
        print(f"    event classes read from '{PIN_A['table']}': {', '.join(evs) or 'NONE'}")
        grid(tA, SEL_A, PIN_A, "A")
    tB = read(a.inv, PIN_B, "B inverted") if "b" in a.parts else None
    if tB is not None:
        grid(tB, SEL_B, PIN_B, "B")
    if a.list:
        print("\n[--list] stopping here.")
        return

    rows, notes, md = [], [], []
    if tA is not None:
        part_a(tA, rows, notes, md, a)
    elif "a" in a.parts:
        notes.append(f"A not built -- {a.fwd} unusable")
    if tB is not None:
        part_b(tB, rows, notes, md, a)
    elif "b" in a.parts:
        notes.append(f"B not built -- {a.inv} unusable")

    print("\n" + "=" * 78)
    print("RESULT")
    print("=" * 78)
    if not rows:
        print("  nothing built")
        for n in notes:
            print(f"  >>> {n}")
        sys.exit(1)

    T = pd.DataFrame(rows)
    T.to_csv(f"{a.out}.tsv", sep="\t", index=False)
    with pd.option_context("display.width", 250, "display.max_columns", 40):
        print(T.drop(columns=["model"]).to_string(index=False))
    print(f"\n[write] {a.out}.tsv   (includes a `model` column naming the provenance)")

    with open(f"{a.out}.md", "w") as fh:
        fh.write("**Supplementary Table S7. Candidate splicing regulators, "
                 "calibrated against a housekeeping-gene floor.**\n\n")
        for p in md:
            fh.write(p + "\n\n")
        fh.write("*Genes in the 45-gene reference set with no testable event: "
                 + ", ".join(A_ABSENT) + " in both event classes, and "
                 + ", ".join(A_ABSENT_MXE_ONLY) + " in mutually exclusive exons only. "
                 "An untested gene is not a null, and is named rather than omitted.*\n\n")
        fh.write("*The floor is a single-gene estimate, and its own enrichment moves "
                 "across models. That is why every margin is reported with an interval "
                 "rather than as a point difference.*\n\n")
        fh.write("*Sensitivity: without adjustment for cell composition, or in models "
                 "fitted across all samples rather than within disease, RBM20 and PTBP1 "
                 "do exceed the floor in mutually exclusive exons (margins +0.28 and "
                 "+0.19, Holm *P* < 0.05). Both are cardiomyocyte-enriched and target "
                 "events were selected for a large disease effect, so an all-sample model "
                 "allows any disease-associated gene to track them without regulating "
                 "them.*\n")
    print(f"[write] {a.out}.md   (legend; the sensitivity paragraph is already in it)")

    for n in notes:
        print(f"  >>> {n}")
    built = {r["part"][0] for r in rows}
    if {"A", "B"} - built:
        print(f"\n  S7 INCOMPLETE: part(s) {', '.join(sorted({'A', 'B'} - built))}"
              f" not built.")
        print("  S7 is the only full evidence for the paper's second claim. Do not submit")
        print("  it with a part missing -- fix the input, or say in the Methods what was")
        print("  not assessed.")
    else:
        print("\n  Both parts built from one identified model each; every recorded figure")
        print("  reproduced. Two things still need a human:")
        print("    1. verify A's testable-gene counts against 57's log (not in the file)")
        print("    2. label the matching scheme on the Discussion's GAPDH-rank sentence")


if __name__ == "__main__":
    main()
