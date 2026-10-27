#!/usr/bin/env python3
"""
151_jc_measurable_and_direction_split.py -- two manuscript numbers, from files
that already exist. No cluster job.

WHAT IT MEASURES
----------------
**A2** -- the manuscript's Methods [027]:

    "...agrees in direction for 97.1% of the 864 reported events JC can measure."

That conditions on the **JCEC-significant** set and asks how JC moves on it.
Script 63 conditions the other way round (direction agreement among
**JC**-significant events: 101/102 at 50 bp, 92/93 at 125 bp), so 63's figure is
NOT a substitute and must not be quoted for this.

**A6** -- Figure 2A's direction split: "584 with increased and 342 with decreased
inclusion in DCM", which sums to the 926 skipped exons. The 125 bp split of 813
is not recorded anywhere.

THE GATE
--------
The script runs 50 bp FIRST and refuses to print any 125 bp number unless the
50 bp pass reproduces all six published values:

    JCEC   926 significant of 62,971
    JC     715 significant of 50,124
    A2     864 measurable, 97.1% agreeing
    A6     584 up / 342 down

If the 50 bp pass does not reproduce those, the reader is wrong and the 125 bp
numbers it would print are worthless. That is the same discipline every patch
this month has used, and it is the whole value of running 50 bp at all.

CONVENTION -- COLLAPSE FIRST
----------------------------
Identical to `62_build_s5.py`'s `read_events`: truncate `key` to NKEEP[et]
fields, keep the lowest-FDR row per collapsed key, and only THEN apply
FDR < 0.05 and |dPSI| >= 0.10. Filtering before collapsing is the 926-vs-991
error (and 813-vs-864 at 125 bp -- see the warning below).

⚠ A COINCIDENCE THAT IS A TRAP
------------------------------
`864` appears twice in this project with two different meanings:

  * the manuscript's 864 = JCEC-significant events JC can measure, at **50 bp**
  * `06_sensitivity_concordance.py`'s 864 = the **filter-first** cassette-exon
    count of the **125 bp** base table (NUMBERS_125bp_REDERIVATION.md section 3,
    which records `by full coordinate key : 944` and `by cassette exon : 864`)

They are unrelated. Do not let one validate the other.

USAGE
-----
  cd /data/project/chureeratp/RNA/HEART
  python3 151_jc_measurable_and_direction_split.py              # 50 bp then 125 bp
  python3 151_jc_measurable_and_direction_split.py --self-test
  python3 151_jc_measurable_and_direction_split.py --only 50    # calibration only

Overrides, if a published value is ever re-recorded:
  CAL_JCEC_SIG=926 CAL_JCEC_TESTED=62971 CAL_JC_SIG=715 CAL_JC_TESTED=50124 \\
  CAL_MEASURABLE=864 CAL_AGREE_PCT=97.1 CAL_UP=584 CAL_DOWN=342
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

ET = "SE"
NKEEP = {"SE": 3, "MXE": 5, "A3SS": 5, "A5SS": 5, "RI": 3}
FDR_CUT = 0.05
DPSI_CUT = 0.10

# The published 50 bp values this must reproduce before 125 bp is believed.
CAL = {
    "jcec_sig": int(os.environ.get("CAL_JCEC_SIG", 926)),
    "jcec_tested": int(os.environ.get("CAL_JCEC_TESTED", 62971)),
    "jc_sig": int(os.environ.get("CAL_JC_SIG", 715)),
    "jc_tested": int(os.environ.get("CAL_JC_TESTED", 50124)),
    "measurable": int(os.environ.get("CAL_MEASURABLE", 864)),
    "agree_pct": float(os.environ.get("CAL_AGREE_PCT", 97.1)),
    "up": int(os.environ.get("CAL_UP", 584)),
    "down": int(os.environ.get("CAL_DOWN", 342)),
}

ARMS = {
    "50": dict(
        label="50 bp (as published)",
        jcec=["SE_v2filtered.all_events.tsv"],
        jc=["jccheck/SE_v2filtered_JC_dedup.all_events.tsv",
            "jccheck/SE_v2filtered_JC.all_events.tsv"],
    ),
    "125": dict(
        label="125 bp (primary)",
        jcec=["SE_v2filtered_rl125.all_events.tsv",
              "SE_rl125.all_events.tsv"],
        jc=["jccheck_rl125/SE_v2filtered_JC_dedup.all_events.tsv",
            "jccheck_rl125/SE_v2filtered_JC.all_events.tsv"],
    ),
}


def first_existing(paths):
    for p in paths:
        if os.path.exists(p) and os.path.getsize(p):
            return p
    return None


def read_events(path, et=ET):
    """Collapse-first, exactly as 62_build_s5.py's read_events does."""
    t = pd.read_csv(path, sep="\t", low_memory=False)
    if "key" not in t.columns or "FDR" not in t.columns:
        return None, f"no key/FDR column; have {list(t.columns)[:8]}"
    d = next((c for c in ("dPSI_adj", "dPSI") if c in t.columns), None)
    if d is None:
        cand = [c for c in t.columns if c.lower().startswith("dpsi")]
        d = next((c for c in cand if "adj" in c.lower()),
                 cand[0] if cand else None)
    if d is None:
        return None, "no dPSI column"
    t["ekey"] = [":".join(str(k).split(":")[:NKEEP[et]]) for k in t["key"]]
    t = t.loc[t.groupby("ekey")["FDR"].idxmin()]
    return t.rename(columns={d: "dpsi"})[["ekey", "dpsi", "FDR"]], None


def significant(T):
    return T[(T.FDR < FDR_CUT) & (T.dpsi.abs() >= DPSI_CUT)]


def measure(arm):
    """Everything for one read length, or (None, reason)."""
    spec = ARMS[arm]
    pe = first_existing(spec["jcec"])
    pj = first_existing(spec["jc"])
    if pe is None:
        return None, f"no JCEC discovery table; looked for {spec['jcec']}"
    if pj is None:
        return None, (f"no JC discovery table; looked for {spec['jc']}. "
                      f"Script 20 writes it")

    E, err = read_events(pe)
    if E is None:
        return None, f"JCEC table: {err}"
    J, err = read_events(pj)
    if J is None:
        return None, f"JC table: {err}"

    Es, Js = significant(E), significant(J)

    # ---- A2: the JCEC-significant events JC can MEASURE (present at all,
    # significant or not), and whether they move the same way.
    m = Es.merge(J, on="ekey", suffixes=("_ec", "_jc"))
    # A sign needs a direction to agree about; drop exact zeros rather than
    # scoring them, as 63 does.
    mm = m[(m.dpsi_ec != 0) & (m.dpsi_jc != 0)]
    same = int((np.sign(mm.dpsi_ec) == np.sign(mm.dpsi_jc)).sum())

    # ---- A6: the direction split of the JCEC-significant set
    up = int((Es.dpsi > 0).sum())
    down = int((Es.dpsi < 0).sum())

    return dict(
        jcec_path=pe, jc_path=pj,
        jcec_tested=len(E), jcec_sig=len(Es),
        jc_tested=len(J), jc_sig=len(Js),
        measurable=len(m), scored=len(mm), same=same,
        agree_pct=100.0 * same / len(mm) if len(mm) else float("nan"),
        unmeasurable=len(Es) - len(m),
        up=up, down=down, zero=len(Es) - up - down,
    ), None


def report(arm, r):
    print()
    print("-" * 78)
    print(f"  {ARMS[arm]['label']}")
    print("-" * 78)
    print(f"    JCEC table  {r['jcec_path']}")
    print(f"    JC   table  {r['jc_path']}")
    print()
    print(f"    JCEC  {r['jcec_sig']:>6} significant of {r['jcec_tested']:>7}")
    print(f"    JC    {r['jc_sig']:>6} significant of {r['jc_tested']:>7}")
    gain = (100.0 * (r['jcec_sig'] - r['jc_sig']) / r['jc_sig']
            if r['jc_sig'] else float('nan'))
    print(f"    JCEC advantage in significant events: {gain:+.1f}%")
    print()
    print("    A2  of the JCEC-significant events, how many can JC measure,")
    print("        and do they move the same way?")
    print(f"          measurable in JC      {r['measurable']:>6}"
          f"   (of {r['jcec_sig']}; {r['unmeasurable']} absent from the JC table)")
    print(f"          scored for direction  {r['scored']:>6}"
          f"   (exact zeros dropped)")
    print(f"          same direction        {r['same']:>6}"
          f"   = {r['agree_pct']:.1f}%")
    print()
    print("    A6  direction split of the JCEC-significant set (Figure 2A)")
    print(f"          increased inclusion   {r['up']:>6}")
    print(f"          decreased inclusion   {r['down']:>6}")
    if r["zero"]:
        print(f"          exactly zero          {r['zero']:>6}  ⚠ unexpected")
    print(f"          sum                   {r['up'] + r['down']:>6}"
          f"   (must equal {r['jcec_sig']})")


def calibrate(r):
    """Every published 50 bp value, or the list of mismatches."""
    bad = []
    for k, label in (("jcec_sig", "JCEC significant"),
                     ("jcec_tested", "JCEC tested"),
                     ("jc_sig", "JC significant"),
                     ("jc_tested", "JC tested"),
                     ("measurable", "measurable in JC"),
                     ("up", "increased inclusion"),
                     ("down", "decreased inclusion")):
        if r[k] != CAL[k]:
            bad.append(f"{label}: {r[k]} here, {CAL[k]} published")
    # the percentage is published to one decimal
    if abs(r["agree_pct"] - CAL["agree_pct"]) > 0.05:
        bad.append(f"direction agreement: {r['agree_pct']:.1f}% here, "
                   f"{CAL['agree_pct']}% published")
    if r["up"] + r["down"] != r["jcec_sig"]:
        bad.append(f"the direction split does not sum to the significant "
                   f"count: {r['up']} + {r['down']} != {r['jcec_sig']}")
    return bad


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", choices=["50", "125"], default=None,
                    help="run one arm only; 125 still requires the 50 bp gate "
                         "unless --skip-gate is given")
    ap.add_argument("--skip-gate", action="store_true",
                    help="report 125 bp without calibrating (say so if used)")
    ap.add_argument("--self-test", action="store_true")
    a = ap.parse_args()
    if a.self_test:
        return self_test()

    self_test()
    print()
    print("=" * 78)
    print("151 — the JC-measurable direction figure, and Figure 2A's split")
    print("=" * 78)
    print("  convention: COLLAPSE FIRST (key[:3] for SE, lowest FDR), then")
    print(f"              FDR < {FDR_CUT} and |dPSI| >= {DPSI_CUT}")

    cal_ok = None
    if a.only != "125":
        r50, err = measure("50")
        if r50 is None:
            print(f"\n  [fatal] 50 bp: {err}", file=sys.stderr)
            return 2
        report("50", r50)
        bad = calibrate(r50)
        print()
        if bad:
            cal_ok = False
            print("  >>> CALIBRATION FAILED. The 50 bp pass does not reproduce")
            print("      the published values:")
            for b in bad:
                print(f"        ✗ {b}")
            print("      The reader is wrong, so every 125 bp number it would")
            print("      print is suspect. Fix the reader first.")
        else:
            cal_ok = True
            print("  ✓ CALIBRATED — all eight published 50 bp values reproduce:")
            print(f"      JCEC {CAL['jcec_sig']}/{CAL['jcec_tested']}, "
                  f"JC {CAL['jc_sig']}/{CAL['jc_tested']},")
            print(f"      {CAL['measurable']} measurable at "
                  f"{CAL['agree_pct']}%, split {CAL['up']}/{CAL['down']}")

    if a.only == "50":
        return 0 if cal_ok else 3

    if cal_ok is False and not a.skip_gate:
        print("\n  REFUSING to report 125 bp: the 50 bp gate did not pass.")
        print("  Re-run with --skip-gate only if you intend to report an")
        print("  uncalibrated number, and say so wherever it is quoted.")
        return 3
    if cal_ok is None and not a.skip_gate:
        print("\n  REFUSING: --only 125 skips the gate. Use --skip-gate and")
        print("  say so, or drop --only.")
        return 3

    r125, err = measure("125")
    if r125 is None:
        print(f"\n  [fatal] 125 bp: {err}", file=sys.stderr)
        return 2
    report("125", r125)

    print()
    print("=" * 78)
    print("FOR THE MANUSCRIPT")
    print("=" * 78)
    if a.only != "125":
        print(f"  [027]  '{CAL['jc_sig']} significant skipped exons against "
              f"{CAL['jcec_sig']}'")
        print(f"         -> {r125['jc_sig']} against {r125['jcec_sig']}")
        print(f"         '(.{CAL['jc_tested']} against {CAL['jcec_tested']} "
              f"exons)'")
        print(f"         -> ({r125['jc_tested']} against "
              f"{r125['jcec_tested']} exons)")
        print(f"         'agrees in direction for {CAL['agree_pct']}% of the "
              f"{CAL['measurable']} reported events JC can measure'")
    print(f"         -> agrees in direction for {r125['agree_pct']:.1f}% of "
          f"the {r125['measurable']}")
    print()
    print(f"  Fig 2A 'increased ... decreased inclusion in DCM'")
    print(f"         -> {r125['up']} increased, {r125['down']} decreased")
    print()
    print("  ⚠ Do NOT quote script 63's 92/93 (98.92%) for the [027] figure.")
    print("    63 conditions on JC-significant events; [027] conditions on")
    print("    JCEC-significant ones. Different denominators, different claim.")
    print()
    print("  ⚠ '864' also appears in 06_sensitivity_concordance.py as the")
    print("    FILTER-FIRST cassette count of the 125 bp base table. Unrelated")
    print("    to the number above. Do not let one validate the other.")
    return 0


def self_test():
    """Synthetic tables with a known answer, through the real functions."""
    import tempfile

    d = tempfile.mkdtemp()
    rng = np.random.default_rng(3)

    # 40 events. Build the truth explicitly so every reported count is known.
    #   ekeys e00..e39, two rows each (so the collapse has work to do)
    #   JCEC significant: e00..e19  (12 up, 8 down)
    #   present in JC:    e00..e14  -> 15 measurable
    #   of those, JC flips the sign on e13 and e14 -> 13/15 agree
    rows_e, rows_j = [], []
    for i in range(40):
        sig = i < 20
        up = sig and i < 12
        dp = (0.4 if up else -0.4) if sig else 0.02
        fdr = 0.001 if sig else 0.5
        # two rows per collapsed key; the lower-FDR one must win
        rows_e.append(dict(key=f"chr1:{i}:{i+9}:a:b", FDR=fdr, dPSI_adj=dp))
        rows_e.append(dict(key=f"chr1:{i}:{i+9}:c:d", FDR=0.9, dPSI_adj=-dp))
        if i < 15:
            flip = i in (13, 14)
            rows_j.append(dict(key=f"chr1:{i}:{i+9}:a:b", FDR=fdr,
                               dPSI_adj=(-dp if flip else dp)))
    pe = os.path.join(d, "SE_v2filtered.all_events.tsv")
    pd.DataFrame(rows_e).to_csv(pe, sep="\t", index=False)
    os.makedirs(os.path.join(d, "jccheck"), exist_ok=True)
    pj = os.path.join(d, "jccheck/SE_v2filtered_JC.all_events.tsv")
    pd.DataFrame(rows_j).to_csv(pj, sep="\t", index=False)

    cwd = os.getcwd()
    os.chdir(d)
    try:
        r, err = measure("50")
        assert r is not None, err
        # collapse-first: 40 collapsed keys from 80 rows, not 80
        assert r["jcec_tested"] == 40, r
        assert r["jcec_sig"] == 20, r
        assert r["up"] == 12 and r["down"] == 8, r
        assert r["up"] + r["down"] == r["jcec_sig"], r
        # A2: 15 of the 20 significant events exist in the JC table
        assert r["measurable"] == 15, r
        assert r["unmeasurable"] == 5, r
        assert r["same"] == 13, r
        assert abs(r["agree_pct"] - 100 * 13 / 15) < 1e-9, r
        # the JC table's own significant count is independent of the above
        assert r["jc_sig"] == 15, r

        # the lowest-FDR row must be the one kept: if the collapse took the
        # other row, every sign would be inverted and `same` would be 2
        assert r["same"] != 2, "the collapse kept the wrong row"

        # the gate must FAIL on this fixture, since it is not the real data
        bad = calibrate(r)
        assert bad, "calibration passed on synthetic data"
        assert any("JCEC significant" in b for b in bad), bad
    finally:
        os.chdir(cwd)

    # a table with no dPSI column is reported, not crashed on
    p2 = os.path.join(d, "bad.tsv")
    pd.DataFrame([dict(key="chr1:1:2:a:b", FDR=0.01)]).to_csv(
        p2, sep="\t", index=False)
    T, err = read_events(p2)
    assert T is None and "dPSI" in err, (T, err)

    # exact zeros are dropped from the direction score, not counted as agreeing
    rows_z = [dict(key="chr1:1:2:a:b", FDR=0.001, dPSI_adj=0.4),
              dict(key="chr1:3:4:a:b", FDR=0.001, dPSI_adj=0.4)]
    rows_zj = [dict(key="chr1:1:2:a:b", FDR=0.001, dPSI_adj=0.0),
               dict(key="chr1:3:4:a:b", FDR=0.001, dPSI_adj=0.4)]
    d2 = tempfile.mkdtemp()
    pd.DataFrame(rows_z).to_csv(
        os.path.join(d2, "SE_v2filtered.all_events.tsv"), sep="\t", index=False)
    os.makedirs(os.path.join(d2, "jccheck"), exist_ok=True)
    pd.DataFrame(rows_zj).to_csv(
        os.path.join(d2, "jccheck/SE_v2filtered_JC.all_events.tsv"),
        sep="\t", index=False)
    os.chdir(d2)
    try:
        r2, err = measure("50")
        assert r2["measurable"] == 2, r2
        assert r2["scored"] == 1, r2          # the zero is dropped
        assert r2["same"] == 1, r2
        assert r2["agree_pct"] == 100.0, r2
    finally:
        os.chdir(cwd)

    print("[ok] self-test passed: collapse-first keeps the lowest-FDR row (80 "
          "rows -> 40 keys, and a deliberately inverted high-FDR duplicate "
          "does not win); the JCEC-significant set is 20 with a 12/8 split "
          "summing to it; A2 counts 15 of those 20 as measurable in JC with 5 "
          "absent and 13 agreeing; the JC table's own significant count is "
          "computed independently; exact zeros are dropped from the direction "
          "score rather than scored as agreeing; a missing dPSI column is "
          "reported not crashed on; and the calibration gate FAILS on "
          "synthetic data instead of passing vacuously")
    return 0


if __name__ == "__main__":
    sys.exit(main())
