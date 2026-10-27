#!/usr/bin/env python3
"""
46_colorado_icm_lookup.py -- settle the Colorado DCM-vs-ICM counts.

The one-pager says "0 of 20,027 SE, from 37 DCM hearts giving 460 SE events against NF".
Figure 3 says "476 against 2, on 19,953 / 19,998 tested", from 36 DCM. Post-exclusion
Colorado is 36 DCM / 13 ICM / 13 NF, so the figure is almost certainly right and the
text stale -- but "almost certainly" is not a citation.

Does it twice and compares:

  A  reads the frozen master_numbers_flagged.tsv, which is what Figure 3 draws
  B  recounts from the _dedup event tables using the project's own rules --
     exon-level collapse, FDR<0.05 and |dPSI|>=0.10

If A and B agree, the number is settled. If they disagree, the frozen file and the
tables have drifted, which matters far more than this one lookup.

Usage:
  export PATH=/data/apps/anaconda3/envs/mamba/envs/rnaseq_tools/bin:$PATH
  cd /data/project/chureeratp/RNA/HEART
  python 46_colorado_icm_lookup.py
"""
import argparse
import glob
import os

import pandas as pd

NKEEP = {"SE": 3, "MXE": 5, "A3SS": 5, "A5SS": 5, "RI": 3}
# Colorado MXE does not reproduce and is excluded project-wide
TYPES = ["SE", "A3SS", "A5SS", "RI"]


def exon_level(path, et, fdr, mind):
    """Recount a contrast the way the project counts: collapse, then threshold."""
    t = pd.read_csv(path, sep="\t")
    d = next((c for c in ("dPSI_adj", "dPSI") if c in t.columns), None)
    if d is None:
        cand = [c for c in t.columns if c.lower().startswith("dpsi")]
        d = next((c for c in cand if "adj" in c.lower()), cand[0] if cand else None)
    if d is None or "key" not in t.columns or "FDR" not in t.columns:
        return None, None
    n = NKEEP[et]
    t["ekey"] = [":".join(str(k).split(":")[:n]) for k in t["key"]]
    t = t.loc[t.groupby("ekey")["FDR"].idxmin()]
    sig = t[(t.FDR < fdr) & (t[d].abs() >= mind)]
    return len(sig), len(t)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--master", default="master_numbers_flagged.tsv")
    ap.add_argument("--fdr", type=float, default=0.05)
    ap.add_argument("--min-dpsi", type=float, default=0.10)
    a = ap.parse_args()

    print("=" * 74)
    print("A — from the frozen master file (what Figure 3 draws)")
    print("=" * 74)
    fromA = {}
    if os.path.exists(a.master) and os.path.getsize(a.master):
        M = pd.read_csv(a.master, sep="\t")
        C = M[M.cohort.astype(str).str.contains("Colorado", case=False, na=False)]
        if len(C):
            cols = [c for c in ("cohort", "contrast", "event_type", "sig_exon",
                                "tested_exon", "significant", "tested", "reliability")
                    if c in C.columns]
            print(C[cols].to_string(index=False))
            sc = "sig_exon" if "sig_exon" in C.columns else "significant"
            tc = "tested_exon" if "tested_exon" in C.columns else "tested"
            for _, r in C.iterrows():
                fromA[(str(r.get("contrast", "?")), str(r.get("event_type", "?")))] = \
                    (int(r[sc]), int(r[tc]))
        else:
            print("  no Colorado rows found")
    else:
        print(f"  {a.master} not found")

    print("\n" + "=" * 74)
    print("B — recounted from the _dedup tables (exon-level, "
          f"FDR<{a.fdr}, |dPSI|>={a.min_dpsi})")
    print("=" * 74)
    print(f"  {'contrast':<16}{'type':<7}{'significant':>12}{'tested':>10}   file")
    fromB = {}
    for label, pats in (("DCM vs ICM", ["{ET}_Colorado_DCMvsICM", "Colorado_{ET}_DCMvsICM"]),
                        ("DCM vs NF", ["{ET}_Colorado_DCMvsNF", "Colorado_{ET}_DCMvsNF"])):
        for et in TYPES:
            hit = None
            for pat in pats:
                for suf in ("_dedup.all_events.tsv", ".all_events.tsv"):
                    p = pat.format(ET=et) + suf
                    if os.path.exists(p) and os.path.getsize(p):
                        hit = p
                        break
                if hit:
                    break
            if not hit:
                continue
            s, n = exon_level(hit, et, a.fdr, a.min_dpsi)
            if s is None:
                print(f"  {label:<16}{et:<7}{'unreadable':>12}{'':>10}   {hit}")
                continue
            fromB[(label, et)] = (s, n)
            flag = "" if "_dedup" in hit else "   <- NOT the _dedup table"
            print(f"  {label:<16}{et:<7}{s:>12,}{n:>10,}   {os.path.basename(hit)}{flag}")
    if not fromB:
        print("  no Colorado contrast tables found. Look for what exists:")
        for f in sorted(glob.glob("*Colorado*all_events.tsv"))[:12]:
            print(f"    {f}")
        return

    print("\n" + "=" * 74)
    print("THE NUMBERS FOR THE TEXT")
    print("=" * 74)
    icm = fromB.get(("DCM vs ICM", "SE"))
    nf = fromB.get(("DCM vs NF", "SE"))
    if icm and nf:
        print(f"  Colorado SE, 36 DCM hearts, exon-level, FDR<{a.fdr} and "
              f"|dPSI|>={a.min_dpsi}:")
        print(f"     vs NF   {nf[0]:>6,} significant of {nf[1]:,} tested")
        print(f"     vs ICM  {icm[0]:>6,} significant of {icm[1]:,} tested")
        if icm[0]:
            print(f"     ratio   {nf[0] / icm[0]:.0f}-fold, same hearts, same pipeline")
        print("\n  Lead with the ratio. The one-pager's '0 of 20,027 from 37 DCM hearts")
        print("  giving 460' is pre-exclusion and should be replaced by these.")
    else:
        print("  one of the two SE contrasts is missing above — cannot state the ratio")

    if fromA and fromB:
        print("\n  A vs B, where both have a value:")
        bad = 0
        for k, (s, n) in fromB.items():
            m = next((v for kk, v in fromA.items()
                      if k[1] == kk[1] and k[0].replace(" ", "").lower()
                      in kk[0].replace(" ", "").lower()), None)
            if m is None:
                continue
            ok = (s, n) == m
            bad += 0 if ok else 1
            print(f"   [{'ok ' if ok else '***'}] {k[0]:<12}{k[1]:<6}"
                  f"master {m[0]:,}/{m[1]:,}   recount {s:,}/{n:,}")
        if bad:
            print(f"\n  *** {bad} disagreement(s). The frozen master file and the event")
            print("      tables have drifted. That is a bigger problem than this lookup:")
            print("      Figure 3 draws the master file, and the recount is the raw truth.")
        else:
            print("\n  [ok ] the frozen file and the tables agree — the number is settled.")


if __name__ == "__main__":
    main()
