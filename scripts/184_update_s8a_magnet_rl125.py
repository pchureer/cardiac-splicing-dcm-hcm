#!/usr/bin/env python3
"""
184_update_s8a_magnet_rl125.py -- Supplementary Table 8(a) at 125 bp.

S8a_permutations.tsv (74, 7 Oct) holds 20 Colorado analyses and 2 MAGNet
permutations. Only MAGNet was re-quantified, so only its two rows change; the
Colorado rows are the same runs at the same read length. Job 411394
(183_magnet_permutations_rl125.sbatch) re-ran the two MAGNet shuffles on the
125 bp table and wrote inflate_rl125/magnet_permutations_rl125.summary.txt.

This script copies S8a_permutations.tsv to S8a_permutations_rl125.tsv with:
  * the two MAGNet permutation rows replaced by the 125 bp run (job id, script,
    log, n_tested, n_significant, source line), same columns, same convention
  * one added row, run_kind=real, for the MAGNet positive control from the same
    job (944 of 80,973 full-coordinate), because a permutation null is only
    interpretable beside the real model from the same invocation -- which is
    exactly what 74 did for Colorado (478/475/1,596/1,605)

It refuses unless the original holds exactly the two MAGNet rows it expects
(seeds 101/102, 50 bp tested counts 103,054/103,073) and the summary holds two
permutations at sig=0 and a real model at 944. The original file is not written.

  python 184_update_s8a_magnet_rl125.py            # dry run: prints the rows
  python 184_update_s8a_magnet_rl125.py --apply    # writes S8a_permutations_rl125.tsv
"""

import argparse
import csv
import os
import re
import sys

SRC = "S8a_permutations.tsv"
OUT = "S8a_permutations_rl125.tsv"
SUMMARY = "inflate_rl125/magnet_permutations_rl125.summary.txt"
JOB = "411394"
SCRIPT = "183_magnet_permutations_rl125.sbatch"
LOG = f"magnet_perm_rl125_{JOB}.out"
EXPECT_50 = {"101": "103054", "102": "103073"}


def die(m):
    print("\n[REFUSE] " + m + "\n", file=sys.stderr)
    sys.exit(2)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--apply", action="store_true")
    a = ap.parse_args()
    for f in (SRC, SUMMARY):
        if not os.path.isfile(f):
            die(f"{f} missing")

    # ---- the 125 bp run
    txt = open(SUMMARY).read()
    real = re.search(r"real model\s+tested=(\d+)\s+sig=(\d+)", txt)
    perms = re.findall(r"permutation seed (\d+)\s+tested=(\d+)\s+sig=(\d+)", txt)
    exon = re.findall(r"\(exon level: tested=(\d+) sig=(\d+)\)", txt)
    if not real or len(perms) != 2:
        die(f"summary does not hold one real model and two permutations: {SUMMARY}")
    if real.group(2) != "944":
        die(f"real model sig={real.group(2)}, expected 944; do not update on a failed calibration")
    if any(s != "0" for _, _, s in perms):
        die("a permutation is non-zero; do not update the table until that is understood")
    print(f"125 bp: real {real.group(2)}/{real.group(1)}; "
          + "; ".join(f"seed {s} {sig}/{n}" for s, n, sig in perms)
          + (f"; exon level " + ", ".join(f"{sg}/{n}" for n, sg in exon) if exon else ""))

    # ---- the original
    with open(SRC, newline="") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        cols = rd.fieldnames
        rows = list(rd)
    need = ["job_id", "script", "log_file", "cohort", "contrast", "event_type", "run_kind",
            "seed", "n_tested", "n_significant", "pct_significant", "counting_convention",
            "significance_rule", "source_line"]
    miss = [c for c in need if c not in cols]
    if miss:
        die(f"{SRC} lacks columns {miss}")
    mag = [r for r in rows if r["cohort"] == "MAGNet"]
    if len(mag) != 2 or sorted(r["seed"] for r in mag) != ["101", "102"]:
        die(f"{SRC}: expected exactly two MAGNet rows with seeds 101/102, found "
            f"{[(r['seed'], r['run_kind']) for r in mag]}")
    for r in mag:
        if r["n_tested"] != EXPECT_50[r["seed"]] or r["n_significant"] != "0":
            die(f"{SRC}: MAGNet seed {r['seed']} reads {r['n_significant']}/{r['n_tested']}, "
                f"expected 0/{EXPECT_50[r['seed']]} (the 50 bp record). Not the file I expect.")
    print(f"50 bp rows to replace: " + "; ".join(f"seed {r['seed']} {r['n_significant']}/{r['n_tested']}" for r in mag))

    # ---- build the 125 bp rows
    lines = txt.splitlines()
    def src_line(pattern):
        for i, l in enumerate(lines, 1):
            if re.search(pattern, l):
                return str(i)
        return ""
    new_rows = []
    for r in rows:
        if r["cohort"] != "MAGNet":
            new_rows.append(r); continue
        n, sig = next((n, sg) for s, n, sg in perms if s == r["seed"])
        q = dict(r)
        q.update(job_id=JOB, script=SCRIPT, log_file=LOG, n_tested=n, n_significant=sig,
                 pct_significant="0.0", source_line=src_line(rf"permutation seed {r['seed']}\b"))
        new_rows.append(q)
    proto = dict(mag[0])
    proto.update(job_id=JOB, script=SCRIPT, log_file=LOG, contrast="MAGNet_DCMvsNF",
                 run_kind="real", seed="", n_tested=real.group(1), n_significant=real.group(2),
                 pct_significant=f"{100*int(real.group(2))/int(real.group(1)):.2f}",
                 source_line=src_line(r"real model"))
    new_rows.append(proto)

    print("\nMAGNet rows in the 125 bp table:")
    for r in new_rows:
        if r["cohort"] == "MAGNet":
            print("  " + "\t".join(r.get(c, "") for c in ("run_kind", "seed", "n_tested", "n_significant",
                                                        "counting_convention", "job_id", "source_line")))
    if not a.apply:
        print("\nDRY RUN. --apply writes", OUT)
        return
    if os.path.exists(OUT) and not os.environ.get("S8A_OVERWRITE"):
        die(f"{OUT} exists; remove it deliberately first, or set S8A_OVERWRITE=1 "
            f"(the 172 regeneration check does this: it backs the file up, then rebuilds in place)")
    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(new_rows)
    print(f"\n[write] {OUT}  ({len(new_rows)} rows; {SRC} untouched)")
    print("Next: 181 item for Supp Table 8(a) -> S8a_permutations_rl125.tsv; caption marker out.")


if __name__ == "__main__":
    main()
