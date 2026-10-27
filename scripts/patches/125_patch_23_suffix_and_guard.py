#!/usr/bin/env python3
"""
125_patch_23_suffix_and_guard.py -- make 23 read-length consistent, and stop it
overwriting the 50 bp comparators.

WHAT IS WRONG
-------------
`23_signature_rbp_scan.sbatch` is currently a THREE-WAY MIX:

    line  58  TGT=compadj_multi_event_type.txt            <- 50 bp, 134 genes
    line 126  M=rmats_magnet_rl125/${ET}.MATS.JCEC.txt    <- 125 bp
    line 138  --out-prefix signature_rbp_${ET}            <- UNSUFFIXED

110 renamed the input directory and left the gene list and the output prefix
alone. Submitted as it stands, the job would:

  1. score 125 bp events against the 50 bp gene list, and
  2. **overwrite** the 12 September 50 bp outputs:

        signature_rbp_SE.events.tsv    155 MB
        signature_rbp_SE.report.txt
        signature_rbp_MXE.events.tsv   187 MB
        signature_rbp_MXE.report.txt

Those four are the comparators `58_signature_rbp_margins.py` consumed on
5 October to build `signature_rbp_margins.tsv`. Overwriting them leaves that
table unreproducible from its own inputs, and **nothing errors.** Same shape as
`11_other_event_types.sbatch` writing the frozen comparators.

WHAT THIS PATCH DOES
--------------------
  * `SUF="${TBL_SUFFIX:-}"`, and derives all three from it:
        RMATS_DIR="rmats_magnet${SUF}"
        TGT="compadj${SUF}_multi_event_type.txt"
        --out-prefix signature_rbp_${ET}${SUF}
    so a run is consistent at ONE read length or it does not run.
  * default empty = the original 50 bp behaviour, reproducing the Sep 12 output
    (the envify convention; note 110 had already broken it for line 126).
  * a hard check that `$RMATS_DIR` exists — the existing per-class guard only
    `continue`s, which turns a wrong directory into a job that completes having
    done nothing.
  * **a don't-clobber guard**: refuses if any output it is about to write
    already exists, naming the file and its date, unless `OVERWRITE_OK=1`.

The expression side needs no change. `100_requant_magnet_125.sh` line 44 states
"The ONLY intended difference is --readLength 125", run on the raw STAR
alignments — so `magnet_expr.per_sample.tsv` (featureCounts, same BAMs) is valid
at both read lengths.

USAGE
-----
  python 125_patch_23_suffix_and_guard.py
  python 125_patch_23_suffix_and_guard.py --apply
  python 125_patch_23_suffix_and_guard.py --self-test

Then:
  TBL_SUFFIX=_rl125 sbatch 23_signature_rbp_scan.sbatch
"""
import argparse
import difflib
import os
import shutil
import sys
import time

TARGET = "23_signature_rbp_scan.sbatch"

A_CAND = ('CANDIDATES="QKI,RBM5,ZRANB2,HNRNPF,HNRNPA0,MBNL1,CELF1,SRSF1,'
          'PTBP1,GAPDH"\n')
B_CAND = '''# --------------------------------------------------------------------------
# Read length. Every path below derives from SUF, so a run is consistent at ONE
# read length or it does not run. Empty = the original 50 bp behaviour, which
# reproduces the 12 Sep output; _rl125 = the re-quantified tables.
#
# 110 renamed the rMATS directory on the M= line and left the gene list and the
# --out-prefix alone, which made this file a three-way mix whose 125 bp run
# would have overwritten the 50 bp comparators.
SUF="${TBL_SUFFIX:-}"
RMATS_DIR="${RMATS_DIR:-rmats_magnet${SUF}}"
# --------------------------------------------------------------------------
''' + A_CAND

A_TGT = "TGT=compadj_multi_event_type.txt"
B_TGT = 'TGT="${TGT:-compadj${SUF}_multi_event_type.txt}"'

A_M = "  M=rmats_magnet_rl125/${ET}.MATS.JCEC.txt\n"
B_M = "  M=${RMATS_DIR}/${ET}.MATS.JCEC.txt\n"

A_OUT = "    --out-prefix signature_rbp_${ET}\n"
B_OUT = "    --out-prefix signature_rbp_${ET}${SUF}\n"

A_LOOP = "for ET in SE MXE; do\n"
B_LOOP = '''# --------------------------------------------------------------------------
# A wrong rMATS directory must be an ERROR, not a skip. The per-class guard
# below only `continue`s, which turns a mistyped directory into a job that
# completes having done nothing and reports no failure.
if [ ! -d "$RMATS_DIR" ]; then
  echo "*** rMATS directory not found: $RMATS_DIR"
  echo "*** TBL_SUFFIX='$SUF'. Set RMATS_DIR explicitly if the directory for"
  echo "*** this read length has another name."
  exit 1
fi

# DON'T-CLOBBER GUARD.
# signature_rbp_{SE,MXE}.{events.tsv,report.txt} from 12 Sep are the 50 bp
# comparators that 58_signature_rbp_margins.py consumed to build
# signature_rbp_margins.tsv. Overwriting them makes that table unreproducible
# from its own inputs, and nothing errors.
for ET in SE MXE; do
  for EXT in events.tsv report.txt; do
    F="signature_rbp_${ET}${SUF}.${EXT}"
    if [ -e "$F" ] && [ -z "${OVERWRITE_OK:-}" ]; then
      echo "*** REFUSING — $F already exists"
      echo "***   modified: $(date -r "$F" '+%Y-%m-%d %H:%M' 2>/dev/null)"
      echo "***   size:     $(stat -c %s "$F" 2>/dev/null) bytes"
      echo "***"
      echo "*** Writing over it would destroy a comparator and nothing would"
      echo "*** error. If you mean to replace it, move it aside or re-run with"
      echo "***   OVERWRITE_OK=1 TBL_SUFFIX='$SUF' sbatch $0"
      exit 1
    fi
  done
done
echo "[guard] no output of signature_rbp_*${SUF} exists — safe to write"
# --------------------------------------------------------------------------

''' + A_LOOP

EDITS = [
    ("SUF and RMATS_DIR", A_CAND, B_CAND, 1),
    ("the target gene list, derived from SUF", A_TGT, B_TGT, 1),
    ("the rMATS input path, derived from RMATS_DIR", A_M, B_M, 1),
    ("the output prefix, suffixed", A_OUT, B_OUT, 1),
    ("the directory check and the don't-clobber guard", A_LOOP, B_LOOP, 1),
]


def apply_edits(text):
    problems, out = [], text
    for name, anchor, repl, want in EDITS:
        got = out.count(anchor)
        if got != want:
            problems.append(f"{name}: anchor occurs {got}x, expected {want}x\n"
                            f"      anchor: {anchor.strip().splitlines()[0][:80]}")
            continue
        if 'SUF="${TBL_SUFFIX:-}"' in out and name.startswith("SUF"):
            problems.append("SUF is already defined — already patched")
            continue
        out = out.replace(anchor, repl, want)
    return out, problems


def verify(text):
    bad = []
    if 'SUF="${TBL_SUFFIX:-}"' not in text:
        bad.append("SUF is not defined")
    if 'RMATS_DIR="${RMATS_DIR:-rmats_magnet${SUF}}"' not in text:
        bad.append("RMATS_DIR is not derived from SUF")
    if "rmats_magnet_rl125/" in text.replace("rmats_magnet${SUF}", ""):
        bad.append("a hard-coded rmats_magnet_rl125 path survives")
    if "--out-prefix signature_rbp_${ET}${SUF}" not in text:
        bad.append("the output prefix is not suffixed")
    if "TGT=compadj_multi_event_type.txt" in text:
        bad.append("the unparameterised gene-list name survives")
    if "OVERWRITE_OK" not in text:
        bad.append("the don't-clobber guard is missing")
    if 'if [ ! -d "$RMATS_DIR" ]; then' not in text:
        bad.append("the rMATS directory check is missing")
    # ordering: SUF defined before every consumer
    i = text.find('SUF="${TBL_SUFFIX:-}"')
    for m in ('RMATS_DIR="${RMATS_DIR:-rmats_magnet${SUF}}"',
              'TGT="${TGT:-compadj${SUF}_multi_event_type.txt}"',
              "--out-prefix signature_rbp_${ET}${SUF}",
              "OVERWRITE_OK"):
        j = text.find(m)
        if j < 0 or i > j:
            bad.append(f"SUF is not defined above its use: {m[:34]}")
    # the guard must precede the work loop that writes
    i_g = text.find("OVERWRITE_OK")
    i_w = text.find("--out-prefix signature_rbp_")
    if i_g < 0 or i_w < 0 or i_g > i_w:
        bad.append("the guard does not precede the writing loop")
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
    print("125 — make 23 read-length consistent, and guard its outputs")
    print("=" * 78)
    if not os.path.exists(a.target):
        print(f"[fatal] {a.target} not found.", file=sys.stderr)
        return 2
    src = open(a.target, errors="replace").read()
    new, problems = apply_edits(src)
    if problems:
        print("REFUSING — anchors do not match:")
        for p in problems:
            print(f"  ✗ {p}")
        print("\n  Nothing written. Send me lines 55-60, 118-142 and I will")
        print("  re-anchor.")
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
        print("Then check the guard fires BEFORE submitting anything:")
        print()
        print("  bash -n 23_signature_rbp_scan.sbatch      # syntax only")
        print("  sbatch 23_signature_rbp_scan.sbatch       # SUF empty -> must")
        print("                                            # REFUSE on the Sep")
        print("                                            # 12 files")
        print()
        print("  # then the real 125 bp run")
        print("  TBL_SUFFIX=_rl125 sbatch 23_signature_rbp_scan.sbatch")
        print()
        print("  It writes signature_rbp_{SE,MXE}_rl125.*, so the Sep 12")
        print("  comparators and signature_rbp_margins.tsv are untouched.")
        print()
        print("⚠ The comment on the TGT line still says '134 ... genes'. The")
        print("  count is computed into $NT and echoed, so the comment is the")
        print("  only stale copy — at 125 bp the list is 131. Edit by hand:")
        print()
        print("  grep -n 'the 134 multi-event-type' 23_signature_rbp_scan.sbatch")
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
    print("✓ verified on disk: one SUF drives the gene list, the rMATS")
    print("  directory and the output prefix; the guard precedes the loop.")
    return 0


def self_test():
    import subprocess
    import tempfile

    fixture = (
        "#!/bin/bash\n"
        + A_CAND
        + 'MAT=magnet_expr.log2cpm.tsv.gz\n'
        + A_TGT + "      # the 134 multi-event-type signature genes\n"
        + 'echo "TGT=$TGT"\n'
        + A_LOOP
        + A_M
        + '  echo "M=$M"\n'
        + '  echo stub \\\n'
        + A_OUT
        + "done\n")

    new, probs = apply_edits(fixture)
    assert not probs, probs
    assert not verify(new), verify(new)

    d = tempfile.mkdtemp()
    p = os.path.join(d, "fx.sh")
    open(p, "w").write(new)
    os.makedirs(os.path.join(d, "rmats_magnet"))
    os.makedirs(os.path.join(d, "rmats_magnet_rl125"))

    def run(env=None, args=()):
        return subprocess.run(["bash", p, *args], capture_output=True,
                              text=True, cwd=d,
                              env=dict(os.environ, **(env or {})))

    # 1. unsuffixed: 50 bp paths throughout
    r = run()
    assert r.returncode == 0, r.stdout + r.stderr
    assert "TGT=compadj_multi_event_type.txt" in r.stdout, r.stdout
    assert "M=rmats_magnet/SE.MATS.JCEC.txt" in r.stdout, r.stdout
    assert "--out-prefix signature_rbp_SE\n" in r.stdout, r.stdout
    assert "[guard]" in r.stdout, r.stdout

    # 2. suffixed: 125 bp paths throughout, including the OUTPUT
    r2 = run({"TBL_SUFFIX": "_rl125"})
    assert r2.returncode == 0, r2.stdout + r2.stderr
    assert "TGT=compadj_rl125_multi_event_type.txt" in r2.stdout, r2.stdout
    assert "M=rmats_magnet_rl125/SE.MATS.JCEC.txt" in r2.stdout, r2.stdout
    assert "--out-prefix signature_rbp_SE_rl125" in r2.stdout, r2.stdout

    # 3. an existing output must REFUSE
    open(os.path.join(d, "signature_rbp_SE.events.tsv"), "w").write("x")
    r3 = run()
    assert r3.returncode == 1, r3.stdout + r3.stderr
    assert "REFUSING" in r3.stdout, r3.stdout
    assert "signature_rbp_SE.events.tsv" in r3.stdout, r3.stdout
    # ...and must NOT have reached the loop
    assert "--out-prefix" not in r3.stdout, r3.stdout
    # the suffixed run is unaffected by the unsuffixed file
    r3b = run({"TBL_SUFFIX": "_rl125"})
    assert r3b.returncode == 0, r3b.stdout + r3b.stderr

    # 4. OVERWRITE_OK lets it through
    r4 = run({"OVERWRITE_OK": "1"})
    assert r4.returncode == 0, r4.stdout + r4.stderr
    assert "--out-prefix signature_rbp_SE\n" in r4.stdout, r4.stdout

    # 5. a missing rMATS directory is an ERROR, not a skip
    r5 = run({"TBL_SUFFIX": "_nope"})
    assert r5.returncode == 1, r5.stdout + r5.stderr
    assert "rMATS directory not found" in r5.stdout, r5.stdout

    # 6. double application refused
    _n, p6 = apply_edits(new)
    assert p6, "second application was not refused"

    print("[ok] self-test passed: one SUF drives gene list, rMATS directory and "
          "output prefix together; an existing output refuses before the loop "
          "runs; the suffixed run is unaffected by the unsuffixed file; "
          "OVERWRITE_OK releases it; a missing rMATS directory errors instead "
          "of skipping; double application refused")
    return 0


if __name__ == "__main__":
    sys.exit(main())
