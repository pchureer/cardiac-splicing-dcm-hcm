#!/usr/bin/env python3
"""
155_patch_65_90_readlength_records.py -- S7 and Figure 5 at 125 bp.

TWO FILES, ALL OR NOTHING: `65_build_s7.py` and `90_build_f5.py`.

WHY
---
65 already has `--fwd`, `--inv` and `--out`, so pointing it at the 125 bp
deposits needs no patch. It still **refused**:

    >>> B NOT BUILT. None of the 8 combinations of model, matching reproduces
        the recorded figures.
            3 x  SE: top candidate by margin is RBM20, not QKI
            2 x  SE: top candidate by margin is PTBP1, not QKI
        Either this is a different file, or the record is wrong. Do not
        pick a combination by eye to make the table appear.

That is correct behaviour. 65 does not merely read the file — it **identifies
which of eight model combinations is the reported one** by reproducing recorded
headline figures. Those records are 50 bp. At 125 bp they all move.

Tenth and eleventh instances of the hard-coded-statistic pattern, after 68, 95,
96, 36, 93, 56's G45 list, 49, 58 and 56's reproduction gate.

⚠ THE ONE FAILURE MODE WORSE THAN REFUSING
------------------------------------------
`CHECK_B` keys on a gene **NAME**:

    CHECK_B = {"SE":  dict(top="QKI",   margin=+0.033, ci=(-0.025, 0.082)),
               "MXE": dict(top="RBM20", margin=+0.009, ci=(-0.053, 0.088))}

At 125 bp the identity changes — **QKI -> PTBP1** in SE and **RBM20 -> PTBP1**
in MXE. A gate told "top is QKI" and handed PTBP1's numbers would **PASS while
naming the wrong gene**, and S7 would then report a model identified on a false
premise. So the overrides are **MANDATORY when a suffixed input is given**,
never defaulted — the same design as 132's and 138's gates, and the same
`_rl<digits>` detection as 131's.

THE 125 bp RECORDS, AND WHERE EACH COMES FROM
---------------------------------------------
| | 50 bp (default) | 125 bp | source |
|---|---|---|---|
| CHECK_A SE | -0.005, (-0.119, 0.080), p 0.76, n 663 | **-0.051, (-0.159, 0.047), p 0.311, n 641** | 57's rl125 deposit; n from RESULT_56 |
| CHECK_A MXE | +0.000, (-0.085, 0.104), p 0.88, n 1393 | **-0.019, (-0.108, 0.072), p 0.690, n 1399** | same |
| CHECK_B SE | top QKI, +0.033, (-0.025, 0.082) | **top PTBP1, +0.062, (-0.025, 0.122)** | 152, SE DCM/cell 3-D |
| CHECK_B MXE | top RBM20, +0.009, (-0.053, 0.088) | **top PTBP1, +0.113, (+0.016, +0.187)** | 152, MXE DCM/cell 3-D |
| B_NTARGET_GENES | 134 / 130 / 124 | **131 / 123 / 119** | RESULT_23_58 |
| REC_FLOOR_RANK | SE 5, MXE 2 | **SE 6, MXE 3** | 152 |
| A_TESTABLE | SE 40, MXE 39 | **unchanged** | RESULT_56 |

Every 50 bp default was verified today: 65 reproduces `S7.tsv` byte-identically,
57 reproduces its own 5 Oct deposit byte-identically, and 152 reproduces the
published floor ranks (MXE 2 of 11, SE 5 of 11).

90 ALSO NEEDS ITS INPUT PATHS FREED
-----------------------------------
`SCAN` and `FWD` are module literals used via `os.path.join(d, SCAN)` (line 149)
and `os.path.join(a.dir, FWD)` (205). `--dir` cannot separate the read lengths
because both deposits live in the same directory. Rather than touch those call
sites, the constants themselves read the environment — no argparse change, no
change to either function.

USAGE
-----
  python 155_patch_65_90_readlength_records.py [--apply] [--self-test]
"""
import argparse
import difflib
import os
import shutil
import subprocess
import sys
import time

T65 = "65_build_s7.py"
T90 = "90_build_f5.py"

# ------------------------------------------------------------------ 65
A_65_CHECK = '''CHECK_A = {"SE":  dict(margin=-0.005, ci=(-0.119, 0.080), p=0.76, n_target=663),
           "MXE": dict(margin=+0.000, ci=(-0.085, 0.104), p=0.88, n_target=1393)}
CHECK_B = {"SE":  dict(top="QKI",   margin=+0.033, ci=(-0.025, 0.082)),
           "MXE": dict(top="RBM20", margin=+0.009, ci=(-0.053, 0.088))}
'''

B_65_CHECK = '''# READ LENGTH, added 22 Oct. Everything below is a 50 bp record. At 125 bp every
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
    if not _re.search(r"_rl\\d+", _argv):
        return
    _miss = [_n for _n in ("S7_CHECK_A", "S7_CHECK_B", "S7_B_NGENES")
             if not os.environ.get(_n, "").strip()]
    if _miss:
        sys.exit(
            "REFUSING: an input path carries _rl<digits> but %s %s unset.\\n"
            "  These are CHECKS that identify the reported model, and the 50 bp\\n"
            "  records would either refuse every combination or -- worse --\\n"
            "  match one while naming the wrong top gene (QKI/RBM20 at 50 bp,\\n"
            "  PTBP1 at 125 bp). Supply all three:\\n"
            "    S7_CHECK_A=\\"SE=-0.051/-0.159/0.047/0.311/641,"
            "MXE=-0.019/-0.108/0.072/0.690/1399\\" \\\\\\n"
            "    S7_CHECK_B=\\"SE=PTBP1/0.062/-0.025/0.122,"
            "MXE=PTBP1/0.113/0.016/0.187\\" \\\\\\n"
            "    S7_B_NGENES=\\"total=131,SE=123,MXE=119\\""
            % (", ".join(_miss), "is" if len(_miss) == 1 else "are"))


_refuse_unparameterised()
'''

A_65_NG = 'B_NTARGET_GENES = {"total": 134, "SE": 130, "MXE": 124}\n'
B_65_NG = '''_G_ENV = _s7_env("S7_B_NGENES")
B_NTARGET_GENES = ({_k: int(_v[0]) for _k, _v in _G_ENV.items()} if _G_ENV
                   else {"total": 134, "SE": 130, "MXE": 124})
'''

# ------------------------------------------------------------------ 90
A_90_FILES = '''SCAN = "signature_rbp_margins.tsv"
FWD = "rbm20_margin_resolve.tsv"
'''
B_90_FILES = '''# READ LENGTH, added 22 Oct. These are joined with --dir at lines ~149 and ~205,
# and --dir cannot separate the read lengths because both deposits sit in the
# same directory. Read from the environment so neither call site changes:
#   F5_SCAN=signature_rbp_margins_rl125.tsv
#   F5_FWD=rbm20_margin_resolve_rl125.tsv
SCAN = os.environ.get("F5_SCAN", "signature_rbp_margins.tsv")
FWD = os.environ.get("F5_FWD", "rbm20_margin_resolve.tsv")
'''

A_90_REC = '''REC_SCAN = {("SE", "QKI"): (0.033, -0.025, 0.082),
            ("MXE", "RBM20"): (0.009, -0.053, 0.088)}
REC_FLOOR_RANK = {"SE": 5, "MXE": 2}          # of 11 genes, by delta descending
REC_FWD = {"SE": (-0.005, -0.119, 0.080), "MXE": (0.000, -0.085, 0.104)}
'''

B_90_REC = '''# READ LENGTH. All three are 50 bp records, and REC_SCAN is keyed by
# (class, GENE) -- the gene changes identity at 125 bp (QKI -> PTBP1 in SE,
# RBM20 -> PTBP1 in MXE), so a gate handed 125 bp numbers under 50 bp keys would
# name the wrong gene in a figure. Mandatory when F5_SCAN/F5_FWD are suffixed.
#   F5_REC_SCAN="SE=PTBP1/0.062/-0.025/0.122,MXE=PTBP1/0.113/0.016/0.187"
#   F5_FLOOR_RANK="SE=6,MXE=3"
#   F5_REC_FWD="SE=-0.051/-0.159/0.047,MXE=-0.019/-0.108/0.072"
def _f5_env(_name):
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


_S_ENV = _f5_env("F5_REC_SCAN")
REC_SCAN = ({(_et, _v[0]): (float(_v[1]), float(_v[2]), float(_v[3]))
             for _et, _v in _S_ENV.items()} if _S_ENV else
            {("SE", "QKI"): (0.033, -0.025, 0.082),
             ("MXE", "RBM20"): (0.009, -0.053, 0.088)})

_R_ENV = _f5_env("F5_FLOOR_RANK")
REC_FLOOR_RANK = ({_k: int(_v[0]) for _k, _v in _R_ENV.items()} if _R_ENV
                  else {"SE": 5, "MXE": 2})   # of 11 genes, by delta descending

_W_ENV = _f5_env("F5_REC_FWD")
REC_FWD = ({_et: (float(_v[0]), float(_v[1]), float(_v[2]))
            for _et, _v in _W_ENV.items()} if _W_ENV else
           {"SE": (-0.005, -0.119, 0.080), "MXE": (0.000, -0.085, 0.104)})

if ("_rl" in SCAN or "_rl" in FWD):
    _miss = [_n for _n in ("F5_REC_SCAN", "F5_FLOOR_RANK", "F5_REC_FWD")
             if not os.environ.get(_n, "").strip()]
    if _miss:
        sys.exit(
            "REFUSING: F5_SCAN/F5_FWD name a suffixed deposit but %s unset.\\n"
            "  REC_SCAN is keyed by GENE and the top gene changes at 125 bp\\n"
            "  (QKI/RBM20 -> PTBP1), so the figure would be labelled with the\\n"
            "  wrong gene rather than fail. Supply all three."
            % ", ".join(_miss))
'''

EDITS_65 = [("65: the recorded headlines", A_65_CHECK, B_65_CHECK, 1),
            ("65: B_NTARGET_GENES", A_65_NG, B_65_NG, 1)]
EDITS_90 = [("90: the input paths", A_90_FILES, B_90_FILES, 1),
            ("90: the record", A_90_REC, B_90_REC, 1)]


def apply_edits(text, edits, already):
    if already in text:
        return text, [f"already carries 155 ({already[:40]}...)"]
    problems, out = [], text
    for name, anchor, repl, want in edits:
        got = out.count(anchor)
        if got != want:
            problems.append(f"{name}: anchor occurs {got}x, expected {want}x\n"
                            f"      first line: {anchor.splitlines()[0][:66]}")
            continue
        out = out.replace(anchor, repl, want)
    return out, problems


def verify_65(text):
    bad = []
    for w in ('def _s7_env(', '_A_ENV = _s7_env("S7_CHECK_A")',
              '_B_ENV = _s7_env("S7_CHECK_B")',
              '_G_ENV = _s7_env("S7_B_NGENES")',
              'def _refuse_unparameterised():', '_refuse_unparameterised()',
              'naming the wrong top gene'):
        if w not in text:
            bad.append(f"missing: {w[:52]}")
    # the 50 bp defaults must survive verbatim
    for w in ('dict(margin=-0.005, ci=(-0.119, 0.080), p=0.76, n_target=663)',
              'dict(top="QKI",   margin=+0.033, ci=(-0.025, 0.082))',
              'dict(top="RBM20", margin=+0.009, ci=(-0.053, 0.088))',
              '{"total": 134, "SE": 130, "MXE": 124}'):
        if w not in text:
            bad.append(f"a 50 bp default was lost: {w[:52]}")
    # A_TESTABLE must NOT be parameterised — it does not change
    if 'A_TESTABLE = {"SE": 40, "MXE": 39}' not in text:
        bad.append("A_TESTABLE was altered; RESULT_56 has it unchanged at "
                   "40/39 and it must stay a literal")
    return bad


def verify_90(text):
    bad = []
    for w in ('SCAN = os.environ.get("F5_SCAN"',
              'FWD = os.environ.get("F5_FWD"',
              'def _f5_env(', '_S_ENV = _f5_env("F5_REC_SCAN")',
              '_R_ENV = _f5_env("F5_FLOOR_RANK")',
              '_W_ENV = _f5_env("F5_REC_FWD")',
              'wrong gene rather than fail'):
        if w not in text:
            bad.append(f"missing: {w[:52]}")
    for w in ('("SE", "QKI"): (0.033, -0.025, 0.082)',
              '("MXE", "RBM20"): (0.009, -0.053, 0.088)',
              '{"SE": 5, "MXE": 2}',
              '{"SE": (-0.005, -0.119, 0.080), "MXE": (0.000, -0.085, 0.104)}'):
        if w not in text:
            bad.append(f"a 50 bp default was lost: {w[:52]}")
    # CANDIDATES is the panel's gene list, not a record — must not move
    if 'CANDIDATES = ["RBM20", "QKI", "RBM5", "ZRANB2", "HNRNPF", "HNRNPA0",' \
            not in text:
        bad.append("CANDIDATES was altered; it is the panel's gene list, not a "
                   "recorded statistic")
    return bad


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--apply", action="store_true")
    ap.add_argument("--self-test", action="store_true")
    a = ap.parse_args()
    if a.self_test:
        return self_test()

    self_test()
    print()
    print("=" * 78)
    print("155 — S7 and Figure 5: the records follow the read length")
    print("=" * 78)

    for t in (T65, T90):
        if not os.path.exists(t):
            print(f"[fatal] {t} not found.", file=sys.stderr)
            return 2

    s65 = open(T65, errors="replace").read()
    s90 = open(T90, errors="replace").read()
    n65, p65 = apply_edits(s65, EDITS_65, "def _s7_env(")
    n90, p90 = apply_edits(s90, EDITS_90, "def _f5_env(")
    if p65 or p90:
        print("REFUSING — anchors do not match. NEITHER file is written:")
        for p in p65 + p90:
            print(f"  ✗ {p}")
        return 3
    b65, b90 = verify_65(n65), verify_90(n90)
    if b65 or b90:
        print("REFUSING — patched text fails verification:")
        for b in b65 + b90:
            print(f"  ✗ {b}")
        return 3
    import tempfile
    d = tempfile.mkdtemp()
    for nm, txt in ((T65, n65), (T90, n90)):
        p = os.path.join(d, nm)
        open(p, "w").write(txt)
        r = subprocess.run([sys.executable, "-c",
                            f"import ast;ast.parse(open({p!r}).read())"],
                           capture_output=True, text=True)
        if r.returncode != 0:
            print(f"REFUSING — {nm} is not valid python:\n{r.stderr}")
            return 3

    for nm, old, new in ((T65, s65, n65), (T90, s90, n90)):
        print("".join(difflib.unified_diff(
            old.splitlines(True), new.splitlines(True),
            fromfile=nm, tofile=nm + " (patched)", n=2)))
    print("=" * 78)

    if not a.apply:
        print("DRY RUN — nothing written. Re-run with --apply, then:")
        print()
        print("  # 1. 65 at 50 bp must stay byte-identical to the published S7")
        print("  python3 65_build_s7.py --out S7_55check 2>&1 | tail -5")
        print("  cmp S7.tsv S7_55check.tsv && echo IDENTICAL")
        print()
        print("  # 2. the refusal, on purpose: suffixed input, no records")
        print("  python3 65_build_s7.py --fwd rbm20_margin_resolve_rl125.tsv \\")
        print("      --inv signature_rbp_margins_rl125.tsv --out /tmp/x"
              "  # must REFUSE")
        print()
        print("  # 3. S7 at 125 bp")
        print('  S7_CHECK_A="SE=-0.051/-0.159/0.047/0.311/641,'
              'MXE=-0.019/-0.108/0.072/0.690/1399" \\')
        print('  S7_CHECK_B="SE=PTBP1/0.062/-0.025/0.122,'
              'MXE=PTBP1/0.113/0.016/0.187" \\')
        print('  S7_B_NGENES="total=131,SE=123,MXE=119" \\')
        print("  python3 65_build_s7.py --fwd rbm20_margin_resolve_rl125.tsv \\")
        print("      --inv signature_rbp_margins_rl125.tsv --out S7_rl125 "
              "2>&1 | tail -30")
        print()
        print("  # 4. Figure 5 at 50 bp, then 125 bp")
        print("  python3 90_build_f5.py --list 2>&1 | tail -20")
        print("  F5_SCAN=signature_rbp_margins_rl125.tsv \\")
        print("  F5_FWD=rbm20_margin_resolve_rl125.tsv \\")
        print('  F5_REC_SCAN="SE=PTBP1/0.062/-0.025/0.122,'
              'MXE=PTBP1/0.113/0.016/0.187" \\')
        print('  F5_FLOOR_RANK="SE=6,MXE=3" \\')
        print('  F5_REC_FWD="SE=-0.051/-0.159/0.047,'
              'MXE=-0.019/-0.108/0.072" \\')
        print("  python3 90_build_f5.py --out legend_check/f5_rl125 "
              "2>&1 | tail -25")
        print()
        print("⚠ STEP 2 IS NOT OPTIONAL. It is the one case worse than a")
        print("  failure: CHECK_B and REC_SCAN key on a gene NAME, and at")
        print("  125 bp the top gene changes (QKI/RBM20 -> PTBP1). Without")
        print("  the refusal a run could match a model while labelling it")
        print("  with the wrong gene.")
        return 0

    stamp = time.strftime("%Y%m%d_%H%M%S")
    baks = []
    for nm, txt in ((T65, n65), (T90, n90)):
        bak = f"{nm}.bak.{stamp}"
        shutil.copy2(nm, bak)
        baks.append((nm, bak))
        open(nm, "w").write(txt)
        print(f"[write] {nm}   (backup: {bak})")
    ok = True
    for nm, _bak in baks:
        disk = open(nm, errors="replace").read()
        bad = verify_65(disk) if nm == T65 else verify_90(disk)
        if bad:
            ok = False
            print(f"✗ {nm} failed on-disk verification:")
            for b in bad:
                print(f"    {b}")
    if not ok:
        for nm, bak in baks:
            shutil.copy2(bak, nm)
        print("restored both files")
        return 3
    print("✓ verified on disk: both scripts' records follow the read length,")
    print("  every 50 bp default survives verbatim, A_TESTABLE and CANDIDATES")
    print("  are untouched, and a suffixed input without records REFUSES")
    print("  rather than naming the wrong gene.")
    return 0


def self_test():
    import tempfile

    fx65 = ("import os\nimport sys\n"
            'SEL_A = ["table_variant", "model", "gene_set", "matching"]\n'
            + A_65_CHECK +
            "CHECK_B_NCAND = 10\n"
            + 'A_TESTABLE = {"SE": 40, "MXE": 39}\n'
            + A_65_NG +
            "TOL_MARGIN, TOL_CI, TOL_P = 0.004, 0.010, 0.03\n"
            'print("A", CHECK_A["SE"]["margin"], CHECK_A["SE"]["n_target"])\n'
            'print("B", CHECK_B["SE"]["top"], CHECK_B["MXE"]["top"])\n'
            'print("G", B_NTARGET_GENES)\n')

    fx90 = ("import os\nimport sys\n"
            + A_90_FILES +
            'FLOOR_GENE = "GAPDH"\n'
            'CANDIDATES = ["RBM20", "QKI", "RBM5", "ZRANB2", "HNRNPF", '
            '"HNRNPA0",\n              "MBNL1", "CELF1", "SRSF1", "PTBP1"]\n'
            + A_90_REC +
            "TOL = 0.004\n"
            'print("SCAN", SCAN)\nprint("FWD", FWD)\n'
            'print("REC", sorted(REC_SCAN))\nprint("RANK", REC_FLOOR_RANK)\n'
            'print("FWDREC", REC_FWD["SE"])\n')

    n65, p65 = apply_edits(fx65, EDITS_65, "def _s7_env(")
    assert not p65, p65
    assert not verify_65(n65), verify_65(n65)
    n90, p90 = apply_edits(fx90, EDITS_90, "def _f5_env(")
    assert not p90, p90
    assert not verify_90(n90), verify_90(n90)
    compile(n65, "fx65", "exec")
    compile(n90, "fx90", "exec")

    d = tempfile.mkdtemp()

    def run(src, name, env=None, argv=()):
        p = os.path.join(d, name)
        open(p, "w").write(src)
        e = dict(os.environ)
        for k in ("S7_CHECK_A", "S7_CHECK_B", "S7_B_NGENES",
                  "F5_SCAN", "F5_FWD", "F5_REC_SCAN", "F5_FLOOR_RANK",
                  "F5_REC_FWD"):
            e.pop(k, None)
        e.update(env or {})
        return subprocess.run([sys.executable, "-I", p, *argv],
                              capture_output=True, text=True, cwd=d, env=e)

    # --- 65 unsuffixed: identical to pre-patch
    r0 = run(fx65, "p65.py")
    r1 = run(n65, "n65.py")
    assert r0.returncode == 0 and r1.returncode == 0, (r0.stderr, r1.stderr)
    assert r0.stdout == r1.stdout, (r0.stdout, r1.stdout)
    assert "B QKI RBM20" in r1.stdout, r1.stdout

    # --- 65 with records supplied
    r2 = run(n65, "n65.py", {
        "S7_CHECK_A": "SE=-0.051/-0.159/0.047/0.311/641,"
                      "MXE=-0.019/-0.108/0.072/0.690/1399",
        "S7_CHECK_B": "SE=PTBP1/0.062/-0.025/0.122,"
                      "MXE=PTBP1/0.113/0.016/0.187",
        "S7_B_NGENES": "total=131,SE=123,MXE=119"})
    assert r2.returncode == 0, r2.stderr
    assert "A -0.051 641" in r2.stdout, r2.stdout
    assert "B PTBP1 PTBP1" in r2.stdout, r2.stdout
    assert "'total': 131" in r2.stdout and "'SE': 123" in r2.stdout, r2.stdout

    # --- THE GUARD: suffixed argv, no records -> refuse
    r3 = run(n65, "n65.py", argv=("--fwd", "rbm20_margin_resolve_rl125.tsv"))
    assert r3.returncode != 0, r3.stdout
    assert "REFUSING" in r3.stderr + r3.stdout, (r3.stdout, r3.stderr)
    assert "naming the wrong top gene" in r3.stderr + r3.stdout, r3.stderr

    # --- a PARTIAL override is still refused
    r4 = run(n65, "n65.py", {"S7_CHECK_A": "SE=-0.051/-0.159/0.047/0.311/641"},
             argv=("--fwd", "x_rl125.tsv"))
    assert r4.returncode != 0, r4.stdout
    assert "S7_CHECK_B" in r4.stderr + r4.stdout, r4.stderr

    # --- unsuffixed argv with no records must NOT refuse
    r5 = run(n65, "n65.py", argv=("--fwd", "rbm20_margin_resolve.tsv"))
    assert r5.returncode == 0, (r5.stdout, r5.stderr)

    # --- 90 unsuffixed: identical to pre-patch
    q0 = run(fx90, "p90.py")
    q1 = run(n90, "n90.py")
    assert q0.returncode == 0 and q1.returncode == 0, (q0.stderr, q1.stderr)
    assert q0.stdout == q1.stdout, (q0.stdout, q1.stdout)
    assert "RANK {'SE': 5, 'MXE': 2}" in q1.stdout, q1.stdout

    # --- 90 with records
    q2 = run(n90, "n90.py", {
        "F5_SCAN": "signature_rbp_margins_rl125.tsv",
        "F5_FWD": "rbm20_margin_resolve_rl125.tsv",
        "F5_REC_SCAN": "SE=PTBP1/0.062/-0.025/0.122,"
                       "MXE=PTBP1/0.113/0.016/0.187",
        "F5_FLOOR_RANK": "SE=6,MXE=3",
        "F5_REC_FWD": "SE=-0.051/-0.159/0.047,MXE=-0.019/-0.108/0.072"})
    assert q2.returncode == 0, q2.stderr
    assert "_rl125" in q2.stdout, q2.stdout
    assert "('MXE', 'PTBP1')" in q2.stdout and "('SE', 'PTBP1')" in q2.stdout, \
        q2.stdout
    assert "RANK {'SE': 6, 'MXE': 3}" in q2.stdout, q2.stdout
    assert "FWDREC (-0.051, -0.159, 0.047)" in q2.stdout, q2.stdout

    # --- 90's guard: suffixed file, no records
    q3 = run(n90, "n90.py", {"F5_SCAN": "signature_rbp_margins_rl125.tsv"})
    assert q3.returncode != 0, q3.stdout
    assert "wrong gene rather than fail" in q3.stderr + q3.stdout, q3.stderr

    # --- all-or-nothing: a bad 90 anchor must not let 65 through
    _n, pbad = apply_edits("nothing here\n", EDITS_90, "def _f5_env(")
    assert pbad and len(pbad) == 2, pbad

    # --- double application refused for both
    _x, pd65 = apply_edits(n65, EDITS_65, "def _s7_env(")
    assert pd65 and "already carries 155" in pd65[0], pd65
    _y, pd90 = apply_edits(n90, EDITS_90, "def _f5_env(")
    assert pd90 and "already carries 155" in pd90[0], pd90

    print("[ok] self-test passed: unsuffixed output is byte-identical to "
          "pre-patch for BOTH scripts (QKI/RBM20 and ranks 5/2 preserved); "
          "supplying the records gives PTBP1/PTBP1, ranks 6/3, 131/123/119 and "
          "the 125 bp forward triple; a suffixed input with NO records REFUSES "
          "in both scripts, naming the wrong-gene risk; a PARTIAL override is "
          "also refused and names what is missing; an unsuffixed input without "
          "records does not refuse; A_TESTABLE and CANDIDATES are untouched; "
          "a bad anchor in one file blocks both; double application refused")
    return 0


if __name__ == "__main__":
    sys.exit(main())
