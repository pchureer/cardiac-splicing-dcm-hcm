#!/usr/bin/env python3
"""
172_verify_regeneration.py -- re-run every producer and prove the submitted
figures and tables are what the deposited code makes.

⚠ THE RISK THIS SCRIPT EXISTS TO MANAGE
---------------------------------------
Re-running producers days before submission can overwrite a good figure with a
broken one. Five of the six figure producers in this project passed all of their
own internal gates while drawing a stale constant, and several now REFUSE to run
unless an override environment variable is set. So:

  * every artefact is copied to a backup directory BEFORE anything runs;
  * a producer that exits non-zero, or writes nothing, has its original
    RESTORED immediately;
  * nothing is ever deleted;
  * --restore puts a whole backup directory back in one command.

Three phases
------------
  1. --discover   find which script writes each artefact, by grepping the
                  scripts for the artefact's own filename. Writes a spec TSV.
                  It does NOT guess a command it cannot see; those are left
                  blank for you to fill, and phase 3 skips them loudly.
  2. (you edit the spec: fill `command` and `env` for each row)
  3. --verify     back up, run, compare, report.

Comparison
----------
  .png          decoded pixels, not bytes. Matplotlib writes a timestamp into
                PNG metadata, so byte-identity is the wrong test and would
                report every figure as changed. Reported as:
                  IDENTICAL       bytes equal
                  PIXEL-IDENTICAL bytes differ, every pixel equal  <- expected
                  DIFFERS         pixels differ (with how many, and by how much)
  .tsv .csv     SHA-256, then a line diff if they differ
  .md .txt      SHA-256, then a unified diff

Usage
-----
  python 172_verify_regeneration.py --discover \
      --root /data/project/chureeratp/RNA/HEART \
      --artefacts F1_discovery_rl125.png,F2_...png,Table1.tsv \
      --spec verify_spec.tsv

  # or let it find every _rl125 figure and every table named in the spec seed
  python 172_verify_regeneration.py --discover --root ... --auto --spec verify_spec.tsv

  # review verify_spec.tsv, fill in `command` and `env`, then
  python 172_verify_regeneration.py --verify --root ... --spec verify_spec.tsv
  python 172_verify_regeneration.py --verify --root ... --spec verify_spec.tsv --run

  # if anything went wrong
  python 172_verify_regeneration.py --restore verify_backup_20261024_143000
"""

import argparse
import difflib
import hashlib
import os
import re
import shlex
import shutil
import subprocess
import sys
import time
from pathlib import Path

# A pipe closed early (`| head`) otherwise ends the run in a BrokenPipeError
# traceback that reads like a failure in the work rather than in the output.
try:
    import signal as _signal
    _signal.signal(_signal.SIGPIPE, _signal.SIG_DFL)
except (AttributeError, ValueError, ImportError):
    pass

SPEC_COLS = ["artefact", "producer", "patched_by", "command", "env",
             "env_optional", "outputs", "notes"]

# Driver scripts that invoke producers. Their command lines are the canonical
# invocations, so they are harvested to pre-fill `command`.
DRIVER_HINTS = ("make_all", "run_all", "_driver", "regenerate", "rebuild")

# Overrides this project's producers now demand. Listed so --discover can tell
# you which one a given script will refuse without.
KNOWN_OVERRIDES = {
    "F1_SLOPES": "Figure 1's cohort slopes (43_, patched by 169_)",
    "F6_COUNTS": "Figure 6's g:Profiler counts (41_, patched by 165_)",
    "F6_TOTALS": "Figure 6's enrichment totals (41_)",
}

SKIP_DIRS = {".git", "deposit", "node_modules", "__pycache__", "probe_rl125_tmp",
             "..."}

# ⚠ The verification and deposit tooling must not be read as producers.
# 173_fill_spec.py holds every artefact's filename as a literal in its KNOWN
# table AND writes files, so it classifies as a writer of all of them and
# displaces the real producers. 172 skipped only itself; one script over was
# enough to break it. Any tool that TABULATES artefacts belongs here.
TOOLING = {"171_build_deposit.py", "172_verify_regeneration.py",
           "173_fill_spec.py", "174_diff_images.py",
           "175_measure_label_ink.py"}


def is_tooling(path):
    """The explicit set, plus any script that names the harness itself.

    174 landed in the project directory and became F4's producer: its usage
    example holds `F4_replication_rl125.png` and it writes files. That is the
    FOURTH tool to do this. Any script that refers to verify_spec.tsv or to
    172 by name is part of the harness, whatever its number.
    """
    if path.name in TOOLING:
        return True
    # Read ONLY the head. read_text() pulls the whole file into memory before
    # slicing, and this was being called on multi-GB tables.
    try:
        with open(path, "rb") as fh:
            head = fh.read(6000).decode("utf-8", "replace")
    except OSError:
        return False
    return ("verify_spec" in head or "172_verify_regeneration" in head
            or "verify_backup_" in head)
# Directory name prefixes to skip. A previous run's backup holds copies of every
# artefact under its ORIGINAL name, so walking it finds each figure twice and
# the second copy is not a deliverable.
SKIP_PREFIXES = ("verify_backup_",)


def keep_dir(name):
    return name not in SKIP_DIRS and not name.startswith(SKIP_PREFIXES)


def say(*a):
    print(*a, flush=True)


def sha256(path, chunk=1 << 20):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while True:
            b = fh.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


# ------------------------------------------------------------------ discovery

def all_scripts(root, extra_skip=()):
    skip = set(TOOLING) | {Path(__file__).name} | set(extra_skip)
    out = []
    for dirpath, dirnames, filenames in os.walk(root):
        dirnames[:] = [d for d in dirnames if keep_dir(d)]
        for fn in filenames:
            # ⚠ Cheap tests first. The content test was running on EVERY file
            # under --root -- thousands of per-sample VCFs and the 1 GB rMATS
            # tables -- and the discover step hung for minutes reading them to
            # inspect their first 6000 characters.
            if not (fn.endswith((".py", ".sh", ".sbatch", ".R"))
                    and re.match(r"^\d", fn)):
                continue
            if fn in skip or is_tooling(Path(dirpath) / fn):
                continue      # tooling: tabulates artefacts, does not make them
            out.append(Path(dirpath) / fn)
    return sorted(out, key=lambda p: p.name)



def clean_command(line):
    """Turn a harvested line into a runnable command, or return "".

    A numbered list inside a docstring looks exactly like an invocation:
        1. python 87_build_table3.py            # G2 should now pass
    That exact line was pre-filled into the spec for Table 3. A command must
    START with an interpreter or a script, after list markers and prompts are
    stripped, or it is prose and is dropped.
    """
    line = line.strip()
    line = re.sub(r"^[-*>]\s+", "", line)                 # bullets, quotes
    line = re.sub(r"^\(?\d+[.)]\s+", "", line)            # 1.  1)  (1)
    line = re.sub(r"^\$\s*", "", line)                    # shell prompt
    line = re.sub(r"^\s*(srun|sbatch)\s+(-\S+\s+|--\S+\s+)*", "", line)
    line = re.sub(r"\s+#.*$", "", line)                   # trailing comment
    line = re.sub(r"\s*[;&|]+\s*$", "", line).strip()
    if not re.match(r"^(python\d?|bash|sh|Rscript)\s+\S|^\./?\d+[a-z]?_|"
                    r"^\d+[a-z]?_[\w.-]+\.(py|sh|sbatch|R)\b", line):
        return ""
    return line


def harvest_invocations(producer, texts):
    """Command lines in OTHER scripts that invoke `producer`.

    A driver such as 89_make_all_figures.sh holds the canonical invocation,
    flags and all. Reading it beats reconstructing the flags from argparse.
    """
    found = []
    pat = re.compile(r"^[^\n#]*\b" + re.escape(producer) + r"\b[^\n]*$", re.M)
    for s, txt in texts.items():
        if s.name == producer:
            continue
        for m in pat.finditer(txt):
            line = m.group(0).strip()
            if not re.search(r"(python|bash|sh|sbatch|srun|Rscript)\b", line):
                continue
            line = clean_command(line)
            if not line:
                continue
            is_driver = any(h in s.name for h in DRIVER_HINTS) or s.suffix in (".sh",)
            found.append((line, s.name, is_driver))
    # prefer driver scripts, then dedupe on the command text
    found.sort(key=lambda t: (not t[2], t[1]))
    seen, keep = set(), []
    for line, who, drv in found:
        if line in seen:
            continue
        seen.add(line)
        keep.append((line, who, drv))
    return keep



# Environment variables every process has; reading one is not an override.
ENV_NOISE = {
    "PATH", "HOME", "USER", "LOGNAME", "PWD", "OLDPWD", "SHELL", "TERM",
    "TMPDIR", "TMP", "TEMP", "LANG", "LC_ALL", "LANGUAGE", "PYTHONPATH",
    "PYTHONHASHSEED", "MPLBACKEND", "MPLCONFIGDIR", "CONDA_PREFIX",
    "CONDA_DEFAULT_ENV", "VIRTUAL_ENV", "OMP_NUM_THREADS",
    "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS", "OPENBLAS_NUM_THREADS",
    "SLURM_JOB_ID", "SLURM_ARRAY_TASK_ID", "SLURM_ARRAY_JOB_ID",
    "SLURM_CPUS_PER_TASK", "SLURM_MEM_PER_NODE", "SLURM_NTASKS",
    "HOSTNAME", "DISPLAY", "COLUMNS", "NO_COLOR",
}


def env_prefix(txt):
    """A prefix the helper prepends to the name it is given, or "".

    90's helper is called `_f5_env("REC_SCAN", ...)` but reads
    `F5_REC_SCAN` from the environment. Reporting the INNER name puts a
    variable in the spec that the producer never reads, and the run then
    refuses for a reason that looks like the value being wrong.
    """
    for pat in (r"""environ(?:\.get)?\(?\s*\[?\s*f["']([A-Z][A-Z0-9_]*_)\{""",
                r"""environ(?:\.get)?\(?\s*\[?\s*["']([A-Z][A-Z0-9_]*_)["']\s*\+""",
                r"""getenv\(\s*f["']([A-Z][A-Z0-9_]*_)\{""",
                r"""getenv\(\s*["']([A-Z][A-Z0-9_]*_)["']\s*\+"""):
        m = re.search(pat, txt)
        if m:
            return m.group(1)
    return ""


def env_overrides(txt):
    """Every environment variable the producer READS, and whether it REFUSES.

    A hard-coded list was wrong: it knew F1_SLOPES, F6_COUNTS and F6_TOTALS and
    missed F3_SUFFIX, F3_HCM_BASE, F3_PANEL_A, F6_MULTI and T3_RECORD. Reading
    them out of the producer cannot go stale when the next patch adds one.

    "Refuses" is decided by FOLLOWING THE LOCAL VARIABLE to a guard, not by
    proximity to a `sys.exit`. A window of a few lines marked F3_DEBUG as
    mandatory because three other variables' guards happened to sit under it --
    and a variable wrongly called mandatory makes the row refuse for no reason.
    """
    # The last pattern catches a HELPER: 90_build_f5.py reads its seven
    # overrides through `_f5_env("F5_SCAN", ...)`, so the names never appear
    # inside an os.environ call and a direct-only scan returned nothing for
    # Figure 5 -- the figure with the most dangerous override of all, because
    # panel D is plotted from constants no gate checks.
    pats = (r"""os\.environ\.get\(\s*["'](\w+)["']""",
            r"""os\.environ\[\s*["'](\w+)["']""",
            r"""os\.getenv\(\s*["'](\w+)["']""",
            r"""\benviron\.get\(\s*["'](\w+)["']""",
            r"""\b\w*(?:env|getenv|override)\w*\(\s*["']"""
            r"""([A-Z][A-Z0-9_]{2,})["']""")
    lines = txt.splitlines()
    found, local = {}, {}            # var -> hard ; var -> {local names}
    for line in lines:
        for pat in pats:
            for m in re.finditer(pat, line):
                v = m.group(1)
                if v in ENV_NOISE:
                    continue
                found.setdefault(v, False)
                # A helper hides the guard: 90 calls
                # `_f5_env("F5_TECH45", required=True)` and the sys.exit lives
                # inside _f5_env, where following the local variable cannot
                # reach it. The CALL says whether it is mandatory.
                call = line[m.start():]
                cut = call.find(")")
                call = call[:cut] if cut > 0 else call
                if re.search(r"(required|mandatory|must|strict)\s*=\s*True",
                             call):
                    found[v] = True
                elif re.search(r"""^\s*\w*(?:env|getenv|override)\w*\(\s*"""
                               r"""["']""" + re.escape(v) + r"""["']\s*$""",
                               call):
                    # name alone, no default supplied -> nothing to fall back to
                    found[v] = True
                asg = re.match(r"\s*([A-Za-z_]\w*)\s*=", line)
                if asg:
                    local.setdefault(v, set()).add(asg.group(1))
                # a guard written inline: if not os.environ.get("V"): exit
                if re.search(r"\b(if|assert)\b", line) and \
                        re.search(r"(sys\.exit|SystemExit|raise)", line):
                    found[v] = True

    EXIT = re.compile(r"(sys\.exit|SystemExit|raise\s+SystemExit|"
                      r"parser\.error|ap\.error)")
    for v, names in local.items():
        for nm in names:
            word = re.compile(r"\b" + re.escape(nm) + r"\b")
            for i, line in enumerate(lines):
                if not word.search(line):
                    continue
                if not re.match(r"\s*(if|elif|assert|while)\b", line):
                    continue
                # the guard's own line, or the block it governs
                if EXIT.search(line):
                    found[v] = True
                    break
                indent = len(line) - len(line.lstrip())
                for nxt in lines[i + 1:i + 5]:
                    if not nxt.strip():
                        continue
                    if len(nxt) - len(nxt.lstrip()) <= indent:
                        break
                    if EXIT.search(nxt):
                        found[v] = True
                        break
                if found[v]:
                    break
    # A third idiom: a LIST OF NAMES tested for presence, then one exit.
    #   _miss = [_n for _n in ("F5_REC_SCAN", ..., "F5_D_LOG")
    #            if not os.environ.get(_n)]
    #   if _miss: sys.exit("REFUSING: ...")
    # F5_D_LOG is mandatory this way AND has a fallback at its own read site,
    # so following the local variable classified it as optional. It is not.
    for m in re.finditer(r"""[\(\[]((?:\s*["'][A-Z][A-Z0-9_]{2,}["']\s*,?\s*)+)"""
                         r"""[\)\]]""", txt):
        tail = txt[m.end():m.end() + 400]
        head = txt[max(0, m.start() - 200):m.start()]
        if not re.search(r"(environ|getenv)", head + tail):
            continue
        if not re.search(r"(sys\.exit|SystemExit|REFUSING|\braise\b)", tail):
            continue
        for q in re.finditer(r"""["']([A-Z][A-Z0-9_]{2,})["']""", m.group(1)):
            v = q.group(1)
            if v not in ENV_NOISE:
                found[v] = True

    pre = env_prefix(txt)
    if pre:
        found = {(k if k.startswith(pre) else pre + k): v
                 for k, v in found.items()}
    return sorted(found.items())


def argparse_flags(txt):
    """Flags with their defaults -- a flag alone does not tell you what to pass.

    `--out=F2_discovery.png` says the 50 bp name is the default, so the 125 bp
    run must override it. That is the difference between a command you can write
    and one you have to go and read the script for.
    """
    out = []
    for m in re.finditer(r"""add_argument\(\s*["'](--[\w-]+)["']([^;]*?)\)\s*$""",
                         txt, re.M):
        flag, rest = m.group(1), m.group(2)
        lab = flag
        if re.search(r"required\s*=\s*True", rest):
            lab += " (REQUIRED)"
        else:
            act = re.search(r"""action\s*=\s*["'](store_true|store_false|append|count)["']""",
                            rest)
            d = re.search(r"""default\s*=\s*(["'][^"']*["']|[\w.+-]+)""", rest)
            if act:
                lab += f" [{act.group(1)}]"
            elif d:
                lab += f"={d.group(1)}"
        out.append(lab)
    return sorted(set(out))


def fragment_search(stem, texts):
    """For a filename assembled at runtime: find the longest literal prefix,
    and the longest literal suffix, that some script actually contains."""
    hits = []
    parts = stem.split("_")
    for cut in range(len(parts), 0, -1):
        frag = "_".join(parts[:cut])
        if len(frag) < 6:
            break
        who = [s.name for s, t in texts.items() if frag in t]
        if who:
            hits.append(("prefix", frag, who))
            break
    for cut in range(0, len(parts)):
        frag = "_".join(parts[cut:])
        if len(frag) < 6:
            break
        who = [s.name for s, t in texts.items() if frag in t]
        if who:
            hits.append(("suffix", frag, who))
            break
    return hits


def find_artefacts(root, patterns):
    found = []
    for dirpath, dirnames, filenames in os.walk(root):
        dirnames[:] = [d for d in dirnames if keep_dir(d)]
        for fn in filenames:
            p = Path(dirpath) / fn
            for pat in patterns:
                if re.search(pat, fn):
                    found.append(p)
                    break
    return sorted(set(found))


def patch_targets(patch_name, patch_text, script_names):
    """Which script does this patch edit?

    This project's patch filenames name their own targets:
      169_patch_43_figure1_slopes.py        -> 43
      162_patch_figures_37_38_readlength.py -> 37, 38
      128_patch_61_73_table2.py             -> 61, 73
      165_patch_41_figure6_rl125.py         -> 41   (not 6, not 125)
    Tokens are split on '_' and only PURELY NUMERIC tokens are taken, so
    'figure6' and 'rl125' do not contribute. The body is also read for an
    explicit target filename, which wins when present.
    """
    # explicit, from the body
    explicit = sorted(set(
        m.group(1) for m in
        re.finditer(r"""["'](\d+[a-z]?_[\w.-]+\.(?:py|sh|sbatch|R))["']""",
                    patch_text)
        if m.group(1) != patch_name))
    if explicit:
        return [e for e in explicit if e in script_names], "its body"

    rest = re.sub(r"^\d+[a-z]?_patch_", "", patch_name)
    rest = re.sub(r"\.(py|sh|sbatch|R)$", "", rest)
    nums = [t for t in rest.split("_") if t.isdigit()]
    hits = []
    for n in nums:
        for nm in script_names:
            if re.match(r"^0*" + n + r"[a-z]?_", nm) and "_patch_" not in nm:
                hits.append(nm)
    return sorted(set(hits)), "its filename"


def classify_mentions(token, texts):
    """Which scripts WRITE `token`, which merely read it, which patch it.

    A grep for a filename finds readers, writers, patches and drivers alike.
    Five of six figures in this project were attributed to a patch for exactly
    that reason. So each mentioning script is classified by WHAT IT DOES on the
    lines that name the file.
    """
    W_PY = re.compile(
        r"(savefig|to_csv|to_excel|to_markdown|to_json|to_parquet|savetxt|"
        r"write_text|write_bytes|\.write\(|\.writelines\(|csv\.writer|"
        r"ExcelWriter|open\s*\([^)]*[\"']\s*[wax])")
    W_SH = re.compile(r"(>>?\s*\S|\btee\b|\s-o\s|--out(put)?[=\s])")
    R_ANY = re.compile(
        r"(read_csv|read_table|read_excel|read_json|read_parquet|loadtxt|"
        r"genfromtxt|imread|read_text|\bcat\b|\bawk\b|\bgrep\b|\bhead\b|"
        r"\bwc\b|\bsort\b|open\s*\([^)]*[\"']\s*r|<\s*\S)")
    INVOKE = re.compile(r"(python\d?|bash|sh|Rscript|sbatch|srun)\s+"
                        r"[\w./-]*\d+[a-z]?_[\w.-]+")

    writers, readers, patches, callers, ambiguous = [], [], [], [], []
    evidence = {}
    for s, txt in texts.items():
        if token not in txt:
            continue
        lines = [l for l in txt.splitlines()
                 if token in l and not l.lstrip().startswith("#")]
        if not lines:
            continue
        if "_patch_" in s.name:
            patches.append(s.name)
            continue
        if lines and all(INVOKE.search(l) for l in lines):
            callers.append(s.name)
            continue

        is_sh = s.suffix in (".sh", ".sbatch")
        wrote = any(W_PY.search(l) or (is_sh and W_SH.search(l)) for l in lines)
        if not wrote:
            # the name may be assigned to a variable that is written later
            for l in lines:
                m = re.match(r"\s*([A-Za-z_]\w{1,})\s*=\s*([^=].*)$", l)
                # the RHS must be a STRING EXPRESSION holding the token, or
                # this matches `df = pd.read_csv("Table1.tsv")` and calls the
                # reader a writer.
                if not m or token not in m.group(2):
                    continue
                if not ('"' in m.group(2) or "'" in m.group(2)):
                    continue
                if re.search(r"read_|loadtxt|genfromtxt|imread", m.group(2)):
                    continue
                v = re.escape(m.group(1))
                if re.search(r"(savefig|to_csv|to_markdown|savetxt|open)\s*\(\s*"
                             + v + r"\b", txt) or \
                   re.search(r"\b" + v + r"\s*\)\s*\.write", txt) or \
                   re.search(r"open\s*\(\s*" + v + r"\s*,\s*[\"']\s*[wax]", txt):
                    wrote = True
                    break
        if wrote:
            ev = next((l.strip()[:120] for l in lines
                       if W_PY.search(l) or (is_sh and W_SH.search(l))),
                      lines[0].strip()[:120])
            writers.append(s.name)
            evidence[s.name] = ev
        elif any(R_ANY.search(l) for l in lines):
            readers.append(s.name)
        else:
            ambiguous.append(s.name)
    return (sorted(writers), sorted(readers), sorted(patches),
            sorted(callers), sorted(ambiguous), evidence)


def discover(root, artefacts, scripts):
    texts = {}
    for s in scripts:
        try:
            texts[s] = s.read_text(encoding="utf-8", errors="replace")
        except OSError:
            continue

    rows = []
    for a in artefacts:
        base, stem = a.name, a.stem
        producer, command, notes = "", "", []

        # Try the full basename first. If no script holds it, retry with the
        # STEM: a producer that does savefig(OUT + ".png") or takes --out
        # without the extension never contains the full filename.
        token = base
        w, r, pat, call, amb, ev = classify_mentions(base, texts)
        if not (w or r or pat or call or amb):
            token = stem
            w, r, pat, call, amb, ev = classify_mentions(stem, texts)
            if w or r or pat or call or amb:
                notes.append(f"matched on the stem '{stem}', not the full "
                             f"filename -- the extension is appended by the "
                             f"producer or passed on the command line")

        if len(w) == 1:
            producer = w[0]
            notes.append("writes it directly")
        elif len(w) > 1:
            producer = "|".join(w)
            notes.append(f"{len(w)} scripts look like writers -- pick the "
                         f"one that makes the submitted copy. The line that "
                         f"made each look like a writer: " +
                         " /// ".join(f"[{k}] {ev.get(k,'?')}" for k in w))
        elif amb:
            producer = "|".join(amb)
            notes.append(f"names it but with no recognisable read or write "
                         f"idiom on that line; confirm by opening it")
        elif call:
            inv_lines, invoked = [], []
            for s, txt in texts.items():
                if s.name not in call:
                    continue
                for l in txt.splitlines():
                    if token not in l or l.lstrip().startswith("#"):
                        continue
                    m = re.search(r"(?:python\d?|bash|sh|Rscript|sbatch|srun)\s+"
                                  r"(?:-\S+\s+)*([\w./-]*?(\d+[a-z]?_[\w.-]+))", l)
                    if not m:
                        continue
                    cmd = re.sub(r"^\s*(srun|sbatch)\s+(-\S+\s+|--\S+\s+)*", "",
                                 l.strip())
                    inv_lines.append((re.sub(r"\s*[;&|]+\s*$", "", cmd), s.name))
                    invoked.append(m.group(2))
            producer = "|".join(sorted(set(invoked))) or "|".join(call)
            if len(set(c for c, _ in inv_lines)) == 1:
                command = inv_lines[0][0]
                notes.append(f"the producer never holds this name -- it is "
                             f"passed in. Command read from {inv_lines[0][1]}.")
            elif inv_lines:
                command = inv_lines[0][0]
                notes.append(f"{len(set(c for c,_ in inv_lines))} invocations; "
                             f"pre-filled the first. Others: " +
                             " /// ".join(sorted(set(c for c, _ in inv_lines))[1:4]))
        elif pat:
            # Only patches hold the name. That is not a dead end: the patch
            # names its target, and after patching the producer builds the
            # name at runtime from a read-length argument, which is why no
            # literal survives in it. Resolve the target and confirm it by
            # looking for the longest literal prefix of the stem inside it.
            names = {k.name for k in texts}
            tgt, how = [], ""
            for pn in pat:
                t = next((v for k, v in texts.items() if k.name == pn), "")
                got, how = patch_targets(pn, t, names)
                tgt.extend(got)
            tgt = sorted(set(tgt))
            if tgt:
                producer = "|".join(tgt)
                notes.append(f"no non-patch script holds this name; the "
                             f"producer builds it at runtime. Target resolved "
                             f"from the patch ({how}): {', '.join(tgt)}.")
                # confirmation: does the target hold a prefix of the stem?
                conf, confirmed = [], []
                for t_name in tgt:
                    txt = next((v for k, v in texts.items()
                                if k.name == t_name), "")
                    parts, got = stem.split("_"), None
                    # a `break` on the too-short fragment would skip a for/else,
                    # so the outcome is tracked explicitly
                    for cut in range(len(parts), 0, -1):
                        frag = "_".join(parts[:cut])
                        if len(frag) < 4:
                            break
                        if frag in txt:
                            got = frag
                            break
                    if got:
                        confirmed.append((len(got), t_name, got))
                        conf.append(f"confirmed: {t_name} holds '{got}'")
                    else:
                        conf.append(
                            f"{t_name} could not be confirmed: '{stem}' has no "
                            f"searchable literal fragment of 4+ characters"
                            if len(stem.split("_")[0]) < 4 else
                            f"{t_name} holds no literal fragment of '{stem}'")
                notes.append("; ".join(conf))
                # A patch can edit several producers at once -- 162 patches both
                # 37 and 38, 128 patches both 61 and 73. Narrow by the LONGEST
                # confirmed fragment: 61 holds 'Table2_rl125' and 73 only
                # 'Table2', so 61 writes this file and 73 writes the 50 bp one.
                # Uniqueness alone is not enough -- both confirmed.
                if len(tgt) > 1 and confirmed:
                    confirmed.sort(reverse=True)
                    best = [c for c in confirmed if c[0] == confirmed[0][0]]
                    if len(best) == 1:
                        producer = best[0][1]
                        notes.append(
                            f"narrowed to {producer}: of the {len(tgt)} scripts "
                            f"this patch edits, it holds the longest literal "
                            f"fragment of this artefact's name "
                            f"('{best[0][2]}', {best[0][0]} chars)"
                            + (f" against '{confirmed[1][2]}' in "
                               f"{confirmed[1][1]}" if len(confirmed) > 1 else ""))
                    else:
                        notes.append(
                            f"⚠ {len(best)} targets tie on the longest fragment "
                            f"'{best[0][2]}' -- open them and find which writes "
                            f"it")
                elif len(tgt) > 1:
                    notes.append(f"⚠ none of the {len(tgt)} targets holds any "
                                 f"literal fragment of the name -- open them "
                                 f"and find which one writes it")
            else:
                notes.append(f"ONLY patch scripts name this ({', '.join(pat)}) "
                             f"and the target could not be resolved from the "
                             f"patch name or body. Open one and read its "
                             f"target. Never put a patch in `command` -- it is "
                             f"already applied.")
        elif r:
            notes.append(f"only READ by {', '.join(r)}; nothing writes it. "
                         f"Either the producer is gone or the file was made by "
                         f"hand. An artefact nobody can rebuild cannot be "
                         f"verified.")
        else:
            notes.append("NO SCRIPT names this file or its stem. Find out what "
                         "made it before it is submitted.")

        if r and producer:
            notes.append(f"read-only: {', '.join(r)}")
        if pat and producer:
            pass        # recorded in the patched_by column

        cands = [c for c in producer.split("|") if c]
        if not command and len(cands) == 1:
            inv = harvest_invocations(cands[0], texts)
            if inv:
                command = inv[0][0]
                notes.append(f"command read from {inv[0][1]}"
                             + (f" ({len(inv)} invocations found; check the "
                                f"others)" if len(inv) > 1 else ""))

        # The flags and their DEFAULTS are reported whether or not a command was
        # found, because a harvested command can be a 50 bp invocation and the
        # default is what tells you so.
        if len(cands) == 1:
            t = next((v for k, v in texts.items() if k.name == cands[0]), "")
            fl = argparse_flags(t)
            notes.append(("flags: " + " ".join(fl)) if fl
                         else "takes no arguments")
            missing = [f.split()[0] for f in fl if "(REQUIRED)" in f
                       and f.split()[0] not in command]
            if command and missing:
                notes.append(f"⚠ the harvested command omits REQUIRED "
                             f"{', '.join(missing)} -- it will refuse. Complete "
                             f"it before running.")
            elif not command:
                notes.append("no script invokes it, so no command to copy -- "
                             "build one from the flags above")

        need, soft = [], []
        for cand in cands:
            t = next((v for k, v in texts.items() if k.name == cand), "")
            for var, hard in env_overrides(t):
                if hard and var not in need:
                    need.append(var)
                elif not hard and var not in soft:
                    soft.append(var)
        soft = [v for v in soft if v not in need]
        if soft:
            notes.append("also reads (no refusal seen, so a default applies — "
                         "check whether the default is the 50 bp one): "
                         + ", ".join(soft))
        if need:
            notes.append("REFUSES without: " + ", ".join(need))
            if len(cands) > 1:
                notes.append(f"⚠ `env` is the UNION over {len(cands)} candidate "
                             f"producers. Narrow `producer` first, then DELETE "
                             f"the variables that are not its own -- an override "
                             f"belonging to another script makes this row refuse "
                             f"for no reason.")

        rows.append({
            "artefact": str(a),
            "producer": producer,
            "patched_by": "|".join(pat),
            "command": command,
            "env": ";".join(f"{v}=" for v in need),
            # ⚠ An override WITH a default is more dangerous than one that
            # refuses: it falls back silently, and the fallback is whatever the
            # 50 bp run used. F4_CEILINGS defaults to 0.584/0.988 -- the 50 bp
            # ceilings -- and omitting it produced a 50 bp figure that passed
            # every check. These used to live only in `notes`, which is easy to
            # drop when reading the spec. They now have a column.
            "env_optional": ";".join(soft),
            "outputs": str(a),
            "notes": " | ".join(notes),
        })
    return rows



def write_spec(path, rows):
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\t".join(SPEC_COLS) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(c, "")).replace("\t", " ")
                               for c in SPEC_COLS) + "\n")


def read_spec(path):
    rows = []
    with open(path, encoding="utf-8") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        if hdr[:len(SPEC_COLS)] != SPEC_COLS:
            sys.exit(f"[fatal] {path}: header is {hdr}, expected {SPEC_COLS}. "
                     f"Keep the columns and their order.")
        for ln, line in enumerate(fh, 2):
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            parts += [""] * (len(SPEC_COLS) - len(parts))
            rows.append(dict(zip(SPEC_COLS, parts), _line=ln))
    return rows


# ------------------------------------------------------------------ compare

def load_pixels(path):
    try:
        from PIL import Image
        import numpy as np
    except ImportError:
        return None, ("Pillow and numpy are needed to compare figures by pixel. "
                      "`pip install pillow numpy` or `conda install pillow numpy`. "
                      "Without them a PNG can only be compared byte-wise, which "
                      "reports every matplotlib figure as changed.")
    try:
        with Image.open(path) as im:
            return np.asarray(im.convert("RGBA")), None
    except Exception as exc:
        return None, f"cannot decode {path}: {exc}"


def pdf_strip_dates(path):
    """PDF bytes with the timestamps removed.

    matplotlib writes /CreationDate into every PDF, so two PDFs from identical
    runs are never byte-identical and a hash comparison reports DIFFERS for
    every figure, every time. Stripping the dates makes the comparison about
    the drawing.
    """
    try:
        raw = Path(path).read_bytes()
    except OSError:
        return None
    for key in (b"CreationDate", b"ModDate"):
        out = bytearray()
        i = 0
        while True:
            j = raw.find(b"/" + key, i)
            if j < 0:
                out += raw[i:]
                break
            out += raw[i:j]
            k = raw.find(b")", j)
            i = (k + 1) if k > 0 else j + len(key) + 1
        raw = bytes(out)
    return raw


def compare_pdf(old_path, new_path):
    """Verdict for a pair of PDFs, ignoring embedded timestamps."""
    a, b = pdf_strip_dates(old_path), pdf_strip_dates(new_path)
    if a is None or b is None:
        return "DIFFERS (bytes)", "could not read one of the PDFs"
    if a == b:
        return "IDENTICAL", "byte-identical once /CreationDate is removed"
    # fall back to the text, then to a render if poppler is present
    try:
        ta = subprocess.run(["pdftotext", str(old_path), "-"],
                            capture_output=True, text=True, timeout=120)
        tb = subprocess.run(["pdftotext", str(new_path), "-"],
                            capture_output=True, text=True, timeout=120)
        if ta.returncode == 0 and tb.returncode == 0:
            if ta.stdout == tb.stdout:
                note = "text identical; the drawing differs"
            else:
                d = [l for l in difflib.unified_diff(
                    ta.stdout.splitlines(), tb.stdout.splitlines(),
                    "submitted", "regenerated", lineterm="", n=0)][:20]
                return "DIFFERS", ("TEXT CHANGED. " + "\n      "
                                   + "\n      ".join(d))
        else:
            note = "pdftotext unavailable"
    except (FileNotFoundError, subprocess.SubprocessError):
        note = "pdftotext unavailable"
    return "DIFFERS", (note + " — render both with `pdftoppm -r 300 -png "
                       "-singlefile` and compare with 174_diff_images.py")


def compare(old_path, new_path):
    """Returns (verdict, detail)."""
    if not Path(new_path).exists():
        return "NOT PRODUCED", "the command exited cleanly but wrote no file here"
    if sha256(old_path) == sha256(new_path):
        return "IDENTICAL", ""

    ext = Path(new_path).suffix.lower()
    if ext == ".pdf":
        return compare_pdf(old_path, new_path)
    if ext in {".png", ".jpg", ".jpeg", ".tif", ".tiff"}:
        a, err = load_pixels(old_path)
        if a is None:
            return "DIFFERS (bytes)", err
        b, err = load_pixels(new_path)
        if b is None:
            return "DIFFERS (bytes)", err
        if a.shape != b.shape:
            return "DIFFERS", f"size changed: {a.shape[1]}x{a.shape[0]} -> " \
                              f"{b.shape[1]}x{b.shape[0]}"
        import numpy as np
        d = np.abs(a.astype(int) - b.astype(int))
        nz = int((d.sum(axis=2) > 0).sum())
        if nz == 0:
            return "PIXEL-IDENTICAL", "bytes differ (PNG metadata only)"
        tot = a.shape[0] * a.shape[1]
        ys, xs = np.nonzero(d.sum(axis=2))
        return "DIFFERS", (f"{nz} of {tot} pixels ({100*nz/tot:.3f}%), max "
                           f"channel delta {int(d.max())}, bounding box "
                           f"x {xs.min()}-{xs.max()}, y {ys.min()}-{ys.max()}")

    # text / table
    try:
        o = Path(old_path).read_text(encoding="utf-8", errors="replace").splitlines()
        n = Path(new_path).read_text(encoding="utf-8", errors="replace").splitlines()
    except OSError as exc:
        return "DIFFERS (bytes)", str(exc)
    if len(o) != len(n):
        head = f"{len(o)} lines -> {len(n)} lines. "
    else:
        head = ""
    diff = [l for l in difflib.unified_diff(o, n, "submitted", "regenerated",
                                            lineterm="", n=0)][:24]
    return "DIFFERS", head + "\n      " + "\n      ".join(diff)


# ------------------------------------------------------------------ verify

def verify(root, rows, do_run, backup_dir):
    results = []
    runnable = [r for r in rows if r["command"].strip()]
    skipped = [r for r in rows if not r["command"].strip()]

    say(f"\n{len(runnable)} row(s) with a command, {len(skipped)} without.")
    if skipped:
        say("\n⚠ NO COMMAND -- these artefacts will NOT be verified:")
        for r in skipped:
            say(f"    {Path(r['artefact']).name:38s} producer: "
                f"{r['producer'] or '(unknown)'}")
        say("  An unverified figure is not a passing figure. Fill `command` "
            "for each, or accept in writing that it was not checked.")

    # ⚠ a patch is not a regeneration
    bad_cmd = [r for r in runnable if "_patch_" in r["command"]]
    if bad_cmd:
        say("\n[fatal] `command` names a patch script. A patch edits a producer "
            "and has already been applied; running it again edits the producer "
            "a second time and verifies nothing:")
        for r in bad_cmd:
            say(f"    {Path(r['artefact']).name}: {r['command']}")
        say("  Put the PRODUCER in `command` and leave the patch in "
            "`patched_by`.")
        sys.exit(1)

    # env completeness
    bad_env = []
    for r in runnable:
        for kv in [x for x in r["env"].split(";") if x.strip()]:
            if "=" not in kv or not kv.split("=", 1)[1].strip():
                bad_env.append((r["artefact"], kv))
    if bad_env:
        say("\n[fatal] env values left empty -- the producer will refuse and "
            "the restore path will fire for no reason:")
        for a, kv in bad_env:
            say(f"    {Path(a).name}: {kv}")
        sys.exit(1)

    if not do_run:
        say("\n--- would run (in order) ---")
        for r in runnable:
            # ⚠ Quote the values. F3_PANEL_A is
            #   MAGNet/DCM vs NF=813/53456,MAGNet/HCM vs NF=431/39607
            # which contains SPACES. The real run passes env through
            # subprocess's env dict and is unaffected, but the printed line was
            # not copy-pasteable: a shell would read `vs` as a command and
            # `NF=813/53456,...` as a second assignment, and the figure would
            # be drawn from the 50 bp fallback.
            parts = []
            for kv in r["env"].split(";"):
                kv = kv.strip()
                if not kv or "=" not in kv:
                    continue
                k, v = kv.split("=", 1)
                parts.append(f"{k}={shlex.quote(v)}" if v != shlex.quote(v)
                             else f"{k}={v}")
            e = " ".join(parts)
            say(f"  {' + '.join(p.name for p in map(Path, r['outputs'].split(','))) if r['outputs'] else Path(r['artefact']).name}")
            say(f"    {(e + ' ') if e else ''}{r['command']}")
        say("\nDRY RUN -- nothing backed up, nothing run. Add --run.")
        return []

    backup_dir.mkdir(parents=True, exist_ok=False)
    say(f"\nbackup: {backup_dir}")

    # ⚠ CHAIN OF CUSTODY. A backup is only "the submitted copy" if nothing
    # replaced the file since the last time it was backed up. F5 was restored
    # from one backup, replaced on disk by its regenerated version, and the
    # next run backed THAT up as the original and reported IDENTICAL against
    # itself. So each file about to be backed up is hashed against the same
    # path in every earlier verify_backup_* directory; a mismatch is reported
    # before the run, naming both hashes, so a self-comparison cannot pass as
    # a verification.
    earlier = sorted(d for d in root.glob("verify_backup_*")
                     if d.is_dir() and d != backup_dir)
    custody = []

    # back up EVERY output of EVERY runnable row, before anything executes
    backed = {}
    for r in runnable:
        for o in [x.strip() for x in r["outputs"].split(",") if x.strip()]:
            p = Path(o)
            if not p.is_absolute():
                p = root / p
            if not p.exists():
                say(f"  [warn] {p} does not exist yet; it will be treated as new")
                continue
            dest = backup_dir / p.relative_to(root) if root in p.parents \
                else backup_dir / p.name
            dest.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(p, dest)
            backed[str(p)] = dest
            rel = dest.relative_to(backup_dir)
            # Compare against the EARLIEST backup that holds this path, and
            # only that one -- it is the submitted copy. The first version
            # kept looking after an equal match and found a LATER backup that
            # differed, flagging F5 as replaced when it had just been restored.
            # And a PDF embeds /CreationDate, so its raw hash differs on every
            # run; content is compared with the dates stripped, or every PDF
            # reads as custody-broken forever.
            first_dir = next((d for d in earlier if (d / rel).is_file()), None)
            if first_dir is not None:
                first = first_dir / rel
                if p.suffix.lower() == ".pdf":
                    same = pdf_strip_dates(first) == pdf_strip_dates(p)
                else:
                    same = sha256(first) == sha256(p)
                if not same:
                    custody.append((rel, first_dir.name,
                                    sha256(first)[:12], sha256(p)[:12]))
    say(f"  backed up {len(backed)} artefact(s). Nothing was moved or deleted.")
    if custody:
        say(f"\n  ⚠⚠ {len(custody)} file(s) DIFFER from the copy an earlier run "
            f"backed up. The file on disk is NOT the same 'original' as "
            f"before -- it was replaced between runs, and a verdict against it "
            f"is a comparison with whatever replaced it, not with the "
            f"submitted artefact:")
        for rel, dname, h0, h1 in custody:
            say(f"      {rel}")
            say(f"          {dname}: {h0}…   now on disk: {h1}…")
        say("      The earliest backup is the submitted copy unless you know "
            "otherwise. Restore from it, or state which copy is authoritative, "
            "before trusting an IDENTICAL on these rows.")

    before = snapshot(root)
    unprotected = []

    for r in runnable:
        name = Path(r["artefact"]).name
        say(f"\n=== {name}")
        env = dict(os.environ)
        for col in ("env", "env_optional"):
            for kv in [x for x in (r.get(col) or "").split(";") if x.strip()]:
                if "=" not in kv:
                    continue
                k, v = kv.split("=", 1)
                if v.strip():
                    env[k] = v
        say(f"  $ {r['command']}")
        t0 = time.time()
        proc = subprocess.run(r["command"], shell=True, cwd=str(root), env=env,
                              capture_output=True, text=True)
        dt = time.time() - t0
        tail = (proc.stdout or "").strip().splitlines()[-6:]
        for l in tail:
            say(f"    | {l[:150]}")
        if proc.returncode != 0:
            err = (proc.stderr or "").strip().splitlines()[-6:]
            for l in err:
                say(f"    ! {l[:150]}")
            say(f"  [FAILED] exit {proc.returncode} after {dt:.1f}s -- RESTORING")
            for o in [x.strip() for x in r["outputs"].split(",") if x.strip()]:
                p = Path(o) if Path(o).is_absolute() else root / o
                if str(p) in backed:
                    shutil.copy2(backed[str(p)], p)
            results.append((name, "FAILED", f"exit {proc.returncode}"))
            continue

        for o in [x.strip() for x in r["outputs"].split(",") if x.strip()]:
            p = Path(o) if Path(o).is_absolute() else root / o
            old = backed.get(str(p))
            if old is None:
                results.append((Path(o).name, "NEW",
                                "no prior copy existed to compare against"))
                continue
            verdict, detail = compare(old, p)
            rel_o = Path(old).relative_to(backup_dir)
            if any(c[0] == rel_o for c in custody) and verdict in (
                    "IDENTICAL", "PIXEL-IDENTICAL"):
                verdict = "IDENTICAL (custody broken)"
                detail = ("matches the file that was on disk, but that file "
                          "differs from an earlier backup -- see the custody "
                          "warning above")
            results.append((Path(o).name, verdict, detail))
            say(f"  [{verdict}] {Path(o).name}")
            if detail:
                say(f"      {detail}")
            if verdict in ("DIFFERS", "DIFFERS (bytes)"):
                # ⚠ Keep the REGENERATED file too. --restore puts the submitted
                # copies back over it, and then there is nothing left to
                # compare -- a diff run after a restore compares each file with
                # itself and reports 0 pixels differ. Preserving it here means
                # the comparison survives the restore.
                keep = Path(str(old) + ".regenerated")
                try:
                    shutil.copy2(p, keep)
                except OSError as exc:
                    say(f"      [warn] could not preserve the regenerated "
                        f"copy: {exc}")
                    keep = None
                say(f"      submitted copy  : {old}")
                if keep:
                    say(f"      regenerated copy: {keep}")
                    say(f"      compare them with:")
                    say(f"        python 174_diff_images.py --old {old} "
                        f"--new {keep} --out {Path(o).stem}_diff")
            elif verdict == "NOT PRODUCED":
                say(f"      submitted copy kept at {old}")

    extra = unexpected_writes(before, root, set(backed))
    extra = [e for e in extra
             if not Path(e).name.startswith("verify_backup_")
             and ".bak." not in Path(e).name]
    if extra:
        say(f"\n⚠ {len(extra)} file(s) were written by these commands but are "
            f"NOT in any row's `outputs`, so NO BACKUP was taken and a "
            f"--restore will not put them back:")
        for e in extra[:40]:
            say(f"    {Path(e).relative_to(root) if root in Path(e).parents else e}")
        if len(extra) > 40:
            say(f"    ... and {len(extra) - 40} more")
        say("  Add them to the row's `outputs` and run again from a clean "
            "state, or accept that these are now the regenerated versions "
            "while the declared outputs may have been restored -- which leaves "
            "a MISMATCHED PAIR on disk.")

    return results


def snapshot(root):
    """mtime of every file under root, so a command's writes can be detected.

    The spec declares what a producer writes; the producer decides. 38 writes a
    .pdf beside every .png and the spec listed only the .png, so a verify run
    overwrote six PDFs with no backup -- and nothing reported it. A snapshot
    before and after each command turns that silence into a warning.
    """
    seen = {}
    for dirpath, dirnames, filenames in os.walk(root):
        dirnames[:] = [d for d in dirnames if keep_dir(d)]
        for fn in filenames:
            q = Path(dirpath) / fn
            try:
                seen[str(q)] = q.stat().st_mtime_ns
            except OSError:
                pass
    return seen


def unexpected_writes(before, root, backed):
    """Files whose mtime moved but which were never backed up."""
    out = []
    after = snapshot(root)
    for k, v in after.items():
        if k in backed:
            continue
        if before.get(k) != v:
            out.append(k)
    return sorted(out)


def restore(backup_dir, root):
    backup_dir = Path(backup_dir)
    if not backup_dir.is_dir():
        sys.exit(f"[fatal] not a directory: {backup_dir}")
    say("⚠ This replaces whatever is in place now with the submitted copies. "
        "If a row reported DIFFERS, its regenerated file is kept here as "
        "*.regenerated and is NOT restored, so it stays available to compare.")
    n = 0
    for dirpath, _, filenames in os.walk(backup_dir):
        for fn in filenames:
            if fn.endswith(".regenerated"):
                continue      # a preserved comparison copy, not an original
            src = Path(dirpath) / fn
            dest = root / src.relative_to(backup_dir)
            dest.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(src, dest)
            say(f"  restored {dest}")
            n += 1
    say(f"{n} file(s) restored from {backup_dir}. The backup was not deleted.")


# ------------------------------------------------------------------ main

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=".", help="project root")
    ap.add_argument("--spec", default="verify_spec.tsv")
    ap.add_argument("--discover", action="store_true")
    ap.add_argument("--verify", action="store_true")
    ap.add_argument("--run", action="store_true",
                    help="with --verify: actually back up and execute")
    ap.add_argument("--restore", default=None, metavar="BACKUP_DIR")
    ap.add_argument("--artefacts", default="",
                    help="comma-separated filenames or regexes to verify")
    ap.add_argument("--auto", action="store_true",
                    help="with --discover: find every *_rl125.png and every "
                         "Table*.tsv/.md under --root")
    ap.add_argument("--skip-script", action="append", default=[],
                    help="a script name to exclude from producer discovery, in "
                         "addition to the 171/172/173 tooling; repeatable. Use "
                         "for anything that tabulates artefact names rather "
                         "than writing them.")
    ap.add_argument("--only", default=None, metavar="FILE",
                    help="with --discover: a file of artefact names or paths, "
                         "one per line (# comments allowed). Restricts the set "
                         "to exactly the artefacts being submitted, so the "
                         "spec holds no 50 bp leftovers or *_BEFORE_* backups.")
    a = ap.parse_args()
    root = Path(a.root).resolve()

    # ⚠ a literal "..." copied out of a usage example is a real, creatable
    # directory name, and everything downstream then runs in the wrong place.
    if not root.is_dir():
        sys.exit(f"[fatal] --root is not a directory: {root}")
    if root.name == "...":
        sys.exit(f"[fatal] --root is literally '...': {root}\n"
                 f"        That was a placeholder in the usage example, not a "
                 f"path. Pass the real project root. If a '...' directory was "
                 f"created by an earlier run, remove it.")
    if not a.restore and not any(root.glob("[0-9]*_*.py")):
        sys.exit(f"[fatal] no numbered script (NN_*.py) directly under {root}. "
                 f"That is almost certainly the wrong --root -- the producers "
                 f"live beside the artefacts.")

    if a.restore:
        restore(a.restore, root)
        return

    if a.discover:
        pats = [re.escape(x.strip()) for x in a.artefacts.split(",") if x.strip()]
        if a.auto:
            pats += [r"_rl125\.png$", r"^Table\d.*\.(tsv|md)$",
                     r"^f\d.*_rl125\.png$"]
        if not pats:
            sys.exit("[fatal] --discover needs --artefacts or --auto.")
        arts = find_artefacts(root, pats)
        if a.only:
            only_p = Path(a.only)
            if not only_p.is_file():
                sys.exit(f"[fatal] --only file not found: {only_p.resolve()}\n"
                         f"        It is a plain text list, one artefact name "
                         f"per line, '#' for comments. Without it, drop --only "
                         f"and --auto finds everything, including the 50 bp "
                         f"copies and the *_BEFORE_* backups.")
            want = [l.strip() for l in only_p.read_text(encoding="utf-8").splitlines()
                    if l.strip() and not l.lstrip().startswith("#")]
            if not want:
                sys.exit(f"[fatal] {only_p} has no artefact lines, only "
                         f"comments or blanks.")
            keep, unmatched = [], []
            for w in want:
                # Anchored. A bare name means "directly under --root"; a
                # name with a slash means that exact relative path. A plain
                # endswith() let `t3_rl125/Table1.tsv` in on the line
                # `Table1.tsv` -- the very file this list excludes by name,
                # and the same same-basename-different-artefact trap as the
                # *_bams_dedup lists.
                wn = w.strip("/")
                if "/" in wn:
                    hit = [p for p in arts
                           if str(p.relative_to(root)) == wn]
                else:
                    hit = [p for p in arts
                           if p.name == wn and p.parent == root]
                if hit:
                    keep.extend(hit)
                elif "/" not in wn and any(p.name == wn for p in arts):
                    where = sorted({str(p.relative_to(root)) for p in arts
                                    if p.name == wn})
                    unmatched.append(
                        f"{w}  -- a bare name matches only files directly "
                        f"under --root; this one exists at: "
                        f"{', '.join(where)}")
                else:
                    unmatched.append(w)
            if unmatched:
                sys.exit("[fatal] --only lists artefacts that are not on disk "
                         "under --root:\n  " + "\n  ".join(unmatched) +
                         "\n        A name in the submitted set that does not "
                         "exist is the defect, not a filter problem. Two false "
                         "alarms in this project came from a guessed filename.")
            dropped = [p for p in arts if p not in keep]
            arts = sorted(set(keep))
            say(f"--only kept {len(arts)}, dropped {len(dropped)}")
        if not arts:
            sys.exit(f"[fatal] no artefact under {root} matched {pats}. "
                     f"Check the names against what is actually on disk -- "
                     f"guessing a filename is how two false alarms happened "
                     f"in this project already.")
        say(f"{len(arts)} artefact(s) found:")
        for p in arts:
            say(f"  {p.relative_to(root) if root in p.parents else p}")
        scripts = all_scripts(root, a.skip_script)
        say(f"\ngrepping {len(scripts)} script(s) for each filename ...")
        rows = discover(root, arts, scripts)
        write_spec(a.spec, rows)
        say(f"\nwrote {a.spec}")
        ready = [r for r in rows if r["command"]]
        one = [r for r in rows if r["producer"] and "|" not in r["producer"]]
        multi = [r for r in rows if "|" in r["producer"]]
        patch_only = [r for r in rows if not r["producer"] and r["patched_by"]]
        unknown = [r for r in rows if not r["producer"] and not r["patched_by"]]
        needs = [r for r in rows if r["env"]]

        say(f"  {len(one)} with a single non-patch producer")
        if ready:
            say(f"  {len(ready)} already have a `command`, read out of a driver "
                f"script -- check each one before running it")
        if multi:
            say(f"  {len(multi)} name more than one non-patch script -- pick "
                f"the writer:")
            for r in multi:
                say(f"      {Path(r['artefact']).name:30s} {r['producer']}")
        if patch_only:
            say(f"  ⚠ {len(patch_only)} are named ONLY by patch scripts. A patch "
                f"is not a producer -- it contains the filename because that is "
                f"the string it edits in. Open the patch and read its target:")
            for r in patch_only:
                say(f"      {Path(r['artefact']).name:30s} patched_by "
                    f"{r['patched_by']}")
        if unknown:
            say(f"  ⚠ {len(unknown)} with no producer at all (name assembled at "
                f"runtime, or written by hand):")
            for r in unknown:
                say(f"      {Path(r['artefact']).name}")
                if r["notes"]:
                    say(f"        {r['notes'][:180]}")
        if needs:
            say(f"  {len(needs)} will refuse without an override; `env` is "
                f"pre-filled with the variable names and empty values.")
        softs = [r for r in rows if r.get("env_optional")]
        if softs:
            say("")
            say(f"  ⚠⚠ {len(softs)} producer(s) read an override that HAS A "
                f"DEFAULT. These do NOT refuse -- they fall back silently, and "
                f"the fallback is whatever the 50 bp run used:")
            for r in softs:
                say(f"      {Path(r['artefact']).name:30s} "
                    f"{r['env_optional'].replace(';', ', ')}")
            say("      Decide for EACH whether its default is the 125 bp value. "
                "An override that refuses protects you; one that falls back "
                "does not.")
        say("\nNext: fill `command` and `env`, then --verify (dry) and "
            "--verify --run.")
        return

    if a.verify:
        if not Path(a.spec).exists():
            sys.exit(f"[fatal] no spec at {a.spec}. Run --discover first.")
        rows = read_spec(a.spec)
        if a.only:
            # --only restricts a verify/run to the listed artefacts. Until 26 Oct
            # it was honoured only by --discover, so a run meant to cover nine
            # supplement rows silently re-ran all 23 main ones instead.
            only_p = Path(a.only)
            if not only_p.is_file():
                sys.exit(f"[fatal] --only file not found: {only_p.resolve()}")
            want = {Path(l.strip()).name for l in only_p.read_text(encoding="utf-8").splitlines()
                    if l.strip() and not l.lstrip().startswith("#")}
            def names(r):
                # a merged row's artefact reads "<a> + <b>"; outputs lists every file
                out = {Path(x.strip()).name for x in r["outputs"].split(",") if x.strip()}
                out |= {Path(x.strip()).name for x in r["artefact"].split(" + ")}
                return out
            kept = [r for r in rows if names(r) & want]
            found = set().union(*(names(r) for r in kept)) & want if kept else set()
            missing = sorted(want - found)
            say(f"--only {only_p.name}: {len(kept)} of {len(rows)} spec row(s) match "
                f"{len(want)} listed artefact(s)")
            if missing:
                say("  not in the spec: " + ", ".join(missing))
            if not kept:
                sys.exit("[fatal] --only matched no spec row. Append the rows first "
                         "(tail -n +2 verify_spec_supplement.tsv >> verify_spec.tsv).")
            rows = kept
        stamp = time.strftime("%Y%m%d_%H%M%S")
        results = verify(root, rows, a.run, root / f"verify_backup_{stamp}")
        if not results:
            return
        say("\n" + "=" * 62)
        say("SUMMARY")
        order = ["IDENTICAL", "PIXEL-IDENTICAL", "IDENTICAL (custody broken)",
                 "NEW", "DIFFERS", "DIFFERS (bytes)", "NOT PRODUCED", "FAILED"]
        counts = {}
        for _, v, _ in results:
            counts[v] = counts.get(v, 0) + 1
        for v in order:
            if v in counts:
                say(f"  {counts[v]:3d}  {v}")
        for v in counts:
            if v not in order:
                say(f"  {counts[v]:3d}  {v}")
        skipped_rows = [r for r in rows if not r["command"].strip()]
        bad = [r for r in results
               if r[1] not in ("IDENTICAL", "PIXEL-IDENTICAL")]
        if skipped_rows:
            say(f"\n⚠ {len(skipped_rows)} artefact(s) had NO COMMAND and were "
                f"NOT run. They are not verified:")
            for r in skipped_rows:
                say(f"    {Path(r['artefact']).name}")
        if bad or skipped_rows:
            if not bad:
                say("\nEvery artefact that RAN reproduced -- but the set is "
                    "incomplete. Fill the missing command(s) and run again.")
                sys.exit(1)
            say("\nNot reproduced:")
            for n, v, d in bad:
                say(f"  [{v}] {n}")
            say(f"\n⚠ Compare the differing files BEFORE restoring -- a "
                f"restore overwrites the regenerated copy in place. The "
                f"regenerated versions are preserved in the backup as "
                f"*.regenerated, so the comparison survives either way.")
            say(f"\nEvery submitted copy is still in verify_backup_{stamp}. "
                f"To put them all back:")
            say(f"  python {Path(__file__).name} --restore "
                f"{root / ('verify_backup_' + stamp)} --root {root}")
            sys.exit(1)
        say("\nEvery artefact reproduced. The deposited code makes the "
            "submitted figures and tables.")
        return

    ap.print_help()


if __name__ == "__main__":
    main()
