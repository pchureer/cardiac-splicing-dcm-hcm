#!/usr/bin/env python3
"""
173_fill_spec.py -- fill verify_spec.tsv from the project record.

Every command below is quoted from the record of the run that actually produced
the submitted artefact, with the document named. Nothing here is reconstructed
from argparse defaults, because a default is what the 50 bp run used.

  python 173_fill_spec.py --spec verify_spec.tsv              # dry run
  python 173_fill_spec.py --spec verify_spec.tsv --apply
  python 173_fill_spec.py --spec verify_spec.tsv --apply --merge

--merge joins rows whose command is identical into one row with several
outputs, so a producer that writes both a .tsv and a .md is run once instead of
twice. Running it twice is not harmful, but the second run compares a file the
first run just rewrote.

It REFUSES to overwrite a non-empty field unless --force, so a value you
corrected by hand is not silently replaced.
"""

import argparse
import re
import shutil
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

# basename -> (producer, command, env, provenance)
# env is ';'-separated KEY=VALUE, exactly as 172 expects.
KNOWN = {
    "F1_design_rl125.png": (
        "43_figure1_design.py",
        "python3 43_figure1_design.py --table Table1.tsv --out F1_design_rl125",
        "F1_SLOPES=0.690/0.968",
        "⚠ SLASH-separated, not comma. 43 line 113 is "
        "`_F1_PLATFORM, _F1_PEOPLE = (float(_x) for _x in _F1S.split(\"/\"))`, "
        "so a comma raises ValueError: could not convert string to float: "
        "'0.690,0.968'. The record prints the pair as '0.690 and 0.968' and "
        "the separator was mine to get wrong. "
        "RESULT_f1_rl125_and_all_figures_complete.md -- panel B draws slope "
        "0.690 and 0.968. 43 prints '0.584, 0.988 (50 bp defaults)' when "
        "F1_SLOPES is unset, and 169's guard arms on 'rl125' appearing in argv, "
        "so --out F1_design_rl125 is what makes the guard fire.",
    ),
    "F2_discovery_rl125.png": (
        "39_figure2_discovery.py",
        "python3 39_figure2_discovery.py --event-type SE "
        "--master master_f2_panelB_rl125.tsv --discovery SE_rl125 "
        "--mats 'rmats_magnet_rl125/{ET}.MATS.JCEC.txt' --out F2_discovery_rl125",
        "",
        "RESULT_f2_rl125_regenerated.md, verbatim. 39 needed NO patch -- panels "
        "A and C were read-length aware already; only panel B was frozen, in a "
        "FILE. 164 must have written master_f2_panelB_rl125.tsv first.",
    ),
    "F3_etiology_rl125.png": (
        "37_figure3_etiology.py",
        "python3 37_figure3_etiology.py --event-type SE "
        "--dcm-vs-hcm SE_DCMvsHCM_rl125 --out F3_etiology_rl125",
        "F3_SUFFIX=_rl125;F3_HCM_BASE=SE_HCMvsNF_rl125;"
        "F3_PANEL_A=MAGNet/DCM vs NF=813/53456,MAGNet/HCM vs NF=431/39607",
        "RESULT_f3_f4_rl125_regenerated_and_legend_audit.md, verbatim. The run "
        "is confirmed by the log line '[info] direct DCM vs HCM: 0 of 38,273'. "
        "An earlier attempt printed 45,380 -- the 50 bp table -- under the same "
        "_rl125 output name. THE DENOMINATOR, NOT THE FILENAME, SAYS WHICH READ "
        "LENGTH A FIGURE IS BUILT FROM.",
    ),
    "F4_replication_rl125.png": (
        "38_figure4_replication.py",
        "python3 38_figure4_replication.py --event-type SE "
        "--conc direction_concordance_dedup_rl125.tsv --out F4_replication_rl125",
        "F4_CEILINGS=0.690/0.968;F4_SUFFIX=_rl125",
        "⚠ BOTH overrides are REQUIRED IN PRACTICE AND NEITHER REFUSES. "
        "38 line 59: _C4 = os.environ.get('F4_CEILINGS', '').strip(); with it "
        "unset, lines 63-64 fall back to CEIL_PEOPLE=0.988 and "
        "CEIL_PLATFORM=0.584 -- THE 50 bp CEILINGS. The 125 bp pair is "
        "0.690/0.968 (RESULT_s4_rl125.md), slash-separated like F1_SLOPES. "
        "F4_SUFFIX=_rl125 re-bases the MAGNet discovery arm only; line 69 says "
        "every REPL cohort must stay unsuffixed because none was re-quantified, "
        "so do NOT add a suffix there. "
        "--conc carries the 125 bp concordance table "
        "(DEFECT_33_overwrote_direction_concordance.md). "
        "Measured consequence of omitting F4_CEILINGS: the orange ceiling line "
        "moves 87 px left and the blue one 16 px right, both labels still four "
        "characters wide, so NOTHING in the text changes -- "
        "RESULT_f4_submitted_is_correct_my_command_drew_50bp_ceilings.md.",
    ),
    "F6_enrichment_rl125.png": (
        "41_figure6_enrichment.py",
        "python3 41_figure6_enrichment.py "
        "--correct compadj_rl125_gprofiler.tsv "
        "--genome compadj_rl125_gprofiler_WHOLEGENOME.tsv "
        "--out F6_enrichment_rl125",
        "F6_COUNTS=12564/805/550/255/73/131;F6_TOTALS=1/85;F6_MULTI=0/17",
        "DEFECT_f6_panelA_counts_every_query_list.md ('Expected after 166'), "
        "verbatim. Provenance must end 'terms: 1 vs 85' with NO *** warning. "
        "Panel A's bars are 1 and 85; panel C is 131 / 17 / 0.",
    ),
    "Table1.tsv": (
        "42_table1.py",
        "python3 42_table1.py --out Table1",
        "",
        "42's --out default IS 'Table1' -- Table 1 is the one artefact with no "
        "_rl125 suffix, because 149 patched 42 in place and kept the 50 bp copy "
        "as Table1_BEFORE_149.tsv. Passed explicitly anyway, so the spec says "
        "what it writes. F1_SLOPES was removed: it came from 43 being a "
        "candidate, and 43 READS Table 1.",
    ),
    "Table1.md": (
        "42_table1.py",
        "python3 42_table1.py --out Table1",
        "",
        "same run as Table1.tsv -- 42 writes the pair. Use --merge.",
    ),
    "f5_rl125.png": (
        "90_build_f5.py",
        "python3 90_build_f5.py --out legend_check/f5_rl125 "
        "--model 'DCM/cell' --matching '3-D match (+|disease beta|)' "
        "--expr-model composition --b1 dcm_bams.txt --b2 nf_bams.txt",
        "F5_SCAN=signature_rbp_margins_rl125.tsv;"
        "F5_FWD=rbm20_margin_resolve_rl125.tsv;"
        "F5_REC_SCAN=SE=PTBP1/0.062/-0.025/0.122,MXE=PTBP1/0.113/0.016/0.187;"
        "F5_FLOOR_RANK=SE=6,MXE=3;"
        "F5_REC_FWD=SE=-0.051/-0.159/0.047,MXE=-0.019/-0.108/0.072;"
        "F5_TECH45=SE=0.96/0.87/1.02/0.92,MXE=1.98/1.24/2.65/1.63;"
        "F5_TECH9=SE=1.03/0.95/1.08/1.00,MXE=1.95/1.29/2.52/1.67;"
        "F5_NTARGET=SE=641,MXE=1399;"
        "F5_D_LOG=rbm20_g45_407951.out",
        "values read out of 90's own comment block (lines 91-203); the eight "
        "TECH pairs, NTARGET and D_LOG all match "
        "RESULT_f5_panelD_and_provenance_rl125.md exactly. Flags from 89 "
        "(--model, --matching, --b1, --b2) and --expr-model composition is the "
        "default and the model the Results quote. ⚠ F5_REC_SCAN SAYS THE TOP "
        "GENE IS PTBP1 AT BOTH SE AND MXE -- at 50 bp it was QKI/RBM20. The "
        "gate is keyed by gene, so a 50 bp record over 125 bp data labels the "
        "wrong gene and PASSES. Check the figure, the legend and the "
        "Discussion for 'top-ranked'. ⚠ F5_SCAN/F5_FWD do NOT refuse -- they "
        "fall back to the 50 bp tables.",
    ),
    "Table2_rl125.tsv": (
        "61_table2_lead_genes.py",
        "python3 61_table2_lead_genes.py --disc SE_rl125 "
        "--compadj '{ET}_v2compadj_rl125' "
        "--prefix Table2_rl125 --types SE --top 20",
        "FROZEN_SE=813",
        "61 line 62: FROZEN_SE = int(os.environ.get('FROZEN_SE', '926')); "
        "line 61 says '926 at 50 bp, 813 at 125'. It is the ONLY self-check in "
        "61 that aborts, and it compares len(sig) for SE -- so it is also an "
        "independent confirmation that the right discovery table was read. "
        "⚠ --compadj MUST be passed: its default '{ET}_v2compadj' is "
        "UNSUFFIXED, so the CompAdj column would be built from the 50 bp table "
        "while FROZEN_SE passes on the discovery count alone. "
        "SE_v2compadj_rl125.all_events.tsv exists and is the right one: the "
        "125 bp compadj arm is 777 significant (498 of the 813). "
        "SE_v2cmfrac_rl125 also exists and gives 939 -- that is the "
        "cardiomyocyte/fibroblast FRACTION model the Results quote separately, "
        "a different adjustment. 939 is not 777; do not substitute one for the "
        "other. The {ET} form mirrors the default so --types is still what "
        "picks the event class.",
    ),
    "Table2_rl125.md": (
        "61_table2_lead_genes.py",
        "python3 61_table2_lead_genes.py --disc SE_rl125 "
        "--compadj '{ET}_v2compadj_rl125' "
        "--prefix Table2_rl125 --types SE --top 20",
        "FROZEN_SE=813",
        "same run as Table2_rl125.tsv. Use --merge.",
    ),
    "Table4_rl125.tsv": (
        "98_build_replication_table.py",
        "python3 98_build_replication_table.py --apply",
        "TBL_SUFFIX=_rl125;"
        "REC_KN=Colorado=311/385,Brigham=298/364,Stanford=345/455,Fuwai=422/596;"
        "REC_SLOPE=Colorado=0.5999,Brigham=0.5099,Stanford=0.4008,Fuwai=0.2948;"
        "REC_POOLED=1376,1800;REC_POOLED_P=2.51e-117;"
        "REC_MEAN_SLOPE=0.451;REC_LBNL=0.5726/78.2",
        "RESULT_table4_rl125.md, section 'THE RECORD, FOR RE-RUNS' -- verbatim. "
        "⚠ 98's DEFAULTS BUILD THE 50 bp TABLE: --list with no environment "
        "reads cohort_type_concordance_dedup.tsv and cohort_slopes_dedup.tsv "
        "and would write Table4.tsv / Table4.md (pooled 1636/2163 = 75.6%, "
        "mean slope 0.406). The 125 bp table is pooled 1376/1800 = 76.4%, mean "
        "0.451. TBL_SUFFIX is what moves both the inputs and the output names.",
    ),
    "Table4_rl125.md": (
        "98_build_replication_table.py",
        "python3 98_build_replication_table.py --apply",
        "TBL_SUFFIX=_rl125;"
        "REC_KN=Colorado=311/385,Brigham=298/364,Stanford=345/455,Fuwai=422/596;"
        "REC_SLOPE=Colorado=0.5999,Brigham=0.5099,Stanford=0.4008,Fuwai=0.2948;"
        "REC_POOLED=1376,1800;REC_POOLED_P=2.51e-117;"
        "REC_MEAN_SLOPE=0.451;REC_LBNL=0.5726/78.2",
        "same run as Table4_rl125.tsv. Use --merge.",
    ),
    "Table3_rl125.tsv": (
        "87_build_table3.py",
        "python3 87_build_table3.py --dir t3_rl125 --out Table3_rl125",
        "T3_RECORD=rl125",
        "RESULT_table3_rl125_built_and_caption_fixed.md -- 167 wrote the input "
        "directory t3_rl125/ and 168 gave 87 a 125 bp copy of its internal "
        "record. 'rl125 on the command line without T3_RECORD=rl125 exits "
        "before G1, and T3_RECORD=rl125 without a 125 bp path also exits.' "
        "⚠ CHECK whether --out is joined to --dir: the submitted file is "
        "t3_rl125/Table3_rl125.tsv. If 87 writes to the cwd, pass "
        "--out t3_rl125/Table3_rl125 instead.",
    ),
    "Table3_rl125.md": (
        "87_build_table3.py",
        "python3 87_build_table3.py --dir t3_rl125 --out Table3_rl125",
        "T3_RECORD=rl125",
        "same run as Table3_rl125.tsv. Use --merge.",
    ),
}

# Deliberately absent, with the reason. 173 prints these so a blank row is a
# known gap rather than an oversight.
UNKNOWN = {}

# The same 61/73 run also writes these. They are NOT in submitted_artefacts.txt;
# if they are submitted, add them to this row's `outputs` so one run verifies
# all of them, rather than giving them rows of their own.
ALSO_WRITTEN = {
    "Table2_rl125.tsv": ["Table2_recount_rl125.tsv",
                         "Table2_rl125_top40.tsv", "Table2_rl125_top40.md"],
}




# Artefacts where the project record beats discovery, with the reason. 164 is a
# real writer -- of master_f2_panelB_rl125.tsv, which 39 READS via --master --
# so discovery is not wrong about the idiom, only about which file.
RECORD_WINS = {
    "F2_discovery_rl125.png":
        "164 writes the panel-B master that 39 reads; it is one step upstream",
}


# 172 writes this into `notes` when it copied a command out of another script.
HARVESTED = "command read from"


def is_blank(col, val, notes=""):
    """A 172 placeholder counts as blank, so the record can fill it.

    `F1_SLOPES=` and `F6_COUNTS=;F6_TOTALS=` are names with no values, written
    by discovery for you to complete. Treating them as set would make --force
    necessary for the normal case, and --force also overwrites real edits.

    A HARVESTED command counts as blank too. 172 copied
    `python 87_build_table3.py` out of another script; it is syntactically a
    command and it is WRONG -- no `--dir t3_rl125`, no `--out Table3_rl125`, so
    it would build the 50 bp table in the working directory. A harvested
    command is a guess; the record is evidence. Only a command you wrote
    yourself is protected.
    """
    v = (val or "").strip()
    if col == "command" and v and HARVESTED in (notes or ""):
        return True
    if not v:
        return True
    if col == "env":
        return all(not kv.split("=", 1)[1].strip()
                   for kv in v.split(";") if "=" in kv) or "=" not in v
    if col == "producer":
        return "|" in v                     # an unnarrowed candidate list
    if col == "command":
        # prose, or a bare script with no interpreter, is not a command
        return not re.match(r"^(python\d?|bash|sh|Rscript)\s+\S", v)
    return False


def parse_env(v):
    out = {}
    for kv in (v or "").split(";"):
        kv = kv.strip()
        if not kv:
            continue
        k, _, val = kv.partition("=")
        out[k.strip()] = val.strip()
    return out


def merge_env(existing, record):
    """Union the record's variables with what the row already has.

    Replacing wholesale loses a value set by hand; refusing loses the variables
    the record adds. 172 pre-fills `TBL_SUFFIX=` as a placeholder and the
    record supplies seven variables -- so a partial env must GROW, not block.
    A non-empty existing value wins for a key both define, and is reported.
    """
    e, r = parse_env(existing), parse_env(record)
    kept = []
    out = dict(r)
    for k, v in e.items():
        if v and k in r and v != r[k]:
            out[k] = v
            kept.append(k)
        elif k not in r and v:
            out[k] = v
            kept.append(k)
    order = [k for k in r] + [k for k in out if k not in r]
    return ";".join(f"{k}={out[k]}" for k in order), kept


def order_rows(rows):
    """Put a producer before anything that reads its output.

    F1's command is `43_figure1_design.py --table Table1.tsv`, and Table 1's own
    run rewrites Table1.tsv. In spec order F1 ran FIRST, so it was drawn from
    the copy that was about to be replaced -- and if 42's output had moved, the
    figure would have been built on the old table and still reported IDENTICAL,
    because its own backup was taken from the same stale file.

    Dependencies are read off the commands: a row depends on another when its
    command names one of that row's outputs.
    """
    produces = {}
    for i, r in enumerate(rows):
        for o in r["outputs"].split(","):
            o = Path(o.strip()).name
            if o:
                produces.setdefault(o, i)

    deps = {i: set() for i in range(len(rows))}
    for i, r in enumerate(rows):
        cmd = r["command"] + " " + r["env"]
        for name, j in produces.items():
            if j != i and re.search(r"(?<![\w./-])" + re.escape(name)
                                    + r"(?![\w])", cmd):
                deps[i].add(j)

    out, seen, guard = [], set(), 0
    while len(out) < len(rows) and guard <= len(rows) + 1:
        guard += 1
        for i in range(len(rows)):
            if i in seen or not deps[i] <= seen:
                continue
            out.append(i)
            seen.add(i)
    if len(out) != len(rows):
        # a cycle: keep the original order rather than invent one
        return rows, None
    moved = [(rows[i], sorted(Path(rows[j]["outputs"].split(",")[0]).name
                              for j in deps[i]))
             for i in out if deps[i]]
    return [rows[i] for i in out], moved


def row_names(artefact):
    """Basenames a row covers.

    --merge rewrites `artefact` to "<path> + <name> + ...", so a second run saw
    "Table1.md + Table1.tsv", found it in neither KNOWN nor UNKNOWN, and
    reported three already-filled rows as still blank.
    """
    out = []
    for part in str(artefact).split(" + "):
        part = part.strip()
        if part:
            out.append(Path(part).name)
    return out


def read_spec(path):
    rows = []
    with open(path, encoding="utf-8") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        if hdr[:len(SPEC_COLS)] != SPEC_COLS:
            sys.exit(f"[fatal] {path}: header is\n  {hdr}\nexpected\n  "
                     f"{SPEC_COLS}\nRe-run 172 --discover with the current "
                     f"version -- the `env_optional` column was added so an "
                     f"override that falls back silently cannot be read past.")
        for line in fh:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            parts += [""] * (len(SPEC_COLS) - len(parts))
            rows.append(dict(zip(SPEC_COLS, parts[:len(SPEC_COLS)])))
    return rows


def write_spec(path, rows):
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\t".join(SPEC_COLS) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(c, "")).replace("\t", " ")
                               for c in SPEC_COLS) + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--spec", default="verify_spec.tsv")
    ap.add_argument("--apply", action="store_true")
    ap.add_argument("--merge", action="store_true",
                    help="join rows with an identical command into one row "
                         "with several outputs")
    ap.add_argument("--force", action="store_true",
                    help="overwrite fields that already have a value")
    a = ap.parse_args()

    spec = Path(a.spec)
    if not spec.is_file():
        sys.exit(f"[fatal] no spec at {spec.resolve()}. Run 172 --discover first.")
    rows = read_spec(spec)
    print(f"{spec}: {len(rows)} row(s)")

    filled, kept, missing, surprises = [], [], [], []
    for r in rows:
        names = row_names(r["artefact"])
        base = next((n for n in names if n in KNOWN),
                    next((n for n in names if n in UNKNOWN), names[0]))
        if base in KNOWN:
            prod, cmd, env, why = KNOWN[base]
            if (r["producer"] and r["producer"] != prod
                    and not is_blank("producer", r["producer"], r["notes"])):
                if base in RECORD_WINS:
                    surprises.append(
                        f"{base}: discovery said {r['producer']}, replaced with "
                        f"{prod} -- {RECORD_WINS[base]}")
                elif not a.force:
                    surprises.append(
                        f"{base}: discovery says {r['producer']}, the record "
                        f"says {prod}. NOT changed -- resolve this first, one "
                        f"of the two is wrong.")
                    kept.append(base)
                    continue
            if env:
                opt = parse_env(r.get("env_optional", ""))
                if opt:
                    rec = parse_env(env)
                    moved = [k for k in opt if k in rec]
                    if moved:
                        r["env_optional"] = ";".join(
                            f"{k}={v}" for k, v in opt.items()
                            if k not in rec)
                        surprises.append(
                            f"{base}.env: {', '.join(moved)} moved from "
                            f"env_optional to env — the record supplies a "
                            f"value, so the default must NOT be used")
                newenv, kept = merge_env(r["env"], env)
                if newenv != r["env"]:
                    r["env"] = newenv
                if kept:
                    surprises.append(
                        f"{base}.env: kept your value(s) for "
                        f"{', '.join(kept)} and added the rest from the record")
                env = ""                 # handled; do not run the generic path
            for col, val in (("producer", prod), ("command", cmd), ("env", env)):
                if not val:
                    # An empty `env` in the record is a POSITIVE statement: this
                    # producer needs no overrides. It must be able to CLEAR a
                    # placeholder, or Table1.tsv keeps the `F1_SLOPES=` that 172
                    # inherited from 43 being a candidate -- and 172 then refuses
                    # the row for an empty env value.
                    if col == "env" and r[col] and is_blank("env", r[col],
                                                             r["notes"]):
                        r[col] = ""
                    continue
                if r[col] and r[col] != val \
                        and not is_blank(col, r[col], r["notes"]) \
                        and not a.force and base not in RECORD_WINS:
                    surprises.append(f"{base}.{col}: already set to "
                                     f"{r[col][:60]!r}; --force to replace")
                    continue
                r[col] = val
            r["notes"] = (r["notes"] + " | FILLED FROM RECORD: " + why).strip(" |")
            filled.append(base)
        elif base in UNKNOWN:
            r["notes"] = (r["notes"] + " | NOT FILLED: " + UNKNOWN[base]).strip(" |")
            missing.append(base)
        else:
            missing.append(base)

    # ⚠ EVERY figure producer writes a .pdf beside the .png. 172 backs up only
    # the paths in `outputs`, so a verify run overwrote six PDFs with no copy
    # kept -- and for F4 that left the .pdf regenerated while the .png was
    # restored to the submitted version, a mismatched pair on disk. The PDF is
    # the file a journal publishes; the PNG is the one you look at.
    SIBLINGS = {".png": (".pdf", ".svg", ".eps"),
                ".pdf": (".png",),
                ".tsv": (".md", ".xlsx"),
                ".md": (".tsv",)}
    for r in rows:
        outs = [o.strip() for o in r["outputs"].split(",") if o.strip()]
        found = []
        for o in list(outs):
            for ext in SIBLINGS.get(Path(o).suffix.lower(), ()):
                sib = str(Path(o).with_suffix(ext))
                if sib not in outs and Path(sib).is_file():
                    outs.append(sib)
                    found.append(Path(sib).name)
        if found:
            r["outputs"] = ",".join(outs)
            r["notes"] = (r["notes"] + " | sibling output(s) added to the "
                          "backup set: " + ", ".join(found)).strip(" |")
            surprises.append(
                f"{row_names(r['artefact'])[0]}: added sibling(s) "
                f"{', '.join(found)} — the producer writes them too and they "
                f"were not being backed up")

    for r in rows:
        extra = []
        for n in row_names(r["artefact"]):
            extra += ALSO_WRITTEN.get(n, [])
        if not extra:
            continue
        base = str(Path(r["outputs"].split(",")[0].strip()).parent)
        outs = [o.strip() for o in r["outputs"].split(",") if o.strip()]
        added = []
        for e in extra:
            cand = str(Path(base) / e) if base not in ("", ".") else e
            if cand not in outs:
                outs.append(cand)
                added.append(e)
        if added:
            r["outputs"] = ",".join(outs)
            r["notes"] = (r["notes"] + " | also backed up and compared: "
                          + ", ".join(added)).strip(" |")
            surprises.append(
                f"{row_names(r['artefact'])[0]}: added {len(added)} file(s) this "
                f"run also writes to `outputs` so they are backed up — "
                f"{', '.join(added)}")

    if a.merge:
        by_cmd, order = {}, []
        for r in rows:
            k = (r["command"], r["env"])
            if r["command"] and k in by_cmd:
                tgt = by_cmd[k]
                outs = [o.strip() for o in tgt["outputs"].split(",") if o.strip()]
                for o in [o.strip() for o in r["outputs"].split(",") if o.strip()]:
                    if o not in outs:
                        outs.append(o)
                tgt["outputs"] = ",".join(outs)
                tgt["artefact"] = tgt["artefact"] + " + " + Path(r["artefact"]).name
            else:
                by_cmd[k] = r
                order.append(r)
        if len(order) != len(rows):
            print(f"\n--merge: {len(rows)} row(s) -> {len(order)} run(s)")
        rows = order

    # Report on what the row HAS, not on whether its name is in a table.
    have = [r for r in rows if r["command"].strip()]
    blank = [r for r in rows if not r["command"].strip()]
    rows, moved = order_rows(rows)
    if moved:
        print("\nordered so a producer runs before its reader:")
        for r, on in moved:
            nm = " + ".join(row_names(r["artefact"]))
            print(f"    {nm}  after  {', '.join(on)}")
    elif moved is None:
        print("\n[warn] the commands form a dependency cycle; spec order kept. "
              "Check by hand which artefact feeds which.")

    print(f"\nfilled from the record : {len(filled)}")
    for b in filled:
        print(f"    {b}")
    if blank:
        print(f"\nstill blank            : {len(blank)}")
        for r in blank:
            nm = " + ".join(row_names(r["artefact"]))
            print(f"    {nm}")
            why = next((UNKNOWN[n] for n in row_names(r["artefact"])
                        if n in UNKNOWN), None)
            if why:
                print(f"      {why[:200]}")
    if len(have) + len(blank) != len(rows):
        print("  [warn] row accounting does not add up")
    if surprises:
        print(f"\n⚠ {len(surprises)} conflict(s) -- NOT changed:")
        for s in surprises:
            print(f"    {s}")

    ready = [r for r in rows if r["command"]]
    # ⚠ WRITE BEFORE REPORTING. `173 --apply | head` closes stdout, SIGPIPE
    # kills the process mid-report, and the write never happens -- while the
    # part of the report already printed says the rows were filled. Writing
    # first makes the output cosmetic and the effect certain.
    bak = None
    if a.apply:
        bak = spec.with_suffix(spec.suffix + ".bak."
                               + time.strftime("%Y%m%d_%H%M%S"))
        shutil.copy2(spec, bak)
        write_spec(spec, rows)

    print(f"\n{len(ready)} of {len(rows)} row(s) now have a command.")
    if not a.apply:
        print("\nDRY RUN -- nothing written. Re-run with --apply.")
        return

    print(f"\nwrote {spec}   (previous version kept at {bak})")
    print("\nNext:")
    print(f"  python 172_verify_regeneration.py --root . --verify --spec {spec}")
    print("  # read the command list, then add --run")


if __name__ == "__main__":
    main()
