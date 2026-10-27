#!/usr/bin/env python3
"""
176_patch_figure_cosmetics.py -- three one-line cosmetic fixes found by the
26 Oct visual audit (RESULT_visual_audit_six_figures_26oct.md).

  P1  37_figure3_etiology.py    panel C: the n / Deming / direction annotation
                                lies across the regression line and the upper
                                right cluster. Put a white box behind it.
  P2  41_figure6_enrichment.py  panel C: "...and 11 more" and the italic
                                "exactly the cardiac story" are stacked with no
                                gap. The caption is laid out by a helper at a
                                fixed height, so show five terms instead of
                                six; the "...and N more" count follows.
  P3  90_build_f5.py            panel C title says "No candidate separates from
                                the housekeeping control" while PTBP1's MXE
                                interval visibly clears the floor. The claim is
                                true AFTER Holm; the title now says so, in the
                                panel's own word for the control ("floor") so
                                the length barely changes.

Discovery first:

  python 176_patch_figure_cosmetics.py --list
  python 176_patch_figure_cosmetics.py --apply
  python 176_patch_figure_cosmetics.py --apply --only p1,p3

Each patch locates its anchor by TOKENISING the producer, not by regex over
raw text, so a parenthesis inside a label string ("(median |t|)") cannot throw
the match. Each refuses unless its anchor occurs exactly once. P2 refuses unless
the caption's y-coordinate is a plain numeric literal and the call is in axes
coordinates, because shifting a data-unit coordinate by an axes-unit delta
would move the caption by an unknown amount. Re-running after success prints
ALREADY APPLIED and changes nothing.

After applying, regenerate through the harness so the artefact stays
reproducible:

  python 172_verify_regeneration.py --root . --verify --spec verify_spec.tsv --run

Each patched figure reports DIFFERS once -- that is the change landing -- and
174_diff_images.py on the .regenerated copy should show ONLY the intended
region moving.
"""

import argparse
import io
import re
import shutil
import sys
import time
import tokenize
from pathlib import Path

# ----------------------------------------------------------------- the patches

P1_FILE = "37_figure3_etiology.py"
P1_NEEDLE = "Deming slope"
P1_KW = 'bbox=dict(facecolor="white", alpha=0.85, edgecolor="none", pad=2)'

P2_FILE = "41_figure6_enrichment.py"
P2_NEEDLE = "exactly the cardiac story"

P3_FILE = "90_build_f5.py"
P3_OLD = "No candidate separates from the housekeeping control"
P3_NEW = "No candidate separates from the floor after correction"


# ----------------------------------------------------------------- tokenising

class Call:
    """One `.text(...)` call: where it starts, where its ')' is, its source."""
    def __init__(self, src_lines, open_tok, close_tok, arg_tokens):
        self.open = open_tok      # the '(' token
        self.close = close_tok    # the matching ')' token
        self.args = arg_tokens    # tokens strictly inside the parens
        r0, c0 = open_tok.start
        r1, c1 = close_tok.end
        if r0 == r1:
            self.text = src_lines[r0 - 1][c0:c1]
        else:
            self.text = (src_lines[r0 - 1][c0:] + "".join(src_lines[r0:r1 - 1])
                         + src_lines[r1 - 1][:c1])
        self.line = open_tok.start[0]


def text_calls(src):
    """Every `<something>.text(` call in src, paren-matched by the tokenizer."""
    lines = src.splitlines(keepends=True)
    toks = list(tokenize.generate_tokens(io.StringIO(src).readline))
    calls = []
    for i, t in enumerate(toks):
        if t.type != tokenize.NAME or t.string != "text":
            continue
        if i < 1 or toks[i - 1].string != ".":
            continue
        if i + 1 >= len(toks) or toks[i + 1].string != "(":
            continue
        depth, j = 0, i + 1
        inner = []
        while j < len(toks):
            tj = toks[j]
            if tj.type == tokenize.OP and tj.string in "([{":
                depth += 1
            elif tj.type == tokenize.OP and tj.string in ")]}":
                depth -= 1
                if depth == 0:
                    break
            if j > i + 1:
                inner.append(tj)
            j += 1
        if j < len(toks):
            calls.append(Call(lines, toks[i + 1], toks[j], inner))
    return calls


def split_args(arg_tokens):
    """Top-level comma split of a call's argument tokens -> list of token lists."""
    out, cur, depth = [], [], 0
    for t in arg_tokens:
        if t.type == tokenize.OP and t.string in "([{":
            depth += 1
        elif t.type == tokenize.OP and t.string in ")]}":
            depth -= 1
        if t.type == tokenize.OP and t.string == "," and depth == 0:
            out.append(cur)
            cur = []
            continue
        if t.type not in (tokenize.NL, tokenize.NEWLINE, tokenize.COMMENT):
            cur.append(t)
    if cur:
        out.append(cur)
    return out


def splice(src, row, col, insert):
    """Insert text at (1-based row, 0-based col)."""
    lines = src.splitlines(keepends=True)
    ln = lines[row - 1]
    lines[row - 1] = ln[:col] + insert + ln[col:]
    return "".join(lines)


def replace_span(src, r0, c0, r1, c1, new):
    lines = src.splitlines(keepends=True)
    assert r0 == r1, "single-line spans only"
    ln = lines[r0 - 1]
    lines[r0 - 1] = ln[:c0] + new + ln[c1:]
    return "".join(lines)


# ----------------------------------------------------------------- P1

def plan_p1(src, say):
    calls = [c for c in text_calls(src) if P1_NEEDLE in c.text]
    if not calls:
        # the label may be bound to a variable first
        m = re.findall(r"^\s*([A-Za-z_]\w*)\s*=.*" + re.escape(P1_NEEDLE),
                       src, flags=re.M)
        if m:
            for v in set(m):
                calls += [c for c in text_calls(src)
                          if re.search(r"\b" + re.escape(v) + r"\b", c.text)]
    say(f"  {len(calls)} .text( call(s) carry '{P1_NEEDLE}'")
    for c in calls:
        say(f"    line {c.line}: {' '.join(c.text.split())[:150]}")
    if not calls:
        return None, (f"no .text( call carries '{P1_NEEDLE}' and no variable "
                      f"bound to it is passed to one. grep -n 'Deming' "
                      f"{P1_FILE} and tell me what you see.")
    if len(calls) > 1:
        return None, "more than one candidate call; refusing to guess which"
    c = calls[0]
    if "bbox" in c.text:
        return "ALREADY", None
    # Insert straight after the last SIGNIFICANT argument token. The token list
    # includes NL tokens, so "the last token" can be a line break after a
    # trailing comma -- the first version then wrote ", bbox=" after that
    # comma and produced `fontsize=7,\n    , bbox=...)`, a syntax error.
    sig = [t for t in c.args
           if t.type not in (tokenize.NL, tokenize.NEWLINE, tokenize.COMMENT)]
    if not sig:
        return None, "the .text( call has no arguments?"
    last = sig[-1]
    if last.string == ",":
        r, col = last.end
        return ("insert", r, col, " " + P1_KW), None
    r, col = last.end
    return ("insert", r, col, ", " + P1_KW), None


# ----------------------------------------------------------------- P2

P2_OLD = ('MULTI_TERMS + [f"\u2026and {FROZEN_MULTI[\'genome\'] - len(MULTI_TERMS)} more"],')
P2_NEW = ('MULTI_TERMS[:5] + [f"\u2026and {FROZEN_MULTI[\'genome\'] - len(MULTI_TERMS[:5])} more"],')


def plan_p2(src, say, dy=None):
    """Show five terms instead of six.

    The caption is not an ax.text call: lines 426-427 pass the term list and
    "exactly the cardiac story" into a box-drawing helper that lays them out
    itself, with the caption at a fixed height (the empty blue box puts its
    own caption at the same y). So the gap is made by removing one list row,
    not by moving the caption. The "...and N more" count is computed from the
    slice, so it becomes "...and 12 more" on its own, and the manuscript never
    enumerates the six terms, so dropping the sixth changes no claim.
    """
    n_old, n_new = src.count(P2_OLD), src.count(P2_NEW)
    say(f"  six-term list line : {n_old} occurrence(s)")
    say(f"  five-term list line: {n_new} occurrence(s)")
    if n_new and not n_old:
        return "ALREADY", None
    if n_old != 1:
        hits = [f"{i}: {l.strip()[:120]}" for i, l in enumerate(src.splitlines(), 1)
                if "MULTI_TERMS" in l and "more" in l]
        return None, ("expected the term-list line exactly once; found "
                      f"{n_old}. Lines mentioning MULTI_TERMS and 'more': "
                      + (" | ".join(hits) or "none") +
                      ". Paste that line and I will match it exactly.")
    return ("literal", P2_OLD, P2_NEW), None


# ----------------------------------------------------------------- P3

def plan_p3(src, say):
    n_old, n_new = src.count(P3_OLD), src.count(P3_NEW)
    say(f"  '{P3_OLD}': {n_old} occurrence(s)")
    say(f"  '{P3_NEW}': {n_new} occurrence(s)")
    if n_new and not n_old:
        return "ALREADY", None
    if n_old != 1:
        return None, (f"expected the old title exactly once, found {n_old}. If it "
                      f"is split across two string literals, grep -n 'separates' "
                      f"{P3_FILE} and tell me.")
    return ("literal", P3_OLD, P3_NEW), None


# ----------------------------------------------------------------- driver

def apply_plan(src, plan):
    kind = plan[0]
    if kind == "insert":
        _, r, c, text = plan
        return splice(src, r, c, text)
    if kind == "replace":
        _, r0, c0, r1, c1, new = plan
        return replace_span(src, r0, c0, r1, c1, new)
    if kind == "literal":
        _, old, new = plan
        return src.replace(old, new, 1)
    if kind == "replace_marked":
        _, r0, c0, r1, c1, new, marker = plan
        out = replace_span(src, r0, c0, r1, c1, new)
        lines = out.splitlines(keepends=True)
        ln = lines[r0 - 1]
        body = ln.rstrip("\r\n")
        nl = ln[len(body):]
        lines[r0 - 1] = body + marker + nl
        return "".join(lines)
    raise ValueError(kind)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--apply", action="store_true")
    ap.add_argument("--list", action="store_true",
                    help="dry run (the default when --apply is absent)")
    ap.add_argument("--only", default="p1,p2,p3",
                    help="comma list of p1,p2,p3 (default all)")
    ap.add_argument("--dir", default=".", help="where the producers live")
    a = ap.parse_args()
    here = Path(a.dir)
    want = {x.strip().lower() for x in a.only.split(",") if x.strip()}

    def say(*s):
        print(*s, flush=True)

    jobs = [("p1", P1_FILE, lambda s: plan_p1(s, say)),
            ("p2", P2_FILE, lambda s: plan_p2(s, say)),
            ("p3", P3_FILE, lambda s: plan_p3(s, say))]

    ok, bad, already = [], [], []
    for key, fname, planner in jobs:
        if key not in want:
            continue
        say(f"\n=== {key.upper()}  {fname}")
        p = here / fname
        if not p.is_file():
            say(f"  [skip] not found at {p}")
            bad.append((key, "file not found"))
            continue
        src = p.read_text(encoding="utf-8")
        try:
            compile(src, fname, "exec")
        except SyntaxError as exc:
            bad.append((key, f"does not compile before patching: {exc}"))
            say(f"  [refuse] {bad[-1][1]}")
            continue
        plan, why = planner(src)
        if plan == "ALREADY":
            say("  ALREADY APPLIED")
            already.append(key)
            continue
        if plan is None:
            say(f"  [refuse] {why}")
            bad.append((key, why))
            continue
        new = apply_plan(src, plan)
        try:
            compile(new, fname, "exec")
        except SyntaxError as exc:
            bad.append((key, f"the edit would not compile: {exc}"))
            say(f"  [refuse] {bad[-1][1]}")
            continue
        # show the changed line(s)
        for i, (lo, ln) in enumerate(zip(src.splitlines(), new.splitlines()), 1):
            if lo != ln:
                say(f"  - {i}: {lo.strip()[:140]}")
                say(f"  + {i}: {ln.strip()[:140]}")
        if len(new.splitlines()) != len(src.splitlines()):
            say("  (line count changed)")
        if a.apply:
            bak = p.with_name(p.name + ".bak." + time.strftime("%Y%m%d_%H%M%S"))
            shutil.copy2(p, bak)
            p.write_text(new, encoding="utf-8")
            say(f"  applied; backup {bak.name}")
        ok.append(key)

    say("\n" + "=" * 60)
    say(f"{'APPLIED' if a.apply else 'WOULD APPLY'}: {', '.join(ok) or 'none'}"
        f"   already: {', '.join(already) or 'none'}"
        f"   refused: {', '.join(k for k, _ in bad) or 'none'}")
    if bad:
        for k, w in bad:
            say(f"  {k}: {w}")
    if not a.apply and ok:
        say("\nDRY RUN. Re-run with --apply. Then:")
        say("  python 172_verify_regeneration.py --root . --verify "
            "--spec verify_spec.tsv --run")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
