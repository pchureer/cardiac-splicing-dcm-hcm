#!/usr/bin/env python3
"""
163_patch_38_ceiling_labels.py -- Figure 4's ceiling labels collide at 125 bp.

At 50 bp the two ceilings sat at 0.584 and 0.988, far enough apart that two
centred labels cleared each other. At 125 bp they are 0.690 and 0.968, and the
rendered figure reads

    platform 0.69people 0.97

Both labels are drawn with `ha="center"` on the same baseline, so the gap
between them is the gap between the lines. Re-basing the ceilings was correct;
the layout assumed a spacing that the correct values do not have.

THE FIX
-------
Anchor the lower ceiling's label to the RIGHT of its text (so it grows
leftward) and the higher one to the LEFT (so it grows rightward), each nudged
a little away from its line. The labels then move apart as the lines move
together, which is the opposite of the current behaviour.

This is layout only. No number changes, and 162's work is untouched.

⚠ Its own script says so: "NOW LOOK AT THE PNG. The validator checks colour,
not layout -- check for label collisions, clipped text and overlapping
annotations." The validator passed on the collided figure.

USAGE
-----
  python3 163_patch_38_ceiling_labels.py [--apply] [--self-test]
"""
import argparse
import difflib
import os
import shutil
import sys
import time

TARGET = "38_figure4_replication.py"

A = '''        ax.text(xv, ypooled, lab, fontsize=6.3, color=col, ha="center", va="center",'''
B = '''        # 163: at 125 bp the ceilings are 0.69 and 0.97 and two centred labels
        # overlap. Anchor the lower one right and the upper one left so they
        # separate as the lines converge.
        _lo = min(CEIL_PLATFORM, CEIL_PEOPLE)
        _ha = "right" if xv == _lo else "left"
        _dx = -0.015 if _ha == "right" else 0.015
        ax.text(xv + _dx, ypooled, lab, fontsize=6.3, color=col, ha=_ha,
                va="center",'''


def apply_edits(text):
    if "163:" in text:
        return text, ["38 already carries 163 (restore the .bak to re-apply)"]
    n = text.count(A)
    if n != 1:
        return text, [f"anchor occurs {n}x, expected 1x: {A.strip()[:60]}"]
    return text.replace(A, B, 1), []


def verify(text):
    bad = []
    for w in ('_lo = min(CEIL_PLATFORM, CEIL_PEOPLE)',
              '_ha = "right" if xv == _lo else "left"',
              'ax.text(xv + _dx, ypooled, lab'):
        if w not in text:
            bad.append(f"missing: {w[:50]}")
    if 'ax.text(xv, ypooled, lab, fontsize=6.3, color=col, ha="center"' in text:
        bad.append("the centred draw call survived")
    # 162 must be intact
    for w in ('_C4 = os.environ.get("F4_CEILINGS", "").strip()',
              'f"platform {CEIL_PLATFORM:.2f}"',
              'f"people {CEIL_PEOPLE:.2f}"'):
        if w not in text:
            bad.append(f"162 was damaged: {w[:46]}")
    return bad


FIX = '''import os
ORANGE="#eb6834"; BLUE="#2a78d6"
_C4 = os.environ.get("F4_CEILINGS", "").strip()
if _C4:
    CEIL_PLATFORM, CEIL_PEOPLE = (float(_x) for _x in _C4.split("/"))
else:
    CEIL_PEOPLE = 0.988
    CEIL_PLATFORM = 0.584


def draw(ax, ypooled):
    for xv, lab, col in ((CEIL_PLATFORM, f"platform {CEIL_PLATFORM:.2f}", ORANGE),
                         (CEIL_PEOPLE, f"people {CEIL_PEOPLE:.2f}", BLUE)):
        ax.axvline(xv, color=col, linewidth=1.3, linestyle=(0, (4, 3)), zorder=1)
        ax.text(xv, ypooled, lab, fontsize=6.3, color=col, ha="center", va="center",
                zorder=5)
'''


def self_test():
    fails = []
    new, probs = apply_edits(FIX)
    if probs:
        fails.append(f"fixture did not patch: {probs}")
    else:
        bad = verify(new)
        if bad:
            fails.append(f"verify rejected its own output: {bad}")
        try:
            compile(new, "x", "exec")
        except SyntaxError as e:
            fails.append(f"does not compile: {e}")
        # the lower ceiling must anchor right, the upper left, at BOTH read lengths
        for ceil, lo, hi in ((None, 0.584, 0.988), ("0.690/0.968", 0.690, 0.968)):
            if ceil: os.environ["F4_CEILINGS"] = ceil
            else: os.environ.pop("F4_CEILINGS", None)
            ns = {}
            exec(compile(new, "x", "exec"), ns)
            got = []
            class Ax:
                def axvline(self, *a, **k): pass
                def text(self, x, y, lab, **k): got.append((round(x, 3), k["ha"]))
            ns["draw"](Ax(), 0)
            want = [(round(lo - 0.015, 3), "right"), (round(hi + 0.015, 3), "left")]
            if got != want:
                fails.append(f"ceilings {lo}/{hi}: got {got}, want {want}")
        os.environ.pop("F4_CEILINGS", None)
        again, p2 = apply_edits(new)
        if not p2 or "already carries 163" not in p2[0]:
            fails.append("re-application did not refuse")
    print("163 self-test")
    for f in fails:
        print("  ✗ " + f)
    if not fails:
        print("  ✓ patches and compiles, the lower ceiling anchors right and")
        print("    the upper left at BOTH read lengths, 162 is intact, and")
        print("    re-application refuses")
    return 1 if fails else 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("target", nargs="?", default=TARGET)
    ap.add_argument("--apply", action="store_true")
    ap.add_argument("--self-test", action="store_true")
    a = ap.parse_args()
    if a.self_test:
        return self_test()
    if self_test():
        print("REFUSING: 163's own self-test failed.")
        return 3
    print()
    if not os.path.exists(a.target):
        print(f"REFUSING: {a.target} not found")
        return 3
    src = open(a.target, errors="replace").read()
    new, probs = apply_edits(src)
    print("=" * 74)
    print("163 — Figure 4's ceiling labels")
    print("=" * 74)
    if probs:
        print("REFUSING:")
        for p in probs:
            print("    " + p)
        return 3
    print("".join(difflib.unified_diff(src.splitlines(True), new.splitlines(True),
                                       fromfile=a.target,
                                       tofile=a.target + " (patched)", n=3)))
    if not a.apply:
        print("DRY RUN. Re-run with --apply, then:")
        print()
        print("  # 50 bp control: must still be byte-identical to the published")
        print("  #   figure, because at 0.58/0.99 the labels never collided and")
        print("  #   a 0.015 nudge is invisible at that separation... it is NOT.")
        print("  #   Expect a DIFFERENCE here, and check it is only the labels.")
        print("  python3 38_figure4_replication.py --out /tmp/f4_cal2")
        print("  md5sum /tmp/f4_cal2.png F4_replication.png")
        print()
        print("  # 125 bp")
        print("  F4_SUFFIX=_rl125 F4_CEILINGS=0.690/0.968 \\")
        print("  python3 38_figure4_replication.py \\")
        print("      --conc direction_concordance_dedup_rl125.tsv \\")
        print("      --out F4_replication_rl125")
        print()
        print("⚠ The 50 bp figure WILL change by 0.015 of an axis unit in two")
        print("  label positions. If the published 50 bp figure must stay")
        print("  byte-identical, regenerate it from the .bak instead.")
        return 0
    bak = f"{a.target}.bak.{time.strftime('%Y%m%d_%H%M%S')}"
    shutil.copy2(a.target, bak)
    open(a.target, "w").write(new)
    print(f"[write] {a.target}   (backup: {bak})")
    disk = open(a.target, errors="replace").read()
    bad = verify(disk)
    try:
        compile(disk, a.target, "exec")
    except SyntaxError as e:
        bad.append(f"does not compile: {e}")
    if bad or disk != new:
        print("✗ ON-DISK VERIFICATION FAILED — restoring:")
        for b in bad:
            print(f"    {b}")
        shutil.copy2(bak, a.target)
        return 3
    print("✓ verified: the lower ceiling's label anchors right, the upper's")
    print("  left, and 162's parameterisation is untouched.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
