#!/usr/bin/env python3
"""
175_measure_label_ink.py -- read a coloured label's SIZE out of a PNG, when the
vector original is gone and the image cannot be looked at.

Figure 4's ceiling labels are orange and blue. Patch 163 changed only their
ANCHORING, so if that is the whole difference the labels must have the same
WIDTH and the same INK COUNT in both versions and differ only in position. A
different width or a different ink count means the TEXT differs -- a longer or
shorter number.

That distinction is what separates "a label moved" from "the figure carries a
different value", and it needs no OCR and no looking.

  python 175_measure_label_ink.py F4_replication_rl125.png f4_regen.png
  python 175_measure_label_ink.py --box 2000,180,2900,920 a.png b.png

Prints, per colour, the total ink, the row range, and each horizontal block
with its width. Compare the WIDTHS between the two files:

  same widths, different x        -> the label moved; 163's repositioning
  different widths or ink counts  -> the label's TEXT is different
"""

import argparse
import sys


def need(mod):
    try:
        return __import__(mod)
    except ImportError:
        sys.exit(f"[fatal] {mod} is required: conda install -c conda-forge "
                 f"pillow numpy")


def blocks(cols, gap):
    out, start, prev = [], None, None
    for x in cols:
        if start is None:
            start = prev = x
            continue
        if x - prev >= gap:
            out.append((start, prev))
            start = x
        prev = x
    if start is not None:
        out.append((start, prev))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("files", nargs="+")
    ap.add_argument("--box", default=None,
                    help="x0,y0,x1,y1 to restrict to a panel or a label band")
    ap.add_argument("--gap", type=int, default=12,
                    help="blank columns that separate two blocks (default 12)")
    a = ap.parse_args()

    need("PIL")
    np = need("numpy")
    from PIL import Image

    box = None
    if a.box:
        try:
            box = tuple(int(v) for v in a.box.split(","))
            assert len(box) == 4
        except Exception:
            sys.exit("[fatal] --box wants x0,y0,x1,y1")

    rows = {}
    for path in a.files:
        try:
            im = Image.open(path).convert("RGB")
        except Exception as exc:
            print(f"{path}: cannot open: {exc}")
            continue
        if box:
            im = im.crop(box)
        arr = np.asarray(im).astype(int)
        r, g, b = arr[:, :, 0], arr[:, :, 1], arr[:, :, 2]
        masks = {
            # tab:orange is (255,127,14); tab:blue is (31,119,180)
            "orange": (r > g) & (g > b) & (r - b > 60) & (r > 120),
            "blue": (b > g) & (g >= r) & (b - r > 60) & (b > 100),
            "red": (r > g) & (r > b) & (r - g > 80) & (g < 90),
            "green": (g > r) & (g > b) & (g - r > 50),
        }
        print(f"\n{path}   {arr.shape[1]} x {arr.shape[0]}"
              + (f"   (cropped to {box})" if box else ""))
        here = {}
        for nm, m in masks.items():
            n = int(m.sum())
            if n < 20:
                continue
            ys, xs = np.nonzero(m)
            cols = sorted(set(xs.tolist()))
            bl = blocks(cols, a.gap)
            print(f"  {nm:7s} ink {n:7d}   rows {ys.min()}-{ys.max()}   "
                  f"{len(bl)} block(s)")
            widths = []
            for x0, x1 in bl:
                sel = (xs >= x0) & (xs <= x1)
                widths.append(x1 - x0 + 1)
                print(f"            x {x0:5d}-{x1:<5d} width {x1-x0+1:5d}   "
                      f"ink {int(sel.sum()):6d}")
            here[nm] = (n, widths)
        rows[path] = here

    if len(rows) < 2:
        return
    names = list(rows)
    print("\n" + "=" * 62)
    print("COMPARISON")
    for nm in ("orange", "blue", "red", "green"):
        a0 = rows[names[0]].get(nm)
        a1 = rows[names[1]].get(nm)
        if not a0 or not a1:
            continue
        print(f"\n  {nm}")
        print(f"    ink    {a0[0]:>7} vs {a1[0]:>7}"
              f"   ({a1[0] - a0[0]:+d})")
        print(f"    widths {a0[1]} vs {a1[1]}")
        if sorted(a0[1]) == sorted(a1[1]) and a0[0] == a1[0]:
            print("    -> SAME widths and SAME ink: the label MOVED. The text "
                  "is unchanged.")
        elif sorted(a0[1]) == sorted(a1[1]):
            print("    -> same widths, different ink count: antialiasing or a "
                  "sub-pixel shift. Text almost certainly unchanged.")
        else:
            print("    -> DIFFERENT widths: the label's TEXT IS DIFFERENT, not "
                  "just its position. A wider block is a longer number.")
            print("       At ~9 pt a digit is roughly 20-25 px at 300 dpi, so "
                  "divide the width difference by that for the number of "
                  "characters gained or lost.")


if __name__ == "__main__":
    main()
