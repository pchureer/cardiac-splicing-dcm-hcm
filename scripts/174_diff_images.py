#!/usr/bin/env python3
"""
174_diff_images.py -- show WHAT moved between two versions of a figure.

A bounding box tells you where the difference is in pixels. It does not tell
you whether a label shifted, a bar changed height, or a whole panel was added.
This writes images you can look at.

  python 174_diff_images.py \
      --old verify_backup_20261026_111321/F4_replication_rl125.png \
      --new F4_replication_rl125.png \
      --out f4_diff

Writes <out>_sidebyside.png, <out>_heat.png, <out>_crop.png and prints a
per-region breakdown so the change can be attributed to a panel.

Nothing is modified: both inputs are opened read-only.
"""

import argparse
import signal
import sys
from pathlib import Path

# `174 ... | head` closes the pipe early; without this the script ends in a
# BrokenPipeError traceback that looks like a failure in the comparison.
try:
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)
except (AttributeError, ValueError):
    pass


def need(mod):
    try:
        return __import__(mod)
    except ImportError:
        sys.exit(f"[fatal] {mod} is required. "
                 f"conda install -c conda-forge pillow numpy")


def classify_colour(rgb):
    """A rough name for a mean colour, so a component can be recognised.

    Figure 4's ceiling labels are orange and blue (the record: "platform 0.69
    anchors right of its ORANGE line and people 0.97 left of its BLUE line"),
    so naming the hue is often enough to identify what moved.
    """
    r, g, b = [float(v) for v in rgb]
    mx, mn = max(r, g, b), min(r, g, b)
    if mx > 235 and mn > 235:
        return "white/background"
    if mx - mn < 22:
        return "black/grey" if mx < 120 else "light grey"
    if r >= g >= b and r - b > 40:
        return "orange/red" if g < 0.72 * r else "yellow/olive"
    if b >= g >= r and b - r > 40:
        return "blue"
    if g >= r and g >= b:
        return "green"
    if b >= r and r >= g:
        return "purple/magenta"
    return f"mixed({int(r)},{int(g)},{int(b)})"


def components(mask, limit=400000):
    """Connected components of a boolean mask, 8-connectivity, no scipy."""
    import numpy as np
    from collections import deque
    H, W = mask.shape
    seen = np.zeros_like(mask, dtype=bool)
    ys, xs = np.nonzero(mask)
    if len(ys) > limit:
        return None
    out = []
    for sy, sx in zip(ys, xs):
        if seen[sy, sx]:
            continue
        q = deque([(sy, sx)])
        seen[sy, sx] = True
        pts = []
        while q:
            y, x = q.popleft()
            pts.append((y, x))
            for dy in (-1, 0, 1):
                for dx in (-1, 0, 1):
                    ny, nx = y + dy, x + dx
                    if 0 <= ny < H and 0 <= nx < W and mask[ny, nx] \
                            and not seen[ny, nx]:
                        seen[ny, nx] = True
                        q.append((ny, nx))
        out.append(pts)
    out.sort(key=len, reverse=True)
    return out


def report_components(mask, A, B, top=12):
    import numpy as np
    comps = components(mask)
    if comps is None:
        print("\n[skip] too many differing pixels for a component breakdown; "
              "read the crop instead.")
        return
    print(f"\n{len(comps)} connected region(s) of difference. "
          f"Largest {min(top, len(comps))}:")
    print(f"  {'px':>7}  {'box (x0,y0)-(x1,y1)':<30} {'submitted':<18} "
          f"{'regenerated':<18} verdict")
    for pts in comps[:top]:
        yy = np.array([p[0] for p in pts]); xx = np.array([p[1] for p in pts])
        ca = A[yy, xx].mean(axis=0)
        cb = B[yy, xx].mean(axis=0)
        na, nb = classify_colour(ca), classify_colour(cb)
        if na == "white/background" and nb != "white/background":
            verdict = "ADDED in the regenerated figure"
        elif nb == "white/background" and na != "white/background":
            verdict = "REMOVED from the regenerated figure"
        else:
            verdict = "changed colour or shifted"
        print(f"  {len(pts):>7}  "
              f"({xx.min()},{yy.min()})-({xx.max()},{yy.max()})".ljust(34)
              + f" {na:<18} {nb:<18} {verdict}")
    added = sum(len(p) for p in comps
                if classify_colour(A[[q[0] for q in p], [q[1] for q in p]]
                                   .mean(axis=0)) == "white/background")
    removed = sum(len(p) for p in comps
                  if classify_colour(B[[q[0] for q in p], [q[1] for q in p]]
                                     .mean(axis=0)) == "white/background")
    print(f"\n  ink only in the regenerated figure : {added} px")
    print(f"  ink only in the submitted figure   : {removed} px")
    if added and removed:
        print("  both -> something MOVED or was redrawn differently.")
    elif added:
        print("  only additions -> the regenerated figure draws something the "
              "submitted one does not.")
    elif removed:
        print("  only removals -> the submitted figure draws something the "
              "regenerated one does not.")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--old", required=True, help="the submitted copy (backup)")
    ap.add_argument("--new", required=True, help="the regenerated copy")
    ap.add_argument("--out", default="diff", help="output prefix")
    ap.add_argument("--grid", type=int, default=3,
                    help="split into GRID x GRID regions for the breakdown "
                         "(default 3, which maps onto a 2x2 or 2x3 panel "
                         "figure well enough to name the panel)")
    ap.add_argument("--top", type=int, default=12,
                    help="how many of the largest changed regions to describe")
    ap.add_argument("--pad", type=int, default=40,
                    help="pixels of context around the crop")
    a = ap.parse_args()

    need("PIL")
    np = need("numpy")
    from PIL import Image, ImageChops

    for p in (a.old, a.new):
        if not Path(p).is_file():
            sys.exit(f"[fatal] not a file: {p}")

    with Image.open(a.old) as i1, Image.open(a.new) as i2:
        im1 = i1.convert("RGB")
        im2 = i2.convert("RGB")
        if im1.size != im2.size:
            print(f"⚠ SIZE CHANGED: {im1.size} -> {im2.size}")
            print("  A size change means the figure's geometry moved, not just "
                  "its content -- a panel added or removed, or a different "
                  "figsize/dpi. Compare the two side by side and stop looking "
                  "for a moved label.")
            w = max(im1.width, im2.width)
            h = max(im1.height, im2.height)
            c1 = Image.new("RGB", (w, h), "white"); c1.paste(im1, (0, 0))
            c2 = Image.new("RGB", (w, h), "white"); c2.paste(im2, (0, 0))
            im1, im2 = c1, c2

        A = np.asarray(im1).astype(int)
        B = np.asarray(im2).astype(int)
        d = np.abs(A - B).sum(axis=2)
        mask = d > 0
        n = int(mask.sum())
        tot = mask.size
        print(f"{n} of {tot} pixels differ ({100 * n / tot:.3f}%), "
              f"max summed delta {int(d.max())}")
        if n == 0:
            print("pixels are identical; the files differ only in metadata.")
            return

        ys, xs = np.nonzero(mask)
        x0, x1, y0, y1 = xs.min(), xs.max(), ys.min(), ys.max()
        print(f"bounding box  x {x0}-{x1}  y {y0}-{y1}  "
              f"({x1 - x0 + 1} x {y1 - y0 + 1} px)")

        # --- per-region breakdown, so the change can be named
        g = a.grid
        H, W = mask.shape
        print(f"\nper-region share (grid {g}x{g}, row-major; "
              f"figure is {W}x{H}px):")
        rows = []
        for r in range(g):
            cells = []
            for c in range(g):
                sub = mask[r * H // g:(r + 1) * H // g,
                           c * W // g:(c + 1) * W // g]
                pct = 100 * sub.sum() / max(1, sub.size)
                cells.append(pct)
            rows.append(cells)
        for r, cells in enumerate(rows):
            band = "top" if r == 0 else ("bottom" if r == g - 1 else "middle")
            print("  " + f"{band:7s}" +
                  "".join(f"{v:9.2f}%" for v in cells))
        print("           " + "".join(f"{lab:>10s}" for lab in
                                      (["left", "centre", "right"] if g == 3
                                       else [f"col{i}" for i in range(g)])))
        hot = max(((v, r, c) for r, row in enumerate(rows)
                   for c, v in enumerate(row)))
        print(f"  busiest region: row {hot[1]}, col {hot[2]} "
              f"({hot[0]:.2f}% of its pixels)")

        report_components(mask, A, B, top=a.top)

        # --- heat map of the difference
        heat = (255 * (d / max(1, d.max()))).astype("uint8")
        Image.fromarray(255 - heat).save(f"{a.out}_heat.png")

        # --- side by side, old | new
        sbs = Image.new("RGB", (im1.width * 2 + 24, im1.height), "white")
        sbs.paste(im1, (0, 0))
        sbs.paste(im2, (im1.width + 24, 0))
        sbs.save(f"{a.out}_sidebyside.png")

        # --- crop of the changed region from both, stacked
        px = a.pad
        box = (max(0, x0 - px), max(0, y0 - px),
               min(im1.width, x1 + px), min(im1.height, y1 + px))
        c1, c2 = im1.crop(box), im2.crop(box)
        crop = Image.new("RGB", (c1.width, c1.height * 2 + 16), "white")
        crop.paste(c1, (0, 0))
        crop.paste(c2, (0, c1.height + 16))
        crop.save(f"{a.out}_crop.png")

        print(f"\nwrote {a.out}_sidebyside.png   old on the left, new on the right")
        print(f"      {a.out}_heat.png          dark where the pixels differ")
        print(f"      {a.out}_crop.png          the changed region only, "
              f"submitted on top, regenerated below")
        print("\nOpen the crop first -- it is the one that answers what changed.")


if __name__ == "__main__":
    main()
