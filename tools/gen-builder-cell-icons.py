#!/usr/bin/env python3
"""Draw the icons of the Cell and Slab rows of the GUI builder window.

Usage: python3 tools/gen-builder-cell-icons.py

Writes 64x64 white-on-transparent PNGs to dat/assets/icons (the GUI
tints them at draw time). Each icon is drawn at 512x512 and
downsampled. Requires Pillow.
"""
import os
from PIL import Image, ImageDraw

BIG = 512  # drawing canvas side (pixels)
OUT = 64  # icon side (pixels)
STROKE = 44  # line width on the drawing canvas
WHITE = (255, 255, 255, 255)
CLEAR = (255, 255, 255, 0)  # transparent white: no dark fringe when filtered
OUTDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "dat", "assets", "icons")


def canvas():
    img = Image.new("RGBA", (BIG, BIG), CLEAR)
    return img, ImageDraw.Draw(img)


def dot(dr, x, y, r):
    dr.ellipse((x - r, y - r, x + r, y + r), fill=WHITE)


def polyline(dr, pts, closed=False, width=STROKE):
    """A line through pts with round joints and ends."""
    if closed:
        pts = list(pts) + [pts[0]]
    dr.line(pts, fill=WHITE, width=width, joint="curve")
    for x, y in pts:
        dot(dr, x, y, width / 2)


def save(img, name):
    img.resize((OUT, OUT), Image.LANCZOS).save(os.path.join(OUTDIR, name))


def cellconv():
    """Conventional cell: a square cell with a lattice point at its center."""
    img, dr = canvas()
    polyline(dr, [(76, 76), (436, 76), (436, 436), (76, 436)], closed=True)
    dot(dr, 256, 256, 58)
    save(img, "ui_cellconv.png")


def cellprim():
    """Primitive cell: a skewed cell with nothing inside."""
    img, dr = canvas()
    polyline(dr, [(196, 96), (456, 96), (316, 416), (56, 416)], closed=True)
    save(img, "ui_cellprim.png")


def slabcenter():
    """Center slab: a cell with the slab across its middle."""
    img, dr = canvas()
    polyline(dr, [(106, 56), (406, 56), (406, 456), (106, 456)], closed=True)
    dr.rectangle((106, 206, 406, 306), fill=WHITE)
    save(img, "ui_slabcenter.png")


if __name__ == "__main__":
    cellconv()
    cellprim()
    slabcenter()
