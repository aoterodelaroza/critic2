#!/usr/bin/env python3
"""Draw the icons of the Cell and Slab rows of the GUI builder window.

Usage: python3 tools/gen-builder-cell-icons.py [name ...]

Writes 64x64 white-on-transparent PNGs to dat/assets/icons (the GUI
tints them at draw time): the icons named (cellconv, slaborigin, ...)
or, with no names, all of them. Each icon is drawn at 512x512 and
downsampled. Requires Pillow.
"""
import os
import sys
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


def molecule(dr):
    """Two bonded atoms, in the middle of the canvas."""
    dr.line([(186, 326), (326, 186)], fill=WHITE, width=34)
    dot(dr, 186, 326, 56)
    dot(dr, 326, 186, 56)


def tocrystal():
    """To crystal: a molecule in a cell."""
    img, dr = canvas()
    polyline(dr, [(56, 56), (456, 56), (456, 456), (56, 456)], closed=True, width=36)
    molecule(dr)
    save(img, "ui_tocrystal.png")


def tomolecule():
    """To molecule: a molecule between the corners that are left of a cell."""
    img, dr = canvas()
    for sx, sy in ((1, 1), (-1, 1), (-1, -1), (1, -1)):
        x, y = 256 - 200 * sx, 256 - 200 * sy
        polyline(dr, [(x + 110 * sx, y), (x, y), (x, y + 110 * sy)], width=36)
    molecule(dr)
    save(img, "ui_tomolecule.png")


def slaborigin():
    """Origin: two edges of the cell, with an atom at the corner they share."""
    img, dr = canvas()
    polyline(dr, [(126, 66), (126, 386), (446, 386)])
    dot(dr, 126, 386, 88)
    save(img, "ui_slaborigin.png")


def slabcenter():
    """Center slab: a cell with the slab across its middle."""
    img, dr = canvas()
    polyline(dr, [(106, 56), (406, 56), (406, 456), (106, 456)], closed=True)
    dr.rectangle((106, 206, 406, 306), fill=WHITE)
    save(img, "ui_slabcenter.png")


ICONS = {f.__name__: f for f in (cellconv, cellprim, tocrystal, tomolecule, slaborigin, slabcenter)}

if __name__ == "__main__":
    for name in sys.argv[1:] or ICONS:
        if name not in ICONS:
            sys.exit(f"unknown icon {name!r}; known: {', '.join(ICONS)}")
        ICONS[name]()
