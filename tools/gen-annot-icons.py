#! /usr/bin/env python3
# Generate the icons of the annotation row of the GUI view: the button
# that shows the row (ui_annotate: a square, a circle and a small
# triangle pointing down) and the
# button that opens the editor of the current object (ui_objprops:
# three sliders). White on transparent, 64x64, drawn at 256x256 and
# downsampled, like the rest of dat/assets/icons.
#
# Usage: gen-annot-icons.py [outdir]  (default: dat/assets/icons)

import sys
from PIL import Image, ImageDraw

N = 256  # drawing size
S = 64   # icon size
W = (255, 255, 255, 255)
T = (0, 0, 0, 0)

def new():
    im = Image.new("RGBA", (N, N), T)
    return im, ImageDraw.Draw(im)

def save(im, name, outdir):
    im = im.resize((S, S), Image.LANCZOS)
    # the transparent texels carry the ink color, so that linear
    # filtering does not darken the stroke edges
    px = im.load()
    for j in range(S):
        for i in range(S):
            r, g, b, a = px[i, j]
            px[i, j] = (255, 255, 255, a)
    im.save(outdir + "/" + name)

def annotate(outdir):
    # a square behind a filled circle, the usual "shapes" glyph, with a
    # small triangle pointing down in the corner (the row drops down)
    im, d = new()
    lw = 20
    d.rectangle([12, 76, 140, 204], outline=W, width=lw)
    d.ellipse([62, 0, 202, 140], fill=T)
    d.ellipse([76, 14, 188, 126], fill=W)
    d.polygon([(160, 178), (252, 178), (206, 236)], fill=W)
    save(im, "ui_annotate.png", outdir)

def objprops(outdir):
    im, d = new()
    lw = 16
    knob = 26
    rows = [(60, 170), (128, 80), (196, 140)]
    for y, xk in rows:
        d.line([(24, y), (232, y)], fill=W, width=lw)
        d.ellipse([xk - knob - 10, y - knob - 10, xk + knob + 10, y + knob + 10], fill=T)
        d.ellipse([xk - knob, y - knob, xk + knob, y + knob], fill=W)
    save(im, "ui_objprops.png", outdir)

if __name__ == "__main__":
    outdir = sys.argv[1] if len(sys.argv) > 1 else "dat/assets/icons"
    annotate(outdir)
    objprops(outdir)
