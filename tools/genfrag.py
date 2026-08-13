#!/usr/bin/env python
"""Draw the fragment diagrams the GUI builder shows for its library.

Every diagram in dat/assets/fragments comes from here, so that they share one
style: same font, same line weight, same colours, and the same red disc for
the attachment point. Run it after adding or changing a fragment in
dat/lib/fragment.dat, and commit the images it writes.

  Usage:  genfrag.py [outdir [name ...]]

with outdir defaulting to dat/assets/fragments and, with no names, every
fragment in the library. Needs RDKit and Pillow, neither of which critic2
itself uses; in a conda environment that has them,

  conda run -n <env> python tools/genfrag.py

The structures are the ones in dat/lib/fragment.dat: 3D geometry with explicit
hydrogens, an anchor atom and an attach placeholder. RDKit perceives the bond
orders and lays the depiction out; the attachment point is the end of the
attach bond for a substituent, and the metal site for a ligand (the centre of
the donors, stepped aside if that would land on the skeleton).

The ligands are stored as they coordinate, with the donor hydrogens removed,
so most of them only perceive as anions. That charge belongs to the free
ligand and not to the picture of it bound to a metal, so it is dropped from
the donors before drawing -- after kekulization, which needs it.
"""
import io
import math
import os
import sys

from PIL import Image, ImageDraw
from rdkit import Chem, RDLogger
from rdkit.Chem import rdDepictor, rdDetermineBonds
from rdkit.Chem.Draw import rdMolDraw2D

RDLogger.DisableLog("rdApp.*")

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DAT = os.path.join(ROOT, "dat", "lib", "fragment.dat")
OUT = os.path.join(ROOT, "dat", "assets", "fragments")
SIDE = 192          # final png side, px
SUPER = 3           # render this many times larger and downsample
DOTFRAC = 0.075     # red disc radius, as a fraction of the image side
DOTCOL = (206, 42, 42)
DONORTOL = 1.25     # a ligand donor is this much farther than the nearest one
DONORCUT = 3.2      # ... and no farther than this from the metal site, angstrom

# darker than the RDKit default for the elements that wash out on white
PALETTE = {9: (0.15, 0.55, 0.20), 17: (0.05, 0.55, 0.05), 16: (0.70, 0.55, 0.00),
           34: (0.60, 0.40, 0.10), 53: (0.45, 0.15, 0.60), 5: (0.60, 0.35, 0.35),
           14: (0.40, 0.40, 0.45), 15: (0.85, 0.45, 0.00)}

# charges to try, in order; the ligands need the negative ones
CHARGES = [0, -1, -2, -3, -4, 1, -5, -6]


def parse(path):
    """The fragments in fragment.dat: name, atoms, coords, anchor, attach."""
    frags, cur, inmol = [], None, False
    for line in open(path):
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        w = s.split()
        if w[0] == "structure":
            cur = dict(name=w[1], sym=[], xyz=[], anchor=0, attach=0, radius=None)
        elif w[0] == "molecule":
            inmol = True
        elif w[0] == "endmolecule":
            inmol = False
        elif w[0] == "anchor":
            cur["anchor"] = int(w[1])
        elif w[0] == "attach":
            cur["attach"] = int(w[1])
        elif w[0] == "radius":
            cur["radius"] = float(w[1])
        elif w[0] == "endstructure":
            frags.append(cur)
            cur = None
        elif inmol and cur is not None:
            cur["sym"].append(w[0])
            cur["xyz"].append([float(x) for x in w[1:4]])
    return frags


def build(frag):
    """Perceive the structure. Returns the mol, the index of the attachment
    hydrogen (substituent) or None, and the indices of the donor atoms
    (ligand), all in the mol's own numbering."""
    sym, xyz = frag["sym"], frag["xyz"]
    iat = frag["attach"] - 1
    isligand = sym[iat].lower() == "bq"

    keep = [i for i in range(len(sym)) if sym[i].lower() != "bq"]
    lines = ["%d" % len(keep), frag["name"]]
    for i in keep:
        lines.append("%-3s %14.8f %14.8f %14.8f" % (sym[i], *xyz[i]))
    raw = Chem.MolFromXYZBlock("\n".join(lines) + "\n")
    if raw is None:
        raise RuntimeError("xyz parse failed")

    err = None
    for q in CHARGES:
        mol = Chem.Mol(raw)
        try:
            rdDetermineBonds.DetermineBonds(mol, charge=q)
            Chem.SanitizeMol(mol)
            break
        except Exception as e:                                  # noqa: BLE001
            err, mol = e, None
    if mol is None:
        raise RuntimeError("bond perception failed: %s" % err)

    # mark the atoms the drawing needs, so RemoveHs cannot renumber them away
    if isligand:
        m = xyz[iat]
        dist = [(math.dist(xyz[i], m), k) for k, i in enumerate(keep)
                if sym[i].lower() != "h"]
        dmin = min(d for d, _ in dist)
        cut = min(max(DONORTOL * dmin, dmin + 0.3), DONORCUT)
        near = [k for d, k in dist if d <= cut]
        # a heteroatom donor beats a carbon at the same distance; only the
        # haptic ligands, where there is no heteroatom at all, bind by carbon
        het = [k for k in near
               if mol.GetAtomWithIdx(k).GetSymbol() in ("N", "O", "P", "S", "Se", "As")]
        near = het or near
        ndon = 0
        for k in near:
            mol.GetAtomWithIdx(k).SetBoolProp("donor", True)
            ndon += 1
        if ndon == 0:
            raise RuntimeError("no donor near the metal site")
        return mol, True, q
    mol.GetAtomWithIdx(keep.index(iat)).SetBoolProp("attach", True)
    return mol, False, q


def tagged(mol, prop):
    return [a.GetIdx() for a in mol.GetAtoms() if a.HasProp(prop)]


def draw(mol, isligand, qtot, path):
    """Render the fragment and stamp the red attachment disc on it."""
    mol = Chem.RWMol(mol)
    Chem.RemoveStereochemistry(mol)
    if isligand:
        # the donors' free valences are taken by the metal, so the charge
        # the free ligand needed is not part of this picture. Kekulize
        # first: an aromatic anion cannot be kekulized once it is neutral
        Chem.Kekulize(mol, clearAromaticFlags=True)
        don = tagged(mol, "donor")
        drop = []
        for i in don:
            a = mol.GetAtomWithIdx(i)
            a.SetFormalCharge(0)
            a.SetNumRadicalElectrons(0)
            # the hydrogens this donor really has go into its label; they
            # are taken out here because RemoveHs will not touch an atom
            # that carries NoImplicit, and would count them twice if it did
            nb = [x.GetIdx() for x in a.GetNeighbors() if x.GetAtomicNum() == 1]
            a.SetNumExplicitHs(len(nb))
            a.SetNoImplicit(True)
            drop += nb
        for i in sorted(drop, reverse=True):
            mol.RemoveAtom(i)
        don = tagged(mol, "donor")
        # a carboxylate that binds through one of its oxygens is drawn
        # C(=O)-O-[M], whichever way round the perception put the double
        # bond: without this the two halves of oxalate come out different
        for c in mol.GetAtoms():
            if c.GetSymbol() != "C":
                continue
            ox = [nb for nb in c.GetNeighbors()
                  if nb.GetSymbol() == "O" and nb.GetDegree() == 1]
            odon = [o for o in ox if o.HasProp("donor")]
            ofree = [o for o in ox if not o.HasProp("donor")]
            if len(odon) == 1 and len(ofree) == 1:
                mol.GetBondBetweenAtoms(c.GetIdx(), odon[0].GetIdx()).SetBondType(
                    Chem.BondType.SINGLE)
                mol.GetBondBetweenAtoms(c.GetIdx(), ofree[0].GetIdx()).SetBondType(
                    Chem.BondType.DOUBLE)
                for o in (odon[0], ofree[0]):
                    o.SetFormalCharge(0)
                    o.SetNumRadicalElectrons(0)
                    o.SetNoImplicit(True)

        # the rest of the ligand is drawn as the free ligand would be: no
        # leftover halves of charge-separated pairs, and the oxygens and
        # nitrogens that are not bound to the metal take their hydrogen
        for a in mol.GetAtoms():
            if a.GetFormalCharge() != 0 and not a.HasProp("donor"):
                a.SetFormalCharge(0)
                a.SetNumRadicalElectrons(0)
                a.SetNoImplicit(not (a.GetSymbol() == "O" and a.GetDegree() == 1))
        if len(don) == 1:
            # one donor: the metal would sit on top of its label, so give
            # it a bond of its own, as a substituent has
            idum = mol.AddAtom(Chem.Atom(0))
            mol.AddBond(don[0], idum, Chem.BondType.SINGLE)
            mol.GetAtomWithIdx(idum).SetBoolProp("attach", True)
            mol.GetAtomWithIdx(don[0]).ClearProp("donor")
    else:
        # the attachment hydrogen becomes the point the disc marks
        a = mol.GetAtomWithIdx(tagged(mol, "attach")[0])
        a.SetAtomicNum(0)
        a.SetNoImplicit(True)
        if qtot < 0:
            # stored without its acidic hydrogen: that charge belongs to
            # the fragment as stored, not to the group being drawn, so the
            # oxygen takes its hydrogen back. Groups that are charged in
            # their own right (nitro, azido, an ammonium) keep theirs
            for b in mol.GetAtoms():
                if b.GetFormalCharge() < 0 and b.GetSymbol() == "O" \
                   and b.GetDegree() == 1:
                    b.SetFormalCharge(0)
                    b.SetNoImplicit(False)
    for a in mol.GetAtoms():
        if a.HasProp("attach"):
            a.SetProp("atomLabel", " ")

    mol = Chem.RemoveHs(mol.GetMol(), sanitize=False)
    mol.RemoveAllConformers()
    rdDepictor.SetPreferCoordGen(True)
    rdDepictor.Compute2DCoords(mol)

    px = SIDE * SUPER
    d = rdMolDraw2D.MolDraw2DCairo(px, px)
    o = d.drawOptions()
    o.dummiesAreAttachments = False
    o.addStereoAnnotation = False
    o.clearBackground = True
    o.padding = 0.10
    o.bondLineWidth = 2 * SUPER
    o.minFontSize = 13 * SUPER
    o.maxFontSize = 16 * SUPER
    o.updateAtomPalette(PALETTE)
    d.DrawMolecule(rdMolDraw2D.PrepareMolForDrawing(mol, not isligand, False))
    d.FinishDrawing()

    im = Image.open(io.BytesIO(d.GetDrawingText())).convert("RGB")
    mark = tagged(mol, "attach") or tagged(mol, "donor")
    pts = [d.GetDrawCoords(i) for i in mark]
    cx = sum(p.x for p in pts) / len(pts)
    cy = sum(p.y for p in pts) / len(pts)
    if len(mark) > 1:
        # the centroid of the donors is where the metal belongs, but in a
        # chelate it can land on the skeleton: step away from it until the
        # disc clears the atoms, which a macrocycle never needs
        pos = [d.GetDrawCoords(i) for i in range(mol.GetNumAtoms())]
        bl = [math.dist((pos[b.GetBeginAtomIdx()].x, pos[b.GetBeginAtomIdx()].y),
                        (pos[b.GetEndAtomIdx()].x, pos[b.GetEndAtomIdx()].y))
              for b in mol.GetBonds()]
        bl = sorted(bl)[len(bl) // 2] if bl else px / 8.
        # a labelled atom takes more room than a plain vertex
        rad = [(0.50 if a.GetSymbol() != "C" else 0.28) * bl for a in mol.GetAtoms()]
        clear = DOTFRAC * px
        near = lambda x, y: min(math.dist((x, y), (q.x, q.y)) - rad[i]
                                for i, q in enumerate(pos))
        if near(cx, cy) < clear:
            best = None
            for step in (0.7, 1.1, 1.5):
                for k in range(24):
                    t = k * math.pi / 12.
                    x = cx + step * bl * math.cos(t)
                    y = cy + step * bl * math.sin(t)
                    if not (clear < x < px - clear and clear < y < px - clear):
                        continue
                    cand = (near(x, y), -step, x, y)
                    if best is None or cand > best:
                        best = cand
                if best is not None and best[0] >= clear:
                    break
            if best is not None:
                cx, cy = best[2], best[3]
    r = DOTFRAC * px
    ImageDraw.Draw(im).ellipse([cx - r, cy - r, cx + r, cy + r], fill=DOTCOL)
    im.resize((SIDE, SIDE), Image.LANCZOS).save(path)


if __name__ == "__main__":
    out = sys.argv[1] if len(sys.argv) > 1 else OUT
    only = sys.argv[2:]
    os.makedirs(out, exist_ok=True)
    nok, bad = 0, []
    for f in parse(DAT):
        if only and f["name"] not in only:
            continue
        try:
            mol, isligand, q = build(f)
            draw(mol, isligand, q, os.path.join(out, f["name"] + ".png"))
            nok += 1
        except Exception as e:                                  # noqa: BLE001
            bad.append((f["name"], str(e).replace("\n", " ")[:80]))
    print("drawn %d, failed %d" % (nok, len(bad)))
    for n, e in bad:
        print("  FAIL %-24s %s" % (n, e))
