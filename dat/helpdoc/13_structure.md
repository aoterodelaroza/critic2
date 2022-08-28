
## Relabel the Atoms in the Structure (ATOMLABEL) {#c2-atomlabel}

The ATOMLABEL keyword can be used to change the atomic labels for the
atoms in the current structure:
~~~
ATOMLABEL template.s
~~~
The template string (`template.s`) is used to build the new atomic
names. The format specifiers for this template are:

* `%aid`: the index for the atom in the non-equivalent atom list.

* `%id`: the index for the atom in the non-equivalent atom list,
   counting only the atoms of the same type.

* `%S`: the atomic symbol, derived from the current atomic number.

* `%s`: same as `%S`, but lowercase.

* `%l`: the current atom label.

## Powder Diffraction (POWDER) {#c2-powder}

The keyword POWDER calculates the powder diffraction pattern for the
current crystal structure:
~~~
POWDER [TH2INI t2i.r] [TH2END t2e.r] [{L|LAMBDA} l.r]
       [FPOL fpol.r] [NPTS npts.i] [SIGMA sigma.r]
       [ROOT root.s] [HARD|SOFT]
~~~
Only the $$2\theta$$ range from `t2i.r` (default: 5 degrees) to
`t2e.r` (default: 90 degrees) will be plotted. The wavelength of the
incident radiation is given by `l.r` (in angstrom, default: 1.5406
angstrom, corresponding to Cu-K&alpha;). The polarization of the X-ray
radiation affects the treatment of the resulting intensities. The
default is `fpol.r = 0`, unpolarized light. For synchrotron radiation,
use `fpol.r = 0.95`. `npts.i` is the number of points in the generated
spectrum (default: 10001). Gaussian broadening is used on the observed
peaks, with width parameter `sigma.r` (default: 0.05 degrees).

By default, two files are generated: `<root>_xrd.dat`, containing the
$$2\theta$$ versus intensity data, and `<root>_xrd.gnu`, the gnuplot
script to plot it. The name of these files can be changed using the
ROOT keyword. The Miller indices of the peaks are written to the
standard output.

The `HARD` and `SOFT` keywords control how peaks right outside the
plot range are treated. In a `HARD` powder diffraction pattern, the
peaks outside the plot range are not computed, and therefore their
tails will not appear on the plot even if the corresponding peak is so
close that the contribution to the RDF would be significant. In a
`SOFT` powder diffraction pattern, all peaks are computed and
represented, even if the maximum is outside the plot range. Default:
`SOFT`.

## Radial Distribution Function (RDF) {#c2-rdf}

The RDF keyword calculates the radial distribution function (RDF) for
the current crystal or molecular structure:
~~~
RDF [RINI t2i.r] [REND t2e.r] [SIGMA sigma.r] [NPTS npts.i]
    [ROOT root.s] [PAIR is1.s is2.s [PAIR is1.s is2.s ...]]
    [HARD|SOFT]
~~~
The definition of RDF is similar to the one found in
[Willighagen et al., Acta Cryst. B 61 (2005) 29](https://doi.org/10.1107/S0108768104028344),
but where the atomic charges are replaced by the square root of the
atomic number. Each pair of atoms A and B at a distance $$d_{AB}$$ in
the system will create a peak in the RDF with a Gaussian shape:
\begin{equation}
I_{AB}(r) = \sqrt{Z_AZ_B}\times \exp\left(-\frac{(r-d_{AB})^2}{2\sigma^2}\right)
\end{equation}
The RDF is constructed from the sum of all $$I_{AB}(r)$$ functions from all
unique atom pairs by doing:
\begin{equation}
{\rm RDF}(r) = \sum_{A<B} \frac{I_{AB}(r)}{r^2N_{\rm cell}}
\end{equation}
where $$N_{\rm cell}$$ is the number of atoms in the unit cell.

The RDF is plotted from an initial distance `t21.r` (default: 0) up to
a maximum distance `t2e.r` (default: 25 bohr) using `npts.i` points in
that interval (default: 10001). Gaussian broadening is used with sigma
equal to `sigma.r` (default: 0.05 bohr). The default units of `RINI`,
`REND`, and `SIGMA` are bohr in crystals and angstrom in molecules
(this can be changed with the
[UNITS](/critic2/manual/inputoutput/#c2-units) keyword).

Two files are generated: `<root>_rdf.dat`, containing the RDF versus
distance data, and `<root>_rdf.gnu`, the gnuplot script to plot
it. The name of these files can be changed using the ROOT keyword. If
`PAIR` is given, only the distances between atoms of type `is1.s` and
`is2.s` will contribute to the RDF. Multiple `PAIR` keywords can be
used. The `is1.s` and `is2.s` must be the name of an atomic species in
the system.

The `HARD` and `SOFT` keywords control how peaks right outside the
plot range are treated. In a `HARD` RDF, the peaks outside the plot
range are not computed, and therefore their tails will not appear on
the plot even if the corresponding peak is so close that the
contribution to the RDF would be significant. In a `SOFT` RDF, all
peaks are computed and represented, even if the maximum is outside the
plot range. Default: `SOFT`.

## Compare Crystal and Molecular Structures (COMPARE) {#c2-compare}

The COMPARE keyword compares two or more structures:
~~~
COMPARE {.|file1.s} {.|file2.s} [{.|file3.s} ...]
COMPARE ... [MOLECULE|CRYSTAL] [REDUCE eps.r]
COMPARE ... [POWDER|RDF] [XEND xend.r] [SIGMA sigma.r]  ## crystals
COMPARE ... [SORTED|RDF|ULLMANN|UMEYAMA]  ## molecules
~~~
At least two structures are required for the comparison.
The structures can be given as external files (`file1.s`,
`file2.s`,...). The behavior is the same as in CRYSTAL and MOLECULE:
the file format is identified using the file extension or its contents
if the extension is not enough. If a dot (".") is used instead of a
file name, the current structure (previously loaded with
CRYSTAL/MOLECULE) is used.

There are two distinct modes of operation in COMPARE, depending on
whether a molecular or crystal comparison is carried out. If the
structures are all molecules or if the MOLECULE keyword is used, then
the structures are compared as molecules. If any one of the structures
is a crystal or if the CRYSTAL keyword is used, a crystal comparison
is carried out.

There are two ways of calculating a comparison between crystals, based
on the radial distribution functions (RDF keyword) or the powder
diffraction patterns (POWDER keyword). The default is POWDER.
In both cases, COMPARE finds the measure of similarity (DIFF) based on
the corresponding functions (RDF or diffractogram).
Two crystal structures are exactly equal if `DIFF = 0`. Maximum
dissimilarity occurs when `DIFF = 1`.  The crystal similarity
measure is calculated using the cross-correlation functions defined in
[de Gelder et al., J. Comput. Chem., 22 (2001) 273](https://doi.org/10.1002/1096-987X(200102)22:3%3C273::AID-JCC1001%3E3.0.CO;2-0),
with the triangle weight. Powder diffraction patterns are calculated
from $$2\theta = 5$$ up to `xend.r` (XEND keyword, default:
50). Radial distribution functions are calculated from zero up to
`xend.r` bohr (XEND keyword, default: 25 bohr). SIGMA is the Gaussian
broadening parameter for the powder diffraction or RDF peaks.

For the molecular comparison, there are several options. If the SORTED
keyword is used, the atomic sequence in each molecule is assumed to be
the same. In this case, COMPARE finds the translation and rotation
that brings the molecules into closest agreement with each other and
reports the resulting root-mean-square (RMS) of the atomic positions.
The molecular rotation is calculated using Walker et al.'s quaternion
algorithm
([Walker et al., CVGIP-Imag. Understan. 54 (1991) 358](https://doi.org/10.1016/1049-9660(91)90036-O)).
For the comparison to work correctly, it is necessary that the two
molecules have the same number of atoms and that the atoms are in the
same sequence.

If the molecules have the same atomic connectivity (the same molecular
diagram) but the atomic sequences are not equivalent, i.e. the atoms
are disordered, then the ULLMANN or UMEYAMA methods can be used.  The
Umeyama approach (keyword: UMEYAMA) uses a weighted graph
matching method based on the algorithm proposed in
[Umeyama, S., IEEE PAMI, 10 (1988) 695-703](http://dx.doi.org/10.1109/34.6778).
Umeyama's method is non-iterative, but approximate. If there are
significant differences between the structures being compared, or if
they are highly symmetric, UMEYAMA may find the wrong atomic
permutation, so it is a good idea to always check that the RMS in the
output is reasonable.  The Ullmann approach (keyword: ULLMANN) uses a
modified version of Ullmann's subgraph matching algorithm
([Ullmann, J. R., J. ACM 23 (1976) 31-42](https://doi.org/10.1145/321921.321925)).
This method finds all possible matching sequences based on graph
connectivity alone, then selects the sequence with lowest RMS upon
molecular rotation. It is more reliable than UMEYAMA but may become
expensive for larger molecules.

The ULLMANN method is the default if the number and types of atoms in
the molecules being compared are the same. If this is not the case, or
if the RDF keyword is used, then radial distribution functions are
employed and the comparison is similar to how RDF works in crystals.

Lastly, if more than two structures are used in COMPARE, critic2 will
compare each pair of structures and present the resulting similarity
matrix. Alternatively, if the REDUCE keyword is used, then a threshold
(`eps.r`) is applied to determine whether two structures are equal or
not. Critic2 then prints a list of unique structures and repeated
structures in the output.

The COMPARE keyword does not require a previous CRYSTAL or MOLECULE
keyword. Hence, valid critic2 inputs would be:
~~~
COMPARE bleh1.scf.in bleh2.cif
COMPARE bleh1.xyz bleh2.wfx
~~~
provided the files exist.

## Transform the Unit Cell (NEWCELL) {#c2-newcell}

The NEWCELL keyword transforms the unit cell used to describe the
current crystal structure to a new cell:
~~~
NEWCELL {x1.r y1.r z1.r x2.r y2.r z2.r x3.r y3.r z3.r|n1.i n2.i n3.i} [INV|INVERSE]
        [ORIGIN x0.r y0.r z0.r]
NEWCELL [{PRIMSTD|STANDARD|PRIMITIVE} [REFINE]]
NEWCELL [NIGGLI|DELAUNAY]
~~~
The new unit cell is given by the vectors (`x1.r` `y1.r` `z1.r`),
(`x2.r` `y2.r` `z2.r`), and (`x3.r` `y3.r` `z3.r`) in
crystallographic coordinates relative to the old unit cell. The x1,
x2, x3 vectors must be pure translations of the old cell; either
lattice vectors, centering vectors, or combinations of the
two. Alternatively, if three integers are given (`n1.i` `n2.i` `n3.i`),
NEWCELL builds a supercell with `n1.i` cells in the a direction,
`n2.i` cells in the b direction, and `n3.i` cells in the c direction.

NEWCELL unloads all fields (except the promolecular density) and
clears the critical point list. If the INV (or INVERSE) keyword is
used, the input vectors correspond to the crystallographic coordinates
of the old cell in the new coordinate system. A NEWCELL transformation
is the inverse of the same transformation using the INV keyword.
Optionally, if an ORIGIN vector is given, (`x0.r` `y0.r` `z0.r`), the
cell origin is translated to x0. The units of x0 are crystallographic
coordinates of the original cell.

The NEWCELL keyword is useful for building supercells or for
performing routine but tedious crystallographic transformations. For
instance, given a face-centered cubic lattice and the conventional
cubic cell one can find the primitive (rhombohedral) cell by doing:
~~~
CRYSTAL LIBRARY mgo
NEWCELL 1/2 1/2 0 1/2 0 1/2 0 1/2 1/2
~~~
Likewise, if the current cell is rhombohedral, the same NEWCELL order
but including the INVERSE keyword transforms to the cubic. That is:
~~~
CRYSTAL LIBRARY mgo
NEWCELL 1/2 1/2 0 1/2 0 1/2 0 1/2 1/2
NEWCELL 1/2 1/2 0 1/2 0 1/2 0 1/2 1/2 INVERSE
~~~
gives a unit cell and crystal structure description that is equivalent
to the initial one read from the library.

NEWCELL also admits specific keywords that perform common
transformations to certain cells of interest. The cell can be
transformed to:

* STANDARD: standard (canonical) unit cell.

* PRIMITIVE: standard primitive unit cell. Does not transform the cell
  if the unit cell is already primitive.

* PRIMSTD: standard primitive unit cell. Does the transformation even
  if the current unit cell is primitive.

* NIGGLI: Niggli-reduced cell for the current lattice. Use a NEWCELL
  PRIMITIVE first to get the primitive Niggli cell.

* DELAUNAY: Delaunay-reduced cell for the current lattice. Use a
  NEWCELL PRIMITIVE first to get the primitive Delaunay cell.

The origin is not translated by any of these keywords. The REFINE
keyword can be used in combination with STANDARD, PRIMITIVE, or
PRIMSTD. If REFINE is used, then the atomic positions are idealized
according to the space group symmetry operations. For instance, an
atomic coordinate of 0.333 may become 0.3333333... in a hexagonal
space group. These transformations use the
[spglib](https://atztogo.github.io/spglib/) library. Please consult
the spglib manual for more information.

## Calculate Atomic Environments (ENVIRON) {#c2-environ}

The ENVIRON keyword prints lists of neighbor atoms:
~~~
ENVIRON [DIST dist.r] [POINT x0.r y0.r z0.r|ATOM at.s/iat.i|CELATOM iat.i]
[BY by.s/iby.i] [SHELL|SHELLS]
~~~
If POINT is given, print the atomic neighbors around the point with
coordinates (`x0.r` `y0.r` `z0.r`) in crystallographic coordinates
(crystal) or or molecular Cartesian coordinates (molecule, default
units: angstrom). Instead, if ATOM is given, print the neighbors
around atom `iat.i` from the non-equivalent atom list or around every
atom with atomic symbol `at.s` (converted internally to atomic
number). If CELATOM is used, then print the environment around atom
iat.i from the complete list.  If neither POINT nor ATOM nor CELATOM
are given, print the environments of all non-equivalent atoms in the
unit cell.

By default, the environments extend up to 5 angstrom from the central
point. The DIST keyword can be used to change this value (by default,
`dist.r` is in bohr in crystals and angstrom in molecules). The BY
keyword allows filtering the neighbor list to print only certain kinds
of atoms. If `iby.i` is given, print only atoms whose non-equivalent ID
is the same as `iby.i`. If `by.s` is given, print only atoms with the same
atomic symbol as `by.s` (converted internally to atomic number). If
SHELLS (or SHELL) is given, group the neighbors in shells by distance
(1e-2 atomic distance threshold for atoms in the same shell) and
non-equivalent ID.

## Calculate Coordination Polyhedra (POLYHEDRA) {#c2-polyhedra}

The POLYHEDRA keyword calculates the coordination polyhedra in a
periodic solid:
~~~
POLYHEDRA atcenter.s atvertex.s [[rmin.r] rmax.r]
~~~
The polyhedra are built with atom `atcenter.s` as the center and
`atvertex.s` as the vertices. These strings can refer to an atomic
species, if such species exists in the current structure. If not, then
the symbols are converted to atomic numbers, and these are used
instead. By default, the distance range to consider the central and
vertex atom are coordinated is from zero to the sum of covalent radii
times the [BONDFACTOR](/critic2/manual/misc/#c2-bondfactor). If a
single number (`rmax.r`) is indicated at the end of the POLYHEDRA
command, then the distance range goes from zero to `rmax.r`. If two
numbers are indicated, then the distance range goes from `rmax.r` to
`rmin.r`.

The output of the POLYHEDRA keyword contains a list of non-equivalent
atoms of type `atcenter.s` that have a coordination polyhedron. The
number of vertices, minimum and maximum vertex-center distance, number
of faces, and polyhedron volume are printed.

## Effective coordination number (ECON) {#c2-econ}

The coordination number of an atom is typically defined as the number
of neighbors that are closest to that atom. This definition can be
unsatisfactory for situations where there is a range of bond lengths
around the central atom. To address these complex cases, the
calculation of an effective coordination number (ECoN) was introduced
by Hoppe in 1979. The ECoN of a given atom is calculated by assigning
to each atom around it a weight based on their distance. The original
procedure is described in
[Hoppe, Z. Kristallogr. 150 (1979) 23](http://dx.doi.org/10.1524/zkri.1979.150.14.23)
and examined in detail in
[Nespolo, Acta Cryst. B, 72 (2016) 51](http://dx.doi.org/10.1107/S2052520615019472).

The implementation of ECoN in critic2 is slightly different from those
two references in that we do not require the user to define
coordination polyhedra for the calculation. Instead, the ECoN in
critic2 is calculated as a formula that takes into account the
distances from the central atom to all other atoms in the
crystal. There are two variants of ECoN: iterative, calculated using a
self-consistent procedure, and non-iterative.

For a given atom, the ECoN is defined as:
\begin{equation}
{\rm ECoN} = \sum_{i} w_{i}
\end{equation}
where the sum runs over all the atoms in the environment of the chosen
central atom (or a subset of atoms belonging to a certain species, see
below). The weight of the ith atom ($$w_{i}$$) is defined as:
\begin{equation}
w_{i} = \exp\left[1-\left(\frac{d_{i}}{d_{\rm av}}\right)^6\right]
\end{equation}
with $$d_{i}$$ the distance to atom $$i$$. In the non-iterative
variant of ECoN, $$d_{\rm av}$$ is the weighted average distance,
defined as:
\begin{equation}
d_{\rm av} = \frac{\sum_{i} d_{i}\exp\left[1-\left(\frac{d_{i}}{d_{\rm min}}\right)^6\right]}
{\sum_{i}\exp\left[1-\left(\frac{d_{i}}{d_{\rm min}}\right)^6\right]}
\end{equation}
where $$d_{\rm min}$$ is the shortest distance to the central atom for
the considered atomic species. In contrast, the iterative variant of
ECoN calculates the weighted average distance by solving the
non-linear equation:
\begin{equation}
d_{\rm av}^{\rm it} = \frac{\sum_{i} d_{i}\exp\left[1-\left(\frac{d_{i}}{d_{\rm av}^{\rm it}}\right)^6\right]}
{\sum_{i}\exp\left[1-\left(\frac{d_{i}}{d_{\rm av}^{\rm it}}\right)^6\right]}
\end{equation}
which is done using a self-consistent iterative procedure.

A typical output for the ECoN keyword is:
~~~
# nid->spc   name(nid)->name(spc)     econ      1econ       nd         1nd
   1 -> *          Ti -> *           5.9754     5.9743     3.6805     3.6804
   1 -> 1          Ti -> Ti          4.8730     4.7197     5.8587     5.8302
   1 -> 2          Ti -> O           5.9754     5.9742     3.6805     3.6804
   2 -> *           O -> *           3.2704     3.2170     3.7153     3.7056
   2 -> 1           O -> Ti          2.9877     2.9871     3.6805     3.6804
   2 -> 2           O -> O           7.9476     6.4604     5.1286     4.9802
~~~
where `nid` represents the non-equivalent atom ID of the central atom
for which the ECoN is calculated and the values corresponding to
different `spc` are obtained by considering the distances to atoms
belonging to those atomic species only (`spc=1`, `2`, etc.) or to all
atoms regardless of species (`spc=*`). Critic2 reports the
iterative ECoN (`econ`), non-iterative ECoN (`1econ`), iterative
weighted average distance (`nd`, $$d_{\rm av}^{\rm it}$$), and
non-iterative weighted average distance (`1nd`, $$d_{\rm av}$$).

## Pair and Triplet Coordination Numbers (COORD) {#c2-coord}

The COORD keyword calculates pair and triplet coordination numbers:
~~~
COORD [DIST dist.r] [FAC fac.r]
~~~
By default two atoms are coordinated if they are within `fac.r` times
the sum of their covalent radii. By default, `fac.r` is equal to the
[BONDFACTOR](/critic2/manual/misc/#c2-bondfactor) and the internal
list of covalent radii is used (see the
[RADII](/critic2/manual/misc/#c2-radii) keyword). The value of `fac.r`
can be changed with the FAC keyword. The atomic radii for atomic
species can be changed with using the [RADII keyword](/critic2/manual/misc/#c2-radii).
If the DIST keyword is used, all atoms within a distance `dist.r` are
coordinated.

On output, COORD will list the number of coordinated pairs per atom in
the unit cell and per atomic species. In addition, it will also list
all coordinated triplets X-Y-Z, where Y runs over all atoms in the
unit cell and over all atomic species.

## Packing Ratio (PACKING) {#c2-packing}
The PACKING keyword computes the packing ratio of the crystal.
~~~
PACKING {COV|VDW|} [PREC prec.r]
~~~
With VDW, use the van der Waals radii. With COV, use the covalent
radii. If neither VDW nor COV are used, use half of the nearest
neighbor distance (in this case, the spheres would not overlap).

In the VDW and COV cases, the calculation is done using a Monte-Carlo
sampling of the unit cell. The PREC keyword allows controlling the
precision of this calculation. PREC corresponds to the standard
deviation in the van der Waals volume divided by the volume
itself. The default `prec.r` is 0.01. The van der Waals and covalent
radii can be changed using the [RADII](/critic2/manual/misc/#c2-radii)
keyword.

## Van der Waals volume (VDW) {#c2-vdw}
The VDW keyword calculates the van der Waals volume of a crystal or
molecule.
~~~
VDW [PREC prec.r]
~~~
The calculation is done using a Monte-Carlo sampling of the unit cell
or the molecular volume. The PREC keyword allows controlling the
precision of this calculation. PREC corresponds to the standard
deviation in the van der Waals volume divided by the volume
itself. The default `prec.r` is 0.01. The van der Waals radii are
taken from the internal tables, and they can be changed using the
[RADII](/critic2/manual/misc/#c2-radii) keyword.

## Identify Atoms in the Structure Given Their Coordinates (IDENTIFY) {#c2-identify}

The IDENTIFY keyword identifies the coordinates in the user input by
matching them against the internal list of atoms and critical points:
~~~
IDENTIFY [ANG|ANGSTROM|BOHR|AU|CRYST]
 x.r y.r z.r [ANG|ANGSTROM|BOHR|AU|CRYST]
 ...
 file.xyz
ENDIDENTIFY/END
IDENTIFY file.xyz
~~~
If a coordinate is close (1e-4 bohr) to an atom or CP, the
corresponding indices as well as the coordinates are written to the
output. The input can come as either the coordinates of the points
themselves or a filename pointing to an `.xyz` file. IDENTIFY can be
used in environment mode (IDENTIFY/ENDIDENTIFY, with several lines of
input) or as a single command when applying it to a single `.xyz`
file.

The default units are crystallographic coordinates in crystals and
molecular Cartesian coordinates in molecules (default unit:
angstrom). However, these units can be modified with one of the
keywords that follow the IDENTIFY command. For specific points, the
unit can be changed by specifying the same keywords after the three
coordinates. The units in the `.xyz` file are angstrom (the xyz file
also has to have the usual syntax, with the number of atoms in the
first line and the title, or a blank line, in the second line).

In addition to the identity of the points, if any, critic2 also
provides the vertices of the cube that encompasses all the points in
the list that did match an atom or CP.

## Point-Charge Electrostatic Energy (EWALD) {#c2-ewald}

The EWALD keyword calculates the electrostatic energy of the lattice
of point charges using Ewald's method:
~~~
EWALD
~~~
The atomic charges are defined using the Q keyword.

## Reorder Atoms in a Molecule or Molecular Crystal (MOLREORDER) {#c2-molreorder}

The MOLREORDER keyword reorders the atoms in a target molecule or all
the molecules in a target molecular crystal to have the same atomic
sequence as in a template molecule. The syntax is:
~~~
MOLREORDER {.|template.s} {.|target.s} [WRITE file.s] [MOVEATOMS] [INV] [UMEYAMA|ULLMANN]
~~~
The template molecule is in file `template.s`. The target structure
file `target.s` must contain either a molecule or a molecular
crystal. If `target.s` is a molecule, it must contain the same number
and types of atoms as the template. If `target.s` is a molecular
crystal, its asymmetric unit must have an integer number of molecules
(Z', see the WHOLEMOLS option in
[SYMM/SYM](/critic2/manual/crystal/#c2-symm)). In addition, all
molecular fragments in the target crystal must have the same number
and atom types as the template. On output, the atoms in structure
`target.s` are reordered such that they are in the same sequence as
`template.s`. A dot (".") given as either the template (`template.s`)
or target (`target.s`) files instructs critic2 to use the currently
loaded structure (with a previous CRYSTAL/MOLECULE command).

MOLREORDER can use two different algorithms to solve the atomic
sequence reordering problem:

* The Umeyama approach (keyword: UMEYAMA) uses a weighted graph
  matching method based on the algorithm proposed in
  [Umeyama, S., IEEE PAMI, 10 (1988) 695-703](http://dx.doi.org/10.1109/34.6778).
  Umeyama's method is non-iterative, but approximate. If there are
  significant differences between the structures being compared, or if
  they are highly symmetric, UMEYAMA may find the wrong atomic
  permutation, so it is a good idea to always check that the RMS in the
  output is small.

* The Ullmann approach (keyword: ULLMANN) uses a modified version of
  Ullmann's subgraph matching algorithm ([Ullmann, J. R., J. ACM 23 (1976) 31-42](https://doi.org/10.1145/321921.321925)).
  This method finds all possible matching sequences based on graph
  connectivity alone, then selects the sequence with lowest RMS upon
  molecular rotation. It is more reliable than UMEYAMA but may become
  expensive for larger molecules. This is the default.

Once the atoms are in the correct sequence, the RMS is obtained by
calculating the rotation matrix that yields the best match, in the
least squares sense, between the template and the reordered
structure(s). The MOLREORDER keyword does not require having any
molecular or crystal structure loaded.

If WRITE is used followed by a file name `file.s`, the output
structure with its atoms in the same order as the template is written
to that file. In addition, if MOVEATOMS is used, the atoms in the
structure written to `file.s` are moved to match the structure of the
template molecule. If INV is used, inversion operations are allowed.

## Move Atoms in a Molecular Crystal to Match Molecules (MOLMOVE) {#c2-molmove}

The `MOLMOVE` keyword is used to modify the atomic positions in a
molecular crystal in such a way that the structures of its fragments
are identical to a set of molecules provided by the user. This is
useful when the molecules in a molecular crystal are extracted (with
the `NMER 1 ONEMOTIF` options of the [WRITE](/critic2/manual/write/)
command), then their structure is modified, most commonly by running a
geometry relaxation in the gas-phase, and then the new molecular
structures need to be re-packed back into the crystal structure. The
syntax for the `MOLMOVE` keyword is:
~~~
MOLMOVE mol1.s mol2.s ... moln.s target.s new.s
~~~
The molecular crystal that needs to be modified (`target.s`) must be
composed of `n` fragments, and `n` molecular structures have to be
given first (`mol1.s` to `moln.s`). The molecular structures can be
given in any order, but their centers of mass and atomic sequence must
be the same as the fragments in the target molecular crystal. Using
the `NMER 1 ONEMOTIF` options to `WRITE` the molecules in a molecular
crystal ensures that the order and centers of mass of the molecular
fragments are correct. The file name where the repacked molecular
structure will be written must be given at the end of the command
(`new.s`).

## Calculate Dimensions of Uniform k-Point Grids (KPOINTS) {#c2-kpoints}

The KPOINTS keyword:
~~~
KPOINTS [rk.r] [RKMAX rkmax.r]
~~~
calculates the number of k-points in each direction (nk1, nk2, nk3) of
a uniform k-point grid with density given by length parameter $$R_k$$
(`rk.r`). More k-points are allocated in short directions of the
crystal according to [VASP's formula](https://www.vasp.at/wiki/index.php/KPOINTS):
\begin{equation}
n_i = \text{int}\left[\text{max}\left(1,R_k |{\rm b}_i| + \frac{1}{2}\right)\right]
\end{equation}
where $$n_i$$ is the number of k-points in direction i and $$|{\rm
b}_i|$$ is the length of the ith reciprocal lattice vector.

If the KPOINTS keyword is used without any further input, critic2
lists all the distinct uniform k-point grids for values of `rk.r` from
zero up to a value equal to `rkmax.r`, which can be modified with the
optional RKMAX keyword (default: 100 bohr). If a number `rk.r` follows
KPOINTS, the k-point grid dimensions are calculated only for that value.

## Print the Brillouin Zone Geometry (BZ) {#c2-bz}

The BZ keyword:
~~~
BZ
~~~
calculates and prints the geometry of the Brillouin cell for the
currently loaded crystal.
