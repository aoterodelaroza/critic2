
## Automatic Determination of Critical Points (AUTO) {#c2-auto}

~~~
AUTO [GRADEPS eps.r] [CPEPS eps.r] [NUCEPS neps.r] [NUCEPSH nepsh.r]
     [EPSDEGEN edeg.r] [DISCARD expr.s] [CHK] [DRY] [SEEDOBJ] ...
AUTO ... [CLIP CUBE x0.r y0.r z0.r x1.r y1.r z1.r]
AUTO ... [CLIP SPHERE x0.r y0.r z0.r rad.r]
AUTO ... [SEED ...] [SEED ...] ...
AUTO SEED WS [DEPTH depth.i] [X0 x0.r y0.r z0.r]
             [RADIUS rad.r]
AUTO SEED OH [DEPTH depth.i]  [X0 x0.r y0.r z0.r]
             [RADIUS rad.r] [NR nr.r]
AUTO SEED SPHERE [X0 x0.r y0.r z0.r] [RADIUS rad.r]
                 [NTHETA ntheta.i] [NPHI nphi.i] [NR nr.r]
AUTO SEED PAIR [DIST dist.r] [NPTS n.i]
AUTO SEED TRIPLET [DIST dist.r]
AUTO SEED LINE [X0 x0.r y0.r z0.r] [X1 x0.r y0.r z0.r]
               [NPTS n.i]
AUTO SEED POINT [X0 x0.r y0.r z0.r]
AUTO SEED MESH
~~~

The search for the critical points (CP) of a scalar field (the points
where the gradient of the field vanishes) is a basic task in the
quantum theory of atoms in molecules (QTAIM). In critic2, this search
is almost always conducted using the automatic CP localization
algorithm implemented in the AUTO keyword.

The automatic search for critical points has two steps: seeding and
searching. In the seeding step, a collection of points are selected in
the space that spans the crystal (the unit cell) or molecular
space. In the search step, a Newton-Raphson algorithm is launched at
each of the seeds in order to find nearby critical points.

The default seeding behavior in critic2 depends on whether the
geometry under study is a crystal (loaded with the
[CRYSTAL](/critic2/manual/crystal/#c2-crystal) keyword) or
a molecule ([MOLECULE](/critic2/manual/molecule/#c2-molecule)
keyword):

* In a crystal, the default behavior of AUTO is to calculate the
  Wigner-Seitz (WS) cell and its irreducible part (the smallest piece
  of the crystal that reproduces the WS cell by symmetry). Once the
  irreducible WS (IWS) cell is found, seed points are chosen by
  subdividing the edges, faces and interior of the tetrahedra forming
  the IWS comprises up to a certain subdivision level (the DEPTH).

* In a molecule, a single seed is planted at the midpoint between
  every pair of atoms less than 15 bohr apart.

Sometimes the default behavior will fail to locate some critical
points. For those situations, AUTO provides multiple seeding
strategies that can be employed by the user via the SEED
keyword. These include searches between pairs of atoms (PAIR), atomic
triplets (TRIPLET), uniform seeding in a sphere (SPHERE), a recursive
subdivision of an octahedron (OH), seeding along lines (LINE) and at
points (POINT), and seeding at a molecular integration mesh
(MESH). The seed list built using these "seeding actions" is pruned to
remove duplicates before the search is performed. Optionally, a
portion of the unit cell can be selected to restrict the search in
real space using the CLIP keyword. Once the seed list is built,
Newton-Raphson is applied at each of the seeds on the list, making
full use of shared-memory parallelization and crystal symmetry.

The default seeding procedure can be changed using one or more SEED
keywords. When a SEED keyword is used, the default seeding strategy is
forgotten by critic2 and manual control of the seeding is used
instead. Several SEED keywords can be used at the same time, each one
specifying a single seeding action that is determined by the keyword
immediately following SEED. This keyword can be:

* WS: a recursive subdivision of the Wigner-Seitz cell to level
  `depth.i` (keyword: DEPTH, default: 1). The WS cell can be displaced
  and centered somewhere in the unit cell using X0 (default: (0,0,0))
  and scaled down to a radius of `rad.r` (keyword: RADIUS, default:
  not used). The units for x0 are crystallographic coordinates in
  crystals and molecular Cartesian coordinates in molecules (default:
  angstrom unless changed by
  [UNITS](/critic2/manual/inputoutput/#c2-units)). For `rad.r`, the
  default units are bohr in crystals and angstrom in molecules.

* OH: a sphere of radius `rad.r` (keyword RADIUS, mandatory) is built
  around (`x0.r` `y0.r` `z0.r`) (keyword X0, mandatory). On the
  surface of that sphere, points are set according to a recursive
  subdivision algorithm that starts from a single octahedron and uses
  `depth.r` recursion levels (keyword: DEPTH, default: 1). This is the
  same procedure as the TRIANG keyword in
  [BASINPLOT](/critic2/manual/basinplot/#c2-basinplot). Each resulting
  point on the sphere surface determines a single ray along which
  `nr.r` seeds are uniformly distributed (keyword: NR, mandatory),
  from X0 to the surface of the sphere. The units for x0 are
  crystallographic coordinates in crystals and molecular Cartesian
  coordinates in molecules (default: angstrom unless changed by
  [UNITS](/critic2/manual/inputoutput/#c2-units)). For `rad.r`, the
  default units are bohr in crystals and angstrom in molecules.

* SPHERE: a sphere of radius `rad.r` (keyword RADIUS, mandatory) is
  built around (`x0.r` `y0.r` `z0.r`) (keyword X0, mandatory). On the
  surface of that sphere, points are uniformly distributed, with a
  placement algorithm that uses `nphi.i` points in the azimuthal angle
  (keyword: NPHI, mandatory) and `ntheta.i` points in the polar angle
  (keyword: NTHETA, mandatory). Each resulting point on the sphere
  surface determines a single ray along which `nr.r` seeds are
  uniformly distributed (keyword: NR, mandatory), from X0 to the
  surface of the sphere. The units for x0 are crystallographic
  coordinates in crystals and molecular Cartesian coordinates in
  molecules (default: angstrom unless changed by
  [UNITS](/critic2/manual/inputoutput/#c2-units)). For `rad.r`, the
  default units are bohr in crystals and angstrom in molecules.

* PAIR: seeds are placed on the interatomic lines, for all atom pairs
  at a distance less than `dist.r` (keyword: DIST, default: 15). The
  number of seeds per line is `n.i` (keyword: NPTS, default: 5). The
  default units for dist.r are bohr in crystals, angstrom in
  molecules.

* TRIPLET: seeds are placed at the barycenter of every atomic triplet
  in which the three atoms are at a distance from each other less than
  `dist.r` (keyword: DIST, default: 15). The default units for
  `dist.r` are bohr in crystals, angstrom in molecules.

* LINE: place `n.i` seeds (keyword: NPTS, default: 5) along a line
  between (`x0.r` `y0.r` `z0.r`) (keyword: X0, default: origin) and
  (`x1.r` `y1.r` `z1.r`) (keyword: X1, mandatory). The units for x0
  and x1 are crystallographic coordinates in crystals and molecular
  Cartesian coordinates in molecules (default: angstrom unless changed
  by [UNITS](/critic2/manual/inputoutput/#c2-units)).

* POINT: place a single seed at (`x0.r` `y0.r` `z0.r`) (keyword: X0,
  mandatory). The units for x0 are crystallographic coordinates in
  crystals and molecular Cartesian coordinates in molecules (default:
  angstrom unless changed by
  [UNITS](/critic2/manual/inputoutput/#c2-units)).

* MESH: place seeds at the nodes of a molecular integration mesh. The
  type of mesh can be controlled with the
  [MESHTYPE](/critic2/manual/misc/#c2-meshtype) keyword.

Multiple SEED keywords can be given in the same AUTO command. For
instance:
~~~
AUTO SEED PAIR SEED WS SEED POINT 1/4 1/4 1/4
~~~
executes three seeding actions: a search between all atoms pairs (1
seed per pair), a recursive subdivision of the WS cell (one level),
and a single seed at (0.25 0.25 0.25). The seed placement can be
visualized using the optional SEEDOBJ keyword. SEEDOBJ writes an OBJ
file (`<root>_seeds.obj`) containing the unit cell and all the seed
positions.

The AUTO search can be restricted to a portion of the unit cell using
the CLIP keyword. The CLIP keyword specifies a region of real
space. Only the seeds inside that region are used, and only the CPs
found inside that region are accepted (although, in crystals, symmetry
can replicate the CPs and send them outside the CLIP region; use
[MYM](/critic2/manual/crystal/#c2-symm) with the CLEAR keyword or
[NOSYMM](/critic2/manual/crystal/#c2-symm) to deactivate symmetry if
necessary). There are two possible region shapes in CLIP: a box (CUBE)
and a sphere (SPHERE). The box is specified by giving two opposite
corners: x0 and x1. The sphere requires a center (x0) and a
radius. The units for x0 and x1 are crystallographic coordinates in
crystals and molecular Cartesian coordinates in molecules (default:
angstrom unless changed by
[UNITS](/critic2/manual/inputoutput/#c2-units)). For `rad.r`, the
default units are bohr in crystals and angstrom in molecules.

A number of additional optional keywords control the behavior of
AUTO. GRADEPS is the gradient norm threshold for the optimization: if
a CP is found with gradient norm less than GRADEPS (default: 1e-12),
then it is accepted as CP.

The DISCARD keyword can be used to reduce the list of critical
points. If the expression `expr.s` evaluated at the critical point is
non-zero, the critical point is discarded. A typical use for this
keyword is when the system has a vacuum region. The (spurious)
critical points in the vacuum region can be eliminated from the list
by doing, for instance, `DISCARD "$rho < 1e-7"` to remove the critical
points with density lower than 1e-7. The arithmetic expression can
involve any number of fields, not just the reference field.

If DRY (dry run) is used, then the seeding is done but the actual CP
search is skipped. This is useful to examine the seed placement (in
combination with SEEDOBJ) and also to print the current list of CPs at
zero computational cost.

CPEPS controls the minimum distance between CPs to consider them
equivalent. The default units for CPEPS are bohr in crystals and
angstrom in molecules. Default: 0.2 bohr. Similarly, NUCEPS controls
the distance to consider a CP the same as a nucleus. If a CP is found
at a distance less than `neps.r` from the closest nucleus, critic2
considers the two are the same. The default value of NUCEPS is 0.1
bohr, except for fields defined on a grid, where it defaults to 2
times the maximum of the grid step in each direction. Hydrogens are
particular in that the maximum of the electron density may be
significantly displaced from the atomic position. A separate NUCEPS
distance criterion is provided for hydrogens (NUCEPSH), which defaults
to 0.2 bohr.

EPSDEGEN controls how degenerate critical points (i.e. one or more
eigenvalues of the Hessian is zero) are discarded. This is important
to get rid of critical points that appear in vacuum regions. A
critical point is considered degenerate if any of the elements of the
diagonal of the Hessian is less than `edeg.r` in absolute value.

Because finding the CPs can be an expensive task in large structures,
the CP list for the current field can be saved to a checkpoint file
using the CHK keyword. This keyword generates a `<root>.chk_cps` file
where the list of critical points is stored. It can be accessed in
subsequent critic2 runs by using the CHK keyword in AUTO. For
instance, to read the CPs from the checkpoint file and skip any
calculation in AUTO, you can do:
~~~
AUTO DRY CHK
~~~
By default, checkpoint files are not used.

### Example Output for AUTO in a Crystal

In crystals, critic2 writes all the critical points found by AUTO to
two internal lists: the "non-equivalent" list, containing only those
CPs not equivalent by symmetry, and the "complete" list, which
contains all the CPs in the unit cell. See the
[input and output notation](/critic2/manual/inputoutput/#c2-notation).
In molecules, symmetry is not used, so both lists are the same.

The most important part of the AUTO output is the "final report",
which gives the non-equivalent CP list and some other useful
information. Its appearance is (the table has been simplified):
~~~
* Critical point list, final report (non-equivalent cps)
  Topological class (n|b|r|c):   2(8) 1(16) 1(16) 2( 8)
  Morse sum :   0
# ncp pg type    position  mult name  f     |grad|     lap
1 Td  (3,-3) n 0.00 0.00 0.00  4 B  7.19E+1 0.00E+00 -3.00E+15
2 Td  (3,-3) n 0.25 0.25 0.25  4 P  2.36E+3 0.00E+00 -3.00E+15
3 C3v (3,-1) b 0.09 0.59 0.59 16 b1 1.25E-1 4.43E-17 -2.25E-01
4 C3v (3, 1) r 0.88 0.88 0.61 16 r1 1.08E-2 1.09E-16  3.45E-02
5 Td  (3, 3) c 0.75 0.75 0.75  4 c1 5.63E-3 1.34E-16  2.22E-02
6 Td  (3, 3) c 1.00 0.00 0.50  4 c2 7.28E-3 1.32E-16  2.58E-02
~~~
The non-equivalent CP list gives the type and position of all the
non-equivalent CPs found, their associated name, rank and signature,
and multiplicity. Critic2 also provides the site-symmetry (pg), and
the values of the reference field (f), its gradient (grad) and its
Laplacian (lap) evaluated at the CPs. The 'topological class' gives
the number of non-equivalent ncp, bcp, rcp, and ccp found. The values
in parentheses correspond to the total number of CPs of each class in
the cell. If the list is complete, then the Morse sum is zero (in a
crystal) or one (in a molecule).

Following the final report is the "analysis of system bonds" and
"analysis of system rings":
~~~
* Analysis of system bonds
# ncp end1 end2 r1(bohr) r2(bohr) r1/r2  r1-B-r2(degree)
  3    Mg    O   1.7789   2.1999 0.80864  179.9999

* Analysis of system rings
# ncp end1 end2 r1(bohr) r2(bohr) r1/r2  r1-R-r2(degree)
  4    c01  c01  1.9894   1.9894 1.00000  179.9999
~~~
These contain a list of all the bond (resp. ring) critical points
found and the corresponding nuclei (resp. cages) they connect. The
connected nuclei are found by tracing an upwards gradient path from
the bonds. The cages are found using a downwards gradient path from
the ring. In addition, it is also shown the geometric distance to the
connected critical point, the ratio between the distances, and the
angle formed by the bond path (ring path) at the bcp (rcp).

The **complete CP list** is the list of all CPs in the unit cell, and
comes next in the output. The identifiers from this list are used as
input for other keywords (for instance,
[GRDVEC](/critic2/manual/gradientpath/#c2-grdvec) or
[FLUXPRINT](/critic2/manual/gradientpath/#c2-fluxprint)) as they
specify a particular position in the crystal. The entries are similar
to the non-equivalent CP list (the table has been simplified):
~~~
* Complete CP list
  (x symbols are the non-equivalent representative atoms)
  cp ncp typ | x  | y  | z  | op.   (lvec+cvec)
x 1  1  n    1.00 0.50 0.32  1   1.0   0.0   0.0
  2  1  n    0.50 0.00 0.67  2   0.0   0.0   1.0
x 3  2  n    1.00 0.50 0.59  1   1.0   0.0   0.0
  4  2  n    0.50 0.00 0.40  2   0.0   0.0   1.0
[...]
~~~
The columns are, in order, the CP identifier (cp), the identifier for
the same CP in the non-equivalent CP list (ncp), the type of CP, the
position in crystallographic coordinates, and the symmetry operation
that transforms the CP from the non-equivalent CP list into the listed
CP. "Op." corresponds to one of the symmetry operations listed in the
output and "lvec+cvec" is a translation. Application of the operation
"op." to the non-equivalent CP followed by the translation recovers
the position of the CP in the complete list. The CPs that are at
exactly the same position as the corresponding CPs in the
non-equivalent list are listed first and marked with an "x". Note that
for these, the operation is always 1 (the identity) and the
translation vector is always a lattice translation.

The complete list is followed by a list of all the bcps and rcps in
the unit cell, together with the particular nuclei linked by the bond
path they represent (bcp) or with the cages linked by the ring path
(rcp):
~~~
* Complete CP list, bcp and rcp connectivity table
# (cp(end)+lvec connected to bcp/rcp)
#cp ncp typ   position    end1 (l) end2 (l)
1    1   n 0.00 0.00 0.00
[...]
9    3   b 0.50 0.50 0.77 4 (0 0 1) 5 (0 0 0)
10   3   b 0.50 1.00 0.27 3 (0 1 0) 6 (0 1 0)
[...]
33   4   r 0.75 0.25 1.00 63 (0 0 1) 57 (0 0 0)
34   4   r 0.75 0.75 0.50 64 (0 0 0) 58 (0 0 0)
~~~
This list is similar to the complete CP list, but the two entries at
the end (end1 and end2) give the two identifiers from the complete CP
list that the bond or ring is connected to, as well as the lattice
vector by which it should be translated to regenerate their actual
position.

More information about the atomic connectivity through bcps can be
found in the "attractor connectivity matrix":
~~~
* Attractor connectivity matrix
              n(1)  n(2)
               Mg     O
 n(1)    Mg     0     6
 n(2)    O      6     0
~~~
This list gives the number of bond paths between the different types
of non-equivalent atoms in the cell. For instance, in the example
above, each Mg (row 1) is bonded to six adjacent oxygens through six
bcps (second entry in the matrix) but there are no Mg-Mg bond
paths. Likewise, every oxygen is bonded to six Mg, but there are no
O-O bond paths.

The final part of the output from AUTO contains a detailed list of all
the non-equivalent critical points found, together with an exahustive
list of properties calculated at those points. More properties can be
calculated by using the
[POINTPROP](/critic2/manual/cpsearch/#c2-pointprop) keyword:
~~~
* Additional properties at the critical points
[...]
+ Critical point no. 4
  Crystallogrpahic coordinates: 0.75 0.25 1.00
  Cartesian coordinates (bohr): 5.96 1.98 7.95
  Type : (3,1)
  Field value (f): 2.350841441E-02
  Field value, valence (fval): 2.350841441E-02
  Gradient (grad f): 3.77e-16 -3.76e-16 2.65e-17
  Gradient norm (|grad f|): 5.340413578E-16
  Laplacian (del2 f): 2.889105008E-02
  Laplacian, valence (del2 fval): 2.889105008E-02
  Hessian eigenvalues: -7.02e-3  3.42e-3  3.24e-2
  Hessian:
     1.795930240E-02   1.453912685E-02  -1.850420411E-19
     1.453912685E-02   1.795930240E-02  -4.127816220E-19
    -1.850420411E-19  -4.127816220E-19  -7.027554725E-03
  Ellipticity (l_1/l_2 - 1): -3.054735092E+00
[...]
+ Flatness (rho_min / rho_b,max): 0.462267
~~~
The list of properties calculated by default contains the
crystallographic and Cartesian coordinates, the rank and signature of
the critical point, the value of the reference field at that point and
its derivatives (gradient, gradient norm, Laplacian, Hessian). The
ellipticity is the quotient between Hessian eigenvalues, only relevant
at bcps or rcps. The final entry, flatness, given at the end of the
list, is the quotient between the density minimum and the maximum
density at the bond critical points.

### Example Output for AUTO in a Molecule

The output for AUTO in a molecule is similar to a crystal but simpler,
since there is no distinction between the non-equivalent critical
point list and the complete CP list. The list of CPs is given first,
together with the class, and the Poincare-Hopf sum (which is 1 for a
complete list of CPs):
~~~
* Critical point list, final report (non-equivalent cps)
  Topological class (n|b|r|c): 12(12) 12(12) 1(1) 0(0)
  Poincare-Hopf sum: 1
# ncp type      position (ang_) name   f    |grad|   lap
  1  (3,-3) n  0.00  1.20  0.69  C  1.1e+2 0.0e+00 -4.3e+5
  [...]
  24 (3,-1) b -0.00  1.81 -1.04 b06 2.7e-1 2.8e-13 -9.6e-1
  25 (3,1 ) r -0.00 -0.00  0.00 r01 2.0e-2 1.0e-15 1.6e-1
~~~
Each row contains the CP identfier, the rank and signature, the type
of critical point, the position in Cartesian coordinates (angstrom
by default), the name, and the value, gradient and Laplacian of the
reference field at that point.

The CP list is followed by the analysis of the connectivity between
atoms via bonds and between cages via rings:
~~~
* Analysis of system bonds
# ncp end1 end2 r1(ang) r2(ang) r1/r2 r1-B-r2(degree)
 13    C    C   0.6979  0.6979  1.000   179.8198
 [...]
 24    C    H   0.6959  0.3907  1.781   179.9998

* Analysis of system rings
# ncp end1 end2 r1(ang) r2(ang) r1/r2 r1-R-r2(degree)
 25    ??  ??   6.7788  6.7788  1.000   180.0000
~~~
followed by the attractor connectivity matrix:
~~~
* Attractor connectivity matrix
              n(1)  n(2)  n(3)  n(4)  n(5)  n(6)
               C     H     C     H     C     H
 n(1)    C      0     1     1     0     0     0
 [...]
n(12)    H      0     0     0     0     0     0
              n(7)  n(8)  n(9) n(10) n(11) n(12)
               C     H     C     H     C     H
 n(1)    C      0     0     0     0     1     0
 [...]
n(12)    H      0     0     0     0     1     0
~~~
These lists have the same meaning as in the crsytal example.

Finally, the output gives an exhaustive list of properties at the
critical points. As in the case of crystals, the list of properties
calculated at the critical points can be varied using the
[POINTPROP](/critic2/manual/cpsearch/#c2-pointprop) keyword:
~~~
* Additional properties at the critical points
[...]
+ Critical point no. 4
  Coordinates (ang_): 0.0000000000 2.1497999997 -1.2411900006
  Type : nucleus
  Field value (f): 4.254133728E-01
  Field value, valence (fval): 4.254133728E-01
  Gradient (grad f): 0.00e+00  0.00e+00  0.00e+00
  Gradient norm (|grad f|): 0.000000000E+00
  Laplacian (del2 f): -1.951050589E+01
  Laplacian, valence (del2 fval): -1.951050589E+01
  Hessian eigenvalues: -6.64e+00  -6.63e+00  -6.22e+00
  Hessian:
    -6.639684283E+00  -4.772491252E-15   2.755444887E-15
    -4.772491252E-15  -6.329875105E+00  -1.827977265E-01
     2.755444887E-15  -1.827977265E-01  -6.540946498E+00
  Field 0 (f,|grad|,lap): 3.16e-01 0.00e+00 -7.26e+00
[...]
~~~
The only difference with the output in a crystal is that the
coordinates are given in Cartesian (angstrom by default) and referred
to the molecular origin. Also, the flatness is missing, because it is
meaningless in a molecule.

### Problems Finding Critical Points

Sometimes the zero (one) sum condition between the number of critical
points of each type is not fulfilled. This is usually caused by
problems in the localization of critical points that originate from
the (often unavoidable) numerical shortcomings of the scalar field.

In WIEN2k and elk densities, there might be spurious CPs at the
surface of the muffin tin, where the density is discontinuous. These
spurious CPs show up in the final report list. Every time a FPLAPW
field is loaded, every atomic muffin is checked for discontinuities
and the report printed to the output.

In DFTB+ and other fields with missing core electrons (e.g. plane-wave
calculations), the use of core-augmentation (ZPSP and CORE keywords)
is recommended to prevent the appearance of spurious critical points
close to the nuclei.

In scalar fields with extreme variations in value (e.g. the Laplacian
of the electron density), it is unlikely that AUTO will find all the
core CPs, since they may be very close to each other or the CP may
have a very small attraction basin in the Newton-Raphson algorithm. A
previous version of critic2 does that (available upon request) but
since the core CPs are not all that interesting, we have decided to
remove that code from this version.

By far, the most problematic type of field for CP localization is a
grid of values. The difficulty is in finding a way to calculate the
first and second derivatives of the field at arbitrary points in space
(i.e. doing an interpolation) in a way that is both accurate and
smooth. These problems have been considered in a 2022 article
([Otero-de-la-Roza, J. Chem. Phys. 156, 224116 (2022)](https://doi.org/10.1063/5.0090232))
and in [this example](/critic2/examples/example_08_01_gridcps/). In
short, if you are looking for the critical points of an all-electron
density given on a grid, then the SMOOTHRHO interpolation is your best
option.

Regions where the value of the field is very small are subject to the
appearance of many spurious critical points, due to the very small
value of the gradient. These regions typically appear in the vacuum
region of a slab, or far away from a molecular system. To eliminate
the critical points in these regions, use the DISCARD keyword combined
with a threshold for the field value. For instance, to discard all
critical points with density less (loaded in field `$rho`) than 1e-5,
use `DISCARD "$rho < 1e-5"`.

### Visualization of Critical Points

Critic2 does not come with a graphical interface to visualize the
results of your calculations (yet). This can be a problem when the
list of critical points is very large, if you want to identify one
particular critical point among many, or if you need to calculate
distances between atoms and/or critical points. However, there are
some molecular visualization programs that can be used to read the
output of critic2, and this section describes the most comfortable
procedures to do so at present.

The key to critical point visualization is the CPREPORT keyword. When
used with one of the available output file formats, critic2 will write
the list of atoms plus the critical points in that format. The
critical points are labeled as "Xn" (nuclear CP), "Xb" (bonds), "Xr"
(rings), and "Xc" (cages). The list of critical points is written to
the output in the same order as in critic2's complete critical point
list, unless one of the special keywords to modify the plotting region
is used (e.g. MOLMOTIF). This means that, if your graphical program
gives you the atom labels, atom number n in the GUI will correspond
exactly to critical point number n in critic2.

The challenege, therefore, is to make external programs understand
that we have objects (critical points) that should be treated like
atoms, but are not atoms. An easy (and sometimes very convenient) way
of doing this is to replace the critical point labels ("Xb", "Xr",...)
with atoms we know are not present in our system (e.g. H).

#### Deprecated visualization method using avogadro

Previously, the way of visualizing the
critical points this was using [avogadro](http://avogadro.cc), an
open-source visualizer for crystals and molecules. Avogadro uses
[openbabel](http://openbabel.org) as the underlying format
converter. In order to make avogadro understand the critical point
types, all you need to do is modify the `element.txt` file that comes
with openbabel, and add the critical point types at the end of the
atomic species list:
~~~
[...]
118     Uuh     0.00    1.60    0.00    0.00    6       294     0.00    0       0       0.99    0.00    0.06    Ununoctium
119     Xn      0.00    0.50    0.00    0.50    0       0       0.00    0       0       0.2805  0.6223  0.0164  nCP
120     Xb      0.00    0.50    0.00    0.50    0       0       0.00    0       0       1.0000  0.8512  0.2380  bCP
121     Xr      0.00    0.50    0.00    0.50    0       0       0.00    0       0       0.5851  0.5349  1.0000  rCP
122     Xc      0.00    0.50    0.00    0.50    0       0       0.00    0       0       1.0000  0.3998  0.3411  cCP
123     Xz      0.00    0.20    0.00    0.20    0       0       0.00    0       0       0.1709  1.0000  0.0000  xCP
~~~
With these changes, atoms and critical points will be represented, and
the latter will use the color scheme above. The recommended file
format both for molecules and crystals is cml, which has several
advantages over xyz and cif: i) it prevents avogadro from calculating
the molecular connectivity, which can be expensive in a very large
system, and ii) it allows showing a fragment of the crystal along with
the unit cell. To access the critical point labels (which correspond
to the critic2 labels), go to Display types, mark "Label", click on
the Label options, and select "Text: atom number".

The capability of avogadro to handle periodic systems is limited, and
more recent versions of avogadro have taken the program in a somewhat
confusing direction, so the above may not apply in your case.

#### Recommended visualization method using vmd

The current best way of visualizing the critical points in your system
is using a vmd script. Write the critical points using:
~~~
CPREPORT name.vmd [GRAPH]
~~~
where the GRAPH keyword can be used to write the bond paths. This
command writes a vmd script (`name.vmd`) as well as the geometry of
the molecule or crystal, including the critical point and gradient
paths (`name.xyz`). To open this file with vmd, execute the command:
~~~
vmd -e name.vmd
~~~
or open vmd and select the "File" menu followed by "Load visualization
state".

In this representation, all atoms in the unit cell (in crystals) or in
the molecule are shown, colored by element. Bond critical points are
in orange, rigns are cyan, cages are red and non-nuclear maxima are
green. Bond  paths are shown in pink. In the case of a
crystal, the unit cell and lattice vectors are also shown. The lattice
vectors are also known to vmd so you can repeat the unit cell in the
"Periodic" tab of the "Graphics -> Representations" menu. The blue
labels appearing on the plot are the atoms in the system, in the order
in which the appear in the complete list.

If you open the Representations window ("Graphics ->
Representations"), you will see there are at most 6 representations
(you may have fewer if the system lacks certain types of critical
points). Types 1, 2, 3, and 4 correspond to nuclear, bond, ring, and
cage critical points, respectively. Type 5 are the bond
paths. The representation with "Element" coloring method are the
atoms. You can double click any of these representations to hide some
of these objects.

Lastly, it is possible to identify a particular atom or critical
point. To do this, select "Mouse -> Pick" and click on an atom or
critical point. A message will be written in the vmd terminal
containing, for instance:
~~~
...
Info) name: Xr138
...
Info) index: 137
...
~~~
which indicates that this is the ring critical point number 138 in the
complete CP list (the critical point ID is always the index minus
1). When picking atoms and critical points, it is a good idea to hide
the bond path representations, to avoid misclicks.

Of course, vmd gives almost infinite variety when it comes to
tailoring your plot to suit your needs, so changing the script or the
contents of the xyz file is encouraged. There is, in fact, no
information about the system in the `.vmd` file other than the lattice
vectors and unit cell information when representing a periodic
crystal.

## Requesting More Information About the Critical Point List (CPREPORT) {#c2-cpreport}

~~~
CPREPORT {SHORT|LONG|VERYLONG|SHELLS [n.i]}
CPREPORT file.{xyz,gjf,cml,vmd} [SPHERE rad.r [x0.r y0.r z0.r]]
         [CUBE side.r [x0.r y0.r z0.r]] [BORDER] [ix.i iy.i iz.i]
         [MOLMOTIF] [ONEMOTIF] [ENVIRON dist.r]
         [NMER nmer.i]
CPREPORT file.{obj,ply,off} [SPHERE rad.r [x0.r y0.r z0.r]]
         [CUBE side.r [x0.r y0.r z0.r]] [BORDER] [ix.i iy.i iz.i]
         [MOLMOTIF] [ONEMOTIF] [CELL] [MOLCELL]
CPREPORT file.scf.in
CPREPORT file.tess
CPREPORT file.cri|file.incritic
CPREPORT {[file.]POSCAR|[file.]CONTCAR}
CPREPORT file.abin
CPREPORT file.elk
CPREPORT file.gau
CPREPORT file.cif
CPREPORT file.m
CPREPORT file.gin
CPREPORT file.lammps
CPREPORT file.fdf
CPREPORT file.STRUCT_IN
CPREPORT file.hsd
CPREPORT file.gen
CPREPORT file.json
CPREPORT file.test
CPREPORT [...] [GRAPH]
~~~

CPREPORT prints additional information about the critical points to
the output.

* SHORT: prints the list of non-equivalent critical points.
* LONG: in a crystal, print the complete list of critical points in
  the unit cell and the connectivity in the case of bcp and rcp (when
  the graph is calculated); in a molecule, it is the same as SHORT.
* VERYLONG: detailed information at every critical point, including
  the derivatives of the reference field, the evaluation of all other
  fields, and the flatness (in a crystal only).
* SHELLS: local neighbor environment of every critical point (up to
  `n.i` shells, default 10).

In addition, any of the file formats available in the WRITE command
can also be used in CPREPORT to write the molecular or crystal
structure plus the critical point list in the relevant region. The
behavior of these options is analogous to
[WRITE](/critic2/manual/write/) (except in that the critical points
are written, in addition to the atoms). The critical points are
written using special symbols: Xn for a nuclear CP, Xb for a bond, Xr
for a ring, and Xc for cage. Miscellaneous or unassigned critical
points, and points along a gradient path are labeled Xz. Unless
special options are used to change the region being writteh by
CPREPORT (such as MOLMOTIF, BURST, etc.), the atoms and critical
points in the output file are in the same order as in critic2's
complete critical point list.

Critic2 can also write the structure and critical point information
for the reference field to a JavaScript Object Notation (JSON) file by
indicating a file with extension `.json`. The generated JSON file
contains the system geometry, the field details, and essentialy the
same information as in the SHORT, LONG, and VERYLONG reports.

The `.test` file format is used in critic2's internal tests. This file
format is intended for debugging purposes only.

The optional GRAPH keyword can be used in combination with any of the
file formats mentioned above. When GRAPH is used, the bond paths are
calculated and represented.

## List of Properties Calculated at Points (POINTPROP) {#c2-pointprop}

The default output for AUTO contains a list of detailed information
at the critical points. This is an example of what this section looks
like:
~~~
* Additional properties at the critical points
[...]
+ Critical point no. 9
  Crystallogrpahic coordinates: 0.65183 0.90869 0.72597
  Cartesian coordinates: 12.89298 22.10026 33.60159
  Type : (3,-1)
  Field value (f): 4.939109713E-03
  Field value, valence (fval): 4.939109713E-03
  Gradient (grad f): -9.5115E-18 3.8415E-17 4.1112E-18
  Gradient norm (|grad f|): 3.978845324E-17
  Laplacian (del2 f): 1.312332087E-02
  Laplacian, valence (del2 fval): 1.312332087E-02
  Hessian eigenvalues: -2.11139E-03 -1.48096E-03 1.67156E-02
  Hessian:
    -2.068242883E-04  -5.185807347E-03   1.179735612E-03
    -5.185807347E-03   1.373314257E-02  -4.510661176E-03
     1.179735612E-03  -4.510661176E-03  -4.029974096E-04
~~~
For each critical point, the coordinates (Cartesian and
crystallographic in a crystal; only Cartesian in a molecule), the
type, and the evaluation of the reference field and its derivatives is
given. In many cases, it is of interest to calculate the value of a
different field, or an arithmetic expression involving several
fields, at those critical points.

To obtain more information at the critical points of the reference
field, the procedure in critic2 is to register a field or an
arithmetic expression involving known fields in the "properties
list", accesible using the POINTPROP keyword, with syntax:
~~~
POINTPROP name.s "expr.s"
POINTPROP shorthand.s
POINTPROP CLEAR
POINTPROP LIST
~~~
The POINTPROP keyword associates the expression `expr.s` with the name
`name.s` and registers that name in a list of properties. When AUTO is
run (or [CPREPORT](/critic2/manual/cpsearch/#c2-cpreport),
if the POINTPROP order comes after AUTO), those
arithmetic expression will be applied to each of the CPs and the
result printed in the output. For instance, if one does:
~~~
POINTPROP MYGTF gtf(1)
POINTPROP STH log($1^2)
~~~
Then the result of AUTO for the critical point above becomes:
~~~
+ Critical point no. 9
  Crystallogrpahic coordinates: 0.65183 0.90869 0.72597
  Cartesian coordinates: 12.89298 22.10026 33.60159
  Type : (3,-1)
  Field value (f): 4.939109713E-03
  Field value, valence (fval): 4.939109713E-03
  Gradient (grad f): -9.5115E-18 3.8415E-17 4.1112E-18
  Gradient norm (|grad f|): 3.978845324E-17
  Laplacian (del2 f): 1.312332087E-02
  Laplacian, valence (del2 fval): 1.312332087E-02
  Hessian eigenvalues: -2.11139E-03 -1.48096E-03 1.67156E-02
  Hessian:
    -2.068242883E-04  -5.185807347E-03   1.179735612E-03
    -5.185807347E-03   1.373314257E-02  -4.510661176E-03
     1.179735612E-03  -4.510661176E-03  -4.029974096E-04
  mygtf (gtf(0)): 4.112914774E-04
  sth (log($0^2)): -1.062114037E+01
~~~
The properties in the list are calculated at the end. The properties
list is also used in other parts of critic2, notably in the output of
[POINT](/critic2/manual/graphics/#c2-point).

Any arithmetic expression can be used in POINTPROP, but it is common
to use one of the
[chemical functions from the critic2 function library](/critic2/manual/arithmetics/#availchemfun).
The shorthand names for the chemical functions can also be used to
apply those functions to the reference field. For instance:
~~~
POINTPROP GTF
~~~
activates the calculation of the Thomas-Fermi kinetic energy density
(`gtf` function) on the reference field. POINTPROP can only be used
with arithmetic expressions involving known fields. The keyword
[CLEAR](/critic2/manual/arithmetics/#c2-clear)
deletes all the properties in the list. The list of properties can be
accesed at any time using POINTPROP LIST.

## Examples

- AUTO:

  + [Locating Critical Points of Densities Given as Grids](/critic2/examples/example_08_01_gridcps/)
