
## Notation Used in this Manual {#c2-notation}

The input for critic2 is free-format and case-insensitive. Lines
preceded by `#` are treated as comments. The input syntax is
keyword-driven: the first word in any (non-blank) input line
determines the task to be carried out.

In this manual, keywords are written in CAPS. Input variables are
denoted using a suffix to indicate their type: a real number (`.r`),
an integer (`.i`) or a string (`.s`). Almost anywhere that a number is
expected, it is possible to use an
[arithmetic expression](/critic2/manual/arithmetics/). If an arithmetic expression
is required (and not merely optional), then it should be given in
either double or single quotes are used (for instance, "$1+$2"). When
several alternative keywords are possible, the or symbol (|) is
used. Square brackets ([]) denote optional keywords and curly braces
({}) are used for grouping.

Some of the sections in the rest of this manual contain a subsection
describing additional options. These provide some independent keywords
that control the behavior of critic2 and are meant to be used either
before or after the keywords in the section in which they appear. For
instance, NOSYMM can be used before CRYSTAL to deactivate the
automatic calculation of the crystal symmetry:
~~~
NOSYMM
CRYSTAL benzene.cif
~~~
Hence, NOSYMM appears in this manual as an additional option to
CRYSTAL.

The critic2 output is mostly self-explanatory, although there are a
number of key concepts that need to be understood. In the case of
crystals, the crystal motif is represented in critic2 by two lists of
atoms. The **non-equivalent atom list** contains the atoms in the
asymmetric unit, that is, the minimal list of atoms that generate all
the atomic positions in the crystal by symmetry. The **cell atom
list**, equivalently called the **complete atom list**, contains all
atoms in the unit cell. The non-equivalent atom list reproduces the
complete atom list by applying all symmetry operations (other than
pure lattice translations) known to critic2.

The atoms in each of those two lists are numbered. The integer
identifier for an atom in the non-equivalent atom list is symbolized
in this manual as `nat.i`. The integer identifier for an atom in the
complete atom list is `at.i`. The
[`syntax.txt` file](/critic2/syntax/) contains a summary of all
available keywords and follows the same notation. The distinction
between atoms from the complete list and from the non-equivalent list
is irrelevant in the case of molecules, because symmetry is not used
(hence both lists are the same).

In exact parallel to the atomic lists, critic2 also maintains a list
of non-equivalent critical points (CP) and a complete (or cell) list
of critical points found for a given scalar field. The CPs in the
non-equivalent list reproduce all the CPs in the complete list by
symmetry, and both lists are the exact same in molecules because
symmetry is not used. The keyword definitions use `ncp.i` for the
integer indices in the non-equivalent CP list, and `cp.i` for the
integer identifiers in the complete CP list.

In critic2 atoms are considered critical points always (this may
change at some point). Therefore, the non-equivalent (complete) atom
list is a subset of the non-equivalent (complete) CP list. Critic2
makes sure that the integer identifier for all atoms in the atom lists
are the same as in the corresponding CP lists. For instance, if an
atom has index `nat.i = 2` and `at.i = 5`, then necessarily `ncp.i =
2` and `cp.i = 5`, regardless of how many additional critical points
have been found for this system.

In most critic2 keywords, atoms can be selected by their atomic symbol
(`at.s` in the syntax definitions), in which case the keyword applies
to all atoms with the same atomic number unless otherwise
stated. Atoms can also be selected by an integer identifier from the
non-equivalent atom list (`nat.i` in the definitions) in those cases
in which symmetry makes it irrelevant which of the symmetry-equivalent
atoms in the cell are used. For example, the non-equivalent atom
identifier can be used to instruct critic2 to calculate the charge of
a certain atom, since all symmetry-equivalent atoms have the same
charge.

Some additional notation and terms that are used in the manual:

* By **scalar field** or **field** we mean a numerical or analytical
  representation of a function that associates a scalar value to every
  point in space. Often, this function is the electron density, for
  which special techniques are provided (for instance, core
  augmentation in the case of valence densities, see
  [ZPSP](/critic2/manual/crystal/#c2-charge)). However, critic2 can
  deal with any scalar field, and examples other than the density are
  the ELF, the Laplacian, etc.

* The **promolecular density** is the scalar field built by using the
  sum of the in-vacuo atomic densities. This object comes up in a
  number of contexts. For instance,
  [NCIPLOT](/critic2/manual/nciplot/) and
  [HIRSHFELD](/critic2/manual/misc/#c2-hirshfeld) use it. The
  promolecular density does not require any input from the user other
  than the crystal or molecular structure, and is always available
  under field identifier $0 (or $rho0).

* We denote by \<root\> the root of the input file, i.e., the name of
  the file minus its extension. If no input file is known (for
  instance, because critic2 is being run interactively), then
  the root defaults to "stdin". The default root can be changed with
  the keyword [ROOT](/critic2/manual/misc/#c2-root).

* The critical points of a field can be classified by their **rank**
  (r) and **signature** (s). The rank is the number of non-zero
  eigenvalues of the Hessian. In the vast majority of cases, r is 3.
  The signature is the number of positive eigenvalues minus the
  number of negative eigenvalues. s = -3 is a maximum, s = -1 is a
  first-order saddle point, s = +1 is a second-order saddle point and
  s = 3 is a minimum. These four types of critical points receive
  special names: nuclear CP, bond CP, ring CP, and cage CP,
  respectively. The abbreviations ncp, bcp, rcp, and ccp are also used
  throughout the manual and in the output. Note that a maximum is
  always a "nuclear critical point" even though it may not be
  associated to any nucleus.

## Input and Output Units {#c2-units}

The default input and output units in critic2 are bohr for crystals
(if the structure is loaded using the
[CRYSTAL](/critic2/manual/crystal/#c2-crystal) keyword) and angstrom
for molecules (if the
[MOLECULE](/critic2/manual/molecule/#c2-molecule) keyword is
used). In the particular case of molecules, the origin is also placed
at the same point as in the structure provided by the user.

This default behavior can be changed with the UNITS keyword:
~~~
UNITS {BOHR|AU|A.U.|ANG|ANGSTROM}
~~~
This command changes the units of all distances in input and output to
either bohr or angstrom.

## Simple Input and Output for a Crystal {#c2-simplecrystal}

As an example, let us consider an input for the conventional cell of
the fluorite (CaF2) crystal. Fluorite is cubic with space group
Fm-3m. Ca forms a face-centered cubic lattice and F occupies all the
tetrahedral voids. The input is:
~~~
CRYSTAL
 SPG F m -3 m
 CELL 5.463 5.463 5.463 90 90 90 ang
 Ca 0   0   0
 F  1/4 1/4 1/4
ENDCRYSTAL
~~~
When this structure is read, the non-equivalent atom list contains two
atoms: Ca at (0,0,0) with multiplicity 4 and F at (1/4,1/4,1/4) with
multiplicity 8. The cell atom list contains twelve atoms: four Ca
atoms at (0,0,0), (1/2,1/2,0), etc. and eight F atoms at
(1/4,1/4,1/4), (3/4,1/4,1/4), etc.

The output for this example follows. First, the output gives the
header with some information about the system, the version (the commit
number), and the location of the relevant library and density files:
~~~
* CRITIC2: analysis of real-space scalar fields in solids and molecules.
  (c) 1996-2019 A. Otero-de-la-Roza, A. Martin-Pendas, V. Lua~na
  Distributed under GNU GPL v.3 (see COPYING for details)
  Bugs, requests, and rants: aoterodelaroza@gmail.com
  Website: https://aoterodelaroza.github.io/critic2/
  If you find this software useful, please cite:
  A. Otero-de-la-Roza et al., Comput. Phys. Commun. 185 (2014) 1007-1018.
  A. Otero-de-la-Roza et al., Comput. Phys. Commun. 180 (2009) 157-166.

+ critic2 (development), version 1.0 (git:1.1dev-152-gb8eb182)
		 host: Linux-5.16.0-6-amd64
		 date: Wed 27 Apr 2022 10:20:16
 compiled dat: /usr/local/share/critic2
	  datadir: /home/alberto/git/critic2/dat
	 dic file: /home/alberto/git/critic2/dat/cif/cif_core.dic
...was found?: T
	   spglib: 1.13.0
		libxc: <unavailable>

CRITIC2--2022/4/29, 07:03:03.280
~~~

After the CRYSTAL keyword is read, critic2 first lists the basic
information about the crystal (note that the input lines read are
copied to the output preceded by the "%%" prefix). The output
starts with the cell parameters and the number of atoms in the crystal
motif:
~~~
%% CRYSTAL
%% SPG f m -3 m
%% CELL 5.463 5.463 5.463 90 90 90 ANG
%% NEQ 0 0 0 ca
%% NEQ 1/4 1/4 1/4 f
%% ENDCRYSTAL
* Crystal structure
  From: <input>
  Lattice parameters (bohr): 10.323574  10.323574  10.323574
  Lattice parameters (ang): 5.463000  5.463000  5.463000
  Lattice angles (degrees): 90.000  90.000  90.000
  Empirical formula:
	ca(1) f(2)
  Number of non-equivalent atoms in the unit cell: 2
  Number of atoms in the unit cell: 12
  Number of atomic species: 2
  Number of electrons (with zero atomic charge): 152
~~~

Next is the list of **atomic species**. Internally, critic2 keeps a
list of all the types of atoms present in the crystal or
molecule. Normally, each atomic species corresponds to a different
element but in some cases, for instance magnetic systems, it may
be useful to differentiate between two different atomic types with
the same atomic number.
~~~
+ List of atomic species:
# spc  Z   name    Q   ZPSP
   1  20    ca     0.0  --
   2   9     f     0.0  --
~~~
In this case, however, we have two species corresponding to Ca
and F, each with their corresponding atomic number.

Next comes the non-equivalent atom list. In this case, the whole
crystal is generated by replicating two atoms: one Ca and one F. The
positions, multiplicities, and the atomic numbers are indicated:
~~~
+ List of non-equivalent atoms in the unit cell (cryst. coords.):
# nat       x              y              z        spc  name   mult  Z
   1   0.0000000000   0.0000000000   0.0000000000   1    ca     4  20
   2   0.2500000000   0.2500000000   0.2500000000   2     f     8   9
~~~

The next table is the complete atom list. Here, critic2 lists all the
atoms in the unit cell: four Ca and eight F. The exact same list is
repeated in Cartesian coordinates, referred to the internal coordinate
system used in critic2. The output also indicates the matrix of
lattice vectors in Cartesian coordinates (repeated below).
~~~
+ List of atoms in the unit cell (cryst. coords.):
# at        x              y              z        spc  name    Z
   1   0.0000000000   0.0000000000   0.0000000000   1    ca    20
   2   0.0000000000   0.5000000000   0.5000000000   1    ca    20
   3   0.5000000000   0.0000000000   0.5000000000   1    ca    20
   4   0.5000000000   0.5000000000   0.0000000000   1    ca    20
   5   0.2500000000   0.2500000000   0.2500000000   2     f     9
   6   0.2500000000   0.7500000000   0.7500000000   2     f     9
   7   0.7500000000   0.2500000000   0.7500000000   2     f     9
   8   0.7500000000   0.7500000000   0.2500000000   2     f     9
   9   0.7500000000   0.7500000000   0.7500000000   2     f     9
  10   0.7500000000   0.2500000000   0.2500000000   2     f     9
  11   0.2500000000   0.7500000000   0.2500000000   2     f     9
  12   0.2500000000   0.2500000000   0.7500000000   2     f     9

+ Lattice vectors (bohr)
	a:   10.3235738640     0.0000000000     0.0000000000
	b:    0.0000000000    10.3235738640     0.0000000000
	c:    0.0000000000     0.0000000000    10.3235738640

+ List of atoms in Cartesian coordinates (bohr):
# at         x                y                z         spc  name    Z     dnn
   1     0.0000000000     0.0000000000     0.0000000000   1    ca    20    4.4702
   2     0.0000000000     5.1617869320     5.1617869320   1    ca    20    4.4702
   3     5.1617869320     0.0000000000     5.1617869320   1    ca    20    4.4702
   4     5.1617869320     5.1617869320     0.0000000000   1    ca    20    4.4702
   5     2.5808934660     2.5808934660     2.5808934660   2     f     9    4.4702
   6     2.5808934660     7.7426803980     7.7426803980   2     f     9    4.4702
   7     7.7426803980     2.5808934660     7.7426803980   2     f     9    4.4702
   8     7.7426803980     7.7426803980     2.5808934660   2     f     9    4.4702
   9     7.7426803980     7.7426803980     7.7426803980   2     f     9    4.4702
  10     7.7426803980     2.5808934660     2.5808934660   2     f     9    4.4702
  11     2.5808934660     7.7426803980     2.5808934660   2     f     9    4.4702
  12     2.5808934660     2.5808934660     7.7426803980   2     f     9    4.4702
~~~

Following this information comes the cell volume, in atomic units and
in angstrom^3:
~~~
+ Cell volume (bohr^3): 1100.24704
+ Cell volume (ang^3): 163.03979
~~~

And then the list of symmetry operations and the space group and point
group information:
~~~
+ List of symmetry operations (48):
  Operation 1:
	 1.000000  0.000000  0.000000  0.000000
	 0.000000  1.000000  0.000000  0.000000
	 0.000000  0.000000  1.000000  0.000000
  Operation 2:
	 0.000000  0.000000 -1.000000  0.000000
	-1.000000  0.000000  0.000000  0.000000
	 0.000000 -1.000000  0.000000  0.000000
[...]
  Operation 48:
	 0.000000  1.000000  0.000000  0.000000
	 0.000000  0.000000  1.000000  0.000000
	-1.000000  0.000000  0.000000  0.000000

+ List of symmetry operations in crystallographic notation:
   1: x,y,z
   2: -z,-x,-y
   3: -y,x,z
[...]
   47: x,-z,-y
   48: y,z,-x

+ List of centering vectors (4):
  Vector 1: 0.000000  0.000000  0.000000
  Vector 2: 0.000000  0.500000  0.500000
  Vector 3: 0.500000  0.000000  0.500000
  Vector 4: 0.500000  0.500000  0.000000

+ Crystal symmetry information
  Space group (Hermann-Mauguin): Fm-3m (number 225)
  Space group (Hall): -F 4 2 3 (number 523)
  Point group (Hermann-Mauguin): m-3m
  Point group (Schoenflies): Oh
  Holohedry: cubic
  Laue class: m-3m
~~~

The Cartesian to crystallographic ("car to crys") transformation
matrices are the transformation operations between the vector basis
formed by the cell vectors (crystallographic coordinates) and the
internal Cartesian axes used in critic2 (Cartesian coordinates). The
crystallographic to Cartesian matrix ("crys to car") gives the cell
vectors in Cartesian axes. The metric tensor is the transpose of "crys
to car" times "crys to car", and contains the scalar products between
lattice vectors.
~~~
+ Car/crys coordinate transformation matrices:
  A = car to crys (xcrys = A * xcar, bohr^-1)
	   0.0968656798    -0.0000000000    -0.0000000000
	   0.0000000000     0.0968656798    -0.0000000000
	   0.0000000000     0.0000000000     0.0968656798
  B = crys to car (xcar = B * xcrys, bohr)
	  10.3235738640     0.0000000000     0.0000000000
	   0.0000000000    10.3235738640     0.0000000000
	   0.0000000000     0.0000000000    10.3235738640
  G = metric tensor (B'*B, bohr^2)
	 106.5761773245     0.0000000000     0.0000000000
	   0.0000000000   106.5761773245     0.0000000000
	   0.0000000000     0.0000000000   106.5761773245
~~~

The next block in the output gives the calculation of the discrete
molecular units in this crystal. By default, critic2 calculates the
full atomic connectivity graph in any given input. Based on this
graph, critic2 determines whether the crystal is composed of molecules
(0D), polymers (1D), slabs (2D) or if it is a three-dimensional
periodic structure. In the case of a molecular crystal, critic2 gives
the list of all molecules in the unit cell and assigns to them an
identifying integer. In this case, the crystal is not molecular:
~~~
+ List of fragments in the system (1)
# Id = fragment ID. nat = number of atoms in fragment. C-o-m = center of mass (bohr).
# Discrete = is this fragment finite?
# Id  nat           Center of mass            Discrete
  1    12     0.871668    0.871668    0.871668   No

+ This is a 3D periodic structure.
~~~

In order to perform efficient calculations of distances in the
crystal, critic2 sets up an **atomic environment**, which is
the collection of all atoms in a number of unit cells surrounding the
cell at the origin. These atoms are placed into bins, which are then
used to speed up the calculation of distances and atomic contributions
to scalar fields. The information about the atomic environments is
given next:
~~~
+ Atomic environment
  Number of atoms (reduced cell/environment): 12 / 3056
  Radius of (unit cell/environment) circumscribed sphere (bohr): 8.9405 / 67.0536
  Maximum interaction distance (bohr): 31.4559
  Covering regions:
	Total number of regions: 216 (6 6 6)
	Minimum region ID: -3 -3 -3
	Maximum region ID: 2 2 2
	Region side (bohr): 13.3037
	Transformation origin (bohr): 5.1618,5.1618,5.1618
	Search offsets: 2197
	Maximum search offset: 6
	Average number of atoms per region: 14.1481
	Maximum number of atoms in a region: 40
~~~

Next the Wigner-Seitz (WS) cell and the reduced Delaunay cell are
determined. This information is used in
the calculation of shortest distances between lattice translations of
two atoms, and for the
[AUTO](/critic2/manual/cpsearch/#c2-auto) and
[YT](/critic2/manual/integrate/#c2-yt) tasks, among other things.
~~~
+ Vertex of the WS cell in cryst. coords. (8)
# id = vertex ID. xyz = vertex cryst. coords. d = vertex distance to origin (bohr).
   id       x            y            z          d (bohr)
	1    0.500000     0.500000    -0.500000     8.94047722
	2    0.500000    -0.500000     0.500000     8.94047722
	3    0.500000     0.500000     0.500000     8.94047722
	4   -0.500000     0.500000     0.500000     8.94047722
	5    0.500000    -0.500000    -0.500000     8.94047722
	6   -0.500000     0.500000    -0.500000     8.94047722
	7   -0.500000    -0.500000    -0.500000     8.94047722
	8   -0.500000    -0.500000     0.500000     8.94047722

+ Faces of the WS cell (6)
# Face ID: vertexID1 vertexID2 ...
   1: 3  4  8  2
   2: 3  4  6  1
   3: 3  1  5  2
   4: 1  6  7  5
   5: 2  8  7  5
   6: 4  8  7  6

+ Lattice vectors for the Wigner-Seitz neighbors
# FaceID: Voronoi lattice vector (cryst. coords.)
   1:  0  0  1
   2:  0  1  0
   3:  1  0  0
   4:  0  0 -1
   5:  0 -1  0
   6: -1  0  0

+ Lattice vectors for the Delaunay reduced cell (cryst. coords.)
  a:  1  0  0
  b:  0  1  0
  c:  0  0  1
  Delaunay reduced cell lengths: 10.323574 10.323574 10.323574
  Delaunay reduced cell angles: 90.000 90.000 90.000

+ Is the cell orthogonal? T
+ Is the reduced cell orthogonal? T
~~~

Critic2 always has a "reference" scalar field defined. The reference
field is the primary target for most keywords. For instance, it
provides the attraction basins integrated when calculating atomic
charges, and it is the field whose critical points are determined by
[AUTO](/critic2/manual/cpsearch/#c2-auto).  In absence of any external
field loaded by the user, critic2 defaults to using the promolecular
density (the sum of atomic densities) as reference.  The promolecular
density is made available to the user through the field identifier $0
(also, $rho0).
~~~
* List of scalar fields
+ Field number 0
  Name: <promolecular>
  Source: <generated>
  Type: promolecular
  Atoms in the environment: 3056
  Use core densities? F
  Numerical derivatives? F
  Nuclear CP signature: -3
  Number of non-equivalent critical points: 2
  Number of critical points in the unit cell: 12
  Alias for this field (2): $0, $rho0
  This is the REFERENCE field.
~~~

A list of the current integrable properties is
given next. This is the list of properties that would be integrated in the
attraction basins if the user runs
[INTEGRALS](/critic2/manual/integrate/#c2-integrals) or any of the other
basin integration methods. The list of integrable properties can be
queried and modified with the [INTEGRABLE](/critic2/manual/integrate/#c2-integrable)
keyword. Our output shows the default integrable properties, which are
the atomic volume and the value and Laplacian of the reference field
~~~
* List of integrable properties (3)
#  Id  Type  Field  Name
	1   v        0  Volume
	2  fval      0  Pop
	3  lval      0  Lap
~~~

Next is the list of additional properties to be calculated at the
critical points. These are used at the end of the automatic
CP search with [AUTO](/critic2/manual/cpsearch/#c2-auto), and can
be modified using the
[POINTPROP](/critic2/manual/cpsearch/#c2-pointprop) keyword. In our
example, there are no additional properties.
~~~
* List of additional properties at critical points (0)
~~~

Each scalar field can be augmented with a core contribution. This can
be useful in cases when the source program only provides the valence
density, and is activated with the
[ZPSP](/critic2/manual/crystal/#c2-charge) keyword. Next in the output
is the list of core and pseudopotential charges for all known fields
(in this case, the promolecular density, which has neither):
~~~
* List of core and pseudopotential charges for each field
# id  type   core?  ZPSP
  0  promol   no
~~~

The execution finishes with a report of the warnings found and the
timestamp. It is always a good idea to check for warnings in the
output:
~~~
CRITIC2 ended successfully (0 WARNINGS, 0 COMMENTS)

Elapsed wall time: 0s
Elapsed CPU time: 0s
CRITIC2--2019/1/31, 15:45:19.603
~~~

## Simple Input and Output for a Molecule {#c2-simplemol}

Molecular structures are read in critic2 using the MOLECULE keyword. A
simple input file for a water molecule is:
~~~
MOLECULE
  O 0.000000 0.000000 0.118882
  H 0.000000 0.756653 -0.475529
  H 0.000000 -0.756653 -0.475529
ENDMOLECULE
~~~
Unlike in CRYSTAL, the coordinates in the MOLECULE environment are
Cartesian coordinates. The default units in and after a MOLECULE
keyword are angstrom.

The output starts off with the same header as in CRYSTAL, and then:
~~~
%% MOLECULE
%% O 0.000000 0.000000 0.118882
%% H 0.000000 0.756653 -0.475529
%% H 0.000000 -0.756653 -0.475529
%% ENDMOLECULE
* Molecular structure
  From: <input>
  Encompassing cell dimensions (bohr): 37.794523  40.654257  38.917797
  Encompassing cell dimensions (ang): 20.000000  21.513306  20.594411
  Empirical formula:
	o(1) h(2)
  Number of atoms: 3
  Number of atomic species: 2
  Number of electrons (with zero atomic charge): 10
~~~

The output shows a copy of the input lines (after the "%%" prefix),
and some general information about the structure. Critic2 works under
periodic boundary conditions, even when dealing with molecular
structures. The molecule is placed inside a very large unit cell
to mimic gas-phase conditions but critic2 treats the molecule in
the same way as a crystal, converting the atomic coordinates to
"crystallographic" coordinates inside the supercell. In the output,
critic2 shows the dimension of this cell
in bohr and angstrom, and the number of atoms and electrons in the
molecule. Keywords are available in the
[MOLECULE](/critic2/manual/molecule/#c2-molecule) keyword and
environment for changing the size and shape of the encompassing
cell.

After that comes the list of atomic species (same as in a crystal) and
the list of atoms in Cartesian coordinates (angstrom and referred to
the same origin as in the input):
~~~
+ List of atomic species:
# spc  Z   name    Q   ZPSP
   1   8     o     0.0  --
   2   1     h     0.0  --

+ List of atoms in Cartesian coordinates (ang_):
# at         x                y                z         spc  name    Z     dnn
   1     0.0000000000     0.0000000000     0.1188820000   1     o     8    0.9622
   2     0.0000000000     0.7566530000    -0.4755290000   2     h     1    0.9622
   3     0.0000000000    -0.7566530000    -0.4755290000   2     h     1    0.9622
~~~
The list of atoms in crystallographic coordinates is not given when
the structure is a molecule. Likewise, symmetry is not used in a
molecular system and hence there is no need for a list of atoms in the
asymmetric unit. All atoms in a molecule have multiplicity 1.

The molecule is placed inside a big cell that
tries to model the empty space around the molecule. This, however, may
lead to some problems with critic2's methods. For instance, the
critical point search will find that at the border of the supercell
the density is discontinuous (because critic2 uses periodic boundary
conditions) and report spurious CPs. Likewise, the gradient path
tracing routines can become trapped at the border of the cell. To
prevent this, critic2 defines by default a second cell, slightly
smaller than the encompassing cell defined above. This
[molecular cell](/critic2/manual/molecule/#c2-molcell) represents the valid
molecular space for the current structure. Regions outside the molecular cell
cannot be traversed by gradient paths and can not hold any critical
points. Essentially, the outside border of the encompassing cell becomes a
representation of infinity for the molecule under study.

The dimensions of the molecular cell are given next in the output:
~~~
+ Limits of the molecular cell (in fractions of the unit cell).
# The part of the unit cell outside the molecular cell represents
# infinity (no CPs or gradient paths in it).
  x-axis: 0.1000 -> 0.9000
  y-axis: 0.0930 -> 0.9070
  z-axis: 0.0971 -> 0.9029
~~~
where the limits are given in fractional coordinates of the
encompassing cell. That is, the molecular cell goes from 0.1 to 0.9 of
the encompassing cell (given above) in the x direction, etc. The
remaining 10% of the cell in each direction becomes the forbidden zone
for this structure.

After this, the output is very similar to CRYSTAL, except the output
related to the crystal symmetry is not present. The atomic
environments and list of molecular fragments are
shown next:
~~~
+ List of fragments in the system (1)
# Id = fragment ID. nat = number of atoms in fragment. C-o-m = center of mass (ang_).
# Discrete = is this fragment finite?
# Id  nat           Center of mass            Discrete
  1    3      0.000000    0.000000    0.052368   Yes

+ Atomic environment
  Number of atoms (reduced cell/environment): 3 / 3
  Radius of (unit cell/environment) circumscribed sphere (ang_): 17.9371 / 0.8129
  Maximum interaction distance (ang_): 11.3820
  Covering regions:
	Total number of regions: 4 (1 2 2)
	Minimum region ID: 0 -1 -1
	Maximum region ID: 0 0 0
	Region side (ang_): 2.6400
	Transformation origin (ang_): 10.0000,10.7567,10.2972
	Search offsets: 1331
	Maximum search offset: 5
	Average number of atoms per region: 0.7500
	Maximum number of atoms in a region: 1

~~~
As in the case of a crystal, critic2 calculates how many discrete
fragments there is in this molecule. In this case, just one. The
environment comprises the whole molecule (3 atoms).

The rest of the output is completely equivalent to the crystal case
(see the discussion above):
~~~
* List of scalar fields
+ Field number 0
  Name: <promolecular>
  Source: <generated>
  Type: promolecular
  Atoms in the environment: 3
  Use core densities? F
  Numerical derivatives? F
  Nuclear CP signature: -3
  Number of non-equivalent critical points: 3
  Number of critical points in the unit cell: 3
  Alias for this field (2): $0, $rho0
  This is the REFERENCE field.

* List of integrable properties (3)
#  Id  Type  Field  Name
	1   v        0  Volume
	2  fval      0  Pop
	3  lval      0  Lap

* List of additional properties at critical points (0)

* List of core and pseudopotential charges for each field
# id  type   core?  ZPSP
  0  promol   no

CRITIC2 ended successfully (0 WARNINGS, 0 COMMENTS)

Elapsed wall time: 0s
Elapsed CPU time: 0s
CRITIC2--2019/1/31, 16:10:37.890
~~~
