
## Introduction

Critic2 can be used as a converter between different crystal structure
file formats. For instance, to convert a cif file to a QE input, we
can use:
~~~
CRYSTAL myfile.cif
WRITE myfile.scf.in
~~~
Sometimes, it is also necessary to create a finite representation of a
crystal by taking the crystal motif, perhaps extended with some of
atoms in the neighboring unit cells. The simplest way of doing this is
by writing an xyz file:
~~~
CRYSTAL myfile_DEN
WRITE myfile.xyz
~~~
The MOLMOTIF keyword writes all atoms in the unit cell and completes
the molecules by using atoms in the neighboring cells. One of the
important aims of WRITE is to write **template** input files for
different programs. The particular keywords used in those templates
(calculation level, basis set, etc.) are mostly meaningless but the
structure is correct. It is up to the user to adapt these templates
to suit their needs.

In critic2, the main keyword to export the molecular or crystal
structure, or a transformation of those structures, is the WRITE
keyword.

## Writing the Structure (WRITE) {#c2-write}
~~~
WRITE file.{xyz,gjf,cml} [ix.i iy.i iz.i] [BORDER]
      [SPHERE rad.r [x0.r y0.r z0.r]] [CUBE side.r [x0.r y0.r z0.r]]
      [MOLMOTIF] [ONEMOTIF] [ENVIRON dist.r] [NMER nmer.i]
WRITE file.{obj,ply,off} [ix.i iy.i iz.i] [BORDER]
      [SPHERE rad.r [x0.r y0.r z0.r]] [CUBE side.r [x0.r y0.r z0.r]]
      [MOLMOTIF] [ONEMOTIF] [CELL] [MOLCELL]
WRITE file.scf.in [rklength.r]
WRITE file.tess
WRITE file.cri|file.incritic
WRITE {[file.]POSCAR|[file.]CONTCAR}
WRITE file.abin
WRITE file.elk
WRITE file.gau
WRITE file.cif [NOSYM|NOSYMM]
WRITE file.d12 [NOSYM|NOSYMM]
WRITE file.m
WRITE file.db
WRITE file.gin
WRITE file.lammps
WRITE file.fdf
WRITE file.STRUCT_IN
WRITE file.hsd
WRITE file.gen
WRITE file.pyscf
WRITE file.in [rklength.r] ## FHIaims geometry.in
WRITE file.frac
~~~
The WRITE keyword writes the currently loaded structure to a file.
A number of file formats can be written by critic2. As in CRYSTAL and
MOLECULE, the type of file is detected by the extension (`.xyz`,
`.in`, `.cri`, etc.).

### Molecular File Formats (xyz, gjf, cml)

With this command:
~~~
WRITE file.{xyz,gjf,cml} [ix.i iy.i iz.i] [BORDER]
      [SPHERE rad.r [x0.r y0.r z0.r]] [CUBE side.r [x0.r y0.r z0.r]]
      [MOLMOTIF] [ONEMOTIF] [ENVIRON dist.r] [NMER nmer.i]
~~~
WRITE generates an xyz file containing a finite
piece of the crystal (if the structure was loaded with CRYSTAL) or the
molecule (resp. MOLECULE). Alternatively, if the `.gjf`
extension is used, a template for a Gaussian input file is
written. If cml is used, a Chemical Markup Language file (xml-style)
is created, containing the same molecular fragment (see below). The number of
cells used in each direction is given by `ix.i`,
`iy.i`, and `iz.i` (default: 1, 1, 1). For the purpose of its graphical
representation, it is sometimes convenient to include atoms that are
almost exactly at the edge of the cell. For instance, the NaCl crystal
is:
~~~
CRYSTAL
  SPG f m -3 m
  CELL 5.64 5.64 5.64 90 90 90 ANG
  NEQ 0. 0. 0. na
  NEQ 1/2 1/2 1/2 cl
ENDCRYSTAL
WRITE nacl.xyz
~~~
Critic2 will (correctly) generate a list of 4 Na and 4 Cl atoms,
representing 1/8th of the conventional cell, because the atoms at (1,
0, 0), (1, 1/2, 0), etc. are repetitions of the atoms in the main
cell. However, this does not look good when the unit cell is
represented because many of the atoms in the cubic cell
are "missing". The BORDER keyword instructs critic2 to include atoms
that are (almost) exactly at the edge of the cell.

The SPHERE keyword writes all atoms inside a sphere of radius `rad.r`
(bohr) and centered around the crystallographic coordinates (`x0.r`,
`y0.r`, `z0.r`). In molecules, the default units for both the center
and radius of the sphere are Cartesian in angstrom. If no center is
given, (0,0,0) is used in both cases. The similar keyword CUBE writes
all atoms inside a cube of side side.r centered around (`x0.r`,
`y0.r`, `z0.r`) (default: (0,0,0)).

The keyword MOLMOTIF is used in molecular crystals. All atoms in the
requested crystal fragment (indicated by the optional
`ix.i`,... integers) are written to the xyz file. Then, the molecules
in the fragment are completed by including atoms from outside the
fragment. Critic2 detects whether the atomic connectivity using a
distance criterion (see the
[BONDFACTOR](/critic2/manual/misc/#c2-bondfactor) and
[RADII](/critic2/manual/misc/#c2-radii) keywords).

The ONEMOTIF and ENVIRON keywords are also used in molecular
crystals. ONEMOTIF writes all atoms in the unit cell, translated by
lattice vectors so that the resulting fragment has whole molecules
only. ENVIRON writes the molecular environment of the unit cell origin
up to a distance `renv.r` (bohr in crystals, angstrom in
molecules). All molecules whose center of mass is at a distance less
than renv.r are written in their entirety, even if some of their atoms
exceed the `renv.r` distance from the origin.

The NMER keyword is used in molecular systems as well, and in
combination with ONEMOTIF, MOLMOTIF, or ENVIRON. When NMER is given,
the fragment of the system selected with either of those three
keywords is split into its component molecules. Then, all monomers,
dimers, trimers,... are written to separate files. All n-mers are
written from monomers up to n-mers, where n is equal to `nmer.i`. In
NMER is used with ENVIRON, the first molecule in all n-mers except for
those with n = `nmer.i` is always part of the Wigner-Seitz cell, which
is useful when generating molecular environments of the crystal for
calculations using incremental methods.

There is an important application of the xyz-format WRITE keyword: the
coordinates written to the xyz file are consistent with the
transformation to Cartesian coordinates in critic2, so it is possible
to bring back all or part of these coordinates to critic2 in order to
represent a subset of the atoms in a crystal. This is very useful
when generating fragments for an [NCIPLOT](/critic2/manual/nciplot/)
calculation (see the FRAGMENT keyword) and in some
[LOAD](/critic2/manual/fields/#c2-load) options to obtain the
promolecular density of a subset of the atoms. Whether the contents of
an xyz file are recognized by critic2 as atoms belonging to the
current system or not can be determined using the
[IDENTIFY](/critic2/manual/structure/#c2-identify) keyword.

The
[CML (Chemical Markup Language)](https://en.wikipedia.org/wiki/Chemical_Markup_Language)
format has the same options as the xyz output format. In the CML
format, an XML-style file is written containing the selected crystal
fragment. If the system is a crystal (loaded with the CRYSTAL
keyword), then the cell geometry is written to the CML file as
well. The CML output format is specially tailored for being easy to
read by [avogadro](https://avogadro.cc) and its underlying engine,
[openbabel](http://openbabel.org).

### Graphical File Formats (obj, ply, off)

The following keyword also writes finite molecular representations of
the structure:
~~~
WRITE file.{obj,ply,off} [ix.i iy.i iz.i] [BORDER]
      [SPHERE rad.r [x0.r y0.r z0.r]] [CUBE side.r [x0.r y0.r z0.r]]
      [MOLMOTIF] [ONEMOTIF] [CELL] [MOLCELL]
~~~
In this case, however, the generated files are graphical
representations.

The OBJ output is the Wavefront OBJ format. The OBJ format is a
three-dimensional model representation, that is, it uses vertices and
faces instead of atoms. This file format is understood by many
visualizers such as view3dscene, meshlab, blender, and others. The
keywords have the same meaning as in the xyz format. The additional
CELL keyword instructs critic2 to write a stick representation of the
unit cell. In a molecular structure, the MOLCELL keyword can be used
to represent the
[molecular cell](/critic2/manual/molecule/#c2-molcell). The
similarly popular PLY (polygon file format or Stanford triangle
format) and OFF (Geomview) file formats can be used as well, with the
same options.

### Quantum ESPRESSO (scf.in)

Quantum ESPRESSO inputs can be written using the extension
`.scf.in`. This conversion is especially useful in the case of
low-symmetry crystals (e.g. monoclinic in a non-conventional setting)
where the conversion from other formats, such as CIF, can be
tricky. The QE input generation works by first determining the Bravais
lattice from the symmetry operations. Critic2 uses 'ibrav=0' always,
and writes a `CELL_PARAMETERS` block containing the
crystallographic-to-Cartesian transformation matrix. QE is particular
about how this matrix should written in order for its own symmetry
module to work. If the crystal setting matches any of those covered in
the QE manual, then that particular matrix is used. Otherwise, critic2
uses its own internal `CELL_PARAMETERS` matrix, which may result in
Quantum ESPRESSO failing to recognize the crystal symmetry. By
default, the crystal cell used by critic2 is written to the QE input
template. To reduce the cell to a primitive, use NEWCELL with the
PRIMITIVE keyword before writing the file. The optional parameter to
the QE input writer (`rklength.r`) is the length parameter (rk-length)
that determines the density of the k-point grid. See the
[KPOINTS](/critic2/manual/structure/#c2-kpoints) keyword.

### Tessel (tess)

A tessel input file (extension `.tess`) can be written. Tessel is a
program for graphical representations of crystals, available upon
request.

### Critic2 (cri, incritic)

A critic2 input file can be written using the .cri or .incritic
extensions. The input file contains the description of the molecular
or crystal structure using the CRYSTAL/MOLECULE environment.

### VASP (POSCAR, CONTCAR)

A VASP `POSCAR` (or `CONTCAR`) is generated by using the `POSCAR` or
`CONTCAR` extension or name. The list of atomic types is written to
the critic2 output. This list is necessary to build the corresponding
`POTCAR`. The atoms are always ordered in increasing atomic number.

### abinit (abin)

An abinit input file containing the input structure can be written by
using the `.abin` extension.

### elk (elk)

An elk input template can be written using the `.elk` extension.

### Crystallographic Information File (cif)

A simple cif file is generated if the `.cif` extension is used. If
NOSYM (or NOSYMM) is used, the cif file is written without symmetry
(i.e. in the P1 space group).

### SHELX (res)

The `.res` extension writes a SHELX file in res format. If NOSYM (or
NOSYMM) is used, the cif file is written without symmetry (i.e. in the
P1 space group).

### Crystal (d12)

A template input file for crystal14 and crystal17 can be written with
the extension `.d12`. If the crystal structure is written without
symmetry (by using the `NOSYM` or `NOSYMM` keywords), then the `.d12`
file contains the complete structure specification. Otherwise, the
`.d12` file uses the `EXTERNAL` keyword to specify the geometry, and
an additional `.fort.34` file with the same root is created. This file
contains the symmetry operations as well as the crystal geometry, and
it should be renamed to just `fort.34` before the crystal calculation
is run. Note that the `.d12` file contains only a keyword to run the
geometry test - there is no basis set specification or any of the other
relevant keywords.

Crystal17 (and probably earlier versions) has different tolerance
values for detecting and accepting symmetry operations than
critic2. It is strongly recommenteded that you use
[SYM REFINE](/critic2/manual/crystal/#c2-symm) before writing a `.d12`
file with symmetry in order to recalculcate the atomic positions in
the detected space group.

### Gaussian (periodic, gau)

A Gaussian input file for calculations under periodic boundary
conditions can be written using the `.gau` extension. For a template
corresponding to a finite molecule, use `.gjf` (see above).

### Escher/Octave (m)

The octave script file (extension `.m`) contains the structure in
octave format, to be read using the
[escher library](https://github.com/aoterodelaroza/escher).

### DCP database (db)

The db file format is intended for a set of automated input generation
octave scripts, the [dcp package](https://github.com/aoterodelaroza/dcp).

### Gulp (gin)

A simple GULP template input file containing the structure (and EEM as
the first line) can be written using the `.gin` extension. No resonant
carbon atoms are detected. The GULP writer is experimental, so please
exercise care and double-check the templates.

### LAMMPS (lammps)

For file names with an extension `.lammps`, critic2 writes a simple
LAMMPS data file containing one unit cell (length units are
angstrom). Only orthogonal cells are supported for now. The LAMMPS
writer is experimental, so please exercise care and double-check the
templates.

### SIESTA (fdf, STRUCT_IN)

Two types of siesta inputs can be generated. The `.fdf` extension
writes a template for a proper functional siesta input template
containing the crystal structure. The `STRUCT_IN` extension or name
writes files that can be read using the `MD.UseStructFile` option.

### DFTB+ (gen, hsd)

Two inputs types for DFTB+ may be written. The `.gen` format contains
only the structure and is meant to be used with the `GenFormat` method
in `Geometry`. The `.hsd` writes a full input template, including the
structure.

### pyscf (pyscf)

The `.pyscf` extension writes a template python script for using with
pyscf. The script creates an instance of a Mole object (in the case of
a molecule) or a Cell object (if it is a crystal).

### FHIaims (in)

A `geometry.in` input file for FHIaims can be generated using the
`.in` extension (note that using `scf.in` generates QE inputs). The
file then needs to be renamed to `geometry.in`. The appropriate
selection of keywords is used depending on whether the current
structure is a crystal or a molecule.

If an additional real number is given (`rklength.r`), an additional
file (`file.in_control`) containing a template `control.in` for the
FHIaims calculation is written. The template `control.in` specifies a
uniform k-point grid in the `k_grid` keyword with a density given by
the length parameter `rklength.r`. See the
[KPOINTS](/critic2/manual/structure/#c2-kpoints) keyword for details.

### TINKER frac files (frac)

A `.frac` file can be generated in TINKER format. The `.frac` file
always contains the cell lengths and angles (even if a molecule is
written) and the atomic connectivity calculated by critic2. The force
field atom types in the file correspond to the "tiny" force field
parameters in the TINKER distribution (`tiny.prm`).

## Writing a .mols File for DMACRYS/NEIGHCRYS (MAKEMOLSNC) {#c2-makemolsnc}
~~~
MAKEMOLSNC file_fort.21.s file_mols.s
~~~
The MAKEMOLSNC keyword is used to write .mols files for NEIGHCRYS, the
companion program to
DMACRYS. [DMACRYS](http://www.chem.ucl.ac.uk/cposs/dmacrys/index.html)
is a program for the calculation of lattice energies and their
derivatives using rigid molecules and the distributed multipole
approach. The atomic multipoles are calculated from a gas-phase
molecular wavefunction using the [GDMA program](http://www-stone.ch.cam.ac.uk/pub/gdma/)

The NEIGHCRYS program is used to help in the preparation of input
files for DMACRYS. Sometimes, NEIGHCRYS has difficulty generating the
.mols file, a file that contains the local frame of reference for the
calculation of the multipoles for every molecule in the crystal. (See
also [SYM WHOLEMOLS](/critic2/manual/crystal/#c2-symm).)

The MAKEMOLSNC keyword is intended to be a replacement for this
function in NEIGHCRYS. MAKEMOLSNC reads the NEIGHCRYS output (a
`fort.21` file) and generates the corresponding mols file with name
`file_mols.s`. MAKEMOLSNC does not work if the molecules are linear
because NEIGHCRYS does not work with them either. The use of
MAKEMOLSNC does not require a molecular or crystal structure being
loaded.

