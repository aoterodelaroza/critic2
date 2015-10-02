Tools
-----

charges.sh
~~~~~~~~~~

The script charges.sh (in the src/ directory) can be used to calculate
atomic charges given an abinit, vasp or cube (QE) grid by applying the
YT algorithm. The usage is simply: 

::

  #ascii#
  charges.sh file_DEN     # (abinit)
  charges.sh CHGCAR       # (vasp)
  charges.sh file.cube    # (qe and cubes)

For more info,

:: 

  #ascii#
  charges.sh -h  

Off2off
~~~~~~~

Off2off is a program to transform and combine OFF and COFF
three-dimensional object description files. To run the program, use:

::

  #ascii#
  off2off [OPTIONS] input_file(s)

where input_files are one or more files written in OFF or COFF
format. The output, also written in OFF or COFF format is written
to the standard output unit, as a default.

The options for off2off are:

* -c : give (solid) color to OFF files (default).

* -u : remove color from COFF files.

* -z : translate coordinates origin to the geometric center of files.

* -m : translate coordinates to the mass center of input points.

* -s : scale point (x,y,z) coordinates to [-1,1].

* -b : read and write BASIN files instead of the standard OFF.

* -d : read and write DBASIN files.

* -o output : write results to given file instead to stdout.

* -t matrix : read in a group of transformation matrices from the
  given file.

Basin2off
~~~~~~~~~

Basin2off uses the data contained in one or more BASIN data files and
prepares OFF and COFF files to plot colormap renderings of a scalar
property on the atomic basins surface. It is run:

::

  #ascii#
  basin2off [-c] [-g] [-l] [-a] [-L] [-A] [-n num] input_files....

The list of input_files will be analyzed together to produce a
common color scale. For each input file a COFF or OFF file will
be produced, its name being that of the input with the suffix
'.basin' stripped and substituted by '.coff' or '.off'. The list of
options is:

* -c : used together with -n 0 will produce colored OFF files. Each
  basin input file will be given a different solid color.

* -g : use a simple scale based on gray values, instead of the
  default color scale.

* -l : assign the colors according to a logarithmic scale rather
  than the default linear scale.

* -a : use a atan() mapping function instead.

* -L : experimental: atan(log()) map

* -A : experimental: atan(atan()) map

* -n num : the scalar property used to produce the COFF files will
  be the num-th one contained in the BASIN files. (Default: num =
  1). (Special case: num=0 will produce OFF files with no color
  defined for the basin faces).

Basinmerge
~~~~~~~~~~
Basinmerge is a script to copy by symmetry and combine the
BASIN/DBASIN/OFF/COFF datafiles of several atoms in the crystal. It
uses the tessel, off2off and basin2off programs. The program
interprets a collection of orders/instructions written according to
the next language:

- **FILETYPE [OFF|COFF|BASIN|DBASIN]**

  Type of file to be processed.

- **NEQION name.s file.s x.r y.r z.r**

  Enter basin description for each non-equivalent ion

- **CRYS2CART**

  ::

          c11.r  c12.r  c13.r  c14.r
          c21.r  c22.r  c23.r  c24.r
          c31.r  c32.r  c33.r  c34.r
          c41.r  c42.r  c43.r  c44.r

  Matrix that transforms from crystal to Cartesian coordinates.
  It is written in the critic2 output.

- **CART2CRYST**

  ::

          d11.r  d12.r  d13.r  d14.r
          d21.r  d22.r  d23.r  d24.r
          d31.r  d32.r  d33.r  d34.r
          d41.r  d42.r  d43.r  d44.r

  Matrix that transforms from Cartesian to crystal coordinates.
  It is written in the critic2 output.

- **COPY fromname.s toname.s fileto.s**

  ::

          t11.r  t12.r  t13.r  t14.r
          t21.r  t22.r  t23.r  t24.r
          t31.r  t32.r  t33.r  t34.r
          t41.r  t42.r  t43.r  t44.r

  Create a copy of the basin file of a NEQ ion into an equivalent
  position. The transformation matrix that is entered should act
  on the crystallographic coordinates, and it is produced by
  critic2.

- **MERGETO  destfile.s  file1.s  [ file2.s ... ]**

  Merge a collection of files into a single one.

- **BASIN2OFF  order.s**

  Run directly the basin2off program.

- **OFF2OFF  order.s**

  Run directly the off2off program.

Scripts
~~~~~~~

For now, just one script (in tools/scripts), cif2scf.sh, that converts
from a cif file to a quantum espresso input. To use it, 

::

  #ascii#
  cif2scf.sh file.cif 

which creates file.scf.in using critic2.

Elk mod for the ELF and Coulomb potential
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The tools/elk_mode directory contains three source files to modify elk
(1.3.2) and make it print a description of the ELF and the Coulomb
potential. See the README in the subdirectory.

