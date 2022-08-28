
## Points (POINT) {#c2-point}

~~~
POINT [x.r y.r z.r|file.s] [ALL] [FIELD {id.s|"expr.s"}]
~~~
Calculates the value of the reference field or an arithmetic
expression at point (`x.r`, `y.r`, `z.r`) in crystallographic
coordinates (if the structure is a CRYSTAL) or molecular Cartesian
coordinates (if a MOLECULE). For the latter, the default units
are angstrom unless changed by the
[UNITS](/critic2/manual/inputoutput/#c2-units) keyword. If a file name
is passed instead (`file.s`), calculate the same quantities at all the
points specified by the file. All non-blank lines in the file that do
not start with a comment symbol (`#`) represent a point, and only the
first three real numbers from each line are read.

If ALL is used, all loaded fields are evaluated. In addition, all
arithmetic expressions that have been registered using the
[POINTPROP](/critic2/manual/cpsearch/#c2-pointprop) keyword are also
calculated. The POINTPROP keyword combined with POINT is useful to
evaluate chemical functions at arbitrary points in space.

If FIELD is used and followed by an integer or field identifier
(`id.s`), then only that field is evaluated. FIELD followed by an
arithmetic expression calculates the value of that expression at the
point.

## Lines (LINE) {#c2-line}

~~~
LINE x0.r y0.r z0.r x1.r y1.r z1.r npts.i [FILE file.s]
     [FIELD id.s|"expr.s"] [GX|GY|GZ|GMOD|HXX|HXY|HXZ|HYX|HYY|
     HYZ|HZX|HZY|HZZ|LAP]
~~~
Calculate a line from (`x0.r`, `y0.r`, `z0.r`) to (`x1.r`, `y1.r`,
`z1.r`) with `npts.i` points. The units for the two endpoints (x0 and
x1) are crystallographic coordinates in crystals and molecular Cartesian
coordinates in molecules. The latter are angstrom by default (unless
[UNITS](/critic2/manual/inputoutput/#c2-units) is used).

By default, the result is written to the standard output, but it can
be redirected to a file using FILE. The reference field is used unless
a FIELD keyword is given, in which case the field `id.s` or the
expression `expr.s` are evaluated. Together with the value of the field,
an additional quantity can be evaluated: the components of the
gradient (GX,GY,GZ), the norm of the gradient (GMOD), the components
of the Hessian (HXX,...) and the Laplacian of the reference (or the
`id.s`) field.

## Planes and Contour Plots (PLANE) {#c2-plane}

~~~
PLANE x0.r y0.r z0.r x1.r y1.r z1.r x2.r y2.r z2.r nx.i ny.i
      [SCALE sx.r sy.r] [EXTENDX zx0.r zx1.r] [EXTENDY zy0.r zy1.r]
      [FILE file.s] [FIELD id.s/"expr"]
      [F,GX,GY,GZ,GMOD,HXX,HXY,HXZ,HYY,HYZ,HZZ,LAP]
      [CONTOUR {LOG niso.i [zmin.r zmax.r]|ATAN niso.i [zmin.r zmax.r]|
      BADER|LIN niso.i [rini.r rend.r]|i1.r i2.r ...}] [COLORMAP [LOG|ATAN]]
      [RELIEF zmin.r zmax.r] [LABELZ labelz.r]
~~~
Calculate the value (or derivatives) of the reference field on a
plane. The results are written to a file, with default name
`<root>_plane.dat`. The geometry of the plane is specified by three
points: the origin (`x0.r`, `y0.r`, `z0.r`), the end of the x-axis
(`x1.r`, `y1.r`, `z1.r`), and the end of the y-axis (`x2.r`, `y2.r`,
`z2.r`). The number of points calculated on each axis are given by
`nx.i` (x-axis) and `ny.i` (y-axis). The units are
crystallographic coordinates in a crystal, and molecular Cartesian
coordinates in a molecule (default: angstrom unless the
[UNITS](/critic2/manual/inputoutput/#c2-units) keyword
is used). The two axes of
the plane can be scaled using the SCALE keyword. If `sx.r` (`sy.r`) is
given, the total length of the x-axis (y-axis) is scaled by `sx.r`
(`sy.r`). If EXTENDX is used, extend the x-axis by `zx0.r` (initial point
of the x-axis) and `zx1.r` (end point). The keyword EXTENDY performs the
equivalent operation on the y-axis. The units for EXTENDX and EXTENDY
are bohr (crystals) or angstrom (molecules) unless changed by the
[UNITS](/critic2/manual/inputoutput/#c2-units) keyword.

The name of the output file can be changed with FILE. Using FIELD, one
of the loaded fields (`id.s`) or an expression (`expr.s`) can be
evaluated. In addition to the field value, a second property can be
evaluated: the field again (F), its derivatives (Gx), its second
derivatives (Hxx), the gradient norm (GMOD) or the Laplacian (LAP).

The keyword CONTOUR writes a contour map representation of the plane:
two contour line files (`.iso` and `.neg.iso`) and a gnuplot script
(`.gnu`). The distribution of contour values can be: logarithmic (LOG,
with `niso.i` contours), arctangent (ATAN, with `niso.i` contours),
same as in the aimpac program (BADER, {1,2,4,8}x10^{-3,-2,-1,0,1}),
linear (LIN, `niso.i` contours from `r0.r` to `r1.r`), or the user can
specify the contour values manually (no keyword). In LOG, ATAN, and
LIN, the default contour values range from the minimum to the maximum
value of the field in the plot. These quantities can be changed by
passing the optional zmin.r and zmax.r parameters to LOG/ATAN. The
field or any of its derivatives, selected with the [F|GX|...] keyword,
is used for the contour plot.  The
[GRDVEC](/critic2/manual/gradientpath/#c2-grdvec) keyword performs the
same task as PLANE with the CONTOUR option, and more (e.g. tracing
gradient paths), but is more complex to use.

The RELIEF keyword writes a gnuplot template for a three-dimensional
relief plot using the data calculated by PLANE. The default suffix is
`-relief.gnu`. The mandatory arguments `zmin.r` and `zmax.r` set the
range of the z axis in the plot.

The COLORMAP keyword writes a template for a colormap plot of the
field on the plane. If the LOG or ATAN keywords are given, the
logarithm or the arctangent of the field are represented in the
colormap.

For the plots that display atomic or critical point labels, LABELZ
controls how many labels are represented. Any atom or critical point
that is at a distance less than `labelz.r` (default: 0.1 bohr) is
shown as a label in the plot.

## Grids (CUBE) {#c2-cube}

~~~
CUBE x0.r y0.r z0.r x1.r y1.r z1.r nx.i ny.i nz.i [FILE file.s] [FIELD id.s/"expr"]
     [F,GX,GY,GZ,GMOD,HXX,HXY,HXZ,HYY,HYZ,HZZ,LAP] [HEADER] [ORTHO]
CUBE x0.r y0.r z0.r x1.r y1.r z1.r lpp.r ...
CUBE CELL {lpp.r|nx.i ny.i nz.i} ...
CUBE GRID [SHIFT ix.i iy.i iz.i] ...
CUBE MLWF ibnd.i nRx.i nRy.i nRz.i [SPIN ispin.i] ...
CUBE WANNIER ibnd.i nRx.i nRy.i nRz.i [SPIN ispin.i] ...
CUBE UNK ibnd.i ik.i [SPIN ispin.i] ...
CUBE PSINK ibnd.i ik.i nRx.i nRy.i nRz.i [SPIN ispin.i] ...
CUBE ... FILE CHGCAR
CUBE ... FILE bleh.cube
CUBE ... FILE bleh.bincube
CUBE ... FILE bleh.xsf
~~~
The CUBE keyword writes a three-dimensional grid in Gaussian cube,
binary cube, VASP CHGCAR, and xsf formats. The limits of the grid can
be set in three ways. By giving the end-points (`x0.r`, `y0.r`,
`z0.r`) and (`x1.r`, `y1.r`, `z1.r`) it is possible to build a grid
from a fragment of the system (this is only possible using
the cube and bincube formats). In crystals, this fragment has the same
shape as the unit cell. To write an orthogonal fragment, use the
additional ORTHO keyword. In molecules, it is always an orthogonal
fragment. The CELL keyword calculates a grid spanning the entire unit
cell. GRID, MLWF, WANNIER, UNK, and PSIK has the same effect as
CELL regarding the output grid geometry.

If the end-points are given, they must be in crystallographic
coordinates if the system is a periodic crystal (the structure was
read using the CRYSTAL keyword) or molecular Cartesian coordinates if
the system is a molecule (read with the MOLECULE keyword). The units
in the latter default to angstrom unless changed using the
[UNITS](/critic2/manual/inputoutput/#c2-units) keyword.

The number of points in the grid can also be controlled in several
ways. If the grid limits are given explicitly or using CELL, then the
number of points on each axis can be indicated by giving three
integers (`nx.i`, `ny.i`, and `nz.i`) corresponding to the number of
points in the x-, y-, and z-axis respectively. If a single number
(`lpp.r`) is found, then the number of points is the length of the
axis divided by `lpp.r` (lpp means "length per point"). The 
units for `lpp.r` are angstrom for molecules and bohr for crystals
unless the [UNITS](/critic2/manual/inputoutput/#c2-units) keyword is
used.

The GRID keyword can be used to write a given grid to a grid file
directly. If FIELD is used in combination with GRID, then the
indicated field or expression is used; otherwise, the reference field
is used. The GRID keyword is useful when combined with LOAD and
arithmetic operations to read, manipulate, and then save grids to an
external file. If GRID is used, both the geometry of the grid and the
number of points are taken from the parent grid field. The SHIFT
keyword is used for shifting the origin of the grid to a different
point. The origin of the shifted grid is the position of point `ix.i`,
`iy.i`, `iz.i` in the old grid.

The MLWF, WANNIER, UNK, and PSINK keywords are similar to GRID in that
they dump a scalar field on a grid to a file directly. These keywords
only work with [Quantum ESPRESSO pwc files](/critic2/manual/fields/c2-qepwc),
which contain the information about the Bloch states in a periodic
solid. The meaning of these keywords is:

- MLWF writes the maximally-localized Wannier function for band
  `ibnd.i` and lattice vector given by the integers `nRx.i`, `nRy.i`,
  and `nRz.i` ($$w_{nR}({\mathbf r})$$). Requires the Wannier
  checkpoint file for the Bloch coefficient rotation.

- WANNIER writes a Wannier function calculated without rotation of the
  Bloch coefficients, for band `ibnd.i` and lattice vector given by
  the integers `nRx.i`, `nRy.i`, and `nRz.i` ($$w_{nR}({\mathbf
  r})$$). Requires that the k-point list corresponds to a uniform
  Monkhorst-Pack grid (i.e. no symmetry).

- UNK writes the periodic part of the Bloch state with band index
  `ibnd.i` and k-point `ik.i` ($$u_{nk}({\mathbf r})$$). The k-point
  identifier can be found from the output when the `.pwc` file is
  loaded. Does not require a Wannier checkpoint file.

- PSINK writes the Bloch state with band index `ibnd.i` and k-point
  `ik.i` at lattice vector `nRx.i`, `nRy.i`, and `nRz.i`
  ($$\psi_{nk}({\mathbf r}-{\mathbf R}) = u_{nk}({\mathbf r}) e^{i{\mathbf k} ({\mathbf r}-{\mathbf R})}$$).
  The k-point identifier can be found from the output when the `.pwc`
  file is loaded. Does not require a Wannier checkpoint file.

In a spin-polarized calculation, the spin channel can be selected
using the SPIN keyword (1 for spin-up and 2 for spin-down).

Independently on how the grid is set up, several options control the
behavior of CUBE. FILE sets the name of the output file (default:
`<root>.cube`). critic2 uses different formats depending on the
extension: binary cube (`.bincube`), cube (`.cube`), xsf (`.xsf`), and
VASP-style CHGCAR (everything else). The HEADER keyword writes only
the header to the output file. FIELD sets the field to be used by the
field number or alias (`id.s`). Alternatively, an arithmetic
expression can be used. Finally, a derivative of the scalar field
(gradient, Hessian, Laplacian) can be selected instead of the value of
the field itself (F) to build the grid.
