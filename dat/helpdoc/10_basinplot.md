
## Attractor Basin Plots (BASINPLOT) {#c2-basinplot}

~~~
BASINPLOT [CUBE [lvl.i] | TRIANG [lvl.i] | 
  SPHERE [ntheta.i nphi.i]]
  [OFF|OBJ|PLY|BASIN|DBASIN [npts.i]}]
  [CP cp.i] [PREC delta.r] [VERBOSE] [MAP id.s|"expr"]
~~~
The BASINPLOT keyword plots the attraction (Bader) basin of the CP
`cp.i` from the complete list. If CP is not given, all the
non-equivalent attractors are used. BASINPLOT works by tracing rays
starting at the CP in question and then doing bisection until the
limit of the basin is found. At each point in the bisection procedure,
a gradient path is traced. If the gradient path ends in the CP we are
plotting then the point is considered to be inside the
basin. Otherwise, the point is outside.

There are several keywords that control the number and position of the
bisection rays for BASINPLOT. With CUBE, a cube is selected as the
starting polyhedron, and recursively subdivided lvl.i times. The final
(convex) polyhedron is placed on the attractor and the zero-flux
surface limit for the rays corresponding to the polyhedron vertices is
determined. TRIANG follows the same process but starting from an
octahedron. SPHERE sets the rays by doing a triangulation of the unit
sphere with `nphi.r` and `ntheta.i` angles. The total number of points
is given by the formula $$2n_{\phi}*(2^{n_{\theta}}-1)+2$$.

The output keyword selects the output format for the basin plot: OFF
(geomview), OBJ (Wavefront obj), PLY (Standford ply), BASIN (tessel),
and DBASIN. The DBASIN file format also contains information about
scalar fields measured along the basin rays.

The naming scheme for the output files is `\<root\>-cp.ext` where root
is the root of the run (the name of the input file up to the first dot
unless changed by the [ROOT](/critic2/manual/misc/#c2-root)
keyword), `cp` is the complete CP list identifier of the attractor and
`ext` is the extension selected with one of the format keywords
above.

The precision of the bisection is set using the PREC keyword, and
equal to `delta.r`. VERBOSE gives more information in the output about
the bisection process, useful if you are impatient.

If a 3D model format is used (OFF, OBJ, PLY), the MAP keyword can be
used to create a color map of a given field (given by the field number
or identifier `id.s`) or a field-containing expression ("expr") onto
the surface. The color scale limits are the minimum and the maximum
value of the field or expression on all the points of the surface. The
mapping function is the same as in gnuplot: 
($$r =\sqrt{x}$$, $$g=x^3$$, $$b=\sin(360\times x)$$, with $$x$$
ranging from 0 to 1).

The default is the TRIANG method, `lvl.i = 3`, 
`ntheta.i = nphi.i = 5`, OBJ output, and one basin plot for
all the non-equivalent attractors found in AUTO.

## Primary Bundle Plots (BUNDLEPLOT) {#c2-bundleplot}

~~~
BUNDLEPLOT x.r y.r z.r
  [CUBE [lvl.i] | TRIANG [lvl.i] | SPHERE [ntheta.i nphi.i]]
  [OFF|OBJ|PLY|BASIN|DBASIN [npts.i]}]
  [ROOT root.s] [PREC delta.r] [VERBOSE] [MAP id.s|"expr"]
~~~
The BUNDLEPLOT keyword plots a primary bundle starting from a point in
its interior, given by `x.r`, `y.r`, and `z.r` in crystallographic
coordinates (crystal) or molecular Cartesian coordinates (molecule,
default units: angstrom). The syntax of the BUNDLEPLOT keyword is
essentially the same as BASINPLOT.

The same bisection algorithm described in BASINPLOT is used in
BUNDLEPLOT. The precision can be set with the PREC keyword to
`delta.r` (default: 1d-5 bohr). The set of rays to be traced are
obtained by a recursive subdivision (`lvl.i` cycles) of a cube (CUBE),
an octahedron (TRIANG) or using a uniform distribution of `ntheta.i`
times `nphi.i` points on the unit sphere (SPHERE). The output file has
root `root.s` (ROOT keyword), and its format may be OFF, OBJ, PLY,
BASIN or DBASIN with `npts.i` points sampled along each ray.

If a 3D model format is used (OFF, OBJ, PLY), the MAP keyword can be
used to create a color map of a given field (given by the field number
or identifier `id.s`) or a field-containing expression ("expr") onto
the surface. The color scale limits are the minimum and the maximum
value of the field or expression on all the points of the surface. The
mapping function is the same as in gnuplot: 
($$r =\sqrt{x}$$, $$g=x^3$$, $$b=\sin(360\times x)$$, with $$x$$
ranging from 0 to 1).

Default values: TRIANG, `lvl.i = 3`, `ntheta.i = nphi.i = 5`, OBJ output,
`root.s = <root>-bundle`.

