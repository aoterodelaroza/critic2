
## Overview

Gradient paths are the solution of the differential equation
$$x^\prime = {\mathbf \nabla}f({\mathbf x})$$ where $$f({\mathbf x})$$
is a scalar field. They play an important role in QTAIM theory because
gradient paths of the electron density cannot cross the boundary
between atomic regions. In consequence, a gradient path plot is a
simple way to investigate the shape and properties of a
basin. Gradient paths originate at maxima (if the field is the density
they are usually the nuclei) and end at the minima (the crystal voids)
or at infinity in case of a gas-phase molecule.

There are two keywords for making gradient path representations in
critic2: GRDVEC (2D gradient paths plus a contour plot) and FLUXPRINT
(3D).

## Gradient Path Representations in a Plane (GRDVEC) {#c2-grdvec}

~~~
GRDVEC
   {FILES|ROOT|ONAME} rootname.s
   PLANE x0.r y0.r z0.r x1.r y1.r z1.r x2.r y2.r z2.r
   SCALE sx.r sy.r
   EXTENDX zx0.r zx1.r
   EXTENDY zy0.r zy1.r
   OUTCP sx.r sy.r
   HMAX hmax.r
   ORIG x.r y.r z.r atr.i up.i down.i
   CP cp.i up.i down.i
   CPALL
   BCPALL up.i down.i
   RBCPALL bup.i bdown.i rup.i rdown.i
   CHECK
        x.r y.r z.r
        ...
   ENDCHECK/END
   CONTOUR {F,GX,GY,GZ,GMOD,HXX,HXY,HXZ,HYY,HYZ,HZZ,LAP} 
     nptsu.i nptsv.i {LIN niso.i [cini.r cend.r]|
     LOG niso.i [zmin.r zmax.r]|ATAN niso.i [zmin.r zmax.r]|
     BADER|i1.i i2.i...}
ENDGRDVEC/END
~~~
GRDVEC makes a plot in a plane. The plot contains the gradient
paths originating from a set of points, critical or otherwise, inside
that plane, in addition to a contour representation of a scalar
field. The GRDVEC syntax consists of a GRDVEC...ENDGRDVEC environment
that accepts a set of input lines (in any order) that control the
characteristics of the plot. The syntax of GRDVEC originated from (and
is similar to) Bader's AIMPAC suite of programs.

By using the FILES keyword (equivalently, ROOT or ONAME), the user
sets the root name of the output files containing the information for
the plot (default: \<root\>). These files include:

* \<root\>.grd : gradient path data.

* \<root\>.dat : values of the reference field on the plane.

* \<root\>.iso, \<root\>.neg.iso : positive (green) and negative (blue)
  contour lines. 

* \<root\>.gnu : gnuplot script file that generates the merged
  gradient/contour plot.

* \<root\>-label.gnu : gnuplot script file loaded in \<root\>.gnu
  containing the information for the position of the CPs in the plot
  plane.

PLANE specifies the plane for the plot using three points: `x0` is the
origin, `x1` the end of the x-axis and `x2` the end of the y-axis. By
default, these points are in crystallographic coordinates in a
crystal, and in molecular Cartesian coordinates in a molecule (default
units: angstrom). The two axes of the plane can be scaled using the
SCALE keyword. If `sx.r` (`sy.r`) is given, the total length of the
x-axis (y-axis) is scaled by `sx.r` (`sy.r`). If EXTENDX is used,
extend the x-axis by `zx0.r` (initial point of the x-axis) and `zx1.r`
(end point). The keyword EXTENDY performs the equivalent operation on
the y-axis. The units for EXTENDX and EXTENDY are bohr (crystals) or
angstrom (molecules) unless changed by the 
[UNITS](/critic2/manual/inputoutput/#c2-units) keyword.

The plot plane may contain regions that are traversed by gradient
lines originating at critical points located inside the plane but
outside the plot region. If this is the case, the OUTCP option allows
the user to extend the plane for the CP labels. The `sx.r` and `sy.r`
in OUTCP are scale parameters for the plane, same as in SCALE, but
only apply to CP labels. The x-axis extends $$(s_x-1)\times l_x$$ in
each direction, where $$s_x$$ is `sx.r` and $$l_x$$ is the length of
the x-axis. The `sy.r` variable works the same way. The plane
determined by the vectors given in PLANE acts as a clipping plane
while the scaled plane determines the gradient path origins.

With HMAX, you can set the maximum distance from a CP to the plane
to be included in the plot (units: bohr in crystals, angstrom in
molecules). Default: 1d-4 bohr. 

The ORIG keyword adds a source of gradient lines to the plot. Its
coordinates are `x.r`, `y.r` and `z.r`. Crystallographic coordinates
are used in a crystal, and molecular Cartesian coordinates in a
molecule (default units: angstrom). `atr.i` is 1 if the point is to be
treated as a ncp or ccp (the up and down trajectories start from
points located on a sphere centered on the origin) and it is 0 if the
point is to be treated as a bcp or ccp (a circle is built around the
CP in the plane determined by two eigenvectors whose eigenvalues have
equal sign. The remaining eigenvector determines a unique
direction). `up.i` and `down.i` are the number of gradient paths to be
started in the upwards and downwards direction, respectively.

The CP keyword accepts a critical point identifier from the complete
CP list (or the complete atom list). The number of upwards and
downwards gradient paths must be given. A special case is the CPALL
keyword, which adds as origins every critical point in the CP list on
the selected plane. The default number of gradient paths is 36 down
for ncps and 36 up for ccps, and 2 up and 2 down for bcps and
rcps. The BCPALL keyword is similar to CPALL, except that only the
bond critical points are included as origins. If BCPALL is used, the
user must supply the number of gradient lines in the upwards and
downwards directions. In a similar way, RBCPALL includes bond and ring
critical points, and the user must give the number of upwards and
downwards gradient paths for bonds (`bup.i`, `bown.i`) and rings
(`rup.i`, `rdown.i`).

The CHECK environment allows the user to enter the crystallographic
coordinates of a CP of the scalar field to add it as an origin. If the
point given is not a CP or if it lies outside the selected plane, it
is excluded from the list of points that are sources of gradient
paths. The valid CPs in the CHECK list are identified and an adequate
number of gradient paths are started according to its character: for a
ncp and ccp, 36 upwards or downwards and for a bcp or rcp, 2 upwards
and 2 downwards.

The CONTOUR keyword makes critic2 generate a plot in which the
gradient paths calculated in GRDVEC are merged with a contour plot, in
the spirit of the CONTOUR option to
[PLANE](/critic2/manual/graphics/#c2-plane). The scalar field for the
contour plot can be selected with 
F (current reference field), GX, GY,... The syntax for the derivatives
is the same as in PLANE. After this, the user must specify the
number of points in each direction of the plane (`nptsu.i` and
`nptsv.i`) and the contour values and type of mapping. The contour
distribution can be: logarithmic (LOG, with `niso.i` contours),
arctangent (ATAN, with `niso.i` contours), same as in the aimpac
program (BADER, {1,2,4,8}x10^{-3,-2,-1,0,1}), linear (LIN, `niso.i`
contours from `r0.r` to `r1.r`), or the user can specify the contour
values manually (no keyword). In LOG and ATAN, the default contours
range from the minimum to the maximum value of the field in the
plot. These quantities can be changed by passing the optional `zmin.r`
and `zmax.r` parameters to LOG/ATAN.

Note that GRDVEC is able to handle non-orthogonal axes for the plot
plane. If the two plane axes determined in the PLANE keyword are
non-orthogonal, the final graph will correctly reflect the actual
appearance of the plane by conserving the original angle between the
x- and y- axis. Also, note that at most 2 gradient lines may be traced
from bcps and rcps, either upwards or downwards. Thus, for example,
`BCPALL 2 2` is equivalent to `BCPALL 2 100` or `BCPALL 100 100`.

## Three-Dimensional Gradient Path Representations (FLUXPRINT) {#c2-fluxprint}

~~~
FLUXPRINT
  POINT {1|-1|0} x.r y.r z.r [step.r epsi.r]
  NCP cp.i ntheta.i nphi.i [step.r epsi.r]
     [LVEC x.i y.i z.i]
  BCP cp.i 1 [step.r epsi.r] [LVEC x.i y.i z.i]
  BCP cp.i {0|-1} n.i [step.r epsi.r]
     [LVEC x.i y.i z.i]
     [BRAINDEAD|QUOTIENT|DYNAMICAL]
  RCP cp.i -1 [step.r epsi.r] [LVEC x.i y.i z.i]
  RCP cp.i {0|1} n.i [step.r epsi.r]
     [LVEC x.i y.i z.i]
     [BRAINDEAD|QUOTIENT|DYNAMICAL]
  CCP cp.i ntheta.i nphi.i [step.r epsi.r]
     [LVEC x.i y.i z.i]
  GRAPH igraph.i [step.r epsi.r]
  COLOR r.i g.i b.i
  TEXT|TESSEL|TESS|OBJ|PLY|OFF|CML
  SHELLS ishl.i
  NOSYM
ENDFLUXPRINT/END
~~~

The FLUXPRINT keyword prints three-dimensional gradient paths. There
are several plotting commands:

- POINT: plot a gradient path starting at point (`x.r` `y.r` `z.r`)
  in crystallographic coordinates (crystals) or in molecular Cartesian
  coordinates (molecules, default unit: angstrom). `step.r` is the
  maximum step (Cartesian coordinates) for the gradient path tracing
  algorithm. If `step.r` > 0, use it also as the initial step. If
  `step.r` < 0, use a small step as initial the initial step
  (1d-3). `epsi.r` is the  gradient norm stop criterion. The default
  values are 0.1 for step and 1e-9 for `epsi.r`. The {1|-1|0} field
  controls the direction of the path. An ascending gradient path is
  obtained with 1 while -1 issues a descending path. 0 = -1 + 1 makes
  FLUXPRINT represent both ascending and descending paths.

- NCP: print gradient paths starting from a (small) sphere centered on
  the nuclear CP identified by `cp.i` (this identifier comes from the
  complete CP list). The number of points is controlled by `ntheta.i`
  (number of points sampling the azimuthal angle) and `nphi.i` (number
  of points sampling the polar angle). `cp.i` specifies a ncp in the
  main cell up to a lattice translation. The LVEC optional keyword
  allows the user to enter a lattice vector to displace the
  represented gradient paths from their initial position, which is
  given by the complete CP list written by AUTO.

- BCP: print gradient paths starting at the vicinity of a bond CP,
  identified by `cp.i`. If the gradient path is ascending (1 in the
  fourth field), the (unique) bond path associated to the bcp is
  represented. If -1 is given instead, the IAS associated to the bcp
  is sampled starting from a small circle surrounding the bcp, with
  `n.i` points on it. With a 0 value, both tasks are performed.

  The three keywords BRAINDEAD, QUOTIENT and DYNAMICAL establish the
  method employed in generating the starting angular grid. With
  BRAINDEAD, critic2 uses a uniform angular grid. Using QUOTIENT, the
  uniform grid is remapped by $$x^{l_1/l_2}$$ where $$l_1$$ and
  $$l_2$$ are the two negative eigenvalues at the bcp. This way, the
  points get accumulated around the bcp with lowest eigenvalue
  (highest if absolute value is taken). DYNAMICAL uses a linearized
  model of the interatomic surface and predicts the initial angles
  critic2 has to take in order to generate a uniform distribution of
  points a given distance away. This distance is calculated as 90% of
  the distance to the nearest ccp found in a coarse exploration of the
  IAS omega-limits. Unfortunately, this algorithm works only in cases
  where the bcp has significant but not too large ellipticty. Also,
  there is no gain in using this method in cases where the number of
  ccps is different than 4.

  By default, BRAINDEAD is used. H1 is experimental.

- RCP: print gradient paths starting at the neighbourhood of a ring
  CP. The situation is analogous to that of the bcps.

- CCP: print gradient paths starting at the vicinity of a cage
  CP. Again, the situation is symmetric to the ncp case.

- GRAPH: represent the complete graph in the unit cell. This means:

  * All the bond paths for which both ncps and the bcp lay inside
    the main unit cell.

  * All the ring paths for which both ccps and the rcp lay inside
    the main unit cell.

  The critical points on the boundary of the main cell are also
  represented.

  The `igraph.i` value represents the amount of information that is
  to be printed. It is the sum of these two values, each representing
  an element to plot:

  * 1 : print ring paths associated to the rcps.

  * 2 : print bond paths associated to the bcps.

Several parameters of the plot generated by FLUXPRINT can be
changed. The color to be applied to the paths resulting from FLUXPRINT
commands can be changed using the COLOR keyword (three integers from 0
to 255 for the red, green, and blue components). By default, NCP uses
the color corresponding to the originating atom, and a gold color is
used otherwise. The format of the output file can be controlled with
the following keywords: TEXT (plain text file), TESS or TESSEL
(tessel), OBJ (Wavefront obj), PLY (ply format), OFF (Geomview's off),
and CML (Chemical Markup Language). The CML format, which can be read
with the avogadro program, is used by default. The color
specifications are not passed to the CML format. To keep the number of
points manageable, critic2 writes one gradient path point every
certain distance along the path.

Finally, the SHELLS keyword applies only to graph and graphcp. It
represents the number of unit cell shells where the graph is going
to be plotted. Thus, 0 represents the main unit cell; 1, the main
unit cell and its 26 neighbours; and so on. By default, SHELLS
adopts the -1 value, which is equivalent to 0 for the graph
keyword and means that the partial graph generated in GRAPHCP is
not replicated by symmetry.

The NOSYM keyword instructs critic2 to write the complete list of CPs
in the unit cell, together with the identity matrix as the only
operation in the space group, to the output.

