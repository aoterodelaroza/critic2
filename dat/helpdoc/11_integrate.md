
## Overview

Critic2 provides several methods to integrate the attractor (Bader)
basins associated to the maxima of a field. In QTAIM theory, this
field is the electron density, the attractors are (usually) nuclei and
the basins are the atomic regions. The integrated properties are
atomic properties (e.g. atomic charges, volumes, moments, etc.). The
attractor basins are defined by a zero-flux condition of the electron
density: no gradient paths cross the boundary between attractor
regions. This makes the basins local to each attractor, but their
definition is a relatively complex algorithmic problem.

The simplest way of integrating an attractor basin is bisection. A
number of points distributed in a small sphere around the atom are
chosen, each of them determining a ray. On each ray, a process of
bisection is started. A point belongs to the basin if the gradient
path traced upwards ends up at the position of the attractor we are
considering. If the end-point is a different attractor, then the point
is not in the basin. By using bisection, it is possible to determine
the basin limit (called the interatomic surface, IAS). The bisection
algorithm is implemented in critic2, and can be accessed with the
INTEGRALS keyword. Bisection works best with analytical fields such as
molecular wavefunctions or WIEN2k, but can be used with any field.

The [qtree algorithm](https://doi.org/10.1002/jcc.21620)
is based on the recursive subdivision of the
irreducible Wigner-Seitz (IWS). In qtree, the smallest
symmetry-irreducible portion of space is considered and a tetrahedral
mesh of points is superimposed on it. The gradient path is traced from
all those points and the points are assigned to different atoms (the
points are "colored"). The integration is performed by quadrature over
the points belonging to a given basin. The qtree algorithm is accessed
through the QTREE keyword, and is suited for small crystals and
analytical fields (e.g. WIEN2k or elk).

Lastly, integration algorithms based on grid discretization are very
popular nowadays thanks to the widespread use of
pseudopotential/plane-waves DFT methods. Critic2 provides the
integration [method of Yu and Trinkle (YT)](https://doi.org/10.1063/1.3553716).
The algorithm is based on the assignment of integration
weights to each point in the numerical grid by evaluating the flow of
the gradient using the neighboring points. This algorithm is extremely
efficient and robust and is strongly recommended in the case of fields
on a grid. The associated keyword is YT.

An alternative to the YT method is the method proposed by Henkelman et
al. ([Comput. Mater. Sci. 36, 254-360 (2006)](https://doi.org/10.1016/j.commatsci.2005.04.010),
[J. Comput. Chem. 28, 899-908 (2007)](https://doi.org/10.1002/jcc.20575),
[J. Phys.: Condens. Matter 21, 084204 (2009)](https://doi.org/10.1088/0953-8984/21/8/084204)),
which is implemented in critic2 through the keyword BADER.

The field that determines the basins being calculated is always the
reference field (see
[REFERENCE](/critic2/manual/fields/#c2-reference)). In general, it is
necessary to define one or more properties to be integrated inside the
basins that use other scalar fields. For instance, to calculate
the electron population inside an ELF basin, the ELF would be the
reference field and the electron density would be an integrable
property.

## List of Properties Integrated in the Attractor Basins (INTEGRABLE) {#c2-integrable}

~~~
INTEGRABLE id.s {F|FVAL|GMOD|LAP|LAPVAL} [NAME name.s]
INTEGRABLE id.s {MULTIPOLE|MULTIPOLES} [lmax.i]
INTEGRABLE id.s DELOC [WANNIER] [PSINK] [NOU] [NOSIJCHK] [NOFACHK] [WANCUT wancut.r]
INTEGRABLE "expr.s"
INTEGRABLE DELOC_SIJCHK file-sij.s
INTEGRABLE DELOC_FACHK file-fa.s
INTEGRABLE CLEAR
INTEGRABLE ... [NAME name.s]
~~~
Critic2 uses an internal list of all properties that will be
integrated in the attraction basins. This list can be modified by the
user with the INTEGRABLE keyword. This keyword has a syntax similar to
the list of properties calculated at the critical points,
[POINTPROP](/critic2/manual/cpsearch/#c2-pointprop).

A single INTEGRABLE command assigns a new quantity to be integrated in
the atomic basins. The new integrable property is related to field
`id.s` (given as field number or identifier). This quantity can be the
field value itself (F), its valence component (if the field is
core-augmented, FVAL), the gradient norm (GMOD), the Laplacian (LAP),
or the valence-component of the Laplacian (LAPVAL). If no keyword is
given after `id.s`, F is used by default.

With the MULTIPOLES (or MULTIPOLE) keyword, the multipole moments of
the field are calculated up to an l equal to `lmax.i` (default:
5). This keyword only applies to the BADER and YT integration methods.
For the others, it is equivalent to the field value (same as F). The
units for all calculated multipoles are atomic units.

The keyword DELOC activates the calculation of the delocalization
indices (DIs) using field `id.s`. The calculation of DIs requires that
the field is a Quantum ESPRESSO [pwc
file](/critic2/manual/fields/#c2-qepwc), which contains the Bloch
states. There are two modes of operation with DELOC:

- Calculate the DIs using maximally localized Wannier functions
  (MLWF). This requires wannier90 checkpoint files to be provided with
  the pwc file, and is the default behavior in this case. More details
  are [given below](/critic2/manual/integrate/#c2-intwandi). The use
  of Wannier functions can be forced using the WANNIER keyword in
  INTEGRABLE DELOC.

- Calculate the DIs using Bloch states. This is much less efficient
  than using Wannier functions but it is simpler because it does not
  require an additional wannier90 calculation and allows calculating
  DIs in metals. More details are [given
  below](/critic2/manual/integrate/#c2-intblodi). The use of Bloch
  states can be forced using the PSINK keyword in INTEGRABLE DELOC.

Also, see the
[example](/critic2/examples/example_11_10_deloc-indices/) for the
calculation of DIs in solids.

Both ways of calculating the DIs generate by default two checkpoint
files containing the (complex) overlap matrices (sij) or the matrix of
atomic integrals of the exchange-correlation density (fa). The
delocalization indices can be calculated from these files bypassing
having to load the pwc and wannier90 checkpoint files with the
DELOC_SIJCHK and DELOC_FACHK keywords. This is useful because the pwc
files can be quite large, so it may be interesting to delete them
after the sij or fa matrices are computed. Both keywords accept a
string pointing to the corresponding file. The DI calculation from
these matrices will typically take only a few seconds.

In addition, it is possible to define an integrable property using an
expression involving more than one field (`expr.s`). For instance, if
the spin-up density is in field 1 and the spin-down density is in
field 2, the atomic moments can be obtained using:
~~~
LOAD AS "$1+$2"
REFERENCE 3
INTEGRABLE "$1-$2"
~~~
Note that the quotation marks are required.

The additional keyword NAME can be used with any of the options above
to change the name of the integrable property, for easy identification
in the output.

The keyword CLEAR resets the list to its initial state (volume,
electron population, and Laplacian in crystals; electron population
and Laplacian in molecules). Using the INTEGRABLE keyword will print a
report on the list of integrable properties.

The default integrable properties are:

* Volume (`1`), in crystals only.

* Pop (`fval`): the value of the reference field is integrated. If the
  reference field is the density, then this is the number of electrons
  in the basin. If core augmentation is active for this field, only
  the valence contribution is integrated.

* Lap (`lap(fval)`): the Laplacian of the reference field. The
  integrated Laplacian has been traditionally used as a check of the
  quality of the integration because the exact integral is zero
  regardless of the basin (because of the divergence
  theorem). However, it is difficult to obtain a zero in the Laplacian
  integral in critic2 in some cases because of numerical inaccuracies:

  - In fields based on a grid, the numerical interpolation gives a
    noisy Laplacian.

  - In FPLAPW fields (WIEN2k and elk), the discontinuity at the muffin
    surface introduce a non-zero contribution to the integral.

  If f is a core-augmented field, only the valence Laplacian is
  integrated.

### Integrating Delocalization Indices in a Solid With Maximally Localized Wannier Functions {#c2-intwandi}

The keyword DELOC activates the calculation of localization and
delocalization indices (DIs) in a crystal using the procedure
described in the
[literature](http://dx.doi.org/10.1021/acs.jctc.8b00549). DIs can be
calculated only if the loaded field contains information about
individual Kohn-Sham states and the orbital rotation that leads to the
maximally localized Wannier functions (MLWF). This is done by using a
field loaded from a `.pwc` file (generated by the `pw2critic.x`
utility in Quantum ESPRESSO) together with a checkpoint (chk) file
from wannier90. The former contains the
[electronic wavefunctions](/critic2/manual/fields/#c2-qepwc) and
the latter the orbital rotation.  For maximum consistency, the pwc
file can also be used to provide the structural information for the
run via the [CRYSTAL keyword](/critic2/manual/crystal/#c2-crypwc).

In addition to these data, the calculation of DIs has a few
requirements: the grid must be consistent with that of the reference
field, and the DIs can be calculated using YT or BADER only.

A typical delocalization index calculation comprises the following
steps:

- Run a PAW calculation, then obtain an all-electron density using
  `pp.x` with `plot_num=21`. This creates a cube file (rhoae.cube)
  that gives the Bader basins for the calculation (the pseudo-valence
  density is not valid for this purpose).

- Run an SCF calculation with norm-conserving pseudopotentials and the
  same ecutrho as the calculation in the first step, so the two grids
  have the same size.

- Use the `open_grid.x` utility in Quantum ESPRESSO to unpack the
  symmetry of the k-point grid, in preparation for the wannier90 run
  (this would normally be accomplished by a non-SCF calculation, but
  with `open_grid.x` it is easier and much faster).

- Use `pw2critic.x` on the output of `open_grid.x` to generate the pwc
  file.

- Run wannier90 on the result of `open_grid.x` to generate the chk
  file.

- Load the all-electron density and the pwc and chk files as two
  fields in critic2. Set the former as the reference density and the
  latter as INTEGRABLE DELOC. If the system is spin-polarized, two
  checkpoint files will be necessary, one for each spin component.

- Run YT or BADER. YT is usually more accurate but takes longer than
  BADER.

Additional options for the DI calculation follow. The NOU option
disables the use of the U rotation matrices to calculate the
MLWFs. This makes critic2 calculate the DI using Wannier functions
obtained by using a straight Wannier transformation from the Bloch
states. This is naturally much slower than the maximally localized
version, since overlaps cannot be discarded, and should be used only
if wannier90 failed to converge for the particular case under study.

By default, two checkpoint files are generated during a DI calculation
run. These files have the same name as the pwc file but with `-sij`
and `-fa` suffixes. The former checkpoint file contains the atomic
overlap matrices, and the latter, the exchange-correlation density
($$F_{AB}$$) integrals required for the DI calculation. The presence
of any of these two files makes critic2 read the information from the
files and bypass the corresponding calculations, which are quite time
consuming in general. The keywords NOSIJCHK and NOFACHK deactivate
reading and writing these checkpoint files. You can calculate the DIs
from these files directly without the pwc file using the DELOC_SIJCHK
and DELOC_FACHK options to INTEGRABLE.

By default, the overlap between two MLWFs whose centers are a certain
distance away are discarded. The WANCUT keyword controls this
distance: overlaps are discarded if the centers are `wancut.r` times
the sum of their spreads away. By default, `wancut.r = 4.0`. A very
large `wancut.r` will prevent critic2 from discarding any overlaps. The
appropriateness of the chosen WANCUT can be checked a posteriori by
comparing the integrated electron population obtained by sum of the
localization and delocalization indices to the value obtained from a
straight integration of the electron density.

See the [example](/critic2/examples/example_11_10_deloc-indices/) for
some sample input files and more details.

### Integrating Delocalization Indices in a Solid With Bloch States {#c2-intblodi}

The DIs can also be calculated without doing the transformation to
Wannier functions, using the Bloch states provided by the
calculation. This is typically far less efficient than using Wannier
functions but has two advantages:

- It is simpler because it does not require the additional wannier90
  step.

- Using Bloch states, DIs can be calculated in metals (the Wannier
  transformation is ill-defined if some bands are partially
  occupied).

The calculation of DIs using Bloch states requires a `.pwc` file only,
generated by the `pw2critic.x` utility in Quantum ESPRESSO. The pwc
file contains [electronic
wavefunctions](/critic2/manual/fields/#c2-qepwc) and the crystal
structure. For maximum consistency, the pwc file should be used to
provide the structural information for the run via the [CRYSTAL
keyword](/critic2/manual/crystal/#c2-crypwc).

In addition to these data, the calculation of DIs with Bloch states
has a few requirements: the grid must be consistent with that of the
reference field, and the DIs can be calculated using YT or BADER only.
A typical delocalization index with Bloch states calculation comprises
the following steps:

- Run a PAW calculation, then obtain an all-electron density using
  `pp.x` with `plot_num=21`. This creates a cube file (rhoae.cube)
  that gives the Bader basins for the calculation (the pseudo-valence
  density is not valid for this purpose).

- Run an SCF calculation with norm-conserving pseudopotentials and the
  same ecutrho as the calculation in the first step, so the two grids
  have the same size. This calculation must be run without k-point
  symmetry, i.e., with `nosym=.true`, and `noinv=.true.`.

- Use `pw2critic.x` to generate the pwc file.

- Load the all-electron density and the pwc files as two fields in
  critic2. Set the former as the reference density and the latter as
  INTEGRABLE DELOC.

- Run YT or BADER. YT is usually more accurate but takes longer than
  BADER.

By default, two checkpoint files are generated during a DI calculation
run. These files have the same name as the pwc file but with `-sij`
and `-fa` suffixes. The former checkpoint file contains the atomic
overlap matrices, and the latter, the exchange-correlation density
($$F_{AB}$$) integrals required for the DI calculation. The presence
of any of these two files makes critic2 read the information from the
files and bypass the corresponding calculations, which are quite time
consuming in general. The keywords NOSIJCHK and NOFACHK deactivate
reading and writing these checkpoint files. You can calculate the DIs
from these files directly without the pwc file using the DELOC_SIJCHK
and DELOC_FACHK options to INTEGRABLE.

The keywords NOU and WANCUT have no effect on the calculation of DIs
using Bloch states.

See the [example](/critic2/examples/example_11_10_deloc-indices/) for
some sample input files and more details.

## Bisection (INTEGRALS and SPHEREINTEGRALS) {#c2-integrals}

~~~
INTEGRALS {GAULEG ntheta.i nphi.i|LEBEDEV nleb.i}
          [CP ncp.i] [RWINT] [VERBOSE]
~~~
The BISECTION keyword integrates the attractor basins using
bisection. If the Gauss-Legendre quadrature is used (GAULEG keyword),
`ntheta.i` and `nphi.i`  are the number of $$\theta$$ (polar
angle) and $$\phi$$ (azimuthal angle) points.

In the case of a Lebedev-Laikov quadrature, selected via the LEBEDEV
keyword, only the total number of points in the spherical quadrature
is needed. The actual value of `nleb.i` is the smallest number larger
than the one given by the user that is included in the list: 6, 14,
26, 38, 50, 74, 86, 110, 146, 170, 194, 230, 266, 302, 350, 434, 590,
770, 974, 1202, 1454, 1730, 2030, 2354, 2702, 3074, 3470, 3890, 4334,
4802, 5294, 5810.

By using the CP keyword, a single non-equivalent CP (`ncp.i`) is
integrated. Otherwise, all the CPs of the correct type (found using
AUTO) are integrated. If RWINT is present, read (if they exist) and
write the `.int` files containing the interatomic surface limit for
the rays associated to the chosen quadrature method.

Defaults: `ntheta.i = nphi.i = 50`, `nleb.i = 4802`.

~~~
SPHEREINTEGRALS {GAULEG ntheta.i nphi.i|
                 LEBEDEV nleb.i} [CP ncp.i] [NR npts.i]
                [R0 r0.r] [REND rend.r]
~~~
The SPHEREINTEGRALS keyword integrates the volume, field and Laplacian
of the reference field in successive spheres centered around each of
the attractor CPs. The meaning of the GAULEG and LEBEDEV keywords is
the same as in INTEGRALS.

A total number of `npts.i` spheres are integrated per nucleus. The
grid is logarithmic, so that the region near the nucleus has a higher
population of points. The grid starts at the radius `r0.r` and ends at
`rend.r` (bohr in crystals, angstrom in molecules). If `rend.r` < 0 then
the final radius is taken as half the nearest neighbor distance for
each atom times `abs(rend.r)`.

Default: `npts.i = 100`. In GAULEG, `ntheta.i = 20` and `nphi.i = 20`. In
LEBEDEV, `nleb.i = 770`. `r0.r = 1d-3` bohr. `rend.r = rnn/2` for each
CP. `id.i = 0` (all attractors).

## Qtree (QTREE) {#c2-qtree}

### General Syntax

The [QTREE integration method](https://doi.org/10.1002/jcc.21620)
integrates QTAIM atomic properties by
discretization of the smallest part of the crystal that reproduces the
whole system by symmetry. QTREE is specific for periodic crystals,
and its use is recommended with fields not given on a grid. For fields
on a grid, either
[BADER](/critic2/manual/integrate/#c2-bader) or
[YT](/critic2/manual/integrate/#c2-yt) are better alternatives. QTREE
is based on a hierarchical subdivision of the irreducible part of
the WS cell, employing a tetrahedral grid. The integration region is
selected so as to maximize the use of symmetry, and partitioned into
tetrahedra. The syntax of QTREE consists of the maximum subdivision
level (`maxlevel.i`), perhaps followed by the pre-splitting level
(`plevel.i`), and then a series of optional keywords that control the
behavior and the options for the integration.
~~~
QTREE [maxlevel.i [plevel.i]] [MINL minl.i] [GRADIENT_MODE gmode.i]
      [QTREE_ODE_MODE omode.i] [STEPSIZE step.r] [ODE_ABSERR abserr.r]
      [INTEG_MODE level.i imode.i] [INTEG_SCHEME ischeme.i] [KEASTNUM k.i]
      [PLOT_MODE plmode.i] [PROP_MODE prmode.i] [MPSTEP inistep.i]
      [QTREEFAC f.r] [CUB_ABS abs.r] [CUB_REL rel.r] [CUB_MPTS mpts.i]
      [SPHFACTOR {ncp.i fac.r|at.s fac.r}] [SPHINTFACTOR atom.i fac.r]
      [DOCONTACTS] [WS_ORIGIN x.r y.r z.r] [WS_SCALE scale.r]
      [NOKILLEXT] [AUTOSPH {1|2}] [CHECKBETA] [NOPLOTSTICKS]
      [COLOR_ALLOCATE {0|1}] [SETSPH_LVL lvl.i] [VCUTOFF vcutoff.r]
~~~
In QTREE, the tetrahedra that comprise the IWS enter a recursive
subdivision process in which each of the tetrahedra is divided in 8 at
each level, up to a level given by the user. This subdivision level is
controlled by the `maxlevel.i` argument given after the QTREE keyword
(default, 6). Every tetrahedron vertex is assigned to a non-equivalent
atom in the unit cell by tracing a gradient path and identifying its
endpoint. Once all the points in the tetrahedral mesh are assigned,
the tetrahedra are integrated and the properties assigned to the
corresponding atoms. The space near the atoms is integrated using a
beta-sphere, which improves the accuracy of the integration.

In the simplest (and most common) approach, qtree can be executed using:
~~~
QTREE [maxlevel.i [plevel.i]]
~~~
where `maxlevel.i` is the level of subdivision. The optional `plevel.i`
value corresponds to the pre-splitting level of the tetrahedra. The
initial tetrahedra list is split into smaller tetrahedra `plevel.i`
times. This can be useful in cases where a very high accuracy (and
therefore a very high level) is required, but there is not enough
memory available to advance to higher `maxlevel.i`. However, using
a relatively high `plevel.i` incurs an overhead, because the
atom assigning procedure is not as efficient when smaller tetrahedra
are used.

### Steps of the QTREE Algorithm

1. The WS cell is constructed and split into tetrahedra, all of which
   have in common, at least, the origin of the WS cell. Then, the site
   symmetry of the origin is calculated and the tetrahedra that are
   unique under the operations of this group are found. This is what
   we call the irreducible Wigner-Seitz cell (IWS). Note, however,
   that it is only "irreducible" in the local site symmetry of the
   origin, not in the full set of space group operations. The IWS is
   the region that is integrated in later steps of QTREE. We
   will refer to a IWS tetrahedra as an IWST.

   It is possible, through the WS\_ORIGIN keyword, to shift the origin
   of the WS cell away from the (0 0 0) position. Because the symmetry
   of the WS cell is determined by the site symmetry of the origin,
   the number and shape of IWST change depending on the origin
   chosen. A general position (with no symmetry) will make the IWS
   exactly equal to the WS.

   Also, for large systems, the user can choose to shrink the size of
   the original WS cell in order to integrate a smaller region, using
   the WS\_SCALE keyword (most likely in combination with WS\_ORIGIN,
   to move the integration region around). If a value is given to
   WS\_SCALE (for instance, `rws.r`), then all the vectors connecting
   the origin of the WS cell with the vertex are shrunk by a factor
   `rws.r`, and the volume of the integration region is decreased by
   a factor equal to the cube of `rws.r`. The IWS is calculated using
   the smaller WS cell and integrated in the same way. Note that this
   integration region is non-periodic: it does not fill the volume of
   the solid and it does not integrate to the total number of
   electrons per cell.

2. Non-overlapping spheres are chosen centered on each of the atoms in
   the cell, the so-called beta-spheres. Atoms equivalent by symmetry
   share the same beta-sphere radius ($$\beta_i$$ for atom $$i$$). The
   beta-sphere has two roles in QTREE:

   -  The atomic properties are integrated inside the beta-spheres
      using a 2D cubature. The cubature can be a product of two 1D
      Gauss-Legendre quadratures or a Lebedev quadrature of the
      sphere. Both methods, and the number of nodes can be selected
      using the INT\_SPHEREQUAD\_* keywords shown below. The radial
      quadrature can be any of the available in critic2, and is
      controlled by the INT\_RADQUAD\_* options. The default values,
      however, are usually fine, integrating the beta-spheres in a
      matter of seconds with a precision that is orders of magnitude
      better than the overall QTREE performance.

      This beta-sphere integration removes the error of the
      finite-elements integration of a region where the integrated
      scalar fields present the steeper variations in value. By
      removing the high-error regions from the grid integration, the
      accuracy of QTREE is enhanced. In particular, this increase in
      precision outweighs the loss of precision caused by creating an
      additional interface between the grid and the sphere.

   - The space inside the beta-sphere of an atom is assumed to be
     inside the basin of that atom. The terminus of any gradient path
     that reaches the interior of the beta-sphere $$i$$ is assumed to
     be the atom $$i$$. It is known that most of the steps in the
     integration of the gradient of the electron density are spent in
     the close vicinity of the terminus. Therefore, this modification
     saves precious function evaluations.

   The default beta-sphere radius is set to $$0.80$$ times half the
   nearest-neighbor distance. Both i) and ii) above assume that the
   beta-sphere is completely contained inside the basin of the
   atom. This may turn out not to be true for the default beta-sphere
   radius (specially for cations in ionic systems). In these cases, the
   keyword SPHFACTOR is used. An example of the use of this keyword
   is:
~~~
SPHFACTOR 1 0.70   # Make b_1 = 0.70 * rNN2(1) (atom type 1)
~~~
   If SPHFACTOR < 0, use the [method by Rodriguez et al.](https://doi.org/10.1002/jcc.21134)
   to determine the beta-sphere radii: a collection of
   points around the atom are selected and the angle between the
   gradient and the radial direction is determined. If all the angles
   are less than 45 degrees, the sphere is accepted. In solids, this
   strategy usually yields spheres that are too large.

   In the case that any atomic SPHFACTOR is zero (the default value
   for all atoms), then a pre-computation at a lower level is done to
   ensure that all beta-spheres lie within the desired basins. There
   are two methods for this, and the method can be chosen using
   the AUTOSPH keyword.

   Method number one involves using a reduced version of QTREE. The
   pre-computation usually takes no longer than some minutes (and
   usually only few seconds) and the spheres are guaranteed to be
   inside the basins. The keyword SETSPH\_LVL controls the level of the
   pre-computation, that must not be higher than 7. The default value
   is 6 or `maxlevel.i`, whichever is smaller.

   The second method (default) traces gradient path on a coarse sphere
   around each nucleus, and reduces the sphere until all of the points
   are inside the basin. NOCHECKBETA is used in this case.

   An additional factor the user can define is the SPHINTFACTOR. It is
   possible to consider the sphere where GPs terminate different from
   the one that is integrated. If SPHINTFACTOR is defined, for
   example as:
~~~
SPHINTFACTOR 1 0.75
~~~
   then the sphere associated to atom 1 where the integration is done
   has a radius which is 0.75 times that of the sphere where GP
   terminate.

   The CHECKBETA and NOCHECKBETA keywords activate and deactivate the
   check that ensures that the beta-spheres is completely contained
   inside the basin.

   If a beta-sphere is not strictly contained in the basin, QTREE
   detects it and stops immediately (specifically, QTREE checks that
   every tetrahedron that is partly contained in a beta-sphere has
   vertex termini that are all assigned to the same atom as the
   atom that owns the beta-sphere).

   For a new system, it is always a good idea to start with a
   low-level QTREE (say, level 4) to check if the default beta-spheres
   are adequate. If one of the beta-spheres is too large, the error
   message looks like:
~~~
An undecided tetrahedron is overlapping with a beta-sphere.
Make beta-spheres smaller for this system.
terms:            1          -1           2           1
~~~
   which indicates that there is a tetrahedron that is partly
   contained in the sphere of the first atom (terminus -1) and that
   has a vertex corresponding to atom 2. Modifying the sphfactor
   solves this problem:
~~~
SPHFACTOR 1 0.70
~~~
   At lower levels, QTREE is reasonably fast, so a trial-and-error
   selection of beta-spheres is acceptable. Note that the
   beta-spheres used in QTREE have no relation to the ones reported
   after an AUTO calculation.

3. If the cell is periodic (which means that WS\_SCALE was not set),
   the contacts between the faces of the IWST are found. These
   contacts are used in a later step to copy the termini information
   between tetrahedron faces.

   The determination of the tetrahedra contacts in a periodic
   integration region is deactivated if the NOCONTACTS keyword is
   used. The DOCONTACTS keyword does the opposite thing: it activate
   the calculation of contacts. By default, the contacts are not
   calculated.

4. A grid is built for each of the IWST. For each of the grid points,
   the termini of the gradient paths starting at them is
   calculated. The positions of the grid points are determined by
   subdividing the IWST `maxlevel.i` times. In each subdivision step,
   a parent tetrahedron is divided in 8 smaller tetrahedron (all with
   the same volume, $$V/8$$) by splitting each edge of the parent
   tetrahedron in two. There are two possible ways of doing this, the
   choice being irrelevant to the performance of QTREE.

   For a given IWST, the size of the grid is given by
$$
\begin{equation}
S = n*(n+1)*(n+2)/6
\end{equation}
$$
   where $$n = 2^l + 1$$ and $$ll$$ is the `maxlevel.i`, the
   approximate scaling being as $$8^l$$. The termini information on
   the grid is saved to the array `trm` of type `integer*1`, with size
   $$n_t * S$$, where $$n_t$$ is the number of IWST (in fact, the
   integer type is that which is the result of
   `selected_int_kind(2)`).

   The subdivision level is the main input parameter for QTREE,
   controlling the accuracy (and cost!) of the integration. For
   small-medium sized systems $$4-5$$ are low cost integrations
   (seconds), $$6-7$$ are medium cost (minutes) and $$8-9$$ are the
   slowest and most accurate (hours). The level is input in the call
   to the QTREE integration:
~~~
QTREE 6   # maxlevel.i is 6
~~~
   By default, the integration level is 5.

   In addition to `trm`, more work space can be allocated if the
   integration is restricted to the volume and charge or to the
   volume, charge and Laplacian. The number and type of properties to
   integrate is controlled by the PROP\_MODE keyword. The following
   values are allowed:

   - 0: Only volume is integrated. This amounts to canceling the
     finite elements integration of tetrahedra and is equivalent to
     INTEG\_MODE 0 (see below).

   - 1: Only charge and volume. If the integration uses the value of
     the density at the grid points (INTEG\_MODE 11, see below) In
     addition to `trm`, another `real*8` array, `fgr`, is allocated
     (strictly it is `selected_real_kind(14)`, not `real*8`). In fgr,
     the value of the density at the grid points is stored.

   - 2: Charge, volume and Laplacian. In a similar way to 1, if the
     information on the grid points is used during the integration
     (INTEG\_MODE 11), an additional `real*8` array, `lapgr`, is
     allocated. It contains the value of the Laplacian of the electron
     density at the grid points.

   - 3: All the integrable properties calculated by the module. The
     number of properties varies with the interface being used. No
     `fgr` or `lapgr` are allocated, as the grid points need to be
     recomputed during integration.

   The default value for PROP\_MODE is 2.

   The termini of the grid points contained in a beta-sphere is marked
   previous to the beginning of the subdivision.

5. Each tetrahedron is subdivided recursively up to a level
   `maxlevel.i`, and integrated at the same time. The IWST integration
   is relatively independent of one another, so for the moment we will
   focus on just one IWST, which we will call the base tetrahedron.

   The result of the integration of the different IWST can be used for
   other IWST:

   - When the integration ends, the termini of the four faces of a
     base tetrahedron are copied to its neighbors' `trm`, according to
     the contacts determined previously.

   - Depending of the method chosen (see GRADIENT\_MODE below), the
     gradient path integration may be aware of the neighboring grid
     points, that may belong to other IWST. In particular, the
     gradient mode number 3 integrates a gradient path following grid
     points. When the endpoint is reached, all the grid points that
     have been traversed by the path are assigned the same common
     terminus. Therefore, there is the possibility that gradient paths
     starting inside a given base tetrahedron write the `trm` of other
     IWST.

6. A tetrahedra stack is built and initialized with one element: the
   base tetrahedron. An iterator works on the stack, performing at
   each step the following tasks:

   - Pop a tetrahedron from the stack.

   - The termini of the vertex of the tetrahedron are calculated, if
     they are not already known.

   - If all the termini of the tetrahedron correspond to the same
     atom, the tetrahedron is "painted". This means that all the grid
     points that are in the interior or border of the tetrahedron are
     assigned the same color as its vertex, thereby saving the tracing
     of the gradient paths.

     "Painting" can be dangerous whenever a (curved) IAS crosses the
     face of the tetrahedron being painted. To this end, a minimum
     level is defined, using the keyword QTREE\_MINL. If the
     subdivision level of the tetrahedron is lower or equal than
     `minl.i`, the tetrahedron is not painted. (Note: the base
     tetrahedron corresponds to level 0).

     Furthermore, if all the termini correspond to the same atom and
     are located outside of the beta-sphere region, the tetrahedron
     is integrated and does not enter another subdivision
     process. Once more, this only happens to tetrahedra with a
     level of subdivision strictly greater than `minl.i`. The "inner
     integration method" is the quadrature method used to integrate
     these tetrahedra.

     If all the termini are located inside the beta-sphere region, the
     tetrahedron does not subdivide, but the properties are not
     integrated, because this region corresponds to the sphere
     integration addressed in point 2.  If the tetrahedron is across
     the border of a beta-sphere, it is divided further.

   - If the tetrahedron is at subdivision level equal to `maxlevel.i`,
     then it does not subdivide, it is integrated and the properties
     are assigned to the atoms. There are several possibilities
     depending of the nature of its vertex' termini:

     + If all the termini correspond to the same atom, and the
       tetrahedron is completely inside or outside of this atom's
       sphere, it corresponds to the case in 6.3.

     + If it is completely inside an atom basin, but on the border of
       a beta-sphere, the part of the tetrahedron that is outside of
       the sphere is integrated and assigned to the atom. Another
       integration method is required for this, different from the
       "inner integration". We will refer to this method as "border,
       same-color integration".

     + If the tetrahedron has termini corresponding to different
       atoms, its properties are integrated and split into
       contribution to atoms, according to the number of termini each
       atom has. These tetrahedra are located on the IAS, and require
       a third class of integration, "border, diff-color integration".

   - A tetrahedron that has not been integrated continues the
     subdivision. In this step, 8 new tetrahedra are pushed into the
     stack. The edges of the parent tetrahedron are split in two. By
     construction, the newly generated points also correspond to grid
     points.

     The subdivision scheme is:

     + 1, 1-2, 1-3, 1-4
     + 2, 1-2, 2-3, 2-4
     + 3, 1-3, 2-3, 3-4
     + 4, 1-4, 2-4, 3-4
     + 2-3, 1-2, 1-3, 1-4
     + 1-4, 1-2, 2-3, 2-4
     + 1-4, 1-3, 2-3, 3-4
     + 2-3, 1-4, 2-4, 3-4

     where "a" represents a vertex of the parent tetrahedron and "a-b"
     the midpoint of both vertex. Each of the 8 child tetrahedra
     enclose the same volume, equal to $$V/8^l$$, where $$l$$ is the
     subdivision level and $$V$$ is the volume of the base
     tetrahedron.

   - When the stack is empty, the work on the base tetrahedron is
     finished.

7. "Inner integration". The inner integration is a quadrature that is
   applied to tetrahedra that are completely contained in the
   non-beta-sphere region of a basin. It can apply to a tetrahedron of
   any level, as long as this level is greater than `minl.i`. The
   integrated properties are assigned to a single atom.

   In the current implementation of QTREE, several integration methods
   are possible, and are controlled by the INTEG\_MODE keyword. The
   possible values of INTEG\_MODE are:

   - 11 : use the information of the density, Laplacian and properties
     at the vertex of the tetrahedron to integrate. The integral is
     approximated by a quadrature of four terms, each corresponding to
     a volume that is 1/4 of the volume of the tetrahedron and
     multiplied by the value of the properties at the vertex. This
     integration method is useful if only charge or charge and
     Laplacian are being integrated, because the information gained
     during the gradient path tracing, and saved in the `fgr` and
     `lapgr` arrays, is used. It is not very accurate for large
     tetrahedra.

   - 12 : use the CUBPACK routines. CUBPACK provides an adaptive
     tetrahedron integration method based on recursive subdivision and
     an integration rule with 43 nodes (degree 8), that is equivalent
     to the DCUTET library by Bernsten et al. The integration rule is
     fully symmetric under the Th group operations. The error
     estimation is compared to the error requested by the user, that
     is controlled using the CUB\_ABS (absolute error), CUB\_REL
     (relative error) and CUB\_MPTS (maximum number of function
     evaluations) keywords. If CUB\_MPTS is exceeded, an error message
     is output, but the QTREE integration continues.

     Note that, no matter how low the error requirements are, the
     CUBPACK integration spends, at least, 43 function evaluations per
     tetrahedron, so it is quite expensive if compared to other
     integration modes. This should be used for large tetrahedra (see
     below) or for really accurate calculations.

   - 1...10 : use a non-adaptive rule from the KEAST library (Keast et
     al., 1986), the number corresponding to:

     + 1  -- order =  1,  degree =  0
     + 2  -- order =  4,  degree =  1
     + 3  -- order =  5,  degree =  2
     + 4  -- order = 10,  degree =  3
     + 5  -- order = 11,  degree =  4
     + 6  -- order = 14,  degree =  4
     + 7  -- order = 15,  degree =  5
     + 8  -- order = 24,  degree =  6
     + 9  -- order = 31,  degree =  7
     + 10 -- order = 45,  degree =  8

     In particular, the first KEAST rule uses the barycenter of the
     tetrahedron.

   The syntax of the INTEG\_MODE keyword is:
~~~
INTEG_MODE lvl.i mode.i
~~~
   where `mode.i` is one of the modes above and `lvl.i` is the level
   to which it applies. This means that, if a tetrahedron of a given
   level is integrated, the value of INTEG\_MODE(level) is checked to
   decide on the method.

   Another INTEG\_MODE value is possible:

   - -1 : do not integrate and force the tetrahedron into the
     subdivision process. This value of INTEG\_MODE can be combined
     with a positive value at higher levels, amounting to a recursive
     integration in the style of CUBPACK. Of course, -1 is not an
     acceptable value of INTEG\_MODE for the last level, `maxlevel.i`.

   As setting these INTEG\_MODE by hand could be confusing, QTREE
   provides sets of INTEG\_MODE values,
   which we will call "integration schemes". An integration scheme
   is a full set of INTEG\_MODEs for all levels. Integration schemes
   are selected with the INTEG_SCHEME keyword, that can have the
   following values:

   - 0: do not integrate, only calculate volume and plot (see
     below). This is equivalent to setting PROP\_MODE to 0.

   - 1: subdivide each tetrahedron up to the highest level and then
     integrate using the vertex information. This is most useful if
     PROP\_MODE is 1 (only charge and volume) or 2 (charge, volume and
     Laplacian) because the information of the gradient path tracing
     (`fgr` and `lapgr` arrays) are used:
~~~
INTEG_MODE = -1 -1 ... -1 11
!            ^^           ^^
!        QTREE_MINL    `maxlevel.i`
~~~
   - 2: subdivide each tetrahedron up to the highest level and then
     integrate using the barycenter.
~~~
INTEG_MODE = -1 -1 ... -1 1.
~~~
   - 3: barycentric integration at all levels of subdivision. Less
     accurate but faster than 2.
~~~
INTEG_MODE = 1 1 ... 1 1.
~~~
   - 4: one of the Keast rules (given by the KEASTNUM keyword) is
     used at all levels. If KEASTNUM is n,
~~~
INTEG_MODE = n n ... n n.
~~~
   - 5: CUBPACK, at all levels. Reserve this one for special
     occasions.
~~~
INTEG_MODE = 12 12 ... 12 12.
~~~
   - 6: this scheme and the next are (poor) attempts at trying an
     adaptive integration scheme. They are not more reliable or
     efficient than, for instance, scheme 2. Integration scheme 6
     calculates levels 4, 5, and 6 using CUBPACK, and the rest with
     subdivision up to the highest level and vertex-based
     integration.
~~~
INTEG_MODE = 12 12 12 -1 ... -1 11.
~~~
   - 7: same as 6 but the final integration uses only the barycenter.
~~~
INTEG_MODE = 12 12 12 -1 ... -1 1.
~~~
   - -1: let the user enter the INTEG_MODEs by hand.

   The default integration scheme is 2, suitable for low and
   medium-accuracy calculations.

8. "Border, same-color integration". This integration method applies
   to tetrahedra that have reached the maximum subdivision level and
   sit on the interface between a beta-sphere and the atomic
   basin. Some of the vertex are inside the sphere and some are
   outside. The out-of-sphere part is integrated and added to the
   atomic properties, while the in-the-sphere part is ignored because
   it has already been integrated.

   The integration works by assuming that the sphere radius is much
   larger than the tetrahedron characteristic lengths and, therefore,
   that the sphere surface can be considered a plane that intersects
   the tetrahedron. The intersection points of the sphere with the
   tetrahedron edges are easily calculated and, for simplicity, we
   refer to them as the "middle" of the edges. There are three cases:

   - One vertex is outside, three inside. The tetrahedron formed by
     the vertex outside and the three middle points of the edges that
     stem from it form a tetrahedron by itself, which is integrated
     and added to the atom properties.

   - Three vertex are outside, one inside. The difference between the
     whole tetrahedron integration and the small tetrahedron inside
     the sphere is added to the atom properties. The small tetrahedron
     is formed by the vertex that is inside the sphere and the three
     edges connected to it.

   - Two vertex are inside, two outside. The region outside of the
     sphere is a "triangular prism", that is split in three tetrahedra
     and integrated.

   Note that the INTEG\_MODE of the maximum subdivision level
   (`maxlevel.i`) applies to all the sub-integrations of the border,
   same-color integration.

9. "Border, diff-color integration". As in the case of "border,
   same-color integration", this method only applies to tetrahedra
   which are at their maximum subdivision level. In this case, the
   termini of the vertex corresponding to, at least, two different
   atoms.

   In the current implementation of QTREE, the tetrahedron is
   integrated as a whole. Then, the properties are equally
   assigned to each of the termini atoms. For instance, if the termini
   are (1 1 1 3), the properties of the tetrahedron are integrated,
   then 3/4 of them assigned to atom 1 and 1/4 to atom 3.

10. Gradient path tracing. The gradient path start always at grid
    points, and are traced using one of three methods, controlled by
    the GRADIENT\_MODE keyword, that can assume the following values:

    - 1: "full gradient". This method is plain ODE integration.  It
      is carried out ignoring the grid information. The gradient path
      is terminated whenever it enters a beta-sphere region.

    - 2: "color gradient". At each point of the gradient path, the
      neighboring grid points are checked. If all of them correspond
      to the same atom, then the terminus of the gradient path is
      assigned to that atom.

      In a tetrahedral mesh, the meaning of "neighboring grid points"
      is not as clear as in a cubic mesh. For a given point x, the
      neighbors are calculated by first converting x to convex
      coordinates, that range from 0 to $$2^l$$, restricted to
      $$x_1 + x_2 + x_3 \leq 2^l$$, where $$l$$ is the subdivision
      levels. The neighboring points are $$(x_1\pm 1, x_2\pm 1, x_3\pm 1)$$
      If any of these neighbors are not valid points in the
      tetrahedron, they are discarded and not checked. This is the
      default, except in the grid module.

    - 3: "qtree gradient". This method behaves much like the "full
      gradient", but whenever the gradient path steps near a grid
      point, it is projected to it. When a projection occurs, the grid
      point is pushed onto a stack. At the end of the gradient path,
      when the terminus is known, all the grid points in the stack are
      popped and assigned the terminus.

      The projection regions are spheres located around each grid
      point, whose radius is controlled by the QTREEFAC keyword. The
      radius of these spheres is $$m_l/2^{l_{\rm max}}/{\rm
      qtreefac}$$, where $$m_l$$ is the smallest edge length of all
      IWST and $$l_{\rm max}$$ is the maximum subdivision level
      (`maxlevel.i`). Note that QTREEFAC equal to 1 is the maximum
      value allowed, and corresponds to touching spheres along at
      least one tetrahedron edge. By default, QTREEFAC is 2. Lower
      levels of QTREEFAC tend to give errors when assigning the grid
      points that lie on the IAS of two atoms (although only
      there). With higher levels, the time saving diminishes, and
      "qtree gradient" is equivalent to "full gradient".

      Additionally, the projection can be started only after a certain
      number of initial steps. The MPSTEP keyword controls this value.
      The default MPSTEP is 0.

    - -1, -2, -3: these correspond to the same as their positive
      values, but each gradient path terminus is compared to their
      "full gradient" version, using the best available ODE
      integration method. Information about the results of the
      comparison are output to stdout, and a .tess file is generated
      (`difftermxx.tess`, where `xx` is the subdivision level)
      containing the position of the points where both termini
      differ.

    If the integration region is not periodic, then methods "color"
    and "qtree" are not defined. There are two possible options,
    controlled by the "KILLEXT" and "NOKILLEXT" keywords. If KILLEXT
    is active (the default behavior), the gradient path tracing is
    killed whenever it leaves the integration region, independently of
    the GRADIENT\_MODE being used. The terminus is then assigned to an
    "unknown" state, and the tetrahedra it generates are not
    integrated. If NOKILLEXT is active, the gradient path is continued
    as a "full gradient", until the terminus is found.

    The default is KILLEXT because, if the integration region is not
    periodic, the integral over atoms that are partially contained in
    it is most likely not meaningful to the user.

    The ODE integration method can be chosen using the
    QTREE\_ODE\_MODE keyword, that can assume the following values:

    - 1 : Euler method, fixed step, 1st order.

    - 2 : Heun method, fixed step, 2nd order.

    - 3 : Kutta method, fixed step, 3rd order.

    - 4 : Runge-Kutta method, fixed step, 4th order.

    - 5 : Euler-Heun embedded method, adaptive step. 1st order with
      2nd order error estimation. 2 evaluations per step.

    - 6 : Bogacki-Shampine embedded method, adaptive step. 3rd order
      with 5th order error estimation. The FSAL (first step also
      last) allows only 4 evaluations per step. Local
      extrapolation.

    - 7 : Runge-Kutta Cash-Karp embedded method, adaptive step. 4th
      order with 5th order error estimation. 6 evaluations per
      step.

    - 8 : Dormand-Prince 4-5 embedded method, adaptive step. 4th order
      with 5th order error estimation. 6 evaluation per step, with
      FSAL. Local extrapolation.

    For embedded methods (4-8), the absolute error requested from the
    method can be set using the ODE\_ABSERR keyword. The default of
    this variable is chosen so that reasonable stepsizes are
    kept. This default is 1d-3 for Euler-Heun and 1d-4 for the rest.
    Typically, methods with greater accuracy (7 and 8) save
    evaluations by increasing stepsize to values much larger than
    their lower accuracy counterparts.

    The step size of the fixed step methods (1-4) is controlled with
    the STEPSIZE keyword (bohr). In the variable step methods (5-8)
    the value of STEPSIZE is the length of the starting step.

    The default QTREE\_ODE\_MODE is Dormand-Prince (8).

11. When the integration of the base tetrahedron is finished, the
    termini of the grid points located at each of its four faces are
    copied to the corresponding neighboring IWST, if DOCONTACTS is
    active.

Once the integration of the IWST is completed, the atomic properties
are summed. The final result is output, together with an analysis of
the contribution of each subdivision level to the total integrated
properties.

It is possible to plot the basins obtained by QTREE using the
PLOT\_MODE keyword. It can assume the values:

- 0: no plotting is done.

- 1: a single tess file is written containing a description of the
  unit cell CPs, the IWS, and balls corresponding to all the grid
  points that have been sampled.

- 2: same as 1 but only the balls that are either on the face of an
  IWST or close to a IAS are output.

- 3: the full WS cell

- 4: a file for the full WS cell and several files, containing a
  description of each of the integrated basins. Note that the basins
  need not be connected.

- 5: same as 4 but only balls belonging to faces of IWST and IAS are
  output.

The default value is 0. If PLOT\_MODE is > 0, then the sticks that
form the tetrahedra are written to `.stick` files. The PLOTSTICKS and
NOPLOTSTICKS keywords control this behavior.

### Additional Considerations

* The integration of the volume is not done using the beta-sphere /
  basin separation because the volume of each tetrahedron is exactly
  known. The integrated cell volume for a periodic integration region
  will always be exact (if it is not, then it is an error). The
  integrated cell charge, on the contrary, is a measure of how well
  the tetrahedra are being integrated, but not of how well the IAS
  is being determined.

* For very high levels of QTREE (say 10-11, depending on the amount of
  memory your computer has), memory usage may be a problem. The
  COLOR\_ALLOCATE keyword controls the amount of memory allocated for
  the color and property arrays. The syntax is:
~~~
COLOR_ALLOCATE {0|1}
~~~
  Using a zero value, the color vector (and possibly the properties
  vectors, depending on PROP\_MODE) is allocated only for the current
  IWST. This saves memory but makes the computation slower, especially
  if the GRADIENT\_MODE is 2 or 3. In addition, setting COLOR\_ALLOCATE
  to 0 deactivates the passing of colors through the contacts
  (DOCONTACTS and NOCONTACTS keywords) and the plotting (sets
  PLOT\_MODE to 0). If COLOR\_ALLOCATE is 1, the color (and optionally
  the properties) of all the IWST are saved. By default,
  COLOR\_ALLOCATE is 1 if `maxlevel.i` is less than 9 and 0 if the
  maximum level is higher.

## Yu and Trinkle Grid Atomic Integration Method (YT) {#c2-yt}

~~~
YT [NNM] [NOATOMS] [WCUBE] [BASINS [OBJ|PLY|OFF] [ibasin.i]] [RATOM ratom.r]
   [DISCARD expr.s] [JSON file.json] [ONLY iat1.i iat2.i ...]
~~~
The [Yu and Trinkle (YT) method](https://doi.org/10.1063/1.3553716)
calculates the attraction basins of the reference fields and computes
QTAIM integrals in them. The reference field must be defined on a
grid. Hence YT will not work directly with wien2k, elk,
aiPI,... densities, although those can be transformed into a grid by
appropriate use of the LOAD keyword. The algorithm proceeds by running
over grid nodes in decreasing order of density. If a point has no
neighboring points with higher density, then it's a local maximum. If
it does, but all of them belong to the same basin then that point
belongs to the interior of that basin as well. Otherwise, it is
sitting on top of the interatomic surface. The actual fraction of a
grid point on an IAS belonging to a particular basin is calculated by
evaluating the trajectory flow to neighboring points.

The YT algorithm is described in
[J. Chem. Phys. 134 (2011) 064111](https://doi.org/10.1063/1.3553716)
which should be consulted for further details. Please, cite this
reference if you use this keyword in your work.

The located maxima in the field are identified by default with the
closest nucleus. If non-nuclear maxima are expected, use the NNM
keyword to assign only maxima that are only within 1 bohr of the
closest atom. This distance can be changed using the RATOM keyword
(`ratom.r` in bohr (crystals) or angstrom (molecules)), which also
controls the distance below which two maxima are considered
the same maximum. Changing the default `ratom.r` using the RATOM
keyword automatically activates the detection of NNM. The NOATOMS
option is appropriate for scalar fields where the maxima are not
expected to be at the atomic positions (or at least not all of
them). If NOATOMS is used, all the maxima found are given as NNM. This
is useful for fields such as the ELF, the Laplacian, etc.

The WCUBE option makes critic2 write cube files for the integration
weights of each attractor. In YT, these weights are values zero
(outside the basin), one (inside the basin), or some intermediate
value near the atomic basin boundary. The generated cube files have
names `<root>_wcube_xx.cube`, where `xx` is the attractor number.

Use the BASINS option to write a graphical representation of the
calculated basins. The format can be chosen using the OBJ, PLY, and
OFF keywords (default: OBJ). If an integer is given after the format
selector (`ibasin.i`), then plot only the basin for that
attractor. Otherwise, plot all of them. The basin surfaces are colored
by the value of the reference field, in the default gnuplot scale.

Any maxima that is not assigned to an existing atom or non-nuclear
critical point is automatically added to the
[critical point list](/critic2/manual/cpsearch/). It is possible to
get more information about these maxima by using the
[CPREPORT](/critic2/manual/cpsearch/#c2-cpreport) keyword.

In some cases, particularly if there is a vacuum region in your system
(for instance, if your system is a molecule or a surface), multiple
spurious maxima will appear due to numerical noise in the grid
values. The number of spurious attractors will increase the
computational cost and serve little purpose, as the vacuum region will
integrate to zero anyway. In these cases, the DISCARD keyword can be
used to make critic2 ignore any attractor that matches the expression
`expr.s` when it is evaluated at that point. For instance, if the
electron density is given by field `$rho`, and we want to discard
low-density critical points, we could use:
~~~
DISCARD "$rho < 1e-7"
~~~
The arithmetic expression may include any field, not just the
reference field.

By using the JSON keyword a JavaScript Object Notation (JSON) file is
created containing the molecular or crystal structure, information
about the reference field and the results of the integration.

The ONLY keyword restricts the integration to only the atoms given by
the user. The integers `iat1.i`,... are maxima identifiers from the
complete critical point list, which contains all the critical points
in the unit cell. For atoms, the identifiers in the complete CP list
coincide with those in the complete atom list, so `iat.1` for an atom
is also the integer ID from the list of atoms in the unit cell (see
the [notation](/critic2/manual/inputoutput/#c2-notation)). The ONLY
keyword is useful for delocalization index calculations in large
systems, where restricting the integration to only a handful of atoms
saves computing time.

Not all the properties defined by the
[INTEGRABLE](/critic2/manual/integrate/#c2-integrable) keyword are
integrated inside the basins. Only the subset of those properties that
are grids, have F or FVAL as the integrand and are congruent with the
reference grid are considered. This limitation can be circumvented by
using LOAD AS to define grid from fields or expressions that are not
given on a grid. In addition, no core is used even if the CORE keyword
is active. The volume is always integrated. A xyz file (`<root>_yt.xyz`)
is always written, containing the unit cell description (with border,
see WRITE) and the position of the maxima, labeled as XX.

Note that in the output ("List of basins and local properties"),
"Charge" refers not to the integrated electron density (because
critic2 does not know whether a given field is an electron density or
not) but to the value of the integral of the reference field in its
own basins (which may not make much sense if you are integrating, for
instance, the ELF or the Laplacian). Loading a second field and using
INTEGRABLE and the field number is the way to go in such cases.

Usage of the YT algorithm for grid fields is strongly recommended, as
it is much more efficient, robust and accurate than the alternatives.

For examples see the [calculation of Bader atomic
properties](/critic2/examples/example_11_01_simple-integration/#c2-bader).

## Henkelman et al. Grid Atomic Integration Method (BADER) {#c2-bader}

The algorithm by Henkelman et al. is implemented in critic2, and can
be used with the BADER keyword:
~~~
BADER [NNM] [NOATOMS] [WCUBE] [BASINS [OBJ|PLY|OFF] [ibasin.i]] [RATOM ratom.r]
      [DISCARD expr.s] [JSON file.json] [ONLY iat1.i iat2.i ...]
~~~
The BADER algorithm uses the reference field to calculate the QTAIM
basins. This field must be defined on a grid. BADER assigns grid nodes
to basins using the near-grid method incrementally described in
[Comput. Mater. Sci. 36, 354-360 (2006)](https://doi.org/10.1016/j.commatsci.2005.04.010),
[J. Comput. Chem. 28, 899-908 (2007)](https://doi.org/10.1002/jcc.20575), and
[J. Phys.: Condens. Matter 21, 084204 (2009)](https://doi.org/10.1088/0953-8984/21/8/084204).
Please, cite these references if you use this method. The same method
is implemented in the
[bader program](http://theory.cm.utexas.edu/henkelman/code/bader/)
available from the Henkelman research group's website.

The output and the optional keywords have the same meaning as in
YT, except in the following cases:

- The WCUBE option makes critic2 write cube files for the integration
  weights of each attractor. In BADER, these weights are either zero
  (outside the basin) or one (inside the basin). The generated cube
  files have names `<root>_wcube_xx.cube`, where `xx` is the attractor
  number. In addition, a file `<root>_wcube_all.cube` is written, in
  which the value at every grid point corresponds to the integer ID of
  the attractor for the corresponding basin.

Using BADER as an alternative to YT is recommended in very large
grids because of its more efficient memory usage, but in general it
gives less accurate integrations (at least in my experience).

For examples see the [calculation of Bader atomic
properties](/critic2/examples/example_11_01_simple-integration/#c2-bader).

## Isosurface Grid Integration (ISOSURFACE) {#c2-isosurface}

The ISOSURFACE keyword is used to integrate the volume and various
scalar fields in regions delimited by isosurfaces of a scalar
field. The isosurfaces are determined by the reference field and a
contour (isosurface) value set by the user. This keyword can be used
only if the reference field is given on a grid. The syntax is:
~~~
ISOSURFACE {HIGHER|LOWER} isov.r [WCUBE] [BASINS [OBJ|PLY|OFF] [ibasin.i]]
           [DISCARD expr.s]
~~~
The ISOSURFACE keyword must be followed by either HIGHER or LOWER and
a real number (`isov.r`). The regions to be integrated are bound by
isosurfaces with contour value equal to `isov.r`. If HIGHER, integrate
the regions with reference field value higher than `isov.r`. If LOWER,
integrate the regions with value lower than `isov.r`.

The WCUBE option is similar to the same option in
[BADER](/critic2/manual/integrate/#c2-bader). WCUBE writes cube files
containing the isosurface regions. The points in cube file
`<root>_wcube_xx.cube` have a value of 1, if the point is inside a
domain for isosurface `xx`, or 0 otherwise. Each point in the
`<root>_wcube_all.cube` file has grid value equal to the integer ID of
the isosurface domain in which it is contained, or 0 if it is not
inside any isosurface domain.

The BASINS keyword writes a graphical representation of the calculated
isosurface domains. The format can be chosen using the OBJ, PLY, and
OFF keywords (default: OBJ). If an integer is given after the format
selector (`ibasin.i`), then plot only the domain for that isosurface
ID. Otherwise, plot all of them.

The DISCARD keyword accepts an arithmetic expression that is evaluated
at every grid point. If the expression is true at a given grid point,
then that point is not considered in the construction of the
isosurface domains. This is useful in combination with critic2's
[structural variables](/critic2/manual/arithmetics/#c2-structvar). For
instance, doing:
~~~
DISCARD "@idnuc != 3"
~~~
represents only the isosurface associated with atom number 3.

## Hirshfeld Atomic Properties (HIRSHFELD) {#c2-hirshfeld}

The [Hirshfeld property](https://doi.org/10.1063/1.2831900)
associated with atom A for scalar field $$f({\bf r})$$ is calculated
as:

$$
\begin{equation}
F_{\rm A} = \int \frac{\rho_{\rm A}({\bf r})}{\rho_{\rm pro}({\bf r})}
\times f({\bf r}) d{\bf r}
\end{equation}
$$

where $$\rho_{\rm pro}$$ is the promolecular density (the sum of
atomic densities) and $$\rho_{\rm A}$$ is in-vacuo atomic density for
atom A. For instance, the Hirshfeld charge of atom A would be
$$Z_{\rm A} - N_{\rm A}$$, where $$N_{\rm A}$$ is the Hirshfeld atomic
integral for the all-electron density.

Hirshfeld atomic properties can be calculated using the HIRSHFELD
keyword:
~~~
HIRSHFELD [WCUBE] [ONLY iat1.i iat2.i ...]
~~~
There are two ways in which this keyword operates. If the reference
field is a grid, then a grid integration of all properties defined as
[INTEGRABLE](/critic2/manual/integrate/#c2-integrable) is carried
out. In the case of a grid integration, the input and output of
HIRSHFELD resembles that of [YT](/critic2/manual/integrate/#c2-yt) and
[BADER](/critic2/manual/integrate/#c2-bader). For this grid
integration to be numerically sensible, the integrable scalar field
must be reasonably smooth (for instance, the pseudo-density coming
from a plane-wave calculation).

The WCUBE option writes cube files for the Hirshfeld weights of each
nucleus, with file names `<root>_wcube_xx.cube`, where `xx` is the
atom identifier from the complete atom list. The ONLY keyword
restricts the integration to only certain atoms, given by their
identifiers from the complete atom list.

If the reference field is not a grid, HIRSHFELD assumes the field
contains the electron density and carries out the integration using
the selected [molecular mesh](/critic2/manual/misc/#c2-meshtype). In
this case, WCUBE and ONLY cannot be used, and only the volume and the
Hirshfeld atomic electron populations are calculated.

For an example see the [calculation of Hirshfeld
charges](/critic2/examples/example_11_01_simple-integration/#c2-hirshfeld).

## Voronoi Atomic Properties (VORONOI) {#c2-voronoi}

The [Voronoi property](https://doi.org/10.1002/jcc.10351)
associated with atom A for scalar field $$f({\bf r})$$ is calculated
as the integral of $$f$$ over all points in space that are closer to A
than to any other atom (the Voronoi region). For instance, the
[Voronoi deformation density](https://doi.org/10.1002/jcc.10351) (VDD)
of an atom is minus the integral of the electron density minus the
promolecular density over its Voronoi region.

Voronoi atomic properties can be calculated using the VORONOI keyword:
~~~
VORONOI [BASINS [OBJ|PLY|OFF] [ibasin.i]] [ONLY iat1.i iat2.i ...]
~~~
Only the grid integration case has been implemented, so the reference
field must be a grid. All properties defined as
[INTEGRABLE](/critic2/manual/integrate/#c2-integrable) are integrated
in the Voronoi regions, and the input and output of
VORONOI resembles that of [YT](/critic2/manual/integrate/#c2-yt) and
[BADER](/critic2/manual/integrate/#c2-bader). For this grid
integration to be numerically sensible, the integrable scalar field
must be reasonably smooth (for instance, the pseudo-density coming
from a plane-wave calculation).

The BASINS option writes a graphical representation of the Voronoi
regions. The format can be chosen using the OBJ, PLY, and OFF keywords
(default: OBJ). If an integer is given after the format selector
(`ibasin.i`), then plot only the basin for that atom (identifier from
the complete atom list). Otherwise, plot all of them. The basin
surfaces are colored by the value of the reference field, in the
default gnuplot scale. The ONLY keyword restricts the integration to
only certain atoms, given by their identifiers from the complete atom
list.

For an example see the [calculation of Voronoi deformation density
(VDD) charges](/critic2/examples/example_11_01_simple-integration/#c2-voronoi).

## Examples

+ [Integration of atomic charges and volumes with grids](/critic2/examples/example_11_01_simple-integration/)

+ [Calculation of delocalization indices in solids](/critic2/examples/example_11_10_deloc-indices/)
