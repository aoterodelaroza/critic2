
## Creating STM plots (STM)

The STM keyword generates plots comparable to those obtained in
scanning tunneling microscopy (STM) experiments:
~~~
STM [CURRENT [curr.r]|HEIGHT [hei.r]] [TOP top.r]
    [{CELL|CELLS} nx.i ny.i] [NPTS n1.i n2.i]
    [LINE x0.r y0.r x1.r y1.r npts.i]
~~~
The STM keyword should be applied only to systems meeting the
following requirements:

* The system is a slab, possibly with molecules adsorbed onto one or
  both faces.

* The slab and the corresponding vacuum are perpendicular to one of
  the crystallographic axes (a, b, or c). The two angles related to
  the perpendicular axis must be right. For instance, a typical
  situation is one where a and b form a hexagonal cell in the plane,
  and c is perpendicular to the slab. The angles related to c (alpha
  and beta) must both have 90 degrees.

* The slab has two faces: one with vacuum above and one below. Only
  the face with vacuum above it will be used (i.e. critic2 will scan
  in the positive vacuum direction).

* The STM plot calculation uses the reference field. 

The STM keyword uses the popular 
[Tersoff-Hamann approximation](https://doi.org/10.1103/PhysRevB.31.805), in
which the observed current is proportional to the local density of
states at the Fermi level:

$$
\begin{equation}
I \approx V \rho_{\rm loc}({\bf r},V)
\end{equation}
$$

where $$I$$ is the current, $$V$$ is the bias voltage, and 
$$\rho_{\rm loc}$$ is the local density of states at $$E_F$$:

$$
\begin{equation}
\rho_{\rm loc}({\bf r},V) = \sum_{ {\bf k},n}^{E_F-eV\to E_F} |\psi_{ {\bf k},n}({\bf r})|^2
\end{equation}
$$

where the sum runs only over the one-electron states that have
energies between $$E_F-eV$$ and $$E_F$$.

To use the STM keyword, the reference field must be the local density
of states at the Fermi level $$\rho_{\rm loc}$$ calculated by your
quantum chemistry package of choice. For instance, in Quantum
ESPRESSO, a Gaussian cube file containing the LDOS can be obtained by
using `plot_num=5` and a `sample_bias` equal to the bias voltage (in
Ry!).

There are two main modes of operation in STM: constant current and
constant height. In constant current mode, the tip is allowed to move
vertically in a way that the current through it is constant. This
results in a map of tip height as a function of displacement over the
surface. In the constant height mode, the height is kept fixed and the
STM plot represents the current across the tip.

Critic2 can generate plots for both STM operation modes. In the
constant current mode, selected with the keyword CURRENT, the plot
represents the height of the tip relative to the top of the surface,
in angstrom. The height corresponds to a constant current given by an
LDOS equal to `curr.r` (default: 0.001 a.u.), so the plot is
essentially the height of the LDOS=`curr.r` isosurface as a function
of the displacement over the surface.

In CURRENT plots, the height is calculated relative to the top of the
surface, which is taken as the coordinate in the vacuum direction of
the last atom before the vacuum. The position of the top of the
surface, which affects only the length scale on the plot and not the
shape of the plot itself, can be changed with the TOP keyword. The
`top.r` value corresponds to the fractional coordinate in the vacuum
direction of the height to which the constant-current plot will be
referred. For instance, to refer the constant-current plot to the
center of the slab, use `TOP 0.5`. The CURRENT keyword is the default
if no CURRENT or HEIGHT is given.

The HEIGHT keyword selects the constant-height mode plot. In this mode
of operation the LDOS, which is proportional to the current across the
tip, is plotted at a constant height over the surface. The
`hei.r` option to HEIGHT is the fractional coordinate along the
perpendicular axis of the plot plane. By default, critic2 uses the
fractional coordinate corresponding to the last atom before the
vacuum plus one bohr into the vacuum.

The CELL (or CELLS) keyword controls the number of unit cells plotted
in each in-plane crystallographic direction (default: `1 1`). The
number of points plotted in each in-plane unit cell is given by the
NPTS keyword (`n1.i` and `n2.i`). By default, critic2 uses the number
of grid points in each of those directions if the reference field is a
grid, or 51x51 points if not. If a grid is used, letting critic2 use
the grid geometry (by not using NPTS) is strongly recommended because
it is cheaper and accurate.

On output, two files will be written. The `<root>_stm.dat` file
contains the 2D data for the STM signal on the selected plane. The
`<root>_stm.gnu` file is an example script that generates a plot
similar to those found in the literature.

The LINE keyword makes critic2 plot a line containing the STM
information (either height or current) instead of a plane. The line
goes from (`x0.r` `y0.r`) to (`x1.r` `y1.r`), in fractional
coordinates for the in-plane crystallographic axes. `npts.i` points
along the line are calculated. The output file is
`<root>_stm_line.dat`.

The current implementation of STM has benefited from the code and the
guidance kindly provided by Enrico Benassi (see the THANKS file).

## Examples

- [Making STM plots with Quantum ESPRESSO and critic2](/critic2/examples/example_14_01_stmqe/)

