<p align="center">
  <img src="https://github.com/aoterodelaroza/critic2/blob/master/dat/logo/critic2_logo-big-with-letters.png?raw=true" alt="critic2 logo" height="30%" width="30%"/>
</p>

**Critic2** is a program for the manipulation and analysis of structural
and chemical information in molecules and periodic solids. Critic2 can
be used to:

* Read, transform, and create molecular and crystal structures.

* Carry out crystallographic computations (environments, coordination
  numbers, structural comparison,...).

* Read, analyze, manipulate, combine, and create scalar fields such as
  the electron density or the ELF.

* Carry out calculations using Bader's atoms in molecules theory
  (finding critical points, integrating atomic basins, plotting
  gradient path manifolds,...).

* Make non-covalent interaction (NCI) and other similar plots.

Critic2 is provides an abstraction layer on top of the underlying
quantum chemical calculation. Critic2 interfaces with many electronic
structure programs: WIEN2k, elk, PI, Quantum ESPRESSO, abinit, VASP,
DFTB+, Gaussian, psi4, siesta, and more.

## Compilation and installation

Critic2 can be compiled on Linux and macOS. For this, you
will need a Fortran and a C compiler, and a build system (either
autotools or cmake, and make). Detailed
[installation instructions](https://aoterodelaroza.github.io/critic2/installation/)
can be found in the manual. Critic2 uses some fairly modern Fortran
features, which may not be implemented on all current Fortran
compilers. Please, check out the
[relevant section of the manual](https://aoterodelaroza.github.io/critic2/installation/#whichcompilerswork)
to see if your compiler is listed.

Critic2 is parallelized for shared-memory architectures (unless
specifically deactivated during the build process). You change the
number of parallel threads by setting the <code>OMP_NUM_THREADS</code>
environment variable.

The environment variable CRITIC_HOME is necessary if critic2 was not
installed with `make install`. It must point to the root directory of
the distribution:

	export CRITIC_HOME=/home/alberto/programs/critic2dir

This variable is necessary for critic2 to find the atomic densities,
the cif dictionary and the library data. These should be in
`${CRITIC_HOME}/dat/`.

Lastly, a number of
[external libraries](https://aoterodelaroza.github.io/critic2/installation/#external-libraries)
can be used to extend critic2's capabilities, including readline
(shell-like features in the critic2 command line),
[Libxc](https://gitlab.com/libxc/libxc) (exchange-correlation energies
and potentials), and [Libcint](https://github.com/sunqm/libcint)
(molecular integrals).

## Using critic2

All of critic2's features are documented in the  [reference manual](https://aoterodelaroza.github.io/critic2/).
For a text version of the manual, please clone the [repository](https://github.com/aoterodelaroza/aoterodelaroza.github.io).
Some examples are provided in the [examples section](https://aoterodelaroza.github.io/critic2/examples/)
of the documentation.

## References and citation

The basic references for critic2 are:

* [A. Otero-de-la-Roza, E. R. Johnson and V. Luaña, Comput. Phys. Commun. **185**, 1007-1018 (2014)](http://dx.doi.org/10.1016/j.cpc.2013.10.026)
* [A. Otero-de-la-Roza, M. A. Blanco, A. Martín Pendás and V. Luaña, Comput. Phys. Commun. **180**, 157–166 (2009)](http://dx.doi.org/10.1016/j.cpc.2008.07.018)

The output and the manual may contain additional references pertaining
to methods employed by particular keywords.

## License

Critic2 is made available under the GNU/GPL v3 license. See the
[LICENSE](https://github.com/aoterodelaroza/critic2/blob/master/LICENSE)
file in the root of the critic2 distribution for more details.
