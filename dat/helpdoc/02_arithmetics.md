
## Basic Usage {#c2-arithbasic}

In critic2, an arithmetic expression can be used almost everywhere in
the input where a real or integer number is expected. Arithmetic
expressions that appear in the input (without an associated keyword)
are evaluated and their result is written to the output. For instance,
you can start critic2 and write:
~~~
3+2*sin(pi/4)
%% 3+2*sin(pi/4)
   4.4142135623731
~~~
Similarly, variables can be defined and utilized in any
expression. Variable names must start with a letter and are composed
only of letters, numbers, and the underscore character. Also, they
cannot have the same name as a known constant (pi, e, eps) or a
function, regardless of case. Variables in expressions are
case-sensitive. To use a variable, first you need to assign it. For
instance:
~~~
a = 20+10
a/7 + 1
%% a/7 + 1
   5.2857142857143
~~~
In fact, by using the `-q`
[command-line option](/critic2/manual/#c2-commandline), critic2 can be
used as a simple calculator:
~~~
$ echo "1+erf($RANDOM/100000)" | critic2 -q
1.2578249310340
~~~
When used in combination with other keywords, arithmetic expressions
must be enclosed in either double quotes ("), single quotes ('), or
parentheses, or they must form a single word (i.e. no spaces). For
instance, this is valid critic2 input:
~~~
a = 0.12
[...]
Be 1/3 "2 / 3" 1/4+a
~~~
but this is not:
~~~
Be 1/3 2 /3 1/4+a
~~~
Arithmetic expressions can contain:

* Operators: `+`, `-`, `*`, `/`, `**` or `^`, `%` (modulo)

* [Functions](/critic2/manual/arithmetics/#availchemfun): the usual
  mathematical functions (`sin`, `exp`, etc.) as well as "chemical
  functions".

* Constants: `pi`, `e` and `eps` (the machine precision).

* Variables defined by the user as above.

* [Structural variables](/critic2/manual/arithmetics/#c2-structvar)

Parentheses can be used and the usual rules of associativity and
precedence apply.

In some cases, arithmetic expressions can be applied to make new
scalar fields by transforming existing fields. Scalar fields are
denoted by a dollar sign (`$`) followed by an identifier and,
optionally, a modifier separated by a colon (`$id:modifier`). The
default identifier for a field is the order in which the field was
loaded in the input. For instance, if the first field (`$1`) is the
spin-up density and the second field (`$2`) is the spin-down density,
the total density can be calculated with the expression `$1+$2` and
the spin density can be calculated with `$1-$2`. Field number zero
(`$0` or `$rho0`) represent the promolecular density, which is always
available once the crystal or molecule structure is known. In our
example, `$1+$2-$rho0` would represent the density difference between
the actual density and the sum of atomic densities. Fields can also be
referred by a name, if the
[ID keyword](/critic2/manual/fields/#c2-addload) is used. Named
fields simplify work when you have multiple fields.

It is possible to specify a field modifier right after the number or
name for that field in order to access its derivatives and other
properties related to it. The modifier is separated by the field name
using the ":" character and is case insensitive. It may be one of:

+ `v`: valence-only value of the field (it is usually employed to access
  the valence density in a grid field in which core augmentation is
  active).
+ `c`: core-only value of the field.
+ `x`, `y`, `z`: first derivatives: $$f_x$$, $$f_y$$, $$f_z$$.
+ `xx`, `xy`, `yx`, `xz`, `zx`, `yy`, `yz`, `zy`, `zz`: second
  derivatives: $$f_{xx}$$, $$f_{xy}$$, etc.
+ `g`: norm of the gradient, $$\lvert\nabla f\rvert$$.
+ `l`: Laplacian, $$\nabla^2 f$$.
+ `lv`: valence Laplacian (Laplacian without core augmentation).
+ `lc`: core Laplacian.

For instance, `$2:l` is the Laplacian of field 2 and `$rho0:xy` is the
xy-component of the Hessian of the promolecular density. In molecular
wavefunctions (wfn/wfx/fchk/molden), the value of particular orbitals
can be selected with the following field modifiers:

+ an integer (`<n>`): selects molecular orbital number `n`. For
  instance, `$1:3` refers to the value of molecular orbital number 3
  in field 1.
+ `HOMO`: the highest-occupied MO in an RHF wavefunction.
+ `LUMO`: the lowest-unoccupied MO in an RHF wavefunction.
+ `AHOMO`: the alpha highest-occupied MO in a UHF wavefunction.
+ `ALUMO`: the alpha lowest-unoccupied MO in a UHF wavefunction.
+ `BHOMO`: the beta highest-occupied MO in a UHF wavefunction.
+ `BLUMO`: the beta lowest-unoccupied MO in a UHF wavefunction.
+ `A<n>`: alpha MO number `n` in a UHF wavefunction.
+ `B<n>`: beta MO number `n` in a UHF wavefunction.

Using these flags to select virtual orbitals requires a file format
that contains information about them (wfn and wfx do not work). In
addition, the molecular wavefunction field has to be loaded using the
[READVIRTUAL keyword](/critic2/manual/fields/#c2-load). Also, in
unrestricted molecular wavefunction fields, one can access the spin
densities:

+ `up`: up-spin density.
+ `dn`: down-spin density.
+ `sp`: (difference) spin density.

In critic2, there is a distinction between expressions that reference
fields and those that do not and, for certain keywords, critic2 will
decide what to do with an expression based on this distinction. For
instance:
~~~
a = 2
CUBE CELL FIELD "a*$1"
~~~
calculates a grid spanning the entire cell for a scalar field that is
built as two times the value of field number `1`, but:
~~~
a = 2
CUBE CELL FIELD "a*1"
~~~
uses field number two (`a*1 = 2`) to calculate the same grid because
this expression does not reference any field. Expressions involving
fields can also be used in [LOAD](/critic2/manual/fields/#c2-load)
as well as in many other keywords
([POINT](/critic2/manual/graphics/#c2-point), [LINE](/critic2/manual/graphics/#c2-line),
[INTEGRABLE](/critic2/manual/integrate/#c2-integrable), etc.).
When using arithmetic expressions to create new fields, it is also possible
to refer to coordinates in real space in those expressions by using
[Structural variables](/critic2/manual/arithmetics/#c2-structvar).

## Clear Variables (CLEAR) {#c2-clear}

The value of a variable can be cleared using the CLEAR keyword:
~~~
CLEAR var1.s var2.s ...
CLEAR ALL
~~~
This keyword deletes the variables var1.s, var2.s, etc. or all the
variables (ALL).

## List Variables (LIST) {#c2-list}

At any moment, the internal list of variables can be printed to the
output using the keyword LIST:
~~~
LIST
~~~
The LIST keyword lists all named variables and fields.

## Special Fields

Some special fields are defined from the crystal (or molecular)
structure alone. For now, the only available special field is `ewald`,
that can be accessed using `$ewald` in arithmetic expressions, and
gives the value of the Ewald potential using the current existing
charges. Atomic charges can be set with the
[Q keyword](/critic2/manual/crystal/#c2-charge). For instance:
~~~
cube cell field "2 * $ewald"
~~~
calculates a grid using 2 times the value of the Ewald potential.

## List of Available Functions {#availchemfun}

Arithmetic expressions can use any of the functions in the critic2
function library. The list of functions includes the usual
mathematical functions (like `exp` or `sin`) but also functions that
are meant to be applied to scalar fields of a certain type (e.g., the
Thomas-Fermi kinetic energy density, `gtf`). The latter are called
"chemical functions", since they carry chemical information. It is
very important to distinguish whether a function expects a numerical
argument (e.g. `sin(x)`), or a field identifier (e.g. `gtf(1)` or
`gtf(rho0)`).

The list of arithmetic functions is: `abs`, `exp`, `sqrt`, `floor`,
`ceil`, `ceiling`, `round`, `log`, `log10`, `sin`, `asin`, `cos`,
`acos`, `tan`, `atan`, `atan2`, `sinh`, `cosh`, `erf`, `erfc`, `min`,
`max`. All these functions apply to numbers (or other arithmetic
expressions), and their behavior is the usual. For instance,
`sin(2*pi)`, `max($1,0)`, and `atan2(y,x)` are all valid expressions.

The chemical functions in the lists below accept one or more field
identifiers as their arguments. Their purpose is to provide shorthands
to build fields from other fields using physically-relevant
formulas. For instance, `gtf(1)` is the Thomas-Fermi kinetic energy
density calculated using the electron density in field 1. In all
instances, `gtf(1)` is equivalent to writing the formula in full
(`3/10*(3*pi^2)^(2d0/3d0)*$1^(5/3)`) but, naturally, much more
convenient. Some of the chemical functions like, for instance, those
that require having access to the one-electron wavefunctions (e.g. the
ELF), can only be used with fields of a certain type. The name in
square brackets, if available, is a shorthand for applying the
chemical function to the reference field (case-insensitive) in the
POINTPROP keyword.

The following list of chemical functions can be used with any field
type. In all cases, field `id` should correspond to the system's
electron density (it is up to the user to make sure this is the
case).

* `gtf(id)` [`GTF`]: Thomas-Fermi kinetic energy density. The kinetic
  energy density for a uniform electron gas with its density given by
  the value of field id at the point ($$g^{\text{TF}}$$).[^yangparr]

* `vtf(id)` [`VTF`]: the potential energy density calculated using the
  Thomas-Fermi kinetic energy density and the local virial theorem
  ($$v^{\text{TF}}({\bf r}) = \frac{1}{4}\nabla^2\rho({\bf r}) - 2g^{\text{TF}}({\bf r})$$ in au).[^yangparr]

* `htf(id)` [`HTF`]: the total energy density calculated using the
  Thomas-Fermi kinetic energy density and the local virial theorem
  ($$h^{\text{TF}}({\bf r}) = g^{\text{TF}}({\bf r}) + v^{\text{TF}}({\bf r})$$).
  The field id must contain the electron density of the system.[^yangparr]

* `gtf_kir(id)` [`GTF_KIR`]: Thomas-Fermi kinetic energy density with the
  semiclassical gradient correction proposed by Kirzhnits for the
  not-so-homogeneous electron gas. The electron density and its
  derivatives are those of field id at every point in space.($$g^{\text{kir}}$$)
  [^kirzhnits1][^kirzhnits2][^abramov][^zhurova][^espinosa]

* `vtf_kir(id)` [`VTF_KIR`]: the potential energy density calculated using
  gtf_kir(id) and the local virial theorem ($$v^{\text{kir}}({\bf r})
  = \frac{1}{4}\nabla^2\rho({\bf r}) - 2g^{\text{kir}}({\bf r})$$ in au).
  [^kirzhnits1][^kirzhnits2][^abramov][^zhurova][^espinosa]

* `htf_kir(id)` [`HTF_KIR`]: the total energy density calculated using
  `gtf_kir(id)` and the local virial theorem
  ($$h^{\text{kir}}({\bf r}) = g^{\text{kir}}({\bf r}) +
  v^{\text{kir}}({\bf r})$$).
  [^kirzhnits1][^kirzhnits2][^abramov][^zhurova][^espinosa]

* `lag(id)` [`LAG`]: the Lagrangian density ($$-\frac{1}{4}\nabla^2\rho$$).

* `lol_kir(id)` [`LOL_KIR`]: the localized-orbital locator (LOL) with
  the kinetic energy density calculated using the Thomas-Fermi
  approximation with Kirzhnits gradient correction.[^tsirelson]

* `rdg(id)` [`RDG`]: the reduced density gradient,

\begin{equation}
s({\bf r}) = \frac{\lvert{\bf \nabla}\rho({\bf r})\rvert}{2(3\pi^2)^{1/3}\rho({\bf r})^{4/3}}
\end{equation}

The following functions require the kinetic energy density, and
therefore can only be used with fields that provide the one-electron
wavefunctions. At present, this is only available for molecular
wavefunction fields.

* `gkin(id)` [`GKIN`]: the kinetic energy density, G-version
  ($$\sum_i{\bf \nabla}\psi_i\cdot{\bf \nabla}\psi_i /2$$).[^bader1][^bader2]

* `kkin(id)` [`KKIN`]: the kinetic energy density, K-version
  ($$\sum_i\psi_i\nabla^2\psi_i /2$$).[^bader1][^bader2]

* `vir(id)` [`VIR`]: the electronic potential energy density, also
  called the virial field.[^keith]

* `he(id)` [`HE`]: the electronic energy density, `vir(id) + gkin(id)`.

* `elf(id)` [`ELF`]: the electron localization function (ELF).[^becke1]

* `lol(id)` [`LOL`]: the localized-orbital locator (LOL).[^schmider1][^schmider2]

* `brhole_a1(id)`, `brhole_a2(id)`, `brhole_a(id)`: the $$A$$
  prefactor of the spherically averaged hole model proposd by Becke
  and Roussel (spin up, down, and average, respectively). The BR hole
  is an exponential $$Ae^{-\alpha r}$$ at a distance $$b$$ from the
  reference point.[^becke2]

* `brhole_b1(id)`, `brhole_b2(id)`, `brhole_b(id)`: the $$b$$ parameter of
  the BR hole model (spin up, down, and average). $$b$$ is distance from
  the exponential center to the reference point.[^becke2]

* `brhole_alf1(id)`, `brhole_alf2(id)`, `brhole_alf(id)`: the exponent
  of the BR hole model (spin up, down, and average).[^becke2]

* `xhcurv1(id)`, `xhcurv2(id)`, `xhcurv(id)`: the curvature of the
  exchange hole at the reference point (spin up, down, and
  average). $$Q_{\sigma}$$ in the literature.[^becke2]

* `dsigs1(id)`, `dsigs2(id)`, `dsigs(id)`: the leading coefficient of
  the same-spin pair density (spin up, down, and
  average). $$D_\sigma$$ in the literature.[^becke2]

The following chemical functions require both a molecular wavefunction
and basis set information (at present, this can only be read from a
Gaussian fchk file). In addition, it is necessary to have critic2
compiled with the [libcint library](/critic2/installation/#c2-libcint) to
calculate the molecular integrals involved.

* `mep(id)`: molecular electrostatic potential.

* `uslater(id)`: Slater potential $$U_x$$. The HF exchange energy is
  $$\int\rho({\bf r})U_x({\bf r})d{\bf r}$$.[^becke3]

* `nheff(id)`: reverse BR efefctive hole normalization.[^becke3]

* `xhole(id,x,y,z)`: Exchange hole with reference point at (x,y,z). The
  coordinates are Cartesian in angstrom referred to the molecular
  origin if the system is a molecule or crystallographic coordinates
  if the system is a crystal.

Other special labels can be used, that activate the calculation of
properties for the reference field. These are:

* `stress`: calculate the Schrodinger stress tensor of the reference
  field. The virial field is the trace of this tensor.[^keith]

A particular case of chemical function is `xc()`, that allows the user
to access the [libxc library](/critic2/installation/#c2-libxc). This is
only possible if the libxc library was linked in the compilation
of critic2.

## Use of LIBXC in Arithmetic Expressions {#libxc}

If critic2 is [linked to the libxc library](/critic2/installation/#c2-libxc)
then the xc() chemical function can be used in arithmetic
expressions. xc() calculates the exchange and/or
correlation energy density for one of the functionals in the libxc
library. The number of arguments to xc() depends on the type of
functional invoked, which is selected using an integer index. The list
of functionals available and their corresponding indices should be
consulted in the libxc documentation. The integer index that selects
the functional always appears **last** in the calling sequence of
xc().

The number and type of arguments to `xc(...,idx)` depend on the type
of functional specified by the `idx` integer. This integer depends on
the libxc version used but a full list can be obtained within critic2
using the [LIBXC](#c2-libxc) keyword. The call to `xc()` can be one
of:

* LDA functional: `xc(rho,idx)`

* GGA functional: `xc(rho,grad,idx)`

* meta-GGA functional: `xc(rho,grad,lapl,tau,idx)`

where `rho` is the electron density ($$\rho({\bf r})$$), `grad` is its
gradient ($$\lvert\nabla\rho({\bf r})\rvert$$), `lapl` is its Laplacian
($$\nabla^2\rho({\bf r})$$) and `tau` is the kinetic energy density
($$1/2\sum_i\lvert\nabla\psi_i\rvert^2$$). Note that
`rho`, `grad`, `lapl`, and `tau` are **expressions**, not field
identifiers as in the chemical functions above or in the `idx`
argument. For instance, the expression for LDA using the electron
density loaded in field number 1 would be:
~~~
xc($1,1)+xc($1,9)
~~~
because `idx=1` is Slater's exchange and `idx=9` is Perdew-Zunger
correlation. PBE (a GGA functional) would be:
~~~
xc($1,$2,101)+xc($1,$2,130)
~~~
Here `idx=101` is PBE exchange and `idx=130` is PBE correlation. Field
`$1` contains the electron density and `$2` is its gradient. In fields
given on a grid spanning the crysatl unit cell, the gradient field can
be calculated using:
~~~
LOAD AS GRAD 1
~~~
which is calculated using fast Fourier transform. Alternatively, if
the field is not a grid, the gradient is best calculated directly
using the `:g` field modifier:
~~~
xc($1,$1:g,101)+xc($1,$1:g,130)
~~~

### List available LIBXC functionals (LIBXC) {#c2-libxc}

When critic2 is compiled with [libxc](/critic2/installation/#c2-libxc)
support, the LIBXC keyword can be used to query the library for a
list of available functoinals:
~~~
LIBXC [REF|REFS] [NAME|NAMES] [FLAGS] [ALL]
~~~
By default, the LIBXC keyword gives a list of the functional IDs
in the first column, followed by the functional name, kind, and
family. If FLAGS is used, then flags containing additional information
about the functional (e.g. if it provides the exchange-correlation
energy, etc.) are given. If NAME or NAMES is used, the full name of
the functional is printed. If REF or REFS is used, the literature
references are also printed. The ALL option is equivalent to using
REFS, NAMES, and FLAGS all at once.

## List of Structural Variables {#c2-structvar}

When creating new fields as transformations of existing fields, it is
possible to use the crystal or molecular structure (the nearest atom,
the coordinates, etc.) as part of the transformation using structural
variables. Structural variables start with the symbol "@" followed by
an identifier that selects the type of variable. For some of the
structural variables, an additional modifier can be applied using the
":" symbol after the variable identifier. For instance, `@dnuc` gives
the distance to the nearest nucleus. `@rho0nuc:3` is the atomic
density at the given point of the nearest atom number 3 (from the
complete list). The following structural variables are accepted in
critic2:

* `dnuc`: Distance to the closest nucleus. By default, the distance is
  in angstrom for molecules and in bohr for crystals, unless changed
  using the [UNITS](/critic2/manual/inputoutput/#c2-units) keyword.[^struvar]

* `xnucx`, `ynucx`, `znucx`: x,y,z coordinates of the nearest nucleus
  in crystallographic coordinates.[^struvar]

* `xnucc`, `ynucc`, `znucc`: x,y,z coordinates of the nearest nucleus
  in Cartesian coordinates. By default, the coordinates have units of
  bohr in crystals, and are referred to the molecular center and have
  units of angstrom in molecules.[^struvar]

* `xx`, `yx`, `zx`: the x,y,z coordinates of the point where the
  arithmetic expression is being evaluated (crystallographic
  coordinates).

* `xc`, `yc`, `zc`: the x,y,z Cartesian coordinates of the point where
  the arithmetic expression is being evaluated. Units are bohr.

* `xm`, `ym`, `zm`: the x,y,z Cartesian coordinates of the point where
  the arithmetic expression is being evaluated. By default, the
  coordinates have units of bohr in crystals, and are referred to the
  molecular center and have units of angstrom in molecules.

* `xxr`, `yxr`, `zxr`: the x,y,z coordinates of the point where the
  arithmetic expression is being evaluated (crystallographic
  coordinates in the reduced unit cell).

* `idnuc`: complete-list ID of the closest nucleus.[^struvar]

* `nidnuc`: non-equivalent list ID of the closest nucleus.[^struvar]

* `rho0nuc`: atomic density contribution from the nearest nucleus.[^struvar]

* `spcnuc`: species ID of the nearest nucleus.[^struvar]

* `zatnuc`: atomic number of the closest nucleus.[^struvar]

For instance, if we have a crystal structure and its electron density
is grid field `$1`, then:
~~~
LOAD AS "(@idnuc == 1) * $1" ID voronoi
SUM VORONOI
~~~
calculates the Voronoi charge of atom 1. Likewise,
~~~
LOAD AS "@rho0nuc:1/$0 * $1" ID hirsh
SUM hirsh
~~~
calculates the Hirshfeld charge (although it is probably more
convenient to use the [HIRSHFELD](/critic2/manual/misc/#c2-hirshfeld)
keyword. Structural variables are also useful
in molecules in combination with the
[MOLCALC](/critic2/manual/misc/#c2-molcalc)
keyword. For instance, to calculate the dipole moment of a neutral
molecule (in units of electrons*angstrom):
~~~
molcalc "$wfx * @xc"
molcalc "$wfx * @yc"
molcalc "$wfx * @zc"
~~~
Similar expressions can also be used to create new scalar fields by
restricting or modifying the values of a scalar field only in certain
areas of the system.

[^yangparr]: Yang and Parr, Density-Functional Theory of Atoms and Molecules.
[^kirzhnits1]: Kirzhnits, (1957). Sov. Phys. JETP, 5, 64-72.
[^kirzhnits2]: Kirzhnits, Field Theoretical Methods in Many-body Systems (Pergamon, New York, 1967).
[^abramov]: Abramov, Y. A. Acta Cryst. A (1997) 264-272.
[^zhurova]: Zhurova and Tsirelson, Acta Cryst. B (2002) 58, 567-575.
[^espinosa]: Espinosa et al., Chem. Phys. Lett. 285 (1998) 170-173.
[^tsirelson]: Tsirelson and Stash, Acta Cryst. (2002) B58, 780.
[^bader1]: Bader and Beddall, J. Chem. Phys. (1972) 56, 3320.
[^bader2]: Bader and Essen, J. Chem. Phys. (1984) 80, 1943.
[^keith]: Keith et al. Int. J. Quantum Chem. (1996) 57, 183-198.
[^becke1]: Becke and Edgecombe J. Chem. Phys. (1990) 92, 5397-5403.
[^schmider1]: Schmider and Becke, J. Mol. Struct. (Theochem) (2000) 527, 51-61.
[^schmider2]: Schmider and Becke, J. Chem. Phys. (2002) 116, 3184-3193.
[^becke2]: A.D. Becke and M.R. Roussel, Phys. Rev. A 39 (1989) 3761.
[^becke3]: A.D. Becke, J. Chem. Phys. 138 (2013) 074109.
[^struvar]: This structural variables accepts a modifier. A modifier is a colon (`:`) followed by a number `id.i`. If given, the modifier restricts the structural variable to atoms with integer ID `id.i` (from the complete list).
