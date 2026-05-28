# Iterative Hirshfeld (`HIRSHFELD_I`)

`HIRSHFELD_I` performs an iterative Hirshfeld (Hirshfeld-I) partitioning
of the reference field on a grid, following Bultinck, Van Alsenoy,
Ayers and Carbó-Dorca, *J. Chem. Phys.* **126**, 144111 (2007)
([10.1063/1.2715563](https://doi.org/10.1063/1.2715563)).

The motivation is well known: stockholder Hirshfeld weights built from
neutral free-atom densities systematically underestimate chemical
polarity, because they do not adapt to the actual charge state of the
atom in the molecule. Hirshfeld-I cures this by making the
reference density of each atom *self-consistent* with the population
it ends up carrying.

## Syntax

```
HIRSHFELD_I [WCUBE] [ONLY iat1.i iat2.i ...]
            [WFCDIR dir.s] [HITOL tol.r] [HIMAXIT maxit.i]
```

| Option         | Default | Meaning                                                                                                     |
| -------------- | ------- | ----------------------------------------------------------------------------------------------------------- |
| `WCUBE`        | off     | Write per-atom weight cubes (same as `HIRSHFELD`).                                                          |
| `ONLY i1 i2…`  | all     | Restrict atomic-property integration to the listed cell atoms.                                              |
| `WFCDIR dir`   | unset   | Look in `dir` for user-supplied charged `.wfc` files named `<elem>_q+N.wfc` / `<elem>_q-N.wfc`. See below.  |
| `HITOL tol`    | `1e-4`  | Convergence threshold on `max|ΔQ|` between successive SCF iterations.                                       |
| `HIMAXIT n`    | `60`    | Maximum number of SCF iterations.                                                                           |

The keyword is otherwise a drop-in replacement for `HIRSHFELD`: it
needs a grid-type reference field, integrates the same scalar
properties through `INTEGRABLE`, and reports the same volumes /
populations / Hirshfeld overlap populations.

## Algorithm

Let `ρ(r)` be the reference (molecular) electron density and let
`ρ_A^{q}(r)` be the spherically averaged density of atom *A* in
integer charge state *q*. The Hirshfeld-I weights are

```
            ρ_A^{q_A}(r)
w_A(r) = ───────────────────
          Σ_B ρ_B^{q_B}(r)
```

where each `q_A` is itself a function of the population that *A* picks
up under those weights:

```
N_A = ∫ w_A(r) ρ(r) dr        Q_A = Z_A − N_A
```

The reference density for a fractional `Q_A` is interpolated linearly
between the two flanking integer-charge densities,

```
ρ_A^ref(r) = (1−f) ρ_A^{⌊Q_A⌋}(r) + f ρ_A^{⌈Q_A⌉}(r)        f = Q_A − ⌊Q_A⌋
```

SCF loop:

1. Initialise every `Q_A = 0` (neutral references).
2. Pre-load the radial densities for every `(Z, q)` pair currently in
   use (see the resolution order below). This step is serial so that
   the next step can be embarrassingly parallel.
3. In a single OpenMP grid sweep, simultaneously rebuild `Σ_B ρ_B^ref`
   at every grid point and accumulate `N_A = ∫ w_A ρ dV`.
4. Update each `Q_A`, derive the bracketing `(⌊Q⌋, ⌈Q⌉)` and the
   mixing fraction `f`.
5. If `max|ΔQ|` < `HITOL` exit, else go to step 2.

The same iteration data structure is then handed to
`intgrid_hirshfeld_fields` and `intgrid_hirshfeld_overlap`, which use
the converged `ρ_A^ref` (instead of the neutral free-atom density)
wherever they previously called `agrid(z)%interp`.

## Charged atomic densities

Critic2's built-in `.wfc` tables only contain neutral atoms. Hirshfeld-I
needs both anionic and cationic references. They are resolved in this
order, per `(Z, q)`:

1. **User-supplied file (`WFCDIR`).** If `WFCDIR <dir>` is given, the
   code looks for
   - `<dir>/<elem>_q+N.wfc` for `q > 0`
   - `<dir>/<elem>_q-N.wfc` for `q < 0`

   where `<elem>` is the lower-case element symbol. The file format is
   the standard critic2 `.wfc` (LD1/QE Slater-orbital export). This is
   the rigorous route — drop in publication-quality charged densities
   when you have them.

2. **Built-in cation table.** For `q > 0` the existing `read_db` path
   produces a cation by zeroing trailing orbital occupations.

3. **Anion extrapolation.** For `q < 0`, `read_critic` has been
   extended: when the requested electron count exceeds the neutral
   orbital occupations, electrons are added to the highest-occupied
   orbital up to its angular-momentum capacity `2(2L+1)`. If the
   capacity is insufficient (which only happens for closed-shell
   anions beyond the orbitals stored in the standard `.wfc`), the
   routine emits a warning and uses the best-effort density it
   managed to build.

4. **Fall back to neutral** if all of the above fail, with a warning.

Iteration `Q_A` is clamped to the half-open interval `[-5, Z)`. The
upper limit `Q = Z` corresponds to a fully stripped atom (zero density
everywhere), which is correctly handled by the radial interpolator
returning zero.

## Worked example — water

Single-point B3LYP/6-31G(d) on H₂O, density cube generated with
`cubegen 4 fdensity=scf h2o.fchk h2o_rho_full.cube 160 h`.

Critic2 input:

```
molecule h2o_rho_full.cube
load    h2o_rho_full.cube
integrable 1
hirshfeld
hirshfeld_i
```

SCF trace:

```
# iter  max|dQ|       sum(Q)        atom-charges
  1     3.318E-01    0.00524    -0.3318  0.1685  0.1685
  2     1.735E-01    0.00524    -0.5053  0.2552  0.2552
  …
  16    1.045E-04    0.00524    -0.7292  0.3672  0.3672
  17    6.257E-05    0.00524    -0.7292  0.3672  0.3672
+ Hirshfeld-I SCF converged in 17 iterations
```

Results:

| Method        | Q(O)     | Q(H)     |
| ------------- | -------- | -------- |
| Hirshfeld     | −0.378   | +0.192   |
| Hirshfeld-I   | **−0.729** | **+0.367** |
| Literature¹   | −0.73 … −0.78 | +0.37 … +0.39 |

¹ Bultinck *et al.* 2007 and successors at comparable functionals/basis.
The Hirshfeld-I numbers reproduce the expected ~2× polarity
amplification over plain Hirshfeld.

## What changed in the source tree

* `src/integration.f90` — new `imtype_hirshfeld_i = 7` constant.
* `src/integration@proc.f90` — parses `HIRSHFELD_I` and the new
  `WFCDIR / HITOL / HIMAXIT` sub-options, dispatches to the SCF
  driver, and makes every code path that previously branched on
  `imtype_hirshfeld` also accept the iterative variant. Both the
  per-atom field integration and the bond-order overlap path now use
  the iterative density evaluator when active.
* `src/global@proc.F90` — the top-level command parser now
  recognises `hirshfeld_i` alongside `hirshfeld`.
* `src/hirshfeld.f90` / `src/hirshfeld@proc.f90` — new public
  procedures `hirsh_i_driver`, `hirsh_i_eval`, `hirsh_i_active`,
  `hirsh_i_cleanup`; module-level cache for `(Z, q)` radial
  densities; serial pre-load step that keeps the parallel grid
  sweep race-free.
* `src/grid1mod.f90` / `src/grid1mod@proc.f90` —
  - `read_critic` now augments occupations for anions
    (electron count > neutral occ): it fills the outermost orbital
    up to its `2(2L+1)` capacity.
  - `read_db` no longer rejects `q < 0`.
  - New public wrapper `grid1_read_file(g, file, z, q)` so a caller
    can load a `.wfc` from an absolute path (used by `WFCDIR`).
* `src/types.f90` — `basindat` gained the input fields `hi_wfcdir`,
  `hi_tol`, `hi_maxit` and the per-integration SCF state
  `hi_isactive`, `hi_qlo`, `hi_qhi`, `hi_frac`, `hi_qfinal`. The state
  lives on `bas` (mirroring how YT keeps `bas%luw`); only the shared
  `(Z, q) → grid1` memoisation cache is module-level in
  `hirshfeld@proc.f90`, where it can stay warm across multiple
  HIRSHFELD_I invocations like the neutral `agrid` table.
* `dat/helpdoc/hirshfeld_i.keyw` — keyword syntax shown by the GUI
  help.

## Limitations and follow-ups

* The GUI keyword catalog (`src/gui/templates@proc.f90`) is not
  updated in this change; the keyword is recognised by the parser but
  does not yet appear in the GUI keyword tree. The GUI build is
  off-by-default in critic2; this can be added later.
* Anion extrapolation only adds electrons to orbitals already present
  in the `.wfc` file. Anions that would need a new shell (e.g.
  closed-shell second-row anions like O²⁻ in some contexts) cannot be
  represented this way and trigger a warning. The `WFCDIR` route is
  the supported workaround.
* The `wcube` weight-cube output for `HIRSHFELD_I` currently calls the
  same `hirsh_weights` routine as plain Hirshfeld (neutral
  references). The integration itself uses the iterative weights; the
  visualisation cubes do not. A future change can have `hirsh_weights`
  consult `hirsh_i_eval` in active mode.
