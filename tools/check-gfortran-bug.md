# Checking a Fortran compiler for the polymorphic-dispatch codegen bug

Some gfortran versions miscompile a pattern the critic2 GUI relies on (interactive
graphics and dynamics): **passing an array section of a _polymorphic_ object's
allocatable component as the argument of a _type-bound_ (dispatch) call**, e.g.

```fortran
xf = c%c2x(md%r(:,i))   ! c is class(crystal); md is class(mdrun)
```

On an affected compiler the array descriptor for the component (`md%r`) is
emitted incorrectly, so the access faults. In the GUI this shows up as a crash
when loading any structure (atoms/bonds rendering) or when starting the
water-cluster / interactive dynamics.

## How to test a compiler

`repro.f90` (this directory) reproduces the bug in isolation — it is the exact
shape of `src/dynamics@proc.f90:157`. Build and run it with the compiler you
intend to use:

```sh
# native (MSYS2 / Linux)
gfortran -g -O0 -fcheck=all repro.f90 -o repro && ./repro

# cross-compiling from Linux to Windows
x86_64-w64-mingw32-gfortran -g -O0 -fcheck=all repro.f90 -o repro.exe
wine ./repro.exe          # or copy repro.exe to Windows and run it
```

- **AFFECTED** — the program aborts, with either
  `Fortran runtime error: Index '1' of dimension 2 of array 'md%r' outside of
  expected range (0:0)` or a raw `SIGSEGV`.
- **GOOD** — the program prints `COMPILER OK`.

This is a *code-generation* bug: it reproduces at `-O0`, so it depends only on
the compiler version and the target, **not** on native-vs-cross. A cross-compiler
built on an affected gfortran version behaves exactly like the native one.

## Known results

| gfortran version                              | target                | result       |
|-----------------------------------------------|-----------------------|--------------|
| **16.1.0** (MSYS2 UCRT; GCC-16 dev snapshot)  | x86_64 Windows (UCRT) | **AFFECTED** |
| **15.3.0** (WinLibs UCRT; latest GCC release) | x86_64 Windows (UCRT) | OK           |

So the bug is a **GCC-16 (development) regression**; the latest stable release
(15.3.0) is unaffected. GCC 16 is not a released series — MSYS2 ships trunk
snapshots. Prefer a stable GCC (<= 15.x) to build the Windows GUI.

If you must use an affected compiler, the source-level workaround is to copy the
section into a plain local variable before the dispatch call
(`ri = md%r(:,i); xf = c%c2x(ri)`), or make the object non-polymorphic. Run this
check whenever you change compilers (native upgrade, a new MSYS2 snapshot, or a
cross-compiler on a different distro).
