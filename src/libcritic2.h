/*
Copyright (c) 2015-2022 Alberto Otero de la Roza <aoterodelaroza@gmail.com>,
Ángel Martín Pendás <angel@fluor.quimica.uniovi.es> and Víctor Luaña
<victor@fluor.quimica.uniovi.es>.

critic2 is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at
your option) any later version.

critic2 is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef LIBCRITIC2_H
#define LIBCRITIC2_H

#ifdef __cplusplus
extern "C" {
#endif

//// Types ////
typedef void crystal;
typedef void xrpd_peaklist;

//// Functions ////

//// CRYSTAL ////

// Create a crystal structure object from a file. The format of the
// file is detected from the extension. Allocate space for the crystal
// structure and return a pointer to it or NULL if there was an
// error. Requires destruction of the object after its use.
crystal *c2_crystal_from_file(const char *file);

// Create a crystal structure from the lattice parameters (lattice,
// in bohr), number of atoms (natom), atomic positions (position,
// fractional coords), and atomic numbers (zat). Allocate space for
// the crystal structure and return a pointer to it or NULL if there
// was an error. Requires destruction of the object after its use.
crystal *c2_crystal_from_lattice(const int natom,const double lattice[3][3],
                                 const double position[][3],const int zat[]);

// Create a crystal structure from the cell lengths (cel, in bohr),
// cell angles (ang, degrees), number of atoms (natom), atomic
// positions (position, fractional coords), and atomic numbers
// (zat). Allocate space for the crystal structure and return a
// pointer to it or NULL if there was an error. Requires destruction
// of the object after its use.
crystal *c2_crystal_from_cellpar(const int natom,const double cel[3], const double ang[3],
                                 const double position[][3],const int zat[]);

// Write the report about the input crystal structure to standard output.
void c2_describe_crystal(crystal *cr);

// Write the crystal structure to a file. The format of the file is
// detected from the extension.
crystal *c2_write_crystal(crystal *cr,const char *file);

// Destroy the input crystal structure object and free the memory.
void c2_destroy_crystal(crystal *cr);

//// PEAKS ////

// Calculate the XRPD peak positions from the crystal structure.
// Limit the peaks to the 2theta range th2ini to th2end. lambda is
// the wavelength in angstrom. fpol is the polarization correction
// factor (0 = unpolarized, 0.95 = synchrotron). If th2ini, th2end,
// lambda, fpol < 0, use default values.
xrpd_peaklist *c2_peaks_from_crystal(crystal *cr,double th2ini,double th2end,double lambda,
                                     double fpol);

// Read xrpd peaks from a file. Allocate space for the crystal
// structure and return a pointer to it or NULL if there was an
// error.
xrpd_peaklist *c2_peaks_from_file(const char *file);

// Destroy the input XRPD peaks structure object and free the memory.
void c2_destroy_peaks(xrpd_peaklist *pk);

//// (VC-)GPWDF ////

// Compare crystal c1 and set of XRPD peaks p2 using GPWDF. alpha =
// Gaussian triangle width. lambda = wavelength in angstrom. fpol =
// polarization correction factor (0 = unpolarized, 0.95 =
// synchrotron). If alpha, lambda, fpol < 0, use default values.
double c2_compare_gpwdf(crystal *c1,xrpd_peaklist *p2,double alpha,double lambda,double fpol);

// Compare crystal c1 and set of XRPD peaks p2 using variable-cell
// (VC-)GPWDF. If global, use the global minimization for the best
// GPWDF; otherwise use only a local minimization. If verbose, print
// search progress to stdout. alpha = Gaussian triangle
// width. lambda = wavelength in angstrom. fpol = polarization
// correction factor (0 = unpolarized, 0.95 = synchrotron). maxfeval
// = maximum number of function evaluations. besteps = do not
// restart the feval count if a diff lower than the best diff is
// found within besteps. max_elong = maximum cell length elongation
// (%). max_ang = maximum cell angle deformation (degrees). Returns
// the VC-GPDWF score and the deformed c1 structure in crout.
double c2_compare_vcgpwdf(crystal *c1,xrpd_peaklist *p2,crystal **crout,bool global,bool verbose,
                          double alpha,double lambda,double fpol,int maxfeval,
                          double besteps,double max_elong,double max_ang);

// Compare crystal c1 and set of XRPD peaks p2 using variable-cell
// (VC-)GPWDF, global minimization. Use "safe" default settings. If
// verbose, print search progress to stdout. lambda = wavelength in
// angstrom. fpol = polarization correction factor (0 = unpolarized,
// 0.95 = synchrotron). Returns the VC-GPDWF score and the deformed
// c1 structure in crout.
double c2_compare_vcgpwdf_global_safe(crystal *c1,xrpd_peaklist *p2,crystal **crout,bool verbose,
                                      double lambda,double fpol);

// Compare crystal c1 and set of XRPD peaks p2 using variable-cell
// (VC-)GPWDF, global minimization. Use "quick" default settings. If
// verbose, print search progress to stdout. lambda = wavelength in
// angstrom. fpol = polarization correction factor (0 = unpolarized,
// 0.95 = synchrotron). Returns the VC-GPDWF score and the deformed
// c1 structure in crout.
double c2_compare_vcgpwdf_global_quick(crystal *c1,xrpd_peaklist *p2,crystal **crout,bool verbose,
                                       double lambda,double fpol);

#ifdef __cplusplus
}
#endif
#endif
