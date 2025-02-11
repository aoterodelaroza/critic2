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

// Return the number of atoms from crystal cr.
int c2_crystal_get_natom(crystal *cr);

// Get the number of atoms (nat), cell lengths (cel, bohr), cell
// angles (ang, degrees), lattice vectors (lattice, bohr),
// atomic positions (position, fractional coords), and atomic numbers
// (zat) from the crystal structure cr. The input variables must point
// to enough space to contain the output (cel[3], ang[3], lattice[3][3],
// position[natom][3], zat[natom]). If the pointer is NULL, do not
// return that variable.
void c2_crystal_get_structure(crystal *cr,int *natom,double *cel,double *ang,
                              double *lattice,double *position, int *zat);

// Write the crystal structure to a file. The format of the file is
// detected from the extension.
crystal *c2_write_crystal(crystal *cr,const char *file);

// Destroy the input crystal structure object and free the memory.
void c2_destroy_crystal(crystal *cr);

//// PEAKS ////

// Calculate the XRPD peak positions from the crystal structure.
// Limit the peaks to the 2theta range th2ini (default: 5) to th2end
// (default: 50. lambda is the wavelength in angstrom (default:
// 1.5406). fpol is the polarization correction factor (0 =
// unpolarized, 0.95 = synchrotron, default: 0). If th2ini, th2end,
// lambda, fpol < 0, use default values. Returns NULL on error.
xrpd_peaklist *c2_peaks_from_crystal(crystal *cr,double th2ini,double th2end,double lambda,
                                     double fpol);

// Read xrpd peaks from a peaks file. Allocate space for the crystal
// structure and return a pointer to it or NULL if there was an error.
xrpd_peaklist *c2_peaks_from_file(const char *file);

// Fit XRPD profile data from an xy file (xyfile) and return the peak
// list and the RMS of the fit (rms), or NULL if there was an
// error. Other arguments:
// - verbose: print fitting info to the stdout.
// - ymax_detect: maxima are added to the initial model as candidate
// peaks if their intensities are higher than this value.
// - def_ymax_detect: use the default ymax_detect (the median of the
// profile intensities). If this is true, ymax_detect is ignored.
// - nadj: a candidate peak is added if the maximum is surrounded by
//   nadj points on either side with smaller intensity. If negative,
//   use the default value (2).
// - pkinput: if non-NULL use a peak list as the initial model (nadj0
//   and ymax_detect0 are not used if pkinput is given).
xrpd_peaklist *c2_peaks_from_profile(const char *xyfile,double *rms,bool verbose,
				     bool def_ymax_detect,double ymax_detect,
				     int nadj,xrpd_peaklist *pkinput);

// Write the peak list to a file.
void c2_write_peaks(xrpd_peaklist *pk, const char *file);

// Calculate the diffraction profile from peak list pk with n points
// between 2*theta values th2ini and th2end (if th2ini and th2end are
// negative, use the th2ini and th2end in the peak list
// structure). Return the 2*theta in x[] and the profile intensity in
// y[]. x and y must be freed after use. If error, return NULL in x
// and y.
void c2_peaks_calculate_profile(xrpd_peaklist *pk,int n,double th2ini,double th2end,
				double **x,double **y);

// Destroy the input XRPD peaks structure object and free the memory.
void c2_destroy_peaks(xrpd_peaklist *pk);

//// Background estimation ////
// Read a XRPD pattern from file xyfile and estimate the background
// contribution using David/Sivia's method with nknot knots in the
// cubic spline (if nknot <= 0, defaults to 20). Write the resulting
// background to file xyback (if non-NULL) and the original pattern
// minus the background to file xyclean (if non-NULL).
void c2_profile_background(const char *xyfile,const char *xyback,const char *xyclean,
			   int nknot);

//// (VC-)GPWDF ////

// Compare crystal c1 and set of XRPD peaks p2 using GPWDF. alpha =
// Gaussian triangle width. lambda = wavelength in angstrom. fpol =
// polarization correction factor (0 = unpolarized, 0.95 =
// synchrotron). If alpha, lambda, fpol < 0, use default values.
// Returns -1 on error.
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
// Returns -1 on error.
double c2_compare_vcgpwdf(crystal *c1,xrpd_peaklist *p2,crystal **crout,bool global,bool verbose,
                          double alpha,double lambda,double fpol,int maxfeval,
                          double besteps,double max_elong,double max_ang);

// Compare crystal c1 and set of XRPD peaks p2 using variable-cell
// (VC-)GPWDF, global minimization. Use "safe" default settings. If
// verbose, print search progress to stdout. lambda = wavelength in
// angstrom. fpol = polarization correction factor (0 = unpolarized,
// 0.95 = synchrotron). Returns the VC-GPDWF score and the deformed
// c1 structure in crout. Returns -1 on error.
double c2_compare_vcgpwdf_global_safe(crystal *c1,xrpd_peaklist *p2,crystal **crout,bool verbose,
                                      double lambda,double fpol);

// Compare crystal c1 and set of XRPD peaks p2 using variable-cell
// (VC-)GPWDF, global minimization. Use "quick" default settings. If
// verbose, print search progress to stdout. lambda = wavelength in
// angstrom. fpol = polarization correction factor (0 = unpolarized,
// 0.95 = synchrotron). Returns the VC-GPDWF score and the deformed
// c1 structure in crout. Returns -1 on error.
double c2_compare_vcgpwdf_global_quick(crystal *c1,xrpd_peaklist *p2,crystal **crout,bool verbose,
                                       double lambda,double fpol);

#ifdef __cplusplus
}
#endif
#endif
