// Copyright (C) 2012 Atsushi Togo
// This file is part of spglib.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef __spin_H__
#define __spin_H__

#include "cell.h"
#include "mathfunc.h"
#include "symmetry.h"

/**
 * @brief
 *
 * @param[out] equivalent_atoms
 * @param[out] permutations such that the p-th operation in `magnetic_symmetry`
 * maps site-`i` to site-`permutations[p * cell->size + i]`.
 * @param[out] prim_lattice
 * @param[in] sym_nonspin Symmetry operations with ignoring spin
 * @param[in] cell
 * @param[in] with_time_reversal true if consider time reversal operation
 * @param[in] is_axial true if site tensors are axial w.r.t. time-reversal
 * operations
 * @param[in] symprec
 * @param[in] angle_tolerance
 * @param[in] mag_symprec if mag_sympprec < 0, use symprec instead
 * @return Return NULL if failed.
 */
MagneticSymmetry *spn_get_operations_with_site_tensors(
    int **equivalent_atoms, int **permutations, double prim_lattice[3][3],
    Symmetry const *sym_nonspin, Cell const *cell, int const with_time_reversal,
    int const is_axial, double const symprec, double const angle_tolerance,
    double const mag_symprec);
VecDBL *spn_collect_pure_translations_from_magnetic_symmetry(
    MagneticSymmetry const *sym_msg);
Cell *spn_get_idealized_cell(int const *permutations, Cell const *cell,
                             MagneticSymmetry const *magnetic_symmetry,
                             int const with_time_reversal, int const is_axial);
double *spn_alloc_site_tensors(int const num_atoms, int const tensor_rank);

#endif
