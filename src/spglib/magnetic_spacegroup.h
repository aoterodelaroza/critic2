// Copyright (C) 2012 Atsushi Togo
// This file is part of spglib.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef __msg_H__
#define __msg_H__

#include "msg_database.h"
#include "spacegroup.h"
#include "symmetry.h"

typedef struct {
    /* Magnetic space-group type */
    int uni_number;
    int msg_type;
    int hall_number;
    /* Transformation to standardized setting */
    double transformation_matrix[3][3];
    double origin_shift[3];
    /* Rigid rotation to standardized lattice */
    double std_rotation_matrix[3][3];
} MagneticDataset;

MagneticDataset *msg_identify_magnetic_space_group_type(
    double const lattice[3][3], MagneticSymmetry const *magnetic_symmetry,
    double const symprec);
Cell *msg_get_transformed_cell(Cell const *cell, double const tmat[3][3],
                               double const origin_shift[3],
                               double const rigid_rot[3][3],
                               MagneticSymmetry const *magnetic_symmetry,
                               double const symprec,
                               double const angle_tolerance);

#endif /*__msg_H__ */
