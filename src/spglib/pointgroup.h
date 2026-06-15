// Copyright (C) 2008 Atsushi Togo
// This file is part of spglib.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef __pointgroup_H__
#define __pointgroup_H__

#include "mathfunc.h"
#include "symmetry.h"

typedef enum {
    HOLOHEDRY_NONE,
    TRICLI,
    MONOCLI,
    ORTHO,
    TETRA,
    TRIGO,
    HEXA,
    CUBIC,
} Holohedry;

typedef enum {
    LAUE_NONE,
    LAUE1,
    LAUE2M,
    LAUEMMM,
    LAUE4M,
    LAUE4MMM,
    LAUE3,
    LAUE3M,
    LAUE6M,
    LAUE6MMM,
    LAUEM3,
    LAUEM3M,
} Laue;

typedef struct {
    int number;
    char symbol[6];
    char schoenflies[4];
    Holohedry holohedry;
    Laue laue;
} Pointgroup;

// @brief Return pointgroup.number = 0 if failed.
// @param[out] transform_mat
// @param[in] rotations
// @param[in] num_rotations
// @param[in] aperiodic_axis Use `aperiodic_axis=-1` for space group.
Pointgroup ptg_get_transformation_matrix(int transform_mat[3][3],
                                         int const rotations[][3][3],
                                         int const num_rotations,
                                         int const aperiodic_axis);
Pointgroup ptg_get_pointgroup(int const pointgroup_number);
PointSymmetry ptg_get_pointsymmetry(int const rotations[][3][3],
                                    int const num_rotations);
#endif
