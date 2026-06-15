// Copyright (C) 2008 Atsushi Togo
// This file is part of spglib.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef __cell_H__
#define __cell_H__

#include "mathfunc.h"

typedef enum {
    NOSPIN = -1,
    COLLINEAR = 0,
    NONCOLLINEAR = 1,
} SiteTensorType;

typedef struct {
    /* Number of atoms */
    int size;
    /* Used for layer group. Set -1 for space-group search */
    int aperiodic_axis;
    /* 3x3 matrix */
    double (*lattice)[3];
    /* Atomic types with length (size, ) */
    int *types;
    /* Scaled positions with length (size, 3) */
    double (*position)[3];
    /* Rank of site tensors. Set COLLINEAR for scalar, */
    /* and NONCOLLINEAR for vector. */
    /* If no site tensors, set SiteTensorType.NOSPIN. */
    SiteTensorType tensor_rank;
    /* For tensor_rank=COLLINEAR, site tensors with (size, ).*/
    /* For tensor_rank=NONCOLLINEAR, site tensors with (size * 3, ).*/
    double *tensors;
} Cell;

SPG_API_TEST Cell *cel_alloc_cell(int const size,
                                  SiteTensorType const tensor_rank);
SPG_API_TEST void cel_free_cell(Cell *cell);
void cel_set_cell(Cell *cell, double const lattice[3][3],
                  double const position[][3], int const types[]);
SPG_API_TEST void cel_set_layer_cell(Cell *cell, double const lattice[3][3],
                                     double const position[][3],
                                     int const types[],
                                     int const aperiodic_axis);
void cel_set_cell_with_tensors(Cell *cell, double const lattice[3][3],
                               double const position[][3], int const types[],
                               double const *tensors);
Cell *cel_copy_cell(Cell const *cell);
int cel_is_overlap(double const a[3], double const b[3],
                   double const lattice[3][3], double const symprec);
int cel_is_overlap_with_same_type(double const a[3], double const b[3],
                                  int const type_a, int const type_b,
                                  double const lattice[3][3],
                                  double const symprec);
int cel_any_overlap(Cell const *cell, double const symprec);
int cel_any_overlap_with_same_type(Cell const *cell, double const symprec);
Cell *cel_trim_cell(int *mapping_table, double const trimmed_lattice[3][3],
                    Cell const *cell, double const symprec);
int cel_layer_is_overlap(double const a[3], double const b[3],
                         double const lattice[3][3], int const periodic_axes[2],
                         double const symprec);
int cel_layer_is_overlap_with_same_type(double const a[3], double const b[3],
                                        int const type_a, int const type_b,
                                        double const lattice[3][3],
                                        int const periodic_axes[2],
                                        double const symprec);
int cel_layer_any_overlap_with_same_type(Cell const *cell,
                                         int const periodic_axes[2],
                                         double const symprec);

#endif
