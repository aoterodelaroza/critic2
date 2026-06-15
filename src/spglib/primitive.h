// Copyright (C) 2008 Atsushi Togo
// This file is part of spglib.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef __primitive_H__
#define __primitive_H__

#include "cell.h"
#include "mathfunc.h"
#include "symmetry.h"

typedef struct {
    Cell *cell;
    int *mapping_table;
    int size;
    double tolerance;
    double angle_tolerance;
    double (*orig_lattice)[3]; /* 3x3 matrix */
} Primitive;

Primitive *prm_alloc_primitive(int const size);
void prm_free_primitive(Primitive *primitive);
Primitive *prm_get_primitive(Cell const *cell, double const symprec,
                             double const angle_tolerance);
int prm_get_primitive_with_pure_trans(Primitive *primitive, Cell const *cell,
                                      VecDBL const *pure_trans,
                                      double const symprec,
                                      double const angle_tolerance);
Symmetry *prm_get_primitive_symmetry(double t_mat[3][3],
                                     Symmetry const *symmetry,
                                     double const symprec);
int prm_get_primitive_lattice_vectors(double prim_lattice[3][3],
                                      Cell const *cell,
                                      VecDBL const *pure_trans,
                                      double const symprec,
                                      double const angle_tolerance);

#endif
