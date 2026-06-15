// Copyright (C) 2008 Atsushi Togo
// This file is part of spglib.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef __symmetry_H__
#define __symmetry_H__

#include "cell.h"
#include "mathfunc.h"

typedef struct {
    int size;
    int (*rot)[3][3];
    double (*trans)[3];
} Symmetry;

typedef struct {
    int size;
    int (*rot)[3][3];
    double (*trans)[3];
    int *timerev;
} MagneticSymmetry;

typedef struct {
    int rot[48][3][3];
    int size;
} PointSymmetry;

Symmetry *sym_alloc_symmetry(int const size);
SPG_API_TEST void sym_free_symmetry(Symmetry *symmetry);
MagneticSymmetry *sym_alloc_magnetic_symmetry(int const size);
void sym_free_magnetic_symmetry(MagneticSymmetry *symmetry);
SPG_API_TEST Symmetry *sym_get_operation(Cell const *primitive,
                                         double const symprec,
                                         double const angle_tolerance);
Symmetry *sym_reduce_operation(Cell const *primitive, Symmetry const *symmetry,
                               double const symprec,
                               double const angle_tolerance);
VecDBL *sym_get_pure_translation(Cell const *cell, double const symprec);
VecDBL *sym_reduce_pure_translation(Cell const *cell, VecDBL const *pure_trans,
                                    double const symprec,
                                    double const angle_tolerance);

#endif
