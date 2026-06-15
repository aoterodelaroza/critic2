// Copyright (C) 2010 Atsushi Togo
// This file is part of spglib.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef __spacegroup_H__
#define __spacegroup_H__

#include "cell.h"
#include "mathfunc.h"
#include "primitive.h"
#include "symmetry.h"

typedef struct {
    int number;
    int hall_number;
    int pointgroup_number;
    char schoenflies[7];
    char hall_symbol[17];
    char international[32];
    char international_long[20];
    char international_short[11];
    char choice[6];
    double bravais_lattice[3][3];
    double origin_shift[3];
} Spacegroup;

typedef enum {
    CENTERING_ERROR,
    PRIMITIVE,
    BODY,
    FACE,
    A_FACE,
    B_FACE,
    C_FACE,
    BASE,
    R_CENTER,
} Centering;

Spacegroup *spa_search_spacegroup(Primitive const *primitive,
                                  int const hall_number, double const symprec,
                                  double const angle_tolerance);
Spacegroup *spa_search_spacegroup_with_symmetry(Symmetry const *symmetry,
                                                double const prim_lat[3][3],
                                                double const symprec);
Cell *spa_transform_to_primitive(int *mapping_table, Cell const *cell,
                                 double const trans_mat[3][3],
                                 Centering const centering,
                                 double const symprec);
Cell *spa_transform_from_primitive(Cell const *primitive,
                                   Centering const centering,
                                   double const symprec);
void spa_copy_spacegroup(Spacegroup *dst, Spacegroup const *src);

#endif
