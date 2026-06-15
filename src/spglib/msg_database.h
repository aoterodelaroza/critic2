// Copyright (C) 2010 Atsushi Togo
// This file is part of spglib.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef __msg_database_H__
#define __msg_database_H__

#include "symmetry.h"

typedef struct {
    int uni_number;
    int litvin_number;
    char bns_number[8];
    char og_number[12];
    int number;
    int type;
} MagneticSpacegroupType;

MagneticSpacegroupType msgdb_get_magnetic_spacegroup_type(int const uni_number);
MagneticSymmetry *msgdb_get_spacegroup_operations(int const uni_number,
                                                  int const hall_number);
void msgdb_get_uni_candidates(int uni_number_range[2], int const hall_number);
Symmetry *msgdb_get_std_transformations(int const uni_number,
                                        int const hall_number);

#endif /* __msg_database_H__ */
