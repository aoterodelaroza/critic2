// Copyright (C) 2008 Atsushi Togo
// This file is part of spglib.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef __kpoint_H__
#define __kpoint_H__

#include <stddef.h>

#include "mathfunc.h"

int kpt_get_irreducible_reciprocal_mesh(int grid_address[][3],
                                        int ir_mapping_table[],
                                        int const mesh[3],
                                        int const is_shift[3],
                                        MatINT const *rot_reciprocal);
size_t kpt_get_dense_irreducible_reciprocal_mesh(int grid_address[][3],
                                                 size_t ir_mapping_table[],
                                                 int const mesh[3],
                                                 int const is_shift[3],
                                                 MatINT const *rot_reciprocal);
int kpt_get_stabilized_reciprocal_mesh(
    int grid_address[][3], int ir_mapping_table[], int const mesh[3],
    int const is_shift[3], int const is_time_reversal, MatINT const *rotations,
    size_t const num_q, double const qpoints[][3]);
size_t kpt_get_dense_stabilized_reciprocal_mesh(
    int grid_address[][3], size_t ir_mapping_table[], int const mesh[3],
    int const is_shift[3], int const is_time_reversal, MatINT const *rotations,
    size_t const num_q, double const qpoints[][3]);
void kpt_get_dense_grid_points_by_rotations(size_t rot_grid_points[],
                                            int const address_orig[3],
                                            int const (*rot_reciprocal)[3][3],
                                            int const num_rot,
                                            int const mesh[3],
                                            int const is_shift[3]);
void kpt_get_dense_BZ_grid_points_by_rotations(
    size_t rot_grid_points[], int const address_orig[3],
    int const (*rot_reciprocal)[3][3], int const num_rot, int const mesh[3],
    int const is_shift[3], size_t const bz_map[]);
int kpt_relocate_BZ_grid_address(int bz_grid_address[][3], int bz_map[],
                                 int const grid_address[][3], int const mesh[3],
                                 double const rec_lattice[3][3],
                                 int const is_shift[3]);
size_t kpt_relocate_dense_BZ_grid_address(
    int bz_grid_address[][3], size_t bz_map[], int const grid_address[][3],
    int const mesh[3], double const rec_lattice[3][3], int const is_shift[3]);
MatINT *kpt_get_point_group_reciprocal(MatINT const *rotations,
                                       int const is_time_reversal);
MatINT *kpt_get_point_group_reciprocal_with_q(MatINT const *rot_reciprocal,
                                              double const symprec,
                                              size_t const num_q,
                                              double const qpoints[][3]);

#endif
