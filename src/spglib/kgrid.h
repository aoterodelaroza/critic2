// Copyright (C) 2015 Atsushi Togo
// This file is part of spglib.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef __kgrid_H__
#define __kgrid_H__

#include <stddef.h>

/* #define GRID_ORDER_XYZ */
/* This changes behaviour of index order of address. */
/* Without GRID_ORDER_XYZ, left most element of address runs first. */
/* grid_address (e.g. 4x4x4 mesh, unless GRID_ORDER_XYZ is defined) */
/*    [[ 0  0  0]                                                   */
/*     [ 1  0  0]                                                   */
/*     [ 2  0  0]                                                   */
/*     [-1  0  0]                                                   */
/*     [ 0  1  0]                                                   */
/*     [ 1  1  0]                                                   */
/*     [ 2  1  0]                                                   */
/*     [-1  1  0]                                                   */
/*     ....      ]                                                  */
/*                                                                  */
/* With GRID_ORDER_XYZ, right most element of address runs first.   */
/* grid_address (e.g. 4x4x4 mesh, if GRID_ORDER_XYZ is defined)     */
/*    [[ 0  0  0]                                                   */
/*     [ 0  0  1]                                                   */
/*     [ 0  0  2]                                                   */
/*     [ 0  0 -1]                                                   */
/*     [ 0  1  0]                                                   */
/*     [ 0  1  1]                                                   */
/*     [ 0  1  2]                                                   */
/*     [ 0  1 -1]                                                   */
/*     ....      ]                                                  */

/* #define GRID_BOUNDARY_AS_NEGATIVE */
/* This changes the behaviour of address elements on the surface of  */
/* parallelepiped. */
/* For odd mesh number, this affects nothing, e.g., [-2, -1, 0, 1, 2]. */
/* regardless of with and without GRID_BOUNDARY_AS_NEGATIVE. */
/* For even mesh number, this affects as follows: */
/* without GRID_BOUNDARY_AS_NEGATIVE, e.g., [-2, -1, 0, 1, 2, 3]. */
/* with GRID_BOUNDARY_AS_NEGATIVE, e.g., [-3, -2, -1, 0, 1, 2]. */

void kgd_get_all_grid_addresses(int grid_address[][3], int const mesh[3]);
int kgd_get_grid_point_double_mesh(int const address_double[3],
                                   int const mesh[3]);
size_t kgd_get_dense_grid_point_double_mesh(int const address_double[3],
                                            int const mesh[3]);
void kgd_get_grid_address_double_mesh(int address_double[3],
                                      int const address[3], int const mesh[3],
                                      int const is_shift[3]);

#endif
