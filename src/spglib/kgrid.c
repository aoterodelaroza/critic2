// Copyright (C) 2015 Atsushi Togo
// This file is part of spglib.
// SPDX-License-Identifier: BSD-3-Clause

#include "kgrid.h"

#include <assert.h>
#include <stddef.h>

static void get_all_grid_addresses(int grid_address[][3], int const mesh[3]);
static size_t get_grid_point_double_mesh(int const address_double[3],
                                         int const mesh[3]);
static size_t get_grid_point_single_mesh(int const address[3],
                                         int const mesh[3]);
static void modulo_i3(int v[3], int const m[3]);
static void reduce_grid_address(int address[3], int const mesh[3]);
static void reduce_grid_address_double(int address[3], int const mesh[3]);

void kgd_get_all_grid_addresses(int grid_address[][3], int const mesh[3]) {
    get_all_grid_addresses(grid_address, mesh);
}

int kgd_get_grid_point_double_mesh(int const address_double[3],
                                   int const mesh[3]) {
    return get_grid_point_double_mesh(address_double, mesh);
}

size_t kgd_get_dense_grid_point_double_mesh(int const address_double[3],
                                            int const mesh[3]) {
    return get_grid_point_double_mesh(address_double, mesh);
}

void kgd_get_grid_address_double_mesh(int address_double[3],
                                      int const address[3], int const mesh[3],
                                      int const is_shift[3]) {
    int i;

    for (i = 0; i < 3; i++) {
        address_double[i] = address[i] * 2 + (is_shift[i] != 0);
    }
    reduce_grid_address_double(address_double, mesh);
}

static void get_all_grid_addresses(int grid_address[][3], int const mesh[3]) {
    int i, j, k;
    size_t grid_point;
    int address[3];

    for (i = 0; i < mesh[0]; i++) {
        address[0] = i;
        for (j = 0; j < mesh[1]; j++) {
            address[1] = j;
            for (k = 0; k < mesh[2]; k++) {
                address[2] = k;
                grid_point = get_grid_point_single_mesh(address, mesh);

                assert((size_t)(mesh[0] * mesh[1] * mesh[2]) > grid_point);

                grid_address[grid_point][0] = address[0];
                grid_address[grid_point][1] = address[1];
                grid_address[grid_point][2] = address[2];
                reduce_grid_address(grid_address[grid_point], mesh);
            }
        }
    }
}

static size_t get_grid_point_double_mesh(int const address_double[3],
                                         int const mesh[3]) {
    int i;
    int address[3];

    for (i = 0; i < 3; i++) {
        if (address_double[i] % 2 == 0) {
            address[i] = address_double[i] / 2;
        } else {
            address[i] = (address_double[i] - 1) / 2;
        }
    }
    modulo_i3(address, mesh);

    return get_grid_point_single_mesh(address, mesh);
}

static size_t get_grid_point_single_mesh(int const address[3],
                                         int const mesh[3]) {
#ifndef GRID_ORDER_XYZ
    return (address[2] * mesh[0] * (size_t)(mesh[1]) + address[1] * mesh[0] +
            address[0]);
#else
    return (address[0] * mesh[1] * (size_t)(mesh[2]) + address[1] * mesh[2] +
            address[2]);
#endif
}

static void modulo_i3(int v[3], int const m[3]) {
    int i;

    for (i = 0; i < 3; i++) {
        v[i] = v[i] % m[i];

        if (v[i] < 0) {
            v[i] += m[i];
        }
    }
}

static void reduce_grid_address(int address[3], int const mesh[3]) {
    int i;

    for (i = 0; i < 3; i++) {
#ifndef GRID_BOUNDARY_AS_NEGATIVE
        address[i] -= mesh[i] * (address[i] > mesh[i] / 2);
#else
        address[i] -= mesh[i] * (address[i] > (mesh[i] - 1) / 2);
#endif
    }
}

static void reduce_grid_address_double(int address[3], int const mesh[3]) {
    int i;

    for (i = 0; i < 3; i++) {
#ifndef GRID_BOUNDARY_AS_NEGATIVE
        address[i] -= 2 * mesh[i] * (address[i] > mesh[i]);
#else
        address[i] -= 2 * mesh[i] * (address[i] > mesh[i] - 1);
#endif
    }
}
