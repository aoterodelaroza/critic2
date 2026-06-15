// Copyright (C) 2017 Atsushi Togo
// This file is part of spglib.
// SPDX-License-Identifier: BSD-3-Clause

#include "cell.h"
#include "mathfunc.h"

/* Contains pre-allocated memory and precomputed data for check_total_overlap.
 */
typedef struct {
    /* Number of atoms. */
    int size;

    /* Pre-allocated memory for various things. */
    void *argsort_work;
    void *blob;

    /* Temp areas for writing stuff. (points into blob) */
    double (*pos_temp_1)[3];
    double (*pos_temp_2)[3];

    /* Temp area for writing lattice point distances. (points into blob) */
    double *distance_temp; /* for lattice point distances */
    int *perm_temp;        /* for permutations during sort */

    /* Sorted data of original cell. (points into blob)*/
    double (*lattice)[3];
    double (*pos_sorted)[3];
    int *types_sorted;

    /* Using array reference to avoid redundant loop */
    int *periodic_axes;
} OverlapChecker;

OverlapChecker *ovl_overlap_checker_init(Cell const *cell);

int ovl_check_total_overlap(OverlapChecker *checker, double const test_trans[3],
                            int const rot[3][3], double const symprec,
                            int const is_identity);

int ovl_check_layer_total_overlap(OverlapChecker *checker,
                                  double const test_trans[3],
                                  int const rot[3][3], double const symprec,
                                  int const is_identity);

void ovl_overlap_checker_free(OverlapChecker *checker);
