// Copyright (C) 2017 Atsushi Togo
// This file is part of spglib.
// SPDX-License-Identifier: BSD-3-Clause

#include "determination.h"

#include <stdlib.h>

#include "cell.h"
#include "debug.h"
#include "primitive.h"
#include "refinement.h"
#include "spacegroup.h"

#define REDUCE_RATE_OUTER 0.9
#define NUM_ATTEMPT_OUTER 10
#define REDUCE_RATE 0.95
#define ANGLE_REDUCE_RATE 0.95
#define NUM_ATTEMPT 20

static DataContainer *get_spacegroup_and_primitive(Cell const *cell,
                                                   int const hall_number,
                                                   double const symprec,
                                                   double const angle_symprec);

DataContainer *det_determine_all(Cell const *cell, int const hall_number,
                                 double const symprec,
                                 double const angle_symprec) {
    int attempt;
    double tolerance;
    DataContainer *container;

    container = NULL;

    if (hall_number > 530) {
        return NULL;
    }

    tolerance = symprec;
    for (attempt = 0; attempt < NUM_ATTEMPT_OUTER; attempt++) {
        if ((container = get_spacegroup_and_primitive(
                 cell, hall_number, tolerance, angle_symprec)) != NULL) {
            if ((container->exact_structure =
                     ref_get_exact_structure_and_symmetry(
                         container->spacegroup, container->primitive->cell,
                         cell, container->primitive->mapping_table,
                         container->primitive->tolerance)) != NULL) {
                goto found;
            }
            debug_print(
                "spglib: ref_get_exact_structure_and_symmetry failed.\n");
            det_free_container(container);
            container = NULL;
        }
        tolerance *= REDUCE_RATE_OUTER;
    }

found:
    return container;
}

void det_free_container(DataContainer *container) {
    if (container != NULL) {
        if (container->spacegroup != NULL) {
            free(container->spacegroup);
            container->spacegroup = NULL;
        }
        if (container->primitive != NULL) {
            prm_free_primitive(container->primitive);
            container->primitive = NULL;
        }
        if (container->exact_structure != NULL) {
            ref_free_exact_structure(container->exact_structure);
            container->exact_structure = NULL;
        }
        free(container);
    }
}

/* NULL is returned if failed */
static DataContainer *get_spacegroup_and_primitive(Cell const *cell,
                                                   int const hall_number,
                                                   double const symprec,
                                                   double const angle_symprec) {
    int attempt;
    double tolerance, angle_tolerance;
    DataContainer *container;

    debug_print("get_spacegroup_and_primitive (tolerance = %f):\n", symprec);

    container = NULL;

    if ((container = (DataContainer *)malloc(sizeof(DataContainer))) == NULL) {
        warning_memory("container");
        return NULL;
    }

    container->primitive = NULL;
    container->spacegroup = NULL;
    container->exact_structure = NULL;

    tolerance = symprec;
    angle_tolerance = angle_symprec;

    for (attempt = 0; attempt < NUM_ATTEMPT; attempt++) {
        if ((container->primitive =
                 prm_get_primitive(cell, tolerance, angle_tolerance)) != NULL) {
            debug_print("primitive lattice\n");
            debug_print_matrix_d3(container->primitive->cell->lattice);

            if ((container->spacegroup = spa_search_spacegroup(
                     container->primitive, hall_number,
                     container->primitive->tolerance,
                     container->primitive->angle_tolerance)) != NULL) {
                goto found;
            }

            prm_free_primitive(container->primitive);
            container->primitive = NULL;
        }

        debug_print("spglib: Attempt %d tolerance = %f failed.\n", attempt,
                    tolerance);

        tolerance *= REDUCE_RATE;
        if (angle_tolerance > 0) {
            angle_tolerance *= ANGLE_REDUCE_RATE;
        }
    }

    det_free_container(container);
    container = NULL;

    return NULL;

found:
    return container;
}
