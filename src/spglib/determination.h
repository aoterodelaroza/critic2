// Copyright (C) 2017 Atsushi Togo
// This file is part of spglib.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef __determination_H__
#define __determination_H__

#include "cell.h"
#include "primitive.h"
#include "refinement.h"
#include "spacegroup.h"

typedef struct {
    Primitive *primitive;
    Spacegroup *spacegroup;
    ExactStructure *exact_structure;
} DataContainer;

DataContainer *det_determine_all(Cell const *cell, int const hall_number,
                                 double const symprec,
                                 double const angle_symprec);
void det_free_container(DataContainer *container);

#endif
