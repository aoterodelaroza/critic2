// Copyright (C) 2010 Atsushi Togo
// This file is part of spglib.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef __hall_symbol_H__
#define __hall_symbol_H__

#include "mathfunc.h"
#include "spacegroup.h"
#include "symmetry.h"

int hal_match_hall_symbol_db(double origin_shift[3],
                             double const bravais_lattice[3][3],
                             int const hall_number, Centering const centering,
                             Symmetry const *symmetry, double const symprec);

#endif
