// Copyright (C) 2010 Atsushi Togo
// This file is part of spglib.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef __delaunay_H__
#define __delaunay_H__

#include "base.h"
#include "mathfunc.h"

SPG_API_TEST int del_delaunay_reduce(double lattice_new[3][3],
                                     double const lattice[3][3],
                                     double const symprec);
SPG_API_TEST int del_layer_delaunay_reduce(double min_lattice[3][3],
                                           double const lattice[3][3],
                                           int const aperiodic_axis,
                                           double const symprec);

// @brief Delaunay reduction for monoclinic/oblique or monoclinic/rectangular
// @param[out] red_lattice
// @param[in] lattice
// @param[in] unique_axis
//            Two-fold axis or mirror-plane-perpendicular axis
// @param[in] aperiodic_axis
// @param[in] symprec
// @note For Monoclinic/oblique, the unique axis is also the aperiodic axis.
//       Axes are {j, k, unique_axis(=aperiodic_axis)}.
//       For Monoclinic/rectangular, axes are {unique_axis, j,
//       k(=aperiodic_axis)}. j and k are delaunay reduced, which can be
//       incomplete for Monoclinic/Rectangular
int del_layer_delaunay_reduce_2D(double min_lattice[3][3],
                                 double const lattice[3][3],
                                 int const unique_axis,
                                 int const aperiodic_axis,
                                 double const symprec);

#endif
