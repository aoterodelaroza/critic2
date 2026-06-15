// Copyright (C) 2011 Atsushi Togo
// This file is part of spglib.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef __sitesym_database_H__
#define __sitesym_database_H__

int ssmdb_get_coordinate(int rot[3][3], double trans[3], int const index);
void ssmdb_get_wyckoff_indices(int indices[2], int const index);
void ssmdb_get_site_symmetry_symbol(char symbol[7], int const index);

#endif
