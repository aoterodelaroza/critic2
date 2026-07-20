/*
Copyright (c) 2007-2022 Alberto Otero de la Roza <aoterodelaroza@gmail.com>

critic2 is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at
your option) any later version.

critic2 is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

/* Smoke test for the libcritic2 C API (BUILD_LIBRARY=ON). This checks
   that the shared library links, that its symbols are exported, and
   that a crystal structure survives a build/query round trip. It uses
   no data from dat/, so it can carry the "nodata" label. */

#include <stdio.h>
#include <math.h>

#include "libcritic2.h"

static int nfail = 0;

static void check_double(const char *what, double got, double expect, double tol) {
  if (!(fabs(got - expect) <= tol)) {
    printf("FAIL %s: got %.10f, expected %.10f (tol %g)\n", what, got, expect, tol);
    nfail++;
  }
}

static void check_int(const char *what, int got, int expect) {
  if (got != expect) {
    printf("FAIL %s: got %d, expected %d\n", what, got, expect);
    nfail++;
  }
}

int main(void) {
  /* CsCl-like cubic cell: a = 5 bohr, one Na at the origin and one Cl
     at the body centre. Chosen so the structure is unambiguous and
     needs no external data. */
  const double cel[3] = {5.0, 5.0, 5.0};
  const double ang[3] = {90.0, 90.0, 90.0};
  const double position[2][3] = {{0.0, 0.0, 0.0}, {0.5, 0.5, 0.5}};
  const int zat[2] = {11, 17};

  crystal *cr = c2_crystal_from_cellpar(2, cel, ang, position, zat);
  if (!cr) {
    printf("FAIL c2_crystal_from_cellpar returned NULL\n");
    return 1;
  }

  check_int("c2_crystal_get_natom", c2_crystal_get_natom(cr), 2);

  /* read the structure back and compare with the input */
  int natom = 0;
  double ocel[3], oang[3], olattice[3][3], oposition[2][3];
  int ozat[2];
  c2_crystal_get_structure(cr, &natom, ocel, oang, (double *) olattice,
                           (double *) oposition, ozat);

  check_int("get_structure natom", natom, 2);
  for (int i = 0; i < 3; i++) {
    char what[64];
    sprintf(what, "cell length %d", i);
    check_double(what, ocel[i], cel[i], 1e-10);
    sprintf(what, "cell angle %d", i);
    check_double(what, oang[i], ang[i], 1e-10);
  }

  /* for this cubic cell the lattice vectors are a*identity */
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      char what[64];
      sprintf(what, "lattice[%d][%d]", i, j);
      check_double(what, olattice[i][j], (i == j) ? 5.0 : 0.0, 1e-10);
    }
  }

  /* the atoms may come back in either order, so match by atomic number */
  for (int i = 0; i < 2; i++) {
    int iref = -1;
    for (int j = 0; j < 2; j++)
      if (ozat[i] == zat[j]) iref = j;
    if (iref < 0) {
      printf("FAIL atom %d has unexpected atomic number %d\n", i, ozat[i]);
      nfail++;
      continue;
    }
    for (int k = 0; k < 3; k++) {
      char what[64];
      sprintf(what, "atom %d (Z=%d) position %d", i, ozat[i], k);
      check_double(what, oposition[i][k], position[iref][k], 1e-10);
    }
  }

  c2_destroy_crystal(cr);

  if (nfail > 0) {
    printf("libcritic2 C API test: %d check(s) failed\n", nfail);
    return 1;
  }
  printf("libcritic2 C API test: all checks passed\n");
  return 0;
}
