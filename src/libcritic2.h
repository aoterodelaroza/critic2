/*
Copyright (c) 2015-2022 Alberto Otero de la Roza <aoterodelaroza@gmail.com>,
Ángel Martín Pendás <angel@fluor.quimica.uniovi.es> and Víctor Luaña
<victor@fluor.quimica.uniovi.es>.

critic2 is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at
your option) any later version.

critic2 is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef LIBCRITIC2_H
#define LIBCRITIC2_H

#ifdef __cplusplus
extern "C" {
#endif

//// Types ////
typedef void crystal;

//// Functions ////
// Read structure from a file and return a pointer to the crystal
// structure object (requires destruction after use).
crystal *c2_crystal_from_file(const char *file);
crystal *c2_crystal_from_lattice(const int natom,const double lattice[3][3],
				 const double position[][3],const int zat[]);
crystal *c2_crystal_from_cellpar(const int natom,const double cel[3], const double ang[3],
				 const double position[][3],const int zat[]);

// Write a report about the crystal structure to standard output.
void c2_describe_crystal(crystal *cr);

// Destroy a crystal structure object.
void c2_destroy_crystal(crystal *cr);

#ifdef __cplusplus
}
#endif
#endif
