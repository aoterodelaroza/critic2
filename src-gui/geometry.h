// Copyright (c) 2017 Alberto Otero de la Roza <aoterodelaroza@gmail.com>,
// Ángel Martín Pendás <angel@fluor.quimica.uniovi.es> and Víctor Luaña
// <victor@fluor.quimica.uniovi.es>. 
//
// critic2 is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or (at
// your option) any later version.
// 
// critic2 is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include <GL/gl3w.h>

// Function declarations
void CreateAndFillBuffers();

// Pointers and array sizes for sphere and cylinder models
extern GLfloat *isphv[4];
extern unsigned int *isphi[4];
extern unsigned int nsphv[4];
extern unsigned int nsphi[4];
extern GLfloat *icylv[4];
extern unsigned int *icyli[4];
extern unsigned int ncylv[4];
extern unsigned int ncyli[4];

// Buffer indices
extern GLuint bufsphi[4], bufsphv[4];
extern GLuint bufcyli[4], bufcylv[4];

