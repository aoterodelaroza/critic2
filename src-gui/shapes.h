// -*-c++-*-
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

#ifndef SHAPES_H
#define SHAPES_H

#include "imgui/gl3w.h"

// public vertex buffer and attribute objects
const int nmaxsph = 4;
extern GLuint sphVAO[];
extern const GLuint sphnel[];
const int nmaxcyl = 3;
extern GLuint cylVAO[];
extern const GLuint cylnel[];

// shape functions
void CreateAndFillBuffers();
void DeleteBuffers();

#endif
