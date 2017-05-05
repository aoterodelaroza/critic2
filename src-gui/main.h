/*
Copyright (c) 2017 Alberto Otero de la Roza
<aoterodelaroza@gmail.com>, Robin Myhr <x@example.com>, Isaac
Visintainer <x@example.com>, Richard Greaves <x@example.com>, Ángel
Martín Pendás <angel@fluor.quimica.uniovi.es> and Víctor Luaña
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

#include "imgui.h"
#include "matrix_math.h"

// Bond and atom resolutions (0 = coarse -> 3 = smooth)
extern const char bondresolution;
extern const char atomresolution;

// Bond thickness and atom/CP size
extern const float bondthickness;
extern const float atomsize;
extern const float cpsize;

// Tooltipdelay
extern const float ttipdelay;

// Show/hide elements of the interface
extern bool show_bonds;
extern bool show_cps;
extern bool show_atoms;
extern bool show_cell;

// Quit flag
extern bool want_quit;

// Current state of the camera
extern CameraInfo cam;

