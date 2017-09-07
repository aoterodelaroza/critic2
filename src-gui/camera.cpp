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

#include "camera.h"
#include <stdio.h>

void Camera::processMouseEvents(MouseState *m){
  static bool firstpass = true;
  if (firstpass){
    firstpass = false;
    pmstate = *m;
    return;
  }

  // save mouse positions for possible drag
  if ((m->rclick || m->lclick) && !(m->ldrag || m->rdrag))
    posdrag = m->pos;

  if (m->rdrag){
    vec3 Right = normalize(cross(LookAt-Position,Up));
    vec2 offset = mousesens_pan * (m->pos - pmstate.pos);
    LookAt += - offset[0] * Right + offset[1] * Up;
    Position += - offset[0] * Right + offset[1] * Up;
  }

  // zoom in and out (requires using box_xmaxlen to scale)
  if (abs(m->scroll) > 1e-6)
    Position += normalize(LookAt-Position) * m->scroll * mousesens_zoom;

  pmstate = *m;
}

