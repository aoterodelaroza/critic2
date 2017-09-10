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

#include "imgui/gl3w.h"
#include <GLFW/glfw3.h>
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
  if ((m->rclick || m->lclick) && !(m->ldrag || m->rdrag)){
    posdrag = m->pos;
    viewdrag = view;
  }

  // rotation with the left button
  if (m->ldrag){
    // float arc = length(m->pos - pmstate.pos);
    // printf("%f\n",arc);
    // vec3 axis = vec3(m->pos.x - pmstate.pos.x,m->pos.y - pmstate.pos.y,0.f);
    // axis = cross(axis,vec3(0.f,0.f,1.f));
    // float angle = length(axis) * mousesens_rot * srad;
    // axis = normalize(axis);
    // mat4 rot = rotate(angle,axis);
    // view = rot * viewdrag;
  }

  // translation with right mouse button
  if (m->rdrag){
    vec3 vr = row(view,0);
    vec3 vu = row(view,1);
    vec2 offset = mousesens_pan * srad * (m->pos - pmstate.pos);
    view = translate(view,offset[0] * vr + offset[1] * vu);
  }

  // zoom by scrolling
  if (abs(m->scroll) > 1e-6){
    vec3 dir = row(view,2);
    view = translate(view,m->scroll * srad * mousesens_zoom * dir);
  }

  // save into the previous mouse state
  pmstate = *m;
}

void Camera::applyMatrices(GLuint shader_id){
  glUniformMatrix4fv(glGetUniformLocation(shader_id, "projection"), 1, GL_FALSE, &projection[0][0]);
  glUniformMatrix4fv(glGetUniformLocation(shader_id, "view"), 1, GL_FALSE, &view[0][0]);
}

