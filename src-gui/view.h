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

#ifndef VIEW_H
#define VIEW_H

#include "imgui/imgui_dock.h"
#include "imgui/mouse.h"
#include "shader.h"
#include "imgui/gl3w.h"

#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/glm.hpp>

using namespace std;
using namespace glm;

struct View
{
  // view methods
  void Draw();
  void Update();
  void Delete();
  bool processMouseEvents();

  // camera methods
  void updateMProjection();
  void updateMView();
  void updateMWVP();

  // camera matrices and vectors
  vec3 v_pos; // position vector
  vec3 v_front; // front vector
  vec3 v_up; // up vector
  mat4 m_projection; // projection
  mat4 m_view; // view
  mat4 m_world; // world
  mat4 m_wvp; // projection * view * world

  // assoicated objects
  GLuint FBO; // framebuffer object
  GLuint FBOtex; // framebuffer object texture
  Shader *shader; // pointer to the current shader
  char *title; // title
  int iscene; // integer identifier of the associated scene
  MouseState *mstate = nullptr; // mouse
  ImGui::Dock *dock; // dock
};
  
// Create a new view linked to scene iscene (0 for no scene).
void CreateView(char *title, Shader *shader, int iscene = 0);

// Draw all available views
void DrawAllViews();

#endif
