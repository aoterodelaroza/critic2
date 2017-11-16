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
#include "settings.h"

#include <glm/gtx/rotate_vector.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/glm.hpp>

using namespace std;
using namespace glm;

struct View
{
  // mouse behavior enum
  enum MouseBehavior_{MB_Navigation,MB_Pointer,MB_Angle,MB_Ruler,MB_Builder,MB_Alignment};

  // view methods
  void Draw();
  void Update();
  void Delete();
  bool processMouseEvents(bool hover);
  bool navigate(bool hover);
  bool updateTexSize();

  // camera methods
  void updateProjection();
  void updateView();
  void updateWorld();
  float getDepth(vec2 ntexpos);
  vec3 sphereProject(vec2 ntexpos);
  vec4 cam_world_coords();

  // projection and unprojection
  void pos_to_ntexpos(vec2 &pos); // screen to normalized texture position
  void ntexpos_to_pos(vec2 &pos); // normalized texture to screen position
  void pos_to_texpos(vec2 &pos); // screen to texture position
  void texpos_to_pos(vec2 &pos); // texture to screen position
  void texpos_to_ntexpos(vec2 &pos); // texture to normalized texture position
  void ntexpos_to_texpos(vec2 &pos); // normalized texture to texture position
  vec2 world_to_texpos(vec3 pos); // world to texture position

  // draw shapes
  void drawSphere(float r0[3],float rad,float rgb[4]);
  void drawCylinder(float r1[3],float r2[3],float rad,float rgb[4]);

  // camera matrices and vectors
  bool iswire = false; // use wire
  bool isortho = true; // is ortho or perspective?
  vec3 v_pos = {}; // position vector
  vec3 v_front = {}; // front vector
  vec3 v_up = {}; // up vector
  mat4 m_projection = mat4(1.0); // projection
  mat4 m_view = mat4(1.0); // view
  mat4 m_world = mat4(1.0); // world

  // saved states
  MouseBehavior_ mousebehavior = MB_Navigation; // mouse behavior
  bool rlock = false; // rmb is dragging
  bool llock = false; // lmb is dragging
  vec3 mpos0; // saved mouse position in screen coords (0 -> 1024)
  vec3 cpos0; // saved camera position in world coords
  mat4 crot0; // saved camera rotation

  // associated objects
  int icurtex = -1; // current texture in use
  GLuint FBO[nmaxtex]; // framebuffer object
  GLuint FBOtex[nmaxtex]; // framebuffer object textures
  GLuint FBOdepth[nmaxtex]; // framebuffer object depth buffers
  Shader *shader = nullptr; // pointer to the current shader
  char *title; // title
  int iscene = -1; // integer identifier of the associated scene
  ImGui::Dock *dock = nullptr; // dock
  ImRect vrect; // rectangle for the current view
};
  
// Create a new view linked to scene iscene (0 for no scene).
View *CreateView(char *title, Shader *shader, int iscene = 0);

// Draw all available views
void DrawAllViews();

#endif
