// -*- c++ -*-
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
#include "imgui/gl3w.h"
#include "settings.h"

#include <glm/gtx/rotate_vector.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/glm.hpp>

struct View
{
  // mouse behavior enum
  enum MouseBehavior_{MB_Navigation,MB_Pointer,MB_Angle,MB_Ruler,
		      MB_Builder,MB_Alignment,MB_Query};

  // view methods
  void SetDefaults();
  void Draw();
  void Update();
  void Delete();
  bool processMouseEvents(bool hover);
  bool Navigate(bool hover);

  // texture methods
  void createTex(float atex);
  void deleteTex();
  bool updateTexSize();

  // camera methods
  void resetView();
  void updateProjection();
  void updateView();
  void updateWorld();
  glm::vec3 cam_world_coords();
  glm::vec3 cam_view_coords();

  // coordinate transformations
  // pos: mouse coordinates
  // texpos: texture coordinates (bl:0,0 | ur:texsize.x,texsize.y)
  // ntexpos: normalized texture coordinates (bl:-1,-1 | ur:1,1)
  // world: world (scene) coordinates
  // view: view (camera) coordinates
  void pos_to_ntexpos(glm::vec2 &pos); 
  void ntexpos_to_pos(glm::vec2 &pos);
  void pos_to_texpos(glm::vec2 &pos);
  void texpos_to_pos(glm::vec2 &pos);
  void texpos_to_ntexpos(glm::vec2 &pos);
  void ntexpos_to_texpos(glm::vec2 &pos);
  glm::vec2 world_to_texpos(glm::vec3 pos);
  glm::vec3 texpos_to_world(glm::vec2 pos, float dist=-1.f); // dist=0, znear; dist<0, scene origin plane; dist>0, distance from camera
  glm::vec2 world_to_ntexpos(glm::vec3 pos);
  glm::vec3 ntexpos_to_world(glm::vec2 pos, float dist=-1.f);  // dist=0, znear; dist<0, scene origin plane; dist>0, distance from camera
  glm::vec2 view_to_texpos(glm::vec3 pos, float *depth); // returns depth in the second argument
  glm::vec3 texpos_to_view(glm::vec2 pos, float depth); 
  float texpos_viewdepth(glm::vec2 texpos); // depth from current view at texpos or 1.0 if background pixel found
  
  // draw shapes
  void drawSphere(glm::vec3 r0, float rad, glm::vec4 rgb, int res, bool blend);
  void drawCylinder(glm::vec3 r1, glm::vec3 r2, float rad, glm::vec4 rgb, int res, bool blend);

  // draw settings
  bool isucell; // draw the unit cell?
  bool isborder; // draw the atoms on the border?
  bool ismotif; // draw the molecular motif
  int ncell[3]; // number of unit cells in each direction

  // view settings
  float resetd; // reset distance (scenerad)
  float zfov; // field of view angle (degrees)
  float bgrgb[4]; // background color
  int isphres; // atom resolution
  int icylres; // bond resolution

  // camera matrices and vectors
  bool iswire = false; // use wire
  bool isortho = false; // is ortho or perspective?
  glm::vec3 v_pos = {}; // position vector
  glm::vec3 v_front = {}; // front vector
  glm::vec3 v_up = {}; // up vector
  glm::mat4 m_projection = glm::mat4(1.0); // projection
  glm::mat4 m_view = glm::mat4(1.0); // view
  glm::mat4 m_world = glm::mat4(1.0); // world

  // saved states for the mouse interaction
  MouseBehavior_ mousebehavior = MB_Navigation; // mouse behavior
  glm::vec2 mposlast; // last mouse position
  bool rlock = false; // dragging
  glm::vec3 mpos0_r; // saved mouse position (for dragging)
  glm::vec3 cpos0_r; // saved camera position in view coords (for dragging)

  bool llock = false; // lmb is dragging
  glm::vec3 mpos0_l; // saved mouse position (for rotating)
  glm::vec3 cpos0_l; // saved camera position in view coords (for rotating)
  glm::mat4 crot0_l; // saved world matrix (for rotating)

  bool slock = false; // scroll dragging
  float mpos0_s; // mouse y position (for scrolling using keys)

  bool updatescene = false; // update the scene next pass

  // associated objects
  GLuint FBO; // framebuffer object
  GLuint FBOtex; // framebuffer object textures
  GLuint FBOdepth; // framebuffer object depth buffers
  float FBO_atex; // side of the texture (pixels)
  float FBO_a; // side of the texture square used for rendering (pixels)
  char *title; // title
  int iscene = -1; // integer identifier of the associated scene
  ImGui::Dock *dock = nullptr; // dock
  ImRect vrect; // rectangle for the current view
};
  
// Create a new view linked to scene iscene (0 for no scene).
View *CreateView(char *title, int iscene = 0);

// Draw all available views
void DrawAllViews();

// Force-update all views
void ForceUpdateAllViews();

// Set all views settings to default value 
void SetDefaultAllViews();

#endif
