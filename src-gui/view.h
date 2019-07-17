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
#include "scene.h"

#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/rotate_vector.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/glm.hpp>

#include <map>

struct View
{
  // mouse behavior enum
  enum MouseBehavior_{MB_Navigation,MB_Pointer,MB_Angle,MB_Ruler,
		      MB_Builder,MB_Alignment,MB_Query};

  enum Variable_{V_ALL, V_lightpos, V_lightcolor, V_ambient, V_diffuse, V_specular, V_shininess,
                 V_wireframe, V_orthogonal, V_fov, V_resetdistance, V_bgrgb, V_show_atoms,
                 V_isphres, V_show_bonds, V_icylres, V_show_labels, V_format_labels, V_lat_labels,
                 V_scale_labels, V_rgb_labels
  };

  // constructor
  View(char *title_, float atex, int iscene=0);
  ~View();

  // view methods
  void changeScene(int isc);
  void setDefaultAllScenes(Variable_ var);
  void setDefault(Scene *sc_=nullptr, Variable_ var=V_ALL);
  void Draw();
  void Update();
  void Delete();
  bool processMouseEvents(bool hover);
  bool Navigate(bool hover);

  // texture methods
  void createTex(float atex);
  void deleteTex();
  bool updateTexSize();

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
  
  // draw atomic labels (requires setting the c2 scene pointers first)
  void drawAtomLabel(glm::vec3 x0, int iatom, int ix, int iy, int iz);

  // draw shapes
  void drawSphere(glm::vec3 r0, float rad, glm::vec4 rgb, int res);
  void drawCylinder(glm::vec3 r1, glm::vec3 r2, float rad, glm::vec4 rgb, int res);
  void drawUnitCell(glm::vec3 &v0, glm::vec3 &vx, glm::vec3 &vy, glm::vec3 &vz, bool colors);

  // map of scenes
  std::map<int,Scene *> scmap; // map of the known scenes
  Scene *sc = nullptr; // pointer to the current scene
  int iscene = -1; // integer identifier of the current scene

  // saved states for the mouse interaction
  MouseBehavior_ mousebehavior = MB_Navigation; // mouse behavior
  glm::vec2 mposlast; // last mouse position
  bool rlock = false; // dragging
  glm::vec3 mpos0_r; // saved mouse position (for dragging)
  glm::vec3 cpos0_r; // saved camera position in view coords (for dragging)
  bool llock = false; // lmb is dragging
  glm::vec3 mpos0_l; // saved mouse position (for rotating)
  glm::vec3 cpos0_l; // saved camera position in view coords (for rotating)
  bool slock = false; // scroll dragging
  float mpos0_s; // mouse y position (for scrolling using keys)

  // special effects for this view
  float SE_pixelated = 1.0f;

  // associated textures
  GLuint FBO; // framebuffer object
  GLuint FBOtex; // framebuffer object texture
  GLuint FBOdepth; // framebuffer object depth buffer
  float FBO_atex; // side of the render texture (pixels)
  float arender; // side of the square used for rendering (pixels), arender <= FBO_atex

  // dock
  float bgrgb[4]; // background color
  char *title; // title
  ImGui::Dock *dock = nullptr; // dock
  ImRect vrect; // rectangle for the current view
};
  
// Create a new view linked to scene iscene (0 for no scene).
View *CreateView(char *title, int iscene = 0);

// Draw all available views
void DrawAllViews();

// Force-update all views
void ForceUpdateAllViews();

// Set default values for all views
void SetDefaultAllViews(View::Variable_ var=View::V_ALL);

// Main view
extern View *mainview;

#endif
