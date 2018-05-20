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

#ifndef SCENE_H
#define SCENE_H

#include "settings.h"
#include <glm/glm.hpp>

#include "critic2.h"

struct Scene{
  // basic scene info
  int iscene = -1; // id for the associated scene

  // from critic2
  bool ismolecule; // true if this is a molecular system, false if crystal
  float scenerad; // scene radius
  float avec[3][3];

  // draw settings
  bool isucell; // draw the unit cell?
  bool ismolcell; // draw the molecular cell?
  bool isborder; // draw the atoms on the border?
  bool ismotif; // draw the molecular motif
  int ncell[3]; // number of unit cells in each direction

  // view settings
  float resetd; // reset distance (scenerad)
  float zfov; // field of view angle (degrees)
  bool show_atoms; // show atoms?
  float scale_atoms; // global scale atoms
  int isphres; // atom resolution
  bool show_bonds; // show bonds?
  float scale_bonds; // global scale bonds
  int icylres; // bond resolution
  bool iswire = false; // use wire
  bool isortho = false; // is ortho or perspective?

  // camera matrices and vectors
  glm::vec3 v_pos = {}; // position vector
  glm::vec3 v_front = {}; // front vector
  glm::vec3 v_up = {}; // up vector
  glm::mat4 m_projection = glm::mat4(1.0); // projection
  glm::mat4 m_view = glm::mat4(1.0); // view
  glm::mat4 m_world = glm::mat4(1.0); // world

  // True if the scene needs updating
  bool updatescene = false; // update the scene next pass
  bool setshaderp = false; // send the projection matrix to the shader next pass
  bool setshaderv = false; // send the view matrix to the shader next pass
  bool setshaderw = false; // send the world matrix to the shader next pass

  // Constructor
  Scene(int isc){
    iscene = isc;
    grabFromC2();
    setDefaults();
    resetView();
    updateAll();  
  };

  void grabFromC2(); // Get the scene parameters from the critic2 interface
  void setDefaults(); // Set the scene to the default settings (does not update)
  void resetView(); // Reset the view parameters (does not update)
  void updateAll(); // Update all matrices
  void updateProjection(); // Update the projection matrix
  void updateView(); // Update the view matrix
  void updateWorld(); // Update the world matrix

  // camera methods
  bool alignViewAxis(int iaxis);
  glm::vec3 cam_world_coords();
  glm::vec3 cam_view_coords();

  
};

#endif
