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
/*
Copyright 2010 Etay Meiri

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef MATRIX_MATH_H
#define MATRIX_MATH_H

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <settings.h>
#include <imgui/mouse.h>

#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/glm.hpp>

using namespace std;
using namespace glm;

class Camera
{
public:
  Camera(float scenerad){
    m_scenerad = scenerad;
    SetCamera(0.f,0.f,-4.0f*scenerad, 0.f,0.f,1.f, 1.f,0.f,0.f);
    pers_FOV = zfov;
    pers_Width = FBO_tex_x;
    pers_Height = FBO_tex_y;
    pers_zNear = znear;
    pers_zFar = zfar;
    UpdateProjection();
    UpdateView();
    UpdateWVP();
  }

  void SetCamera(float Pos[3], float Target[3], float Up[3]){
    memcpy(&m_camera_pos, Pos, sizeof(float)*3);
    memcpy(&m_camera_target, Target, sizeof(float)*3);
    memcpy(&m_camera_up, Up, sizeof(float)*3);
  }

  void SetCamera(float px,float py,float pz,float tx,float ty,float tz,float ux,float uy,float uz){
    m_camera_pos[0] = px; m_camera_pos[1] = py; m_camera_pos[2] = pz;
    m_camera_target[0] = tx; m_camera_target[1] = ty; m_camera_target[2] = tz;
    m_camera_up[0] = ux; m_camera_up[1] = uy; m_camera_up[2] = uz;
  }

  void SetScenerad(float scenerad){
    m_scenerad = scenerad;
  }

  void SetOrthoProjInfo(float l, float r, float b, float t, float n, float f){
    ortho_Left = l;
    ortho_Right = r;
    ortho_Bottom = b;
    ortho_Top = t;
    ortho_zNear = n;
    ortho_zFar = f;
  }

  void UpdateProjection();
  void UpdateView();
  void UpdateWVP();

  void ProcessMouseEvents(MouseState *mstate);

  mat4 m_projection;
  mat4 m_view;
  mat4 m_world;
  mat4 m_wvp;
private:
  float m_scenerad;

  vec2 mpos0;
  vec2 cpos0;
  mat4 rot0;

  float pers_FOV;
  float pers_Width;
  float pers_Height;
  float pers_zNear;
  float pers_zFar;
  float ortho_Right;
  float ortho_Left;
  float ortho_Bottom;
  float ortho_Top;
  float ortho_zNear;
  float ortho_zFar;

  float m_camera_pos[3];
  float m_camera_target[3];
  float m_camera_up[3];
};

#endif
