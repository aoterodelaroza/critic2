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

class Pipeline
{
public:
  Pipeline(float scenerad){
    m_scale[0] = 1.f; m_scale[1] = 1.f, m_scale[2] = 1.f;
    m_pos[0] = 0.f; m_pos[1] = 0.f, m_pos[2] = 0.f;
    m_rotate[0] = 0.f; m_rotate[1] = 0.f, m_rotate[2] = 0.f;
    m_post_rotate[0] = 0.f; m_post_rotate[1] = 0.f, m_post_rotate[2] = 0.f;
    m_scenerad = scenerad;
    SetCamera(0.f,0.f,-4.0f*scenerad, 0.f,0.f,1.f, 0.f,1.f,0.f);
    pers_FOV = zfov;
    pers_Width = FBO_tex_x;
    pers_Height = FBO_tex_y;
    pers_zNear = znear;
    pers_zFar = zfar;
  }

  void Scale(float x, float y, float z){
    m_scale[0] = x; m_scale[1] = y; m_scale[2] = z;
  }

  void Translate(float x, float y, float z){
    m_pos[0] = x; m_pos[1] = y; m_pos[2] = z;
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

  const glm::mat4 * GetProjTrans();
  const glm::mat4 * GetViewTrans();
  const glm::mat4 * GetWorldTrans();

  void ProcessMouseEvents(MouseState *mstate);

private:
  float m_scenerad;
  float m_scale[3];
  float m_pos[3];
  float m_rotate[3];
  glm::mat4 m_rotate_trans;
  float m_post_rotate[3];
  glm::mat4 m_post_rotate_trans;

  glm::vec2 mpos0;
  glm::vec2 cpos0;
  glm::mat4 rot0;

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

  glm::mat4 m_ProjTransformation;
  glm::mat4 m_Vtransformation;
  glm::mat4 m_Wtransformation;
};

#endif
