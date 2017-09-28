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
#ifdef WIN32

#define _USE_MATH_DEFINES
#endif // WIN32

#include <stdlib.h>
#include <matrix_math.h>

using namespace std;

void Pipeline::ProcessMouseEvents(MouseState *mstate){
  // scroll
  if (abs(mstate->scroll) > 1e-5){
    m_camera_pos[2] += mstate->scroll * mousesens_zoom * m_scenerad;
    m_camera_pos[2] = fmin(m_camera_pos[2],-znear);
  }

  // drag
  if (mstate->rclick){ 
    mpos0 = mstate->pos;
    cpos0 = {m_camera_pos[0],m_camera_pos[1]};
  } else if (mstate->rdrag){
    glm::vec2 dx = mstate->pos - mpos0;
    m_camera_pos[0] = cpos0.x + mousesens_pan * m_scenerad * dx.x;
    m_camera_pos[1] = cpos0.y - mousesens_pan * m_scenerad * dx.y;
  }

  // rotate
  if (mstate->lclick){ 
    mpos0 = mstate->pos;
    rot0 = m_post_rotate_trans;
  } else if (mstate->ldrag) {
    glm::vec3 curRotAxis = glm::vec3((float)(mpos0.x-mstate->pos.x), (float)(mstate->pos.y-mpos0.y), 0.f);
    curRotAxis = glm::cross(curRotAxis,glm::vec3(0, 0, 1));
    float curRotAng = glm::length(curRotAxis) * mousesens_rot * m_scenerad;
    curRotAxis = glm::normalize(curRotAxis);
 
    glm::mat4 curRot = glm::rotate(glm::mat4(),curRotAng,glm::vec3(curRotAxis.x,curRotAxis.y,curRotAxis.z));
    m_post_rotate_trans = curRot * rot0;
  }
}

const glm::mat4 * Pipeline::GetProjTrans(){

  m_ProjTransformation = glm::infinitePerspective(pers_FOV,pers_Width/pers_Height,pers_zNear);
  // m_ProjTransformation = glm::ortho(-10.f,10.f,-10.f,10.f,0.1f,1000.f);

  return &m_ProjTransformation;
}

const glm::mat4 * Pipeline::GetViewTrans(){
  glm::vec3 pos = glm::vec3(m_camera_pos[0],m_camera_pos[1],m_camera_pos[2]);
  glm::vec3 target = pos + glm::vec3(m_camera_target[0],m_camera_target[1],m_camera_target[2]);
  glm::vec3 up = glm::vec3(m_camera_up[0],m_camera_up[1],m_camera_up[2]);

  m_Vtransformation = glm::lookAt(pos,target,up);

  return &m_Vtransformation;
}

const glm::mat4 * Pipeline::GetWorldTrans(){
  glm::mat4 scale;
  scale = glm::scale(scale,glm::vec3(m_scale[0],m_scale[1],m_scale[2]));

  glm::mat4 translate;
  translate = glm::translate(translate,glm::vec3(m_pos[0],m_pos[1],m_pos[2]));

  m_Wtransformation = m_post_rotate_trans * translate * m_rotate_trans * scale;

  return &m_Wtransformation;
}
