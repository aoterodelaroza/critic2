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

#include <list>

#include "imgui/imgui.h"
#include "imgui/imgui_widgets.h"
#include "imgui/mouse.h"

#include "settings.h"
#include "shapes.h"
#include "view.h"
#include "critic2.h"

using namespace ImGui;
using namespace std;

// A linked list for all current views.
static list<View*> viewlist;

void CreateView(char *title, Shader *shader, int iscene/*=0*/){
  View *aview = new View;

  // save the scene id and set the pointers to that scene
  aview->iscene = iscene;
  aview->title = title;
  aview->shader = shader;

  // generate the texture buffer
  glGenFramebuffers(1, &aview->FBO);
  glGenTextures(1, &aview->FBOtex);
  glBindFramebuffer(GL_FRAMEBUFFER, aview->FBO);
  glBindTexture(GL_TEXTURE_2D, aview->FBOtex);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, FBO_tex_x, FBO_tex_y, 0, GL_RGB, GL_UNSIGNED_BYTE, NULL);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);  
  glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, aview->FBOtex, 0); 
  if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
    exit(EXIT_FAILURE);
  glBindFramebuffer(GL_FRAMEBUFFER, 0);

  // initialize the camera vectors
  aview->v_pos    = {0.f,0.f,-10.f};
  aview->v_front  = {0.f,0.f,1.f};
  aview->v_up     = {1.f,0.f,0.f};

  // new mouse state
  aview->mstate = new MouseState;

  if (iscene > 0){
    // set the camera to fit the scene size
    c2::set_scene_pointers(iscene);
    aview->v_pos = {0.f,0.f,-4.f*c2::scenerad};
  }

  // initialize the camera matrices
  aview->updateMProjection();  
  aview->updateMView();
  aview->updateMWVP();

  // plot the scene to the texture
  aview->Update();

  // add the view to the list
  viewlist.push_back(aview);
}

void DrawAllViews(){
  for (auto iv : viewlist)
    iv->Draw();
}

void View::Draw(){
  ImGuiContext *g = GetCurrentContext();
  PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(0.0,g->Style.WindowPadding.y));
  PushStyleColor(ImGuiCol_WindowBg,ImVec4(bgrgb[0],bgrgb[1],bgrgb[2],bgrgb[3]));
  // xxxx //
  SetNextWindowSize(ImVec2(300.f,300.f));
  // xxxx //
  if (BeginDock("Main view") && this->iscene > 0){
    // set the pointers to the current scene
    c2::set_scene_pointers(this->iscene);

    glBindFramebuffer(GL_FRAMEBUFFER, this->FBO);
    if (ImageInteractive((void *) this->FBOtex,this->mstate) && this->mstate->hover)
      if (processMouseEvents())
        this->Update();

    glBindFramebuffer(GL_FRAMEBUFFER, 0);

    // // mouse to world coordinates
    // if (this->mstate->hover && this->mstate->ldown){
    //   glm::vec4 vndc = {0.f,0.f,0.f,1.f};
    //   vndc = this->cam->m_wvp * vndc;
    //   vndc.x = this->mstate->ndpos.x * vndc.w;
    //   vndc.y = this->mstate->ndpos.y * vndc.w;
    //   vndc = inverse(this->cam->m_wvp) * vndc;
    //   vndc /= vndc.w;
    //   float rgb[3] = {0.f,1.f,0.f};
    //   mat4 m_model;
    //   m_model = translate(m_model,vec3(vndc.x,vndc.y,vndc.z));
    //   glUniform4fv(glGetUniformLocation(this->shader->id, "vColor"), 1, (const GLfloat *)rgb);
    //   glUniformMatrix4fv(glGetUniformLocation(this->shader->id, "model"), 1, GL_FALSE, value_ptr(m_model));
    //   glDrawElements(GL_TRIANGLES, spherenel[0], GL_UNSIGNED_INT, 0);
    // }
    
    // // world to mouse coordinates
    // glm::vec4 vndc = {10.f,10.0f,10.f,1.f};
    // printf("ball1: %f %f %f\n",vndc.x,vndc.y,vndc.z);
    // float rgb[3] = {0.f,1.f,0.f};
    // mat4 m_model;
    // m_model = translate(m_model,vec3(vndc.x,vndc.y,vndc.z));
    // glUniform4fv(glGetUniformLocation(shader.id, "vColor"), 1, (const GLfloat *)rgb);
    // glUniformMatrix4fv(glGetUniformLocation(shader.id, "model"), 1, GL_FALSE, value_ptr(m_model));
    // glDrawElements(GL_TRIANGLES, spherenel[0], GL_UNSIGNED_INT, 0);
    
    // glm::vec4 outpos = p.m_wvp * vndc;
    // outpos /= outpos.w;

    // printf("screen1 : %f %f %f %f\n",outpos.x,outpos.y,outpos.z,outpos.w);
  }
  this->dock = GetCurrentDock();
  EndDock();
  PopStyleColor();
  PopStyleVar();
  
  // dviewcont->newDock(dviewdock);
  // dviewdock->setDetachedDockSize(300.f,300.f);
}

void View::Update(){

  glBindFramebuffer(GL_FRAMEBUFFER, this->FBO);
  glViewport(0.,0.,FBO_tex_x,FBO_tex_y);

  glClearColor(bgrgb[0],bgrgb[1],bgrgb[2],bgrgb[3]);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
      
  if (this->iscene > 0){
    c2::set_scene_pointers(this->iscene);
    glBindVertexArray(sphereVAO[isphres]);

    glUniformMatrix4fv(glGetUniformLocation(this->shader->id, "wvp"), 1, GL_FALSE, value_ptr(this->m_wvp));
    for (int i=0;i<c2::nat;i++){
      mat4 m_model;
      m_model = translate(m_model,vec3(c2::at[i].r[0],c2::at[i].r[1],c2::at[i].r[2]));
      m_model = scale(m_model,vec3(c2::at[i].rad,c2::at[i].rad,c2::at[i].rad));

      glUniform4fv(glGetUniformLocation(this->shader->id, "vColor"), 1, (const GLfloat *)c2::at[i].rgb);
      glUniformMatrix4fv(glGetUniformLocation(this->shader->id, "model"), 1, GL_FALSE, value_ptr(m_model));
      glDrawElements(GL_TRIANGLES, spherenel[isphres], GL_UNSIGNED_INT, 0);
    }
  }
  glBindFramebuffer(GL_FRAMEBUFFER, 0);
}

void View::Delete(){
  for (auto it = viewlist.begin(); it != viewlist.end(); it++) {
    if (*it == this){
      glDeleteTextures(1, &this->FBOtex);
      glDeleteFramebuffers(1, &this->FBO); 
      if (this->mstate)
	delete this->mstate;
      viewlist.erase(it);
      break;
    }
  }
}

bool View::processMouseEvents(){
  c2::set_scene_pointers(iscene);

  bool updateview = false, updatewvp = false;

  // mouse scroll = zoom
  if (abs(mstate->scroll) > 1e-5){
    v_pos[2] += mstate->scroll * mousesens_zoom * c2::scenerad;
    v_pos[2] = fmin(v_pos[2],-znear);
    updateview = true;
    updatewvp = true;
  }

  // // drag
  // if (mstate->rclick){ 
  //   mpos0 = mstate->pos;
  //   cpos0 = {m_camera_pos[1],m_camera_pos[0]};
  // } else if (mstate->rdrag){
  //   vec2 dx = mstate->pos - mpos0;
  //   m_camera_pos[0] = cpos0.y + mousesens_pan * m_scenerad * dx.y;
  //   m_camera_pos[1] = cpos0.x - mousesens_pan * m_scenerad * dx.x;
  //   UpdateView();
  //   UpdateWVP();
  // }

  // // rotate
  // if (mstate->lclick){ 
  //   mpos0 = mstate->pos;
  //   rot0 = m_world;
  // } else if (mstate->ldrag) {
  //   vec3 curRotAxis = vec3((float)(mpos0.y-mstate->pos.y), (float)(mstate->pos.x-mpos0.x), 0.f);
  //   curRotAxis = cross(curRotAxis,vec3(0, 0, 1));
  //   float curRotAng = length(curRotAxis) * mousesens_rot * m_scenerad;
  //   curRotAxis = normalize(curRotAxis);
 
  //   mat4 curRot = rotate(mat4(),curRotAng,vec3(curRotAxis.x,curRotAxis.y,curRotAxis.z));
  //   m_world = curRot * rot0;
  //   UpdateWVP();
  // }

  if (updateview)
    updateMView();
  if (updateview || updatewvp)
    updateMWVP();

  return updateview || updatewvp;
}

void View::updateMProjection(){
  m_projection = infinitePerspective(zfov,FBO_tex_x/FBO_tex_y,znear);
}
void View::updateMView(){
  m_view = lookAt(v_pos,v_pos+v_front,v_up);
}
void View::updateMWVP(){
  m_wvp = m_projection * m_view * m_world;
}

