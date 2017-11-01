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

  if (iscene > 0){
    // initialize the camera
    c2::set_scene_pointers(iscene);
    aview->cam = new Camera();
    aview->mstate = new MouseState;
    aview->cam->SetSceneRad(c2::scenerad);
    aview->Update();
  }

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

    // this->Update();

    glBindFramebuffer(GL_FRAMEBUFFER, this->FBO);
    if (ImageInteractive((void *) this->FBOtex,this->mstate) && this->mstate->hover){
      // p.ProcessMouseEvents(&mstate);
    }
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

    glUniformMatrix4fv(glGetUniformLocation(this->shader->id, "wvp"), 1, GL_FALSE, value_ptr(this->cam->m_wvp));
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
      if (this->cam)
	delete this->cam;
      viewlist.erase(it);
      break;
    }
  }
}
