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

// The viewport vector
const vec4 viewport = {0.f,0.f,FBO_tex_a,FBO_tex_a};

void CreateView(char *title, Shader *shader, int iscene/*=0*/){
  View *aview = new View;

  // save the scene id and set the pointers to that scene
  aview->iscene = iscene;
  aview->title = title;
  aview->shader = shader;

  // texture
  glGenTextures(1, &aview->FBOtex);
  glBindTexture(GL_TEXTURE_2D, aview->FBOtex);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, FBO_tex_a, FBO_tex_a, 0, GL_RGBA, GL_UNSIGNED_BYTE, 0);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);  
  glBindTexture(GL_TEXTURE_2D, 0);

  // render buffer
  glGenRenderbuffers(1, &aview->FBOdepth);
  glBindRenderbuffer(GL_RENDERBUFFER, aview->FBOdepth);
  glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT, FBO_tex_a, FBO_tex_a);
  glBindRenderbuffer(GL_RENDERBUFFER, 0);

  // frame buffer
  glGenFramebuffers(1, &aview->FBO);
  glBindFramebuffer(GL_FRAMEBUFFER, aview->FBO);
  glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, aview->FBOtex, 0); 
  glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, aview->FBOdepth);

  if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
    exit(EXIT_FAILURE);
  glBindFramebuffer(GL_FRAMEBUFFER, 0);

  // initialize the camera vectors
  aview->v_pos    = {0.f,0.f,10.f};
  aview->v_front  = {0.f,0.f,-1.f};
  aview->v_up     = {0.f,1.f,0.f};

  // new mouse state
  aview->mstate = new MouseState;

  if (iscene > 0){
    // set the camera to fit the scene size
    c2::set_scene_pointers(iscene);
    aview->v_pos = {0.f,0.f,4.f*c2::scenerad};
  }

  // initialize the camera matrices
  aview->updateProjection();  
  aview->updateView();
  aview->updateWorld();

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
  SetNextWindowSize(ImVec2(300.f,300.f),ImGuiSetCond_Once);
  if (BeginDock("Main view") && iscene > 0){
    ImageInteractive((void *) FBOtex,mstate);
    if (processMouseEvents())
      Update();
  }
  dock = GetCurrentDock();
  EndDock();
  PopStyleColor();
  PopStyleVar();
}

void View::Update(){

  glBindFramebuffer(GL_FRAMEBUFFER, FBO);
  glViewport(0.,0.,FBO_tex_a,FBO_tex_a);

  glClearColor(bgrgb[0],bgrgb[1],bgrgb[2],bgrgb[3]);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
      
  if (iscene > 0){
    c2::set_scene_pointers(iscene);
    glBindVertexArray(sphereVAO[isphres]);

    // scene
    for (int i=0;i<c2::nat;i++){
      mat4 m_model = mat4(1.0f);
      m_model = translate(m_model,vec3(c2::at[i].r[0],c2::at[i].r[1],c2::at[i].r[2]));
      m_model = scale(m_model,vec3(c2::at[i].rad,c2::at[i].rad,c2::at[i].rad));

      shader->setVec4("vColor",(const GLfloat *)c2::at[i].rgb);
      shader->setMat4("model",value_ptr(m_model));
      glDrawElements(GL_TRIANGLES, spherenel[isphres], GL_UNSIGNED_INT, 0);
    }
    // coordinate axes
    mat4 m_model = mat4(1.0f);
    vec3 rgb = vec3(1.f,0.f,0.f);
    m_model = translate(m_model,vec3(c2::scenerad*1.2,0.f,0.f));
    m_model = scale(m_model,vec3(1.f,1.f,1.f));
    shader->setVec4("vColor",value_ptr(rgb));
    shader->setMat4("model",value_ptr(m_model));
    glDrawElements(GL_TRIANGLES, spherenel[isphres], GL_UNSIGNED_INT, 0);

    m_model = mat4(1.0f);
    rgb = vec3(0.f,1.f,0.f);
    m_model = translate(m_model,vec3(0.f,c2::scenerad*1.2,0.f));
    m_model = scale(m_model,vec3(1.f,1.f,1.f));
    shader->setVec4("vColor",value_ptr(rgb));
    shader->setMat4("model",value_ptr(m_model));
    glDrawElements(GL_TRIANGLES, spherenel[isphres], GL_UNSIGNED_INT, 0);

    m_model = mat4(1.0f);
    rgb = vec3(0.f,0.f,1.f);
    m_model = translate(m_model,vec3(0.f,0.f,c2::scenerad*1.2));
    m_model = scale(m_model,vec3(1.f,1.f,1.f));
    shader->setVec4("vColor",value_ptr(rgb));
    shader->setMat4("model",value_ptr(m_model));
    glDrawElements(GL_TRIANGLES, spherenel[isphres], GL_UNSIGNED_INT, 0);
  }

  glBindFramebuffer(GL_FRAMEBUFFER, 0);
}

void View::Delete(){
  for (auto it = viewlist.begin(); it != viewlist.end(); it++) {
    if (*it == this){
      glDeleteTextures(1, &FBOtex);
      glDeleteFramebuffers(1, &FBO); 
      if (mstate)
	delete mstate;
      viewlist.erase(it);
      break;
    }
  }
}

bool View::processMouseEvents(){
  c2::set_scene_pointers(iscene);

  const float eps = 1e-8;
  bool updateview = false, updateworld = false, updateprojection = false;

  // mouse scroll = zoom
  if (mstate->hover && abs(mstate->scroll) > eps){
    float ratio = fmin(mousesens_zoom * mstate->scroll,0.5f);
    v_pos = v_pos - ratio * v_pos;
    if (length(v_pos) < min_zoom)
      v_pos = v_pos / length(v_pos) * min_zoom;
    if (isortho)
      updateprojection = true;
    updateview = true;
  }

  // drag
  if (mstate->hover && mstate->rclick){ 
    float depth = getDepth(mstate->ndpos);
    if (depth < 1.0){
      mpos0 = {mstate->ndpos.x*FBO_tex_a,mstate->ndpos.y*FBO_tex_a,depth};
    }else{
      vec3 origin = {0.f,0.f,0.f};
      origin = project(origin,m_view,m_projection,viewport);
      mpos0 = {mstate->ndpos.x*FBO_tex_a,mstate->ndpos.y*FBO_tex_a,origin.z};
    }
    cpos0 = {v_pos[0],v_pos[1],0.f};
    rlock = true;
  } else if (rlock){
    if (mstate->rdown){
      vec3 vnew = {mstate->ndpos.x*FBO_tex_a,mstate->ndpos.y*FBO_tex_a,mpos0.z};
      vec3 vold = mpos0;
      vnew = unProject(vnew,m_view,m_projection,viewport);
      vold = unProject(vold,m_view,m_projection,viewport);
      v_pos.x = cpos0.x - (vnew.x - vold.x);
      v_pos.y = cpos0.y - (vnew.y - vold.y);
      updateview = true;
    } else {
      rlock = false;
    }
  }

  // rotate
  if (mstate->hover && mstate->lclick){
    mpos0 = {mstate->ndpos.x,mstate->ndpos.y,0.f};
    cpos0 = sphereProject(mstate->ndpos);
    crot0 = m_world;
    llock = true;
  } else if (llock) {
    if (mstate->ldown){
      vec3 cpos = sphereProject(mstate->ndpos);
      vec3 axis = cross(cpos0,cpos);
      if (length(axis) > 1e-10f){
        vec2 mpos = {mstate->ndpos.x-mpos0.x,mstate->ndpos.y-mpos0.y};
        float ang = 2.0f * length(mpos) * mousesens_rot;
        m_world = rotate(crot0,ang,axis);
        updateworld = true;
      }
    } else { 
      llock = false;
    }
  }

  if (updateworld)
    updateWorld();
  if (updateview)
    updateView();
  if (updateprojection)
    updateProjection();

  return updateview || updateworld || updateprojection;
}

void View::updateProjection(){
  if (isortho){
    float hw2 = tan(0.5f*zfov) * v_pos[2];
    m_projection = ortho(-hw2,hw2,-hw2,hw2,znear,1000.f);
  } else {
    m_projection = infinitePerspective(radians(zfov),1.0f,znear);
  }
  shader->setMat4("projection",value_ptr(m_projection));
}

void View::updateView(){
  m_view = lookAt(v_pos,v_pos+v_front,v_up);
  shader->setMat4("view",value_ptr(m_view));
}

void View::updateWorld(){
  shader->setMat4("world",value_ptr(m_world));
}

float View::getDepth(vec2 ndpos){
    float depth;
    glBindFramebuffer(GL_FRAMEBUFFER, FBO);
    glReadPixels(mstate->ndpos.x*FBO_tex_a,mstate->ndpos.y*FBO_tex_a, 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &depth);
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    return depth;
}

vec3 View::sphereProject(vec2 ndpos){
  vec2 xs = {(clamp(mstate->ndpos.x,0.f,1.f)-0.5f), (clamp(mstate->ndpos.y,0.f,1.f)-0.5f)};
  float a = 2.0f * fmin(length(xs),0.5f);
  float b = atan2f(xs.y,xs.x);
  return vec3(cosf(b) * sinf(a), sinf(b) * sinf(a), cosf(a));
}

