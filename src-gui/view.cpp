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

#include "imgui/font_glyphs.h"
#include "imgui/imgui.h"
#include "imgui/imgui_widgets.h"
#include "imgui/mouse.h"

#include "critic2.h"
#include "shapes.h"
#include "view.h"
#include "settings.h"

using namespace ImGui;
using namespace std;

// A linked list for all current views.
static list<View*> viewlist;

View *CreateView(char *title, Shader *shader, int iscene/*=0*/){
  View *aview = new View;

  // save the scene id and set the pointers to that scene
  aview->iscene = iscene;
  aview->title = title;
  aview->shader = shader;

  // preparation of the textures
  glGenTextures(1, &(aview->FBOtex));
  glGenRenderbuffers(1, &(aview->FBOdepth));
  glGenFramebuffers(1, &(aview->FBO));

  // create the texture
  glBindTexture(GL_TEXTURE_2D, aview->FBOtex);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, FBO_tex_a, FBO_tex_a, 0, GL_RGBA, GL_UNSIGNED_BYTE, 0);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);  
  glBindTexture(GL_TEXTURE_2D, 0);

  // render buffer
  glBindRenderbuffer(GL_RENDERBUFFER, aview->FBOdepth);
  glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT, FBO_tex_a, FBO_tex_a);
  glBindRenderbuffer(GL_RENDERBUFFER, 0);

  // frame buffer
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

  return aview;
}

void DrawAllViews(){
  for (auto iv : viewlist){
    iv->Draw();
  }
}

void View::Draw(){
  ImGuiContext *g = GetCurrentContext();

  // Variables for associated dialogs
  static bool drawprefs = false;

  // The dock window
  SetNextWindowSize(ImVec2(300.f,300.f),ImGuiSetCond_Once);
  PushStyleColor(ImGuiCol_WindowBg,ImVec4(bgrgb[0],bgrgb[1],bgrgb[2],bgrgb[3]));
  if (BeginDock(title)){
    // save cursor position at the top left
    ImVec2 cpos = GetCursorPos();

    // variables for the buttons
    const int nicon = 11;
    const int ihfill = 6;
    static bool hovered[nicon] = {};
    static bool held[nicon] = {};
    bool pressed[nicon] = {};
    bool usegray[nicon] = {};
    char *buttonchar[nicon] = {ICON_SM_ARROWS,ICON_SM_MOUSE_POINTER,
			       ICON_SM_COMPASS_ANGLE,ICON_SM_RULER,ICON_SM_PENCIL,
			       ICON_SM_ALIGNMENT,ICON_SM_COG,ICON_SM_TAG,
			       ICON_SM_FLOPPY_O,ICON_SM_QUESTION,ICON_SM_TIMES};
    char *buttontip[nicon] = {
      "Navigate\n\nLeft: Rotate\nRight: Translate\nMouse wheel: Zoom\nDouble click: Reset the view",
      "Select atoms and bonds",
      "Manipulate angles",
      "Measure",
      "Build",
      "Align",
      "Preferences",
      "Atom labels",
      "Save",
      "Query",
      "Close",
    };
    ImVec2 buttonsize = ImVec2(fonticon_size+1,fonticon_size+1);
    ImColor heldcolor = ImColor(0.7216f,0.5254,0.04314f);
    ImColor hovercolor = ImColor(0.8549f,0.6471f,0.1255f);
    ImColor graycolor = ImColor(0.5f,0.5f,0.5f);

    // Interactive part of the buttons; save variables for later
    PushStyleVar(ImGuiStyleVar_FramePadding, ImVec2(0.f,0.f));
    PushStyleVar(ImGuiStyleVar_ItemSpacing, ImVec2(0.f,0.f));
    PushFont(fonticon);
    for (int i = 0; i < nicon; i++){
      if (i == ihfill){
	SameLine(); Dummy(ImVec2(0.25f*fonticon_size,0.f));
	SameLine(); Text(ICON_SM_VBAR);
	SameLine(); Dummy(ImVec2(0.25f*fonticon_size,0.f));
      }
      PushID(buttonchar[i]);
      SameLine();
      pressed[i] = InvisibleButtonEx(buttonchar[i],buttonsize,&hovered[i],&held[i]); 
      AttachTooltip(buttontip[i],tooltip_delay,tooltip_maxwidth,fontdefault);
      PopID();
    }
    PopFont();
    PopStyleVar(2);

    // Overlay the image 
    bool hover = false;
    SetCursorPos(cpos);
    if (iscene > 0){
      ImageInteractive((void *) FBOtex,FBO_a/FBO_tex_a,&hover,&vrect);
    } else {
      hover = false;
      vrect = dock->window->Rect();
    }

    // process mouse events
    if (iscene > 0){
      c2::set_scene_pointers(iscene);
      if (processMouseEvents(hover) || updateTexSize())
	Update();
    }

    // 0: ICON_SM_ARROWS
    // 1: ICON_SM_MOUSE_POINTER,
    // 2: ICON_SM_COMPASS_ANGLE
    // 3: ICON_SM_RULER
    // 4: ICON_SM_PENCIL,
    // 5: ICON_SM_ALIGNMENT
    // 6: ICON_SM_COG
    // 7: ICON_SM_TAG,
    // 8: ICON_SM_FLOPPY_O
    // 9: ICON_SM_QUESTION ok
    // 10: ICON_SM_TIMES
    // xx: ICON_SM_VBAR

    // process button interactions
    usegray[0] = usegray[1] = usegray[2] = usegray[3] = usegray[4] = usegray[5] = true;
    if (pressed[0]) mousebehavior = MB_Navigation;
    if (pressed[1]) mousebehavior = MB_Pointer;
    if (pressed[2]) mousebehavior = MB_Angle;
    if (pressed[3]) mousebehavior = MB_Ruler;
    if (pressed[4]) mousebehavior = MB_Builder;
    if (pressed[5]) mousebehavior = MB_Alignment;
    if (pressed[6]) drawprefs = true;
    usegray[mousebehavior] = false;

    // Render the buttons on top of the image
    SetCursorPos(cpos);
    PushStyleColor(ImGuiCol_Button, ImColor(0, 0, 0, 0));
    PushStyleColor(ImGuiCol_ButtonHovered, ImColor(0, 0, 0, 0));
    PushStyleColor(ImGuiCol_ButtonActive, ImColor(0, 0, 0, 0));
    PushStyleVar(ImGuiStyleVar_FramePadding, ImVec2(0.f,0.f));
    PushStyleVar(ImGuiStyleVar_ItemSpacing, ImVec2(0.f,0.f));
    PushFont(fonticon);
    for (int i = 0; i < nicon; i++){
      if (i == ihfill){
	Dummy(ImVec2(0.25f*fonticon_size,0.f)); SameLine();
	Text(ICON_SM_VBAR); SameLine();
	Dummy(ImVec2(0.25f*fonticon_size,0.f)); SameLine();
      }
      PushID(buttonchar[i]);
      if (held[i])
	PushStyleColor(ImGuiCol_Text, heldcolor);
      else if (hovered[i])
	PushStyleColor(ImGuiCol_Text, hovercolor);
      else if (usegray[i])
	PushStyleColor(ImGuiCol_Text, graycolor);
      Button(buttonchar[i],buttonsize);
      if (held[i] || hovered[i] || usegray[i])
	PopStyleColor();
      SameLine();
      PopID();
    }
    PopFont();
    PopStyleVar(2);
    PopStyleColor(3);
  }
  dock = GetCurrentDock();
  EndDock();
  PopStyleColor();

  // Preferences dialog
  if (drawprefs){
    char prefstr[strlen(title) + 15];
    ImFormatString(prefstr,IM_ARRAYSIZE(prefstr),"Prefrences (%s)",title);
    if (BeginDock(prefstr,&drawprefs)){
      Text("Blah!");
    }
    EndDock();
  }
}

void View::Update(){

  if (icurtex < 0)
    updateTexSize();
  glBindFramebuffer(GL_FRAMEBUFFER, FBO);
  glViewport(0.,0.,FBO_a,FBO_a);

  glClearColor(bgrgb[0],bgrgb[1],bgrgb[2],bgrgb[3]);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  if (iswire)
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

  if (iscene > 0){
    // scene atoms
    for (int i=0;i<c2::nat;i++)
      drawSphere(make_vec3(c2::at[i].r),c2::at[i].rad,make_vec4(c2::at[i].rgb),isphres,false);

    // scene bonds
    for (int i=0;i<c2::nbond;i++){
      float rgb[4] = {0.5f,0.f,0.f,1.f};
      drawCylinder(make_vec3(c2::bond[i].r1),make_vec3(c2::bond[i].r2),c2::bond[i].rad,make_vec4(rgb),icylres,false);
    }

    // grid and Cartesian axes
    // int nmaxgrid = (int) (ceil(c2::scenerad * AUTOANG)+1e-5);
    // vec3 r1 = {}, r2 = {};
    // vec4 rgb = {1.f,1.f,1.f,1.f};
    // rgb.x = 1.f; rgb.y = 0.f; rgb.z = 0.f;
    // r1.y = 0.f; r1.x = -nmaxgrid/AUTOANG;
    // r2.y = 0.f; r2.x =  nmaxgrid/AUTOANG;
    // drawCylinder(r1,r2,radgrid,rgb,0,false);
    // rgb.x = 0.f; rgb.y = 1.f; rgb.z = 0.f;
    // r1.x = 0.f; r1.y = -nmaxgrid/AUTOANG;
    // r2.x = 0.f; r2.y =  nmaxgrid/AUTOANG;
    // drawCylinder(r1,r2,radgrid,rgb,0,false);
    // rgb.x = 0.f; rgb.y = 0.f; rgb.z = 1.f;
    // r1.y = 0.f; r1.z = -nmaxgrid/AUTOANG;
    // r2.y = 0.f; r2.z =  nmaxgrid/AUTOANG;
    // drawCylinder(r1,r2,radgrid,rgb,0,false);
    // for (int i = 1; i <= nmaxgrid; i++){
    //   r1.z = 0.f; r2.z = 0.f;
    //   rgb.x = 1.f; rgb.y = 1.f; rgb.z = 1.f; rgb.w = 1.f; 
    //   r1.x = i/AUTOANG; r1.y = -nmaxgrid/AUTOANG;
    //   r2.x = i/AUTOANG; r2.y =  nmaxgrid/AUTOANG;
    //   drawCylinder(r1,r2,radgrid,rgb);
    //   r1.x = -i/AUTOANG; r1.y = -nmaxgrid/AUTOANG;
    //   r2.x = -i/AUTOANG; r2.y =  nmaxgrid/AUTOANG;
    //   drawCylinder(r1,r2,radgrid,rgb);
    //   r1.y =  i/AUTOANG; r1.x = -nmaxgrid/AUTOANG;
    //   r2.y =  i/AUTOANG; r2.x =  nmaxgrid/AUTOANG;
    //   drawCylinder(r1,r2,radgrid,rgb);
    //   r1.y = -i/AUTOANG; r1.x = -nmaxgrid/AUTOANG;
    //   r2.y = -i/AUTOANG; r2.x =  nmaxgrid/AUTOANG;
    //   drawCylinder(r1,r2,radgrid,rgb);
    // }

    // rotation sphere
    // if (llock){
    //   vec3 v0;
    //   float rad0 = 0.25f * v_pos[2];
    //   float rad1;
    //   vec4 rgb;
    //   vec3 uvec = vec3(0.f,0.f,1.f);

    //   v0 = inverse(mat3(m_world)) * uvec;
    //   v0 = v0 * rad0;
    //   rgb = {0.6f,0.2f,0.8f,0.8f};
    //   rad1 = 0.05f * rad0;
    //   drawSphere(v0,rad1,rgb,3,false);

    //   v0 = inverse(mat3(crot0_l)) * uvec;
    //   v0 = v0 * rad0;
    //   rgb = {1.0f,0.8f,0.0f,0.8f};
    //   rad1 = 0.05f * rad0;
    //   drawSphere(v0,rad1,rgb,3,false);

    //   v0 = vec3(0.f,0.f,0.f);
    //   rgb = {1.0f,1.0f,1.0f,0.4f};
    //   drawSphere(v0,rad0,rgb,3,true);
    // }
  }

  glBindFramebuffer(GL_FRAMEBUFFER, 0);
  if (iswire)
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
}

void View::Delete(){
  for (auto it = viewlist.begin(); it != viewlist.end(); it++) {
    if (*it == this){
      glDeleteTextures(1, &FBOtex);
      glDeleteRenderbuffers(1, &FBOdepth);
      glDeleteFramebuffers(1, &FBO);
      viewlist.erase(it);
      break;
    }
  }
}

bool View::processMouseEvents(bool hover){
  if (mousebehavior == MB_Navigation)
    return Navigate(hover);
  else if (mousebehavior == MB_Pointer)
    return false;
  else if (mousebehavior == MB_Angle)
    return false;
  else if (mousebehavior == MB_Ruler)
    return false;
  else if (mousebehavior == MB_Builder)
    return false;
  else if (mousebehavior == MB_Alignment)
    return false;
}

bool View::Navigate(bool hover){
  const float eps = 1e-8;
  const vec4 viewport = {0.f,0.f,FBO_a,FBO_a};
  bool updateview = false, updateworld = false, updateprojection = false;
  bool updatenone = false;

  // calculate the texture coordinates
  vec2 texpos = mstate.pos;
  pos_to_texpos(texpos);

  // mouse scroll = zoom
  if (hover && abs(mstate.scroll) > eps && !rlock && !llock){
    float ratio = fmin(mousesens_zoom * mstate.scroll,0.5f);
    v_pos = v_pos - ratio * v_pos;
    if (length(v_pos) < min_zoom)
      v_pos = v_pos / length(v_pos) * min_zoom;
    if (isortho)
      updateprojection = true;
    updateview = true;
  }

  // drag
  if (hover && mstate.rclick && !llock){ 
    float depth = texpos_viewdepth(texpos);
    if (depth < 1.0){
      mpos0_r = {texpos.x,texpos.y,depth};
    }else{
      mpos0_r = {texpos.x,texpos.y,0.f};
      view_to_texpos({0.f,0.f,0.f},&mpos0_r.z);
    }
    cpos0_r = {v_pos[0],v_pos[1],0.f};
    rlock = true;
  } else if (rlock) {
    if (mstate.rdown){
      vec3 vnew = texpos_to_view(texpos,mpos0_r.z);
      vec3 vold = texpos_to_view(vec2(mpos0_r),mpos0_r.z);
      v_pos.x = cpos0_r.x - (vnew.x - vold.x);
      v_pos.y = cpos0_r.y - (vnew.y - vold.y);
      updateview = true;
    } else {
      rlock = false;
    }
    // updatenone = true; // only if we draw some guiding element
  }

  // rotate
  if (hover && mstate.lclick && !rlock){
    mpos0_l = {texpos.x, texpos.y, 0.f};
    view_to_texpos({0.f,0.f,0.f},&mpos0_l.z);
    cpos0_l = texpos_to_view(texpos,mpos0_l.z);
    crot0_l = m_world;
    llock = true;
  } else if (llock) {
    if (mstate.ldown){
      vec3 cpos1 = texpos_to_view(texpos,mpos0_l.z);
      vec3 axis = cross(vec3(0.f,0.f,1.f),cpos1-cpos0_l);
      float lax = length(axis);
      if (lax > 1e-10f){
	axis = inverse(mat3(crot0_l)) * normalize(axis);
	vec2 mpos = {texpos.x-mpos0_l.x, texpos.y-mpos0_l.y};
	float ang = 2.0f * length(mpos) * mousesens_rot / FBO_a;
	m_world = rotate(crot0_l,ang,axis);
	updateworld = true;
      }
    } else { 
      llock = false;
    }
    // updatenone = true; // only if we draw some guiding element
  }

  // double click
  if (hover && mstate.ldclick){
    if (iscene > 0)
      v_pos = {0.f,0.f,4.f*c2::scenerad};
    else
      v_pos = {0.f,0.f,10.f};
    v_front = {0.f,0.f,-1.f};
    v_up    = {0.f,1.f,0.f};
    m_world = mat4(1.0f);
    crot0_l = mat4(1.0f);
    llock = false;
    rlock = false;
    if (isortho)
      updateprojection = true;
    updateview = true;
    updateworld = true;
  }

  if (updateworld)
    updateWorld();
  if (updateview)
    updateView();
  if (updateprojection)
    updateProjection();

  return updateworld || updateview || updateprojection || updatenone;
}


bool View::updateTexSize(){
  float amax = ((dock && dock->size.x > 0.f && dock->size.y > 0.f)? fmax(dock->size.x,dock->size.y) : 200.f);
  if (FBO_a != amax){
    FBO_a = amax;
    return true;
  }
  return false;
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

vec3 View::cam_world_coords(){
  vec4 pos4 = inverse(m_world) * vec4(v_pos,1.0f);
  return vec3(pos4.x/pos4.w,pos4.y/pos4.w,pos4.z/pos4.w);
}

vec3 View::cam_view_coords(){
  return v_pos;
}

void View::pos_to_ntexpos(vec2 &pos){
  float x = vrect.Max.x - vrect.Min.x;
  float y = vrect.Max.y - vrect.Min.y;
  float xratio = 2.f * x / fmax(x,y);
  float yratio = 2.f * y / fmax(x,y);
  
  pos.x = ((pos.x - vrect.Min.x) / x - 0.5f) * xratio;
  pos.y = (0.5f - (pos.y - vrect.Min.y) / y) * yratio;
}

void View::ntexpos_to_pos(vec2 &pos){
  float x = vrect.Max.x - vrect.Min.x;
  float y = vrect.Max.y - vrect.Min.y;
  float xratio1 = 0.5f * fmax(x,y) / x;
  float yratio1 = 0.5f * fmax(x,y) / y;
  
  pos.x = vrect.Min.x + x * (0.5f + xratio1 * pos.x);
  pos.y = vrect.Min.y + y * (0.5f - yratio1 * pos.y);
}

void View::pos_to_texpos(vec2 &pos){
  pos_to_ntexpos(pos);
  ntexpos_to_texpos(pos);
}

void View::texpos_to_pos(vec2 &pos){
  texpos_to_ntexpos(pos);
  ntexpos_to_pos(pos);
}

void View::texpos_to_ntexpos(vec2 &pos){
  pos = (pos / FBO_a) * 2.f - 1.f;
}

void View::ntexpos_to_texpos(vec2 &pos){
  pos = (0.5f * pos + 0.5f) * FBO_a;
}

vec2 View::world_to_texpos(vec3 pos){
  const vec4 viewport = {0.f,0.f,FBO_a,FBO_a};
  vec3 pos3 = project(pos,m_view * m_world,m_projection,viewport);
  return vec2(pos3);
}

// dist=0, znear; dist<0, scene origin plane; dist>0, distance from camera
vec3 View::texpos_to_world(vec2 pos, float dist/*=-1.f*/){
  const vec4 viewport = {0.f,0.f,FBO_a,FBO_a};
  vec3 wpos = {};
  if (dist < 0.f){
    // Set the point on the plane parallel to the z-plane that passes through
    // the origin of the scene
    wpos = project(wpos,m_view * m_world,m_projection,viewport);
    wpos.x = pos.x; wpos.y = pos.y;
    wpos = unProject(wpos,m_view * m_world,m_projection,viewport);
  } else if (dist == 0.f) {
    // Set the point on the near plane
    wpos.x = pos.x; wpos.y = pos.y;
    wpos = unProject(wpos,m_view * m_world,m_projection,viewport);
  } else {
    // Set the point at a distance dist from the camera
    vec3 origpos = texpos_to_world(pos,-1.f);
    vec3 nearpos = texpos_to_world(pos,0.f);
    vec3 dir = normalize(origpos - nearpos);
    wpos = nearpos + fmax(dist - znear,0.f) * dir;
  }
  return wpos;
}

vec2 View::world_to_ntexpos(vec3 pos){
  vec2 pos2 = world_to_texpos(pos);
  texpos_to_ntexpos(pos2);
  return pos2;
}

// dist=0, znear; dist<0, scene origin plane; dist>0, distance from camera
vec3 View::ntexpos_to_world(vec2 pos, float dist/*=-1.f*/){
  ntexpos_to_texpos(pos);
  return texpos_to_world(pos,dist);
}

vec2 View::view_to_texpos(vec3 pos, float *depth){
  const vec4 viewport = {0.f,0.f,FBO_a,FBO_a};
  vec3 pos3 = project(pos,m_view,m_projection,viewport);
  *depth = pos3.z;
  return vec2(pos3);
}

vec3 View::texpos_to_view(vec2 pos, float depth){
  const vec4 viewport = {0.f,0.f,FBO_a,FBO_a};
  vec3 wpos = {pos.x,pos.y,depth};
  wpos = unProject(wpos,m_view,m_projection,viewport);
  return wpos;
}

float View::texpos_viewdepth(vec2 texpos){
    float depth;
    glBindFramebuffer(GL_FRAMEBUFFER, FBO);
    glReadPixels(texpos.x,texpos.y, 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &depth);
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    return depth;
}

void View::drawSphere(vec3 r0, float rad, vec4 rgb, int res, bool blend){
  if (blend){
    glEnable(GL_BLEND);
    glEnable(GL_CULL_FACE);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glDepthMask(0);
  }
  glBindVertexArray(sphVAO[res]);

  mat4 m_model = mat4(1.0f);
  m_model = translate(m_model,r0);
  m_model = scale(m_model,vec3(rad,rad,rad));
  mat3 m_normrot = transpose(inverse(mat3(m_view) * mat3(m_world) * mat3(m_model)));

  shader->setVec4("vColor",value_ptr(rgb));
  shader->setMat4("model",value_ptr(m_model));
  shader->setMat3("normrot",value_ptr(m_normrot));
  glDrawElements(GL_TRIANGLES, 3*sphnel[res], GL_UNSIGNED_INT, 0);
  if (blend){
    glDisable(GL_BLEND);
    glDisable(GL_CULL_FACE);
    glDepthMask(1);
  }
}

void View::drawCylinder(vec3 r1, vec3 r2, float rad, vec4 rgb, int res, bool blend){
  if (blend){
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glDepthMask(0);
  }
  glBindVertexArray(cylVAO[res]);
  vec3 xmid = 0.5f * (r1 + r2);
  vec3 xdif = r2 - r1;
  float blen = length(xdif);
  
  mat4 m_model = mat4(1.0f);
  m_model = translate(m_model,xmid);
  m_model = m_model * orientation(xdif/blen,vec3(0.f,0.f,1.f));
  m_model = scale(m_model,vec3(rad,rad,blen));
  mat3 m_normrot = transpose(inverse(mat3(m_view) * mat3(m_world) * mat3(m_model)));
  shader->setMat3("normrot",value_ptr(m_normrot));
  shader->setVec4("vColor",value_ptr(rgb));
  shader->setMat4("model",value_ptr(m_model));
  glDrawElements(GL_TRIANGLES, 3*cylnel[icylres], GL_UNSIGNED_INT, 0);
  if (blend){
    glDisable(GL_BLEND);
    glDepthMask(1);
  }
}

