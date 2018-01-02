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
#include <algorithm>
#include <sstream>

#include "imgui/icon_glyphs.h"
#include "imgui/imgui.h"
#include "imgui/imgui_widgets.h"

#include "critic2.h"
#include "shapes.h"
#include "view.h"
#include "settings.h"
#include "keybinding.h"

using namespace ImGui;

// A linked list for all current views.
static list<View*> viewlist;

// Create a tooltip label
static string view_tooltip_label(int id){
  stringstream str;
  
  switch(id){
  case 0:
    str << "Navigate\n\n" 
	<< BindKeyName(BIND_NAV_ROTATE) << ": Rotate\n"
        << BindKeyName(BIND_NAV_TRANSLATE) << ": Translate\n"
	<< BindKeyName(BIND_NAV_ZOOM) << ": Zoom\n"
	<< BindKeyName(BIND_NAV_RESET) << ": Reset the view\n"
	<< BindKeyName(BIND_VIEW_ALIGN_A_AXIS) << "/"
	<< BindKeyName(BIND_VIEW_ALIGN_B_AXIS) << "/"
	<< BindKeyName(BIND_VIEW_ALIGN_C_AXIS) << ": Align to crystallographic axes\n"
	<< BindKeyName(BIND_VIEW_ALIGN_X_AXIS) << "/"
	<< BindKeyName(BIND_VIEW_ALIGN_Y_AXIS) << "/"
	<< BindKeyName(BIND_VIEW_ALIGN_Z_AXIS) << ": Align to Cartesian axes";
    break;
  case 1:  str << "Select atoms and bonds"; break;
  case 2:  str << "Manipulate angles"; break;
  case 3:  str << "Measure"; break;
  case 4:  str << "Build"; break;
  case 5:  str << "Align"; break;
  case 6:  str << "Query"; break;
  case 7:  str << "Preferences"; break;
  };
  return str.str();
}

void View::SetDefaults(){
  iswire = view_wireframe;
  isortho = view_orthogonal;
  zfov = view_fov;
  resetd = view_resetdistance; 
  for (int i=0; i<4; i++)
    bgrgb[i] = view_bgrgb[i];
  show_atoms = view_show_atoms;
  scale_atoms = view_scale_atoms;
  isphres = view_isphres;
  show_bonds = view_show_bonds;
  scale_bonds = view_scale_bonds;
  icylres = view_icylres;
  ncell[0] = ncell[1] = ncell[2] = 1;
  isucell = false;
  isborder = false;
  ismotif = false;
  ismolcell = false;
  if (iscene > 0 && !c2::ismolecule){
    isucell = true;
    isborder = true;
    ismotif = true;
  }
}

void View::Draw(){
  ImGuiContext *g = GetCurrentContext();
  ImGuiIO& io = GetIO();

  // Variables for associated dialogs
  bool drawprefs = false;
  PushStyleColor(ImGuiCol_WindowBg,ImVec4(bgrgb[0],bgrgb[1],bgrgb[2],bgrgb[3]));
  if (BeginDock(title)){
    // save cursor position at the top left
    ImVec2 cpos = GetCursorPos();

    // variables for the buttons
    const int nicon = 8;
    const int ihfill = 7;
    static bool hovered[nicon] = {};
    static bool held[nicon] = {};
    bool pressed[nicon] = {};
    bool usegray[nicon] = {};
    char *buttonchar[nicon] = {ICON_SM_ARROWS,ICON_SM_MOUSE_POINTER,
			       ICON_SM_COMPASS_ANGLE,ICON_SM_RULER,ICON_SM_PENCIL,
			       ICON_SM_ALIGNMENT,ICON_SM_QUESTION,
			       ICON_SM_COG};
    ImVec2 buttonsize = ImVec2(ImGuiStyleUI.FontSizeIcon + 1,ImGuiStyleUI.FontSizeIcon + 1);
    ImVec4 color = ImGuiStyleUI.Colors[ImGuiColUI_ViewIcon];
    ImVec4 heldcolor = ImGuiStyleUI.Colors[ImGuiColUI_ViewIconActive];
    ImVec4 hovercolor = ImGuiStyleUI.Colors[ImGuiColUI_ViewIconHovered];
    ImVec4 graycolor = ImGuiStyleUI.Colors[ImGuiColUI_ViewIconInactive];
    float fpad = 0.f;

    // Interactive part of the buttons; save variables for later
    PushStyleVar(ImGuiStyleVar_FramePadding, ImVec2(0.f,0.f));
    PushStyleVar(ImGuiStyleVar_ItemSpacing, ImVec2(0.f,0.f));
    PushFont(fonticon);
    for (int i = 0; i < nicon; i++){
      if (i == ihfill){
	fpad = GetContentRegionAvailWidth()-ImGuiStyleUI.FontSizeIcon-1;
	SameLine(0, fpad);
      }
      PushID(buttonchar[i]);
      pressed[i] = InvisibleButtonEx(buttonchar[i],buttonsize,&hovered[i],&held[i]); 
      if (ImGuiStyleUI.TooltipEnabled)
	AttachTooltip(view_tooltip_label(i).c_str(),ImGuiStyleUI.TooltipDelay,ImGuiStyleUI.TooltipMaxwidth,fontdefault);
      SameLine();
      PopID();
    }
    PopFont();
    PopStyleVar(2);

    // Set the scene pointers
    if (iscene > 0)
      c2::set_scene_pointers(iscene);

    // Maybe we need to change the texture sizes
    updatescene = updateTexSize() || updatescene;

    // Update the scene
    if (updatescene){
      Update();
      updatescene = false;
    }

    // Overlay the image
    bool hover = false;
    SetCursorPos(cpos);
    if (iscene > 0){
      ImageInteractive((void *) FBOtex0,FBO_a/FBO_atex,&hover,&vrect);
    } else {
      hover = false;
      vrect = dock->window->Rect();
    }

    // Process mouse events
    if (iscene > 0)
      updatescene = processMouseEvents(hover) || updatescene;

    // process button interactions
    usegray[0] = usegray[1] = usegray[2] = usegray[3] = usegray[4] = usegray[5] = usegray[6] = true;
    if (pressed[0]) mousebehavior = MB_Navigation;
    if (pressed[1]) mousebehavior = MB_Pointer;
    if (pressed[2]) mousebehavior = MB_Angle;
    if (pressed[3]) mousebehavior = MB_Ruler;
    if (pressed[4]) mousebehavior = MB_Builder;
    if (pressed[5]) mousebehavior = MB_Alignment;
    if (pressed[6]) mousebehavior = MB_Query;
    if (pressed[7]) drawprefs = true;
    usegray[mousebehavior] = false;

    // Render the buttons on top of the image
    SetCursorPos(cpos);
    PushStyleColor(ImGuiCol_Button, ImVec4(0, 0, 0, 0));
    PushStyleColor(ImGuiCol_ButtonHovered, ImVec4(0, 0, 0, 0));
    PushStyleColor(ImGuiCol_ButtonActive, ImVec4(0, 0, 0, 0));
    PushStyleVar(ImGuiStyleVar_FramePadding, ImVec2(0.f,0.f));
    PushStyleVar(ImGuiStyleVar_ItemSpacing, ImVec2(0.f,0.f));
    PushFont(fonticon);
    for (int i = 0; i < nicon; i++){
      if (i == ihfill)
	SameLine(0, fpad);
      PushID(buttonchar[i]);
      if (held[i])
	PushStyleColor(ImGuiCol_Text, heldcolor);
      else if (hovered[i])
	PushStyleColor(ImGuiCol_Text, hovercolor);
      else if (usegray[i])
	PushStyleColor(ImGuiCol_Text, graycolor);
      else
	PushStyleColor(ImGuiCol_Text, color);
      Button(buttonchar[i],buttonsize);
      PopStyleColor();
      SameLine();
      PopID();
    }
    PopFont();
    PopStyleVar(2);
    PopStyleColor(3);

    // Set a reasonable content size for the window
    NewLine();
    Dummy(ImVec2(0.4f*io.DisplaySize.x,0.4f*io.DisplaySize.y));
  }
  dock = GetCurrentDock();
  EndDock();
  PopStyleColor();

  // Preferences dialog
  if (drawprefs)
    ImGui::OpenPopup("prefview");
  if (BeginPopup("prefview")){
    bool changed = false;
    float itemwidth = 5.f * g->FontSize;
    updatescene |= Checkbox("Unit cell", &isucell);
    if (iscene > 0 && c2::ismolecule)
      updatescene |= Checkbox("Molecular cell", &ismolcell);
    PushItemWidth(itemwidth);
    if (iscene > 0){
      if (!c2::ismolecule){
	AlignTextToFramePadding();
	Text("Number of cells:");
	SameLine();
	if (Button("Reset")){
	  updatescene = true;
	  ncell[0] = ncell[1] = ncell[2] = 1;
	}
	Indent();
	AlignTextToFramePadding();
	Text("a:"); SameLine(0.f,0.f);
	updatescene |= InputInt("##aaxis", &ncell[0]);
	SameLine();
	Text("b:"); SameLine(0.f,0.f);
	updatescene |= InputInt("##baxis", &ncell[1]);
	SameLine();
	Text("c:"); SameLine(0.f,0.f);
	updatescene |= InputInt("##caxis", &ncell[2]);
	for (int i=0; i<3; i++)
	  ncell[i] = std::max(1,ncell[i]);
	Unindent();
      }

      Text("View/Axis Alignment:");
      Indent();
      int ialign = 0;
      const float aside = g->FontSize + g->Style.FramePadding.y * 2.0f;
      const ImVec2 siz = {aside,aside};
      if (!c2::ismolecule){
	if (Button("a",siz)) ialign = 1; SameLine();
	if (Button("b",siz)) ialign = 2; SameLine();
	if (Button("c",siz)) ialign = 3; SameLine();
      }
      if (Button("x",siz)) ialign = -1; SameLine();
      if (Button("y",siz)) ialign = -2; SameLine();
      if (Button("z",siz)) ialign = -3;
      Unindent();
      if (ialign != 0){
	alignViewAxis(ialign);
	updatescene = true;
      }

      if (!c2::ismolecule){
	updatescene |= Checkbox("Crystal packing", &isborder);
	updatescene |= Checkbox("Molecular motif", &ismotif);
      }
    }
    Separator();
    changed |= Checkbox("Wireframe rendering", &view_wireframe);
    changed |= Checkbox("Orthgonal projection", &view_orthogonal);
    if (!view_orthogonal){
      Indent();
      changed |= DragFloat("Field of view (degrees)", &view_fov, 2.5f, 0.0f, 180.0f, "%.1f", 1.0f); 
      Unindent();
    }
    changed |= DragFloat("Reset distance (scene radius)", &view_resetdistance, 0.02f, 0.0f, 10.f, "%.2f", 1.0f); 
    changed |= ColorEdit4("Background color", view_bgrgb, coloreditflags);
    changed |= Checkbox("Show atoms", &view_show_atoms);
    if (view_show_atoms){
      Indent();
      changed |= SliderInt("Atom resolution", &view_isphres, 0, nmaxsph-1); 
      changed |= DragFloat("Atom size", &view_scale_atoms, 0.01f, 0.0f, 5.f, "%.2f", 1.0f);
      Unindent();
    }
    changed |= Checkbox("Show bonds", &view_show_bonds);
    if (view_show_bonds){
      Indent();
      changed |= SliderInt("Bond resolution", &view_icylres, 0, nmaxcyl-1); 
      changed |= DragFloat("Bond size", &view_scale_bonds, 0.01f, 0.0f, 5.f, "%.2f", 1.0f);
      Unindent();
    }
    PopItemWidth();
    
    // For now...
    if (changed)
      SetDefaultAllViews();
    if (IsMouseClicked(1) || IsBindEvent(BIND_CLOSE_LAST_DIALOG,false))
      CloseCurrentPopup();
    EndPopup();
  }
}

void View::Update(){
  glEnable(GL_MULTISAMPLE);
  glBindFramebuffer(GL_FRAMEBUFFER, FBO);
  glViewport(0.,0.,FBO_a,FBO_a);

  glClearColor(bgrgb[0],bgrgb[1],bgrgb[2],bgrgb[3]);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  if (iswire)
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

  if (iscene > 0 && c2::isinit == 2){
    shader->setInt("uselighting",1);
    const float rthr = 0.01f;

    glm::vec3 v0 = {0.f,0.f,0.f};
    glm::vec3 vx = {c2::avec[0][0],c2::avec[0][1],c2::avec[0][2]};
    glm::vec3 vy = {c2::avec[1][0],c2::avec[1][1],c2::avec[1][2]};
    glm::vec3 vz = {c2::avec[2][0],c2::avec[2][1],c2::avec[2][2]};

    // prepare the lattice vector limits for each atom
    int imin[c2::nat][3], imax[c2::nat][3];
    int iminf[c2::nmol][3] = {}, imaxf[c2::nmol][3] = {};
    for (int i=0; i<c2::nat; i++){
      for (int j=0; j<3; j++){
	iminf[c2::at[i].ifrag][j] = std::min(iminf[c2::at[i].ifrag][j],-c2::at[i].flvec[j]);
	imaxf[c2::at[i].ifrag][j] = std::max(imaxf[c2::at[i].ifrag][j],-c2::at[i].flvec[j]);
      }
    }
    for (int i=0;i<c2::nat;i++){
      for (int j=0; j<3; j++){
	imin[i][j] = 0;
	imax[i][j] = ncell[j];
	if (isborder){
	  if (ismotif && c2::moldiscrete[c2::at[i].ifrag]){
	    imin[i][j] = iminf[c2::at[i].ifrag][j];
	    imax[i][j] = ncell[j] + imaxf[c2::at[i].ifrag][j];
	  } else {
	    if (c2::at[i].x[j] < rthr)
	      imax[i][j] = ncell[j] + 1;
	    else if (c2::at[i].x[j] > 1.f - rthr)
	      imin[i][j] = -1;
	  }
	}
	if (ismotif && c2::moldiscrete[c2::at[i].ifrag]){
	  imin[i][j] = imin[i][j] + c2::at[i].flvec[j];
	  imax[i][j] = imax[i][j] + c2::at[i].flvec[j];
	}
      }
    }

    // calculate the center of the scene (right now: center of the middle cell)
    glm::vec3 center = v0 + 0.5f * ((float) ncell[0] * vx + (float) ncell[1] * vy + (float) ncell[2] * vz);
    v0 = v0 - center;

    // scene atoms and bonds
    int (*idcon_)[c2::mncon] = (int (*)[c2::mncon]) c2::idcon;
    int (*lcon_)[c2::mncon][3] = (int (*)[c2::mncon][3]) c2::lcon;

    for (int i=0;i<c2::nat;i++){
      glm::vec3 r0 = glm::make_vec3(c2::at[i].r) - center;
      glm::vec4 rgb = glm::make_vec4(c2::at[i].rgb);

      for (int ix=imin[i][0]; ix<imax[i][0]; ix++){
	for (int iy=imin[i][1]; iy<imax[i][1]; iy++){
	  for (int iz=imin[i][2]; iz<imax[i][2]; iz++){
	    glm::vec3 x0 = r0 + (float) ix * vx + (float) iy * vy + (float) iz * vz;
	    if (show_atoms)
	      drawSphere(x0,scale_atoms * c2::at[i].rad,rgb,isphres,false);

	    if (!show_bonds) continue;
	    for (int j=0;j<c2::at[i].ncon;j++){
	      int ineigh = idcon_[i][j];
	      int ixn = ix + lcon_[i][j][0];
	      int iyn = iy + lcon_[i][j][1];
	      int izn = iz + lcon_[i][j][2];
	      if (ineigh < i || ineigh == i && (ixn < ix || ixn == ix && (iyn < iy || iyn == iy && izn < iz))) continue;
	      if (ixn < imin[ineigh][0] || ixn >= imax[ineigh][0] || iyn < imin[ineigh][1] || iyn >= imax[ineigh][1] ||
		  izn < imin[ineigh][2] || izn >= imax[ineigh][2]) continue;

	      glm::vec3 x1 = glm::make_vec3(c2::at[ineigh].r) - center + (float) ixn * vx + (float) iyn * vy + (float) izn * vz;
	      glm::vec3 xmid = x0 + 0.5f * (x1 - x0);
	      glm::vec4 rgbn = glm::make_vec4(c2::at[ineigh].rgb);

	      const float rad = 0.2f;
	      drawCylinder(x0,xmid,scale_bonds * rad,rgb,icylres,false);
	      drawCylinder(xmid,x1,scale_bonds * rad,rgbn,icylres,false);
	    }
	  }
	}
      }
    }

    // unit cell
    shader->setInt("uselighting",0);
    if (isucell)
      drawUnitCell(v0,vx,vy,vz,true);
    if (ismolcell){
      glm::vec3 v0_ = v0 + c2::molborder[0] * vx + c2::molborder[1] * vy + c2::molborder[2] * vz;
      glm::vec3 vx_ = vx * (1.f - 2.f * c2::molborder[0]);
      glm::vec3 vy_ = vy * (1.f - 2.f * c2::molborder[1]);
      glm::vec3 vz_ = vz * (1.f - 2.f * c2::molborder[2]);
      drawUnitCell(v0_,vx_,vy_,vz_,false);
    }

    // // the scenerad spehre, for testing
    // vec3 v0 = vec3(0.f,0.f,0.f);
    // vec4 rgb = {1.0f,1.0f,1.0f,0.4f};
    // drawSphere(v0,c2::scenerad,rgb,3,true);
  }

  glDisable(GL_MULTISAMPLE);
  glBindFramebuffer(GL_READ_FRAMEBUFFER, FBO);
  glBindFramebuffer(GL_DRAW_FRAMEBUFFER, FBO0);
  glBlitFramebuffer(0, 0, FBO_atex, FBO_atex, 0, 0, FBO_atex, FBO_atex, GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT, GL_NEAREST);
  glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
  glBindFramebuffer(GL_READ_FRAMEBUFFER, 0);

  glBindFramebuffer(GL_FRAMEBUFFER, 0);
  if (iswire)
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

}

void View::Delete(){
  for (auto it = viewlist.begin(); it != viewlist.end(); it++) {
    if (*it == this){
      deleteTex();
      viewlist.erase(it);
      break;
    }
  }
}

bool View::processMouseEvents(bool hover){
  bool outb = false;

  // common to all modes
  if (hover && IsBindEvent(BIND_VIEW_ALIGN_A_AXIS,false))
    outb |= alignViewAxis(1);
  else if (hover && IsBindEvent(BIND_VIEW_ALIGN_B_AXIS,false))
    outb |= alignViewAxis(2);
  else if (hover && IsBindEvent(BIND_VIEW_ALIGN_C_AXIS,false))
    outb |= alignViewAxis(3);
  else if (hover && IsBindEvent(BIND_VIEW_ALIGN_X_AXIS,false))
    outb |= alignViewAxis(-1);
  else if (hover && IsBindEvent(BIND_VIEW_ALIGN_Y_AXIS,false))
    outb |= alignViewAxis(-2);
  else if (hover && IsBindEvent(BIND_VIEW_ALIGN_Z_AXIS,false))
    outb |= alignViewAxis(-3);

  // process keybindings specific to each mode
  if (mousebehavior == MB_Navigation)
    outb |= Navigate(hover);
  else if (mousebehavior == MB_Pointer)
    outb |= false;
  else if (mousebehavior == MB_Angle)
    outb |= false;
  else if (mousebehavior == MB_Ruler)
    outb |= false;
  else if (mousebehavior == MB_Builder)
    outb |= false;
  else if (mousebehavior == MB_Query)
    outb |= false;
  else if (mousebehavior == MB_Alignment)
    outb |= false;
  return outb;
}

bool View::Navigate(bool hover){
  const float eps = 1e-8;
  const glm::vec4 viewport = {0.f,0.f,FBO_a,FBO_a};
  bool updateview = false, updateworld = false, updateprojection = false;
  bool updatenone = false;

  // calculate the texture coordinates
  ImGuiContext *g = GetCurrentContext();
  glm::vec2 mousepos = {g->IO.MousePos.x,g->IO.MousePos.y};
  glm::vec2 texpos = mousepos;
  pos_to_texpos(texpos);

  // Zoom the view. There are two behaviors: mouse scroll and hold key
  // and translate the mouse.
  float ratio = 0.f;
  if (hover && IsBindEvent(BIND_NAV_ZOOM,false) && !rlock && !llock){
    if (keybind[BIND_NAV_ZOOM] == GLFW_MOUSE_SCROLL){
      ratio = mousesens_zoom0 * view_mousesens_zoom * g->IO.MouseWheel;
    } else {
      mpos0_s = mousepos.y;
      slock = true;
    }
  } else if (slock) {
    if (IsBindEvent(BIND_NAV_ZOOM,true)){
      // 10/a to make it adimensional
      ratio = mousesens_zoom0 * view_mousesens_zoom * (mpos0_s-mousepos.y) * (10.f / FBO_a); 
      mpos0_s = mousepos.y;
    } else {
      slock = false;
    }
    // updatenone = true; // only if we draw some guiding element
  }
  if (ratio != 0.f){
    ratio = glm::clamp(ratio,-0.999f,0.999f);
    v_pos = v_pos - ratio * v_pos;
    if (glm::length(v_pos) < min_zoom)
      v_pos = v_pos / glm::length(v_pos) * min_zoom;
    if (glm::length(v_pos) > max_zoom * c2::scenerad)
      v_pos = v_pos / glm::length(v_pos) * (max_zoom * c2::scenerad);
    if (isortho)
      updateprojection = true;
    updateview = true;
  }

  // drag
  if (hover && IsBindEvent(BIND_NAV_TRANSLATE,false) && !llock && !slock){ 
    float depth = texpos_viewdepth(texpos);
    if (depth < 1.0){
      mpos0_r = {texpos.x,texpos.y,depth};
    }else{
      mpos0_r = {texpos.x,texpos.y,0.f};
      view_to_texpos({0.f,0.f,0.f},&mpos0_r.z);
    }
    cpos0_r = {v_pos[0],v_pos[1],0.f};
    rlock = true;
    mposlast = mousepos;
  } else if (rlock) {
    SetMouseCursor(ImGuiMouseCursor_Move);
    if (IsBindEvent(BIND_NAV_TRANSLATE,true)){
      if (mousepos.x != mposlast.x || mousepos.y != mposlast.y){
	glm::vec3 vnew = texpos_to_view(texpos,mpos0_r.z);
	glm::vec3 vold = texpos_to_view(glm::vec2(mpos0_r),mpos0_r.z);
	v_pos.x = cpos0_r.x - (vnew.x - vold.x);
	v_pos.y = cpos0_r.y - (vnew.y - vold.y);
	mposlast = mousepos;
	updateview = true;
      }
    } else {
      rlock = false;
    }
    // updatenone = true; // only if we draw some guiding element
  }

  // rotate
  if (hover && IsBindEvent(BIND_NAV_ROTATE,false) && !rlock && !slock){
    mpos0_l = {texpos.x, texpos.y, 0.f};
    view_to_texpos({0.f,0.f,0.f},&mpos0_l.z);
    cpos0_l = texpos_to_view(texpos,mpos0_l.z);
    llock = true;
  } else if (llock) {
    SetMouseCursor(ImGuiMouseCursor_Move);
    if (IsBindEvent(BIND_NAV_ROTATE,true)){
      if (texpos.x != mpos0_l.x || texpos.y != mpos0_l.y){
	glm::vec3 cpos1 = texpos_to_view(texpos,mpos0_l.z);
	glm::vec3 axis = glm::cross(glm::vec3(0.f,0.f,1.f),cpos1-cpos0_l);
	float lax = glm::length(axis);
	if (lax > 1e-10f){
	  axis = glm::inverse(glm::mat3(m_world)) * glm::normalize(axis);
	  glm::vec2 mpos = {texpos.x-mpos0_l.x, texpos.y-mpos0_l.y};
	  float ang = 2.0f * glm::length(mpos) * mousesens_rot0 * view_mousesens_rot / FBO_a;
	  m_world = glm::rotate(m_world,ang,axis);
	  updateworld = true;
	}
	mpos0_l = {texpos.x, texpos.y, 0.f};
	cpos0_l = texpos_to_view(texpos,mpos0_l.z);
      }
    } else { 
      llock = false;
    }
    // updatenone = true; // only if we draw some guiding element
  }

  // double click
  if (hover && IsBindEvent(BIND_NAV_RESET,false)){
    resetView();
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

void View::createTex(float atex){
  const int nsample = 4;

  // Create the multisampling FBO and buffers
  glGenTextures(1, &(FBOtex));
  glGenRenderbuffers(1, &(FBOdepth));
  glGenFramebuffers(1, &(FBO));

  // texture
  glBindTexture(GL_TEXTURE_2D_MULTISAMPLE, FBOtex);
  glTexImage2DMultisample(GL_TEXTURE_2D_MULTISAMPLE, nsample, GL_RGBA, atex, atex, GL_TRUE);
  glTexParameteri(GL_TEXTURE_2D_MULTISAMPLE, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D_MULTISAMPLE, GL_TEXTURE_MAG_FILTER, GL_NEAREST);  
  glBindTexture(GL_TEXTURE_2D_MULTISAMPLE, 0);

  // render buffer
  glBindRenderbuffer(GL_RENDERBUFFER, FBOdepth);
  glRenderbufferStorageMultisample(GL_RENDERBUFFER, nsample, GL_DEPTH_COMPONENT, atex, atex);
  glBindRenderbuffer(GL_RENDERBUFFER, 0);

  // frame buffer
  glBindFramebuffer(GL_FRAMEBUFFER, FBO);
  glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D_MULTISAMPLE, FBOtex, 0); 
  glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, FBOdepth);

  if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
    exit(EXIT_FAILURE);
  glBindFramebuffer(GL_FRAMEBUFFER, 0);

  // non-multisampled FBO and buffers
  glGenTextures(1, &(FBOtex0));
  glGenRenderbuffers(1, &(FBOdepth0));
  glGenFramebuffers(1, &(FBO0));

  // texture
  glBindTexture(GL_TEXTURE_2D, FBOtex0);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, atex, atex, 0, GL_RGBA, GL_UNSIGNED_BYTE, 0);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);  
  glBindTexture(GL_TEXTURE_2D, 0);

  // render buffer
  glBindRenderbuffer(GL_RENDERBUFFER, FBOdepth0);
  glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT, atex, atex);
  glBindRenderbuffer(GL_RENDERBUFFER, 0);
  
  // frame buffer
  glBindFramebuffer(GL_FRAMEBUFFER, FBO0);
  glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, FBOtex0, 0); 
  glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, FBOdepth0);

  if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
    exit(EXIT_FAILURE);
  glBindFramebuffer(GL_FRAMEBUFFER, 0);

  FBO_atex = atex;
}

void View::deleteTex(){
  glDeleteTextures(1, &FBOtex);
  glDeleteRenderbuffers(1, &FBOdepth);
  glDeleteFramebuffers(1, &FBO);
  glDeleteTextures(1, &FBOtex0);
  glDeleteRenderbuffers(1, &FBOdepth0);
  glDeleteFramebuffers(1, &FBO0);
}

bool View::updateTexSize(){
  float amax = ((dock && dock->size.x > 0.f && dock->size.y > 0.f)? fmax(dock->size.x,dock->size.y) : 200.f);
  bool redraw = false;

  if (amax >= FBO_atex){
    deleteTex();
    createTex(ceil(1.5f * amax));
    redraw = true;
  }

  if (FBO_a != amax){
    FBO_a = amax;
    redraw = true;
  }
  return redraw;
}

// Align the view with a given scene axis. a,b,c = 1,2,3 and x,y,z =
// -1,-2,-3. Returns true if the view was updated.
bool View::alignViewAxis(int iaxis){
  if (iscene == 0 || c2::isinit < 2) return false;

  glm::vec3 oaxis = {0.f,0.f,1.f};
  glm::vec3 naxis;
  switch(iaxis){
  case 1: // a
    naxis = {c2::avec[0][0],c2::avec[0][1],c2::avec[0][2]};
    naxis = glm::normalize(naxis);
    break;
  case 2: // b
    naxis = {c2::avec[1][0],c2::avec[1][1],c2::avec[1][2]};
    naxis = glm::normalize(naxis);
    break;
  case 3: // c
    naxis = {c2::avec[2][0],c2::avec[2][1],c2::avec[2][2]};
    naxis = glm::normalize(naxis);
    break;
  case -1: // x
    naxis = {1.0f,0.f,0.f}; break;
  case -2: // y
    naxis = {0.f,1.0f,0.f}; break;
  case -3: // z
    naxis = {0.f,0.f,1.0f}; break;
  default:
    return false;
  }
  glm::vec3 raxis = glm::cross(oaxis,naxis);
  float angle = std::asin(glm::length(raxis));
  resetView();
  if (angle < 1e-10f){
    m_world = glm::mat4(1.0f);
  } else {
    raxis = glm::normalize(raxis);
    m_world = glm::rotate(-angle,raxis);
  }
  updateWorld();
  updateView();
  if (isortho)
    updateProjection();
  return true;
}

void View::resetView(){
  v_front  = {0.f,0.f,-1.f};
  v_up     = {0.f,1.f,0.f};
  v_pos[0] = v_pos[1] = 0.f;
  float scaledsrad = c2::scenerad * std::sqrt(std::max(ncell[0],std::max(ncell[1],ncell[2])));
  v_pos[2] = iscene > 0? resetd * scaledsrad / (tan(0.5f*glm::radians(zfov))):10.f;
  m_world = glm::mat4(1.0f);
  llock = false;
  rlock = false;
}

void View::updateProjection(){
  if (isortho){
    float hw2 = tan(0.5f*glm::radians(zfov)) * v_pos[2];
    m_projection = glm::ortho(-hw2,hw2,-hw2,hw2,znear,1000.f);
  } else {
    m_projection = glm::infinitePerspective(glm::radians(zfov),1.0f,znear);
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

glm::vec3 View::cam_world_coords(){
  glm::vec4 pos4 = glm::inverse(m_world) * glm::vec4(v_pos,1.0f);
  return glm::vec3(pos4.x/pos4.w,pos4.y/pos4.w,pos4.z/pos4.w);
}

glm::vec3 View::cam_view_coords(){
  return v_pos;
}

void View::pos_to_ntexpos(glm::vec2 &pos){
  float x = vrect.Max.x - vrect.Min.x;
  float y = vrect.Max.y - vrect.Min.y;
  float xratio = 2.f * x / fmax(x,y);
  float yratio = 2.f * y / fmax(x,y);
  
  pos.x = ((pos.x - vrect.Min.x) / x - 0.5f) * xratio;
  pos.y = (0.5f - (pos.y - vrect.Min.y) / y) * yratio;
}

void View::ntexpos_to_pos(glm::vec2 &pos){
  float x = vrect.Max.x - vrect.Min.x;
  float y = vrect.Max.y - vrect.Min.y;
  float xratio1 = 0.5f * fmax(x,y) / x;
  float yratio1 = 0.5f * fmax(x,y) / y;
  
  pos.x = vrect.Min.x + x * (0.5f + xratio1 * pos.x);
  pos.y = vrect.Min.y + y * (0.5f - yratio1 * pos.y);
}

void View::pos_to_texpos(glm::vec2 &pos){
  pos_to_ntexpos(pos);
  ntexpos_to_texpos(pos);
}

void View::texpos_to_pos(glm::vec2 &pos){
  texpos_to_ntexpos(pos);
  ntexpos_to_pos(pos);
}

void View::texpos_to_ntexpos(glm::vec2 &pos){
  pos = (pos / FBO_a) * 2.f - 1.f;
}

void View::ntexpos_to_texpos(glm::vec2 &pos){
  pos = (0.5f * pos + 0.5f) * FBO_a;
}

glm::vec2 View::world_to_texpos(glm::vec3 pos){
  const glm::vec4 viewport = {0.f,0.f,FBO_a,FBO_a};
  glm::vec3 pos3 = project(pos,m_view * m_world,m_projection,viewport);
  return glm::vec2(pos3);
}

// dist=0, znear; dist<0, scene origin plane; dist>0, distance from camera
glm::vec3 View::texpos_to_world(glm::vec2 pos, float dist/*=-1.f*/){
  const glm::vec4 viewport = {0.f,0.f,FBO_a,FBO_a};
  glm::vec3 wpos = {};
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
    glm::vec3 origpos = texpos_to_world(pos,-1.f);
    glm::vec3 nearpos = texpos_to_world(pos,0.f);
    glm::vec3 dir = glm::normalize(origpos - nearpos);
    wpos = nearpos + fmax(dist - znear,0.f) * dir;
  }
  return wpos;
}

glm::vec2 View::world_to_ntexpos(glm::vec3 pos){
  glm::vec2 pos2 = world_to_texpos(pos);
  texpos_to_ntexpos(pos2);
  return pos2;
}

// dist=0, znear; dist<0, scene origin plane; dist>0, distance from camera
glm::vec3 View::ntexpos_to_world(glm::vec2 pos, float dist/*=-1.f*/){
  ntexpos_to_texpos(pos);
  return texpos_to_world(pos,dist);
}

glm::vec2 View::view_to_texpos(glm::vec3 pos, float *depth){
  const glm::vec4 viewport = {0.f,0.f,FBO_a,FBO_a};
  glm::vec3 pos3 = project(pos,m_view,m_projection,viewport);
  *depth = pos3.z;
  return glm::vec2(pos3);
}

glm::vec3 View::texpos_to_view(glm::vec2 pos, float depth){
  const glm::vec4 viewport = {0.f,0.f,FBO_a,FBO_a};
  glm::vec3 wpos = {pos.x,pos.y,depth};
  wpos = unProject(wpos,m_view,m_projection,viewport);
  return wpos;
}

float View::texpos_viewdepth(glm::vec2 texpos){
    float depth;
    glBindFramebuffer(GL_FRAMEBUFFER, FBO0);
    glReadPixels(texpos.x,texpos.y, 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &depth);
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    return depth;
}

void View::drawSphere(glm::vec3 r0, float rad, glm::vec4 rgb, int res, bool blend){
  if (blend){
    glEnable(GL_BLEND);
    glEnable(GL_CULL_FACE);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glDepthMask(0);
  }
  glBindVertexArray(sphVAO[res]);

  glm::mat4 m_model = glm::mat4(1.0f);
  m_model = glm::translate(m_model,r0);
  m_model = glm::scale(m_model,glm::vec3(rad,rad,rad));
  glm::mat3 m_normrot = glm::transpose(glm::inverse(glm::mat3(m_view) * glm::mat3(m_world) * glm::mat3(m_model)));

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

void View::drawCylinder(glm::vec3 r1, glm::vec3 r2, float rad, glm::vec4 rgb, int res, bool blend){
  if (blend){
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glDepthMask(0);
  }
  glBindVertexArray(cylVAO[res]);
  glm::vec3 xmid = 0.5f * (r1 + r2);
  glm::vec3 xdif = r2 - r1;
  float blen = glm::length(xdif);
  xdif = xdif/blen;
  glm::vec3 up = {0.f,0.f,1.f};
  glm::vec3 crs = glm::cross(up,xdif);
  
  glm::mat4 m_model = glm::mat4(1.0f);
  m_model = glm::translate(m_model,xmid);
  if (glm::length(crs) > 1.e-12)
    m_model = m_model * glm::rotate(acos(dot(xdif,up)),crs);
  m_model = glm::scale(m_model,glm::vec3(rad,rad,blen));
  glm::mat3 m_normrot = glm::transpose(glm::inverse(glm::mat3(m_view) * glm::mat3(m_world) * glm::mat3(m_model)));
  shader->setMat3("normrot",value_ptr(m_normrot));
  shader->setVec4("vColor",value_ptr(rgb));
  shader->setMat4("model",value_ptr(m_model));
  glDrawElements(GL_TRIANGLES, 3*cylnel[icylres], GL_UNSIGNED_INT, 0);
  if (blend){
    glDisable(GL_BLEND);
    glDepthMask(1);
  }
}

void View::drawUnitCell(glm::vec3 &v0, glm::vec3 &vx, glm::vec3 &vy, glm::vec3 &vz, bool colors){
  const float cellthick = 0.1f;
  glm::vec4 ucellrgbx, ucellrgby, ucellrgbz, ucellrgbo;
  if (colors){
    ucellrgbx = {1.f,0.f,0.f,1.f};
    ucellrgby = {0.f,1.f,0.f,1.f};
    ucellrgbz = {0.f,0.f,1.f,1.f};
  } else {
    ucellrgbx = {0.5f,0.5f,0.5f,1.f};
    ucellrgby = {0.5f,0.5f,0.5f,1.f};
    ucellrgbz = {0.5f,0.5f,0.5f,1.f};
  }
  ucellrgbo = {1.0f,1.0f,1.0f,1.f};
  
  drawCylinder(v0,v0+vx,cellthick,ucellrgbx,0,false);
  drawCylinder(v0,v0+vy,cellthick,ucellrgby,0,false);
  drawCylinder(v0,v0+vz,cellthick,ucellrgbz,0,false);
  for (int ix=0; ix<ncell[0]; ix++){
    for (int iy=0; iy<ncell[1]; iy++){
      for (int iz=0; iz<ncell[2]; iz++){
	glm::vec3 lvec = (float) ix * vx + (float) iy * vy + (float) iz * vz;

	if (ix == 0 && iy == 0 && iz > 0)
	  drawCylinder(v0+lvec,v0+vz+lvec,cellthick,ucellrgbo,0,false);
	if (ix == 0 && iz == 0 && iy > 0)
	  drawCylinder(v0+lvec,v0+vy+lvec,cellthick,ucellrgbo,0,false);
	if (iy == 0 && iz == 0 && ix > 0)
	  drawCylinder(v0+lvec,v0+vx+lvec,cellthick,ucellrgbo,0,false);
	if (iz == 0) { 
	  drawCylinder(v0+vx+lvec,v0+vx+vy+lvec,cellthick,ucellrgbo,0,false);
	  drawCylinder(v0+vy+lvec,v0+vy+vx+lvec,cellthick,ucellrgbo,0,false);
	}
	if (ix == 0){
	  drawCylinder(v0+vy+lvec,v0+vy+vz+lvec,cellthick,ucellrgbo,0,false);
	  drawCylinder(v0+vz+lvec,v0+vz+vy+lvec,cellthick,ucellrgbo,0,false);
	}
	if (iy == 0){
	  drawCylinder(v0+vx+lvec,v0+vx+vz+lvec,cellthick,ucellrgbo,0,false);
	  drawCylinder(v0+vz+lvec,v0+vz+vx+lvec,cellthick,ucellrgbo,0,false);
	}
	drawCylinder(v0+vx+vy+lvec,v0+vx+vy+vz+lvec,cellthick,ucellrgbo,0,false);
	drawCylinder(v0+vx+vz+lvec,v0+vx+vy+vz+lvec,cellthick,ucellrgbo,0,false);
	drawCylinder(v0+vy+vz+lvec,v0+vx+vy+vz+lvec,cellthick,ucellrgbo,0,false);
      }
    }
  }
}

View *CreateView(char *title, int iscene/*=0*/){
  View *aview = new View;

  // save the scene id and set the pointers to that scene
  aview->iscene = iscene;
  aview->title = title;

  // create the texture
  aview->createTex(FBO_tex_a);

  // scene pointers
  if (iscene > 0)
    c2::set_scene_pointers(iscene);

  // set default settings
  aview->SetDefaults();

  // initialize the camera vectors
  aview->resetView();

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

void ForceUpdateAllViews(){
  for (auto iv : viewlist){
    iv->updatescene = true;
  }
}

void SetDefaultAllViews(){
  for (auto iv : viewlist){
    iv->iswire = view_wireframe;
    iv->isortho = view_orthogonal;
    iv->zfov = view_fov;
    iv->resetd = view_resetdistance; 
    for (int i=0; i<4; i++)
      iv->bgrgb[i] = view_bgrgb[i];
    iv->show_atoms = view_show_atoms;
    iv->scale_atoms = view_scale_atoms;
    iv->isphres = view_isphres;
    iv->show_bonds = view_show_bonds;
    iv->scale_bonds = view_scale_bonds;
    iv->icylres = view_icylres;

    iv->updateProjection();
    iv->updatescene = true;
  }
}

