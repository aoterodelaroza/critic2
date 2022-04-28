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

#include "scene.h"
#include "shader.h"
#include "critic2.h"

#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/rotate_vector.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/glm.hpp>

const float ofov = 45.0f; // "field of view" for orthogonal projection

Scene::Scene(int isc, float atex_){
  loadShader(false);
  iscene = isc;
  grabFromC2();
  setDefaults();
  resetView();
  updateAllMatrix();  
  setTextureSize(atex_);
}

void Scene::loadShader(bool reload){
  std::string vpath = std::string(c2::c2home) + std::string("/shaders/phong.vs");
  std::string fpath = std::string(c2::c2home) + std::string("/shaders/phong.fs");
  if (reload)
    shphong->loadShader(vpath.c_str(),fpath.c_str(),false);
  else
    shphong = new Shader(vpath.c_str(),fpath.c_str());

  vpath = std::string(c2::c2home) + std::string("/shaders/text.vs");
  fpath = std::string(c2::c2home) + std::string("/shaders/text.fs");
  if (reload)
    shtext->loadShader(vpath.c_str(),fpath.c_str(),false);
  else
    shtext = new Shader(vpath.c_str(),fpath.c_str());
  shcurrent = nullptr;
  if (reload)
    updateAllShaderUniforms();
}

void Scene::usephong(){
  if (shphong && shcurrent != shphong){
    shphong->use();
    shcurrent = shphong;
  }
}

void Scene::usetext(){
  if (shtext && shcurrent != shtext){
    shtext->use();
    shcurrent = shtext;
  }
}

void Scene::grabFromC2(){
  c2::set_scene_pointers(iscene);
  ismolecule = c2::ismolecule;
  ismol3d = c2::ismol3d;
  scenerad = c2::scenerad;
  for (int i=0; i<=2; i++)
    for (int j=0; j<=2; j++)
      avec[i][j] = c2::avec[i][j];
}

void Scene::setDefaults(){
  isucell = (iscene > 0 && !ismolecule);
  ismolcell = (iscene > 0 && ismolecule);
  isborder = (iscene > 0 && !ismolecule && !ismol3d);
  ismotif = (iscene > 0 && !ismolecule);
  ncell[0] = ncell[1] = ncell[2] = 1;

  resetd = view_resetdistance; 
  zfov = view_fov;
  show_atoms = view_show_atoms;
  scale_atoms = view_scale_atoms;
  isphres = view_isphres;
  show_bonds = view_show_bonds;
  scale_bonds = view_scale_bonds;
  show_labels = view_show_labels;
  format_labels = view_format_labels;
  lat_labels = view_lat_labels;
  scale_labels = view_scale_labels;
  icylres = view_icylres;
  iswire = view_wireframe;
  isortho = view_orthogonal;

  setLightpos(view_lightpos);
  setLightcolor(view_lightcolor);
  setAmbient(view_ambient);
  setDiffuse(view_diffuse);
  setSpecular(view_specular);
  setShininess(view_shininess);
  setTextColor(glm::vec3(view_rgb_labels[0],view_rgb_labels[1],view_rgb_labels[2]));
}

void Scene::setLightpos(glm::vec3 lightpos_){
  lightpos = lightpos_;
  usephong();
  shphong->setVec3("lightPos",value_ptr(lightpos));
}
void Scene::setLightcolor(glm::vec3 lightcolor_){
  lightcolor = lightcolor_;
  usephong();
  shphong->setVec3("lightColor",value_ptr(lightcolor));
}
void Scene::setAmbient(float ambient_){
  ambient = ambient_;
  usephong();
  shphong->setFloat("ambient",ambient);
}
void Scene::setDiffuse(float diffuse_){
  diffuse = diffuse_;
  usephong();
  shphong->setFloat("diffuse",diffuse);
}
void Scene::setSpecular(float specular_){
  specular = specular_;
  usephong();
  shphong->setFloat("specular",specular);
}
void Scene::setShininess(int shininess_){
  shininess = shininess_;
  usephong();
  shphong->setInt("shininess",shininess);
}
void Scene::setTextColor(glm::vec3 textcolor_){
  textcolor = textcolor_;
  usetext();
  shtext->setVec3("textColor",value_ptr(textcolor));
}
void Scene::setTextureSize(float atex_){
  atex = atex_;
  usetext();
  glm::mat4 proj = glm::ortho(0.0f, atex, 0.0f, atex);
  shtext->setMat4("projection",value_ptr(proj));
}

void Scene::resetView(){
  v_front  = {0.f,0.f,-1.f};
  v_up     = {0.f,1.f,0.f};
  v_pos[0] = v_pos[1] = 0.f;
  float scaledsrad = scenerad * std::sqrt(std::max(ncell[0],std::max(ncell[1],ncell[2])));
  if (isortho)
    v_pos[2] = iscene > 0? resetd * scaledsrad / (tan(0.5f*glm::radians(ofov))):10.f;
  else
    v_pos[2] = iscene > 0? resetd * scaledsrad / (tan(0.5f*glm::radians(zfov))):10.f;
  m_world = glm::mat4(1.0f);
}

void Scene::updateAllShaderUniforms(){
  usephong();
  updateAllMatrix();
  shphong->setVec3("lightPos",value_ptr(lightpos));
  shphong->setVec3("lightColor",value_ptr(lightcolor));
  shphong->setFloat("ambient",ambient);
  shphong->setFloat("diffuse",diffuse);
  shphong->setFloat("specular",specular);
  shphong->setInt("shininess",shininess);

  usetext();
  shtext->setVec3("textColor",value_ptr(textcolor));
  glm::mat4 proj = glm::ortho(0.0f, atex, 0.0f, atex);
  shtext->setMat4("projection",value_ptr(proj));
}

void Scene::updateAllMatrix(){
  usephong();
  updateProjection();
  updateView();
  updateWorld();
}

void Scene::updateProjection(){
  if (isortho){
    float hw2 = tan(0.5f*glm::radians(ofov)) * v_pos[2];
    m_projection = glm::ortho(-hw2,hw2,-hw2,hw2,znear,1000.f);
  } else {
    m_projection = glm::infinitePerspective(glm::radians(zfov),1.0f,znear);
  }
  shphong->setMat4("projection",value_ptr(m_projection));
  updatescene = true;
}

void Scene::updateView(){
  m_view = lookAt(v_pos,v_pos+v_front,v_up);
  shphong->setMat4("view",value_ptr(m_view));
  updatescene = true;
}

void Scene::updateWorld(){
  shphong->setMat4("world",value_ptr(m_world));
  updatescene = true;
}

// Align the view with a given scene axis. a,b,c = 1,2,3 and x,y,z =
// -1,-2,-3. Returns true if the view was updated.
bool Scene::alignViewAxis(int iaxis){
  glm::vec3 oaxis = {0.f,0.f,1.f};
  glm::vec3 naxis;
  switch(iaxis){
  case 1: // a
    naxis = {avec[0][0],avec[0][1],avec[0][2]};
    naxis = glm::normalize(naxis);
    break;
  case 2: // b
    naxis = {avec[1][0],avec[1][1],avec[1][2]};
    naxis = glm::normalize(naxis);
    break;
  case 3: // c
    naxis = {avec[2][0],avec[2][1],avec[2][2]};
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
  updateAllMatrix();
  return true;
}

glm::vec3 Scene::cam_world_coords(){
  glm::vec4 pos4 = glm::inverse(m_world) * glm::vec4(v_pos,1.0f);
  return glm::vec3(pos4.x/pos4.w,pos4.y/pos4.w,pos4.z/pos4.w);
}

glm::vec3 Scene::cam_view_coords(){
  return v_pos;
}

