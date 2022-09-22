! Copyright (c) 2019-2022 Alberto Otero de la Roza <aoterodelaroza@gmail.com>,
! Ángel Martín Pendás <angel@fluor.quimica.uniovi.es> and Víctor Luaña
! <victor@fluor.quimica.uniovi.es>.
!
! critic2 is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or (at
! your option) any later version.
!
! critic2 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

! Scene object and GL rendering utilities
submodule (scenes) proc
  implicit none

  real*8 :: min_scenerad = 10.d0

contains

  !> Initialize a scene object associated with system isys.
  module subroutine scene_init(s,isys)
    use gui_main, only: nsys, sysc, sys_init, sys
    class(scene), intent(inout), target :: s
    integer, intent(in) :: isys

    real*8 :: xmin(3), xmax(3), x(3)
    integer :: i, i1, i2, i3

    if (isys < 1 .or. isys > nsys) return
    if (sysc(isys)%status /= sys_init) return

    ! basic variables
    s%id = isys
    s%showcell = .not.(sys(isys)%c%ismolecule)
    s%ncell = 1
    s%isinit = .true.

    ! default imotif
    if (sys(isys)%c%ismolecule) then
       s%imotif = 0
    elseif (sys(isys)%c%ismol3d) then
       s%imotif = 3
    else
       s%imotif = 1
    end if

    ! calculate the initial scene radius
    xmin = 0d0
    xmax = 0d0
    if (sys(isys)%c%ncel > 0) then
       xmin = sys(isys)%c%atcel(1)%r
       xmax = sys(isys)%c%atcel(1)%r
       do i = 2, sys(isys)%c%ncel
          xmax = max(sys(isys)%c%atcel(i)%r,xmax)
          xmin = min(sys(isys)%c%atcel(i)%r,xmin)
       end do

       if (.not.sys(isys)%c%ismolecule) then
          do i1 = 0, 1
             do i2 = 0, 1
                do i3 = 0, 1
                   x = real((/i1,i2,i3/),8)
                   x = sys(isys)%c%x2c(x)
                   xmax = max(x,xmax)
                   xmin = min(x,xmin)
                end do
             end do
          end do
       end if
    end if
    s%scenerad = max(norm2(xmax-xmin),min_scenerad)

  end subroutine scene_init

  !> Terminate a scene object
  module subroutine scene_end(s)
    class(scene), intent(inout), target :: s

    s%isinit = .false.
    s%id = 0

  end subroutine scene_end

  !> Draw the scene
  module subroutine scene_render(s)
    use interfaces_opengl3
    use gui_main, only: sysc, sys_init
    use shaders, only: shader_phong, useshader
    use shapes, only: testVAO
    class(scene), intent(inout), target :: s

    if (.not.sysc(s%id)%status == sys_init) return

    ! set up the shader and the uniforms
    call useshader(shader_phong)

    ! call glBindVertexArray(testVAO)
    ! call glDrawArrays(GL_TRIANGLES, 0, 3)
    ! call glBindVertexArray(0)

  end subroutine scene_render

end submodule proc

! Scene::Scene(int isc, float atex_){
!   setDefaults();
!   resetView();
!   updateAllMatrix();
!   setTextureSize(atex_);
! }

!   resetd = view_resetdistance;
!   zfov = view_fov;
!   isortho = view_orthogonal;

!   setLightpos(view_lightpos);
!   setLightcolor(view_lightcolor);
!   setAmbient(view_ambient);
!   setDiffuse(view_diffuse);
!   setSpecular(view_specular);
!   setShininess(view_shininess);
!   setTextColor(glm::vec3(view_rgb_labels[0],view_rgb_labels[1],view_rgb_labels[2]));
! }

!   shphong->setVec3("lightPos",value_ptr(lightpos));
!   shphong->setVec3("lightColor",value_ptr(lightcolor));
!   shphong->setFloat("ambient",ambient);
!   shphong->setFloat("diffuse",diffuse);
!   shphong->setFloat("specular",specular);
!   shphong->setInt("shininess",shininess);

! void Scene::setTextureSize(float atex_){
!   atex = atex_;
!   usetext();
!   glm::mat4 proj = glm::ortho(0.0f, atex, 0.0f, atex);
!   shtext->setMat4("projection",value_ptr(proj));
! }

! void Scene::resetView(){
!   v_front  = {0.f,0.f,-1.f};
!   v_up     = {0.f,1.f,0.f};
!   v_pos[0] = v_pos[1] = 0.f;
!   float scaledsrad = scenerad * std::sqrt(std::max(ncell[0],std::max(ncell[1],ncell[2])));
!   if (isortho)
!     v_pos[2] = iscene > 0? resetd * scaledsrad / (tan(0.5f*glm::radians(ofov))):10.f;
!   else
!     v_pos[2] = iscene > 0? resetd * scaledsrad / (tan(0.5f*glm::radians(zfov))):10.f;
!   m_world = glm::mat4(1.0f);
! }

! void Scene::updateAllMatrix(){
!   updateProjection();
!   updateView();
!   updateWorld();
! }

! void Scene::updateProjection(){
!   if (isortho){
!     float hw2 = tan(0.5f*glm::radians(ofov)) * v_pos[2];
!     m_projection = glm::ortho(-hw2,hw2,-hw2,hw2,znear,1000.f);
!   } else {
!     m_projection = glm::infinitePerspective(glm::radians(zfov),1.0f,znear);
!   }
!   shphong->setMat4("projection",value_ptr(m_projection));
!   updatescene = true;
! }

! void Scene::updateView(){
!   m_view = lookAt(v_pos,v_pos+v_front,v_up);
!   shphong->setMat4("view",value_ptr(m_view));
!   updatescene = true;
! }

! void Scene::updateWorld(){
!   shphong->setMat4("world",value_ptr(m_world));
!   updatescene = true;
! }

! // Align the view with a given scene axis. a,b,c = 1,2,3 and x,y,z =
! // -1,-2,-3. Returns true if the view was updated.
! bool Scene::alignViewAxis(int iaxis){
!   glm::vec3 oaxis = {0.f,0.f,1.f};
!   glm::vec3 naxis;
!   switch(iaxis){
!   case 1: // a
!     naxis = {avec[0][0],avec[0][1],avec[0][2]};
!     naxis = glm::normalize(naxis);
!     break;
!   case 2: // b
!     naxis = {avec[1][0],avec[1][1],avec[1][2]};
!     naxis = glm::normalize(naxis);
!     break;
!   case 3: // c
!     naxis = {avec[2][0],avec[2][1],avec[2][2]};
!     naxis = glm::normalize(naxis);
!     break;
!   case -1: // x
!     naxis = {1.0f,0.f,0.f}; break;
!   case -2: // y
!     naxis = {0.f,1.0f,0.f}; break;
!   case -3: // z
!     naxis = {0.f,0.f,1.0f}; break;
!   default:
!     return false;
!   }
!   glm::vec3 raxis = glm::cross(oaxis,naxis);
!   float angle = std::asin(glm::length(raxis));
!   resetView();
!   if (angle < 1e-10f){
!     m_world = glm::mat4(1.0f);
!   } else {
!     raxis = glm::normalize(raxis);
!     m_world = glm::rotate(-angle,raxis);
!   }
!   updateAllMatrix();
!   return true;
! }

! glm::vec3 Scene::cam_world_coords(){
!   glm::vec4 pos4 = glm::inverse(m_world) * glm::vec4(v_pos,1.0f);
!   return glm::vec3(pos4.x/pos4.w,pos4.y/pos4.w,pos4.z/pos4.w);
! }

! glm::vec3 Scene::cam_view_coords(){
!   return v_pos;
! }
