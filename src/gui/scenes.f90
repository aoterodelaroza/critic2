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
module scenes
  use iso_c_binding
  implicit none

  private

  type scene
     logical :: isinit = .false. ! whether the scene has been initialized
     integer :: id ! system ID
     logical :: showcell = .false. ! show the unit cell?
     integer :: ncell(3) ! number of cells show in each direction
     integer :: imotif = 0 ! 0 = none, 1 = border, 2 = onemotif, 3 = molmotif
     real*8 :: scenerad = 1d0 ! scene radius
     real*8 :: scenecenter(3) ! scene center
     ! phong shader settings
     real(c_float) :: lightpos(3) ! light position
     real(c_float) :: lightcolor(3) ! light color
     real(c_float) :: ambient ! ambient light coefficent
     real(c_float) :: diffuse ! diffuse light coefficent
     real(c_float) :: specular ! specular light coefficent
     integer(c_int) :: shininess ! shininess parameter
     ! scene transformation matrices
     real(c_float) :: ortho_fov ! orthographic field of view
     real(c_float) :: persp_fov ! perspective field of view
     real(c_float) :: znear ! position of the near plane
     real(c_float) :: zfar  ! position of the far plane
     real(c_float) :: campos(3) ! position of the camera
     real(c_float) :: world(4,4) ! world transform matrix
     real(c_float) :: view(4,4) ! view transform matrix
     real(c_float) :: projection(4,4) ! projection transform matrix
   contains
     procedure :: init => scene_init
     procedure :: end => scene_end
     procedure :: reset => scene_reset
     procedure :: render => scene_render
     procedure :: update_projection_matrix
     procedure :: update_view_matrix
     procedure :: draw_sphere
     procedure :: draw_cylinder
  end type scene
  public :: scene

  ! module procedure interfaces
  interface
     module subroutine scene_init(s,isys)
       class(scene), intent(inout), target :: s
       integer, intent(in) :: isys
     end subroutine scene_init
     module subroutine scene_end(s)
       class(scene), intent(inout), target :: s
     end subroutine scene_end
     module subroutine scene_reset(s)
       class(scene), intent(inout), target :: s
     end subroutine scene_reset
     module subroutine scene_render(s)
       class(scene), intent(inout), target :: s
     end subroutine scene_render
     module subroutine update_projection_matrix(s)
       class(scene), intent(inout), target :: s
     end subroutine update_projection_matrix
     module subroutine update_view_matrix(s)
       class(scene), intent(inout), target :: s
     end subroutine update_view_matrix
     module subroutine draw_sphere(s,x0,rad,rgba)
       class(scene), intent(inout), target :: s
       real(c_float), intent(in) :: x0(3)
       real(c_float), intent(in) :: rad
       real(c_float), intent(in) :: rgba(4)
     end subroutine draw_sphere
     module subroutine draw_cylinder(s,x1,x2,rad,rgba)
       class(scene), intent(inout), target :: s
       real(c_float), intent(in) :: x1(3)
       real(c_float), intent(in) :: x2(3)
       real(c_float), intent(in) :: rad
       real(c_float), intent(in) :: rgba(4)
     end subroutine draw_cylinder
  end interface

end module scenes

