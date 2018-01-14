! Copyright (c) 2015 Alberto Otero de la Roza <aoterodelaroza@gmail.com>,
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

! Output in several common graphics formats
module graphics
  implicit none
  
  private

  type grhandle
     integer :: lu = 0 !< Logical unit for the graphics file
     integer :: lumtl = 0 !< Logical unit for the matierlas file (obj format)
     integer :: fmt = 0 !< File format
     character*(255) :: file = "" !< File name
     integer :: nball = 0 !< Number of balls
     integer :: nstick = 0 !< Number of sticks
     integer :: nsurf = 0 !< Number of surfaces
     integer :: nface = 0 !< Number of faces
     integer :: nmtl = 0 !< Number of materials
     integer :: nv = 0 !< Number of vertices
     integer :: nf = 0 !< Number of faces
     integer, allocatable :: mtlrgb(:,:) !< Material definitions
   contains
     procedure :: open => graphics_open
     procedure :: close => graphics_close
     procedure :: ball => graphics_ball
     procedure :: polygon => graphics_polygon
     procedure :: stick => graphics_stick
     procedure :: surf => graphics_surf
  end type grhandle
  public :: grhandle

  private :: graphics_init
  private :: graphics_open
  private :: graphics_close
  private :: graphics_ball
  private :: graphics_polygon
  private :: graphics_stick
  private :: graphics_surf
  private :: obj_open
  private :: obj_close
  private :: obj_ball
  private :: obj_polygon
  private :: obj_stick
  private :: obj_surf
  private :: register_texture
  private :: ply_open
  private :: ply_close
  private :: ply_ball
  private :: ply_polygon
  private :: ply_stick
  private :: ply_surf
  private :: off_open
  private :: off_close
  private :: off_ball
  private :: off_polygon
  private :: off_stick
  private :: off_surf
  
  interface
     module subroutine graphics_init()
     end subroutine graphics_init
     module subroutine graphics_open(g,fmt,file)
       class(grhandle), intent(inout) :: g
       character*3, intent(in) :: fmt
       character*(*), intent(in) :: file
     end subroutine graphics_open
     module subroutine graphics_close(g)
       class(grhandle), intent(inout) :: g
     end subroutine graphics_close
     module subroutine graphics_ball(g,x,rgb,r)
       class(grhandle), intent(inout) :: g
       real*8, intent(in) :: x(3)
       integer, intent(in) :: rgb(3)
       real*8, intent(in) :: r
     end subroutine graphics_ball
     module subroutine graphics_polygon(g,x,rgb)
       class(grhandle), intent(inout) :: g
       real*8, intent(in) :: x(:,:)
       integer, intent(in) :: rgb(3)
     end subroutine graphics_polygon
     module subroutine graphics_stick(g,x1,x2,rgb,r)
       class(grhandle), intent(inout) :: g
       real*8, intent(in) :: x1(3), x2(3)
       integer, intent(in) :: rgb(3)
       real*8, intent(in) :: r
     end subroutine graphics_stick
     module subroutine graphics_surf(g,srf,fsurf)
       use surface, only: minisurf
       class(grhandle), intent(inout) :: g
       type(minisurf), intent(in) :: srf
       real*8, intent(in), optional :: fsurf(:)
     end subroutine graphics_surf
     module subroutine obj_open(g)
       type(grhandle), intent(inout) :: g
     end subroutine obj_open
     module subroutine obj_close(g)
       type(grhandle), intent(inout) :: g
     end subroutine obj_close
     module subroutine obj_ball(g,x,rgb,r)
       type(grhandle), intent(inout) :: g
       real*8, intent(in) :: x(3)
       integer, intent(in) :: rgb(3)
       real*8, intent(in) :: r
     end subroutine obj_ball
     module subroutine obj_polygon(g,x,rgb)
       type(grhandle), intent(inout) :: g
       real*8, intent(in) :: x(:,:)
       integer, intent(in) :: rgb(3)
     end subroutine obj_polygon
     module subroutine obj_stick(g,x1,x2,rgb,r)
       type(grhandle), intent(inout) :: g
       real*8, intent(in) :: x1(3), x2(3)
       integer, intent(in) :: rgb(3)
       real*8, intent(in) :: r
     end subroutine obj_stick
     module subroutine obj_surf(g,srf,fsurf)
       use surface, only: minisurf
       type(grhandle), intent(inout) :: g
       type(minisurf), intent(in) :: srf
       real*8, intent(in), optional :: fsurf(:)
     end subroutine obj_surf
     module function register_texture(g,rgb) result(imtl)
       type(grhandle), intent(inout) :: g
       integer, intent(in) :: rgb(3)
       integer :: imtl
     end function register_texture
     module subroutine ply_open(g)
       type(grhandle), intent(inout) :: g
     end subroutine ply_open
     module subroutine ply_close(g)
       type(grhandle), intent(inout) :: g
     end subroutine ply_close
     module subroutine ply_ball(g,x,rgb,r)
       type(grhandle), intent(inout) :: g
       real*8, intent(in) :: x(3)
       integer, intent(in) :: rgb(3)
       real*8, intent(in) :: r
     end subroutine ply_ball
     module subroutine ply_polygon(g,x,rgb)
       type(grhandle), intent(inout) :: g
       real*8, intent(in) :: x(:,:)
       integer, intent(in) :: rgb(3)
     end subroutine ply_polygon
     module subroutine ply_stick(g,x1,x2,rgb,r)
       type(grhandle), intent(inout) :: g
       real*8, intent(in) :: x1(3), x2(3)
       integer, intent(in) :: rgb(3)
       real*8, intent(in) :: r
     end subroutine ply_stick
     module subroutine ply_surf(g,srf,fsurf)
       use surface, only: minisurf
       type(grhandle), intent(inout) :: g
       type(minisurf), intent(in) :: srf
       real*8, intent(in), optional :: fsurf(:)
     end subroutine ply_surf
     module subroutine off_open(g)
       type(grhandle), intent(inout) :: g
     end subroutine off_open
     module subroutine off_close(g)
       type(grhandle), intent(inout) :: g
     end subroutine off_close
     module subroutine off_ball(g,x,rgb,r)
       type(grhandle), intent(inout) :: g
       real*8, intent(in) :: x(3)
       integer, intent(in) :: rgb(3)
       real*8, intent(in) :: r
     end subroutine off_ball
     module subroutine off_polygon(g,x,rgb)
       type(grhandle), intent(inout) :: g
       real*8, intent(in) :: x(:,:)
       integer, intent(in) :: rgb(3)
     end subroutine off_polygon
     module subroutine off_stick(g,x1,x2,rgb,r)
       type(grhandle), intent(inout) :: g
       real*8, intent(in) :: x1(3), x2(3)
       integer, intent(in) :: rgb(3)
       real*8, intent(in) :: r
     end subroutine off_stick
     module subroutine off_surf(g,srf,fsurf)
       use surface, only: minisurf
       type(grhandle), intent(inout) :: g
       type(minisurf), intent(in) :: srf
       real*8, intent(in), optional :: fsurf(:)
     end subroutine off_surf
  end interface

end module graphics
