! Copyright (c) 2015-2022 Alberto Otero de la Roza <aoterodelaroza@gmail.com>,
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
!

!> Basic plotting capabilities: contour diagrams, 1d, 2d and 3d representations.
module rhoplot
  implicit none

  private

  public :: rhoplot_point
  public :: rhoplot_line
  public :: rhoplot_line_values
  public :: rhoplot_cube
  public :: rhoplot_plane
  public :: rhoplot_plane_values
  public :: rhoplot_contour_lines
  public :: rhoplot_grdvec

  interface
     module subroutine rhoplot_point(line)
       character*(*), intent(in) :: line
     end subroutine rhoplot_point
     module subroutine rhoplot_line(line)
       character*(*), intent(in) :: line
     end subroutine rhoplot_line
     module subroutine rhoplot_line_values(s,x0,x1,np,id,expr,nti,f,g,errmsg)
       use systemmod, only: system
       type(system), intent(inout) :: s
       real*8, intent(in) :: x0(3), x1(3)
       integer, intent(in) :: np, id, nti
       character*(*), intent(in) :: expr
       real*8, intent(out) :: f(np), g(np)
       character(len=:), allocatable, intent(out) :: errmsg
     end subroutine rhoplot_line_values
     module subroutine rhoplot_cube(line)
       character*(*), intent(in) :: line
     end subroutine rhoplot_cube
     module subroutine rhoplot_plane(line)
       character*(*), intent(in) :: line
     end subroutine rhoplot_plane
     module subroutine rhoplot_plane_values(s,x0,uu,vv,nx,ny,id,expr,nti,ff,errmsg,progress)
       use systemmod, only: system
       type(system), intent(inout) :: s
       real*8, intent(in) :: x0(3), uu(3), vv(3)
       integer, intent(in) :: nx, ny, id, nti
       character*(*), intent(in) :: expr
       real*8, intent(out) :: ff(nx,ny)
       character(len=:), allocatable, intent(out) :: errmsg
       logical, intent(in), optional :: progress
     end subroutine rhoplot_plane_values
     module subroutine rhoplot_contour_lines(ff,nx,ny,r01,r02,cosalfa,niso,ziso,nline,lptr,xy,liso)
       integer, intent(in) :: nx, ny
       real*8, intent(in) :: ff(nx,ny)
       real*8, intent(in) :: r01, r02, cosalfa
       integer, intent(in) :: niso
       real*8, intent(in) :: ziso(niso)
       integer, intent(out) :: nline
       integer, allocatable, intent(inout) :: lptr(:)
       real*8, allocatable, intent(inout) :: xy(:,:)
       integer, allocatable, intent(inout) :: liso(:)
     end subroutine rhoplot_contour_lines
     module subroutine rhoplot_grdvec()
     end subroutine rhoplot_grdvec
  end interface

end module rhoplot

