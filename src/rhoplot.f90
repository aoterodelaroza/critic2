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
!

!> Basic plotting capabilities: contour diagrams, 1d, 2d and 3d representations.
module rhoplot
  implicit none

  private

  public :: rhoplot_point
  public :: rhoplot_line
  public :: rhoplot_cube
  public :: rhoplot_plane
  private :: contour
  private :: relief
  private :: colormap
  private :: hallarpuntos
  private :: ordenarpuntos
  private :: linea
  public :: rhoplot_grdvec
  private :: plotvec
  private :: wrtpath
  private :: autochk
  private :: write_fichlabel
  private :: write_fichgnu
  
  interface
     module subroutine rhoplot_point(line)
       character*(*), intent(in) :: line
     end subroutine rhoplot_point
     module subroutine rhoplot_line(line)
       character*(*), intent(in) :: line
     end subroutine rhoplot_line
     module subroutine rhoplot_cube(line)
       character*(*), intent(in) :: line
     end subroutine rhoplot_cube
     module subroutine rhoplot_plane(line)
       character*(*), intent(in) :: line
     end subroutine rhoplot_plane
     module subroutine contour(ff,r0,r1,r2,nx,ny,niso,ziso,rootname,dognu,dolabels)
       integer, intent(in) :: nx, ny
       real*8, intent(in) :: ff(nx,ny)
       real*8, intent(in) :: r0(3), r1(3), r2(3)
       integer, intent(in) :: niso
       real*8, intent(in) :: ziso(niso)
       character*(*), intent(in) :: rootname
       logical, intent(in) :: dognu, dolabels
     end subroutine contour
     module subroutine relief(rootname,outfile,zmin,zmax)
       real*8, intent(in) :: zmin, zmax
       character*(*), intent(in) :: rootname, outfile
     end subroutine relief
     module subroutine colormap(rootname,outfile,cmopt)
       character*(*), intent(in) :: rootname, outfile
       integer, intent(in) :: cmopt
     end subroutine colormap
     module subroutine hallarpuntos(ff,zc,x,y,nx,ny)
       integer, intent(in) :: nx, ny
       real*8, intent(in) :: ff(nx,ny)
       real*8, intent(in) :: zc
       real*8, intent(in) :: x(:), y(:)
     end subroutine hallarpuntos
     module subroutine ordenarpuntos (luw,calpha,ziso)
       integer, intent(in) :: luw
       real*8, intent(in) :: calpha, ziso
     end subroutine ordenarpuntos
     module subroutine linea (x,y,n,luw,ziso)
       integer, intent(in) :: n
       real*8, dimension(n), intent(in) :: x, y
       integer, intent(in) :: luw
       real*8, intent(in) :: ziso
     end subroutine linea
     module subroutine rhoplot_grdvec()
     end subroutine rhoplot_grdvec
     module subroutine plotvec (r0,r1,r2,autocheck,udat)
       integer, intent(in) :: udat
       logical, intent(in) :: autocheck
       real*8, dimension(3), intent(in) :: r0, r1, r2
     end subroutine plotvec
     module subroutine wrtpath(xpath,nptf,udat,rp0,r01,r02,cosalfa,sinalfa)
       use types, only: gpathp
       type(gpathp) :: xpath(*)
       integer, intent(in) :: nptf
       integer, intent(in) :: udat
       real*8, intent(in) :: rp0(3), r01, r02, cosalfa, sinalfa
     end subroutine wrtpath
     module subroutine autochk(rp0)
       real*8 :: rp0(3)
     end subroutine autochk
     module subroutine write_fichlabel(rootname)
       character*(*), intent(in) :: rootname
     end subroutine write_fichlabel
     module subroutine write_fichgnu(rootname,dolabels,docontour,dograds)
       character*(*), intent(in) :: rootname
       logical, intent(in) :: dolabels, docontour, dograds
     end subroutine write_fichgnu
  end interface

end module rhoplot

