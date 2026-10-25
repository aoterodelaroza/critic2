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
  public :: rhoplot_cube_values
  public :: rhoplot_cube_gridcheck
  public :: rhoplot_cube_gridvalues
  public :: rhoplot_cube_write
  public :: rhoplot_cube_outform
  public :: rhoplot_plane
  public :: rhoplot_plane_values
  public :: rhoplot_contour_lines
  public :: rhoplot_nti_request
  public :: rhoplot_nti_value
  public :: rhoplot_grdvec

  ! formats of the files written by rhoplot_cube_write
  integer, parameter, public :: cubefmt_cube = 1 ! Gaussian cube
  integer, parameter, public :: cubefmt_bincube = 2 ! binary cube
  integer, parameter, public :: cubefmt_vasp = 3 ! VASP CHGCAR
  integer, parameter, public :: cubefmt_vaspnov = 4 ! VASP ELFCAR (no volume factor)
  integer, parameter, public :: cubefmt_xsf = 5 ! xcrysden xsf

  ! wavefunctions of rhoplot_cube_gridvalues
  integer, parameter, public :: cubepsi_none = 0 ! the grid itself
  integer, parameter, public :: cubepsi_mlwf = 1 ! maximally localized Wannier function
  integer, parameter, public :: cubepsi_wannier = 2 ! Wannier function (no localization)
  integer, parameter, public :: cubepsi_unk = 3 ! periodic part of the Bloch function
  integer, parameter, public :: cubepsi_psink = 4 ! Bloch function

  ! parts of a wavefunction in rhoplot_cube_gridvalues
  integer, parameter, public :: cubepart_abs = 0
  integer, parameter, public :: cubepart_real = 1
  integer, parameter, public :: cubepart_imag = 2

  interface
     module subroutine rhoplot_point(line)
       character*(*), intent(in) :: line
     end subroutine rhoplot_point
     module subroutine rhoplot_line(line)
       character*(*), intent(in) :: line
     end subroutine rhoplot_line
     pure module function rhoplot_nti_request(nti) result(request)
       use types, only: field_evaluation_avail
       integer, intent(in) :: nti
       type(field_evaluation_avail) :: request
     end function rhoplot_nti_request
     pure module function rhoplot_nti_value(res,nti) result(v)
       use types, only: scalar_value
       type(scalar_value), intent(in) :: res
       integer, intent(in) :: nti
       real*8 :: v
     end function rhoplot_nti_value
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
     module subroutine rhoplot_cube_values(s,x0,xd,nn,id,expr,nti,oncell,lf,errmsg)
       use systemmod, only: system
       type(system), intent(inout), target :: s
       real*8, intent(in) :: x0(3), xd(3,3)
       integer, intent(in) :: nn(3), id, nti
       character*(*), intent(in) :: expr
       logical, intent(in) :: oncell
       real*8, allocatable, intent(inout) :: lf(:,:,:)
       character(len=:), allocatable, intent(out) :: errmsg
     end subroutine rhoplot_cube_values
     module function rhoplot_cube_gridcheck(s,id,ipsi,ibnd,supercell) result(errmsg)
       use systemmod, only: system
       type(system), intent(in) :: s
       integer, intent(in) :: id, ipsi, ibnd
       logical, intent(in) :: supercell
       character(len=:), allocatable :: errmsg
     end function rhoplot_cube_gridcheck
     module subroutine rhoplot_cube_gridvalues(s,id,ia,na,ipsi,ibnd,ik,inr,ispin,ipart,f)
       use systemmod, only: system
       type(system), intent(inout) :: s
       integer, intent(in) :: id, ia(3), na(3), ipsi, ibnd, ik, inr(3), ispin, ipart
       real*8, allocatable, intent(inout) :: f(:,:,:)
     end subroutine rhoplot_cube_gridvalues
     module subroutine rhoplot_cube_write(cr,outform,f,file,header,xd,x0,ishift,errmsg)
       use crystalmod, only: crystal
       type(crystal), intent(in) :: cr
       integer, intent(in) :: outform
       real*8, allocatable, intent(in) :: f(:,:,:)
       character*(*), intent(in) :: file
       logical, intent(in) :: header
       real*8, intent(in) :: xd(3,3), x0(3)
       integer, intent(in) :: ishift(3)
       character(len=:), allocatable, intent(out) :: errmsg
     end subroutine rhoplot_cube_write
     module function rhoplot_cube_outform(file) result(outform)
       character*(*), intent(in) :: file
       integer :: outform
     end function rhoplot_cube_outform
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

