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

!> Qtree, work with tetrahedra.
module qtree_tetrawork
  implicit none

  private
  public :: tetrah_subdivide
  public :: term_rec
  public :: tetrah_paint
  public :: integ_inner_keast
  public :: integ_inner_cubpack
  public :: integ_border_keast
  public :: integ_border_cubpack
  public :: integ_corner
  public :: integ_corner_deferred
  public :: integ_corner_sum
  public :: paint_inside_spheres
  private :: cubpack_f
  
  interface
     module subroutine tetrah_subdivide(base_t,iiv,il,acum_atprop,trm,fgr,lapgr,vgr)
       use qtree_basic, only: qtreei, qtreer
       integer, intent(in) :: base_t
       integer, intent(in) :: iiv(3,4)
       integer, intent(in) :: il
       integer(qtreei), intent(inout) :: trm(:,:)
       real(qtreer), intent(inout) :: fgr(:,:), lapgr(:,:), vgr(:)
       real*8, intent(inout) :: acum_atprop(:,:)
     end subroutine tetrah_subdivide
     module function term_rec(base_t,iver,l,trm,fgr,lapgr)
       use qtree_basic, only: qtreei, qtreer
       integer, intent(in) :: base_t, iver(3), l
       integer(qtreei), intent(inout) :: trm(:,:)
       real(qtreer), intent(inout) :: fgr(:,:), lapgr(:,:)
       integer :: term_rec
     end function term_rec
     module subroutine tetrah_paint(base_t,iv,l,color,trm)
       use qtree_basic, only: qtreei
       integer, intent(in) :: base_t
       integer, intent(in) :: iv(3,4)
       integer, intent(in) :: l
       integer, intent(in) :: color
       integer(qtreei), intent(inout) :: trm(:,:)
     end subroutine tetrah_paint
     module subroutine integ_inner_keast(base_t,iv,l,color,klvl,acum_atprop)
       integer, intent(in) :: base_t
       integer, intent(in) :: iv(3,4)
       integer, intent(in) :: l
       integer, intent(in) :: color
       integer, intent(in) :: klvl
       real*8, intent(inout) :: acum_atprop(:,:)
     end subroutine integ_inner_keast
     module subroutine integ_inner_cubpack(base_t,iv,l,color,acum_atprop)
       integer, intent(in) :: base_t
       integer, intent(in) :: iv(3,4)
       integer, intent(in) :: l
       integer, intent(in) :: color
       real*8, intent(inout) :: acum_atprop(:,:)
     end subroutine integ_inner_cubpack
     module subroutine integ_border_keast(base_t,iv,l,ts,klvl,acum_atprop)
       integer, intent(in) :: base_t
       integer, intent(in) :: iv(3,4)
       integer, intent(in) :: l
       integer, intent(in) :: ts(4)
       integer, intent(in) :: klvl
       real*8, intent(inout) :: acum_atprop(:,:)
     end subroutine integ_border_keast
     module subroutine integ_border_cubpack(base_t,iv,l,ts,acum_atprop)
       integer, intent(in) :: base_t
       integer, intent(in) :: iv(3,4)
       integer, intent(in) :: l
       integer, intent(in) :: ts(4)
       real*8, intent(inout) :: acum_atprop(:,:)
     end subroutine integ_border_cubpack
     module subroutine integ_corner(base_t,iv,l,ts,acum_atprop,fgr,lapgr)
       use qtree_basic, only: qtreer
       integer, intent(in) :: base_t
       integer, intent(in) :: iv(3,4)
       integer, intent(in) :: l
       integer, intent(in) :: ts(4)
       real(qtreer), intent(inout) :: fgr(:,:), lapgr(:,:)
       real*8, intent(inout) :: acum_atprop(:,:)
     end subroutine integ_corner
     module subroutine integ_corner_deferred(base_t,iv,l,ts,vgr)
       use qtree_basic, only: qtreer
       integer, intent(in) :: base_t
       integer, intent(in) :: iv(3,4)
       integer, intent(in) :: l
       integer, intent(in) :: ts(4)
       real(qtreer), intent(inout) :: vgr(:)
     end subroutine integ_corner_deferred
     module subroutine integ_corner_sum(base_t,trm,vgr,acum_atprop)
       use qtree_basic, only: qtreei, qtreer
       integer, intent(in) :: base_t
       integer(qtreei), intent(inout) :: trm(:,:)
       real(qtreer), intent(inout) :: vgr(:)
       real*8, intent(inout) :: acum_atprop(:,:)
     end subroutine integ_corner_sum
     module subroutine paint_inside_spheres(tt,tto,trm)
       use qtree_basic, only: qtreei
       integer, intent(in) :: tt, tto
       integer(qtreei), intent(inout) :: trm(:,:)
     end subroutine paint_inside_spheres
     module function cubpack_f(numfun,x) result(value)
       use precision_model, only: stnd
       integer, intent(in) :: numfun
       real(kind=stnd), dimension(:), intent(in) :: x
       real(kind=stnd), dimension(numfun) :: value
     end function cubpack_f
  end interface

end module qtree_tetrawork
