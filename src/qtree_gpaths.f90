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

!> Qtree, gradient path tracing.
module qtree_gpaths
  implicit none

  private
  public :: gradient_qtree
  public :: gradient_color
  public :: gradient_full

  interface
     module subroutine gradient_qtree(xp,base_t,res,iidx,cpout,ier,usestack,trm,fgr,lapgr)
       use qtree_basic, only: qtreeidx, qtreei, qtreer
       use types, only: scalar_value
       real*8, intent(inout) :: xp(3)
       integer, intent(inout) :: base_t
       type(scalar_value), intent(inout) :: res
       integer(qtreeidx), intent(in) :: iidx
       integer, intent(out) :: cpout
       integer, intent(out) :: ier
       logical, intent(in) :: usestack
       integer(qtreei), intent(inout) :: trm(:,:)
       real(qtreer), intent(inout) :: fgr(:,:), lapgr(:,:)
     end subroutine gradient_qtree
     module subroutine gradient_color(xp,base_t,rver,res,cpout,ier,trm)
       use qtree_basic, only: qtreei
       use types, only: scalar_value
       real*8, intent(inout) :: xp(3)
       integer, intent(inout) :: base_t
       real*8, intent(inout) :: rver(3)
       type(scalar_value), intent(inout) :: res
       integer, intent(out) :: cpout
       integer, intent(out) :: ier
       integer(qtreei), intent(inout) :: trm(:,:)
     end subroutine gradient_color
     module subroutine gradient_full(xp,base_t,rver,res,cpout,ier)
       use types, only: scalar_value
       real*8, intent(inout) :: xp(3)
       integer, intent(inout) :: base_t
       real*8, intent(inout) :: rver(3)
       type(scalar_value), intent(inout) :: res
       integer, intent(out) :: cpout
       integer, intent(out) :: ier
     end subroutine gradient_full
  end interface

end module qtree_gpaths
