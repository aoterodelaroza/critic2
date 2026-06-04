! Copyright (c) 2007-2022 Alberto Otero de la Roza <aoterodelaroza@gmail.com>,
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

! Molecular point-group determination. The point_group type holds the
! symmetry operations and symbol of a finite molecule; calc_point_group
! computes them from a list of atomic positions and atomic numbers. This
! module depends only on low-level modules, so it can be used both by the
! crystal class and by the molecular fragment class.
module molsymmod
  use types, only: molsymop
  implicit none

  private
  public :: point_group
  public :: calc_point_group

  !> Molecular point group: symmetry operations and symbol.
  type point_group
     logical :: avail = .false. ! the point group has been calculated
     logical :: isatom ! whether the system is a single atom
     logical :: islinear ! whether the system is a linear molecule
     logical :: isplanar ! whether the system is a planar molecule
     real*8 :: xcm(3) ! center of mass
     integer :: nop ! number of symmetry operations
     type(molsymop), allocatable :: op(:) ! symmetry operations
     character(len=:), allocatable :: symbol ! symbol for the point group
   contains
     procedure :: clear => point_group_clear
     procedure :: report => point_group_report
     procedure :: init_as_c1 => point_group_init_as_c1
  end type point_group

  interface
     module subroutine calc_point_group(nat,xin,z,pg,errmsg)
       integer, intent(in) :: nat
       real*8, intent(in) :: xin(3,nat)
       integer, intent(in) :: z(nat)
       type(point_group), intent(inout) :: pg
       character(len=:), allocatable, intent(out) :: errmsg
     end subroutine calc_point_group
     module subroutine point_group_clear(p)
       class(point_group), intent(inout) :: p
     end subroutine point_group_clear
     module subroutine point_group_report(p)
       class(point_group), intent(inout) :: p
     end subroutine point_group_report
     module subroutine point_group_init_as_c1(p)
       class(point_group), intent(inout) :: p
     end subroutine point_group_init_as_c1
  end interface

end module molsymmod
