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

!> Tools for automatic or manual critical point determination.
module autocp
  implicit none

  private 

  public :: autocritic
  public :: cpreport
  private :: critshell
  private :: writechk
  private :: readchk
  private :: atomic_connect_report
  private :: cp_short_report
  private :: barycentric
  private :: barycentric_divide
  private :: seed_from_simplex
  private :: cp_long_report
  private :: cp_vlong_report
  private :: graph_short_report
  private :: makegraph
  private :: scale_ws
  
  interface
     module subroutine autocritic(line)
       character*(*), intent(in) :: line
     end subroutine autocritic
     module subroutine cpreport(line)
       character*(*), intent(in) :: line
     end subroutine cpreport
     module subroutine critshell(shmax)
       integer, intent(in) :: shmax
     end subroutine critshell
     module subroutine writechk()
     end subroutine writechk
     module subroutine readchk()
     end subroutine readchk
     module subroutine atomic_connect_report()
     end subroutine atomic_connect_report
     module subroutine cp_short_report()
     end subroutine cp_short_report
     module subroutine barycentric(iniv,depthmax,nn,xseed)
       real*8, dimension(4,3), intent(in) :: iniv
       integer, intent(in) :: depthmax
       integer, intent(inout) :: nn
       real*8, intent(inout), allocatable :: xseed(:,:)
     end subroutine barycentric
     recursive module subroutine barycentric_divide(n,simp,depth)
       integer, intent(in) :: n
       integer, intent(in) :: simp(:,:)
       integer, intent(in) :: depth
     end subroutine barycentric_divide
     module subroutine seed_from_simplex(simp, dim, nn, xseed)
       real*8, intent(in) :: simp(4,3)
       integer, intent(in) :: dim
       integer, intent(inout) :: nn
       real*8, intent(inout), allocatable :: xseed(:,:)
     end subroutine seed_from_simplex
     module subroutine cp_long_report()
     end subroutine cp_long_report
     module subroutine cp_vlong_report()
     end subroutine cp_vlong_report
     module subroutine graph_short_report()
     end subroutine graph_short_report
     module subroutine makegraph()
     end subroutine makegraph
     module subroutine scale_ws(rad,wso,ntetrag,tetrag)
       real*8, intent(in) :: rad
       real*8, intent(in) :: wso(3)
       integer, intent(inout) :: ntetrag
       real*8, intent(inout) :: tetrag(:,:,:)
     end subroutine scale_ws
  end interface

end module autocp

