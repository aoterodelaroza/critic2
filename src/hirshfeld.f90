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

! Hirshfeld integration
module hirshfeld
  implicit none

  private

  public :: hirsh_grid
  public :: hirsh_nogrid
  public :: hirsh_weights
  public :: hirsh_cutoffs
  public :: hirsh_rho
  public :: hirsh_option
  public :: voronoi_grid

  ! Hirshfeld-I defaults: convergence threshold on max|dN| and maximum iterations
  real*8, parameter, public :: hirsh_tol_def = 1d-7
  integer, parameter, public :: hirsh_maxit_def = 200

  interface
     module subroutine hirsh_grid(s,bas)
       use systemmod, only: system
       use types, only: basindat
       type(system), intent(inout) :: s
       type(basindat), intent(inout) :: bas
     end subroutine hirsh_grid
     module subroutine hirsh_weights(s,bas,idb,w)
       use systemmod, only: system
       use types, only: basindat
       type(system), intent(inout) :: s
       type(basindat), intent(in) :: bas
       integer, intent(in) :: idb
       real*8, intent(out) :: w(:,:,:)
     end subroutine hirsh_weights
     module function hirsh_rho(c,xn,i,r)
       use crystalmod, only: crystal
       type(crystal), intent(in) :: c
       real*8, intent(in) :: xn(:)
       integer, intent(in) :: i
       real*8, intent(in) :: r
       real*8 :: hirsh_rho
     end function hirsh_rho
     module subroutine hirsh_cutoffs(c,xn,rcut)
       use crystalmod, only: crystal
       type(crystal), intent(in) :: c
       real*8, intent(in) :: xn(:)
       real*8, allocatable, intent(inout) :: rcut(:,:)
     end subroutine hirsh_cutoffs
     module subroutine hirsh_option(word,line,lp,iter,tol,maxit,found)
       character*(*), intent(in) :: word
       character*(*), intent(in) :: line
       integer, intent(inout) :: lp
       logical, intent(inout) :: iter
       real*8, intent(inout) :: tol
       integer, intent(inout) :: maxit
       logical, intent(out) :: found
     end subroutine hirsh_option
     module subroutine voronoi_grid(s,bas)
       use systemmod, only: system
       use types, only: basindat
       type(system), intent(inout) :: s
       type(basindat), intent(inout) :: bas
     end subroutine voronoi_grid
     module subroutine hirsh_nogrid(line)
       character*(*), intent(in) :: line
     end subroutine hirsh_nogrid
  end interface

end module hirshfeld
