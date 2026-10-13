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

!> 3d plotting of gradient paths.
module flux
  use types, only: gpathp
  implicit none

  private

  !> A gradient path to trace: its starting point, the direction of
  !> the trace and, if it leaves a critical point, the CP (the first
  !> point of the path).
  type flux_seed
     real*8 :: x0(3) = 0d0 ! starting point (Cartesian, bohr)
     real*8 :: xorig(3) = 0d0 ! origin of the path (Cartesian, bohr): the CP if haso, else x0
     logical :: haso = .false. ! whether the path leaves the CP at xorig
     real*8 :: dir(3) = 0d0 ! starting direction (Cartesian unit vector, zero if none)
     integer :: iup = 1 ! direction of the path (+1 up, -1 down)
  end type flux_seed
  public :: flux_seed

  !> A traced gradient path: its seed, the points (p, crystallographic
  !> coordinates in %x continuous along the path) and the cell CP at the
  !> end of the path (icpend, 0 = the path does not end at a CP) and its
  !> distance to the last point (dend, bohr).
  type flux_path
     type(flux_seed) :: seed ! the seed of this path
     integer :: ier = 0 ! termination code of the gradient routine
     integer :: icpend = 0 ! cell CP at the end of the path (0 = none)
     real*8 :: dend = 0d0 ! distance between the end of the path and that CP (bohr)
     type(gpathp), allocatable :: p(:) ! the points of the path
  end type flux_path
  public :: flux_path

  ! methods for the starting angles around a BCP/RCP (flux_seeds_axis)
  integer, parameter, public :: flux_method_braindead = 0 ! uniform angles
  integer, parameter, public :: flux_method_quotient = 1 ! remapped with the eigenvalue quotient
  integer, parameter, public :: flux_method_dynamical = 2 ! linearized dynamical system

  !xx! list of routines
  public  :: fluxprint
  public  :: flux_seeds_point
  public  :: flux_seeds_sphere
  public  :: flux_seeds_axis
  public  :: flux_trace

  interface
     module subroutine fluxprint()
     end subroutine fluxprint
     module subroutine flux_seeds_point(x,iup,seeds)
       real*8, intent(in) :: x(3)
       integer, intent(in) :: iup
       type(flux_seed), allocatable, intent(inout) :: seeds(:)
     end subroutine flux_seeds_point
     module subroutine flux_seeds_sphere(x,iup,ntheta,nphi,seeds)
       real*8, intent(in) :: x(3)
       integer, intent(in) :: iup, ntheta, nphi
       type(flux_seed), allocatable, intent(inout) :: seeds(:)
     end subroutine flux_seeds_sphere
     module subroutine flux_seeds_axis(f,x,ircp,iup,npts,method,seeds,prune,errmsg)
       use fieldmod, only: field
       type(field), intent(inout) :: f
       real*8, intent(in) :: x(3)
       integer, intent(in) :: ircp, iup, npts, method
       type(flux_seed), allocatable, intent(inout) :: seeds(:)
       real*8, intent(in) :: prune
       character(len=:), allocatable, intent(out) :: errmsg
     end subroutine flux_seeds_axis
     module subroutine flux_trace(f,seeds,paths,prune,abortable)
       use fieldmod, only: field
       type(field), intent(inout) :: f
       type(flux_seed), intent(in) :: seeds(:)
       type(flux_path), allocatable, intent(inout) :: paths(:)
       real*8, intent(in) :: prune
       logical, intent(in), optional :: abortable
     end subroutine flux_trace
  end interface

end module flux
