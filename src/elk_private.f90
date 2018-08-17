!> Part of the following code has been adapted from the elk distribution, 
!> version 1.3.2.
!> Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
!> Distributed under the terms of the GNU General Public License.

! J.L. Casals Sainz contributed the adaptation to elk 2.1.25

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

!> Interface to ELK densities.
module elk_private
  use environmod, only: environ
  implicit none

  private

  type elkwfn
     real*8 :: x2c(3,3)
     real*8 :: c2x(3,3)
     integer :: lmaxvr
     real*8, allocatable :: spr(:,:)
     real*8, allocatable :: spr_a(:)
     real*8, allocatable :: spr_b(:)
     integer, allocatable :: nrmt(:)
     integer :: ngvec
     real*8, allocatable :: vgc(:,:) 
     integer, allocatable :: igfft(:)
     integer :: ncel0
     real*8, allocatable :: xcel(:,:)
     integer, allocatable :: iesp(:)
     real*8, allocatable :: rhomt(:,:,:)
     complex*16, allocatable :: rhok(:)
     real*8, allocatable :: rmt(:)
     integer :: n(3)
     real*8 :: maxrmt = 0d0
     logical :: isealloc = .false.
     type(environ), pointer :: e
   contains
     procedure :: end => elkwfn_end !< Deallocate all data
     procedure :: rmt_atom !< RMT from the closest atom
     procedure :: read_out !< Read the elkwfn data from an elk's OUT file
     procedure :: rho2 !< Calculate the density and derivatives at a point
     procedure :: tolap !< Convert an elkwfn into its Laplacian
  end type elkwfn
  public :: elkwfn

  interface
     module subroutine elkwfn_end(f)
       class(elkwfn), intent(inout) :: f
     end subroutine elkwfn_end
     module function rmt_atom(f,x)
       class(elkwfn), intent(in) :: f
       real*8, intent(in) :: x(3)
       real*8 :: rmt_atom
     end function rmt_atom
     module subroutine read_out(f,env,file,file2,file3)
       class(elkwfn), intent(inout) :: f
       type(environ), intent(in), target :: env
       character*(*), intent(in) :: file, file2
       character*(*), intent(in), optional :: file3
     end subroutine read_out
     module subroutine rho2(f,vpl,nder,frho,gfrho,hfrho)
       class(elkwfn), intent(in) :: f
       real(8), intent(in) :: vpl(3)
       real(8), intent(out) :: frho, gfrho(3), hfrho(3,3)
       integer, intent(in) :: nder
     end subroutine rho2
     module subroutine tolap(f)
       class(elkwfn), intent(inout) :: f
     end subroutine tolap
  end interface

end module elk_private

