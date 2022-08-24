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

!> Interface to ELK densities and structures.
module elk_private
  use environmod, only: environ
  use types, only: thread_info
  implicit none

  private

  ! elk wavefunction type
  type elkwfn
     real*8 :: x2c(3,3) ! cryst to cart matrix (elk's may be different from critic's)
     integer :: lmaxvr ! max LM in muffin tins
     real*8, allocatable :: spr(:,:) ! radial grid
     real*8, allocatable :: spr_a(:) ! radial coef - prefactor
     real*8, allocatable :: spr_b(:) ! radial coef - exponent
     integer, allocatable :: nrmt(:) ! number of radial points
     integer :: ngvec ! number of G-vectors
     real*8, allocatable :: vgc(:,:) ! G-vector coordinates
     integer, allocatable :: igfft(:) ! G-vector mapping
     real*8, allocatable :: rhomt(:,:,:) ! muffin tin coeffs
     complex*16, allocatable :: rhok(:) ! interstitial rho coeffs
     real*8, allocatable :: rmt(:) ! muffin tin radii
     integer :: n(3) ! interstitial grid size
     type(environ), pointer :: e ! pointer to the environment
   contains
     procedure :: end => elkwfn_end !< Deallocate all data
     procedure :: read_out !< Read the elkwfn data from an elk's OUT file
     procedure :: rho2 !< Calculate the density and derivatives at a point
     procedure :: tolap !< Convert an elkwfn into its Laplacian
  end type elkwfn
  public :: elkwfn

  interface
     module subroutine elkwfn_end(f)
       class(elkwfn), intent(inout) :: f
     end subroutine elkwfn_end
     module subroutine read_out(f,env,file,file2,file3,ti)
       class(elkwfn), intent(inout) :: f
       type(environ), intent(in), target :: env
       character*(*), intent(in) :: file, file2
       character*(*), intent(in), optional :: file3
       type(thread_info), intent(in), optional :: ti
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

