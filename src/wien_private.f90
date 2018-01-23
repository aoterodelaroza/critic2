!> This module contains routines  adapted (with permission) from the WIEN2k
!> source by P. Blaha, K. Schwarz et al. 

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

!> Interface to WIEN2k densities.
module wien_private
  implicit none

  private

  type wienwfn
     ! field quantities
     logical :: cnorm = .true.
     integer, allocatable :: lm(:,:,:)
     integer, allocatable :: lmmax(:)
     real*8, allocatable :: slm(:,:,:)
     real*8, allocatable :: sk(:)
     real*8, allocatable :: ski(:)
     real*8, allocatable :: tauk(:)
     real*8, allocatable :: tauki(:)
     integer :: lastind
     integer :: lastniz
     ! structural
     logical :: ishlat
     integer :: nat
     real*8, allocatable  :: rotloc(:,:,:)
     real*8, allocatable :: rnot(:)
     real*8, allocatable :: rmt(:)
     real*8, allocatable :: dx(:)
     integer, allocatable :: jri(:)
     integer, allocatable :: multw(:)
     real*8 :: br1(3,3)
     real*8 :: br2(3,3)
     real*8 :: br3(3,3)
     logical :: ortho
     integer :: ndat
     integer, allocatable :: iatnr(:)
     real*8, allocatable :: pos(:,:)
     integer, allocatable :: iop(:)
     integer :: niord
     integer, allocatable :: iz(:,:,:)
     real*8, allocatable :: tau(:,:)
     integer :: npos 
     real*8  :: atp(3,343)
     integer :: nwav
     real*8, allocatable :: krec(:,:)
     integer, allocatable :: inst(:)
     logical :: cmpl
     real*8 :: a(3)
   contains
     procedure :: end => wien_end
     procedure :: rmt_atom
     procedure :: read_clmsum
     procedure :: rho2
     procedure :: tolap
  end type wienwfn
  public :: wienwfn
  
  interface
     module subroutine wien_end(f)
       class(wienwfn), intent(inout) :: f
     end subroutine wien_end
     module function rmt_atom(f,x)
       class(wienwfn), intent(in) :: f
       real*8, intent(in) :: x(3)
       real*8 :: rmt_atom
     end function rmt_atom
     module subroutine read_clmsum(f,file,file2)
       class(wienwfn), intent(inout) :: f
       character*(*), intent(in) :: file, file2
     end subroutine read_clmsum
     module subroutine rho2(f,v0,rho,grad,h)
       class(wienwfn), intent(in) :: f
       real*8, dimension(3), intent(in) :: v0
       real*8, intent(out) :: rho, grad(3), h(3,3)
     end subroutine rho2
     module subroutine tolap(f)
       class(wienwfn), intent(inout) :: f
     end subroutine tolap
  end interface

end module wien_private

