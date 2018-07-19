!> Parts of code taken from pi7, by VLC, EFM, AMP, MFA, MBV, MAB.
!> (c) Victor Lua~na Cabal, Universidad de Oviedo, 1987--

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

!> Interface to aiPI (pi7) densities.
module pi_private
  use grid1mod, only: grid1
  implicit none

  private

  type piwfn
     ! atomic orbitals
     character*(20), allocatable :: piname(:)
     logical, allocatable :: pi_used(:)
     integer, allocatable :: nsym(:)
     integer, allocatable :: naos(:,:)
     integer, allocatable :: naaos(:,:)
     integer, allocatable :: nsto(:,:)
     integer, allocatable :: nasto(:,:)
     integer, allocatable :: nn(:,:)
     real*8, allocatable :: z(:,:)
     real*8, allocatable :: xnsto(:,:)
     real*8, allocatable :: c(:,:,:)
     real*8, allocatable :: nelec(:,:)
     type(grid1), allocatable :: pgrid(:)
     ! structural info
     real*8, allocatable :: renv(:,:)
     integer, allocatable :: idxenv(:)
     integer, allocatable :: zenv(:)
     integer :: nenv
   contains
     procedure :: end => pi_end
     procedure :: register_struct
     procedure :: read_ion
     procedure :: rho2
     procedure :: fillinterpol
  end type piwfn
  public :: piwfn

  interface
     module subroutine pi_end(f)
       class(piwfn), intent(inout) :: f
     end subroutine pi_end
     module subroutine register_struct(f,nenv,spc,atenv)
       use types, only: anyatom, species
       class(piwfn), intent(inout) :: f
       integer, intent(in) :: nenv
       type(species), intent(in) :: spc(*)
       type(anyatom), intent(in) :: atenv(nenv)
     end subroutine register_struct
     module subroutine read_ion(f,fichero,ni)
       class(piwfn), intent(inout) :: f
       character*(*) :: fichero
       integer :: ni
     end subroutine read_ion
     module subroutine rho2(f,xpos,exact,rho,grad,h)
       class(piwfn), intent(in) :: f
       real*8, intent(in) :: xpos(3)
       logical, intent(in) :: exact
       real*8, intent(out) :: rho, grad(3), h(3,3)
     end subroutine rho2
     module subroutine fillinterpol(f)
       class(piwfn), intent(inout) :: f
     end subroutine fillinterpol
  end interface

end module pi_private

