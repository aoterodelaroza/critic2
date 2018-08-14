! Copyright (c) 2016 Alberto Otero de la Roza <aoterodelaroza@gmail.com>,
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

! Interface to DFTB+ wavefunctions.
module dftb_private
  use grid1mod, only: grid1
  use environmod, only: environ
  implicit none
  
  private
  
  !> atomic basis set information in dftb fields
  type dftbatom
     character*10 :: name  !< name of the atom
     integer :: z !< atomic number
     integer :: norb !< number of orbitals
     integer :: l(4) !< orbital angular momentum
     real*8 :: occ(4) !< orbital occupations
     real*8 :: cutoff(4) !< orbital cutoffs
     integer :: nexp(4) !< number of exponents 
     real*8 :: eexp(5,4) !< exponents
     integer :: ncoef(5,4) !< number of coefficients
     real*8 :: coef(5,5,4) !< coefficients (icoef,iexp,iorb)
     type(grid1), allocatable :: orb(:) !< radial component of the orbitals (grid)
  end type dftbatom
  
  type dftbwfn
     logical :: isreal
     integer :: nkpt
     integer :: nspin
     integer :: nstates
     integer :: norb
     real*8, allocatable :: docc(:,:,:)
     real*8, allocatable :: dkpt(:,:)
     real*8, allocatable :: evecr(:,:,:)
     complex*16, allocatable :: evecc(:,:,:,:)
     integer, allocatable :: ispec(:)
     integer, allocatable :: idxorb(:)
     integer :: midxorb
     integer :: maxnorb
     integer :: maxlm
     type(dftbatom), allocatable :: bas(:)
     ! structural info
     real*8 :: globalcutoff = 0d0
     logical :: isealloc = .false.
     type(environ), pointer :: e
   contains
     procedure :: end => dftb_end
     procedure :: read => dftb_read
     procedure :: rho2
     procedure :: register_struct
  end type dftbwfn
  public :: dftbwfn

  interface
     module subroutine dftb_end(f)
       class(dftbwfn), intent(inout) :: f
     end subroutine dftb_end
     module subroutine dftb_read(f,filexml,filebin,filehsd,atcel,spc)
       use types, only: anyatom, species
       class(dftbwfn), intent(inout) :: f
       character*(*), intent(in) :: filexml
       character*(*), intent(in) :: filebin
       character*(*), intent(in) :: filehsd
       class(anyatom), intent(in) :: atcel(:)
       type(species), intent(in) :: spc(:)
     end subroutine dftb_read
     module subroutine rho2(f,xpos,exact,nder,rho,grad,h,gkin)
       class(dftbwfn), intent(inout) :: f
       real*8, intent(in) :: xpos(3)
       logical, intent(in) :: exact
       integer, intent(in) :: nder 
       real*8, intent(out) :: rho
       real*8, intent(out) :: grad(3)
       real*8, intent(out) :: h(3,3)
       real*8, intent(out) :: gkin
     end subroutine rho2
     module subroutine register_struct(f,e)
       use types, only: anyatom, species
       class(dftbwfn), intent(inout) :: f
       type(environ), intent(in), target :: e
     end subroutine register_struct
  end interface

end module dftb_private
