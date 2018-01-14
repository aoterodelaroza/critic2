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
     integer :: nenv
     real*8, allocatable :: renv(:,:)
     integer, allocatable :: lenv(:,:)
     integer, allocatable :: idxenv(:)
     integer, allocatable :: zenv(:)
     real*8 :: globalcutoff = 0d0
   contains
     procedure :: end => dftb_end
     procedure :: read => dftb_read
     procedure :: rho2
     procedure :: register_struct
  end type dftbwfn
  public :: dftbwfn

  private :: next_logical
  private :: next_integer
  private :: read_kpointsandweights
  private :: read_occupations
  private :: dftb_read_reals1
  private :: next_hsd_atom
  private :: build_interpolation_grid1
  private :: calculate_rl
  
  interface
     module subroutine dftb_end(f)
       class(dftbwfn), intent(inout) :: f
     end subroutine dftb_end
     module subroutine dftb_read(f,filexml,filebin,filehsd,atcel,spc)
       use types, only: celatom, species
       class(dftbwfn), intent(inout) :: f
       character*(*), intent(in) :: filexml
       character*(*), intent(in) :: filebin
       character*(*), intent(in) :: filehsd
       type(celatom), intent(in) :: atcel(:)
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
     module subroutine register_struct(f,rmat,atenv,spc)
       use types, only: celatom, species
       class(dftbwfn), intent(inout) :: f
       real*8, intent(in) :: rmat(3,3)
       type(celatom), intent(in) :: atenv(:)
       type(species), intent(in) :: spc(:)
     end subroutine register_struct
     module function next_logical(lu,line0,key0) result(next)
       integer, intent(in) :: lu
       character*(*), intent(in) :: line0, key0
       logical :: next
     end function next_logical
     module function next_integer(lu,line0,key0) result(next)
       integer, intent(in) :: lu
       character*(*), intent(in) :: line0, key0
       integer :: next
     end function next_integer
     module subroutine read_kpointsandweights(lu,kpts,w)
       integer, intent(in) :: lu
       real*8, intent(out) :: kpts(:,:)
       real*8, intent(out) :: w(:)
     end subroutine read_kpointsandweights
     module subroutine read_occupations(lu,occ)
       integer, intent(in) :: lu
       real*8, intent(out) :: occ(:,:,:)
     end subroutine read_occupations
     module function dftb_read_reals1(lu,n) result(x)
       integer, intent(in) :: lu, n
       real*8 :: x(n)
     endfunction dftb_read_reals1
     module function next_hsd_atom(lu,at) result(ok)
       integer, intent(in) :: lu
       type(dftbatom), intent(out) :: at
       logical :: ok
     end function next_hsd_atom
     module subroutine build_interpolation_grid1(ff)
       class(dftbwfn), intent(inout) :: ff
     end subroutine build_interpolation_grid1
     module subroutine calculate_rl(ff,it,iorb,r0,f,fp,fpp)
       class(dftbwfn), intent(in) :: ff
       integer, intent(in) :: it, iorb
       real*8, intent(in) :: r0
       real*8, intent(out) :: f, fp, fpp
     end subroutine calculate_rl
     module subroutine realloc_dftbatom(a,nnew)
       type(dftbatom), intent(inout), allocatable :: a(:)
       integer, intent(in) :: nnew
     end subroutine realloc_dftbatom
  end interface

end module dftb_private
