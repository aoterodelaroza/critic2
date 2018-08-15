!> Part of the following code has been adapted from postg
!> Copyright (c) 2015 Alberto Otero de la Roza <aoterodelaroza@ucmerced.edu>,
!> Felix Kannemann <felix.kannemann@dal.ca>, Erin R. Johnson <ejohnson29@ucmerced.edu>,
!> Ross M. Dickson <ross.dickson@dal.ca>, Hartmut Schmider <hs7@post.queensu.ca>,
!> and Axel D. Becke <axel.becke@dal.ca>

! Copyright (c) 2009-2017 Alberto Otero de la Roza <aoterodelaroza@gmail.com>,
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

! molecular wavefunction readers and tools
module wfn_private
  implicit none
  
  private

  ! Molecular wavefunction type
  ! Order of the orbitals:
  ! Restricted wavefunctions (wfntyp = wfn_rhf):
  !
  !     alpha and beta           alpha and beta
  ! | ..................... | ..................... |
  ! 1                     nmoocc                  nmoall
  !
  ! Unrestricted wavefunctions (wfntyp = wfn_uhf):
  !
  !     alpha       beta        alpha       beta
  ! | ......... | ......... | ......... | ......... |
  ! 1         nalpha      nmoocc     nmoocc+      nmoall
  !                                nalpha_virt
  ! 
  ! Fractional wavefunctions (wfntyp = wfn_frac): 
  ! 
  ! | ............................................. |
  ! 1                                         nmoall = nmoocc
  ! 
  ! In rhf and uhf wavefunctions, the virtual orbitals (nmoocc to
  ! nmoall) are only available if hasvirtual = .true.
  type molwfn
     logical :: useecp !< this wavefunction comes from a calc with ECPs
     integer :: nmoall !< number of molecular orbitals (occupied and virtual)
     integer :: nmoocc !< number of occupied molecular orbitals
     integer :: nalpha !< number of alpha electrons
     integer :: nalpha_virt !< number of virtual alpha orbitals
     integer :: npri !< number of primitives
     integer :: wfntyp !< type of wavefunction (rhf, uhf, fractional occ)
     logical :: issto !< are the primitives GTOs or STOs?
     logical :: hasvirtual !< are the virtual orbitals known?
     integer :: ixmaxsto(4) !< maximum exponent for x, y, z, and r in STOs
     integer, allocatable :: icenter(:) !< primitive center
     integer, allocatable :: itype(:) !< primitive type (see li(:,:) array)
     real*8, allocatable :: d2ran(:) !< maximum d^2 (GTO) or d (STO) to discard the primitive
     real*8, allocatable :: e(:) !< primitive exponents
     real*8, allocatable :: occ(:) !< MO occupation numbers
     real*8, allocatable :: cmo(:,:) !< MO coefficients
     integer :: nedf !< number of EDFs (electron density functions - core density for ECPs)
     integer, allocatable :: icenter_edf(:) !< EDF centers
     integer, allocatable :: itype_edf(:) !< EDF types
     real*8, allocatable :: e_edf(:) !< EDF exponents
     real*8, allocatable :: c_edf(:) !< EDF coefficients
     ! atomic positions, internal copy
     integer :: nat !< number of atoms
     real*8, allocatable :: xat(:,:) !< atomic coordinates (Cartesian, bohr)
   contains
     procedure :: end => wfn_end !< deallocate all arrays in wfn object
     procedure :: read_wfn !< read wavefunction from a Gaussian wfn file
     procedure :: read_wfx !< read wavefunction from a Gaussian wfx file
     procedure :: read_fchk !< read wavefunction from a Gaussian formatted checkpoint file
     procedure :: read_molden !< read wavefunction from a molden file
     procedure :: register_struct !< pass the atomic number and positions to the object
     procedure :: rho2 !< calculate the density, derivatives, and other properties
     procedure :: calculate_mo !< calculate the MO values at a point (driver)
     procedure :: calculate_mo_sto !< calculate the MO values at a point (STO version)
     procedure :: calculate_mo_gto !< calculate the MO values at a point (GTO version)
  end type molwfn
  public :: molwfn

  public :: wfn_read_xyz_geometry
  public :: wfn_read_wfn_geometry
  public :: wfn_read_wfx_geometry
  public :: wfn_read_fchk_geometry
  public :: wfn_read_molden_geometry
  public :: wfn_read_log_geometry

  ! wfn type identifier
  integer, parameter, public :: wfn_rhf = 0
  integer, parameter, public :: wfn_uhf = 1
  integer, parameter, public :: wfn_frac = 2
  
  interface
     module subroutine wfn_end(f)
       class(molwfn), intent(inout) :: f
     end subroutine wfn_end
     module subroutine wfn_read_xyz_geometry(file,n,x,z,name,errmsg)
       character*(*), intent(in) :: file
       integer, intent(out) :: n
       real*8, allocatable, intent(inout) :: x(:,:)
       integer, allocatable, intent(inout) :: z(:)
       character*(10), allocatable, intent(inout) :: name(:)
       character(len=:), allocatable, intent(out) :: errmsg
     end subroutine wfn_read_xyz_geometry
     module subroutine wfn_read_wfn_geometry(file,n,x,z,name,errmsg)
       character*(*), intent(in) :: file
       integer, intent(out) :: n
       real*8, allocatable, intent(inout) :: x(:,:)
       integer, allocatable, intent(inout) :: z(:)
       character*(10), allocatable, intent(inout) :: name(:)
       character(len=:), allocatable, intent(out) :: errmsg
     end subroutine wfn_read_wfn_geometry
     module subroutine wfn_read_wfx_geometry(file,n,x,z,name,errmsg)
       character*(*), intent(in) :: file
       integer, intent(out) :: n
       real*8, allocatable, intent(inout) :: x(:,:)
       integer, allocatable, intent(inout) :: z(:)
       character*(10), allocatable, intent(inout) :: name(:)
       character(len=:), allocatable, intent(out) :: errmsg
     end subroutine wfn_read_wfx_geometry
     module subroutine wfn_read_fchk_geometry(file,n,x,z,name,errmsg)
       character*(*), intent(in) :: file
       integer, intent(out) :: n
       real*8, allocatable, intent(inout) :: x(:,:)
       integer, allocatable, intent(inout) :: z(:)
       character*(10), allocatable, intent(inout) :: name(:)
       character(len=:), allocatable, intent(out) :: errmsg
     end subroutine wfn_read_fchk_geometry
     module subroutine wfn_read_molden_geometry(file,n,x,z,name,errmsg)
       character*(*), intent(in) :: file
       integer, intent(out) :: n
       real*8, allocatable, intent(inout) :: x(:,:)
       integer, allocatable, intent(inout) :: z(:)
       character*(10), allocatable, intent(inout) :: name(:)
       character(len=:), allocatable, intent(out) :: errmsg
     end subroutine wfn_read_molden_geometry
     module subroutine wfn_read_log_geometry(file,n,x,z,name,errmsg)
       character*(*), intent(in) :: file
       integer, intent(out) :: n
       real*8, allocatable, intent(inout) :: x(:,:)
       integer, allocatable, intent(inout) :: z(:)
       character*(10), allocatable, intent(inout) :: name(:)
       character(len=:), allocatable, intent(out) :: errmsg
     end subroutine wfn_read_log_geometry
     module subroutine read_wfn(f,file)
       class(molwfn), intent(inout) :: f
       character*(*), intent(in) :: file
     end subroutine read_wfn
     module subroutine read_wfx(f,file)
       class(molwfn), intent(inout) :: f
       character*(*), intent(in) :: file
     end subroutine read_wfx
     module subroutine read_fchk(f,file,readvirtual)
       class(molwfn), intent(inout) :: f
       character*(*), intent(in) :: file
       logical, intent(in) :: readvirtual
     end subroutine read_fchk
     module subroutine read_molden(f,file,readvirtual)
       class(molwfn), intent(inout) :: f
       character*(*), intent(in) :: file
       logical, intent(in) :: readvirtual
     end subroutine read_molden
     module subroutine register_struct(f,ncel,atcel)
       use types, only: celatom
       class(molwfn), intent(inout) :: f
       integer, intent(in) :: ncel
       type(celatom), intent(in) :: atcel(:)
     end subroutine register_struct
     module subroutine rho2(f,xpos,nder,rho,grad,h,gkin,vir,stress,xmo)
       use tools_io, only: ferror, faterr
       class(molwfn), intent(in) :: f
       real*8, intent(in) :: xpos(3)
       integer, intent(in) :: nder 
       real*8, intent(out) :: rho
       real*8, intent(out) :: grad(3)
       real*8, intent(out) :: h(3,3)
       real*8, intent(out) :: gkin
       real*8, intent(out) :: vir
       real*8, intent(out) :: stress(3,3)
       real*8, allocatable, intent(out), optional :: xmo(:)
     end subroutine rho2
     module subroutine calculate_mo(f,xpos,phi,fder)
       class(molwfn), intent(in) :: f
       real*8, intent(in) :: xpos(3)
       real*8, intent(out) :: phi
       character*(*), intent(in) :: fder
     end subroutine calculate_mo
     module subroutine calculate_mo_sto(f,xpos,phi,philb,imo0,imo1,nder)
       class(molwfn), intent(in) :: f
       real*8, intent(in) :: xpos(3)
       real*8, intent(inout) :: phi(:,:)
       integer, intent(in) :: philb
       integer, intent(in) :: imo0
       integer, intent(in) :: imo1
       integer, intent(in) :: nder
     end subroutine calculate_mo_sto
     module subroutine calculate_mo_gto(f,xpos,phi,philb,imo0,imo1,nder)
       class(molwfn), intent(in) :: f
       real*8, intent(in) :: xpos(3)
       real*8, intent(inout) :: phi(:,:)
       integer, intent(in) :: philb
       integer, intent(in) :: imo0
       integer, intent(in) :: imo1
       integer, intent(in) :: nder
     end subroutine calculate_mo_gto
  end interface

end module wfn_private
