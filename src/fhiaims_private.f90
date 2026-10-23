! Copyright (c) 2016-2022 Alberto Otero de la Roza <aoterodelaroza@gmail.com>,
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

! FHI-aims wavefunctions, read from the dgrid_aims.dat file (output dgrid).
module fhiaims_private
  use iso_c_binding, only: c_ptr
  use types, only: thread_info
  implicit none

  private

  !> Header of a dgrid_aims.dat file (first three lines). The file
  !> layout is: header; grid sizes per species; atoms; two text lines
  !> (energy, chemical potential); radial functions; spline
  !> coefficients; eigenvalues and occupations; eigenvectors; and,
  !> for crystals, the lattice vectors, the reciprocal vectors and
  !> the k-points with their weights.
  type dgrid_header
     integer :: nspc !< number of species
     integer :: nat !< number of atoms
     integer :: nstates !< number of states
     integer :: nbas !< number of basis functions
     integer :: nfn !< number of radial functions (all species)
     integer :: nper !< number of periodic dimensions (0 or 3)
     integer :: nspin !< number of spin channels
     integer :: nkpt !< number of k-points
     integer :: nkxyz(3) !< k-point grid
  end type dgrid_header
  public :: dgrid_header

  !> Radial function of the FHI-aims basis set: phi = u(r)/r * S_lm,
  !> where u is a cubic spline in the index of the logarithmic grid
  !> of its species, r_i = rmin * exp((i-1)*linc).
  type fhiaimsradial
     integer :: l !< angular momentum
     real*8 :: rcut !< outer radius (bohr), the function is zero beyond it
     real*8, allocatable :: spl(:,:) !< spline coefficients (4,ngrid)
  end type fhiaimsradial

  !> Species of the FHI-aims calculation
  type fhiaimsspecies
     character*20 :: name !< species name
     integer :: ngrid !< number of points in the logarithmic grid
     real*8 :: rmin !< first point of the logarithmic grid (bohr)
     real*8 :: linc !< log of the grid increment
     real*8 :: rcut !< maximum outer radius over the radial functions
     integer :: lmax !< maximum angular momentum
     integer :: nfn !< number of radial functions
     integer :: nbas !< number of basis functions (sum over 2l+1)
     type(fhiaimsradial), allocatable :: fn(:) !< radial functions, in file order
  end type fhiaimsspecies

  !> FHI-aims wavefunction, stored as a block-sparse real-space
  !> density matrix (summed over spins and k-points):
  !>   rho(r) = sum_{aT,bT'} phi_aT(r) D_ab(T-T') phi_bT'(r)
  !> for basis functions on atoms a, b (FHI-aims order) in cells T,
  !> T'. Only one of each pair of transposed blocks D_ab(dT) =
  !> D_ba(-dT)^T is stored: those with a < b, or a == b and dT
  !> lexicographically non-negative. The blocks of a pair (a,b) are
  !> contiguous in the block list, in increasing order of their
  !> lattice-vector key (see find_block).
  type fhiaimswfn
     logical :: isreal !< real eigenvectors in the file
     logical :: ismol !< molecular calculation (no periodicity)
     integer :: nspin !< number of spin channels
     integer :: nkpt !< number of k-points
     integer :: nstates !< number of states
     integer :: nbas !< number of basis functions
     integer :: nat !< number of atoms
     integer :: nspc !< number of species
     integer :: lmax !< maximum angular momentum over all species
     type(fhiaimsspecies), allocatable :: spc(:) !< species
     integer, allocatable :: ispc(:) !< species of each FHI-aims atom
     ! mapping to the critic2 cell atoms
     integer, allocatable :: iaims(:) !< FHI-aims atom for each critic2 cell atom
     integer, allocatable :: ioff(:,:) !< aims position = critic2 position + ioff (lattice)
     real*8, allocatable :: spcutoff(:,:) !< cutoff for list_near_atoms (c%nspc,2)
     ! density matrix blocks
     integer :: nblk !< number of blocks
     integer :: tmax(3) !< lattice vectors of the blocks run from -tmax to tmax
     integer, allocatable :: pair0(:,:) !< first block of pair (a,b), a <= b
     integer, allocatable :: pairn(:,:) !< number of blocks of pair (a,b), a <= b
     integer, allocatable :: blkt(:,:) !< lattice vector dT of each block (3,nblk)
     integer, allocatable :: blkkey(:) !< integer key of dT, increasing within a pair
     integer*8, allocatable :: blkoff(:) !< offset of each block in dm
     real*8, allocatable :: dm(:) !< density matrix blocks, column-major (nbas_a,nbas_b)
     type(c_ptr) :: cptr !< pointer to the crystal structure
   contains
     procedure :: end => fhiaims_end
     procedure :: read => fhiaims_read
     procedure :: rho2
  end type fhiaimswfn
  public :: fhiaimswfn

  public :: dgrid_is_file
  public :: dgrid_read_header
  public :: dgrid_read_atoms
  public :: dgrid_read_tail

  interface
     module function dgrid_is_file(file,ti) result(ok)
       character*(*), intent(in) :: file
       type(thread_info), intent(in), optional :: ti
       logical :: ok
     end function dgrid_is_file
     module subroutine dgrid_read_header(lu,hdr,errmsg)
       integer, intent(in) :: lu
       type(dgrid_header), intent(out) :: hdr
       character(len=:), allocatable, intent(out) :: errmsg
     end subroutine dgrid_read_header
     module subroutine dgrid_read_atoms(lu,hdr,ngrid,ispc,xat,spcname,errmsg)
       integer, intent(in) :: lu
       type(dgrid_header), intent(in) :: hdr
       integer, allocatable, intent(inout) :: ngrid(:)
       integer, allocatable, intent(inout) :: ispc(:)
       real*8, allocatable, intent(inout) :: xat(:,:)
       character*20, allocatable, intent(inout) :: spcname(:)
       character(len=:), allocatable, intent(out) :: errmsg
     end subroutine dgrid_read_atoms
     module subroutine dgrid_read_tail(file,nkpt,lat,kpt,wk,errmsg)
       character*(*), intent(in) :: file
       integer, intent(in) :: nkpt
       real*8, intent(out) :: lat(3,3)
       real*8, intent(out) :: kpt(3,nkpt)
       real*8, intent(out) :: wk(nkpt)
       character(len=:), allocatable, intent(out) :: errmsg
     end subroutine dgrid_read_tail
     module subroutine fhiaims_end(f)
       class(fhiaimswfn), intent(inout) :: f
     end subroutine fhiaims_end
     module subroutine fhiaims_read(f,cptr,file,errmsg,ti)
       class(fhiaimswfn), intent(inout) :: f
       type(c_ptr), intent(in) :: cptr
       character*(*), intent(in) :: file
       character(len=:), allocatable, intent(out) :: errmsg
       type(thread_info), intent(in), optional :: ti
     end subroutine fhiaims_read
     module subroutine rho2(f,xpos,nder,rho,grad,h,gkin)
       class(fhiaimswfn), intent(inout) :: f
       real*8, intent(in) :: xpos(3)
       integer, intent(in) :: nder
       real*8, intent(out) :: rho
       real*8, intent(out) :: grad(3)
       real*8, intent(out) :: h(3,3)
       real*8, intent(out) :: gkin
     end subroutine rho2
  end interface

end module fhiaims_private
