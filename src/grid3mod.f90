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

! Class for 3d grids and related tools.
module grid3mod
  use iso_c_binding, only: c_ptr
  use types, only: thread_info
  use param, only: mlen
  implicit none

  private

  ! grid interpolation modes
  integer, parameter, public :: mode_nearest = 1 !< nearest grid node
  integer, parameter, public :: mode_trilinear = 2 !< trilinear
  integer, parameter, public :: mode_trispline = 3 !< trispline
  integer, parameter, public :: mode_tricubic = 4 !< tricubic
  integer, parameter, public :: mode_smr = 5 !< smoothrho
  integer, parameter, public :: mode_default = mode_tricubic

  !> Information for QE Kohn-Sham states, maybe plus wannier functions
  type qedat
     integer :: nks !< Number of k-points (nk1*nk2*nk3)
     integer :: nk(3) !< Number of k-points in each grid direction (0 if not a grid)
     integer :: nbnd !< Number of bands
     integer :: nspin !< Number of spins
     logical :: gamma_only !< Whether gamma-only tricks are being used
     character(len=mlen) :: fpwc !< pwc file name
     real*8, allocatable :: kpt(:,:) !< k-points in cryst. coords. (3,k)
     real*8, allocatable :: wk(:) !< k-point weights (k)
     real*8, allocatable :: ek(:,:,:) !< band energies (bnd,k,spin)
     real*8, allocatable :: occ(:,:,:) !< band occupations (bnd,k,spin)
     integer, allocatable :: ngk(:) !< number of plane-waves for each k-point (k)
     integer, allocatable :: igk_k(:,:) !< fft reorder (npw,nk)
     integer, allocatable :: nl(:) !< fft reorder (ngms)
     integer, allocatable :: nlm(:) !< fft reorder (ngms) (gamma only)
     ! Wannier transformation
     integer :: nbndw(2) !< Number of occuped bands (for the Wannier transform)
     real*8, allocatable :: center(:,:,:) !< wannier function centers (cryst)
     real*8, allocatable :: spread(:,:) !< wannier function spreads (bohr)
     complex*16, allocatable :: u(:,:,:,:) !< u matrix
  end type qedat
  public :: qedat

  !> Three-dimensional grid class
  type grid3
     ! initialization and info flags
     logical :: isinit = .false. !< is the grid initialized?
     logical :: isqe = .false. !< does it have qe ks states info?
     logical :: iswan = .false. !< does it have wannier info?
     ! basic grid info
     integer :: mode = mode_default !< interpolation mode
     integer :: n(3) !< number of grid points in each direction
     type(c_ptr) :: cptr !< pointer to the crystal structure
     ! information for partial grids (that do not span the whole cell)
     logical :: partial = .false. !< true if the grid only spans a part of the cell
     real*8 :: x0(3) !< origin of the grid (cryst. coords.)
     real*8 :: x2cl(3,3) !< local x2c: lattice vectors for the grid domain
     real*8 :: c2xl(3,3) !< grid Car. to cryst. matrix, inverse of x2cl
     real*8 :: x2cg(3,3) !< grid step vectors: x2cl / n
     ! geometry of the crystal (also of the grid, if partial = .false.)
     real*8 :: x2c(3,3) !< crystallographic to Cartesian matrix (crystal)
     real*8 :: c2x(3,3) !< Cartesian to crystallographic matrix (crystal)
     ! grid point neighbor topology
     real*8 :: dmax !< largest grid step
     integer :: nvec !< number of neighbor grid points
     integer, allocatable :: vec(:,:) !< grid coordinates of neighbor grid points
     real*8, allocatable :: area(:) !< area of the Voronoi facets
     ! the actual grid data
     real*8, allocatable :: f(:,:,:) !< grid values
     ! trispline interpolation
     real*8, allocatable :: c2(:,:,:,:) !< cubic coefficients
     ! smoothrho interpolation
     integer :: smr_nlist ! number of nodes in the stencil
     integer, allocatable :: smr_ilist(:,:) ! integer offsets of the stencil nodes
     real*8, allocatable :: smr_xlist(:,:) ! positions of the stencil nodes
     real*8, allocatable :: smr_phiinv(:,:) ! inverse of the phi matrix
     real*8, allocatable :: smr_rho0(:,:,:) ! promolecular density on the grid (for smoothing)
     integer :: smr_nenv ! target number of nodes in the stencil
     real*8 :: smr_fdmax ! dmax factor for the continuity smoothing
     ! QE band states and Wannier function transformation
     type(qedat) :: qe
   contains
     procedure :: end => grid_end !< deallocate all arrays and uninitialize
     procedure :: setmode !< set the interpolation mode of a grid
     procedure :: normalize !< normalize the grid to a given value
     procedure :: from_array3 !< build a grid3 from a 3d array of real numbers
     procedure :: read_cube !< grid3 from a Gaussian cube file
     procedure :: read_bincube !< grid3 from a binary cube file
     procedure :: read_siesta !< grid3 from siesta RHO file
     procedure :: read_fplo !< grid3 from FPLO 001 file
     procedure :: read_abinit !< grid3 from abinit binary file
     procedure :: read_vasp !< grid3 from VASP file (CHG, CHGCAR, etc.)
     procedure :: read_qub !< grid3 from aimpac qub format
     procedure :: read_xsf !< grid3 from xsf (xcrysden) file
     procedure :: read_fmt !< grid3 from fmt (CASTEP) file
     procedure :: read_txt !< grid3 from txt (BAND) file
     procedure :: read_pwc !< read a pwc file created by pw2critic.x
     procedure :: read_elk !< grid3 from elk file format
     procedure :: read_wannier_chk !< qe/wannier info from chk file
     procedure :: interp !< interpolate the grid at an arbitrary point
     procedure :: fft !< grid3 as the FFT of another grid3
     procedure :: resample !< grid3 as a Fourier resampling of another grid3
     procedure :: rotate_qe_evc !< write U-rotated scratch files using QE evc file
     procedure :: get_qe_wnr !< build a Wannier function from Bloch coeffs (pre-open, parallel)
     procedure :: get_qe_wnr_standalone !< build a Wannier function from Bloch coeffs (standalone)
     procedure :: get_qe_psink_standalone !< build a psink(r) function from Bloch coeffs (standalone)
     procedure :: new_eval !< grid3 from an arithmetic expression
  end type grid3
  public :: grid3

  interface
     module subroutine new_eval(f,sptr,cptr,n,expr,x2c)
       use iso_c_binding, only: c_ptr
       class(grid3), intent(inout) :: f
       type(c_ptr), intent(in) :: sptr
       type(c_ptr), intent(in) :: cptr
       integer, intent(in) :: n(3)
       character(*), intent(in) :: expr
       real*8, intent(in) :: x2c(3,3)
     end subroutine new_eval
     module subroutine grid_end(f)
       class(grid3), intent(inout) :: f
     end subroutine grid_end
     module subroutine setmode(f,mode)
       class(grid3), intent(inout) :: f
       character*(*), intent(in) :: mode
     end subroutine setmode
     module subroutine normalize(f,norm)
       class(grid3), intent(inout) :: f
       real*8, intent(in) :: norm
     end subroutine normalize
     module subroutine from_array3(f,g,x2c,cptr)
       class(grid3), intent(inout) :: f
       real*8, intent(in) :: g(:,:,:)
       real*8, intent(in) :: x2c(3,3)
       type(c_ptr), intent(in) :: cptr
     end subroutine from_array3
     module subroutine read_cube(f,cptr,file,x2c,molx0,errmsg,ti)
       class(grid3), intent(inout) :: f
       type(c_ptr), intent(in) :: cptr
       character*(*), intent(in) :: file
       real*8, intent(in) :: x2c(3,3), molx0(3)
       character(len=:), allocatable, intent(out) :: errmsg
       type(thread_info), intent(in), optional :: ti
     end subroutine read_cube
     module subroutine read_bincube(f,cptr,file,x2c,molx0,errmsg,ti)
       class(grid3), intent(inout) :: f
       type(c_ptr), intent(in) :: cptr
       character*(*), intent(in) :: file
       real*8, intent(in) :: x2c(3,3), molx0(3)
       character(len=:), allocatable, intent(out) :: errmsg
       type(thread_info), intent(in), optional :: ti
     end subroutine read_bincube
     module subroutine read_siesta(f,cptr,file,x2c,errmsg,ti)
       class(grid3), intent(inout) :: f
       type(c_ptr), intent(in) :: cptr
       character*(*), intent(in) :: file
       real*8, intent(in) :: x2c(3,3)
       character(len=:), allocatable, intent(out) :: errmsg
       type(thread_info), intent(in), optional :: ti
     end subroutine read_siesta
     module subroutine read_fplo(f,cptr,file,x2c,errmsg,ti)
       class(grid3), intent(inout) :: f
       type(c_ptr), intent(in) :: cptr
       character*(*), intent(in) :: file
       real*8, intent(in) :: x2c(3,3)
       character(len=:), allocatable, intent(out) :: errmsg
       type(thread_info), intent(in), optional :: ti
     end subroutine read_fplo
     module subroutine read_abinit(f,cptr,file,x2c,errmsg,ti)
       class(grid3), intent(inout) :: f
       type(c_ptr), intent(in) :: cptr
       character*(*), intent(in) :: file
       real*8, intent(in) :: x2c(3,3)
       character(len=:), allocatable, intent(out) :: errmsg
       type(thread_info), intent(in), optional :: ti
     end subroutine read_abinit
     module subroutine read_vasp(f,cptr,file,x2c,vscal,ibl,errmsg,ti)
       class(grid3), intent(inout) :: f
       type(c_ptr), intent(in) :: cptr
       character*(*), intent(in) :: file
       real*8, intent(in) :: x2c(3,3)
       logical, intent(in) :: vscal
       integer, intent(in), optional :: ibl
       character(len=:), allocatable, intent(out) :: errmsg
       type(thread_info), intent(in), optional :: ti
     end subroutine read_vasp
     module subroutine read_qub(f,cptr,file,x2c,errmsg,ti)
       class(grid3), intent(inout) :: f
       type(c_ptr), intent(in) :: cptr
       character*(*), intent(in) :: file
       real*8, intent(in) :: x2c(3,3)
       character(len=:), allocatable, intent(out) :: errmsg
       type(thread_info), intent(in), optional :: ti
     end subroutine read_qub
     module subroutine read_xsf(f,cptr,file,x2c,errmsg,ti)
       class(grid3), intent(inout) :: f
       type(c_ptr), intent(in) :: cptr
       character*(*), intent(in) :: file
       real*8, intent(in) :: x2c(3,3)
       character(len=:), allocatable, intent(out) :: errmsg
       type(thread_info), intent(in), optional :: ti
     end subroutine read_xsf
     module subroutine read_fmt(f,cptr,file,x2c,errmsg,ti)
       class(grid3), intent(inout) :: f
       type(c_ptr), intent(in) :: cptr
       character*(*), intent(in) :: file
       real*8, intent(in) :: x2c(3,3)
       character(len=:), allocatable, intent(out) :: errmsg
       type(thread_info), intent(in), optional :: ti
     end subroutine read_fmt
     module subroutine read_txt(f,cptr,file,x2c,errmsg,ti)
       class(grid3), intent(inout) :: f
       type(c_ptr), intent(in) :: cptr
       character*(*), intent(in) :: file
       real*8, intent(in) :: x2c(3,3)
       character(len=:), allocatable, intent(out) :: errmsg
       type(thread_info), intent(in), optional :: ti
     end subroutine read_txt
     module subroutine read_pwc(f,cptr,fpwc,ispin,ikpt,ibnd,emin,emax,x2c,errmsg,ti)
       class(grid3), intent(inout) :: f
       type(c_ptr), intent(in) :: cptr
       character*(*), intent(in) :: fpwc
       integer, intent(in) :: ispin
       integer, intent(in), allocatable :: ikpt(:)
       integer, intent(in), allocatable :: ibnd(:)
       real*8, intent(in) :: emin, emax
       real*8, intent(in) :: x2c(3,3)
       character(len=:), allocatable, intent(out) :: errmsg
       type(thread_info), intent(in), optional :: ti
     end subroutine read_pwc
     module subroutine read_elk(f,cptr,file,x2c,errmsg,ti)
       class(grid3), intent(inout) :: f
       type(c_ptr), intent(in) :: cptr
       character*(*), intent(in) :: file
       real*8, intent(in) :: x2c(3,3)
       character(len=:), allocatable, intent(out) :: errmsg
       type(thread_info), intent(in), optional :: ti
     end subroutine read_elk
     module subroutine read_wannier_chk(f,fileup,filedn,errmsg,ti)
       class(grid3), intent(inout) :: f
       character*(*), intent(in) :: fileup
       character*(*), intent(in), optional :: filedn
       character(len=:), allocatable, intent(out), optional :: errmsg
       type(thread_info), intent(in), optional :: ti
     end subroutine read_wannier_chk
     module subroutine interp(f,xi,y,yp,ypp,valid)
       class(grid3), intent(inout) :: f !< Grid to interpolate
       real*8, intent(in) :: xi(3) !< Target point (cryst. coords.)
       real*8, intent(out) :: y !< Interpolated value
       real*8, intent(out) :: yp(3) !< First derivative
       real*8, intent(out) :: ypp(3,3) !< Second derivative
       logical, intent(out) :: valid !< whether the point was in the grid domain
     end subroutine interp
     module subroutine fft(fnew,fold,iff)
       class(grid3), intent(inout) :: fnew
       type(grid3), intent(in) :: fold
       integer, intent(in) :: iff
     end subroutine fft
     module subroutine resample(frs,frho,n2)
       class(grid3), intent(inout) :: frs
       type(grid3), intent(in) :: frho
       integer, intent(in) :: n2(3)
     end subroutine resample
     module subroutine rotate_qe_evc(f,luevc,luevc_ibnd,useu,ti)
       class(grid3), intent(inout) :: f
       integer, intent(out) :: luevc(2)
       integer, intent(out) :: luevc_ibnd(2)
       logical, intent(in) :: useu
       type(thread_info), intent(in), optional :: ti
     end subroutine rotate_qe_evc
     module subroutine get_qe_wnr(f,omega,ibnd,ispin,luevc,luevc_ibnd,fout)
       class(grid3), intent(in) :: f
       real*8, intent(in) :: omega
       integer, intent(in) :: ibnd
       integer, intent(in) :: ispin
       integer, intent(in) :: luevc(2)
       integer, intent(inout) :: luevc_ibnd(2)
       complex*16, intent(out) :: fout(:,:,:,:)
     end subroutine get_qe_wnr
     module subroutine get_qe_wnr_standalone(f,omega,ibnd,ispin,inr,rotate,fout,ti)
       class(grid3), intent(in) :: f
       real*8, intent(in) :: omega
       integer, intent(in) :: ibnd
       integer, intent(in) :: ispin
       integer, intent(in) :: inr(3)
       logical, intent(in) :: rotate
       complex*16, intent(out) :: fout(:,:,:)
       type(thread_info), intent(in), optional :: ti
     end subroutine get_qe_wnr_standalone
     module subroutine get_qe_psink_standalone(f,omega,ibnd,ik,ispin,usephase,inr,fout,ti)
       class(grid3), intent(in) :: f
       real*8, intent(in) :: omega
       integer, intent(in) :: ibnd
       integer, intent(in) :: ik
       integer, intent(in) :: ispin
       logical :: usephase
       integer, intent(in) :: inr(3)
       complex*16, intent(out) :: fout(:,:,:)
       type(thread_info), intent(in), optional :: ti
     end subroutine get_qe_psink_standalone
  end interface

end module grid3mod
