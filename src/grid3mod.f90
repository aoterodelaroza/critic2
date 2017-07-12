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

! Class for 3d grids and related tools.
module grid3mod
  implicit none

  private

  private :: grinterp_nearest
  private :: grinterp_trilinear
  private :: grinterp_trispline
  private :: grinterp_tricubic
  private :: grid_near
  private :: grid_floor
  private :: init_trispline
  private :: pop_grid
  
  !> Information for wannier functions
  type wandat
     integer :: nks !< Number of k-points (lattice vectors)
     integer :: nwan(3) !< Number of lattice vectors
     integer :: nbnd !< Number of bands
     integer :: nspin !< Number of spins
     logical :: sijchk = .true. !< Save the sij to a checkpoint (chk-sij) file
     logical :: fachk = .true. !< Save the fa to a checkpoint (chk-fa) file
     logical :: haschk !< A chk-sij file is currently available
     logical :: useu = .true. !< Use the U transformation to get MLWF
     real*8 :: cutoff !< Cutoff for atomic overlaps
     character*(255) :: fevc !< evc file name
     real*8, allocatable :: kpt(:,:) !< k-points in fract. coords.
     real*8, allocatable :: center(:,:,:) !< wannier function centers (cryst)
     real*8, allocatable :: spread(:,:) !< wannier function spreads (bohr)
     integer, allocatable :: ngk(:) !< number of plane-waves for each k-point
     integer, allocatable :: igk_k(:,:) !< fft reorder
     integer, allocatable :: nls(:) !< fft reorder
     complex*16, allocatable :: u(:,:,:) !< u matrix
  end type wandat

  !> Three-dimensional grid class
  type grid3
     logical :: isinit = .false. !< is the grid initialized?
     logical :: iswan = .false. !< does it have wannier info?
     integer :: mode !< interpolation mode
     integer :: n(3) !< number of grid points in each direction
     real*8, allocatable :: f(:,:,:) !< grid values
     real*8, allocatable :: c2(:,:,:,:) !< cubic coefficients for spline interpolation
     type(wandat) :: wan !< Wannier functions and related information
   contains
     procedure :: end => grid_end !< deallocate all arrays and uninitialize
     procedure :: setmode !< set the interpolation mode of a grid
     procedure :: normalize !< normalize the grid to a given value
     procedure :: from_array3 !< build a grid3 from a 3d array of real numbers
     procedure :: read_cube !< grid3 from a Gaussian cube file
     procedure :: read_siesta !< grid3 from siesta RHO file
     procedure :: read_abinit !< grid3 from abinit binary file
     procedure :: read_vasp !< grid3 from VASP file (CHG, CHGCAR, etc.)
     procedure :: read_qub !< grid3 from aimpac qub format
     procedure :: read_xsf !< grid3 from xsf (xcrysden) file
     procedure :: read_unk !< read a UNK file by QE for wannier90 
     procedure :: read_unkgen !< read a unkgen file created by pw2wannier.x
     procedure :: read_elk !< grid3 from elk file format
     procedure :: interp !< interpolate the grid at an arbitrary point
     procedure :: laplacian !< grid3 as the Laplacian of another grid3
     procedure :: gradrho !< grid3 as the gradrho of another grid3
     procedure :: hxx !< grid3 as the Hessian components of another grid3
     procedure :: get_qe_wnr !< build a Wannier function from the wannier info/files
     procedure :: new_eval !< grid3 from an arithmetic expression
  end type grid3
  public :: grid3

  integer, parameter :: mode_nearest = 1 !< interpolation mode: nearest grid node
  integer, parameter :: mode_trilinear = 2 !< interpolation mode: trilinear
  integer, parameter :: mode_trispline = 3 !< interpolation mode: trispline
  integer, parameter :: mode_tricubic = 4 !< interpolation mode: tricubic
  integer, parameter :: mode_default = mode_tricubic

  ! The 64x64 matrix for tricubic interpolation
  real*8, parameter, private :: c(64,64) = reshape((/&                      ! values for c(i,j), with...  (i,  j)
    1d0,  0d0, -3d0,  2d0,  0d0,  0d0,  0d0,  0d0, -3d0,  0d0,   9d0,  -6d0,  2d0,  0d0,  -6d0,   4d0,&  ! 1-16, 1
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&	 ! 17-32, 1
   -3d0,  0d0,  9d0, -6d0,  0d0,  0d0,  0d0,  0d0,  9d0,  0d0, -27d0,  18d0, -6d0,  0d0,  18d0, -12d0,&	 ! 33-48, 1
    2d0,  0d0, -6d0,  4d0,  0d0,  0d0,  0d0,  0d0, -6d0,  0d0,  18d0, -12d0,  4d0,  0d0, -12d0,   8d0,&	 ! 49-64, 1
    0d0,  0d0,  3d0, -2d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  -9d0,   6d0,  0d0,  0d0,   6d0,  -4d0,&	 ! 1-16, 2
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&	 ! 17-32, 2
    0d0,  0d0, -9d0,  6d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  27d0, -18d0,  0d0,  0d0, -18d0,  12d0,&	 ! 33-48, 2
    0d0,  0d0,  6d0, -4d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0, -18d0,  12d0,  0d0,  0d0,  12d0,  -8d0,&	 ! 49-64, 2
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  3d0,  0d0,  -9d0,   6d0, -2d0,  0d0,   6d0,  -4d0,&	 ! 1-16, 3
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&	 ! 17-32, 3
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0, -9d0,  0d0,  27d0, -18d0,  6d0,  0d0, -18d0,  12d0,&	 ! 33-48, 3
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  6d0,  0d0, -18d0,  12d0, -4d0,  0d0,  12d0,  -8d0,&	 ! 49-64, 3
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   9d0,  -6d0,  0d0,  0d0,  -6d0,   4d0,&	 ! 1-16, 4
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&	 ! 17-32, 4
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0, -27d0,  18d0,  0d0,  0d0,  18d0, -12d0,&	 ! 33-48, 4
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  18d0, -12d0,  0d0,  0d0, -12d0,   8d0,&	 ! 49-64, 4
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&	 ! 1-16, 5
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&	 ! 17-32, 5
    3d0,  0d0, -9d0,  6d0,  0d0,  0d0,  0d0,  0d0, -9d0,  0d0,  27d0, -18d0,  6d0,  0d0, -18d0,  12d0,&	 ! 33-48, 5
   -2d0,  0d0,  6d0, -4d0,  0d0,  0d0,  0d0,  0d0,  6d0,  0d0, -18d0,  12d0, -4d0,  0d0,  12d0,  -8d0,&	 ! 49-64, 5
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&	 ! 1-16, 6
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&	 ! 17-32, 6
    0d0,  0d0,  9d0, -6d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0, -27d0,  18d0,  0d0,  0d0,  18d0, -12d0,&	 ! 33-48, 6
    0d0,  0d0, -6d0,  4d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  18d0, -12d0,  0d0,  0d0, -12d0,   8d0,&	 ! 49-64, 6
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&	 ! 1-16, 7
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&	 ! 17-32, 7
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  9d0,  0d0, -27d0,  18d0, -6d0,  0d0,  18d0, -12d0,&	 ! 33-48, 7
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0, -6d0,  0d0,  18d0, -12d0,  4d0,  0d0, -12d0,   8d0,&	 ! 49-64, 7
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&	 ! 1-16, 8
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&	 ! 17-32, 8
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  27d0, -18d0,  0d0,  0d0, -18d0,  12d0,&	 ! 33-48, 8
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0, -18d0,  12d0,  0d0,  0d0,  12d0,  -8d0,&	 ! 49-64, 8
    0d0,  1d0, -2d0,  1d0,  0d0,  0d0,  0d0,  0d0,  0d0, -3d0,   6d0,  -3d0,  0d0,  2d0,  -4d0,   2d0,&	 ! 1-16, 9
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&	 ! 17-32, 9
    0d0, -3d0,  6d0, -3d0,  0d0,  0d0,  0d0,  0d0,  0d0,  9d0, -18d0,   9d0,  0d0, -6d0,  12d0,  -6d0,&	 ! 33-48, 9
    0d0,  2d0, -4d0,  2d0,  0d0,  0d0,  0d0,  0d0,  0d0, -6d0,  12d0,  -6d0,  0d0,  4d0,  -8d0,   4d0,&	 ! 49-64, 9
    0d0,  0d0, -1d0,  1d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   3d0,  -3d0,  0d0,  0d0,  -2d0,   2d0,&	 ! 1-16, 10
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&	 ! 17-32, 10
    0d0,  0d0,  3d0, -3d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  -9d0,   9d0,  0d0,  0d0,   6d0,  -6d0,&	 ! 33-48, 10
    0d0,  0d0, -2d0,  2d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   6d0,  -6d0,  0d0,  0d0,  -4d0,   4d0,&	 ! 49-64, 10
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  3d0,  -6d0,   3d0,  0d0, -2d0,   4d0,  -2d0,&	 ! 1-16, 11
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&	 ! 17-32, 11
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0, -9d0,  18d0,  -9d0,  0d0,  6d0, -12d0,   6d0,&	 ! 33-48, 11
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  6d0, -12d0,   6d0,  0d0, -4d0,   8d0,  -4d0,&	 ! 49-64, 11
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  -3d0,   3d0,  0d0,  0d0,   2d0,  -2d0,&	 ! 1-16, 12
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&	 ! 17-32, 12
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   9d0,  -9d0,  0d0,  0d0,  -6d0,   6d0,&	 ! 33-48, 12
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  -6d0,   6d0,  0d0,  0d0,   4d0,  -4d0,&	 ! 49-64, 12
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&	 ! 1-16, 13
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&	 ! 17-32, 13
    0d0,  3d0, -6d0,  3d0,  0d0,  0d0,  0d0,  0d0,  0d0, -9d0,  18d0,  -9d0,  0d0,  6d0, -12d0,   6d0,&	 ! 33-48, 13
    0d0, -2d0,  4d0, -2d0,  0d0,  0d0,  0d0,  0d0,  0d0,  6d0, -12d0,   6d0,  0d0, -4d0,   8d0,  -4d0,&	 ! 49-64, 13
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&	 ! 1-16, 14
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&	 ! 17-32, 14
    0d0,  0d0, -3d0,  3d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   9d0,  -9d0,  0d0,  0d0,  -6d0,   6d0,&	 ! 33-48, 14
    0d0,  0d0,  2d0, -2d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  -6d0,   6d0,  0d0,  0d0,   4d0,  -4d0,&	 ! 49-64, 14
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&	 ! 1-16, 15
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&	 ! 17-32, 15
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  9d0, -18d0,   9d0,  0d0, -6d0,  12d0,  -6d0,&	 ! 33-48, 15
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0, -6d0,  12d0,  -6d0,  0d0,  4d0,  -8d0,   4d0,&	 ! 49-64, 15
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&	 ! 1-16, 16
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&	 ! 17-32, 16
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  -9d0,   9d0,  0d0,  0d0,   6d0,  -6d0,&	 ! 33-48, 16
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   6d0,  -6d0,  0d0,  0d0,  -4d0,   4d0,&	 ! 49-64, 16
    0d0,  0d0,  0d0,  0d0,  1d0,  0d0, -3d0,  2d0, -2d0,  0d0,   6d0,  -4d0,  1d0,  0d0,  -3d0,   2d0,&	 ! 1-16, 17
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&	 ! 17-32, 17
    0d0,  0d0,  0d0,  0d0, -3d0,  0d0,  9d0, -6d0,  6d0,  0d0, -18d0,  12d0, -3d0,  0d0,   9d0,  -6d0,&	 ! 33-48, 17
    0d0,  0d0,  0d0,  0d0,  2d0,  0d0, -6d0,  4d0, -4d0,  0d0,  12d0,  -8d0,  2d0,  0d0,  -6d0,   4d0,&	 ! 49-64, 17
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  3d0, -2d0,  0d0,  0d0,  -6d0,   4d0,  0d0,  0d0,   3d0,  -2d0,&	 ! 1-16, 18
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&	 ! 17-32, 18
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0, -9d0,  6d0,  0d0,  0d0,  18d0, -12d0,  0d0,  0d0,  -9d0,   6d0,&	 ! 33-48, 18
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  6d0, -4d0,  0d0,  0d0, -12d0,   8d0,  0d0,  0d0,   6d0,  -4d0,&	 ! 49-64, 18
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0, -1d0,  0d0,   3d0,  -2d0,  1d0,  0d0,  -3d0,   2d0,&	 ! 1-16, 19
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&	 ! 17-32, 19
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  3d0,  0d0,  -9d0,   6d0, -3d0,  0d0,   9d0,  -6d0,&	 ! 33-48, 19
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0, -2d0,  0d0,   6d0,  -4d0,  2d0,  0d0,  -6d0,   4d0,&	 ! 49-64, 19
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  -3d0,   2d0,  0d0,  0d0,   3d0,  -2d0,&	 ! 1-16, 20
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&	 ! 17-32, 20
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   9d0,  -6d0,  0d0,  0d0,  -9d0,   6d0,&	 ! 33-48, 20
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  -6d0,   4d0,  0d0,  0d0,   6d0,  -4d0,&	 ! 49-64, 20
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&	 ! 1-16, 21
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&	 ! 17-32, 21
    0d0,  0d0,  0d0,  0d0,  3d0,  0d0, -9d0,  6d0, -6d0,  0d0,  18d0, -12d0,  3d0,  0d0,  -9d0,   6d0,&	 ! 33-48, 21
    0d0,  0d0,  0d0,  0d0, -2d0,  0d0,  6d0, -4d0,  4d0,  0d0, -12d0,   8d0, -2d0,  0d0,   6d0,  -4d0,&	 ! 49-64, 21
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&	 ! 1-16, 22
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&	 ! 17-32, 22
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  9d0, -6d0,  0d0,  0d0, -18d0,  12d0,  0d0,  0d0,   9d0,  -6d0,&	 ! 33-48, 22
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0, -6d0,  4d0,  0d0,  0d0,  12d0,  -8d0,  0d0,  0d0,  -6d0,   4d0,&	 ! 49-64, 22
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&	 ! 1-16, 23
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&	 ! 17-32, 23
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0, -3d0,  0d0,   9d0,  -6d0,  3d0,  0d0,  -9d0,   6d0,&	 ! 33-48, 23
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  2d0,  0d0,  -6d0,   4d0, -2d0,  0d0,   6d0,  -4d0,&	 ! 49-64, 23
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&	 ! 1-16, 24
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&	 ! 17-32, 24
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  -9d0,   6d0,  0d0,  0d0,   9d0,  -6d0,&	 ! 33-48, 24
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   6d0,  -4d0,  0d0,  0d0,  -6d0,   4d0,&	 ! 49-64, 24
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&	 ! 1-16, 25
    1d0,  0d0, -3d0,  2d0,  0d0,  0d0,  0d0,  0d0, -3d0,  0d0,   9d0,  -6d0,  2d0,  0d0,  -6d0,   4d0,&	 ! 17-32, 25
   -2d0,  0d0,  6d0, -4d0,  0d0,  0d0,  0d0,  0d0,  6d0,  0d0, -18d0,  12d0, -4d0,  0d0,  12d0,  -8d0,&	 ! 33-48, 25
    1d0,  0d0, -3d0,  2d0,  0d0,  0d0,  0d0,  0d0, -3d0,  0d0,   9d0,  -6d0,  2d0,  0d0,  -6d0,   4d0,&	 ! 49-64, 25
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&	 ! 1-16, 26
    0d0,  0d0,  3d0, -2d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  -9d0,   6d0,  0d0,  0d0,   6d0,  -4d0,&	 ! 17-32, 26
    0d0,  0d0, -6d0,  4d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  18d0, -12d0,  0d0,  0d0, -12d0,   8d0,&	 ! 33-48, 26
    0d0,  0d0,  3d0, -2d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  -9d0,   6d0,  0d0,  0d0,   6d0,  -4d0,&	 ! 49-64, 26
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&	 ! 1-16, 27
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  3d0,  0d0,  -9d0,   6d0, -2d0,  0d0,   6d0,  -4d0,&	 ! 17-32, 27
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0, -6d0,  0d0,  18d0, -12d0,  4d0,  0d0, -12d0,   8d0,&	 ! 33-48, 27
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  3d0,  0d0,  -9d0,   6d0, -2d0,  0d0,   6d0,  -4d0,&	 ! 49-64, 27
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&	 ! 1-16, 28
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   9d0,  -6d0,  0d0,  0d0,  -6d0,   4d0,&	 ! 17-32, 28
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0, -18d0,  12d0,  0d0,  0d0,  12d0,  -8d0,&	 ! 33-48, 28
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   9d0,  -6d0,  0d0,  0d0,  -6d0,   4d0,&	 ! 49-64, 28
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&	 ! 1-16, 29
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&	 ! 17-32, 29
   -1d0,  0d0,  3d0, -2d0,  0d0,  0d0,  0d0,  0d0,  3d0,  0d0,  -9d0,   6d0, -2d0,  0d0,   6d0,  -4d0,&	 ! 33-48, 29
    1d0,  0d0, -3d0,  2d0,  0d0,  0d0,  0d0,  0d0, -3d0,  0d0,   9d0,  -6d0,  2d0,  0d0,  -6d0,   4d0,&	 ! 49-64, 29
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&	 ! 1-16, 30
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&	 ! 17-32, 30
    0d0,  0d0, -3d0,  2d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   9d0,  -6d0,  0d0,  0d0,  -6d0,   4d0,&	 ! 33-48, 30
    0d0,  0d0,  3d0, -2d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  -9d0,   6d0,  0d0,  0d0,   6d0,  -4d0,&	 ! 49-64, 30
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&	 ! 1-16, 31
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&	 ! 17-32, 31
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0, -3d0,  0d0,   9d0,  -6d0,  2d0,  0d0,  -6d0,   4d0,&	 ! 33-48, 31
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  3d0,  0d0,  -9d0,   6d0, -2d0,  0d0,   6d0,  -4d0,&	 ! 49-64, 31
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&	 ! 1-16, 32
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&	 ! 17-32, 32
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  -9d0,   6d0,  0d0,  0d0,   6d0,  -4d0,&	 ! 33-48, 32
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   9d0,  -6d0,  0d0,  0d0,  -6d0,   4d0,&	 ! 49-64, 32
    0d0,  0d0,  0d0,  0d0,  0d0,  1d0, -2d0,  1d0,  0d0, -2d0,   4d0,  -2d0,  0d0,  1d0,  -2d0,   1d0,&	 ! 1-16, 33
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&	 ! 17-32, 33
    0d0,  0d0,  0d0,  0d0,  0d0, -3d0,  6d0, -3d0,  0d0,  6d0, -12d0,   6d0,  0d0, -3d0,   6d0,  -3d0,&	 ! 33-48, 33
    0d0,  0d0,  0d0,  0d0,  0d0,  2d0, -4d0,  2d0,  0d0, -4d0,   8d0,  -4d0,  0d0,  2d0,  -4d0,   2d0,&	 ! 49-64, 33
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0, -1d0,  1d0,  0d0,  0d0,   2d0,  -2d0,  0d0,  0d0,  -1d0,   1d0,&	 ! 1-16, 34
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&	 ! 17-32, 34
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  3d0, -3d0,  0d0,  0d0,  -6d0,   6d0,  0d0,  0d0,   3d0,  -3d0,&	 ! 33-48, 34
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0, -2d0,  2d0,  0d0,  0d0,   4d0,  -4d0,  0d0,  0d0,  -2d0,   2d0,&	 ! 49-64, 34
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0, -1d0,   2d0,  -1d0,  0d0,  1d0,  -2d0,   1d0,&	 ! 1-16, 35
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&	 ! 17-32, 35
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  3d0,  -6d0,   3d0,  0d0, -3d0,   6d0,  -3d0,&	 ! 33-48, 35
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0, -2d0,   4d0,  -2d0,  0d0,  2d0,  -4d0,   2d0,&	 ! 49-64, 35
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   1d0,  -1d0,  0d0,  0d0,  -1d0,   1d0,&	 ! 1-16, 36
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&	 ! 17-32, 36
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  -3d0,   3d0,  0d0,  0d0,   3d0,  -3d0,&	 ! 33-48, 36
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   2d0,  -2d0,  0d0,  0d0,  -2d0,   2d0,&	 ! 49-64, 36
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&	 ! 1-16, 37
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&	 ! 17-32, 37
    0d0,  0d0,  0d0,  0d0,  0d0,  3d0, -6d0,  3d0,  0d0, -6d0,  12d0,  -6d0,  0d0,  3d0,  -6d0,   3d0,&	 ! 33-48, 37
    0d0,  0d0,  0d0,  0d0,  0d0, -2d0,  4d0, -2d0,  0d0,  4d0,  -8d0,   4d0,  0d0, -2d0,   4d0,  -2d0,&	 ! 49-64, 37
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&	 ! 1-16, 38
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&	 ! 17-32, 38
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0, -3d0,  3d0,  0d0,  0d0,   6d0,  -6d0,  0d0,  0d0,  -3d0,   3d0,&	 ! 33-48, 38
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  2d0, -2d0,  0d0,  0d0,  -4d0,   4d0,  0d0,  0d0,   2d0,  -2d0,&	 ! 49-64, 38
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&	 ! 1-16, 39
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&	 ! 17-32, 39
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0, -3d0,   6d0,  -3d0,  0d0,  3d0,  -6d0,   3d0,&	 ! 33-48, 39
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  2d0,  -4d0,   2d0,  0d0, -2d0,   4d0,  -2d0,&	 ! 49-64, 39
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&	 ! 1-16, 40
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&	 ! 17-32, 40
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   3d0,  -3d0,  0d0,  0d0,  -3d0,   3d0,&	 ! 33-48, 40
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  -2d0,   2d0,  0d0,  0d0,   2d0,  -2d0,&	 ! 49-64, 40
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&	 ! 1-16, 41
    0d0,  1d0, -2d0,  1d0,  0d0,  0d0,  0d0,  0d0,  0d0, -3d0,   6d0,  -3d0,  0d0,  2d0,  -4d0,   2d0,&	 ! 17-32, 41
    0d0, -2d0,  4d0, -2d0,  0d0,  0d0,  0d0,  0d0,  0d0,  6d0, -12d0,   6d0,  0d0, -4d0,   8d0,  -4d0,&	 ! 33-48, 41
    0d0,  1d0, -2d0,  1d0,  0d0,  0d0,  0d0,  0d0,  0d0, -3d0,   6d0,  -3d0,  0d0,  2d0,  -4d0,   2d0,&	 ! 49-64, 41
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&	 ! 1-16, 42
    0d0,  0d0, -1d0,  1d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   3d0,  -3d0,  0d0,  0d0,  -2d0,   2d0,&	 ! 17-32, 42
    0d0,  0d0,  2d0, -2d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  -6d0,   6d0,  0d0,  0d0,   4d0,  -4d0,&	 ! 33-48, 42
    0d0,  0d0, -1d0,  1d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   3d0,  -3d0,  0d0,  0d0,  -2d0,   2d0,&	 ! 49-64, 42
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&	 ! 1-16, 43
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  3d0,  -6d0,   3d0,  0d0, -2d0,   4d0,  -2d0,&	 ! 17-32, 43
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0, -6d0,  12d0,  -6d0,  0d0,  4d0,  -8d0,   4d0,&	 ! 33-48, 43
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  3d0,  -6d0,   3d0,  0d0, -2d0,   4d0,  -2d0,&	 ! 49-64, 43
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&	 ! 1-16, 44
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  -3d0,   3d0,  0d0,  0d0,   2d0,  -2d0,&	 ! 17-32, 44
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   6d0,  -6d0,  0d0,  0d0,  -4d0,   4d0,&	 ! 33-48, 44
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  -3d0,   3d0,  0d0,  0d0,   2d0,  -2d0,&	 ! 49-64, 44
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&	 ! 1-16, 45
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&	 ! 17-32, 45
    0d0, -1d0,  2d0, -1d0,  0d0,  0d0,  0d0,  0d0,  0d0,  3d0,  -6d0,   3d0,  0d0, -2d0,   4d0,  -2d0,&	 ! 33-48, 45
    0d0,  1d0, -2d0,  1d0,  0d0,  0d0,  0d0,  0d0,  0d0, -3d0,   6d0,  -3d0,  0d0,  2d0,  -4d0,   2d0,&	 ! 49-64, 45
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&	 ! 1-16, 46
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&	 ! 17-32, 46
    0d0,  0d0,  1d0, -1d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  -3d0,   3d0,  0d0,  0d0,   2d0,  -2d0,&	 ! 33-48, 46
    0d0,  0d0, -1d0,  1d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   3d0,  -3d0,  0d0,  0d0,  -2d0,   2d0,&	 ! 49-64, 46
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&	 ! 1-16, 47
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&	 ! 17-32, 47
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0, -3d0,   6d0,  -3d0,  0d0,  2d0,  -4d0,   2d0,&	 ! 33-48, 47
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  3d0,  -6d0,   3d0,  0d0, -2d0,   4d0,  -2d0,&	 ! 49-64, 47
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&	 ! 1-16, 48
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&	 ! 17-32, 48
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   3d0,  -3d0,  0d0,  0d0,  -2d0,   2d0,&	 ! 33-48, 48
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  -3d0,   3d0,  0d0,  0d0,   2d0,  -2d0,&	 ! 49-64, 48
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&	 ! 1-16, 49
    0d0,  0d0,  0d0,  0d0,  1d0,  0d0, -3d0,  2d0, -2d0,  0d0,   6d0,  -4d0,  1d0,  0d0,  -3d0,   2d0,&	 ! 17-32, 49
    0d0,  0d0,  0d0,  0d0, -2d0,  0d0,  6d0, -4d0,  4d0,  0d0, -12d0,   8d0, -2d0,  0d0,   6d0,  -4d0,&	 ! 33-48, 49
    0d0,  0d0,  0d0,  0d0,  1d0,  0d0, -3d0,  2d0, -2d0,  0d0,   6d0,  -4d0,  1d0,  0d0,  -3d0,   2d0,&	 ! 49-64, 49
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&	 ! 1-16, 50
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  3d0, -2d0,  0d0,  0d0,  -6d0,   4d0,  0d0,  0d0,   3d0,  -2d0,&	 ! 17-32, 50
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0, -6d0,  4d0,  0d0,  0d0,  12d0,  -8d0,  0d0,  0d0,  -6d0,   4d0,&	 ! 33-48, 50
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  3d0, -2d0,  0d0,  0d0,  -6d0,   4d0,  0d0,  0d0,   3d0,  -2d0,&	 ! 49-64, 50
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&	 ! 1-16, 51
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0, -1d0,  0d0,   3d0,  -2d0,  1d0,  0d0,  -3d0,   2d0,&	 ! 17-32, 51
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  2d0,  0d0,  -6d0,   4d0, -2d0,  0d0,   6d0,  -4d0,&	 ! 33-48, 51
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0, -1d0,  0d0,   3d0,  -2d0,  1d0,  0d0,  -3d0,   2d0,&	 ! 49-64, 51
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&	 ! 1-16, 52
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  -3d0,   2d0,  0d0,  0d0,   3d0,  -2d0,&	 ! 17-32, 52
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   6d0,  -4d0,  0d0,  0d0,  -6d0,   4d0,&	 ! 33-48, 52
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  -3d0,   2d0,  0d0,  0d0,   3d0,  -2d0,&	 ! 49-64, 52
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&	 ! 1-16, 53
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&	 ! 17-32, 53
    0d0,  0d0,  0d0,  0d0, -1d0,  0d0,  3d0, -2d0,  2d0,  0d0,  -6d0,   4d0, -1d0,  0d0,   3d0,  -2d0,&	 ! 33-48, 53
    0d0,  0d0,  0d0,  0d0,  1d0,  0d0, -3d0,  2d0, -2d0,  0d0,   6d0,  -4d0,  1d0,  0d0,  -3d0,   2d0,&	 ! 49-64, 53
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&	 ! 1-16, 54
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&	 ! 17-32, 54
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0, -3d0,  2d0,  0d0,  0d0,   6d0,  -4d0,  0d0,  0d0,  -3d0,   2d0,&	 ! 33-48, 54
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  3d0, -2d0,  0d0,  0d0,  -6d0,   4d0,  0d0,  0d0,   3d0,  -2d0,&	 ! 49-64, 54
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&	 ! 1-16, 55
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&	 ! 17-32, 55
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  1d0,  0d0,  -3d0,   2d0, -1d0,  0d0,   3d0,  -2d0,&	 ! 33-48, 55
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0, -1d0,  0d0,   3d0,  -2d0,  1d0,  0d0,  -3d0,   2d0,&	 ! 49-64, 55
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&	 ! 1-16, 56
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&	 ! 17-32, 56
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   3d0,  -2d0,  0d0,  0d0,  -3d0,   2d0,&	 ! 33-48, 56
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  -3d0,   2d0,  0d0,  0d0,   3d0,  -2d0,&	 ! 49-64, 56
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&	 ! 1-16, 57
    0d0,  0d0,  0d0,  0d0,  0d0,  1d0, -2d0,  1d0,  0d0, -2d0,   4d0,  -2d0,  0d0,  1d0,  -2d0,   1d0,&	 ! 17-32, 57
    0d0,  0d0,  0d0,  0d0,  0d0, -2d0,  4d0, -2d0,  0d0,  4d0,  -8d0,   4d0,  0d0, -2d0,   4d0,  -2d0,&	 ! 33-48, 57
    0d0,  0d0,  0d0,  0d0,  0d0,  1d0, -2d0,  1d0,  0d0, -2d0,   4d0,  -2d0,  0d0,  1d0,  -2d0,   1d0,&	 ! 49-64, 57
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&	 ! 1-16, 58
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0, -1d0,  1d0,  0d0,  0d0,   2d0,  -2d0,  0d0,  0d0,  -1d0,   1d0,&	 ! 17-32, 58
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  2d0, -2d0,  0d0,  0d0,  -4d0,   4d0,  0d0,  0d0,   2d0,  -2d0,&	 ! 33-48, 58
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0, -1d0,  1d0,  0d0,  0d0,   2d0,  -2d0,  0d0,  0d0,  -1d0,   1d0,&	 ! 49-64, 58
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&	 ! 1-16, 59
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0, -1d0,   2d0,  -1d0,  0d0,  1d0,  -2d0,   1d0,&	 ! 17-32, 59
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  2d0,  -4d0,   2d0,  0d0, -2d0,   4d0,  -2d0,&	 ! 33-48, 59
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0, -1d0,   2d0,  -1d0,  0d0,  1d0,  -2d0,   1d0,&	 ! 49-64, 59
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&	 ! 1-16, 60
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   1d0,  -1d0,  0d0,  0d0,  -1d0,   1d0,&	 ! 17-32, 60
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  -2d0,   2d0,  0d0,  0d0,   2d0,  -2d0,&	 ! 33-48, 60
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   1d0,  -1d0,  0d0,  0d0,  -1d0,   1d0,&	 ! 49-64, 60
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&	 ! 1-16, 61
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&	 ! 17-32, 61
    0d0,  0d0,  0d0,  0d0,  0d0, -1d0,  2d0, -1d0,  0d0,  2d0,  -4d0,   2d0,  0d0, -1d0,   2d0,  -1d0,&	 ! 33-48, 61
    0d0,  0d0,  0d0,  0d0,  0d0,  1d0, -2d0,  1d0,  0d0, -2d0,   4d0,  -2d0,  0d0,  1d0,  -2d0,   1d0,&	 ! 49-64, 61
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&	 ! 1-16, 62
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&	 ! 17-32, 62
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  1d0, -1d0,  0d0,  0d0,  -2d0,   2d0,  0d0,  0d0,   1d0,  -1d0,&	 ! 33-48, 62
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0, -1d0,  1d0,  0d0,  0d0,   2d0,  -2d0,  0d0,  0d0,  -1d0,   1d0,&	 ! 49-64, 62
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&	 ! 1-16, 63
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&	 ! 17-32, 63
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  1d0,  -2d0,   1d0,  0d0, -1d0,   2d0,  -1d0,&	 ! 33-48, 63
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0, -1d0,   2d0,  -1d0,  0d0,  1d0,  -2d0,   1d0,&	 ! 49-64, 63
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&	 ! 1-16, 64
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&	 ! 17-32, 64
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  -1d0,   1d0,  0d0,  0d0,   1d0,  -1d0,&	 ! 33-48, 64
    0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   1d0,  -1d0,  0d0,  0d0,  -1d0,   1d0&	 ! 49-64, 64
    /),shape(c))

contains

  subroutine grid_end(f)
    class(grid3), intent(inout) :: f

    f%isinit = .false.
    f%iswan = .false.
    if (allocated(f%f)) deallocate(f%f)
    if (allocated(f%c2)) deallocate(f%c2)
    if (allocated(f%wan%kpt)) deallocate(f%wan%kpt)
    if (allocated(f%wan%center)) deallocate(f%wan%center)
    if (allocated(f%wan%spread)) deallocate(f%wan%spread)
    if (allocated(f%wan%ngk)) deallocate(f%wan%ngk)
    if (allocated(f%wan%igk_k)) deallocate(f%wan%igk_k)
    if (allocated(f%wan%nls)) deallocate(f%wan%nls)
    if (allocated(f%wan%u)) deallocate(f%wan%u)
    
  end subroutine grid_end

  !> Set the interpolation mode. The possible modes arenearest,
  !> trilinear, trispline, and tricubic (lowercase).
  subroutine setmode(f,mode)
    use tools_io, only: equal, lower
    class(grid3), intent(inout) :: f
    character*(*), intent(in) :: mode
    
    character(len=:), allocatable :: lmode
    
    lmode = lower(mode)
    if (equal(lmode,'tricubic')) then
       f%mode = mode_tricubic
    else if (equal(lmode,'trispline')) then
       f%mode = mode_trispline
       if (allocated(f%c2)) deallocate(f%c2)
    else if (equal(lmode,'trilinear')) then
       f%mode = mode_trilinear
    else if (equal(lmode,'nearest')) then
       f%mode = mode_nearest
    else if (equal(lmode,'default')) then
       f%mode = mode_default
    end if

  end subroutine setmode

  !> Normalize the grid to a given value. omega is the cell volume.
  subroutine normalize(f,norm,omega)
    class(grid3), intent(inout) :: f
    real*8, intent(in) :: norm, omega

    f%f = f%f / (sum(f%f) * omega / real(product(f%n),8)) * norm
    if (allocated(f%c2)) deallocate(f%c2)
    
  end subroutine normalize

  !> Build a grid field from a three-dimensional array
  subroutine from_array3(f,g)
    use tools_io, only: ferror, faterr
    class(grid3), intent(inout) :: f
    real*8, intent(in) :: g(:,:,:)

    integer :: istat, n(3)

    call f%end()
    f%isinit = .true.
    f%iswan = .false.
    f%mode = mode_default
    n = ubound(g) - lbound(g) + 1
    f%n(:) = n
    allocate(f%f(n(1),n(2),n(3)),stat=istat)
    if (istat /= 0) &
       call ferror('from_array3','Error allocating grid',faterr)
    f%f = g

  end subroutine from_array3

  !> Read a grid in Gaussian CUBE format
  subroutine read_cube(f,file)
    use tools_io, only: fopen_read, ferror, faterr, fclose
    class(grid3), intent(inout) :: f
    character*(*), intent(in) :: file !< Input file

    integer :: luc
    integer :: nat
    integer :: istat, n(3), i, j, k

    call f%end()
    luc = fopen_read(file)

    read (luc,*) 
    read (luc,*) 
    read (luc,*,iostat=istat) nat
    if (istat /= 0) &
       call ferror('read_cube','Error reading nat',faterr,file)
    do i = 1, 3
       read (luc,*,iostat=istat) n(i)
       if (istat /= 0) &
          call ferror('read_cube','Error reading nx, ny, nz',faterr,file)
    end do
    do i = 1, nat
       read (luc,*)
    end do

    f%isinit = .true.
    f%iswan = .false.
    f%mode = mode_default
    f%n(:) = n
    allocate(f%f(n(1),n(2),n(3)),stat=istat)
    if (istat /= 0) &
       call ferror('read_cube','Error allocating grid',faterr,file)
    read(luc,*,iostat=istat) (((f%f(i,j,k),k=1,n(3)),j=1,n(2)),i=1,n(1))
    if (istat /= 0) &
       call ferror('read_cube','Error reading grid',faterr,file)

    call fclose(luc)

  end subroutine read_cube

  !> Read a grid in siesta RHO format
  subroutine read_siesta(f,file)
    use tools_io, only: fopen_read, faterr, ferror, fclose
    class(grid3), intent(inout) :: f
    character*(*), intent(in) :: file !< Input file

    integer :: luc, nspin, istat
    integer :: i, iy, iz, n(3)
    real*8 :: r(3,3)
    real*4, allocatable :: g(:)

    ! initialize
    call f%end()
    f%isinit = .true.
    f%mode = mode_default

    ! open file
    luc = fopen_read(file,'unformatted')

    ! assume unformatted
    read (luc) r
    read (luc) n, nspin
    f%n = n

    allocate(f%f(n(1),n(2),n(3)),stat=istat)
    if (istat /= 0) &
       call ferror('read_siesta','Error allocating grid',faterr,file)
    allocate(g(n(1)),stat=istat)
    if (istat /= 0) &
       call ferror('read_siesta','Error allocating auxiliary grid',faterr,file)
    f%f = 0d0
    do i = 1, nspin
       do iz = 1, n(3)
          do iy = 1, n(2)
             read (luc) g
             f%f(:,iy,iz) = f%f(:,iy,iz) + g
          end do
       end do
    end do
    deallocate(g)

    call fclose(luc)

  end subroutine read_siesta

  !> Read a grid in abinit format
  subroutine read_abinit(f,file)
    use tools_io, only: fopen_read, ferror, faterr, fclose
    use abinit_private, only: hdr_type, hdr_io
    class(grid3), intent(inout) :: f
    character*(*), intent(in) :: file !< Input file

    integer :: luc
    integer :: fform0, istat, n(3)
    type(hdr_type) :: hdr
    real*8, allocatable :: g(:,:,:)

    call f%end()
    luc = fopen_read(file,'unformatted')

    ! read the header
    call hdr_io(fform0,hdr,1,luc)

    f%isinit = .true.
    f%mode = mode_default
    f%n(:) = hdr%ngfft(:)
    n = f%n
    allocate(f%f(n(1),n(2),n(3)),stat=istat)
    if (istat /= 0) &
       call ferror('read_abinit','Error allocating grid',faterr,file)
    allocate(g(n(1),n(2),n(3)))
    read(luc,iostat=istat) g
    f%f = g
    deallocate(g)
    if (istat /= 0) &
       call ferror('read_abinit','Error reading grid',faterr,file)

    call fclose(luc)

  end subroutine read_abinit

  !> Read a grid in VASP format
  subroutine read_vasp(f,file,omega)
    use tools_io, only: fopen_read, getline_raw, faterr, ferror, fclose
    class(grid3), intent(inout) :: f
    character*(*), intent(in) :: file !< Input file
    real*8, intent(in) :: omega !< Cell volume

    integer :: luc
    integer :: istat, n(3), i, j, k
    character(len=:), allocatable :: line
    logical :: ok

    call f%end()
    luc = fopen_read(file)

    do while(.true.)
       ok = getline_raw(luc,line,.true.)
       if (len(trim(line)) == 0) exit
    end do

    read (luc,*,iostat=istat) n
    if (istat /= 0) &
       call ferror('read_vasp','Error reading nx, ny, nz',faterr,file)

    f%isinit = .true.
    f%mode = mode_default
    f%n(:) = n
    allocate(f%f(n(1),n(2),n(3)),stat=istat)
    if (istat /= 0) &
       call ferror('read_vasp','Error allocating grid',faterr,file)
    read(luc,*,iostat=istat) (((f%f(i,j,k),i=1,n(1)),j=1,n(2)),k=1,n(3))
    if (istat /= 0) &
       call ferror('read_vasp','Error reading grid',faterr,file)
    f%f(:,:,:) = f%f(:,:,:) / omega

    call fclose(luc)

  end subroutine read_vasp

  !> Read a grid in aimpac qub format
  subroutine read_qub(f,file)
    use tools_io, only: fopen_read, ferror, faterr, fclose
    class(grid3), intent(inout) :: f
    character*(*), intent(in) :: file !< Input file

    integer :: luc
    integer :: istat, n(3), i, j, k

    call f%end()
    luc = fopen_read(file)

    read (luc,*,iostat=istat) f%n
    if (istat /= 0) &
       call ferror('read_qub','Error reading n1, n2, n3',faterr,file)

    f%isinit = .true.
    f%mode = mode_default
    n = f%n(:)
    allocate(f%f(n(1),n(2),n(3)),stat=istat)
    if (istat /= 0) &
       call ferror('read_qub','Error allocating grid',faterr,file)
    read(luc,*,iostat=istat) (((f%f(i,j,k),i=1,n(1)),j=1,n(2)),k=1,n(3))
    if (istat /= 0) &
       call ferror('read_qub','Error reading grid',faterr,file)

    call fclose(luc)

  end subroutine read_qub

  !> Read a grid in xcrysden xsf format -- only first 3d grid in first 3d block
  subroutine read_xsf(f,file)
    use tools_io, only: fopen_read, getline_raw, lgetword, equal, ferror, faterr, &
       fclose
    use types, only: realloc
    class(grid3), intent(inout) :: f
    character*(*), intent(in) :: file !< Input file

    integer :: luc
    integer :: istat, n(3), lp, i, j, k
    character(len=:), allocatable :: line, word
    logical :: found, ok
    real*8, dimension(3) :: x0, x1, x2, x3
    real*8 :: pmat(3,3)
    real*8, allocatable :: ggloc(:,:,:)

    call f%end()

    ! open file for reading
    luc = fopen_read(file)

    ! position at the beginning of the first grid, ignore the rest
    found = .false.
    ok = .false.
    do while (getline_raw(luc,line))
       lp = 1
       word = lgetword(line,lp)
       if (equal(word,'primvec'))then
          ok = .true.
          read(luc,*,iostat=istat) pmat
          if (istat /= 0) &
             call ferror('read_xsf','Error PRIMVEC',faterr,file)
       else if (equal(word,'begin_block_datagrid_3d').or.equal(word,'begin_block_datagrid3d'))then
          found = .true.
          exit
       end if
    end do
    if (.not.found) call ferror('read_xsf','BEGIN_BLOCK_DATAGRID_3D not found',faterr,file)
    if (.not.ok) call ferror('read_xsf','PRIMVEC not found',faterr,file)

    found = .false.
    do while (getline_raw(luc,line))
       lp = 1
       word = lgetword(line,lp)
       if (equal(word(1:min(17,len(word))),'begin_datagrid_3d') .or. equal(word(1:min(11,len(word))),'datagrid_3d').or.&
           equal(word(1:min(16,len(word))),'begin_datagrid3d') .or. equal(word(1:min(10,len(word))),'datagrid3d')) then
          found = .true.
          exit
       end if
    end do
    if (.not.found) call ferror('read_xsf','BEGIN_DATAGRID_3D... not found',faterr,file)

    ! grid dimension
    read (luc,*,iostat=istat) n 
    if (istat /= 0) &
       call ferror('read_xsf','Error reading n1, n2, n3',faterr,file)
    f%n = n - 1

    ! origin and edge vectors
    read (luc,*,iostat=istat) x0, x1, x2, x3

    f%isinit = .true.
    f%mode = mode_default
    allocate(ggloc(n(1),n(2),n(3)),stat=istat)
    if (istat /= 0) &
       call ferror('read_xsf','Error allocating grid',faterr,file)
    if (istat /= 0) &
       call ferror('read_xsf','Error allocating grid',faterr,file)
    read(luc,*,iostat=istat) (((ggloc(i,j,k),i=1,n(1)),j=1,n(2)),k=1,n(3))
    if (istat /= 0) &
       call ferror('read_xsf','Error reading grid',faterr,file)

    allocate(f%f(f%n(1),f%n(2),f%n(3)),stat=istat)
    f%f = ggloc(1:n(1)-1,1:n(2)-1,1:n(3)-1)
    deallocate(ggloc)
    n = f%n

    call fclose(luc)

  end subroutine read_xsf

  !> Read the UNK generated by QE for Wannier and the wannier90
  !> checkpoint file. The part that reads the checkpoint file has been
  !> adapted from wannier90,
  !> Copyright (C) 2007-13 Jonathan Yates, Arash Mostofi,       
  !>                Giovanni Pizzi, Young-Su Lee,               
  !>                Nicola Marzari, Ivo Souza, David Vanderbilt 
  !> Distributed under GNU/GPL v2.
  subroutine read_unk(f,file,filedn,omega,nou,dochk)
    use tools_math, only: det, matinv
    use tools_io, only: fopen_read, getline_raw, lgetword, equal, ferror, faterr, &
       fclose, string, fopen_write, uout
    use types, only: realloc
    use param, only: bohrtoa
    class(grid3), intent(inout) :: f
    character*(*), intent(in) :: file !< Input file (spin up or total)
    character*(*), intent(in) :: filedn !< Input file (spin down)
    real*8, intent(in) :: omega
    logical, intent(in) :: nou
    logical, intent(in) :: dochk

    integer :: luc, luw
    integer :: nspin, ispin, ibnd, nbnd, jbnd, idum, nall(3)
    integer :: n(3), ik, ik1, ik2, ik3, naux(3), ikk
    real*8 :: fspin
    integer :: i, j, k, l
    complex*16, allocatable :: raux(:,:,:), raux2(:,:,:), rseq(:)
    integer :: nk1, nk2, nk3, nk
    character(len=:), allocatable :: fname, oname, sijfname
    logical :: haschk
    ! for the wannier checkpoint (wannier90, 2.0.1)
    character(len=33) :: header
    real*8 :: rlatt(3,3), rclatt(3,3), rlatti(3,3)
    character(len=20) :: chkpt1
    logical :: have_disentangled
    complex*16, allocatable :: u_matrix(:,:,:)
    complex*16 :: cdum

    call f%end()

    ! spin
    if (len_trim(filedn) < 1) then
       nspin = 1
       fspin = 2d0
    else
       nspin = 2
       fspin = 1d0
    end if

    ! checkpoint files
    haschk = .false.
    if (dochk) then
       sijfname = trim(file) // "-sij"
       inquire(file=sijfname,exist=haschk)
    end if

    ! read the grid size from the unk or the checkpoint
    if (haschk) then
       luc = fopen_read(sijfname,form="unformatted")
       read (luc) n
       call fclose(luc)
    else
       luc = fopen_read("UNK00001.1",form="unformatted")
       read(luc) n, idum, jbnd
       call fclose(luc)
    endif
    f%n = n

    ! allocate the density
    if (allocated(f%f)) deallocate(f%f)
    allocate(f%f(f%n(1),f%n(2),f%n(3)))
    f%f = 0d0

    ! read the chk file/s
    do ispin = 1, nspin
       if (ispin == 1) then
          luc = fopen_read(file,form="unformatted")
       else
          luc = fopen_read(filedn,form="unformatted")
       end if

       ! header and number of bands
       read(luc) header
       read(luc) nbnd 
       read(luc) jbnd 
       if (jbnd > 0) &
          call ferror("read_unk","number of excluded bands /= 0",faterr)
       read(luc) (idum,i=1,jbnd) 

       ! real and reciprocal lattice
       read(luc) ((rlatt(i,j),i=1,3),j=1,3)
       if (abs(det(rlatt) / bohrtoa**3 - omega) / omega > 1d-2) & 
          call ferror("read_unk","wannier and current structure's volumes differ by more than 1%",faterr)
       read(luc) ((rclatt(i,j),i=1,3),j=1,3)
    
       ! number of k-points
       read(luc) nk 
       read(luc) nk1, nk2, nk3
       if (nk == 0 .or. nk1 == 0 .or. nk2 == 0 .or. nk3 == 0 .or. nk /= (nk1*nk2*nk3)) &
          call ferror("read_unk","no monkhorst-pack grid or inconsistent k-point number",faterr)

       ! k-points
       allocate(f%wan%kpt(3,nk))
       read(luc) ((f%wan%kpt(i,j),i=1,3),j=1,nk) 
       do i = 1, nk
          ik1 = nint(f%wan%kpt(1,i) * nk1)
          ik2 = nint(f%wan%kpt(2,i) * nk2)
          ik3 = nint(f%wan%kpt(3,i) * nk3)
          if (abs(f%wan%kpt(1,i) * nk1 - ik1) > 1d-5 .or.abs(f%wan%kpt(2,i) * nk2 - ik2) > 1d-5 .or.&
             abs(f%wan%kpt(3,i) * nk3 - ik3) > 1d-5) then
             write (uout,*) f%wan%kpt(:,i)
             write (uout,*) f%wan%kpt(1,i)*nk1,f%wan%kpt(1,i)*nk2,f%wan%kpt(1,i)*nk3
             write (uout,*) ik1, ik2, ik3
             call ferror("read_unk","not a (uniform) monkhorst-pack grid or shifted grid",faterr)
          end if
       end do

       read(luc) idum ! number of nearest k-point neighbours
       read(luc) jbnd ! number of wannier functions
       if (jbnd /= nbnd) &
          call ferror("read_unk","number of wannier functions /= number of bands",faterr)
       
       ! checkpoint positon and disentanglement
       read(luc) chkpt1
       read(luc) have_disentangled
       if (have_disentangled) & 
          call ferror("read_unk","can not handle disentangled wannier functions",faterr)

       ! u and m matrices
       if (allocated(u_matrix)) deallocate(u_matrix)
       allocate(u_matrix(nbnd,nbnd,nk))
       read(luc) (((u_matrix(i,j,k),i=1,nbnd),j=1,nbnd),k=1,nk)
       read(luc) ((((cdum,i=1,nbnd),j=1,nbnd),k=1,idum),l=1,nk) ! m matrix

       ! wannier centers and spreads
       if (ispin == 1) then
          if (allocated(f%wan%center)) deallocate(f%wan%center)
          if (allocated(f%wan%spread)) deallocate(f%wan%spread)
          allocate(f%wan%center(3,nbnd,nspin),f%wan%spread(nbnd,nspin))
       else
          if (nbnd > size(f%wan%spread,1)) then
             call realloc(f%wan%center,3,nbnd,nspin)
             call realloc(f%wan%spread,nbnd,nspin)
          end if
       end if
       read(luc) ((f%wan%center(i,j,ispin),i=1,3),j=1,nbnd)
       read(luc) (f%wan%spread(i,ispin),i=1,nbnd)

       ! end of wannier checkpoint
       call fclose(luc)

       ! dimensions for the supercell
       f%wan%nwan = (/nk1,nk2,nk3/)
       f%wan%nks = nk1*nk2*nk3
       nall = f%n * f%wan%nwan

       ! convert centers to crystallographic and spread to bohr
       rlatti = matinv(rlatt)
       do i = 1, nbnd
          f%wan%center(:,i,ispin) = matmul(f%wan%center(:,i,ispin),rlatti)
          do j = 1, 3
             if (f%wan%center(j,i,ispin) > f%wan%nwan(j)) &
                f%wan%center(j,i,ispin) = f%wan%center(j,i,ispin) - f%wan%nwan(j)
             if (f%wan%center(j,i,ispin) < 0d0) &
                f%wan%center(j,i,ispin) = f%wan%center(j,i,ispin) + f%wan%nwan(j)
          end do
          f%wan%spread(i,ispin) = sqrt(f%wan%spread(i,ispin)) / bohrtoa
       end do

       if (.not.haschk) then
          ! allocate arrays and fill some info
          if (.not.nou) allocate(raux(n(1),n(2),n(3)))
          allocate(raux2(n(1),n(2),n(3)))

          do ik = 1, nk
             ! unk file name
             fname = "UNK" // string(ik,5,pad0=.true.) // "." // string(ispin)

             ! open file for reading
             luc = fopen_read(fname,form="unformatted")

             ! read the header
             read(luc) naux, ikk, idum

             ! read the bands and pass them to the WNK file open file for writing
             do ibnd = 1, nbnd
                oname = "WNK." // string(ik) // "." // string(ibnd) // "." // string(ispin)
                luw = fopen_write(oname,form="unformatted")

                if (.not.nou) then
                   ! apply the transformation
                   raux2 = 0d0
                   rewind(luc)
                   read(luc) naux, ikk, idum
                   do jbnd = 1, nbnd
                      read(luc) raux
                      raux2 = raux2 + u_matrix(jbnd,ibnd,ik) * raux
                   end do
                else
                   ! transfer this band to the new file
                   read(luc) raux2
                end if

                ! write and close
                write(luw) raux2
                call fclose(luw)

                f%f = f%f + conjg(raux2) * raux2
             end do

             ! close the unk file
             call fclose(luc)
          end do
          if (allocated(raux)) deallocate(raux)
          deallocate(raux2)
       end if
    end do

    if (haschk) then
       ! read the density from the checkpoint file
       luc = fopen_read(sijfname,form="unformatted")
       read (luc) n
       read (luc) f%f
       call fclose(luc)
    else
       f%f = f%f * fspin / omega / real(nk,8)
    end if

    f%wan%sijchk = dochk
    f%wan%haschk = haschk
    f%isinit = .true.
    f%mode = mode_default
    f%iswan = .true.
    f%wan%nbnd = nbnd
    f%wan%nspin = nspin
    if (allocated(f%wan%ngk)) deallocate(f%wan%ngk)
    if (allocated(f%wan%igk_k)) deallocate(f%wan%igk_k)
    if (allocated(f%wan%nls)) deallocate(f%wan%nls)
    if (allocated(f%wan%u)) deallocate(f%wan%u)

  end subroutine read_unk

  !> Read unkgen file created by pw2wannier.x.
  subroutine read_unkgen(f,fchk,fchkdn,funkgen,fevc,omega,dochk)
    use tools_math, only: det, matinv
    use tools_io, only: fopen_read, getline_raw, lgetword, equal, ferror, faterr, &
       fclose, string, fopen_write, uout
    use types, only: realloc
    use param, only: bohrtoa
    class(grid3), intent(inout) :: f
    character*(*), intent(in) :: fchk !< Input file (chk file from wannier90)
    character*(*), intent(in) :: fchkdn !< Input file (chk for the down spin)
    character*(*), intent(in) :: funkgen !< unkgen file (unkgen file from wannier90)
    character*(*), intent(in) :: fevc !< unkgen file (evc file from pw2wannier)
    real*8, intent(in) :: omega !< unit cell
    logical, intent(in) :: dochk !< use checkpoint file

    integer :: luc, luw
    integer :: nspin, ispin, ibnd, nbnd, jbnd, idum, nall(3)
    integer :: n(3), ik, ik1, ik2, ik3, naux(3), ikk
    real*8 :: fspin
    integer :: i, j, k, l
    complex*16, allocatable :: raux(:,:,:), raux2(:,:,:), rseq(:), evc(:)
    integer :: nk1, nk2, nk3, nk
    character(len=:), allocatable :: fname, oname, sijfname
    ! for the wannier checkpoint (wannier90, 2.0.1)
    character(len=33) :: header
    real*8 :: rlatt(3,3), rclatt(3,3), rlatti(3,3)
    character(len=20) :: chkpt1
    logical :: have_disentangled, ok1, ok2, haschk
    complex*16 :: cdum

    call f%end()

    ! spin
    if (len_trim(fchkdn) < 1) then
       nspin = 1
       fspin = 2d0
    else
       nspin = 2
       fspin = 1d0
    end if
    
    ! checkpoint files
    haschk = .false.
    if (dochk) then
       sijfname = trim(fchk) // "-sij"
       inquire(file=sijfname,exist=haschk)
    end if
    if (.not.haschk) then
       inquire(file=funkgen,exist=ok1)
       inquire(file=fevc,exist=ok2)
       if (.not.ok1.or..not.ok2) &
          call ferror("read_unkgen","unkgen/evc and chk-sij files not found",faterr)
    end if

    ! read the grid size from the unkgen or the checkpoint
    if (haschk) then
       luc = fopen_read(sijfname,form="unformatted")
       read (luc) n
       call fclose(luc)
    else
       luc = fopen_read(funkgen,form="unformatted")
       read (luc) ikk, idum, ibnd, ispin
       read (luc) n
       call fclose(luc)
    end if
    f%n = n

    ! read the chk file/s
    do ispin = 1, nspin
       if (ispin == 1) then
          luc = fopen_read(fchk,form="unformatted")
       else
          luc = fopen_read(fchkdn,form="unformatted")
       end if

       ! header and number of bands
       read(luc) header
       read(luc) nbnd 
       read(luc) jbnd 
       if (jbnd > 0) &
          call ferror("read_unkgen","number of excluded bands /= 0",faterr)
       read(luc) (idum,i=1,jbnd) 

       ! real and reciprocal lattice
       read(luc) ((rlatt(i,j),i=1,3),j=1,3)
       if (abs(det(rlatt) / bohrtoa**3 - omega) / omega > 1d-2) & 
          call ferror("read_unkgen","wannier and current structure's volumes differ by more than 1%",faterr)
       read(luc) ((rclatt(i,j),i=1,3),j=1,3)
    
       ! number of k-points
       read(luc) nk 
       read(luc) nk1, nk2, nk3
       if (nk == 0 .or. nk1 == 0 .or. nk2 == 0 .or. nk3 == 0 .or. nk /= (nk1*nk2*nk3)) &
          call ferror("read_unkgen","no monkhorst-pack grid or inconsistent k-point number",faterr)

       ! k-points
       if (allocated(f%wan%kpt)) deallocate(f%wan%kpt)
       allocate(f%wan%kpt(3,nk))
       read(luc) ((f%wan%kpt(i,j),i=1,3),j=1,nk) 
       do i = 1, nk
          ik1 = nint(f%wan%kpt(1,i) * nk1)
          ik2 = nint(f%wan%kpt(2,i) * nk2)
          ik3 = nint(f%wan%kpt(3,i) * nk3)
          if (abs(f%wan%kpt(1,i) * nk1 - ik1) > 1d-5 .or.abs(f%wan%kpt(2,i) * nk2 - ik2) > 1d-5 .or.&
             abs(f%wan%kpt(3,i) * nk3 - ik3) > 1d-5) then
             write (uout,*) f%wan%kpt(:,i)
             write (uout,*) f%wan%kpt(1,i)*nk1,f%wan%kpt(1,i)*nk2,f%wan%kpt(1,i)*nk3
             write (uout,*) ik1, ik2, ik3
             call ferror("read_unkgen","not a (uniform) monkhorst-pack grid or shifted grid",faterr)
          end if
       end do

       read(luc) idum ! number of nearest k-point neighbours
       read(luc) jbnd ! number of wannier functions
       if (jbnd /= nbnd) &
          call ferror("read_unkgen","number of wannier functions /= number of bands",faterr)
       
       ! checkpoint positon and disentanglement
       read(luc) chkpt1
       read(luc) have_disentangled
       if (have_disentangled) & 
          call ferror("read_unkgen","can not handle disentangled wannier functions",faterr)

       ! u and m matrices
       if (allocated(f%wan%u)) deallocate(f%wan%u)
       allocate(f%wan%u(nbnd,nbnd,nk))
       read(luc) (((f%wan%u(i,j,k),i=1,nbnd),j=1,nbnd),k=1,nk)
       read(luc) ((((cdum,i=1,nbnd),j=1,nbnd),k=1,idum),l=1,nk) ! m matrix

       ! wannier centers and spreads
       if (ispin == 1) then
          if (allocated(f%wan%center)) deallocate(f%wan%center)
          if (allocated(f%wan%spread)) deallocate(f%wan%spread)
          allocate(f%wan%center(3,nbnd,nspin),f%wan%spread(nbnd,nspin))
       else
          if (nbnd > size(f%wan%spread,1)) then
             call realloc(f%wan%center,3,nbnd,nspin)
             call realloc(f%wan%spread,nbnd,nspin)
          end if
       end if
       read(luc) ((f%wan%center(i,j,ispin),i=1,3),j=1,nbnd)
       read(luc) (f%wan%spread(i,ispin),i=1,nbnd)

       ! end of wannier checkpoint
       call fclose(luc)

       ! dimensions for the supercell
       f%wan%nwan = (/nk1,nk2,nk3/)
       nall = f%n * f%wan%nwan

       ! convert centers to crystallographic and spread to bohr
       rlatti = matinv(rlatt)
       do i = 1, nbnd
          f%wan%center(:,i,ispin) = matmul(f%wan%center(:,i,ispin),rlatti)
          do j = 1, 3
             if (f%wan%center(j,i,ispin) > f%wan%nwan(j)) &
                f%wan%center(j,i,ispin) = f%wan%center(j,i,ispin) - f%wan%nwan(j)
             if (f%wan%center(j,i,ispin) < 0d0) &
                f%wan%center(j,i,ispin) = f%wan%center(j,i,ispin) + f%wan%nwan(j)
          end do
          f%wan%spread(i,ispin) = sqrt(f%wan%spread(i,ispin)) / bohrtoa
       end do
    end do
    f%wan%nks = nk
    f%wan%nbnd = nbnd
    f%wan%nspin = nspin

    ! save the evc file name
    f%wan%fevc = fevc

    if (.not.haschk) then
       ! read the unkgen info
       luc = fopen_read(funkgen,form="unformatted")
       read (luc) ikk, idum, ibnd, ispin
       if (ikk /= 1 .or. idum /= nk .or. ibnd /= nbnd .or. ispin /= nspin) &
          call ferror("read_unkgen","inconsistent unkgen",faterr)

       read (luc) n
       if (any(n /= f%n)) &
          call ferror("read_unkgen","inconsistent unkgen",faterr)

       ! allocate and read integer indices
       read(luc) naux
       if (allocated(f%wan%ngk)) deallocate(f%wan%ngk)
       allocate(f%wan%ngk(nk))
       if (allocated(f%wan%igk_k)) deallocate(f%wan%igk_k)
       allocate(f%wan%igk_k(naux(2),nk))
       if (allocated(f%wan%nls)) deallocate(f%wan%nls)
       allocate(f%wan%nls(naux(3)))
       read (luc) f%wan%ngk
       read (luc) f%wan%igk_k
       read (luc) f%wan%nls
       call fclose(luc)
    end if

    ! allocate the density
    if (allocated(f%f)) deallocate(f%f)
    allocate(f%f(f%n(1),f%n(2),f%n(3)))
    f%f = 0d0

    if (.not.haschk) then
       ! open the evc file
       luc = fopen_read(fevc,form="unformatted")
       allocate(evc(maxval(f%wan%ngk(1:nk))))

       ! calculate the electron density
       allocate(raux(n(1),n(2),n(3)),rseq(n(1)*n(2)*n(3)))
       do ispin = 1, nspin
          do ik = 1, nk
             do ibnd = 1, nbnd
                read (luc) evc(1:f%wan%ngk(ik))
                rseq = 0d0
                rseq(f%wan%nls(f%wan%igk_k(1:f%wan%ngk(ik),ik))) = evc(1:f%wan%ngk(ik))
                raux = reshape(rseq,shape(raux))
                call cfftnd(3,n,+1,raux)
                f%f = f%f + conjg(raux) * raux
             end do
          end do
       end do
       deallocate(raux,rseq,evc)
       call fclose(luc)
       f%f = f%f * fspin / omega / real(nk,8)
    else
       luc = fopen_read(sijfname,form="unformatted")
       read (luc) n
       read (luc) f%f
       call fclose(luc)
    end if

    f%wan%sijchk = dochk
    f%wan%haschk = haschk
    f%isinit = .true.
    f%mode = mode_default
    f%iswan = .true.
    f%wan%nbnd = nbnd
    f%wan%nspin = nspin

  end subroutine read_unkgen

  !> Read a grid in elk format -- only first 3d grid in first 3d block
  subroutine read_elk(f,file)
    use tools_io, only: fopen_read, ferror, faterr, fclose
    class(grid3), intent(inout) :: f
    character*(*), intent(in) :: file !< Input file

    integer :: luc, ios
    integer :: n(3), i, j, k
    real*8 :: dum(3)

    call f%end()

    ! open file for reading
    luc = fopen_read(file)

    ! grid dimension
    read (luc,*,iostat=ios) n
    if (ios /= 0) &
       call ferror('read_elk','Error reading n1, n2, n3',faterr,file)

    f%n = n
    allocate(f%f(n(1),n(2),n(3)),stat=ios)
    if (ios /= 0) &
       call ferror('read_elk','Error allocating grid',faterr,file)
    do k = 1, n(3)
       do j = 1, n(2)
          do i = 1, n(1)
             read (luc,*) dum, f%f(i,j,k)
          end do
       end do
    end do

    call fclose(luc)

    n = f%n
    f%isinit = .true.
    f%mode = mode_default

  end subroutine read_elk

  !> Interpolate the function value, first and second derivative at
  !> point x0 (crystallographic coords.) using the grid g.  This
  !> routine is thread-safe.
  subroutine interp(f,xi,y,yp,ypp) 
    class(grid3), intent(inout) :: f !< Grid to interpolate
    real*8, intent(in) :: xi(3) !< Target point (cryst. coords.)
    real*8, intent(out) :: y !< Interpolated value
    real*8, intent(out) :: yp(3) !< First derivative
    real*8, intent(out) :: ypp(3,3) !< Second derivative

    real*8 :: x0(3)

    x0 = xi - floor(xi)
    y = 0d0
    yp = 0d0
    ypp = 0d0

    if (f%mode == mode_nearest) then
       call grinterp_nearest(f,x0,y)
       yp = 0d0
       ypp = 0d0
    else if (f%mode == mode_trilinear) then
       call grinterp_trilinear(f,x0,y,yp)
       ypp = 0d0
    else if (f%mode == mode_trispline) then
       call grinterp_trispline(f,x0,y,yp,ypp)
    else if (f%mode == mode_tricubic) then
       call grinterp_tricubic(f,x0,y,yp,ypp)
    end if

  end subroutine interp

  !> Interpolate by giving the value of the grid at the nearest point
  !> (or the would-be nearest, it is nearest only in orthogonal
  !> grids). f is the grid field to interpolate. xi is the point in
  !> crystallographic coordinates. y is the interpolated value at xi.
  !> This routine is thread-safe.
  subroutine grinterp_nearest(f,x0,y)
    class(grid3), intent(in) :: f !< Input grid
    real*8, intent(in) :: x0(3) !< Target point (cryst. coords.)
    real*8, intent(out) :: y !< Interpolated value

    real*8 :: x(3)
    integer :: idx(3)

    x = modulo(x0,1d0)
    idx = grid_near(f,x)
    y = f%f(idx(1),idx(2),idx(3))

  end subroutine grinterp_nearest

  !> Interpolate using a trilinear interpolant. f is the grid field to
  !> interpolate. xi is the point in crystallographic coordinates. y
  !> and yp are the value and gradient at xi.  This routine is
  !> thread-safe.
  subroutine grinterp_trilinear(f,x0,y,yp)
    class(grid3), intent(in) :: f !< Input grid
    real*8, intent(in) :: x0(3) !< Target point (cryst. coords.)
    real*8, intent(out) :: y !< Interpolated value
    real*8, intent(out) :: yp(3) !< First derivative

    integer :: idx(3), iidx(3), i, j, k
    real*8 :: ff(0:2,0:2,0:2), r(3), s(3), x(3)

    ! compute value at the cube vertices
    x = modulo(x0,1d0)
    idx = grid_floor(f,x)
    do i = 0, 1
       do j = 0, 1
          do k = 0, 1
             iidx = modulo(idx+(/i,j,k/)-1,f%n)+1
             ff(i,j,k) = f%f(iidx(1),iidx(2),iidx(3))
          end do
       end do
    end do

    ! x and 1-x
    r = f%n * x - idx + 1
    s = 1d0 - r

    ! trilinear interpolation
    do i = 0, 1
       do j = 0, 1
          ff(i,j,2) = ff(i,j,0) * s(3) + ff(i,j,1) * r(3)
          ff(2,i,j) = ff(0,i,j) * s(1) + ff(1,i,j) * r(1)
       end do
       ff(i,2,2) = ff(i,0,2) * s(2) + ff(i,1,2) * r(2)
       ff(2,i,2) = ff(0,i,2) * s(1) + ff(1,i,2) * r(1)
       ff(2,2,i) = ff(2,0,i) * s(2) + ff(2,1,i) * r(2)
    end do
    ff(2,2,2) = ff(0,2,2) * s(1) + ff(1,2,2) * r(1)
    y = ff(2,2,2)
    yp(1) = ff(1,2,2) - ff(0,2,2)
    yp(2) = ff(2,1,2) - ff(2,0,2)
    yp(3) = ff(2,2,1) - ff(2,2,0)

    ! from cell to crystallographic coordinates
    yp = yp * f%n

  end subroutine grinterp_trilinear

  !> Using the global cubic spline information (c2 array) with
  !> periodic boundary conditions, calculate the density and c2
  !> coeffs.  of the 6 points on the cube faces, and then average the
  !> spline prediction. This is a modified version of the abinit
  !> density interpolation subroutine. f is the grid field to
  !> interpolate. xi is the point in crystallographic coordinates. y,
  !> yp, and ypp are the value, gradient, and Hessian at xi.  This
  !> routine is thread-safe.
  subroutine grinterp_trispline(f,x0,y,yp,ypp)
    class(grid3), intent(inout), target :: f !< Input grid
    real*8, intent(in) :: x0(3) !< Target point
    real*8, intent(out) :: y !< Interpolated value
    real*8, intent(out) :: yp(3) !< First derivative
    real*8, intent(out) :: ypp(3,3) !< Second derivative

    integer :: i, ii, jj, kk, ll, nn, indx(3), oii, onn, omm
    integer :: inii(4,3)
    real*8 :: bbb, ddu(2), hrh(2), hh(4,2), grd(4), lder(4)
    real*8 :: cof(2,3), ddstar(6), rhstar(6), sqder(6,4), sqvlr(6,4)
    real*8 :: pomsq(2,3), pom2sq(2,3)
    real*8,pointer :: ptddx(:,:,:),ptddy(:,:,:),ptddz(:,:,:),ptrho(:,:,:)
    real*8 :: xx(3), dix(3)

    real*8, parameter :: eps = 1d-12

    !$omp critical (checkalloc)
    if (.not.allocated(f%c2)) then
       call init_trispline(f)
    end if
    !$omp end critical (checkalloc)

    nullify(ptddx,ptddy,ptddz,ptrho)

    xx = modulo(x0,1d0)

    do i=1,3
       dix(i)=1d0/f%n(i)
    end do

    ! determine the index in the grid
    do ii=1,3
       indx(ii)=int(xx(ii)*f%n(ii))
       bbb=(xx(ii)-indx(ii)*dix(ii))*f%n(ii)
       if (indx(ii)==f%n(ii)) then
          indx(ii)=1
          xx(ii)=0.d0
       else
          indx(ii)=indx(ii)+1
       end if
       !  Explicit handling to avoid numeric problems
       if (bbb > 1.d0+eps ) then
          cof(1,ii)=0.d0
          cof(2,ii)=1.d0
       elseif (bbb < -eps ) then
          cof(1,ii)=1.d0
          cof(2,ii)=0.d0
       else
          cof(1,ii)=1.d0-bbb
          cof(2,ii)=bbb
       end if
    end do

    ! 3d interpolation of the valence density

    ! determination of the values of density and of its second derivative
    ! at the "star" = constructed at vv with primitive directions
    ! To interpolation the values at the faces of the grid cell are needed

    rhstar(:)=0.d0
    sqder(:,:)=0.d0
    sqvlr(:,:)=0.d0
    ddstar(:)=0.d0
    pomsq(:,:)=0.d0
    pom2sq(:,:)=0.d0

    oii=1; onn=1; omm=1
    if (indx(1)==f%n(1)) oii=1-f%n(1)
    if (indx(2)==f%n(2)) onn=1-f%n(2)
    if (indx(3)==f%n(3)) omm=1-f%n(3)

    ! the values in the corners of the grid cell

    ptddx=>f%c2(indx(1):indx(1)+oii:oii,indx(2):indx(2)+onn:onn,indx(3):indx(3)+omm:omm,1)
    ptddy=>f%c2(indx(1):indx(1)+oii:oii,indx(2):indx(2)+onn:onn,indx(3):indx(3)+omm:omm,2)
    ptddz=>f%c2(indx(1):indx(1)+oii:oii,indx(2):indx(2)+onn:onn,indx(3):indx(3)+omm:omm,3)
    ptrho=>f%f(indx(1):indx(1)+oii:oii,indx(2):indx(2)+onn:onn,indx(3):indx(3)+omm:omm)

    ! the coefficients for spline interpolation of density and its derivation
    do ii=1,3
       do jj=1,2
          pomsq(jj,ii)=(cof(jj,ii)*cof(jj,ii)*cof(jj,ii)-cof(jj,ii))/6.d0*dix(ii)*dix(ii)
          pom2sq(jj,ii)=(3.d0*cof(jj,ii)*cof(jj,ii)-1.d0)/6.d0*dix(ii)
          if (jj==1) pom2sq(jj,ii)=-pom2sq(jj,ii)
       end do
    end do


    do ii=1,2
       do jj=1,2
          do kk=1,2
             ddstar(ii)=ddstar(ii)+cof(jj,2)*cof(kk,3)*ptddx(ii,jj,kk)
             ddstar(ii+2)=ddstar(ii+2)+cof(jj,3)*cof(kk,1)*ptddy(kk,ii,jj)
             ddstar(ii+4)=ddstar(ii+4)+cof(jj,1)*cof(kk,2)*ptddz(jj,kk,ii)
             sqder(ii,jj)=sqder(ii,jj)+cof(kk,2)*ptddz(ii,kk,jj)
             sqder(ii,jj+2)=sqder(ii,jj+2)+cof(kk,3)*ptddy(ii,jj,kk)
             sqder(ii+2,jj)=sqder(ii+2,jj)+cof(kk,3)*ptddx(jj,ii,kk)
             sqder(ii+2,jj+2)=sqder(ii+2,jj+2)+cof(kk,1)*ptddz(kk,ii,jj)
             sqder(ii+4,jj)=sqder(ii+4,jj)+cof(kk,1)*ptddy(kk,jj,ii)
             sqder(ii+4,jj+2)=sqder(ii+4,jj+2)+cof(kk,2)*ptddx(jj,kk,ii)
             sqvlr(ii,jj)=sqvlr(ii,jj)+cof(kk,2)*ptrho(ii,kk,jj)+pomsq(kk,2)*ptddy(ii,kk,jj)
             sqvlr(ii,jj+2)=sqvlr(ii,jj+2)+cof(kk,3)*ptrho(ii,jj,kk)+pomsq(kk,3)*ptddz(ii,jj,kk)
             sqvlr(ii+2,jj+2)=sqvlr(ii+2,jj+2)+cof(kk,1)*ptrho(kk,ii,jj)+pomsq(kk,1)*ptddx(kk,ii,jj)
          end do
       end do
    end do

    do ii=1,2
       do jj=1,2
          sqvlr(ii+2,jj)=sqvlr(jj,ii+2)
          sqvlr(ii+4,jj)=sqvlr(jj+2,ii+2)
          sqvlr(ii+4,jj+2)=sqvlr(jj,ii)
       end do
    end do

    do ii=1,2
       do jj=1,2
          rhstar(ii)=rhstar(ii)+cof(jj,3)*sqvlr(ii,jj)+pomsq(jj,3)*sqder(ii,jj)+&
             &    cof(jj,2)*sqvlr(ii,jj+2)+pomsq(jj,2)*sqder(ii,jj+2)
          rhstar(ii+2)=rhstar(ii+2)+cof(jj,1)*sqvlr(ii+2,jj)+pomsq(jj,1)*sqder(ii+2,jj)+&
             &    cof(jj,3)*sqvlr(ii+2,jj+2)+pomsq(jj,3)*sqder(ii+2,jj+2)
          rhstar(ii+4)=rhstar(ii+4)+cof(jj,2)*sqvlr(ii+4,jj)+pomsq(jj,2)*sqder(ii+4,jj)+&
             &    cof(jj,1)*sqvlr(ii+4,jj+2)+pomsq(jj,1)*sqder(ii+4,jj+2)
       end do
    end do
    rhstar(:)=rhstar(:)/2.d0

    y = 0d0
    yp = 0d0
    ypp = 0d0
    kk=1; nn=1
    do ii=1,5,2
       do jj=1,2
          nn=-nn
          y=y+cof(jj,kk)*rhstar(ii+jj-1)+pomsq(jj,kk)*ddstar(ii+jj-1)
          yp(kk)=yp(kk)+pom2sq(jj,kk)*ddstar(ii+jj-1)
          ypp(kk,kk)=ypp(kk,kk)+cof(jj,kk)*ddstar(ii+jj-1)
          yp(kk)=yp(kk)+nn*rhstar(ii+jj-1)/dix(kk)
       end do
       kk=kk+1
    end do
    y=y/3.d0

    ! Off-diagonal elements of the hessian

    ! for the speed reasons the polynomial interpolation
    ! for second derivation fields is used in this case
    ! but the last step is always done by spline interpolation.


    do ii=1,3
       do jj=-1,2
          inii(jj+2,ii)=indx(ii)+jj
          if (inii(jj+2,ii) < 1) inii(jj+2,ii)=inii(jj+2,ii)+f%n(ii)
          if (inii(jj+2,ii) > f%n(ii)) inii(jj+2,ii)=inii(jj+2,ii)-f%n(ii)
       end do
    end do

    ! Not very nice

    do ii=1,3
       select case (ii)
       case (1)
          do jj=1,4
             ddu(1)=cof(1,2)*f%c2(inii(jj,1),inii(2,2),inii(2,3),3)+cof(2,2)*f%c2(inii(jj,1),inii(3,2),inii(2,3),3)
             ddu(2)=cof(1,2)*f%c2(inii(jj,1),inii(2,2),inii(3,3),3)+cof(2,2)*f%c2(inii(jj,1),inii(3,2),inii(3,3),3)
             hrh(1)=cof(1,2)*f%f(inii(jj,1),inii(2,2),inii(2,3))+cof(2,2)*f%f(inii(jj,1),inii(3,2),inii(2,3))+&
                &      pomsq(1,2)*f%c2(inii(jj,1),inii(2,2),inii(2,3),2)+pomsq(2,2)*f%c2(inii(jj,1),inii(3,2),inii(2,3),2)
             hrh(2)=cof(1,2)*f%f(inii(jj,1),inii(2,2),inii(3,3))+cof(2,2)*f%f(inii(jj,1),inii(3,2),inii(3,3))+&
                &      pomsq(1,2)*f%c2(inii(jj,1),inii(2,2),inii(3,3),2)+pomsq(2,2)*f%c2(inii(jj,1),inii(3,2),inii(3,3),2)
             hh(jj,2)=(hrh(2)-hrh(1))/dix(3)+pom2sq(1,3)*ddu(1)+pom2sq(2,3)*ddu(2)

             ddu(1)=cof(1,3)*f%c2(inii(jj,1),inii(2,2),inii(2,3),2)+cof(2,3)*f%c2(inii(jj,1),inii(2,2),inii(3,3),2)
             ddu(2)=cof(1,3)*f%c2(inii(jj,1),inii(3,2),inii(2,3),2)+cof(2,3)*f%c2(inii(jj,1),inii(3,2),inii(3,3),2)
             hrh(1)=cof(1,3)*f%f(inii(jj,1),inii(2,2),inii(2,3))+cof(2,3)*f%f(inii(jj,1),inii(2,2),inii(3,3))+&
                &      pomsq(1,3)*f%c2(inii(jj,1),inii(2,2),inii(2,3),3)+pomsq(2,3)*f%c2(inii(jj,1),inii(2,2),inii(3,3),3)
             hrh(2)=cof(1,3)*f%f(inii(jj,1),inii(3,2),inii(2,3))+cof(2,3)*f%f(inii(jj,1),inii(3,2),inii(3,3))+&
                &      pomsq(1,3)*f%c2(inii(jj,1),inii(3,2),inii(2,3),3)+pomsq(2,3)*f%c2(inii(jj,1),inii(3,2),inii(3,3),3)
             hh(jj,1)=(hrh(2)-hrh(1))/dix(2)+pom2sq(1,2)*ddu(1)+pom2sq(2,2)*ddu(2)
          end do
       case (2)
          do jj=1,4
             ddu(1)=cof(1,3)*f%c2(inii(2,1),inii(jj,2),inii(2,3),1)+cof(2,3)*f%c2(inii(2,1),inii(jj,2),inii(3,3),1)
             ddu(2)=cof(1,3)*f%c2(inii(3,1),inii(jj,2),inii(2,3),1)+cof(2,3)*f%c2(inii(3,1),inii(jj,2),inii(3,3),1)
             hrh(1)=cof(1,3)*f%f(inii(2,1),inii(jj,2),inii(2,3))+cof(2,3)*f%f(inii(2,1),inii(jj,2),inii(3,3))+&
                &      pomsq(1,3)*f%c2(inii(2,1),inii(jj,2),inii(2,3),3)+pomsq(2,3)*f%c2(inii(2,1),inii(jj,2),inii(3,3),3)
             hrh(2)=cof(1,3)*f%f(inii(3,1),inii(jj,2),inii(2,3))+cof(2,3)*f%f(inii(3,1),inii(jj,2),inii(3,3))+&
                &      pomsq(1,3)*f%c2(inii(3,1),inii(jj,2),inii(2,3),3)+pomsq(2,3)*f%c2(inii(3,1),inii(jj,2),inii(3,3),3)
             hh(jj,2)=(hrh(2)-hrh(1))/dix(1)+pom2sq(1,1)*ddu(1)+pom2sq(2,1)*ddu(2)

             ddu(1)=cof(1,1)*f%c2(inii(2,1),inii(jj,2),inii(2,3),3)+cof(2,1)*f%c2(inii(3,1),inii(jj,2),inii(2,3),3)
             ddu(2)=cof(1,1)*f%c2(inii(2,1),inii(jj,2),inii(3,3),3)+cof(2,1)*f%c2(inii(3,1),inii(jj,2),inii(3,3),3)
             hrh(1)=cof(1,1)*f%f(inii(2,1),inii(jj,2),inii(2,3))+cof(2,1)*f%f(inii(3,1),inii(jj,2),inii(2,3))+&
                &      pomsq(1,1)*f%c2(inii(2,1),inii(jj,2),inii(2,3),1)+pomsq(2,1)*f%c2(inii(3,1),inii(jj,2),inii(2,3),1)
             hrh(2)=cof(1,1)*f%f(inii(2,1),inii(jj,2),inii(3,3))+cof(2,1)*f%f(inii(3,1),inii(jj,2),inii(3,3))+&
                &      pomsq(1,1)*f%c2(inii(2,1),inii(jj,2),inii(3,3),1)+pomsq(2,1)*f%c2(inii(3,1),inii(jj,2),inii(3,3),1)
             hh(jj,1)=(hrh(2)-hrh(1))/dix(3)+pom2sq(1,3)*ddu(1)+pom2sq(2,3)*ddu(2)
          end do
       case (3)
          do jj=1,4
             ddu(1)=cof(1,1)*f%c2(inii(2,1),inii(2,2),inii(jj,3),2)+cof(2,1)*f%c2(inii(3,1),inii(2,2),inii(jj,3),2)
             ddu(2)=cof(1,1)*f%c2(inii(2,1),inii(3,2),inii(jj,3),2)+cof(2,1)*f%c2(inii(3,1),inii(3,2),inii(jj,3),2)
             hrh(1)=cof(1,1)*f%f(inii(2,1),inii(2,2),inii(jj,3))+cof(2,1)*f%f(inii(3,1),inii(2,2),inii(jj,3))+&
                &      pomsq(1,1)*f%c2(inii(2,1),inii(2,2),inii(jj,3),1)+pomsq(2,1)*f%c2(inii(3,1),inii(2,2),inii(jj,3),1)
             hrh(2)=cof(1,1)*f%f(inii(2,1),inii(3,2),inii(jj,3))+cof(2,1)*f%f(inii(3,1),inii(3,2),inii(jj,3))+&
                &      pomsq(1,1)*f%c2(inii(2,1),inii(3,2),inii(jj,3),1)+pomsq(2,1)*f%c2(inii(3,1),inii(3,2),inii(jj,3),1)
             hh(jj,2)=(hrh(2)-hrh(1))/dix(2)+pom2sq(1,2)*ddu(1)+pom2sq(2,2)*ddu(2)

             ddu(1)=cof(1,2)*f%c2(inii(2,1),inii(2,2),inii(jj,3),1)+cof(2,2)*f%c2(inii(2,1),inii(3,2),inii(jj,3),1)
             ddu(2)=cof(1,2)*f%c2(inii(3,1),inii(2,2),inii(jj,3),1)+cof(2,2)*f%c2(inii(3,1),inii(3,2),inii(jj,3),1)
             hrh(1)=cof(1,2)*f%f(inii(2,1),inii(2,2),inii(jj,3))+cof(2,2)*f%f(inii(2,1),inii(3,2),inii(jj,3))+&
                &      pomsq(1,2)*f%c2(inii(2,1),inii(2,2),inii(jj,3),2)+pomsq(2,2)*f%c2(inii(2,1),inii(3,2),inii(jj,3),2)
             hrh(2)=cof(1,2)*f%f(inii(3,1),inii(2,2),inii(jj,3))+cof(2,2)*f%f(inii(3,1),inii(3,2),inii(jj,3))+&
                &      pomsq(1,2)*f%c2(inii(3,1),inii(2,2),inii(jj,3),2)+pomsq(2,2)*f%c2(inii(3,1),inii(3,2),inii(jj,3),2)
             hh(jj,1)=(hrh(2)-hrh(1))/dix(1)+pom2sq(1,1)*ddu(1)+pom2sq(2,1)*ddu(2)
          end do
       end select
       do jj=-2,1
          grd(jj+3)=(indx(ii)+jj)*dix(ii)
       end do

       !  write(6,'("hh: ",/,4F16.8,/,4F16.8)') ((hh(kk,jj),kk=1,4),jj=1,2)
       !  write(6,'("grad: ",3F16.8)') (grad(kk),kk=1,3)
       !  write(6,'("dix: ",3F16.8)') (dix(kk),kk=1,3)
       !  write(6,'("grd: ",4F16.8)') (grd(kk),kk=1,4)
       !  write(6,'("inii: ",4I4)') (inii(kk,ii),kk=1,4)

       do jj=1,2

          !   polynomial interpolation

          do kk=1,3
             do ll=4,kk+1,-1
                hh(ll,jj)=(hh(ll,jj)-hh(ll-1,jj))/(grd(ll)-grd(ll-1))
             end do
          end do
          lder(4)=hh(4,jj)
          do kk=3,1,-1
             lder(kk)=hh(kk,jj)+(xx(ii)-grd(kk))*lder(kk+1)
          end do
          do kk=1,2
             do ll=3,kk+1,-1
                lder(ll)=lder(ll)+(xx(ii)-grd(ll-kk))*lder(ll+1)
             end do
          end do
          nn=ii+jj
          if (nn > 3) nn=nn-3
          ypp(ii,nn)=ypp(ii,nn)+lder(2)
          ypp(nn,ii)=ypp(nn,ii)+lder(2)
       end do
    end do

    ! averaging of the mixed derivations obtained in different order
    do ii=1,3
       do jj=1,3
          if (ii /= jj) ypp(ii,jj)=ypp(ii,jj)/2.d0
       end do
    end do

    nullify(ptddx,ptddy,ptddz,ptrho)

  end subroutine grinterp_trispline

  !> Tricubic interpolation based on: Lekien and Marsden,
  !> Int. J. Numer. Meth. Engng, 63 (2005) 455-471.  This
  !> interpolation is C^1 and local (the interpolant uses information
  !> of the grid points close to the point).  This subroutine has been
  !> adapted from the likely code, by David Kirkby, University of
  !> California, Irvine (https://github.com/deepzot/likely). f is the
  !> grid field to interpolate. xi is the point in crystallographic
  !> coordinates. y, yp, and ypp are the value, gradient, and Hessian
  !> at xi.  This routine is thread-safe.
  subroutine grinterp_tricubic(f,xi,y,yp,ypp)
    class(grid3), intent(inout), target :: f !< Input grid
    real*8, intent(in) :: xi(3) !< Target point
    real*8, intent(out) :: y !< Interpolated value
    real*8, intent(out) :: yp(3) !< First derivative
    real*8, intent(out) :: ypp(3,3) !< Second derivative

    integer :: idx(3), iidx(3), i, j, k, l
    real*8 :: g(-1:2,-1:2,-1:2), x(3)
    real*8 :: a(64), b(64)
    real*8, dimension(0:3) :: aa, bb, bbx, bbxx, aax, aay, aaxy
    real*8, dimension(0:3) :: aaxx, aayy

    ! initialize
    y = 0d0
    yp = 0d0
    ypp = 0d0

    ! fetch values at the cube vertices
    x = modulo(xi,1d0)
    idx = grid_floor(f,x)
    do i = -1, 2
       do j = -1, 2
          do k = -1, 2
             iidx = modulo(idx+(/i,j,k/)-1,f%n)+1
             g(i,j,k) = f%f(iidx(1),iidx(2),iidx(3))
          end do
       end do
    end do

    ! Fill the values in the b-vector (see Lekien and Marsden)
    ! f
    b(1) = g(0,0,0)
    b(2) = g(1,0,0)
    b(3) = g(0,1,0)
    b(4) = g(1,1,0)
    b(5) = g(0,0,1)
    b(6) = g(1,0,1)
    b(7) = g(0,1,1)
    b(8) = g(1,1,1)

    ! fx
    b(9)  = 0.5d0*(g(1,0,0)-g(-1,0,0))
    b(10) = 0.5d0*(g(2,0,0)-g(0,0,0))
    b(11) = 0.5d0*(g(1,1,0)-g(-1,1,0))
    b(12) = 0.5d0*(g(2,1,0)-g(0,1,0))
    b(13) = 0.5d0*(g(1,0,1)-g(-1,0,1))
    b(14) = 0.5d0*(g(2,0,1)-g(0,0,1))
    b(15) = 0.5d0*(g(1,1,1)-g(-1,1,1))
    b(16) = 0.5d0*(g(2,1,1)-g(0,1,1))

    ! fy
    b(17) = 0.5d0*(g(0,1,0)-g(0,-1,0))
    b(18) = 0.5d0*(g(1,1,0)-g(1,-1,0))
    b(19) = 0.5d0*(g(0,2,0)-g(0,0,0))
    b(20) = 0.5d0*(g(1,2,0)-g(1,0,0))
    b(21) = 0.5d0*(g(0,1,1)-g(0,-1,1))
    b(22) = 0.5d0*(g(1,1,1)-g(1,-1,1))
    b(23) = 0.5d0*(g(0,2,1)-g(0,0,1))
    b(24) = 0.5d0*(g(1,2,1)-g(1,0,1))

    ! fz
    b(25) = 0.5d0*(g(0,0,1)-g(0,0,-1))
    b(26) = 0.5d0*(g(1,0,1)-g(1,0,-1))
    b(27) = 0.5d0*(g(0,1,1)-g(0,1,-1))
    b(28) = 0.5d0*(g(1,1,1)-g(1,1,-1))
    b(29) = 0.5d0*(g(0,0,2)-g(0,0,0))
    b(30) = 0.5d0*(g(1,0,2)-g(1,0,0))
    b(31) = 0.5d0*(g(0,1,2)-g(0,1,0))
    b(32) = 0.5d0*(g(1,1,2)-g(1,1,0))

    ! fxy
    b(33) = 0.25d0*(g(1,1,0)-g(-1,1,0)-g(1,-1,0)+g(-1,-1,0))
    b(34) = 0.25d0*(g(2,1,0)-g(0,1,0)-g(2,-1,0)+g(0,-1,0))
    b(35) = 0.25d0*(g(1,2,0)-g(-1,2,0)-g(1,0,0)+g(-1,0,0))
    b(36) = 0.25d0*(g(2,2,0)-g(0,2,0)-g(2,0,0)+g(0,0,0))
    b(37) = 0.25d0*(g(1,1,1)-g(-1,1,1)-g(1,-1,1)+g(-1,-1,1))
    b(38) = 0.25d0*(g(2,1,1)-g(0,1,1)-g(2,-1,1)+g(0,-1,1))
    b(39) = 0.25d0*(g(1,2,1)-g(-1,2,1)-g(1,0,1)+g(-1,0,1))
    b(40) = 0.25d0*(g(2,2,1)-g(0,2,1)-g(2,0,1)+g(0,0,1))

    ! fxz
    b(41) = 0.25d0*(g(1,0,1)-g(-1,0,1)-g(1,0,-1)+g(-1,0,-1))
    b(42) = 0.25d0*(g(2,0,1)-g(0,0,1)-g(2,0,-1)+g(0,0,-1))
    b(43) = 0.25d0*(g(1,1,1)-g(-1,1,1)-g(1,1,-1)+g(-1,1,-1))
    b(44) = 0.25d0*(g(2,1,1)-g(0,1,1)-g(2,1,-1)+g(0,1,-1))
    b(45) = 0.25d0*(g(1,0,2)-g(-1,0,2)-g(1,0,0)+g(-1,0,0))
    b(46) = 0.25d0*(g(2,0,2)-g(0,0,2)-g(2,0,0)+g(0,0,0))
    b(47) = 0.25d0*(g(1,1,2)-g(-1,1,2)-g(1,1,0)+g(-1,1,0))
    b(48) = 0.25d0*(g(2,1,2)-g(0,1,2)-g(2,1,0)+g(0,1,0))

    ! fyz
    b(49) = 0.25d0*(g(0,1,1)-g(0,-1,1)-g(0,1,-1)+g(0,-1,-1))
    b(50) = 0.25d0*(g(1,1,1)-g(1,-1,1)-g(1,1,-1)+g(1,-1,-1))
    b(51) = 0.25d0*(g(0,2,1)-g(0,0,1)-g(0,2,-1)+g(0,0,-1))
    b(52) = 0.25d0*(g(1,2,1)-g(1,0,1)-g(1,2,-1)+g(1,0,-1))
    b(53) = 0.25d0*(g(0,1,2)-g(0,-1,2)-g(0,1,0)+g(0,-1,0))
    b(54) = 0.25d0*(g(1,1,2)-g(1,-1,2)-g(1,1,0)+g(1,-1,0))
    b(55) = 0.25d0*(g(0,2,2)-g(0,0,2)-g(0,2,0)+g(0,0,0))
    b(56) = 0.25d0*(g(1,2,2)-g(1,0,2)-g(1,2,0)+g(1,0,0))

    ! fxyz
    b(57) = 0.125d0*(g(1,1,1)-g(-1,1,1)-g(1,-1,1)+g(-1,-1,1)-g(1,1,-1)+g(-1,1,-1)+g(1,-1,-1)-g(-1,-1,-1))
    b(58) = 0.125d0*(g(2,1,1)-g(0,1,1)-g(2,-1,1)+g(0,-1,1)-g(2,1,-1)+g(0,1,-1)+g(2,-1,-1)-g(0,-1,-1))
    b(59) = 0.125d0*(g(1,2,1)-g(-1,2,1)-g(1,0,1)+g(-1,0,1)-g(1,2,-1)+g(-1,2,-1)+g(1,0,-1)-g(-1,0,-1))
    b(60) = 0.125d0*(g(2,2,1)-g(0,2,1)-g(2,0,1)+g(0,0,1)-g(2,2,-1)+g(0,2,-1)+g(2,0,-1)-g(0,0,-1))
    b(61) = 0.125d0*(g(1,1,2)-g(-1,1,2)-g(1,-1,2)+g(-1,-1,2)-g(1,1,0)+g(-1,1,0)+g(1,-1,0)-g(-1,-1,0))
    b(62) = 0.125d0*(g(2,1,2)-g(0,1,2)-g(2,-1,2)+g(0,-1,2)-g(2,1,0)+g(0,1,0)+g(2,-1,0)-g(0,-1,0))
    b(63) = 0.125d0*(g(1,2,2)-g(-1,2,2)-g(1,0,2)+g(-1,0,2)-g(1,2,0)+g(-1,2,0)+g(1,0,0)-g(-1,0,0))
    b(64) = 0.125d0*(g(2,2,2)-g(0,2,2)-g(2,0,2)+g(0,0,2)-g(2,2,0)+g(0,2,0)+g(2,0,0)-g(0,0,0))

    ! calculate the coefficient vector
    a = matmul(c,b)

    ! interpolation in integer coordinates
    x = (x * f%n - (idx-1))

    l = 1 ! packed coefficient vector index
    do k = 0, 3
       do j = 0, 3
          ! horner's rule on x
          bb(j) = a(l) + x(1) * (a(l+1) + x(1) * (a(l+2) + x(1) * a(l+3)))
          bbx(j) = a(l+1) + x(1) * (2d0 * a(l+2) +  x(1) * 3d0 * a(l+3))
          bbxx(j) = 2d0 * a(l+2) + 6d0 * x(1) * a(l+3)

          ! advance
          l = l + 4
       end do
       ! horner's rule on y
       aa(k) = bb(0) + x(2) * (bb(1) + x(2) * (bb(2) + x(2) * bb(3)))

       aax(k) = bbx(0) + x(2) * (bbx(1) + x(2) * (bbx(2) + x(2) * bbx(3)))
       aay(k) = bb(1) + x(2) * (2d0 * bb(2) + x(2) * 3d0 * bb(3))
       aaxy(k) = bbx(1) + x(2) * (2d0 * bbx(2) + x(2) * 3d0 * bbx(3))

       aaxx(k) = bbxx(0) + x(2) * (bbxx(1) + x(2) * (bbxx(2) + x(2) * bbxx(3)))
       aayy(k) = 2d0 * bb(2) + 6d0 * x(2) * bb(3)
    end do

    ! field value
    y = aa(0) + x(3) * (aa(1) + x(3) * (aa(2) + x(3) * aa(3)))

    ! gradient
    yp(1) = aax(0) + x(3) * (aax(1) + x(3) * (aax(2) + x(3) * aax(3)))
    yp(2) = aay(0) + x(3) * (aay(1) + x(3) * (aay(2) + x(3) * aay(3)))
    yp(3) = aa(1) + x(3) * (2d0 * aa(2) + x(3) * 3d0 * aa(3))

    ! hessian
    ypp(1,1) = aaxx(0) + x(3) * (aaxx(1) + x(3) * (aaxx(2) + x(3) * aaxx(3)))
    ypp(1,2) = aaxy(0) + x(3) * (aaxy(1) + x(3) * (aaxy(2) + x(3) * aaxy(3)))
    ypp(1,3) = aax(1) + x(3) * (2d0 * aax(2) + x(3) * 3d0 * aax(3))
    ypp(2,2) = aayy(0) + x(3) * (aayy(1) + x(3) * (aayy(2) + x(3) * aayy(3)))
    ypp(2,3) = aay(1) + x(3) * (2d0 * aay(2) + x(3) * 3d0 * aay(3))
    ypp(3,3) = 2d0 * aa(2) + 6d0 * x(3) * aa(3)

    ! transform back to fractional coordinates and fill the Hessian
    do i = 1, 3
       yp(i) = yp(i) * f%n(i)
       do j = i, 3
          ypp(i,j) = ypp(i,j) * f%n(i) * f%n(j)
          ypp(j,i) = ypp(i,j)
       end do
    end do

  end subroutine grinterp_tricubic

  !> Pseudo-nearest grid point of a x (crystallographic) (only nearest in 
  !> orthogonal grids).
  function grid_near(f,x)

    class(grid3), intent(in) :: f !< Input grid
    real*8, intent(in) :: x(3) !< Target point (cryst. coords.)
    integer :: grid_near(3) 

    grid_near = modulo(nint(x * f%n),f%n)+1

  end function grid_near

  !> Floor grid point of a point x in crystallographic coords.
  function grid_floor(f,x)

    class(grid3), intent(in) :: f !< Input grid
    real*8, intent(in) :: x(3) !< Target point (cryst. coords.)
    integer :: grid_floor(3)

    grid_floor = modulo(floor(x * f%n),f%n)+1

  end function grid_floor

  !> Initialize the grid for the trispline interpolation. This is a
  !> modified version of the corresponding subroutine in abinit.
  subroutine init_trispline(f)
    use tools_io, only: ferror, faterr
    class(grid3), intent(inout) :: f !< Input grid

    integer :: istat
    integer :: d, nmax, i, i1, i2
    real*8, allocatable :: l(:,:), fg(:)
    real*8 :: fprev, fnext, fuse, fone

    if (allocated(f%c2)) return

    allocate(f%c2(f%n(1),f%n(2),f%n(3),3),stat=istat)
    if (istat /= 0) &
       call ferror('init_trispline','Error allocating c2 grids',faterr)

    nmax = maxval(f%n)
    allocate(l(nmax,nmax))
    allocate(fg(nmax))

    ! cholesky decomposition of the matrix: 
    ! A = 
    !  ( 4 1 0 ... 0 1 )
    !  ( 1 4 1 ... 0 0 )
    !  ( 0 1 4 ... 0 0 )
    !      .......
    !  ( 0 0 ... 1 4 1 )
    !  ( 1 0 ... 0 1 4 )
    ! that is the coefficient matrix of the system that determines
    ! the x^2 coefficients of the cubic spline (array c2 / 2).
    ! L is a lower-triangular matrix, A = L * L^t

    ! direction x->y->z
    do d = 1, 3
       nmax = f%n(d)
       l = 0d0
       l(1,1) = 2d0
       l(2,1) = 1d0 / l(1,1)
       l(2,2) = sqrt(15d0) / 2d0
       l(3,2) = 1d0 / l(2,2)
       l(nmax,1) = 0.5d0
       l(nmax,2) = - 0.5d0 / sqrt(15d0)
       do i = 3, nmax-1
          l(i,i) = sqrt(4d0 - 1d0 / l(i-1,i-1)**2)
          l(i+1,i) = 1d0 / l(i,i)
          l(nmax,i) = - l(nmax,i-1) * l(i,i-1) / l(i,i)
       end do
       l(nmax,nmax-1) = (1d0 - l(nmax,nmax-2)*l(nmax-1,nmax-2)) /&
          l(nmax-1,nmax-1)
       l(nmax,nmax) = sqrt(4d0 - sum(l(nmax,1:nmax-1)**2))
       l = l / sqrt(6d0 * nmax**2)

       ! for each of the grid points in the plane, find the spline c2
       ! (c2*x^2) coefficients
       do i1 = 1,f%n(mod(d,3)+1)
          do i2 = 1,f%n(mod(d+1,3)+1)
             select case(d)
             case(1)
                fg(1:nmax) = f%f(:,i1,i2)
             case(2)
                fg(1:nmax) = f%f(i2,:,i1)
             case(3)
                fg(1:nmax) = f%f(i1,i2,:)
             end select

             ! constant terms for the system
             fone = fg(1)
             fprev = fg(nmax)
             do i = 1, nmax-1
                fnext = fg(i+1)
                fuse = fprev + fnext - fg(i) - fg(i)
                fprev = fg(i)
                fg(i) = fuse
             end do
             fg(nmax) = fprev + fone - fg(nmax) - fg(nmax)

             ! solve by direct substitution, L^t b = c
             fg(1) = fg(1) / l(1,1)
             do i = 2, nmax-1
                fg(i) = (fg(i) - l(i,i-1) * fg(i-1)) / l(i,i)
             end do
             fg(nmax) = (fg(nmax) - &
                dot_product(l(nmax,1:nmax-1),fg(1:nmax-1))) / l(nmax,nmax)

             ! again, L a = b
             fg(nmax) = fg(nmax) / l(nmax,nmax)
             fg(nmax-1) = (fg(nmax-1) - l(nmax,nmax-1) * fg(nmax)) / l(nmax-1,nmax-1)
             do i = nmax-2, 1, -1
                fg(i) = (fg(i) - l(i+1,i) * fg(i+1) - l(nmax,i)*fg(nmax)) / l(i,i)
             end do

             ! write down the second derivatives
             select case(d)
             case(1)
                f%c2(:,i1,i2,d) = fg(1:nmax)
             case(2)
                f%c2(i2,:,i1,d) = fg(1:nmax)
             case(3)
                f%c2(i1,i2,:,d) = fg(1:nmax)
             end select

          end do
       end do

    end do

    deallocate(l,fg)

  end subroutine init_trispline

  !> Given the electron density in the isrho slot, calculate the laplacian 
  !> in islap using FFT. x2c is the crystallographic to Cartesian matrix
  !> for this grid.
  subroutine laplacian(flap,frho,x2c)
    use tools_io, only: ferror, faterr
    use tools_math, only: cross, det
    use param, only: pi
    class(grid3), intent(inout) :: flap
    type(grid3), intent(in) :: frho
    real*8, intent(in) :: x2c(3,3)

    integer :: n(3), i1, i2, i3

    complex*16 :: zaux(frho%n(1)*frho%n(2)*frho%n(3))
    real*8 :: bvec(3,3), vol
    integer :: ig, ntot
    integer :: j1, j2, j3
    integer, allocatable :: ivg(:,:), igfft(:)
    real*8, allocatable :: vgc(:,:)

    call flap%end()
    if (.not.frho%isinit) &
       call ferror('grid_laplacian','no density grid',faterr)

    ! allocate slot
    n = frho%n    
    flap%n = n
    flap%isinit = .true.
    flap%mode = mode_default
    ntot = n(1) * n(2) * n(3)

    ! bvec
    bvec(:,1) = cross(x2c(:,3),x2c(:,2))
    bvec(:,2) = cross(x2c(:,1),x2c(:,3))
    bvec(:,3) = cross(x2c(:,2),x2c(:,1))
    vol = abs(det(x2c))
    bvec = 2d0 * pi / vol * bvec

    allocate(ivg(3,ntot))
    ig = 0
    do i1 = n(1)/2-n(1)+1, n(1)/2
       do i2 = n(2)/2-n(2)+1, n(2)/2
          do i3 = n(3)/2-n(3)+1, n(3)/2
             ig = ig + 1
             ivg(1,ig)=i1
             ivg(2,ig)=i2
             ivg(3,ig)=i3
          end do
       end do
    end do
    allocate(igfft(ntot),vgc(3,ntot))
    do ig = 1, ntot
       i1=ivg(1,ig)
       i2=ivg(2,ig)
       i3=ivg(3,ig)
       if (i1.ge.0) then
          j1=i1
       else
          j1=n(1)+i1
       end if
       if (i2.ge.0) then
          j2=i2
       else
          j2=n(2)+i2
       end if
       if (i3.ge.0) then
          j3=i3
       else
          j3=n(3)+i3
       end if
       igfft(ig)=j3*n(2)*n(1)+j2*n(1)+j1+1
       vgc(:,ig)=dble(i1)*bvec(:,1)+dble(i2)*bvec(:,2)+dble(i3)*bvec(:,3)
    end do

    zaux = 0d0
    zaux = reshape(frho%f,shape(zaux))
    call cfftnd(3,n,-1,zaux)

    do ig = 1, ntot
       zaux(igfft(ig)) = -dot_product(vgc(:,ig),vgc(:,ig)) * zaux(igfft(ig))
    end do

    call cfftnd(3,n,+1,zaux)
    allocate(flap%f(n(1),n(2),n(3)))
    flap%f = real(reshape(real(zaux,8),shape(flap%f)),8)

    deallocate(igfft,vgc,ivg)

  end subroutine laplacian

  !> Calculate the gradient norm of a scalar field using FFT. x2c is
  !> the crystallographic to Cartesian matrix for this grid.
  subroutine gradrho(fgrho,frho,x2c)
    use tools_io, only: faterr, ferror
    use tools_math, only: det, cross
    use param, only: pi
    class(grid3), intent(inout) :: fgrho
    type(grid3), intent(in) :: frho
    real*8, intent(in) :: x2c(3,3)

    integer :: n(3), i, i1, i2, i3, nr1, nr2, iq1, iq2, n12
    complex*16 :: zaux(frho%n(1)*frho%n(2)*frho%n(3))
    real*8 :: bvec(3,3), vol
    real*8 :: vgc
    integer :: ig, ntot
    integer :: j1, j2, j3, igfft

    call fgrho%end()
    if (.not.frho%isinit) &
       call ferror('grid_gradgrho','no density grid',faterr)

    ! allocate slot
    n = frho%n    
    fgrho%n = n
    fgrho%isinit = .true.
    fgrho%mode = mode_default
    ntot = n(1) * n(2) * n(3)

    ! bvec
    bvec(:,1) = cross(x2c(:,3),x2c(:,2))
    bvec(:,2) = cross(x2c(:,1),x2c(:,3))
    bvec(:,3) = cross(x2c(:,2),x2c(:,1))
    vol = abs(det(x2c))
    bvec = 2d0 * pi / vol * bvec

    n12 = n(1)*n(2)
    allocate(fgrho%f(n(1),n(2),n(3)))
    fgrho%f = 0d0
    do i = 1, 3
       zaux = reshape(frho%f,shape(zaux))
       call cfftnd(3,n,-1,zaux)

       do ig = 1, ntot
          iq1 = (ig-1) / (n(2)*n(3))
          i1 = iq1 - n(1)/2 + 1
          nr1 = (ig-1) - n(2)*n(3) * iq1
          iq2 = nr1 / n(3)
          i2 = iq2 - n(2)/2 + 1
          nr2 = nr1 - n(3) * iq2
          i3 = nr2 - n(3)/2 + 1
          if (i1.ge.0) then
             j1=i1
          else
             j1=n(1)+i1
          end if
          if (i2.ge.0) then
             j2=i2
          else
             j2=n(2)+i2
          end if
          if (i3.ge.0) then
             j3=i3
          else
             j3=n(3)+i3
          end if
          igfft=j3*n(2)*n(1)+j2*n(1)+j1+1
          vgc = dble(i1)*bvec(i,1)+dble(i2)*bvec(i,2)+dble(i3)*bvec(i,3)
          zaux(igfft) = vgc * cmplx(-aimag(zaux(igfft)),dble(zaux(igfft)),8)
       end do

       call cfftnd(3,n,+1,zaux)
       fgrho%f = fgrho%f + (real(reshape(real(zaux,8),shape(fgrho%f)),8))**2
    end do
    fgrho%f = sqrt(fgrho%f)

  end subroutine gradrho

  !> Given the electron density in the isrho slot, calculate the
  !> diagonal component ix (x=1,y=2,z=3) of the Hessian using FFT. x2c
  !> is the crystallographic to Cartesian matrix for this grid.
  subroutine hxx(fxx,frho,ix,x2c)
    use tools_io, only: ferror, faterr
    use tools_math, only: det, cross
    use param, only: pi
    class(grid3), intent(inout) :: fxx
    type(grid3), intent(in) :: frho
    integer, intent(in) :: ix
    real*8, intent(in) :: x2c(3,3)

    integer :: n(3), i1, i2, i3

    complex*16 :: zaux(frho%n(1)*frho%n(2)*frho%n(3))
    real*8 :: bvec(3,3), vol
    integer :: ig, ntot
    integer :: j1, j2, j3
    integer, allocatable :: ivg(:,:), igfft(:)
    real*8, allocatable :: vgc(:,:)

    call fxx%end()
    if (.not.frho%isinit) &
       call ferror('hxx','no density grid',faterr)

    ! allocate slot
    n = frho%n    
    fxx%n = n
    fxx%isinit = .true.
    fxx%mode = mode_default
    ntot = n(1) * n(2) * n(3)

    ! bvec
    bvec(:,1) = cross(x2c(:,3),x2c(:,2))
    bvec(:,2) = cross(x2c(:,1),x2c(:,3))
    bvec(:,3) = cross(x2c(:,2),x2c(:,1))
    vol = abs(det(x2c))
    bvec = 2d0 * pi / vol * bvec

    allocate(ivg(3,ntot))
    ig = 0
    do i1 = n(1)/2-n(1)+1, n(1)/2
       do i2 = n(2)/2-n(2)+1, n(2)/2
          do i3 = n(3)/2-n(3)+1, n(3)/2
             ig = ig + 1
             ivg(1,ig)=i1
             ivg(2,ig)=i2
             ivg(3,ig)=i3
          end do
       end do
    end do
    allocate(igfft(ntot),vgc(3,ntot))
    do ig = 1, ntot
       i1=ivg(1,ig)
       i2=ivg(2,ig)
       i3=ivg(3,ig)
       if (i1.ge.0) then
          j1=i1
       else
          j1=n(1)+i1
       end if
       if (i2.ge.0) then
          j2=i2
       else
          j2=n(2)+i2
       end if
       if (i3.ge.0) then
          j3=i3
       else
          j3=n(3)+i3
       end if
       igfft(ig)=j3*n(2)*n(1)+j2*n(1)+j1+1
       vgc(:,ig)=dble(i1)*bvec(:,1)+dble(i2)*bvec(:,2)+dble(i3)*bvec(:,3)
    end do

    zaux = 0d0
    zaux = reshape(frho%f,shape(zaux))
    call cfftnd(3,n,-1,zaux)

    do ig = 1, ntot
       zaux(igfft(ig)) = -vgc(ix,ig) * vgc(ix,ig) * zaux(igfft(ig))
    end do

    call cfftnd(3,n,+1,zaux)
    allocate(fxx%f(n(1),n(2),n(3)))
    fxx%f = real(reshape(real(zaux,8),shape(fxx%f)),8)

    deallocate(igfft,vgc,ivg)

  end subroutine hxx

  !> Build a Wannier function from QEs unk(r) functions in the UNK files.
  subroutine get_qe_wnr(f,ibnd,ispin,fout)
    use tools_io, only: string, ferror, faterr, fopen_read, fopen_write,&
       fclose
    use param, only: tpi, img
    class(grid3), intent(in) :: f
    integer, intent(in) :: ibnd
    integer, intent(in) :: ispin
    complex*16, intent(out) :: fout(f%wan%nwan(1)*f%n(1),f%wan%nwan(2)*f%n(2),f%wan%nwan(3)*f%n(3))

    integer :: luc, i, j, k, nk1, nk2, nk3
    integer :: nk, ik, ik1, ik2, ik3, nall(3), naux(3), ikk, nbnd
    integer :: jbnd, imax(3)
    complex*16 :: raux(f%n(1),f%n(2),f%n(3)), rseq(f%n(1)*f%n(2)*f%n(3))
    complex*16 :: raux2(f%n(1),f%n(2),f%n(3)), raux3(f%n(1),f%n(2),f%n(3))
    character(len=:), allocatable :: fname
    complex*16 :: tnorm, ph
    complex*16, allocatable :: evc(:)
    integer :: n(3)

    n = f%n
    nk1 = f%wan%nwan(1)
    nk2 = f%wan%nwan(2)
    nk3 = f%wan%nwan(3)
    fout = 0d0
    nk = nk1 * nk2 * nk3
    nall(1) = n(1) * nk1
    nall(2) = n(2) * nk2
    nall(3) = n(3) * nk3

    if (allocated(f%wan%ngk)) then
       allocate(evc(maxval(f%wan%ngk(1:nk))))
       luc = fopen_read(f%wan%fevc,form="unformatted")
       if (ispin == 2) then
          do ik = 1, nk
             do jbnd = 1, f%wan%nbnd
                read (luc) evc(1:f%wan%ngk(ik))
             end do
          end do
       end if
    end if

    ! run over k-points
    do ik = 1, nk
       ! phase factor
       do k = 1, n(3)
          do j = 1, n(2)
             do i = 1, n(1)
                raux2(i,j,k) = exp(tpi*img*(f%wan%kpt(1,ik)*real(i-1,8)/real(n(1),8)+&
                   f%wan%kpt(2,ik)*real(j-1,8)/real(n(2),8)+f%wan%kpt(3,ik)*real(k-1,8)/real(n(3),8)))
             end do
          end do
       end do

       if (allocated(f%wan%ngk)) then
          ! from a unkgen
          if (f%wan%useu) then
             rseq = 0d0
             do jbnd = 1, f%wan%nbnd
                read (luc) evc(1:f%wan%ngk(ik))
                rseq(f%wan%nls(f%wan%igk_k(1:f%wan%ngk(ik),ik))) = &
                   rseq(f%wan%nls(f%wan%igk_k(1:f%wan%ngk(ik),ik))) + &
                   f%wan%u(jbnd,ibnd,ik) * evc(1:f%wan%ngk(ik))
             end do
          else
             rseq = 0d0
             do jbnd = 1, f%wan%nbnd
                read (luc) evc(1:f%wan%ngk(ik))
                if (jbnd == ibnd) &
                   rseq(f%wan%nls(f%wan%igk_k(1:f%wan%ngk(ik),ik))) = evc(1:f%wan%ngk(ik))
             end do
          end if
          raux = reshape(rseq,shape(raux))
          call cfftnd(3,n,+1,raux)
       else
          ! from the WNK created using the UNK files
          fname = "WNK." // string(ik) // "." // string(ibnd) // "." // string(ispin)
          luc = fopen_read(fname,form="unformatted")
          read(luc) raux
          call fclose(luc)
       end if

       ! add the contribution from this k-point
       raux2 = raux2 * raux
       do ikk = 1, nk
          ik1 = nint(f%wan%kpt(1,ikk) * nk1)
          ik2 = nint(f%wan%kpt(2,ikk) * nk2)
          ik3 = nint(f%wan%kpt(3,ikk) * nk3)
          ph = exp(tpi*img*(f%wan%kpt(1,ik)*ik1+f%wan%kpt(2,ik)*ik2+f%wan%kpt(3,ik)*ik3))
          fout(ik1*n(1)+1:(ik1+1)*n(1),ik2*n(2)+1:(ik2+1)*n(2),ik3*n(3)+1:(ik3+1)*n(3)) = &
             fout(ik1*n(1)+1:(ik1+1)*n(1),ik2*n(2)+1:(ik2+1)*n(2),ik3*n(3)+1:(ik3+1)*n(3)) + &
             raux2 * ph
       end do
    end do

    if (allocated(f%wan%ngk)) then
       deallocate(evc)
       call fclose(luc)
    end if

    ! normalize
    imax = maxloc(abs(fout))
    tnorm = fout(imax(1),imax(2),imax(3))
    tnorm = tnorm / abs(tnorm) * nk
    fout = fout / tnorm

  end subroutine get_qe_wnr

  !> Build a 3d grid using an arithmetic expression.
  subroutine new_eval(f,sptr,n,expr,fh,field_cube)
    use hashmod, only: hash
    use arithmetic, only: token, tokenize, token_undef, token_num, token_fun,&
       token_op, token_lpar, token_rpar, token_comma, token_field, iprec, iassoc,&
       istype, fun_openpar, fun_xc
    use types, only: realloc
    use iso_c_binding, only: c_ptr
    class(grid3), intent(inout) :: f
    type(c_ptr), intent(in) :: sptr
    integer, intent(in) :: n(3)
    character(*), intent(in) :: expr
    type(hash), intent(in) :: fh

    interface
       function field_cube(sptr,n,id,fder,dry,ifail) result(q)
         use iso_c_binding, only: c_ptr
         type(c_ptr), intent(in) :: sptr
         character*(*), intent(in) :: id
         integer, intent(in) :: n(3)
         character*4, intent(in) :: fder
         logical, intent(in) :: dry
         logical, intent(out) :: ifail
         real*8 :: q(n(1),n(2),n(3))
       end function field_cube
    end interface

    integer :: i, ntok, lp
    integer :: c, s(100)
    logical :: again, ok, ifail
    integer :: nq, ns, nder
    real*8, allocatable :: q(:,:,:,:)
    type(token), allocatable :: toklist(:)

    call f%end()
    f%n = n
    f%mode = mode_default
    f%isinit = .true.

    ! tokenize the expression in input
    lp = 1
    ok = tokenize(expr,ntok,toklist,lp,fh)
    if (.not.ok) then
       goto 999
       return
    end if

    ! initialize
    nq = 0
    ns = 0
    allocate(q(n(1),n(2),n(3),1))

    ! the grid version of the arithmetic evaluator does not support
    ! certain types of operators (xc, chemfunction).
    do i = 1, ntok
       if (toklist(i)%ival == fun_xc) goto 999
       if (istype(toklist(i)%ival,'chemfunction')) goto 999
       if (toklist(i)%type == token_field) then
          q(:,:,:,1) = field_cube(sptr,n,toklist(i)%sval,toklist(i)%fder,.true.,ifail)
          if (ifail) goto 999
       end if
    end do

    ! run over tokens
    do i = 1, ntok
       if (toklist(i)%type == token_num) then
          ! a number
          nq = nq + 1
          if (nq > size(q,4)) call realloc(q,n(1),n(2),n(3),nq)
          q(:,:,:,nq) = toklist(i)%fval
       elseif (toklist(i)%type == token_fun) then
          ! a function
          ns = ns + 1
          s(ns) = toklist(i)%ival
       elseif (toklist(i)%type == token_op) then
          ! a binary operator
          c = toklist(i)%ival
          again = .true.
          do while (again)
             again = .false.
             if (ns > 0) then
                if (iprec(c) < iprec(s(ns)) .or. iassoc(c)==-1 .and. iprec(c)<=iprec(s(ns))) then
                   call pop_grid(q,nq,s,ns,ifail)
                   if (ifail) goto 999
                   again = .true.
                end if
             end if
          end do
          ns = ns + 1
          s(ns) = c
       elseif (toklist(i)%type == token_lpar) then
          ! left parenthesis
          ns = ns + 1
          s(ns) = fun_openpar
       elseif (toklist(i)%type == token_rpar) then
          ! right parenthesis
           do while (ns > 0)
              if (s(ns) == fun_openpar) exit
              call pop_grid(q,nq,s,ns,ifail)
              if (ifail) goto 999
           end do
           if (ns == 0) goto 999
           ns = ns - 1
           ! if the top of the stack is a function, pop it
           if (ns > 0) then
              c = s(ns)
              if (istype(c,'function')) then
                 call pop_grid(q,nq,s,ns,ifail)
                 if (ifail) goto 999
              end if
           end if
        elseif (toklist(i)%type == token_comma) then
           ! a comma
           do while (ns > 0)
              if (s(ns) == fun_openpar) exit
              call pop_grid(q,nq,s,ns,ifail)
              if (ifail) goto 999
           end do
           if (s(ns) /= fun_openpar) goto 999
        elseif (toklist(i)%type == token_field) then
           ! a field
           nq = nq + 1
           if (nq > size(q,4)) call realloc(q,n(1),n(2),n(3),nq)
           q(:,:,:,nq) = field_cube(sptr,n,toklist(i)%sval,toklist(i)%fder,.false.,ifail)
           if (ifail) goto 999
       else
          goto 999
       end if
    end do

    ! unwind the stack
    do while (ns > 0)
       call pop_grid(q,nq,s,ns,ifail)
       if (ifail) goto 999
    end do
    allocate(f%f(n(1),n(2),n(3)))
    f%f = q(:,:,:,1)
    if (allocated(q)) deallocate(q)

    return
999 continue
    if (allocated(q)) deallocate(q)
    call f%end()

  end subroutine new_eval

  !> Pop from the stack and operate on the queue.  This routine is
  !> thread-safe. Grid version
  subroutine pop_grid(q,nq,s,ns,fail)
    use arithmetic, only: istype, fun_plus, fun_minus, fun_prod, fun_div,&
       fun_power, fun_modulo, fun_atan2, fun_min, fun_max, fun_great,&
       fun_great, fun_less, fun_leq, fun_geq, fun_equal, fun_neq, fun_and,&
       fun_or, fun_uplus, fun_uminus, fun_abs, fun_exp, fun_sqrt, fun_floor,&
       fun_ceiling, fun_round, fun_log, fun_log10, fun_sin, fun_asin, fun_cos,&
       fun_acos, fun_tan, fun_atan, fun_sinh, fun_cosh, fun_erf, fun_erfc,&
       fun_xc

    real*8, allocatable, intent(inout) :: q(:,:,:,:)
    integer, intent(inout) :: s(:)
    integer, intent(inout) :: nq, ns
    logical, intent(out) :: fail

    integer :: c

    ! pop from the stack
    fail = .false.
    if (ns == 0) goto 999
    c = s(ns)
    ns = ns - 1

    ! And apply to the queue
    if (c == fun_xc) then
       goto 999
    elseif (istype(c,'binary')) then
       ! a binary operator or function
       if (nq < 2) goto 999
       nq = nq - 1
       select case(c)
       case (fun_plus)
          q(:,:,:,nq) = q(:,:,:,nq) + q(:,:,:,nq+1)
       case (fun_minus)
          q(:,:,:,nq) = q(:,:,:,nq) - q(:,:,:,nq+1)
       case (fun_prod)
          q(:,:,:,nq) = q(:,:,:,nq) * q(:,:,:,nq+1)
       case (fun_div)
          q(:,:,:,nq) = q(:,:,:,nq) / q(:,:,:,nq+1)
       case (fun_power)
          q(:,:,:,nq) = q(:,:,:,nq) ** q(:,:,:,nq+1)
       case (fun_modulo)
          q(:,:,:,nq) = modulo(q(:,:,:,nq),q(:,:,:,nq+1))
       case (fun_atan2)
          q(:,:,:,nq) = atan2(q(:,:,:,nq),q(:,:,:,nq+1))
       case (fun_min)
          q(:,:,:,nq) = min(q(:,:,:,nq),q(:,:,:,nq+1))
       case (fun_max)
          q(:,:,:,nq) = max(q(:,:,:,nq),q(:,:,:,nq+1))
       case (fun_great)
          where (q(:,:,:,nq) > q(:,:,:,nq+1))
             q(:,:,:,nq) = 1d0
          elsewhere
             q(:,:,:,nq) = 0d0
          end where
       case (fun_less)
          where (q(:,:,:,nq) < q(:,:,:,nq+1))
             q(:,:,:,nq) = 1d0
          elsewhere
             q(:,:,:,nq) = 0d0
          end where
       case (fun_leq)
          where (q(:,:,:,nq) <= q(:,:,:,nq+1))
             q(:,:,:,nq) = 1d0
          elsewhere
             q(:,:,:,nq) = 0d0
          end where
       case (fun_geq)
          where (q(:,:,:,nq) >= q(:,:,:,nq+1))
             q(:,:,:,nq) = 1d0
          elsewhere
             q(:,:,:,nq) = 0d0
          end where
       case (fun_equal)
          where (q(:,:,:,nq) == q(:,:,:,nq+1))
             q(:,:,:,nq) = 1d0
          elsewhere
             q(:,:,:,nq) = 0d0
          end where
       case (fun_neq)
          where (q(:,:,:,nq) /= q(:,:,:,nq+1))
             q(:,:,:,nq) = 1d0
          elsewhere
             q(:,:,:,nq) = 0d0
          end where
       case (fun_and)
          where (.not.(q(:,:,:,nq)==0d0).and..not.(q(:,:,:,nq+1)==0d0))
             q(:,:,:,nq) = 1d0
          elsewhere
             q(:,:,:,nq) = 0d0
          end where
       case (fun_or)
          where (.not.(q(:,:,:,nq)==0d0).or..not.(q(:,:,:,nq+1)==0d0))
             q(:,:,:,nq) = 1d0
          elsewhere
             q(:,:,:,nq) = 0d0
          end where
       end select
    elseif (istype(c,'unary')) then
       ! a unary operator or function
       if (nq < 1) goto 999
       select case(c)
       case (fun_uplus)
          q(:,:,:,nq) = +q(:,:,:,nq)
       case (fun_uminus)
          q(:,:,:,nq) = -q(:,:,:,nq)
       case (fun_abs)
          q(:,:,:,nq) = abs(q(:,:,:,nq))
       case (fun_exp)
          q(:,:,:,nq) = exp(q(:,:,:,nq))
       case (fun_sqrt)
          q(:,:,:,nq) = sqrt(q(:,:,:,nq))
       case (fun_floor)
          q(:,:,:,nq) = floor(q(:,:,:,nq))
       case (fun_ceiling)
          q(:,:,:,nq) = ceiling(q(:,:,:,nq))
       case (fun_round)
          q(:,:,:,nq) = nint(q(:,:,:,nq))
       case (fun_log)
          q(:,:,:,nq) = log(q(:,:,:,nq))
       case (fun_log10)
          q(:,:,:,nq) = log10(q(:,:,:,nq))
       case (fun_sin)
          q(:,:,:,nq) = sin(q(:,:,:,nq))
       case (fun_asin)
          q(:,:,:,nq) = asin(q(:,:,:,nq))
       case (fun_cos)
          q(:,:,:,nq) = cos(q(:,:,:,nq))
       case (fun_acos)
          q(:,:,:,nq) = acos(q(:,:,:,nq))
       case (fun_tan)
          q(:,:,:,nq) = tan(q(:,:,:,nq))
       case (fun_atan)
          q(:,:,:,nq) = atan(q(:,:,:,nq))
       case (fun_sinh)
          q(:,:,:,nq) = sinh(q(:,:,:,nq))
       case (fun_cosh)
          q(:,:,:,nq) = cosh(q(:,:,:,nq))
       case (fun_erf)
          q(:,:,:,nq) = erf(q(:,:,:,nq))
       case (fun_erfc)
          q(:,:,:,nq) = erfc(q(:,:,:,nq))
       end select
    else
       goto 999
    end if

    return
999 continue
    fail = .true.

  end subroutine pop_grid

end module grid3mod
