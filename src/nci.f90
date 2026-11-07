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
!

! non-covalent interaction index plots
module nci
  use fragmentmod, only: fragment
  implicit none

  private
  public :: nciplot
  public :: nci_options, nci_result
  public :: nci_options_init, nci_options_set_field, nci_dopromol, nci_compute, nci_masks
  public :: nci_compute_promol, nci_frag_region
  public :: nci_in_dat, nci_domains
  public :: nci_write_cubes, nci_write_dat, nci_write_vmd, nci_write_files
  public :: nci_chk_read, nci_chk_write

  ! value of the reduced density gradient at the grid points removed
  ! from the isosurface plot by the cutoffs
  real*8, parameter, public :: nci_rdg_masked = 100d0

  !> Options of an NCI calculation: the density, the grid, the
  !> fragments, and the cutoffs. nci_options_init sets the defaults
  !> for a density; the grid (x0, xmat, nstep, periodic) is set by the
  !> caller before nci_compute. The cutoffs are not used by
  !> nci_compute: they can be changed and the masks and files redone
  !> without recalculating.
  type nci_options
     ! the density
     integer :: ifield = 0 ! field whose density is used (0 = promolecular)
     logical :: fourierint = .false. ! derivatives of a grid field by Fourier transform (vs. differentiating the interpolant)
     ! the grid: nstep points along each column of xmat, starting at x0 (Cartesian, cell frame)
     logical :: periodic = .true. ! the grid spans whole unit cells (crystals; the plot wraps around)
     real*8 :: x0(3) = 0d0 ! origin
     real*8 :: xmat(3,3) = 0d0 ! step vectors (columns)
     integer :: nstep(3) = 0 ! number of points along each step vector (the files list the atoms inside
                              ! the box around the grid, non-periodic grids only)
     ! the fragments: only the regions between them are plotted
     integer :: nfrag = 0 ! number of fragments
     type(fragment), allocatable :: fr(:) ! the fragments
     logical :: frcell = .false. ! fragments are sets of cell atoms (by cidx), with their lattice translations
     real*8 :: rhoparam = 0.95d0 ! no fragment contributes more than this fraction of the fragments' promolecular density
     real*8 :: rhoparam2 = 0.75d0 ! the fragments contribute at least this fraction of the total promolecular density
     logical :: dovoid = .false. ! plot only the voids: where the promolecular density is below rho_void
     real*8 :: rho_void = 0.01d0 ! promolecular density below which a point is in a void (dovoid)
     ! the cutoffs
     real*8 :: rhocut = 0.2d0 ! density cutoff of the s vs. rho points (dat file)
     real*8 :: dimcut = 2d0 ! RDG cutoff of the s vs. rho points (dat file)
     real*8 :: rhoplot = 0.1d0 ! density cutoff of the isosurface plot (and half-width of its color scale)
     real*8 :: dimplot = 0.5d0 ! RDG isovalue of the isosurface plot
     real*8 :: srhorange(2) = (/-1d30,1d30/) ! range of the signed density kept in the isosurface plot
     logical :: onlyneg = .false. ! plot only where the second eigenvalue of the Hessian is negative
  end type nci_options

  !> Grids of an NCI calculation (nci_compute) and the density and
  !> grid they were calculated on, in grid order: f(i,j,k) is the
  !> value at x0 + (i-1)*xmat(:,1) + (j-1)*xmat(:,2) + (k-1)*xmat(:,3).
  type nci_result
     integer :: ifield = -1 ! field whose density was used (-1 = no result)
     logical :: fourierint = .false. ! whether the derivatives were calculated by Fourier transform
     real*8 :: x0(3) = 0d0 ! grid origin
     real*8 :: xmat(3,3) = 0d0 ! grid step vectors (columns)
     integer :: n(3) = 0 ! grid dimensions
     real*8, allocatable :: srho(:,:,:) ! density times the sign of the second Hessian eigenvalue
     real*8, allocatable :: rdg(:,:,:) ! reduced density gradient
     real*8, allocatable :: rhoat(:,:,:) ! promolecular density (if nci_dopromol)
     real*8, allocatable :: rhofrag(:,:,:,:) ! promolecular density of each fragment (nfrag > 0 only)
  end type nci_result

  interface
     module subroutine nciplot()
     end subroutine nciplot
     module subroutine nci_options_init(opts,s,ifield)
       use systemmod, only: system
       type(nci_options), intent(out) :: opts
       type(system), intent(in) :: s
       integer, intent(in) :: ifield
     end subroutine nci_options_init
     module subroutine nci_options_set_field(opts,s,ifield)
       use systemmod, only: system
       type(nci_options), intent(inout) :: opts
       type(system), intent(in) :: s
       integer, intent(in) :: ifield
     end subroutine nci_options_set_field
     pure module function nci_dopromol(opts) result(dopromol)
       type(nci_options), intent(in) :: opts
       logical :: dopromol
     end function nci_dopromol
     module subroutine nci_compute(opts,s,res,errmsg,pre)
       use systemmod, only: system
       type(nci_options), intent(in) :: opts
       type(system), intent(inout), target :: s
       type(nci_result), intent(out) :: res
       character(len=:), allocatable, intent(out) :: errmsg
       type(nci_result), intent(inout), optional :: pre
     end subroutine nci_compute
     module subroutine nci_compute_promol(opts,s,res,errmsg)
       use systemmod, only: system
       type(nci_options), intent(in) :: opts
       type(system), intent(inout), target :: s
       type(nci_result), intent(out) :: res
       character(len=:), allocatable, intent(out) :: errmsg
     end subroutine nci_compute_promol
     module subroutine nci_frag_region(opts,res,g)
       type(nci_options), intent(in) :: opts
       type(nci_result), intent(in) :: res
       real*8, allocatable, intent(out) :: g(:,:,:)
     end subroutine nci_frag_region
     module subroutine nci_masks(opts,res,inter,splot,voidq,voidqp,voidv,errmsg)
       type(nci_options), intent(in) :: opts
       type(nci_result), intent(in) :: res
       logical(1), allocatable, intent(out), optional :: inter(:,:,:)
       real*8, allocatable, intent(out) :: splot(:,:,:)
       real*8, intent(out) :: voidq, voidqp, voidv
       character(len=:), allocatable, intent(out) :: errmsg
     end subroutine nci_masks
     pure module function nci_in_dat(opts,res,inter,i,j,k) result(ok)
       type(nci_options), intent(in) :: opts
       type(nci_result), intent(in) :: res
       logical(1), intent(in) :: inter(:,:,:)
       integer, intent(in) :: i, j, k
       logical :: ok
     end function nci_in_dat
     module subroutine nci_domains(opts,res,splot,ndom,lbl,npts,vol,rmin,rmax,errmsg)
       type(nci_options), intent(in) :: opts
       type(nci_result), intent(in) :: res
       real*8, intent(in) :: splot(:,:,:)
       integer, intent(out) :: ndom
       integer, allocatable, intent(out) :: lbl(:,:,:)
       integer, allocatable, intent(out) :: npts(:)
       real*8, allocatable, intent(out) :: vol(:), rmin(:), rmax(:)
       character(len=:), allocatable, intent(out) :: errmsg
     end subroutine nci_domains
     module subroutine nci_write_cubes(oname,opts,s,res,splot,errmsg)
       use systemmod, only: system
       character*(*), intent(in) :: oname
       type(nci_options), intent(in) :: opts
       type(system), intent(in) :: s
       type(nci_result), intent(in) :: res
       real*8, intent(in) :: splot(:,:,:)
       character(len=:), allocatable, intent(out) :: errmsg
     end subroutine nci_write_cubes
     module subroutine nci_write_dat(file,opts,res,inter,errmsg)
       character*(*), intent(in) :: file
       type(nci_options), intent(in) :: opts
       type(nci_result), intent(in) :: res
       logical(1), intent(in) :: inter(:,:,:)
       character(len=:), allocatable, intent(out) :: errmsg
     end subroutine nci_write_dat
     module subroutine nci_write_vmd(oname,opts,s,domolmotif,errmsg,refroot)
       use systemmod, only: system
       character*(*), intent(in) :: oname
       type(nci_options), intent(in) :: opts
       type(system), intent(inout) :: s
       logical, intent(in) :: domolmotif
       character(len=:), allocatable, intent(out) :: errmsg
       character*(*), intent(in), optional :: refroot
     end subroutine nci_write_vmd
     module subroutine nci_write_files(root,opts,s,res,inter,splot,domolmotif,errmsg,refroot)
       use systemmod, only: system
       character*(*), intent(in) :: root
       type(nci_options), intent(in) :: opts
       type(system), intent(inout) :: s
       type(nci_result), intent(in) :: res
       logical(1), intent(in) :: inter(:,:,:)
       real*8, intent(in) :: splot(:,:,:)
       logical, intent(in) :: domolmotif
       character(len=:), allocatable, intent(out) :: errmsg
       character*(*), intent(in), optional :: refroot
     end subroutine nci_write_files
     module function nci_chk_read(file,opts,s,res) result(ok)
       use systemmod, only: system
       character*(*), intent(in) :: file
       type(nci_options), intent(in) :: opts
       type(system), intent(in) :: s
       type(nci_result), intent(out) :: res
       logical :: ok
     end function nci_chk_read
     module subroutine nci_chk_write(file,opts,s,res,errmsg)
       use systemmod, only: system
       character*(*), intent(in) :: file
       type(nci_options), intent(in) :: opts
       type(system), intent(in) :: s
       type(nci_result), intent(in) :: res
       character(len=:), allocatable, intent(out) :: errmsg
     end subroutine nci_chk_write
  end interface

end module nci
