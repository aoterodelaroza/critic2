!> Abinit structure and scalar field reader.
!> Most routines in this module were adapted from abinit.
! COPYRIGHT
! Copyright (C) 1998-2017 ABINIT group (DCA, XG, GMR, MKV, MT)
! This file is distributed under the terms of the
! GNU General Public License, see ~abinit/COPYING
! or http://www.gnu.org/copyleft/gpl.txt .
! For the initials of contributors, see ~abinit/doc/developers/contributor

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

! abinit tools
module abinit_private
  implicit none

  private
  public :: hdr_type, hdr_io ! abinit

  integer, parameter :: dp=kind(1.0d0)
  character(len=1), parameter :: ch10 = char(10)
  real(dp), parameter :: tol6= 0.000001_dp
  real(dp), parameter :: Ha_eV=27.21138386_dp
  integer, parameter :: md5_slen = 32

  type pawrhoij_type
     integer :: cplex
     integer :: itypat
     integer :: lmn_size
     integer :: lmn2_size
     integer :: lmnmix_sz
     integer :: ngrhoij = 0
     integer :: nrhoijsel
     integer :: nspden
     integer :: nspinor
     integer :: nsppol
     integer :: use_rhoij_
     integer :: use_rhoijp = 0
     integer :: use_rhoijres = 0
     real(dp), allocatable :: rhoijp (:,:)
     integer, allocatable :: rhoijselect(:)
  end type pawrhoij_type

  ! The hdr_type to be exported
  type hdr_type
     real(dp) :: rprimd(3,3)
     integer :: natom
     integer :: ntypat
     real(dp), allocatable :: xred(:,:)
     real(dp), allocatable :: znucltypat(:)
     integer, allocatable :: typat(:)
     integer :: ngfft(3)
  end type hdr_type

  ! Old version of the hdr_type (headform = 23,34,40,41,42,44,53,56,57)
  type hdr_type_1
     integer :: bantot        ! total number of bands (sum of nband on all kpts and spins)
     integer :: date          ! starting date
     integer :: headform      ! format of the header
     integer :: intxc,ixc,natom,nkpt,npsp,nspden    ! input variables
     integer :: nspinor,nsppol,nsym,ntypat,occopt   ! input variables
     integer :: pertcase      ! the index of the perturbation, 0 if GS calculation
     integer :: usepaw        ! input variable (0=norm-conserving psps, 1=paw)
     integer :: usewvl        ! input variable (0=plane-waves, 1=wavelets)
     integer :: ngfft(3)      ! input variable
     integer, allocatable :: istwfk(:)    ! input variable istwfk(nkpt)
     integer, allocatable :: lmn_size(:)  ! lmn_size(npsp) from psps
     integer, allocatable :: nband(:)     ! input variable nband(nkpt*nsppol)
     integer, allocatable :: npwarr(:)    ! npwarr(nkpt) array holding npw for each k point
     integer          :: nwvlarr(2)   ! nwvlarr(2) array holding the number of wavelets
     integer, allocatable :: pspcod(:)    ! pscod(npsp) from psps
     integer, allocatable :: pspdat(:)    ! psdat(npsp) from psps
     integer, allocatable :: pspso(:)     ! pspso(npsp) from psps
     integer, allocatable :: pspxc(:)     ! pspxc(npsp) from psps
     integer, allocatable :: so_psp(:)    ! input variable so_psp(npsp)
     integer, allocatable :: symafm(:)    ! input variable symafm(nsym)
     integer, allocatable :: symrel(:,:,:)! input variable symrel(3,3,nsym)
     integer, allocatable :: typat(:)     ! input variable typat(natom)
     real(dp) :: ecut                  ! input variable
     real(dp) :: ecutdg                ! input variable (ecut for NC psps, pawecutdg for paw)
     real(dp) :: ecutsm                ! input variable
     real(dp) :: ecut_eff              ! ecut*dilatmx**2 (dilatmx is an input variable)
     real(dp) :: etot,fermie,residm    ! EVOLVING variables
     real(dp) :: qptn(3)               ! the wavevector, in case of a perturbation
     real(dp) :: rprimd(3,3)           ! EVOLVING variables
     real(dp) :: stmbias               ! input variable
     real(dp) :: tphysel               ! input variable
     real(dp) :: tsmear                ! input variable
     real(dp), allocatable :: kptns(:,:)   ! input variable kptns(3,nkpt)
     real(dp), allocatable :: occ(:)       ! EVOLVING variable occ(bantot)
     real(dp), allocatable :: tnons(:,:)   ! input variable tnons(3,nsym)
     real(dp), allocatable :: wtk(:)       ! weight of kpoints wtk(nkpt)
     real(dp), allocatable :: xred(:,:)    ! EVOLVING variable xred(3,natom)
     real(dp), allocatable :: zionpsp(:)   ! zionpsp(npsp) from psps
     real(dp), allocatable :: znuclpsp(:)  ! znuclpsp(npsp) from psps
     real(dp), allocatable :: znucltypat(:)! znucltypat(ntypat) from alchemy
     character(len=6) :: codvsn              ! version of the code
     character(len=132), allocatable :: title(:) ! title(npsp) from psps
     type(pawrhoij_type), allocatable :: pawrhoij(:) ! EVOLVING variable, only for paw
  end type hdr_type_1

  type hdr_type_2
     integer :: bantot        ! total number of bands (sum of nband on all kpts and spins)
     integer :: date          ! starting date
     integer :: headform      ! format of the header
     integer :: intxc         ! input variable
     integer :: ixc           ! input variable
     integer :: mband         ! maxval(hdr%nband)
     integer :: natom         ! input variable
     integer :: nkpt          ! input variable
     integer :: npsp          ! input variable
     integer :: nspden        ! input variable
     integer :: nspinor       ! input variable
     integer :: nsppol        ! input variable
     integer :: nsym          ! input variable
     integer :: ntypat        ! input variable
     integer :: occopt        ! input variable
     integer :: pertcase      ! the index of the perturbation, 0 if GS calculation
     integer :: usepaw        ! input variable (0=norm-conserving psps, 1=paw)
     integer :: usewvl        ! input variable (0=plane-waves, 1=wavelets)
     integer :: kptopt        ! input variable (defines symmetries used for k-point sampling)
     integer :: pawcpxocc     ! input variable
     integer :: nshiftk_orig=1 ! original number of shifts given in input (changed in inkpts, the actual value is nshiftk)
     integer :: nshiftk=1       ! number of shifts after inkpts.
     integer :: icoulomb        ! input variable.
     real(dp) :: ecut         ! input variable
     real(dp) :: ecutdg       ! input variable (ecut for NC psps, pawecutdg for paw)
     real(dp) :: ecutsm       ! input variable
     real(dp) :: ecut_eff     ! ecut*dilatmx**2 (dilatmx is an input variable)
     real(dp) :: etot         ! EVOLVING variable
     real(dp) :: fermie       ! EVOLVING variable
     real(dp) :: residm       ! EVOLVING variable
     real(dp) :: stmbias      ! input variable
     real(dp) :: tphysel      ! input variable
     real(dp) :: tsmear       ! input variable
     real(dp) :: nelect       ! number of electrons (computed from pseudos and charge)
     real(dp) :: charge       ! input variable
     real(dp) :: qptn(3)
     real(dp) :: rprimd(3,3)
     integer :: ngfft(3)
     integer :: nwvlarr(2)
     integer :: kptrlatt_orig(3,3)
     integer :: kptrlatt(3,3)
     integer, allocatable :: istwfk(:)
     integer, allocatable :: lmn_size(:)
     integer, allocatable :: nband(:)
     integer, allocatable :: npwarr(:)
     integer, allocatable :: pspcod(:)
     integer, allocatable :: pspdat(:)
     integer, allocatable :: pspso(:)
     integer, allocatable :: pspxc(:)
     integer, allocatable :: so_psp(:)
     integer, allocatable :: symafm(:)
     integer, allocatable :: symrel(:,:,:)
     integer, allocatable :: typat(:)
     real(dp), allocatable :: kptns(:,:)
     real(dp), allocatable :: occ(:)
     real(dp), allocatable :: tnons(:,:)
     real(dp), allocatable :: wtk(:)
     real(dp),allocatable :: shiftk_orig(:,:)
     real(dp),allocatable :: shiftk(:,:)
     real(dp),allocatable :: amu(:)
     real(dp), allocatable :: xred(:,:)
     real(dp), allocatable :: zionpsp(:)
     real(dp), allocatable :: znuclpsp(:)
     real(dp), allocatable :: znucltypat(:)
     character(len=6) :: codvsn
     character(len=132), allocatable :: title(:)
     character(len=md5_slen),allocatable :: md5_pseudos(:)
     type(pawrhoij_type), allocatable :: pawrhoij(:)
  end type hdr_type_2

  interface
     module subroutine hdr_io(fform,hdr,rdwr,unitfi,errmsg)
       integer, intent(inout) :: fform
       integer, intent(in) :: rdwr,unitfi
       type(hdr_type), intent(inout) :: hdr
       character(len=:), allocatable, intent(out) :: errmsg
     end subroutine hdr_io
     module subroutine hdr_io_1(fform,hdr,unitfi,errmsg)
       integer,intent(inout) :: fform
       integer,intent(in) :: unitfi
       type(hdr_type_1),intent(inout) :: hdr
       character(len=:), allocatable, intent(out) :: errmsg
     end subroutine hdr_io_1
     module subroutine hdr_io_2(fform,hdr,unit,errmsg)
       integer, intent(out) :: fform
       integer, intent(in) :: unit
       type(hdr_type_2), intent(out) :: hdr
       character(len=:), allocatable, intent(out) :: errmsg
     end subroutine hdr_io_2
  end interface

end module abinit_private
