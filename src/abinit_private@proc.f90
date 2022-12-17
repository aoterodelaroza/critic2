!> Abinit structure and scalar field reader.
!> Most routines in this module were adapted from abinit.
! COPYRIGHT
! Copyright (C) 1998-2017 ABINIT group (DCA, XG, GMR, MKV, MT)
! This file is distributed under the terms of the
! GNU General Public License, see ~abinit/COPYING
! or http://www.gnu.org/copyleft/gpl.txt .
! For the initials of contributors, see ~abinit/doc/developers/contributor

! Copyright (c) 2007-2022 Alberto Otero de la Roza <aoterodelaroza@gmail.com>,
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

! abinit tools and readers
submodule (abinit_private) proc
  implicit none

  ! private parameters
  integer, parameter :: dp=kind(1.0d0)
  character(len=1), parameter :: ch10 = char(10)
  real(dp), parameter :: tol6= 0.000001_dp
  real(dp), parameter :: Ha_eV=27.21138386_dp
  integer, parameter :: md5_slen = 32

  ! paw info
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

  !xx! private subroutines
  ! subroutine hdr_io_1(fform,hdr,unitfi,headform_1,errmsg)
  ! subroutine hdr_io_2(fform,hdr,unit,headform_2,errmsg)
  ! subroutine pawrhoij_io(pawrhoij,unitfi,nspden_in,typat,headform,errmsg)

contains

  ! Driver to choose which version of hdr_io to use.
  module subroutine hdr_io(fform,hdr,rdwr,unitfi,errmsg)
    use tools_io, only: string
    integer, intent(inout) :: fform
    integer, intent(in) :: rdwr,unitfi
    type(hdr_type), intent(inout) :: hdr
    character(len=:), allocatable, intent(out) :: errmsg

    type(hdr_type_1) :: hdr_1
    type(hdr_type_2) :: hdr_2
    integer :: headform_1, headform_2

    errmsg = ""
    ! read the headform
    if (rdwr /= 1 .and. rdwr /= 5) then
       errmsg = "Only read operations allowed in abinit files."
       return
    end if
    if (rdwr == 1) then
       rewind(unitfi)
    end if

    call hdr_io_1(fform,hdr_1,unitfi,headform_1,errmsg)
    if (len_trim(errmsg) > 0) return
    if (headform_1 /= 0) then
       rewind(unitfi)
       call hdr_io_2(fform,hdr_2,unitfi,headform_2,errmsg)
       if (headform_2 /= 0) then
          errmsg = "Can not handle headform: " // string(headform_2) // &
             ". The version of abinit you are using is not supported &
             &by critic2. Please, e-mail the critic2 developer this &
             &message and the version of abinit you are using."
          return
       end if
       if (len_trim(errmsg) > 0) return
       if (allocated(hdr%xred)) deallocate(hdr%xred)
       if (allocated(hdr%znucltypat)) deallocate(hdr%znucltypat)
       if (allocated(hdr%typat)) deallocate(hdr%typat)
       allocate(hdr%xred(3,hdr_2%natom))
       allocate(hdr%znucltypat(hdr_2%ntypat))
       allocate(hdr%typat(hdr_2%natom))
       hdr%rprimd = hdr_2%rprimd
       hdr%ntypat = hdr_2%ntypat
       hdr%natom = hdr_2%natom
       hdr%xred = hdr_2%xred
       hdr%znucltypat = hdr_2%znucltypat
       hdr%typat = hdr_2%typat
       hdr%ngfft = hdr_2%ngfft
    else
       if (allocated(hdr%xred)) deallocate(hdr%xred)
       if (allocated(hdr%znucltypat)) deallocate(hdr%znucltypat)
       if (allocated(hdr%typat)) deallocate(hdr%typat)
       allocate(hdr%xred(3,hdr_1%natom))
       allocate(hdr%znucltypat(hdr_1%ntypat))
       allocate(hdr%typat(hdr_1%natom))
       hdr%rprimd = hdr_1%rprimd
       hdr%ntypat = hdr_1%ntypat
       hdr%natom = hdr_1%natom
       hdr%xred = hdr_1%xred
       hdr%znucltypat = hdr_1%znucltypat
       hdr%typat = hdr_1%typat
       hdr%ngfft = hdr_1%ngfft
    end if

  end subroutine hdr_io

  !xx! private subroutines

  ! The following functions are from ABINIT
  ! Copyright (C) 2007-2009 ABINIT group (MT)
  ! This file is distributed under the terms of the
  ! GNU General Public License, see ~ABINIT/Infos/copyright
  ! or http://www.gnu.org/copyleft/gpl.txt .

  !> The hdr_io subroutine from abinit.
  subroutine hdr_io_1(fform,hdr,unitfi,headform_1,errmsg)
    use tools_io, only: string
    !! This subroutine deals with the I/O of the hdr_type
    !! structured variables (read/write/echo).
    !! According to the value of rdwr, it reads the header
    !! of a file, writes it, or echo the value of the structured
    !! variable to a file.
    !! Note that, when reading, different records of hdr
    !! are allocated here, according to the values of the
    !! read variables. Records of hdr should be deallocated
    !! correctly by a call to hdr_clean when hdr is not used anymore.
    !! Two instances of the hdr_io routines are defined :
    !!  hdr_io_int to which only the unit number is given
    !!  hdr_io_wfftype to which a wffil datatype is given
    !!
    !! COPYRIGHT
    !! Copyright (C) 2002-2009 ABINIT group (XG,MB)
    !! This file is distributed under the terms of the
    !! GNU General Public License, see ~abinit/COPYING
    !! or http://www.gnu.org/copyleft/gpl.txt .
    !! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
    !!
    !! INPUTS
    !!  rdwr= if 1, read the hdr structured variable from the header of the file,
    !!        if 2, write the header to unformatted file
    !!        if 3, echo part of the header to formatted file (records 1 and 2)
    !!        if 4, echo the header to formatted file
    !!        if 5, read the hdr without rewinding (unformatted)
    !!        if 6, write the hdr without rewinding (unformatted)
    !!  unitfi=unit number of the file (unformatted if rdwr=1, 2, 5 or 6 formatted if rdwr=3,4)
    !!
    !! OUTPUT
    !!  (see side effects)
    !!
    !! SIDE EFFECTS
    !!  The following variables are both input or output :
    !!  fform=kind of the array in the file
    !!   if rdwr=1,5 : will be output ; if the reading fail, return fform=0
    !!   if rdwr=2,3,4,6 : should be input, will be written or echo to file
    !!  hdr <type(hdr_type)>=the header structured variable
    !!   if rdwr=1,5 : will be output
    !!   if rdwr=2,3,4,6 : should be input, will be written or echo to file
    !!
    !! NOTES
    !! In all cases, the file is supposed to be open already
    !! When reading (rdwr=1) or writing (rdwr=2), rewind the file
    !! When echoing (rdwr=3) does not rewind the file.
    !! When reading (rdwr=5) or writing (rdwr=6), DOES NOT rewind the file
    !!
    !! PARENTS
    !!      conducti,cut3d,initaim,inwffil3,ioarr,macroave,newsp,outkss,outwf,print_ij
    !!      rdem1,rdkss,testepsm1,testlda,uderiv,vtorho3,wrem1
    !!
    !! CHILDREN
    !!      leave_new,wrtout
    !!
    !! SOURCE

    !    use defs_basis
    !    use defs_datatypes
    !
    !    use m_IO_tools, only : flush_unit

    !This section has been created automatically by the script Abilint (TD).
    !Do not modify the following lines by hand.
    !    use interfaces_14_hidewrite
    !    use interfaces_16_hideleave
    !    use interfaces_32_util
    !End of the abilint section

    !Arguments ------------------------------------
    integer,intent(inout) :: fform
    integer,intent(in) :: unitfi
    type(hdr_type_1),intent(inout) :: hdr
    integer, intent(out) :: headform_1
    character(len=:), allocatable, intent(out) :: errmsg

    !Local variables-------------------------------
    integer :: bantot,bsize,cplex,headform,iatom,ierr,ipsp,ispden
    integer :: lloc,lmax,mmax,natom,nkpt,npsp,nspden,nsppol
    integer :: nsym,ntypat
    character(len=6) :: codvsn
    integer, allocatable :: ibuffer(:),nsel44(:,:),nsel56(:)
    real(dp) :: acell(3)
    real(dp), allocatable :: buffer(:)
    ! *************************************************************************

    headform_1 = 0
    ! Reading the first record of the file ------------------------------------

    read(unitfi,iostat=ierr)codvsn,fform
    if (ierr /=0) then
       fform=0
       headform = -1
       errmsg = "Error reading file."
       return   ! This is to allow treatment of old epsm1 format
    end if

    if(fform==1   .or. &
       fform==2   .or. &
       fform==51  .or. &
       fform==52  .or. &
       fform==101 .or. &
       fform==102       )then
       !  This is the old format
       headform=22
    else

       !  Format beyond 22 have a different first line, so need reading again the first line
       backspace (unitfi)
       read (unitfi)   codvsn,headform,fform

       if(headform/=23 .and. &
          headform/=34 .and. &
          headform/=40 .and. &
          headform/=41 .and. &
          headform/=42 .and. &
          headform/=44 .and. &
          headform/=53 .and. &
          headform/=56 .and. &
          headform/=57)then
          headform_1 = headform
          return
       end if
    end if

    hdr%codvsn=codvsn
    hdr%headform=headform
    ! fform is not a record of hdr_type

    !DEBUG
    ! write(6,*)' hdr_io : debug '
    ! write(6,*)' hdr_io : codvsn,headform,fform',codvsn,headform,fform
    !ENDDEBUG

    ! Reading the second record of the file ------------------------------------

    ! Initialize the values that are not present for all versions (exception : npsp)
    hdr%nspden=1
    hdr%nspinor=1
    hdr%occopt=1
    hdr%pertcase=1
    hdr%usepaw=0
    hdr%usewvl=0
    hdr%ecut=0d0
    hdr%ecutdg=0d0
    hdr%ecutsm=0d0
    hdr%qptn(1:3)=0d0
    hdr%stmbias=0d0
    hdr%tphysel=0d0
    hdr%tsmear=0d0

    errmsg = "Error reading file."
    if(headform==22)then

       read(unitfi,err=999,end=999) bantot, hdr%date, hdr%intxc, hdr%ixc, natom, hdr%ngfft(1:3),&
          nkpt, nsppol, nsym, ntypat,&
          acell, hdr%ecut_eff, hdr%rprimd
       npsp=ntypat

    else if(headform==23)then

       !   Compared to v2.2, add nspden, nspinor, occopt
       read(unitfi,err=999,end=999) bantot, hdr%date, hdr%intxc, hdr%ixc, natom, hdr%ngfft(1:3),&
          nkpt, hdr%nspden, hdr%nspinor, nsppol, nsym, ntypat, hdr%occopt,&
          acell, hdr%ecut_eff, hdr%rprimd
       npsp=ntypat

    else if(headform==34)then

       !   Compared to v2.3, subtract acell, and add npsp
       read(unitfi,err=999,end=999) bantot, hdr%date, hdr%intxc, hdr%ixc, natom, hdr%ngfft(1:3),&
          nkpt, hdr%nspden, hdr%nspinor, nsppol, nsym, npsp, ntypat, hdr%occopt,&
          hdr%ecut_eff, hdr%rprimd

    else if(headform==40)then

       !   Compared to v3.4, add ecut, ecutsm, tphysel, tsmear
       read(unitfi,err=999,end=999) bantot, hdr%date, hdr%intxc, hdr%ixc, natom, hdr%ngfft(1:3),&
          nkpt, hdr%nspden, hdr%nspinor, nsppol, nsym, npsp, ntypat, hdr%occopt,&
          hdr%ecut, hdr%ecutsm, hdr%ecut_eff, hdr%rprimd, hdr%tphysel, hdr%tsmear

    else if(headform==41)then

       !   Compared to v4.0, add pertcase and qptn(3)
       read(unitfi,err=999,end=999) bantot, hdr%date, hdr%intxc, hdr%ixc, natom, hdr%ngfft(1:3),&
          nkpt, hdr%nspden, hdr%nspinor, nsppol, nsym, npsp, ntypat, hdr%occopt, hdr%pertcase,&
          hdr%ecut, hdr%ecutsm, hdr%ecut_eff, hdr%qptn(1:3), hdr%rprimd, hdr%tphysel, hdr%tsmear

    else if(headform==42)then

       !   Compared to v4.1, add stmbias
       read(unitfi,err=999,end=999) bantot, hdr%date, hdr%intxc, hdr%ixc, natom, hdr%ngfft(1:3),&
          nkpt, hdr%nspden, hdr%nspinor, nsppol, nsym, npsp, ntypat, hdr%occopt, hdr%pertcase,&
          hdr%ecut, hdr%ecutsm, hdr%ecut_eff, hdr%qptn(1:3), hdr%rprimd,&
          hdr%stmbias, hdr%tphysel, hdr%tsmear

    else if(headform>=44 .and. headform<57)then

       !   Compared to v4.2, add usepaw and ecutdg
       read(unitfi,err=999,end=999) bantot, hdr%date, hdr%intxc, hdr%ixc, natom, hdr%ngfft(1:3),&
          nkpt, hdr%nspden, hdr%nspinor, nsppol, nsym, npsp, ntypat, hdr%occopt, hdr%pertcase,&
          hdr%usepaw, hdr%ecut, hdr%ecutdg, hdr%ecutsm, hdr%ecut_eff, hdr%qptn(1:3), hdr%rprimd,&
          hdr%stmbias, hdr%tphysel, hdr%tsmear

    else if(headform>=57)then

       !   Compared to v4.4, add usewvl
       read(unitfi,err=999,end=999) bantot, hdr%date, hdr%intxc, hdr%ixc, natom, hdr%ngfft(1:3),&
          nkpt, hdr%nspden, hdr%nspinor, nsppol, nsym, npsp, ntypat, hdr%occopt, hdr%pertcase,&
          hdr%usepaw, hdr%ecut, hdr%ecutdg, hdr%ecutsm, hdr%ecut_eff, hdr%qptn(1:3), hdr%rprimd,&
          hdr%stmbias, hdr%tphysel, hdr%tsmear, hdr%usewvl

    end if

    hdr%bantot=bantot
    hdr%natom =natom
    hdr%nkpt  =nkpt
    hdr%npsp  =npsp
    hdr%nsppol=nsppol
    hdr%nsym  =nsym
    hdr%ntypat =ntypat

    if(hdr%ecutsm>tol6 .and. headform<44 .and. .not.(fform==51.or.fform==52.or.fform==101.or.fform==102))then
       errmsg = "The value of ecutsm is "//string(hdr%ecutsm,'f',8,2)//", while the file has been produced prior to v4.4 ."
       return
    end if

    !DEBUG
    !  write(6,*)' hdr_io : before allocate '
    !  write(6,*)' hdr_io : bantot,natom,nkpt,npsp,nsppol,nsym,ntypat',&
    !&  bantot,natom,nkpt,npsp,nsppol,nsym,ntypat
    !ENDDEBUG

    ! Allocate all parts of hdr that need to be --------------------------------
    allocate(hdr%istwfk(nkpt))
    allocate(hdr%kptns(3,nkpt))
    allocate(hdr%lmn_size(npsp))
    allocate(hdr%nband(nkpt*nsppol))
    allocate(hdr%npwarr(nkpt)) ! Warning : npwarr here has only one dim
    allocate(hdr%occ(bantot))
    allocate(hdr%pspcod(npsp))
    allocate(hdr%pspdat(npsp))
    allocate(hdr%pspso(npsp))
    allocate(hdr%pspxc(npsp))
    allocate(hdr%so_psp(npsp))
    allocate(hdr%symafm(nsym))
    allocate(hdr%symrel(3,3,nsym))
    allocate(hdr%title(npsp))
    allocate(hdr%tnons(3,nsym))
    allocate(hdr%typat(natom))
    allocate(hdr%wtk(nkpt))
    allocate(hdr%xred(3,natom))
    allocate(hdr%zionpsp(npsp))
    allocate(hdr%znuclpsp(npsp))
    allocate(hdr%znucltypat(ntypat))
    if(hdr%usepaw==1) allocate(hdr%pawrhoij(natom))

    !DEBUG
    ! write(6,*)' hdr_io : after allocate '
    !ENDDEBUG

    ! Reading the third record of the file ------------------------------------

    ! Initialize the values that are not present for all versions
    hdr%istwfk(:)=1
    hdr%so_psp(:)=1
    hdr%symafm(:)=1

    if(headform==22 .and. (fform==1 .or. fform==51 .or. fform==101))then

       !   This is very old (pre-2.0) format !
       read(unitfi,err=999,end=999) hdr%nband(:), hdr%npwarr(:), hdr%symrel(:,:,:), &
          hdr%typat(:), hdr%kptns(:,:), hdr%occ(:), &
          hdr%tnons(:,:), hdr%znucltypat(:)

    else if(headform==22 .or. headform==23 .or. headform==34)then

       !   Compared to pre v2.0, add istwfk
       read(unitfi,err=999,end=999) hdr%nband(:), hdr%npwarr(:), hdr%symrel(:,:,:), &
          hdr%typat(:), hdr%istwfk(:), hdr%kptns(:,:), hdr%occ(:), &
          hdr%tnons(:,:), hdr%znucltypat(:)

    else if(headform>=40 .and. headform < 50)then

       !   Compared to pre v4.0, add so_psp and symafm, and switch istwfk

       read(unitfi,err=999,end=999)  hdr%istwfk(:), hdr%nband(:), hdr%npwarr(:), &
          hdr%so_psp(:), hdr%symafm(:), hdr%symrel(:,:,:), &
          hdr%typat(:), hdr%kptns(:,:), hdr%occ(:), &
          hdr%tnons(:,:), hdr%znucltypat(:)

    else if(headform>=50)then

       !   Compared to pre v5.0, add wtk
       read(unitfi,err=999,end=999)  hdr%istwfk(:), hdr%nband(:), hdr%npwarr(:), &
          hdr%so_psp(:), hdr%symafm(:), hdr%symrel(:,:,:), &
          hdr%typat(:), hdr%kptns(:,:), hdr%occ(:), &
          hdr%tnons(:,:), hdr%znucltypat(:), hdr%wtk(:)

    end if

    ! Reading the records with psp information ---------------------------------

    ! Initialize the values that are not present for all versions
    hdr%pspso(:)=1
    hdr%lmn_size(:)=0

    if(headform==22)then

       do ipsp=1,npsp
          read(unitfi,err=999,end=999) hdr%title(ipsp), hdr%znuclpsp(ipsp), &
             hdr%zionpsp(ipsp), hdr%pspdat(ipsp), hdr%pspcod(ipsp), &
             hdr%pspxc(ipsp), lmax, lloc, mmax
       end do

    else if(headform==23)then

       !  Compared to 2.2, add pspso
       do ipsp=1,npsp
          read(unitfi,err=999,end=999) hdr%title(ipsp), hdr%znuclpsp(ipsp), &
             hdr%zionpsp(ipsp), hdr%pspso(ipsp), hdr%pspdat(ipsp), &
             hdr%pspcod(ipsp), hdr%pspxc(ipsp), lmax, lloc, mmax
       end do

    else if(headform==34 .or. headform==40 .or. headform==41 &
       .or. headform==42)then

       !  Compared to 2.3, suppress lmax, lloc, mmax
       do ipsp=1,npsp
          read(unitfi,err=999,end=999) hdr%title(ipsp), hdr%znuclpsp(ipsp), &
             hdr%zionpsp(ipsp), hdr%pspso(ipsp), hdr%pspdat(ipsp), &
             hdr%pspcod(ipsp), hdr%pspxc(ipsp)
       end do

    else if(headform>=44)then

       !  Compared to 4.2, add lmn_size
       do ipsp=1,npsp
          read(unitfi,err=999,end=999) hdr%title(ipsp), hdr%znuclpsp(ipsp), &
             hdr%zionpsp(ipsp), hdr%pspso(ipsp), hdr%pspdat(ipsp), &
             hdr%pspcod(ipsp), hdr%pspxc(ipsp), hdr%lmn_size(ipsp)
       end do

    end if

    ! Reading the final record of the header  ---------------------------------

    ! Initialize the values that are not present for all versions
    hdr%fermie=0d0

    if(headform==22)then
       read(unitfi,err=999,end=999) hdr%residm, hdr%xred(:,:), hdr%etot
    else if(headform==23 .or. headform==34 .or. headform>=40)then
       read(unitfi,err=999,end=999) hdr%residm, hdr%xred(:,:), hdr%etot, hdr%fermie
    end if

    !DEBUG
    !  write(6,*)' hdr_io : read mode, hdr%so_psp(:), hdr%symafm(:)=',&
    !&                                 hdr%so_psp(:), hdr%symafm(:)
    !ENDDEBUG

    ! Reading the Rhoij tab if the PAW method was used  -----------------------
    if (hdr%usepaw==1) then
       if ((headform>=44).and.(headform<56)) then
          allocate(nsel44(hdr%nspden,hdr%natom))
          read(unitfi,err=999,end=999) ((nsel44(ispden,iatom),ispden=1,hdr%nspden),iatom=1,hdr%natom)
          bsize=sum(nsel44)
          allocate(ibuffer(bsize),buffer(bsize));
          read(unitfi,err=999,end=999) ibuffer(:),buffer(:)
          deallocate(ibuffer,buffer,nsel44)
       else if (headform>=56) then
          allocate(nsel56(hdr%natom))
          if (headform==56) then
             read(unitfi,err=999,end=999) (nsel56(iatom),iatom=1,hdr%natom),cplex
             nspden=hdr%nspden
          else
             read(unitfi,err=999,end=999) (nsel56(iatom),iatom=1,hdr%natom),cplex,nspden
          end if
          bsize=sum(nsel56)
          allocate(ibuffer(bsize),buffer(bsize*nspden*cplex))
          read(unitfi,err=999,end=999) ibuffer(:),buffer(:)
          deallocate(ibuffer,buffer,nsel56)
       end if
    end if

    errmsg = ""
999 continue

  end subroutine hdr_io_1

  subroutine hdr_io_2(fform,hdr,unit,headform_2,errmsg)
    implicit none
    integer, intent(out) :: fform
    integer, intent(in) :: unit
    type(hdr_type_2), intent(out) :: hdr
    integer, intent(out) :: headform_2
    character(len=:), allocatable, intent(out) :: errmsg

    !Local variables-------------------------------
    integer :: ipsp
    character(len=500) :: errmsgl
    real(dp),allocatable :: occ3d(:,:,:)
    integer :: ii,band,ikpt,spin
    character(len=:), allocatable :: errmsg2

    !*************************************************************************

    errmsg = "Error reading file."
    ! Reading the first record of the file ------------------------------------
    ! fform is not a record of hdr_type
    read(unit, err=10, iomsg=errmsgl) hdr%codvsn,hdr%headform,fform

    headform_2 = 0
    if (hdr%headform < 80) then
       headform_2 = hdr%headform
       errmsg = ""
       return
    end if

    !Reading the second record of the file ------------------------------------
    read(unit, err=10, iomsg=errmsgl) &
       hdr%bantot, hdr%date, hdr%intxc, hdr%ixc, hdr%natom, hdr%ngfft(1:3),&
       hdr%nkpt, hdr%nspden, hdr%nspinor, hdr%nsppol, hdr%nsym, hdr%npsp, hdr%ntypat, hdr%occopt, hdr%pertcase,&
       hdr%usepaw, hdr%ecut, hdr%ecutdg, hdr%ecutsm, hdr%ecut_eff, hdr%qptn(1:3), hdr%rprimd,&
       hdr%stmbias, hdr%tphysel, hdr%tsmear, hdr%usewvl, hdr%nshiftk_orig, hdr%nshiftk, hdr%mband


    !Allocate all parts of hdr that need to be --------------------------------
    allocate(hdr%istwfk(hdr%nkpt))
    allocate(hdr%nband(hdr%nkpt*hdr%nsppol))
    allocate(hdr%npwarr(hdr%nkpt))
    allocate(hdr%pspcod(hdr%npsp))
    allocate(hdr%pspdat(hdr%npsp))
    allocate(hdr%pspso(hdr%npsp))
    allocate(hdr%pspxc(hdr%npsp))
    allocate(hdr%lmn_size(hdr%npsp))
    allocate(hdr%so_psp(hdr%npsp))
    allocate(hdr%symafm(hdr%nsym))
    allocate(hdr%symrel(3,3,hdr%nsym))
    allocate(hdr%typat(hdr%natom))
    allocate(hdr%kptns(3,hdr%nkpt))
    allocate(hdr%occ(hdr%bantot))
    allocate(hdr%tnons(3,hdr%nsym))
    allocate(hdr%wtk(hdr%nkpt))
    allocate(hdr%xred(3,hdr%natom))
    allocate(hdr%zionpsp(hdr%npsp))
    allocate(hdr%znuclpsp(hdr%npsp))
    allocate(hdr%znucltypat(hdr%ntypat))
    allocate(hdr%title(hdr%npsp))
    allocate(hdr%shiftk(3,hdr%nshiftk))
    allocate(hdr%shiftk_orig(3,hdr%nshiftk_orig))
    allocate(hdr%md5_pseudos(hdr%npsp))
    allocate(hdr%amu(hdr%ntypat))

    if (hdr%usepaw==1)  then
       allocate(hdr%pawrhoij(hdr%natom))
    end if

    ! Reading the third record of the file ------------------------------------

    ! Take into account future migration to occ(:,:,:) in the Format
    ! read 3d matrix with stride and transfer to (stupid) 1d hdr%occ in packed form.
    allocate(occ3d(hdr%mband,hdr%nkpt,hdr%nsppol))

    read(unit, err=10, iomsg=errmsgl) &
       hdr%istwfk(:), hdr%nband(:), hdr%npwarr(:), &
       hdr%so_psp(:), hdr%symafm(:), hdr%symrel(:,:,:), &
       hdr%typat(:), hdr%kptns(:,:), occ3d, &
       hdr%tnons(:,:), hdr%znucltypat(:), hdr%wtk(:)

    ii = 0
    do spin=1,hdr%nsppol
       do ikpt=1,hdr%nkpt
          do band=1,hdr%nband(ikpt + (spin-1) * hdr%nkpt)
             ii = ii +1
             hdr%occ(ii) = occ3d(band,ikpt,spin)
          end do
       end do
    end do
    deallocate(occ3d)

    ! Reading the final record of the header  ---------------------------------
    read(unit, err=10, iomsg=errmsgl) hdr%residm, hdr%xred(:,:), hdr%etot, hdr%fermie, hdr%amu(:)

    read(unit, err=10, iomsg=errmsgl)&
       hdr%kptopt,hdr%pawcpxocc,hdr%nelect,hdr%charge,hdr%icoulomb,&
       hdr%kptrlatt,hdr%kptrlatt_orig, hdr%shiftk_orig,hdr%shiftk

    ! Reading the records with psp information ---------------------------------
    do ipsp=1,hdr%npsp
       read(unit, err=10, iomsg=errmsgl) &
          hdr%title(ipsp), hdr%znuclpsp(ipsp), hdr%zionpsp(ipsp), hdr%pspso(ipsp), hdr%pspdat(ipsp), &
          hdr%pspcod(ipsp), hdr%pspxc(ipsp), hdr%lmn_size(ipsp), hdr%md5_pseudos(ipsp)
    end do

    if (hdr%usepaw==1) then ! Reading the Rhoij tab if the PAW method was used.
       call pawrhoij_io(hdr%pawrhoij,unit,hdr%nspden,hdr%typat,hdr%headform,errmsg2)
       if (len_trim(errmsg2) > 0) then
          errmsg = errmsg2
          goto 10
       end if
    end if

    errmsg = ""
10  continue

  end subroutine hdr_io_2

  subroutine pawrhoij_io(pawrhoij,unitfi,nspden_in,typat,headform,errmsg)
    implicit none
    !Arguments ------------------------------------
    !scalars
    integer,intent(in) :: unitfi,headform,nspden_in
    !arrays
    integer,intent(in) :: typat(:)
    type(pawrhoij_type),intent(inout) :: pawrhoij(:)
    character(len=:), allocatable, intent(out) :: errmsg

    !Local variables-------------------------------
    !scalars
    integer :: iatom,natom,ispden,bsize
    integer :: my_cplex,my_natom,my_nspden
    !arrays
    integer,allocatable :: ibuffer(:),nsel44(:,:),nsel56(:)
    real(dp), allocatable :: buffer(:)

    ! *************************************************************************

    errmsg = "Error reading file."
    my_natom=SIZE(pawrhoij)
    if (my_natom==0) goto 999
    my_nspden=nspden_in
    natom=size(typat)


    if ((headform>=44).and.(headform<56)) then
       allocate(nsel44(nspden_in,natom))
       read(unitfi,err=999,end=999) ((nsel44(ispden,iatom),ispden=1,nspden_in),iatom=1,natom)
       bsize=sum(nsel44)
       allocate(ibuffer(bsize),buffer(bsize))
       read(unitfi,err=999,end=999) ibuffer(:),buffer(:)
       deallocate(ibuffer,buffer,nsel44)
    else if (headform>=56) then
       allocate(nsel56(natom))
       if (headform==56) then
          read(unitfi,err=999,end=999) (nsel56(iatom),iatom=1,natom),my_cplex
       else
          read(unitfi,err=999,end=999) (nsel56(iatom),iatom=1,natom),my_cplex,my_nspden
       end if
       bsize=sum(nsel56)
       allocate(ibuffer(bsize),buffer(bsize*my_nspden*my_cplex))
       read(unitfi,err=999,end=999) ibuffer(:),buffer(:)
       deallocate(ibuffer,buffer,nsel56)
    end if

    errmsg = ""
999 continue

  end subroutine pawrhoij_io

end submodule proc
