!> Some routines in this module were adapted from abinit.
!>   COPYRIGHT
!>   Copyright (C) 2007-2009 ABINIT group (MT)
!>   This file is distributed under the terms of the
!>   GNU General Public License, see ~ABINIT/Infos/copyright
!>   or http://www.gnu.org/copyleft/gpl.txt .

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

  public :: hdr_type, hdr_io ! abinit

  integer, parameter :: dp=kind(1.0d0)
  character(len=1), parameter :: ch10 = char(10)
  real(dp), parameter :: tol6= 0.000001_dp
  real(dp), parameter :: Ha_eV=27.21138386_dp ! 1 Hartree, in eV
  ! defs_datatypes.F90
  type pawrhoij_type
     !Integer scalars
     integer :: cplex
     ! cplex=1 if rhoij are real, 2 if rhoij are complex
     integer :: lmn_size
     ! Number of (l,m,n) elements for the paw basis
     integer :: lmn2_size
     ! lmn2_size=lmn_size*(lmn_size+1)/2
     ! where lmn_size is the number of (l,m,n) elements for the paw basis
     integer :: lmnmix_sz
     ! lmnmix_sz=number of (lmn,lmn_prime) verifying l<=lmix and l_prime<=lmix
     !           i.e. number of rhoij elements being mixed during SCF cycle
     ! lmnmix_sz=0 if mixing data are note used
     integer :: ngrhoij
     ! First dimension of array grhoij
     integer :: nrhoijsel
     ! nrhoijsel
     ! Number of non-zero value of rhoij
     ! This is the size of rhoijp(:,:) (see below in this datastructure)
     integer :: nspden
     ! Number of spin-density components for rhoij (may be different from nspden for density)
     integer :: nsppol
     ! Number of independant spin-components
     integer :: use_rhoij_
     ! 1 if pawrhoij%rhoij_ is allocated
     integer :: use_rhoijres
     ! 1 if pawrhoij%rhoijres is allocated
     !Integer arrays
     integer, pointer :: kpawmix(:)
     ! kpawmix(lmnmix_sz)
     ! Indirect array selecting the elements of rhoij
     ! being mixed during SCF cycle
     integer, pointer :: rhoijselect(:)
     ! rhoijselect(lmn2_size)
     ! Indirect array selecting the non-zero elements of rhoij:
     ! rhoijselect(isel,ispden)=klmn if rhoij(klmn,ispden) is non-zero
     !Real (real(dp)) arrays
     real(dp), pointer :: grhoij (:,:,:)
     ! grhoij(ngrhoij,cplex*lmn2_size,nspden)
     ! Gradients of Rho_ij wrt xred, strains, ... (non-packed storage)
     real(dp), pointer :: rhoij_ (:,:)
     ! rhoij_(cplex*lmn2_size,nspden)
     ! Array used to (temporary) store Rho_ij in a non-packed storage mode
     real(dp), pointer :: rhoijp (:,:)
     ! rhoijp(cplex*lmn2_size,nspden)
     ! Augmentation waves occupancies Rho_ij
     ! in PACKED STORAGE (only non-zero elements are stored)
     real(dp), pointer :: rhoijres (:,:)
     ! rhoijres(cplex*lmn2_size,nspden)
     ! Rho_ij residuals during SCF cycle (non-packed storage)
  end type pawrhoij_type

  type hdr_type
     integer :: bantot        ! total number of bands (sum of nband on all kpts and spins)
     integer :: date          ! starting date
     integer :: headform      ! format of the header
     integer :: intxc,ixc,natom,nkpt,npsp,nspden    ! input variables
     integer :: nspinor,nsppol,nsym,ntypat,occopt   ! input variables
     integer :: pertcase      ! the index of the perturbation, 0 if GS calculation
     integer :: usepaw        ! input variable (0=norm-conserving psps, 1=paw)
     integer :: usewvl        ! input variable (0=plane-waves, 1=wavelets)
     integer :: ngfft(3)      ! input variable

     ! This record is not a part of the hdr_type, although it is present in the
     ! header of the files. This is because it depends on the kind of file
     ! that is written, while all other information does not depend on it.
     ! It was preferred to let it be initialized or defined outside of hdr_type.
     ! integer :: fform         ! file descriptor (or file format)

     integer, pointer :: istwfk(:)    ! input variable istwfk(nkpt)
     integer, pointer :: lmn_size(:)  ! lmn_size(npsp) from psps
     integer, pointer :: nband(:)     ! input variable nband(nkpt*nsppol)
     integer, pointer :: npwarr(:)    ! npwarr(nkpt) array holding npw for each k point
     integer          :: nwvlarr(2)   ! nwvlarr(2) array holding the number of wavelets
     ! for each resolution.
     integer, pointer :: pspcod(:)    ! pscod(npsp) from psps
     integer, pointer :: pspdat(:)    ! psdat(npsp) from psps
     integer, pointer :: pspso(:)     ! pspso(npsp) from psps
     integer, pointer :: pspxc(:)     ! pspxc(npsp) from psps
     integer, pointer :: so_psp(:)    ! input variable so_psp(npsp)
     integer, pointer :: symafm(:)    ! input variable symafm(nsym)
     integer, pointer :: symrel(:,:,:)! input variable symrel(3,3,nsym)
     integer, pointer :: typat(:)     ! input variable typat(natom)

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
     real(dp), pointer :: kptns(:,:)   ! input variable kptns(3,nkpt)
     real(dp), pointer :: occ(:)       ! EVOLVING variable occ(bantot)
     real(dp), pointer :: tnons(:,:)   ! input variable tnons(3,nsym)
     real(dp), pointer :: wtk(:)       ! weight of kpoints wtk(nkpt)
     real(dp), pointer :: xred(:,:)    ! EVOLVING variable xred(3,natom)
     real(dp), pointer :: zionpsp(:)   ! zionpsp(npsp) from psps
     real(dp), pointer :: znuclpsp(:)  ! znuclpsp(npsp) from psps
     ! Note the difference between znucl and znuclpsp !!
     real(dp), pointer :: znucltypat(:)! znucltypat(ntypat) from alchemy

     character(len=6) :: codvsn              ! version of the code
     character(len=132), pointer :: title(:) ! title(npsp) from psps

     type(pawrhoij_type), pointer :: pawrhoij(:) ! EVOLVING variable, only for paw

     !Should make a list of supplementary infos
     ! MG: For postprocessing purposes, it is quite useful to
     !  have kptrlatt as well as nshiftk and shiftk. also kptopt is useful
     !  to know if time reversal can be employed
  end type hdr_type

contains

  ! private
  
  ! The following functions are from ABINIT
  ! Copyright (C) 2007-2009 ABINIT group (MT)
  ! This file is distributed under the terms of the
  ! GNU General Public License, see ~ABINIT/Infos/copyright
  ! or http://www.gnu.org/copyleft/gpl.txt .
  
  !> The hdr_io subroutine from abinit.
  subroutine hdr_io(fform,hdr,rdwr,unitfi)
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
    !!      leave_new,rhoij_alloc,wrtout
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
    integer,intent(in) :: rdwr,unitfi
    type(hdr_type),intent(inout) :: hdr

    !Local variables-------------------------------
    integer :: bantot,bsize,cplex,headform,iatom,ierr,ii,ikpt,ipsp,ispden,isym
    integer :: jj,lloc,lmax,mmax,natinc,natom,nkpt,npsp,nselect,nspden,nsppol
    integer :: nsym,ntypat
    character(len=500) :: message
    character(len=6) :: codvsn
    integer, allocatable :: ibuffer(:),nsel44(:,:),nsel56(:)
    real(dp) :: acell(3)
    real(dp), allocatable :: buffer(:)

    ! *************************************************************************

    !DEBUG
    !write(6,*)' hdr_io : enter hdr_io_int '
    ! call flush(6)
    !ENDDEBUG

    ! -------------------------------------------------------------------------
    ! Reading the header of an unformatted file
    ! -------------------------------------------------------------------------

    if(rdwr==1 .or. rdwr==5)then

       if (rdwr==1) then
          rewind(unitfi)
       end if

       ! Reading the first record of the file ------------------------------------

       read(unitfi,iostat=ierr)codvsn,fform
       if (ierr /=0) then
          fform=0
          return   ! This is to allow treatment of old epsm1 format
       end if

       if(fform==1   .or. &
          &    fform==2   .or. &
          &    fform==51  .or. &
          &    fform==52  .or. &
          &    fform==101 .or. &
          &    fform==102       )then
          !  This is the old format
          headform=22

       else

          !  Format beyond 22 have a different first line, so need reading again the first line
          backspace (unitfi)
          read (unitfi)   codvsn,headform,fform

          if(headform/=23 .and. &
             &     headform/=34 .and. &
             &     headform/=40 .and. &
             &     headform/=41 .and. &
             &     headform/=42 .and. &
             &     headform/=44 .and. &
             &     headform/=53 .and. &
             &     headform/=56 .and. &
             &     headform/=57         )then
             write(message, '(4a,i3,3a,i8,3a)' ) ch10,&
                &    ' hdr_io : ERROR -',ch10,&
                &    '  The first line of the (WF, DEN or POT) file read in unit ',unitfi,' is erroneous.',ch10,&
                &    '  headform is ',headform,', while it should be 23, 34, 40, 41, 42, 44, 53 or 56 or 57.',ch10,&
                &    '  Action : check the correctness of your file.'
             call wrtout(message)
             stop 1
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

       if(headform==22)then

          read(unitfi) bantot, hdr%date, hdr%intxc, hdr%ixc, natom, hdr%ngfft(1:3),&
             &    nkpt, nsppol, nsym, ntypat,&
             &    acell, hdr%ecut_eff, hdr%rprimd
          npsp=ntypat

       else if(headform==23)then

          !   Compared to v2.2, add nspden, nspinor, occopt
          read(unitfi) bantot, hdr%date, hdr%intxc, hdr%ixc, natom, hdr%ngfft(1:3),&
             &    nkpt, hdr%nspden, hdr%nspinor, nsppol, nsym, ntypat, hdr%occopt,&
             &    acell, hdr%ecut_eff, hdr%rprimd
          npsp=ntypat

       else if(headform==34)then

          !   Compared to v2.3, subtract acell, and add npsp
          read(unitfi) bantot, hdr%date, hdr%intxc, hdr%ixc, natom, hdr%ngfft(1:3),&
             &    nkpt, hdr%nspden, hdr%nspinor, nsppol, nsym, npsp, ntypat, hdr%occopt,&
             &    hdr%ecut_eff, hdr%rprimd

       else if(headform==40)then

          !   Compared to v3.4, add ecut, ecutsm, tphysel, tsmear
          read(unitfi) bantot, hdr%date, hdr%intxc, hdr%ixc, natom, hdr%ngfft(1:3),&
             &    nkpt, hdr%nspden, hdr%nspinor, nsppol, nsym, npsp, ntypat, hdr%occopt,&
             &    hdr%ecut, hdr%ecutsm, hdr%ecut_eff, hdr%rprimd, hdr%tphysel, hdr%tsmear

       else if(headform==41)then

          !   Compared to v4.0, add pertcase and qptn(3)
          read(unitfi) bantot, hdr%date, hdr%intxc, hdr%ixc, natom, hdr%ngfft(1:3),&
             &    nkpt, hdr%nspden, hdr%nspinor, nsppol, nsym, npsp, ntypat, hdr%occopt, hdr%pertcase,&
             &    hdr%ecut, hdr%ecutsm, hdr%ecut_eff, hdr%qptn(1:3), hdr%rprimd, hdr%tphysel, hdr%tsmear

       else if(headform==42)then

          !   Compared to v4.1, add stmbias
          read(unitfi) bantot, hdr%date, hdr%intxc, hdr%ixc, natom, hdr%ngfft(1:3),&
             &    nkpt, hdr%nspden, hdr%nspinor, nsppol, nsym, npsp, ntypat, hdr%occopt, hdr%pertcase,&
             &    hdr%ecut, hdr%ecutsm, hdr%ecut_eff, hdr%qptn(1:3), hdr%rprimd,&
             &    hdr%stmbias, hdr%tphysel, hdr%tsmear

       else if(headform>=44 .and. headform<57)then

          !   Compared to v4.2, add usepaw and ecutdg
          read(unitfi) bantot, hdr%date, hdr%intxc, hdr%ixc, natom, hdr%ngfft(1:3),&
             &    nkpt, hdr%nspden, hdr%nspinor, nsppol, nsym, npsp, ntypat, hdr%occopt, hdr%pertcase,&
             &    hdr%usepaw, hdr%ecut, hdr%ecutdg, hdr%ecutsm, hdr%ecut_eff, hdr%qptn(1:3), hdr%rprimd,&
             &    hdr%stmbias, hdr%tphysel, hdr%tsmear

       else if(headform>=57)then

          !   Compared to v4.4, add usewvl
          read(unitfi) bantot, hdr%date, hdr%intxc, hdr%ixc, natom, hdr%ngfft(1:3),&
             &    nkpt, hdr%nspden, hdr%nspinor, nsppol, nsym, npsp, ntypat, hdr%occopt, hdr%pertcase,&
             &    hdr%usepaw, hdr%ecut, hdr%ecutdg, hdr%ecutsm, hdr%ecut_eff, hdr%qptn(1:3), hdr%rprimd,&
             &    hdr%stmbias, hdr%tphysel, hdr%tsmear, hdr%usewvl

       end if

       hdr%bantot=bantot
       hdr%natom =natom
       hdr%nkpt  =nkpt
       hdr%npsp  =npsp
       hdr%nsppol=nsppol
       hdr%nsym  =nsym
       hdr%ntypat =ntypat

       if(hdr%ecutsm>tol6 .and. headform<44 .and. .not.(fform==51.or.fform==52.or.fform==101.or.fform==102))then
          write(message, '(4a,es16.6,9a)' ) ch10,&
             &   ' hdr_io : ERROR -',ch10,&
             &   '  The value of ecutsm is',hdr%ecutsm,', while the file has been produced prior to v4.4 .',ch10,&
             &   '  The definition of the smearing function has changed, so that you are not allowed',ch10,&
             &   '  to restart from a old wavefunction file. By contrast, you can restart from an old',ch10,&
             &   '  potential or density file, and perform a self-consistent cycle with a new ABINIT version.',ch10,&
             &   '  Action : produce a density or potential file using the old version of ABINIT, and restart from it.'
          call wrtout(message)
          stop 1
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
          read(unitfi) hdr%nband(:), hdr%npwarr(:), hdr%symrel(:,:,:), &
             &    hdr%typat(:), hdr%kptns(:,:), hdr%occ(:), &
             &    hdr%tnons(:,:), hdr%znucltypat(:)

       else if(headform==22 .or. headform==23 .or. headform==34)then

          !   Compared to pre v2.0, add istwfk
          read(unitfi) hdr%nband(:), hdr%npwarr(:), hdr%symrel(:,:,:), &
             &    hdr%typat(:), hdr%istwfk(:), hdr%kptns(:,:), hdr%occ(:), &
             &    hdr%tnons(:,:), hdr%znucltypat(:)

       else if(headform>=40 .and. headform < 50)then

          !   Compared to pre v4.0, add so_psp and symafm, and switch istwfk

          read(unitfi)  hdr%istwfk(:), hdr%nband(:), hdr%npwarr(:), &
             &    hdr%so_psp(:), hdr%symafm(:), hdr%symrel(:,:,:), &
             &    hdr%typat(:), hdr%kptns(:,:), hdr%occ(:), &
             &    hdr%tnons(:,:), hdr%znucltypat(:)

       else if(headform>=50)then

          !   Compared to pre v5.0, add wtk
          read(unitfi)  hdr%istwfk(:), hdr%nband(:), hdr%npwarr(:), &
             &    hdr%so_psp(:), hdr%symafm(:), hdr%symrel(:,:,:), &
             &    hdr%typat(:), hdr%kptns(:,:), hdr%occ(:), &
             &    hdr%tnons(:,:), hdr%znucltypat(:), hdr%wtk(:)

       end if

       ! Reading the records with psp information ---------------------------------

       ! Initialize the values that are not present for all versions
       hdr%pspso(:)=1
       hdr%lmn_size(:)=0

       if(headform==22)then

          do ipsp=1,npsp
             read(unitfi) hdr%title(ipsp), hdr%znuclpsp(ipsp), &
                &    hdr%zionpsp(ipsp), hdr%pspdat(ipsp), hdr%pspcod(ipsp), &
                &    hdr%pspxc(ipsp), lmax, lloc, mmax
          end do

       else if(headform==23)then

          !  Compared to 2.2, add pspso
          do ipsp=1,npsp
             read(unitfi) hdr%title(ipsp), hdr%znuclpsp(ipsp), &
                &    hdr%zionpsp(ipsp), hdr%pspso(ipsp), hdr%pspdat(ipsp), &
                &    hdr%pspcod(ipsp), hdr%pspxc(ipsp), lmax, lloc, mmax
          end do

       else if(headform==34 .or. headform==40 .or. headform==41 &
          &                      .or. headform==42)then

          !  Compared to 2.3, suppress lmax, lloc, mmax
          do ipsp=1,npsp
             read(unitfi) hdr%title(ipsp), hdr%znuclpsp(ipsp), &
                &    hdr%zionpsp(ipsp), hdr%pspso(ipsp), hdr%pspdat(ipsp), &
                &    hdr%pspcod(ipsp), hdr%pspxc(ipsp)
          end do

       else if(headform>=44)then

          !  Compared to 4.2, add lmn_size
          do ipsp=1,npsp
             read(unitfi) hdr%title(ipsp), hdr%znuclpsp(ipsp), &
                &    hdr%zionpsp(ipsp), hdr%pspso(ipsp), hdr%pspdat(ipsp), &
                &    hdr%pspcod(ipsp), hdr%pspxc(ipsp), hdr%lmn_size(ipsp)
          end do

       end if

       ! Reading the final record of the header  ---------------------------------

       ! Initialize the values that are not present for all versions
       hdr%fermie=0d0

       if(headform==22)then
          read(unitfi) hdr%residm, hdr%xred(:,:), hdr%etot
       else if(headform==23 .or. headform==34 .or. headform>=40)then
          read(unitfi) hdr%residm, hdr%xred(:,:), hdr%etot, hdr%fermie
       end if

       !DEBUG
       !  write(6,*)' hdr_io : read mode, hdr%so_psp(:), hdr%symafm(:)=',&
       !&                                 hdr%so_psp(:), hdr%symafm(:)
       !ENDDEBUG

       ! Reading the Rhoij tab if the PAW method was used  -----------------------
       if (hdr%usepaw==1) then


          if ((headform>=44).and.(headform<56)) then
             allocate(nsel44(hdr%nspden,hdr%natom))
             read(unitfi) ((nsel44(ispden,iatom),ispden=1,hdr%nspden),iatom=1,hdr%natom)
             call rhoij_alloc(1,hdr%lmn_size,hdr%nspden,hdr%nsppol,hdr%pawrhoij,hdr%typat)
             do iatom=1,hdr%natom
                hdr%pawrhoij(iatom)%nrhoijsel=nsel44(1,iatom)
             end do
             bsize=sum(nsel44);allocate(ibuffer(bsize),buffer(bsize));ii=0
             read(unitfi) ibuffer(:),buffer(:)
             do iatom=1,hdr%natom
                nselect=nsel44(1,iatom)
                hdr%pawrhoij(iatom)%rhoijselect(1:nselect)=ibuffer(ii+1:ii+nselect)
                do ispden=1,hdr%nspden
                   hdr%pawrhoij(iatom)%rhoijp(1:nselect,ispden)=buffer(ii+1:ii+nselect)
                   ii=ii+nselect
                end do
             end do
             deallocate(ibuffer,buffer,nsel44)

          else if (headform>=56) then
             allocate(nsel56(hdr%natom))
             if (headform==56) then
                read(unitfi) (nsel56(iatom),iatom=1,hdr%natom),cplex
                nspden=hdr%nspden
             else
                read(unitfi) (nsel56(iatom),iatom=1,hdr%natom),cplex,nspden
             end if
             call rhoij_alloc(cplex,hdr%lmn_size,nspden,hdr%nsppol,hdr%pawrhoij,hdr%typat)
             do iatom=1,hdr%natom
                hdr%pawrhoij(iatom)%nrhoijsel=nsel56(iatom)
             end do
             bsize=sum(nsel56);allocate(ibuffer(bsize),buffer(bsize*nspden*cplex))
             ii=0;jj=0
             read(unitfi) ibuffer(:),buffer(:)
             do iatom=1,hdr%natom
                nselect=nsel56(iatom)
                hdr%pawrhoij(iatom)%rhoijselect(1:nselect)=ibuffer(ii+1:ii+nselect)
                ii=ii+nselect
                do ispden=1,nspden
                   hdr%pawrhoij(iatom)%rhoijp(1:cplex*nselect,ispden)=buffer(jj+1:jj+cplex*nselect)
                   jj=jj+cplex*nselect
                end do
             end do
             deallocate(ibuffer,buffer,nsel56)
          end if

       end if


       ! -------------------------------------------------------------------------
       ! Writing the header of an unformatted file
       ! -------------------------------------------------------------------------

    else if(rdwr==2 .or. rdwr==6)then

       ! natom,nkpt,npsp,ntypat... are not defined in this section :
       ! always address them from hdr

       if(rdwr==2) then
          rewind(unitfi)
       end if

       ! Writing always use last format version
       headform=57
       write(unitfi) hdr%codvsn, headform, fform

       write(unitfi) hdr%bantot, hdr%date, hdr%intxc, hdr%ixc, &
          &  hdr%natom, hdr%ngfft(1:3), hdr%nkpt, &
          &  hdr%nspden, hdr%nspinor, &
          &  hdr%nsppol, hdr%nsym, hdr%npsp, hdr%ntypat, hdr%occopt, hdr%pertcase,&
          &  hdr%usepaw, hdr%ecut, hdr%ecutdg, hdr%ecutsm, hdr%ecut_eff, &
          &  hdr%qptn, hdr%rprimd, hdr%stmbias, hdr%tphysel, hdr%tsmear, &
          &  hdr%usewvl

       write(unitfi) hdr%istwfk(:), hdr%nband(:), hdr%npwarr(:),&
          &   hdr%so_psp(:), hdr%symafm(:), hdr%symrel(:,:,:), &
          &   hdr%typat(:), hdr%kptns(:,:), hdr%occ(:), &
          &   hdr%tnons(:,:), hdr%znucltypat(:), hdr%wtk(:)

       !DEBUG
       !   write(6,*)' hdr_io : write psp record, headform= ',hdr%headform
       !ENDDEBUG

       do ipsp=1,hdr%npsp

          write(unitfi) hdr%title(ipsp), hdr%znuclpsp(ipsp), &
             &   hdr%zionpsp(ipsp), hdr%pspso(ipsp), hdr%pspdat(ipsp), &
             &   hdr%pspcod(ipsp), hdr%pspxc(ipsp), hdr%lmn_size(ipsp)

       end do

       !DEBUG
       !  write(6,*)' hdr_io : write mode, hdr%so_psp(:), hdr%symafm(:)=',&
       !&                                  hdr%so_psp(:), hdr%symafm(:)
       !ENDDEBUG

       write(unitfi) hdr%residm, hdr%xred(:,:), hdr%etot, hdr%fermie

       if (hdr%usepaw==1) then
          allocate(nsel56(hdr%natom))
          cplex=hdr%pawrhoij(1)%cplex
          nspden=hdr%pawrhoij(1)%nspden
          do iatom=1,hdr%natom
             nsel56(iatom)=hdr%pawrhoij(iatom)%nrhoijsel
          end do
          write(unitfi) (nsel56(iatom),iatom=1,hdr%natom),cplex,nspden
          bsize=sum(nsel56);allocate(ibuffer(bsize),buffer(bsize*nspden*cplex))
          ii=0;jj=0
          do iatom=1,hdr%natom
             nselect=nsel56(iatom)
             ibuffer(ii+1:ii+nselect)=hdr%pawrhoij(iatom)%rhoijselect(1:nselect)
             ii=ii+nselect
             do ispden=1,nspden
                buffer(jj+1:jj+cplex*nselect)=hdr%pawrhoij(iatom)%rhoijp(1:cplex*nselect,ispden)
                jj=jj+cplex*nselect
             end do
          end do
          write(unitfi) ibuffer(:),buffer(:)
          deallocate(ibuffer,buffer,nsel56)
       end if

       ! -------------------------------------------------------------------------
       ! Writing the header of a formatted file
       ! -------------------------------------------------------------------------

    else if(rdwr==3 .or. rdwr==4)then

       write(unitfi, '(a)' )&
          &  ' ==============================================================================='
       if(rdwr==3)write(unitfi, '(a)' ) ' ECHO of part of the ABINIT file header '
       if(rdwr==4)write(unitfi, '(a)' ) ' ECHO of the ABINIT file header '
       write(unitfi, '(a)' ) ' '
       write(unitfi, '(a)' ) ' First record :'
       write(unitfi, '(a,a6,2i5)' )  '.codvsn,headform,fform = ',&
          &  hdr%codvsn, hdr%headform, fform     ! Do not worry about 22 format

       write(unitfi, '(a)' ) ' '
       write(unitfi, '(a)' ) ' Second record :'
       write(unitfi, '(a,4i6)') ' bantot,intxc,ixc,natom  =',&
          &                            hdr%bantot, hdr%intxc, hdr%ixc, hdr%natom
       write(unitfi, '(a,4i6)') ' ngfft(1:3),nkpt         =',&
          &                            hdr%ngfft(1:3), hdr%nkpt

       if(hdr%headform>=23)then
          write(unitfi, '(a,2i6)') ' nspden,nspinor          =',&
             &                            hdr%nspden, hdr%nspinor
       end if

       if(hdr%headform<=23)then
          write(unitfi, '(a,4i6)' ) ' nsppol,nsym,ntypat,occopt=',&
             &                            hdr%nsppol,hdr%nsym,hdr%ntypat,hdr%occopt
       else if(hdr%headform<=40)then
          write(unitfi, '(a,5i6)' ) ' nsppol,nsym,npsp,ntypat,occopt=',&
             &                            hdr%nsppol,hdr%nsym,hdr%npsp,hdr%ntypat,hdr%occopt

       else if(hdr%headform==41 .or. hdr%headform==42)then
          write(unitfi, '(a,6i6)' ) ' nsppol,nsym,npsp,ntypat,occopt,pertcase=',&
             &                            hdr%nsppol,hdr%nsym,hdr%npsp,hdr%ntypat,hdr%occopt,hdr%pertcase
       else if(hdr%headform>=44)then
          write(unitfi, '(a,4i6)' ) ' nsppol,nsym,npsp,ntypat =',&
             &                            hdr%nsppol,hdr%nsym,hdr%npsp,hdr%ntypat
          write(unitfi, '(a,3i6)' ) ' occopt,pertcase,usepaw  =',&
             &                            hdr%occopt,hdr%pertcase,hdr%usepaw
       end if

       if(hdr%headform==40 .or. hdr%headform==41 .or. hdr%headform==42)then
          write(unitfi, '(a,2es18.10)') ' ecut,ecutsm             =',hdr%ecut, hdr%ecutsm
       else if(hdr%headform>=44)then
          write(unitfi, '(a,3es18.10)') ' ecut,ecutdg,ecutsm      =',&
             &   hdr%ecut, hdr%ecutdg, hdr%ecutsm
       end if

       write(unitfi, '(a, es18.10)' ) ' ecut_eff                =',hdr%ecut_eff

       if(hdr%headform>=41)then
          write(unitfi, '(a,3es18.10)') ' qptn(1:3)               =',hdr%qptn(1:3)
       end if

       write(unitfi, '(a,3es18.10)' ) ' rprimd(1:3,1)           =',hdr%rprimd(1:3,1)
       write(unitfi, '(a,3es18.10)' ) ' rprimd(1:3,2)           =',hdr%rprimd(1:3,2)
       write(unitfi, '(a,3es18.10)' ) ' rprimd(1:3,3)           =',hdr%rprimd(1:3,3)

       if(hdr%headform==40.or.hdr%headform==41)then
          write(unitfi, '(a,2es18.10)') ' tphysel,tsmear          =',hdr%tphysel, hdr%tsmear
       else if(hdr%headform>=42)then
          write(unitfi, '(a,3es18.10)') ' stmbias,tphysel,tsmear  =',&
             &   hdr%stmbias,hdr%tphysel, hdr%tsmear
       end if

       write(unitfi, '(a)' )
       if(rdwr==3)then
          write(unitfi, '(a,i3,a)' ) ' The header contain ',hdr%npsp+2,' additional records.'
       else

          write(unitfi, '(a)' ) ' Third record :'
          write(unitfi, '(a,(12i4,8x))') ' istwfk=',hdr%istwfk(:)
          write(unitfi, '(a,(12i4,8x))') ' nband =',hdr%nband(:)
          write(unitfi, '(a,(10i5,8x))') ' npwarr=',hdr%npwarr(:)

          if(hdr%headform>=40)then
             write(unitfi, '(a,(12i4,8x))') ' so_psp=',hdr%so_psp(:)
          end if

          if(hdr%headform>=40)then
             write(unitfi, '(a)') ' symafm='
             write(unitfi, '(8x,24i3,8x)') hdr%symafm(:)
          end if

          write(unitfi, '(a)' ) ' symrel='
          do isym=1,hdr%nsym/2
             write(unitfi, '(a,9i4,a,9i4)' ) '        ',hdr%symrel(:,:,2*isym-1),&
                &                                         '  ',hdr%symrel(:,:,2*isym)
          end do
          if(2*(hdr%nsym/2)/=hdr%nsym)write(unitfi, '(a,9i4)' ) '        ',hdr%symrel(:,:,hdr%nsym)

          write(unitfi, '(a,(12i4,8x))') ' type  =',hdr%typat(:)
          write(unitfi, '(a)' ) ' kptns =                 (max 50 k-points will be written)'
          do ikpt=1,min(hdr%nkpt,50)
             write(unitfi, '(a,3es16.6)' ) '        ',hdr%kptns(:,ikpt)
          end do
          write(unitfi, '(a)' ) ' wtk ='
          do ikpt=1,hdr%nkpt,10
             write(unitfi, '(a,10f6.2)' ) '        ',hdr%wtk(ikpt:min(hdr%nkpt,ikpt + 10 - 1))
          end do
          write(unitfi, '(a)' ) '   occ ='
          do ii=1,hdr%bantot,10
             write(unitfi, '(a,10f6.2)') '        ',hdr%occ(ii:min(hdr%bantot,ii+10-1))
          end do
          write(unitfi, '(a)' ) ' tnons ='
          do isym=1,hdr%nsym/2
             write(unitfi, '(a,3f10.6,a,3f10.6)' ) '        ',hdr%tnons(:,2*isym-1),&
                &                                               '  ',hdr%tnons(:,2*isym)
          end do
          if(2*(hdr%nsym/2)/=hdr%nsym)write(unitfi, '(a,3f10.6)' ) '        ',hdr%tnons(:,hdr%nsym)
          write(unitfi, '(a,(10f6.2,8x))') '  znucl=',hdr%znucltypat(:)
          write(unitfi,'(a)')

          write(unitfi, '(a)' ) ' Pseudopotential info :'
          do ipsp=1,hdr%npsp
             write(unitfi,'(a,a)' ) ' title=',trim(hdr%title(ipsp))
             write(unitfi,'(a,f6.2,a,f6.2,a,i3,a,i6,a,i3,a,i3)' ) &
                &    '  znuclpsp=',hdr%znuclpsp(ipsp),    ', zionpsp=',  hdr%zionpsp(ipsp),    &
                &    ', pspso=' , hdr%pspso(ipsp),  ', pspdat=',hdr%pspdat(ipsp),  &
                &    ', pspcod=', hdr%pspcod(ipsp), ', pspxc=', hdr%pspxc(ipsp)
             if(hdr%headform>=44)then
                if(hdr%usepaw==1)then
                   write(unitfi,'(a,i3)' ) '  lmn_size=', hdr%lmn_size(ipsp)
                else
                   write(unitfi,'(a,i3)' ) '  lmnmax  =', hdr%lmn_size(ipsp)
                end if
             end if

          end do

          write(unitfi, '(a)' ) ' '
          write(unitfi, '(a)' ) ' Last record :'
          write(unitfi, '(a,es16.6,es22.12,es16.6)' ) &
             &   ' residm,etot,fermie=',hdr%residm, hdr%etot, hdr%fermie
          write(unitfi, '(a)' ) ' xred ='
          do iatom=1,hdr%natom
             write(unitfi, '(a,3es16.6)' ) '        ',hdr%xred(:,iatom)
          end do

          if(hdr%usepaw==1)then
             natinc=1;if(hdr%natom>1) natinc=hdr%natom-1
             do iatom=1,hdr%natom,natinc
                do ispden=1,hdr%pawrhoij(iatom)%nspden
                   write(unitfi, '(a,i4,a,i1,a)' ) ' rhoij(',iatom,',',ispden,&
                      &           ')=  (max 12 non-zero components will be written)'
                   call print_ij(hdr%pawrhoij(iatom)%rhoijp(:,ispden),&
                      &                   hdr%pawrhoij(iatom)%nrhoijsel,&
                      &                   hdr%pawrhoij(iatom)%cplex,&
                      &                   hdr%pawrhoij(iatom)%lmn_size,-1,ibuffer,1,0,&
                      &                   hdr%pawrhoij(iatom)%rhoijselect(:),-1.d0,1)
                end do
             end do
          end if

          if(rdwr==3)write(unitfi, '(a)' ) ' End the ECHO of part of the ABINIT file header '
          if(rdwr==4)write(unitfi, '(a)' ) ' End the ECHO of the ABINIT file header '
          write(unitfi, '(a)' )&
             &  ' ==============================================================================='

       end if ! rdwr is 3 or 4
       call flush(unitfi)
    end if ! choice read/write/echo

    return

  end subroutine hdr_io

  !> The wrtout subroutine from abinit.
  subroutine wrtout(message)
    !Arguments ------------------------------------
    !scalars
    character(len=500),intent(inout) :: message

    !Local variables-------------------------------
    !scalars
    integer,save :: iexit=0,ncomment=0,nwarning=0
    integer :: lenmessage,rtnpos
    character(len=500) :: messtmp

    !******************************************************************
    !BEGIN EXECUTABLE SECTION

    if(message/=' ') then
       messtmp=message
       lenmessage=len(message)
       ! Here, split the message, according to the char(10)
       ! characters (carriage return). This technique is
       ! portable accross different OS.
       rtnpos=index(messtmp,ch10)
       ! do while(rtnpos/=0)
       do while(rtnpos > 1)
          write(uout, '(a)' ) trim(messtmp(1:rtnpos-1))
          messtmp=messtmp(rtnpos+1:lenmessage)
          lenmessage=lenmessage-rtnpos
          rtnpos=index(messtmp,ch10)
       end do
       write(uout, '(a)' ) trim(messtmp)
    else
       write(uout,*)
    end if

    if( index(trim(message),'BUG') /= 0 )then
       write(uout, '(a)' ) '  Action : contact ABINIT group.'
       write(uout,*)
    end if

    if( index(trim(message),'BUG') /= 0   .or. &
       & index(trim(message),'Calculation completed') /= 0 )then
       if(nwarning<10000 .and. ncomment<1000)then
          write(uout, '(a,i5,a,i4,a)' ) &
             &   '.Delivered',nwarning,' WARNINGs and',ncomment,' COMMENTs to log file.'
       else
          write(uout, '(a,i6,a,i6,a)' ) &
             &   '.Delivered',nwarning,' WARNINGs and',ncomment,' COMMENTs to log file.'
       end if
       if(iexit/=0)then
          write(uout, '(a)' ) ' Note : exit requested by the user.'
       end if
    end if

    if( index(trim(message),'Exit') /= 0 )then
       iexit=1
    end if

    !Count the number of warnings and comments. Only take into
    !account unit 6, in order not to duplicate these numbers.
    if( index(trim(message),'WARNING') /= 0 .and. uout==6 )then
       nwarning=nwarning+1
    end if
    if( index(trim(message),'COMMENT') /= 0 .and. uout==6 )then
       ncomment=ncomment+1
    end if

  end subroutine wrtout

  !> The print_ij subroutine from abinit.
  subroutine print_ij(a_ij,adim,cplex,ndim,opt_l,opt_l_index,opt_pack,opt_prtvol,pack2ij,test_value,unt, &
     &                   opt_sym,asym_ij)    !Optional arguments

    !Arguments ---------------------------------------------
    !scalars
    integer,intent(in) :: adim,cplex,ndim,opt_l,opt_pack,opt_prtvol,unt
    integer,intent(in),optional :: opt_sym
    real(dp),intent(in) :: test_value
    !arrays
    integer,intent(in) :: opt_l_index(ndim*min(1+opt_l,1)),pack2ij(adim*opt_pack)
    real(dp),intent(in) :: a_ij(cplex*adim)
    real(dp),intent(in),optional :: asym_ij(cplex*adim)

    !Local variables ---------------------------------------
    ! Adjust format bellow according to maxprt
    !scalars
    integer,parameter :: maxprt_default=12
    integer :: dplex,ilmn,ilmn1,j0lmn,jlmn,jlmn1,klmn,klmn1,klmn2,maxprt,nhigh
    integer :: nmin,optsym
    real(dp) :: testval
    logical :: use_asym
    character(len=500) :: message
    !arrays
    real(dp),parameter :: fact_re(4)=(/1d0,1d0,-1d0,-1d0/),fact_im(4)=(/1d0,-1d0,-1d0,1d0/)
    real(dp) :: tabmax(cplex),tabmin(cplex)
    real(dp),allocatable :: b_ij(:),bsym_ij(:),prtab(:,:,:)

    ! *************************************************************************

10  format(100(1x,f9.5))
11  format(12(1x,f9.5),a) !Change this format according to variable "maxprt"

    !DEBUG
    !write(6,*)' print_ij : enter '
    !ENDDEBUG

    !Optional arguments
    use_asym=present(asym_ij)
    if (present(opt_sym)) then
       optsym=opt_sym
    else
       optsym=2
    end if

    !Define size of square matrix
    if (opt_prtvol>=0) then
       maxprt=maxprt_default
    else
       maxprt=ndim
    end if
    nmin=min(ndim,maxprt)

    if (opt_l>=0) nmin=count(opt_l_index(:)==opt_l)
    allocate(prtab(cplex,nmin,nmin))
    dplex=cplex-1

    !Eventually unpack input matrix(es)
    allocate(b_ij(cplex*ndim*(ndim+1)/2))
    if (opt_pack==0) then
       b_ij=a_ij
    else if (opt_pack==1) then
       b_ij=0d0
       do klmn=1,adim
          klmn1=cplex*klmn-dplex
          klmn2=cplex*pack2ij(klmn)-dplex
          b_ij(klmn2:klmn2+dplex)=a_ij(klmn1:klmn1+dplex)
       end do
    end if
    if (opt_prtvol<0.and.opt_l<0) then
       if (cplex==1) then
          tabmax(1)=maxval(abs(b_ij))
          tabmin(1)=minval(abs(b_ij))
       else
          tabmax(1:2)=0d0;tabmin(1:2)=1.d20
          do klmn=1,ndim
             klmn2=2*klmn
             tabmax(1)=max(tabmax(1),b_ij(klmn2-1))
             tabmin(1)=min(tabmin(1),b_ij(klmn2-1))
             tabmax(2)=max(tabmax(2),b_ij(klmn2  ))
             tabmin(2)=min(tabmin(2),b_ij(klmn2  ))
          end do
       end if
    end if
    if (use_asym) then
       allocate(bsym_ij(cplex*ndim*(ndim+1)/2))
       if (opt_pack==0) then
          bsym_ij=asym_ij
       else if (opt_pack==1) then
          bsym_ij=0d0
          do klmn=1,adim
             klmn1=cplex*klmn-dplex
             klmn2=cplex*pack2ij(klmn)-dplex
             bsym_ij(klmn2:klmn2+dplex)=asym_ij(klmn1:klmn1+dplex)
          end do
       end if
       if (opt_prtvol<0.and.opt_l<0) then
          if (cplex==1) then
             tabmax(1)=max(tabmax(1),maxval(abs(bsym_ij)))
             tabmin(1)=min(tabmin(1),minval(abs(bsym_ij)))
          else
             do klmn=1,ndim
                klmn2=2*klmn
                tabmax(1)=max(tabmax(1),bsym_ij(klmn2-1))
                tabmin(1)=min(tabmin(1),bsym_ij(klmn2-1))
                tabmax(2)=max(tabmax(2),bsym_ij(klmn2  ))
                tabmin(2)=min(tabmin(2),bsym_ij(klmn2  ))
             end do
          end if
       end if
    end if

    !Transfer triangular matrix to rectangular one
    jlmn1=0
    do jlmn=1,ndim
       if (opt_l<0) then
          jlmn1=jlmn;if (jlmn1>nmin) cycle
       else if (opt_l_index(jlmn)==opt_l) then
          jlmn1=jlmn1+1
       else
          cycle
       end if
       ilmn1=0;j0lmn=jlmn*(jlmn-1)/2
       do ilmn=1,jlmn
          if (opt_l<0) then
             ilmn1=ilmn
          else if (opt_l_index(ilmn)==opt_l) then
             ilmn1=ilmn1+1
          else
             cycle
          end if
          klmn=j0lmn+ilmn
          if (cplex==1) then
             prtab(1,ilmn1,jlmn1)=b_ij(klmn)
             if (use_asym) then
                prtab(1,jlmn1,ilmn1)=fact_re(optsym)*bsym_ij(klmn)
             else
                prtab(1,jlmn1,ilmn1)=fact_re(optsym)*b_ij(klmn)
             end if
          else
             klmn=2*klmn
             prtab(1:2,ilmn1,jlmn1)=b_ij(klmn-1:klmn)
             if (use_asym) then
                prtab(1,jlmn1,ilmn1)=fact_re(optsym)*bsym_ij(klmn-1)
                prtab(2,jlmn1,ilmn1)=fact_im(optsym)*bsym_ij(klmn  )
             else
                prtab(1,jlmn1,ilmn1)=fact_re(optsym)*b_ij(klmn-1)
                prtab(2,jlmn1,ilmn1)=fact_im(optsym)*b_ij(klmn  )
             end if
          end if
       end do
    end do
    deallocate(b_ij)

    if (use_asym) deallocate(bsym_ij)

    if (unt==2) then
       prtab=prtab*Ha_eV
       if (opt_prtvol<0.and.opt_l<0) then
          tabmax=tabmax*Ha_eV
          tabmin=tabmin*Ha_eV
       end if
    end if

    if (cplex==2) then
       write(message,'(3x,a)') '=== REAL PART:'
       call wrtout(message)
    end if

    if (ndim<=maxprt.or.opt_l>=0) then
       do ilmn=1,nmin
          write(message,fmt=10) prtab(1,1:nmin,ilmn)
          call wrtout(message)
       end do
    else
       do ilmn=1,nmin
          write(message,fmt=11) prtab(1,1:nmin,ilmn),' ...'
          call wrtout(message)
       end do
       write(message,'(3x,a,i2,a)') '...  only ',maxprt,'  components have been written...'
       call wrtout(message)
    end if
    if (opt_prtvol<0.and.opt_l<0) then
       write(message,'(3x,2(a,es9.2))') 'max. value= ',tabmax(1),', min. value= ',tabmin(1)
       call wrtout(message)
    end if

    if (cplex==2) then
       write(message,'(3x,a)') '=== IMAGINARY PART:'
       call wrtout(message)
       if (ndim<=maxprt.or.opt_l>=0) then
          do ilmn=1,nmin
             write(message,fmt=10) prtab(2,1:nmin,ilmn)
             call wrtout(message)
          end do
       else
          do ilmn=1,nmin
             write(message,fmt=11) prtab(2,1:nmin,ilmn),' ...'
             call wrtout(message)
          end do
          write(message,'(3x,a,i2,a)') '...  only ',maxprt,'  components have been written...'
          call wrtout(message)
       end if
       if (opt_prtvol<0.and.opt_l<0) then
          write(message,'(3x,2(a,es9.2))') 'max. value= ',tabmax(2),', min. value= ',tabmin(2)
          call wrtout(message)
       end if
    end if

    if (test_value>0d0) then
       testval=test_value;if (unt==2) testval=testval*Ha_eV
       nhigh=0;nhigh=count(abs(prtab(:,:,:))>=testval)
       if (nhigh>0) then
          write(message,'(5a,i3,a,f6.1,7a)')&
             &   ' print_ij: WARNING -',ch10,&
             &   '  The matrix seems to have high value(s) !',ch10,&
             &   '  (',nhigh,' components have a value greater than ',testval,').',ch10,&
             &   '  It can cause instabilities during SCF convergence.',ch10,&
             &   '  Action: you should check your atomic dataset (psp file)',ch10,&
             &   '          and look for "high" projector functions...'
          call wrtout(message)
       end if
    end if

    deallocate(prtab)

    !DEBUG
    !write(6,*)' print_ij : exit '
    !ENDDEBUG

  end subroutine print_ij

  !> The rhoij_alloc subroutine from abinit.
  subroutine rhoij_alloc(cplex,nlmn,nspden,nsppol,pawrhoij,typat,&      ! Mandatory arguments
     &                      ngrhoij,nlmnmix,use_rhoij_,use_rhoijres) ! Optional arguments
    !! NAME
    !! rhoij_alloc
    !!
    !! FUNCTION
    !! Initialize and allocate a pawrhoij datastructure
    !!
    !! COPYRIGHT
    !! Copyright (C) 2007-2009 ABINIT group (MT)
    !! This file is distributed under the terms of the
    !! GNU General Public License, see ~ABINIT/Infos/copyright
    !! or http://www.gnu.org/copyleft/gpl.txt .
    !!
    !! INPUTS
    !! cplex=1 if rhoij is REAL,2 if COMPLEX
    !! nlmn(:)=array of (l,m,n) sizes for rhoij
    !! nspden=number of spin-components for rhoij
    !! nsppol=number of independant spin-components for rhoij
    !! typat(:)=types of atoms
    !! ngrhoij=number of gradients to be allocated (OPTIONAL, default=0)
    !! nlmnmix=number of rhoij elements to be mixed during SCF cycle (OPTIONAL, default=0)
    !! use_rhoij_=1 if pawrhoij(:)%rhoij_ has to be allocated (OPTIONAL, default=0)
    !! use_rhoijres=1 if pawrhoij(:)%rhoijres has to be allocated (OPTIONAL, default=0)
    !!
    !! SIDE EFFECTS
    !! pawrhoij(:)<type(pawrhoij_type)>= rhoij datastructure

    !Arguments ------------------------------------
    !scalars
    integer,intent(in) :: cplex,nspden,nsppol
    integer,intent(in),optional :: ngrhoij,nlmnmix,use_rhoij_,use_rhoijres
    !arrays
    integer,intent(in) :: nlmn(:),typat(:)
    type(pawrhoij_type),intent(inout) :: pawrhoij(:)

    !Local variables-------------------------------
    !scalars
    integer :: irhoij,itypat,lmn2_size,nn1,nn2,nrhoij

    ! *************************************************************************

    nrhoij=size(pawrhoij);nn1=size(nlmn);nn2=size(typat)
    if (nrhoij/=nn2.or.maxval(typat)>nn1) stop "Error in rhoij_alloc: wrong sizes ! "

    do irhoij=1,nrhoij

       itypat=typat(irhoij)
       lmn2_size=nlmn(itypat)*(nlmn(itypat)+1)/2

       ! Scalars initializations
       pawrhoij(irhoij)%cplex=cplex
       pawrhoij(irhoij)%lmn_size=nlmn(itypat)
       pawrhoij(irhoij)%lmn2_size=lmn2_size
       pawrhoij(irhoij)%nspden=nspden
       pawrhoij(irhoij)%nsppol=nsppol
       pawrhoij(irhoij)%nrhoijsel=0
       pawrhoij(irhoij)%lmnmix_sz=0
       pawrhoij(irhoij)%ngrhoij=0
       pawrhoij(irhoij)%use_rhoij_=0
       pawrhoij(irhoij)%use_rhoijres=0

       ! Mandatory pointers allocations
       allocate(pawrhoij(irhoij)%rhoijselect(lmn2_size))
       allocate(pawrhoij(irhoij)%rhoijp(cplex*lmn2_size,nspden))

       ! Optional pointers allocations
       if (present(ngrhoij)) then
          if (ngrhoij>0) then
             pawrhoij(irhoij)%ngrhoij=ngrhoij
             allocate(pawrhoij(irhoij)%grhoij(ngrhoij,cplex*lmn2_size,nspden))
          end if
       end if
       if (present(nlmnmix)) then
          if (nlmnmix>0) then
             pawrhoij(irhoij)%lmnmix_sz=nlmnmix
             allocate(pawrhoij(irhoij)%kpawmix(nlmnmix))
          end if
       end if
       if (present(use_rhoij_)) then
          if (use_rhoij_>0) then
             pawrhoij(irhoij)%use_rhoij_=use_rhoij_
             allocate(pawrhoij(irhoij)%rhoij_(cplex*lmn2_size,nspden))
          end if
       end if
       if (present(use_rhoijres)) then
          if (use_rhoijres>0) then
             pawrhoij(irhoij)%use_rhoijres=use_rhoijres
             allocate(pawrhoij(irhoij)%rhoijres(cplex*lmn2_size,nspden))
          end if
       end if

       ! Intializations to zero (to avoid overflow)
       pawrhoij(irhoij)%rhoijselect(:)=0
       pawrhoij(irhoij)%rhoijp(:,:)=0d0

    end do

  end subroutine rhoij_alloc

end module abinit_private
