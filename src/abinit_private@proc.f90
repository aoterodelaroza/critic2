!> Abinit structure and scalar field reader.
!> Most routines in this module were adapted from abinit.
! COPYRIGHT
! Copyright (C) 1998-2017 ABINIT group (DCA, XG, GMR, MKV, MT)
! This file is distributed under the terms of the
! GNU General Public License, see ~abinit/COPYING
! or http://www.gnu.org/copyleft/gpl.txt .
! For the initials of contributors, see ~abinit/doc/developers/contributor

! Copyright (c) 2007-2018 Alberto Otero de la Roza <aoterodelaroza@gmail.com>,
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

submodule (abinit_private) proc
  implicit none

contains
  
  ! Driver to choose which version of hdr_io to use.
  module subroutine hdr_io(fform,hdr,rdwr,unitfi)
    use tools_io, only: ferror, faterr, uout
    integer, intent(inout) :: fform
    integer, intent(in) :: rdwr,unitfi
    type(hdr_type), intent(inout) :: hdr

    type(hdr_type_1) :: hdr_1
    type(hdr_type_2) :: hdr_2

    ! read the headform
    if (rdwr /= 1 .and. rdwr /= 5) then
       call ferror("hdr_io","only read operations allowed in abinit files",faterr)
    end if
    if (rdwr == 1) then
       rewind(unitfi)
    end if

    call hdr_io_1(fform,hdr_1,unitfi)
    if (headform_1 /= 0) then
       rewind(unitfi)
       call hdr_io_2(fform,hdr_2,unitfi)
       if (headform_2 /= 0) then
          write (uout,'("Can not handle headform: ",I3)') headform_2
          write (uout,'("The version of abinit you are using is not supported by critic2.")') 
          write (uout,'("Please, e-mail the critic2 developer this message and the version")') 
          write (uout,'("of abinit you are using.")') 
          call ferror("hdr_io","unknown abinit headform",faterr)
       end if
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

  ! private
  
  ! The following functions are from ABINIT
  ! Copyright (C) 2007-2009 ABINIT group (MT)
  ! This file is distributed under the terms of the
  ! GNU General Public License, see ~ABINIT/Infos/copyright
  ! or http://www.gnu.org/copyleft/gpl.txt .
  
  !> The hdr_io subroutine from abinit.
  module subroutine hdr_io_1(fform,hdr,unitfi)
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
    integer,intent(in) :: unitfi
    type(hdr_type_1),intent(inout) :: hdr

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
    
    headform_1 = 0
    ! Reading the first record of the file ------------------------------------

    read(unitfi,iostat=ierr)codvsn,fform
    if (ierr /=0) then
       fform=0
       headform = -1
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

    if(headform==22)then

       read(unitfi) bantot, hdr%date, hdr%intxc, hdr%ixc, natom, hdr%ngfft(1:3),&
          nkpt, nsppol, nsym, ntypat,&
          acell, hdr%ecut_eff, hdr%rprimd
       npsp=ntypat

    else if(headform==23)then

       !   Compared to v2.2, add nspden, nspinor, occopt
       read(unitfi) bantot, hdr%date, hdr%intxc, hdr%ixc, natom, hdr%ngfft(1:3),&
          nkpt, hdr%nspden, hdr%nspinor, nsppol, nsym, ntypat, hdr%occopt,&
          acell, hdr%ecut_eff, hdr%rprimd
       npsp=ntypat

    else if(headform==34)then

       !   Compared to v2.3, subtract acell, and add npsp
       read(unitfi) bantot, hdr%date, hdr%intxc, hdr%ixc, natom, hdr%ngfft(1:3),&
          nkpt, hdr%nspden, hdr%nspinor, nsppol, nsym, npsp, ntypat, hdr%occopt,&
          hdr%ecut_eff, hdr%rprimd

    else if(headform==40)then

       !   Compared to v3.4, add ecut, ecutsm, tphysel, tsmear
       read(unitfi) bantot, hdr%date, hdr%intxc, hdr%ixc, natom, hdr%ngfft(1:3),&
          nkpt, hdr%nspden, hdr%nspinor, nsppol, nsym, npsp, ntypat, hdr%occopt,&
          hdr%ecut, hdr%ecutsm, hdr%ecut_eff, hdr%rprimd, hdr%tphysel, hdr%tsmear

    else if(headform==41)then

       !   Compared to v4.0, add pertcase and qptn(3)
       read(unitfi) bantot, hdr%date, hdr%intxc, hdr%ixc, natom, hdr%ngfft(1:3),&
          nkpt, hdr%nspden, hdr%nspinor, nsppol, nsym, npsp, ntypat, hdr%occopt, hdr%pertcase,&
          hdr%ecut, hdr%ecutsm, hdr%ecut_eff, hdr%qptn(1:3), hdr%rprimd, hdr%tphysel, hdr%tsmear

    else if(headform==42)then

       !   Compared to v4.1, add stmbias
       read(unitfi) bantot, hdr%date, hdr%intxc, hdr%ixc, natom, hdr%ngfft(1:3),&
          nkpt, hdr%nspden, hdr%nspinor, nsppol, nsym, npsp, ntypat, hdr%occopt, hdr%pertcase,&
          hdr%ecut, hdr%ecutsm, hdr%ecut_eff, hdr%qptn(1:3), hdr%rprimd,&
          hdr%stmbias, hdr%tphysel, hdr%tsmear

    else if(headform>=44 .and. headform<57)then

       !   Compared to v4.2, add usepaw and ecutdg
       read(unitfi) bantot, hdr%date, hdr%intxc, hdr%ixc, natom, hdr%ngfft(1:3),&
          nkpt, hdr%nspden, hdr%nspinor, nsppol, nsym, npsp, ntypat, hdr%occopt, hdr%pertcase,&
          hdr%usepaw, hdr%ecut, hdr%ecutdg, hdr%ecutsm, hdr%ecut_eff, hdr%qptn(1:3), hdr%rprimd,&
          hdr%stmbias, hdr%tphysel, hdr%tsmear

    else if(headform>=57)then

       !   Compared to v4.4, add usewvl
       read(unitfi) bantot, hdr%date, hdr%intxc, hdr%ixc, natom, hdr%ngfft(1:3),&
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
       write(message, '(4a,es16.6,9a)' ) ch10,&
          ' hdr_io : ERROR -',ch10,&
          '  The value of ecutsm is',hdr%ecutsm,', while the file has been produced prior to v4.4 .',ch10,&
          '  The definition of the smearing function has changed, so that you are not allowed',ch10,&
          '  to restart from a old wavefunction file. By contrast, you can restart from an old',ch10,&
          '  potential or density file, and perform a self-consistent cycle with a new ABINIT version.',ch10,&
          '  Action : produce a density or potential file using the old version of ABINIT, and restart from it.'
       write (*,*) message
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
          hdr%typat(:), hdr%kptns(:,:), hdr%occ(:), &
          hdr%tnons(:,:), hdr%znucltypat(:)

    else if(headform==22 .or. headform==23 .or. headform==34)then

       !   Compared to pre v2.0, add istwfk
       read(unitfi) hdr%nband(:), hdr%npwarr(:), hdr%symrel(:,:,:), &
          hdr%typat(:), hdr%istwfk(:), hdr%kptns(:,:), hdr%occ(:), &
          hdr%tnons(:,:), hdr%znucltypat(:)

    else if(headform>=40 .and. headform < 50)then

       !   Compared to pre v4.0, add so_psp and symafm, and switch istwfk

       read(unitfi)  hdr%istwfk(:), hdr%nband(:), hdr%npwarr(:), &
          hdr%so_psp(:), hdr%symafm(:), hdr%symrel(:,:,:), &
          hdr%typat(:), hdr%kptns(:,:), hdr%occ(:), &
          hdr%tnons(:,:), hdr%znucltypat(:)

    else if(headform>=50)then

       !   Compared to pre v5.0, add wtk
       read(unitfi)  hdr%istwfk(:), hdr%nband(:), hdr%npwarr(:), &
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
          read(unitfi) hdr%title(ipsp), hdr%znuclpsp(ipsp), &
             hdr%zionpsp(ipsp), hdr%pspdat(ipsp), hdr%pspcod(ipsp), &
             hdr%pspxc(ipsp), lmax, lloc, mmax
       end do

    else if(headform==23)then

       !  Compared to 2.2, add pspso
       do ipsp=1,npsp
          read(unitfi) hdr%title(ipsp), hdr%znuclpsp(ipsp), &
             hdr%zionpsp(ipsp), hdr%pspso(ipsp), hdr%pspdat(ipsp), &
             hdr%pspcod(ipsp), hdr%pspxc(ipsp), lmax, lloc, mmax
       end do

    else if(headform==34 .or. headform==40 .or. headform==41 &
       .or. headform==42)then

       !  Compared to 2.3, suppress lmax, lloc, mmax
       do ipsp=1,npsp
          read(unitfi) hdr%title(ipsp), hdr%znuclpsp(ipsp), &
             hdr%zionpsp(ipsp), hdr%pspso(ipsp), hdr%pspdat(ipsp), &
             hdr%pspcod(ipsp), hdr%pspxc(ipsp)
       end do

    else if(headform>=44)then

       !  Compared to 4.2, add lmn_size
       do ipsp=1,npsp
          read(unitfi) hdr%title(ipsp), hdr%znuclpsp(ipsp), &
             hdr%zionpsp(ipsp), hdr%pspso(ipsp), hdr%pspdat(ipsp), &
             hdr%pspcod(ipsp), hdr%pspxc(ipsp), hdr%lmn_size(ipsp)
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

  end subroutine hdr_io_1

  !> The rhoij_alloc subroutine from abinit.
  module subroutine rhoij_alloc(cplex,nlmn,nspden,nsppol,pawrhoij,typat)
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
    !!
    !! SIDE EFFECTS
    !! pawrhoij(:)<type(pawrhoij_type)>= rhoij datastructure

    !Arguments ------------------------------------
    !scalars
    integer,intent(in) :: cplex,nspden,nsppol
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

       ! Mandatory pointers allocations
       allocate(pawrhoij(irhoij)%rhoijselect(lmn2_size))
       allocate(pawrhoij(irhoij)%rhoijp(cplex*lmn2_size,nspden))

       ! Intializations to zero (to avoid overflow)
       pawrhoij(irhoij)%rhoijselect(:)=0
       pawrhoij(irhoij)%rhoijp(:,:)=0d0
    end do

  end subroutine rhoij_alloc

  module subroutine hdr_io_2(fform,hdr,unit)
    use tools_io, only: ferror, faterr
    implicit none
    integer, intent(out) :: fform
    integer, intent(in) :: unit
    type(hdr_type_2), intent(out) :: hdr

    !Local variables-------------------------------
    integer :: ipsp
    character(len=500) :: msg,errmsg
    real(dp),allocatable :: occ3d(:,:,:)
    integer :: ii,band,ikpt,spin

    !*************************************************************************

    ! Reading the first record of the file ------------------------------------
    ! fform is not a record of hdr_type
    read(unit, err=10, iomsg=errmsg) hdr%codvsn,hdr%headform,fform

    if (hdr%headform < 80) then
       headform_2 = hdr%headform
       return
    end if

    !Reading the second record of the file ------------------------------------
    read(unit, err=10, iomsg=errmsg) &
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

    read(unit, err=10, iomsg=errmsg) &
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
    read(unit, err=10, iomsg=errmsg) hdr%residm, hdr%xred(:,:), hdr%etot, hdr%fermie, hdr%amu(:)

    read(unit, err=10, iomsg=errmsg)&
       hdr%kptopt,hdr%pawcpxocc,hdr%nelect,hdr%charge,hdr%icoulomb,&
       hdr%kptrlatt,hdr%kptrlatt_orig, hdr%shiftk_orig,hdr%shiftk

    ! Reading the records with psp information ---------------------------------
    do ipsp=1,hdr%npsp
       read(unit, err=10, iomsg=errmsg) &
          hdr%title(ipsp), hdr%znuclpsp(ipsp), hdr%zionpsp(ipsp), hdr%pspso(ipsp), hdr%pspdat(ipsp), &
          hdr%pspcod(ipsp), hdr%pspxc(ipsp), hdr%lmn_size(ipsp), hdr%md5_pseudos(ipsp)
    end do

    if (hdr%usepaw==1) then ! Reading the Rhoij tab if the PAW method was used.
       call pawrhoij_io(hdr%pawrhoij,unit,hdr%nsppol,hdr%nspinor,hdr%nspden,hdr%lmn_size,hdr%typat,hdr%headform)
    end if

    return

    ! Handle IO-error: write warning and let the caller handle the exception.
10  continue
    call ferror("hdr_io_2","error reading abinit file",faterr)

  end subroutine hdr_io_2

  subroutine pawrhoij_io(pawrhoij,unitfi,nsppol_in,nspinor_in,nspden_in,nlmn_type,typat,headform)
    implicit none
    !Arguments ------------------------------------
    !scalars
    integer,intent(in) :: unitfi,headform,nspden_in,nspinor_in,nsppol_in
    !arrays
    integer,intent(in) :: typat(:),nlmn_type(:)
    type(pawrhoij_type),intent(inout) :: pawrhoij(:)

    !Local variables-------------------------------
    !scalars
    integer,parameter :: fort_formatted=1,fort_binary=2,netcdf_io=3
    integer :: cplex,i1,i2,iatom,iatom_tot,natom,ispden,bsize,ii,jj,lmn2_size
    integer :: nselect,my_cplex,my_natinc,my_natom,my_nspden,ngrhoijmx,size_rhoij2
    integer :: iomode,ncid,natom_id,cplex_id,nspden_id,nsel56_id,buffer_id,ibuffer_id,ncerr
    integer :: bsize_id,bufsize_id
    logical :: paral_atom
    character(len=500) :: msg
    !arrays
    integer,allocatable :: ibuffer(:),nsel44(:,:),nsel56(:)
    real(dp), allocatable :: buffer(:)

    ! *************************************************************************

    my_natom=SIZE(pawrhoij);if (my_natom==0) return
    my_nspden=nspden_in
    natom=size(typat)
    paral_atom=(my_natom/=natom)

    iomode = fort_binary

    ncid = unitfi

    if ((headform>=44).and.(headform<56)) then
       allocate(nsel44(nspden_in,natom))
       read(unitfi  ) ((nsel44(ispden,iatom),ispden=1,nspden_in),iatom=1,natom)
       call pawrhoij_alloc(pawrhoij,1,nspden_in,nspinor_in,nsppol_in,typat,nlmn_type)
       do iatom=1,natom
          pawrhoij(iatom)%nrhoijsel=nsel44(1,iatom)
       end do
       bsize=sum(nsel44)
       allocate(ibuffer(bsize))
       allocate(buffer(bsize))
       read(unitfi  ) ibuffer(:),buffer(:)
       ii=0
       do iatom=1,natom
          nselect=nsel44(1,iatom)
          pawrhoij(iatom)%rhoijselect(1:nselect)=ibuffer(ii+1:ii+nselect)
          do ispden=1,nspden_in
             pawrhoij(iatom)%rhoijp(1:nselect,ispden)=buffer(ii+1:ii+nselect)
             ii=ii+nselect
          end do
       end do
       deallocate(ibuffer)
       deallocate(buffer)
       deallocate(nsel44)
    else if (headform>=56) then
       allocate(nsel56(natom))
       if (headform==56) then
          read(unitfi  ) (nsel56(iatom),iatom=1,natom),my_cplex
       else
          read(unitfi  ) (nsel56(iatom),iatom=1,natom),my_cplex,my_nspden
       end if
       call pawrhoij_alloc(pawrhoij,my_cplex,my_nspden,nspinor_in,nsppol_in,typat,nlmn_type)
       do iatom=1,natom
          pawrhoij(iatom)%nrhoijsel=nsel56(iatom)
       end do
       bsize=sum(nsel56)
       allocate(ibuffer(bsize))
       allocate(buffer(bsize*my_nspden*my_cplex))
       read(unitfi  ) ibuffer(:),buffer(:)
       ii=0;jj=0
       do iatom=1,natom
          nselect=nsel56(iatom)
          pawrhoij(iatom)%rhoijselect(1:nselect)=ibuffer(ii+1:ii+nselect)
          ii=ii+nselect
          do ispden=1,my_nspden
             pawrhoij(iatom)%rhoijp(1:my_cplex*nselect,ispden)=buffer(jj+1:jj+my_cplex*nselect)
             jj=jj+my_cplex*nselect
          end do
       end do
       deallocate(ibuffer)
       deallocate(buffer)
       deallocate(nsel56)
    end if

  end subroutine pawrhoij_io

  module subroutine pawrhoij_alloc(pawrhoij,cplex,nspden,nspinor,nsppol,typat,lmnsize)
    implicit none

    !Arguments ------------------------------------
    !scalars
    integer,intent(in) :: cplex,nspden,nspinor,nsppol
    !arrays
    integer,intent(in) :: typat(:)
    integer,target,intent(in) :: lmnsize(:)
    type(pawrhoij_type),intent(inout) :: pawrhoij(:)

    !Local variables-------------------------------
    !scalars
    integer :: irhoij,irhoij_,itypat,lmn2_size,my_comm_atom,nn1,natom,nrhoij
    logical :: has_rhoijp
    character(len=500) :: msg

    nrhoij=size(pawrhoij)
    natom=size(typat)
    if (nrhoij>natom) then
       msg=' wrong sizes (1) !'
       write (*,*) msg
       stop 1
    end if

    !Select lmn_size for each atom type
    nn1=size(lmnsize)
    if (maxval(typat)>nn1) then
       msg=' wrong sizes (3) !'
       write (*,*) msg
       stop 1
    end if

    if (nrhoij>0) then
       do irhoij=1,nrhoij
          irhoij_=irhoij
          itypat=typat(irhoij_)

          lmn2_size=lmnsize(itypat)*(lmnsize(itypat)+1)/2

          !    Scalars initializations
          pawrhoij(irhoij)%cplex=cplex
          pawrhoij(irhoij)%itypat=itypat
          pawrhoij(irhoij)%lmn_size=lmnsize(itypat)
          pawrhoij(irhoij)%lmn2_size=lmn2_size
          pawrhoij(irhoij)%nspden=nspden
          pawrhoij(irhoij)%nspinor=nspinor
          pawrhoij(irhoij)%nsppol=nsppol
          pawrhoij(irhoij)%nrhoijsel=0
          pawrhoij(irhoij)%lmnmix_sz=0
          pawrhoij(irhoij)%ngrhoij=0
          pawrhoij(irhoij)%use_rhoij_=0
          pawrhoij(irhoij)%use_rhoijres=0

          !    Arrays allocations
          has_rhoijp=.true.
          if (has_rhoijp) then
             pawrhoij(irhoij)%use_rhoijp=1
             allocate(pawrhoij(irhoij)%rhoijselect(lmn2_size))
             allocate(pawrhoij(irhoij)%rhoijp(cplex*lmn2_size,nspden))
             pawrhoij(irhoij)%rhoijselect(:)=0
             pawrhoij(irhoij)%rhoijp(:,:) = 0
          end if
       end do
    end if

  end subroutine pawrhoij_alloc

end submodule proc
