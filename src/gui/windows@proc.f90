! Copyright (c) 2019-2022 Alberto Otero de la Roza <aoterodelaroza@gmail.com>,
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

! The class to handle ImGui windows, general routines.
submodule (windows) proc
  use interfaces_cimgui
  implicit none

  ! initial side for the view texture
  integer(c_int), parameter :: initial_texture_side = 1024_c_int

  ! strings and permutations for the open file dialogs
  integer, allocatable :: isperm_openfiles(:), isperm_inv_openfiles(:)
  character(kind=c_char,len=:), allocatable, target :: combostr_openfiles, dialogstr_openfiles
  integer, allocatable :: isperm_openfieldfile(:), isperm_inv_openfieldfile(:)
  character(kind=c_char,len=:), allocatable, target :: combostr_openfieldfile, dialogstr_openfieldfile
  integer, allocatable :: isperm_openvibfile(:), isperm_inv_openvibfile(:)
  character(kind=c_char,len=:), allocatable, target :: combostr_openvibfile, dialogstr_openvibfile

  !xx! private procedures
  ! subroutine dialog_user_callback(vFilter, vUserData, vCantContinue)

contains

  !> Initialize data required in the windows@proc submodule
  module subroutine windows_init()
    use param

    integer :: i, nn

    ! open structure file
    !! for the dialog: list of extensions "seen"
    dialogstr_openfiles = &
       "&
       &All files (*.*){*.*},&
       &ABINIT (DEN...){(DEN|ELF|POT|VHA\VHXC|VXC|GDEN1|GDEN2|GDEN3|LDEN|KDEN|PAWDEN|VCLMB|VPSP)},&
       &ADF (molden){.molden},&
       &cif (cif){.cif},&
       &CASTEP (cell|geom){.cell,.geom},&
       &CRYSTAL (out){.out},&
       &cube (cube|bincube){.cube,.bincube},&
       &DFTB+ (gen){.gen},&
       &DMACRYS (dmain|16|21){.dmain,.16,.21},&
       &elk (OUT){.OUT},&
       &FHIaims (in|in.next_step|out|own){.in,.next_step,.out,.own},&
       &FPLO (out){.out},&
       &Gaussian (com|gjf|zmat|log|wfn|wfx|fchk|cube){.com,.gjf,.zmat,.log,.wfn,.wfx,.fchk,.cube},&
       &ORCA (out|molden|molden.input){.out,.molden,.input},&
       &pdb (pdb){.pdb},&
       &postg (pgout){.pgout},&
       &psi4 (molden|dat){.molden,.dat},&
       &Quantum ESPRESSO (out|in|pwi|pwo|cube|pwc) {.out,.in,.pwi,.pwo,.cube,.pwc},&
       &SHELX (res|ins){.res,.ins.16},&
       &SIESTA (STRUCT_IN|STRUCT_OUT) {.STRUCT_IN,.STRUCT_OUT},&
       &TINKER (frac) {.frac},&
       &Tripos MOL2 (mol2){.mol2},&
       &VASP (POSCAR|CONTCAR|...){(CONTCAR|CHGCAR|ELFCAR|CHG|AECCAR0|AECCAR1|AECCAR2|POSCAR)},&
       &WIEN2k (struct){.struct},&
       &Xcrysden (xsf|axsf) {.xsf,.axsf},&
       &xyz (xyz){.xyz},&
       &"// c_null_char

    !! for the dialog: options to force a particular format
    combostr_openfiles = "" &
       // "Auto-detect" // c_null_char &             ! isformat_unknown
       // "Abinit DEN-style file" // c_null_char &   ! isformat_abinit
       // "Binary cube" // c_null_char &             ! isformat_bincube
       // "CASTEP cell" // c_null_char &             ! isformat_castepcell
       // "CASTEP geom" // c_null_char &             ! isformat_castepgeom
       // "CIF file" // c_null_char &                ! isformat_cif
       // "CRYSTAL output" // c_null_char &          ! isformat_crystal
       // "Cube" // c_null_char &                    ! isformat_cube
       // "DFTB+ gen file" // c_null_char &          ! isformat_gen
       // "DMACRYS .21 file" // c_null_char &        ! isformat_f21
       // "DMACRYS dmain file" // c_null_char &      ! isformat_dmain
       // "elk GEOMETRY.OUT" // c_null_char &        ! isformat_elk
       // "FHIaims input" // c_null_char &           ! isformat_aimsin
       // "FHIaims output" // c_null_char &          ! isformat_aimsout
       // "FPLO output" // c_null_char &             ! isformat_fploout
       // "Gaussian fchk" // c_null_char &           ! isformat_fchk
       // "Gaussian input" // c_null_char &          ! isformat_gjf
       // "Gaussian output" // c_null_char &         ! isformat_gaussian
       // "Gaussian wfn" // c_null_char &            ! isformat_wfn
       // "Gaussian wfx" // c_null_char &            ! isformat_wfx
       // "Gaussian zmat file" // c_null_char &      ! isformat_zmat
       // "mol2 file" // c_null_char &               ! isformat_mol2
       // "Molden-style file" // c_null_char &       ! isformat_molden
       // "ORCA output file" // c_null_char &        ! isformat_orca
       // "pdb file" // c_null_char &                ! isformat_pdb
       // "postg output" // c_null_char &            ! isformat_pgout
       // "psi4 output" // c_null_char &             ! isformat_dat
       // "Quantum ESPRESSO input" // c_null_char &  ! isformat_qein
       // "Quantum ESPRESSO output" // c_null_char & ! isformat_qeout
       // "Quantum ESPRESSO pwc" // c_null_char &    ! isformat_pwc
       // "SHELX" // c_null_char &                   ! isformat_shelx
       // "SIESTA IN/OUT file" // c_null_char &      ! isformat_siesta
       // "TINKER frac file" // c_null_char &        ! isformat_tinkerfrac
       // "VASP" // c_null_char &                    ! isformat_vasp
       // "WIEN2k struct file" // c_null_char &      ! isformat_struct
       // "Xcrysden axsf file" // c_null_char &      ! isformat_axsf
       // "Xcrysden xsf file" // c_null_char &       ! isformat_xsf
       // "xyz file" // c_null_char                  ! isformat_xyz
    nn = 0
    do i = 1,len(combostr_openfiles)
       if (combostr_openfiles(i:i) == c_null_char) nn = nn + 1
    end do

    !! Permutations for interpreting the format coming out of the dialog combo
    !! (same sequence as in combostr)
    allocate(isperm_openfiles(0:nn-1))
    isperm_openfiles(0:nn-1) = (/isformat_unknown,isformat_abinit,isformat_bincube,isformat_castepcell,&
       isformat_castepgeom,isformat_cif,isformat_crystal,isformat_cube,isformat_gen,isformat_f21,isformat_dmain,&
       isformat_elk,isformat_aimsin,isformat_aimsout,isformat_fploout,isformat_fchk,isformat_gjf,&
       isformat_gaussian,isformat_wfn,isformat_wfx,isformat_zmat,isformat_mol2,isformat_molden,isformat_orca,&
       isformat_pdb,isformat_pgout,isformat_dat,isformat_qein,isformat_qeout,isformat_pwc,isformat_shelx,&
       isformat_siesta,isformat_tinkerfrac,isformat_vasp,isformat_struct,isformat_axsf,&
       isformat_xsf,isformat_xyz/)

    !! Inverse of the permutation above
    allocate(isperm_inv_openfiles(0:maxval(isperm_openfiles)))
    isperm_inv_openfiles = 0
    do i = 0, nn-1
       isperm_inv_openfiles(isperm_openfiles(i)) = i
    end do

    ! open field file
    !! for the dialog: list of extensions "seen"
    dialogstr_openfieldfile =&
                "&
                &All files (*.*){*.*},&
                &ABINIT (DEN...){(DEN|ELF|POT|VHA\VHXC|VXC|GDEN1|GDEN2|GDEN3|LDEN|KDEN|PAWDEN|VCLMB|VPSP)},&
                &Aimpac (qub){.qub},&
                &Binary cube file (bincube){.bincube},&
                &CASTEP grid file (fmt){.fmt},&
                &Cube file (cube){.cube},&
                &DFTB+ (detailed.xml){.xml},&
                &elk (grid){.grid},&
                &elk STATE.OUT (OUT){.OUT},&
                &FPLO grid (001){.001},&
                &Gaussian wavefunction (wfn|wfx|fchk){.wfn,.wfx,.fchk},&
                &Molden-style file (molden){.molden},&
                &Quantum ESPRESSO pwc file (pwc){.pwc},&
                &SIESTA (RHO...) {(RHO|BADER|DRHO|LDOS|VT|VH)},&
                &VASP ELFCAR{(ELFCAR)},&
                &VASP (CHGCAR...){(CONTCAR|CHGCAR|ELFCAR|CHG|AECCAR0|AECCAR1|AECCAR2|POSCAR)},&
                &WIEN2k (clmsum...){.clmsum,.clmup,.clmdn},&
                &Xcrysden (xsf|axsf) {.xsf,.axsf},&
                &"// c_null_char

    combostr_openfieldfile = "" &
       // "Auto-detect" // c_null_char &              ! ifformat_unknown
       // "ABINIT DEN-style file" // c_null_char &    ! ifformat_abinit
       // "Aimpac qub" // c_null_char &               ! ifformat_qub
       // "Binary cube" // c_null_char &              ! ifformat_bincube
       // "CASTEP grid file" // c_null_char &         ! ifformat_fmt
       // "Cube file" // c_null_char &                ! ifformat_cube
       // "DFTB+ detailed.xml" // c_null_char &       ! ifformat_dftb
       // "elk grid" // c_null_char &                 ! ifformat_elkgrid
       // "elk STATE.OUT" // c_null_char &            ! ifformat_elk
       // "FPLO grid" // c_null_char &                ! ifformat_fplogrid
       // "Gaussian wfn" // c_null_char &             ! ifformat_wfn
       // "Gaussian wfx" // c_null_char &             ! ifformat_wfx
       // "Gaussian fchk" // c_null_char &            ! ifformat_fchk
       // "Molden file" // c_null_char &              ! ifformat_molden
       // "Quantum ESPRESSO pwc" // c_null_char &     ! ifformat_pwc
       // "SIESTA RHO-style file" // c_null_char &    ! ifformat_siestagrid
       // "VASP ELFCAR-style file" // c_null_char &   ! ifformat_vaspnov
       // "VASP CHGCAR-style file" // c_null_char &   ! ifformat_vasp
       // "WIEN2k clmsum-style file" // c_null_char & ! ifformat_wien
       // "Xcrysden xsf" // c_null_char               ! ifformat_xsf
    nn = 0
    do i = 1,len(combostr_openfieldfile)
       if (combostr_openfieldfile(i:i) == c_null_char) nn = nn + 1
    end do

    !! Permutations for interpreting the format coming out of the dialog combo
    !! (same sequence as in combostr)
    allocate(isperm_openfieldfile(0:nn-1))
    isperm_openfieldfile(0:nn-1) = (/ifformat_unknown,ifformat_abinit,ifformat_qub,ifformat_bincube,&
       ifformat_fmt,ifformat_cube,ifformat_dftb,ifformat_elkgrid,ifformat_elk,ifformat_fplogrid,&
       ifformat_wfn,ifformat_wfx,ifformat_fchk,ifformat_molden,ifformat_pwc,ifformat_siestagrid,&
       ifformat_vaspnov,ifformat_vasp,ifformat_wien,ifformat_xsf/)

    !! Inverse of the permutation above
    allocate(isperm_inv_openfieldfile(0:maxval(isperm_openfieldfile)))
    isperm_inv_openfieldfile = 0
    do i = 0, nn-1
       isperm_inv_openfieldfile(isperm_openfieldfile(i)) = i
    end do

    ! open vibration file
    dialogstr_openvibfile = &
       "&
       &All files (*.*){*.*},&
       &Quantum ESPRESSO matdyn.x modes file (modes){.modes},&
       &Quantum ESPRESSO eig file (eig){.eig},&
       &Quantum ESPRESSO ph.x dyn file (dyn){.dyn},&
       &phonopy ascii file (ascii){.ascii},&
       &phonopy yaml file (yaml){.yaml},&
       &phonopy hdf5 file (hdf5){.hdf5},&
       &crystal output file (out){.out},&
       &Gaussian output file (log){.log},&
       &Gaussian fchk file (fchk){.fchk},&
       &"// c_null_char
    combostr_openvibfile = "" &
       // "Auto-detect" // c_null_char &                 ! ivformat_unknown
       // "Quantum ESPRESSO modes file" // c_null_char & ! ivformat_matdynmodes
       // "Quantum ESPRESSO eig file" // c_null_char   & ! ivformat_matdyneig
       // "Quantum ESPRESSO dyn file" // c_null_char   & ! ivformat_qedyn
       // "phonopy ascii file" // c_null_char          & ! ivformat_phonopy_ascii
       // "phonopy yaml file" // c_null_char           & ! ivformat_phonopy_yaml
       // "phonopy hdf5 file" // c_null_char           & ! ivformat_phonopy_hdf5
       // "crystal output file" // c_null_char         & ! ivformat_crystal_out
       // "Gaussian output file" // c_null_char        & ! ivformat_gaussian_log
       // "Gaussian fchk file" // c_null_char            ! ivformat_gaussian_fchk
    nn = 0
    do i = 1,len(combostr_openvibfile)
       if (combostr_openvibfile(i:i) == c_null_char) nn = nn + 1
    end do

    allocate(isperm_openvibfile(0:nn-1))
    isperm_openvibfile(0:nn-1) = (/ivformat_unknown,ivformat_matdynmodes,&
       ivformat_matdyneig,ivformat_qedyn,ivformat_phonopy_ascii,ivformat_phonopy_yaml,&
       ivformat_phonopy_hdf5,ivformat_crystal_out,ivformat_gaussian_log,&
       ivformat_gaussian_fchk/)

    allocate(isperm_inv_openvibfile(0:maxval(isperm_openvibfile)))
    isperm_inv_openvibfile = 0
    do i = 0, nn-1
       isperm_inv_openvibfile(isperm_openvibfile(i)) = i
    end do

  end subroutine windows_init

  !> Deallocate all data in a command and reset it to emtpy
  module subroutine command_end(c)
    class(command_inout), intent(inout) :: c

    c%id = 0
    c%status = command_inout_empty
    c%size = 0
    c%scrolly = 0._c_float
    if (allocated(c%tooltipinfo)) deallocate(c%tooltipinfo)
    if (allocated(c%input)) deallocate(c%input)
    if (allocated(c%output)) deallocate(c%output)

  end subroutine command_end

  !xx! Module functions and subroutines

  !> Reallocate the window stack if there are only a few windows
  !> left. This must be done *outside* any window method because it
  !> will change the location of all window pointers.
  module subroutine stack_realloc_maybe()
    type(window), allocatable :: winaux(:)

    integer, parameter :: iroom = 5

    if (nwin+iroom > size(win,1)) then
       allocate(winaux(2*(nwin+iroom)))
       winaux(1:size(win,1)) = win
       call move_alloc(winaux,win)
       call regenerate_window_pointers()
    end if

  end subroutine stack_realloc_maybe

  !> Create a window in the window stack with the given type. Returns
  !> the window ID. If isopen, initialize the window as open. purpose =
  !> for dialogs, purpose of the dialog. isys = associated system.
  !> irep = associated representation. idparent = window ID of the
  !> caller. permanent = do not kill the window when closed.
  !> orraise = if this is a valid window ID, raise the window instead
  !> of creating a new one; if orraise < 0, raise the first window
  !> from the stack with the same type as type.
  module function stack_create_window(type,isopen,purpose,isys,irep,idparent,itoken,permanent,orraise)
    use tools_io, only: ferror, faterr
    use windows, only: window, nwin, win
    integer, intent(in) :: type
    logical, intent(in) :: isopen
    integer, intent(in), optional :: purpose
    integer, intent(in), optional :: isys
    integer, intent(in), optional :: irep
    integer, intent(in), optional :: idparent
    integer, intent(in), optional :: itoken
    logical, intent(in), optional :: permanent
    integer, intent(in), optional :: orraise

    integer :: stack_create_window

    integer :: i, id, raiseid
    logical :: ok
    integer, parameter :: maxwin = 40

    ! If orraise and the window exists, raise it, update parameters, and exit
    if (present(orraise)) then
       if (orraise < 0) then
          raiseid = 0
          do i = 1, nwin
             ok = win(i)%type == type .and. win(i)%isopen
             ! specific tests according to type
             if (ok.and.type == wintype_dialog.and.present(purpose)) ok = (win(i)%dialog_data%purpose == purpose)
             if (ok.and.type == wintype_editrep.and.present(isys).and.present(irep).and.present(idparent)) then
                ok = (win(i)%isys == isys .and. win(i)%irep == irep .and. win(i)%idparent == idparent)
                if (ok.and.present(itoken)) &
                   ok = (win(i)%itoken == itoken)
             end if
             if (ok.and.type == wintype_scfplot.and.present(isys)) ok = (win(i)%isys == isys)
             if (ok.and.type == wintype_geometry.and.present(isys)) ok = (win(i)%isys == isys)
             if (ok.and.type == wintype_rebond.and.present(isys)) ok = (win(i)%isys == isys)
             if (ok) then
                raiseid = i
                exit
             end if
          end do
       else
          raiseid = orraise
       end if
       if (raiseid > 0 .and. raiseid <= nwin) then
          stack_create_window = raiseid
          if (present(isys)) win(raiseid)%isys = isys
          if (present(irep)) win(raiseid)%irep = irep
          if (present(idparent)) win(raiseid)%idparent = idparent
          if (present(itoken)) win(raiseid)%itoken = itoken
          call igSetWindowFocus_Str(c_loc(win(raiseid)%name))
          return
       end if
    end if

    ! find the first unused window or create a new one
    id = 0
    do i = 1, nwin
       if (.not.win(i)%isinit.and..not.win(i)%permanent) then
          id = i
          exit
       end if
    end do
    if (id == 0) then
       nwin = nwin + 1
       id = nwin
    end if

    ! reallocate if necessary
    if (.not.allocated(win)) allocate(win(maxwin))

    if (nwin > size(win,1)) &
       call ferror('stack_create_window','too many windows',faterr)

    ! initialize the new window
    call win(id)%init(type,isopen,id,purpose,isys,irep,idparent,itoken)
    if (present(permanent)) then
       win(id)%permanent = permanent
    else
       win(id)%permanent = .false.
    end if
    stack_create_window = id

  end function stack_create_window

  !> This routine regenerates all pointers to the widows in the win(:)
  !> structure and its components. It is used when an array size is
  !> exceeded and move_alloc needs to be used to allocate more memory.
  module subroutine regenerate_window_pointers()
    use gui_main, only: sysc, sys_init, ok_system

    integer :: i, iv

    do i = 1, nwin
       if (win(i)%type == wintype_editrep .and. win(i)%irep > 0) then
          win(i)%rep => win(win(i)%idparent)%sc%rep(win(i)%irep)
       elseif (win(i)%type == wintype_view) then
          win(i)%sc => null()
          iv = win(i)%view_selected
          if (ok_system(iv,sys_init)) then
             if (win(i)%ismain) win(i)%sc => sysc(iv)%sc
          end if
       end if
    end do

  end subroutine regenerate_window_pointers

  !> Initialize a window of the given type. If isiopen, initialize it
  !> as open.
  module subroutine window_init(w,type,isopen,id,purpose,isys,irep,idparent,itoken)
    use interfaces_opengl3
    use gui_main, only: ColorDialogDir, ColorDialogFile
    use tools_io, only: ferror, faterr
    use param, only: bohrtoa
    class(window), intent(inout), target :: w
    integer, intent(in) :: type
    logical, intent(in) :: isopen
    integer, intent(in) :: id
    integer, intent(in), optional :: purpose
    integer, intent(in), optional :: isys
    integer, intent(in), optional :: irep
    integer, intent(in), optional :: idparent
    integer, intent(in), optional :: itoken

    character(kind=c_char,len=:), allocatable, target :: str1

    ! initialization of the state
    w%isinit = .true.
    w%firstpass = .true.
    w%isopen = isopen
    w%type = type
    w%id = -id
    w%name = "" // c_null_char
    w%errmsg = ""
    w%tree_selected = 1
    w%tree_sortcid = 0
    w%tree_sortdir = 1
    w%forceresize = .false.
    w%forcesort = .false.
    w%forceupdate = .false.
    w%forceinit = .false.
    w%forceselect = 0
    w%inpcon_selected = 1
    w%okfile = ""
    w%okfile_set = .false. ! whether the library file has been set by the user
    w%okfile_read = .false. ! whether the structure list should be re-read from the lib
    w%view_selected = 1
    nullify(w%sc)
    w%view_mousebehavior = MB_navigation
    w%forcerender = .true.
    if (allocated(w%iord)) deallocate(w%iord)
    w%dialog_data%dptr = c_null_ptr
    w%dialog_data%mol = -1
    w%dialog_data%showhidden = .false._c_bool
    w%dialog_data%isformat = isformat_unknown
    w%dialog_data%readlastonly = .false._c_bool
    w%dialog_data%purpose = wpurp_unknown
    w%dialog_data%molcubic = .false.
    w%dialog_data%rborder = rborder_def*bohrtoa
    w%plotn = 0
    if (present(isys)) w%isys = isys
    if (present(irep)) w%irep = irep
    if (present(idparent)) w%idparent = idparent
    if (present(itoken)) w%itoken = itoken
    if (present(purpose)) w%dialog_purpose = purpose

    ! type-specific initialization
    if (type == wintype_dialog) then
       ! dialog
       w%dptr = IGFD_Create()
       str1 = "+" // c_null_char
       call IGFD_SetFileStyle(w%dptr,IGFD_FileStyleByTypeDir,c_null_ptr,ColorDialogDir,c_loc(str1),c_null_ptr)
       str1 = " " // c_null_char
       call IGFD_SetFileStyle(w%dptr,IGFD_FileStyleByTypeFile,c_null_ptr,ColorDialogFile,c_loc(str1),c_null_ptr)
       if (.not.present(purpose)) &
          call ferror('window_init','dialog requires a purpose',faterr)
    elseif (type == wintype_load_field) then
       ! dialog: load field window
       if (.not.present(isys)) &
          call ferror('window_init','load_field requires isys',faterr)
    elseif (type == wintype_scfplot) then
       ! SCF plot window
       if (.not.present(isys)) &
          call ferror('window_init','scfplot requires isys',faterr)
    elseif (type == wintype_editrep) then
       ! edit representation window
       if (.not.present(isys)) &
          call ferror('window_init','editrep requires isys',faterr)
       if (.not.present(irep)) &
          call ferror('window_init','editrep requires irep',faterr)
       if (.not.present(idparent)) &
          call ferror('window_init','editrep requires idparent',faterr)
       w%rep => win(idparent)%sc%rep(irep)
    elseif (type == wintype_exportimage) then
       ! export image window
       if (.not.present(idparent)) &
          call ferror('window_init','exportimage requires idparent',faterr)
    elseif (type == wintype_exportimage) then
       ! vibrations window
       if (.not.present(idparent)) &
          call ferror('window_init','vibrations requires idparent',faterr)
    elseif (type == wintype_view) then
       ! view window
       if (.not.present(purpose)) &
          call ferror('window_init','view requires purpose',faterr)
       if (purpose == wpurp_view_main) then
          w%ismain = .true.
       elseif (purpose == wpurp_view_alternate) then
          w%ismain = .false.
          allocate(w%sc)
       else
          call ferror('window_init','unknown view purpose',faterr)
       end if
       call w%create_texture_view(initial_texture_side)
    elseif (type == wintype_rebond) then
       ! recalculate bonds window
       if (.not.present(isys)) &
          call ferror('window_init','rebond requires isys',faterr)
    elseif (type == wintype_geometry) then
       ! geometry window
       if (.not.present(isys)) &
          call ferror('window_init','geometry requires isys',faterr)
    end if

  end subroutine window_init

  !> End a window and deallocate the data.
  module subroutine window_end(w)
    use gui_main, only: ok_system, sysc, sys_init
    use interfaces_opengl3
    class(window), intent(inout), target :: w

    ! window-specific destruction
    if (w%isinit) then
       if (w%type == wintype_dialog .and. c_associated(w%dptr)) then
          ! destroy the dialog
          call IGFD_Destroy(w%dptr)
       elseif (w%type == wintype_view) then
          ! delete the texture and deallocate the scene if not the main view
          call w%delete_texture_view()
          if (.not.w%ismain) then
             if (associated(w%sc)) deallocate(w%sc)
          end if
       elseif (w%type == wintype_vibrations) then
          ! reset the animation status of the parent
          if (w%idparent > 0 .and. w%idparent <= nwin) then
             if (associated(win(w%idparent)%sc)) then
                win(w%idparent)%sc%iqpt_selected = 0
                win(w%idparent)%sc%ifreq_selected = 0
                win(w%idparent)%sc%animation = 0
                win(w%idparent)%forcerender = .true.
             end if
          end if
       elseif (w%type == wintype_geometry) then
          ! remove all highlights
          if (allocated(w%geometry_selected)) w%geometry_selected = .false.
          if (ok_system(w%isys,sys_init)) &
             call sysc(w%isys)%highlight_clear(.false.)
       end if
    end if

    ! deallocate the rest of the data
    w%firstpass = .true.
    w%isinit = .false.
    w%isopen = .false.
    w%id = -1
    w%idparent = 0
    w%itoken = 0
    w%isys = 0
    w%irep = 0
    w%tied_to_tree = .false.
    w%name = "" // c_null_char
    if (allocated(w%iord)) deallocate(w%iord)
    if (allocated(w%forceremove)) deallocate(w%forceremove)
    w%plotn = 0
    if (allocated(w%plotx)) deallocate(w%plotx)
    if (allocated(w%ploty)) deallocate(w%ploty)
    nullify(w%rep)
    w%irep = 0

  end subroutine window_end

  !> Return true if the root of this window is focused
  module function window_focused(w)
    use gui_main, only: g
    class(window), intent(inout) :: w
    logical :: window_focused

    type(ImGuiWindow), pointer :: wptr, navwin, navwinroot

    window_focused = .false.
    if (.not.c_associated(w%ptr).or..not.c_associated(g%NavWindow)) return
    call c_f_pointer(w%ptr,wptr)
    call c_f_pointer(g%NavWindow,navwin)
    if (.not.c_associated(navwin%RootWindow)) return
    call c_f_pointer(navwin%RootWindow,navwinroot)
    window_focused = associated(wptr,navwinroot)

  end function window_focused

  !> Draw an ImGui window.
  module subroutine window_draw(w)
    use gui_main, only: fontsize
    use utils, only: iw_text, get_nice_next_window_pos
    use tools_io, only: string, ferror, faterr
    class(window), intent(inout), target :: w

    character(kind=c_char,len=:), allocatable, target :: str1, str2, str3
    type(ImVec2) :: inisize, pos, pivot
    type(ImGuiWindow), pointer :: wptr

    if (.not.w%isinit) return
    if (.not.w%isopen.and..not.w%permanent) then
       ! the window may have been closed in imgui (x button)
       call w%end()
       return
    end if

    ! First pass on creation: assign ID, name, and flags
    if (w%id < 0) then
       w%firstpass = .true.
       w%id = abs(w%id)
       w%pos = 0._c_float
       w%isdocked = .false.
       if (w%type == wintype_tree) then
          w%name = "Tree##" // string(w%id) // c_null_char
          w%flags = ImGuiWindowFlags_None
       elseif (w%type == wintype_view) then
          w%flags = ImGuiWindowFlags_None
          if (w%ismain) then
             w%name = "Main View##" // string(w%id)  // c_null_char
          else
             w%name = "Alternate View##" // string(w%id)  // c_null_char
             inisize%x = 90 * fontsize%x
             inisize%y = 30 * fontsize%y
             call igSetNextWindowSize(inisize,ImGuiCond_FirstUseEver)
             call get_nice_next_window_pos(pos)
             pivot%x = 0._c_float
             pivot%y = 0._c_float
             call igSetNextWindowPos(pos,ImGuiCond_None,pivot)
          end if
       elseif (w%type == wintype_console_input) then
          w%name = "Input Console##" // string(w%id)  // c_null_char
          w%flags = ImGuiWindowFlags_None
       elseif (w%type == wintype_console_output) then
          w%name = "Output Console##" // string(w%id)  // c_null_char
          w%flags = ImGuiWindowFlags_None
       elseif (w%type == wintype_about) then
          w%name = "About##" // string(w%id)  // c_null_char
          w%flags = ImGuiWindowFlags_None
          inisize%x = 52 * fontsize%x
          inisize%y = 17 * fontsize%y
          call igSetNextWindowSize(inisize,ImGuiCond_FirstUseEver)
       elseif (w%type == wintype_dialog) then
          w%dialog_data%dptr = w%dptr
          w%dialog_data%purpose = w%dialog_purpose
          if (w%dialog_data%showhidden) then
             w%flags = ImGuiFileDialogFlags_None
          else
             w%flags = ImGuiFileDialogFlags_DontShowHiddenFiles
          end if
          str2 = "" // c_null_char ! default path
          inisize%x = 90 * fontsize%x
          inisize%y = 30 * fontsize%y
          call igSetNextWindowSize(inisize,ImGuiCond_FirstUseEver)

          str1 = "All files (*.*){*.*}" // c_null_char
          if (w%dialog_purpose == wpurp_dialog_openfiles) then
             ! open dialog
             w%name = "Open File(s)##" // string(w%id)  // c_null_char
             call IGFD_OpenPaneDialog2(w%dptr,c_loc(w%name),c_loc(w%name),c_loc(dialogstr_openfiles),c_loc(str2),&
                c_funloc(dialog_user_callback),280._c_float,0_c_int,c_loc(w%dialog_data),w%flags)
          elseif (w%dialog_purpose == wpurp_dialog_savelogfile) then
             w%name = "Save Log File##" // string(w%id)  // c_null_char
             str2 = "file.log" // c_null_char
             str3 = "./" // c_null_char
             call IGFD_OpenPaneDialog(w%dptr,c_loc(w%name),c_loc(w%name),c_loc(str1),c_loc(str3),c_loc(str2),&
                c_funloc(dialog_user_callback),280._c_float,1_c_int,c_loc(w%dialog_data),w%flags)
          elseif (w%dialog_purpose == wpurp_dialog_openlibraryfile) then
             w%name = "Open Library File##" // string(w%id)  // c_null_char
             call IGFD_OpenPaneDialog2(w%dptr,c_loc(w%name),c_loc(w%name),c_loc(str1),c_loc(str2),&
                c_funloc(dialog_user_callback),280._c_float,1_c_int,c_loc(w%dialog_data),w%flags)
          elseif (w%dialog_purpose == wpurp_dialog_openfieldfile) then
             w%name = "Open Field File(s)##" // string(w%id)  // c_null_char
             call IGFD_OpenPaneDialog2(w%dptr,c_loc(w%name),c_loc(w%name),c_loc(dialogstr_openfieldfile),&
                c_loc(str2),c_funloc(dialog_user_callback),280._c_float,1_c_int,c_loc(w%dialog_data),w%flags)
          elseif (w%dialog_purpose == wpurp_dialog_openvibfile) then
             w%name = "Open Vibration Data File(s)##" // string(w%id)  // c_null_char
             call IGFD_OpenPaneDialog2(w%dptr,c_loc(w%name),c_loc(w%name),c_loc(dialogstr_openvibfile),&
                c_loc(str2),c_funloc(dialog_user_callback),280._c_float,0_c_int,c_loc(w%dialog_data),w%flags)
          elseif (w%dialog_purpose == wpurp_dialog_openonefilemodal) then
             w%name = "Open File##" // string(w%id)  // c_null_char
             call IGFD_OpenPaneDialog2(w%dptr,c_loc(w%name),c_loc(w%name),c_loc(str1),c_loc(str2),&
                c_funloc(dialog_user_callback),280._c_float,1_c_int,c_loc(w%dialog_data),w%flags)
          elseif (w%dialog_purpose == wpurp_dialog_saveimagefile) then
             w%name = "Save Image File##" // string(w%id) // c_null_char
             str1 = "PNG (*.png) {.png},BMP (*.bmp) {.bmp},TGA (*.tga) {.tga},JPEG (*.jpg) {.jpg}"// c_null_char
             str2 = "image.png" // c_null_char
             str3 = "./" // c_null_char
             call IGFD_OpenPaneDialog(w%dptr,c_loc(w%name),c_loc(w%name),c_loc(str1),c_loc(str3),c_loc(str2),&
                c_funloc(dialog_user_callback),280._c_float,1_c_int,c_loc(w%dialog_data),w%flags)
          else
             call ferror('window_draw','unknown dialog purpose',faterr)
          end if
       elseif (w%type == wintype_new_struct) then
          w%name = "New Structure##" // string(w%id)  // c_null_char
          w%flags = ImGuiWindowFlags_None
          inisize%x = 90 * fontsize%x
          inisize%y = 40 * fontsize%y
          call igSetNextWindowSize(inisize,ImGuiCond_FirstUseEver)
       elseif (w%type == wintype_new_struct_library) then
          w%name = "New Structure from Library##" // string(w%id)  // c_null_char
          w%flags = ImGuiWindowFlags_None
          inisize%x = 90 * fontsize%x
          inisize%y = 30 * fontsize%y
          call igSetNextWindowSize(inisize,ImGuiCond_FirstUseEver)
       elseif (w%type == wintype_load_field) then
          w%name = "Load Field##" // string(w%id) // c_null_char
          w%flags = ImGuiWindowFlags_None
          inisize%x = 65 * fontsize%x
          inisize%y = 16 * fontsize%y
          call igSetNextWindowSize(inisize,ImGuiCond_FirstUseEver)
       elseif (w%type == wintype_scfplot) then
          w%name = "SCF Iterations (" // string(w%isys) // ")##" // string(w%id)  // c_null_char
          w%flags = ImGuiWindowFlags_None
          inisize%x = 45 * fontsize%x
          inisize%y = inisize%x
          call igSetNextWindowSize(inisize,ImGuiCond_FirstUseEver)
       elseif (w%type == wintype_editrep) then
          w%name = "Object [" // string(w%rep%name) // "]##" // string(w%id) // c_null_char
          w%flags = ImGuiWindowFlags_None
          inisize%x = 63 * fontsize%x
          inisize%y = 46 * fontsize%y
          call igSetNextWindowSize(inisize,ImGuiCond_FirstUseEver)
       elseif (w%type == wintype_exportimage) then
          w%name = "Export to Image" // "##" // string(w%id) // c_null_char
          w%flags = ImGuiWindowFlags_None
          inisize%x = 50 * fontsize%x
          inisize%y = 17 * fontsize%y
          call igSetNextWindowSize(inisize,ImGuiCond_FirstUseEver)
       elseif (w%type == wintype_vibrations) then
          w%name = "Vibrations" // "##" // string(w%id) // c_null_char
          w%flags = ImGuiWindowFlags_None
          inisize%x = 62 * fontsize%x
          inisize%y = 25 * fontsize%y
          call igSetNextWindowSize(inisize,ImGuiCond_FirstUseEver)
       elseif (w%type == wintype_rebond) then
          w%name = "Recalculate Bonds###rebond"  // string(w%id) // c_null_char
          w%flags = ImGuiWindowFlags_None
          inisize%x = 55 * fontsize%x
          inisize%y = 23 * fontsize%y
          call igSetNextWindowSize(inisize,ImGuiCond_FirstUseEver)
       elseif (w%type == wintype_geometry) then
          w%name = "View/Edit Geometry##"  // string(w%id) // c_null_char
          w%flags = ImGuiWindowFlags_None
          inisize%x = 55 * fontsize%x
          inisize%y = 29 * fontsize%y
          call igSetNextWindowSize(inisize,ImGuiCond_FirstUseEver)
       elseif (w%type == wintype_preferences) then
          w%name = "Preferences##" // string(w%id)  // c_null_char
          w%flags = ImGuiWindowFlags_None
          inisize%x = 65 * fontsize%x
          inisize%y = 23 * fontsize%y
          call igSetNextWindowSize(inisize,ImGuiCond_FirstUseEver)
       elseif (w%type == wintype_treeplot) then
          w%name = "Plot Tree Data##" // string(w%id)  // c_null_char
          w%flags = ImGuiWindowFlags_None
          inisize%x = 55 * fontsize%x
          inisize%y = 25 * fontsize%y
          call igSetNextWindowSize(inisize,ImGuiCond_FirstUseEver)
       end if
    end if

    ! Window pre-processing: some windows (scf plot, edit rep)
    ! require having a valid system/representation/etc.
    ! associated
    if (w%isopen) then
       if (w%type == wintype_scfplot) then
          call w%update_scfplot()
       elseif (w%type == wintype_editrep) then
          call w%update_editrep()
       elseif (w%type == wintype_geometry) then
          call w%update_geometry()
       elseif (w%type == wintype_rebond) then
          call w%update_rebond()
       end if
    end if

    ! Draw the window contents, depending on type
    if (w%isopen) then
       ! assign the pointer ID for the window, if not a dialog
       if (w%type == wintype_dialog) then
          w%ptr = IGFD_GetCurrentWindow(w%dptr)
          call w%draw_dialog()
       else
          if (igBegin(c_loc(w%name),w%isopen,w%flags)) then
             w%ptr = igGetCurrentWindow()
             if (w%type == wintype_tree) then
                call w%draw_tree()
             elseif (w%type == wintype_view) then
                call w%draw_view()
             elseif (w%type == wintype_console_input) then
                call w%draw_ci()
             elseif (w%type == wintype_console_output) then
                call w%draw_co()
             elseif (w%type == wintype_about) then
                call w%draw_about()
             elseif (w%type == wintype_new_struct) then
                call w%draw_new_struct()
             elseif (w%type == wintype_new_struct_library) then
                call w%draw_new_struct_from_library()
             elseif (w%type == wintype_load_field) then
                call w%draw_load_field()
             elseif (w%type == wintype_scfplot) then
                call w%draw_scfplot()
             elseif (w%type == wintype_editrep) then
                call w%draw_editrep()
             elseif (w%type == wintype_exportimage) then
                call w%draw_exportimage()
             elseif (w%type == wintype_vibrations) then
                call w%draw_vibrations()
             elseif (w%type == wintype_rebond) then
                call w%draw_rebond()
             elseif (w%type == wintype_geometry) then
                call w%draw_geometry()
             elseif (w%type == wintype_preferences) then
                call w%draw_preferences()
             elseif (w%type == wintype_treeplot) then
                call w%draw_treeplot()
             end if
          end if
          call igEnd()
       end if
       w%firstpass = .false.

       ! write down window info
       if (c_associated(w%ptr)) then
          call c_f_pointer(w%ptr,wptr)
          w%pos(1) = wptr%Pos%x
          w%pos(2) = wptr%Pos%y
          w%isdocked = wptr%DockIsActive
       end if
    else
       w%firstpass = .true.
    end if

  end subroutine window_draw

  !> Draw the about window
  module subroutine draw_about(w)
    use config, only: istring_version, getstring
    use gui_main, only: g
    use utils, only: iw_text, iw_button, iw_calcwidth
    use keybindings, only: is_bind_event, BIND_CLOSE_FOCUSED_DIALOG, BIND_CLOSE_ALL_DIALOGS,&
       BIND_OK_FOCUSED_DIALOG
    class(window), intent(inout), target :: w

    real(c_float) :: wwidth, twidth
    type(ImVec2) :: szavail

    call iw_text("  ---- critic2 GUI ----",danger=.true.,centered=.true.)
    call iw_text("version " // getstring(istring_version),centered=.true.)
    call igNewLine()
    call iw_text("Critic2 is a program for visualization and analysis",centered=.true.)
    call iw_text("of structures and data in computational chemistry",centered=.true.)
    call igNewLine()
    call iw_text("Copyright (c) Alberto Otero de la Roza",centered=.true.)
    call iw_text("Distributed under GNU/GPL license, version 3",centered=.true.)
    call iw_text("Contact: aoterodelaroza@gmail.com",highlight=.true.,centered=.true.)
    call igNewLine()

    ! bottom-align for the rest of the contents
    call igGetContentRegionAvail(szavail)
    if (szavail%y > igGetTextLineHeightWithSpacing() + g%Style%WindowPadding%y) &
       call igSetCursorPosY(igGetCursorPosY() + szavail%y - igGetTextLineHeightWithSpacing() - g%Style%WindowPadding%y)
    wwidth = igGetWindowWidth()
    twidth = iw_calcwidth(5,1)
    call igSetCursorPosX((wwidth - twidth) * 0.5_c_float)
    if (iw_button("Close")) w%isopen = .false.

    ! exit if focused and received the close keybinding
    if ((w%focused() .and. (is_bind_event(BIND_CLOSE_FOCUSED_DIALOG) .or. is_bind_event(BIND_OK_FOCUSED_DIALOG))) .or.&
       is_bind_event(BIND_CLOSE_ALL_DIALOGS))&
       w%isopen = .false.

  end subroutine draw_about

  !> Draw the preferences window
  module subroutine draw_preferences(w)
    use gui_main, only: g, tooltip_enabled, tooltip_delay, tooltip_wrap_factor,&
       tree_select_updates_inpcon, tree_select_updates_view, io, fontsize,&
       set_default_ui_settings, ColorTableCellBg
    use interfaces_cimgui
    use keybindings
    use utils, only: iw_tooltip, iw_button, iw_text, iw_calcwidth, iw_clamp_color4,&
       iw_checkbox
    class(window), intent(inout), target :: w

    character(kind=c_char,len=:), allocatable, target :: str, str2, zeroc
    character(len=:), allocatable, target :: strf
    logical(c_bool) :: ldum
    logical :: doquit
    type(ImVec2) :: sz, szero
    integer :: i, newkey
    integer(c_int) :: flags

    logical, save :: ttshown = .false. ! tooltip flag
    type(c_ptr), save :: cfilter = c_null_ptr ! filter object (allocated first pass, never destroyed)
    integer(c_int), save :: catid = 0 ! category ID (from left panel) 0=interface,1=keybinding,2=colors
    integer(c_int), save :: getbind = -1 ! get binding flag

    real(c_float), parameter :: wleft = 120._c_float

    ! first pass when opened, reset the state
    if (w%firstpass) call init_state()

    ! initialize
    doquit = .false.
    szero%x = 0
    szero%y = 0

    ! text filter
    call igAlignTextToFramePadding();
    call iw_text("Filter")
    call igSameLine(0._c_float,-1._c_float)
    if (.not.c_associated(cfilter)) &
       cfilter = ImGuiTextFilter_ImGuiTextFilter(c_loc(zeroc))
    str = "##preferencesfilter" // c_null_char
    ldum = ImGuiTextFilter_Draw(cfilter,c_loc(str),0._c_float)
    call iw_tooltip("Filter UI settings by name in the list below. Use comma-separated fields&
       & and - for excluding. Example: inc1,inc2,-exc includes all settings&
       & with inc1 or inc2 and excludes settings with exc.",ttshown)
    if (iw_button("Clear",sameline=.true.)) then
       if (c_associated(cfilter)) &
          call ImGuiTextFilter_Clear(cfilter)
    end if
    call iw_tooltip("Clear the filter",ttshown)

    ! left panel with selectables
    str = "leftpanel" // c_null_char
    sz%x = wleft
    sz%y = 0._c_float
    if (igBeginChild_Str(c_loc(str),sz,.true._c_bool,ImGuiWindowFlags_None)) then
       str = "Interface" // c_null_char
       if (igSelectable_Bool(c_loc(str),logical(catid == 0,c_bool),ImGuiSelectableFlags_None,szero)) catid = 0
       call iw_tooltip("Global settings for the graphical user interface",ttshown)
       str = "Key bindings" // c_null_char
       if (igSelectable_Bool(c_loc(str),logical(catid == 1,c_bool),ImGuiSelectableFlags_None,szero)) catid = 1
       call iw_tooltip("User interface key bindings",ttshown)
       str = "Tree colors" // c_null_char
       if (igSelectable_Bool(c_loc(str),logical(catid == 2,c_bool),ImGuiSelectableFlags_None,szero)) catid = 2
       call iw_tooltip("Background colors for the systems in the tree window",ttshown)
    end if
    call igEndChild()
    call igSameLine(0._c_float,-1._c_float)

    ! right panel
    call igBeginGroup()
    str = "rightpanel" // c_null_char
    sz%x = 0._c_float
    sz%y = -igGetFrameHeightWithSpacing() - g%Style%ItemSpacing%y
    if (igBeginChild_Str(c_loc(str),sz,.false._c_bool,ImGuiWindowFlags_None)) then

       call igPushItemWidth(6._c_float * fontsize%x)
       if (catid == 0) then
          !! Interface
          call iw_text("Interface",highlight=.true.)
          call igSeparator()

          str = "Font scale" // c_null_char
          str2 = "%.2f" // c_null_char
          if (ImGuiTextFilter_PassFilter(cfilter,c_loc(str),c_null_ptr)) then
             ldum = igDragFloat(c_loc(str),io%FontGlobalScale,0.005_c_float,0.3_c_float,3.0_c_float,c_loc(str2),&
                ImGuiSliderFlags_AlwaysClamp)
             call iw_tooltip("Scale factor for the user interface font",ttshown)
          end if

          str = "Enable tooltips" // c_null_char
          if (ImGuiTextFilter_PassFilter(cfilter,c_loc(str),c_null_ptr)) then
             ldum = iw_checkbox(str,tooltip_enabled)
             call iw_tooltip("Show/hide the tooltips when hovering interface elements with the mouse",ttshown)
          end if

          str = "Tooltip delay (s)" // c_null_char
          str2 = "%.1f" // c_null_char
          if (ImGuiTextFilter_PassFilter(cfilter,c_loc(str),c_null_ptr)) then
             ldum = igDragFloat(c_loc(str),tooltip_delay,0.1_c_float,0._c_float,5._c_float,c_loc(str2),&
                ImGuiSliderFlags_AlwaysClamp)
             call iw_tooltip("Delay for showing the tooltips",ttshown)
          end if

          str = "Tooltip maximum width (pixels)" // c_null_char
          str2 = "%.1f" // c_null_char
          if (ImGuiTextFilter_PassFilter(cfilter,c_loc(str),c_null_ptr)) then
             ldum = igDragFloat(c_loc(str),tooltip_wrap_factor,1._c_float,0._c_float,1000._c_float,&
                c_loc(str2),ImGuiSliderFlags_AlwaysClamp)
             call iw_tooltip("Width of the interface tooltips",ttshown)
          end if

          str = "Tree selects input console system" // c_null_char
          if (ImGuiTextFilter_PassFilter(cfilter,c_loc(str),c_null_ptr)) then
             ldum = iw_checkbox(str,tree_select_updates_inpcon)
             call iw_tooltip("Selecting a system on the tree changes the system in the input console",ttshown)
          end if

          str = "Tree selects view system" // c_null_char
          if (ImGuiTextFilter_PassFilter(cfilter,c_loc(str),c_null_ptr)) then
             ldum = iw_checkbox(str,tree_select_updates_view)
             call iw_tooltip("Selecting a system on the tree changes the system in the view window",ttshown)
          end if

       elseif (catid == 1) then
          !! key bindings
          call iw_text("Key bindings",highlight=.true.)
          call iw_text("(?)",sameline=.true.)
          call iw_tooltip("Left click to assign a new binding. Right-click to toggle double click behavior&
             & (only for mouse input). Middle click to erase the binding.")
          call igSeparator()

          str = "##keybindtable" // c_null_char
          flags = ImGuiTableFlags_NoSavedSettings
          flags = ior(flags,ImGuiTableFlags_SizingFixedFit)
          flags = ior(flags,ImGuiTableFlags_NoBordersInBody)
          sz%x = 0
          sz%y = 0
          if (igBeginTable(c_loc(str),2,flags,sz,0._c_float)) then
             do i = 1, BIND_NUM
                str2 = trim(bindnames(i)) // c_null_char
                if (.not.ImGuiTextFilter_PassFilter(cfilter,c_loc(str2),c_null_ptr)) cycle

                call igTableNextRow(ImGuiTableRowFlags_None, 0._c_float)
                if (igTableSetColumnIndex(0)) then
                   call iw_text(str2)
                end if
                if (igTableSetColumnIndex(1)) then
                   strf = trim(get_bind_keyname(i))

                   call igPushID_Int(int(i,c_int))
                   if (iw_button(strf)) getbind = i
                   call igPopID()

                   if (igIsItemHovered(ImGuiHoveredFlags_None)) then
                      ! right click to toggle
                      if (igIsMouseClicked(ImGuiPopupFlags_MouseButtonRight,.false._c_bool)) then
                         newkey = ImGuiKey_None
                         if (keybind(i) == ImGuiKey_MouseLeft) then
                            newkey = ImGuiKey_MouseLeftDouble
                         elseif (keybind(i) == ImGuiKey_MouseLeftDouble) then
                            newkey = ImGuiKey_MouseLeft
                         elseif (keybind(i) == ImGuiKey_MouseRight) then
                            newkey = ImGuiKey_MouseRightDouble
                         elseif (keybind(i) == ImGuiKey_MouseRightDouble) then
                            newkey = ImGuiKey_MouseRight
                         elseif (keybind(i) == ImGuiKey_MouseMiddle) then
                            newkey = ImGuiKey_MouseMiddleDouble
                         elseif (keybind(i) == ImGuiKey_MouseMiddleDouble) then
                            newkey = ImGuiKey_MouseMiddle
                         end if
                         if (newkey /= ImGuiKey_None) then
                            call set_bind(i,newkey,modbind(i))
                         end if
                      end if

                      ! middle click to erase
                      if (igIsMouseClicked(ImGuiPopupFlags_MouseButtonMiddle,.false._c_bool)) &
                         call set_bind(i,ImGuiKey_None,modbind(i))
                   end if
                end if
             end do
             call igEndTable()
          end if

          ! popup to read a key
          if (getbind /= -1) then
             use_keybindings = .false.
             str2 = "Choosekey" // c_null_char
             flags = ImGuiWindowFlags_AlwaysAutoResize
             flags = ior(flags,ImGuiWindowFlags_NoTitleBar)
             flags = ior(flags,ImGuiWindowFlags_NoMove)
             call igOpenPopup_Str(c_loc(str2),ImGuiPopupFlags_None)
             if (igBeginPopupModal(c_loc(str2),logical(.true.,c_bool),flags)) then
                call iw_text("Please press a key or mouse button.")
                if (set_bind_from_user_input(getbind)) then
                   getbind = -1
                   use_keybindings = .true.
                   call igCloseCurrentPopup()
                end if
                call igEndPopup()
             end if
          end if
       elseif (catid == 2) then
          !! tree colors
          call iw_text("Tree colors",highlight=.true.)
          call igSeparator()

          ! color editors
          call color_edit4("Crystal","A molecular system with a single molecule",ColorTableCellBg(:,0))
          call color_edit4("Layered crystal","A layered crystal, with two-dimensional periodic bonding networks",&
             ColorTableCellBg(:,1))
          call color_edit4("Crystal with 1D bonded chains",&
             "A crystal, with a one-dimensional periodic bonding networks",&
             ColorTableCellBg(:,2))
          call color_edit4("Molecular crystal","A molecular crystal",ColorTableCellBg(:,3))
          call color_edit4("Surface","A 2D surface",ColorTableCellBg(:,4))
          call color_edit4("Chain","An one-dimensional atomic structure (e.g. a polymer)",ColorTableCellBg(:,5))
          call color_edit4("Molecule in a box","A molecule in a large periodic box",ColorTableCellBg(:,6))
          call color_edit4("Molecule","A molecule",ColorTableCellBg(:,7))
          call color_edit4("Molecular cluster","A molecular cluster",ColorTableCellBg(:,8))
       end if
       call igPopItemWidth()
    end if
    call igEndChild()

    ! child for final buttons
    str = "finalbuttons" // c_null_char
    sz%x = 0._c_float
    sz%y = 0._c_float
    if (igBeginChild_Str(c_loc(str),sz,.false._c_bool,ImGuiWindowFlags_None)) then
       call igSetCursorPosX(iw_calcwidth(10,2,from_end=.true.) - g%Style%ScrollbarSize)
       if (iw_button("Reset",danger=.true.)) &
          call set_default_ui_settings()
       call iw_tooltip("Reset to the default settings",ttshown)
       if (iw_button("Close",sameline=.true.)) doquit = .true.
       call iw_tooltip("Close this window",ttshown)
    end if
    call igEndChild()
    call igEndGroup()

    ! exit if focused and received the close keybinding
    if (w%focused() .and. is_bind_event(BIND_OK_FOCUSED_DIALOG)) doquit = .true.
    if ((w%focused() .and. is_bind_event(BIND_CLOSE_FOCUSED_DIALOG)).or.&
       is_bind_event(BIND_CLOSE_ALL_DIALOGS)) doquit = .true.

    ! quit the window
    if (doquit) then
       call end_state()
       call w%end()
    end if

  contains
    ! initialize the state for this window
    subroutine init_state()
       if (c_associated(cfilter)) call ImGuiTextFilter_Clear(cfilter)
       catid = 0_c_int
       getbind = -1_c_int
    end subroutine init_state

    ! terminate the state for this window
    subroutine end_state()
       if (c_associated(cfilter)) call ImGuiTextFilter_Clear(cfilter)
    end subroutine end_state

    ! display the color editor
    subroutine color_edit4(str,strtt,rgba)
      character*(*), intent(in) :: str, strtt
      real(c_float), intent(inout) :: rgba(4)

      character(len=:), allocatable, target :: str_
      logical(c_bool) :: ldum

      str_ = str // c_null_char
      if (ImGuiTextFilter_PassFilter(cfilter,c_loc(str_),c_null_ptr)) then
         ldum = igColorEdit4(c_loc(str_),rgba,ImGuiColorEditFlags_NoInputs)
         call iw_tooltip(strtt,ttshown)
         call iw_clamp_color4(rgba)
      end if

    end subroutine color_edit4

  end subroutine draw_preferences

  !xx! private procedures

  ! the callback for the right-hand-side pane of the dialog
  subroutine dialog_user_callback(vFilter, vUserData, vCantContinue)
    use gui_main, only: g
    use utils, only: igIsItemHovered_delayed, iw_tooltip, iw_text, iw_radiobutton,&
       iw_combo_simple, iw_checkbox
    use interfaces_cimgui
    use tools_io, only: string
    type(c_ptr), intent(in), value :: vFilter ! const char *
    type(c_ptr), value :: vUserData ! void *
    logical(c_bool) :: vCantContinue ! bool *

    character(kind=c_char,len=:), allocatable, target :: str, stropt, strex
    type(dialog_userdata), pointer :: data
    logical(c_bool) :: ldum
    type(ImVec2) :: sz
    integer(c_int) :: flags

    logical, save :: ttshown = .false. ! tooltip flag

    ! generate the data pointer
    call c_f_pointer(vUserData,data)

    !! common options !!

    ! header
    call iw_text("Options",highlight=.true.)

    ! show hidden files
    if (iw_checkbox("Show hidden files",data%showhidden)) then
       flags = IGFD_GetFlags(data%dptr)
       if (data%showhidden) then
          flags = iand(flags,not(ImGuiFileDialogFlags_DontShowHiddenFiles))
       else
          flags = ior(flags,ImGuiFileDialogFlags_DontShowHiddenFiles)
       end if
       call IGFD_SetFlags(data%dptr,flags)
    end if
    call iw_tooltip("Show the hidden files and directories",ttshown)

    !! options specific to the open files dialog !!
    if (data%purpose == wpurp_dialog_openfiles) then
       ldum = iw_checkbox("Read last structure only",data%readlastonly)
       call iw_tooltip("Read only the last structure in the file",ttshown)
       call igNewLine()

       ! radio buttons for auto/crystal/molecule
       call iw_text("Read structures as...",highlight=.true.)
       ldum = iw_radiobutton("Auto-detect",int=data%mol,intval=-1_c_int)
       call iw_tooltip("Auto-detect whether new structures are read as crystals or molecules",ttshown)
       ldum = iw_radiobutton("Crystal",int=data%mol,intval=0_c_int)
       call iw_tooltip("Force new structures to be read as crystals",ttshown)
       ldum = iw_radiobutton("Molecule",int=data%mol,intval=1_c_int)
       call iw_tooltip("Force new structures to be read as molecules",ttshown)

       ! molecular options
       call igIndent(0._c_float)
       str = "Cell border (Å)" // c_null_char
       stropt = "%.3f" // c_null_char
       strex = string(data%rborder,'f',decimal=3) // c_null_char
       call igCalcTextSize(sz,c_loc(strex),c_null_ptr,.false._c_bool,-1._c_float)
       call igPushItemWidth(sz%x + 2 * g%Style%FramePadding%x)
       ldum = igInputFloat(c_loc(str),data%rborder,0._c_float,0._c_float,&
          c_loc(stropt),ImGuiInputTextFlags_None)
       call iw_tooltip("Size of the periodic cell border around new molecules",ttshown)
       call igPopItemWidth()

       ldum = iw_checkbox("Cubic cell",data%molcubic)
       call iw_tooltip("Read new molecules inside cubic periodic cell",ttshown)
       call igUnindent(0._c_float)
       call igNewLine()

       ! Input structure format (isformat)
       call iw_text("Read file format",highlight=.true.)
       data%isformat = isperm_inv_openfiles(data%isformat)
       call iw_combo_simple("##formatcombo",combostr_openfiles,data%isformat)
       call iw_tooltip("Force the new structure to be read with this file format, or auto-detect&
          &from the extension",ttshown)
       call igNewLine()
       data%isformat = isperm_openfiles(data%isformat)

    elseif (data%purpose == wpurp_dialog_openfieldfile) then
       ! Input field format (ifformat)
       call iw_text("Read file format",highlight=.true.)
       data%isformat = isperm_inv_openfieldfile(data%isformat)
       call iw_combo_simple("##formatcombo",combostr_openfieldfile,data%isformat)
       call iw_tooltip("Force the new field to be read with the given file format, or auto-detect&
          &from the extension",ttshown)
       call igNewLine()
       data%isformat = isperm_openfieldfile(data%isformat)

    elseif (data%purpose == wpurp_dialog_openvibfile) then
       ! Input field format (ifformat)
       call iw_text("Read vibration data format",highlight=.true.)
       data%isformat = isperm_inv_openvibfile(data%isformat)
       call iw_combo_simple("##formatcombo",combostr_openvibfile,data%isformat)
       call iw_tooltip("Force the new field to be read with the given file format, or auto-detect&
          &from the extension",ttshown)
       call igNewLine()
       data%isformat = isperm_openvibfile(data%isformat)
    end if

  end subroutine dialog_user_callback

end submodule proc
