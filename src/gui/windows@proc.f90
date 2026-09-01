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

  !> Return the ID of the view window that the view keybindings act
  !> upon: the focused view, or the parent view of a focused builder or
  !> dynamics window; if no window has the focus, the main view. If some
  !> other window is focused, return the main view, or zero if strict.
  module function view_target_window(strict)
    logical, intent(in), optional :: strict
    integer :: view_target_window

    integer :: i
    logical :: ok, strict_

    strict_ = .false.
    if (present(strict)) strict_ = strict

    view_target_window = iwin_view
    do i = 1, nwin
       if (.not.win(i)%isinit) cycle
       if (.not.win(i)%focused()) cycle
       if (win(i)%type == wintype_view) then
          view_target_window = i
       elseif (win(i)%type == wintype_builder .or. win(i)%type == wintype_dynamics) then
          ! these windows drive the view they were opened from
          ok = win(i)%idparent >= 1 .and. win(i)%idparent <= nwin
          if (ok) ok = win(win(i)%idparent)%isinit
          if (ok) ok = win(win(i)%idparent)%type == wintype_view
          if (ok) view_target_window = win(i)%idparent
       elseif (strict_) then
          view_target_window = 0
       end if
       exit
    end do

  end function view_target_window

  !> Initialize data required in the windows@proc submodule
  module subroutine windows_init()
    use fieldseedmod, only: field_format_string
    use param

    integer :: i, nn
    integer, allocatable :: itmp(:)
    character(len=:), allocatable :: fmtkey, fmtdesc
    logical :: isg

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
       // "Auto-detect" // c_null_char &             ! isformat_r_unknown
       // "Abinit DEN-style file" // c_null_char &   ! isformat_r_abinit
       // "Binary cube" // c_null_char &             ! isformat_r_bincube
       // "CASTEP cell" // c_null_char &             ! isformat_r_castepcell
       // "CASTEP geom" // c_null_char &             ! isformat_r_castepgeom
       // "CIF file" // c_null_char &                ! isformat_r_cif
       // "CRYSTAL output" // c_null_char &          ! isformat_r_crystal
       // "Cube" // c_null_char &                    ! isformat_r_cube
       // "DFTB+ gen file" // c_null_char &          ! isformat_r_gen
       // "DMACRYS .21 file" // c_null_char &        ! isformat_r_f21
       // "DMACRYS dmain file" // c_null_char &      ! isformat_r_dmain
       // "elk GEOMETRY.OUT" // c_null_char &        ! isformat_r_elk
       // "FHIaims input" // c_null_char &           ! isformat_r_aimsin
       // "FHIaims output" // c_null_char &          ! isformat_r_aimsout
       // "FPLO output" // c_null_char &             ! isformat_r_fploout
       // "Gaussian fchk" // c_null_char &           ! isformat_r_fchk
       // "Gaussian input" // c_null_char &          ! isformat_r_gjf
       // "Gaussian output" // c_null_char &         ! isformat_r_gaussian
       // "Gaussian wfn" // c_null_char &            ! isformat_r_wfn
       // "Gaussian wfx" // c_null_char &            ! isformat_r_wfx
       // "Gaussian zmat file" // c_null_char &      ! isformat_r_zmat
       // "mol2 file" // c_null_char &               ! isformat_r_mol2
       // "Molden-style file" // c_null_char &       ! isformat_r_molden
       // "ORCA output file" // c_null_char &        ! isformat_r_orca
       // "pdb file" // c_null_char &                ! isformat_r_pdb
       // "postg output" // c_null_char &            ! isformat_r_pgout
       // "psi4 output" // c_null_char &             ! isformat_r_dat
       // "Quantum ESPRESSO input" // c_null_char &  ! isformat_r_qein
       // "Quantum ESPRESSO output" // c_null_char & ! isformat_r_qeout
       // "Quantum ESPRESSO pwc" // c_null_char &    ! isformat_r_pwc
       // "SHELX" // c_null_char &                   ! isformat_r_shelx
       // "SIESTA IN/OUT file" // c_null_char &      ! isformat_r_siesta
       // "TINKER frac file" // c_null_char &        ! isformat_r_tinkerfrac
       // "VASP" // c_null_char &                    ! isformat_r_vasp
       // "WIEN2k struct file" // c_null_char &      ! isformat_r_struct
       // "Xcrysden axsf file" // c_null_char &      ! isformat_r_axsf
       // "Xcrysden xsf file" // c_null_char &       ! isformat_r_xsf
       // "xyz file" // c_null_char                  ! isformat_r_xyz
    nn = 0
    do i = 1,len(combostr_openfiles)
       if (combostr_openfiles(i:i) == c_null_char) nn = nn + 1
    end do

    !! Permutations for interpreting the format coming out of the dialog combo
    !! (same sequence as in combostr_openfiles)
    allocate(isperm_openfiles(0:nn-1))
    isperm_openfiles(0:nn-1) = (/isformat_r_unknown,isformat_r_abinit,isformat_r_bincube,&
       isformat_r_castepcell,isformat_r_castepgeom,isformat_r_cif,isformat_r_crystal,&
       isformat_r_cube,isformat_r_gen,isformat_r_f21,isformat_r_dmain,isformat_r_elk,&
       isformat_r_aimsin,isformat_r_aimsout,isformat_r_fploout,isformat_r_fchk,isformat_r_gjf,&
       isformat_r_gaussian,isformat_r_wfn,isformat_r_wfx,isformat_r_zmat,isformat_r_mol2,&
       isformat_r_molden,isformat_r_orca,isformat_r_pdb,isformat_r_pgout,isformat_r_dat,&
       isformat_r_qein,isformat_r_qeout,isformat_r_pwc,isformat_r_shelx,isformat_r_siesta,&
       isformat_r_tinkerfrac,isformat_r_vasp,isformat_r_struct,isformat_r_axsf,&
       isformat_r_xsf,isformat_r_xyz/)

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
                &DAT grid (DAT){.DAT},&
                &DFTB+ (detailed.xml){.xml},&
                &elk (grid){.grid},&
                &elk STATE.OUT (OUT){.OUT},&
                &FPLO grid (001){.001},&
                &Gaussian wavefunction (wfn|wfx|fchk){.wfn,.wfx,.fchk},&
                &Molden-style file (molden){.molden},&
                &Plain text grid (txt){.txt},&
                &Quantum ESPRESSO pwc file (pwc){.pwc},&
                &SIESTA (RHO...) {(RHO|BADER|DRHO|LDOS|VT|VH)},&
                &VASP ELFCAR{(ELFCAR)},&
                &VASP (CHGCAR...){(CONTCAR|CHGCAR|ELFCAR|CHG|AECCAR0|AECCAR1|AECCAR2|POSCAR)},&
                &WIEN2k (clmsum...){.clmsum,.clmup,.clmdn},&
                &Xcrysden (xsf|axsf) {.xsf,.axsf},&
                &"// c_null_char

    !! Formats offered in the dialog combo, in display order (Auto-detect
    !! first). This list is the only owner of the order; the display names
    !! come from field_format_string, so a new format needs only the
    !! fieldseedmod table, an entry here, and (optionally) the extension
    !! filter above
    itmp = (/ifformat_unknown,ifformat_abinit,ifformat_qub,ifformat_bincube,&
       ifformat_fmt,ifformat_cube,ifformat_dat,ifformat_dftb,ifformat_elkgrid,ifformat_elk,&
       ifformat_fplogrid,ifformat_wfn,ifformat_wfx,ifformat_fchk,ifformat_molden,ifformat_txt,&
       ifformat_pwc,ifformat_siestagrid,ifformat_vaspnov,ifformat_vasp,ifformat_wien,ifformat_xsf/)
    nn = size(itmp,1)
    allocate(isperm_openfieldfile(0:nn-1))
    isperm_openfieldfile(0:nn-1) = itmp
    deallocate(itmp)
    combostr_openfieldfile = "Auto-detect" // c_null_char
    do i = 1, nn-1
       call field_format_string(isperm_openfieldfile(i),fmtkey,fmtdesc,isg)
       combostr_openfieldfile = combostr_openfieldfile // fmtdesc // c_null_char
    end do

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
       &CASTEP phonon file (phonon){.phonon},&
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
       // "Gaussian fchk file" // c_null_char          & ! ivformat_gaussian_fchk
       // "CASTEP phonon file" // c_null_char            ! ivformat_castep_phonon
    nn = 0
    do i = 1,len(combostr_openvibfile)
       if (combostr_openvibfile(i:i) == c_null_char) nn = nn + 1
    end do

    allocate(isperm_openvibfile(0:nn-1))
    isperm_openvibfile(0:nn-1) = (/ivformat_unknown,ivformat_matdynmodes,&
       ivformat_matdyneig,ivformat_qedyn,ivformat_phonopy_ascii,ivformat_phonopy_yaml,&
       ivformat_phonopy_hdf5,ivformat_crystal_out,ivformat_gaussian_log,&
       ivformat_gaussian_fchk,ivformat_castep_phonon/)

    allocate(isperm_inv_openvibfile(0:maxval(isperm_openvibfile)))
    isperm_inv_openvibfile = 0
    do i = 0, nn-1
       isperm_inv_openvibfile(isperm_openvibfile(i)) = i
    end do

  end subroutine windows_init

  !xx! Module functions and subroutines

  !> Reallocate the window stack if there are only a few windows
  !> left. This must be done *outside* any window method because it
  !> will change the location of all window pointers.
  module subroutine stack_realloc_maybe()
    type(window), allocatable :: winaux(:)

    integer :: i

    integer, parameter :: iroom = 5

    if (nwin+iroom > size(win,1)) then
       ! the window array is deep-copied below, so drop the caches that
       ! can be rebuilt: an MO window can be holding hundreds of MB of
       ! sampling grids, and copying them to grow the stack is pure cost
       do i = 1, nwin
          call win(i)%drop_caches()
       end do
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
  module function stack_create_window(type,isopen,purpose,isys,irep,idparent,itoken,permanent,orraise,&
     dialog_filter)
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
    character(len=*,kind=c_char), intent(in), optional :: dialog_filter

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
             if (ok.and.type == wintype_dialog.and.present(purpose)) ok = (win(i)%purpose == purpose)
             if (ok.and.type == wintype_editrep.and.present(isys).and.present(irep).and.present(idparent)) then
                ok = (win(i)%isys == isys .and. win(i)%irep == irep .and. win(i)%idparent == idparent)
                if (ok.and.present(itoken)) &
                   ok = (win(i)%itoken == itoken)
             end if
             if (ok.and.type == wintype_scfplot.and.present(isys)) ok = (win(i)%isys == isys)
             if (ok.and.type == wintype_geometry.and.present(idparent)) ok = (win(i)%idparent == idparent)
             if (ok.and.type == wintype_vibrations.and.present(idparent)) ok = (win(i)%idparent == idparent)
             if (ok.and.type == wintype_mo.and.present(idparent)) ok = (win(i)%idparent == idparent)
             if (ok.and.type == wintype_dynamics.and.present(idparent)) ok = (win(i)%idparent == idparent)
             if (ok.and.type == wintype_builder.and.present(idparent)) ok = (win(i)%idparent == idparent)
             if (ok.and.type == wintype_exportimage.and.present(idparent)) ok = (win(i)%idparent == idparent)
             if (ok.and.type == wintype_saveas.and.present(idparent)) ok = (win(i)%idparent == idparent)
             if (ok.and.type == wintype_save_multiple.and.present(idparent)) ok = (win(i)%idparent == idparent)
             if (ok.and.type == wintype_extract.and.present(idparent)) ok = (win(i)%idparent == idparent)
             if (ok.and.type == wintype_rattle.and.present(idparent)) ok = (win(i)%idparent == idparent)
             if (ok.and.type == wintype_water_cluster.and.present(idparent)) ok = (win(i)%idparent == idparent)
             if (ok.and.type == wintype_load_field.and.present(isys)) ok = (win(i)%isys == isys)
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

    ! allocate the stack if this is the first window
    if (.not.allocated(win)) allocate(win(maxwin))

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

    ! The stack is grown by stack_realloc_maybe once per frame; it cannot be
    ! grown here because that would move all the window pointers (see the
    ! comment in stack_realloc_maybe). Overflowing it means more windows were
    ! created in a single frame than the slack allows: a bug, not a user error.
    if (nwin > size(win,1)) &
       call ferror('stack_create_window','too many windows',faterr)

    ! initialize the new window
    call win(id)%init(type,isopen,id,purpose,isys,irep,idparent,itoken,dialog_filter)
    if (present(permanent)) then
       win(id)%permanent = permanent
    else
       win(id)%permanent = .false.
    end if
    stack_create_window = id

  end function stack_create_window

  !> Show a danger-colored warning if the window's okfile exists on
  !> disk. The existence check is cached and re-run only when okfile
  !> changes.
  module subroutine okfile_warn_overwrite(w)
    use utils, only: iw_text
    class(window), intent(inout) :: w

    if (w%okfile /= w%okfile_lastcheck) then
       w%okfile_lastcheck = w%okfile
       if (len_trim(w%okfile) > 0) then
          inquire(file=w%okfile,exist=w%okfile_exists)
       else
          w%okfile_exists = .false.
       end if
    end if
    if (w%okfile_exists) &
       call iw_text("File exists and will be overwritten",danger=.true.,wrap=.true.)

  end subroutine okfile_warn_overwrite

  !> Default file name for a save window operating on system isys: the
  !> root of the system's source file (defname if unavailable) plus
  !> extension ext. Unless uselastdir is present and false, put the
  !> file in the last directory recorded by okfile_save_dir (if any)
  !> instead of the source file's directory.
  module function okfile_default(isys,defname,ext,uselastdir) result(file)
    use systems, only: sysc, ok_system, sys_loaded_not_init
    use utils, only: get_current_working_dir, file_name_root
    use param, only: dirsep
    integer, intent(in) :: isys
    character(len=*), intent(in) :: defname
    character(len=*), intent(in) :: ext
    logical, intent(in), optional :: uselastdir
    character(kind=c_char,len=:), allocatable :: file

    character(len=:), allocatable :: root
    integer :: idx
    logical :: uselast

    uselast = .true.
    if (present(uselastdir)) uselast = uselastdir
    if (.not.allocated(okfile_lastdir)) okfile_lastdir = ""

    ! root of the system's source file, or defname in the cwd
    root = ""
    if (ok_system(isys,sys_loaded_not_init)) then
       if (len_trim(sysc(isys)%seed%file) > 0) &
          root = file_name_root(sysc(isys)%seed%file)
    end if
    if (len_trim(root) == 0) &
       root = get_current_working_dir() // dirsep // defname

    ! replace the directory with the last-used one
    if (uselast .and. len_trim(okfile_lastdir) > 0) then
       idx = index(root,dirsep,back=.true.)
       root = okfile_lastdir // dirsep // root(idx+1:)
    end if

    file = root // "." // trim(ext)

  end function okfile_default

  !> Record the directory of file as the last one used by the save
  !> windows (used by okfile_default). No-op if file has no path.
  module subroutine okfile_save_dir(file)
    use param, only: dirsep
    character(len=*), intent(in) :: file

    integer :: idx

    idx = index(file,dirsep,back=.true.)
    if (idx > 1) okfile_lastdir = file(1:idx-1)

  end subroutine okfile_save_dir

  !> Build the option string for the write-format combos, unless it
  !> has been built already. The first entry is "Auto-detect" (used
  !> only by the save-as window) and the rest are the write formats in
  !> fmtperm order, starting at write_format_combostr(icombo_fmt1:).
  module subroutine build_write_format_combo()

    integer :: i

    if (allocated(write_format_combostr)) return
    write_format_combostr = "Auto-detect" // c_null_char
    icombo_fmt1 = len(write_format_combostr) + 1
    do i = 1, isformat_w_max
       write_format_combostr = write_format_combostr // trim(fmtnames(fmtperm(i))) // c_null_char
    end do

  end subroutine build_write_format_combo

  !> This routine regenerates all pointers to the widows in the win(:)
  !> structure and its components. It is used when an array size is
  !> exceeded and move_alloc needs to be used to allocate more memory.
  module subroutine regenerate_window_pointers()
    use systems, only: sysc, sys_init, ok_system

    integer :: i, iv

    do i = 1, nwin
       if (win(i)%type == wintype_editrep .and. win(i)%irep > 0) then
          win(i)%rep => win(win(i)%idparent)%sc%rep(win(i)%irep)
       elseif (win(i)%type == wintype_view) then
          win(i)%sc => null()
          iv = win(i)%isys
          if (ok_system(iv,sys_init)) then
             if (win(i)%ismain) win(i)%sc => sysc(iv)%sc
          end if
       end if
    end do

  end subroutine regenerate_window_pointers

  !> Initialize a window of the given type. If isiopen, initialize it
  !> as open.
  module subroutine window_init(w,type,isopen,id,purpose,isys,irep,idparent,itoken,dialog_filter)
    use interfaces_opengl3
    use gui_main, only: ColorDialogDir, ColorDialogFile
    use systems, only: always_read_virtuals
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
    character(len=*,kind=c_char), intent(in), optional :: dialog_filter

    character(kind=c_char,len=:), allocatable, target :: str1

    ! initialization of the state
    w%isinit = .true.
    w%firstpass = .true.
    w%isopen = isopen
    w%type = type
    w%id = -id
    w%name = "" // c_null_char
    w%errmsg = ""
    w%okmsg = ""
    w%isys = 1
    w%growtofit = .false.
    w%needheight = 0._c_float
    w%needwidth = 0._c_float
    w%sortcid = 0
    w%sortdir = 1
    w%okfile = ""
    ! only the one-file modal dialog (wpurp_dialog_openonefilemodal)
    ! honors dialog_filter; other purposes use their own filter lists
    if (allocated(w%dialog_filter)) deallocate(w%dialog_filter)
    if (present(dialog_filter)) w%dialog_filter = dialog_filter
    w%okfile_lastcheck = ""
    w%okfile_set = .false. ! whether the library file has been set by the user
    w%okfile_read = .false. ! whether the structure list should be re-read from the lib
    nullify(w%sc)
    w%forcerender = .true.
    if (allocated(w%iord)) deallocate(w%iord)
    w%lastselected = 0
    w%tabselected = ""
    w%dialog_data%dptr = c_null_ptr
    w%dialog_data%mol = -1
    w%dialog_data%showhidden = .false._c_bool
    w%dialog_data%isformat = isformat_r_unknown
    w%dialog_data%readlastonly = .false._c_bool
    w%dialog_data%purpose = wpurp_unknown
    w%dialog_data%molcubic = .false.
    w%dialog_data%rborder = rborder_def*bohrtoa
    w%plotn = 0
    if (present(isys)) w%isys = isys
    if (present(irep)) w%irep = irep
    if (present(idparent)) w%idparent = idparent
    if (present(itoken)) w%itoken = itoken
    if (present(purpose)) w%purpose = purpose
    w%geometry_expression = ""

    ! Type-specific initialization
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
       w%lf = loadfield_state() ! the slot may carry a previous window's form
       w%lf%readvirtual = always_read_virtuals
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
    elseif (type == wintype_saveas) then
       ! save structure as window
       if (.not.present(idparent)) &
          call ferror('window_init','saveas requires idparent',faterr)
       if (allocated(w%saveas_lastcontrol)) deallocate(w%saveas_lastcontrol)
       w%saveas_control_exists = .false.
    elseif (type == wintype_save_multiple) then
       ! save multiple structures window
       if (.not.present(idparent)) &
          call ferror('window_init','save_multiple requires idparent',faterr)
       w%sm = savemult_state() ! the slot may carry a previous window's list
       ! a collapsed window consumes firstpass without drawing its body, so
       ! these have to be usable before that block ever runs
       w%savemult_pattern = ""
       w%savemult_scope = 0
       w%savemult_format = 1
       w%savemult_exist = 0
       w%savemult_rk = 50._c_float
       w%savemult_nosym = .false.
       w%savemult_cartesian = .false.
       w%savemult_docell = .true.
       if (allocated(w%savemult_lastcheck)) deallocate(w%savemult_lastcheck)
       w%savemult_direxists = .false.
    elseif (type == wintype_extract) then
       ! extract cluster as molecule window
       if (.not.present(idparent)) &
          call ferror('window_init','extract requires idparent',faterr)
    elseif (type == wintype_rattle) then
       ! rattle structure window
       if (.not.present(idparent)) &
          call ferror('window_init','rattle requires idparent',faterr)
    elseif (type == wintype_vibrations) then
       ! vibrations window
       if (.not.present(idparent)) &
          call ferror('window_init','vibrations requires idparent',faterr)
    elseif (type == wintype_mo) then
       ! molecular orbitals window
       if (.not.present(idparent)) &
          call ferror('window_init','mo requires idparent',faterr)
    elseif (type == wintype_dynamics) then
       ! dynamics window
       if (.not.present(idparent)) &
          call ferror('window_init','dynamics requires idparent',faterr)
       w%isys = 0
    elseif (type == wintype_water_cluster) then
       ! water cluster demonstration window
       if (.not.present(idparent)) &
          call ferror('window_init','water_cluster requires idparent',faterr)
       w%isys = 0
       ! per-participant state: never inherit the previous player's name
       ! from a recycled window slot
       w%wc_name = ""
    elseif (type == wintype_builder) then
       ! builder window, tied to the parent view
       if (.not.present(idparent)) &
          call ferror('window_init','builder requires idparent',faterr)
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
    elseif (type == wintype_geometry) then
       ! geometry window: anchored to the view whose system it edits
       if (.not.present(idparent)) &
          call ferror('window_init','geometry requires idparent',faterr)
    end if

  end subroutine window_init

  !> Release the caches this window can rebuild on demand, without
  !> otherwise changing what it shows. Used before the window stack is
  !> deep-copied to grow it, where carrying the caches along is pure
  !> cost.
  module subroutine window_drop_caches(w)
    class(window), intent(inout) :: w

    if (w%type == wintype_mo) then
       ! the per-orbital sampling grids and the diagram's level list,
       ! both rebuilt on demand
       w%mo_cache = mo_cache_state()
       w%mo_diag = mo_diagram_state()
    end if

  end subroutine window_drop_caches

  !> End a window and deallocate the data.
  module subroutine window_end(w)
    use systems, only: ok_system, sysc, sys, sys_init, lastchange_geometry, remove_system
    use interfaces_opengl3
    use tools_io, only: ferror, warning
    class(window), intent(inout), target :: w

    integer :: isysd
    character(len=:), allocatable :: errmsg

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
       elseif (w%type == wintype_mo) then
          ! release the cached orbital sampling grids and the level list
          w%mo_cache = mo_cache_state()
          w%mo_diag = mo_diagram_state()
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
       elseif (w%type == wintype_dynamics) then
          ! stop the dynamics run, rebuild the crystal, free the state
          isysd = w%isys
          if (ok_system(isysd,sys_init)) then
             if (sysc(isysd)%md%ready) then
                ! the window is going away, so its errmsg cannot be shown:
                ! report a failed rebuild in the output console instead
                call sys(isysd)%c%rebuild_after_move(errmsg=errmsg)
                if (len_trim(errmsg) > 0) &
                   call ferror('window_end','could not rebuild the structure after the '//&
                   'dynamics run: '//trim(errmsg),warning)
                call sysc(isysd)%post_event(lastchange_geometry)
             end if
             sysc(isysd)%md_run = .false.
             call sysc(isysd)%md%free()
             if (w%idparent > 0 .and. w%idparent <= nwin) then
                if (associated(win(w%idparent)%sc)) &
                   win(w%idparent)%forcerender = .true.
             end if
          end if
       elseif (w%type == wintype_rattle) then
          ! a sampling run belongs to this window: closing it must
          ! stop the dynamics and put the sampled-from structure back
          if (w%rattle_running) then
             isysd = w%rattle_isys
             w%rattle_running = .false.
             if (ok_system(isysd,sys_init)) then
                if (sysc(isysd)%md%ready) call sysc(isysd)%md%reset(sys(isysd)%c)
                call sysc(isysd)%md_stop(errmsg)
                if (len_trim(errmsg) > 0) &
                   call ferror('window_end','could not rebuild the structure after the '//&
                   'sampling run: '//trim(errmsg),warning)
             end if
          end if
          if (allocated(w%rattle_seed)) deallocate(w%rattle_seed)
       elseif (w%type == wintype_save_multiple) then
          if (allocated(w%savemult_pattern)) deallocate(w%savemult_pattern)
          w%sm = savemult_state() ! the file list is this window's, not the session's
       elseif (w%type == wintype_water_cluster) then
          ! the demo owns its generated cluster; remove it on close so it does not linger
          isysd = w%isys
          if (ok_system(isysd,sys_init)) then
             sysc(isysd)%md_run = .false.
             if (sysc(isysd)%md%ready) call sysc(isysd)%md%free()
             call remove_system(isysd)
          end if
       elseif (w%type == wintype_builder) then
          ! release the forced builder mode on the parent view, if still active
          if (w%builder_vm /= 0 .and. w%idparent >= 1 .and. w%idparent <= nwin) then
             if (win(w%idparent)%isinit) &
                call win(w%idparent)%viewmode_release_forced(w%id,w%builder_vm)
          end if
          ! stop the edit session
          call w%edit_stop()
       elseif (w%type == wintype_geometry) then
          ! remove all highlights
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
    w%ptr = c_null_ptr
    w%timelast_focused = 0d0
    w%name = "" // c_null_char
    if (allocated(w%iord)) deallocate(w%iord)
    w%lastselected = 0
    w%tabselected = ""
    w%plotn = 0
    if (allocated(w%plotx)) deallocate(w%plotx)
    if (allocated(w%ploty)) deallocate(w%ploty)
    nullify(w%rep)
    w%irep = 0
    w%viewmode = vm_navigate
    w%viewmode_transient = .false.
    w%vmdata = viewmode_data() ! drop the stale forced-mode owner and pick data
    w%geometry_expression = ""
    w%geometry_expression_ok = .false.
    w%geometry_expr_error = ""
    call w%geometry_addbond%clear()
    w%editrep_pick_item = 0
    w%editrep_pick_slot = 0
    w%editrep_isopick = -1
    w%editrep_isoline = 0
    call w%editrep_pick%clear()
    w%wc_started = .false.
    w%builder_tool = it_none
    w%builder_vm = 0
    w%builder_isys = 0
    ! (edit-session markers are cleared by the edit_stop call above, and
    ! the session payload is re-seeded by the start routines)
    w%edit_pending = .false.
    w%edit_time = 0d0
    if (allocated(w%edit_frag)) deallocate(w%edit_frag)

  end subroutine window_end

  !> Draw the energy-backend combo for system isys.
  module subroutine draw_ff_backend_combo(isys,strid,nchars)
    use systems, only: sys, sysc
    use energy, only: ff_list, ff_backend_applicable, ff_backend_label, ff_backend_default
    use utils, only: iw_combo_simple, iw_calcwidth
    integer, intent(in) :: isys, nchars
    character(len=*), intent(in) :: strid

    integer :: i, nback, icombo, backids(size(ff_list))
    character(len=:), allocatable :: str

    ! resolve the automatic default and inapplicable backends
    if (sysc(isys)%md_backend < 0 .or.&
       .not.ff_backend_applicable(sysc(isys)%md_backend,sys(isys)%c)) &
       sysc(isys)%md_backend = ff_backend_default(sys(isys)%c)

    ! build the list of applicable backends
    str = ""
    nback = 0
    do i = 1, size(ff_list)
       if (.not.ff_backend_applicable(ff_list(i),sys(isys)%c)) cycle
       nback = nback + 1
       backids(nback) = ff_list(i)
       str = str // trim(ff_backend_label(ff_list(i))) // c_null_char
    end do

    ! translate the stored backend id to its position in the filtered list
    icombo = 0
    do i = 1, nback
       if (backids(i) == sysc(isys)%md_backend) icombo = i - 1
    end do
    call igSameLine(0._c_float,-1._c_float)
    call igPushItemWidth(iw_calcwidth(nchars,1))
    call iw_combo_simple(strid,str,icombo)
    call igPopItemWidth()

    ! persist the selection
    sysc(isys)%md_backend = backids(icombo+1)

  end subroutine draw_ff_backend_combo

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

  !> Index in of the live view window this window is anchored to,
  !> taken from w%idparent, or zero if the anchor is gone or was never
  !> a view window.
  module function window_anchor_view(w) result(iview)
    class(window), intent(in) :: w
    integer :: iview

    iview = 0
    if (w%idparent < 1 .or. w%idparent > nwin) return
    if (.not.win(w%idparent)%isinit) return
    if (win(w%idparent)%type /= wintype_view) return
    iview = w%idparent

  end function window_anchor_view

  !> Resolve this window's anchor: the view window it is bound to
  !> (w%idparent) and the system that view is currently
  !> showing. Returns .false. if the anchor is gone or is no longer a
  !> view window, in which case the caller should close
  !> itself. changed is .true. on the frame the anchor view switched
  !> to a different system, so that the caller can drop the state it
  !> had cached for the previous one; it also covers a window that was
  !> raised onto a different view by stack_create_window. The resolved
  !> system is cached in w%isys.
  module function window_anchor(w,iview,isys,changed) result(ok)
    class(window), intent(inout) :: w
    integer, intent(out) :: iview
    integer, intent(out) :: isys
    logical, intent(out) :: changed
    logical :: ok

    iview = w%anchor_view()
    ok = (iview > 0)
    isys = 0
    changed = .false.
    if (.not.ok) return

    isys = win(iview)%isys
    changed = (isys /= w%isys)
    w%isys = isys

  end function window_anchor

  !> Point this window's anchor view at system isys. Does nothing if
  !> the anchor is gone.
  module subroutine window_retarget(w,isys)
    class(window), intent(inout) :: w
    integer, intent(in) :: isys

    integer :: iview

    iview = w%anchor_view()
    if (iview > 0) &
       call win(iview)%select_view(isys)

  end subroutine window_retarget

  !> Draw an ImGui window.
  module subroutine window_draw(w)
    use gui_main, only: fontsize, io, g
    use interfaces_glfw, only: glfwGetTime
    use utils, only: iw_text, get_nice_next_window_pos, iw_calcwidth
    use tools_io, only: string, ferror, faterr
    class(window), intent(inout), target :: w

    character(kind=c_char,len=:), allocatable, target :: str1, str2, str3
    real(c_float) :: panewidth, hneed
    type(ImVec2) :: inisize, pos, pivot, szmin, szmax, szwant
    type(ImGuiWindow), pointer :: wptr

    if (.not.w%isinit) return
    if (.not.w%isopen.and..not.w%permanent) then
       ! the window may have been closed in imgui (x button)
       call w%end()
       return
    end if

    ! record the last time this window was focused (close-dialog bind)
    if (w%focused()) w%timelast_focused = glfwGetTime()

    ! First pass on creation: assign ID, name, and flags
    if (w%id < 0) then
       w%firstpass = .true.
       w%id = abs(w%id)
       w%pos = 0._c_float
       w%isdocked = .false.
       if (w%type == wintype_tree) then
          call init_window("Tree")
       elseif (w%type == wintype_view) then
          if (w%ismain) then
             call init_window("Main View")
          else
             call init_window("Alternate View",90,30)
             call get_nice_next_window_pos(pos)
             pivot%x = 0._c_float
             pivot%y = 0._c_float
             call igSetNextWindowPos(pos,ImGuiCond_None,pivot)
          end if
       elseif (w%type == wintype_console_input) then
          call init_window("Input")
       elseif (w%type == wintype_console_output) then
          call init_window("Output")
       elseif (w%type == wintype_about) then
          call init_window("About",52)
       elseif (w%type == wintype_dialog) then
          w%dialog_data%dptr = w%dptr
          w%dialog_data%purpose = w%purpose
          if (w%dialog_data%showhidden) then
             w%flags = ImGuiFileDialogFlags_None
          else
             w%flags = ImGuiFileDialogFlags_DontShowHiddenFiles
          end if
          str2 = "" // c_null_char ! default path
          inisize%x = 90 * fontsize%x
          inisize%y = 30 * fontsize%y
          call clamp_to_display(inisize)
          call igSetNextWindowSize(inisize,ImGuiCond_FirstUseEver)
          ! a file dialog is the widest window in the GUI: center it, or the
          ! cascade offset pushes its side pane past the edge of the display
          pos%x = max(0.5_c_float * (io%DisplaySize%x - inisize%x),0._c_float)
          pos%y = max(0.5_c_float * (io%DisplaySize%y - inisize%y),0._c_float)
          pivot%x = 0._c_float
          pivot%y = 0._c_float
          call igSetNextWindowPos(pos,ImGuiCond_FirstUseEver,pivot)

          ! width of the dialog side pane: wide enough for its longest
          ! option label at the current font size, so that raising the UI
          ! scale does not clip the pane contents
          panewidth = max(280._c_float,iw_calcwidth(28,0,1))

          str1 = "All files (*.*){*.*}" // c_null_char
          if (w%purpose == wpurp_dialog_openfiles) then
             ! open dialog
             w%name = "Open File(s)##" // string(w%id)  // c_null_char
             call IGFD_OpenPaneDialog2(w%dptr,c_loc(w%name),c_loc(w%name),c_loc(dialogstr_openfiles),c_loc(str2),&
                c_funloc(dialog_user_callback),panewidth,0_c_int,c_loc(w%dialog_data),w%flags)
          elseif (w%purpose == wpurp_dialog_savelogfile) then
             w%name = "Save Log File##" // string(w%id)  // c_null_char
             str2 = "file.log" // c_null_char
             str3 = "./" // c_null_char
             call IGFD_OpenPaneDialog(w%dptr,c_loc(w%name),c_loc(w%name),c_loc(str1),c_loc(str3),c_loc(str2),&
                c_funloc(dialog_user_callback),panewidth,1_c_int,c_loc(w%dialog_data),w%flags)
          elseif (w%purpose == wpurp_dialog_openlibraryfile) then
             w%name = "Open Library File##" // string(w%id)  // c_null_char
             call IGFD_OpenPaneDialog2(w%dptr,c_loc(w%name),c_loc(w%name),c_loc(str1),c_loc(str2),&
                c_funloc(dialog_user_callback),panewidth,1_c_int,c_loc(w%dialog_data),w%flags)
          elseif (w%purpose == wpurp_dialog_openfieldfile) then
             w%name = "Open Field File(s)##" // string(w%id)  // c_null_char
             call IGFD_OpenPaneDialog2(w%dptr,c_loc(w%name),c_loc(w%name),c_loc(dialogstr_openfieldfile),&
                c_loc(str2),c_funloc(dialog_user_callback),panewidth,1_c_int,c_loc(w%dialog_data),w%flags)
          elseif (w%purpose == wpurp_dialog_openvibfile) then
             w%name = "Open Vibration Data File(s)##" // string(w%id)  // c_null_char
             call IGFD_OpenPaneDialog2(w%dptr,c_loc(w%name),c_loc(w%name),c_loc(dialogstr_openvibfile),&
                c_loc(str2),c_funloc(dialog_user_callback),panewidth,0_c_int,c_loc(w%dialog_data),w%flags)
          elseif (w%purpose == wpurp_dialog_openonefilemodal) then
             w%name = "Open File##" // string(w%id)  // c_null_char
             if (allocated(w%dialog_filter)) &
                str1 = w%dialog_filter // ",All files (*.*){*.*}" // c_null_char
             call IGFD_OpenPaneDialog2(w%dptr,c_loc(w%name),c_loc(w%name),c_loc(str1),c_loc(str2),&
                c_funloc(dialog_user_callback),panewidth,1_c_int,c_loc(w%dialog_data),w%flags)
          elseif (w%purpose == wpurp_dialog_saveimagefile) then
             w%name = "Save Image File##" // string(w%id) // c_null_char
             str1 = "PNG (*.png) {.png},BMP (*.bmp) {.bmp},TGA (*.tga) {.tga},JPEG (*.jpg) {.jpg}"// c_null_char
             call dialog_initial_file(w%idparent,"image.png",str2,str3)
             call IGFD_OpenPaneDialog(w%dptr,c_loc(w%name),c_loc(w%name),c_loc(str1),c_loc(str3),c_loc(str2),&
                c_funloc(dialog_user_callback),panewidth,1_c_int,c_loc(w%dialog_data),w%flags)
          elseif (w%purpose == wpurp_dialog_savefile) then
             w%name = "Save Structure File##" // string(w%id) // c_null_char
             str1 = "&
                &FHIaims input (*.in) {.in},&
                &CIF (*.cif) {.cif},&
                &Quantum ESPRESSO input (*.pwi) {.pwi},&
                &VASP (*.POSCAR) {.POSCAR},&
                &CASTEP cell (*.cell) {.cell},&
                &SHELX res (*.res) {.res},&
                &CRYSTAL input (*.d12) {.d12},&
                &abinit input (*.abin) {.abin},&
                &elk GEOMETRY.OUT (*.elk) {.elk},&
                &Gaussian input (*.gjf) {.gjf},&
                &Gaussian input, periodic (*.gau) {.gau},&
                &xyz (*.xyz) {.xyz},&
                &CML (*.cml) {.cml},&
                &critic2 input (*.cri) {.cri},&
                &escher/octave (*.m) {.m},&
                &dcp database (*.db) {.db},&
                &GULP input (*.gin) {.gin},&
                &LAMMPS data (*.lammps) {.lammps},&
                &SIESTA fdf (*.fdf) {.fdf},&
                &SIESTA STRUCT_IN (*.struct_in) {.struct_in},&
                &DFTB+ hsd (*.hsd) {.hsd},&
                &DFTB+ gen (*.gen) {.gen},&
                &pyscf script (*.pyscf) {.pyscf},&
                &TINKER frac (*.frac) {.frac},&
                &tessel (*.tess) {.tess},&
                &pdb (*.pdb) {.pdb},&
                &Wavefront obj (*.obj) {.obj},&
                &PLY (*.ply) {.ply},&
                &OFF (*.off) {.off},&
                &All files (*.*){*.*}"// c_null_char
             call dialog_initial_file(w%idparent,"structure.in",str2,str3)
             ! the overwrite flag goes only to the file dialog, not to w%flags: it
             ! would be reused as ImGuiWindowFlags in IGFD_DisplayDialog
             call IGFD_OpenPaneDialog(w%dptr,c_loc(w%name),c_loc(w%name),c_loc(str1),c_loc(str3),c_loc(str2),&
                c_funloc(dialog_user_callback),panewidth,1_c_int,c_loc(w%dialog_data),&
                ior(w%flags,ImGuiFileDialogFlags_ConfirmOverwrite))
          elseif (w%purpose == wpurp_dialog_selectdir) then
             w%name = "Select Directory##" // string(w%id) // c_null_char
             str2 = "" // c_null_char
             str3 = "./" // c_null_char
             if (w%idparent >= 1 .and. w%idparent <= nwin) then
                if (allocated(win(w%idparent)%okfile)) then
                   if (len_trim(win(w%idparent)%okfile) > 0) &
                      str3 = trim(win(w%idparent)%okfile) // c_null_char
                end if
             end if
             ! a null filter list puts the dialog in directory-selection mode
             call IGFD_OpenPaneDialog(w%dptr,c_loc(w%name),c_loc(w%name),c_null_ptr,c_loc(str3),c_loc(str2),&
                c_funloc(dialog_user_callback),panewidth,1_c_int,c_loc(w%dialog_data),w%flags)
          else
             call ferror('window_draw','unknown dialog purpose: ' // string(w%purpose),faterr)
          end if
       elseif (w%type == wintype_new_struct) then
          call init_window("New Structure",90,40)
       elseif (w%type == wintype_new_struct_library) then
          call init_window("New Structure from Library",90,30)
       elseif (w%type == wintype_load_field) then
          call init_window("Load Field",65)
       elseif (w%type == wintype_scfplot) then
          call init_window("SCF Iterations (" // string(w%isys) // ")",45,45,square=.true.)
       elseif (w%type == wintype_editrep) then
          call init_window("Object [" // string(w%rep%name) // "]",63,46)
       elseif (w%type == wintype_exportimage) then
          call init_window("Export to Image",52)
       elseif (w%type == wintype_saveas) then
          call init_window("Save As",52)
       elseif (w%type == wintype_save_multiple) then
          call init_window("Save Multiple Structures",60,36)
       elseif (w%type == wintype_extract) then
          call init_window("Extract as Molecule(s)",52)
       elseif (w%type == wintype_rattle) then
          call init_window("Rattle Structure",52)
       elseif (w%type == wintype_vibrations) then
          call init_window("Vibrations",62)
       elseif (w%type == wintype_mo) then
          ! a starting point only: draw_mo asks for a width that fits
          ! the diagram and however many spin tables this wavefunction
          ! turns out to need, as soon as it knows
          call init_window("Molecular Orbitals",70)
       elseif (w%type == wintype_dynamics) then
          call init_window("Dynamics",55)
       elseif (w%type == wintype_water_cluster) then
          call init_window("Water Cluster Demonstration",55,30)
       elseif (w%type == wintype_geometry) then
          call init_window("View/Edit Geometry",70,30)
       elseif (w%type == wintype_preferences) then
          call init_window("Preferences",81,33)
       elseif (w%type == wintype_treeplot) then
          call init_window("Plot Tree Data",55,25)
       elseif (w%type == wintype_builder) then
          ! wide enough for the tool palette on a single row
          call init_window("Builder",63,30)
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
       elseif (w%type == wintype_console_output) then
          call w%update_co()
       end if
    end if

    ! Draw the window contents, depending on type
    if (w%isopen) then
       ! assign the pointer ID for the window, if not a dialog
       if (w%type == wintype_dialog) then
          w%ptr = IGFD_GetCurrentWindow(w%dptr)
          call w%draw_dialog()
       else
          ! keep floating windows within a usable size range: at least a
          ! few glyphs, at most the display (docked windows follow their
          ! dock node and are left alone). In a grow-to-fit window whose
          ! content overflowed last frame, raise the minimum to the
          ! height of the content, which grows the window to show all of
          ! it; manual resizing works otherwise, so the window can never
          ! stay shorter than its content but keeps its size when the
          ! content shrinks. Content taller than the display grows the
          ! window to the display and scrolls for the rest: giving up on
          ! the grow instead would leave the window at whatever size it
          ! happened to have, which on a small display is the bare
          ! minimum -- a couple of header lines and a scroll bar.
          if (.not.w%isdocked) then
             szmax%x = huge(1._c_float)
             szmax%y = huge(1._c_float)
             call clamp_to_display(szmax)
             szmin%x = min(25._c_float * fontsize%x,szmax%x)
             szmin%y = min(6._c_float * fontsize%y,szmax%y)
             if (w%needheight > 0._c_float) then
                hneed = min(w%needheight,szmax%y)
                szmin%y = max(szmin%y,hneed)
                ! the buttons are anchored to the window's bottom: raise
                ! the window if the grown window would stick out under
                ! the bottom edge of the display
                call c_f_pointer(w%ptr,wptr)
                if (wptr%Pos%y + hneed > szmax%y) then
                   pos%x = wptr%Pos%x
                   pos%y = max(szmax%y - hneed,0._c_float)
                   pivot%x = 0._c_float
                   pivot%y = 0._c_float
                   call igSetNextWindowPos(pos,ImGuiCond_Always,pivot)
                end if
             end if
             call igSetNextWindowSizeConstraints(szmin,szmax,c_null_funptr,c_null_ptr)

             ! a one-shot width, asked for by the window's own content
             ! once it knows how wide it wants to be (the initial width
             ! passed to init_window is fixed before the content is
             ! known). Applied once and cleared, so the user can resize
             ! the window afterwards and keep the size.
             if (w%needwidth > 0._c_float .and. c_associated(w%ptr)) then
                call c_f_pointer(w%ptr,wptr)
                szwant%x = min(max(w%needwidth,szmin%x),szmax%x)
                szwant%y = wptr%Size%y
                call igSetNextWindowSize(szwant,ImGuiCond_Always)
             end if
          end if
          ! cleared whether or not it could be applied: a docked window
          ! must not keep the request and spring to that width the first
          ! time it is undocked, long after the user has sized it
          w%needwidth = 0._c_float
          if (igBegin(c_loc(w%name),w%isopen,w%flags)) then
             w%ptr = igGetCurrentWindow()
             w%heightslack = 0._c_float ! only a body that stretches to fill sets this
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
                call w%draw_new_struct_library()
             elseif (w%type == wintype_load_field) then
                call w%draw_load_field()
             elseif (w%type == wintype_scfplot) then
                call w%draw_scfplot()
             elseif (w%type == wintype_editrep) then
                call w%draw_editrep()
             elseif (w%type == wintype_exportimage) then
                call w%draw_exportimage()
             elseif (w%type == wintype_saveas) then
                call w%draw_saveas()
             elseif (w%type == wintype_save_multiple) then
                call w%draw_save_multiple()
             elseif (w%type == wintype_extract) then
                call w%draw_extract()
             elseif (w%type == wintype_rattle) then
                call w%draw_rattle()
             elseif (w%type == wintype_vibrations) then
                call w%draw_vibrations()
             elseif (w%type == wintype_mo) then
                call w%draw_mo()
             elseif (w%type == wintype_dynamics) then
                call w%draw_dynamics()
             elseif (w%type == wintype_water_cluster) then
                call w%draw_water_cluster()
             elseif (w%type == wintype_geometry) then
                call w%draw_geometry()
             elseif (w%type == wintype_preferences) then
                call w%draw_preferences()
             elseif (w%type == wintype_treeplot) then
                call w%draw_treeplot()
             elseif (w%type == wintype_builder) then
                call w%draw_builder()
             end if

             ! The window height that shows all of this window's content,
             ! for the grow-to-fit constraint of the next frame. Measured
             ! from the laid-out content rather than from the scrollbar:
             ! a scrollbar only exists while the window is too short, so
             ! deriving the height from it makes the constraint switch
             ! itself off as soon as it has worked, and the window then
             ! flips between the two sizes on alternate frames.
             ! A body that stretches part of itself to fill the window
             ! reports how much it added: measured with the stretch in it,
             ! the content would grow with the window, the minimum would
             ! ratchet up to whatever height the user last set, and the
             ! window could never be made smaller again
             if (w%growtofit .and. .not.w%isdocked) &
                w%needheight = igGetCursorPosY() + g%Style%WindowPadding%y - w%heightslack
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

  contains
    !> Name the window being created after title (plus the window ID, which
    !> keeps the ImGui identifier unique), give it the default flags and, if nx
    !> and ny are present, an initial size of nx by ny characters. If only nx
    !> is present, fit the initial height to the first frame's content and grow
    !> the window whenever its content overflows it (grow-to-fit).
    subroutine init_window(title,nx,ny,square)
      character(len=*,kind=c_char), intent(in) :: title
      integer, intent(in), optional :: nx, ny
      logical, intent(in), optional :: square

      logical :: square_

      square_ = .false.
      if (present(square)) square_ = square

      w%name = title // "##" // string(w%id) // c_null_char
      w%flags = ImGuiWindowFlags_None
      if (present(nx)) then
         inisize%x = nx * fontsize%x
         if (present(ny)) then
            ! square windows measure the height in glyph widths too, so that
            ! they come out actually square and not one line-height per
            ! character
            if (square_) then
               inisize%y = ny * fontsize%x
            else
               inisize%y = ny * fontsize%y
            end if
         else
            ! no height given: fit the height to the content on the first
            ! frame (a zero axis means auto-fit) and grow the window when
            ! its content overflows it (see the size constraints before
            ! igBegin)
            w%growtofit = .true.
            inisize%y = 0._c_float
         end if
         call clamp_to_display(inisize)
         call igSetNextWindowSize(inisize,ImGuiCond_FirstUseEver)
      end if

    end subroutine init_window

    !> Clamp an initial window size to what fits in the application
    !> window, leaving room for the title bar and the drop shadow. A
    !> window sized in characters overflows the display when the user
    !> raises the font scale, and whatever sits on its right (the side
    !> pane of a file dialog, say) is then cut off.
    subroutine clamp_to_display(siz)
      use gui_main, only: io
      type(ImVec2), intent(inout) :: siz

      real(c_float), parameter :: margin = 2._c_float

      if (io%DisplaySize%x > 0._c_float) &
         siz%x = min(siz%x,io%DisplaySize%x - margin * fontsize%x)
      if (io%DisplaySize%y > 0._c_float) &
         siz%y = min(siz%y,io%DisplaySize%y - margin * fontsize%y)

    end subroutine clamp_to_display

    !> Initial file name and path for a save dialog, taken from the
    !> caller window's current okfile if available; defname in ./
    !> otherwise. Both are returned null-terminated.
    subroutine dialog_initial_file(idparent,defname,file,path)
      use param, only: dirsep
      integer, intent(in) :: idparent
      character(len=*), intent(in) :: defname
      character(kind=c_char,len=:), allocatable, target, intent(out) :: file, path

      integer :: idx
      logical :: ok

      file = defname // c_null_char
      path = "./" // c_null_char
      ok = (idparent >= 1 .and. idparent <= nwin)
      if (ok) ok = allocated(win(idparent)%okfile)
      if (ok) ok = (len_trim(win(idparent)%okfile) > 0)
      if (ok) then
         idx = index(win(idparent)%okfile,dirsep,back=.true.)
         if (idx > 1) path = win(idparent)%okfile(1:idx-1) // c_null_char
         if (idx < len_trim(win(idparent)%okfile)) &
            file = trim(win(idparent)%okfile(idx+1:)) // c_null_char
      end if

    end subroutine dialog_initial_file

  end subroutine window_draw

  !xx! private procedures

  ! the callback for the right-hand-side pane of the dialog
  subroutine dialog_user_callback(vFilter, vUserData, vCantContinue)
    use gui_main, only: g
    use utils, only: igIsItemHovered_delayed, iw_tooltip, iw_text, iw_radiobutton, iw_combo_simple,&
       iw_checkbox, iw_inputfloat
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
       ldum = iw_inputfloat(str,data%rborder)
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

  !xx! pairpick: a staged atom pick for two-atom operations

  !> Stamp the pick time without staging an atom (the pick session
  !> starts now; earlier geometry changes do not invalidate it).
  module subroutine pairpick_arm(pp)
    use interfaces_glfw, only: glfwGetTime
    class(pairpick), intent(inout) :: pp

    pp%idx = 0
    pp%time = glfwGetTime()

  end subroutine pairpick_arm

  !> Stage an atom (cell index + lattice vector) and stamp the time.
  module subroutine pairpick_stage(pp,idx)
    use interfaces_glfw, only: glfwGetTime
    class(pairpick), intent(inout) :: pp
    integer, intent(in) :: idx(4)

    pp%idx = idx
    pp%time = glfwGetTime()

  end subroutine pairpick_stage

  !> Forget the staged atom and the time stamp.
  module subroutine pairpick_clear(pp)
    class(pairpick), intent(inout) :: pp

    pp%idx = 0
    pp%time = 0d0

  end subroutine pairpick_clear

  !> Whether an atom is staged.
  pure module function pairpick_is_staged(pp)
    class(pairpick), intent(in) :: pp
    logical :: pairpick_is_staged

    pairpick_is_staged = (pp%idx(1) > 0)

  end function pairpick_is_staged

  !> Whether the stored pick is stale: the geometry changed after the
  !> stamp, so any stored cell-atom index may point elsewhere.
  pure module function pairpick_is_stale(pp,timelastchange)
    class(pairpick), intent(in) :: pp
    real*8, intent(in) :: timelastchange
    logical :: pairpick_is_stale

    pairpick_is_stale = (timelastchange > pp%time)

  end function pairpick_is_stale

  !> Whether the staged atom (index and lattice vector) equals idx.
  pure module function pairpick_same(pp,idx)
    class(pairpick), intent(in) :: pp
    integer, intent(in) :: idx(4)
    logical :: pairpick_same

    pairpick_same = all(pp%idx == idx)

  end function pairpick_same

  !> Color of the atom with index iat in the atom list of type itype
  !> (system isys) as drawn by the first shown atoms representation in
  !> view window iview. Returns .true. and fills rgb if found. Used to
  !> paint atom buttons in the windows attached to a view.
  module function atom_view_rgb(iview,isys,itype,iat,rgb) result(have)
    use representations, only: reptype_atoms
    use systems, only: sysc
    integer, intent(in) :: iview, isys, itype, iat
    real(c_float), intent(out) :: rgb(3)
    logical :: have

    integer :: jrep, idd

    have = .false.
    rgb = 0._c_float
    if (iview < 1 .or. iview > nwin) return
    if (.not.win(iview)%isinit) return
    ! The style index is computed from isys but indexes the view's scene, so the
    ! two must agree. They can disagree for one frame after a window retargets
    ! its view, and the style arrays are sized per system.
    if (win(iview)%isys /= isys) return
    if (.not.associated(win(iview)%sc)) return
    do jrep = 1, win(iview)%sc%nrep
       if (win(iview)%sc%rep(jrep)%type == reptype_atoms .and. win(iview)%sc%rep(jrep)%isinit .and.&
          win(iview)%sc%rep(jrep)%shown) then
          idd = sysc(isys)%attype_type_id_to_id(itype,iat,win(iview)%sc%rep(jrep)%atoms%style%type)
          if (idd /= 0) then
             have = .true.
             rgb = win(iview)%sc%rep(jrep)%atoms%style%rgb(:,idd)
             return
          end if
       end if
    end do

  end function atom_view_rgb

  !> Short atom-anchor label: atom name + "#" + cell-atom index, or `notset`
  !> when idx does not name a valid cell atom of system isys. With species,
  !> the species name and a space are used instead, which is the shorter form
  !> the atom buttons want. In a crystal, an atom outside the (0,0,0) cell
  !> carries its lattice vector.
  module function anchor_label(isys,idx,notset,species) result(s)
    use systems, only: sys
    use tools_io, only: string
    integer, intent(in) :: isys
    integer(c_int), intent(in) :: idx(4)
    character(len=*), intent(in) :: notset
    logical, intent(in), optional :: species
    character(len=:), allocatable :: s

    logical :: species_

    species_ = .false.
    if (present(species)) species_ = species

    if (idx(1) < 1 .or. idx(1) > sys(isys)%c%ncel) then
       s = notset
    else
       if (species_) then
          s = trim(sys(isys)%c%spc(sys(isys)%c%atcel(idx(1))%is)%name) // " " // string(idx(1))
       else
          s = trim(sys(isys)%c%at(sys(isys)%c%atcel(idx(1))%idx)%name) // "#" // string(idx(1))
       end if
       if (any(idx(2:4) /= 0)) &
          s = s // "+(" // string(idx(2)) // "," // string(idx(3)) // "," // string(idx(4)) // ")"
    end if

  end function anchor_label

end submodule proc
