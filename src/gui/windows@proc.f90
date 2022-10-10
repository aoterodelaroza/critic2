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

  ! Count unique IDs for keeping track of windows and widget
  integer :: idcount = 0

  ! initial side for the view texture
  integer(c_int), parameter :: initial_texture_side = 1024_c_int

  !xx! private procedures
  ! subroutine dialog_user_callback(vFilter, vUserData, vCantContinue)

contains

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
    end if

  end subroutine stack_realloc_maybe

  !> Create a window in the window stack with the given type. Returns
  !> the window ID.
  module function stack_create_window(type,isopen,purpose,isys)
    use tools_io, only: ferror, faterr
    use windows, only: window, nwin, win
    integer, intent(in) :: type
    logical, intent(in) :: isopen
    integer, intent(in), optional :: purpose
    integer, intent(in), optional :: isys

    integer :: stack_create_window

    integer :: i, id
    integer, parameter :: maxwin = 40

    ! find the first unused window or create a new one
    id = 0
    do i = 1, nwin
       if (.not.win(i)%isinit) then
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
    call win(id)%init(type,isopen,purpose,isys)
    stack_create_window = id

  end function stack_create_window

  !> Check whether the window with the given id is still open. If it is
  !> not, or id points to an invalid window, set it to id = 0. If changed
  !> is present, set it to the old id if the id has been changed.
  module subroutine update_window_id(id,changed)
    integer, intent(inout) :: id
    integer, intent(out), optional :: changed

    integer :: oid

    oid = id
    id = 0
    if (present(changed)) changed = oid
    if (oid < 1 .or. oid > nwin) return
    if (.not.win(oid)%isinit .or. .not.win(oid)%isopen) return
    id = oid
    if (present(changed)) changed = 0

  end subroutine update_window_id

  !> Initialize a window of the given type. If isiopen, initialize it
  !> as open.
  module subroutine window_init(w,type,isopen,purpose,isys)
    use interfaces_opengl3
    use gui_main, only: ColorDialogDir, ColorDialogFile
    use tools_io, only: ferror, faterr
    use param, only: bohrtoa
    class(window), intent(inout), target :: w
    integer, intent(in) :: type
    logical, intent(in) :: isopen
    integer, intent(in), optional :: purpose
    integer, intent(in), optional :: isys

    character(kind=c_char,len=:), allocatable, target :: str1

    ! initialization of the state
    w%isinit = .true.
    w%firstpass = .true.
    w%isopen = isopen
    w%type = type
    w%id = -1
    w%name = "" // c_null_char
    w%table_selected = 1
    w%table_sortcid = 0
    w%table_sortdir = 1
    w%forceresize = .false.
    w%forcesort = .false.
    w%forceupdate = .false.
    w%forceinit = .false.
    w%inpcon_selected = 1
    w%okfile_set = .false. ! whether the library file has been set by the user
    w%okfile_read = .false. ! whether the structure list should be re-read from the lib
    w%view_selected = 1
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
    w%forcequitdialog = .false.

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
       w%dialog_purpose = purpose
    elseif (type == wintype_load_field) then
       ! dialog: load field window
       if (.not.present(isys)) &
          call ferror('window_init','load_field requires isys',faterr)
       w%loadfield_isys = isys
    elseif (type == wintype_scfplot) then
       ! SCF plot window
       if (.not.present(isys)) &
          call ferror('window_init','scfplot requires isys',faterr)
       w%scfplot_isys = isys
    elseif (type == wintype_view) then
       ! view window
       call w%create_texture_view(initial_texture_side)
    end if

  end subroutine window_init

  !> End a window and deallocate the data.
  module subroutine window_end(w)
    use interfaces_opengl3
    class(window), intent(inout), target :: w

    ! window-specific destruction
    if (w%isinit) then
       if (w%type == wintype_dialog .and. c_associated(w%dptr)) then
          call IGFD_Destroy(w%dptr)
       elseif (w%type == wintype_view) then
          call w%delete_texture_view()
       end if
    end if

    ! deallocate the rest of the data
    w%isinit = .false.
    w%isopen = .false.
    w%id = -1
    w%name = "" // c_null_char
    if (allocated(w%iord)) deallocate(w%iord)
    if (allocated(w%forceremove)) deallocate(w%forceremove)

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
    use utils, only: iw_text
    use tools_io, only: string, ferror, faterr
    class(window), intent(inout), target :: w

    character(kind=c_char,len=:), allocatable, target :: str1, str2
    type(ImVec2) :: inisize

    if (.not.w%isinit) return
    if (.not.w%isopen) return

    ! First pass on creation: assign ID, name, and flags
    if (w%id < 0) then
       w%firstpass = .true.
       idcount = idcount + 1
       w%id = idcount
       if (w%type == wintype_tree) then
          w%name = "Tree" // c_null_char
          w%flags = ImGuiWindowFlags_None
       elseif (w%type == wintype_view) then
          w%name = "View" // c_null_char
          w%flags = ImGuiWindowFlags_None
       elseif (w%type == wintype_console_input) then
          w%name = "Input Console" // c_null_char
          w%flags = ImGuiWindowFlags_None
       elseif (w%type == wintype_console_output) then
          w%name = "Output Console" // c_null_char
          w%flags = ImGuiWindowFlags_None
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
             w%name = "Open File(s)..." // c_null_char
             str1 = &
                "&
                &All files (*.*){*.*},&
                &ABINIT (DEN...){(DEN|ELF|POT|VHA\VHXC|VXC|GDEN1|GDEN2|GDEN3|LDEN|KDEN|PAWDEN|VCLMB|VPSP)},&
                &ADF (molden){.molden},&
                &cif (cif){.cif},&
                &CRYSTAL (out){.out},&
                &cube (cube|bincube){.cube,.bincube},&
                &DFTB+ (gen){.gen},&
                &DMACRYS (dmain|16|21){.dmain,.16,.21},&
                &elk (OUT){.OUT},&
                &FHIaims (in|in.next_step|out|own){.in,.next_step,.out,.own},&
                &Gaussian (log|wfn|wfx|fchk|cube){.log,.wfn,.wfx,.fchk,.cube},&
                &ORCA (molden|molden.input){.molden,.input},&
                &postg (pgout){.pgout},&
                &psi4 (molden|dat){.molden,.dat},&
                &Quantum ESPRESSO (out|in|cube|pwc) {.out,.in,.cube,.pwc},&
                &SHELX (res|ins){.res,.ins.16},&
                &SIESTA (STRUCT_IN|STRUCT_OUT) {.STRUCT_IN,.STRUCT_OUT},&
                &TINKER (frac) {.frac},&
                &VASP (POSCAR|CONTCAR|...){(CONTCAR|CHGCAR|ELFCAR|CHG|AECCAR0|AECCAR1|AECCAR2|POSCAR)},&
                &WIEN2k (struct){.struct},&
                &Xcrysden (xsf|axsf) {.xsf,.axsf},&
                &xyz (xyz){.xyz},&
                &"// c_null_char
             call IGFD_OpenPaneDialog2(w%dptr,c_loc(w%name),c_loc(w%name),c_loc(str1),c_loc(str2),&
                c_funloc(dialog_user_callback),280._c_float,0_c_int,c_loc(w%dialog_data),w%flags)
          elseif (w%dialog_purpose == wpurp_dialog_savelogfile) then
             w%name = "Save Log File..." // c_null_char
             call IGFD_OpenPaneDialog2(w%dptr,c_loc(w%name),c_loc(w%name),c_loc(str1),c_loc(str2),&
                c_funloc(dialog_user_callback),280._c_float,1_c_int,c_loc(w%dialog_data),w%flags)
          elseif (w%dialog_purpose == wpurp_dialog_openlibraryfile) then
             w%name = "Open Library File..." // c_null_char
             call IGFD_OpenPaneDialog2(w%dptr,c_loc(w%name),c_loc(w%name),c_loc(str1),c_loc(str2),&
                c_funloc(dialog_user_callback),280._c_float,1_c_int,c_loc(w%dialog_data),w%flags)
          elseif (w%dialog_purpose == wpurp_dialog_openfieldfile) then
             w%name = "Open Field File(s)..." // c_null_char
             str1 = &
                "&
                &All files (*.*){*.*},&
                &ABINIT (DEN...){(DEN|ELF|POT|VHA\VHXC|VXC|GDEN1|GDEN2|GDEN3|LDEN|KDEN|PAWDEN|VCLMB|VPSP)},&
                &Aimpac (qub){.qub},&
                &Binary cube file (bincube){.bincube},&
                &Cube file (cube){.cube},&
                &DFTB+ (detailed.xml){.xml},&
                &elk (grid){.grid},&
                &elk STATE.OUT (OUT){.OUT},&
                &Gaussian wavefunction (wfn|wfx|fchk){.wfn,.wfx,.fchk},&
                &Molden-style file (molden){.molden},&
                &Quantum ESPRESSO pwc file (pwc){.pwc},&
                &SIESTA (RHO...) {(RHO|BADER|DRHO|LDOS|VT|VH)},&
                &VASP ELFCAR{(ELFCAR)},&
                &VASP (CHGCAR...){(CONTCAR|CHGCAR|ELFCAR|CHG|AECCAR0|AECCAR1|AECCAR2|POSCAR)},&
                &WIEN2k (clmsum...){.clmsum,.clmup,.clmdn},&
                &Xcrysden (xsf|axsf) {.xsf,.axsf},&
                &"// c_null_char
             call IGFD_OpenPaneDialog2(w%dptr,c_loc(w%name),c_loc(w%name),c_loc(str1),c_loc(str2),&
                c_funloc(dialog_user_callback),280._c_float,1_c_int,c_loc(w%dialog_data),w%flags)
          elseif (w%dialog_purpose == wpurp_dialog_openonefilemodal) then
             w%name = "Open File..." // c_null_char
             call IGFD_OpenPaneDialog2(w%dptr,c_loc(w%name),c_loc(w%name),c_loc(str1),c_loc(str2),&
                c_funloc(dialog_user_callback),280._c_float,1_c_int,c_loc(w%dialog_data),w%flags)
          else
             call ferror('window_draw','unknown dialog purpose',faterr)
          end if
       elseif (w%type == wintype_new_struct) then
          w%name = "New Structure..." // c_null_char
          w%flags = ImGuiWindowFlags_None
          inisize%x = 90 * fontsize%x
          inisize%y = 40 * fontsize%y
          call igSetNextWindowSize(inisize,ImGuiCond_FirstUseEver)
       elseif (w%type == wintype_new_struct_library) then
          w%name = "New Structure from Library..." // c_null_char
          w%flags = ImGuiWindowFlags_None
          inisize%x = 90 * fontsize%x
          inisize%y = 30 * fontsize%y
          call igSetNextWindowSize(inisize,ImGuiCond_FirstUseEver)
       elseif (w%type == wintype_load_field) then
          w%name = "Load Field..." // c_null_char
          w%flags = ImGuiWindowFlags_None
          inisize%x = 60 * fontsize%x
          inisize%y = 15 * fontsize%y
          call igSetNextWindowSize(inisize,ImGuiCond_FirstUseEver)
       elseif (w%type == wintype_scfplot) then
          w%name = "SCF Iterations (" // string(w%scfplot_isys) // ")" // c_null_char
          w%flags = ImGuiWindowFlags_None
          inisize%x = 45 * fontsize%x
          inisize%y = inisize%x
          call igSetNextWindowSize(inisize,ImGuiCond_FirstUseEver)
       end if
    end if

    if (w%isopen) then
       ! assign the pointer ID for the window, if not a dialog
       if (w%type == wintype_dialog) then
          w%ptr = IGFD_GetCurrentWindow(w%dptr)
          call w%draw_dialog()
       else
          if (igBegin(c_loc(w%name),w%isopen,w%flags)) then
             ! draw the window contents, depending on type
             w%ptr = igGetCurrentWindow()
             if (w%type == wintype_tree) then
                call w%draw_tree()
             elseif (w%type == wintype_view) then
                call w%draw_view()
             elseif (w%type == wintype_console_input) then
                call w%draw_ci()
             elseif (w%type == wintype_console_output) then
                call w%draw_co()
             elseif (w%type == wintype_new_struct) then
                call w%draw_new_struct()
             elseif (w%type == wintype_new_struct_library) then
                call w%draw_new_struct_from_library()
             elseif (w%type == wintype_load_field) then
                call w%draw_load_field()
             elseif (w%type == wintype_scfplot) then
                call w%draw_scfplot()
             end if
          end if
          call igEnd()
       end if
       w%firstpass = .false.
    else
       w%firstpass = .true.
    end if

  end subroutine window_draw

  !xx! private procedures

  ! the callback for the right-hand-side pane of the dialog
  subroutine dialog_user_callback(vFilter, vUserData, vCantContinue) bind(c)
    use gui_main, only: g
    use utils, only: igIsItemHovered_delayed, iw_tooltip, iw_text, iw_radiobutton,&
       iw_combo_simple
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
    str = "Show hidden files" // c_null_char
    if (igCheckbox(c_loc(str),data%showhidden)) then
       flags = IGFD_GetFlags(data%dptr)
       if (data%showhidden) then
          flags = iand(flags,not(ImGuiFileDialogFlags_DontShowHiddenFiles))
       else
          flags = ior(flags,ImGuiFileDialogFlags_DontShowHiddenFiles)
       end if
       call IGFD_SetFlags(data%dptr,flags)
    end if
    call iw_tooltip("Show the OS hidden files and directories in this dialog",ttshown)

    !! options specific to the open files dialog !!
    if (data%purpose == wpurp_dialog_openfiles) then
       str = "Read last structure only" // c_null_char
       ldum = igCheckbox(c_loc(str),data%readlastonly)
       call iw_tooltip("Read only the last structure in the file",ttshown)
       call igNewLine()

       ! radio buttons for auto/crystal/molecule
       call iw_text("Read structures as...",highlight=.true.)
       ldum = iw_radiobutton("Auto-detect",int=data%mol,intval=-1_c_int)
       call iw_tooltip("Auto-detect whether new structures are read as crystals or molecules",ttshown)
       ldum = iw_radiobutton("Crystal",int=data%mol,intval=0_c_int)
       call iw_tooltip("Force new structures to be read as crystals",ttshown)
       ldum = iw_radiobutton("Molecule",int=data%mol,intval=2_c_int)
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
       call iw_tooltip("Periodic cell border around new molecules",ttshown)
       call igPopItemWidth()

       str = "Cubic cell" // c_null_char
       ldum = igCheckbox(c_loc(str),data%molcubic)
       call iw_tooltip("Read new molecules inside cubic periodic cell",ttshown)
       call igUnindent(0._c_float)
       call igNewLine()

       ! Input structure format (isformat)
       call iw_text("Read file format",highlight=.true.)
       stropt = "" &
          // "Auto-detect" // c_null_char &             ! isformat_unknown = 0
          // "Abinit DEN-style file" // c_null_char &   ! isformat_abinit = 7
          // "Binary cube" // c_null_char &             ! isformat_bincube = 5
          // "CIF file" // c_null_char &                ! isformat_cif = 1
          // "CRYSTAL output" // c_null_char &          ! isformat_crystal = 11
          // "Cube" // c_null_char &                    ! isformat_cube = 4
          // "DFTB+ gen file" // c_null_char &          ! isformat_gen = 20
          // "DMACRYS .21 file" // c_null_char &        ! isformat_f21 = 3
          // "DMACRYS dmain file" // c_null_char &      ! isformat_dmain = 27
          // "elk GEOMETRY.OUT" // c_null_char &        ! isformat_elk = 8
          // "FHIaims input" // c_null_char &           ! isformat_aimsin = 28
          // "FHIaims output" // c_null_char &          ! isformat_aimsout = 29
          // "Gaussian fchk" // c_null_char &           ! isformat_fchk = 15
          // "Gaussian output" // c_null_char &         ! isformat_gaussian = 17
          // "Gaussian wfn" // c_null_char &            ! isformat_wfn = 13
          // "Gaussian wfx" // c_null_char &            ! isformat_wfx = 14
          // "ORCA molden file" // c_null_char &        ! isformat_orca = 26
          // "postg output" // c_null_char &            ! isformat_pgout = 25
          // "psi4 molden file" // c_null_char &        ! isformat_molden = 16
          // "psi4 output" // c_null_char &             ! isformat_dat = 24
          // "Quantum ESPRESSO input" // c_null_char &  ! isformat_qein = 9
          // "Quantum ESPRESSO output" // c_null_char & ! isformat_qeout = 10
          // "Quantum ESPRESSO pwc" // c_null_char &    ! isformat_pwc = 22
          // "SHELX" // c_null_char &                   ! isformat_shelx = 2
          // "SIESTA IN/OUT file" // c_null_char &      ! isformat_siesta = 18
          // "TINKER frac file" // c_null_char &        ! isformat_tinkerfrac = 30
          // "VASP" // c_null_char &                    ! isformat_vasp = 21
          // "WIEN2k struct file" // c_null_char &      ! isformat_struct = 6
          // "Xcrysden axsf file" // c_null_char &      ! isformat_axsf = 23
          // "Xcrysden xsf file" // c_null_char &       ! isformat_xsf = 19
          // "xyz file" // c_null_char                  ! isformat_xyz = 12
       call iw_combo_simple("##formatcombo",stropt,data%isformat)
       call iw_tooltip("Force the new structure to be read with this file format, or auto-detect&
          &from the extension",ttshown)
       call igNewLine()
    elseif (data%purpose == wpurp_dialog_openfieldfile) then
       ! Input field format (ifformat)
       call iw_text("Read file format",highlight=.true.)
       stropt = "" &
          // "Auto-detect" // c_null_char &              ! ifformat_unknown = 0
          // "ABINIT DEN-style file" // c_null_char &    ! ifformat_abinit = 6
          // "Aimpac qub" // c_null_char &               ! ifformat_qub = 9
          // "Binary cube" // c_null_char &              ! ifformat_bincube = 5
          // "Cube file" // c_null_char &                ! ifformat_cube = 4
          // "DFTB+ detailed.xml" // c_null_char &       ! ifformat_dftb = 13
          // "elk grid" // c_null_char &                 ! ifformat_elkgrid = 11
          // "elk STATE.OUT" // c_null_char &            ! ifformat_elk = 2
          // "Gaussian wfn" // c_null_char &             ! ifformat_wfn = 15
          // "Gaussian wfx" // c_null_char &             ! ifformat_wfx = 16
          // "Gaussian fchk" // c_null_char &            ! ifformat_fchk = 17
          // "Molden file" // c_null_char &              ! ifformat_molden = 18
          // "Quantum ESPRESSO pwc" // c_null_char &     ! ifformat_pwc = 14
          // "SIESTA RHO-style file" // c_null_char &    ! ifformat_siestagrid = 12
          // "VASP ELFCAR-style file" // c_null_char &   ! ifformat_vaspnov = 8
          // "VASP CHGCAR-style file" // c_null_char &   ! ifformat_vasp = 7
          // "WIEN2k clmsum-style file" // c_null_char & ! ifformat_wien = 1
          // "Xcrysden xsf" // c_null_char               ! ifformat_xsf = 10
       call iw_combo_simple("##formatcombo",stropt,data%isformat)
       call iw_tooltip("Force the new field to be read with the given file format, or auto-detect&
          &from the extension",ttshown)
       call igNewLine()
    end if

  end subroutine dialog_user_callback

end submodule proc
