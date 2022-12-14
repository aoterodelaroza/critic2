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
  module function stack_create_window(type,isopen,purpose,isys,irep,idcaller)
    use tools_io, only: ferror, faterr
    use windows, only: window, nwin, win
    integer, intent(in) :: type
    logical, intent(in) :: isopen
    integer, intent(in), optional :: purpose
    integer, intent(in), optional :: isys
    integer, intent(in), optional :: irep
    integer, intent(in), optional :: idcaller

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
    call win(id)%init(type,isopen,id,purpose,isys,irep,idcaller)
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
  module subroutine window_init(w,type,isopen,id,purpose,isys,irep,idcaller)
    use interfaces_opengl3
    use gui_main, only: ColorDialogDir, ColorDialogFile, sysc
    use tools_io, only: ferror, faterr
    use param, only: bohrtoa
    class(window), intent(inout), target :: w
    integer, intent(in) :: type
    logical, intent(in) :: isopen
    integer, intent(in) :: id
    integer, intent(in), optional :: purpose
    integer, intent(in), optional :: isys
    integer, intent(in), optional :: irep
    integer, intent(in), optional :: idcaller

    character(kind=c_char,len=:), allocatable, target :: str1

    ! initialization of the state
    w%isinit = .true.
    w%firstpass = .true.
    w%isopen = isopen
    w%type = type
    w%id = -id
    w%name = "" // c_null_char
    w%errmsg = ""
    w%table_selected = 1
    w%table_sortcid = 0
    w%table_sortdir = 1
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
    w%view_mousebehavior = MB_navigation
    w%idexportwin = 0
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
    w%idsave = 0

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
       w%isys = isys
    elseif (type == wintype_scfplot) then
       ! SCF plot window
       if (.not.present(isys)) &
          call ferror('window_init','scfplot requires isys',faterr)
       w%isys = isys
    elseif (type == wintype_editrep) then
       ! edit representation window
       if (.not.present(isys)) &
          call ferror('window_init','editrep requires isys',faterr)
       if (.not.present(irep)) &
          call ferror('window_init','editrep requires irep',faterr)
       if (.not.present(idcaller)) &
          call ferror('window_init','editrep requires idcaller',faterr)
       w%isys = isys
       w%rep => sysc(isys)%sc%rep(irep)
       w%idparent = idcaller
    elseif (type == wintype_exportimage) then
       ! export image window
       if (.not.present(idcaller)) &
          call ferror('window_init','exportimage requires idcaller',faterr)
       w%idparent = idcaller
    elseif (type == wintype_view) then
       ! view window
       call w%create_texture_view(initial_texture_side)
    elseif (type == wintype_rebond) then
       ! recalculate bonds window
       if (.not.present(isys)) &
          call ferror('window_init','rebond requires isys',faterr)
       w%isys = isys
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
    if (allocated(w%plotx)) deallocate(w%plotx)
    if (allocated(w%ploty)) deallocate(w%ploty)
    nullify(w%rep)

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

    character(kind=c_char,len=:), allocatable, target :: str1, str2, str3
    type(ImVec2) :: inisize

    if (.not.w%isinit) return
    if (.not.w%isopen) return

    ! First pass on creation: assign ID, name, and flags
    if (w%id < 0) then
       w%firstpass = .true.
       w%id = abs(w%id)
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
       elseif (w%type == wintype_about) then
          w%name = "About" // c_null_char
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
             w%name = "Open File(s)" // c_null_char
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
             w%name = "Save Log File" // c_null_char
             str2 = "file.log" // c_null_char
             str3 = "./" // c_null_char
             call IGFD_OpenPaneDialog(w%dptr,c_loc(w%name),c_loc(w%name),c_loc(str1),c_loc(str3),c_loc(str2),&
                c_funloc(dialog_user_callback),280._c_float,1_c_int,c_loc(w%dialog_data),w%flags)
          elseif (w%dialog_purpose == wpurp_dialog_openlibraryfile) then
             w%name = "Open Library File" // c_null_char
             call IGFD_OpenPaneDialog2(w%dptr,c_loc(w%name),c_loc(w%name),c_loc(str1),c_loc(str2),&
                c_funloc(dialog_user_callback),280._c_float,1_c_int,c_loc(w%dialog_data),w%flags)
          elseif (w%dialog_purpose == wpurp_dialog_openfieldfile) then
             w%name = "Open Field File(s)" // c_null_char
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
             w%name = "Open File" // c_null_char
             call IGFD_OpenPaneDialog2(w%dptr,c_loc(w%name),c_loc(w%name),c_loc(str1),c_loc(str2),&
                c_funloc(dialog_user_callback),280._c_float,1_c_int,c_loc(w%dialog_data),w%flags)
          elseif (w%dialog_purpose == wpurp_dialog_saveimagefile) then
             w%name = "Save Image File" // c_null_char
             str1 = "PNG (*.png) {.png},BMP (*.bmp) {.bmp},TGA (*.tga) {.tga},JPEG (*.jpg) {.jpg}"// c_null_char
             str2 = "file.png" // c_null_char
             str3 = "./" // c_null_char
             call IGFD_OpenPaneDialog(w%dptr,c_loc(w%name),c_loc(w%name),c_loc(str1),c_loc(str3),c_loc(str2),&
                c_funloc(dialog_user_callback),280._c_float,1_c_int,c_loc(w%dialog_data),w%flags)
          else
             call ferror('window_draw','unknown dialog purpose',faterr)
          end if
       elseif (w%type == wintype_new_struct) then
          w%name = "New Structure" // c_null_char
          w%flags = ImGuiWindowFlags_None
          inisize%x = 90 * fontsize%x
          inisize%y = 40 * fontsize%y
          call igSetNextWindowSize(inisize,ImGuiCond_FirstUseEver)
       elseif (w%type == wintype_new_struct_library) then
          w%name = "New Structure from Library" // c_null_char
          w%flags = ImGuiWindowFlags_None
          inisize%x = 90 * fontsize%x
          inisize%y = 30 * fontsize%y
          call igSetNextWindowSize(inisize,ImGuiCond_FirstUseEver)
       elseif (w%type == wintype_load_field) then
          w%name = "Load Field" // c_null_char
          w%flags = ImGuiWindowFlags_None
          inisize%x = 60 * fontsize%x
          inisize%y = 15 * fontsize%y
          call igSetNextWindowSize(inisize,ImGuiCond_FirstUseEver)
       elseif (w%type == wintype_scfplot) then
          w%name = "SCF Iterations (" // string(w%isys) // ")" // c_null_char
          w%flags = ImGuiWindowFlags_None
          inisize%x = 45 * fontsize%x
          inisize%y = inisize%x
          call igSetNextWindowSize(inisize,ImGuiCond_FirstUseEver)
       elseif (w%type == wintype_editrep) then
          w%name = "Representation [" // string(w%rep%name) // ", " // &
             string(w%isys) // "]##" // string(w%rep%idrep) // c_null_char
          w%flags = ImGuiWindowFlags_None
          inisize%x = 60 * fontsize%x
          inisize%y = 44 * fontsize%y
          call igSetNextWindowSize(inisize,ImGuiCond_FirstUseEver)
       elseif (w%type == wintype_exportimage) then
          w%name = "Export Image" // "##" // string(w%idparent) // c_null_char
          w%flags = ImGuiWindowFlags_None
          inisize%x = 50 * fontsize%x
          inisize%y = 17 * fontsize%y
          call igSetNextWindowSize(inisize,ImGuiCond_FirstUseEver)
       elseif (w%type == wintype_rebond) then
          w%name = "Recalculate Bonds" // c_null_char
          w%flags = ImGuiWindowFlags_None
          inisize%x = 55 * fontsize%x
          inisize%y = 19 * fontsize%y
          call igSetNextWindowSize(inisize,ImGuiCond_FirstUseEver)
       elseif (w%type == wintype_preferences) then
          w%name = "Preferences..." // c_null_char
          w%flags = ImGuiWindowFlags_None
          inisize%x = 55 * fontsize%x
          inisize%y = 19 * fontsize%y
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
             elseif (w%type == wintype_rebond) then
                call w%draw_rebond()
             elseif (w%type == wintype_preferences) then
                call w%draw_preferences()
             end if
          end if
          call igEnd()
       end if
       w%firstpass = .false.
    else
       w%firstpass = .true.
    end if

  end subroutine window_draw

  !> Draw the about window
  module subroutine draw_about(w)
    use gui_main, only: g
    use utils, only: iw_text, iw_button, iw_calcwidth
    use keybindings, only: is_bind_event, BIND_CLOSE_FOCUSED_DIALOG, BIND_CLOSE_ALL_DIALOGS,&
       BIND_OK_FOCUSED_DIALOG
    class(window), intent(inout), target :: w

    real(c_float) :: wwidth, twidth
    type(ImVec2) :: szavail

    call iw_text("  ---- critic2 GUI ----",danger=.true.,centered=.true.)
    call iw_text("[development version]",centered=.true.)
    call igNewLine()
    call iw_text("Critic2 is a program for visualization and analysis",centered=.true.)
    call iw_text("of structures and data in computational chemistry",centered=.true.)
    call igNewLine()
    call iw_text("Copyright (c) 2022- Alberto Otero de la Roza",centered=.true.)
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
    use gui_main, only: g, tooltip_enabled, tooltip_delay, tooltip_wrap_factor
    use interfaces_cimgui
    use keybindings, only: is_bind_event, BIND_CLOSE_FOCUSED_DIALOG, BIND_CLOSE_ALL_DIALOGS,&
       BIND_OK_FOCUSED_DIALOG
    use utils, only: iw_tooltip, iw_button, iw_text, iw_calcwidth
    class(window), intent(inout), target :: w

    character(kind=c_char,len=:), allocatable, target :: str, str2, zeroc
    logical(c_bool) :: ldum
    logical :: doquit
    type(ImVec2) :: sz, szero

    logical, save :: ttshown = .false. ! tooltip flag
    type(c_ptr), save :: cfilter = c_null_ptr ! filter object (allocated first pass, never destroyed)
    integer(c_int), save :: catid = 0 ! category ID (from left panel)

    real(c_float), parameter :: wleft = 150._c_float
    real(c_float), parameter :: wright = 400._c_float

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
       call iw_tooltip("Settings for the graphical user interface",ttshown)
       str = "Key bindings" // c_null_char
       if (igSelectable_Bool(c_loc(str),logical(catid == 1,c_bool),ImGuiSelectableFlags_None,szero)) catid = 1
       call iw_tooltip("User interface key bindings",ttshown)
       call igEndChild()
    end if
    call igSameLine(0._c_float,-1._c_float)

    ! right panel
    call igBeginGroup()
    str = "rightpanel" // c_null_char
    sz%x = 0._c_float
    sz%y = -igGetFrameHeightWithSpacing() - g%Style%ItemSpacing%y
    if (igBeginChild_Str(c_loc(str),sz,.true._c_bool,ImGuiWindowFlags_None)) then

       call igPushItemWidth(4._c_float * g%FontSize)
       if (catid == 0) then
          !! Interface
          str = "Enable tooltips" // c_null_char
          if (ImGuiTextFilter_PassFilter(cfilter,c_loc(str),c_null_ptr)) then
             ldum = igCheckbox(c_loc(str), tooltip_enabled)
             call iw_tooltip("Show/hide the tooltips in the user interface elements",ttshown)
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

       elseif (catid == 1) then
          !! key bindings

  !       static int getbind = -1;
  !       TextDisabled("(Right-click on the button to toggle double-click)");

  !       // Key bindings
  !       BeginChild("keybindselector",ImVec2(wright,0.f),false);
  !       Columns(2,"keybindingcolumns",false);
  !       for (int i = 0; i < BIND_MAX; i++){
  !         if (!filter.PassFilter(BindNames[i]))
  !           continue;
  !         string result = BindKeyName(i);
  !         if (result.length() == 0)
  !           result = "<not bound>";
  !         Text(BindNames[i]);
  !         NextColumn();

  !         PushID(i);
  !         if (Button(result.c_str())){
  !           getbind = i;
  !         }
  !         PopID();

  !         if (IsMouseClicked(1) && IsItemHovered()){
  !           int newkey;
  !           switch (keybind[i]){
  !           case GLFW_MOUSE_LEFT:
  !             newkey = GLFW_MOUSE_LEFT_DOUBLE; break;
  !           case GLFW_MOUSE_LEFT_DOUBLE:
  !             newkey = GLFW_MOUSE_LEFT; break;
  !           case GLFW_MOUSE_RIGHT:
  !             newkey = GLFW_MOUSE_RIGHT_DOUBLE; break;
  !           case GLFW_MOUSE_RIGHT_DOUBLE:
  !             newkey = GLFW_MOUSE_RIGHT; break;
  !           case GLFW_MOUSE_MIDDLE:
  !             newkey = GLFW_MOUSE_MIDDLE_DOUBLE; break;
  !           case GLFW_MOUSE_MIDDLE_DOUBLE:
  !             newkey = GLFW_MOUSE_MIDDLE; break;
  !           case GLFW_MOUSE_BUTTON3:
  !             newkey = GLFW_MOUSE_BUTTON3_DOUBLE; break;
  !           case GLFW_MOUSE_BUTTON3_DOUBLE:
  !             newkey = GLFW_MOUSE_BUTTON3; break;
  !           case GLFW_MOUSE_BUTTON4:
  !             newkey = GLFW_MOUSE_BUTTON4_DOUBLE; break;
  !           case GLFW_MOUSE_BUTTON4_DOUBLE:
  !             newkey = GLFW_MOUSE_BUTTON4; break;
  !           default:
  !             newkey = NOKEY;
  !           }
  !           if (newkey)
  !             SetBind(i,newkey,modbind[i]);
  !         }
  !         NextColumn();
  !       }
  !       Columns(1);
  !       EndChild();

  !       // popup to read a key
  !       if (getbind != -1){
  !         SetBindEventLevel(1);
  !         OpenPopup("Choosekey");
  !         if (BeginPopupModal("Choosekey", NULL, ImGuiWindowFlags_AlwaysAutoResize | ImGuiWindowFlags_NoTitleBar |
  !                           ImGuiWindowFlags_NoMove)){
  !           Text("Please press a key or mouse button.");
  !           if (SetBindFromUserInput(getbind)){
  !             getbind = -1;
  !             SetBindEventLevel();
  !             CloseCurrentPopup();
  !           }
  !           EndPopup();
  !         }
  !       }
  !     }
       end if
       call igPopItemWidth()

  !     if (updateviews)
  !       ForceUpdateAllViews();
  !   } // BeginDock()

       !! copy tooltips from the previous UI
       !! copy UI options from gui_main and other modules
       !! segfault if the group is too small??
       !! -- implement reset --

       call igEndChild()
    end if

    ! child for final buttons
    str = "finalbuttons" // c_null_char
    sz%x = 0._c_float
    sz%y = 0._c_float
    if (igBeginChild_Str(c_loc(str),sz,.false._c_bool,ImGuiWindowFlags_None)) then
       call igSetCursorPosX(iw_calcwidth(10,2,from_end=.true.) - g%Style%ScrollbarSize)
       if (iw_button("Reset",danger=.true.)) then
          write (*,*) "Bleh!"
       end if
       call iw_tooltip("Reset to the default UI preferences",ttshown)
       if (iw_button("Close",sameline=.true.)) doquit = .true.
       call iw_tooltip("Close this window",ttshown)
       call igEndChild()
    end if
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
       catid = 0
    end subroutine init_state

    ! terminate the state for this window
    subroutine end_state()
       if (c_associated(cfilter)) call ImGuiTextFilter_Clear(cfilter)
    end subroutine end_state

  end subroutine draw_preferences

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
