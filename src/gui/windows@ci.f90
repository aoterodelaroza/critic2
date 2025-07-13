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

! Windows, console input.
submodule (windows) ci
  use interfaces_cimgui
  implicit none

contains

  !> Draw the contents of the input console
  module subroutine draw_ci(w)
    use keybindings, only: BIND_INPCON_RUN, get_bind_keyname, is_bind_event
    use systems, only: sys, sysc, nsys, sys_init
    use gui_main, only: g, force_run_commands, tree_select_updates_inpcon
    use templates, only: draw_keyword_context_menu
    use utils, only: igIsItemHovered_delayed, iw_tooltip, iw_button, iw_text, iw_menuitem
    use systemmod, only: sy
    use tools_io, only: string
    use param, only: newline
    class(window), intent(inout), target :: w

    character(kind=c_char,len=:), allocatable, target :: str1, str2
    type(ImVec2) :: sz, szero, szavail
    logical(c_bool) :: ldum, is_selected
    real(c_float) :: combowidth
    logical :: ok, ok2
    integer :: i, idx

    logical, save :: ttshown = .false. ! tooltip flag

    ! initialize
    szero%x = 0
    szero%y = 0

    ! allocate the input buffer if not already done
    if (.not.allocated(inputb)) then
       allocate(character(len=maxlib+1) :: inputb)
       inputb(1:1) = c_null_char
    end if

    ! update based on time signals between dependent windows
    if (w%timelast_inpcon_assign < win(iwin_tree)%timelast_tree_assign.and.&
       tree_select_updates_inpcon) then
       call assign_system(win(iwin_tree)%tree_selected)
    end if

    ! first line: text
    call iw_text("Input",highlight=.true.)

    ! first line: clear button
    if (iw_button("Clear",sameline=.true.)) &
       inputb(1:1) = c_null_char
    call iw_tooltip("Clear the text from the input console",ttshown)

    ! first line: template button
    ldum = iw_button("Template",sameline=.true.,popupcontext=ok,popupflags=ImGuiPopupFlags_MouseButtonLeft)
    if (ok) &
       call draw_keyword_context_menu(.true.)
    call iw_tooltip("Insert a template for a critic2 command",ttshown)

    ! first line: help button
    ldum = iw_button("Help",sameline=.true.,popupcontext=ok,popupflags=ImGuiPopupFlags_MouseButtonLeft)
    if (ok) &
       call draw_keyword_context_menu(.false.)
    call iw_tooltip("Bring up the critic2 command reference",ttshown)

    ! second line: calculate size of the RUN button
    sz%x = 2 * (igGetTextLineHeight() + 2 * g%Style%FramePadding%y) + g%Style%ItemSpacing%y
    sz%y = sz%x

    ! second line: system selector
    call igBeginGroup()
    call iw_text("System")
    call igSameLine(0._c_float,-1._c_float)

    call igGetContentRegionAvail(szavail)
    combowidth = max(szavail%x - sz%x - g%Style%ItemSpacing%x,0._c_float)

    !! set the system pointer and determine the preview string
    str2 = "" // c_null_char
    sy => null()
    if (w%inpcon_selected >= 1 .and. w%inpcon_selected <= nsys) then
       if (sysc(w%inpcon_selected)%status == sys_init) then
          str2 = string(w%inpcon_selected) // ": " // trim(sysc(w%inpcon_selected)%seed%name) // c_null_char
          sy => sys(w%inpcon_selected)
       end if
    end if
    str1 = "##systemcombo" // c_null_char
    call igSetNextItemWidth(combowidth)
    if (igBeginCombo(c_loc(str1),c_loc(str2),ImGuiComboFlags_None)) then
       do i = 1, nsys
          if (sysc(i)%status == sys_init) then
             is_selected = (w%inpcon_selected == i)
             str2 = string(i) // ": " // trim(sysc(i)%seed%name) // c_null_char
             if (igSelectable_Bool(c_loc(str2),is_selected,ImGuiSelectableFlags_None,szero)) &
                call assign_system(i)
             if (is_selected) &
                call igSetItemDefaultFocus()
          end if
       end do
       call igEndCombo()
    end if
    call iw_tooltip("Set the current system (input commands are applied to it)",ttshown)

    ! third line: field selector
    call iw_text("Field ")
    call igSameLine(0._c_float,-1._c_float)

    str1 = "##fieldcombo" // c_null_char
    str2 = "" // c_null_char
    if (associated(sy)) &
       str2 = string(sy%iref) // ": " // trim(sy%f(sy%iref)%name) // c_null_char
    call igSetNextItemWidth(combowidth)
    if (igBeginCombo(c_loc(str1),c_loc(str2),ImGuiComboFlags_None)) then
       if (associated(sy)) then
          do i = 0, sy%nf
             if (.not.sy%f(i)%isinit) cycle
             is_selected = (sy%iref == i)
             str2 = string(i) // ": " // trim(sy%f(i)%name) // c_null_char
             if (igSelectable_Bool(c_loc(str2),is_selected,ImGuiSelectableFlags_None,szero)) &
                call sy%set_reference(i,.false.)
             if (is_selected) &
                call igSetItemDefaultFocus()
          end do
       end if
       call igEndCombo()
    end if
    call iw_tooltip("Set the reference field (input commands are applied to it)",ttshown)
    call igEndGroup()

    ! right-hand-side of lines 2 and 3: RUN button
    ok = iw_button("RUN",danger=.true.,sameline=.true.,siz=(/sz%x,sz%y/),&
       popupcontext=ok2,popupflags=ImGuiPopupFlags_MouseButtonRight)
    ok = ok .or. (igIsWindowFocused(ImGuiFocusedFlags_None) .and. is_bind_event(BIND_INPCON_RUN))
    if (ok) then
       idx = index(inputb,c_null_char)
       if (idx > 1) then
          if (associated(sy)) force_run_commands = 1
       end if
    end if
    if (ok2) then
       if (iw_menuitem("Run on all systems")) &
          force_run_commands = 2
       call iw_tooltip("Run these commands on all loaded systems")
       call igEndPopup()
    end if
    call iw_tooltip("Run the commands (" // trim(get_bind_keyname(BIND_INPCON_RUN)) // "). "//newline//&
       "Right-click to run commands on all systems.",ttshown)

    ! calculate sizes and draw the multiline
    call igGetContentRegionAvail(sz)
    str1 = "##inputmultiline" // c_null_char
    ldum = igInputTextMultiline(c_loc(str1),c_loc(inputb),maxlib,sz,&
       ImGuiInputTextFlags_AllowTabInput,c_null_funptr,c_null_ptr)

  contains
    subroutine assign_system(idx)
      use interfaces_glfw, only: glfwGetTime
      integer, intent(in) :: idx

      w%inpcon_selected = idx
      w%timelast_inpcon_assign = glfwGetTime()

    end subroutine assign_system

  end subroutine draw_ci

  !> Run the commands from the console input
  module subroutine run_commands_ci(w)
    use systemmod, only: sy
    use systems, only: launch_initialization_thread, kill_initialization_thread, are_threads_running,&
       sysc, sys_init, nsys, sys, lastchange_geometry
    use global, only: critic_main
    use tools_io, only: falloc, uin, fclose, ferror, faterr
    use iso_fortran_env, only: input_unit
    class(window), intent(inout), target :: w

    integer :: idx
    integer :: ios
    logical :: reinit, ldum

    ! if no system selected, return
    if (w%inpcon_selected < 1 .or. w%inpcon_selected > nsys) return
    if (sysc(w%inpcon_selected)%status /= sys_init) return

    ! check we have some input
    idx = index(inputb,c_null_char)
    if (idx <= 1) return

    ! connect the system
    sy => sys(w%inpcon_selected)

    ! if the initialization is happening, stop it
    reinit = are_threads_running()
    if (reinit) call kill_initialization_thread()

    ! connect a scratch file to uin, write the commands, rewind, and run
    uin = falloc()
    open(unit=uin,status='scratch',form='formatted',access='stream',iostat=ios)
    if (ios /= 0) &
       call ferror("run_commands","cannot open buffer for critic2 input",faterr)
    write (uin,'(A)') inputb(1:idx-1)
    rewind(uin)
    call critic_main()

    ! read the output
    ldum = read_output_uout(.true.)

    ! reinitialize the threads
    if (reinit) call launch_initialization_thread()

    ! set the time to rebuild lists
    call sysc(w%inpcon_selected)%set_timelastchange(lastchange_geometry)

    ! clean up
    call fclose(uin)
    uin = input_unit

  end subroutine run_commands_ci

  !> Block the GUI by dimming the background and showing the current
  !> console input. Useful for preparing the screen for running the
  !> commands in the console input. If allsys, the commands will apply
  !> to all loaded systems.
  module subroutine block_gui_ci(w,allsys)
    use gui_main, only: mainvwp, io, g, ColorWaitBg
    use utils, only: iw_text
    use param, only: newline
    class(window), intent(inout), target :: w
    logical, intent(in) :: allsys

    integer(c_int) :: flags
    type(ImVec2) :: sz, pivot
    character(kind=c_char,len=:), allocatable, target :: str1, text
    character(kind=c_char,len=:), allocatable, target :: csystem, cfield
    logical(c_bool) :: ldum

    !! blank the background
    flags = ImGuiWindowFlags_NoDecoration
    flags = ior(flags,ImGuiWindowFlags_NoMove)
    flags = ior(flags,ImGuiWindowFlags_NoSavedSettings)
    sz%x = 0._c_float
    sz%y = 0._c_float
    call igSetNextWindowPos(mainvwp%WorkPos,0,sz)
    call igSetNextWindowSize(mainvwp%WorkSize,0)
    call igPushStyleColor_Vec4(ImGuiCol_WindowBg,ColorWaitBg)
    str1 = "##blankbackground" // c_null_char
    ldum = .true.
    call igSetNextWindowFocus()
    ldum = igBegin(c_loc(str1), ldum, flags)
    call igEnd()
    call igPopStyleColor(1)

    !! overlay
    ! set window position at the center
    sz%x = io%DisplaySize%x * 0.5_c_float
    sz%y = io%DisplaySize%y * 0.5_c_float
    pivot%x = 0.5_c_float
    pivot%y = 0.5_c_float
    call igSetNextWindowPos(sz,0,pivot)

    ! get the input details
    call w%get_input_details_ci(csystem,cfield)

    ! set window size
    if (allsys) then
       text = "...Running critic2 input..." // newline //&
          "System: <all systems>" // newline //&
          "Field:  <all reference fields>" // newline //&
          "Input:  " // newline // inputb
    else
       text = "...Running critic2 input..." // newline //&
          "System: " // csystem // newline //&
          "Field:  " // cfield // newline //&
          "Input:  " // newline // inputb
    end if
    call igCalcTextSize(sz,c_loc(text),c_null_ptr,.false._c_bool,-1._c_float)
    sz%y = sz%y + 2 * g%Style%WindowPadding%y
    sz%x = sz%x + 2 * g%Style%WindowPadding%x
    call igSetNextWindowSize(sz,0)

    ! draw the window
    flags = ImGuiWindowFlags_NoDecoration
    flags = ior(flags,ImGuiWindowFlags_NoDocking)
    flags = ior(flags,ImGuiWindowFlags_AlwaysAutoResize)
    flags = ior(flags,ImGuiWindowFlags_NoSavedSettings)
    flags = ior(flags,ImGuiWindowFlags_NoFocusOnAppearing)
    flags = ior(flags,ImGuiWindowFlags_NoNav)
    ldum = .true.
    str1 = "##popupwait" // c_null_char
    call igSetNextWindowFocus()
    if (igBegin(c_loc(str1), ldum, flags)) then
       sz%x = 0
       sz%y = 0
       call igPushStyleVar_Vec2(ImGuiStyleVar_ItemSpacing,sz)
       call iw_text("...Running critic2 input...",highlight=.true.)
       call iw_text("System: ",highlight=.true.)
       if (allsys) then
          call iw_text("<all systems>",sameline=.true.)
       else
          call iw_text(csystem,sameline=.true.)
       end if
       call iw_text("Field:  ",highlight=.true.)
       if (allsys) then
          call iw_text("<all reference fields>",sameline=.true.)
       else
          call iw_text(cfield,sameline=.true.)
       end if
       call iw_text("Input:  ",highlight=.true.)
       call igIndent(0._c_float)
       call igText(c_loc(inputb))
       call igPopStyleVar(1_c_int)
       call igUnindent(0._c_float)
    end if
    call igEnd()

  end subroutine block_gui_ci

  !> Fill the input buffer with the given string.
  module subroutine fill_input_ci(w,str)
    class(window), intent(inout), target :: w
    character(len=*), intent(in) :: str

    integer :: idx

    idx = index(str,c_null_char)
    if (idx == 0) then
       idx = len_trim(str)
    else
       idx = idx - 1
    end if
    inputb(1:idx+1) = trim(str(1:idx)) // c_null_char

  end subroutine fill_input_ci

  !> Get the system and field strings for current input (without null char).
  module subroutine get_input_details_ci(w,csystem,cfield)
    use systems, only: nsys, sysc, sys_init, sys
    use tools_io, only: string
    class(window), intent(inout), target :: w
    character(len=:), allocatable, intent(inout) :: csystem, cfield

    integer :: iref

    ! system name and field name
    csystem = "<unknown>"
    cfield = "<unknown>"
    if (w%inpcon_selected >= 1 .and. w%inpcon_selected <= nsys) then
       if (sysc(w%inpcon_selected)%status == sys_init) then
          csystem = "(" // string(w%inpcon_selected) // ") " // trim(sysc(w%inpcon_selected)%seed%name)
          iref = sys(w%inpcon_selected)%iref
          cfield = "(" // string(iref) // ") " // trim(sys(w%inpcon_selected)%f(iref)%name)
       end if
    end if

  end subroutine get_input_details_ci

end submodule ci
