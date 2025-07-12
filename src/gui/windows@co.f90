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

! Windows, console output.
submodule (windows) co
  use interfaces_cimgui
  implicit none

contains

  !> Update tasks for the rebond window, before the window is
  !> created.
  module subroutine update_co(w)
    use tools_io, only: string
    class(window), intent(inout), target :: w

    if (w%timelast_outcon_focused < timelast_output_written) then
       w%flags = ImGuiWindowFlags_UnsavedDocument
    else
       w%flags = ImGuiWindowFlags_None
    end if

  end subroutine update_co

  !> Draw the contents of the output console
  module subroutine draw_co(w)
    use interfaces_glfw, only: glfwGetTime
    use commands, only: com, clear_all_commands, nicom, icom, idcom
    use windows, only: stack_create_window
    use gui_main, only: g, ColorDangerButton, ColorFrameBgAlt
    use utils, only: igIsItemHovered_delayed, iw_tooltip, iw_button, iw_text,&
       iw_setposx_fromend, iw_calcheight, iw_calcwidth, iw_menuitem
    use tools_io, only: string
    class(window), intent(inout), target :: w

    integer :: i, curline, ndrawn, idx, idum
    character(kind=c_char,len=:), allocatable, target :: str1
    type(ImVec2) :: sz, szero, szavail
    logical(c_bool) :: ldum
    logical :: setscroll, skip, pushed, ok
    real(c_float) :: itemspacing, xavail, xavail1, xx, ally
    integer :: navail, navail1

    real(c_float), save :: maxallscrolly = 0._c_float ! max scroll value of the All pane
    real(c_float), save :: allscrolly = 0._c_float ! scroll value of the All pane
    logical, save :: ttshown = .false. ! tooltip flag

    ! initialize
    szero%x = 0._c_float
    szero%y = 0._c_float
    setscroll = .false.
    if (w%firstpass) &
       w%timelast_outcon_focused = glfwGetTime()

    ! if focused, update the time
    if (w%focused()) &
       w%timelast_outcon_focused = glfwGetTime()

    ! read new output, if available
    ldum = read_output_uout(.false.)

    ! get the available x
    call igGetContentRegionAvail(szavail)
    xavail = szavail%x

    ! first line: text
    call iw_text("Output",highlight=.true.)

    ! first line: clear button
    if (idcom == 0) then
       if (iw_button("Clear",sameline=.true.)) then
          outputb(1:1) = c_null_char
          timelast_output_written = glfwGetTime()
          lob = 0
       end if
       call iw_tooltip("Clear the output log",ttshown)
    end if

    ! first line: copy button
    if (iw_button("Copy",sameline=.true.)) then
       if (idcom == 0) then
          call igSetClipboardText(c_loc(outputb))
       else
          call igSetClipboardText(c_loc(com(icom(idcom))%output))
       end if
    end if
    call iw_tooltip("Copy the active output log to clipboard",ttshown)

    ! first line: save button
    if (iw_button("Save",sameline=.true.)) &
       idum = stack_create_window(wintype_dialog,.true.,wpurp_dialog_savelogfile,orraise=-1)
    call iw_tooltip("Save the active output log to a file",ttshown)

    ! first line: remove all button
    call iw_setposx_fromend(10,1)

    if (iw_button("Remove All",danger=.true.)) &
       call clear_all_commands()
    call iw_tooltip("Remove the output logs from all commands",ttshown)

    ! second line: all button
    xavail1 = max(xavail - iw_calcwidth(3,1) + g%Style%ItemSpacing%x,0._c_float)
    if (idcom == 0) then
       call igPushStyleColor_Vec4(ImGuiCol_Button,g%Style%Colors(ImGuiCol_ButtonActive+1))
    else
       call igPushStyleColor_Vec4(ImGuiCol_Button,ColorDangerButton)
    end if
    if (iw_button("All")) then
       idcom = 0
       setscroll = .true.
    end if
    call iw_tooltip("Show outputs from all commands in the console",ttshown)
    call igPopStyleColor(1)

    !! second line: list of command i/os
    ! make the itemspacing zero
    sz%x = 0
    sz%y = g%Style%ItemSpacing%y
    itemspacing = g%Style%ItemSpacing%x
    call igPushStyleVar_Vec2(ImGuiStyleVar_ItemSpacing,sz)

    ! calculate the button size and the number of buttons that fit
    xx = iw_calcwidth(ceiling(log10(max(maxval(com(icom(1:nicom))%id),1) + 0.1)),1)
    xavail = xavail / xx
    xavail1 = xavail1 / xx
    navail = max(floor(xavail),1)
    navail1 = max(floor(xavail1),1)

    ! render the buttons
    curline = 1
    ndrawn = 0
    do i = nicom, 1,  -1
       ! skip if no more buttons can be shown
       if (curline == 1) then
          skip = (ndrawn >= navail1)
       else
          skip = (ndrawn >= navail)
       end if

       ! skip or sameline
       if (skip) then
          ndrawn = 0
          curline = curline + 1
       else
          if (i == nicom) then
             call igSameLine(0._c_float,itemspacing)
          else
             call igSameLine(0._c_float,-1._c_float)
          end if
       end if

       ! render the button in alternate colors
       pushed = .true.
       if (idcom == i) then
          call igPushStyleColor_Vec4(ImGuiCol_Button,g%Style%Colors(ImGuiCol_ButtonActive+1))
       elseif (mod(i,2) == 0) then
          call igPushStyleColor_Vec4(ImGuiCol_Button,ColorFrameBgAlt)
       else
          pushed = .false.
       end if
       if (iw_button(string(com(icom(i))%id),popupcontext=ok,popupflags=ImGuiPopupFlags_MouseButtonRight,&
          siz=(/xx,iw_calcheight(1,0)/))) then
          idcom = i
          setscroll = .true.
       end if
       if (pushed) call igPopStyleColor(1)

       ! tooltip
       call iw_tooltip(com(icom(i))%tooltipinfo)

       ! context menu
       if (ok) then
          if (iw_menuitem("Edit Input")) then
             idx = index(com(icom(i))%input,c_null_char)
             if (idx > 0) then
                call w%fill_input_ci(com(icom(i))%input(1:idx))
             end if
             call igFocusWindow(win(iwin_console_input)%ptr)
          end if

          if (iw_menuitem("Remove")) then
             call com(icom(i))%end()
             icom(i:nicom-1) = icom(i+1:nicom)
             nicom = nicom - 1
             if (idcom == i) then
                idcom = 0
             elseif (idcom > i) then
                idcom = idcom - 1
             end if
          end if

          call igEndPopup()
       end if

       ndrawn = ndrawn + 1
    end do
    call igPopStyleVar(1_c_int)

    ! calculate sizes and draw the multiline (with dark background and border)
    call igGetContentRegionAvail(sz)
    call igPushStyleColor_Vec4(ImGuiCol_FrameBg,g%Style%Colors(ImGuiCol_WindowBg+1))
    call igPushStyleVar_Float(ImGuiStyleVar_FrameBorderSize,1._c_float)
    str1 = "##outputmultiline" // c_null_char
    if (idcom == 0) then
       ldum = igInputTextMultiline(c_loc(str1),c_loc(outputb),lob,sz,&
          ImGuiInputTextFlags_ReadOnly,c_null_funptr,c_null_ptr)
    else
       ldum = igInputTextMultiline(c_loc(str1),c_loc(com(icom(idcom))%output),&
          com(icom(idcom))%size+1,sz,ImGuiInputTextFlags_ReadOnly,&
          c_null_funptr,c_null_ptr)
    end if
    call igPopStyleVar(1)
    call igPopStyleColor(1)

    ! Set scroll to previous value if changing output, the All view scrolls to the
    ! end on new output
    ldum = igBeginChild_Str(c_loc(str1),szero,.false._c_bool,ImGuiWindowFlags_None)
    if (setscroll) then
       if (idcom > 0) then
          call igSetScrollY_Float(com(icom(idcom))%scrolly)
       else
          call igSetScrollY_Float(allscrolly)
       end if
    endif

    ! if there is new output in all, set the scroll to the beginning
    ! of the new output
    if (idcom == 0) then
       ally = igGetScrollMaxY()
       if (abs(ally - maxallscrolly) > 1d-5) then
          if (maxallscrolly > 0._c_float) then
             call igGetItemRectSize(sz)
             allscrolly = min(maxallscrolly + sz%y - 2 * igGetTextLineHeight(),ally)
             call igSetScrollY_Float(allscrolly)
          end if
          maxallscrolly = ally
       end if
    end if

    ! write down the current scroll y
    if (idcom > 0) then
       com(icom(idcom))%scrolly = igGetScrollY()
    else
       allscrolly = igGetScrollY()
    end if
    call igEndChild()

  end subroutine draw_co

  !> Read new output from the scratch LU uout and append it to the
  !> output buffer.  If iscom, this output corresponds to a command,
  !> so create the command I/O object. If cominfo, use the provided
  !> string as command info instead of the contents of the input buffer.
  !> Return true if output has been read.
  module function read_output_uout(iscom,cominfo)
    use interfaces_glfw, only: glfwGetTime
    use commands, only: command_inout, com, icom, command_inout_empty, command_inout_used,&
       maxcomout, ncom, ncomid, nicom
    use utils, only: get_time_string
    use gui_main, only: are_threads_running
    use tools_io, only: uout, getline_raw, string, ferror, faterr
    use types, only: realloc
    use param, only: newline
    logical, intent(in) :: iscom
    character*(*), intent(in), optional :: cominfo
    logical :: read_output_uout

    character(kind=c_char,len=:), allocatable, target :: csystem, cfield
    type(command_inout), allocatable :: aux(:)
    character(len=:), allocatable :: line, commonstr, showinp
    integer(c_size_t) :: pos, lshift, ll, olob, total, newsize
    integer :: idx, ithis, i
    logical :: ok

    ! allocate the output buffer if not allocated
    read_output_uout = .false.
    if (.not.allocated(outputb)) then
       allocate(character(len=maxlob+1) :: outputb)
       outputb(1:1) = c_null_char
       lob = 0
    end if

    ! allocate the command input/output stack, if not allocated
    if (.not.allocated(com)) then
       ncom = 0
       allocate(com(10))
    end if
    if (.not.allocated(icom)) then
       nicom = 0
       allocate(icom(10))
    end if

    ! do not read if initialization threads are still running (may generate output)
    if (are_threads_running()) return

    ! get size of the new output
    inquire(uout,pos=pos)
    if (pos > 1) then
       ! there is new output, rewind, read it, and rewind again

       ! I am going to need pos-1 new characters
       ll = pos-1
       rewind(uout)

       ! do not have enough room to accomodate this much output = discard new text
       do while (ll > maxlob)
          ok = getline_raw(uout,line)
          if (.not.ok) return
          ll = ll - (len(line)+1)
       end do

       ! check whether we can accomodate the new text = discard from the beginning
       if (lob + ll > maxlob) then
          lshift = ll + lob - maxlob
          outputb(1:lob-lshift) = outputb(lshift+1:lob)
          lob = lob - lshift
          ! adjust to the next newline
          idx = index(outputb,newline)
          if (idx == 0) then
             lob = 0
          else
             outputb(1:lob-idx) = outputb(idx+1:lob)
             lob = lob - idx
          end if
       end if

       ! read the new output and rewind
       olob = lob
       do while(getline_raw(uout,line))
          ll = len(line)
          if (ll > 0) then
             outputb(lob+1:lob+ll) = line(1:ll)
             lob = lob + ll
          end if
          outputb(lob+1:lob+1) = newline
          lob = lob + 1
          inquire(uout,pos=pos)
       end do
       outputb(lob+1:lob+1) = c_null_char
       rewind(uout)

       if (iscom) then
          ! try to remove commands if we exceed the size
          newsize = lob - olob + 1
          if (newsize > maxcomout) &
             call ferror('read_output_ci','exceeded output buffer for command',faterr)
          total = newsize
          do i = 1, nicom
             total = total + com(icom(i))%size
          end do
          i = 0
          do while (total > maxcomout)
             i = i + 1
             total = total - com(icom(i))%size
             call com(icom(i))%end()
          end do
          icom(1:nicom-i) = icom(i+1:nicom)
          nicom = nicom - i

          ! get an unusued command ID
          ithis = 0
          do i = 1, ncom
             if (com(i)%status == command_inout_empty) then
                ithis = i
                exit
             end if
          end do

          ! add a new command, reallocate if necessary
          if (ithis == 0) then
             ncom = ncom + 1
             if (ncom > size(com,1)) then
                allocate(aux(2*ncom))
                aux(1:size(com,1)) = com
                call move_alloc(aux,com)
             end if
             ithis = ncom
          end if

          ! fill the new command info
          ncomid = ncomid + 1
          call win(iwin_console_input)%get_input_details_ci(csystem,cfield)
          com(ithis)%id = ncomid
          com(ithis)%status = command_inout_used
          com(ithis)%size = lob - olob + 1
          if (present(cominfo)) then
             showinp = trim(cominfo)
             com(ithis)%input = c_null_char
          else
             idx = index(inputb,c_null_char)
             showinp = inputb(1:idx-1)
             com(ithis)%input = inputb(1:idx)
          end if

          commonstr = "### Command: " // string(ncomid) // " (" // get_time_string() // ")" // newline
          com(ithis)%tooltipinfo = commonstr // &
             "System: " // csystem // newline //&
             "Field: " // cfield // newline //&
             "Input: " // newline // showinp // newline //&
             "#########" // newline // newline // "[Right-click for options]"
          com(ithis)%output = commonstr // &
             "## System: " // csystem // newline //&
             "## Field: " // cfield // newline //&
             "## Input: " // newline // showinp // newline //&
             "#########" // newline // newline // outputb(olob+1:lob+1)

          ! add it to the list of active commands
          nicom = nicom + 1
          if (nicom > size(icom,1)) &
             call realloc(icom,2*nicom)
          icom(nicom) = ithis
       end if

       ! we have new data
       read_output_uout = .true.

       ! update the time
       timelast_output_written = glfwGetTime()
    end if

  end function read_output_uout

end submodule co
