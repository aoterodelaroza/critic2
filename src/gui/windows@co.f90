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

  !> Draw the contents of the output console
  module subroutine draw_co(w)
    use windows, only: stack_create_window
    use gui_main, only: g, ColorDangerButton, ColorFrameBgAlt
    use utils, only: igIsItemHovered_delayed, iw_tooltip, iw_button, iw_text,&
       iw_setposx_fromend, iw_calcheight, iw_calcwidth, iw_menuitem
    use tools_io, only: string
    class(window), intent(inout), target :: w

    integer :: i, curline, ndrawn, idx, idum
    character(kind=c_char,len=:), allocatable, target :: str1, strpop
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

    ! read new output, if available
    ldum = w%read_output_ci(.false.)

    ! get the available x
    call igGetContentRegionAvail(szavail)
    xavail = szavail%x

    ! first line: text
    call iw_text("Output",highlight=.true.)

    ! first line: clear button
    if (idcom == 0) then
       if (iw_button("Clear",sameline=.true.)) then
          outputb(1:1) = c_null_char
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

    if (iw_button("Remove All",danger=.true.)) then
       ! remove all command i/o information
       ncomid = 0
       ncom = 0
       nicom = 0
       idcom = 0
       icom = 0
       do i = 1, ncom
          call com(i)%end()
       end do
    end if
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

end submodule co
