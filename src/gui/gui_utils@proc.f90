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

! Some utilities for building the GUI (e.g. wrappers around ImGui routines).
submodule (gui_utils) proc
  use iso_c_binding
  implicit none

contains

  !> Calculate the height of nline text lines and npadline padded lines:
  !>   frame-line1-frame-itemspace-frame-line2-frame-windowpad
  !> If endpad, add the end window padding
  module function iw_calcheight(npadline,nline,endpad)
    use gui_interfaces_cimgui
    use gui_main, only: g
    integer, intent(in) :: npadline
    integer, intent(in) :: nline
    logical, intent(in), optional :: endpad
    real(c_float) :: iw_calcheight

    logical :: endpad_

    endpad_ = .false.
    if (present(endpad)) endpad_ = endpad

    iw_calcheight = nline * igGetTextLineHeight() + &
       npadline * (igGetTextLineHeight() + 2 * g%Style%FramePadding%y) + &
       (npadline+nline - 1) * g%Style%ItemSpacing%y
    if (endpad_) iw_calcheight = iw_calcheight + g%Style%WindowPadding%y

  end function iw_calcheight

  !> Calculate the width of ntext characters and nbutton buttons.
  !>
  !> SomeText    |  Button1  |    |  Button2  |  ||--end of window
  !>             ^--^     ^--^    ^--^     ^--^    <-- FramePadding.x
  !>        ^----^           ^----^                <-- ItemSpacing.x
  !>                                          ^--^ <-- WindowPadding.x
  !>
  !> If from_end, return the position to place the cursor at the
  !> calculated distance from the end of the window, or the current
  !> position if it is negative.
  module function iw_calcwidth(ntext,nbutton,from_end)
    use gui_interfaces_cimgui
    use gui_main, only: g
    integer, intent(in) :: ntext
    integer, intent(in) :: nbutton
    logical, intent(in), optional :: from_end
    real(c_float) :: iw_calcwidth

    type(ImVec2) :: sz
    character(len=:,kind=c_char), allocatable, target :: strc

    ! text size
    allocate(character(len=ntext+1,kind=c_char) :: strc)
    strc(1:ntext) = ""
    strc(ntext+1:ntext+1) = c_null_char
    call igCalcTextSize(sz,c_loc(strc),c_null_ptr,.false._c_bool,-1._c_float)

    ! calculate width
    iw_calcwidth = sz%x
    if (nbutton > 0) then
       iw_calcwidth = iw_calcwidth + nbutton * (2 * g%Style%FramePadding%x) + &
          (nbutton - 1) * g%Style%ItemSpacing%x
    end if
    if (present(from_end)) then
       if (from_end) then
          iw_calcwidth = igGetWindowWidth() - g%style%WindowPadding%x - iw_calcwidth
       end if
    end if

  end function iw_calcwidth

  !> Draw text. If highlight, use the highlight color. If disabled,
  !> use the disabled font. If sameline, draw it in the same line as
  !> the previous widget.
  module subroutine iw_text(str,highlight,disabled,sameline)
    use gui_interfaces_cimgui
    use gui_main, only: ColorHighlightText
    character(len=*,kind=c_char), intent(in) :: str
    logical, intent(in), optional :: highlight
    logical, intent(in), optional :: sameline
    logical, intent(in), optional :: disabled

    character(len=:,kind=c_char), allocatable, target :: str1

    logical :: highlight_, disabled_, sameline_

    highlight_ = .false.
    sameline_ = .false.
    disabled_ = .false.
    if (present(highlight)) highlight_ = highlight
    if (present(sameline)) sameline_ = sameline
    if (present(disabled)) disabled_ = disabled

    if (sameline_) call igSameLine(0._c_float,-1._c_float)
    str1 = str // c_null_char
    if (disabled_) then
       call igTextDisabled(c_loc(str1))
    elseif (highlight_) then
       call igTextColored(ColorHighlightText,c_loc(str1))
    else
       call igText(c_loc(str1))
    end if

  end subroutine iw_text

  !> Draw a button. If danger, use the danger color. If sameline, draw
  !> the button in the same line as the preceding widgets.  If
  !> disabled, disable the button. If siz, use this size for the
  !> button. If popupcontext and poupflags, open a popup context with
  !> the given flags and return the resulting bool in popupcontext.
  module function iw_button(str,danger,sameline,disabled,siz,popupcontext,popupflags)
    use gui_interfaces_cimgui
    use gui_main, only: ColorDangerButton
    character(len=*,kind=c_char), intent(in) :: str
    logical, intent(in), optional :: danger
    logical, intent(in), optional :: sameline
    logical, intent(in), optional :: disabled
    real(c_float), intent(in), optional :: siz(2)
    logical, intent(inout), optional :: popupcontext
    integer(c_int), intent(in), optional :: popupflags
    logical :: iw_button

    character(len=:,kind=c_char), allocatable, target :: str1
    logical :: danger_, sameline_, disabled_
    type(ImVec2) :: sz

    if (present(siz)) then
       sz%x = siz(1)
       sz%y = siz(2)
    else
       sz%x = 0._c_float
       sz%y = 0._c_float
    end if
    danger_ = .false.
    sameline_ = .false.
    disabled_ = .false.
    if (present(danger)) danger_ = danger
    if (present(sameline)) sameline_ = sameline
    if (present(disabled)) disabled_ = disabled

    if (sameline_) call igSameLine(0._c_float,-1._c_float)
    str1 = trim(str) // c_null_char
    if (danger_) &
       call igPushStyleColor_Vec4(ImGuiCol_Button,ColorDangerButton)
    call igBeginDisabled(logical(disabled_,c_bool))
    iw_button = logical(igButton(c_loc(str1),sz))
    call igEndDisabled()
    if (danger_) &
       call igPopStyleColor(1)
    if (present(popupcontext) .and. present(popupflags)) &
       popupcontext = igBeginPopupContextItem(c_loc(str1),popupflags)

  end function iw_button

  !> Create a wrapped tooltip, maybe with a delay to show. ttshown
  !> activates the delay, and is the show flag for the delayed tooltip.
  module subroutine iw_tooltip(str,ttshown)
    use gui_interfaces_cimgui
    use gui_main, only: tooltip_wrap_factor, tooltip_delay
    character(len=*,kind=c_char), intent(in) :: str
    logical, intent(inout), optional :: ttshown

    character(len=:,kind=c_char), allocatable, target :: strloc

    if (present(ttshown)) then
       if (igIsItemHovered_delayed(ImGuiHoveredFlags_None,tooltip_delay,ttshown)) call show_tooltip()
    else
       if (igIsItemHovered(ImGuiHoveredFlags_None)) call show_tooltip()
    end if

  contains
    subroutine show_tooltip()
      strloc = trim(str) // c_null_char
      call igBeginTooltip()
      call igPushTextWrapPos(tooltip_wrap_factor * igGetFontSize())
      call igTextWrapped(c_loc(strloc))
      call igPopTextWrapPos()
      call igEndTooltip()
    end subroutine show_tooltip
  end subroutine iw_tooltip

  ! Returns true if the last item has been hovered for at least thr
  ! seconds. If already_shown (the tooltip has already been displayed),
  ! do not use the delay.
  module function igIsItemHovered_delayed(flags,thr,already_shown)
    use gui_main, only: g
    use gui_interfaces_cimgui, only: igIsItemHovered
    integer(c_int), value :: flags
    real(c_float), intent(in) :: thr
       logical, intent(in) :: already_shown
    logical(c_bool) :: igIsItemHovered_delayed

    igIsItemHovered_delayed = igIsItemHovered(flags) .and. (already_shown .or. g%HoveredIdTimer >= thr)

  end function igIsItemHovered_delayed

  ! Get a date/time string in a format adequate for the GUI
  module function get_time_string() result(output)
    use tools_io, only: string
    character(len=:), allocatable :: output

    integer :: values(8)

    call date_and_time(values=values)

    output = string(values(5),2,pad0=.true.) // ":" // string(values(6),2,pad0=.true.) //&
       ":" // string(values(7),2,pad0=.true.) // ", " // string(values(1)) // "/" // &
       string(values(2)) // "/" // string(values(3))

  end function get_time_string

end submodule proc
