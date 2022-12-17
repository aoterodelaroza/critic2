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
submodule (utils) proc
  use iso_c_binding
  implicit none

contains

  !> Set the cursor X position a distance from the end of the content
  !> region corresponding to ntext characters and nbutton buttons.
  module subroutine iw_setposx_fromend(ntext,nbutton)
    use interfaces_cimgui
    integer, intent(in) :: ntext
    integer, intent(in) :: nbutton

    real(c_float) :: posx

    call igSameLine(0._c_float,-1._c_float)
    posx = iw_calcwidth(ntext,nbutton,from_end=.true.)
    if (posx > igGetCursorPosX()) &
       call igSetCursorPosX(posx)

  end subroutine iw_setposx_fromend

  !> Calculate the height of nline text lines and npadline padded lines:
  !>   frame-line1-frame-itemspace-frame-line2-frame-windowpad
  !> If endpad, add the end window padding
  module function iw_calcheight(npadline,nline,endpad)
    use interfaces_cimgui
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
    use interfaces_cimgui
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
    iw_calcwidth = max(iw_calcwidth,0._c_float)

  end function iw_calcwidth

  !> Simple combo with title str. stropt contains the options
  !> separated by \0 and terminated by \0\0. ival is the current value
  !> of the combo.  sameline = place it in the same line as the last
  !> item. changed = returns true if the combo option changed.
  module subroutine iw_combo_simple(str,stropt,ival,sameline,changed)
    use interfaces_cimgui
    character(len=*,kind=c_char), intent(in) :: str
    character(len=*,kind=c_char), intent(in) :: stropt
    integer, intent(inout) :: ival
    logical, intent(in), optional :: sameline
    logical(c_bool), intent(out), optional :: changed

    character(len=:,kind=c_char), allocatable, target :: str1
    character(len=:,kind=c_char), allocatable, target :: stropt1
    logical :: ldum, sameline_
    integer :: ll, maxlen, nidx, idx

    sameline_ = .false.
    if (present(sameline)) sameline_ = sameline

    str1 = str // c_null_char
    stropt1 = stropt // c_null_char // c_null_char

    if (sameline_) &
       call igSameLine(0._c_float,-1._c_float)

    maxlen = 0
    idx = 0
    ll = len(stropt1)-1
    do while (idx < ll)
       nidx = idx + index(stropt1(idx+1:),c_null_char)
       maxlen = max(maxlen,nidx-1-idx)
       idx = nidx
    end do
    call igSetNextItemWidth(iw_calcwidth(maxlen+4,0))

    ldum = igCombo_Str(c_loc(str1), ival, c_loc(stropt1), -1_c_int)
    if (present(changed)) changed = ldum

  end subroutine iw_combo_simple

  !> Draw a radio button with title str. If bool and boolval are given
  !> associated with logical bool and with value boolval. If int and
  !> intval are presented, associated with integer int and with value
  !> intval
  module function iw_radiobutton(str,bool,boolval,int,intval,sameline)
    use interfaces_cimgui
    character(len=*,kind=c_char), intent(in) :: str
    logical, intent(inout), optional :: bool
    logical, intent(in), optional :: boolval
    integer(c_int), intent(inout), optional :: int
    integer(c_int), intent(in), optional :: intval
    logical, intent(in), optional :: sameline
    logical :: iw_radiobutton

    character(len=:,kind=c_char), allocatable, target :: str1
    logical :: sameline_

    iw_radiobutton = .false.
    sameline_ = .false.
    if (present(sameline)) sameline_ = sameline

    if (sameline_) &
       call igSameLine(0._c_float,-1._c_float)
    str1 = str // c_null_char
    if (present(bool).and.present(boolval)) then
       if (igRadioButton_Bool(c_loc(str1),logical(bool.eqv.boolval,c_bool))) then
          bool = boolval
          iw_radiobutton = .true.
       end if
    elseif (present(int).and.present(intval)) then
       if (igRadioButton_IntPtr(c_loc(str1),int,intval)) iw_radiobutton = .true.
    end if

  end function iw_radiobutton

  !> Draw text. If highlight, use the highlight color. If disabled,
  !> use the disabled font. If sameline, draw it in the same line as
  !> the previous widget. If sameline_nospace, draw it in the same
  !> line adjacent to the previous item. If noadvance, do not advance
  !> the cursor after writing. If copy_to_output, write the text to
  !> uout as well (without advancing to a new line and with a comma
  !> after the string). If centered, center the text in the window.
  module subroutine iw_text(str,highlight,danger,disabled,sameline,sameline_nospace,&
     noadvance,copy_to_output,centered)
    use interfaces_cimgui
    use gui_main, only: ColorHighlightText, ColorDangerText
    use tools_io, only: uout
    character(len=*,kind=c_char), intent(in) :: str
    logical, intent(in), optional :: highlight
    logical, intent(in), optional :: danger
    logical, intent(in), optional :: sameline
    logical, intent(in), optional :: sameline_nospace
    logical, intent(in), optional :: disabled
    logical, intent(in), optional :: noadvance
    logical, intent(in), optional :: copy_to_output
    logical, intent(in), optional :: centered

    character(len=:,kind=c_char), allocatable, target :: str1

    logical :: highlight_, danger_, disabled_, sameline_, sameline_nospace_
    logical :: noadvance_,copy_to_output_, centered_
    real(c_float) :: pos, wwidth, twidth
    type(ImVec2) :: sz

    highlight_ = .false.
    danger_ = .false.
    sameline_ = .false.
    sameline_nospace_ = .false.
    disabled_ = .false.
    noadvance_ = .false.
    copy_to_output_ = .false.
    centered_ = .false.
    if (present(highlight)) highlight_ = highlight
    if (present(danger)) danger_ = danger
    if (present(sameline)) sameline_ = sameline
    if (present(sameline_nospace)) sameline_nospace_ = sameline_nospace
    if (present(disabled)) disabled_ = disabled
    if (present(noadvance)) noadvance_ = noadvance
    if (present(copy_to_output)) copy_to_output_ = copy_to_output
    if (present(centered)) centered_ = centered

    str1 = str // c_null_char
    if (noadvance_) pos = igGetCursorPosX()
    if (sameline_) call igSameLine(0._c_float,-1._c_float)
    if (sameline_nospace_) call igSameLine(0._c_float,0._c_float)
    if (centered_) then
       wwidth = igGetWindowWidth()
       call igCalcTextSize(sz,c_loc(str1),c_null_ptr,.false._c_bool,-1.0_c_float)
       twidth = sz%x
       call igSetCursorPosX((wwidth - twidth) * 0.5_c_float)
    end if
    if (disabled_) then
       call igTextDisabled(c_loc(str1))
    elseif (highlight_) then
       call igTextColored(ColorHighlightText,c_loc(str1))
    elseif (danger_) then
       call igTextColored(ColorDangerText,c_loc(str1))
    else
       call igText(c_loc(str1))
    end if
    if (noadvance_) then
       call igSameLine(0._c_float,0._c_float)
       call igSetCursorPosX(pos)
    end if
    if (copy_to_output_) &
       write (uout,'(A,",")',advance='no') str

  end subroutine iw_text

  !> Draw a button. If danger, use the danger color. If sameline, draw
  !> the button in the same line as the preceding widgets.  If
  !> disabled, disable the button. If siz, use this size for the
  !> button. If popupcontext and poupflags, open a popup context with
  !> the given flags and return the resulting bool in popupcontext.
  module function iw_button(str,danger,sameline,disabled,siz,popupcontext,popupflags)
    use interfaces_cimgui
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
    use interfaces_cimgui
    use gui_main, only: tooltip_wrap_factor, tooltip_delay, tooltip_enabled, fontsize
    character(len=*,kind=c_char), intent(in) :: str
    logical, intent(inout), optional :: ttshown

    character(len=:,kind=c_char), allocatable, target :: strloc
    integer :: flags

    if (.not.tooltip_enabled) return
    flags = ImGuiHoveredFlags_AllowWhenBlockedByPopup
    if (present(ttshown)) then
       if (igIsItemHovered_delayed(flags,tooltip_delay,ttshown)) call show_tooltip()
    else
       if (igIsItemHovered(flags)) call show_tooltip()
    end if

  contains
    subroutine show_tooltip()
      strloc = trim(str) // c_null_char
      call igBeginTooltip()
      call igPushTextWrapPos(tooltip_wrap_factor * fontsize%x)
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
    use interfaces_cimgui, only: igIsItemHovered
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
       ":" // string(values(7),2,pad0=.true.) // " - " // string(values(1)) // "/" // &
       string(values(2)) // "/" // string(values(3))

  end function get_time_string

  !> Read the buffer buf line by line up to the null character. Write
  !> the read lines to logical unit LU
  module subroutine buffer_to_string_array(buf,lu,prefix,suffix)
    use param, only: newline
    character*(*), intent(in) :: buf
    integer, intent(in) :: lu
    character*(*), intent(in), optional :: prefix
    character*(*), intent(in), optional :: suffix

    integer :: idx
    character(kind=c_char,len=:), allocatable, target :: left, pre, suf
    logical :: exloop

    pre = ""
    suf = ""
    if (present(prefix)) pre = prefix
    if (present(suffix)) suf = suffix

    idx = 0
    left = buf(1:index(buf,c_null_char)-1)
    exloop = .false.
    do while (.not.exloop)
       idx = index(left,newline)
       if (idx == 0) then
          idx = len_trim(left) + 1
          exloop = .true.
       end if
       if (len_trim(left(:idx-1)) > 0) then
          write (lu,'(A,A,A)') pre, left(:idx-1), suf
       end if
       if(.not.exloop) left = left(idx+1:)
    end do

  end subroutine buffer_to_string_array

end submodule proc
