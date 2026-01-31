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

  !> Drag float button for 1, 2, 3, and 4 floating point numbers. The
  !> number of buttons shown is determined by whether the x1, x2, x3,
  !> or x4 argument is passed.  Speed = step for the dragfloat. Min
  !> and max = minimum and maximum values. scale = scale the numbers
  !> by this value before and after the drag. sformat = C format of the
  !> label. flags = combination of ImGuiSliderFlags_* flags. Version
  !> for real(c_float) type.
  module function iw_dragfloat_realc(str,x1,x2,x3,x4,speed,min,max,scale,sformat,flags)
    use interfaces_cimgui
    character(len=*,kind=c_char), intent(in) :: str
    real(c_float), intent(inout), optional :: x1
    real(c_float), intent(inout), optional :: x2(2)
    real(c_float), intent(inout), optional :: x3(3)
    real(c_float), intent(inout), optional :: x4(4)
    real(c_float), intent(in), optional :: speed, min, max, scale
    character(len=*,kind=c_char), intent(in), optional :: sformat
    integer(c_int), intent(in), optional :: flags
    logical :: iw_dragfloat_realc

    real(c_float) :: speed_, min_, max_, scale_
    character(len=:,kind=c_char), allocatable, target :: str_, sformat_
    integer(c_int) :: flags_
    real(c_float) :: x1_, x2_(2), x3_(3), x4_(4)

    str_ = trim(str) // c_null_char
    speed_ = 1._c_float
    if (present(speed)) speed_ = speed
    min_ = -FLT_MAX
    if (present(min)) min_ = min
    max_ = FLT_MAX
    if (present(max)) max_ = max
    scale_ = 1._c_float
    if (present(scale)) scale_ = scale
    sformat_ = "%.3f" // c_null_char
    if (present(sformat)) sformat_ = trim(sformat) // c_null_char
    flags_ = 0_c_int
    if (present(flags)) flags_ = flags

    if (present(x1)) then
       x1_ = x1 * scale_
       iw_dragfloat_realc = igDragFloat(c_loc(str_),x1_,speed_,min_,max_,c_loc(sformat_),flags_)
       if (iw_dragfloat_realc) x1 = x1_ / scale_
    elseif (present(x2)) then
       x2_ = x2 * scale_
       iw_dragfloat_realc = igDragFloat2(c_loc(str_),x2_,speed_,min_,max_,c_loc(sformat_),flags_)
       if (iw_dragfloat_realc) x2 = x2_ / scale_
    elseif (present(x3)) then
       x3_ = x3 * scale_
       iw_dragfloat_realc = igDragFloat3(c_loc(str_),x3_,speed_,min_,max_,c_loc(sformat_),flags_)
       if (iw_dragfloat_realc) x3 = x3_ / scale_
    elseif (present(x4)) then
       x4_ = x4 * scale_
       iw_dragfloat_realc = igDragFloat4(c_loc(str_),x4_,speed_,min_,max_,c_loc(sformat_),flags_)
       if (iw_dragfloat_realc) x4 = x4_ / scale_
    else
       iw_dragfloat_realc = .false.
    end if

  end function iw_dragfloat_realc

  !> Drag float button for 1, 2, 3, and 4 floating point numbers. The
  !> number of buttons shown is determined by whether the x1, x2, x3,
  !> or x4 argument is passed.  Speed = step for the dragfloat. Min
  !> and max = minimum and maximum values. scale = scale the numbers
  !> by this value before and after the drag. sformat = C format of the
  !> label. flags = combination of ImGuiSliderFlags_* flags. Version
  !> for real*8 type.
  module function iw_dragfloat_real8(str,x1,x2,x3,x4,speed,min,max,scale,sformat,flags)
    use interfaces_cimgui
    character(len=*,kind=c_char), intent(in) :: str
    real*8, intent(inout), optional :: x1
    real*8, intent(inout), optional :: x2(2)
    real*8, intent(inout), optional :: x3(3)
    real*8, intent(inout), optional :: x4(4)
    real*8, intent(in), optional :: speed, min, max, scale
    character(len=*,kind=c_char), intent(in), optional :: sformat
    integer(c_int), intent(in), optional :: flags
    logical :: iw_dragfloat_real8

    real(c_float) :: speed_, min_, max_, scale_, x1_, x2_(2), x3_(3), x4_(4)
    character(len=:,kind=c_char), allocatable, target :: str_, sformat_
    integer(c_int) :: flags_

    str_ = trim(str) // c_null_char
    speed_ = 1._c_float
    if (present(speed)) speed_ = real(speed,c_float)
    min_ = -FLT_MAX
    if (present(min)) min_ = real(min,c_float)
    max_ = FLT_MAX
    if (present(max)) max_ = real(max,c_float)
    scale_ = 1._c_float
    if (present(scale)) scale_ = scale
    sformat_ = "%.3f" // c_null_char
    if (present(sformat)) sformat_ = trim(sformat) // c_null_char
    flags_ = 0_c_int
    if (present(flags)) flags_ = flags

    if (present(x1)) then
       x1_ = real(x1 * scale_,c_float)
       iw_dragfloat_real8 = igDragFloat(c_loc(str_),x1_,speed_,min_,max_,c_loc(sformat_),flags_)
       if (iw_dragfloat_real8) x1 = x1_ / scale_
    elseif (present(x2)) then
       x2_ = real(x2 * scale_,c_float)
       iw_dragfloat_real8 = igDragFloat2(c_loc(str_),x2_,speed_,min_,max_,c_loc(sformat_),flags_)
       if (iw_dragfloat_real8) x2 = x2_ / scale_
    elseif (present(x3)) then
       x3_ = real(x3 * scale_,c_float)
       iw_dragfloat_real8 = igDragFloat3(c_loc(str_),x3_,speed_,min_,max_,c_loc(sformat_),flags_)
       if (iw_dragfloat_real8) x3 = x3_ / scale_
    elseif (present(x4)) then
       x4_ = real(x4 * scale_,c_float)
       iw_dragfloat_real8 = igDragFloat4(c_loc(str_),x4_,speed_,min_,max_,c_loc(sformat_),flags_)
       if (iw_dragfloat_real8) x4 = x4_ / scale_
    else
       iw_dragfloat_real8 = .false.
    end if

  end function iw_dragfloat_real8

  !> Clamp a 3-color to the 0->1 interval
  module subroutine iw_clamp_color3(rgb)
    real(c_float), intent(inout) :: rgb(3)

    rgb = min(rgb,1._c_float)
    rgb = max(rgb,0._c_float)

  end subroutine iw_clamp_color3

  !> Clamp a 4-color to the 0->1 interval
  module subroutine iw_clamp_color4(rgba)
    real(c_float), intent(inout) :: rgba(4)

    rgba = min(rgba,1._c_float)
    rgba = max(rgba,0._c_float)

  end subroutine iw_clamp_color4

  !> Draw a dragcolor3 widget with title str and color rgb. If
  !> sameline, draw it in the same line as the last object. If
  !> nolabel, do not show the widget label. Returns true if the color
  !> changed.
  module function iw_coloredit(str,rgb,rgba,sameline,nolabel,nointeraction)
    use interfaces_cimgui
    character(len=*,kind=c_char), intent(in) :: str
    real(c_float), intent(inout), optional :: rgb(3)
    real(c_float), intent(inout), optional :: rgba(4)
    logical, intent(in), optional :: sameline, nolabel, nointeraction
    logical :: iw_coloredit

    character(len=:,kind=c_char), allocatable, target :: str1
    logical :: sameline_, nolabel_, nointeraction_
    integer(c_int) :: flags

    sameline_ = .false.
    if (present(sameline)) sameline_ = sameline
    nolabel_ = .false.
    if (present(nolabel)) nolabel_ = nolabel
    nointeraction_ = .false.
    if (present(nointeraction)) nointeraction_ = nointeraction

    if (sameline_) &
       call igSameLine(0._c_float,-1._c_float)
    str1 = str // c_null_char

    flags = ImGuiColorEditFlags_NoInputs
    if (nolabel_) flags = ImGuiColorEditFlags_NoLabel
    if (nointeraction_) then
       flags = ior(flags,ImGuiColorEditFlags_NoPicker)
       flags = ior(flags,ImGuiColorEditFlags_NoOptions)
       flags = ior(flags,ImGuiColorEditFlags_NoTooltip)
       flags = ior(flags,ImGuiColorEditFlags_NoSidePreview)
       flags = ior(flags,ImGuiColorEditFlags_NoDragDrop)
    end if
    if (present(rgb)) then
       iw_coloredit = igColorEdit3(c_loc(str1),rgb,flags)
       call iw_clamp_color3(rgb)
    elseif (present(rgba)) then
       iw_coloredit = igColorEdit4(c_loc(str1),rgba,flags)
       call iw_clamp_color4(rgba)
    end if

  end function iw_coloredit

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
  !> separated by \0 and terminated by \0. ival is the current value
  !> of the combo. sameline = place it in the same line as the last
  !> item. samline_nospace = like sameline, but no extra
  !> space. changed = returns true if the combo option
  !> changed. noarrow = hide the arrow on the right of the combo.
  module subroutine iw_combo_simple(str,stropt,ival,sameline,sameline_nospace,changed,noarrow)
    use interfaces_cimgui
    use types, only: realloc
    character(len=*,kind=c_char), intent(in) :: str
    character(len=*,kind=c_char), intent(in) :: stropt
    integer, intent(inout) :: ival
    logical, intent(in), optional :: sameline
    logical, intent(in), optional :: sameline_nospace
    logical(c_bool), intent(out), optional :: changed
    logical, intent(in), optional :: noarrow

    type(ImVec2) :: szero
    character(len=:,kind=c_char), allocatable, target :: str1, str2, preview
    character(len=:,kind=c_char), allocatable, target :: stropt1
    logical :: sameline_, sameline_nospace_, noarrow_
    integer :: ll, maxlen, nword, i, iselect
    integer(c_int) :: flags
    integer, allocatable :: idx(:)
    logical(c_bool) :: selected

    ! process input options
    noarrow_ = .false.
    if (present(noarrow)) noarrow_ = noarrow
    sameline_ = .false.
    if (present(sameline)) sameline_ = sameline
    sameline_nospace_ = .false.
    if (present(sameline_nospace)) sameline_nospace_ = sameline_nospace
    szero%x = 0
    szero%y = 0

    ! strip null chars from the options string
    str1 = str // c_null_char
    stropt1 = stropt
    do while (.true.)
       ll = len(stropt1)
       if (index(stropt1,c_null_char,back=.true.) == ll) then
          stropt1 = stropt1(1:ll-1)
       else
          exit
       end if
    end do
    stropt1 = stropt1 // c_null_char

    ! number of words and indices of the nulls
    allocate(idx(20))
    preview = c_null_char
    nword = 0
    idx(1) = 0
    maxlen = 0
    do i = 1, len(stropt1)
       if (stropt1(i:i) == c_null_char) then
          nword = nword + 1
          if (nword+1 > size(idx,1)) call realloc(idx,2*nword)
          idx(nword+1) = i
          if (ival+1 == nword) preview = stropt1(idx(nword)+1:i-1) // c_null_char
          maxlen = max(maxlen,i-1-idx(nword))
       end if
    end do

    ! same line
    if (sameline_) call igSameLine(0._c_float,-1._c_float)
    if (sameline_nospace_) call igSameLine(0._c_float,0._c_float)

    ! display the combo
    if (noarrow_) then
       call igSetNextItemWidth(iw_calcwidth(maxlen+1,0))
    else
       call igSetNextItemWidth(iw_calcwidth(maxlen+4,0))
    end if
    iselect = ival
    flags = ImGuiComboFlags_None
    if (noarrow_) flags = ior(flags,ImGuiComboFlags_NoArrowButton)
    if (igBeginCombo(c_loc(str1),c_loc(preview),flags)) then
       do i = 1, nword
          str2 = stropt1(idx(i)+1:idx(i+1)-1) // c_null_char
          selected = (i == ival+1)
          if (igSelectable_Bool(c_loc(str2),selected,ImGuiSelectableFlags_None,szero)) &
             iselect = i-1
          if (selected) &
             call igSetItemDefaultFocus()
       end do
       call igEndCombo()
    end if
    if (present(changed)) changed = (ival /= iselect)
    ival = iselect

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

  !> Draw a checkbox with title str. The value of the checkbox is
  !> associated with bool. sameline = draw it in the same line as the
  !> previous widget.
  module function iw_checkbox(str,bool,sameline,highlight)
    use interfaces_cimgui
    use gui_main, only: ColorHighlightText
    character(len=*,kind=c_char), intent(in) :: str
    logical, intent(inout) :: bool
    logical, intent(in), optional :: sameline
    logical, intent(in), optional :: highlight
    logical :: iw_checkbox

    logical :: sameline_
    character(len=:,kind=c_char), allocatable, target :: str1
    logical(c_bool) :: bool_
    logical :: highlight_

    iw_checkbox = .false.
    bool_ = logical(bool,c_bool)
    sameline_ = .false.
    highlight_ = .false.
    if (present(sameline)) sameline_ = sameline
    if (present(highlight)) highlight_ = highlight

    if (sameline_) &
       call igSameLine(0._c_float,-1._c_float)
    if (highlight_) &
       call igPushStyleColor_Vec4(ImGuiCol_Text,ColorHighlightText)

    str1 = str // c_null_char
    iw_checkbox = igCheckbox(c_loc(str1),bool_)
    bool = logical(bool_)

    if (highlight_) &
       call igPopStyleColor(1)

  end function iw_checkbox

  !> Draw text. highlight = use the highlight color. danger = use the
  !> danger color. disabled = use the disabled font. sameline = draw
  !> it in the same line as the previous widget. sameline_nospace =
  !> draw it in the same line adjacent to the previous item. noadvance
  !> = do not advance the cursor after writing. copy_to_output = write
  !> the text to uout as well (without advancing to a new line and with
  !> a comma after the string). centered = center the text in the window.
  !> rgba = use this color for the text.
  module subroutine iw_text(str,highlight,danger,disabled,sameline,sameline_nospace,&
     noadvance,copy_to_output,centered,rgb,rgba)
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
    real(c_float), intent(in), optional :: rgb(3)
    real(c_float), intent(in), optional :: rgba(4)

    character(len=:,kind=c_char), allocatable, target :: str1

    logical :: highlight_, danger_, disabled_, sameline_, sameline_nospace_
    logical :: noadvance_,copy_to_output_, centered_
    real(c_float) :: pos, wwidth, twidth
    type(ImVec2) :: sz
    type(ImVec4) :: col

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
    elseif (present(rgb)) then
       col%x = rgb(1)
       col%y = rgb(2)
       col%z = rgb(3)
       col%w = 1._c_float
       call igTextColored(col,c_loc(str1))
    elseif (present(rgba)) then
       col%x = rgba(1)
       col%y = rgba(2)
       col%z = rgba(3)
       col%w = rgba(4)
       call igTextColored(col,c_loc(str1))
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

  !> Create a menu item with the given label. If keybind is present
  !> use the key bind text as shortcut. If selected, mark the menu
  !> item as selected (default = false). If enabled, mark the menu
  !> item as enabled (default = true).
  module function iw_menuitem(label,keybind,selected,enabled,shortcut_text)
    use interfaces_cimgui
    use keybindings, only: get_bind_keyname
    character(len=*,kind=c_char), intent(in) :: label
    integer, intent(in), optional :: keybind
    logical, intent(in), optional :: selected
    logical, intent(in), optional :: enabled
    character(len=*,kind=c_char), intent(in), optional :: shortcut_text
    logical :: iw_menuitem

    character(len=:,kind=c_char), allocatable, target :: str1, str2
    type(c_ptr) :: shortcutptr
    logical(c_bool) :: selected_, enabled_

    str1 = trim(label) // c_null_char
    if (present(keybind)) then
       str2 = trim(get_bind_keyname(keybind)) // c_null_char
       shortcutptr = c_loc(str2)
    elseif (present(shortcut_text)) then
       str2 = trim(shortcut_text) // c_null_char
       shortcutptr = c_loc(str2)
    else
       shortcutptr = c_null_ptr
    end if
    selected_ = .false._c_bool
    if (present(selected)) selected_ = selected
    enabled_ = .true._c_bool
    if (present(enabled)) enabled_ = enabled

    iw_menuitem = igMenuItem_Bool(c_loc(str1),shortcutptr,selected_,enabled_)

  end function iw_menuitem

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
  module subroutine iw_tooltip(str,ttshown,rgba,nowrap)
    use interfaces_cimgui
    use gui_main, only: tooltip_wrap_factor, tooltip_delay, tooltip_enabled, fontsize
    character(len=*,kind=c_char), intent(in) :: str
    logical, intent(inout), optional :: ttshown
    real(c_float), intent(in), optional :: rgba(4)
    logical, intent(in), optional :: nowrap

    character(len=:,kind=c_char), allocatable, target :: strloc
    integer :: flags
    type(ImVec4) :: col
    logical :: nowrap_

    nowrap_ = .false.
    if (present(nowrap)) nowrap_ = nowrap
    if (present(rgba)) then
       col%x = rgba(1)
       col%y = rgba(2)
       col%z = rgba(3)
       col%w = rgba(4)
       call igPushStyleColor_Vec4(ImGuiCol_Text,col)
    end if

    if (.not.tooltip_enabled) return
    flags = ImGuiHoveredFlags_AllowWhenBlockedByPopup
    if (present(ttshown)) then
       if (igIsItemHovered_delayed(flags,tooltip_delay,ttshown)) call show_tooltip()
    else
       if (igIsItemHovered(flags)) call show_tooltip()
    end if

    if (present(rgba)) &
       call igPopStyleColor(1)

  contains
    subroutine show_tooltip()
      strloc = trim(str) // c_null_char
      call igBeginTooltip()
      if (nowrap_) then
         call igText(c_loc(strloc))
      else
         call igPushTextWrapPos(tooltip_wrap_factor * fontsize%x)
         call igTextWrapped(c_loc(strloc))
         call igPopTextWrapPos()
      end if
      call igEndTooltip()
    end subroutine show_tooltip
  end subroutine iw_tooltip

  !> Create a selectable that highlights the current row. Return true
  !> if the selectable is hovered. If clicked is present, return .true.
  !> if clicked.
  module function iw_highlight_selectable(str,clicked)
    use interfaces_cimgui
    use gui_main, only: g, ColorTableHighlightRow
    character(len=*,kind=c_char), intent(in) :: str
    logical, intent(out), optional :: clicked
    logical :: iw_highlight_selectable

    type(ImVec2) :: sz1, szero
    real(c_float) :: pos
    integer(c_int) :: flags
    character(kind=c_char,len=:), allocatable, target :: str2
    logical :: ldum

    szero%x = 0
    szero%y = 0
    sz1%x = g%Style%ItemSpacing%x
    sz1%y = g%Style%ItemSpacing%y + g%Style%CellPadding%y + g%Style%FramePadding%y
    call igPushStyleVar_Vec2(ImGuiStyleVar_ItemSpacing,sz1)
    call igPushStyleColor_Vec4(ImGuiCol_HeaderHovered,ColorTableHighlightRow)
    pos = igGetCursorPosX()
    flags = ImGuiSelectableFlags_SpanAllColumns
    flags = ior(flags,ImGuiSelectableFlags_AllowItemOverlap)
    flags = ior(flags,ImGuiSelectableFlags_SelectOnNav)
    str2 = trim(str) // c_null_char
    call igSameLine(0._c_float,0._c_float)
    ldum = igSelectable_Bool(c_loc(str2),.false._c_bool,flags,szero)
    if (present(clicked)) clicked = ldum
    call igSetCursorPosX(pos)
    iw_highlight_selectable = .false.
    if (igIsItemHovered(ImGuiHoveredFlags_None)) then
       iw_highlight_selectable = &
          igIsMouseHoveringRect(g%LastItemData%NavRect%min,g%LastItemData%NavRect%max,.false._c_bool)
    end if
    call igPopStyleColor(1_c_int)
    call igPopStyleVar(1_c_int)

  end function iw_highlight_selectable

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

  !> Calculate a nice starting position for the next window based on
  !> the existing windows.
  module subroutine get_nice_next_window_pos(pos)
    use windows, only: nwin, win
    use interfaces_cimgui
    type(ImVec2), intent(out) :: pos

    integer :: i
    real(c_float) :: step
    logical :: found

    step = igGetTextLineHeightWithSpacing()
    pos%x = step
    pos%y = step

    found = .true.
    do while (found)
       pos%x = pos%x + step
       pos%y = pos%y + step

       found = .false.
       do i = 1, nwin
          if (win(i)%isopen.and..not.win(i)%isdocked) then
             if (abs(win(i)%pos(1) - pos%x) < step .and. abs(win(i)%pos(2) - pos%y) < step) then
                found = .true.
                exit
             end if
          end if
       end do
    end do

  end subroutine get_nice_next_window_pos

  !> Get the current working directory from imgui. No trailing /.
  module function get_current_working_dir()
    use interfaces_cimgui, only: getCurrentWorkDir
    use param, only: dirsep
    character(len=:), allocatable :: get_current_working_dir

    character(kind=c_char,len=:), allocatable, target :: strc

    integer :: idum, in

    allocate(character(len=10241) :: strc)
    idum = getCurrentWorkDir(c_loc(strc),10240_c_size_t)
    in = index(strc,c_null_char)-1
    if (strc(in:in) == dirsep.and.in > 0) in = in - 1
    get_current_working_dir = strc(1:in)

  end function get_current_working_dir

end submodule proc
