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

  ! deferred-commit state for the iw_inputtext/iw_inputint/iw_inputint3/iw_inputfloat/
  ! iw_inputfloat3 widgets (notlive): the ImGui ID and in-progress value of the input
  ! currently being edited. Only one item is active at a time, so the edit ID is shared
  ! across all of them; only the value buffer differs by type.
  integer(c_int) :: input_editid = 0_c_int
  integer(c_int) :: inputint_editbuf = 0_c_int
  integer(c_int) :: inputint3_editbuf(3) = 0_c_int
  real(c_float) :: inputfloat_editbuf = 0._c_float
  real(c_float) :: inputfloat3_editbuf(3) = 0._c_float
  character(kind=c_char,len=:), allocatable :: inputtext_editbuf

contains

  !> Draw a periodic table and return the selected atomic number.
  module function iw_periodictable()
    use interfaces_cimgui
    use tools_io, only: nameguess, string
    use gui_main, only: g, ColorElement, ColorBlack
    use param, only: elemname
    integer :: iw_periodictable

    integer :: irow, icol, iz
    real(c_float) :: width
    type(ImVec2) :: szero
    real(c_float) :: ssquare(2)
    type(ImVec4) :: color = ColorBlack
    integer, parameter :: ptable(18,9) = reshape((/&
        1, -1, -1, -1, -1, -1, -1,  0, -1, -1, -1, -1, -1, -1, -1, -1, -1,  2,&
        3,  4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  5,  6,  7,  8,  9, 10,&
       11, 12, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 13, 14, 15, 16, 17, 18,&
       19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36,&
       37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54,&
       55, 56, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86,&
       87, 88,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,&
       -1, -1, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, -1, -1,&
       -1, -1, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99,100,101,102, -1, -1 &
       /),(/18,9/))

    ! initialize
    szero%x = 0._c_float
    szero%y = 0._c_float
    width = iw_calcwidth(2,1)
    ssquare = width

    ! set some padding variables to zero
    call igPushStyleVar_Vec2(ImGuiStyleVar_FramePadding,szero)

    ! draw the buttons
    iw_periodictable = -1
    do irow = 1, size(ptable,2)
       do icol = 1, size(ptable,1)
          iz = ptable(icol,irow)
          if (iz >= 0) then
             color%x = ColorElement(1,iz)
             color%y = ColorElement(2,iz)
             color%z = ColorElement(3,iz)
             call igPushStyleColor_Vec4(ImGuiCol_Button,color)
             call igPushStyleColor_Vec4(ImGuiCol_Text,ColorBlack)

             call igSameLine((icol-1)*width + g%Style%WindowPadding%y,1._c_float)
             if (iw_button(nameguess(iz,.true.),siz=ssquare)) then
                iw_periodictable = iz
                call igPopStyleColor(2)
                call igPopStyleVar(1_c_int)
                return
             end if
             call igPopStyleColor(2)
             call iw_tooltip(trim(elemname(iz)) // " (Z=" // string(iz) // ")")
          end if
       end do
       call igNewLine()
    end do

    ! reset the padding
    call igPopStyleVar(1_c_int)

  end function iw_periodictable

  !> Create an input text widget with the given label and the given
  !> width, and buffer size bufsize. The text is given as either texta
  !> (allocatable string) or textf (fixed string). grabfocus = grab focus
  !> when drawn. sameline = place in the same line as the previous widget.
  !> flags = flags from ImGuiInputTextFlags_*.
  module function iw_inputtext(label,bufsize,texta,textf,width,grabfocus,sameline,notlive,flags,&
     nlines)
    use interfaces_cimgui
    character(len=*), intent(in) :: label
    integer, intent(in) :: bufsize
    character(len=:), allocatable, intent(inout), optional :: texta
    character(len=*), intent(inout), optional :: textf
    integer, intent(in), optional :: width
    logical, intent(in), optional :: grabfocus
    logical, intent(in), optional :: sameline
    logical, intent(in), optional :: notlive
    integer(c_int), intent(in), optional :: flags
    integer, intent(in), optional :: nlines
    logical :: iw_inputtext

    integer(c_int) :: flags_, myid
    character(kind=c_char,len=:), allocatable, target :: label_, text_
    character(kind=c_char,len=:), allocatable :: src
    integer :: i, ll
    logical :: sameline_, notlive_, editing
    type(ImVec2) :: szml

    ! process input options
    flags_ = ImGuiInputTextFlags_None
    if (present(flags)) flags_ = flags
    sameline_ = .false.
    if (present(sameline)) sameline_ = sameline
    notlive_ = .false.
    if (present(notlive)) notlive_ = notlive

    ! must have exactly one of the text outputs
    if (.not.present(texta) .and. .not.present(textf)) return
    label_ = trim(label) // c_null_char
    myid = igGetID_Str(c_loc(label_))

    ! choose the source string used to fill the edit buffer this frame. With deferred
    ! commit, keep editing a persistent buffer for the active item so texta/textf are
    ! left untouched while typing and only updated when the value is committed.
    if (notlive_ .and. myid == input_editid) then
       src = inputtext_editbuf
    elseif (present(texta)) then
       src = texta
    else
       src = trim(textf)
    end if

    ! set up the C buffer from the source string
    allocate(character(len=bufsize+1) :: text_)
    ll = min(len(src),bufsize)
    if (ll > 0) text_(1:ll) = src(1:ll)
    text_(ll+1:) = c_null_char

    ! push width
    if (present(width)) &
       call igPushItemWidth(iw_calcwidth(width,1))

    ! same line
    if (sameline_) &
       call igSameLine(0._c_float,-1._c_float)

    ! grab focus
    if (present(grabfocus)) then
       if (grabfocus) call igSetKeyboardFocusHere(0_c_int)
    end if

    ! call inputtext (multiline if a number of lines was given)
    if (present(nlines)) then
       szml%x = 0._c_float
       szml%y = iw_calcheight(nlines,0,.false.)
       iw_inputtext = logical(igInputTextMultiline(c_loc(label_),c_loc(text_),int(bufsize,c_size_t),&
          szml,flags_,c_null_funptr,c_null_ptr))
    else
       iw_inputtext = logical(igInputText(c_loc(label_),c_loc(text_),int(bufsize,c_size_t),flags_,c_null_funptr,c_null_ptr))
    end if

    ! pop the width
    if (present(width)) &
       call igPopItemWidth()

    ! text currently held in the C buffer
    ll = index(text_,c_null_char)-1
    if (ll < 0) ll = len(text_)

    if (notlive_) then
       ! deferred commit: edit a persistent buffer so texta/textf are left untouched while
       ! typing and only updated when the value is committed (Enter, Tab, or click-away). We do
       ! not use ImGuiInputTextFlags_EnterReturnsTrue, which would discard the text on Tab/click.
       call input_edit_update(myid,editing,iw_inputtext)
       if (editing) inputtext_editbuf = text_(1:ll)
    end if

    ! write the text out: on commit if deferred, on every change if live
    if (iw_inputtext) then
       if (present(texta)) then
          texta = text_(1:ll)
       else
          textf = text_(1:ll)
       end if
    end if

  end function iw_inputtext

  !> Drag float button for 1, 2, 3, and 4 floating point numbers. The
  !> number of buttons shown is determined by whether the x1, x2, x3,
  !> or x4 argument is passed.  Speed = step for the dragfloat. Min
  !> and max = minimum and maximum values. scale = scale the numbers
  !> by this value before and after the drag. decimal = the number of
  !> decimal places on the label. flags = combination of
  !> ImGuiSliderFlags_* flags. Version for real(c_float) type.
  module function iw_dragfloat_realc(str,x1,x2,x3,x4,speed,min,max,scale,decimal,&
     sameline,notlive,flags)
    use interfaces_cimgui
    use tools_io, only: string
    character(len=*,kind=c_char), intent(in) :: str
    real(c_float), intent(inout), optional :: x1
    real(c_float), intent(inout), optional :: x2(2)
    real(c_float), intent(inout), optional :: x3(3)
    real(c_float), intent(inout), optional :: x4(4)
    real(c_float), intent(in), optional :: speed, min, max, scale
    integer, intent(in), optional :: decimal
    logical, intent(in), optional :: sameline
    logical, intent(in), optional :: notlive
    integer(c_int), intent(in), optional :: flags
    logical :: iw_dragfloat_realc

    real(c_float) :: speed_, min_, max_, scale_
    character(len=:,kind=c_char), allocatable, target :: str_, sformat_
    integer(c_int) :: flags_
    logical :: notlive_
    real(c_float) :: x1_, x2_(2), x3_(3), x4_(4)
    integer :: n, decimal_
    real(c_float) :: width
    logical :: sameline_

    ! process options
    str_ = trim(str) // c_null_char
    speed_ = 1._c_float
    if (present(speed)) speed_ = speed
    min_ = -FLT_MAX
    if (present(min)) min_ = min
    max_ = FLT_MAX
    if (present(max)) max_ = max
    scale_ = 1._c_float
    if (present(scale)) scale_ = scale
    decimal_ = 3
    if (present(decimal)) decimal_ = decimal
    sformat_ = "%." // string(decimal_) // "f" // c_null_char
    sameline_ = .false.
    if (present(sameline)) sameline_ = sameline
    flags_ = 0_c_int
    if (present(flags)) flags_ = flags
    notlive_ = .false.
    if (present(notlive)) notlive_ = notlive

    ! same line
    if (sameline_) &
       call igSameLine(0._c_float,-1._c_float)

    ! calculate width
    if (present(x1)) then
       n = 1
    elseif (present(x2)) then
       n = 2
    elseif (present(x3)) then
       n = 3
    elseif (present(x4)) then
       n = 4
    end if
    width = dragfloat_width(n,decimal_,min_,max_,present(min),present(max))

    ! draw the float
    call igPushItemWidth(width)
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
    call igPopItemWidth()

    if (notlive_) then
       if (iw_dragfloat_realc .and. igIsItemActive()) &
          iw_dragfloat_realc = igIsMouseDragging(ImGuiMouseButton_Left,-1._c_float)
       ! also commit when the user leaves the field by Enter, Tab, or click-away
       iw_dragfloat_realc = iw_dragfloat_realc .or. logical(igIsItemDeactivatedAfterEdit())
    end if

  end function iw_dragfloat_realc

  !> Drag float button for 1, 2, 3, and 4 floating point numbers. The
  !> number of buttons shown is determined by whether the x1, x2, x3,
  !> or x4 argument is passed.  Speed = step for the dragfloat. Min
  !> and max = minimum and maximum values. scale = scale the numbers
  !> by this value before and after the drag. sformat = C format of the
  !> label. flags = combination of ImGuiSliderFlags_* flags. Version
  !> for real*8 type.
  module function iw_dragfloat_real8(str,x1,x2,x3,x4,speed,min,max,scale,decimal,&
     sameline,notlive,committed,flags)
    use interfaces_cimgui
    use tools_io, only: string
    character(len=*,kind=c_char), intent(in) :: str
    real*8, intent(inout), optional :: x1
    real*8, intent(inout), optional :: x2(2)
    real*8, intent(inout), optional :: x3(3)
    real*8, intent(inout), optional :: x4(4)
    real*8, intent(in), optional :: speed, min, max, scale
    integer, intent(in), optional :: decimal
    logical, intent(in), optional :: sameline
    logical, intent(in), optional :: notlive
    logical, intent(out), optional :: committed
    integer(c_int), intent(in), optional :: flags
    logical :: iw_dragfloat_real8

    real*8 :: scale_
    real(c_float) :: speed_, min_, max_, x1_, x2_(2), x3_(3), x4_(4)
    character(len=:,kind=c_char), allocatable, target :: str_, sformat_
    integer(c_int) :: flags_
    logical :: notlive_
    integer :: n, decimal_
    real(c_float) :: width
    logical :: sameline_

    ! process options
    str_ = trim(str) // c_null_char
    speed_ = 1._c_float
    if (present(speed)) speed_ = real(speed,c_float)
    min_ = -FLT_MAX
    if (present(min)) min_ = real(min,c_float)
    max_ = FLT_MAX
    if (present(max)) max_ = real(max,c_float)
    scale_ = 1._c_float
    if (present(scale)) scale_ = scale
    decimal_ = 3
    if (present(decimal)) decimal_ = decimal
    sformat_ = "%." // string(decimal_) // "f" // c_null_char
    sameline_ = .false.
    if (present(sameline)) sameline_ = sameline
    flags_ = 0_c_int
    if (present(flags)) flags_ = flags
    notlive_ = .false.
    if (present(notlive)) notlive_ = notlive

    ! same line
    if (sameline_) &
       call igSameLine(0._c_float,-1._c_float)

    ! calculate width
    if (present(x1)) then
       n = 1
    elseif (present(x2)) then
       n = 2
    elseif (present(x3)) then
       n = 3
    elseif (present(x4)) then
       n = 4
    end if
    width = dragfloat_width(n,decimal_,min_,max_,present(min),present(max))

    ! draw the float
    call igPushItemWidth(width)
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
    call igPopItemWidth()

    if (present(committed)) committed = .false.
    if (notlive_) then
       if (iw_dragfloat_real8 .and. igIsItemActive()) &
          iw_dragfloat_real8 = igIsMouseDragging(ImGuiMouseButton_Left,-1._c_float)
       ! also commit when the user leaves the field by Enter, Tab, or click-away;
       ! report that final commit separately so callers can distinguish it from a
       ! drag still in progress
       if (present(committed)) committed = logical(igIsItemDeactivatedAfterEdit())
       iw_dragfloat_real8 = iw_dragfloat_real8 .or. logical(igIsItemDeactivatedAfterEdit())
    end if

  end function iw_dragfloat_real8

  !> Integer input box around igInputInt. step and step_fast are the
  !> increments for the +/- buttons (default 0 = hide the buttons).
  !> width = field width in number of characters. sameline = draw in the
  !> same line as the previous widget. flags = ImGuiInputTextFlags_*. If
  !> notlive, return .true. only when the value is committed (Enter,
  !> Tab to the next item, or clicking away), instead of on every keystroke.
  module function iw_inputint(str,ival,step,step_fast,width,sameline,notlive,flags)
    use interfaces_cimgui
    character(len=*,kind=c_char), intent(in) :: str
    integer(c_int), intent(inout) :: ival
    integer(c_int), intent(in), optional :: step, step_fast
    integer, intent(in), optional :: width
    logical, intent(in), optional :: sameline
    logical, intent(in), optional :: notlive
    integer(c_int), intent(in), optional :: flags
    logical :: iw_inputint

    character(len=:,kind=c_char), allocatable, target :: str_
    integer(c_int) :: flags_, step_, step_fast_, v, myid
    logical :: sameline_, notlive_, editing

    ! process options
    str_ = trim(str) // c_null_char
    step_ = 0_c_int
    if (present(step)) step_ = step
    step_fast_ = 0_c_int
    if (present(step_fast)) step_fast_ = step_fast
    sameline_ = .false.
    if (present(sameline)) sameline_ = sameline
    flags_ = 0_c_int
    if (present(flags)) flags_ = flags
    notlive_ = .false.
    if (present(notlive)) notlive_ = notlive

    ! same line
    if (sameline_) &
       call igSameLine(0._c_float,-1._c_float)

    ! width
    if (present(width)) &
       call igPushItemWidth(iw_calcwidth(width,1))

    if (notlive_) then
       ! deferred commit: edit a persistent buffer so ival is left untouched while typing and
       ! only updated when the value is committed (Enter, Tab, or click-away). We do not use
       ! ImGuiInputTextFlags_EnterReturnsTrue, which would discard the typed value on Tab/click.
       myid = igGetID_Str(c_loc(str_))
       if (myid == input_editid) then
          v = inputint_editbuf
       else
          v = ival
       end if
       iw_inputint = logical(igInputInt(c_loc(str_),v,step_,step_fast_,flags_))
       call input_edit_update(myid,editing,iw_inputint)
       if (editing) inputint_editbuf = v
       if (iw_inputint) ival = v
    else
       ! live: ival updates on every change
       iw_inputint = logical(igInputInt(c_loc(str_),ival,step_,step_fast_,flags_))
    end if

    ! pop the width
    if (present(width)) &
       call igPopItemWidth()

  end function iw_inputint

  !> Input box for three integers, wrapping igInputInt3. width = field width
  !> in number of characters. sameline = draw in the same line as the previous
  !> widget. flags = ImGuiInputTextFlags_*. If notlive, ival is left untouched
  !> while typing and updated (and .true. returned) only when the value is
  !> committed (Enter, Tab, or click-away), instead of on every keystroke.
  module function iw_inputint3(str,ival,width,sameline,notlive,flags)
    use interfaces_cimgui
    character(len=*,kind=c_char), intent(in) :: str
    integer(c_int), intent(inout) :: ival(3)
    integer, intent(in), optional :: width
    logical, intent(in), optional :: sameline
    logical, intent(in), optional :: notlive
    integer(c_int), intent(in), optional :: flags
    logical :: iw_inputint3

    character(len=:,kind=c_char), allocatable, target :: str_
    integer(c_int) :: flags_, v(3), myid
    logical :: sameline_, notlive_, editing

    ! process options
    str_ = trim(str) // c_null_char
    sameline_ = .false.
    if (present(sameline)) sameline_ = sameline
    flags_ = 0_c_int
    if (present(flags)) flags_ = flags
    notlive_ = .false.
    if (present(notlive)) notlive_ = notlive

    ! same line
    if (sameline_) &
       call igSameLine(0._c_float,-1._c_float)

    ! width
    if (present(width)) &
       call igPushItemWidth(iw_calcwidth(width,1))

    if (notlive_) then
       ! deferred commit: edit a persistent buffer so ival is left untouched
       ! while typing (see iw_inputint)
       myid = igGetID_Str(c_loc(str_))
       if (myid == input_editid) then
          v = inputint3_editbuf
       else
          v = ival
       end if
       iw_inputint3 = logical(igInputInt3(c_loc(str_),v,flags_))
       call input_edit_update(myid,editing,iw_inputint3)
       if (editing) inputint3_editbuf = v
       if (iw_inputint3) ival = v
    else
       ! live: ival updates on every change
       iw_inputint3 = logical(igInputInt3(c_loc(str_),ival,flags_))
    end if

    ! pop the width
    if (present(width)) &
       call igPopItemWidth()

  end function iw_inputint3

  !> Input box for a floating-point number, wrapping igInputFloat. step and
  !> step_fast are the increments for the +/- buttons (default 0 = hide them).
  !> decimal = number of decimals shown (default 3). width = field width in
  !> number of characters. sameline = draw in the same line as the previous
  !> widget. flags = ImGuiInputTextFlags_*. If notlive, val is left untouched
  !> while typing and updated (and .true. returned) only when the value is
  !> committed (Enter, Tab, or click-away), instead of on every keystroke.
  module function iw_inputfloat(str,val,step,step_fast,decimal,width,sameline,notlive,flags)
    use interfaces_cimgui
    use tools_io, only: string
    character(len=*,kind=c_char), intent(in) :: str
    real(c_float), intent(inout) :: val
    real(c_float), intent(in), optional :: step, step_fast
    integer, intent(in), optional :: decimal
    integer, intent(in), optional :: width
    logical, intent(in), optional :: sameline
    logical, intent(in), optional :: notlive
    integer(c_int), intent(in), optional :: flags
    logical :: iw_inputfloat

    character(len=:,kind=c_char), allocatable, target :: str_, fmt_
    integer(c_int) :: flags_, myid
    real(c_float) :: step_, step_fast_, v
    logical :: sameline_, notlive_, editing

    ! process options
    str_ = trim(str) // c_null_char
    step_ = 0._c_float
    if (present(step)) step_ = step
    step_fast_ = 0._c_float
    if (present(step_fast)) step_fast_ = step_fast
    fmt_ = "%.3f" // c_null_char
    if (present(decimal)) fmt_ = "%." // string(decimal) // "f" // c_null_char
    sameline_ = .false.
    if (present(sameline)) sameline_ = sameline
    flags_ = 0_c_int
    if (present(flags)) flags_ = flags
    notlive_ = .false.
    if (present(notlive)) notlive_ = notlive

    ! same line
    if (sameline_) &
       call igSameLine(0._c_float,-1._c_float)

    ! width
    if (present(width)) &
       call igPushItemWidth(iw_calcwidth(width,1))

    if (notlive_) then
       ! deferred commit: edit a persistent buffer so val is left untouched
       ! while typing (see iw_inputint)
       myid = igGetID_Str(c_loc(str_))
       if (myid == input_editid) then
          v = inputfloat_editbuf
       else
          v = val
       end if
       iw_inputfloat = logical(igInputFloat(c_loc(str_),v,step_,step_fast_,c_loc(fmt_),flags_))
       call input_edit_update(myid,editing,iw_inputfloat)
       if (editing) inputfloat_editbuf = v
       if (iw_inputfloat) val = v
    else
       ! live: val updates on every change
       iw_inputfloat = logical(igInputFloat(c_loc(str_),val,step_,step_fast_,c_loc(fmt_),flags_))
    end if

    ! pop the width
    if (present(width)) &
       call igPopItemWidth()

  end function iw_inputfloat

  !> Input box for three floating-point numbers, wrapping igInputFloat3.
  !> decimal = number of decimals shown (default 3). width = field width in
  !> number of characters. sameline = draw in the same line as the previous
  !> widget. flags = ImGuiInputTextFlags_*. If notlive, val is left untouched
  !> while typing and updated (and .true. returned) only when the value is
  !> committed (Enter, Tab, or click-away), instead of on every keystroke.
  module function iw_inputfloat3(str,val,decimal,width,sameline,notlive,flags)
    use interfaces_cimgui
    use tools_io, only: string
    character(len=*,kind=c_char), intent(in) :: str
    real(c_float), intent(inout) :: val(3)
    integer, intent(in), optional :: decimal
    integer, intent(in), optional :: width
    logical, intent(in), optional :: sameline
    logical, intent(in), optional :: notlive
    integer(c_int), intent(in), optional :: flags
    logical :: iw_inputfloat3

    character(len=:,kind=c_char), allocatable, target :: str_, fmt_
    integer(c_int) :: flags_, myid
    real(c_float) :: v(3)
    logical :: sameline_, notlive_, editing

    ! process options
    str_ = trim(str) // c_null_char
    fmt_ = "%.3f" // c_null_char
    if (present(decimal)) fmt_ = "%." // string(decimal) // "f" // c_null_char
    sameline_ = .false.
    if (present(sameline)) sameline_ = sameline
    flags_ = 0_c_int
    if (present(flags)) flags_ = flags
    notlive_ = .false.
    if (present(notlive)) notlive_ = notlive

    ! same line
    if (sameline_) &
       call igSameLine(0._c_float,-1._c_float)

    ! width
    if (present(width)) &
       call igPushItemWidth(iw_calcwidth(width,1))

    if (notlive_) then
       ! deferred commit (see iw_inputint)
       myid = igGetID_Str(c_loc(str_))
       if (myid == input_editid) then
          v = inputfloat3_editbuf
       else
          v = val
       end if
       iw_inputfloat3 = logical(igInputFloat3(c_loc(str_),v,c_loc(fmt_),flags_))
       call input_edit_update(myid,editing,iw_inputfloat3)
       if (editing) inputfloat3_editbuf = v
       if (iw_inputfloat3) val = v
    else
       ! live: val updates on every change
       iw_inputfloat3 = logical(igInputFloat3(c_loc(str_),val,c_loc(fmt_),flags_))
    end if

    ! pop the width
    if (present(width)) &
       call igPopItemWidth()

  end function iw_inputfloat3

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

  !> Move the cursor to the bottom of the window and, unless centered, to the
  !> right, leaving room for the trailing widget row: ntext characters of text
  !> and nbutton buttons (as in iw_calcwidth). Used for the final button row of
  !> a dialog, which sits in the bottom-right corner regardless of how much
  !> content is above it. If centered, the row is centered horizontally instead.
  module subroutine iw_setpos_bottomright(ntext,nbutton,ncheck,centered)
    use interfaces_cimgui
    use gui_main, only: g
    integer, intent(in) :: ntext
    integer, intent(in) :: nbutton
    integer, intent(in), optional :: ncheck
    logical, intent(in), optional :: centered

    logical :: centered_
    integer :: ncheck_
    type(ImVec2) :: szavail

    centered_ = .false.
    if (present(centered)) centered_ = centered
    ncheck_ = 0
    if (present(ncheck)) ncheck_ = ncheck

    call igGetContentRegionAvail(szavail)

    ! horizontal: right-aligned (leaving room for the scrollbar) or centered
    if (centered_) then
       call igSetCursorPosX(0.5_c_float * (igGetWindowWidth() - iw_calcwidth(ntext,nbutton,ncheck_)))
    else
       call igSetCursorPosX(iw_calcwidth(ntext,nbutton,ncheck_,from_end=.true.) - g%Style%ScrollbarSize)
    end if

    ! vertical: skip the remaining space, if there is room for one more line
    if (szavail%y > igGetTextLineHeightWithSpacing() + g%Style%WindowPadding%y) &
       call igSetCursorPosY(igGetCursorPosY() + szavail%y - igGetTextLineHeightWithSpacing() - &
       g%Style%WindowPadding%y)

  end subroutine iw_setpos_bottomright

  !> Set up one column of the current table, wrapping igTableSetupColumn. The
  !> column is identified either by id, or by icol, a counter that is
  !> incremented and then used as the id (so successive calls number the columns
  !> from the caller's starting value). If sortid and icolsort are given, record
  !> in icolsort the sorting key sortid associated with this column. flags =
  !> ImGuiTableColumnFlags_* (default: None). width = initial width or weight
  !> (default: 0, automatic).
  module subroutine iw_table_column(label,id,icol,sortid,icolsort,flags,width)
    use interfaces_cimgui
    character(len=*,kind=c_char), intent(in) :: label
    integer(c_int), intent(in), optional :: id
    integer, intent(inout), optional :: icol
    integer, intent(in), optional :: sortid
    integer, intent(inout), optional :: icolsort(0:)
    integer(c_int), intent(in), optional :: flags
    real(c_float), intent(in), optional :: width

    character(len=:,kind=c_char), allocatable, target :: str1
    integer(c_int) :: flags_, id_
    real(c_float) :: width_

    ! process options
    flags_ = ImGuiTableColumnFlags_None
    if (present(flags)) flags_ = flags
    width_ = 0._c_float
    if (present(width)) width_ = width

    ! column id: advance the counter if one was given
    id_ = 0_c_int
    if (present(icol)) then
       icol = icol + 1
       id_ = int(icol,c_int)
    elseif (present(id)) then
       id_ = id
    end if

    ! record the sorting key for this column
    if (present(sortid) .and. present(icolsort)) &
       icolsort(id_) = sortid

    ! not trimmed: some column labels carry a significant trailing blank
    str1 = label // c_null_char
    call igTableSetupColumn(c_loc(str1),flags_,width_,id_)

  end subroutine iw_table_column

  !> Begin a menu with the given label, wrapping igBeginMenu. enabled = whether
  !> the menu can be opened (default: true). If it returns .true., the caller
  !> must close the menu with igEndMenu.
  module function iw_beginmenu(label,enabled)
    use interfaces_cimgui
    character(len=*,kind=c_char), intent(in) :: label
    logical, intent(in), optional :: enabled
    logical :: iw_beginmenu

    character(len=:,kind=c_char), allocatable, target :: str1
    logical(c_bool) :: enabled_

    enabled_ = .true._c_bool
    if (present(enabled)) enabled_ = logical(enabled,c_bool)

    str1 = trim(label) // c_null_char
    iw_beginmenu = logical(igBeginMenu(c_loc(str1),enabled_))

  end function iw_beginmenu

  !> Begin a tab item with the given label, wrapping igBeginTabItem. flags =
  !> ImGuiTabItemFlags_* (default: None). If it returns .true., the caller must
  !> close the tab item with igEndTabItem.
  module function iw_begintabitem(label,flags)
    use interfaces_cimgui
    character(len=*,kind=c_char), intent(in) :: label
    integer(c_int), intent(in), optional :: flags
    logical :: iw_begintabitem

    character(len=:,kind=c_char), allocatable, target :: str1
    integer(c_int) :: flags_

    flags_ = ImGuiTabItemFlags_None
    if (present(flags)) flags_ = flags

    str1 = trim(label) // c_null_char
    iw_begintabitem = logical(igBeginTabItem(c_loc(str1),c_null_ptr,flags_))

  end function iw_begintabitem

  !> Whether the keybinding to close a dialog fired this frame: the OK or close
  !> binding while the dialog is focused (pass w%focused() in focused), or the
  !> close-all binding regardless of focus. If okcloses is false, the OK binding
  !> is not a close event; use it in dialogs whose OK action does not always
  !> close the window, and that handle the OK binding themselves.
  module function iw_close_event(focused,okcloses)
    use keybindings, only: is_bind_event, BIND_OK_FOCUSED_DIALOG, BIND_CLOSE_FOCUSED_DIALOG,&
       BIND_CLOSE_ALL_DIALOGS
    logical, intent(in) :: focused
    logical, intent(in), optional :: okcloses
    logical :: iw_close_event

    logical :: okcloses_

    okcloses_ = .true.
    if (present(okcloses)) okcloses_ = okcloses

    iw_close_event = (focused .and. ((okcloses_ .and. is_bind_event(BIND_OK_FOCUSED_DIALOG)) .or.&
       is_bind_event(BIND_CLOSE_FOCUSED_DIALOG))) .or. is_bind_event(BIND_CLOSE_ALL_DIALOGS)

  end function iw_close_event

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

  !> Calculate the width of ntext characters, nbutton buttons, and ncheck
  !> checkboxes (a checkbox is a tick box one frame high plus the gap to its
  !> label, so it is wider than a button with the same label).
  !>
  !> SomeText    |  Button1  |    |  Button2  |  ||--end of window
  !>             ^--^     ^--^    ^--^     ^--^    <-- FramePadding.x
  !>        ^----^           ^----^                <-- ItemSpacing.x
  !>                                          ^--^ <-- WindowPadding.x
  !>
  !> If from_end, return the position to place the cursor at the
  !> calculated distance from the end of the window, or the current
  !> position if it is negative.
  module function iw_calcwidth(ntext,nbutton,ncheck,from_end)
    use interfaces_cimgui
    use gui_main, only: g
    integer, intent(in) :: ntext
    integer, intent(in) :: nbutton
    integer, intent(in), optional :: ncheck
    logical, intent(in), optional :: from_end
    real(c_float) :: iw_calcwidth

    type(ImVec2) :: sz
    character(len=:,kind=c_char), allocatable, target :: strc
    integer :: ncheck_, nitem

    ncheck_ = 0
    if (present(ncheck)) ncheck_ = ncheck

    ! text size
    allocate(character(len=ntext+1,kind=c_char) :: strc)
    strc(1:ntext) = ""
    strc(ntext+1:ntext+1) = c_null_char
    call igCalcTextSize(sz,c_loc(strc),c_null_ptr,.false._c_bool,-1._c_float)

    ! calculate width: the text, plus the frame padding of each button, plus
    ! the tick box and label gap of each checkbox, plus the spacing between
    ! consecutive items
    iw_calcwidth = sz%x
    nitem = nbutton + ncheck_
    if (nbutton > 0) &
       iw_calcwidth = iw_calcwidth + nbutton * (2 * g%Style%FramePadding%x)
    if (ncheck_ > 0) &
       iw_calcwidth = iw_calcwidth + ncheck_ * (igGetTextLineHeight() + &
          2 * g%Style%FramePadding%y + g%Style%ItemInnerSpacing%x)
    if (nitem > 0) &
       iw_calcwidth = iw_calcwidth + (nitem - 1) * g%Style%ItemSpacing%x
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
  module subroutine iw_combo_simple(str,stropt,ival,sameline,sameline_nospace,changed,noarrow,startsatone)
    use interfaces_cimgui
    use types, only: realloc
    character(len=*,kind=c_char), intent(in) :: str
    character(len=*,kind=c_char), intent(in) :: stropt
    integer, intent(inout) :: ival
    logical, intent(in), optional :: sameline
    logical, intent(in), optional :: sameline_nospace
    logical, intent(out), optional :: changed
    logical, intent(in), optional :: noarrow
    logical, intent(in), optional :: startsatone

    type(ImVec2) :: szero
    character(len=:,kind=c_char), allocatable, target :: str1, str2, preview
    character(len=:,kind=c_char), allocatable, target :: stropt1
    logical :: sameline_, sameline_nospace_, noarrow_
    integer :: ll, maxlen, nword, i, iselect, off, ival0
    integer(c_int) :: flags
    integer, allocatable :: idx(:)
    logical(c_bool) :: selected

    ! process input options
    noarrow_ = .false.
    if (present(noarrow)) noarrow_ = noarrow
    ! offset for 1-based indexing (ival=1 selects the first element); work internally with 0-based ival0
    off = 0
    if (present(startsatone)) then
       if (startsatone) off = 1
    end if
    ival0 = ival - off
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
          if (ival0+1 == nword) preview = stropt1(idx(nword)+1:i-1) // c_null_char
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
    iselect = ival0
    flags = ImGuiComboFlags_None
    if (noarrow_) flags = ior(flags,ImGuiComboFlags_NoArrowButton)
    if (igBeginCombo(c_loc(str1),c_loc(preview),flags)) then
       do i = 1, nword
          str2 = stropt1(idx(i)+1:idx(i+1)-1) // c_null_char
          selected = (i == ival0+1)
          if (igSelectable_Bool(c_loc(str2),selected,ImGuiSelectableFlags_None,szero)) &
             iselect = i-1
          if (selected) &
             call igSetItemDefaultFocus()
       end do
       call igEndCombo()
    end if
    if (present(changed)) changed = (ival0 /= iselect)
    ival = iselect + off

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

  !> Draw an integer stepper widget: an optional left label, a "-"
  !> button (decrement by one), an integer input field, and a "+" button
  !> (increment by one). ival is the value, modified in place. The
  !> built-in input steppers are suppressed; the field width is sized to
  !> the number of digits (ndigit, or computed from ival if absent). The
  !> value is clamped to [minval,maxval] (whichever is present). sameline
  !> draws the widget in the same line as the previous one. If notlive
  !> is present and true, the input commits only on Enter. If tooltip is
  !> present, it is shown when hovering over the input field or either
  !> button. Returns .true. if the value changed.
  module function iw_intstepper(str,ival,label,minval,maxval,ndigit,sameline,notlive,flags,tooltip)
    use interfaces_cimgui
    use gui_main, only: g
    character(len=*,kind=c_char), intent(in) :: str
    integer(c_int), intent(inout) :: ival
    character(len=*,kind=c_char), intent(in), optional :: label
    integer(c_int), intent(in), optional :: minval, maxval
    integer, intent(in), optional :: ndigit
    logical, intent(in), optional :: sameline
    logical, intent(in), optional :: notlive
    integer(c_int), intent(in), optional :: flags
    character(len=*,kind=c_char), intent(in), optional :: tooltip
    logical :: iw_intstepper

    logical :: sameline_, ldum
    integer :: nd
    integer(c_int) :: old
    logical, save :: ttshown = .false.

    old = ival
    sameline_ = .false.
    if (present(sameline)) sameline_ = sameline

    ! field width (number of digits)
    if (present(ndigit)) then
       nd = ndigit
    else
       nd = ceiling(log10(max(ival,1_c_int) + 0.1))
    end if

    ! optional left label (tight against the "-" button)
    if (sameline_) call igSameLine(0._c_float,-1._c_float)
    if (present(label)) then
       call iw_text(label,alignframe=.true.)
       call igSameLine(0._c_float,g%Style%FramePadding%x)
    end if

    ! "-" button / input field / "+" button (buttons auto-repeat when held). The input field
    ! reuses iw_inputint so it inherits the notlive/flags handling (and its no-loss commit).
    call igPushButtonRepeat(.true._c_bool)
    if (iw_button("-##" // str)) ival = ival - 1
    if (present(tooltip)) call iw_tooltip(tooltip,ttshown)
    call igSameLine(0._c_float,0.5_c_float*g%Style%FramePadding%x)
    ldum = iw_inputint("##" // str,ival,width=nd,notlive=notlive,flags=flags)
    if (present(tooltip)) call iw_tooltip(tooltip,ttshown)
    call igSameLine(0._c_float,0.5_c_float*g%Style%FramePadding%x)
    if (iw_button("+##" // str)) ival = ival + 1
    if (present(tooltip)) call iw_tooltip(tooltip,ttshown)
    call igPopButtonRepeat()

    ! clamp and report change
    if (present(minval)) ival = max(ival,minval)
    if (present(maxval)) ival = min(ival,maxval)
    iw_intstepper = (ival /= old)

  end function iw_intstepper

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
     noadvance,copy_to_output,centered,alignframe,rgb,rgba,wrap)
    use interfaces_cimgui
    use gui_main, only: g, ColorHighlightText, ColorDangerText
    use tools_io, only: uout
    character(len=*,kind=c_char), intent(in) :: str
    logical, intent(in), optional :: highlight
    logical, intent(in), optional :: danger
    logical, intent(in), optional :: disabled
    logical, intent(in), optional :: sameline
    logical, intent(in), optional :: sameline_nospace
    logical, intent(in), optional :: noadvance
    logical, intent(in), optional :: copy_to_output
    logical, intent(in), optional :: centered
    logical, intent(in), optional :: alignframe
    real(c_float), intent(in), optional :: rgb(3)
    real(c_float), intent(in), optional :: rgba(4)
    logical, intent(in), optional :: wrap

    character(len=:,kind=c_char), allocatable, target :: str1

    logical :: highlight_, danger_, disabled_, sameline_, sameline_nospace_
    logical :: noadvance_,copy_to_output_, centered_, alignframe_, wrap_
    logical :: pushedcolor
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
    alignframe_ = .false.
    wrap_ = .false.
    if (present(wrap)) wrap_ = wrap
    if (present(highlight)) highlight_ = highlight
    if (present(danger)) danger_ = danger
    if (present(sameline)) sameline_ = sameline
    if (present(sameline_nospace)) sameline_nospace_ = sameline_nospace
    if (present(disabled)) disabled_ = disabled
    if (present(noadvance)) noadvance_ = noadvance
    if (present(copy_to_output)) copy_to_output_ = copy_to_output
    if (present(centered)) centered_ = centered
    if (present(alignframe)) alignframe_ = alignframe

    str1 = str // c_null_char
    if (alignframe_) call igAlignTextToFramePadding()
    if (noadvance_) pos = igGetCursorPosX()
    if (sameline_) call igSameLine(0._c_float,-1._c_float)
    if (sameline_nospace_) call igSameLine(0._c_float,0._c_float)
    if (centered_) then
       wwidth = igGetWindowWidth()
       call igCalcTextSize(sz,c_loc(str1),c_null_ptr,.false._c_bool,-1.0_c_float)
       twidth = sz%x
       call igSetCursorPosX((wwidth - twidth) * 0.5_c_float)
    end if
    pushedcolor = .true.
    if (disabled_) then
       col = g%Style%Colors(ImGuiCol_TextDisabled+1)
    elseif (highlight_) then
       col = ColorHighlightText
    elseif (danger_) then
       col = ColorDangerText
    elseif (present(rgb)) then
       col = ImVec4(rgb(1),rgb(2),rgb(3),1._c_float)
    elseif (present(rgba)) then
       col = ImVec4(rgba(1),rgba(2),rgba(3),rgba(4))
    else
       pushedcolor = .false.
    end if
    if (pushedcolor) call igPushStyleColor_Vec4(ImGuiCol_Text,col)
    ! wrap at the window edge if requested
    if (wrap_) call igPushTextWrapPos(0._c_float)
    ! unformatted: the text is data, not a printf format string
    call igTextUnformatted(c_loc(str1),c_null_ptr)
    if (wrap_) call igPopTextWrapPos()
    if (pushedcolor) call igPopStyleColor(1)
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
  module function iw_menuitem(label,keybind,selected,enabled,shortcut_text,danger)
    use interfaces_cimgui
    use keybindings, only: get_bind_keyname
    use gui_main, only: ColorDangerText
    character(len=*,kind=c_char), intent(in) :: label
    integer, intent(in), optional :: keybind
    logical, intent(in), optional :: selected
    logical, intent(in), optional :: enabled
    character(len=*,kind=c_char), intent(in), optional :: shortcut_text
    logical, intent(in), optional :: danger
    logical :: iw_menuitem

    character(len=:,kind=c_char), allocatable, target :: str1, str2
    type(c_ptr) :: shortcutptr
    logical(c_bool) :: selected_, enabled_
    logical :: danger_

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
    danger_ = .false.
    if (present(danger)) danger_ = danger

    if (danger_) call igPushStyleColor_Vec4(ImGuiCol_Text,ColorDangerText)
    iw_menuitem = igMenuItem_Bool(c_loc(str1),shortcutptr,selected_,enabled_)
    if (danger_) call igPopStyleColor(1_c_int)

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

  !> Draw a button labeled str whose background is the color rgb of the
  !> atom it stands for, with the label in black or white, whichever
  !> reads better on it. havergb = the color is known; if false, the
  !> button keeps the default colors. disabled = grayed out and
  !> unresponsive. inert = does not react to the mouse and never
  !> returns pressed, for buttons that only display an atom. Returns
  !> .true. when clicked.
  module function iw_atom_button(str,rgb,havergb,sameline,disabled,inert) result(pressed)
    use interfaces_cimgui
    use gui_main, only: g, ColorButtonHoverFactor, ColorButtonActiveFactor, lumweights,&
       ColorBlack, ColorWhite
    character(len=*,kind=c_char), intent(in) :: str
    real(c_float), intent(in) :: rgb(3)
    logical, intent(in) :: havergb
    logical, intent(in), optional :: sameline
    logical, intent(in), optional :: disabled
    logical, intent(in), optional :: inert
    logical :: pressed

    logical :: inert_
    integer(c_int) :: npop
    type(ImVec4) :: cbase, col4

    inert_ = .false.
    if (present(inert)) inert_ = inert

    ! the background color: the atom's, or the default button color
    if (havergb) then
       cbase = ImVec4(rgb(1),rgb(2),rgb(3),1._c_float)
    else
       cbase = g%Style%Colors(ImGuiCol_Button+1)
    end if

    ! the three button states. An inert button keeps the same color in
    ! all of them, so it does not light up under the mouse. Nothing
    ! needs pushing if the color is the default one and it does react.
    npop = 0
    if (havergb .or. inert_) then
       call igPushStyleColor_Vec4(ImGuiCol_Button,cbase)
       if (inert_) then
          col4 = cbase
       else
          col4 = ImVec4(min(cbase%x*ColorButtonHoverFactor,1._c_float),&
             min(cbase%y*ColorButtonHoverFactor,1._c_float),&
             min(cbase%z*ColorButtonHoverFactor,1._c_float),cbase%w)
       end if
       call igPushStyleColor_Vec4(ImGuiCol_ButtonHovered,col4)
       if (.not.inert_) &
          col4 = ImVec4(cbase%x*ColorButtonActiveFactor,cbase%y*ColorButtonActiveFactor,&
          cbase%z*ColorButtonActiveFactor,cbase%w)
       call igPushStyleColor_Vec4(ImGuiCol_ButtonActive,col4)
       npop = 3
    end if

    ! readable label: black on light atoms, white on dark ones
    if (havergb) then
       if (dot_product(lumweights,rgb) > 0.5_c_float) then
          col4 = ColorBlack
       else
          col4 = ColorWhite
       end if
       call igPushStyleColor_Vec4(ImGuiCol_Text,col4)
       npop = npop + 1
    end if

    pressed = iw_button(str,sameline=sameline,disabled=disabled)
    if (npop > 0) call igPopStyleColor(npop)
    if (inert_) pressed = .false.

  end function iw_atom_button

  !> Draw a clickable icon button of side fontsize%y: an invisible button
  !> (id from strid) with the tinted texture tex drawn on top, brightened
  !> while hovered. If the texture is unavailable (tex==0), the fallback
  !> text glyph is drawn instead, so the control keeps working with no
  !> icon assets. Returns .true. when clicked.
  module function iw_icon_button(strid,tex,rgba,fallback) result(pressed)
    use interfaces_cimgui
    use gui_main, only: fontsize
    character(len=*,kind=c_char), intent(in) :: strid
    integer(c_int), intent(in) :: tex
    real(c_float), intent(in) :: rgba(4)
    character(len=*,kind=c_char), intent(in) :: fallback
    logical :: pressed

    type(ImVec2) :: sz, uv0, uv1
    type(ImVec4) :: tintcol, nobord
    real(c_float) :: posx, col(4)
    character(kind=c_char,len=:), allocatable, target :: strl

    sz = ImVec2(fontsize%y,fontsize%y)
    uv0 = ImVec2(0._c_float,0._c_float)
    uv1 = ImVec2(1._c_float,1._c_float)
    nobord = ImVec4(0._c_float,0._c_float,0._c_float,0._c_float)
    call igAlignTextToFramePadding()
    posx = igGetCursorPosX()

    ! invisible button first (captures the click and hover), then draw the
    ! icon (or a text glyph, if the texture is unavailable) on top
    strl = strid // c_null_char
    pressed = logical(igInvisibleButton(c_loc(strl),sz,ImGuiButtonFlags_None))
    ! brighten the tint toward white while hovered
    col = rgba
    if (igIsItemHovered(ImGuiHoveredFlags_None)) &
       col(1:3) = rgba(1:3) + (1._c_float-rgba(1:3))*0.4_c_float
    call igSameLine(0._c_float,-1._c_float)
    call igSetCursorPosX(posx)
    if (tex /= 0) then
       tintcol = ImVec4(col(1),col(2),col(3),col(4))
       call igImage(int(tex,c_intptr_t),sz,uv0,uv1,tintcol,nobord)
    else
       call iw_text(fallback,rgba=col)
    end if

  end function iw_icon_button

  !> Draw an icon toggle button: the tinted texture tex on a button
  !> whose background is highlighted while state is on and transparent
  !> while off. The icon side is iconscale text line heights (a little
  !> over one, because the artwork is line art and reads as spindly at
  !> exactly the height of a glyph). If the texture is unavailable
  !> (tex==0), a text button with the fallback glyph is drawn instead,
  !> so the control keeps working with no icon assets.
  !> Clicking flips state; returns .true. when that happens. If state
  !> is absent, the button is a momentary flat button (always drawn in
  !> the off style). If popupcontext and popupflags are given, open a
  !> context popup attached to the button and return whether it is
  !> open. If danger, tint the icon red (fallback: danger button
  !> color). scale further multiplies the icon side (default 1), for
  !> palettes where the glyph is the only label.
  module function iw_icon_togglebutton(strid,tex,fallback,state,disabled,sameline,&
     popupcontext,popupflags,danger,scale) result(changed)
    use interfaces_cimgui
    use gui_main, only: g, fontsize, iconscale, ColorDangerButton
    character(len=*,kind=c_char), intent(in) :: strid
    integer(c_int), intent(in) :: tex
    character(len=*,kind=c_char), intent(in) :: fallback
    logical, intent(inout), optional :: state
    logical, intent(in), optional :: disabled
    logical, intent(in), optional :: sameline
    logical, intent(inout), optional :: popupcontext
    integer(c_int), intent(in), optional :: popupflags
    logical, intent(in), optional :: danger
    real(c_float), intent(in), optional :: scale
    logical :: changed

    logical :: disabled_, danger_, state_
    real(c_float) :: side
    type(ImVec2) :: sz, uv0, uv1
    type(ImVec4) :: bgcol, tintcol
    character(kind=c_char,len=:), allocatable, target :: strl

    disabled_ = .false.
    if (present(disabled)) disabled_ = disabled
    danger_ = .false.
    if (present(danger)) danger_ = danger
    state_ = .false.
    if (present(state)) state_ = state
    side = iconscale * fontsize%y
    if (present(scale)) side = scale * iconscale * fontsize%y
    if (present(sameline)) then
       if (sameline) call igSameLine(0._c_float,-1._c_float)
    end if

    ! button background: pressed-button color when on, transparent when off
    if (state_) then
       call igPushStyleColor_Vec4(ImGuiCol_Button,g%Style%Colors(ImGuiCol_ButtonActive+1))
    else
       bgcol = ImVec4(0._c_float,0._c_float,0._c_float,0._c_float)
       call igPushStyleColor_Vec4(ImGuiCol_Button,bgcol)
    end if
    call igBeginDisabled(logical(disabled_,c_bool))
    if (tex /= 0) then
       sz = ImVec2(side,side)
       uv0 = ImVec2(0._c_float,0._c_float)
       uv1 = ImVec2(1._c_float,1._c_float)
       bgcol = ImVec4(0._c_float,0._c_float,0._c_float,0._c_float)
       if (danger_) then
          tintcol = ImVec4(0.86_c_float,0.20_c_float,0.18_c_float,1._c_float)
       else
          tintcol = ImVec4(1._c_float,1._c_float,1._c_float,1._c_float)
       end if
       strl = strid // c_null_char
       changed = logical(igImageButtonEx(igGetID_Str(c_loc(strl)),int(tex,c_intptr_t),sz,uv0,uv1,&
          g%Style%FramePadding,bgcol,tintcol))
    else
       ! same footprint as the icon branch, so a row of these lays out
       ! identically whether or not the assets loaded (a zero height
       ! would give the default frame height, ignoring the scale)
       sz = ImVec2(side + 2._c_float*g%Style%FramePadding%x,&
          side + 2._c_float*g%Style%FramePadding%y)
       strl = fallback // "###" // strid // c_null_char
       if (danger_) &
          call igPushStyleColor_Vec4(ImGuiCol_Button,ColorDangerButton)
       changed = logical(igButton(c_loc(strl),sz))
       if (danger_) &
          call igPopStyleColor(1)
    end if
    call igEndDisabled()
    call igPopStyleColor(1)
    if (changed .and. present(state)) state = .not.state
    ! attach the context popup to the last item (the button in either branch)
    if (present(popupcontext) .and. present(popupflags)) &
       popupcontext = igBeginPopupContextItem(c_null_ptr,popupflags)

  end function iw_icon_togglebutton

  !> Height of the frame that iw_icon_togglebutton draws, for laying
  !> out whatever has to line up with a row of them (a row label, say).
  !> scale is the same optional multiplier that function takes.
  module function iw_iconbutton_height(scale) result(h)
    use interfaces_cimgui
    use gui_main, only: g, fontsize, iconscale
    real(c_float), intent(in), optional :: scale
    real(c_float) :: h

    h = iconscale * fontsize%y
    if (present(scale)) h = scale * h
    h = h + 2._c_float * g%Style%FramePadding%y

  end function iw_iconbutton_height

  !> Draw the standard close button (a red X icon) with the given id.
  !> Returns .true. when clicked.
  module function iw_close_button(strid) result(pressed)
    use icons, only: icon_tex, icon_ui_close
    character(len=*,kind=c_char), intent(in) :: strid
    logical :: pressed

    real(c_float), parameter :: rgba_close(4) = (/0.86_c_float,0.20_c_float,0.18_c_float,1.00_c_float/)

    pressed = iw_icon_button(strid,icon_tex(icon_ui_close),rgba_close,"✕")

  end function iw_close_button

  !> Create a wrapped tooltip, maybe with a delay to show. ttshown
  !> activates the delay, and is the show flag for the delayed tooltip.
  !> whendisabled shows the tooltip also for a disabled widget, which is
  !> where the user most needs to be told why it cannot be used.
  module subroutine iw_tooltip(str,ttshown,rgba,nowrap,whendisabled)
    use interfaces_cimgui
    use gui_main, only: tooltip_wrap_factor, tooltip_delay, tooltip_enabled, fontsize
    character(len=*,kind=c_char), intent(in) :: str
    logical, intent(inout), optional :: ttshown
    real(c_float), intent(in), optional :: rgba(4)
    logical, intent(in), optional :: nowrap
    logical, intent(in), optional :: whendisabled

    character(len=:,kind=c_char), allocatable, target :: strloc
    integer :: flags
    type(ImVec4) :: col
    logical :: nowrap_

    if (.not.tooltip_enabled) return
    nowrap_ = .false.
    if (present(nowrap)) nowrap_ = nowrap
    if (present(rgba)) then
       col%x = rgba(1)
       col%y = rgba(2)
       col%z = rgba(3)
       col%w = rgba(4)
       call igPushStyleColor_Vec4(ImGuiCol_Text,col)
    end if

    flags = ImGuiHoveredFlags_AllowWhenBlockedByPopup
    if (present(whendisabled)) then
       if (whendisabled) flags = ior(flags,ImGuiHoveredFlags_AllowWhenDisabled)
    end if
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
      if (.not.nowrap_) call igPushTextWrapPos(tooltip_wrap_factor * fontsize%x)
      ! unformatted: the text is data, not a printf format string
      call igTextUnformatted(c_loc(strloc),c_null_ptr)
      if (.not.nowrap_) call igPopTextWrapPos()
      call igEndTooltip()
    end subroutine show_tooltip
  end subroutine iw_tooltip

  !> Draw a "(?)" help marker (in the same line as the previous widget
  !> unless sameline is present and false) and attach the tooltip str to
  !> it. The tooltip has no delay (it is shown without a ttshown timer).
  module subroutine iw_helpermark(str,sameline)
    character(len=*,kind=c_char), intent(in) :: str
    logical, intent(in), optional :: sameline

    logical :: sameline_

    sameline_ = .true.
    if (present(sameline)) sameline_ = sameline
    call iw_text("(?)",sameline=sameline_)
    call iw_tooltip(str)

  end subroutine iw_helpermark

  !> A button that opens the arithmetic expressions page of the manual,
  !> for the widgets that take an expression from the user. strid is the
  !> ImGui id suffix ("##...") that tells this button from the others in
  !> the window. Companion of the iw_arith_help text.
  module subroutine iw_arith_help_button(strid,ttshown)
    use interfaces_cimgui
    character(len=*,kind=c_char), intent(in) :: strid
    logical, intent(inout), optional :: ttshown

    character(len=:,kind=c_char), allocatable, target :: str1

    if (iw_button("Help" // strid,sameline=.true.)) then
       str1 = "https://aoterodelaroza.github.io/critic2/manual/arithmetics" // c_null_char
       call openLink(c_loc(str1))
    end if
    call iw_tooltip("Open the manual page regarding arithmetic expressions.",ttshown)

  end subroutine iw_arith_help_button

  !> Create a selectable that highlights the current row. Return true
  !> if the selectable is hovered. If clicked is present, return .true.
  !> if clicked.
  module function iw_highlight_selectable(str,clicked,selected)
    use interfaces_cimgui
    use gui_main, only: g, ColorTableHighlightRow
    character(len=*,kind=c_char), intent(in) :: str
    logical, intent(out), optional :: clicked
    logical, intent(in), optional :: selected
    logical :: iw_highlight_selectable

    type(ImVec2) :: sz1, szero
    real(c_float) :: pos
    integer(c_int) :: flags
    character(kind=c_char,len=:), allocatable, target :: str2
    logical :: ldum
    logical(c_bool) :: sel_

    sel_ = .false._c_bool
    if (present(selected)) sel_ = logical(selected,c_bool)
    szero%x = 0
    szero%y = 0
    sz1%x = g%Style%ItemSpacing%x
    sz1%y = g%Style%ItemSpacing%y + g%Style%CellPadding%y + g%Style%FramePadding%y
    call igPushStyleVar_Vec2(ImGuiStyleVar_ItemSpacing,sz1)
    call igPushStyleColor_Vec4(ImGuiCol_HeaderHovered,ImVec4(ColorTableHighlightRow(1),&
       ColorTableHighlightRow(2),ColorTableHighlightRow(3),ColorTableHighlightRow(4)))
    pos = igGetCursorPosX()
    flags = ImGuiSelectableFlags_SpanAllColumns
    flags = ior(flags,ImGuiSelectableFlags_AllowItemOverlap)
    flags = ior(flags,ImGuiSelectableFlags_SelectOnNav)
    str2 = trim(str) // c_null_char
    call igSameLine(0._c_float,0._c_float)
    ldum = igSelectable_Bool(c_loc(str2),sel_,flags,szero)
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

  !> A duration of t seconds, in whatever unit reads best for its size
  !> (us/ms/s/min/h).
  module function duration_string(t) result(str)
    use tools_io, only: string
    real*8, intent(in) :: t
    character(len=:), allocatable :: str

    if (t < 1d-3) then
       str = string(t*1d6,'f',decimal=0) // " µs"
    elseif (t < 1d0) then
       str = string(t*1d3,'f',decimal=1) // " ms"
    elseif (t < 60d0) then
       str = string(t,'f',decimal=1) // " s"
    elseif (t < 3600d0) then
       str = string(t/60d0,'f',decimal=1) // " min"
    else
       str = string(t/3600d0,'f',decimal=1) // " h"
    end if

  end function duration_string

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

  !> Return file with the extension removed from its base name (the
  !> path, if any, is kept). For espresso/alamode input files
  !> (x.scf.in, x.alm.in), remove the double extension (similar to the
  !> rule in struct_detect_write_format in crystalmod). If the base
  !> name has no extension, return file unchanged (trimmed).
  module function file_name_root(file) result(root)
    use tools_io, only: equal, lower
    use param, only: dirsep
    character(len=*), intent(in) :: file
    character(len=:), allocatable :: root

    character(len=:), allocatable :: ext, ext0
    integer :: idir, idx

    root = trim(file)
    idir = index(root,dirsep,back=.true.)
    idx = index(root(idir+1:),'.',back=.true.)
    if (idx == 0) return
    ext = lower(root(idir+idx+1:))
    root = root(:idir+idx-1)
    ! compound extensions come off whole: .scf.in and .scf.out (quantum
    ! espresso), .alm.in (alamode)
    if (equal(ext,'in') .or. equal(ext,'out')) then
       ext0 = ext
       idx = index(root(idir+1:),'.',back=.true.)
       if (idx > 0) then
          ext = lower(root(idir+idx+1:))
          if (equal(ext,'scf') .or. (equal(ext,'alm') .and. equal(ext0,'in'))) &
             root = root(:idir+idx-1)
       end if
    end if

  end function file_name_root

  !xx! private procedures !xx!

  !> Deferred-commit bookkeeping for the notlive iw_input* widgets; call
  !> immediately after the widget has been drawn. Takes over the edit state if
  !> the widget was just activated, then reports in editing whether this item is
  !> the one being edited (so its persistent buffer must be refreshed with the
  !> value shown this frame) and in commit whether the edit has just been
  !> finished (Enter, Tab, or click-away), in which case the edit state is
  !> released and the caller writes the value back to its own variable. Holding
  !> this in a single place is what enforces the invariant shared by all the
  !> notlive widgets: only one item is being edited at any time, so a single
  !> edit ID is enough and only the value buffer differs by type.
  subroutine input_edit_update(myid,editing,commit)
    use interfaces_cimgui, only: igIsItemActivated, igIsItemDeactivatedAfterEdit
    integer(c_int), intent(in) :: myid
    logical, intent(out) :: editing
    logical, intent(out) :: commit

    if (logical(igIsItemActivated())) input_editid = myid
    editing = (myid == input_editid)
    commit = (editing .and. logical(igIsItemDeactivatedAfterEdit()))
    if (commit) input_editid = 0_c_int

  end subroutine input_edit_update

  !> Width of a drag-float widget holding n numbers with decimal decimals.
  !> Reserve room for the integer part of the field's range (taken from the
  !> bounds vmin/vmax, active only if havemin/havemax) plus the decimals, so a
  !> large maximum (e.g. a temperature up to 10000) is not clipped; small-valued
  !> fields keep the previous default of (4+decimal) characters.
  function dragfloat_width(n,decimal,vmin,vmax,havemin,havemax) result(width)
    use interfaces_cimgui
    use gui_main, only: g
    integer, intent(in) :: n
    integer, intent(in) :: decimal
    real(c_float), intent(in) :: vmin, vmax
    logical, intent(in) :: havemin, havemax
    real(c_float) :: width

    integer :: nintdig, nchar
    real(c_float) :: vabs
    type(ImVec2) :: sz
    character(len=:,kind=c_char), allocatable, target :: strc

    ! width of a single digit
    strc = "0" // c_null_char
    call igCalcTextSize(sz,c_loc(strc),c_null_ptr,.false._c_bool,-1._c_float)

    ! largest absolute value the field can take
    vabs = 1._c_float
    if (havemax .and. abs(vmax) > vabs) vabs = abs(vmax)
    if (havemin .and. abs(vmin) > vabs) vabs = abs(vmin)

    ! integer digits (biased up near powers of ten so 10000 gets 5, never fewer
    ! than 3 to preserve the old (4+decimal) default for small-valued fields)
    nintdig = int(log10(vabs) + 1e-4_c_float) + 1
    if (nintdig < 3) nintdig = 3

    ! characters in the widest value: sign, integer digits, decimal
    ! point (only if there are decimals) and decimals
    nchar = 1 + nintdig + decimal
    if (decimal > 0) nchar = nchar + 1
    width = nchar * sz%x * n + n * (2 * g%Style%FramePadding%x) + (n-1) * g%Style%ItemInnerSpacing%x

  end function dragfloat_width

  !> Draw the periodicity controls: the "Periodicity" label, the
  !> None/Automatic/Manual radio buttons over pertype when it is given, and
  !> the a/b/c cell-count steppers over nc with a Reset button. Without
  !> pertype the steppers are always drawn (a scene, which always has an
  !> explicit cell count); with it they appear only in manual mode.
  !> Returns true if the periodicity changed.
  module function iw_periodicity_widget(nc,ttshown,pertype) result(changed)
    integer(c_int), intent(inout) :: nc(3)
    logical, intent(inout), optional :: ttshown
    integer(c_int), intent(inout), optional :: pertype
    logical :: changed

    logical :: ldum, dosteppers
    integer :: ipad
    integer(c_int) :: nc_(3)

    changed = .false.

    call iw_text("Periodicity",highlight=.true.,alignframe=.true.)

    ! the periodicity type, where the caller has one
    dosteppers = .true.
    if (present(pertype)) then
       changed = changed .or. iw_radiobutton("None",int=pertype,intval=0_c_int,sameline=.true.)
       call iw_tooltip("This object is represented only in the main cell and not repeated by translation",&
          ttshown)
       changed = changed .or. iw_radiobutton("Automatic",int=pertype,intval=1_c_int,sameline=.true.)
       call iw_tooltip("Number of periodic cells controlled by the +/- options in the view window",ttshown)
       changed = changed .or. iw_radiobutton("Manual",int=pertype,intval=2_c_int,sameline=.true.)
       call iw_tooltip("Manually set the number of periodic cells",ttshown)
       dosteppers = (pertype == 2_c_int)
    end if
    if (.not.dosteppers) return

    ! number of cells along each axis (all three share a digit width, so
    ! they do not jump around as the counts grow), then the reset
    ipad = ceiling(log10(max(maxval(nc),1) + 0.1))
    nc_ = nc
    ldum = iw_intstepper("aaxis",nc_(1),label="a:",minval=1_c_int,ndigit=ipad,notlive=.true.,&
       tooltip="Number of cells represented along the a crystallographic axis")
    ldum = iw_intstepper("baxis",nc_(2),label="b:",minval=1_c_int,ndigit=ipad,notlive=.true.,&
       sameline=.true.,tooltip="Number of cells represented along the b crystallographic axis")
    ldum = iw_intstepper("caxis",nc_(3),label="c:",minval=1_c_int,ndigit=ipad,notlive=.true.,&
       sameline=.true.,tooltip="Number of cells represented along the c crystallographic axis")
    if (any(nc_ /= nc)) then
       nc = nc_
       changed = .true.
    end if

    if (iw_button("Reset##periodicity",sameline=.true.)) then
       nc = 1
       changed = .true.
    end if
    call iw_tooltip("Reset the number of cells to one",ttshown)

  end function iw_periodicity_widget

  !> Combo that selects one of the fields of system isys (including the
  !> promolecular field 0), listing only the initialized ones. ifield is
  !> updated and the result is true when the selection changed. width is
  !> the widget width (the available content region if absent, which is
  !> what a table cell wants); nonestr is the preview shown when ifield is
  !> not a field of this system. If onlygrid, list only grid fields.
  module function iw_field_combo(strid,isys,ifield,width,nonestr,onlygrid) result(ch)
    use interfaces_cimgui
    use systems, only: sys
    use fieldmod, only: type_grid
    use tools_io, only: string
    character(len=*,kind=c_char), intent(in) :: strid
    integer, intent(in) :: isys
    integer, intent(inout) :: ifield
    real(c_float), intent(in), optional :: width
    character(len=*), intent(in), optional :: nonestr
    logical, intent(in), optional :: onlygrid
    logical :: ch

    integer :: i
    type(ImVec2) :: szero, szavail
    logical(c_bool) :: issel
    logical :: onlygrid_
    character(kind=c_char,len=:), allocatable, target :: str1, str2

    ch = .false.
    onlygrid_ = .false.
    if (present(onlygrid)) onlygrid_ = onlygrid
    szero%x = 0
    szero%y = 0
    if (sys(isys)%goodfield(ifield)) then
       str2 = string(ifield) // ": " // trim(sys(isys)%f(ifield)%name) // c_null_char
    elseif (present(nonestr)) then
       str2 = nonestr // c_null_char
    else
       str2 = "<none>" // c_null_char
    end if
    str1 = strid // c_null_char
    if (present(width)) then
       call igSetNextItemWidth(width)
    else
       call igGetContentRegionAvail(szavail)
       call igSetNextItemWidth(szavail%x)
    end if
    if (igBeginCombo(c_loc(str1),c_loc(str2),ImGuiComboFlags_None)) then
       do i = 0, sys(isys)%nf
          if (.not.sys(isys)%f(i)%isinit) cycle
          if (onlygrid_ .and. .not.sys(isys)%goodfield(i,type=type_grid)) cycle
          issel = (ifield == i)
          str2 = string(i) // ": " // trim(sys(isys)%f(i)%name) // c_null_char
          if (igSelectable_Bool(c_loc(str2),issel,ImGuiSelectableFlags_None,szero)) then
             if (ifield /= i) then
                ifield = i
                ch = .true.
             end if
          end if
          if (issel) &
             call igSetItemDefaultFocus()
       end do
       call igEndCombo()
    end if

  end function iw_field_combo

  !> Fill lut (3,n) with n colors evenly sampled from colormap icmap (an
  !> index into iw_cmap_name), from the low end of the map to the high end.
  !> Sampling into a table once is much cheaper than asking implot for every
  !> value to be colored, and the table is what the surface and the legend
  !> share.
  module subroutine iw_colormap_lut(icmap,lut)
    use interfaces_cimgui
    integer, intent(in) :: icmap
    real(c_float), intent(out) :: lut(:,:)

    integer :: i, n
    real(c_float) :: t
    type(ImVec4) :: col

    n = size(lut,2)
    do i = 1, n
       t = 0._c_float
       if (n > 1) t = real(i-1,c_float) / real(n-1,c_float)
       call ipSampleColormap(t,iw_colormap_id(icmap),col)
       lut(1,i) = col%x
       lut(2,i) = col%y
       lut(3,i) = col%z
    end do

  end subroutine iw_colormap_lut

  !> The implot id of colormap icmap (an index into iw_cmap_name).
  module function iw_colormap_id(icmap) result(id)
    use interfaces_cimgui
    integer, intent(in) :: icmap
    integer(c_int) :: id

    ! the implot constants are bind(c) variables, so this cannot be a
    ! parameter table; keep the order the same as iw_cmap_name
    select case (icmap)
    case (2)
       id = ImPlotColormap_Plasma
    case (3)
       id = ImPlotColormap_Hot
    case (4)
       id = ImPlotColormap_Cool
    case (5)
       id = ImPlotColormap_Jet
    case (6)
       id = ImPlotColormap_RdBu
    case (7)
       id = ImPlotColormap_Spectral
    case (8)
       id = ImPlotColormap_Greys
    case default
       id = ImPlotColormap_Viridis
    end select

  end function iw_colormap_id

end submodule proc
