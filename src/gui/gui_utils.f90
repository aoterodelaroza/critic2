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
module gui_utils
  use iso_c_binding
  implicit none

  private

  public :: iw_calcheight
  public :: iw_calcwidth
  public :: iw_combo_simple
  public :: iw_radiobutton
  public :: iw_text
  public :: iw_button
  public :: iw_tooltip
  public :: igIsItemHovered_delayed
  public :: get_time_string
  public :: buffer_to_string_array

  ! module procedure interfaces
  interface
     module function iw_calcheight(npadline,nline,endpad)
       integer, intent(in) :: npadline
       integer, intent(in) :: nline
       logical, intent(in), optional :: endpad
       real(c_float) :: iw_calcheight
     end function iw_calcheight
     module function iw_calcwidth(ntext,nbutton,from_end)
       integer, intent(in) :: ntext
       integer, intent(in) :: nbutton
       logical, intent(in), optional :: from_end
       real(c_float) :: iw_calcwidth
     end function iw_calcwidth
     module subroutine iw_combo_simple(str,stropt,ival,sameline)
       character(len=*,kind=c_char), intent(in) :: str
       character(len=*,kind=c_char), intent(in) :: stropt
       integer, intent(inout) :: ival
       logical, intent(in), optional :: sameline
     end subroutine iw_combo_simple
     module function iw_radiobutton(str,bool,boolval,int,intval,sameline)
       character(len=*,kind=c_char), intent(in) :: str
       logical, intent(inout), optional :: bool
       logical, intent(in), optional :: boolval
       integer(c_int), intent(inout), optional :: int
       integer(c_int), intent(in), optional :: intval
       logical, intent(in), optional :: sameline
       logical :: iw_radiobutton
     end function iw_radiobutton
     module subroutine iw_text(str,highlight,disabled,sameline)
       character(len=*,kind=c_char), intent(in) :: str
       logical, intent(in), optional :: highlight
       logical, intent(in), optional :: sameline
       logical, intent(in), optional :: disabled
     end subroutine iw_text
     module function iw_button(str,danger,sameline,disabled,siz,&
        popupcontext,popupflags)
       character(len=*,kind=c_char), intent(in) :: str
       logical, intent(in), optional :: danger
       logical, intent(in), optional :: sameline
       logical, intent(in), optional :: disabled
       real(c_float), intent(in), optional :: siz(2)
       logical, intent(inout), optional :: popupcontext
       integer(c_int), intent(in), optional :: popupflags
       logical :: iw_button
     end function iw_button
     module subroutine iw_tooltip(str,ttshown)
       character(len=*,kind=c_char), intent(in) :: str
       logical, intent(inout), optional :: ttshown
     end subroutine iw_tooltip
     module function igIsItemHovered_delayed(flags,thr,already_shown)
       integer(c_int), value :: flags
       real(c_float), intent(in) :: thr
       logical, intent(in) :: already_shown
       logical(c_bool) :: igIsItemHovered_delayed
     end function igIsItemHovered_delayed
     module function get_time_string()
       character(len=:), allocatable :: get_time_string
     end function get_time_string
     module subroutine buffer_to_string_array(buf,lu,prefix,suffix)
       character*(*), intent(in) :: buf
       integer, intent(in) :: lu
       character*(*), intent(in), optional :: prefix
       character*(*), intent(in), optional :: suffix
     end subroutine buffer_to_string_array
  end interface

end module gui_utils

