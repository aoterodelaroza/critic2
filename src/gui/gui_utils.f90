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
  public :: iw_text
  public :: iw_button
  public :: iw_tooltip
  public :: igIsItemHovered_delayed
  public :: get_time_string

  ! module procedure interfaces
  interface
     module function iw_calcheight(nline)
       use gui_interfaces_cimgui
       integer, intent(in) :: nline
       real(c_float) :: iw_calcheight
     end function iw_calcheight
     module function iw_calcwidth(ntext,nbutton,from_end)
       integer, intent(in) :: ntext
       integer, intent(in) :: nbutton
       logical, intent(in), optional :: from_end
       real(c_float) :: iw_calcwidth
     end function iw_calcwidth
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
  end interface

end module gui_utils

