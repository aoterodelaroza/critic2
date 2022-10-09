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

! Some utilities for building the GUI: wrappers around ImGui routines,
! math utilities adapted from the glm library,...
module utils
  use iso_c_binding
  implicit none

  private

  !xx! proc submodule !xx!
  public :: iw_setposx_fromend
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
  !xx! math submodule !xx!
  public :: infiniteperspective
  public :: ortho
  public :: lookat
  public :: project
  public :: unproject

  ! module procedure interfaces
  interface
     !xx! proc submodule !xx!
     module subroutine iw_setposx_fromend(ntext,nbutton)
       integer, intent(in) :: ntext
       integer, intent(in) :: nbutton
     end subroutine iw_setposx_fromend
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
     module subroutine iw_text(str,highlight,disabled,sameline,sameline_nospace,noadvance,&
        copy_to_output)
       character(len=*,kind=c_char), intent(in) :: str
       logical, intent(in), optional :: highlight
       logical, intent(in), optional :: disabled
       logical, intent(in), optional :: sameline
       logical, intent(in), optional :: sameline_nospace
       logical, intent(in), optional :: noadvance
       logical, intent(in), optional :: copy_to_output
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
     !xx! math submodule !xx!
     module function infiniteperspective(fovy,aspect,znear)
       use iso_c_binding, only: c_float
       real(c_float), intent(in) :: fovy, aspect, znear
       real(c_float) :: infiniteperspective(4,4)
     end function infiniteperspective
     module function ortho(left,right,bottom,top,znear,zfar)
       use iso_c_binding, only: c_float
       real(c_float), intent(in) :: left, right, bottom, top, znear, zfar
       real(c_float) :: ortho(4,4)
     end function ortho
     module function lookat(eye,center,up)
       use iso_c_binding, only: c_float
       real(c_float) :: eye(3)
       real(c_float) :: center(3)
       real(c_float) :: up(3)
       real(c_float) :: lookat(4,4)
     end function lookat
     module function project(pos,mview,proj,viewport_a)
       use iso_c_binding, only: c_float, c_int
       real(c_float), intent(in) :: pos(3)
       real(c_float), intent(in) :: mview(4,4)
       real(c_float), intent(in) :: proj(4,4)
       integer(c_int), intent(in) :: viewport_a
       real(c_float) :: project(3)
     end function project
     module function unproject(pos,mview,proj,viewport_a)
       use iso_c_binding, only: c_float, c_int
       real(c_float), intent(in) :: pos(3)
       real(c_float), intent(in) :: mview(4,4)
       real(c_float), intent(in) :: proj(4,4)
       integer(c_int), intent(in) :: viewport_a
       real(c_float) :: unproject(3)
     end function unproject
  end interface

end module utils

