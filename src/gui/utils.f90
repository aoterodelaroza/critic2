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
  public :: iw_clamp_color3
  public :: iw_clamp_color4
  public :: iw_coloredit
  public :: iw_setposx_fromend
  public :: iw_calcheight
  public :: iw_calcwidth
  public :: iw_combo_simple
  public :: iw_radiobutton
  public :: iw_checkbox
  public :: iw_text
  public :: iw_button
  public :: iw_menuitem
  public :: iw_tooltip
  public :: iw_highlight_selectable
  public :: igIsItemHovered_delayed
  public :: get_time_string
  public :: buffer_to_string_array
  public :: get_nice_next_window_pos
  public :: get_current_working_dir
  !xx! math submodule !xx!
  public :: infiniteperspective
  public :: ortho
  public :: lookat
  public :: project
  public :: unproject
  public :: translate
  public :: rotate
  public :: mult
  public :: invmult

  ! module procedure interfaces
  interface
     !xx! proc submodule !xx!
     module subroutine iw_clamp_color3(rgb)
       real(c_float), intent(inout) :: rgb(3)
     end subroutine iw_clamp_color3
     module subroutine iw_clamp_color4(rgba)
       real(c_float), intent(inout) :: rgba(4)
     end subroutine iw_clamp_color4
     module function iw_coloredit(str,rgb,rgba,sameline,nolabel,nointeraction)
       character(len=*,kind=c_char), intent(in) :: str
       real(c_float), intent(inout), optional :: rgb(3)
       real(c_float), intent(inout), optional :: rgba(4)
       logical, intent(in), optional :: sameline, nolabel, nointeraction
       logical :: iw_coloredit
     end function iw_coloredit
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
     module subroutine iw_combo_simple(str,stropt,ival,sameline,sameline_nospace,changed,&
        noarrow)
       character(len=*,kind=c_char), intent(in) :: str
       character(len=*,kind=c_char), intent(in) :: stropt
       integer, intent(inout) :: ival
       logical, intent(in), optional :: sameline
       logical, intent(in), optional :: sameline_nospace
       logical(c_bool), intent(out), optional :: changed
       logical, intent(in), optional :: noarrow
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
     module function iw_checkbox(str,bool,sameline,highlight)
       character(len=*,kind=c_char), intent(in) :: str
       logical, intent(inout) :: bool
       logical, intent(in), optional :: sameline
       logical, intent(in), optional :: highlight
       logical :: iw_checkbox
     end function iw_checkbox
     module subroutine iw_text(str,highlight,danger,disabled,sameline,sameline_nospace,&
        noadvance,copy_to_output,centered,rgb,rgba)
       character(len=*,kind=c_char), intent(in) :: str
       logical, intent(in), optional :: highlight
       logical, intent(in), optional :: danger
       logical, intent(in), optional :: disabled
       logical, intent(in), optional :: sameline
       logical, intent(in), optional :: sameline_nospace
       logical, intent(in), optional :: noadvance
       logical, intent(in), optional :: copy_to_output
       logical, intent(in), optional :: centered
       real(c_float), intent(in), optional :: rgb(3)
       real(c_float), intent(in), optional :: rgba(4)
     end subroutine iw_text
     module function iw_menuitem(label,keybind,selected,enabled,shortcut_text)
       character(len=*,kind=c_char), intent(in) :: label
       integer, intent(in), optional :: keybind
       logical, intent(in), optional :: selected
       logical, intent(in), optional :: enabled
       character(len=*,kind=c_char), intent(in), optional :: shortcut_text
       logical :: iw_menuitem
     end function iw_menuitem
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
     module subroutine iw_tooltip(str,ttshown,rgba,nowrap)
       character(len=*,kind=c_char), intent(in) :: str
       logical, intent(inout), optional :: ttshown
       real(c_float), intent(in), optional :: rgba(4)
       logical, intent(in), optional :: nowrap
     end subroutine iw_tooltip
     module function iw_highlight_selectable(str,selected,clicked)
       character(len=*,kind=c_char), intent(in) :: str
       logical, intent(in), optional :: selected
       logical, intent(out), optional :: clicked
       logical :: iw_highlight_selectable
     end function iw_highlight_selectable
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
     module subroutine get_nice_next_window_pos(pos)
       use interfaces_cimgui, only: ImVec2
       type(ImVec2), intent(out) :: pos
     end subroutine get_nice_next_window_pos
     module function get_current_working_dir()
       character(len=:), allocatable :: get_current_working_dir
     end function get_current_working_dir
     !xx! math submodule !xx!
     module subroutine infiniteperspective(m,fovy,aspect,znear)
       use iso_c_binding, only: c_float
       real(c_float), intent(out) :: m(4,4)
       real(c_float), intent(in) :: fovy, aspect, znear
     end subroutine infiniteperspective
     module subroutine ortho(m,left,right,bottom,top,znear,zfar)
       use iso_c_binding, only: c_float
       real(c_float), intent(out) :: m(4,4)
       real(c_float), intent(in) :: left, right, bottom, top, znear, zfar
     end subroutine ortho
     module subroutine lookat(m,eye,center,up)
       use iso_c_binding, only: c_float
       real(c_float), intent(out) :: m(4,4)
       real(c_float), intent(in) :: eye(3)
       real(c_float), intent(in) :: center(3)
       real(c_float), intent(in) :: up(3)
     end subroutine lookat
     module subroutine project(pos,mview,proj,viewport_a)
       use iso_c_binding, only: c_float, c_int
       real(c_float), intent(inout) :: pos(3)
       real(c_float), intent(in) :: mview(4,4)
       real(c_float), intent(in) :: proj(4,4)
       integer(c_int), intent(in) :: viewport_a
     end subroutine project
     module subroutine unproject(pos,mview,proj,viewport_a)
       use iso_c_binding, only: c_float, c_int
       real(c_float), intent(inout) :: pos(3)
       real(c_float), intent(in) :: mview(4,4)
       real(c_float), intent(in) :: proj(4,4)
       integer(c_int), intent(in) :: viewport_a
     end subroutine unproject
     module subroutine translate(mat,v)
       real(c_float), intent(inout) :: mat(4,4)
       real(c_float), intent(in) :: v(3)
     end subroutine translate
     module subroutine rotate(mat,angle,axis)
       real(c_float), intent(inout) :: mat(4,4)
       real(c_float), intent(in) :: angle
       real(c_float), intent(in) :: axis(3)
     end subroutine rotate
     module subroutine mult(m,mat,v,notrans)
       real(c_float), intent(out) :: m(3)
       real(c_float), intent(in) :: mat(4,4)
       real(c_float), intent(in) :: v(3)
       logical, intent(in), optional :: notrans
     end subroutine mult
     module subroutine invmult(v,mat,notrans)
       real(c_float), intent(inout) :: v(3)
       real(c_float), intent(in) :: mat(4,4)
       logical, intent(in), optional :: notrans
     end subroutine invmult
  end interface

end module utils

