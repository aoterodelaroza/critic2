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

! Icon textures for the GUI
module icons
  use iso_c_binding, only: c_int, c_float
  use param, only: isformat_r_max
  implicit none

  private

  ! icon IDs: system properties
  integer, parameter, public :: icon_prop_fields = 1 ! fields loaded (nested rings)
  integer, parameter, public :: icon_prop_vib = 2    ! vibration data (zigzag)
  integer, parameter, public :: icon_prop_occ = 3    ! partial occupancies (half-filled circle)
  ! icon IDs: UI controls
  integer, parameter, public :: icon_ui_close = 4    ! close button (X)
  integer, parameter, public :: icon_ui_expand = 5   ! expand button (right triangle)
  integer, parameter, public :: icon_ui_collapse = 6 ! collapse button (down triangle)
  integer, parameter, public :: icon_NUM = 6

  ! range of the format icons, indexed by the isformat_r_* constants
  ! in param.F90
  integer, parameter, public :: icon_fmt_MAX = isformat_r_max

  ! OpenGL texture IDs for the icons (0 = unavailable, skip drawing)
  integer(c_int), target, public :: icon_tex(icon_NUM) = 0        ! property icons
  integer(c_int), target, public :: icon_tex_fmt(0:icon_fmt_MAX) = 0 ! format icons

  ! tint colors for the format icons (second index = isformat_r_*);
  ! filled in icons_init
  real(c_float), public :: rgba_icon_fmt(4,0:icon_fmt_MAX)

  public :: icons_init
  public :: icons_end
  public :: format_name

  interface
     module subroutine icons_init()
     end subroutine icons_init
     module subroutine icons_end()
     end subroutine icons_end
     module function format_name(isformat) result(name)
       integer, intent(in) :: isformat
       character(len=:), allocatable :: name
     end function format_name
  end interface

end module icons
