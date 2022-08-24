! Copyright (c) 2015-2022 Alberto Otero de la Roza <aoterodelaroza@gmail.com>,
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

! Module for variables configured by the build system & associated tools.
module config
  use param, only: dirsep
  implicit none

  public

  ! enum for the strings from the config flag
  integer, parameter :: istring_package = 1
  integer, parameter :: istring_version = 2
  integer, parameter :: istring_atarget = 3
  integer, parameter :: istring_adate = 4
  integer, parameter :: istring_datadir = 5

  interface
     module function getstring(id) result(str)
       integer, intent(in) :: id
       character(len=:), allocatable :: str
     end function getstring
  end interface

end module config
