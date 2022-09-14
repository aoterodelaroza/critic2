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

! Fortran interfaces for opengl3
module gui_interfaces_opengl3
  use iso_c_binding
  implicit none

  public

  !xx! OpenGL constants
  include "interfaces/opengl3_constants.inc"

  interface
     !xx! GL3W
     include "interfaces/gl3w_proc.inc"

     !xx! OpenGL
     include "interfaces/opengl3_proc.inc"
  end interface

end module gui_interfaces_opengl3
