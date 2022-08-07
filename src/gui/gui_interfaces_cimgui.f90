! Copyright (c) 2019 Alberto Otero de la Roza <aoterodelaroza@gmail.com>,
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

! Fortran interfaces for cimgui
module gui_interfaces_cimgui
  use iso_c_binding
  implicit none

  public

  !xx! cimgui constants
  include "interfaces/cimgui_constants.inc"

  !xx! cimgui user-defined types
  include "interfaces/cimgui_types.inc"

  interface
     !xx! ImGui GLFW platform backend (imgui_impl_glfw)
     include "interfaces/imgui_glfw_proc.inc"

     !xx! ImGui OpenGL3 renderer backend (imgui_impl_opengl3)
     include "interfaces/imgui_opengl3_proc.inc"

     !xx! cimgui
     include "interfaces/cimgui_proc.inc"
  end interface

end module gui_interfaces_cimgui
