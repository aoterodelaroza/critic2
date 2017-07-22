! Copyright (c) 2017 Alberto Otero de la Roza
! <aoterodelaroza@gmail.com>, Robin Myhr <x@example.com>, Isaac
! Visintainer <x@example.com>, Richard Greaves <x@example.com>, Ángel
! Martín Pendás <angel@fluor.quimica.uniovi.es> and Víctor Luaña
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

!> Fortran bindings for the GLFW library API
module gui_glfw
  use iso_c_binding, only: c_ptr, c_int
  implicit none

  public

  ! These constants from the GLFW3 API on my computer (version
  ! 3.2.1-1). They may need to be changed or extracted at runtime if
  ! they change for other platforms.

  interface
     ! void glfwGetFramebufferSize(GLFWwindow* window,int* width,int* height);
     subroutine glfwGetFramebufferSize(window,width,height) bind(c,name="glfwGetFramebufferSize")
       import :: c_ptr, c_int
       type(c_ptr), value :: window
       integer(c_int) :: width, height
     end subroutine glfwGetFramebufferSize
  end interface
  
contains

end module gui_glfw
