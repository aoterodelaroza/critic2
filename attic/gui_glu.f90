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

!> Fortran bindings for the OpenGL utility library (GLU) API
module gui_glu
  use gui_gl, only: GLdouble, GLint
  use iso_c_binding, only: c_ptr
  implicit none

  public

  ! These constants from the mesa implementation of the OpenGL API on
  ! my computer (version 17.1.4-1). They may need to be changed or
  ! extracted at runtime if they change for other platforms.

  interface
     ! GLUquadric* gluNewQuadric (void);
     function gluNewQuadric() bind(c,name="gluNewQuadric")
       import :: c_ptr
       type(c_ptr) :: gluNewQuadric
     end function gluNewQuadric
     ! void gluSphere(GLUquadric* quad,GLdouble radius,GLint slices,GLint stacks)
     subroutine gluSphere(quad,radius,slices,stacks) bind(c,name="gluSphere")
       import :: GLdouble, GLint, c_ptr
       type(c_ptr), value :: quad
       real(GLdouble), value :: radius
       integer(GLint), value :: slices, stacks
     end subroutine gluSphere
  end interface

contains

end module gui_glu
