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

!> Fortran bindings for the OpenGL API
module gui_gl
  use iso_c_binding, only: c_int, c_float, c_double, c_signed_char, &
     c_short
  implicit none

  public

  ! These constants from the mesa implementation of the OpenGL API on
  ! my computer (version 17.1.4-1). They may need to be changed or
  ! extracted at runtime if they change for other platforms.

  ! Datatypes
  integer, parameter :: GLenum     = c_int         ! unsigned int
  integer, parameter :: GLboolean  = c_signed_char ! unsigned char
  integer, parameter :: GLbitfield = c_int         ! unsigned int
  ! integer, parameter :: GLvoid     = void        ! void
  integer, parameter :: GLbyte     = c_signed_char ! signed char
  integer, parameter :: GLshort    = c_short       ! short
  integer, parameter :: GLint      = c_int         ! int
  integer, parameter :: GLubyte    = c_signed_char ! unsigned char
  integer, parameter :: GLushort   = c_short       ! unsigned short
  integer, parameter :: GLuint     = c_int         ! unsigned int
  integer, parameter :: GLsizei    = c_int         ! int
  integer, parameter :: GLfloat    = c_float       ! float
  integer, parameter :: GLclampf   = c_float       ! float
  integer, parameter :: GLdouble   = c_double      ! double
  integer, parameter :: GLclampd   = c_double      ! double

  ! glPush/PopAttrib bits
  integer(GLbitfield), parameter :: GL_CURRENT_BIT         = z'00000001'
  integer(GLbitfield), parameter :: GL_POINT_BIT           = z'00000002'
  integer(GLbitfield), parameter :: GL_LINE_BIT            = z'00000004'
  integer(GLbitfield), parameter :: GL_POLYGON_BIT         = z'00000008'
  integer(GLbitfield), parameter :: GL_POLYGON_STIPPLE_BIT = z'00000010'
  integer(GLbitfield), parameter :: GL_PIXEL_MODE_BIT      = z'00000020'
  integer(GLbitfield), parameter :: GL_LIGHTING_BIT        = z'00000040'
  integer(GLbitfield), parameter :: GL_FOG_BIT             = z'00000080'
  integer(GLbitfield), parameter :: GL_DEPTH_BUFFER_BIT    = z'00000100'
  integer(GLbitfield), parameter :: GL_ACCUM_BUFFER_BIT    = z'00000200'
  integer(GLbitfield), parameter :: GL_STENCIL_BUFFER_BIT  = z'00000400'
  integer(GLbitfield), parameter :: GL_VIEWPORT_BIT        = z'00000800'
  integer(GLbitfield), parameter :: GL_TRANSFORM_BIT       = z'00001000'
  integer(GLbitfield), parameter :: GL_ENABLE_BIT          = z'00002000'
  integer(GLbitfield), parameter :: GL_COLOR_BUFFER_BIT    = z'00004000'
  integer(GLbitfield), parameter :: GL_HINT_BIT            = z'00008000'
  integer(GLbitfield), parameter :: GL_EVAL_BIT            = z'00010000'
  integer(GLbitfield), parameter :: GL_LIST_BIT            = z'00020000'
  integer(GLbitfield), parameter :: GL_TEXTURE_BIT         = z'00040000'
  integer(GLbitfield), parameter :: GL_SCISSOR_BIT         = z'00080000'
  integer(GLbitfield), parameter :: GL_ALL_ATTRIB_BITS     = z'000FFFFF' ! original: 0xffffffff - this is the OR() of all the above

  ! Matrix Mode
  integer(GLenum), parameter :: GL_MATRIX_MODE = z'0BA0'
  integer(GLenum), parameter :: GL_MODELVIEW   = z'1700'
  integer(GLenum), parameter :: GL_PROJECTION  = z'1701'
  integer(GLenum), parameter :: GL_TEXTURE     = z'1702'

  ! Primitives
  integer(GLenum), parameter :: GL_POINTS         = z'0000'
  integer(GLenum), parameter :: GL_LINES          = z'0001'
  integer(GLenum), parameter :: GL_LINE_LOOP      = z'0002'
  integer(GLenum), parameter :: GL_LINE_STRIP     = z'0003'
  integer(GLenum), parameter :: GL_TRIANGLES      = z'0004'
  integer(GLenum), parameter :: GL_TRIANGLE_STRIP = z'0005'
  integer(GLenum), parameter :: GL_TRIANGLE_FAN   = z'0006'
  integer(GLenum), parameter :: GL_QUADS          = z'0007'
  integer(GLenum), parameter :: GL_QUAD_STRIP     = z'0008'
  integer(GLenum), parameter :: GL_POLYGON        = z'0009'

  !! GL functions in alphabetical order !!
  interface
     ! void glBegin(GLenum mode);
     subroutine glBegin(mode) bind(c,name="glBegin")
       import :: GLenum
       integer(GLenum), value :: mode
     end subroutine glBegin
     ! void glClear(GLbitfield mask);
     subroutine glClear(mask) bind(c,name="glClear")
       import :: GLbitfield
       integer(GLbitfield), value :: mask
     end subroutine glClear
     ! void glColor3f(GLfloat r,GLfloat g,GLfloat b);
     subroutine glColor3f(r,g,b) bind(c,name="glColor3f")
       import :: GLfloat
       real(GLfloat), value :: r,g,b
     end subroutine glColor3f
     ! void glClearColor(GLclampf r,GLclampf g,GLclampf b,GLclampf a);
     subroutine glClearColor(r,g,b,a) bind(c,name="glClearColor")
       import :: GLclampf
       real(GLclampf), value :: r,g,b,a
     end subroutine glClearColor
     ! void glEnd(void);
     subroutine glEnd() bind(c,name="glEnd")
     end subroutine glEnd
     ! void glFlush(void);
     subroutine glFlush() bind(c,name="glFlush")
     end subroutine glFlush
     ! void glLoadIdentity(void);
     subroutine glLoadIdentity() bind(c,name="glLoadIdentity")
     end subroutine glLoadIdentity
     ! void glMatrixMode(GLenum mode);
     subroutine glMatrixMode(mode) bind(c,name="glMatrixMode")
       import :: GLenum
       integer(GLenum), value :: mode
     end subroutine glMatrixMode
     ! void glOrtho(GLdouble left, GLdouble right, GLdouble bottom,
     !              GLdouble top, GLdouble near_val, GLdouble far_val)
     subroutine glOrtho(left,right,bottom,top,near_val,far_val) bind(c,name="glOrtho")
       import :: GLdouble
       real(GLdouble), value :: left,right,bottom,top,near_val,far_val
     end subroutine glOrtho
     ! void glVertex3f(GLfloat x,GLfloat y,GLfloat z);
     subroutine glVertex3f(x,y,z) bind(c,name="glVertex3f")
       import :: GLfloat
       real(GLfloat), value :: x,y,z
     end subroutine glVertex3f
     ! void glViewport(GLint x,GLint y,GLsizei width,GLsizei height);
     subroutine glViewport(x,y,w,h) bind(c,name="glViewport")
       import :: GLint, GLsizei
       integer(GLint), value :: x, y
       integer(GLsizei), value :: w, h
     end subroutine glViewport
  end interface

contains

end module gui_gl
