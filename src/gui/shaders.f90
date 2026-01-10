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

! OpenGL shader programs
module shaders
  use iso_c_binding
  implicit none

  private

  integer, parameter, public :: shader_test = 1
  integer, parameter, public :: shader_simple = 2
  integer, parameter, public :: shader_phong = 3
  integer, parameter, public :: shader_text_direct = 4
  integer, parameter, public :: shader_text_onscene = 5
  integer, parameter, public :: shader_pickindex = 6
  integer, parameter, public :: shader_NUM = 6

  public :: shaders_init
  public :: shaders_end
  public :: useshader
  public :: get_uniform_location
  public :: setuniform_int
  public :: setuniform_float
  public :: setuniform_vec3
  public :: setuniform_vec4
  public :: setuniform_mat3
  public :: setuniform_mat4

  interface setuniform_float
     module procedure setuniform_float_realc
     module procedure setuniform_float_real8
  end interface setuniform_float
  interface setuniform_vec3
     module procedure setuniform_vec3_realc
     module procedure setuniform_vec3_real8
  end interface setuniform_vec3
  interface setuniform_vec4
     module procedure setuniform_vec4_realc
     module procedure setuniform_vec4_real8
  end interface setuniform_vec4
  interface setuniform_mat3
     module procedure setuniform_mat3_realc
     module procedure setuniform_mat3_real8
  end interface setuniform_mat3
  interface setuniform_mat4
     module procedure setuniform_mat4_realc
     module procedure setuniform_mat4_real8
  end interface setuniform_mat4

  ! module procedure interfaces
  interface
     module subroutine shaders_init()
     end subroutine shaders_init
     module subroutine shaders_end()
     end subroutine shaders_end
     module subroutine useshader(id)
       integer, intent(in) :: id
     end subroutine useshader
     module function get_uniform_location(name)
       character*(*), intent(in) :: name
       integer(c_int) :: get_uniform_location
     end function get_uniform_location
     module subroutine setuniform_int(x,name,idxi)
       integer(c_int), intent(in) :: x
       character*(*), intent(in), optional :: name
       integer(c_int), intent(in), optional :: idxi
     end subroutine setuniform_int
     module subroutine setuniform_float_realc(x,name,idxi)
       real(c_float), intent(in) :: x
       character*(*), intent(in), optional :: name
       integer(c_int), intent(in), optional :: idxi
     end subroutine setuniform_float_realc
     module subroutine setuniform_float_real8(x,name,idxi)
       real*8, intent(in) :: x
       character*(*), intent(in), optional :: name
       integer(c_int), intent(in), optional :: idxi
     end subroutine setuniform_float_real8
     module subroutine setuniform_vec3_realc(x,name,idxi)
       real(c_float), intent(in) :: x(3)
       character*(*), intent(in), optional :: name
       integer(c_int), intent(in), optional :: idxi
     end subroutine setuniform_vec3_realc
     module subroutine setuniform_vec3_real8(x,name,idxi)
       real*8, intent(in) :: x(3)
       character*(*), intent(in), optional :: name
       integer(c_int), intent(in), optional :: idxi
     end subroutine setuniform_vec3_real8
     module subroutine setuniform_vec4_realc(x,name,idxi)
       real(c_float), intent(in) :: x(4)
       character*(*), intent(in), optional :: name
       integer(c_int), intent(in), optional :: idxi
     end subroutine setuniform_vec4_realc
     module subroutine setuniform_vec4_real8(x,name,idxi)
       real*8, intent(in) :: x(4)
       character*(*), intent(in), optional :: name
       integer(c_int), intent(in), optional :: idxi
     end subroutine setuniform_vec4_real8
     module subroutine setuniform_mat3_realc(x,name,idxi)
       real(c_float), intent(in) :: x(3,3)
       character*(*), intent(in), optional :: name
       integer(c_int), intent(in), optional :: idxi
     end subroutine setuniform_mat3_realc
     module subroutine setuniform_mat3_real8(x,name,idxi)
       real*8, intent(in) :: x(3,3)
       character*(*), intent(in), optional :: name
       integer(c_int), intent(in), optional :: idxi
     end subroutine setuniform_mat3_real8
     module subroutine setuniform_mat4_realc(x,name,idxi)
       real(c_float), intent(in) :: x(4,4)
       character*(*), intent(in), optional :: name
       integer(c_int), intent(in), optional :: idxi
     end subroutine setuniform_mat4_realc
     module subroutine setuniform_mat4_real8(x,name,idxi)
       real*8, intent(in) :: x(4,4)
       character*(*), intent(in), optional :: name
       integer(c_int), intent(in), optional :: idxi
     end subroutine setuniform_mat4_real8
  end interface

end module shaders

