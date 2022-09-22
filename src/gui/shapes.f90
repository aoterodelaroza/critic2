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

! OpenGL buffers for simple shapes
module shapes
  use iso_c_binding
  implicit none

  private

  public :: shapes_init
  public :: shapes_end

  ! sphere objects
  integer, parameter, public :: nmaxsph = 1
  integer(c_int), target, public :: sphVAO(nmaxsph) ! sphere: vertex array object
  integer(c_int), target, public :: sphVBO          ! sphere: vertex buffer object
  integer(c_int), target, public :: sphEBO(nmaxsph) ! sphere: element buffer object

  ! test objects
  integer(c_int), target, public :: testVAO ! test: vertex array object
  integer(c_int), target, public :: testVBO ! test: vertex buffer object

  ! module procedure interfaces
  interface
     module subroutine shapes_init()
     end subroutine shapes_init
     module subroutine shapes_end()
     end subroutine shapes_end
  end interface

end module shapes

