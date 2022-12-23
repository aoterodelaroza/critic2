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
  integer, parameter, public :: nmaxsph = 5
  integer(c_int), target, public :: sphVAO(nmaxsph) ! sphere: vertex array object
  integer(c_int), target, public :: sphVBO          ! sphere: vertex buffer object
  integer(c_int), target, public :: sphEBO(nmaxsph) ! sphere: element buffer object

  ! sphere dimensions
  integer(c_int), parameter, public :: sphnve(nmaxsph) = (/12,42,162,642,2562/)
  integer(c_int), parameter, public :: sphnel(nmaxsph) = (/20,80,320,1280,5120/)
  integer(c_int), parameter, public :: sphneladd(0:nmaxsph) = (/0,20,100,420,1700,6820/)

  ! cylinder objects
  integer, parameter, public :: nmaxcyl = 3
  integer(c_int), target, public :: cylVAO(nmaxcyl) ! cylinder: vertex array object
  integer(c_int), target, public :: cylVBO(nmaxcyl) ! cylinder: vertex buffer object
  integer(c_int), target, public :: cylEBO(nmaxcyl) ! cylinder: element buffer object

  ! cylinder dimensions
  integer(c_int), parameter, public :: cylnve(nmaxcyl) = (/14,26,38/)
  integer(c_int), parameter, public :: cylnveadd(0:nmaxcyl) = (/0,14,40,78/)
  integer(c_int), parameter, public :: cylnel(nmaxcyl) = (/24,48,72/)
  integer(c_int), parameter, public :: cylneladd(0:nmaxcyl) = (/0,24,72,144/)

  ! text objects
  integer(c_int), target, public :: textVAO
  integer(c_int), target, public :: textVBO
  integer(c_int), target, public :: textVAOos
  integer(c_int), target, public :: textVBOos
  integer(c_int), parameter, public :: text_maxvert = 6 * 128

  ! module procedure interfaces
  interface
     module subroutine shapes_init()
     end subroutine shapes_init
     module subroutine shapes_end()
     end subroutine shapes_end
  end interface

end module shapes

