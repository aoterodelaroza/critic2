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

! Fortran interfaces for tinycthread
module gui_interfaces_threads
  use iso_c_binding
  implicit none

  public

  interface
     ! int mtx_init(mtx_t *mtx, int type);
     function mtx_init(mtx, type) bind(c,name="mtx_init")
       import c_int, c_ptr
       type(c_ptr), value :: mtx
       integer(c_int), value :: type
       integer(c_int) :: mtx_init
     end function mtx_init
  end interface

end module gui_interfaces_threads
