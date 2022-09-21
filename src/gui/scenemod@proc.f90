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

! Scene object and GL rendering utilities
submodule (scenemod) proc
  implicit none

contains

  ! Initialize a scene object associated with system isys.
  module subroutine scene_init(s,isys)
    use gui_main, only: nsys, sysc, sys_init, sys
    class(scene), intent(inout), target :: s
    integer, intent(in) :: isys

    if (isys < 1 .or. isys > nsys) return
    if (sysc(isys)%status /= sys_init) return

    s%id = isys
    s%showcell = .not.(sys(isys)%c%ismolecule)
    s%ncell = 1
    s%isinit = .true.

  end subroutine scene_init

  ! Terminate a scene object
  module subroutine scene_end(s)
    class(scene), intent(inout), target :: s

    s%isinit = .false.
    s%id = 0

  end subroutine scene_end

end submodule proc
