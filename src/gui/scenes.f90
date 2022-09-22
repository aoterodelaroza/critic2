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
module scenes
  use iso_c_binding
  implicit none

  private

  type scene
     logical :: isinit = .false. ! whether the scene has been initialized
     integer :: id ! system ID
     logical :: showcell = .false. ! show the unit cell?
     integer :: ncell(3) ! number of cells show in each direction
     integer :: imotif = 0 ! 0 = none, 1 = border, 2 = onemotif, 3 = molmotif
     real*8 :: scenerad = 1d0 ! scene radius
   contains
     procedure :: init => scene_init
     procedure :: end => scene_end
     procedure :: render => scene_render
  end type scene
  public :: scene

  ! module procedure interfaces
  interface
     module subroutine scene_init(s,isys)
       class(scene), intent(inout), target :: s
       integer, intent(in) :: isys
     end subroutine scene_init
     module subroutine scene_end(s)
       class(scene), intent(inout), target :: s
     end subroutine scene_end
     module subroutine scene_render(s)
       class(scene), intent(inout), target :: s
     end subroutine scene_render
  end interface

end module scenes

