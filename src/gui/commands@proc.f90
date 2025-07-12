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

! The class to handle ImGui windows, general routines.
submodule (commands) proc
  implicit none

contains

  !> Deallocate all data in a command and reset it to emtpy
  module subroutine command_end(c)
    class(command_inout), intent(inout) :: c

    c%id = 0
    c%status = command_inout_empty
    c%size = 0
    c%scrolly = 0._c_float
    if (allocated(c%tooltipinfo)) deallocate(c%tooltipinfo)
    if (allocated(c%input)) deallocate(c%input)
    if (allocated(c%output)) deallocate(c%output)

  end subroutine command_end

  !> Clear all commands
  module subroutine clear_all_commands()

    integer :: i

    do i = 1, ncom
       call com(i)%end()
    end do
    ncomid = 0
    ncom = 0
    nicom = 0
    idcom = 0
    icom = 0

  end subroutine clear_all_commands

end submodule proc
