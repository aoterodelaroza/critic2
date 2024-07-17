! Copyright (c) 2007-2022 Alberto Otero de la Roza <aoterodelaroza@gmail.com>,
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

! Routines for reading and handling molecular and crystal vibrations.
submodule (crystalmod) edit
  implicit none

contains

  !> Read a file containing the vibrational information for this
  !> system. Populates c%vib.
  module subroutine read_vibrations_file(c,file,errmsg,ti)
    use param, only: isformat_unknown
    use crystalseedmod, only: vibrations_detect_format
    class(crystal), intent(inout) :: c
    character*(*), intent(in) :: file
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    integer :: ivformat

    ! detect the format
    call vibrations_detect_format(file,ivformat,ti)
    if (ivformat == isformat_unknown) then
       errmsg = "Unknown vibration file format: " // trim(file)
       return
    end if

    write (*,*) "bleh!", ivformat
    stop 1

  end subroutine read_vibrations_file

end submodule edit
