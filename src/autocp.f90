! Copyright (c) 2015-2022 Alberto Otero de la Roza <aoterodelaroza@gmail.com>,
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

!> Tools for automatic or manual critical point determination.
module autocp
  implicit none

  private

  public :: autocritic
  public :: autocritic_graph
  public :: cpreport

  interface
     module subroutine autocritic(line,success,clear,seeds)
       character*(*), intent(in) :: line
       logical, intent(out), optional :: success
       logical, intent(in), optional :: clear
       real*8, allocatable, intent(inout), optional :: seeds(:,:)
     end subroutine autocritic
     module subroutine autocritic_graph()
     end subroutine autocritic_graph
     module subroutine cpreport(line)
       character*(*), intent(in) :: line
     end subroutine cpreport
  end interface

end module autocp

