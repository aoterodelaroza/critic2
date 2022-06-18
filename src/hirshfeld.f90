! Copyright (c) 2015 Alberto Otero de la Roza <aoterodelaroza@gmail.com>,
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

! hirshfeld integration
module hirshfeld
  implicit none

  private

  public :: hirsh_driver
  public :: hirsh_grid
  public :: hirsh_weights

  interface
     module subroutine hirsh_driver()
     end subroutine hirsh_driver
     module subroutine hirsh_grid(s,bas)
       use systemmod, only: system
       use types, only: basindat
       type(system), intent(inout) :: s
       type(basindat), intent(inout) :: bas
     end subroutine hirsh_grid
     module subroutine hirsh_weights(s,bas,idb,w)
       use systemmod, only: system
       use types, only: basindat
       type(system), intent(inout) :: s
       type(basindat), intent(in) :: bas
       integer, intent(in) :: idb
       real*8, intent(out) :: w(:,:,:)
     end subroutine hirsh_weights
  end interface

end module hirshfeld
