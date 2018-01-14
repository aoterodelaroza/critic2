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

!> Integration by recursive division of the IWS.
module qtree
  use global, only: INT_lebedev
  implicit none
  
  private

  public :: qtree_integration
  private :: qtree_setsph1
  private :: qtree_setsph2
  public :: qtree_setsphfactor
  
  interface
     module subroutine qtree_integration(lvl, plvl)
       integer, intent(in) :: lvl
       integer, intent(in) :: plvl
     end subroutine qtree_integration
     module subroutine qtree_setsph1(lvl,verbose)
       integer, intent(in) :: lvl
       logical, intent(in) :: verbose
     end subroutine qtree_setsph1
     module subroutine qtree_setsph2(verbose)
       logical, intent(in) :: verbose
     end subroutine qtree_setsph2
     module subroutine qtree_setsphfactor(line)
       character*(*), intent(in) :: line
     end subroutine qtree_setsphfactor
  end interface

end module qtree

