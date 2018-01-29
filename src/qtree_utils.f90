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

!> Qtree, utilities.
module qtree_utils
  use param, only: mlen
  implicit none

  private
  public :: small_writetess
  public :: open_difftess
  public :: close_difftess
  public :: getkeast
  
  interface
     module subroutine small_writetess(roottess,otrm,trm)
       use qtree_basic, only: qtreei
       use param, only: mlen
       character*(*), intent(in) :: roottess
       integer, intent(in) :: otrm
       integer(qtreei), intent(in) :: trm(:,:)
     end subroutine small_writetess
     module subroutine open_difftess(roottess)
       use param, only: mlen
       character*(*), intent(in) :: roottess
     end subroutine open_difftess
     module subroutine close_difftess(roottess)
       use param, only: mlen
       character*(*), intent(in) :: roottess
     end subroutine close_difftess
     module subroutine getkeast()
     end subroutine getkeast
  end interface

end module qtree_utils
