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

!> merge-sort routines contributed by:
!> Jose Luis Casals Sainz <joseluiscasalssainz@gmail.com>

!Toolbox with misc. tools.
module tools
  implicit none

  private

  !xx! sorting and uniqing
  public :: mergesort, imergesort, qcksort, iqcksort, uniqc
  !xx! system
  public :: unlink
  
  interface
     module subroutine qcksort (arr, iord, first, last)
       real*8, dimension(:), intent(in) :: arr
       integer, dimension(:), intent(inout) :: iord
       integer, intent(in) :: first
       integer, intent(in) :: last 
     end subroutine qcksort
     module subroutine iqcksort (iarr, iord, first, last)
       integer, dimension(:), intent(in) :: iarr
       integer, dimension(:), intent(inout) :: iord
       integer, intent(in) :: first
       integer, intent(in) :: last 
     end subroutine iqcksort
     module subroutine uniqc(v, nini, nfin, eps)
       integer, intent(in) :: nini
       integer, intent(inout) :: nfin
       real*8, intent(inout), dimension(3,*) :: v
       real*8, intent(in) :: eps
     end subroutine uniqc
     module subroutine mergesort(xvalt, irngt)
       real*8, dimension (:), intent (in) :: xvalt
       integer, dimension (:), intent (out) :: irngt
     end subroutine mergesort
     module subroutine imergesort(xvalt, irngt)
       integer, dimension (:), intent (in) :: xvalt
       integer, dimension (:), intent (out) :: irngt
     end subroutine imergesort
     module subroutine unlink(file)
       character(len=*), intent(in) :: file
     end subroutine unlink
  end interface

end module tools
