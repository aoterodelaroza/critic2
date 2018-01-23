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

! This module is based on the code kindly provided by Enrico Benassi,
! Scuola Normale Superiore di Pisa, Italy <ebenassi3@gmail.com>

! Calculate the STM map of the surface slab in the input using the
! Tersoff-Hamann approximation and either a constant height or a constant
! current. The reference field is assumed to be the DOS at the Fermi level
! (although the density will also give reasonable results). 
module stm
  implicit none

  private

  public :: stm_driver
  
  interface
     module subroutine stm_driver(line)
       character*(*), intent(inout) :: line
     end subroutine stm_driver
  end interface
  
end module stm
