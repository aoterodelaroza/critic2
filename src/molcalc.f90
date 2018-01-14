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

! Calculations using molecular wavefunctions.
module molcalc
  implicit none

  private

  public :: molcalc_driver
  private :: molcalc_nelec
  private :: molcalc_peach
  private :: molcalc_integral
  
  interface
     module subroutine molcalc_driver(line)
       character*(*), intent(inout) :: line
     end subroutine molcalc_driver
     module subroutine molcalc_nelec()
     end subroutine molcalc_nelec
     module subroutine molcalc_peach()
     end subroutine molcalc_peach
     module subroutine molcalc_integral()
     end subroutine molcalc_integral
  end interface

end module molcalc
