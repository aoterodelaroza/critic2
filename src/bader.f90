! This module contains Bader integration-on-a-grid as proposed by
! Henkelman and collaborators. The code in this module has been
! adapted from the 'bader' program (version 0.28a, 07/12/12), that
! can be located at the Henkelman's group webpage: 
!    http://theory.cm.utexas.edu/
!    http://theory.cm.utexas.edu/vasp/bader/
! The authors of bader are Wenjie Tang, Andri Arnaldsson, Samuel T.
! Chill, and Graeme Henkelman.
! Also, see:
!   Comput. Mater. Sci. 36, 254-360 (2006).
!   J. Comput. Chem. 28, 899-908 (2007).
!   J. Phys.: Condens. Matter 21, 084204 (2009)

! Copyright 2009
! Wenjie Tang, Andri Arnaldsson, Samuel T. Chill, and Graeme Henkelman
!
! Bader is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.

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
module bader
  use systemmod, only: system
  implicit none

  private

  public :: bader_integrate
  
  interface
     module subroutine bader_integrate(s,bas)
       use types, only: basindat
       type(system), intent(inout) :: s
       type(basindat), intent(inout) :: bas
     end subroutine bader_integrate
  end interface

end module bader
