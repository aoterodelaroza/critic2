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

! Calculations using molecular wavefunctions.
module molcalc
  implicit none

  private

  public :: molcalc_driver
  private :: molcalc_nelec

contains

  !> Driver for molecular calculations
  subroutine molcalc_driver(line)
    use systemmod, only: sy
    use tools_io, only: ferror, faterr, uout, lgetword, equal

    character*(*), intent(inout) :: line

    character(len=:), allocatable :: word
    integer :: lp

    if (.not.sy%c%ismolecule) then
       call ferror("molcalc_driver","MOLCALC can not be used with crystals",faterr,syntax=.true.)
       return
    end if

    write (uout,'("* MOLCALC: calculations using molecular wavefunctions ")')

    lp = 1
    do while(.true.)
       word = lgetword(line,lp)
       if (equal(word,"nelec")) then
          call molcalc_nelec()
       elseif (len_trim(word) > 0) then
          call ferror('molcalc_driver','Wrong syntax in MOLCALC',faterr,syntax=.true.)
          return
       else
          exit
       end if
    end do

  end subroutine molcalc_driver

  subroutine molcalc_nelec()
    use systemmod, only: sy
    use meshmod, only: mesh
    use tools_io, only: string, uout
    use param, only: im_rho
    
    type(mesh) :: m
    integer :: prop(1), id(1)

    call m%gen(sy%c)

    write (uout,'("+ Molecular integrals")')
    write (uout,'(2X,"mesh size      ",A)') string(m%n)
    if (m%type == 0) then
       write (uout,'(2X,"mesh type      Becke")')
    elseif (m%type == 1) then
       write (uout,'(2X,"mesh type      Franchini (small)")')
    elseif (m%type == 2) then
       write (uout,'(2X,"mesh type      Franchini (normal)")')
    elseif (m%type == 3) then
       write (uout,'(2X,"mesh type      Franchini (good)")')
    elseif (m%type == 4) then
       write (uout,'(2X,"mesh type      Franchini (very good)")')
    elseif (m%type == 5) then
       write (uout,'(2X,"mesh type      Franchini (excellent)")')
    end if

    id(1) = 1
    prop(1) = im_rho
    call m%fill(sy%f(sy%iref),id,prop,.false.)

    write (uout,'("+ Number of electrons (NELEC) = ",A/)') string(sum(m%f(:,1) * m%w),'f',14,8)

  end subroutine molcalc_nelec

end module molcalc
