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

! #ifdef HAVE_CINT
!     write (*,*) "bleh!"
! #else
!     call ferror("molcalc_peach","PEACH requires the CINT library",faterr)
! #endif

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
       elseif (equal(word,"peach")) then
          call molcalc_peach()
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
    integer :: prop(1)

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

    prop(1) = im_rho
    call m%fill(sy%f(sy%iref),prop,.false.)

    write (uout,'("+ Number of electrons (NELEC) = ",A/)') string(sum(m%f(:,1) * m%w),'f',14,8)

  end subroutine molcalc_nelec

  subroutine molcalc_peach()
    use systemmod, only: sy
    use meshmod, only: mesh
    use fieldmod, only: type_wfn
    use tools_io, only: ferror, faterr, getline, uin, ucopy, string, isinteger, isreal,&
       lgetword, equal, uout
    use types, only: realloc
    use param, only: im_rho

    type(mesh) :: m
    integer :: i, n, lp
    logical :: ok
    real*8 :: lam, dden, oia
    character(len=:), allocatable :: line, word
    integer, allocatable :: imo1(:), imo2(:), prop(:)
    real*8, allocatable :: kk(:)

    if (.not.sy%c%ismolecule) then
       call ferror("molcalc_driver","MOLCALC can not be used with crystals",faterr,syntax=.true.)
       return
    end if
    if (sy%f(sy%iref)%type /= type_wfn) then
       call ferror("molcalc_driver","PEACH can be used with molecular wavefunctions only",faterr,syntax=.true.)
       return
    end if

    write (uout,'("  Measure of overlap between orbitals in an excitation. ")')
    write (uout,'("  Please cite: Peach et al., J. Chem. Phys. 128 (2008) 044118.")')
    allocate(imo1(10),imo2(10),kk(10))
    n = 0
    do while (getline(uin,line,ucopy=ucopy))
       lp = 1
       word = lgetword(line,lp)
       if (equal(word,'endmolcalc').or.equal(word,'end')) then
          exit
       end if
       n = n + 1
       if (n > ubound(imo1,1)) then
          call realloc(imo1,2*n)
          call realloc(imo2,2*n)
          call realloc(kk,2*n)
       end if

       ! first orbital
       ok = isinteger(imo1(n),word)
       if (.not.ok) goto 999

       ! second orbital
       word = lgetword(line,lp)
       if (equal(word,"->")) &
          word = lgetword(line,lp)
       ok = isinteger(imo2(n),word)
       if (.not.ok) goto 999

       ! coefficient
       word = lgetword(line,lp)
       ok = isreal(kk(n),word)
       if (.not.ok) goto 999
    enddo
    if (n == 0) then
       call ferror("molcalc_driver","No MOs in PEACH",faterr,syntax=.true.)
       return
    end if

    ! generate and fill the mesh
    call m%gen(sy%c)
    allocate(prop(2*n))
    do i = 1, n
       prop(i) = 100 + imo1(i)
       prop(n+i) = 100 + imo2(i)
    end do
    call m%fill(sy%f(sy%iref),prop,.false.)
    deallocate(prop)

    lam = 0d0
    dden = 0d0
    do i = 1, n
       oia = sum(abs(m%f(:,i)) * abs(m%f(:,n+i)) * m%w)
       lam = lam + kk(i) * kk(i) * oia
       dden = dden + kk(i) * kk(i)
    end do
    lam = lam / dden

    write (uout,'("+ PEACH = ",A)') string(lam,'f',8,3)
    write (uout,*)

    deallocate(imo1,imo2,kk)
    return
999 continue
    call ferror("molcalc_peach","error reading line: " // string(line),faterr)

  end subroutine molcalc_peach

end module molcalc
