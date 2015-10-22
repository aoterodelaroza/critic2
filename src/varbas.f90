! Copyright (c) 2015 Alberto Otero de la Roza <alberto@fluor.quimica.uniovi.es>,
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

!> External density-independent global variables and tools.
module varbas
  use global, only: mneq
  use types, only: cp_type
  implicit none

  public

  ! is the list initialized?
  logical :: cplist_init = .false.
  
  ! Non-equivalent critical points
  integer, parameter :: mneqcp0 = 100
  integer :: ncp = 0
  type(cp_type), allocatable :: cp(:)

  ! Compelte list of critical points
  integer, parameter :: mcp0 = 1000
  integer :: ncpcel = 0
  type(cp_type), allocatable :: cpcel(:)

  real*8, parameter :: Rbetadef  = 0.5d0 !< default beta sphere radius

contains

  !> Deallocate the arrays for the critical points
  subroutine varbas_end()
    if (allocated(cp)) deallocate(cp)
    if (allocated(cpcel)) deallocate(cpcel)
    ncp = 0
    ncpcel = 0
    cplist_init = .false.
  end subroutine varbas_end

  !> Given the point xp in crystallographic coordinates, calculates
  !> the nearest CP of type 'type' or non-equivalent index 'idx'. In the
  !> output, nid represents the id (complete CP list), dist is the
  !> distance. If nozero is used, skip zero distance CPs.
  subroutine nearestcp(xp,nid,dist,type,idx,nozero)
    use struct_basic
    implicit none

    real*8, intent(in) :: xp(:)
    integer, intent(out) :: nid
    real*8, intent(out) :: dist
    integer, intent(in), optional :: type
    integer, intent(in), optional :: idx
    logical, intent(in), optional :: nozero

    real*8, parameter :: eps2 = 1d-10 * 1d-10

    real*8 :: temp(3), d2, d2min
    integer :: j

    ! check if it is a known cp
    d2min = 1d30
    do j = 1, ncpcel
       if (present(type)) then
          if (cpcel(j)%typ /= type) cycle
       end if
       if (present(idx)) then
          if (cpcel(j)%idx /= idx) cycle
       end if
       temp = cpcel(j)%x - xp
       call cr%shortest(temp,d2)
       if (present(nozero)) then
          if (d2 < eps2) cycle
       end if
       if (d2 < d2min) then
          nid = j
          d2min = d2
       end if
    end do
    dist = sqrt(d2min)

  end subroutine nearestcp

  !> Given a position xp (usually the end-point of a gradient path),
  !> identify the critical point it represents and return the
  !> complete CP list id. If not found, 0. eps is the distance
  !> criterion (bohr). xp is in crystallographic coordinates.
  function whichcp(xp,eps)
    use struct_basic
    use tools_io

    real*8, intent(in) :: xp(3)
    real*8, intent(in) :: eps
    integer :: whichcp

    real*8 :: x(3), eps2, d2
    integer :: j

    ! check if it is a known cp
    eps2 = eps * eps
    do j = 1, ncpcel
       x = cpcel(j)%x - xp
       call cr%shortest(x,d2)
       ! is it this cp?
       if (d2 < eps2) then
          whichcp = j
          return 
       end if
    end do
    whichcp = 0

  end function whichcp

  subroutine varbas_identify(line0,lp)
    use struct_basic
    use global
    use tools_io
    use types
    use param
    
    character*(*), intent(in) :: line0
    integer, intent(inout) :: lp

    logical :: ok
    character(len=:), allocatable :: line, word
    real*8 :: x0(3), xmin(3), xmax(3)
    real*8, allocatable :: pointlist(:,:)
    logical, allocatable :: isrec(:)
    integer :: i, j, n, lu, nat, idx
    logical :: found

    integer :: ldunit, unit, mm
    integer, parameter :: unit_au = 1
    integer, parameter :: unit_ang = 2
    integer, parameter :: unit_x = 3
    integer, parameter :: unit_rec = 4

    real*8, parameter :: eps = 1d-4

    ! parse the first word
    ldunit = unit_x
    word = lgetword(line0,lp)
    if (equal(word,'angstrom') .or.equal(word,'ang')) then
       ldunit = unit_ang
    elseif (equal(word,'bohr') .or.equal(word,'au')) then
       ldunit = unit_au
    elseif (equal(word,'cryst')) then
       ldunit = unit_x
    elseif (equal(word,'reciprocal')) then
       ldunit = unit_rec
    elseif (len_trim(word) > 0) then
       call ferror('varbas_identify','Unkwnon extra keyword',faterr,line0)
    endif
    word = lgetword(line0,lp)
    if (len_trim(word) > 0) &
       call ferror('varbas_identify','Unkwnon extra keyword',faterr,line0)

    lp = 1
    n = 0
    allocate(pointlist(3,10),isrec(10))
    ok = getline(uin,line,ucopy=ucopy)
    if (ok) then
       word = lgetword(line,lp)
    else
       word = ""
       lp = 1
    end if
    do while (ok.and..not.equal(word,'endidentify').and..not.equal(word,'end'))
       lp = 1
       ok = eval_next (x0(1), line, lp)
       ok = ok .and. eval_next (x0(2), line, lp)
       ok = ok .and. eval_next (x0(3), line, lp)
       if (ok) then
          ! this is a point, parse the last word
          word = lgetword(line,lp)
          if (equal(word,'angstrom') .or.equal(word,'ang')) then
             unit = unit_ang
          elseif (equal(word,'bohr') .or.equal(word,'au')) then
             unit = unit_au
          elseif (equal(word,'cryst')) then
             unit = unit_x
          elseif (equal(word,'reciprocal')) then
             unit = unit_rec
          else
             unit = ldunit
          endif

          if (unit == unit_ang) then
             x0 = cr%c2x(x0 / bohrtoa)
          elseif (unit == unit_au) then
             x0 = cr%c2x(x0)
          endif
          n = n + 1
          if (n > size(pointlist,2)) then
             call realloc(pointlist,3,2*n)
             call realloc(isrec,2*n)
          end if
          pointlist(:,n) = x0
          isrec(n) = (unit == unit_rec)
       else
          ! this is an xyz file
          lu = fopen_read(word)
          read(lu,*) nat
          read(lu,*) 
          do i = 1, nat
             read(lu,*) word, x0
             x0 = cr%c2x(x0 / bohrtoa)
             n = n + 1
             if (n > size(pointlist,2)) then
                call realloc(pointlist,3,2*n)
                call realloc(isrec,2*n)
             end if
             pointlist(:,n) = x0
             isrec(n) = (unit == unit_rec)
          end do
          call fclose(lu)
       endif

       ! read next line
       lp = 1
       ok = getline(uin,line,ucopy=ucopy) 
       if (ok) then
          word = lgetword(line,lp)
       else
          line = ""
          lp = 1
       end if
    enddo
    word = lgetword(line,lp)
    if (len_trim(word) > 0) &
       call ferror('varbas_identify','Unkwnon extra keyword',faterr,line)

    xmin = 1d40
    xmax = -1d40
    found = .false.
    ! identify the atoms
    write(uout,'("* IDENTIFY: match input coordinates to atoms or CPs in the crystal")')
    write(uout,'("# (x,y,z) is the position in crystallographic coordinates ")')
    write(uout,'("# Sym is the atomic symbol/CP identifier (n=ncp,b=bcp,r=rcp,c=ccp) ")')
    write(uout,'("# N-neq is the index for the atom/CP in the non-equivalent list ")')
    write(uout,'("# Cel is the index for the atom/CP in the unit cell list ")')
    write(uout,'("# n        x             y             z     mult name nneq  cel")')

    do i = 1, n
       x0 = pointlist(:,i)
       if (.not.isrec(i)) then
          idx = identify_cp(x0,eps)
          mm = cr%get_mult(x0)
          if (idx > 0) then
             write (uout,'(99(A,X))') string(i,length=3,justify=ioj_left), &
                (string(x0(j),'f',length=13,decimal=8,justify=4),j=1,3), &
                string(mm,length=3,justify=ioj_center), &
                string(cpcel(idx)%name,length=5,justify=ioj_center), &
                string(cpcel(idx)%idx,length=4,justify=ioj_center), &
                string(idx,length=4,justify=ioj_center)
             do j = 1, 3
                xmin(j) = min(xmin(j),x0(j))
                xmax(j) = max(xmax(j),x0(j))
                found = .true.
             end do
          else
             write (uout,'(99(A,X))') string(i,length=3,justify=ioj_left), &
                (string(x0(j),'f',length=13,decimal=8,justify=4),j=1,3), &
                string(mm,length=3,justify=ioj_center), &
                string(" --- not found --- ")
          endif
       else
          mm = cr%get_mult_reciprocal(x0)
          write (uout,'(99(A,X))') string(i,length=3,justify=ioj_left), &
             (string(x0(j),'f',length=13,decimal=8,justify=4),j=1,3), &
             string(mm,length=3,justify=ioj_center), &
             string(" --- not found --- ")
       endif
    end do
    deallocate(pointlist,isrec)

    if (found) then
       write(uout,'("#")')
       write(uout,'("+ Cube, x0 (cryst): ",3(A,X))') (string(xmin(j),'f',decimal=8),j=1,3)
       write(uout,'("+ Cube, x1 (cryst): ",3(A,X))') (string(xmax(j),'f',decimal=8),j=1,3)
       xmin = cr%c2x(cr%x2c(xmin) - 2)
       xmax = cr%c2x(cr%x2c(xmax) + 2)
       write(uout,'("+ Cube + 2 bohr, x0 (cryst): ",3(A,X))') (string(xmin(j),'f',decimal=8),j=1,3)
       write(uout,'("+ Cube + 2 bohr, x1 (cryst): ",3(A,X))') (string(xmax(j),'f',decimal=8),j=1,3)
       xmin = cr%c2x(cr%x2c(xmin) - 3)
       xmax = cr%c2x(cr%x2c(xmax) + 3)
       write(uout,'("+ Cube + 5 bohr, x0 (cryst): ",3(A,X))') (string(xmin(j),'f',decimal=8),j=1,3)
       write(uout,'("+ Cube + 5 bohr, x1 (cryst): ",3(A,X))') (string(xmax(j),'f',decimal=8),j=1,3)
       write(uout,*)
    end if

  end subroutine varbas_identify

  !> Identify a CP in the unit cell. Input: cryst coords. Output:
  !> the cell CP index.
  function identify_cp(x0,eps)
    use struct_basic
    use tools_io

    integer :: identify_cp
    real*8, intent(in) :: x0(3)
    real*8, intent(in) :: eps

    real*8 :: x(3), dist2, eps2
    integer :: i

    identify_cp = 0
    eps2 = eps*eps
    do i = 1, ncpcel
       x = x0 - cpcel(i)%x
       call cr%shortest(x,dist2)
       if (dist2 < eps2) then
          identify_cp = i
          return
       end if
    end do

  end function identify_cp

end module varbas
