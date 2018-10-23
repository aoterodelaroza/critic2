! Copyright (c) 2007-2018 Alberto Otero de la Roza <aoterodelaroza@gmail.com>,
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

submodule (molcalc) proc
  implicit none

  !xx! private procedures
  ! subroutine molcalc_nelec()
  ! subroutine molcalc_peach()
  ! subroutine molcalc_integral()
  ! subroutine molcalc_expression(expr)

contains

  !> Driver for molecular calculations
  module subroutine molcalc_driver(line)
    use tools_io, only: ferror, faterr, uout, isexpression_or_word, lower, equal
    character*(*), intent(inout) :: line

    character(len=:), allocatable :: word, lword
    integer :: lp
    logical :: ok

    write (uout,'("* MOLCALC: calculations using molecular meshes and wavefunctions ")')

    lp = 1

    ok =  isexpression_or_word(word,line,lp)
    lword = lower(word)
    if (.not.ok) then
       call ferror('molcalc_driver','Wrong syntax in MOLCALC',faterr,syntax=.true.)
       return
    else if (equal(lword,"nelec")) then
       call molcalc_nelec()
    elseif (equal(lword,"peach")) then
       call molcalc_peach()
    elseif (equal(lword,"integral")) then
       call molcalc_integral()
    else
       call molcalc_expression(word)
    end if

  end subroutine molcalc_driver

  !xx! private procedures

  subroutine molcalc_nelec()
    use systemmod, only: sy
    use meshmod, only: mesh
    use tools_io, only: string, uout
    use param, only: im_volume, im_rho
    
    type(mesh) :: m
    integer :: prop(2)

    call m%gen(sy%c)

    write (uout,'("+ Simple molecular integrals (NELEC)")')
    call m%report()

    prop(1) = im_volume
    prop(2) = im_rho
    call m%fill(sy%f(sy%iref),prop,.not.sy%c%ismolecule)

    write (uout,'("+ Volume (bohr^3) = ",A)') string(sum(m%f(:,1) * m%w),'f',14,8)
    write (uout,'("+ Number of electrons = ",A)') string(sum(m%f(:,2) * m%w),'f',14,8)
    write (uout,*)

  end subroutine molcalc_nelec

  subroutine molcalc_peach()
    use systemmod, only: sy
    use meshmod, only: mesh
    use fieldmod, only: type_wfn
    use tools_io, only: ferror, faterr, getline, uin, ucopy, string, isinteger, isreal,&
       lgetword, equal, uout
    use types, only: realloc

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

    write (uout,'("+ Measure of overlap between orbitals in an excitation (PEACH). ")')
    write (uout,'("  Please cite: Peach et al., J. Chem. Phys. 128 (2008) 044118.")')
    call m%report()
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

  subroutine molcalc_integral()
    use tools_io, only: ferror, faterr
#ifdef HAVE_CINT

    write (*,*) "bleh!"

#else
    call ferror("molcalc_integral","INTEGRAL requires the CINT library",faterr)
#endif
  end subroutine molcalc_integral

  subroutine molcalc_expression(expr)
    use systemmod, only: sy
    use meshmod, only: mesh
    use tools_io, only: string, uout
    use types, only: scalar_value
    use param, only: im_volume, im_rho
    character*(*), intent(in) :: expr
    
    type(mesh) :: m
    real*8, allocatable :: ff(:)
    real*8 :: fval
    integer :: i
    logical :: ok

    call m%gen(sy%c)

    write (uout,'("+ Molecular integral calculation")')
    write (uout,'("  Expression: ",A)') trim(expr)
    call m%report()

    allocate(ff(m%n))
    !$omp parallel do private(fval,ok)
    do i = 1, m%n
       fval = sy%eval(expr,.true.,ok,m%x(:,i))
       !$omp critical (save)
       ff(i) = fval
       !$omp end critical (save)
    end do
    !$omp end parallel do
    write (uout,'("+ Integral(",A,") = ",A/)') string(expr), string(sum(ff * m%w),'f',14,8)
    deallocate(ff)

  end subroutine molcalc_expression

end submodule proc
