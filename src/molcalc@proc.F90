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
    use tools_io, only: ferror, faterr, uout, isexpression_or_word, lower, equal,&
       lgetword, getword
    character*(*), intent(inout) :: line

    character(len=:), allocatable :: word, lword, expr, savevar
    integer :: lp
    logical :: ok

    write (uout,'("* MOLCALC: calculations using molecular meshes and wavefunctions ")')

    savevar = ""
    lp = 1
    ok =  isexpression_or_word(expr,line,lp)
    if (.not.ok) then
       call ferror('molcalc_driver','Wrong syntax in MOLCALC',faterr,syntax=.true.)
       return
    end if

    do while (.true.)
       word = lgetword(line,lp)
       if (equal(word,'assign')) then
          savevar = getword(line,lp)
          if (len_trim(savevar) == 0) then
             call ferror('molcalc_driver','Zero-length variable name in MOLCALC/ASSIGN',faterr,line,syntax=.true.)
             return
          end if
       else if (len_trim(word) > 0) then
          call ferror('molcalc_driver','Unknown keyword',faterr,line,syntax=.true.)
          return
       else
          exit
       end if
    end do

    lword = lower(expr)
    if (equal(lword,"nelec")) then
       call molcalc_nelec()
    elseif (equal(lword,"peach")) then
       call molcalc_peach()
    elseif (equal(lword,"integral")) then
       call molcalc_integral()
    else
       call molcalc_expression(expr,savevar)
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

  !> Compute an expression in the molecular mesh. Save the result in variable
  !> savevar.
  subroutine molcalc_expression(expr,savevar)
    use systemmod, only: sy
    use meshmod, only: mesh
    use arithmetic, only: setvariable
    use tools_io, only: string, uout
    use types, only: scalar_value
    character*(*), intent(in) :: expr, savevar
    
    type(mesh) :: m
    real*8, allocatable :: ff(:)
    real*8 :: fval, fsum
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
    fsum = sum(ff * m%w)
    write (uout,'("+ Integral(",A,") = ",A/)') string(expr), string(fsum,'f',14,8)
    deallocate(ff)
    if (len_trim(savevar) > 0) call setvariable(trim(savevar),fsum)

  end subroutine molcalc_expression

  subroutine molcalc_integral()
    use tools_io, only: ferror, faterr
#ifdef HAVE_CINT

    integer :: off
    integer :: natm, atm(6,2), nbas, bas(8,2)
    integer :: i, j, di, dj, is0, shls(4)
    real*8 :: env(1000)
    real*8, allocatable :: buf1e(:,:,:)
    real*8, external :: CINTgto_norm
    integer, external :: CINTcgto_cart, CINT1e_kin_cart, CINT1e_ovlp_cart
    real*8 :: tmn(2,2), pmn(2,2), tkin

    off = 0
    natm = 2
    atm = 0
    env = 0d0
    nbas = 2
    bas = 0
    env = 0d0

    atm(1,1) = 1
    atm(2,1) = off
    atm(3,1) = 1
    atm(4,1) = off+3
    atm(5:6,1) = 0
    env(off+1:off+3) = (/0d0,0d0,0.5d0/) / 0.52917720859d0 ! in bohr!
    env(off+4) = 0d0
    off = off + 4

    atm(1,2) = 1
    atm(2,2) = off
    atm(3,2) = 1
    atm(4,2) = off+3
    atm(5:6,2) = 0
    env(off+1:off+3) = (/0d0,0d0,-0.5d0/) / 0.52917720859d0 ! in bohr!
    env(off+4) = 0d0
    off = off + 4

    ! CINTgto_norm is same as my gnorm but ---without the 4*pi---!
    ! S   3   1.00
    ! 3.42525091             0.15432897
    ! 0.62391373             0.53532814
    ! 0.16885540             0.44463454
    env(off+1) = 3.42525091d0
    env(off+2) = 0.62391373d0
    env(off+3) = 0.16885540d0
    env(off+4) = 0.15432897d0 * CINTgto_norm(0,env(off+1))
    env(off+5) = 0.53532814d0 * CINTgto_norm(0,env(off+2))
    env(off+6) = 0.44463454d0 * CINTgto_norm(0,env(off+3))

    bas(1,1) = 0
    bas(2,1) = 0
    bas(3,1) = 3
    bas(4,1) = 1
    bas(5,1) = 0
    bas(6,1) = off
    bas(7,1) = off+3
    bas(8,1) = 0

    bas(1,2) = 1
    bas(2,2) = 0
    bas(3,2) = 3
    bas(4,2) = 1
    bas(5,2) = 0
    bas(6,2) = off
    bas(7,2) = off+3
    bas(8,2) = 0

    do i = 1, nbas
       di = CINTcgto_cart(i-1, bas)
       do j = 1, i
          dj = CINTcgto_cart(j-1, bas)
          allocate(buf1e(di,di,1))
          shls(1) = i-1
          shls(2) = j-1
          is0 = CINT1e_kin_cart(buf1e,shls,atm,natm,bas,nbas,env)          
          ! is0 = CINT1e_ovlp_cart(buf1e,shls,atm,natm,bas,nbas,env)          
          tmn(i,j) = buf1e(1,1,1)
          deallocate(buf1e)
       end do
    end do
    tmn(1,2) = tmn(2,1)

    pmn = (5.78027982d-01 * sqrt(2d0))**2 ! this matches pyscf's make_rdm1

    tkin = pmn(1,1)*tmn(1,1) + pmn(1,2)*tmn(1,2) + pmn(2,1)*tmn(2,1) + pmn(2,2)*tmn(2,2)
    write (*,*) pmn(1,1), tmn(1,1), pmn(1,1)*tmn(1,1)
    write (*,*) pmn(1,2), tmn(1,2), pmn(1,2)*tmn(1,2)
    write (*,*) pmn(2,1), tmn(2,1), pmn(2,1)*tmn(2,1)
    write (*,*) pmn(2,2), tmn(2,2), pmn(2,2)*tmn(2,2)
    write (*,*) tkin

    stop 1

#else
    call ferror("molcalc_integral","INTEGRAL requires the CINT library",faterr)
#endif
  end subroutine molcalc_integral

end submodule proc
