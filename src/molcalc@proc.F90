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
    use systemmod, only: sy
    use tools_io, only: ferror, faterr
    use param, only: pi
#ifdef HAVE_CINT

    integer :: off, ioff, joff
    integer :: natm, nbas, nbast
    integer :: i, j, di, dj, is0, shls(4)
    real*8 :: env(1000)
    integer, allocatable :: atm(:,:), bas(:,:)
    real*8, allocatable :: buf1e(:,:,:)
    real*8, external :: CINTgto_norm
    integer, external :: CINTcgto_cart, CINT1e_kin_cart, CINT1e_ovlp_cart
    real*8, allocatable :: tmn(:,:), pmn(:,:), smn(:,:), mocoef(:)
    real*8, allocatable :: aux(:,:)
    integer, allocatable :: imap(:), itype(:)
    real*8 :: norm1, norm2, aexp

    off = 20
    natm = 3
    env = 0d0
    nbas = 10
    env = 0d0

    allocate(atm(6,natm),bas(8,nbas))
    atm = 0
    bas = 0

    ! atoms
    atm(1,1) = 8
    atm(2,1) = off
    atm(3,1) = 1
    atm(4,1) = off+3
    atm(5:6,1) = 0
    env(off+1:off+3) = (/0.d0,0.d0,0.1188819994d0/) / 0.52917720859d0 ! in bohr!
    env(off+4) = 0d0
    off = off + 4

    atm(1,2) = 1
    atm(2,2) = off
    atm(3,2) = 1
    atm(4,2) = off+3
    atm(5:6,2) = 0
    env(off+1:off+3) = (/0.d0,0.7566529960d0,-0.4755289978d0/) / 0.52917720859d0 ! in bohr!
    env(off+4) = 0d0
    off = off + 4

    atm(1,3) = 1
    atm(2,3) = off
    atm(3,3) = 1
    atm(4,3) = off+3
    atm(5:6,2) = 0
    env(off+1:off+3) = (/0.0000000000d0,-0.7566529960d0,-0.4755289978d0/) / 0.52917720859d0 ! in bohr!
    env(off+4) = 0d0
    off = off + 4
    
    ! coefficients and exponents
    ! CINTgto_norm is same as my gnorm but ---without the 4*pi---!
    ! -- oxygen --
    ! S   6   1.00
    !    5484.6717000              0.0018311
    !     825.2349500              0.0139501
    !     188.0469600              0.0684451
    !      52.9645000              0.2327143
    !      16.8975700              0.4701930
    !       5.7996353              0.3585209
    env(off+1)  = 5484.6717000d0
    env(off+2)  =  825.2349500d0
    env(off+3)  =  188.0469600d0
    env(off+4)  =   52.9645000d0
    env(off+5)  =   16.8975700d0
    env(off+6)  =    5.7996353d0
    env(off+7)  = 0.0018311d0 * CINTgto_norm(0,env(off+1))
    env(off+8)  = 0.0139501d0 * CINTgto_norm(0,env(off+2))
    env(off+9)  = 0.0684451d0 * CINTgto_norm(0,env(off+3))
    env(off+10) = 0.2327143d0 * CINTgto_norm(0,env(off+4))
    env(off+11) = 0.4701930d0 * CINTgto_norm(0,env(off+5))
    env(off+12) = 0.3585209d0 * CINTgto_norm(0,env(off+6))
    ! SP   3   1.00
    !      15.5396160             -0.1107775              0.0708743
    !       3.5999336             -0.1480263              0.3397528
    !       1.0137618              1.1307670              0.7271586
    env(off+13) = 15.5396160d0
    env(off+14) =  3.5999336d0
    env(off+15) =  1.0137618d0
    env(off+16) = -0.1107775d0 * CINTgto_norm(0,env(off+13))
    env(off+17) = -0.1480263d0 * CINTgto_norm(0,env(off+14))
    env(off+18) =  1.1307670d0 * CINTgto_norm(0,env(off+15))
    env(off+19) =  0.0708743d0 * CINTgto_norm(1,env(off+13))
    env(off+20) =  0.3397528d0 * CINTgto_norm(1,env(off+14))
    env(off+21) =  0.7271586d0 * CINTgto_norm(1,env(off+15))
    ! SP   1   1.00
    !       0.2700058              1.0000000              1.0000000
    env(off+22) = 0.2700058d0
    env(off+23) = 1d0 * CINTgto_norm(0,env(off+22))
    env(off+24) = 1d0 * CINTgto_norm(1,env(off+22))
    ! D   1   1.00
    !       0.8000000              1.0000000
    env(off+25) = 0.8d0
    env(off+26) = 1d0 * CINTgto_norm(2,env(off+25))
    ! -- hydrogen --
    ! S   3   1.00
    !      18.7311370              0.03349460
    !       2.8253937              0.23472695
    !       0.6401217              0.81375733
    env(off+27) = 18.7311370d0
    env(off+28) =  2.8253937d0
    env(off+29) =  0.6401217d0
    env(off+30) =  0.03349460d0 * CINTgto_norm(0,env(off+27))
    env(off+31) =  0.23472695d0 * CINTgto_norm(0,env(off+28))
    env(off+32) =  0.81375733d0 * CINTgto_norm(0,env(off+29))
    ! S   1   1.00
    !       0.1612778              1.0000000
    env(off+33) = 0.1612778d0
    env(off+34) = 1d0 * CINTgto_norm(0,env(off+33))

    ! shells
    bas(1,1) = 0
    bas(2,1) = 0
    bas(3,1) = 6
    bas(4,1) = 1
    bas(5,1) = 0
    bas(6,1) = off
    bas(7,1) = off+6
    bas(8,1) = 0

    bas(1,2) = 0
    bas(2,2) = 0
    bas(3,2) = 3
    bas(4,2) = 1
    bas(5,2) = 0
    bas(6,2) = off+12
    bas(7,2) = off+15
    bas(8,2) = 0

    bas(1,3) = 0
    bas(2,3) = 1
    bas(3,3) = 3
    bas(4,3) = 1
    bas(5,3) = 0
    bas(6,3) = off+12
    bas(7,3) = off+18
    bas(8,3) = 0

    bas(1,4) = 0
    bas(2,4) = 0
    bas(3,4) = 1
    bas(4,4) = 1
    bas(5,4) = 0
    bas(6,4) = off+21
    bas(7,4) = off+22
    bas(8,4) = 0

    bas(1,5) = 0
    bas(2,5) = 1
    bas(3,5) = 1
    bas(4,5) = 1
    bas(5,5) = 0
    bas(6,5) = off+21
    bas(7,5) = off+23
    bas(8,5) = 0

    bas(1,6) = 0
    bas(2,6) = 2
    bas(3,6) = 1
    bas(4,6) = 1
    bas(5,6) = 0
    bas(6,6) = off+24
    bas(7,6) = off+25
    bas(8,6) = 0

    bas(1,7) = 1
    bas(2,7) = 0
    bas(3,7) = 3
    bas(4,7) = 1
    bas(5,7) = 0
    bas(6,7) = off+26
    bas(7,7) = off+29
    bas(8,7) = 0

    bas(1,8) = 1
    bas(2,8) = 0
    bas(3,8) = 1
    bas(4,8) = 1
    bas(5,8) = 0
    bas(6,8) = off+32
    bas(7,8) = off+33
    bas(8,8) = 0

    bas(1,9) = 2
    bas(2,9) = 0
    bas(3,9) = 3
    bas(4,9) = 1
    bas(5,9) = 0
    bas(6,9) = off+26
    bas(7,9) = off+29
    bas(8,9) = 0

    bas(1,10) = 2
    bas(2,10) = 0
    bas(3,10) = 1
    bas(4,10) = 1
    bas(5,10) = 0
    bas(6,10) = off+32
    bas(7,10) = off+33
    bas(8,10) = 0

    nbast = 0
    do i = 1, nbas
       di = CINTcgto_cart(i-1, bas)
       nbast = nbast + di
    end do
    allocate(tmn(nbast,nbast),smn(nbast,nbast),pmn(nbast,nbast))

    ioff = 0
    tmn = 0d0
    smn = 0d0
    do i = 1, nbas
       di = CINTcgto_cart(i-1, bas)
       joff = 0
       do j = 1, nbas
          dj = CINTcgto_cart(j-1, bas)
          allocate(buf1e(di,dj,1))
          shls(1) = i-1
          shls(2) = j-1

          is0 = CINT1e_kin_cart(buf1e,shls,atm,natm,bas,nbas,env)          
          tmn(ioff+1:ioff+di,joff+1:joff+dj) = buf1e(:,:,1)

          is0 = CINT1e_ovlp_cart(buf1e,shls,atm,natm,bas,nbas,env)          
          smn(ioff+1:ioff+di,joff+1:joff+dj) = buf1e(:,:,1)

          deallocate(buf1e)
          joff = joff + dj
       end do
       ioff = ioff + di
    end do
    
    ! remap the moc matrix
    allocate(aux(size(sy%f(1)%wfn%bas%moc,1),nbast),imap(nbast),itype(nbast))
    ! fchk -> libcint
    ! 1 -> 1
    ! 2 -> 2
    ! 3-5 -> 3-5
    ! 6 -> 6
    ! 7-9 -> 7-9
    ! 10 xx -> 10 xx
    ! 11 yy -> 13 yy
    ! 12 zz -> 15 zz
    ! 13 xy -> 11 xy
    ! 14 xz -> 12 xz
    ! 15 yz -> 14 yz
    ! 16 -> 16
    ! 17 -> 17
    ! 18 -> 18
    ! 19 -> 19
    ! imap(fchk) = cint
    imap  = (/1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 13, 15, 11, 12, 14, 16, 17, 18, 19/)
    itype = (/1, 1, 2, 3, 4, 1, 2, 3, 4, 5,   6,  7,  8,  9, 10,  1,  1,  1,  1/)
    do i = 1, nbast
       aexp = 0.8d0
       if (itype(i) >= 5 .and. itype(i) <= 7) then
          norm1 = 2**(11d0/4d0) * aexp**(7d0/4d0) / pi**(3d0/4d0) / sqrt(3d0) ! critic
          norm2 = CINTgto_norm(2,aexp) ! cint
       else if (itype(i) >= 8 .and. itype(i) <= 10) then
          norm1 = 2**(11d0/4d0) * aexp**(7d0/4d0) / pi**(3d0/4d0) ! critic
          norm2 = CINTgto_norm(2,aexp) ! cint
       else
          norm1 = 1d0
          norm2 = 1d0
       end if
       aux(:,imap(i)) = sy%f(1)%wfn%bas%moc(:,i) * norm1 / norm2
    end do

    ! make the 1-dm
    pmn = matmul(transpose(aux),aux) * 2d0

    ! write (*,*) pmn(:,1)
    ! write (*,*) smn(:,1)
    ! write (*,*) tmn(:,1)

    write (*,*) sum(pmn * tmn)
    
    ! tkin = pmn(1,1)*tmn(1,1) + pmn(1,2)*tmn(1,2) + pmn(2,1)*tmn(2,1) + pmn(2,2)*tmn(2,2)
    ! write (*,*) pmn(1,1), tmn(1,1), pmn(1,1)*tmn(1,1)
    ! write (*,*) pmn(1,2), tmn(1,2), pmn(1,2)*tmn(1,2)
    ! write (*,*) pmn(2,1), tmn(2,1), pmn(2,1)*tmn(2,1)
    ! write (*,*) pmn(2,2), tmn(2,2), pmn(2,2)*tmn(2,2)
    ! write (*,*) tkin

    stop 1

#else
    call ferror("molcalc_integral","INTEGRAL requires the CINT library",faterr)
#endif
  end subroutine molcalc_integral

end submodule proc
