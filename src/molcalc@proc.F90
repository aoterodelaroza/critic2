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
    use tools_io, only: ferror, faterr, uout, string
#ifdef HAVE_CINT

    integer :: ioff, joff, koff, loff
    integer :: nbas, nbast
    integer :: i, j, k, l, di, dj, dk, dl, is0, shls(4), i1, j1
    real*8, allocatable :: buf1e(:,:,:), buf2e(:,:,:,:,:)
    real*8, external :: CINTgto_norm
    integer, external :: CINTcgto_cart, CINT1e_kin_cart, CINT1e_ovlp_cart, CINT1e_nuc_cart
    integer, external :: CINT2e_cart
    real*8, allocatable :: hmn(:,:), pmn(:,:), smn(:,:), eri(:,:,:,:), jmn(:,:), kmn(:,:)
    real*8, allocatable :: vmn(:,:)
    real*8 :: ee, enuc, etot, dij

    nbas = sy%f(1)%wfn%cint%nbas
    nbast = sy%f(1)%wfn%cint%nbast
    allocate(hmn(nbast,nbast),smn(nbast,nbast),pmn(nbast,nbast))
    allocate(jmn(nbast,nbast),kmn(nbast,nbast),vmn(nbast,nbast))
    ! allocate(eri(nbast,nbast,nbast,nbast)) ! in-core - cannot be done for most systems

    ioff = 0
    hmn = 0d0
    smn = 0d0
    do i = 1, nbas
       di = CINTcgto_cart(i-1, sy%f(1)%wfn%cint%bas)
       joff = 0
       do j = 1, nbas
          dj = CINTcgto_cart(j-1, sy%f(1)%wfn%cint%bas)
          if (j >= i) then
             allocate(buf1e(di,dj,1))
             shls(1) = i-1
             shls(2) = j-1

             ! kinetic energy
             is0 = CINT1e_kin_cart(buf1e,shls,sy%f(1)%wfn%cint%atm,sy%f(1)%wfn%cint%natm,sy%f(1)%wfn%cint%bas,sy%f(1)%wfn%cint%nbas,sy%f(1)%wfn%cint%env)
             hmn(ioff+1:ioff+di,joff+1:joff+dj) = hmn(ioff+1:ioff+di,joff+1:joff+dj) + buf1e(:,:,1)

             ! nuclear attraction
             is0 = CINT1e_nuc_cart(buf1e,shls,sy%f(1)%wfn%cint%atm,sy%f(1)%wfn%cint%natm,sy%f(1)%wfn%cint%bas,sy%f(1)%wfn%cint%nbas,sy%f(1)%wfn%cint%env)
             hmn(ioff+1:ioff+di,joff+1:joff+dj) = hmn(ioff+1:ioff+di,joff+1:joff+dj) + buf1e(:,:,1)

             ! overlap
             is0 = CINT1e_ovlp_cart(buf1e,shls,sy%f(1)%wfn%cint%atm,sy%f(1)%wfn%cint%natm,sy%f(1)%wfn%cint%bas,sy%f(1)%wfn%cint%nbas,sy%f(1)%wfn%cint%env)
             smn(ioff+1:ioff+di,joff+1:joff+dj) = buf1e(:,:,1)

             ! propagate to the upper half
             hmn(joff+1:joff+dj,ioff+1:ioff+di) = transpose(hmn(ioff+1:ioff+di,joff+1:joff+dj))
             smn(joff+1:joff+dj,ioff+1:ioff+di) = transpose(smn(ioff+1:ioff+di,joff+1:joff+dj))

             deallocate(buf1e)
          end if
          joff = joff + dj
       end do
       ioff = ioff + di
    end do

    ! make the 1-dm
    pmn = matmul(transpose(sy%f(1)%wfn%cint%moc),sy%f(1)%wfn%cint%moc) * 2d0

    ! eri = 0d0
    jmn = 0d0
    kmn = 0d0
    ioff = 0
    do i = 1, nbas
       di = CINTcgto_cart(i-1, sy%f(1)%wfn%cint%bas)
       joff = 0
       do j = 1, nbas
          dj = CINTcgto_cart(j-1, sy%f(1)%wfn%cint%bas)
          koff = 0
          do k = 1, nbas
             dk = CINTcgto_cart(k-1, sy%f(1)%wfn%cint%bas)
             loff = 0
             do l = 1, nbas
                dl = CINTcgto_cart(l-1, sy%f(1)%wfn%cint%bas)

                allocate(buf2e(di,dj,dk,dl,1))
                shls = (/i-1,j-1,k-1,l-1/)
                is0 = CINT2e_cart(buf2e,shls,sy%f(1)%wfn%cint%atm,sy%f(1)%wfn%cint%natm,sy%f(1)%wfn%cint%bas,sy%f(1)%wfn%cint%nbas,sy%f(1)%wfn%cint%env,0_8)

                ! eri(ioff+1:ioff+di,joff+1:joff+dj,koff+1:koff+dk,loff+1:loff+dl) = buf2e(:,:,:,:,1)
                do j1 = loff+1, loff+dl
                   do i1 = joff+1, joff+dj
                      kmn(i1,j1) = kmn(i1,j1) + sum(pmn(ioff+1:ioff+di,koff+1:koff+dk) * buf2e(:,i1-joff,:,j1-loff,1))
                   end do
                   do i1 = koff+1, koff+dk
                      jmn(i1,j1) = jmn(i1,j1) + sum(pmn(ioff+1:ioff+di,joff+1:joff+dj) * buf2e(:,:,i1-koff,j1-loff,1))
                   end do
                end do

                deallocate(buf2e)

                loff = loff + dl
             end do
             koff = koff + dk
          end do
          joff = joff + dj
       end do
       ioff = ioff + di
    end do

    ! calculate V
    vmn = jmn - 0.5d0 * kmn

    ! calculate energies
    ee = sum(pmn * (hmn + 0.5d0 * vmn))
    enuc = 0d0
    do i = 1, sy%c%ncel
       do j = i+1, sy%c%ncel
          dij = norm2(sy%c%at(i)%r - sy%c%at(j)%r)
          enuc = enuc + sy%c%spc(sy%c%at(i)%is)%z * sy%c%spc(sy%c%at(j)%is)%z / dij
       end do
    end do
    etot = enuc + ee

    ! energies
    write (uout,'("+ Total energy = ",A," Hartree")') string(etot,'f',decimal=10)
    write (uout,*)

#else
    call ferror("molcalc_integral","INTEGRAL requires the CINT library",faterr)
#endif
  end subroutine molcalc_integral

end submodule proc
