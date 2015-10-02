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

!> Ewald summations
module ewald
  implicit none

  private

  public :: ewald_energy, ewald_pot

contains

  !> Calculates the Ewald electrostatic energy, using the input charges.
  subroutine ewald_energy(ewe) 
    use struct_basic
    use tools_io

    real*8, intent(out) :: ewe

    real*8 :: qsum, q2sum
    real*8 :: rcut, hcut, eta
    integer :: lrmax(3), lhmax(3)
    real*8 :: x(3)
    integer :: i

    call sum_charges(qsum,q2sum)
    call calculate_cutoff(q2sum,rcut,hcut,eta,lrmax,lhmax)
    
    ewe = 0d0
    do i = 1, cr%nneq
       x = cr%at(i)%x
       ewe = ewe + cr%at(i)%mult * cr%at(i)%qat * &
          ewald_pot(x,.true.,rcut,hcut,eta,lrmax,lhmax,qsum,q2sum)
    end do
    ewe = ewe / 2d0

  end subroutine ewald_energy

  !> calculate the Ewald electrostatic potential at 
  !> an arbitrary position x (crystallographic coords.)
  !> If x is the nucleus j, return pot - q_j / |r-rj| at rj.
  function ewald_pot(x,isnuc,rcut,hcut,eta,lrmax,lhmax,qsum,q2sum)
    use struct_basic
    use param

    real*8, intent(in) :: x(3)
    logical, intent(in) :: isnuc
    real*8, intent(in) :: rcut, hcut, eta
    integer, intent(in) :: lrmax(3), lhmax(3)
    real*8, intent(in) :: qsum, q2sum
    real*8 :: ewald_pot

    real*8 :: nuc_cutoff2 = 1d-10

    real*8 :: rcut2
    integer :: i, i1, i2, i3, qnuc
    real*8 :: px(3), lvec(3), d2, d, dh
    real*8 :: sfac_c, sfacp, bbarg
    real*8 :: sum_real, sum_rec, sum0, sum_back

    ! is this a nuclear position? -> get charge
    qnuc = 0
    if (isnuc) then
       do i = 1, cr%ncel
          px = cr%atcel(i)%x - x
          d2 = dot_product(px,matmul(cr%gtensor,px))
          if (d2 < nuc_cutoff2) then
             qnuc = cr%at(cr%atcel(i)%idx)%qat
             exit
          end if
       end do
    end if

    ! real space sum
    rcut2 = rcut * rcut
    sum_real = 0
    do i1 = -lrmax(1),lrmax(1)
       do i2 = -lrmax(2),lrmax(2)
          do i3 = -lrmax(3),lrmax(3)
             lvec = real((/i1,i2,i3/),8)
             do i = 1,cr%ncel
                px = x - cr%atcel(i)%x - lvec
                d2 = dot_product(px,matmul(cr%gtensor,px))
                if (d2 < 1d-12 .or. d2 > rcut2) cycle
                d = sqrt(d2) / eta
                sum_real = sum_real + cr%at(cr%atcel(i)%idx)%qat * erfc(d) / d
             end do
          end do
       end do
    end do
    sum_real = sum_real / eta

    ! reciprocal space sum
    sum_rec = 0
    do i1 = -lhmax(1),lhmax(1)
       do i2 = -lhmax(2),lhmax(2)
          do i3 = -lhmax(3),lhmax(3)
             lvec = tpi * (/i1,i2,i3/)
             dh = sqrt(dot_product(lvec,matmul(cr%grtensor,lvec)))
             if (dh < 1d-12 .or. dh > hcut) cycle
             bbarg = 0.5d0 * dh * eta

             sfac_c = 0
             do i = 1, cr%ncel
                sfac_c = sfac_c + cr%at(cr%atcel(i)%idx)%qat * &
                   cos(dot_product(lvec,x-cr%atcel(i)%x))
             end do
             sfacp = 2d0 * sfac_c

             sum_rec = sum_rec + sfacp / dh**2 * exp(-bbarg**2)
          end do
       end do
    end do
    sum_rec = sum_rec * 2d0 * pi / cr%omega
    
    ! h = 0 term, apply only at the nucleus
    if (isnuc) then
       sum0 = - 2d0 * qnuc / sqpi / eta
    else
       sum0 = 0d0
    end if

    ! compensating background charge term
    sum_back = -qsum * eta**2 * pi / cr%omega 

    ! sum up and exit
    ewald_pot = sum_real + sum_rec + sum0 + sum_back

  end function ewald_pot

  !> Calculate real and reciprocal space sum cutoffs
  subroutine calculate_cutoff(q2sum,rcut,hcut,eta,lrmax,lhmax)
    use struct_basic
    use param

    real*8, intent(in) :: q2sum
    real*8, intent(out) :: rcut, hcut, eta
    integer, intent(out) :: lrmax(3), lhmax(3)

    real*8, parameter :: sgrow = 1.4d0
    real*8, parameter :: epscut = 1d-5
    real*8, parameter :: eeps = 1d-12

    integer :: i
    real*8 :: aux
    integer :: ia, ib, ic
    real*8 :: alrmax(3)
    real*8 :: rcut1, rcut2, err_real
    real*8 :: hcut1, hcut2, err_rec

    ! determine shortest vector in real space
    aux = 0d0
    do i = 1, 3
       if (cr%aa(i) > aux) then
          aux = cr%aa(i)
          ia = i
       end if
    end do
    ! determine shortest vector in reciprocal space, dif. from ia
    aux = 0d0
    do i = 1, 3
       if (cr%ar(i) > aux .and. i /= ia) then
          aux = cr%ar(i)
          ic = i
       end if
    end do
    ! the remaining vector is ib
    ib = 1
    do i = 1, 3
       if (i /= ia .and. i /= ic) ib = i
    end do

    ! convergence parameter
    eta = sqrt(cr%omega / pi / cr%aa(ib) / sin(cr%bb(ic)*rad))

    ! real space cutoff
    rcut1 = 1d0
    rcut2 = 2d0 / sgrow
    err_real = 1d30
    do while (err_real >= eeps)
       rcut2 = rcut2 * sgrow
       err_real = pi * cr%ncel**2 * q2sum / cr%omega * eta**2 * erfc(rcut2 / eta)
    end do
    do while (rcut2-rcut1 >= epscut)
       rcut = 0.5*(rcut1+rcut2)
       err_real = pi * cr%ncel**2 * q2sum / cr%omega * eta**2 * erfc(rcut / eta)
       if (err_real > eeps) then
          rcut1 = rcut
       else
          rcut2 = rcut
       endif
    end do
    rcut = 0.5*(rcut1+rcut2)
    ! real space cells to explore
    alrmax = 0d0
    alrmax(1) = cr%aa(2) * cr%aa(3) * sin(cr%bb(1)*rad)
    alrmax(2) = cr%aa(1) * cr%aa(3) * sin(cr%bb(2)*rad)
    alrmax(3) = cr%aa(1) * cr%aa(2) * sin(cr%bb(3)*rad)
    lrmax = ceiling(rcut * alrmax / cr%omega)

    ! reciprocal space cutoff
    hcut1 = 1d0
    hcut2 = 2d0 / sgrow
    err_rec = 1d30
    do while(err_rec >= eeps)
       hcut2 = hcut2 * sgrow
       err_rec = cr%ncel**2 * q2sum / sqpi / eta * erfc(eta * hcut2 / 2)
    end do
    do while(hcut2-hcut1 > epscut)
       hcut = 0.5*(hcut1+hcut2)
       err_rec = cr%ncel**2 * q2sum / sqpi / eta * erfc(eta * hcut / 2)
       if (err_rec > eeps) then
          hcut1 = hcut
       else
          hcut2 = hcut
       endif
    end do
    hcut = 0.5*(hcut1+hcut2)
    ! reciprocal space cells to explore
    lhmax = ceiling(cr%aa(ia) / tpi * hcut)

  end subroutine calculate_cutoff

  !> Sum the cell charges and squares
  subroutine sum_charges(qsum,q2sum)
    use struct_basic
    use tools_io

    real*8, intent(out) :: qsum, q2sum

    integer :: i

    ! calculate sum of charges and charges**2
    qsum = 0d0
    q2sum = 0d0
    do i = 1, cr%nneq
       if (cr%at(i)%qat == 0) &
          call ferror('ewald_energy','Some of the charges are 0',faterr)
       qsum = qsum + real(cr%at(i)%mult * cr%at(i)%qat,8)
       q2sum = q2sum + real(cr%at(i)%mult * cr%at(i)%qat**2,8)
    end do

  end subroutine sum_charges

end module ewald
