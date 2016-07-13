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

    real*8 :: x(3)
    integer :: i

    if (.not.cr%ewald_ready) &
       call cr%calculate_ewald_cutoffs()
    
    ewe = 0d0
    do i = 1, cr%nneq
       x = cr%at(i)%x
       ewe = ewe + cr%at(i)%mult * cr%at(i)%qat * &
          ewald_pot(x,.true.)
    end do
    ewe = ewe / 2d0

  end subroutine ewald_energy

  !> calculate the Ewald electrostatic potential at 
  !> an arbitrary position x (crystallographic coords.)
  !> If x is the nucleus j, return pot - q_j / |r-rj| at rj.
  function ewald_pot(x,isnuc)
    use struct_basic
    use param

    real*8, intent(in) :: x(3)
    logical, intent(in) :: isnuc
    real*8 :: ewald_pot

    real*8 :: nuc_cutoff2 = 1d-14

    real*8 :: rcut2
    integer :: i, i1, i2, i3, qnuc
    real*8 :: px(3), lvec(3), d2, d, dh
    real*8 :: sfac_c, sfacp, bbarg
    real*8 :: sum_real, sum_rec, sum0, sum_back

    !$omp critical (fill_ewald)
    if (.not.cr%ewald_ready) &
       call cr%calculate_ewald_cutoffs()
    !$omp end critical (fill_ewald)
    
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
    rcut2 = cr%rcut * cr%rcut
    sum_real = 0
    do i1 = -cr%lrmax(1),cr%lrmax(1)
       do i2 = -cr%lrmax(2),cr%lrmax(2)
          do i3 = -cr%lrmax(3),cr%lrmax(3)
             lvec = real((/i1,i2,i3/),8)
             do i = 1,cr%ncel
                px = x - cr%atcel(i)%x - lvec
                d2 = dot_product(px,matmul(cr%gtensor,px))
                if (d2 < 1d-12 .or. d2 > rcut2) cycle
                d = sqrt(d2) / cr%eta
                sum_real = sum_real + cr%at(cr%atcel(i)%idx)%qat * erfc(d) / d
             end do
          end do
       end do
    end do
    sum_real = sum_real / cr%eta

    ! reciprocal space sum
    sum_rec = 0
    do i1 = -cr%lhmax(1),cr%lhmax(1)
       do i2 = -cr%lhmax(2),cr%lhmax(2)
          do i3 = -cr%lhmax(3),cr%lhmax(3)
             lvec = tpi * (/i1,i2,i3/)
             dh = sqrt(dot_product(lvec,matmul(cr%grtensor,lvec)))
             if (dh < 1d-12 .or. dh > cr%hcut) cycle
             bbarg = 0.5d0 * dh * cr%eta

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
       sum0 = - 2d0 * qnuc / sqpi / cr%eta
    else
       sum0 = 0d0
    end if

    ! compensating background charge term
    sum_back = -cr%qsum * cr%eta**2 * pi / cr%omega 

    ! sum up and exit
    ewald_pot = sum_real + sum_rec + sum0 + sum_back

  end function ewald_pot

end module ewald
