! Copyright (c) 2007-2022 Alberto Otero de la Roza <aoterodelaroza@gmail.com>,
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

! Routines for editing molecular and crystal structures
submodule (crystalmod) complex
  implicit none

contains

  !> Calculate real and reciprocal space sum cutoffs
  module subroutine calculate_ewald_cutoffs(c)
    use tools_io, only: ferror, faterr
    use param, only: pi, rad, sqpi, tpi
    class(crystal), intent(inout) :: c

    real*8, parameter :: sgrow = 1.4d0
    real*8, parameter :: epscut = 1d-5
    real*8, parameter :: eeps = 1d-12

    integer :: i
    real*8 :: aux
    integer :: ia, ib, ic
    real*8 :: alrmax(3), qq
    real*8 :: rcut1, rcut2, err_real
    real*8 :: hcut1, hcut2, err_rec

    if (c%isewald) return

    ! calculate sum of charges and charges**2
    c%qsum = 0d0
    c%q2sum = 0d0
    do i = 1, c%nneq
       qq = c%spc(c%at(i)%is)%qat
       if (abs(qq) < 1d-6) &
          call ferror('calculate_ewald_cutoffs','Some of the charges are 0',faterr)
       c%qsum = c%qsum + real(c%at(i)%mult * qq,8)
       c%q2sum = c%q2sum + real(c%at(i)%mult * qq**2,8)
    end do

    ! determine shortest vector in real space
    aux = 0d0
    do i = 1, 3
       if (c%aa(i) > aux) then
          aux = c%aa(i)
          ia = i
       end if
    end do
    ! determine shortest vector in reciprocal space, dif. from ia
    aux = 0d0
    do i = 1, 3
       if (c%ar(i) > aux .and. i /= ia) then
          aux = c%ar(i)
          ic = i
       end if
    end do
    ! the remaining vector is ib
    ib = 1
    do i = 1, 3
       if (i /= ia .and. i /= ic) ib = i
    end do

    ! convergence parameter
    c%eta = sqrt(c%omega / pi / c%aa(ib) / sin(c%bb(ic)*rad))

    ! real space cutoff
    rcut1 = 1d0
    rcut2 = 2d0 / sgrow
    err_real = 1d30
    do while (err_real >= eeps)
       rcut2 = rcut2 * sgrow
       err_real = pi * c%ncel**2 * c%q2sum / c%omega * c%eta**2 * erfc(rcut2 / c%eta)
    end do
    do while (rcut2-rcut1 >= epscut)
       c%rcut = 0.5*(rcut1+rcut2)
       err_real = pi * c%ncel**2 * c%q2sum / c%omega * c%eta**2 * erfc(c%rcut / c%eta)
       if (err_real > eeps) then
          rcut1 = c%rcut
       else
          rcut2 = c%rcut
       endif
    end do
    c%rcut = 0.5*(rcut1+rcut2)
    ! real space cells to explore
    alrmax = 0d0
    alrmax(1) = c%aa(2) * c%aa(3) * sin(c%bb(1)*rad)
    alrmax(2) = c%aa(1) * c%aa(3) * sin(c%bb(2)*rad)
    alrmax(3) = c%aa(1) * c%aa(2) * sin(c%bb(3)*rad)
    c%lrmax = ceiling(c%rcut * alrmax / c%omega)

    ! reciprocal space cutoff
    hcut1 = 1d0
    hcut2 = 2d0 / sgrow
    err_rec = 1d30
    do while(err_rec >= eeps)
       hcut2 = hcut2 * sgrow
       err_rec = c%ncel**2 * c%q2sum / sqpi / c%eta * erfc(c%eta * hcut2 / 2)
    end do
    do while(hcut2-hcut1 > epscut)
       c%hcut = 0.5*(hcut1+hcut2)
       err_rec = c%ncel**2 * c%q2sum / sqpi / c%eta * erfc(c%eta * c%hcut / 2)
       if (err_rec > eeps) then
          hcut1 = c%hcut
       else
          hcut2 = c%hcut
       endif
    end do
    c%hcut = 0.5*(hcut1+hcut2)
    ! reciprocal space cells to explore
    c%lhmax = ceiling(c%aa(ia) / tpi * c%hcut)
    c%isewald = .true.

  end subroutine calculate_ewald_cutoffs

  !> Calculates the Ewald electrostatic energy, using the input charges.
  module function ewald_energy(c) result(ewe)
    use tools_io, only: uout, string
    class(crystal), intent(inout) :: c
    real*8 :: ewe

    real*8 :: x(3), pot
    integer :: i

    if (.not.c%isewald) &
       call c%calculate_ewald_cutoffs()

    write (uout,'("+ Electrostatic potential at atomic positions")')
    write (uout,'("#id name mult    charge         Vel(Ha/e)")')
    ewe = 0d0
    do i = 1, c%nneq
       x = c%at(i)%x
       pot = c%ewald_pot(x)
       ewe = ewe + c%at(i)%mult * c%spc(c%at(i)%is)%qat * pot
       write (uout,'(99(A," "))') string(i,4), string(c%spc(c%at(i)%is)%name,4),&
          string(c%at(i)%mult,4),&
          string(c%spc(c%at(i)%is)%qat,'e',14,6,3), string(pot,'e',18,10,3)
    end do
    write (uout,*)
    ewe = ewe / 2d0

  end function ewald_energy

  !> Calculate the Ewald electrostatic potential at an arbitrary
  !> position x (crystallographic coords.)  If x is the nucleus j,
  !> return pot - q_j / |r-rj| at rj.
  module function ewald_pot(c,x)
    use param, only: tpi, pi, sqpi, icrd_crys
    class(crystal), intent(inout) :: c
    real*8, intent(in) :: x(3)
    real*8 :: ewald_pot

    real*8 :: rcut2, qnuc
    integer :: i, i1, i2, i3, idnuc
    real*8 :: px(3), lvec(3), d2, d, dh
    real*8 :: sfac_c, sfacp, bbarg
    real*8 :: sum_real, sum_rec, sum0, sum_back

    if (.not.c%isewald) then
       !$omp critical (fill_ewald)
       call c%calculate_ewald_cutoffs()
       !$omp end critical (fill_ewald)
    end if

    ! is this a nuclear position? -> get charge
    idnuc = c%identify_atom_env(x,icrd_crys)
    if (idnuc > 0) then
       qnuc = c%spc(c%atcel(idnuc)%is)%qat
    else
       qnuc = 0d0
    endif

    ! real space sum
    rcut2 = c%rcut * c%rcut
    sum_real = 0
    do i1 = -c%lrmax(1),c%lrmax(1)
       do i2 = -c%lrmax(2),c%lrmax(2)
          do i3 = -c%lrmax(3),c%lrmax(3)
             lvec = real((/i1,i2,i3/),8)
             do i = 1,c%ncel
                px = x - c%atcel(i)%x - lvec
                d2 = dot_product(px,matmul(c%gtensor,px))
                if (d2 < 1d-12 .or. d2 > rcut2) cycle
                d = sqrt(d2) / c%eta
                sum_real = sum_real + c%spc(c%atcel(i)%is)%qat * erfc(d) / d
             end do
          end do
       end do
    end do
    sum_real = sum_real / c%eta

    ! reciprocal space sum
    sum_rec = 0
    do i1 = -c%lhmax(1),c%lhmax(1)
       do i2 = -c%lhmax(2),c%lhmax(2)
          do i3 = -c%lhmax(3),c%lhmax(3)
             lvec = tpi * (/i1,i2,i3/)
             dh = sqrt(dot_product(lvec,matmul(c%grtensor,lvec)))
             if (dh < 1d-12 .or. dh > c%hcut) cycle
             bbarg = 0.5d0 * dh * c%eta

             sfac_c = 0
             do i = 1, c%ncel
                sfac_c = sfac_c + c%spc(c%atcel(i)%is)%qat * &
                   cos(dot_product(lvec,x-c%atcel(i)%x))
             end do
             sfacp = 2d0 * sfac_c

             sum_rec = sum_rec + sfacp / dh**2 * exp(-bbarg**2)
          end do
       end do
    end do
    sum_rec = sum_rec * 2d0 * pi / c%omega

    ! h = 0 term, applied only at the nucleus
    sum0 = - 2d0 * qnuc / sqpi / c%eta

    ! compensating background charge term
    sum_back = -c%qsum * c%eta**2 * pi / c%omega

    ! sum up and exit
    ewald_pot = sum_real + sum_rec + sum0 + sum_back

  end function ewald_pot

  !> Calculate the core or promolecular densities on a grid with n(:)
  !> points. If a fragment is given, then only the atoms in it
  !> contribute.  This routine is thread-safe.
  module subroutine promolecular_grid(c,f,n,zpsp,fr)
    use grid3mod, only: grid3
    use grid1mod, only: grid1
    use fragmentmod, only: fragment
    use param, only: icrd_crys
    class(crystal), intent(inout) :: c
    type(grid3), intent(out) :: f
    integer, intent(in) :: n(3)
    integer, intent(in), optional :: zpsp(:)
    type(fragment), intent(in), optional :: fr

    integer :: i, j, k
    real*8 :: x(3), xdelta(3,3), rdum1(3), rdum2(3,3), rho

    ! initialize
    call f%end()
    f%isinit = .true.
    f%n = n
    allocate(f%f(n(1),n(2),n(3)))

    do i = 1, 3
       xdelta(:,i) = 0d0
       xdelta(i,i) = 1d0 / real(n(i),8)
    end do

    !$omp parallel do private(x,rho,rdum1,rdum2)
    do k = 1, n(3)
       do j = 1, n(2)
          do i = 1, n(1)
             x = (i-1) * xdelta(:,1) + (j-1) * xdelta(:,2) + (k-1) * xdelta(:,3)
             call c%promolecular_env(x,icrd_crys,rho,rdum1,rdum2,0,zpsp,fr)

             !$omp critical(write)
             f%f(i,j,k) = rho
             !$omp end critical(write)
          end do
       end do
    end do
    !$omp end parallel do

  end subroutine promolecular_grid

  !> Calculate the packing ratio (in %) using the nearest-neighbor
  !> information. Each atom is assigned a ratio equal to half the distance
  !> to its nearest neighbor.
  module function get_pack_ratio(c) result(px)
    use param, only: pi
    class(crystal), intent(inout) :: c
    real*8 :: px

    integer :: i

    px = 0d0
    do i = 1, c%nneq
       px = px + c%at(i)%mult * 4d0/3d0 * pi * c%at(i)%rnn2**3
    end do
    px = px / c%omega * 100d0

  end function get_pack_ratio

  !> Calculate the vdw volume in a molecule or crystal by Monte-Carlo
  !> sampling.  relerr = use enough points to obtain a standard
  !> deviation divided by the volume equal to this value. If
  !> rtable(1:maxzat0) is present, use those radii instead of the
  !> van der walls radii
  module function vdw_volume(c,relerr,rtable) result(vvdw)
    use param, only: VBIG, atmvdw, icrd_cart
    class(crystal), intent(inout) :: c
    real*8, intent(in) :: relerr
    real*8, intent(in), optional :: rtable(:)
    real*8 :: vvdw

    real*8 :: xmin(3), xmax(3), x(3), vtot, svol, pp
    integer :: i, nat
    real*8, allocatable :: rvdw(:,:)
    integer*8 :: nin, ntot
    logical :: again

    ! build the list of atomic radii
    allocate(rvdw(c%nspc,2))
    do i = 1, c%nspc
       rvdw(i,1) = 0d0
       if (present(rtable)) then
          rvdw(i,2) = rtable(c%spc(i)%z)
       else
          rvdw(i,2) = atmvdw(c%spc(i)%z)
       end if
    end do

    ! calculate the encompassing box
    if (c%ismolecule) then
       xmin = VBIG
       xmax = -VBIG
       do i = 1, c%ncel
          xmin = min(xmin,c%atcel(i)%r)
          xmax = max(xmax,c%atcel(i)%r)
       end do
       xmin = xmin - maxval(rvdw(:,2))
       xmax = xmax + maxval(rvdw(:,2))
       vtot = product(xmax-xmin)
    else
       vtot = c%omega
    end if

    ! use Monte-Carlo to determine the volume
    again = .true.
    ntot = 0
    nin = 0
    do while (again)
       ntot = ntot + 1
       call random_number(x)

       if (c%ismolecule) then
          x = xmin + x * (xmax - xmin)
       else
          x = c%x2c(x)
       end if
       call c%list_near_atoms(x,icrd_cart,.false.,nat,up2dsp=rvdw)
       if (nat > 0) then
          nin = nin + 1
       end if

       pp = real(nin,8) / real(ntot,8)
       svol = vtot * sqrt(pp * (1-pp)/real(ntot,8))
       vvdw = vtot * pp
       if (ntot > 100) then
          again = (svol > relerr * vvdw)
       end if
    end do

  end function vdw_volume

  !> Calculate the number of k-points (nk) for a given rk-length. Uses
  !> the VASP formula.
  module subroutine get_kpoints(c,rk,nk)
    class(crystal), intent(in) :: c
    real*8, intent(in) :: rk
    integer, intent(out) :: nk(3)

    nk = int(max(1d0,rk * c%ar + 0.5d0))

  end subroutine get_kpoints

end submodule complex
