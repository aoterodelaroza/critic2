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

  ! OpenMP is used in the Ewald routines only above these sizes.
  integer, parameter :: ewald_omp_nmin = 48  !< min atoms to thread the real-space loops
  integer, parameter :: ewald_omp_kmin = 256 !< min k-vectors to thread the reciprocal loop

contains

  !> Calculate real and reciprocal space sum cutoffs. qfrac contains the
  !> atomic charges (one per species).
  module subroutine calculate_ewald_cutoffs(c,qfrac,rcut,hcut,eta,qsum,lrmax,lhmax)
    use tools_io, only: ferror, faterr
    use param, only: pi, rad, sqpi, tpi
    class(crystal), intent(inout) :: c

    real*8, intent(in) :: qfrac(:)
    real*8, intent(out) :: rcut, hcut, eta, qsum
    integer, intent(out) :: lrmax(3), lhmax(3)

    real*8, parameter :: sgrow = 1.4d0
    real*8, parameter :: epscut = 1d-5
    real*8, parameter :: eeps = 1d-12

    integer :: i
    real*8 :: aux, q2sum
    integer :: ia, ib, ic
    real*8 :: alrmax(3), qq
    real*8 :: rcut1, rcut2, err_real
    real*8 :: hcut1, hcut2, err_rec

    ! calculate sum of charges and charges**2
    qsum = 0d0
    q2sum = 0d0
    do i = 1, c%nneq
       qq = qfrac(c%at(i)%is)
       if (abs(qq) < 1d-6) &
          call ferror('calculate_ewald_cutoffs','Some of the charges are 0',faterr)
       qsum = qsum + real(c%at(i)%mult * qq,8)
       q2sum = q2sum + real(c%at(i)%mult * qq**2,8)
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
    eta = sqrt(c%omega / pi / c%aa(ib) / sin(c%bb(ic)*rad))

    ! real space cutoff
    rcut1 = 1d0
    rcut2 = 2d0 / sgrow
    err_real = 1d30
    do while (err_real >= eeps)
       rcut2 = rcut2 * sgrow
       err_real = pi * c%ncel**2 * q2sum / c%omega * eta**2 * erfc(rcut2 / eta)
    end do
    do while (rcut2-rcut1 >= epscut)
       rcut = 0.5*(rcut1+rcut2)
       err_real = pi * c%ncel**2 * q2sum / c%omega * eta**2 * erfc(rcut / eta)
       if (err_real > eeps) then
          rcut1 = rcut
       else
          rcut2 = rcut
       endif
    end do
    rcut = 0.5*(rcut1+rcut2)
    ! real space cells to explore
    alrmax = 0d0
    alrmax(1) = c%aa(2) * c%aa(3) * sin(c%bb(1)*rad)
    alrmax(2) = c%aa(1) * c%aa(3) * sin(c%bb(2)*rad)
    alrmax(3) = c%aa(1) * c%aa(2) * sin(c%bb(3)*rad)
    lrmax = ceiling(rcut * alrmax / c%omega)

    ! reciprocal space cutoff
    hcut1 = 1d0
    hcut2 = 2d0 / sgrow
    err_rec = 1d30
    do while(err_rec >= eeps)
       hcut2 = hcut2 * sgrow
       err_rec = c%ncel**2 * q2sum / sqpi / eta * erfc(eta * hcut2 / 2)
    end do
    do while(hcut2-hcut1 > epscut)
       hcut = 0.5*(hcut1+hcut2)
       err_rec = c%ncel**2 * q2sum / sqpi / eta * erfc(eta * hcut / 2)
       if (err_rec > eeps) then
          hcut1 = hcut
       else
          hcut2 = hcut
       endif
    end do
    hcut = 0.5*(hcut1+hcut2)
    ! reciprocal space cells to explore
    lhmax = ceiling(c%aa(ia) / tpi * hcut)

  end subroutine calculate_ewald_cutoffs

  !> Calculates the Ewald electrostatic energy, using the input charges
  !> qfrac (one per species).
  module function ewald_energy(c,qfrac) result(ewe)
    use tools_io, only: uout, string
    class(crystal), intent(inout) :: c
    real*8, intent(in) :: qfrac(:)
    real*8 :: ewe

    real*8 :: x(3), pot
    integer :: i
    real*8 :: rcut, hcut, eta, qsum
    integer :: lrmax(3), lhmax(3)
    real*8, allocatable :: qat(:)

    call c%calculate_ewald_cutoffs(qfrac,rcut,hcut,eta,qsum,lrmax,lhmax)

    ! per-atom charges for the generalized potential routine
    allocate(qat(c%ncel))
    do i = 1, c%ncel
       qat(i) = qfrac(c%atcel(i)%is)
    end do

    write (uout,'("+ Electrostatic potential at atomic positions")')
    write (uout,'("#id name mult    charge         Vel(Ha/e)")')
    ewe = 0d0
    do i = 1, c%nneq
       x = c%at(i)%x
       pot = c%ewald_pot(qat,x,rcut,hcut,eta,lrmax,lhmax)
       ewe = ewe + c%at(i)%mult * qfrac(c%at(i)%is) * pot
       write (uout,'(99(A," "))') string(i,4), string(c%at(i)%name,4),&
          string(c%at(i)%mult,4),&
          string(qfrac(c%at(i)%is),'e',14,6,3), string(pot,'e',18,10,3)
    end do
    write (uout,*)
    ewe = ewe / 2d0

  end function ewald_energy

  !> Ewald electrostatic potential at an arbitrary position x (crystallographic
  !> coords.) for the per-atom charges qat (length c%ncel). If x is the nucleus
  !> j, return pot - q_j/|r-rj| at rj (the singular self term is regularized by
  !> the h=0 self term). Atomic units.
  module function ewald_pot(c,qat,x,rcut,hcut,eta,lrmax,lhmax) result(pot)
    use param, only: pi, sqpi, icrd_crys
    class(crystal), intent(inout) :: c
    real*8, intent(in) :: qat(:)
    real*8, intent(in) :: x(3)
    real*8, intent(in) :: rcut, hcut, eta
    integer, intent(in) :: lrmax(3), lhmax(3)
    real*8 :: pot

    real*8 :: qnuc, g(c%ncel)
    integer :: idnuc

    ! is this a nuclear position? -> get its charge for the h=0 self term
    idnuc = c%identify_atom(x,icrd_crys)
    if (idnuc > 0) then
       qnuc = qat(idnuc)
    else
       qnuc = 0d0
    end if

    ! pot = sum_j q_j g_j(x) - self(nucleus) - background
    call ewald_gval(c,x,eta,rcut,hcut,lrmax,lhmax,g)
    pot = dot_product(qat,g) - 2d0*qnuc/(sqpi*eta) - sum(qat)*eta**2*pi/c%omega

  end function ewald_pot

  !> Build the periodic Ewald 1/r interaction matrix aew(n,n) (n =
  !> c%ncel, atomic units): aew(i,j) = potential at atom i due to a
  !> unit charge on atom j and all its periodic images, including the
  !> self term on the diagonal and the neutralizing-background
  !> term. The Ewald energy for a set of per-atom charges q is then
  !> 0.5*dot_product(q, matmul(aew,q)).
  module subroutine ewald_matrix(c,aew,rcut,hcut,eta,lrmax,lhmax)
    use param, only: pi, sqpi
    class(crystal), intent(inout) :: c
    real*8, intent(out) :: aew(:,:)
    real*8, intent(in) :: rcut, hcut, eta
    integer, intent(in) :: lrmax(3), lhmax(3)

    integer :: i
    real*8 :: g(c%ncel)

    ! each row is independent (ewald_gval is read-only in c); no reduction needed
    !$omp parallel do private(i,g) if(c%ncel >= ewald_omp_nmin)
    do i = 1, c%ncel
       call ewald_gval(c,c%atcel(i)%x,eta,rcut,hcut,lrmax,lhmax,g)
       aew(i,:) = g
    end do
    !$omp end parallel do
    ! self term (diagonal) and neutralizing background (uniform)
    do i = 1, c%ncel
       aew(i,i) = aew(i,i) - 2d0/(sqpi*eta)
    end do
    aew = aew - pi*eta*eta/c%omega

  end subroutine ewald_matrix

  !> Periodic Ewald 1/r energy (eene) and, accumulated into grad/vir,
  !> its forces and virial, for the fixed per-atom charges qat (length
  !> c%ncel).
  module subroutine ewald_energy_grad(c,qat,rcut,hcut,eta,lrmax,lhmax,eene,grad,vir,dostress,&
     nbptr,nbj,nblv)
    use param, only: pi, sqpi, tpi
    class(crystal), intent(inout) :: c
    real*8, intent(in) :: qat(:)
    real*8, intent(in) :: rcut, hcut, eta
    integer, intent(in) :: lrmax(3), lhmax(3)
    real*8, intent(out) :: eene
    real*8, intent(inout) :: grad(:,:)
    real*8, intent(inout) :: vir(3,3)
    logical, intent(in) :: dostress
    integer, intent(in), optional :: nbptr(:), nbj(:), nblv(:,:)

    integer :: n, i, j, i1, i2, i3, a, b, p
    real*8 :: cosk(c%ncel), sink(c%ncel) ! per-atom structure-factor phases (automatic)
    real*8 :: r, dh, bbarg
    real*8 :: px(3), kv(3), kc(3), d(3), g(3), gp, fmag, gr
    real*8 :: sc, ss, brec, s2, fac, qtot
    real*8 :: minvt(3,3)
    logical :: uselist, dorecip, dopar

    n = c%ncel
    eene = 0d0
    minvt = transpose(c%m_c2x)
    uselist = present(nbptr) .and. present(nbj) .and. present(nblv)
    dopar = n >= ewald_omp_nmin

    ! --- real space: 0.5 sum q_i q_j erfc(r/eta)/r within rcut. Self image
    ! (i==j, L=0) is excluded by r<1d-10; own images (i==j, L/=0) contribute
    ! energy and stress but give no net force. Uses the caller's precomputed
    ! CSR pair list when supplied (fast path for MD/relaxation); otherwise a
    ! generic image-cell triple loop. ---
    if (uselist) then
       !$omp parallel do private(i,p,j,px,d,r,gr,gp,fmag,g,a,b) reduction(+:eene,grad,vir) if(dopar)
       do i = 1, n
          do p = nbptr(i), nbptr(i+1)-1
             j = nbj(p)
             px = c%atcel(i)%x - c%atcel(j)%x - real(nblv(:,p),8)
             d = matmul(c%m_x2c,px)
             r = norm2(d)
             if (r < 1d-10 .or. r > rcut) cycle
             gr = erfc(r/eta)/r
             eene = eene + 0.5d0*qat(i)*qat(j)*gr
             gp = -erfc(r/eta)/(r*r) - (2d0/(eta*sqpi))*exp(-(r/eta)**2)/r
             fmag = 0.5d0*qat(i)*qat(j)*gp/r
             g = fmag*d
             if (i /= j) then
                grad(:,i) = grad(:,i) + g
                grad(:,j) = grad(:,j) - g
             end if
             if (dostress) then
                do a = 1, 3
                   do b = 1, 3
                      vir(a,b) = vir(a,b) - d(a)*g(b)
                   end do
                end do
             end if
          end do
       end do
       !$omp end parallel do
    else
       !$omp parallel do private(i,j,i1,i2,i3,px,d,r,gr,gp,fmag,g,a,b) reduction(+:eene,grad,vir) if(dopar)
       do i = 1, n
          do j = 1, n
             do i1 = -lrmax(1), lrmax(1)
                do i2 = -lrmax(2), lrmax(2)
                   do i3 = -lrmax(3), lrmax(3)
                      px = c%atcel(i)%x - c%atcel(j)%x - real((/i1,i2,i3/),8)
                      d = matmul(c%m_x2c,px)
                      r = norm2(d)
                      if (r < 1d-10 .or. r > rcut) cycle
                      gr = erfc(r/eta)/r
                      eene = eene + 0.5d0*qat(i)*qat(j)*gr
                      gp = -erfc(r/eta)/(r*r) - (2d0/(eta*sqpi))*exp(-(r/eta)**2)/r
                      fmag = 0.5d0*qat(i)*qat(j)*gp/r
                      g = fmag*d
                      if (i /= j) then
                         grad(:,i) = grad(:,i) + g
                         grad(:,j) = grad(:,j) - g
                      end if
                      if (dostress) then
                         do a = 1, 3
                            do b = 1, 3
                               vir(a,b) = vir(a,b) - d(a)*g(b)
                            end do
                         end do
                      end if
                   end do
                end do
             end do
          end do
       end do
       !$omp end parallel do
    end if

    ! --- self energy (diagonal Ewald term) ---
    do i = 1, n
       eene = eene - qat(i)*qat(i)/(sqpi*eta)
    end do
    ! --- neutralizing background (zero for a neutral cell) ---
    qtot = sum(qat(1:n))
    eene = eene - 0.5d0*pi*eta*eta/c%omega*qtot*qtot

    ! --- reciprocal space: sum_k brec |S(k)|^2 (dominant cost; parallel over k) ---
    dorecip = product(2*lhmax+1) >= ewald_omp_kmin
    !$omp parallel do collapse(3) private(kv,kc,dh,bbarg,brec,sc,ss,i,cosk,sink,s2,fac,a,b) &
    !$omp    reduction(+:eene,grad,vir) if(dorecip)
    do i1 = -lhmax(1), lhmax(1)
       do i2 = -lhmax(2), lhmax(2)
          do i3 = -lhmax(3), lhmax(3)
             kv = tpi*real((/i1,i2,i3/),8)
             kc = tpi*matmul(minvt,real((/i1,i2,i3/),8))  ! cartesian k-vector
             dh = norm2(kc)
             if (dh < 1d-12 .or. dh > hcut) cycle
             bbarg = 0.5d0*dh*eta
             brec = (tpi/c%omega)/dh**2*exp(-bbarg**2)
             sc = 0d0
             ss = 0d0
             do i = 1, n
                cosk(i) = cos(dot_product(kv,c%atcel(i)%x))
                sink(i) = sin(dot_product(kv,c%atcel(i)%x))
                sc = sc + qat(i)*cosk(i)
                ss = ss + qat(i)*sink(i)
             end do
             s2 = sc*sc + ss*ss
             eene = eene + brec*s2
             do i = 1, n
                grad(:,i) = grad(:,i) + brec*2d0*qat(i)*(ss*cosk(i) - sc*sink(i))*kc
             end do
             if (dostress) then
                fac = 2d0*(1d0/dh**2 + 0.25d0*eta*eta)
                do a = 1, 3
                   do b = 1, 3
                      vir(a,b) = vir(a,b) + brec*s2*(merge(1d0,0d0,a==b) - fac*kc(a)*kc(b))
                   end do
                end do
             end if
          end do
       end do
    end do
    !$omp end parallel do

  end subroutine ewald_energy_grad

  !> Per-source Ewald Green's function at field point x (crystallographic): for
  !> each cell atom j returns g(j) = sum over j's images of erfc(r/eta)/r plus
  !> the reciprocal-space sum, i.e. the (regularized, self-image-at-L=0 excluded)
  !> potential at x from a unit charge on atom j and its periodic images, WITHOUT
  !> the singular self and background terms (added by the caller).
  subroutine ewald_gval(c,x,eta,rcut,hcut,lrmax,lhmax,g)
    use param, only: tpi
    class(crystal), intent(in) :: c
    real*8, intent(in) :: x(3), eta, rcut, hcut
    integer, intent(in) :: lrmax(3), lhmax(3)
    real*8, intent(out) :: g(:)

    integer :: i, i1, i2, i3
    real*8 :: rcut2, px(3), d(3), d2, dd, kc(3), dh, bbarg, ktrm, arg

    g = 0d0
    rcut2 = rcut*rcut

    ! real space
    do i1 = -lrmax(1), lrmax(1)
       do i2 = -lrmax(2), lrmax(2)
          do i3 = -lrmax(3), lrmax(3)
             do i = 1, c%ncel
                px = x - c%atcel(i)%x - real((/i1,i2,i3/),8)
                d = matmul(c%m_x2c,px)
                d2 = dot_product(d,d)
                if (d2 < 1d-12 .or. d2 > rcut2) cycle
                dd = sqrt(d2)
                g(i) = g(i) + erfc(dd/eta)/dd
             end do
          end do
       end do
    end do

    ! reciprocal space
    do i1 = -lhmax(1), lhmax(1)
       do i2 = -lhmax(2), lhmax(2)
          do i3 = -lhmax(3), lhmax(3)
             kc = tpi*matmul(transpose(c%m_c2x),real((/i1,i2,i3/),8))
             dh = norm2(kc)
             if (dh < 1d-12 .or. dh > hcut) cycle
             bbarg = 0.5d0*dh*eta
             ktrm = (tpi/c%omega)*2d0/dh**2*exp(-bbarg**2)
             do i = 1, c%ncel
                arg = tpi*dot_product(real((/i1,i2,i3/),8), x - c%atcel(i)%x)
                g(i) = g(i) + ktrm*cos(arg)
             end do
          end do
       end do
    end do

  end subroutine ewald_gval

  !> Calculate the core or promolecular densities on a grid with n(:)
  !> points. If a fragment is given, then only the atoms in it
  !> contribute.  This routine is thread-safe.
  module subroutine promolecular_array3(c,f,n,zpsp,fr)
    use grid1mod, only: grid1
    use fragmentmod, only: fragment
    use param, only: icrd_crys
    class(crystal), intent(inout) :: c
       real*8, intent(inout), allocatable :: f(:,:,:)
    integer, intent(in) :: n(3)
    integer, intent(in), optional :: zpsp(:)
    type(fragment), intent(in), optional :: fr

    integer :: i, j, k
    real*8 :: x(3), xdelta(3,3), rdum1(3), rdum2(3,3), rho

    if (allocated(f)) deallocate(f)
    allocate(f(n(1),n(2),n(3)))

    do i = 1, 3
       xdelta(:,i) = 0d0
       xdelta(i,i) = 1d0 / real(n(i),8)
    end do

    !$omp parallel do private(x,rho,rdum1,rdum2)
    do k = 1, n(3)
       do j = 1, n(2)
          do i = 1, n(1)
             x = (i-1) * xdelta(:,1) + (j-1) * xdelta(:,2) + (k-1) * xdelta(:,3)
             call c%promolecular_atom(x,icrd_crys,rho,rdum1,rdum2,0,zpsp,fr)

             f(i,j,k) = rho
          end do
       end do
    end do
    !$omp end parallel do

  end subroutine promolecular_array3

  !> Calculate the packing ratio (in %) using the nearest-neighbor
  !> information. Each atom is assigned a ratio equal to half the distance
  !> to its nearest neighbor.
  module function get_pack_ratio(c) result(px)
    use param, only: pi
    class(crystal), intent(inout) :: c
    real*8 :: px

    integer :: i
    real*8 :: rnn2

    px = 0d0
    do i = 1, c%nneq
       rnn2 = c%get_rnn2(i)
       px = px + c%at(i)%mult * 4d0/3d0 * pi * rnn2**3
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
  pure module subroutine get_kpoints(c,rk,nk)
    class(crystal), intent(in) :: c
    real*8, intent(in) :: rk
    integer, intent(out) :: nk(3)

    nk = int(max(1d0,rk * c%ar + 0.5d0))

  end subroutine get_kpoints

end submodule complex
