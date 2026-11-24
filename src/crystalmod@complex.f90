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

  ! The ideal coordination shapes of coordination_shape: number of
  ! vertices and name (vertices in shape_vertices).
  integer, parameter :: nshape = 29
  integer, parameter :: shape_nv(nshape) = (/2,2,2,3,3,3,4,4,4,5,5,6,6,6,6,7,7,7,8,8,8,8,&
     9,9,10,12,12,12,14/)
  character(len=29), parameter :: shape_name(nshape) = (/character(len=29) ::&
     "linear","V-shaped","L-shaped","trigonal planar","trigonal pyramidal","T-shaped",&
     "tetrahedral","square planar","seesaw","trigonal bipyramidal","square pyramidal",&
     "octahedral","trigonal prismatic","hexagonal planar","pentagonal pyramidal",&
     "pentagonal bipyramidal","capped octahedral","capped trigonal prismatic","cubic",&
     "square antiprismatic","dodecahedral","hexagonal bipyramidal",&
     "tricapped trigonal prismatic","capped square antiprismatic",&
     "bicapped square antiprismatic","cuboctahedral","icosahedral","anticuboctahedral",&
     "rhombic dodecahedral"/)

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
    !$omp    reduction(+:eene,grad,vir) if(dorecip) schedule(dynamic)
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
    use global, only: abort_poll, progress_start, progress_step
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

    call progress_start(n(3)*n(2),"lines")
    !$omp parallel do private(x,rho,rdum1,rdum2) collapse(2)
    do k = 1, n(3)
       do j = 1, n(2)
          ! a cancelled calculation (the GUI): the caller stops
          if (abort_poll()) cycle
          do i = 1, n(1)
             x = (i-1) * xdelta(:,1) + (j-1) * xdelta(:,2) + (k-1) * xdelta(:,3)
             call c%promolecular_atom(x,icrd_crys,rho,rdum1,rdum2,0,zpsp,fr)

             f(i,j,k) = rho
          end do
          call progress_step()
       end do
    end do
    !$omp end parallel do

  end subroutine promolecular_array3

  !> Find the voids in a scalar field sampled over the unit cell,
  !> understood as the connected regions where the field is lower than
  !> isoval. f is the field on a uniform grid spanning the cell (from
  !> promolecular_array3, say); two of its points belong to the same
  !> void if they are neighbors along one of the three lattice
  !> directions. Returns the total void volume (vtot, bohr^3) and the
  !> number of voids found (nvoid). For each void, in order of
  !> decreasing volume: its volume (vol, bohr^3), the position of its
  !> deepest point (xdeep, crystallographic coordinates), and the
  !> value of the field at that point (rhodeep). errmsg is non-empty
  !> in case of error.
  module subroutine void_domains(c,f,isoval,vtot,nvoid,vol,xdeep,rhodeep,errmsg,ilbl)
    use grid3mod, only: grid_domains
    class(crystal), intent(inout) :: c
    real*8, intent(in) :: f(:,:,:)
    real*8, intent(in) :: isoval
    real*8, intent(out) :: vtot
    integer, intent(out) :: nvoid
    real*8, allocatable, intent(out) :: vol(:)
    real*8, allocatable, intent(out) :: xdeep(:,:)
    real*8, allocatable, intent(out) :: rhodeep(:)
    character(len=:), allocatable, intent(out) :: errmsg
    integer, allocatable, intent(inout), optional :: ilbl(:,:,:)

    integer :: i, j, k, id, n(3)
    real*8 :: dvol
    integer, allocatable :: lbl(:,:,:), ncount(:), imin(:,:)
    real*8, allocatable :: rmin(:)

    vtot = 0d0
    n = shape(f)

    ! the domains of the grid under the isovalue (wrapping around the cell)
    call grid_domains(f,isoval,.true.,nvoid,lbl,ncount,errmsg)
    if (len_trim(errmsg) > 0) then
       allocate(vol(0),xdeep(3,0),rhodeep(0))
       if (present(ilbl)) then
          if (allocated(ilbl)) deallocate(ilbl)
       end if
       return
    end if

    ! the deepest point of each void
    allocate(imin(3,nvoid),rmin(nvoid))
    rmin = huge(1d0)
    imin = 1
    do k = 1, n(3)
       do j = 1, n(2)
          do i = 1, n(1)
             id = lbl(i,j,k)
             if (id == 0) cycle
             if (f(i,j,k) < rmin(id)) then
                rmin(id) = f(i,j,k)
                imin(:,id) = (/i,j,k/)
             end if
          end do
       end do
    end do

    ! the volumes, largest first (grid_domains numbers them that way)
    dvol = c%omega / real(product(int(n,8)),8)
    vtot = real(sum(ncount),8) * dvol
    allocate(vol(nvoid),xdeep(3,nvoid),rhodeep(nvoid))
    do i = 1, nvoid
       vol(i) = real(ncount(i),8) * dvol
       xdeep(:,i) = real(imin(:,i)-1,8) / real(n,8)
       rhodeep(i) = rmin(i)
    end do

    ! return the labels if requested
    if (present(ilbl)) then
       if (allocated(ilbl)) deallocate(ilbl)
       call move_alloc(lbl,ilbl)
    end if

  end subroutine void_domains

  !> Calculate the coordination polyhedron centered on point x0
  !> (crystallographic coordinates). The vertices of the polyhedron are
  !> the atoms of species is0 (or, if is0 is zero, the atoms with atomic
  !> number iz0; or, if ispc is present, the atoms of any species flagged
  !> in it) at a distance between rmin and rmax (bohr) from the center;
  !> the atom at the center itself, if there is one, is not a vertex. Returns the number of vertices (nat), the shortest and
  !> longest vertex distance (dmin and dmax, bohr), the number of faces
  !> (nf), and the volume of the polyhedron (vol, bohr^3). If dev is
  !> present, it returns the maximum distance from a vertex to the best-fit
  !> plane of all of them (bohr), zero if they are coplanar and the
  !> polyhedron is therefore flat; dev is also zero if nat <= 3, when it is
  !> not calculated. ier is non-zero if the triangulation failed. Fewer than
  !> four vertices in range is not an error: it gives nat <= 3 and a zero
  !> volume, and no polyhedron.
  module subroutine coord_polyhedron(c,x0,is0,iz0,rmin,rmax,nat,dmin,dmax,nf,vol,ier,dev,ispc)
    use param, only: icrd_crys
    class(crystal), intent(inout) :: c
    real*8, intent(in) :: x0(3)
    integer, intent(in) :: is0
    integer, intent(in) :: iz0
    real*8, intent(in) :: rmin
    real*8, intent(in) :: rmax
    integer, intent(out) :: nat
    real*8, intent(out) :: dmin
    real*8, intent(out) :: dmax
    integer, intent(out) :: nf
    real*8, intent(out) :: vol
    integer, intent(out) :: ier
    real*8, intent(out), optional :: dev
    logical, intent(in), optional :: ispc(:)

    integer :: i, nat0
    integer, allocatable :: eid(:), lvec(:,:)
    real*8, allocatable :: dist(:), xstar(:,:), up2dsp(:,:)

    nat = 0
    dmin = 0d0
    dmax = 0d0

    ! the atoms that make up the vertices
    if (present(ispc)) then
       ! several species at once: a distance window per species, zero for the
       ! ones that are not vertices (the same way the GUI collects its corners)
       allocate(up2dsp(c%nspc,2))
       up2dsp = 0d0
       do i = 1, min(c%nspc,size(ispc,1))
          if (ispc(i)) up2dsp(i,2) = rmax
       end do
       call c%list_near_atoms(x0,icrd_crys,.true.,nat0,eid=eid,dist=dist,lvec=lvec,&
          up2dsp=up2dsp,nozero=.true.)
       deallocate(up2dsp)
    elseif (is0 /= 0) then
       call c%list_near_atoms(x0,icrd_crys,.true.,nat0,eid=eid,dist=dist,lvec=lvec,&
          up2d=rmax,ispc0=is0,nozero=.true.)
    else
       call c%list_near_atoms(x0,icrd_crys,.true.,nat0,eid=eid,dist=dist,lvec=lvec,&
          up2d=rmax,iz0=iz0,nozero=.true.)
    end if

    ! discard the atoms closer than rmin (the list comes sorted by distance)
    allocate(xstar(3,nat0))
    do i = 1, nat0
       if (dist(i) < rmin) cycle
       nat = nat + 1
       xstar(:,nat) = c%x2c(c%atcel(eid(i))%x + lvec(:,i))
       if (nat == 1) dmin = dist(i)
       dmax = dist(i)
    end do

    call polyhedron_from_points(c%x2c(x0),xstar(:,1:nat),nf,vol,ier,dev)

  end subroutine coord_polyhedron

  !> Convex polyhedron with vertices xv(:,1:n) around the center xc
  !> (Cartesian, any units). Returns the number of faces (nf) and the
  !> volume (vol, units cubed), as the sum of the tetrahedra joining
  !> the faces to the center. If dev is present, it returns the
  !> maximum distance from a vertex to the best-fit plane of all of
  !> them, zero if they are coplanar and the polyhedron is therefore
  !> flat. A solid needs four vertices: with three or fewer, nf, vol
  !> and dev are zero, and ier is zero too. ier is non-zero if the
  !> triangulation failed. A vertex on the center is ignored by the
  !> triangulation (but not by dev).
  module subroutine polyhedron_from_points(xc,xv,nf,vol,ier,dev)
    use tools_math, only: mixed, plane_from_points
    use iso_c_binding, only: c_int, c_double, c_ptr
    real*8, intent(in) :: xc(3)
    real*8, intent(in) :: xv(:,:)
    integer, intent(out) :: nf
    real*8, intent(out) :: vol
    integer, intent(out) :: ier
    real*8, intent(out), optional :: dev

    integer :: j, n
    integer, allocatable :: iface(:,:)
    real*8 :: xp1(3), xp2(3), xp3(3), xcen(3), xnor(3)
    type(c_ptr) :: ctx

    interface
       subroutine runqhull_basintriangulate_step1(n,x0,xvert,nf,ctx,ier) bind(c)
         import c_int, c_double, c_ptr
         integer(c_int), value :: n
         real(c_double) :: x0(3)
         real(c_double) :: xvert(3,n)
         integer(c_int) :: nf
         type(c_ptr) :: ctx
         integer(c_int) :: ier
       end subroutine runqhull_basintriangulate_step1
       subroutine runqhull_basintriangulate_step2(nf,iface,ctx) bind(c)
         import c_int, c_double, c_ptr
         integer(c_int), value :: nf
         integer(c_int) :: iface(3,nf)
         type(c_ptr), value :: ctx
       end subroutine runqhull_basintriangulate_step2
    end interface

    n = size(xv,2)
    nf = 0
    vol = 0d0
    ier = 0
    if (present(dev)) dev = 0d0

    ! three or fewer points are not passed to the triangulation, which
    ! cannot build an initial simplex out of them
    if (n <= 3) return

    ! how far the vertices stray from a common plane
    if (present(dev)) call plane_from_points(xv,n,xcen,xnor,dev)

    ! project on a sphere and triangulate the convex polyhedron. The context
    ! is allocated by step1 even if it fails, and freed by step2
    call runqhull_basintriangulate_step1(n,xc,xv,nf,ctx,ier)
    allocate(iface(3,max(nf,1)))
    call runqhull_basintriangulate_step2(nf,iface,ctx)
    if (ier /= 0) return

    ! the volume of the polyhedron, as the sum of the tetrahedra that join
    ! its faces to the center
    do j = 1, nf
       xp1 = xv(:,iface(1,j)) - xc
       xp2 = xv(:,iface(2,j)) - xc
       xp3 = xv(:,iface(3,j)) - xc
       vol = vol + abs(mixed(xp1,xp2,xp3)) / 6d0
    end do

  end subroutine polyhedron_from_points

  !> Closest ideal shape to the coordination environment with vertices
  !> xv(:,1:n) around the center xc (Cartesian, any units). The angles
  !> between the vertices at the center are sorted and compared with
  !> those of the ideal shapes with the same number of vertices (all
  !> at the same distance from the center). Returns the closest shape
  !> in ish(1) and the root mean square deviation of the angles in
  !> dev(1) (degrees), and the second closest in ish(2) and dev(2).
  !> ish is zero if there is no shape (no ideal shape with n vertices,
  !> or only one for ish(2), or a vertex on the center). The names of
  !> the shapes come from coordination_shape_name.
  module subroutine coordination_shape(xc,xv,ish,dev)
    use tools, only: qcksort
    real*8, intent(in) :: xc(3)
    real*8, intent(in) :: xv(:,:)
    integer, intent(out) :: ish(2)
    real*8, intent(out) :: dev(2)

    integer :: i, n
    real*8 :: d
    real*8, allocatable :: x(:,:), ang(:), angref(:)

    ish = 0
    dev = 0d0
    n = size(xv,2)
    if (.not.any(shape_nv == n)) return

    ! the directions from the center
    allocate(x(3,n),ang(n*(n-1)/2),angref(n*(n-1)/2))
    do i = 1, n
       x(:,i) = xv(:,i) - xc
       d = norm2(x(:,i))
       if (d < 1d-10) return
       x(:,i) = x(:,i) / d
    end do
    call sorted_angles(x,ang)

    ! compare with the ideal shapes
    dev = huge(1d0)
    do i = 1, nshape
       if (shape_nv(i) /= n) cycle
       call shape_vertices(i,x)
       call sorted_angles(x,angref)
       d = sqrt(sum((ang - angref)**2) / size(ang))
       if (d < dev(1)) then
          ish(2) = ish(1)
          dev(2) = dev(1)
          ish(1) = i
          dev(1) = d
       elseif (d < dev(2)) then
          ish(2) = i
          dev(2) = d
       end if
    end do
    where (ish == 0) dev = 0d0

  contains
    !> The sorted angles (degrees) between the unit vectors x(:,1:n).
    subroutine sorted_angles(x,ang)
      use param, only: rad
      real*8, intent(in) :: x(:,:)
      real*8, intent(out) :: ang(:)

      integer :: i, j, k, n

      n = size(x,2)
      k = 0
      do i = 1, n
         do j = i+1, n
            k = k + 1
            ang(k) = acos(max(min(dot_product(x(:,i),x(:,j)),1d0),-1d0)) / rad
         end do
      end do
      call qcksort(ang)

    end subroutine sorted_angles
  end subroutine coordination_shape

  !> Name of the ideal coordination shape with index ish (see
  !> coordination_shape). Empty if ish is out of range.
  module function coordination_shape_name(ish) result(name)
    integer, intent(in) :: ish
    character(len=:), allocatable :: name

    if (ish < 1 .or. ish > nshape) then
       name = ""
    else
       name = trim(shape_name(ish))
    end if

  end function coordination_shape_name

  !> Unit vectors from the center to the vertices of the ideal
  !> coordination shape ish, in x(:,1:shape_nv(ish)).
  subroutine shape_vertices(ish,x)
    use tools_io, only: ferror, faterr
    use param, only: rad
    integer, intent(in) :: ish
    real*8, intent(inout) :: x(:,:)

    ! dodecahedron, A and B sites (most favorable polyhedron of Hoard
    ! and Silverton, Inorg. Chem. 2 (1963) 235)
    real*8, parameter :: tha = 35.2d0
    real*8, parameter :: thb = 73.5d0
    ! square pyramid, apex-center-base angle: the center is above the
    ! base, as in most of them (90 degrees is a vacant octahedron)
    real*8, parameter :: spyang = 100d0
    integer :: nput
    real*8 :: ttet, tcube, tp, ta, tc, ti

    nput = 0
    ttet = acos(-1d0/3d0) / rad ! tetrahedral angle
    tcube = 0.5d0 * ttet ! cube, body diagonal to edge
    tp = atan(2d0/sqrt(3d0)) / rad ! trigonal prism with equal edges
    ta = asin(sqrt(4d0/(4d0+sqrt(2d0)))) / rad ! square antiprism with equal edges
    tc = acos(sqrt(2d0/3d0)) / rad ! cuboctahedra, out of the hexagon plane
    ti = atan(2d0) / rad ! icosahedron

    ! the polar angle (from z) and the azimuth of each vertex, in degrees
    select case(ish)
    case(1) ! linear
       call put(1,0d0,0d0)
       call put(2,180d0,0d0)
    case(2) ! V-shaped
       call put(1,0d0,0d0)
       call put(2,ttet,0d0)
    case(3) ! L-shaped
       call put(1,0d0,0d0)
       call put(2,90d0,0d0)
    case(4) ! trigonal planar
       call ring(1,3,90d0,0d0)
    case(5) ! trigonal pyramidal
       call ring(1,3,ttet,0d0)
    case(6) ! T-shaped
       call put(1,90d0,0d0)
       call put(2,90d0,180d0)
       call put(3,90d0,90d0)
    case(7) ! tetrahedral
       call put(1,0d0,0d0)
       call ring(2,3,ttet,0d0)
    case(8) ! square planar
       call ring(1,4,90d0,0d0)
    case(9) ! seesaw
       call put(1,0d0,0d0)
       call put(2,180d0,0d0)
       call put(3,90d0,0d0)
       call put(4,90d0,120d0)
    case(10) ! trigonal bipyramidal
       call bipyramid(1,3)
    case(11) ! square pyramidal
       call put(1,0d0,0d0)
       call ring(2,4,spyang,0d0)
    case(12) ! octahedral
       call bipyramid(1,4)
    case(13) ! trigonal prismatic
       call prism(1)
    case(14) ! hexagonal planar
       call ring(1,6,90d0,0d0)
    case(15) ! pentagonal pyramidal
       call put(1,0d0,0d0)
       call ring(2,5,90d0,0d0)
    case(16) ! pentagonal bipyramidal
       call bipyramid(1,5)
    case(17) ! capped octahedral, over the (1,1,1) face
       call bipyramid(1,4)
       call put(7,tcube,45d0)
    case(18) ! capped trigonal prismatic, over a square face
       call prism(1)
       call put(7,90d0,60d0)
    case(19) ! cubic
       call cube(1)
    case(20) ! square antiprismatic
       call antiprism(1)
    case(21) ! dodecahedral (D2d)
       call put(1,tha,0d0)
       call put(2,tha,180d0)
       call put(3,180d0-tha,90d0)
       call put(4,180d0-tha,270d0)
       call put(5,thb,90d0)
       call put(6,thb,270d0)
       call put(7,180d0-thb,0d0)
       call put(8,180d0-thb,180d0)
    case(22) ! hexagonal bipyramidal
       call bipyramid(1,6)
    case(23) ! tricapped trigonal prismatic
       call prism(1)
       call ring(7,3,90d0,60d0)
    case(24) ! capped square antiprismatic
       call antiprism(1)
       call put(9,0d0,0d0)
    case(25) ! bicapped square antiprismatic
       call antiprism(1)
       call put(9,0d0,0d0)
       call put(10,180d0,0d0)
    case(26) ! cuboctahedral (fcc)
       call ring(1,6,90d0,0d0)
       call ring(7,3,tc,30d0)
       call ring(10,3,180d0-tc,90d0)
    case(27) ! icosahedral
       call put(1,0d0,0d0)
       call put(2,180d0,0d0)
       call ring(3,5,ti,0d0)
       call ring(8,5,180d0-ti,36d0)
    case(28) ! anticuboctahedral (hcp)
       call ring(1,6,90d0,0d0)
       call ring(7,3,tc,30d0)
       call ring(10,3,180d0-tc,30d0)
    case(29) ! rhombic dodecahedral (bcc, first and second neighbors)
       call cube(1)
       call bipyramid(9,4)
    end select

    ! the tables of the shapes and this routine have to agree
    if (nput /= shape_nv(ish) .or. nput /= size(x,2)) &
       call ferror('shape_vertices','inconsistent number of vertices in an ideal shape',faterr)

  contains
    subroutine put(k,th,ph)
      integer, intent(in) :: k
      real*8, intent(in) :: th, ph
      x(:,k) = (/sin(th*rad)*cos(ph*rad),sin(th*rad)*sin(ph*rad),cos(th*rad)/)
      nput = nput + 1
    end subroutine put
    subroutine ring(k0,n,th,ph0)
      integer, intent(in) :: k0, n
      real*8, intent(in) :: th, ph0
      integer :: k
      do k = 0, n-1
         call put(k0+k,th,ph0+360d0*k/n)
      end do
    end subroutine ring
    subroutine bipyramid(k0,n)
      integer, intent(in) :: k0, n
      call put(k0,0d0,0d0)
      call put(k0+1,180d0,0d0)
      call ring(k0+2,n,90d0,0d0)
    end subroutine bipyramid
    subroutine cube(k0)
      integer, intent(in) :: k0
      call ring(k0,4,tcube,45d0)
      call ring(k0+4,4,180d0-tcube,45d0)
    end subroutine cube
    subroutine prism(k0)
      integer, intent(in) :: k0
      call ring(k0,3,tp,0d0)
      call ring(k0+3,3,180d0-tp,0d0)
    end subroutine prism
    subroutine antiprism(k0)
      integer, intent(in) :: k0
      call ring(k0,4,ta,0d0)
      call ring(k0+4,4,180d0-ta,45d0)
    end subroutine antiprism
  end subroutine shape_vertices

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
  !> sampling. relerr = use enough points to obtain a standard
  !> deviation divided by the volume equal to this value. If
  !> rtable(0:maxzat0) is present, use those radii instead of the
  !> van der walls radii. If ratom(1:c%ncel) is present, use those
  !> radii, one per atom in the complete list (takes precedence over
  !> rtable; for radii that cannot be written as a function of the
  !> atomic number, like half the nearest-neighbor distance).
  module function vdw_volume(c,relerr,rtable,ratom) result(vvdw)
    use global, only: abort_poll, progress_start, progress_set
    use param, only: VBIG, atmvdw, icrd_cart, maxzat0
    class(crystal), intent(inout) :: c
    real*8, intent(in) :: relerr
    real*8, intent(in), optional :: rtable(0:maxzat0)
    real*8, intent(in), optional :: ratom(1:c%ncel)
    real*8 :: vvdw

    real*8 :: xmin(3), xmax(3), x(3), vtot, svol, pp, rmax, nneed
    integer :: i, nat

    integer, parameter :: npollstep = 1024 ! samples between cancel checks
    real*8, allocatable :: rvdw(:,:)
    integer*8 :: nin, ntot
    logical :: again, peratom

    ! build the list of atomic radii
    peratom = present(ratom)
    if (peratom) then
       rmax = maxval(ratom)
    else
       allocate(rvdw(c%nspc,2))
       do i = 1, c%nspc
          rvdw(i,1) = 0d0
          if (present(rtable)) then
             rvdw(i,2) = rtable(c%spc(i)%z)
          else
             rvdw(i,2) = atmvdw(c%spc(i)%z)
          end if
       end do
       rmax = maxval(rvdw(:,2))
    end if

    ! calculate the encompassing box
    if (c%ismolecule) then
       xmin = VBIG
       xmax = -VBIG
       do i = 1, c%ncel
          xmin = min(xmin,c%atcel(i)%r)
          xmax = max(xmax,c%atcel(i)%r)
       end do
       xmin = xmin - rmax
       xmax = xmax + rmax
       vtot = product(xmax-xmin)
    else
       vtot = c%omega
    end if

    ! use Monte-Carlo to determine the volume
    call progress_start(0,"samples")
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
       if (peratom) then
          call c%list_near_atoms(x,icrd_cart,.false.,nat,up2dcidx=ratom)
       else
          call c%list_near_atoms(x,icrd_cart,.false.,nat,up2dsp=rvdw)
       end if
       if (nat > 0) then
          nin = nin + 1
       end if

       pp = real(nin,8) / real(ntot,8)
       svol = vtot * sqrt(pp * (1-pp)/real(ntot,8))
       vvdw = vtot * pp
       if (ntot > 100) then
          again = (svol > relerr * vvdw)
       end if

       ! a cancelled calculation (the GUI): the caller discards the
       ! volume. The progress is the samples done of those the current
       ! fraction pp needs to reach relerr (an estimate)
       if (mod(ntot,int(npollstep,8)) == 0) then
          if (abort_poll()) exit
          nneed = min(max(vdw_volume_npoints(pp,relerr),real(ntot,8)),2d9)
          call progress_set(nint(min(real(ntot,8),nneed)),nint(nneed))
       end if
    end do

  end function vdw_volume

  !> Number of Monte Carlo samples vdw_volume needs to reach relative
  !> error relerr when a fraction p of the samples falls inside the
  !> spheres: it stops when the standard deviation of the volume
  !> divided by the volume falls below relerr, after (1-p)/(p*relerr^2)
  !> samples, and not before 100.
  pure module function vdw_volume_npoints(p,relerr) result(n)
    real*8, intent(in) :: p, relerr
    real*8 :: n

    if (p <= 0d0) then
       n = 100d0
    else
       n = max((1d0 - p) / (p * relerr**2),100d0)
    end if

  end function vdw_volume_npoints

  !> Calculate the number of k-points (nk) for a given rk-length. Uses
  !> the VASP formula.
  pure module subroutine get_kpoints(c,rk,nk)
    class(crystal), intent(in) :: c
    real*8, intent(in) :: rk
    integer, intent(out) :: nk(3)

    nk = int(max(1d0,rk * c%ar + 0.5d0))

  end subroutine get_kpoints

  !> Ranges of the rklength (bohr) that give the same k-point grid
  !> (get_kpoints), scanning rk = kpoints_rkstep, 2*kpoints_rkstep, ...,
  !> rkmax. The i-th range goes from rk1(i) to rk2(i) with grid nk(:,i).
  pure module subroutine get_kpoints_ranges(c,rkmax,rk1,rk2,nk)
    class(crystal), intent(in) :: c
    real*8, intent(in) :: rkmax
    real*8, allocatable, intent(out) :: rk1(:), rk2(:)
    integer, allocatable, intent(out) :: nk(:,:)

    integer, parameter :: maxstep = 100000000 ! longest scan (avoids integer overflow)

    integer :: j, n, nstep, ipass, nkj(3), nkold(3)

    ! the scan counts steps, so rounding does not accumulate along it;
    ! count the ranges in the first pass, fill them in the second
    nstep = nint(max(min(rkmax / kpoints_rkstep,real(maxstep,8)),1d0))
    do ipass = 1, 2
       call c%get_kpoints(kpoints_rkstep,nkold)
       n = 1
       if (ipass == 2) then
          rk1(1) = kpoints_rkstep
          nk(:,1) = nkold
       end if
       do j = 2, nstep
          call c%get_kpoints(j*kpoints_rkstep,nkj)
          if (all(nkj == nkold)) cycle
          nkold = nkj
          if (ipass == 2) rk2(n) = (j-1)*kpoints_rkstep
          n = n + 1
          if (ipass == 2) then
             rk1(n) = j*kpoints_rkstep
             nk(:,n) = nkj
          end if
       end do
       if (ipass == 1) allocate(rk1(n),rk2(n),nk(3,n))
    end do
    rk2(n) = nstep*kpoints_rkstep

  end subroutine get_kpoints_ranges

  !> Rattle the atoms of crystal c: append nstruct new seeds to seed(:), each
  !> a copy of c with every atom displaced by mag (bohr) in a random direction.
  !> nseed is the number of seeds already in the list, and is updated.
  module subroutine bulk_rattle_seeds(c,nstruct,mag,seed,nseed)
    use crystalseedmod, only: crystalseed, realloc_crystalseed
    use tools_math, only: random_unit_vector
    type(crystal), intent(in) :: c
    integer, intent(in) :: nstruct
    real*8, intent(in) :: mag
    type(crystalseed), allocatable, intent(inout) :: seed(:)
    integer, intent(inout) :: nseed

    type(crystalseed) :: seed0
    real*8 :: xdelta(3)
    integer :: i, j

    if (nstruct < 1) return
    if (.not.allocated(seed)) then
       allocate(seed(nseed + nstruct))
    elseif (nseed + nstruct > size(seed,1)) then
       call realloc_crystalseed(seed,2*(nseed + nstruct))
    end if

    ! the template seed, from the structure as it is now
    call c%makeseed(seed0,.false.)

    do i = nseed+1, nseed+nstruct
       seed(i) = seed0
       do j = 1, c%ncel
          call random_unit_vector(xdelta)
          seed(i)%x(:,j) = c%c2x(c%atcel(j)%r + mag * xdelta)
       end do
    end do
    nseed = nseed + nstruct

  end subroutine bulk_rattle_seeds

end submodule complex
