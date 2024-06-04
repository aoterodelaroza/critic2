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

! Routines for powder diffraction and related techniques
submodule (crystalmod) powderproc
  implicit none

  !xx! private procedures
  ! function loglike_bayes(t)

contains

  !> Calculate the powder diffraction pattern of the crystal structure.
  !> Modes of operation(mode) =
  !>   0 : return a simulated diffractogram using Gaussian functions
  !>       at every peak.
  !>   1 : return the list of 2*theta vs. intensity values.
  !> Common input:
  !>   th2ini0 = initial 2*theta (in degrees).
  !>   th2end0 = final 2*theta (in degrees).
  !>   lambda0 = wavelength of the radiation (in angstrom).
  !>   fpol = polarization correction factor (0 = unpolarized, 0.95 = syncrhotron).
  !> Input, mode 0 (simulated diffractogram):
  !>   npts = number of points in the diffractoram.
  !>   sigma =  parameter for the Gaussian broadening.
  !>   ishard = if false represent the tails of the peaks outside the plot range.
  !> Output, global:
  !>   th2p = 2*theta for the located maxima (in degrees, optional)
  !>   ip = itensities at the peak maxima (optional)
  !>   hvecp= reciprocal lattice vector corresponding to the peaks (optional)
  !> Output, mode 0 (simulated diffractogram):
  !>   t = 2*theta uniform grid with npts points between th2ini0 and th2end0 (degrees, optional)
  !>   ih = intensity on the 2*theta grid (optional)
  module subroutine powder(c,mode,& ! mode
     th2ini0,th2end0,lambda0,fpol,& ! global input
     npts,sigma,ishard,& ! mode 0 input
     th2p,ip,hvecp,discardp,& ! global output
     t,ih) ! mode 0 output
    use param, only: pi, bohrtoa, cscatt, c2scatt
    use tools_io, only: ferror, faterr
    use tools, only: qcksort
    use types, only: realloc
    class(crystal), intent(in) :: c
    integer, intent(in) :: mode
    real*8, intent(in) :: th2ini0, th2end0
    real*8, intent(in) :: lambda0
    real*8, intent(in) :: fpol
    integer, intent(in), optional :: npts
    real*8, intent(in), optional :: sigma
    logical, intent(in), optional :: ishard
    real*8, allocatable, intent(inout), optional :: th2p(:)
    real*8, allocatable, intent(inout), optional :: ip(:)
    integer, allocatable, intent(inout), optional :: hvecp(:,:)
    real*8, intent(in), optional :: discardp
    real*8, allocatable, intent(inout), optional :: t(:)
    real*8, allocatable, intent(inout), optional :: ih(:)

    integer :: i, j, np, hcell, h, k, l, iz, idx
    real*8 :: th2ini, th2end, lambda, hvec(3), kvec(3), th, sth, th2, cth, cth2
    real*8 :: sigma2, smax, dh2, dh, dh3, sthlam, cterm, sterm
    real*8 :: ffac, as(4), bs(4), cs, c2s(4), int, mcorr, afac
    real*8 :: ipmax, ihmax, tshift
    integer :: hmax
    logical :: again
    integer, allocatable :: io(:)
    real*8, allocatable :: isum(:)

    integer, parameter :: mp = 20
    real*8, parameter :: ieps = 1d-5
    real*8, parameter :: theps = 1d-5
    real*8, parameter :: iepscont = 1d-10
    real*8, parameter :: intmax = 1d15

    ! check input consistency
    if (mode == 0 .and..not.(present(npts).and.present(sigma).and.present(ishard))) &
       call ferror('struct_powder','mode=0 requires npts, sigma, and ishard',faterr)
    if (mode == 1 .and.(present(npts).or.present(sigma).or.present(ishard))) &
       call ferror('struct_powder','mode=1 incompatible with npts, sigma, and ishard',faterr)
    if (mode == 1 .and.(present(t).or.present(ih))) &
       call ferror('struct_powder','mode=1 incompatible with t and ih',faterr)
    if (mode == 1 .and..not.(present(th2p).and.present(ip))) &
       call ferror('struct_powder','mode=1 requires th2p and ip as output',faterr)

    ! calculate tshift
    tshift = 0d0
    if (present(ishard)) then
       if (.not.ishard) then
          tshift = sigma * sqrt(abs(-2d0 * log(iepscont/intmax)))
          tshift = tshift * pi / 180d0
       end if
    end if

    ! prepare the grid limits
    if (present(t)) then
       if (allocated(t)) deallocate(t)
       allocate(t(npts))
       do i = 1, npts
          t(i) = th2ini0 + real(i-1,8) / real(npts-1,8) * (th2end0-th2ini0)
       end do
    end if
    if (present(ih)) then
       if (allocated(ih)) deallocate(ih)
       allocate(ih(npts),isum(npts))
       ih = 0d0
    end if
    th2ini = th2ini0 * pi / 180d0
    th2end = th2end0 * pi / 180d0

    ! allocate for peak list
    if (present(th2p)) then
       if (allocated(th2p)) deallocate(th2p)
       allocate(th2p(mp))
    end if
    if (present(ip)) then
       if (allocated(ip)) deallocate(ip)
       allocate(ip(mp))
    end if
    if (present(hvecp)) then
       if (allocated(hvecp)) deallocate(hvecp)
       allocate(hvecp(3,mp))
    end if

    ! cell limits, convert lambda to bohr
    lambda = lambda0 / bohrtoa
    smax = sin((th2end+tshift)/2d0)
    hmax = 2*ceiling(2*smax/lambda/minval(c%ar))
    ! broadening -> gaussian
    if (present(sigma)) sigma2 = sigma * sigma

    ! calculate the intensities
    np = 0
    hcell = 0
    again = .true.
    do while (again)
       hcell = hcell + 1
       again = (hcell <= hmax)
       do h = -hcell, hcell
          do k = -hcell, hcell
             do l = -hcell, hcell
                if (abs(h)/=hcell.and.abs(k)/=hcell.and.abs(l)/=hcell) cycle

                ! reciprocal lattice vector length
                hvec = real((/h,k,l/),8)
                dh2 = dot_product(hvec,matmul(c%grtensor,hvec))
                dh = sqrt(dh2)
                dh3 = dh2 * dh

                ! the theta is not outside the spectrum range
                sth = 0.5d0 * lambda * dh
                if (abs(sth) > smax) cycle
                th = asin(sth)
                th2 = 2d0 * th
                if (th2 < th2ini-tshift .or. th2 > th2end+tshift) cycle

                ! more stuff we need
                sthlam = dh / bohrtoa / 2d0
                kvec = 2 * pi * hvec

                ! calculate the raw intensity for this (hkl)
                cterm = 0d0
                sterm = 0d0
                do i = 1, c%ncel
                   iz = c%spc(c%atcel(i)%is)%z
                   if (iz < 1 .or. iz > size(cscatt,2)) &
                      call ferror('struct_powder','invalid Z -> no atomic scattering factors',faterr)
                   as = (/cscatt(1,iz),cscatt(3,iz),cscatt(5,iz),cscatt(7,iz)/)
                   bs = (/cscatt(2,iz),cscatt(4,iz),cscatt(6,iz),cscatt(8,iz)/)
                   cs = cscatt(9,iz)
                   if (dh < 2d0) then
                      ffac = as(1)*exp(-bs(1)*dh2)+as(2)*exp(-bs(2)*dh2)+&
                         as(3)*exp(-bs(3)*dh2)+as(4)*exp(-bs(4)*dh2)+cs
                   elseif (iz == 1) then
                      ffac = 0d0
                   else
                      c2s = c2scatt(:,iz)
                      ffac = exp(c2s(1)+c2s(2)*dh+c2s(3)*dh2/10d0+c2s(4)*dh3/100d0)
                   end if
                   ffac = ffac * exp(-sthlam**2)
                   cterm = cterm + ffac * cos(dot_product(kvec,c%atcel(i)%x))
                   sterm = sterm + ffac * sin(dot_product(kvec,c%atcel(i)%x))
                end do
                int = cterm**2 + sterm**2

                ! lorentz-polarization correction.
                ! this formula is compatible with Pecharsky, IuCR
                ! website, and the Fox, gdis, and dioptas implementations
                cth = cos(th)
                cth2 = cth*cth-sth*sth
                afac = (1-fpol) / (1+fpol)
                mcorr = (1 + afac * cth2 * cth2) / (1 + afac) / cth / (sth*sth)
                int = int * mcorr

                ! sum the peak
                if (int > ieps) then
                   again = .true.
                   ! use a gaussian profile, add the intensity
                   if (present(ih)) then
                      isum = int * exp(-(t-th2*180/pi)**2 / 2d0 / sigma2)
                      ih = ih + isum
                   end if

                   ! identify the new peak
                   if (present(th2p).and.present(ip)) then
                      if (th2 > th2ini .and. th2 < th2end) then
                         idx = 0
                         do i = 1, np
                            if (abs(th2p(i)-th2) <= theps) then
                               idx = i
                               exit
                            end if
                         end do

                         if (idx == 0) then
                            np = np + 1
                            if (np > size(th2p)) then
                               call realloc(th2p,2*np)
                               call realloc(ip,2*np)
                               if (present(hvecp)) call realloc(hvecp,3,2*np)
                            end if
                            th2p(np) = th2
                            ip(np) = int
                            if (present(hvecp)) hvecp(:,np) = (/h,k,l/)
                         else
                            ip(idx) = ip(idx) + int
                            if (present(hvecp)) hvecp(:,idx) = (/h,k,l/)
                         endif
                      end if
                   end if

                end if
             end do
          end do
       end do
    end do
    if (present(th2p)) call realloc(th2p,np)
    if (present(ip)) call realloc(ip,np)
    if (present(hvecp)) call realloc(hvecp,3,np)

    ! normalize the intensities to 100
    if (present(ip)) then
       ipmax = maxval(ip)
       if (ipmax > 1d-10) ip = ip / ipmax * 100
    end if
    if (present(ih)) then
       ihmax = maxval(ih)
       if (ihmax > 1d-10) ih = ih / ihmax * 100
    end if

    ! discard if necessary
    if (present(discardp) .and. present(ip)) then
       j = 0
       do i = 1, np
          if (ip(i) < discardp) cycle
          j = j + 1
          ip(j) = ip(i)
          if (present(th2p)) th2p(j) = th2p(i)
          if (present(hvecp)) hvecp(:,j) = hvecp(:,i)
       end do
       np = j
       call realloc(ip,np)
       if (present(th2p)) call realloc(th2p,np)
       if (present(hvecp)) call realloc(hvecp,3,np)
    end if

    ! sort the peaks
    if (present(th2p).and.present(ip)) then
       if (np == 0) &
          call ferror('struct_powder','no peaks found in the 2theta range',faterr)

       allocate(io(np))
       do i = 1, np
          io(i) = i
       end do
       call qcksort(th2p,io,1,np)
       th2p = th2p(io)
       ip = ip(io)
       if (present(hvecp)) hvecp = hvecp(:,io)
       deallocate(io)
    end if

  end subroutine powder

  !> Calculate the radial distribution function.  On input, npts is
  !> the number of bins from the initial (0) to the final (rend)
  !> distance and sigma is the Gaussian broadening of the peaks. If
  !> ishard=.false., represent the tails of the peaks outside the plot
  !> range.  On output, t is the distance grid, and ih is the value of
  !> the RDF. This routine is based on:
  !>   Willighagen et al., Acta Cryst. B 61 (2005) 29.
  !> except using the sqrt of the atomic numbers instead of the
  !> charges.
  !>
  !> Optionally, if npairs0/ipairs0 are given, return the RDF of only
  !> the pairs of species given in the ipairs0 array. npairs0 is the
  !> number of selected pairs.
  !>
  !> If ihat is present, return the RDF for non-equivalent atom i in
  !> ihat(:,i).
  !>
  !> If intpeak is present, use intpeak(i) for the intensity of a peak
  !> associated with non-equivalent atom i, instead of its atomic
  !> number.
  module subroutine rdf(c,rini,rend,sigma,ishard,npts,t,ih,npairs0,ipairs0,ihat,intpeak)
    use global, only: atomeps
    use param, only: icrd_cart
    class(crystal), intent(inout) :: c
    real*8, intent(in) :: rini
    real*8, intent(in) :: rend
    real*8, intent(in) :: sigma
    logical, intent(in) :: ishard
    integer, intent(in) :: npts
    real*8, allocatable, intent(inout) :: t(:)
    real*8, allocatable, intent(inout) :: ih(:)
    integer, intent(in), optional :: npairs0
    integer, intent(in), optional :: ipairs0(:,:)
    real*8, allocatable, intent(inout), optional :: ihat(:,:)
    real*8, intent(in), optional :: intpeak(:)

    integer :: i, j, k, nat, iz, jz, kz, npairs, iaux
    integer :: idx, mmult, mza
    real*8 :: int, sigma2, fac, tshift, imax
    logical :: found
    integer, allocatable :: eid(:), ipairs(:,:)
    real*8, allocatable :: dist(:), ihaux(:)

    real*8, parameter :: ieps = 1d-10

    ! set up the chosen pairs (ipairs(1) < ipairs(2))
    npairs = 0
    if (present(npairs0) .and. present(ipairs0)) then
       npairs = npairs0
       if (npairs > 0) then
          allocate(ipairs(size(ipairs0,1),npairs))
          ipairs = ipairs0(:,1:npairs)
          do i = 1, npairs
             if (ipairs(2,i) < ipairs(1,i)) then
                iaux = ipairs(1,i)
                ipairs(1,i) = ipairs(2,i)
                ipairs(2,i) = iaux
             end if
          end do
       end if
    end if

    ! set up the space for atomic RDFs
    if (present(ihat)) then
       if (allocated(ihat)) deallocate(ihat)
       allocate(ihat(npts,c%nneq))
       ihat = 0d0
    end if

    ! sigma2 and tshift for soft RDFs
    sigma2 = sigma * sigma
    if (.not.ishard) then
       mmult = 0
       mza = 0
       do i = 1, c%nneq
          mmult = max(c%at(i)%mult,mmult)
          mza = max(c%spc(c%at(i)%is)%z,mza)
       end do
       imax = real(mmult,8) * real(mza,8)
       tshift = sigma * sqrt(abs(-2d0 * log(ieps/imax)))
    else
       tshift = 0d0
    end if

    ! prepare the grid limits
    if (allocated(t)) deallocate(t)
    if (allocated(ih)) deallocate(ih)
    allocate(t(npts),ih(npts),ihaux(npts))
    do i = 1, npts
       t(i) = rini + real(i-1,8) / real(npts-1,8) * (rend-rini)
    end do
    ih = 0d0
    ihaux = 0d0

    ! calculate the radial distribution function for the crystal
    ! RDF(r) = sum_i=1...c%nneq sum_j=1...c%env%n sqrt(Zi*Zj) / c%nneq / rij * delta(r-rij)
    do i = 1, c%nneq
       if (npairs > 0) then
          found = .false.
          do k = 1, npairs
             if (ipairs(1,k) == c%at(i)%is) then
                found = .true.
                exit
             end if
          end do
          if (.not.found) cycle
       end if
       iz = c%spc(c%at(i)%is)%z
       call c%list_near_atoms(c%at(i)%r,icrd_cart,.true.,nat,eid,dist,up2d=rend+tshift,nozero=.true.)

       do j = 1, nat
          ! skip if at a distance lower than rini or higher than rend
          if (dist(j) < rini-tshift .or. dist(j) > rend+tshift) cycle

          ! get info for this atom from the environment
          kz = c%atcel(eid(j))%is
          jz = c%spc(kz)%z
          idx = c%atcel(eid(j))%idx

          ! only the chosen atomic pairs generate a peak
          if (npairs > 0) then
             found = .false.
             do k = 1, npairs
                if (ipairs(2,k) == kz) then
                   found = .true.
                   exit
                end if
             end do
             if (.not.found) cycle
          end if

          ! factor of 1/2 if this pair will be repeated twice
          if (i == idx) then
             fac = 0.5d0
          else
             fac = 1d0
          end if

          if (present(intpeak)) then
             int = fac * sqrt(intpeak(i) * intpeak(idx)) * c%at(i)%mult
          else
             int = fac * sqrt(real(iz * jz,8)) * c%at(i)%mult
          end if
          ihaux = int * exp(-(t - dist(j))**2 / 2d0 / sigma2)
          ih = ih + ihaux
          if (present(ihat)) &
             ihat(:,i) = ihat(:,i) + ihaux
       end do
    end do
    if (.not.c%ismolecule) then
       do i = 1, npts
          if (abs(t(i)) < atomeps) cycle
          ih(i) = ih(i) / t(i)**2
          if (present(ihat)) &
             ihat(i,:) = ihat(i,:) / t(i)**2
       end do
       ih = ih / c%ncel
       if (present(ihat)) &
          ihat = ihat / c%ncel
    end if

  end subroutine rdf

  !> Calcualte the average minimum distances up to imax nearest
  !> neighbors.  Return the AMD vector (1:imax) in res.
  !> Widdowson et al., "Average Minimum Distances of Periodic Point Sets-Foundational Invariants for Mapping Periodic Crystals"
  !> Match. Commun. Math. Comput. Chem., 87 (2022) 529, doi:10.46793/match.87-3.529W
  !> Reference implementation: https://github.com/dwiddo/average-minimum-distance
  module subroutine amd(c,imax,res)
    use tools_io, only: string
    use param, only: icrd_crys
    class(crystal), intent(inout) :: c
    integer, intent(in) :: imax
    real*8, intent(out) :: res(imax)

    integer :: i, jini
    integer :: nat, imax_
    real*8, allocatable :: dist(:)
    logical :: ok

    ! initialize, check input and allocate
    imax_ = imax
    if (c%ismolecule .and. imax > c%ncel-1) imax_ = c%ncel - 1
    res = 0d0

    ! calculate the amd
    jini = 1
    do i = 1, c%nneq
       ok = .false.
       call c%list_near_atoms(c%at(i)%x,icrd_crys,.true.,nat,dist=dist,up2n=imax_,nozero=.true.)
       res = res + c%at(i)%mult * dist
    end do
    res = res / real(c%ncel,8)

  end subroutine amd

  !> Terminate an xrpd_peaklist object.
  module subroutine xrpd_peaklist_end(p)
    class(xrpd_peaklist), intent(inout) :: p

    p%npeak = 0
    p%haveth2limits = .false.
    p%haveradiation = .false.
    p%havehvec = .false.
    p%havegradients = .false.
    p%havepeakshape = .false.
    p%lambda = 0d0
    p%fpol = 0d0
    p%th2ini = 0d0
    p%th2end = 0d0
    if (allocated(p%th2)) deallocate(p%th2)
    if (allocated(p%ip)) deallocate(p%ip)
    if (allocated(p%hvec)) deallocate(p%hvec)
    if (allocated(p%th2g)) deallocate(p%th2g)
    if (allocated(p%ipg)) deallocate(p%ipg)
    if (allocated(p%fwhm)) deallocate(p%fwhm)
    if (allocated(p%cgau)) deallocate(p%cgau)

  end subroutine xrpd_peaklist_end

  !> Return the list of powder peaks between the initial (th2ini0) and
  !> final (th2end0) 2*theta (in angles), with wavelength lambda0 (in
  !> angstrom) and polarization factor fpol. If usehvecp, calculate
  !> only the reflections in p%hvec (requires p%npeaks and p%hvec). If
  !> calcderivs, calculate derivatives and set p%th2g and p%ipg. If gg
  !> is present, use gg instead of the metric tensor from c. Return
  !> non-zero errmsg on error.
  module subroutine xrpd_peaks_from_crystal_powder(p,c,th2ini0,th2end0,lambda0,fpol,usehvecp,calcderivs,&
     errmsg,gg)
    use tools_math, only: matinv
    use tools, only: qcksort
    use types, only: realloc
    use param, only: pi, bohrtoa
    class(xrpd_peaklist), intent(inout) :: p
    type(crystal), intent(in) :: c
    real*8, intent(in) :: th2ini0, th2end0
    real*8, intent(in) :: lambda0
    real*8, intent(in) :: fpol
    logical, intent(in) :: usehvecp
    logical, intent(in) :: calcderivs
    character(len=:), allocatable, intent(out) :: errmsg
    real*8, intent(in), optional :: gg(3,3)

    integer :: i
    real*8 :: lambda, hvec(3)
    real*8 :: int, intg(6), th2, th2g(6), th2ini, th2end, smax, imax
    integer :: hmax, hcell, h, k, l
    logical :: again
    integer, allocatable :: io(:)
    real*8 :: grtensor(3,3), ar(3)
    integer, allocatable :: hvecini(:,:)
    integer :: nhvecini

    integer, parameter :: mp = 200
    real*8, parameter :: ieps = 1d-5

    ! consistency checks
    errmsg = ""
    if (usehvecp) then
       if (.not.p%havehvec.or..not.allocated(p%hvec).or.p%npeak<=0) then
          errmsg = 'requested usehvecp without hvec information'
          return
       end if
       nhvecini = p%npeak
       hvecini = p%hvec
    end if

    ! initialize angles and radiation
    call p%end()
    p%haveth2limits = .true.
    p%th2ini = th2ini0
    p%th2end = th2end0
    th2ini = th2ini0 * pi / 180d0
    th2end = th2end0 * pi / 180d0
    p%haveradiation = .true.
    p%lambda = lambda
    p%fpol = fpol

    ! grtensor and ar
    if (present(gg)) then
       grtensor = gg
       call matinv(grtensor,3)
       do i = 1, 3
          ar(i) = sqrt(grtensor(i,i))
       end do
    else
       grtensor = c%grtensor
       ar = c%ar
    end if

    ! allocate peak list and initialize output
    p%havehvec = .true.
    p%havegradients = calcderivs
    if (usehvecp) then
       p%npeak = nhvecini
       allocate(p%th2(p%npeak),p%ip(p%npeak))
       if (calcderivs) allocate(p%th2g(6,p%npeak),p%ipg(6,p%npeak))
    else
       if (allocated(p%hvec)) deallocate(p%hvec)
       allocate(p%th2(mp),p%ip(mp),p%hvec(3,mp))
       if (calcderivs) allocate(p%th2g(6,mp),p%ipg(6,mp))
    end if

    ! metric tensor, cell limits, convert lambda to bohr
    lambda = lambda0 / bohrtoa
    smax = sin(th2end/2d0)
    hmax = 2*ceiling(2*smax/lambda/minval(ar))

    ! calculate the intensities
    if (.not.usehvecp) then
       p%npeak = 0
       hcell = 0
       again = .true.
       do while (again)
          hcell = hcell + 1
          again = (hcell <= hmax)
          do h = -hcell, hcell
             do k = -hcell, hcell
                do l = -hcell, hcell
                   if (abs(h)/=hcell.and.abs(k)/=hcell.and.abs(l)/=hcell) cycle

                   ! calculate this reciprocal lattice vector
                   hvec = real((/h,k,l/),8)
                   call run_function_body()
                   if (len_trim(errmsg) > 0) return

                   ! sum the peak
                   if (int > ieps) then
                      again = .true.
                      p%npeak = p%npeak + 1
                      if (p%npeak > size(p%th2,1)) then
                         call realloc(p%th2,2*p%npeak)
                         call realloc(p%ip,2*p%npeak)
                         call realloc(p%hvec,3,2*p%npeak)
                         if (calcderivs) then
                            call realloc(p%th2g,6,2*p%npeak)
                            call realloc(p%ipg,6,2*p%npeak)
                         end if
                      end if
                      p%th2(p%npeak) = th2
                      p%ip(p%npeak) = int
                      p%hvec(1,p%npeak) = h
                      p%hvec(2,p%npeak) = k
                      p%hvec(3,p%npeak) = l
                      if (calcderivs) then
                         p%th2g(:,p%npeak) = th2g
                         p%ipg(:,p%npeak) = intg
                      end if
                   end if
                end do
             end do
          end do
       end do
       call realloc(p%th2,p%npeak)
       call realloc(p%ip,p%npeak)
       call realloc(p%hvec,3,p%npeak)
       if (calcderivs) then
          call realloc(p%th2g,6,p%npeak)
          call realloc(p%ipg,6,p%npeak)
       end if
    else
       p%hvec = hvecini
       do i = 1, p%npeak
          hvec = real(p%hvec(:,i),8)
          call run_function_body()
          if (len_trim(errmsg) > 0) return
          p%th2(i) = th2
          p%ip(i) = int
          if (calcderivs) then
             p%th2g(:,i) = th2g
             p%ipg(:,i) = intg
          end if
       end do
    end if

    ! sort the peaks
    allocate(io(p%npeak))
    do i = 1, p%npeak
       io(i) = i
    end do
    call qcksort(p%th2,io,1,p%npeak)
    p%th2 = p%th2(io)
    p%ip = p%ip(io)
    p%hvec = p%hvec(:,io)
    if (calcderivs) then
       p%th2g = p%th2g(:,io)
       p%ipg = p%ipg(:,io)
    end if
    deallocate(io)

    ! scale the angles
    p%th2 = p%th2 * 180d0 / pi
    if (calcderivs) p%th2g = p%th2g * 180d0 / pi

    ! normalize the intensities (highest peak is 100)
    imax = maxval(p%ip)
    p%ip = p%ip / imax * 100d0
    if (calcderivs) p%ipg = p%ipg / imax

  contains
    subroutine run_function_body()
      use param, only: cscatt, c2scatt
      real*8 :: xfac, ffacg(6), esthlam, afac, mcorr, mcorrg(6)
      real*8 :: cterm, sterm, ctermg(6), stermg(6), ffac
      real*8 :: kx, ckx, skx
      real*8 :: as(4), bs(4), cs, c2s(4)
      real*8 :: dhv(3), dhm(6), sthlam, kvec(3), ebs(4)
      real*8 :: th, dh, dh2, dh3, sth, cth, cth2
      integer :: i, iz

      ! initialize
      int = 0d0
      intg = 0d0
      th = 0d0
      th2 = 0d0

      ! plane distance and derivatives
      dh2 = dot_product(hvec,matmul(grtensor,hvec))
      dh = sqrt(dh2)
      dh3 = dh2 * dh
      if (calcderivs) then
         dhv = matmul(grtensor,hvec)
         dhm(1) = dhv(1) * dhv(1)
         dhm(2) = dhv(1) * dhv(2)
         dhm(3) = dhv(1) * dhv(3)
         dhm(4) = dhv(2) * dhv(2)
         dhm(5) = dhv(2) * dhv(3)
         dhm(6) = dhv(3) * dhv(3)
         dhm = -0.5d0 * dhm / dh
      end if

      ! the theta is not outside the spectrum range
      sth = 0.5d0 * lambda * dh
      if (abs(sth) > 1d0) return
      if (.not.usehvecp .and. abs(sth) > smax) return
      th = asin(sth)
      th2 = 2d0 * th
      if (.not.usehvecp .and. (th2 < th2ini .or. th2 > th2end)) return
      cth = cos(th)
      if (calcderivs) &
         th2g = 2d0 * sth / (cth * dh) * dhm

      ! more stuff we need
      sthlam = dh / bohrtoa / 2d0
      kvec = 2 * pi * hvec

      ! calculate the raw intensity for this (hkl)
      cterm = 0d0
      sterm = 0d0
      if (calcderivs) then
         ctermg = 0d0
         stermg = 0d0
      end if
      do i = 1, c%ncel
         iz = c%spc(c%atcel(i)%is)%z
         if (iz < 1 .or. iz > size(cscatt,2)) then
            errmsg = 'invalid Z -> no atomic scattering factors'
            return
         end if
         as = (/cscatt(1,iz),cscatt(3,iz),cscatt(5,iz),cscatt(7,iz)/)
         bs = (/cscatt(2,iz),cscatt(4,iz),cscatt(6,iz),cscatt(8,iz)/)
         cs = cscatt(9,iz)
         if (dh < 2d0) then
            ebs = exp(-bs * dh2)
            ffac = as(1) * ebs(1) + as(2) * ebs(2) + as(3) * ebs(3) + as(4) * ebs(4) + cs
            if (calcderivs) then
               xfac = as(1)*bs(1)*ebs(1) + as(2)*bs(2)*ebs(2) + as(3)*bs(3)*ebs(3) + as(4)*bs(4)*ebs(4)
               ffacg = -2d0 * dh * xfac * dhm
            end if
         elseif (iz == 1) then
            ffac = 0d0
            if (calcderivs) ffacg = 0d0
         else
            c2s = c2scatt(:,iz)
            ffac = exp(c2s(1)+c2s(2)*dh+c2s(3)*dh2/10d0+c2s(4)*dh3/100d0)
            if (calcderivs) then
               xfac = c2s(2) + 2d0 * c2s(3) * dh / 10d0 + 3d0 * c2s(4) * dh2 / 100d0
               ffacg = ffac * xfac * dhm
            end if
         end if
         esthlam = exp(-sthlam**2)
         if (calcderivs) &
            ffacg = (ffacg - sthlam / bohrtoa * ffac * dhm) * esthlam
         ffac = ffac * esthlam

         kx = dot_product(kvec,c%atcel(i)%x)
         ckx = cos(kx)
         skx = sin(kx)
         cterm = cterm + ffac * ckx
         sterm = sterm + ffac * skx
         if (calcderivs) then
            ctermg = ctermg + ffacg * ckx
            stermg = stermg + ffacg * skx
         end if
      end do
      int = cterm*cterm + sterm*sterm
      if (calcderivs) &
         intg = 2d0 * (cterm * ctermg + sterm * stermg)

      ! lorentz-polarization correction.
      ! this formula is compatible with Pecharsky, IuCR
      ! website, and the Fox, gdis, and dioptas implementations
      cth2 = cth*cth-sth*sth
      afac = (1-fpol) / (1+fpol)
      mcorr = (1 + afac * cth2 * cth2) / (1 + afac) / cth / (sth*sth)
      if (calcderivs) then
         mcorrg = 0.5d0 * (- 8 * afac * cth2 / ((1+afac) * sth) - mcorr * cth / sth - mcorr * cth2 / (sth * cth)) * th2g
         intg = intg * mcorr + int * mcorrg
      end if
      int = int * mcorr

    end subroutine run_function_body

  end subroutine xrpd_peaks_from_crystal_powder

  !> Read the XRPD peaks from a peaks file. The peaks are sorted by th2
  !> on output. If error, return non-zero-length errmsg.
  module subroutine xrpd_peaks_from_peaks_file(p,file,errmsg)
    use tools, only: qcksort
    use tools_io, only: fopen_read, getline, isreal, fclose, getword, lower, isreal
    use types, only: realloc
    class(xrpd_peaklist), intent(inout) :: p
    character*(*), intent(in) :: file
    character(len=:), allocatable, intent(out) :: errmsg

    integer :: i, lu, lp
    character(len=:), allocatable :: line, word
    real*8 :: r1, r2
    logical :: ok
    logical :: readat(4)

    integer, parameter :: mpeak = 1000
    integer, allocatable :: io(:)

    ! initialize output
    errmsg = ""
    call p%end()

    ! read the pattern from the peaks file
    readat = .false.
    p%havepeakshape = .true.
    allocate(p%th2(mpeak),p%ip(mpeak),p%fwhm(mpeak),p%cgau(mpeak))
    lu = fopen_read(file)
    do while (getline(lu,line))
       if (line(1:1) == "@") then
          ! interpret special keywords
          lp = 1
          word = getword(line,lp)
          word = lower(word)
          if (word == "@th2ini") then
             ok = isreal(p%th2ini,line,lp)
             readat(1) = .true.
          elseif (word == "@th2end") then
             ok = isreal(p%th2end,line,lp)
             readat(2) = .true.
          elseif (word == "@lambda") then
             ok = isreal(p%lambda,line,lp)
             readat(3) = .true.
          elseif (word == "@fpol") then
             ok = isreal(p%fpol,line,lp)
             readat(4) = .true.
          end if
          if (.not.ok) then
             errmsg = "Error reading keyword symbol: " // word
             return
          end if
          cycle
       end if

       lp = 1
       ok = isreal(r1,line,lp)
       ok = ok .and. isreal(r2,line,lp)
       if (.not.ok) then
          errmsg = "invalid input in peaks file"
          return
       end if
       p%npeak = p%npeak + 1
       if (p%npeak > size(p%th2,1)) then
          call realloc(p%th2,2*p%npeak)
          call realloc(p%ip,2*p%npeak)
          if (p%havepeakshape) then
             call realloc(p%fwhm,2*p%npeak)
             call realloc(p%cgau,2*p%npeak)
          end if
       end if
       p%th2(p%npeak) = r1
       p%ip(p%npeak) = r2

       if (p%havepeakshape) then
          ok = isreal(r1,line,lp)
          ok = ok .and. isreal(r2,line,lp)
          if (.not.ok) then
             p%havepeakshape = .false.
             if (allocated(p%fwhm)) deallocate(p%fwhm)
             if (allocated(p%cgau)) deallocate(p%cgau)
          else
             p%fwhm(p%npeak) = r1
             p%cgau(p%npeak) = r2
          end if
       end if
    end do
    call fclose(lu)
    call realloc(p%th2,p%npeak)
    call realloc(p%ip,p%npeak)
    if (allocated(p%fwhm)) call realloc(p%fwhm,p%npeak)
    if (allocated(p%cgau)) call realloc(p%cgau,p%npeak)

    ! process at symbols
    p%haveth2limits = (readat(1) .and. readat(2))
    p%haveradiation = (readat(3) .and. readat(4))

    ! sort the peaks
    allocate(io(p%npeak))
    do i = 1, p%npeak
       io(i) = i
    end do
    call qcksort(p%th2,io,1,p%npeak)
    p%th2 = p%th2(io)
    p%ip = p%ip(io)
    if (allocated(p%fwhm)) p%fwhm = p%fwhm(io)
    if (allocated(p%cgau)) p%cgau = p%cgau(io)
    deallocate(io)

  end subroutine xrpd_peaks_from_peaks_file

  !> Fit XRPD profile data from an xy file (xyfile) and return the
  !> peak list in p and the RMS of the fit in rms. The profile model
  !> uses a sum of pseudo-Voigt functions. If an error occurs, return
  !> the a non-zero errmsg with the error description. Optinal arguments:
  !> - verbose0: if present and true, write progress messages to stdout.
  !> - ymax_detect0: a candidate peak is added to the initial model if
  !>   the maximum is higher than ymax_detect0. If -huge(1d0), use the
  !>   default value (median of the profile intensities).
  !> - nadj0: a candidate peak is added if the maximum is surrounded by
  !>   nadj0 points on either side with smaller intensity. If negative,
  !>   use the default value (2).
  !> - pkinput: use a peak list as the initial model (nadj0 and ymax_detect0
  !>   are not used if pkinput is given).
  !> - xorig/yorig: return the pattern profile read from the file.
  !> - ycalc: return the fitted pattern.
  module subroutine xrpd_peaks_from_profile_file(p,xyfile,rms,errmsg,verbose0,&
     ymax_detect0,nadj0,pkinput,xorig,yorig,ycalc)
#ifdef HAVE_NLOPT
    use param, only: pi
    use tools, only: qcksort
    use tools_io, only: uout, file_read_xy, string
    use types, only: realloc
#endif
    class(xrpd_peaklist), intent(inout) :: p
    character*(*), intent(in) :: xyfile
    real*8, intent(out) :: rms
    character(len=:), allocatable, intent(out) :: errmsg
    logical, intent(in), optional :: verbose0
    real*8, intent(in), optional :: ymax_detect0
    integer, intent(in), optional :: nadj0
    type(xrpd_peaklist), intent(in), optional :: pkinput
    real*8, intent(inout), allocatable, optional :: xorig(:)
    real*8, intent(inout), allocatable, optional :: yorig(:)
    real*8, intent(inout), allocatable, optional :: ycalc(:)

    ! bail out if NLOPT is not available
#ifndef HAVE_NLOPT
    errmsg = "gaussian_compare can only be used if nlopt is available"
    rms = huge(1d0)
    return
#else

    integer :: ii, i, ip, n, nprm, nprml
    integer, allocatable :: pid(:), io(:)
    real*8, allocatable :: x(:), y(:), pth2(:), phei(:), prm(:), lb(:), ub(:), yfit(:), ysum(:)
    real*8, allocatable :: yread(:)
    real*8 :: fac, minx, maxx, maxy, ssq, maxa, xdif, xshift
    logical :: ok
    integer :: npeaks, npeaks_
    integer*8 :: opt
    integer :: ires
    logical :: verbose
    real*8 :: ymax_detect
    integer :: nadj

    include 'nlopt.f'

    integer, parameter :: main_algorithm = NLOPT_LD_CCSAQ
    integer, parameter :: fallback_algorithm = NLOPT_LD_MMA
    real*8, parameter :: ftol_eps = 1d-8
    real*8, parameter :: xtol_eps_abs(4) = (/1d-4,1d-8,1d-8,1d-4/) ! 2th, fwhm, eta, Int
    real*8, parameter :: fwhm_max = 0.6 ! maximum peak FWHM
    real*8, parameter :: fwhm_max_prefit = 0.1 ! maximum peak FWHM in prefit
    real*8, parameter :: area_peak_filter = 1d-4
    real*8, parameter :: gamma_default = 0.1d0
    real*8, parameter :: eta_default = 0.0d0
    integer, parameter :: prefit_percentile = 4

    ! process input
    errmsg = ""
    verbose = .false.
    if (present(verbose0)) verbose = verbose0
    ymax_detect = -huge(1d0)
    if (present(ymax_detect0)) ymax_detect = ymax_detect0
    nadj = 2
    if (present(nadj0)) nadj = nadj0

    ! read the pattern
    if (verbose) &
       write (uout,'("+ Reading the pattern from xy file: ",A)') trim(xyfile)
    call file_read_xy(xyfile,n,x,y,errmsg)
    if (len_trim(errmsg) > 0) return
    yread = y
    minx = minval(x)
    maxx = maxval(x)
    maxy = maxval(y)

    ! if default ymax_detect, use the median of the data
    if (ymax_detect == -huge(1d0)) then
       ! straightforward calculation of the median
       ysum = y
       allocate(io(n))
       do i = 1, n
          io(i) = i
       end do
       call qcksort(ysum,io,1,n)
       if (mod(n,2) == 0) then
          ymax_detect = 0.5d0 * (ysum(io(n/2)) + ysum(io(n/2+1)))
       else
          ymax_detect = ysum(io(n/2+1))
       end if
       deallocate(io,ysum)
    end if

    if (.not.present(pkinput)) then
       ! auto-detect the peaks
       if (verbose) &
          write (uout,'("+ Auto-detecting the peaks")')
       npeaks = 0
       allocate(pth2(100),pid(100),phei(100))
       if (verbose) then
          write (uout,'("  Minimum intensity cutoff for peak detection: ",A)') &
             string(ymax_detect,'f',decimal=4)
       end if
       do i = 3, n-2
          if (nadj == 2) then
             ok = (y(i) > ymax_detect .and. y(i) > y(i-1) .and. y(i) > y(i-2) .and. y(i) > y(i+1) .and.&
                y(i) > y(i+2))
          else
             ok = (y(i) > ymax_detect .and. y(i) > y(i-1) .and. y(i) > y(i+1))
          end if
          if (ok) then
             npeaks = npeaks + 1
             if (npeaks > size(pth2,1)) then
                call realloc(pth2,2*npeaks)
                call realloc(pid,2*npeaks)
                call realloc(phei,2*npeaks)
             end if
             pth2(npeaks) = x(i)
             phei(npeaks) = y(i)
             pid(npeaks) = i
          end if
       end do
       call realloc(pth2,npeaks)
       call realloc(pid,npeaks)
       call realloc(phei,npeaks)
       if (verbose) &
          write (uout,'("  Peaks detected: ",A)') string(npeaks)

       ! sort the peaks from highest to lowest
       allocate(io(npeaks))
       do i = 1, npeaks
          io(i) = i
       end do
       call qcksort(phei,io,1,npeaks)
       pth2 = pth2(io(npeaks:1:-1))
       phei = phei(io(npeaks:1:-1))
       pid = pid(io(npeaks:1:-1))
       deallocate(io)

       ! prepare the parameter list
       fac = eta_default * 2 * sqrt(log(2d0)) / (sqrt(pi) * gamma_default) + &
          (1-eta_default) * 2 / (pi * gamma_default)
       allocate(prm(4*npeaks),lb(4*npeaks),ub(4*npeaks))
       nprm = 0
       do i = 1, npeaks
          ! peak position (2*theta)
          nprm = nprm + 1
          prm(nprm) = pth2(i)
          if (pid(i) == 1) then
             lb(nprm) = x(1) - 1d-10
             ub(nprm) = pth2(i) + 2 * (x(pid(i)+1) - pth2(i))
          elseif (pid(i) == n) then
             lb(nprm) = pth2(i) - 2 * (pth2(i) - x(pid(i)-1))
             ub(nprm) = x(n) + 1d-10
          else
             lb(nprm) = pth2(i) - 2 * (pth2(i) - x(pid(i)-1))
             ub(nprm) = pth2(i) + 2 * (x(pid(i)+1) - pth2(i))
          end if

          ! peak FWHM (gamma)
          nprm = nprm + 1
          prm(nprm) = gamma_default
          lb(nprm) = 1d-5
          ub(nprm) = fwhm_max
          if (i <= npeaks/prefit_percentile) then
             ub(nprm) = fwhm_max_prefit
          else
             ub(nprm) = fwhm_max
          end if

          ! Gaussian/Lorentz coefficient (eta)
          nprm = nprm + 1
          prm(nprm) = eta_default
          lb(nprm) = 0d0
          ub(nprm) = 1d0

          ! peak area
          nprm = nprm + 1
          prm(nprm) = phei(i) / fac
          lb(nprm) = 0d0
          ub(nprm) = (maxx-minx) * maxy
       end do

       ! pre-fit the peaks
       if (verbose)&
          write (uout,'(/"+++ Pre-fitting peaks (this make take some time...)")')
       allocate(yfit(n),ysum(n))
       ysum = 0d0
       do ip = 1, npeaks
          ! number of parameters
          nprml = 4

          ! prepare fitting data and last peak
          y = yread - ysum
          prm(nprm) = max(y(pid(ip)) / fac,1d-20)
          lb(nprm) = 0d0
          ub(nprm) = (maxx-minx) * maxy

          ! run the minimization (main algorithm)
          call nlo_create(opt, main_algorithm, nprml)
          call nlo_set_xtol_abs(ires, opt, xtol_eps_abs)
          call nlo_set_lower_bounds(ires, opt, lb(4*(ip-1)+1:4*ip))
          call nlo_set_upper_bounds(ires, opt, ub(4*(ip-1)+1:4*ip))
          call nlo_set_min_objective(ires, opt, ffit, 0)
          call nlo_optimize(ires, opt, prm(4*(ip-1)+1:4*ip), ssq)
          call nlo_destroy(opt)
          if (ires < 0) then
             ! run the minimization (fallback algorithm)
             call nlo_create(opt, fallback_algorithm, nprml)
             call nlo_set_xtol_abs(ires, opt, xtol_eps_abs)
             call nlo_set_lower_bounds(ires, opt, lb(4*(ip-1)+1:4*ip))
             call nlo_set_upper_bounds(ires, opt, ub(4*(ip-1)+1:4*ip))
             call nlo_set_min_objective(ires, opt, ffit, 0)
             call nlo_optimize(ires, opt, prm(4*(ip-1)+1:4*ip), ssq)
             call nlo_destroy(opt)
          end if

          ! output message
          if (ires < 0) then
             if (ires == -1) then
                errmsg = "FAILURE"
             elseif (ires == -2) then
                errmsg = "FAILURE: invalid arguments"
             elseif (ires == -3) then
                errmsg = "FAILURE: out of memory"
             elseif (ires == -4) then
                errmsg = "FAILURE: roundoff errors limited progress"
             elseif (ires == -5) then
                errmsg = "FAILURE: termination forced by user"
             end if
             if (verbose) &
                write (uout,'(A)') errmsg
             return
          end if

          ! output peak
          rms = sqrt(ssq / real(n,8))
          if (verbose) then
             write (uout,'("Peak ",A,": center=",A," fwhm=",A," gaussian_ratio=",A," area=",A," rms=",A)') &
                string(ip), string(prm(4*(ip-1)+1),'f',decimal=4), string(prm(4*(ip-1)+2),'f',decimal=4),&
                string(prm(4*(ip-1)+3),'f',decimal=4), string(prm(4*(ip-1)+4),'f',decimal=4),&
                string(rms,'e',decimal=4)
          end if

          ! calculate final profile and sum it
          call ffit(ssq,4,prm(4*(ip-1)+1:4*ip),prm(4*(ip-1)+1:4*ip),0,0.)
          ysum = ysum + yfit

          ! lu = fopen_write("fit.dat")
          ! write (lu,'("## x y yfit std-resid")')
          ! do i = 1, n
          !    write (lu,'(3(A," "))') string(x(i),'f',decimal=10), string(yread(i),'f',decimal=10),&
          !       string(ysum(i),'f',decimal=10)
          ! end do
          ! call fclose(lu)
       end do
       y = yread
       rms = sqrt(sum((ysum - yread)**2) / real(n,8))
       if (ires < 0) then
          errmsg = "FAILURE in the prefitting stage"
          return
       end if
       write (uout,'("+ Finished pre-fitting peaks.")')
       write (*,'("+ RMS of the fit (after prefit) = ",A)') string(rms,'f',decimal=4)

       ! prune the peaks
       maxa = maxval(prm(4:nprm:4))
       npeaks_ = 0
       do i = 1, npeaks
          if (abs(prm(4*(i-1)+4) / maxa * 100d0) > area_peak_filter) then
             npeaks_ = npeaks_ + 1
             prm(4*(npeaks_-1)+1:4*npeaks_) = prm(4*(i-1)+1:4*i)
             lb(4*(npeaks_-1)+1:4*npeaks_) = lb(4*(i-1)+1:4*i)
             ub(4*(npeaks_-1)+1:4*npeaks_) = ub(4*(i-1)+1:4*i)
          end if
       end do
       npeaks = npeaks_
       nprm = 4*npeaks
       call realloc(prm,nprm)
       call realloc(lb,nprm)
       call realloc(ub,nprm)
       if (verbose) &
          write (uout,'("+ Number of peaks after pruning: ",A/)') string(npeaks)

       ! put back the maximum FWHM for all peaks and set the minimum at
       ! one average distance between x-points
       xdif = (maxx - minx) / real(n-1,8)
       do i = 1, npeaks
          ub(4*(i-1)+2) = fwhm_max
          lb(4*(i-1)+2) = min(xdif,0.99d0*fwhm_max)
          prm(4*(i-1)+2) = min(max(prm(4*(i-1)+2),lb(4*(i-1)+2)),ub(4*(i-1)+2))
       end do
    else
       ! sort the peaks from highest to lowest
       npeaks = pkinput%npeak
       allocate(io(npeaks),phei(npeaks))
       do i = 1, npeaks
          io(i) = i
          phei(i) = -pkinput%ip(i)
       end do
       call qcksort(phei,io,1,npeaks)
       deallocate(phei)

       ! prepare the parameter list
       fac = eta_default * 2 * sqrt(log(2d0)) / (sqrt(pi) * gamma_default) + &
          (1-eta_default) * 2 / (pi * gamma_default)
       allocate(prm(4*npeaks),lb(4*npeaks),ub(4*npeaks))
       nprm = 0
       xshift = sum(x(2:n) - x(1:n-1)) / real(n-1,8)
       do ii = 1, npeaks
          i = io(ii)
          ! peak position (2*theta)
          nprm = nprm + 1
          prm(nprm) = pkinput%th2(i)
          lb(nprm) = prm(nprm) - 2*xshift
          ub(nprm) = prm(nprm) + 2*xshift

          ! peak FWHM (gamma)
          nprm = nprm + 1
          if (pkinput%havepeakshape) then
             prm(nprm) = pkinput%fwhm(i)
          else
             prm(nprm) = gamma_default
          end if
          lb(nprm) = 1d-5
          ub(nprm) = fwhm_max

          ! Gaussian/Lorentz coefficient (eta)
          nprm = nprm + 1
          if (pkinput%havepeakshape) then
             prm(nprm) = pkinput%cgau(i)
          else
             prm(nprm) = 0d0
          end if
          lb(nprm) = 0d0
          ub(nprm) = 1d0

          ! peak area
          nprm = nprm + 1
          prm(nprm) = pkinput%ip(i)
          lb(nprm) = 0d0
          ub(nprm) = (maxx-minx) * maxy
       end do
       deallocate(io)
    end if

    ! fitting the whole pattern
    if (verbose) &
       write (uout,'("+++ Fitting the whole pattern")')
    ! run the minimization (main algorithm)
    call nlo_create(opt, main_algorithm, nprm)
    call nlo_set_ftol_rel(ires, opt, ftol_eps)
    call nlo_set_lower_bounds(ires, opt, lb)
    call nlo_set_upper_bounds(ires, opt, ub)
    call nlo_set_min_objective(ires, opt, ffit, 0)
    call nlo_optimize(ires, opt, prm, ssq)
    call nlo_destroy(opt)
    if (ires < 0) then
       call nlo_create(opt, fallback_algorithm, nprm)
       call nlo_set_ftol_rel(ires, opt, ftol_eps)
       call nlo_set_lower_bounds(ires, opt, lb)
       call nlo_set_upper_bounds(ires, opt, ub)
       call nlo_set_min_objective(ires, opt, ffit, 0)
       call nlo_optimize(ires, opt, prm, ssq)
       call nlo_destroy(opt)
    end if

    ! output message
    if (ires < 0) then
       call errmsg_from_ires()
       if (verbose) &
          write (uout,'(A)') errmsg
       return
    end if

    ! prune the peaks again
    maxa = maxval(prm(4:nprm:4))
    npeaks_ = 0
    do i = 1, npeaks
       if (abs(prm(4*(i-1)+4) / maxa * 100d0) > area_peak_filter) then
          npeaks_ = npeaks_ + 1
          prm(4*(npeaks_-1)+1:4*npeaks_) = prm(4*(i-1)+1:4*i)
          lb(4*(npeaks_-1)+1:4*npeaks_) = lb(4*(i-1)+1:4*i)
          ub(4*(npeaks_-1)+1:4*npeaks_) = ub(4*(i-1)+1:4*i)
       end if
    end do
    npeaks = npeaks_
    nprm = 4*npeaks
    call realloc(prm,nprm)
    call realloc(lb,nprm)
    call realloc(ub,nprm)
    if (verbose) &
       write (uout,'("+ Number of peaks after final pruning: ",A)') string(npeaks)

    ! peaks to standard output
    if (verbose) then
       write (uout,'("+ Final list of peaks:")')
       maxa = maxval(prm(4:nprm:4))
       do ip = 1, npeaks
          write (uout,'("Peak ",A,": center=",A," fwhm=",A," gaussian_ratio=",A," area=",A," norm_area=",A)') &
             string(ip), string(prm(4*(ip-1)+1),'f',decimal=4), string(prm(4*(ip-1)+2),'f',decimal=4),&
             string(prm(4*(ip-1)+3),'f',decimal=4), string(prm(4*(ip-1)+4),'f',decimal=4),&
             string(prm(4*(ip-1)+4)/maxa*100d0,'f',decimal=4)
       end do
    end if

    ! calculate final profile and rms
    call ffit(ssq,nprm,prm,prm,0,0.)
    rms = sqrt(sum((yfit - yread)**2) / real(n,8))
    if (verbose) then
       write (uout,'("+ Finished fitting pattern.")')
       write (*,'("+ RMS of the fit (final) = ",A/)') string(rms,'f',decimal=4)
    end if

    ! output the XRPD peak list
    call p%end()
    p%npeak = npeaks
    p%haveth2limits = .true.
    p%th2ini = x(1)
    p%th2end = x(n)
    p%havepeakshape = .true.
    allocate(p%th2(npeaks))
    allocate(p%ip(npeaks))
    allocate(p%fwhm(npeaks))
    allocate(p%cgau(npeaks))
    do ip = 1, npeaks
       p%th2(ip) = prm(4*(ip-1)+1)
       p%fwhm(ip) = prm(4*(ip-1)+2)
       p%cgau(ip) = prm(4*(ip-1)+3)
       p%ip(ip) = prm(4*(ip-1)+4)
    end do

    ! output the pattern
    if (present(xorig)) xorig = x
    if (present(yorig)) yorig = yread
    if (present(ycalc)) ycalc = yfit

  contains
    !> Helper routine: get the errormsg from the NLOPT error code.
    subroutine errmsg_from_ires()
       if (ires == -1) then
          errmsg = "FAILURE"
       elseif (ires == -2) then
          errmsg = "FAILURE: invalid arguments"
       elseif (ires == -3) then
          errmsg = "FAILURE: out of memory"
       elseif (ires == -4) then
          errmsg = "FAILURE: roundoff errors limited progress"
       elseif (ires == -5) then
          errmsg = "FAILURE: termination forced by user"
       end if
    end subroutine errmsg_from_ires

    !> Routine for NLOPT optimization. Returns sum of the squares of
    !> the difference in val for a model with nprm parameters in the
    !> vector prm. If need_gradient is non-zero, return the gradient
    !> of the parameters in grad. f_data is ignored. Side-effect:
    !> Writes the pattern corresponding to prm in the yfit vector
    !> from the parent routine.
    subroutine ffit(val, nprm, prm, grad, need_gradient, f_data)
      use tools_math, only: gaussian, lorentzian
      use param, only: pi
      real*8 :: val, prm(nprm), grad(nprm)
      integer :: nprm, need_gradient
      real :: f_data

      integer :: i
      real*8 :: x0, gamma, eta, int, s, g2
      real*8, allocatable :: gau(:), lor(:), yfitg(:,:), ydif(:)

      if (.not.allocated(yfit)) allocate(yfit(n))
      allocate(gau(n), lor(n), yfitg(n,nprm), ydif(n))

      yfit = 0d0
      do i = 1, nprm, 4
         x0 = prm(i)
         gamma = max(prm(i+1),1d-80)
         eta = prm(i+2)
         int = prm(i+3)

         gau = gaussian(x,x0,gamma)
         lor = lorentzian(x,x0,gamma)

         yfit = yfit + int * (eta * gau + (1-eta) * lor)
         if (need_gradient /= 0) then
            g2 = gamma / 2d0
            s = g2 / sqrt(2d0 * log(2d0))

            yfitg(:,i) = int * (eta * (gau * (x-x0) / (s * s)) +&
               (1-eta) * (lor * lor * 4 * pi * (x-x0) / gamma))
            yfitg(:,i+1) = int * (eta * (-gau / gamma + gau * (x-x0)*(x-x0) / (s*s) / gamma) +&
               (1-eta) * (lor / gamma - pi * lor * lor))
            yfitg(:,i+2) = int * (gau - lor)
            yfitg(:,i+3) = eta * gau + (1-eta) * lor
         end if
      end do

      ydif = y - yfit
      val = sum(ydif * ydif)
      if (need_gradient /= 0) then
         do i = 1, nprm
            grad(i) = -2d0 * sum(ydif * yfitg(:,i))
         end do
      end if

    end subroutine ffit
#endif
  end subroutine xrpd_peaks_from_profile_file

  !> Write the peak list to a file.
  module subroutine xrpd_write_to_file(p,file)
    use tools_io, only: fopen_write, fclose, string
    class(xrpd_peaklist), intent(in) :: p
    character*(*), intent(in) :: file

    integer :: lu, ip
    character(len=:), allocatable :: str

    lu = fopen_write(file)
    write (lu,'("## List of peaks")')
    write (lu,'("## 2*theta   Area   FWHM   gau/lor")')
    if (p%haveth2limits) then
       write (lu,'("@th2ini ",A)') string(p%th2ini,'f',decimal=8)
       write (lu,'("@th2end ",A)') string(p%th2end,'f',decimal=8)
    end if
    if (p%haveradiation) then
       write (lu,'("@lambda ",A)') string(p%lambda,'f',decimal=8)
       write (lu,'("@fpol ",A)') string(p%fpol,'f',decimal=8)
    end if
    do ip = 1, p%npeak
       str = string(p%th2(ip),'f',decimal=10) // " " // string(p%ip(ip),'f',decimal=10)
       if (p%havepeakshape) then
          str = str // " " // string(p%fwhm(ip),'f',decimal=10) // " " // string(p%cgau(ip),'f',decimal=10)
          if (p%havehvec) then
             str = str // string(p%hvec(1,ip)) // " " // string(p%hvec(2,ip)) // " " // string(p%hvec(3,ip))
          end if
       end if
       write (lu,'(A)') str
    end do
    call fclose(lu)

  end subroutine xrpd_write_to_file

  !> Calculate the XRPD profile corresponding to the peak list in
  !> p. Requires having the peak shapes. The profile contains n evenly
  !> distributed points: the angles (th2, degrees) are returned in x
  !> and the intensity is returned in y. If the pattern contains
  !> no th2ini and th2end info or if the user wants to override them,
  !> use th2ini and th2end to set the profile limits. Return non-zero
  !> errmsg on error.
  module subroutine xrpd_calculate_profile(p,n,x,y,errmsg,th2ini,th2end)
    use tools_math, only: gaussian, lorentzian
    class(xrpd_peaklist), intent(inout) :: p
    integer, intent(in) :: n
    real*8, allocatable, intent(inout) :: x(:), y(:)
    character(len=:), allocatable, intent(out) :: errmsg
    real*8, intent(in), optional :: th2ini, th2end

    integer :: i
    real*8 :: xini, xend, x0, gamma, eta, int

    ! consistency checks
    errmsg = ""
    if (n <= 0) then
       errmsg = "must have positive number of points for the grid"
       return
    end if

    ! initialize grid limits
    if (.not.present(th2ini).and..not.present(th2end)) then
       if (.not.p%haveth2limits) then
          errmsg = "must have th2ini and th2end from the procedure arguments"
          return
       end if
    end if
    if (p%haveth2limits) then
       xini = p%th2ini
       xend = p%th2end
    end if
    if (present(th2ini)) xini = th2ini
    if (present(th2end)) xend = th2end

    ! initialize output arrays
    if (allocated(x)) deallocate(x)
    if (allocated(y)) deallocate(y)
    allocate(x(n))
    allocate(y(n))
    do i = 1, n
       x(i) = xini + real(i-1,8) / real(n-1,8) * (xend-xini)
    end do

    ! calculate the profile
    y = 0d0
    do i = 1, p%npeak
       x0 = p%th2(i)
       if (p%havepeakshape) then
          gamma = p%fwhm(i)
          eta = p%cgau(i)
       else
          gamma = 2d0 * sqrt(2d0 * log(2d0)) * xrpd_sigma_def
          eta = 1d0
       end if
       int = p%ip(i)
       y = y + int * (eta * gaussian(x,x0,gamma) + (1-eta) * lorentzian(x,x0,gamma))
    end do

  end subroutine xrpd_calculate_profile

  ! Calculate the Gaussian cross-correlation between the peak patterns
  ! p1 and p2 with Gaussian triangle function width alpha and
  ! Gaussian width sigma, and output the value in d12. If calcderivs
  ! is present and true, calculate the derivatives wrt the metric
  ! tensor of p1 in d12g (requires p1%havegradients).
  module subroutine crosscorr_gaussian(p1,p2,alpha,sigma,d12,errmsg,calcderivs,d12g)
    use param, only: pi
    type(xrpd_peaklist), intent(in) :: p1, p2
    real*8, intent(in) :: alpha, sigma
    real*8, intent(out) :: d12
    character(len=:), allocatable, intent(out) :: errmsg
    logical, intent(in), optional :: calcderivs
    real*8, intent(out), optional :: d12g(6)

    logical :: calcg
    real*8 :: z, zp, z2p, zsq, thdeps, expt12, thdif
    integer :: imin, i, j

    real*8, parameter :: eps_discard = 1d-10

    ! consistency checks
    calcg = present(calcderivs)
    if (calcg) calcg = calcderivs
    if (calcg .and..not.p1%havegradients) then
       errmsg = 'calcderivs requires p1%havegradients'
       return
    end if
    if (calcg .and..not.present(d12g)) then
       errmsg = 'calcderivs requires d12g being present'
       return
    end if

    z = 1d0 / (alpha**2 + 4d0 * pi * sigma**2)
    zp = pi * z
    z2p = 2d0 * zp
    zsq = sqrt(z)

    thdeps = sqrt(abs(-log(eps_discard) / zp))
    d12 = 0d0
    imin = 1
    if (calcg) d12g = 0d0
    do j = 1, p2%npeak
       iloop: do i = imin, p1%npeak
          thdif = p1%th2(i) - p2%th2(j)
          if (abs(thdif) > thdeps) then
             if (thdif > thdeps) then
                exit iloop
             else
                imin = i+1
                cycle
             end if
          end if
          expt12 = exp(-zp * thdif * thdif)
          d12 = d12 + p1%ip(i) * p2%ip(j) * expt12
          if (calcg) d12g = d12g + p2%ip(j) * expt12 * (p1%ipg(:,i) - z2p * p1%ip(i) * thdif * p1%th2g(:,i))
       end do iloop
    end do
    d12 = d12 * zsq
    if (calcg) d12g = d12g * zsq

  end subroutine crosscorr_gaussian

  ! Compare crystal structure c1 against XRPD pattern p2 using
  ! Gaussian PWDF (GPWDF) and variants. imode must be one of 0 (no
  ! change to c1), 1 (local minization of diff), 2 (global
  ! minimization of diff).  The difference (diff) between the two is
  ! output in diff. In the case of local or global minimization, the
  ! seed for the deformed c1 structure is returned in seedout.
  ! Optional input arguments:
  ! - verbose0 = write messages to output
  ! - alpha0 = width of the triangle Gaussian function
  ! - lambda0 = wavefunction of the light in angstrom
  ! - fpol0 = polarization factor
  ! - maxfeval0 = maximum number of evaluations before stopping
  ! - besteps0 = in global minimization, do not reset the eval
  !   counter if difference wrt the lowest diff is lower than this
  !   value.
  ! - max_elong_def0 = maximum cell length deformation
  ! - max_ang_def0 = maximum angle deformation (degrees)
  module subroutine gaussian_compare(c1,p2,imode,diff,errmsg,seedout,verbose0,alpha0,&
     lambda0,fpol0,maxfeval0,besteps0,max_elong_def0,max_ang_def0)
    use crystalseedmod, only: crystalseed
#ifdef HAVE_NLOPT
    use tools_io, only: uout, string
    use tools_math, only: m_x2c_from_cellpar, det3
#endif
    type(crystal), intent(in) :: c1
    type(xrpd_peaklist), intent(in) :: p2
    integer, intent(in) :: imode
    real*8, intent(out) :: diff
    character(len=:), allocatable, intent(out) :: errmsg
    type(crystalseed), intent(out), optional :: seedout
    logical, intent(in), optional :: verbose0
    real*8, intent(in), optional :: alpha0
    real*8, intent(in), optional :: lambda0
    real*8, intent(in), optional :: fpol0
    integer, intent(in), optional :: maxfeval0
    real*8, intent(in), optional :: besteps0
    real*8, intent(in), optional :: max_elong_def0
    real*8, intent(in), optional :: max_ang_def0

    ! bail out if NLOPT is not available
#ifndef HAVE_NLOPT
    errmsg = "gaussian_compare can only be used if nlopt is available"
    diff = 1d0
    if (present(seedout)) call seedout%end()
    return
#else

    ! local block
    integer :: i
    real*8 :: dfg22, x(6), xorig(6), vorig, grad(6), x2c(3,3), omega
    integer*8 :: opt, lopt
    integer :: ires, maxfeval
    real*8 :: lb(6), ub(6), th2ini, th2end, alpha, lambda, besteps
    real*8 :: max_elong_def, max_ang_def, fpol
    logical :: iresok
    logical :: verbose
    type(xrpd_peaklist) :: p1

    ! global block
    integer :: neval
    real*8 :: lastval
    real*8 :: bestval
    integer :: nbesteval
    real*8, parameter :: ftol_eps = 1d-5

    ! parameters
    integer, parameter :: imode_sp = 0
    integer, parameter :: imode_local = 1
    integer, parameter :: imode_global = 2

    include 'nlopt.f'

    ! consistency checks
    errmsg = ""
    if (p2%npeak == 0) then
       errmsg = 'No peaks found in the pattern'
       return
    end if
    if (imode /= imode_sp .and. imode /= imode_local .and. imode /= imode_global) then
       errmsg = 'Incorrect imode (must be one of 0, 1, 2)'
       return
    end if

    ! process the input options
    alpha = xrpd_alpha_def
    if (present(alpha0)) alpha = alpha0
    lambda = xrpd_lambda_def
    if (present(lambda0)) lambda = lambda0
    maxfeval = xrpd_maxfeval_def_safe
    if (present(maxfeval0)) maxfeval = maxfeval0
    besteps = xrpd_besteps_def_safe
    if (present(besteps0)) besteps = besteps0
    max_elong_def = xrpd_max_elong_def_safe
    if (present(max_elong_def0)) max_elong_def = max_elong_def0
    max_ang_def = xrpd_max_ang_def_safe
    if (present(max_ang_def0)) max_ang_def = max_ang_def0
    fpol = xrpd_fpol_def
    if (present(fpol0)) fpol = fpol0
    verbose = .true.
    if (present(verbose0)) verbose = verbose0

    ! initialize
    neval = 0
    lastval = -1d0
    bestval = 1.1d0
    nbesteval = 0
    th2ini = p2%th2(1) - 1d-2
    th2end = p2%th2(p2%npeak) + 1d-2
    if (verbose) &
       write (uout,'("# step    DIFF        -- cell parameters --")')

    call p1%from_crystal(c1,th2ini,th2end,lambda,fpol,.false.,.false.,errmsg)
    if (len_trim(errmsg) > 0) return
    call crosscorr_gaussian(p2,p2,alpha,xrpd_sigma_def,dfg22,errmsg,.false.)
    if (len_trim(errmsg) > 0) return

    x(1:3) = c1%aa
    x(4:6) = c1%bb
    xorig = x
    vorig = c1%omega
    if (imode /= imode_sp) then
       if (imode == imode_global) then
          ! global minimization
          call nlo_create(lopt, NLOPT_LD_SLSQP, 6)
          call nlo_set_ftol_rel(ires, lopt, ftol_eps)

          call nlo_create(opt, NLOPT_G_MLSL_LDS, 6)
          call nlo_set_local_optimizer(ires, opt, lopt)
       elseif (imode == imode_local) then
          ! local minimization
          call nlo_create(opt, NLOPT_LD_SLSQP, 6)
          call nlo_set_ftol_rel(ires, opt, ftol_eps)
       end if

       lb(1:3) = max(x(1:3) * (1 - max_elong_def),0d0)
       ub(1:3) = x(1:3) * (1 + max_elong_def)
       lb(4:6) = max(x(4:6) - max_ang_def,10d0)
       ub(4:6) = min(x(4:6) + max_ang_def,170d0)
       call nlo_set_lower_bounds(ires, opt, lb)
       call nlo_set_upper_bounds(ires, opt, ub)
       call nlo_set_min_objective(ires, opt, diff_fun, 0)

       ! run the minimization
       iresok = .false.
       call nlo_optimize(ires, opt, x, diff)

       ! final message
       if (verbose) then
          if (ires == 1) then
             write (uout,'("+ SUCCESS")')
          elseif (ires == 2) then
             write (uout,'("+ SUCCESS: maximum iterations reached")')
          elseif (ires == 3) then
             write (uout,'("+ SUCCESS: ftol_rel or ftol_abs was reached")')
          elseif (ires == 4) then
             write (uout,'("+ SUCCESS: xtol_rel or xtol_abs was reached")')
          elseif (ires == 5) then
             write (uout,'("+ SUCCESS: maxeval was reached")')
          elseif (ires == 6) then
             write (uout,'("+ SUCCESS: maxtime was reached")')
          elseif (ires == -1) then
             write (uout,'("+ FAILURE")')
          elseif (ires == -2) then
             write (uout,'("+ FAILURE: invalid arguments")')
          elseif (ires == -3) then
             write (uout,'("+ FAILURE: out of memory")')
          elseif (ires == -4) then
             write (uout,'("+ FAILURE: roundoff errors limited progress")')
          elseif (ires == -5) then
             if (iresok) then
                write (uout,'("+ SUCCESS? maximum number of evaluations reached")')
             else
                write (uout,'("+ FAILURE: termination forced by user")')
             end if
          end if
       end if
       if (.not.iresok .and. ires < 0) then
          errmsg = "Error during the minimization process"
          return
       end if

       ! clean up
       if (imode == imode_global) &
          call nlo_destroy(lopt)
       call nlo_destroy(opt)

       ! write message to output
       if (verbose) then
          write (uout,'("+ Lattice parameters: ")')
          write (uout,'("  Initial (1): ",6(A," "))') &
             (string(xorig(i),'f',length=11,decimal=8),i=1,3), &
             (string(xorig(i),'f',length=10,decimal=6),i=4,6)
          write (uout,'("  Final (1):   ",6(A," "))') &
             (string(x(i),'f',length=11,decimal=8),i=1,3), &
             (string(x(i),'f',length=10,decimal=6),i=4,6)
          write (uout,'("  Relative length deformations:   ",3(A," "))') &
             (string(abs(xorig(i)-x(i))/xorig(i),'f',length=11,decimal=8),i=1,3)
          write (uout,'("  Angle displacements:   ",3(A," "))') &
             (string(abs(xorig(i)-x(i)),'f',length=7,decimal=4),i=4,6)

          x2c = m_x2c_from_cellpar(x(1:3),x(4:6))
          omega = det3(x2c)
          write (uout,'("  Initial volume (bohr3): ",A)') string(vorig,'f',decimal=4)
          write (uout,'("  Final volume (bohr3): ",A)') string(omega,'f',decimal=4)
          write (uout,'("  Volume deformation: ",A)') string(abs(omega-vorig)/vorig,'f',decimal=8)
       end if

       ! make final structure seed
       if (present(seedout)) then
          call c1%makeseed(seedout,.false.)
          seedout%useabr = 1
          seedout%aa = x(1:3)
          seedout%bb = x(4:6)
       end if
    else
       call diff_fun(diff,6,x,grad,0,1.)
    end if

  contains
    subroutine diff_fun(val, n, x, grad, need_gradient, f_data)
      use param, only: pi
      use tools_io, only: string, uout, ioj_left, ioj_right
      use tools_math, only: det3sym
      real*8 :: val, x(n), grad(n)
      integer :: n, need_gradient
      real :: f_data

      real*8 :: calp, cbet, cgam, a, b, c, salp, sbet, sgam, gg(3,3)
      real*8 :: dd, ddg(6)
      real*8 :: dfg11, dfgg11(6), dfg12, dfgg12(6)
      integer :: i, ires
      character(len=:), allocatable :: str, errmsg

      neval = neval + 1
      a = x(1)
      b = x(2)
      c = x(3)
      calp = cos(x(4) * pi / 180d0)
      cbet = cos(x(5) * pi / 180d0)
      cgam = cos(x(6) * pi / 180d0)
      salp = sin(x(4) * pi / 180d0)
      sbet = sin(x(5) * pi / 180d0)
      sgam = sin(x(6) * pi / 180d0)
      gg(1,1) = a * a
      gg(1,2) = a * b * cgam
      gg(2,1) = gg(1,2)
      gg(1,3) = a * c * cbet
      gg(3,1) = gg(1,3)
      gg(2,2) = b * b
      gg(2,3) = b * c * calp
      gg(3,2) = gg(2,3)
      gg(3,3) = c * c

      if (det3sym(gg) < 0d0) then
         val = huge(1d0)
         if (need_gradient /= 0) grad = 0d0
         return
      end if

      ! only recompute peak pattern for crystal 1
      call p1%from_crystal(c1,th2ini,th2end,lambda,fpol,.true.,.true.,errmsg,gg)
      if (len_trim(errmsg) > 0) goto 999
      call crosscorr_gaussian(p1,p1,alpha,xrpd_sigma_def,dfg11,errmsg,.true.,dfgg11)
      if (len_trim(errmsg) > 0) goto 999
      call crosscorr_gaussian(p1,p2,alpha,xrpd_sigma_def,dfg12,errmsg,.true.,dfgg12)
      if (len_trim(errmsg) > 0) goto 999

      dd = dfg12 / sqrt(dfg11 * dfg22)

      ! write output
      val = 1d0 - dd
      if (need_gradient /= 0) then
         ! derivatives wrt Gij (the 0.5 in the second term is missing
         ! because there are two dfgg11, one from each "1"
         ddg = -dd * (dfgg12 / dfg12 - dfgg11 / dfg11)

         ! Off-diagonal components (2,3,5) are half what they should
         ! be because we did not impose symmetric matrix
         grad(1) = a * ddg(1) + b * cgam * ddg(2) + c * cbet * ddg(3)
         grad(2) = a * cgam * ddg(2) + b * ddg(4) + c * calp * ddg(5)
         grad(3) = a * cbet * ddg(3) + b * calp * ddg(5) + c * ddg(6)
         grad(4) = - b * c * salp * ddg(5) * pi / 180d0
         grad(5) = - a * c * sbet * ddg(3) * pi / 180d0
         grad(6) = - a * b * sgam * ddg(2) * pi / 180d0
         grad = 2 * grad
      end if

      ! message
      lastval = val
      str = ""
      if (val < bestval * (1d0 - besteps)) then
         bestval = val
         nbesteval = neval
      elseif (verbose) then
         str = " | " // string(bestval,'f',12,8)  // "best for last " // string(neval - nbesteval) // " iterations"
      end if
      if (verbose) then
         write (uout,'(A," ",A," at ",7(A," "))') string(neval,8,ioj_left),&
            string(val,'f',length=12,decimal=8),&
            (string(x(i),'f',length=7,decimal=4,justify=ioj_right),i=1,3), &
            (string(x(i),'f',length=6,decimal=2,justify=ioj_right),i=4,6), &
            str
      end if

      ! force termination?
      if (imode == imode_global .and. neval - nbesteval > maxfeval) then
         iresok = .true.
         call nlo_set_force_stop(ires, opt, 2)
      end if
      return
999   continue

      iresok = .false.
      call nlo_set_force_stop(ires, opt, -1)

    end subroutine diff_fun
#endif
  end subroutine gaussian_compare

  !> Calculate the background for a diffraction pattern with n points
  !> given by 2*theta (x) and profile intensity (y) and return it in
  !> yb(n). The x must be in strict increasing order. The background
  !> is estimated using robust Bayesian analysis, as proposed in David
  !> and Sivia, J. Appl. Cryst. 34 (2001) 318.
  module function david_sivia_calculate_background(n,x,y,errmsg) result(yb)
    use tools_math, only: splinefit, splineval
    use tools_io, only: fopen_write, fclose
    integer, intent(in) :: n
    real*8, intent(in) :: x(n)
    real*8, intent(in) :: y(n)
    character(len=:), allocatable, intent(out) :: errmsg
    real*8 :: yb(n)
#ifndef HAVE_NLOPT
    errmsg = "david_sivia_calculate_background can only be used if nlopt is available"
    yb = 0d0
    return
#else
    integer :: i, nk, j, lu
    real*8, allocatable :: xk(:), yk(:), c(:,:)
    real*8 :: ssq
    integer*8 :: opt
    integer :: ires

    integer, parameter :: nknot0 = 25
    real*8, parameter :: ftol_eps = 1d-8

    include 'nlopt.f'

    ! set up the initial knots
    errmsg = ""
    nk = max(min(nknot0,n/2),2)
    nk = 20 ! xxxxx

    allocate(xk(nk),yk(nk))
    do i = 1, nk
       xk(i) = x(1) + real(i - 1,8) / real(nk - 1,8) * (x(n) - x(1))
    end do
    yk = sum(y) / real(n,8)

    ! optimize
    call nlo_create(opt, NLOPT_LN_PRAXIS, nk)
    call nlo_set_ftol_rel(ires, opt, ftol_eps)
    call nlo_set_min_objective(ires, opt, lsfun, 0)
    call nlo_optimize(ires, opt, yk, ssq)
    call nlo_destroy(opt)

    ! output message and destroy
    if (ires < 0) then
       if (ires == -1) then
          errmsg = "FAILURE"
       elseif (ires == -2) then
          errmsg = "FAILURE: invalid arguments"
       elseif (ires == -3) then
          errmsg = "FAILURE: out of memory"
       elseif (ires == -4) then
          errmsg = "FAILURE: roundoff errors limited progress"
       elseif (ires == -5) then
          errmsg = "FAILURE: termination forced by user"
       end if
       return
    end if

    ! output the final background pattern
    c = splinefit(nk,xk,yk)
    do i = 1, n
       yb(i) = splineval(nk,c,x(i))
    end do

  contains
    subroutine lsfun(val, nin, yin, grad, need_gradient, f_data)
      real*8 :: val, yin(nin), grad(nin)
      integer :: nin, need_gradient
      real :: f_data

      integer :: i
      real*8 :: t, diff

      c = splinefit(nk,xk,yin)
      val = 0d0
      do i = 1, n
         t = (y(i) - splineval(nk,c,x(i))) / sqrt(y(i))
         diff = loglike_bayes(t)
         val = val - diff
      end do

    end subroutine lsfun
#endif
  end function david_sivia_calculate_background

  !xx! private procedures

  !> Returns the log(likelihood) function for the robust estimation of
  !> background in a diffraction pattern, from David and Sivia
  !> J. Appl. Cryst. 34 (2001) 318. The input t must be equal to
  !> (y_bckgrnd - y_obs) / (sqrt(2)*sigma). This routine implements
  !> equation 12 in the article:
  !>   f(t) =
  !>      -t^2, if t <= 0
  !>      A + ln( int_eps^infty exp(-(u-t)^2)/u du ), if 0 < t <= 10
  !>      B + ln(t), if t > 10
  !> The values in the t <= 0 and t > 10 ranges correspond to the
  !> asymptotic behavior of the log-integral function. The A and B
  !> values are chosen to make the function continuous.
  !> Note: the integrand is ill-behaved in the u -> 0 limit, and the
  !> shape of the log(integral) is heavily dependent on a) the value
  !> of eps and b) the way the numerical integration is carried
  !> out. The spline coefficients in this routine were determined with
  !> a) eps = 4e-6 and b) a trapezoidal rule with 100001 points between
  !> eps and t + 10.
  function loglike_bayes(t)
    use tools_math, only: splineval
    real*8, intent(in) :: t
    real*8 :: loglike_bayes

    integer :: i, j

    real*8, allocatable :: caux(:,:)

    integer, parameter :: nknot = 41
    real*8, parameter :: c(5,nknot) = reshape((/&
       -0.00000000000000E+00,1.34920863570493E-01,-9.73650027881392E-01,4.44282155173727E-02,0.00000000000000E+00, &
       -2.64287199825048E-02,-3.43573859960696E-01,-9.40328866243362E-01,4.44282155173710E-02,2.50000000000000E-01,&
       -1.70398548245430E-01,-8.05408002672870E-01,-9.07007704605332E-01,9.05752872736656E-02,5.00000000000000E-01,&
       -4.27023291587830E-01,-1.24192898861172E+00,-8.39076239150080E-01,1.67533612034369E-01,7.50000000000000E-01,&
       -7.87330090999603E-01,-1.63005455593032E+00,-7.13426030124303E-01,3.13952473145982E-01,1.00000000000000E+00,&
       -1.23452734947205E+00,-1.92790148227760E+00,-4.77961675264813E-01,5.45173798848488E-01,1.25000000000000E+00,&
       -1.73785698413849E+00,-2.06466223262591E+00,-6.90813261284493E-02,7.57202459633888E-01,1.50000000000000E+00,&
       -2.24650883674622E+00,-1.95722743450878E+00,4.98820518596971E-01,5.91490300604228E-01,1.75000000000000E+00, &
       -2.69539737701416E+00,-1.59691274384701E+00,9.42438244050142E-01,-2.94975853906578E-02,2.00000000000000E+00,&
       -3.03618407249451E+00,-1.13122441908268E+00,9.20315055007149E-01,-4.40281575252538E-01,2.25000000000000E+00,&
       -3.26834988594055E+00,-7.53619686938960E-01,5.90103873567745E-01,-3.81281631177316E-01,2.50000000000000E+00,&
       -3.42583084106445E+00,-5.30058056000834E-01,3.04142650184758E-01,-2.00927959608507E-01,2.75000000000000E+00,&
       -3.54247593879700E+00,-4.15660723335050E-01,1.53446680478378E-01,-8.67350460136542E-02,3.00000000000000E+00,&
       -3.63815593719482E+00,-3.55200204223421E-01,8.83953959681358E-02,-3.65892537978070E-02,3.25000000000000E+00,&
       -3.72200298309326E+00,-3.17862991326442E-01,6.09534556197804E-02,-1.78234661388643E-02,3.50000000000000E+00,&
       -3.79793763160706E+00,-2.90728163417589E-01,4.75858560156328E-02,-1.07997429748607E-02,3.75000000000000E+00,&
       -3.86781430244446E+00,-2.68960187217559E-01,3.94860487844879E-02,-7.64211274294979E-03,4.00000000000000E+00,&
       -3.93270587921143E+00,-2.50650058964618E-01,3.37544642272753E-02,-5.91879335802759E-03,4.25000000000000E+00,&
       -3.99335122108459E+00,-2.34882600605610E-01,2.93153692087539E-02,-4.72397359056131E-03,4.50000000000000E+00,&
       -4.05031347274780E+00,-2.21110661049464E-01,2.57723890158322E-02,-3.85657692815933E-03,4.75000000000000E+00,&
       -4.10404062271118E+00,-2.08947574715578E-01,2.28799563197128E-02,-3.24144232961121E-03,5.00000000000000E+00,&
       -4.15489816665649E+00,-1.98115366992523E-01,2.04488745725047E-02,-2.70890375339761E-03,5.25000000000000E+00,&
       -4.20319128036499E+00,-1.88398849160033E-01,1.84171967574561E-02,-2.31088210992247E-03,5.50000000000000E+00,&
       -4.24917602539062E+00,-1.79623541176915E-01,1.66840351750149E-02,-1.99410101003927E-03,5.75000000000000E+00,&
       -4.29307031631470E+00,-1.71655417528790E-01,1.51884594174855E-02,-1.73663963117310E-03,6.00000000000000E+00,&
       -4.33506202697754E+00,-1.64386807750892E-01,1.38859796941057E-02,-1.49635218401833E-03,6.25000000000000E+00,&
       -4.37531423568726E+00,-1.57724383938343E-01,1.27637155560920E-02,-1.34167233587945E-03,6.50000000000000E+00,&
       -4.41396856307983E+00,-1.51594089723274E-01,1.17574613041821E-02,-1.16308151933708E-03,6.75000000000000E+00,&
       -4.45115041732788E+00,-1.45933436856059E-01,1.08851501646792E-02,-1.05556213364544E-03,7.00000000000000E+00,&
       -4.48696994781494E+00,-1.40688779673778E-01,1.00934785644456E-02,-9.31808617958829E-04,7.25000000000000E+00,&
       -4.52152585983276E+00,-1.35816754507423E-01,9.39462210097652E-03,-8.32437769521022E-04,7.50000000000000E+00,&
       -4.55490589141846E+00,-1.31275525538719E-01,8.77029377383598E-03,-7.73840694582972E-04,7.75000000000000E+00,&
       -4.58718872070312E+00,-1.27035473782036E-01,8.18991325289864E-03,-6.49836170897089E-04,8.00000000000000E+00,&
       -4.61844587326050E+00,-1.23062361437630E-01,7.70253612472593E-03,-6.55134934329116E-04,8.25000000000000E+00,&
       -4.64874029159546E+00,-1.19333931175453E-01,7.21118492397932E-03,-5.44321357412780E-04,8.50000000000000E+00,&
       -4.67813158035278E+00,-1.15830398967979E-01,6.80294390591984E-03,-5.24513229771095E-04,8.75000000000000E+00,&
       -4.70667219161987E+00,-1.12527273245601E-01,6.40955898359147E-03,-4.70418692253283E-04,9.00000000000000E+00,&
       -4.73441076278687E+00,-1.09410697258603E-01,6.05674496440178E-03,-4.31946766841662E-04,9.25000000000000E+00,&
       -4.76139163970947E+00,-1.06463314795184E-01,5.73278488927076E-03,-3.95788381006845E-04,9.50000000000000E+00,&
       -4.78765535354614E+00,-1.03671132671988E-01,5.43594360351551E-03,-3.95788381006401E-04,9.75000000000000E+00,&
       -4.81323957443237E+00,-1.01027371191669E-01,2.50000000000000E-01,1.16025403784439E-01,1.00000000000000E+01/),&
       shape(c))
    real*8, parameter :: ll10 = c(1,nknot)

    if (t <= 0d0) then
       loglike_bayes = -t * t
    elseif (t <= 10d0) then
       loglike_bayes = splineval(nknot,c,t)
    else
       loglike_bayes = -log(t) + log(10d0) + ll10
    end if

  end function loglike_bayes

end submodule powderproc
