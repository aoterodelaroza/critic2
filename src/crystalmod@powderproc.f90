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

  !> Return the list of powder peaks between the initial (th2ini0) and
  !> final (th2end0) 2*theta (in angles), with wavelength lambda0 (in
  !> angstrom) and polarization factor fpol. If usehvecp, calculate
  !> only the reflections in p%hvec (requires p%npeaks and p%hvec). If
  !> calcderivs, calculate derivatives and set p%th2g and p%ipg. If gg
  !> is present, use gg instead of the metric tensor from c.
  module subroutine powder_peaks(c,p,th2ini0,th2end0,lambda0,fpol,usehvecp,calcderivs,gg)
    use tools_math, only: matinv
    use tools_io, only: ferror, faterr
    use tools, only: qcksort
    use types, only: realloc
    use param, only: pi, bohrtoa
    class(crystal), intent(in) :: c
    type(xrpd_peaklist), intent(inout) :: p
    real*8, intent(in) :: th2ini0, th2end0
    real*8, intent(in) :: lambda0
    real*8, intent(in) :: fpol
    logical, intent(in) :: usehvecp
    logical, intent(in) :: calcderivs
    real*8, intent(in), optional :: gg(3,3)

    integer :: i
    real*8 :: th2ini, th2end, lambda, hvec(3)
    real*8 :: smax, imax
    real*8 :: int, intg(6), th2, th2g(6)
    integer :: hmax, hcell, h, k, l
    logical :: again
    integer, allocatable :: io(:)
    real*8 :: grtensor(3,3), ar(3)

    integer, parameter :: mp = 200
    real*8, parameter :: ieps = 1d-5

    ! consistency checks
    if (usehvecp .and.(.not.p%havehvec.or..not.allocated(p%hvec).or.p%npeak<=0)) &
       call ferror('powder_peaks','requested usehvecp without hvec information',faterr)

    ! initialize angles
    th2ini = th2ini0 * pi / 180d0
    th2end = th2end0 * pi / 180d0

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
    p%havepeakshape = .false.
    p%havehvec = .true.
    p%havegradients = calcderivs
    if (allocated(p%th2)) deallocate(p%th2)
    if (allocated(p%ip)) deallocate(p%ip)
    if (allocated(p%th2g)) deallocate(p%th2g)
    if (allocated(p%ipg)) deallocate(p%ipg)
    if (allocated(p%fwhm)) deallocate(p%fwhm)
    if (allocated(p%cgau)) deallocate(p%cgau)
    if (usehvecp) then
       allocate(p%th2(p%npeak),p%ip(p%npeak))
       if (calcderivs) allocate(p%th2g(6,p%npeak),p%ipg(6,p%npeak))
    else
       if (allocated(p%hvec)) deallocate(p%hvec)
       allocate(p%th2(mp),p%ip(mp),p%hvec(3,mp))
       if (calcderivs) allocate(p%th2g(6,mp),p%ipg(6,mp))
    end if

    ! metric tensor, cell limits, convert lambda to bohr
    lambda = lambda0 / bohrtoa
    smax = sin((th2end)/2d0)
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
       do i = 1, p%npeak
          hvec = real(p%hvec(:,i),8)
          call run_function_body()
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

    ! normalize the intensities (highest peak is 1)
    imax = maxval(p%ip)
    p%ip = p%ip / imax
    if (calcderivs) p%ipg = p%ipg / imax

  contains
    subroutine run_function_body()
      use tools_io, only: ferror, faterr
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
         if (iz < 1 .or. iz > size(cscatt,2)) &
            call ferror('struct_powder','invalid Z -> no atomic scattering factors',faterr)
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

  end subroutine powder_peaks

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

  !> Read the XRPD peaks from a peaks file. The peaks are sorted by th2
  !> on output.
  module subroutine xrpd_peaks_from_file(p,file)
    use tools, only: qcksort
    use tools_io, only: fopen_read, getline, isreal, ferror, faterr, fclose
    use types, only: realloc
    type(xrpd_peaklist), intent(inout) :: p
    character*(*), intent(in) :: file

    integer :: i, lu, lp
    character(len=:), allocatable :: line
    real*8 :: r1, r2
    logical :: ok

    integer, parameter :: mpeak = 1000
    integer, allocatable :: io(:)

    ! initialize output
    p%npeak = 0
    if (allocated(p%th2)) deallocate(p%th2)
    if (allocated(p%ip)) deallocate(p%ip)
    if (allocated(p%hvec)) deallocate(p%hvec)
    if (allocated(p%th2g)) deallocate(p%th2g)
    if (allocated(p%ipg)) deallocate(p%ipg)
    if (allocated(p%fwhm)) deallocate(p%fwhm)
    if (allocated(p%cgau)) deallocate(p%cgau)
    p%havehvec = .false.
    p%havegradients = .false.
    p%havepeakshape = .true.

    ! read the pattern from the peaks file
    allocate(p%th2(mpeak),p%ip(mpeak),p%fwhm(mpeak),p%cgau(mpeak))
    lu = fopen_read(file)
    do while (getline(lu,line))
       lp = 1
       ok = isreal(r1,line,lp)
       ok = ok .and. isreal(r2,line,lp)
       if (.not.ok) &
          call ferror("xrpd_peaks_from_file","invalid input in peaks file",faterr)
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

  end subroutine xrpd_peaks_from_file

  ! Calculate the Gaussian cross-correlation between the peak patterns
  ! p1 and p2 with Gaussian triangle function width alpha and
  ! Gaussian width sigma, and output the value in d12. If calcderivs
  ! is present and true, calculate the derivatives wrt the metric
  ! tensor of p1 in d12g (requires p1%havegradients).
  module subroutine crosscorr_gaussian(p1,p2,alpha,sigma,d12,calcderivs,d12g)
    use tools_io, only: ferror, faterr
    use param, only: pi
    type(xrpd_peaklist), intent(in) :: p1, p2
    real*8, intent(in) :: alpha, sigma
    real*8, intent(out) :: d12
    logical, intent(in), optional :: calcderivs
    real*8, intent(out), optional :: d12g(6)

    logical :: calcg
    real*8 :: z, zp, z2p, zsq, thdeps, expt12, thdif
    integer :: imin, i, j

    real*8, parameter :: eps_discard = 1d-10

    ! consistency checks
    calcg = present(calcderivs)
    if (calcg) calcg = calcderivs
    if (calcg .and..not.p1%havegradients) &
       call ferror('crosscorr_gaussian','calcderivs requires p1%havegradients',faterr)
    if (calcg .and..not.present(d12g)) &
       call ferror('crosscorr_gaussian','calcderivs requires d12g being present',faterr)

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
  module subroutine gaussian_compare(c1,p2,imode,diff,seedout,verbose0,alpha0,&
     lambda0,fpol0,maxfeval0,besteps0,max_elong_def0,max_ang_def0)
    use crystalseedmod, only: crystalseed
    use tools_io, only: ferror
#ifndef HAVE_NLOPT
    use tools_io, only: warning
#else
    use tools_io, only: faterr, uout, string
    use tools_math, only: m_x2c_from_cellpar, det3
#endif
    type(crystal), intent(in) :: c1
    type(xrpd_peaklist), intent(in) :: p2
    integer, intent(in) :: imode
    real*8, intent(out) :: diff
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
    call ferror("gaussian_compare","gaussian_compare can only be used if nlopt is available",warning)
    diff = 1d0
    if (present(seedout)) call seedout%end()
    return
#else

    ! local block
    integer :: i
    real*8 :: dfg22, x(6), xorig(6), vorig, grad(6), x2c(3,3), omega
    integer*8 :: opt, lopt
    integer :: ires, imode, maxfeval
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
    real*8, parameter :: max_elong_defdef = 0.1d0
    real*8, parameter :: max_ang_defdef = 5d0
    integer, parameter :: maxfeval_def = 15000
    real*8, parameter :: besteps_def = 1d-4

    real*8, parameter :: sigma = 0.05d0
    real*8, parameter :: lambda_def = 1.5406d0
    real*8, parameter :: fpol_def = 0d0
    real*8, parameter :: alpha_def = 1.0d0 ! width of the "triangle"

    include 'nlopt.f'

    ! consistency checks
    if (p2%npeak == 0) &
       call ferror('gaussian_compare','no peaks found',faterr)
    if (imode /= imode_sp .and. imode /= imode_local .and.imode /= imode_global) &
       call ferror('gaussian_compare','imode must be one of 0, 1, 2',faterr)

    ! process the input options
    alpha = alpha_def
    if (present(alpha0)) alpha = alpha0
    lambda = lambda_def
    if (present(lambda0)) lambda = lambda0
    maxfeval = maxfeval_def
    if (present(maxfeval0)) maxfeval = maxfeval0
    besteps = besteps_def
    if (present(besteps0)) besteps = besteps0
    max_elong_def = max_elong_defdef
    if (present(max_elong_def0)) max_elong_def = max_elong_def0
    max_ang_def = max_ang_defdef
    if (present(max_ang_def0)) max_ang_def = max_ang_def0
    fpol = fpol_def
    if (present(fpol0)) fpol = fpol0
    verbose = .true.
    if (present(verbose0)) verbose = verbose0

    ! initialize
    neval = 0
    lastval = -1d0
    bestval = 1.1d0
    nbesteval = 0
    th2ini = p2%th2(1) - 2 * sigma
    th2end = p2%th2(p2%npeak) + 2 * sigma

    call c1%powder_peaks(p1,th2ini,th2end,lambda,fpol,.false.,.false.)
    call crosscorr_gaussian(p2,p2,alpha,sigma,dfg22,.false.)

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
      character(len=:), allocatable :: str

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
      call c1%powder_peaks(p1,th2ini,th2end,lambda,fpol,.true.,.true.,gg)
      call crosscorr_gaussian(p1,p1,alpha,sigma,dfg11,.true.,dfgg11)
      call crosscorr_gaussian(p1,p2,alpha,sigma,dfg12,.true.,dfgg12)

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

    end subroutine diff_fun
#endif
  end subroutine gaussian_compare

end submodule powderproc
