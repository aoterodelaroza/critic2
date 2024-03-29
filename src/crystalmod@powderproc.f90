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

  !> Calculate the radial distribution function.  On input, npts is
  !> the number of bins points from the initial (0) to the final
  !> (rend) distance and sigma is the Gaussian broadening of the
  !> peaks. If ishard=.false., represent the tails of the peaks
  !> outside the plot range.
  !> On output, t is the distance grid, and ih is the value of the
  !> RDF. This routine is based on:
  !>   Willighagen et al., Acta Cryst. B 61 (2005) 29.
  !> except using the sqrt of the atomic numbers instead of the
  !> charges.
  !>
  !> Optionally, if npairs0/ipairs0 are given, return the RDF of only
  !> the pairs of species given in the ipairs0 array. npairs0 is the
  !> number of selected pairs. If dp, ip, and isp are provided, on
  !> output they will contain the list of peaks that contribute to the
  !> RDF. Specifically, the distances (dp, in bohr), intensities (ip),
  !> and the atomic species (isp(2,*)).
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

    ! set up the chosen pairs
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
                if (ipairs(1,k) == kz) then
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

end submodule powderproc
