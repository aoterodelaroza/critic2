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

! Routines for the atomic environment class
submodule (environmod) proc
  implicit none

  integer, parameter :: menv0 = 10

  real*8 :: boxsize_default = 4.0 ! length of the region side (bohr)

  !xx! private procedures
  ! subroutine calculate_regions(e)
  ! function packoffset(i1,i2,i3,nr,nr21) result(ip)
  ! function unpackoffset(ip,nr,nr21) result(ioff)

contains
  
  !xx! environ class methods
  !> Free allocated arrays and nullify variables
  module subroutine environ_end(e)
    class(environ), intent(inout) :: e

    if (allocated(e%spc)) deallocate(e%spc)
    if (allocated(e%at)) deallocate(e%at)
    if (allocated(e%imap)) deallocate(e%imap)
    if (allocated(e%nrlo)) deallocate(e%nrlo)
    if (allocated(e%nrhi)) deallocate(e%nrhi)
    if (allocated(e%rs_ioffset)) deallocate(e%rs_ioffset)
    if (allocated(e%rs_rcut)) deallocate(e%rs_rcut)
    e%n = 0
    e%nspc = 0
    e%ncell = 0
    e%nreg = 0
    e%nmin = 0
    e%nmax = 0
    e%nregion = 0
    e%rs_nreg = 0

  end subroutine environ_end
  
  !> Build an environment from molecule data. nspc = number of species.
  !> spc = species. n = number of atoms in the cell. at = atoms in the
  !> cell. m_xr2c = reduced crystallographic to Cartesian matrix. 
  !> m_x2xr = crystallographic to reduced crystallographic matrix.
  module subroutine environ_build_from_molecule(e,nspc,spc,n,at,m_xr2c,m_x2xr,m_x2c)
    use tools_math, only: matinv
    use types, only: realloc, celatom
    use param, only: atmcov
    class(environ), intent(inout) :: e
    integer, intent(in) :: nspc
    type(species), intent(in) :: spc(nspc)
    integer, intent(in) :: n
    type(celatom), intent(in) :: at(n)
    real*8, intent(in) :: m_xr2c(3,3)
    real*8, intent(in) :: m_x2xr(3,3)
    real*8, intent(in) :: m_x2c(3,3)

    integer :: i 
    real*8 :: sphmax, rmin(3), rmax(3), dmax

    e%ismolecule = .true.
    e%nspc = nspc
    if (allocated(e%spc)) deallocate(e%spc)
    e%spc = spc

    ! calculate the boxsize, maximum bondfactor of 2.0
    dmax = 0d0
    do i = 1, nspc
       if (spc(i)%z > 0) then
          dmax = max(dmax,atmcov(spc(i)%z))
       end if
    end do
    e%boxsize = max(boxsize_default,4d0*dmax + 1d-10)

    ! fill the matrices
    e%m_xr2c = m_xr2c
    e%m_c2xr = matinv(m_xr2c)
    e%m_x2xr = m_x2xr
    e%m_xr2x = matinv(m_x2xr)
    e%m_x2c = m_x2c
    e%m_c2x = matinv(m_x2c)

    ! calculate the maximum diagonal half-length (sphmax)
    sphmax = norm2(e%xr2c((/0d0,0d0,0d0/) - (/0.5d0,0.5d0,0.5d0/)))
    sphmax = max(sphmax,norm2(e%xr2c((/1d0,0d0,0d0/) - (/0.5d0,0.5d0,0.5d0/))))
    sphmax = max(sphmax,norm2(e%xr2c((/0d0,1d0,0d0/) - (/0.5d0,0.5d0,0.5d0/))))
    sphmax = max(sphmax,norm2(e%xr2c((/0d0,0d0,1d0/) - (/0.5d0,0.5d0,0.5d0/))))
    e%rsph_uc = sphmax

    rmin = 1d40
    rmax = -1d40
    e%n = n
    if (allocated(e%at)) deallocate(e%at)
    allocate(e%at(e%n))
    do i = 1, n
       e%at(i)%x = at(i)%x
       e%at(i)%r = at(i)%r
       rmin = min(e%at(i)%r,rmin)
       rmax = max(e%at(i)%r,rmax)
       e%at(i)%idx = at(i)%idx
       e%at(i)%cidx = i
       e%at(i)%lvec = 0
       e%at(i)%is = at(i)%is
    end do
    e%ncell = n
    e%rsph_env = max(0.5d0 * norm2(rmax-rmin),e%rsph_uc)
    e%dmax0 = huge(1d0)

    call calculate_regions(e)

  end subroutine environ_build_from_molecule

  !> Build an environment from crystal data. nspc = number of species.
  !> spc = species. n = number of atoms in the cell. at = atoms in the
  !> cell. m_xr2c = reduced crystallographic to Cartesian matrix. 
  !> m_x2xr = crystallographic to reduced crystallographic matrix.
  !> dmax0 = the environment will contain all atoms within a distance 
  !> dmax0 of any point in the unit cell.
  module subroutine environ_build_from_crystal(e,nspc,spc,n,at,m_xr2c,m_x2xr,m_x2c,dmax0)
    use global, only: cutrad
    use tools_math, only: matinv
    use types, only: realloc, species, celatom
    use param, only: atmcov, bohrtoa
    class(environ), intent(inout) :: e
    integer, intent(in) :: nspc
    type(species), intent(in) :: spc(nspc)
    integer, intent(in) :: n
    type(celatom), intent(in) :: at(n)
    real*8, intent(in) :: m_xr2c(3,3)
    real*8, intent(in) :: m_x2xr(3,3)
    real*8, intent(in) :: m_x2c(3,3)
    real*8, intent(in), optional :: dmax0

    logical :: dorepeat
    real*8 :: sphmax, dmax, x(3), xhalf(3)
    integer :: i1, i2, i3, i, imax, i3min

    e%ismolecule = .false.
    e%nspc = nspc
    if (allocated(e%spc)) deallocate(e%spc)
    e%spc = spc

    ! calculate the boxsize
    dmax = 0d0
    do i = 1, nspc
       if (spc(i)%z > 0) then
          dmax = max(dmax,atmcov(spc(i)%z))
       end if
    end do
    e%boxsize = max(boxsize_default,4d0*dmax + 1d-10)

    ! fill the matrices
    e%m_xr2c = m_xr2c
    e%m_c2xr = matinv(m_xr2c)
    e%m_x2xr = m_x2xr
    e%m_xr2x = matinv(m_x2xr)
    e%m_x2c = m_x2c
    e%m_c2x = matinv(m_x2c)

    ! calculate the maximum diagonal half-length (sphmax)
    sphmax = norm2(e%xr2c((/0d0,0d0,0d0/) - (/0.5d0,0.5d0,0.5d0/)))
    sphmax = max(sphmax,norm2(e%xr2c((/1d0,0d0,0d0/) - (/0.5d0,0.5d0,0.5d0/))))
    sphmax = max(sphmax,norm2(e%xr2c((/0d0,1d0,0d0/) - (/0.5d0,0.5d0,0.5d0/))))
    sphmax = max(sphmax,norm2(e%xr2c((/0d0,0d0,1d0/) - (/0.5d0,0.5d0,0.5d0/))))
    e%rsph_uc = sphmax

    ! calculate the dmax if not given
    if (present(dmax0)) then
       dmax = dmax0
    else
       dmax = 0d0
       do i = 1, n
          if (spc(at(i)%is)%z > 0) dmax = max(dmax,cutrad(spc(at(i)%is)%z))
       end do
    end if
    e%dmax0 = dmax

    ! add all atoms in the main reduced cell (around 0,0,0)
    e%n = n
    if (allocated(e%at)) deallocate(e%at)
    allocate(e%at(e%n))
    do i = 1, n
       e%at(i)%x = e%x2xr(at(i)%x)
       e%at(i)%x = e%at(i)%x - floor(e%at(i)%x)
       x = e%xr2x(e%at(i)%x)
       e%at(i)%r = e%xr2c(e%at(i)%x)
       e%at(i)%idx = at(i)%idx
       e%at(i)%cidx = i
       e%at(i)%lvec = floor(x - at(i)%x)
       e%at(i)%is = at(i)%is
    end do
    e%ncell = n

    ! include all atoms at a distance sphmax+dmax from the center of the cell
    dorepeat = .true.
    imax = 0
    xhalf = 0.5d0
    xhalf = e%xr2c(xhalf)
    do while (dorepeat)
       dorepeat = .false.
       imax = imax + 1
       do i1 = -imax, imax
          do i2 = -imax, imax
             if (abs(i1) == imax .or. abs(i2) == imax) then
                i3min = 0
             else
                i3min = imax
             end if
             do i3 = i3min, imax
                call addcell(i1,i2,i3)
                if (i3 /= 0) &
                   call addcell(i1,i2,-i3)
             end do
          end do
       end do
    end do
    call realloc(e%at,e%n)
    e%rsph_env = sphmax+dmax

    call calculate_regions(e)
    
  contains
    subroutine addcell(i1,i2,i3)
      integer, intent(in) :: i1, i2, i3

      integer :: p(3), px(3)
      real*8 :: xc(3), dd

      do i = 1, n
         p = (/i1, i2, i3/)
         px = nint(e%xr2x(real(p,8)))

         x = e%at(i)%x + p
         xc = e%xr2c(x)
         dd = norm2(xc - xhalf)
         if (dd < sphmax+dmax) then
            dorepeat = .true.
            e%n = e%n + 1
            if (e%n > size(e%at,1)) call realloc(e%at,2*e%n)

            e%at(e%n)%x = x
            e%at(e%n)%r = xc
            e%at(e%n)%idx = e%at(i)%idx
            e%at(e%n)%cidx = e%at(i)%cidx
            e%at(e%n)%lvec = e%at(i)%lvec + px
            e%at(e%n)%is = e%at(i)%is
         end if
      end do

    end subroutine addcell
  end subroutine environ_build_from_crystal

  !> Reduced crystallographic to Cartesian transform
  pure module function xr2c(e,xx) result(res)
    class(environ), intent(in) :: e
    real*8, intent(in)  :: xx(3)
    real*8 :: res(3)
    res = matmul(e%m_xr2c,xx)
  end function xr2c

  !> Cartesian to reduced crystallographic transform
  pure module function c2xr(e,xx) result(res)
    class(environ), intent(in) :: e
    real*8, intent(in)  :: xx(3)
    real*8 :: res(3)
    res = matmul(e%m_c2xr,xx)
  end function c2xr

  !> Reduced crystallographic to crystallographic transform
  pure module function xr2x(e,xx) result(res)
    class(environ), intent(in) :: e
    real*8, intent(in)  :: xx(3)
    real*8 :: res(3)
    res = matmul(e%m_xr2x,xx)
  end function xr2x

  !> Crystallographic to reduced crystallographic transform
  pure module function x2xr(e,xx) result(res)
    class(environ), intent(in) :: e
    real*8, intent(in)  :: xx(3)
    real*8 :: res(3)
    res = matmul(e%m_x2xr,xx)
  end function x2xr

  !> Cartesian to crystallographic transform
  pure module function c2x(e,xx) result(res)
    class(environ), intent(in) :: e
    real*8, intent(in)  :: xx(3)
    real*8 :: res(3)
    res = matmul(e%m_c2x,xx)
  end function c2x

  !> Crystallographic to Cartesian transform
  pure module function x2c(e,xx) result(res)
    class(environ), intent(in) :: e
    real*8, intent(in)  :: xx(3)
    real*8 :: res(3)
    res = matmul(e%m_x2c,xx)
  end function x2c

  !> Convert between coordinate type icrd (cartesian, cryst., reduced
  !> cryst.) and type ocrd. x2c and c2x not available.
  pure module function y2z(e,xx,icrd,ocrd) result(res)
    use param, only: icrd_cart, icrd_crys, icrd_rcrys
    class(environ), intent(in) :: e
    real*8, intent(in)  :: xx(3)
    integer, intent(in) :: icrd, ocrd
    real*8 :: res(3)

    if (icrd == icrd_rcrys .and. ocrd == icrd_cart) then
       res = matmul(e%m_xr2c,xx)
    else if (icrd == icrd_cart .and. ocrd == icrd_rcrys) then
       res = matmul(e%m_c2xr,xx)
    else if (icrd == icrd_rcrys .and. ocrd == icrd_crys) then
       res = matmul(e%m_xr2x,xx)
    else if (icrd == icrd_crys .and. ocrd == icrd_rcrys) then
       res = matmul(e%m_x2xr,xx)
    else if (icrd == icrd_cart .and. ocrd == icrd_crys) then
       res = matmul(e%m_c2x,xx)
    else if (icrd == icrd_crys .and. ocrd == icrd_cart) then
       res = matmul(e%m_x2c,xx)
    else
       res = xx
    end if

  end function y2z

  !> Given a point xx in icrd coordinates, transform to ocrd
  !> coordinates. If the environment belongs to a crystal, also
  !> translate to the main cell. If lvec is present, returns the
  !> lattice translation in reduced crystallographic coordinates.
  pure module subroutine y2z_center(e,xx,icrd,ocrd,lvec)
    use param, only: icrd_rcrys
    class(environ), intent(in) :: e
    real*8, intent(inout) :: xx(3)
    integer, intent(in) :: icrd
    integer, intent(in) :: ocrd
    integer, intent(out), optional :: lvec(3)

    xx = e%y2z(xx,icrd,icrd_rcrys)
    if (present(lvec)) lvec = floor(xx)
    xx = xx - floor(xx)
    xx = e%y2z(xx,icrd_rcrys,ocrd)

  end subroutine y2z_center

  !> Cartesian to region transform
  pure module function c2p(e,xx) result(res)
    class(environ), intent(in) :: e
    real*8, intent(in)  :: xx(3)
    integer :: res(3)

    res = floor((xx - e%x0) / e%boxsize)

  end function c2p

  !> Region to integer region transform
  pure module function p2i(e,xx) result(res)
    class(environ), intent(in) :: e
    integer, intent(in)  :: xx(3)
    integer :: res

    integer :: ix(3)

    ix = min(max(xx,e%nmin),e%nmax) - e%nmin
    res = 1 + ix(1) + e%nreg(1) * (ix(2) + e%nreg(2) * ix(3))

  end function p2i

  !> Cartesian to integer region transform
  pure module function c2i(e,xx) result(res)
    class(environ), intent(in) :: e
    real*8, intent(in)  :: xx(3)
    integer :: res

    integer :: ix(3)

    ix = min(max(floor((xx - e%x0) / e%boxsize),e%nmin),e%nmax) - e%nmin
    res = 1 + ix(1) + e%nreg(1) * (ix(2) + e%nreg(2) * ix(3))

  end function c2i

  !> Translate point x0 (with icrd input coordinates) to the main cell
  !> and, if it corresponds to an atomic position (to within atomeps),
  !> return the ID of the atom. Otherwise, return 0. If lncel is
  !> .false. or not present, the ID is for the non-equivalent atom
  !> list. Otherwise, it is for the complete list. This routine is
  !> thread-safe.
  module function identify_atom(e,x0,icrd,lncel)
    use global, only: atomeps
    use param, only: icrd_cart
    class(environ), intent(in) :: e
    real*8, intent(in) :: x0(3)
    integer, intent(in) :: icrd
    logical, intent(in), optional :: lncel
    integer :: identify_atom
    
    real*8 :: xp(3), x0r(3), x1r(3), dist, distmin
    logical :: ln
    integer :: ireg0(3), idx0, imin(3), imax(3), i, j, k, i1, i2, i3
    integer :: ireg(3), idx, kmin

    ! initialize
    ln = .false.
    if (present(lncel)) ln = lncel

    ! bring the atom to the main cell, convert, identify the region
    xp = x0
    call e%y2z_center(xp,icrd,icrd_cart)
    ireg0 = e%c2p(x0)
    idx0 = e%p2i(ireg0)

    ! calculate the regions to explore
    imin = 0
    imax = 0
    x0r = e%x0 + e%boxsize * ireg0
    x1r = e%x0 + e%boxsize * (ireg0+1)
    do i = 1, 3
       if (xp(i)-x0r(i) <= atomeps) imin(i) = -1
       if (x1r(i)-xp(i) <= atomeps) imax(i) = 1
    end do

    ! Identify the atom
    kmin = 0
    distmin = 1d40
    do i1 = imin(1), imax(1)
       do i2 = imin(2), imax(2)
          do i3 = imin(3), imax(3)
             ireg = ireg0 + (/i1,i2,i3/)
             if (any(ireg < e%nmin) .or. any(ireg > e%nmax)) cycle
             idx = e%p2i(ireg)
             if (e%nrhi(idx) == 0) cycle
             do j = e%nrlo(idx), e%nrhi(idx)
                k = e%imap(j)
                dist = norm2(e%at(k)%r - xp)
                if (dist < distmin .and. dist < atomeps) then
                   kmin = k
                end if
             end do
          end do
       end do
    end do

    ! return the atom index
    identify_atom = 0
    if (kmin > 0) then
       if (ln) then
          identify_atom = e%at(kmin)%cidx
       else
          identify_atom = e%at(kmin)%idx
       end if
    end if

  end function identify_atom

  !> Given the point xp (in icrd coordinates), translates to the main
  !> cell and calculates the nearest atom.  The nearest atom has ID
  !> nid from the complete list (atcel) and is at a distance dist, or
  !> nid=0 and dist=0d0 if the search did not produce any atoms.  On
  !> output, the optional argument lvec contains the lattice vector to
  !> the nearest atom (i.e. its position is atcel(nid)%x + lvec). If
  !> nid0, consider only atoms with index nid0 from the non-equivalent
  !> list. If id0, consider only atoms with index id0 from the
  !> complete list. If nozero, disregard zero-distance atoms. This
  !> routine is thread-safe.
  module subroutine nearest_atom(e,xp,icrd,nid,dist,lvec,nid0,id0,nozero)
    use param, only: icrd_cart
    class(environ), intent(in) :: e
    real*8, intent(in) :: xp(3)
    integer, intent(in) :: icrd
    integer, intent(out) :: nid
    real*8, intent(out) :: dist
    integer, intent(out), optional :: lvec(3)
    integer, intent(in), optional :: nid0
    integer, intent(in), optional :: id0
    logical, intent(in), optional :: nozero

    real*8, parameter :: eps = 1d-10

    real*8 :: x0(3), dist, distmin
    integer :: ireg0(3), ireg(3), idxreg
    integer :: i, j, k, kmin, lvec0(3)

    ! Find the integer region for the main cell copy of the input point
    x0 = xp
    call e%y2z_center(x0,icrd,icrd_cart,lvec0)
    ireg0 = e%c2p(x0)
    
    ! run over regions sorted by distance
    distmin = 1d40
    kmin = 0
    main: do i = 1, e%rs_nreg
       ireg = ireg0 + unpackoffset(e%rs_ioffset(i),e%rs_imax,e%rs_2imax1)
       if (any(ireg < e%nmin) .or. any(ireg > e%nmax)) cycle
       idxreg = e%p2i(ireg)
       if (e%nrhi(idxreg) == 0) cycle

       do j = e%nrlo(idxreg), e%nrhi(idxreg)
          k = e%imap(j)
          if (distmin < e%rs_rcut(i)) exit main
          if (present(nid0)) then
             if (e%at(k)%idx /= nid0) cycle
          end if
          if (present(id0)) then
             if (e%at(k)%cidx /= id0) cycle
          end if

          dist = norm2(e%at(k)%r - x0)
          if (present(nozero)) then
             if (dist < eps) cycle
          end if
          if (dist < distmin) then
             distmin = dist
             kmin = k
          end if
       end do
    end do main

    if (kmin > 0) then
       nid = e%at(kmin)%cidx
       dist = distmin
       if (present(lvec)) then
          lvec = e%at(kmin)%lvec + nint(e%xr2x(real(lvec0,8)))
       end if
    else
       nid = 0
       dist = 0d0
       if (present(lvec)) lvec = 0
    end if

  end subroutine nearest_atom

  !> Given the point xp (in icrd coordinates), center the point in the
  !> main cell and calculate the list of nearest atoms. If sorted is
  !> true, the list is sorted by distance and shell (if up2n or up2sh
  !> is used, the list is always sorted). The output list eid contains
  !> nat atoms with IDs nid(1:nat) from the environment, distances to
  !> the input point equal to dist(1:nat) and lattice vector lvec in
  !> cryst. coords. The position of atom i in cryst. coords. is
  !> e%at(nid(i))%x + lvec. Optionally, ishell0(i) contains the shell
  !> ID for atom i in the output list. One or more of four cutoff
  !> criteria must be chosen: list up to a distance up2d, between a
  !> minimum and a maximum species-dependent distance (up2dsp), up to
  !> a number of shells up2sh or up to a number of atoms up2n. If
  !> nid0, consider only atoms with index nid0 from the non-equivalent
  !> list. If id0, consider only atoms with index id0 from the
  !> complete list. If nozero, disregard zero-distance atoms. The
  !> output error condition ierr is 0 if the search was successful or
  !> non-zero if the input search conditions could not be met by this
  !> environment.
  module subroutine list_near_atoms(e,xp,icrd,sorted,nat,eid,dist,lvec,ierr,ishell0,up2d,up2dsp,up2sh,up2n,nid0,id0,nozero)
    use global, only: atomeps
    use tools_io, only: ferror, faterr
    use tools, only: mergesort
    use types, only: realloc
    use param, only: icrd_cart, ctsq32, ctsq3
    class(environ), intent(in) :: e
    real*8, intent(in) :: xp(3)
    integer, intent(in) :: icrd
    logical, intent(in) :: sorted
    integer, intent(out) :: nat
    integer, allocatable, intent(inout) :: eid(:)
    real*8, allocatable, intent(inout) :: dist(:)
    integer, intent(out) :: lvec(3)
    integer, intent(out) :: ierr
    integer, allocatable, intent(inout), optional :: ishell0(:)
    real*8, intent(in), optional :: up2d
    real*8, intent(in), optional :: up2dsp(:,:)
    integer, intent(in), optional :: up2sh
    integer, intent(in), optional :: up2n
    integer, intent(in), optional :: nid0
    integer, intent(in), optional :: id0
    logical, intent(in), optional :: nozero

    real*8, parameter :: eps = 1d-10

    real*8 :: x0(3), dist0, rcutshel, rcutn, up2rmax, rext, xhalf(3), dmaxenv
    integer :: ireg0(3), ireg(3), idxreg, nats
    integer :: i, j, k, imax, i1, i2, i3, ii1, ii2, ii3, isign(3), i0loop
    integer, allocatable :: iord(:), iiord(:), ishell(:)
    integer :: nshel
    real*8, allocatable :: rshel(:)
    integer, allocatable :: idxshel(:)
    logical :: doshell, enough

    if (.not.present(up2d).and..not.present(up2dsp).and..not.present(up2sh).and..not.present(up2n)) &
       call ferror("list_near_atoms","must give one of up2d, up2dsp, up2sh, or up2n",faterr)
    doshell = present(ishell0) .or. present(up2sh)

    ! Find the integer region for the main cell copy of the input point
    x0 = xp
    call e%y2z_center(x0,icrd,icrd_cart,lvec)
    ireg0 = e%c2p(x0)
    lvec = nint(e%xr2x(real(lvec,8)))

    ! calculate the maximum distance for which we guarantee we will
    ! get all the atoms
    xhalf = 0.5d0
    xhalf = e%xr2c(xhalf)
    dmaxenv = max(e%rsph_env - norm2(x0 - xhalf),0d0)

    ! Initialize the output and auxiliary arrays
    nat = 0
    nshel = 0
    rcutshel = -1d0
    rcutn = -1d0
    if (allocated(eid)) deallocate(eid)
    if (allocated(dist)) deallocate(dist)
    allocate(eid(10),dist(10))
    if (doshell) then
       allocate(ishell(10),rshel(10),idxshel(10))
    end if
    if (present(up2dsp)) up2rmax = maxval(up2dsp(1:e%nspc,2))

    ! Run over search regions around ireg0 and over atoms belonging to
    ! those regions. Use the sorted search data in the rs_* variables.
    enough = .false.
    i = 0
    nats = 0
    do while (.not.enough .and. i < e%rs_nreg)
       i = i + 1
       ireg = ireg0 + unpackoffset(e%rs_ioffset(i),e%rs_imax,e%rs_2imax1)

       if (all(ireg >= e%nmin) .and. all(ireg <= e%nmax)) then
          idxreg = e%p2i(ireg)
          do j = e%nrlo(idxreg), e%nrhi(idxreg)
             k = e%imap(j)
             nats = nats + 1

             ! apply nid0 and id0 conditions
             if (present(nid0)) then
                if (e%at(k)%idx /= nid0) cycle
             end if
             if (present(id0)) then
                if (e%at(k)%cidx /= id0) cycle
             end if

             ! calculate the distance to this atom
             dist0 = norm2(e%at(k)%r - x0)

             ! apply the nozero and up2x conditions
             if (present(nozero)) then
                if (dist0 < eps) cycle
             end if
             if (present(up2d)) then
                if (dist0 > up2d) cycle
             end if
             if (present(up2dsp)) then
                if (dist0 < up2dsp(e%at(k)%is,1) .or. dist0 > up2dsp(e%at(k)%is,2)) cycle
             end if

             ! Add this atom to the list (contains section subroutine) and update rcutn.
             call add_atom_to_output_list()

             ! Add this shell to the list (contains section subroutine) and update rcutshel.
             call add_shell_to_output_list()
          end do
       end if

       ! We have enough if...
       ! - this is a molecule and we went over all the atoms
       enough = (e%n == nats .and. e%ismolecule)
       if (present(up2d)) then
          ! ...or we captured all the atoms up to distance up2d, which happens if:
          ! - Current region has all its points at a distance > up2d from any point in the reference cell (e%rs_rcut(i) > up2d).
          ! - The regions not considered by our offsets all have rs_rcut(i) > rsph_env > dmaxenv > up2d or there are no
          !   regions not considered by our offsets because this is a molecule
          !   (note: a molecule has no atoms outside the search offsets)
          enough = enough .or. (e%rs_rcut(i) > up2d .and. (dmaxenv > up2d .or. e%ismolecule))
       end if
       if (present(up2dsp)) then
          ! Same as above, but with maxval(up2dsp).
          enough = enough .or. (e%rs_rcut(i) > up2rmax .and. (dmaxenv > up2rmax .or. e%ismolecule))
       end if
       if (present(up2sh)) then
          ! ...or we captured all the atoms up to shell up2sh if:
          ! - We know at least up2sh shells and the farthest known shell is at rcutshel distance
          ! - Current region has all its points at a distance > farthest shell (rs_rcut(i) > rcutshel > rcutshel(ordered list))
          ! - The external regions all have rs_rcut(i) > rsph_env > dmaxenv > rcutshel > rcutshel(ordered list)
          !   (note: a molecule has no atoms outside the search offsets)
          enough = enough .or. (nshel >= up2sh .and. e%rs_rcut(i) > rcutshel .and. (dmaxenv > rcutshel .or. e%ismolecule))
       end if
       if (present(up2n)) then
          ! ...or we captured all the atoms up to number up2n if:
          ! - We know at least up2n atoms, and the farthest known atom is at rcutn distance
          ! - Current region has all its points at a distance > farthest atom (rs_rcut(i) > rcutn > rcutn(ordered list))
          ! - The external regions all have rs_rcut(i) > rsph_env > dmaxenv > rcutn > rcutn(ordered list)
          !   (note: a molecule has no atoms outside the search offsets)
          enough = enough .or. (nat >= up2n .and. e%rs_rcut(i) > rcutn .and. (dmaxenv > rcutn .or. e%ismolecule))
       end if
    end do

    ! Rearrange the arrays 
    if (nat > 0) then
       ! First reallocation
       call realloc(eid,nat)
       call realloc(dist,nat)
       if (doshell) call realloc(ishell,nat)

       ! Re-order shell labels first by distance, then by atom ID
       if (doshell) then
          allocate(iord(nshel),iiord(nshel))
          do i = 1, nshel
             iord(i) = i
          end do
          call mergesort(idxshel,iord,1,nshel)
          call mergesort(rshel,iord,1,nshel)
          do i = 1, nshel
             iiord(iord(i)) = i
          end do
          deallocate(iord)
          do i = 1, nat
             ishell(i) = iiord(ishell(i))
          end do
          deallocate(iiord)
       end if

       ! Sort output atoms by distance
       if (sorted .or. present(up2sh) .or. present(up2n)) then
          ! First by shell (if available) then by distance
          allocate(iord(nat))
          do i = 1, nat
             iord(i) = i
          end do
          if (doshell) call mergesort(ishell,iord,1,nat)
          call mergesort(dist,iord,1,nat)
          eid = eid(iord)
          dist = dist(iord)
          if (doshell) ishell = ishell(iord)
          deallocate(iord)

          ! Prune the extra atoms and reallocate
          if (present(up2sh) .or. present(up2n)) then
             if (present(up2sh)) then
                do i = 1, nat
                   if (ishell(i) > up2sh) then
                      nat = i-1
                      exit
                   end if
                end do
             end if
             if (present(up2n)) nat = min(up2n,nat)
             call realloc(eid,nat)
             call realloc(dist,nat)
             call realloc(ishell,nat)
          end if
       end if
    end if
    if (allocated(rshel)) deallocate(rshel)
    if (allocated(idxshel)) deallocate(idxshel)
    if (present(ishell0)) then
       if (allocated(ishell0)) deallocate(ishell0)
       call move_alloc(ishell,ishell0)
    end if

    ! Write the output error condition
    ierr = 0
    if (.not.enough) ierr = 1

  contains
    subroutine add_atom_to_output_list()
      nat = nat + 1
      if (nat > size(eid,1)) then
         call realloc(eid,2*nat)
         call realloc(dist,2*nat)
         if (doshell) call realloc(ishell,2*nat)
      end if
      eid(nat) = k
      dist(nat) = dist0

      ! Update the up2n rcut, distance to farthest known atom in the initial 1->up2n list.  
      ! The final (ordered) list will have its farthest atom at a distance less than rcutn.
      if (present(up2n)) then
         if (nat <= up2n) rcutn = max(dist0,rcutn)
      end if
    end subroutine add_atom_to_output_list
    
    subroutine add_shell_to_output_list()
      use tools_io, only: string
      integer :: lthis, l

      if (.not.doshell) return

      ! See if this shell is already known. (Make a better hash to improve this.)
      lthis = 0
      do l = 1, nshel
         if ((abs(dist0 - rshel(l)) < atomeps) .and. (idxshel(l) == e%at(k)%idx)) then
            lthis = l
            exit
         end if
      end do

      ! If not known, create the new shell
      if (lthis == 0) then
         nshel = nshel + 1
         if (nshel > size(rshel,1)) then
            call realloc(rshel,2*nshel)
            call realloc(idxshel,2*nshel)
         end if
         lthis = nshel
         idxshel(nshel) = e%at(k)%idx
         rshel(nshel) = dist0
      end if

      ! Update the up2sh rcut, distance to farthest known shell in the initial 1->up2sh list.
      ! The final (ordered) list will have its farthest shell at a distance less than rcutshel.
      if (present(up2sh)) then
         if (nshel <= up2sh) rcutshel = max(dist0,rcutshel)
      end if

      ! Write down the assigned shell for this atom
      ishell(nat) = lthis

    end subroutine add_shell_to_output_list

  end subroutine list_near_atoms
    
  !> Translate the point x0 to the main cell and calculate the core
  !> (if zpsp is present) or promolecular densities at a point x0
  !> (coord format given by icrd) using atomic radial grids up to a
  !> number of derivatives nder (max: 2). Returns the density (f),
  !> gradient (fp, nder >= 1), and Hessian (fpp, nder >= 2). If a
  !> fragment (fr) is given, then only the atoms in it
  !> contribute. This routine is thread-safe.
  module subroutine promolecular(e,x0,icrd,f,fp,fpp,nder,zpsp,fr)
    use grid1mod, only: cgrid, agrid, grid1
    use global, only: cutrad
    use fragmentmod, only: fragment
    use tools_io, only: ferror, faterr
    use param, only: icrd_cart, maxzat
    class(environ), intent(in) :: e
    real*8, intent(in) :: x0(3) !< Point in cryst. coords.
    integer, intent(in) :: icrd !< Input coordinates
    real*8, intent(out) :: f !< Density
    real*8, intent(out) :: fp(3) !< Density gradient
    real*8, intent(out) :: fpp(3,3) !< Density hessian
    integer, intent(in) :: nder !< Number of derivatives to calculate
    integer, intent(in), optional :: zpsp(:) !< core charges
    type(fragment), intent(in), optional :: fr !< Fragment contributing to the density

    integer :: i, j, k, ii, iz, ierr
    real*8 :: xc(3), xx(3), rlvec(3), r, rinv1, rinv1rp
    real*8 :: rho, rhop, rhopp, rfac, radd
    logical :: iscore
    type(grid1), pointer :: g
    integer :: nat, lvec(3)
    integer, allocatable :: nid(:)
    real*8, allocatable :: dist(:), rcutmax(:,:)
    logical, allocatable :: isinfr(:)

    f = 0d0
    fp = 0d0
    fpp = 0d0
    iscore = present(zpsp)
    if (iscore) then
       if (all(zpsp <= 0)) return
    end if
    if (iscore.and..not.allocated(cgrid)) then
       call ferror("promolecular","cgrid not allocated",faterr)
    elseif (.not.iscore.and..not.allocated(agrid)) then
       call ferror("promolecular","agrid not allocated",faterr)
    end if

    ! convert to Cartesian and move to the main cell
    xc = x0
    call e%y2z_center(xc,icrd,icrd_cart)

    ! compute the list of atoms that contribute to the point
    allocate(rcutmax(e%nspc,2))
    rcutmax = 0d0
    do i = 1, e%nspc
       iz = e%spc(i)%z
       if (iz == 0 .or. iz > maxzat) cycle
       if (iscore) then
          g => cgrid(iz,zpsp(iz))
       else
          g => agrid(iz)
       end if
       rcutmax(i,2) = min(cutrad(iz),g%rmax)
    end do
    call e%list_near_atoms(xc,icrd_cart,.false.,nat,nid,dist,lvec,ierr,up2dsp=rcutmax)
    rlvec = lvec
    rlvec = e%x2c(rlvec)
    deallocate(rcutmax)

    ! if fragment is provided, use only the atoms in it (by their cell index)
    if (present(fr)) then
       allocate(isinfr(e%ncell))
       isinfr = .false.
       do i = 1, fr%nat
          isinfr(fr%at(i)%cidx) = .true.
       end do
    end if

    ! Do the density and derivatives sum
    do ii = 1, nat
       i = nid(ii)
       if (present(fr)) then
          if (.not.isinfr(e%at(i)%cidx)) cycle
       end if

       iz = e%spc(e%at(i)%is)%z
       if (iz == 0 .or. iz > maxzat) cycle
       r = dist(ii)
       if (r > cutrad(iz)) cycle

       if (iscore) then
          g => cgrid(iz,zpsp(iz))
       else
          g => agrid(iz)
       end if
       r = max(max(r,g%r(1)),1d-14)
       call g%interp(r,rho,rhop,rhopp)
       rho = max(rho,0d0)
       f = f + rho

       if (nder < 1) cycle
       xx = xc - (e%at(i)%r + rlvec)
       rinv1 = 1d0 / r
       rinv1rp = rinv1 * rhop
       fp = fp + xx * rinv1rp

       if (nder < 2) cycle
       rfac = (rhopp - rinv1rp) * rinv1 * rinv1
       do j = 1, 3
          fpp(j,j) = fpp(j,j) + rinv1rp + rfac * xx(j) * xx(j)
          do k = 1, j-1
             radd = rfac * xx(j) * xx(k)
             fpp(j,k) = fpp(j,k) + radd
             fpp(k,j) = fpp(k,j) + radd
          end do
       end do
    end do

    if (allocated(isinfr)) deallocate(isinfr)

  end subroutine promolecular

  !> Find the covalent bond connectivity and return the bonds in the
  !> nstar array. To use this routine it is necessary that the boxsize
  !> is larger than the maximum covalent distance plus the tolerance.
  subroutine find_asterisms_covalent(e,nstar)
    use global, only: atomeps, bondfactor
    use tools_io, only: ferror, faterr
    use types, only: realloc
    use param, only: icrd_cart, atmcov
    class(environ), intent(in) :: e
    type(neighstar), allocatable, intent(inout) :: nstar(:)
    
    integer :: i, j, imin(3), imax(3)
    real*8 :: xmin(3), xmax(3), x0(3), dist2, ri, rj, rij2, r2
    integer :: i1, i2, i3, p0(3), p1(3), idx0, idx1
    integer :: j1, j2, j3, ki, kj, iz, is, js
    real*8, allocatable :: rij2(:,:,:)

    ! allocate the asterism arrays
    if (allocated(nstar)) deallocate(nstar)
    allocate(nstar(e%ncell))
    do i = 1, e%ncell
       nstar(i)%ncon = 0
       if (allocated(nstar(i)%idcon)) deallocate(nstar(i)%idcon)
       if (allocated(nstar(i)%lcon)) deallocate(nstar(i)%lcon)
       allocate(nstar(i)%idcon(20))
       allocate(nstar(i)%lcon(3,20))
    end do
    
    ! pre-calculate the distance^2 matrix
    allocate(rij2(e%nspc,e%nspc,2))
    rij2 = 0d0
    do i = 1, e%nspc
       ri = atmcov(e%spc(i)%z)
       do j = 1, e%nspc
          rj = atmcov(e%spc(j)%z)
          
          r2 = (ri+rj) * bondfactor
          rij2(i,j,2) = r2*r2
          rij2(j,i,2) = rij2(i,j,2)

          r2 = ri+rj / bondfactor
          rij2(i,j,1) = r2*r2
          rij2(j,i,1) = rij2(i,j,1)
       end do
    end do
    if (sqrt(maxval(rij2)) > e%boxsize) &
       call ferror("find_asterisms_covalent","boxsize too small for find_asterisms_covalent",faterr)

    ! find the first and last region that cover the unit cell
    xmin = 1d40
    xmax = -1d40
    do i = 1, e%ncell
       xmin = min(xmin,e%at(i)%r)
       xmax = max(xmax,e%at(i)%r)
    end do
    imin = e%c2p(xmin)
    imax = e%c2p(xmax)

    ! run over atoms in the unit cell and build the connectivity star
    do ki = 1, e%ncell
       p0 = e%c2p(e%at(ki)%r)
       idx0 = e%p2i(p0)
       is = e%at(ki)%is
       if (e%spc(is)%z == 0) cycle

       do j1 = -1, 1
          do j2 = -1, 1
             do j3 = -1, 1
                p1 = p0 + (/j1,j2,j3/)
                if (any(p1 < e%nmin .or. p1 > e%nmax)) cycle
                idx1 = e%p2i(p1)

                do j = e%nrlo(idx1), e%nrhi(idx1)
                   kj = e%imap(j)
                   js = e%at(kj)%is
                   if (e%spc(js)%z == 0) cycle

                   x0 = e%at(ki)%r - e%at(kj)%r
                   dist2 = x0(1)*x0(1)+x0(2)*x0(2)+x0(3)*x0(3)
                   if (dist2 >= rij2(is,js,1) .and. dist2 <= rij2(is,js,2)) then
                      nstar(ki)%ncon = nstar(ki)%ncon + 1
                      if (nstar(ki)%ncon > size(nstar(ki)%idcon,1)) then
                         call realloc(nstar(ki)%idcon,2*nstar(ki)%ncon)
                         call realloc(nstar(ki)%lcon,3,2*nstar(ki)%ncon)
                      end if
                      nstar(ki)%idcon(nstar(ki)%ncon) = e%at(kj)%cidx
                      nstar(ki)%lcon(:,nstar(ki)%ncon) = e%at(kj)%lvec
                   end if
                end do

             end do
          end do
       end do
    end do

    ! reallocate the asterism arrays
    do i = 1, e%ncell
       if (nstar(i)%ncon > 0) then
          call realloc(nstar(i)%idcon,nstar(i)%ncon)
          call realloc(nstar(i)%lcon,3,nstar(i)%ncon)
       end if
    end do

  end subroutine find_asterisms_covalent

  !> Find the bond connectivity using the radii in rtable (bohr) and a
  !> tolerance factor. rtable is indexed with the atomic number.
  !> Return the bonds in the nstar array. This routine can be used
  !> regardless of boxsize but is slower than find_asterisms_covalent.
  subroutine find_asterisms_listatoms(e,nstar,rtable,factor)
    use param, only: icrd_cart
    class(environ), intent(in) :: e
    type(neighstar), allocatable, intent(inout) :: nstar(:)
    real*8, intent(in) :: rtable(:)
    real*8, intent(in) :: factor
    
    integer :: i, j, iz, ierr, lvec(3), nat
    real*8, allocatable :: d0sp(:), dsp(:,:), dist(:)
    integer, allocatable :: eid(:)

    if (allocated(nstar)) deallocate(nstar)
    allocate(nstar(e%ncell))
    
    ! build the table of radii for the species
    allocate(d0sp(e%nspc),dsp(e%nspc,2))
    d0sp = 0d0
    do i = 1, e%nspc
       iz = e%spc(i)%z
       if (iz > 0) then
          d0sp(i) = rtable(iz)
       else
          d0sp(i) = -1d40
       end if
    end do

    ! list the neighbors for each atom in the unit cell and build the stars
    do i = 1, e%ncell
       iz = e%spc(e%at(i)%is)%z
       nstar(i)%ncon = 0
       if (iz > 0) then
          dsp(:,1) = (d0sp + rtable(iz)) * factor
          dsp(:,2) = (d0sp + rtable(iz)) / factor
          call e%list_near_atoms(e%at(i)%r,icrd_cart,.false.,nat,eid,dist,lvec,ierr,up2dsp=dsp,nozero=.true.)

          if (allocated(nstar(i)%idcon)) deallocate(nstar(i)%idcon)
          if (allocated(nstar(i)%lcon)) deallocate(nstar(i)%lcon)

          if (nat > 0) then
             allocate(nstar(i)%idcon(nat),nstar(i)%lcon(3,nat))
             nstar(i)%ncon = nat
             do j = 1, nat
                nstar(i)%idcon(j) = e%at(eid(j))%cidx
                nstar(i)%lcon(:,j) = e%at(eid(j))%lvec + lvec
             end do
          end if
       end if
    end do

  end subroutine find_asterisms_listatoms

  !> Write a report about the environment to stdout
  module subroutine environ_report(e)
    use global, only: iunitname0, dunit0, iunit
    use tools_io, only: uout, string
    class(environ), intent(in) :: e

    integer :: j
    integer :: idmax, idavg

    write (uout,'("+ Atomic environment")')
    write (uout,'("  Number of atoms (reduced cell/environment): ",A," / ",A)') string(e%ncell), string(e%n)
    write (uout,'("  Radius of (unit cell/environment) circumscribed sphere (",A,"): ",A," / ",A)') &
       iunitname0(iunit), trim(string(e%rsph_uc*dunit0(iunit),'f',8,4)), trim(string(e%rsph_env*dunit0(iunit),'f',8,4))
    if (e%dmax0 == huge(1d0)) then
       write (uout,'("  Maximum interaction distance: all atoms in the environment.")')
    else
       write (uout,'("  Maximum interaction distance: ",A)') string(e%dmax0,'f',8,4)
    end if
    write (uout,'("  Covering regions: ")')
    write (uout,'("    Total number of regions: ",A," (",2(A,X),A,")")') string(e%nregion), (string(e%nreg(j)),j=1,3)
    write (uout,'("    Minimum region ID: ",3(A,X))') (string(e%nmin(j)),j=1,3)
    write (uout,'("    Maximum region ID: ",3(A,X))') (string(e%nmax(j)),j=1,3)
    write (uout,'("    Region side (",A,"): ",A)') iunitname0(iunit), trim(string(e%boxsize * dunit0(iunit),'f',8,4))
    write (uout,'("    Transformation origin (",A,"): ",A,",",A,",",A)') iunitname0(iunit), &
       (trim(string(e%x0(j) * dunit0(iunit),'f',8,4)),j=1,3)
    write (uout,'("    Search offsets: ",A)') string(e%rs_nreg)
    write (uout,'("    Maximum search offset: ",A)') string(e%rs_imax)
    
    idavg = 0
    idmax = 0
    do j = 1, e%nregion
       idavg = idavg + (e%nrhi(j) - e%nrlo(j) + 1)
       idmax = max(idmax,e%nrhi(j) - e%nrlo(j) + 1,idmax)
    end do
    write (uout,'("    Average number of atoms per region: ",A)') trim(string(real(idavg,8)/real(e%nregion,8),'f',8,4))
    write (uout,'("    Maximum number of atoms in a region: ",A)') string(idmax)

    write (uout,*)

  end subroutine environ_report

  !xx! private procedures

  !> Calculate regions associated with the current environment and
  !> assign atoms to each region.
  subroutine calculate_regions(e)
    use grid1mod, only: agrid
    use global, only: cutrad
    use tools, only: qcksort
    use types, only: realloc
    use param, only: ctsq3, ctsq32, maxzat
    type(environ), intent(inout) :: e
    
    integer :: i, m
    integer, allocatable :: iord(:)
    integer :: i1, i2, i3, imax, nreg, iz, nregc(3)
    real*8 :: x0(3), x1(3), dist, rcut0
    real*8 :: xmin(3), xmax(3), xminc(3), xmaxc(3)

    ! find the encompassing boxes, for the main cell
    xminc = 1d40
    xmaxc = -1d40
    do i1 = 0, 1
       do i2 = 0, 1
          do i3 = 0, 1
             x0 = (/i1,i2,i3/)
             x0 = e%xr2c(x0)
             xminc = min(xminc,x0)
             xmaxc = max(xmaxc,x0)
          end do
       end do
    end do

    ! for the whole environment
    xmin = xminc
    xmax = xmaxc
    do i = e%ncell+1,e%n
       xmin = min(xmin,e%at(i)%r)
       xmax = max(xmax,e%at(i)%r)
    end do

    ! cap the box size at around 100^3 boxes
    e%boxsize = max(boxsize_default,2d0*e%rsph_env/100d0)

    ! calculate the position of the origin and the region partition
    nregc = ceiling(max((xmaxc - xminc) / e%boxsize,1d-14))
    e%x0 = xminc - 0.5d0 * (nregc * e%boxsize - (xmaxc - xminc))
    e%nmin = floor((xmin - e%x0) / e%boxsize)
    e%nmax = floor((xmax - e%x0) / e%boxsize)
    e%nreg = e%nmax - e%nmin + 1
    e%nregion = product(e%nreg)

    ! build the ordered list of atoms
    if (allocated(e%imap)) deallocate(e%imap)
    allocate(iord(e%n),e%imap(e%n))
    do i = 1, e%n
       iord(i) = e%c2i(e%at(i)%r)
       e%imap(i) = i
    end do
    call qcksort(iord,e%imap,1,e%n)

    ! limits for each region
    if (allocated(e%nrlo)) deallocate(e%nrlo)
    if (allocated(e%nrhi)) deallocate(e%nrhi)
    allocate(e%nrlo(e%nregion),e%nrhi(e%nregion))
    e%nrlo = 1
    e%nrhi = 0
    m = 0
    do i = 1, e%n
       if (iord(e%imap(i)) /= m) then
          if (i > 1) e%nrhi(m) = i-1
          e%nrlo(iord(e%imap(i))) = i
          e%nrhi(iord(e%imap(i))) = i
       end if
       m = iord(e%imap(i))
    end do
    e%nrhi(iord(e%imap(e%n))) = e%n
    deallocate(iord)

    ! Calculate the limits of the search region
    e%rs_imax = ceiling(e%rsph_env / e%boxsize - 0.5d0)
    e%rs_2imax1 = 2 * e%rs_imax + 1
    e%rs_nreg = e%rs_2imax1**3

    ! Take a reference point and calculate regions around it. The
    ! regions are sorted by distance to the reference point. This will
    ! become useful for nearest neighbor searches.
    if (allocated(e%rs_ioffset)) deallocate(e%rs_ioffset)
    if (allocated(e%rs_rcut)) deallocate(e%rs_rcut)
    allocate(e%rs_ioffset(e%rs_nreg),e%rs_rcut(e%rs_nreg))

    nreg = 0
    do i1 = -e%rs_imax, e%rs_imax
       do i2 = -e%rs_imax, e%rs_imax
          do i3 = -e%rs_imax, e%rs_imax
             nreg = nreg + 1
             x0 = max( ((/i1,i2,i3/) - 0.5d0),0d0)
             x1 = max(-((/i1,i2,i3/) + 0.5d0),0d0)

             ! dist = minimum distance from the origin to the (x0,x1) cube
             ! rcut = this cube has to be included in all searches where the maximum interaction
             !        distance is rcut or higher
             dist = sqrt(x0(1)*x0(1)+x0(2)*x0(2)+x0(3)*x0(3)+x1(1)*x1(1)+x1(2)*x1(2)+x1(3)*x1(3))
             rcut0 = max(dist - ctsq32,0d0) * e%boxsize

             e%rs_ioffset(nreg) = packoffset(i1,i2,i3,e%rs_imax,e%rs_2imax1)
             e%rs_rcut(nreg) = rcut0
          end do
       end do
    end do

    ! sort and populate the environ arrays
    allocate(iord(e%rs_nreg))
    do i = 1, e%rs_nreg
       iord(i) = i
    end do
    call qcksort(e%rs_rcut,iord,1,e%rs_nreg)
    e%rs_ioffset = e%rs_ioffset(iord)
    e%rs_rcut = e%rs_rcut(iord)
    deallocate(iord)

  end subroutine calculate_regions

  !> Pack an offset into a single integer index. nr are the offset
  !> limits (i1,i2,i3 = -(nr-1) .. nr-1). nr21 = 2 * nr + 1.
  function packoffset(i1,i2,i3,nr,nr21) result(ip)
    integer, intent(in) :: i1, i2, i3, nr, nr21
    integer :: ip

    ip = i1 + nr + nr21 * (i2 + nr + nr21 * (i3 + nr))

  end function packoffset

  !> Unpack an offset from a single integer index. nr are the offset
  !> limits (i1,i2,i3 = -(nr-1) .. nr-1). nr21 = 2 * nr + 1.
  function unpackoffset(ip,nr,nr21) result(ioff)
    integer, intent(in) :: ip, nr, nr21
    integer :: ioff(3)

    integer :: iaux, idd

    iaux = ip

    idd = modulo(iaux,nr21)
    ioff(1) = idd - nr
    iaux = (iaux - idd) / nr21

    idd = modulo(iaux,nr21)
    ioff(2) = idd - nr
    iaux = (iaux - idd) / nr21

    ioff(3) = iaux - nr

  end function unpackoffset

end submodule proc
