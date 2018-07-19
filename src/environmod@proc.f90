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

  real*8 :: boxsize_default = 4 ! length of the region side (bohr)

  !xx! private procedures
  ! subroutine calculate_regions(e)

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
    if (allocated(e%iaddregs)) deallocate(e%iaddregs)
    if (allocated(e%rcutregs)) deallocate(e%rcutregs)
    e%n = 0
    e%nspc = 0
    e%ncell = 0
    e%nregc = 0
    e%nreg = 0
    e%nmin = 0
    e%nmax = 0
    e%nregion = 0
    e%nregs = 0

  end subroutine environ_end
  
  !> Build an environment from molecule data. nspc = number of species.
  !> spc = species. n = number of atoms in the cell. at = atoms in the
  !> cell. m_xr2c = reduced crystallographic to Cartesian matrix. 
  !> m_x2xr = crystallographic to reduced crystallographic matrix.
  module subroutine environ_build_from_molecule(e,nspc,spc,n,at,m_xr2c,m_x2xr,m_x2c)
    use tools_math, only: matinv
    use types, only: realloc, celatom
    class(environ), intent(inout) :: e
    integer, intent(in) :: nspc
    type(species), intent(in) :: spc(nspc)
    integer, intent(in) :: n
    type(celatom), intent(in) :: at(n)
    real*8, intent(in) :: m_xr2c(3,3)
    real*8, intent(in) :: m_x2xr(3,3)
    real*8, intent(in) :: m_x2c(3,3)

    integer :: i 
    real*8 :: sphmax

    e%ismolecule = .true.
    e%nspc = nspc
    if (allocated(e%spc)) deallocate(e%spc)
    e%spc = spc

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
    e%sphmax = sphmax
    e%dmax0 = sphmax

    e%n = n
    if (allocated(e%at)) deallocate(e%at)
    allocate(e%at(e%n))
    do i = 1, n
       e%at(i)%x = at(i)%x
       e%at(i)%r = at(i)%r
       e%at(i)%idx = at(i)%idx
       e%at(i)%cidx = i
       e%at(i)%lvec = 0
       e%at(i)%is = at(i)%is
    end do
    e%ncell = n

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
    real*8 :: sphmax, dmax, x(3), xc(3), dd, xhalf(3)
    integer :: i1, i2, i3, i, imax, p(3), px(3)

    e%ismolecule = .false.
    e%nspc = nspc
    if (allocated(e%spc)) deallocate(e%spc)
    e%spc = spc

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
    e%sphmax = sphmax

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
             do i3 = -imax, imax
                if (abs(i1) /= imax .and. abs(i2) /= imax .and. abs(i3) /= imax) cycle

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
             end do
          end do
       end do
    end do
    call realloc(e%at,e%n)

    call calculate_regions(e)

  contains
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

    if (e%ismolecule) then
       xx = e%y2z(xx,icrd,ocrd)
       if (present(lvec)) lvec = 0
    else
       xx = e%y2z(xx,icrd,icrd_rcrys)
       if (present(lvec)) lvec = floor(xx)
       xx = xx - floor(xx)
       xx = e%y2z(xx,icrd_rcrys,ocrd)
    end if

  end subroutine y2z_center

  !> Cartesian to region transform
  pure module function c2p(e,xx) result(res)
    class(environ), intent(in) :: e
    real*8, intent(in)  :: xx(3)
    integer :: res(3)

    res = min(max(floor((xx - e%x0) / e%boxsize),e%nmin),e%nmax)

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

  !> Given the point xp (in icrd coordinates), calculates the nearest
  !> atom.  The nearest atom has ID nid from the complete list (atcel)
  !> and is at a distance dist. On output, the optional argument lvec
  !> contains the lattice vector to the nearest atom (i.e. its
  !> position is atcel(nid)%x + lvec). If nid0, consider only atoms
  !> with index nid0 from the non-equivalent list. If id0, consider
  !> only atoms with index id0 from the complete list.  If nozero,
  !> disregard zero-distance atoms. This routine is thread-safe.
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
    integer :: i, j, k, kmin, lvec0(3), nreg21(3)

    ! Find the integer region for the main cell copy of the input point
    x0 = xp
    call e%y2z_center(x0,icrd,icrd_cart,lvec0)
    ireg0 = e%c2p(x0)
    
    ! run over regions sorted by distance
    nreg21 = 2 * e%nreg + 1
    distmin = 1d40
    kmin = 0
    main: do i = 1, e%nregs
       ireg = ireg0 + unpackoffset(e%iaddregs(i),e%nreg,nreg21)
       if (any(ireg < e%nmin) .or. any(ireg > e%nmax)) cycle
       idxreg = e%p2i(ireg)
       if (e%nrhi(idxreg) == 0) cycle

       do j = e%nrlo(idxreg), e%nrhi(idxreg)
          k = e%imap(j)
          if (distmin < e%rcutregs(i)) exit main
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

    nid = e%at(kmin)%cidx
    dist = distmin
    if (present(lvec)) then
       lvec = e%at(kmin)%lvec + nint(e%xr2x(real(lvec0,8)))
    end if

  end subroutine nearest_atom

  !> Given the point xp (in icrd coordinates), calculates the list of
  !> nearest atoms. The list is sorted by distance and shell if sorted
  !> is true (using up2n or up2sh makes sorted = true). The output
  !> list contains nat atoms with IDs nid(1:nat) from the environment,
  !> distances to the input point equal to dist(1:nat) and lattice
  !> vector lvec in cryst. coords. The position of atom i in
  !> cryst. coords. is e%at(nid(i))%x + lvec. Optionally, ishell(i)
  !> contains the shell ID for atom i in the output list. One or more
  !> of four cutoff criteria must be chosen: list up to a distance
  !> up2d, list up to a species-dependent distance (up2dsp), up to a
  !> number of shells up2sh or up to a number of atoms up2n. If nid0,
  !> consider only atoms with index nid0 from the non-equivalent
  !> list. If id0, consider only atoms with index id0 from the
  !> complete list. If nozero, disregard zero-distance atoms.
  module subroutine list_near_atoms(e,xp,icrd,sorted,nat,eid,dist,lvec,ishell0,up2d,up2dsp,up2sh,up2n,nid0,id0,nozero)
    use global, only: atomeps
    use tools_io, only: ferror, faterr
    use tools, only: qcksort, iqcksort
    use types, only: realloc
    use param, only: icrd_cart
    class(environ), intent(in) :: e
    real*8, intent(in) :: xp(3)
    integer, intent(in) :: icrd
    logical, intent(in) :: sorted
    integer, intent(out) :: nat
    integer, allocatable, intent(inout) :: eid(:)
    real*8, allocatable, intent(inout) :: dist(:)
    integer, intent(out) :: lvec(3)
    integer, allocatable, intent(inout), optional :: ishell0(:)
    real*8, intent(in), optional :: up2d
    real*8, intent(in), optional :: up2dsp(:)
    integer, intent(in), optional :: up2sh
    integer, intent(in), optional :: up2n
    integer, intent(in), optional :: nid0
    integer, intent(in), optional :: id0
    logical, intent(in), optional :: nozero

    real*8, parameter :: eps = 1d-10

    real*8 :: x0(3), dist0, rcutshel, rcutn
    integer :: ireg0(3), ireg(3), idxreg
    integer :: i, j, k, l, lthis
    integer, allocatable :: iord(:), iiord(:), ishell(:)
    integer :: nshel, nreg21(3)
    real*8, allocatable :: rshel(:)
    integer, allocatable :: idxshel(:)
    logical :: doshell

    if (.not.present(up2d).and..not.present(up2dsp).and..not.present(up2sh).and..not.present(up2n)) &
       call ferror("list_near_atoms","must give one of up2d, up2dsp, up2sh, or up2n",faterr)
    doshell = present(ishell0) .or. present(up2sh)

    ! Find the integer region for the main cell copy of the input point
    x0 = xp
    call e%y2z_center(x0,icrd,icrd_cart,lvec)
    ireg0 = e%c2p(x0)
    lvec = nint(e%xr2x(real(lvec,8)))

    ! initialize the output arrays
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

    nreg21 = 2 * e%nreg + 1
    ! run over regions around ireg0 and over atoms belonging to those regions
    main: do i = 1, e%nregs
       ireg = ireg0 + unpackoffset(e%iaddregs(i),e%nreg,nreg21)
       if (any(ireg < e%nmin) .or. any(ireg > e%nmax)) cycle
       idxreg = e%p2i(ireg)
       if (e%nrhi(idxreg) == 0) cycle

       do j = e%nrlo(idxreg), e%nrhi(idxreg)
          k = e%imap(j)

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
             if (e%rcutregs(i) > up2d) exit main
             if (dist0 > up2d) cycle
          end if
          if (present(up2dsp)) then
             if (e%rcutregs(i) > maxval(up2dsp)) exit main
             if (dist0 > up2dsp(e%at(k)%is)) cycle
          end if
          if (present(up2sh)) then
             if (nshel >= up2sh .and. e%rcutregs(i) > rcutshel) &
                exit main
          end if
          if (present(up2n)) then
             if (nat >= up2n .and. e%rcutregs(i) > rcutn) &
                exit main
          end if

          ! ! add this atom
          nat = nat + 1
          if (nat > size(eid,1)) then
             call realloc(eid,2*nat)
             call realloc(dist,2*nat)
             if (doshell) call realloc(ishell,2*nat)
          end if
          eid(nat) = k
          dist(nat) = dist0

          ! update the up2n rcut. All remaining atoms will have
          ! dist > than the current region's rcut. rcutn will be >
          ! than all atom distances from 1 to up2n.
          if (present(up2n)) then
             if (nat <= up2n) rcutn = max(dist0,rcutn)
          end if

          ! process the shells
          if (doshell) then
             ! see if this shell is already known
             lthis = 0
             do l = 1, nshel
                if ((abs(dist0 - rshel(l)) < atomeps) .and. (idxshel(l) == e%at(k)%idx)) then
                   lthis = l
                   exit
                end if
             end do

             ! if not known, create the new shell
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

             ! update the up2sh rcut. All remaining atoms will have
             ! dist > than the current region's rcut. rcutsh will be >
             ! than all atom distances from 1 to shell up2sh.
             if (present(up2sh)) then
                if (up2sh <= nshel) rcutshel = max(dist0,rcutshel)
             end if

             ! write down the shell for this atom
             ishell(nat) = lthis
          end if
       end do
    end do main

    ! rearrange the arrays 
    if (nat > 0) then
       ! first reallocation
       call realloc(eid,nat)
       call realloc(dist,nat)
       if (doshell) call realloc(ishell,nat)

       ! re-order shell labels
       if (doshell) then
          allocate(iord(nshel),iiord(nshel))
          do i = 1, nshel
             iord(i) = i
          end do
          call iqcksort(idxshel,iord,1,nshel)
          call qcksort(rshel,iord,1,nshel)
          do i = 1, nshel
             iiord(iord(i)) = i
          end do
          deallocate(iord)
          do i = 1, nat
             ishell(i) = iiord(ishell(i))
          end do
          deallocate(iiord)
       end if

       if (sorted .or. present(up2sh) .or. present(up2n)) then
          ! sort by distance and by shell
          allocate(iord(nat))
          do i = 1, nat
             iord(i) = i
          end do
          if (doshell) call iqcksort(ishell,iord,1,nat)
          call qcksort(dist,iord,1,nat)
          eid = eid(iord)
          dist = dist(iord)
          if (doshell) ishell = ishell(iord)
          deallocate(iord)

          ! prune the extra atoms
          if (present(up2sh) .or. present(up2n)) then
             if (present(up2sh)) then
                do i = 1, nat
                   if (ishell(i) > up2sh) then
                      nat = i-1
                      exit
                   end if
                end do
             end if
             if (present(up2n)) nat = up2n
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

  end subroutine list_near_atoms
    
  !> Calculate the core (if zpsp is present) or promolecular densities
  !> at a point x0 (coord format given by icrd) using atomic radial
  !> grids up to a number of derivatives nder (max: 2). Returns the
  !> density (f), gradient (fp, nder >= 1), and Hessian (fpp, nder >=
  !> 2). If a fragment (fr) is given, then only the atoms in it
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

    integer :: i, j, k, ii, iz
    real*8 :: xc(3), xx(3), rlvec(3), r, rinv1, rinv1rp
    real*8 :: rho, rhop, rhopp, rfac, radd
    logical :: iscore
    type(grid1), pointer :: g
    integer :: nat, lvec(3)
    integer, allocatable :: nid(:)
    real*8, allocatable :: dist(:), rcutmax(:)

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
    allocate(rcutmax(e%nspc))
    rcutmax = 0d0
    do i = 1, e%nspc
       iz = e%spc(i)%z
       if (iz == 0 .or. iz > maxzat) cycle
       if (iscore) then
          g => cgrid(iz,zpsp(iz))
       else
          g => agrid(iz)
       end if
       rcutmax(i) = min(cutrad(iz),g%rmax)
    end do
    call e%list_near_atoms(xc,icrd_cart,.false.,nat,nid,dist,lvec,up2dsp=rcutmax)
    rlvec = lvec
    rlvec = e%x2c(rlvec)
    deallocate(rcutmax)

    ! Do the density and derivatives sum
    do ii = 1, nat
       i = nid(ii)
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

  end subroutine promolecular

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
       iunitname0(iunit), trim(string(e%sphmax*dunit0(iunit),'f',8,4)), trim(string((e%sphmax+e%dmax0)*dunit0(iunit),'f',8,4))
    write (uout,'("  Maximum interaction distance: ",A)') string(e%dmax0,'f',8,4)
    write (uout,'("  Covering regions: ")')
    write (uout,'("    Unit cell: ",3(A,X))') (string(e%nregc(j)),j=1,3)
    write (uout,'("    Environment: ",3(A,X))') (string(e%nreg(j)),j=1,3)
    write (uout,'("    Search offsets: ",A)') string(e%nregs)
    write (uout,'("    Total: ",A)') string(e%nregion)
    write (uout,'("    Region side (",A,"): ",A)') iunitname0(iunit), trim(string(e%boxsize * dunit0(iunit),'f',8,4))
    write (uout,'("    Transformation origin (",A,"): ",A,",",A,",",A)') iunitname0(iunit), &
       (trim(string(e%x0(j) * dunit0(iunit),'f',8,4)),j=1,3)
    write (uout,'("    Minimum region ID: ",3(A,X))') (string(e%nmin(j)),j=1,3)
    write (uout,'("    Maximum region ID: ",3(A,X))') (string(e%nmax(j)),j=1,3)
    
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
    use tools, only: iqcksort, qcksort
    use types, only: realloc
    type(environ), intent(inout) :: e
    
    integer :: i, m
    integer, allocatable :: iord(:)
    logical :: dorepeat
    integer :: i1, i2, i3, imax, nreg, nreg21(3)
    real*8 :: x0(3), x1(3), dist, rcut0

    ! find the encompassing boxes, for the main cell
    e%xminc = 1d40
    e%xmaxc = -1d40
    do i = 1, e%ncell
       e%xminc = min(e%xminc,e%at(i)%r)
       e%xmaxc = max(e%xmaxc,e%at(i)%r)
    end do

    ! for the whole environment
    e%xmin = e%xminc
    e%xmax = e%xmaxc
    do i = e%ncell+1,e%n
       e%xmin = min(e%xmin,e%at(i)%r)
       e%xmax = max(e%xmax,e%at(i)%r)
    end do

    ! calculate the position of the origin and the region partition
    e%boxsize = boxsize_default
    e%nregc = ceiling((e%xmaxc - e%xminc) / e%boxsize)
    e%x0 = e%xminc - 0.5d0 * (e%nregc * e%boxsize - (e%xmaxc - e%xminc))
    e%nmin = floor((e%xmin - e%x0) / e%boxsize)
    e%nmax = floor((e%xmax - e%x0) / e%boxsize)
    e%nreg = e%nmax - e%nmin + 1
    e%nregion = product(e%nreg)

    ! build the ordered list of atoms
    if (allocated(e%imap)) deallocate(e%imap)
    allocate(iord(e%n),e%imap(e%n))
    do i = 1, e%n
       iord(i) = e%c2i(e%at(i)%r)
       e%imap(i) = i
    end do
    call iqcksort(iord,e%imap,1,e%n)

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

    ! Take a reference point and calculate regions around it. The
    ! regions are sorted by distance to the reference point. This will
    ! become useful for nearest neighbor searches.
    if (allocated(e%iaddregs)) deallocate(e%iaddregs)
    if (allocated(e%rcutregs)) deallocate(e%rcutregs)
    allocate(e%iaddregs(e%nregion),e%rcutregs(e%nregion))
    nreg21 = 2 * e%nreg + 1
    dorepeat = .true.
    imax = 0
    nreg = 1
    e%iaddregs(1) = packoffset(0,0,0,e%nreg,nreg21)
    e%rcutregs(1) = 0d0
    do while (dorepeat)
       dorepeat = .false.
       imax = imax + 1
       if (all(imax >= e%nreg)) exit

       do i1 = -imax, imax
          do i2 = -imax, imax
             do i3 = -imax, imax
                if (abs(i1) /= imax .and. abs(i2) /= imax .and. abs(i3) /= imax) cycle
                if (abs(i1) >= e%nreg(1) .or. abs(i2) >= e%nreg(2) .or. abs(i3) >= e%nreg(3)) cycle
                x0 = e%boxsize * ((/i1,i2,i3/) - 0.5d0)
                x1 = e%boxsize * ((/i1,i2,i3/) + 0.5d0)

                ! dist = minimum distance from the origin to the (x0,x1) cube
                ! rcut = this cube has to be included in all searches where the maximum interaction
                !        distance is rcut or higher
                dist = sqrt(max(x0(1),0d0)**2 + max(-x1(1),0d0)**2 + max(x0(2),0d0)**2 + max(-x1(2),0d0)**2 + &
                   max(x0(3),0d0)**2 + max(-x1(3),0d0)**2)
                rcut0 = max(dist - e%boxsize * sqrt(3d0) / 2d0,0d0)

                if (rcut0 < 1.5d0 * e%dmax0) then
                   dorepeat = .true. 

                   nreg = nreg + 1
                   if (nreg > size(e%iaddregs,1)) then
                      call realloc(e%iaddregs,2*nreg)
                      call realloc(e%rcutregs,2*nreg)
                   end if
                   ! because the p -> i transformation is linear, we can write the offset as a single
                   ! integer
                   e%iaddregs(nreg) = packoffset(i1,i2,i3,e%nreg,nreg21)
                   e%rcutregs(nreg) = rcut0
                end if
             end do
          end do
       end do
    end do
    call realloc(e%iaddregs,nreg)
    call realloc(e%rcutregs,nreg)

    ! sort and populate the environ arrays
    allocate(iord(nreg))
    do i = 1, nreg
       iord(i) = i
    end do
    call qcksort(e%rcutregs,iord,1,nreg)
    e%iaddregs = e%iaddregs(iord)
    e%rcutregs = e%rcutregs(iord)
    e%nregs = nreg
    deallocate(iord)

  end subroutine calculate_regions

  !> Pack an offset into a single integer index. nr are the offset
  !> limits (i1,i2,i3 = -(nr-1) .. nr-1). nr21 = 2 * nr + 1.
  function packoffset(i1,i2,i3,nr,nr21) result(ip)
    integer, intent(in) :: i1, i2, i3, nr(3), nr21(3)
    integer :: ip

    ip = i1 + nr(1) + nr21(1) * (i2 + nr(2) + nr21(2) * (i3 + nr(3)))

  end function packoffset

  !> Unpack an offset from a single integer index. nr are the offset
  !> limits (i1,i2,i3 = -(nr-1) .. nr-1). nr21 = 2 * nr + 1.
  function unpackoffset(ip,nr,nr21) result(ioff)
    integer, intent(in) :: ip, nr(3), nr21(3)
    integer :: ioff(3)

    integer :: iaux, idd

    iaux = ip

    idd = modulo(iaux,nr21(1))
    ioff(1) = idd - nr(1)
    iaux = (iaux - idd) / nr21(1)

    idd = modulo(iaux,nr21(2))
    ioff(2) = idd - nr(2)
    iaux = (iaux - idd) / nr21(2)

    ioff(3) = iaux - nr(3)

  end function unpackoffset

end submodule proc
