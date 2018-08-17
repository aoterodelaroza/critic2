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
  
  !> Build an environment from molecular or crystal data. ismol = true
  !> if the environment is for a molecule (false = crystal). nspc =
  !> number of species.  spc = species. n = number of atoms in the
  !> cell. at = atoms in the cell. m_xr2c = reduced crystallographic
  !> to Cartesian matrix.  m_x2xr = crystallographic to reduced
  !> crystallographic matrix. dmax0 = the environment will contain all
  !> atoms within a distance dmax0 of any atom in the unit
  !> cell/molecule. If at_in_xr is present, it contains the unit cell
  !> atoms in reduced crystallographic coords with anyatom type, and
  !> at(:) is ignored (used to load from another environment).
  module subroutine environ_build(e,ismol,nspc,spc,n,at,m_xr2c,m_x2xr,m_x2c,dmax0,at_in_xr)
    use global, only: cutrad
    use tools, only: qcksort
    use tools_math, only: matinv
    use types, only: realloc, species, celatom
    use param, only: atmcov, ctsq32
    class(environ), intent(inout) :: e
    logical, intent(in) :: ismol
    integer, intent(in) :: nspc
    type(species), intent(in) :: spc(nspc)
    integer, intent(in) :: n
    type(celatom), intent(in) :: at(n)
    real*8, intent(in) :: m_xr2c(3,3)
    real*8, intent(in) :: m_x2xr(3,3)
    real*8, intent(in) :: m_x2c(3,3)
    real*8, intent(in), optional :: dmax0
    type(anyatom), intent(in), optional :: at_in_xr(:)

    logical :: dorepeat
    real*8 :: sphmax, dmax, x(3), rmin(3), rmax(3)
    real*8 :: x0(3), x1(3), dist, rcut0
    integer :: i1, i2, i3, i, m, imax, i3min, nreg
    integer, allocatable :: iord(:)

    ! fill the matrices
    e%m_xr2c = m_xr2c
    e%m_c2xr = matinv(m_xr2c)
    e%m_x2xr = m_x2xr
    e%m_xr2x = matinv(m_x2xr)
    e%m_x2c = m_x2c
    e%m_c2x = matinv(m_x2c)

    ! do nothing else if there are no atoms
    e%n = n
    e%nspc = nspc
    if (n == 0 .or. nspc == 0) return

    ! species and ismolecule
    e%ismolecule = ismol
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

    ! calculate the maximum diagonal half-length (sphmax)
    sphmax = norm2(e%xr2c((/0d0,0d0,0d0/) - (/0.5d0,0.5d0,0.5d0/)))
    sphmax = max(sphmax,norm2(e%xr2c((/1d0,0d0,0d0/) - (/0.5d0,0.5d0,0.5d0/))))
    sphmax = max(sphmax,norm2(e%xr2c((/0d0,1d0,0d0/) - (/0.5d0,0.5d0,0.5d0/))))
    sphmax = max(sphmax,norm2(e%xr2c((/0d0,0d0,1d0/) - (/0.5d0,0.5d0,0.5d0/))))
    e%rsph_uc = sphmax

    ! calculate the e%dmax0 if not given
    if (present(dmax0)) then
       e%dmax0 = dmax0 + 1d-2
    else
       e%dmax0 = 0d0
       do i = 1, nspc
          if (spc(i)%z > 0) e%dmax0 = max(e%dmax0,cutrad(spc(i)%z))
       end do
       e%dmax0 = e%dmax0 + 1d-2
    end if

    ! parepare to write down the first n atoms
    if (allocated(e%at)) deallocate(e%at)
    allocate(e%at(e%n))

    ! origin is the center of the unit cell
    e%x0 = 0.5d0
    e%x0 = e%xr2c(e%x0)

    if (ismol) then
       ! A molecule: add all atoms in the molecule to the environment
       rmin = 1d40
       rmax = -1d40
       do i = 1, n
          if (present(at_in_xr)) then
             e%at(i)%x = at_in_xr(i)%x
             e%at(i)%r = at_in_xr(i)%r
             e%at(i)%idx = at_in_xr(i)%idx
             e%at(i)%is = at_in_xr(i)%is
          else
             e%at(i)%x = at(i)%x
             e%at(i)%r = at(i)%r
             e%at(i)%idx = at(i)%idx
             e%at(i)%is = at(i)%is
          end if
          rmin = min(e%at(i)%r,rmin)
          rmax = max(e%at(i)%r,rmax)
          e%at(i)%cidx = i
          e%at(i)%lvec = 0
       end do
       e%ncell = n
    else
       ! A crystal
       ! add all atoms in the main reduced cell (coords: 0 -> 1)
       rmin = 1d40
       rmax = -1d40
       do i = 1, n
          if (present(at_in_xr)) then
             e%at(i)%x = at_in_xr(i)%x
             e%at(i)%x = e%at(i)%x - floor(e%at(i)%x)
             e%at(i)%lvec = nint(e%xr2x(e%at(i)%x - at_in_xr(i)%x))
             e%at(i)%idx = at_in_xr(i)%idx
             e%at(i)%is = at_in_xr(i)%is
          else
             e%at(i)%x = e%x2xr(at(i)%x)
             e%at(i)%x = e%at(i)%x - floor(e%at(i)%x)
             e%at(i)%lvec = nint(e%xr2x(e%at(i)%x) - at(i)%x)
             e%at(i)%idx = at(i)%idx
             e%at(i)%is = at(i)%is
          end if
          e%at(i)%r = e%xr2c(e%at(i)%x)
          rmin = min(e%at(i)%r,rmin)
          rmax = max(e%at(i)%r,rmax)
          e%at(i)%cidx = i
       end do
       e%ncell = n

       ! include all atoms at a distance sphmax+dmax from the center of the cell
       dorepeat = .true.
       imax = 0
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
    end if

    ! sphere that circumscribes the environment
    e%rsph_env = 0.5d0 * norm2(rmax-rmin)

    ! cap the box size at around 100^3 boxes
    e%boxsize = max(e%boxsize,2d0*e%rsph_env/100d0)

    ! minimum and maximum region
    e%nmin = floor((rmin - e%x0) / e%boxsize)
    e%nmax = floor((rmax - e%x0) / e%boxsize)

    ! number of regions
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

    ! Calculate the limits of the search region. This choice of
    ! rs_imax ensures that rs_rcut(i) > e%rsph_env and e%dmax0 for all
    ! the regions no covered by the offsets.
    e%rs_imax = ceiling(max(e%dmax0,e%rsph_env) / e%boxsize)
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
         dd = norm2(xc - e%x0)
         if (dd < sphmax+e%dmax0) then
            dorepeat = .true.
            e%n = e%n + 1
            if (e%n > size(e%at,1)) call realloc(e%at,2*e%n)

            e%at(e%n)%x = x
            e%at(e%n)%r = xc
            rmin = min(e%at(e%n)%r,rmin)
            rmax = max(e%at(e%n)%r,rmax)
            e%at(e%n)%idx = e%at(i)%idx
            e%at(e%n)%cidx = e%at(i)%cidx
            e%at(e%n)%lvec = e%at(i)%lvec + px
            e%at(e%n)%is = e%at(i)%is
         end if
      end do

    end subroutine addcell
  endsubroutine environ_build

  !> Extend an environment e0 to a new interaction distance equal to dmax0.
  module subroutine environ_extend(e,e0,dmax0)
    use types, only: celatom
    class(environ), intent(inout) :: e
    type(environ), intent(in) :: e0
    real*8, intent(in) :: dmax0

    type(celatom) :: at(1)

    call e%end()
    call e%build(e0%ismolecule,e0%nspc,e0%spc(1:e0%nspc),e0%ncell,at,e0%m_xr2c,e0%m_x2xr,e0%m_x2c,dmax0,e0%at(1:e0%ncell))

  end subroutine environ_extend

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

  !> Given point x0 (with icrd input coordinates), translate to the
  !> main cell if the environment is from a crystal. Then if x0
  !> corresponds to an atomic position (to within distmax or atomeps
  !> if distmax is not given), return the complete-list ID of the
  !> atom. Otherwise, return 0. Optionally, return the lattice vector
  !> translation (lvec) and the distance (dist) to the closest atom.
  !> This routine is thread-safe.
  module function identify_atom(e,x0,icrd,lvec,dist,distmax)
    use global, only: atomeps
    use tools_io, only: ferror, faterr
    class(environ), intent(in) :: e
    real*8, intent(in) :: x0(3)
    integer, intent(in) :: icrd
    integer, intent(out), optional :: lvec(3)
    real*8, intent(out), optional :: dist
    real*8, intent(in), optional :: distmax
    integer :: identify_atom
    
    real*8 :: eps, dist0
    integer :: ierr, lvec0(3)

    eps = atomeps
    if (present(distmax)) eps = distmax

    call e%nearest_atom_short(x0,icrd,eps,identify_atom,lvec0,dist0,ierr)
    if (ierr == 0 .or. ierr == 1) then
       if (present(lvec)) lvec = lvec0
       if (present(dist)) dist = dist0
       return
    end if

    call e%nearest_atom_long(x0,icrd,eps,identify_atom,lvec0,dist0,ierr)
    if (ierr == 2) then
       call ferror("identify_atom","distmax too large for environment boxsize",faterr)
    end if
    if (present(lvec)) lvec = lvec0
    if (present(dist)) dist = dist0

  end function identify_atom

  !> Given the point xp (in icrd coordinates), translate to the main
  !> cell if the environment is from a crystal. Then, calculate the
  !> nearest atom up to a distance distmax. The nearest atom has ID
  !> nid from the complete list (atcel) and is at a distance dist, or
  !> nid=0 and dist=distmax if the search did not produce any atoms.
  !> The default distmax is the environment's dmax0.  On output, the
  !> optional argument lvec contains the lattice vector to the nearest
  !> atom (i.e. its position is atcel(nid)%x + lvec). If cidx,
  !> consider only atoms with index cidx0 from the complete list. If
  !> idx0, consider only atoms with index id0 from the non-equivalent
  !> list. If nozero, disregard zero-distance atoms. This routine is
  !> thread-safe.
  module subroutine nearest_atom(e,xp,icrd,nid,dist,distmax,lvec,cidx0,idx0,nozero)
    use param, only: icrd_cart
    class(environ), intent(in) :: e
    real*8, intent(in) :: xp(3)
    integer, intent(in) :: icrd
    integer, intent(out) :: nid
    real*8, intent(out) :: dist
    real*8, intent(in), optional :: distmax
    integer, intent(out), optional :: lvec(3)
    integer, intent(in), optional :: cidx0
    integer, intent(in), optional :: idx0
    logical, intent(in), optional :: nozero

    real*8 :: eps
    integer :: ierr, lvec0(3)

    eps = e%dmax0
    if (present(distmax)) eps = distmax

    call e%nearest_atom_short(xp,icrd,eps,nid,lvec0,dist,ierr,cidx0,idx0,nozero)
    if (ierr == 0 .or. ierr == 1) then
       if (present(lvec)) lvec = lvec0
       return
    end if

    call e%nearest_atom_long(xp,icrd,eps,nid,lvec0,dist,ierr,cidx0,idx0,nozero)
    if (present(lvec)) lvec = lvec0
    if (ierr == 2) then
       nid = 0
       dist = eps
    end if

  end subroutine nearest_atom
  
  !> Given point x0 (with icrd input coordinates), translate to the
  !> main cell if the environment is from a crystal. Then, calculate
  !> the environment atom closest to x0 up to a distance of
  !> distmax. If an atom was found within distmax, return the
  !> complete list ID (cidx), the lattice vector translation to the
  !> closest point (lvec), the distance (dist), and ierr = 0. If no
  !> atom was found, return eid = 0 and: a) ierr = 1 if it is
  !> guaranteed that there are no atoms up to distmax or b) ierr = 2
  !> if distmax exceeds half the environment's boxsize and therefore
  !> we can be sure there are no atoms only up to distmax=boxisze/2.
  !> If cidx0, consider only atoms with ID cidx0 from the complete
  !> list. If idx0, consider only atoms with ID idx0 from the
  !> non-equivalent atom list. If nozero, discard atoms very close to
  !> x0 (1d-10 default). This routine is limited to a distmax equal to
  !> half the environment's boxsize, but it is significantly faster
  !> than nearest_atom_long. Thread-safe.
  module subroutine nearest_atom_short(e,x0,icrd,distmax,cidx,lvec,dist,ierr,cidx0,idx0,nozero)
    use param, only: icrd_cart
    class(environ), intent(in) :: e
    real*8, intent(in) :: x0(3)
    integer, intent(in) :: icrd
    real*8, intent(in) :: distmax
    integer, intent(out) :: cidx
    integer, intent(out) :: lvec(3)
    real*8, intent(out) :: dist
    integer, intent(out) :: ierr
    integer, intent(in), optional :: cidx0
    integer, intent(in), optional :: idx0
    logical, intent(in), optional :: nozero

    real*8, parameter :: epsmall2 = 1d-20
    real*8 :: xp(3), x0r(3), x1r(3), dist, d2, d2min, eps2, x(3), d0
    integer :: ireg0(3), imin(3), imax(3), i, j, k, i1, i2, i3
    integer :: ireg(3), idx, kmin
    
    ! initialize
    d0 = min(distmax,0.5d0 * e%boxsize - 1d-10)
    dist = 0d0
    lvec = 0
    ierr = 0
    eps2 = d0 * d0

    ! bring the atom to the main cell, convert, identify the region
    xp = x0
    call e%y2z_center(xp,icrd,icrd_cart,lvec)
    ireg0 = e%c2p(xp)

    ! calculate the regions to explore
    imin = 0
    imax = 0
    x0r = e%x0 + e%boxsize * ireg0
    x1r = e%x0 + e%boxsize * (ireg0+1)
    do i = 1, 3
       if (xp(i)-x0r(i) <= d0) imin(i) = -1
       if (x1r(i)-xp(i) <= d0) imax(i) = 1
    end do

    ! Identify the atom
    kmin = 0
    d2min = 1d40
    do i1 = imin(1), imax(1)
       do i2 = imin(2), imax(2)
          do i3 = imin(3), imax(3)
             ireg = ireg0 + (/i1,i2,i3/)
             if (any(ireg < e%nmin) .or. any(ireg > e%nmax)) cycle
             idx = e%p2i(ireg)
             do j = e%nrlo(idx), e%nrhi(idx)
                k = e%imap(j)
                if (present(cidx0)) then
                   if (e%at(k)%cidx /= cidx0) cycle
                end if
                if (present(idx0)) then
                   if (e%at(k)%idx /= idx0) cycle
                end if

                x = e%at(k)%r - xp
                d2 = x(1)*x(1)+x(2)*x(2)+x(3)*x(3)
                if (present(nozero)) then
                   if (nozero .and. d2 < epsmall2) cycle
                end if

                if (d2 < d2min .and. d2 < eps2) then
                   kmin = k
                   d2min = d2
                end if
             end do
          end do
       end do
    end do

    ! return the atom index and other info
    if (kmin > 0 .and. d2min < eps2) then
       cidx = e%at(kmin)%cidx
       dist = sqrt(d2min)
       lvec = e%at(kmin)%lvec + nint(e%xr2x(real(lvec,8)))
       ierr = 0
    else
       cidx = 0
       dist = 0d0
       lvec = 0
       if (d0 == distmax) then
          ierr = 1
       else
          ierr = 2
       end if
    end if

  end subroutine nearest_atom_short

  !> Given point x0 (with icrd input coordinates), translate to the
  !> main cell if the environment is from a crystal. Then, calculate
  !> the environment atom closest to x0 up to a distance of
  !> distmax. If an atom was found within distmax, return the
  !> complete list ID (cidx), the lattice vector translation to the
  !> closest point (lvec), the distance (dist), and ierr = 0. If no
  !> atom was found, return eid = 0 and: a) ierr = 1 if it is
  !> guaranteed that there are no atoms up to distmax or b) ierr = 2
  !> if distmax exceeds half the environment's boxsize and therefore
  !> we can be sure there are no atoms only up to distmax=boxisze/2.
  !> If cidx0, consider only atoms with ID cidx0 from the complete
  !> list. If idx0, consider only atoms with ID idx0 from the
  !> non-equivalent atom list. If nozero, discard atoms very close to
  !> x0 (1d-10 default). This routine is limited to a distmax equal to
  !> the dmax0 of the environment. It is slower than the short-range
  !> version, nearest_atom_short. Thread-safe.
  module subroutine nearest_atom_long(e,x0,icrd,distmax,cidx,lvec,dist,ierr,cidx0,idx0,nozero)
    use param, only: icrd_cart
    class(environ), intent(in) :: e
    real*8, intent(in) :: x0(3)
    integer, intent(in) :: icrd
    real*8, intent(in) :: distmax
    integer, intent(out) :: cidx
    integer, intent(out) :: lvec(3)
    real*8, intent(out) :: dist
    integer, intent(out) :: ierr
    integer, intent(in), optional :: cidx0
    integer, intent(in), optional :: idx0
    logical, intent(in), optional :: nozero

    real*8, parameter :: epsmall2 = 1d-20
    real*8 :: x(3), xp(3)
    integer :: ireg0(3), ireg(3), kmin, idxreg
    integer :: i, j, k
    real*8 :: d0, eps2, d2min, d2, rcut2

    d0 = min(distmax,e%dmax0)
    dist = 0d0
    lvec = 0
    ierr = 0
    eps2 = d0 * d0

    ! Find the integer region for the main cell copy of the input point
    xp = x0
    call e%y2z_center(xp,icrd,icrd_cart,lvec)
    ireg0 = e%c2p(xp)
    
    ! run over regions sorted by distance
    kmin = 0
    d2min = 1d40
    do i = 1, e%rs_nreg
       ireg = ireg0 + unpackoffset(e%rs_ioffset(i),e%rs_imax,e%rs_2imax1)
       if (any(ireg < e%nmin) .or. any(ireg > e%nmax)) cycle

       rcut2 = e%rs_rcut(i) * e%rs_rcut(i)
       if (d2min < rcut2 .or. eps2 < rcut2) exit
       idxreg = e%p2i(ireg)
       
       do j = e%nrlo(idxreg), e%nrhi(idxreg)
          k = e%imap(j)
          if (present(cidx0)) then
             if (e%at(k)%cidx /= cidx0) cycle
          end if
          if (present(idx0)) then
             if (e%at(k)%idx /= idx0) cycle
          end if

          x = e%at(k)%r - xp
          d2 = x(1)*x(1)+x(2)*x(2)+x(3)*x(3)
          if (present(nozero)) then
             if (nozero .and. d2 < epsmall2) cycle
          end if

          if (d2 < d2min .and. d2 < eps2) then
             kmin = k
             d2min = d2
          end if
       end do
    end do

    ! return the atom index and other info
    if (kmin > 0 .and. d2min < eps2) then
       cidx = e%at(kmin)%cidx
       dist = sqrt(d2min)
       lvec = e%at(kmin)%lvec + nint(e%xr2x(real(lvec,8)))
       ierr = 0
    else
       cidx = 0
       dist = 0d0
       lvec = 0
       if (d0 == distmax) then
          ierr = 1
       else
          ierr = 2
       end if
    end if

  end subroutine nearest_atom_long

  !> Given the point xp (in icrd coordinates), center the point in the
  !> main cell if the environment is from a crystal. Then, calculate
  !> the list of nearest atoms. If sorted is true, the list is sorted
  !> by distance and shell (if up2n or up2sh is used, the list is
  !> always sorted). The output list eid contains nat atoms with IDs
  !> nid(1:nat) from the environment, distances to the input point
  !> equal to dist(1:nat) and lattice vector lvec in
  !> cryst. coords. The position of atom i in cryst. coords. is
  !> e%at(nid(i))%x + lvec. Optionally, ishell0(i) contains the shell
  !> ID for atom i in the output list. One or more of four cutoff
  !> criteria must be chosen: list up to a distance up2d, between a
  !> minimum and a maximum species-dependent distance (up2dsp), up
  !> to a distance to each atom in the complete list (up2dcidx), up to
  !> a number of shells up2sh or up to a number of atoms up2n. If
  !> nid0, consider only atoms with index nid0 from the non-equivalent
  !> list. If id0, consider only atoms with index id0 from the
  !> complete list. If iz0, consider only atoms with atomic number
  !> iz0.  If nozero, disregard zero-distance atoms. The output error
  !> condition ierr is 0 if the search was successful or non-zero if
  !> the input search conditions could not be met by this environment.
  module subroutine list_near_atoms(e,xp,icrd,sorted,nat,eid,dist,lvec,ierr,ishell0,up2d,up2dsp,up2dcidx,up2sh,up2n,nid0,id0,iz0,nozero)
    use global, only: atomeps
    use tools_io, only: ferror, faterr
    use tools, only: mergesort
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
    integer, intent(out) :: ierr
    integer, allocatable, intent(inout), optional :: ishell0(:)
    real*8, intent(in), optional :: up2d
    real*8, intent(in), optional :: up2dsp(:,:)
    real*8, intent(in), optional :: up2dcidx(:)
    integer, intent(in), optional :: up2sh
    integer, intent(in), optional :: up2n
    integer, intent(in), optional :: nid0
    integer, intent(in), optional :: id0
    integer, intent(in), optional :: iz0
    logical, intent(in), optional :: nozero

    real*8, parameter :: eps = 1d-10

    real*8 :: x0(3), dist0, rcutshel, rcutn, up2rmax, dmaxenv
    integer :: ireg0(3), ireg(3), idxreg, nats
    integer :: i, j, k
    integer, allocatable :: iord(:), iiord(:), ishell(:), iaux(:,:)
    integer :: nshel
    real*8, allocatable :: rshel(:)
    integer, allocatable :: idxshel(:)
    logical :: doshell, enough

    if (.not.present(up2d).and..not.present(up2dsp).and..not.present(up2dcidx).and.&
        .not.present(up2sh).and..not.present(up2n)) &
       call ferror("list_near_atoms","must give one of up2d, up2dsp, up2dcidx, up2sh, or up2n",faterr)
    doshell = present(ishell0) .or. present(up2sh)

    ! Find the integer region for the main cell copy of the input point
    x0 = xp
    call e%y2z_center(x0,icrd,icrd_cart,lvec)
    ireg0 = e%c2p(x0)
    lvec = nint(e%xr2x(real(lvec,8)))

    ! calculate the maximum distance for which we guarantee we will
    ! get all the atoms
    dmaxenv = max(e%rsph_env - norm2(x0 - e%x0),0d0)

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
    if (present(up2dcidx)) up2rmax = maxval(up2dcidx(1:e%ncell))

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
             if (present(iz0)) then
                if (e%spc(e%at(k)%is)%z /= iz0) cycle
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
             if (present(up2dcidx)) then
                if (dist0 > up2dcidx(e%at(k)%cidx)) cycle
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
       if (present(up2dcidx)) then
          ! Same as above, but with maxval(up2dcidx).
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
          allocate(iord(nat))
          do i = 1, nat
             iord(i) = i
          end do
          if (doshell) then
             ! First by shell (if available) then by distance within each shell
             ! (otherwise the atomeps threshold may screw up the shell order)
             call mergesort(ishell,iord,1,nat)

             allocate(iaux(2,nshel))
             do i = 1, nat
                if (i == 1) then
                   iaux(1,ishell(iord(i))) = 1
                else if (ishell(iord(i)) /= ishell(iord(i-1))) then
                   iaux(2,ishell(iord(i-1))) = i-1
                   iaux(:,ishell(iord(i))) = i
                else if (i == nat) then
                   iaux(2,ishell(iord(i))) = i
                end if
             end do
             do i = 1, nshel
                call mergesort(dist,iord,iaux(1,i),iaux(2,i))
             end do
             deallocate(iaux)
          else
             ! We do not do shell so by distance only
             call mergesort(dist,iord,1,nat)
          end if
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
    
  !> Translate the point x0 to the main cell if the environment is
  !> from a crystal. Then, calculate the core (if zpsp is present) or
  !> promolecular densities at a point x0 (coord format given by icrd)
  !> using atomic radial grids up to a number of derivatives nder
  !> (max: 2). Returns the density (f), gradient (fp, nder >= 1), and
  !> Hessian (fpp, nder >= 2). If a fragment (fr) is given, then only
  !> the atoms in it contribute. This routine is thread-safe.
  module subroutine promolecular(e,x0,icrd,f,fp,fpp,nder,zpsp,fr)
    use grid1mod, only: cgrid, agrid, grid1
    use global, only: cutrad
    use fragmentmod, only: fragment
    use tools_io, only: ferror, faterr
    use param, only: icrd_cart, maxzat
    class(environ), intent(in) :: e
    real*8, intent(in) :: x0(3) !< Point in cryst. coords.
    integer, intent(in) :: icrd !< Input coordinate format
    real*8, intent(out) :: f !< Density
    real*8, intent(out) :: fp(3) !< Density gradient
    real*8, intent(out) :: fpp(3,3) !< Density hessian
    integer, intent(in) :: nder !< Number of derivatives to calculate
    integer, intent(in), optional :: zpsp(:) !< core charges
    type(fragment), intent(in), optional :: fr !< Fragment contributing to the density

    integer :: i, j, k, ii, iz, ierr
    real*8 :: xc(3), xx(3), rlvec(3), r, rinv1, rinv1rp
    real*8 :: rho, rhop, rhopp, rfac, radd, rmax
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

    ! convert to Cartesian and move to the main cell (crystals)
    xc = x0
    call e%y2z_center(xc,icrd,icrd_cart)

    ! calculate the cutoffs
    allocate(rcutmax(e%nspc,2))
    rcutmax = 0d0
    rmax = 0d0
    do i = 1, e%nspc
       iz = e%spc(i)%z
       if (iz == 0 .or. iz > maxzat) cycle
       if (iscore) then
          g => cgrid(iz,zpsp(iz))
       else
          g => agrid(iz)
       end if
       rcutmax(i,2) = min(cutrad(iz),g%rmax)
       rmax = max(rcutmax(i,2),rmax)
    end do
    if (rmax >= e%dmax0) &
       call ferror("promolecular","dmax0 not large enough for promolecular",faterr)

    ! compute the list of atoms that contribute to the point
    call e%list_near_atoms(xc,icrd_cart,.false.,nat,nid,dist,lvec,ierr,up2dsp=rcutmax)
    deallocate(rcutmax)
    rlvec = lvec
    rlvec = e%x2c(rlvec)
    if (nat == 0) return

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
  !> If present, return half the nearest-neighbor distance for each
  !> atom in rnn2, or 0.0 if not found.
  subroutine find_asterisms_covalent(e,nstar,rnn2)
    use global, only: bondfactor, atomeps
    use tools_io, only: ferror, faterr, uout, string
    use types, only: realloc
    use param, only: atmcov
    class(environ), intent(in) :: e
    type(neighstar), allocatable, intent(inout) :: nstar(:)
    real*8, allocatable, intent(inout), optional :: rnn2(:)
    
    integer :: i, j
    real*8 :: x0(3), dist2, ri, rj, rij2, r2
    integer :: p0(3), p1(3), idx1
    integer :: j1, j2, j3, ki, kj, is, js
    real*8, allocatable :: rij2(:,:,:)

    ! return if there are no atoms
    if (e%n == 0 .or. e%nspc == 0) return

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
    
    ! allocate the rnn2 array
    if (present(rnn2)) then
       if (allocated(rnn2)) deallocate(rnn2)
       allocate(rnn2(e%ncell))
       rnn2 = atomeps * atomeps
    end if

    ! pre-calculate the distance^2 matrix
    allocate(rij2(e%nspc,e%nspc,2))
    rij2 = 0d0
    do i = 1, e%nspc
       ri = atmcov(e%spc(i)%z)
       do j = i, e%nspc
          rj = atmcov(e%spc(j)%z)
          
          r2 = (ri+rj) * bondfactor
          rij2(i,j,2) = r2*r2
          rij2(j,i,2) = rij2(i,j,2)

          r2 = (ri+rj) / bondfactor
          rij2(i,j,1) = r2*r2
          rij2(j,i,1) = rij2(i,j,1)
       end do
    end do
    if (sqrt(maxval(rij2)) > e%boxsize) then
       write (uout,'("rij2    = ",A)') string(sqrt(maxval(rij2)),'f',10,4)
       write (uout,'("boxsize = ",A)') string(e%boxsize,'f',10,4)
       call ferror("find_asterisms_covalent","boxsize too small for find_asterisms_covalent",faterr)
    end if

    ! run over atoms in the unit cell and build the connectivity star
    do ki = 1, e%ncell
       p0 = e%c2p(e%at(ki)%r)
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
                      nstar(ki)%lcon(:,nstar(ki)%ncon) = e%at(kj)%lvec- e%at(ki)%lvec
                      if (present(rnn2)) then
                         if (rnn2(ki) < 1d-10 .or. dist2 < rnn2(ki)) then
                            rnn2(ki) = dist2
                         end if
                      end if
                   end if
                end do ! j = e%nrlo(idx1), e%nrhi(idx1)

             end do ! j3 = -1, 1
          end do ! j2 = -1, 1
       end do ! j1 = -1, 1
    end do ! ki = 1, e%ncell

    ! reallocate the asterism arrays and convert to half nearest
    ! neighbor distance
    do i = 1, e%ncell
       if (nstar(i)%ncon > 0) then
          call realloc(nstar(i)%idcon,nstar(i)%ncon)
          call realloc(nstar(i)%lcon,3,nstar(i)%ncon)
       end if
    end do
    if (present(rnn2)) then
       rnn2 = sqrt(rnn2)
    end if

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
       write (uout,'("  Maximum interaction distance (",A,"): ",A)') iunitname0(iunit), string(e%dmax0*dunit0(iunit),'f',8,4)
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

  !> Check the consistency of an environment. For debug
  module subroutine environ_check(e)
    use tools_io, only: ferror, faterr
    use param, only: ctsq32, atmcov
    class(environ), intent(in) :: e

    integer :: i, ireg(3), i1, i2, i3, imax, i3min
    real*8 :: mcov, x0(3), x1(3), dist, rcut0

    ! all atoms in the covered regions
    do i = 1, e%n
       ireg = e%c2p(e%at(i)%r)
       if (any(ireg < e%nmin) .or. any(ireg > e%nmax)) then
          call ferror("environ_check","atoms outside covered regions",faterr)
       end if
    end do

    ! calculate the boxsize, maximum bondfactor of 2.0
    mcov = 0d0
    do i = 1, e%nspc
       if (e%spc(i)%z > 0) then
          mcov = max(mcov,atmcov(e%spc(i)%z))
       end if
    end do
    if (e%boxsize < 4d0 * mcov) then
       call ferror("environ_check","boxsize too small",faterr)
    end if

    ! all atoms not in the search region have rs_rcut(i) > rsph_env
    imax = e%rs_imax + 1
    do i1 = -imax, imax
       do i2 = -imax, imax
          if (abs(i1) == imax .or. abs(i2) == imax) then
             i3min = 0
          else
             i3min = imax
          end if
          do i3 = i3min, imax
             x0 = max( ((/i1,i2,i3/) - 0.5d0),0d0)
             x1 = max(-((/i1,i2,i3/) + 0.5d0),0d0)
             dist = sqrt(x0(1)*x0(1)+x0(2)*x0(2)+x0(3)*x0(3)+x1(1)*x1(1)+x1(2)*x1(2)+x1(3)*x1(3))
             rcut0 = max(dist - ctsq32,0d0) * e%boxsize
             if (rcut0 <= e%rsph_env) then
                call ferror("environ_check","one instance of rs_rcut(i) <= rsph_env",faterr)
             end if
          end do
       end do
    end do
    
    write (*,*) "all tests passed!"
    stop 1

  end subroutine environ_check

  !xx! private procedures

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
