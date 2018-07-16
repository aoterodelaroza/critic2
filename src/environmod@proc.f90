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

    if (allocated(e%at)) deallocate(e%at)
    if (allocated(e%imap)) deallocate(e%imap)
    if (allocated(e%nrlo)) deallocate(e%nrlo)
    if (allocated(e%nrhi)) deallocate(e%nrhi)
    e%n = 0
    e%ncell = 0
    e%nregc = 0
    e%nreg = 0
    e%nmin = 0
    e%nmax = 0
    e%nregion = 0

  end subroutine environ_end
  
  !> Build an environment from molecule data
  module subroutine environ_build_from_molecule(e,n,at)
    use types, only: realloc
    class(environ), intent(inout) :: e
    integer, intent(in) :: n
    type(celatom), intent(in) :: at(n)

    integer :: i 

    e%n = n
    if (allocated(e%at)) deallocate(e%at)
    allocate(e%at(e%n))
    do i = 1, n
       e%at(i)%x = at(i)%x
       e%at(i)%x = e%at(i)%x - nint(e%at(i)%x)
       e%at(i)%r = at(i)%r
       e%at(i)%idx = at(i)%idx
       e%at(i)%cidx = i
       e%at(i)%ir = at(i)%ir
       e%at(i)%ic = at(i)%ic
       e%at(i)%lvec = 0
       e%at(i)%lenv = 0
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
  module subroutine environ_build_from_crystal(e,nspc,spc,n,at,m_xr2c,m_x2xr,dmax0)
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
    real*8, intent(in), optional :: dmax0

    logical :: dorepeat
    real*8 :: sphmax, dmax, m_xr2x(3,3), x(3), xc(3), dd
    integer :: i1, i2, i3, i, imax, p(3), px(3)

    ! prepare
    m_xr2x = matinv(m_x2xr)

    ! calculate the maximum diagonal half-length (sphmax)
    sphmax = norm2(xr2c((/0d0,0d0,0d0/) - (/0.5d0,0.5d0,0.5d0/)))
    sphmax = max(sphmax,norm2(xr2c((/1d0,0d0,0d0/) - (/0.5d0,0.5d0,0.5d0/))))
    sphmax = max(sphmax,norm2(xr2c((/0d0,1d0,0d0/) - (/0.5d0,0.5d0,0.5d0/))))
    sphmax = max(sphmax,norm2(xr2c((/0d0,0d0,1d0/) - (/0.5d0,0.5d0,0.5d0/))))
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
       e%at(i)%x = x2xr(at(i)%x)
       e%at(i)%x = e%at(i)%x - nint(e%at(i)%x)
       x = xr2x(e%at(i)%x)
       e%at(i)%r = xr2c(e%at(i)%x)
       e%at(i)%idx = at(i)%idx
       e%at(i)%cidx = i
       e%at(i)%ir = at(i)%ir
       e%at(i)%ic = at(i)%ic
       e%at(i)%lenv = nint(x - at(i)%x)
       e%at(i)%lvec = at(i)%lvec + e%at(i)%lenv
       e%at(i)%is = at(i)%is
    end do
    e%ncell = n

    ! include all atoms at a distance sphmax+dmax from the origin
    dorepeat = .true.
    imax = 0
    do while (dorepeat)
       dorepeat = .false.
       imax = imax + 1
       do i1 = -imax, imax
          do i2 = -imax, imax
             do i3 = -imax, imax
                if (abs(i1) /= imax .and. abs(i2) /= imax .and. abs(i3) /= imax) cycle

                do i = 1, n
                   p = (/i1, i2, i3/)
                   px = nint(xr2x(real(p,8)))

                   x = e%at(i)%x + p
                   xc = xr2c(x)
                   dd = norm2(xc)
                   if (dd < sphmax+dmax) then
                      dorepeat = .true.
                      e%n = e%n + 1
                      if (e%n > size(e%at,1)) call realloc(e%at,2*e%n)

                      e%at(e%n)%x = x
                      e%at(e%n)%r = xc
                      e%at(e%n)%idx = e%at(i)%idx
                      e%at(e%n)%cidx = e%at(i)%cidx
                      e%at(e%n)%ir = e%at(i)%ir
                      e%at(e%n)%ic = e%at(i)%ic
                      e%at(e%n)%lenv = e%at(i)%lenv + px
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
    pure function xr2c(xx) result(res)
      real*8, intent(in)  :: xx(3)
      real*8 :: res(3)
      res = matmul(m_xr2c,xx)
    end function xr2c
    pure function x2xr(xx) result(res)
      real*8, intent(in)  :: xx(3)
      real*8 :: res(3)
      res = matmul(m_x2xr,xx)
    end function x2xr
    pure function xr2x(xx) result(res)
      real*8, intent(in)  :: xx(3)
      real*8 :: res(3)
      res = matmul(m_xr2x,xx)
    end function xr2x
  end subroutine environ_build_from_crystal

  !> Cartesian to region transform
  pure module function c2p(e,xx) result(res)
    class(environ), intent(in) :: e
    real*8, intent(in)  :: xx(3)
    integer :: res(3)

    res = floor((max(min(xx,e%xmax),e%xmin) - e%x0) / e%boxsize)

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
    integer :: i1, i2, i3, imax, iadd, nreg
    real*8 :: x0(3), x1(3), dist, rcut0
    integer, allocatable :: iadd(:)
    real*8, allocatable :: rcut(:)

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
    allocate(iadd(e%nregion),rcut(e%nregion))
    dorepeat = .true.
    imax = 0
    nreg = 1
    iadd(1) = 0
    rcut(1) = 0d0
    do while (dorepeat)
       dorepeat = .false.
       imax = imax + 1
       do i1 = -imax, imax
          do i2 = -imax, imax
             do i3 = -imax, imax
                if (abs(i1) /= imax .and. abs(i2) /= imax .and. abs(i3) /= imax) cycle
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
                   if (nreg > size(iadd,1)) then
                      call realloc(iadd,2*nreg)
                      call realloc(rcut,2*nreg)
                   end if
                   ! because the p -> i transformation is linear, we can write the offset as a single
                   ! integer
                   iadd(nreg) = i1 + e%nreg(1) * (i2 + e%nreg(2) * i3)
                   rcut(nreg) = rcut0
                end if
             end do
          end do
       end do
    end do
    call realloc(iadd,nreg)
    call realloc(rcut,nreg)

    ! sort and populate the environ arrays
    allocate(iord(nreg))
    do i = 1, nreg
       iord(i) = i
    end do
    call qcksort(rcut,iord,1,nreg)
    allocate(e%iaddregs(nreg),e%rcutregs(nreg))
    do i = 1, nreg
       e%iaddregs(i) = iadd(iord(i))
       e%rcutregs(i) = rcut(iord(i))
    end do
    e%nregs = nreg
    deallocate(iord,iadd,rcut)

  end subroutine calculate_regions

end submodule proc
