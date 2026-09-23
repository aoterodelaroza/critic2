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

!> One-dimensional grid class
submodule (grid1mod) proc
  implicit none

  !xx! private procedures
  ! subroutine read_fit_block(file,key,n,np,al,co,found,navail,ti)
  ! subroutine tabulate(g,z,np,al,co)

  ! logarithmic grid for the tabulated densities, r_i = exp(xmin + (i-1)*dx) / Z
  real*8, parameter :: tab_xmin = -10d0
  real*8, parameter :: tab_dx = 0.005d0
  real*8, parameter :: tab_rcap = 200d0 !< Never tabulate beyond this radius (bohr)

  ! cutoffs
  real*8, parameter :: core_cutdens = 1d-08 !< Cutoff contribution for core radial grids

contains

  !> Deallocate arrays and uninitialize
  module subroutine grid1_end(g)
    class(grid1), intent(inout) :: g

    g%isinit = .false.
    if (allocated(g%r)) deallocate(g%r)
    if (allocated(g%f)) deallocate(g%f)
    if (allocated(g%fp)) deallocate(g%fp)
    if (allocated(g%fpp)) deallocate(g%fpp)

  end subroutine grid1_end

  !> Build the radial density of atom z from the fitted analytical
  !> densities in the database (critic_home/atomdens/fit_ZZZ_Sym.dat).
  !> With q = 0, this is the all-electron density of the neutral atom
  !> (block STATE N=z). With 0 < q < z, it is the frozen-core density
  !> for a pseudopotential with q valence electrons (block CORE
  !> N=z-q). The fit, rho(r) = sum_i c_i
  !> r^n_i exp(-alpha_i r), is tabulated with its exact derivatives on
  !> a logarithmic grid. If the density is not available, the grid is
  !> left uninitialized and a warning is issued, or the reason is
  !> returned in errmsg if present.
  module subroutine read_db(g,z,q,ti,errmsg)
    use global, only: critic_home
    use tools_io, only: nameguess, ferror, warning, string
    use param, only: dirsep
    class(grid1), intent(inout) :: g !< Output radial grid
    integer, intent(in) :: z !< Atomic number
    integer, intent(in) :: q !< Atomic pseudopotential charge (0 = all-electron)
    type(thread_info), intent(in), optional :: ti
    character(len=:), allocatable, intent(out), optional :: errmsg

    character(len=:), allocatable :: file, msg
    integer, allocatable :: np(:), navail(:)
    real*8, allocatable :: al(:), co(:)
    logical :: found
    integer :: i

    call g%grid1_end()
    if (present(errmsg)) errmsg = ""

    ! build the file name and read the block
    file = trim(critic_home) // dirsep // "atomdens" // dirsep // "fit_" //&
       string(z,3,pad0=.true.) // "_" // trim(nameguess(z,.true.)) // ".dat"
    if (q == 0) then
       call read_fit_block(file,"STATE",z,np,al,co,found,navail,ti)
       msg = 'Atomic density for Z = ' // string(z) // ' not found in: ' // file
    else
       call read_fit_block(file,"CORE",z-q,np,al,co,found,navail,ti)
       msg = 'No core density with ' // string(z-q) // ' electrons (ZPSP = ' // string(q) //&
          ') for Z = ' // string(z) // '. Available ZPSP:'
       do i = 1, size(navail)
          msg = msg // " " // string(z - navail(i))
       end do
    end if
    if (found) found = (size(np) > 0)
    if (.not.found) then
       if (present(errmsg)) then
          errmsg = msg
       else
          call ferror('read_db',msg,warning)
       end if
       return
    end if

    ! tabulate
    call tabulate(g,z,np,al,co)
    g%z = z
    g%qat = q

  end subroutine read_db

  !> Interpolate the radial grid g at distance r0, and obtain the value,
  !> first derivative and second derivative.
  module subroutine interp(g,r0,f,fp,fpp)
    class(grid1), intent(in) :: g !< The radial grid.
    real*8, intent(in) :: r0 !< Value of the radial coordinate.
    real*8, intent(out) :: f !< Interpolated value
    real*8, intent(out) :: fp !< Interpolated first derivative
    real*8, intent(out) :: fpp !< Interpolated second derivative

    integer :: ir, i, j, ii
    real*8 :: r, prod, rr(4), dr1(4), x1dr12(4,4)

    f = 0d0
    fp = 0d0
    fpp = 0d0

    if (.not.g%isinit) return
    if (r0 >= g%rmax) return

    ! careful with grid limits.
    if (r0 <= g%r(1)) then
       ir = 1
       r = g%r(1)
    else
       ir = 1 + floor(log(r0/g%a)/g%b)
       r = r0
    end if

    x1dr12 = 0d0
    do i = 1, 4
       ii = min(max(ir,2),g%ngrid-2) - 2 + i
       rr(i) = g%r(ii)
       dr1(i) = r - rr(i)
       do j = 1, i-1
          x1dr12(i,j) = 1d0 / (rr(i) - rr(j))
          x1dr12(j,i) = -x1dr12(i,j)
       end do
    end do

    ! interpolate, lagrange 3rd order, 4 nodes
    do i = 1, 4
       ii = min(max(ir,2),g%ngrid-2) - 2 + i
       prod = 1.d0
       do j = 1 ,4
          if (i == j) cycle
          prod = prod * dr1(j) * x1dr12(i,j)
       end do
       f = f + g%f(ii) * prod
       fp = fp + g%fp(ii) * prod
       fpp = fpp + g%fpp(ii) * prod
    end do

  end subroutine interp

  !> Read the core density from the internal density tables for atom
  !> with Z = iz and ZPSP = iq. A pseudopotential with all electrons
  !> in valence (iq = iz) has no core, and its grid is left empty.
  !> If the core is not available, the reason is returned in errmsg
  !> (empty otherwise).
  module subroutine grid1_register_core(iz,iq,errmsg)
    use param, only: maxzat0, maxzat
    integer, intent(in) :: iz, iq
    character(len=:), allocatable, intent(out) :: errmsg

    errmsg = ""
    if (.not.allocated(cgrid)) allocate(cgrid(maxzat0,maxzat0))
    if (iz <= 0 .or. iz > maxzat) return
    if (iq <= 0 .or. iq >= iz) return
    if (cgrid(iz,iq)%isinit) return
    call cgrid(iz,iq)%read_db(iz,iq,errmsg=errmsg)

  end subroutine grid1_register_core

  !> Read the all-electron density from the internal density tables
  !> for atom with Z = iz.
  module subroutine grid1_register_ae(iz)
    use param, only: maxzat0, maxzat
    integer, intent(in) :: iz

    integer :: i

    if (iz <= 0 .or. iz > maxzat) return
    if (.not.allocated(agrid)) then
       allocate(agrid(maxzat0))
       do i = 1, maxzat0
          agrid(i)%isinit = .false.
          agrid(i)%z = 0
          agrid(i)%qat = 0
       end do
    end if

    if (agrid(iz)%isinit) then
       if (agrid(iz)%z == iz) return
    end if

    call agrid(iz)%read_db(iz,0)

  end subroutine grid1_register_ae

  !> Deallocate the core and all-electron density grid
  module subroutine grid1_clean_grids()

    if (allocated(agrid)) deallocate(agrid)
    if (allocated(cgrid)) deallocate(cgrid)

  end subroutine grid1_clean_grids

  !xx! private procedures

  !> Read one block of an analytical density file (fit_*.dat). The
  !> block has the header "STATE <N> <q> <nterm>
  !> <sym>" (key = STATE) or "CORE <N> <nterm> <sym>" (key = CORE),
  !> followed by nterm lines "n alpha c". Returns the powers (np),
  !> exponents (al), and coefficients (co) of the block with N = n,
  !> dropping the terms with c = 0. found is false if the file or the
  !> block do not exist. navail is the list of N of the blocks in the
  !> file, for error messages.
  subroutine read_fit_block(file,key,n,np,al,co,found,navail,ti)
    use tools_io, only: fopen_read, fclose, getline_raw, getword, isinteger, isreal, &
       equal
    use types, only: realloc
    character*(*), intent(in) :: file
    character*(*), intent(in) :: key
    integer, intent(in) :: n
    integer, allocatable, intent(inout) :: np(:)
    real*8, allocatable, intent(inout) :: al(:), co(:)
    logical, intent(out) :: found
    integer, allocatable, intent(out) :: navail(:)
    type(thread_info), intent(in), optional :: ti

    integer :: lu, lp, i, nblock, nterm, nn, nav
    character(len=:), allocatable :: line, word
    real*8 :: rn, ra, rc
    logical :: exist, ok

    found = .false.
    nav = 0
    allocate(navail(10))
    inquire(file=file,exist=exist)
    if (exist) then
       lu = fopen_read(file,ti=ti)
       do while (getline_raw(lu,line,.false.))
          lp = 1
          word = getword(line,lp)
          if (.not.equal(word,key)) cycle

          ! block header: N, (q,) nterm
          ok = isinteger(nblock,line,lp)
          if (key == "STATE") ok = ok .and. isinteger(i,line,lp)
          ok = ok .and. isinteger(nterm,line,lp)
          if (.not.ok) cycle
          nav = nav + 1
          if (nav > size(navail)) call realloc(navail,2*nav)
          navail(nav) = nblock
          if (nblock /= n .or. found) cycle

          ! read the terms, skip zero coefficients
          if (allocated(np)) deallocate(np)
          if (allocated(al)) deallocate(al)
          if (allocated(co)) deallocate(co)
          allocate(np(nterm),al(nterm),co(nterm))
          nn = 0
          do i = 1, nterm
             ok = getline_raw(lu,line,.true.)
             lp = 1
             ok = isreal(rn,line,lp)
             ok = ok .and. isreal(ra,line,lp)
             ok = ok .and. isreal(rc,line,lp)
             if (.not.ok) exit
             if (rc == 0d0) cycle
             nn = nn + 1
             np(nn) = nint(rn)
             al(nn) = ra
             co(nn) = rc
          end do
          if (.not.ok) exit
          np = np(1:nn)
          al = al(1:nn)
          co = co(1:nn)
          found = .true.
       end do
       call fclose(lu)
    end if
    call realloc(navail,nav)

  end subroutine read_fit_block

  !> Tabulate the analytical density rho(r) = sum_i co_i r^np_i
  !> exp(-al_i r) and its exact first and second derivatives on the
  !> logarithmic grid r_i = exp(tab_xmin + (i-1)*tab_dx) / z, up to the
  !> first node where the density drops below core_cutdens (and at
  !> least four nodes, the interpolation stencil).
  subroutine tabulate(g,z,np,al,co)
    use types, only: realloc
    class(grid1), intent(inout) :: g
    integer, intent(in) :: z
    integer, intent(in) :: np(:)
    real*8, intent(in) :: al(:), co(:)

    integer :: i, k, nmax
    real*8 :: r, t, u

    g%a = exp(tab_xmin) / real(z,8)
    g%b = tab_dx
    nmax = ceiling((log(tab_rcap * z) - tab_xmin) / tab_dx) + 1
    allocate(g%r(nmax),g%f(nmax),g%fp(nmax),g%fpp(nmax))

    ! d/dr r^n e^(-ar) = r^n e^(-ar) (n/r - a)
    ! d2/dr2 r^n e^(-ar) = r^n e^(-ar) ((n/r - a)^2 - n/r^2)
    do i = 1, nmax
       r = g%a * exp(g%b * (i-1))
       g%r(i) = r
       g%f(i) = 0d0
       g%fp(i) = 0d0
       g%fpp(i) = 0d0
       do k = 1, size(np)
          t = co(k) * r**np(k) * exp(-al(k) * r)
          u = np(k) / r - al(k)
          g%f(i) = g%f(i) + t
          g%fp(i) = g%fp(i) + t * u
          g%fpp(i) = g%fpp(i) + t * (u * u - np(k) / (r * r))
       end do
       if (g%f(i) < core_cutdens .and. i >= 4) exit
    end do
    g%ngrid = min(i,nmax)
    call realloc(g%r,g%ngrid)
    call realloc(g%f,g%ngrid)
    call realloc(g%fp,g%ngrid)
    call realloc(g%fpp,g%ngrid)
    g%rmax = g%r(g%ngrid)
    g%rmax2 = g%rmax * g%rmax
    g%isinit = .true.

  end subroutine tabulate

end submodule proc
