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

!> One-dimensional grid class
submodule (grid1mod) proc
  implicit none

  !xx! private procedures
  ! subroutine read_critic(g,file,n,abspath,ti)

  ! radial grid derivation formulas
  integer, parameter :: noef(6,3) = reshape((/&
     0,  1,  2,  3,  4,  5,&
     -2, -1,  0,  1,  2,  3,&
     -5, -4, -3, -2, -1,  0/),shape(noef)) !< Node offsets for 6-point derivation formulas.
  real*8, parameter :: coef1(6,3) = reshape((/&
     -274,  600, -600,  400,  -150,  24,&
     6,  -60,  -40,  120,   -30,   4,&
     -24, 150, -400,  600,  -600, 274 /),shape(coef1)) !< Coefficients for first derivative.
  real*8, parameter :: coef2(6,3) = reshape((/&
     225, -770,  1070,  -780,   305,   -50,&
     -5,   80,  -150,    80,    -5,     0,&
     -50,  305,  -780 ,  1070,  -770,   225/),shape(coef2)) !< Coefficients for second derivative.
  real*8, parameter :: fac1=1d0/120d0 !< Prefactor for first derivative.
  real*8, parameter :: fac2=2d0/120d0 !< Prefactor for second derivative.

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

  !> Read a density core file from the database.
  module subroutine read_db(g,z,q,ti)
    use global, only: critic_home
    use tools_io, only: nameguess, warning, lower, ferror
    use param, only: dirsep
    class(grid1), intent(inout) :: g !< Output radial grid
    integer, intent(in) :: z !< Atomic number
    integer, intent(in) :: q !< Atomic pseudopotential charge
    type(thread_info), intent(in), optional :: ti

    character(len=:), allocatable :: file

    ! build the file name
    file = trim(critic_home) // dirsep // "wfc" // dirsep // lower(nameguess(z)) // "_pbe.wfc"

    ! anions not supported
    if (q < 0) &
       call ferror('read_db','Anions not supported: neutral atomic density used instead',warning)

    ! do it
    call read_critic(g,file,z-q,.true.,ti=ti)
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
    if (g%z - g%qat <= 0) return

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
  !> with Z = iz and ZPSP = iq.
  module subroutine grid1_register_core(iz,iq)
    use param, only: maxzat0, maxzat
    integer, intent(in) :: iz, iq

    integer :: i, j

    if (iz <= 0 .or. iz > maxzat) return
    if (iq <= 0 .or. iq > iz) return

    if (.not.allocated(cgrid)) then
       allocate(cgrid(maxzat0,maxzat0))
       do i = 1, maxzat0
          do j = 1, maxzat0
             cgrid(i,j)%z = 0
             cgrid(i,j)%qat = 0
          end do
       end do
    end if

    if (cgrid(iz,iq)%isinit) then
       if(cgrid(iz,iq)%z == iz .and. cgrid(iz,iq)%qat == iq) return
    end if
    call cgrid(iz,iq)%read_db(iz,iq)

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

  !> Read grid in critic format. This format is adapted from the wfc files
  !> of the ld1 program in the quantum espresso distribution. Only n electrons
  !> out of the total Z are used to build the grid.
  subroutine read_critic(g,file,n,abspath,ti)
    use types, only: realloc
    use tools_io, only: uout, warning, ferror, fopen_read, string, fclose
    use param, only: pi
    class(grid1), intent(inout) :: g !< One-dimensional grid on output
    character*(*), intent(in) :: file !< File with the grid description
    integer, intent(in) :: n !< Number of electrons read
    logical, intent(in) :: abspath !< Absolute path?
    type(thread_info), intent(in), optional :: ti

    integer :: i, j, lu, nn
    integer, allocatable :: occ(:)
    real*8, allocatable :: wfcin(:), rr(:,:)
    character*2, allocatable :: wfcl(:)
    real*8 :: xmin, zz, dx, r, r1, r2, r3 ,r4, delta, delta2
    integer :: ngrid, ns, ic
    logical :: exist

    ! check that the file exists
    inquire(file=file,exist=exist)
    if (.not.exist) then
       write (uout,'("File: ",A)') trim(file)
       call ferror("grid1_read_critic","Atomic density file not found",warning)
       g%isinit = .false.
       return
    end if

    ! Read header and allocate arrays
    lu = fopen_read(file,abspath0=abspath,ti=ti)
    read (lu,*) nn

    allocate(wfcin(nn),wfcl(nn),occ(nn))
    occ = 0
    read (lu,*) (wfcl(i),i=1,nn)
    read (lu,*) (occ(i),i=1,nn)
    read (lu,*) xmin, zz, dx, ngrid

    if (sum(occ) /= n) then
       ns = 0
       do i = 1, nn
          if (ns + occ(i) > n) then
             occ(i) = n - ns
             occ(i+1:nn) = 0
             exit
          else
             ns = ns + occ(i)
          end if
       end do
    end if

    ! Read the grid and build the density
    allocate(g%r(ngrid),rr(ngrid,0:2))
    rr = 0d0
    do i = 1, ngrid
       read (lu,*) r, (wfcin(j),j=1,nn)
       g%r(i) = r
       rr(i,0) = dot_product(occ(1:nn),wfcin(1:nn)**2)
       if (rr(i,0)/(4d0*pi*r**2) < core_cutdens .and. i > 1) then
          ngrid = i
          call realloc(g%r,ngrid)
          exit
       end if
    end do

    ! ! laplacian transformation
    ! ! 1/r^2 * d/dr ( r^2 * df/dr)
    ! do i = 1, ngrid
    !    if (i <= 2) then
    !       ic = 1
    !    else if (i >= ngrid-2) then
    !       ic = 3
    !    else
    !       ic = 2
    !    end if
    !    do j = 1, 6
    !       rr(i,1) = rr(i,1) + coef1(j,ic) * rr(i+noef(j,ic),0)
    !       rr(i,2) = rr(i,2) + coef2(j,ic) * rr(i+noef(j,ic),0)
    !    end do
    !    rr(i,1) = rr(i,1) * fac1
    !    rr(i,2) = rr(i,2) * fac2

    !    r = g%r(i)
    !    r1 = 1d0 / r
    !    r2 = r1 * r1
    !    r3 = r2 * r1
    !    r4 = r3 * r1
    !    delta=1.d0/g%b
    !    delta2=delta*delta

    !    g%f(i) = rr(i,0) * r2
    !    g%fp(i) = (rr(i,1) * delta - 2.d0 * rr(i,0)) * r3
    !    g%fpp(i) = (rr(i,2) * delta2 - 5.d0 * rr(i,1) * delta + 6.d0 * rr(i,0)) * r4
    ! end do

    ! fill grid info
    g%isinit = .true.
    g%a = exp(xmin) / zz
    g%b = dx
    g%ngrid = ngrid
    g%rmax = g%r(ngrid)
    g%rmax2 = g%r(ngrid)**2

    ! calculate derivatives
    allocate(g%f(ngrid),g%fp(ngrid),g%fpp(ngrid))
    do i = 1, ngrid
       if (i <= 2) then
          ic = 1
       else if (i >= ngrid-2) then
          ic = 3
       else
          ic = 2
       end if
       do j = 1, 6
          rr(i,1) = rr(i,1) + coef1(j,ic) * rr(i+noef(j,ic),0)
          rr(i,2) = rr(i,2) + coef2(j,ic) * rr(i+noef(j,ic),0)
       end do
       rr(i,1) = rr(i,1) * fac1
       rr(i,2) = rr(i,2) * fac2

       r = g%r(i)
       r1 = 1d0 / r
       r2 = r1 * r1
       r3 = r2 * r1
       r4 = r3 * r1
       delta=1.d0/g%b
       delta2=delta*delta

       g%f(i) = rr(i,0) * r2
       g%fp(i) = (rr(i,1) * delta - 2.d0 * rr(i,0)) * r3
       g%fpp(i) = (rr(i,2) * delta2 - 5.d0 * rr(i,1) * delta + 6.d0 * rr(i,0)) * r4
    end do
    g%f = g%f / (4d0*pi)
    g%fp = g%fp / (4d0*pi)
    g%fpp = g%fpp / (4d0*pi)

    ! close the density file
    call fclose(lu)

    ! ! check the normalization of the density file and output
    ! if (verbose) then
    !    econf = ""
    !    do i = 1, nn
    !       econf = econf // string(wfcl(i),2)
    !       econf = econf // "(" // string(occ(i)) // ")"
    !    end do

    !    write (uout,'("+ Read density file: ", A)') string(file)
    !    write (uout,'("  Log grid (r = a*e^(b*x)) with a = ",A,", b = ",A)') &
    !       string(g%a,'e',length=10,decimal=4), string(g%b,'e',length=10,decimal=4)
    !    write (uout,'("  Num. grid points = ",A,", rmax (bohr) = ",A)') &
    !       string(g%ngrid), string(g%rmax,'f',decimal=7)
    !    write (uout,'("  Integrated charge = ",A)') &
    !       string(sum(g%f * g%r**3 * g%b * 4d0 * pi),'f',decimal=10)
    !    write (uout,'("  El. conf.: ",A)') string(econf)
    ! end if

    ! cleanup
    deallocate(rr,wfcin,wfcl,occ)

  end subroutine read_critic

end submodule proc
