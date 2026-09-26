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

! meshmod: becke-style meshes for molecular integrals.
submodule (meshmod) proc
  implicit none

  !xx! private procedures
  ! function table_filename(lvl) result(file)
  ! subroutine load_table(lvl)
  ! function find_table_entry(lvl,iz,zeff) result(ie)
  ! subroutine atomic_grid_spec(iz,zeff,type,lvl,r,wr,nang,fromtable)
  ! subroutine rmesh_power(rmin,rmax,n,r,wintr)
  ! subroutine partition_becke(m)
  ! subroutine partition_promolecular(m,c)
  ! subroutine rmesh_postg(n,iz,r,wintr)
  ! subroutine rmesh_franchini(n,iz,r,wintr)
  ! function z2nr(z,lvl) result(nr)
  ! function z2nang(z,lvl) result(nang)

  !> One entry (element/effective charge) of a tuned atomic grid table
  type table_entry
     integer :: z = 0 !< atomic number
     integer :: zeff = 0 !< effective charge (electrons treated explicitly)
     real*8 :: rmin = 0d0 !< first radius of the power transform (bohr)
     real*8 :: rmax = 0d0 !< last radius of the power transform (bohr)
     integer, allocatable :: nang(:) !< (number of shells) Lebedev order on each shell
  end type table_entry

  !> A tuned atomic grid table (dat/meshes/tv-13.7-N.txt), one per
  !> level. e is allocated only if the file was found and parsed.
  type table
     logical :: isinit = .false. !< tried to load the file
     type(table_entry), allocatable :: e(:) !< entries
  end type table

  !> Accuracy tier (number of digits) of the table used for each mesh level
  integer, parameter :: level2tier(mesh_level_small:mesh_level_amazing) = (/3,4,5,6,7/)

  !> Tuned atomic grid tables, loaded on first use
  type(table), save :: tab(mesh_level_small:mesh_level_amazing)

contains

  !> Deallocate and uninitalize
  pure module subroutine endmesh(m)
    class(mesh), intent(inout) :: m

    m%n = 0
    m%nat = 0
    if (allocated(m%w)) deallocate(m%w)
    if (allocated(m%x)) deallocate(m%x)
    if (allocated(m%f)) deallocate(m%f)
    if (allocated(m%idat)) deallocate(m%idat)
    if (allocated(m%zat)) deallocate(m%zat)
    if (allocated(m%zeffat)) deallocate(m%zeffat)
    if (allocated(m%fromtable)) deallocate(m%fromtable)
    if (allocated(m%xat)) deallocate(m%xat)
    if (allocated(m%ishoff)) deallocate(m%ishoff)
    if (allocated(m%rsh)) deallocate(m%rsh)
    if (allocated(m%wrsh)) deallocate(m%wrsh)
    if (allocated(m%ipsh)) deallocate(m%ipsh)
    if (allocated(m%wa)) deallocate(m%wa)
    if (allocated(m%wp)) deallocate(m%wp)

  end subroutine endmesh

  !> Generate a molecular/crystal integration mesh for crystal c. The
  !> mesh is the union of atomic grids multiplied by a partition
  !> function of the given type (mesh_type_becke, only for molecules,
  !> or mesh_type_franchini, promolecular weights, always used for
  !> crystals). lvl selects the accuracy (mesh_level_*, default:
  !> good). The atomic grids are taken from the tuned tables in
  !> dat/meshes for the element and effective charge (zpsp, per
  !> species; <= 0 or absent means all-electron); atoms not covered
  !> by the tables use the legacy grids (constant Lebedev order,
  !> postg or Franchini radial grids).
  module subroutine genmesh(m,c,type,lvl,zpsp)
    use crystalmod, only: crystal
    use tools_math, only: select_lebedev
    use types, only: realloc
    use param, only: maxzat, fourpi
    class(mesh), intent(inout) :: m
    type(crystal), intent(inout) :: c
    integer, intent(in), optional :: type
    integer, intent(in), optional :: lvl
    integer, intent(in), optional :: zpsp(:)

    integer :: tmesh, lmesh
    integer :: i, k, iz, is, zeff, nsh, np, ish, il, ip, nang, mang, nr
    real*8 :: xnuc(3)
    real*8, allocatable :: r(:), wr(:), xang(:), yang(:), zang(:), wang(:)
    integer, allocatable :: nangr(:)
    logical :: ft

    if (.not.c%ismolecule) then
       tmesh = mesh_type_franchini
    else if (present(type)) then
       tmesh = type
    else
       tmesh = mesh_type_becke
    end if
    if (present(lvl)) then
       lmesh = lvl
    else
       lmesh = mesh_level_good
    end if

    ! reset the arrays
    call m%end()
    m%type = tmesh
    m%lvl = lmesh

    ! atoms with a grid
    m%nat = 0
    do i = 1, c%ncel
       iz = c%spc(c%atcel(i)%is)%z
       if (iz < 1 .or. iz > maxzat) cycle
       m%nat = m%nat + 1
    end do
    allocate(m%idat(m%nat),m%zat(m%nat),m%zeffat(m%nat),m%fromtable(m%nat))
    allocate(m%xat(3,m%nat),m%ishoff(m%nat+1))
    allocate(m%rsh(100),m%wrsh(100),m%ipsh(101))

    ! atomic grid specifications and offsets (serial: loads the tables)
    nsh = 0
    np = 0
    k = 0
    do i = 1, c%ncel
       is = c%atcel(i)%is
       iz = c%spc(is)%z
       if (iz < 1 .or. iz > maxzat) cycle
       k = k + 1
       zeff = iz
       if (present(zpsp)) then
          if (is <= size(zpsp)) then
             if (zpsp(is) > 0) zeff = zpsp(is)
          end if
       end if
       if (c%ismolecule) then
          xnuc = c%atcel(i)%r
       else
          xnuc = c%x2xr(c%atcel(i)%x)
          xnuc = xnuc - floor(xnuc)
          xnuc = c%xr2c(xnuc)
       end if
       m%idat(k) = i
       m%zat(k) = iz
       m%zeffat(k) = zeff
       m%xat(:,k) = xnuc

       call atomic_grid_spec(iz,zeff,tmesh,lmesh,r,wr,nangr,ft)
       m%fromtable(k) = ft
       nr = size(r)
       if (nsh + nr > size(m%rsh)) then
          call realloc(m%rsh,2*(nsh+nr))
          call realloc(m%wrsh,2*(nsh+nr))
          call realloc(m%ipsh,2*(nsh+nr)+1)
       end if
       m%ishoff(k) = nsh + 1
       do ish = 1, nr
          nsh = nsh + 1
          m%rsh(nsh) = r(ish)
          m%wrsh(nsh) = wr(ish)
          m%ipsh(nsh) = np + 1
          np = np + nangr(ish)
       end do
    end do
    m%ishoff(m%nat+1) = nsh + 1
    m%ipsh(nsh+1) = np + 1
    call realloc(m%rsh,nsh)
    call realloc(m%wrsh,nsh)
    call realloc(m%ipsh,nsh+1)
    m%n = np
    allocate(m%x(3,m%n),m%w(m%n),m%wa(m%n),m%wp(m%n))
    if (m%n == 0) return

    ! points and atomic weights, each atom writes into its own slice
    mang = maxval(m%ipsh(2:nsh+1) - m%ipsh(1:nsh))
    !$omp parallel do private(xang,yang,zang,wang,ish,nang,il,ip) schedule(dynamic)
    do k = 1, m%nat
       allocate(xang(mang),yang(mang),zang(mang),wang(mang))
       do ish = m%ishoff(k), m%ishoff(k+1)-1
          nang = m%ipsh(ish+1) - m%ipsh(ish)
          call select_lebedev(nang,xang,yang,zang,wang)
          do il = 1, nang
             ip = m%ipsh(ish) + il - 1
             m%x(:,ip) = m%xat(:,k) + m%rsh(ish) * (/xang(il),yang(il),zang(il)/)
             m%wa(ip) = m%wrsh(ish) * wang(il) / fourpi
          end do
       end do
       deallocate(xang,yang,zang,wang)
    end do
    !$omp end parallel do

    ! partition weights
    if (tmesh == mesh_type_becke) then
       call partition_becke(m)
    else
       call partition_promolecular(m,c)
    end if
    m%w = m%wa * m%wp

  end subroutine genmesh

  !> Calculate one or more scalar fields on the molecular mesh (m)
  !> using field f. prop is a one-dimensional array containing an
  !> integer label for the property (param module).  If prop is 100+n,
  !> then calculate molecular orbital number n (only for if ff is a
  !> molecular wavefunction).  If periodic, assume the mesh is for a
  !> crystal (calculate the properties by moving the points back to
  !> the main cell.
  module subroutine fillmesh(m,ff,prop,periodic)
    use fieldmod, only: field
    use tools_io, only: faterr, ferror
    use tools_math, only: bhole
    use types, only: scalar_value, realloc, field_evaluation_avail,&
       fieldeval_category_f, fieldeval_category_fder1, fieldeval_category_fder2,&
       fieldeval_category_gkin, fieldeval_category_mo, id_mo_id
    use param, only: im_volume, im_rho, im_gradrho, im_gkin, im_b, im_null
    class(mesh), intent(inout) :: m
    type(field), intent(inout) :: ff
    integer, intent(in) :: prop(:)
    logical, intent(in) :: periodic

    type(scalar_value) :: res
    integer :: i, j, n, nder
    real*8 :: fval, rhos, laps, tau, drhos2, dsigs, quads, br_b, br_alf, br_a
    character*10 :: fder
    type(field_evaluation_avail) :: request, requestmo

    if (.not.ff%isinit) &
       call ferror("fillmesh","field not initialized",faterr)

    n = size(prop)
    if (allocated(m%f)) then
       if (size(m%f,2) /= n) call realloc(m%f,m%n,n)
    else
       allocate(m%f(m%n,n))
    end if
    call request%clear()
    call requestmo%clear()
    requestmo%avail(fieldeval_category_mo) = .true.
    requestmo%moini = id_mo_id

    ! calculate the maximum derivative
    nder = -1
    do j = 1, n
       select case(prop(j))
       case(im_rho)
          nder = 0
          request%avail(fieldeval_category_f) = .true.
       case(im_gradrho)
          nder = 1
          request%avail(fieldeval_category_fder1) = .true.
       case(im_gkin)
          nder = 1
          request%avail(fieldeval_category_gkin) = .true.
       case(im_b)
          nder = 2
          request%avail(fieldeval_category_fder1) = .true.
          request%avail(fieldeval_category_fder2) = .true.
          request%avail(fieldeval_category_gkin) = .true.
       end select
    end do

    !$omp parallel do private(fval,res,fder,rhos,laps,tau,drhos2,dsigs,quads,br_b,br_alf,br_a) schedule(dynamic)
    do i = 1, m%n
       if (nder >= 0) then
          call ff%grd(m%x(:,i),request,res,periodic=periodic)
          if (.not.res%satisfied) &
             call ferror("fillmesh","field unable to provide requested properties",faterr)
       end if

       do j = 1, n
          if (prop(j) == im_volume) then
             fval = 1d0
          else if (prop(j) == im_rho) then
             fval = res%f
          else if (prop(j) == im_gradrho) then
             fval = res%gfmod
          else if (prop(j) == im_gkin) then
             fval = res%gkin
          else if (prop(j) == im_b) then
             rhos = 0.5d0 * res%f
             laps = 0.5d0 * res%del2f
             tau = res%gkin
             drhos2 = (0.5d0 * res%gfmod)
             drhos2 = drhos2 * drhos2
             dsigs = tau - 0.25d0 * drhos2 / max(rhos,1d-30)
             quads = (laps - 2d0 * dsigs) / 6d0
             call bhole(rhos,quads,1d0,br_b,br_alf,br_a)
             fval = br_b
          else if (prop(j) > 100) then
             requestmo%moend = prop(j) - 100
             call ff%grd(m%x(:,i),requestmo,res,periodic=periodic)
             if (.not.res%satisfied) &
                call ferror("fillmesh","field unable to provide individual MO values",faterr)
             fval = res%fspc
          else if (prop(j) == im_null) then
             fval = 0d0
          else
             call ferror("fillmesh","unknown quantity to calculate in fillmesh",faterr)
          end if
          m%f(i,j) = fval
       end do
    end do
    !$omp end parallel do

  end subroutine fillmesh

  !> Spherical average of the function f (given on the mesh points)
  !> around the atomic grid iat: for each radial shell of that atom,
  !> the average of f over the shell with the shell's angular weights.
  !> The radii of the shells are m%rsh(m%ishoff(iat):m%ishoff(iat+1)-1).
  module function spherical_average(m,iat,f) result(fr)
    class(mesh), intent(in) :: m
    integer, intent(in) :: iat
    real*8, intent(in) :: f(:)
    real*8, allocatable :: fr(:)

    integer :: ish, i0, i1, j

    allocate(fr(m%ishoff(iat+1)-m%ishoff(iat)))
    j = 0
    do ish = m%ishoff(iat), m%ishoff(iat+1)-1
       j = j + 1
       i0 = m%ipsh(ish)
       i1 = m%ipsh(ish+1) - 1
       fr(j) = sum(m%wa(i0:i1) * f(i0:i1)) / m%wrsh(ish)
    end do

  end function spherical_average

  !> Write information about the mesh to the standard output.
  module subroutine report(m)
    use tools_io, only: uout, string
    class(mesh), intent(inout) :: m

    integer :: i, nfb
    character(len=:), allocatable :: str

    write (uout,'("  Mesh size      ",A)') string(m%n)
    if (m%type == mesh_type_becke) then
       write (uout,'("  Mesh type      Becke")')
    elseif (m%type == mesh_type_franchini) then
       write (uout,'("  Mesh type      Franchini")')
    end if
    if (m%lvl == mesh_level_small) then
       write (uout,'("  Mesh level     small")')
    elseif (m%lvl == mesh_level_normal) then
       write (uout,'("  Mesh level     normal")')
    elseif (m%lvl == mesh_level_good) then
       write (uout,'("  Mesh level     good")')
    elseif (m%lvl == mesh_level_vgood) then
       write (uout,'("  Mesh level     very good")')
    elseif (m%lvl == mesh_level_amazing) then
       write (uout,'("  Mesh level     amazing")')
    end if
    if (m%nat > 0) then
       nfb = count(.not.m%fromtable(1:m%nat))
       if (nfb == m%nat) then
          write (uout,'("  Atomic grids   legacy scheme (",A," atoms, ",A," shells)")') &
             string(m%nat), string(size(m%rsh))
       else
          write (uout,'("  Atomic grids   tv-13.7-",A," (",A," atoms, ",A," shells)")') &
             string(level2tier(m%lvl)), string(m%nat), string(size(m%rsh))
       end if
       if (nfb > 0 .and. nfb < m%nat) then
          str = ""
          do i = 1, m%nat
             if (m%fromtable(i)) cycle
             str = str // " " // string(m%zat(i))
             if (m%zeffat(i) /= m%zat(i)) str = str // "(" // string(m%zeffat(i)) // ")"
          end do
          write (uout,'("  Legacy grids   ",A," atoms not covered by the table, Z(Zeff):",A)') &
             string(nfb), str
       end if
    end if
    write (uout,'("  Please cite: ")')
    if (m%type == mesh_type_becke) then
       write (uout,'("    A. D. Becke, J. Chem. Phys. 88, 2547 (1988). (10.1063/1.454033)")')
    end if
    write (uout,'("    V. I. Lebedev and D. N. Laikov, Dokl. Math. 59, 477 (1999).")')
    if (m%nat > 0) then
       if (nfb < m%nat) then
          write (uout,'("    Tuned atomic grids from HORTON 2.1.0: T. Verstraelen, P. Tecmer, F. Heidar-Zadeh,")')
          write (uout,'("       C. E. Gonzalez-Espinoza, M. Chan, T. D. Kim, K. Boguslawski, S. Fias,")')
          write (uout,'("       S. Vandenbrande, D. Berrocal, and P. W. Ayers, http://theochem.github.com/horton/ (2017).")')
       end if
    end if

  end subroutine report

  !xx! private procedures

  !> Name of the tuned atomic grid table file for mesh level lvl.
  function table_filename(lvl) result(file)
    use global, only: critic_home
    use tools_io, only: string
    use param, only: dirsep
    integer, intent(in) :: lvl
    character(len=:), allocatable :: file

    file = trim(critic_home) // dirsep // "meshes" // dirsep // "tv-13.7-" //&
       string(level2tier(lvl)) // ".txt"

  end function table_filename

  !> Load the tuned atomic grid table for mesh level lvl from
  !> dat/meshes/tv-13.7-N.txt. Only the first call for a given level
  !> reads the file. If the file does not exist or can not be parsed,
  !> a warning is written and the table is left unallocated. Each
  !> entry is three lines: "Z [Zeff]", "PowerRTransform rmin rmax
  !> npt", and the npt Lebedev orders; comments (#) and blank lines
  !> are skipped by getline.
  subroutine load_table(lvl)
    use tools_io, only: fopen_read, fclose, getline, lgetword, isinteger, isreal,&
       equal, ferror, warning
    integer, intent(in) :: lvl

    character(len=:), allocatable :: file, line, word
    integer :: lu, lp, i, n, nline, npt
    logical :: exist, ok
    type(table_entry), allocatable :: e(:)

    if (tab(lvl)%isinit) return
    tab(lvl)%isinit = .true.

    file = table_filename(lvl)
    inquire(file=file,exist=exist)
    if (.not.exist) then
       call ferror('load_table','Atomic grid table not found (using legacy grids): ' // file,warning)
       return
    end if

    ! count the entries (three lines each) and read them
    lu = fopen_read(file)
    nline = 0
    do while (getline(lu,line))
       nline = nline + 1
    end do
    rewind(lu)
    allocate(e(nline/3))
    ok = (mod(nline,3) == 0)
    n = 0
    do while (ok .and. n < size(e))
       n = n + 1
       ! Z [Zeff]
       ok = getline(lu,line)
       lp = 1
       ok = ok .and. isinteger(e(n)%z,line,lp)
       if (.not.ok) exit
       e(n)%zeff = e(n)%z
       ok = isinteger(e(n)%zeff,line,lp) ! optional; zeff is left as z if absent
       ok = .true.

       ! PowerRTransform rmin rmax npt
       ok = getline(lu,line)
       lp = 1
       word = lgetword(line,lp)
       ok = ok .and. equal(word,"powerrtransform")
       ok = ok .and. isreal(e(n)%rmin,line,lp)
       ok = ok .and. isreal(e(n)%rmax,line,lp)
       ok = ok .and. isinteger(npt,line,lp)
       if (.not.ok) exit
       ok = (npt >= 8) .and. (e(n)%rmin > 0d0) .and. (e(n)%rmax > e(n)%rmin)
       if (.not.ok) exit

       ! Lebedev order on each shell
       ok = getline(lu,line)
       lp = 1
       allocate(e(n)%nang(npt))
       do i = 1, npt
          ok = ok .and. isinteger(e(n)%nang(i),line,lp)
       end do
    end do
    call fclose(lu)

    if (.not.ok) then
       call ferror('load_table','Error reading atomic grid table (using legacy grids): ' // file,warning)
       return
    end if
    call move_alloc(e,tab(lvl)%e)

  end subroutine load_table

  !> Find the entry in the table for level lvl for atomic number iz
  !> and effective charge zeff: among the entries for that element,
  !> the one with the smallest table Zeff that is >= zeff. Returns
  !> zero if there is none.
  function find_table_entry(lvl,iz,zeff) result(ie)
    integer, intent(in) :: lvl, iz, zeff
    integer :: ie

    integer :: i

    ie = 0
    call load_table(lvl)
    if (.not.allocated(tab(lvl)%e)) return
    do i = 1, size(tab(lvl)%e)
       if (tab(lvl)%e(i)%z /= iz) cycle
       if (tab(lvl)%e(i)%zeff < zeff) cycle
       if (ie == 0) then
          ie = i
       elseif (tab(lvl)%e(i)%zeff < tab(lvl)%e(ie)%zeff) then
          ie = i
       end if
    end do

  end function find_table_entry

  !> Atomic grid for element iz with effective charge zeff, mesh type
  !> and level: radii r(:), radial weights wr(:) (including 4*pi*r^2)
  !> and Lebedev order on each shell nang(:). Taken from the tuned
  !> table if the element is covered (fromtable = .true.); otherwise,
  !> from the legacy scheme (constant angular order from z2nang, and
  !> the postg (Becke type) or Franchini radial grids; the postg grid
  !> is markedly more accurate for heavy all-electron atoms in
  !> molecules).
  subroutine atomic_grid_spec(iz,zeff,type,lvl,r,wr,nang,fromtable)
    use tools_math, only: good_lebedev
    integer, intent(in) :: iz, zeff, type, lvl
    real*8, allocatable, intent(inout) :: r(:), wr(:)
    integer, allocatable, intent(inout) :: nang(:)
    logical, intent(out) :: fromtable

    integer :: ie, nr, i, nang1

    if (allocated(r)) deallocate(r)
    if (allocated(wr)) deallocate(wr)
    if (allocated(nang)) deallocate(nang)

    ie = find_table_entry(lvl,iz,zeff)
    fromtable = (ie > 0)
    if (fromtable) then
       nr = size(tab(lvl)%e(ie)%nang)
       allocate(r(nr),wr(nr),nang(nr))
       call rmesh_power(tab(lvl)%e(ie)%rmin,tab(lvl)%e(ie)%rmax,nr,r,wr)
       nang = tab(lvl)%e(ie)%nang
       do i = 1, nr
          call good_lebedev(nang(i))
       end do
    else
       nr = z2nr(iz,lvl)
       allocate(r(nr),wr(nr),nang(nr))
       if (type == mesh_type_becke) then
          call rmesh_postg(nr,iz,r,wr)
       else
          call rmesh_franchini(nr,iz,r,wr)
       end if
       nang1 = z2nang(iz,lvl)
       call good_lebedev(nang1)
       nang = nang1
    end if

  end subroutine atomic_grid_spec

  !> Radial grid from the power transform used in HORTON's tuned
  !> tables: r_i = rmin * (i+1)^p for i = 0,...,n-1, with p =
  !> ln(rmax/rmin)/ln(n), integrated with the composite Simpson rule
  !> in i. The weights include the 4*pi*r^2 factor. n >= 8.
  subroutine rmesh_power(rmin,rmax,n,r,wintr)
    use param, only: fourpi
    real*8, intent(in) :: rmin, rmax
    integer, intent(in) :: n
    real*8, intent(out) :: r(n), wintr(n)

    real*8, parameter :: wend(4) = (/17d0/48d0, 59d0/48d0, 43d0/48d0, 49d0/48d0/)
    real*8 :: p, drdt, w1d
    integer :: i

    p = log(rmax/rmin) / log(real(n,8))
    do i = 1, n
       r(i) = rmin * real(i,8)**p
       drdt = p * rmin * real(i,8)**(p-1d0)
       if (i <= 4) then
          w1d = wend(i)
       elseif (i > n-4) then
          w1d = wend(n-i+1)
       else
          w1d = 1d0
       end if
       wintr(i) = fourpi * r(i)**2 * drdt * w1d
    end do

  end subroutine rmesh_power

  !> Becke partition weights (A. D. Becke, J. Chem. Phys. 88 (1988)
  !> 2547) for the mesh points, with k=3 iterations of the smoothing
  !> polynomial and the heteronuclear size adjustment of the appendix
  !> (Pyykko covalent radii from param, a_ij clamped to +-0.45 as in
  !> HORTON; the fixed table atmcov0 is used so that the RADII keyword
  !> does not change the integration weights). Only for molecules.
  subroutine partition_becke(m)
    use param, only: atmcov0
    type(mesh), intent(inout) :: m

    integer :: i, j, k, ip, it
    real*8 :: chi, u, a, mu, nu, s, vpsum
    real*8, allocatable :: rr(:,:), aij(:,:), d(:), p(:)

    ! interatomic distances and size-adjustment parameters
    allocate(rr(m%nat,m%nat),aij(m%nat,m%nat))
    rr = 0d0
    aij = 0d0
    do i = 1, m%nat
       do j = i+1, m%nat
          rr(i,j) = norm2(m%xat(:,i) - m%xat(:,j))
          rr(j,i) = rr(i,j)
          chi = atmcov0(m%zat(i)) / atmcov0(m%zat(j))
          u = (chi - 1d0) / (chi + 1d0)
          a = u / (u*u - 1d0)
          a = max(min(a,0.45d0),-0.45d0)
          aij(i,j) = a
          aij(j,i) = -a
       end do
    end do

    ! cell functions: s(j,i) = -s(i,j), so each pair is done once
    !$omp parallel do private(ip,d,p,i,j,mu,nu,s,it,vpsum) schedule(dynamic)
    do k = 1, m%nat
       allocate(d(m%nat),p(m%nat))
       do ip = m%ipsh(m%ishoff(k)), m%ipsh(m%ishoff(k+1))-1
          do i = 1, m%nat
             d(i) = norm2(m%x(:,ip) - m%xat(:,i))
          end do
          p = 1d0
          do i = 1, m%nat
             do j = i+1, m%nat
                mu = (d(i) - d(j)) / rr(i,j)
                nu = mu + aij(i,j) * (1d0 - mu*mu)
                s = nu
                do it = 1, 3
                   s = 1.5d0 * s - 0.5d0 * s**3
                end do
                p(i) = p(i) * 0.5d0 * (1d0 - s)
                p(j) = p(j) * 0.5d0 * (1d0 + s)
             end do
          end do
          vpsum = sum(p)
          if (vpsum > 0d0) then
             m%wp(ip) = p(k) / vpsum
          else
             m%wp(ip) = 0d0
          end if
       end do
       deallocate(d,p)
    end do
    !$omp end parallel do

  end subroutine partition_becke

  !> Promolecular partition weights for the mesh points, using
  !> exp(-2r)/r^3 as the atomic weight function (scaled by 0.3 for
  !> hydrogen) and a neighbor list, so the cost does not grow with
  !> the number of atoms and the mesh works for crystals. The list
  !> reaches the largest grid radius so that every atom closer to the
  !> point than its owner is included.
  subroutine partition_promolecular(m,c)
    use crystalmod, only: crystal
    use param, only: maxzat, icrd_cart
    type(mesh), intent(inout) :: m
    type(crystal), intent(inout) :: c

    real*8, parameter :: rthres = 12d0 ! contribution to weight: 2e-14

    integer :: k, ish, ip, j, iz, nat
    real*8 :: r, vp0, vpsum, rmax
    integer, allocatable :: eid(:)
    real*8, allocatable :: dist(:)

    rmax = max(maxval(m%rsh),rthres)

    !$omp parallel do private(ish,ip,r,vp0,vpsum,nat,eid,dist,j,iz) schedule(dynamic)
    do k = 1, m%nat
       do ish = m%ishoff(k), m%ishoff(k+1)-1
          r = m%rsh(ish)
          vp0 = wfun(m%zat(k),r)
          do ip = m%ipsh(ish), m%ipsh(ish+1)-1
             ! all atoms within rmax of the mesh point
             call c%list_near_atoms(m%x(:,ip),icrd_cart,.false.,nat,eid=eid,dist=dist,up2d=rmax)
             vpsum = 0d0
             do j = 1, nat
                iz = c%spc(c%atcel(eid(j))%is)%z
                if (iz < 1 .or. iz > maxzat) cycle
                vpsum = vpsum + wfun(iz,dist(j))
             end do
             vpsum = max(vp0,vpsum)
             if (vpsum > 0d0) then
                m%wp(ip) = vp0 / vpsum
             else
                m%wp(ip) = 0d0
             end if
          end do
       end do
    end do
    !$omp end parallel do

  contains
    function wfun(iz,r)
      integer, intent(in) :: iz
      real*8, intent(in) :: r
      real*8 :: wfun

      if (iz == 1) then
         wfun = 0.3d0 * exp(-2d0 * r) / max(r,1d-10)**3
      else
         wfun = exp(-2d0 * r) / max(r,1d-10)**3
      end if

    end function wfun
  end subroutine partition_promolecular

  !> Radial mesh for integration of exponential-like functions.
  !> From postg.
  subroutine rmesh_postg(n,iz,r,wintr)
    ! Radial mesh and integration weights,
    ! derivatives of variable r with respect to variable q.
    ! The q-mesh is uniform on the interval (0,+1). Transformation is
    ! r = rmid * q / (1-q)
    use param, only: third, fourpi
    integer, intent(in) :: n
    integer, intent(in) :: iz
    real*8, intent(out) :: r(n), wintr(n)

    real*8 :: h, q, rmid
    integer :: i

    rmid = 1d0/real(iz,8)**third
    h = 1d0/real(n+1,8)
    do i=1,n
       q=h*i
       r(i) = rmid * q / (1.d0-q)
       wintr(i) = fourpi * h * r(i)**2 * rmid / (1.d0-q)**2
    enddo

  end subroutine rmesh_postg

  !> Radial mesh for integration of exponential-like functions.
  !> After Franchini et al., J. Comput. Chem. 34 (2013) 1819.
  !> The transformation is:
  !>   r = zeta/ln(2) * (1+q) * ln(2/(1-q))
  !> where zeta is a function of the atomic number.
  subroutine rmesh_franchini(n,iz,r,wintr)
    use tools_math, only: gauleg
    use param, only: log2, fourpi
    integer, intent(in) :: n
    integer, intent(in) :: iz
    real*8, intent(out) :: r(n), wintr(n)

    !> From the SI in Franchini's.
    real*8, parameter :: zeta(103) = (/&
       0.8, 0.9, 1.8, 1.4, 1.3, 1.1, 0.9, 0.9, 0.9, 0.9,& !  1-10
       1.4, 1.3, 1.3, 1.2, 1.1, 1.0, 1.0, 1.0, 1.5, 1.4,& ! 11-20
       1.3, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.1, 1.1, 1.1,& ! 21-30
       1.1, 1.0, 0.9, 0.9, 0.9, 0.9, 1.4, 1.4, 1.1, 1.3,& ! 31-40
       1.0, 1.2, 0.9, 0.9, 0.9, 1.0, 0.9, 1.0, 1.0, 1.3,& ! 41-50
       1.2, 1.2, 0.9, 1.0, 1.7, 1.5, 1.5, 1.3, 1.3, 1.4,& ! 51-60
       1.8, 1.4, 1.2, 1.3, 1.3, 1.4, 1.1, 1.1, 1.2, 1.6,& ! 61-70
       1.4, 1.3, 1.2, 1.0, 1.0, 0.9, 1.3, 1.2, 1.2, 1.0,& ! 71-80
       1.2, 1.2, 1.1, 1.2, 1.1, 2.1, 2.2, 1.8, 1.7, 1.3,& ! 81-90
       1.4, 1.2, 1.2, 1.3, 1.4, 1.4, 1.7, 1.9, 1.9, 2.0,& ! 91-100
       2.0, 1.6, 2.0/) ! 101-103

    real*8 :: h, q
    integer :: i

    call gauleg(-1d0,1d0,r,wintr,n)
    h = 1d0/real(n+1,8)
    do i = 1, n
       q = r(i)
       r(i) = zeta(iz)/log2 * (1+q) * log(2/(1-q))
       wintr(i) = fourpi * r(i)**2 * wintr(i) * zeta(iz)/log2 * (log(2/(1-q)) + (1+q)/(1-q))
    enddo

  end subroutine rmesh_franchini

  !> Atomic number to radial point number.
  !> After Franchini et al., J. Comput. Chem. 34 (2013) 1819.
  !> lvl = 1 (small), 2 (normal), 3 (good), 4(very good), 5 (excellent)
  !> Some changes to the choice of ni0
  function z2nr(z,lvl) result(nr)
    use tools_io, only: faterr, ferror
    integer, intent(in) :: z
    integer, intent(in) :: lvl
    integer :: nr

    nr = 15
    if (z > 2) nr = 20
    if (z > 10) nr = 25
    if (z > 18) nr = 35
    if (z > 36) nr = 60
    if (z > 54) nr = 85
    if (z > 86) nr = 110

    if (lvl == mesh_level_small) then
       nr = ceiling(nr * 2.37)
    elseif (lvl == mesh_level_normal) then
       nr = ceiling(nr * 3.08)
    elseif (lvl == mesh_level_good) then
       nr = ceiling(nr * 3.42)
    elseif (lvl == mesh_level_vgood) then
       nr = ceiling(nr * 4.27)
    elseif (lvl == mesh_level_amazing) then
       nr = ceiling(nr * 6.72)
    else
       call ferror("z2nr","unknown accuracy level",faterr)
    end if

  endfunction z2nr

  !> Select the number of points for the angular (Lebedev) quadrature.
  !> After Franchini et al., J. Comput. Chem. 34 (2013) 1819.
  !> lvl = 1 (small), 2 (normal), 3 (good), 4(very good), 5 (excellent)
  function z2nang(z,lvl) result(nang)
    use tools_io, only: ferror, faterr
    integer, intent(in) :: z
    integer, intent(in) :: lvl
    integer :: nang

    if (lvl == mesh_level_small) then
       nang = 110
    elseif (lvl == mesh_level_normal) then
       nang = 194
    elseif (lvl == mesh_level_good) then
       nang = 302
    elseif (lvl == mesh_level_vgood) then
       nang = 590
    elseif (lvl == mesh_level_amazing) then
       nang = 770
    else
       call ferror("z2nang","unknown accuracy level",faterr)
    end if

  endfunction z2nang

end submodule proc
