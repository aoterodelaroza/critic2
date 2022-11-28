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

! Structure class and routines for basic crystallography computations
submodule (crystalmod) proc
  implicit none

  !xx! private procedures
  ! subroutine typeop(rot,type,vec,order)
  ! function equiv_tetrah(c,x0,t1,t2,leqv,lrotm,eps)
  ! function perm3(p,r,t) result(res)

  ! symmetry operation symbols
  integer, parameter :: ident=0 !< identifier for sym. operations
  integer, parameter :: inv=1 !< identifier for sym. operations
  integer, parameter :: c2=2 !< identifier for sym. operations
  integer, parameter :: c3=3 !< identifier for sym. operations
  integer, parameter :: c4=4 !< identifier for sym. operations
  integer, parameter :: c6=5 !< identifier for sym. operations
  integer, parameter :: s3=6 !< identifier for sym. operations
  integer, parameter :: s4=7 !< identifier for sym. operations
  integer, parameter :: s6=8 !< identifier for sym. operations
  integer, parameter :: sigma=9 !< identifier for sym. operations

  ! array initialization values
  integer, parameter :: mspc0 = 4
  integer, parameter :: mneq0 = 4
  integer, parameter :: mcel0 = 10

contains

  !xx! crystal class methods
  !> Initialize the struct arrays
  module subroutine struct_init(c)
    class(crystal), intent(inout) :: c

    integer :: i

    ! deallocate all the arrays and reinitialize
    call c%end()

    ! allocate space for atoms
    if (.not.allocated(c%spc)) allocate(c%spc(mspc0))
    if (.not.allocated(c%at)) allocate(c%at(mneq0))
    if (.not.allocated(c%atcel)) allocate(c%atcel(mcel0))
    c%nspc = 0
    c%nneq = 0
    c%ncel = 0

    ! nullify metrics
    c%aa = 0d0
    c%bb = 0d0
    c%ar = 0d0
    c%omega = 0d0
    c%gtensor = 0d0
    c%grtensor = 0d0
    c%m_x2c = 0d0
    c%m_c2x = 0d0
    c%n2_x2c = 0d0
    c%n2_c2x = 0d0

    ! no symmetry
    call c%clearsym()

    ! no molecule
    c%ismolecule = .false.
    c%molx0 = 0d0
    c%molborder = 0d0

    ! no ws
    c%ws_nv = 0
    c%ws_nf = 0
    c%ws_mnfv = 0
    if (allocated(c%ws_ineighx)) deallocate(c%ws_ineighx)
    if (allocated(c%ws_ineighxr)) deallocate(c%ws_ineighxr)
    if (allocated(c%ws_ineighc)) deallocate(c%ws_ineighc)
    if (allocated(c%ws_nside)) deallocate(c%ws_nside)
    if (allocated(c%ws_iside)) deallocate(c%ws_iside)
    if (allocated(c%ws_x)) deallocate(c%ws_x)
    c%isortho = .false.
    c%isortho_del = .false.

    ! initialize species
    do i = 1, mspc0
       c%spc(i)%name = ""
       c%spc(i)%z = 0
       c%spc(i)%qat = 0d0
    end do

    ! initialize atoms
    do i = 1, mneq0
       c%at(i)%rnn2 = 0d0
    end do

    ! no molecular fragments
    c%nmol = 0
    if (allocated(c%nstar)) deallocate(c%nstar)
    if (allocated(c%mol)) deallocate(c%mol)
    c%nlvac = 0
    c%lvac = 0
    c%lcon = 0
    if (allocated(c%idatcelmol)) deallocate(c%idatcelmol)

    ! no 3d molecular crystals
    c%ismol3d = .false.

    ! core charges
    c%zpsp = -1

    ! the crystal is not initialized until struct_new is run
    c%file = ""
    c%isinit = .false.
    c%havesym = 0
    c%isewald = .false.

  end subroutine struct_init

  !> Terminate allocated arrays
  module subroutine struct_end(c)
    class(crystal), intent(inout) :: c

    c%isinit = .false.
    if (allocated(c%spc)) deallocate(c%spc)
    if (allocated(c%at)) deallocate(c%at)
    if (allocated(c%atcel)) deallocate(c%atcel)
    if (allocated(c%cen)) deallocate(c%cen)
    if (allocated(c%ws_iside)) deallocate(c%ws_iside)
    if (allocated(c%ws_x)) deallocate(c%ws_x)
    if (allocated(c%nstar)) deallocate(c%nstar)
    if (allocated(c%mol)) deallocate(c%mol)
    if (allocated(c%idatcelmol)) deallocate(c%idatcelmol)
    call c%env%end()
    c%isinit = .false.
    c%havesym = 0
    c%isewald = .false.
    c%file = ""
    c%nspc = 0
    c%nneq = 0
    c%ncel = 0
    c%neqv = 0
    c%ncv = 0
    c%ismolecule = .false.
    c%molx0 = 0d0
    c%molborder = 0d0
    c%ws_nv = 0
    c%ws_nf = 0
    c%ws_mnfv = 0
    c%nmol = 0
    c%nlvac = 0
    c%ismol3d = .false.

  end subroutine struct_end

  !> Create a new, complete crystal/molecule from a crystal seed. If
  !> failed and crashfail is true, crash the program. Otherwise,
  !> return a the error status through c%isinit. If noenv is present
  !> and true, do not load the atomic grids or the environment.
  module subroutine struct_new(c,seed,crashfail,noenv,ti)
    use crystalseedmod, only: crystalseed
    use grid1mod, only: grid1_register_ae
    use global, only: crsmall, atomeps
    use tools_math, only: m_x2c_from_cellpar, m_c2x_from_cellpar, matinv, &
       det3, mnorm2
    use tools_io, only: ferror, faterr, zatguess, string
    use tools, only: wscell
    use types, only: realloc
    use param, only: pi, eyet, icrd_cart
    class(crystal), intent(inout) :: c
    type(crystalseed), intent(in) :: seed
    logical, intent(in) :: crashfail
    logical, intent(in), optional :: noenv
    type(thread_info), intent(in), optional :: ti

    real*8 :: g(3,3), xmax(3), xmin(3), xcm(3), dist, border, xx(3)
    logical :: good, clearsym, doenv
    integer :: i, j, k, l, iat, newmult
    real*8, allocatable :: atpos(:,:), area(:)
    integer, allocatable :: irotm(:), icenv(:)
    logical, allocatable :: useatom(:)
    character(len=:), allocatable :: errmsg
    logical :: haveatoms

    if (.not.seed%isused) then
       if (crashfail) then
          call ferror("struct_new","uninitialized seed",faterr)
       else
          return
       end if
    end if

    ! initialize the structure
    call c%init()
    c%ismolecule = seed%ismolecule

    !! first the cell, then the atoms !!

    ! if this is a molecule, calculate the center and encompassing cell
    border = max(seed%border,1d-6)
    if (seed%ismolecule .and. seed%nat > 0) then
       xmax = -1d40
       xmin =  1d40
       do i = 1, seed%nat
          do j = 1, 3
             xmax(j) = max(seed%x(j,i)+border,xmax(j))
             xmin(j) = min(seed%x(j,i)-border,xmin(j))
          end do
       end do
       if (seed%cubic) then
          xmin = minval(xmin)
          xmax = maxval(xmax)
       end if
       xcm = xmin + 0.5d0 * (xmax-xmin)
    else
       xmax = 1d0
       xcm = 0.5d0
       xmin = 0d0
    end if

    ! basic cell and centering information
    if (seed%useabr == 0) then
       if (.not.seed%ismolecule) then
          if (crashfail) then
             call ferror("struct_new","cell data unavailable",faterr)
          else
             return
          end if
       end if
       ! this is a molecule, for which no cell has been given
       c%aa = xmax - xmin
       c%bb = 90d0

       c%m_x2c = m_x2c_from_cellpar(c%aa,c%bb)
       c%m_c2x = c%m_x2c
       call matinv(c%m_c2x,3)
       g = matmul(transpose(c%m_x2c),c%m_x2c)
    elseif (seed%useabr == 1) then
       ! use aa and bb
       c%aa = seed%aa
       c%bb = seed%bb
       c%m_x2c = m_x2c_from_cellpar(c%aa,c%bb)
       c%m_c2x = c%m_x2c
       call matinv(c%m_c2x,3)
       g = matmul(transpose(c%m_x2c),c%m_x2c)
    elseif (seed%useabr == 2) then
       ! use m_x2c
       c%m_x2c = seed%m_x2c
       c%m_c2x = c%m_x2c
       call matinv(c%m_c2x,3)
       g = matmul(transpose(c%m_x2c),c%m_x2c)
       do i = 1, 3
          c%aa(i) = sqrt(g(i,i))
       end do
       c%bb(1) = acos(g(2,3) / c%aa(2) / c%aa(3)) * 180d0 / pi
       c%bb(2) = acos(g(1,3) / c%aa(1) / c%aa(3)) * 180d0 / pi
       c%bb(3) = acos(g(1,2) / c%aa(1) / c%aa(2)) * 180d0 / pi
    else
       if (crashfail) then
          call ferror("struct_new","unknown useabr",faterr)
       else
          return
       end if
    end if

    ! rest of the cell metrics
    c%gtensor = g
    c%omega = sqrt(max(det3(g),0d0))
    c%grtensor = g
    call matinv(c%grtensor,3)
    do i = 1, 3
       c%ar(i) = sqrt(c%grtensor(i,i))
    end do
    c%n2_x2c = mnorm2(c%m_x2c)
    c%n2_c2x = mnorm2(c%m_c2x)

    ! calculate the wigner-seitz cell
    call wscell(c%m_x2c,.not.c%ismolecule,c%ws_nf,c%ws_nv,c%ws_mnfv,c%ws_iside,c%ws_nside,&
       c%ws_x,c%ws_ineighc,c%ws_ineighx,area,c%isortho,c%m_xr2x,c%m_x2xr,c%m_xr2c,c%m_c2xr,&
       c%n2_xr2x,c%n2_x2xr,c%n2_xr2c,c%n2_c2xr,c%ws_ineighxr,c%isortho_del)

    !! now the atoms, we can use shortest() and related functions past this !!

    ! check for repeats in the seed (for cif & shelx formats)
    allocate(useatom(seed%nat))
    useatom = .true.
    if (seed%havesym > 0 .and. .not.seed%ismolecule .and. seed%nat > 0 .and. seed%checkrepeats) then
       do i = 1, seed%nat
          if (.not.useatom(i)) cycle
          do j = 1, seed%neqv
             do k = 1, seed%ncv
                xx = matmul(seed%rotm(1:3,1:3,j),seed%x(:,i)) + seed%rotm(:,4,j) + seed%cen(:,k)
                do l = i+1, seed%nat
                   if (.not.useatom(l)) cycle
                   if (c%eql_distance(seed%x(:,l),xx) < atomeps) useatom(l) = .false.
                end do
             enddo
          enddo
       end do
    end if

    ! copy the atomic information
    c%nspc = seed%nspc
    c%nneq = count(useatom)
    if (c%nneq > 0) then
       if (allocated(c%at)) deallocate(c%at)
       allocate(c%at(c%nneq))
       if (allocated(c%spc)) deallocate(c%spc)
       allocate(c%spc(c%nspc))
       c%spc = seed%spc(1:seed%nspc)
       i = 0
       do j = 1, seed%nat
          if (useatom(j)) then
             i = i + 1
             c%at(i)%x = seed%x(:,j)
             c%at(i)%is = seed%is(j)
          end if
       end do
    end if
    deallocate(useatom)

    ! check if we have any atoms
    haveatoms = .false.
    if (c%nneq > 0) then
       do i = 1, c%nspc
          if (c%spc(i)%z > 0) then
             haveatoms = .true.
             exit
          end if
       end do
    end if

    ! transform the atomic coordinates in the case of a molecule, and fill
    ! the remaining molecular fields
    if (seed%ismolecule) then
       if (seed%useabr == 0) then
          ! a cell has not been given

          ! center in the cell and convert to crystallographic coordinates
          do i = 1, c%nneq
             c%at(i)%x = (c%at(i)%x-xcm+0.5d0*c%aa) / c%aa
          end do

          ! Keep the (1/2,1/2,1/2) translation applied
          c%molx0 = -(/0.5d0, 0.5d0, 0.5d0/) * c%aa + xcm

          ! Set up the molecular cell. c%molborder is in fractional coordinates
          ! and gives the position of the molecular cell in each axis. By default,
          ! choose the molecular cell as the minimal encompassing cell for the molecule
          ! plus 80% of the border or 2 bohr, whichever is larger. The molecular cell
          ! can not exceed the actual unit cell
          c%molborder = max(border - max(2d0,0.8d0 * border),0d0) / (xmax - xmin)
       else
          if (any(abs(c%bb - 90d0) > 1d-3)) then
             if (crashfail) then
                call ferror("struct_new","MOLECULE does not allow non-orthogonal cells",faterr)
             else
                return
             end if
          end if
          ! a cell has been given, save the origin
          if (seed%havex0) then
             c%molx0 = seed%molx0
          else
             c%molx0 = -(/0.5d0, 0.5d0, 0.5d0/) * c%aa
          endif

          ! calculate the molecular cell
          if (c%nneq > 0) then
             xmin = 1d40
             do i = 1, c%nneq
                do j = 1, 3
                   xmin(j) = min(c%at(i)%x(j),xmin(j))
                   xmin(j) = min(1d0-max(c%at(i)%x(j),1d0-xmin(j)),xmin(j))
                end do
             end do
          end if
          c%molborder = max(xmin - max(0.8d0 * xmin,2d0/c%aa),0d0)
       end if
    end if

    ! move the crystallographic coordinates to the main cell, calculate the
    ! Cartesian coordinates. Set the default wyckoff letter.
    do i = 1, c%nneq
       c%at(i)%x = c%at(i)%x - floor(c%at(i)%x)
       c%at(i)%r = c%x2c(c%at(i)%x)
       c%at(i)%wyc = "?"
    end do

    ! copy the symmetry information, if available
    if (seed%havesym > 0 .and..not.seed%ismolecule) then
       c%havesym = 1
       c%neqv = seed%neqv
       c%ncv = seed%ncv
       if (allocated(c%cen)) deallocate(c%cen)
       allocate(c%cen(3,seed%ncv))
       c%cen = seed%cen(:,1:seed%ncv)
       c%rotm = 0d0
       c%rotm(:,:,1:seed%neqv) = seed%rotm(:,:,1:seed%neqv)

       ! permute the symmetry operations to make the identity the first
       if (c%neqv > 1) then
          if (.not.all(abs(eyet - c%rotm(:,:,1)) < 1d-12)) then
             good = .false.
             do i = 1, c%neqv
                if (all(abs(eyet - c%rotm(:,:,i)) < 1d-12)) then
                   c%rotm(:,:,i) = c%rotm(:,:,1)
                   c%rotm(:,:,1) = eyet
                   good = .true.
                   exit
                end if
             end do
             if (.not.good) then
                if (crashfail) then
                   call ferror('struct_new','identity operation not found',faterr)
                else
                   return
                end if
             end if
          end if
       end if

       ! generate the complete atom list
       if (c%nneq > 0) then
          if (allocated(c%atcel)) deallocate(c%atcel)
          allocate(c%atcel(c%nneq*c%neqv*c%ncv))
          iat = 0
          do i = 1, c%nneq
             call c%symeqv(c%at(i)%x,c%at(i)%mult,atpos,irotm,icenv,atomeps)

             newmult = 0
             jloop: do j = 1, c%at(i)%mult
                iat = iat + 1
                newmult = newmult + 1
                c%atcel(iat)%x = atpos(:,j)
                c%atcel(iat)%r = c%x2c(atpos(:,j))
                c%atcel(iat)%idx = i
                c%atcel(iat)%cidx = iat
                c%atcel(iat)%ir = irotm(j)
                c%atcel(iat)%ic = icenv(j)
                c%atcel(iat)%lvec = nint(atpos(:,j) - &
                   (matmul(c%rotm(1:3,1:3,irotm(j)),c%at(i)%x) + &
                   c%rotm(1:3,4,irotm(j)) + c%cen(:,icenv(j))))
                c%atcel(iat)%is = c%at(i)%is
             end do jloop

             if (newmult /= c%at(i)%mult) then
                if (crashfail) then
                   call ferror('struct_new','inconsistent multiplicity for atom ' // string(i),faterr)
                else
                   return
                end if
             end if
          end do
          c%ncel = iat
          if (allocated(atpos)) deallocate(atpos)
          if (allocated(irotm)) deallocate(irotm)
          if (allocated(icenv)) deallocate(icenv)
          call realloc(c%atcel,c%ncel)
       else
          c%ncel = 0
       end if
    else
       call c%clearsym()
    end if

    ! symmetry from spglib
    clearsym = .true.
    if (.not.seed%ismolecule .and. seed%havesym == 0) then
       if ((seed%findsym == 1 .or. seed%findsym == -1 .and. seed%nat <= crsmall) .and. haveatoms) then
          ! symmetry was not available, and I want it
          ! this operation fills the symmetry info, at(i)%mult, and ncel/atcel
          call c%calcsym(.true.,errmsg,ti=ti)
          if (len_trim(errmsg) > 0) then
             call ferror("struct_new","spglib: "//errmsg,faterr)
          else
             clearsym = .false.
          end if
       end if

    else if (.not.seed%ismolecule .and. c%havesym > 0) then
       ! symmetry was already available, but I still want the space group details
       call c%spglib_wrap(c%spg,.false.,errmsg,ti=ti)
       if (len_trim(errmsg) > 0) then
          call ferror("struct_new","spglib: "//errmsg,faterr)
       else
          c%spgavail = .true.
          call c%spgtowyc(c%spg)
          clearsym = .false.
       end if
    end if

    if (clearsym) then
       ! symmetry was not available or there was an error, and I do not
       ! want symmetry - make a copy of at() to atcel() and set P1
       call c%clearsym(neq2cel=.true.)
    end if

    doenv = .true.
    if (present(noenv)) then
       if (noenv) doenv = .false.
    end if

    if (doenv) then
       ! load the atomic density grids
       do i = 1, c%nspc
          call grid1_register_ae(c%spc(i)%z)
       end do

       if (haveatoms) then
          ! Build the atomic environments
          call c%env%build(c%ismolecule,c%nspc,c%spc(1:c%nspc),c%ncel,c%atcel(1:c%ncel),c%m_x2c)

          ! Find the atomic connectivity and the molecular fragments
          call c%env%find_asterisms_covalent(c%nstar)
          call c%fill_molecular_fragments()
          call c%calculate_molecular_equivalence()

          ! Write the half nearest-neighbor distance
          do i = 1, c%nneq
             if (.not.c%ismolecule .or. c%ncel > 1) then
                call c%env%nearest_atom(c%at(i)%r,icrd_cart,iat,dist,nozero=.true.)
                c%at(i)%rnn2 = 0.5d0 * dist
             else
                c%at(i)%rnn2 = 0d0
             end if
          end do
       end if
    end if

    ! the initialization is done - this crystal is ready to use
    c%file = seed%file
    c%isinit = .true.

  end subroutine struct_new

  !> Identify a species in the crystal from a string. Step 1: check if
  !> str is equal to the name of any species (case sensitive). 2: same
  !> but case insensitive.  3: transform to Z and check if any species
  !> has that Z (if several, return the first one). 4: return 0 (not
  !> found).
  module function identify_spc(c,str) result(res)
    use crystalseedmod, only: crystalseed
    use tools_io, only: equal, lower, zatguess
    class(crystal), intent(inout) :: c
    character*(*), intent(in) :: str
    integer :: res

    integer :: i, j, iz
    character*(len_trim(str)) :: stri

    stri = trim(str)

    ! steps 1 and 2
    do j = 1, 2
       do i = 1, c%nspc
          if (equal(c%spc(i)%name,stri)) then
             res = i
             return
          end if
       end do
       if (j == 2) exit
       stri = lower(stri)
    end do

    ! step 3
    iz = zatguess(stri)
    if (iz > 0) then
       do i = 1, c%nspc
          if (c%spc(i)%z == iz) then
             res = i
             return
          end if
       end do
    end if

    ! not found: return zero
    res = 0

  end function identify_spc

  !> Convert input cryst. -> cartesian. This routine is thread-safe.
  pure module function x2c(c,xx) result(res)
    class(crystal), intent(in) :: c
    real*8, intent(in) :: xx(3)
    real*8 :: res(3)

    res = matmul(c%m_x2c,xx)

  end function x2c

  !> Convert input cartesian -> cryst. This routine is thread-safe.
  pure module function c2x(c,xx) result(res)
    class(crystal), intent(in) :: c
    real*8, intent(in)  :: xx(3)
    real*8 :: res(3)

    res = matmul(c%m_c2x,xx)

  end function c2x

  !> Convert reduced cryst. -> cartesian. This routine is thread-safe.
  pure module function xr2c(c,xx) result(res)
    class(crystal), intent(in) :: c
    real*8, intent(in)  :: xx(3)
    real*8 :: res(3)

    res = matmul(c%m_xr2c,xx)

  end function xr2c

  !> Convert cartesian -> reduced cryst. This routine is thread-safe.
  pure module function c2xr(c,xx) result(res)
    class(crystal), intent(in) :: c
    real*8, intent(in)  :: xx(3)
    real*8 :: res(3)

    res = matmul(c%m_c2xr,xx)

  end function c2xr

  !> Convert reduced cryst. -> input cryst. This routine is thread-safe.
  pure module function xr2x(c,xx) result(res)
    class(crystal), intent(in) :: c
    real*8, intent(in)  :: xx(3)
    real*8 :: res(3)

    res = matmul(c%m_xr2x,xx)

  end function xr2x

  !> Convert input cryst. -> reduced cryst. This routine is thread-safe.
  pure module function x2xr(c,xx) result(res)
    class(crystal), intent(in) :: c
    real*8, intent(in)  :: xx(3)
    real*8 :: res(3)

    res = matmul(c%m_x2xr,xx)

  end function x2xr

  !> Compute the distance between points in crystallographic.  This
  !> routine is thread-safe.
  pure module function distance(c,x1,x2)
    class(crystal), intent(in) :: c !< Input crystal
    real*8, intent(in) :: x1(3) !< First point in cryst. coordinates
    real*8, intent(in) :: x2(3) !< Second point in cryst. coordinates
    real*8 :: distance

    real*8 :: xd(3)

    xd = c%x2c(x1 - x2)
    distance = norm2(xd)

  end function distance

  !> Compute the shortest distance between a point x1 and all
  !> lattice translations of another point x2. Input points in cryst.
  !> coordinates. This routine is thread-safe.
  pure module function eql_distance(c,x1,x2)
    class(crystal), intent(in) :: c !< Input crystal
    real*8, intent(in) :: x1(3) !< First point in cryst. coordinates
    real*8, intent(in) :: x2(3) !< Second point in cryst. coordinates
    real*8 :: eql_distance

    real*8 :: xd(3), dist

    xd = x1 - x2
    call c%shortest(xd,dist)
    eql_distance = dist

  end function eql_distance

  !> Given a point in crystallographic coordinates (x), find the
  !> lattice-translated copy of x with the shortest length. Returns
  !> the shortest-length vector in Cartesian coordinates and the
  !> distance. This routine is thread-safe.
  pure module subroutine shortest(c,x,dist)
    class(crystal), intent(in) :: c
    real*8, intent(inout) :: x(3)
    real*8, intent(out) :: dist

    integer :: i
    real*8 :: xtry(3), dvws, x0(3)

    if (c%isortho) then
       x = x - nint(x)
       x = matmul(c%m_x2c,x)
       dist = norm2(x)
    else
       x = matmul(c%m_x2xr,x)
       x = x - nint(x)
       x = matmul(c%m_xr2c,x)
       dist = norm2(x)

       if (.not.c%isortho_del) then
          x0 = x
          do i = 1, c%ws_nf
             xtry = x0 + c%ws_ineighc(:,i)
             dvws = norm2(xtry)
             if (dvws < dist) then
                x = xtry
                dist = dvws
             endif
          end do
       end if
    endif

  end subroutine shortest

  !> Determine if two points x0 and x1 (cryst.) are at a distance less
  !> than eps. Logical version of c%distance(). If d2 is present and
  !> are_close is .true., return the square of the distance in that
  !> argument.  This routine is thread-safe.
  module function are_close(c,x0,x1,eps,dd)
    use param, only: ctsq3
    class(crystal), intent(in) :: c
    real*8, intent(in) :: x0(3), x1(3)
    real*8, intent(in) :: eps
    real*8, intent(out), optional :: dd
    logical :: are_close

    real*8 :: x(3), dbound, dist

    are_close = .false.
    x = x0 - x1
    dbound = minval(abs(x)) * ctsq3 / c%n2_c2x
    if (dbound > eps) return
    x = matmul(c%m_x2c,x)
    if (any(abs(x) > eps)) return
    dist = norm2(x)
    are_close = (dist < eps)
    if (present(dd) .and. are_close) dd = dist

  end function are_close

  !> Determine if a points x0 is at a distance less than eps from x1
  !> or any of its lattice translations. x0 and x1 are in cryst.
  !> coords. Logical version of c%ldistance(). If d2 is present and
  !> are_close is .true., return the square of the distance in that
  !> argument. This routine is thread-safe.
  module function are_lclose(c,x0,x1,eps,dd)
    class(crystal), intent(in) :: c
    real*8, intent(in) :: x0(3), x1(3)
    real*8, intent(in) :: eps
    real*8, intent(out), optional :: dd
    logical :: are_lclose

    real*8 :: x(3), dist

    are_lclose = .false.
    x = x0 - x1
    call c%shortest(x,dist)
    are_lclose = (dist < eps)
    if (present(dd) .and. are_lclose) dd = dist

  end function are_lclose

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
  !> list. If is0, consider only atoms of species is0. If nozero,
  !> disregard zero-distance atoms. This routine is a wrapper for the
  !> environment's nearest_atom. Thread-safe.
  module subroutine nearest_atom(c,xp,icrd,nid,dist,distmax,lvec,cidx0,idx0,is0,nozero)
    class(crystal), intent(in) :: c
    real*8, intent(in) :: xp(3)
    integer, intent(in) :: icrd
    integer, intent(out) :: nid
    real*8, intent(out) :: dist
    real*8, intent(in), optional :: distmax
    integer, intent(out), optional :: lvec(3)
    integer, intent(in), optional :: cidx0
    integer, intent(in), optional :: idx0
    integer, intent(in), optional :: is0
    logical, intent(in), optional :: nozero

    call c%env%nearest_atom(xp,icrd,nid,dist,distmax,lvec,cidx0,idx0,is0,nozero)

  end subroutine nearest_atom

  !> Given point x0 (with icrd input coordinates), translate to the
  !> main cell if the environment is from a crystal. Then if x0
  !> corresponds to an atomic position (to within distmax or atomeps
  !> if distmax is not given), return the complete-list ID of the
  !> atom. Otherwise, return 0. Optionally, return the lattice vector
  !> translation (lvec) and the distance (dist) to the closest atom.
  !> This routine is a wrapper for the environment's identify_atom
  !> function. Thread-safe.
  module function identify_atom(c,x0,icrd,lvec,dist,distmax)
    class(crystal), intent(in) :: c
    real*8, intent(in) :: x0(3)
    integer, intent(in) :: icrd
    integer, intent(out), optional :: lvec(3)
    real*8, intent(out), optional :: dist
    real*8, intent(in), optional :: distmax
    integer :: identify_atom

    identify_atom = c%env%identify_atom(x0,icrd,lvec,dist,distmax)

  endfunction identify_atom

  !> Identify a fragment in the unit cell. Input: cartesian coords. Output:
  !> A fragment object. This routine is thread-safe.
  module function identify_fragment(c,nat,x0) result(fr)
    use types, only: realloc
    use param, only: icrd_cart
    class(crystal), intent(in) :: c
    integer, intent(in) :: nat
    real*8, intent(in) :: x0(3,nat)
    type(fragment) :: fr

    integer :: id, i, n

    fr%nat = nat
    allocate(fr%at(nat))
    fr%nspc = c%nspc
    allocate(fr%spc(c%nspc))
    fr%spc = c%spc

    n = 0
    do i = 1, nat
       id = c%identify_atom(x0(:,i),icrd_cart)
       if (id > 0) then
          n = n + 1
          fr%at(n)%r = x0(:,i)
          fr%at(n)%x = c%c2x(x0(:,i))
          fr%at(n)%cidx = id
          fr%at(n)%idx = c%atcel(id)%idx
          fr%at(n)%lvec = nint(fr%at(n)%x - c%atcel(id)%x)
          fr%at(n)%is = c%atcel(id)%is
       end if
    end do
    fr%nat = n
    call realloc(fr%at,n)

  end function identify_fragment

  !> Identify a fragment in the unit cell from an external xyz
  !> file. An instance of a fragment object is returned. If any of the
  !> atoms is not correctly identified, return 0.
  module function identify_fragment_from_xyz(c,file,ti) result(fr)
    use tools_io, only: fopen_read, string, ferror, faterr, fclose
    use param, only: bohrtoa, icrd_cart, mlen
    use types, only: realloc

    class(crystal), intent(in) :: c
    character*(*) :: file
    type(thread_info), intent(in), optional :: ti
    type(fragment) :: fr

    integer :: lu, nat
    integer :: id, i
    real*8 :: x0(3)
    character(len=mlen) :: word

    lu = fopen_read(file,ti=ti)
    read(lu,*,err=999) nat
    read(lu,*,err=999)

    fr%nat = nat
    allocate(fr%at(nat))
    fr%nspc = c%nspc
    allocate(fr%spc(fr%nspc))
    fr%spc = c%spc
    do i = 1, nat
       read(lu,*,err=999) word, x0
       x0 = x0 / bohrtoa - c%molx0
       id = c%identify_atom(x0,icrd_cart)
       if (id == 0) then
          fr%nat = 0
          deallocate(fr%at)
          fr%nspc = 0
          deallocate(fr%spc)
          return
       endif
       fr%at(i)%r = x0
       fr%at(i)%x = c%c2x(x0)
       fr%at(i)%cidx = id
       fr%at(i)%idx = c%atcel(id)%idx
       fr%at(i)%lvec = nint(fr%at(i)%x - c%atcel(id)%x)
       fr%at(i)%is = c%atcel(id)%is
    end do
    call fclose(lu)
    call realloc(fr%at,fr%nat)

    return
999 continue
    call ferror('identify_fragment_from_xyz','error reading xyz file: '//string(file),faterr)

  end function identify_fragment_from_xyz

  !> Obtain symmetry equivalent positions of xp0 and write them to
  !> vec. Write the multiplicity of the xp0 position to mmult. xp0 and
  !> vec are in crystallographic coordinates. irotm and icenv contain
  !> the index of the rotation matrix and centering vectors
  !> responsible for the transformation of xp into the corresponding
  !> vector in vec. eps is the minimum distance to consider two points
  !> equivalent (in bohr). vec, irotm, icenv, and eps0 are optional.
  module subroutine symeqv(c,xp0,mmult,vec,irotm,icenv,eps0)
    use types, only: realloc
    class(crystal), intent(in) :: c !< Input crystal
    real*8, intent(in) :: xp0(3) !< input position (crystallographic)
    integer, intent(out) :: mmult !< multiplicity
    real*8, allocatable, intent(inout), optional :: vec(:,:) !< sym-eq positions (crystallographic)
    integer, allocatable, intent(inout), optional :: irotm(:) !< index of the operation
    integer, allocatable, intent(inout), optional :: icenv(:) !< index of the cent. vector
    real*8, intent(in), optional :: eps0 !< Minimum distance to consider two vectors different (bohr)

    real*8 :: avec(3,c%neqv*c%ncv)
    integer :: arotm(c%neqv*c%ncv)
    integer :: acenv(c%neqv*c%ncv)
    integer :: i, j
    integer :: mrot, mrot0
    real*8 :: xp(3)
    real*8 :: eps

    real*8, parameter :: eps_default = 1d-2

    ! the eps
    if (present(eps0)) then
       eps = eps0
    else
       eps = eps_default
    end if

    !.Translate position to main cell
    xp = xp0 - floor(xp0)

    !.Run over symmetry operations and create the list of copies
    mrot0 = 0
    do i = 1, c%neqv
       do  j = 1, c%ncv
          mrot0 = mrot0 + 1
          avec(:,mrot0) = matmul(c%rotm(1:3,1:3,i),xp) + c%rotm(:,4,i) + c%cen(:,j)
          avec(:,mrot0) = avec(:,mrot0) - floor(avec(:,mrot0))
          arotm(mrot0) = i
          acenv(mrot0) = j
       enddo
    enddo

    ! calculate distances and (possibly) write vec
    if (present(vec)) call realloc(vec,3,mrot0)
    if (present(irotm)) call realloc(irotm,mrot0)
    if (present(icenv)) call realloc(icenv,mrot0)
    mrot = 0
    d: do i = 1, mrot0
       do j = 1,i-1
          if (c%are_lclose(avec(:,i),avec(:,j),eps)) cycle d
       end do
       mrot = mrot + 1
       if (present(vec)) vec(:,mrot) = avec(:,i)
       if (present(irotm)) irotm(mrot) = arotm(i)
       if (present(icenv)) icenv(mrot) = acenv(i)
    end do d

    mmult=mrot
    if(present(vec)) call realloc(vec,3,mmult)
    if(present(irotm)) call realloc(irotm,mmult)
    if(present(icenv)) call realloc(icenv,mmult)

  end subroutine symeqv

  !> Calculate the multiplicity of the point x0 (cryst. coord.)
  module function get_mult(c,x0) result (mult)
    class(crystal), intent(in) :: c
    real*8, intent(in) :: x0(3)
    integer :: mult

    call c%symeqv(x0,mult)

  end function get_mult

  !> Calculate the number of k-points (nk) for a given rk-length. Uses
  !> the VASP formula.
  module subroutine get_kpoints(c,rk,nk)
    class(crystal), intent(in) :: c
    real*8, intent(in) :: rk
    integer, intent(out) :: nk(3)

    nk = int(max(1d0,rk * c%ar + 0.5d0))

  end subroutine get_kpoints

  !> Calculate the distance matrix (molecules only). If inverse,
  !> invert the distances. If conn, return the distance if the atoms
  !> are bonded, and the inverse distance if they are not, and zero
  !> along the diagonals.
  module subroutine distmatrix(c,d,inverse,conn)
    class(crystal), intent(in) :: c
    real*8, allocatable, intent(inout) :: d(:,:)
    logical, intent(in), optional :: inverse
    logical, intent(in), optional :: conn

    integer :: i, j
    logical :: inv_, conn_

    inv_ = .false.
    conn_ = .false.
    if (present(inverse)) inv_ = inverse
    if (present(conn)) conn_ = conn

    if (.not.c%ismolecule) return
    if (allocated(d)) deallocate(d)
    allocate(d(c%ncel,c%ncel))

    if (conn_) then
       d = 0.d0
       do i = 1, c%ncel
          do j = i+1, c%ncel
             d(i,j) = 1d0/norm2(c%atcel(i)%r - c%atcel(j)%r)
             d(j,i) = d(i,j)
          end do
       end do
       do i = 1, c%ncel
          do j = 1, c%nstar(i)%ncon
             d(i,c%nstar(i)%idcon(j)) = norm2(c%atcel(i)%r - c%atcel(c%nstar(i)%idcon(j))%r)
             d(c%nstar(i)%idcon(j),i) = d(i,c%nstar(i)%idcon(j))
          end do
       end do
    else
       d = 0d0
       do i = 1, c%ncel
          do j = i+1, c%ncel
             if (inv_) then
                d(i,j) = 1d0/norm2(c%atcel(i)%r - c%atcel(j)%r)
             else
                d(i,j) = norm2(c%atcel(i)%r - c%atcel(j)%r)
             end if
             d(j,i) = d(i,j)
          end do
       end do
    end if

  end subroutine distmatrix

  !> List atoms in a number of cells around the main cell (nx cells),
  !> possibly with border (doborder).
  module function listatoms_cells(c,nx,doborder) result(fr)
    use types, only: realloc
    class(crystal), intent(in) :: c
    integer, intent(in) :: nx(3)
    logical, intent(in) :: doborder
    type(fragment) :: fr

    real*8, parameter :: rthr = 0.01d0
    real*8, parameter :: rthr1 = 1-rthr

    integer :: ix, iy, iz, i
    logical :: if1

    fr%nspc = c%nspc
    allocate(fr%spc(c%nspc))
    fr%spc(1:fr%nspc) = c%spc(1:c%nspc)
    allocate(fr%at(1))
    fr%nat = 0

    ! All atoms in these cells
    fr%nat = 0
    do ix = 0,nx(1)-1
       do iy = 0,nx(2)-1
          do iz = 0,nx(3)-1
             do i = 1, c%ncel
                fr%nat = fr%nat + 1
                if (fr%nat > size(fr%at)) call realloc(fr%at,2*fr%nat)
                fr%at(fr%nat)%x = c%atcel(i)%x + (/ix,iy,iz/)
                fr%at(fr%nat)%r = c%x2c(fr%at(fr%nat)%x)
                fr%at(fr%nat)%cidx = i
                fr%at(fr%nat)%idx = c%atcel(i)%idx
                fr%at(fr%nat)%lvec = (/ix,iy,iz/)
                fr%at(fr%nat)%is = c%atcel(i)%is
             end do
          end do
       end do
    end do

    ! border: pick atoms
    if (doborder) then
       do ix = -1,nx(1)
          do iy = -1,nx(2)
             do iz = -1,nx(3)
                if (ix > -1 .and. ix < nx(1) .and. iy > -1 .and. iy < nx(2) .and.&
                   iz > -1 .and. iz < nx(3)) cycle
                do i = 1, c%ncel
                   ! border
                   if1 = (ix == -1 .and. c%atcel(i)%x(1)<rthr1 .or. &
                      ix == nx(1) .and. c%atcel(i)%x(1)>rthr .or. &
                      iy == -1 .and. c%atcel(i)%x(2)<rthr1 .or. &
                      iy == nx(2) .and. c%atcel(i)%x(2)>rthr .or. &
                      iz == -1 .and. c%atcel(i)%x(3)<rthr1 .or. &
                      iz == nx(3) .and. c%atcel(i)%x(3)>rthr)
                   if (.not.if1) then
                      fr%nat = fr%nat + 1
                      if (fr%nat > size(fr%at)) call realloc(fr%at,2*fr%nat)
                      fr%at(fr%nat)%x = c%atcel(i)%x + (/ix,iy,iz/)
                      fr%at(fr%nat)%r = c%x2c(fr%at(fr%nat)%x)
                      fr%at(fr%nat)%cidx = i
                      fr%at(fr%nat)%idx = c%atcel(i)%idx
                      fr%at(fr%nat)%lvec = (/ix,iy,iz/)
                      fr%at(fr%nat)%is = c%atcel(i)%is
                      cycle
                   end if
                end do
             end do
          end do
       end do
    end if
    call realloc(fr%at,fr%nat)

  end function listatoms_cells

  !> List atoms inside a sphere of radius rsph and center xsph
  !> (cryst.)  or a cube of side rcub and center xcub (cryst.). Return
  !> the list of atomic positions (Cartesian) in x, the atomic numbers
  !> in z and the number of atoms in nat.
  module function listatoms_sphcub(c,rsph,xsph,rcub,xcub) result(fr)
    use tools_io, only: ferror, faterr
    use types, only: realloc
    class(crystal), intent(in) :: c
    real*8, intent(in), optional :: rsph, xsph(3)
    real*8, intent(in), optional :: rcub, xcub(3)
    type(fragment) :: fr

    integer :: ix, iy, iz, i, nn
    real*8 :: x0(3), d
    logical :: doagain, dosph

    if (.not.(present(rsph).and.present(xsph)).and..not.(present(rcub).and.present(xcub))) &
       call ferror("listatoms_sphcub","Need sphere or cube input",faterr)
    dosph = present(rsph)

    allocate(fr%at(1))
    fr%nat = 0
    allocate(fr%spc(c%nspc))
    fr%nspc = c%nspc
    fr%spc = c%spc

    ! all atoms in a sphere
    doagain = .true.
    nn = -1
    do while(doagain)
       doagain = .false.
       nn = nn + 1
       do ix = -nn, nn
          do iy = -nn, nn
             do iz = -nn, nn
                if (abs(ix)/=nn .and. abs(iy)/=nn .and. abs(iz)/=nn) cycle
                do i = 1, c%ncel
                   if (dosph) then
                      x0 = c%x2c(c%atcel(i)%x + (/ix,iy,iz/) - xsph)
                      if (all(abs(x0) > rsph)) cycle
                      d = norm2(x0)
                      if (d >= rsph) cycle
                   else
                      x0 = c%x2c(c%atcel(i)%x + (/ix,iy,iz/) - xcub)
                      if (any(abs(x0) > rcub)) cycle
                   endif

                   ! add this atom
                   fr%nat = fr%nat + 1
                   if (fr%nat > size(fr%at)) call realloc(fr%at,2*fr%nat)
                   fr%at(fr%nat)%x = c%atcel(i)%x + (/ix,iy,iz/)
                   fr%at(fr%nat)%r = c%x2c(fr%at(fr%nat)%x)
                   fr%at(fr%nat)%cidx = i
                   fr%at(fr%nat)%idx = c%atcel(i)%idx
                   fr%at(fr%nat)%lvec = (/ix,iy,iz/)
                   fr%at(fr%nat)%is = c%atcel(i)%is
                   doagain = .true.
                end do
             end do
          end do
       end do
    end do
    call realloc(fr%at,fr%nat)

  end function listatoms_sphcub

  !> Using the calculated asterisms for each atom determine the
  !> molecules in the system by finding the connected components of
  !> the graph. This routine also determines whether the crystal is
  !> extended or molecular. Fills c%nmol and c%mol.
  module subroutine fill_molecular_fragments(c)
    use fragmentmod, only: realloc_fragment
    use tools_math, only: lattice_direction
    use tools_io, only: ferror, faterr
    use tools, only: qcksort
    use types, only: realloc
    class(crystal), intent(inout) :: c

    integer :: i, j, k, l, jid, newid, newl(3)
    integer :: nat, nlvec, lwork, info
    logical :: found, fdisc
    integer, allocatable :: id(:), lvec(:,:), iord(:), icidx(:)
    logical, allocatable :: ldone(:), used(:)
    real*8, allocatable :: rlvec(:,:), sigma(:), uvec(:,:), vvec(:,:), work(:)
    real*8 :: xcm(3)

    if (.not.allocated(c%nstar)) &
       call ferror('fill_molecular_fragments','no asterisms found',faterr)
    if (allocated(c%mol)) deallocate(c%mol)
    if (allocated(c%idatcelmol)) deallocate(c%idatcelmol)

    ! initizialize
    allocate(used(c%ncel))
    used = .false.
    c%nmol = 0
    allocate(c%mol(1),id(10),lvec(3,10),ldone(10))

    ! run over atoms in the unit cell
    do i = 1, c%ncel
       if (used(i)) cycle

       ! increment the fragment counter
       c%nmol = c%nmol + 1
       if (c%nmol > size(c%mol)) then
          call realloc_fragment(c%mol,2*c%nmol)
       end if
       c%mol(c%nmol)%discrete = .true.
       c%mol(c%nmol)%nlvec = 0

       ! initialize the stack with atom i in the seed
       nat = 1
       id(1) = i
       lvec(:,1) = 0d0
       ldone(1) = .false.
       ! run the stack
       do while (.not.all(ldone(1:nat)))
          ! find the next atom that is not done
          do j = 1, nat
             if (.not.ldone(j)) exit
          end do
          ldone(j) = .true.
          jid = id(j)

          ! run over all neighbors of j
          do k = 1, c%nstar(jid)%ncon
             ! id for the atom and lattice vector
             newid = c%nstar(jid)%idcon(k)
             newl = c%nstar(jid)%lcon(:,k) + lvec(:,j)

             ! Is this atom in the fragment already? -> skip it. If it
             ! has a different lattice vector, mark the fragment as
             ! not discrete.
             found = .false.
             do l = 1, nat
                found = (newid == id(l))
                fdisc = all(newl == lvec(:,l))
                if (found) exit
             end do
             if (found) then
                if (.not.fdisc.and..not.c%ismolecule) then
                   if (.not.allocated(c%mol(c%nmol)%lvec)) allocate(c%mol(c%nmol)%lvec(3,1))

                   newl = newl - lvec(:,l)
                   found = .false.
                   do l = 1, c%mol(c%nmol)%nlvec
                      if (all(c%mol(c%nmol)%lvec(:,l) == newl) .or. all(c%mol(c%nmol)%lvec(:,l) == -newl)) then
                         found = .true.
                         exit
                      end if
                   end do
                   if (.not.found) then
                      c%mol(c%nmol)%nlvec = c%mol(c%nmol)%nlvec + 1
                      if (c%mol(c%nmol)%nlvec > size(c%mol(c%nmol)%lvec,2)) then
                         call realloc(c%mol(c%nmol)%lvec,3,2*c%mol(c%nmol)%nlvec)
                      end if
                      c%mol(c%nmol)%lvec(:,c%mol(c%nmol)%nlvec) = newl
                   end if

                   c%mol(c%nmol)%discrete = .false.
                end if
                cycle
             end if

             ! Have we used this atom already?
             if (used(newid)) cycle

             ! Add the atom to the stack and mark it as used.
             nat = nat + 1
             if (nat > size(ldone)) then
                call realloc(id,2*nat)
                call realloc(lvec,3,2*nat)
                call realloc(ldone,2*nat)
             end if
             id(nat) = newid
             lvec(:,nat) = newl
             used(newid) = .true.
             ldone(nat) = .false.
          end do
       end do

       ! add this fragment to the list
       used(i) = .true.
       allocate(c%mol(c%nmol)%spc(c%nspc))
       c%mol(c%nmol)%nspc = c%nspc
       c%mol(c%nmol)%spc = c%spc(1:c%nspc)
       allocate(c%mol(c%nmol)%at(nat))
       c%mol(c%nmol)%nat = nat
       do j = 1, nat
          c%mol(c%nmol)%at(j)%x = c%atcel(id(j))%x + lvec(:,j)
          c%mol(c%nmol)%at(j)%r = c%x2c(c%mol(c%nmol)%at(j)%x)
          c%mol(c%nmol)%at(j)%cidx = id(j)
          c%mol(c%nmol)%at(j)%idx = c%atcel(id(j))%idx
          c%mol(c%nmol)%at(j)%lvec = lvec(:,j)
          c%mol(c%nmol)%at(j)%is = c%atcel(id(j))%is
       end do
       if (c%mol(c%nmol)%nlvec > 0) &
          call realloc(c%mol(c%nmol)%lvec,3,c%mol(c%nmol)%nlvec)

       ! reorder the fragment in order of increasing cidx
       allocate(iord(nat),icidx(nat))
       do j = 1, nat
          iord(j) = j
          icidx(j) = c%mol(c%nmol)%at(j)%cidx
       end do
       call qcksort(icidx,iord,1,nat)
       c%mol(c%nmol)%at = c%mol(c%nmol)%at(iord)
       deallocate(iord,icidx)
    end do
    call realloc_fragment(c%mol,c%nmol)
    deallocate(used)

    ! translate all fragments to the main cell
    if (.not.c%ismolecule) then
       do i = 1, c%nmol
          xcm = c%mol(i)%cmass()
          newl = floor(c%c2x(xcm))
          do j = 1, c%mol(i)%nat
             c%mol(i)%at(j)%x = c%mol(i)%at(j)%x - newl
             c%mol(i)%at(j)%r = c%x2c(c%mol(i)%at(j)%x)
             c%mol(i)%at(j)%lvec = c%mol(i)%at(j)%lvec - newl
          end do
       end do
    end if

    ! accumulate all the lattice vectors and calculate the vacuum lattice vectors
    nlvec = 0
    c%nlvac = 0
    do i = 1, c%nmol
       nlvec = nlvec + c%mol(i)%nlvec
    end do
    if (nlvec > 0) then
       allocate(rlvec(3,nlvec))
       k = 0
       do i = 1, c%nmol
          do j = 1, c%mol(i)%nlvec
             k = k + 1
             rlvec(:,k) = c%mol(i)%lvec(:,j)
          end do
       end do

       lwork = 5*3 + max(nlvec,3)
       allocate(sigma(3),uvec(3,3),vvec(1,1),work(lwork))
       sigma = 0d0
       call dgesvd('A','N',3,k,rlvec,3,sigma,uvec,3,vvec,1,work,lwork,info)
       if (info /= 0) &
          call ferror("fill_molecular_fragments","dgesvd failed!",faterr)

       c%nlvac = count(abs(sigma) < 1d-12)
       c%lvac = 0d0
       c%lcon = 0d0
       if (c%nlvac == 2) then
          c%lvac(:,1) = lattice_direction(uvec(:,3),.true.)
          c%lvac(:,2) = lattice_direction(uvec(:,2),.true.)
          c%lcon(:,1) = lattice_direction(uvec(:,1),.true.)
       elseif (c%nlvac == 1) then
          c%lvac(:,1) = lattice_direction(uvec(:,3),.true.)
          c%lcon(:,1) = lattice_direction(uvec(:,2),.true.)
          c%lcon(:,2) = lattice_direction(uvec(:,1),.true.)
       end if
       deallocate(rlvec,sigma,uvec,vvec,work)
    end if

    ! fill the mapping between atoms and molecules
    if (allocated(c%idatcelmol)) deallocate(c%idatcelmol)
    allocate(c%idatcelmol(c%ncel))
    c%idatcelmol = 0
    do i = 1, c%nmol
       do j = 1, c%mol(i)%nat
          c%idatcelmol(c%mol(i)%at(j)%cidx) = i
       end do
    end do

  end subroutine fill_molecular_fragments

  !> Calculate whether the crystal is a molecular crystal. If it is,
  !> for each molecule, calculate whether it is contained as a whole
  !> in the asymmetric unit (idxmol=-1) or not (idxmol >= 0). If it
  !> is not, then the molecule can be symmetry-unique idxmol = 0 or
  !> equivalent to the molecule with index idxmol>0. Fills ismol3d
  !> and idxmol.
  module subroutine calculate_molecular_equivalence(c)
    use tools, only: qcksort
    class(crystal), intent(inout) :: c

    integer :: mmax, i, j, nat
    integer, allocatable :: midx(:,:)

    ! 3d molecular crystals calculations
    c%ismol3d = all(c%mol(1:c%nmol)%discrete).and..not.c%ismolecule
    if (c%ismol3d) then
       ! assign fractional (-1)/symmetry unique (0)/symmetry equivalent (>0)
       if (allocated(c%idxmol)) deallocate(c%idxmol)
       allocate(c%idxmol(c%nmol))
       c%idxmol = 0

       ! find the maximum number of atoms
       mmax = 0
       do i = 1, c%nmol
          mmax = max(mmax,c%mol(i)%nat)
       end do

       ! sort the atomic indices
       allocate(midx(mmax,c%nmol))
       do i = 1, c%nmol
          nat = c%mol(i)%nat
          do j = 1, nat
             midx(j,i) = c%mol(i)%at(j)%idx
          end do
          call qcksort(midx(1:nat,i))

          ! check if the molecule is fractional
          do j = 1, c%mol(i)%nat-1
             if (midx(j,i) == midx(j+1,i)) then
                c%idxmol(i) = -1
                exit
             end if
          end do
       end do

       ! assign idxmol based on the atomic sequence
       do i = 1, c%nmol
          nat = c%mol(i)%nat
          if (c%idxmol(i) < 0) cycle
          do j = 1, i-1
             if (c%idxmol(j) < 0) cycle
             if (all(midx(1:nat,i) == midx(1:nat,j))) then
                c%idxmol(i) = j
                exit
             end if
          end do
       end do
       deallocate(midx)
    end if

  end subroutine calculate_molecular_equivalence

  !> List all molecules resulting from completing the initial fragment
  !> fri by adding adjacent atoms that are covalently bonded. Return
  !> the number of fragment (nfrag), the fragments themselves (fr),
  !> and whether the fragments are discrete (not connected to copies
  !> of themselves in a different cell).
  module subroutine listmolecules(c,fri,nfrag,fr,isdiscrete)
    use fragmentmod, only: realloc_fragment
    use types, only: realloc
    class(crystal), intent(inout) :: c
    type(fragment), intent(in) :: fri
    integer, intent(out) :: nfrag
    type(fragment), intent(out), allocatable :: fr(:)
    logical, intent(out), allocatable :: isdiscrete(:)

    integer :: i, j, k, l, newid, newl(3), jid
    integer :: nat
    integer, allocatable :: id(:), lvec(:,:)
    logical, allocatable :: ldone(:)
    logical :: found, ldist
    integer :: nseed
    integer, allocatable :: idseed(:), lseed(:,:)
    logical, allocatable :: fseed(:)

    ! unwrap the input fragment
    nseed = fri%nat
    allocate(idseed(nseed),lseed(3,nseed),fseed(nseed))
    do j = 1, nseed
       idseed(j) = fri%at(j)%cidx
       lseed(:,j) = fri%at(j)%lvec
    end do
    fseed = .false.

    ! allocate stuff
    nfrag = 0
    allocate(fr(1),isdiscrete(1),id(10),lvec(3,10),ldone(10))
    isdiscrete = .true.

    do i = 1, nseed
        if (fseed(i)) cycle

       ! initialize the stack with atom i in the seed
       nat = 1
       id(1) = idseed(i)
       lvec(:,1) = lseed(:,i)
       ldone(1) = .false.
       ldist = .true.
       ! run the stack
       do while (.not.all(ldone(1:nat)))
          ! find the next atom that is not done
          do j = 1, nat
             if (.not.ldone(j)) exit
          end do
          ldone(j) = .true.
          jid = id(j)

          ! run over all neighbors of j
          do k = 1, c%nstar(jid)%ncon
             ! id for the atom and lattice vector
             newid = c%nstar(jid)%idcon(k)
             newl = c%nstar(jid)%lcon(:,k) + lvec(:,j)

             ! is this atom in the fragment already? -> skip it
             found = .false.
             do l = 1, nat
                found = (newid == id(l)) .and. all(newl == lvec(:,l))
                if (found) exit
             end do
             if (found) cycle

             ! is this atom in the fragment already with a different
             ! lattice vector?  -> add it to the list but not to the
             ! stack, and mark the fragment as non-discrete
             found = .false.
             do l = 1, nat
                found = (newid == id(l))
                if (found) exit
             end do
             nat = nat + 1
             if (nat > size(ldone)) then
                call realloc(id,2*nat)
                call realloc(lvec,3,2*nat)
                call realloc(ldone,2*nat)
             end if
             id(nat) = newid
             lvec(:,nat) = newl

             if (found) then
                ldone(nat) = .true.
                ldist = .false.
             else
                ! if it wasn't found, then add the atom to the stack
                ldone(nat) = .false.
             end if
          end do
       end do

       ! add this fragment to the list
       fseed(i) = .true.
       nfrag = nfrag + 1
       if (nfrag > size(fr)) then
          call realloc_fragment(fr,2*nfrag)
          call realloc(isdiscrete,2*nfrag)
       end if
       fr(nfrag)%nat = nat
       allocate(fr(nfrag)%at(nat))
       if (allocated(fr(nfrag)%spc)) deallocate(fr(nfrag)%spc)
       fr(nfrag)%nspc = c%nspc
       allocate(fr(nfrag)%spc(c%nspc))
       fr(nfrag)%spc = c%spc
       isdiscrete(nfrag) = ldist
       do j = 1, nat
          fr(nfrag)%at(j)%x = c%atcel(id(j))%x + lvec(:,j)
          fr(nfrag)%at(j)%r = c%x2c(fr(nfrag)%at(j)%x)
          fr(nfrag)%at(j)%cidx = id(j)
          fr(nfrag)%at(j)%idx = c%atcel(id(j))%idx
          fr(nfrag)%at(j)%lvec = lvec(:,j)
          fr(nfrag)%at(j)%is = c%atcel(id(j))%is
       end do

       ! run over all atoms in the new fragment and mark those atoms in the seed
       do j = 1, nat
          do k = 1, nseed
             if (fseed(k)) cycle
             if (id(j) == idseed(k) .and. all(lvec(:,j) == lseed(:,k))) &
                fseed(k) = .true.
          end do
       end do
    end do
    call realloc_fragment(fr,nfrag)
    call realloc(isdiscrete,nfrag)

  end subroutine listmolecules

  !> Determines the site symmetry for a point x0 in cryst. coords.
  !> Two points are the same if their distance is less than eps0.
  !> Returns the site symmetry group symbol (sitesymm), the
  !> number of operations in this group (leqv) and the rotation
  !> operations (lrotm)
  module function sitesymm(c,x0,eps0,leqv,lrotm)
    use tools_io, only: string
    use param, only: eye
    class(crystal), intent(in) :: c !< Input crystal
    real*8, intent(in) :: x0(3) !< Input point in cryst. coords.
    real*8, intent(in), optional :: eps0 !< Two points are different if distance is > eps
    character*3 :: sitesymm !< point group symbol
    integer, optional :: leqv !< Number of operations in the group
    real*8, optional :: lrotm(3,3,48) !< Point group operations

    integer :: i, m
    real*8 :: dumy(3), dist2, eps, vec(3)
    integer :: type
    logical :: ok
    integer :: highest, highests
    integer :: nnsym, order
    integer :: ordersym(c%neqv), masksym(0:9)

    real*8, parameter :: eps_default = 1d-2

    if (c%ismolecule .or. c%havesym == 0 .or. (c%neqv == 1 .and. c%ncv == 1)) then
       sitesymm='C1'
       if (present(leqv).and.present(lrotm)) then
          leqv = 1
          lrotm(:,:,1) = eye
       end if
       return
    end if

    if (present(eps0)) then
       eps = eps0
    else
       eps = eps_default
    end if

    ! Run over all proper symmetry elements of symmetry
    nnsym = 0
    masksym = 0
    ordersym = 0
    do i = 1, c%neqv
       ok = .false.
       do m = 1, c%ncv
          dumy = matmul(c%rotm(1:3,1:3,i),x0)
          dumy = - dumy + x0 - c%rotm(1:3,4,i) - c%cen(:,m)
          call c%shortest(dumy,dist2)
          ok = (dist2 < eps)
          if (ok) exit
       end do
       if (ok) then
          ! A symmetry operation at location x has been found.
          nnsym = nnsym+1
          call typeop(c%rotm(:,:,i),type,vec,order)
          ordersym(nnsym) = order
          masksym(type) = masksym(type) + 1
          if (present(leqv).and.present(lrotm)) then
             leqv = nnsym
             lrotm(:,:,leqv) = c%rotm(1:3,1:3,i)
          end if
       endif
    enddo

    ! calculate the point group
    sitesymm = ""
    if (masksym(c3) > 2) then
       ! cubic groups
       if (masksym(c4) /= 0) then
          if (masksym(inv) /= 0) then
             sitesymm='Oh'
          else
             sitesymm='O'
          endif
       elseif (masksym(s4).ne.0) then
          sitesymm= 'Td'
       elseif (masksym(inv).ne.0) then
          sitesymm = 'Th'
       else
          sitesymm = 'T'
       endif
    else
       !Compute highest order proper axis.
       highest=0
       highests=0
       if (masksym(c2) /= 0) highest=2
       if (masksym(c3) /= 0) highest=3
       if (masksym(c4) /= 0) highest=4
       if (masksym(c6) /= 0) highest=6
       if (masksym(s3) /= 0) highests=3
       if (masksym(s4) /= 0) highests=4
       if (masksym(s6) /= 0) highests=6
       if (highest == 0) then
          if (masksym(inv) /= 0) then
             sitesymm='i'
          elseif (masksym(sigma) /= 0) then
             sitesymm='Cs'
          else
             sitesymm='C1'
          endif
       elseif (masksym(c2) >= highest) then
          if (masksym(sigma) == 0) then
             sitesymm='D' // string(highest)
          elseif (masksym(inv) .eq. 1) then
             if (highest == 3) then
                sitesymm= 'D3d'
             else
                sitesymm= 'D' // string(highest) // 'h'
             endif
          else
             if (highest .eq. 3) then
                sitesymm= 'D3h'
             else
                sitesymm= 'D' // string(highest) // 'd'
             endif
          endif
       elseif (masksym(sigma) == 0) then
          if (highests /= 0) then
             sitesymm= 'S' // string(highests/2)
          else
             sitesymm= 'C' // string(highest)
          endif
       elseif (masksym(sigma) .lt. highest) then
          sitesymm= 'C' // string(highest) // 'h'
       else
          sitesymm= 'C' // string(highest) // 'v'
       endif
    endif

  end function sitesymm

  !> Calculate the packing ratio (in %) using the nearest-neighbor
  !> information. Each atom is assigned a ratio equal to half the distance
  !> to its nearest neighbor.
  module function get_pack_ratio(c) result(px)
    use param, only: pi
    class(crystal), intent(inout) :: c
    real*8 :: px

    integer :: i

    px = 0d0
    do i = 1, c%nneq
       px = px + c%at(i)%mult * 4d0/3d0 * pi * c%at(i)%rnn2**3
    end do
    px = px / c%omega * 100d0

  end function get_pack_ratio

  !> Calculate the vdw volume in a molecule or crystal by Monte-Carlo
  !> sampling.  relerr = use enough points to obtain a standard
  !> deviation divided by the volume equal to this value. If
  !> rtable(1:maxzat0) is present, use those radii instead of the
  !> van der walls radii
  module function vdw_volume(c,relerr,rtable) result(vvdw)
    use tools_io, only: ferror, faterr
    use param, only: VBIG, atmvdw, icrd_cart
    class(crystal), intent(inout) :: c
    real*8, intent(in) :: relerr
    real*8, intent(in), optional :: rtable(:)
    real*8 :: vvdw

    real*8 :: xmin(3), xmax(3), x(3), vtot, svol, pp
    integer :: i, nat
    real*8, allocatable :: rvdw(:,:)
    integer :: ierr
    integer*8 :: nin, ntot
    logical :: again

    ! build the list of atomic radii
    allocate(rvdw(c%nspc,2))
    do i = 1, c%nspc
       rvdw(i,1) = 0d0
       if (present(rtable)) then
          rvdw(i,2) = rtable(c%spc(i)%z)
       else
          rvdw(i,2) = atmvdw(c%spc(i)%z)
       end if
    end do

    ! calculate the encompassing box
    if (c%ismolecule) then
       xmin = VBIG
       xmax = -VBIG
       do i = 1, c%ncel
          xmin = min(xmin,c%atcel(i)%r)
          xmax = max(xmax,c%atcel(i)%r)
       end do
       xmin = xmin - maxval(rvdw(:,2))
       xmax = xmax + maxval(rvdw(:,2))
       vtot = product(xmax-xmin)
    else
       vtot = c%omega
    end if

    ! use Monte-Carlo to determine the volume
    again = .true.
    ntot = 0
    nin = 0
    do while (again)
       ntot = ntot + 1
       call random_number(x)

       if (c%ismolecule) then
          x = xmin + x * (xmax - xmin)
       else
          x = c%x2c(x)
       end if
       call c%env%list_near_atoms(x,icrd_cart,.false.,nat,ierr,up2dsp=rvdw)
       if (ierr > 0) then
          call ferror('vdw_volume','error determining the list of near atoms',faterr)
       elseif (nat > 0) then
          nin = nin + 1
       end if

       pp = real(nin,8) / real(ntot,8)
       svol = vtot * sqrt(pp * (1-pp)/real(ntot,8))
       vvdw = vtot * pp
       if (ntot > 100) then
          again = (svol > relerr * vvdw)
       end if
    end do

  end function vdw_volume

  !> Calculate the powder diffraction pattern.
  !> On input, npts is the number of 2*theta points from the initial
  !> (th2ini0) to the final (th2end0) 2*theta values. Both angles are
  !> in degrees. lambda0 is the wavelength of the radiation (in
  !> angstrom). sigma is the parameter for the Gaussian broadening.
  !> fpol is the polarization correction factor (0 = unpolarized, 0.95
  !> = syncrhotron). If ishard=.false., represent the tails of the peaks
  !> outside the plot range.
  !> On output, t is the 2*theta grid, ih is the intensity on the
  !> 2*theta grid, th2p is the 2*theta for the located maxima, ip is
  !> the list of maxima itensities, and hvecp is the reciprocal
  !> lattice vector corresponding to the peaks.
  module subroutine powder(c,th2ini0,th2end0,ishard,npts,lambda0,fpol,&
     sigma,t,ih,th2p,ip,hvecp)
    use param, only: pi, bohrtoa, cscatt, c2scatt
    use tools_io, only: ferror, faterr
    use tools, only: qcksort
    use types, only: realloc
    class(crystal), intent(in) :: c
    real*8, intent(in) :: th2ini0, th2end0
    logical, intent(in) :: ishard
    integer, intent(in) :: npts
    real*8, intent(in) :: lambda0
    real*8, intent(in) :: fpol
    real*8, intent(in) :: sigma
    real*8, allocatable, intent(inout) :: t(:)
    real*8, allocatable, intent(inout) :: ih(:)
    real*8, allocatable, intent(inout) :: th2p(:)
    real*8, allocatable, intent(inout) :: ip(:)
    integer, allocatable, intent(inout) :: hvecp(:,:)

    integer :: i, np, hcell, h, k, l, iz, idx
    real*8 :: th2ini, th2end, lambda, hvec(3), kvec(3), th, sth, th2
    real*8 :: sigma2, smax, dh2, dh, dh3, sthlam, cterm, sterm
    real*8 :: ffac, as(4), bs(4), cs, c2s(4), int, mcorr, afac
    real*8 :: ipmax, ihmax, tshift
    integer :: hmax
    logical :: again
    integer, allocatable :: multp(:)
    integer, allocatable :: io(:)
    real*8, allocatable :: isum(:)

    integer, parameter :: mp = 20
    real*8, parameter :: ieps = 1d-5
    real*8, parameter :: theps = 1d-5
    real*8, parameter :: iepscont = 1d-10
    real*8, parameter :: intmax = 1d15

    ! calculate tshift
    if (.not.ishard) then
       tshift = sigma * sqrt(abs(-2d0 * log(iepscont/intmax)))
    else
       tshift = 0d0
    end if
    tshift = tshift * pi / 180d0

    ! prepare the grid limits
    if (allocated(t)) deallocate(t)
    if (allocated(ih)) deallocate(ih)
    allocate(t(npts),ih(npts),isum(npts))
    do i = 1, npts
       t(i) = th2ini0 + real(i-1,8) / real(npts-1,8) * (th2end0-th2ini0)
    end do
    ih = 0d0
    th2ini = th2ini0 * pi / 180d0
    th2end = th2end0 * pi / 180d0

    ! allocate for peak list
    if (allocated(th2p)) deallocate(th2p)
    if (allocated(ip)) deallocate(ip)
    if (allocated(hvecp)) deallocate(hvecp)
    allocate(th2p(mp),ip(mp),multp(mp),hvecp(3,mp))

    ! cell limits, convert lambda to bohr
    lambda = lambda0 / bohrtoa
    smax = sin((th2end+tshift)/2d0)
    hmax = 2*ceiling(2*smax/lambda/minval(c%ar))
    ! broadening -> gaussian
    sigma2 = sigma * sigma

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

                ! profile correction
                ! Yinghua J. Appl. Cryst. 20 (1987) 258
                ! lpf = (1 + cos(th2)**2) / sin(th)**2
                ! int = int * lpf

                ! gdis lorentz-polarization correction; checked in diamond;
                ! criticized by Yinghua because it is a correction for the integrated
                ! intensity
                ! lpf = 0.5*(1+cos(th2)**2) / sin(th2) / sin(th)
                ! int = int * lpf

                ! FoX-compatible
                ! lorentz correction
                mcorr = 1d0 / sin(th2)
                ! slit aperture
                mcorr = mcorr / sin(th)
                ! polarization
                afac = (1-fpol) / (1+fpol)
                mcorr = mcorr * (1+afac*(0.5d0+0.5d0*cos(2*th2))) / (1+afac)
                int = int * mcorr

                ! sum the peak
                if (int > ieps) then
                   again = .true.
                   ! use a gaussian profile, add the intensity
                   isum = int * exp(-(t-th2*180/pi)**2 / 2d0 / sigma2)
                   ih = ih + isum

                   ! identify the new peak
                   if (th2 > th2ini .and. th2 < th2end) then
                      if (all(abs(th2p(1:np)-th2) > theps)) then
                         np = np + 1
                         if (np > size(th2p)) then
                            call realloc(th2p,2*np)
                            call realloc(ip,2*np)
                            call realloc(multp,2*np)
                            call realloc(hvecp,3,2*np)
                         end if
                         th2p(np) = th2
                         ip(np) = int
                         multp(np) = 1
                         hvecp(:,np) = (/h,k,l/)
                      else
                         do idx = 1, np
                            if (abs(th2p(idx)-th2) <= theps) exit
                         end do
                         multp(idx) = multp(idx) + 1
                         ! usually the hvec with the most positive indices is the last one
                         hvecp(:,idx) = (/h,k,l/)
                      endif
                   end if
                end if
             end do
          end do
       end do
    end do
    call realloc(th2p,np)
    call realloc(ip,np)
    call realloc(multp,np)
    call realloc(hvecp,3,np)

    ! normalize the intensities to 100
    if (np == 0) &
       call ferror('struct_powder','no peaks found in the 2theta range',faterr)
    ip = ip * multp
    ipmax = maxval(ip)
    ihmax = maxval(ih)
    ih = ih / ihmax * 100
    ip = ip / ipmax * 100

    ! deallocate the multiplicities
    deallocate(multp)

    ! sort the peaks
    allocate(io(np))
    do i = 1, np
       io(i) = i
    end do
    call qcksort(th2p,io,1,np)
    th2p = th2p(io)
    ip = ip(io)
    hvecp = hvecp(:,io)
    deallocate(io)

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
    use environmod, only: environ
    class(crystal), intent(in) :: c
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

    integer :: i, j, k, nat, lvec(3), ierr, iz, jz, kz, npairs, iaux
    integer :: idx, mmult, mza
    real*8 :: int, sigma2, fac, tshift, imax
    logical :: localenv, found
    type(environ) :: le
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

    ! prepare the environment
    localenv = .true.
    if (c%ismolecule .or. c%env%dmax0 > rend+tshift) then
       localenv = .false.
    end if
    if (localenv) then
       call le%build(c%ismolecule,c%nspc,c%spc(1:c%nspc),c%ncel,c%atcel(1:c%ncel),c%m_x2c)
    end if

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
       if (localenv) then
          call le%list_near_atoms(c%at(i)%r,icrd_cart,.true.,nat,ierr,eid,dist,lvec,up2d=rend+tshift,nozero=.true.)
       else
          call c%env%list_near_atoms(c%at(i)%r,icrd_cart,.true.,nat,ierr,eid,dist,lvec,up2d=rend+tshift,nozero=.true.)
       end if

       do j = 1, nat
          ! skip if at a distance lower than rini or higher than rend
          if (dist(j) < rini-tshift .or. dist(j) > rend+tshift) cycle

          ! get info for this atom from the environment
          if (localenv) then
             kz = le%at(eid(j))%is
             jz = le%spc(kz)%z
             idx = le%at(eid(j))%idx
          else
             kz = c%env%at(eid(j))%is
             jz = c%env%spc(kz)%z
             idx = c%env%at(eid(j))%idx
          end if

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
    use environmod, only: environ
    use tools_io, only: ferror, faterr, string
    use param, only: icrd_crys
    class(crystal), intent(in) :: c
    integer, intent(in) :: imax
    real*8, intent(out) :: res(imax)

    integer :: i, j, jini
    integer :: nat, ierr, imax_
    real*8, allocatable :: dist(:)
    real*8 :: dmax0
    logical :: ok
    type(environ) :: le
    logical :: localenv

    integer, parameter :: maxtries = 5 ! maximum number of environment max. distance doubling before crashing

    ! initialize, check input and allocate
    localenv = .false.
    dmax0 = c%env%dmax0
    imax_ = imax
    if (c%ismolecule .and. imax > c%ncel-1) imax_ = c%ncel - 1
    res = 0d0

    ! calculate the amd
    jini = 1
    do i = 1, c%nneq
       ok = .false.
       do j = jini, maxtries
          jini = j
          if (localenv) then
             call le%list_near_atoms(c%at(i)%x,icrd_crys,.true.,nat,ierr,dist=dist,up2n=imax_,nozero=.true.)
          else
             call c%env%list_near_atoms(c%at(i)%x,icrd_crys,.true.,nat,ierr,dist=dist,up2n=imax_,nozero=.true.)
          end if

          ! failed to get the requested number of atoms
          if (ierr /= 0 .or. nat < imax_) then
             ! create an environment with double the maximum distance and point to it
             localenv = .true.
             dmax0 = 2d0 * dmax0
             call le%build(c%ismolecule,c%nspc,c%spc(1:c%nspc),c%ncel,c%atcel(1:c%ncel),c%m_x2c,dmax0=dmax0)
          else
             ok = .true.
             exit
          end if
       end do
       if (.not.ok) &
          call ferror('struct_amd','Error calculating environment for atom: ' // string(i),faterr)
       res = res + c%at(i)%mult * dist
    end do
    res = res / real(c%ncel,8)

  end subroutine amd

  !> Calculate real and reciprocal space sum cutoffs
  module subroutine calculate_ewald_cutoffs(c)
    use tools_io, only: ferror, faterr
    use param, only: pi, rad, sqpi, tpi
    class(crystal), intent(inout) :: c

    real*8, parameter :: sgrow = 1.4d0
    real*8, parameter :: epscut = 1d-5
    real*8, parameter :: eeps = 1d-12

    integer :: i
    real*8 :: aux
    integer :: ia, ib, ic
    real*8 :: alrmax(3), qq
    real*8 :: rcut1, rcut2, err_real
    real*8 :: hcut1, hcut2, err_rec

    if (c%isewald) return

    ! calculate sum of charges and charges**2
    c%qsum = 0d0
    c%q2sum = 0d0
    do i = 1, c%nneq
       qq = c%spc(c%at(i)%is)%qat
       if (abs(qq) < 1d-6) &
          call ferror('calculate_ewald_cutoffs','Some of the charges are 0',faterr)
       c%qsum = c%qsum + real(c%at(i)%mult * qq,8)
       c%q2sum = c%q2sum + real(c%at(i)%mult * qq**2,8)
    end do

    ! determine shortest vector in real space
    aux = 0d0
    do i = 1, 3
       if (c%aa(i) > aux) then
          aux = c%aa(i)
          ia = i
       end if
    end do
    ! determine shortest vector in reciprocal space, dif. from ia
    aux = 0d0
    do i = 1, 3
       if (c%ar(i) > aux .and. i /= ia) then
          aux = c%ar(i)
          ic = i
       end if
    end do
    ! the remaining vector is ib
    ib = 1
    do i = 1, 3
       if (i /= ia .and. i /= ic) ib = i
    end do

    ! convergence parameter
    c%eta = sqrt(c%omega / pi / c%aa(ib) / sin(c%bb(ic)*rad))

    ! real space cutoff
    rcut1 = 1d0
    rcut2 = 2d0 / sgrow
    err_real = 1d30
    do while (err_real >= eeps)
       rcut2 = rcut2 * sgrow
       err_real = pi * c%ncel**2 * c%q2sum / c%omega * c%eta**2 * erfc(rcut2 / c%eta)
    end do
    do while (rcut2-rcut1 >= epscut)
       c%rcut = 0.5*(rcut1+rcut2)
       err_real = pi * c%ncel**2 * c%q2sum / c%omega * c%eta**2 * erfc(c%rcut / c%eta)
       if (err_real > eeps) then
          rcut1 = c%rcut
       else
          rcut2 = c%rcut
       endif
    end do
    c%rcut = 0.5*(rcut1+rcut2)
    ! real space cells to explore
    alrmax = 0d0
    alrmax(1) = c%aa(2) * c%aa(3) * sin(c%bb(1)*rad)
    alrmax(2) = c%aa(1) * c%aa(3) * sin(c%bb(2)*rad)
    alrmax(3) = c%aa(1) * c%aa(2) * sin(c%bb(3)*rad)
    c%lrmax = ceiling(c%rcut * alrmax / c%omega)

    ! reciprocal space cutoff
    hcut1 = 1d0
    hcut2 = 2d0 / sgrow
    err_rec = 1d30
    do while(err_rec >= eeps)
       hcut2 = hcut2 * sgrow
       err_rec = c%ncel**2 * c%q2sum / sqpi / c%eta * erfc(c%eta * hcut2 / 2)
    end do
    do while(hcut2-hcut1 > epscut)
       c%hcut = 0.5*(hcut1+hcut2)
       err_rec = c%ncel**2 * c%q2sum / sqpi / c%eta * erfc(c%eta * c%hcut / 2)
       if (err_rec > eeps) then
          hcut1 = c%hcut
       else
          hcut2 = c%hcut
       endif
    end do
    c%hcut = 0.5*(hcut1+hcut2)
    ! reciprocal space cells to explore
    c%lhmax = ceiling(c%aa(ia) / tpi * c%hcut)
    c%isewald = .true.

  end subroutine calculate_ewald_cutoffs

  !> Calculates the Ewald electrostatic energy, using the input charges.
  module function ewald_energy(c) result(ewe)
    use tools_io, only: uout, string
    class(crystal), intent(inout) :: c
    real*8 :: ewe

    real*8 :: x(3), pot
    integer :: i

    if (.not.c%isewald) &
       call c%calculate_ewald_cutoffs()

    write (uout,'("+ Electrostatic potential at atomic positions")')
    write (uout,'("#id name mult    charge         Vel(Ha/e)")')
    ewe = 0d0
    do i = 1, c%nneq
       x = c%at(i)%x
       pot = c%ewald_pot(x)
       ewe = ewe + c%at(i)%mult * c%spc(c%at(i)%is)%qat * pot
       write (uout,'(99(A," "))') string(i,4), string(c%spc(c%at(i)%is)%name,4),&
          string(c%at(i)%mult,4),&
          string(c%spc(c%at(i)%is)%qat,'e',14,6,3), string(pot,'e',18,10,3)
    end do
    write (uout,*)
    ewe = ewe / 2d0

  end function ewald_energy

  !> Calculate the Ewald electrostatic potential at an arbitrary
  !> position x (crystallographic coords.)  If x is the nucleus j,
  !> return pot - q_j / |r-rj| at rj.
  module function ewald_pot(c,x)
    use param, only: tpi, pi, sqpi, icrd_crys
    class(crystal), intent(inout) :: c
    real*8, intent(in) :: x(3)
    real*8 :: ewald_pot

    real*8 :: rcut2, qnuc
    integer :: i, i1, i2, i3, idnuc
    real*8 :: px(3), lvec(3), d2, d, dh
    real*8 :: sfac_c, sfacp, bbarg
    real*8 :: sum_real, sum_rec, sum0, sum_back

    if (.not.c%isewald) then
       !$omp critical (fill_ewald)
       call c%calculate_ewald_cutoffs()
       !$omp end critical (fill_ewald)
    end if

    ! is this a nuclear position? -> get charge
    idnuc = c%identify_atom(x,icrd_crys)
    if (idnuc > 0) then
       qnuc = c%spc(c%atcel(idnuc)%is)%qat
    else
       qnuc = 0d0
    endif

    ! real space sum
    rcut2 = c%rcut * c%rcut
    sum_real = 0
    do i1 = -c%lrmax(1),c%lrmax(1)
       do i2 = -c%lrmax(2),c%lrmax(2)
          do i3 = -c%lrmax(3),c%lrmax(3)
             lvec = real((/i1,i2,i3/),8)
             do i = 1,c%ncel
                px = x - c%atcel(i)%x - lvec
                d2 = dot_product(px,matmul(c%gtensor,px))
                if (d2 < 1d-12 .or. d2 > rcut2) cycle
                d = sqrt(d2) / c%eta
                sum_real = sum_real + c%spc(c%atcel(i)%is)%qat * erfc(d) / d
             end do
          end do
       end do
    end do
    sum_real = sum_real / c%eta

    ! reciprocal space sum
    sum_rec = 0
    do i1 = -c%lhmax(1),c%lhmax(1)
       do i2 = -c%lhmax(2),c%lhmax(2)
          do i3 = -c%lhmax(3),c%lhmax(3)
             lvec = tpi * (/i1,i2,i3/)
             dh = sqrt(dot_product(lvec,matmul(c%grtensor,lvec)))
             if (dh < 1d-12 .or. dh > c%hcut) cycle
             bbarg = 0.5d0 * dh * c%eta

             sfac_c = 0
             do i = 1, c%ncel
                sfac_c = sfac_c + c%spc(c%atcel(i)%is)%qat * &
                   cos(dot_product(lvec,x-c%atcel(i)%x))
             end do
             sfacp = 2d0 * sfac_c

             sum_rec = sum_rec + sfacp / dh**2 * exp(-bbarg**2)
          end do
       end do
    end do
    sum_rec = sum_rec * 2d0 * pi / c%omega

    ! h = 0 term, applied only at the nucleus
    sum0 = - 2d0 * qnuc / sqpi / c%eta

    ! compensating background charge term
    sum_back = -c%qsum * c%eta**2 * pi / c%omega

    ! sum up and exit
    ewald_pot = sum_real + sum_rec + sum0 + sum_back

  end function ewald_pot

  !> Make a crystal seed (seed) from a crystal structure
  module subroutine makeseed(c,seed,copysym)
    use crystalseedmod, only: crystalseed
    class(crystal), intent(in) :: c
    type(crystalseed), intent(out) :: seed
    logical, intent(in) :: copysym

    integer :: i

    ! general
    seed%isused = .true.
    seed%name = ""
    seed%file = c%file

    ! atoms
    if (copysym .and. c%spgavail) then
       seed%nat = c%nneq
       allocate(seed%x(3,c%nneq),seed%is(c%nneq))
       do i = 1, c%nneq
          seed%x(:,i) = c%at(i)%x
          seed%is(i) = c%at(i)%is
       end do
    else
       seed%nat = c%ncel
       allocate(seed%x(3,c%ncel),seed%is(c%ncel))
       do i = 1, c%ncel
          seed%x(:,i) = c%atcel(i)%x
          seed%is(i) = c%atcel(i)%is
       end do
    end if

    ! species
    seed%nspc = c%nspc
    allocate(seed%spc(c%nspc))
    do i = 1, c%nspc
       seed%spc(i) = c%spc(i)
    end do

    ! cell
    seed%useabr = 2
    seed%aa = 0d0
    seed%bb = 0d0
    seed%m_x2c = c%m_x2c

    ! symmetry
    seed%findsym = -1
    seed%checkrepeats = .false.
    if (copysym .and. c%spgavail) then
       seed%havesym = 1
       seed%neqv = c%neqv
       seed%ncv = c%ncv
       allocate(seed%rotm(3,4,c%neqv),seed%cen(3,c%ncv))
       seed%rotm = c%rotm(:,:,1:c%neqv)
       seed%cen = c%cen(:,1:c%ncv)
    else
       seed%havesym = 0
       seed%neqv = 0
       seed%ncv = 0
    end if

    ! molecular fields
    seed%ismolecule = c%ismolecule
    seed%cubic = (abs(c%aa(1)-c%aa(2)) < 1d-5).and.(abs(c%aa(1)-c%aa(3)) < 1d-5).and.&
       all(abs(c%bb-90d0) < 1d-3)
    seed%border = 0d0
    seed%havex0 = .true.
    seed%molx0 = c%molx0

  end subroutine makeseed

  !> Create a new structure by reordering the atoms in the current
  !> structure. iperm is the permutation vector (atom i in the new
  !> structure is iperm(i) in the old structure).
  module subroutine reorder_atoms(c,iperm,ti)
    use crystalseedmod, only: crystalseed
    class(crystal), intent(inout) :: c
    integer, intent(in) :: iperm(:)
    type(thread_info), intent(in), optional :: ti

    type(crystalseed) :: seed
    real*8, allocatable :: x(:,:)
    integer, allocatable :: is(:)
    integer :: i

    ! make the new seed
    call c%makeseed(seed,.false.)

    ! reorder the atoms
    allocate(x(3,seed%nat),is(seed%nat))
    do i = 1, seed%nat
       x(:,i) = seed%x(:,iperm(i))
       is(i) = seed%is(iperm(i))
    end do
    seed%x = x
    seed%is = is
    deallocate(x,is)

    ! reload the crystal
    call c%struct_new(seed,.true.,ti=ti)

  end subroutine reorder_atoms

  !> Calculate the nearest atom on a grid of points whose dimensions
  !> are given by n. Return the complete atom list IDs in the idg
  !> array.
  module subroutine nearest_atom_grid(c,n,idg)
    use param, only: icrd_crys
    class(crystal), intent(inout) :: c
    integer, intent(in) :: n(3)
    integer, allocatable, intent(inout) :: idg(:,:,:)

    integer :: i, j, k, nid
    real*8 :: xdelta(3,3), dist, x(3)

    ! initialize
    if (allocated(idg)) deallocate(idg)
    allocate(idg(n(1),n(2),n(3)))
    do i = 1, 3
       xdelta(:,i) = 0d0
       xdelta(i,i) = 1d0 / real(n(i),8)
    end do

    !$omp parallel do private(x,nid,dist)
    do k = 1, n(3)
       do j = 1, n(2)
          do i = 1, n(1)
             x = (i-1) * xdelta(:,1) + (j-1) * xdelta(:,2) + (k-1) * xdelta(:,3)
             call c%env%nearest_atom(x,icrd_crys,nid,dist)
             !$omp critical(write)
             idg(i,j,k) = nid
             !$omp end critical(write)
          end do
       end do
    end do
    !$omp end parallel do

  end subroutine nearest_atom_grid

  !> Given a crystal structure (c) and three lattice vectors in cryst.
  !> coords (x0(:,1), x0(:,2), x0(:,3)), build the same crystal
  !> structure using the unit cell given by those vectors. If nnew,
  !> xnew, and isnew are given, replace the nnew atoms with the
  !> new positions (xnew) and species (isnew). If noenv is present
  !> and true, do not load the atomic grids or the environments in
  !> the new cell.
  module subroutine newcell(c,x00,t0,nnew,xnew,isnew,noenv,ti)
    use crystalseedmod, only: crystalseed
    use tools_math, only: det3, matinv, mnorm2
    use tools_io, only: ferror, faterr, warning, string, uout
    use types, only: realloc
    class(crystal), intent(inout) :: c
    real*8, intent(in) :: x00(3,3)
    real*8, intent(in), optional :: t0(3)
    integer, intent(in), optional :: nnew
    real*8, intent(in), optional :: xnew(:,:)
    integer, intent(in), optional :: isnew(:)
    logical, intent(in), optional :: noenv
    type(thread_info), intent(in), optional :: ti

    type(crystalseed) :: ncseed
    logical :: ok, found, atgiven
    real*8 :: x0(3,3), x0inv(3,3), fvol
    real*8 :: x(3), dx(3), dd, t(3)
    integer :: i, j, k, l, m
    integer :: nn
    integer :: nlat
    real*8, allocatable :: xlat(:,:)

    real*8, parameter :: eps = 1d-2

    if (c%ismolecule) &
       call ferror('newcell','NEWCELL incompatible with molecules',faterr)

    ! initialize
    atgiven = (present(nnew).and.present(xnew).and.present(isnew))
    x0 = x00
    dd = det3(x0)

    ! check the new lattice vectors are sane
    if (abs(dd) < eps) then
       call ferror('newcell','invalid input vectors',faterr)
    elseif (dd < 0d0) then
       ! flip the cell
       x0 = -x0
       dd = det3(x0)
       call ferror('newcell','det < 0; vectors flipped',warning)
    endif

    ! set the origin translation
    if (present(t0)) then
       t = t0
    else
       t = 0d0
    end if

    ! species list for the new seed
    allocate(ncseed%spc(c%nspc))
    ncseed%nspc = c%nspc
    ncseed%spc = c%spc

    ! metrics of the new cell
    ncseed%m_x2c = matmul(c%m_x2c,x0)
    ncseed%useabr = 2

    ! check if we have new atomic positions
    if (atgiven) then
       ! sanity check
       if (size(xnew,1) /= 3 .or. size(xnew,2) /= nnew .or. size(isnew) /= nnew) &
          call ferror('newcell','NEWCELL: error in size of new atomic positions array',faterr)

       ! build the atomic info in the new seed from the info in input
       ncseed%nat = nnew
       allocate(ncseed%x(3,nnew),ncseed%is(nnew))
       do i = 1, nnew
          ncseed%x(:,i) = xnew(:,i)
          ncseed%is(i) = isnew(i)
       end do
    else
       ! check that the vectors are pure translations
       do i = 1, 3
          ok = .false.
          do j = 1, c%ncv
             ok = (c%are_lclose(x0(:,i),c%cen(:,j),1d-4))
             if (ok) exit
          end do
          if (.not.ok) then
             write (uout,'("! Error: cell vectors are not lattice translations. This can happen for")')
             write (uout,'("! several reasons but a common one is that the symmetry in your crystal")')
             write (uout,'("! was not correctly calculated, or read from an external file (a cif")')
             write (uout,'("! file) that does not have all the symmetry operations. If this is the")')
             write (uout,'("! case, please try to run SYM RECALC before before attempting the cell")')
             write (uout,'("! transformation.")')
             call ferror("newcell","Cell vector number " // string(i) // &
                " is not a pure translation",faterr)
          end if
       end do

       ! inverse matrix
       x0inv = x0
       call matinv(x0inv,3)

       ! find a star of lattice vectors and supercell centering vectors, if any
       ! first lattice vector is (0 0 0)
       allocate(xlat(3,10))
       xlat = 0d0
       nlat = 1
       do i = minval(floor(x0(1,:))),maxval(ceiling(x0(1,:)))
          do j = minval(floor(x0(2,:))),maxval(ceiling(x0(2,:)))
             do k = minval(floor(x0(3,:))),maxval(ceiling(x0(3,:)))
                x = matmul((/i, j, k/),transpose(x0inv))
                if (any(abs(x-nint(x)) > eps)) then
                   ! this is a new candidate for supercell centering vector
                   ! check if we have it already
                   x = x - floor(x)
                   found = .false.
                   do l = 1, nlat
                      if (all(abs(xlat(:,l) - x) < eps)) then
                         found = .true.
                         exit
                      end if
                   end do
                   if (.not.found) then
                      nlat = nlat + 1
                      if (nlat > size(xlat,2)) call realloc(xlat,3,2*nlat)
                      xlat(:,nlat) = x
                   end if
                endif
             end do
          end do
       end do

       ! build the new atom list
       ncseed%nat = 0
       fvol = abs(det3(x0))
       nn = nint(c%ncel * fvol)
       if (abs(nn - (c%ncel*fvol)) > eps) &
          call ferror('newcell','inconsistent number of atoms in newcell',faterr)
       allocate(ncseed%x(3,nn),ncseed%is(nn))
       do i = 1, nlat
          do j = 1, c%ncel
             ! candidate atom
             x = matmul(c%atcel(j)%x-t,transpose(x0inv)) + xlat(:,i)
             x = x - floor(x)

             ! check if we have it already
             ok = .true.
             do m = 1, ncseed%nat
                dx = x - ncseed%x(:,m)
                dx = abs(dx - nint(dx))
                if (all(dx < eps)) then
                   ok = .false.
                   exit
                end if
             end do
             if (ok) then
                ! add it to the list
                ncseed%nat = ncseed%nat + 1
                if (ncseed%nat > size(ncseed%x,2)) then
                   call realloc(ncseed%x,3,2*ncseed%nat)
                   call realloc(ncseed%is,2*ncseed%nat)
                end if
                ncseed%x(:,ncseed%nat) = x
                ncseed%is(ncseed%nat) = c%atcel(j)%is
             end if
          end do
       end do
       call realloc(ncseed%x,3,ncseed%nat)
       call realloc(ncseed%is,ncseed%nat)
       deallocate(xlat)
    end if

    ! rest of the seed information
    ncseed%isused = .true.
    ncseed%file = c%file
    ncseed%havesym = 0
    ncseed%findsym = -1
    ncseed%ismolecule = .false.

    ! initialize the structure
    call c%struct_new(ncseed,.true.,noenv,ti=ti)

  end subroutine newcell

  !> Transform to the standard cell. If toprim, convert to the
  !> primitive standard cell. If doforce = .true., force the
  !> transformation to the primitive even if it does not lead to a
  !> smaller cell. Return the transformation matrix, or a matrix
  !> of zeros if no change was done.
  module function cell_standard(c,toprim,doforce,refine,noenv,ti) result(x0)
    use iso_c_binding, only: c_double
    use spglib, only: spg_standardize_cell
    use global, only: symprec
    use tools_math, only: det3, matinv
    use tools_io, only: ferror, faterr
    use types, only: realloc
    use param, only: eye
    class(crystal), intent(inout) :: c
    logical, intent(in) :: toprim
    logical, intent(in) :: doforce
    logical, intent(in) :: refine
    logical, intent(in), optional :: noenv
    real*8 :: x0(3,3)
    type(thread_info), intent(in), optional :: ti

    integer :: ntyp, nat
    integer :: i, nnew, iprim, inorefine
    real(c_double), allocatable :: x(:,:)
    integer, allocatable :: types_(:)
    real*8 :: rmat(3,3)

    x0 = 0d0

    ! ignore molecules
    if (c%ismolecule) return

    ! use spglib transformation to the standard cell
    rmat = transpose(c%m_x2c)
    nat = c%ncel
    ntyp = c%nspc
    allocate(x(3,4*c%ncel),types_(4*c%ncel)) ! to make room for spglib
    do i = 1, c%ncel
       x(:,i) = c%atcel(i)%x
       types_(i) = c%atcel(i)%is
    end do

    ! parse options
    iprim = 0
    if (toprim) iprim = 1
    inorefine = 1
    if (refine) inorefine = 0

    nnew = spg_standardize_cell(rmat,x,types_,nat,iprim,inorefine,symprec)
    if (nnew == 0) &
       call ferror("cell_standard","could not find primitive cell",faterr)
    call realloc(x,3,nnew)
    call realloc(types_,nnew)
    rmat = transpose(rmat)
    do i = 1, 3
       rmat(:,i) = c%c2x(rmat(:,i))
    end do

    ! flip the cell?
    if (det3(rmat) < 0d0) rmat = -rmat

    ! rmat = transpose(matinv(c%spg%transformation_matrix))
    if (refine) then
       call c%newcell(rmat,nnew=nnew,xnew=x,isnew=types_,noenv=noenv,ti=ti)
    else
       ! if a primitive is wanted but det is not less than 1, do not make the change
       if (all(abs(rmat - eye) < symprec)) return
       if (toprim .and. .not.(det3(rmat) < 1d0-symprec) .and..not.doforce) return
       call c%newcell(rmat,noenv=noenv,ti=ti)
    end if
    x0 = rmat

  end function cell_standard

  !> Transform to the Niggli cell. Return the transformation matrix,
  !> or a matrix of zeros if no change was done.
  module function cell_niggli(c,noenv,ti) result(x0)
    use spglib, only: spg_niggli_reduce
    use global, only: symprec
    use tools_io, only: ferror, faterr
    use tools_math, only: det3
    use param, only: eye
    class(crystal), intent(inout) :: c
    logical, intent(in), optional :: noenv
    type(thread_info), intent(in), optional :: ti
    real*8 :: x0(3,3)

    real*8 :: rmat(3,3)
    integer :: id, i

    x0 = 0d0

    ! ignore molecules
    if (c%ismolecule) return

    ! use spglib niggli reduction
    rmat = transpose(c%m_x2c)
    id = spg_niggli_reduce(rmat,symprec)
    if (id == 0) &
       call ferror("cell_niggli","could not find Niggli reduction",faterr)
    rmat = transpose(rmat)
    do i = 1, 3
       rmat(:,i) = c%c2x(rmat(:,i))
    end do

    ! flip the cell?
    if (det3(rmat) < 0d0) rmat = -rmat

    ! transform
    if (any(abs(rmat - eye) > symprec)) then
       call c%newcell(rmat,noenv=noenv,ti=ti)
       x0 = rmat
    end if

  end function cell_niggli

  !> Transform to the Delaunay cell. Return the transformation matrix,
  !> or a matrix of zeros if no change was done.
  module function cell_delaunay(c,noenv,ti) result(x0)
    use spglib, only: spg_delaunay_reduce
    use global, only: symprec
    use tools_io, only: ferror, faterr
    use tools_math, only: det3
    use param, only: eye
    class(crystal), intent(inout) :: c
    logical, intent(in), optional :: noenv
    type(thread_info), intent(in), optional :: ti
    real*8 :: x0(3,3)

    real*8 :: rmat(3,3)
    integer :: id, i

    x0 = 0d0

    ! ignore molecules
    if (c%ismolecule) return

    ! use spglib delaunay reduction
    rmat = transpose(c%m_x2c)
    id = spg_delaunay_reduce(rmat,symprec)
    if (id == 0) &
       call ferror("cell_delaunay","could not find Delaunay reduction",faterr)
    rmat = transpose(rmat)
    do i = 1, 3
       rmat(:,i) = c%c2x(rmat(:,i))
    end do

    ! flip the cell?
    if (det3(rmat) < 0d0) rmat = -rmat

    ! transform
    if (any(abs(rmat - eye) > symprec)) then
       call c%newcell(rmat,noenv=noenv,ti=ti)
       x0 = rmat
    end if

  end function cell_delaunay

  !> Write information about the crystal structure to the output. lcrys =
  !> information about the structure. lq = charges.
  module subroutine struct_report(c,lcrys,lq)
    use global, only: iunitname0, dunit0, iunit
    use tools_math, only: gcd
    use tools_io, only: uout, string, ioj_center, ioj_left, ioj_right
    use param, only: bohrtoa, maxzat, pi, atmass, pcamu,bohr2cm
    class(crystal), intent(in) :: c
    logical, intent(in) :: lcrys
    logical, intent(in) :: lq

    integer :: i, j, k, iz, is
    integer :: nelec
    real*8 :: maxdv, xcm(3), x0(3), xlen(3), xang(3), xred(3,3)
    real*8 :: dens, mass
    character(len=:), allocatable :: str1
    integer, allocatable :: nis(:)
    integer :: izp0

    character*1, parameter :: lvecname(3) = (/"a","b","c"/)

    if (lcrys) then
       ! Header
       if (.not.c%ismolecule) then
          write (uout,'("* Crystal structure")')
          write (uout,'("  From: ",A)') string(c%file)
          write (uout,'("  Lattice parameters (bohr): ",3(A,"  "))') &
             string(c%aa(1),'f',decimal=6), string(c%aa(2),'f',decimal=6), string(c%aa(3),'f',decimal=6)
          write (uout,'("  Lattice parameters (ang): ",3(A,"  "))') &
             string(c%aa(1)*bohrtoa,'f',decimal=6), string(c%aa(2)*bohrtoa,'f',decimal=6), string(c%aa(3)*bohrtoa,'f',decimal=6)
          write (uout,'("  Lattice angles (degrees): ",3(A,"  "))') &
             string(c%bb(1),'f',decimal=3), string(c%bb(2),'f',decimal=3), string(c%bb(3),'f',decimal=3)
       else
          write (uout,'("* Molecular structure")')
          write (uout,'("  From: ",A)') string(c%file)
          write (uout,'("  Encompassing cell dimensions (bohr): ",3(A,"  "))') &
             string(c%aa(1),'f',decimal=6), string(c%aa(2),'f',decimal=6), string(c%aa(3),'f',decimal=6)
          write (uout,'("  Encompassing cell dimensions (ang): ",3(A,"  "))') &
             string(c%aa(1)*bohrtoa,'f',decimal=6), string(c%aa(2)*bohrtoa,'f',decimal=6), string(c%aa(3)*bohrtoa,'f',decimal=6)
       endif

       ! Compute unit formula, and z
       allocate(nis(c%nspc))
       nis = 0
       do i = 1, c%nneq
          nis(c%at(i)%is) = nis(c%at(i)%is) + c%at(i)%mult
       end do
       maxdv = gcd(nis,c%nspc)
       write (uout,'("  Empirical formula: ",999(10(A,"(",A,") ")))') &
          (string(c%spc(i)%name), string(nint(nis(i)/maxdv)), i=1,c%nspc)
       deallocate(nis)
       if (.not.c%ismolecule) then
          write (uout,'("  Number of non-equivalent atoms in the unit cell: ",A)') string(c%nneq)
          write (uout,'("  Number of atoms in the unit cell: ",A)') string(c%ncel)
       else
          write (uout,'("  Number of atoms: ",A)') string(c%ncel)
       endif
       write (uout,'("  Number of atomic species: ",A)') string(c%nspc)
       nelec = 0
       mass = 0d0
       do i = 1, c%nneq
          iz = c%spc(c%at(i)%is)%z
          if (iz >= maxzat .or. iz <= 0) cycle
          nelec = nelec + iz * c%at(i)%mult
          mass = mass + atmass(iz) * c%at(i)%mult
       end do
       write (uout,'("  Number of electrons (with zero atomic charge): ",A)') string(nelec)
       if (.not.c%ismolecule) then
          write (uout,'("  Molar mass (amu, per unit cell): ",A)') string(mass,'f',decimal=3)
       else
          write (uout,'("  Molar mass (amu): ",A)') string(mass,'f',decimal=3)
       end if

       ! Cell volume and density, space group short report
       if (.not.c%ismolecule) then
          dens = (mass*pcamu) / (c%omega*bohr2cm**3)
          write (uout,'("  Cell volume (bohr^3): ",A)') string(c%omega,'f',decimal=5)
          write (uout,'("  Cell volume (ang^3): ",A)') string(c%omega * bohrtoa**3,'f',decimal=5)
          write (uout,'("  Density (g/cm^3): ",A)') string(dens,'f',decimal=5)
          ! space group, very short report
          if (c%havesym > 0 .and. c%spgavail) then
             write(uout,'("  Space group (H-M): ",A, " (",A,")")') &
                string(c%spg%international_symbol), string(c%spg%spacegroup_number)
          else
             write(uout,'("  Space group (H-M): ---")')
          end if
       end if

       write (uout,*)
    end if

    if (lq) then
       write (uout,'("+ List of atomic species: ")')
       write (uout,'("# spc = atomic species. Z = atomic number. name = atomic name (symbol).")')
       write (uout,'("# Q = charge. ZPSP = pseudopotential charge.")')
       write (uout,'("# ",99(A," "))') string("spc",3,ioj_center), &
          string("Z",3,ioj_center), string("name",7,ioj_center),&
          string("Q",length=7,justify=ioj_center),&
          string("ZPSP",length=4,justify=ioj_right)
       do i = 1, c%nspc
          str1 = " -- "
          if (c%spc(i)%z > 0) then
             if (c%zpsp(c%spc(i)%z) > 0) &
                str1 = string(c%zpsp(c%spc(i)%z))
          end if
          write (uout,'("  ",99(A," "))') string(i,3,ioj_center), &
             string(c%spc(i)%z,3,ioj_center), string(c%spc(i)%name,7,ioj_center),&
             string(c%spc(i)%qat,'f',length=7,decimal=4,justify=ioj_right),&
             str1
       end do
       write (uout,*)
    end if

    if (lcrys) then
       ! List of atoms in crystallographic coordinates
       if (.not.c%ismolecule) then
          write (uout,'("+ List of non-equivalent atoms in the unit cell (cryst. coords.): ")')
          write (uout,'("# at = complete list atomic ID. xyz = Cartesian coordinates. spc = atomic species.")')
          write (uout,'("# wyck = wyckoff position. name = atomic name (symbol). mult = multiplicity.")')
          write (uout,'("# Z = atomic number.")')

          write (uout,'("# ",99(A," "))') string("nat",3,ioj_center), &
             string("x",14,ioj_center), string("y",14,ioj_center),&
             string("z",14,ioj_center), string("spc",3,ioj_center), string("wyck",4,ioj_center), &
             string("name",7,ioj_center), string("mult",4,ioj_center), string("Z",3,ioj_center)
          do i=1, c%nneq
             is = c%at(i)%is
             write (uout,'("  ",99(A," "))') string(i,3,ioj_center),&
                (string(c%at(i)%x(j),'f',length=14,decimal=10,justify=3),j=1,3),&
                string(is,3,ioj_center), string(c%at(i)%mult,3,ioj_right) // c%at(i)%wyc, &
                string(c%spc(is)%name,7,ioj_center), &
                string(c%at(i)%mult,4,ioj_center), string(c%spc(is)%z,3,ioj_center)
          enddo
          write (uout,*)

          write (uout,'("+ List of atoms in the unit cell (cryst. coords.): ")')
          write (uout,'("# at = complete list atomic ID. xyz = Cartesian coordinates. spc = atomic species.")')
          write (uout,'("# name = atomic name (symbol). Z = atomic number. nat = non-equivalent atom id.")')
          if (allocated(c%mol) .and. c%nmol > 0) then
             write (uout,'("# mol = molecular fragment.")')
             str1 = string("mol",3,ioj_center)
          else
             str1 = ""
          end if
          write (uout,'("# ",99(A," "))') string("at",3,ioj_center),&
             string("x",14,ioj_center), string("y",14,ioj_center),&
             string("z",14,ioj_center), string("spc",3,ioj_center), string("name",7,ioj_center),&
             string("Z",3,ioj_center), string("nat",3,ioj_center), str1
          do i=1,c%ncel
             is = c%atcel(i)%is
             if (allocated(c%mol) .and. c%nmol > 0) then
                str1 = string(c%idatcelmol(i),3,ioj_center)
             else
                str1 = ""
             end if
             write (uout,'("  ",99(A," "))') &
                string(i,3,ioj_center),&
                string(c%atcel(i)%x(1),'f',length=14,decimal=10,justify=3),&
                string(c%atcel(i)%x(2),'f',length=14,decimal=10,justify=3),&
                string(c%atcel(i)%x(3),'f',length=14,decimal=10,justify=3),&
                string(is,3,ioj_center),&
                string(c%spc(is)%name,7,ioj_center),&
                string(c%spc(is)%z,3,ioj_center),&
                string(c%atcel(i)%idx,3,ioj_center), str1
          enddo
          write (uout,*)

          write (uout,'("+ Lattice vectors (",A,")")') iunitname0(iunit)
          do i = 1, 3
             write (uout,'("    ",A,": ",3(A," "))') lvecname(i),&
                (string(c%m_x2c(j,i)*dunit0(iunit),'f',length=16,decimal=10,justify=5),j=1,3)
          end do
          write (uout,*)
       end if

       ! List of atoms in Cartesian coordinates
       write (uout,'("+ List of atoms in Cartesian coordinates (",A,"): ")') iunitname0(iunit)
       write (uout,'("# at = complete list atomic ID. xyz = Cartesian coordinates. spc = atomic species.")')
       write (uout,'("# name = atomic name (symbol). Z = atomic number. dnn = nearest-neighbor distance.")')
       if (allocated(c%mol) .and. c%nmol > 0) then
          write (uout,'("# nat = non-equivalent atom id. mol = molecular fragment.")')
          str1 = string("mol",3,ioj_center)
       else
          write (uout,'("# nat = non-equivalent atom id. ")')
          str1 = ""
       end if
       write (uout,'("# ",99(A," "))') string("at",3,ioj_center), &
          string("x",16,ioj_center), string("y",16,ioj_center),&
          string("z",16,ioj_center), string("spc",3,ioj_center), string("name",7,ioj_center),&
          string("Z",3,ioj_center), string("dnn",10,ioj_center), string("nat",3,ioj_center), str1
       do i=1,c%ncel
          is = c%atcel(i)%is
          if (allocated(c%mol) .and. c%nmol > 0) then
             str1 = string(c%idatcelmol(i),3,ioj_center)
          else
             str1 = ""
          end if
          write (uout,'("  ",99(A," "))') &
             string(i,3,ioj_center),&
             (string((c%atcel(i)%r(j)+c%molx0(j))*dunit0(iunit),'f',length=16,decimal=10,justify=5),j=1,3),&
             string(is,3,ioj_center),string(c%spc(is)%name,7,ioj_center), string(c%spc(is)%z,3,ioj_center),&
             string(2d0*c%at(c%atcel(i)%idx)%rnn2*dunit0(iunit),'f',length=10,decimal=4,justify=4), &
             string(c%atcel(i)%idx,3,ioj_center), str1
       enddo
       write (uout,*)

       ! Encompassing region for the molecule
       if (c%ismolecule) then
          write (uout,'("+ Limits of the molecular cell (in fractions of the unit cell).")')
          write (uout,'("# The part of the unit cell outside the molecular cell represents")')
          write (uout,'("# infinity (no CPs or gradient paths in it).")')
          write (uout,'("  x-axis: ",A," -> ",A)') trim(string(c%molborder(1),'f',10,4)), trim(string(1d0-c%molborder(1),'f',10,4))
          write (uout,'("  y-axis: ",A," -> ",A)') trim(string(c%molborder(2),'f',10,4)), trim(string(1d0-c%molborder(2),'f',10,4))
          write (uout,'("  z-axis: ",A," -> ",A)') trim(string(c%molborder(3),'f',10,4)), trim(string(1d0-c%molborder(3),'f',10,4))
          write (uout,*)
       end if

       ! Write symmetry operations
       if (.not.c%ismolecule) then
          write(uout,'("+ List of symmetry operations (",A,"):")') string(c%neqv)
          do k = 1, c%neqv
             write (uout,'("  Operation ",A,":")') string(k)
             write (uout,'(2("    ",4(A," ")/),"    ",4(A," "))') &
                ((string(c%rotm(i,j,k),'f',length=9,decimal=6,justify=3), j = 1, 4), i = 1, 3)
          enddo
          write (uout,*)

          call c%struct_report_symxyz(doaxes=.true.)

          write(uout,'("+ List of centering vectors (",A,"):")') string(c%ncv)
          do k = 1, c%ncv
             write (uout,'("  Vector ",A,": ",3(A," "))') string(k), &
                (string(c%cen(i,k),'f',length=9,decimal=6), i = 1, 3)
          enddo
          write (uout,*)

          call c%struct_report_symmetry()

          write (uout,'("+ Cartesian/crystallographic coordinate transformation matrices:")')
          write (uout,'("  A = car to crys (xcrys = A * xcar, ",A,"^-1)")') iunitname0(iunit)
          do i = 1, 3
             write (uout,'("    ",3(A," "))') (string(c%m_c2x(i,j)/dunit0(iunit),'f',length=16,decimal=10,justify=5),j=1,3)
          end do
          write (uout,'("  B = crys to car (xcar = B * xcrys, ",A,")")') iunitname0(iunit)
          do i = 1, 3
             write (uout,'("    ",3(A," "))') (string(c%m_x2c(i,j)*dunit0(iunit),'f',length=16,decimal=10,justify=5),j=1,3)
          end do
          write (uout,'("  G = metric tensor (B''*B, ",A,"^2)")') iunitname0(iunit)
          do i = 1, 3
             write (uout,'("    ",3(A," "))') (string(c%gtensor(i,j)*dunit0(iunit)**2,'f',length=16,decimal=10,justify=5),j=1,3)
          end do
          write (uout,*)
       end if

       ! Discrete molecules, if available
       if (allocated(c%nstar) .and. allocated(c%mol) .and. c%nmol > 0) then
          write (uout,'("+ List of fragments in the system (",A,")")') string(c%nmol)
          write (uout,'("# Id = fragment ID. nat = number of atoms in fragment. C-o-m = center of mass (",A,").")')&
             iunitname0(iunit)
          write (uout,'("# Discrete = is this fragment finite?")')

          if (.not.c%ismolecule.and.c%ismol3d.and.allocated(c%idxmol)) then
             ! the table with the molecular equivalence integers
             write (uout,'("# idx = 0 (molecule is symmetry-unique), n>0 (equivalent to molecule n), -1 (symmetric to self)")')
             write (uout,'("# Id  nat           Center of mass            Discrete idx")')
             do i = 1, c%nmol
                if (c%ismolecule) then
                   xcm = (c%mol(i)%cmass()+c%molx0) * dunit0(iunit)
                else
                   xcm = c%c2x(c%mol(i)%cmass())
                end if
                write (uout,'(99("  ",A))') string(i,3,ioj_left), string(c%mol(i)%nat,4,ioj_left),&
                   (string(xcm(j),'f',10,6,3),j=1,3), string(c%mol(i)%discrete), string(c%idxmol(i),3,ioj_right)
             end do
          else
             ! the simple table with just the list of fragments
             write (uout,'("# Id  nat           Center of mass            Discrete")')
             do i = 1, c%nmol
                if (c%ismolecule) then
                   xcm = (c%mol(i)%cmass()+c%molx0) * dunit0(iunit)
                else
                   xcm = c%c2x(c%mol(i)%cmass())
                end if
                write (uout,'(99("  ",A))') string(i,3,ioj_left), string(c%mol(i)%nat,4,ioj_left),&
                   (string(xcm(j),'f',10,6,3),j=1,3), string(c%mol(i)%discrete)
             end do
          end if

          if (.not.c%ismolecule) then
             if (c%ismol3d .or. c%nlvac == 3) then
                write (uout,'(/"+ This is a molecular crystal.")')
                write (uout,'("  Number of molecules per cell (Z) = ",A)') string(c%nmol)
                izp0 = 0
                do i = 1, c%nmol
                   if (c%idxmol(i) < 0) then
                      izp0 = -1
                      exit
                   elseif (c%idxmol(i) == 0) then
                      izp0 = izp0 + 1
                   end if
                end do
                if (izp0 > 0) then
                   write (uout,'("  Number of molecules in the asymmetric unit (Z'') = ",A)') string(izp0)
                else
                   write (uout,'("  Number of molecules in the asymmetric unit (Z'') < 1")')
                end if

             else if (c%nlvac == 2) then
                write (uout,'(/"+ This is a 1D periodic (polymer) structure.")')
                write (uout,'("  Vacuum lattice vectors: (",2(A," "),A,"), (",2(A," "),A,")")') &
                   (string(c%lvac(j,1)),j=1,3), (string(c%lvac(j,2)),j=1,3)
                write (uout,'("  Connected lattice vectors: (",2(A," "),A,")")') &
                   (string(c%lcon(j,1)),j=1,3)
             else if (c%nlvac == 1) then
                write (uout,'(/"+ This is a 2D periodic (layered) structure.")')
                write (uout,'("  Vacuum lattice vectors: (",2(A," "),A,")")') &
                   (string(c%lvac(j,1)),j=1,3)
                write (uout,'("  Connected lattice vectors: (",2(A," "),A,"), (",2(A," "),A,")")') &
                   (string(c%lcon(j,1)),j=1,3), (string(c%lcon(j,2)),j=1,3)
             else
                write (uout,'(/"+ This is a 3D periodic structure.")')
             end if
          end if
          write (uout,*)
       end if

       ! Number of atoms in the atomic environment
       call c%env%report()

       ! Wigner-Seitz cell
       if (.not.c%ismolecule) then
          write (uout,'("+ Vertex of the WS cell in cryst. coords. (",A,")")') string(c%ws_nv)
          write (uout,'("# id = vertex ID. xyz = vertex cryst. coords. d = vertex distance to origin (",A,").")') iunitname0(iunit)
          write (uout,'(5("  ",A))') string("id",length=3,justify=ioj_right),&
             string("x",length=11,justify=ioj_center),&
             string("y",length=11,justify=ioj_center),&
             string("z",length=11,justify=ioj_center),&
             string("d ("//iunitname0(iunit)//")",length=14,justify=ioj_center)
          do i = 1, c%ws_nv
             x0 = c%x2c(c%ws_x(:,i))
             write (uout,'(5("  ",A))') string(i,length=3,justify=ioj_right), &
                (string(c%ws_x(j,i),'f',length=11,decimal=6,justify=4),j=1,3), &
                string(norm2(x0)*dunit0(iunit),'f',length=14,decimal=8,justify=4)
          enddo
          write (uout,*)

          write (uout,'("+ Faces of the WS cell (",A,")")') string(c%ws_nf)
          write (uout,'("# Face ID: vertexID1 vertexID2 ...")')
          do i = 1, c%ws_nf
             write (uout,'("  ",A,": ",999(A," "))') string(i,length=2,justify=ioj_right), &
                (string(c%ws_iside(j,i),length=2),j=1,c%ws_nside(i))
          end do
          write (uout,*)

          write (uout,'("+ Lattice vectors for the Wigner-Seitz neighbors")')
          write (uout,'("# FaceID: Voronoi lattice vector (cryst. coords.)")')
          do i = 1, c%ws_nf
             write (uout,'("  ",A,": ",99(A," "))') string(i,length=2,justify=ioj_right), &
                (string(c%ws_ineighx(j,i),length=2,justify=ioj_right),j=1,3)
          end do
          write (uout,*)

          write (uout,'("+ Lattice vectors for the Delaunay reduced cell (cryst. coords.)")')
          do i = 1, 3
             write (uout,'("  ",A,": ",99(A," "))') lvecname(i), &
                (string(nint(c%m_xr2x(j,i)),length=2,justify=ioj_right),j=1,3)
          end do

          do i = 1, 3
             x0 = c%m_xr2x(:,i)
             xred(:,i) = c%x2c(x0)
             xlen(i) = norm2(xred(:,i))
          end do
          xang(1) = acos(dot_product(xred(:,2),xred(:,3)) / xlen(2) / xlen(3)) * 180d0 / pi
          xang(2) = acos(dot_product(xred(:,1),xred(:,3)) / xlen(1) / xlen(3)) * 180d0 / pi
          xang(3) = acos(dot_product(xred(:,1),xred(:,2)) / xlen(1) / xlen(2)) * 180d0 / pi

          write (uout,'("  Delaunay reduced cell lengths: ",99(A," "))') &
             (string(xlen(j),'f',decimal=6,justify=ioj_right),j=1,3)
          write (uout,'("  Delaunay reduced cell angles: ",99(A," "))') &
             (string(xang(j),'f',decimal=3,justify=ioj_right),j=1,3)
          write (uout,*)

          write (uout,'("+ Is the cell orthogonal? ",L1)') c%isortho
          write (uout,'("+ Is the reduced cell orthogonal? ",L1/)') c%isortho_del
       end if
    end if

  end subroutine struct_report

  !> Write information about the crystal symmetry.
  module subroutine struct_report_symmetry(c)
    use tools_io, only: string, uout
    class(crystal), intent(in) :: c

    character(len=3) :: schpg
    integer :: holo, laue

    if (.not.c%ismolecule .and. c%havesym > 0) then
       write(uout,'("+ Crystal symmetry information")')
       if (c%spgavail) then
          write(uout,'("  Space group (Hermann-Mauguin): ",A, " (number ",A,")")') &
             string(c%spg%international_symbol), string(c%spg%spacegroup_number)
          write(uout,'("  Point group (Hermann-Mauguin): ",A)') string(c%spg%pointgroup_symbol)

          call pointgroup_info(c%spg%pointgroup_symbol,schpg,holo,laue)
          write(uout,'("  Point group (Schoenflies): ",A)') string(schpg)
          write(uout,'("  Holohedry: ",A)') string(holo_string(holo))
          write(uout,'("  Laue class: ",A)') string(laue_string(laue))
       else
          write(uout,'("  -- not available --")')
       end if
       write (uout,*)
    end if

  end subroutine struct_report_symmetry

  !> Write the list of symmetry operations to stdout, using
  !> crystallographic notation (if possible). If strfin is present,
  !> return the strings in that variable instead of writing them to
  !> uout.
  module subroutine struct_report_symxyz(c,strfin,doaxes)
    use tools_math, only: eig, det3
    use tools_io, only: uout, string, ioj_right
    use global, only: symprec
    use param, only: mlen, pi
    class(crystal), intent(in) :: c
    character(len=mlen), intent(out), optional :: strfin(c%neqv*c%ncv)
    logical, intent(in), optional :: doaxes

    real*8, parameter :: rfrac(25) = (/-12d0/12d0,-11d0/12d0,-10d0/12d0,&
       -9d0/12d0,-8d0/12d0,-7d0/12d0,-6d0/12d0,-5d0/12d0,-4d0/12d0,-3d0/12d0,&
       -2d0/12d0,-1d0/12d0,0d0/12d0,1d0/12d0,2d0/12d0,3d0/12d0,4d0/12d0,&
       5d0/12d0,6d0/12d0,7d0/12d0,8d0/12d0,9d0/12d0,10d0/12d0,11d0/12d0,12d0/12d0/)
    character*6, parameter :: sfrac(25) = (/"      ","-11/12","-5/6  ",&
       "-3/4  ","-2/3  ","-7/12 ","-1/2  ","-5/12 ","-1/3  ","-1/4  ","-1/6  ",&
       "-1/12 ","      ","1/12  ","1/6   ","1/4   ","1/3   ","5/12  ","1/2   ",&
       "7/12  ","2/3   ","3/4   ","5/6   ","11/12 ","      "/)
    character*1, parameter :: xyz(3) = (/"x","y","z"/)
    real*8, parameter :: eps = 1d-5

    logical :: ok, iszero, doax
    integer :: i1, i2, i, j, k, idx, rotnum
    character(len=mlen) :: strout(c%neqv*c%ncv)
    real*8 :: xtrans, rmat(3,3), eval(3), evali(3), rotaxis(3)
    real*8 :: trace, det, ang, ridx
    character(len=mlen), allocatable :: rotchar(:)

    doax = .false.
    if (present(doaxes)) doax = doaxes

    i = 0
    do i1 = 1, c%ncv
       do i2 = 1, c%neqv
          i = i + 1
          strout(i) = "<not found>"
       end do
    end do

    ! classify the rotations
    if (doax) then
       allocate(rotchar(c%neqv))
       do i2 = 1, c%neqv
          rmat = c%rotm(:,1:3,i2)
          trace = rmat(1,1)+rmat(2,2)+rmat(3,3)
          det = det3(rmat)

          if (abs(trace - 3d0) < 1d-5) then
             rotchar(i2) = " 1"
          elseif (abs(trace + 3d0) < 1d-5) then
             rotchar(i2) = "-1"
          else
             ! determine the angle of rotation
             ang = 0.5d0*(trace-det)
             if (abs(ang) > 1d0) ang = sign(1d0,ang)
             ang = acos(ang)
             if (abs(ang) < eps) then
                rotnum = 1
             else
                rotnum = nint((2d0*pi) / ang)
             end if
             if (det > 0d0) then
                rotchar(i2) = string(rotnum,2,ioj_right)
             else
                if (rotnum == 2) then
                   rotchar(i2) = " m"
                else
                   rotchar(i2) = string(-rotnum,2,ioj_right)
                end if
             end if

             ! determine the axis of rotation
             call eig(rmat,3,eval,evali)
             idx = 0
             do j = 1, 3
                if (abs(evali(j)) < eps .and. abs(eval(j)-det) < eps) then
                   idx = j
                   exit
                end if
             end do
             if (idx > 0) then
                rotaxis = rmat(:,idx)
                ! divide by the smallest non-zero element in the vector
                ridx = 1d40
                do j = 1, 3
                   if (abs(rotaxis(j)) > eps .and. abs(rotaxis(j)) < ridx) ridx = abs(rotaxis(j))
                end do
                rotaxis = rotaxis / ridx
             endif

             if (all(abs(rotaxis - nint(rotaxis)) < eps)) then
                rotchar(i2) = rotchar(i2)(1:2) // " [" // string(nint(rotaxis(1))) // "," //&
                   string(nint(rotaxis(2))) // "," // string(nint(rotaxis(3))) // "]"
             else
                rotchar(i2) = rotchar(i2)(1:2) // " [" // string(rotaxis(1),'f',decimal=1) // "," //&
                   string(rotaxis(2),'f',decimal=1) // "," // string(rotaxis(3),'f',decimal=1) // "]"
             end if

             rotaxis = c%x2c(rotaxis)
             rotaxis = rotaxis / norm2(rotaxis)
             rotchar(i2) = trim(rotchar(i2)) // "; [" // string(rotaxis(1),'f',decimal=3) // "," //&
                string(rotaxis(2),'f',decimal=3) // "," // string(rotaxis(3),'f',decimal=3) // "]"
          end if
       end do
    end if

    i = 0
    main: do i1 = 1, c%ncv
       loopi: do i2 = 1, c%neqv
          i = i + 1
          strout(i) = ""

          do j = 1, 3
             ! translation
             ok = .false.
             do k = 1, 25
                xtrans = c%rotm(j,4,i2)+c%cen(j,i1) - rfrac(k)
                xtrans = xtrans - nint(xtrans)
                if (abs(xtrans) < symprec) then
                   ok = .true.
                   strout(i) = trim(strout(i)) // sfrac(k)
                   iszero = (k == 13) .or. (k == 1) .or. (k == 25)
                   exit
                end if
             end do
             if (.not.ok) then
                strout(i) = "<not found>"
                cycle loopi
             end if

             ! rotation
             do k = 1, 3
                if (abs(c%rotm(j,k,i2) - 1d0) < symprec) then
                   if (iszero) then
                      strout(i) = trim(strout(i)) // xyz(k)
                   else
                      strout(i) = trim(strout(i)) // "+" // xyz(k)
                   end if
                   iszero = .false.
                elseif (abs(c%rotm(j,k,i2) + 1d0) < symprec) then
                   strout(i) = trim(strout(i)) // "-" // xyz(k)
                   iszero = .false.
                elseif (abs(c%rotm(j,k,i2)) > symprec) then
                   strout(i) = "<not found>"
                   cycle loopi
                end if
             end do

             ! the comma
             if (j < 3) &
                strout(i) = trim(strout(i)) // ","
          end do
          if (doax) then
             strout(i) = string(strout(i),30) // " ## " // trim(rotchar(i2))
          end if
       end do loopi
    end do main
    if (doax) deallocate(rotchar)

    if (present(strfin)) then
       strfin = strout
    else
       write(uout,'("+ List of symmetry operations in crystallographic notation:")')
       if (doax) &
          write(uout,'("# number: operation ... ## order of rotation [axis]; [axis in Cartesian]")')
       do k = 1, c%neqv*c%ncv
          write (uout,'("   ",A,": ",A)') string(k,length=3,justify=ioj_right), string(strout(k))
       enddo
       write (uout,*)
    end if

  end subroutine struct_report_symxyz

  !> Write the field info to a JSON object. The structure object to
  !> json with root p.
  module subroutine struct_write_json(c,json,p)
    use json_module, only: json_value, json_core
    class(crystal), intent(in) :: c
    type(json_core), intent(inout) :: json
    type(json_value), pointer, intent(inout) :: p

    type(json_value), pointer :: s, ap, arr

    integer :: i
    character(len=mlen), allocatable :: strfin(:)

    if (.not.c%isinit) return
    call json%create_object(s,'structure')
    call json%add(p,s)

    call json%add(s,'cell_lengths',c%aa)
    call json%add(s,'cell_angles',c%bb)
    call json%add(s,'cell_volume',c%omega)
    call json%add(s,'crys_to_cart_matrix',reshape(c%m_x2c,(/9/)))
    call json%add(s,'cart_to_crys_matrix',reshape(c%m_c2x,(/9/)))

    call json%add(s,'is_molecule',c%ismolecule)
    call json%add(s,'molecule_centering_vector',c%molx0)
    call json%add(s,'molecular_cell_border',c%molborder)
    call json%add(s,'periodicity',3-c%nlvac)

    call json%add(s,'number_of_species',c%nspc)
    call json%create_array(arr,'species')
    call json%add(s,arr)
    do i = 1, c%nspc
       call json%create_object(ap,'')
       call json%add(ap,'id',i)
       call json%add(ap,'name',trim(c%spc(i)%name))
       call json%add(ap,'atomic_number',c%spc(i)%z)
       call json%add(arr,ap)
       nullify(ap)
    end do
    nullify(arr)

    call json%add(s,'number_of_nonequivalent_atoms',c%nneq)
    call json%create_array(arr,'nonequivalent_atoms')
    call json%add(s,arr)
    do i = 1, c%nneq
       call json%create_object(ap,'')
       call json%add(ap,'id',i)
       call json%add(ap,'species',c%at(i)%is)
       call json%add(ap,'fractional_coordinates',c%at(i)%x(:))
       call json%add(ap,'cartesian_coordinates',c%at(i)%r(:))
       call json%add(ap,'multiplicity',c%at(i)%mult)
       call json%add(ap,'wyckoff_letter',c%at(i)%wyc)
       call json%add(ap,'half_nn_distance',c%at(i)%rnn2)
       call json%add(arr,ap)
       nullify(ap)
    end do
    nullify(arr)

    call json%add(s,'number_of_cell_atoms',c%ncel)
    call json%create_array(arr,'cell_atoms')
    call json%add(s,arr)
    do i = 1, c%ncel
       call json%create_object(ap,'')
       call json%add(ap,'id',i)
       call json%add(ap,'species',c%atcel(i)%is)
       call json%add(ap,'fractional_coordinates',c%atcel(i)%x(:))
       call json%add(ap,'cartesian_coordinates',c%atcel(i)%r(:))
       call json%add(ap,'nonequivalent_id',c%atcel(i)%idx)
       call json%add(ap,'symop_to_nneq',c%atcel(i)%ir)
       call json%add(ap,'centering_vector_to_nneq',c%atcel(i)%ic)
       call json%add(ap,'lattice_vector_to_nneq',c%atcel(i)%lvec)
       call json%add(arr,ap)
       nullify(ap)
    end do
    nullify(arr)

    if (.not.c%ismolecule) then
       call json%add(s,'have_symmetry',c%havesym > 0)
       if (c%spgavail) then
          call json%add(s,'space_group_hm',trim(c%spg%international_symbol))
          call json%add(s,'space_group_ita_number',c%spg%spacegroup_number)
       end if
       call json%add(s,'number_of_symops',c%neqv)

       allocate(strfin(c%neqv*c%ncv))
       call c%struct_report_symxyz(strfin)
       call json%create_array(arr,'symops')
       call json%add(s,arr)
       do i = 1, c%neqv
          call json%create_object(ap,'')
          call json%add(ap,'id',i)
          call json%add(ap,'operation',trim(strfin(i)))
          call json%add(ap,'rotation',reshape(c%rotm(1:3,1:3,i),(/9/)))
          call json%add(ap,'translation',c%rotm(:,4,i))
          call json%add(arr,ap)
          nullify(ap)
       end do
       nullify(arr)

       call json%add(s,'number_of_centering_vectors',c%ncv)
       call json%add(s,'centering_vectors',reshape(c%cen(1:3,1:c%ncv),(/3*c%ncv/)))
    endif
    call json%add(s,'number_of_molecular_fragments',c%nmol)
    nullify(s)

  end subroutine struct_write_json

  !> Use the spg library to find information about the space group,
  !> encapsulated in the spg user-defined type (spg).
  !> Input: cell vectors (m_x2c), ncel, atcel(:), at(:).
  !> If usenneq, use nneq and at(:) instead of ncel and atcel(:).
  !> If error, return the error in errmsg. Otherwise, return
  !> a zero-length string.
  module subroutine spglib_wrap(c,spg,usenneq,errmsg,ti)
    use iso_c_binding, only: c_double
    use spglib, only: spg_get_dataset, spg_get_error_message
    use global, only: symprec
    use tools_io, only: string, equal
    use param, only: maxzat0
    use types, only: realloc
    class(crystal), intent(in) :: c
    type(SpglibDataset), intent(inout) :: spg
    logical, intent(in) :: usenneq
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    real(c_double) :: lattice(3,3)
    real(c_double), allocatable :: x(:,:)
    integer, allocatable :: typ(:)
    integer :: ntyp, nat
    integer :: i, iz(maxzat0)
    character(len=32) :: error

    ! get the dataset from spglib
    errmsg = ""
    lattice = transpose(c%m_x2c)
    iz = 0
    ntyp = 0
    if (usenneq) then
       nat = c%nneq
       ntyp = c%nspc
       allocate(x(3,c%nneq),typ(c%nneq))
       do i = 1, c%nneq
          x(:,i) = c%at(i)%x
          typ(i) = c%at(i)%is
       end do
    else
       nat = c%ncel
       ntyp = c%nspc
       allocate(x(3,c%ncel),typ(c%ncel))
       do i = 1, c%ncel
          x(:,i) = c%atcel(i)%x
          typ(i) = c%atcel(i)%is
       end do
    end if

    spg = spg_get_dataset(lattice,x,typ,nat,symprec)
    deallocate(x,typ)

    ! check error messages (only if this is not threaded, because
    ! this routine is not thread-safe)
    if (.not.present(ti)) then
       error = trim(spg_get_error_message(c%spg%spglib_error))
       if (.not.equal(error,"no error")) then
          errmsg = string(error)
          return
       end if
    end if

  end subroutine spglib_wrap

  !> Set the Wyckoff positions in crystal c from the information in spg.
  !> Assumes the atoms are in the same order. If spg is not present
  !> or not initialized, then fill with "?".
  module subroutine spgtowyc(c,spg)
    class(crystal), intent(inout) :: c
    type(SpglibDataset), intent(inout), optional :: spg

    integer :: i, idx, ilet
    logical :: doit

    character(len=26), parameter :: wycklet = "abcdefghijklmnopqrstuvwxyz"

    doit = .false.
    if (present(spg)) then
       if (spg%n_atoms > 0) doit = .true.
    end if

    ! fill wyckoff letters
    if (doit) then
       do i = 1, c%ncel
          idx = c%atcel(i)%idx
          if (c%at(idx)%wyc == "?") then
             ilet = c%spg%wyckoffs(i)+1
             c%at(idx)%wyc = wycklet(ilet:ilet)
          end if
       end do
    else
       do i = 1, c%nneq
          c%at(i)%wyc = "?"
       end do
    end if

  end subroutine spgtowyc

  !> Calculate the crystal symmetry operations.
  !> Input: cell vectors (m_x2c), ncel, atcel(:), at(:)
  !> Output: neqv, rotm, spg, at()%mult, ncel, and atcel().
  !> If usenneq, use nneq and at(:) instead of ncel and atcel(:).
  !> If error, errmsg in output has length > 0 and contains the error
  !> message.
  module subroutine calcsym(c,usenneq,errmsg,ti)
    use iso_c_binding, only: c_double
    use spglib, only: spg_get_error_message
    use global, only: symprec
    use tools_io, only: string, equal
    use param, only: eyet, eye
    use types, only: realloc
    class(crystal), intent(inout) :: c
    logical, intent(in) :: usenneq
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    integer, allocatable :: iidx(:)
    integer :: i, j, k, iat, idx
    logical :: found
    real*8 :: rotm(3,3), x0(3)

    c%spgavail = .false.
    call c%spglib_wrap(c%spg,usenneq,errmsg,ti=ti)
    if (len_trim(errmsg) > 0) return
    c%spgavail = .true.

    ! make a copy of nneq into ncel, if appropriate
    if (usenneq) then
       c%ncel = c%nneq
       if (allocated(c%atcel)) deallocate(c%atcel)
       allocate(c%atcel(c%ncel))
       do i = 1, c%ncel
          c%atcel(i)%x = c%at(i)%x
          c%atcel(i)%r = c%at(i)%r
          c%atcel(i)%is = c%at(i)%is
       end do
    end if

    ! re-write nneq list based on the information from spg
    allocate(iidx(c%ncel))
    iidx = 0
    c%nneq = 0
    do i = 1, c%spg%n_atoms
       idx = c%spg%equivalent_atoms(i) + 1
       if (iidx(idx) == 0) then
          c%nneq = c%nneq + 1
          iidx(idx) = c%nneq
          if (c%nneq > size(c%at,1)) call realloc(c%at,2*c%nneq)
          c%at(c%nneq)%x = c%atcel(idx)%x
          c%at(c%nneq)%r = c%atcel(idx)%r
          c%at(c%nneq)%is = c%atcel(idx)%is
          c%at(c%nneq)%mult = 1
          c%at(c%nneq)%rnn2 = 0d0
          c%at(c%nneq)%wyc = "?"
       else
          c%at(iidx(idx))%mult = c%at(iidx(idx))%mult + 1
       end if
       c%atcel(i)%idx = iidx(idx)
    end do
    call realloc(c%at,c%nneq)

    ! unpack spglib's output into pure translations and symops
    c%neqv = 1
    c%rotm(:,:,1) = eyet
    c%ncv = 1
    if (.not.allocated(c%cen)) allocate(c%cen(3,1))
    c%cen = 0d0
    do i = 1, c%spg%n_operations
       rotm = transpose(c%spg%rotations(:,:,i))

       ! is this a pure translation?
       if (all(abs(rotm - eye) < symprec)) then
          found = .false.
          do j = 1, c%ncv
             if (all(abs(c%spg%translations(:,i) - c%cen(:,j)) < symprec)) then
                found = .true.
                exit
             end if
          end do
          if (.not.found) then
             c%ncv = c%ncv + 1
             if (c%ncv > size(c%cen,2)) &
                call realloc(c%cen,3,2*c%ncv)
             c%cen(:,c%ncv) = c%spg%translations(:,i)
          end if
          cycle
       end if

       ! Do I have this rotation already?
       found = .false.
       do j = 1, c%neqv
          if (all(abs(rotm - c%rotm(1:3,1:3,j)) < symprec)) then
             found = .true.
          end if
       end do
       if (.not.found) then
          c%neqv = c%neqv + 1
          c%rotm(1:3,1:3,c%neqv) = rotm
          c%rotm(1:3,4,c%neqv) = c%spg%translations(:,i)
       end if
    end do
    call realloc(c%cen,3,c%ncv)
    c%havesym = 1

    ! generate symmetry operation info for the complete atom list
    ! use the spg%equivalent_atoms information
    do k = 1, c%ncel
       idx = c%spg%equivalent_atoms(k) + 1
       iat = iidx(idx)

       found = .false.
       loopi: do i = 1, c%neqv
          do j = 1, c%ncv
             x0 = matmul(c%rotm(1:3,1:3,i),c%at(iat)%x) + c%rotm(:,4,i) + c%cen(:,j)
             if (c%are_lclose(x0,c%atcel(k)%x,2d0*symprec)) then
                c%atcel(k)%ir = i
                c%atcel(k)%ic = j
                c%atcel(k)%lvec = nint(c%atcel(k)%x - x0)
                found = .true.
                exit loopi
             end if
          end do
       end do loopi
       if (.not.found) then
          errmsg = "error building rotation and center information for atom list"
          return
       end if
    end do
    deallocate(iidx)

    ! write down the Wyckoff positions from the spg
    call c%spgtowyc(c%spg)

  end subroutine calcsym

  !> Clear the symmetry information from the crystal and rebuild the
  !> atom list. If cel2neq, copy the atom information from atcel(:) to
  !> at(:). If neq2cel, copy from at(:) to atcel(:). Note that if
  !> neq2cel is used, then the environment and asterisms become
  !> obsolete and need to be recalculated.
  module subroutine clearsym(c,cel2neq,neq2cel)
    use types, only: neqatom, realloc
    use param, only: eyet
    class(crystal), intent(inout) :: c
    logical, intent(in), optional :: cel2neq
    logical, intent(in), optional :: neq2cel

    type(neqatom), allocatable :: aux(:)
    integer :: idir, i

    idir = 0
    if (present(cel2neq)) then
       if (cel2neq) idir = 1
    end if
    if (present(neq2cel)) then
       if (neq2cel) idir = -1
    end if

    ! nullify the space group and all the symmetry info
    c%spgavail = .false.
    c%spg%n_atoms = 0
    c%havesym = 0
    c%neqv = 1
    c%rotm = 0d0
    c%rotm(:,:,1) = eyet
    c%ncv = 1
    if (allocated(c%cen)) deallocate(c%cen)
    allocate(c%cen(3,4))
    c%cen = 0d0
    if (allocated(c%at) .and. c%nneq > 0) then
       do i = 1, c%nneq
          c%at(i)%mult = 1
          c%at(i)%wyc = "?"
       end do
    end if

    if (idir == 1) then
       ! convert ncel to nneq
       if (c%nneq /= c%ncel) then
          allocate(aux(c%nneq))
          aux = c%at(1:c%nneq)
          c%nneq = c%ncel
          call realloc(c%at,c%ncel)
          do i = 1, c%ncel
             c%at(i) = aux(c%atcel(i)%idx)
             c%at(i)%x = c%atcel(i)%x
             c%at(i)%is = c%atcel(i)%is
             c%at(i)%mult = 1
             c%at(i)%wyc = "?"
          end do
          deallocate(aux)
       end if
    elseif (idir == -1) then
       ! convert nneq to ncel
       c%ncel = c%nneq
       if (allocated(c%atcel)) deallocate(c%atcel)
       allocate(c%atcel(c%ncel))
       do i = 1, c%ncel
          c%atcel(i)%x = c%at(i)%x
          c%atcel(i)%r = c%at(i)%r
          c%atcel(i)%idx = i
          c%atcel(i)%cidx = i
          c%atcel(i)%ir = 1
          c%atcel(i)%ic = 1
          c%atcel(i)%lvec = 0
          c%atcel(i)%is = c%at(i)%is
       end do
    end if

  end subroutine clearsym

  module subroutine checkgroup(c)
    use tools_math, only: matinv
    use tools_io, only: uout, string
    class(crystal), intent(inout) :: c

    integer :: i, j, i1, j1, ifound(2)
    real*8 :: op(4,4), op1(4,4)
    real*8, parameter :: eye4(4,4) = reshape((/1d0,0d0,0d0,0d0,0d0,1d0,0d0,0d0,0d0,0d0,1d0,0d0,0d0,0d0,0d0,1d0/),shape(eye4))

    write (uout,'("* Checking the consistency of the space group")')

    ! list of operations
    write (uout,'("+ List of operations")')
    do i = 1, c%neqv
       do  j = 1, c%ncv
          call assignop(i,j,op)
          write (uout,'("  Operation ",A,"/",A)') string(i), string(j)
          op(1:3,4) = c%x2c(op(1:3,4)) * 0.52917720859d0
          write (uout,'(3("    ",4(A," ")/),"    ",4(A," "))') &
             ((string(op(i1,j1),'f',length=9,decimal=6,justify=3), j1=1,4), i1=1,4)
       enddo
    end do

    ! identity
    call checkexists(eye4,ifound)
    if (all(ifound > 0)) then
       write (uout,'("+ Is there an identity operation? ... yes")')
    else
       write (uout,'("+ Is there an identity operation? ... no")')
       goto 999
    end if

    ! inverse
    write (uout,'("+ Inverse operations")')
    do i = 1, c%neqv
       do  j = 1, c%ncv
          call assignop(i,j,op)
          call matinv(op,4)
          call checkexists(op,ifound)
          if (all(ifound > 0)) then
             write (uout,'("  Op. ",A,"/",A," has inverse: ",A,"/",A)') string(i), string(j),&
                string(ifound(1)), string(ifound(2))
          else
             write (uout,'("  Op. ",A,"/",A," has inverse: not found!")') string(i), string(j)
             goto 999
          end if
       end do
    end do

    ! closure
    write (uout,'("+ Closure")')
    do i = 1, c%neqv
       do  j = 1, c%ncv
          call assignop(i,j,op)
          do i1 = 1, c%neqv
             do  j1 = 1, c%ncv
                call assignop(i1,j1,op1)
                op1 = matmul(op,op1)
                call checkexists(op1,ifound)
                if (all(ifound > 0)) then
                   write (uout,'("  Op. ",A,"/",A," times ",A,"/",A," is in the closure: ",A,"/",A)') &
                      string(i), string(j), string(i1), string(j1), string(ifound(1)), string(ifound(2))
                else
                   write (uout,'("  Op. ",A,"/",A," times ",A,"/",A," is not in the closure")') &
                      string(i), string(j), string(i1), string(j1)
                   goto 999
                end if
             end do
          end do
       end do
    end do

999 continue
    write (uout,*)

  contains
    subroutine assignop(ieqv,icv,op)
      real*8, intent(out) :: op(4,4)
      integer, intent(in) :: ieqv, icv

      op = 0d0
      op(1:3,1:3) = c%rotm(:,1:3,ieqv)
      op(1:3,4) = c%rotm(:,4,ieqv) + c%cen(:,icv)
      op(4,4) = 1d0

    end subroutine assignop
    subroutine checkexists(op,ifound)
      real*8, intent(in) :: op(4,4)
      integer, intent(out) :: ifound(2)

      real*8 :: op1(4,4), x(3)
      integer :: i, j

      ifound = 0
      do i = 1, c%neqv
         do  j = 1, c%ncv
            call assignop(i,j,op1)
            x = op(1:3,4) - op1(1:3,4)
            x = x - nint(x)
            if (all(abs(op(1:3,1:3) - op1(1:3,1:3)) < 1d-4) .and. all(abs(x) < 1d-4)) then
               ifound(1) = i
               ifound(2) = j
               return
            end if
         end do
      end do

    end subroutine checkexists
  end subroutine checkgroup

  !> Re-assign atomic types to have an asymmetric unit with whole molecules
  module subroutine wholemols(c,ti)
    use crystalseedmod, only: crystalseed
    use tools, only: qcksort
    use tools_io, only: ferror, faterr
    use types, only: realloc
    use param, only: icrd_crys
    class(crystal), intent(inout) :: c
    type(thread_info), intent(in), optional :: ti

    logical, allocatable :: ismap(:,:), ldone(:)
    integer :: i, j, k, l, id
    real*8 :: x0(3)
    type(crystalseed) :: ncseed
    logical, allocatable :: sgroup(:,:), agroup(:), isuse(:)
    integer, allocatable :: pr(:,:), iorder(:), idxi(:)
    real*8 :: rot(3,4), xd(3), xxd(3)
    logical :: found, ok, again
    integer :: ig, ngroup, ngroupold

    real*8, parameter :: eps = 1d-4

    ! this is only appropriate for molecular crystals
    if (c%ismolecule) return
    if (.not.c%ismol3d) return

    ! compute the product table for the rotational parts of this group
    allocate(pr(c%neqv,c%neqv))
    do i = 1, c%neqv
       do j = 1, c%neqv
          if (i == 1) then
             pr(i,j) = j
          else if (j == 1) then
             pr(i,j) = i
          else
             ! compute the product of i and j
             rot(:,1:3) = matmul(c%rotm(1:3,1:3,i),c%rotm(1:3,1:3,j))
             rot(:,4) = matmul(c%rotm(1:3,1:3,i),c%rotm(1:3,4,j)) + c%rotm(1:3,4,i)

             ! identify the product operation
             found = .false.
             loopk: do k = 1, c%neqv
                ok = all(abs(rot(1:3,1:3) - c%rotm(1:3,1:3,k)) < eps)
                if (ok) then
                   xd = rot(1:3,4) - c%rotm(1:3,4,k)
                   do l = 1, c%ncv
                      xxd = xd - c%cen(:,l)
                      if (all(abs(xxd - nint(xxd)) < eps)) then
                         found = .true.
                         pr(i,j) = k
                         exit loopk
                      end if
                   end do
                end if
             end do loopk
             if (.not.found) &
                call ferror("wholemols","symmetry operations do not form a group",faterr)
          end if
       end do
    end do

    ! initialize the subgroup list with the trivial subgroups
    allocate(sgroup(c%neqv,2),agroup(c%neqv))
    ngroup = 2
    sgroup = .true.
    sgroup(2:,2) = .false.

    ! construct all the other subgroups
    ngroupold = 1
    do while (ngroupold < ngroup)
       ! consider the next subgroup
       ngroupold = ngroupold + 1
       ig = ngroupold

       ! try adding operations one by one
       do i = 2, c%neqv
          agroup = sgroup(:,ig)
          if (agroup(i)) cycle
          agroup(i) = .true.

          ! check if it is a subset of a known non-trivial subgroup
          found = .false.
          do j = 3, ngroup
             if (all((agroup .and. sgroup(:,j)) .eqv. agroup)) then
                found = .true.
                exit
             end if
          end do
          if (found) cycle

          ! complete the subgroup
          again = .true.
          do while (again)
             again = .false.
             do k = 1, c%neqv
                if (.not.agroup(k)) cycle
                do j = 1, c%neqv
                   if (.not.agroup(j)) cycle
                   if (.not.agroup(pr(k,j))) then
                      agroup(pr(k,j)) = .true.
                      again = .true.
                   end if
                end do
             end do
          end do

          ! check it is not the full group
          if (all(agroup)) cycle

          ! add it to the list
          ngroup = ngroup + 1
          if (ngroup > size(sgroup,2)) call realloc(sgroup,c%neqv,2*ngroup)
          sgroup(:,ngroup) = agroup
       end do
    end do
    call realloc(sgroup,c%neqv,ngroup)
    deallocate(agroup,pr)

    ! sort the subgroups by order
    allocate(iorder(ngroup),idxi(ngroup))
    do i = 1, ngroup
       iorder(i) = count(sgroup(:,i))
       idxi(i) = i
    end do
    call qcksort(iorder,idxi,1,ngroup)
    sgroup = sgroup(:,idxi)
    deallocate(iorder,idxi)

    ! identify the largest known subgroup with whole molecules in the asymmetric unit
    allocate(ismap(c%ncv,c%neqv),ldone(c%nneq))
    ig = 0
    do l = ngroup,1,-1
       ismap = .false.
       do i = 2, c%neqv
          if (.not.sgroup(i,l)) cycle
          do j = 1, c%ncv

             ldone = .false.
             do k = 1, c%ncel
                if (ldone(c%atcel(k)%idx)) cycle

                x0 = matmul(c%rotm(1:3,1:3,i),c%atcel(k)%x) + c%rotm(:,4,i) + c%cen(:,j)
                id = c%identify_atom(x0,icrd_crys)

                if (id == 0) &
                   call ferror('wholemols','error identifying rotated atom',faterr)

                if (c%idatcelmol(k) == c%idatcelmol(id)) &
                   ismap(j,i) = .true.
                ldone(c%atcel(k)%idx) = .true.
             end do
          end do
       end do

       if (all(.not.ismap)) then
          ig = l
          exit
       end if
    end do
    deallocate(ldone,ismap)

    ! calculate the asymmetric unit for the new group, write it down in the seed
    allocate(ncseed%x(3,c%ncel),ncseed%is(c%ncel))
    ncseed%nat = 0
    allocate(isuse(c%ncel))
    isuse = .false.
    do i = 1, c%ncel
       if (isuse(i)) cycle
       isuse(i) = .true.
       do j = 1, c%neqv
          if (.not.sgroup(j,ig)) cycle
          do k = 1, c%ncv
             x0 = matmul(c%rotm(1:3,1:3,j),c%atcel(i)%x) + c%rotm(:,4,j) + c%cen(:,k)
             id = c%identify_atom(x0,icrd_crys)
             if (id == 0) &
                call ferror('wholemols','error identifying rotated atom',faterr)
             isuse(id) = .true.
          end do
       end do

       ncseed%nat = ncseed%nat + 1
       ncseed%x(:,ncseed%nat) = c%atcel(i)%x
       ncseed%is(ncseed%nat) = c%atcel(i)%is
    end do
    deallocate(isuse)

    ! write down the new symmetry operations in the seed
    ncseed%ncv = c%ncv
    allocate(ncseed%cen(3,c%ncv))
    ncseed%cen = c%cen(1:3,1:c%ncv)
    ncseed%neqv = count(sgroup(:,ig))
    allocate(ncseed%rotm(3,4,ncseed%neqv))
    k = 0
    do i = 1, c%neqv
       if (sgroup(i,ig)) then
          k = k + 1
          ncseed%rotm(:,:,k) = c%rotm(:,:,i)
       end if
    end do
    deallocate(sgroup)

    ! write the rest of the seed
    allocate(ncseed%spc(c%nspc))
    ncseed%nspc = c%nspc
    ncseed%spc = c%spc
    ncseed%useabr = 2
    ncseed%m_x2c = c%m_x2c
    ncseed%ismolecule = c%ismolecule
    ncseed%isused = .true.
    ncseed%file = c%file
    ncseed%havesym = 1
    ncseed%findsym = 0
    ncseed%checkrepeats = .false.

    ! build the new crystal
    call c%struct_new(ncseed,.true.,ti=ti)

  end subroutine wholemols

  !> Calculate the irreducible WS wedge around point xorigin (cryst coords)
  !> and partition it into tetrahedra.
  module subroutine getiws(c,xorigin,ntetrag,tetrag)
    use tools_math, only: mixed
    use types, only: realloc
    class(crystal), intent(in) :: c !< the crystal structure
    real*8, intent(in) :: xorigin(3)
    integer, intent(out), optional :: ntetrag
    real*8, allocatable, intent(inout), optional :: tetrag(:,:,:)

    integer :: leqv
    real*8 :: lrotm(3,3,48)
    logical, allocatable :: active(:)
    character*3 :: pg
    integer :: i, j, n, m
    real*8 :: x0(3), xp1(3), xp2(3), xp3(3), xoriginc(3)

    real*8, parameter :: eps_wspesca = 1d-5 !< Criterion for tetrahedra equivalence
    real*8, parameter :: ws_eps_vol = 1d-5 !< Reject tetrahedra smaller than this.

    ! origin in Cartesian
    xoriginc = c%x2c(xorigin)

    ! local symmetry group
    pg = c%sitesymm(xorigin,leqv=leqv,lrotm=lrotm)

    ! count tetrahedra
    ntetrag = 0
    do i = 1, c%ws_nf
       ntetrag = ntetrag + 2 * c%ws_nside(i)
    end do

    ! build all the tetrahedra
    m = 0
    if (allocated(tetrag)) deallocate(tetrag)
    allocate(tetrag(3,4,ntetrag))
    do i = 1, c%ws_nf
       n = c%ws_nside(i)
       ! calculate middle point of the face
       x0 = 0d0
       do j = 1, n
          x0 = x0 + c%ws_x(:,c%ws_iside(j,i))
       end do
       x0 = c%x2c(x0) / n

       do j = 1, n
          xp1 = c%x2c(c%ws_x(:,c%ws_iside(j,i)))
          xp2 = c%x2c(c%ws_x(:,c%ws_iside(mod(j,n)+1,i)))
          ! tetrah 1
          m = m + 1
          tetrag(:,1,m) = (/ 0d0, 0d0, 0d0 /) + xoriginc
          tetrag(:,2,m) = x0 + xoriginc
          tetrag(:,3,m) = xp1 + xoriginc
          tetrag(:,4,m) = 0.5d0 * (xp1 + xp2) + xoriginc
          ! tetrah 2
          m = m + 1
          tetrag(:,1,m) = (/ 0d0, 0d0, 0d0 /) + xoriginc
          tetrag(:,2,m) = x0 + xoriginc
          tetrag(:,3,m) = xp2 + xoriginc
          tetrag(:,4,m) = 0.5d0 * (xp1 + xp2) + xoriginc
       end do
    end do

    ! check volumes and convert to crystallographic
    allocate(active(ntetrag))
    active = .true.
    do i = 1, ntetrag
       ! calculate volume
       xp1 = tetrag(:,2,i) - tetrag(:,1,i)
       xp2 = tetrag(:,3,i) - tetrag(:,1,i)
       xp3 = tetrag(:,4,i) - tetrag(:,1,i)
       if (abs(mixed(xp1,xp2,xp3)) / 6d0 < ws_eps_vol) then
          active(i) = .false.
          cycle
       end if
       do j = 1, 4
          tetrag(:,j,i) = matmul(c%m_c2x,tetrag(:,j,i))
       end do
    end do

    ! reduce tetrahedra equivalent by symmetry
    do i = 1, ntetrag
       if (.not.active(i)) cycle
       do j = i+1, ntetrag
          if (equiv_tetrah(c,xorigin,tetrag(:,:,i),tetrag(:,:,j),leqv,lrotm,eps_wspesca)) then
             active(j) = .false.
          end if
       end do
    end do

    ! rebuild the tetrahedra list
    in: do i = 1, ntetrag
       if (active(i)) cycle
       do j = i+1, ntetrag
          if (active(j)) then
             tetrag(:,:,i) = tetrag(:,:,j)
             active(j) = .false.
             cycle in
          end if
       end do
       ntetrag = i - 1
       exit
    end do in
    call realloc(tetrag,3,4,ntetrag)

    deallocate(active)

  end subroutine getiws

  !> Calculate the number of lattice vectors in each direction in
  !> order to be sure that the main cell is surrounded by a shell at
  !> least rmax thick. x2r is the cryst-to-car matrix for the
  !> lattice. The heuristic method has been adapted from gulp.
  module subroutine search_lattice(x2r,rmax,imax,jmax,kmax)
    use param, only: pi

    real*8, intent(in) :: x2r(3,3), rmax
    integer, intent(out) :: imax, jmax, kmax

    integer :: i, nadd
    real*8 :: acell(3), bcell(3), gmat(3,3)

    gmat = matmul(transpose(x2r),x2r)
    do i = 1, 3
       acell(i) = sqrt(gmat(i,i))
    end do
    bcell(1) = acos(gmat(2,3) / acell(2) / acell(3)) * 180d0 / pi
    bcell(2) = acos(gmat(1,3) / acell(1) / acell(3)) * 180d0 / pi
    bcell(3) = acos(gmat(1,2) / acell(1) / acell(2)) * 180d0 / pi

    ! determine cell list, trick adapted from gulp
    nadd = 0
    if (bcell(1)<30 .or. bcell(2)<30 .or. bcell(3)<30 .or. bcell(1)>150 .or. bcell(2)>150 .or. bcell(3)>150) then
       nadd = 3
    else if (bcell(1)<50 .or. bcell(2)<50 .or. bcell(3)<50 .or. bcell(1)>130 .or. bcell(2)>130 .or. bcell(3)>130) then
       nadd = 2
    else if (bcell(1)<70 .or. bcell(2)<70 .or. bcell(3)<70 .or. bcell(1)>110 .or. bcell(2)>110 .or. bcell(3)>110) then
       nadd = 1
    else
       nadd = 1
    end if
    imax = ceiling(rmax / acell(1)) + nadd
    jmax = ceiling(rmax / acell(2)) + nadd
    kmax = ceiling(rmax / acell(3)) + nadd

  end subroutine search_lattice

  !> Write the structure to a file. Use the format derived from the
  !> extension of file and use default values for all options.
  module subroutine write_simple_driver(c,file,ti)
    use tools_io, only: lower, ferror, faterr, equal
    class(crystal), intent(inout) :: c
    character*(*), intent(in) :: file
    type(thread_info), intent(in), optional :: ti

    character(len=:), allocatable :: wext, wext2, wroot
    integer :: idx

    wext = lower(file(index(file,'.',.true.)+1:))
    wroot = file(:index(file,'.',.true.)-1)

    if (equal(wext,'xyz').or.equal(wext,'gjf').or.equal(wext,'cml')) then
       call c%write_mol(file,wext,ti=ti)
    else if(equal(wext,'obj').or.equal(wext,'ply').or.equal(wext,'off')) then
       call c%write_3dmodel(file,wext,ti=ti)
    elseif (equal(wext,'gau')) then
       call c%write_gaussian(file,ti=ti)
    elseif (equal(wext,'in')) then
       idx = index(wroot,'.',.true.)
       if (idx > 0) then
          wext2 = lower(wroot(idx+1:))
          if (equal(wext2,'scf')) then
             call c%write_espresso(file,ti=ti)
          else
             idx = 0
          end if
       end if
       if (idx == 0) &
            call c%write_fhi(file,.true.,ti=ti)
    elseif (equal(wext,'poscar') .or. equal(wext,'contcar')) then
       call c%write_vasp(file,.false.,ti=ti)
    elseif (equal(wext,'abin')) then
       call c%write_abinit(file,ti=ti)
    elseif (equal(wext,'elk')) then
       call c%write_elk(file,ti=ti)
    elseif (equal(wext,'tess')) then
       call c%write_tessel(file,ti=ti)
    elseif (equal(wext,'incritic').or.equal(wext,'cri')) then
       call c%write_critic(file,ti=ti)
    elseif (equal(wext,'cif')) then
       call c%write_cif(file,.true.,ti=ti)
    elseif (equal(wext,'d12').or.equal(wext,'34')) then
       call c%write_d12(file,.true.,ti=ti)
    elseif (equal(wext,'res')) then
       call c%write_res(file,-1,ti=ti)
    elseif (equal(wext,'m')) then
       call c%write_escher(file,ti=ti)
    elseif (equal(wext,'db')) then
       call c%write_db(file,ti=ti)
    elseif (equal(wext,'gin')) then
       call c%write_gulp(file,ti=ti)
    elseif (equal(wext,'lammps')) then
       call c%write_lammps(file,ti=ti)
    elseif (equal(wext,'fdf')) then
       call c%write_siesta_fdf(file,ti=ti)
    elseif (equal(wext,'struct_in')) then
       call c%write_siesta_in(file,ti=ti)
    elseif (equal(wext,'hsd')) then
       call c%write_dftbp_hsd(file,ti=ti)
    elseif (equal(wext,'gen')) then
       call c%write_dftbp_gen(file,ti=ti)
    elseif (equal(wext,'pyscf')) then
       call c%write_pyscf(file,ti=ti)
    elseif (equal(wext,'fhi')) then
       call c%write_fhi(file,.true.,ti=ti)
    elseif (equal(wext,'frac')) then
       call c%write_tinkerfrac(file,ti=ti)
    else
       call ferror('struct_write','unrecognized file format',faterr)
    end if

  end subroutine write_simple_driver

  !> Write a xyz/gjf/cml file containing a finite piece of the crystal
  !> structure. fmt can be one of xyz, gjf, or cml. ix is the number
  !> of unit cells to plot.  If doborder is .true., add all atoms at
  !> the border. If onemotif is .true., write all molecules in the
  !> unit cell.  If molmotif is .true., complete molecules with atoms
  !> in adjacent cells. If docell, add sticks for the unit cell
  !> limits. If environ is true, write all molecules up to a distance
  !> renv (bohr) from the origin. If lnmer, partition the resulting
  !> list of molecules into n-mers, up to nth order. If rsph (bohr) is
  !> positive, then use all atoms in a sphere around xsph (cryst.). If
  !> rcub (bohr) is positive, use all atoms in a cube around xcub
  !> (cryst.). If luout is present, return the LU in that argument
  !> and do not close the file.
  module subroutine write_mol(c,file,fmt,ix0,doborder0,onemotif0,molmotif0,&
     environ0,renv0,lnmer0,nmer0,rsph0,xsph0,rcub0,xcub0,usenames0,luout,ti)
    use global, only: dunit0, iunit
    use tools_math, only: nchoosek, comb
    use tools_io, only: ferror, faterr, uout, string, ioj_left, string, ioj_right,&
       equal
    use types, only: realloc
    use fragmentmod, only: fragment, realloc_fragment
    class(crystal), intent(inout) :: c
    character*(*), intent(in) :: file
    character*3, intent(in) :: fmt
    integer, intent(in), optional :: ix0(3)
    logical, intent(in), optional :: doborder0, onemotif0, molmotif0, environ0
    real*8, intent(in), optional :: renv0
    logical, intent(in), optional :: lnmer0
    integer, intent(in), optional :: nmer0
    real*8, intent(in), optional :: rsph0, xsph0(3)
    real*8, intent(in), optional :: rcub0, xcub0(3)
    logical, intent(in), optional :: usenames0
    integer, intent(out), optional :: luout
    type(thread_info), intent(in), optional :: ti

    type(fragment) :: fr
    type(fragment), allocatable :: fr0(:)
    logical, allocatable :: isdiscrete(:)
    integer :: i, j, k, l, m, nmol, icel, lvec(3), ncm
    integer :: ncomb, nlimi, nlimj, icount
    integer, allocatable :: icomb(:), origmol(:)
    character(len=:), allocatable :: wroot, file0, aux
    logical :: doagain
    real*8, allocatable :: cmlist(:,:)
    real*8 :: xcm(3), dist
    integer :: ix(3)
    logical :: doborder, onemotif, molmotif, environ
    real*8 :: renv
    logical :: lnmer, usenames
    integer :: nmer
    real*8 :: rsph, xsph(3)
    real*8 :: rcub, xcub(3)

    ! set the default options
    ix = 1
    doborder = .false.
    onemotif = .false.
    molmotif = .false.
    environ = .false.
    renv = 0d0
    lnmer = .false.
    nmer = 1
    rsph = -1d0
    xsph = 0d0
    rcub = -1d0
    xcub = 0d0
    usenames = .false.
    if (present(ix0)) ix = ix0
    if (present(doborder0)) doborder = doborder0
    if (present(onemotif0)) onemotif = onemotif0
    if (present(molmotif0)) molmotif = molmotif0
    if (present(environ0)) environ = environ0
    if (present(renv0)) renv = renv0
    if (present(lnmer0)) lnmer = lnmer0
    if (present(nmer0)) nmer = nmer0
    if (present(rsph0)) rsph = rsph0
    if (present(xsph0)) xsph = xsph0
    if (present(rcub0)) rcub = rcub0
    if (present(xcub0)) xcub = xcub0
    if (present(usenames0)) usenames = usenames0

    ! determine the fragments
    if (onemotif) then
       call fr%merge_array(c%mol(1:c%nmol),.false.)
       allocate(fr0(c%nmol))
       fr0 = c%mol
       nmol = c%nmol
    elseif (environ) then
       ! calculate the centers of mass for all fragments in the molecular motif
       allocate(cmlist(3,c%nmol),fr0(c%nmol),origmol(c%nmol))
       nmol = c%nmol
       fr0 = c%mol
       do i = 1, c%nmol
          cmlist(:,i) = c%mol(i)%cmass()
          origmol(i) = i
       end do
       ncm = c%nmol

       doagain = .true.
       icel = 0
       do while(doagain)
          doagain = .false.
          icel = icel + 1
          do k = 1, 6
             if (k == 1 .or. k == 2) then
                nlimi = icel
                nlimj = icel
             elseif (k == 3 .or. k == 4) then
                nlimi = icel
                nlimj = icel-1
             else
                nlimi = icel-1
                nlimj = icel-1
             end if
             do i = -nlimi, nlimi
                do j = -nlimj, nlimj
                   lvec = icelcomb(k,i,j,icel)
                   do l = 1, c%nmol
                      xcm = c%c2x(cmlist(:,l)) + lvec
                      xcm = c%x2c(xcm)
                      dist = norm2(xcm)
                      if (dist <= renv) then
                         ncm = ncm + 1
                         if (ncm > size(cmlist,2)) then
                            call realloc(cmlist,3,2*ncm)
                            call realloc(origmol,2*ncm)
                            call realloc_fragment(fr0,2*ncm)
                         end if
                         cmlist(:,ncm) = xcm
                         fr0(ncm) = c%mol(l)
                         do m = 1, fr0(ncm)%nat
                            fr0(ncm)%at(m)%x = fr0(ncm)%at(m)%x + lvec
                            fr0(ncm)%at(m)%r = c%x2c(fr0(ncm)%at(m)%x)
                            fr0(ncm)%at(m)%lvec = fr0(ncm)%at(m)%lvec + lvec
                         end do
                         origmol(ncm) = l
                         doagain = .true.
                      end if
                   end do
                end do
             end do
          end do
       end do

       call realloc_fragment(fr0,ncm)
       call fr%merge_array(fr0,.false.)
       nmol = ncm
    else
       if (rcub > 0) then
          fr = c%listatoms_sphcub(rcub=rcub,xcub=xcub)
       elseif (rsph > 0) then
          fr = c%listatoms_sphcub(rsph=rsph,xsph=xsph)
       else
          fr = c%listatoms_cells(ix,doborder)
       endif
       if (molmotif) then
          call c%listmolecules(fr,nmol,fr0,isdiscrete)
          call fr%merge_array(fr0,.false.)
       end if
    end if

    if (lnmer .and..not.allocated(fr0)) &
       call ferror('write_mol','ONEMOTIF, MOLMOTIF, or ENVIRON are necessary with NMER',faterr)

    ! If environ, report the identities of all the molecules in the environment
    if (environ) then
       write (uout,'("+ List of fragments in the molecular environment")')
       write (uout,'("  Number of fragments: ",A)') string(nmol)
       write (uout,'("  Fragment number Id with nat atoms at center-of-mass comes from")')
       write (uout,'("  from fragment idmol in the Wigner-Seitz cell translated by")')
       write (uout,'("  lattice vector lvec.")')
       write (uout,'("# Id nat               Center of mass          idmol     lvec")')
       do i = 1, nmol
          if (c%ismolecule) then
             xcm = cmlist(:,i) * dunit0(iunit)
          else
             xcm = c%c2x(cmlist(:,i))
          end if
          write (uout,'(99("  ",A))') string(i,3,ioj_left), string(fr0(i)%nat,4,ioj_left),&
             (string(xcm(j),'f',10,6,3),j=1,3), string(origmol(i)), &
             (string(nint(fr0(i)%at(1)%x(j)-c%mol(origmol(i))%at(1)%x(j)),3,ioj_right),j=1,3)
       end do
       write (uout,*)
    end if

    ! if this is a molecule, translate to the proper origin
    if (c%ismolecule) then
       do i = 1, fr%nat
          fr%at(i)%r = fr%at(i)%r + c%molx0
       end do
       if (allocated(fr0)) then
          do i = 1, nmol
             do j = 1, fr0(i)%nat
                fr0(i)%at(j)%r = fr0(i)%at(j)%r + c%molx0
             end do
          end do
       end if
    end if

    if (.not.lnmer) then
       ! normal write
       call dowrite(file,fr,ti=ti)
    else
       wroot = file(:index(file,'.',.true.)-1)
       do i = 1, nmer
          if (i == 1) then
             if (nmer == 1) then
                nlimj = nmol
             else
                nlimj = c%nmol
             end if
             ! monomers
             do j = 1, nlimj
                file0 = trim(wroot) // "_" // string(j) // "." // fmt
                call dowrite(file0,fr0(j),ti=ti)
             end do
             write (uout,'("+ Written ",A," ",A,"-mers")') string(c%nmol), string(i)
          elseif (i == nmer) then
             ! n-mers
             allocate(icomb(i-1))
             icount = 0
             do l = 1, c%nmol
                ncomb = nchoosek(nmol,i-1)
                do j = 1, ncomb
                   call comb(nmol,i-1,j,icomb)
                   if (any(icomb == l)) cycle
                   file0 = trim(wroot) // "_" // string(l)
                   fr = fr0(l)
                   do k = 1, i-1
                      aux = trim(file0) // "_" // string(icomb(k))
                      file0 = aux
                      call fr%append(fr0(icomb(k)))
                   end do
                   aux = trim(file0) // "." // fmt
                   file0 = aux
                   call dowrite(file0,fr,ti=ti)
                   icount = icount + 1
                end do
             end do
             deallocate(icomb)
             write (uout,'("+ Written ",A," ",A,"-mers")') string(icount), string(i)
          else
             ! everything in between
             ncomb = nchoosek(nmol,i)
             allocate(icomb(i))
             icount = 0
             do j = 1, ncomb
                call comb(nmol,i,j,icomb)
                file0 = wroot
                call fr%init()
                do k = 1, i
                   aux = trim(file0) // "_" // string(icomb(k))
                   file0 = aux
                   call fr%append(fr0(icomb(k)))
                end do
                aux = trim(file0) // "." // fmt
                file0 = aux
                call dowrite(file0,fr,ti=ti)
                icount = icount + 1
             end do
             deallocate(icomb)
             write (uout,'("+ Written ",A," ",A,"-mers")') string(icount), string(i)
          end if
       end do
    end if

  contains
    subroutine dowrite(fileo,fro,ti)
      character*(*) :: fileo
      type(fragment) :: fro
      type(thread_info), intent(in), optional :: ti

      if (equal(fmt,"xyz")) then
         call fro%writexyz(fileo,usenames,ti=ti)
      elseif (equal(fmt,"gjf")) then
         call fro%writegjf(fileo,ti=ti)
      elseif (equal(fmt,"cml")) then
         if (c%ismolecule) then
            call fro%writecml(fileo,luout=luout,ti=ti)
         else
            call fro%writecml(fileo,c%m_x2c,luout=luout,ti=ti)
         end if
      else
         call ferror("write_mol","Unknown format",faterr)
      endif

    end subroutine dowrite
    function icelcomb(idx,i,j,icel)
      integer, intent(in) :: idx, i, j, icel
      integer :: icelcomb(3)

      icelcomb = 0
      if (idx == 1) then
         icelcomb = (/i,j,icel/)
      elseif (idx == 2) then
         icelcomb = (/i,j,-icel/)
      elseif (idx == 3) then
         icelcomb = (/i,icel,j/)
      elseif (idx == 4) then
         icelcomb = (/i,-icel,j/)
      elseif (idx == 5) then
         icelcomb = (/icel,i,j/)
      elseif (idx == 6) then
         icelcomb = (/-icel,i,j/)
      endif

    end function icelcomb
  end subroutine write_mol

  !> Write an obj/ply/off file containing the crystal structure. fmt
  !> can be one of obj, ply, or off. ix is the number of unit cells to
  !> plot. If doborder is .true., add all atoms at the border. If
  !> onemotif is .true., write all molecules in the unit cell. If
  !> molmotif is .true., complete molecules with atoms in adjacent
  !> cells. If docell, add sticks for the unit cell limits. If
  !> domolcell, add sticks tfor the molecular cell. If rsph (bohr) is
  !> positive, then use all atoms in a sphere around xsph (cryst.). If
  !> rcub (bohr) is positive, use all atoms in a cube around xcub
  !> (cryst.). If gr0 is present, then return the graphics handle and
  !> do not close the files.
  module subroutine write_3dmodel(c,file,fmt,ix0,doborder0,onemotif0,molmotif0,&
     docell0,domolcell0,rsph0,xsph0,rcub0,xcub0,gr0,ti)
    use graphics, only: grhandle
    use fragmentmod, only: fragment
    use tools_io, only: equal
    use param, only: maxzat, atmcov, jmlcol
    class(crystal), intent(inout) :: c
    character*(*), intent(in) :: file
    character*3, intent(in) :: fmt
    integer, intent(in), optional :: ix0(3)
    logical, intent(in), optional :: doborder0, onemotif0, molmotif0
    logical, intent(in), optional :: docell0, domolcell0
    real*8, intent(in), optional :: rsph0, xsph0(3)
    real*8, intent(in), optional :: rcub0, xcub0(3)
    type(grhandle), intent(out), optional :: gr0
    type(thread_info), intent(in), optional :: ti

    integer :: i, j
    real*8 :: d, xd(3), x0(3), x1(3), rr
    type(fragment) :: fr
    type(fragment), allocatable :: fr0(:)
    logical, allocatable :: isdiscrete(:)
    integer :: nmol
    type(grhandle) :: gr
    integer :: ix(3)
    logical :: doborder, onemotif, molmotif
    logical :: docell, domolcell
    real*8 :: rsph, xsph(3)
    real*8 :: rcub, xcub(3)

    real*8, parameter :: rfac = 1.4d0
    real*8, parameter :: x0cell(3,2,12) = reshape((/&
       0d0, 0d0, 0d0,   0d0, 0d0, 1d0,&
       0d0, 0d0, 0d0,   0d0, 1d0, 0d0,&
       0d0, 0d0, 0d0,   1d0, 0d0, 0d0,&
       1d0, 0d0, 0d0,   1d0, 1d0, 0d0,&
       1d0, 1d0, 0d0,   0d0, 1d0, 0d0,&
       0d0, 1d0, 0d0,   0d0, 1d0, 1d0,&
       0d0, 1d0, 1d0,   0d0, 0d0, 1d0,&
       0d0, 0d0, 1d0,   1d0, 0d0, 1d0,&
       1d0, 0d0, 1d0,   1d0, 0d0, 0d0,&
       1d0, 1d0, 1d0,   0d0, 1d0, 1d0,&
       1d0, 1d0, 1d0,   1d0, 0d0, 1d0,&
       1d0, 1d0, 1d0,   1d0, 1d0, 0d0/),shape(x0cell))

    ! set default values
    ix = 1
    doborder = .false.
    onemotif = .false.
    molmotif = .false.
    docell = .false.
    domolcell = .false.
    rsph = -1d0
    xsph = 0d0
    rcub = -1d0
    xcub = 0d0
    if (present(ix0)) ix = ix0
    if (present(doborder0)) doborder = doborder0
    if (present(onemotif0)) onemotif = onemotif0
    if (present(molmotif0)) molmotif = molmotif0
    if (present(docell0)) docell = docell0
    if (present(domolcell0)) domolcell = domolcell0
    if (present(rsph0)) rsph = rsph0
    if (present(xsph0)) xsph = xsph0
    if (present(rcub0)) rcub = rcub0
    if (present(xcub0)) xcub = xcub0

    ! open and get the atom list
    if (onemotif) then
       call fr%merge_array(c%mol(1:c%nmol),.false.)
    else
       if (rcub > 0) then
          fr = c%listatoms_sphcub(rcub=rcub,xcub=xcub)
       elseif (rsph > 0) then
          fr = c%listatoms_sphcub(rsph=rsph,xsph=xsph)
       else
          fr = c%listatoms_cells(ix,doborder)
       endif
       if (molmotif) then
          call c%listmolecules(fr,nmol,fr0,isdiscrete)
          call fr%merge_array(fr0,.false.)
       end if
    end if

    ! if this is a molecule, translate to the proper origin
    if (c%ismolecule) then
       do i = 1, fr%nat
          fr%at(i)%r = fr%at(i)%r + c%molx0
       end do
    end if

    call gr%open(fmt,file,ti=ti)

    ! add the balls
    do i = 1, fr%nat
       if (fr%spc(fr%at(i)%is)%z > maxzat) then
          rr = 0.21d0
       else
          rr = 0.6d0*atmcov(fr%spc(fr%at(i)%is)%z)
       endif
       call gr%ball(fr%at(i)%r,JMLcol(:,fr%spc(fr%at(i)%is)%z),rr)
    end do

    ! add the sticks
    do i = 1, fr%nat
       do j = i+1, fr%nat
          if (fr%spc(fr%at(i)%is)%z > maxzat .or. fr%spc(fr%at(j)%is)%z > maxzat) cycle
          xd = fr%at(i)%r - fr%at(j)%r
          d = norm2(xd)
          if (d < (atmcov(fr%spc(fr%at(i)%is)%z) + atmcov(fr%spc(fr%at(j)%is)%z)) * rfac) then
             xd = fr%at(i)%r + 0.5d0 * (fr%at(j)%r - fr%at(i)%r)
             call gr%stick(fr%at(i)%r,xd,JMLcol(:,fr%spc(fr%at(i)%is)%z),0.05d0)
             call gr%stick(fr%at(j)%r,xd,JMLcol(:,fr%spc(fr%at(j)%is)%z),0.05d0)
          end if
       end do
    end do

    ! add the cell
    if (docell) then
       do i = 1, 12
          x0 = c%x2c(x0cell(:,1,i)) + c%molx0
          x1 = c%x2c(x0cell(:,2,i)) + c%molx0
          call gr%stick(x0,x1,(/255,0,0/),0.03d0)
       end do
    end if

    ! add the molecular cell
    if (domolcell .and. c%ismolecule) then
       do i = 1, 12
          x0 = x0cell(:,1,i)
          x1 = x0cell(:,2,i)
          do j = 1, 3
             if (abs(x0(j)) < 1d-12) x0(j) = c%molborder(j)
             if (abs(x0(j)-1d0) < 1d-12) x0(j) = 1d0-c%molborder(j)
             if (abs(x1(j)) < 1d-12) x1(j) = c%molborder(j)
             if (abs(x1(j)-1d0) < 1d-12) x1(j) = 1d0-c%molborder(j)
          end do
          x0 = c%x2c(x0) + c%molx0
          x1 = c%x2c(x1) + c%molx0
          call gr%stick(x0,x1,(/0,0,255/),0.03d0)
       end do
    end if

    ! close or give the handles to the calling routine, cleanup
    if (present(gr0)) then
       gr0 = gr
    else
       call gr%close(ti=ti)
    end if

    if (allocated(fr0)) deallocate(fr0)

  end subroutine write_3dmodel

  !> Write a quantum espresso input template
  module subroutine write_espresso(c,file,rklength,ti)
    use tools_io, only: fopen_write, lower, fclose, string
    use param, only: atmass
    class(crystal), intent(in) :: c
    character*(*), intent(in) :: file
    real*8, intent(in), optional :: rklength
    type(thread_info), intent(in), optional :: ti

    integer :: i, lu, nk(3)
    real*8 :: rk

    rk = 40d0
    if (present(rklength)) rk = rklength

    lu = fopen_write(file,ti=ti)
    write (lu,'("&control")')
    write (lu,'(" title=''crystal'',")')
    write (lu,'(" prefix=''crystal'',")')
    write (lu,'(" pseudo_dir=''.'',")')
    write (lu,'(" calculation=''vc-relax'',")')
    write (lu,'("/")')
    write (lu,'("&system")')
    write (lu,'(" ibrav=0,")')
    write (lu,'(" nat=",I6,",")') c%ncel
    write (lu,'(" ntyp=",I3,",")') c%nspc
    write (lu,'(" ecutwfc=60.0,")')
    write (lu,'(" ecutrho=600.0,")')
    write (lu,'(" xdm=.true.,")')
    write (lu,'("/")')
    write (lu,'("&electrons"/" conv_thr = 1d-8,"/"/")')
    write (lu,'("&ions"/"/")')
    write (lu,'("&cell"/"/")')
    write (lu,'("ATOMIC_SPECIES")')
    do i = 1, c%nspc
       write (lu,'(A," ",F12.6," ",A,".UPF")') trim(c%spc(i)%name), atmass(c%spc(i)%z), trim(lower(c%spc(i)%name))
    end do
    write (lu,'(/"ATOMIC_POSITIONS crystal")')
    do i = 1, c%ncel
       write (lu,'(A,3(" ",F13.8," "))') trim(c%spc(c%atcel(i)%is)%name), c%atcel(i)%x
    end do

    call c%get_kpoints(rk,nk)
    write (lu,'(/"K_POINTS automatic"/3(A," ")" 1 1 1"/)') (string(nk(i)),i=1,3)

    write (lu,'("CELL_PARAMETERS bohr")')
    do i = 1, 3
       write (lu,'(3(F18.12," "))') c%m_x2c(:,i)
    end do
    call fclose(lu)

  end subroutine write_espresso

  !> Write a VASP POSCAR file. If verbose, write the atom sequence
  !> to the output. If append, append to an existing file.
  module subroutine write_vasp(c,file,verbose,append,ti)
    use tools_io, only: fopen_write, fopen_append, string, uout, fclose, nameguess
    use param, only: bohrtoa
    class(crystal), intent(in) :: c
    character*(*), intent(in) :: file
    logical, intent(in) :: verbose
    logical, intent(in), optional :: append
    type(thread_info), intent(in), optional :: ti

    character(len=:), allocatable :: lbl1, lbl2, aux, auxname
    integer :: i, j, lu, ntyp
    logical :: append_, ok

    ! whether to append
    append_ = .false.
    if (present(append)) append_ = append

    ! Cell
    if (append_) then
       inquire(file=file,exist=ok)
       if (ok) then
          lu = fopen_append(file,ti=ti)
       else
          lu = fopen_write(file,ti=ti)
       end if
    else
       lu = fopen_write(file,ti=ti)
    end if
    write (lu,'("critic2 | ",A)') string(c%file)
    write (lu,'("1.0")')
    do i = 1, 3
       write (lu,'(3(F15.10," "))') c%m_x2c(:,i) * bohrtoa
    end do

    ! Number of atoms per type and Direct
    lbl1 = ""
    lbl2 = ""
    do i = 1, c%nspc
       ntyp = 0
       do j = 1, c%ncel
          if (c%atcel(j)%is == i) ntyp = ntyp + 1
       end do
       aux = lbl1 // " " // string(ntyp)
       lbl1 = aux
       auxname = lbl2 // " " // trim(nameguess(c%spc(i)%z,.true.))
       lbl2 = auxname
    end do
    write (lu,'(A)') lbl2
    write (lu,'(A)') lbl1
    write (lu,'("Direct")')

    ! Atomic positions
    do i = 1, c%nspc
       do j = 1, c%ncel
          if (c%atcel(j)%is == i) then
             write (lu,'(3(F13.8," "))') c%atcel(j)%x
          end if
       end do
    end do
    call fclose(lu)

    if (verbose) &
       write (uout,'("+ Atom type sequence: ",999(A," "))') (string(c%spc(j)%name),j=1,c%nspc)

  end subroutine write_vasp

  !> Write an abinit input template
  module subroutine write_abinit(c,file,ti)
    use tools_io, only: fopen_write, string, fclose
    use param, only: pi
    class(crystal), intent(in) :: c
    character*(*), intent(in) :: file
    type(thread_info), intent(in), optional :: ti

    character(len=:), allocatable :: lbl1, aux
    integer :: ntyp
    real*8 :: aap(3), bbp(3), gpq(3,3)
    integer :: i, j, lu

    ! Find the lengths and angles of the cell
    gpq = matmul(transpose(c%m_x2c),c%m_x2c)
    do i = 1, 3
       aap(i) = sqrt(gpq(i,i))
    end do
    bbp(1)=acos(gpq(2,3)/sqrt(gpq(2,2)*gpq(3,3)))/pi*180d0
    bbp(2)=acos(gpq(1,3)/sqrt(gpq(1,1)*gpq(3,3)))/pi*180d0
    bbp(3)=acos(gpq(1,2)/sqrt(gpq(1,1)*gpq(2,2)))/pi*180d0

    ! Write input
    lu = fopen_write(file,ti=ti)
    write (lu,'("acell ",3(F14.10," "))') aap
    write (lu,'("angdeg ",3(F14.10," "))') bbp
    write (lu,'("ntypat ",I3)') c%nspc

    lbl1 = ""
    do i = 1, c%nspc
       aux = lbl1 // " " // string(c%spc(i)%z)
       lbl1 = aux
    end do
    write (lu,'("znucl ",A)') lbl1
    write (lu,'("natom ",I5)') c%ncel

    lbl1 = ""
    do i = 1, c%nspc
       ntyp = 0
       do j = 1, c%ncel
          if (c%atcel(j)%is == i) ntyp = ntyp + 1
       end do
       aux = lbl1 // " " // string(ntyp) // "*" // string(i)
       lbl1 = aux
    end do
    write (lu,'("typat ",A)') lbl1

    write (lu,'("xred ")')
    do i = 1, c%nspc
       do j = 1, c%ncel
          if (c%atcel(j)%is == i) then
             write (lu,'(" ",3(F15.10," "))') c%atcel(j)%x
          end if
       end do
    end do
    write (lu,*)

    ! template for the rest
    write (lu,'("#Definition of the planewave basis set")')
    write (lu,'("ecut 15")')
    write (lu,'("")')
    write (lu,'("# k-grid")')
    write (lu,'("kptopt 1")')
    write (lu,'("nshiftk 4")')
    write (lu,'("shiftk  0.5 0.5 0.5  ")')
    write (lu,'("        0.5 0.0 0.0")')
    write (lu,'("        0.0 0.5 0.0")')
    write (lu,'("        0.0 0.0 0.5")')
    write (lu,'("ngkpt   1 1 1")')
    write (lu,'("")')
    write (lu,'("#Definition of the SCF procedure")')
    write (lu,'("nstep  20      ")')
    write (lu,'("toldfe 1.0d-10")')
    write (lu,'("diemac 12.0   ")')
    call fclose(lu)

  end subroutine write_abinit

  !> Write an elk input template
  module subroutine write_elk(c,file,ti)
    use tools_io, only: fopen_write, fclose
    class(crystal), intent(in) :: c
    character*(*), intent(in) :: file
    type(thread_info), intent(in), optional :: ti

    integer :: ntyp
    integer :: i, j, lu

    ! Write input
    lu = fopen_write(file,ti=ti)
    write (lu,'("tasks"/,"0"/)')
    write (lu,'("xctype"/,"20"/)')
    write (lu,'("avec")')
    do i = 1, 3
       write (lu,'("  ",3(F15.10," "))') c%m_x2c(:,i)
    end do
    write (lu,*)

    write (lu,'("sppath"/,"''./''"/)')

    write (lu,'("atoms")')
    write (lu,'("  ",I4)') c%nspc
    do i = 1, c%nspc
       write (lu,'("  ","''",A,".in''")') trim(c%spc(i)%name)
       ntyp = 0
       do j = 1, c%ncel
          if (c%atcel(j)%is == i) ntyp = ntyp + 1
       end do
       write (lu,'("  ",I3)') ntyp
       do j = 1, c%ncel
          if (c%atcel(j)%is == i) then
             write (lu,'("  ",3(F14.10," "),"0.0 0.0 0.0")') c%atcel(j)%x
          end if
       end do
    end do
    write (lu,*)

    write (lu,'("ngridk"/,"  4 4 4"/)')
    write (lu,'("rgkmax"/,"  7.0"/)')
    write (lu,'("highq"/,"  .true."/)')
    call fclose(lu)

  end subroutine write_elk

  !> Write a Gaussian template input (periodic).
  module subroutine write_gaussian(c,file,ti)
    use tools_io, only: fopen_write, string, nameguess, ioj_left, fclose
    use param, only: bohrtoa
    class(crystal), intent(in) :: c
    character*(*), intent(in) :: file
    type(thread_info), intent(in), optional :: ti

    character(len=:), allocatable :: wroot
    integer :: lu, i, j

    wroot = file(:index(file,'.',.true.)-1)

    lu = fopen_write(file,ti=ti)
    write (lu,'("%chk=",A,".chk")') wroot
    write (lu,'("%nprocs=8")')
    write (lu,'("%mem=2GB")')
    write (lu,'("#p pbepbe/sto-3g int(grid=ultrafine) pop=regular fmm=(print)")')
    write (lu,'("   iop1=timestamp iop(5/13=1,5/33=1,5/181=10,5/184=186)")')
    write (lu,'("   scf=(novaracc,noincfock,tight)")')
    write (lu,*)
    write (lu,'("title")')
    write (lu,*)
    write (lu,'("0 1")')
    do i = 1, c%ncel
       write (lu,'(99(A," "))') string(nameguess(c%spc(c%atcel(i)%is)%z,.true.),2,ioj_left),&
          (string(c%atcel(i)%r(j)*bohrtoa,'f',14,8,ioj_left),j=1,3)
    end do
    do i = 1, 3
       write (lu,'(99(A," "))') string("Tv",2,ioj_left),&
          (string(c%m_x2c(j,i)*bohrtoa,'f',14,8,ioj_left),j=1,3)
    end do
    write (lu,*)

    call fclose(lu)

  end subroutine write_gaussian

  !> Write a tessel input template
  module subroutine write_tessel(c,file,ti)
    use global, only: fileroot
    use tools_io, only: fopen_write, fclose
    class(crystal), intent(in) :: c
    character*(*), intent(in) :: file
    type(thread_info), intent(in), optional :: ti

    integer :: lu, i

    lu = fopen_write(file,ti=ti)

    write (lu,'("set camangle 75 -10 45")')
    write (lu,'("set background background {color rgb <1,1,1>}")')
    write (lu,'("set use_planes .false.")')
    write (lu,'("set ball_texture finish{specular 0.2 roughness 0.1 reflection 0.1}")')
    write (lu,'("set equalscale noscale")')
    write (lu,'("molecule")')
    write (lu,'("  crystal")')
    write (lu,'("    symmatrix seitz")')
    do i = 1, c%ncv
       write (lu,'("     ",A,3(F15.12," "))') "cen ",c%cen(:,i)
    end do
    write (lu,'("     ","#")')
    do i = 1, c%neqv
       write (lu,'("     ",3(F5.2," "),F15.12)') c%rotm(1,:,i)
       write (lu,'("     ",3(F5.2," "),F15.12)') c%rotm(2,:,i)
       write (lu,'("     ",3(F5.2," "),F15.12)') c%rotm(3,:,i)
       write (lu,'("     ","#")')
    end do
    write (lu,'("     ","endsymmatrix")')
    write (lu,'("     ",A,6(F12.8," "))') "cell", c%aa, c%bb
    write (lu,'("     ","crystalbox  -2.30 -2.30 -2.30 2.30 2.30 2.30")')
    write (lu,'("     ",A,6(F6.3," "))') "clippingbox ",-0.02,-0.02,-0.02,+1.02,+1.02,+1.02
    do i = 1, c%nneq
       write (lu,'("     ","neq ",3(F12.8," "),A10)') c%at(i)%x, trim(c%spc(c%at(i)%is)%name)
    end do
    write (lu,'("  endcrystal")')
    write (lu,'(A)') "  unitcell radius 0.01 rgb 1.0 0.5 0.5 many"
    write (lu,'("  molmotif allmaincell jmol")')
    write (lu,'("  off ",A,".off")') trim(adjustl(fileroot))
    write (lu,'("  vrml ",A,".wrl")') trim(adjustl(fileroot))
    write (lu,'("  povray ",A,".pov")') trim(adjustl(fileroot))
    write (lu,'("endmolecule")')
    write (lu,'("# run povray -D -UV +I",A,".pov +O",A,".png +W2000 +H2000 +A")') &
       trim(adjustl(fileroot)), trim(adjustl(fileroot))
    write (lu,'("end")')

    call fclose(lu)

  end subroutine write_tessel

  !> Write a critic2 input template
  module subroutine write_critic(c,file,ti)
    use tools_io, only: fopen_write, fclose, string
    class(crystal), intent(in) :: c
    character*(*), intent(in) :: file
    type(thread_info), intent(in), optional :: ti

    integer :: lu, i, j

    lu = fopen_write(file,ti=ti)

    write (lu,'("crystal")')
    write (lu,'("  cell ",6(A," "))') (string(c%aa(i),'f',decimal=10),i=1,3),&
       (string(c%bb(i),'f',decimal=6),i=1,3)
    do i = 1, c%ncel
       write (lu,'("  neq ",3(A," "),A)') (string(c%atcel(i)%x(j),'f',decimal=10),j=1,3),&
          string(c%spc(c%atcel(i)%is)%name)
    end do
    write (lu,'("endcrystal")')
    write (lu,'("end")')
    call fclose(lu)

  end subroutine write_critic

  !> Write a simple cif file (filename = file) with the c crystal
  !> structure. If usesym0, write symmetry to the cif file; otherwise
  !> use P1.
  module subroutine write_cif(c,file,usesym0,ti)
    use global, only: fileroot, testing
    use tools_io, only: fopen_write, fclose, string, nameguess, deblank, nameguess,&
       ferror, faterr
    use param, only: bohrtoa, maxzat
    use tools_math, only: gcd
    class(crystal), intent(in) :: c
    character*(*), intent(in) :: file
    logical, intent(in) :: usesym0
    type(thread_info), intent(in), optional :: ti

    integer :: i, j, iz, lu, idx, gcdz
    character(len=mlen), allocatable :: strfin(:)
    character*2 :: sym
    character*3 :: schpg
    character(len=:), allocatable :: str
    integer :: holo, laue, natmol
    logical :: usesym, doz
    integer :: datvalues(8)
    integer, allocatable :: atc(:,:)

    ! Hill order for chemical formula. First C, then H, then all the other
    ! elements in alphabetical order.
    integer, parameter :: hillord(maxzat) = &
       (/6,1,89,47,13,95,18,33,85,79,5,56,4,107,83,97,35,20,48,58,98,17,&
       96,112,27,24,55,29,105,110,66,68,99,63,9,26,114,100,87,31,64,32,2,72,80,&
       67,108,53,49,77,19,36,57,3,103,71,116,115,101,12,25,42,109,7,11,41,60,10,&
       113,28,102,93,8,118,76,15,91,82,46,61,84,59,78,94,88,37,75,104,111,45,86,&
       44,16,51,21,34,106,14,62,50,38,73,65,43,52,90,22,81,69,117,92,23,74,54,&
       39,70,30,40/)

    ! use symmetry?
    usesym = usesym0 .and. c%spgavail

    ! open output file
    lu = fopen_write(file,ti=ti)

    ! header (date and time mucks up testing)
    write (lu,'("data_",A)') string(deblank(fileroot))
    write (lu,'("_audit_creation_method ''critic2''")')
    if (.not.testing) then
       call date_and_time(values=datvalues)
       write (lu,'("_audit_creation_date ",A,"-",A,"-",A)') &
          (string(datvalues(i)),i=1,3)
    end if

    ! formula: count the number of element types
    allocate(atc(maxzat,c%nmol))
    doz = c%ismol3d
    atc = 0
    do i = 1, c%nmol
       do j = 1, c%mol(i)%nat
          if (c%mol(i)%spc(c%mol(i)%at(j)%is)%z > 0 .and. c%mol(i)%spc(c%mol(i)%at(j)%is)%z <= maxzat) then
             atc(c%mol(i)%spc(c%mol(i)%at(j)%is)%z,i) = atc(c%mol(i)%spc(c%mol(i)%at(j)%is)%z,i) + 1
          end if
       end do
       if (i > 1) then
          if (any(atc(:,i) - atc(:,1) /= 0)) then
             doz = .false.
          end if
       end if
    end do

    ! formula
    if (.not.doz) then
       ! not a molecular crystal or different types of molecules,
       ! collect all atc then calculate gdc of all non-zero
       ! numbers. atc(:,1) is now the formula unit and gcdz = Z.
       do i = 2, c%nmol
          atc(:,1) = atc(:,1) + atc(:,i)
       end do
       gcdz = -1
       do i = 1, maxzat
          if (atc(i,1) > 0) then
             if (gcdz < 0) then
                gcdz = atc(i,1)
             else
                gcdz = gcd(gcdz,atc(i,1))
             end if
          end if
       end do
       atc(:,1) = atc(:,1) / gcdz
    else
       ! a molecular crystal with always the same molecule: gcdz = Z
       natmol = sum(atc(:,1))
       if (abs(real(c%ncel,8)/real(natmol,8) - c%ncel/natmol) > 1d-10) &
          call ferror('write_cif','inconsistent number of atoms in fragment',faterr)
       gcdz = c%ncel / natmol
    end if

    ! formula: build the molecular formula/formula unit and write to cif
    str = ""
    do i = 1, maxzat
       idx = hillord(i)
       if (atc(idx,1) > 0) then
          sym = nameguess(idx,.true.)
          if (atc(idx,1) > 1) then
             str = str // trim(sym) // string(atc(idx,1)) // " "
          else
             str = str // trim(sym) //  " "
          end if
       end if
    end do
    str = trim(str)
    write (lu,'("_chemical_formula_sum ''",A,"''")') str
    write (lu,'("_cell_formula_units_Z ",A)') string(gcdz)
    deallocate(atc)

    ! cell dimensions
    write (lu,'("_cell_length_a ",F20.10)') c%aa(1)*bohrtoa
    write (lu,'("_cell_length_b ",F20.10)') c%aa(2)*bohrtoa
    write (lu,'("_cell_length_c ",F20.10)') c%aa(3)*bohrtoa
    write (lu,'("_cell_angle_alpha ",F14.4)') c%bb(1)
    write (lu,'("_cell_angle_beta ",F14.4)') c%bb(2)
    write (lu,'("_cell_angle_gamma ",F14.4)') c%bb(3)
    write (lu,'("_cell_volume ",F20.6)') c%omega * bohrtoa**3

    if (usesym) then
       allocate(strfin(c%neqv*c%ncv))
       call c%struct_report_symxyz(strfin)
       do i = 1, c%neqv*c%ncv
          if (index(strfin(i),"not found") > 0) then
             usesym = .false.
             exit
          end if
       end do
    end if

    ! write the symmetry, if applicable
    if (usesym) then
       call pointgroup_info(c%spg%pointgroup_symbol,schpg,holo,laue)
       write (lu,'("_space_group_crystal_system ",A)') string(holo_string(holo))
       write (lu,'("_space_group_IT_number ",A)') string(c%spg%spacegroup_number)
       write (lu,'("_space_group_name_H-M_alt ''",A,"''")') string(c%spg%international_symbol)
       write (lu,'("_space_group_name_Hall ''",A,"''")') string(c%spg%hall_symbol)

       write (lu,'("loop_")')
       write (lu,'("_symmetry_equiv_pos_site_id")')
       write (lu,'("_symmetry_equiv_pos_as_xyz")')
       do i = 1, c%neqv*c%ncv
          write (lu,'(" ",A," ''",A,"''")') string(i), string(strfin(i))
       end do
    else
       write (lu,'("_space_group_crystal_system triclinic")')
       write (lu,'("_space_group_IT_number 1")')
       write (lu,'("_space_group_name_H-M_alt ''P 1''")')
       write (lu,'("_space_group_name_Hall ''P 1''")')

       write (lu,'("loop_")')
       write (lu,'("_symmetry_equiv_pos_site_id")')
       write (lu,'("_symmetry_equiv_pos_as_xyz")')
       write (lu,'(" 1 ''x,y,z''")')
    end if
    if (allocated(strfin)) deallocate(strfin)

    write (lu,'("loop_")')
    write (lu,'("_atom_site_label")')
    write (lu,'("_atom_site_type_symbol")')
    write (lu,'("_atom_site_fract_x")')
    write (lu,'("_atom_site_fract_y")')
    write (lu,'("_atom_site_fract_z")')
    if (usesym) then
       do i = 1, c%nneq
          iz = c%at(i)%is
          write (lu,'(A5," ",A3," ",3(F20.14," "))') c%spc(iz)%name, nameguess(c%spc(iz)%z,.true.), c%at(i)%x
       end do
    else
       do i = 1, c%ncel
          iz = c%atcel(i)%is
          write (lu,'(A5," ",A3," ",3(F20.14," "))') c%spc(iz)%name, nameguess(c%spc(iz)%z,.true.), c%atcel(i)%x
       end do
    end if
    call fclose(lu)

  end subroutine write_cif

  !> Write a simple d12 file
  module subroutine write_d12(c,file,dosym,ti)
    use tools_io, only: fopen_write, fclose, string
    use param, only: bohrtoa
    class(crystal), intent(in) :: c
    character*(*), intent(in) :: file
    logical, intent(in) :: dosym
    type(thread_info), intent(in), optional :: ti

    character(len=:), allocatable :: file34
    character(len=3) :: schpg
    integer :: lu, holo, laue
    integer :: i, j, k, l
    real*8 :: x(3)
    real*8 :: dum(3,3)

    lu = fopen_write(file,ti=ti)
    write (lu,'("Title")')
    if (c%ismolecule) then
       write (lu,'("MOLECULE")')
       write (lu,'("1")')
       write (lu,'(A)') string(c%ncel)
       do i = 1, c%ncel
          write (lu,'(4(A," "))') string(c%spc(c%atcel(i)%is)%z), &
             (string((c%atcel(i)%r(j)+c%molx0(j))*bohrtoa,'f',15,8),j=1,3)
       end do
    elseif (dosym) then
       write (lu,'("EXTERNAL")')
    else
       write (lu,'("CRYSTAL")')
       write (lu,'("0 0 0")')
       write (lu,'("1")')
       write (lu,'(6(A," "))') (string(c%aa(i)*bohrtoa,'f',15,8),i=1,3), (string(c%bb(j),'f',15,8),j=1,3)
       write (lu,'(A)') string(c%ncel)
       do i = 1, c%ncel
          write (lu,'(4(A," "))') string(c%spc(c%atcel(i)%is)%z), (string(c%atcel(i)%x(j),'f',15,8),j=1,3)
       end do
    end if
    write (lu,'("TESTGEOM")')
    write (lu,'("END")')
    write (lu,'("END")')
    write (lu,'("END")')
    call fclose(lu)

    if (.not.c%ismolecule.and.dosym) then
       file34 = file(:index(file,'.',.true.)-1) // ".fort.34"
       lu = fopen_write(file34,ti=ti)

       ! for a crystal, if symmetry is available
       ! header: dimensionality, centring type, and crystal holohedry
       call pointgroup_info(c%spg%pointgroup_symbol,schpg,holo,laue)
       write (lu,'("3 1 ",A)') string(holo-1)
       do i = 1, 3
          write (lu,'(3(A," "))') (string(c%m_x2c(j,i)*bohrtoa,'f',decimal=10),j=1,3)
       end do

       ! symmetry operations
       write (lu,'(A)') string(c%neqv*c%ncv)
       do i = 1, c%neqv
          do j = 1, c%ncv
             dum = transpose(matmul(matmul(c%m_x2c,c%rotm(1:3,1:3,i)),c%m_c2x))
             do k = 1, 3
                write (lu,'(3(" ",E19.12))') (dum(l,k),l=1,3)
             end do
             x = c%rotm(:,4,i)+c%cen(:,j)
             x = c%x2c(x)
             write (lu,'(3(" ",E19.12))') (x(l)*bohrtoa,l=1,3)
          end do
       end do

       ! atoms
       write (lu,'(A)') string(c%nneq)
       do i = 1, c%nneq
          write (lu,'(4(A," "))') string(c%spc(c%at(i)%is)%z), (string(c%at(i)%r(l)*bohrtoa,'f',decimal=10),l=1,3)
       end do

       call fclose(lu)
    end if

  end subroutine write_d12

  !> Write a shelx res file (filename = file) with the c crystal
  !> structure. If usesym0, write symmetry to the cif file; otherwise
  !> use P1. dosym = 0 (do not use symmetry), 1 (use symmetry), or
  !> -1 (use symmetry only if possible, do not emit warnings)
  module subroutine write_res(c,file,dosym,ti)
    use tools_io, only: fopen_write, fclose, string, ferror, warning, nameguess
    use tools_math, only: det3
    use param, only: bohrtoa, eye
    class(crystal), intent(in) :: c
    character*(*), intent(in) :: file
    integer, intent(in) :: dosym
    type(thread_info), intent(in), optional :: ti

    integer :: i, j, lu, ilatt
    character(len=mlen), allocatable :: strfin(:)
    character(len=:), allocatable :: str
    logical :: usesym, ok3(3), ok2(2)
    real*8 :: dd

    real*8, parameter :: cen_i(3)  = (/0.5d0,0.5d0,0.5d0/)
    real*8, parameter :: cen_a(3)  = (/0.0d0,0.5d0,0.5d0/)
    real*8, parameter :: cen_b(3)  = (/0.5d0,0.0d0,0.5d0/)
    real*8, parameter :: cen_c(3)  = (/0.5d0,0.5d0,0.0d0/)
    real*8, parameter :: cen_r1(3) = (/2d0,1d0,1d0/) / 3d0
    real*8, parameter :: cen_r2(3) = (/1d0,2d0,2d0/) / 3d0
    real*8, parameter :: eps = 1d-5

    ! use symmetry?
    usesym = (dosym==1 .or. dosym==-1) .and. c%spgavail

10  continue

    ! open output file
    lu = fopen_write(file,ti=ti)

    ! header
    write (lu,'("TITL critic2 | ",A)') trim(c%file)
    write (lu,'("CELL 0.71073 ",6(A," "))') (string(c%aa(i)*bohrtoa,'f',12,8),i=1,3), &
       (string(c%bb(j),'f',10,6),j=1,3)
    if (usesym) then
       write (lu,'("ZERR ",A," 0.0001 0.0001 0.0001 0.0001 0.0001 0.0001")') string(c%ncel/c%nneq)
    else
       write (lu,'("ZERR 1 0.0001 0.0001 0.0001 0.0001 0.0001 0.0001")')
    end if
    if (usesym) then
       ! identify lattice type
       ilatt = 0
       if (c%ncv == 1) then
          ilatt = -1 ! P
       elseif (c%ncv == 2) then
          if (all(abs(c%cen(:,2) - cen_i - nint(c%cen(:,2) - cen_i)) < eps)) then
             ilatt = -2 ! I
          elseif (all(abs(c%cen(:,2) - cen_a - nint(c%cen(:,2) - cen_a)) < eps)) then
             ilatt = -5 ! A
          elseif (all(abs(c%cen(:,2) - cen_b - nint(c%cen(:,2) - cen_b)) < eps)) then
             ilatt = -6 ! B
          elseif (all(abs(c%cen(:,2) - cen_c - nint(c%cen(:,2) - cen_c)) < eps)) then
             ilatt = -7 ! C
          end if
       elseif (c%ncv == 3) then
          ok2 = .false.
          do i = 2, c%ncv
             if (.not.ok2(1)) ok2(1) = all(abs(c%cen(:,i) - cen_r1 - nint(c%cen(:,i) - cen_r1)) < eps)
             if (.not.ok2(2)) ok2(2) = all(abs(c%cen(:,i) - cen_r2 - nint(c%cen(:,i) - cen_r2)) < eps)
          end do
          if (all(ok2)) ilatt = -3
       elseif (c%ncv == 4) then
          ok3 = .false.
          do i = 2, c%ncv
             if (.not.ok3(1)) ok3(1) = all(abs(c%cen(:,i) - cen_a - nint(c%cen(:,i) - cen_a)) < eps)
             if (.not.ok3(2)) ok3(2) = all(abs(c%cen(:,i) - cen_b - nint(c%cen(:,i) - cen_b)) < eps)
             if (.not.ok3(3)) ok3(3) = all(abs(c%cen(:,i) - cen_c - nint(c%cen(:,i) - cen_c)) < eps)
          end do
          if (all(ok3)) ilatt = -4
       end if
       if (ilatt == 0) then
          if (dosym == 1) &
             call ferror('write_res','unknown set of centering vectors',warning)
          usesym = .false.
          call fclose(lu)
          goto 10
       end if

       ! identify centrosymmetry
       do i = 1, c%neqv
          if (all(abs(c%rotm(1:3,1:3,i) + eye) < eps) .and. all(abs(c%rotm(1:3,4,i)) < eps)) then
             ilatt = -ilatt
             exit
          end if
       end do
    else
       ilatt = -1
    end if
    write (lu,'("LATT ",A)') string(ilatt)

    if (usesym) then
       allocate(strfin(c%neqv*c%ncv))
       call c%struct_report_symxyz(strfin)
       do i = 2, c%neqv ! skip the identity
          dd = det3(c%rotm(1:3,1:3,i))
          if (dd > 0d0 .or. ilatt < 0) then
             if (index(strfin(i),"not found") > 0) then
                if (dosym == 1) &
                   call ferror('write_res','unknown set of centering vectors',warning)
                usesym = .false.
                call fclose(lu)
                goto 10
             end if
             write (lu,'("SYMM ",A)') trim(strfin(i))
          end if
       end do
    end if

    ! atomic species
    write (lu,'("SFAC ",999(A," "))') (trim(nameguess(c%spc(i)%z,.true.)),i=1,c%nspc)

    ! number of atoms of each type
    str = ""
    do i = 1, c%nspc
       str = str // " " // string(count(c%atcel(:)%is == i))
    end do
    write (lu,'("UNIT ",A)') str
    write (lu,'("FVAR 1.00")')

    ! list of atoms
    if (usesym) then
       do i = 1, c%nneq
          write (lu,'(999(A," "))') trim(c%spc(c%at(i)%is)%name) // string(i), string(c%at(i)%is), &
             (string(c%at(i)%x(j),'f',12,8),j=1,3), string(real(c%at(i)%mult,8)/(c%neqv*c%ncv),'f',12,8), &
             "0.05"
       end do
    else
       do i = 1, c%ncel
          write (lu,'(999(A," "))') trim(c%spc(c%atcel(i)%is)%name) // string(i), string(c%atcel(i)%is), &
             (string(c%atcel(i)%x(j),'f',12,8),j=1,3), "1.0", "0.05"
       end do
    end if

    ! close the file
    write (lu,'("END")')
    call fclose(lu)

  end subroutine write_res

  !> Write an escher octave script
  module subroutine write_escher(c,file,ti)
    use global, only: fileroot
    use tools_io, only: fopen_write, string, fclose
    use param, only: pi
    class(crystal), intent(in) :: c
    character*(*), intent(in) :: file
    type(thread_info), intent(in), optional :: ti

    character(len=:), allocatable :: lbl1, aux
    integer :: lu, i, n

    lu = fopen_write(file,ti=ti)

    ! count number of atoms per type
    write (lu,'("cr = struct();")')
    write (lu,'("cr.name = """,A,""";")') trim(adjustl(fileroot))
    write (lu,'("cr.a = [",1p,3(E22.14," "),"];")') c%aa
    write (lu,'("cr.b = [",1p,3(E22.14," "),"];")') c%bb * pi / 180d0
    write (lu,'("cr.nat = ",I6,";")') c%ncel
    write (lu,'("cr.ntyp = ",I6,";")') c%nspc
    write (lu,'("cr.r = [")')
    do i = 1, 3
       write (lu,'("  ",1p,3(E22.14," "))') c%m_x2c(:,i)
    end do
    write (lu,'("  ","];")')
    write (lu,'("cr.g = [")')
    do i = 1, 3
       write (lu,'("  ",1p,3(E22.14," "))') c%gtensor(:,i)
    end do
    write (lu,'("  ","];")')
    write (lu,'("cr.omega = ",1p,E22.14,";")') c%omega

    lbl1 = "cr.ztyp = ["
    do i = 1, c%nspc
       aux = lbl1 // " " // string(c%spc(i)%z)
       lbl1 = aux
    end do
    aux = lbl1 // "];"
    lbl1 = aux
    write (lu,'(A)') lbl1

    lbl1 = "cr.attyp = {"
    do i = 1, c%nspc
       if (i > 1) then
          aux = lbl1 // ","
          lbl1 = aux
       end if
       aux = lbl1 // '"' // string(c%spc(i)%name) // '"'
       lbl1 = aux
    end do
    aux = lbl1 // "};"
    lbl1 = aux
    write (lu,'(A)') lbl1

    lbl1 = "cr.typ = ["
    do i = 1, c%ncel
       aux = lbl1 // " " // string(c%atcel(i)%is)
       lbl1 = aux
    end do
    aux = lbl1 // "];"
    lbl1 = aux
    write (lu,'(A)') lbl1

    write (lu,'("cr.x = [")')
    n = 0
    do i = 1, c%ncel
       write (lu,'("  ",1p,3(E22.14," "))') c%atcel(i)%x
    end do
    write (lu,'("  ];")')

    call fclose(lu)

  end subroutine write_escher

  !> Write a db file for the dcp automatic input generator
  module subroutine write_db(c,file,ti)
    use tools_io, only: fopen_write, string, fclose, nameguess
    use param, only: bohrtoa
    class(crystal), intent(in) :: c
    character*(*), intent(in) :: file
    type(thread_info), intent(in), optional :: ti

    integer :: lu, i, j

    lu = fopen_write(file,ti=ti)
    write (lu,'("type crystal_energy")')
    write (lu,'("kpts 4")')
    write (lu,'("crys")')
    write (lu,'(6(A," "))') (string(c%aa(i)*bohrtoa,'f',18,10),i=1,3), (string(c%bb(i),'f',18,10),i=1,3)
    do i = 1, c%ncel
       write (lu,'(A," ",1p,3(A," "))') string(nameguess(c%spc(c%atcel(i)%is)%z,.true.)), &
          (string(c%atcel(i)%x(j),'f',18,10),j=1,3)
    end do
    write (lu,'("end")')
    call fclose(lu)

  end subroutine write_db

  !> Write a gulp input script
  module subroutine write_gulp(c,file,ti)
    use tools_io, only: fopen_write, nameguess, fclose, string
    use param, only: bohrtoa
    class(crystal), intent(inout) :: c
    character*(*), intent(in) :: file
    type(thread_info), intent(in), optional :: ti

    integer :: lu, i, j

    lu = fopen_write(file,ti=ti)
    write (lu,'("eem")')
    write (lu,'("cell ",6(A," "))') (string(c%aa(j) * bohrtoa,'f',13,9),j=1,3), &
       (string(c%bb(j),'f',10,5),j=1,3)
    write (lu,'("fractional")')
    do i = 1, c%ncel
       write (lu,'(A5," ",3(A," "))') trim(c%spc(c%atcel(i)%is)%name),&
          (string(c%atcel(i)%x(j),'f',15,9),j=1,3)
    end do

    call fclose(lu)

  end subroutine write_gulp

  !> Write a lammps data file
  module subroutine write_lammps(c,file,ti)
    use tools_io, only: fopen_write, ferror, faterr, fclose
    use tools_math, only: m_x2c_from_cellpar
    use param, only: bohrtoa, atmass
    class(crystal), intent(in) :: c
    character*(*), intent(in) :: file
    type(thread_info), intent(in), optional :: ti

    integer :: i, j, l
    integer :: lu
    real*8 :: rnew(3,3)

    lu = fopen_write(file,ti=ti)

    ! header
    write (lu,'("LAMMPS data file created by critic2. (experimental)",/)')
    write (lu,'(I9," atoms")') c%ncel
    write (lu,'(I9," atom types")') c%nspc
    write (lu,*)

    ! metrics of the cell --> this needs more testing
    rnew = m_x2c_from_cellpar(c%aa,c%bb)
    if (abs(c%m_x2c(1,2)) > 1d-12 .or. abs(c%m_x2c(1,3)) > 1d-12 .or.&
       abs(c%m_x2c(2,3)) > 1d-12) then
       call ferror ('write_lammps','non-orthogonal cells not implemented',faterr)
    end if
    write (lu,'(2(F18.10," ")," xlo xhi")') 0d0, c%m_x2c(1,1)*bohrtoa
    write (lu,'(2(F18.10," ")," ylo yhi")') 0d0, c%m_x2c(2,2)*bohrtoa
    write (lu,'(2(F18.10," ")," zlo zhi")') 0d0, c%m_x2c(3,3)*bohrtoa
    write (lu,'(3(F18.10," ")," xy xz yz")') 0d0, 0d0, 0d0
    write (lu,*)

    write (lu,'("Masses"/)')
    do i = 1, c%nspc
       write (lu,'(I3," ",F10.4)') i, atmass(c%spc(i)%z)
    end do
    write (lu,*)

    write (lu,'("Atoms"/)')
    l = 0
    do i = 1, c%nspc
       do j = 1, c%ncel
          if (c%atcel(j)%is /= i) cycle
          l = l + 1
          write (lu,'(I7," ",I3," ",F4.1,3(F15.8," "))') l, i, 0d0, c%atcel(j)%r*bohrtoa
       end do
    end do

    call fclose(lu)

  end subroutine write_lammps

  !> Write a siesta fdf data file
  module subroutine write_siesta_fdf(c,file,ti)
    use tools_io, only: fopen_write, nameguess, lower, fclose
    use param, only: bohrtoa
    class(crystal), intent(in) :: c
    character*(*), intent(in) :: file
    type(thread_info), intent(in), optional :: ti

    integer :: i, j
    integer :: lu

    lu = fopen_write(file,ti=ti)

    ! header
    write (lu,'("# fdf file created by critic2.",/)')
    write (lu,'("SystemName crystal")')
    write (lu,'("SystemLabel crystal")')
    write (lu,*)

    write (lu,'("NumberOfSpecies ",I3)') c%nspc
    write (lu,'("NumberOfAtoms ", I6)') c%ncel
    write (lu,'("%block Chemical_Species_Label")')
    do i = 1, c%nspc
       write (lu,'(I3,I3," ",A2)') i, c%spc(i)%z, lower(nameguess(c%spc(i)%z,.true.))
    end do
    write (lu,'("%endblock Chemical_Species_Label")')
    write (lu,*)

    write (lu,'("LatticeConstant 1.0 ang")')
    write (lu,'("%block LatticeParameters")')
    write (lu,'(3(F16.10," "),3(F16.8," "))') c%aa*bohrtoa, c%bb
    write (lu,'("%endblock LatticeParameters")')
    write (lu,'("AtomicCoordinatesFormat Fractional")')
    write (lu,'("%block AtomicCoordinatesAndAtomicSpecies")')
    do i = 1, c%nspc
       do j = 1, c%ncel
          if (c%atcel(j)%is == i) &
             write (lu,'(3(F18.12," "),I3)') c%atcel(j)%x, i
       end do
    end do
    write (lu,'("%endblock AtomicCoordinatesAndAtomicSpecies")')
    write (lu,*)

    write (lu,'("XC.functional GGA")')
    write (lu,'("XC.authors PBE")')
    write (lu,'("SpinPolarized .false.")')
    write (lu,'("MaxSCFIterations 100")')
    write (lu,'("MeshCutoff 100. Ry")')
    write (lu,'("DM.NumberPulay 3")')
    write (lu,*)

    write (lu,'("PAO.BasisSize DZP")')
    write (lu,*)

    write (lu,'("kgrid_cutoff 10.0 ang")')
    write (lu,*)

    write (lu,'("ElectronicTemperature 5 K")')
    write (lu,*)

    write (lu,'("# options")')
    write (lu,'("LongOutput")')
    write (lu,'("SaveRho")')
    write (lu,'("SaveBaderCharge")')
    write (lu,'("DM.UseSaveDM")')
    write (lu,'("WriteDenchar")')
    write (lu,'("WriteCoorXmol")')

    call fclose(lu)

  end subroutine write_siesta_fdf

  !> Write a siesta STRUCT_IN data file
  module subroutine write_siesta_in(c,file,ti)
    use tools_io, only: fopen_write, uout, nameguess, string, fclose
    use param, only: bohrtoa
    class(crystal), intent(in) :: c
    character*(*), intent(in) :: file
    type(thread_info), intent(in), optional :: ti

    integer :: lu
    real*8 :: r(3,3)
    integer :: i, j, k

    lu = fopen_write(file,ti=ti)

    ! lattice vectors
    r = transpose(c%m_x2c) * bohrtoa
    do i = 1, 3
       write (lu,'(3(F20.12," "))') r(i,:)
    end do

    ! atoms
    write (lu,*) c%ncel
    j = 0
    do i = 1, c%nspc
       do k = 1, c%ncel
          if (c%atcel(k)%is == i) &
             write (lu,'(I3," ",I3," ",3(F20.12," "))') i, c%spc(i)%z, c%atcel(k)%x
       end do
    end do

    call fclose(lu)

    ! Write the chemical species block to standard output
    write (uout,'("%block Chemical_Species_Label")')
    do i = 1, c%nspc
       write (uout,'(3("  ",A))') string(i), string(c%spc(i)%z), &
          string(nameguess(c%spc(i)%z,.true.))
    end do
    write (uout,'("%endblock Chemical_Species_Label")')
    write (uout,*)

  end subroutine write_siesta_in

  !> Write a DFTB+ human-friendly structured data format (hsd) file
  module subroutine write_dftbp_hsd(c,file,ti)
    use tools_io, only: fopen_write, string, nameguess, fclose
    use param, only: maxzat0
    class(crystal), intent(in) :: c
    character*(*), intent(in) :: file
    type(thread_info), intent(in), optional :: ti

    real*8, parameter :: hderiv(maxzat0) = (/&
     -0.1857d0,      0.d0,      0.d0,    0.d0,      0.d0, -0.1492d0,& ! 1:6   (H-C)
     -0.1535d0, -0.1575d0, -0.1623d0,    0.d0, -0.0454d0,   -0.02d0,& ! 7:12  (N-Mg)
          0.d0,      0.d0,   -0.14d0, -0.11d0, -0.0697d0,      0.d0,& ! 13:18 (Al-Ar)
     -0.0339d0, -0.0340d0,      0.d0,    0.d0,      0.d0,      0.d0,& ! 19:24 (K-Cr)
          0.d0,      0.d0,      0.d0,    0.d0,      0.d0,   -0.03d0,& ! 25:30 (Mn-Zn)
          0.d0,      0.d0,      0.d0,    0.d0, -0.0573d0,      0.d0,& ! 31:36 (Ga-Kr)
          0.d0,      0.d0,      0.d0,    0.d0,      0.d0,      0.d0,& ! 37:42 (Rb-Mo)
          0.d0,      0.d0,      0.d0,    0.d0,      0.d0,      0.d0,& ! 43:48 (Tc-Cd)
          0.d0,      0.d0,      0.d0,    0.d0, -0.0433d0,      0.d0,& ! 49:54 (In-Xe)
          0.d0,      0.d0,      0.d0,    0.d0,      0.d0,      0.d0,& ! 55:60 (Cs-Nd)
          0.d0,      0.d0,      0.d0,    0.d0,      0.d0,      0.d0,& ! 61:66 (Pm-Dy)
          0.d0,      0.d0,      0.d0,    0.d0,      0.d0,      0.d0,& ! 67:72 (Ho-Hf)
          0.d0,      0.d0,      0.d0,    0.d0,      0.d0,      0.d0,& ! 73:78 (Ta-Pt)
          0.d0,      0.d0,      0.d0,    0.d0,      0.d0,      0.d0,& ! 79:84 (Au-Po)
          0.d0,      0.d0,      0.d0,    0.d0,      0.d0,      0.d0,& ! 85:90 (At-Th)
          0.d0,      0.d0,      0.d0,    0.d0,      0.d0,      0.d0,& ! 91:96 (Pa-Cm)
          0.d0,      0.d0,      0.d0,    0.d0,      0.d0,      0.d0,& ! 97:102 (Bk-No)
          0.d0,      0.d0,      0.d0,    0.d0,      0.d0,      0.d0,& ! 103:108 (Lr-Hs)
          0.d0,      0.d0,      0.d0,    0.d0,      0.d0,      0.d0,& ! 109:114 (Mt-Uuq)
          0.d0,      0.d0,      0.d0,    0.d0,&                       ! 115:118 (Uup-Uuo)
          0.d0,      0.d0,      0.d0,    0.d0,      0.d0/)            ! 119:123

    character*1, parameter :: maxang(maxzat0) = (/&
       "s", "x", "x", "x", "x", "p",& ! 1:6   (H-C)
       "p", "p", "p", "x", "p", "p",& ! 7:12  (N-Mg)
       "x", "x", "d", "d", "d", "x",& ! 13:18 (Al-Ar)
       "p", "p", "x", "x", "x", "x",& ! 19:24 (K-Cr)
       "x", "x", "x", "x", "x", "d",& ! 25:30 (Mn-Zn)
       "x", "x", "x", "x", "d", "x",& ! 31:36 (Ga-Kr)
       "x", "x", "x", "x", "x", "x",& ! 37:42 (Rb-Mo)
       "x", "x", "x", "x", "x", "x",& ! 43:48 (Tc-Cd)
       "x", "x", "x", "x", "d", "x",& ! 49:54 (In-Xe)
       "x", "x", "x", "x", "x", "x",& ! 55:60 (Cs-Nd)
       "x", "x", "x", "x", "x", "x",& ! 61:66 (Pm-Dy)
       "x", "x", "x", "x", "x", "x",& ! 67:72 (Ho-Hf)
       "x", "x", "x", "x", "x", "x",& ! 73:78 (Ta-Pt)
       "x", "x", "x", "x", "x", "x",& ! 79:84 (Au-Po)
       "x", "x", "x", "x", "x", "x",& ! 85:90 (At-Th)
       "x", "x", "x", "x", "x", "x",& ! 91:96 (Pa-Cm)
       "x", "x", "x", "x", "x", "x",& ! 97:102 (Bk-No)
       "x", "x", "x", "x", "x", "x",& ! 103:108 (Lr-Hs)
       "x", "x", "x", "x", "x", "x",& ! 109:114 (Mt-Uuq)
       "x", "x", "x", "x",&           ! 115:118 (Uup-Uuo)
       "x", "x", "x", "x", "x"/)      ! 119:123

    integer :: lu, i

    lu = fopen_write(file,ti=ti)
    write (lu,'("Geometry = GenFormat {")')
    call c%write_dftbp_gen(file,lu,ti=ti)
    write(lu,'("}")')
    write(lu,'("")')
    write(lu,'("Driver = ConjugateGradient {")')
    write(lu,'("       MovedAtoms = 1:-1")')
    write(lu,'("       MaxForceComponent = 1e-5")')
    write(lu,'("       MaxSteps = 3000")')
    write(lu,'("       LatticeOpt = Yes")')
    write(lu,'("       OutputPrefix = ""geo_end""")')
    write(lu,'("}")')
    write(lu,'("")')
    write(lu,'("Hamiltonian = DFTB{")')
    write(lu,'("  ThirdOrderFull = Yes")')
    write(lu,'("  SCC = Yes")')
    write(lu,'("  SCCTolerance = 1e-7")')
    write(lu,'("  MaxSCCIterations = 125")')
    write(lu,'("  MaxAngularMomentum = {")')
    do i = 1, c%nspc
       write (lu,'("    ",A," = ",A)') string(nameguess(c%spc(i)%z,.true.)), &
          string(maxang(c%spc(i)%z))
    end do
    write(lu,'("  }")')
    write(lu,'("  SlaterKosterFiles = Type2FileNames {")')
    write(lu,'("    Prefix = ""xxx""")')
    write(lu,'("    Separator = ""-""")')
    write(lu,'("    Suffix = "".skf""")')
    write(lu,'("    LowerCaseTypeName = No")')
    write(lu,'("  }")')
    if (.not.c%ismolecule) then
       write(lu,'("  KPointsAndWeights = SupercellFolding {")')
       write(lu,'("    4 0 0 ")')
       write(lu,'("    0 4 0")')
       write(lu,'("    0 0 4")')
       write(lu,'("    0.5 0.5 0.5")')
       write(lu,'("  }")')
    end if
    write(lu,'("  DampXH = Yes")')
    write(lu,'("  DampXHExponent = 4.2")')
    write(lu,'("  HubbardDerivs {")')
    do i = 1, c%nspc
       write (lu,'("    ",A," = ",A)') string(nameguess(c%spc(i)%z,.true.)), &
          string(hderiv(c%spc(i)%z),'f',decimal=4)
    end do
    write(lu,'("  }")')
    write(lu,'("}")')
    write(lu,'("")')
    write(lu,'("Options {")')
    write(lu,'("  WriteDetailedXML = Yes")')
    write(lu,'("}")')
    write(lu,'("")')
    write(lu,'("ParserOptions {")')
    write(lu,'("  ParserVersion = 4")')
    write(lu,'("}")')
    write(lu,'("")')
    call fclose(lu)

  end subroutine write_dftbp_hsd

  !> Write the crystal structure in pyscf format (python script)
  module subroutine write_pyscf(c,file,ti)
    use tools_io, only: fopen_write, fclose, string, nameguess
    class(crystal), intent(in) :: c
    character*(*), intent(in) :: file
    type(thread_info), intent(in), optional :: ti

    integer :: lu, i, j

    lu = fopen_write(file,ti=ti)
    if (c%ismolecule) then
       write (lu,'("from pyscf import gto")')
       write (lu,'("")')
       write (lu,'("mol = gto.Mole()")')
       write (lu,'("mol.atom = ''''''")')
       do i = 1, c%ncel
          write (lu,'(A,3(" ",A))') string(nameguess(c%spc(c%atcel(i)%is)%z,.true.)),&
             (string(c%atcel(i)%r(j),'f',18,10),j=1,3)
       end do
       write (lu,'("''''''")')
       write (lu,'("mol.basis = ''aug-cc-pvtz''")')
       write (lu,'("mol.verbose = 4")')
       write (lu,'("mol.spin = 0")')
       write (lu,'("mol.charge = 0")')
       write (lu,'("mol.build()")')
    else
       write (lu,'("from pyscf.pbc import gto")')
       write (lu,'("")')
       write (lu,'("cell = gto.Cell()")')
       write (lu,'("cell.atom = ''''''")')
       do i = 1, c%ncel
          write (lu,'(A,3(" ",A))') string(nameguess(c%spc(c%atcel(i)%is)%z,.true.)),&
             (string(c%atcel(i)%r(j),'f',18,10),j=1,3)
       end do
       write (lu,'("''''''")')
       write (lu,'("cell.a = ''''''")')
       do i = 1, 3
          write (lu,'("  ",1p,3(E22.14," "))') c%m_x2c(:,i)
       end do
       write (lu,'("''''''")')

       write (lu,'("cell.unit = ''Bohr''")')
       write (lu,'("cell.basis = ''gth-szv''")')
       write (lu,'("cell.pseudo = ''gth-pade''")')
       write (lu,'("cell.verbose = 4")')
       write (lu,'("cell.build()")')
    end if
    call fclose(lu)

  end subroutine write_pyscf

  !> Write the crystal or molecualr structure in FHIaims geometry.in
  !> format.  If frac = .true., use atom_frac instead of atom if c is
  !> a crystal.
  module subroutine write_fhi(c,file,frac,rklength,ti)
    use tools_io, only: fopen_write, fclose, string, nameguess
    use param, only: bohrtoa
    class(crystal), intent(in) :: c
    character*(*), intent(in) :: file
    logical, intent(in) :: frac
    real*8, intent(in), optional :: rklength
    type(thread_info), intent(in), optional :: ti

    integer :: lu, i, j, nk(3)
    character*2 :: name

    nk = 0
    if (present(rklength)) &
       call c%get_kpoints(rklength,nk)

    lu = fopen_write(file,ti=ti)
    write (lu,'("## FHIaims input file generated by critic2.")')
    if (.not.c%ismolecule) then
       do i = 1, 3
          write (lu,'("lattice_vector ",3(A," "))') (string(c%m_x2c(j,i) * bohrtoa,'f',18,10),j=1,3)
       end do
    end if
    do i = 1, c%ncel
       name = nameguess(c%spc(c%atcel(i)%is)%z,.true.)
       if (frac .and. .not.c%ismolecule) then
          write (lu,'("atom_frac ",4(" ",A))') (string(c%atcel(i)%x(j),'f',18,10),j=1,3),&
             name
       else
          write (lu,'("atom ",4(" ",A))') (string(c%atcel(i)%r(j)*bohrtoa,'f',18,10),j=1,3),&
             name
       end if
    end do
    call fclose(lu)

    if (present(rklength) .and. all(nk > 0)) then
       lu = fopen_write(trim(file) // "_control",ti=ti)
       write (lu,'("xc b86bpbe")')
       write (lu,'("spin none")')
       write (lu,'("charge 0")')
       write (lu,'("output_level MD_light")')
       write (lu,'("xdm")')
       write (lu,'("relativistic atomic_zora scalar")')
       write (lu,'("")')
       write (lu,'("sc_accuracy_rho 1e-7")')
       write (lu,'("")')
       write (lu,'("relax_geometry bfgs 0.005")')
       write (lu,'("relax_unit_cell full")')
       write (lu,'("max_relaxation_steps 200")')
       write (lu,'("")')
       write (lu,'("k_grid",3(" ",A))') (string(nk(i)),i=1,3)
       call fclose(lu)
    end if

  end subroutine write_fhi

  !> Write a TINKER frac file with tiny FF parametrization.
  module subroutine write_tinkerfrac(c,file,ti)
    use tools_io, only: fopen_write, string, fclose, nameguess, ioj_right, ioj_left, ferror, faterr,&
       warning
    use tools, only: tiny_atom_type
    use types, only: realloc
    use param, only: bohrtoa
    use fragmentmod, only: fragment
    class(crystal), intent(in) :: c
    character*(*), intent(in) :: file
    type(thread_info), intent(in), optional :: ti

    integer :: i, j, nn, lu
    integer :: iz, nneig, ityp
    type(fragment) :: fr
    integer, allocatable :: imap(:), ineig(:)
    logical :: dowarn

    ! open file
    lu = fopen_write(file,ti=ti)

    ! line 1: number of atoms and title (the spaces at the end are there for the USPEX reader)
    write (lu,'(A,"          ")') string(c%ncel)

    ! line 2: cell parameters (the spaces at the end are there for the USPEX reader)
    write (lu,'(6(A," "),"          ")') (string(c%aa(i)*bohrtoa,'f',decimal=8),i=1,3), &
       (string(c%bb(i),'f',decimal=8),i=1,3)

    ! merge all molecules into the onemotif
    call fr%merge_array(c%mol(1:c%nmol),.false.)

    ! build the atom map
    allocate(imap(c%ncel))
    imap = 0
    nn = 0
    do i = 1, fr%nat
       imap(fr%at(i)%cidx) = i
    end do
    if (any(imap == 0)) &
       call ferror ('write_tinkerfrac',&
       'some atoms are not mapped into molecular fragments; perhaps not a molecular crystal?',faterr)

    ! do the atom typing and write the atom list
    dowarn = .false.
    allocate(ineig(10))
    do i = 1, fr%nat
       ! type of atom
       iz = fr%spc(fr%at(i)%is)%z
       nneig = c%nstar(fr%at(i)%cidx)%ncon
       ityp = tiny_atom_type(iz,nneig)
       dowarn = dowarn .or. (ityp == 0)

       ! neighbor list
       if (nneig > size(ineig,1)) call realloc(ineig,nneig)
       do j = 1, nneig
          ineig(j) = imap(c%nstar(fr%at(i)%cidx)%idcon(j))
       end do

       write (lu,'(999(A," "))') string(i,5,ioj_right), &
          string(nameguess(iz,.true.),5,ioj_left), &
          (string(fr%at(i)%x(j),'f',12,8,ioj_right),j=1,3), string(ityp),&
          (string(ineig(j)),j=1,nneig)
    end do
    deallocate(imap,ineig)

    ! close and clean up
    if (dowarn) &
       call ferror('write_tinkerfrac','some atom types could not be assigned',warning)
    call fclose(lu)

  end subroutine write_tinkerfrac

  !> Write a DFTB+ human-friendly gen structure file
  module subroutine write_dftbp_gen(c,file,lu0,ti)
    use tools_io, only: fopen_write, nameguess, string, fclose
    use param, only: bohrtoa
    class(crystal), intent(in) :: c
    character*(*), intent(in) :: file
    integer, intent(in), optional :: lu0
    type(thread_info), intent(in), optional :: ti

    integer :: lu, i, j, k
    real*8 :: r(3,3)
    character(len=:), allocatable :: strtyp, aux

    ! open file
    if (present(lu0)) then
       lu = lu0
    else
       lu = fopen_write(file,ti=ti)
    end if

    ! atom types
    strtyp = ""
    do i = 1, c%nspc
       aux = strtyp // " " // string(nameguess(c%spc(i)%z,.true.))
       strtyp = aux
    end do

    if (c%ismolecule) then
       ! molecule
       write (lu,'(A," C")') string(c%ncel)
       write (lu,'(A)') strtyp

       ! Cartesian coordinates
       do i = 1, c%nspc
          do k = 1, c%ncel
             if (c%atcel(k)%is == i) &
                write (lu,'(99(A," "))') string(k), string(i), &
                (string(c%atcel(k)%r(j)*bohrtoa,'f',20,12),j=1,3)
          end do
       end do
    else
       ! crystal
       write (lu,'(A," F")') string(c%ncel)
       write (lu,'(A)') strtyp

       ! fractional coordinates
       do i = 1, c%nspc
          do k = 1, c%ncel
             if (c%atcel(k)%is == i) &
                write (lu,'(99(A," "))') string(k), string(i), &
                (string(c%atcel(k)%x(j),'f',20,12),j=1,3)
          end do
       end do

       ! lattice vectors
       r = c%m_x2c * bohrtoa
       write (lu,'(3(A," "))') (string(0d0,'f',20,12),j=1,3)
       do i = 1, 3
          write (lu,'(3(A," "))') (string(r(j,i),'f',20,12),j=1,3)
       end do
    endif

    ! close file
    if (.not.present(lu0)) call fclose(lu)

  end subroutine write_dftbp_gen

  !> Write a grid to a cube file. The input is the crystal (c),
  !> the grid in 3D array form (g), the filename (file), and whether
  !> to write the whole cube or only the header (onlyheader). If xd0
  !> is given, use it as the metric of the cube; otherwise, use the
  !> unit cell. If x00 is given, use it as the origin of the cube
  !> (in bohr). Otherwise, use the crystal's molx0.
  module subroutine writegrid_cube(c,g,file,onlyheader,binary,xd0,x00,ishift0,ti)
    use global, only: precisecube
    use tools_io, only: fopen_write, fclose
    use param, only: eye
    class(crystal), intent(in) :: c
    real*8, intent(in), allocatable :: g(:,:,:)
    character*(*), intent(in) :: file
    logical, intent(in) :: onlyheader
    logical, intent(in) :: binary
    real*8, intent(in), optional :: xd0(3,3)
    real*8, intent(in), optional :: x00(3)
    integer, intent(in), optional :: ishift0(3)
    type(thread_info), intent(in), optional :: ti

    integer :: n(3), i, ix, iy, lu, ishift(3)
    integer :: iix, iiy, iiz
    real*8 :: xd(3,3), x0(3), rshift(3)

    ! grid size
    if (onlyheader.and..not.allocated(g)) then
       n = 1
    else
       do i = 1, 3
          n(i) = size(g,i)
       end do
    end if

    ! process optional arguments
    ishift = 0
    if (present(ishift0)) then
       ishift = ishift0
    end if
    rshift = real(ishift,8) / n
    rshift = c%x2c(rshift)
    if (present(xd0)) then
       xd = xd0
    else
       xd = eye
       do i = 1, 3
          xd(:,i) = c%x2c(xd(:,i))
          xd(:,i) = xd(:,i) / real(n(i),8)
       end do
    endif
    if (present(x00)) then
       x0 = x00
    else
       x0 = c%molx0
    endif

    ! write
    if (binary) then
       lu = fopen_write(file,form="unformatted",ti=ti)
       write(lu) c%ncel, x0
       write(lu) n, xd
       do i = 1, c%ncel
          write(lu) c%spc(c%atcel(i)%is)%z, 0d0, c%atcel(i)%r(:) + c%molx0 - rshift
       end do
       write (lu) g
    else
       lu = fopen_write(file,ti=ti)
       write(lu,'("critic2-cube")')
       write(lu,'("critic2-cube")')
       if (precisecube) then
          write(lu,'(I5," ",3(E22.14," "))') c%ncel, x0
          do i = 1, 3
             write(lu,'(I5," ",3(E22.14," "))') n(i), xd(:,i)
          end do
          do i = 1, c%ncel
             write(lu,'(I4,F5.1,3(E22.14," "))') c%spc(c%atcel(i)%is)%z, 0d0, &
                c%atcel(i)%r(:) + c%molx0 - rshift
          end do
       else
          write(lu,'(I5,3(F12.6))') c%ncel, x0
          do i = 1, 3
             write(lu,'(I5,3(F12.6))') n(i), xd(:,i)
          end do
          do i = 1, c%ncel
             write(lu,'(I4,F5.1,F11.6,F11.6,F11.6)') c%spc(c%atcel(i)%is)%z, 0d0, &
                c%atcel(i)%r(:) + c%molx0 - rshift
          end do
       end if
       if (.not.onlyheader) then
          do iix = 0, n(1)-1
             ix = modulo(iix + ishift(1),n(1)) + 1
             do iiy = 0, n(2)-1
                iy = modulo(iiy + ishift(2),n(2)) + 1
                if (precisecube) then
                   write (lu,'(6(" ",e22.14))') (g(ix,iy,modulo(iiz+ishift(3),n(3))+1),iiz=0,n(3)-1)
                else
                   write (lu,'(1p,6(" ",e12.5))') (g(ix,iy,modulo(iiz+ishift(3),n(3))+1),iiz=0,n(3)-1)
                end if
             enddo
          enddo
       end if
    end if
    call fclose(lu)

  end subroutine writegrid_cube

  !> Write a grid to a VASP CHGCAR file. The input is the crystal (c),
  !> the grid in 3D array form (g), the filename (file), and whether
  !> to write the whole cube or only the header (onlyheader).
  module subroutine writegrid_vasp(c,g,file,onlyheader,ishift0,ti)
    use tools_io, only: fopen_write, string, nameguess, fclose
    use param, only: bohrtoa
    class(crystal), intent(in) :: c
    real*8, intent(in), allocatable :: g(:,:,:)
    character*(*), intent(in) :: file
    logical, intent(in) :: onlyheader
    integer, intent(in), optional :: ishift0(3)
    type(thread_info), intent(in), optional :: ti

    integer :: n(3), i, j, ix, iy, iz, lu, ishift(3), nat
    character(len=:), allocatable :: line0, aux
    real*8 :: xshift(3)

    ! grid size
    if (onlyheader.and..not.allocated(g)) then
       n = 1
    else
       do i = 1, 3
          n(i) = size(g,i)
       end do
    end if

    ! process optional arguments
    ishift = 0
    if (present(ishift0)) then
       ishift = ishift0
    end if
    xshift = real(ishift,8) / n

    lu = fopen_write(file,ti=ti)
    write (lu,'("CHGCAR generated by critic2")')
    write (lu,'("1.0")')
    do i = 1, 3
       write (lu,'(1p,3(E22.14," "))') c%m_x2c(:,i) * bohrtoa
    end do

    ! species line
    line0 = ""
    do i = 1, c%nspc
       aux = line0 // " " // string(nameguess(c%spc(i)%z,.true.))
       line0 = aux
    end do
    write (lu,'(A)') line0

    ! number of atoms line
    line0 = ""
    do i = 1, c%nspc
       nat = 0
       do j = 1, c%ncel
          if (c%atcel(j)%is == i) nat = nat + 1
       end do
       aux = line0 // " " // string(nat)
       line0 = aux
    end do
    write (lu,'(A)') line0

    ! positions
    write (lu,'("Direct")')
    do i = 1, c%nspc
       do j = 1, c%ncel
          if (c%atcel(j)%is == i) &
             write (lu,'(1p,3(E22.14," "))') c%atcel(j)%x - xshift
       end do
    end do
    write (lu,*)
    write (lu,'(3(I5," "))') n
    if (.not.onlyheader) then
       write (lu,'(5(" ",e22.14))') &
          (((g(modulo(ix+ishift(1),n(1))+1,modulo(iy+ishift(2),n(2))+1,modulo(iz+ishift(3),n(3))+1)*c%omega,&
          ix=0,n(1)-1),iy=0,n(2)-1),iz=0,n(3)-1)
    end if
    call fclose(lu)

  end subroutine writegrid_vasp

  !> Write a grid to a xsf file. The input is the crystal (c), the
  !> grid in 3D array form (g), the filename (file), and whether to
  !> write the whole xsf or only the structure (onlyheader).
  module subroutine writegrid_xsf(c,g,file,onlyheader,ishift0,ti)
    use tools_io, only: fopen_write, fclose, string, ferror, faterr
    use param, only: bohrtoa
    class(crystal), intent(in) :: c
    real*8, intent(in), allocatable :: g(:,:,:)
    character*(*), intent(in) :: file
    logical, intent(in) :: onlyheader
    integer, intent(in), optional :: ishift0(3)
    type(thread_info), intent(in), optional :: ti

    real*8 :: x(3), rshift(3)
    integer :: n(3), i, j, ix, iy, iz, lu, ishift(3)

    if (c%ismolecule) then
       call ferror('writegrid_xsf','molecular xsf files not supported yet',faterr)
    end if

    ! grid size
    if (onlyheader .and..not.allocated(g)) then
       n = 1
    else
       do i = 1, 3
          n(i) = size(g,i)
       end do
    end if

    ! process optional arguments
    ishift = 0
    if (present(ishift0)) then
       ishift = ishift0
    end if
    rshift = real(ishift,8) / n
    rshift = c%x2c(rshift)

    lu = fopen_write(file,ti=ti)
    write (lu,'("## xsf file generated by critic2")')
    write (lu,'("CRYSTAL")')
    write (lu,'("PRIMVEC")')
    do i = 1, 3
       write (lu,'(1p,3(E22.14," "))') c%m_x2c(:,i) * bohrtoa
    end do
    write (lu,'("CONVVEC")')
    do i = 1, 3
       write (lu,'(1p,3(E22.14," "))') c%m_x2c(:,i) * bohrtoa
    end do

    write (lu,'("PRIMCOORD")')
    write (lu,'(A," 1")') string(c%ncel)
    do j = 1, c%ncel
       write (lu,'(A," ",1p,3(E22.14," "))') string(c%spc(c%atcel(j)%is)%z), &
          (c%atcel(j)%r - rshift) * bohrtoa
    end do
    if (.not.onlyheader) then
       write (lu,'("BEGIN_BLOCK_DATAGRID3D")')
       write (lu,'("critic2:grid")')
       write (lu,'("DATAGRID_3D_GRID")')
       write (lu,'(3(A," "))') (string(n(j)+1),j=1,3)
       x = 0d0
       write (lu,'(1p,3(" ",e22.14))') x
       do i = 1, 3
          x = 0d0
          x(i) = 1d0
          x = c%x2c(x) * bohrtoa
          write (lu,'(1p,3(" ",e22.14))') x
       end do
       write (lu,'(1p,5(" ",e22.14))') (((g(modulo(ix+ishift(1),n(1))+1,modulo(iy+ishift(2),n(2))+1,&
          modulo(iz+ishift(3),n(3))+1),ix=0,n(1)),iy=0,n(2)),iz=0,n(3))
       write (lu,'("END_DATAGRID_3D")')
       write (lu,'("END_BLOCK_DATAGRID3D")')
    end if

    call fclose(lu)

  end subroutine writegrid_xsf

  !> Translate the point x0 to the main cell and calculate the core
  !> (if zpsp is present) or promolecular densities at a point x0
  !> (coord format given by icrd) using atomic radial grids up to a
  !> number of derivatives nder (max: 2). Returns the density (f),
  !> gradient (fp, nder >= 1), and Hessian (fpp, nder >= 2). If a
  !> fragment (fr) is given, then only the atoms in it
  !> contribute. This routine is a wrapper for the environment's
  !> promolecular. Thread-safe.
  module subroutine promolecular(c,x0,icrd,f,fp,fpp,nder,zpsp,fr)
    use fragmentmod, only: fragment
    class(crystal), intent(in) :: c
    real*8, intent(in) :: x0(3) !< Point in cryst. coords.
    integer, intent(in) :: icrd !< Input coordinate format
    real*8, intent(out) :: f !< Density
    real*8, intent(out) :: fp(3) !< Density gradient
    real*8, intent(out) :: fpp(3,3) !< Density hessian
    integer, intent(in) :: nder !< Number of derivatives to calculate
    integer, intent(in), optional :: zpsp(:)
    type(fragment), intent(in), optional :: fr !< Fragment contributing to the density

    call c%env%promolecular(x0,icrd,f,fp,fpp,nder,zpsp,fr)

  end subroutine promolecular

  !> Calculate the core or promolecular densities on a grid with n(:)
  !> points. If a fragment is given, then only the atoms in it
  !> contribute.  This routine is thread-safe.
  module subroutine promolecular_grid(c,f,n,zpsp,fr)
    use grid3mod, only: grid3
    use grid1mod, only: grid1
    use fragmentmod, only: fragment
    use param, only: icrd_crys
    class(crystal), intent(in) :: c
    type(grid3), intent(out) :: f
    integer, intent(in) :: n(3)
    integer, intent(in), optional :: zpsp(:)
    type(fragment), intent(in), optional :: fr

    integer :: i, j, k
    real*8 :: x(3), xdelta(3,3), rdum1(3), rdum2(3,3), rho

    ! initialize
    call f%end()
    f%isinit = .true.
    f%n = n
    allocate(f%f(n(1),n(2),n(3)))

    do i = 1, 3
       xdelta(:,i) = 0d0
       xdelta(i,i) = 1d0 / real(n(i),8)
    end do

    !$omp parallel do private(x,rho,rdum1,rdum2)
    do k = 1, n(3)
       do j = 1, n(2)
          do i = 1, n(1)
             x = (i-1) * xdelta(:,1) + (j-1) * xdelta(:,2) + (k-1) * xdelta(:,3)
             call c%promolecular(x,icrd_crys,rho,rdum1,rdum2,0,zpsp,fr)

             !$omp critical(write)
             f%f(i,j,k) = rho
             !$omp end critical(write)
          end do
       end do
    end do
    !$omp end parallel do

  end subroutine promolecular_grid

  !> Get the holohedry and the Laue class from the Hermann-Mauguin
  !> point group label. Adapted from spglib, takes spglib HM point
  !> group labels.
  module subroutine pointgroup_info(hmpg,schpg,holo,laue)
    use tools_io, only: equal
    character*(*), intent(in) :: hmpg
    character(len=3), intent(out) :: schpg
    integer, intent(out) :: holo
    integer, intent(out) :: laue

    if (equal(hmpg,"")) then
       schpg = ""
       holo = holo_unk
       laue = laue_unk
    elseif (equal(hmpg,"1")) then
       schpg = "C1"
       holo = holo_tric
       laue = laue_1
    elseif (equal(hmpg,"-1")) then
       schpg = "Ci"
       holo = holo_tric
       laue = laue_1
    elseif (equal(hmpg,"2")) then
       schpg = "C2"
       holo = holo_mono
       laue = laue_2m
    elseif (equal(hmpg,"m")) then
       schpg = "Cs"
       holo = holo_mono
       laue = laue_2m
    elseif (equal(hmpg,"2/m")) then
       schpg = "C2h"
       holo = holo_mono
       laue = laue_2m
    elseif (equal(hmpg,"222")) then
       schpg = "D2"
       holo = holo_ortho
       laue = laue_mmm
    elseif (equal(hmpg,"mm2")) then
       schpg = "C2v"
       holo = holo_ortho
       laue = laue_mmm
    elseif (equal(hmpg,"mmm")) then
       schpg = "D2h"
       holo = holo_ortho
       laue = laue_mmm
    elseif (equal(hmpg,"4")) then
       schpg = "C4"
       holo = holo_tetra
       laue = laue_4m
    elseif (equal(hmpg,"-4")) then
       schpg = "S4"
       holo = holo_tetra
       laue = laue_4m
    elseif (equal(hmpg,"4/m")) then
       schpg = "C4h"
       holo = holo_tetra
       laue = laue_4m
    elseif (equal(hmpg,"422")) then
       schpg = "D4"
       holo = holo_tetra
       laue = laue_4mmm
    elseif (equal(hmpg,"4mm")) then
       schpg = "C4v"
       holo = holo_tetra
       laue = laue_4mmm
    elseif (equal(hmpg,"-42m")) then
       schpg = "D2d"
       holo = holo_tetra
       laue = laue_4mmm
    elseif (equal(hmpg,"4/mmm")) then
       schpg = "D4h"
       holo = holo_tetra
       laue = laue_4mmm
    elseif (equal(hmpg,"3")) then
       schpg = "C3"
       holo = holo_trig
       laue = laue_3
    elseif (equal(hmpg,"-3")) then
       schpg = "C3i"
       holo = holo_trig
       laue = laue_3
    elseif (equal(hmpg,"32")) then
       schpg = "D3"
       holo = holo_trig
       laue = laue_3m
    elseif (equal(hmpg,"3m")) then
       schpg = "C3v"
       holo = holo_trig
       laue = laue_3m
    elseif (equal(hmpg,"-3m")) then
       schpg = "D3d"
       holo = holo_trig
       laue = laue_3m
    elseif (equal(hmpg,"6")) then
       schpg = "C6"
       holo = holo_hex
       laue = laue_6m
    elseif (equal(hmpg,"-6")) then
       schpg = "C3h"
       holo = holo_hex
       laue = laue_6m
    elseif (equal(hmpg,"6/m")) then
       schpg = "C6h"
       holo = holo_hex
       laue = laue_6m
    elseif (equal(hmpg,"622")) then
       schpg = "D6"
       holo = holo_hex
       laue = laue_6mmm
    elseif (equal(hmpg,"6mm")) then
       schpg = "C6v"
       holo = holo_hex
       laue = laue_6mmm
    elseif (equal(hmpg,"-6m2")) then
       schpg = "D3h"
       holo = holo_hex
       laue = laue_6mmm
    elseif (equal(hmpg,"6/mmm")) then
       schpg = "D6h"
       holo = holo_hex
       laue = laue_6mmm
    elseif (equal(hmpg,"23")) then
       schpg = "T"
       holo = holo_cub
       laue = laue_m3
    elseif (equal(hmpg,"m-3")) then
       schpg = "Th"
       holo = holo_cub
       laue = laue_m3
    elseif (equal(hmpg,"432")) then
       schpg = "O"
       holo = holo_cub
       laue = laue_m3m
    elseif (equal(hmpg,"-43m")) then
       schpg = "Td"
       holo = holo_cub
       laue = laue_m3m
    elseif (equal(hmpg,"m-3m")) then
       schpg = "Oh"
       holo = holo_cub
       laue = laue_m3m
    end if

  end subroutine pointgroup_info

  !xx! private procedures

  !> For the symmetry operation rot, compute the order and
  !> characteristic direction of the operation. The output is given in
  !> type (the type of operation: ident, c2, c3,...)  and vec (the
  !> direction of the axis or the normal to the plane, 0 for a point
  !> symmetry element). Used in the structure initialization.
  subroutine typeop(rot,type,vec,order)
    use tools_math, only: eig
    use tools_io, only: ferror, faterr
    use param, only: tpi, eye

    real*8, intent(in) :: rot(3,4) !< rotm operation
    integer, intent(out) :: type !< output type
    real*8, dimension(3), intent(out) :: vec !< output eigenvector
    integer, intent(out) :: order

    real*8, parameter :: eps = 1e-3

    integer, dimension(2,2) :: iord
    real*8 :: trace, tone
    real*8, dimension(3) :: eigen, eigeni
    real*8, dimension(3,3) :: mat, mat3
    integer :: nm, nones, nminusones
    integer :: i

    trace = rot(1,1) + rot(2,2) + rot(3,3)
    eigen = 0d0
    eigeni = 0d0
    mat = rot(1:3,1:3)
    iord = 0
    vec = 0d0

    !.Determine type of operation
    order = 1
    if (abs(trace+3) .lt. eps) then
       ! inversion
       order = 2
       type = inv
       return
    elseif (abs(trace-3) .lt.eps) then
       ! identity
       type = ident
       return
    else
       nm=3
       ! It is an axis or a plane. Diagonalize.
       call eig(mat,3,eigen,eigeni)

       ! Determine actual operation.
       nones = 0
       nminusones = 0
       do i = 1, 3
          if (abs(eigen(i)-1) < eps .and. abs(eigeni(i)) < eps) then
             nones = nones + 1
             iord(1,nones) = i
          elseif (abs(eigen(i)+1) < eps .and. abs(eigeni(i)) < eps) then
             nminusones = nminusones +1
             iord(2,nminusones) = i
          endif
       enddo

       ! What is it?
       if (nones .eq. 1) then
          !.A proper axis. Order?
          tone = min(max((trace-1) / 2,-1d0),1d0)
          order = nint(abs(tpi/(acos(tone))))
          do i = 1,3
             vec(i) = mat(i,iord(1,1))
          enddo
          if (order .eq. 2) then
             type = c2
          elseif (order .eq. 3) then
             type = c3
          elseif (order .eq. 4) then
             type = c4
          elseif (order .eq. 6) then
             type = c6
          else
             call ferror ('typeop','Axis unknown',faterr)
          endif
       elseif (nones .eq. 2) then
          !A plane, Find directions.
          order = 2
          type = sigma
          do i = 1,3
             vec(i) = mat(i,iord(2,1))
          enddo
       elseif (nminusones .eq. 1) then
          !.An improper axes. Order?
          tone = min(max((trace+1) / 2,-1d0),1d0)
          order = nint(abs(tpi/(acos(tone))))
          do i = 1,3
             vec(i) = mat(i,iord(2,1))
          enddo
          if (order .eq. 3) then
             type = s3
          elseif (order .eq. 4) then
             type = s4
          elseif (order .eq. 6) then
             mat3 = rot(1:3,1:3)
             mat3 = matmul(matmul(mat3,mat3),mat3)
             if (all(abs(mat3+eye) < 1d-5)) then
                type = s3
             else
                type = s6
             end if
          else
             call ferror ('typeop','Axis unknown',faterr)
          endif
       else
          call ferror ('typeop', 'Sym. Element unknown', faterr)
       endif
    endif
    vec = vec / norm2(vec)

  end subroutine typeop

  !> Private for wigner. Determines if two tetrahedra are equivalent.
  function equiv_tetrah(c,x0,t1,t2,leqv,lrotm,eps)
    logical :: equiv_tetrah
    type(crystal), intent(in) :: c
    real*8, intent(in) :: x0(3)
    real*8, dimension(3,0:3), intent(in) :: t1, t2
    integer, intent(in) :: leqv
    real*8, intent(in) :: lrotm(3,3,48), eps

    integer :: i, k, p
    real*8 :: r1(3,0:3), d1(3,0:3), dist2(0:3), xdum(3)

    equiv_tetrah = .false.

    do i = 1, leqv
       do k = 0, 3
          r1(:,k) = matmul(lrotm(:,1:3,i),t1(:,k)-x0) + x0
       end do

       do p = 1, 6
          d1 = perm3(p,r1,t2)
          do k = 0, 3
             xdum = d1(:,k)
             d1(:,k) = c%x2c(xdum)
             dist2(k) = norm2(d1(:,k))
          end do
          if (all(dist2 < eps)) then
             equiv_tetrah = .true.
             return
          end if
       end do
    end do

  end function equiv_tetrah

  !> Private for equiv_tetrah, wigner. 3! permutations.
  function perm3(p,r,t) result(res)
    integer, intent(in) :: p
    real*8, intent(in) :: r(3,0:3), t(3,0:3)
    real*8 :: res(3,0:3)

    select case(p)
    case(1)
       res(:,0) = r(:,0) - t(:,0)
       res(:,1) = r(:,1) - t(:,1)
       res(:,2) = r(:,2) - t(:,2)
       res(:,3) = r(:,3) - t(:,3)
    case(2)
       res(:,0) = r(:,0) - t(:,0)
       res(:,1) = r(:,1) - t(:,1)
       res(:,2) = r(:,2) - t(:,3)
       res(:,3) = r(:,3) - t(:,2)
    case(3)
       res(:,0) = r(:,0) - t(:,0)
       res(:,1) = r(:,1) - t(:,2)
       res(:,2) = r(:,2) - t(:,1)
       res(:,3) = r(:,3) - t(:,3)
    case(4)
       res(:,0) = r(:,0) - t(:,0)
       res(:,1) = r(:,1) - t(:,2)
       res(:,2) = r(:,2) - t(:,3)
       res(:,3) = r(:,3) - t(:,1)
    case(5)
       res(:,0) = r(:,0) - t(:,0)
       res(:,1) = r(:,1) - t(:,3)
       res(:,2) = r(:,2) - t(:,1)
       res(:,3) = r(:,3) - t(:,2)
    case(6)
       res(:,0) = r(:,0) - t(:,0)
       res(:,1) = r(:,1) - t(:,3)
       res(:,2) = r(:,2) - t(:,2)
       res(:,3) = r(:,3) - t(:,1)
    end select

  end function perm3

end submodule proc
