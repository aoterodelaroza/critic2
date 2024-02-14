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

! Structure class and routines for molecular and crystal structure computations
submodule (crystalmod) proc
  implicit none

contains

  !> Initialize the struct arrays
  module subroutine struct_init(c)
    class(crystal), intent(inout) :: c

    integer :: i

    ! array initialization values
    integer, parameter :: mspc0 = 4
    integer, parameter :: mneq0 = 4
    integer, parameter :: mcel0 = 10

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
    c%vaclength = 0d0
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
    use global, only: crsmall, atomeps_structnew
    use tools_math, only: m_x2c_from_cellpar, m_c2x_from_cellpar, matinv, &
       det3, mnorm2
    use tools_io, only: ferror, faterr, zatguess, string
    use tools, only: wscell, qcksort
    use types, only: realloc
    use param, only: pi, eyet, eye, icrd_cart
    class(crystal), intent(inout) :: c
    type(crystalseed), intent(in) :: seed
    logical, intent(in) :: crashfail
    logical, intent(in), optional :: noenv
    type(thread_info), intent(in), optional :: ti

    real*8 :: g(3,3), xmax(3), xmin(3), xcm(3), dist, border, xx(3)
    logical :: good, good2, clearsym, doenv
    integer :: i, j, k, l, iat, newmult
    real*8, allocatable :: atpos(:,:), area(:), xcoord(:)
    integer, allocatable :: irotm(:), icenv(:), iord(:)
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
                   if (c%eql_distance(seed%x(:,l),xx) < atomeps_structnew) useatom(l) = .false.
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
                if (all(abs(eye - c%rotm(1:3,1:3,i)) < 1d-12)) then
                   good2 = .false.
                   do j = 1, c%ncv
                      if (all(abs(c%rotm(:,4,i) - c%cen(:,j)) < 1d-12) .or.&
                          all(abs(c%rotm(:,4,i) + c%cen(:,j)) < 1d-12)) then
                         good2 = .true.
                      end if
                   end do
                   if (.not.good2) then
                      if (crashfail) then
                         call ferror('struct_new','identity operation in rotm with unknown translation vector',faterr)
                      else
                         return
                      end if
                   end if
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
             call c%symeqv(c%at(i)%x,c%at(i)%mult,atpos,irotm,icenv,atomeps_structnew)

             newmult = 0
             jloop: do j = 1, c%at(i)%mult
                iat = iat + 1
                newmult = newmult + 1
                c%atcel(iat)%x = atpos(:,j)
                c%atcel(iat)%r = c%x2c(atpos(:,j))

                c%atcel(iat)%rxc = c%x2xr(atpos(:,j))
                c%atcel(iat)%rxc = c%atcel(iat)%rxc - floor(c%atcel(iat)%rxc)
                c%atcel(iat)%rxc = c%xr2c(c%atcel(iat)%rxc)

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

       ! set up the block environments
       call c%build_env()

       ! find atomic connectivity and molecular fragments
       call c%find_asterisms_covalent()
       call c%fill_molecular_fragments()
       call c%calculate_molecular_equivalence()

       ! calculate vacuum lengths
       allocate(xcoord(2*c%ncel),iord(2*c%ncel))
       do i = 1, 3
          do j = 1, c%ncel
             xcoord(j) = (c%atcel(j)%x(i) - floor(c%atcel(j)%x(i))) * c%aa(i)
             xcoord(c%ncel+j) = xcoord(j) + c%aa(i)
             iord(j) = j
             iord(c%ncel+j) = c%ncel + j
          end do
          call qcksort(xcoord,iord,1,2*c%ncel)
          c%vaclength(i) = 0d0
          do j = 2, 2*c%ncel
             c%vaclength(i) = max(c%vaclength(i),xcoord(iord(j))-xcoord(iord(j-1)))
          end do
       end do
       deallocate(xcoord,iord)

       ! Write the half nearest-neighbor distance
       do i = 1, c%nneq
          if (.not.c%ismolecule .or. c%ncel > 1) then
             call c%nearest_atom_env(c%at(i)%r,icrd_cart,iat,dist,nozero=.true.)
             c%at(i)%rnn2 = 0.5d0 * dist
          else
             c%at(i)%rnn2 = 0d0
          end if
       end do
    end if

    ! the initialization is done - this crystal is ready to use
    c%file = seed%file
    c%isinit = .true.

  end subroutine struct_new

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
             call c%nearest_atom_env(x,icrd_crys,nid,dist)
             !$omp critical(write)
             idg(i,j,k) = nid
             !$omp end critical(write)
          end do
       end do
    end do
    !$omp end parallel do

  end subroutine nearest_atom_grid

  !> Identify a species in the structure from a string. Step 1: check
  !> if str is equal to the name of any species (case sensitive). 2:
  !> same but case insensitive.  3: transform to Z and check if any
  !> species has that Z (if several, return the first one). 4: return
  !> 0 (not found).
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

end submodule proc
