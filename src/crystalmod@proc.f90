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

submodule (crystalmod) proc
  implicit none

  !xx! private procedures
  ! subroutine lattpg(rmat,ncen,xcen,nn,rot)
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
    use param, only: eyet
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
    c%neqv = 1
    c%rotm = 0d0
    c%rotm(:,:,1) = eyet
    c%ncv = 1
    if (allocated(c%cen)) deallocate(c%cen)
    allocate(c%cen(3,4))
    c%cen = 0d0
    c%spg%n_atoms = 0

    ! no molecule
    c%ismolecule = .false.
    c%molx0 = 0d0
    c%molborder = 0d0

    ! no ws
    c%ws_nv = 0
    c%ws_nf = 0
    c%ws_mnfv = 0
    c%ws_ineighx = 0
    c%ws_ineighc = 0d0
    c%ws_nside = 00
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

  end subroutine struct_end

  !> Create a new, complete crystal/molecule from a crystal seed. If
  !> failed and crashfail is true, crash the program. Otherwise,
  !> return a the error status through c%isinit.
  module subroutine struct_new(c,seed,crashfail)
    use crystalseedmod, only: crystalseed
    use grid1mod, only: grid1_register_ae
    use global, only: crsmall, atomeps
    use tools_math, only: m_x2c_from_cellpar, m_c2x_from_cellpar, matinv, &
       det, mnorm2
    use tools_io, only: ferror, faterr, zatguess, string
    use types, only: realloc
    use param, only: pi, eyet, icrd_cart
    class(crystal), intent(inout) :: c
    type(crystalseed), intent(in) :: seed
    logical, intent(in) :: crashfail
    
    real*8 :: g(3,3), xmax(3), xmin(3), xcm(3), dist
    logical :: good, hasspg
    integer :: i, j, iat
    real*8, allocatable :: atpos(:,:), rnn2(:)
    integer, allocatable :: irotm(:), icenv(:)

    if (.not.seed%isused) then
       if (crashfail) then
          call ferror("struct_new","uninitialized seed",faterr)
       else
          return
       end if
    end if

    ! initialize the structure
    call c%init()

    ! copy the atomic information
    c%nspc = seed%nspc
    c%nneq = seed%nat
    if (c%nneq > 0) then
       if (allocated(c%at)) deallocate(c%at)
       allocate(c%at(c%nneq))
       if (allocated(c%spc)) deallocate(c%spc)
       allocate(c%spc(c%nspc))
       c%spc = seed%spc(1:seed%nspc)
       do i = 1, c%nneq
          c%at(i)%x = seed%x(:,i)
          c%at(i)%is = seed%is(i)
       end do

       ! if this is a molecule, calculate the center and encompassing cell
       if (seed%ismolecule) then
          xmax = -1d40
          xmin =  1d40
          xcm = 0d0
          do i = 1, seed%nat
             do j = 1, 3
                xmax(j) = max(seed%x(j,i)+seed%border,xmax(j))
                xmin(j) = min(seed%x(j,i)-seed%border,xmin(j))
             end do
             xcm = xcm + seed%x(:,i)
          end do
          xcm = xcm / seed%nat
          if (seed%cubic) then
             xmin = minval(xmin)
             xmax = maxval(xmax)
          end if
       end if
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
       c%m_c2x = matinv(c%m_x2c)
       g = matmul(transpose(c%m_x2c),c%m_x2c)
    elseif (seed%useabr == 1) then
       ! use aa and bb
       c%aa = seed%aa
       c%bb = seed%bb
       c%m_x2c = m_x2c_from_cellpar(c%aa,c%bb)
       c%m_c2x = matinv(c%m_x2c)
       g = matmul(transpose(c%m_x2c),c%m_x2c)
    elseif (seed%useabr == 2) then
       ! use m_x2c
       c%m_x2c = seed%m_x2c
       c%m_c2x = matinv(c%m_x2c)
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

    ! transform the atomic coordinates in the case of a molecule, and fill 
    ! the remaining molecular fields
    if (seed%ismolecule) then
       c%ismolecule = seed%ismolecule

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
          c%molborder = max(seed%border - max(2d0,0.8d0 * seed%border),0d0) / (xmax - xmin)
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
    ! Cartesian coordinates
    do i = 1, c%nneq
       c%at(i)%x = c%at(i)%x - floor(c%at(i)%x)
       c%at(i)%r = c%x2c(c%at(i)%x)
    end do

    ! rest of the cell metrics
    c%gtensor = g
    c%omega = sqrt(max(det(g),0d0))
    c%grtensor = matinv(g)
    do i = 1, 3
       c%ar(i) = sqrt(c%grtensor(i,i))
    end do
    c%n2_x2c = mnorm2(c%m_x2c)
    c%n2_c2x = mnorm2(c%m_c2x)

    ! calculate the wigner-seitz cell
    call c%wigner()

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

             do j = 1, c%at(i)%mult
                iat = iat + 1
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
             end do
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
       c%havesym = 0
       c%neqv = 1
       c%rotm = 0d0
       c%rotm(:,:,1) = eyet
       c%ncv = 1
       if (.not.allocated(c%cen)) allocate(c%cen(3,4))
       c%cen = 0d0
    end if

    ! symmetry from spglib
    hasspg = .false.
    if (.not.seed%ismolecule .and. seed%havesym == 0 .and. (seed%findsym == 1 .or. seed%findsym == -1 .and. seed%nat <= crsmall)) then
       ! symmetry was not available, and I want it
       ! this operation fills the symmetry info, at(i)%mult, and ncel/atcel
       call c%spglib_wrap(.true.,.false.)
       hasspg = .true.
    else if (c%havesym > 0 .and..not.hasspg) then
       ! symmetry was already available, but I still want the space group details
       call c%spglib_wrap(.false.,.true.)
    else
       ! symmetry was not available, and I do not want it - make a copy of at() to atcel()
       c%ncel = c%nneq
       if (allocated(c%atcel)) deallocate(c%atcel)
       allocate(c%atcel(c%ncel))
       do i = 1, c%ncel
          c%at(i)%mult = 1
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

    ! load the atomic density grids
    do i = 1, c%nspc
       call grid1_register_ae(c%spc(i)%z)
    end do

    if (c%ncel > 0) then
       ! Build the atomic environments
       call c%env%build(c%ismolecule,c%nspc,c%spc(1:c%nspc),c%ncel,c%atcel(1:c%ncel),c%m_xr2c,c%m_x2xr,c%m_x2c)

       ! Find the atomic connectivity and the molecular fragments
       call c%env%find_asterisms_covalent(c%nstar,rnn2)
       call c%fill_molecular_fragments()
       
       ! Write the half nearest-neighbor distance
       do i = 1, c%nneq
          if (rnn2(i) > 0d0) then
             c%at(i)%rnn2 = 0.5d0 * rnn2(i)
          elseif (.not.c%ismolecule .or. c%ncel > 1) then
             call c%env%nearest_atom(c%at(i)%r,icrd_cart,iat,dist,nozero=.true.)
             c%at(i)%rnn2 = 0.5d0 * dist 
          end if
       end do
       if (allocated(rnn2)) deallocate(rnn2)
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
  module function identify_fragment_from_xyz(c,file) result(fr)
    use tools_io, only: fopen_read, string, ferror, faterr, fclose
    use param, only: bohrtoa, icrd_cart
    use types, only: realloc

    class(crystal), intent(in) :: c
    character*(*) :: file
    type(fragment) :: fr

    integer :: lu, nat
    integer :: id, i
    real*8 :: x0(3)
    character(len=:), allocatable :: word

    lu = fopen_read(file)
    read(lu,*,err=999) nat
    read(lu,*,err=999) 
    
    fr%nat = nat
    allocate(fr%at(nat))
    fr%nspc = c%nspc
    allocate(fr%spc(fr%nspc))
    fr%spc = c%spc
    do i = 1, nat
       word = ""
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
    use tools_io, only: ferror, faterr, warning
    use types, only: realloc
    class(crystal), intent(inout) :: c

    integer :: i, j, k, l, jid, newid, newl(3)
    integer :: nat, nlvec, lwork, info
    logical :: found, fdisc
    integer, allocatable :: id(:), lvec(:,:)
    logical, allocatable :: ldone(:), used(:)
    real*8, allocatable :: rlvec(:,:), sigma(:), uvec(:,:), vvec(:,:), work(:)
    real*8 :: xcm(3), x(3)

    if (.not.allocated(c%nstar)) &
       call ferror('fill_molecular_fragments','no asterisms found',faterr)
    if (allocated(c%mol)) deallocate(c%mol)

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

  end subroutine fill_molecular_fragments

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

  !> Calculate the powder diffraction pattern. 
  !> On input, npts is the number of 2*theta points from the initial
  !> (th2ini0) to the final (th2end0) 2*theta values. Both angles are
  !> in degrees. lambda0 is the wavelength of the radiation (in
  !> angstrom). sigma is the parameter for the Gaussian broadening.
  !> fpol is the polarization correction factor (0 = unpolarized, 0.95
  !> = syncrhotron).
  !> On output, t is the 2*theta grid, ih is the intensity on the
  !> 2*theta grid, th2p is the 2*theta for the located maxima, ip is
  !> the list of maxima itensities, and hvecp is the reciprocal
  !> lattice vector corresponding to the peaks.
  module subroutine powder(c,th2ini0,th2end0,npts,lambda0,fpol,&
     sigma,t,ih,th2p,ip,hvecp)
    use param, only: pi, bohrtoa, cscatt, c2scatt
    use tools_io, only: ferror, faterr
    use tools, only: qcksort
    use types, only: realloc
    class(crystal), intent(in) :: c
    real*8, intent(in) :: th2ini0, th2end0
    integer, intent(in) :: npts
    real*8, intent(in) :: lambda0
    real*8, intent(in) :: fpol
    real*8, intent(in) :: sigma
    real*8, allocatable, intent(inout) :: t(:)
    real*8, allocatable, intent(inout) :: ih(:)
    real*8, allocatable, intent(inout) :: th2p(:)
    real*8, allocatable, intent(inout) :: ip(:)
    integer, allocatable, intent(inout) :: hvecp(:,:)

    integer :: i, ii, np, hcell, h, k, l, iz, idx
    real*8 :: th2ini, th2end, lambda, hvec(3), kvec(3), th, sth, th2
    real*8 :: sigma2, smax, dh2, dh, dh3, sthlam, cterm, sterm
    real*8 :: ffac, as(4), bs(4), cs, c2s(4), int, mcorr, afac
    real*8 :: ipmax, ihmax
    integer :: hmax
    integer, allocatable :: multp(:)
    integer, allocatable :: io(:)
    real*8, allocatable :: th2p_(:), ip_(:)
    integer, allocatable :: hvecp_(:,:)

    integer, parameter :: mp = 20
    real*8, parameter :: ieps = 1d-5
    real*8, parameter :: theps = 1d-5

    ! prepare the grid limits
    if (allocated(t)) deallocate(t)
    if (allocated(ih)) deallocate(ih)
    allocate(t(npts),ih(npts))
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
    smax = sin(th2end/2d0)
    hmax = 2*ceiling(2*smax/lambda/minval(c%ar))
    ! broadening -> gaussian
    sigma2 = sigma * sigma

    ! calculate the intensities
    np = 0
    do hcell = 1, hmax
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
                if (th2 < th2ini .or. th2 > th2end) cycle

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
                   ! use a gaussian profile, add the intensity
                   ih = ih + int * exp(-(t-th2*180/pi)**2 / 2d0 / sigma2)

                   ! identify the new peak
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
    allocate(io(np),th2p_(np),ip_(np),hvecp_(3,np))
    do i = 1, np
       io(i) = i
    end do
    call qcksort(th2p,io,1,np)
    do ii = 1, np
       i = io(ii)
       th2p_(ii) = th2p(i)
       ip_(ii) = ip(i)
       hvecp_(:,ii) = hvecp(:,i)
    end do
    th2p = th2p_
    ip = ip_
    hvecp = hvecp_
    deallocate(th2p_,ip_,hvecp_,io)

  end subroutine powder

  !> Calculate the radial distribution function.  On input, npts is
  !> the number of bins points from the initial (0) to the final
  !> (rend) distance and sigma is the Gaussian broadening of the
  !> peaks. On output, t is the distance grid, and ih is the value of
  !> the RDF. This routine is based on:
  !>   Willighagen et al., Acta Cryst. B 61 (2005) 29.
  !> except using the sqrt of the atomic numbers instead of the 
  !> charges.
  module subroutine rdf(c,rend,sigma,npts,t,ih)
    use param, only: icrd_cart
    use environmod, only: environ
    class(crystal), intent(in) :: c
    real*8, intent(in) :: rend
    real*8, intent(in) :: sigma
    integer, intent(in) :: npts
    real*8, allocatable, intent(inout) :: t(:)
    real*8, allocatable, intent(inout) :: ih(:)

    integer :: i, j, nat, lvec(3), ierr, iz, jz
    real*8 :: hfac, int, sigma2
    logical :: localenv
    type(environ) :: le
    integer, allocatable :: eid(:)
    real*8, allocatable :: dist(:)

    sigma2 = sigma * sigma

    ! prepare the grid limits
    if (allocated(t)) deallocate(t)
    if (allocated(ih)) deallocate(ih)
    allocate(t(npts),ih(npts))
    do i = 1, npts
       t(i) = real(i-1,8) / real(npts-1,8) * rend
    end do
    ih = 0d0

    ! prepare the environment
    localenv = .true.
    if (c%ismolecule .or. c%env%dmax0 > rend) then
       localenv = .false.
    end if
    if (localenv) then
       call le%build(c%ismolecule,c%nspc,c%spc(1:c%nspc),c%ncel,c%atcel(1:c%ncel),c%m_xr2c,c%m_x2xr,c%m_x2c)
    end if

    ! calculate the radial distribution function for the crystal
    ! RDF(r) = sum_i=1...c%nneq sum_j=1...c%env%n sqrt(Zi*Zj) / c%nneq / rij * delta(r-rij)
    hfac = (npts-1) / rend
    do i = 1, c%nneq
       iz = c%spc(c%at(i)%is)%z
       if (localenv) then
          call le%list_near_atoms(c%at(i)%r,icrd_cart,.false.,nat,eid,dist,lvec,ierr,up2d=rend,nozero=.true.)
       else
          call c%env%list_near_atoms(c%at(i)%r,icrd_cart,.false.,nat,eid,dist,lvec,ierr,up2d=rend,nozero=.true.)
       end if

       do j = 1, nat
          if (localenv) then
             jz = le%spc(le%at(eid(j))%is)%z
          else
             jz = c%env%spc(c%env%at(eid(j))%is)%z
          end if
          int = sqrt(real(iz * jz,8)) * c%at(i)%mult
          ih = ih + int * exp(-(t - dist(j))**2 / 2d0 / sigma2)
       end do
    end do
    if (.not.c%ismolecule) then
       do i = 2, npts
          ih(i) = ih(i) / t(i)**2
       end do
       ih = ih / c%ncel
    end if

  end subroutine rdf

  !> Calculate real and reciprocal space sum cutoffs
  module subroutine calculate_ewald_cutoffs(c)
    use tools_io, only: ferror, faterr
    use param, only: pi, rad, sqpi, tpi
    class(crystal), intent(inout) :: c

    real*8, parameter :: sgrow = 1.4d0
    real*8, parameter :: epscut = 1d-5
    real*8, parameter :: eeps = 1d-12

    integer :: i
    real*8 :: aux, q2sum
    integer :: ia, ib, ic
    real*8 :: alrmax(3), qq
    real*8 :: rcut1, rcut2, err_real
    real*8 :: hcut1, hcut2, err_rec

    if (c%isewald) return

    ! calculate sum of charges and charges**2
    c%qsum = 0d0
    q2sum = 0d0
    do i = 1, c%nneq
       qq = c%spc(c%at(i)%is)%qat
       if (abs(qq) < 1d-6) &
          call ferror('ewald_energy','Some of the charges are 0',faterr)
       c%qsum = c%qsum + real(c%at(i)%mult * qq,8)
       q2sum = q2sum + real(c%at(i)%mult * qq**2,8)
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
       err_real = pi * c%ncel**2 * q2sum / c%omega * c%eta**2 * erfc(rcut2 / c%eta)
    end do
    do while (rcut2-rcut1 >= epscut)
       c%rcut = 0.5*(rcut1+rcut2)
       err_real = pi * c%ncel**2 * q2sum / c%omega * c%eta**2 * erfc(c%rcut / c%eta)
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
       err_rec = c%ncel**2 * q2sum / sqpi / c%eta * erfc(c%eta * hcut2 / 2)
    end do
    do while(hcut2-hcut1 > epscut)
       c%hcut = 0.5*(hcut1+hcut2)
       err_rec = c%ncel**2 * q2sum / sqpi / c%eta * erfc(c%eta * c%hcut / 2)
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
    class(crystal), intent(inout) :: c
    real*8 :: ewe

    real*8 :: x(3)
    integer :: i

    if (.not.c%isewald) &
       call c%calculate_ewald_cutoffs()
    
    ewe = 0d0
    do i = 1, c%nneq
       x = c%at(i)%x
       ewe = ewe + c%at(i)%mult * c%spc(c%at(i)%is)%qat * &
          c%ewald_pot(x,.true.)
    end do
    ewe = ewe / 2d0

  end function ewald_energy

  !> Calculate the Ewald electrostatic potential at an arbitrary
  !> position x (crystallographic coords.)  If x is the nucleus j,
  !> return pot - q_j / |r-rj| at rj.
  module function ewald_pot(c,x,isnuc)
    use param, only: tpi, pi, sqpi
    class(crystal), intent(inout) :: c
    real*8, intent(in) :: x(3)
    logical, intent(in) :: isnuc
    real*8 :: ewald_pot

    real*8 :: nuc_cutoff2 = 1d-14

    real*8 :: rcut2, qnuc
    integer :: i, i1, i2, i3
    real*8 :: px(3), lvec(3), d2, d, dh
    real*8 :: sfac_c, sfacp, bbarg
    real*8 :: sum_real, sum_rec, sum0, sum_back

    if (.not.c%isewald) then
       !$omp critical (fill_ewald)
       call c%calculate_ewald_cutoffs()
       !$omp end critical (fill_ewald)
    end if

    ! is this a nuclear position? -> get charge
    qnuc = 0d0
    if (isnuc) then
       do i = 1, c%ncel
          px = c%atcel(i)%x - x
          d2 = dot_product(px,matmul(c%gtensor,px))
          if (d2 < nuc_cutoff2) then
             qnuc = c%spc(c%atcel(i)%is)%qat
             exit
          end if
       end do
    end if

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
    
    ! h = 0 term, apply only at the nucleus
    if (isnuc) then
       sum0 = - 2d0 * qnuc / sqpi / c%eta
    else
       sum0 = 0d0
    end if

    ! compensating background charge term
    sum_back = -c%qsum * c%eta**2 * pi / c%omega 

    ! sum up and exit
    ewald_pot = sum_real + sum_rec + sum0 + sum_back

  end function ewald_pot

  !> Given a crystal structure (c) and three lattice vectors in cryst.
  !> coords (x0(:,1), x0(:,2), x0(:,3)), build the same crystal
  !> structure using the unit cell given by those vectors. 
  module subroutine newcell(c,x00,t0,verbose0)
    use crystalseedmod, only: crystalseed
    use tools_math, only: det, matinv, mnorm2
    use tools_io, only: ferror, faterr, warning, string, uout
    use types, only: realloc
    class(crystal), intent(inout) :: c
    real*8, intent(in) :: x00(3,3)
    real*8, intent(in), optional :: t0(3)
    logical, intent(in), optional :: verbose0

    type(crystalseed) :: ncseed
    logical :: ok, found, verbose
    real*8 :: x0(3,3), x0inv(3,3), fvol
    real*8 :: x(3), dx(3), dd, t(3)
    integer :: i, j, k, l, m
    integer :: nr, nn
    integer :: nlat
    real*8, allocatable :: xlat(:,:)

    real*8, parameter :: eps = 1d-6

    if (c%ismolecule) &
       call ferror('newcell','NEWCELL incompatible with molecules',faterr)

    ! initialize
    x0 = x00
    dd = det(x0)

    if (abs(dd) < eps) then
       call ferror('newcell','invalid input vectors',faterr)
    elseif (dd < 0d0) then
       ! flip the cell
       x0 = -x0
       dd = det(x0)
       call ferror('newcell','det < 0; vectors flipped',warning)
    endif
    if (present(t0)) then
       t = t0
    else
       t = 0d0
    end if
    verbose = .false.
    if (present(verbose0)) verbose = verbose0

    ! check that the vectors are pure translations
    do i = 1, 3
       ok = .false.
       do j = 1, c%ncv
          ok = (c%are_lclose(x0(:,i),c%cen(:,j),1d-4))
          if (ok) exit
       end do
       if (.not.ok) &
          call ferror("newcell","Cell vector number " // string(i) // &
          " is not a pure translation",faterr)
    end do

    ! is this a smaller or a larger cell? Arrange vectors.
    if (abs(dd-1d0) < eps) then
       nr = 1
    elseif (dd > 1d0) then
       nr = nint(dd)
       if (abs(nr-dd) > eps) &
          call ferror('newcell','inconsistent determinant of lat. vectors',faterr)
    else
       nr = -nint(1d0/dd)
       if (abs(-nr-1d0/dd) > eps) &
          call ferror('newcell','inconsistent determinant of lat. vectors',faterr)
    end if

    if (verbose) then
       write (uout,'("* Transformation to a new unit cell (NEWCELL)")')
       write (uout,'("  Lattice vectors of the new cell in the old setting (cryst. coord.):")')
       write (uout,'(4X,3(A,X))') (string(x0(i,1),'f',12,7,4),i=1,3)
       write (uout,'(4X,3(A,X))') (string(x0(i,2),'f',12,7,4),i=1,3)
       write (uout,'(4X,3(A,X))') (string(x0(i,3),'f',12,7,4),i=1,3)
       write (uout,'("  Origin translation: ",3(A,X))') (string(t(i),'f',12,7,4),i=1,3)
       write (uout,*)
    end if

    ! inverse matrix
    x0inv = matinv(x0)

    ! metrics of the new cell
    ncseed%m_x2c = matmul(c%m_x2c,x0)
    ncseed%useabr = 2
    fvol = abs(det(x0))
    if (abs(nint(fvol)-fvol) > eps .and. abs(nint(1d0/fvol)-1d0/fvol) > eps) &
       call ferror("newcell","Inconsistent newcell volume",faterr)

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

    ! species list
    allocate(ncseed%spc(c%nspc))
    ncseed%nspc = c%nspc
    ncseed%spc = c%spc

    ! build the new atom list
    ncseed%nat = 0
    nn = nint(c%ncel * fvol)
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

    if (nr > 0) then
       if (ncseed%nat / c%ncel /= nr) then
          write (uout,*) "c%nneq = ", c%ncel
          write (uout,*) "ncseed%nat = ", ncseed%nat
          write (uout,*) "nr = ", nr
          call ferror('newcell','inconsistent cell # of atoms (nr > 0)',faterr)
       end if
    else
       if (c%ncel / ncseed%nat /= -nr) then
          write (uout,*) "c%nneq = ", c%ncel
          write (uout,*) "ncseed%nat = ", ncseed%nat
          write (uout,*) "nr = ", nr
          call ferror('newcell','inconsistent cell # of atoms (nr < 0)',faterr)
       end if
    endif

    ! rest of the seed information
    ncseed%isused = .true.
    ncseed%file = "<derived>"
    ncseed%havesym = 0
    ncseed%findsym = -1
    ncseed%ismolecule = .false.

    ! initialize the structure
    call c%struct_new(ncseed,.true.)
    if (verbose) call c%report(.true.,.true.)

  end subroutine newcell

  !> Transform to the standard cell. If toprim, convert to the
  !> primitive standard cell. If verbose, write
  !> information about the new crystal. If doforce = .true.,
  !> force the transformation to the primitive even if it does
  !> not lead to a smaller cell.
  module subroutine cell_standard(c,toprim,doforce,verbose)
    use iso_c_binding, only: c_double
    use spglib, only: spg_standardize_cell, spg_get_dataset
    use global, only: symprec
    use tools_math, only: det, matinv
    use tools_io, only: ferror, faterr, uout
    use param, only: eye
    class(crystal), intent(inout) :: c
    logical, intent(in) :: toprim
    logical, intent(in) :: doforce
    logical, intent(in) :: verbose
    
    integer :: ntyp, nat
    integer :: i, id, iprim
    real(c_double), allocatable :: x(:,:)
    integer, allocatable :: types(:)
    real*8 :: rmat(3,3)

    ! ignore molecules
    if (c%ismolecule) return

    ! use spglib transformation to the standard cell
    rmat = transpose(c%m_x2c)
    nat = c%ncel
    ntyp = c%nspc
    allocate(x(3,c%ncel),types(c%ncel))
    do i = 1, c%ncel
       x(:,i) = c%atcel(i)%x
       types(i) = c%atcel(i)%is
    end do

    iprim = 0
    if (toprim) iprim = 1
    id = spg_standardize_cell(rmat,x,types,nat,iprim,1,symprec)
    if (id == 0) &
       call ferror("cell_standard","could not find primitive cell",faterr)
    rmat = transpose(rmat)
    do i = 1, 3
       rmat(:,i) = c%c2x(rmat(:,i))
    end do

    ! flip the cell?
    if (det(rmat) < 0d0) rmat = -rmat

    ! if a primitive is wanted but det is not less than 1, do not make the change
    if (all(abs(rmat - eye) < symprec)) then
       if (verbose) &
          write (uout,'("+ Cell transformation leads to the same cell: skipping."/)')
       return
    end if
    if (toprim .and. .not.(det(rmat) < 1d0-symprec) .and..not.doforce) then
       if (verbose) &
          write (uout,'("+ Cell transformation does not lead to a smaller cell: skipping."/)')
       return
    end if

    ! rmat = transpose(matinv(c%spg%transformation_matrix))
    call c%newcell(rmat,verbose0=verbose)

  end subroutine cell_standard

  !> Transform to the Niggli cell. If verbose, write information
  !> about the new crystal.
  module subroutine cell_niggli(c,verbose)
    use spglib, only: spg_niggli_reduce
    use global, only: symprec
    use tools_io, only: ferror, faterr
    use tools_math, only: det
    class(crystal), intent(inout) :: c
    logical, intent(in) :: verbose
    
    real*8 :: rmat(3,3)
    integer :: id, i

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
    if (det(rmat) < 0d0) rmat = -rmat

    ! transform
    call c%newcell(rmat,verbose0=verbose)

  end subroutine cell_niggli

  !> Transform to the Delaunay cell. If verbose, write information
  !> about the new crystal.
  module subroutine cell_delaunay(c,verbose)
    use spglib, only: spg_delaunay_reduce
    use global, only: symprec
    use tools_io, only: ferror, faterr
    use tools_math, only: det
    class(crystal), intent(inout) :: c
    logical, intent(in) :: verbose

    real*8 :: rmat(3,3)
    integer :: id, i

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
    if (det(rmat) < 0d0) rmat = -rmat

    ! transform
    call c%newcell(rmat,verbose0=verbose)

  end subroutine cell_delaunay

  !> Transforms the current basis to the Delaunay reduced basis.
  !> Return the four Delaunay vectors in crystallographic coordinates
  !> (rmat) cell, see 9.1.8 in ITC. If rbas is present, it contains
  !> the three shortest of the seven Delaunay lattice vectors that
  !> form a cell (useful to transform to one of the delaunay reduced
  !> cells).
  module subroutine delaunay_reduction(c,rmat,rbas)
    use tools, only: qcksort
    use tools_math, only: det
    use tools_io, only: faterr, ferror
    class(crystal), intent(in) :: c
    real*8, intent(out) :: rmat(3,4)
    real*8, intent(out), optional :: rbas(3,3)
    
    integer :: i, j, k, iord(7)
    real*8 :: sc(4,4), xstar(3,7), xlen(7), dd
    logical :: again, ok

    real*8, parameter :: eps = 1d-10

    ! build the four Delaunay vectors
    rmat = 0d0
    do i = 1, 3
       rmat(i,i) = 1d0
       rmat(:,i) = c%x2c(rmat(:,i))
    end do
    rmat(:,4) = -(rmat(:,1)+rmat(:,2)+rmat(:,3))

    ! reduce until all the scalar products are negative or zero
    again = .true.
    sc = -1d0
    do while(again)
       do i = 1, 4
          do j = i+1, 4
             sc(i,j) = dot_product(rmat(:,i),rmat(:,j))
             sc(j,i) = sc(i,j)
          end do
       end do

       if (any(sc > eps)) then
          ai: do i = 1, 4
             aj: do j = i+1, 4
                if (sc(i,j) > eps) exit ai
             end do aj
          end do ai
          do k = 1, 4
             if (i == k .or. j == k) cycle
             rmat(:,k) = rmat(:,i) + rmat(:,k)
          end do
          rmat(:,i) = -rmat(:,i)
       else
          again = .false.
       end if
    end do

    if (present(rbas)) then
       xstar(:,1)  = rmat(:,1)
       xstar(:,2)  = rmat(:,2)
       xstar(:,3)  = rmat(:,3)
       xstar(:,4)  = rmat(:,4)
       xstar(:,5)  = rmat(:,1)+rmat(:,2)
       xstar(:,6)  = rmat(:,1)+rmat(:,3)
       xstar(:,7)  = rmat(:,2)+rmat(:,3)
       do i = 1, 7
          xlen(i) = norm2(xstar(:,i))
          iord(i) = i
       end do
       call qcksort(xlen,iord,1,7)
       rbas(:,1) = xstar(:,iord(1))
       ok = .false.
       iloop: do i = 2, 7
          rbas(:,2) = xstar(:,iord(i))
          do j = i+1, 7
             rbas(:,3) = xstar(:,iord(j))
             dd = det(rbas)
             if (abs(dd) > eps) then
                ok = .true.
                exit iloop
             end if
          end do
       end do iloop
       if (.not.ok) &
          call ferror("delaunay_reduction","could not find reduced basis",faterr)
       if (dd < 0d0) rbas = -rbas
       do i = 1, 3
          rbas(:,i) = c%c2x(rbas(:,i))
       end do
    end if

    do i = 1, 4
       rmat(:,i) = c%c2x(rmat(:,i))
    end do

  end subroutine delaunay_reduction

  !> Write information about the crystal structure to the output. lcrys = 
  !> information about the structure. lq = charges.
  module subroutine struct_report(c,lcrys,lq)
    use global, only: iunitname0, dunit0, iunit
    use tools_math, only: gcd
    use tools_io, only: uout, string, ioj_center, ioj_left, ioj_right
    use param, only: bohrtoa, maxzat, pi
    class(crystal), intent(in) :: c
    logical, intent(in) :: lcrys
    logical, intent(in) :: lq

    integer :: i, j, k, iz, is
    integer :: nelec, holo, laue
    real*8 :: maxdv, xcm(3), x0(3), xlen(3), xang(3), xred(3,3)
    character(len=:), allocatable :: str1
    character(len=3) :: schpg
    integer, allocatable :: nis(:)

    character*1, parameter :: lvecname(3) = (/"a","b","c"/)

    if (lcrys) then
       ! Header
       if (.not.c%ismolecule) then
          write (uout,'("* Crystal structure")')
          write (uout,'("  From: ",A)') string(c%file)
          write (uout,'("  Lattice parameters (bohr): ",3(A,2X))') &
             string(c%aa(1),'f',decimal=6), string(c%aa(2),'f',decimal=6), string(c%aa(3),'f',decimal=6)
          write (uout,'("  Lattice parameters (ang): ",3(A,2X))') &
             string(c%aa(1)*bohrtoa,'f',decimal=6), string(c%aa(2)*bohrtoa,'f',decimal=6), string(c%aa(3)*bohrtoa,'f',decimal=6)
          write (uout,'("  Lattice angles (degrees): ",3(A,2X))') &
             string(c%bb(1),'f',decimal=3), string(c%bb(2),'f',decimal=3), string(c%bb(3),'f',decimal=3)
       else
          write (uout,'("* Molecular structure")')
          write (uout,'("  From: ",A)') string(c%file)
          write (uout,'("  Encompassing cell dimensions (bohr): ",3(A,2X))') &
             string(c%aa(1),'f',decimal=6), string(c%aa(2),'f',decimal=6), string(c%aa(3),'f',decimal=6)
          write (uout,'("  Encompassing cell dimensions (ang): ",3(A,2X))') &
             string(c%aa(1)*bohrtoa,'f',decimal=6), string(c%aa(2)*bohrtoa,'f',decimal=6), string(c%aa(3)*bohrtoa,'f',decimal=6)
       endif

       ! Compute unit formula, and z
       allocate(nis(c%nspc))
       nis = 0
       do i = 1, c%nneq
          nis(c%at(i)%is) = nis(c%at(i)%is) + c%at(i)%mult
       end do
       maxdv = gcd(nis,c%nspc)
       write (uout,'("  Empirical formula: ",999(/4X,10(A,"(",A,") ")))') &
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
       do i = 1, c%nneq
          iz = c%spc(c%at(i)%is)%z
          if (iz >= maxzat) cycle
          nelec = nelec + iz * c%at(i)%mult
       end do
       write (uout,'("  Number of electrons (with zero atomic charge): ",A/)') string(nelec)
    end if

    if (lq) then
       write (uout,'("+ List of atomic species: ")')
       write (uout,'("# spc = atomic species. Z = atomic number. name = atomic name (symbol).")')
       write (uout,'("# Q = charge. ZPSP = pseudopotential charge.")')
       write (uout,'("# ",99(A,X))') string("spc",3,ioj_center), &
          string("Z",3,ioj_center), string("name",7,ioj_center),&
          string("Q",length=4,justify=ioj_center),&
          string("ZPSP",length=4,justify=ioj_right)
       do i = 1, c%nspc
          if (c%zpsp(c%spc(i)%z) > 0) then
             str1 = string(c%zpsp(c%spc(i)%z))
          else
             str1 = " -- "
          end if
          write (uout,'("  ",99(A,X))') string(i,3,ioj_center), &
             string(c%spc(i)%z,3,ioj_center), string(c%spc(i)%name,7,ioj_center),&
             string(c%spc(i)%qat,'f',length=4,decimal=1,justify=ioj_right),&
             str1
       end do
       write (uout,*)
    end if

    if (lcrys) then
       ! List of atoms in crystallographic coordinates
       if (.not.c%ismolecule) then
          write (uout,'("+ List of non-equivalent atoms in the unit cell (cryst. coords.): ")')
          write (uout,'("# at = complete list atomic ID. xyz = Cartesian coordinates. spc = atomic species.")')
          write (uout,'("# name = atomic name (symbol). mult = multiplicity. Z = atomic number.")')
          write (uout,'("# ",8(A,X))') string("nat",3,ioj_center), &
             string("x",14,ioj_center), string("y",14,ioj_center),&
             string("z",14,ioj_center), string("spc",3,ioj_center), string("name",7,ioj_center), &
             string("mult",3,ioj_center), string("Z",3,ioj_center)
          do i=1, c%nneq
             is = c%at(i)%is
             write (uout,'(2x,8(A,X))') &
                string(i,3,ioj_center),&
                string(c%at(i)%x(1),'f',length=14,decimal=10,justify=3),&
                string(c%at(i)%x(2),'f',length=14,decimal=10,justify=3),&
                string(c%at(i)%x(3),'f',length=14,decimal=10,justify=3),& 
                string(is,3,ioj_center), &
                string(c%spc(is)%name,7,ioj_center), &
                string(c%at(i)%mult,3,ioj_center), string(c%spc(is)%z,3,ioj_center)
          enddo
          write (uout,*)

          write (uout,'("+ List of atoms in the unit cell (cryst. coords.): ")')
          write (uout,'("# at = complete list atomic ID. xyz = Cartesian coordinates. spc = atomic species.")')
          write (uout,'("# name = atomic name (symbol). Z = atomic number.")')
          write (uout,'("# ",7(A,X))') string("at",3,ioj_center),&
             string("x",14,ioj_center), string("y",14,ioj_center),&
             string("z",14,ioj_center), string("spc",3,ioj_center), string("name",7,ioj_center),&
             string("Z",3,ioj_center)
          do i=1,c%ncel
             is = c%atcel(i)%is
             write (uout,'(2x,7(A,X))') &
                string(i,3,ioj_center),&
                string(c%atcel(i)%x(1),'f',length=14,decimal=10,justify=3),&
                string(c%atcel(i)%x(2),'f',length=14,decimal=10,justify=3),&
                string(c%atcel(i)%x(3),'f',length=14,decimal=10,justify=3),& 
                string(is,3,ioj_center),&
                string(c%spc(is)%name,7,ioj_center),&
                string(c%spc(is)%z,3,ioj_center)
          enddo
          write (uout,*)

          write (uout,'("+ Lattice vectors (",A,")")') iunitname0(iunit)
          do i = 1, 3
             write (uout,'(4X,A,": ",3(A,X))') lvecname(i), (string(c%m_x2c(j,i)*dunit0(iunit),'f',length=16,decimal=10,justify=5),j=1,3)
          end do
          write (uout,*)
       end if

       ! List of atoms in Cartesian coordinates
       write (uout,'("+ List of atoms in Cartesian coordinates (",A,"): ")') iunitname0(iunit)
       write (uout,'("# at = complete list atomic ID. xyz = Cartesian coordinates. spc = atomic species.")')
       write (uout,'("# name = atomic name (symbol). Z = atomic number. dnn = nearest-neighbor distance.")')
       write (uout,'("# ",99(A,X))') string("at",3,ioj_center), &
          string("x",16,ioj_center), string("y",16,ioj_center),&
          string("z",16,ioj_center), string("spc",3,ioj_center), string("name",7,ioj_center),&
          string("Z",3,ioj_center), string("dnn",10,ioj_center)
       do i=1,c%ncel
          is = c%atcel(i)%is
          write (uout,'(2x,99(A,X))') &
             string(i,3,ioj_center),&
             (string((c%atcel(i)%r(j)+c%molx0(j))*dunit0(iunit),'f',length=16,decimal=10,justify=5),j=1,3),&
             string(is,3,ioj_center),string(c%spc(is)%name,7,ioj_center),&
             string(c%spc(is)%z,3,ioj_center), string(2d0*c%at(c%atcel(i)%idx)%rnn2*dunit0(iunit),'f',length=10,decimal=4,justify=4)
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

       ! Cell volume
       if (.not.c%ismolecule) then
          write (uout,'("+ Cell volume (bohr^3): ",A)') string(c%omega,'f',decimal=5)
          write (uout,'("+ Cell volume (ang^3): ",A)') string(c%omega * bohrtoa**3,'f',decimal=5)
          write (uout,*)
       end if

       ! Write symmetry operations 
       if (.not.c%ismolecule) then
          write(uout,'("+ List of symmetry operations (",A,"):")') string(c%neqv)
          do k = 1, c%neqv
             write (uout,'(2X,"Operation ",A,":")') string(k)
             write (uout,'(2(4X,4(A,X)/),4X,4(A,X))') ((string(c%rotm(i,j,k),'f',length=9,decimal=6,justify=3), j = 1, 4), i = 1, 3)
          enddo
          write (uout,*)

          call c%struct_report_symxyz()

          write(uout,'("+ List of centering vectors (",A,"):")') string(c%ncv)
          do k = 1, c%ncv
             write (uout,'(2X,"Vector ",A,": ",3(A,X))') string(k), &
                (string(c%cen(i,k),'f',length=9,decimal=6), i = 1, 3)
          enddo
          write (uout,*)

          if (c%havesym > 0) then
             write(uout,'("+ Crystal symmetry information")')
             if (c%spg%n_atoms > 0) then
                if (len_trim(c%spg%choice) > 0) then
                   write(uout,'("  Space group (Hermann-Mauguin): ",A, " (number ",A,", setting ",A,")")') &
                      string(c%spg%international_symbol), string(c%spg%spacegroup_number),&
                      string(c%spg%choice)
                else
                   write(uout,'("  Space group (Hermann-Mauguin): ",A, " (number ",A,")")') &
                      string(c%spg%international_symbol), string(c%spg%spacegroup_number)
                end if
                write(uout,'("  Space group (Hall): ",A, " (number ",A,")")') &
                   string(c%spg%hall_symbol), string(c%spg%hall_number)
                write(uout,'("  Point group (Hermann-Mauguin): ",A)') string(c%spg%pointgroup_symbol)

                call pointgroup_info(c%spg%pointgroup_symbol,schpg,holo,laue)
                write(uout,'("  Point group (Schoenflies): ",A)') string(schpg)
                write(uout,'("  Holohedry: ",A)') string(holo_string(holo))
                write(uout,'("  Laue class: ",A)') string(laue_string(laue))
             else
                write(uout,'("  Unavailable because symmetry read from external file")')
             end if
             write (uout,*)
          end if

          write (uout,'("+ Cartesian/crystallographic coordinate transformation matrices:")')
          write (uout,'("  A = car to crys (xcrys = A * xcar, ",A,"^-1)")') iunitname0(iunit)
          do i = 1, 3
             write (uout,'(4X,3(A,X))') (string(c%m_c2x(i,j)/dunit0(iunit),'f',length=16,decimal=10,justify=5),j=1,3)
          end do
          write (uout,'("  B = crys to car (xcar = B * xcrys, ",A,")")') iunitname0(iunit)
          do i = 1, 3
             write (uout,'(4X,3(A,X))') (string(c%m_x2c(i,j)*dunit0(iunit),'f',length=16,decimal=10,justify=5),j=1,3)
          end do
          write (uout,'("  G = metric tensor (B''*B, ",A,"^2)")') iunitname0(iunit)
          do i = 1, 3
             write (uout,'(4X,3(A,X))') (string(c%gtensor(i,j)*dunit0(iunit)**2,'f',length=16,decimal=10,justify=5),j=1,3)
          end do
          write (uout,*)
       end if

       ! Discrete molecules, if available
       if (allocated(c%nstar) .and. allocated(c%mol) .and. c%nmol > 0) then
          write (uout,'("+ List of fragments in the system (",A,")")') string(c%nmol)
          write (uout,'("# Id = fragment ID. nat = number of atoms in fragment. C-o-m = center of mass (",A,").")') iunitname0(iunit)
          write (uout,'("# Discrete = is this fragment finite?")')
          write (uout,'("# Id  nat           Center of mass            Discrete  ")')
          do i = 1, c%nmol
             if (c%ismolecule) then
                xcm = (c%mol(i)%cmass()+c%molx0) * dunit0(iunit)
             else
                xcm = c%c2x(c%mol(i)%cmass())
             end if
             write (uout,'(99(2X,A))') string(i,3,ioj_left), string(c%mol(i)%nat,4,ioj_left),&
                (string(xcm(j),'f',10,6,3),j=1,3), string(c%mol(i)%discrete)
          end do
          if (.not.c%ismolecule) then
             if (all(c%mol(1:c%nmol)%discrete) .or. c%nlvac == 3) then
                write (uout,'(/"+ This is a molecular crystal.")')
             else if (c%nlvac == 2) then
                write (uout,'(/"+ This is a 1D periodic (polymer) structure.")')
                write (uout,'("  Vacuum lattice vectors: (",2(A,X),A,"), (",2(A,X),A,")")') &
                   (string(c%lvac(j,1)),j=1,3), (string(c%lvac(j,2)),j=1,3)
                write (uout,'("  Connected lattice vectors: (",2(A,X),A,")")') &
                   (string(c%lcon(j,1)),j=1,3)
             else if (c%nlvac == 1) then
                write (uout,'(/"+ This is a 2D periodic (layered) structure.")')
                write (uout,'("  Vacuum lattice vectors: (",2(A,X),A,")")') &
                   (string(c%lvac(j,1)),j=1,3)
                write (uout,'("  Connected lattice vectors: (",2(A,X),A,"), (",2(A,X),A,")")') &
                   (string(c%lcon(j,1)),j=1,3), (string(c%lcon(j,2)),j=1,3)
             else 
                write (uout,'(/"+ This is a 3D periodic structure.")')
             end if
             write (uout,*)
          end if
       end if

       ! Number of atoms in the atomic environment
       call c%env%report()

       ! Wigner-Seitz cell
       if (.not.c%ismolecule) then
          write (uout,'("+ Vertex of the WS cell in cryst. coords. (",A,")")') string(c%ws_nv)
          write (uout,'("# id = vertex ID. xyz = vertex cryst. coords. d = vertex distance to origin (",A,").")') iunitname0(iunit)
          write (uout,'(5(2X,A))') string("id",length=3,justify=ioj_right),&
             string("x",length=11,justify=ioj_center),&
             string("y",length=11,justify=ioj_center),&
             string("z",length=11,justify=ioj_center),&
             string("d ("//iunitname0(iunit)//")",length=14,justify=ioj_center)
          do i = 1, c%ws_nv
             x0 = c%x2c(c%ws_x(:,i))
             write (uout,'(5(2X,A))') string(i,length=3,justify=ioj_right), &
                (string(c%ws_x(j,i),'f',length=11,decimal=6,justify=4),j=1,3), &
                string(norm2(x0)*dunit0(iunit),'f',length=14,decimal=8,justify=4)
          enddo
          write (uout,*)

          write (uout,'("+ Faces of the WS cell (",A,")")') string(c%ws_nf)
          write (uout,'("# Face ID: vertexID1 vertexID2 ...")')
          do i = 1, c%ws_nf
             write (uout,'(2X,A,": ",999(A,X))') string(i,length=2,justify=ioj_right), &
                (string(c%ws_iside(j,i),length=2),j=1,c%ws_nside(i))
          end do
          write (uout,*)

          write (uout,'("+ Lattice vectors for the Wigner-Seitz neighbors")')
          write (uout,'("# FaceID: Voronoi lattice vector (cryst. coords.)")')
          do i = 1, c%ws_nf
             write (uout,'(2X,A,": ",99(A,X))') string(i,length=2,justify=ioj_right), &
                (string(c%ws_ineighx(j,i),length=2,justify=ioj_right),j=1,3)
          end do
          write (uout,*)

          write (uout,'("+ Lattice vectors for the Delaunay reduced cell (cryst. coords.)")')
          do i = 1, 3
             write (uout,'(2X,A,": ",99(A,X))') lvecname(i), &
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

          write (uout,'("  Delaunay reduced cell lengths: ",99(A,X))') &
             (string(xlen(j),'f',decimal=6,justify=ioj_right),j=1,3)
          write (uout,'("  Delaunay reduced cell angles: ",99(A,X))') &
             (string(xang(j),'f',decimal=3,justify=ioj_right),j=1,3)
          write (uout,*)

          write (uout,'("+ Is the cell orthogonal? ",L1)') c%isortho
          write (uout,'("+ Is the reduced cell orthogonal? ",L1/)') c%isortho_del
       end if
    end if

  end subroutine struct_report

  !> Write the list of symmetry operations to stdout, using crystallographic
  !> notation (if possible).
  module subroutine struct_report_symxyz(c,strfin)
    use tools_io, only: uout, string
    use global, only: symprec
    use param, only: mlen
    class(crystal), intent(in) :: c
    character(len=mlen), intent(out), optional :: strfin(c%neqv)

    real*8, parameter :: rfrac(25) = (/-12d0/12d0,-11d0/12d0,-10d0/12d0,&
       -9d0/12d0,-8d0/12d0,-7d0/12d0,-6d0/12d0,-5d0/12d0,-4d0/12d0,-3d0/12d0,&
       -2d0/12d0,-1d0/12d0,0d0/12d0,1d0/12d0,2d0/12d0,3d0/12d0,4d0/12d0,&
       5d0/12d0,6d0/12d0,7d0/12d0,8d0/12d0,9d0/12d0,10d0/12d0,11d0/12d0,12d0/12d0/)
    character*6, parameter :: sfrac(25) = (/"      ","-11/12","-5/6  ",&
       "-3/4  ","-2/3  ","-7/12 ","-1/2  ","-5/12 ","-1/3  ","-1/4  ","-1/6  ",&
       "-1/12 ","      ","1/12  ","1/6   ","1/4   ","1/3   ","5/12  ","1/2   ",&
       "7/12  ","2/3   ","3/4   ","5/6   ","11/12 ","      "/)
    character*1, parameter :: xyz(3) = (/"x","y","z"/)

    logical :: ok, iszero
    integer :: i, j, k
    character*255 :: strout(c%neqv)

    do i = 1, c%neqv
       strout(i) = ""
       do j = 1, 3
          ! translation
          ok = .false.
          do k = 1, 25
             if (abs(c%rotm(j,4,i) - rfrac(k)) < symprec) then
                ok = .true.
                strout(i) = trim(strout(i)) // sfrac(k)
                iszero = (k == 13) .or. (k == 1) .or. (k == 25)
                exit
             end if
          end do
          if (.not.ok) return

          ! rotation
          do k = 1, 3
             if (abs(c%rotm(j,k,i) - 1d0) < symprec) then
                if (iszero) then
                   strout(i) = trim(strout(i)) // xyz(k)
                else
                   strout(i) = trim(strout(i)) // "+" // xyz(k)
                end if
                iszero = .false.
             elseif (abs(c%rotm(j,k,i) + 1d0) < symprec) then
                strout(i) = trim(strout(i)) // "-" // xyz(k)
                iszero = .false.
             elseif (abs(c%rotm(j,k,i)) > symprec) then
                return
             end if
          end do

          ! the comma
          if (j < 3) &
             strout(i) = trim(strout(i)) // ","
       end do
    end do

    if (present(strfin)) then
       strfin = strout
    else
       write(uout,'("+ List of symmetry operations in crystallographic notation:")')
       do k = 1, c%neqv
          write (uout,'(3X,A,": ",A)') string(k), string(strout(k))
       enddo
       write (uout,*)
    end if

  end subroutine struct_report_symxyz

  !> Use the spg library to find information about the space group.
  !> In: cell vectors (m_x2c), ncel, atcel(:), at(:) Out: neqv, rotm,
  !> spg, at()%mult, ncel, and atcel(). If usenneq is .true., use nneq
  !> and at(:) instead of ncel and atcel. If onlyspg is .true., fill
  !> only the spg field and leave the others unchanged.
  module subroutine spglib_wrap(c,usenneq,onlyspg)
    use iso_c_binding, only: c_double
    use spglib, only: spg_get_dataset, spg_get_error_message
    use global, only: symprec, atomeps
    use tools_io, only: string, ferror, equal, faterr
    use param, only: maxzat0, eyet, eye
    use types, only: realloc
    class(crystal), intent(inout) :: c
    logical, intent(in) :: usenneq
    logical, intent(in) :: onlyspg

    real(c_double) :: lattice(3,3)
    real(c_double), allocatable :: x(:,:)
    integer, allocatable :: typ(:), iidx(:)
    integer :: ntyp, nat
    integer :: i, j, k, iat, iz(maxzat0), idx
    character(len=32) :: error
    logical :: found
    real*8 :: rotm(3,3), x0(3)
    logical, allocatable :: used(:)

    ! get the dataset from spglib
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
    c%spg = spg_get_dataset(lattice,x,typ,nat,symprec)
    deallocate(x,typ)

    ! check error messages
    error = trim(spg_get_error_message(c%spg%spglib_error))
    if (.not.equal(error,"no error")) &
       call ferror("spglib_wrap","error from spglib: "//string(error),faterr)

    if (onlyspg) return

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
       else
          c%at(iidx(idx))%mult = c%at(iidx(idx))%mult + 1
       end if
       c%atcel(i)%idx = iidx(idx)
    end do
    deallocate(iidx)
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
    allocate(used(c%ncel))
    used = .false.
    main: do iat = 1, c%nneq
       do i = 1, c%neqv
          do j = 1, c%ncv
             x0 = matmul(c%rotm(1:3,1:3,i),c%at(iat)%x) + c%rotm(:,4,i) + c%cen(:,j)
             found = .false.
             do k = 1, c%ncel
                if (used(k)) cycle
                if (c%are_lclose(x0,c%atcel(k)%x,atomeps)) then
                   c%atcel(k)%ir = i
                   c%atcel(k)%ic = j
                   c%atcel(k)%lvec = nint(c%atcel(k)%x - x0)
                   used(k) = .true.
                   exit
                end if
             end do
             if (all(used)) exit main
          end do
       end do
    end do main
    if (.not.all(used)) &
       call ferror("spglib_wrap","error building rotation and center information for atom list",faterr)

  end subroutine spglib_wrap

  !xx! Wigner-Seitz cell tools and cell partition

  !> Builds the Wigner-Seitz cell. Writes the WS components of c.
  !> Also writes the Delaunay reduction parameters and, if area is
  !> present, calculates the areas of the WS facets. Sets the matrices
  !> for the Delaunay transformations.
  module subroutine wigner(c,area)
    use, intrinsic :: iso_c_binding, only: c_char, c_null_char, c_int, c_double
    use tools_math, only: mixed, cross, matinv, mnorm2, det
    use tools_io, only: string, fopen_write, fopen_read,&
       ferror, faterr, fclose
    use types, only: realloc
    use param, only: pi, eye
    class(crystal), intent(inout) :: c
    real*8, intent(out), optional :: area(14) !< area of the WS faces 

    interface
       ! The definitions and documentation for these functions are in doqhull.c
       subroutine runqhull_voronoi_step1(n,xstar,nf,nv,mnfv) bind(c)
         use, intrinsic :: iso_c_binding, only: c_int, c_double
         integer(c_int), value :: n
         real(c_double) :: xstar(3,n)
         integer(c_int) :: nf, nv, mnfv
       end subroutine runqhull_voronoi_step1
       subroutine runqhull_voronoi_step2(nf,nv,mnfv,ivws,xvws,nfvws,fvws) bind(c)
         use, intrinsic :: iso_c_binding, only: c_int, c_double
         integer(c_int), value :: nf, nv, mnfv
         integer(c_int) :: ivws(nf)
         real(c_double) :: xvws(3,nv)
         integer(c_int) :: nfvws(mnfv)
         integer(c_int) :: fvws(mnfv)
       end subroutine runqhull_voronoi_step2
    end interface

    real*8, parameter :: eps_dnorm = 1d-5 !< minimum lattice vector length

    integer :: i, j, k
    real*8 :: av(3), bary(3), rmat(3,4)
    integer(c_int) :: n
    real(c_double) :: xstar(3,14)
    integer(c_int), allocatable :: ivws(:)
    real(c_double), allocatable :: xvws(:,:)
    real*8 :: rdel(3,3)
    real*8 :: xlen(3), xang(3), x0(3), xred(3,3)

    ! delaunay reduction
    call c%delaunay_reduction(rmat,rbas=rdel)

    ! construct star of lattice vectors -> use Delaunay reduction
    ! see 9.1.8 in ITC.
    n = 14
    xstar(:,1)  = rmat(:,1)
    xstar(:,2)  = rmat(:,2)
    xstar(:,3)  = rmat(:,3)
    xstar(:,4)  = rmat(:,4)
    xstar(:,5)  = rmat(:,1)+rmat(:,2)
    xstar(:,6)  = rmat(:,1)+rmat(:,3)
    xstar(:,7)  = rmat(:,2)+rmat(:,3)
    xstar(:,8)  = -(rmat(:,1))
    xstar(:,9)  = -(rmat(:,2))
    xstar(:,10) = -(rmat(:,3))
    xstar(:,11) = -(rmat(:,4))
    xstar(:,12) = -(rmat(:,1)+rmat(:,2))
    xstar(:,13) = -(rmat(:,1)+rmat(:,3))
    xstar(:,14) = -(rmat(:,2)+rmat(:,3))
    do i = 1, 14
       xstar(:,i) = matmul(c%m_x2c,xstar(:,i))
       if (norm2(xstar(:,i)) < eps_dnorm) &
          call ferror("wigner","Lattice vector too short. Please, check the unit cell definition.",faterr)
    end do

    call runqhull_voronoi_step1(n,xstar,c%ws_nf,c%ws_nv,c%ws_mnfv)
    allocate(ivws(c%ws_nf),c%ws_iside(c%ws_mnfv,c%ws_nf),xvws(3,c%ws_nv))
    c%ws_nside = 0
    call runqhull_voronoi_step2(c%ws_nf,c%ws_nv,c%ws_mnfv,ivws,xvws,c%ws_nside(1:c%ws_nf),c%ws_iside)

    if (allocated(c%ws_x)) deallocate(c%ws_x)
    allocate(c%ws_x(3,c%ws_nv))
    do i = 1, c%ws_nv
       c%ws_x(:,i) = matmul(c%m_c2x,xvws(:,i))
    end do

    c%ws_ineighc = 0
    c%ws_ineighx = 0
    do i = 1, c%ws_nf
       c%ws_ineighc(:,i) = xstar(:,ivws(i))
       c%ws_ineighx(:,i) = nint(matmul(c%m_c2x,xstar(:,ivws(i))))
    end do
    if (present(area)) then
       do i = 1, c%ws_nf
          ! lattice point
          bary = 0d0
          do j = 1, c%ws_nside(i)
             bary = bary + xvws(:,c%ws_iside(j,i))
          end do
          bary = 2d0 * bary / c%ws_nside(i)

          ! area of a convex polygon
          av = 0d0
          do j = 1, c%ws_nside(i)
             k = mod(j,c%ws_nside(i))+1
             av = av + cross(xvws(:,c%ws_iside(j,i)),xvws(:,c%ws_iside(k,i)))
          end do
          area(i) = 0.5d0 * abs(dot_product(bary,av) / norm2(bary))
       end do
    end if

    deallocate(ivws,xvws)
    
    ! is the cell orthogonal?
    c%isortho = (c%ws_nf <= 6)
    if (c%isortho) then
       do i = 1, c%ws_nf
          c%isortho = c%isortho .and. (count(abs(c%ws_ineighx(:,i)) == 1) == 1) .and.&
             (count(abs(c%ws_ineighx(:,i)) == 0) == 2)
       end do
    end if

    ! calculate the delaunay reduction parameters for shortest vector search
    if (c%ismolecule) then
       c%m_xr2x = eye
       c%m_x2xr = eye
    else
       c%m_xr2x = rdel
       c%m_x2xr = matinv(c%m_xr2x)
    end if
    c%m_xr2c = matmul(c%m_x2c,c%m_xr2x)
    c%m_c2xr = matinv(c%m_xr2c)

    do i = 1, 3
       x0 = c%m_xr2x(i,:)
       xred(:,i) = c%x2c(x0)
       xlen(i) = norm2(xred(:,i))
    end do
    xang(1) = acos(dot_product(xred(:,2),xred(:,3)) / xlen(2) / xlen(3)) * 180d0 / pi
    xang(2) = acos(dot_product(xred(:,1),xred(:,3)) / xlen(1) / xlen(3)) * 180d0 / pi
    xang(3) = acos(dot_product(xred(:,1),xred(:,2)) / xlen(1) / xlen(2)) * 180d0 / pi

    c%n2_xr2x = mnorm2(c%m_xr2x)
    c%n2_x2xr = mnorm2(c%m_x2xr)
    c%n2_xr2c = mnorm2(c%m_xr2c)
    c%n2_c2xr = mnorm2(c%m_c2xr)

    do i = 1, c%ws_nf
       c%ws_ineighxr(:,i) = nint(matmul(c%ws_ineighx(:,i),c%m_x2xr))
    end do

    c%isortho_del = .true.
    do i = 1, c%ws_nf
       c%isortho_del = c%isortho_del .and. (count(abs(c%ws_ineighxr(:,i)) == 1) == 1) .and.&
          (count(abs(c%ws_ineighxr(:,i)) == 0) == 2)
    end do

  end subroutine wigner

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
  module subroutine write_mol(c,file,fmt,ix,doborder,onemotif,molmotif,&
     environ,renv,lnmer,nmer,rsph,xsph,rcub,xcub,luout)
    use global, only: dunit0, iunit
    use tools_math, only: nchoosek, comb
    use tools_io, only: ferror, faterr, uout, string, ioj_left, string, ioj_right,&
       equal
    use types, only: realloc
    use fragmentmod, only: fragment, realloc_fragment
    class(crystal), intent(inout) :: c
    character*(*), intent(in) :: file
    character*3, intent(in) :: fmt
    integer, intent(in) :: ix(3)
    logical, intent(in) :: doborder, onemotif, molmotif, environ
    real*8, intent(in) :: renv
    logical, intent(in) :: lnmer
    integer, intent(in) :: nmer
    real*8, intent(in) :: rsph, xsph(3)
    real*8, intent(in) :: rcub, xcub(3)
    integer, intent(out), optional :: luout

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
          write (uout,'(99(2X,A))') string(i,3,ioj_left), string(fr0(i)%nat,4,ioj_left),&
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
       call dowrite(file,fr)
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
                call dowrite(file0,fr0(j))
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
                   call dowrite(file0,fr)
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
                call dowrite(file0,fr)
                icount = icount + 1
             end do
             deallocate(icomb)
             write (uout,'("+ Written ",A," ",A,"-mers")') string(icount), string(i)
          end if
       end do
    end if

  contains
    subroutine dowrite(fileo,fro)
      character*(*) :: fileo
      type(fragment) :: fro
      
      if (equal(fmt,"xyz")) then
         call fro%writexyz(fileo)
      elseif (equal(fmt,"gjf")) then
         call fro%writegjf(fileo)
      elseif (equal(fmt,"cml")) then
         if (c%ismolecule) then
            call fro%writecml(fileo,luout=luout)
         else
            call fro%writecml(fileo,c%m_x2c,luout=luout)
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
  module subroutine write_3dmodel(c,file,fmt,ix,doborder,onemotif,molmotif,&
     docell,domolcell,rsph,xsph,rcub,xcub,gr0)
    use graphics, only: grhandle
    use fragmentmod, only: fragment
    use tools_io, only: equal
    use param, only: maxzat, atmcov, jmlcol
    class(crystal), intent(inout) :: c
    character*(*), intent(in) :: file
    character*3, intent(in) :: fmt
    integer, intent(in) :: ix(3)
    logical, intent(in) :: doborder, onemotif, molmotif
    logical, intent(in) :: docell, domolcell
    real*8, intent(in) :: rsph, xsph(3)
    real*8, intent(in) :: rcub, xcub(3)
    type(grhandle), intent(out), optional :: gr0

    integer :: i, j
    real*8 :: d, xd(3), x0(3), x1(3), rr
    type(fragment) :: fr
    type(fragment), allocatable :: fr0(:)
    logical, allocatable :: isdiscrete(:)
    integer :: nmol
    type(grhandle) :: gr

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

    call gr%open(fmt,file)

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
       call gr%close()
    end if

    if (allocated(fr0)) deallocate(fr0)

  end subroutine write_3dmodel

  !> Write a quantum espresso input template
  module subroutine write_espresso(c,file)
    use tools_io, only: fopen_write, lower, fclose
    use param, only: atmass
    class(crystal), intent(in) :: c
    character*(*), intent(in) :: file

    integer :: i, lu

    lu = fopen_write(file)
    write (lu,'("&control")')
    write (lu,'(" title=''crystal'',")')
    write (lu,'(" prefix=''crystal'',")')
    write (lu,'(" pseudo_dir=''.'',")')
    write (lu,'(" calculation=''vc-relax'',")')
    write (lu,'("/")')
    write (lu,'("&system")')
    write (lu,'(" ibrav=0,")')
    write (lu,'(" celldm(1)=1.0,")')
    write (lu,'(" nat=",I6,",")') c%ncel
    write (lu,'(" ntyp=",I3,",")') c%nspc
    write (lu,'(" ecutwfc=60.0,")')
    write (lu,'(" ecutrho=600.0,")')
    write (lu,'(" xdm=.true.,")')
    write (lu,'(" xdm_a1=0.4073,")')
    write (lu,'(" xdm_a2=2.4150,")')
    write (lu,'("/")')
    write (lu,'("&electrons"/" conv_thr = 1d-8,"/"/")')
    write (lu,'("&ions"/"/")')
    write (lu,'("&cell"/"/")')
    write (lu,'("ATOMIC_SPECIES")')
    do i = 1, c%nspc
       write (lu,'(A,X,F12.6,X,A,".UPF")') trim(c%spc(i)%name), atmass(c%spc(i)%z), trim(lower(c%spc(i)%name))
    end do
    write (lu,'(/"ATOMIC_POSITIONS crystal")')
    do i = 1, c%ncel
       write (lu,'(A,3(X,F13.8,X))') trim(c%spc(c%atcel(i)%is)%name), c%atcel(i)%x
    end do
    write (lu,'(/"K_POINTS automatic"/"2 2 2 1 1 1"/)')
    write (lu,'("CELL_PARAMETERS")')
    do i = 1, 3
       write (lu,'(3(F18.12,X))') c%m_x2c(:,i)
    end do
    call fclose(lu)

  end subroutine write_espresso

  !> Write a VASP POSCAR template
  module subroutine write_vasp(c,file,verbose)
    use tools_io, only: fopen_write, string, uout, fclose
    use param, only: bohrtoa
    class(crystal), intent(in) :: c
    character*(*), intent(in) :: file
    logical, intent(in) :: verbose

    character(len=:), allocatable :: lbl1, aux
    integer :: i, j, lu, ntyp

    ! Cell
    lu = fopen_write(file)
    write (lu,'("crystal")')
    write (lu,'("1.0")')
    do i = 1, 3
       write (lu,'(3(F15.10,X))') c%m_x2c(:,i) * bohrtoa
    end do

    ! Number of atoms per type and Direct
    lbl1 = ""
    do i = 1, c%nspc
       ntyp = 0
       do j = 1, c%ncel
          if (c%atcel(j)%is == i) ntyp = ntyp + 1
       end do
       aux = lbl1 // " " // string(ntyp)
       lbl1 = aux
    end do
    write (lu,'(A)') lbl1
    write (lu,'("Direct")')

    ! Atomic positions
    do i = 1, c%nspc
       do j = 1, c%ncel
          if (c%atcel(j)%is == i) then
             write (lu,'(3(F13.8,X))') c%atcel(j)%x
          end if
       end do
    end do
    call fclose(lu)

    if (verbose) &
       write (uout,'("+ Atom type sequence: ",999(A,X))') (string(c%spc(j)%name),j=1,c%nspc)

  end subroutine write_vasp

  !> Write an abinit input template
  module subroutine write_abinit(c,file)
    use tools_io, only: fopen_write, string, fclose
    use param, only: pi
    class(crystal), intent(in) :: c
    character*(*), intent(in) :: file

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
    lu = fopen_write(file)
    write (lu,'("acell ",3(F14.10,X))') aap
    write (lu,'("angdeg ",3(F14.10,X))') bbp
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
             write (lu,'(X,3(F15.10,X))') c%atcel(j)%x
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
  module subroutine write_elk(c,file)
    use tools_io, only: fopen_write, fclose
    class(crystal), intent(in) :: c
    character*(*), intent(in) :: file

    integer :: ntyp
    integer :: i, j, lu

    ! Write input
    lu = fopen_write(file)
    write (lu,'("tasks"/,"0"/)')
    write (lu,'("xctype"/,"20"/)')
    write (lu,'("avec")')
    do i = 1, 3
       write (lu,'(2X,3(F15.10,X))') c%m_x2c(:,i)
    end do
    write (lu,*)

    write (lu,'("sppath"/,"''./''"/)')

    write (lu,'("atoms")')
    write (lu,'(2X,I4)') c%nspc
    do i = 1, c%nspc
       write (lu,'(2X,"''",A,".in''")') trim(c%spc(i)%name)
       ntyp = 0
       do j = 1, c%ncel
          if (c%atcel(j)%is == i) ntyp = ntyp + 1
       end do
       write (lu,'(2X,I3)') ntyp
       do j = 1, c%ncel
          if (c%atcel(j)%is == i) then
             write (lu,'(2X,3(F14.10,X),"0.0 0.0 0.0")') c%atcel(j)%x
          end if
       end do
    end do
    write (lu,*)

    write (lu,'("ngridk"/,"  4 4 4"/)')
    write (lu,'("rgkmax"/,"  7.0"/)')
    call fclose(lu)

  end subroutine write_elk

  !> Write a Gaussian template input (periodic).
  module subroutine write_gaussian(c,file)
    use tools_io, only: fopen_write, string, nameguess, ioj_left, fclose
    use param, only: bohrtoa
    class(crystal), intent(in) :: c
    character*(*), intent(in) :: file

    character(len=:), allocatable :: wroot
    integer :: lu, i, j

    wroot = file(:index(file,'.',.true.)-1)

    lu = fopen_write(file)
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
       write (lu,'(99(A,X))') string(nameguess(c%spc(c%atcel(i)%is)%z,.true.),2,ioj_left),&
          (string(c%atcel(i)%r(j)*bohrtoa,'f',14,8,ioj_left),j=1,3)
    end do
    do i = 1, 3
       write (lu,'(99(A,X))') string("Tv",2,ioj_left),&
          (string(c%m_x2c(j,i)*bohrtoa,'f',14,8,ioj_left),j=1,3)
    end do
    write (lu,*)

    call fclose(lu)

  end subroutine write_gaussian

  !> Write a tessel input template
  module subroutine write_tessel(c,file)
    use global, only: fileroot
    use tools_io, only: fopen_write, fclose
    class(crystal), intent(in) :: c
    character*(*), intent(in) :: file

    integer :: lu, i

    lu = fopen_write(file)

    write (lu,'("set camangle 75 -10 45")')
    write (lu,'("set background background {color rgb <1,1,1>}")')
    write (lu,'("set use_planes .false.")')
    write (lu,'("set ball_texture finish{specular 0.2 roughness 0.1 reflection 0.1}")')
    write (lu,'("set equalscale noscale")')
    write (lu,'("molecule")')
    write (lu,'("  crystal")')
    write (lu,'("    symmatrix seitz")')
    do i = 1, c%ncv
       write (lu,'(5X,A,3(F15.12,X))') "cen ",c%cen(:,i)
    end do
    write (lu,'(5X,"#")')
    do i = 1, c%neqv
       write (lu,'(5X,3(F5.2,X),F15.12)') c%rotm(1,:,i)
       write (lu,'(5X,3(F5.2,X),F15.12)') c%rotm(2,:,i)
       write (lu,'(5X,3(F5.2,X),F15.12)') c%rotm(3,:,i)
       write (lu,'(5X,"#")')
    end do
    write (lu,'(5X,"endsymmatrix")')
    write (lu,'(5X,A,6(F12.8,X))') "cell", c%aa, c%bb
    write (lu,'(5X,"crystalbox  -2.30 -2.30 -2.30 2.30 2.30 2.30")')
    write (lu,'(5X,A,6(F6.3,X))') "clippingbox ",-0.02,-0.02,-0.02,+1.02,+1.02,+1.02
    do i = 1, c%nneq
       write (lu,'(5X,"neq ",3(F12.8," "),A10)') c%at(i)%x, trim(c%spc(c%at(i)%is)%name)
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
  module subroutine write_critic(c,file)
    use tools_io, only: fopen_write, fclose, string
    class(crystal), intent(in) :: c
    character*(*), intent(in) :: file

    integer :: lu, i, j

    lu = fopen_write(file)

    write (lu,'("crystal")')
    write (lu,'("  cell ",6(A,X))') (string(c%aa(i),'f',20,10),i=1,3),&
       (string(c%bb(i),'f',20,10),i=1,3)
    do i = 1, c%ncel
       write (lu,'("  neq ",3(A," "),A10)') (string(c%atcel(i)%x(j),'f',20,10),j=1,3),&
          string(c%spc(c%atcel(i)%is)%name)
    end do
    write (lu,'("endcrystal")')
    write (lu,'("end")')
    call fclose(lu)

  end subroutine write_critic

  !> Write a simple cif file
  module subroutine write_cif(c,file)
    use global, only: fileroot
    use tools_io, only: fopen_write, fclose, string, nameguess
    use param, only: bohrtoa
    class(crystal), intent(in) :: c
    character*(*), intent(in) :: file

    integer :: i, iz, lu

    lu = fopen_write(file)

    write (lu,'("data_",A)') string(fileroot)
    write (lu,'("_cell_volume ",F20.6)') c%omega * bohrtoa**3
    write (lu,'("_symmetry_space_group_name_H-M ''P 1''")');
    write (lu,'("_symmetry_Int_Tables_number 1")');
    write (lu,'("loop_")');
    write (lu,'("_symmetry_equiv_pos_site_id")');
    write (lu,'("_symmetry_equiv_pos_as_xyz")');
    write (lu,'("1 x,y,z")');
    write (lu,'("_cell_length_a ",F20.10)') c%aa(1)*bohrtoa
    write (lu,'("_cell_length_b ",F20.10)') c%aa(2)*bohrtoa
    write (lu,'("_cell_length_c ",F20.10)') c%aa(3)*bohrtoa
    write (lu,'("_cell_angle_alpha ",F14.4)') c%bb(1)
    write (lu,'("_cell_angle_beta ",F14.4)') c%bb(2)
    write (lu,'("_cell_angle_gamma ",F14.4)') c%bb(3)
    write (lu,'("_cell_formula_units_Z 1")')
    write (lu,'("loop_")');
    write (lu,'("_atom_site_label")');
    write (lu,'("_atom_site_type_symbol")');
    write (lu,'("_atom_site_fract_x")');
    write (lu,'("_atom_site_fract_y")');
    write (lu,'("_atom_site_fract_z")');
    do i = 1, c%ncel
       iz = c%atcel(i)%is
       write (lu,'(A5,X,A3,X,3(F20.14,X))') c%spc(iz)%name, &
          nameguess(c%spc(iz)%z,.true.), c%atcel(i)%x
    end do
    call fclose(lu)

  end subroutine write_cif

  !> Write a simple cif file
  module subroutine write_d12(c,file,dosym)
    use tools_io, only: fopen_write, fclose, string, ferror, faterr
    use global, only: symprec
    use param, only: bohrtoa, pi
    class(crystal), intent(in) :: c
    character*(*), intent(in) :: file
    logical, intent(in) :: dosym

    integer :: holo, laue, i, j, nmin
    character(len=3) :: schpg
    type(crystal) :: nc
    integer :: lu, spgnum, irhomb, count90, count120
    real*8 :: xmin(6)
    logical :: ok

    ! we need symmetry for this
    if (dosym) then
       ! This is to address a bug in gfortran 4.9 (and possibly earlier versions)
       ! regarding assignment of user-defined types with allocatable components. 
       nc = c
       call nc%cell_standard(.false.,.false.,.false.)
       call pointgroup_info(nc%spg%pointgroup_symbol,schpg,holo,laue)
       xmin = 0d0
       irhomb = 0
       if (holo == holo_unk) then
          call ferror("write_d12","unknown holohedry",faterr)
       elseif (holo == holo_tric) then
          nmin = 6
          xmin(1:3) = nc%aa * bohrtoa
          xmin(4:6) = nc%bb
       elseif (holo == holo_mono) then
          nmin = 4
          xmin(1:3) = nc%aa * bohrtoa
          ok = .false.
          do i = 1, 3
             if (abs(nc%bb(i) - 90d0) > symprec) then
                xmin(4) = nc%bb(i)
                ok = .true.
                exit
             end if
          end do
          if (.not.ok) &
             xmin(4) = 90d0
       elseif (holo == holo_ortho) then
          nmin = 3
          xmin(1:3) = nc%aa * bohrtoa
       elseif (holo == holo_tetra) then
          nmin = 2
          xmin(1) = nc%aa(1) * bohrtoa
          xmin(2) = nc%aa(3) * bohrtoa
       elseif (holo == holo_trig) then
          count90 = count(abs(nc%bb - 90d0) < 1d-1)
          count120 = count(abs(sin(nc%bb * pi / 180d0) - sqrt(3d0)/2d0) < 1d-2)
          nmin = 2
          xmin(1) = nc%aa(1) * bohrtoa
          if (count90 == 2 .and. count120 == 1) then
             ! hexagonal axes
             xmin(2) = nc%aa(3) * bohrtoa
             irhomb = 0
          else
             ! rhombohedral axes
             xmin(2) = nc%bb(1)
             irhomb = 1
          end if
       elseif (holo == holo_hex) then
          nmin = 2
          xmin(1) = nc%aa(1) * bohrtoa
          xmin(2) = nc%aa(3) * bohrtoa
       elseif (holo == holo_cub) then
          nmin = 1
          xmin(1) = nc%aa(1) * bohrtoa
       end if
       spgnum = nc%spg%spacegroup_number
    else
       spgnum = 1
       irhomb = 0
       nmin = 6
       xmin(1:3) = c%aa * bohrtoa
       xmin(4:6) = c%bb
    end if

    lu = fopen_write(file)
    write (lu,'("Title")')
    write (lu,'("CRYSTAL")')
    write (lu,'("0 ",A," 0")') string(irhomb)
    write (lu,'(A)') string(spgnum)
    write (lu,'(6(A,X))') (string(xmin(i),'f',15,8),i=1,nmin)
    if (dosym) then
       write (lu,'(A)') string(nc%nneq)
       do i = 1, nc%nneq
          write (lu,'(4(A,X))') string(nc%spc(nc%at(i)%is)%z), (string(nc%at(i)%x(j),'f',15,8),j=1,3)
       end do
    else
       write (lu,'(A)') string(c%ncel)
       do i = 1, c%ncel
          write (lu,'(4(A,X))') string(c%spc(c%atcel(i)%is)%z), (string(c%atcel(i)%x(j),'f',15,8),j=1,3)
       end do
    end if
    write (lu,'("SETPRINT")')
    write (lu,'("1")')
    write (lu,'("3 1")')
    write (lu,'("END")')
    write (lu,'("xx basis xx")')
    write (lu,'("99 0")')
    write (lu,'("END")')
    write (lu,'("SHRINK")')
    write (lu,'("4 4")')
    write (lu,'("TOLDEE")')
    write (lu,'("7")')
    write (lu,'("END")')

    call fclose(lu)

  end subroutine write_d12

  !> Write an escher octave script
  module subroutine write_escher(c,file)
    use global, only: fileroot
    use tools_io, only: fopen_write, string, fclose
    use param, only: pi
    class(crystal), intent(in) :: c
    character*(*), intent(in) :: file

    character(len=:), allocatable :: lbl1, aux
    integer :: lu, i, n

    lu = fopen_write(file)

    ! count number of atoms per type
    write (lu,'("cr = struct();")')
    write (lu,'("cr.name = """,A,""";")') trim(adjustl(fileroot))
    write (lu,'("cr.a = [",1p,3(E22.14,X),"];")') c%aa
    write (lu,'("cr.b = [",1p,3(E22.14,X),"];")') c%bb * pi / 180d0
    write (lu,'("cr.nat = ",I6,";")') c%ncel
    write (lu,'("cr.ntyp = ",I6,";")') c%nspc
    write (lu,'("cr.r = [")')
    do i = 1, 3
       write (lu,'(2X,1p,3(E22.14,X))') c%m_x2c(:,i)
    end do
    write (lu,'(2X,"];")')
    write (lu,'("cr.g = [")')
    do i = 1, 3
       write (lu,'(2X,1p,3(E22.14,X))') c%gtensor(:,i)
    end do
    write (lu,'(2X,"];")')
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
       write (lu,'(2X,1p,3(E22.14,X))') c%atcel(i)%x
    end do
    write (lu,'("  ];")')

    call fclose(lu)

  end subroutine write_escher

  !> Write a db file for the dcp automatic input generator
  module subroutine write_db(c,file)
    use tools_io, only: fopen_write, string, fclose, nameguess
    use param, only: bohrtoa
    class(crystal), intent(in) :: c
    character*(*), intent(in) :: file

    integer :: lu, i, j

    lu = fopen_write(file)
    write (lu,'("type crystal_energy")')
    write (lu,'("kpts 4")')
    write (lu,'("crys")')
    write (lu,'(6(A,X))') (string(c%aa(i)*bohrtoa,'f',18,10),i=1,3), (string(c%bb(i),'f',18,10),i=1,3)
    do i = 1, c%ncel
       write (lu,'(A,X,1p,3(A,X))') string(nameguess(c%spc(c%atcel(i)%is)%z,.true.)), &
          (string(c%atcel(i)%x(j),'f',18,10),j=1,3)
    end do
    write (lu,'("end")')
    call fclose(lu)

  end subroutine write_db

  !> Write a gulp input script
  module subroutine write_gulp(c,file)
    use tools_io, only: fopen_write, nameguess, fclose, string
    use param, only: bohrtoa
    class(crystal), intent(inout) :: c
    character*(*), intent(in) :: file

    integer :: lu, i, j

    lu = fopen_write(file)
    write (lu,'("eem")')
    write (lu,'("cell ",6(A,X))') (string(c%aa(j) * bohrtoa,'f',13,9),j=1,3), &
       (string(c%bb(j),'f',10,5),j=1,3)
    write (lu,'("fractional")')
    do i = 1, c%ncel
       write (lu,'(A5,X,3(A,X))') trim(c%spc(c%atcel(i)%is)%name),&
          (string(c%atcel(i)%x(j),'f',15,9),j=1,3)
    end do

    call fclose(lu)

  end subroutine write_gulp

  !> Write a lammps data file
  module subroutine write_lammps(c,file)
    use tools_io, only: fopen_write, ferror, faterr, fclose
    use tools_math, only: m_x2c_from_cellpar
    use param, only: bohrtoa, atmass
    class(crystal), intent(in) :: c
    character*(*), intent(in) :: file

    integer :: i, j, l
    integer :: lu
    real*8 :: rnew(3,3)

    lu = fopen_write(file)

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
    write (lu,'(2(F18.10,X)," xlo xhi")') 0d0, c%m_x2c(1,1)*bohrtoa
    write (lu,'(2(F18.10,X)," ylo yhi")') 0d0, c%m_x2c(2,2)*bohrtoa
    write (lu,'(2(F18.10,X)," zlo zhi")') 0d0, c%m_x2c(3,3)*bohrtoa
    write (lu,'(3(F18.10,X)," xy xz yz")') 0d0, 0d0, 0d0
    write (lu,*)

    write (lu,'("Masses"/)')
    do i = 1, c%nspc
       write (lu,'(I3,X,F10.4)') i, atmass(c%spc(i)%z)
    end do
    write (lu,*)

    write (lu,'("Atoms"/)')
    l = 0
    do i = 1, c%nspc
       do j = 1, c%ncel
          if (c%atcel(j)%is /= i) cycle
          l = l + 1
          write (lu,'(I7,X,I3,X,F4.1,3(F15.8,X))') l, i, 0d0, c%atcel(j)%r*bohrtoa
       end do
    end do

    call fclose(lu)

  end subroutine write_lammps

  !> Write a siesta fdf data file
  module subroutine write_siesta_fdf(c,file)
    use tools_io, only: fopen_write, nameguess, lower, fclose
    use param, only: bohrtoa
    class(crystal), intent(in) :: c
    character*(*), intent(in) :: file

    integer :: i, j
    integer :: lu

    lu = fopen_write(file)

    ! header
    write (lu,'("# fdf file created by critic2.",/)')
    write (lu,'("SystemName crystal")') 
    write (lu,'("SystemLabel crystal")') 
    write (lu,*)

    write (lu,'("NumberOfSpecies ",I3)') c%nspc
    write (lu,'("NumberOfAtoms ", I6)') c%ncel
    write (lu,'("%block Chemical_Species_Label")') 
    do i = 1, c%nspc
       write (lu,'(I3,I3,X,A2)') i, c%spc(i)%z, lower(nameguess(c%spc(i)%z,.true.))
    end do
    write (lu,'("%endblock Chemical_Species_Label")') 
    write (lu,*)

    write (lu,'("LatticeConstant 1.0 ang")')
    write (lu,'("%block LatticeParameters")')
    write (lu,'(3(F16.10,X),3(F16.8,X))') c%aa*bohrtoa, c%bb
    write (lu,'("%endblock LatticeParameters")')
    write (lu,'("AtomicCoordinatesFormat Fractional")')
    write (lu,'("%block AtomicCoordinatesAndAtomicSpecies")')
    do i = 1, c%nspc
       do j = 1, c%ncel
          if (c%atcel(j)%is == i) &
             write (lu,'(3(F18.12,X),I3)') c%atcel(j)%x, i
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
  module subroutine write_siesta_in(c,file)
    use tools_io, only: fopen_write, uout, nameguess, string, fclose
    use param, only: bohrtoa
    class(crystal), intent(in) :: c
    character*(*), intent(in) :: file

    integer :: lu
    real*8 :: r(3,3)
    integer :: i, j, k

    lu = fopen_write(file)

    ! lattice vectors
    r = transpose(c%m_x2c) * bohrtoa
    do i = 1, 3
       write (lu,'(3(F20.12,X))') r(i,:)
    end do

    ! atoms
    write (lu,*) c%ncel
    j = 0
    do i = 1, c%nspc
       do k = 1, c%ncel
          if (c%atcel(k)%is == i) &
             write (lu,'(I3,X,I3,X,3(F20.12,X))') i, c%spc(i)%z, c%atcel(k)%x
       end do
    end do

    call fclose(lu)

    ! Write the chemical species block to standard output
    write (uout,'("%block Chemical_Species_Label")') 
    do i = 1, c%nspc
       write (uout,'(3(2X,A))') string(i), string(c%spc(i)%z), &
          string(nameguess(c%spc(i)%z,.true.))
    end do
    write (uout,'("%endblock Chemical_Species_Label")') 
    write (uout,*)

  end subroutine write_siesta_in

  !> Write a DFTB+ human-friendly structured data format (hsd) file
  module subroutine write_dftbp_hsd(c,file)
    use tools_io, only: fopen_write, string, nameguess, fclose
    use param, only: maxzat0
    class(crystal), intent(in) :: c
    character*(*), intent(in) :: file

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

    lu = fopen_write(file)
    write (lu,'("Geometry = GenFormat {")')
    call c%write_dftbp_gen(file,lu)
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
       write (lu,'(4X,A," = ",A)') string(nameguess(c%spc(i)%z,.true.)), &
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
       write (lu,'(4X,A," = ",A)') string(nameguess(c%spc(i)%z,.true.)), &
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

  !> Write a DFTB+ human-friendly gen structure file
  module subroutine write_dftbp_gen(c,file,lu0)
    use tools_io, only: fopen_write, nameguess, string, fclose
    use param, only: bohrtoa
    class(crystal), intent(in) :: c
    character*(*), intent(in) :: file
    integer, intent(in), optional :: lu0

    integer :: lu, i, j, k
    real*8 :: r(3,3)
    character(len=:), allocatable :: strtyp, aux

    ! open file
    if (present(lu0)) then
       lu = lu0
    else
       lu = fopen_write(file)
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
                write (lu,'(99(A,X))') string(k), string(i), &
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
                write (lu,'(99(A,X))') string(k), string(i), &
                (string(c%atcel(k)%x(j),'f',20,12),j=1,3)
          end do
       end do

       ! lattice vectors
       r = c%m_x2c * bohrtoa
       write (lu,'(3(A,X))') (string(0d0,'f',20,12),j=1,3)
       do i = 1, 3
          write (lu,'(3(A,X))') (string(r(j,i),'f',20,12),j=1,3)
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
  module subroutine writegrid_cube(c,g,file,onlyheader,binary,xd0,x00)
    use global, only: precisecube
    use tools_io, only: fopen_write, fclose
    use param, only: eye
    class(crystal), intent(in) :: c
    real*8, intent(in) :: g(:,:,:)
    character*(*), intent(in) :: file
    logical, intent(in) :: onlyheader
    logical, intent(in) :: binary
    real*8, intent(in), optional :: xd0(3,3)
    real*8, intent(in), optional :: x00(3)

    integer :: n(3), i, ix, iy, iz, lu
    real*8 :: xd(3,3), x0(3)

    do i = 1, 3
       n(i) = size(g,i)
    end do
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

    if (binary) then
       lu = fopen_write(file,form="unformatted")
       write(lu) c%ncel, x0
       write(lu) n, xd
       do i = 1, c%ncel
          write(lu) c%spc(c%atcel(i)%is)%z, 0d0, c%atcel(i)%r(:) + c%molx0
       end do
       write (lu) g
    else
       lu = fopen_write(file)
       write(lu,'("critic2-cube")')
       write(lu,'("critic2-cube")')
       if (precisecube) then
          write(lu,'(I5,1x,3(E22.14,1X))') c%ncel, x0
          do i = 1, 3
             write(lu,'(I5,1x,3(E22.14,1X))') n(i), xd(:,i)
          end do
          do i = 1, c%ncel
             write(lu,'(I4,F5.1,3(E22.14,1X))') c%spc(c%atcel(i)%is)%z, 0d0, c%atcel(i)%r(:) + c%molx0
          end do
       else
          write(lu,'(I5,3(F12.6))') c%ncel, x0
          do i = 1, 3
             write(lu,'(I5,3(F12.6))') n(i), xd(:,i)
          end do
          do i = 1, c%ncel
             write(lu,'(I4,F5.1,F11.6,F11.6,F11.6)') c%spc(c%atcel(i)%is)%z, 0d0, c%atcel(i)%r(:) + c%molx0
          end do
       end if
       if (.not.onlyheader) then
          do ix = 1, n(1)
             do iy = 1, n(2)
                if (precisecube) then
                   write (lu,'(6(1x,e22.14))') (g(ix,iy,iz),iz=1,n(3))
                else
                   write (lu,'(1p,6(1x,e12.5))') (g(ix,iy,iz),iz=1,n(3))
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
  module subroutine writegrid_vasp(c,g,file,onlyheader)
    use tools_io, only: fopen_write, string, nameguess, fclose
    use param, only: bohrtoa
    class(crystal), intent(in) :: c
    real*8, intent(in) :: g(:,:,:)
    character*(*), intent(in) :: file
    logical :: onlyheader

    integer :: n(3), i, j, ix, iy, iz, lu
    character(len=:), allocatable :: line0, aux

    do i = 1, 3
       n(i) = size(g,i)
    end do

    lu = fopen_write(file)
    write (lu,'("CHGCAR generated by critic2")')
    write (lu,'("1.0")')
    do i = 1, 3
       write (lu,'(1p,3(E22.14,X))') c%m_x2c(:,i) * bohrtoa
    end do

    line0 = ""
    do i = 1, c%nspc
       aux = line0 // " " // string(nameguess(c%spc(i)%z,.true.))
       line0 = aux
    end do
    write (lu,'(A)') line0
    write (lu,'("Direct")')
    do i = 1, c%nspc
       do j = 1, c%ncel
          if (c%atcel(j)%is == i) cycle
          write (lu,'(1p,3(E22.14,X))') c%atcel(j)%x
       end do
    end do
    write (lu,*)
    write (lu,'(3(I5,X))') n
    if (.not.onlyheader) then
       write (lu,'(5(1x,e22.14))') (((g(ix,iy,iz)*c%omega,ix=1,n(1)),iy=1,n(2)),iz=1,n(3))
    end if
    call fclose(lu)

  end subroutine writegrid_vasp


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

    !$omp parallel do private(x,rho,rdum1,rdum2) schedule(dynamic)
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

  !xx! private procedures

  !> Calculate the point group of a lattice. rmat is the matrix of
  !> lattice vectors (G = rmat * rmat'). verbose activates output to
  !> uout. ncen and xcen are the centering vectors for the lattice in
  !> crystallographic coordinates. nn and rot, if present, receive the
  !> symmetry operations for the lattice.
  subroutine lattpg(rmat,ncen,xcen,nn,rot)
    use sympg, only: nopsym, opsym, sym3d
    use tools_math, only: matinv
    use types, only: realloc
    real*8, intent(in) :: rmat(3,3)
    integer, intent(in) :: ncen
    real*8, intent(in) :: xcen(3,ncen)
    integer, intent(out), optional :: nn
    real*8, intent(out), optional :: rot(3,3,48)

    real*8 :: rmati(3,3), aal(3), gmat(3,3)
    integer :: i, na, nb, nc, npos, ia, ib, ic, it, op
    real*8  :: amax, amax2e, x(3), t(3), d2
    real*8, allocatable :: ax(:,:)
    integer, allocatable :: atZmol(:)

    ! the reciprocal-space matrix is the transpose of the inverse
    rmati = matinv(rmat)
    gmat = matmul(rmat,transpose(rmat))
    do i = 1, 3
       aal(i) = sqrt(gmat(i,i))
    end do

    ! every lattice vector must be represented in the molecule
    amax = 2d0*maxval(aal)
    amax2e = amax + 1d-5
    na = int((amax/aal(1)) - 1d-7) + 2
    nb = int((amax/aal(2)) - 1d-7) + 2
    nc = int((amax/aal(3)) - 1d-7) + 2

    ! build the lattice point molecule
    npos = 0
    allocate(ax(3,10))
    do ia = -na, na
       do ib = -nb, nb
          do ic = -nc, nc
             do it = 1, ncen
                x = real((/ia,ib,ic/),8) + xcen(:,it)
                t = matmul(x,rmat)
                d2 = norm2(t)
                if (d2 <= amax2e) then
                   npos = npos + 1
                   if (npos > size(ax,2)) call realloc(ax,3,2*npos)
                   ax(:,npos) = t
                end if
             end do
          end do
       end do
    end do
    call realloc(ax,3,npos)
    allocate(atZmol(npos))
    atzmol = 1

    ! Use symmetry to determine the point group. Default options
    call sym3d(rmat,ax,atZmol,npos,.false.)

    ! fill the reciprocal space matrices and clean up
    if (present(nn) .and. present(rot)) then
       nn = nopsym
       do op = 1, nn
          rot(:,:,op) = matmul(transpose(rmati),matmul(opsym(op,:,:),transpose(rmat)))
       end do
    end if
    deallocate(atZmol,ax)

  end subroutine lattpg

  !> For the symmetry operation rot, compute the order and
  !> characteristic direction of the operation. The output is given in
  !> type (the type of operation: ident, c2, c3,...)  and vec (the
  !> direction of the axis or the normal to the plane, 0 for a point
  !> symmetry element). Used in the structure initialization.
  subroutine typeop(rot,type,vec,order)
    use tools_math, only: eigns
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
       call eigns(mat,eigen,eigeni)

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
