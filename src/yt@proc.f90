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

! This module contains code for the Yu and Trinkle Bader integration-on-a-grid.
! See:
!    Min Yu, Dallas Trinkle, J. Chem. Phys. 134 (2011) 064111.
! for details on the algorithm.
submodule (yt) proc
  implicit none

contains

  !> Do the YT integration on system s and field iref. Return the
  !> number of basins (bas%nattr), their coordinates (crystallographic
  !> coordinates, bas%xattr), the integer id that gives the basin for
  !> each grid point (bas%idg) and the logical unit of an open scratch
  !> file containing the weights (bas%luw).  If the arithmetic
  !> expression bas%expr is not empty, then apply that expression to
  !> the basin attractors. If the expression is non-zero, discard the
  !> attractor. If bas%atexist is true, then the code is aware of the
  !> presence of atoms, which are added as attractors at the beginning
  !> of the run. Two attractors are considered equal if they are
  !> within a ditsance of bas%ratom (bohr).
  module subroutine yt_integrate(s,bas)
    use systemmod, only: system
    use tools_math, only: m_x2c_from_cellpar, matinv
    use tools_io, only: ferror, faterr, fopen_scratch
    use param, only: vsmall, icrd_crys
    use tools, only: qcksort
    use types, only: realloc, basindat
    type(system), intent(inout) :: s
    type(basindat), intent(inout) :: bas

    real*8, allocatable :: g(:)
    integer, allocatable :: io(:), iio(:), vec(:,:)
    integer :: i, ii, j, n(3), nn, nvec, ib(3), jb(3), jj, k, kk
    real*8 :: csum
    integer :: nhi
    integer, allocatable :: ibasin(:), ihi(:), inear(:,:), nlo(:)
    real*8, allocatable :: chi(:), fnear(:,:), al(:)
    logical :: isias, isassigned, ok
    integer :: nid
    real*8 :: dv(3), fval, x(3)
    character(len=:), allocatable :: errmsg

    if (.not.s%isinit) &
       call ferror("yt_integrate","system not initialized",faterr)
    if (.not.allocated(s%c)) &
       call ferror("yt_integrate","system does not have crystal",faterr)

    ! Pre-allocate atoms and nnm as maxima
    allocate(bas%xattr(3,s%f(s%iref)%ncpcel))
    bas%xattr = 0d0
    if (bas%atexist) then
       bas%nattr = s%f(s%iref)%ncpcel
       do i = 1, s%f(s%iref)%ncpcel
          bas%xattr(:,i) = s%f(s%iref)%cpcel(i)%x
       end do
    else
       bas%nattr = 0
    end if

    ! Copy the field onto a one-dimensional array
    do i = 1, 3
       n(i) = size(bas%f,i)
    end do
    nn = n(1)*n(2)*n(3)
    allocate(g(nn))
    g = reshape(bas%f,shape(g))

    ! sort g, from smaller to larger field value
    allocate(io(nn),iio(nn))
    do i = 1, nn
       io(i) = i
    end do
    call qcksort(g,io,1,nn)
    do i = 1, nn
       iio(io(i)) = i
    end do

    ! copy the voronoi-relevant vectors
    nvec = s%f(s%iref)%grid%nvec
    vec = s%f(s%iref)%grid%vec(:,1:nvec)
    al = s%f(s%iref)%grid%area(1:nvec)

    ! run over grid points in order of decreasing density
    allocate(ibasin(nn),ihi(nvec),chi(nvec),inear(nvec,nn),fnear(nvec,nn),nlo(nn))
    ibasin = 0
    ihi = 0
    chi = 0d0
    nlo = 0
    inear = 0
    fnear = 0d0
    do ii = nn, 1, -1
       ! find the number of points with higher density
       nhi = 0
       i = io(ii)
       ib = to3(i)
       csum = 0d0
       chi = 0d0
       do k = 1, nvec
          jb = ib + vec(:,k)
          j = to1(jb)
          jj = iio(j)
          if (jj > ii) then
             nhi = nhi + 1
             ihi(nhi) = jj
             chi(nhi) = max(al(k) * (g(j) - g(i)),vsmall)
             csum = csum + chi(nhi)
          end if
       end do

       ! classify the point
       nlo(ii) = 0
       if (nhi == 0) then
          ! this is a local maximum
          dv = real(ib-1,8) / n

          ! check if it is an atom (use ratom)
          isassigned = .false.
          if (bas%atexist) then
             nid = s%c%identify_atom(dv,icrd_crys,distmax=bas%ratom)
             if (nid > 0) then
                ibasin(ii) = nid
                isassigned = .true.
             end if
          end if
          ! check if it is a known nnm
          if (.not.isassigned .and. bas%ratom > vsmall) then
             do k = 1, bas%nattr
                if (s%c%are_lclose(dv,bas%xattr(:,k),bas%ratom)) then
                   ibasin(ii) = k
                   isassigned = .true.
                   exit
                end if
             end do
          end if
          ! well, it must be a new attractor then
          if (.not.isassigned) then
             ok = .true.
             if (len_trim(bas%expr) > 0) then
                x = s%c%x2c(dv)
                fval = s%eval(bas%expr,errmsg,x)
                if (len_trim(errmsg) > 0) &
                   call ferror("yt","Invalid DISCARD expression: "//trim(errmsg),faterr)
                ok = (abs(fval) < vsmall)
             end if
             if (ok) then
                bas%nattr = bas%nattr + 1
                if (bas%nattr > size(bas%xattr,2)) call realloc(bas%xattr,3,2*bas%nattr)
                ibasin(ii) = bas%nattr
                bas%xattr(:,bas%nattr) = dv
             end if
          end if
       else
          isias = (ibasin(ihi(1)) == 0)
          do k = 1, nhi
             isias = isias .or. (ibasin(ihi(k)) /= ibasin(ihi(1)))
          end do
          if (.not.isias) then
             ! interior point
             ibasin(ii) = ibasin(ihi(1))
          else
             ! ias point, save the neighbor info for later
             ibasin(ii) = 0
             do k = 1, nhi
                kk = ihi(k)
                nlo(kk) = nlo(kk) + 1
                inear(nlo(kk),kk) = ii
                fnear(nlo(kk),kk) = chi(k) / max(csum,vsmall)
             end do
          end if
       end if
    end do

    ! clean up a little, we calculate the weights now
    deallocate(io,ihi,chi)

    ! Save the necessary information to reproduce the weights to an
    ! external file. Saving the weights directly gives files
    ! that are too large.
    bas%luw = fopen_scratch()
    write (bas%luw) bas%nattr, nn, nvec
    write (bas%luw) nlo
    write (bas%luw) ibasin
    write (bas%luw) iio
    write (bas%luw) inear
    write (bas%luw) fnear
    call flush(bas%luw)
    rewind(bas%luw)
    deallocate(inear,fnear,nlo,g,iio)

    ! clean up and output
    allocate(bas%idg(n(1),n(2),n(3)))
    bas%idg = reshape(ibasin,shape(bas%idg))
    deallocate(ibasin)
    call realloc(bas%xattr,3,bas%nattr)

  contains
    function to3(k)
      integer :: to3(3), k
      to3(1) = modulo(k-1,n(1))+1
      to3(2) = modulo((k-1) / n(1),n(2)) + 1
      to3(3) = modulo((k-1) / (n(1)*n(2)),n(3)) + 1
    end function to3
    function to1(k)
      integer :: k(3), to1
      to1 = modulo(k(1)-1,n(1)) + n(1) * (modulo(k(2)-1,n(2)) + n(2) * (modulo(k(3)-1,n(3)))) + 1
    end function to1
  end subroutine yt_integrate

  !> Calculate the isosurface regions associated with a value higher
  !> than the contour value bas%isov. On output, return the number of
  !> regions (bas%nattr), the coordinate of the highest-function
  !> maximum inside (crystallographic coordinates, bas%xattr), and the
  !> integer id that gives the region for each grid point (bas%idg) If
  !> the arithmetic expression bas%expr is not empty, then discard the
  !> grid points where the expression is true.
  module subroutine yt_isosurface(s,bas)
    use types, only: basindat
    use crystalmod, only: crystal
    use tools_math, only: m_x2c_from_cellpar, matinv
    use tools_io, only: ferror, faterr
    use tools, only: qcksort
    use types, only: realloc
    use param, only: vsmall
    type(system), intent(inout) :: s
    type(basindat), intent(inout) :: bas

    real*8, allocatable :: g(:)
    integer :: nn, n(3)
    integer, allocatable :: io(:), iio(:), vec(:,:)
    integer :: i, ii, j, jj, k, ib(3), jb(3), nvec
    integer :: nhi, imin
    integer, allocatable :: ibasin(:), ihi(:)
    real*8 :: fval, x(3), dv(3)
    logical :: interior
    integer, allocatable :: imap(:)
    character(len=:), allocatable :: errmsg

    if (.not.s%isinit) &
       call ferror("yt_isosurface","system not initialized",faterr)
    if (.not.allocated(s%c)) &
       call ferror("yt_isosurface","system does not have crystal",faterr)

    ! initialize
    allocate(bas%xattr(3,10),imap(10))
    bas%xattr = 0d0
    bas%nattr = 0
    imap = 0

    ! Copy the field onto a one-dimensional array
    do i = 1, 3
       n(i) = size(bas%f,i)
    end do
    nn = n(1)*n(2)*n(3)
    allocate(g(nn))
    g = reshape(bas%f,shape(g))

    ! sort g, from smaller to larger field value
    allocate(io(nn),iio(nn))
    do i = 1, nn
       io(i) = i
    end do
    call qcksort(g,io,1,nn)
    do i = 1, nn
       iio(io(i)) = i
    end do

    ! copy the voronoi-relevant vectors
    nvec = s%f(s%iref)%grid%nvec
    vec = s%f(s%iref)%grid%vec(:,1:nvec)

    ! run over grid points in order of decreasing density
    allocate(ibasin(nn),ihi(nvec))
    ibasin = 0
    ihi = 0
    do ii = nn, 1, -1
       i = io(ii)
       ib = to3(i)

       ! skip the points lower than the contour value
       if (g(i) < bas%isov) cycle

       ! if the discard expression is present and true, cycle
       dv = real(ib-1,8) / n
       if (len_trim(bas%expr) > 0) then
          x = s%c%x2c(dv)
          fval = s%eval(bas%expr,errmsg,x)
          if (len_trim(errmsg) > 0) &
             call ferror("yt","Invalid DISCARD expression: " // trim(errmsg),faterr)
          if (abs(fval) > vsmall) cycle
       end if

       ! find the number of points with higher density
       nhi = 0
       do k = 1, nvec
          jb = ib + vec(:,k)
          j = to1(jb)
          jj = iio(j)
          if (jj > ii) then
             nhi = nhi + 1
             ihi(nhi) = jj
          end if
       end do

       ! classify the point
       if (nhi == 0) then
          ! this is a maximum
          bas%nattr = bas%nattr + 1
          if (bas%nattr > size(bas%xattr,2)) then
             call realloc(bas%xattr,3,2*bas%nattr)
             call realloc(imap,2*bas%nattr)
             imap(bas%nattr+1:) = 0
          end if
          ibasin(ii) = bas%nattr
          bas%xattr(:,bas%nattr) = dv
       else
          interior = .true.
          do k = 2, nhi
             interior = interior .and. (ibasin(ihi(k)) == ibasin(ihi(1)))
          end do
          if (interior) then
             ! interior point
             ibasin(ii) = ibasin(ihi(1))
          else
             ! contact between isosurfaces, reassign indices
             imin = minval(ibasin(ihi(1:nhi)))
             ibasin(ii) = imin
             do k = 1, nhi
                if (ibasin(ihi(k)) /= imin) imap(ibasin(ihi(k))) = imin
             end do
          end if
       end if
    end do

    ! clean up
    deallocate(io,ihi,g)

    ! reassing maxima
    nn = 0
    do i = 1, bas%nattr
       if (imap(i) == 0) then
          nn = nn + 1
       else
          ii = i
          do while (imap(ii) /= 0)
             ii = imap(ii)
          end do
          where (ibasin == i)
             ibasin = ii
          end where
       end if
    end do
    bas%nattr = nn
    deallocate(imap)

    ! clean up and output
    allocate(bas%idg(n(1),n(2),n(3)))
    bas%idg = reshape(ibasin(iio),shape(bas%idg))
    deallocate(ibasin,iio)
    call realloc(bas%xattr,3,bas%nattr)

  contains
    function to3(k)
      integer :: to3(3), k
      to3(1) = modulo(k-1,n(1))+1
      to3(2) = modulo((k-1) / n(1),n(2)) + 1
      to3(3) = modulo((k-1) / (n(1)*n(2)),n(3)) + 1
    end function to3
    function to1(k)
      integer :: k(3), to1
      to1 = modulo(k(1)-1,n(1)) + n(1) * (modulo(k(2)-1,n(2)) + n(2) * (modulo(k(3)-1,n(3)))) + 1
    end function to1

  end subroutine yt_isosurface

  !> Read or use the neighbor and fractions and generate the YT
  !> weights for the given input basin. The input for the weight has
  !> to come from a file (logical unit luw), or from the data passed
  !> in the din argument. On output, this routine gives the read
  !> data (dout) and the weights for basin idb. All arguments are
  !> optional. One of luw and din has to be given (but not both).
  !> One of dout and (w,idb) is allowed (but not both).
  module subroutine yt_weights(luw,din,idb,w,dout)
    use tools_io, only: ferror, faterr

    integer, intent(in), optional :: luw !< Logical unit for the YT checkpoint
    type(ytdata), intent(in), optional :: din !< Data for generating weights (input)
    integer, intent(in), optional :: idb !< Basin to calculate the weights for
    type(ytdata), intent(inout), optional :: dout !< Data for generating weights (output)
    real*8, intent(inout), optional :: w(:,:,:) !< Output weights for this basin

    integer :: nbasin, nn, nvec
    integer, allocatable :: nlo(:), ibasin(:), iio(:), inear(:,:)
    real*8, allocatable :: fnear(:,:), waux(:)
    logical :: opened, exist
    integer :: i, j, k

    if (present(luw) .and. present(din)) &
       call ferror('yt_weights','luw and din both present in yt_weights',faterr)
    if (.not.present(luw) .and. .not.present(din)) &
       call ferror('yt_weights','luw and din both absent in yt_weights',faterr)
    if (present(dout) .and. present(w)) &
       call ferror('yt_weights','dout and w/idb both present in yt_weights',faterr)
    if (.not.present(dout) .and. .not.present(w)) &
       call ferror('yt_weights','dout and w/idb both absent in yt_weights',faterr)

    ! read the yt checkpoint file header
    if (present(luw)) then
       inquire(luw,opened=opened,exist=exist)
       if (.not.opened.or..not.exist) &
          call ferror('yt_weights','YT checkpoint is not open',faterr)
       rewind(luw)
       read (luw) nbasin, nn, nvec

       ! allocate and rest of the yt checkpoint
       allocate(ibasin(nn),nlo(nn),inear(nvec,nn),fnear(nvec,nn),iio(nn))
       read (luw) nlo
       read (luw) ibasin
       read (luw) iio
       read (luw) inear
       read (luw) fnear
    else
       nbasin = din%nbasin
       nn = din%nn
       nvec = din%nvec
    end if

    ! write dout
    if (present(dout)) then
       if (present(luw)) then
          dout%nbasin = nbasin
          dout%nn = nn
          dout%nvec = nvec
          if (allocated(dout%nlo)) deallocate(dout%nlo)
          call move_alloc(nlo,dout%nlo)
          if (allocated(dout%ibasin)) deallocate(dout%ibasin)
          call move_alloc(ibasin,dout%ibasin)
          if (allocated(dout%iio)) deallocate(dout%iio)
          call move_alloc(iio,dout%iio)
          if (allocated(dout%inear)) deallocate(dout%inear)
          call move_alloc(inear,dout%inear)
          if (allocated(dout%fnear)) deallocate(dout%fnear)
          call move_alloc(fnear,dout%fnear)
       else
          dout = din
       end if
       return
    end if

    ! Calculate the weight for this attractor
    if (.not.present(w).or..not.present(idb)) &
       call ferror('yt_weights','one of w/idb is missing',faterr)
    if (idb < 0 .or. idb > nn) &
       call ferror('yt_weights','unknown basin',faterr)

    ! prepare the output grid
    allocate(waux(nn))

    ! recompute w for the ias points
    if (present(luw)) then
       where (ibasin == idb)
          waux = 1d0
       elsewhere
          waux = 0d0
       end where

       do j = nn, 1, -1
          if (abs(waux(j)) > 0d0) then
             do k = 1, nlo(j)
                waux(inear(k,j)) = waux(inear(k,j)) + fnear(k,j) * waux(j)
             end do
          end if
       end do

       nn = 0
       do k = lbound(w,3), ubound(w,3)
          do j = lbound(w,2), ubound(w,2)
             do i = lbound(w,1), ubound(w,1)
                nn = nn + 1
                w(i,j,k) = waux(iio(nn))
             end do
          end do
       end do
       deallocate(nlo,ibasin,iio,inear,fnear)
    else
       where (din%ibasin == idb)
          waux = 1d0
       elsewhere
          waux = 0d0
       end where

       do j = nn, 1, -1
          if (abs(waux(j)) > 0d0) then
             do k = 1, din%nlo(j)
                waux(din%inear(k,j)) = waux(din%inear(k,j)) + din%fnear(k,j) * waux(j)
             end do
          end if
       end do

       nn = 0
       do k = lbound(w,3), ubound(w,3)
          do j = lbound(w,2), ubound(w,2)
             do i = lbound(w,1), ubound(w,1)
                nn = nn + 1
                w(i,j,k) = waux(din%iio(nn))
             end do
          end do
       end do
    end if

    ! clean up
    deallocate(waux)

  end subroutine yt_weights

  !> Remap the attractors from a yt calculation for DIs.
  module subroutine yt_remap(s,bas,dat,nattn,ilvec,iatt)
    use types, only: basindat, realloc
    type(system), intent(in) :: s
    type(basindat), intent(in) :: bas
    type(ytdata), intent(in) :: dat
    integer, intent(out) :: nattn
    integer, allocatable, intent(inout) :: iatt(:), ilvec(:,:)

    integer :: i, j, m1, m2, m3, p(3)
    real*8 :: x(3), xs(3), d2
    logical :: found
    real*8, allocatable :: w(:,:,:)

    if (allocated(iatt)) deallocate(iatt)
    allocate(iatt(bas%nattr))
    nattn = bas%nattr
    do i = 1, bas%nattr
       iatt(i) = i
    enddo
    if (allocated(ilvec)) deallocate(ilvec)
    allocate(ilvec(3,bas%nattr))
    ilvec = 0

    allocate(w(bas%n(1),bas%n(2),bas%n(3)))
    do i = 1, bas%nattr
       call yt_weights(din=dat,idb=i,w=w)
       do m3 = 1, bas%n(3)
          do m2 = 1, bas%n(2)
             do m1 = 1, bas%n(1)
                if (abs(w(m1,m2,m3)) < 1d-15) cycle
                p = (/m1,m2,m3/)
                x = real(p-1,8) / bas%n - bas%xattr(:,i)
                xs = x
                call s%c%shortest(xs,d2)
                p = nint(x - s%c%c2x(xs))
                if (any(p /= 0)) then
                   found = .false.
                   do j = bas%nattr+1, nattn
                      if (iatt(j) == i .and. all(p == ilvec(:,j))) then
                         found = .true.
                         exit
                      end if
                   end do
                   if (.not.found) then
                      nattn = nattn + 1
                      if (nattn > size(ilvec,2)) then
                         call realloc(ilvec,3,2*nattn)
                         call realloc(iatt,2*nattn)
                      end if
                      ilvec(:,nattn) = p
                      iatt(nattn) = i
                   end if
                end if
             end do
          end do
       end do
    end do
    deallocate(w)
    call realloc(ilvec,3,nattn)
    call realloc(iatt,nattn)

  endsubroutine yt_remap

  !> Deallocate a ytdata object
  module subroutine ytdata_clean(dat)
    type(ytdata), intent(inout) :: dat

    if (allocated(dat%nlo)) deallocate(dat%nlo)
    if (allocated(dat%ibasin)) deallocate(dat%ibasin)
    if (allocated(dat%iio)) deallocate(dat%iio)
    if (allocated(dat%inear)) deallocate(dat%inear)
    if (allocated(dat%fnear)) deallocate(dat%fnear)

  end subroutine ytdata_clean

end submodule proc
