! Copyright (c) 2015 Alberto Otero de la Roza <aoterodelaroza@gmail.com>,
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
module yt
  implicit none

  private

  public :: yt_integrate
  public :: yt_weights

contains

  !> Do the YT integration on crystal c and field f. Return the number
  !> of basins (nbasin), their coordinates (cryst. coordinates,
  !> xcoord), the integer id that gives the basin for each grid point
  !> (idg) and the logical unit of an open scratch file containing the
  !> weights (luw). If the arithmetic expression discexpr is not
  !> empty, then apply that expression to the basin attractors. If the
  !> expression is non-zero, discard the attractor. If atexist is
  !> true, then the code is aware of the presence of atoms, which are
  !> added as attractors at the beginning of the run. Two attractors
  !> are considered equal if they are within a ditsance of ratom
  !> (bohr).
  subroutine yt_integrate(c,ff,discexpr,atexist,ratom,nbasin,xcoord,idg,luw)
    use struct_basic, only: crystal
    use fields, only: fields_fcheck, fields_feval
    use tools_io, only: ferror, faterr, fopen_scratch
    use arithmetic, only: eval
    use param, only: vsmall
    use tools, only: qcksort
    use types, only: realloc
    type(crystal), intent(in) :: c
    real*8, intent(in) :: ff(:,:,:)
    character*(*), intent(in) :: discexpr
    logical, intent(in) :: atexist
    real*8, intent(in) :: ratom
    integer, intent(out) :: nbasin
    real*8, allocatable, intent(inout) :: xcoord(:,:)
    integer, allocatable, intent(inout) :: idg(:,:,:)
    integer, intent(out) :: luw

    real*8, allocatable :: g(:)
    integer, allocatable :: io(:), iio(:)
    integer :: i, ii, j, n(3), nn, nvec, vec(3,16), ib(3), jb(3), jj, k, kk
    real*8 :: al(40), csum
    integer :: nhi
    integer, allocatable :: ibasin(:), ihi(:), inear(:,:), nlo(:)
    real*8, allocatable :: chi(:), fnear(:,:)
    logical :: isias, isassigned, ok
    integer :: nid, lvec(3)
    real*8 :: dist, dv(3), fval, x(3)
    type(crystal) :: caux

    ! Pre-allocate atoms as maxima
    allocate(xcoord(3,c%ncel))
    xcoord = 0d0
    if (atexist) then
       nbasin = c%ncel
       do i = 1, c%ncel
          xcoord(:,i) = c%atcel(i)%x
       end do
    else
       nbasin = 0
    end if

    ! Copy the field onto a one-dimensional array
    do i = 1, 3
       n(i) = size(ff,i)
    end do
    nn = n(1)*n(2)*n(3)
    allocate(g(nn))
    g = reshape(ff,shape(g))

    ! sort g, from smaller to larger field value
    allocate(io(nn),iio(nn))
    do i = 1, nn
       io(i) = i
    end do
    call qcksort(g,io,1,nn)
    do i = 1, nn
       iio(io(i)) = i
    end do

    ! calculate areas*lengths and grid vectors. Use a smaller crystal
    ! where the "lattice" corresponds to the cube grid points.
    caux%isinit = .true.
    caux%aa = c%aa / real(n,8)
    caux%bb = c%bb
    call caux%set_cryscar()
    call caux%wigner((/0d0,0d0,0d0/),nvec=nvec,vec=vec,area0=al)

    ! run over grid points in order of decreasing density
    allocate(ibasin(nn),ihi(nvec),chi(nvec),inear(nvec,nn),fnear(nvec,nn),nlo(nn))
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
             chi(nhi) = al(k) * (g(j) - g(i))
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
          if (atexist) then
             nid = 0
             call c%nearest_atom(dv,nid,dist,lvec)
             if (dist < ratom) then
                ibasin(ii) = nid
                isassigned = .true.
             end if
          end if
          ! check if it is a known nnm
          if (.not.isassigned .and. ratom > vsmall) then
             do k = 1, nbasin
                if (c%are_lclose(dv,xcoord(:,k),ratom)) then
                   ibasin(ii) = k
                   isassigned = .true.
                   exit
                end if
             end do
          end if
          ! well, it must be a new attractor then
          if (.not.isassigned) then
             ok = .true.
             if (len_trim(discexpr) > 0) then
                x = c%x2c(dv)
                fval = eval(discexpr,.false.,ok,x,fields_fcheck,fields_feval)
                if (.not.ok) &
                   call ferror("yt","invalid DISCARD expression",faterr)
                ok = (abs(fval) < 1d-30)
             end if
             if (ok) then
                nbasin = nbasin + 1
                if (nbasin > size(xcoord,2)) call realloc(xcoord,3,2*nbasin)
                ibasin(ii) = nbasin
                xcoord(:,nbasin) = dv
             end if
          end if
       else
          isias = ibasin(ihi(1)) == 0
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
                fnear(nlo(kk),kk) = chi(k) / max(csum,1d-40)
             end do
          end if
       end if
    end do

    ! clean up a little, we calculate the weights now
    deallocate(io,ihi,chi)

    ! Save the necessary information to reproduce the weights to an
    ! external file. Saving the weights directly gives files
    ! that are too large
    luw = fopen_scratch()
    write (luw) nbasin, nn, nvec
    write (luw) nlo
    write (luw) ibasin
    write (luw) iio
    write (luw) inear
    write (luw) fnear
    rewind(luw)
    deallocate(inear,fnear,nlo,g,iio)

    ! clean up and output 
    allocate(idg(n(1),n(2),n(3)))
    idg = reshape(ibasin,shape(idg))
    deallocate(ibasin)
    call realloc(xcoord,3,nbasin)

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

  !> Read the neighbor and fractions from the external file and
  !> generate the YT weights for the given input basin.
  subroutine yt_weights(luw,idb,w)
    use tools_io, only: ferror, faterr

    integer, intent(in) :: luw !< Logical unit for the YT checkpoint
    integer, intent(in) :: idb !< Basin to integrate
    real*8, intent(inout) :: w(:,:,:) !< Output weights for this basin

    integer :: nbasin, nn, nvec
    integer, allocatable :: nlo(:), ibasin(:), iio(:), inear(:,:)
    real*8, allocatable :: fnear(:,:), waux(:)
    logical :: opened, exist
    integer :: i, j, k

    ! read the yt checkpoint file header
    inquire(luw,opened=opened,exist=exist)
    if (.not.opened.or..not.exist) &
       call ferror('yt_weights','YT checkpoint is not open',faterr)
    rewind(luw)
    read (luw) nbasin, nn, nvec
    if (idb < 0 .or. idb > nn) &
       call ferror('yt_weights','unknown basin',faterr)

    ! allocate and rest of the yt checkpoint
    allocate(ibasin(nn),nlo(nn),inear(nvec,nn),fnear(nvec,nn),iio(nn))
    read (luw) nlo
    read (luw) ibasin
    read (luw) iio
    read (luw) inear
    read (luw) fnear

    ! prepare the output grid
    allocate(waux(nn))

    ! which nodes belong to this basin?
    where (ibasin == idb)
       waux = 1d0
    elsewhere
       waux = 0d0
    end where

    ! recompute w for the ias points
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

    ! clean up
    deallocate(nlo,ibasin,iio,inear,fnear,waux)
  
  end subroutine yt_weights

end module yt
