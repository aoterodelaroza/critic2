! Copyright (c) 2015 Alberto Otero de la Roza <alberto@fluor.quimica.uniovi.es>,
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

contains

  !> Do the YT integration. Return the number of basins, their coordinates
  !> (integer coordinates), the integer id of each grid point and the logical
  !> unit of an open scratch file containing the weights.
  subroutine yt_integrate(nbasin,xcoord,idatt,idg,ratom,luw)
    use struct
    use tools_io
    use tools
    use fields
    use global
    use struct_basic
    use types
    
    integer, intent(out) :: nbasin
    real*8, allocatable, intent(inout) :: xcoord(:,:)
    integer, allocatable, intent(inout) :: idatt(:)
    integer, allocatable, intent(inout) :: idg(:,:,:)
    real*8, intent(in) :: ratom
    integer, intent(out) :: luw

    real(kind=gk), allocatable :: g(:)
    integer, allocatable :: io(:), iio(:)
    integer :: i, ii, j, n(3), nn, nvec, vec(3,40), ib(3), jb(3), jj, k, kk
    real*8 :: al(40), csum
    integer :: nhi, nbas(cr%ncel)
    integer, allocatable :: ibasin(:), ihi(:), inear(:,:), nlo(:), neqb(:)
    real(kind=gk), allocatable :: chi(:), fnear(:,:), w(:)
    logical :: isias
    integer :: nid, lvec(3), bat(cr%ncel), nnnm
    real*8 :: dist, dv(3)

    ! Copy the field onto a grid
    n = f(refden)%n(:)
    nn = n(1)*n(2)*n(3)
    allocate(g(nn))
    g = max(reshape(f(refden)%f,shape(g)),0d0)

    ! sort g, from smaller to larger
    allocate(io(nn),iio(nn))
    do i = 1, nn
       io(i) = i
    end do
    call qcksort(g,io,1,nn)
    do i = 1, nn
       iio(io(i)) = i
    end do

    ! calculate areas*lengths and grid vectors
    call cr%wigner((/0d0,0d0,0d0/),.false.,n,nvec,vec,al)

    ! run over grid points in order of decreasing density
    allocate(neqb(cr%ncel),idatt(cr%ncel))
    allocate(ibasin(nn),ihi(nvec),chi(nvec),inear(nvec,nn),fnear(nvec,nn),nlo(nn))
    nbasin = 0
    nbas = 0
    neqb = 0
    idatt = 0
    bat = 0
    nnnm = 0
    do ii = nn, 1, -1
       ! find the number of points with higher density
       nhi = 0
       i = io(ii)
       ib = to3(i)
       csum = 0d0
       do k = 1, nvec
          jb = ib + vec(:,k)
          j = to1(jb)
          jj = iio(j)
          if (jj > ii) then
             nhi = nhi + 1
             ihi(nhi) = jj
             chi(nhi) = max(al(k) * (g(j) - g(i)),1d-30)
             csum = csum + chi(nhi)
          end if
       end do

       ! classify the point
       nlo(ii) = 0
       if (nhi == 0) then
          ! this is a local maximum
          dv = real(ib-1,8) / n
          nid = 0
          call cr%nearest_atom(dv,nid,dist,lvec)
          if (dist <= ratom) then
             if (bat(nid) > 0) then
                ibasin(ii) = bat(nid)
             else
                nbasin = nbasin + 1
                if (nbasin > size(neqb)) call realloc(neqb,2*size(neqb))
                if (nbasin > size(idatt)) call realloc(idatt,2*size(neqb))
                bat(nid) = nbasin
                ibasin(ii) = nbasin
                neqb(nbasin) = i
                idatt(nbasin) = nid
             endif
          else
             nbasin = nbasin + 1
             nnnm = nnnm + 1
             if (nbasin > size(neqb)) call realloc(neqb,2*size(neqb))
             if (nbasin > size(idatt)) call realloc(idatt,2*size(neqb))
             ibasin(ii) = nbasin
             neqb(nbasin) = i
             idatt(nbasin) = -nnnm
          endif
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
                fnear(nlo(kk),kk) = chi(k) / csum
             end do
          end if
       end if
    end do

    ! clean up a little, we calculate the weights now
    deallocate(io,ihi,chi)
    allocate(w(nn))

    ! compute weights and save them to an external file, in order
    ! reorder the weights in the output
    luw = fopen_scratch()
    do i = 1, nbasin
       ! which nodes belong to this basin?
       where (ibasin == i)
          w = 1d0
       elsewhere
          w = 0d0
       end where

       ! recompute w for the ias points
       do j = nn, 1, -1
          if (abs(w(j)) > 0d0) then
             do k = 1, nlo(j)
                w(inear(k,j)) = w(inear(k,j)) + fnear(k,j) * w(j)
             end do
          end if
       end do
       
       write (luw) (w(iio(j)),j=1,nn)
    end do
    deallocate(w,inear,fnear,nlo,g,iio)
    rewind(luw)

    ! output 
    allocate(xcoord(3,nbasin))
    do i = 1, nbasin
       if (idatt(i) > 0) then
          xcoord(:,i) = cr%atcel(idatt(i))%x
       else
          xcoord(:,i) = real(to3(neqb(i))-1,8) / n
       end if
    end do
    deallocate(neqb)
    allocate(idg(n(1),n(2),n(3)))
    idg = reshape(ibasin,shape(idg))
    deallocate(ibasin)
    call realloc(idatt,nbasin)

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

end module yt
