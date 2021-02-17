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

submodule (tools) proc
  implicit none

contains

  !xx! uniqing
  !> Returns all diferent vectors v(i,3) (i = 1, n) existing in the
  !> initial set v(k,3) (k=1, n, cartesian coordinates). Only
  !> v(nini,3) through v(nfin,3) are sorted. In the output, nfin
  !> returns the index for the last element. The comparison is
  !> performed in cartesian coordinates.
  module subroutine uniqc(v, nini, nfin, eps)
    integer, intent(in) :: nini !< Starting index
    integer, intent(inout) :: nfin !< End index
    real*8, intent(inout), dimension(3,*) :: v !< The initial and final array
    real*8, intent(in) :: eps

    real*8 :: eps2, x(3), dist2
    integer :: i, j, nfin0
    logical :: found(nini:nfin)

    if (nfin < nini) return

    eps2 = eps * eps

    found = .false.
    !$omp parallel do private(x,dist2)
    do i = nini, nfin
       !$omp flush(found)
       if (found(i)) cycle
       do j = i+1, nfin
          !$omp flush(found)
          if (found(j)) cycle
          x = abs(v(:,i) - v(:,j))
          if (any(x > eps)) cycle
          dist2 = x(1)*x(1) + x(2)*x(2) + x(3)*x(3)
          if (dist2 < eps2) then
             !$omp critical (write)
             found(j) = .true.
             !$omp end critical (write)
             exit
          end if
       end do
    end do
    !$omp end parallel do

    ! rearrange
    nfin0 = nini-1
    do i = nini, nfin
       if (.not.found(i)) then
          nfin0 = nfin0 + 1
          v(:,nfin0) = v(:,i)
       end if
    end do
    nfin = nfin0

  end subroutine uniqc

  !xx! sorting
  !> Sort the elements of the array arr in ascending order using the
  !> quicksort algorithm (real*8 version). iord is the initial order
  !> of data in arr and first and last the intervals for elements to
  !> be analyzed. In the output, iord contains the final order in the
  !> array arr. This sort is not stable.
  module subroutine qcksort_r8(arr, iord, first, last)
    use tools_io, only: ferror, faterr
    !.....Maximum number of elements to be sorted depends on nstack value:
    !     nstack...... ~2log2(last-first+1)
    real*8, dimension(:), intent(in) :: arr !< Array to be sorted
    integer, dimension(:), intent(inout) :: iord !< Permutation array
    integer, intent(in) :: first !< First index
    integer, intent(in) :: last  !< Last index

    integer           m, nstack
    parameter         (m=7, nstack=50)
    real*8            fm, fa, fc, fmi
    parameter         (fm=7875d0, fa=211d0, fc=1663d0, fmi=1.2698413d-4)
    integer           istack(nstack), jstack, l, ir, j, na, i, iq
    real*8            fx, a

    jstack=0
    l=first
    ir=last
    fx=0d0
10  if(ir-l.lt.m)then
       do j=l+1,ir
          na=iord(j)
          a=arr(na)
          do i=j-1,first,-1
             if(arr(iord(i)).le.a) go to 12
             iord(i+1)=iord(i)
          enddo
          i=first-1
12        iord(i+1)=na
       enddo
       if(jstack.eq.0) then
          return
       end if
       ir=istack(jstack)
       l=istack(jstack-1)
       jstack=jstack-2
    else
       i=l
       j=ir
       fx=mod(fx*fa+fc,fm)
       iq=l+int((ir-l+1)*(fx*fmi))
       na=iord(iq)
       a=arr(na)
       iord(iq)=iord(l)
20     continue
21     if (j.ge.first) then
          if (a.lt.arr(iord(j))) then
             j=j-1
             goto 21
          endif
       endif
       if(j.le.i)then
          iord(i)=na
          go to 30
       endif
       iord(i)=iord(j)
       i=i+1
22     if (i.le.last) then
          if (a.gt.arr(iord(i))) then
             i=i+1
             goto 22
          endif
       endif
       if(j.le.i)then
          iord(j)=na
          i=j
          go to 30
       endif
       iord(j)=iord(i)
       j=j-1
       go to 20
30     jstack=jstack+2
       if(jstack.gt.nstack) call ferror ('qcksort','Increase nstack',faterr)
       if(ir-i.ge.i-l)then
          istack(jstack)=ir
          istack(jstack-1)=i+1
          ir=i-1
       else
          istack(jstack)=i-1
          istack(jstack-1)=l
          l=i+1
       endif
    endif
    go to 10
  end subroutine qcksort_r8

  !> Sort the elements of the array arr in ascending order using the
  !> quicksort algorithm (integer*4 version). iord is the initial
  !> order of data in arr and first and last the intervals for
  !> elements to be analyzed. In the output, iord contains the final
  !> order in the array arr. This sort is not stable.
  module subroutine qcksort_i4(iarr, iord, first, last)
    use tools_io, only: ferror, faterr

    integer, dimension(:), intent(in) :: iarr !< Array to be sorted
    integer, dimension(:), intent(inout) :: iord !< Permutation array
    integer, intent(in) :: first !< First index
    integer, intent(in) :: last  !< Last index

    integer, parameter :: m=7, nstack=100
    real*8, parameter :: fm=7875.d0, fa=211.d0, fc=1663.d0, fmi=1.2698413d-4
    integer :: istack(nstack)
    integer :: jstack, l, ir, j, na, a, i, iq
    real*8 :: fx

    jstack=0
    l=first
    ir=last
    fx=0
10  if(ir-l.lt.m)then
       do j=l+1,ir
          na=iord(j)
          a=iarr(na)
          do i=j-1,first,-1
             if(iarr(iord(i)).le.a) go to 12
             iord(i+1)=iord(i)
          enddo
          i=first-1
12        iord(i+1)=na
       enddo
       if(jstack.eq.0)return
       ir=istack(jstack)
       l=istack(jstack-1)
       jstack=jstack-2
    else
       i=l
       j=ir
       fx=mod(fx*fa+fc,fm)
       iq=l+int((ir-l+1)*(fx*fmi))
       na=iord(iq)
       a=iarr(na)
       iord(iq)=iord(l)
20     continue
21     if (j.ge.first) then
          if (a.lt.iarr(iord(j))) then
             j=j-1
             goto 21
          endif
       endif
       if(j.le.i)then
          iord(i)=na
          go to 30
       endif
       iord(i)=iord(j)
       i=i+1
22     if (i.le.last) then
          if (a.gt.iarr(iord(i))) then
             i=i+1
             goto 22
          endif
       endif
       if(j.le.i)then
          iord(j)=na
          i=j
          go to 30
       endif
       iord(j)=iord(i)
       j=j-1
       go to 20
30     jstack=jstack+2
       if(jstack.gt.nstack) call ferror ('qcksort','Increase nstack',faterr)
       if(ir-i.ge.i-l)then
          istack(jstack)=ir
          istack(jstack-1)=i+1
          ir=i-1
       else
          istack(jstack)=i-1
          istack(jstack-1)=l
          l=i+1
       endif
    endif
    go to 10
  end subroutine qcksort_i4

  !> Apply quicksort in-place, real*8 version
  module subroutine qcksort_r8_inplace(arr)
    real*8, intent(inout) :: arr(:)

    integer :: i, n
    integer, allocatable :: iord(:)

    n = size(arr,1)
    allocate(iord(n))
    do i = 1, n
       iord(i) = i
    end do
    call qcksort_r8(arr,iord,1,n)
    arr = arr(iord)
    deallocate(iord)

  end subroutine qcksort_r8_inplace

  !> Apply quicksort in-place, integer*4 version
  module subroutine qcksort_i4_inplace(arr)
    integer, intent(inout) :: arr(:)

    integer :: i, n
    integer, allocatable :: iord(:)

    n = size(arr,1)
    allocate(iord(n))
    do i = 1, n
       iord(i) = i
    end do
    call qcksort_i4(arr,iord,1,n)
    arr = arr(iord)
    deallocate(iord)

  end subroutine qcksort_i4_inplace

  !> Sort a real*8 array (arr(ini:n)) in ascending order by
  !> mergesort. Returns the ordering permutation in array iord
  !> (normally, iord = ini:n). This sort is stable and uses n
  !> additional work space.
  module subroutine mergesort_r8(arr,iord,ini,n)
    real*8, dimension(:), intent(in) :: arr
    integer, dimension(:), intent(inout) :: iord
    integer, intent(in) :: ini, n

    integer :: w, i
    integer :: ileft, iright, iend
    integer :: i0, i1, k
    integer, allocatable :: iaux(:)
    logical :: doleft, useaux

    if (n <= ini) return
    w = 1
    allocate(iaux(ini:n))

    useaux = .false.
    do while (w < n)
       useaux = .not.useaux
       do i = ini, n, 2*w
          ileft = i
          iright = min(i+w,n)
          iend = min(i+2*w,n+1)

          i0 = ileft
          i1 = iright
          do k = ileft, iend-1
             doleft = (i1 >= iend)
             if (useaux) then
                if (.not.doleft) doleft = (arr(iord(i0)) <= arr(iord(i1)))
                if (i0 < iright .and. doleft) then
                   iaux(k) = iord(i0)
                   i0 = i0 + 1
                else
                   iaux(k) = iord(i1)
                   i1 = i1 + 1
                end if
             else
                if (.not.doleft) doleft = (arr(iaux(i0)) <= arr(iaux(i1)))
                if (i0 < iright .and. doleft) then
                   iord(k) = iaux(i0)
                   i0 = i0 + 1
                else
                   iord(k) = iaux(i1)
                   i1 = i1 + 1
                end if
             end if
          end do
       end do
       w = 2 * w
    end do
    if (useaux) iord(ini:n) = iaux
    deallocate(iaux)

  end subroutine mergesort_r8

  !> Sort an integer*4 array (arr(ini:n)) in ascending order by
  !> mergesort. Returns the ordering permutation in array iord
  !> (normally, iord = ini:n). This sort is stable and uses n
  !> additional work space.
  module subroutine mergesort_i4(arr, iord, ini, n)
    integer, dimension(:), intent(in) :: arr
    integer, dimension(:), intent(inout) :: iord
    integer, intent(in) :: ini, n

    integer :: w, i
    integer :: ileft, iright, iend
    integer :: i0, i1, k
    integer, allocatable :: iaux(:)
    logical :: doleft, useaux

    if (n <= ini) return
    w = 1
    allocate(iaux(ini:n))

    useaux = .false.
    do while (w < n)
       useaux = .not.useaux
       do i = ini, n, 2*w
          ileft = i
          iright = min(i+w,n)
          iend = min(i+2*w,n+1)

          i0 = ileft
          i1 = iright
          do k = ileft, iend-1
             doleft = (i1 >= iend)
             if (useaux) then
                if (.not.doleft) doleft = (arr(iord(i0)) <= arr(iord(i1)))
                if (i0 < iright .and. doleft) then
                   iaux(k) = iord(i0)
                   i0 = i0 + 1
                else
                   iaux(k) = iord(i1)
                   i1 = i1 + 1
                end if
             else
                if (.not.doleft) doleft = (arr(iaux(i0)) <= arr(iaux(i1)))
                if (i0 < iright .and. doleft) then
                   iord(k) = iaux(i0)
                   i0 = i0 + 1
                else
                   iord(k) = iaux(i1)
                   i1 = i1 + 1
                end if
             end if
          end do
       end do
       w = 2 * w
    end do
    if (useaux) iord(ini:n) = iaux
    deallocate(iaux)

  end subroutine mergesort_i4

end submodule proc
