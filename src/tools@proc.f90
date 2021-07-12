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

  !> Return the atom type for TINKER's tiny force field (see tiny.prm
  !> in the TINKER distribution). iz is the atomic number and nn is
  !> the number of neighbors. Return 0 if the atom is unknown.
  module function tiny_atom_type(iz,nn)
    integer, intent(in) :: iz, nn
    integer :: tiny_atom_type

    if (iz == 1 .and. nn == 0) then
       tiny_atom_type = 10 ! atom 10 H "H Atom" 1 1.000 0
    elseif (iz == 1 .and. nn == 1) then
       tiny_atom_type = 11 ! atom 11 H "H Monovalent" 1 1.000 1
    elseif (iz == 1 .and. nn == 2) then
       tiny_atom_type = 12 ! atom 12 H "H Divalent" 1 1.000 2
    elseif (iz == 2 .and. nn == 0) then
       tiny_atom_type = 20 ! atom 20 He "He Atom" 2 4.000 0
    elseif (iz == 5 .and. nn == 3) then
       tiny_atom_type = 53 ! atom 53 B "B Trivalent" 5 10.800 3
    elseif (iz == 5 .and. nn == 4) then
       tiny_atom_type = 54 ! atom 54 B "B Tetravalent" 5 10.800 4
    elseif (iz == 6 .and. nn == 2) then
       tiny_atom_type = 62 ! atom 62 C "C Divalent" 6 12.000 2
    elseif (iz == 6 .and. nn == 3) then
       tiny_atom_type = 63 ! atom 63 C "C Trivalent" 6 12.000 3
    elseif (iz == 6 .and. nn == 4) then
       tiny_atom_type = 64 ! atom 64 C "C Tetravalent" 6 12.000 4
    elseif (iz == 7 .and. nn == 1) then
       tiny_atom_type = 71 ! atom 71 N "N Monovalent" 7 14.000 1
    elseif (iz == 7 .and. nn == 2) then
       tiny_atom_type = 72 ! atom 72 N "N Divalent" 7 14.000 2
    elseif (iz == 7 .and. nn == 3) then
       tiny_atom_type = 73 ! atom 73 N "N Trivalent" 7 14.000 3
    elseif (iz == 7 .and. nn == 4) then
       tiny_atom_type = 74 ! atom 74 N "N+ Tetravalent" 7 14.000 4
    elseif (iz == 8 .and. nn == 1) then
       tiny_atom_type = 81 ! atom 81 O "O Monovalent" 8 16.000 1
    elseif (iz == 8 .and. nn == 2) then
       tiny_atom_type = 82 ! atom 82 O "O Divalent" 8 16.000 2
    elseif (iz == 8 .and. nn == 3) then
       tiny_atom_type = 83 ! atom 83 O "O+ Trivalent" 8 16.000 3
    elseif (iz == 9 .and. nn == 0) then
       tiny_atom_type = 90 ! atom 90 F "F- Ion" 9 18.000 0
    elseif (iz == 9 .and. nn == 1) then
       tiny_atom_type = 91 ! atom 91 F "F Monovalent" 9 18.000 1
    elseif (iz == 10 .and. nn == 0) then
       tiny_atom_type = 100 ! atom 100 Ne "Ne Atom" 10 20.200 0
    elseif (iz == 11 .and. nn == 0) then
       tiny_atom_type = 110 ! atom 110 Na "Na+ Ion" 11 23.000 0
    elseif (iz == 12 .and. nn == 0) then
       tiny_atom_type = 120 ! atom 120 Mg "Mg+2 Ion" 12 24.300 0
    elseif (iz == 13 .and. nn == 0) then
       tiny_atom_type = 130 ! atom 130 Al "Al+3 Ion" 13 27.000 0
    elseif (iz == 14 .and. nn == 4) then
       tiny_atom_type = 144 ! atom 144 Si "Si Tetravalent" 14 28.100 4
    elseif (iz == 15 .and. nn == 3) then
       tiny_atom_type = 153 ! atom 153 P "P Trivalent" 15 31.000 3
    elseif (iz == 15 .and. nn == 4) then
       tiny_atom_type = 154 ! atom 154 P "P Tetravalent" 15 31.000 4
    elseif (iz == 15 .and. nn == 5) then
       tiny_atom_type = 155 ! atom 155 P "P Pentavalent" 15 31.000 5
    elseif (iz == 16 .and. nn == 1) then
       tiny_atom_type = 161 ! atom 161 S "S Monovalent" 16 32.000 1
    elseif (iz == 16 .and. nn == 2) then
       tiny_atom_type = 162 ! atom 162 S "S Divalent" 16 32.000 2
    elseif (iz == 16 .and. nn == 3) then
       tiny_atom_type = 163 ! atom 163 S "S Trivalent" 16 32.000 3
    elseif (iz == 16 .and. nn == 4) then
       tiny_atom_type = 164 ! atom 164 S "S Tetravalent" 16 32.000 4
    elseif (iz == 16 .and. nn == 6) then
       tiny_atom_type = 166 ! atom 166 S "S Hexavalent" 16 32.000 6
    elseif (iz == 17 .and. nn == 0) then
       tiny_atom_type = 170 ! atom 170 Cl "Cl- Ion" 17 35.500 0
    elseif (iz == 17 .and. nn == 1) then
       tiny_atom_type = 171 ! atom 171 Cl "Cl Monovalent" 17 35.500 1
    elseif (iz == 18 .and. nn == 0) then
       tiny_atom_type = 180 ! atom 180 Ar "Ar Atom" 18 39.900 0
    elseif (iz == 19 .and. nn == 0) then
       tiny_atom_type = 190 ! atom 190 K "K Ion" 19 39.100 0
    elseif (iz == 20 .and. nn == 0) then
       tiny_atom_type = 200 ! atom 200 Ca "Ca+2 Ion" 20 40.100 0
    elseif (iz == 35 .and. nn == 0) then
       tiny_atom_type = 350 ! atom 350 Br "Br- Ion" 35 79.000 0
    elseif (iz == 35 .and. nn == 1) then
       tiny_atom_type = 351 ! atom 351 Br "Br Monovalent" 35 79.000 1
    elseif (iz == 36 .and. nn == 0) then
       tiny_atom_type = 360 ! atom 360 Kr "Kr Atom" 36 83.800 0
    elseif (iz == 37 .and. nn == 0) then
       tiny_atom_type = 370 ! atom 370 Rb "Rb+ Ion" 37 85.500 0
    elseif (iz == 38 .and. nn == 0) then
       tiny_atom_type = 380 ! atom 380 Sr "Sr+2 Ion" 38 87.600 0
    elseif (iz == 53 .and. nn == 0) then
       tiny_atom_type = 530 ! atom 530 I "I- Ion" 53 126.900 0
    elseif (iz == 53 .and. nn == 1) then
       tiny_atom_type = 531 ! atom 531 I "I Monovalent" 53 126.900 1
    elseif (iz == 54 .and. nn == 0) then
       tiny_atom_type = 540 ! atom 540 Xe "Xe Atom" 54 131.300 0
    elseif (iz == 55 .and. nn == 0) then
       tiny_atom_type = 550 ! atom 550 Cs "Cs+ Ion" 55 132.900 0
    elseif (iz == 56 .and. nn == 0) then
       tiny_atom_type = 560 ! atom 560 Ba "Ba+2 Ion" 56 137.300 0
    else
       tiny_atom_type = 0
    endif

  end function tiny_atom_type

  !> Transforms the basis in x2c to the Delaunay reduced basis.
  !> x2c(:,i) are the Cartesian coordinates of lattice vector i.
  !> Return the four Delaunay vectors in crystallographic coordinates
  !> (rmat) cell, see 9.1.8 in ITC. If rbas is present, it contains
  !> the three shortest of the seven Delaunay lattice vectors that
  !> form a cell (useful to transform to one of the delaunay reduced
  !> cells).
  module subroutine delaunay_reduction(x2c,rmat,rbas)
    use tools_math, only: det3, matinv
    use tools_io, only: faterr, ferror
    real*8, intent(in) :: x2c(3,3)
    real*8, intent(out) :: rmat(3,4)
    real*8, intent(out), optional :: rbas(3,3)

    integer :: i, j, k, iord(7)
    real*8 :: c2x(3,3), sc(4,4), xstar(3,7), xlen(7), dd
    logical :: again, ok

    real*8, parameter :: eps = 1d-10

    ! calculate the c2x matrix
    c2x = x2c
    call matinv(c2x,3)

    ! build the four Delaunay vectors
    rmat = 0d0
    do i = 1, 3
       rmat(i,i) = 1d0
       rmat(:,i) = matmul(x2c,rmat(:,i))
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
             dd = det3(rbas)
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
          rbas(:,i) = matmul(c2x,rbas(:,i))
       end do
    end if

    do i = 1, 4
       rmat(:,i) = matmul(c2x,rmat(:,i))
    end do

  end subroutine delaunay_reduction

end submodule proc
