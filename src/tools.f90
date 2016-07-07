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

!> merge-sort routines contributed by:
!> Jose Luis Casals Sainz <joseluiscasalssainz@gmail.com>

!Toolbox with misc. tools.
module tools
  implicit none

  private

  !xx! sorting and uniqing
  public :: mergesort, imergesort, qcksort, iqcksort, uniqc
  !xx! system
  public :: unlink

contains

  !xx! sorting
  !> Sort the elements of the array arr in ascending order using the
  !> quicksort algorithm. iord is the initial order of data in arr and 
  !> first and last the intervals for elements to be analyzed. In the output,
  !> iord contains the final order in the array arr.
  !> WARNING! It is not stable!?!?
  subroutine qcksort (arr, iord, first, last)
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
  end subroutine qcksort

  !> Same as qcksort, for integer arrays.
  !> WARNING! It is not stable!?!?
  subroutine iqcksort (iarr, iord, first, last)
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
  end subroutine iqcksort

  !> Returns all diferent vectors v(i,3) (i = 1, n) existing in the
  !> initial set v(k,3) (k=1, n, cartesian coordinates). Only
  !> v(nini,3) through v(nfin,3) are sorted. In the output, nfin
  !> returns the index for the last element. The comparison is
  !> performed in cartesian coordinates.
  subroutine uniqc(v, nini, nfin, eps)
    use tools_io, only: ferror, faterr
    
    integer, intent(in) :: nini !< Starting index
    integer, intent(inout) :: nfin !< End index
    real*8, intent(inout), dimension(3,*) :: v !< The initial and final array
    real*8, intent(in) :: eps
  
    real*8 :: eps2, x(3), dist2
    integer :: i, j, nfin0
    logical :: found
  
    if (nfin < nini) return
  
    eps2 = eps * eps
  
    nfin0 = nini-1
    do i = nini, nfin
       found = .false.
       do j = nini, nfin0
          x = v(:,i) - v(:,j)
          dist2 = x(1)*x(1) + x(2)*x(2) + x(3)*x(3)
          if (dist2 < eps2) then
             found = .true.
             exit
          end if
       end do
       if (.not.found) then
          nfin0 = nfin0 + 1
          v(:,nfin0) = v(:,i)
       end if
    end do
    nfin = nfin0
  
  end subroutine uniqc

  !xx! crystal work

  !>.equalcrysv -- determines if to crystallographic positions are
  !> equivalent by means of a lattice translation. v1 and v2 are given
  !> in crystallographic coordinates.
  function equalcrysv(v1,v2)

    logical :: equalcrysv
    real*8, intent(in) :: v1(3) !< First vector
    real*8, intent(in) :: v2(3) !< Second vector

    real*8, parameter :: pusheps = 1d-14
    real*8, parameter :: TOLx = 1d-5
    
    real*8  :: del(3)
    integer :: i

    do i = 1, 3
       del(i) = v1(i) - v2(i)
       if (del(i) .gt. -pusheps) then
          del(i) = del(i) - int(del(i)+pusheps)
       else 
          del(i) = del(i) - int(del(i)+pusheps) + 1d0
       end if
    end do

    if ((del(1) .lt. TOLx .or. del(1) .gt. 1d0-TOLx) .and.&
        (del(2) .lt. TOLx .or. del(2) .gt. 1d0-TOLx) .and.&
        (del(3) .lt. TOLx .or. del(3) .gt. 1d0-TOLx)) then
       equalcrysv = .true.
    else
       equalcrysv = .false.
    end if

  end function equalcrysv

  !xx! arrays

  !> Contributed by:
  !> Jose Luis Casals Sainz <joseluiscasalssainz@gmail.com>
  !> merge-sort sorting of an array
  !> for performance reasons, the first 2 passes are taken
  !> out of the standard loop, and use dedicated coding.
  subroutine mergesort(xvalt, irngt)

    use tools_io

    real*8, dimension (:), intent (in) :: xvalt !< array to order
    integer, dimension (:), intent (out) :: irngt !< ordered index

    integer, dimension (size(irngt)) :: jwrkt
    integer :: lmtna, lmtnc, irng1, irng2
    integer :: nval, iind, iwrkd, iwrk, iwrkf, jinda, iinda, iindb
    real(kind(xvalt)) :: xvala, xvalb
 
    nval = min (size(xvalt), size(irngt))
    select case (nval)
    case (:0)
      return
    case (1)
      irngt (1) = 1
      return
    case default
      continue
    end select

    ! fill-in the index array, creating ordered couples
    do iind = 2, nval, 2
      if (xvalt(iind-1) <= xvalt(iind)) then
        irngt (iind-1) = iind - 1
        irngt (iind) = iind
      else
        irngt (iind-1) = iind
        irngt (iind) = iind - 1
      end if
    end do
    if (mod(nval, 2) /= 0) then
      irngt (nval) = nval
    end if

    ! we will now have ordered subsets a - b - a - b - ...
    ! and merge a and b couples into     c   -   c   - ...
    lmtna = 2
    lmtnc = 4
    
    ! first iteration. the length of the ordered subsets goes from 2 to 4
    do
      if (nval <= 2) exit
      ! loop on merges of a and b into c
      do iwrkd = 0, nval - 1, 4
        if ((iwrkd+4) > nval) then
          if ((iwrkd+2) >= nval) exit
          ! 1 2 3
          if (xvalt(irngt(iwrkd+2)) <= xvalt(irngt(iwrkd+3))) exit
            !   1 3 2
            if (xvalt(irngt(iwrkd+1)) <= xvalt(irngt(iwrkd+3))) then
              irng2 = irngt (iwrkd+2)
              irngt (iwrkd+2) = irngt (iwrkd+3)
              irngt (iwrkd+3) = irng2
              !   3 1 2
            else
              irng1 = irngt (iwrkd+1)
              irngt (iwrkd+1) = irngt (iwrkd+3)
              irngt (iwrkd+3) = irngt (iwrkd+2)
              irngt (iwrkd+2) = irng1
            end if
            exit
          end if
          ! 1 2 3 4
          if (xvalt(irngt(iwrkd+2)) <= xvalt(irngt(iwrkd+3))) cycle
          ! 1 3 x x
          if (xvalt(irngt(iwrkd+1)) <= xvalt(irngt(iwrkd+3))) then
            irng2 = irngt (iwrkd+2)
            irngt (iwrkd+2) = irngt (iwrkd+3)
            if (xvalt(irng2) <= xvalt(irngt(iwrkd+4))) then
              ! 1 3 2 4
              irngt (iwrkd+3) = irng2
            else
              ! 1 3 4 2
              irngt (iwrkd+3) = irngt (iwrkd+4)
              irngt (iwrkd+4) = irng2
            end if
            ! 3 x x x
          else
            irng1 = irngt (iwrkd+1)
            irng2 = irngt (iwrkd+2)
            irngt (iwrkd+1) = irngt (iwrkd+3)
            if (xvalt(irng1) <= xvalt(irngt(iwrkd+4))) then
              irngt (iwrkd+2) = irng1
            if (xvalt(irng2) <= xvalt(irngt(iwrkd+4))) then
              ! 3 1 2 4
              irngt (iwrkd+3) = irng2
            else
              ! 3 1 4 2
              irngt (iwrkd+3) = irngt (iwrkd+4)
              irngt (iwrkd+4) = irng2
            end if
          else
            ! 3 4 1 2
            irngt (iwrkd+2) = irngt (iwrkd+4)
            irngt (iwrkd+3) = irng1
            irngt (iwrkd+4) = irng2
          end if
        end if
      end do
      ! the cs become as and bs
      lmtna = 4
      exit
    end do
    ! iteration loop. each time, the length of the ordered subsets
    ! is doubled.
    do
      if (lmtna >= nval) exit
      iwrkf = 0
      lmtnc = 2 * lmtnc
      ! loop on merges of a and b into c
      do
        iwrk = iwrkf
        iwrkd = iwrkf + 1
        jinda = iwrkf + lmtna
        iwrkf = iwrkf + lmtnc
        if (iwrkf >= nval) then
          if (jinda >= nval) exit
          iwrkf = nval
        end if
        iinda = 1
        iindb = jinda + 1
        ! shortcut for the case when the max of a is smaller
        ! than the min of b. this line may be activated when the
        ! initial set is already close to sorted.
        !
        ! if (xvalt(irngt(jinda)) <= xvalt(irngt(iindb))) cycle
        !
        ! one steps in the c subset, that we build in the final rank array
        ! make a copy of the rank array for the merge iteration
        jwrkt (1:lmtna) = irngt (iwrkd:jinda)
        xvala = xvalt (jwrkt(iinda))
        xvalb = xvalt (irngt(iindb))
        do
          iwrk = iwrk + 1
          ! we still have unprocessed values in both a and b
          if (xvala > xvalb) then
            irngt (iwrk) = irngt (iindb)
            iindb = iindb + 1
            if (iindb > iwrkf) then
              ! only a still with unprocessed values
              irngt (iwrk+1:iwrkf) = jwrkt (iinda:lmtna)
              exit
            end if
            xvalb = xvalt (irngt(iindb))
          else
            irngt (iwrk) = jwrkt (iinda)
            iinda = iinda + 1
            if (iinda > lmtna) exit! only b still with unprocessed values
              xvala = xvalt (jwrkt(iinda))
            end if
        end do
      end do
      ! the cs become as and bs
      lmtna = 2 * lmtna
    end do

  end subroutine mergesort

  !> Contributed by:
  !> Jose Luis Casals Sainz <joseluiscasalssainz@gmail.com>
  !> merge-sort sorting of an array
  !> for performance reasons, the first 2 passes are taken
  !> out of the standard loop, and use dedicated coding.
  subroutine imergesort(xvalt, irngt)

    use tools_io

    integer, dimension (:), intent (in) :: xvalt !< array to order
    integer, dimension (:), intent (out) :: irngt !< ordered index

    integer, dimension (size(irngt)) :: jwrkt
    integer :: lmtna, lmtnc, irng1, irng2
    integer :: nval, iind, iwrkd, iwrk, iwrkf, jinda, iinda, iindb
    integer(kind(xvalt)) :: xvala, xvalb
 
    nval = min (size(xvalt), size(irngt))
    select case (nval)
    case (:0)
      return
    case (1)
      irngt (1) = 1
      return
    case default
      continue
    end select

    ! fill-in the index array, creating ordered couples
    do iind = 2, nval, 2
      if (xvalt(iind-1) <= xvalt(iind)) then
        irngt (iind-1) = iind - 1
        irngt (iind) = iind
      else
        irngt (iind-1) = iind
        irngt (iind) = iind - 1
      end if
    end do
    if (mod(nval, 2) /= 0) then
      irngt (nval) = nval
    end if

    ! we will now have ordered subsets a - b - a - b - ...
    ! and merge a and b couples into     c   -   c   - ...
    lmtna = 2
    lmtnc = 4
    
    ! first iteration. the length of the ordered subsets goes from 2 to 4
    do
      if (nval <= 2) exit
      ! loop on merges of a and b into c
      do iwrkd = 0, nval - 1, 4
        if ((iwrkd+4) > nval) then
          if ((iwrkd+2) >= nval) exit
          ! 1 2 3
          if (xvalt(irngt(iwrkd+2)) <= xvalt(irngt(iwrkd+3))) exit
            !   1 3 2
            if (xvalt(irngt(iwrkd+1)) <= xvalt(irngt(iwrkd+3))) then
              irng2 = irngt (iwrkd+2)
              irngt (iwrkd+2) = irngt (iwrkd+3)
              irngt (iwrkd+3) = irng2
              !   3 1 2
            else
              irng1 = irngt (iwrkd+1)
              irngt (iwrkd+1) = irngt (iwrkd+3)
              irngt (iwrkd+3) = irngt (iwrkd+2)
              irngt (iwrkd+2) = irng1
            end if
            exit
          end if
          ! 1 2 3 4
          if (xvalt(irngt(iwrkd+2)) <= xvalt(irngt(iwrkd+3))) cycle
          ! 1 3 x x
          if (xvalt(irngt(iwrkd+1)) <= xvalt(irngt(iwrkd+3))) then
            irng2 = irngt (iwrkd+2)
            irngt (iwrkd+2) = irngt (iwrkd+3)
            if (xvalt(irng2) <= xvalt(irngt(iwrkd+4))) then
              ! 1 3 2 4
              irngt (iwrkd+3) = irng2
            else
              ! 1 3 4 2
              irngt (iwrkd+3) = irngt (iwrkd+4)
              irngt (iwrkd+4) = irng2
            end if
            ! 3 x x x
          else
            irng1 = irngt (iwrkd+1)
            irng2 = irngt (iwrkd+2)
            irngt (iwrkd+1) = irngt (iwrkd+3)
            if (xvalt(irng1) <= xvalt(irngt(iwrkd+4))) then
              irngt (iwrkd+2) = irng1
            if (xvalt(irng2) <= xvalt(irngt(iwrkd+4))) then
              ! 3 1 2 4
              irngt (iwrkd+3) = irng2
            else
              ! 3 1 4 2
              irngt (iwrkd+3) = irngt (iwrkd+4)
              irngt (iwrkd+4) = irng2
            end if
          else
            ! 3 4 1 2
            irngt (iwrkd+2) = irngt (iwrkd+4)
            irngt (iwrkd+3) = irng1
            irngt (iwrkd+4) = irng2
          end if
        end if
      end do
      ! the cs become as and bs
      lmtna = 4
      exit
    end do
    ! iteration loop. each time, the length of the ordered subsets
    ! is doubled.
    do
      if (lmtna >= nval) exit
      iwrkf = 0
      lmtnc = 2 * lmtnc
      ! loop on merges of a and b into c
      do
        iwrk = iwrkf
        iwrkd = iwrkf + 1
        jinda = iwrkf + lmtna
        iwrkf = iwrkf + lmtnc
        if (iwrkf >= nval) then
          if (jinda >= nval) exit
          iwrkf = nval
        end if
        iinda = 1
        iindb = jinda + 1
        ! shortcut for the case when the max of a is smaller
        ! than the min of b. this line may be activated when the
        ! initial set is already close to sorted.
        !
        ! if (xvalt(irngt(jinda)) <= xvalt(irngt(iindb))) cycle
        !
        ! one steps in the c subset, that we build in the final rank array
        ! make a copy of the rank array for the merge iteration
        jwrkt (1:lmtna) = irngt (iwrkd:jinda)
        xvala = xvalt (jwrkt(iinda))
        xvalb = xvalt (irngt(iindb))
        do
          iwrk = iwrk + 1
          ! we still have unprocessed values in both a and b
          if (xvala > xvalb) then
            irngt (iwrk) = irngt (iindb)
            iindb = iindb + 1
            if (iindb > iwrkf) then
              ! only a still with unprocessed values
              irngt (iwrk+1:iwrkf) = jwrkt (iinda:lmtna)
              exit
            end if
            xvalb = xvalt (irngt(iindb))
          else
            irngt (iwrk) = jwrkt (iinda)
            iinda = iinda + 1
            if (iinda > lmtna) exit! only b still with unprocessed values
              xvala = xvalt (jwrkt(iinda))
            end if
        end do
      end do
      ! the cs become as and bs
      lmtna = 2 * lmtna
    end do

  end subroutine imergesort

  !> Delete a file
  subroutine unlink(file)
    use, intrinsic :: iso_c_binding
    character(len=*), intent(in) :: file

    interface 
       subroutine unlink_wrap(file) bind(c,name="unlink_wrap")
         use, intrinsic :: iso_c_binding
         character(kind=c_char) :: file
       end subroutine unlink_wrap
    end interface

    character(len=len_trim(file)+1) :: filec

    filec = trim(file) // C_NULL_CHAR
    call unlink_wrap(filec)

  end subroutine unlink

end module tools
