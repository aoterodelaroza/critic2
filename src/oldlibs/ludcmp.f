!>
!> @file   ludcmp.f
!>
!> @brief
!> Part of the ??? library. 
!>
      subroutine ludcmp (a, n, np, indx, d)
      integer n,np,indx(n),nmax
      real*8 d,a(np,np),tiny
c
c.....ludcmp - LU decomposition of a matrix a.
c
c.....Input parameters:
c
c     a(,) ....... matrix to decompose
c     n, np ...... real and declared dimensions of the matrix
c
c.....Output parameters:
c
c     a(,) ....... LU decomposed matrix
c     indx() ..... index vector of the LU decomposed rows
c     d .......... sign of the number of row interchanges
c
c.....Scratch space:
c
c     vv() ....... scratch space vector used for pivoting
c
c
      parameter (tiny    = 1.0d-20)
      parameter (maxkpot = 15)
      parameter (maxnf   = 85)
      parameter (nmax    = maxkpot+maxnf)
      
      integer i,imax,j,k
      real*8 aamax,dum,sum,vv(nmax)
c
c.....get the scaling information in vv
c
      d = 1d0
      do i = 1, n
         aamax = 0d0
         do j = 1, n
            if (abs(a(i,j)).gt.aamax) aamax = abs(a(i,j))
         enddo
         if (aamax.eq.0d0) stop 'ludcmp : Singular matrix !! '
         vv(i) = 1d0 / aamax
      enddo
c
c.....loop over columns
c
      do j = 1, n
         if (j.gt.1) then
            do i = 1, j-1
               sum = a(i,j)
               if (i.gt.1) then
                  do k = 1, i-1
                     sum = sum - a(i,k) * a(k,j)
                  enddo
                  a(i,j) = sum
               endif
            enddo
         endif
c
c.....search for a pivot and row loop
c
         aamax = 0d0
         imax = 0
         do i=j,n
            sum = a(i,j)
            if (j.gt.1) then
               do k = 1, j-1
                  sum = sum - a(i,k) * a(k,j)
               enddo
               a(i,j) = sum
            endif
            dum = vv(i) * abs(sum)
            if (dum.ge.aamax) then
               imax = i
               aamax = dum
            endif
         enddo
c
c.....interchange rows
c
         if (j.ne.imax) then
            do k = 1, n
               dum = a(imax,k)
               a(imax,k) = a(j,k)
               a(j,k) = dum
            enddo
            d = - d
            vv(imax) = vv(j)
         endif
         indx(j) = imax
c
c........Divide by the pivot, avoiding singular matrices
c
         if (a(j,j).eq.0d0) a(j,j) = tiny
         if (j.ne.n) then
            dum = 1d0 / a(j,j)
            do i = j+1, n
               a(i,j) = a(i,j) * dum
            enddo
         endif
      enddo
      return
      end
