! Part of the following code has been adapted from the elk distribution, 
! version 1.3.2.
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! Distributed under the terms of the GNU General Public License.
! See each specific routine for individual licensing terms.

! Contains code from the fftpack5 library
! Licensed under the GNU General Public License (GPL).
! Copyright (C) 1995-2004, Scientific Computing Division,
! University Corporation for Atmospheric Research
! by Paul Swarztrauber and Richard Valent

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

subroutine cfftnd(nd,n,sgn,c)
!
! DESCRIPTION:
!  In-place fast Fourier transform for complex arrays in $n_d$ dimensions. The
!  forward transform is scaled by one over the size of the array. Uses a
!  modified version of the FFTPACK5 library.
!
! INPUT/OUTPUT PARAMETERS:
!   nd  : number of dimensions (in,integer)
!   n   : mesh size (in,integer(nd))
!   sgn : FFT direction, -1: forward, 1: backward (in,integer)
!   c   : array to transform (inout,complex(n(1)*n(2)*...*n(nd)))
!
!  Copyright (C) 2005 J. K. Dewhurst
!  Distributed under the terms of the GNU General Public License.
!  See the file COPYING for license details.
!
implicit none
! arguments
integer, intent(in) :: nd
integer, intent(in) :: n(nd)
integer, intent(in) :: sgn
complex(8), intent(inout) :: c(*)
! local variables
integer i,j,k,l,p,q,iw,iw1,ni
integer lensav,lenwrk
! allocatable arrays
real(8), allocatable :: wsave(:)
real(8), allocatable :: work(:)
if (nd.le.0) then
  write(*,*)
  write(*,'("Error(cfftnd): invalid number of dimensions : ",I8)') nd
  write(*,*)
  stop
end if
p=1
lensav=1
do i=1,nd
  if (n(i).le.0) then
    write(*,*)
    write(*,'("Error(cfftnd): invalid n : ",I8)') n(i)
    write(*,'(" for dimension ",I4)') i
    write(*,*)
    stop
  end if
  p=p*n(i)
  lensav=max(lensav,2*n(i)+int(log(real(n(i),8)))+4)
end do
lenwrk=2*p
allocate(wsave(lensav))
allocate(work(lenwrk))
if (sgn.gt.0) then
  q=1
  do i=1,nd
    ni=n(i)
    if (ni.gt.1) then
      iw=ni+ni+1
      iw1=iw+1
      p=p/ni
      call cfftmi(ni,wsave,lensav)
      j=1
      k=q*ni
      do l=1,p
        call cmfm1b(q,1,ni,q,c(j),work,wsave,wsave(iw),wsave(iw1))
        j=j+k
      end do
      q=k
    end if
  end do
else
  q=1
  do i=1,nd
    ni=n(i)
    if (ni.gt.1) then
      iw=ni+ni+1
      iw1=iw+1
      p=p/ni
      call cfftmi(ni,wsave,lensav)
      j=1
      k=q*ni
      do l=1,p
        call cmfm1f(q,1,ni,q,c(j),work,wsave,wsave(iw),wsave(iw1))
        j=j+k
      end do
      q=k
    end if
  end do
end if
deallocate(wsave,work)
return
end subroutine

subroutine cfftmi ( n, wsave, lensav )

!*******************************************************************************
!
!! CFFTMI: initialization for CFFTMB and CFFTMF.
!
!  Discussion:
!
!    CFFTMI initializes array WSAVE for use in its companion routines 
!    CFFTMB and CFFTMF.  CFFTMI must be called before the first call 
!    to CFFTMB or CFFTMF, and after whenever the value of integer N changes. 
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    24 March 2005
!
!  Author:
!
!    Paul Swarztrauber
!    Richard Valent
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
!    Input, integer N, the length of each sequence to be transformed. 
!    The transform is most efficient when N is a product of small primes. 
!
!    Input, integer LENSAV, the dimension of the WSAVE array.  LENSAV must be
!    at least 2*N + INT(LOG(REAL(N))) + 4. 
!
!    Output, real WSAVE(LENSAV), containing the prime factors of N and
!    also containing certain trigonometric values which will be used in
!    routines CFFTMB or CFFTMF.
!
!
  implicit none

  integer lensav

  integer iw1
  integer n
  real(8) wsave(lensav)

  if ( n == 1 ) then
    return
  end if

  iw1 = n + n + 1
  call mcfti1 ( n, wsave, wsave(iw1), wsave(iw1+1) )

  return
end

subroutine mcfti1 ( n, wa, fnf, fac )

!*******************************************************************************
!
!! MCFTI1 is an FFTPACK5 auxiliary routine.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Author:
!
!    Paul Swarztrauber
!    Richard Valent
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  real(8) fac(*)
  real(8) fnf
  integer ido
  integer ip
  integer iw
  integer k1
  integer l1
  integer l2
  integer n
  integer nf
  real(8) wa(*)
!
!  Get the factorization of N.
!
  call factor ( n, nf, fac )
  fnf = real ( nf, 8 )
  iw = 1
  l1 = 1
!
!  Set up the trigonometric tables.
!
  do k1 = 1, nf
    ip = int ( fac(k1) )
    l2 = l1 * ip
    ido = n / l2
    call tables ( ido, ip, wa(iw) )
    iw = iw + ( ip - 1 ) * ( ido + ido )
    l1 = l2
  end do

  return
end

subroutine factor ( n, nf, fac )

!*******************************************************************************
!
!! FACTOR determines the factors of an integer.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    28 March 2005
!
!  Author:
!
!    Paul Swarztrauber
!    Richard Valent
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
!    Input, integer N, the number for which factorization and other information
!    is needed.
!
!    Output, integer NF, the number of factors.
!
!    Output, real FAC(*), a list of factors of N.
!
  implicit none

  real(8) fac(*)
  integer j
  integer n
  integer nf
  integer nl
  integer nq
  integer nr
  integer ntry

  nl = n
  nf = 0
  j = 0

  do while ( 1 < nl )

    j = j + 1

    if ( j == 1 ) then
      ntry = 4
    else if ( j == 2 ) then
      ntry = 2
    else if ( j == 3 ) then
      ntry = 3
    else if ( j == 4 ) then
      ntry = 5
    else
      ntry = ntry + 2
    end if

    do

      nq = nl / ntry
      nr = nl - ntry * nq

      if ( nr /= 0 ) then
        exit
      end if

      nf = nf + 1
      fac(nf) = ntry
      nl = nq

    end do

  end do

  return
end

subroutine tables ( ido, ip, wa )

!*******************************************************************************
!
!! TABLES computes trigonometric tables needed by the FFT routines.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Author:
!
!    Paul Swarztrauber
!    Richard Valent
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  integer ido
  integer ip

  real(8) arg1
  real(8) arg2
  real(8) arg3
  real(8) arg4
  real(8) argz
  integer i
  integer j
  real(8) tpi
  real(8) wa(ido,ip-1,2)

  tpi = 8.0E+00_8 * atan ( 1.0E+00_8 )
  argz = tpi / real ( ip, 8 )
  arg1 = tpi / real ( ido * ip, 8 )

  do j = 2, ip

    arg2 = real ( j - 1, 8 ) * arg1

    do i = 1, ido
      arg3 = real ( i - 1, 8 ) * arg2
      wa(i,j-1,1) = cos ( arg3 )
      wa(i,j-1,2) = sin ( arg3 )
    end do

    if ( 5 < ip ) then
      arg4 = real ( j - 1, 8 ) * argz
      wa(1,j-1,1) = cos ( arg4 )
      wa(1,j-1,2) = sin ( arg4 )
    end if

  end do

  return
  end

subroutine cmfm1b ( lot, jump, n, inc, c, ch, wa, fnf, fac )

!*******************************************************************************
!
!! CMFM1B is an FFTPACK5 auxiliary routine.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Author:
!
!    Paul Swarztrauber
!    Richard Valent
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  complex(8) c(*)
  real(8) ch(*)
  real(8) fac(*)
  real(8) fnf
  integer ido
  integer inc
  integer ip
  integer iw
  integer jump
  integer k1
  integer l1
  integer l2
  integer lid
  integer lot
  integer n
  integer na
  integer nbr
  integer nf
  real(8) wa(*)

  nf = int ( fnf )
  na = 0
  l1 = 1
  iw = 1

  do k1 = 1, nf

    ip = int ( fac(k1) )
    l2 = ip * l1
    ido = n / l2
    lid = l1 * ido
    nbr = 1 + na + 2 * min ( ip - 2, 4 )

    if ( nbr == 1 ) then
      call cmf2kb ( lot, ido, l1, na, c, jump, inc, ch, 1, lot, wa(iw) )
    else if ( nbr == 2 ) then
      call cmf2kb ( lot, ido, l1, na, ch, 1, lot, c, jump, inc, wa(iw) )
    else if ( nbr == 3 ) then
      call cmf3kb ( lot, ido, l1, na, c, jump, inc, ch, 1, lot, wa(iw) )
    else if ( nbr == 4 ) then
      call cmf3kb ( lot, ido, l1, na, ch, 1, lot, c, jump, inc, wa(iw) )
    else if ( nbr == 5 ) then
      call cmf4kb ( lot, ido, l1, na, c, jump, inc, ch, 1, lot, wa(iw) )
    else if ( nbr == 6 ) then
      call cmf4kb ( lot, ido, l1, na, ch, 1, lot, c, jump, inc, wa(iw) )
    else if ( nbr == 7 ) then
      call cmf5kb ( lot, ido, l1, na, c, jump, inc, ch, 1, lot, wa(iw) )
    else if ( nbr == 8 ) then
      call cmf5kb ( lot, ido, l1, na, ch, 1, lot, c, jump, inc, wa(iw) )
    else if ( nbr == 9 ) then
      call cmfgkb ( lot, ido, ip, l1, lid, na, c, c, jump, inc, ch, ch, &
        1, lot, wa(iw) )
    else if ( nbr == 10 ) then
      call cmfgkb ( lot, ido, ip, l1, lid, na, ch, ch, 1, lot, c, c, &
        jump, inc, wa(iw) )
    end if

    l1 = l2
    iw = iw + ( ip - 1 ) * ( ido + ido )

    if ( ip <= 5 ) then
      na = 1 - na
    end if

  end do

  return
end

subroutine cmfm1f ( lot, jump, n, inc, c, ch, wa, fnf, fac )

!*******************************************************************************
!
!! CMFM1F is an FFTPACK5 auxiliary routine.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Author:
!
!    Paul Swarztrauber
!    Richard Valent
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  complex(8) c(*)
  real(8) ch(*)
  real(8) fac(*)
  real(8) fnf
  integer ido
  integer inc
  integer ip
  integer iw
  integer jump
  integer k1
  integer l1
  integer l2
  integer lid
  integer lot
  integer n
  integer na
  integer nbr
  integer nf
  real(8) wa(*)

  nf = int ( fnf )
  na = 0
  l1 = 1
  iw = 1

  do k1 = 1, nf

    ip = int ( fac(k1) )
    l2 = ip * l1
    ido = n / l2
    lid = l1 * ido
    nbr = 1 + na + 2 * min ( ip - 2, 4 )

    if ( nbr == 1 ) then
      call cmf2kf ( lot, ido, l1, na, c, jump, inc, ch, 1, lot, wa(iw) )
    else if ( nbr == 2 ) then
      call cmf2kf ( lot, ido, l1, na, ch, 1, lot, c, jump, inc, wa(iw) )
    else if ( nbr == 3 ) then
      call cmf3kf ( lot, ido, l1, na, c, jump, inc, ch, 1, lot, wa(iw) )
    else if ( nbr == 4 ) then
      call cmf3kf ( lot, ido, l1, na, ch, 1, lot, c, jump, inc, wa(iw) )
    else if ( nbr == 5 ) then
      call cmf4kf ( lot, ido, l1, na, c, jump, inc, ch, 1, lot, wa(iw) )
    else if ( nbr == 6 ) then
      call cmf4kf ( lot, ido, l1, na, ch, 1, lot, c, jump, inc, wa(iw) )
    else if ( nbr == 7 ) then
      call cmf5kf ( lot, ido, l1, na, c, jump, inc, ch, 1, lot, wa(iw) )
    else if ( nbr == 8 ) then
      call cmf5kf ( lot, ido, l1, na, ch, 1, lot, c, jump, inc, wa(iw) )
    else if ( nbr == 9 ) then
      call cmfgkf ( lot, ido, ip, l1, lid, na, c, c, jump, inc, ch, ch, &
        1, lot, wa(iw) )
    else if ( nbr == 10 ) then
      call cmfgkf ( lot, ido, ip, l1, lid, na, ch, ch, 1, lot, c, c, &
        jump, inc, wa(iw) )
    end if

    l1 = l2
    iw = iw + ( ip - 1 ) * ( ido + ido )

    if ( ip <= 5 ) then
      na = 1 - na
    end if

  end do

  return
end

subroutine cmf2kb ( lot, ido, l1, na, cc, im1, in1, ch, im2, in2, wa )

!*******************************************************************************
!
!! CMF2KB is an FFTPACK5 auxiliary routine.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Author:
!
!    Paul Swarztrauber
!    Richard Valent
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  integer ido
  integer in1
  integer in2
  integer l1

  real(8) cc(2,in1,l1,ido,2)
  real(8) ch(2,in2,l1,2,ido)
  real(8) chold1
  real(8) chold2
  integer i
  integer im1
  integer im2
  integer k
  integer lot
  integer m1
  integer m1d
  integer m2
  integer m2s
  integer na
  real(8) ti2
  real(8) tr2
  real(8) wa(ido,1,2)

  m1d = ( lot - 1 ) * im1 + 1
  m2s = 1 - im2

  if ( 1 < ido .or. na == 1 ) then

    do k = 1, l1
      m2 = m2s
      do m1 = 1, m1d, im1
        m2 = m2 + im2
        ch(1,m2,k,1,1) = cc(1,m1,k,1,1) + cc(1,m1,k,1,2)
        ch(1,m2,k,2,1) = cc(1,m1,k,1,1) - cc(1,m1,k,1,2)
        ch(2,m2,k,1,1) = cc(2,m1,k,1,1) + cc(2,m1,k,1,2)
        ch(2,m2,k,2,1) = cc(2,m1,k,1,1) - cc(2,m1,k,1,2)
      end do
    end do

    do i = 2, ido
      do k = 1, l1
        m2 = m2s
        do m1 = 1, m1d, im1
          m2 = m2 + im2
          ch(1,m2,k,1,i) = cc(1,m1,k,i,1)+cc(1,m1,k,i,2)
          tr2 = cc(1,m1,k,i,1)-cc(1,m1,k,i,2)
          ch(2,m2,k,1,i) = cc(2,m1,k,i,1)+cc(2,m1,k,i,2)
          ti2 = cc(2,m1,k,i,1)-cc(2,m1,k,i,2)
          ch(2,m2,k,2,i) = wa(i,1,1)*ti2+wa(i,1,2)*tr2
          ch(1,m2,k,2,i) = wa(i,1,1)*tr2-wa(i,1,2)*ti2
        end do
      end do
    end do

  else

    do k = 1, l1
      do m1 = 1, m1d, im1

        chold1         = cc(1,m1,k,1,1) + cc(1,m1,k,1,2)
        cc(1,m1,k,1,2) = cc(1,m1,k,1,1) - cc(1,m1,k,1,2)
        cc(1,m1,k,1,1) = chold1

        chold2         = cc(2,m1,k,1,1) + cc(2,m1,k,1,2)
        cc(2,m1,k,1,2) = cc(2,m1,k,1,1) - cc(2,m1,k,1,2)
        cc(2,m1,k,1,1) = chold2

      end do
    end do

  end if

  return
end

subroutine cmf2kf ( lot, ido, l1, na, cc, im1, in1, ch, im2, in2, wa )

!*******************************************************************************
!
!! CMF2KF is an FFTPACK5 auxiliary routine.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Author:
!
!    Paul Swarztrauber
!    Richard Valent
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  integer ido
  integer in1
  integer in2
  integer l1

  real(8) cc(2,in1,l1,ido,2)
  real(8) ch(2,in2,l1,2,ido)
  real(8) chold1
  real(8) chold2
  integer i
  integer im1
  integer im2
  integer k
  integer lot
  integer m1
  integer m1d
  integer m2
  integer m2s
  integer na
  real(8) sn
  real(8) ti2
  real(8) tr2
  real(8) wa(ido,1,2)

  m1d = ( lot - 1 ) * im1 + 1
  m2s = 1 - im2

  if ( 1 < ido ) then

    do k = 1, l1
      m2 = m2s
      do m1 = 1, m1d, im1
        m2 = m2 + im2
        ch(1,m2,k,1,1) = cc(1,m1,k,1,1)+cc(1,m1,k,1,2)
        ch(1,m2,k,2,1) = cc(1,m1,k,1,1)-cc(1,m1,k,1,2)
        ch(2,m2,k,1,1) = cc(2,m1,k,1,1)+cc(2,m1,k,1,2)
        ch(2,m2,k,2,1) = cc(2,m1,k,1,1)-cc(2,m1,k,1,2)
      end do
    end do

    do i = 2, ido
      do k = 1, l1
        m2 = m2s
        do m1 = 1, m1d, im1
          m2 = m2 + im2
          ch(1,m2,k,1,i) = cc(1,m1,k,i,1)+cc(1,m1,k,i,2)
          tr2 = cc(1,m1,k,i,1)-cc(1,m1,k,i,2)
          ch(2,m2,k,1,i) = cc(2,m1,k,i,1)+cc(2,m1,k,i,2)
          ti2 = cc(2,m1,k,i,1)-cc(2,m1,k,i,2)
          ch(2,m2,k,2,i) = wa(i,1,1)*ti2-wa(i,1,2)*tr2
          ch(1,m2,k,2,i) = wa(i,1,1)*tr2+wa(i,1,2)*ti2
        end do
      end do
    end do

  else if ( na == 1 ) then

    sn = 1.0E+00_8 / real ( 2 * l1, 8 )

    do k = 1, l1
      m2 = m2s
      do m1 = 1, m1d, im1
        m2 = m2 + im2
        ch(1,m2,k,1,1) = sn * ( cc(1,m1,k,1,1) + cc(1,m1,k,1,2) )
        ch(1,m2,k,2,1) = sn * ( cc(1,m1,k,1,1) - cc(1,m1,k,1,2) )
        ch(2,m2,k,1,1) = sn * ( cc(2,m1,k,1,1) + cc(2,m1,k,1,2) )
        ch(2,m2,k,2,1) = sn * ( cc(2,m1,k,1,1) - cc(2,m1,k,1,2) )
      end do
    end do

  else

    sn = 1.0E+00_8 / real ( 2 * l1, 8 )

    do k = 1, l1
      do m1 = 1, m1d, im1

        chold1         = sn * ( cc(1,m1,k,1,1) + cc(1,m1,k,1,2) )
        cc(1,m1,k,1,2) = sn * ( cc(1,m1,k,1,1) - cc(1,m1,k,1,2) )
        cc(1,m1,k,1,1) = chold1

        chold2         = sn * ( cc(2,m1,k,1,1) + cc(2,m1,k,1,2) )
        cc(2,m1,k,1,2) = sn * ( cc(2,m1,k,1,1) - cc(2,m1,k,1,2) )
        cc(2,m1,k,1,1) = chold2

      end do
    end do

  end if

  return
end

subroutine cmf3kb ( lot, ido, l1, na, cc, im1, in1, ch, im2, in2, wa )

!*******************************************************************************
!
!! CMF3KB is an FFTPACK5 auxiliary routine.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Author:
!
!    Paul Swarztrauber
!    Richard Valent
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  integer ido
  integer in1
  integer in2
  integer l1

  real(8) cc(2,in1,l1,ido,3)
  real(8) ch(2,in2,l1,3,ido)
  real(8) ci2
  real(8) ci3
  real(8) cr2
  real(8) cr3
  real(8) di2
  real(8) di3
  real(8) dr2
  real(8) dr3
  integer i
  integer im1
  integer im2
  integer k
  integer lot
  integer m1
  integer m1d
  integer m2
  integer m2s
  integer na
  real(8), parameter :: taui =  0.866025403784439E+00_8
  real(8), parameter :: taur = -0.5E+00_8
  real(8) ti2
  real(8) tr2
  real(8) wa(ido,2,2)

  m1d = ( lot - 1 ) * im1 + 1
  m2s = 1 - im2

  if ( 1 < ido .or. na == 1 ) then
    do k = 1, l1
      m2 = m2s
      do m1 = 1, m1d, im1
        m2 = m2 + im2
        tr2 = cc(1,m1,k,1,2)+cc(1,m1,k,1,3)
        cr2 = cc(1,m1,k,1,1)+taur*tr2
        ch(1,m2,k,1,1) = cc(1,m1,k,1,1)+tr2
        ti2 = cc(2,m1,k,1,2)+cc(2,m1,k,1,3)
        ci2 = cc(2,m1,k,1,1)+taur*ti2
        ch(2,m2,k,1,1) = cc(2,m1,k,1,1)+ti2
        cr3 = taui * (cc(1,m1,k,1,2)-cc(1,m1,k,1,3))
        ci3 = taui * (cc(2,m1,k,1,2)-cc(2,m1,k,1,3))
        ch(1,m2,k,2,1) = cr2-ci3
        ch(1,m2,k,3,1) = cr2+ci3
        ch(2,m2,k,2,1) = ci2+cr3
        ch(2,m2,k,3,1) = ci2-cr3
      end do
    end do

    do i = 2, ido
      do k = 1, l1
        m2 = m2s
        do m1 = 1, m1d, im1
          m2 = m2 + im2
          tr2 = cc(1,m1,k,i,2)+cc(1,m1,k,i,3)
          cr2 = cc(1,m1,k,i,1)+taur*tr2
          ch(1,m2,k,1,i) = cc(1,m1,k,i,1)+tr2
          ti2 = cc(2,m1,k,i,2)+cc(2,m1,k,i,3)
          ci2 = cc(2,m1,k,i,1)+taur*ti2
          ch(2,m2,k,1,i) = cc(2,m1,k,i,1)+ti2
          cr3 = taui*(cc(1,m1,k,i,2)-cc(1,m1,k,i,3))
          ci3 = taui*(cc(2,m1,k,i,2)-cc(2,m1,k,i,3))
          dr2 = cr2-ci3
          dr3 = cr2+ci3
          di2 = ci2+cr3
          di3 = ci2-cr3
          ch(2,m2,k,2,i) = wa(i,1,1)*di2+wa(i,1,2)*dr2
          ch(1,m2,k,2,i) = wa(i,1,1)*dr2-wa(i,1,2)*di2
          ch(2,m2,k,3,i) = wa(i,2,1)*di3+wa(i,2,2)*dr3
          ch(1,m2,k,3,i) = wa(i,2,1)*dr3-wa(i,2,2)*di3
        end do
      end do
    end do

  else

    do k = 1, l1
      do m1 = 1, m1d, im1
        tr2 = cc(1,m1,k,1,2)+cc(1,m1,k,1,3)
        cr2 = cc(1,m1,k,1,1)+taur*tr2
        cc(1,m1,k,1,1) = cc(1,m1,k,1,1)+tr2
        ti2 = cc(2,m1,k,1,2)+cc(2,m1,k,1,3)
        ci2 = cc(2,m1,k,1,1)+taur*ti2
        cc(2,m1,k,1,1) = cc(2,m1,k,1,1)+ti2
        cr3 = taui*(cc(1,m1,k,1,2)-cc(1,m1,k,1,3))
        ci3 = taui*(cc(2,m1,k,1,2)-cc(2,m1,k,1,3))
        cc(1,m1,k,1,2) = cr2-ci3
        cc(1,m1,k,1,3) = cr2+ci3
        cc(2,m1,k,1,2) = ci2+cr3
        cc(2,m1,k,1,3) = ci2-cr3
      end do
    end do

  end if

  return
end

subroutine cmf3kf ( lot, ido, l1, na, cc, im1, in1, ch, im2, in2, wa )

!*******************************************************************************
!
!! CMF3KF is an FFTPACK5 auxiliary routine.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Author:
!
!    Paul Swarztrauber
!    Richard Valent
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  integer ido
  integer in1
  integer in2
  integer l1

  real(8) cc(2,in1,l1,ido,3)
  real(8) ch(2,in2,l1,3,ido)
  real(8) ci2
  real(8) ci3
  real(8) cr2
  real(8) cr3
  real(8) di2
  real(8) di3
  real(8) dr2
  real(8) dr3
  integer i
  integer im1
  integer im2
  integer k
  integer lot
  integer m1
  integer m1d
  integer m2
  integer m2s
  integer na
  real(8) sn
  real(8), parameter :: taui = -0.866025403784439E+00_8
  real(8), parameter :: taur = -0.5E+00_8
  real(8) ti2
  real(8) tr2
  real(8) wa(ido,2,2)

  m1d = ( lot - 1 ) * im1 + 1
  m2s = 1 - im2

  if ( 1 < ido ) then

    do k = 1, l1
      m2 = m2s
      do m1 = 1, m1d, im1
        m2 = m2 + im2
        tr2 = cc(1,m1,k,1,2)+cc(1,m1,k,1,3)
        cr2 = cc(1,m1,k,1,1)+taur*tr2
        ch(1,m2,k,1,1) = cc(1,m1,k,1,1)+tr2
        ti2 = cc(2,m1,k,1,2)+cc(2,m1,k,1,3)
        ci2 = cc(2,m1,k,1,1)+taur*ti2
        ch(2,m2,k,1,1) = cc(2,m1,k,1,1)+ti2
        cr3 = taui*(cc(1,m1,k,1,2)-cc(1,m1,k,1,3))
        ci3 = taui*(cc(2,m1,k,1,2)-cc(2,m1,k,1,3))
        ch(1,m2,k,2,1) = cr2-ci3
        ch(1,m2,k,3,1) = cr2+ci3
        ch(2,m2,k,2,1) = ci2+cr3
        ch(2,m2,k,3,1) = ci2-cr3
      end do
    end do

    do i = 2, ido
      do k = 1, l1
        m2 = m2s
        do m1 = 1, m1d, im1
          m2 = m2 + im2
          tr2 = cc(1,m1,k,i,2)+cc(1,m1,k,i,3)
          cr2 = cc(1,m1,k,i,1)+taur*tr2
          ch(1,m2,k,1,i) = cc(1,m1,k,i,1)+tr2
          ti2 = cc(2,m1,k,i,2)+cc(2,m1,k,i,3)
          ci2 = cc(2,m1,k,i,1)+taur*ti2
          ch(2,m2,k,1,i) = cc(2,m1,k,i,1)+ti2
          cr3 = taui*(cc(1,m1,k,i,2)-cc(1,m1,k,i,3))
          ci3 = taui*(cc(2,m1,k,i,2)-cc(2,m1,k,i,3))
          dr2 = cr2-ci3
          dr3 = cr2+ci3
          di2 = ci2+cr3
          di3 = ci2-cr3
          ch(2,m2,k,2,i) = wa(i,1,1)*di2-wa(i,1,2)*dr2
          ch(1,m2,k,2,i) = wa(i,1,1)*dr2+wa(i,1,2)*di2
          ch(2,m2,k,3,i) = wa(i,2,1)*di3-wa(i,2,2)*dr3
          ch(1,m2,k,3,i) = wa(i,2,1)*dr3+wa(i,2,2)*di3
        end do
      end do
    end do

  else if ( na == 1 ) then

    sn = 1.0E+00_8 / real ( 3 * l1, 8 )

    do k = 1, l1
      m2 = m2s
      do m1 = 1, m1d, im1
        m2 = m2 + im2
        tr2 = cc(1,m1,k,1,2)+cc(1,m1,k,1,3)
        cr2 = cc(1,m1,k,1,1)+taur*tr2
        ch(1,m2,k,1,1) = sn*(cc(1,m1,k,1,1)+tr2)
        ti2 = cc(2,m1,k,1,2)+cc(2,m1,k,1,3)
        ci2 = cc(2,m1,k,1,1)+taur*ti2
        ch(2,m2,k,1,1) = sn*(cc(2,m1,k,1,1)+ti2)
        cr3 = taui*(cc(1,m1,k,1,2)-cc(1,m1,k,1,3))
        ci3 = taui*(cc(2,m1,k,1,2)-cc(2,m1,k,1,3))
        ch(1,m2,k,2,1) = sn*(cr2-ci3)
        ch(1,m2,k,3,1) = sn*(cr2+ci3)
        ch(2,m2,k,2,1) = sn*(ci2+cr3)
        ch(2,m2,k,3,1) = sn*(ci2-cr3)
      end do
    end do

  else

    sn = 1.0E+00_8 / real ( 3 * l1, 8 )

    do k = 1, l1
      do m1 = 1, m1d, im1
        tr2 = cc(1,m1,k,1,2)+cc(1,m1,k,1,3)
        cr2 = cc(1,m1,k,1,1)+taur*tr2
        cc(1,m1,k,1,1) = sn*(cc(1,m1,k,1,1)+tr2)
        ti2 = cc(2,m1,k,1,2)+cc(2,m1,k,1,3)
        ci2 = cc(2,m1,k,1,1)+taur*ti2
        cc(2,m1,k,1,1) = sn*(cc(2,m1,k,1,1)+ti2)
        cr3 = taui*(cc(1,m1,k,1,2)-cc(1,m1,k,1,3))
        ci3 = taui*(cc(2,m1,k,1,2)-cc(2,m1,k,1,3))
        cc(1,m1,k,1,2) = sn*(cr2-ci3)
        cc(1,m1,k,1,3) = sn*(cr2+ci3)
        cc(2,m1,k,1,2) = sn*(ci2+cr3)
        cc(2,m1,k,1,3) = sn*(ci2-cr3)
      end do
    end do

  end if

  return
end

subroutine cmf4kb ( lot, ido, l1, na, cc, im1, in1, ch, im2, in2, wa )

!*******************************************************************************
!
!! CMF4KB is an FFTPACK5 auxiliary routine.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Author:
!
!    Paul Swarztrauber
!    Richard Valent
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  integer ido
  integer in1
  integer in2
  integer l1

  real(8) cc(2,in1,l1,ido,4)
  real(8) ch(2,in2,l1,4,ido)
  real(8) ci2
  real(8) ci3
  real(8) ci4
  real(8) cr2
  real(8) cr3
  real(8) cr4
  integer i
  integer im1
  integer im2
  integer k
  integer lot
  integer m1
  integer m1d
  integer m2
  integer m2s
  integer na
  real(8) ti1
  real(8) ti2
  real(8) ti3
  real(8) ti4
  real(8) tr1
  real(8) tr2
  real(8) tr3
  real(8) tr4
  real(8) wa(ido,3,2)

  m1d = ( lot - 1 ) * im1 + 1
  m2s = 1 - im2

  if ( 1 < ido .or. na == 1 ) then

    do k = 1, l1
      m2 = m2s
      do m1 = 1, m1d, im1
        m2 = m2 + im2
        ti1 = cc(2,m1,k,1,1)-cc(2,m1,k,1,3)
        ti2 = cc(2,m1,k,1,1)+cc(2,m1,k,1,3)
        tr4 = cc(2,m1,k,1,4)-cc(2,m1,k,1,2)
        ti3 = cc(2,m1,k,1,2)+cc(2,m1,k,1,4)
        tr1 = cc(1,m1,k,1,1)-cc(1,m1,k,1,3)
        tr2 = cc(1,m1,k,1,1)+cc(1,m1,k,1,3)
        ti4 = cc(1,m1,k,1,2)-cc(1,m1,k,1,4)
        tr3 = cc(1,m1,k,1,2)+cc(1,m1,k,1,4)
        ch(1,m2,k,1,1) = tr2+tr3
        ch(1,m2,k,3,1) = tr2-tr3
        ch(2,m2,k,1,1) = ti2+ti3
        ch(2,m2,k,3,1) = ti2-ti3
        ch(1,m2,k,2,1) = tr1+tr4
        ch(1,m2,k,4,1) = tr1-tr4
        ch(2,m2,k,2,1) = ti1+ti4
        ch(2,m2,k,4,1) = ti1-ti4
      end do
    end do

    do i = 2, ido
      do k = 1, l1
        m2 = m2s
        do m1 = 1, m1d, im1
          m2 = m2 + im2
          ti1 = cc(2,m1,k,i,1)-cc(2,m1,k,i,3)
          ti2 = cc(2,m1,k,i,1)+cc(2,m1,k,i,3)
          ti3 = cc(2,m1,k,i,2)+cc(2,m1,k,i,4)
          tr4 = cc(2,m1,k,i,4)-cc(2,m1,k,i,2)
          tr1 = cc(1,m1,k,i,1)-cc(1,m1,k,i,3)
          tr2 = cc(1,m1,k,i,1)+cc(1,m1,k,i,3)
          ti4 = cc(1,m1,k,i,2)-cc(1,m1,k,i,4)
          tr3 = cc(1,m1,k,i,2)+cc(1,m1,k,i,4)
          ch(1,m2,k,1,i) = tr2+tr3
          cr3 = tr2-tr3
          ch(2,m2,k,1,i) = ti2+ti3
          ci3 = ti2-ti3
          cr2 = tr1+tr4
          cr4 = tr1-tr4
          ci2 = ti1+ti4
          ci4 = ti1-ti4
          ch(1,m2,k,2,i) = wa(i,1,1)*cr2-wa(i,1,2)*ci2
          ch(2,m2,k,2,i) = wa(i,1,1)*ci2+wa(i,1,2)*cr2
          ch(1,m2,k,3,i) = wa(i,2,1)*cr3-wa(i,2,2)*ci3
          ch(2,m2,k,3,i) = wa(i,2,1)*ci3+wa(i,2,2)*cr3
          ch(1,m2,k,4,i) = wa(i,3,1)*cr4-wa(i,3,2)*ci4
          ch(2,m2,k,4,i) = wa(i,3,1)*ci4+wa(i,3,2)*cr4
        end do
      end do
    end do

  else

    do k = 1, l1
      do m1 = 1, m1d, im1
        ti1 = cc(2,m1,k,1,1)-cc(2,m1,k,1,3)
        ti2 = cc(2,m1,k,1,1)+cc(2,m1,k,1,3)
        tr4 = cc(2,m1,k,1,4)-cc(2,m1,k,1,2)
        ti3 = cc(2,m1,k,1,2)+cc(2,m1,k,1,4)
        tr1 = cc(1,m1,k,1,1)-cc(1,m1,k,1,3)
        tr2 = cc(1,m1,k,1,1)+cc(1,m1,k,1,3)
        ti4 = cc(1,m1,k,1,2)-cc(1,m1,k,1,4)
        tr3 = cc(1,m1,k,1,2)+cc(1,m1,k,1,4)
        cc(1,m1,k,1,1) = tr2+tr3
        cc(1,m1,k,1,3) = tr2-tr3
        cc(2,m1,k,1,1) = ti2+ti3
        cc(2,m1,k,1,3) = ti2-ti3
        cc(1,m1,k,1,2) = tr1+tr4
        cc(1,m1,k,1,4) = tr1-tr4
        cc(2,m1,k,1,2) = ti1+ti4
        cc(2,m1,k,1,4) = ti1-ti4
      end do
    end do

  end if

  return
end

subroutine cmf4kf ( lot, ido, l1, na, cc, im1, in1, ch, im2, in2, wa )

!*******************************************************************************
!
!! CMF4KF is an FFTPACK5 auxiliary routine.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Author:
!
!    Paul Swarztrauber
!    Richard Valent
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  integer ido
  integer in1
  integer in2
  integer l1

  real(8) cc(2,in1,l1,ido,4)
  real(8) ch(2,in2,l1,4,ido)
  real(8) ci2
  real(8) ci3
  real(8) ci4
  real(8) cr2
  real(8) cr3
  real(8) cr4
  integer i
  integer im1
  integer im2
  integer k
  integer lot
  integer m1
  integer m1d
  integer m2
  integer m2s
  integer na
  real(8) sn
  real(8) ti1
  real(8) ti2
  real(8) ti3
  real(8) ti4
  real(8) tr1
  real(8) tr2
  real(8) tr3
  real(8) tr4
  real(8) wa(ido,3,2)

  m1d = ( lot - 1 ) * im1 + 1
  m2s = 1 - im2

  if ( 1 < ido ) then

    do k = 1, l1
      m2 = m2s
      do m1 = 1, m1d, im1
        m2 = m2 + im2
        ti1 = cc(2,m1,k,1,1)-cc(2,m1,k,1,3)
        ti2 = cc(2,m1,k,1,1)+cc(2,m1,k,1,3)
        tr4 = cc(2,m1,k,1,2)-cc(2,m1,k,1,4)
        ti3 = cc(2,m1,k,1,2)+cc(2,m1,k,1,4)
        tr1 = cc(1,m1,k,1,1)-cc(1,m1,k,1,3)
        tr2 = cc(1,m1,k,1,1)+cc(1,m1,k,1,3)
        ti4 = cc(1,m1,k,1,4)-cc(1,m1,k,1,2)
        tr3 = cc(1,m1,k,1,2)+cc(1,m1,k,1,4)
        ch(1,m2,k,1,1) = tr2+tr3
        ch(1,m2,k,3,1) = tr2-tr3
        ch(2,m2,k,1,1) = ti2+ti3
        ch(2,m2,k,3,1) = ti2-ti3
        ch(1,m2,k,2,1) = tr1+tr4
        ch(1,m2,k,4,1) = tr1-tr4
        ch(2,m2,k,2,1) = ti1+ti4
        ch(2,m2,k,4,1) = ti1-ti4
      end do
    end do

    do i = 2, ido
      do k = 1, l1
        m2 = m2s
        do m1 = 1, m1d, im1
          m2 = m2 + im2
          ti1 = cc(2,m1,k,i,1)-cc(2,m1,k,i,3)
          ti2 = cc(2,m1,k,i,1)+cc(2,m1,k,i,3)
          ti3 = cc(2,m1,k,i,2)+cc(2,m1,k,i,4)
          tr4 = cc(2,m1,k,i,2)-cc(2,m1,k,i,4)
          tr1 = cc(1,m1,k,i,1)-cc(1,m1,k,i,3)
          tr2 = cc(1,m1,k,i,1)+cc(1,m1,k,i,3)
          ti4 = cc(1,m1,k,i,4)-cc(1,m1,k,i,2)
          tr3 = cc(1,m1,k,i,2)+cc(1,m1,k,i,4)
          ch(1,m2,k,1,i) = tr2+tr3
          cr3 = tr2-tr3
          ch(2,m2,k,1,i) = ti2+ti3
          ci3 = ti2-ti3
          cr2 = tr1+tr4
          cr4 = tr1-tr4
          ci2 = ti1+ti4
          ci4 = ti1-ti4
          ch(1,m2,k,2,i) = wa(i,1,1)*cr2+wa(i,1,2)*ci2
          ch(2,m2,k,2,i) = wa(i,1,1)*ci2-wa(i,1,2)*cr2
          ch(1,m2,k,3,i) = wa(i,2,1)*cr3+wa(i,2,2)*ci3
          ch(2,m2,k,3,i) = wa(i,2,1)*ci3-wa(i,2,2)*cr3
          ch(1,m2,k,4,i) = wa(i,3,1)*cr4+wa(i,3,2)*ci4
          ch(2,m2,k,4,i) = wa(i,3,1)*ci4-wa(i,3,2)*cr4
        end do
      end do
    end do

  else if ( na == 1 ) then

    sn = 1.0E+00_8 / real ( 4 * l1, 8 )

    do k = 1, l1
      m2 = m2s
      do m1 = 1, m1d, im1
        m2 = m2 + im2
        ti1 = cc(2,m1,k,1,1)-cc(2,m1,k,1,3)
        ti2 = cc(2,m1,k,1,1)+cc(2,m1,k,1,3)
        tr4 = cc(2,m1,k,1,2)-cc(2,m1,k,1,4)
        ti3 = cc(2,m1,k,1,2)+cc(2,m1,k,1,4)
        tr1 = cc(1,m1,k,1,1)-cc(1,m1,k,1,3)
        tr2 = cc(1,m1,k,1,1)+cc(1,m1,k,1,3)
        ti4 = cc(1,m1,k,1,4)-cc(1,m1,k,1,2)
        tr3 = cc(1,m1,k,1,2)+cc(1,m1,k,1,4)
        ch(1,m2,k,1,1) = sn*(tr2+tr3)
        ch(1,m2,k,3,1) = sn*(tr2-tr3)
        ch(2,m2,k,1,1) = sn*(ti2+ti3)
        ch(2,m2,k,3,1) = sn*(ti2-ti3)
        ch(1,m2,k,2,1) = sn*(tr1+tr4)
        ch(1,m2,k,4,1) = sn*(tr1-tr4)
        ch(2,m2,k,2,1) = sn*(ti1+ti4)
        ch(2,m2,k,4,1) = sn*(ti1-ti4)
      end do
    end do

  else

    sn = 1.0E+00_8 / real ( 4 * l1, 8 )

    do k = 1, l1
      do m1 = 1, m1d, im1
        ti1 = cc(2,m1,k,1,1)-cc(2,m1,k,1,3)
        ti2 = cc(2,m1,k,1,1)+cc(2,m1,k,1,3)
        tr4 = cc(2,m1,k,1,2)-cc(2,m1,k,1,4)
        ti3 = cc(2,m1,k,1,2)+cc(2,m1,k,1,4)
        tr1 = cc(1,m1,k,1,1)-cc(1,m1,k,1,3)
        tr2 = cc(1,m1,k,1,1)+cc(1,m1,k,1,3)
        ti4 = cc(1,m1,k,1,4)-cc(1,m1,k,1,2)
        tr3 = cc(1,m1,k,1,2)+cc(1,m1,k,1,4)
        cc(1,m1,k,1,1) = sn*(tr2+tr3)
        cc(1,m1,k,1,3) = sn*(tr2-tr3)
        cc(2,m1,k,1,1) = sn*(ti2+ti3)
        cc(2,m1,k,1,3) = sn*(ti2-ti3)
        cc(1,m1,k,1,2) = sn*(tr1+tr4)
        cc(1,m1,k,1,4) = sn*(tr1-tr4)
        cc(2,m1,k,1,2) = sn*(ti1+ti4)
        cc(2,m1,k,1,4) = sn*(ti1-ti4)
      end do
    end do

  end if

  return
end

subroutine cmf5kb ( lot, ido, l1, na, cc, im1, in1, ch, im2, in2, wa )

!*******************************************************************************
!
!! CMF5KB is an FFTPACK5 auxiliary routine.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Author:
!
!    Paul Swarztrauber
!    Richard Valent
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  integer ido
  integer in1
  integer in2
  integer l1

  real(8) cc(2,in1,l1,ido,5)
  real(8) ch(2,in2,l1,5,ido)
  real(8) chold1
  real(8) chold2
  real(8) ci2
  real(8) ci3
  real(8) ci4
  real(8) ci5
  real(8) cr2
  real(8) cr3
  real(8) cr4
  real(8) cr5
  real(8) di2
  real(8) di3
  real(8) di4
  real(8) di5
  real(8) dr2
  real(8) dr3
  real(8) dr4
  real(8) dr5
  integer i
  integer im1
  integer im2
  integer k
  integer lot
  integer m1
  integer m1d
  integer m2
  integer m2s
  integer na
  real(8) ti2
  real(8) ti3
  real(8) ti4
  real(8) ti5
  real(8), parameter :: ti11 =  0.9510565162951536E+00_8
  real(8), parameter :: ti12 =  0.5877852522924731E+00_8
  real(8) tr2
  real(8) tr3
  real(8) tr4
  real(8) tr5
  real(8), parameter :: tr11 =  0.3090169943749474E+00_8
  real(8), parameter :: tr12 = -0.8090169943749474E+00_8
  real(8) wa(ido,4,2)

  m1d = ( lot - 1 ) * im1 + 1
  m2s = 1 - im2

  if ( 1 < ido .or. na == 1 ) then

    do k = 1, l1
      m2 = m2s
      do m1 = 1, m1d, im1
        m2 = m2 + im2
        ti5 = cc(2,m1,k,1,2)-cc(2,m1,k,1,5)
        ti2 = cc(2,m1,k,1,2)+cc(2,m1,k,1,5)
        ti4 = cc(2,m1,k,1,3)-cc(2,m1,k,1,4)
        ti3 = cc(2,m1,k,1,3)+cc(2,m1,k,1,4)
        tr5 = cc(1,m1,k,1,2)-cc(1,m1,k,1,5)
        tr2 = cc(1,m1,k,1,2)+cc(1,m1,k,1,5)
        tr4 = cc(1,m1,k,1,3)-cc(1,m1,k,1,4)
        tr3 = cc(1,m1,k,1,3)+cc(1,m1,k,1,4)
        ch(1,m2,k,1,1) = cc(1,m1,k,1,1)+tr2+tr3
        ch(2,m2,k,1,1) = cc(2,m1,k,1,1)+ti2+ti3
        cr2 = cc(1,m1,k,1,1)+tr11*tr2+tr12*tr3
        ci2 = cc(2,m1,k,1,1)+tr11*ti2+tr12*ti3
        cr3 = cc(1,m1,k,1,1)+tr12*tr2+tr11*tr3
        ci3 = cc(2,m1,k,1,1)+tr12*ti2+tr11*ti3
        cr5 = ti11*tr5+ti12*tr4
        ci5 = ti11*ti5+ti12*ti4
        cr4 = ti12*tr5-ti11*tr4
        ci4 = ti12*ti5-ti11*ti4
        ch(1,m2,k,2,1) = cr2-ci5
        ch(1,m2,k,5,1) = cr2+ci5
        ch(2,m2,k,2,1) = ci2+cr5
        ch(2,m2,k,3,1) = ci3+cr4
        ch(1,m2,k,3,1) = cr3-ci4
        ch(1,m2,k,4,1) = cr3+ci4
        ch(2,m2,k,4,1) = ci3-cr4
        ch(2,m2,k,5,1) = ci2-cr5
      end do
    end do

    do i = 2, ido
      do k = 1, l1
        m2 = m2s
        do m1 = 1, m1d, im1
          m2 = m2 + im2
          ti5 = cc(2,m1,k,i,2)-cc(2,m1,k,i,5)
          ti2 = cc(2,m1,k,i,2)+cc(2,m1,k,i,5)
          ti4 = cc(2,m1,k,i,3)-cc(2,m1,k,i,4)
          ti3 = cc(2,m1,k,i,3)+cc(2,m1,k,i,4)
          tr5 = cc(1,m1,k,i,2)-cc(1,m1,k,i,5)
          tr2 = cc(1,m1,k,i,2)+cc(1,m1,k,i,5)
          tr4 = cc(1,m1,k,i,3)-cc(1,m1,k,i,4)
          tr3 = cc(1,m1,k,i,3)+cc(1,m1,k,i,4)
          ch(1,m2,k,1,i) = cc(1,m1,k,i,1)+tr2+tr3
          ch(2,m2,k,1,i) = cc(2,m1,k,i,1)+ti2+ti3
          cr2 = cc(1,m1,k,i,1)+tr11*tr2+tr12*tr3
          ci2 = cc(2,m1,k,i,1)+tr11*ti2+tr12*ti3
          cr3 = cc(1,m1,k,i,1)+tr12*tr2+tr11*tr3
          ci3 = cc(2,m1,k,i,1)+tr12*ti2+tr11*ti3
          cr5 = ti11*tr5+ti12*tr4
          ci5 = ti11*ti5+ti12*ti4
          cr4 = ti12*tr5-ti11*tr4
          ci4 = ti12*ti5-ti11*ti4
          dr3 = cr3-ci4
          dr4 = cr3+ci4
          di3 = ci3+cr4
          di4 = ci3-cr4
          dr5 = cr2+ci5
          dr2 = cr2-ci5
          di5 = ci2-cr5
          di2 = ci2+cr5
          ch(1,m2,k,2,i) = wa(i,1,1) * dr2 - wa(i,1,2) * di2
          ch(2,m2,k,2,i) = wa(i,1,1) * di2 + wa(i,1,2) * dr2
          ch(1,m2,k,3,i) = wa(i,2,1) * dr3 - wa(i,2,2) * di3
          ch(2,m2,k,3,i) = wa(i,2,1) * di3 + wa(i,2,2) * dr3
          ch(1,m2,k,4,i) = wa(i,3,1) * dr4 - wa(i,3,2) * di4
          ch(2,m2,k,4,i) = wa(i,3,1) * di4 + wa(i,3,2) * dr4
          ch(1,m2,k,5,i) = wa(i,4,1) * dr5 - wa(i,4,2) * di5
          ch(2,m2,k,5,i) = wa(i,4,1) * di5 + wa(i,4,2) * dr5
        end do
      end do
    end do

  else

    do k = 1, l1
      do m1 = 1, m1d, im1
        ti5 = cc(2,m1,k,1,2)-cc(2,m1,k,1,5)
        ti2 = cc(2,m1,k,1,2)+cc(2,m1,k,1,5)
        ti4 = cc(2,m1,k,1,3)-cc(2,m1,k,1,4)
        ti3 = cc(2,m1,k,1,3)+cc(2,m1,k,1,4)
        tr5 = cc(1,m1,k,1,2)-cc(1,m1,k,1,5)
        tr2 = cc(1,m1,k,1,2)+cc(1,m1,k,1,5)
        tr4 = cc(1,m1,k,1,3)-cc(1,m1,k,1,4)
        tr3 = cc(1,m1,k,1,3)+cc(1,m1,k,1,4)

        chold1 = cc(1,m1,k,1,1) + tr2 + tr3
        chold2 = cc(2,m1,k,1,1) + ti2 + ti3

        cr2 = cc(1,m1,k,1,1) + tr11 * tr2 + tr12 * tr3
        ci2 = cc(2,m1,k,1,1) + tr11 * ti2 + tr12 * ti3
        cr3 = cc(1,m1,k,1,1) + tr12 * tr2 + tr11 * tr3
        ci3 = cc(2,m1,k,1,1) + tr12 * ti2 + tr11 * ti3

        cc(1,m1,k,1,1) = chold1
        cc(2,m1,k,1,1) = chold2

        cr5 = ti11*tr5 + ti12*tr4
        ci5 = ti11*ti5 + ti12*ti4
        cr4 = ti12*tr5 - ti11*tr4
        ci4 = ti12*ti5 - ti11*ti4
        cc(1,m1,k,1,2) = cr2-ci5
        cc(1,m1,k,1,5) = cr2+ci5
        cc(2,m1,k,1,2) = ci2+cr5
        cc(2,m1,k,1,3) = ci3+cr4
        cc(1,m1,k,1,3) = cr3-ci4
        cc(1,m1,k,1,4) = cr3+ci4
        cc(2,m1,k,1,4) = ci3-cr4
        cc(2,m1,k,1,5) = ci2-cr5
      end do
    end do

  end if

  return
end

subroutine cmf5kf ( lot, ido, l1, na, cc, im1, in1, ch, im2, in2, wa )

!*******************************************************************************
!
!! CMF5KF is an FFTPACK5 auxiliary routine.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Author:
!
!    Paul Swarztrauber
!    Richard Valent
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  integer ido
  integer in1
  integer in2
  integer l1

  real(8) cc(2,in1,l1,ido,5)
  real(8) ch(2,in2,l1,5,ido)
  real(8) chold1
  real(8) chold2
  real(8) ci2
  real(8) ci3
  real(8) ci4
  real(8) ci5
  real(8) cr2
  real(8) cr3
  real(8) cr4
  real(8) cr5
  real(8) di2
  real(8) di3
  real(8) di4
  real(8) di5
  real(8) dr2
  real(8) dr3
  real(8) dr4
  real(8) dr5
  integer i
  integer im1
  integer im2
  integer k
  integer lot
  integer m1
  integer m1d
  integer m2
  integer m2s
  integer na
  real(8) sn
  real(8) ti2
  real(8) ti3
  real(8) ti4
  real(8) ti5
  real(8), parameter :: ti11 = -0.9510565162951536E+00_8
  real(8), parameter :: ti12 = -0.5877852522924731E+00_8
  real(8) tr2
  real(8) tr3
  real(8) tr4
  real(8) tr5
  real(8), parameter :: tr11 =  0.3090169943749474E+00_8
  real(8), parameter :: tr12 = -0.8090169943749474E+00_8
  real(8) wa(ido,4,2)

  m1d = ( lot - 1 ) * im1 + 1
  m2s = 1 - im2

  if ( 1 < ido ) then

    do k = 1, l1
      m2 = m2s
      do m1 = 1, m1d, im1
        m2 = m2 + im2
        ti5 = cc(2,m1,k,1,2)-cc(2,m1,k,1,5)
        ti2 = cc(2,m1,k,1,2)+cc(2,m1,k,1,5)
        ti4 = cc(2,m1,k,1,3)-cc(2,m1,k,1,4)
        ti3 = cc(2,m1,k,1,3)+cc(2,m1,k,1,4)
        tr5 = cc(1,m1,k,1,2)-cc(1,m1,k,1,5)
        tr2 = cc(1,m1,k,1,2)+cc(1,m1,k,1,5)
        tr4 = cc(1,m1,k,1,3)-cc(1,m1,k,1,4)
        tr3 = cc(1,m1,k,1,3)+cc(1,m1,k,1,4)
        ch(1,m2,k,1,1) = cc(1,m1,k,1,1)+tr2+tr3
        ch(2,m2,k,1,1) = cc(2,m1,k,1,1)+ti2+ti3
        cr2 = cc(1,m1,k,1,1)+tr11*tr2+tr12*tr3
        ci2 = cc(2,m1,k,1,1)+tr11*ti2+tr12*ti3
        cr3 = cc(1,m1,k,1,1)+tr12*tr2+tr11*tr3
        ci3 = cc(2,m1,k,1,1)+tr12*ti2+tr11*ti3
        cr5 = ti11*tr5+ti12*tr4
        ci5 = ti11*ti5+ti12*ti4
        cr4 = ti12*tr5-ti11*tr4
        ci4 = ti12*ti5-ti11*ti4
        ch(1,m2,k,2,1) = cr2-ci5
        ch(1,m2,k,5,1) = cr2+ci5
        ch(2,m2,k,2,1) = ci2+cr5
        ch(2,m2,k,3,1) = ci3+cr4
        ch(1,m2,k,3,1) = cr3-ci4
        ch(1,m2,k,4,1) = cr3+ci4
        ch(2,m2,k,4,1) = ci3-cr4
        ch(2,m2,k,5,1) = ci2-cr5
      end do
    end do

    do i = 2, ido
      do k = 1, l1
        m2 = m2s
        do m1 = 1, m1d, im1
          m2 = m2 + im2
          ti5 = cc(2,m1,k,i,2)-cc(2,m1,k,i,5)
          ti2 = cc(2,m1,k,i,2)+cc(2,m1,k,i,5)
          ti4 = cc(2,m1,k,i,3)-cc(2,m1,k,i,4)
          ti3 = cc(2,m1,k,i,3)+cc(2,m1,k,i,4)
          tr5 = cc(1,m1,k,i,2)-cc(1,m1,k,i,5)
          tr2 = cc(1,m1,k,i,2)+cc(1,m1,k,i,5)
          tr4 = cc(1,m1,k,i,3)-cc(1,m1,k,i,4)
          tr3 = cc(1,m1,k,i,3)+cc(1,m1,k,i,4)
          ch(1,m2,k,1,i) = cc(1,m1,k,i,1)+tr2+tr3
          ch(2,m2,k,1,i) = cc(2,m1,k,i,1)+ti2+ti3
          cr2 = cc(1,m1,k,i,1)+tr11*tr2+tr12*tr3
          ci2 = cc(2,m1,k,i,1)+tr11*ti2+tr12*ti3
          cr3 = cc(1,m1,k,i,1)+tr12*tr2+tr11*tr3
          ci3 = cc(2,m1,k,i,1)+tr12*ti2+tr11*ti3
          cr5 = ti11*tr5+ti12*tr4
          ci5 = ti11*ti5+ti12*ti4
          cr4 = ti12*tr5-ti11*tr4
          ci4 = ti12*ti5-ti11*ti4
          dr3 = cr3-ci4
          dr4 = cr3+ci4
          di3 = ci3+cr4
          di4 = ci3-cr4
          dr5 = cr2+ci5
          dr2 = cr2-ci5
          di5 = ci2-cr5
          di2 = ci2+cr5
          ch(1,m2,k,2,i) = wa(i,1,1)*dr2+wa(i,1,2)*di2
          ch(2,m2,k,2,i) = wa(i,1,1)*di2-wa(i,1,2)*dr2
          ch(1,m2,k,3,i) = wa(i,2,1)*dr3+wa(i,2,2)*di3
          ch(2,m2,k,3,i) = wa(i,2,1)*di3-wa(i,2,2)*dr3
          ch(1,m2,k,4,i) = wa(i,3,1)*dr4+wa(i,3,2)*di4
          ch(2,m2,k,4,i) = wa(i,3,1)*di4-wa(i,3,2)*dr4
          ch(1,m2,k,5,i) = wa(i,4,1)*dr5+wa(i,4,2)*di5
          ch(2,m2,k,5,i) = wa(i,4,1)*di5-wa(i,4,2)*dr5
        end do
      end do
    end do

  else if ( na == 1 ) then

    sn = 1.0E+00_8 / real ( 5 * l1, 8 )
 
    do k = 1, l1
      m2 = m2s
      do m1 = 1, m1d, im1
        m2 = m2 + im2
        ti5 = cc(2,m1,k,1,2)-cc(2,m1,k,1,5)
        ti2 = cc(2,m1,k,1,2)+cc(2,m1,k,1,5)
        ti4 = cc(2,m1,k,1,3)-cc(2,m1,k,1,4)
        ti3 = cc(2,m1,k,1,3)+cc(2,m1,k,1,4)
        tr5 = cc(1,m1,k,1,2)-cc(1,m1,k,1,5)
        tr2 = cc(1,m1,k,1,2)+cc(1,m1,k,1,5)
        tr4 = cc(1,m1,k,1,3)-cc(1,m1,k,1,4)
        tr3 = cc(1,m1,k,1,3)+cc(1,m1,k,1,4)
        ch(1,m2,k,1,1) = sn*(cc(1,m1,k,1,1)+tr2+tr3)
        ch(2,m2,k,1,1) = sn*(cc(2,m1,k,1,1)+ti2+ti3)
        cr2 = cc(1,m1,k,1,1)+tr11*tr2+tr12*tr3
        ci2 = cc(2,m1,k,1,1)+tr11*ti2+tr12*ti3
        cr3 = cc(1,m1,k,1,1)+tr12*tr2+tr11*tr3
        ci3 = cc(2,m1,k,1,1)+tr12*ti2+tr11*ti3
        cr5 = ti11*tr5+ti12*tr4
        ci5 = ti11*ti5+ti12*ti4
        cr4 = ti12*tr5-ti11*tr4
        ci4 = ti12*ti5-ti11*ti4
        ch(1,m2,k,2,1) = sn*(cr2-ci5)
        ch(1,m2,k,5,1) = sn*(cr2+ci5)
        ch(2,m2,k,2,1) = sn*(ci2+cr5)
        ch(2,m2,k,3,1) = sn*(ci3+cr4)
        ch(1,m2,k,3,1) = sn*(cr3-ci4)
        ch(1,m2,k,4,1) = sn*(cr3+ci4)
        ch(2,m2,k,4,1) = sn*(ci3-cr4)
        ch(2,m2,k,5,1) = sn*(ci2-cr5)
      end do
    end do

  else

    sn = 1.0E+00_8 / real ( 5 * l1, 8 )

    do k = 1, l1
      do m1 = 1, m1d, im1

        ti5 = cc(2,m1,k,1,2) - cc(2,m1,k,1,5)
        ti2 = cc(2,m1,k,1,2) + cc(2,m1,k,1,5)
        ti4 = cc(2,m1,k,1,3) - cc(2,m1,k,1,4)
        ti3 = cc(2,m1,k,1,3) + cc(2,m1,k,1,4)
        tr5 = cc(1,m1,k,1,2) - cc(1,m1,k,1,5)
        tr2 = cc(1,m1,k,1,2) + cc(1,m1,k,1,5)
        tr4 = cc(1,m1,k,1,3) - cc(1,m1,k,1,4)
        tr3 = cc(1,m1,k,1,3) + cc(1,m1,k,1,4)

        chold1 = sn * ( cc(1,m1,k,1,1) + tr2 + tr3 )
        chold2 = sn * ( cc(2,m1,k,1,1) + ti2 + ti3 )

        cr2 = cc(1,m1,k,1,1) + tr11 * tr2 + tr12 * tr3
        ci2 = cc(2,m1,k,1,1) + tr11 * ti2 + tr12 * ti3
        cr3 = cc(1,m1,k,1,1) + tr12 * tr2 + tr11 * tr3
        ci3 = cc(2,m1,k,1,1) + tr12 * ti2 + tr11 * ti3

        cc(1,m1,k,1,1) = chold1
        cc(2,m1,k,1,1) = chold2

        cr5 = ti11 * tr5 + ti12 * tr4
        ci5 = ti11 * ti5 + ti12 * ti4
        cr4 = ti12 * tr5 - ti11 * tr4
        ci4 = ti12 * ti5 - ti11 * ti4

        cc(1,m1,k,1,2) = sn * ( cr2 - ci5 )
        cc(1,m1,k,1,5) = sn * ( cr2 + ci5 )
        cc(2,m1,k,1,2) = sn * ( ci2 + cr5 )
        cc(2,m1,k,1,3) = sn * ( ci3 + cr4 )
        cc(1,m1,k,1,3) = sn * ( cr3 - ci4 )
        cc(1,m1,k,1,4) = sn * ( cr3 + ci4 )
        cc(2,m1,k,1,4) = sn * ( ci3 - cr4 )
        cc(2,m1,k,1,5) = sn * ( ci2 - cr5 )

      end do
    end do

  end if

  return
end

subroutine cmfgkb ( lot, ido, ip, l1, lid, na, cc, cc1, im1, in1, &
  ch, ch1, im2, in2, wa )

!*******************************************************************************
!
!! CMFGKB is an FFTPACK5 auxiliary routine.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Author:
!
!    Paul Swarztrauber
!    Richard Valent
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  integer ido
  integer in1
  integer in2
  integer ip
  integer l1
  integer lid

  real(8) cc(2,in1,l1,ip,ido)
  real(8) cc1(2,in1,lid,ip)
  real(8) ch(2,in2,l1,ido,ip)
  real(8) ch1(2,in2,lid,ip)
  real(8) chold1
  real(8) chold2
  integer i
  integer idlj
  integer im1
  integer im2
  integer ipp2
  integer ipph
  integer j
  integer jc
  integer k
  integer ki
  integer l
  integer lc
  integer lot
  integer m1
  integer m1d
  integer m2
  integer m2s
  integer na
  real(8) wa(ido,ip-1,2)
  real(8) wai
  real(8) war

  m1d = ( lot - 1 ) * im1 + 1
  m2s = 1 - im2
  ipp2 = ip + 2
  ipph = ( ip + 1 ) / 2

  do ki = 1, lid
    m2 = m2s
    do m1 = 1, m1d, im1
      m2 = m2 + im2
      ch1(1,m2,ki,1) = cc1(1,m1,ki,1)
      ch1(2,m2,ki,1) = cc1(2,m1,ki,1)
    end do
  end do

  do j = 2, ipph
    jc = ipp2 - j
    do ki = 1, lid
      m2 = m2s
      do m1 = 1, m1d, im1
        m2 = m2 + im2
        ch1(1,m2,ki,j) =  cc1(1,m1,ki,j) + cc1(1,m1,ki,jc)
        ch1(1,m2,ki,jc) = cc1(1,m1,ki,j) - cc1(1,m1,ki,jc)
        ch1(2,m2,ki,j) =  cc1(2,m1,ki,j) + cc1(2,m1,ki,jc)
        ch1(2,m2,ki,jc) = cc1(2,m1,ki,j) - cc1(2,m1,ki,jc)
      end do
    end do
  end do

  do j = 2, ipph
    do ki = 1, lid
      m2 = m2s
      do m1 = 1, m1d, im1
        m2 = m2 + im2
        cc1(1,m1,ki,1) = cc1(1,m1,ki,1) + ch1(1,m2,ki,j)
        cc1(2,m1,ki,1) = cc1(2,m1,ki,1) + ch1(2,m2,ki,j)
      end do
    end do
  end do

  do l = 2, ipph

    lc = ipp2 - l
    do ki = 1, lid
      m2 = m2s
      do m1 = 1, m1d, im1
        m2 = m2 + im2

        cc1(1,m1,ki,l)  = ch1(1,m2,ki,1) + wa(1,l-1,1) * ch1(1,m2,ki,2)
        cc1(1,m1,ki,lc) =                  wa(1,l-1,2) * ch1(1,m2,ki,ip)
        cc1(2,m1,ki,l)  = ch1(2,m2,ki,1) + wa(1,l-1,1) * ch1(2,m2,ki,2)
        cc1(2,m1,ki,lc) =                  wa(1,l-1,2) * ch1(2,m2,ki,ip)

      end do
    end do

    do j = 3, ipph
      jc = ipp2 - j
      idlj = mod ( ( l - 1 ) * ( j - 1 ), ip )
      war = wa(1,idlj,1)
      wai = wa(1,idlj,2)
      do ki = 1, lid
        m2 = m2s
        do m1 = 1, m1d, im1
          m2 = m2 + im2
          cc1(1,m1,ki,l)  = cc1(1,m1,ki,l)  + war * ch1(1,m2,ki,j)
          cc1(1,m1,ki,lc) = cc1(1,m1,ki,lc) + wai * ch1(1,m2,ki,jc)
          cc1(2,m1,ki,l)  = cc1(2,m1,ki,l)  + war * ch1(2,m2,ki,j)
          cc1(2,m1,ki,lc) = cc1(2,m1,ki,lc) + wai * ch1(2,m2,ki,jc)
        end do
      end do
    end do

  end do

  if( 1 < ido .or. na == 1 ) then

    do ki = 1, lid
      m2 = m2s
      do m1 = 1, m1d, im1
        m2 = m2 + im2
        ch1(1,m2,ki,1) = cc1(1,m1,ki,1)
        ch1(2,m2,ki,1) = cc1(2,m1,ki,1)
      end do
    end do

    do j = 2, ipph
      jc = ipp2 - j
      do ki = 1, lid
        m2 = m2s
        do m1 = 1, m1d, im1
          m2 = m2 + im2
          ch1(1,m2,ki,j)  = cc1(1,m1,ki,j) - cc1(2,m1,ki,jc)
          ch1(1,m2,ki,jc) = cc1(1,m1,ki,j) + cc1(2,m1,ki,jc)
          ch1(2,m2,ki,jc) = cc1(2,m1,ki,j) - cc1(1,m1,ki,jc)
          ch1(2,m2,ki,j)  = cc1(2,m1,ki,j) + cc1(1,m1,ki,jc)
        end do
      end do
    end do

    if ( ido == 1 ) then
      return
    end if

    do i = 1, ido
      do k = 1, l1
        m2 = m2s
        do m1 = 1, m1d, im1
          m2 = m2 + im2
          cc(1,m1,k,1,i) = ch(1,m2,k,i,1)
          cc(2,m1,k,1,i) = ch(2,m2,k,i,1)
        end do
      end do
    end do

    do j = 2, ip
      do k = 1, l1
        m2 = m2s
        do m1 = 1, m1d, im1
          m2 = m2 + im2
          cc(1,m1,k,j,1) = ch(1,m2,k,1,j)
          cc(2,m1,k,j,1) = ch(2,m2,k,1,j)
        end do
      end do
    end do

    do j = 2, ip
      do i = 2, ido
        do k = 1, l1
          m2 = m2s
          do m1 = 1, m1d, im1
            m2 = m2 + im2
            cc(1,m1,k,j,i) = wa(i,j-1,1) * ch(1,m2,k,i,j) &
                           - wa(i,j-1,2) * ch(2,m2,k,i,j)
            cc(2,m1,k,j,i) = wa(i,j-1,1) * ch(2,m2,k,i,j) &
                           + wa(i,j-1,2) * ch(1,m2,k,i,j)
          end do
        end do
      end do
    end do

  else

    do j = 2, ipph
      jc = ipp2 - j
      do ki = 1, lid
        do m1 = 1, m1d, im1

          chold1         = cc1(1,m1,ki,j) - cc1(2,m1,ki,jc)
          chold2         = cc1(1,m1,ki,j) + cc1(2,m1,ki,jc)
          cc1(1,m1,ki,j) = chold1

          cc1(2,m1,ki,jc) = cc1(2,m1,ki,j) - cc1(1,m1,ki,jc)
          cc1(2,m1,ki,j)  = cc1(2,m1,ki,j) + cc1(1,m1,ki,jc)
          cc1(1,m1,ki,jc) = chold2

        end do
      end do
    end do

  end if

  return
end

subroutine cmfgkf ( lot, ido, ip, l1, lid, na, cc, cc1, im1, in1, &
  ch, ch1, im2, in2, wa )

!*******************************************************************************
!
!! CMFGKF is an FFTPACK5 auxiliary routine.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Author:
!
!    Paul Swarztrauber
!    Richard Valent
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  integer ido
  integer in1
  integer in2
  integer ip
  integer l1
  integer lid

  real(8) cc(2,in1,l1,ip,ido)
  real(8) cc1(2,in1,lid,ip)
  real(8) ch(2,in2,l1,ido,ip)
  real(8) ch1(2,in2,lid,ip)
  real(8) chold1
  real(8) chold2
  integer i
  integer idlj
  integer im1
  integer im2
  integer ipp2
  integer ipph
  integer j 
  integer jc
  integer k
  integer ki
  integer l
  integer lc
  integer lot
  integer m1
  integer m1d
  integer m2
  integer m2s
  integer na
  real(8) sn
  real(8) wa(ido,ip-1,2)
  real(8) wai
  real(8) war

  m1d = ( lot - 1 ) * im1 + 1
  m2s = 1 - im2
  ipp2 = ip + 2
  ipph = ( ip + 1 ) / 2

  do ki = 1, lid
    m2 = m2s
    do m1 = 1, m1d, im1
      m2 = m2 + im2
      ch1(1,m2,ki,1) = cc1(1,m1,ki,1)
      ch1(2,m2,ki,1) = cc1(2,m1,ki,1)
    end do
  end do

  do j = 2, ipph
    jc = ipp2 - j
    do ki = 1, lid
      m2 = m2s
      do m1 = 1, m1d, im1
        m2 = m2 + im2
        ch1(1,m2,ki,j) =  cc1(1,m1,ki,j) + cc1(1,m1,ki,jc)
        ch1(1,m2,ki,jc) = cc1(1,m1,ki,j) - cc1(1,m1,ki,jc)
        ch1(2,m2,ki,j) =  cc1(2,m1,ki,j) + cc1(2,m1,ki,jc)
        ch1(2,m2,ki,jc) = cc1(2,m1,ki,j) - cc1(2,m1,ki,jc)
      end do
    end do
  end do

  do j = 2, ipph
    do ki = 1, lid
      m2 = m2s
      do m1 = 1, m1d, im1
        m2 = m2 + im2
        cc1(1,m1,ki,1) = cc1(1,m1,ki,1) + ch1(1,m2,ki,j)
        cc1(2,m1,ki,1) = cc1(2,m1,ki,1) + ch1(2,m2,ki,j)
      end do
    end do
  end do

  do l = 2, ipph

    lc = ipp2 - l

    do ki = 1, lid
      m2 = m2s
      do m1 = 1, m1d, im1
        m2 = m2 + im2
        cc1(1,m1,ki,l)  = ch1(1,m2,ki,1) + wa(1,l-1,1) * ch1(1,m2,ki,2)
        cc1(1,m1,ki,lc) =                - wa(1,l-1,2) * ch1(1,m2,ki,ip)
        cc1(2,m1,ki,l)  = ch1(2,m2,ki,1) + wa(1,l-1,1) * ch1(2,m2,ki,2)
        cc1(2,m1,ki,lc) =                - wa(1,l-1,2) * ch1(2,m2,ki,ip)
      end do
    end do

    do j = 3, ipph
      jc = ipp2 - j
      idlj = mod ( ( l - 1 ) * ( j - 1 ), ip )
      war = wa(1,idlj,1)
      wai = -wa(1,idlj,2)
      do ki = 1, lid
        m2 = m2s
        do m1 = 1, m1d, im1
          m2 = m2 + im2
          cc1(1,m1,ki,l)  = cc1(1,m1,ki,l)  + war * ch1(1,m2,ki,j)
          cc1(1,m1,ki,lc) = cc1(1,m1,ki,lc) + wai * ch1(1,m2,ki,jc)
          cc1(2,m1,ki,l)  = cc1(2,m1,ki,l)  + war * ch1(2,m2,ki,j)
          cc1(2,m1,ki,lc) = cc1(2,m1,ki,lc) + wai * ch1(2,m2,ki,jc)
        end do
      end do
    end do

  end do

  if ( 1 < ido ) then

    do ki = 1, lid
      m2 = m2s
      do m1 = 1, m1d, im1
        m2 = m2 + im2
        ch1(1,m2,ki,1) = cc1(1,m1,ki,1)
        ch1(2,m2,ki,1) = cc1(2,m1,ki,1)
      end do
    end do

    do j = 2, ipph
      jc = ipp2 - j
      do ki = 1, lid
        m2 = m2s
        do m1 = 1, m1d, im1
          m2 = m2 + im2
          ch1(1,m2,ki,j)  = cc1(1,m1,ki,j) - cc1(2,m1,ki,jc)
          ch1(2,m2,ki,j)  = cc1(2,m1,ki,j) + cc1(1,m1,ki,jc)
          ch1(1,m2,ki,jc) = cc1(1,m1,ki,j) + cc1(2,m1,ki,jc)
          ch1(2,m2,ki,jc) = cc1(2,m1,ki,j) - cc1(1,m1,ki,jc)
        end do
      end do
    end do

    do i = 1, ido
      do k = 1, l1
        m2 = m2s
        do m1 = 1, m1d, im1
          m2 = m2 + im2
          cc(1,m1,k,1,i) = ch(1,m2,k,i,1)
          cc(2,m1,k,1,i) = ch(2,m2,k,i,1)
        end do
      end do
    end do

    do j = 2, ip
      do k = 1, l1
        m2 = m2s
        do m1 = 1, m1d, im1
          m2 = m2 + im2
          cc(1,m1,k,j,1) = ch(1,m2,k,1,j)
          cc(2,m1,k,j,1) = ch(2,m2,k,1,j)
        end do
      end do
    end do

    do j = 2, ip
      do i = 2, ido
        do k = 1, l1
          m2 = m2s
          do m1 = 1, m1d, im1
            m2 = m2 + im2
            cc(1,m1,k,j,i) = wa(i,j-1,1) * ch(1,m2,k,i,j) &
                           + wa(i,j-1,2) * ch(2,m2,k,i,j)
            cc(2,m1,k,j,i) = wa(i,j-1,1) * ch(2,m2,k,i,j) &
                           - wa(i,j-1,2) * ch(1,m2,k,i,j)
          end do
        end do
      end do
    end do

  else if ( na == 1 ) then

    sn = 1.0E+00_8 / real ( ip * l1, 8 )

    do ki = 1, lid
      m2 = m2s
      do m1 = 1, m1d, im1
        m2 = m2 + im2
        ch1(1,m2,ki,1) = sn * cc1(1,m1,ki,1)
        ch1(2,m2,ki,1) = sn * cc1(2,m1,ki,1)
      end do
    end do

    do j = 2, ipph
      jc = ipp2 - j
      do ki = 1, lid
        m2 = m2s
        do m1 = 1, m1d, im1
          m2 = m2 + im2
          ch1(1,m2,ki,j)  = sn * ( cc1(1,m1,ki,j) - cc1(2,m1,ki,jc) )
          ch1(2,m2,ki,j)  = sn * ( cc1(2,m1,ki,j) + cc1(1,m1,ki,jc) )
          ch1(1,m2,ki,jc) = sn * ( cc1(1,m1,ki,j) + cc1(2,m1,ki,jc) )
          ch1(2,m2,ki,jc) = sn * ( cc1(2,m1,ki,j) - cc1(1,m1,ki,jc) )
        end do
      end do
    end do

  else

    sn = 1.0E+00_8 / real ( ip * l1, 8 )

    do ki = 1, lid
      m2 = m2s
      do m1 = 1, m1d, im1
        m2 = m2 + im2
        cc1(1,m1,ki,1) = sn * cc1(1,m1,ki,1)
        cc1(2,m1,ki,1) = sn * cc1(2,m1,ki,1)
      end do
    end do

    do j = 2, ipph
      jc = ipp2 - j
      do ki = 1, lid
        do m1 = 1, m1d, im1
          chold1 = sn * ( cc1(1,m1,ki,j) - cc1(2,m1,ki,jc) )
          chold2 = sn * ( cc1(1,m1,ki,j) + cc1(2,m1,ki,jc) )
          cc1(1,m1,ki,j) = chold1
          cc1(2,m1,ki,jc) = sn * ( cc1(2,m1,ki,j) - cc1(1,m1,ki,jc) )
          cc1(2,m1,ki,j)  = sn * ( cc1(2,m1,ki,j) + cc1(1,m1,ki,jc) )
          cc1(1,m1,ki,jc) = chold2
        end do
      end do
    end do

  end if

  return
end


