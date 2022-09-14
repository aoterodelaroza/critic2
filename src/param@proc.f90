! Copyright (c) 2015-2022 Alberto Otero de la Roza <aoterodelaroza@gmail.com>,
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

! Mathematical and physical constants, color specifications and factorials.
submodule (param) proc
  implicit none

  !xx! private procedures
  ! double precision function D1MACH(I)
  ! subroutine i1mcry(A, A1, B, C, D)

contains

  !> Initialize basic variables. Use at the beginning of the run.
  module subroutine param_init()

    integer :: i
    integer :: n, clock
    integer, dimension(:), allocatable :: seed

    ! factorial matrix
    fact(0)=1d0
    do i=1,mfact
       fact(i)=i*fact(i-1)
    end do

    ! random seed
    call random_seed(size = n)
    allocate(seed(n))
    call system_clock(count=clock)
    seed = clock + 37 * (/ (i - 1, i = 1, n) /)
    call random_seed(put = seed)
    deallocate(seed)

    ! machine constants
    do i = 1, 5
       d1mach_(i) = d1mach(i)
    end do

    ! wien
    ! define site-symmetric transformation coefficients
    c_kub=0.0d0
    c_kub(0,0)=1.d0
    c_kub(3,2)=1.d0
    c_kub(4,0)=.5d0*sqrt(7.d0/3.d0)
    c_kub(4,4)=.5*sqrt(5.d0/3.d0)
    c_kub(6,0)=.5d0*sqrt(.5d0)
    c_kub(6,2)=.25d0*sqrt(11.d0)
    c_kub(6,4)=-.5d0*sqrt(7.d0/2.d0)
    c_kub(6,6)=-.25d0*sqrt(5.d0)
    c_kub(7,2)=.5d0*sqrt(13.d0/6.d0)
    c_kub(7,6)=.5d0*sqrt(11.d0/6.d0)
    c_kub(8,0)=.125d0*sqrt(33.d0)
    c_kub(8,4)=.25d0*sqrt(7.d0/3.d0)
    c_kub(8,8)=.125d0*sqrt(65.d0/3.d0)
    c_kub(9,2)=.25d0*sqrt(3.d0)
    c_kub(9,4)=.5d0*sqrt(17.d0/6.d0)
    c_kub(9,6)=-.25d0*sqrt(13.d0)
    c_kub(9,8)=-.5d0*sqrt(7.d0/6.d0)
    c_kub(10,0)=.125*sqrt(65.D0/6.D0)
    c_kub(10,2)=.125*sqrt(247.D0/6.D0)
    c_kub(10,4)=-.25*sqrt(11.D0/2.D0)
    c_kub(10,6)=0.0625d0*sqrt(19.D0/3.D0)
    c_kub(10,8)=-.125*sqrt(187.D0/6.D0)
    c_kub(10,10)=-.0625d0*sqrt(85.d0)

    ! initialize the constants in the variable hash
    call vh%put("pi",pi)
    call vh%put("e",cte)
    call vh%put("eps",epsilon(1d0))

  end subroutine param_init

  ! from quadpack.
  double precision function D1MACH(I)
      INTEGER I
!
!  DOUBLE-PRECISION MACHINE CONSTANTS
!  D1MACH( 1) = B**(EMIN-1), THE SMALLEST POSITIVE MAGNITUDE.
!  D1MACH( 2) = B**EMAX*(1 - B**(-T)), THE LARGEST MAGNITUDE.
!  D1MACH( 3) = B**(-T), THE SMALLEST RELATIVE SPACING.
!  D1MACH( 4) = B**(1-T), THE LARGEST RELATIVE SPACING.
!  D1MACH( 5) = LOG10(B)
!
      INTEGER SMALL(2)
      INTEGER LARGE(2)
      INTEGER RIGHT(2)
      INTEGER DIVER(2)
      INTEGER LOG10(2)
      INTEGER SC, CRAY1(38), J
      COMMON /D9MACH/ CRAY1
      SAVE SMALL, LARGE, RIGHT, DIVER, LOG10, SC
      DOUBLE PRECISION DMACH(5)
      EQUIVALENCE (DMACH(1),SMALL(1))
      EQUIVALENCE (DMACH(2),LARGE(1))
      EQUIVALENCE (DMACH(3),RIGHT(1))
      EQUIVALENCE (DMACH(4),DIVER(1))
      EQUIVALENCE (DMACH(5),LOG10(1))
!  THIS VERSION ADAPTS AUTOMATICALLY TO MOST CURRENT MACHINES.
!  R1MACH CAN HANDLE AUTO-DOUBLE COMPILING, BUT THIS VERSION OF
!  D1MACH DOES NOT, BECAUSE WE DO NOT HAVE QUAD CONSTANTS FOR
!  MANY MACHINES YET.
!  TO COMPILE ON OLDER MACHINES, ADD A C IN COLUMN 1
!  ON THE NEXT LINE
      DATA SC/0/
!  AND REMOVE THE C FROM COLUMN 1 IN ONE OF THE SECTIONS BELOW.
!  CONSTANTS FOR EVEN OLDER MACHINES CAN BE OBTAINED BY
!          mail netlib@research.bell-labs.com
!          send old1mach from blas
!  PLEASE SEND CORRECTIONS TO dmg OR ehg@bell-labs.com.
!
!     MACHINE CONSTANTS FOR THE HONEYWELL DPS 8/70 SERIES.
!      DATA SMALL(1),SMALL(2) / O402400000000, O000000000000 /
!      DATA LARGE(1),LARGE(2) / O376777777777, O777777777777 /
!      DATA RIGHT(1),RIGHT(2) / O604400000000, O000000000000 /
!      DATA DIVER(1),DIVER(2) / O606400000000, O000000000000 /
!      DATA LOG10(1),LOG10(2) / O776464202324, O117571775714 /, SC/987/
!
!     MACHINE CONSTANTS FOR PDP-11 FORTRANS SUPPORTING
!     32-BIT INTEGERS.
!      DATA SMALL(1),SMALL(2) /    8388608,           0 /
!      DATA LARGE(1),LARGE(2) / 2147483647,          -1 /
!      DATA RIGHT(1),RIGHT(2) /  612368384,           0 /
!      DATA DIVER(1),DIVER(2) /  620756992,           0 /
!      DATA LOG10(1),LOG10(2) / 1067065498, -2063872008 /, SC/987/
!
!     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES.
!      DATA SMALL(1),SMALL(2) / O000040000000, O000000000000 /
!      DATA LARGE(1),LARGE(2) / O377777777777, O777777777777 /
!      DATA RIGHT(1),RIGHT(2) / O170540000000, O000000000000 /
!      DATA DIVER(1),DIVER(2) / O170640000000, O000000000000 /
!      DATA LOG10(1),LOG10(2) / O177746420232, O411757177572 /, SC/987/
!
!     ON FIRST CALL, IF NO DATA UNCOMMENTED, TEST MACHINE TYPES.
      IF (SC .NE. 987) THEN
         DMACH(1) = 1.D13
         IF (      SMALL(1) .EQ. 1117925532 &
            .AND. SMALL(2) .EQ. -448790528) THEN
!           *** IEEE BIG ENDIAN ***
            SMALL(1) = 1048576
            SMALL(2) = 0
            LARGE(1) = 2146435071
            LARGE(2) = -1
            RIGHT(1) = 1017118720
            RIGHT(2) = 0
            DIVER(1) = 1018167296
            DIVER(2) = 0
            LOG10(1) = 1070810131
            LOG10(2) = 1352628735
         ELSE IF ( SMALL(2) .EQ. 1117925532 &
            .AND. SMALL(1) .EQ. -448790528) THEN
!           *** IEEE LITTLE ENDIAN ***
            SMALL(2) = 1048576
            SMALL(1) = 0
            LARGE(2) = 2146435071
            LARGE(1) = -1
            RIGHT(2) = 1017118720
            RIGHT(1) = 0
            DIVER(2) = 1018167296
            DIVER(1) = 0
            LOG10(2) = 1070810131
            LOG10(1) = 1352628735
         ELSE IF ( SMALL(1) .EQ. -2065213935 &
            .AND. SMALL(2) .EQ. 10752) THEN
!               *** VAX WITH D_FLOATING ***
            SMALL(1) = 128
            SMALL(2) = 0
            LARGE(1) = -32769
            LARGE(2) = -1
            RIGHT(1) = 9344
            RIGHT(2) = 0
            DIVER(1) = 9472
            DIVER(2) = 0
            LOG10(1) = 546979738
            LOG10(2) = -805796613
         ELSE IF ( SMALL(1) .EQ. 1267827943 &
            .AND. SMALL(2) .EQ. 704643072) THEN
!               *** IBM MAINFRAME ***
            SMALL(1) = 1048576
            SMALL(2) = 0
            LARGE(1) = 2147483647
            LARGE(2) = -1
            RIGHT(1) = 856686592
            RIGHT(2) = 0
            DIVER(1) = 873463808
            DIVER(2) = 0
            LOG10(1) = 1091781651
            LOG10(2) = 1352628735
         ELSE IF ( SMALL(1) .EQ. 1120022684 &
            .AND. SMALL(2) .EQ. -448790528) THEN
!           *** CONVEX C-1 ***
            SMALL(1) = 1048576
            SMALL(2) = 0
            LARGE(1) = 2147483647
            LARGE(2) = -1
            RIGHT(1) = 1019215872
            RIGHT(2) = 0
            DIVER(1) = 1020264448
            DIVER(2) = 0
            LOG10(1) = 1072907283
            LOG10(2) = 1352628735
         ELSE IF ( SMALL(1) .EQ. 815547074 &
            .AND. SMALL(2) .EQ. 58688) THEN
!           *** VAX G-FLOATING ***
            SMALL(1) = 16
            SMALL(2) = 0
            LARGE(1) = -32769
            LARGE(2) = -1
            RIGHT(1) = 15552
            RIGHT(2) = 0
            DIVER(1) = 15568
            DIVER(2) = 0
            LOG10(1) = 1142112243
            LOG10(2) = 2046775455
         ELSE
            DMACH(2) = 1.D27 + 1
            DMACH(3) = 1.D27
            LARGE(2) = LARGE(2) - RIGHT(2)
            IF (LARGE(2) .EQ. 64 .AND. SMALL(2) .EQ. 0) THEN
               CRAY1(1) = 67291416
               DO 10 J = 1, 20
                  CRAY1(J+1) = CRAY1(J) + CRAY1(J)
 10               CONTINUE
               CRAY1(22) = CRAY1(21) + 321322
               DO 20 J = 22, 37
                  CRAY1(J+1) = CRAY1(J) + CRAY1(J)
 20               CONTINUE
               IF (CRAY1(38) .EQ. SMALL(1)) THEN
!                  *** CRAY ***
                  CALL I1MCRY(SMALL(1), J, 8285, 8388608, 0)
                  SMALL(2) = 0
                  CALL I1MCRY(LARGE(1), J, 24574, 16777215, 16777215)
                  CALL I1MCRY(LARGE(2), J, 0, 16777215, 16777214)
                  CALL I1MCRY(RIGHT(1), J, 16291, 8388608, 0)
                  RIGHT(2) = 0
                  CALL I1MCRY(DIVER(1), J, 16292, 8388608, 0)
                  DIVER(2) = 0
                  CALL I1MCRY(LOG10(1), J, 16383, 10100890, 8715215)
                  CALL I1MCRY(LOG10(2), J, 0, 16226447, 9001388)
               ELSE
                  WRITE(*,9000)
                  STOP 779
                  END IF
            ELSE
               WRITE(*,9000)
               STOP 779
               END IF
            END IF
         SC = 987
         END IF
!    SANITY CHECK
      IF (DMACH(4) .GE. 1.0D0) STOP 778
      IF (I .LT. 1 .OR. I .GT. 5) THEN
         WRITE(*,*) 'D1MACH(I): I =',I,' is out of bounds.'
         STOP
         END IF
      D1MACH = DMACH(I)
      RETURN
 9000 FORMAT(/' Adjust D1MACH by uncommenting data statements'/ &
         ' appropriate for your machine.')
! /* Standard C source for D1MACH -- remove the * in column 1 */
!#include <stdio.h>
!#include <float.h>
!#include <math.h>
!double d1mach_(long *i)
!{
!	switch(*i){
!	  case 1: return DBL_MIN;
!	  case 2: return DBL_MAX;
!	  case 3: return DBL_EPSILON/FLT_RADIX;
!	  case 4: return DBL_EPSILON;
!	  case 5: return log10((double)FLT_RADIX);
!	  }
!	fprintf(stderr, "invalid argument: d1mach(%ld)\n", *i);
!	exit(1); return 0; /* some compilers demand return values */
!}
  end function d1mach

  SUBROUTINE I1MCRY(A, A1, B, C, D)
!*** SPECIAL COMPUTATION FOR OLD CRAY MACHINES ****
      INTEGER A, A1, B, C, D
      A1 = 16777216*B + C
      A = 16777216*A1 + D
  end subroutine i1mcry

end submodule proc
