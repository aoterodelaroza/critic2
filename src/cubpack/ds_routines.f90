! This file is F-compatible, except for upper/lower case conventions.
!--------------------------------------------------------------------
MODULE DS_ROUTINES

USE Precision_Model, ONLY: stnd
Implicit NONE
private

public :: DSINIT, DSGET, DSSPUT, DSUPUT, DSSTAT, DSPINT, DSCOPY
public :: DSSUM, DSFREE, DSUSED


CONTAINS
!-----------------------------------------------------------------------
      SUBROUTINE DSINIT(DIMENS,NRVERT,NIINFO,NRINFO,NRFUNC,NRVACA,      &
                        BOTTIS,BOTTRS,ISTORE,IFAIL)
!***BEGIN PROLOGUE DSINIT
!***DATE WRITTEN   900612   (YYMMDD)
!***REVISION DATE  970612   (YYMMDD)
!***AUTHOR
!          Ronald Cools, Dept. of Computer Science,
!          Katholieke Universiteit Leuven, Celestijnenlaan 200A,
!          B-3001 Heverlee, Belgium
!          Email:  Ronald.Cools@cs.kuleuven.ac.be
!
!***PURPOSE  To initialise a data structure that stores information
!            about the subregions in an adaptive integrator.
!***DESCRIPTION
!
! Organization of data structure:
! -------------------------------
!
!  For each subregion in the data structure a record with the
!  following information is stored:
!  - a sortkey ( 1 double precision number )
!  - the vertices to describe the subregion
!             ( NRVERT*DIMENS double precision numbers )
!  - approximation for the integrals over the subregion
!             ( NRFUNC double precision numbers )
!  - estimates for the errors of the approximations
!             ( NRFUNC double precision numbers )
!  - additional information about the subregion
!             ( NIINFO integers and NRINFO double precision numbers)
!
!  These records appear in a partially sorted binary tree
!  with the element with the largest sortkey on top or in an
!  unsorted pool.
!  All information is kept in 2 arrays: one for the integers
!  and one for the double precision numbers. The records have
!  a fixed place. Only the pointers to these records are
!  modified if a new element is added or if the top is removed.
!  If the topelement is removed from the tree, the place
!  were the record was stored can be reused.
!  The pointer to this vacant place is saved. NRVACA integers
!  are reserved for this purpose.
!
!  The following figure shows how all information is distributed
!  over the heaps.
!
!                    ISTORE                 RSTORE
!                ---------------        ---------------
!                | 1) DIMENS   |        | ----------- |
!                | 2) NRVERT   |        | | sortkey | |
!                | 3) NIINFO   |        | | err est | |
!                | 4) NRINFO   |        | | int appr| |
!                | 5) NRFUNC   |        | | vertices| |
!                | 6) NRVACA   |        | |  info   | |
!                | 7) BOTTIS   |        | ----------- |
!                | 8) BOTTRS   |        | | sortkey | |
!                | 9) OFFSET   |        | | err est | |
!                | 10) START   |        | | int appr| |
!                | 11) BLOCK   |        | | vertices| |
!                | 12) INTREE  |        | |  info   | |
!                | 13) INPOOL  |        | ----------- |
!                | 14) HOLES   |        |      .      |
!                | 15) LOST    |        |      .      |
!                |=============|        |      .      |
!                |       .     |        |             |
!                |       .     |        |             |
!                | ----------- |        |             |
!                | |vacancies| |        |             |
!     offset   ->| ----------- |        |             |
!                |=============|        |             |
!  offset + 1  ->| ----------- |        |             |
!                | | pointer | |        |             |
!                | | to tree | |        |             |
!                | ----------- |        |             |
!                |      .      |        |             |
!                |      .      |        |             |
!                |      .      |        |             |
!                |             |        |             |
!                |             |        |             |
!                |             |        |             |
!                |             |        |             |
!                |      .      |        |             |
!                |      .      |        |             |
!                |      .      |        |             |
!                | ----------- |        |             |
!                | |   pool  | |        |             |
!                | ----------- |        |             |
!                |=============|        |             |
!     start    ->|      .      |        |             |
!                |      .      |        |             |
!                |      .      |        |             |
!                | ----------- |        |             |
!                | |  region | |        |             |
!                | |  info   | |        |             |
!       bottis ->| ----------- |        |             |<- bottrs
!                ---------------        ---------------
!
!  The maximum number of subregions that can be stored in this
!  data structure is
!  min( BOTTRS/BLOCK , (BOTTIS-CONST-NRVACA)/(1+NRINFO) )
!
! Input parameters
! ----------------
!
!  DIMENS = dimension of (sub-)regions
!           DIMENS > 0
!  NRVERT = number of vertices to describe a subregion
!           NRVERT > DIMENS
!  NIINFO = number of integers used to save additional information
!           for each subregion
!           NIINFO > 0
!  NRINFO = number of double precision numbers used to save additional
!           information for each subregion
!           NRINFO > 0
!  NRFUNC = number of integrand functions for which information must
!           be stored
!           NRFUNC > 0
!  NRVACA = maximum number of pointers to empty spaces in the heap
!           that must be stored for later use
!           NRVACA >= 0
!  BOTTIS = length of integer array ISTORE.
!           Needed because this heap is also filled up starting from
!           the bottom
!           bottis > 0
!  BOTTRS = length of the double precision heap
!           Only needed for checks
!           bottrs > 0
!  ISTORE = integer array of dimension (BOTTIS)
!           Used to store integer part of records, pointers to
!           records and information about the data structure.
!
! Output parameter
! ----------------
!
!  IFAIL = integer to indicate success or failure
!          IFAIL = 0 for normal exit
!          IFAIL = 10001 if DIMENS is less than 1
!          IFAIL = 10002 if NRVERT is less than or equal to dimens
!          IFAIL = 10003 if NIINFO is less than 1
!          IFAIL = 10004 if NRINFO is less than 1
!          IFAIL = 10005 if NRFUNC is less than 1
!          IFAIL = 10006 if NRVACA is negative
!          IFAIL = 10007 if the array of integers cannot contain one record
!          IFAIL = 10008 if the array of reals cannot contain one record
!  ISTORE = contains initialisation information about the data structure
!
!***LONG DESCRIPTION
!
!  Other procedures available to work on this data structure are:
!  SUBROUTINE DSGET
!  SUBROUTINE DSSPUT
!  SUBROUTINE DSUPUT
!  SUBROUTINE DSSTAT
!  SUBROUTINE DSPINT
!  SUBROUTINE DSSUM
!  INTEGER FUNCTION DSFREE
!  INTEGER FUNCTION DSUSED
!
!***END PROLOGUE DSINIT
!
!  Global variables
!
      INTEGER, INTENT(IN) :: DIMENS,NRVERT,NIINFO,NRINFO,NRFUNC,NRVACA,  &
                             BOTTIS,BOTTRS
      INTEGER, INTENT(OUT) :: IFAIL
      INTEGER, DIMENSION(:), INTENT(IN OUT) :: ISTORE
!
!  Local variables and constants
!
!  CONST = the number of constants and information variables during
!          the existence of the data structure
!  INTREE = number of subregions in the sorted tree
!  INPOOL = number of subregions in the unsorted pool
!  HOLES  = number of pointers to holes in the heap
!           0 <= HOLES <= NRVACA
!  LOST = number of records that cannot be accessed any more
!  OFFSET : The first pointer to the tree is stored in ISTORE(OFFSET+1)
!  START : The first pointer to the pool is stored in ISTORE(START-1)
!  BLOCK = number of double precision numbers in a record
!
      INTEGER, PARAMETER ::  CONST = 15, INTREE = 0, INPOOL = 0, &
                             HOLES = 0,  LOST = 0
      INTEGER ::    OFFSET,BLOCK,START
!***FIRST EXECUTABLE STATEMENT
!
! Check input
!
      BLOCK = NRINFO + NRVERT*DIMENS + 2*NRFUNC
      IF (NRFUNC /= 1) THEN
          BLOCK = BLOCK + 1
      END IF
      IF (DIMENS <= 0) THEN
          IFAIL = 10001
      ELSE IF (NRVERT <= DIMENS) THEN
          IFAIL = 10002
      ELSE IF (NIINFO <= 0) THEN
          IFAIL = 10003
      ELSE IF (NRINFO <= 0) THEN
          IFAIL = 10004
      ELSE IF (NRFUNC <= 0) THEN
          IFAIL = 10005
      ELSE IF (NRVACA < 0) THEN
          IFAIL = 10006
      ELSE IF (BOTTIS-NRVACA-CONST < 1+NIINFO) THEN
          IFAIL = 10007
      ELSE IF (BOTTRS < BLOCK) THEN
          IFAIL = 10008
      ELSE
          ISTORE(1) = DIMENS
          ISTORE(2) = NRVERT
          ISTORE(3) = NIINFO
          ISTORE(4) = NRINFO
          ISTORE(5) = NRFUNC
          ISTORE(6) = NRVACA
          ISTORE(7) = BOTTIS
          ISTORE(8) = BOTTRS
          OFFSET = CONST + NRVACA
          ISTORE(9) = OFFSET
          START = BOTTIS + 1 - ((BOTTIS-NRVACA-CONST)/ (1+NIINFO))*NIINFO
          ISTORE(10) = START
          ISTORE(11) = BLOCK
          ISTORE(12) = INTREE
          ISTORE(13) = INPOOL
          ISTORE(14) = HOLES
          ISTORE(15) = LOST
          IFAIL = 0
      END IF
!
!***END DSINIT
!
      RETURN
      END SUBROUTINE DSINIT
!-----------------------------------------------------------------------
 SUBROUTINE DSSPUT(VERTIC,INTAPP,ERREST,IRGINF,RRGINF,ISTORE,RSTORE,IFAIL)
!
!***BEGIN PROLOGUE DSSPUT
!***DATE WRITTEN   900612   (YYMMDD)
!***REVISION DATE  970612   (YYMMDD)
!***AUTHOR
!          Ronald Cools, Dept. of Computer Science,
!          Katholieke Universiteit Leuven, Celestijnenlaan 200A,
!          B-3001 Heverlee, Belgium
!          Email:  Ronald.Cools@cs.kuleuven.ac.be
!
!***PURPOSE  To add a record to a data structure that stores information
!            about the subregions in an adaptive integrator.
!            The record is added to a sorted tree.
!***DESCRIPTION
!   See the description of SUBROUTINE DSINIT
!
!  Input parameters
!  ----------------
!
!  VERTIC = double precision array of dimension (DIMENS,NRVERT)
!           Contains the vertices that describe a subregion.
!           ( VERTIC(1,i),...,VERTIC(DIMENS,i) ) are the coordinates
!           of the i-th vertex.
!  INTAPP = double precision array of dimension (NRFUNC)
!           Contains approximations to the integrals.
!  ERREST = double precision array of dimension (NRFUNC)
!           Contains error estimates.
!  IRGINF = integer array of dimension (NIINFO).
!           Contains additional information about the subregion
!  RRGINF = double precision array of dimension (NRINFO).
!           Contains additional information about the subregion
!  ISTORE = integer array of dimension (BOTTIS).
!           Contains the integers of the records.
!  RSTORE = double precision array of dimension (BOTTRS).
!           Contains the double precision numbers of the records.
!
!  Output parameters
!  -----------------
!
!  IFAIL = integer to indicate success or failure
!          IFAIL = 0 for normal exit
!          IFAIL = 10009 if the integer array ISTORE is full
!          IFAIL = 10010 if the double precision array RSTORE is full
!  ISTORE = integer array of dimension (BOTTIS).
!           Contains the integers of the records.
!  RSTORE = double precision array of dimension (BOTTRS).
!           Contains the double precision numbers of the records.
!
!***END PROLOGUE DSSPUT
!
! Variables and constants used by DS-procedures only
! --------------------------------------------------
!
!  DIMENS = dimension of (sub-)regions
!  NRVERT = number of vertices to describe a subregion
!  NIINFO = number of integers used to save additional information
!           for each subregion
!  NRINFO = number of double precision numbers used to save additional
!           information for each subregion
!  NRFUNC = number of integrand functions for which information must
!           be stored
!  NRVACA = maximum number of pointers to empty spaces in the heap
!           that must be stored for later use
!  BOTTIS = length of integer array ISTORE.
!           Needed because this heap is also filled up starting from
!           the bottom
!  BOTTRS = length of the double precision heap
!           Only needed for checks
!  INTREE = number of subregions in the sorted tree
!  INPOOL = number of subregions in the unsorted pool
!  HOLES  = number of pointers to holes in the heap
!           0 <= HOLES <= NRVACA
!  LOST = number of records that cannot be accessed any more
!  OFFSET : The first pointer to the tree is stored in ISTORE(OFFSET+1)
!  START : The first pointer to the pool is stored in ISTORE(START-1)
!  BLOCK = number of double precision numbers in a record
!
!  Global variables
!
      INTEGER, DIMENSION(:), INTENT(IN) :: IRGINF
      INTEGER, DIMENSION(:), INTENT(IN OUT) :: ISTORE
      INTEGER, INTENT(OUT)   :: IFAIL
      REAL(kind=stnd), DIMENSION(:), INTENT(IN) :: INTAPP,ERREST,RRGINF
      REAL(kind=stnd), DIMENSION(:), INTENT(IN OUT) :: RSTORE
      REAL(kind=stnd), DIMENSION(:,:), INTENT(IN) :: VERTIC
!
! Local variables
!
! SORKEY = sortkey for maintaining the partially sorted tree
!        = maximum(errest(1),...,errest(nrfunc))
! SPACE  = index to the place were a new record can be
!          inserted in the heap
!
      INTEGER :: DIMENS,NRVERT,NIINFO,NRINFO,NRFUNC,BOTTIS,BOTTRS,  &
                 INTREE,INPOOL,HOLES,LOST,BLOCK,OFFSET,START
      REAL(kind=stnd) :: SORKEY
      INTEGER :: I,SPACE,POINT,SUBRGN,SUBTMP
!***FIRST EXECUTABLE STATEMENT
!
!  The following initialisation statements are included
!  to make the code readable
!
!  Initialise data structure constants
      DIMENS = ISTORE(1)
      NRVERT = ISTORE(2)
      NIINFO = ISTORE(3)
      NRINFO = ISTORE(4)
      NRFUNC = ISTORE(5)
      BOTTIS = ISTORE(7)
      BOTTRS = ISTORE(8)
      OFFSET = ISTORE(9)
      START = ISTORE(10)
      BLOCK = ISTORE(11)
!  Initialise data structure variables
      INTREE = ISTORE(12)
      INPOOL = ISTORE(13)
      HOLES = ISTORE(14)
      LOST = ISTORE(15)
!
! Check if enough space is left in the arrays to put in a subregion
!
      IF (HOLES <= 0) THEN
          IF ((BOTTIS-START+1- (INTREE+INPOOL+LOST)*NIINFO <          &
              NIINFO) .OR. (START-OFFSET-1-INTREE-INPOOL < 1)) THEN
              IFAIL = 10009
              RETURN
          ELSE IF (BOTTRS- (INTREE+INPOOL+LOST)*BLOCK < BLOCK) THEN
              IFAIL = 10010
              RETURN
          END IF
      END IF
      INTREE = INTREE + 1
!
! Determine index for new record
!
      IF (HOLES <= 0) THEN
          SPACE = INTREE + INPOOL + LOST
      ELSE
          SPACE = ISTORE(OFFSET+1-HOLES)
          HOLES = HOLES - 1
      END IF
!
! Compute sortkey
!
      SORKEY = MAXVAL(ERREST(1:NRFUNC))
!
! Put the new record in the heap
!
      POINT = (SPACE-1)*BLOCK
      IF (NRFUNC > 1) THEN
          POINT = POINT + 1
          RSTORE(POINT) = SORKEY
      END IF
      RSTORE(POINT+1:POINT+NRFUNC) = ERREST(1:NRFUNC)
      RSTORE(POINT+NRFUNC+1 : POINT+NRFUNC*2) = INTAPP(1:NRFUNC)
      POINT = POINT + NRFUNC*2
      DO I = 1,DIMENS
         RSTORE(POINT+1:POINT+NRVERT) = VERTIC(I,1:NRVERT)
         POINT = POINT + NRVERT
      END DO
      RSTORE(POINT+1:POINT+NRINFO) = RRGINF(1:NRINFO)

      POINT = BOTTIS - NIINFO*SPACE
      ISTORE(POINT+1:POINT+NIINFO) = IRGINF(1:NIINFO)
!
! Insert the index in the tree
!
      SUBRGN = INTREE
      DO
         SUBTMP = SUBRGN/2
         IF (SUBTMP >= 1) THEN
!
!          Compare max. child with parent.
!          If parent is max, then done.
!
             IF (SORKEY > RSTORE(1+ (ISTORE(OFFSET+SUBTMP)-1)*BLOCK)) THEN
!
!                Move the pointer at position subtmp down the heap.
!
                 ISTORE(OFFSET+SUBRGN) = ISTORE(OFFSET+SUBTMP)
                 SUBRGN = SUBTMP
                 CYCLE
             END IF
         END IF
         EXIT
      END DO
!
!  Set the pointer to the new index in the heap.
!
      ISTORE(OFFSET+SUBRGN) = SPACE
!
!  Save data structure variables
!
      ISTORE(12) = INTREE
      ISTORE(14) = HOLES
      IFAIL = 0
!
!***END DSSPUT
!
      RETURN
      END SUBROUTINE DSSPUT
!-----------------------------------------------------------------------
      SUBROUTINE DSGET(VERTIC,INTAPP,ERREST,IRGINF,RRGINF,ISTORE,RSTORE,&
                       IFAIL)
!***BEGIN PROLOGUE DSGET
!***DATE WRITTEN   900612   (YYMMDD)
!***REVISION DATE  961206   (YYMMDD)
!***AUTHOR
!          Ronald Cools, Dept. of Computer Science,
!          Katholieke Universiteit Leuven, Celestijnenlaan 200A,
!          B-3001 Heverlee, Belgium
!          Email:  Ronald.Cools@cs.kuleuven.ac.be
!
!***PURPOSE  To get and delete the record at the root of a sorted tree
!            that stores information about the subregions in an
!            adaptive integrator.
!***DESCRIPTION
!   See the description of SUBROUTINE DSINIT
!
!  Output parameters
!  -----------------
!
!  VERTIC = double precision array of dimension (DIMENS,NRVERT)
!           Contains the vertices that describe a subregion.
!           ( VERTIC(1,i),...,VERTIC(DIMENS,i) ) are the coordinates
!           of the i-th vertex.
!  INTAPP = double precision array of dimension (NRFUNC)
!           Contains approximations to the integrals.
!  ERREST = double precision array of dimension (NRFUNC)
!           Contains error estimates.
!  IRGINF = integer array of dimension (NIINFO).
!           Contains additional information about the subregion
!  RRGINF = double precision array of dimension (NRINFO).
!           Contains additional information about the subregion
!  ISTORE = integer array of dimension (BOTTIS).
!           Contains the integers of the records.
!  RSTORE = double precision array of dimension (BOTTRS).
!           Contains the double precision numbers of the records.
!  IFAIL = integer to indicate success or failure
!          IFAIL = 0 for normal exit
!          IFAIL = 10011 : attempt to get something out of empty tree
!          IFAIL = 9999 if a hole is created that cannot be saved.
!                     This is a warning, not an error !
!
!***END PROLOGUE DSGET
!
! Variables and constants used by DS-procedures only
! --------------------------------------------------
!
!  DIMENS = dimension of (sub-)regions
!  NRVERT = number of vertices to describe a subregion
!  NIINFO = number of integers used to save additional information
!           for each subregion
!  NRINFO = number of double precision numbers used to save additional
!           information for each subregion
!  NRFUNC = number of integrand functions for which information must
!           be stored
!  NRVACA = maximum number of pointers to empty spaces in the heap
!           that must be stored for later use
!  BOTTIS = length of integer array ISTORE.
!           Needed because this heap is also filled up starting from
!           the bottom
!  INTREE = number of subregions in the sorted tree
!  HOLES  = number of pointers to holes in the heap
!           0 <= HOLES <= NRVACA
!  LOST =   number of records that cannot be accessed any more
!  OFFSET : The first pointer to the tree is stored in ISTORE(OFFSET+1)
!  BLOCK =  number of double precision numbers in a record
!
!  Global variables
!
      INTEGER, INTENT(OUT) :: IFAIL
      INTEGER, DIMENSION(:), INTENT(IN OUT) :: ISTORE
      INTEGER, DIMENSION(:), INTENT(OUT) :: IRGINF
      REAL(kind=stnd), DIMENSION(:), INTENT(OUT) :: INTAPP,ERREST,RRGINF
      REAL(kind=stnd), DIMENSION(:), INTENT(IN OUT) :: RSTORE
      REAL(kind=stnd), DIMENSION(:,:), INTENT(OUT) :: VERTIC
!
! Local variables
!
! SORKEY = sortkey for maintaining the partially sorted tree
!        = maximum(errest(1),...,errest(nrfunc)
! SPACE  = index to the place were a new record can be
!          inserted in the heap
!
      INTEGER :: DIMENS,NRVERT,NIINFO,NRINFO,NRFUNC,NRVACA,BOTTIS,INTREE,  &
                 INPOOL,HOLES,LOST,BLOCK,OFFSET
      REAL(kind=stnd):: SORKEY
      INTEGER :: SUBRGN,SUBTMP,I,POINT,SPACE
!***FIRST EXECUTABLE STATEMENT
!  Initialise data structure constants
      DIMENS = ISTORE(1)
      NRVERT = ISTORE(2)
      NIINFO = ISTORE(3)
      NRINFO = ISTORE(4)
      NRFUNC = ISTORE(5)
      NRVACA = ISTORE(6)
      BOTTIS = ISTORE(7)
!     bottrs = istore(8)
      OFFSET = ISTORE(9)
!     start = istore(10)
      BLOCK = ISTORE(11)
      INTREE = ISTORE(12)
      INPOOL = ISTORE(13)
      HOLES = ISTORE(14)
      LOST = ISTORE(15)
!
! Check if something is in the sorted tree
!
      IF (INTREE <= 0) THEN
          IFAIL = 10011

      ELSE
!
! Get the top-record out of the tree
!
          SPACE = ISTORE(OFFSET+1)
          IF (NRFUNC == 1) THEN
              POINT = (SPACE-1)*BLOCK
          ELSE
              POINT = (SPACE-1)*BLOCK + 1
          END IF

          ERREST(1:NRFUNC) = RSTORE(POINT+1 : POINT+NRFUNC)
          INTAPP(1:NRFUNC) = RSTORE(POINT+NRFUNC+1 : POINT+NRFUNC*2)
          POINT = POINT + 2*NRFUNC
          DO I = 1,DIMENS
             VERTIC(I,1:NRVERT) = RSTORE(POINT+1:POINT+NRVERT)
             POINT = POINT + NRVERT
          END DO
          RRGINF(1:NRINFO) = RSTORE(POINT+1:POINT+NRINFO)

          POINT = BOTTIS - NIINFO*SPACE
          IRGINF(1:NIINFO) = ISTORE(POINT+1:POINT+NIINFO)
!
          IF (SPACE /= INTREE+HOLES+LOST+INPOOL) THEN
!
!             Fill the record with zeros.
!
                  POINT = (SPACE-1)*BLOCK
                  RSTORE(POINT + 1 : POINT + BLOCK) = 0
                  POINT = BOTTIS - NIINFO*SPACE
                  ISTORE(POINT + 1 : POINT + NIINFO) = 0
!
              IF (HOLES == NRVACA) THEN
                  LOST = LOST + 1

              ELSE
!
!             Save the place of the hole in the heap
!
                  ISTORE(OFFSET-HOLES) = SPACE
                  HOLES = HOLES + 1
              END IF

          END IF
!
! Rearrange the tree
!
          SORKEY = RSTORE(1+BLOCK* (ISTORE(OFFSET+INTREE)-1))
          INTREE = INTREE - 1
          SUBRGN = 1
          DO
          SUBTMP = 2*SUBRGN
          IF (SUBTMP <= INTREE) THEN
              IF (SUBTMP /= SUBRGN) THEN
!
!          Find max .of left and right child
!
                  IF (RSTORE(1+BLOCK* (ISTORE(OFFSET+SUBTMP)-1)) <   &
                      RSTORE(1+BLOCK* (ISTORE(OFFSET+SUBTMP+1)- 1))) THEN
                          SUBTMP = SUBTMP + 1
                  END IF
              END IF
!
!        Compare max .child with parent
!        If parent is max., then done
!
              IF (SORKEY < RSTORE(1+BLOCK* (ISTORE(OFFSET+SUBTMP)- 1)))&
                 THEN
!
!           Move the pointer at position subtmp up the heap.
!
                  ISTORE(OFFSET+SUBRGN) = ISTORE(OFFSET+SUBTMP)
                  SUBRGN = SUBTMP
                  CYCLE

              END IF

          END IF
          EXIT
          END DO
!
!  Update the pointer
!
          IF (INTREE > 0) THEN
             ISTORE(OFFSET+SUBRGN) = ISTORE(OFFSET+INTREE+1)
          END IF
!
!  Save data structure variables
!
          ISTORE(12) = INTREE
          ISTORE(14) = HOLES
          IF ( ISTORE(15) /= LOST ) THEN
             ISTORE(15) = LOST
             IFAIL = 9999
          ELSE
             IFAIL = 0
          END IF
!
!***END DSGET
!
      END IF
      RETURN

      END  SUBROUTINE DSGET
!-----------------------------------------------------------------------
      FUNCTION DSFREE(ISTORE) RESULT(FREE)    ! INTEGER RESULT
!***BEGIN PROLOGUE DSFREE
!***DATE WRITTEN   900612   (YYMMDD)
!***REVISION DATE  950829   (YYMMDD)
!***AUTHOR
!          Ronald Cools, Dept. of Computer Science,
!          Katholieke Universiteit Leuven, Celestijnenlaan 200A,
!          B-3001 Heverlee, Belgium
!          Email:  Ronald.Cools@cs.kuleuven.ac.be
!
!***PURPOSE  To return the number of records that may be added to the
!            data structure before all space is used.
!***DESCRIPTION
!   See the description of SUBROUTINE DSINIT
!
!***END PROLOGUE DSFREE
!
! Variables and constants used by DS-procedures only
! --------------------------------------------------
!
!  NIINFO = number of integers used to save additional information
!           for each subregion
!  BOTTIS = length of integer array ISTORE.
!           Needed because this heap is also filled up starting from
!           the bottom
!  BOTTRS = length of the double precision heap
!           Only needed for checks
!  INTREE = number of subregions in the sorted tree
!  INPOOL = number of subregions in the unsorted pool
!  LOST = number of records that cannot be accessed any more
!  OFFSET : The first pointer to the tree is stored in ISTORE(OFFSET+1)
!  BLOCK = number of double precision numbers in a record
!
!  Global variables
!
      INTEGER, DIMENSION(:), INTENT(IN) :: ISTORE
!
! Local variables
!
!   freerh = number of free records in the double precision heap
!   freeih = number of free records in the integer heap
!
      INTEGER :: FREE
      INTEGER :: NIINFO,BOTTIS,BOTTRS,INTREE,INPOOL,LOST,BLOCK,OFFSET
      INTEGER :: FREERS,FREEIS
!***FIRST EXECUTABLE STATEMENT
!  Initialise data structure constants
      NIINFO = ISTORE(3)
      BOTTIS = ISTORE(7)
      BOTTRS = ISTORE(8)
      OFFSET = ISTORE(9)
      BLOCK = ISTORE(11)
!  Initialise data structure variables
      INTREE = ISTORE(12)
      INPOOL = ISTORE(13)
      LOST = ISTORE(15)
!
      FREERS = BOTTRS/BLOCK - INTREE - INPOOL - LOST
      FREEIS = (BOTTIS-OFFSET)/ (NIINFO+1) - INTREE - INPOOL - LOST
      FREE = MIN(FREERS,FREEIS)
!
!***END DSFREE
!
      RETURN
      END FUNCTION DSFREE
!-----------------------------------------------------------------------
      FUNCTION DSUSED(ISTORE) RESULT(USED)   ! INTEGER RESULT
!***BEGIN PROLOGUE DSUSED
!***DATE WRITTEN   900612   (YYMMDD)
!***REVISION DATE  950829   (YYMMDD)
!***AUTHOR
!          Ronald Cools, Dept. of Computer Science,
!          Katholieke Universiteit Leuven, Celestijnenlaan 200A,
!          B-3001 Heverlee, Belgium
!          Email:  Ronald.Cools@cs.kuleuven.ac.be
!
!***PURPOSE  To return the number of records in the data structure.
!
!***DESCRIPTION
!   See the description of SUBROUTINE DSINIT
!
!***END PROLOGUE DSUSED
!
! Variables and constants used by DS-procedures only
! --------------------------------------------------
!
!  INTREE = number of subregions in the sorted tree
!  INPOOL = number of subregions in the unsorted pool
!
!  Global variables
!
      INTEGER, DIMENSION(:), INTENT(IN) :: ISTORE
!
!  Local variables
!
      INTEGER :: USED
      INTEGER :: INTREE,INPOOL
!***FIRST EXECUTABLE STATEMENT
!  Initialise data structure variables
      INTREE = ISTORE(12)
      INPOOL = ISTORE(13)
!
      USED = INTREE + INPOOL
!
!***END DSUSED
!
      RETURN
      END FUNCTION DSUSED
!-----------------------------------------------------------------------
      SUBROUTINE DSSUM(VALUE,ABSERR,ISTORE,RSTORE,IFAIL)
!***BEGIN PROLOGUE DSSUM
!***DATE WRITTEN   901129   (YYMMDD)
!***REVISION DATE  950829   (YYMMDD)
!***AUTHOR
!          Ronald Cools, Dept. of Computer Science,
!          Katholieke Universiteit Leuven, Celestijnenlaan 200A,
!          B-3001 Heverlee, Belgium
!          Email:  Ronald.Cools@cs.kuleuven.ac.be
!
!***PURPOSE To compute more accurate values of VALUE and ABSERR.
!   This is done to reduce the effect of roundoff on final results.
!   Large intermediate sums in the computation may course large,
!   unnecessary roundoff errors. Thus recomputing the sums of errors
!   and estimates and in addition grouping the sums in three groups
!   should remove this problem.
!
!***DESCRIPTION
!   See the description of SUBROUTINE DSINIT
!
!  Input parameters
!  ----------------
!
!  ISTORE = integer array of dimension (BOTTIS).
!           Contains the integers of the records.
!  RSTORE = double precision array of dimension (BOTTRS).
!           Contains the double precision numbers of the records.
!
!  Output parameters
!  -----------------
!
!  VALUE = double precision array of dimension (NRFUNC)
!           Contains approximations to the integrals.
!  ABSERR = double precision array of dimension (NRFUNC)
!           Contains error estimates.
!  IFAIL  = integer to indicate success or failure
!           IFAIL = 0 for normal termination.
!old        IFAIL = 10012 if HOLES is not equal to zero.
!
!***END PROLOGUE DSSUM
!
! Variables and constants used by DS-procedures only
! --------------------------------------------------
!
!  DIMENS = dimension of (sub-)regions
!  NRVERT = number of vertices to describe a subregion
!  NIINFO = number of integers used to save additional information
!           for each subregion
!  NRINFO = number of double precision numbers used to save additional
!           information for each subregion
!  NRFUNC = number of integrand functions for which information must
!           be stored
!  NRVACA = maximum number of pointers to empty spaces in the heap
!           that must be stored for later use
!  BOTTIS = length of integer array ISTORE.
!           Needed because this heap is also filled up starting from
!           the bottom
!  BOTTRS = length of the double precision heap
!           Only needed for checks
!  INTREE = number of subregions in the sorted tree
!  INPOOL = number of subregions in the unsorted pool
!  HOLES  = number of pointers to holes in the heap
!           0 <= HOLES <= NRVACA
!  LOST = number of records that cannot be accessed any more
!  OFFSET : The first pointer to the tree is stored in ISTORE(OFFSET+1)
!  START : The first pointer to the pool is stored in ISTORE(START-1)
!  BLOCK = number of double precision numbers in a record
!
!  Global variables
!
      INTEGER, INTENT(OUT) :: IFAIL
      INTEGER, DIMENSION(:), INTENT(IN) :: ISTORE
      REAL(kind=stnd), DIMENSION(:), INTENT(IN) :: RSTORE
      REAL(kind=stnd), DIMENSION(:), INTENT(OUT) :: ABSERR,VALUE
!
! Local variables
!
!
      INTEGER :: NRFUNC,INTREE,INPOOL,HOLES,LOST,BLOCK
      REAL(kind=stnd), DIMENSION(2) :: SUMVAL, SUMERR
      INTEGER :: SBRGNS,I,J,POINT,I1,I2,GROUP,LEFT
!***FIRST EXECUTABLE STATEMENT
!  Initialise data structure constants
      NRFUNC = ISTORE(5)
      BLOCK = ISTORE(11)
!  Initialise data structure variables
      INTREE = ISTORE(12)
      INPOOL = ISTORE(13)
      HOLES = ISTORE(14)
      LOST = ISTORE(15)
!
! Assume there are no holes
!
!     IF (HOLES /= 0) THEN
!         IFAIL = 10012
!         RETURN
!     END IF
      SBRGNS = INTREE + INPOOL + LOST + HOLES
      DO J = 1,NRFUNC
          IF (NRFUNC == 1) THEN 
                                  POINT = 1
                           ELSE 
                                  POINT = J + 1
          END IF
          VALUE(J) = 0
          ABSERR(J) = 0
          GROUP = SBRGNS**(1.0/3.0)     ! This is accurate enough
          IF (GROUP**3 /= SBRGNS) THEN
              GROUP = GROUP + 1
          END IF
          LEFT = SBRGNS
          DO I2 = 0,SBRGNS - 1,GROUP**2
              SUMVAL(2) = 0
              SUMERR(2) = 0
              DO I1 = 0,MIN(LEFT,GROUP**2) - 1,GROUP
                  SUMVAL(1) = 0
                  SUMERR(1) = 0
                  DO I = 1 + I1 + I2,MIN(LEFT,GROUP) + I1 + I2
                      SUMVAL(1) = SUMVAL(1) + RSTORE(POINT+NRFUNC)
                      SUMERR(1) = SUMERR(1) + RSTORE(POINT)
                      POINT = POINT + BLOCK
                  END DO
                  LEFT = LEFT - MIN(LEFT,GROUP)
                  SUMVAL(2) = SUMVAL(2) + SUMVAL(1)
                  SUMERR(2) = SUMERR(2) + SUMERR(1)
              END DO
              VALUE(J) = VALUE(J) + SUMVAL(2)
              ABSERR(J) = ABSERR(J) + SUMERR(2)
          END DO
      END DO
      IFAIL = 0
!
!***END DSSUM
!
      RETURN
      END SUBROUTINE DSSUM
!-----------------------------------------------------------------------
      SUBROUTINE DSSTAT(ISTORE,RSTORE)
!***DATE WRITTEN   901129   (YYMMDD)
!***REVISION DATE  970320   (YYMMDD)
!
!  Global variables
!
      INTEGER, DIMENSION(:), INTENT(IN) :: ISTORE
      REAL(kind=stnd), DIMENSION(:), INTENT(IN) :: RSTORE
!
!  Global functions
!
!     INTEGER :: DSFREE,DSUSED
!
!  Local variables

      INTEGER :: I,J,POINT,NIINFO,START
      INTEGER, DIMENSION(15) :: HIST
!
      write(unit=*,fmt=*) "DATA STRUCTURE CONSTANTS:"
      write(unit=*,fmt=*) "DIMENS = ",ISTORE(1)
      write(unit=*,fmt=*) "NRVERT = ",ISTORE(2)
      write(unit=*,fmt=*) "NRFUNC = ",ISTORE(5)
      write(unit=*,fmt=*) "NIINFO = ",ISTORE(3)
      write(unit=*,fmt=*) "NRINFO = ",ISTORE(4)
      write(unit=*,fmt=*) "NRVACA = ",ISTORE(6)
      write(unit=*,fmt=*) "BLOCK  = ",ISTORE(11)
      write(unit=*,fmt=*) "BOTTIS = ",ISTORE(7)
      write(unit=*,fmt=*) "BOTTRS = ",ISTORE(8)
      write(unit=*,fmt=*) "OFFSET = ",ISTORE(9)
      write(unit=*,fmt=*) "START  = ",ISTORE(10)
!
      write(unit=*,fmt=*) "DATA STRUCTURE VARIABLES:"
      write(unit=*,fmt=*) "INTREE = ",ISTORE(12)
      write(unit=*,fmt=*) "INPOOL = ",ISTORE(13)
      write(unit=*,fmt=*) "HOLES  = ",ISTORE(14)
      write(unit=*,fmt=*) "LOST   = ",ISTORE(15)
      write(unit=*,fmt=*) "FREE   = ",DSFREE(ISTORE)
!
      DO I = 1,15
          HIST(I) = 0
      END DO
      NIINFO = ISTORE(3)
      START = ISTORE(7) - niinfo + 2
      DO I = START,START - (DSUSED(ISTORE)-1)*NIINFO,-NIINFO
          IF (ISTORE(I) < 15) THEN
              HIST(ISTORE(I)) = HIST(ISTORE(I)) + 1
          ELSE
              HIST(15) = HIST(15) + 1
          END IF
      END DO
      write(unit=*,fmt=*) "SUBDIVISION INFORMATION:"
      DO I = 1,14
          write(unit=*,fmt=*) "RATIO = ",I," -> ",HIST(I)," REGIONS"
      END DO
      write(unit=*,fmt=*) "RATIO > 14 -> ",HIST(15)," REGIONS"
!
!     dump RSTORE to fort.1
!
      open (unit=1,status="replace",action="write",file="tmp_realstore")
      point = 1
      do i=1,istore(12)+istore(13)+istore(14)+istore(15)
          write (unit=1,fmt="( "" :"",i3,"":"",30(""-"") )") i
          do j = 1,istore(11)
           write(unit=1,fmt="(4x,i3,"" -> "",es24.16)") j,rstore(point)
             point = point+1
         end do
      end do
!
!     dump ISTORE to fort.2
!
       open (unit=2,status="replace",action="write",file="tmp_integerstore")
       point = istore(7)
       do  i=1,istore(12)+istore(13)+istore(14)+istore(15)
         write(unit=2,fmt="("" :"",i3,"":"",30(""-""))") i
           do  j = 1,istore(3)
             write(unit=2,fmt="(4x,i3,"" -> "",i10)") j,istore(point)
              point = point-1
          end do
       end do

       write(unit=2,fmt=*) "The pointers in the TREE:"
       do i = istore(9)+1,istore(9)+istore(12)
          write(unit=2,fmt=*) istore(i)
       end do
       write(unit=2,fmt=*) "The pointers in the POOL:"
       do i = istore(10)-1,istore(10)-istore(13),-1
           write(unit=2,fmt=*) istore(i)
       end do
       close(unit=1)
       write(unit=*,fmt=*) "Real    part of data structure is dumped in file tmp_realstore."
       write(unit=*,fmt=*) "Integer part of data structure is dumped in file tmp_integerstore."
       close(unit=2)
!
      RETURN
      END SUBROUTINE DSSTAT
!-----------------------------------------------------------------------
      SUBROUTINE DSCOPY(ISTOREIN,RSTOREIN,ISTOREOUT,RSTOREOUT)
!***BEGIN PROLOGUE DSINIT
!***DATE WRITTEN   950612   (YYMMDD)
!***REVISION DATE  990607   (YYMMDD)   (bug fix)
!***AUTHOR
!          Ronald Cools, Dept. of Computer Science,
!          Katholieke Universiteit Leuven, Celestijnenlaan 200A,
!          B-3001 Heverlee, Belgium
!          Email:  Ronald.Cools@cs.kuleuven.ac.be
!
!***PURPOSE  To copy one region collection into the other.
      INTEGER, DIMENSION(:), INTENT(IN)  :: ISTOREIN
      INTEGER, DIMENSION(:), INTENT(OUT) :: ISTOREOUT
      REAL(kind=stnd), DIMENSION(:), INTENT(IN)  :: RSTOREIN
      REAL(kind=stnd), DIMENSION(:), INTENT(OUT) :: RSTOREOUT

      INTEGER :: RECORDS, NIINFO, NRVACA, BOTTISIN, OFFSET, STARTIN, &
              BLOCK, INTREE, INPOOL, HOLES,LOST,BOTTISOUT,BOTTRSOUT,&
              STARTOUT

      INTEGER, PARAMETER :: CONST=15

! Constants of incoming data structure
      NIINFO = ISTOREIN(3)
      NRVACA = ISTOREIN(6)
      BOTTISIN = ISTOREIN(7)
      OFFSET = ISTOREIN(9)
      STARTIN = ISTOREIN(10)
      BLOCK = ISTOREIN(11)
      INTREE = ISTOREIN(12)
      INPOOL = ISTOREIN(13)
      HOLES = ISTOREIN(14)
      LOST = ISTOREIN(15)
      RECORDS = INTREE + INPOOL + HOLES + LOST
! Constants of outcoming data structure
      BOTTISOUT = size(ISTOREOUT)
      BOTTRSOUT = size(RSTOREOUT)
      STARTOUT = BOTTISOUT + 1 - ((BOTTISOUT-NRVACA-CONST)/ (1+NIINFO))*NIINFO

! Copy integer part of data structure
      ISTOREOUT(1:CONST) = ISTOREIN(1:CONST)
      ISTOREOUT(7) = BOTTISOUT
      ISTOREOUT(8) = BOTTRSOUT
      ISTOREOUT(10) = STARTOUT
      ISTOREOUT(OFFSET-HOLES+1:OFFSET) = ISTOREIN(OFFSET-HOLES+1:OFFSET)
      ISTOREOUT(OFFSET+1:OFFSET+INTREE) = ISTOREIN(OFFSET+1:OFFSET+INTREE)
      ISTOREOUT(STARTOUT-INPOOL:STARTOUT-1) = ISTOREIN(STARTIN-INPOOL:STARTIN-1)
      ISTOREOUT(BOTTISOUT-records*niinfo:BOTTISOUT) = ISTOREIN(BOTTISIN-records*niinfo:BOTTISIN)

! Copy real part of data structure
      RSTOREOUT(1:BLOCK*records) = rstorein(1:BLOCK*records)

      return 
      end subroutine dscopy
!-----------------------------------------------------------------------

SUBROUTINE DSUPUT(VERTIC,INTAPP,ERREST,IRGINF,RRGINF,ISTORE,RSTORE,IFAIL)                                   
!***BEGIN PROLOGUE DSUPUT                                               
!***DATE WRITTEN   901126   (YYMMDD)                                    
!***REVISION DATE  961206   (YYMMDD)                                    
!***AUTHOR                                                              
!          Ronald Cools, Dept. of Computer Science,                     
!          Katholieke Universiteit Leuven, Celestijnenlaan 200A,        
!          B-3001 Heverlee, Belgium                                     
!          Email:  Ronald.Cools@cs.kuleuven.ac.be                            
!                                                                       
!***PURPOSE  To add a record to a data structure that stores information
!            about the subregions in an adaptive integrator.            
!            The record is added to an unsorted pool.                   
!***DESCRIPTION                                                         
!   See the description of SUBROUTINE DSINIT                            
!                                                                       
!  Input parameters                                                     
!  ----------------                                                     
!                                                                       
!  VERTIC = double precision array of dimension (DIMENS,NRVERT)         
!           Contains the vertices that describe a subregion.            
!           ( VERTIC(1,i),...,VERTIC(DIMENS,i) ) are the coordinates    
!           of the i-th vertex.                                         
!  INTAPP = double precision array of dimension (NRFUNC)                
!           Contains approximations to the integrals.                   
!  ERREST = double precision array of dimension (NRFUNC)                
!           Contains error estimates.                                   
!  IRGINF = integer array of dimension (NIINFO).                        
!           Contains additional information about the subregion         
!  RRGINF = double precision array of dimension (NRINFO).               
!           Contains additional information about the subregion         
!  ISTORE = integer array of dimension (BOTTIS).                        
!           Contains the integers of the records.                       
!  RSTORE = double precision array of dimension (BOTTRS).               
!           Contains the double precision numbers of the records.       
!                                                                       
!  Output parameters                                                    
!  -----------------                                                    
!                                                                       
!  IFAIL = integer to indicate success or failure                       
!          IFAIL = 0 for normal exit                                    
!          IFAIL = 10009 if the integer array ISTORE is full                
!          IFAIL = 10010 if the double precision array RSTORE is full      
!  ISTORE = integer array of dimension (BOTTIS).                        
!           Contains the integers of the records.                       
!  RSTORE = double precision array of dimension (BOTTRS).               
!           Contains the double precision numbers of the records.       
!                                                                       
!***END PROLOGUE DSUPUT                                                 
!                                                                       
! Variables and constants used by DS-procedures only                    
! --------------------------------------------------                    
!                                                                       
!  DIMENS = dimension of (sub-)regions                                  
!  NRVERT = number of vertices to describe a subregion                  
!  NIINFO = number of integers used to save additional information      
!           for each subregion                                          
!  NRINFO = number of double precision numbers used to save additional  
!           information for each subregion                              
!  NRFUNC = number of integrand functions for which information must    
!           be stored                                                   
!  NRVACA = maximum number of pointers to empty spaces in the heap      
!           that must be stored for later use                           
!  BOTTIS = length of integer array ISTORE.                             
!           Needed because this heap is also filled up starting from    
!           the bottom                                                  
!  BOTTRS = length of the double precision heap                         
!           Only needed for checks                                      
!  INTREE = number of subregions in the sorted tree                     
!  INPOOL = number of subregions in the unsorted pool                   
!  HOLES  = number of pointers to holes in the heap                     
!           0 <= HOLES <= NRVACA                                        
!  LOST = number of records that cannot be accessed any more            
!  OFFSET : The first pointer to the tree is stored in ISTORE(OFFSET+1) 
!  START : The first pointer to the pool is stored in ISTORE(START-1)   
!  BLOCK = number of double precision numbers in a record               
!                                                                       
!  Global variables                                                     
!                                                                       
      INTEGER, INTENT(OUT) :: IFAIL 
      INTEGER, DIMENSION(:), INTENT(IN OUT) :: ISTORE
      INTEGER, DIMENSION(:), INTENT(IN) :: IRGINF
      REAL(kind=stnd), DIMENSION(:), INTENT(IN) :: INTAPP,ERREST,RRGINF
      REAL(kind=stnd), DIMENSION(:), INTENT(IN OUT) :: RSTORE 
      REAL(kind=stnd), DIMENSION(:,:), INTENT(IN) :: VERTIC
!                                                                       
! Local variables                                                       
!                                                                       
! SORKEY = sortkey for maintaining the partially sorted tree            
!        = maximum(errest(1),...,errest(nrfunc)                         
! SPACE  = index to the place were a new record can be                  
!          inserted in the heap                                         
!                                                                       
      INTEGER :: DIMENS,NRVERT,NIINFO,NRINFO,NRFUNC,BOTTIS,BOTTRS,INTREE,  &
              INPOOL,HOLES,LOST,BLOCK,OFFSET,START                      
      REAL(kind=stnd):: SORKEY 
      INTEGER :: I,SPACE,POINT 
!***FIRST EXECUTABLE STATEMENT                                          
!                                                                       
!  The following initialisation statements are included                 
!  to make the code readable                                            
!                                                                       
!  Initialise data structure constants                                  
      DIMENS = ISTORE(1) 
      NRVERT = ISTORE(2) 
      NIINFO = ISTORE(3) 
      NRINFO = ISTORE(4) 
      NRFUNC = ISTORE(5) 
      BOTTIS = ISTORE(7) 
      BOTTRS = ISTORE(8) 
      OFFSET = ISTORE(9) 
      START = ISTORE(10) 
      BLOCK = ISTORE(11) 
!  Initialise data structure variables                                  
      INTREE = ISTORE(12) 
      INPOOL = ISTORE(13) 
      HOLES = ISTORE(14) 
      LOST = ISTORE(15) 
!                                                                       
! Check if enough space is left in the arrays to put in a subregion     
!                                                                       
      IF (HOLES <= 0) THEN 
         IF ((BOTTIS-START+1- (INTREE+INPOOL+LOST)*NIINFO <           &
             NIINFO) .OR. (START-OFFSET-1-INTREE-INPOOL < 1)) THEN    
             IFAIL = 10009 
          ELSE IF (BOTTRS- (INTREE+INPOOL+LOST)*BLOCK < BLOCK) THEN 
              IFAIL = 10010 
              RETURN 
          END IF 
      END IF 
      INPOOL = INPOOL + 1 
!                                                                       
! Determine index for new record                                        
!                                                                       
      IF (HOLES <= 0) THEN 
          SPACE = INTREE + INPOOL + LOST 
      ELSE 
          SPACE = ISTORE(OFFSET+1-HOLES) 
          HOLES = HOLES - 1 
      END IF 
!                                                                       
! Compute sortkey (in case this record has to be sorted later)          
!                                                                       
      SORKEY = MAXVAL(ERREST(1:NRFUNC))
!                                                                       
! Put the new record in the heap                                        
!                                                                       
      POINT = (SPACE-1)*BLOCK 
      IF (NRFUNC > 1) THEN 
          POINT = POINT + 1 
          RSTORE(POINT) = SORKEY 
      END IF 
      RSTORE(POINT+1:POINT+NRFUNC) = ERREST(1:NRFUNC)
      RSTORE(POINT+NRFUNC+1 : POINT+NRFUNC*2) = INTAPP(1:NRFUNC)
      POINT = POINT + NRFUNC * 2
      DO I = 1,DIMENS
         RSTORE(POINT+1:POINT+NRVERT) = VERTIC(I,1:NRVERT)
         POINT = POINT + NRVERT
      END DO
      RSTORE(POINT+1:POINT+NRINFO) = RRGINF(1:NRINFO)

      POINT = BOTTIS - NIINFO*SPACE 
      ISTORE(POINT+1:POINT+NIINFO) = IRGINF(1:NIINFO)
!                                                                       
! Insert the index in the pool                                          
!                                                                       
      ISTORE(START-INPOOL) = SPACE 
!                                                                       
!  Save data structure variables                                        
!                                                                       
      ISTORE(13) = INPOOL 
      ISTORE(14) = HOLES 
      IFAIL = 0 
!                                                                       
!***END DSUPUT                                                          
!                                                                       
      RETURN 
      END SUBROUTINE DSUPUT               


!---------------------------------------------------------------------------
      SUBROUTINE DSPINT(ISTORE,RSTORE) 
!***BEGIN PROLOGUE DSPINT                                               
!***DATE WRITTEN   901129   (YYMMDD)                                    
!***REVISION DATE  961206   (YYMMDD)                                    
!***AUTHOR                                                              
!          Ronald Cools, Dept. of Computer Science,                     
!          Katholieke Universiteit Leuven, Celestijnenlaan 200A,        
!          B-3001 Heverlee, Belgium                                     
!          Email:  Ronald.Cools@cs.kuleuven.ac.be                            
!                                                                       
!***PURPOSE  All records in the unsorted part of the data structure     
!          are put in the sorted part of the data structure.            
!                                                                       
!***DESCRIPTION                                                         
!   See the description of SUBROUTINE DSINIT                            
!                                                                       
!  Input parameters                                                     
!  ----------------                                                     
!                                                                       
!  ISTORE = integer array of dimension (BOTTIS).                        
!           Contains the integers of the records.                       
!  RSTORE = double precision array of dimension (BOTTRS).               
!           Contains the double precision numbers of the records.       
!                                                                       
!  Output parameters                                                    
!  -----------------                                                    
!                                                                       
!  ISTORE = integer array of dimension (BOTTIS).                        
!           Contains the integers of the records.                       
!  RSTORE = double precision array of dimension (BOTTRS).               
!           Contains the double precision numbers of the records.       
!                                                                       
!***END PROLOGUE DSPINT                                                 
!                                                                       
! Variables and constants used by DS-procedures only                    
! --------------------------------------------------                    
!                                                                       
!  INTREE = number of subregions in the sorted tree                     
!  INPOOL = number of subregions in the unsorted pool                   
!  OFFSET : The first pointer to the tree is stored in ISTORE(OFFSET+1) 
!  START :  The first pointer to the pool is stored in ISTORE(START-1)  
!  BLOCK =  number of double precision numbers in a record              
!                                                                       
!  Global variables                                                     
!                                                                       
      INTEGER, DIMENSION(:), INTENT(IN OUT) :: ISTORE
      REAL(kind=stnd), DIMENSION(:), INTENT(IN OUT) :: RSTORE
!                                                                       
! Local variables                                                       
!                                                                       
! SORKEY = sortkey for maintaining the partially sorted tree            
!        = maximum(errest(1),...,errest(nrfunc)                         
! SPACE  = index to the place were a new record can be                  
!          inserted in the heap                                         
!                                                                       
      INTEGER :: INTREE,INPOOL,BLOCK,OFFSET,START
      REAL(kind=stnd):: SORKEY 
      INTEGER :: SPACE,SUBRGN,SUBTMP,NR 
!***FIRST EXECUTABLE STATEMENT                                          
!                                                                       
!  The following initialisation statements are included                 
!  to make the code readable                                            
!                                                                       
!  Initialise data structure constants                                  
      OFFSET = ISTORE(9) 
      START = ISTORE(10) 
      BLOCK = ISTORE(11) 
!  Initialise data structure variables                                  
      INTREE = ISTORE(12) 
      INPOOL = ISTORE(13) 
!                                                                       
!                                                                       
!                                                                       
      DO NR = INPOOL,1,-1 
          INTREE = INTREE + 1 
!                                                                       
! extract sortkey                                                       
!                                                                       
          SPACE = ISTORE(START-NR) 
          SORKEY = RSTORE((SPACE-1)*BLOCK+1) 
!                                                                       
! Insert the index in the tree                                          
!                                                                       
          SUBRGN = INTREE 
          DO
             SUBTMP = SUBRGN/2 
             IF (SUBTMP >= 1) THEN 
!                                                                       
!          Compare max. child with parent.                                 
!          If parent is max, then done.                                    
!                                                                       
                 IF (SORKEY > RSTORE(1+ (ISTORE(OFFSET+SUBTMP)-1)*BLOCK)) THEN                                       
!                                                                       
!                Move the pointer at position subtmp down the heap.        
!                                                                       
                     ISTORE(OFFSET+SUBRGN) = ISTORE(OFFSET+SUBTMP) 
                     SUBRGN = SUBTMP 
                     CYCLE
                 END IF 
             END IF 
             EXIT
          END DO
!                                                                       
!  Set the pointer to the new index in the heap.                        
!                                                                       
          ISTORE(OFFSET+SUBRGN) = SPACE 
      END DO 
      INPOOL = 0 
!                                                                       
!  Save data structure variables                                        
!                                                                       
      ISTORE(12) = INTREE 
      ISTORE(13) = INPOOL 
!                                                                       
!***END DSPINT                                                          
!                                                                       
      RETURN 
      END SUBROUTINE DSPINT

!----------------------------------------------------------------------------

END MODULE DS_ROUTINES
