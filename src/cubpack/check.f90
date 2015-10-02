! This file is F-compatible, except for upper/lower case conventions.
!--------------------------------------------------------------------
Module Check_Input

USE Precision_Model, ONLY: stnd
USE internal_types

Implicit NONE

PRIVATE

PUBLIC :: CHECK
PRIVATE :: CHECK_RESTART, CHECK_NORESTART

INTERFACE CHECK
   MODULE PROCEDURE CHECK_RESTART, CHECK_NORESTART
END INTERFACE

CONTAINS
      SUBROUTINE CHECK_RESTART(DIMENS,NUMFUN,BOTTRH,BOTTIH,&
                       ISTORE,INFORM,JOB,IFAIL,EPSABS,EPSREL, &
                       MINPTS,MAXPTS,TUNE)
!***BEGIN PROLOGUE CHECK_RESTART
!***REVISION DATE  950503   (YYMMDD)
!***REVISION DATE  980407   (YYMMDD)
!***REVISION DATE  990531   (YYMMDD)
!***REVISION DATE  990624   (YYMMDD)
!***REVISION DATE  010719   (YYMMDD)
!***AUTHOR
!          Ronald Cools, Dept. of Computer Science,
!          Katholieke Universiteit Leuven, Celestijnenlaan 200A,
!          B-3001 Heverlee, Belgium
!          Email:  Ronald.Cools@cs.kuleuven.ac.be

!***PURPOSE  CHECK_RESTART checks the validity of the
!            input parameters to CUBATR on restart.
!***DESCRIPTION
!
!
!   ON ENTRY
!
!     DIMENS Integer.
!            Number of variables.
!     NUMFUN Integer.
!            Number of components of the integral.
!     NUMRGN Integer.
!            NUMRGN should contain the initial number of regions.
!     MINPTS Integer.
!            Minimum number of FUNSUB calls.
!     MAXPTS Integer.
!            Maximum number of FUNSUB calls.
!     TUNE   Real.
!            Requested reliability.
!     JOB    Integer.
!            Describes what algorithm CUBATR must use.
!     EPSABS Real.
!            Requested absolute accuracy.
!     EPSREL Real.
!            Requested relative accuracy.
!     BOTTRH Integer.
!            Defines the length of the real array containing
!            the region collection.
!     BOTTIH Integer.
!            Defines the length of the integer array containing
!            the region collection.
!     RGTYPE Integer array of dimension (NUMRGN).
!           RGTYPE(L) describes the type of region L.
!
!   ON RETURN
!
!     IFAIL  Integer.
!            If IFAIL has an illegal value on entry, it is reset to 0.
!
!     INFORM Integer.
!            INFORM = 0 for normal exit.
!            INFORM = 8 if DIMENS <= 1                ( RESTART = FALSE )
!                          DIMENS used inconsistently ( RESTAR = TRUE ).
!            INFORM = 16 if NUMFUN < 1.               ( RESTART = FALSE )
!                          NUMFUN used inconsistently ( RESTART = TRUE ).
!            INFORM = 32 if NUMRGN < 1.               ( RESTART = FALSE )
!                          BOTTIH used inconsistently ( RESTART = TRUE ).
!            INFORM = 64 if there is no support for this region (RESTART = FALSE )
!                          BOTTRH used inconsistently           ( RESTART = TRUE ).
!            INFORM = 128 if ISTORE is strangely small (RESTART = TRUE)
!            INFORM = 256 if MAXPTS <1 or MINPTS > MAXPTS.
!            INFORM = 512 if TUNE > 1 OR TUNE < 0.
!            INFORM = 1024 if EPSABS < 0 or EPSREL < 0
!            INFORM = 2048 if IFAIL not in {-1,0,1}
!            INFORM = 4096 if ABS(JOB) not in {0,1,2,11,12}
!            If more than one input parameter is wrong,
!            then INFORM is set to the sum of the values mentioned above.
!            If errors occured during a restart, INFORM = INFORM + 8192.
!
!***END PROLOGUE CHECK_RESTART
!
!   Global variables.
!
      INTEGER, INTENT(IN)  :: DIMENS,NUMFUN,BOTTRH,BOTTIH,JOB
      INTEGER, INTENT(OUT) :: INFORM
      INTEGER, DIMENSION(:), INTENT(IN) :: ISTORE
      INTEGER, OPTIONAL, INTENT(IN OUT)  :: IFAIL
      INTEGER, OPTIONAL, INTENT(IN)     :: MINPTS,MAXPTS
      REAL(kind=stnd), OPTIONAL, INTENT(IN)  :: EPSABS,EPSREL,TUNE
!
!***FIRST EXECUTABLE STATEMENT CHECK_RESTART
!

      INFORM = 0
          IF ( SIZE(ISTORE) <= 15) THEN
              INFORM = 128
          ELSE
              ! We assume that istore is present and that BOTTIH and
              ! BORTRH have the right values.
              !
              ! Check valid DIMENS.
              !
              IF (DIMENS /= ISTORE(1)) THEN 
                                            INFORM = 8
              END IF
              !
              ! Check positive NUMFUN.
              !
              IF (NUMFUN /= ISTORE(5)) THEN
                                            INFORM = INFORM + 16
              END IF
              !
              ! Check workspace.
              !
              IF ((BOTTIH /= ISTORE(7)) .OR. (SIZE(ISTORE) /= BOTTIH)) THEN
                 INFORM = INFORM + 32
              END IF
              IF (BOTTRH /= ISTORE(8)) THEN
                                            INFORM = INFORM + 64
              END IF
          END IF
          IF (INFORM /= 0) THEN
                                INFORM = INFORM + 8192
          END IF
!
!     Check valid limits on allowed number of function evaluations.
!
      IF (PRESENT(MAXPTS)) THEN
         IF (PRESENT(MINPTS)) THEN
            IF (MINPTS > MAXPTS) THEN
                                      INFORM = INFORM + 256
            END IF
         ELSE
            IF (MAXPTS < 1) THEN 
                                 INFORM = INFORM + 256
                            END IF
         END IF
      END IF
!
!     Check valid requested reliablitiy.
!
      IF (PRESENT(TUNE)) THEN
         IF ((TUNE < 0) .OR. (TUNE > 1)) THEN
                                              INFORM = INFORM + 512
         END IF
      END IF
!
!     Check valid accuracy requests.
!
      IF (PRESENT(EPSABS)) THEN
         IF (EPSABS < 0) THEN
                              INFORM = INFORM + 1024
         END IF
      ELSE IF (PRESENT(EPSREL)) THEN
         IF (EPSREL < 0) THEN
                              INFORM = INFORM + 1024
         END IF
      END IF
!
!     Check valid IFAIL
!
      IF (PRESENT(IFAIL)) THEN
         IF ((IFAIL < -1) .OR. (IFAIL > 1)) THEN
              INFORM = INFORM + 2048
              IFAIL = 0
         END IF
      END IF
!
!     Check valid JOB. JOB = 0 cannot occur at this stage
!    
      IF ((ABS(JOB) /= 1) .AND. (ABS(JOB) /= 2) .AND. (ABS(JOB) /= 11) .AND.  (ABS(JOB) /= 12)) THEN
                                                      INFORM = INFORM + 4096
      END IF
!
      RETURN
      END SUBROUTINE CHECK_RESTART

      SUBROUTINE CHECK_NORESTART(DIMENS,NUMFUN,NUMRGN,RGTYPE, &
                       INFORM,JOB,IFAIL,EPSABS,EPSREL, &
                       MINPTS,MAXPTS,TUNE)
!***BEGIN PROLOGUE CHECK_NORESTART
!***REVISION DATE  950503   (YYMMDD)
!***REVISION DATE  980407   (YYMMDD)
!***REVISION DATE  990531   (YYMMDD)
!***AUTHOR
!          Ronald Cools, Dept. of Computer Science,
!          Katholieke Universiteit Leuven, Celestijnenlaan 200A,
!          B-3001 Heverlee, Belgium
!          Email:  Ronald.Cools@cs.kuleuven.ac.be

!***PURPOSE  CHECK_NORESTART checks the validity of the
!            input parameters to CUBATR.
!***DESCRIPTION
!
!
!   ON ENTRY
!
!     DIMENS Integer.
!            Number of variables.
!     NUMFUN Integer.
!            Number of components of the integral.
!     NUMRGN Integer.
!            NUMRGN should contain the initial number of regions.
!     MINPTS Integer.
!            Minimum number of FUNSUB calls.
!     MAXPTS Integer.
!            Maximum number of FUNSUB calls.
!     TUNE   Real.
!            Requested reliability.
!     JOB    Integer.
!            Describes what algorithm CUBATR must use.
!     EPSABS Real.
!            Requested absolute accuracy.
!     EPSREL Real.
!            Requested relative accuracy.
!     BOTTRH Integer.
!            Defines the length of the real array containing
!            the region collection.
!     BOTTIH Integer.
!            Defines the length of the integer array containing
!            the region collection.
!     RGTYPE Integer array of dimension (NUMRGN).
!           RGTYPE(L) describes the type of region L.
!
!   ON RETURN
!
!     IFAIL  Integer.
!            If IFAIL has an illegal value on entry, it is reset to 0.
!
!     INFORM Integer.
!            INFORM = 0 for normal exit.
!            INFORM = 8 if DIMENS <= 1                ( RESTART = FALSE )
!                          DIMENS used inconsistently ( RESTAR = TRUE ).
!            INFORM = 16 if NUMFUN < 1.               ( RESTART = FALSE )
!                          NUMFUN used inconsistently ( RESTART = TRUE ).
!            INFORM = 32 if NUMRGN < 1.               ( RESTART = FALSE )
!                          BOTTIH used inconsistently ( RESTART = TRUE ).
!            INFORM = 64 if there is no support for this region (RESTART = FALSE )
!                          BOTTRH used inconsistently           ( RESTART = TRUE ).
!            INFORM = 128 if JOB={2,11} and DIMENS > 3
!            INFORM = 256 if MAXPTS <1 or MINPTS > MAXPTS.
!            INFORM = 512 if TUNE > 1 OR TUNE < 0.
!            INFORM = 1024 if EPSABS < 0 or EPSREL < 0
!            INFORM = 2048 if IFAIL not in {-1,0,1}
!            INFORM = 4096 if ABS(JOB) not in {0,1,2,11,12}
!            If more than one input parameter is wrong,
!            then INFORM is set to the sum of the values mentioned above.
!            If errors occured during a restart, INFORM = INFORM + 8192.
!
!***END PROLOGUE CHECK_NORESTART
!
!   Global variables.
!
      INTEGER, INTENT(IN)  :: DIMENS,NUMFUN,NUMRGN,JOB
      INTEGER, INTENT(OUT) :: INFORM
      INTEGER, DIMENSION(:), INTENT(IN) :: RGTYPE
      INTEGER, OPTIONAL, INTENT(IN OUT)  :: IFAIL
      INTEGER, OPTIONAL, INTENT(IN)     :: MINPTS,MAXPTS
      REAL(kind=stnd), OPTIONAL, INTENT(IN)  :: EPSABS,EPSREL,TUNE
!
!   Local variables
!
      INTEGER :: I
!
!***FIRST EXECUTABLE STATEMENT CHECK_NORESTART
!

      INFORM = 0
!
!     Check valid DIMENS.
!
          IF (DIMENS <= 0) THEN
                                INFORM = 8
          END IF
!
!     Check positive NUMFUN.
!
          IF (NUMFUN < 1) THEN
                               INFORM = INFORM + 16
          END IF
!
!     Check valid NUMRGN.
!
          IF (NUMRGN <= 0) THEN
                                INFORM = INFORM + 32
          END IF
!
!     Check valid region type
!
          DO I = 1,NUMRGN
              IF ((RGTYPE(I) /= Simplex) .and. (RGTYPE(I) /= Hyperrectangle)) THEN
                  INFORM = INFORM + 64
                  EXIT
              END IF
          END DO
!
!     Check if DIMENS and JOB or in agreement for special values of JOB
!
      IF ( ((JOB == 2) .OR. (JOB == 11)) .AND. (DIMENS > 3)) THEN
         INFORM = INFORM + 128
      END IF
!
!     Check valid limits on allowed number of function evaluations.
!
      IF (PRESENT(MAXPTS)) THEN
         IF (PRESENT(MINPTS)) THEN
            IF (MINPTS > MAXPTS) THEN
                                      INFORM = INFORM + 256
            END IF
         ELSE
            IF (MAXPTS < 1) THEN 
                                 INFORM = INFORM + 256
                            END IF
         END IF
      END IF
!
!     Check valid requested reliablitiy.
!
      IF (PRESENT(TUNE)) THEN
         IF ((TUNE < 0) .OR. (TUNE > 1)) THEN
                                              INFORM = INFORM + 512
         END IF
      END IF
!
!     Check valid accuracy requests.
!
      IF (PRESENT(EPSABS)) THEN
         IF (EPSABS < 0) THEN
                              INFORM = INFORM + 1024
         END IF
      ELSE IF (PRESENT(EPSREL)) THEN
         IF (EPSREL < 0) THEN
                              INFORM = INFORM + 1024
         END IF
      END IF
!
!     Check valid IFAIL
!
      IF (PRESENT(IFAIL)) THEN
         IF ((IFAIL < -1) .OR. (IFAIL > 1)) THEN
              INFORM = INFORM + 2048
              IFAIL = 0
         END IF
      END IF
!
!     Check valid JOB. JOB = 0 cannot occur at this stage
!    
      IF ((ABS(JOB) /= 1) .AND. (ABS(JOB) /= 2) .AND. (ABS(JOB) /= 11) .AND.  (ABS(JOB) /= 12)) THEN
                                                      INFORM = INFORM + 4096
      END IF
!
      RETURN
      END SUBROUTINE CHECK_NORESTART

END MODULE Check_Input
