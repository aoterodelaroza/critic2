! This file is F-compatible, except for upper/lower case conventions.
!--------------------------------------------------------------------
Module Region_Processor

USE Precision_Model, ONLY: stnd
USE internal_types
USE Subdivisions
USE CubatureRule_General

Implicit NONE

PRIVATE
PUBLIC :: Process_Region

CONTAINS
      SUBROUTINE Process_Region(DIMENS,CINFO,NRVERT,MAXSUB,NUMFUN,       &
                                Integrand,VEROLD,INFOLD,RINFOL,       &
                                MAXPTS,OUTSUB,VERNEW,INFNEW,RINFNE,   &
                                VALNEW,ERRNEW,NREVAL,IFAIL)
!***BEGIN PROLOGUE Process_Region
!***DATE WRITTEN   910507   (YYMMDD)
!***REVISION DATE  970605   (YYMMDD)
!***REVISION DATE  990527   (YYMMDD) (F conversion)
!***REVISION DATE  990624   (YYMMDD) (2/4/8 division of C3 activated)
!***REVISION DATE  010919   (YYMMDD) (and now this can be switched off)
!***AUTHOR
!          Ronald Cools, Dept. of Computer Science,
!          Katholieke Universiteit Leuven, Celestijnenlaan 200A,
!          B-3001 Heverlee, Belgium
!          Email:  Ronald.Cools@cs.kuleuven.ac.be
!
!***PURPOSE  Improve the estimate of an integral over a region.
!
!***DESCRIPTION
!
!   Input parameters
!   ----------------
!
!   DIMENS Integer.
!          The dimension of the region of integration.
!   NRVERT Integer.
!          The number of vertices to determine a region.
!   MAXSUB  Integer.
!          A region is divided into at most MAXSUB subregions.
!   NUMFUN Integer.
!          Number of components of the integral.
!   NIINFO Integer = size(INFOLD)
!          The number of integers used to save information about
!          the region.
!          Conventions for info-record:
!          info-record(5) = 1 if there was asymptotic behaviour when the
!                           region was processed before.
!                         = 0 otherwise
!          info-record(4) = suggestions for subdivision
!          info-record(3) = number of the original region where this
!                           region is a part of
!          info-record(2) = (volume of orinal region)/
!                           (volume of this region)
!          info-record(1) = type of region
!   NRINFO Integer = size(RINFOL)
!          The number of reals used to save information about
!          the region.
!          Conventions for info-record:
!          info-record(1) = volume
!
!   MAXPTS Integer.
!          The maximum number of function evaluations that Process_Region
!          may use.
!   Integrand Externally declared subroutine for computing
!          all components of the integrand at the given
!          evaluation point.
!          It must have parameters (X,NUMFUN,FUNVLS)
!          Input parameters:
!            X(1)   The x-coordinate of the evaluation point.
!            X(2)   The y-coordinate of the evaluation point.
!            ...
!            X(DIMENS) The z-coordinate of the evaluation point.
!            NUMFUN Integer that defines the number of
!                   components of I.
!          Output parameter:
!            FUNVLS Real array of dimension NUMFUN
!                   that defines NUMFUN components of the integrand.
!
!
!   Output parameters
!   -----------------
!
!   OUTSUB Integer.
!          The number of regions returned by this routine.
!   NREVAL Integer.
!          Number of function evaluations used by Process_Region.
!   IFAIL  Integer.
!          IFAIL = 0 for normal exit.
!          IFAIL = 1 if nothing was done because nothing can be
!                    done with MAXPTS function evaluations.
!          IFAIL = 6 if an unsupported type of region is given.
!          IFAIL = 7 or 6 if DIVIDE returns with an error
!           If IFAIL not equal to 0 then OUTSUB = 0
!
!***ROUTINES CALLED Rule_General,DIVIDE
!***END PROLOGUE Process_Region

!   Global variables
!
      INTERFACE 
         FUNCTION Integrand(NUMFUN,X) RESULT(Value)
            USE Precision_Model
            INTEGER, INTENT(IN) :: NUMFUN
            REAL(kind=stnd), DIMENSION(:), INTENT(IN) :: X
            REAL(kind=stnd), DIMENSION(NUMFUN) :: Value
         END FUNCTION Integrand
      END INTERFACE
      INTEGER, INTENT(IN) ::  DIMENS,NRVERT,MAXSUB,NUMFUN,MAXPTS
      INTEGER, DIMENSION(:), INTENT(IN) :: INFOLD
      INTEGER, INTENT(OUT)::  NREVAL,IFAIL,OUTSUB
      INTEGER, DIMENSION(:,:), INTENT(OUT):: INFNEW
      REAL(kind=stnd), DIMENSION(:,:), INTENT(IN) :: VEROLD
      REAL(kind=stnd), DIMENSION(:), INTENT(IN) :: RINFOL
      REAL(kind=stnd), DIMENSION(:,:), INTENT(OUT):: RINFNE,VALNEW,ERRNEW
      REAL(kind=stnd), DIMENSION(:,:,:), INTENT(OUT):: VERNEW
      TYPE(integrator_info), INTENT(IN) :: CINFO
!
!   Local variables
!
      INTEGER :: L,NUM,REQSUB
!
      NREVAL = 0
      REQSUB = MAXSUB
      IF ( (.NOT. CINFO%UNIFORM_SUBDIV) .AND.  (INFOLD(4) /= 0) ) THEN
         ! INFOLD(4) is set by Rule_C* ; not by Rule_T*
         IF (INFOLD(4) < 0) THEN
             ! Rule_C3 can return negative numbers if 4-division is
             ! suggested.
             REQSUB = MIN(REQSUB,4)
         ELSE IF (INFOLD(4) < 90) THEN
             ! Rule_Cn and Rule_C2 can only returns positive numbers
             !   and this suggests 2-division.
             ! Rule_C3 returns a number >90 if uniform subdivision is
             !   recommended.
             REQSUB = MIN(REQSUB,2)
         END IF
      END IF
!
!  Divide the region in (at most) REQSUB subregions.
!
      DO
        CALL DIVIDE(DIMENS,NRVERT,REQSUB,CINFO%UNIFORM_SUBDIV,NUMFUN,VEROLD,INFOLD, &
                    RINFOL,Integrand,OUTSUB,NUM,VERNEW,INFNEW,RINFNE,IFAIL)
        NREVAL = NREVAL + NUM
        IF (IFAIL == 0) THEN
            EXIT
        END IF
        IF ( CINFO%UNIFORM_SUBDIV .OR. REQSUB == 2) THEN 
            ! Now IFAIL /= 0 and we have no alternative
            OUTSUB = 0
            RETURN
        END IF
        REQSUB = 2  ! Retry, asking for a 2-division
      END DO

!
!  Apply the basic rule to each subregion.
!
      NUM = 0
      DO L = 1,OUTSUB
         NUM = NUM + Rule_Cost( DIMENS, INFNEW(1,L), CINFO%Key)
      END DO
      IF ( NREVAL + NUM > MAXPTS ) THEN
              IFAIL = 1
              OUTSUB = 0
              RETURN
      END IF
      DO L = 1,OUTSUB
          CALL Rule_General(DIMENS,CINFO,VERNEW(:,:,L),INFNEW(:,L),  & 
                            RINFNE(:,L),NUMFUN,Integrand,VALNEW(:,L),&
                            ERRNEW(:,L),NUM,IFAIL)
          NREVAL = NREVAL + NUM
          IF (IFAIL /= 0) THEN
              OUTSUB = 0
              RETURN
          END IF
      END DO
!
      RETURN
      END SUBROUTINE Process_Region

END MODULE Region_Processor
