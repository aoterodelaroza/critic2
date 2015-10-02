Module CubatureRule_General

USE Precision_Model, ONLY: stnd
USE internal_types
USE QuadratureRule
USE CubatureRule_T2
USE CubatureRule_T3
USE CubatureRule_Tn
USE CubatureRule_C2
USE CubatureRule_C3
USE CubatureRule_Cn

Implicit NONE

PRIVATE

PUBLIC :: Rule_General, Rule_Cost

CONTAINS
   SUBROUTINE Rule_General(DIMENS,CINFO,VER,IINFO,RINFO,NUMFUN,Integrand, &
                           BASVAL,RGNERR,NUM,IFAIL)

!***BEGIN PROLOGUE Rule_General
!***DATE WRITTEN   901114   (YYMMDD)
!***REVISION DATE  970507   (YYMMDD)
!***REVISION DATE  980330   (YYMMDD) (1D added)
!***REVISION DATE  980406   (YYMMDD) (dcuhre added)
!***REVISION DATE  990527   (YYMMDD) (F conversion)
!***REVISION DATE  000814   (YYMMDD) (rule selection for Cn changed)
!***REVISION DATE  010829   (YYMMDD) (init IFAIL)
!***REVISION DATE  020716   (YYMMDD) (rule selection for Cn changed))
!***AUTHOR
!          Ronald Cools, Dept. of Computer Science,
!          Katholieke Universiteit Leuven, Celestijnenlaan 200A,
!          B-3001 Heverlee, Belgium
!          Email:  Ronald.Cools@cs.kuleuven.ac.be
!
!***PURPOSE  To compute basic integration rule values and
!            corresponding error estimates.
!***DESCRIPTION Rule_General selects a basic integration rule
!            suitable for the given region.
!
!   Input parameters
!   ----------------
!
!   DIMENS Integer
!          The dimension of the region
!   VER    Real array of dimension (DIMENS,NRVERT).
!          The coordinates of the vertices of the region.
!   NUMFUN Integer
!          Number of components of the vector integrand.
!   CINFO  Type integrator_info
!          Paramaters to select proper integration rule.
!   IINFO  Integer array
!   RINFO  Real array
!          The 2 arrays contain additional information about the
!          subregion. See MODULE Global_Adaptive_Algorithm
!          for details.
!   Integrand Externally declared subroutine for computing
!            all components of the integrand at the given
!            evaluation point.
!            It must have parameters (X,NUMFUN,FUNVLS)
!            Input parameters:
!              X(1)      The x-coordinate of the evaluation point.
!              X(2)      The y-coordinate of the evaluation point.
!              NUMFUN Integer that defines the number of
!                     components of I.
!            Output parameter:
!              FUNVLS Real array of dimension NUMFUN
!                     that defines NUMFUN components of the integrand.
!
!
!   Output parameters
!   -----------------
!
!   BASVAL Real array of dimension NUMFUN.
!          The values for the basic rule for each component
!          of the integrand.
!   RGNERR Real array of dimension NUMFUN.
!          The error estimates for each component of the integrand.
!   NUM    Integer.
!          The number of function evaluations used.
!   IFAIL  Integer.
!          IFAIL = 0 on normal exit.
!          IFAIL = 5 if a given type of region is not supported.
!   IINFO  Integer array.
!   RINFO  Real array.
!          The 2 arrays contain additional information about the
!          subregion. See MODULE Global_Adaptive_Algorithm
!          for details.
!
!***ROUTINES CALLED Integrand
!***END PROLOGUE Rule_General
!
!   Global variables.
!
   INTERFACE 
      FUNCTION Integrand(NUMFUN,X) RESULT(Value)
         USE Precision_Model
         INTEGER, INTENT(IN) :: NUMFUN
         REAL(kind=stnd), DIMENSION(:), INTENT(IN) :: X
         REAL(kind=stnd), DIMENSION(NUMFUN) :: Value
      END FUNCTION Integrand
   END INTERFACE
    INTEGER, INTENT(IN)               :: NUMFUN,DIMENS
    INTEGER, INTENT(OUT)              :: NUM,IFAIL
    INTEGER, DIMENSION(:), INTENT(IN OUT) :: IINFO
    REAL(kind=STND), DIMENSION(:,:), INTENT(IN) :: VER
    REAL(kind=STND), DIMENSION(:), INTENT(OUT) :: BASVAL, RGNERR
    REAL(kind=STND), DIMENSION(:), INTENT(IN OUT) :: RINFO
    TYPE(integrator_info), INTENT(IN) :: cinfo
!
!  Local variables
!
!   KEY    Integer.
!          Rule selection parameter for DREPRO
!   TUNE   Real.
!          Requested reliability of DREPRO: 0 <= TUNE <= 1
!   RGTYPE  Integer
!           Indicates the type of region:
!           RGTYPE = 1 => the region is a simplex
!           RGTYPE = 2 => the region is a hyperrectangle
!           RGTYPE = 3 => the region is an octahedron
!           Note: If DIMENS=1 then RGTYPE=1 and RGTYPE=2 are equivalent
!
   INTEGER         :: RGTYPE,KEY
   REAL(kind=STND) :: TUNE
!
   KEY = CINFO%KEY
   TUNE = CINFO%TUNE
   RGTYPE = IINFO(1)
   IFAIL = 0
   SELECT CASE (RGTYPE)
   CASE (Simplex)
       IF (DIMENS == 1) THEN
          CALL Dqk_drv(KEY,VER,NUMFUN,Integrand,BASVAL,RGNERR,NUM)
       ELSE IF (DIMENS > 3) THEN
          CALL Rule_Tn(TUNE,DIMENS,VER,RINFO(1),NUMFUN,Integrand,KEY,BASVAL,RGNERR,NUM)
       ELSE IF ((KEY >= 1) .AND. (KEY <= 4)) THEN
          CALL Rule_Tn(TUNE,DIMENS,VER,RINFO(1),NUMFUN,Integrand,KEY,BASVAL,RGNERR,NUM)
       ELSE IF (DIMENS == 2) THEN
          CALL Rule_T2a(VER,RINFO(1),NUMFUN,Integrand,BASVAL,RGNERR,NUM)
       ELSE ! IF (DIMENS == 3) THEN
          CALL Rule_T3a(VER,RINFO(1),NUMFUN,Integrand,BASVAL,RGNERR,NUM)
       END IF
   CASE (Hyperrectangle)
       IF (DIMENS == 1) THEN
          CALL Dqk_drv(KEY,VER,NUMFUN,Integrand,BASVAL,RGNERR,NUM)
       ELSE IF (DIMENS > 3) THEN
          CALL Rule_Cn(KEY,DIMENS,VER,IINFO,RINFO(1),NUMFUN,Integrand,BASVAL,RGNERR,NUM)
       ELSE IF ((KEY >= 3) .AND. (KEY <= 4)) THEN
          CALL Rule_Cn(KEY,DIMENS,VER,IINFO,RINFO(1),NUMFUN,Integrand,BASVAL,RGNERR,NUM)
       ELSE IF (DIMENS == 2) THEN
          CALL Rule_C2a(VER,IINFO,RINFO(1),NUMFUN,Integrand,BASVAL,RGNERR,NUM)
       ELSE ! IF (DIMENS == 3) THEN
          CALL Rule_C3a(VER,IINFO,RINFO(1),NUMFUN,Integrand,BASVAL,RGNERR,NUM)
       END IF
   CASE DEFAULT
       IFAIL = 5
       NUM = 0
   END SELECT
   RETURN
   END SUBROUTINE Rule_General

   FUNCTION Rule_Cost( DIMENS, RGTYPE, KEY ) RESULT(RULCLS)
!
!     Integer function for computing the number of function values
!     needed by the local integration rule.
!     Input parameters are NOT checked!
!
! Global variables
!
!     DIMENS Integer number of dimensions.
!     RGTYPE Integer type of integration region.
!     KEY    Integer type of integration rule.
!
   INTEGER, INTENT(IN) :: DIMENS, RGTYPE, KEY
   INTEGER :: RULCLS
!
!     Local Variables
!
   INTEGER :: NKEY
 
   SELECT CASE (RGTYPE)

   CASE (Simplex) 
      IF ( DIMENS == 1 ) THEN
         SELECT CASE (KEY)
            CASE(:1)
            RULCLS = 15
            CASE(2)
            RULCLS = 21
            CASE(3)
            RULCLS = 31
            CASE(4)
            RULCLS = 41
            CASE(5)
            RULCLS = 51
            CASE(6:)
            RULCLS = 61
         END SELECT
      ELSE IF ( (DIMENS > 3) .OR. ((KEY >= 1) .AND. (KEY <= 4)) ) THEN
         !
         !  Compute RULCLS for DIMENS-simplex rules.
         !
         IF ( KEY == 0 ) THEN 
                               NKEY = 3
                         ELSE 
                               NKEY = KEY
         END IF
         ! First count the Grundmann and Moller points.
         RULCLS = DIMENS + 2
         IF ( NKEY > 1 ) THEN
              RULCLS = ( DIMENS + 3 )*RULCLS/2
         END IF
         IF ( NKEY > 2 ) THEN
              RULCLS = ( DIMENS + 4 )*RULCLS/3
         END IF
         IF ( NKEY > 3 ) THEN
              RULCLS = ( DIMENS + 5 )*RULCLS/4
         END IF
         ! Add those from the degree 5 Stroud rule.
         RULCLS = RULCLS + DIMENS + 1
         IF ( NKEY > 1 ) THEN
              RULCLS = RULCLS +     DIMENS + 1
         END IF
         IF ( NKEY > 2 ) THEN
              RULCLS = RULCLS +   ( DIMENS + 1 )*DIMENS 
         END IF
         ! Add those of the degree 7 Mysovskikh rule.
         IF ( NKEY > 3 ) THEN
              RULCLS = RULCLS + 3*( DIMENS + 1 )*( DIMENS + 2 )/2
         END IF
         ! subtract a generator if DIMENS == 3
         IF ((DIMENS == 3) .AND. (NKEY > 2)) THEN
              RULCLS = RULCLS - ((DIMENS+1)*DIMENS)/2
         END IF
      ELSE IF ( DIMENS == 2 ) THEN
         RULCLS = 37
      ELSE ! IF ( DIMENS == 3 ) THEN
         RULCLS = 43
      END IF

   CASE (Hyperrectangle) 
      !
      ! Compute RULCLS for DIMENS-hyperrectangle rules.
      !
      IF ( DIMENS == 1) THEN
         SELECT CASE( KEY )
            CASE(:1)
            RULCLS = 15
            CASE(2)
            RULCLS = 21
            CASE(3)
            RULCLS = 31
            CASE(4)
            RULCLS = 41
            CASE(5)
            RULCLS = 51
            CASE(6:)
            RULCLS = 61
         END SELECT
      ELSE IF ((DIMENS > 3) .OR. (KEY == 3) .OR. (KEY == 4)) THEN
         IF ( KEY == 4 ) THEN
            RULCLS = 1 + 4*2*DIMENS + 2*DIMENS* (DIMENS-1) + 4*DIMENS* (DIMENS-1) +     &
                  4*DIMENS* (DIMENS-1)* (DIMENS-2)/3 + 2**DIMENS
         ELSE
            RULCLS = 1 + 3*2*DIMENS + 2*DIMENS* (DIMENS-1) + 2**DIMENS
         END IF
      ELSE IF (DIMENS == 2) THEN
         RULCLS = 37
      ELSE ! IF (DIMENS == 3) THEN
         RULCLS = 89
      END IF
   END SELECT
   RETURN
   END FUNCTION Rule_Cost

END MODULE CubatureRule_General
