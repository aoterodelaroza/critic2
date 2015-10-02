! This file is F-compatible, except for upper/lower case conventions.
!--------------------------------------------------------------------
Module CubatureRule_T3

USE Precision_Model, ONLY: stnd

Implicit NONE

PRIVATE

PUBLIC :: Rule_T3a
PRIVATE :: OrbT3_Sum

CONTAINS
    SUBROUTINE Rule_T3a(VER,VOLUME,NUMFUN,Integrand,BASVAL,RGNERR,NUM)
!
!***BEGIN PROLOGUE Rule_T3a
!***REFER TO DCUTET
!***REVISION DATE  970324   (YYMMDD) 
!***REVISION DATE  990528   (YYMMDD)  (F conversion)
!***PURPOSE  To compute basic integration rule values and
!            corresponding error estimates.
!***DESCRIPTION Rule_T3a computes basic integration rule values
!            for a vector of integrands over a tetrahedron.
!            Rule_T3a also computes estimates for the errors by
!            using several null rule approximations.
!   ON ENTRY
!
!   VER    Real array of dimension (3,4).
!          The coordinates of the vertices of the tetrahedron.
!          vertex i -> ( ver(1,i),ver(2,i),ver(3,i) )
!   NUMFUN Integer.
!          Number of components of the vector integrand.
!   Integrand Externally declared subroutine for computing
!            all components of the integrand at the given
!            evaluation point.
!            It must have parameters (X)
!            Input parameters:
!              X(1)   The x-coordinate of the evaluation point.
!              X(2)   The y-coordinate of the evaluation point.
!              X(3)   The z-coordinate of the evaluation point.
!
!   ON RETURN
!
!   BASVAL Real array of dimension NUMFUN.
!          The values for the basic rule for each component
!          of the integrand.
!   RGNERR Real array of dimension NUMFUN.
!          The error estimates for each component of the integrand.
!   NUM    Integer
!          The number of function evaluations used.
!
!***REFERENCES M. Beckers and A. Haegemans,
!              The construction of cubature formula for the tetrahedron
!              Report TW128, K.U. Leuven (1990).
!***ROUTINES CALLED
!                   OrbT3_Sum,Integrand
!***END PROLOGUE Rule_T3a
!
!   Parameters
!
!   ORBITS Integer
!          The number of orbits of the cubature formula and null rules
!   CRIVAL Real
!          The decision to choose the optimistic part of the error
!          estimator is based on CRIVAL
!   FACMED Real
!          FACMED is the safety coefficient used in the non-optimistic
!          part of the error estimator. FACMED is related to CRIVAL
!          and FACOPT.
!   FACOPT Real
!          FACOPT is the safety coefficient used in the optimistic part
!          of the error estimator.
!   K      Integer array of dimension (0:3) that contains the structure
!          parameters. K(I) = number of orbits of type I.
!   TYPE1  Real array of dimension (K(1)).
!          Contains the first homogeneous coordinate of the generators
!          of type 1
!   TYPE2  Real array of dimension (K(2)).
!          Contains the first homogeneous coordinate of the generators
!          of type 2
!   TYPE3  Real array of dimension (2,K(2)).
!          Contains the first two homogeneous coordinates of
!          the generators of type 3.
!   WEIGHT Real array of dimension (9,ORBITS).
!          The weights of the cubature formula and the null rules.
!          WEIGHT(1,1) ,..., WEIGHT(1,ORBITS) are the weights of the
!                cubature formula
!          WEIGHT(I,1) ,..., WEIGHT(I,ORBITS) for I > 1, are the weights
!                of the null rules
!
!
!   Global variables.
!
      INTEGER, INTENT(IN)                         :: NUMFUN
      INTEGER, INTENT(OUT)                        :: NUM
      REAL(kind=stnd), INTENT(IN)                 :: VOLUME
      REAL(kind=stnd), DIMENSION(:,:), INTENT(IN) :: VER
      REAL(kind=stnd), DIMENSION(:), INTENT(OUT)  :: BASVAL, RGNERR
      INTERFACE 
         FUNCTION Integrand(NUMFUN,X) RESULT(Value)
            USE Precision_Model
            INTEGER, INTENT(IN) :: NUMFUN
            REAL(kind=stnd), DIMENSION(:), INTENT(IN) :: X
            REAL(kind=stnd), DIMENSION(NUMFUN) :: Value
         END FUNCTION Integrand
      END INTERFACE
!
!   Constants
!
      INTEGER, PARAMETER :: ORBITS = 7
      REAL(kind=stnd), PARAMETER:: CRIVAL=0.5_stnd,         &
                                   FACMED=5,                &
                                   FACOPT=FACMED/CRIVAL,    &
                                   TRES=50*EPSILON(crival)
!
!  Cubature formula of degree 8 with 43 points
!
      INTEGER, DIMENSION(0:3), PARAMETER ::                         &
              K = (/1,3,1,2/)         ! Rule structure parameters
!
!  Information for the generators
!
      REAL(kind=stnd), DIMENSION(1:3), PARAMETER ::                 &
              TYPE1 = (/ 0.379510205167980387748057300876_stnd,     &
                         0.753689235068359830728182577696_stnd,     &
                         0.982654148484406008240470085259_stnd/)
      REAL(kind=stnd), DIMENSION(1:1), PARAMETER ::                 &
              TYPE2 = (/ 0.449467259981105775574375471447_stnd/)
      REAL(kind=stnd), DIMENSION(1:2,1:2), PARAMETER ::             &
              TYPE3 = RESHAPE( SOURCE=                              &
                   (/ 0.506227344977843677082264893876_stnd,        &
                      0.356395827885340437169173969841E-1_stnd,     &
                      0.736298458958971696943019005441_stnd,        &
                      0.190486041934633455699433285302_stnd/),      &
                   SHAPE=(/2,2/), ORDER=(/1,2/) )
!
!  Weights of the cubature formula
!
      REAL(kind=stnd), DIMENSION(1:ORBITS), PARAMETER ::            &
           W1 = (/                                                  &
          -0.123001131951839495043519102752_stnd,                   &
           0.855018349372014074906384482699E-1_stnd,                &
           0.118021998788034059253768205083E-1_stnd,                &
           0.101900465455732427902646736855E-2_stnd,                &
           0.274781029468036908044610867719E-1_stnd,                &
           0.342269148520915110408153517904E-1_stnd,                &
           0.128431148469725555789001180031E-1_stnd/)
!
!  Weights of the null rule of degree 5
!
      REAL(kind=stnd), DIMENSION(1:ORBITS), PARAMETER ::            &
           W2 = (/                                                  &
           0.211921237628032658308230999090_stnd,                   &
          -0.660207516445726284649283745987E-1_stnd,                &
           0.225058824086711710443385047042E-1_stnd,                &
          -0.375962972067425589765730699401E-3_stnd,                &
           0.710066020561055159657284834784E-2_stnd,                &
           0.156515256061747694921427149028E-2_stnd,                &
          -0.814530839643584660306807872526E-2_stnd/)
!
!  Weights of null rule of degree 4
!
      REAL(kind=stnd), DIMENSION(1:ORBITS), PARAMETER ::            &
           W3 = (/                                                  &
          -0.508105488137100551376844924797E-1_stnd,                &
           0.104596681151665328209751420525E-1_stnd,                &
           0.927471438532788763594989973184E-1_stnd,                &
           0.210489990008917994323967321174E-2_stnd,                &
           0.379184172251962722213408547663E-1_stnd,                &
          -0.111747242913563605790923001557E-1_stnd,                &
          -0.386541758762774673113423570465E-1_stnd/)
!
!  Weights of first null rule of degree 3
!
      REAL(kind=stnd), DIMENSION(1:ORBITS), PARAMETER ::            &
           W4 = (/                                                  &
          -0.775992773232808462404390159802E-1_stnd,                &
          -0.527453289659022924847298408064E-1_stnd,                &
           0.145876238555932704488677626554E-1_stnd,                &
           0.739374873393616192857532718429E-2_stnd,                &
          -0.374618791364332892611678523428E-1_stnd,                &
           0.538502846550653076078817013885E-1_stnd,                &
          -0.183980865177843057548322735665E-1_stnd/)
!
!  Weights of second null rule of degree 3
!
      REAL(kind=stnd), DIMENSION(1:ORBITS), PARAMETER ::            &
           W5 = (/                                                  &
           0.181767621501470154602720474731E-1_stnd,                &
           0.179938831310058580533178529022E-1_stnd,                &
           0.713210362750414891598257378898E-1_stnd,                &
          -0.443935688958258805893448212636E-1_stnd,                &
          -0.657639036547720234169662790056E-1_stnd,                &
          -0.101551807522541414699808460583E-1_stnd,                &
           0.265486188970540796821750584204E-1_stnd/)
!
!  Weights of null rule of degree 2
!
      REAL(kind=stnd), DIMENSION(1:ORBITS), PARAMETER ::            &
           W6 = (/                                                  &
          -0.867629853722843888927184699428E-1_stnd,                &
          -0.715881271235661902772072127812E-1_stnd,                &
           0.886720767790426261677273459523E-2_stnd,                &
          -0.577885573028655167063092577589E-1_stnd,                &
           0.430310167581202031805055255554E-1_stnd,                &
          -0.606467834856775537069463817445E-2_stnd,                &
           0.319492443333738343104163265406E-1_stnd/)
!
!   Weights of null rule of degree 1
!
      REAL(kind=stnd), DIMENSION(1:ORBITS), PARAMETER ::            &
           W7 = (/                                                  &
           0.510374015624925451319499382594E-1_stnd,                &
           0.463998830432033721597269299429E-1_stnd,                &
          -0.191086148397852799983451475821E-1_stnd,                &
          -0.973768821003670776204287367278E-1_stnd,                &
           0.180352562073914141268335496511E-1_stnd,                &
           0.277129527093489643801598303110E-1_stnd,                &
          -0.176218263109360550515567818653E-1_stnd/)
!
      REAL(kind=stnd), DIMENSION(1:7,1:ORBITS), PARAMETER ::        &
           WEIGHT =  RESHAPE( SOURCE= (/ W1,W2,W3,W4,W5,W6,W7/),    &
                            SHAPE=(/7,ORBITS/), ORDER=(/2,1/) )
!
!   Local variables.
!
      INTEGER :: J,NR,P,GENTYPE
      REAL(kind=stnd):: NOISE,DEG4,DEG3,DEG1,R2,R1,R
      REAL(kind=stnd), DIMENSION(NUMFUN) :: SUMVAL
      REAL(kind=stnd), DIMENSION(NUMFUN,6) :: NullRule
      REAL(kind=stnd), DIMENSION(3) :: Z
!
!***FIRST EXECUTABLE STATEMENT Rule_T3a
!
!  The number of points used by the cubature formula is
!     NUM    = K(0) + 4*K(1) + 6*K(2) + 12*K(3) = 43
      NUM = 43
!
!  Initialise BASVAL and NullRule
!
      BASVAL = 0 
      NullRule = 0
!
!  Compute contributions from orbits with 1, 4, 6 and 12  points
!
      P = 1
      DO GENTYPE = 0,3
          DO NR = 1,K(GENTYPE)
              SELECT CASE (GENTYPE)
              CASE (1) !       Generator ( z(1) , z(2), z(2) , z(2) )
                  Z(1) = TYPE1(NR)
                  Z(2) = (1-Z(1))/3
              CASE (2) !       Generator ( z(1) , z(1), z(2) , z(2) )
                  Z(1) = TYPE2(NR)
                  Z(2) = (1-2*Z(1))/2
              CASE (3) !       Generator ( z(1) , z(2), z(3) , z(3) )
                  Z(1:2) = TYPE3(1:2,NR)
                  Z(3) = (1-Z(1)-Z(2))/2
              END SELECT
              CALL OrbT3_Sum(GENTYPE,Z,VER,NUMFUN,Integrand,SUMVAL)
              BASVAL = BASVAL + WEIGHT(1,P)*SUMVAL
              DO J = 1,NUMFUN
                 NullRule(J,1:6) = NullRule(J,1:6) + WEIGHT(2:7,P)*SUMVAL(J)
              END DO
              P = P + 1
          END DO
      END DO
!
!  Compute error estimates
!
      DO J = 1,NUMFUN
          NOISE = ABS(BASVAL(J))*TRES
          DEG4 = SQRT(NullRule(J,1)**2+NullRule(J,2)**2)
          DEG3 = SQRT(NullRule(J,3)**2+NullRule(J,4)**2)
          IF (DEG4 <= NOISE) THEN
              RGNERR(J) = NOISE
          ELSE
              DEG1 = SQRT(NullRule(J,5)**2+NullRule(J,6)**2)
              IF (DEG3 /= 0) THEN 
                                    R1 = (DEG4/DEG3)**2
                             ELSE 
                                    R1 = 1
              END IF
              IF (DEG1 /= 0) THEN 
                                    R2 = DEG3/DEG1
                             ELSE 
                                    R2 = 1
              END IF
              R = MAX(R1,R2)
              IF (R >= CRIVAL) THEN 
                                    RGNERR(J) = FACMED*R*DEG4
                               ELSE 
                                    RGNERR(J) = FACOPT*(R**2)*DEG4
              END IF
              RGNERR(J) = MAX(NOISE,RGNERR(J))
          END IF
          RGNERR(J) = VOLUME*RGNERR(J)
          BASVAL(J) = VOLUME*BASVAL(J)
      END DO
      RETURN
      END SUBROUTINE Rule_T3a


      SUBROUTINE OrbT3_Sum(GENTYPE,GENER,VER,NUMFUN,Integrand,SUMVAL)
!***BEGIN PROLOGUE OrbT3_Sum
!***PURPOSE  To compute the sum of function values over all points
!            of an orbit.
!***DESCRIPTION
!   ON ENTRY
!
!   GENTYPE Integer
!          The type of the orbit.
!   GENER  Integer array of dimension (3).
!          The generator for the orbit in homogeneous coordinates.
!   VER    Real array of dimension (3,4).
!          The coordinates of the vertices of the tetrahedron.
!          vertex i -> ( ver(1,i),ver(2,i),ver(3,i) )
!   NUMFUN Integer.
!          Number of components of the vector integrand.
!   Integrand Externally declared subroutine for computing
!            all components of the integrand at the given
!            evaluation point.
!            It must have parameters (DIM,X,NUMFUN,FUNVLS)
!            Input parameters:
!              DIM = 3
!              X(1)   The x-coordinate of the evaluation point.
!              X(2)   The y-coordinate of the evaluation point.
!              X(3)   The z-coordinate of the evaluation point.
!              NUMFUN Integer that defines the number of
!                     components of the vector integrand.
!            Output parameter:
!              FUNVLS Real array of dimension NUMFUN
!                     that defines NUMFUN components of the integrand.

!   ON RETURN
!
!   SUMVAL Real array of dimension (NUMFUN).
!          The sum of function values over all points
!          of the given orbit.
!
!***END PROLOGUE  OrbT3_Sum
!
!  Global variables
!
      INTEGER, INTENT(IN) :: NUMFUN,GENTYPE
      REAL(kind=stnd), DIMENSION(:), INTENT(IN)   :: GENER
      REAL(kind=stnd), DIMENSION(:), INTENT(OUT)  :: SUMVAL
      REAL(kind=stnd), DIMENSION(:,:), INTENT(IN) :: VER
      INTERFACE 
         FUNCTION Integrand(NUMFUN,X) RESULT(Value)
            USE Precision_Model
            INTEGER, INTENT(IN) :: NUMFUN
            REAL(kind=stnd), DIMENSION(:), INTENT(IN) :: X
            REAL(kind=stnd), DIMENSION(NUMFUN) :: Value
         END FUNCTION Integrand
      END INTERFACE
!
!  Local variables
!
      INTEGER ::  J,NUMBER
      REAL(kind=stnd):: Z1,Z2,Z3
      REAL(kind=stnd), DIMENSION(3,12) :: X
      REAL(kind=stnd), DIMENSION(NUMFUN) :: WORK
!***FIRST EXECUTABLE STATEMENT  OrbT3_Sum
      SELECT CASE (GENTYPE)
!
!  Generator with homogeneous coordinates (1/4,1/4,1/4,1/4)
!
      CASE (0)
          NUMBER = 1
          X(:,1) = SUM( VER, DIM=2 )/4
!
!  Generator with homogeneous coordinates (z1,z2,z2,z2)
!
      CASE (1)
          NUMBER = 4
          Z1 = GENER(1)
          Z2 = GENER(2)
          X(:,1) = Z1*VER(:,1) + Z2* (VER(:,2)+VER(:,3)+VER(:,4))
          X(:,2) = Z1*VER(:,2) + Z2* (VER(:,1)+VER(:,3)+VER(:,4))
          X(:,3) = Z1*VER(:,3) + Z2* (VER(:,2)+VER(:,1)+VER(:,4))
          X(:,4) = Z1*VER(:,4) + Z2* (VER(:,2)+VER(:,3)+VER(:,1))
!
!  Generator with homogeneous coordinates (z1,z1,z2,z2)
!
      CASE (2)
        NUMBER = 6
        Z1 = GENER(1)
        Z2 = GENER(2)
        X(:,1) = Z1* (VER(:,1)+VER(:,2)) + Z2* (VER(:,3)+VER(:,4))
        X(:,2) = Z1* (VER(:,1)+VER(:,3)) + Z2* (VER(:,2)+VER(:,4))
        X(:,3) = Z1* (VER(:,1)+VER(:,4)) + Z2* (VER(:,3)+VER(:,2))
        X(:,4) = Z1* (VER(:,2)+VER(:,3)) + Z2* (VER(:,1)+VER(:,4))
        X(:,5) = Z1* (VER(:,2)+VER(:,4)) + Z2* (VER(:,1)+VER(:,3))
        X(:,6) = Z1* (VER(:,3)+VER(:,4)) + Z2* (VER(:,1)+VER(:,2))
!
!  Generator with homogeneous coordinates (z1,z2,z3,z3)
!
      CASE (3)
        NUMBER = 12
        Z1 = GENER(1)
        Z2 = GENER(2)
        Z3 = GENER(3)
        X(:,1) = Z1*VER(:,1) + Z2*VER(:,2) + Z3* (VER(:,3)+VER(:,4))
        X(:,2) = Z1*VER(:,1) + Z2*VER(:,3) + Z3* (VER(:,2)+VER(:,4))
        X(:,3) = Z1*VER(:,1) + Z2*VER(:,4) + Z3* (VER(:,2)+VER(:,3))
        X(:,4) = Z1*VER(:,2) + Z2*VER(:,1) + Z3* (VER(:,3)+VER(:,4))
        X(:,5) = Z1*VER(:,2) + Z2*VER(:,3) + Z3* (VER(:,1)+VER(:,4))
        X(:,6) = Z1*VER(:,2) + Z2*VER(:,4) + Z3* (VER(:,1)+VER(:,3))
        X(:,7) = Z1*VER(:,3) + Z2*VER(:,1) + Z3* (VER(:,2)+VER(:,4))
        X(:,8) = Z1*VER(:,3) + Z2*VER(:,2) + Z3* (VER(:,1)+VER(:,4))
        X(:,9) = Z1*VER(:,3) + Z2*VER(:,4) + Z3* (VER(:,1)+VER(:,2))
        X(:,10) = Z1*VER(:,4) + Z2*VER(:,1) + Z3* (VER(:,2)+VER(:,3))
        X(:,11) = Z1*VER(:,4) + Z2*VER(:,2) + Z3* (VER(:,1)+VER(:,3))
        X(:,12) = Z1*VER(:,4) + Z2*VER(:,3) + Z3* (VER(:,1)+VER(:,2))
      END SELECT
      SUMVAL = Integrand(NUMFUN,X(:,1))
      DO J = 2,NUMBER
         WORK = Integrand(NUMFUN,X(:,J))
         SUMVAL(1:NUMFUN) = SUMVAL(1:NUMFUN) + WORK(1:NUMFUN)
      END DO
      RETURN
      END SUBROUTINE OrbT3_Sum

END MODULE CubatureRule_T3
