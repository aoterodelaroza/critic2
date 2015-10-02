! This file is F-compatible, except for upper/lower case conventions.
!--------------------------------------------------------------------
Module CubatureRule_C2

USE Precision_Model, ONLY: stnd

Implicit NONE

PRIVATE

PUBLIC :: Rule_C2a

CONTAINS
SUBROUTINE Rule_C2a(VER,INFOLD,AREA,NUMFUN,Integrand,BASVAL,RGNERR,NUM)
!
!***BEGIN PROLOGUE Rule_C2a
!***PURPOSE  To compute basic integration rule values and
!            corresponding error estimates.
! ***REVISION DATE  950531   (YYMMDD) (Fortran90 transformation)
! ***REVISION DATE  990527   (YYMMDD) (F transformation)
! ***AUTHOR
!          Ronald Cools, Dept. of Computer Science,
!          Katholieke Universiteit Leuven, Celestijnenlaan 200A,
!          B-3001 Heverlee, Belgium
!          Email:  ronald@cs.kuleuven.ac.be
!
! ***REFERENCES
!          The cubature formula of degree 13 with 37 points is from
!          Rabinowitz & Richter. The tuning of the error estimator
!          is described in:
!          R. Cools.
!          "The subdivision strategy and reliablity in adaptive
!           integration revisited."
!          Report TW 213, Dept. of Computer Science, K.U.Leuven, 1994.
!
!***DESCRIPTION Rule_C2a computes basic integration rule values
!            for a vector of integrands over a rectangular region.
!            Rule_C2a also computes estimates for the errors by
!            using several null rule approximations.
!   ON ENTRY
!
!   VER    Real array of dimension (2,3).
!          The coordinates of the vertices of the parallellogram.
!   NUMFUN Integer.
!          Number of components of the vector integrand.
!   INFOLD Integer array
!   Integrand Externally declared subroutine for computing
!            all components of the integrand at the given
!            evaluation point.
!            It must have parameters (DIM,X,NUMFUN,FUNVLS)
!            Input parameters:
!              DIM = 2
!              X(1)      The x-coordinate of the evaluation point.
!              X(2)      The y-coordinate of the evaluation point.
!              NUMFUN Integer that defines the number of
!                     components of I.
!            Output parameter:
!              FUNVLS Real array of dimension NUMFUN
!                     that defines NUMFUN components of the integrand.
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
!   INFOLD Integer array
!
!***ROUTINES CALLED Integrand
!***END PROLOGUE Rule_C2a
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
      INTEGER, INTENT(IN) :: NUMFUN
      INTEGER, INTENT(OUT) :: NUM
      INTEGER, DIMENSION(:), INTENT(IN OUT) :: INFOLD
      REAL(kind=stnd), INTENT(IN) :: AREA
      REAL(kind=stnd), DIMENSION(:,:), INTENT(IN) :: VER
      REAL(kind=stnd), DIMENSION(:), INTENT(OUT) :: BASVAL, RGNERR
!
!   Parameters
!
      INTEGER, DIMENSION(0:3), PARAMETER :: &
              K = (/1,2,3,2/)                ! Rule structure parameters

      INTEGER, PARAMETER ::                 &
              ORBITS = 8                     ! Number of orbits in rule

      REAL(kind=stnd), PARAMETER ::         &
              HALF = 0.5_stnd,              &
              FOUR = 4.0_stnd,              &
              CRIVAL = 0.4_stnd,            &
              FACMED = 8.0_stnd,            &
              FACOPT = FACMED/CRIVAL**2,    &
              TRES = 50*EPSILON(HALF),      &
              CUTOFF = 1.0E-4_stnd ,        &
              DFCLEV  = 0.55_stnd

      REAL(kind=stnd), DIMENSION(0:2), PARAMETER ::     &
              DFC = (/2.97397430397053625382_stnd,      &
                      1.0_stnd,                         &
                     -2.48698715198526812691_stnd  /)
!
!  Cubature formula of degree 13 with 37 points (Rabinowitz & Richter)
!
!
!  Information for the generators
!
      INTEGER :: I
      REAL(kind=stnd), DIMENSION(1:2), PARAMETER ::                    &
              TYPE1 = (/ 0.9909890363004326469792722978603_stnd,  &
                         0.6283940712305315063814483471116_stnd  /)

      REAL(kind=stnd), DIMENSION(1:3), PARAMETER ::                    &
              TYPE2 = (/ 0.9194861553393073086142137772149_stnd,  &
                         0.6973201917871173078084506730937_stnd,  &
                         0.3805687186904854497424188074662_stnd  /)
      REAL(kind=stnd), DIMENSION(1:2,1:2), PARAMETER ::                &
              TYPE3 = RESHAPE( SOURCE=                            &
                       (/ 0.9708504361720225062147290554088_stnd, &
                          0.6390348393207252159077623446225_stnd, &
                          0.8623637916722781475018696425693_stnd, &
                          0.3162277660168700033875075593701_stnd /),&
                          SHAPE=(/2,2/) )

!       The weights of the basic rule and the null rules.
!       WEIGHT(1,1),...,WEIGHT(1,ORBITS) are weights for the basic rule.
!       WEIGHT(I,1),...,WEIGHT(I,ORBITS) for I>1 are null rule weights.
!
! 
!  Weights of the cubature formula.
!
       REAL(kind=stnd), DIMENSION(ORBITS), PARAMETER ::           &
         W1 = (/                                                  &
          2.995235559387052215463143056692E-1_stnd ,              &
          3.311006686692356205977471655046E-2_stnd ,              &
          1.802214941550624038355347399683E-1_stnd ,              &
          3.916727896035153300761243260674E-2_stnd ,              &
          1.387748348777288706306435595057E-1_stnd ,              &
          2.268881207335707037147066705814E-1_stnd ,              &
          3.657395765508995601240002438981E-2_stnd ,              &
          1.169047000557533546701746277951E-1_stnd /)
!
!  Weights of the rules of degree 7, 7, 5 , 5 , 3 , 3  and 1.
!
       REAL(kind=stnd), DIMENSION(ORBITS), PARAMETER ::           &
         W2 = (/                                                  &
          7.610781847149629154716409791983E-2_stnd ,              &
          1.486101247399760261471935168346E-1_stnd ,              &
         -2.077685631717747007172983323970E-1_stnd ,              &
          6.850758313011924198538315395405E-2_stnd ,              &
          2.024205813317813585572881715385E-1_stnd ,              &
          1.108627473745508429879249169864E-1_stnd ,              &
         -1.187411393304862640859204217487E-1_stnd ,              &
         -5.208857468077715683772080394959E-2_stnd /)
!
       REAL(kind=stnd), DIMENSION(ORBITS), PARAMETER ::           &
         W3 = (/                                                  &
          4.016494861405949747097510013162E-2_stnd ,              &
         -1.093132962444079541048635452881E-1_stnd ,              &
         -2.270251673633777452624380129694E-1_stnd ,              &
          1.231674163356097016086203579325E-2_stnd ,              &
         -1.420402526499201540699111172200E-1_stnd ,              &
          1.189080551229557928776504129312E-1_stnd ,              &
         -4.482039658150474743804189300793E-3_stnd ,              &
          1.730383808319875827592824151609E-1_stnd /)
!
       REAL(kind=stnd), DIMENSION(ORBITS), PARAMETER ::           &
         W4 = (/                                                  &
         -5.643905795781771973971259866415E-1_stnd ,              &
          2.878418073676293225652331648545E-2_stnd ,              &
          1.159354231997583294689565314470E-1_stnd ,              &
          1.376081498690624477894043101438E-1_stnd ,              &
         -7.909780225340130915490382973570E-2_stnd ,              &
          1.174335441429478112778176601234E-1_stnd ,              &
         -1.107251942334134124782600707843E-1_stnd ,              &
          2.094226883312045633400182488252E-2_stnd /)
!
       REAL(kind=stnd), DIMENSION(ORBITS), PARAMETER ::           &
         W5 = (/                                                  &
         -2.269001713589584730602581694579E-1_stnd ,              &
          2.976190892690301120078774620049E-2_stnd ,              &
         -7.440193483272787588251423144751E-2_stnd ,              &
         -1.224665989043784131260454301280E-1_stnd ,              &
         -4.857910454732976198562745578156E-2_stnd ,              &
          2.228157325962656425537280474671E-1_stnd ,              &
          1.459764751457503859063666414952E-1_stnd ,              &
         -1.211789553452468781539987084682E-1_stnd /)
!
       REAL(kind=stnd), DIMENSION(ORBITS), PARAMETER ::           &
         W6 = (/                                                  &
         -3.326760468009974589269134283992E-1_stnd ,              &
          1.796655319904795478676993902115E-1_stnd ,              &
         -4.389976396805911868560791966472E-2_stnd ,              &
         -2.295841771339316497310760908889E-1_stnd ,              &
          6.182618387692816082856552878852E-2_stnd ,              &
         -1.202703885325137746461829140891E-1_stnd ,              &
          5.109536580363550180208564374234E-3_stnd ,              &
          1.126062761533095493689566169969E-1_stnd /)
!
       REAL(kind=stnd), DIMENSION(ORBITS), PARAMETER ::           &
         W7 = (/                                                  &
          2.290638530086106512999345512401E-1_stnd ,              &
          2.702070398116919449911037051753E-1_stnd ,              &
         -9.078047988731123605988441792069E-3_stnd ,              &
          4.618480310858703283999169489655E-2_stnd ,              &
         -2.598231009547631799096616255056E-1_stnd ,              &
         -2.518433931146441037986247681820E-2_stnd ,              &
         -1.257796993152456033984707367389E-2_stnd ,              &
         -2.720818902721190304043617320910E-2_stnd /)
!
       REAL(kind=stnd), DIMENSION(ORBITS), PARAMETER ::           &
         W8 = (/                                                  &
          2.746908885094872977794932213372E-1_stnd ,              &
         -1.149427039769738298032807785523E-2_stnd ,              &
          1.596178537820019535731955591283E-1_stnd ,              &
         -2.180626972663360443142752377527E-1_stnd ,              &
         -8.711748038292630173597899697063E-3_stnd ,              &
          1.902786182960269617633915869710E-1_stnd ,              &
         -1.189840649092108827784089292890E-1_stnd ,              &
          2.883382565767354162177931122471E-2_stnd /)

      REAL(kind=stnd), DIMENSION(1:8,1:ORBITS), PARAMETER ::      &
          WEIGHT = RESHAPE( SOURCE= (/ W1,W2,W3,W4,W5,W6,W7,W8 /),&
              SHAPE=(/8,ORBITS/), ORDER=(/2,1/) )
!
!   Local variables.
!
      INTEGER :: J,NUMBER,GENTYPE,NR,P
      REAL(kind=stnd):: R1,R2,R3,R,NOISE,DEG7,DEG5,DEG3,DEG1,   &
                        DIFFX,DIFFY,Z1,Z2
      REAL(kind=stnd), DIMENSION(2,8) :: X
      REAL(kind=stnd), DIMENSION(NUMFUN,7) :: NullRule
      REAL(kind=stnd), DIMENSION(NUMFUN) :: SUMVAL
!
!***FIRST EXECUTABLE STATEMENT Rule_C2a
!
!  The number of points used by the cubature formula is
!     NUM    = K(0) + 4*K(1) + 4*K(2) + 8*K(3)
      NUM = 37
!
!
!  Initialise BASVAL and NullRule
!
      BASVAL = 0
      NullRule = 0
      P = 1
!
!  Compute contributions from orbits with 1, 4 and 8 points
!
      DO GENTYPE = 0,3
          DO NR = 1,K(GENTYPE)
              SELECT CASE (GENTYPE)
              CASE (0) !                 Generator ( 0 , 0 )
                  NUMBER = 1
                  X(:,1) = (VER(:,2)+VER(:,3))*HALF

              CASE (1) !                 Generator ( z1 , 0 )
                  Z1 = TYPE1(NR)
                  NUMBER = 4
                  Z1 = Z1*HALF
                  X(:,1) = -VER(:,1)*Z1 + VER(:,2)*HALF +       &
                               VER(:,3)* (Z1+HALF)
                  X(:,2) = VER(:,1)*Z1 + VER(:,2)*HALF +        &
                               VER(:,3)* (-Z1+HALF)
                  X(:,3) = VER(:,1)*Z1 + VER(:,2)* (-Z1+HALF) + &
                               VER(:,3)*HALF
                  X(:,4) = -VER(:,1)*Z1 + VER(:,2)* (Z1+HALF) + &

                               VER(:,3)*HALF
              CASE (2) !                 Generator ( z(1) , z(1) )
                  Z1 = TYPE2(NR)
                  NUMBER = 4
                  Z1 = Z1*HALF
                  X(:,1) = -2*VER(:,1)*Z1 + VER(:,2)* (HALF+Z1) +&
                               VER(:,3)* (Z1+HALF)
                  X(:,2) = VER(:,2)* (HALF+Z1) +  VER(:,3)* (-Z1+HALF)
                  X(:,3) = VER(:,2)* (HALF-Z1) + VER(:,3)* (Z1+HALF)
                  X(:,4) = 2*VER(:,1)*Z1 + VER(:,2)* (HALF-Z1) + &
                               VER(:,3)* (-Z1+HALF)

              CASE (3) !                 Generator ( z(1) , z(2) )
                  Z1 = TYPE3(1,NR)*HALF
                  Z2 = TYPE3(2,NR)*HALF
                  NUMBER = 8
                  X(:,1) = VER(:,1)* (-Z1-Z2) +                     &
                               VER(:,2)* (HALF+Z2) + VER(:,3)* (HALF+Z1)
                  X(:,2) = VER(:,1)* (+Z1-Z2) +                     &
                               VER(:,2)* (HALF+Z2) + VER(:,3)* (HALF-Z1)
                  X(:,3) = VER(:,1)* (-Z1+Z2) +                     &
                               VER(:,2)* (HALF-Z2) + VER(:,3)* (HALF+Z1)
                  X(:,4) = VER(:,1)* (+Z1+Z2) +                     &
                               VER(:,2)* (HALF-Z2) + VER(:,3)* (HALF-Z1)
                  X(:,5) = VER(:,1)* (-Z1-Z2) +                     &
                               VER(:,2)* (HALF+Z1) + VER(:,3)* (HALF+Z2)
                  X(:,6) = VER(:,1)* (+Z2-Z1) +                     &
                               VER(:,2)* (HALF+Z1) + VER(:,3)* (HALF-Z2)
                  X(:,7) = VER(:,1)* (-Z2+Z1) +                     &
                               VER(:,2)* (HALF-Z1) + VER(:,3)* (HALF+Z2)
                  X(:,8) = VER(:,1)* (+Z1+Z2) +                     &
                               VER(:,2)* (HALF-Z1) + VER(:,3)* (HALF-Z2)
              END SELECT
!             CALL Integrand(2,X(1,1),NUMFUN,SUMVAL)
              SUMVAL = Integrand(NUMFUN,X(:,1))
              SELECT CASE (GENTYPE)
              CASE (0)
                  DIFFy = SUMVAL(1)*DFC(0)
                  DIFFx = DIFFy
              CASE (1)
                  DIFFy = DIFFy + SUMVAL(1)*DFC(NR)
              END SELECT

              DO J = 2,NUMBER
                  RGNERR = Integrand(NUMFUN,X(:,J))
!                 CALL Integrand(2,X(1,J),NUMFUN,RGNERR)
                  IF (GENTYPE == 1) THEN
                      IF (J <= 2) THEN
                          DIFFy = DIFFy + RGNERR(1)*DFC(NR)
                      ELSE
                          DIFFx = DIFFx + RGNERR(1)*DFC(NR)
                      END IF
                  END IF
                  DO I = 1,NUMFUN
                      SUMVAL(I) = SUMVAL(I) + RGNERR(I)
                  END DO
              END DO
              DO J = 1,NUMFUN
                  BASVAL(J) = BASVAL(J) + WEIGHT(1,P)*SUMVAL(J)
                  DO I = 1,7
                      NullRule(J,I) = NullRule(J,I) + WEIGHT(I+1,P)*SUMVAL(J)
                  END DO
              END DO
              P = P + 1
          END DO
      END DO
!
!  Decide on future subdivision direction
!
      DIFFy = ABS(DIFFy)
      DIFFx = ABS(DIFFx)
      IF (MAX(DIFFy,DIFFx) < CUTOFF) THEN
          INFOLD(4) = 0
      ELSE IF (DIFFy < DFCLEV*DIFFx) THEN
          INFOLD(4) = 1
      ELSE IF (DIFFx < DFCLEV*DIFFy) THEN
          INFOLD(4) = 2
      ELSE
          INFOLD(4) = 0
      END IF
!
!    Compute errors.
!
      DO J = 1,NUMFUN
          NOISE = ABS(BASVAL(J))*TRES
          DEG7 = SQRT(NullRule(J,1)**2+NullRule(J,2)**2)
          IF (DEG7 <= NOISE) THEN
              RGNERR(J) = NOISE
          ELSE
              DEG5 = SQRT(NullRule(J,3)**2+NullRule(J,4)**2)
              DEG3 = SQRT(NullRule(J,5)**2+NullRule(J,6)**2)
              DEG1 = SQRT(NullRule(J,7)**2+NullRule(J,6)**2)
              IF (DEG5 /= 0) THEN
                  R1 = DEG7/DEG5
              ELSE
                  R1 = 1
              END IF
              IF (DEG3 /= 0) THEN
                  R2 = DEG5/DEG3
              ELSE
                  R2 = 1
              END IF
              IF (DEG1 /= 0) THEN
                  R3 = DEG3/DEG1
              ELSE
                  R3 = 1
              END IF
              R = MAX(R1,R2,R3)
              IF (R >= 1) THEN
                  INFOLD(5) = 0
                  RGNERR(J) = FACMED*DEG7
              ELSE IF (R >= CRIVAL) THEN
                  INFOLD(5) = 0
                  RGNERR(J) = FACMED*DEG7*R
              ELSE
                  INFOLD(5) = 1
                  RGNERR(J) = FACOPT* (R**3)*DEG7
              END IF
              RGNERR(J) = MAX(NOISE,RGNERR(J))
          END IF
          RGNERR(J) = AREA*RGNERR(J)/FOUR
          BASVAL(J) = AREA*BASVAL(J)/FOUR
      END DO
      RETURN
      END SUBROUTINE Rule_C2a

END Module CubatureRule_C2
