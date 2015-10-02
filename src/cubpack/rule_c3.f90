! This file is F-compatible, except for upper/lower case conventions.
!--------------------------------------------------------------------
Module CubatureRule_C3

USE Precision_Model, ONLY: stnd

Implicit NONE

PRIVATE

PUBLIC :: Rule_C3a

CONTAINS
      SUBROUTINE Rule_C3a(VER,INFOLD,VOLUME,NUMFUN,Integrand,BASVAL,RGNERR,NUM)
!
!***BEGIN PROLOGUE Rule_C3a
!***DATE WRITTEN   970430   (YYMMDD)
!***REVISION DATE  990528   (YYMMDD)  (F conversion)
!***REVISION DATE  990604   (YYMMDD)  (divisions removed)
!***REVISION DATE  010919   (YYMMDD)  (subdivision information changed)
!***AUTHOR
!          Erwin Goor & Ronald Cools
!***PURPOSE  To compute basic integration rule values and
!            corresponding error estimates.
!***DESCRIPTION Rule_C3a computes basic integration rule values
!            for a vector of integrands over a cube.
!            Rule_C3a also computes estimates for the errors by
!            using several null rule approximations.
!   ON ENTRY
!
!   VER    Real array of dimension (3,4).
!          The coordinates of the vertices of the cube.
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
!   INFOLD(4) contains useful information for future subdivisions.
!          This is a 2 digit number.
!          The least significant digit contains info for 2-division.
!          The most significant digit contains info for 2/4/8-division.
!
!***REFERENCES Terje O. Espelid,
!              On the construction of good fully symmetric integration rules
!              SIAM J. NUMER. ANAL., Vol 24, No 4, August 1987.
!***ROUTINES CALLED
!                   OrbC3_Sum,Integrand
!***END PROLOGUE Rule_C3a
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
!   K      Integer array of dimension (0:4) that contains the structure
!          parameters. K(I) = number of orbits of type I.
!   TYPE0  Real array of dimension (K(0)).
!          Contains the first homogeneous coordinate of the generators
!          of type 0
!   TYPE1  Real array of dimension (K(1)).
!          Contains the first homogeneous coordinate of the generators
!          of type 1
!   TYPE2  Real array of dimension (K(2)).
!          Contains the first homogeneous coordinate of the generators
!          of type 2
!   TYPE3  Real array of dimension (K(3)).
!          Contains the first homogeneous coordinate of the generators
!          of type 3
!   TYPE4  Real array of dimension (2,K(4)).
!          Contains the first two homogeneous coordinates of
!          the generators of type 4.
!   WEIGHT Real array of dimension (8,ORBITS).
!          The weights of the cubature formula and the null rules.
!          WEIGHT(1,1) ,..., WEIGHT(1,ORBITS) are the weights of the
!                cubature formula
!          WEIGHT(I,1) ,..., WEIGHT(I,ORBITS) for I > 1, are the weights
!                of the null rules
!
!
!   Global variables.
!
      INTEGER, INTENT(IN)  :: NUMFUN
      INTEGER, INTENT(OUT) ::  NUM
      INTEGER, DIMENSION(:), INTENT(IN OUT) :: INFOLD
      REAL(kind=stnd), DIMENSION(:,:), INTENT(IN) :: VER
      REAL(kind=stnd), DIMENSION(:), INTENT(OUT) :: BASVAL, RGNERR
      REAL(kind=stnd), INTENT(IN) :: VOLUME

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
  INTEGER, PARAMETER ::   ORBITS = 8
  REAL(kind=stnd), PARAMETER:: HALF=0.5_stnd,                      &
                               TRES=50*EPSILON(HALF),              &
                               EIGHT=8.0_stnd,                     &
                               CRIVAL=0.5_stnd,                    &
                               FACMED=16.0_stnd,                   &
                               FACOPT=FACMED/CRIVAL,               &
                               CUTOFF=0.0001_stnd,                 &
                               BOUND1=0.55_stnd,                   &
                               BOUND2=0.1_stnd          
  REAL(kind=stnd),DIMENSION(0:2),PARAMETER ::                      &
                            DFC = (/-1.888242615950863158_stnd,    &
                                     1.0_stnd,                     &      
                                    -0.055878692024568421_stnd/)
!
!  Cubature formula of degree 11 with 89 points
!
          INTEGER, DIMENSION(0:4), PARAMETER ::                    &
                      K = (/1,2,1,2,2/)    ! Rule structure parameters
!
!  Information for the generators
!
     INTEGER :: I
     REAL(kind=stnd), DIMENSION(1:2), PARAMETER ::                 &
             TYPE1 = (/ 0.18052075573249470822058794973169_stnd,   &
                        0.76366700531389881917491285487593_stnd/)     
     REAL(kind=stnd), DIMENSION(1:1), PARAMETER ::                 &
             TYPE2 = (/ 0.78830647250493547382701809890320_stnd/)    
     REAL(kind=stnd), DIMENSION(1:2), PARAMETER ::                 &
             TYPE3 = (/ 0.52144095618907883093466563244662_stnd,   &
                        0.82617221731829521680258049987187_stnd/)
     REAL(kind=stnd), DIMENSION(1:2,1:2), PARAMETER ::             &
              TYPE4 = RESHAPE( SOURCE= (/                          &
                         0.97242481190902569341735290339642_stnd,  &
                         0.42110799536126494561573004212340_stnd,  &
                         0.46744847718063717386066856549398_stnd,  &
                         0.95930890177116312068682678508629_stnd/),&
                        SHAPE=(/2,2/), ORDER=(/2,1/) )
!
!  Weights of the cubature formula
!
     REAL(kind=stnd), DIMENSION(1:ORBITS), PARAMETER ::            &
          W1 = (/                                                  &
                       -3.8953757160532611282600158759761_stnd,    &
                       0.81856660649734492124488620848997_stnd,    &
                       0.21531990811243713389570303707084_stnd,    &
                       0.13294125786605972133694036186057_stnd,    &
                       0.20860736817139167236507600668966_stnd,    &
                       0.58773830247770211247176280672459E-1_stnd, &
                       0.14162520571796920029532663140195E-1_stnd, &
                       0.67408810538559624656933076917621E-1_stnd /)
!
!  Weights of the null rule of degree 7
!
      REAL(kind=stnd), DIMENSION(1:ORBITS), PARAMETER ::           &
           W2 = (/                                                 &
                      -4.09716713830698487528731726354974_stnd,    &
                       0.74677969141515975595448607337805_stnd,    &
                      -0.46417363966526945568711797402935E-1_stnd, &
                       0.17727896857214604523364317327611E-1_stnd, &
                      -0.52470267108342327292624832196503E-1_stnd, &
                       0.11664449246716569286439872549707E-1_stnd, &
                      -0.53791657604525012052256857235844E-2_stnd, &
                       0.57418721863532854861331639295054E-2_stnd /)
!
!  Weights of null rule of degree 5
!
      REAL(kind=stnd), DIMENSION(1:ORBITS), PARAMETER ::           &
           W3 = (/                                                 &
                      -0.38418257346076321180942534161643E-1_stnd, &
                      -0.13809393819461658172313782493931_stnd,    &
                      -0.36966067391537395134033779422780_stnd,    &
                      -0.84000216250042784132468820210803_stnd,    &
                       0.98945444020489415899976989198423_stnd,    &
                      -0.41711084444475038186621253136889_stnd,    &
                       0.16359771673332304281428383670012_stnd,    &
                       0.19416157968042709878528265453063_stnd /)
!
!  Weights of first null rule of degree 5
!
      REAL(kind=stnd), DIMENSION(1:ORBITS), PARAMETER ::           &
           W4 = (/                                                 &
                      -0.14339024357997965865651041863131_stnd,    &
                      -0.66574999845489902738481549374176_stnd,    &
                       1.41839687939449181019382815244241_stnd,    &
                      -0.38081150117982085429476186080312E-1_stnd, &
                       0.34239475207312939941374219725925_stnd,    &
                      -0.23914717080898300739086398648709_stnd,    &
                       0.17286429544022414400537142491791_stnd,    &
                      -0.37042670755501427522315796603400_stnd/)
!
!  Weights of second null rule of degree 5
!
      REAL(kind=stnd), DIMENSION(1:ORBITS), PARAMETER ::           &
           W5 = (/                                                 &
                       0.22799738734760429615753528373102E-2_stnd, &
                       0.11905590372687096358500385350270E-1_stnd, &
                       0.33061061773693696525589383114487_stnd,    &
                      -0.85087326218125136673746262977113_stnd,    &
                      -0.36662650042420411259687111438947_stnd,    &
                       0.67323896068825610599714357264483_stnd,    &
                      -0.21824874793205523689197441637244_stnd,    &
                       0.45575717466252940526695071801422_stnd/)
!
!  Weights of null rule of degree 3
!
      REAL(kind=stnd), DIMENSION(1:ORBITS), PARAMETER ::           &
           W6 = (/                                                 &
                       0.90256803909235587533523794574584E-1_stnd, &
                       0.48637869392568312594620915328833_stnd,    &
                       0.22605864607138151215790621639873_stnd,    &
                      -0.39733053275758221554551132969977_stnd,    &
                      -0.70260525791346956325766210839408_stnd,    &
                      -0.99463956560911321939850485360746_stnd,    &
                       0.40396337255003509831702095024785_stnd,    &
                       0.17858013317413262800119803474017_stnd/)
!
!   Weights of null rule of degree 3
!
      REAL(kind=stnd), DIMENSION(1:ORBITS), PARAMETER ::           &
           W7 = (/                                                 &
                      -0.14921225605249988383647637370117_stnd,    &
                      -0.85591937094309731035699238944552_stnd,    &
                      -0.19099496525254529541694749463933_stnd,    &
                       0.19040522291618896750255618030982_stnd,    &
                      -0.27115051668249444743700654045145_stnd,    &
                      -0.35979905900549990374887346342547_stnd,    &
                      -0.30255473948150702142743236951006_stnd,    &
                       0.68561441463717546800811910057288_stnd/)
!
!   Weights of null rule of degree 1
!
      REAL(kind=stnd), DIMENSION(1:ORBITS), PARAMETER ::           &
           W8 = (/                                                 &
                      -0.08621714823499486708026259436068_stnd,    &
                      -0.50649036609789532778219912403455_stnd,    &
                      -0.32380295374003822370848233043689_stnd,    &
                      -0.20985459608863301976909856362672_stnd,    &
                      -0.32887361592661612467795805597359_stnd,    &
                       0.21614962128074085470205806229955_stnd,    &
                       0.73079630784536188048918935662702_stnd,    &
                      -0.37712863378314543994499210520609_stnd/)
!
     REAL(kind=stnd), DIMENSION(1:8,1:ORBITS), PARAMETER ::        &
          WEIGHT =  RESHAPE( SOURCE= (/ W1,W2,W3,W4,W5,W6,W7,W8/), &
                          SHAPE=(/8,ORBITS/), ORDER=(/2,1/) )
!
!   Local variables.
!
      INTEGER :: J,NR,P,REGTYPE,NUMBER
      INTEGER, DIMENSION(1:3) :: DIFFINDEX
      REAL(kind=stnd):: NOISE,DEG7,DEG5,DEG3,DEG1,R3,R2,R1,R,      &
                Z1,Z2
      REAL(kind=stnd), DIMENSION(NUMFUN,7) :: NullRule
      REAL(kind=stnd), DIMENSION(NUMFUN) :: SUMVAL
      REAL(kind=stnd), DIMENSION(1:3) :: DIFF
      REAL(kind=stnd), DIMENSION(3,24) :: X
!
!***FIRST EXECUTABLE STATEMENT Rule_C3a
!

                                          
!  The number of points used by the cubature formula is
!     NUM    = 1*K(0) + 6*K(1) + 12*K(2) + 8*K(3) +24*K(4) = 89
          NUM = 89
!
!  Initialise BASVAL and NullRule
!
          BASVAL = 0 
          NullRule = 0
!
!  Compute contributions from orbits with 1, 6, 12, 8 and 24  points
!
          P = 1
          DO REGTYPE = 0,4
                DO NR = 1,K(REGTYPE)
                      SELECT CASE (REGTYPE)
                      CASE (0) !       Generator ( 0, 0, 0 )
                              NUMBER = 1
                              X(:,1) = (VER(:,2) + VER(:,3))*HALF +                 & 
                                             ( VER(:,4) - VER(:,1) )*HALF
                      CASE (1) !       Generator ( z1 , 0, 0 )
                              Z1 = TYPE1(NR)*HALF
                              NUMBER = 6
                              X(:,1) = -VER(:,1)*(z1+HALF) + VER(:,2)*HALF +        &
                                   VER(:,3)*HALF + VER(:,4)*(Z1+HALF)         
                              X(:,2) = VER(:,1)*(z1-HALF) + VER(:,2)*HALF +         &
                                   VER(:,3)*HALF + VER(:,4)*(-Z1+HALF)                
                              X(:,3) = -VER(:,1)*(z1+HALF) + VER(:,2)*HALF +        &
                                   VER(:,4)*HALF + VER(:,3)*(Z1+HALF)
                              X(:,4) = VER(:,1)*(z1-HALF) + VER(:,2)*HALF +         &
                                   VER(:,4)*HALF + VER(:,3)*(-Z1+HALF)
                              X(:,5) = -VER(:,1)*(z1+HALF) + VER(:,4)*HALF +        &
                                   VER(:,3)*HALF + VER(:,2)*(Z1+HALF)
                              X(:,6) = VER(:,1)*(z1-HALF) + VER(:,4)*HALF +         &
                                   VER(:,3)*HALF + VER(:,2)*(-Z1+HALF)
                      CASE (2) !       Generator ( z1 , z1, 0 )
                              Z1 = TYPE2(NR)*HALF
                              NUMBER = 12
                              X(:,1) = -VER(:,1)*(2*z1+HALF) + VER(:,4)*HALF +      &
                                   (VER(:,2) + VER(:,3))*(Z1+HALF) 
                              X(:,2) = VER(:,1)*(2*z1-HALF) + VER(:,4)*HALF +       &
                                   (VER(:,2) + VER(:,3))*(-Z1+HALF) 
                              X(:,3) = -VER(:,1)*(HALF) + VER(:,4)*HALF +           &
                                   VER(:,2)*(-z1+HALF) + VER(:,3)*(z1+HALF) 
                              X(:,4) = -VER(:,1)*(HALF) + VER(:,4)*HALF +           &
                                   VER(:,3)*(-z1+HALF) + VER(:,2)*(z1+HALF) 
                              X(:,5) = -VER(:,1)*(2*z1+HALF) + VER(:,3)*HALF +      &
                                   (VER(:,2) + VER(:,4))*(Z1+HALF) 
                              X(:,6) = VER(:,1)*(2*z1-HALF) + VER(:,3)*HALF +       &
                                   (VER(:,2) + VER(:,4))*(-Z1+HALF) 
                              X(:,7) = -VER(:,1)*(HALF) + VER(:,3)*HALF +           &
                                   VER(:,2)*(-z1+HALF) + VER(:,4)*(z1+HALF) 
                              X(:,8) = -VER(:,1)*(HALF) + VER(:,3)*HALF +           &
                                   VER(:,4)*(-z1+HALF) + VER(:,2)*(z1+HALF) 
                              X(:,9) = -VER(:,1)*(2*z1+HALF) + VER(:,2)*HALF +      &
                                   (VER(:,4) + VER(:,3))*(Z1+HALF) 
                              X(:,10) = VER(:,1)*(2*z1-HALF) + VER(:,2)*HALF +      &
                                    (VER(:,4) + VER(:,3))*(-Z1+HALF) 
                              X(:,11) = -VER(:,1)*(HALF) + VER(:,2)*HALF +          &
                                    VER(:,4)*(-z1+HALF) + VER(:,3)*(z1+HALF) 
                              X(:,12) = -VER(:,1)*(HALF) + VER(:,2)*HALF +          &
                                    VER(:,3)*(-z1+HALF) + VER(:,4)*(z1+HALF) 
                      CASE (3) !       Generator ( z1 , z1, z1 )          
                              Z1 = TYPE3(NR)*HALF
                              NUMBER = 8
                              X(:,1) = -VER(:,1)*(3*z1 + HALF) + VER(:,2)*(HALF + z1) + &
                                   VER(:,3)*(HALF + z1) + VER(:,4)*(HALF + z1)
                              X(:,2) = -VER(:,1)*(z1 + HALF) + VER(:,2)*(HALF - z1) +   &
                                   VER(:,3)*(HALF + z1) + VER(:,4)*(HALF + z1)
                              X(:,3) = -VER(:,1)*(z1 + HALF) + VER(:,2)*(HALF + z1) +   &
                                   VER(:,3)*(HALF - z1) + VER(:,4)*(HALF + z1)
                              X(:,4) = -VER(:,1)*(z1 + HALF) + VER(:,2)*(HALF + z1) +   &
                                   VER(:,3)*(HALF + z1) + VER(:,4)*(HALF - z1)
                              X(:,5) = -VER(:,1)*(-z1 + HALF) + VER(:,2)*(HALF - z1) +  &
                                   VER(:,3)*(HALF - z1) + VER(:,4)*(HALF + z1)
                              X(:,6) = -VER(:,1)*(-z1 + HALF) + VER(:,2)*(HALF - z1) +  &
                                   VER(:,3)*(HALF + z1) + VER(:,4)*(HALF - z1)
                              X(:,7) = -VER(:,1)*(-z1 + HALF) + VER(:,2)*(HALF + z1) +  &
                                   VER(:,3)*(HALF - z1) + VER(:,4)*(HALF - z1)
                              X(:,8) = VER(:,1)*(3*z1 - HALF) + VER(:,2)*(HALF - z1) +  &
                                   VER(:,3)*(HALF - z1) + VER(:,4)*(HALF - z1)
                      CASE (4) !       Generator ( z1 , z1, z2 )   
                                                Z1 = TYPE4(1,NR)*HALF
                              Z2 = TYPE4(2,NR)*HALF
                              NUMBER = 24
                              X(:,1) = -VER(:,1)*(2*z1+z2+HALF) + VER(:,4)*(HALF + z2) +   &
                                   VER(:,2)*(HALF + z1)+ VER(:,3)*(HALF + z1) 
                              X(:,2) = -VER(:,1)*(2*z1-z2+HALF) + VER(:,4)*(HALF - z2) +   &
                                   VER(:,2)*(HALF + z1)+ VER(:,3)*(HALF + z1) 
                              X(:,3) = -VER(:,1)*(z2+HALF) + VER(:,4)*(HALF + z2) +        &
                                   VER(:,2)*(HALF - z1)+ VER(:,3)*(HALF + z1) 
                              X(:,4) = -VER(:,1)*(z2+HALF) + VER(:,4)*(HALF + z2) +        &
                                   VER(:,2)*(HALF + z1)+ VER(:,3)*(HALF - z1) 
                              X(:,5) = -VER(:,1)*(-2*z1+z2+HALF) + VER(:,4)*(HALF + z2) +  &
                                   VER(:,2)*(HALF - z1)+ VER(:,3)*(HALF - z1) 
                              X(:,6) = -VER(:,1)*(-z2+HALF) + VER(:,4)*(HALF - z2) +       &
                                   VER(:,2)*(HALF + z1)+ VER(:,3)*(HALF - z1) 
                              X(:,7) = -VER(:,1)*(-z2+HALF) + VER(:,4)*(HALF - z2) +       &
                                   VER(:,2)*(HALF - z1)+ VER(:,3)*(HALF + z1) 
                              X(:,8) = -VER(:,1)*(-2*z1-z2+HALF) + VER(:,4)*(HALF - z2) +  &
                                   VER(:,2)*(HALF - z1)+ VER(:,3)*(HALF - z1)            
                              X(:,9) = -VER(:,1)*(2*z1+z2+HALF) + VER(:,2)*(HALF + z2) +   &
                                   VER(:,4)*(HALF + z1)+ VER(:,3)*(HALF + z1) 
                              X(:,10) = -VER(:,1)*(2*z1-z2+HALF) + VER(:,2)*(HALF - z2) +  &
                                    VER(:,4)*(HALF + z1)+ VER(:,3)*(HALF + z1) 
                              X(:,11) = -VER(:,1)*(z2+HALF) + VER(:,2)*(HALF + z2) +       &
                                    VER(:,4)*(HALF - z1)+ VER(:,3)*(HALF + z1) 
                              X(:,12) = -VER(:,1)*(z2+HALF) + VER(:,2)*(HALF + z2) +       &
                                    VER(:,4)*(HALF + z1)+ VER(:,3)*(HALF - z1) 
                              X(:,13) = -VER(:,1)*(-2*z1+z2+HALF) + VER(:,2)*(HALF + z2) + &
                                    VER(:,4)*(HALF - z1)+ VER(:,3)*(HALF - z1) 
                              X(:,14) = -VER(:,1)*(-z2+HALF) + VER(:,2)*(HALF - z2) +      &
                                    VER(:,4)*(HALF + z1)+ VER(:,3)*(HALF - z1)           
                              X(:,15) = -VER(:,1)*(-z2+HALF) + VER(:,2)*(HALF - z2) +      &
                                    VER(:,4)*(HALF - z1)+ VER(:,3)*(HALF + z1) 
                              X(:,16) = -VER(:,1)*(-2*z1-z2+HALF) + VER(:,2)*(HALF - z2) + &
                                    VER(:,4)*(HALF - z1)+ VER(:,3)*(HALF - z1) 
                              X(:,17) = -VER(:,1)*(2*z1+z2+HALF) + VER(:,3)*(HALF + z2) +  &
                                    VER(:,2)*(HALF + z1)+ VER(:,4)*(HALF + z1) 
                              X(:,18) = -VER(:,1)*(2*z1-z2+HALF) + VER(:,3)*(HALF - z2) +  &
                                    VER(:,2)*(HALF + z1)+ VER(:,4)*(HALF + z1) 
                              X(:,19) = -VER(:,1)*(z2+HALF) + VER(:,3)*(HALF + z2) +       &
                                    VER(:,2)*(HALF - z1)+ VER(:,4)*(HALF + z1) 
                              X(:,20) = -VER(:,1)*(z2+HALF) + VER(:,3)*(HALF + z2) +       &
                                    VER(:,2)*(HALF + z1)+ VER(:,4)*(HALF - z1) 
                              X(:,21) = -VER(:,1)*(-2*z1+z2+HALF) + VER(:,3)*(HALF + z2) + &
                                    VER(:,2)*(HALF - z1)+ VER(:,4)*(HALF - z1) 
                              X(:,22) = -VER(:,1)*(-z2+HALF) + VER(:,3)*(HALF - z2) +      &
                                    VER(:,2)*(HALF + z1)+ VER(:,4)*(HALF - z1) 
                              X(:,23) = -VER(:,1)*(-z2+HALF) + VER(:,3)*(HALF - z2) +      &
                                    VER(:,2)*(HALF - z1)+ VER(:,4)*(HALF + z1) 
                              X(:,24) = -VER(:,1)*(-2*z1-z2+HALF) + VER(:,3)*(HALF - z2) + &
                                    VER(:,2)*(HALF - z1)+ VER(:,4)*(HALF - z1) 
                      END SELECT
                      SUMVAL = Integrand(NUMFUN,X(:,1))
                      SELECT CASE (REGTYPE)
                      CASE (0)
                           DIFF(2) = SUMVAL(1)*DFC(0)
                           DIFF(1) = DIFF(2)
                           DIFF(3) = DIFF(2)
                      CASE (1)
                           DIFF(3) = DIFF(3) + SUMVAL(1)*DFC(NR)
                      END SELECT
                      DO J = 2,NUMBER
                           RGNERR = Integrand(NUMFUN,X(:,J))
                           IF ( REGTYPE == 1 ) THEN
                                IF ( J <= 2 ) THEN
                                     DIFF(3) = DIFF(3) + RGNERR(1)*DFC(NR)
                                ELSE IF ( (J == 3) .OR. (J == 4) ) THEN
                                     DIFF(2) = DIFF(2) + RGNERR(1)*DFC(NR)
                                ELSE
                                     DIFF(1) = DIFF(1) + RGNERR(1)*DFC(NR)
                                END IF
                           END IF
                           SUMVAL(1:NUMFUN) = SUMVAL(1:NUMFUN) + RGNERR(1:NUMFUN)
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
          DIFF = ABS(DIFF)
          ! Sort the fourth differences   
          DIFFINDEX = (/ 1,2,3 /)
          DO I = 1,2
               DO J = 1,(3-I)
                    IF ( DIFF(J+1) > DIFF(J) ) THEN
                         R = DIFF(J)  ! dummy
                         DIFF(J) = DIFF(J+1)
                         DIFF(J+1) = R
                         P = DIFFINDEX(J)  ! dummy
                         DIFFINDEX(J) = DIFFINDEX(J+1)
                         DIFFINDEX(J+1) = P
                    END IF
               END DO
          END DO
          if ( diff(1) < cutoff ) then   ! recommend uniform division
             !!  infold(4) = 0
             infold(4) = 90+diffindex(1)
          else
             if ( diff(2) > bound1*diff(1) ) then
                if ( diff(3) > bound2*diff(2) ) then
                   ! recommend uniform division
                   !! infold(4) = 0
                   infold(4) = 90+diffindex(1)
                else
                   ! recommend division in 4 
                   select case (diffindex(3))
                   case(1)
                          infold(4) = -3
                   case(2)
                          infold(4) = -1  
                   case(3)
                          infold(4) = -2
                   end select
                   infold(4) = 10*infold(4) - diffindex(1)
                end if
             else
                ! recommend division in 2
                infold(4) = diffindex(1)
             end if
          end if
!
!  Compute error estimates
!
          DO J = 1,NUMFUN
                NOISE = ABS(BASVAL(J))*TRES
                DEG7 = ABS( NullRule(J,1) )
                IF (DEG7 <= NOISE) THEN
                      RGNERR(J) = NOISE
                ELSE
                         DEG5 = SQRT(NullRule(J,2)**2+NullRule(J,3)**2)
                      DEG3 = SQRT(NullRule(J,4)**2+NullRule(J,5)**2)
                      DEG1 = SQRT(NullRule(J,6)**2+NullRule(J,7)**2)
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
                           RGNERR(J) = FACMED*R*DEG7
                      ELSE  
                           INFOLD(5) = 1
                           RGNERR(J) = FACOPT*(R**2)*DEG7
                      END IF
                      RGNERR(J) = MAX(NOISE,RGNERR(J))
                END IF
                RGNERR(J) = VOLUME*RGNERR(J)/EIGHT
                BASVAL(J) = VOLUME*BASVAL(J)/EIGHT
          END DO
          RETURN
          END SUBROUTINE Rule_C3a

END Module CubatureRule_C3
