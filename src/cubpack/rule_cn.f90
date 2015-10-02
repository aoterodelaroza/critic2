! This file is F-compatible, except for upper/lower case conventions.
!--------------------------------------------------------------------

Module CubatureRule_Cn

USE Precision_Model, ONLY: stnd

Implicit NONE

PRIVATE

PUBLIC  :: Rule_Cn
PRIVATE :: Rule_Deg7, Rule_Deg9, SymCub_Sum

CONTAINS
  SUBROUTINE Rule_Cn(KEY,N,VERTEX,INFOLD,VOLUME,NF,Integrand,BASVAL,RGNERR,NUM)
!
!***BEGIN PROLOGUE Rule_Cn
!***DATE WRITTEN   990701   (YYMMDD)
!***REVISION DATE  000814   (Handling of KEY now similar to Tn; changed by RC)
!***AUTHOR
!
!            Alan Genz
!            Department of Mathematics
!            Washington State University
!            Pullman, WA 99164-3113, USA
!
!
!***PURPOSE  To compute basic integration rule values and
!            corresponding error estimates.
!***DESCRIPTION Rule_Cn7computes basic integration rule values
!            for a vector of integrands over a dim hyperrectangle
!            Rule_Cn7 also computes estimates for the errors by
!            using several null rule approximations.
!            We use a degree 7 integration rule,
!            two degree 5 null rules, one degree 3 null rule and one
!            degree 1 null rule for the hypercube.
!            RESTRICTION : this routine will only give correct results
!            for dim > 2.
!   ON ENTRY
!
!   KEY       Integer.
!             If Key > 2 and Key < 5 then a rule of degree 2*Key + 1
!             is used; otherwise a default rule of degree 7 is used.
!   N         Integer, dimension of the integration problem
!   VERTEX    Real array of dimension (N,0:N).
!             The coordinates of the vertices of the parallelepiped.
!              vertex i -> ( vertex(1,i),vertex(2,i),...,vertex(N,i) )
!   NF        Integer, number of components of the vector integrand.
!   Integrand Real vector function of length NF for computing components of
!              the integrand at X.
!              It must have parameters ( NF, X ); see interface below
!               Input parameters:
!                 X      Real array of length N, the evaluation point.
!                 NF     Integer number of components of Integrand.
!
!   ON RETURN
!
!   BASVAL Real array of dimension NF.
!          The values for the basic rule for each component
!          of the integrand.
!   RGNERR Real array of dimension NF.
!          The error estimates for each component of the integrand.
!   NUM    Integer, number of function evaluations used.
!
!***REFERENCES A. Genz and A. Malik,
!             "An Imbedded Family of Fully Symmetric Numerical
!              Integration Rules",
!              SIAM J Numer. Anal. 20 (1983), pp. 580-588.
!***ROUTINES CALLED
!                   Integrand, SymCub_Sum, Rule_Deg7, Rule_Deg9
!***END PROLOGUE Rule_Cn
!
   INTEGER,                          INTENT(IN)    :: KEY, N, NF
   INTEGER,                          INTENT(OUT)   :: NUM
   INTEGER, DIMENSION(:),            INTENT(INOUT) :: INFOLD
   REAL(KIND=STND), DIMENSION(:,0:), INTENT(IN)    :: VERTEX
   REAL(KIND=STND), DIMENSION(:),    INTENT(OUT)   :: RGNERR, BASVAL
   REAL(KIND=STND),                  INTENT(IN)    :: VOLUME
   INTERFACE 
      FUNCTION INTEGRAND(NF,X) RESULT(VALUE)
         USE PRECISION_MODEL
         INTEGER,                       INTENT(IN) :: NF
         REAL(KIND=STND), DIMENSION(:), INTENT(IN) :: X
         REAL(KIND=STND), DIMENSION(NF)            :: VALUE
      END FUNCTION INTEGRAND
   END INTERFACE
!
   INTEGER, PARAMETER                        :: MXW = 9, MXG = 4, RLS = 5
   REAL(KIND=STND), DIMENSION(MXW,RLS), SAVE :: W
   REAL(KIND=STND), DIMENSION(MXG,MXW), SAVE :: G
   INTEGER,                             SAVE :: OLDKEY = -1, OLDN = 0
   INTEGER,                             SAVE :: WTS, NUMR
!
   INTEGER                                   :: I, DIVAXN
   REAL(KIND=STND), PARAMETER                :: ONE = 1
   REAL(KIND=STND), PARAMETER                :: SMALL = 100*EPSILON(ONE) 
   REAL(KIND=STND), DIMENSION(N)             :: CENTER, DFS, GTEMP
   REAL(KIND=STND), DIMENSION(N,N)           :: VERDIF
   REAL(KIND=STND), DIMENSION(NF,RLS)        :: RULE
   REAL(KIND=STND), DIMENSION(NF,3)          :: TEMP
   REAL(KIND=STND), DIMENSION(NF)            :: FRTHDF
   REAL(KIND=STND)                           :: RATIO
!
   IF ( KEY /= OLDKEY .OR. OLDKEY == -1 .OR. N /= OLDN ) THEN
      OLDKEY = KEY
      OLDN = N
      SELECT CASE ( KEY )
      CASE (4)
         CALL Rule_Deg9( N, W, G, WTS, NUMR )
      CASE DEFAULT
         CALL Rule_Deg7( N, W, G, WTS, NUMR )
      END SELECT
   END IF 
   NUM = NUMR
   VERDIF = ( VERTEX(:,1:N) - SPREAD( VERTEX(:,0), 2, N ) )/2
   CENTER = VERTEX(:,0) + SUM( VERDIF, 2 )
   DIVAXN = SUM( MAXLOC( SUM( ABS( VERDIF ), 1 ) ) )
   TEMP(:,1) = INTEGRAND( NF, CENTER )
   RULE = MATMUL( TEMP(:,1:1), W(1:1,1:RLS)  ) 
   RATIO = ( G(1,3)/G(1,2) )**2
   DO I = 1, N
      TEMP(:,2) = INTEGRAND( NF, CENTER - G(1,2)*VERDIF(:,I) )               &
                + INTEGRAND( NF, CENTER + G(1,2)*VERDIF(:,I) )
      TEMP(:,3) = INTEGRAND( NF, CENTER - G(1,3)*VERDIF(:,I) )               &
                + INTEGRAND( NF, CENTER + G(1,3)*VERDIF(:,I) )
      RULE = RULE + MATMUL( TEMP(:,2:3), W(2:3,1:RLS) )
      FRTHDF = ABS( 2*(1-RATIO)*TEMP(:,1) + RATIO*TEMP(:,2) - TEMP(:,3) )/4
      DFS(I) = SUM( FRTHDF, MASK = ABS(TEMP(:,1)) + FRTHDF > ABS(TEMP(:,1)) )
   END DO
   IF ( MAXVAL( DFS ) > 0 ) THEN
      DIVAXN = SUM( MAXLOC( DFS ) )
   END IF
   INFOLD(4) = DIVAXN
!
!    Finish computing the rule values.
!
   DO I = 4, WTS
      GTEMP( 1 : MIN(N,MXG-1) ) = G( 1 : MIN(N,MXG-1) , I )
      IF ( N >= MXG ) THEN
         GTEMP(MXG:N) = G(MXG,I)
      END IF
      TEMP(:,1) = SymCub_Sum( N, VERTEX, GTEMP, NF, INTEGRAND )
      RULE = RULE + MATMUL( TEMP(:,1:1), W(I:I,1:RLS) )
   END DO
!
!    Compute errors.
!
      RULE(:,2:5) = ABS( RULE(:,2:5) )
      RULE(:,3) = MAX( RULE(:,2), RULE(:,3) )
      RULE(:,2) = ABS( RULE(:,1) )
      DO I = 3, 5
         WHERE ( RULE(:,2) + RULE(:,I)/NUM <= RULE(:,2) )        
            RULE(:,I) = 0
         END WHERE
      END DO
      RATIO = 5 + 8*KEY
      WHERE ( RATIO*RULE(:,3) <= RULE(:,4) .AND. RATIO*RULE(:,4) <= RULE(:,5) )
          RGNERR = VOLUME*RULE(:,3)/2
      ELSEWHERE 
          RGNERR = VOLUME*MAXVAL( RULE(:,3:5), 2 )
      END WHERE
      BASVAL = VOLUME*RULE(:,1)
!
END SUBROUTINE Rule_Cn
!
SUBROUTINE Rule_Deg7( N, W, G, WTSR, NUMR )
!
!***BEGIN PROLOGUE Rule_Deg7
!***KEYWORDS basic integration rule, degree 7
!***PURPOSE  To initialize a degree 7 basic rule, and null rules.
!***AUTHOR   
!
!            Alan Genz
!            Department of Mathematics
!            Washington State University
!            Pullman, WA 99164-3113, USA
!
!***DATE WRITTEN  990701   (YYMMDD)
!***DESCRIPTION  Rule_Deg7 initializes a degree 7 integration rule,
!            two degree 5 null rules, one degree 3 null rule and one
!            degree 1 null rule for the hypercube [-1,1]**N.
!
!   ON ENTRY
!
!   N      Integer, number of variables.
!
!   ON RETURN
!   W      Real array of dimension (WTS,5).
!          The weights for the basic and null rules.
!          W(1,1),...,W(WTS,1) are weights for the basic rule.
!          W(I,1),...,W(WTS,I), for I > 1 are null rule weights.
!   G      Real array of dimension (N, WTS).
!          The fully symmetric sum generators for the rules.
!          G(1, J), ..., G(N, J) are the are the generators for the
!          points associated with the Jth weights.
!   WTSR   Integer, WTS
!   NUMR   Integer, number of points for Rule_Deg7
!
!***REFERENCES A. Genz and A. Malik,
!             "An Imbedded Family of Fully Symmetric Numerical
!              Integration Rules",
!              SIAM J Numer. Anal. 20 (1983), pp. 580-588.
!***ROUTINES CALLED-NONE
!***END PROLOGUE Rule_Deg7
!
!   Global variables
!
!   WTS Integer, PARAMETER.
!          The number of weights in each of the rules : 6
!
      INTEGER,                         INTENT(IN)  :: N
      REAL(KIND=STND), DIMENSION(:,:), INTENT(OUT) :: W, G
      INTEGER,                         INTENT(OUT) :: WTSR, NUMR
!
!   Constant, the size of RULPTS
!
      INTEGER, PARAMETER              :: WTS = 6  
!
!   Local Variables
!
      INTEGER                         :: K
      REAL(KIND=STND)                 :: TEMP, LAM0, LAM1, LAM2, LAMP, TWONDM
      REAL(KIND=STND), DIMENSION(WTS) :: RULPTS
      REAL(KIND=STND), DIMENSION(3)   :: ALPHA
!
!     Initialize generators, weights and RULPTS
!
      G = 0
      W = 0
      TWONDM = 2**N
      RULPTS(1) = 1
      RULPTS(2:WTS-2) = 2*N
      RULPTS(WTS-1) = 2*N*(N-1)
      RULPTS(WTS) = TWONDM
!
!     Compute squared generator parameters
!
      LAM0 = 0.4707_STND
      LAMP = 0.5625_STND
      LAM1 = 4/( 15 - 5/LAM0 )
      TEMP = (1 - LAM1/LAM0 )/27
      LAM2 = ( 5 - 7*LAM1 - 35*TEMP )/( 7 - 35*LAM1/3 - 35*TEMP/LAM0 )
!
!     Compute degree 7 rule weights
!
      W(6,1) = 1/(3*LAM0)**3/TWONDM
      W(5,1) = ( 1 - 5*LAM0/3 )/( 60*(LAM1-LAM0)*LAM1**2 )
      W(3,1) = ( 1 - 5*LAM2/3 - 5*TWONDM*W(6,1)*LAM0*(LAM0-LAM2) )          &
              /( 10*LAM1*(LAM1-LAM2) ) - 2*(N-1)*W(5,1)
      W(2,1) = ( 1 - 5*LAM1/3 - 5*TWONDM*W(6,1)*LAM0*(LAM0-LAM1) )          &
              /( 10*LAM2*(LAM2-LAM1) )
!
!     Compute weights for 2 degree 5, 1 degree 3 and 1 degree 1 rules
!
      W(6,2) = 1/( 36*LAM0**3 )/TWONDM
      W(5,2) = ( 1 - 9*TWONDM*W(6,2)*LAM0**2 )/( 36*LAM1**2 )
      W(3,2) = ( 1 - 5*LAM2/3 - 5*TWONDM*W(6,2)*LAM0*(LAM0-LAM2) )          &
              /( 10*LAM1* (LAM1-LAM2) ) - 2* (N-1)*W(5,2)
      W(2,2) = ( 1 - 5*LAM1/3 - 5*TWONDM*W(6,2)*LAM0*(LAM0-LAM1) )          &
              /( 10*LAM2* (LAM2-LAM1) )
      W(6,3) = 5/(108*LAM0**3)/TWONDM
      W(5,3) = ( 1 - 9*TWONDM*W(6,3)*LAM0**2 )/( 36*LAM1**2 )
      W(4,3) = ( 1 - 5*LAM1/3 - 5*TWONDM*W(6,3)*LAM0*(LAM0-LAM1) )          &
              /( 10*LAMP*(LAMP-LAM1) )
      W(3,3) = ( 1 - 5*LAMP/3 - 5*TWONDM*W(6,3)*LAM0* (LAM0-LAMP) )         &
              /( 10*LAM1*(LAM1-LAMP) ) - 2*(N-1)*W(5,3)
      W(6,4) = 1/( 54*LAM0**3 )/TWONDM
      W(5,4) = ( 1 - 18*TWONDM*W(6,4)*LAM0**2 )/( 72*LAM1**2 )
      W(3,4) = ( 1 - 10*LAM2/3 - 10*TWONDM*W(6,4)*LAM0*(LAM0-LAM2) )        &
              /( 20*LAM1*(LAM1-LAM2) ) - 2*(N-1)*W(5,4)
      W(2,4) = ( 1 - 10*LAM1/3 - 10*TWONDM*W(6,4)*LAM0*(LAM0-LAM1) )        &
              /( 20*LAM2*(LAM2-LAM1) )
!
!     Set generator values
!
      LAM0 = SQRT(LAM0)
      LAM1 = SQRT(LAM1)
      LAM2 = SQRT(LAM2)
      LAMP = SQRT(LAMP)
      G(  :,6) = LAM0
      G(1:2,5) = LAM1
      G(  1,2) = LAM2
      G(  1,3) = LAM1
      G(  1,4) = LAMP
!
!     Compute constant weight values.
!
      W(1,1:5) = 1 - MATMUL( RULPTS(2:WTS), W(2:WTS,1:5) ) 
!
!     Compute final weight values; null rule weights are computed as 
!     differences between weights from highest degree and lower degree rules.
!
      W(1:WTS,2:5) = W(1:WTS,2:5) - SPREAD( W(1:WTS,1), 2, 4 )
!
!        Orthogonalize and normalize null rules.
!
      TEMP = SUM( RULPTS*W(1:WTS,1)*W(1:WTS,1) )
      W(1:WTS,2) = W(1:WTS,2)*SQRT( TEMP/SUM( RULPTS*W(1:WTS,2)*W(1:WTS,2) ) )
      DO K = 3, 5
         ALPHA(1:K-2) = -MATMUL( TRANSPOSE(W(1:WTS,2:K-1)), RULPTS*W(1:WTS,K) )
         W(1:WTS,K) = W(1:WTS,K) + MATMUL( W(1:WTS,2:K-1), ALPHA(1:K-2) )/TEMP
         W(1:WTS,K) = W(1:WTS,K)*SQRT(TEMP/SUM(RULPTS*W(1:WTS,K)*W(1:WTS,K)))
      END DO
      WTSR = WTS
      NUMR = SUM( RULPTS )
!
END SUBROUTINE Rule_Deg7
!
SUBROUTINE Rule_Deg9( N, W, G, WTSR, NUMR )
!***BEGIN PROLOGUE Rule_Deg9
!***KEYWORDS basic integration rule, degree 9
!***PURPOSE  To initialize a degree 9 basic rule and null rules.
!***AUTHOR   
!
!            Alan Genz
!            Department of Mathematics
!            Washington State University
!            Pullman, WA 99164-3113, USA
!
!***DATE WRITTEN  990701   (YYMMDD)
!***DESCRIPTION  Rule_Deg9 initializes a degree 9 integration rule,
!            two degree 7 null rules, one degree 5 null rule and one
!            degree 3 null rule for the hypercube [-1,1]**N.
!
!   ON ENTRY
!
!   N      Integer, number of variables.
!
!   ON RETURN
!   W      Real array of dimension (WTS,5).
!          The weights for the basic and null rules.
!          W(1,1),...,W(WTS,1) are weights for the basic rule.
!          W(I,1),...,W(WTS,I), for I > 1 are null rule weights.
!   G      Real array of dimension (N, WTS).
!          The fully symmetric sum generators for the rules.
!          G(1, J), ..., G(N, J) are the are the generators for the
!          points associated with the Jth weights.
!   WTSR   Integer, WTS
!   NUMR   Integer, number of points for Rule_Deg7
!
!***REFERENCES A. Genz and A. Malik,
!             "An Imbedded Family of Fully Symmetric Numerical
!              Integration Rules",
!              SIAM J Numer. Anal. 20 (1983), pp. 580-588.
!***ROUTINES CALLED-NONE
!***END PROLOGUE Rule_Deg9
!
!   Global variables
!
      INTEGER,                         INTENT(IN)  :: N
      REAL(KIND=STND), DIMENSION(:,:), INTENT(OUT) :: W, G
      INTEGER,                         INTENT(OUT) :: WTSR, NUMR
!
      INTEGER, PARAMETER              :: WTS = 9
!
!   Local Variables
!
      INTEGER                         :: K
      REAL(KIND=STND), DIMENSION(WTS) :: RULPTS
      REAL(KIND=STND), DIMENSION(3)   :: ALPHA
      REAL(KIND=STND)                 :: TEMP, LAM0,LAM1,LAM2,LAM3,LAMP, TWONDM
!
!***FIRST EXECUTABLE STATEMENT Rule_Deg9
!
!     Initialize generators, weights and RULPTS
!
      G = 0
      W = 0
      TWONDM = 2**N
      RULPTS(1) = 1
      RULPTS(2:5) = 2*N
      RULPTS(6) = 2*N*(N-1)
      RULPTS(7) = 4*N*(N-1)
      RULPTS(8) = (4*N*(N-1)*(N-2))/3
      RULPTS(9) = TWONDM
!
!     Compute squared generator parameters
!
      LAM0 = 0.4707_STND
      LAMP = 0.0625_STND
      LAM1 = 4/( 15 - 5/LAM0 )
      TEMP = ( 1 - LAM1/LAM0 )/27
      LAM2 = ( 5 - 7*LAM1 - 35*TEMP )/( 7 - 35*LAM1/3 - 35*TEMP/LAM0 )
      TEMP = TEMP*( 1 - LAM2/LAM0 )/3
      LAM3 = ( 7  - 9*(LAM2+LAM1)   + 63*LAM2*LAM1/5 - 63*TEMP )              &
            /( 9 - 63*(LAM2+LAM1)/5 + 21*LAM2*LAM1 - 63*TEMP/LAM0 )
!
!     Compute degree 9 rule weights
!
      W(9,1) = 1/(3*LAM0)**4/TWONDM
      W(8,1) = ( 1 - 1/(3*LAM0) )/(6*LAM1)**3
      W(7,1) = ( 1 - 7*(LAM0+LAM1)/5 + 7*LAM0*LAM1/3 )                        &
              /( 84*LAM1*LAM2*(LAM2-LAM0)*(LAM2-LAM1) )
      W(6,1) = ( 1 - 7*(LAM0+LAM2)/5 + 7*LAM0*LAM2/3 )                        &
              /( 84*LAM1*LAM1*(LAM1-LAM0)*(LAM1-LAM2) )                       &
              - W(7,1)*LAM2/LAM1 - 2*(N-2)*W(8,1)
      W(4,1) = ( 1 - 9*( (LAM0+LAM1+LAM2)/7                                   &
                       - (LAM0*LAM1+LAM0*LAM2+LAM1*LAM2)/5 )                  &
                   - 3*LAM0*LAM1*LAM2 )                                       &
              /( 18*LAM3*(LAM3-LAM0)*(LAM3-LAM1)*(LAM3-LAM2) )
      W(3,1) = ( 1 - 9*( (LAM0+LAM1+LAM3)/7                                   &
                       - (LAM0*LAM1+LAM0*LAM3+LAM1*LAM3)/5 )                  &
                    -3*LAM0*LAM1*LAM3 )                                       &
              /( 18*LAM2*(LAM2-LAM0)*(LAM2-LAM1)*(LAM2-LAM3) )                &
              - 2*(N-1)*W(7,1)
      W(2,1) = ( 1 - 9*( (LAM0+LAM2+LAM3)/7                                   &
                       - (LAM0*LAM2+LAM0*LAM3+LAM2*LAM3)/5 )                  &
                   - 3*LAM0*LAM2*LAM3 )                                       &
              /( 18*LAM1*(LAM1-LAM0)*(LAM1-LAM2)*(LAM1-LAM3) )                &
              - 2*(N-1)*( W(6,1) + W(7,1) + (N-2)*W(8,1) )
!
!     Compute weights for 2 degree 7, 1 degree 5 and 1 degree 3 rules
!
      W(9,2) = 1/( 108*LAM0**4 )/TWONDM
      W(8,2) = ( 1 - 27*TWONDM*W(9,2)*LAM0**3 )/(6*LAM1)**3
      W(7,2) = ( 1 - 5*LAM1/3 - 15*TWONDM*W(9,2)*LAM0**2*(LAM0-LAM1) )        &
              /( 60*LAM1*LAM2*(LAM2-LAM1) )
      W(6,2) = ( 1 - 9*( 8*LAM1*LAM2*W(7,2) + TWONDM*W(9,2)*LAM0**2 ) )       &
              /(36*LAM1*LAM1) - 2*W(8,2)*(N-2)
      W(4,2) = ( 1 - 7*( (LAM1+LAM2)/5 - LAM1*LAM2/3                          &
                       + TWONDM*W(9,2)*LAM0*(LAM0-LAM1)*(LAM0-LAM2) ) )       &
              /( 14*LAM3*(LAM3-LAM1)*(LAM3-LAM2) )
      W(3,2) = ( 1 - 7*( (LAM1+LAM3)/5 - LAM1*LAM3/3                          &
                       + TWONDM*W(9,2)*LAM0*(LAM0-LAM1)*(LAM0-LAM3) ) )       &
              /( 14*LAM2*(LAM2-LAM1)*(LAM2-LAM3) ) - 2*(N-1)*W(7,2)
      W(2,2) = ( 1 - 7*( (LAM2+LAM3)/5 - LAM2*LAM3/3                          &
                       + TWONDM*W(9,2)*LAM0*(LAM0-LAM2)*(LAM0-LAM3) ) )       &
              /( 14*LAM1*(LAM1-LAM2)*(LAM1-LAM3) )                            &
              - 2*(N-1)*( W(6,2) + W(7,2) + (N-2)*W(8,2) )
      W(9,3) = 5/( 324*LAM0**4 )/TWONDM
      W(8,3) = ( 1 - 27*TWONDM*W(9,3)*LAM0**3 )/(6*LAM1)**3
      W(7,3) = ( 1 - 5*LAM1/3 - 15*TWONDM*W(9,3)*LAM0**2*(LAM0-LAM1) )        &
              /( 60*LAM1*LAM2* (LAM2-LAM1) )
      W(6,3) = ( 1 - 9*( 8*LAM1*LAM2*W(7,3) + TWONDM*W(9,3)*LAM0**2 ) )       &
              /( 36*LAM1*LAM1) - 2*W(8,3)*(N-2)
      W(5,3) = ( 1 - 7*( (LAM1+LAM2)/5 - LAM1*LAM2/3                          &
                       + TWONDM*W(9,3)*LAM0*(LAM0-LAM1)*(LAM0-LAM2) ) )       &
              /( 14*LAMP* (LAMP-LAM1)*(LAMP-LAM2) )
      W(3,3) = ( 1 - 7*( (LAM1+LAMP)/5 - LAM1*LAMP/3                          &
                       + TWONDM*W(9,3)*LAM0*(LAM0-LAM1)*(LAM0-LAMP) ) )       &
              /( 14*LAM2*(LAM2-LAM1)*(LAM2-LAMP) ) - 2*(N-1)*W(7,3)
      W(2,3) = ( 1 - 7*( (LAM2+LAMP)/5 - LAM2*LAMP/3                          &
                       + TWONDM*W(9,3)*LAM0*(LAM0-LAM2)*(LAM0-LAMP) ) )       &
              /( 14*LAM1*(LAM1-LAM2)*(LAM1-LAMP) )                            &
              - 2*(N-1)*( W(6,3) + W(7,3) + (N-2)*W(8,3) )
      W(9,4) = 2/( 81*LAM0**4 )/TWONDM
      W(8,4) = ( 2 - 27*TWONDM*W(9,4)*LAM0**3)/(6*LAM1)**3
      W(7,4) = ( 2 - 15*LAM1/9-15*TWONDM*W(9,4)*LAM0*(LAM0-LAM1) )            &
              /( 60*LAM1*LAM2*(LAM2-LAM1) )
      W(6,4) = ( 1 - 9*( 8*LAM1*LAM2*W(7,4) + TWONDM*W(9,4)*LAM0**2 ) )       &
              /( 36*LAM1*LAM1 ) - 2*W(8,4)*(N-2)
      W(4,4) = ( 2 - 7*( (LAM1+LAM2)/5 - LAM1*LAM2/3                          &
                       + TWONDM*W(9,4)*LAM0*(LAM0-LAM1)*(LAM0-LAM2) ) )       &
              /( 14*LAM3*(LAM3-LAM1)*(LAM3-LAM2) )
      W(3,4) = ( 2 - 7*( (LAM1+LAM3)/5 - LAM1*LAM3/3                          &
                       + TWONDM*W(9,4)*LAM0*(LAM0-LAM1)*(LAM0-LAM3) ) )       &
              /( 14*LAM2*(LAM2-LAM1)*(LAM2-LAM3) ) - 2*(N-1)*W(7,4)
      W(2,4) = ( 2 - 7*( (LAM2+LAM3)/5 - LAM2*LAM3/3                          &
                       + TWONDM*W(9,4)*LAM0*(LAM0-LAM2)*(LAM0-LAM3) ) )       &
              /( 14*LAM1*(LAM1-LAM2)*(LAM1-LAM3) )                            &
              - 2*(N-1)*( W(6,4) + W(7,4) + (N-2)*W(8,4) )
      W(2,5) = 1/( 6*LAM1 )
!
!     Set generator values
!
      LAM0 = SQRT(LAM0)
      LAM1 = SQRT(LAM1)
      LAM2 = SQRT(LAM2)
      LAM3 = SQRT(LAM3)
      LAMP = SQRT(LAMP)
      G(  :,9) = LAM0
      G(1:3,8) = LAM1
      G(  1,7) = LAM1
      G(  2,7) = LAM2
      G(1:2,6) = LAM1
      G(  1,5) = LAMP
      G(  1,4) = LAM3
      G(  1,3) = LAM2
      G(  1,2) = LAM1
!
!     Compute constant weight values.
!
      W(1,1:5) = 1 - MATMUL( RULPTS(2:WTS), W(2:WTS,1:5) ) 
!
!     Compute final weight values; null rule weights are computed as 
!     differences between weights from highest degree and lower degree rules.
!
      W(1:WTS,2:5) = W(1:WTS,2:5) - SPREAD( W(1:WTS,1), 2, 4 )
!
!        Orthogonalize and normalize null rules.
!
      TEMP = SUM( RULPTS*W(1:WTS,1)*W(1:WTS,1) )
      W(1:WTS,2) = W(1:WTS,2)*SQRT( TEMP/SUM( RULPTS*W(1:WTS,2)*W(1:WTS,2) ) )
      DO K = 3, 5
         ALPHA(1:K-2) = -MATMUL( TRANSPOSE(W(1:WTS,2:K-1)), RULPTS*W(1:WTS,K) )
         W(1:WTS,K) = W(1:WTS,K) + MATMUL( W(1:WTS,2:K-1), ALPHA(1:K-2) )/TEMP
         W(1:WTS,K) = W(1:WTS,K)*SQRT(TEMP/SUM(RULPTS*W(1:WTS,K)*W(1:WTS,K)))
      END DO
      WTSR = WTS
      NUMR = SUM( RULPTS )
!
END SUBROUTINE Rule_Deg9
!
!
FUNCTION SymCub_Sum( N, VERTEX, GIN, NF, Integrand ) RESULT(SymCubSum)
!
!***BEGIN PROLOGUE SymCub_Sum
!***KEYWORDS fully symmetric sum
!***PURPOSE  To compute fully symmetric basic rule sums
!***AUTHOR
!
!            Alan Genz
!            Department of Mathematics
!            Washington State University
!            Pullman, WA 99164-3113, USA
!
!***LAST MODIFICATION 99-06
!***DESCRIPTION SymCub_Sum computes a fully symmetric sum for a vector of
!            integrand values over a parallelepiped. The sum is taken over all
!            sign combinations and permutations of the generators for the sum.
!
!   ON ENTRY
!
!   N       Integer, number of variables.
!   VERTEX  Real array of dimension (N,0:N)
!           The vertices of the simplex, one vertex per column.
!   NF      Integer, number of components for the vector integrand.
! Integrand Real vector function of length NF for computing components of 
!            the integrand at Z.
!            It must have parameters ( NF, Z ); see interface below
!            Input parameters:
!              Z      Real array of length N, the evaluation point.
!              NF     Integer number of components of Integrand.
!   GIN    Real Array of dimension (1:N).
!          The generators for the fully symmetric sum. 

!
!   ON RETURN
!
! SymCub_Sum Real array of length NF, the values for the fully symmetric 
!            sums for each component of the integrand.
!
!***ROUTINES CALLED: Integrand
!
!***END PROLOGUE SymCub_Sum
!
!   Global variables.
!
      INTEGER,                          INTENT(IN) :: N, NF
      INTERFACE 
         FUNCTION Integrand(NF,Z) RESULT(Value)
            USE Precision_Model
            INTEGER,                       INTENT(IN) :: NF
            REAL(KIND=STND), DIMENSION(:), INTENT(IN) :: Z
            REAL(KIND=STND), DIMENSION(NF)            :: Value
         END FUNCTION Integrand
      END INTERFACE
      REAL(KIND=STND), DIMENSION(:,0:), INTENT(IN) :: VERTEX
      REAL(KIND=STND), DIMENSION(:),    INTENT(IN) :: GIN
      REAL(KIND=STND), DIMENSION(NF)               :: SymCubSum
!
!   Local variables.
!
      INTEGER                         :: IX, JX, I, J
      REAL(KIND=STND), DIMENSION(N,N) :: VERDIF
      REAL(KIND=STND), DIMENSION(N)   :: CENTER, G
      REAL(KIND=STND)                 :: GI, GJ
!
!***FIRST PROCESSING STATEMENT SymCub_Sum
!
      SymCubSum = 0
      G = ABS( GIN )      
!
!     Sort input generators if necessary
!
      DO I = 2, N
         IF ( G(I) > G(I-1) ) THEN
            GI = G(I)
            DO J = I-1, 1, -1
               IF ( GI <= G(J) ) THEN
                  EXIT                 
               END IF
               G(J+1) = G(J)
            END DO
            G(J+1) = GI            
         END IF
      END DO
      VERDIF = ( VERTEX(:,1:N) - SPREAD( VERTEX(:,0), 2, N ) )/2
      CENTER = VERTEX(:,0) + SUM( VERDIF, 2 )
!
!     Compute integrand values for sign changes and permutations of G
!
      DO 
         DO 
            SymCubSum = SymCubSum + Integrand( NF, CENTER + MATMUL(VERDIF,G) )
            DO I = 1, N
                G(I) = - G(I)
                IF ( G(I) < 0 ) THEN
                   EXIT
                END IF
            END DO
            IF ( I > N ) THEN
               EXIT
            END IF
         END DO
!
!        Find next distinct permuation of G and loop back for value.
!        Permutations are generated in reverse lexicographic order.
!
         DO I = 2, N
            IF ( G(I-1) > G(I) ) THEN
               GI = G(I)
               IX = I - 1
               DO J = 1, (I-1)/2
                  GJ = G(J)
                  G(J) = G(I-J)
                  G(I-J) = GJ
                  IF (  GJ <= GI ) THEN
                     IX = IX - 1
                  END IF
                  IF ( G(J) > GI ) THEN
                     JX = J
                  END IF
               END DO
               IF ( G(IX) <= GI ) THEN
                  IX = JX
               END IF
               G(I) = G(IX)
               G(IX) = GI
               EXIT
            END IF
         END DO
         IF ( I > N ) THEN
            EXIT
         END IF
      END DO
!
      END Function SymCub_Sum
!
END Module CubatureRule_Cn
