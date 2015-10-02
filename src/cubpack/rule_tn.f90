! This file is F-compatible, except for upper/lower case conventions.
!--------------------------------------------------------------------
MODULE CubatureRule_Tn
USE Precision_Model
IMPLICIT NONE

PRIVATE

PUBLIC :: Rule_Tn
PRIVATE :: RuleParms_Tn, SymSmp_Sum

CONTAINS 
   SUBROUTINE Rule_Tn( TUNE, NDIM, VERTEX, VOLUME, NF, Integrand,             &
                       INKEY, BASVAL, RGNERR, FVALT )
!***BEGIN PROLOGUE Rule_Tn
!***KEYWORDS basic numerical integration rule
!***PURPOSE  To compute basic integration rule values.
!***AUTHOR
!
!            Alan Genz
!            Department of Mathematics
!            Washington State University
!            Pullman, WA 99164-3113, USA
!            AlanGenz@wsu.edu
!
!***LAST MODIFICATION 02-07-09
!***DESCRIPTION Rule_Tn computes basic integration rule values for a
!            vector of integrands over a hyper-rectangular region.
!            These are estimates for the integrals. Rule_Tn also computes
!            estimates for the errors.
!
!   ON ENTRY
!
!     TUNE   Real, tuning parameter.
!     NDIM   Integer, number of variables.
!     VERTEX Real array of dimension (NDIM,0:NDIM).
!            The simplex vertices; vertex J must have components
!            VERTEX(I,J), I = 1, 2, ..., NDIM.
!     VOLUME Real, volume of simplex.
!     NF     Integer, number of components for the vector integrand.
!  Integrand Real vector function of length NF for computing components of 
!            the integrand at Z.
!            It must have parameters ( NF, Z ); see interface below
!            Input parameters:
!              Z      Real array of length NDIM, the evaluation point.
!              NF     Integer number of components of Integrand.
!     INKEY  Integer rule parameter. 
!            If INKEY > 0 and INKEY < 5 then a rule of degree 2*INKEY + 1 
!            is used; otherwise a default rule of degree 7 is used.
!
!   ON RETURN
!
!     BASVAL Real array of length NF, values for the basic rule for 
!            each component of the integrand.
!     RGNERR Real array of length NF, error estimates for BASVAL.
!     FVALT  Integer, number of integrand values USEd by Rule_Tn.
!
!
!***ROUTINES CALLED: RuleParms_Tn, SymSmp_Sum
!
!***END PROLOGUE Rule_Tn
!
!   Global variables.
!
      INTEGER,                             INTENT(IN)  :: NF, NDIM, INKEY
      INTEGER,                             INTENT(OUT) :: FVALT
      INTERFACE 
         FUNCTION Integrand(NF,X) RESULT(Value)
            USE Precision_Model
            INTEGER,                       INTENT(IN)  :: NF
            REAL(KIND=STND), DIMENSION(:), INTENT(IN)  :: X
            REAL(KIND=STND), DIMENSION(NF)             :: Value
         END FUNCTION Integrand
      END INTERFACE
      REAL(KIND=STND), DIMENSION(:,0:),    INTENT(IN)  :: VERTEX
      REAL(KIND=STND),                     INTENT(IN)  :: VOLUME, TUNE
      REAL(KIND=STND), DIMENSION(:),       INTENT(OUT) :: BASVAL, RGNERR
!
!   Local variables.
!
!   WTS    Integer number of weights in the integration rules.
!   W      Real array of dimension (WTS,RLS).
!          The weights for the basic and null rules.
!          W(1,1),...,W(WTS,1) are weights for the basic rule.
!          W(1,I),...,W(WTS,I), for I > 1 are null rule weights.
!   G      Real array of dimension (0:NDIM, WTS).
!          The fully symmetric sum generators for the rules.
!          G(0, J), ..., G(NDIM, J) are the are the generators for the
!          points associated with the Jth weights.
!   X      Real work array of length NDIM.
!   GT     Real work array of length 0:NDIM.
!   RULE   Real work array of dimension (NF,MXRLS).
!   ERROR  Real work array of length NF.
!   RATIO  Real work array of length NF.
!
      INTEGER, PARAMETER                      :: MXW = 21, MXRLS = 7, MXG = 4 
      REAL(KIND=STND), PARAMETER                  :: ONE = 1
      REAL(KIND=STND), PARAMETER                  :: SMALL = 100*EPSILON(ONE) 
      INTEGER, SAVE                               :: OLDKEY = -1, OLDN = 0
      INTEGER, SAVE                               :: RLS, WTS, KEY, FVALS
      INTEGER, DIMENSION(MXW), SAVE               :: PTS
      REAL(KIND=STND), DIMENSION(MXW,MXRLS), SAVE :: W
      REAL(KIND=STND), DIMENSION(0:MXG,MXW), SAVE :: G
      REAL(KIND=STND), DIMENSION(NF,MXRLS)        :: RULE
      REAL(KIND=STND), DIMENSION(MXRLS)           :: ALPHA
      REAL(KIND=STND), DIMENSION(0:NDIM)          :: GTEMP 
      REAL(KIND=STND), DIMENSION(NF)              :: RATIO
      REAL(KIND=STND)                             :: NORMCF, NORMNL, ERRCOF
      INTEGER                                     :: K
!
!***FIRST PROCESSING STATEMENT Rule_Tn
!
      IF ( OLDKEY /= INKEY .OR. OLDN /= NDIM ) THEN
         OLDN = NDIM
         OLDKEY = INKEY
         IF ( INKEY > 0 .AND. INKEY < 5 ) THEN
            KEY = INKEY
         ELSE
            KEY = 3
         END IF
!
!        Compute WTS, RLS, weights, generators, ERRCOF and PTS.
!
         CALL RuleParms_Tn( NDIM, KEY, W, G, WTS, RLS, PTS )
!
!        Orthogonalize and normalize null rules.
!
         NORMCF = DOT_PRODUCT( PTS(1:WTS)*W(1:WTS,1), W(1:WTS,1) )
         DO K = 2, RLS
            ALPHA(2:K-1) = -MATMUL( TRANSPOSE(W(:,2:K-1)), PTS*W(:,K) ) 
            W(:,K) = W(:,K) + MATMUL( W(:,2:K-1), ALPHA(2:K-1) )/NORMCF
            NORMNL = DOT_PRODUCT( PTS*W(:,K), W(:,K) )
            W(:,K) = W(:,K)*SQRT( NORMCF/NORMNL )
         END DO
         FVALS = SUM( PTS(1:WTS) )
      END IF
!
!     Compute the rule values.
!
      RULE = 0
      DO K = 1, WTS
         IF ( PTS(K) > 0 ) THEN
            GTEMP( 0: MIN(NDIM,MXG-1) ) = G( 0: MIN(NDIM,MXG-1) , K )
            IF ( NDIM >= MXG ) THEN
               GTEMP(MXG:NDIM) = G(MXG,K)
            END IF
            BASVAL = SymSmp_Sum( NDIM, VERTEX, NF, Integrand, GTEMP )
            RULE = RULE + MATMUL( RESHAPE( BASVAL, (/ NF, 1 /) ), W(K:K,:) ) 
         END IF
      END DO
      BASVAL = RULE(:,1)
!
!     Scale integral values and compute the error estimates.
!
      ERRCOF = 7*TUNE + 1
      RATIO = 0
      RULE(:,RLS) = MAX ( ABS(RULE(:,RLS)) , ABS(RULE(:,RLS-1)) )
      RGNERR = RULE(:,RLS)
      IF ( KEY > 1 ) THEN
         DO K = RLS-2, 3, -2
            RULE(:,K) = MAX( ABS(RULE(:,K  )) , ABS(RULE(:,K-1)) )
            WHERE ( ABS(BASVAL) + RULE(:,K)/(100*KEY) > ABS(BASVAL) )        
               RATIO = MAX( RULE(:,K)/RULE(:,K+2), RATIO )
            END WHERE
            RGNERR = MAX( RULE(:,K), RGNERR )
         END DO
         RATIO = MAX( ONE/10, RATIO )
         WHERE ( RATIO >= 1 ) 
            RGNERR = TUNE*RGNERR + ( 1 - TUNE )*RULE(:,3)
         ELSEWHERE 
            RGNERR = RATIO*RULE(:,3)
         END WHERE
      END IF
      RGNERR = VOLUME*MAX( ERRCOF*RGNERR, SMALL*ABS( BASVAL ) )
      BASVAL = VOLUME*BASVAL
      FVALT = FVALS
      RETURN
!
!***END Rule_Tn
!
   END Subroutine Rule_Tn
!
   SUBROUTINE RuleParms_Tn( NDIM, KEY, W, G, WTS, RLS, PTS )
!
!***BEGIN PROLOGUE RuleParms_Tn
!***KEYWORDS basic integration rule, degree 2*KEY+1
!***PURPOSE  To initialize a degree 2*KEY+1 basic rule and null rules.
!***AUTHOR
!
!            Alan Genz
!            Department of Mathematics
!            Washington State University
!            Pullman, WA 99164-3113, USA
!            AlanGenz@wsu.edu
!
!            Ronald Cools, Dept. of Computer Science,
!            Katholieke Universiteit Leuven, Celestijnenlaan 200A,
!            B-3001 Heverlee, Belgium
!            Email:  Ronald.Cools@cs.kuleuven.ac.be
!
!***LAST MODIFICATION by Alan 99-05
!***LAST MODIFICATION by Ronald 01-07-19 (cleaning code)
!***DESCRIPTION  RuleParms_Tn initializes a degree 2*KEY+1 rule, and
!                and max(2*KEY,2) lower degree null rules.
!
!   ON ENTRY
!
!   NDIM    Integer, number of variables.
!   KEY    Integer,  < 5 and >= 0, rule parameter.
!          If KEY > 0 a degree 2*KEY+1 rule is initialized.
!          If KEY = 0 a degree 7 rule is initialized.
!
!   ON RETURN
!   RLS    Integer, total number of rules.
!   WTS    Integer, total number of weights in each of the rules.
!   W      Real array of dimension (MXW,*).
!          The weights for the basic and null rules.
!          W(1,1),...,W(WTS,1) are weights for the basic rule.
!          W(I,1),...,W(WTS,I) for I > 1 are null rule weights.
!   G      Real array of dimension (0:MXG,MXW).
!          The fully symmetric sum generators for the rules.
!          G(0,J), ..., G(MXG,J) are the generators for the
!          points associated with the Jth weights.
!   PTS    Integer array of length (MXW). PTS(J) is the number of integrand 
!          values needed for generator J.
!
!***REFERENCES
!
!  Axel Grundmann and H. M. Moller
!  "Invariant Integration Formulas for the n-Simplex by Combinatorial Methods",
!   SIAM J Numer. Anal. 15(1978), 282--290,
! and
!  A. H. Stroud
!  "A Fifth Degree Integration Formula for the n-Simplex
!  SIAM J Numer. Anal. 6(1969), 90--98,
! and           
!  I. P. Mysovskikh
!  "On a cubature formula for the simplex"
!  Vopros. Vycisl. i Prikl. Mat., Tashkent 51(1978), 74--90.
!
!
!***ROUTINES CALLED NONE
!***END PROLOGUE RuleParms_Tn
!
!   Global variables
!
      INTEGER, INTENT(IN)                           :: NDIM, KEY
      INTEGER, INTENT(OUT)                          :: WTS, RLS
      INTEGER, DIMENSION(:),            INTENT(OUT) :: PTS
      REAL(KIND=STND), DIMENSION(:,:),  INTENT(OUT) :: W
      REAL(KIND=STND), DIMENSION(0:,:), INTENT(OUT) :: G
!
!   Local Variables
!
      REAL(KIND=STND), PARAMETER :: ONE = 1, FFTEEN = 15
      REAL(KIND=STND) :: DR, DR2, DR4, DR6, DR8
      REAL(KIND=STND) :: R1, S1, R2, S2, U1, V1, U2, V2, L1, L2, D1, D2
      REAL(KIND=STND) :: A1, A2, A3, P0, P1, P2, P3, U5, U6, U7, SG
      REAL(KIND=STND) :: R, A, P, Q, TH, TP
      INTEGER         :: IW, GMS
!
!***FIRST PROCESSING STATEMENT RuleParms_Tn
!
!
!     Initialize RLS and GMS.
!
      IF ( KEY == 1 ) THEN
         RLS = 3
         GMS = 2
         WTS = 3
      ELSE IF ( KEY == 2 ) THEN
         RLS = 5
         GMS = 4
         WTS = 6
      ELSE IF ( KEY == 3 .OR. KEY == 0 ) THEN
         RLS = 7
         GMS = 7
         WTS = 11
      ELSE IF ( KEY == 4 ) THEN
         RLS = 7
         IF ( NDIM == 2 ) THEN
            GMS = 11
            WTS = 20
         ELSE
            GMS = 12
            WTS = 21
         END IF
      END IF
!
!     Initialize generators, weights and PTS.
!
      W(:,1:RLS) = 0
      PTS = 0

!
!     Compute generator, PTS and weight values for all rules.
!
      DR = NDIM
      DR2 =    ( DR + 1 )*( DR + 2 )
      DR4 = DR2*( DR + 3 )*( DR + 4 )
      DR6 = DR4*( DR + 5 )*( DR + 6 )
      DR8 = DR6*( DR + 7 )*( DR + 8 )
      G(0:,1) = 1/( DR + 1 )
          PTS(1) = 1
      R1 = ( DR + 4 - SQRT(FFTEEN) )/( DR*DR + 8*DR + 1 )
      S1 = 1 - DR*R1
      L1 = S1 - R1
      G(0    ,GMS+1) = S1
      G(1:,GMS+1) = R1
          PTS(GMS+1) = NDIM + 1
      IW = RLS
      IF ( KEY < 4 )  THEN
!
!        Compute weights for special degree 1 rule.
!
         W(1,IW) = 1
         IW = IW - 1
         W(GMS+1,IW) = 1/( DR + 1 )
         IW = IW - 1
      END IF
!
!     Compute weights, generators and PTS for degree 3 rule.
!
      G(0    ,2) = 3/( DR + 3 )
      G(1:,2) = 1/( DR + 3 )
          PTS(2) = NDIM + 1
            W(2,IW) = ( DR + 3 )**3/( 4*DR2*( DR + 3 ) )
      IF ( KEY > 1 ) THEN
         IW = IW - 1
!
!        Compute weights, generators and PTS for degree 3 and degree 5 rules.
!
         IF ( NDIM == 2 ) THEN
!
!           Special degree 3 rule.
!
            L2 = 0.6205464826720063258904603436171006977619_STND
            L1 = -SQRT( ONE/2 - L2**2 )
            R1 = ( 1 - L1 )/3
            S1 = 1 - 2*R1
            G(0    ,GMS+1) = S1
            G(1:,GMS+1) = R1
                PTS(GMS+1) = 3
                  W(GMS+1,IW) = ONE/6
            R2 = ( 1 - L2 )/3
            S2 = 1 - 2*R2
            G(0    ,GMS+2) = S2
            G(1:,GMS+2) = R2
                PTS(GMS+2) = 3
                  W(GMS+2,IW) = ONE/6
         ELSE
!
!           Degree 3 rule using Stroud points.
!
            R2 = ( DR + 4 + SQRT(FFTEEN) )/( DR*DR + 8*DR + 1 )
            S2 = 1 - DR*R2
            L2 = S2 - R2
            G(0    ,GMS+2) = S2
            G(1:,GMS+2) = R2
                PTS(GMS+2) = NDIM + 1
                  W(GMS+2,IW) = ( 2/(DR+3) - L1 )/(DR+1)/(DR+2)/(L2-L1)/L2**2
                  W(GMS+1,IW) = ( 2/(DR+3) - L2 )/(DR+1)/(DR+2)/(L1-L2)/L1**2
         END IF
         IW = IW - 1
!
!        Grundmann-Moller degree 5 rule.
!
         G(0    ,3) = 5/( DR + 5 )
         G(1:,3) = 1/( DR + 5 )
             PTS(3) = NDIM + 1
         G(0:1  ,4) = 3/( DR + 5 )
         G(2:,4) = 1/( DR + 5 )
             PTS(4) = ( ( NDIM + 1 )*NDIM )/2
             W(2,  IW) = -( DR + 3 )**5/( 16*DR4 )
             W(3:4,IW) =  ( DR + 5 )**5/( 16*DR4*( DR + 5 ) )
      END IF
      IF ( KEY > 2 )  THEN
         IW = IW - 1
!
!        Compute weights, generators and PTS for degree 5 and degree 7 rules.
!
!
!        Stroud degree 5 rule.
!
         U1 = ( DR + 7 + 2*SQRT(FFTEEN) )/( DR*DR + 14*DR - 11 )
         V1 = ( 1 - ( DR - 1 )*U1 )/2
         D1 = V1 - U1
         G(0:1  ,GMS+3) = V1
         G(2:,GMS+3) = U1
             PTS(GMS+3) = ( ( NDIM + 1 )*NDIM )/2
         U2 = ( DR + 7 - 2*SQRT(FFTEEN) )/( DR*DR + 14*DR - 11 )
         V2 = ( 1 - ( DR - 1 )*U2 )/2
         D2 = V2 - U2
         G(0:1  ,GMS+4) = V2
         G(2:,GMS+4) = U2
             PTS(GMS+4) = ( ( NDIM + 1 )*NDIM )/2
         IF ( NDIM == 2 ) THEN
               W(GMS+3,IW) = ( 155 - SQRT(FFTEEN) )/1200
               W(GMS+4,IW) = ( 155 + SQRT(FFTEEN) )/1200
               W(1,    IW) = 1 - 3*( W(GMS+3,IW) + W(GMS+4,IW) ) 
         ELSE IF ( NDIM == 3 ) THEN
               W(GMS+1,IW) = ( 2665 + 14*SQRT(FFTEEN) )/37800
               W(GMS+2,IW) = ( 2665 - 14*SQRT(FFTEEN) )/37800
               W(GMS+3,IW) = 2*FFTEEN/567
             PTS(GMS+4) = 0
         ELSE
               W(GMS+1,IW) = ( 2*(27-DR)/(DR+5)-L2*(13-DR) )/L1**4/(L1-L2)/DR4
               W(GMS+2,IW) = ( 2*(27-DR)/(DR+5)-L1*(13-DR) )/L2**4/(L2-L1)/DR4
               W(GMS+3,IW)=( 2/( DR + 5 ) - D2 )/( DR4*( D1 - D2 )*D1**4 )
               W(GMS+4,IW)=( 2/( DR + 5 ) - D1 )/( DR4*( D2 - D1 )*D2**4 )
         END IF
         IW = IW - 1
!
!        Grundmann-Moller degree 7 rule.
!
         G(0    ,5) = 7/( DR + 7 )
         G(1:,5) = 1/( DR + 7 )
             PTS(5) = NDIM + 1 
         G(0    ,6) = 5/( DR + 7 )
         G(1    ,6) = 3/( DR + 7 )
         G(2:,6) = 1/( DR + 7 )
             PTS(6) = ( NDIM + 1 )*NDIM
         G(0:2  ,7) = 3/( DR + 7 )
         G(3:,7) = 1/( DR + 7 )
             PTS(7) = ( ( NDIM + 1 )*NDIM*( NDIM - 1 ) )/6
             W(2,  IW) =  ( DR + 3 )**7/( 2*64*DR4*( DR + 5 ) )
             W(3:4,IW) = -( DR + 5 )**7/(   64*DR6 )
             W(5:7,IW) =  ( DR + 7 )**7/(   64*DR6*( DR + 7 ) )
      END IF
      IF ( KEY == 4 )  THEN
         IW = IW - 1
!
!        Compute weights, generators and PTS for degree 7 and degree 9 rules.
!
!        Mysovskikh degree 7 rule.
!
         SG = 1/( 23328*DR6 )
         U5 = -6**3*SG*( 52212 - DR*( 6353 + DR*( 1934 - DR*27 ) ) )       
         U6 =  6**4*SG*(  7884 - DR*( 1541 - DR*9 ) )
         U7 = -6**5*SG*(  8292 - DR*( 1139 - DR*3 ) )/( DR + 7 )
         P0 = -144*( 142528 + DR*( 23073 - DR*115 ) )
         P1 = -12*( 6690556 + DR*( 2641189 + DR*( 245378 - DR*1495 ) ) )
         P2 = -16*(6503401 + DR*( 4020794+DR*(787281+DR*(47323-DR*385)) ) )
         P3 = -(6386660+DR*(4411997+DR*(951821+DR*(61659-DR*665))))*(DR+7)
!---------------------------------------------------------------------------
!        Compute 3 zeros by 4 Newton iterations (good for about 30 digits)
!        A1 = -2/( DR + 3 )
!        A1 = A1 - ( P0+A1*(P1+A1*(P2+A1*P3)) )/( P1+A1*(2*P2+A1*3*P3) )
!        A1 = A1 - ( P0+A1*(P1+A1*(P2+A1*P3)) )/( P1+A1*(2*P2+A1*3*P3) )
!        A1 = A1 - ( P0+A1*(P1+A1*(P2+A1*P3)) )/( P1+A1*(2*P2+A1*3*P3) )
!        A1 = A1 - ( P0+A1*(P1+A1*(P2+A1*P3)) )/( P1+A1*(2*P2+A1*3*P3) )
!        G(0    ,GMS+5) = ( 1 - DR*A1 )/( DR + 1 )
!        G(1:,GMS+5) =   ( 1 + A1 )/( DR + 1 )
!            PTS(GMS+5) = NDIM + 1
!        P2 = P2 + A1*P3
!        P1 = P1 + A1*P2
!        A2 = ( -P2 - SQRT( P2**2 - 4*P1*P3 ) )/( 2*P3 )
!        G(0    ,GMS+6) = ( 1 - DR*A2 )/( DR + 1 )
!        G(1:,GMS+6) =   ( 1 + A2 )/( DR + 1 )
!            PTS(GMS+6) = NDIM + 1
!        A3 = ( -P2 + SQRT( P2**2 - 4*P1*P3 ) )/( 2*P3 )
!---------------------------------------------------------------------------
!        Compute 3 zeros by closed formula due to Cardan/Tartaglia
         A = P2/( 3*P3 )
         P = A*( P1/P2 - A )
         Q = A*( 2*A*A - P1/P3 ) + P0/P3
         R = SQRT( -P**3 )
         TH = ACOS( -Q/( 2*R ) )/3
         R = 2*R**( ONE/3 )
         TP = 2*ACOS(-ONE)/3
         A1 = -A + R*COS( TH ) 
         G(0    ,GMS+5) = ( 1 - DR*A1 )/( DR + 1 )
         G(1:,GMS+5) =   ( 1 + A1 )/( DR + 1 )
             PTS(GMS+5) = NDIM + 1
         A2 = -A + R*COS( TH + TP + TP )
         G(0    ,GMS+6) = ( 1 - DR*A2 )/( DR + 1 )
         G(1:,GMS+6) =   ( 1 + A2 )/( DR + 1 )
             PTS(GMS+6) = NDIM + 1
         A3 = -A + R*COS( TH + TP )
!---------------------------------------------------------------------------
         G(0    ,GMS+7) = ( 1 - DR*A3 )/( DR + 1 )
         G(1:,GMS+7) =   ( 1 + A3 )/( DR + 1 )
             PTS(GMS+7) = NDIM + 1
               W(GMS+5,IW) =                                                 &
                  ( U7-(A2+A3)*U6+A2*A3*U5 )/( A1**2-(A2+A3)*A1+A2*A3 )/A1**5
               W(GMS+6,IW) =                                                 &
                  ( U7-(A1+A3)*U6+A1*A3*U5 )/( A2**2-(A1+A3)*A2+A1*A3 )/A2**5
               W(GMS+7,IW) =                                                 &
                  ( U7-(A2+A1)*U6+A2*A1*U5 )/( A3**2-(A2+A1)*A3+A2*A1 )/A3**5
         G(0:1  ,GMS+8) = 4/( DR + 7 )
         G(2:,GMS+8) = 1/( DR + 7 )
             PTS(GMS+8) = ( ( NDIM + 1 )*NDIM )/2
               W(GMS+8,IW) = 10*(DR+7)**6/( 729*DR6 )
         G(0    ,GMS+9) = 11/( DR + 7 )/2
         G(1    ,GMS+9) =  5/( DR + 7 )/2
         G(2:,GMS+9) =  1/( DR + 7 )
             PTS(GMS+9) = ( ( NDIM + 1 )*NDIM )
               W(GMS+9,IW) = 64*(DR+7)**6/( 6561*DR6 )
               W(    4,IW) = W(4,IW+1)
               W(    7,IW) = W(7,IW+1)
         IW = IW - 1
!
!        Grundmann-Moller degree 9 rule.
!
         G(0    ,8) = 9/( DR + 9 )
         G(1:,8) = 1/( DR + 9 )
             PTS(8) = NDIM + 1 
         G(0    ,9) = 7/( DR + 9 )
         G(1    ,9) = 3/( DR + 9 )
         G(2:,9) =   1/( DR + 9 )
             PTS(9) = ( NDIM + 1 )*NDIM 
         G(0:1  ,10) = 5/( DR + 9 )
         G(2:,10) = 1/( DR + 9 )
             PTS(10) = ( ( NDIM + 1 )*NDIM )/2
         G(0    ,11) = 5/( DR + 9 )
         G(1:2  ,11) = 3/( DR + 9 )
         G(3:,11) = 1/( DR + 9 )
             PTS(11) = ( ( NDIM + 1 )*NDIM*( NDIM - 1 ) )/2
             W(2   ,IW) = -( DR + 3 )**9/( 6*256*DR6 )
             W(3:4 ,IW) =  ( DR + 5 )**9/( 2*256*DR6*(DR+7) )
             W(5:7 ,IW) = -( DR + 7 )**9/(   256*DR8 )
             W(8:11,IW) =  ( DR + 9 )**9/(   256*DR8*(DR+9) )
         IF ( NDIM > 2 ) THEN
            G(0:3  ,12) = 3/( DR + 9 )
            G(4:,12) = 1/( DR + 9 )
                PTS(12) = ( ( NDIM + 1 )*NDIM*( NDIM - 1 )*( NDIM - 2 ) )/24
                  W(12,IW) =    W(8,IW)
         END IF         
      END IF
!
!     Compute constant weight values.
!
      W(1,1:RLS) = 1 - MATMUL( PTS(2:WTS), W(2:WTS,1:RLS) ) 
!
!     Compute final weight values; null rule weights are computed as 
!     differences between weights from highest degree and lower degree rules.
!
      W(:,2:RLS) = W(:,2:RLS) - SPREAD( W(:,1), 2, RLS-1 )
!
      RETURN
   END Subroutine RuleParms_Tn
!
   FUNCTION SymSmp_Sum( N, VERTEX, NF, Integrand, GIN ) RESULT(SymSmpSum)
!
!***BEGIN PROLOGUE SymSmp_Sum
!***KEYWORDS fully symmetric sum
!***PURPOSE  To compute fully symmetric basic rule sums
!***AUTHOR
!
!            Alan Genz
!            Department of Mathematics
!            Washington State University
!            Pullman, WA 99164-3113, USA
!
!***LAST MODIFICATION 99-05
!***DESCRIPTION SymSmp_Sum computes a fully symmetric sum for a vector
!            of integrand values over a simplex. The sum is taken over
!            all permutations of the generators for the sum.
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
!   GIN    Real Array of dimension (0:N).
!          The generators for the fully symmetric sum. 

!
!   ON RETURN
!
! SymSmp_Sum Real array of length NF, the values for the fully symmetric 
!            sums for each component of the integrand.
!
!***ROUTINES CALLED: Integrand
!
!***END PROLOGUE SymSmp_Sum
!
!   Global variables.
!
      INTEGER,                          INTENT(IN) :: N, NF
      INTERFACE 
         FUNCTION Integrand(NF,Z) RESULT(Value)
            USE Precision_Model
            INTEGER, INTENT(IN) :: NF
            REAL(KIND=STND), DIMENSION(:), INTENT(IN) :: Z
            REAL(KIND=STND), DIMENSION(NF) :: Value
         END FUNCTION Integrand
      END INTERFACE
      REAL(KIND=STND), DIMENSION(:,0:), INTENT(IN) :: VERTEX
      REAL(KIND=STND), DIMENSION(0:),   INTENT(IN) :: GIN
      REAL(KIND=STND), DIMENSION(NF)               :: SymSmpSum
!
!   Local variables.
!
      INTEGER                         :: IX, LX, I, J, K, L
      REAL(KIND=STND), DIMENSION(0:N) :: G
      REAL(KIND=STND)                 :: GL, GI
!
!***FIRST PROCESSING STATEMENT SymSmp_Sum
!
      SymSmpSum = 0
      G = GIN
!
!     Sort input generators if necessary
!
      K = 0
      DO I = 1, N
         IF ( G(I) > G(I-1) ) THEN
            K = 1
      END IF
      END DO
      IF ( K > 0 ) THEN
         DO I = 1, N
            K = I - 1
            DO J = I, N
               IF ( G(J) > G(K) ) THEN
                  K = J
               END IF
            END DO
            IF ( K >= I ) THEN
               GI = G(I-1)
               G(I-1) = G(K)
               G(K) = GI
            END IF
         END DO
      END IF
!
!     Compute integrand value for permutations of G
!
      DO 
         SymSmpSum = SymSmpSum + Integrand( NF, MATMUL( VERTEX, G ) )
!
!        Find next distinct permuation of G and loop back for value.
!        Permutations are generated in reverse lexicographic order.
!
         DO I = 1, N
            IF ( G(I-1) > G(I) ) THEN
               GI = G(I)
               IX = I - 1
               DO L = 0, I/2-1
                  GL = G(L)
                  G(L) = G(I-L-1)
                  G(I-L-1) = GL
                  IF (  GL <= GI ) THEN
                     IX = IX - 1
                  END IF
                  IF ( G(L) > GI ) THEN
                     LX = L
                  END IF
               END DO
               IF ( G(IX) <= GI ) THEN
                  IX = LX
               END IF
               G(I) = G(IX)
               G(IX) = GI
               EXIT
            END IF
         END DO
         IF ( I == N+1 ) THEN
            EXIT
         END IF
      END DO
!
      RETURN
      END Function SymSmp_Sum
!
END MODULE CubatureRule_Tn
