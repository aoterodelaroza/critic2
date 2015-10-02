! This file is F-compatible, except for upper/lower case conventions.
!--------------------------------------------------------------------

Module CubatureRule_T2

USE Precision_Model, ONLY: stnd

Implicit NONE

PRIVATE
PUBLIC :: Rule_T2a

CONTAINS
      SUBROUTINE Rule_T2a(VER,AREA,NUMFUN,Integrand,BASVAL,RGNERR,NUM)


!***BEGIN PROLOGUE Rule_T2a
!***PURPOSE  To compute basic integration rule values and
!            corresponding error estimates.
!***REFER TO Module CubatureRule_General
!      This subroutine is based on DRLTRI, part of DCUTRI. See below.
!***REVISION DATE  950823   (YYMMDD)
!***REVISION DATE  990527   (YYMMDD) (F conversion)
!***AUTHOR
!      Original version
!          Jarle Berntsen, The Computing Centre,
!          University of Bergen, Thormohlens gt. 55,
!          N-5008 Bergen, NORWAY
!          Email:  jarle@eik.ii.uib.no
!          Terje O. Espelid, Department of Informatics,
!          University of Bergen, Thormohlens gt. 55,
!          N-5008 Bergen, NORWAY
!          Email:  terje@eik.ii.uib.no
!      Translation and modification by
!          Ronald Cools, Dept. of Computer Science,
!          Katholieke Universiteit Leuven, Celestijnenlaan 200A,
!          B-3001 Heverlee, Belgium
!          Email:  ronald@cs.kuleuven.ac.be
!
!***REFERENCES
!   J.Berntsen and T.O.Espelid
!     Algorithm 706: DCUTRI: An algorithm for adaptive cubature over a
!     collection of triangles
!   ACM. Trans. Math. Software, Vol. 18 (1992), pp 329-342.
!
!***DESCRIPTION Rule_T2a computes basic integration rule values
!            for a vector of integrands over a triangular region.
!            Rule_T2a also computes estimates for the errors by
!            using several null rule approximations.
!
!   ON ENTRY
!
!   VER    Real array of dimension (2,3).
!          The coordinates of the vertices of the triangle.
!   AREA   Real.
!          The area of the given region.
!   NUMFUN Integer.
!          Number of components of the vector integrand.
!   Integrand Externally declared subroutine for computing
!            all components of the integrand at the given
!            evaluation point.
!            It must have parameters X
!            Input parameters:
!              X(1)      The x-coordinate of the evaluation point.
!              X(2)      The y-coordinate of the evaluation point.
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
!***REFERENCES Berntsen,J. and Espelid,T.O., Degree 13 Symmetric
!              Quadrature Rules for the Triangle, Report
!***ROUTINES CALLED Integrand
!***END PROLOGUE Rule_T2a
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
      REAL(kind=stnd), INTENT(IN) :: AREA
      REAL(kind=stnd), DIMENSION(:,:), INTENT(IN) :: VER
      REAL(kind=stnd), DIMENSION(:), INTENT(OUT) :: BASVAL, RGNERR
!
!   Constants
!
!   G      Real array of dimension (2,ORBITS).
!          The homogeneous coordinates for the generators of
!          the evaluation points.
!          The integration rule is using symmetric evaluation
!          points and has the structure (1,6,3). That is,
!          1 point of multiplicity 1,
!          6 sets of points of multiplicity 3 and
!          3 sets of points of multiplicity 6.
!          This gives totally 37 evaluation points.
!          In order to reduce the number of loops in Rule_T2a,
!          the 3 loops for the sets of multiplicity 6 are split
!          into 6 loops and added to the loops for the sets of
!          multiplicity 3.
!          The number of weights we have to give with
!          this splitting is 13(ORBITS).
!
!   W      Real array of dimension (9,ORBITS).
!          The weights of the basic rule and the null rules.
!          W(1,1),...,W(1,ORBITS) are weights for the basic rule.
!          W(I,1),...,W(I,ORBITS) for I>1 are null rule weights.
!
      INTEGER, PARAMETER :: ORBITS=13 ! The number of orbits in the Cf.
      REAL(kind=stnd), PARAMETER  :: CRIVAL=0.5_stnd,             &
                                TRES=50*EPSILON(crival),          &
                                FACMED = 10,                      &
                                FACOPT = FACMED/(CRIVAL**2)
!
!  The abscissas are given in homogeneous coordinates.
!
      REAL(kind=stnd), DIMENSION(ORBITS), PARAMETER ::            &
           G1 = (/                                                &
           0.333333333333333333333333333333_stnd,                 &
           0.950275662924105565450352089520_stnd,                 &
           0.171614914923835347556304795551_stnd,                 &
           0.539412243677190440263092985511_stnd,                 &
           0.772160036676532561750285570113_stnd,                 &
           0.009085399949835353883572964740_stnd,                 &
           0.062277290305886993497083640527_stnd,                 &
           0.022076289653624405142446876931_stnd,                 &
           0.018620522802520968955913511549_stnd,                 &
           0.096506481292159228736516560903_stnd,                 &
           0.851306504174348550389457672223_stnd,                 &
           0.689441970728591295496647976487_stnd,                 &
           0.635867859433872768286976979827_stnd/)
      REAL(kind=stnd), DIMENSION(ORBITS), PARAMETER ::            &
           G2 = (/                                                &
           0.333333333333333333333333333333_stnd,                 &
           0.024862168537947217274823955239_stnd,                 &
           0.414192542538082326221847602214_stnd,                 &
           0.230293878161404779868453507244_stnd,                 &
           0.113919981661733719124857214943_stnd,                 &
           0.495457300025082323058213517632_stnd,                 &
           0.468861354847056503251458179727_stnd,                 &
           0.851306504174348550389457672223_stnd,                 &
           0.689441970728591295496647976487_stnd,                 &
           0.635867859433872768286976979827_stnd,                 &
           0.022076289653624405142446876931_stnd,                 &
           0.018620522802520968955913511549_stnd,                 &
           0.096506481292159228736516560903_stnd/)
      REAL(kind=stnd), DIMENSION(2,ORBITS), PARAMETER ::          &
           G =  RESHAPE( SOURCE= (/ G1 , G2 /), SHAPE=(/2,ORBITS/), ORDER=(/2,1/) )
!
!   Weights of the degree 13 quadrature rule.
!
      REAL(kind=stnd), DIMENSION(ORBITS), PARAMETER ::            &
           W1 = (/                                                &
           0.051739766065744133555179145422_stnd,                 &
           0.008007799555564801597804123460_stnd,                 &
           0.046868898981821644823226732071_stnd,                 &
           0.046590940183976487960361770070_stnd,                 &
           0.031016943313796381407646220131_stnd,                 &
           0.010791612736631273623178240136_stnd,                 &
           0.032195534242431618819414482205_stnd,                 &
           0.015445834210701583817692900053_stnd,                 &
           0.017822989923178661888748319485_stnd,                 &
           0.037038683681384627918546472190_stnd,                 &
           0.015445834210701583817692900053_stnd,                 &
           0.017822989923178661888748319485_stnd,                 &
           0.037038683681384627918546472190_stnd/)
!
!   Weights of the first null rule of degree 7.
!
      REAL(kind=stnd), DIMENSION(ORBITS), PARAMETER ::            &
           W2 = (/                                                &
          -0.077738051051462052051304462750_stnd,                 &
           0.001640389740236881582083124927_stnd,                 &
           0.078124083459915167386776552733_stnd,                 &
          -0.030706528522391137165581298102_stnd,                 &
           0.010246307817678312345028512621_stnd,                 &
           0.012586300774453821540476193059_stnd,                 &
          -0.043630506151410607808929481439_stnd,                 &
          -0.004567055157220063810223671248_stnd,                 &
           0.003393373439889186878847613140_stnd,                 &
           0.0_stnd,                                              &
          -0.004567055157220063810223671248_stnd,                 &
           0.003393373439889186878847613140_stnd,                 &
           0.0_stnd/)
!
!   Weights of the second null rule of degree 7.
!
      REAL(kind=stnd), DIMENSION(ORBITS), PARAMETER ::            &
           W3 = (/                                                &
          -0.064293709240668260928898888457_stnd,                 &
           0.003134665264639380635175608661_stnd,                 &
           0.007822550509742830478456728602_stnd,                 &
           0.048653051907689492781049400973_stnd,                 &
           0.032883327334384971735434067029_stnd,                 &
          -0.017019508374229390108580829589_stnd,                 &
           0.025973557893399824586684707198_stnd,                 &
          -0.010716753326806275930657622320_stnd,                 &
           0.018315629578968063765722278290_stnd,                 &
          -0.047607080313197299401024682666_stnd,                 &
          -0.010716753326806275930657622320_stnd,                 &
           0.018315629578968063765722278290_stnd,                 &
          -0.047607080313197299401024682666_stnd/)
!
!   Weights of the first degree 5 null rule.
!
!
      REAL(kind=stnd), DIMENSION(ORBITS), PARAMETER ::            &
           W4 = (/                                                &
           0.021363205584741860993131879186_stnd,                 &
           0.022716410154120323440432428315_stnd,                 &
          -0.026366191271182090678117381002_stnd,                 &
           0.029627021479068212693155637482_stnd,                 &
           0.004782834546596399307634111034_stnd,                 &
           0.004178667433984132052378990240_stnd,                 &
          -0.065398996748953861618846710897_stnd,                 &
          -0.033589813176131630980793760168_stnd,                 &
           0.033018320112481615757912576257_stnd,                 &
           0.012241086002709814125707333127_stnd,                 &
          -0.033589813176131630980793760168_stnd,                 &
           0.033018320112481615757912576257_stnd,                 &
           0.012241086002709814125707333127_stnd/)
!
!   Weights of the second degree 5 null rule.
!
      REAL(kind=stnd), DIMENSION(ORBITS), PARAMETER ::            &
           W5 = (/                                                &
          -0.046058756832790538620830792345_stnd,                 &
           0.005284159186732627192774759959_stnd,                 &
           0.009325799301158899112648198129_stnd,                 &
          -0.006101110360950124560783393745_stnd,                 &
          -0.056223328794664871336486737231_stnd,                 &
          -0.062516479198185693171971930698_stnd,                 &
           0.022428226812039547178810743269_stnd,                 &
          -0.000026014926110604563130107142_stnd,                 &
           0.032882099937471182365626663487_stnd,                 &
           0.018721740987705986426812755881_stnd,                 &
          -0.000026014926110604563130107142_stnd,                 &
           0.032882099937471182365626663487_stnd,                 &
           0.018721740987705986426812755881_stnd/)
!
!   Weights of first degree 3 null rule.
!
      REAL(kind=stnd), DIMENSION(ORBITS), PARAMETER ::            &
           W6 = (/                                                &
           0.080867117677405246540283712799_stnd,                 &
          -0.033915806661511608094988607349_stnd,                 &
           0.014813362053697845461526433401_stnd,                 &
           0.001442315416337389214102507204_stnd,                 &
          -0.024309696484708683486455879210_stnd,                 &
          -0.005135085639122398522835391664_stnd,                 &
          -0.034649417896235909885490654650_stnd,                 &
           0.035748423431577326597742956780_stnd,                 &
           0.024548155266816447583155562333_stnd,                 &
          -0.032897267038856299280541675107_stnd,                 &
           0.035748423431577326597742956780_stnd,                 &
           0.024548155266816447583155562333_stnd,                 &
          -0.032897267038856299280541675107_stnd/)
!
!   Weights of second degree 3 null rule.
!
      REAL(kind=stnd), DIMENSION(ORBITS), PARAMETER ::            &
           W7 = (/                                                &
          -0.038457863913548248582247346193_stnd,                 &
          -0.055143631258696406147982448269_stnd,                 &
          -0.021536994314510083845999131455_stnd,                 &
           0.001547467894857008228010564582_stnd,                 &
           0.057409361764652373776043522086_stnd,                 &
          -0.040636938884669694118908764512_stnd,                 &
          -0.020801144746964801777584428369_stnd,                 &
           0.019490770404993674256256421103_stnd,                 &
           0.002606109985826399625043764771_stnd,                 &
           0.023893703367437102825618048130_stnd,                 &
           0.019490770404993674256256421103_stnd,                 &
           0.002606109985826399625043764771_stnd,                 &
           0.023893703367437102825618048130_stnd/)
!
!   Weights of first degree 1 null rule.
!
      REAL(kind=stnd), DIMENSION(ORBITS), PARAMETER ::            &
           W8 = (/                                                &
           0.074839568911184074117081012527_stnd,                 &
          -0.004270103034833742737299816615_stnd,                 &
           0.049352639555084484177095781183_stnd,                 &
           0.048832124609719176627453278550_stnd,                 &
           0.001042698696559292759051590242_stnd,                 &
          -0.044445273029113458906055765365_stnd,                 &
          -0.004670751812662861209726508477_stnd,                 &
          -0.015613390485814379318605247424_stnd,                 &
          -0.030581651696100000521074498679_stnd,                 &
           0.010801113204340588798240297593_stnd,                 &
          -0.015613390485814379318605247424_stnd,                 &
          -0.030581651696100000521074498679_stnd,                 &
           0.010801113204340588798240297593_stnd/)
!
!   Weights of second degree 1 null rule.
!
      REAL(kind=stnd), DIMENSION(ORBITS), PARAMETER ::            &
           W9 = (/                                                &
           0.009373028261842556370231264134_stnd,                 &
          -0.074249368848508554545399978725_stnd,                 &
           0.014709707700258308001897299938_stnd,                 &
           0.009538502545163567494354463302_stnd,                 &
          -0.014268362488069444905870465047_stnd,                 &
           0.040126396495352694403045023109_stnd,                 &
           0.028737181842214741174950928350_stnd,                 &
          -0.031618075834734607275229608099_stnd,                 &
           0.016879961075872039084307382161_stnd,                 &
           0.010878914758683152984395046434_stnd,                 &
          -0.031618075834734607275229608099_stnd,                 &
           0.016879961075872039084307382161_stnd,                 &
           0.010878914758683152984395046434_stnd/)
      REAL(kind=stnd), DIMENSION(9,ORBITS), PARAMETER ::          &
           W =  RESHAPE( SOURCE= (/ W1,W2,W3,W4,W5,W6,W7,W8,W9 /),&
                   SHAPE=(/9,ORBITS/), ORDER=(/2,1/) )
!
!  Local variables
!
!
!   NullRule  Real array of dimension (NUMFUN,8).
!          A work array.
!
      INTEGER :: I,J,L
      REAL(kind=stnd):: Z1,Z2,Z3,R1,R2,R3,R,DEG7,DEG5,DEG3,DEG1,NOISE
      REAL(kind=stnd), DIMENSION(NUMFUN,8) :: NullRule
      REAL(kind=stnd), DIMENSION(2,3)      :: X
!
!***FIRST EXECUTABLE STATEMENT Rule_T2a
!
!   Compute contributions from the center of the triangle.
!
      X(1:2,1) = (VER(1:2,1)+VER(1:2,2)+VER(1:2,3))/3
      RGNERR=Integrand(NUMFUN,X(:,1))
      BASVAL(1:NUMFUN) = W(1,1)*RGNERR(1:NUMFUN)
      DO J = 1,NUMFUN
          NullRule(J,1:8) = W(2:9,1)*RGNERR(J)
      END DO
!
!   Compute contributions from points with
!   multiplicity 3.
!
      DO I = 2,ORBITS
          Z1 = G(1,I)
          Z2 = G(2,I)
          Z3 = 1 - Z1 - Z2
          X(1:2,1) = Z1*VER(1:2,1) + Z2*VER(1:2,2) + Z3*VER(1:2,3)
          X(1:2,2) = Z2*VER(1:2,1) + Z3*VER(1:2,2) + Z1*VER(1:2,3)
          X(1:2,3) = Z3*VER(1:2,1) + Z1*VER(1:2,2) + Z2*VER(1:2,3)
          DO L = 1,3
              RGNERR=Integrand(NUMFUN,X(:,L))
              BASVAL(1:NUMFUN) = BASVAL(1:NUMFUN) + W(1,I)*RGNERR(1:NUMFUN)
              DO J = 1,NUMFUN
                  NullRule(J,1:8) = NullRule(J,1:8) + W(2:9,I)*RGNERR(J)
              END DO
          END DO
      END DO
!
!    Compute errors.
!
      DO J = 1,NUMFUN
          NOISE = ABS(BASVAL(J))*TRES
          DEG7 = SQRT(NullRule(J,1)**2+NullRule(J,2)**2)
          IF ( DEG7 <= NOISE) THEN
              RGNERR(J) = NOISE
          ELSE
          DEG5 = SQRT(NullRule(J,3)**2+NullRule(J,4)**2)
          DEG3 = SQRT(NullRule(J,5)**2+NullRule(J,6)**2)
          DEG1 = SQRT(NullRule(J,7)**2+NullRule(J,8)**2)
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
              RGNERR(J) = 10*MAX(DEG1,DEG3,DEG5,DEG7)
          ELSE IF (R >= CRIVAL) THEN
              RGNERR(J) = facmed*R*DEG7
          ELSE
              RGNERR(J) = facopt* (R**3)*DEG7
          END IF
          RGNERR(J) = MAX(NOISE,RGNERR(J))
          END IF
          RGNERR(J) = AREA*RGNERR(J)
          BASVAL(J) = AREA*BASVAL(J)
          RGNERR(J) = MIN(ABS(BASVAL(J)),RGNERR(J))
      END DO
      NUM = 37
      RETURN
      END SUBROUTINE Rule_T2a

END MODULE CubatureRule_T2
