! http://www.math.wsu.edu/faculty/genz/papers/simplex/node1.html
           !------------------------!
           ! Cubpack User Interface !
           !------------------------!
Module CUI

USE Precision_Model
USE internal_types

Implicit NONE

PRIVATE

PUBLIC :: CUBATR, CUBPACK_INFO
!-----------------------------------------------------------------------
!***BEGIN PROLOGUE CUBATR
!***DATE WRITTEN   901114   (YYMMDD)
!***REVISION DATE  970620   (YYMMDD)
!***REVISION DATE  980406   (YYMMDD) (MDIV removed)
!***REVISION DATE  000809   (YYMMDD)
!***REVISION DATE  010719   (YYMMDD)
!***REVISION DATE  020715   (YYMMDD) (CUBPACK_INFO added)
!***AUTHOR
!          Ronald Cools, Dept. of Computer Science,
!          Katholieke Universiteit Leuven, Celestijnenlaan 200A,
!          B-3001 Heverlee, Belgium
!          Email:  Ronald.Cools@cs.kuleuven.ac.be
!
!***PURPOSE  Computation of integrals over a collection of regions.
!
!***DESCRIPTION
!          CUBATR is the driver routine for CUBPACK and the only
!          routine that a user has to deal with (at the moment).
!
!-----------------------------------------------------------------------
PRIVATE :: CUBATR_X, CUBATR_1, CUBATR_CLEAR

!
! Module variables
!
INTEGER, SAVE, PRIVATE :: PreJob=0
INTEGER, PRIVATE                :: BOTTIH,BOTTRH
INTEGER,         DIMENSION(:), PRIVATE, ALLOCATABLE  :: IWork
REAL(kind=stnd), DIMENSION(:), PRIVATE, ALLOCATABLE  :: RWork
TYPE(EPSALG_MEM), PRIVATE              :: M

INTERFACE CUBATR
   MODULE PROCEDURE CUBATR_X, CUBATR_1, CUBATR_CLEAR
END INTERFACE

CONTAINS

SUBROUTINE CUBATR_X     &
     (DIMENS,NumFun,Integrand,NumRgn,Vertices,RgType,Value,AbsErr,    &
!   and optional parameters
IFAIL,Neval,EpsAbs,EpsRel,Restart,MinPts,MaxPts,Key,Job,Tune)
!-----------------------------------------------------------------------
!   Input parameters
!   ----------------
!
!   DIMENS Integer.
!          The dimension of the region of integration.
!
!   NumFun Integer.
!          Number of components of the integrand.
!
!   Integrand
!          Externally declared function for computing all components
!          of the integrand at the given evaluation point.
!          It must have input parameter X:
!              X(1)   The x-coordinate of the evaluation point.
!              X(2)   The y-coordinate of the evaluation point.
!              ...
!              X(DIMENS) The z-coordinate of the evaluation point.
!         and NumFun, the number of components of the integrand.
!         It must be compatible with the following interface:
!           INTERFACE 
!              FUNCTION Integrand(NUMFUN,X) RESULT(Value)
!                USE Precision_Model
!                INTEGER, INTENT(IN) :: NUMFUN
!                REAL(kind=stnd), DIMENSION(:), INTENT(IN) :: X
!                REAL(kind=stnd), DIMENSION(NUMFUN) :: Value
!              END FUNCTION Integrand
!           END INTERFACE
!
!   NumRgn Integer.
!          The number of given regions.
!
!   Vertices
!          Real array of dimension (DIMENS,DIMENS+1,NumRgn).
!          Vertices(1:DIMENS,K,L) are the x, y, ... coordinates
!          respectively of vertex K of region L, where
!          K = 1,...,DIMENS+1 and L = 1,...,NumRgn.
!
!   RgType Integer array of dimension (NumRgn).
!          RgType(L) describes the type of region L.
!
!   Value  Real array of dimension NumFun.
!          Approximations to all components of the integral if
!          the procedure is restarted.
!
!   AbsErr Real array of dimension NumFun.
!          Estimates of absolute errors if the procedure is restarted.
!
!   IFAIL  Optional integer argument.
!          This follows the NAG convention:
!          IFAIL = 1 : soft silent error
!                      Control returned to calling program.
!          IFAIL = -1: soft noisy error
!                      Error message is printed.
!                      Control returned to calling program.
!          IFAIL = 0 : hard noisy error
!                      Error message is printed and program is stopped.
!          Default IFAIL = -1.
!
!   EpsAbs Optional real argument.
!          Requested absolute error.
!          Default  EpsAbs = 0.
!
!   EpsRel Optional real argument.
!          Requested relative error.
!          Default EpsRel = sqrt(machine precision).
!
!   Restart Optional boolean argument.
!          If Restart = FALSE, this is the first attempt to compute
!                              the integral.
!          If Restart = TRUE, then we restart a previous attempt.
!          In this case the only parameters for CUBATR that may
!          be changed (with respect to the previous call of CUBATR)
!          are MinPts, MaxPts, EpsAbs, EpsRel, Key and Restart.
!          Default Restart = FALSE.
!
!   MinPts Optional integer argument.
!          The minimum allowed number of integrand evaluations.
!          Default MinPts = 0.
!
!   MaxPts Optional integer argument.
!          The maximum allowed number of integrand evaluations.
!          Default MaxPts = enough to do 500 subdivisions.
!
!   Key    Optional integer argument.
!          Can be used by Rule_General to choose between several
!          local integration rules.
!          Default Key = 2 if Dimension=1 and extrapolation is used 
!                                        (This corresponds to QAGS)
!          Default Key = 0 otherwise
!
!   Job    Optional integer argument.
!          If |Job| = 0, then nothing will be done except freeing all
!                        allocated memory.
!                        This is usefull after a call of CUBATR if no
!                        Restart will be done later and memory usage
!                        might become an issue later.
!                        Equivalently, one can call CUBATR()
!                        without any arguments.
!                   = 1, the global adaptive algorithm is called
!                   = 2, extrapolation using the epsilon algorithm is used.
!                   = 11, a region will be divided in 2**DIMENS subregions
!                        and the global adaptive algorithm is called.
!                        In combination with Key=0, this resembles DUCTRI and DCUTET.
!                   = 12, a region will be divided in 2 subregions
!                        and the global adaptive algorithm is called.
!                        In combination with Key=3 or 4, this resembles DCUHRE.
!          If Job < 0, then an overview of the Region Collection is dumped.
!          This will create the files tmp_integerstore and tmp_realstore.
!          Default Job = 1.
!
!   Tune   Optional real argument.
!          Can be used by Global_Adapt or the local error estimators
!          to influence the reliability. 0 <= Tune <= 1.
!          Tune = 1 is the most reliable available.
!          Default Tune = 1.
!          Note that this is an experimental and controversial parameter.
!          In this version, only Tune = 1 is supported for all regions.
!
!   Output parameters
!   -----------------
!
!   Value  Real array of dimension NumFun.
!          Approximations to all components of the integral
!
!   AbsErr Real array of dimension NumFun.
!          Estimates of absolute errors.
!
!   NEval  Optional Integer.
!          Number of integrand evaluations used by CUBATR for this call.
!
!   IFAIL  Optional Integer.
!          IFAIL = 0 for normal exit.
!
!            AbsErr(K) <=  EpsAbs or
!            AbsErr(K) <=  ABS(Value(K))*EpsRel with MaxPts or less
!            function evaluations for all values of K,
!            1 <= K <= NumFun .
!
!          IFAIL = 1 if MaxPts was too small to obtain the required
!            accuracy. In this case Global_Adapt returns values of
!            Value with estimated absolute errors AbsErr.
!
!          IFAIL > 1 in more serious case of trouble.
!-----------------------------------------------------------------------
! MODULES USED
USE Check_Input
USE Error_Handling
USE DS_ROUTINES, ONLY: DSCOPY, DSINIT, DSUSED, DSSTAT, DSSUM, DSPINT
USE CubatureRule_General, ONLY: Rule_Cost
USE Global_Adaptive_Algorithm
!***END PROLOGUE CUBATR
!-----------------------------------------------------------------------
!
! Global variables
!
INTERFACE 
   FUNCTION Integrand(NUMFUN,X) RESULT(Value)
     USE Precision_Model
     INTEGER, INTENT(IN) :: NUMFUN
     REAL(kind=stnd), DIMENSION(:), INTENT(IN) :: X
     REAL(kind=stnd), DIMENSION(NUMFUN) :: Value
   END FUNCTION Integrand
END INTERFACE
INTEGER, INTENT(IN)               :: DIMENS,NumFun,NumRgn
INTEGER, DIMENSION(:), INTENT(IN) :: RgType
LOGICAL, INTENT(IN),    OPTIONAL  :: Restart
INTEGER, INTENT(OUT),   OPTIONAL  :: NEval
INTEGER, INTENT(IN),    OPTIONAL  :: Job,Key,MaxPts,MinPts
INTEGER, INTENT(IN OUT), OPTIONAL  :: IFAIL
REAL(kind=stnd), INTENT(IN), OPTIONAL :: Tune,EpsAbs,EpsRel
REAL(kind=stnd), INTENT(IN), DIMENSION(:,:,:) :: Vertices
REAL(kind=stnd), INTENT(IN OUT), DIMENSION(:)  :: AbsErr, Value
!
! Named constants
!
INTEGER, PARAMETER ::  NRINFO=1, NIINFO=5
!
!  Local variables
!
INTEGER :: BLOCK,i,Inform,Leval,LJob,LMaxPts,LMinPts,MinCost, &
           NRVERT,NrSub,NRVACA,MAXRGN,RULCLS,STATUS,Tmp
LOGICAL                                     :: EpsAlg,LRestart
REAL(kind=stnd)                             :: LEpsAbs, LEpsRel
REAL(kind=stnd), DIMENSION(:), ALLOCATABLE  :: TmpRWork
INTEGER,    DIMENSION(:), ALLOCATABLE       :: TmpIWork
TYPE(INTEGRATOR_INFO)                       :: CINFO
TYPE(USER_INFO)                             :: UINFO
!-----------------------------------------------------------------------
!
! Check array sizes
!      Array size mismatch results in hard error.
   Inform = 0
   IF (size(rgtype) < numrgn) THEN
       write(unit=*,fmt=*) "Error: size(rgtype) < numrgn"
       Inform = Inform + 1
   END IF
   IF (size(abserr) < numfun) THEN
       write(unit=*,fmt=*) "Error: size(abserr) < numfun"
       Inform = Inform + 1
   END IF
   IF (size(Value) < numfun) THEN
       write(unit=*,fmt=*) "Error: size(Value) < numfun"
       Inform = Inform + 1
   END IF
   IF ((size(vertices,1) /= dimens) .or. (size(vertices,2) /= dimens+1) &
          .or. (size(vertices,3) < numrgn)) THEN
       Inform = Inform + 1
       write(unit=*,fmt=*)"Error: size(vertices) /= (/dimens,dimens+1,numrgn/)"
   END IF
   IF (Inform /= 0) THEN
      WRITE(unit=*,fmt=*)  "Array size mismatch results in hard error."
      STOP   ! "Array size mismatch results in hard error."
   END IF

!-----------------------------------------------------------------------
IF (PRESENT(NEval)) THEN
    NEval = 0
END IF
!-----------------------------------------------------------------------
IF (PRESENT(Job)) THEN
    LJob = Job
    IF (Job == 0) THEN
       CALL CUBATR_CLEAR()
       RETURN
    END IF
 ELSE
    LJob = 1 
END IF
!-----------------------------------------------------------------------
!
! Set optional arguments
!
IF ( PRESENT(Restart)) THEN 
                              LRestart = Restart
                       ELSE 
                              LRestart = .FALSE.
END IF
IF ( PRESENT(Key)) THEN 
                        CINFO%Key = Key
                   ELSE 
                       IF ((ABS(LJob) == 2) .AND. (DIMENS == 1)) THEN
                           CINFO%Key = 2   ! simulate QAGS
                       ELSE
                           CINFO%Key = 0
                       END IF
END IF
!-----------------------------------------------------------------------
!
!  Check input parameters
!
IF ( .NOT. LRestart) THEN
         CALL CHECK(DIMENS,NumFun,NumRgn,RgType,Inform, &
                    LJob,IFAIL,EpsAbs,EpsRel,MinPts,MaxPts,Tune)
ELSE IF ( ALLOCATED(IWork) ) THEN
         CALL CHECK(DIMENS,NumFun,BOTTRH,BOTTIH,IWork,  &
                    Inform,LJob,IFAIL,EpsAbs,EpsRel,MinPts,MaxPts,Tune)
ELSE
   Inform = 4096  ! There is nothing to restart from
END IF
IF (Inform /= 0) THEN 
    CALL Handle_Error(Inform,IFAIL)
    RETURN
END IF
!-----------------------------------------------------------------------
RULCLS = Rule_Cost( DIMENS, RgType(1), CINFO%Key )
MinCost = RULCLS
DO i = 2,NumRgn
   Tmp = Rule_Cost( DIMENS, RgType(i), CINFO%Key)
   RULCLS = max(RULCLS,Tmp)
   MinCost = MinCost + Tmp
END DO
!-----------------------------------------------------------------------
!
! Set optional arguments
!
   
IF ( PRESENT(MinPts)) THEN 
                              LMinPts = MinPts
                      ELSE 
                              LMinPts = 0
END IF
IF ( PRESENT(MaxPts)) THEN 
                              LMaxPts = MaxPts
                      ELSE 
                              LMaxPts = 500*RULCLS
END IF
IF ( PRESENT(Tune)) THEN 
                              CINFO%Tune = Tune
                    ELSE 
                              CINFO%Tune = 1
END IF
IF ( PRESENT(EpsAbs)) THEN 
                              LEpsAbs = EpsAbs
                      ELSE 
                              LEpsAbs = 0
END IF
IF ( PRESENT(EpsRel)) THEN 
                              LEpsRel = EpsRel
                      ELSE 
                              LEpsRel = SQRT(EPSILON(LEpsRel))
END IF
!-----------------------------------------------------------------------
!
! Set other parameters of the Global Adaptive algorithm
!      
!
!  NrSub is an upper limit for the number of subregions after subdivision.
!  This influence memory managment, so don't exagerate here.
IF (DIMENS <= 3) THEN 
                              NrSub = 2**DIMENS
                 ELSE 
                              NrSub = 4
END IF
!
EpsAlg = ( ABS(LJob) == 2 )
CINFO%UNIFORM_SUBDIV = EpsAlg
CINFO%NrSub = NrSub

IF (( ABS(LJob) == 11) .AND. (DIMENS <= 3)) THEN
   ! simulate dcutri and dcutet ; NrSub = 2**DIMENS ; EpsAlg = .FALSE.
   CINFO%UNIFORM_SUBDIV = .TRUE.
END IF

IF ( ABS(LJob) == 12 ) THEN
   ! simulate dcuhre; NrSub = 2 ; EpsAlg = .FALSE.
   CINFO%NrSub = 2
END IF

NRVERT = DIMENS + 1   ! Only cubes and simplices are implemented here.

UINFO%NumFun = NumFun  
UINFO%NumRgn = NumRgn
UINFO%MinPts = LMinPts 
UINFO%MaxPts = LMaxPts
UINFO%EpsAbs = LEpsAbs 
UINFO%EpsRel = LEpsRel
!-----------------------------------------------------------------------
IF (LRestart) THEN
  ! This requires allocating larger arrays and copying
  ! the region collection.
  ALLOCATE(TmpRWork(SIZE(RWork)),STAT=status)
  IF (status /=  0) THEN
     WRITE(unit=*,fmt=*)  "Problem allocating real workspace."
     STOP  !  "Problem allocating real workspace."
  END IF
  ALLOCATE(TmpIWork(SIZE(IWork)),STAT=status)
  IF (status /=  0) THEN
     WRITE(unit=*,fmt=*) "Problem allocating integer workspace."
     STOP  ! "Problem allocating integer workspace."
  END IF
  MAXRGN = DSUSED(IWork)
  CALL DSCOPY(IWork,RWork,TmpIWork,TmpRWork)
ELSE
  MAXRGN = NumRgn
END IF
!-----------------------------------------------------------------------
! NRVACA is the number of regions the global adaptive algorithm
!        removes from the data structure for further processing.
!        In some routines for shared memory parallel machines
!        this is the variable MDIV
NRVACA = 1
! MAXRGN depends on the number of function evalutions
MAXRGN = MAXRGN + 1 + (NrSub-1)*(LMaxPts - RULCLS*NumRgn)/(RULCLS*NrSub)
!
! Compute length of workspace needed.
!
BOTTIH = MAXRGN*(1+NIINFO) + 15 + NRVACA
BLOCK = NRINFO+NRVERT*DIMENS+2*NumFun
IF (NumFun > 1) THEN
    BLOCK = BLOCK + 1
END IF
BOTTRH = MAXRGN*BLOCK
!
! Allocate space for the region collection
!
IF (ALLOCATED(RWork)) THEN
    DEALLOCATE(RWork)
END IF
ALLOCATE(RWork(BOTTRH),STAT=status)
IF (status /=  0) THEN
   WRITE(unit=*,fmt=*) "Problem allocating real workspace."
   STOP  ! "Problem allocating real workspace."
END IF
IF (ALLOCATED(IWork)) THEN
   DEALLOCATE(IWork)
END IF
ALLOCATE(IWork(BOTTIH),STAT=status)
IF (status  /=  0) THEN
   WRITE(unit=*,fmt=*) "Problem allocating integer workspace."
   STOP  ! "Problem allocating integer workspace."
END IF
!
! Initialise region collection
!
IF ( LRestart) THEN
    CALL DSCOPY(TmpIWork,TmpRWork,IWork,RWork)
    DEALLOCATE(TmpIWork,TmpRWork)
ELSE
   IF ( MinCost > LMaxPts ) THEN
       Inform = 128    ! Dit nummer werd al gebruikt !
   ELSE
       CALL DSINIT(DIMENS,NRVERT,NIINFO,NRINFO,NumFun,NRVACA,   &
                   BOTTIH,BOTTRH,IWork,Inform)
   END IF
   IF (Inform /= 0) THEN 
       CALL Handle_Error(Inform,IFAIL)
       RETURN
   END IF
END IF
!-----------------------------------------------------------------------
!
!  Call integration routine
!
If (EpsAlg) THEN
   IF ( LRestart .AND. (PreJob /= ABS(LJob))) THEN
      Inform = 3
   ELSE
      ! Observe that only relevant array sections are passed !
      CALL Global_Adapt_Extrap(DIMENS,CINFO,UINFO,NRVERT,NIINFO,    &
                     NRINFO, Vertices(1:DIMENS,1:NRVERT,1:NUMRGN),  &
                     RgType(1:NUMRGN),Integrand,LRestart,           &
                     Value(1:NUMFUN),AbsErr(1:NUMFUN),LEval,Inform, &
                     RWork,IWork,M)
   END IF
ELSE
   IF ( LRestart .AND. (PreJob /= LJob)) THEN
       IF ( ASSOCIATED(M%RESLA)) THEN
           DEALLOCATE(M%RESLA,M%ERLARG,M%RESULT1,M%ABSERR1,M%RCOPY)
       END IF
       CALL DSPINT(IWork,RWork)
       CALL DSSUM(Value,Abserr,IWork,RWork,Inform)
   END IF
   ! Observe that only relevant array sections are passed !
   CALL Global_Adapt(DIMENS,CINFO,UINFO,NRVERT,NIINFO,NRINFO,    &
                  Vertices(1:DIMENS,1:NRVERT,1:NUMRGN),          &
                  RgType(1:NUMRGN),Integrand,LRestart,           &
                  Value(1:NUMFUN),AbsErr(1:NUMFUN),LEval,Inform, &
                  RWork,IWork)
END IF
IF (PRESENT(NEval)) THEN
   NEval = LEval
END IF
!-----------------------------------------------------------------------
 IF (LJob < 0) THEN
    WRITE(unit=*,fmt=*) "Debug mode: dumping region collection overview."
    CALL DSSTAT(IWork(:),RWork(:))          ! For debugging.
 END IF
!-----------------------------------------------------------------------
! IF ((Inform >= 8) .or. (Inform == 3)) THEN
!    Something went wrong but the data structure remains untouched
!    and so this call can be ignored.
IF ((Inform < 8) .AND. (Inform /= 3)) THEN
   PreJob = ABS(LJob)
END IF
CALL Handle_Error(Inform,IFAIL)
RETURN
END SUBROUTINE CUBATR_X


SUBROUTINE CUBATR_1                                              &
     (DIMENS,Integrand,SVertices,SRgType,SValue,SAbsErr,         &
!   and optional parameters                                      & 
      IFAIL,Neval,EpsAbs,EpsRel,Restart,MaxPts,Key,Job)
!-----------------------------------------------------------------------
!   Input parameters
!   ----------------
!
!   DIMENS Integer.
!          The dimension of the region of integration.
!
!   Integrand
!          Externally declared function for computing all components
!          of the integrand at the given evaluation point.
!          It must have input parameter X:
!              X(1)   The x-coordinate of the evaluation point.
!              X(2)   The y-coordinate of the evaluation point.
!              ...
!              X(DIMENS) The z-coordinate of the evaluation point.
!         and NumFun, the number of components of the integrand.
!         It must be compatible with the following interface:
!           INTERFACE 
!              FUNCTION Integrand(NUMFUN,X) RESULT(Value)
!                USE Precision_Model
!                INTEGER, INTENT(IN) :: NUMFUN
!                REAL(kind=stnd), DIMENSION(:), INTENT(IN) :: X
!                REAL(kind=stnd), DIMENSION(NUMFUN) :: Value
!              END FUNCTION Integrand
!           END INTERFACE
!
!   SVertices
!          Real array of dimension (DIMENS,DIMENS+1).
!          Vertices(1:DIMENS,K) are the x, y, ... coordinates
!          respectively of vertex K of the region, where
!          K = 1,...,DIMENS+1.
!
!   SRgType Integer.
!          RgType describes the type of region L.
!
!   SValue Real.
!          Approximation to the integral if the procedure is restarted.
!
!   SAbsErr Real.
!          Estimate of the absolute error if the procedure is restarted.
!
!   IFAIL  Optional integer argument.
!          This follows the NAG convention:
!          IFAIL = 1 : soft silent error
!                      Control returned to calling program.
!          IFAIL = -1: soft noisy error
!                      Error message is printed.
!                      Control returned to calling program.
!          IFAIL = 0 : hard noisy error
!                      Error message is printed and program is stopped.
!          Default IFAIL = -1.
!
!   EpsAbs Optional real argument.
!          Requested absolute error.
!          Default  EpsAbs = 0.
!
!   EpsRel Optional real argument.
!          Requested relative error.
!          Default EpsRel = sqrt(machine precision).
!
!   Restart Optional boolean argument.
!          If Restart = FALSE, this is the first attempt to compute
!                              the integral.
!          If Restart = TRUE, then we restart a previous attempt.
!          In this case the only parameters for CUBATR that may
!          be changed (with respect to the previous call of CUBATR)
!          are MinPts, MaxPts, EpsAbs, EpsRel, Key and Restart.
!          Default Restart = FALSE.
!
!   MaxPts Optional integer argument.
!          The maximum allowed number of integrand evaluations.
!          Default MaxPts = enough to do 500 subdivisions.
!
!   Key    Optional integer argument.
!          Can be used by Rule_General to choose between several
!          local integration rules.
!          Default Key = 2 if Dimension=1 and extrapolation is used 
!                                        (This corresponds to QAGS)
!          Default Key = 0 otherwise
!
!   Job    Optional integer argument.
!          If |Job| = 0, then nothing will be done except freeing all
!                        allocated memory.
!                        This is usefull after a call of CUBATR if no
!                        Restart will be done later and memory usage
!                        might become an issue later.
!                        Equivalently, one can call CUBATR()
!                        without any arguments.
!                   = 1, the global adaptive algorithm is called
!                   = 2, extrapolation using the epsilon algorithm is used.
!                   = 11, a region will be divided in 2**DIMENS subregions
!                        and the global adaptive algorithm is called.
!                        In combination with Key=0, this resembles DUCTRI and DCUTET.
!                   = 12, a region will be divided in 2 subregions
!                        and the global adaptive algorithm is called.
!                        In combination with Key=3 or 4, this resembles DCUHRE.
!          If Job < 0, then an overview of the Region Collection is dumped.
!          This will create the files tmp_integerstore and tmp_realstore.
!          Default Job = 1.
!
!   Output parameters
!   -----------------
!
!   SValue Real.
!          Approximation to the integral
!
!   AbsErr Real.
!          Estimate of the absolute error.
!
!   NEval  Optional Integer.
!          Number of integrand evaluations used by CUBATR for this call.
!
!   IFAIL  Optional Integer.
!          IFAIL = 0 for normal exit.
!
!            AbsErr(K) <=  EpsAbs or
!            AbsErr(K) <=  ABS(Value(K))*EpsRel with MaxPts or less
!            function evaluations for all values of K,
!            1 <= K <= NumFun .
!
!          IFAIL = 1 if MaxPts was too small to obtain the required
!            accuracy. In this case Global_Adapt returns values of
!            Value with estimated absolute errors AbsErr.
!
!          IFAIL > 1 in more serious case of trouble.
!-----------------------------------------------------------------------
!
! Global variables
!
INTERFACE 
   FUNCTION Integrand(NUMFUN,X) RESULT(Value)
      USE Precision_Model
      INTEGER, INTENT(IN) :: NUMFUN
      REAL(kind=stnd), DIMENSION(:), INTENT(IN) :: X
      REAL(kind=stnd), DIMENSION(NUMFUN) :: Value
   END FUNCTION Integrand
END INTERFACE
LOGICAL, OPTIONAL, INTENT(IN)                 :: Restart
INTEGER, INTENT(IN)                           :: DIMENS,SRgType
INTEGER, INTENT(OUT), OPTIONAL                :: NEval
INTEGER, INTENT(IN), OPTIONAL                 :: Key,MaxPts,Job
INTEGER, INTENT(IN OUT), OPTIONAL             :: IFAIL
REAL(kind=stnd), INTENT(IN), OPTIONAL         :: EpsAbs,EpsRel
REAL(kind=stnd), DIMENSION(:,:), INTENT(IN)   :: SVertices
REAL(kind=stnd), INTENT(IN OUT)               :: SValue,SAbsErr
!
! Local variables
!
INTEGER, DIMENSION(1)                         :: RgType
REAL(kind=stnd), DIMENSION(1)                 :: Value, AbsErr
REAL(kind=stnd), DIMENSION(DIMENS,DIMENS+1,1) :: Vertices
!-------------------
RgType(1) = SRgType
Vertices(:,:,1) = SVertices
IF (PRESENT(Restart)) THEN
    IF ( Restart ) THEN
       Value(1) = SValue
       AbsErr(1) = SAbsErr
    END IF
END IF
CALL CUBATR             &
  (DIMENS,1,Integrand,1,Vertices,RgType,Value,AbsErr,        &
   !   and optional parameters                               &
   ifail=IFAIL,neval=Neval,epsabs=EpsAbs,epsrel=EpsRel,      &
   restart=Restart,maxpts=MaxPts,key=key,job=Job)
SValue = Value(1)
SAbsErr = AbsErr(1)
RETURN
END SUBROUTINE CUBATR_1

SUBROUTINE CUBATR_CLEAR()
  IF ( ALLOCATED(Iwork) ) THEN
     DEALLOCATE(RWork,IWork)
  END IF
  IF ( ASSOCIATED(M%RESLA)) THEN
     DEALLOCATE(M%RESLA,M%ERLARG,M%RESULT1,M%ABSERR1,M%RCOPY)
  END IF
  PreJob = 0
  RETURN
END SUBROUTINE CUBATR_CLEAR

SUBROUTINE CUBPACK_INFO(uout)

integer uout

REAL(kind=stnd) :: x=1.0e-30  ! lowest accuracy of cubature formula constants

write (uout,*) " ---------------------------------------------------------------"
write (uout,*) "              CUBPACK information"
write (uout,*) "              -------------------"
write (uout,*) " The model for real numbers in the current installed version,"
write (uout,*) " obtained with the declaration REAL(KIND=stnd), has the"
write (uout,*) " following characteristics:"
write (uout,*) "                base =  ",radix(x)
write (uout,*) " digits in this base = ", digits(x)
! write (uout,*) " highest exponent = ", maxexponent(x)
! write (uout,*) " lowest exponent = "minexponent(x), "(normalized numbers)
write (uout,*) " This implies:"
write (uout,*) "            machine epsilon = ",epsilon(x)
write (uout,*) "        largest real number = ", huge(x)
write (uout,*) " smallest normalized number = ", tiny(x)
write (uout,*)
write (uout,*) " The lowest relative error that may be obtained with this"
write (uout,*) "(""  version is about "",G8.2)",max(50*epsilon(x),x)
write (uout,*) " Asking for lower error will push the routine to use the"
write (uout,*) " maximal number of function evaluations it is allowed."
write (uout,*)
write (uout,*)
write (uout,*) " This version accepts a collection of hyper-rectangles"
write (uout,*) "   (and parallelepipeds) and simplices as integration regions."
write (uout,*) " Extrapolation using the epsilon-algorithm is available"
write (uout,*) "      for dimensions 1, 2 and 3."
write (uout,*) " The following values of KEY give different integration rules:"
write (uout,*) " - finite interval: KEY = 1, 2, 3, 4, 5."
write (uout,*) "      KEY < 1 defaults to 1; KEY > 5 defaults to 5."
write (uout,*) " - n-cube:          KEY = 3, 4       uses rule of degree 2*KEY+1"
write (uout,*) "           otherwise, uses for a square a rule of degree 13"
write (uout,*) "                                 3-cube a rule of degree 11"
write (uout,*) "                                        a rule of degree  7"
write (uout,*) " - n-simplex:       KEY = 1, 2, 3, 4 uses rule of degree 2*KEY+1" 
write (uout,*) "         otherwise, uses for a triangle a rule of degree 13"
write (uout,*) "                            tetrahedron a rule of degree  8"
write (uout,*) "                                        a rule of degree  7"
write (uout,*) " KEY = 0 corresponds to our preferred choice."
write (uout,*) " ---------------------------------------------------------------"

RETURN
END SUBROUTINE CUBPACK_INFO

END Module CUI
