MODULE Global_Adaptive_Algorithm

USE Precision_Model
USE internal_types
USE Volume_Computation
USE Region_Processor

Implicit NONE
PRIVATE

PUBLIC  :: Global_Adapt, Global_Adapt_Extrap
PRIVATE :: Epsalg

CONTAINS
   SUBROUTINE Global_Adapt(DIMENS,CINFO,UINFO,NRVERT,NIINFO,NRINFO,   &
                           VERTIC,RGTYPE,Integrand,RESTART,VALUE, &
                           ABSERR,NEVAL,IFAIL,RSTORE,ISTORE)
USE DS_ROUTINES, ONLY: DSGET, DSSPUT, DSSUM, DSFREE

!***BEGIN PROLOGUE Global_Adapt
!***DATE WRITTEN   901114   (YYMMDD)
!***REVISION DATE  910503   (YYMMDD)
!***REVISION DATE  950503   (YYMMDD) (Fortran90 transformation)
!***REVISION DATE  970611   (YYMMDD) (more Fortran90)
!***REVISION DATE  980324   (YYMMDD) (MDIV removed)
!***AUTHOR
!          Ronald Cools, Dept. of Computer Science,
!          Katholieke Universiteit Leuven, Celestijnenlaan 200A,
!          B-3001 Heverlee, Belgium
!          Email:  Ronald.Cools@cs.kuleuven.ac.be
!
!***PURPOSE  Computation of integrals over a collection of regions.
!
!***DESCRIPTION
!            Global_Adapt repeatedly
!            subdivides the region with greatest estimated error
!            and estimates the integrals and the errors over the
!            new sub-regions until the error request is met or
!            MAXPTS function evaluations have been used.
!
!   Input parameters
!   ----------------
!
!   DIMENS Integer.
!          The dimension of the region of integration.
!   NRVERT Integer.
!          The number of vertices to determine a region.
!   NRSUB  Integer.
!          A region is divided into NRSUB subregions.
!   NUMFUN Integer.
!          Number of components of the integral.
!   NIINFO Integer.
!          The number of integers used to save information about
!          the region.
!          Conventions for info-record:
!          info-record(5) = 1 if there was asymptotic behaviour when the
!                           region was processed before.
!                         = 0 otherwise
!          info-record(4) = information on best direction for future division
!          info-record(3) = number of the original region where this
!                           region is a part of
!          info-record(2) = (volume of orinal region)/
!                           (volume of this region)
!          info-record(1) = type of region
!   NRINFO Integer.
!          The number of reals used to save information about
!          the region.
!   VERTIC Real array of dimension (DIMENS,NRVERT,NUMRGN).
!          VER(1,K,L), VER(2,K,L),..., VER(DIMENS,K,L) are the x, y, ...
!          coordinates respectively of vertex K of region L, where
!          K = 1,...,NRVERT and L = 1,...,NUMRGN.
!   RGTYPE Integer array of dimension (NUMRGN).
!           RGTYPE(L) describes the type of region L.
!   NUMRGN Integer.
!          The number of given regions.
!   MINPTS Integer.
!          The minimum allowed number of function evaluations.
!   MAXPTS Integer.
!          The maximum allowed number of function evaluations.
!   Integrand Externally declared function for computing
!          all components of the integrand at the given
!          evaluation point.
!          It must be compatible with the following interface:
!           INTERFACE
!              FUNCTION Integrand(NUMFUN,X)
!                 USE Precision_Model
!                 INTEGER NUMFUN
!                 REAL(stnd) X(:)
!                 REAL(stnd) Integrand(NUMFUN)
!              END
!           END INTERFACE
!          Input parameters:
!            X(1)   The x-coordinate of the evaluation point.
!            X(2)   The y-coordinate of the evaluation point.
!            ...
!            X(DIMENS) The z-coordinate of the evaluation point.
!            NUMFUN Integer that defines the number of
!                   components of I.
!
!   EPSABS Real.
!          Requested absolute error.
!   EPSREL Real.
!          Requested relative error.
!
!   RESTART Boolean.
!          If RESTART = FALSE, this is the first attempt to compute
!          the integral.
!          If RESTART = TRUE,
!          then we restart a previous attempt.
!          In this case the only parameters for Global_Adapt that may
!          be changed (with respect to the previous call of Global_Adapt)
!          are MINPTS, MAXPTS, EPSABS, EPSREL and RESTART.
!   MINPTS Integer.
!          Minimum number of integrand function evaluations.
!   MAXPTS Integer.
!          Maximum number of integrand function evaluations.
!
!   Output parameters
!   -----------------
!
!   VALUE  Real array of dimension NUMFUN.
!          Approximations to all components of the integral.
!          (It is an input parameter if RESTART=.true.)
!   ABSERR Real array of dimension NUMFUN.
!          Estimates of absolute errors.
!          (It is an input parameter if RESTART=.true.)
!   NEVAL  Integer.
!          Number of function evaluations used by Global_Adapt.
!   IFAIL  Integer.
!          IFAIL = 0 for normal exit.
!
!            ABSERR(K) <=  EPSABS or
!            ABSERR(K) <=  ABS(VALUE (K))*EPSREL with MAXPTS or less
!            function evaluations for all values of K,
!            1 <= K <= NUMFUN .
!
!          IFAIL = 1 if MAXPTS was too small for Global_Adapt
!            to obtain the required accuracy. In this case Global_Adapt
!            returns values of VALUE  with estimated absolute
!            errors ABSERR.
!
!          IFAIL = 2 if the region collection was not large enough
!            to obtain the required accuracy. In this case Global_Adapt
!            returns values of VALUE  with estimated absolute
!            errors ABSERR.
!
!          IFAIL > 10000 : Failure of the heap-maintaining routines.
!            This should never happen !
!            If IFAIL = 1000X this is IFAIL = X of DSINIT.
!
!***ROUTINES CALLED Process_Region,DSSPUT,DSGET,DSFREE
!***END PROLOGUE Global_Adapt
!
!
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
      TYPE(integrator_info), INTENT(IN) :: CINFO
      TYPE(user_info), INTENT(IN) ::       UINFO
      LOGICAL, INTENT(IN) ::    RESTART
      INTEGER, INTENT(IN) ::    DIMENS,NRVERT,NIINFO,NRINFO
      INTEGER, DIMENSION(:), INTENT(IN) ::    RGTYPE
      INTEGER, INTENT(OUT) ::   NEVAL,IFAIL
      INTEGER, DIMENSION(:), INTENT(IN OUT) :: ISTORE
      REAL(kind=stnd), DIMENSION(:,:,:), INTENT(IN) :: VERTIC
      REAL(kind=stnd), DIMENSION(:), INTENT(IN OUT):: VALUE ,ABSERR
      REAL(kind=stnd), DIMENSION(:), INTENT(IN OUT):: RSTORE
!
!   Local automatic variables.
!
      INTEGER, DIMENSION(NIINFO) :: INFOLD
      INTEGER, DIMENSION(NIINFO,CINFO%NRSUB) :: INFNEW
      REAL(kind=stnd), DIMENSION(DIMENS,NRVERT) :: VEROLD
      REAL(kind=stnd), DIMENSION(DIMENS,NRVERT,CINFO%NRSUB) :: VERNEW
      REAL(kind=stnd), DIMENSION(UINFO%NUMFUN) :: VALOLD, ERROLD
      REAL(kind=stnd), DIMENSION(UINFO%NUMFUN,CINFO%NRSUB) :: VALNEW, ERRNEW
      REAL(kind=stnd), DIMENSION(NRINFO) :: RINFOL
      REAL(kind=stnd), DIMENSION(NRINFO,CINFO%NRSUB) :: RINFNE
!
!   Local variables
!
!   NUM    Integer.
!          The number of points used by the basic rule in DRLGIN.
!   MAXSUB Integer.
!          The maximum number of regions that Process_Region may return
!   OUTSUB Integer
!          The number of regions returned by Process_Region
!
      INTEGER :: I,INFORM,L,NUM,MAXSUB,OUTSUB,NUMFUN
!
!***FIRST EXECUTABLE STATEMENT
!
      NUMFUN = UINFO%NUMFUN
      NEVAL = 0
      IF (.NOT.RESTART) THEN
          !
          !  Initialize the set of given regions.
          !  Compute estimates for integrals and errors for these regions.
          !
          VALUE  = 0
          ABSERR = 0

          DO I = 1,UINFO%NUMRGN
              !
              !  Initialise a region record
              !
              RINFOL(1) = VOLUME(DIMENS,RGTYPE(I),VERTIC(:,:,I))
              RINFOL(2:NRINFO) = 0  ! assignment to zero-sized array is legal
              INFOLD(1) = RGTYPE(I)
              INFOLD(2) = 1
              INFOLD(3) = I
              INFOLD(4:NIINFO) = 0
              !
              !  Apply the basic rule to each given region.
              !
              CALL Process_Region(DIMENS,CINFO,NRVERT,1,NUMFUN,Integrand, &
                          VERTIC(:,:,I),INFOLD,RINFOL,                 &
                          UINFO%MAXPTS-NEVAL,OUTSUB,VERNEW,            &
                          INFNEW,RINFNE,VALNEW,ERRNEW,NUM,INFORM)
              IF (INFORM /= 0) THEN
                  IFAIL = INFORM
                  RETURN
              END IF
              NEVAL = NEVAL + NUM
              !
              !  Adjust VALUE and ABSERR.
              !
              VALUE  = VALUE  + VALNEW(:,1)
              ABSERR = ABSERR + ERRNEW(:,1)
              !
              !  Store the results.
              !
              CALL DSSPUT(VERTIC(:,:,I),VALNEW(:,1),ERRNEW(:,1),       &
                          INFNEW(:,1),RINFNE(:,1),ISTORE,RSTORE,IFAIL)
              IF (IFAIL /= 0) THEN
                 RETURN
              END IF
          END DO
      END IF

!
!   Check for termination:
!   IF the number of points used is smaller than a user suplied minimum
!     OR
!   IF the estimated error of one of the approximations is too large,
!   THEN continue the proces.
!
   DO
      IF (NEVAL >= UINFO%MINPTS) THEN
          IF ( ALL( ABSERR < MAX(UINFO%EPSREL*ABS(VALUE),UINFO%EPSABS) ) ) THEN
             IFAIL = 0
             IF ( .NOT. RESTART ) THEN
                CALL DSSUM(VALUE ,ABSERR,ISTORE,RSTORE,IFAIL)
             END IF
             RETURN
          END IF
      END IF
!
!  If there is enough space to do further subdivisions
!  and it is allowed to do more function evaluations ...
!  Determine the maximum number of subregions after subdivsion
!     It should be o.k. to use dsfree()+1 but this is safer.
      MAXSUB = MIN(DSFREE(ISTORE),cinfo%NRSUB)
      IF ((MAXSUB >= 1) .AND. (UINFO%MAXPTS > NEVAL)) THEN
!
!     ... then prepare to apply the basic rule over each subregion
!         produced by dividing the region with greatest error.
!
           !
           !  Pick the region from the collection.
           !
           CALL DSGET(VEROLD,VALOLD,ERROLD,INFOLD,RINFOL,ISTORE,RSTORE,IFAIL)
           IF (IFAIL > 0) THEN
               RETURN
           END IF
           !
           !  Process the region
           !
           CALL Process_Region(DIMENS,CINFO,NRVERT,MAXSUB,NUMFUN,Integrand,&
                       VEROLD,INFOLD,RINFOL,UINFO%MAXPTS-NEVAL,OUTSUB,  &
                       VERNEW,INFNEW,RINFNE,VALNEW,ERRNEW,NUM,INFORM)
           NEVAL = NEVAL + NUM

           IF (INFORM /= 0) THEN
              !
              !  Restore data structure
              !
              CALL DSSPUT(VEROLD,VALOLD,ERROLD,INFOLD,RINFOL,ISTORE,   &
                           RSTORE,IFAIL)
              CALL DSSUM(VALUE ,ABSERR,ISTORE,RSTORE,IFAIL)
              IFAIL = INFORM
              RETURN
           END IF
           !
           !  Adjust VALUE  and ABSERR
           !
           VALUE  = VALUE  - VALOLD + SUM( VALNEW(:,1:OUTSUB), 2 )
           ABSERR = ABSERR - ERROLD + SUM( ERRNEW(:,1:OUTSUB), 2 )
           !
           !  Store the results.
           !
           DO L = 1,OUTSUB
               CALL DSSPUT(VERNEW(:,:,L),VALNEW(:,L),ERRNEW(:,L),      &
                           INFNEW(:,L),RINFNE(:,L),ISTORE,RSTORE,IFAIL)
               IF (IFAIL /= 0) THEN
                  RETURN
               END IF
           END DO
      ELSE
!
!   ... else there was not enough space available to reach the
!       requested accuracy or the maximum number of points
!       allowed is reached.
!
          CALL DSSUM(VALUE ,ABSERR,ISTORE,RSTORE,IFAIL)
          IF (UINFO%MAXPTS <= NEVAL) THEN
                                          IFAIL = 1
                                     ELSE
                                          IFAIL = 2
          END IF
          RETURN
      END IF
    END DO     ! Check for termination
    RETURN

    END SUBROUTINE Global_Adapt


SUBROUTINE Global_Adapt_Extrap(DIMENS,CINFO,UINFO,NRVERT,NIINFO,NRINFO,&
                           VERTIC,RGTYPE,Integrand,RESTART,VALUE,      &
                           ABSERR,NEVAL,IFAIL,RSTORE,ISTORE,MEM)

USE DS_ROUTINES, ONLY: DSGET, DSSPUT, DSSUM, DSFREE, DSUPUT, DSPINT

!***BEGIN PROLOGUE Global_Adapt_Extrap
!***DATE WRITTEN   9xxxxx   (YYMMDD)
!***REVISION DATE  970613   (YYMMDD)
!***REVISION DATE  980324   (YYMMDD)
!***REVISION DATE  990609   (YYMMDD)
!***AUTHOR
!          Ronald Cools, Dept. of Computer Science,
!          Katholieke Universiteit Leuven, Celestijnenlaan 200A,
!          B-3001 Heverlee, Belgium
!          Email:  Ronald.Cools@cs.kuleuven.ac.be
!
!***PURPOSE  Computation of integrals over a collection of regions.
!
!***DESCRIPTION
!            Global_Adapt_Extrap repeatedly subdivides the large
!            region with greatest estimated error, estimates the
!            integrals and the errors over the new sub-regions
!            and sometimes employs the epsilon algorithm to speed
!            up convergence, until the error request is met or
!            MAXPTS function evaluations have been used.
!
!   Input parameters
!   ----------------
!
!   DIMENS Integer.
!          The dimension of the region of integration.
!   NRVERT Integer.
!          The number of vertices to determine a region.
!   NRSUB  Integer.
!          A region is divided into NRSUB subregions.
!   NUMFUN Integer.
!          Number of components of the integral.
!   NIINFO Integer.
!          The number of integers used to save information about
!          the region.
!          Conventions for info-record:
!          info-record(5) = 1 if there was asymptotic behaviour when the
!                           region was processed before.
!                         = 0 otherwise
!          info-record(4) = information on best direction for future division
!          info-record(3) = number of the original region where this
!                           region is a part of
!          info-record(2) = (volume of original region)/
!                           (volume of this region)
!          info-record(1) = type of region
!   NRINFO Integer.
!          The number of reals used to save information about
!          the region.
!   VERTIC Real array of dimension (DIMENS,NRVERT,NUMRGN).
!          VER(1,K,L), VER(2,K,L),..., VER(DIMENS,K,L) are the x, y, ...
!          coordinates respectively of vertex K of region L, where
!          K = 1,...,NRVERT and L = 1,...,NUMRGN.
!   RGTYPE Integer array of dimension (NUMRGN).
!           RGTYPE(L) describes the type of region L.
!   NUMRGN Integer.
!          The number of given regions.
!   MINPTS Integer.
!          The minimum allowed number of function evaluations.
!   MAXPTS Integer.
!          The maximum allowed number of function evaluations.
!   Integrand Externally declared function for computing
!          all components of the integrand at the given
!          evaluation point.
!          It must be compatible with the following interface:
!           INTERFACE
!              FUNCTION Integrand(NUMFUN,X)
!                 USE Precision_Model
!                 INTEGER NUMFUN
!                 REAL(kind=stnd) X(:)
!                 REAL(kind=stnd) Integrand(NUMFUN)
!              END
!           END INTERFACE
!          Input parameters:
!            X(1)   The x-coordinate of the evaluation point.
!            X(2)   The y-coordinate of the evaluation point.
!            ...
!            X(DIMENS) The z-coordinate of the evaluation point.
!            NUMFUN Integer that defines the number of
!                   components of I.
!
!   EPSABS Real.
!          Requested absolute error.
!   EPSREL Real.
!          Requested relative error.
!
!   RESTART Boolean.
!          If RESTART = FALSE, this is the first attempt to compute
!          the integral.
!          If RESTART = TRUE,
!          then we restart a previous attempt.
!          In this case the only parameters for Global_Adapt_Extrap that may
!          be changed (with respect to the previous call of Global_Adapt_Extrap)
!          are MINPTS, MAXPTS, EPSABS, EPSREL and RESTART.
!   MINPTS Integer.
!          Minimum number of integrand function evaluations.
!   MAXPTS Integer.
!          Maximum number of integrand function evaluations.
!
!   Output parameters
!   -----------------
!
!   VALUE  Real array of dimension NUMFUN.
!          Approximations to all components of the integral.
!          (It is an input parameter if RESTART=.true.)
!   ABSERR Real array of dimension NUMFUN.
!          Estimates of absolute errors.
!          (It is an input parameter if RESTART=.true.)
!   NEVAL  Integer.
!          Number of function evaluations used by Global_Adapt_Extrap.
!   IFAIL  Integer.
!          IFAIL = 0 for normal exit.
!
!            ABSERR(K) <=  EPSABS or
!            ABSERR(K) <=  ABS(VALUE (K))*EPSREL with MAXPTS or less
!            function evaluations for all values of K,
!            1 <= K <= NUMFUN .
!
!          IFAIL = 1 if MAXPTS was too small for Global_Adapt_Extrap
!            to obtain the required accuracy. In this case Global_Adapt_Extrap
!            returns values of VALUE  with estimated absolute
!            errors ABSERR.
!
!          IFAIL = 2 if the region collection was not large enough
!            to obtain the required accuracy. In this case Global_Adapt_Extrap
!            returns values of VALUE  with estimated absolute
!            errors ABSERR.
!
!          IFAIL > 10000 : Failure of the heap-maintaining routines.
!            This should never happen !
!            If IFAIL = 1000X this is IFAIL = X of DSINIT.
!
!***ROUTINES CALLED Process_Region,DSSPUT,DSGET,DSFREE
!***END PROLOGUE Global_Adapt_Extrap
!
!
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
      TYPE(integrator_info), INTENT(IN) :: CINFO
      TYPE(user_info), INTENT(IN) ::       UINFO
      TYPE(epsalg_mem), INTENT(IN OUT) ::   MEM
      LOGICAL, INTENT(IN) ::    RESTART
      INTEGER, INTENT(IN) ::    DIMENS,NRVERT,NIINFO,NRINFO
      INTEGER, DIMENSION(:), INTENT(IN) ::    RGTYPE
      INTEGER, INTENT(OUT) ::   NEVAL,IFAIL
      INTEGER, DIMENSION(:), INTENT(IN OUT) :: ISTORE
      REAL(kind=stnd), DIMENSION(:,:,:), INTENT(IN) :: VERTIC
      REAL(kind=stnd), DIMENSION(:), INTENT(IN OUT):: VALUE , ABSERR
      REAL(kind=stnd), DIMENSION(:), INTENT(IN OUT):: RSTORE
!
!   Local automatic variables.
!
      INTEGER, DIMENSION(NIINFO) :: INFOLD
      INTEGER, DIMENSION(NIINFO,CINFO%NRSUB) :: INFNEW
      REAL(kind=stnd), DIMENSION(DIMENS,NRVERT) :: VEROLD
      REAL(kind=stnd), DIMENSION(DIMENS,NRVERT,CINFO%NRSUB) :: VERNEW
      REAL(kind=stnd), DIMENSION(UINFO%NUMFUN) :: VALOLD, ERROLD, &
                                     BNDEC, ERLAST, ERRO4, ERRBND

      REAL(kind=stnd), DIMENSION(UINFO%NUMFUN,CINFO%NRSUB) ::  VALNEW,ERRNEW
      REAL(kind=stnd), DIMENSION(NRINFO) :: RINFOL
      REAL(kind=stnd), DIMENSION(NRINFO,CINFO%NRSUB) ::  RINFNE
!
!   Local variables
!
!   NUM    Integer.
!          The number of points used by the basic rule in DRLGIN.
!   MAXSUB Integer.
!          The maximum number of regions that Process_Region may return
!   OUTSUB Integer
!          The number of regions returned by Process_Region
!   NRRCOPY Integer array.
!          Number of elements in RCOPY
!   ERLARG Real array
!          The sum of the error-estimates of the 'big' regions
!   MAXERRORPOOL Real
!          The maximum error over a 'small' region
!
      INTEGER :: I,INFORM,J,L,OUTSUB,STATUS
      INTEGER :: NUM,MAXSUB,NUMFUN,NRES
      LOGICAL :: ready2extrap
      REAL(kind=stnd) :: MAXHELP
      REAL(kind=stnd), PARAMETER :: EPRN = 5.0_stnd*EPSILON(MAXHELP)
!
!***FIRST EXECUTABLE STATEMENT
!
      NUMFUN = UINFO%NUMFUN
      NEVAL = 0
      NRES = 0
      IF (.NOT.RESTART) THEN
          !
          !  Initialize the set of given regions.
          !  Compute estimates for integrals and errors for these regions.
          !
          IF ( ASSOCIATED(MEM%RCOPY) ) THEN
             DEALLOCATE( MEM%RCOPY,MEM%NRRCOPY,MEM%RESLA,MEM%ERLARG, &
                         MEM%RESULT1,MEM%ABSERR1 )
          END IF
          VALUE = 0
          ABSERR = 0

          DO I = 1,UINFO%NUMRGN
              !
              !  Initialise a region record
              !
              RINFOL(1) = VOLUME(DIMENS,RGTYPE(I),VERTIC(:,:,I))
              RINFOL(2:NRINFO) = 0
              INFOLD(1) = RGTYPE(I)
              INFOLD(2) = 1
              INFOLD(3) = I
              INFOLD(4:NIINFO) = 0
              !
              !  Apply the basic rule to each given region.
              !
              CALL Process_Region(DIMENS,CINFO,NRVERT,1,NUMFUN,Integrand, &
                          VERTIC(:,:,I),INFOLD,RINFOL,UINFO%MAXPTS-NEVAL, &
                          OUTSUB,VERNEW,INFNEW,RINFNE,VALNEW,ERRNEW,NUM,INFORM)
              IF (INFORM /= 0) THEN
                  IFAIL = INFORM
                  RETURN
              END IF
              NEVAL = NEVAL + NUM
              !
              !  Adjust VALUE and ABSERR.
              !
              VALUE = VALUE + VALNEW(:,1)
              ABSERR = ABSERR + ERRNEW(:,1)
              !
              !  Store the results.
              !
              CALL DSSPUT(VERTIC(:,:,I),VALNEW(:,1),ERRNEW(:,1),      &
                          INFNEW(:,1),RINFNE(:,1),ISTORE,RSTORE,IFAIL)
              IF (IFAIL /= 0) THEN
                 write(unit=*,fmt=*) "DSSPUT error 4"
                 STOP
              END IF
          END DO
      END IF
      IF (NEVAL >= UINFO%MINPTS) THEN
           IF ( ALL( ABSERR < MAX(UINFO%EPSREL*ABS(VALUE),UINFO%EPSABS) ) ) THEN
              IFAIL = 0
              IF ( .NOT. RESTART ) THEN
                  CALL DSSUM(VALUE,ABSERR,ISTORE,RSTORE,IFAIL)
              END IF
              RETURN
           END IF
       END IF

IF ( ( .NOT. RESTART ) .OR. ( RESTART .AND. ( .NOT. ASSOCIATED(MEM%RCOPY) ) ) ) THEN
   ALLOCATE( MEM%RCOPY(EPSTABLENGHT,UINFO%NUMFUN),MEM%RESLA(UINFO%NUMFUN,3), &
             MEM%ERLARG(UINFO%NUMFUN), MEM%NRRCOPY(UINFO%NUMFUN),  &
             MEM%RESULT1(UINFO%NUMFUN),MEM%ABSERR1(UINFO%NUMFUN),  &
             STAT=status )
   IF (status /=  0) THEN
           write(unit=*,fmt=*) "Problem allocating real workspace."
           STOP
   END IF
   MEM%RESLA = 0                !  RC 23-7-2001
   MEM%NRRCOPY = 1
   MEM%DIVLEVEL = 2
   MEM%ERRORMAXPOOL = 0
   MEM%RCOPY (1,:) = VALUE
   MEM%ERLARG = ABSERR
   MEM%HEURISTIC_USED = .FALSE.
   MEM%EPSABS = UINFO%EPSABS
   MEM%EPSREL = UINFO%EPSREL
ELSE
   CALL DSSUM(VALUE,ABSERR,ISTORE,RSTORE,IFAIL)
   IF ( ((UINFO%EPSABS < MEM%EPSABS) .OR. (UINFO%EPSREL < MEM%EPSREL)) &
        .AND. MEM%HEURISTIC_USED ) THEN
      !
      ! Restarting with higher precision requests might confuse
      ! the heuristic. Hence a new extrapolation table is started!
      ! This is reliable, but usually very expensive.
      ! The advice is NOT to use this!
      !
      MEM%NRRCOPY = 1
      MEM%DIVLEVEL = 2
      MEM%ERRORMAXPOOL = 0      ! This is redundant
      MEM%RCOPY (1,:) = VALUE
      MEM%ERLARG = ABSERR
      MEM%HEURISTIC_USED = .FALSE.
      MEM%EPSABS = UINFO%EPSABS
      MEM%EPSREL = UINFO%EPSREL
   END IF
END IF
!
! End of preparation. Now we can really start
!
CALL DSGET(VEROLD,VALOLD,ERROLD,INFOLD,RINFOL,ISTORE,RSTORE,IFAIL)
ready2extrap = .false.
DO
!
!  If there is enough space to do further subdivisions
!  and it is allowed to do more function evaluations ...
!  Determine the maximum number of subregions after subdivsion
!     It should be o.k. to use dsfree()+1 but this is safer.
!
      MAXSUB = MIN(DSFREE(ISTORE),cinfo%NRSUB)
      IF ( (CINFO%UNIFORM_SUBDIV .AND. (MAXSUB /= cinfo%NRSUB)) .OR.  &
           ( MAXSUB < 2) .OR. (UINFO%MAXPTS <= NEVAL)) THEN
!
!   ... there was not enough space available to reach the
!       requested accuracy or the maximum number of points
!       allowed is reached.
!
          CALL DSSPUT(VEROLD,VALOLD,ERROLD,INFOLD,RINFOL,ISTORE,RSTORE,IFAIL)
          CALL DSSUM(VALUE,ABSERR,ISTORE,RSTORE,IFAIL)
          IF ( ALL(MEM%NRRCOPY > 3) ) THEN   ! Why 3 in this statement ?
             IF ( ALL(MEM%ABSERR1 < ABSERR) ) THEN  ! Refine
                ABSERR = MEM%ABSERR1
                VALUE = MEM%RESULT1
             END IF
          END IF
          IFAIL = 1
          RETURN
      ELSE
!     ... then prepare to apply the basic rule over each subregion
!         produced by dividing the region with greatest errors.
!
          I = 1  ! While no parallel computing is implemented yet
          IF ( .NOT. ready2extrap ) THEN
             IF (MAXVAL(ABS(ERROLD)) < MEM%ERRORMAXPOOL ) THEN
                !
                !  Process the region
                !
                DO J = 1, UINFO%NUMFUN
                   ERRBND(J) = MAX(UINFO%EPSABS,UINFO%EPSREL*ABS(VALUE(J)),EPRN*ABS(VALUE(J)))
                   BNDEC(J) = MAX( EPRN*(ABS(VALUE(J))),       &
                                   MIN(0.1_stnd*ERRBND(J),0.001_stnd*ABS(VALUE(J))))
                END DO
                ! The following test is the reason why restarting is so difficult.
                ready2extrap = ALL(MEM%ERLARG < BNDEC)
                IF ( ready2extrap ) THEN
                   !
                   ! Ready for an extrapolation step. Restore first.
                   !
                   MEM%HEURISTIC_USED = .TRUE.
                   CALL DSSPUT(VEROLD,VALOLD,ERROLD,     &
                            INFOLD,RINFOL,ISTORE,RSTORE,IFAIL)
                   IF ( IFAIL /= 0 ) THEN
                       write(unit=*,fmt=*) "DSSPUT-1 error ",ifail
                       STOP
                   END IF
                END IF
             END IF
             IF ( .NOT. ready2extrap ) THEN
                CALL Process_Region(DIMENS,CINFO,NRVERT,MAXSUB,NUMFUN,      &
                            Integrand,VEROLD,INFOLD,         &
                            RINFOL,UINFO%MAXPTS-NEVAL,OUTSUB,       &
                            VERNEW,INFNEW,RINFNE,VALNEW,ERRNEW,NUM,INFORM)
                IF (INFORM /= 0) THEN
                   !
                   !  Restore data structure
                   !
                   CALL DSSPUT(VEROLD,VALOLD,ERROLD,INFOLD,RINFOL, & 
                               ISTORE,RSTORE,IFAIL)
                   IF (IFAIL /= 0) THEN
                      write(unit=*,fmt=*) "DSSPUT error 2"
                      STOP
                   END IF
                   CALL DSSUM(VALUE,ABSERR,ISTORE,RSTORE,IFAIL)
                   IF (((ALL(MEM%NRRCOPY > 3)) .AND. (DIMENS > 1) ) .OR.  &
                       ((ALL(MEM%NRRCOPY > 5)) .AND. (DIMENS == 1))) THEN
                      IF ( ALL(MEM%ABSERR1 < ABSERR) ) THEN ! Refine
                         ABSERR = MEM%ABSERR1
                         VALUE = MEM%RESULT1
                      END IF
                   END IF
                   IFAIL = INFORM
                   RETURN
                END IF
                NEVAL = NEVAL + NUM
                ERRO4 = 0
                DO J = 1, OUTSUB
                   ERRO4 = ERRO4 + ERRNEW(:,J)
                END DO
                ERLAST = ERROLD
                MEM%ERLARG = MEM%ERLARG - ERLAST
                IF (INFNEW(2,1) < MEM%DIVLEVEL) THEN
                   MEM%ERLARG = MEM%ERLARG + ERRO4
                END IF
                !
                !  Adjust VALUE and ABSERR
                !
                VALUE = VALUE - VALOLD
                ABSERR = ABSERR - ERROLD
                DO L = 1,OUTSUB
                   VALUE = VALUE + VALNEW(:,L)
                   ABSERR = ABSERR + ERRNEW(:,L)
                END DO
                !
                !  Store the results.
                !
                DO L = 1,OUTSUB
                   IF ( INFNEW(2,1) >= MEM%DIVLEVEL ) THEN
                      MAXHELP = MAXVAL(ERRNEW(:,L))
                      MEM%ERRORMAXPOOL = MAX(MEM%ERRORMAXPOOL,MAXHELP)
                      CALL DSUPUT(VERNEW(:,:,L),VALNEW(:,L),          &
                                  ERRNEW(:,L),INFNEW(:,L),            &
                                  RINFNE(:,L),ISTORE,RSTORE,IFAIL)
                      IF (IFAIL /= 0) THEN
                        WRITE(unit=*,fmt=*) "DSUPUT-1 error ", ifail
                        STOP
                      END IF
                   ELSE
                      CALL DSSPUT(VERNEW(:,:,L),VALNEW(:,L),          &
                                  ERRNEW(:,L),INFNEW(:,L),            &
                                  RINFNE(:,L),ISTORE,RSTORE,IFAIL)
                      IF (IFAIL /= 0) THEN
                         WRITE(unit=*,fmt=*)  "DSSPUT error 3"
                         STOP
                      END IF
                   END IF
                END DO

                ! HIER EVENTUEEL OOK OP TERMINATIE TESTEN

                CALL DSGET(VEROLD,VALOLD,ERROLD,         &
                           INFOLD,RINFOL,ISTORE,RSTORE,IFAIL)
                ready2extrap = (IFAIL /= 0) ! sorted list is empty
             END IF
          END IF
          IF ( ready2extrap ) THEN
             MEM%NRRCOPY = MEM%NRRCOPY + 1
             DO J = 1,NUMFUN
                MEM%RCOPY (MEM%NRRCOPY(J),J) = VALUE(J)
             END DO
             CALL DSPINT(ISTORE,RSTORE)
             MEM%ERRORMAXPOOL = 0
             MEM%DIVLEVEL = MEM%DIVLEVEL + 1
             IF (ALL(MEM%NRRCOPY == 2)) THEN   !Dit is niet echt juist
                MEM%RESLA(:,3) = VALUE
                MEM%ERLARG = ABSERR
             ELSE
                DO J = 1,NUMFUN
                   CALL EPSALG(MEM%NRRCOPY(J),MEM%RCOPY(:,J),MEM%RESULT1(J),  &
                               MEM%ABSERR1(J),MEM%RESLA(J,:),NRES,DIMENS)
                END DO
                MEM%ABSERR1 = MEM%ABSERR1 + MEM%ERLARG
                MEM%ERLARG = ABSERR
                IF ( .NOT. ALL( (ABSERR-MEM%ABSERR1) > 0 )) THEN
                   IF (ALL((UINFO%EPSREL*ABS(VALUE)-ABSERR) > 0) .OR.  &
                       (ALL((UINFO%EPSABS-ABSERR) > 0 )) ) THEN
                      IFAIL = 0
                      RETURN
                   END IF
                ELSE
                   IF (ALL((UINFO%EPSREL*ABS(MEM%RESULT1)-MEM%ABSERR1)>0) &
                           .OR. (ALL((UINFO%EPSABS-MEM%ABSERR1)>0))) THEN
                      ABSERR = MEM%ABSERR1
                      VALUE = MEM%RESULT1
                      IFAIL = 0
                      RETURN
                   END IF
                END IF
             END IF
             CALL DSGET(VEROLD,VALOLD,ERROLD,INFOLD,RINFOL,ISTORE,RSTORE,IFAIL)
             ready2extrap = (IFAIL /= 0) ! sorted list is empty
          END IF
      END IF
   END DO
   RETURN

   END SUBROUTINE Global_Adapt_Extrap

!-----------------------------------------------------------------------


SUBROUTINE Epsalg(n, epstab, value, abserr, res3last, nres, dimens)
!
!-----------------------------------------------------------------------
!***DATE WRITTEN 961024   (YYMMDD)
!***REVISION DATE  990610   (YYMMDD) (Init. value added)
!***AUTHOR
!***BEGIN PROLOGUE epsalg
!***PURPOSE
!            The routine transforms a given sequence of approximations
!            using the epsilon algorithm of P. Wynn.
!            An estimate of the absolute error is also given.
!            the condensed epsilon table is computed. Only those
!            elements needed for the computation of the next diagonal
!            are preserved.
!***DESCRIPTION
!   ON ENTRY
!   N      integer
!          epstab(n) contains the new element in the
!          first column of the epsilon table.
!
!   EPSTAB real ( stnd ) one dimensional array of dimension epstablenght
!          containing the elements of two lower diagonals of the
!          triangular epsilon table. the elements are numbered
!          starting at the right-hand corner of the triangle.
!          the dimension should be at least (limexp+2).
!
!   RES3LAST  real ( stnd ) one dimensional array
!          previous result
!
!   NRES   integer ( only used if DIMENS==1 )
!          number of calls to the routine
!          (should be zero at first call)
!
!   DIMENS integer
!          the dimension of the integration problem
!
!   ON RETURN
!
!   VALUE real ( stnd )
!          resulting approximation to the integral
!
!   ABSERR real ( stnd )
!          estimate of the absolute error
!
!***REFERENCES
!          Algorithm 612 TRIEX:
!          Integration over a triangle using nonlinear extrapolation,
!          E. de Doncker, ACM TOMS, Vol 10, No. 1, March 1984, Pages 17-22
!***ROUTINES CALLED
!***END PROLOGUE Epsalg
!-----------------------------------------------------------------------
!
!           EPMACH   the largest relative space
!           OFLOW    the largest positive magnitude
      integer, intent( in ) :: dimens
      integer, intent( in out ) :: n
!        Changing n will give problems if NUMFUN > 1
      integer, intent( in out ) :: nres
      real (kind=stnd ), dimension(:), intent( in out ) :: res3last
      real (kind=stnd ), intent( out ) :: value, abserr
      real (kind=stnd ), dimension (:), intent( in out ) :: epstab

      real (kind=stnd ), parameter :: one = 1.0_stnd
      real (kind=stnd ), parameter :: epmach = EPSILON ( one ),      &
                                      oflow  = HUGE ( one ),         &
                                      five = 5.0_stnd
      ! The following constant is derived from one defined in
      ! module internal_types:
      integer, parameter           :: limexp = EPSTABLENGHT - 2

      real (kind=stnd ) :: delta1, delta2, delta3, epsinf, error,  &
          err1, err2, err3, e0, e1, e2, e3, e1abs, res, ss, tol1,  &
          tol2, tol3
      integer :: i, ib, ib2, k1, newelm, num
!
!           limexp  is the maximum number of elements the epsilon
!           table can contain. if this number is reached, the upper
!           diagonal of the epsilon table is deleted.
!           ( epstab  is of dimension (limexp+2) at least.)
!
!           list of major variables
!           -----------------------
!
!           E0       the 4 elements on which the
!           E1       computation of a new element in
!           E2       the epsilon table is based
!           E3                 E0
!                        E3    E1    NEW
!                              E2
!           NEWELM   number of elements to be computed in the new
!                    diagonal
!           ERROR    error = abs(e0-e1)+abs(e1-e2)+abs(e2-new)
!           VALUE   the element in the new diagonal with least error
!           NUM     a copy of the original value of N
!
!***FIRST EXECUTABLE STATEMENT
      abserr = oflow
      value = epstab(n)
      IF ( dimens == 1 ) THEN
         !
         ! Since our implementation of epsalg is taken from Triex, for dimens==1
         ! some modifications have to be made to simulate dqelg from
         ! quadpack/dqags
         !
         nres = nres + 1
         IF ( n < 3 ) THEN
            abserr = max(abserr, five*epmach*abs(value))
            return
         END IF
      END IF
      epstab( n + 2 ) = epstab( n )
      epstab( n ) = oflow
      newelm = ( n - 1 ) / 2
      num = n
      k1 = n
      do i = 1, newelm
         res = epstab( k1 + 2 )
         e2 = res
         e1 = epstab( k1 - 1 )
         e0 = epstab( k1 - 2 )
         e1abs = abs( e1 )
         delta2 = e2 - e1
         err2 = abs( delta2 )
         tol2 = max( abs( e2 ), e1abs ) * epmach
         delta3 = e1 - e0
         err3 = abs( delta3 )
         tol3 = max( e1abs, abs( e0 ) ) * epmach
         if ( .NOT. (err2 > tol2 .OR. err3 > tol3 ) ) then
!
!              if e0, e1 and e2 are equal to within machine
!              accuracy, convergence is assumed.
!              value = e2
!              abserr = abs(e1-e0)+abs(e2-e1)
!
            value = res
            abserr = err2 + err3
            abserr = max ( abserr, five * epmach * abs ( value ) )
            return
         else
            e3 = epstab( k1 )
            epstab( k1 ) = e1
            delta1 = e1 - e3
            err1 = abs( delta1 )
            tol1 = max( e1abs, abs( e3 ) ) * epmach
!
!              if two elements are very close to each other, omit
!              a part of the table by adjusting the value of n.
!
            if ( err1 <= tol1 .OR. err2 <= tol2 .OR. err3 <= tol3 ) then
               n = i + i - 1
               exit
            else
               ss = one / delta1 + one / delta2 - one / delta3
               epsinf = abs( ss * e1 )
!
!                 test to detect irregular behaviour in the table, and
!                 eventually omit a part of the table adjusting the value
!                 of n.
!
               if ( epsinf > 0.1e-03_stnd ) then
!
!                 compute a new element and eventually adjust
!                 the value of value
!
                  res = e1 + one / ss
                  epstab( k1 ) = res
                  k1 = k1 - 2
                  error = err2 + abs( res - e2 ) + err3
                  if ( .NOT. ( error > abserr ) )then
                     abserr = error
                     value = res
                  end if
               else
                  n = i + i - 1
                  exit
               end if
            end if
         end if
      end do
!
!           shift the table
!
      if (n == limexp) then
          n = 2 * ( limexp / 2 ) - 1
      end if
      if ( modulo(num,2) == 0 ) then
         ib = 2
      else
         ib = 1
      end if
      do i = 1, newelm + 1
         ib2 = ib + 2
         epstab ( ib ) = epstab ( ib2 )
         ib = ib2
      end do

      if ( num /= n ) then
         epstab(1:n) = epstab(num-n+1:num)
      end if

      SELECT CASE (DIMENS)
      CASE (1)

      if ( nres < 4 ) then
         res3last(nres) = value
         abserr = oflow
         abserr = max( abserr, five*epmach*abs(value) )
         RETURN
      end if
!
!           compute error estimate
!
      abserr = abs(value-res3last(3))+abs(value-res3last(2))+abs(value-res3last(1))

      CASE DEFAULT

      !
      !     compute error estimate
      !
      abserr = abs ( value - res3last(3) )

      END SELECT
      abserr = max ( abserr, five * epmach * abs ( value ) )
      res3last(1:2) = res3last(2:3)
      res3last(3) = value
      RETURN
   END SUBROUTINE Epsalg

END MODULE Global_Adaptive_Algorithm
