! This file is F-compatible, except for upper/lower case conventions.
!--------------------------------------------------------------------
MODULE Subdivisions

USE Precision_Model, ONLY: stnd
USE internal_types

Implicit NONE 

PRIVATE

PUBLIC :: DIVIDE
PRIVATE :: DIVSMP

CONTAINS
      SUBROUTINE DIVIDE(DIMENS,NRVERT,MAXSUB,UNIFORM_SUBDIV,NUMFUN,VEROLD,   &
            INFOLD,RINFOL,Integrand,OUTSUB,NUM,VERNEW,INFNEW,RINFNE,IFAIL)
!***BEGIN PROLOGUE DIVIDE
!***DATE WRITTEN   900615  (YYMMDD)
!***REVISION DATE  910506  (YYMMDD)
!***REVISION DATE  970401  (YYMMDD) (subdivision modifications)
!***REVISION DATE  980331  (YYMMDD) (1D added)
!***REVISION DATE  980406  (YYMMDD) (2div for nCube added)
!***REVISION DATE  980408  (YYMMDD) (F conversion)
!***REVISION DATE  990525  (YYMMDD) (re-organising)
!***REVISION DATE  990602  (YYMMDD) (2/4-division for T2 activated)
!***REVISION DATE  990624  (YYMMDD) (2/4/8-division for C3 activated)
!***REVISION DATE  010919  (YYMMDD) (infold(4) code changed elsewhere)
!***REVISION DATE  020716  (YYMMDD) (replace MOD intrinsic by MODULO)
!***AUTHOR
!          Ronald Cools, Dept. of Computer Science,
!          Katholieke Universiteit Leuven, Celestijnenlaan 200A,
!          B-3001 Heverlee, Belgium
!          Email: Ronald.Cools@cs.kuleuven.ac.be
!
!          Alan Genz
!          Department of Mathematics
!          Washington State University
!          Pullman, WA 99164-3113, USA
!          Email: AlanGenz@wsu.edu

!
!***PURPOSE To divide a given region in OUTSUB subregions
!           with equal volume.
!
!***DESCRIPTION
!   Input parameters
!   ----------------
!   DIMENS Integer, dimension of the regions
!   NRVERT Integer, number of vertices that describe a region
!   MAXSUB Integer.
!          The given region is divided in at most MAXSUB subregions
!   UNIFORM_SUBDIV Logical.
!          If true, this routine does a 2**dim subdivision.
!          If false, this routine may decide to divide into less than
!          MAXSUB regions
!   NUMFUN Integer, number of components of the integral.
!   VEROLD Real array of dimension (dimens,nrvert)
!          Contains the vertices of the given region
!   INFOLD
!   RINFOL
!   Integrand 
!          Externally declared function for computing all components
!          of the integrand at the given evaluation point.
!          It must have input parameter X:
!              X(1)      The x-coordinate of the evaluation point.
!              X(2)      The y-coordinate of the evaluation point.
!              ...
!              X(DIMENS) The z-coordinate of the evaluation point.
!         and NumFun, the number of components of the integrand.
!         It must be compatible with the following interface:
!           INTERFACE 
!              FUNCTION Integrand(NumFun,X)
!                 USE Precision_Model
!                 INTEGER :: NumFun
!                 REAL(kind=stnd) :: X(:)
!                 REAL(kind=stnd) :: Integrand(NumFun)
!              END 
!           END INTERFACE
!
!   Output parameters
!   -----------------
!   NUM    : Integer number of integrand values used to decide subdivision.
!   OUTSUB : Integer number of subregions the given region was divided into.
!   VERNEW : Real array of dimension (dimens,nrvert,MAXSUB)
!            Contains the vertices of the subregions.
!   INFNEW : Integer array of dimension (:,MAXSUM)
!            Contains additional information for each subregion.
!   RINFNE : Real array of dimension (:,MAXSUM)
!            Contains additional information for each subregion.
!   IFAIL  : integer to indicate success or failure
!            IFAIL = 0 on normal exit
!            IFAIL = 7 if the desired subdivision is not implemented
!                       for this type of region
!            IFAIL = 6 if no subdivsions are implemented for this type
!                       of region
!
!***ROUTINES CALLED  DIVSMP
!***END PROLOGUE DIVIDE
!
!  Global variables
      INTEGER, INTENT(IN)    :: DIMENS,NRVERT,MAXSUB,NUMFUN
      INTEGER, DIMENSION(:), INTENT(IN) :: INFOLD
      LOGICAL, INTENT(IN)    :: UNIFORM_SUBDIV
      REAL(kind=stnd), DIMENSION(:,:), INTENT(IN) :: VEROLD
      REAL(kind=stnd), DIMENSION(:), INTENT(IN) :: RINFOL
      INTEGER, INTENT(OUT)   :: NUM,OUTSUB,IFAIL
      INTEGER, DIMENSION(:,:), INTENT(OUT) :: INFNEW
      REAL(kind=stnd), DIMENSION(:,:,:), INTENT(OUT):: VERNEW
      REAL(kind=stnd), DIMENSION(:,:), INTENT(OUT) :: RINFNE
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
!   GEOMETRY : specifies the type of region
!
      INTEGER    :: I,J,GEOMETRY
      REAL(kind=stnd) :: VOLUME
      REAL(kind=stnd), DIMENSION(3) :: HALF_HEIGHT, HEIGHT
!
!***FIRST EXEUTABLE STATEMENT
      IF (MAXSUB == 1) THEN
            VERNEW(:,:,1) = VEROLD
            INFNEW(:,1) = INFOLD
            RINFNE(:,1) = RINFOL
            OUTSUB = 1
            NUM = 0
            IFAIL = 0
            RETURN
      END IF
!
      GEOMETRY = INFOLD(1)
      SELECT CASE (GEOMETRY)
      CASE (Simplex)
          IF (DIMENS == 1 ) THEN
              VERNEW(:,:,1) = VEROLD
              VERNEW(:,:,2) = VEROLD
              VERNEW(:,2,1) = (VEROLD(1,1)+VEROLD(1,2))/2
              VERNEW(:,1,2) = VERNEW(:,2,1)
              OUTSUB = 2
              NUM = 0
!
          ELSE IF ((DIMENS == 2) .AND. UNIFORM_SUBDIV ) THEN
              VERNEW(:,:,1) = VEROLD
              VERNEW(:,:,2) = VEROLD
              VERNEW(:,:,3) = VEROLD
              VERNEW(:,:,4) = VEROLD
              VERNEW(:,2,1) = ( VEROLD(:,1) + VEROLD(:,2) )/2
              VERNEW(:,3,1) = ( VEROLD(:,1) + VEROLD(:,3) )/2
              VERNEW(:,1,2) = VERNEW(:,2,1)
              VERNEW(:,3,2) = ( VEROLD(:,2) + VEROLD(:,3) )/2
              VERNEW(:,1,3) = VERNEW(:,3,1)
              VERNEW(:,2,3) = VERNEW(:,3,2)
              VERNEW(:,1,4) = VERNEW(:,3,2)
              VERNEW(:,2,4) = VERNEW(:,2,1)
              VERNEW(:,3,4) = VERNEW(:,3,1)      
              OUTSUB = 4
              NUM = 0
!
          ELSE IF ((DIMENS == 3) .AND. UNIFORM_SUBDIV ) THEN
              VERNEW(:,1,1) = ( VEROLD(:,1) + VEROLD(:,4) )/2
              VERNEW(:,2,1) = ( VEROLD(:,1) + VEROLD(:,3) )/2
              VERNEW(:,3,1) = ( VEROLD(:,1) + VEROLD(:,2) )/2
              VERNEW(:,4,1) = VEROLD(:,1)
              VERNEW(:,:,2) = VERNEW(:,:,1)
              VERNEW(:,4,2) = ( VEROLD(:,2) + VEROLD(:,4) )/2
              VERNEW(:,:,3) = VERNEW(:,:,2)
              VERNEW(:,3,3) = ( VEROLD(:,3) + VEROLD(:,4) )/2
              VERNEW(:,:,4) = VERNEW(:,:,2)
              VERNEW(:,1,4) = ( VEROLD(:,2) + VEROLD(:,3) )/2
              VERNEW(:,:,5) = VERNEW(:,:,3)
              VERNEW(:,2,5) = VEROLD(:,4)
              VERNEW(:,:,6) = VERNEW(:,:,4)
              VERNEW(:,3,6) = VERNEW(:,3,3)
              VERNEW(:,:,7) = VERNEW(:,:,6)
              VERNEW(:,4,7) = VEROLD(:,3)
              VERNEW(:,:,8) = VERNEW(:,:,4)
              VERNEW(:,2,8) = VEROLD(:,2)
              OUTSUB = 8
              NUM = 0
!
          ELSE IF (.NOT. UNIFORM_SUBDIV ) THEN
              CALL DIVSMP( DIMENS, NUMFUN, MAXSUB, VEROLD, Integrand,  &
                           NUM, OUTSUB, VERNEW )
          ELSE
              IFAIL = 7
              RETURN
          END IF
!
      CASE (Hyperrectangle)
          SELECT CASE (DIMENS)
          CASE (1)
              VERNEW(:,:,1) = VEROLD
              VERNEW(:,:,2) = VEROLD
              VERNEW(:,2,1) = (VEROLD(1,1)+VEROLD(1,2))/2
              VERNEW(:,1,2) = VERNEW(:,2,1)
              OUTSUB = 2
              NUM = 0
          CASE (2)
            IF ((MAXSUB == 4) .OR. UNIFORM_SUBDIV) THEN
!
!          Divide a parallellogram in 4.
!
              VERNEW(:,1,1) = VEROLD(:,1)
              VERNEW(:,2,1) = ( VEROLD(:,1) + VEROLD(:,2) )/2
              VERNEW(:,3,1) = ( VEROLD(:,1) + VEROLD(:,3) )/2
              VERNEW(:,1,2) = ( VEROLD(:,2) + VEROLD(:,3) )/2
              VERNEW(:,2,2) = VEROLD(:,2) + ( VEROLD(:,3) - VEROLD(:,1) )/2
              VERNEW(:,3,2) = VEROLD(:,3) + ( VEROLD(:,2) - VEROLD(:,1) )/2
              VERNEW(:,1,3) = VERNEW(:,2,1)
              VERNEW(:,2,3) = VEROLD(:,2)
              VERNEW(:,3,3) = VERNEW(:,1,2)
              VERNEW(:,1,4) = VERNEW(:,3,1)
              VERNEW(:,2,4) = VERNEW(:,1,2)
              VERNEW(:,3,4) = VEROLD(:,3)
              OUTSUB = 4
              NUM = 0
          ELSE IF (MAXSUB == 2) THEN
!
!          Divide a parallellogram in 2.
!
              IF ( INFOLD(4) == 2 ) THEN
!                Cut orthogonal to the line through vertices 1-3
                 VERNEW(:,1,1) = VEROLD(:,1)
                 VERNEW(:,2,1) = VEROLD(:,2)
                 VERNEW(:,3,1) = ( VEROLD(:,1) + VEROLD(:,3) )/2
                 VERNEW(:,1,2) = VERNEW(:,3,1)
                 VERNEW(:,2,2) = VEROLD(:,2)+( VEROLD(:,3) - VEROLD(:,1) )/2
                 VERNEW(:,3,2) = VEROLD(:,3)
              ELSE
!                Cut orthogonal to the line through vertices 1-2
                 VERNEW(:,1,1) = VEROLD(:,1)
                 VERNEW(:,2,1) = ( VEROLD(:,1)+VEROLD(:,2) )/2
                 VERNEW(:,3,1) = VEROLD(:,3)
                 VERNEW(:,1,2) = VERNEW(:,2,1)
                 VERNEW(:,3,2) = VEROLD(:,2)
                 VERNEW(:,2,2) = VEROLD(:,3)+( VEROLD(:,2) - VEROLD(:,1) )/2
              END IF
              OUTSUB = 2
              NUM = 0
          ELSE
              IFAIL = 7
              RETURN
          END IF
!
          CASE (3)
          IF ((MAXSUB == 8) .OR. UNIFORM_SUBDIV) THEN
!
!          Divide a 3D-Cube in 8.
!
                 HALF_HEIGHT = (verold(1:3,4)-verold(1:3,1))/2                                      
                 vernew(1:3,1,1) = verold(1:3,1)
                 vernew(1:3,2,1) = (verold(1:3,1)+verold(1:3,2))/2
                 vernew(1:3,3,1) = (verold(1:3,1)+verold(1:3,3))/2
                 vernew(1:3,4,1) = (verold(1:3,1)+verold(1:3,4))/2
                 vernew(1:3,1,2) = vernew(1:3,2,1)
                 vernew(1:3,2,2) = verold(1:3,2)
                 vernew(1:3,3,2) = (verold(1:3,2)+verold(1:3,3))/2
                 vernew(1:3,4,2) = (verold(1:3,2)+verold(1:3,4))/2
                 vernew(1:3,1,3) = vernew(1:3,3,1)
                 vernew(1:3,2,3) = vernew(1:3,3,2)
                 vernew(1:3,3,3) = verold(1:3,3)
                 vernew(1:3,4,3) = (verold(1:3,3)+verold(1:3,4))/2
                 vernew(1:3,1,4) = vernew(1:3,3,2)
                 vernew(1:3,2,4) = verold(1:3,2) +                      &
                                   (verold(1:3,3)-verold(1:3,1))/2
                 vernew(1:3,3,4) = verold(1:3,3) +                      &
                                   (verold(1:3,2)-verold(1:3,1))/2
                 vernew(1:3,4,4) = vernew(1:3,3,2) + HALF_HEIGHT
                 vernew(1:3,1,5) = vernew(1:3,4,1)
                 vernew(1:3,2,5) = vernew(1:3,4,2)
                 vernew(1:3,3,5) = vernew(1:3,4,3)
                 vernew(1:3,4,5) = verold(1:3,4)
                 vernew(1:3,1,6) = vernew(1:3,4,2)
                 vernew(1:3,2,6) = verold(1:3,2) +  HALF_HEIGHT
                 vernew(1:3,3,6) = vernew(1:3,4,4)
                 vernew(1:3,4,6) = vernew(1:3,4,2) + HALF_HEIGHT
                 vernew(1:3,1,7) = vernew(1:3,4,3)
                 vernew(1:3,2,7) = vernew(1:3,4,4)
                 vernew(1:3,3,7) = verold(1:3,3) + HALF_HEIGHT
                 vernew(1:3,4,7) = vernew(1:3,4,3) + HALF_HEIGHT
                 vernew(1:3,1,8) = vernew(1:3,4,4)
                 vernew(1:3,2,8) = vernew(1:3,2,4) + HALF_HEIGHT
                 vernew(1:3,3,8) = vernew(1:3,3,4) + HALF_HEIGHT
                 vernew(1:3,4,8) = vernew(1:3,4,4) + HALF_HEIGHT
              OUTSUB = 8
              NUM = 0
          
          ELSE IF (maxsub == 4) THEN
!
!          Divide a 3D-Cube in 4.
!
              IF (infold(4)/10 == -2) THEN 
!     Cut orthogonal to the line through vertices 1-2 and vertices 1-3
                 HEIGHT = (verold(1:3,4) - verold(1:3,1))
                 vernew(1:3,1,1) = verold(1:3,1)
                 vernew(1:3,2,1) = (verold(1:3,1)+verold(1:3,2))/2
                 vernew(1:3,3,1) = (verold(1:3,1)+verold(1:3,3))/2
                 vernew(1:3,4,1) = verold(1:3,4)
                 vernew(1:3,1,2) = vernew(1:3,2,1)
                 vernew(1:3,2,2) = verold(1:3,2)
                 vernew(1:3,3,2) = (verold(1:3,2)+verold(1:3,3))/2
                 vernew(1:3,4,2) = vernew(1:3,2,1) + HEIGHT
                 vernew(1:3,1,3) = vernew(1:3,3,1)
                 vernew(1:3,2,3) = vernew(1:3,3,2)
                 vernew(1:3,3,3) = verold(1:3,3)
                 vernew(1:3,4,3) = vernew(1:3,3,1) + HEIGHT
                 vernew(1:3,1,4) = vernew(1:3,3,2)
                 vernew(1:3,2,4) = vernew(1:3,3,1) +                      &
                                   (verold(1:3,2) - verold(1:3,1))
                 vernew(1:3,3,4) = vernew(1:3,2,1) +                      &
                                   (verold(1:3,3) - verold(1:3,1))
                 vernew(1:3,4,4) = vernew(1:3,3,2) + HEIGHT
              ELSE IF (infold(4)/10 == -3) THEN
!     Cut orthogonal to the line through vertices 1-3 and vertices 1-4
                 HEIGHT = (verold(1:3,2) - verold(1:3,1))
                 HALF_HEIGHT = (verold(1:3,4) - verold(1:3,1))/2
                 vernew(1:3,1,1) = verold(1:3,1)
                 vernew(1:3,2,1) = verold(1:3,2)
                 vernew(1:3,3,1) = (verold(1:3,1)+ verold(1:3,3))/2
                 vernew(1:3,4,1) = verold(1:3,1) + HALF_HEIGHT
                 vernew(1:3,1,2) = vernew(1:3,3,1)
                 vernew(1:3,2,2) = vernew(1:3,3,1) + HEIGHT
                 vernew(1:3,3,2) = verold(1:3,3)
                 vernew(1:3,4,2) = vernew(1:3,3,1) + HALF_HEIGHT
                 vernew(1:3,1,3) = vernew(1:3,4,1)
                 vernew(1:3,2,3) = vernew(1:3,4,1) + HEIGHT
                 vernew(1:3,3,3) = vernew(1:3,4,2)
                 vernew(1:3,4,3) = verold(1:3,4)
                 vernew(1:3,1,4) = vernew(1:3,4,2)
                 vernew(1:3,2,4) = vernew(1:3,4,2) + HEIGHT
                 vernew(1:3,3,4) = verold(1:3,3) + HALF_HEIGHT
                 vernew(1:3,4,4) = vernew(1:3,4,2) + HALF_HEIGHT
              ELSE
!     Cut orthogonal to the line through vertices 1-2 and vertices 1-4
                 HALF_HEIGHT = (verold(1:3,4) - verold(1:3,1))/2
                 vernew(1:3,1,1) = verold(1:3,1)
                 vernew(1:3,2,1) = (verold(1:3,2) + verold(1:3,1))/2
                 vernew(1:3,3,1) = verold(1:3,3)
                 vernew(1:3,4,1) = verold(1:3,1) + HALF_HEIGHT
                 vernew(1:3,1,2) = vernew(1:3,2,1)
                 vernew(1:3,2,2) = verold(1:3,2)
                 vernew(1:3,3,2) = verold(1:3,3) +                       &
                                   (verold(1:3,2) - verold(1:3,1))/2
                 vernew(1:3,4,2) = vernew(1:3,2,1) + HALF_HEIGHT
                 vernew(1:3,1,3) = vernew(1:3,4,1)
                 vernew(1:3,2,3) = (verold(1:3,2) + verold(1:3,4))/2
                 vernew(1:3,3,3) = verold(1:3,3) + HALF_HEIGHT
                 vernew(1:3,4,3) = verold(1:3,4)
                 vernew(1:3,1,4) = vernew(1:3,4,2)
                 vernew(1:3,2,4) = verold(1:3,2) + HALF_HEIGHT
                 vernew(1:3,3,4) = vernew(1:3,3,2) + HALF_HEIGHT
                 vernew(1:3,4,4) = vernew(1:3,4,2) + HALF_HEIGHT

              END IF
              OUTSUB = 4
              NUM = 0

          ELSE IF (maxsub == 2) THEN
!
!          Divide a 3D-Cube in 2.
!
              IF (modulo(abs(infold(4)),10) == 2) THEN 
!     Cut orthogonal to the line through vertices 1-3
                 vernew(1:3,1,1) = verold(1:3,1)
                 vernew(1:3,2,1) = verold(1:3,2)
                 vernew(1:3,3,1) = (verold(1:3,1)+verold(1:3,3))/2
                 vernew(1:3,4,1) = verold(1:3,4)
                 vernew(1:3,1,2) = vernew(1:3,3,1)
                 vernew(1:3,2,2) = vernew(1:3,3,1) +                     &
                                   (verold(1:3,2) - verold(1:3,1))
                 vernew(1:3,3,2) = verold(1:3,3)
                 vernew(1:3,4,2) = vernew(1:3,3,1) +                     &
                                   (verold(1:3,4) - verold(1:3,1))
              ELSE IF (modulo(abs(infold(4)),10) == 3) THEN
!     Cut orthogonal to the line through vertices 1-4
                 vernew(1:3,1,1) = verold(1:3,1)
                 vernew(1:3,2,1) = verold(1:3,2)
                 vernew(1:3,3,1) = verold(1:3,3)
                 vernew(1:3,4,1) = (verold(1:3,1)+verold(1:3,4))/2
                 vernew(1:3,1,2) = vernew(1:3,4,1)
                 vernew(1:3,2,2) = vernew(1:3,4,1) +                     &
                                   (verold(1:3,2) - verold(1:3,1))
                 vernew(1:3,3,2) = vernew(1:3,4,1) +                     &
                                   (verold(1:3,3) - verold(1:3,1))
                 vernew(1:3,4,2) = verold(1:3,4)
              ELSE
!     Cut orthogonal to the line through vertices 1-2
                 vernew(1:3,1,1) = verold(1:3,1)
                 vernew(1:3,2,1) = (verold(1:3,1)+verold(1:3,2))/2
                 vernew(1:3,3,1) = verold(1:3,3)
                 vernew(1:3,4,1) = verold(1:3,4)
                 vernew(1:3,1,2) = vernew(1:3,2,1)
                 vernew(1:3,2,2) = verold(1:3,2)
                 vernew(1:3,3,2) = vernew(1:3,2,1) +                     &
                                   (verold(1:3,3) - verold(1:3,1))
                 vernew(1:3,4,2) = vernew(1:3,2,1) +                     &
                                   (verold(1:3,4) - verold(1:3,1))
              END IF
              OUTSUB = 2
              NUM = 0
          ELSE
              IFAIL = 7
              RETURN
          END IF
          CASE DEFAULT
!
!          IF DIMENS > 3, then divide in 2 according to infold(4).
!
              J = INFOLD(4) 
              IF ( J == 0 ) THEN
                                   J = 1
              END IF
              DO I = 1,NRVERT
                 IF ( J == I-1 ) THEN
                    VERNEW(:,I,1) = (VEROLD(:,1) + VEROLD(:,I))/2
                 ELSE 
                    VERNEW(:,I,1) = VEROLD(:,I)
                 END IF
              END DO
              HALF_HEIGHT(1) = (VEROLD(J,J+1) - VEROLD(J,1))/2
              DO I = 1,NRVERT
                 IF ( J == I-1 ) THEN
                    VERNEW(:,I,2) = VEROLD(:,I)
                 ELSE 
                    VERNEW(:,I,2) = VERNEW(:,I,1)
                    VERNEW(J,I,2) = VERNEW(J,I,2) + HALF_HEIGHT(1)
                 END IF
              END DO
              OUTSUB = 2
              NUM = 0
          END SELECT
      CASE DEFAULT
          IFAIL = 6
          RETURN
      END SELECT
      VOLUME = RINFOL(1)/OUTSUB
!
!  Update/copy information record
!
      INFNEW(1,1:OUTSUB) = GEOMETRY
      INFNEW(2,1:OUTSUB) = INFOLD(2) + 1
      INFNEW(3,1:OUTSUB) = INFOLD(3)
      INFNEW(4:5,1:OUTSUB) = 0
      RINFNE(1,1:OUTSUB) = VOLUME
      DO I = 1,OUTSUB
          INFNEW(6:,I) = INFOLD(6:)
          RINFNE(2:,I) = RINFOL(2:)
      END DO
      IFAIL = 0
      RETURN
      END SUBROUTINE DIVIDE


      SUBROUTINE DIVSMP( DIMENS, NF, MAXSUB, VEROLD, Integrand,       &
                         FUNCLS, OUTSUB, VERNEW )
!***BEGIN PROLOGUE DIVSMP
!***PURPOSE  To compute new subregions
!***AUTHOR
!
!          Alan Genz
!          Department of Mathematics
!          Washington State University
!          Pullman, WA 99164-3113, USA
!          AlanGenz@wsu.edu
!
!          Ronald Cools, Dept. of Computer Science,
!          Katholieke Universiteit Leuven, Celestijnenlaan 200A,
!          B-3001 Heverlee, Belgium
!          Email: Ronald.Cools@cs.kuleuven.ac.be
!
!
!***LAST MODIFICATION 97-03
!***DESCRIPTION DIVSMP computes fourth differences along each edge
!            direction. It uses these differences to determine a
!            subdivision of the orginal subregion into two new subregions.
!
!   ON ENTRY
!
!   DIMENS  Integer number of variables.
!   NF      Integer number of components for the vector integrand.
!   MAXSUB  Integer.
!           The given region is divided into at most MAXSUB subregions.
!   VEROLD  Real array of dimension (N,0:N), orginal subregion vertices.
! Integrand Real vector function of length NF for computing components of 
!            the integrand at Z.
!            It must have parameters ( NF, Z ). See interface below.
!            Input parameters:
!              Z      Real array of length DIMENS, the evaluation point.
!              NF     Integer number of components of Integrand.
!
!   ON RETURN
!
!   OUTSUB Integer number of subregions the given region was divided into.
!   FUNCLS Integer number of Integrand calls used by DIVSMP.
!   VERNEW Real array of dimension (N,0:N,MAXSUB).
!          The vertices of the MAXSUB new subegions.
!
!***ROUTINES CALLED: Integrand
!***END PROLOGUE DIVSMP
      INTERFACE
         FUNCTION Integrand( NF, Z ) RESULT(Value)
             USE Precision_Model
         INTEGER, INTENT(IN) :: NF
         REAL(kind=STND), DIMENSION(:), INTENT(IN) :: Z
         REAL(kind=STND), DIMENSION(NF) :: Value
         END FUNCTION Integrand
      END INTERFACE
      INTEGER, INTENT(IN) :: DIMENS, NF, MAXSUB
      INTEGER, INTENT(OUT) :: OUTSUB, FUNCLS
      REAL(kind=STND), INTENT(IN), DIMENSION(1:,0:) :: VEROLD
      REAL(kind=STND), INTENT(OUT), DIMENSION(1:,0:,1:) :: VERNEW
!
!   Local Arrays
!   X      Real work array of length DIMENS.
!   H      Real work array of length DIMENS.
!   CENTER Real work array of length DIMENS.
!   WORK   Real work array of dimension (NF,5).
!   FRTHDF Real work array of dimension (0:DIMENS,0:DIMENS).
!   EWIDTH Real work array of dimension (0:DIMENS,0:DIMENS).
!
      REAL(kind=STND), DIMENSION(5), PARAMETER :: DIFCON = (/ 1, -4, 6, -4, 1 /)
      REAL(kind=STND), DIMENSION(DIMENS) :: X, H, CENTER
      REAL(kind=STND), DIMENSION(NF,5) :: WORK
      REAL(kind=STND), DIMENSION(0:DIMENS,0:DIMENS) :: FRTHDF, EWIDTH
      REAL(kind=STND) :: DIFMID, DIFMAX, DIFMIN
      INTEGER :: J, K, L, IMX, JMX, KMX, LMX, LTMP
      INTEGER, DIMENSION(2) :: INDX
!
!***FIRST PROCESSING STATEMENT DIVSMP
!
!
!     Compute the differences.
!
      CENTER = MATMUL( VEROLD, SPREAD( 1, 1, DIMENS+1 ) )/(DIMENS+1)
      WORK(:,3) = Integrand( NF, CENTER )
      FRTHDF = 0 
      EWIDTH = 0
      DO L = 0, DIMENS-1
         DO K = L+1, DIMENS
            H = VEROLD(:,K) - VEROLD(:,L)
            EWIDTH(L,K) = SUM( ABS( H ) ) 
            
            H = 2*H/( 5*( DIMENS + 1 ) )
            X = CENTER - 3*H
            DO J = 1, 5
               X = X + H
               IF ( J  /=  3 ) THEN
                                     WORK(:,J) = Integrand( NF, X )
               END IF
            END DO
            DIFMID = SUM( ABS( MATMUL( WORK(:,:), DIFCON ) ) )
!
!           Ignore differences below roundoff
!     
            IF ( SUM( ABS(WORK(:,3)) ) + DIFMID/8 > SUM( ABS(WORK(:,3)) ) ) THEN 
                 FRTHDF(L,K) = DIFMID
            END IF
         END DO
      END DO      
      FRTHDF = FRTHDF*EWIDTH 
      IF ( MAXVAL( FRTHDF ) == 0 ) THEN
                                         FRTHDF = EWIDTH
      END IF
      INDX(1:2) = MAXLOC( FRTHDF ) - 1
      LMX = INDX(1)
      KMX = INDX(2)
      VERNEW(:,:,1) = VEROLD
      VERNEW(:,:,2) = VEROLD
!
      IF ((DIMENS == 2) .AND. (MAXSUB >=4)) THEN
         DIFMAX = FRTHDF(LMX,KMX)
         DIFMIN = MIN(frthdf(0,1),frthdf(0,2),frthdf(1,2))
         IF ( (DIFMAX > 0.001*ewidth(lmx,kmx)) .AND. ( DIFMIN <= 0.45_stnd*DIFMAX )) THEN  ! 2-division
            VERNEW(:,LMX,2) = ( VEROLD(:,KMX) + VEROLD(:,LMX) )/2
            VERNEW(:,KMX,1) =   VERNEW(:,LMX,2)
            OUTSUB = 2
         ELSE ! 4-division
            ! VERNEW(:,:,1) = VEROLD ; VERNEW(:,:,2) = VEROLD
            VERNEW(:,:,3) = VEROLD
            VERNEW(:,:,4) = VEROLD
            VERNEW(:,1,1) = ( VEROLD(:,0) + VEROLD(:,1) )/2
            VERNEW(:,2,1) = ( VEROLD(:,0) + VEROLD(:,2) )/2
            VERNEW(:,0,2) = VERNEW(:,1,1)
            VERNEW(:,2,2) = ( VEROLD(:,1) + VEROLD(:,2) )/2
            VERNEW(:,0,3) = VERNEW(:,2,1)
            VERNEW(:,1,3) = VERNEW(:,2,2)
            VERNEW(:,0,4) = VERNEW(:,2,2)
            VERNEW(:,1,4) = VERNEW(:,1,1)
            VERNEW(:,2,4) = VERNEW(:,2,1)      
            OUTSUB = 4
         END IF
      ELSE IF ( MAXSUB == 2 ) THEN
!
!        Compute two new subregions.
!
         VERNEW(:,LMX,2) = ( VEROLD(:,KMX) + VEROLD(:,LMX) )/2
         VERNEW(:,KMX,1) =   VERNEW(:,LMX,2)
         OUTSUB = 2
      ELSE 
         DIFMAX = FRTHDF(LMX,KMX)
!
         FRTHDF(LMX,KMX) = 0
         INDX(1:2) = MAXLOC( FRTHDF ) - 1
         JMX = INDX(1)
         IMX = INDX(2)
!
         FRTHDF(JMX,IMX) = 0
         INDX(1:2) = MAXLOC( FRTHDF ) - 1
         JMX = INDX(1)
         IMX = INDX(2)
!
         DIFMID = FRTHDF(JMX,IMX)
         IF ( DIFMAX > 2*DIFMID .OR. MAXSUB == 3 ) THEN   
! Tobedone            ^^ tune this parameter                RC
!
!           Compute three new subregions.
!
            VERNEW(:,:,3) = VEROLD
            WHERE ( FRTHDF == 0 )
                  FRTHDF = TRANSPOSE( FRTHDF )
            END WHERE
            INDX(1:1) = MAXLOC( FRTHDF(LMX,:) + FRTHDF(KMX,:) ) - 1
            JMX = INDX(1)
            IF ( FRTHDF(LMX,JMX) > FRTHDF(KMX,JMX) ) THEN
               LTMP = KMX
               KMX = LMX
               LMX = LTMP
            END IF
            DIFMID = FRTHDF(KMX,JMX)
            VERNEW(:,KMX,1) = ( 2*VEROLD(:,LMX) + VEROLD(:,KMX) )/3
            VERNEW(:,LMX,2) =     VERNEW(:,KMX,1)
            VERNEW(:,KMX,2) =     VEROLD(:,JMX)
            IF ( DIFMID > DIFMAX/8 ) THEN
! Tobedone                       ^^ tune this parameter              RC
               VERNEW(:,JMX,2) = ( VEROLD(:,KMX) + VEROLD(:,JMX) )/2
               VERNEW(:,JMX,3) =   VERNEW(:,LMX,2)
            ELSE
               VERNEW(:,JMX,2) = ( VEROLD(:,LMX) + 2*VEROLD(:,KMX) )/3
               VERNEW(:,JMX,3) =   VEROLD(:,JMX)
            END IF
            VERNEW(:,LMX,3) = VERNEW(:,JMX,2)
            OUTSUB = 3
         ELSE
!
!           Compute four new subregions.
!
            VERNEW(:,LMX,2) = ( VEROLD(:,KMX) + VEROLD(:,LMX) )/2
            VERNEW(:,KMX,1) =   VERNEW(:,LMX,2)
            VERNEW(:,:,3) =     VERNEW(:,:,1)
            VERNEW(:,JMX,3) = ( VERNEW(:,IMX,1) + VERNEW(:,JMX,1) )/2
            VERNEW(:,IMX,1) =   VERNEW(:,JMX,3)
            VERNEW(:,:,4) =     VERNEW(:,:,2)
            VERNEW(:,JMX,4) = ( VERNEW(:,IMX,2) + VERNEW(:,JMX,2) )/2
            VERNEW(:,IMX,2) =   VERNEW(:,JMX,4)
            OUTSUB = 4
         END IF
      END IF
      FUNCLS = 1 + 2*DIMENS*(DIMENS+1)
      RETURN
      END SUBROUTINE DIVSMP
END MODULE Subdivisions
