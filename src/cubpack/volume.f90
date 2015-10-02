! This file is F-compatible, except for upper/lower case conventions.
!--------------------------------------------------------------------
Module Volume_Computation

USE Precision_Model, ONLY: stnd

Implicit NONE

PRIVATE
PUBLIC :: VOLUME

CONTAINS
      FUNCTION VOLUME(DIMENS,GEOMETRY,VERTIC) RESULT(Value)
!***BEGIN PROLOGUE VOLUME
!***PURPOSE  To compute the volume of a polytope
!***AUTHORS
!        Ronald Cools                     Alan Genz
!        Dept. of Computer Science        Computer Science Department
!        Katholieke Universiteit Leuven   Washington State University
!        Celestijnenlaan 200A             Pullman, WA 99164-2752
!        B-3001 Heverlee, Belgium         USA
!
!***REVISION DATE  950419   (YYMMDD) (Fortran90 transformation)
!***REVISION DATE  980408   (YYMMDD) (F         transformation)
!***DESCRIPTION VOLUME
!
!   GEOMETRY = 1 : the VERTICes specify a simplex
!   GEOMETRY = 2 : the VERTICes specify a cube
!   GEOMETRY = 3 : the VERTICes specify an octahedron
!
!  WARNING: If a region of an unsupported shape is presented,
!           then this function returns 0.
!           This is the only indication of a failure !
!
!   Global variables.
!
      INTEGER, INTENT(IN) :: DIMENS,GEOMETRY
      REAL(kind=stnd), DIMENSION(:,:), INTENT(IN) :: VERTIC
      REAL(kind=stnd) :: Value
!
!   Local variables.
!
      INTEGER         :: I,J,K,PIVPOS,FACDIM
      REAL(kind=stnd) :: MULT,VOL
      REAL(kind=stnd), DIMENSION(DIMENS) :: TMP
      REAL(kind=stnd), DIMENSION(DIMENS,DIMENS) :: WORK

      IF ((GEOMETRY == 1) .OR. (GEOMETRY == 2) .OR. (GEOMETRY == 3)) THEN
          SELECT CASE (DIMENS)
          CASE (1) !    Compute length of an interval
              Value = ABS(VERTIC(1,2)-VERTIC(1,1))

          CASE (2) !    Compute area of a rectangle.
              Value = ABS((VERTIC(1,2)-VERTIC(1,1))*                 &
                       (VERTIC(2,3)-VERTIC(2,1))-                    &
                       (VERTIC(2,2)-VERTIC(2,1))*                    &
                       (VERTIC(1,3)-VERTIC(1,1)))

          CASE (3) !    Compute the volume of a cube.
              Value = ABS((VERTIC(1,2)-VERTIC(1,1))*                 &
                       ((VERTIC(2,3)-VERTIC(2,1))* (VERTIC(3,        &
                       4)-VERTIC(3,1))- (VERTIC(2,4)-VERTIC(2,       &
                       1))* (VERTIC(3,3)-VERTIC(3,1)))-              &
                       (VERTIC(2,2)-VERTIC(2,1))*                    &
                       ((VERTIC(1,3)-VERTIC(1,1))* (VERTIC(3,        &
                       4)-VERTIC(3,1))- (VERTIC(1,4)-VERTIC(1,       &
                       1))* (VERTIC(3,3)-VERTIC(3,1)))+              &
                       (VERTIC(3,2)-VERTIC(3,1))*                    &
                       ((VERTIC(1,3)-VERTIC(1,1))* (VERTIC(2,        &
                       4)-VERTIC(2,1))- (VERTIC(1,4)-VERTIC(1,       &
                       1))* (VERTIC(2,3)-VERTIC(2,1))))

          CASE DEFAULT !  Compute the volume of a DIMENS-dimensional cube
              DO J = 1,DIMENS
                 WORK(1:DIMENS,J) = VERTIC(1:DIMENS,J+1) - VERTIC(1:DIMENS,1)
              END DO
              VOL = 1
              DO K = 1,DIMENS
                  PIVPOS = K
                  DO J = K + 1,DIMENS
                     IF (ABS(WORK(K,J)) > ABS(WORK(K,PIVPOS))) THEN
                        PIVPOS = J
                     END IF
                  END DO
                  TMP(K:DIMENS) = WORK(K:DIMENS,K)
                  WORK(K:DIMENS,K) = WORK(K:DIMENS,PIVPOS)
                  WORK(K:DIMENS,PIVPOS) = TMP(K:DIMENS)
                  VOL = VOL*WORK(K,K)
                  DO J = K + 1,DIMENS
                      MULT = WORK(K,J)/WORK(K,K)
                      WORK(K+1:DIMENS,J) = WORK(K+1:DIMENS,J) - MULT*WORK(K+1:DIMENS,K)
                  END DO
              END DO
              Value = ABS(VOL)
          END SELECT
          IF ((GEOMETRY == 1) .OR. (GEOMETRY == 3)) THEN
!
!            The volume of an DIMENS-dimensional simplex is the
!            DIMENS! part of the volume of the cube.
!
              FACDIM = DIMENS
              DO I = 2,DIMENS - 1
                  FACDIM = FACDIM*I
              END DO
              Value = Value/FACDIM
          END IF
          IF (GEOMETRY == 3) THEN
!
!            The volume of an DIMENS-dimensional octahedron is
!            2**(DIMENS-1) times the volume of the simplex.
!
              Value = Value*2**(DIMENS-1)
          END IF
      ELSE
          Value = 0
      END IF
      RETURN
      END FUNCTION VOLUME

END MODULE Volume_Computation
