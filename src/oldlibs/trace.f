!>
!> @file   trace.f
!>
!> @brief
!> Part of the ??? library. 
!>
        SUBROUTINE      TRACE (H, E, W, N, IERR)
c
c    trace -- diagonalize the symmetric real matrix h
c             trace, tredig and tqlgrm are not commented.
c
c    input: 
c    h(n,n)                   symmetric real matrix to diagonalize
c    e                        -- used in output --
c    w                        work array, at least w(n)
c    ierr                     -- used in output --
c
c    output:
c    h(n,n)                   orthonormal eigenvectors
c    e                        eigenvalues, smaller first
c    ierr                     error code (not zero means error)
C
c    original info:
C TRACE CALLS TREDIG AND TLQGRM TO DIAGONALIZE A SYMMETRIC REAL MATRIX.
C THE MATRIX IS PASSED DOWN IN H AND IS REPLACED BY THE EIGENVECTORS.
C THE EIGENVALUES IN E ARE STORED SMALLEST FIRST.
C THE WORK STORE W SHOULD BE AT LEAST OF DIMENSION N.
C SKK ==================================================================

        IMPLICIT        DOUBLE PRECISION (A-H, O-Z)
        DIMENSION       H(N,N), E(*), W(*)

        CALL TREDIG     (N, E, W, H)
        CALL TQLGRM     (N, E, W, H, IERR)

        RETURN
        END
