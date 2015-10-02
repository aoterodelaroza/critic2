!>
!> @file   tqlgrm.f
!>
!> @brief
!> Part of the ??? library. 
!>
        SUBROUTINE      TQLGRM        (N, D, E, Z, IERR)
        IMPLICIT        DOUBLE PRECISION (A-H, O-Z)
        DIMENSION       D(*), E(*), Z(N,N)
        PARAMETER (AMACH = 16.0E-13)
        PARAMETER (ZERO = 0.0D0, ONE = 1.0D0)
C
        IERR    = 0
        IF (N .EQ. 1) RETURN
C
        DO 30 I = 2,N
        E(I-1)  = E(I)
30      CONTINUE

        F       = ZERO
        B       = ZERO
        E(N)    = ZERO

        DO 31 L = 1,N
        J       = 0
        H       = AMACH*(DABS(D(L)) + DABS(E(L)))
        IF (B .LT. H) B = H

                DO 32 M = L,N
                IF (DABS(E(M)) .LE. B) GOTO 120
32              CONTINUE

120     IF (M .EQ. L) GOTO 220

130     IF (J .EQ. 30) THEN
                IERR    = L
                RETURN
        END IF

        J       = J + 1
        L1      = L + 1
        G       = D(L)
        P       = (D(L1) - G)/(2*E(L))
        IF (DABS(P*AMACH) .GT. ONE) THEN
                R       = P
        ELSE
                R       = DSQRT(P*P + 1)
        END IF
        D(L)    = E(L)/(P + DSIGN(R,P))
        H       = G - D(L)

                DO 33 I = L1,N
                D(I)    = D(I) - H
33              CONTINUE

        F       = F + H
        P       = D(M)
        C       = ONE
        S       = ZERO
        MML     = M - L

                DO 34 II = 1,MML
                I       = M - II
                G       = C*E(I)
                H       = C*P
                IF (DABS(P) .GE. DABS(E(I))) THEN
                        C       = E(I)/P
                        R       = DSQRT(C*C + 1)
                        E(I+1)  = S*P*R
                        S       = C/R
                        C       = ONE/R
                ELSE
                        C       = P/E(I)
                        R       = DSQRT(C*C + 1)
                        E(I+1)  = S*E(I)*R
                        S       = 1.D0/R
                        C       = C*S
                END IF
                P       = C*D(I) - S*G
                D(I+1)  = H + S*(C*G + S*D(I))

                        DO 35 K = 1,N
                        H       = Z(K,I+1)
                        Z(K,I+1)= S*Z(K,I) + C*H
                        Z(K,I)  = C*Z(K,I) - S*H
35                      CONTINUE

34              CONTINUE

        E(L)    = S*P
        D(L)    = C*P
        IF (DABS(E(L)) .GT. B) GOTO 130

220     D(L)    = D(L) + F
31      CONTINUE

        DO 300 II = 2,N
        I       = II - 1
        K       = I
        P       = D(I)

                DO 260 J = II,N
                IF (D(J) .GE. P) GOTO 260
                K       = J
                P       = D(J)
260             CONTINUE

        IF (K .EQ. I) GOTO 300
        D(K)    = D(I)
        D(I)    = P

                DO 37 J = 1,N
                P       = Z(J,I)
                Z(J,I)  = Z(J,K)
                Z(J,K)  = P
37              CONTINUE

300     CONTINUE
        RETURN
        END
