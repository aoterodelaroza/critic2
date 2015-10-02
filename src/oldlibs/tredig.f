!>
!> @file   tredig.f
!>
!> @brief 
!> Part of ??? library.
!> 
        SUBROUTINE      TREDIG (N, D, E, Z)
        IMPLICIT        DOUBLE PRECISION (A-H, O-Z)
        DIMENSION       D(*), E(*), Z(N,N)
        PARAMETER       (ZERO = 0.0D0, ONE = 1.0D0)

        IF (N .EQ. 1) GOTO 320

        DO 30 II = 2,N
        I       = N + 2 - II
        L       = I - 1
        H       = ZERO
        SCALE   = ZERO

        IF (L .LT. 2) GOTO 130

                DO 31 K = 1,L
                SCALE   = SCALE + DABS(Z(I,K))
31              CONTINUE

        IF (SCALE .NE. ZERO) GOTO 140
130     E(I)    = Z(I,L)
        GOTO 290

140     RSCALE = ONE/SCALE
                DO 32 K = 1,L
                Z(I,K)  = Z(I,K)*RSCALE
                H       = H + Z(I,K)*Z(I,K)
32              CONTINUE
        F       = Z(I,L)
        G       = -DSIGN(DSQRT(H),F)
        E(I)    = SCALE*G
        H       = H - F*G
        Z(I,L)  = F - G
        F       = ZERO
        RH      = ONE/H
        RHSCALE = RH*RSCALE        

                DO 33 J = 1,L
                Z(J,I)  = Z(I,J)*RHSCALE
                G       = ZERO

                        DO 34 K = 1,J
                        G       = G + Z(J,K)*Z(I,K)
34                      CONTINUE

                JP1 = J + 1
                IF (L .LT. JP1) GOTO 220

                        DO 35 K = JP1,L
                        G       = G + Z(K,J)*Z(I,K)
35                      CONTINUE

220             E(J)    = G*RH
                F       = F + E(J)*Z(I,J)
33              CONTINUE

        HH = F/(H + H)

                DO 36 J = 1,L
                F       = Z(I,J)
                G       = E(J) - HH*F
                E(J)    = G
                        DO 37 K = 1,J
                        Z(J,K)  = Z(J,K) - F*E(K) - G*Z(I,K)
37                      CONTINUE
36              CONTINUE

                DO 38 K        = 1,L
                Z(I,K)  =  SCALE*Z(I,K)
38              CONTINUE

290     D(I)    = H
30      CONTINUE

320     D(1)    = ZERO
        E(1)    = ZERO

        DO 500 I = 1,N
        L       = I - 1
        IF (D(I) .EQ. ZERO) GOTO 380

                DO 40 J = 1,L
                G       = ZERO

                        DO 41 K = 1,L
                        G       = G + Z(I,K)*Z(K,J)
41                      CONTINUE

                        DO 42 K = 1,L
                        Z(K,J)  = Z(K,J) - G*Z(K,I)
42                      CONTINUE

40              CONTINUE

380     D(I)    = Z(I,I)
        Z(I,I)  = ONE
        IF(L .LT. 1) GOTO 500

                DO 43 J = 1,L
                Z(J,I)  = ZERO
                Z(I,J)  = ZERO
43              CONTINUE

500     CONTINUE
        RETURN
        END
