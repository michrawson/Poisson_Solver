C FILE: ARRAY.F
      SUBROUTINE SOLVER(N,POT)
C
C     INCREMENT THE FIRST ROW AND DECREMENT THE FIRST COLUMN OF A
C
      INTEGER N
      REAL*8 POT(N,N,N)
Cf2py intent(inplace) POT

      call solve( n, pot )

      END

