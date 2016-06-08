C FILE: ARRAY.F
      SUBROUTINE SOLVER(N,INPUT,RESULT)
C
C     INCREMENT THE FIRST ROW AND DECREMENT THE FIRST COLUMN OF A
C
      INTEGER N
      REAL*8 INPUT(N,N,N)
Cf2py intent(inplace) INPUT
      REAL*8 RESULT(N,N,N)
Cf2py intent(out) RESULT

      call solve( n, input, result )

      END

