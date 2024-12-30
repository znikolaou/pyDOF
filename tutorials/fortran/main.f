      PROGRAM MAIN
      IMPLICIT NONE

      INTEGER I, N=500
      DOUBLE PRECISION PI=4.0*ATAN(1.0)
      DOUBLE PRECISION DX, X(N), Y(N), YBAR(N)

      DX=2.0*PI/(N-1)
      DO I=1, N
       X(I)=I*DX
       Y(I)=SIN(X(I))
      ENDDO

      YBAR=FILTER(Y)


      END
