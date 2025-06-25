      PROGRAM TEST_FILTER1D
      IMPLICIT NONE
      CHARACTER(LEN=*) SRC
      PARAMETER(SRC='coefs.txt')
      INTEGER I,N
      DOUBLE PRECISION PI,DX
      PARAMETER(N=600,PI=4.0*ATAN(1.0))
      DOUBLE PRECISION X(N),Y(N),FY(N),RY(N)

      DX=2.0*PI/(N-1)
      DO I=1,N
       X(I)=(I-1)*DX
       Y(I)=SIN(X(I))
      ENDDO

      CALL PYDOF_INIT(NX,1,1,'P','P','P',SRC)
      CALL PYDOF_FILTER(Y,'F',FY)
      CALL PYDOF_FILTER(Y,'R',RY)

      STOP
      END
