      PROGRAM MAIN
      IMPLICIT NONE

      INTEGER I,N
      DOUBLE PRECISION PI
      PARAMETER(N=500,PI=4.0*ATAN(1.0))
      DOUBLE PRECISION DX,X(N),Y(N),FY(N),RY(N)

      CALL PYDOF_INIT(N,1,1,'P','P','P',1)

      DX=2.0*PI/(N-1)
      DO I=1,N
       X(I)=I*DX
       Y(I)=SIN(X(I))
      ENDDO

      CALL FILTER(Y,FY)
      CALL FILTER_R(FY,RY)

      OPEN(UNI=1,FILE='results.txt',STATUS='OLD',FORM='FORMATTED')
      WRITE(1,*) 'X','Y','FY','RY'
      DO I=1,N 
       WRITE(1,*) X(I),Y(I),FY(I),RY(I)
      ENDDO
      CLOSE(1)

      END
