      PROGRAM TEST_FILTER_1D
      IMPLICIT NONE
      INCLUDE "pydof_common.h"
      INTEGER I,N
      DOUBLE PRECISION PI,LX
      PARAMETER(N=200,PI=4.0*ATAN(1.0),LX=1.0) 
      DOUBLE PRECISION Y(N),X(N),FY(N),FYA(N),RY(N),DX,DELTA,LAMDA,K
      CHARACTER(LEN=2) BCX(2),BCY(2),BCZ(2)
      DOUBLE PRECISION A(3,4,1)

      BCX(1)='CP'
      BCX(2)='CP'
      BCY(1)='CP'
      BCY(2)='CP'
      BCZ(1)='CP'
      BCZ(2)='CP'

      LAMDA=LX/16.0
      DX=LX/(N-1)
      DELTA=8.0*DX
      K=2.0*PI/LAMDA
      DO I=1,N
       X(I)=(I-1)*DX
       Y(I)=SIN(K*X(I))
       FYA(I)=Y(I)*EXP(-DELTA**2.0*K**2.0/24.0)
      ENDDO

      CALL PYDOF_INIT(8,N,1,1,BCX,BCY,BCZ)
      CALL PYDOF_FILTER(1,N,1,1,Y,FY)
      CALL PYDOF_FILTER(2,N,1,1,FY,RY)
     
      OPEN(UNIT=1,FILE='data.txt',STATUS='NEW',FORM='FORMATTED')
      WRITE(1,*) 'x,','y,','fy,','fya,','ry'
      DO I=1,N
       WRITE(1,*) X(I),',',Y(I),',',FY(I),',',FYA(I),',',RY(I)
      ENDDO
      CLOSE(1)

      STOP
      END
