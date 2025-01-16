      PROGRAM TEST_PYDOF_INIT
      IMPLICIT NONE
      INCLUDE "pydof_common.h"
      INTEGER I,N
      DOUBLE PRECISION PI,LX
      PARAMETER(N=400,PI=4.0*ATAN(1.0),LX=1.0) 
      DOUBLE PRECISION Y(N),X(N),FY(N),FYA(N),DX,DELTA,LAMDA,K

      LAMDA=LX/16.0
      DX=LX/(N-1)
      DELTA=8.0*DX
      K=2.0*PI/LAMDA
      DO I=1,N
       X(I)=(I-1)*DX
       Y(I)=SIN(K*X(I))
       FYA(I)=Y(I)*EXP(-DELTA**2.0*K**2.0/24.0)
      ENDDO

      CALL PYDOF_INIT(N,1,1,'P','P','P')
      CALL PYDOF_FILTER(1,N,1,1,Y,FY)
      OPEN(UNIT=1,FILE='data.txt',STATUS='NEW',FORM='FORMATTED')
      WRITE(1,*) 'x,','y,','fy,','fya'
      DO I=1,N
       WRITE(1,*) X(I),',',Y(I),',',FY(I),',',FYA(I)
      ENDDO
      CLOSE(1)

      STOP
      END
