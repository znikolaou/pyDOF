      PROGRAM TEST_FILTER_2D_RADIAL
      IMPLICIT NONE
      INCLUDE "pydof_common.h"
      INTEGER I,J,NX,NY
      DOUBLE PRECISION PI,LX,LY
      PARAMETER(NX=400,NY=NX,PI=4.0*ATAN(1.0),LX=1.0,LY=LX) 
      DOUBLE PRECISION DX,DY,DR,X(NX),Y(NY),KX,KY,KR,LAMDAX,LAMDAY
      DOUBLE PRECISION R,A(NX,NY),FA(NX,NY),RA(NX,NY)
      CHARACTER(LEN=2) BCX(2),BCY(2),BCZ(2)
      
      BCX(1)='CP'
      BCX(2)='CP'
      BCY(1)='CP'
      BCY(2)='CP'
      BCZ(1)='CP'
      BCZ(2)='CP'
      
      DX=2.0*LX/(NX-1)
      DY=2.0*LY/(NY-1)
      DR=SQRT(DX**2.0+DY**2.0)
      DO I=1,NX
       X(I)=(I-1)*DX-LX
      ENDDO
      DO J=1,NY
       Y(J)=(J-1)*DY-LX
      ENDDO

      LAMDAX=LX/16.0
      LAMDAY=LY/16.0
      KX=2.0*PI/LAMDAX
      KY=2.0*PI/LAMDAY
      KR=SQRT(KX**2.0+KY**2.0)
      DO J=1,NY
       DO I=1,NX
        R=SQRT(X(I)**2.0+Y(J)**2.0)
        A(I,J)=SIN(KR*R)
        !Y(I)*EXP(-DELTA**2.0*K**2.0/24.0)
       ENDDO
      ENDDO

      CALL PYDOF_INIT(1,NX,NY,1,BCX,BCY,BCZ)
      CALL PYDOF_FILTER(1,NX,NY,1,A,FA)
      CALL PYDOF_FILTER(2,NX,NY,1,FA,RA)
     
      OPEN(UNIT=1,FILE='data.bin',STATUS='NEW',FORM='UNFORMATTED')
      WRITE(1) NX,NY,A,FA,RA
      CLOSE(1)

      STOP
      END
