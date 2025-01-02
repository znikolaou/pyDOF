      PROGRAM MAIN
      IMPLICIT NONE
      INCLUDE "PYDOF_COMMON.h"

      INTEGER I, N=500
      DOUBLE PRECISION PI=4.0*ATAN(1.0)
      DOUBLE PRECISION DX, X(N), Y(N), YBAR(N)

      !INIT FILTER
      CALL PYDOF_INIT()

      DX=2.0*PI/(N-1)
      DO I=1, N
       X(I)=I*DX
       Y(I)=SIN(X(I))
      ENDDO

      YBAR=FILTER(Y)

      OPEN(UNI=1,FILE='output.txt',STATUS='OLD',FORM='FORMATTED')
      WRITE(1,*) 'X', 'Y', 'YBAR'
      DO I=1,N 
       WRITE(1,*) X(I),Y(I),YBAR(I)
      ENDDO
      CLOSE(1)

      END
