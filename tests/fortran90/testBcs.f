      PROGRAM TEST_PYDOF_INIT
      IMPLICIT NONE
      INCLUDE "pydof_common.h"
      INTEGER I,J,K,NX
      PARAMETER(NX=5)
      CHARACTER(LEN=2) BCX(2),BCY(2),BCZ(2)
      
      BCX(1)='DI'
      BCX(2)='CP'
      BCY(1)='DP'
      BCY(2)='DI'
      BCZ(1)='CP'
      BCZ(2)='DP'
      
      CALL PYDOF_INIT(NX,1,1,BCX,BCY,BCZ)

      DO K=1,NF
        WRITE(*,*) 'FILTER=',K
       DO J=1,NX
        WRITE(*,*) '  POINT=',J
        DO I=-NGMX,NGMX
         WRITE(*,*) '   ',I,IL(I,J,K)
        ENDDO
       ENDDO
      ENDDO

      STOP
      END
