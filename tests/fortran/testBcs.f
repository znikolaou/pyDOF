      PROGRAM TEST_PYDOF_INIT
      IMPLICIT NONE
      INCLUDE "pydof_common.h"
      INTEGER I,J,K,NX
      PARAMETER(NX=20)
      
      WRITE(*,*) 'CONTINUOUS'
      CALL PYDOF_INIT(NX,1,1,'C','C','C')

      DO K=1,NF
        WRITE(*,*) 'FILTER=',K
       DO J=1,NX
        WRITE(*,*) '  POINT=',J
        DO I=-NGMX,NGMX
         WRITE(*,*) '   ',I,IL(I,J,K)
        ENDDO
       ENDDO
      ENDDO

      WRITE(*,*) 'PERIODIC'
      CALL PYDOF_INIT(NX,1,1,'P','P','P')

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
