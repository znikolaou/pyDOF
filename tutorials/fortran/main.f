      PROGRAM MAIN
      IMPLICIT NONE
      INCLUDE "pydof_common.h"
      INTEGER I


      DO I=1,NC(1)
       WRITE(*,*) COEFS(I,1)
      ENDDO

      STOP
      END
