      PROGRAM TEST_PYDOF_INIT
      IMPLICIT NONE
      INCLUDE "pydof_common.h"
      INTEGER I,J,K,NX
      PARAMETER(NX=1000000)
      CHARACTER(LEN=2) BCX(2),BCY(2),BCZ(2)
      
      BCX(1)='DI'
      BCX(2)='CP'
      BCY(1)='DP'
      BCY(2)='DI'
      BCZ(1)='CP'
      BCZ(2)='DP'
      
      CALL PYDOF_INIT(NX,1,1,BCX,BCY,BCZ)

      STOP
      END
