      PROGRAM TEST_PYDOF_HEADER
      IMPLICIT NONE
      INCLUDE "pydof_common.h"
      INTEGER I,J

      WRITE(*,'(AXI3)') 'NFILTERS=',NFILTERS
      WRITE(*,'(AXI3)') 'NCMAX=',NC_MAX
      DO J=1,NFILTERS
       WRITE(*,'(3(AXI3X))') 'FILTER:',J,'NCOEFS=',NC(J),
     & 'FORWARD(1)/INVERSE(0):',ISFORWARD(J)
       DO I=0,NCM1(J)
        WRITE(*,'(E32.16)') COEFS(I,J)
       ENDDO
      ENDDO

      STOP
      END
