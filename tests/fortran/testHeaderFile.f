      PROGRAM TEST_PYDOF_HEADER
      IMPLICIT NONE
      INCLUDE "pydof_common.h"
      INTEGER I,J

      WRITE(*,'(AXI3)') 'NFILTERS=',NF
      WRITE(*,'(AXI3)') 'NCMAX=',NCMX
      DO J=1,NF
       WRITE(*,'(3(AXI3X))') 'FILTER:',J,'NCOEFS=',NCOEF(J),
     & 'FORWARD(1)/INVERSE(0):',ISFORWARD(J)
       DO I=0,NCOEFM1(J)
        WRITE(*,'(E32.16)') COEFS(I,J)
       ENDDO
      ENDDO

      STOP
      END
