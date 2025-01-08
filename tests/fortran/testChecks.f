      PROGRAM TEST_CHECKS
      IMPLICIT NONE 

      CALL PYDOF_INIT(0,1,1,'C','C','C')
      CALL PYDOF_INIT(1,0,1,'C','C','C')
      CALL PYDOF_INIT(1,1,0,'C','C','C')
      CALL PYDOF_INIT(1,1,1,'K','C','C')
      CALL PYDOF_INIT(1,1,1,'C','K','C')
      CALL PYDOF_INIT(1,1,1,'C','C','K')

      STOP
      END
