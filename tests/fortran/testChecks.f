      PROGRAM TEST_CHECKS
      IMPLICIT NONE 
      CHARACTER(LEN=*) SRC
      PARAMETER(SRC='./coeffs.txt')

      CALL PYDOF_INIT(0,1,1,'C','C','C',SRC)
      CALL PYDOF_INIT(1,0,1,'C','C','C',SRC)
      CALL PYDOF_INIT(1,1,0,'C','C','C',SRC)
      CALL PYDOF_INIT(1,1,1,'K','C','C',SRC)
      CALL PYDOF_INIT(1,1,1,'C','K','C',SRC)
      CALL PYDOF_INIT(1,1,1,'C','C','K',SRC)

      STOP
      END
