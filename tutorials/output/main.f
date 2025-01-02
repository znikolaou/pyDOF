      PROGRAM MAIN
      IMPLICIT NONE
      INTEGER I,J,NM(2),IDUM
      CHARACTER(LEN=1) CH
      DOUBLE PRECISION COEFFS(100)

      OPEN(UNIT=1,FILE='coeffs.txt',STATUS='OLD',FORM='FORMATTED')

      J=1
      DO WHILE(1.EQ.1)               
       READ(1,'(A)',END=10) CH
       IF(CH.NE.'!'.AND.CH.NE.'#') THEN
        
        DO I=1,NM(J)
         READ(1,*) IDUM,COEFFS(I)
        ENDDO
        J=2
       ENDIF
      ENDDO
10    CONTINUE

      DO I=1,NM(1)
       WRITE(*,*) COEFFS(I)
      ENDDO
      DO I=NM(1),NM(2)
       WRITE(*,*) COEFFS(I)
      ENDDO

     
      CLOSE(1)

      STOP
      END
