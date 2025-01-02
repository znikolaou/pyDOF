      SUBROUTINE PYDOF_READ_COEFS(SRC,FORI,COEFS)
      IMPLICIT NONE
      INTEGER NM
      CHARACTER(LEN=1000) SRC,LINE 
      CHARACTER(LEN=1) FORI,CH
      DOUBLE PRECISION COEFS(:),GETDOUBLE
      LOGICAL ISCOMMENT, TEST

      !
      LINE='!SFSF'
      WRITE(6,*) GETDOUBLE(4.0)
      !TEST=ISCOMMENT('!sfsd')

      WRITE(*,*) TEST
  
      STOP
      OPEN(UNIT=1,FILE=SRC,STATUS='OLD',FORM='FORMATTED')

      READ(1,'(A)') CH
      DO WHILE(1.EQ.1)
       READ(1,'(A)',END=10) CH
       !WRITE(*,*) ISCOMMENT(LINE)
       !IF(ISCOMM) THEN
       ! WRITE(*,*) CH
       !ENDIF
      ENDDO

10    CONTINUE

      READ(1,*) LINE
      WRITE(*,*) TRIM(LINE)
      CLOSE(1)

      RETURN
      END
      !-----------------------------------------------------------------
      DOUBLE PRECISION FUNCTION GETDOUBLE(A)
      DOUBLE PRECISION A
      TEST=A*2.0
      END FUNCTION

      !-----------------------------------------------------------------
      LOGICAL FUNCTION ISCOMMENT(LINE)
      IMPLICIT NONE
      CHARACTER(LEN=*) LINE
      CHARACTER(LEN=1) CH

      CH=LINE(1)
      IF(CH.EQ.'!'.OR.CH.EQ.'#') THEN 
       ISCOMMENT=.TRUE.
      ELSE 
       ISCOMMENT=.FALSE.
      ENDIF
         
      END FUNCTION
      !-----------------------------------------------------------------
