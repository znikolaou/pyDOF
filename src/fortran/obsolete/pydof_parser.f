      SUBROUTINE PYDOF_READ_COEFS(SRC)
      IMPLICIT NONE
      INCLUDE 'pydof_common.h'
      INTEGER I,J,NM,ISCOMMENT,IU
      PARAMETER(IU=1)
      CHARACTER(LEN=*) SRC
      CHARACTER(LEN=100) CH 
      !

      WRITE(*,*) 'pydof_read_coeffs: Reading coeffs file ...'
      FCOEFS(:)=0.0E0
      OPEN(UNIT=IU,FILE=SRC,STATUS='OLD',FORM='FORMATTED')
       READ(IU,'(100A)') CH
       DO WHILE(TRIM(CH).NE.'#Forward filter')
        READ(IU,'(100A)') CH
       ENDDO
       DO I=1,7
        READ(IU,*) 
       ENDDO
       READ(IU,'(I3.3)') NM
       DO I=1,6
        READ(IU,'(I3,E32.16)') J,FCOEFS(I)
        WRITE(*,*) 'Coeffs',FCOEFS(I)
       ENDDO
      CLOSE(IU)

      WRITE(*,*) 'pydof_read_coeffs: Coeffs loaded!'
 
      RETURN
      END
      !-----------------------------------------------------------------

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
