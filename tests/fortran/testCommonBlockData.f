      PROGRAM MAIN
      INTEGER I
      DOUBLE PRECISION COEFS(3)
      COMMON/COEFS_DATA/ COEFS

      !CALL DATAHOLDER

      DO I=1,3
       WRITE(*,*) COEFS(I), SIZE(COEFS)
      ENDDO

      CALL DATAHOLDERB(COEFS)

      DO I=1,3
       WRITE(*,*) COEFS(I)
      ENDDO

      STOP
      END

      SUBROUTINE DATAHOLDER
      DOUBLE PRECISION COEFS(3)
      DATA COEFS/1.0, 
     &           2.0,
     &           3.0/
      COMMON/COEFS_DATA/COEFS 
 
      END

      SUBROUTINE DATAHOLDERB(COEFS)
      DOUBLE PRECISION COEFS(3)
      COEFS(1)=2.0
      COEFS(2)=3.0
      COEFS(3)=4.0      
      RETURN
      END

