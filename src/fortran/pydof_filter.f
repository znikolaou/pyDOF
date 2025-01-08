      !SUBROUTINE PYDOF_COEFS_DATA
      ! AUTO-GENERATED (PYDOF)
      ! AUTHOR: Z. NIKOLAOU
      !
      !INTEGER NM
      !PARAMETER(NMP1=NM+1)
      !DOUBLE PRECISION COEFS(0:NPM1)
      !DATA COEFS//
      !COMMON/PYDOF_COMM__DATA/NM,COEFS
      !
      !END
      !----------------------------------------------------------------
      !SUBROUTINE SET_DIMENSION(NX,NY,NZ)
      !IMPLICIT NONE
      !INTEGER NX,NY,NZ,ND
      !
      !ND=MIN(NX-1,1)+MIN(NY-1,1)+MIN(NZ-1,1)
      !
      !END
      !-----------------------------------------------------------------
      INTEGER FUNCTION GET_NEIGHBOURS(BC, N, NM)
      IMPLICIT NONE
      CHARACTER(LEN=1) BC
      INTEGER N,NM,I,IC,IS,IL(-NM:NM,N)
  
      IL(-NM:NM,1:N)=0

      IF(N.EQ.1) THEN
       IL(-NM:NM,1:N)=1
       RETURN
      ENDIF

      IS=0
      DO I=1,N
       IS=IS+1
       DO IC=-NM,NM
        IL(IC,IS)=I+IC  
        IF(IL(IC,IS).LT.1) THEN
         IF(BC.EQ.'P') THEN
          IL(IC,IS)=N-MOD(-IL(IC,IS)+1,N)
          !IL(IC,IS)=NX+IL(IC,IS)
         ELSEIF(BC.EQ.'C') THEN
          IL(IC,IS)=1  
         ENDIF
        ENDIF          
        IF(IL(IC,IS).GT.N) THEN
         IF(BC.EQ.'P') THEN        
          IL(IC,IS)=1+MOD(IL(IC,IS),N)
          !IL(IC,IS)=IL(IC,IS)-NX
         ELSEIF(BC.EQ.'C') THEN
          IL(IC,IS)=N  
         ENDIF
        ENDIF              
        IL(IC,IS)=MIN(MAX(IL(IC,IS),1),N)          
       ENDDO
      ENDDO
      
      END FUNCTION
      !-----------------------------------------------------------------
      LOGICAL FUNCTION CHECK_BC(BC)
      IMPLICIT NONE
      CHARACTER(LEN=1) BC
      
      IF(BC.EQ.'C'.OR.BC.EQ.'P') THEN
       CHECK_BC=.TRUE.
      ELSE
       CHECK_BC=.FALSE.
      ENDIF
     
      END FUNCTION
      !-----------------------------------------------------------------
      LOGICAL FUNCTION CHECK_POSITIVE(N)
      IMPLICIT NONE
      INTEGER N

      IF(N<=0) THEN
       CHECK_POSITIVE=.FALSE.
      ELSE
       CHECK_POSITIVE=.TRUE.
      ENDIF

      END FUNCTION 
      !-----------------------------------------------------------------
      SUBROUTINE PYDOF_CHECK_DIMS(NX,NY,NZ)
      IMPLICIT NONE
      CHARACTER(LEN=*) MSG
      PARAMETER(MSG='ERROR! PYDOF_CHECK_DIMS')
      LOGICAL CHECK_POSITIVE
      INTEGER NX,NY,NZ

      IF(.NOT.CHECK_POSITIVE(NX)) THEN
       WRITE(*,*) MSG,'NX<=0, NX=',NX
      ENDIF
      IF(.NOT.CHECK_POSITIVE(NY)) THEN
       WRITE(*,*) MSG,'NY<=0, NY=',NY
      ENDIF
      IF(.NOT.CHECK_POSITIVE(NZ)) THEN
       WRITE(*,*) MSG,'NZ<=0, NZ=',NZ
      ENDIF

      RETURN
      END
      !-----------------------------------------------------------------
      SUBROUTINE PYDOF_CHECK_BCS(BCX,BCY,BCZ)
      IMPLICIT NONE
      CHARACTER(LEN=*) MSG
      PARAMETER(MSG='ERROR! PYDOF_CHECK_BCS:')
      LOGICAL CHECK_BC
      CHARACTER(LEN=1) BCX,BCY,BCZ

      IF(.NOT.CHECK_BC(BCX)) THEN
       WRITE(*,*) MSG,'BCX != C OR P, BCX=',BCX  
      ELSE IF(.NOT.CHECK_BC(BCY)) THEN   
       WRITE(*,*) MSG,'BCY != C OR P, BCY=',BCY
      ELSE IF(.NOT.CHECK_BC(BCZ)) THEN    
       WRITE(*,*) MSG,'BCZ != C OR P, BCZ=',BCZ     
      ENDIF

      RETURN
      END
      !-----------------------------------------------------------------
      SUBROUTINE PYDOF_INIT(NX,NY,NZ,BCX,BCY,BCZ,SRC)
      IMPLICIT NONE
      CHARACTER(LEN=*) SRC
      INTEGER NX,NY,NZ
      CHARACTER(LEN=1) BCX,BCY,BCZ

      WRITE(*,*) 
      WRITE(*,*) '***********PYDOF***********'
      WRITE(*,*) 
      WRITE(*,*) 'Initialising filter ...'
      WRITE(*,*) 'Coefs src:'
      WRITE(*,*) SRC 
      CALL PYDOF_CHECK_DIMS(NX,NY,NZ)
      CALL PYDOF_CHECK_BCS(BCX,BCY,BCZ)
      !CALL PYDOF_READ_COEFS(SRC) 
      !CALL PYDOF_SET_BCS(NX,NY,NZ,BCX,BCY,BCZ,IL,JL,KL)
      !CALL PYDOF_SET_DFILTER(NX,NY,NZ,IL,JL,KL) 

      END
      !-----------------------------------------------------------------
