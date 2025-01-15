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
      SUBROUTINE GET_NEIGH(CBCX,NX,NWRKX,NGX,NG,IL) 
      IMPLICIT NONE
      !GATHERS NEIGHBOURING POINTS OF A FUNCTION GIVEN THE NUMBER OF 
      !POINTS LEFT AND RIGHT OF A POINT. THIS IS DONE AT ALL POINTS OF 
      !THE FUNCTION. BOUNDARY CONDITIONS AS BELOW.
      CHARACTER(LEN=1) CBCX
      INTEGER NX,NWRKX,NG,NGX,I,IC,IS
      INTEGER IL(-NG:NG,NWRKX)
      !
  
      IL(-NG:NG,1:NWRKX)=0

      !IF((NX.NE.1).AND(NG.GT.NX)) THEN
      ! WRITE(*,*) 'GET_NEIGH: ERROR, NG > NX',NG,NX
      ! STOP
      !ENDIF
      IF(NX.EQ.1) THEN
       IL(-NG:NG,1:NX)=1
       RETURN
      ENDIF

      IS=0
      DO I=1,NX
       IS=IS+1
       DO IC=-NGX,NGX
        IL(IC,IS)=I+IC  
        IF(IL(IC,IS).LT.1) THEN
         IF(CBCX.EQ.'P') THEN
          !IL(IC,IS)=NX-MOD(-IL(IC,IS)+1,NX)
          IL(IC,IS)=NX+IL(IC,IS)
         ELSEIF(CBCX.EQ.'C') THEN
          IL(IC,IS)=1  
         ENDIF
        ENDIF          
        IF(IL(IC,IS).GT.NX) THEN
         IF(CBCX.EQ.'P') THEN        
          !IL(IC,IS)=1+MOD(IL(IC,IS),NX)
          IL(IC,IS)=IL(IC,IS)-NX
         ELSEIF(CBCX.EQ.'C') THEN
          IL(IC,IS)=NX  
         ENDIF
        ENDIF              
        IL(IC,IS)=MIN(MAX(IL(IC,IS),1),NX)          
        !WRITE(*,*) 'X',I,IC,IL(IC,IS)
       ENDDO
      ENDDO
      
     
      !
      RETURN
      END 
      !-----------------------------------------------------------------
      SUBROUTINE PYDOF_SET_BCS(NX,NY,NZ,BCX,BCY,BCZ)
      IMPLICIT NONE
      INCLUDE "pydof_common.h"
      INTEGER NX,NY,NZ,IFLT
      CHARACTER(LEN=1) BCX,BCY,BCZ

      IL(-NGMX:NGMX,:,:)=0
      JL(-NGMX:NGMX,:,:)=0
      KL(-NGMX:NGMX,:,:)=0

      DO IFLT=1,NF
       CALL GET_NEIGH(BCX,NX,NXMX,NCOEFM1(IFLT),NGMX,IL(:,:,IFLT))
       CALL GET_NEIGH(BCY,NY,NYMX,NCOEFM1(IFLT),NGMX,JL(:,:,IFLT))
       CALL GET_NEIGH(BCZ,NZ,NZMX,NCOEFM1(IFLT),NGMX,KL(:,:,IFLT))
      ENDDO

      RETURN
      END
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
      SUBROUTINE PYDOF_INIT(NX,NY,NZ,BCX,BCY,BCZ)
      IMPLICIT NONE
      INTEGER NX,NY,NZ
      CHARACTER(LEN=1) BCX,BCY,BCZ

      WRITE(*,*) 
      WRITE(*,*) '***********PYDOF***********'
      WRITE(*,*) 
      WRITE(*,*) ' Initialising filter ...'
       
      CALL PYDOF_CHECK_DIMS(NX,NY,NZ)
      CALL PYDOF_CHECK_BCS(BCX,BCY,BCZ)
      CALL PYDOF_SET_BCS(NX,NY,NZ,BCX,BCY,BCZ)
      !CALL PYDOF_SET_DFILTER(NX,NY,NZ,IL,JL,KL)

       
      WRITE(*,*) '***********PYDOF***********'
 
      END
      !-----------------------------------------------------------------
