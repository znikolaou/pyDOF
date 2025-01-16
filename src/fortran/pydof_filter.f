      SUBROUTINE PYDOF_SET_NCM1OR1(NDIM)
      IMPLICIT NONE
      INCLUDE 'pydof_common.h'
      INTEGER NDIM
      
      NCM1OR1(1,:)=NCOEFM1(:)
      NCM1OR1(2,:)=1
      NCM1OR1(3,:)=1
      IF(NDIM.EQ.2) THEN
       NCM1OR1(2,:)=NCOEFM1(:)
      ELSEIF(NDIM.EQ.3) THEN
       NCM1OR1(2,:)=NCOEFM1(:)
       NCM1OR1(3,:)=NCOEFM1(:)
      ENDIF

      RETURN
      END
      !-----------------------------------------------------------------
      INTEGER FUNCTION PYDOF_GET_DIM(NX,NY,NZ)
      IMPLICIT NONE
      INTEGER NX,NY,NZ
      
      PYDOF_GET_DIM=MIN(NX-1,1)+MIN(NY-1,1)+MIN(NZ-1,1)
      
      END FUNCTION
      !-----------------------------------------------------------------
      SUBROUTINE PYDOF_SET_COEFS_A(NDIM)
      IMPLICIT NONE
      INCLUDE 'pydof_common.h'
      INTEGER NDIM,I,J,K,IFLT,NG
    
      COEFS_A(:,:,:,:)=0.0E0
      
      IF(NDIM.EQ.3) THEN

      DO IFLT=1,NF 
       NG=NCOEFM1(IFLT)
       DO K=-NG,NG
        DO J=-NG,NG
         DO I=-NG,NG
          COEFS_A(I,J,K,IFLT)=COEFS(I,IFLT)*COEFS(J,IFLT)*COEFS(K,IFLT)
         ENDDO
        ENDDO
       ENDDO
      ENDDO
    
      ELSEIF(NDIM.EQ.2) THEN

      DO IFLT=1,NF 
       NG=NCOEFM1(IFLT)
       DO K=0,0
        DO J=-NG,NG
         DO I=-NG,NG
          COEFS_A(I,J,K,IFLT)=COEFS(I,IFLT)*COEFS(J,IFLT)
         ENDDO
        ENDDO
       ENDDO
      ENDDO
    
      ELSEIF(NDIM.EQ.1) THEN

      DO IFLT=1,NF 
       NG=NCOEFM1(IFLT)
       DO K=0,0
        DO J=0,0
         DO I=-NG,NG
          COEFS_A(I,J,K,IFLT)=COEFS(I,IFLT)
         ENDDO
        ENDDO
       ENDDO
      ENDDO

      ENDIF
    
      RETURN
      END      
      !-----------------------------------------------------------------
      SUBROUTINE PYDOF_FILTER(IFLT,NX,NY,NZ,VAR,FVAR)
      IMPLICIT NONE
      !INCLUDE 'omp_lib.h'
      INCLUDE 'pydof_common.h'
      INTEGER IFLT,NX,NY,NZ,I,J,K,TID,NCHECK
      DOUBLE PRECISION VAR(NX,NY,NZ),FVAR(NX,NY,NZ),FV

       
       !$OMP PARALLEL DEFAULT(SHARED) 
       !$OMP& PRIVATE(TID,NCHECK,I,J,K,FV)
      !TID=OMP_GET_THREAD_NUM()
      !NCHECK=OMP_GET_NUM_THREADS()
      !WRITE(*,'(2(A,X,I3,X))') 'THREAD',TID,'NUM_THREADS=',NCHECK
      !$OMP DO SCHEDULE(STATIC),COLLAPSE(3)
      DO K=1,NZ
       DO J=1,NY
        DO I=1,NX                 
         CALL FILTP(IFLT,I,J,K,NX,NY,NZ,VAR,FV)
         WRITE(*,*) 'test',FV
         FVAR(I,J,K)=FV
        ENDDO
       ENDDO
      ENDDO
      !$OMP ENDDO
      !$OMP END PARALLEL
      
      RETURN
      END
      !----------------------------------------------------------------
      SUBROUTINE FILTP(IFLT,I,J,K,NX,NY,NZ,V,FV)
      !
      ! AUTHOR: Z. NIKOLAOU (2023)
      ! FILTERS ON GIVEN POINT I,J,K USING A SYMMETRIC
      ! FILTER IN X,Y,Z. EXPLOITS FILTER SYMMETRY TO SPEED UP
      ! CALCULATIONS.
      ! 
      ! OPERATION: 
      ! \bar{u}_{i,j,k}=
      ! \sum_{-Mz}^{M_z}c_p\sum_{-My}^{My}b_l\sum_{-Mx}^{Mx}a_r
      !  u_{i+r},{j+l},{k+p}
      ! 
      ! IFLT         FILTER IDENTIFIER 
      ! I,J,K        POINT IN ARRAY V TO GET FILTERED VALUE
      ! NX,NY,NZ     ARRAY SIZE 
      ! NMX,NMY,NMZ  HALF-STENCIL FILTER SIZE IN X,Y,Z
      ! V            ARRAY TO BE FILTERED
      ! FV           FILTER VALUE 
      !
      IMPLICIT NONE
      INCLUDE 'pydof_common.h'
      INTEGER IFLT,I,J,K,NX,NY,NZ,NMX,NMY,NMZ,NCOR1(3)
      DOUBLE PRECISION V(NX,NY,NZ),FV,DUM
      INTEGER IC,JC,KC,IP,IM,JP,JM,KP,KM
     
      FV=0.0E0
      FV=COEFS_A(0,0,0,IFLT)*V(I,J,K)
      !1D-X
      DO IC=1,NCM1OR1(1,IFLT)
       IP=IL(IC,I,IFLT)
       IM=IL(-IC,I,IFLT)
       DUM=V(IP,J,K)+V(IM,J,K) 
       FV=FV+COEFS_A(IC,0,0,IFLT)*DUM
      ENDDO
      !1D-Y
      DO JC=1,NCM1OR1(2,IFLT)
       JP=JL(JC,J,IFLT)
       JM=JL(-JC,J,IFLT)
       DUM=V(I,JP,K)+V(I,JM,K)
       FV=FV+COEFS_A(0,JC,0,IFLT)*DUM
      ENDDO
      !1D-Z
      DO KC=1,NCM1OR1(3,IFLT)
       KP=KL(KC,K,IFLT)
       KM=KL(-KC,K,IFLT)
       DUM=V(I,J,KP)+V(I,J,KM)
       FV=FV+COEFS_A(0,0,KC,IFLT)*DUM
      ENDDO
      !XY
      DO JC=1,NCM1OR1(2,IFLT)
       DO IC=1,NCM1OR1(1,IFLT)
        IP=IL(IC,I,IFLT)
        IM=IL(-IC,I,IFLT)
        JP=JL(JC,J,IFLT)
        JM=JL(-JC,J,IFLT)
        DUM=V(IP,JP,K)+V(IP,JM,K)+V(IM,JP,K)+V(IM,JM,K)
        FV=FV+COEFS_A(IC,JC,0,IFLT)*DUM
       ENDDO
      ENDDO
      !XZ
      DO KC=1,NCM1OR1(3,IFLT)
       DO IC=1,NCM1OR1(1,IFLT)
        IP=IL(IC,I,IFLT)
        IM=IL(-IC,I,IFLT)
        KP=KL(KC,K,IFLT)
        KM=KL(-KC,K,IFLT)
        DUM=V(IP,J,KP)+V(IP,J,KM)+V(IM,J,KP)+V(IM,J,KM)
        FV=FV+COEFS_A(IC,0,KC,IFLT)*DUM
       ENDDO
      ENDDO
      !YZ
      DO KC=1,NCM1OR1(3,IFLT)
       DO JC=1,NCM1OR1(2,IFLT)
        KP=KL(KC,K,IFLT)
        KM=KL(-KC,K,IFLT)
        JP=JL(JC,J,IFLT)
        JM=JL(-JC,J,IFLT)
        DUM=V(I,JP,KP)+V(I,JP,KM)+V(I,JM,KP)+V(I,JM,KM)
        FV=FV+COEFS_A(0,JC,KC,IFLT)*DUM
       ENDDO
      ENDDO
      !XYZ
      DO KC=1,NCM1OR1(3,IFLT)
       DO JC=1,NCM1OR1(2,IFLT)
        DO IC=1,NCM1OR1(1,IFLT)
         KP=KL(KC,K,IFLT)
         KM=KL(-KC,K,IFLT)
         JP=JL(JC,J,IFLT)
         JM=JL(-JC,J,IFLT)
         IP=IL(IC,I,IFLT)
         IM=IL(-IC,I,IFLT)
         DUM=V(IP,JP,KP)+V(IP,JP,KM)+V(IP,JM,KP)+V(IP,JM,KM) +
     &       V(IM,JP,KP)+V(IM,JP,KM)+V(IM,JM,KP)+V(IM,JM,KM)
         FV=FV+COEFS_A(IC,JC,KC,IFLT)*DUM
        ENDDO
       ENDDO
      ENDDO

      RETURN
      END
      !----------------------------------------------------------------
      SUBROUTINE GET_NEIGH(CBCX,NX,NWRKX,NGX,NG,IL) 
      !TODO: REFACTOR
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
          IL(IC,IS)=NX-MOD(-IL(IC,IS)+1,NX)
          !IL(IC,IS)=NX+IL(IC,IS)
         ELSEIF(CBCX.EQ.'C') THEN
          IL(IC,IS)=1  
         ENDIF
        ENDIF          
        IF(IL(IC,IS).GT.NX) THEN
         IF(CBCX.EQ.'P') THEN        
          IL(IC,IS)=1+MOD(IL(IC,IS),NX)
          !IL(IC,IS)=IL(IC,IS)-NX
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
      INTEGER NX,NY,NZ,NDIM,PYDOF_GET_DIM
      CHARACTER(LEN=1) BCX,BCY,BCZ

      WRITE(*,*) 
      WRITE(*,*) '***********PYDOF***********'
      WRITE(*,*) 
      WRITE(*,*) ' Initialising filter ...'
       
      CALL PYDOF_CHECK_DIMS(NX,NY,NZ)
      NDIM=PYDOF_GET_DIM(NX,NY,NZ)
      CALL PYDOF_CHECK_BCS(BCX,BCY,BCZ)
      CALL PYDOF_SET_NCM1OR1(NDIM)
      CALL PYDOF_SET_BCS(NX,NY,NZ,BCX,BCY,BCZ)
      CALL PYDOF_SET_COEFS_A(NDIM)
      
      WRITE(*,*) '***********PYDOF***********'
 
      END
      !-----------------------------------------------------------------
