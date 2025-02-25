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
      !DOUBLE PRECISION COEFSL(-NGMX:NGMX,NF)
    
      COEFS_A(:,:,:,:)=0.0D0
           
      IF(NDIM.EQ.3) THEN

      WRITE(*,'(A)') 'PYDOF_SET_COEFS_A: 3D'       
      DO IFLT=1,NF 
      NG=NCOEFM1(IFLT)
       DO K=0,NG
        DO J=0,NG
         DO I=0,NG
          COEFS_A(I,J,K,IFLT)=COEFS(I,IFLT)*COEFS(J,IFLT)*
     &                        COEFS(K,IFLT)
         ENDDO
        ENDDO
       ENDDO
      ENDDO
    
      ELSEIF(NDIM.EQ.2) THEN

      WRITE(*,'(A)') 'PYDOF_SET_COEFS_A: 2D'
      DO IFLT=1,NF 
       NG=NCOEFM1(IFLT)
       DO J=0,NG
        DO I=0,NG
         COEFS_A(I,J,0,IFLT)=COEFS(I,IFLT)*COEFS(J,IFLT)
        ENDDO
       ENDDO
      ENDDO
    
      ELSEIF(NDIM.EQ.1) THEN

      WRITE(*,'(A)') 'PYDOF_SET_COEFS_A: 1D'
      DO IFLT=1,NF 
       NG=NCOEFM1(IFLT)
       DO I=0,NG
        COEFS_A(I,0,0,IFLT)=COEFS(I,IFLT)
       ENDDO
      ENDDO

      ENDIF
    
      RETURN
      END      
      !-----------------------------------------------------------------
      SUBROUTINE PYDOF_FILTER(IFLT,NX,NY,NZ,VAR,FVAR)
      IMPLICIT NONE
      INCLUDE 'pydof_common.h'
      INCLUDE 'omp_lib.h'
      INTEGER IFLT,NX,NY,NZ,I,J,K,TID,NCHECK
      DOUBLE PRECISION VAR(NX,NY,NZ),FVAR(NX,NY,NZ)

!$OMP PARALLEL DEFAULT(SHARED),PRIVATE(TID,NCHECK,I,J,K)
      !TID=OMP_GET_THREAD_NUM()
      !NCHECK=OMP_GET_NUM_THREADS()
      !WRITE(*,'(2(AXI4X))') 'THREAD',TID,'NUM_THREADS=',NCHECK
!$OMP DO SCHEDULE(STATIC),COLLAPSE(3)
      DO K=1,NZ
       DO J=1,NY
        DO I=1,NX                 
         CALL FILTP(IFLT,I,J,K,NX,NY,NZ,VAR,FVAR(I,J,K))
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
      INTEGER IFLT,I,J,K,NX,NY,NZ
      DOUBLE PRECISION V(NX,NY,NZ),FV,DUM
      INTEGER IC,JC,KC,IP,IM,JP,JM,KP,KM
     
      FV=0.0D0
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
      IMPLICIT NONE
      CHARACTER(LEN=2) CBCX
      INTEGER NX,NWRKX,NG,NGX,I,IC,IL(-NG:NG,NWRKX)
      
      IL(-NG:NG,1:NWRKX)=0

      IF(NX.EQ.1) THEN
       IL(-NG:NG,1:NX)=1
       RETURN
      ENDIF

      DO I=1,NX
       DO IC=-NGX,NGX
        IL(IC,I)=I+IC  
        IF(IL(IC,I).LT.1) THEN
         IF(CBCX.EQ.'CP') THEN
          IL(IC,I)=NX-MOD(-IL(IC,I)+1,NX)
         ELSEIF(CBCX.EQ.'DP') THEN
          IL(IC,I)=NX+IL(IC,I)
         ELSEIF(CBCX.EQ.'DI') THEN
          IL(IC,I)=1  
         ENDIF
        ENDIF          
        IF(IL(IC,I).GT.NX) THEN
         IF(CBCX.EQ.'CP') THEN        
          IL(IC,I)=1+MOD(IL(IC,I),NX)
         ELSEIF(CBCX.EQ.'DP') THEN
          IL(IC,I)=IL(IC,I)-NX
         ELSEIF(CBCX.EQ.'DI') THEN
          IL(IC,I)=NX  
         ENDIF
        ENDIF              
        IL(IC,I)=MIN(MAX(IL(IC,I),1),NX)          
       ENDDO
      ENDDO
      
      RETURN
      END
      !-----------------------------------------------------------------
      SUBROUTINE SET_BC_LR(BCL,BCR,NX,NXMX,NGMX,I,IC,IL)
      IMPLICIT NONE
      CHARACTER(LEN=2) BCL,BCR
      INTEGER NX,NXMX,NGMX,I,IC,IL(-NGMX:NGMX,NXMX)
         
      IF(IL(IC,I).GE.1.AND.IL(IC,I).LE.NX) THEN
       RETURN
      ENDIF
      
      IF(IL(IC,I).LT.1) THEN
       IF(BCL.EQ.'CP') THEN
        IL(IC,I)=NX-MOD(-IL(IC,I)+1,NX)
       ELSEIF(BCL.EQ.'DP') THEN
        IL(IC,I)=NX+IL(IC,I)
       ELSEIF(BCL.EQ.'DI') THEN
        IL(IC,I)=1  
       ENDIF
      ENDIF

      IF(IL(IC,I).GT.NX) THEN
       IF(BCR.EQ.'CP') THEN        
        IL(IC,I)=1+MOD(IL(IC,I),NX)
       ELSEIF(BCR.EQ.'DP') THEN
        IL(IC,I)=IL(IC,I)-NX
       ELSEIF(BCR.EQ.'DI') THEN
        IL(IC,I)=NX  
       ENDIF
      ENDIF 

      RETURN
      END     
      !-----------------------------------------------------------------
      SUBROUTINE GET_NEIGH_LR(BC,NX,NXMX,NG,NGMX,IL) 
      IMPLICIT NONE
      CHARACTER(LEN=2) BC(2)
      INTEGER NX,NXMX,NG,NGMX,I,IC,IL(-NGMX:NGMX,NXMX)
      
      IF(NX.EQ.1) THEN
       IL(-NG:NG,1:NX)=1
       RETURN
      ENDIF

      !TODO: OMP 
      DO I=1,NX
       DO IC=-NG,NG
        IL(IC,I)=I+IC 
        CALL SET_BC_LR(BC(1),BC(2),NX,NXMX,NGMX,I,IC,IL)
        !IL(IC,I)=MIN(MAX(IL(IC,I),1),NX)  
       ENDDO
      ENDDO      
         
      RETURN
      END 
      !-----------------------------------------------------------------
      SUBROUTINE PYDOF_SET_BCS(NX,NY,NZ,BCX,BCY,BCZ)
      IMPLICIT NONE
      INCLUDE "pydof_common.h"
      INTEGER NX,NY,NZ,IFLT
      CHARACTER(LEN=2) BCX(2),BCY(2),BCZ(2)

      IL(-NGMX:NGMX,:,:)=0
      JL(-NGMX:NGMX,:,:)=0
      KL(-NGMX:NGMX,:,:)=0

      DO IFLT=1,NF
       CALL GET_NEIGH_LR(BCX,NX,NXMX,NCOEFM1(IFLT),NGMX,IL(:,:,IFLT))
       CALL GET_NEIGH_LR(BCY,NY,NYMX,NCOEFM1(IFLT),NGMX,JL(:,:,IFLT))
       CALL GET_NEIGH_LR(BCZ,NZ,NZMX,NCOEFM1(IFLT),NGMX,KL(:,:,IFLT))
      ENDDO

      RETURN
      END
      !-----------------------------------------------------------------
      SUBROUTINE CHECK_WRK_SIZE(N,NMAX,CN,CNMAX)
      IMPLICIT NONE
      INTEGER N,NMAX
      CHARACTER(LEN=*) CN,CNMAX

      IF(N.GT.NMAX) THEN
       WRITE(*,*) 'ERROR ',CN,'>',CNMAX,':',CN,'=', N,CNMAX,'=',NMAX
       WRITE(*,*) 'INCREASE ARRAY SIZE IN pydof_common.h AND RECOMPILE!'
       WRITE(*,*) 'PYDOF TERMINATING ...'
       STOP
      ENDIF

      END
      !-----------------------------------------------------------------
      SUBROUTINE PYDOF_CHECK_WRK_SIZE(NX,NY,NZ)
      IMPLICIT NONE
      INCLUDE "pydof_common.h"
      INTEGER NX,NY,NZ
      
      CALL CHECK_WRK_SIZE(NX,NXMX,'NX','NXMX')
      CALL CHECK_WRK_SIZE(NY,NYMX,'NY','NYMX')
      CALL CHECK_WRK_SIZE(NZ,NZMX,'NZ','NZMX')
     
      END
      !-----------------------------------------------------------------
      LOGICAL FUNCTION CHECK_BC(BC)
      IMPLICIT NONE
      CHARACTER(LEN=2) BC
      
      IF(BC.EQ.'CP'.OR.BC.EQ.'DP'.OR.BC.EQ.'DI') THEN
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
       WRITE(*,*) 'TERMINATING ...'
       STOP
      ENDIF
      IF(.NOT.CHECK_POSITIVE(NY)) THEN
       WRITE(*,*) MSG,'NY<=0, NY=',NY
       WRITE(*,*) 'TERMINATING ...'
       STOP
      ENDIF
      IF(.NOT.CHECK_POSITIVE(NZ)) THEN
       WRITE(*,*) MSG,'NZ<=0, NZ=',NZ
       WRITE(*,*) 'TERMINATING ...'
       STOP
      ENDIF

      RETURN
      END
      !-----------------------------------------------------------------
      SUBROUTINE PYDOF_CHECK_BCS(BCX,BCY,BCZ)
      IMPLICIT NONE
      CHARACTER(LEN=*) MSG
      PARAMETER(MSG='ERROR! PYDOF_CHECK_BCS:')
      LOGICAL CHECK_BC
      CHARACTER(LEN=2) BCX(2),BCY(2),BCZ(2)
      INTEGER I

      DO I=1,2
       IF(.NOT.CHECK_BC(BCX(I))) THEN
        WRITE(*,*) MSG,'BCX != CP OR DP OR DI: BCX=',BCX(I)  
        WRITE(*,*) 'TERMINATING ...'
        STOP
       ELSE IF(.NOT.CHECK_BC(BCY(I))) THEN   
        WRITE(*,*) MSG,'BCY != CP OR DP OR DI: BCY=',BCY(I)
        WRITE(*,*) 'TERMINATING ...'
        STOP
       ELSE IF(.NOT.CHECK_BC(BCZ(I))) THEN    
        WRITE(*,*) MSG,'BCZ != CP OR DP OR DI: BCZ=',BCZ(I)  
        WRITE(*,*) 'TERMINATING ...'
        STOP
       ENDIF
      ENDDO

      RETURN
      END
      !-----------------------------------------------------------------
      SUBROUTINE PYDOF_SET_OMP_PARAMS(NTHREADS)
      IMPLICIT NONE
      INCLUDE 'omp_lib.h'
      INTEGER NTHREADS,NMX

      NMX=OMP_GET_MAX_THREADS()
      IF(NTHREADS.GT.NMX) THEN
       WRITE(*,*) 'PYDOF_SET_OMP_PARAMS: WARNING'
       WRITE(*,*) 'NTHREADS > SYSTEM MAX THREADS, MAX THREADS SET'
      ENDIF
      CALL OMP_SET_NUM_THREADS(MIN(NMX,NTHREADS))

      RETURN
      END
      !-----------------------------------------------------------------
      SUBROUTINE PYDOF_INIT(NTHREADS,NX,NY,NZ,BCX,BCY,BCZ)
      IMPLICIT NONE
      INCLUDE 'pydof_common.h'
      INTEGER NTHREADS,NX,NY,NZ,NDIM,PYDOF_GET_DIM
      CHARACTER(LEN=2) BCX(2),BCY(2),BCZ(2)
     
      WRITE(*,*) 
      WRITE(*,'(A)') '*********************PYDOF*********************'
      WRITE(*,*) 
      CALL PYDOF_CHECK_DIMS(NX,NY,NZ)
      CALL PYDOF_CHECK_WRK_SIZE(NX,NY,NZ)
      CALL PYDOF_CHECK_BCS(BCX,BCY,BCZ)
      CALL PYDOF_SET_OMP_PARAMS(NTHREADS)
      NDIM=PYDOF_GET_DIM(NX,NY,NZ)

      WRITE(*,'(AXI3)') 'DIMENSION=',NDIM
      WRITE(*,'(AXI3)') 'NO THREADS=',NTHREADS
     
      CALL PYDOF_SET_NCM1OR1(NDIM)
      !TODO
      !CALL PYDOF_CHECK_SUM()
      CALL PYDOF_SET_COEFS_A(NDIM)
      CALL PYDOF_SET_BCS(NX,NY,NZ,BCX,BCY,BCZ)
     
      WRITE(*,'(AXI1)') 'NO OF DISCRETE FILTERS=',NF
      WRITE(*,'(3(AXI5X))') 'MAX ARRAY SIZE',NXMX,'x',NYMX,'x',NZMX
      WRITE(*,*)
      WRITE(*,'(A)') 'DISCRETE FILTER(S) INITIALISED!'
      WRITE(*,*)
      WRITE(*,'(A)') '*********************PYDOF*********************' 
      WRITE(*,*) 
 
      END
      !-----------------------------------------------------------------

