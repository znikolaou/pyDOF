      SUBROUTINE PYDOF_INIT(NX,NY,NZ,BCX,BCY,BCZ)
      IMPLICIT NONE
      INTEGER NX,NY,NZ
      CHARACTER(LEN=1) BCX,BCY,BCZ

      CALL PYDOF_READ_COEFFS()
      CALL PYDOF_SET_DIMESION(NX,NY,NZ)
      CALL PYDOF_SET_BOUNDARY_CONDITIONS(NX,NY,NZ,BCX,BCY,BCZ)
      CALL PYDOF_SET_FILTER(NX,NY,NZ) 

      RETURN 
      END
      !-----------------------------------------------------------------
      SUBROUTINE PYDOF_READ_COEFFS(FILENAME)
      IMPLICIT NONE

      RETURN 
      END
      !-----------------------------------------------------------------
      DOUBLE PRECISION FUNCTION GET_1D_COEFFS
      IMPLICIT NONE
      INTEGER NXM=3
      DOUBLE PRECISION COEFFS(-NXM:NXM)

      COEFFS(0)=
      COEFFS(1)=
      COEFFS(2)=

      DO I=1:NXM
       COEFFS(-I)=COEFFS(I)
      ENDDO
      GET_1D_COEFFS=COEFFS
      
      END FUNCTION
      !-----------------------------------------------------------------
      SUBROUTINE SET_DIMENSION(NX,NY,NZ)
      IMPLICIT NONE
      INTEGER NX,NY,NZ,ND
      COMMON/BLOCK_DIMENSION/ ND

      ND=MIN(NX-1,1)+MIN(NY-1,1)+MIN(NZ-1,1)

      END
      !-----------------------------------------------------------------
      DOUBLE PRECISION FUNCTION GET_COEFFS_TENSOR(NX,NY,NZ)
      IMPLICIT NONE
      DOUBLE PRECISION COEFFS(0:NXM,0:NXM,0:NXM)
      COMMON/COEFFS_TENSOR/COEFFS

      !TODO: BUILD COEFFS TENSOR


      END FUNCTION
      !-----------------------------------------------------------------
      SUBROUTINE INITIALISE_FILTER(NX,NY,NZ,BCX,BCY,BCZ)
      INTEGER NX,NY,NZ
      CHARACTER(LEN=1) BCX,BCY,BCZ

      CALL READ_1DCOEFFS()
      CALL SET_DIMESION(NX,NY,NZ)
      CALL SET_BOUNDARY_CONDITIONS(NX,NY,NZ,BCX,BCY,BCZ)
      CALL SET_DFILTER(NX,NY,NZ) 

      END
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
        IL(IC,IS)=MIN(MAX(IL(IC,IS),1),NX)          
       ENDDO
      ENDDO
      
      END FUNCTION
      !-----------------------------------------------------------------
      DOUBLE PRECISION FUNCTION FILTER(Y,NX,NY,NZ,BCX,BCY,BCZ)
      IMPLICIT NONE
      CHARACTER(LEN=1) BCX,BCY,BCZ
      INTEGER NX,NY,NZ
      DOUBLE PRECISION Y(NX,NY,NZ), YFLT(NX,NY,NZ)
     
      !SET DISCRETE FILTER COEFFS FOR GIVEN DIMS
      COEFFS=GET_ALLDIM_COEFFS(NX,NY,NZ)

      !SET BOUNDARY CONDITIONS
      IL=GET_NEIGHBOURS(BCX,NX)
      JL=GET_NEIGHBOURS(BCY,NY)
      KL=GET_NEIGHBOURS(BCZ,NZ)

      YFLT(:,:,:)=0.0D0
      DO K=1,NZ
       DO J=1,NY
        DO I=1,NX
         DO KC=-NZM,NZM
          DO JC=-NYM,NYM
           DO IC=-NXM,NXM
            YFLT(I,J,K)=COEFFS(IC,JC,KC)*Y(IL(IC,I),JL(JC,J),KL(KC,K)) &
                        +YFLT(I,J,K)
           ENDDO
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ENDDO

      FILTER=YFLT
      END FUNCTION
      !-----------------------------------------------------------------
