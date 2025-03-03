!     ----------------------------------
!     This code has been auto-generated:
!     2025-03-01 20:45:40.423638
!     
!     PYDOF 2025.03
!     Author: Z. Nikolaou
!     ----------------------------------
      INTEGER NF,NCMX,NGMX
      PARAMETER (NF=2)
      PARAMETER (NCMX=2)
      PARAMETER (NGMX=1)
      INTEGER NCOEF(NF)
      INTEGER NCOEFM1(NF)
      INTEGER ISFORWARD(NF)
      DOUBLE PRECISION COEFS(0:NGMX,NF)
      DATA NCOEF/
     & 2,
     & 2/
      DATA NCOEFM1/
     & 1,
     & 1/
      DATA ISFORWARD/
     & 1,
     & 0/
      DATA COEFS(:,1)/
     & 5.0000000000000000e-01,
     & 2.4999999999999997e-01/
      DATA COEFS(:,2)/
     & 1.0000000000022471e+00,
     & -1.1235457009206584e-12/
!     
!     COMMON BLOCK DATA
!     
      INTEGER NXMX,NYMX,NZMX
      PARAMETER (NXMX=10000)
      PARAMETER (NYMX=10000)
      PARAMETER (NZMX=10000)
      INTEGER NCM1OR1(3,NF)
      INTEGER IL(-NGMX:NGMX,NXMX,NF)
      INTEGER JL(-NGMX:NGMX,NYMX,NF)
      INTEGER KL(-NGMX:NGMX,NZMX,NF)
      DOUBLE PRECISION COEFS_A(-NGMX:NGMX,-NGMX:NGMX,-NGMX:NGMX,NF)
      COMMON/PYDOF_COMM/ IL,JL,KL,COEFS_A,NCM1OR1
!     ------END------
