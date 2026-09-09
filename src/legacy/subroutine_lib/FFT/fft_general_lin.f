chelp+
!
! USAGE: CALL FFT(A, B, NSEG, N, NSPN, ISN, err)
!
! ARRAYS A AND B ORIGINALLY HOLD THE REAL AND IMAGINARY
!      COMPONENTS OF THE DATA, AND RETURN THE REAL AND
!      IMAGINARY COMPONENTS OF THE RESULTING FOURIER COEFFICIENTS.
! MULTIVARIATE DATA IS INDEXED ACCORDING TO THE FORTRAN
!      ARRAY ELEMENT SUCCESSOR FUNCTION, WITHOUT LIMIT
!      ON THE NUMBER OF IMPLIED MULTIPLE SUBSCRIPTS.
!      THE SUBROUTINE IS CALLED ONCE FOR EACH VARIATE.
!      THE CALLS FOR A MULTIVARIATE TRANSFORM MAY BE IN ANY ORDER.
!
! N IS THE DIMENSION OF THE CURRENT VARIABLE.
! NSPN IS THE SPACING OF CONSECUTIVE DATA VALUES
!      WHILE INDEXING THE CURRENT VARIABLE.
! NSEG*N*NSPN IS THE TOTAL NUMBER OF COMPLEX DATA VALUES.
! THE SIGN OF ISN DETERMINES THE SIGN OF THE COMPLEX
!      EXPONENTIAL, AND THE MAGNITUDE OF ISN IS NORMALLY ONE.
!      THE MAGNITUDE OF ISN DETERMINES THE INDEXING INCREMENT FOR A&B.
!
! IF FFT IS CALLED TWICE, WITH OPPOSITE SIGNS ON ISN, AN
!      IDENTITY TRANSFORMATION IS DONE...CALLS CAN BE IN EITHER ORDER.
!      THE RESULTS ARE SCALED BY 1/N WHEN THE SIGN OF ISN IS POSITIVE.
!
! A TRI-VARIATE TRANSFORM WITH A(N1,N2,N3), B(N1,N2,N3)
! IS COMPUTED BY
!        CALL FFT(A,B,N2*N3,N1,1,-1)
!        CALL FFT(A,B,N3,N2,N1,-1)
!        CALL FFT(A,B,1,N3,N1*N2,-1)
!
! A SINGLE-VARIATE TRANSFORM OF N COMPLEX DATA VALUES IS COMPUTED BY
!        CALL FFT(A,B,1,N,1,-1)
!
! THE DATA MAY ALTERNATIVELY BE STORED IN A SINGLE COMPLEX
!      ARRAY A, THEN THE MAGNITUDE OF ISN CHANGED TO TWO TO
!      GIVE THE CORRECT INDEXING INCREMENT AND A(2) USED TO
!      PASS THE INITIAL ADDRESS FOR THE SEQUENCE OF IMAGINARY
!      VALUES, E.G.
!        CALL FFT(A,A(2),NSEG,N,NSPN,-2)
!
! ARRAY NFAC IS WORKING STORAGE FOR FACTORING N.  THE SMALLEST
!      NUMBER EXCEEDING THE 15 LOCATIONS PROVIDED IS 12,754,584.
!
chelp-
C-----------------------------------------------------------------------
C BLOCK DATA:  INITIALIZES LABELED COMMON
C-----------------------------------------------------------------------
C
c      BLOCK DATA
C
c      COMMON /CSTAK/ DSTAK(2500)
C
c      DOUBLE PRECISION DSTAK
c      INTEGER ISTAK(5000)
c      INTEGER ISIZE(5)
C
c      EQUIVALENCE (DSTAK(1),ISTAK(1))
c      EQUIVALENCE (ISTAK(1),LOUT)
c      EQUIVALENCE (ISTAK(2),LNOW)
c      EQUIVALENCE (ISTAK(3),LUSED)
c      EQUIVALENCE (ISTAK(4),LMAX)
c      EQUIVALENCE (ISTAK(5),LBOOK)
c      EQUIVALENCE (ISTAK(6),ISIZE(1))
C
c      DATA ISIZE(1), ISIZE(2), ISIZE(3), ISIZE(4), ISIZE(5) /1,1,1,2,2/
c      DATA LOUT, LNOW, LUSED, LMAX, LBOOK /0,10,10,5000,10/
C
c      END
C
C-----------------------------------------------------------------------
C SUBROUTINE:  FFT
C MULTIVARIATE COMPLEX FOURIER TRANSFORM, COMPUTED IN PLACE
C USING MIXED-RADIX FAST FOURIER TRANSFORM ALGORITHM.
C-----------------------------------------------------------------------
C
      SUBROUTINE FFT(A, B, NSEG, N, NSPN, ISN, err)
C
C
      DIMENSION A(1), B(1), NFAC(15)
C
c      COMMON /CSTAK/ DSTAK(2500)
      COMMON /CSTAK/ DSTAK(2**13)
      DOUBLE PRECISION DSTAK
c      INTEGER ISTAK(5000)
      INTEGER ISTAK(2**14)
c      REAL RSTAK(5000)
      REAL RSTAK(2**14)
C
      EQUIVALENCE (DSTAK(1),ISTAK(1))
      EQUIVALENCE (DSTAK(1),RSTAK(1))

      integer err
      integer isize(5)
      equivalence (istak(1),lout)
      equivalence (istak(2),lnow)
      equivalence (istak(3),lused)
      equivalence (istak(4),lmax)
      equivalence (istak(5),lbook)
      equivalence (istak(6),isize(1))
      data isize(1), isize(2), isize(3), isize(4), isize(5) /1,1,1,2,2/
c      data lout, lnow, lused, lmax, lbook /0,10,10,5000,10/
      data lout, lnow, lused, lmax, lbook /0,10,10,16384,10/

C
C DETERMINE THE FACTORS OF N
C
      M = 0
      NF = IABS(N)
      K = NF
      IF (NF.EQ.1) then
        err = 0
        RETURN
      endif
      NSPAN = IABS(NF*NSPN)
      NTOT = IABS(NSPAN*NSEG)
      IF (ISN*NTOT.NE.0) GO TO 20
      IERR = 6  ! I1MACH(4)
c      WRITE (IERR,9999) NSEG, N, NSPN, ISN
c9999  FORMAT (31H ERROR - ZERO IN FFT PARAMETERS, 4I10)
      err = 213 ! zero in fft parameters
      RETURN
C
  10  M = M + 1
      NFAC(M) = 4
      K = K/16
  20  IF (K-(K/16)*16.EQ.0) GO TO 10
      J = 3
      JJ = 9
      GO TO 40
  30  M = M + 1
      NFAC(M) = J
      K = K/JJ
  40  IF (MOD(K,JJ).EQ.0) GO TO 30
      J = J + 2
      JJ = J**2
      IF (JJ.LE.K) GO TO 40
      IF (K.GT.4) GO TO 50
      KT = M
      NFAC(M+1) = K
      IF (K.NE.1) M = M + 1
      GO TO 90
  50  IF (K-(K/4)*4.NE.0) GO TO 60
      M = M + 1
      NFAC(M) = 2
      K = K/4
C ALL SQUARE FACTORS OUT NOW, BUT K .GE. 5 STILL
  60  KT = M
      MAXP = MAX0(KT+KT+2,K-1)
      J = 2
  70  IF (MOD(K,J).NE.0) GO TO 80
      M = M + 1
      NFAC(M) = J
      K = K/J
  80  J = ((J+1)/2)*2 + 1
      IF (J.LE.K) GO TO 70
  90  IF (M.LE.KT+1) MAXP = M + KT + 1
      IF (M+KT.GT.15) GO TO 120
      IF (KT.EQ.0) GO TO 110
      J = KT
 100  M = M + 1
      NFAC(M) = NFAC(J)
      J = J - 1
      IF (J.NE.0) GO TO 100
C
 110  MAXF = M - KT
      MAXF = NFAC(MAXF)
      IF (KT.GT.0) MAXF = MAX0(NFAC(KT),MAXF)
      J = ISTKGT(MAXF*4,3,err)
      if (err.ne.0) return
      JJ = J + MAXF
      J2 = JJ + MAXF
      J3 = J2 + MAXF
      K = ISTKGT(MAXP,2,err)
      if (err.ne.0) return
      CALL FFTMX(A, B, NTOT, NF, NSPAN, ISN, M, KT, RSTAK(J),
     *    RSTAK(JJ), RSTAK(J2), RSTAK(J3), ISTAK(K), NFAC)
      CALL ISTKRL(2,err)
      RETURN
C
 120  IERR = 6  ! I1MACH(4)
c      WRITE (IERR,9998) N
c9998  FORMAT (50H ERROR - FFT PARAMETER N HAS MORE THAN 15 FACTORS-,
c     *    I20)
      err = 215 ! fft parameter N has more than 15 factors
      RETURN
      END
C
C-----------------------------------------------------------------------
C SUBROUTINE:  FFTMX
C CALLED BY SUBROUTINE 'FFT' TO COMPUTE MIXED-RADIX FOURIER TRANSFORM
C-----------------------------------------------------------------------
C
      SUBROUTINE FFTMX(A, B, NTOT, N, NSPAN, ISN, M, KT, AT, CK, BT,
     *    SK, NP, NFAC)
C
c      DIMENSION A(1), B(1), AT(1), CK(1), BT(1), SK(1), NP(1), NFAC(1)
       DIMENSION A(1), B(1), AT(1), CK(1), BT(1), SK(1), NP(2), NFAC(1)
C
      INC = IABS(ISN)
      NT = INC*NTOT
      KS = INC*NSPAN
      RAD = ATAN(1.0)
      S72 = RAD/0.625
      C72 = COS(S72)
      S72 = SIN(S72)
      S120 = SQRT(0.75)
      IF (ISN.GT.0) GO TO 10
      S72 = -S72
      S120 = -S120
      RAD = -RAD
      GO TO 30
C
C SCALE BY 1/N FOR ISN .GT. 0
C
  10  AK = 1.0/FLOAT(N)
      DO 20 J=1,NT,INC
        A(J) = A(J)*AK
        B(J) = B(J)*AK
  20  CONTINUE
C
  30  KSPAN = KS
      NN = NT - INC
      JC = KS/N
C
C SIN, COS VALUES ARE RE-INITIALIZED EACH LIM STEPS
C
      LIM = 32
      KLIM = LIM*JC
      I = 0
      JF = 0
      MAXF = M - KT
      MAXF = NFAC(MAXF)
      IF (KT.GT.0) MAXF = MAX0(NFAC(KT),MAXF)
C
C COMPUTE FOURIER TRANSFORM
C
  40  DR = 8.0*FLOAT(JC)/FLOAT(KSPAN)
      CD = 2.0*SIN(0.5*DR*RAD)**2
      SD = SIN(DR*RAD)
      KK = 1
      I = I + 1
      IF (NFAC(I).NE.2) GO TO 110
C
C TRANSFORM FOR FACTOR OF 2 (INCLUDING ROTATION FACTOR)
C
      KSPAN = KSPAN/2
      K1 = KSPAN + 2
  50  K2 = KK + KSPAN
      AK = A(K2)
      BK = B(K2)
      A(K2) = A(KK) - AK
      B(K2) = B(KK) - BK
      A(KK) = A(KK) + AK
      B(KK) = B(KK) + BK
      KK = K2 + KSPAN
      IF (KK.LE.NN) GO TO 50
      KK = KK - NN
      IF (KK.LE.JC) GO TO 50
      IF (KK.GT.KSPAN) GO TO 350
  60  C1 = 1.0 - CD
      S1 = SD
      MM = MIN0(K1/2,KLIM)
      GO TO 80
  70  AK = C1 - (CD*C1+SD*S1)
      S1 = (SD*C1-CD*S1) + S1
C
C THE FOLLOWING THREE STATEMENTS COMPENSATE FOR TRUNCATION
C ERROR.  IF ROUNDED ARITHMETIC IS USED, SUBSTITUTE
C C1=AK
C
      C1 = 0.5/(AK**2+S1**2) + 0.5
      S1 = C1*S1
      C1 = C1*AK
  80  K2 = KK + KSPAN
      AK = A(KK) - A(K2)
      BK = B(KK) - B(K2)
      A(KK) = A(KK) + A(K2)
      B(KK) = B(KK) + B(K2)
      A(K2) = C1*AK - S1*BK
      B(K2) = S1*AK + C1*BK
      KK = K2 + KSPAN
      IF (KK.LT.NT) GO TO 80
      K2 = KK - NT
      C1 = -C1
      KK = K1 - K2
      IF (KK.GT.K2) GO TO 80
      KK = KK + JC
      IF (KK.LE.MM) GO TO 70
      IF (KK.LT.K2) GO TO 90
      K1 = K1 + INC + INC
      KK = (K1-KSPAN)/2 + JC
      IF (KK.LE.JC+JC) GO TO 60
      GO TO 40
  90  S1 = FLOAT((KK-1)/JC)*DR*RAD
      C1 = COS(S1)
      S1 = SIN(S1)
      MM = MIN0(K1/2,MM+KLIM)
      GO TO 80
C
C TRANSFORM FOR FACTOR OF 3 (OPTIONAL CODE)
C
 100  K1 = KK + KSPAN
      K2 = K1 + KSPAN
      AK = A(KK)
      BK = B(KK)
      AJ = A(K1) + A(K2)
      BJ = B(K1) + B(K2)
      A(KK) = AK + AJ
      B(KK) = BK + BJ
      AK = -0.5*AJ + AK
      BK = -0.5*BJ + BK
      AJ = (A(K1)-A(K2))*S120
      BJ = (B(K1)-B(K2))*S120
      A(K1) = AK - BJ
      B(K1) = BK + AJ
      A(K2) = AK + BJ
      B(K2) = BK - AJ
      KK = K2 + KSPAN
      IF (KK.LT.NN) GO TO 100
      KK = KK - NN
      IF (KK.LE.KSPAN) GO TO 100
      GO TO 290
C
C TRANSFORM FOR FACTOR OF 4
C
 110  IF (NFAC(I).NE.4) GO TO 230
      KSPNN = KSPAN
      KSPAN = KSPAN/4
 120  C1 = 1.0
      S1 = 0
      MM = MIN0(KSPAN,KLIM)
      GO TO 150
 130  C2 = C1 - (CD*C1+SD*S1)
      S1 = (SD*C1-CD*S1) + S1
C
C THE FOLLOWING THREE STATEMENTS COMPENSATE FOR TRUNCATION
C ERROR.  IF ROUNDED ARITHMETIC IS USED, SUBSTITUTE
C C1=C2
C
      C1 = 0.5/(C2**2+S1**2) + 0.5
      S1 = C1*S1
      C1 = C1*C2
 140  C2 = C1**2 - S1**2
      S2 = C1*S1*2.0
      C3 = C2*C1 - S2*S1
      S3 = C2*S1 + S2*C1
 150  K1 = KK + KSPAN
      K2 = K1 + KSPAN
      K3 = K2 + KSPAN
      AKP = A(KK) + A(K2)
      AKM = A(KK) - A(K2)
      AJP = A(K1) + A(K3)
      AJM = A(K1) - A(K3)
      A(KK) = AKP + AJP
      AJP = AKP - AJP
      BKP = B(KK) + B(K2)
      BKM = B(KK) - B(K2)
      BJP = B(K1) + B(K3)
      BJM = B(K1) - B(K3)
      B(KK) = BKP + BJP
      BJP = BKP - BJP
      IF (ISN.LT.0) GO TO 180
      AKP = AKM - BJM
      AKM = AKM + BJM
      BKP = BKM + AJM
      BKM = BKM - AJM
      IF (S1.EQ.0.0) GO TO 190
 160  A(K1) = AKP*C1 - BKP*S1
      B(K1) = AKP*S1 + BKP*C1
      A(K2) = AJP*C2 - BJP*S2
      B(K2) = AJP*S2 + BJP*C2
      A(K3) = AKM*C3 - BKM*S3
      B(K3) = AKM*S3 + BKM*C3
      KK = K3 + KSPAN
      IF (KK.LE.NT) GO TO 150
 170  KK = KK - NT + JC
      IF (KK.LE.MM) GO TO 130
      IF (KK.LT.KSPAN) GO TO 200
      KK = KK - KSPAN + INC
      IF (KK.LE.JC) GO TO 120
      IF (KSPAN.EQ.JC) GO TO 350
      GO TO 40
 180  AKP = AKM + BJM
      AKM = AKM - BJM
      BKP = BKM - AJM
      BKM = BKM + AJM
      IF (S1.NE.0.0) GO TO 160
 190  A(K1) = AKP
      B(K1) = BKP
      A(K2) = AJP
      B(K2) = BJP
      A(K3) = AKM
      B(K3) = BKM
      KK = K3 + KSPAN
      IF (KK.LE.NT) GO TO 150
      GO TO 170
 200  S1 = FLOAT((KK-1)/JC)*DR*RAD
      C1 = COS(S1)
      S1 = SIN(S1)
      MM = MIN0(KSPAN,MM+KLIM)
      GO TO 140
C
C TRANSFORM FOR FACTOR OF 5 (OPTIONAL CODE)
C
 210  C2 = C72**2 - S72**2
      S2 = 2.0*C72*S72
 220  K1 = KK + KSPAN
      K2 = K1 + KSPAN
      K3 = K2 + KSPAN
      K4 = K3 + KSPAN
      AKP = A(K1) + A(K4)
      AKM = A(K1) - A(K4)
      BKP = B(K1) + B(K4)
      BKM = B(K1) - B(K4)
      AJP = A(K2) + A(K3)
      AJM = A(K2) - A(K3)
      BJP = B(K2) + B(K3)
      BJM = B(K2) - B(K3)
      AA = A(KK)
      BB = B(KK)
      A(KK) = AA + AKP + AJP
      B(KK) = BB + BKP + BJP
      AK = AKP*C72 + AJP*C2 + AA
      BK = BKP*C72 + BJP*C2 + BB
      AJ = AKM*S72 + AJM*S2
      BJ = BKM*S72 + BJM*S2
      A(K1) = AK - BJ
      A(K4) = AK + BJ
      B(K1) = BK + AJ
      B(K4) = BK - AJ
      AK = AKP*C2 + AJP*C72 + AA
      BK = BKP*C2 + BJP*C72 + BB
      AJ = AKM*S2 - AJM*S72
      BJ = BKM*S2 - BJM*S72
      A(K2) = AK - BJ
      A(K3) = AK + BJ
      B(K2) = BK + AJ
      B(K3) = BK - AJ
      KK = K4 + KSPAN
      IF (KK.LT.NN) GO TO 220
      KK = KK - NN
      IF (KK.LE.KSPAN) GO TO 220
      GO TO 290
C
C TRANSFORM FOR ODD FACTORS
C
 230  K = NFAC(I)
      KSPNN = KSPAN
      KSPAN = KSPAN/K
      IF (K.EQ.3) GO TO 100
      IF (K.EQ.5) GO TO 210
      IF (K.EQ.JF) GO TO 250
      JF = K
      S1 = RAD/(FLOAT(K)/8.0)
      C1 = COS(S1)
      S1 = SIN(S1)
      CK(JF) = 1.0
      SK(JF) = 0.0
      J = 1
 240  CK(J) = CK(K)*C1 + SK(K)*S1
      SK(J) = CK(K)*S1 - SK(K)*C1
      K = K - 1
      CK(K) = CK(J)
      SK(K) = -SK(J)
      J = J + 1
      IF (J.LT.K) GO TO 240
 250  K1 = KK
      K2 = KK + KSPNN
      AA = A(KK)
      BB = B(KK)
      AK = AA
      BK = BB
      J = 1
      K1 = K1 + KSPAN
 260  K2 = K2 - KSPAN
      J = J + 1
      AT(J) = A(K1) + A(K2)
      AK = AT(J) + AK
      BT(J) = B(K1) + B(K2)
      BK = BT(J) + BK
      J = J + 1
      AT(J) = A(K1) - A(K2)
      BT(J) = B(K1) - B(K2)
      K1 = K1 + KSPAN
      IF (K1.LT.K2) GO TO 260
      A(KK) = AK
      B(KK) = BK
      K1 = KK
      K2 = KK + KSPNN
      J = 1
 270  K1 = K1 + KSPAN
      K2 = K2 - KSPAN
      JJ = J
      AK = AA
      BK = BB
      AJ = 0.0
      BJ = 0.0
      K = 1
 280  K = K + 1
      AK = AT(K)*CK(JJ) + AK
      BK = BT(K)*CK(JJ) + BK
      K = K + 1
      AJ = AT(K)*SK(JJ) + AJ
      BJ = BT(K)*SK(JJ) + BJ
      JJ = JJ + J
      IF (JJ.GT.JF) JJ = JJ - JF
      IF (K.LT.JF) GO TO 280
      K = JF - J
      A(K1) = AK - BJ
      B(K1) = BK + AJ
      A(K2) = AK + BJ
      B(K2) = BK - AJ
      J = J + 1
      IF (J.LT.K) GO TO 270
      KK = KK + KSPNN
      IF (KK.LE.NN) GO TO 250
      KK = KK - NN
      IF (KK.LE.KSPAN) GO TO 250
C
C MULTIPLY BY ROTATION FACTOR (EXCEPT FOR FACTORS OF 2 AND 4)
C
 290  IF (I.EQ.M) GO TO 350
      KK = JC + 1
 300  C2 = 1.0 - CD
      S1 = SD
      MM = MIN0(KSPAN,KLIM)
      GO TO 320
 310  C2 = C1 - (CD*C1+SD*S1)
      S1 = S1 + (SD*C1-CD*S1)
C
C THE FOLLOWING THREE STATEMENTS COMPENSATE FOR TRUNCATION
C ERROR.  IF ROUNDED ARITHMETIC IS USED, THEY MAY
C BE DELETED.
C
      C1 = 0.5/(C2**2+S1**2) + 0.5
      S1 = C1*S1
      C2 = C1*C2
 320  C1 = C2
      S2 = S1
      KK = KK + KSPAN
 330  AK = A(KK)
      A(KK) = C2*AK - S2*B(KK)
      B(KK) = S2*AK + C2*B(KK)
      KK = KK + KSPNN
      IF (KK.LE.NT) GO TO 330
      AK = S1*S2
      S2 = S1*C2 + C1*S2
      C2 = C1*C2 - AK
      KK = KK - NT + KSPAN
      IF (KK.LE.KSPNN) GO TO 330
      KK = KK - KSPNN + JC
      IF (KK.LE.MM) GO TO 310
      IF (KK.LT.KSPAN) GO TO 340
      KK = KK - KSPAN + JC + INC
      IF (KK.LE.JC+JC) GO TO 300
      GO TO 40
 340  S1 = FLOAT((KK-1)/JC)*DR*RAD
      C2 = COS(S1)
      S1 = SIN(S1)
      MM = MIN0(KSPAN,MM+KLIM)
      GO TO 320
C
C PERMUTE THE RESULTS TO NORMAL ORDER---DONE IN TWO STAGES
C PERMUTATION FOR SQUARE FACTORS OF N
C
 350  NP(1) = KS
      IF (KT.EQ.0) GO TO 440
      K = KT + KT + 1
      IF (M.LT.K) K = K - 1
      J = 1
      NP(K+1) = JC
 360  NP(J+1) = NP(J)/NFAC(J)
      NP(K) = NP(K+1)*NFAC(J)
      J = J + 1
      K = K - 1
      IF (J.LT.K) GO TO 360
      K3 = NP(K+1)
      KSPAN = NP(2)
      KK = JC + 1
      K2 = KSPAN + 1
      J = 1
      IF (N.NE.NTOT) GO TO 400
C
C PERMUTATION FOR SINGLE-VARIATE TRANSFORM (OPTIONAL CODE)
C
 370  AK = A(KK)
      A(KK) = A(K2)
      A(K2) = AK
      BK = B(KK)
      B(KK) = B(K2)
      B(K2) = BK
      KK = KK + INC
      K2 = KSPAN + K2
      IF (K2.LT.KS) GO TO 370
 380  K2 = K2 - NP(J)
      J = J + 1
      K2 = NP(J+1) + K2
      IF (K2.GT.NP(J)) GO TO 380
      J = 1
 390  IF (KK.LT.K2) GO TO 370
      KK = KK + INC
      K2 = KSPAN + K2
      IF (K2.LT.KS) GO TO 390
      IF (KK.LT.KS) GO TO 380
      JC = K3
      GO TO 440
C
C PERMUTATION FOR MULTIVARIATE TRANSFORM
C
 400  K = KK + JC
 410  AK = A(KK)
      A(KK) = A(K2)
      A(K2) = AK
      BK = B(KK)
      B(KK) = B(K2)
      B(K2) = BK
      KK = KK + INC
      K2 = K2 + INC
      IF (KK.LT.K) GO TO 410
      KK = KK + KS - JC
      K2 = K2 + KS - JC
      IF (KK.LT.NT) GO TO 400
      K2 = K2 - NT + KSPAN
      KK = KK - NT + JC
      IF (K2.LT.KS) GO TO 400
 420  K2 = K2 - NP(J)
      J = J + 1
      K2 = NP(J+1) + K2
      IF (K2.GT.NP(J)) GO TO 420
      J = 1
 430  IF (KK.LT.K2) GO TO 400
      KK = KK + JC
      K2 = KSPAN + K2
      IF (K2.LT.KS) GO TO 430
      IF (KK.LT.KS) GO TO 420
      JC = K3
 440  IF (2*KT+1.GE.M) RETURN
      KSPNN = NP(KT+1)
C
C PERMUTATION FOR SQUARE-FREE FACTORS OF N
C
      J = M - KT
      NFAC(J+1) = 1
 450  NFAC(J) = NFAC(J)*NFAC(J+1)
      J = J - 1
      IF (J.NE.KT) GO TO 450
      KT = KT + 1
      NN = NFAC(KT) - 1
      JJ = 0
      J = 0
      GO TO 480
 460  JJ = JJ - K2
      K2 = KK
      K = K + 1
      KK = NFAC(K)
 470  JJ = KK + JJ
      IF (JJ.GE.K2) GO TO 460
      NP(J) = JJ
 480  K2 = NFAC(KT)
      K = KT + 1
      KK = NFAC(K)
      J = J + 1
      IF (J.LE.NN) GO TO 470
C
C DETERMINE THE PERMUTATION CYCLES OF LENGTH GREATER THAN 1
C
      J = 0
      GO TO 500
 490  K = KK
      KK = NP(K)
      NP(K) = -KK
      IF (KK.NE.J) GO TO 490
      K3 = KK
 500  J = J + 1
      KK = NP(J)
      IF (KK.LT.0) GO TO 500
      IF (KK.NE.J) GO TO 490
      NP(J) = -J
      IF (J.NE.NN) GO TO 500
      MAXF = INC*MAXF
C
C REORDER A AND B, FOLLOWING THE PERMUTATION CYCLES
C
      GO TO 570
 510  J = J - 1
      IF (NP(J).LT.0) GO TO 510
      JJ = JC
 520  KSPAN = JJ
      IF (JJ.GT.MAXF) KSPAN = MAXF
      JJ = JJ - KSPAN
      K = NP(J)
      KK = JC*K + I + JJ
      K1 = KK + KSPAN
      K2 = 0
 530  K2 = K2 + 1
      AT(K2) = A(K1)
      BT(K2) = B(K1)
      K1 = K1 - INC
      IF (K1.NE.KK) GO TO 530
 540  K1 = KK + KSPAN
      K2 = K1 - JC*(K+NP(K))
      K = -NP(K)
 550  A(K1) = A(K2)
      B(K1) = B(K2)
      K1 = K1 - INC
      K2 = K2 - INC
      IF (K1.NE.KK) GO TO 550
      KK = K2
      IF (K.NE.J) GO TO 540
      K1 = KK + KSPAN
      K2 = 0
 560  K2 = K2 + 1
      A(K1) = AT(K2)
      B(K1) = BT(K2)
      K1 = K1 - INC
      IF (K1.NE.KK) GO TO 560
      IF (JJ.NE.0) GO TO 520
      IF (J.NE.1) GO TO 510
 570  J = K3 + 1
      NT = NT - KSPNN
      I = NT - INC + 1
      IF (NT.GE.0) GO TO 510
      RETURN
      END
C
C-----------------------------------------------------------------------
C SUBROUTINE:  REALS
C USED WITH 'FFT' TO COMPUTE FOURIER TRANSFORM OR INVERSE FOR REAL DATA
C-----------------------------------------------------------------------
C
      SUBROUTINE REALS(A, B, N, ISN, err)
C
C IF ISN=-1, THIS SUBROUTINE COMPLETES THE FOURIER TRANSFORM
C      OF 2*N REAL DATA VALUES, WHERE THE ORIGINAL DATA VALUES ARE
C      STORED ALTERNATELY IN ARRAYS A AND B, AND ARE FIRST
C      TRANSFORMED BY A COMPLEX FOURIER TRANSFORM OF DIMENSION N.
C      THE COSINE COEFFICIENTS ARE IN A(1),A(2),...A(N),A(N+1)
C      AND THE SINE COEFFICIENTS ARE IN B(1),B(2),...B(N),B(N+1).
C      NOTE THAT THE ARRAYS A AND B MUST HAVE DIMENSION N+1.
C      A TYPICAL CALLING SEQUENCE IS
C        CALL FFT(A,B,N,N,N,-1)
C        CALL REALS(A,B,N,-1)
C
C IF ISN=1, THE INVERSE TRANSFORMATION IS DONE, THE FIRST
C      STEP IN EVALUATING A REAL FOURIER SERIES.
C      A TYPICAL CALLING SEQUENCE IS
C        CALL REALS(A,B,N,1)
C        CALL FFT(A,B,N,N,N,1)
C      THE TIME DOMAIN RESULTS ALTERNATE IN ARRAYS A AND B,
C      I.E. A(1),B(1),A(2),B(2),...A(N),B(N).
C
C THE DATA MAY ALTERNATIVELY BE STORED IN A SINGLE COMPLEX
C      ARRAY A, THEN THE MAGNITUDE OF ISN CHANGED TO TWO TO
C      GIVE THE CORRECT INDEXING INCREMENT AND A(2) USED TO
C      PASS THE INITIAL ADDRESS FOR THE SEQUENCE OF IMAGINARY
C      VALUES, E.G.
C        CALL FFT(A,A(2),N,N,N,-2)
C        CALL REALS(A,A(2),N,-2)
C      IN THIS CASE, THE COSINE AND SINE COEFFICIENTS ALTERNATE IN A.
C
      integer err
      DIMENSION A(1), B(1)
      INC = IABS(ISN)
      NF = IABS(N)
      IF (NF*ISN.NE.0) GO TO 10
      IERR = 6  ! I1MACH(4)
c      WRITE (IERR,9999) N, ISN
c9999  FORMAT (33H ERROR - ZERO IN REALS PARAMETERS, 2I10)
      err = 217 ! zero in reals parameters
      RETURN
C
  10  NK = NF*INC + 2
      NH = NK/2
      RAD = ATAN(1.0)
      DR = -4.0/FLOAT(NF)
      CD = 2.0*SIN(0.5*DR*RAD)**2
      SD = SIN(DR*RAD)
C
C SIN,COS VALUES ARE RE-INITIALIZED EACH LIM STEPS
C
      LIM = 32
      MM = LIM
      ML = 0
      SN = 0.0
      IF (ISN.GT.0) GO TO 40
      CN = 1.0
      A(NK-1) = A(1)
      B(NK-1) = B(1)
  20  DO 30 J=1,NH,INC
        K = NK - J
        AA = A(J) + A(K)
        AB = A(J) - A(K)
        BA = B(J) + B(K)
        BB = B(J) - B(K)
        RE = CN*BA + SN*AB
        EM = SN*BA - CN*AB
        B(K) = (EM-BB)*0.5
        B(J) = (EM+BB)*0.5
        A(K) = (AA-RE)*0.5
        A(J) = (AA+RE)*0.5
        ML = ML + 1
        IF (ML.EQ.MM) GO TO 50
        AA = CN - (CD*CN+SD*SN)
        SN = (SD*CN-CD*SN) + SN
C
C THE FOLLOWING THREE STATEMENTS COMPENSATE FOR TRUNCATION
C ERROR.  IF ROUNDED ARITHMETIC IS USED, SUBSTITUTE
C CN=AA
C
        CN = 0.5/(AA**2+SN**2) + 0.5
        SN = CN*SN
        CN = CN*AA
  30  CONTINUE
      err = 0
      RETURN
C
  40  CN = -1.0
      SD = -SD
      GO TO 20
C
  50  MM = MM + LIM
      SN = FLOAT(ML)*DR*RAD
      CN = COS(SN)
      IF (ISN.GT.0) CN = -CN
      SN = SIN(SN)
      err = 0           !
      return            ! GO TO 30
      END
C
C-----------------------------------------------------------------------
C FUNCTION:  ISTKGT(NITEMS,ITYPE)
C ALLOCATES WORKING STORAGE FOR NITEMS OF ITYPE, AS FOLLOWS
C
C 1 - LOGICAL
C 2 - INTEGER
C 3 - REAL
C 4 - DOUBLE PRECISION
C 5 - COMPLEX
C
C-----------------------------------------------------------------------
C
      INTEGER FUNCTION ISTKGT(NITEMS, ITYPE, err)
C
c      COMMON /CSTAK/ DSTAK(2500)
      COMMON /CSTAK/ DSTAK(2**13)
C
      DOUBLE PRECISION DSTAK
c      INTEGER ISTAK(5000)
      INTEGER ISTAK(2**14)
      INTEGER ISIZE(5)
      integer err
C
      EQUIVALENCE (DSTAK(1),ISTAK(1))
      EQUIVALENCE (ISTAK(1),LOUT)
      EQUIVALENCE (ISTAK(2),LNOW)
      EQUIVALENCE (ISTAK(3),LUSED)
      EQUIVALENCE (ISTAK(4),LMAX)
      EQUIVALENCE (ISTAK(5),LBOOK)
      EQUIVALENCE (ISTAK(6),ISIZE(1))
C
      ISTKGT = (LNOW*ISIZE(2)-1)/ISIZE(ITYPE) + 2
      I = ((ISTKGT-1+NITEMS)*ISIZE(ITYPE)-1)/ISIZE(2) + 3
      IF (I.GT.LMAX) GO TO 10
      ISTAK(I-1) = ITYPE
      ISTAK(I) = LNOW
      LOUT = LOUT + 1
      LNOW = I
      LUSED = MAX0(LUSED,LNOW)
      err = 0
      RETURN
C
  10  IERR = 6  ! I1MACH(4)
c      WRITE (IERR,9999) I
c9999  FORMAT (1H , 39HOVERFLOW OF COMMON ARRAY ISTAK --- NEED, I10)
c      WRITE (IERR,9998) (ISTAK(J),J=1,10), ISTAK(LNOW-1), ISTAK(LNOW)
c9998  FORMAT (1X, I5, 11I6)
      err = 219 ! overflow of common array istak
      return
c      STOP
      END
C
C-----------------------------------------------------------------------
C SUBROUTINE:  ISTKRL(K)
C DE-ALLOCATES THE LAST K WORKING STORAGE AREAS
C-----------------------------------------------------------------------
C
      SUBROUTINE ISTKRL(K,err)
C
c      COMMON /CSTAK/ DSTAK(2500)
      COMMON /CSTAK/ DSTAK(2**13)
C
      integer err
      DOUBLE PRECISION DSTAK
c      INTEGER ISTAK(5000)
      INTEGER ISTAK(2**14)
C
      EQUIVALENCE (DSTAK(1),ISTAK(1))
      EQUIVALENCE (ISTAK(1),LOUT)
      EQUIVALENCE (ISTAK(2),LNOW)
      EQUIVALENCE (ISTAK(3),LUSED)
      EQUIVALENCE (ISTAK(4),LMAX)
      EQUIVALENCE (ISTAK(5),LBOOK)
C
      IN = K
C
      IF (LBOOK.LE.LNOW .AND. LNOW.LE.LUSED .AND. LUSED.LE.LMAX) GO TO
     *    10
      IERR = 6  ! I1MACH(4)
c      WRITE (IERR,9999)
c9999  FORMAT (53H WARNING...ISTAK(2),ISTAK(3),ISTAK(4) OR ISTAK(5) HIT)
c      WRITE (IERR,9997) (ISTAK(J),J=1,10), ISTAK(LNOW-1), ISTAK(LNOW)
      err = 221 ! istak(2),istak(3),istak(4) or istak(5) hit
      return
C
  10  IF (IN.LE.0) then
        err = 0
        RETURN
      endif
      IF (LBOOK.GT.ISTAK(LNOW) .OR. ISTAK(LNOW).GE.LNOW-1) GO TO 20
      LOUT = LOUT - 1
      LNOW = ISTAK(LNOW)
      IN = IN - 1
      GO TO 10
C
  20  IERR = 6  ! I1MACH(4)
c      WRITE (IERR,9998)
c9998  FORMAT (45H WARNING...POINTER AT ISTAK(LNOW) OVERWRITTEN/11X,
c     *    27HDE-ALLOCATION NOT COMPLETED)
c      WRITE (IERR,9997) (ISTAK(J),J=1,10), ISTAK(LNOW-1), ISTAK(LNOW)
c9997  FORMAT (1X, I5, 11I6)
      err = 223 ! pointer at istak(lnow) overwritten - allocation incomplete
      RETURN
C
      END
