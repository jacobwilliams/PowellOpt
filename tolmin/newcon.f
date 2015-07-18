      SUBROUTINE NEWCON (N,M,A,IA,IACT,NACT,Z,U,D,RELACC,MDEG,ZZDIAG,
     1  GMNEW,CGRAD)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(IA,*),IACT(*),Z(*),U(*),D(*),ZZDIAG(*),GMNEW(*),
     1  CGRAD(*)
C
C     Initialization.
C
      NP=NACT+1
      KHIGH=MDEG
      IZ=0
      DO 20 I=1,N
      ZZDIAG(I)=0.0
      DO 10 J=NP,N
   10 ZZDIAG(I)=ZZDIAG(I)+Z(IZ+J)**2
   20 IZ=IZ+N
C
C     Calculate the scalar products of D with its constraints.
C
   30 CVMAX=0.0
      DO 50 K=NP,KHIGH
      J=IACT(K)
      IF (J .LE. M) THEN
          SUM=0.0
          SUMABS=0.0
          SUMD=0.0
          DO 40 I=1,N
          TEMP=D(I)*A(I,J)
          SUM=SUM+TEMP
          SUMABS=SUMABS+DABS(TEMP)
   40     SUMD=SUMD+ZZDIAG(I)*A(I,J)**2
      ELSE
          JM=J-M
          IF (JM .LE. N) THEN
              SUM=-D(JM)
          ELSE
              JM=JM-N
              SUM=D(JM)
          END IF
          SUMABS=DABS(SUM)
          SUMD=ZZDIAG(JM)
      END IF
C
C     Pick out the most violated constraint, or return if the
C       violation is negligible.
C
      IF (SUM .GT. RELACC*SUMABS) THEN
          CVIOL=SUM*SUM/SUMD
          IF (CVIOL .GT. CVMAX) THEN
              CVMAX=CVIOL
              IADD=K
              SAVSUM=SUM
              SAVABS=SUMABS
          END IF
      END IF
   50 CONTINUE
      IF (CVMAX .LE. 0.0) GOTO 140
      IF (NACT .EQ. 0) GOTO 120
C
C     Set GMNEW to the gradient of the most violated constraint.
C
      J=IACT(IADD)
      IF (J .LE. M) THEN
          JMV=0
          DO 60 I=1,N
   60     GMNEW(I)=A(I,J)
      ELSE
          JMV=J-M
          DO 70 I=1,N
   70     GMNEW(I)=0.0
          IF (JMV .LE. N) THEN
              GMNEW(JMV)=-1.0
          ELSE
              JMV=JMV-N
              GMNEW(JMV)=1.0
          END IF
      END IF
C
C     Modify GMNEW for the next active constraint.
C
      K=NACT
   80 TEMP=0.0
      IZ=K
      DO 90 I=1,N
      TEMP=TEMP+Z(IZ)*GMNEW(I)
   90 IZ=IZ+N
      TEMP=TEMP*U(K)
      J=IACT(K)
      IF (J .LE. M) THEN
          DO 100 I=1,N
  100     GMNEW(I)=GMNEW(I)-TEMP*A(I,J)
      ELSE
          JM=J-M
          IF (JM .LE. N) THEN
              GMNEW(JM)=GMNEW(JM)+TEMP
          ELSE
              GMNEW(JM-N)=GMNEW(JM-N)-TEMP
          END IF
      END IF
C
C     Revise the values of SAVSUM and SAVABS.
C
      SUM=0.0
      SUMABS=0.0
      DO 110 I=1,N
      TEMP=D(I)*GMNEW(I)
      SUM=SUM+TEMP
  110 SUMABS=SUMABS+DABS(TEMP)
      SAVSUM=DMIN1(SAVSUM,SUM)
      SAVABS=DMAX1(SAVABS,SUMABS)
      K=K-1
      IF (K .GE. 1) GOTO 80
C
C     Add the new constraint to the active set if the constraint
C       violation is still significant.
C
      IF (JMV .GT. 0) D(JMV)=0.0
      IF (SAVSUM .LE. RELACC*SAVABS) GOTO 130
  120 K=NACT
      CALL ADDCON (N,M,A,IA,IACT,NACT,Z,U,RELACC,IADD,GMNEW,CGRAD)
      IF (NACT .GT. K) GOTO 140
C
C     Seek another constraint violation.
C
      IADD=NP
  130 IF (NP .LT. KHIGH) THEN
          K=IACT(KHIGH)
          IACT(KHIGH)=IACT(IADD)
          IACT(IADD)=K
          KHIGH=KHIGH-1
          GOTO 30
      END IF
  140 RETURN
      END
