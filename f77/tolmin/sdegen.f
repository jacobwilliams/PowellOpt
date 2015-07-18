      SUBROUTINE SDEGEN (N,M,A,IA,IACT,NACT,PAR,Z,U,D,ZTG,GM,RELACC,
     1  DDOTGM,MEQL,MDEG,GMNEW,PARNEW,CGRAD)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(IA,*),IACT(*),PAR(*),Z(*),U(*),D(*),ZTG(*),GM(*),
     1  GMNEW(*),PARNEW(*),CGRAD(*)
      MP=MEQL+1
      DTEST=0.0
C
C     Calculate the search direction and branch if it is not downhill.
C
   10 CALL SDIRN (N,NACT,Z,D,ZTG,GM,RELACC,DDOTGM)
      IF (DDOTGM .EQ. 0.0) GOTO 120
C
C     Branch if there is no need to consider any degenerate constraints.
C     The test gives termination if two consecutive additions to the
C       active set fail to increase the predicted new value of F.
C
      IF (NACT .EQ. MDEG) GOTO 120
      NP=NACT+1
      SUM=0.0
      DO 20 J=NP,N
   20 SUM=SUM+ZTG(J)**2
      IF (DTEST .GT. 0.0 .AND. SUM .GE. DTEST) THEN
          IF (ITEST .EQ. 1) GOTO 120
          ITEST=1
      ELSE
          DTEST=SUM
          ITEST=0
      END IF
C
C     Add a constraint to the active set if there are any significant
C       violations of degenerate constraints.
C
      K=NACT
      CALL NEWCON (N,M,A,IA,IACT,NACT,Z,U,D,RELACC,MDEG,GMNEW,PARNEW,
     1  CGRAD)
      IF (NACT .EQ. K) GOTO 120
      PAR(NACT)=0.0
C
C     Calculate the new reduced gradient and Lagrange parameters.
C
   30 DO 40 I=1,N
   40 GMNEW(I)=GM(I)
      K=NACT
   50 TEMP=0.0
      IZ=K
      DO 60 I=1,N
      TEMP=TEMP+Z(IZ)*GMNEW(I)
   60 IZ=IZ+N
      TEMP=TEMP*U(K)
      PARNEW(K)=PAR(K)+TEMP
      IF (K .EQ. NACT) PARNEW(K)=DMIN1(PARNEW(K),0.0D0)
      J=IACT(K)
      IF (J .LE. M) THEN
          DO 70 I=1,N
   70     GMNEW(I)=GMNEW(I)-TEMP*A(I,J)
      ELSE
          JM=J-M
          IF (JM .LE. N) THEN
              GMNEW(JM)=GMNEW(JM)+TEMP
          ELSE
              GMNEW(JM-N)=GMNEW(JM-N)-TEMP
          END IF
      END IF
      K=K-1
      IF (K .GT. MEQL) GOTO 50
C
C     Set RATIO for linear interpolation between PAR and PARNEW.
C
      RATIO=0.0
      IF (MP .LT. NACT) THEN
          KU=NACT-1
          DO 80 K=MP,KU
          IF (PARNEW(K) .GT. 0.0) THEN
              RATIO=PARNEW(K)/(PARNEW(K)-PAR(K))
              IDROP=K
          END IF
   80     CONTINUE
      END IF
C
C     Apply the linear interpolation.
C
      THETA=1.0-RATIO
      DO 90 K=MP,NACT
   90 PAR(K)=DMIN1(THETA*PARNEW(K)+RATIO*PAR(K),0.0D0)
      DO 100 I=1,N
  100 GM(I)=THETA*GMNEW(I)+RATIO*GM(I)
C
C     Drop a constraint if RATIO is positive.
C
      IF (RATIO .GT. 0.0) THEN
          CALL DELCON (N,M,A,IA,IACT,NACT,Z,U,RELACC,IDROP)
          DO 110 K=IDROP,NACT
  110     PAR(K)=PAR(K+1)
          GOTO 30
      END IF
C
C     Return if there is no freedom for a new search direction.
C
      IF (NACT .LT. N) GOTO 10
      DDOTGM=0.0
  120 RETURN
      END
