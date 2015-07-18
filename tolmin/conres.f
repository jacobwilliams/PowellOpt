      SUBROUTINE CONRES (N,M,A,IA,B,XL,XU,X,IACT,NACT,PAR,G,Z,U,XBIG,
     1  BRES,D,ZTG,RELACC,TOL,STEPCB,SUMRES,MEQL,MSAT,MTOT,INDXBD,
     2  GM,GMNEW,PARNEW,CGRAD)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(IA,*),B(*),XL(*),XU(*),X(*),IACT(*),PAR(*),G(*),
     1  Z(*),U(*),XBIG(*),BRES(*),D(*),ZTG(*),GM(*),GMNEW(*),PARNEW(*),
     2  CGRAD(*)
      IDIFF=MTOT-MSAT
C
C     Calculate and partition the residuals of the inactive constraints,
C       and set the gradient vector when seeking feasibility.
C
      IF (IDIFF .GT. 0.0) THEN
          DO 10 I=1,N
   10     G(I)=0.0
          SUMRES=0.0
      END IF
      MSATK=MSAT
      MDEG=NACT
      MSAT=NACT
      KL=MEQL+1
      DO 50 K=KL,MTOT
      J=IACT(K)
C
C     Calculate the residual of the current constraint.
C
      IF (J .LE. M) THEN
          RES=B(J)
          RESABS=DABS(B(J))
          DO 20 I=1,N
          RES=RES-X(I)*A(I,J)
   20     RESABS=RESABS+DABS(XBIG(I)*A(I,J))
      ELSE
          JM=J-M
          IF (JM .LE. N) THEN
              RES=X(JM)-XL(JM)
              RESABS=DABS(XBIG(JM))+DABS(XL(JM))
          ELSE
              JM=JM-N
              RES=XU(JM)-X(JM)
              RESABS=DABS(XBIG(JM))+DABS(XU(JM))
          END IF
      END IF
      BRES(J)=RES
C
C     Set TEMP to the relative residual.
C
      TEMP=0.0
      IF (RESABS .NE. 0.0) TEMP=RES/RESABS
      IF (K .GT. MSATK .AND. TEMP .LT. 0.0) THEN
          IF (TEMP+RELACC .GE. 0.0) THEN
              IF (J .LE. M) THEN
                  SUM=DABS(B(J))
                  DO 30 I=1,N
   30             SUM=SUM+DABS(X(I)*A(I,J))
              ELSE
                  JM=J-M
                  IF (JM .LE. N) THEN
                      SUM=DABS(X(JM))+DABS(XL(JM))
                  ELSE
                      SUM=DABS(X(JM-N))+DABS(XU(JM-N))
                  END IF
              END IF
              IF (DABS(RES) .LE. SUM*RELACC) TEMP=0.0
          END IF
      END IF
C
C     Place the residual in the appropriate position.
C
      IF (K .LE. NACT) GOTO 50
      IF (K .LE. MSATK .OR. TEMP .GE. 0.0) THEN
          MSAT=MSAT+1
          IF (MSAT .LT. K) THEN
              IACT(K)=IACT(MSAT)
          END IF
          IF (TEMP .GT. TOL) THEN
              IACT(MSAT)=J
          ELSE
              MDEG=MDEG+1
              IACT(MSAT)=IACT(MDEG)
              IACT(MDEG)=J
          END IF
C
C     Update the gradient and SUMRES if the constraint is violated when
C       seeking feasibility.
C
      ELSE
          IF (J .LE. M) THEN
              DO 40 I=1,N
   40         G(I)=G(I)+A(I,J)
          ELSE
              J=J-M
              IF (J .LE. N) THEN
                  G(J)=G(J)-1.0
              ELSE
                  G(J-N)=G(J-N)+1.0
              END IF
          END IF
          SUMRES=SUMRES+DABS(RES)
      END IF
   50 CONTINUE
C
C     Seek the next search direction unless CONRES was called from GETFES
C       and feasibility has been achieved.
C
      STEPCB=0.0
      IF (IDIFF .GT. 0 .AND. MSAT .EQ. MTOT) GOTO 60
      CALL GETD (N,M,A,IA,IACT,NACT,PAR,G,Z,U,D,ZTG,RELACC,DDOTG,MEQL,
     1  MDEG,GM,GMNEW,PARNEW,CGRAD)
C
C     Calculate the (bound on the) step-length due to the constraints.
C
      IF (DDOTG .LT. 0.0) THEN
          CALL STEPBD (N,M,A,IA,IACT,BRES,D,STEPCB,DDOTG,MDEG,MSAT,
     1      MTOT,INDXBD)
      END IF
      IF (IDIFF .EQ. 0) SUMRES=DDOTG
   60 RETURN
      END
