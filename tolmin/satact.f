      SUBROUTINE SATACT (N,M,A,IA,B,XL,XU,X,IACT,NACT,INFO,Z,U,XBIG,
     1  RELACC,TOL,MEQL)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(IA,*),B(*),XL(*),XU(*),X(*),IACT(*),Z(*),U(*),
     1  XBIG(*)
      IF (NACT .EQ. 0) GOTO 50
      DO 30 K=1,NACT
C
C     Calculate the next constraint residual.
C
      J=IACT(K)
      IF (J .LE. M) THEN
          RES=B(J)
          RESABS=DABS(B(J))
          RESBIG=RESABS
          DO 10 I=1,N
          TEMPA=A(I,J)
          TEMP=TEMPA*X(I)
          RES=RES-TEMP
          RESABS=RESABS+DABS(TEMP)
   10     RESBIG=RESBIG+DABS(TEMPA)*XBIG(I)
      ELSE
          JX=J-M
          IF (JX .LE. N) THEN
              RES=X(JX)-XL(JX)
              RESABS=DABS(X(JX))+DABS(XL(JX))
              RESBIG=XBIG(JX)+DABS(XL(JX))
              SAVEX=XL(JX)
          ELSE
              JX=JX-N
              RES=XU(JX)-X(JX)
              RESABS=DABS(X(JX))+DABS(XU(JX))
              RESBIG=XBIG(JX)+DABS(XU(JX))
              SAVEX=XU(JX)
          END IF
      END IF
C
C     Shift X if necessary.
C
      IF (RES .NE. 0.0) THEN
          TEMP=RES/RESABS
          IF (K .LE. MEQL) TEMP=-DABS(TEMP)
          IF (TOL .EQ. RELACC .OR. TEMP+RELACC .LT. 0.0) THEN
              INFO=1
              SCALE=RES*U(K)
              IZ=K
              DO 20 I=1,N
              X(I)=X(I)+SCALE*Z(IZ)
              IZ=IZ+N
   20         XBIG(I)=DMAX1(XBIG(I),DABS(X(I)))
              IF (J .GT. M) X(JX)=SAVEX
C
C     Else flag a constraint deletion if necessary.
C
          ELSE IF (RES/RESBIG .GT. TOL) THEN
              IACT(K)=-IACT(K)
          END IF
      END IF
   30 CONTINUE
C
C     Delete any flagged constraints and then return.
C
      IDROP=NACT
  40  IF (IACT(IDROP) .LT. 0) THEN
          IACT(IDROP)=-IACT(IDROP)
          CALL DELCON (N,M,A,IA,IACT,NACT,Z,U,RELACC,IDROP)
      END IF
      IDROP=IDROP-1
      IF (IDROP .GT. MEQL) GOTO 40
   50 RETURN
      END
