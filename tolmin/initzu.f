      SUBROUTINE INITZU (N,M,XL,XU,X,IACT,MEQL,INFO,Z,U,XBIG,RELACC)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION XL(*),XU(*),X(*),IACT(*),Z(*),U(*),XBIG(*)
C
C     Set RELACC.
C
      ZTPAR=100.0
      RELACC=1.0
   10 RELACC=0.5*RELACC
      TEMPA=ZTPAR+0.5*RELACC
      TEMPB=ZTPAR+RELACC
      IF (ZTPAR .LT. TEMPA .AND. TEMPA .LT. TEMPB) GOTO 10
C
C     Seek bound inconsistencies and bound equality constraints.
C
      MEQL=0
      DO 20 I=1,N
      IF (XL(I) .GT. XU(I)) GOTO 50
      IF (XL(I) .EQ. XU(I)) MEQL=MEQL+1
   20 CONTINUE
C
C     Initialize U, Z and XBIG.
C
      JACT=0
      NN=N*N
      DO 30 I=1,NN
   30 Z(I)=0.0
      IZ=0
      DO 40 I=1,N
      IF (XL(I) .EQ. XU(I)) THEN
          X(I)=XU(I)
          JACT=JACT+1
          U(JACT)=1.0
          IACT(JACT)=I+M+N
          J=JACT
      ELSE
          J=I+MEQL-JACT
      END IF
      Z(IZ+J)=1.0
      IZ=IZ+N
   40 XBIG(I)=DABS(X(I))
      INFO=1
   50 RETURN
      END
