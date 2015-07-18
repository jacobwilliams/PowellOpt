      SUBROUTINE GETD (N,M,A,IA,IACT,NACT,PAR,G,Z,U,D,ZTG,RELACC,
     1  DDOTG,MEQL,MDEG,GM,GMNEW,PARNEW,CGRAD)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(IA,*),IACT(*),PAR(*),G(*),Z(*),U(*),D(*),ZTG(*),
     1  GM(*),GMNEW(*),PARNEW(*),CGRAD(*)
C
C     Initialize GM and cycle backwards through the active set.
C
   10 DO 20 I=1,N
   20 GM(I)=G(I)
      K=NACT
   30 IF (K .GT. 0) THEN
C
C     Set TEMP to the next multiplier, but reduce the active set if
C       TEMP has an unacceptable sign.
C
          TEMP=0.0
          IZ=K
          DO 40 I=1,N
          TEMP=TEMP+Z(IZ)*GM(I)
   40     IZ=IZ+N
          TEMP=TEMP*U(K)
          IF (K .GT. MEQL .AND. TEMP .GT. 0.0) THEN
              CALL DELCON (N,M,A,IA,IACT,NACT,Z,U,RELACC,K)
              GOTO 10
          END IF
C
C     Update GM using the multiplier that has just been calculated.
C
          J=IACT(K)
          IF (J .LE. M) THEN
              DO 50 I=1,N
   50         GM(I)=GM(I)-TEMP*A(I,J)
          ELSE
              JM=J-M
              IF (JM .LE. N) THEN
                  GM(JM)=GM(JM)+TEMP
              ELSE
                  GM(JM-N)=GM(JM-N)-TEMP
              END IF
          END IF
          PAR(K)=TEMP
          K=K-1
          GOTO 30
      END IF
C
C     Calculate the search direction and DDOTG.
C
      DDOTG=0.0
      IF (NACT .LT. N) THEN
          CALL SDEGEN (N,M,A,IA,IACT,NACT,PAR,Z,U,D,ZTG,GM,RELACC,
     1      DDOTGM,MEQL,MDEG,GMNEW,PARNEW,CGRAD)
          IF (DDOTGM .LT. 0.0) THEN
              DO 60 I=1,N
   60         DDOTG=DDOTG+D(I)*G(I)
          END IF
      END IF
      RETURN
      END
