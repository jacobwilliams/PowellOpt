      SUBROUTINE ADJTOL (N,M,A,IA,B,XL,XU,X,IACT,NACT,XBIG,RELACC,TOL,
     1  MEQL)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(IA,*),B(*),XL(*),XU(*),X(*),IACT(*),XBIG(*)
C
C     Set VIOL to the greatest relative constraint residual of the first
C       NACT constraints.
C
      VIOL=0.0
      IF (NACT .GT. MEQL) THEN
          KL=MEQL+1
          DO 20 K=KL,NACT
          J=IACT(K)
          IF (J .LE. M) THEN
              RES=B(J)
              RESABS=DABS(B(J))
              DO 10 I=1,N
              RES=RES-A(I,J)*X(I)
   10         RESABS=RESABS+DABS(A(I,J)*XBIG(I))
          ELSE
              JM=J-M
              IF (JM .LE. N) THEN
                  RES=X(JM)-XL(JM)
                  RESABS=XBIG(JM)+DABS(XL(JM))
              ELSE
                  JM=JM-N
                  RES=XU(JM)-X(JM)
                  RESABS=XBIG(JM)+DABS(XU(JM))
              END IF
          END IF
          IF (RES .GT. 0.0) VIOL=DMAX1(VIOL,RES/RESABS)
   20     CONTINUE
      END IF
C
C     Adjust TOL.
C
      TOL=0.1*DMIN1(TOL,VIOL)
      IF (TOL .LE. RELACC+RELACC) THEN
          TOL=RELACC
          DO 30 I=1,N
   30     XBIG(I)=DABS(X(I))
      END IF
      RETURN
      END
