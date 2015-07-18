      SUBROUTINE KTVEC (N,M,A,IA,IACT,NACT,PAR,G,RESKT,Z,U,BRES,RELAXF,
     1  MEQL,SSQKT,PARW,RESKTW)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(IA,*),IACT(*),PAR(*),G(*),RESKT(*),Z(*),U(*),
     1  BRES(*),PARW(*),RESKTW(*)
C
C     Calculate the Lagrange parameters and the residual vector.
C
      DO 10 I=1,N
   10 RESKT(I)=G(I)
      IF (NACT .GT. 0) THEN
          ICASE=0
   20     DO 50 KK=1,NACT
          K=NACT+1-KK
          J=IACT(K)
          TEMP=0.0
          IZ=K
          DO 30 I=1,N
          TEMP=TEMP+Z(IZ)*RESKT(I)
   30     IZ=IZ+N
          TEMP=TEMP*U(K)
          IF (ICASE .EQ. 0) PAR(K)=0.0
          IF (K .LE. MEQL .OR. PAR(K)+TEMP .LT. 0.0) THEN
              PAR(K)=PAR(K)+TEMP
          ELSE
              TEMP=-PAR(K)
              PAR(K)=0.0
          END IF
          IF (TEMP .NE. 0.0) THEN
              IF (J .LE. M) THEN
                  DO 40 I=1,N
   40             RESKT(I)=RESKT(I)-TEMP*A(I,J)
              ELSE
                  JM=J-M
                  IF (JM .LE. N) THEN
                      RESKT(JM)=RESKT(JM)+TEMP
                  ELSE
                      RESKT(JM-N)=RESKT(JM-N)-TEMP
                  END IF
              END IF
          END IF
   50     CONTINUE
C
C     Calculate the sum of squares of the KT residual vector.
C
          SSQKT=0.0
          IF (NACT .EQ. N) GOTO 130
          DO 60 I=1,N
   60     SSQKT=SSQKT+RESKT(I)**2
C
C     Apply iterative refinement to the residual vector.
C
          IF (ICASE .EQ. 0) THEN
              ICASE=1
              DO 70 K=1,NACT
   70         PARW(K)=PAR(K)
              DO 80 I=1,N
   80         RESKTW(I)=RESKT(I)
              SSQKTW=SSQKT
              GOTO 20
          END IF
C
C     Undo the iterative refinement if it does not reduce SSQKT.
C
          IF (SSQKTW .LT. SSQKT) THEN
              DO 90 K=1,NACT
   90         PAR(K)=PARW(K)
              DO 100 I=1,N
  100         RESKT(I)=RESKTW(I)
              SSQKT=SSQKTW
          END IF
C
C     Calculate SSQKT when there are no active constraints.
C
      ELSE
          SSQKT=0.0
          DO 110 I=1,N
  110     SSQKT=SSQKT+G(I)**2
      END IF
C
C     Predict the reduction in F if one corrects any positive residuals
C       of active inequality constraints.
C
      RELAXF=0.0
      IF (MEQL. LT. NACT) THEN
          KL=MEQL+1
          DO 120 K=KL,NACT
          J=IACT(K)
          IF (BRES(J) .GT. 0.0) THEN
              RELAXF=RELAXF-PAR(K)*BRES(J)
          END IF
  120     CONTINUE
      END IF
  130 RETURN
      END
