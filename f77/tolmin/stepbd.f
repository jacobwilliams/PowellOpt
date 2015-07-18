      SUBROUTINE STEPBD (N,M,A,IA,IACT,BRES,D,STEPCB,DDOTG,MDEG,MSAT,
     1  MTOT,INDXBD)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(IA,*),IACT(*),BRES(*),D(*)
C
C     Set steps to constraint boundaries and find the least positive one.
C
      IFLAG=0
      STEPCB=0.0
      INDXBD=0
      K=MDEG
   10 K=K+1
      IF (K .GT. MTOT) GOTO 40
C
C     Form the scalar product of D with the current constraint normal.
C
   20     J=IACT(K)
          IF (J .LE. M) THEN
              SP=0.0
              DO 30 I=1,N
   30         SP=SP+D(I)*A(I,J)
          ELSE
              JM=J-M
              IF (JM .LE. N) THEN
                  SP=-D(JM)
              ELSE
                  SP=D(JM-N)
              END IF
          END IF
C
C     The next branch is taken if label 20 was reached via label 50.
C
          IF (IFLAG .EQ. 1) GOTO 60
C
C     Set BRES(J) to indicate the status of the j-th constraint.
C
          IF (SP*BRES(J) .LE. 0.0) THEN
              BRES(J)=0.0
          ELSE
              BRES(J)=BRES(J)/SP
              IF (STEPCB .EQ. 0.0 .OR. BRES(J) .LT. STEPCB) THEN
                  STEPCB=BRES(J)
                  INDXBD=K
              END IF
          END IF
          GO TO 10
   40 CONTINUE
C
C     Try to pass through the boundary of a violated constraint.
C
   50 IF (INDXBD .LE. MSAT) GOTO 80
          IFLAG=1
          K=INDXBD
          GOTO 20
   60     MSAT=MSAT+1
          IACT(INDXBD)=IACT(MSAT)
          IACT(MSAT)=J
          BRES(J)=0.0
          INDXBD=MSAT
          DDOTG=DDOTG-SP
          IF (DDOTG .LT. 0.0 .AND. MSAT .LT. MTOT) THEN
C
C     Seek the next constraint boundary along the search direction.
C
              TEMP=0.0
              KL=MDEG+1
              DO 70 K=KL,MTOT
              J=IACT(K)
              IF (BRES(J) .GT. 0.0) THEN
                  IF (TEMP .EQ. 0.0 .OR. BRES(J) .LT. TEMP) THEN
                      TEMP=BRES(J)
                      INDXBD=K
                  END IF
              END IF
   70         CONTINUE
              IF (TEMP .GT. 0.0) THEN
                  STEPCB=TEMP
                  GOTO 50
              END IF
          END IF
   80 CONTINUE
      RETURN
      END
