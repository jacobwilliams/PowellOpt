      SUBROUTINE EQCONS (N,M,MEQ,A,IA,B,XU,IACT,MEQL,INFO,Z,U,RELACC,
     1  AM,CGRAD)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(IA,*),B(*),XU(*),IACT(*),Z(*),U(*),AM(*),CGRAD(*)
C
C     Try to add the next equality constraint to the active set.
C
      DO 50 KEQ=1,MEQ
      IF (MEQL .LT. N) THEN
          NP=MEQL+1
          IACT(NP)=KEQ
          CALL ADDCON (N,M,A,IA,IACT,MEQL,Z,U,RELACC,NP,AM,CGRAD)
          IF (MEQL .EQ. NP) GOTO 50
      END IF
C
C     If linear dependence occurs then find the multipliers of the
C       dependence relation and apply them to the right hand sides.
C
      SUM=B(KEQ)
      SUMABS=DABS(B(KEQ))
      IF (MEQL .GT. 0) THEN
          DO 10 I=1,N
   10     AM(I)=A(I,KEQ)
          K=MEQL
   20     VMULT=0.0
          IZ=K
          DO 30 I=1,N
          VMULT=VMULT+Z(IZ)*AM(I)
   30     IZ=IZ+N
          VMULT=VMULT*U(K)
          J=IACT(K)
          IF (J .LE. M) THEN
              DO 40 I=1,N
   40         AM(I)=AM(I)-VMULT*A(I,J)
              RHS=B(J)
          ELSE
              JM=J-M-N
              AM(JM)=AM(JM)-VMULT
              RHS=XU(JM)
          END IF
          SUM=SUM-RHS*VMULT
          SUMABS=SUMABS+DABS(RHS*VMULT)
          K=K-1
          IF (K .GE. 1) GOTO 20
      END IF
C
C     Error return if the constraints are inconsistent.
C
      IF (DABS(SUM) .GT. RELACC*SUMABS) THEN
          INFO=5
          GOTO 60
      END IF
   50 CONTINUE
   60 RETURN
      END
