C
C     The pentagon problem.
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(10,15),B(15),XL(6),XU(6),X(6),IACT(27),PAR(20),
     1  W(1000)
C
C     The two values of ICASE provide two different values of ACC, the latter
C     accuracy being so demanding that a return with INFO=2 occurs. The
C     final values of the objective function in the two cases agree well
C     and constraint violations are negligible, considering the differences
C     between the final values of the variables.
C
      IPRINT=10
      IA=10
      N=6
      DO 100 ICASE=1,2
      ACC=1.0D-6
      IF (ICASE .EQ. 2) ACC=1.0D-14
C
C     Set the components of XL, XU and X.
C
      DO 10 I=1,N
      XL(I)=-1.0D6
      XU(I)=1.0D6
   10 X(I)=0.5D0*DFLOAT(I-3)
      X(2)=0.0D0
      X(4)=-1.0D0
      X(6)=1.0D0
C
C     Set the constraints.
C
      M=0
      MEQ=0
      PI=4.0D0*DATAN(1.0D0)
      DO 30 K=1,5
      THETA=0.4D0*DFLOAT(K-1)*PI
      CTH=DCOS(THETA)
      STH=DSIN(THETA)
      DO 30 J=2,N,2
      M=M+1
      DO 20 I=1,N
   20 A(I,M)=0.0D0
      A(J-1,M)=CTH
      A(J,M)=STH
   30 B(M)=1.0D0
C
C     Call the optimization package.
C
      INFO=0
      PRINT 40, ACC,IPRINT
   40 FORMAT (//5X,'CALL OF GETMIN WITH  ACC =',1PD11.4,
     1  '  AND  IPRINT =',I3)
      CALL GETMIN (N,M,MEQ,A,IA,B,XL,XU,X,ACC,IACT,NACT,PAR,IPRINT,
     1  INFO,W)
      PRINT 50, INFO
   50 FORMAT (/5X,'RETURN FROM TOLMIN WITH INFO =',I2)
      CALL FGCALC (N,X,F,W)
      PRINT 60, F
   60 FORMAT (/5X,'FINAL VALUE OF OBJECTIVE FUNCTION =',1PD20.12)
      PRINT 70, (X(I),I=1,N)
   70 FORMAT (/5X,'FINAL COMPONENTS OF X ='//(4X,1P3D20.12))
      DO 80 K=1,M
      DO 80 I=1,N
   80 B(K)=B(K)-A(I,K)*X(I)
      PRINT 90, (B(K),K=1,M)
   90 FORMAT (/5X,'FINAL CONSTRAINT RESIDUALS ='//(3X,1P6D12.4))
  100 CONTINUE
      STOP
      END
