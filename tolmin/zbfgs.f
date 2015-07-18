      SUBROUTINE ZBFGS (N,X,NACT,G,Z,ZTG,XS,GS,ZZNORM)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(*),G(*),Z(*),ZTG(*),XS(*),GS(*)
C
C     Test if there is sufficient convexity for the update.
C
      DD=0.0
      DG=0.0
      TEMP=0.0
      DO 10 I=1,N
      XS(I)=X(I)-XS(I)
      DD=DD+XS(I)**2
      TEMP=TEMP+GS(I)*XS(I)
      GS(I)=G(I)-GS(I)
   10 DG=DG+GS(I)*XS(I)
      IF (DG .LT. 0.1*DABS(TEMP)) GOTO 90
C
C     Transform the Z matrix.
C
      K=N
   20 KP=K
      K=K-1
      IF (K .GT. NACT) THEN
          IF (ZTG(KP) .EQ. 0.0) GOTO 20
          TEMP=DABS(ZTG(KP))*DSQRT(1.0+(ZTG(K)/ZTG(KP))**2)
          WCOS=ZTG(K)/TEMP
          WSIN=ZTG(KP)/TEMP
          ZTG(K)=TEMP
          IZ=K
          DO 30 I=1,N
          TEMP=WCOS*Z(IZ+1)-WSIN*Z(IZ)
          Z(IZ)=WCOS*Z(IZ)+WSIN*Z(IZ+1)
          Z(IZ+1)=TEMP
   30     IZ=IZ+N
          GOTO 20
      END IF
C
C     Update the value of ZZNORM.
C
      IF (ZZNORM .LT. 0.0) THEN
          ZZNORM=DD/DG
      ELSE
          TEMP=DSQRT(ZZNORM*DD/DG)
          ZZNORM=DMIN1(ZZNORM,TEMP)
          ZZNORM=DMAX1(ZZNORM,0.1D0*TEMP)
      END IF
C
C     Complete the updating of Z.
C
      NP=NACT+1
      TEMP=DSQRT(DG)
      IZ=NP
      DO 40 I=1,N
      Z(IZ)=XS(I)/TEMP
   40 IZ=IZ+N
      IF (NP .LT. N) THEN
          KM=NP+1
          DO 80 K=KM,N
          TEMP=0.0
          IZ=K
          DO 50 I=1,N
          TEMP=TEMP+GS(I)*Z(IZ)
   50     IZ=IZ+N
          TEMP=TEMP/DG
          SUM=0.0
          IZ=K
          DO 60 I=1,N
          Z(IZ)=Z(IZ)-TEMP*XS(I)
          SUM=SUM+Z(IZ)**2
   60     IZ=IZ+N
          IF (SUM .LT. ZZNORM) THEN
              TEMP=DSQRT(ZZNORM/SUM)
              IZ=K
              DO 70 I=1,N
              Z(IZ)=TEMP*Z(IZ)
   70         IZ=IZ+N
          END IF
   80     CONTINUE
      END IF
   90 RETURN
      END
