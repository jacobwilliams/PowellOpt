      SUBROUTINE MINFUN (N,M,A,IA,B,XL,XU,X,ACC,IACT,NACT,PAR,IPRINT,
     1  INFO,G,Z,U,XBIG,RELACC,ZZNORM,TOL,MEQL,MTOT,ITERC,NFVALS,
     2  NFMAX,RESKT,BRES,D,ZTG,GM,XS,GS)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(IA,*),B(*),XL(*),XU(*),X(*),IACT(*),PAR(*),G(*),Z(*),
     1  U(*),XBIG(*),RESKT(*),BRES(*),D(*),ZTG(*),GM(*),XS(*),GS(*)
      SAVE F
C
C     Initialize the minimization calculation.
C
      MSAT=MTOT
      ITERK=ITERC
      NFVALK=NFVALS
      IF (NFVALS .EQ. 0 .OR. INFO .EQ. 1) THEN
          CALL FGCALC (N,X,F,G)
          NFVALS=NFVALS+1
      END IF
      FPREV=DABS(F+F+1.0)
      ITERP=-1
      IF (IPRINT .NE. 0) THEN
          PRINT 1000, TOL
 1000     FORMAT (/5X,'NEW VALUE OF TOL =',1PD13.5)
          ITERP=ITERC+IABS(IPRINT)
          IF (ITERC .EQ. 0) ITERP=0
      END IF
C
C     Calculate the next search direction.
C
   10 CALL CONRES (N,M,A,IA,B,XL,XU,X,IACT,NACT,PAR,G,Z,U,XBIG,BRES,D,
     1  ZTG,RELACC,TOL,STEPCB,DDOTG,MEQL,MSAT,MTOT,INDXBD,GM,RESKT,XS,
     2  GS)
C
C     Calculate the Kuhn Tucker residual vector.
C
      CALL KTVEC (N,M,A,IA,IACT,NACT,PAR,G,RESKT,Z,U,BRES,RELAXF,MEQL,
     1  SSQKT,XS,GS)
C
C     Test for convergence.
C
      IF (SSQKT .LE. ACC*ACC) THEN
          INFO=1
          GOTO 70
      END IF
      IF (DDOTG .GE. 0.0) THEN
          INFO=2
          GOTO 70
      END IF
C
C     Test for termination due to no decrease in F.
C
      IF (F .GE. FPREV) THEN
          IF (TOL .EQ. RELACC .OR. NACT .EQ. 0) THEN
              IF (DIFF .GT. 0.0) GOTO 20
          END IF
          INFO=3
          GOTO 70
      END IF
   20 DIFF=FPREV-F
      FPREV=F
C
C     Test that more calls of FGCALC are allowed.
C
      IF (NFVALS .EQ. NFMAX) THEN
          INFO=8
          GOTO 70
      END IF
C
C     Test whether to reduce TOL and to provide printing.
C
      IF (TOL .GT. RELACC .AND. ITERC .GT. ITERK .AND.
     1  0.1*RELAXF .GE. DMAX1(DIFF,-0.5D0*DDOTG)) GOTO 70
      IF (ITERP .EQ. ITERC) GOTO 80
C
C     Calculate the step along the search direction.
C
   40 ITERC=ITERC+1
      CALL LSRCH (N,X,G,D,XS,GS,RELACC,STEPCB,DDOTG,F,STEP,NFVALS,
     1  NFMAX,BRES)
      IF (STEP .EQ. 0.0) THEN
          INFO=3
          SUM=0.0
          DO 50 I=1,N
   50     SUM=SUM+DABS(D(I)*GS(I))
          IF (DDOTG+RELACC*SUM .GE. 0.0) INFO=2
          GOTO 70
      END IF
C
C     Revise XBIG.
C
      DO 60 I=1,N
   60 XBIG(I)=DMAX1(XBIG(I),DABS(X(I)))
C
C     Revise the second derivative approximation.
C
      CALL ZBFGS (N,X,NACT,G,Z,ZTG,XS,GS,ZZNORM)
C
C     Add a constraint to the active set if it restricts the step.
C
      IF (STEP .EQ. STEPCB) THEN
          K=IACT(INDXBD)
          IF (K .GT. M) THEN
              K=K-M
              IF (K .LE. N) THEN
                  X(K)=XL(K)
              ELSE
                  X(K-N)=XU(K-N)
              END IF
          END IF
          CALL ADDCON (N,M,A,IA,IACT,NACT,Z,U,RELACC,INDXBD,XS,GS)
      END IF
      GOTO 10
C
C     Printing from the subroutine.
C
   70 IF (IPRINT .EQ. 0) GOTO 90
      ITERK=-1
   80 PRINT 1010, ITERC,NFVALS,F
 1010 FORMAT (/5X,'ITERS =',I4,5X,'F.VALS =',I4,5X,'F =',1PD15.7)
      PRINT 1020, (X(I),I=1,N)
 1020 FORMAT ('  X =',(1P5D14.5))
      PRINT 1030, (G(I),I=1,N)
 1030 FORMAT ('  G =',(1P5D14.5))
      IF (IPRINT .LT. 0) THEN
          IF (NACT .EQ. 0) THEN
              PRINT 1050
 1050         FORMAT (5X,'NO ACTIVE CONSTRAINTS')
          ELSE
              PRINT 1060, (IACT(I),I=1,NACT)
 1060         FORMAT (' IA =',(14I5))
              PRINT 1070, (PAR(I),I=1,NACT)
 1070         FORMAT (' LP =',(1P5D14.5))
          END IF
          IF (NACT .EQ. N) THEN
              PRINT 1080
 1080         FORMAT (5X,'KT RESIDUAL VECTOR IS ZERO')
          ELSE
              PRINT 1090,(RESKT(I),I=1,N)
 1090         FORMAT (' KT =',(1P5D14.5))
          END IF
      END IF
      ITERP=ITERC+IABS(IPRINT)
      IF (ITERK .GE. 0) GOTO 40
   90 RETURN
      END
