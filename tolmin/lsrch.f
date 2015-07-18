      SUBROUTINE LSRCH (N,X,G,D,XS,GS,RELACC,STEPCB,DDOTG,F,STEP,
     1  NFVALS,NFMAX,GOPT)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(*),G(*),D(*),XS(*),GS(*),GOPT(*)
C
C     Initialization.
C
      RELINT=0.9
      ICOUNT=0
      RATIO=-1.0
      DO 10 I=1,N
      XS(I)=X(I)
      GS(I)=G(I)
      GOPT(I)=G(I)
      IF (D(I) .NE. 0.0) THEN
          TEMP=DABS(X(I)/D(I))
          IF (RATIO .LT. 0.0 .OR. TEMP .LT. RATIO) RATIO=TEMP
      END IF
   10 CONTINUE
      STEP=DMIN1(1.0D0,STEPCB)
C
C     The following number 1.0D-12 is independent of the working
C       accuracy of the computer arithmetic.
C
      STPMIN=DMAX1(RELACC*RATIO,1.0D-12*STEP)
      STEP=DMAX1(STPMIN,STEP)
      SBASE=0.0
      FBASE=F
      DDOTGB=DDOTG
      STPLOW=0.0
      FLOW=F
      DGLOW=DDOTG
      STPHGH=0.0
      STPOPT=0.0
      FOPT=F
      DGOPT=DABS(DDOTG)
C
C     Calculate another function and gradient value.
C
   20 DO 30 I=1,N
   30 X(I)=XS(I)+STEP*D(I)
      CALL FGCALC (N,X,F,G)
      ICOUNT=ICOUNT+1
      DGMID=0.0
      DO 40 I=1,N
   40 DGMID=DGMID+D(I)*G(I)
      IF (F .LE. FOPT) THEN
          IF (F .LT. FOPT .OR. DABS(DGMID) .LT. DGOPT) THEN
              STPOPT=STEP
              FOPT=F
              DO 50 I=1,N
   50         GOPT(I)=G(I)
              DGOPT=DABS(DGMID)
          END IF
      END IF
      IF (NFVALS+ICOUNT .EQ. NFMAX) GOTO 70
C
C      Modify the bounds on the steplength or convergence.
C
      IF (F .GE. FBASE+0.1*(STEP-SBASE)*DDOTGB) THEN
          IF (STPHGH .GT. 0.0 .OR. F .GT. FBASE .OR. DGMID .GT.
     1      0.5*DDOTG) THEN
              STPHGH=STEP
              FHGH=F
              DGHGH=DGMID
              GOTO 60
          END IF
          SBASE=STEP
          FBASE=F
          DDOTGB=DGMID
      END IF
      IF (DGMID .GE. 0.7*DDOTGB) GOTO 70
      STPLOW=STEP
      FLOW=F
      DGLOW=DGMID
   60 IF (STPHGH .GT. 0.0 .AND. STPLOW .GE. RELINT*STPHGH) GOTO 70
C
C     Calculate the next step length or end the iterations.
C
      IF (STPHGH .EQ. 0.0) THEN
          IF (STEP .EQ. STEPCB) GOTO 70
          TEMP=10.0
          IF (DGMID .GT. 0.9*DDOTG) TEMP=DDOTG/(DDOTG-DGMID)
          STEP=DMIN1(TEMP*STEP,STEPCB)
          GOTO 20
      ELSE IF (ICOUNT .EQ. 1 .OR. STPLOW .GT. 0.0) THEN
          DGKNOT=2.0*(FHGH-FLOW)/(STPHGH-STPLOW)-0.5*(DGLOW+DGHGH)
          IF (DGKNOT .GE. 0.0) THEN
              RATIO=DMAX1(0.1D0,0.5D0*DGLOW/(DGLOW-DGKNOT))
          ELSE
              RATIO=(0.5*DGHGH-DGKNOT)/(DGHGH-DGKNOT)
          END IF
          STEP=STPLOW+RATIO*(STPHGH-STPLOW)
          GOTO 20
      ELSE
          STEP=0.1*STEP
          IF (STEP .GE. STPMIN) GOTO 20
      END IF
C
C     Return from subroutine.
C
   70 IF (STEP .NE. STPOPT) THEN
          STEP=STPOPT
          F=FOPT
          DO 80 I=1,N
          X(I)=XS(I)+STEP*D(I)
   80     G(I)=GOPT(I)
      END IF
      NFVALS=NFVALS+ICOUNT
      RETURN
      END
