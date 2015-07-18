      SUBROUTINE MINFLC (N,M,MEQ,A,IA,B,XL,XU,X,ACC,IACT,NACT,PAR,
     1  IPRINT,INFO,G,Z,U,XBIG,RESKT,BRES,D,ZTG,GM,XS,GS)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(IA,*),B(*),XL(*),XU(*),X(*),IACT(*),PAR(*),G(*),
     1  Z(*),U(*),XBIG(*),RESKT(*),BRES(*),D(*),ZTG(*),GM(*),XS(*),
     2  GS(*)
C
C     Initialize ZZNORM, ITERC, NFVALS and NFMAX.
C
      ZZNORM=-1.0
      ITERC=0
      NFVALS=0
      NFMAX=0
      IF (INFO .GT. 0) NFMAX=INFO
C
C     Check the bounds on N, M and MEQ.
C
      INFO=4
      IF (MAX0(1-N,-M,MEQ*(MEQ-M)) .GT. 0) THEN
          IF (IPRINT .NE. 0) PRINT 1010
 1010     FORMAT (/5X,'ERROR RETURN FROM GETMIN BECAUSE A CONDITION',
     1      ' ON N, M OR MEQ IS VIOLATED')
          GOTO 40
      END IF
C
C     Initialize RELACC, Z, U and TOL.
C
      CALL INITZU (N,M,XL,XU,X,IACT,MEQL,INFO,Z,U,XBIG,RELACC)
      TOL=DMAX1(0.01D0,10.0D0*RELACC)
      IF (INFO .EQ. 4) THEN
          IF (IPRINT .NE. 0) PRINT 1020
 1020     FORMAT (/5X,'ERROR RETURN FROM GETMIN BECAUSE A LOWER',
     1      ' BOUND EXCEEDS AN UPPER BOUND')
          GOTO 40
      END IF
C
C     Add any equality constraints to the active set.
C
      IF (MEQ .GT. 0) THEN
          CALL EQCONS (N,M,MEQ,A,IA,B,XU,IACT,MEQL,INFO,Z,U,RELACC,XS,
     1      GS)
          IF (INFO .EQ. 5) THEN
              IF (IPRINT .NE. 0) PRINT 1030
 1030         FORMAT (/5X,'ERROR RETURN FROM GETMIN BECAUSE THE',
     1          ' EQUALITY CONSTRAINTS ARE INCONSISTENT')
              GOTO 40
          END IF
      END IF
      NACT=MEQL
      MSAT=MEQL
C
C     Add the bounds to the list of constraints.
C
      MTOT=NACT
      DO 10 I=1,N
      IF (XL(I) .LT. XU(I)) THEN
          MTOT=MTOT+2
          IACT(MTOT-1)=M+I
          IACT(MTOT)=M+N+I
      END IF
   10 CONTINUE
C
C     Try to satisfy the bound constraints.
C
      CALL GETFES (N,M,A,IA,B,XL,XU,X,IACT,NACT,PAR,INFO,G,Z,U,XBIG,
     1  RELACC,TOL,MEQL,MSAT,MTOT,BRES,D,ZTG,GM,RESKT,XS,GS)
      IF (MSAT .LT. MTOT) THEN
          IF (IPRINT .NE. 0) PRINT 1040
 1040     FORMAT (/5X,'ERROR RETURN FROM GETMIN BECAUSE THE',
     1      ' EQUALITIES AND BOUNDS ARE INCONSISTENT')
          INFO=6
          GOTO 40
      END IF
C
C     Add the ordinary inequalities to the list of constraints.
C
      IF (M .GT. MEQ) THEN
          MP=MEQ+1
          DO 20 K=MP,M
          MTOT=MTOT+1
   20     IACT(MTOT)=K
      END IF
C
C     Correct any constraint violations.
C
   30 CALL GETFES (N,M,A,IA,B,XL,XU,X,IACT,NACT,PAR,INFO,G,Z,U,XBIG,
     1  RELACC,TOL,MEQL,MSAT,MTOT,BRES,D,ZTG,GM,RESKT,XS,GS)
      IF (MSAT .LT. MTOT) THEN
          IF (IPRINT .NE. 0) PRINT 1050
 1050     FORMAT (/5X,'ERROR RETURN FROM GETMIN BECAUSE THE',
     1      ' CONSTRAINTS ARE INCONSISTENT')
          INFO=7
          GOTO 40
      ELSE IF (MEQL .EQ. N) THEN
          IF (IPRINT .NE. 0) PRINT 1060
 1060     FORMAT (/5X,'GETMIN FINDS THAT THE VARIABLES ARE',
     1      ' DETERMINED BY THE EQUALITY CONSTRAINTS')
          GOTO 40
      END IF
C
C     Minimize the objective function in the case when constraints are
C       treated as degenerate if their residuals are less than TOL.
C
      CALL MINFUN (N,M,A,IA,B,XL,XU,X,ACC,IACT,NACT,PAR,IPRINT,INFO,G,Z,
     1  U,XBIG,RELACC,ZZNORM,TOL,MEQL,MTOT,ITERC,NFVALS,NFMAX,RESKT,
     2  BRES,D,ZTG,GM,XS,GS)
C
C     Reduce TOL if necessary.
C
      IF (TOL .GT. RELACC .AND. NACT .GT. 0) THEN
          IF (NFVALS .NE. NFMAX) THEN
              CALL ADJTOL (N,M,A,IA,B,XL,XU,X,IACT,NACT,XBIG,RELACC,TOL,
     1          MEQL)
              GOTO 30
          ELSE
              INFO=8
          END IF
      END IF
      IF (IPRINT .NE. 0) THEN
          IF (INFO .EQ. 1) PRINT 1070
 1070     FORMAT (/5X,'GETMIN HAS ACHIEVED THE REQUIRED ACCURACY')
          IF (INFO .EQ. 2) PRINT 1080
 1080     FORMAT (/5X,'GETMIN CAN MAKE NO FURTHER PROGRESS BECAUSE',
     1      ' OF ROUNDING ERRORS')
          IF (INFO .EQ. 3) PRINT 1090
 1090     FORMAT (/5X,'GETMIN CAN MAKE NO FURTHER PROGRESS BECAUSE',
     1      ' F WILL NOT DECREASE ANY MORE')
          IF (INFO .EQ. 8) PRINT 1100
 1100     FORMAT (/5X,'GETMIN HAS REACHED THE GIVEN LIMIT ON THE',
     1      ' NUMBER OF CALLS OF FGCALC')
      END IF
   40 RETURN
      END
