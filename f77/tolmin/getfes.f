      SUBROUTINE GETFES (N,M,A,IA,B,XL,XU,X,IACT,NACT,PAR,INFO,G,Z,
     1  U,XBIG,RELACC,TOL,MEQL,MSAT,MTOT,BRES,D,ZTG,GM,GMNEW,PARNEW,
     2  CGRAD)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(IA,*),B(*),XL(*),XU(*),X(*),IACT(*),PAR(*),G(*),Z(*),
     1  U(*),XBIG(*),BRES(*),D(*),ZTG(*),GM(*),GMNEW(*),PARNEW(*),
     2  CGRAD(*)
C
C     Make the correction to X for the active constraints.
C
      INFO=0
   10 CALL SATACT (N,M,A,IA,B,XL,XU,X,IACT,NACT,INFO,Z,U,XBIG,RELACC,
     1  TOL,MEQL)
      IF (INFO .GT. 0) MSAT=NACT
      IF (MSAT .EQ. MTOT) GOTO 60
C
C     Try to correct the infeasibility.
C
   20 MSATK=MSAT
      SUMRSK=0.0
   30 CALL CONRES (N,M,A,IA,B,XL,XU,X,IACT,NACT,PAR,G,Z,U,XBIG,BRES,
     1  D,ZTG,RELACC,TOL,STEPCB,SUMRES,MEQL,MSAT,MTOT,INDXBD,GM,GMNEW,
     2  PARNEW,CGRAD)
C
C     Include the new constraint in the active set.
C
      IF (STEPCB .GT. 0.0) THEN
          DO 40 I=1,N
          X(I)=X(I)+STEPCB*D(I)
   40     XBIG(I)=DMAX1(XBIG(I),DABS(X(I)))
          CALL ADDCON (N,M,A,IA,IACT,NACT,Z,U,RELACC,INDXBD,GMNEW,CGRAD)
      END IF
C
C     Test whether to continue the search for feasibility.
C
      IF (MSAT .LT. MTOT) THEN
          IF (STEPCB .EQ. 0.0) GOTO 50
          IF (MSATK .LT. MSAT) GOTO 20
          IF (SUMRSK .EQ. 0.0 .OR. SUMRES .LT. SUMRSK) THEN
              SUMRSK=SUMRES
              ITEST=0
          END IF
          ITEST=ITEST+1
          IF (ITEST .LE. 2) GOTO 30
C
C     Reduce TOL if it may be too large to allow feasibility.
C
   50     IF (TOL .GT. RELACC) THEN
              CALL ADJTOL (N,M,A,IA,B,XL,XU,X,IACT,NACT,XBIG,RELACC,
     1          TOL,MEQL)
              GOTO 10
          END IF
      END IF
   60 RETURN
      END
