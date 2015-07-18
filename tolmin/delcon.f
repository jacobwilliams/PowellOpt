      SUBROUTINE DELCON (N,M,A,IA,IACT,NACT,Z,U,RELACC,IDROP)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(IA,*),IACT(*),Z(*),U(*)
      NM=NACT-1
      IF (IDROP .EQ. NACT) GOTO 60
      ISAVE=IACT(IDROP)
C
C     Cycle through the constraint exchanges that are needed.
C
      DO 50 J=IDROP,NM
      JP=J+1
      ICON=IACT(JP)
      IACT(J)=ICON
C
C     Calculate the (J,JP) element of R.
C
      IF (ICON .LE. M) THEN
          RJJP=0.0
          IZ=J
          DO 10 I=1,N
          RJJP=RJJP+Z(IZ)*A(I,ICON)
   10     IZ=IZ+N
      ELSE
          IBD=ICON-M
          IF (IBD .LE. N) THEN
              IZBD=IBD*N-N
              RJJP=-Z(IZBD+J)
          ELSE
              IBD=IBD-N
              IZBD=IBD*N-N
              RJJP=Z(IZBD+J)
          END IF
      END IF
C
C     Calculate the parameters of the next rotation.
C
      UJP=U(JP)
      TEMP=RJJP*UJP
      DENOM=DABS(TEMP)
      IF (DENOM*RELACC .LT. 1.0) DENOM=DSQRT(1.0+DENOM*DENOM)
      WCOS=TEMP/DENOM
      WSIN=1.0/DENOM
C
C     Rotate Z when a bound constraint is promoted.
C
      IZ=J
      IF (ICON .GT. M) THEN
          DO 20 I=1,N
          TEMP=WCOS*Z(IZ+1)-WSIN*Z(IZ)
          Z(IZ)=WCOS*Z(IZ)+WSIN*Z(IZ+1)
          Z(IZ+1)=TEMP
   20     IZ=IZ+N
          Z(IZBD+JP)=0.0
C
C     Rotate Z when an ordinary constraint is promoted.
C
      ELSE
          WPIV=0.0
          DO 30 I=1,N
          TEMPA=WCOS*Z(IZ+1)
          TEMPB=WSIN*Z(IZ)
          TEMP=DABS(A(I,ICON))*(DABS(TEMPA)+DABS(TEMPB))
          IF (TEMP .GT. WPIV) THEN
              WPIV=TEMP
              IPIV=I
          END IF
          Z(IZ)=WCOS*Z(IZ)+WSIN*Z(IZ+1)
          Z(IZ+1)=TEMPA-TEMPB
   30     IZ=IZ+N
C
C     Ensure orthogonality to promoted constraint.
C
          SUM=0.0
          IZ=JP
          DO 40 I=1,N
          SUM=SUM+Z(IZ)*A(I,ICON)
   40     IZ=IZ+N
          IF (SUM .NE. 0.0) THEN
              IZ=IPIV*N-N+JP
              Z(IZ)=Z(IZ)-SUM/A(IPIV,ICON)
          END IF
      END IF
C
C     Set the new diagonal elements of U.
C
      U(JP)=-DENOM*U(J)
      U(J)=UJP/DENOM
   50 CONTINUE
C
C     Return.
C
      IACT(NACT)=ISAVE
   60 NACT=NM
      RETURN
      END
