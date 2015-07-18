      SUBROUTINE ADDCON (N,M,A,IA,IACT,NACT,Z,U,RELACC,INDXBD,ZTC,
     1  CGRAD)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(IA,*),IACT(*),Z(*),U(*),ZTC(*),CGRAD(*)
      NP=NACT+1
      ICON=IACT(INDXBD)
      IACT(INDXBD)=IACT(NP)
      IACT(NP)=ICON
C
C     Form ZTC when the new constraint is a bound.
C
      IF (ICON .GT. M) THEN
          INEWBD=ICON-M
          IF (INEWBD .LE. N) THEN
              TEMP=-1.0
          ELSE
              INEWBD=INEWBD-N
              TEMP=1.0
          END IF
          IZNBD=INEWBD*N-N
          DO 10 J=1,N
   10     ZTC(J)=TEMP*Z(IZNBD+J)
C
C     Else form ZTC for an ordinary constraint.
C
      ELSE
          DO 20 I=1,N
   20     CGRAD(I)=A(I,ICON)
          DO 30 J=1,N
          ZTC(J)=0.0
          IZ=J
          DO 30 I=1,N
          ZTC(J)=ZTC(J)+Z(IZ)*CGRAD(I)
   30     IZ=IZ+N
      END IF
C
C     Find any Givens rotations to apply to the last columns of Z.
C
      J=N
   40 JP=J
      J=J-1
      IF (J .GT. NACT) THEN
          IF (ZTC(JP) .EQ. 0.0) GOTO 40
          IF (DABS(ZTC(JP)) .LE. RELACC*DABS(ZTC(J))) THEN
              TEMP=DABS(ZTC(J))
          ELSE IF (DABS(ZTC(J)) .LE. RELACC*DABS(ZTC(JP))) THEN
              TEMP=DABS(ZTC(JP))
          ELSE
              TEMP=DABS(ZTC(JP))*DSQRT(1.0+(ZTC(J)/ZTC(JP))**2)
          END IF
          WCOS=ZTC(J)/TEMP
          WSIN=ZTC(JP)/TEMP
          ZTC(J)=TEMP
C
C     Apply the rotation when the new constraint is a bound.
C
          IZ=J
          IF (ICON .GT. M) THEN
              DO 50 I=1,N
              TEMP=WCOS*Z(IZ+1)-WSIN*Z(IZ)
              Z(IZ)=WCOS*Z(IZ)+WSIN*Z(IZ+1)
              Z(IZ+1)=TEMP
   50         IZ=IZ+N
              Z(IZNBD+JP)=0.0
C
C     Else apply the rotation for an ordinary constraint.
C
          ELSE
              WPIV=0.0
              DO 60 I=1,N
              TEMPA=WCOS*Z(IZ+1)
              TEMPB=WSIN*Z(IZ)
              TEMP=DABS(CGRAD(I))*(DABS(TEMPA)+DABS(TEMPB))
              IF (TEMP .GT. WPIV) THEN
                  WPIV=TEMP
                  IPIV=I
              END IF
              Z(IZ)=WCOS*Z(IZ)+WSIN*Z(IZ+1)
              Z(IZ+1)=TEMPA-TEMPB
   60         IZ=IZ+N
C
C     Ensure orthogonality of Z(.,JP) to CGRAD.
C
              SUM=0.0
              IZ=JP
              DO 70 I=1,N
              SUM=SUM+Z(IZ)*CGRAD(I)
   70         IZ=IZ+N
              IF (SUM .NE. 0.0) THEN
                  IZ=IPIV*N-N+JP
                  Z(IZ)=Z(IZ)-SUM/CGRAD(IPIV)
              END IF
          END IF
          GO TO 40
      END IF
C
C     Test for linear independence in the proposed new active set.
C
      IF (ZTC(NP) .EQ. 0.0) GOTO 90
      IF (ICON .LE. M) THEN
          SUM=0.0
          SUMABS=0.0
          IZ=NP
          DO 80 I=1,N
          TEMP=Z(IZ)*CGRAD(I)
          SUM=SUM+TEMP
          SUMABS=SUMABS+DABS(TEMP)
   80     IZ=IZ+N
          IF (DABS(SUM) .LE. RELACC*SUMABS) GOTO 90
      END IF
C
C     Set the new diagonal element of U and return.
C
      U(NP)=1.0/ZTC(NP)
      NACT=NP
   90 RETURN
      END
