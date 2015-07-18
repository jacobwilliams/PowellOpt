      SUBROUTINE SDIRN (N,NACT,Z,D,ZTG,GM,RELACC,DDOTGM)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION Z(*),D(*),ZTG(*),GM(*)
      DDOTGM=0.0
      IF (NACT .GE. N) GOTO 60
C
C     Premultiply GM by the transpose of Z.
C
      NP=NACT+1
      DO 20 J=NP,N
      SUM=0.0
      SUMABS=0.0
      IZ=J
      DO 10 I=1,N
      TEMP=Z(IZ)*GM(I)
      SUM=SUM+TEMP
      SUMABS=SUMABS+DABS(TEMP)
   10 IZ=IZ+N
      IF (DABS(SUM) .LE. RELACC*SUMABS) SUM=0.0
   20 ZTG(J)=SUM
C
C     Form D by premultiplying ZTG by -Z.
C
      IZ=0
      DO 40 I=1,N
      SUM=0.0
      SUMABS=0.0
      DO 30 J=NP,N
      TEMP=Z(IZ+J)*ZTG(J)
      SUM=SUM-TEMP
   30 SUMABS=SUMABS+DABS(TEMP)
      IF (DABS(SUM) .LE. RELACC*SUMABS) SUM=0.0
      D(I)=SUM
   40 IZ=IZ+N
C
C     Test that the search direction is downhill.
C
      SUMABS=0.0
      DO 50 I=1,N
      TEMP=D(I)*GM(I)
      DDOTGM=DDOTGM+TEMP
   50 SUMABS=SUMABS+DABS(TEMP)
      IF (DDOTGM+RELACC*SUMABS .GE. 0.0) DDOTGM=0.0
   60 RETURN
      END
