      SUBROUTINE FGCALC (N,X,F,G)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(*),G(*)
C
C     Calculate the objective function and its gradient.
C
      WA=(X(1)-X(3))**2+(X(2)-X(4))**2
      WB=(X(3)-X(5))**2+(X(4)-X(6))**2
      WC=(X(5)-X(1))**2+(X(6)-X(2))**2
      F=1.0/(WA**8)+1.0/(WB**8)+1.0/(WC**8)
      G(1)=16.0*((X(3)-X(1))/(WA**9)+(X(5)-X(1))/(WC**9))
      G(2)=16.0*((X(4)-X(2))/(WA**9)+(X(6)-X(2))/(WC**9))
      G(3)=16.0*((X(5)-X(3))/(WB**9)+(X(1)-X(3))/(WA**9))
      G(4)=16.0*((X(6)-X(4))/(WB**9)+(X(2)-X(4))/(WA**9))
      G(5)=16.0*((X(1)-X(5))/(WC**9)+(X(3)-X(5))/(WB**9))
      G(6)=16.0*((X(2)-X(6))/(WC**9)+(X(4)-X(6))/(WB**9))
      RETURN
      END

