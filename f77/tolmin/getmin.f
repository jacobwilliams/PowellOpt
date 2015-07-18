      SUBROUTINE GETMIN (N,M,MEQ,A,IA,B,XL,XU,X,ACC,IACT,NACT,PAR,
     1  IPRINT,INFO,W)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(IA,*),B(*),XL(*),XU(*),X(*),IACT(*),PAR(*),W(*)
C
C  This is the entry point to a package of subroutines that calculate the
C     the least value of a differentiable function of several variables
C     subject to linear constraints on the values of the variables, using
C     the method that is described in the paper "A tolerant algorithm for
C     linearly constrained optimization calculations", Math. Programming B,
C     Vol. 45, pp. 547-566 (1989).
C
C  N is the number of variables and must be set by the user.
C  M is the number of linear constraints (excluding simple bounds) and
C     must be set by the user.
C  MEQ is the number of constraints that are equalities and must be set
C     by the user.
C  A(.,.) is a 2-dimensional array whose columns are the gradients of the
C     M constraint functions.  Its entries must be set by the user and
C     its dimensions must be at least N by M.
C  IA is the actual first dimension of the array A that is supplied by the
C     user, so its value may not be less than N.
C  B(.) is a vector of constraint right hand sides that must also be set
C     by the user.  Specifically the constraints on the variables X(I)
C     I=1(1)N are
C          A(1,K)*X(1)+...+A(N,K)*X(N) .EQ. B(K)  K=1,...,MEQ
C          A(1,K)*X(1)+...+A(N,K)*X(N) .LE. B(K)  K=MEQ+1,...,M  .
C     Note that the data that define the equality constraints come before
C     the data of the inequalities.
C  XL(.) and XU(.) are vectors whose components must be set to lower and
C     upper bounds on the variables.  Choose very large negative and
C     positive entries if a component should be unconstrained, or set
C     XL(I)=XU(I) to freeze the I-th variable.  Specifically these simple
C     bounds are
C          XL(I) .LE. X(I) and X(I) .LE. XU(I)  I=1,...,N  .
C  X(.) is the vector of variables of the optimization calculation.  Its
C     initial elements must be set by the user to an estimate of the
C     required solution.  The subroutines can usually cope with poor
C     estimates, and there is no need for X(.) to be feasible initially.
C     These variables are adjusted automatically and the values that give
C     the least feasible calculated value of the objective function are
C     available in X(.) on the return from GETMIN.
C  ACC is a tolerance on the first order conditions at the calculated
C     solution of the optimization problem.  These first order conditions
C     state that, if X(.) is a solution, then there is a set of active
C     constraints with indices IACT(K) K=1(1)NACT, say, such that X(.) is
C     on the boundaries of these constraints, and the gradient of the
C     objective function can be expressed in the form
C          GRAD(F)=PAR(1)*GRAD(C(IACT(1)))+...
C                        ...+PAR(NACT)*GRAD(C(IACT(NACT)))  .
C     Here PAR(K) K=1(1)NACT are Lagrange multipliers that are nonpositive
C     for inequality constraints, and GRAD(C(IACT(K))) is the gradient of
C     the IACT(K)-th constraint function, so it is A(.,IACT(K)) if IACT(K)
C     .LE. M, and it is minus or plus the J-th coordinate vector if the
C     constraint is the lower or upper bound on X(J) respectively.  The
C     normal return from the calculation occurs when X(.) is feasible and
C     the sum of squares of components of the vector RESKT(.) is at most
C     ACC**2, where RESKT(.) is the N-component vector of residuals of
C     the first order condition that is displayed above.  Sometimes the
C     package cannot satisfy this condition, because noise in the function
C     values can prevent a change to the variables, no line search being
C     allowed to increase the objective function.
C  IACT(.) is a working space array of integer variables that must be
C     provided by the user.  Its length must be at least (M+2*N).  Its
C     leading entries on the return from the subroutine are the indices
C     IACT(K) K=1(1)NACT that are mentioned in the previous paragraph:
C     in other words they are the indices of the final active constraints.
C     Here the indices M+1,...,M+N and M+N+1,...,M+2*N denote the lower
C     and upper bounds respectively.
C  NACT is set automatically to the integer variable of this ilk that has
C     been introduced already.
C  PAR(.) is a one-dimensional array that will hold the Lagrange
C     multipliers PAR(K) K=1(1)NACT on the return from GETMIN, these
C     parameters being defined in the above paragraph on ACC.  The length
C     of PAR should be at least N.
C  IPRINT must be set by the user to specify the frequency of printing
C     during the execution of the optimization package.  There is no
C     printed output if IPRINT=0.  Otherwise, after ensuring feasibility,
C     information is given every IABS(IPRINT) iterations and whenever a
C     parameter called TOL is reduced.  The printing provides the values
C     of X(.), F(.) and G(.)=GRAD(F) if IPRINT is positive, while if
C     IPRINT is negative this information is augmented by the current
C     values of IACT(K) K=1(1)NACT, PAR(K) K=1(1)NACT and RESKT(I)
C     I=1(1)N.  The reason for returning to the calling program is also
C     displayed when IPRINT is nonzero.
C  INFO is an integer variable that should be set to zero initially,
C     unless the user wishes to impose an upper bound on the number of
C     evaluations of the objective function and its gradient, in which
C     case INFO should be set to the value of this bound.  On the exit
C     from GETMIN it will have one of the following integer values to
C     indicate the reason for leaving the optimization package:
C          INFO=1   X(.) is feasible and the condition that depends on
C     ACC is satisfied.
C          INFO=2   X(.) is feasible and rounding errors are preventing
C     further progress.
C          INFO=3   X(.) is feasible but the objective function fails to
C     decrease although a decrease is predicted by the current gradient
C     vector.  If this return occurs and RESKT(.) has large components
C     then the user's calculation of the gradient of the objective
C     function may be incorrect.  One should also question the coding of
C     the gradient when the final rate of convergence is slow.
C          INFO=4   In this case the calculation cannot begin because IA
C     is less than N or because the lower bound on a variable is greater
C     than the upper bound.
C          INFO=5   This value indicates that the equality constraints
C     are inconsistent.   These constraints include any components of
C     X(.) that are frozen by setting XL(I)=XU(I).
C          INFO=6   In this case there is an error return because the
C     equality constraints and the bounds on the variables are found to
C     be inconsistent.
C          INFO=7   This value indicates that there is no vector of
C     variables that satisfies all of the constraints.  Specifically,
C     when this return or an INFO=6 return occurs, the current active
C     constraints (whose indices are IACT(K) K=1(1)NACT) prevent any
C     change in X(.) that reduces the sum of constraint violations,
C     where only bounds are included in this sum if INFO=6.
C          INFO=8   In this case the limit on the number of calls of
C     subroutine FGCALC (see below) has been reached, and there would
C     have been further calculation otherwise.
C  W(.) is a working space array of real variables that must be provided
C     by the user.  Its length must be at least (M+11*N+N**2).  On exit
C     from the package one can find the final components of GRAD(F) and
C     RESKT(.) in W(1),...,W(N) and W(N+1),...,W(2*N) respectively.
C  Note 1.   The variables N, M, MEQ, IA, ACC and IPRINT and the elements
C     of the arrays A(,.,), B(.), XL(.) and XU(.) are not altered by the
C     optimization procedure.  Their values, the value of INFO and the
C     initial components of X(.) must be set on entry to GETMIN.
C  Note 2.   Of course the package needs the objective function and its
C     gradient.  Therefore the user must provide a subroutine called
C     FGCALC, its first two lines being
C          SUBROUTINE FGCALC (N,X,F,G)
C          DIMENSION X(*),G(*)   .
C     It is called automatically with N set as above and with X(.) set
C     to a feasible vector of variables.  It should calculate the value
C     of the objective function and its gradient for this X(.) and should
C     set them in F and G(I) I=1(1)N respectively, without disturbing N
C     or any of the components of X(.).
C
C     Partition the workspace array.
C
      IG=1
      IRESKT=IG+N
      IZ=IRESKT+N
      IU=IZ+N*N
      IXBIG=IU+N
      IBRES=IXBIG+N
      ID=IBRES+M+N+N
      IZTG=ID+N
      IGM=IZTG+N
      IXS=IGM+N
      IGS=IXS+N
C
C     Call the optimization package.
C
      CALL MINFLC (N,M,MEQ,A,IA,B,XL,XU,X,ACC,IACT,NACT,PAR,IPRINT,
     1  INFO,W(IG),W(IZ),W(IU),W(IXBIG),W(IRESKT),W(IBRES),W(ID),
     2  W(IZTG),W(IGM),W(IXS),W(IGS))
      RETURN
      END
