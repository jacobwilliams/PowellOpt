    module tolmin_module
    
    private

    abstract interface
        subroutine func(n,x,f,g)  !! FGCALC interface
        implicit none
        integer :: n
        real * 8 :: x(*)
        real * 8 :: f
        real * 8 :: g(*)
        end subroutine func
    end interface
    
    public :: getmin
    public :: tolmin_test
    
    contains

Subroutine getmin (n, m, meq, a, ia, b, xl, xu, x, acc, iact, nact, &
  par, iprint, info, w, fgcalc)
      Implicit real * 8 (a-h, o-z)
      Dimension a(ia,*), b(*), xl(*), xu(*), x(*), iact(*), par(*), w(*)
      procedure(func) :: fgcalc
!
!  This is the entry point to a package of subroutines that calculate the
!     the least value of a differentiable function of several variables
!     subject to linear constraints on the values of the variables, using
!     the method that is described in the paper "A tolerant algorithm for
!     linearly constrained optimization calculations", Math. Programming B,
!     Vol. 45, pp. 547-566 (1989).
!
!  N is the number of variables and must be set by the user.
!  M is the number of linear constraints (excluding simple bounds) and
!     must be set by the user.
!  MEQ is the number of constraints that are equalities and must be set
!     by the user.
!  A(.,.) is a 2-dimensional array whose columns are the gradients of the
!     M constraint functions.  Its entries must be set by the user and
!     its dimensions must be at least N by M.
!  IA is the actual first dimension of the array A that is supplied by the
!     user, so its value may not be less than N.
!  B(.) is a vector of constraint right hand sides that must also be set
!     by the user.  Specifically the constraints on the variables X(I)
!     I=1(1)N are
!          A(1,K)*X(1)+...+A(N,K)*X(N) .EQ. B(K)  K=1,...,MEQ
!          A(1,K)*X(1)+...+A(N,K)*X(N) .LE. B(K)  K=MEQ+1,...,M  .
!     Note that the data that define the equality constraints come before
!     the data of the inequalities.
!  XL(.) and XU(.) are vectors whose components must be set to lower and
!     upper bounds on the variables.  Choose very large negative and
!     positive entries if a component should be unconstrained, or set
!     XL(I)=XU(I) to freeze the I-th variable.  Specifically these simple
!     bounds are
!          XL(I) .LE. X(I) and X(I) .LE. XU(I)  I=1,...,N  .
!  X(.) is the vector of variables of the optimization calculation.  Its
!     initial elements must be set by the user to an estimate of the
!     required solution.  The subroutines can usually cope with poor
!     estimates, and there is no need for X(.) to be feasible initially.
!     These variables are adjusted automatically and the values that give
!     the least feasible calculated value of the objective function are
!     available in X(.) on the return from GETMIN.
!  ACC is a tolerance on the first order conditions at the calculated
!     solution of the optimization problem.  These first order conditions
!     state that, if X(.) is a solution, then there is a set of active
!     constraints with indices IACT(K) K=1(1)NACT, say, such that X(.) is
!     on the boundaries of these constraints, and the gradient of the
!     objective function can be expressed in the form
!          GRAD(F)=PAR(1)*GRAD(C(IACT(1)))+...
!                        ...+PAR(NACT)*GRAD(C(IACT(NACT)))  .
!     Here PAR(K) K=1(1)NACT are Lagrange multipliers that are nonpositive
!     for inequality constraints, and GRAD(C(IACT(K))) is the gradient of
!     the IACT(K)-th constraint function, so it is A(.,IACT(K)) if IACT(K)
!     .LE. M, and it is minus or plus the J-th coordinate vector if the
!     constraint is the lower or upper bound on X(J) respectively.  The
!     normal return from the calculation occurs when X(.) is feasible and
!     the sum of squares of components of the vector RESKT(.) is at most
!     ACC**2, where RESKT(.) is the N-component vector of residuals of
!     the first order condition that is displayed above.  Sometimes the
!     package cannot satisfy this condition, because noise in the function
!     values can prevent a change to the variables, no line search being
!     allowed to increase the objective function.
!  IACT(.) is a working space array of integer variables that must be
!     provided by the user.  Its length must be at least (M+2*N).  Its
!     leading entries on the return from the subroutine are the indices
!     IACT(K) K=1(1)NACT that are mentioned in the previous paragraph:
!     in other words they are the indices of the final active constraints.
!     Here the indices M+1,...,M+N and M+N+1,...,M+2*N denote the lower
!     and upper bounds respectively.
!  NACT is set automatically to the integer variable of this ilk that has
!     been introduced already.
!  PAR(.) is a one-dimensional array that will hold the Lagrange
!     multipliers PAR(K) K=1(1)NACT on the return from GETMIN, these
!     parameters being defined in the above paragraph on ACC.  The length
!     of PAR should be at least N.
!  IPRINT must be set by the user to specify the frequency of printing
!     during the execution of the optimization package.  There is no
!     printed output if IPRINT=0.  Otherwise, after ensuring feasibility,
!     information is given every IABS(IPRINT) iterations and whenever a
!     parameter called TOL is reduced.  The printing provides the values
!     of X(.), F(.) and G(.)=GRAD(F) if IPRINT is positive, while if
!     IPRINT is negative this information is augmented by the current
!     values of IACT(K) K=1(1)NACT, PAR(K) K=1(1)NACT and RESKT(I)
!     I=1(1)N.  The reason for returning to the calling program is also
!     displayed when IPRINT is nonzero.
!  INFO is an integer variable that should be set to zero initially,
!     unless the user wishes to impose an upper bound on the number of
!     evaluations of the objective function and its gradient, in which
!     case INFO should be set to the value of this bound.  On the exit
!     from GETMIN it will have one of the following integer values to
!     indicate the reason for leaving the optimization package:
!          INFO=1   X(.) is feasible and the condition that depends on
!     ACC is satisfied.
!          INFO=2   X(.) is feasible and rounding errors are preventing
!     further progress.
!          INFO=3   X(.) is feasible but the objective function fails to
!     decrease although a decrease is predicted by the current gradient
!     vector.  If this return occurs and RESKT(.) has large components
!     then the user's calculation of the gradient of the objective
!     function may be incorrect.  One should also question the coding of
!     the gradient when the final rate of convergence is slow.
!          INFO=4   In this case the calculation cannot begin because IA
!     is less than N or because the lower bound on a variable is greater
!     than the upper bound.
!          INFO=5   This value indicates that the equality constraints
!     are inconsistent.   These constraints include any components of
!     X(.) that are frozen by setting XL(I)=XU(I).
!          INFO=6   In this case there is an error return because the
!     equality constraints and the bounds on the variables are found to
!     be inconsistent.
!          INFO=7   This value indicates that there is no vector of
!     variables that satisfies all of the constraints.  Specifically,
!     when this return or an INFO=6 return occurs, the current active
!     constraints (whose indices are IACT(K) K=1(1)NACT) prevent any
!     change in X(.) that reduces the sum of constraint violations,
!     where only bounds are included in this sum if INFO=6.
!          INFO=8   In this case the limit on the number of calls of
!     subroutine FGCALC (see below) has been reached, and there would
!     have been further calculation otherwise.
!  W(.) is a working space array of real variables that must be provided
!     by the user.  Its length must be at least (M+11*N+N**2).  On exit
!     from the package one can find the final components of GRAD(F) and
!     RESKT(.) in W(1),...,W(N) and W(N+1),...,W(2*N) respectively.
!  Note 1.   The variables N, M, MEQ, IA, ACC and IPRINT and the elements
!     of the arrays A(,.,), B(.), XL(.) and XU(.) are not altered by the
!     optimization procedure.  Their values, the value of INFO and the
!     initial components of X(.) must be set on entry to GETMIN.
!  Note 2.   Of course the package needs the objective function and its
!     gradient.  Therefore the user must provide a subroutine called
!     FGCALC, its first two lines being
!          SUBROUTINE FGCALC (N,X,F,G)
!          DIMENSION X(*),G(*)   .
!     It is called automatically with N set as above and with X(.) set
!     to a feasible vector of variables.  It should calculate the value
!     of the objective function and its gradient for this X(.) and should
!     set them in F and G(I) I=1(1)N respectively, without disturbing N
!     or any of the components of X(.).
!
!     Partition the workspace array.
!
      ig = 1
      ireskt = ig + n
      iz = ireskt + n
      iu = iz + n * n
      ixbig = iu + n
      ibres = ixbig + n
      id = ibres + m + n + n
      iztg = id + n
      igm = iztg + n
      ixs = igm + n
      igs = ixs + n
!
!     Call the optimization package.
!
      Call minflc (n, m, meq, a, ia, b, xl, xu, x, acc, iact, nact, &
       par, iprint, info, w(ig), w(iz), w(iu), w(ixbig), w(ireskt), &
       w(ibres), w(id), w(iztg), w(igm), w(ixs), w(igs), fgcalc)

End Subroutine getmin

Subroutine addcon (n, m, a, ia, iact, nact, z, u, relacc, indxbd, ztc, cgrad)
      Implicit real * 8 (a-h, o-z)
      Dimension a (ia,*), iact (*), z (*), u (*), ztc (*), cgrad (*)
      np = nact + 1
      icon = iact (indxbd)
      iact (indxbd) = iact (np)
      iact (np) = icon
!
!     Form ZTC when the new constraint is a bound.
!
      If (icon > m) Then
         inewbd = icon - m
         If (inewbd <= n) Then
            temp = - 1.0
         Else
            inewbd = inewbd - n
            temp = 1.0
         End If
         iznbd = inewbd * n - n
         Do 10 j = 1, n
10       ztc (j) = temp * z (iznbd+j)
!
!     Else form ZTC for an ordinary constraint.
!
      Else
         Do 20 i = 1, n
20       cgrad (i) = a (i, icon)
         Do 30 j = 1, n
            ztc (j) = 0.0
            iz = j
            Do 30 i = 1, n
               ztc (j) = ztc (j) + z (iz) * cgrad (i)
30       iz = iz + n
      End If
!
!     Find any Givens rotations to apply to the last columns of Z.
!
      j = n
40    jp = j
      j = j - 1
      If (j > nact) Then
         If (ztc(jp) == 0.0) Go To 40
         If (dabs(ztc(jp)) <= relacc*dabs(ztc(j))) Then
            temp = dabs (ztc(j))
         Else If (dabs(ztc(j)) <= relacc*dabs(ztc(jp))) Then
            temp = dabs (ztc(jp))
         Else
            temp = dabs (ztc(jp)) * dsqrt (1.0+(ztc(j)/ztc(jp))**2)
         End If
         wcos = ztc (j) / temp
         wsin = ztc (jp) / temp
         ztc (j) = temp
!
!     Apply the rotation when the new constraint is a bound.
!
         iz = j
         If (icon > m) Then
            Do 50 i = 1, n
               temp = wcos * z (iz+1) - wsin * z (iz)
               z (iz) = wcos * z (iz) + wsin * z (iz+1)
               z (iz+1) = temp
50          iz = iz + n
            z (iznbd+jp) = 0.0
!
!     Else apply the rotation for an ordinary constraint.
!
         Else
            wpiv = 0.0
            Do 60 i = 1, n
               tempa = wcos * z (iz+1)
               tempb = wsin * z (iz)
               temp = dabs (cgrad(i)) * (dabs(tempa)+dabs(tempb))
               If (temp > wpiv) Then
                  wpiv = temp
                  ipiv = i
               End If
               z (iz) = wcos * z (iz) + wsin * z (iz+1)
               z (iz+1) = tempa - tempb
60          iz = iz + n
!
!     Ensure orthogonality of Z(.,JP) to CGRAD.
!
            sum = 0.0
            iz = jp
            Do 70 i = 1, n
               sum = sum + z (iz) * cgrad (i)
70          iz = iz + n
            If (sum /= 0.0) Then
               iz = ipiv * n - n + jp
               z (iz) = z (iz) - sum / cgrad (ipiv)
            End If
         End If
         Go To 40
      End If
!
!     Test for linear independence in the proposed new active set.
!
      If (ztc(np) == 0.0) Go To 90
      If (icon <= m) Then
         sum = 0.0
         sumabs = 0.0
         iz = np
         Do 80 i = 1, n
            temp = z (iz) * cgrad (i)
            sum = sum + temp
            sumabs = sumabs + dabs (temp)
80       iz = iz + n
         If (dabs(sum) <= relacc*sumabs) Go To 90
      End If
!
!     Set the new diagonal element of U and return.
!
      u (np) = 1.0 / ztc (np)
      nact = np
90    Return
End Subroutine addcon

Subroutine adjtol (n, m, a, ia, b, xl, xu, x, iact, nact, xbig, relacc, &
& tol, meql)
      Implicit real * 8 (a-h, o-z)
      Dimension a (ia,*), b (*), xl (*), xu (*), x (*), iact (*), xbig &
     & (*)
!
!     Set VIOL to the greatest relative constraint residual of the first
!       NACT constraints.
!
      viol = 0.0
      If (nact > meql) Then
         kl = meql + 1
         Do 20 k = kl, nact
            j = iact (k)
            If (j <= m) Then
               res = b (j)
               resabs = dabs (b(j))
               Do 10 i = 1, n
                  res = res - a (i, j) * x (i)
10             resabs = resabs + dabs (a(i, j)*xbig(i))
            Else
               jm = j - m
               If (jm <= n) Then
                  res = x (jm) - xl (jm)
                  resabs = xbig (jm) + dabs (xl(jm))
               Else
                  jm = jm - n
                  res = xu (jm) - x (jm)
                  resabs = xbig (jm) + dabs (xu(jm))
               End If
            End If
            If (res > 0.0) viol = dmax1 (viol, res/resabs)
20       Continue
      End If
!
!     Adjust TOL.
!
      tol = 0.1 * dmin1 (tol, viol)
      If (tol <= relacc+relacc) Then
         tol = relacc
         Do 30 i = 1, n
30       xbig (i) = dabs (x(i))
      End If
      Return
End Subroutine adjtol

Subroutine conres (n, m, a, ia, b, xl, xu, x, iact, nact, par, g, z, u, &
& xbig, bres, d, ztg, relacc, tol, stepcb, sumres, meql, msat, mtot, &
& indxbd, gm, gmnew, parnew, cgrad)
      Implicit real * 8 (a-h, o-z)
      Dimension a (ia,*), b (*), xl (*), xu (*), x (*), iact (*), par &
     & (*), g (*), z (*), u (*), xbig (*), bres (*), d (*), ztg (*), gm &
     & (*), gmnew (*), parnew (*), cgrad (*)
      idiff = mtot - msat
!
!     Calculate and partition the residuals of the inactive constraints,
!       and set the gradient vector when seeking feasibility.
!
      If (idiff > 0.0) Then
         Do 10 i = 1, n
10       g (i) = 0.0
         sumres = 0.0
      End If
      msatk = msat
      mdeg = nact
      msat = nact
      kl = meql + 1
      Do 50 k = kl, mtot
         j = iact (k)
!
!     Calculate the residual of the current constraint.
!
         If (j <= m) Then
            res = b (j)
            resabs = dabs (b(j))
            Do 20 i = 1, n
               res = res - x (i) * a (i, j)
20          resabs = resabs + dabs (xbig(i)*a(i, j))
         Else
            jm = j - m
            If (jm <= n) Then
               res = x (jm) - xl (jm)
               resabs = dabs (xbig(jm)) + dabs (xl(jm))
            Else
               jm = jm - n
               res = xu (jm) - x (jm)
               resabs = dabs (xbig(jm)) + dabs (xu(jm))
            End If
         End If
         bres (j) = res
!
!     Set TEMP to the relative residual.
!
         temp = 0.0
         If (resabs /= 0.0) temp = res / resabs
         If (k > msatk .And. temp < 0.0) Then
            If (temp+relacc >= 0.0) Then
               If (j <= m) Then
                  sum = dabs (b(j))
                  Do 30 i = 1, n
30                sum = sum + dabs (x(i)*a(i, j))
               Else
                  jm = j - m
                  If (jm <= n) Then
                     sum = dabs (x(jm)) + dabs (xl(jm))
                  Else
                     sum = dabs (x(jm-n)) + dabs (xu(jm-n))
                  End If
               End If
               If (dabs(res) <= sum*relacc) temp = 0.0
            End If
         End If
!
!     Place the residual in the appropriate position.
!
         If (k <= nact) Go To 50
         If (k <= msatk .Or. temp >= 0.0) Then
            msat = msat + 1
            If (msat < k) Then
               iact (k) = iact (msat)
            End If
            If (temp > tol) Then
               iact (msat) = j
            Else
               mdeg = mdeg + 1
               iact (msat) = iact (mdeg)
               iact (mdeg) = j
            End If
!
!     Update the gradient and SUMRES if the constraint is violated when
!       seeking feasibility.
!
         Else
            If (j <= m) Then
               Do 40 i = 1, n
40             g (i) = g (i) + a (i, j)
            Else
               j = j - m
               If (j <= n) Then
                  g (j) = g (j) - 1.0
               Else
                  g (j-n) = g (j-n) + 1.0
               End If
            End If
            sumres = sumres + dabs (res)
         End If
50    Continue
!
!     Seek the next search direction unless CONRES was called from GETFES
!       and feasibility has been achieved.
!
      stepcb = 0.0
      If (idiff > 0 .And. msat == mtot) Go To 60
      Call getd (n, m, a, ia, iact, nact, par, g, z, u, d, ztg, relacc, &
     & ddotg, meql, mdeg, gm, gmnew, parnew, cgrad)
!
!     Calculate the (bound on the) step-length due to the constraints.
!
      If (ddotg < 0.0) Then
         Call stepbd (n, m, a, ia, iact, bres, d, stepcb, ddotg, mdeg, &
        & msat, mtot, indxbd)
      End If
      If (idiff == 0) sumres = ddotg
60    Return
End Subroutine conres

Subroutine delcon (n, m, a, ia, iact, nact, z, u, relacc, idrop)
      Implicit real * 8 (a-h, o-z)
      Dimension a (ia,*), iact (*), z (*), u (*)
      nm = nact - 1
      If (idrop == nact) Go To 60
      isave = iact (idrop)
!
!     Cycle through the constraint exchanges that are needed.
!
      Do 50 j = idrop, nm
         jp = j + 1
         icon = iact (jp)
         iact (j) = icon
!
!     Calculate the (J,JP) element of R.
!
         If (icon <= m) Then
            rjjp = 0.0
            iz = j
            Do 10 i = 1, n
               rjjp = rjjp + z (iz) * a (i, icon)
10          iz = iz + n
         Else
            ibd = icon - m
            If (ibd <= n) Then
               izbd = ibd * n - n
               rjjp = - z (izbd+j)
            Else
               ibd = ibd - n
               izbd = ibd * n - n
               rjjp = z (izbd+j)
            End If
         End If
!
!     Calculate the parameters of the next rotation.
!
         ujp = u (jp)
         temp = rjjp * ujp
         denom = dabs (temp)
         If (denom*relacc < 1.0) denom = dsqrt (1.0+denom*denom)
         wcos = temp / denom
         wsin = 1.0 / denom
!
!     Rotate Z when a bound constraint is promoted.
!
         iz = j
         If (icon > m) Then
            Do 20 i = 1, n
               temp = wcos * z (iz+1) - wsin * z (iz)
               z (iz) = wcos * z (iz) + wsin * z (iz+1)
               z (iz+1) = temp
20          iz = iz + n
            z (izbd+jp) = 0.0
!
!     Rotate Z when an ordinary constraint is promoted.
!
         Else
            wpiv = 0.0
            Do 30 i = 1, n
               tempa = wcos * z (iz+1)
               tempb = wsin * z (iz)
               temp = dabs (a(i, icon)) * (dabs(tempa)+dabs(tempb))
               If (temp > wpiv) Then
                  wpiv = temp
                  ipiv = i
               End If
               z (iz) = wcos * z (iz) + wsin * z (iz+1)
               z (iz+1) = tempa - tempb
30          iz = iz + n
!
!     Ensure orthogonality to promoted constraint.
!
            sum = 0.0
            iz = jp
            Do 40 i = 1, n
               sum = sum + z (iz) * a (i, icon)
40          iz = iz + n
            If (sum /= 0.0) Then
               iz = ipiv * n - n + jp
               z (iz) = z (iz) - sum / a (ipiv, icon)
            End If
         End If
!
!     Set the new diagonal elements of U.
!
         u (jp) = - denom * u (j)
         u (j) = ujp / denom
50    Continue
!
!     Return.
!
      iact (nact) = isave
60    nact = nm
      Return
End Subroutine delcon

Subroutine eqcons (n, m, meq, a, ia, b, xu, iact, meql, info, z, u, &
& relacc, am, cgrad)
      Implicit real * 8 (a-h, o-z)
      Dimension a (ia,*), b (*), xu (*), iact (*), z (*), u (*), am &
     & (*), cgrad (*)
!
!     Try to add the next equality constraint to the active set.
!
      Do 50 keq = 1, meq
         If (meql < n) Then
            np = meql + 1
            iact (np) = keq
            Call addcon (n, m, a, ia, iact, meql, z, u, relacc, np, am, &
           & cgrad)
            If (meql == np) Go To 50
         End If
!
!     If linear dependence occurs then find the multipliers of the
!       dependence relation and apply them to the right hand sides.
!
         sum = b (keq)
         sumabs = dabs (b(keq))
         If (meql > 0) Then
            Do 10 i = 1, n
10          am (i) = a (i, keq)
            k = meql
20          vmult = 0.0
            iz = k
            Do 30 i = 1, n
               vmult = vmult + z (iz) * am (i)
30          iz = iz + n
            vmult = vmult * u (k)
            j = iact (k)
            If (j <= m) Then
               Do 40 i = 1, n
40             am (i) = am (i) - vmult * a (i, j)
               rhs = b (j)
            Else
               jm = j - m - n
               am (jm) = am (jm) - vmult
               rhs = xu (jm)
            End If
            sum = sum - rhs * vmult
            sumabs = sumabs + dabs (rhs*vmult)
            k = k - 1
            If (k >= 1) Go To 20
         End If
!
!     Error return if the constraints are inconsistent.
!
         If (dabs(sum) > relacc*sumabs) Then
            info = 5
            Go To 60
         End If
50    Continue
60    Return
End Subroutine eqcons

!
Subroutine getd (n, m, a, ia, iact, nact, par, g, z, u, d, ztg, relacc, &
& ddotg, meql, mdeg, gm, gmnew, parnew, cgrad)
      Implicit real * 8 (a-h, o-z)
      Dimension a (ia,*), iact (*), par (*), g (*), z (*), u (*), d &
     & (*), ztg (*), gm (*), gmnew (*), parnew (*), cgrad (*)
!
!     Initialize GM and cycle backwards through the active set.
!
10    Do 20 i = 1, n
20    gm (i) = g (i)
      k = nact
30    If (k > 0) Then
!
!     Set TEMP to the next multiplier, but reduce the active set if
!       TEMP has an unacceptable sign.
!
         temp = 0.0
         iz = k
         Do 40 i = 1, n
            temp = temp + z (iz) * gm (i)
40       iz = iz + n
         temp = temp * u (k)
         If (k > meql .And. temp > 0.0) Then
            Call delcon (n, m, a, ia, iact, nact, z, u, relacc, k)
            Go To 10
         End If
!
!     Update GM using the multiplier that has just been calculated.
!
         j = iact (k)
         If (j <= m) Then
            Do 50 i = 1, n
50          gm (i) = gm (i) - temp * a (i, j)
         Else
            jm = j - m
            If (jm <= n) Then
               gm (jm) = gm (jm) + temp
            Else
               gm (jm-n) = gm (jm-n) - temp
            End If
         End If
         par (k) = temp
         k = k - 1
         Go To 30
      End If
!
!     Calculate the search direction and DDOTG.
!
      ddotg = 0.0
      If (nact < n) Then
         Call sdegen (n, m, a, ia, iact, nact, par, z, u, d, ztg, gm, &
        & relacc, ddotgm, meql, mdeg, gmnew, parnew, cgrad)
         If (ddotgm < 0.0) Then
            Do 60 i = 1, n
60          ddotg = ddotg + d (i) * g (i)
         End If
      End If
      Return
End Subroutine getd

Subroutine getfes (n, m, a, ia, b, xl, xu, x, iact, nact, par, info, g, &
& z, u, xbig, relacc, tol, meql, msat, mtot, bres, d, ztg, gm, gmnew, &
& parnew, cgrad)
      Implicit real * 8 (a-h, o-z)
      Dimension a (ia,*), b (*), xl (*), xu (*), x (*), iact (*), par &
     & (*), g (*), z (*), u (*), xbig (*), bres (*), d (*), ztg (*), gm &
     & (*), gmnew (*), parnew (*), cgrad (*)
!
!     Make the correction to X for the active constraints.
!
      info = 0
10    Call satact (n, m, a, ia, b, xl, xu, x, iact, nact, info, z, u, &
     & xbig, relacc, tol, meql)
      If (info > 0) msat = nact
      If (msat == mtot) Go To 60
!
!     Try to correct the infeasibility.
!
20    msatk = msat
      sumrsk = 0.0
30    Call conres (n, m, a, ia, b, xl, xu, x, iact, nact, par, g, z, u, &
     & xbig, bres, d, ztg, relacc, tol, stepcb, sumres, meql, msat, &
     & mtot, indxbd, gm, gmnew, parnew, cgrad)
!
!     Include the new constraint in the active set.
!
      If (stepcb > 0.0) Then
         Do 40 i = 1, n
            x (i) = x (i) + stepcb * d (i)
40       xbig (i) = dmax1 (xbig(i), dabs(x(i)))
         Call addcon (n, m, a, ia, iact, nact, z, u, relacc, indxbd, &
        & gmnew, cgrad)
      End If
!
!     Test whether to continue the search for feasibility.
!
      If (msat < mtot) Then
         If (stepcb == 0.0) Go To 50
         If (msatk < msat) Go To 20
         If (sumrsk == 0.0 .Or. sumres < sumrsk) Then
            sumrsk = sumres
            itest = 0
         End If
         itest = itest + 1
         If (itest <= 2) Go To 30
!
!     Reduce TOL if it may be too large to allow feasibility.
!
50       If (tol > relacc) Then
            Call adjtol (n, m, a, ia, b, xl, xu, x, iact, nact, xbig, &
           & relacc, tol, meql)
            Go To 10
         End If
      End If
60    Return
End Subroutine getfes

Subroutine initzu (n, m, xl, xu, x, iact, meql, info, z, u, xbig, &
& relacc)
      Implicit real * 8 (a-h, o-z)
      Dimension xl (*), xu (*), x (*), iact (*), z (*), u (*), xbig (*)
!
!     Set RELACC.
!
      ztpar = 100.0
      relacc = 1.0
10    relacc = 0.5 * relacc
      tempa = ztpar + 0.5 * relacc
      tempb = ztpar + relacc
      If (ztpar < tempa .And. tempa < tempb) Go To 10
!
!     Seek bound inconsistencies and bound equality constraints.
!
      meql = 0
      Do 20 i = 1, n
         If (xl(i) > xu(i)) Go To 50
         If (xl(i) == xu(i)) meql = meql + 1
20    Continue
!
!     Initialize U, Z and XBIG.
!
      jact = 0
      nn = n * n
      Do 30 i = 1, nn
30    z (i) = 0.0
      iz = 0
      Do 40 i = 1, n
         If (xl(i) == xu(i)) Then
            x (i) = xu (i)
            jact = jact + 1
            u (jact) = 1.0
            iact (jact) = i + m + n
            j = jact
         Else
            j = i + meql - jact
         End If
         z (iz+j) = 1.0
         iz = iz + n
40    xbig (i) = dabs (x(i))
      info = 1
50    Return
End Subroutine initzu

Subroutine ktvec (n, m, a, ia, iact, nact, par, g, reskt, z, u, bres, &
& relaxf, meql, ssqkt, parw, resktw)
      Implicit real * 8 (a-h, o-z)
      Dimension a (ia,*), iact (*), par (*), g (*), reskt (*), z (*), u &
     & (*), bres (*), parw (*), resktw (*)
!
!     Calculate the Lagrange parameters and the residual vector.
!
      Do 10 i = 1, n
10    reskt (i) = g (i)
      If (nact > 0) Then
         icase = 0
20       Do 50 kk = 1, nact
            k = nact + 1 - kk
            j = iact (k)
            temp = 0.0
            iz = k
            Do 30 i = 1, n
               temp = temp + z (iz) * reskt (i)
30          iz = iz + n
            temp = temp * u (k)
            If (icase == 0) par (k) = 0.0
            If (k <= meql .Or. par(k)+temp < 0.0) Then
               par (k) = par (k) + temp
            Else
               temp = - par (k)
               par (k) = 0.0
            End If
            If (temp /= 0.0) Then
               If (j <= m) Then
                  Do 40 i = 1, n
40                reskt (i) = reskt (i) - temp * a (i, j)
               Else
                  jm = j - m
                  If (jm <= n) Then
                     reskt (jm) = reskt (jm) + temp
                  Else
                     reskt (jm-n) = reskt (jm-n) - temp
                  End If
               End If
            End If
50       Continue
!
!     Calculate the sum of squares of the KT residual vector.
!
         ssqkt = 0.0
         If (nact == n) Go To 130
         Do 60 i = 1, n
60       ssqkt = ssqkt + reskt (i) ** 2
!
!     Apply iterative refinement to the residual vector.
!
         If (icase == 0) Then
            icase = 1
            Do 70 k = 1, nact
70          parw (k) = par (k)
            Do 80 i = 1, n
80          resktw (i) = reskt (i)
            ssqktw = ssqkt
            Go To 20
         End If
!
!     Undo the iterative refinement if it does not reduce SSQKT.
!
         If (ssqktw < ssqkt) Then
            Do 90 k = 1, nact
90          par (k) = parw (k)
            Do 100 i = 1, n
100         reskt (i) = resktw (i)
            ssqkt = ssqktw
         End If
!
!     Calculate SSQKT when there are no active constraints.
!
      Else
         ssqkt = 0.0
         Do 110 i = 1, n
110      ssqkt = ssqkt + g (i) ** 2
      End If
!
!     Predict the reduction in F if one corrects any positive residuals
!       of active inequality constraints.
!
      relaxf = 0.0
      If (meql < nact) Then
         kl = meql + 1
         Do 120 k = kl, nact
            j = iact (k)
            If (bres(j) > 0.0) Then
               relaxf = relaxf - par (k) * bres (j)
            End If
120      Continue
      End If
130   Return
End Subroutine ktvec

Subroutine lsrch (n, x, g, d, xs, gs, relacc, stepcb, ddotg, f, step, &
& nfvals, nfmax, gopt, fgcalc)
      Implicit real * 8 (a-h, o-z)
      Dimension x (*), g (*), d (*), xs (*), gs (*), gopt (*)
      procedure(func) :: fgcalc
!
!     Initialization.
!
      relint = 0.9
      icount = 0
      ratio = - 1.0
      Do 10 i = 1, n
         xs (i) = x (i)
         gs (i) = g (i)
         gopt (i) = g (i)
         If (d(i) /= 0.0) Then
            temp = dabs (x(i)/d(i))
            If (ratio < 0.0 .Or. temp < ratio) ratio = temp
         End If
10    Continue
      step = dmin1 (1.0d0, stepcb)
!
!     The following number 1.0D-12 is independent of the working
!       accuracy of the computer arithmetic.
!
      stpmin = dmax1 (relacc*ratio, 1.0d-12*step)
      step = dmax1 (stpmin, step)
      sbase = 0.0
      fbase = f
      ddotgb = ddotg
      stplow = 0.0
      flow = f
      dglow = ddotg
      stphgh = 0.0
      stpopt = 0.0
      fopt = f
      dgopt = dabs (ddotg)
!
!     Calculate another function and gradient value.
!
20    Do 30 i = 1, n
30    x (i) = xs (i) + step * d (i)
      Call fgcalc (n, x, f, g)
      icount = icount + 1
      dgmid = 0.0
      Do 40 i = 1, n
40    dgmid = dgmid + d (i) * g (i)
      If (f <= fopt) Then
         If (f < fopt .Or. dabs(dgmid) < dgopt) Then
            stpopt = step
            fopt = f
            Do 50 i = 1, n
50          gopt (i) = g (i)
            dgopt = dabs (dgmid)
         End If
      End If
      If (nfvals+icount == nfmax) Go To 70
!
!      Modify the bounds on the steplength or convergence.
!
      If (f >= fbase+0.1*(step-sbase)*ddotgb) Then
         If (stphgh > 0.0 .Or. f > fbase .Or. dgmid > &
        & 0.5*ddotg) Then
            stphgh = step
            fhgh = f
            dghgh = dgmid
            Go To 60
         End If
         sbase = step
         fbase = f
         ddotgb = dgmid
      End If
      If (dgmid >= 0.7*ddotgb) Go To 70
      stplow = step
      flow = f
      dglow = dgmid
60    If (stphgh > 0.0 .And. stplow >= relint*stphgh) Go To 70
!
!     Calculate the next step length or end the iterations.
!
      If (stphgh == 0.0) Then
         If (step == stepcb) Go To 70
         temp = 10.0
         If (dgmid > 0.9*ddotg) temp = ddotg / (ddotg-dgmid)
         step = dmin1 (temp*step, stepcb)
         Go To 20
      Else If (icount == 1 .Or. stplow > 0.0) Then
         dgknot = 2.0 * (fhgh-flow) / (stphgh-stplow) - 0.5 * &
        & (dglow+dghgh)
         If (dgknot >= 0.0) Then
            ratio = dmax1 (0.1d0, 0.5d0*dglow/(dglow-dgknot))
         Else
            ratio = (0.5*dghgh-dgknot) / (dghgh-dgknot)
         End If
         step = stplow + ratio * (stphgh-stplow)
         Go To 20
      Else
         step = 0.1 * step
         If (step >= stpmin) Go To 20
      End If
!
!     Return from subroutine.
!
70    If (step /= stpopt) Then
         step = stpopt
         f = fopt
         Do 80 i = 1, n
            x (i) = xs (i) + step * d (i)
80       g (i) = gopt (i)
      End If
      nfvals = nfvals + icount
      Return
End Subroutine lsrch

Subroutine minflc (n, m, meq, a, ia, b, xl, xu, x, acc, iact, nact, &
  par, iprint, info, g, z, u, xbig, reskt, bres, d, ztg, gm, xs, gs, fgcalc)
   Implicit real * 8 (a-h, o-z)
   Dimension a (ia,*), b (*), xl (*), xu (*), x (*), iact (*), par (*), &
    g (*), z (*), u (*), xbig (*), reskt (*), bres (*), d (*), ztg (*), &
    gm (*), xs (*), gs (*)
      procedure(func) :: fgcalc
!
!     Initialize ZZNORM, ITERC, NFVALS and NFMAX.
!
   zznorm = - 1.0
   iterc = 0
   nfvals = 0
   nfmax = 0
   If (info > 0) nfmax = info
!
!     Check the bounds on N, M and MEQ.
!
   info = 4
   If (max0(1-n,-m, meq*(meq-m)) > 0) Then
      If (iprint /= 0) Print 1010
1010  Format (/ 5 x, 'ERROR RETURN FROM GETMIN BECAUSE A CONDITION',&
      ' ON N, M OR MEQ IS VIOLATED')
      Go To 40
   End If
!
!     Initialize RELACC, Z, U and TOL.
!
   Call initzu (n, m, xl, xu, x, iact, meql, info, z, u, xbig, relacc)
   tol = dmax1 (0.01d0, 10.0d0*relacc)
   If (info == 4) Then
      If (iprint /= 0) Print 1020
1020  Format (/ 5 x, 'ERROR RETURN FROM GETMIN BECAUSE A LOWER',&
      ' BOUND EXCEEDS AN UPPER BOUND')
      Go To 40
   End If
!
!     Add any equality constraints to the active set.
!
   If (meq > 0) Then
      Call eqcons (n, m, meq, a, ia, b, xu, iact, meql, info, z, u, &
     & relacc, xs, gs)
      If (info == 5) Then
         If (iprint /= 0) Print 1030
1030     Format (/ 5 x, 'ERROR RETURN FROM GETMIN BECAUSE THE',&
         ' EQUALITY CONSTRAINTS ARE INCONSISTENT')
         Go To 40
      End If
   End If
   nact = meql
   msat = meql
!
!     Add the bounds to the list of constraints.
!
   mtot = nact
   Do 10 i = 1, n
      If (xl(i) < xu(i)) Then
         mtot = mtot + 2
         iact (mtot-1) = m + i
         iact (mtot) = m + n + i
      End If
10 Continue
!
!     Try to satisfy the bound constraints.
!
   Call getfes (n, m, a, ia, b, xl, xu, x, iact, nact, par, info, g, z, &
  & u, xbig, relacc, tol, meql, msat, mtot, bres, d, ztg, gm, reskt, &
  & xs, gs)
   If (msat < mtot) Then
      If (iprint /= 0) Print 1040
1040  Format (/ 5 x, 'ERROR RETURN FROM GETMIN BECAUSE THE',&
      ' EQUALITIES AND BOUNDS ARE INCONSISTENT')
      info = 6
      Go To 40
   End If
!
!     Add the ordinary inequalities to the list of constraints.
!
   If (m > meq) Then
      mp = meq + 1
      Do 20 k = mp, m
         mtot = mtot + 1
20    iact (mtot) = k
   End If
!
!     Correct any constraint violations.
!
30 Call getfes (n, m, a, ia, b, xl, xu, x, iact, nact, par, info, g, z, &
  & u, xbig, relacc, tol, meql, msat, mtot, bres, d, ztg, gm, reskt, &
  & xs, gs)
   If (msat < mtot) Then
      If (iprint /= 0) Print 1050
1050  Format (/ 5 x, 'ERROR RETURN FROM GETMIN BECAUSE THE',&
      ' CONSTRAINTS ARE INCONSISTENT')
      info = 7
      Go To 40
   Else If (meql == n) Then
      If (iprint /= 0) Print 1060
1060  Format (/ 5 x, 'GETMIN FINDS THAT THE VARIABLES ARE',&
      ' DETERMINED BY THE EQUALITY CONSTRAINTS')
      Go To 40
   End If
!
!     Minimize the objective function in the case when constraints are
!       treated as degenerate if their residuals are less than TOL.
!
   Call minfun (n, m, a, ia, b, xl, xu, x, acc, iact, nact, par, &
  & iprint, info, g, z, u, xbig, relacc, zznorm, tol, meql, mtot, &
  & iterc, nfvals, nfmax, reskt, bres, d, ztg, gm, xs, gs, fgcalc)
!
!     Reduce TOL if necessary.
!
   If (tol > relacc .And. nact > 0) Then
      If (nfvals /= nfmax) Then
         Call adjtol (n, m, a, ia, b, xl, xu, x, iact, nact, xbig, &
        & relacc, tol, meql)
         Go To 30
      Else
         info = 8
      End If
   End If
   If (iprint /= 0) Then
      If (info == 1) Print 1070
1070  Format (/ 5 x, 'GETMIN HAS ACHIEVED THE REQUIRED ACCURACY')
      If (info == 2) Print 1080
1080  Format (/ 5 x, 'GETMIN CAN MAKE NO FURTHER PROGRESS BECAUSE',&
      ' OF ROUNDING ERRORS')
      If (info == 3) Print 1090
1090  Format (/ 5 x, 'GETMIN CAN MAKE NO FURTHER PROGRESS BECAUSE',&
      ' F WILL NOT DECREASE ANY MORE')
      If (info == 8) Print 1100
1100  Format (/ 5 x, 'GETMIN HAS REACHED THE GIVEN LIMIT ON THE',&
      ' NUMBER OF CALLS OF FGCALC')
   End If
40 Return
End Subroutine minflc

Subroutine minfun (n, m, a, ia, b, xl, xu, x, acc, iact, nact, par, &
& iprint, info, g, z, u, xbig, relacc, zznorm, tol, meql, mtot, iterc, &
& nfvals, nfmax, reskt, bres, d, ztg, gm, xs, gs, fgcalc)
   Implicit real * 8 (a-h, o-z)
   Dimension a (ia,*), b (*), xl (*), xu (*), x (*), iact (*), par (*), &
  & g (*), z (*), u (*), xbig (*), reskt (*), bres (*), d (*), ztg (*), &
  & gm (*), xs (*), gs (*)
   Save f
   procedure(func) :: fgcalc
!
!     Initialize the minimization calculation.
!
   msat = mtot
   iterk = iterc
   nfvalk = nfvals
   If (nfvals == 0 .Or. info == 1) Then
      Call fgcalc (n, x, f, g)
      nfvals = nfvals + 1
   End If
   fprev = dabs (f+f+1.0)
   iterp = - 1
   If (iprint /= 0) Then
      Print 1000, tol
1000  Format (/ 5 x, 'NEW VALUE OF TOL =', 1 pd13.5)
      iterp = iterc + iabs (iprint)
      If (iterc == 0) iterp = 0
   End If
!
!     Calculate the next search direction.
!
10 Call conres (n, m, a, ia, b, xl, xu, x, iact, nact, par, g, z, u, &
  & xbig, bres, d, ztg, relacc, tol, stepcb, ddotg, meql, msat, mtot, &
  & indxbd, gm, reskt, xs, gs)
!
!     Calculate the Kuhn Tucker residual vector.
!
   Call ktvec (n, m, a, ia, iact, nact, par, g, reskt, z, u, bres, &
  & relaxf, meql, ssqkt, xs, gs)
!
!     Test for convergence.
!
   If (ssqkt <= acc*acc) Then
      info = 1
      Go To 70
   End If
   If (ddotg >= 0.0) Then
      info = 2
      Go To 70
   End If
!
!     Test for termination due to no decrease in F.
!
   If (f >= fprev) Then
      If (tol == relacc .Or. nact == 0) Then
         If (diff > 0.0) Go To 20
      End If
      info = 3
      Go To 70
   End If
20 diff = fprev - f
   fprev = f
!
!     Test that more calls of FGCALC are allowed.
!
   If (nfvals == nfmax) Then
      info = 8
      Go To 70
   End If
!
!     Test whether to reduce TOL and to provide printing.
!
   If (tol > relacc .And. iterc > iterk .And. 0.1*relaxf >= &
  & dmax1(diff,-0.5d0*ddotg)) Go To 70
   If (iterp == iterc) Go To 80
!
!     Calculate the step along the search direction.
!
40 iterc = iterc + 1
   Call lsrch (n, x, g, d, xs, gs, relacc, stepcb, ddotg, f, step, &
    nfvals, nfmax, bres, fgcalc)
   If (step == 0.0) Then
      info = 3
      sum = 0.0
      Do 50 i = 1, n
50    sum = sum + dabs (d(i)*gs(i))
      If (ddotg+relacc*sum >= 0.0) info = 2
      Go To 70
   End If
!
!     Revise XBIG.
!
   Do 60 i = 1, n
60 xbig (i) = dmax1 (xbig(i), dabs(x(i)))
!
!     Revise the second derivative approximation.
!
   Call zbfgs (n, x, nact, g, z, ztg, xs, gs, zznorm)
!
!     Add a constraint to the active set if it restricts the step.
!
   If (step == stepcb) Then
      k = iact (indxbd)
      If (k > m) Then
         k = k - m
         If (k <= n) Then
            x (k) = xl (k)
         Else
            x (k-n) = xu (k-n)
         End If
      End If
      Call addcon (n, m, a, ia, iact, nact, z, u, relacc, indxbd, xs, &
     & gs)
   End If
   Go To 10
!
!     Printing from the subroutine.
!
70 If (iprint == 0) Go To 90
   iterk = - 1
80 Print 1010, iterc, nfvals, f
1010 Format (/ 5 x, 'ITERS =', i4, 5 x, 'F.VALS =', i4, 5 x, 'F =', 1 &
  & pd15.7)
   Print 1020, (x(i), i=1, n)
1020 Format ('  X =', (1 p5d14.5))
   Print 1030, (g(i), i=1, n)
1030 Format ('  G =', (1 p5d14.5))
   If (iprint < 0) Then
      If (nact == 0) Then
         Print 1050
1050     Format (5 x, 'NO ACTIVE CONSTRAINTS')
      Else
         Print 1060, (iact(i), i=1, nact)
1060     Format (' IA =', (14 i5))
         Print 1070, (par(i), i=1, nact)
1070     Format (' LP =', (1 p5d14.5))
      End If
      If (nact == n) Then
         Print 1080
1080     Format (5 x, 'KT RESIDUAL VECTOR IS ZERO')
      Else
         Print 1090, (reskt(i), i=1, n)
1090     Format (' KT =', (1 p5d14.5))
      End If
   End If
   iterp = iterc + iabs (iprint)
   If (iterk >= 0) Go To 40
90 Return
End Subroutine minfun

Subroutine newcon (n, m, a, ia, iact, nact, z, u, d, relacc, mdeg, &
& zzdiag, gmnew, cgrad)
   Implicit real * 8 (a-h, o-z)
   Dimension a (ia,*), iact (*), z (*), u (*), d (*), zzdiag (*), gmnew &
  & (*), cgrad (*)
!
!     Initialization.
!
   np = nact + 1
   khigh = mdeg
   iz = 0
   Do 20 i = 1, n
      zzdiag (i) = 0.0
      Do 10 j = np, n
10    zzdiag (i) = zzdiag (i) + z (iz+j) ** 2
20 iz = iz + n
!
!     Calculate the scalar products of D with its constraints.
!
30 cvmax = 0.0
   Do 50 k = np, khigh
      j = iact (k)
      If (j <= m) Then
         sum = 0.0
         sumabs = 0.0
         sumd = 0.0
         Do 40 i = 1, n
            temp = d (i) * a (i, j)
            sum = sum + temp
            sumabs = sumabs + dabs (temp)
40       sumd = sumd + zzdiag (i) * a (i, j) ** 2
      Else
         jm = j - m
         If (jm <= n) Then
            sum = - d (jm)
         Else
            jm = jm - n
            sum = d (jm)
         End If
         sumabs = dabs (sum)
         sumd = zzdiag (jm)
      End If
!
!     Pick out the most violated constraint, or return if the
!       violation is negligible.
!
      If (sum > relacc*sumabs) Then
         cviol = sum * sum / sumd
         If (cviol > cvmax) Then
            cvmax = cviol
            iadd = k
            savsum = sum
            savabs = sumabs
         End If
      End If
50 Continue
   If (cvmax <= 0.0) Go To 140
   If (nact == 0) Go To 120
!
!     Set GMNEW to the gradient of the most violated constraint.
!
   j = iact (iadd)
   If (j <= m) Then
      jmv = 0
      Do 60 i = 1, n
60    gmnew (i) = a (i, j)
   Else
      jmv = j - m
      Do 70 i = 1, n
70    gmnew (i) = 0.0
      If (jmv <= n) Then
         gmnew (jmv) = - 1.0
      Else
         jmv = jmv - n
         gmnew (jmv) = 1.0
      End If
   End If
!
!     Modify GMNEW for the next active constraint.
!
   k = nact
80 temp = 0.0
   iz = k
   Do 90 i = 1, n
      temp = temp + z (iz) * gmnew (i)
90 iz = iz + n
   temp = temp * u (k)
   j = iact (k)
   If (j <= m) Then
      Do 100 i = 1, n
100   gmnew (i) = gmnew (i) - temp * a (i, j)
   Else
      jm = j - m
      If (jm <= n) Then
         gmnew (jm) = gmnew (jm) + temp
      Else
         gmnew (jm-n) = gmnew (jm-n) - temp
      End If
   End If
!
!     Revise the values of SAVSUM and SAVABS.
!
   sum = 0.0
   sumabs = 0.0
   Do 110 i = 1, n
      temp = d (i) * gmnew (i)
      sum = sum + temp
110 sumabs = sumabs + dabs (temp)
   savsum = dmin1 (savsum, sum)
   savabs = dmax1 (savabs, sumabs)
   k = k - 1
   If (k >= 1) Go To 80
!
!     Add the new constraint to the active set if the constraint
!       violation is still significant.
!
   If (jmv > 0) d (jmv) = 0.0
   If (savsum <= relacc*savabs) Go To 130
120 k = nact
   Call addcon (n, m, a, ia, iact, nact, z, u, relacc, iadd, gmnew, &
  & cgrad)
   If (nact > k) Go To 140
!
!     Seek another constraint violation.
!
   iadd = np
130 If (np < khigh) Then
      k = iact (khigh)
      iact (khigh) = iact (iadd)
      iact (iadd) = k
      khigh = khigh - 1
      Go To 30
   End If
140 Return
End Subroutine newcon

Subroutine satact (n, m, a, ia, b, xl, xu, x, iact, nact, info, z, u, &
  xbig, relacc, tol, meql)
   Implicit real * 8 (a-h, o-z)
   Dimension a (ia,*), b (*), xl (*), xu (*), x (*), iact (*), z (*), u(*), xbig (*)
   If (nact == 0) Go To 50
   Do 30 k = 1, nact
!
!     Calculate the next constraint residual.
!
      j = iact (k)
      If (j <= m) Then
         res = b (j)
         resabs = dabs (b(j))
         resbig = resabs
         Do 10 i = 1, n
            tempa = a (i, j)
            temp = tempa * x (i)
            res = res - temp
            resabs = resabs + dabs (temp)
10       resbig = resbig + dabs (tempa) * xbig (i)
      Else
         jx = j - m
         If (jx <= n) Then
            res = x (jx) - xl (jx)
            resabs = dabs (x(jx)) + dabs (xl(jx))
            resbig = xbig (jx) + dabs (xl(jx))
            savex = xl (jx)
         Else
            jx = jx - n
            res = xu (jx) - x (jx)
            resabs = dabs (x(jx)) + dabs (xu(jx))
            resbig = xbig (jx) + dabs (xu(jx))
            savex = xu (jx)
         End If
      End If
!
!     Shift X if necessary.
!
      If (res /= 0.0) Then
         temp = res / resabs
         If (k <= meql) temp = - dabs (temp)
         If (tol == relacc .Or. temp+relacc < 0.0) Then
            info = 1
            scale = res * u (k)
            iz = k
            Do 20 i = 1, n
               x (i) = x (i) + scale * z (iz)
               iz = iz + n
20          xbig (i) = dmax1 (xbig(i), dabs(x(i)))
            If (j > m) x (jx) = savex
!
!     Else flag a constraint deletion if necessary.
!
         Else If (res/resbig > tol) Then
            iact (k) = - iact (k)
         End If
      End If
30 Continue
!
!     Delete any flagged constraints and then return.
!
   idrop = nact
40 If (iact(idrop) < 0) Then
      iact (idrop) = - iact (idrop)
      Call delcon (n, m, a, ia, iact, nact, z, u, relacc, idrop)
   End If
   idrop = idrop - 1
   If (idrop > meql) Go To 40
50 Return
End Subroutine satact

Subroutine sdegen (n, m, a, ia, iact, nact, par, z, u, d, ztg, gm, &
  relacc, ddotgm, meql, mdeg, gmnew, parnew, cgrad)
   Implicit real * 8 (a-h, o-z)
   Dimension a (ia,*), iact (*), par (*), z (*), u (*), d (*), ztg (*), &
   gm (*), gmnew (*), parnew (*), cgrad (*)
   mp = meql + 1
   dtest = 0.0
!
!     Calculate the search direction and branch if it is not downhill.
!
10 Call sdirn (n, nact, z, d, ztg, gm, relacc, ddotgm)
   If (ddotgm == 0.0) Go To 120
!
!     Branch if there is no need to consider any degenerate constraints.
!     The test gives termination if two consecutive additions to the
!       active set fail to increase the predicted new value of F.
!
   If (nact == mdeg) Go To 120
   np = nact + 1
   sum = 0.0
   Do 20 j = np, n
20 sum = sum + ztg (j) ** 2
   If (dtest > 0.0 .And. sum >= dtest) Then
      If (itest == 1) Go To 120
      itest = 1
   Else
      dtest = sum
      itest = 0
   End If
!
!     Add a constraint to the active set if there are any significant
!       violations of degenerate constraints.
!
   k = nact
   Call newcon (n, m, a, ia, iact, nact, z, u, d, relacc, mdeg, gmnew, &
  & parnew, cgrad)
   If (nact == k) Go To 120
   par (nact) = 0.0
!
!     Calculate the new reduced gradient and Lagrange parameters.
!
30 Do 40 i = 1, n
40 gmnew (i) = gm (i)
   k = nact
50 temp = 0.0
   iz = k
   Do 60 i = 1, n
      temp = temp + z (iz) * gmnew (i)
60 iz = iz + n
   temp = temp * u (k)
   parnew (k) = par (k) + temp
   If (k == nact) parnew (k) = dmin1 (parnew(k), 0.0d0)
   j = iact (k)
   If (j <= m) Then
      Do 70 i = 1, n
70    gmnew (i) = gmnew (i) - temp * a (i, j)
   Else
      jm = j - m
      If (jm <= n) Then
         gmnew (jm) = gmnew (jm) + temp
      Else
         gmnew (jm-n) = gmnew (jm-n) - temp
      End If
   End If
   k = k - 1
   If (k > meql) Go To 50
!
!     Set RATIO for linear interpolation between PAR and PARNEW.
!
   ratio = 0.0
   If (mp < nact) Then
      ku = nact - 1
      Do 80 k = mp, ku
         If (parnew(k) > 0.0) Then
            ratio = parnew (k) / (parnew(k)-par(k))
            idrop = k
         End If
80    Continue
   End If
!
!     Apply the linear interpolation.
!
   theta = 1.0 - ratio
   Do 90 k = mp, nact
90 par (k) = dmin1 (theta*parnew(k)+ratio*par(k), 0.0d0)
   Do 100 i = 1, n
100 gm (i) = theta * gmnew (i) + ratio * gm (i)
!
!     Drop a constraint if RATIO is positive.
!
   If (ratio > 0.0) Then
      Call delcon (n, m, a, ia, iact, nact, z, u, relacc, idrop)
      Do 110 k = idrop, nact
110   par (k) = par (k+1)
      Go To 30
   End If
!
!     Return if there is no freedom for a new search direction.
!
   If (nact < n) Go To 10
   ddotgm = 0.0
120 Return
End Subroutine sdegen

Subroutine sdirn (n, nact, z, d, ztg, gm, relacc, ddotgm)
   Implicit real * 8 (a-h, o-z)
   Dimension z (*), d (*), ztg (*), gm (*)
   ddotgm = 0.0
   If (nact >= n) Go To 60
!
!     Premultiply GM by the transpose of Z.
!
   np = nact + 1
   Do 20 j = np, n
      sum = 0.0
      sumabs = 0.0
      iz = j
      Do 10 i = 1, n
         temp = z (iz) * gm (i)
         sum = sum + temp
         sumabs = sumabs + dabs (temp)
10    iz = iz + n
      If (dabs(sum) <= relacc*sumabs) sum = 0.0
20 ztg (j) = sum
!
!     Form D by premultiplying ZTG by -Z.
!
   iz = 0
   Do 40 i = 1, n
      sum = 0.0
      sumabs = 0.0
      Do 30 j = np, n
         temp = z (iz+j) * ztg (j)
         sum = sum - temp
30    sumabs = sumabs + dabs (temp)
      If (dabs(sum) <= relacc*sumabs) sum = 0.0
      d (i) = sum
40 iz = iz + n
!
!     Test that the search direction is downhill.
!
   sumabs = 0.0
   Do 50 i = 1, n
      temp = d (i) * gm (i)
      ddotgm = ddotgm + temp
50 sumabs = sumabs + dabs (temp)
   If (ddotgm+relacc*sumabs >= 0.0) ddotgm = 0.0
60 Return
End Subroutine sdirn

Subroutine stepbd (n, m, a, ia, iact, bres, d, stepcb, ddotg, mdeg, &
  msat, mtot, indxbd)
   Implicit real * 8 (a-h, o-z)
   Dimension a (ia,*), iact (*), bres (*), d (*)
!
!     Set steps to constraint boundaries and find the least positive one.
!
   iflag = 0
   stepcb = 0.0
   indxbd = 0
   k = mdeg
10 k = k + 1
   If (k > mtot) Go To 40
!
!     Form the scalar product of D with the current constraint normal.
!
20 j = iact (k)
   If (j <= m) Then
      sp = 0.0
      Do 30 i = 1, n
30    sp = sp + d (i) * a (i, j)
   Else
      jm = j - m
      If (jm <= n) Then
         sp = - d (jm)
      Else
         sp = d (jm-n)
      End If
   End If
!
!     The next branch is taken if label 20 was reached via label 50.
!
   If (iflag == 1) Go To 60
!
!     Set BRES(J) to indicate the status of the j-th constraint.
!
   If (sp*bres(j) <= 0.0) Then
      bres (j) = 0.0
   Else
      bres (j) = bres (j) / sp
      If (stepcb == 0.0 .Or. bres(j) < stepcb) Then
         stepcb = bres (j)
         indxbd = k
      End If
   End If
   Go To 10
40 Continue
!
!     Try to pass through the boundary of a violated constraint.
!
50 If (indxbd <= msat) Go To 80
   iflag = 1
   k = indxbd
   Go To 20
60 msat = msat + 1
   iact (indxbd) = iact (msat)
   iact (msat) = j
   bres (j) = 0.0
   indxbd = msat
   ddotg = ddotg - sp
   If (ddotg < 0.0 .And. msat < mtot) Then
!
!     Seek the next constraint boundary along the search direction.
!
      temp = 0.0
      kl = mdeg + 1
      Do 70 k = kl, mtot
         j = iact (k)
         If (bres(j) > 0.0) Then
            If (temp == 0.0 .Or. bres(j) < temp) Then
               temp = bres (j)
               indxbd = k
            End If
         End If
70    Continue
      If (temp > 0.0) Then
         stepcb = temp
         Go To 50
      End If
   End If
80 Continue
   Return
End Subroutine stepbd

Subroutine zbfgs (n, x, nact, g, z, ztg, xs, gs, zznorm)
   Implicit real * 8 (a-h, o-z)
   Dimension x (*), g (*), z (*), ztg (*), xs (*), gs (*)
!
!     Test if there is sufficient convexity for the update.
!
   dd = 0.0
   dg = 0.0
   temp = 0.0
   Do 10 i = 1, n
      xs (i) = x (i) - xs (i)
      dd = dd + xs (i) ** 2
      temp = temp + gs (i) * xs (i)
      gs (i) = g (i) - gs (i)
10 dg = dg + gs (i) * xs (i)
   If (dg < 0.1*dabs(temp)) Go To 90
!
!     Transform the Z matrix.
!
   k = n
20 kp = k
   k = k - 1
   If (k > nact) Then
      If (ztg(kp) == 0.0) Go To 20
      temp = dabs (ztg(kp)) * dsqrt (1.0+(ztg(k)/ztg(kp))**2)
      wcos = ztg (k) / temp
      wsin = ztg (kp) / temp
      ztg (k) = temp
      iz = k
      Do 30 i = 1, n
         temp = wcos * z (iz+1) - wsin * z (iz)
         z (iz) = wcos * z (iz) + wsin * z (iz+1)
         z (iz+1) = temp
30    iz = iz + n
      Go To 20
   End If
!
!     Update the value of ZZNORM.
!
   If (zznorm < 0.0) Then
      zznorm = dd / dg
   Else
      temp = dsqrt (zznorm*dd/dg)
      zznorm = dmin1 (zznorm, temp)
      zznorm = dmax1 (zznorm, 0.1d0*temp)
   End If
!
!     Complete the updating of Z.
!
   np = nact + 1
   temp = dsqrt (dg)
   iz = np
   Do 40 i = 1, n
      z (iz) = xs (i) / temp
40 iz = iz + n
   If (np < n) Then
      km = np + 1
      Do 80 k = km, n
         temp = 0.0
         iz = k
         Do 50 i = 1, n
            temp = temp + gs (i) * z (iz)
50       iz = iz + n
         temp = temp / dg
         sum = 0.0
         iz = k
         Do 60 i = 1, n
            z (iz) = z (iz) - temp * xs (i)
            sum = sum + z (iz) ** 2
60       iz = iz + n
         If (sum < zznorm) Then
            temp = dsqrt (zznorm/sum)
            iz = k
            Do 70 i = 1, n
               z (iz) = temp * z (iz)
70          iz = iz + n
         End If
80    Continue
   End If
90 Return
End Subroutine zbfgs

subroutine tolmin_test()
!
!     The pentagon problem.
!
Implicit real * 8 (a-h, o-z)
Dimension a (10, 15), b (15), xl (6), xu (6), x (6), iact (27), par &
& (20), w (1000)
!
!     The two values of ICASE provide two different values of ACC, the latter
!     accuracy being so demanding that a return with INFO=2 occurs. The
!     final values of the objective function in the two cases agree well
!     and constraint violations are negligible, considering the differences
!     between the final values of the variables.
!
iprint = 10
ia = 10
n = 6
Do 100 icase = 1, 2
   acc = 1.0d-6
   If (icase == 2) acc = 1.0d-14
!
!     Set the components of XL, XU and X.
!
   Do 10 i = 1, n
      xl (i) = - 1.0d6
      xu (i) = 1.0d6
10 x (i) = 0.5d0 * dfloat (i-3)
   x (2) = 0.0d0
   x (4) = - 1.0d0
   x (6) = 1.0d0
!
!     Set the constraints.
!
   m = 0
   meq = 0
   pi = 4.0d0 * datan (1.0d0)
   Do 30 k = 1, 5
      theta = 0.4d0 * dfloat (k-1) * pi
      cth = dcos (theta)
      sth = dsin (theta)
      Do 30 j = 2, n, 2
         m = m + 1
         Do 20 i = 1, n
20       a (i, m) = 0.0d0
         a (j-1, m) = cth
         a (j, m) = sth
30 b (m) = 1.0d0
!
!     Call the optimization package.
!
   info = 0
   Print 40, acc, iprint
40 Format (/ / 5 x, 'CALL OF GETMIN WITH  ACC =', 1 pd11.4, '  AND  IPR&
  &INT =', i3)
   Call getmin (n, m, meq, a, ia, b, xl, xu, x, acc, iact, nact, par, &
  & iprint, info, w, fgcalc)
   Print 50, info
50 Format (/ 5 x, 'RETURN FROM TOLMIN WITH INFO =', i2)
   Call fgcalc (n, x, f, w)
   Print 60, f
60 Format (/ 5 x, 'FINAL VALUE OF OBJECTIVE FUNCTION =', 1 pd20.12)
   Print 70, (x(i), i=1, n)
70 Format (/ 5 x, 'FINAL COMPONENTS OF X =' // (4 x, 1 p3d20.12))
   Do 80 k = 1, m
      Do 80 i = 1, n
80 b (k) = b (k) - a (i, k) * x (i)
   Print 90, (b(k), k=1, m)
90 Format (/ 5 x, 'FINAL CONSTRAINT RESIDUALS =' // (3 x, 1 p6d12.4))
100 Continue

contains

    Subroutine fgcalc (n, x, f, g)
          Implicit real * 8 (a-h, o-z)
          Dimension x (*), g (*)
    !
    !     Calculate the objective function and its gradient.
    !
          wa = (x(1)-x(3)) ** 2 + (x(2)-x(4)) ** 2
          wb = (x(3)-x(5)) ** 2 + (x(4)-x(6)) ** 2
          wc = (x(5)-x(1)) ** 2 + (x(6)-x(2)) ** 2
          f = 1.0 / (wa**8) + 1.0 / (wb**8) + 1.0 / (wc**8)
          g (1) = 16.0 * ((x(3)-x(1))/(wa**9)+(x(5)-x(1))/(wc**9))
          g (2) = 16.0 * ((x(4)-x(2))/(wa**9)+(x(6)-x(2))/(wc**9))
          g (3) = 16.0 * ((x(5)-x(3))/(wb**9)+(x(1)-x(3))/(wa**9))
          g (4) = 16.0 * ((x(6)-x(4))/(wb**9)+(x(2)-x(4))/(wa**9))
          g (5) = 16.0 * ((x(1)-x(5))/(wc**9)+(x(3)-x(5))/(wb**9))
          g (6) = 16.0 * ((x(2)-x(6))/(wc**9)+(x(4)-x(6))/(wb**9))
    End Subroutine fgcalc

End subroutine tolmin_test

end module tolmin_module