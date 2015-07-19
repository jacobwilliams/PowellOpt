    module lincoa_module
    
    use kind_module, only: wp

    private
 
    abstract interface
        subroutine func (N,X,F)  !! calfun interface
        import :: wp
        implicit none
        integer :: n
        real(wp) :: x(*)
        real(wp) :: f
        end subroutine func
    end interface
   
    public :: lincoa
    public :: lincoa_test
    
    contains

Subroutine lincoa (n, npt, m, a, ia, b, x, rhobeg, rhoend, iprint, &
  maxfun, w, calfun)
      Implicit real(wp) (a-h, o-z)
      Dimension a (ia,*), b (*), x (*), w (*)
      procedure(func) :: calfun
!
!     This subroutine seeks the least value of a function of many variables,
!       subject to general linear inequality constraints, by a trust region
!       method that forms quadratic models by interpolation. Usually there
!       is much freedom in each new model after satisfying the interpolation
!       conditions, which is taken up by minimizing the Frobenius norm of
!       the change to the second derivative matrix of the model. One new
!       function value is calculated on each iteration, usually at a point
!       where the current model predicts a reduction in the least value so
!       far of the objective function subject to the linear constraints.
!       Alternatively, a new vector of variables may be chosen to replace
!       an interpolation point that may be too far away for reliability, and
!       then the new point does not have to satisfy the linear constraints.
!       The arguments of the subroutine are as follows.
!
!     N must be set to the number of variables and must be at least two.
!     NPT must be set to the number of interpolation conditions, which is
!       required to be in the interval [N+2,(N+1)(N+2)/2]. Typical choices
!       of the author are NPT=N+6 and NPT=2*N+1. Larger values tend to be
!       highly inefficent when the number of variables is substantial, due
!       to the amount of work and extra difficulty of adjusting more points.
!     M must be set to the number of linear inequality constraints.
!     A is a matrix whose columns are the constraint gradients, which are
!       required to be nonzero.
!     IA is the first dimension of the array A, which must be at least N.
!     B is the vector of right hand sides of the constraints, the J-th
!       constraint being that the scalar product of A(.,J) with X(.) is at
!       most B(J). The initial vector X(.) is made feasible by increasing
!       the value of B(J) if necessary.
!     X is the vector of variables. Initial values of X(1),X(2),...,X(N)
!       must be supplied. If they do not satisfy the constraints, then B
!       is increased as mentioned above. X contains on return the variables
!       that have given the least calculated F subject to the constraints.
!     RHOBEG and RHOEND must be set to the initial and final values of a
!       trust region radius, so both must be positive with RHOEND<=RHOBEG.
!       Typically, RHOBEG should be about one tenth of the greatest expected
!       change to a variable, and RHOEND should indicate the accuracy that
!       is required in the final values of the variables.
!     The value of IPRINT should be set to 0, 1, 2 or 3, which controls the
!       amount of printing. Specifically, there is no output if IPRINT=0 and
!       there is output only at the return if IPRINT=1. Otherwise, the best
!       feasible vector of variables so far and the corresponding value of
!       the objective function are printed whenever RHO is reduced, where
!       RHO is the current lower bound on the trust region radius. Further,
!       each new value of F with its variables are output if IPRINT=3.
!     MAXFUN must be set to an upper bound on the number of calls of CALFUN,
!       its value being at least NPT+1.
!     W is an array used for working space. Its length must be at least
!       M*(2+N) + NPT*(4+N+NPT) + N*(9+3*N) + MAX [ M+3*N, 2*M+N, 2*NPT ].
!       On return, W(1) is set to the final value of F, and W(2) is set to
!       the total number of function evaluations plus 0.5.
!
!     SUBROUTINE CALFUN (N,X,F) has to be provided by the user. It must set
!       F to the value of the objective function for the variables X(1),
!       X(2),...,X(N). The value of the argument F is positive when CALFUN
!       is called if and only if the current X satisfies the constraints
!       to working accuracy.

    integer,dimension(n) :: iact  !to avoid type mismatch error - JW
!
!     Check that N, NPT and MAXFUN are acceptable.
!
      zero = 0.0_wp
      smallx = 1.0e-6_wp * rhoend
      np = n + 1
      nptm = npt - np
      If (n <= 1) Then
         Print 10
10       Format (/ 4 x, 'Return from LINCOA because N is less than 2.')
         Go To 80
      End If
      If (npt < n+2 .Or. npt > ((n+2)*np)/2) Then
         Print 20
20       Format (/ 4 x, 'Return from LINCOA because NPT is not in',&
         ' the required interval.')
         Go To 80
      End If
      If (maxfun <= npt) Then
         Print 30
30       Format (/ 4 x, 'Return from LINCOA because MAXFUN is less',&
         ' than NPT+1.')
         Go To 80
      End If
!
!     Normalize the constraints, and copy the resultant constraint matrix
!       and right hand sides into working space, after increasing the right
!       hand sides if necessary so that the starting point is feasible.
!
      iamat = max0 (m+3*n, 2*m+n, 2*npt) + 1
      ib = iamat + m * n
      iflag = 0
      If (m > 0) Then
         iw = iamat - 1
         Do 60 j = 1, m
            sum = zero
            temp = zero
            Do 40 i = 1, n
               sum = sum + a (i, j) * x (i)
40          temp = temp + a (i, j) ** 2
            If (temp == zero) Then
               Print 50
50             Format (/ 4 x, 'Return from LINCOA because the gradient of',&
               ' a constraint is zero.')
               Go To 80
            End If
            temp = sqrt (temp)
            If (sum-b(j) > smallx*temp) iflag = 1
            w (ib+j-1) = max (b(j), sum) / temp
            Do 60 i = 1, n
               iw = iw + 1
60       w (iw) = a (i, j) / temp
      End If
      If (iflag == 1) Then
         If (iprint > 0) Print 70
70       Format (/ 4 x, 'LINCOA has made the initial X feasible by',&
         ' increasing part(s) of B.')
      End If
!
!     Partition the working space array, so that different parts of it can be
!     treated separately by the subroutine that performs the main calculation.
!
      ndim = npt + n
      ixb = ib + m
      ixp = ixb + n
      ifv = ixp + n * npt
      ixs = ifv + npt
      ixo = ixs + n
      igo = ixo + n
      ihq = igo + n
      ipq = ihq + (n*np) / 2
      ibmat = ipq + npt
      izmat = ibmat + ndim * n
      istp = izmat + npt * nptm
      isp = istp + n
      ixn = isp + npt + npt
      iac = ixn + n
      irc = iac + n
      iqf = irc + m
      irf = iqf + n * n
      ipqw = irf + (n*np) / 2
!
!     The above settings provide a partition of W for subroutine LINCOB.
!
      Call lincob (n, npt, m, w(iamat), w(ib), x, rhobeg, rhoend, &
       iprint, maxfun, w(ixb), w(ixp), w(ifv), w(ixs), w(ixo), w(igo), &
       w(ihq), w(ipq), w(ibmat), w(izmat), ndim, w(istp), w(isp), &
       w(ixn), iact, w(irc), w(iqf), w(irf), w(ipqw), w, calfun)    !--JW mod
       !w(ixn), w(iac), w(irc), w(iqf), w(irf), w(ipqw), w) !--original 
       
80    Return

End Subroutine lincoa

Subroutine lincob (n, npt, m, amat, b, x, rhobeg, rhoend, iprint, &
  maxfun, xbase, xpt, fval, xsav, xopt, gopt, hq, pq, bmat, zmat, ndim, &
  step, sp, xnew, iact, rescon, qfac, rfac, pqw, w, calfun)
      Implicit real(wp) (a-h, o-z)
      Dimension amat (n,*), b (*), x (*), xbase (*), xpt (npt,*), fval &
       (*), xsav (*), xopt (*), gopt (*), hq (*), pq (*), bmat &
       (ndim,*), zmat (npt,*), step (*), sp (*), xnew (*), iact (*), &
       rescon (*), qfac (n,*), rfac (*), pqw (*), w (*)
      procedure(func) :: calfun
!
!     The arguments N, NPT, M, X, RHOBEG, RHOEND, IPRINT and MAXFUN are
!       identical to the corresponding arguments in SUBROUTINE LINCOA.
!     AMAT is a matrix whose columns are the constraint gradients, scaled
!       so that they have unit length.
!     B contains on entry the right hand sides of the constraints, scaled
!       as above, but later B is modified for variables relative to XBASE.
!     XBASE holds a shift of origin that should reduce the contributions
!       from rounding errors to values of the model and Lagrange functions.
!     XPT contains the interpolation point coordinates relative to XBASE.
!     FVAL holds the values of F at the interpolation points.
!     XSAV holds the best feasible vector of variables so far, without any
!       shift of origin.
!     XOPT is set to XSAV-XBASE, which is the displacement from XBASE of
!       the feasible vector of variables that provides the least calculated
!       F so far, this vector being the current trust region centre.
!     GOPT holds the gradient of the quadratic model at XSAV = XBASE+XOPT.
!     HQ holds the explicit second derivatives of the quadratic model.
!     PQ contains the parameters of the implicit second derivatives of the
!       quadratic model.
!     BMAT holds the last N columns of the big inverse matrix H.
!     ZMAT holds the factorization of the leading NPT by NPT submatrix
!       of H, this factorization being ZMAT times Diag(DZ) times ZMAT^T,
!       where the elements of DZ are plus or minus one, as specified by IDZ.
!     NDIM is the first dimension of BMAT and has the value NPT+N.
!     STEP is employed for trial steps from XOPT. It is also used for working
!       space when XBASE is shifted and in PRELIM.
!     SP is reserved for the scalar products XOPT^T XPT(K,.), K=1,2,...,NPT,
!       followed by STEP^T XPT(K,.), K=1,2,...,NPT.
!     XNEW is the displacement from XBASE of the vector of variables for
!       the current calculation of F, except that SUBROUTINE TRSTEP uses it
!       for working space.
!     IACT is an integer array for the indices of the active constraints.
!     RESCON holds useful information about the constraint residuals. Every
!       nonnegative RESCON(J) is the residual of the J-th constraint at the
!       current trust region centre. Otherwise, if RESCON(J) is negative, the
!       J-th constraint holds as a strict inequality at the trust region
!       centre, its residual being at least |RESCON(J)|; further, the value
!       of |RESCON(J)| is at least the current trust region radius DELTA.
!     QFAC is the orthogonal part of the QR factorization of the matrix of
!       active constraint gradients, these gradients being ordered in
!       accordance with IACT. When NACT is less than N, columns are added
!       to QFAC to complete an N by N orthogonal matrix, which is important
!       for keeping calculated steps sufficiently close to the boundaries
!       of the active constraints.
!     RFAC is the upper triangular part of this QR factorization, beginning
!       with the first diagonal element, followed by the two elements in the
!       upper triangular part of the second column and so on.
!     PQW is used for working space, mainly for storing second derivative
!       coefficients of quadratic functions. Its length is NPT+N.
!     The array W is also used for working space. The required number of
!       elements, namely MAX[M+3*N,2*M+N,2*NPT], is set in LINCOA.
!
!     Set some constants.
!
      half = 0.5_wp
      one = 1.0_wp
      tenth = 0.1_wp
      zero = 0.0_wp
      
      np = n + 1
      nh = (n*np) / 2
      nptm = npt - np
!
!     Set the elements of XBASE, XPT, FVAL, XSAV, XOPT, GOPT, HQ, PQ, BMAT,
!       ZMAT and SP for the first iteration. An important feature is that,
!       if the interpolation point XPT(K,.) is not feasible, where K is any
!       integer from [1,NPT], then a change is made to XPT(K,.) if necessary
!       so that the constraint violation is at least 0.2*RHOBEG. Also KOPT
!       is set so that XPT(KOPT,.) is the initial trust region centre.
!
      Call prelim (n, npt, m, amat, b, x, rhobeg, iprint, xbase, xpt, &
       fval, xsav, xopt, gopt, kopt, hq, pq, bmat, zmat, idz, ndim, sp, &
       rescon, step, pqw, w, calfun)
!
!     Begin the iterative procedure.
!
      nf = npt
      fopt = fval (kopt)
      rho = rhobeg
      delta = rho
      ifeas = 0
      nact = 0
      itest = 3
10    knew = 0
      nvala = 0
      nvalb = 0
!
!     Shift XBASE if XOPT may be too far from XBASE. First make the changes
!       to BMAT that do not depend on ZMAT.
!
20    fsave = fopt
      xoptsq = zero
      Do 30 i = 1, n
30    xoptsq = xoptsq + xopt (i) ** 2
      If (xoptsq >= 1.0e4_wp*delta*delta) Then
         qoptsq = 0.25_wp * xoptsq
         Do 50 k = 1, npt
            sum = zero
            Do 40 i = 1, n
40          sum = sum + xpt (k, i) * xopt (i)
            sum = sum - half * xoptsq
            w (npt+k) = sum
            sp (k) = zero
            Do 50 i = 1, n
               xpt (k, i) = xpt (k, i) - half * xopt (i)
               step (i) = bmat (k, i)
               w (i) = sum * xpt (k, i) + qoptsq * xopt (i)
               ip = npt + i
               Do 50 j = 1, i
50       bmat (ip, j) = bmat (ip, j) + step (i) * w (j) + w (i) * step(j)
!
!     Then the revisions of BMAT that depend on ZMAT are calculated.
!
         Do 90 k = 1, nptm
            sumz = zero
            Do 60 i = 1, npt
               sumz = sumz + zmat (i, k)
60          w (i) = w (npt+i) * zmat (i, k)
            Do 80 j = 1, n
               sum = qoptsq * sumz * xopt (j)
               Do 70 i = 1, npt
70             sum = sum + w (i) * xpt (i, j)
               step (j) = sum
               If (k < idz) sum = - sum
               Do 80 i = 1, npt
80          bmat (i, j) = bmat (i, j) + sum * zmat (i, k)
            Do 90 i = 1, n
               ip = i + npt
               temp = step (i)
               If (k < idz) temp = - temp
               Do 90 j = 1, i
90       bmat (ip, j) = bmat (ip, j) + temp * step (j)
!
!     Update the right hand sides of the constraints.
!
         If (m > 0) Then
            Do 110 j = 1, m
               temp = zero
               Do 100 i = 1, n
100            temp = temp + amat (i, j) * xopt (i)
110         b (j) = b (j) - temp
         End If
!
!     The following instructions complete the shift of XBASE, including the
!       changes to the parameters of the quadratic model.
!
         ih = 0
         Do 130 j = 1, n
            w (j) = zero
            Do 120 k = 1, npt
               w (j) = w (j) + pq (k) * xpt (k, j)
120         xpt (k, j) = xpt (k, j) - half * xopt (j)
            Do 130 i = 1, j
               ih = ih + 1
               hq (ih) = hq (ih) + w (i) * xopt (j) + xopt (i) * w (j)
130      bmat (npt+i, j) = bmat (npt+j, i)
         Do 140 j = 1, n
            xbase (j) = xbase (j) + xopt (j)
            xopt (j) = zero
140      xpt (kopt, j) = zero
      End If
!
!     In the case KNEW=0, generate the next trust region step by calling
!       TRSTEP, where SNORM is the current trust region radius initially.
!       The final value of SNORM is the length of the calculated step,
!       except that SNORM is zero on return if the projected gradient is
!       unsuitable for starting the conjugate gradient iterations.
!
      delsav = delta
      ksave = knew
      If (knew == 0) Then
         snorm = delta
         Do 150 i = 1, n
150      xnew (i) = gopt (i)
         Call trstep (n, npt, m, amat, b, xpt, hq, pq, nact, iact, &
          rescon, qfac, rfac, snorm, step, xnew, w, w(m+1), pqw, &
          pqw(np), w(m+np))
!
!     A trust region step is applied whenever its length, namely SNORM, is at
!       least HALF*DELTA. It is also applied if its length is at least 0.1999
!       times DELTA and if a line search of TRSTEP has caused a change to the
!       active set. Otherwise there is a branch below to label 530 or 560.
!
         temp = half * delta
         If (xnew(1) >= half) temp = 0.1999_wp * delta
         If (snorm <= temp) Then
            delta = half * delta
            If (delta <= 1.4_wp*rho) delta = rho
            nvala = nvala + 1
            nvalb = nvalb + 1
            temp = snorm / rho
            If (delsav > rho) temp = one
            If (temp >= half) nvala = zero
            If (temp >= tenth) nvalb = zero
            If (delsav > rho) Go To 530
            If (nvala < 5 .And. nvalb < 3) Go To 530
            If (snorm > zero) ksave = - 1
            Go To 560
         End If
         nvala = zero
         nvalb = zero
!
!     Alternatively, KNEW is positive. Then the model step is calculated
!       within a trust region of radius DEL, after setting the gradient at
!       XBASE and the second derivative parameters of the KNEW-th Lagrange
!       function in W(1) to W(N) and in PQW(1) to PQW(NPT), respectively.
!
      Else
         del = max (tenth*delta, rho)
         Do 160 i = 1, n
160      w (i) = bmat (knew, i)
         Do 170 k = 1, npt
170      pqw (k) = zero
         Do 180 j = 1, nptm
            temp = zmat (knew, j)
            If (j < idz) temp = - temp
            Do 180 k = 1, npt
180      pqw (k) = pqw (k) + temp * zmat (k, j)
         Call qmstep (n, npt, m, amat, b, xpt, xopt, nact, iact, &
          rescon, qfac, kopt, knew, del, step, w, pqw, w(np), w(np+m), &
          ifeas)
      End If
!
!     Set VQUAD to the change to the quadratic model when the move STEP is
!       made from XOPT. If STEP is a trust region step, then VQUAD should be
!       negative. If it is nonnegative due to rounding errors in this case,
!       there is a branch to label 530 to try to improve the model.
!
      vquad = zero
      ih = 0
      Do 190 j = 1, n
         vquad = vquad + step (j) * gopt (j)
         Do 190 i = 1, j
            ih = ih + 1
            temp = step (i) * step (j)
            If (i == j) temp = half * temp
190   vquad = vquad + temp * hq (ih)
      Do 210 k = 1, npt
         temp = zero
         Do 200 j = 1, n
            temp = temp + xpt (k, j) * step (j)
200      sp (npt+k) = temp
210   vquad = vquad + half * pq (k) * temp * temp
      If (ksave == 0 .And. vquad >= zero) Go To 530
!
!     Calculate the next value of the objective function. The difference
!       between the actual new value of F and the value predicted by the
!       model is recorded in DIFF.
!
220   nf = nf + 1
      If (nf > maxfun) Then
         nf = nf - 1
         If (iprint > 0) Print 230
230      Format (/ 4 x, 'Return from LINCOA because CALFUN has been',&
         ' called MAXFUN times.')
         Go To 600
      End If
      xdiff = zero
      Do 240 i = 1, n
         xnew (i) = xopt (i) + step (i)
         x (i) = xbase (i) + xnew (i)
240   xdiff = xdiff + (x(i)-xsav(i)) ** 2
      xdiff = sqrt (xdiff)
      If (ksave ==-1) xdiff = rho
      If (xdiff <= tenth*rho .Or. xdiff >= delta+delta) Then
         ifeas = 0
         If (iprint > 0) Print 250
250      Format (/ 4 x, 'Return from LINCOA because rounding errors',&
         ' prevent reasonable changes to X.')
         Go To 600
      End If
      If (ksave <= 0) ifeas = 1
      f = dfloat (ifeas)
      Call calfun (n, x, f)
      If (iprint == 3) Then
         Print 260, nf, f, (x(i), i=1, n)
260      Format (/ 4 x, 'Function number', i6, '    F =', 1 pd18.10,&
         '    The corresponding X is:' / (2 x, 5d15.6))
      End If
      If (ksave ==-1) Go To 600
      diff = f - fopt - vquad
!
!     If X is feasible, then set DFFALT to the difference between the new
!       value of F and the value predicted by the alternative model.
!
      If (ifeas == 1 .And. itest < 3) Then
         Do 270 k = 1, npt
            pqw (k) = zero
270      w (k) = fval (k) - fval (kopt)
         Do 290 j = 1, nptm
            sum = zero
            Do 280 i = 1, npt
280         sum = sum + w (i) * zmat (i, j)
            If (j < idz) sum = - sum
            Do 290 k = 1, npt
290      pqw (k) = pqw (k) + sum * zmat (k, j)
         vqalt = zero
         Do 310 k = 1, npt
            sum = zero
            Do 300 j = 1, n
300         sum = sum + bmat (k, j) * step (j)
            vqalt = vqalt + sum * w (k)
310      vqalt = vqalt + pqw (k) * sp (npt+k) * (half*sp(npt+k)+sp(k))
         dffalt = f - fopt - vqalt
      End If
      If (itest == 3) Then
         dffalt = diff
         itest = 0
      End If
!
!     Pick the next value of DELTA after a trust region step.
!
      If (ksave == 0) Then
         ratio = (f-fopt) / vquad
         If (ratio <= tenth) Then
            delta = half * delta
         Else If (ratio <= 0.7_wp) Then
            delta = max (half*delta, snorm)
         Else
            temp = sqrt (2.0_wp) * delta
            delta = max (half*delta, snorm+snorm)
            delta = min (delta, temp)
         End If
         If (delta <= 1.4_wp*rho) delta = rho
      End If
!
!     Update BMAT, ZMAT and IDZ, so that the KNEW-th interpolation point
!       can be moved. If STEP is a trust region step, then KNEW is zero at
!       present, but a positive value is picked by subroutine UPDATE.
!
      Call update (n, npt, xpt, bmat, zmat, idz, ndim, sp, step, kopt, &
       knew, pqw, w)
      If (knew == 0) Then
         If (iprint > 0) Print 320
320      Format (/ 4 x, 'Return from LINCOA because the denominator of the updating formula is zero.')
         Go To 600
      End If
!
!     If ITEST is increased to 3, then the next quadratic model is the
!       one whose second derivative matrix is least subject to the new
!       interpolation conditions. Otherwise the new model is constructed
!       by the symmetric Broyden method in the usual way.
!
      If (ifeas == 1) Then
         itest = itest + 1
         If (abs(dffalt) >= tenth*abs(diff)) itest = 0
      End If
!
!     Update the second derivatives of the model by the symmetric Broyden
!       method, using PQW for the second derivative parameters of the new
!       KNEW-th Lagrange function. The contribution from the old parameter
!       PQ(KNEW) is included in the second derivative matrix HQ. W is used
!       later for the gradient of the new KNEW-th Lagrange function.
!
      If (itest < 3) Then
         Do 330 k = 1, npt
330      pqw (k) = zero
         Do 350 j = 1, nptm
            temp = zmat (knew, j)
            If (temp /= zero) Then
               If (j < idz) temp = - temp
               Do 340 k = 1, npt
340            pqw (k) = pqw (k) + temp * zmat (k, j)
            End If
350      Continue
         ih = 0
         Do 360 i = 1, n
            w (i) = bmat (knew, i)
            temp = pq (knew) * xpt (knew, i)
            Do 360 j = 1, i
               ih = ih + 1
360      hq (ih) = hq (ih) + temp * xpt (knew, j)
         pq (knew) = zero
         Do 370 k = 1, npt
370      pq (k) = pq (k) + diff * pqw (k)
      End If
!
!     Include the new interpolation point with the corresponding updates of
!       SP. Also make the changes of the symmetric Broyden method to GOPT at
!       the old XOPT if ITEST is less than 3.
!
      fval (knew) = f
      sp (knew) = sp (kopt) + sp (npt+kopt)
      ssq = zero
      Do 380 i = 1, n
         xpt (knew, i) = xnew (i)
380   ssq = ssq + step (i) ** 2
      sp (npt+knew) = sp (npt+kopt) + ssq
      If (itest < 3) Then
         Do 390 k = 1, npt
            temp = pqw (k) * sp (k)
            Do 390 i = 1, n
390      w (i) = w (i) + temp * xpt (k, i)
         Do 400 i = 1, n
400      gopt (i) = gopt (i) + diff * w (i)
      End If
!
!     Update FOPT, XSAV, XOPT, KOPT, RESCON and SP if the new F is the
!       least calculated value so far with a feasible vector of variables.
!
      If (f < fopt .And. ifeas == 1) Then
         fopt = f
         Do 410 j = 1, n
            xsav (j) = x (j)
410      xopt (j) = xnew (j)
         kopt = knew
         snorm = sqrt (ssq)
         Do 430 j = 1, m
            If (rescon(j) >= delta+snorm) Then
               rescon (j) = snorm - rescon (j)
            Else
               rescon (j) = rescon (j) + snorm
               If (rescon(j)+delta > zero) Then
                  temp = b (j)
                  Do 420 i = 1, n
420               temp = temp - xopt (i) * amat (i, j)
                  temp = max (temp, zero)
                  If (temp >= delta) temp = - temp
                  rescon (j) = temp
               End If
            End If
430      Continue
         Do 440 k = 1, npt
440      sp (k) = sp (k) + sp (npt+k)
!
!     Also revise GOPT when symmetric Broyden updating is applied.
!
         If (itest < 3) Then
            ih = 0
            Do 450 j = 1, n
               Do 450 i = 1, j
                  ih = ih + 1
                  If (i < j) gopt (j) = gopt (j) + hq (ih) * step &
                   (i)
450         gopt (i) = gopt (i) + hq (ih) * step (j)
            Do 460 k = 1, npt
               temp = pq (k) * sp (npt+k)
               Do 460 i = 1, n
460         gopt (i) = gopt (i) + temp * xpt (k, i)
         End If
      End If
!
!     Replace the current model by the least Frobenius norm interpolant if
!       this interpolant gives substantial reductions in the predictions
!       of values of F at feasible points.
!
      If (itest == 3) Then
         Do 470 k = 1, npt
            pq (k) = zero
470      w (k) = fval (k) - fval (kopt)
         Do 490 j = 1, nptm
            sum = zero
            Do 480 i = 1, npt
480         sum = sum + w (i) * zmat (i, j)
            If (j < idz) sum = - sum
            Do 490 k = 1, npt
490      pq (k) = pq (k) + sum * zmat (k, j)
         Do 500 j = 1, n
            gopt (j) = zero
            Do 500 i = 1, npt
500      gopt (j) = gopt (j) + w (i) * bmat (i, j)
         Do 510 k = 1, npt
            temp = pq (k) * sp (k)
            Do 510 i = 1, n
510      gopt (i) = gopt (i) + temp * xpt (k, i)
         Do 520 ih = 1, nh
520      hq (ih) = zero
      End If
!
!     If a trust region step has provided a sufficient decrease in F, then
!       branch for another trust region calculation. Every iteration that
!       takes a model step is followed by an attempt to take a trust region
!       step.
!
      knew = 0
      If (ksave > 0) Go To 20
      If (ratio >= tenth) Go To 20
!
!     Alternatively, find out if the interpolation points are close enough
!       to the best point so far.
!
530   distsq = max (delta*delta, 4.0_wp*rho*rho)
      Do 550 k = 1, npt
         sum = zero
         Do 540 j = 1, n
540      sum = sum + (xpt(k, j)-xopt(j)) ** 2
         If (sum > distsq) Then
            knew = k
            distsq = sum
         End If
550   Continue
!
!     If KNEW is positive, then branch back for the next iteration, which
!       will generate a "model step". Otherwise, if the current iteration
!       has reduced F, or if DELTA was above its lower bound when the last
!       trust region step was calculated, then try a "trust region" step
!       instead.
!
      If (knew > 0) Go To 20
      knew = 0
      If (fopt < fsave) Go To 20
      If (delsav > rho) Go To 20
!
!     The calculations with the current value of RHO are complete.
!       Pick the next value of RHO.
!
560   If (rho > rhoend) Then
         delta = half * rho
         If (rho > 250.0_wp*rhoend) Then
            rho = tenth * rho
         Else If (rho <= 16.0_wp*rhoend) Then
            rho = rhoend
         Else
            rho = sqrt (rho*rhoend)
         End If
         delta = max (delta, rho)
         If (iprint >= 2) Then
            If (iprint >= 3) Print 570
570         Format (5 x)
            Print 580, rho, nf
580         Format (/ 4 x, 'New RHO =', 1 pd11.4, 5 x, 'Number of',&
            ' function values =', i6)
            Print 590, fopt, (xbase(i)+xopt(i), i=1, n)
590         Format (4 x, 'Least value of F =', 1 pd23.15, 9 x,&
            'The corresponding X is:'/(2 x, 5d15.6))
         End If
         Go To 10
      End If
!
!     Return from the calculation, after branching to label 220 for another
!       Newton-Raphson step if it has not been tried before.
!
      If (ksave ==-1) Go To 220
600   If (fopt <= f .Or. ifeas == 0) Then
         Do 610 i = 1, n
610      x (i) = xsav (i)
         f = fopt
      End If
      If (iprint >= 1) Then
         Print 620, nf
620      Format (/ 4 x, 'At the return from LINCOA', 5 x,&
         'Number of function values =', i6)
         Print 590, f, (x(i), i=1, n)
      End If
      w (1) = f
      w (2) = dfloat (nf) + half
      Return
End Subroutine lincob

    
Subroutine getact (n, m, amat, b, nact, iact, qfac, rfac, snorm, &
  resnew, resact, g, dw, vlam, w)
      Implicit real(wp) (a-h, o-z)
      Dimension amat (n,*), b (*), iact (*), qfac (n,*), rfac (*), &
       resnew (*), resact (*), g (*), dw (*), vlam (*), w (*)
!
!     N, M, AMAT, B, NACT, IACT, QFAC and RFAC are the same as the terms
!       with these names in SUBROUTINE LINCOB. The current values must be
!       set on entry. NACT, IACT, QFAC and RFAC are kept up to date when
!       GETACT changes the current active set.
!     SNORM, RESNEW, RESACT, G and DW are the same as the terms with these
!       names in SUBROUTINE TRSTEP. The elements of RESNEW and RESACT are
!       also kept up to date.
!     VLAM and W are used for working space, the vector VLAM being reserved
!       for the Lagrange multipliers of the calculation. Their lengths must
!       be at least N.
!     The main purpose of GETACT is to pick the current active set. It is
!       defined by the property that the projection of -G into the space
!       orthogonal to the active constraint normals is as large as possible,
!       subject to this projected steepest descent direction moving no closer
!       to the boundary of every constraint whose current residual is at most
!       0.2*SNORM. On return, the settings in NACT, IACT, QFAC and RFAC are
!       all appropriate to this choice of active set.
!     Occasionally this projected direction is zero, and then the final value
!       of W(1) is set to zero. Otherwise, the direction itself is returned
!       in DW, and W(1) is set to the square of the length of the direction.
!
!     Set some constants and a temporary VLAM.
!
      one = 1.0_wp
      tiny = 1.0e-60_wp
      zero = 0.0_wp
      tdel = 0.2_wp * snorm
      ddsav = zero
      Do 10 i = 1, n
         ddsav = ddsav + g (i) ** 2
10    vlam (i) = zero
      ddsav = ddsav + ddsav
!
!     Set the initial QFAC to the identity matrix in the case NACT=0.
!
      If (nact == 0) Then
         Do 30 i = 1, n
            Do 20 j = 1, n
20          qfac (i, j) = zero
30       qfac (i, i) = one
         Go To 100
      End If
!
!     Remove any constraints from the initial active set whose residuals
!       exceed TDEL.
!
      iflag = 1
      ic = nact
40    If (resact(ic) > tdel) Go To 800
50    ic = ic - 1
      If (ic > 0) Go To 40
!
!     Remove any constraints from the initial active set whose Lagrange
!       multipliers are nonnegative, and set the surviving multipliers.
!
      iflag = 2
60    If (nact == 0) Go To 100
      ic = nact
70    temp = zero
      Do 80 i = 1, n
80    temp = temp + qfac (i, ic) * g (i)
      idiag = (ic*ic+ic) / 2
      If (ic < nact) Then
         jw = idiag + ic
         Do 90 j = ic + 1, nact
            temp = temp - rfac (jw) * vlam (j)
90       jw = jw + j
      End If
      If (temp >= zero) Go To 800
      vlam (ic) = temp / rfac (idiag)
      ic = ic - 1
      If (ic > 0) Go To 70
!
!     Set the new search direction D. Terminate if the 2-norm of D is zero
!       or does not decrease, or if NACT=N holds. The situation NACT=N
!       occurs for sufficiently large SNORM if the origin is in the convex
!       hull of the constraint gradients.
!
100   If (nact == n) Go To 290
      Do 110 j = nact + 1, n
         w (j) = zero
         Do 110 i = 1, n
110   w (j) = w (j) + qfac (i, j) * g (i)
      dd = zero
      Do 130 i = 1, n
         dw (i) = zero
         Do 120 j = nact + 1, n
120      dw (i) = dw (i) - w (j) * qfac (i, j)
130   dd = dd + dw (i) ** 2
      If (dd >= ddsav) Go To 290
      If (dd == zero) Go To 300
      ddsav = dd
      dnorm = sqrt (dd)
!
!     Pick the next integer L or terminate, a positive value of L being
!       the index of the most violated constraint. The purpose of CTOL
!       below is to estimate whether a positive value of VIOLMX may be
!       due to computer rounding errors.
!
      l = 0
      If (m > 0) Then
         test = dnorm / snorm
         violmx = zero
         Do 150 j = 1, m
            If (resnew(j) > zero .And. resnew(j) <= tdel) Then
               sum = zero
               Do 140 i = 1, n
140            sum = sum + amat (i, j) * dw (i)
               If (sum > test*resnew(j)) Then
                  If (sum > violmx) Then
                     l = j
                     violmx = sum
                  End If
               End If
            End If
150      Continue
         ctol = zero
         temp = 0.01_wp * dnorm
         If (violmx > zero .And. violmx < temp) Then
            If (nact > 0) Then
               Do 170 k = 1, nact
                  j = iact (k)
                  sum = zero
                  Do 160 i = 1, n
160               sum = sum + dw (i) * amat (i, j)
170            ctol = max (ctol, abs(sum))
            End If
         End If
      End If
      w (1) = one
      If (l == 0) Go To 300
      If (violmx <= 10.0_wp*ctol) Go To 300
!
!     Apply Givens rotations to the last (N-NACT) columns of QFAC so that
!       the first (NACT+1) columns of QFAC are the ones required for the
!       addition of the L-th constraint, and add the appropriate column
!       to RFAC.
!
      nactp = nact + 1
      idiag = (nactp*nactp-nactp) / 2
      rdiag = zero
      Do 200 j = n, 1, - 1
         sprod = zero
         Do 180 i = 1, n
180      sprod = sprod + qfac (i, j) * amat (i, l)
         If (j <= nact) Then
            rfac (idiag+j) = sprod
         Else
            If (abs(rdiag) <= 1.0e-20_wp*abs(sprod)) Then
               rdiag = sprod
            Else
               temp = sqrt (sprod*sprod+rdiag*rdiag)
               cosv = sprod / temp
               sinv = rdiag / temp
               rdiag = temp
               Do 190 i = 1, n
                  temp = cosv * qfac (i, j) + sinv * qfac (i, j+1)
                  qfac (i, j+1) = - sinv * qfac (i, j) + cosv * qfac(i, j+1)
190            qfac (i, j) = temp
            End If
         End If
200   Continue
      If (rdiag < zero) Then
         Do 210 i = 1, n
210      qfac (i, nactp) = - qfac (i, nactp)
      End If
      rfac (idiag+nactp) = abs (rdiag)
      nact = nactp
      iact (nact) = l
      resact (nact) = resnew (l)
      vlam (nact) = zero
      resnew (l) = zero
!
!     Set the components of the vector VMU in W.
!
220   w (nact) = one / rfac ((nact*nact+nact)/2) ** 2
      If (nact > 1) Then
         Do 240 i = nact - 1, 1, - 1
            idiag = (i*i+i) / 2
            jw = idiag + i
            sum = zero
            Do 230 j = i + 1, nact
               sum = sum - rfac (jw) * w (j)
230         jw = jw + j
240      w (i) = sum / rfac (idiag)
      End If
!
!     Calculate the multiple of VMU to subtract from VLAM, and update VLAM.
!
      vmult = violmx
      ic = 0
      j = 1
250   If (j < nact) Then
         If (vlam(j) >= vmult*w(j)) Then
            ic = j
            vmult = vlam (j) / w (j)
         End If
         j = j + 1
         Go To 250
      End If
      Do 260 j = 1, nact
260   vlam (j) = vlam (j) - vmult * w (j)
      If (ic > 0) vlam (ic) = zero
      violmx = max (violmx-vmult, zero)
      If (ic == 0) violmx = zero
!
!     Reduce the active set if necessary, so that all components of the
!       new VLAM are negative, with resetting of the residuals of the
!       constraints that become inactive.
!
      iflag = 3
      ic = nact
270   If (vlam(ic) < zero) Go To 280
      resnew (iact(ic)) = max (resact(ic), tiny)
      Go To 800
280   ic = ic - 1
      If (ic > 0) Go To 270
!
!     Calculate the next VMU if VIOLMX is positive. Return if NACT=N holds,
!       as then the active constraints imply D=0. Otherwise, go to label
!       100, to calculate the new D and to test for termination.
!
      If (violmx > zero) Go To 220
      If (nact < n) Go To 100
290   dd = zero
300   w (1) = dd
      Return
!
!     These instructions rearrange the active constraints so that the new
!       value of IACT(NACT) is the old value of IACT(IC). A sequence of
!       Givens rotations is applied to the current QFAC and RFAC. Then NACT
!       is reduced by one.
!
800   resnew (iact(ic)) = max (resact(ic), tiny)
      jc = ic
810   If (jc < nact) Then
         jcp = jc + 1
         idiag = jc * jcp / 2
         jw = idiag + jcp
         temp = sqrt (rfac(jw-1)**2+rfac(jw)**2)
         cval = rfac (jw) / temp
         sval = rfac (jw-1) / temp
         rfac (jw-1) = sval * rfac (idiag)
         rfac (jw) = cval * rfac (idiag)
         rfac (idiag) = temp
         If (jcp < nact) Then
            Do 820 j = jcp + 1, nact
               temp = sval * rfac (jw+jc) + cval * rfac (jw+jcp)
               rfac (jw+jcp) = cval * rfac (jw+jc) - sval * rfac(jw+jcp)
               rfac (jw+jc) = temp
820         jw = jw + j
         End If
         jdiag = idiag - jc
         Do 830 i = 1, n
            If (i < jc) Then
               temp = rfac (idiag+i)
               rfac (idiag+i) = rfac (jdiag+i)
               rfac (jdiag+i) = temp
            End If
            temp = sval * qfac (i, jc) + cval * qfac (i, jcp)
            qfac (i, jcp) = cval * qfac (i, jc) - sval * qfac (i, jcp)
830      qfac (i, jc) = temp
         iact (jc) = iact (jcp)
         resact (jc) = resact (jcp)
         vlam (jc) = vlam (jcp)
         jc = jcp
         Go To 810
      End If
      nact = nact - 1
      Go To (50, 60, 280), iflag
End Subroutine getact

Subroutine prelim (n, npt, m, amat, b, x, rhobeg, iprint, xbase, xpt, &
  fval, xsav, xopt, gopt, kopt, hq, pq, bmat, zmat, idz, ndim, sp, &
  rescon, step, pqw, w, calfun)
   Implicit real(wp) (a-h, o-z)
   Dimension amat (n,*), b (*), x (*), xbase (*), xpt (npt,*), fval &
    (*), xsav (*), xopt (*), gopt (*), hq (*), pq (*), bmat (ndim,*), &
    zmat (npt,*), sp (*), rescon (*), step (*), pqw (*), w (*)
   procedure(func) :: calfun
!
!     The arguments N, NPT, M, AMAT, B, X, RHOBEG, IPRINT, XBASE, XPT, FVAL,
!       XSAV, XOPT, GOPT, HQ, PQ, BMAT, ZMAT, NDIM, SP and RESCON are the
!       same as the corresponding arguments in SUBROUTINE LINCOB.
!     KOPT is set to the integer such that XPT(KOPT,.) is the initial trust
!       region centre.
!     IDZ is going to be set to one, so that every element of Diag(DZ) is
!       one in the product ZMAT times Diag(DZ) times ZMAT^T, which is the
!       factorization of the leading NPT by NPT submatrix of H.
!     STEP, PQW and W are used for working space, the arrays STEP and PQW
!       being taken from LINCOB. The length of W must be at least N+NPT.
!
!     SUBROUTINE PRELIM provides the elements of XBASE, XPT, BMAT and ZMAT
!       for the first iteration, an important feature being that, if any of
!       of the columns of XPT is an infeasible point, then the largest of
!       the constraint violations there is at least 0.2*RHOBEG. It also sets
!       the initial elements of FVAL, XOPT, GOPT, HQ, PQ, SP and RESCON.
!
!     Set some constants.
!
   half = 0.5_wp
   one = 1.0_wp
   zero = 0.0_wp
   nptm = npt - n - 1
   rhosq = rhobeg * rhobeg
   recip = one / rhosq
   reciq = sqrt (half) / rhosq
   test = 0.2_wp * rhobeg
   idz = 1
   kbase = 1
!
!     Set the initial elements of XPT, BMAT, SP and ZMAT to zero.
!
   Do 20 j = 1, n
      xbase (j) = x (j)
      Do 10 k = 1, npt
10    xpt (k, j) = zero
      Do 20 i = 1, ndim
20 bmat (i, j) = zero
   Do 30 k = 1, npt
      sp (k) = zero
      Do 30 j = 1, npt - n - 1
30 zmat (k, j) = zero
!
!     Set the nonzero coordinates of XPT(K,.), K=1,2,...,min[2*N+1,NPT],
!       but they may be altered later to make a constraint violation
!       sufficiently large. The initial nonzero elements of BMAT and of
!       the first min[N,NPT-N-1] columns of ZMAT are set also.
!
   Do 40 j = 1, n
      xpt (j+1, j) = rhobeg
      If (j < npt-n) Then
         jp = n + j + 1
         xpt (jp, j) = - rhobeg
         bmat (j+1, j) = half / rhobeg
         bmat (jp, j) = - half / rhobeg
         zmat (1, j) = - reciq - reciq
         zmat (j+1, j) = reciq
         zmat (jp, j) = reciq
      Else
         bmat (1, j) = - one / rhobeg
         bmat (j+1, j) = one / rhobeg
         bmat (npt+j, j) = - half * rhosq
      End If
40 Continue
!
!     Set the remaining initial nonzero elements of XPT and ZMAT when the
!       number of interpolation points exceeds 2*N+1.
!
   If (npt > 2*n+1) Then
      Do 50 k = n + 1, npt - n - 1
         itemp = (k-1) / n
         ipt = k - itemp * n
         jpt = ipt + itemp
         If (jpt > n) jpt = jpt - n
         xpt (n+k+1, ipt) = rhobeg
         xpt (n+k+1, jpt) = rhobeg
         zmat (1, k) = recip
         zmat (ipt+1, k) = - recip
         zmat (jpt+1, k) = - recip
50    zmat (n+k+1, k) = recip
   End If
!
!     Update the constraint right hand sides to allow for the shift XBASE.
!
   If (m > 0) Then
      Do 70 j = 1, m
         temp = zero
         Do 60 i = 1, n
60       temp = temp + amat (i, j) * xbase (i)
70    b (j) = b (j) - temp
   End If
!
!     Go through the initial points, shifting every infeasible point if
!       necessary so that its constraint violation is at least 0.2*RHOBEG.
!
   Do 150 nf = 1, npt
      feas = one
      bigv = zero
      j = 0
80    j = j + 1
      If (j <= m .And. nf >= 2) Then
         resid = - b (j)
         Do 90 i = 1, n
90       resid = resid + xpt (nf, i) * amat (i, j)
         If (resid <= bigv) Go To 80
         bigv = resid
         jsav = j
         If (resid <= test) Then
            feas = - one
            Go To 80
         End If
         feas = zero
      End If
      If (feas < zero) Then
         Do 100 i = 1, n
100      step (i) = xpt (nf, i) + (test-bigv) * amat (i, jsav)
         Do 110 k = 1, npt
            sp (npt+k) = zero
            Do 110 j = 1, n
110      sp (npt+k) = sp (npt+k) + xpt (k, j) * step (j)
         Call update (n, npt, xpt, bmat, zmat, idz, ndim, sp, step, &
          kbase, nf, pqw, w)
         Do 120 i = 1, n
120      xpt (nf, i) = step (i)
      End If
!
!     Calculate the objective function at the current interpolation point,
!       and set KOPT to the index of the first trust region centre.
!
      Do 130 j = 1, n
130   x (j) = xbase (j) + xpt (nf, j)
      f = feas
      Call calfun (n, x, f)
      If (iprint == 3) Then
         Print 140, nf, f, (x(i), i=1, n)
140      Format (/ 4 x, 'Function number', i6, '    F =', 1 pd18.10,&
         '    The corresponding X is:' / (2 x, 5d15.6))
      End If
      If (nf == 1) Then
         kopt = 1
      Else If (f < fval(kopt) .And. feas > zero) Then
         kopt = nf
      End If
150 fval (nf) = f
!
!     Set PQ for the first quadratic model.
!
   Do 160 j = 1, nptm
      w (j) = zero
      Do 160 k = 1, npt
160 w (j) = w (j) + zmat (k, j) * fval (k)
   Do 170 k = 1, npt
      pq (k) = zero
      Do 170 j = 1, nptm
170 pq (k) = pq (k) + zmat (k, j) * w (j)
!
!     Set XOPT, SP, GOPT and HQ for the first quadratic model.
!
   Do 180 j = 1, n
      xopt (j) = xpt (kopt, j)
      xsav (j) = xbase (j) + xopt (j)
180 gopt (j) = zero
   Do 200 k = 1, npt
      sp (k) = zero
      Do 190 j = 1, n
190   sp (k) = sp (k) + xpt (k, j) * xopt (j)
      temp = pq (k) * sp (k)
      Do 200 j = 1, n
200 gopt (j) = gopt (j) + fval (k) * bmat (k, j) + temp * xpt (k, j)
   Do 210 i = 1, (n*n+n) / 2
210 hq (i) = zero
!
!     Set the initial elements of RESCON.
!
   Do 230 j = 1, m
      temp = b (j)
      Do 220 i = 1, n
220   temp = temp - xopt (i) * amat (i, j)
      temp = max (temp, zero)
      If (temp >= rhobeg) temp = - temp
230 rescon (j) = temp
   Return
End Subroutine prelim

Subroutine qmstep (n, npt, m, amat, b, xpt, xopt, nact, iact, rescon, &
  qfac, kopt, knew, del, step, gl, pqw, rstat, w, ifeas)
   Implicit real(wp) (a-h, o-z)
   Dimension amat (n,*), b (*), xpt (npt,*), xopt (*), iact (*), rescon &
    (*), qfac (n,*), step (*), gl (*), pqw (*), rstat (*), w (*)
!
!     N, NPT, M, AMAT, B, XPT, XOPT, NACT, IACT, RESCON, QFAC, KOPT are the
!       same as the terms with these names in SUBROUTINE LINCOB.
!     KNEW is the index of the interpolation point that is going to be moved.
!     DEL is the current restriction on the length of STEP, which is never
!       greater than the current trust region radius DELTA.
!     STEP will be set to the required step from XOPT to the new point.
!     GL must be set on entry to the gradient of LFUNC at XBASE, where LFUNC
!       is the KNEW-th Lagrange function. It is used also for some other
!       gradients of LFUNC.
!     PQW provides the second derivative parameters of LFUNC.
!     RSTAT and W are used for working space. Their lengths must be at least
!       M and N, respectively. RSTAT(J) is set to -1.0, 0.0, or 1.0 if the
!       J-th constraint is irrelevant, active, or both inactive and relevant,
!       respectively.
!     IFEAS will be set to 0 or 1 if XOPT+STEP is infeasible or feasible.
!
!     STEP is chosen to provide a relatively large value of the modulus of
!       LFUNC(XOPT+STEP), subject to ||STEP|| .LE. DEL. A projected STEP is
!       calculated too, within the trust region, that does not alter the
!       residuals of the active constraints. The projected step is preferred
!       if its value of | LFUNC(XOPT+STEP) | is at least one fifth of the
!       original one, but the greatest violation of a linear constraint must
!       be at least 0.2*DEL, in order to keep the interpolation points apart.
!       The remedy when the maximum constraint violation is too small is to
!       restore the original step, which is perturbed if necessary so that
!       its maximum constraint violation becomes 0.2*DEL.
!
!     Set some constants.
!
   half = 0.5_wp
   one = 1.0_wp
   tenth = 0.1_wp
   zero = 0.0_wp
   test = 0.2_wp * del
!
!     Replace GL by the gradient of LFUNC at the trust region centre, and
!       set the elements of RSTAT.
!
   Do 20 k = 1, npt
      temp = zero
      Do 10 j = 1, n
10    temp = temp + xpt (k, j) * xopt (j)
      temp = pqw (k) * temp
      Do 20 i = 1, n
20 gl (i) = gl (i) + temp * xpt (k, i)
   If (m > 0) Then
      Do 30 j = 1, m
         rstat (j) = one
30    If (abs(rescon(j)) >= del) rstat (j) = - one
      Do 40 k = 1, nact
40    rstat (iact(k)) = zero
   End If
!
!     Find the greatest modulus of LFUNC on a line through XOPT and
!       another interpolation point within the trust region.
!
   iflag = 0
   vbig = zero
   Do 60 k = 1, npt
      If (k == kopt) Go To 60
      ss = zero
      sp = zero
      Do 50 i = 1, n
         temp = xpt (k, i) - xopt (i)
         ss = ss + temp * temp
50    sp = sp + gl (i) * temp
      stp = - del / sqrt (ss)
      If (k == knew) Then
         If (sp*(sp-one) < zero) stp = - stp
         vlag = abs (stp*sp) + stp * stp * abs (sp-one)
      Else
         vlag = abs (stp*(one-stp)*sp)
      End If
      If (vlag > vbig) Then
         ksav = k
         stpsav = stp
         vbig = vlag
      End If
60 Continue
!
!     Set STEP to the move that gives the greatest modulus calculated above.
!       This move may be replaced by a steepest ascent step from XOPT.
!
   gg = zero
   Do 70 i = 1, n
      gg = gg + gl (i) ** 2
70 step (i) = stpsav * (xpt(ksav, i)-xopt(i))
   vgrad = del * sqrt (gg)
   If (vgrad <= tenth*vbig) Go To 220
!
!     Make the replacement if it provides a larger value of VBIG.
!
   ghg = zero
   Do 90 k = 1, npt
      temp = zero
      Do 80 j = 1, n
80    temp = temp + xpt (k, j) * gl (j)
90 ghg = ghg + pqw (k) * temp * temp
   vnew = vgrad + abs (half*del*del*ghg/gg)
   If (vnew > vbig) Then
      vbig = vnew
      stp = del / sqrt (gg)
      If (ghg < zero) stp = - stp
      Do 100 i = 1, n
100   step (i) = stp * gl (i)
   End If
   If (nact == 0 .Or. nact == n) Go To 220
!
!     Overwrite GL by its projection. Then set VNEW to the greatest
!       value of |LFUNC| on the projected gradient from XOPT subject to
!       the trust region bound. If VNEW is sufficiently large, then STEP
!       may be changed to a move along the projected gradient.
!
   Do 110 k = nact + 1, n
      w (k) = zero
      Do 110 i = 1, n
110 w (k) = w (k) + gl (i) * qfac (i, k)
   gg = zero
   Do 130 i = 1, n
      gl (i) = zero
      Do 120 k = nact + 1, n
120   gl (i) = gl (i) + qfac (i, k) * w (k)
130 gg = gg + gl (i) ** 2
   vgrad = del * sqrt (gg)
   If (vgrad <= tenth*vbig) Go To 220
   ghg = zero
   Do 150 k = 1, npt
      temp = zero
      Do 140 j = 1, n
140   temp = temp + xpt (k, j) * gl (j)
150 ghg = ghg + pqw (k) * temp * temp
   vnew = vgrad + abs (half*del*del*ghg/gg)
!
!     Set W to the possible move along the projected gradient.
!
   stp = del / sqrt (gg)
   If (ghg < zero) stp = - stp
   ww = zero
   Do 160 i = 1, n
      w (i) = stp * gl (i)
160 ww = ww + w (i) ** 2
!
!     Set STEP to W if W gives a sufficiently large value of the modulus
!       of the Lagrange function, and if W either preserves feasibility
!       or gives a constraint violation of at least 0.2*DEL. The purpose
!       of CTOL below is to provide a check on feasibility that includes
!       a tolerance for contributions from computer rounding errors.
!
   If (vnew/vbig >= 0.2_wp) Then
      ifeas = 1
      bigv = zero
      j = 0
170   j = j + 1
      If (j <= m) Then
         If (rstat(j) == one) Then
            temp = - rescon (j)
            Do 180 i = 1, n
180         temp = temp + w (i) * amat (i, j)
            bigv = max (bigv, temp)
         End If
         If (bigv < test) Go To 170
         ifeas = 0
      End If
      ctol = zero
      temp = 0.01_wp * sqrt (ww)
      If (bigv > zero .And. bigv < temp) Then
         Do 200 k = 1, nact
            j = iact (k)
            sum = zero
            Do 190 i = 1, n
190         sum = sum + w (i) * amat (i, j)
200      ctol = max (ctol, abs(sum))
      End If
      If (bigv <= 10.0_wp*ctol .Or. bigv >= test) Then
         Do 210 i = 1, n
210      step (i) = w (i)
         Go To 260
      End If
   End If
!
!     Calculate the greatest constraint violation at XOPT+STEP with STEP at
!       its original value. Modify STEP if this violation is unacceptable.
!
220 ifeas = 1
   bigv = zero
   resmax = zero
   j = 0
230 j = j + 1
   If (j <= m) Then
      If (rstat(j) < zero) Go To 230
      temp = - rescon (j)
      Do 240 i = 1, n
240   temp = temp + step (i) * amat (i, j)
      resmax = max (resmax, temp)
      If (temp < test) Then
         If (temp <= bigv) Go To 230
         bigv = temp
         jsav = j
         ifeas = - 1
         Go To 230
      End If
      ifeas = 0
   End If
   If (ifeas ==-1) Then
      Do 250 i = 1, n
250   step (i) = step (i) + (test-bigv) * amat (i, jsav)
      ifeas = 0
   End If
!
!     Return the calculated STEP and the value of IFEAS.
!
260 Return
End Subroutine qmstep

Subroutine trstep (n, npt, m, amat, b, xpt, hq, pq, nact, iact, rescon, &
  qfac, rfac, snorm, step, g, resnew, resact, d, dw, w)
   Implicit real(wp) (a-h, o-z)
   Dimension amat (n,*), b (*), xpt (npt,*), hq (*), pq (*), iact (*), &
    rescon (*), qfac (n,*), rfac (*), step (*), g (*), resnew (*), &
    resact (*), d (*), dw (*), w (*)
!
!     N, NPT, M, AMAT, B, XPT, HQ, PQ, NACT, IACT, RESCON, QFAC and RFAC
!       are the same as the terms with these names in LINCOB. If RESCON(J)
!       is negative, then |RESCON(J)| must be no less than the trust region
!       radius, so that the J-th constraint can be ignored.
!     SNORM is set to the trust region radius DELTA initially. On the
!       return, however, it is the length of the calculated STEP, which is
!       set to zero if the constraints do not allow a long enough step.
!     STEP is the total calculated step so far from the trust region centre,
!       its final value being given by the sequence of CG iterations, which
!       terminate if the trust region boundary is reached.
!     G must be set on entry to the gradient of the quadratic model at the
!       trust region centre. It is used as working space, however, and is
!       always the gradient of the model at the current STEP, except that
!       on return the value of G(1) is set to ONE instead of to ZERO if
!       and only if GETACT is called more than once.
!     RESNEW, RESACT, D, DW and W are used for working space. A negative
!       value of RESNEW(J) indicates that the J-th constraint does not
!       restrict the CG steps of the current trust region calculation, a
!       zero value of RESNEW(J) indicates that the J-th constraint is active,
!       and otherwise RESNEW(J) is set to the greater of TINY and the actual
!       residual of the J-th constraint for the current STEP. RESACT holds
!       the residuals of the active constraints, which may be positive.
!       D is the search direction of each line search. DW is either another
!       search direction or the change in gradient along D. The length of W
!       must be at least MAX[M,2*N].
!
!     Set some numbers for the conjugate gradient iterations.
!
   half = 0.5_wp
   one = 1.0_wp
   tiny = 1.0e-60_wp
   zero = 0.0_wp
   ctest = 0.01_wp
   snsq = snorm * snorm
!
!     Set the initial elements of RESNEW, RESACT and STEP.
!
   If (m > 0) Then
      Do 10 j = 1, m
         resnew (j) = rescon (j)
         If (rescon(j) >= snorm) Then
            resnew (j) = - one
         Else If (rescon(j) >= zero) Then
            resnew (j) = max (resnew(j), tiny)
         End If
10    Continue
      If (nact > 0) Then
         Do 20 k = 1, nact
            resact (k) = rescon (iact(k))
20       resnew (iact(k)) = zero
      End If
   End If
   Do 30 i = 1, n
30 step (i) = zero
   ss = zero
   reduct = zero
   ncall = 0
!
!     GETACT picks the active set for the current STEP. It also sets DW to
!       the vector closest to -G that is orthogonal to the normals of the
!       active constraints. DW is scaled to have length 0.2*SNORM, as then
!       a move of DW from STEP is allowed by the linear constraints.
!
40 ncall = ncall + 1
   Call getact (n, m, amat, b, nact, iact, qfac, rfac, snorm, resnew, &
  & resact, g, dw, w, w(n+1))
   If (w(n+1) == zero) Go To 320
   scale = 0.2_wp * snorm / sqrt (w(n+1))
   Do 50 i = 1, n
50 dw (i) = scale * dw (i)
!
!     If the modulus of the residual of an active constraint is substantial,
!       then set D to the shortest move from STEP to the boundaries of the
!       active constraints.
!
   resmax = zero
   If (nact > 0) Then
      Do 60 k = 1, nact
60    resmax = max (resmax, resact(k))
   End If
   gamma = zero
   If (resmax > 1.0e-4_wp*snorm) Then
      ir = 0
      Do 80 k = 1, nact
         temp = resact (k)
         If (k >= 2) Then
            Do 70 i = 1, k - 1
               ir = ir + 1
70          temp = temp - rfac (ir) * w (i)
         End If
         ir = ir + 1
80    w (k) = temp / rfac (ir)
      Do 90 i = 1, n
         d (i) = zero
         Do 90 k = 1, nact
90    d (i) = d (i) + w (k) * qfac (i, k)
!
!     The vector D that has just been calculated is also the shortest move
!       from STEP+DW to the boundaries of the active constraints. Set GAMMA
!       to the greatest steplength of this move that satisfies the trust
!       region bound.
!
      rhs = snsq
      ds = zero
      dd = zero
      Do 100 i = 1, n
         sum = step (i) + dw (i)
         rhs = rhs - sum * sum
         ds = ds + d (i) * sum
100   dd = dd + d (i) ** 2
      If (rhs > zero) Then
         temp = sqrt (ds*ds+dd*rhs)
         If (ds <= zero) Then
            gamma = (temp-ds) / dd
         Else
            gamma = rhs / (temp+ds)
         End If
      End If
!
!     Reduce the steplength GAMMA if necessary so that the move along D
!       also satisfies the linear constraints.
!
      j = 0
110   If (gamma > zero) Then
         j = j + 1
         If (resnew(j) > zero) Then
            ad = zero
            adw = zero
            Do 120 i = 1, n
               ad = ad + amat (i, j) * d (i)
120         adw = adw + amat (i, j) * dw (i)
            If (ad > zero) Then
               temp = max ((resnew(j)-adw)/ad, zero)
               gamma = min (gamma, temp)
            End If
         End If
         If (j < m) Go To 110
      End If
      gamma = min (gamma, one)
   End If
!
!     Set the next direction for seeking a reduction in the model function
!       subject to the trust region bound and the linear constraints.
!
   If (gamma <= zero) Then
      Do 130 i = 1, n
130   d (i) = dw (i)
      icount = nact
   Else
      Do 140 i = 1, n
140   d (i) = dw (i) + gamma * d (i)
      icount = nact - 1
   End If
   alpbd = one
!
!     Set ALPHA to the steplength from STEP along D to the trust region
!       boundary. Return if the first derivative term of this step is
!       sufficiently small or if no further progress is possible.
!
150 icount = icount + 1
   rhs = snsq - ss
   If (rhs <= zero) Go To 320
   dg = zero
   ds = zero
   dd = zero
   Do 160 i = 1, n
      dg = dg + d (i) * g (i)
      ds = ds + d (i) * step (i)
160 dd = dd + d (i) ** 2
   If (dg >= zero) Go To 320
   temp = sqrt (rhs*dd+ds*ds)
   If (ds <= zero) Then
      alpha = (temp-ds) / dd
   Else
      alpha = rhs / (temp+ds)
   End If
   If (-alpha*dg <= ctest*reduct) Go To 320
!
!     Set DW to the change in gradient along D.
!
   ih = 0
   Do 170 j = 1, n
      dw (j) = zero
      Do 170 i = 1, j
         ih = ih + 1
         If (i < j) dw (j) = dw (j) + hq (ih) * d (i)
170 dw (i) = dw (i) + hq (ih) * d (j)
   Do 190 k = 1, npt
      temp = zero
      Do 180 j = 1, n
180   temp = temp + xpt (k, j) * d (j)
      temp = pq (k) * temp
      Do 190 i = 1, n
190 dw (i) = dw (i) + temp * xpt (k, i)
!
!     Set DGD to the curvature of the model along D. Then reduce ALPHA if
!       necessary to the value that minimizes the model.
!
   dgd = zero
   Do 200 i = 1, n
200 dgd = dgd + d (i) * dw (i)
   alpht = alpha
   If (dg+alpha*dgd > zero) Then
      alpha = - dg / dgd
   End If
!
!     Make a further reduction in ALPHA if necessary to preserve feasibility,
!       and put some scalar products of D with constraint gradients in W.
!
   alphm = alpha
   jsav = 0
   If (m > 0) Then
      Do 220 j = 1, m
         ad = zero
         If (resnew(j) > zero) Then
            Do 210 i = 1, n
210         ad = ad + amat (i, j) * d (i)
            If (alpha*ad > resnew(j)) Then
               alpha = resnew (j) / ad
               jsav = j
            End If
         End If
220   w (j) = ad
   End If
   alpha = max (alpha, alpbd)
   alpha = min (alpha, alphm)
   If (icount == nact) alpha = min (alpha, one)
!
!     Update STEP, G, RESNEW, RESACT and REDUCT.
!
   ss = zero
   Do 230 i = 1, n
      step (i) = step (i) + alpha * d (i)
      ss = ss + step (i) ** 2
230 g (i) = g (i) + alpha * dw (i)
   If (m > 0) Then
      Do 240 j = 1, m
         If (resnew(j) > zero) Then
            resnew (j) = max (resnew(j)-alpha*w(j), tiny)
         End If
240   Continue
   End If
   If (icount == nact .And. nact > 0) Then
      Do 250 k = 1, nact
250   resact (k) = (one-gamma) * resact (k)
   End If
   reduct = reduct - alpha * (dg+half*alpha*dgd)
!
!     Test for termination. Branch to label 40 if there is a new active
!       constraint and if the distance from STEP to the trust region
!       boundary is at least 0.2*SNORM.
!
   If (alpha == alpht) Go To 320
   temp = - alphm * (dg+half*alphm*dgd)
   If (temp <= ctest*reduct) Go To 320
   If (jsav > 0) Then
      If (ss <= 0.64_wp*snsq) Go To 40
      Go To 320
   End If
   If (icount == n) Go To 320
!
!     Calculate the next search direction, which is conjugate to the
!       previous one except in the case ICOUNT=NACT.
!
   If (nact > 0) Then
      Do 260 j = nact + 1, n
         w (j) = zero
         Do 260 i = 1, n
260   w (j) = w (j) + g (i) * qfac (i, j)
      Do 280 i = 1, n
         temp = zero
         Do 270 j = nact + 1, n
270      temp = temp + qfac (i, j) * w (j)
280   w (n+i) = temp
   Else
      Do 290 i = 1, n
290   w (n+i) = g (i)
   End If
   If (icount == nact) Then
      beta = zero
   Else
      wgd = zero
      Do 300 i = 1, n
300   wgd = wgd + w (n+i) * dw (i)
      beta = wgd / dgd
   End If
   Do 310 i = 1, n
310 d (i) = - w (n+i) + beta * d (i)
   alpbd = zero
   Go To 150
!
!     Return from the subroutine.
!
320 snorm = zero
   If (reduct > zero) snorm = sqrt (ss)
   g (1) = zero
   If (ncall > 1) g (1) = one
   Return
End Subroutine trstep

Subroutine update (n, npt, xpt, bmat, zmat, idz, ndim, sp, step, kopt, &
  knew, vlag, w)
   Implicit real(wp) (a-h, o-z)
   Dimension xpt (npt,*), bmat (ndim,*), zmat (npt,*), sp (*), step(*), &
   vlag (*), w (*)
!
!     The arguments N, NPT, XPT, BMAT, ZMAT, IDZ, NDIM ,SP and STEP are
!       identical to the corresponding arguments in SUBROUTINE LINCOB.
!     KOPT is such that XPT(KOPT,.) is the current trust region centre.
!     KNEW on exit is usually positive, and then it is the index of an
!       interpolation point to be moved to the position XPT(KOPT,.)+STEP(.).
!       It is set on entry either to its final value or to 0. In the latter
!       case, the final value of KNEW is chosen to maximize the denominator
!       of the matrix updating formula times a weighting factor.
!     VLAG and W are used for working space, the first NPT+N elements of
!       both of these vectors being required.
!
!     The arrays BMAT and ZMAT with IDZ are updated, the new matrices being
!       the ones that are suitable after the shift of the KNEW-th point to
!       the new position XPT(KOPT,.)+STEP(.). A return with KNEW set to zero
!       occurs if the calculation fails due to a zero denominator in the
!       updating formula, which should never happen.
!
!     Set some constants.
!
   half = 0.5_wp
   one = 1.0_wp
   zero = 0.0_wp
   nptm = npt - n - 1
!
!     Calculate VLAG and BETA for the current choice of STEP. The first NPT
!       elements of VLAG are set to the values of the Lagrange functions at
!       XPT(KOPT,.)+STEP(.). The first NPT components of W_check are held
!       in W, where W_check is defined in a paper on the updating method.
!
   Do 20 k = 1, npt
      w (k) = sp (npt+k) * (half*sp(npt+k)+sp(k))
      sum = zero
      Do 10 j = 1, n
10    sum = sum + bmat (k, j) * step (j)
20 vlag (k) = sum
   beta = zero
   Do 40 k = 1, nptm
      sum = zero
      Do 30 i = 1, npt
30    sum = sum + zmat (i, k) * w (i)
      If (k < idz) Then
         beta = beta + sum * sum
         sum = - sum
      Else
         beta = beta - sum * sum
      End If
      Do 40 i = 1, npt
40 vlag (i) = vlag (i) + sum * zmat (i, k)
   bsum = zero
   dx = zero
   ssq = zero
   Do 70 j = 1, n
      sum = zero
      Do 50 i = 1, npt
50    sum = sum + w (i) * bmat (i, j)
      bsum = bsum + sum * step (j)
      jp = npt + j
      Do 60 k = 1, n
60    sum = sum + bmat (jp, k) * step (k)
      vlag (jp) = sum
      bsum = bsum + sum * step (j)
      dx = dx + step (j) * xpt (kopt, j)
70 ssq = ssq + step (j) ** 2
   beta = dx * dx + ssq * (sp(kopt)+dx+dx+half*ssq) + beta - bsum
   vlag (kopt) = vlag (kopt) + one
!
!     If KNEW is zero initially, then pick the index of the interpolation
!       point to be deleted, by maximizing the absolute value of the
!       denominator of the updating formula times a weighting factor.
!
!
   If (knew == 0) Then
      denmax = zero
      Do 100 k = 1, npt
         hdiag = zero
         Do 80 j = 1, nptm
            temp = one
            If (j < idz) temp = - one
80       hdiag = hdiag + temp * zmat (k, j) ** 2
         denabs = abs (beta*hdiag+vlag(k)**2)
         distsq = zero
         Do 90 j = 1, n
90       distsq = distsq + (xpt(k, j)-xpt(kopt, j)) ** 2
         temp = denabs * distsq * distsq
         If (temp > denmax) Then
            denmax = temp
            knew = k
         End If
100   Continue
   End If
!
!     Apply the rotations that put zeros in the KNEW-th row of ZMAT.
!
   jl = 1
   If (nptm >= 2) Then
      Do 120 j = 2, nptm
         If (j == idz) Then
            jl = idz
         Else If (zmat(knew, j) /= zero) Then
            temp = sqrt (zmat(knew, jl)**2+zmat(knew, j)**2)
            tempa = zmat (knew, jl) / temp
            tempb = zmat (knew, j) / temp
            Do 110 i = 1, npt
               temp = tempa * zmat (i, jl) + tempb * zmat (i, j)
               zmat (i, j) = tempa * zmat (i, j) - tempb * zmat (i, jl)
110         zmat (i, jl) = temp
            zmat (knew, j) = zero
         End If
120   Continue
   End If
!
!     Put the first NPT components of the KNEW-th column of the Z Z^T matrix
!       into W, and calculate the parameters of the updating formula.
!
   tempa = zmat (knew, 1)
   If (idz >= 2) tempa = - tempa
   If (jl > 1) tempb = zmat (knew, jl)
   Do 130 i = 1, npt
      w (i) = tempa * zmat (i, 1)
      If (jl > 1) w (i) = w (i) + tempb * zmat (i, jl)
130 Continue
   alpha = w (knew)
   tau = vlag (knew)
   tausq = tau * tau
   denom = alpha * beta + tausq
   vlag (knew) = vlag (knew) - one
   If (denom == zero) Then
      knew = 0
      Go To 180
   End If
   sqrtdn = sqrt (abs(denom))
!
!     Complete the updating of ZMAT when there is only one nonzero element
!       in the KNEW-th row of the new matrix ZMAT. IFLAG is set to one when
!       the value of IDZ is going to be reduced.
!
   iflag = 0
   If (jl == 1) Then
      tempa = tau / sqrtdn
      tempb = zmat (knew, 1) / sqrtdn
      Do 140 i = 1, npt
140   zmat (i, 1) = tempa * zmat (i, 1) - tempb * vlag (i)
      If (denom < zero) Then
         If (idz == 1) Then
            idz = 2
         Else
            iflag = 1
         End If
      End If
   Else
!
!     Complete the updating of ZMAT in the alternative case.
!
      ja = 1
      If (beta >= zero) ja = jl
      jb = jl + 1 - ja
      temp = zmat (knew, jb) / denom
      tempa = temp * beta
      tempb = temp * tau
      temp = zmat (knew, ja)
      scala = one / sqrt (abs(beta)*temp*temp+tausq)
      scalb = scala * sqrtdn
      Do 150 i = 1, npt
         zmat (i, ja) = scala * (tau*zmat(i, ja)-temp*vlag(i))
150   zmat (i, jb) = scalb * (zmat(i, jb)-tempa*w(i)-tempb*vlag(i))
      If (denom <= zero) Then
         If (beta < zero) Then
            idz = idz + 1
         Else
            iflag = 1
         End If
      End If
   End If
!
!     Reduce IDZ when the diagonal part of the ZMAT times Diag(DZ) times
!       ZMAT^T factorization gains another positive element. Then exchange
!       the first and IDZ-th columns of ZMAT.
!
   If (iflag == 1) Then
      idz = idz - 1
      Do 160 i = 1, npt
         temp = zmat (i, 1)
         zmat (i, 1) = zmat (i, idz)
160   zmat (i, idz) = temp
   End If
!
!     Finally, update the matrix BMAT.
!
   Do 170 j = 1, n
      jp = npt + j
      w (jp) = bmat (knew, j)
      tempa = (alpha*vlag(jp)-tau*w(jp)) / denom
      tempb = (-beta*w(jp)-tau*vlag(jp)) / denom
      Do 170 i = 1, jp
         bmat (i, j) = bmat (i, j) + tempa * vlag (i) + tempb * w (i)
         If (i > npt) bmat (jp, i-npt) = bmat (i, j)
170 Continue
180 Return
End Subroutine update

subroutine lincoa_test()
!     Calculate the tetrahedron of least volume that encloses the points
!       (XP(J),YP(J),ZP(J)), J=1,2,...,NP. Our method requires the origin
!       to be strictly inside the convex hull of these points. There are
!       twelve variables that define the four faces of each tetrahedron
!       that is considered. Each face has the form ALPHA*X+BETA*Y+GAMMA*Z=1,
!       the variables X(3K-2), X(3K-1) and X(3K) being the values of ALPHA,
!       BETA and GAMMA for the K-th face, K=1,2,3,4. Let the set T contain
!       all points in three dimensions that can be reached from the origin
!       without crossing a face. Because the volume of T may be infinite,
!       the objective function is the smaller of FMAX and the volume of T,
!       where FMAX is set to an upper bound on the final volume initially.
!       There are 4*NP linear constraints on the variables, namely that each
!       of the given points (XP(J),YP(J),ZP(J)) shall be in T. Let XS = min
!       XP(J), YS = min YP(J), ZS = min ZP(J) and SS = max XP(J)+YP(J)+ZP(J),
!       where J runs from 1 to NP. The initial values of the variables are
!       X(1)=1/XS, X(5)=1/YS, X(9)=1/ZS, X(2)=X(3)=X(4)=X(6)=X(7) =X(8)=0
!       and X(10)=X(11)=X(12)=1/SS, which satisfy the linear constraints,
!       and which provide the bound FMAX=(SS-XS-YS-ZS)**3/6. Other details
!       of the test calculation are given below, including the choice of
!       the data points (XP(J),YP(J),ZP(J)), J=1,2,...,NP. The smaller final
!       value of the objective function in the case NPT=35 shows that the
!       problem has local minima.
!
    Implicit real(wp) (a-h, o-z)
    Common fmax
    Dimension xp (50), yp (50), zp (50), a (12, 200), b (200), x (12), w(500000)
    !
    !     Set some constants.
    !
    one = 1.0_wp
    two = 2.0_wp
    zero = 0.0_wp
    pi = 4.0_wp * atan (one)
    ia = 12
    n = 12
    !
    !     Set the data points.
    !
    np = 50
    sumx = zero
    sumy = zero
    sumz = zero
    Do 10 j = 1, np
       theta = dfloat (j-1) * pi / dfloat (np-1)
       xp (j) = cos (theta) * cos (two*theta)
       sumx = sumx + xp (j)
       yp (j) = sin (theta) * cos (two*theta)
       sumy = sumy + yp (j)
       zp (j) = sin (two*theta)
10  sumz = sumz + zp (j)
    sumx = sumx / dfloat (np)
    sumy = sumy / dfloat (np)
    sumz = sumz / dfloat (np)
    Do 20 j = 1, np
       xp (j) = xp (j) - sumx
       yp (j) = yp (j) - sumy
20  zp (j) = zp (j) - sumz
    !
    !     Set the linear constraints.
    !
    m = 4 * np
    Do 30 k = 1, m
       b (k) = one
       Do 30 i = 1, n
30  a (i, k) = zero
    Do 40 j = 1, np
       Do 40 i = 1, 4
          k = 4 * j + i - 4
          iw = 3 * i
          a (iw-2, k) = xp (j)
          a (iw-1, k) = yp (j)
40        a (iw, k) = zp (j)
    !
    !     Set the initial vector of variables. The JCASE=1,6 loop gives six
    !       different choices of NPT when LINCOA is called.
    !
    xs = zero
    ys = zero
    zs = zero
    ss = zero
    Do 50 j = 1, np
       xs = min (xs, xp(j))
       ys = min (ys, yp(j))
       zs = min (zs, zp(j))
50     ss = max (ss, xp(j)+yp(j)+zp(j))
    fmax = (ss-xs-ys-zs) ** 3 / 6.0_wp
    Do 80 jcase = 1, 6
       Do 60 i = 2, 8
60     x (i) = zero
       x (1) = one / xs
       x (5) = one / ys
       x (9) = one / zs
       x (10) = one / ss
       x (11) = one / ss
       x (12) = one / ss
    !
    !     Call of LINCOA, which provides the printing given at the end of this
    !       note.
    !
       npt = 5 * jcase + 10
       rhobeg = 1.0_wp
       rhoend = 1.0e-6_wp
       iprint = 1
       maxfun = 10000
       Print 70, npt, rhoend
70     Format (/ / 4 x, 'Output from LINCOA with  NPT =', i4,&
       '  and  RHOEND =', 1 pd12.4)
       Call lincoa (n, npt, m, a, ia, b, x, rhobeg, rhoend, iprint, maxfun,w,calfun)
80     Continue

contains

    Subroutine calfun (n, x, f)
      Implicit real(wp) (a-h, o-z)
      !Common fmax
      Dimension x (*)
      zero = 0.0_wp
      f = fmax
      v12 = x (1) * x (5) - x (4) * x (2)
      v13 = x (1) * x (8) - x (7) * x (2)
      v14 = x (1) * x (11) - x (10) * x (2)
      v23 = x (4) * x (8) - x (7) * x (5)
      v24 = x (4) * x (11) - x (10) * x (5)
      v34 = x (7) * x (11) - x (10) * x (8)
      del1 = v23 * x (12) - v24 * x (9) + v34 * x (6)
      If (del1 <= zero) return
      del2 = - v34 * x (3) - v13 * x (12) + v14 * x (9)
      If (del2 <= zero) return
      del3 = - v14 * x (6) + v24 * x (3) + v12 * x (12)
      If (del3 <= zero) return
      del4 = - v12 * x (9) + v13 * x (6) - v23 * x (3)
      If (del4 <= zero) return
      temp = (del1+del2+del3+del4) ** 3 / (del1*del2*del3*del4)
      f = min (temp/6.0_wp, fmax)
    End Subroutine calfun

end subroutine lincoa_test

end module lincoa_module