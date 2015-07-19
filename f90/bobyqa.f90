    module bobyqa_module

    private

    abstract interface
        subroutine func (N,X,F)  !! calfun interface
        implicit none
        integer :: n
        real * 8 :: x(*)
        real * 8 :: f
        end subroutine func
    end interface

    public :: bobyqa
    public :: bobyqa_test

    contains

Subroutine bobyqa (n, npt, x, xl, xu, rhobeg, rhoend, iprint, maxfun, w, calfun)
      Implicit real * 8 (a-h, o-z)
      Dimension x (*), xl (*), xu (*), w (*)
      procedure(func) :: calfun
!
!     This subroutine seeks the least value of a function of many variables,
!     by applying a trust region method that forms quadratic models by
!     interpolation. There is usually some freedom in the interpolation
!     conditions, which is taken up by minimizing the Frobenius norm of
!     the change to the second derivative of the model, beginning with the
!     zero matrix. The values of the variables are constrained by upper and
!     lower bounds. The arguments of the subroutine are as follows.
!
!     N must be set to the number of variables and must be at least two.
!     NPT is the number of interpolation conditions. Its value must be in
!       the interval [N+2,(N+1)(N+2)/2]. Choices that exceed 2*N+1 are not
!       recommended.
!     Initial values of the variables must be set in X(1),X(2),...,X(N). They
!       will be changed to the values that give the least calculated F.
!     For I=1,2,...,N, XL(I) and XU(I) must provide the lower and upper
!       bounds, respectively, on X(I). The construction of quadratic models
!       requires XL(I) to be strictly less than XU(I) for each I. Further,
!       the contribution to a model from changes to the I-th variable is
!       damaged severely by rounding errors if XU(I)-XL(I) is too small.
!     RHOBEG and RHOEND must be set to the initial and final values of a trust
!       region radius, so both must be positive with RHOEND no greater than
!       RHOBEG. Typically, RHOBEG should be about one tenth of the greatest
!       expected change to a variable, while RHOEND should indicate the
!       accuracy that is required in the final values of the variables. An
!       error return occurs if any of the differences XU(I)-XL(I), I=1,...,N,
!       is less than 2*RHOBEG.
!     The value of IPRINT should be set to 0, 1, 2 or 3, which controls the
!       amount of printing. Specifically, there is no output if IPRINT=0 and
!       there is output only at the return if IPRINT=1. Otherwise, each new
!       value of RHO is printed, with the best vector of variables so far and
!       the corresponding value of the objective function. Further, each new
!       value of F with its variables are output if IPRINT=3.
!     MAXFUN must be set to an upper bound on the number of calls of CALFUN.
!     The array W will be used for working space. Its length must be at least
!       (NPT+5)*(NPT+N)+3*N*(N+5)/2.
!
!     SUBROUTINE CALFUN (N,X,F) has to be provided by the user. It must set
!     F to the value of the objective function for the current values of the
!     variables X(1),X(2),...,X(N), which are generated automatically in a
!     way that satisfies the bounds given in XL and XU.
!
!     Return if the value of NPT is unacceptable.
!
      np = n + 1
      If (npt < n+2 .Or. npt > ((n+2)*np)/2) Then
         Print 10
10       Format (/ 4 x, 'Return from BOBYQA because NPT is not in',&
         ' the required interval')
         Go To 40
      End If
!
!     Partition the working space array, so that different parts of it can
!     be treated separately during the calculation of BOBYQB. The partition
!     requires the first (NPT+2)*(NPT+N)+3*N*(N+5)/2 elements of W plus the
!     space that is taken by the last array in the argument list of BOBYQB.
!
      ndim = npt + n
      ixb = 1
      ixp = ixb + n
      ifv = ixp + n * npt
      ixo = ifv + npt
      igo = ixo + n
      ihq = igo + n
      ipq = ihq + (n*np) / 2
      ibmat = ipq + npt
      izmat = ibmat + ndim * n
      isl = izmat + npt * (npt-np)
      isu = isl + n
      ixn = isu + n
      ixa = ixn + n
      id = ixa + n
      ivl = id + n
      iw = ivl + ndim
!
!     Return if there is insufficient space between the bounds. Modify the
!     initial X if necessary in order to avoid conflicts between the bounds
!     and the construction of the first quadratic model. The lower and upper
!     bounds on moves from the updated X are set now, in the ISL and ISU
!     partitions of W, in order to provide useful and exact information about
!     components of X that become within distance RHOBEG from their bounds.
!
      zero = 0.0d0
      Do 30 j = 1, n
         temp = xu (j) - xl (j)
         If (temp < rhobeg+rhobeg) Then
            Print 20
20          Format (/ 4 x, 'Return from BOBYQA because one of the',&
            ' differences XU(I)-XL(I)' / 6 x, ' is less than 2*RHOBEG.')
            Go To 40
         End If
         jsl = isl + j - 1
         jsu = jsl + n
         w (jsl) = xl (j) - x (j)
         w (jsu) = xu (j) - x (j)
         If (w(jsl) >=-rhobeg) Then
            If (w(jsl) >= zero) Then
               x (j) = xl (j)
               w (jsl) = zero
               w (jsu) = temp
            Else
               x (j) = xl (j) + rhobeg
               w (jsl) = - rhobeg
               w (jsu) = dmax1 (xu(j)-x(j), rhobeg)
            End If
         Else If (w(jsu) <= rhobeg) Then
            If (w(jsu) <= zero) Then
               x (j) = xu (j)
               w (jsl) = - temp
               w (jsu) = zero
            Else
               x (j) = xu (j) - rhobeg
               w (jsl) = dmin1 (xl(j)-x(j),-rhobeg)
               w (jsu) = rhobeg
            End If
         End If
30    Continue
!
!     Make the call of BOBYQB.
!
      Call bobyqb (n, npt, x, xl, xu, rhobeg, rhoend, iprint, maxfun, &
     & w(ixb), w(ixp), w(ifv), w(ixo), w(igo), w(ihq), w(ipq), &
     & w(ibmat), w(izmat), ndim, w(isl), w(isu), w(ixn), w(ixa), w(id), &
     & w(ivl), w(iw), calfun)
40    Return
End Subroutine bobyqa

Subroutine bobyqb (n, npt, x, xl, xu, rhobeg, rhoend, iprint, maxfun, &
 xbase, xpt, fval, xopt, gopt, hq, pq, bmat, zmat, ndim, sl, su, xnew, &
 xalt, d, vlag, w, calfun)
      Implicit real * 8 (a-h, o-z)
      Dimension x (*), xl (*), xu (*), xbase (*), xpt (npt,*), fval &
     & (*), xopt (*), gopt (*), hq (*), pq (*), bmat (ndim,*), zmat &
     & (npt,*), sl (*), su (*), xnew (*), xalt (*), d (*), vlag (*), w &
     & (*)
    procedure(func) :: calfun
!
!     The arguments N, NPT, X, XL, XU, RHOBEG, RHOEND, IPRINT and MAXFUN
!       are identical to the corresponding arguments in SUBROUTINE BOBYQA.
!     XBASE holds a shift of origin that should reduce the contributions
!       from rounding errors to values of the model and Lagrange functions.
!     XPT is a two-dimensional array that holds the coordinates of the
!       interpolation points relative to XBASE.
!     FVAL holds the values of F at the interpolation points.
!     XOPT is set to the displacement from XBASE of the trust region centre.
!     GOPT holds the gradient of the quadratic model at XBASE+XOPT.
!     HQ holds the explicit second derivatives of the quadratic model.
!     PQ contains the parameters of the implicit second derivatives of the
!       quadratic model.
!     BMAT holds the last N columns of H.
!     ZMAT holds the factorization of the leading NPT by NPT submatrix of H,
!       this factorization being ZMAT times ZMAT^T, which provides both the
!       correct rank and positive semi-definiteness.
!     NDIM is the first dimension of BMAT and has the value NPT+N.
!     SL and SU hold the differences XL-XBASE and XU-XBASE, respectively.
!       All the components of every XOPT are going to satisfy the bounds
!       SL(I) .LEQ. XOPT(I) .LEQ. SU(I), with appropriate equalities when
!       XOPT is on a constraint boundary.
!     XNEW is chosen by SUBROUTINE TRSBOX or ALTMOV. Usually XBASE+XNEW is the
!       vector of variables for the next call of CALFUN. XNEW also satisfies
!       the SL and SU constraints in the way that has just been mentioned.
!     XALT is an alternative to XNEW, chosen by ALTMOV, that may replace XNEW
!       in order to increase the denominator in the updating of UPDATE.
!     D is reserved for a trial step from XOPT, which is usually XNEW-XOPT.
!     VLAG contains the values of the Lagrange functions at a new point X.
!       They are part of a product that requires VLAG to be of length NDIM.
!     W is a one-dimensional array that is used for working space. Its length
!       must be at least 3*NDIM = 3*(NPT+N).
!
!     Set some constants.
!
      half = 0.5d0
      one = 1.0d0
      ten = 10.0d0
      tenth = 0.1d0
      two = 2.0d0
      zero = 0.0d0
      np = n + 1
      nptm = npt - np
      nh = (n*np) / 2
!
!     The call of PRELIM sets the elements of XBASE, XPT, FVAL, GOPT, HQ, PQ,
!     BMAT and ZMAT for the first iteration, with the corresponding values of
!     of NF and KOPT, which are the number of calls of CALFUN so far and the
!     index of the interpolation point at the trust region centre. Then the
!     initial XOPT is set too. The branch to label 720 occurs if MAXFUN is
!     less than NPT. GOPT will be updated if KOPT is different from KBASE.
!
      Call prelim (n, npt, x, xl, xu, rhobeg, iprint, maxfun, xbase, &
     & xpt, fval, gopt, hq, pq, bmat, zmat, ndim, sl, su, nf, kopt, calfun)
      xoptsq = zero
      Do 10 i = 1, n
         xopt (i) = xpt (kopt, i)
10    xoptsq = xoptsq + xopt (i) ** 2
      fsave = fval (1)
      If (nf < npt) Then
         If (iprint > 0) Print 390
         Go To 720
      End If
      kbase = 1
!
!     Complete the settings that are required for the iterative procedure.
!
      rho = rhobeg
      delta = rho
      nresc = nf
      ntrits = 0
      diffa = zero
      diffb = zero
      itest = 0
      nfsav = nf
!
!     Update GOPT if necessary before the first iteration and after each
!     call of RESCUE that makes a call of CALFUN.
!
20    If (kopt /= kbase) Then
         ih = 0
         Do 30 j = 1, n
            Do 30 i = 1, j
               ih = ih + 1
               If (i < j) gopt (j) = gopt (j) + hq (ih) * xopt (i)
30       gopt (i) = gopt (i) + hq (ih) * xopt (j)
         If (nf > npt) Then
            Do 50 k = 1, npt
               temp = zero
               Do 40 j = 1, n
40             temp = temp + xpt (k, j) * xopt (j)
               temp = pq (k) * temp
               Do 50 i = 1, n
50          gopt (i) = gopt (i) + temp * xpt (k, i)
         End If
      End If
!
!     Generate the next point in the trust region that provides a small value
!     of the quadratic model subject to the constraints on the variables.
!     The integer NTRITS is set to the number "trust region" iterations that
!     have occurred since the last "alternative" iteration. If the length
!     of XNEW-XOPT is less than HALF*RHO, however, then there is a branch to
!     label 650 or 680 with NTRITS=-1, instead of calculating F at XNEW.
!
60    Call trsbox (n, npt, xpt, xopt, gopt, hq, pq, sl, su, delta, &
     & xnew, d, w, w(np), w(np+n), w(np+2*n), w(np+3*n), dsq, crvmin)
      dnorm = dmin1 (delta, dsqrt(dsq))
      If (dnorm < half*rho) Then
         ntrits = - 1
         distsq = (ten*rho) ** 2
         If (nf <= nfsav+2) Go To 650
!
!     The following choice between labels 650 and 680 depends on whether or
!     not our work with the current RHO seems to be complete. Either RHO is
!     decreased or termination occurs if the errors in the quadratic model at
!     the last three interpolation points compare favourably with predictions
!     of likely improvements to the model within distance HALF*RHO of XOPT.
!
         errbig = dmax1 (diffa, diffb, diffc)
         frhosq = 0.125d0 * rho * rho
         If (crvmin > zero .And. errbig > frhosq*crvmin) Go To &
        & 650
         bdtol = errbig / rho
         Do 80 j = 1, n
            bdtest = bdtol
            If (xnew(j) == sl(j)) bdtest = w (j)
            If (xnew(j) == su(j)) bdtest = - w (j)
            If (bdtest < bdtol) Then
               curv = hq ((j+j*j)/2)
               Do 70 k = 1, npt
70             curv = curv + pq (k) * xpt (k, j) ** 2
               bdtest = bdtest + half * curv * rho
               If (bdtest < bdtol) Go To 650
            End If
80       Continue
         Go To 680
      End If
      ntrits = ntrits + 1
!
!     Severe cancellation is likely to occur if XOPT is too far from XBASE.
!     If the following test holds, then XBASE is shifted so that XOPT becomes
!     zero. The appropriate changes are made to BMAT and to the second
!     derivatives of the current model, beginning with the changes to BMAT
!     that do not depend on ZMAT. VLAG is used temporarily for working space.
!
90    If (dsq <= 1.0d-3*xoptsq) Then
         fracsq = 0.25d0 * xoptsq
         sumpq = zero
         Do 110 k = 1, npt
            sumpq = sumpq + pq (k)
            sum = - half * xoptsq
            Do 100 i = 1, n
100         sum = sum + xpt (k, i) * xopt (i)
            w (npt+k) = sum
            temp = fracsq - half * sum
            Do 110 i = 1, n
               w (i) = bmat (k, i)
               vlag (i) = sum * xpt (k, i) + temp * xopt (i)
               ip = npt + i
               Do 110 j = 1, i
110      bmat (ip, j) = bmat (ip, j) + w (i) * vlag (j) + vlag (i) * w &
        & (j)
!
!     Then the revisions of BMAT that depend on ZMAT are calculated.
!
         Do 150 jj = 1, nptm
            sumz = zero
            sumw = zero
            Do 120 k = 1, npt
               sumz = sumz + zmat (k, jj)
               vlag (k) = w (npt+k) * zmat (k, jj)
120         sumw = sumw + vlag (k)
            Do 140 j = 1, n
               sum = (fracsq*sumz-half*sumw) * xopt (j)
               Do 130 k = 1, npt
130            sum = sum + vlag (k) * xpt (k, j)
               w (j) = sum
               Do 140 k = 1, npt
140         bmat (k, j) = bmat (k, j) + sum * zmat (k, jj)
            Do 150 i = 1, n
               ip = i + npt
               temp = w (i)
               Do 150 j = 1, i
150      bmat (ip, j) = bmat (ip, j) + temp * w (j)
!
!     The following instructions complete the shift, including the changes
!     to the second derivative parameters of the quadratic model.
!
         ih = 0
         Do 170 j = 1, n
            w (j) = - half * sumpq * xopt (j)
            Do 160 k = 1, npt
               w (j) = w (j) + pq (k) * xpt (k, j)
160         xpt (k, j) = xpt (k, j) - xopt (j)
            Do 170 i = 1, j
               ih = ih + 1
               hq (ih) = hq (ih) + w (i) * xopt (j) + xopt (i) * w (j)
170      bmat (npt+i, j) = bmat (npt+j, i)
         Do 180 i = 1, n
            xbase (i) = xbase (i) + xopt (i)
            xnew (i) = xnew (i) - xopt (i)
            sl (i) = sl (i) - xopt (i)
            su (i) = su (i) - xopt (i)
180      xopt (i) = zero
         xoptsq = zero
      End If
      If (ntrits == 0) Go To 210
      Go To 230
!
!     XBASE is also moved to XOPT by a call of RESCUE. This calculation is
!     more expensive than the previous shift, because new matrices BMAT and
!     ZMAT are generated from scratch, which may include the replacement of
!     interpolation points whose positions seem to be causing near linear
!     dependence in the interpolation conditions. Therefore RESCUE is called
!     only if rounding errors have reduced by at least a factor of two the
!     denominator of the formula for updating the H matrix. It provides a
!     useful safeguard, but is not invoked in most applications of BOBYQA.
!
190   nfsav = nf
      kbase = kopt
      Call rescue (n, npt, xl, xu, iprint, maxfun, xbase, xpt, fval, &
     & xopt, gopt, hq, pq, bmat, zmat, ndim, sl, su, nf, delta, kopt, &
     & vlag, w, w(n+np), w(ndim+np), calfun)
!
!     XOPT is updated now in case the branch below to label 720 is taken.
!     Any updating of GOPT occurs after the branch below to label 20, which
!     leads to a trust region iteration as does the branch to label 60.
!
      xoptsq = zero
      If (kopt /= kbase) Then
         Do 200 i = 1, n
            xopt (i) = xpt (kopt, i)
200      xoptsq = xoptsq + xopt (i) ** 2
      End If
      If (nf < 0) Then
         nf = maxfun
         If (iprint > 0) Print 390
         Go To 720
      End If
      nresc = nf
      If (nfsav < nf) Then
         nfsav = nf
         Go To 20
      End If
      If (ntrits > 0) Go To 60
!
!     Pick two alternative vectors of variables, relative to XBASE, that
!     are suitable as new positions of the KNEW-th interpolation point.
!     Firstly, XNEW is set to the point on a line through XOPT and another
!     interpolation point that minimizes the predicted value of the next
!     denominator, subject to ||XNEW - XOPT|| .LEQ. ADELT and to the SL
!     and SU bounds. Secondly, XALT is set to the best feasible point on
!     a constrained version of the Cauchy step of the KNEW-th Lagrange
!     function, the corresponding value of the square of this function
!     being returned in CAUCHY. The choice between these alternatives is
!     going to be made when the denominator is calculated.
!
210   Call altmov (n, npt, xpt, xopt, bmat, zmat, ndim, sl, su, kopt, &
     & knew, adelt, xnew, xalt, alpha, cauchy, w, w(np), w(ndim+1))
      Do 220 i = 1, n
220   d (i) = xnew (i) - xopt (i)
!
!     Calculate VLAG and BETA for the current choice of D. The scalar
!     product of D with XPT(K,.) is going to be held in W(NPT+K) for
!     use when VQUAD is calculated.
!
230   Do 250 k = 1, npt
         suma = zero
         sumb = zero
         sum = zero
         Do 240 j = 1, n
            suma = suma + xpt (k, j) * d (j)
            sumb = sumb + xpt (k, j) * xopt (j)
240      sum = sum + bmat (k, j) * d (j)
         w (k) = suma * (half*suma+sumb)
         vlag (k) = sum
250   w (npt+k) = suma
      beta = zero
      Do 270 jj = 1, nptm
         sum = zero
         Do 260 k = 1, npt
260      sum = sum + zmat (k, jj) * w (k)
         beta = beta - sum * sum
         Do 270 k = 1, npt
270   vlag (k) = vlag (k) + sum * zmat (k, jj)
      dsq = zero
      bsum = zero
      dx = zero
      Do 300 j = 1, n
         dsq = dsq + d (j) ** 2
         sum = zero
         Do 280 k = 1, npt
280      sum = sum + w (k) * bmat (k, j)
         bsum = bsum + sum * d (j)
         jp = npt + j
         Do 290 i = 1, n
290      sum = sum + bmat (jp, i) * d (i)
         vlag (jp) = sum
         bsum = bsum + sum * d (j)
300   dx = dx + d (j) * xopt (j)
      beta = dx * dx + dsq * (xoptsq+dx+dx+half*dsq) + beta - bsum
      vlag (kopt) = vlag (kopt) + one
!
!     If NTRITS is zero, the denominator may be increased by replacing
!     the step D of ALTMOV by a Cauchy step. Then RESCUE may be called if
!     rounding errors have damaged the chosen denominator.
!
      If (ntrits == 0) Then
         denom = vlag (knew) ** 2 + alpha * beta
         If (denom < cauchy .And. cauchy > zero) Then
            Do 310 i = 1, n
               xnew (i) = xalt (i)
310         d (i) = xnew (i) - xopt (i)
            cauchy = zero
            Go To 230
         End If
         If (denom <= half*vlag(knew)**2) Then
            If (nf > nresc) Go To 190
            If (iprint > 0) Print 320
320         Format (/ 5 x, 'Return from BOBYQA because of much',&
            ' cancellation in a denominator.')
            Go To 720
         End If
!
!     Alternatively, if NTRITS is positive, then set KNEW to the index of
!     the next interpolation point to be deleted to make room for a trust
!     region step. Again RESCUE may be called if rounding errors have damaged
!     the chosen denominator, which is the reason for attempting to select
!     KNEW before calculating the next value of the objective function.
!
      Else
         delsq = delta * delta
         scaden = zero
         biglsq = zero
         knew = 0
         Do 350 k = 1, npt
            If (k == kopt) Go To 350
            hdiag = zero
            Do 330 jj = 1, nptm
330         hdiag = hdiag + zmat (k, jj) ** 2
            den = beta * hdiag + vlag (k) ** 2
            distsq = zero
            Do 340 j = 1, n
340         distsq = distsq + (xpt(k, j)-xopt(j)) ** 2
            temp = dmax1 (one, (distsq/delsq)**2)
            If (temp*den > scaden) Then
               scaden = temp * den
               knew = k
               denom = den
            End If
            biglsq = dmax1 (biglsq, temp*vlag(k)**2)
350      Continue
         If (scaden <= half*biglsq) Then
            If (nf > nresc) Go To 190
            If (iprint > 0) Print 320
            Go To 720
         End If
      End If
!
!     Put the variables for the next calculation of the objective function
!       in XNEW, with any adjustments for the bounds.
!
!
!     Calculate the value of the objective function at XBASE+XNEW, unless
!       the limit on the number of calculations of F has been reached.
!
360   Do 380 i = 1, n
         x (i) = dmin1 (dmax1(xl(i), xbase(i)+xnew(i)), xu(i))
         If (xnew(i) == sl(i)) x (i) = xl (i)
         If (xnew(i) == su(i)) x (i) = xu (i)
380   Continue
      If (nf >= maxfun) Then
         If (iprint > 0) Print 390
390      Format (/ 4 x, 'Return from BOBYQA because CALFUN has been',&
         ' called MAXFUN times.')
         Go To 720
      End If
      nf = nf + 1
      Call calfun (n, x, f)
      If (iprint == 3) Then
         Print 400, nf, f, (x(i), i=1, n)
400      Format (/ 4 x, 'Function number', i6, '    F =', 1 pd18.10,&
         '   The corresponding X is:' / (2 x, 5d15.6))
      End If
      If (ntrits ==-1) Then
         fsave = f
         Go To 720
      End If
!
!     Use the quadratic model to predict the change in F due to the step D,
!       and set DIFF to the error of this prediction.
!
      fopt = fval (kopt)
      vquad = zero
      ih = 0
      Do 410 j = 1, n
         vquad = vquad + d (j) * gopt (j)
         Do 410 i = 1, j
            ih = ih + 1
            temp = d (i) * d (j)
            If (i == j) temp = half * temp
410   vquad = vquad + hq (ih) * temp
      Do 420 k = 1, npt
420   vquad = vquad + half * pq (k) * w (npt+k) ** 2
      diff = f - fopt - vquad
      diffc = diffb
      diffb = diffa
      diffa = dabs (diff)
      If (dnorm > rho) nfsav = nf
!
!     Pick the next value of DELTA after a trust region step.
!
      If (ntrits > 0) Then
         If (vquad >= zero) Then
            If (iprint > 0) Print 430
430         Format (/ 4 x, 'Return from BOBYQA because a trust',&
            ' region step has failed to reduce Q.')
            Go To 720
         End If
         ratio = (f-fopt) / vquad
         If (ratio <= tenth) Then
            delta = dmin1 (half*delta, dnorm)
         Else If (ratio <= 0.7d0) Then
            delta = dmax1 (half*delta, dnorm)
         Else
            delta = dmax1 (half*delta, dnorm+dnorm)
         End If
         If (delta <= 1.5d0*rho) delta = rho
!
!     Recalculate KNEW and DENOM if the new F is less than FOPT.
!
         If (f < fopt) Then
            ksav = knew
            densav = denom
            delsq = delta * delta
            scaden = zero
            biglsq = zero
            knew = 0
            Do 460 k = 1, npt
               hdiag = zero
               Do 440 jj = 1, nptm
440            hdiag = hdiag + zmat (k, jj) ** 2
               den = beta * hdiag + vlag (k) ** 2
               distsq = zero
               Do 450 j = 1, n
450            distsq = distsq + (xpt(k, j)-xnew(j)) ** 2
               temp = dmax1 (one, (distsq/delsq)**2)
               If (temp*den > scaden) Then
                  scaden = temp * den
                  knew = k
                  denom = den
               End If
460         biglsq = dmax1 (biglsq, temp*vlag(k)**2)
            If (scaden <= half*biglsq) Then
               knew = ksav
               denom = densav
            End If
         End If
      End If
!
!     Update BMAT and ZMAT, so that the KNEW-th interpolation point can be
!     moved. Also update the second derivative terms of the model.
!
      Call update (n, npt, bmat, zmat, ndim, vlag, beta, denom, knew, &
     & w)
      ih = 0
      pqold = pq (knew)
      pq (knew) = zero
      Do 470 i = 1, n
         temp = pqold * xpt (knew, i)
         Do 470 j = 1, i
            ih = ih + 1
470   hq (ih) = hq (ih) + temp * xpt (knew, j)
      Do 480 jj = 1, nptm
         temp = diff * zmat (knew, jj)
         Do 480 k = 1, npt
480   pq (k) = pq (k) + temp * zmat (k, jj)
!
!     Include the new interpolation point, and make the changes to GOPT at
!     the old XOPT that are caused by the updating of the quadratic model.
!
      fval (knew) = f
      Do 490 i = 1, n
         xpt (knew, i) = xnew (i)
490   w (i) = bmat (knew, i)
      Do 520 k = 1, npt
         suma = zero
         Do 500 jj = 1, nptm
500      suma = suma + zmat (knew, jj) * zmat (k, jj)
         sumb = zero
         Do 510 j = 1, n
510      sumb = sumb + xpt (k, j) * xopt (j)
         temp = suma * sumb
         Do 520 i = 1, n
520   w (i) = w (i) + temp * xpt (k, i)
      Do 530 i = 1, n
530   gopt (i) = gopt (i) + diff * w (i)
!
!     Update XOPT, GOPT and KOPT if the new calculated F is less than FOPT.
!
      If (f < fopt) Then
         kopt = knew
         xoptsq = zero
         ih = 0
         Do 540 j = 1, n
            xopt (j) = xnew (j)
            xoptsq = xoptsq + xopt (j) ** 2
            Do 540 i = 1, j
               ih = ih + 1
               If (i < j) gopt (j) = gopt (j) + hq (ih) * d (i)
540      gopt (i) = gopt (i) + hq (ih) * d (j)
         Do 560 k = 1, npt
            temp = zero
            Do 550 j = 1, n
550         temp = temp + xpt (k, j) * d (j)
            temp = pq (k) * temp
            Do 560 i = 1, n
560      gopt (i) = gopt (i) + temp * xpt (k, i)
      End If
!
!     Calculate the parameters of the least Frobenius norm interpolant to
!     the current data, the gradient of this interpolant at XOPT being put
!     into VLAG(NPT+I), I=1,2,...,N.
!
      If (ntrits > 0) Then
         Do 570 k = 1, npt
            vlag (k) = fval (k) - fval (kopt)
570      w (k) = zero
         Do 590 j = 1, nptm
            sum = zero
            Do 580 k = 1, npt
580         sum = sum + zmat (k, j) * vlag (k)
            Do 590 k = 1, npt
590      w (k) = w (k) + sum * zmat (k, j)
         Do 610 k = 1, npt
            sum = zero
            Do 600 j = 1, n
600         sum = sum + xpt (k, j) * xopt (j)
            w (k+npt) = w (k)
610      w (k) = sum * w (k)
         gqsq = zero
         gisq = zero
         Do 630 i = 1, n
            sum = zero
            Do 620 k = 1, npt
620         sum = sum + bmat (k, i) * vlag (k) + xpt (k, i) * w (k)
            If (xopt(i) == sl(i)) Then
               gqsq = gqsq + dmin1 (zero, gopt(i)) ** 2
               gisq = gisq + dmin1 (zero, sum) ** 2
            Else If (xopt(i) == su(i)) Then
               gqsq = gqsq + dmax1 (zero, gopt(i)) ** 2
               gisq = gisq + dmax1 (zero, sum) ** 2
            Else
               gqsq = gqsq + gopt (i) ** 2
               gisq = gisq + sum * sum
            End If
630      vlag (npt+i) = sum
!
!     Test whether to replace the new quadratic model by the least Frobenius
!     norm interpolant, making the replacement if the test is satisfied.
!
         itest = itest + 1
         If (gqsq < ten*gisq) itest = 0
         If (itest >= 3) Then
            Do 640 i = 1, max0 (npt, nh)
               If (i <= n) gopt (i) = vlag (npt+i)
               If (i <= npt) pq (i) = w (npt+i)
               If (i <= nh) hq (i) = zero
               itest = 0
640         Continue
         End If
      End If
!
!     If a trust region step has provided a sufficient decrease in F, then
!     branch for another trust region calculation. The case NTRITS=0 occurs
!     when the new interpolation point was reached by an alternative step.
!
      If (ntrits == 0) Go To 60
      If (f <= fopt+tenth*vquad) Go To 60
!
!     Alternatively, find out if the interpolation points are close enough
!       to the best point so far.
!
      distsq = dmax1 ((two*delta)**2, (ten*rho)**2)
650   knew = 0
      Do 670 k = 1, npt
         sum = zero
         Do 660 j = 1, n
660      sum = sum + (xpt(k, j)-xopt(j)) ** 2
         If (sum > distsq) Then
            knew = k
            distsq = sum
         End If
670   Continue
!
!     If KNEW is positive, then ALTMOV finds alternative new positions for
!     the KNEW-th interpolation point within distance ADELT of XOPT. It is
!     reached via label 90. Otherwise, there is a branch to label 60 for
!     another trust region iteration, unless the calculations with the
!     current RHO are complete.
!
      If (knew > 0) Then
         dist = dsqrt (distsq)
         If (ntrits ==-1) Then
            delta = dmin1 (tenth*delta, half*dist)
            If (delta <= 1.5d0*rho) delta = rho
         End If
         ntrits = 0
         adelt = dmax1 (dmin1(tenth*dist, delta), rho)
         dsq = adelt * adelt
         Go To 90
      End If
      If (ntrits ==-1) Go To 680
      If (ratio > zero) Go To 60
      If (dmax1(delta, dnorm) > rho) Go To 60
!
!     The calculations with the current value of RHO are complete. Pick the
!       next values of RHO and DELTA.
!
680   If (rho > rhoend) Then
         delta = half * rho
         ratio = rho / rhoend
         If (ratio <= 16.0d0) Then
            rho = rhoend
         Else If (ratio <= 250.0d0) Then
            rho = dsqrt (ratio) * rhoend
         Else
            rho = tenth * rho
         End If
         delta = dmax1 (delta, rho)
         If (iprint >= 2) Then
            If (iprint >= 3) Print 690
690         Format (5 x)
            Print 700, rho, nf
700         Format (/ 4 x, 'New RHO =', 1 pd11.4, 5 x, 'Number of',&
            ' function values =', i6)
            Print 710, fval (kopt), (xbase(i)+xopt(i), i=1, n)
710         Format (4 x, 'Least value of F =', 1 pd23.15, 9 x,&
            'The corresponding X is:'/(2 x, 5d15.6))
         End If
         ntrits = 0
         nfsav = nf
         Go To 60
      End If
!
!     Return from the calculation, after another Newton-Raphson step, if
!       it is too short to have been tried before.
!
      If (ntrits ==-1) Go To 360
720   If (fval(kopt) <= fsave) Then
         Do 730 i = 1, n
            x (i) = dmin1 (dmax1(xl(i), xbase(i)+xopt(i)), xu(i))
            If (xopt(i) == sl(i)) x (i) = xl (i)
            If (xopt(i) == su(i)) x (i) = xu (i)
730      Continue
         f = fval (kopt)
      End If
      If (iprint >= 1) Then
         Print 740, nf
740      Format (/ 4 x, 'At the return from BOBYQA', 5 x,&
        'Number of function values =', i6)
         Print 710, f, (x(i), i=1, n)
      End If
      Return
End Subroutine bobyqb

Subroutine altmov (n, npt, xpt, xopt, bmat, zmat, ndim, sl, su, kopt, &
    knew, adelt, xnew, xalt, alpha, cauchy, glag, hcol, w)

      Implicit real * 8 (a-h, o-z)
      Dimension xpt (npt,*), xopt (*), bmat (ndim,*), zmat (npt,*), sl &
     & (*), su (*), xnew (*), xalt (*), glag (*), hcol (*), w (*)
!
!     The arguments N, NPT, XPT, XOPT, BMAT, ZMAT, NDIM, SL and SU all have
!       the same meanings as the corresponding arguments of BOBYQB.
!     KOPT is the index of the optimal interpolation point.
!     KNEW is the index of the interpolation point that is going to be moved.
!     ADELT is the current trust region bound.
!     XNEW will be set to a suitable new position for the interpolation point
!       XPT(KNEW,.). Specifically, it satisfies the SL, SU and trust region
!       bounds and it should provide a large denominator in the next call of
!       UPDATE. The step XNEW-XOPT from XOPT is restricted to moves along the
!       straight lines through XOPT and another interpolation point.
!     XALT also provides a large value of the modulus of the KNEW-th Lagrange
!       function subject to the constraints that have been mentioned, its main
!       difference from XNEW being that XALT-XOPT is a constrained version of
!       the Cauchy step within the trust region. An exception is that XALT is
!       not calculated if all components of GLAG (see below) are zero.
!     ALPHA will be set to the KNEW-th diagonal element of the H matrix.
!     CAUCHY will be set to the square of the KNEW-th Lagrange function at
!       the step XALT-XOPT from XOPT for the vector XALT that is returned,
!       except that CAUCHY is set to zero if XALT is not calculated.
!     GLAG is a working space vector of length N for the gradient of the
!       KNEW-th Lagrange function at XOPT.
!     HCOL is a working space vector of length NPT for the second derivative
!       coefficients of the KNEW-th Lagrange function.
!     W is a working space vector of length 2N that is going to hold the
!       constrained Cauchy step from XOPT of the Lagrange function, followed
!       by the downhill version of XALT when the uphill step is calculated.
!
!     Set the first NPT components of W to the leading elements of the
!     KNEW-th column of the H matrix.
!
      half = 0.5d0
      one = 1.0d0
      zero = 0.0d0
      const = one + dsqrt (2.0d0)
      Do 10 k = 1, npt
10    hcol (k) = zero
      Do 20 j = 1, npt - n - 1
         temp = zmat (knew, j)
         Do 20 k = 1, npt
20    hcol (k) = hcol (k) + temp * zmat (k, j)
      alpha = hcol (knew)
      ha = half * alpha
!
!     Calculate the gradient of the KNEW-th Lagrange function at XOPT.
!
      Do 30 i = 1, n
30    glag (i) = bmat (knew, i)
      Do 50 k = 1, npt
         temp = zero
         Do 40 j = 1, n
40       temp = temp + xpt (k, j) * xopt (j)
         temp = hcol (k) * temp
         Do 50 i = 1, n
50    glag (i) = glag (i) + temp * xpt (k, i)
!
!     Search for a large denominator along the straight lines through XOPT
!     and another interpolation point. SLBD and SUBD will be lower and upper
!     bounds on the step along each of these lines in turn. PREDSQ will be
!     set to the square of the predicted denominator for each line. PRESAV
!     will be set to the largest admissible value of PREDSQ that occurs.
!
      presav = zero
      Do 80 k = 1, npt
         If (k == kopt) Go To 80
         dderiv = zero
         distsq = zero
         Do 60 i = 1, n
            temp = xpt (k, i) - xopt (i)
            dderiv = dderiv + glag (i) * temp
60       distsq = distsq + temp * temp
         subd = adelt / dsqrt (distsq)
         slbd = - subd
         ilbd = 0
         iubd = 0
         sumin = dmin1 (one, subd)
!
!     Revise SLBD and SUBD if necessary because of the bounds in SL and SU.
!
         Do 70 i = 1, n
            temp = xpt (k, i) - xopt (i)
            If (temp > zero) Then
               If (slbd*temp < sl(i)-xopt(i)) Then
                  slbd = (sl(i)-xopt(i)) / temp
                  ilbd = - i
               End If
               If (subd*temp > su(i)-xopt(i)) Then
                  subd = dmax1 (sumin, (su(i)-xopt(i))/temp)
                  iubd = i
               End If
            Else If (temp < zero) Then
               If (slbd*temp > su(i)-xopt(i)) Then
                  slbd = (su(i)-xopt(i)) / temp
                  ilbd = i
               End If
               If (subd*temp < sl(i)-xopt(i)) Then
                  subd = dmax1 (sumin, (sl(i)-xopt(i))/temp)
                  iubd = - i
               End If
            End If
70       Continue
!
!     Seek a large modulus of the KNEW-th Lagrange function when the index
!     of the other interpolation point on the line through XOPT is KNEW.
!
         If (k == knew) Then
            diff = dderiv - one
            step = slbd
            vlag = slbd * (dderiv-slbd*diff)
            isbd = ilbd
            temp = subd * (dderiv-subd*diff)
            If (dabs(temp) > dabs(vlag)) Then
               step = subd
               vlag = temp
               isbd = iubd
            End If
            tempd = half * dderiv
            tempa = tempd - diff * slbd
            tempb = tempd - diff * subd
            If (tempa*tempb < zero) Then
               temp = tempd * tempd / diff
               If (dabs(temp) > dabs(vlag)) Then
                  step = tempd / diff
                  vlag = temp
                  isbd = 0
               End If
            End If
!
!     Search along each of the other lines through XOPT and another point.
!
         Else
            step = slbd
            vlag = slbd * (one-slbd)
            isbd = ilbd
            temp = subd * (one-subd)
            If (dabs(temp) > dabs(vlag)) Then
               step = subd
               vlag = temp
               isbd = iubd
            End If
            If (subd > half) Then
               If (dabs(vlag) < 0.25d0) Then
                  step = half
                  vlag = 0.25d0
                  isbd = 0
               End If
            End If
            vlag = vlag * dderiv
         End If
!
!     Calculate PREDSQ for the current line search and maintain PRESAV.
!
         temp = step * (one-step) * distsq
         predsq = vlag * vlag * (vlag*vlag+ha*temp*temp)
         If (predsq > presav) Then
            presav = predsq
            ksav = k
            stpsav = step
            ibdsav = isbd
         End If
80    Continue
!
!     Construct XNEW in a way that satisfies the bound constraints exactly.
!
      Do 90 i = 1, n
         temp = xopt (i) + stpsav * (xpt(ksav, i)-xopt(i))
90    xnew (i) = dmax1 (sl(i), dmin1(su(i), temp))
      If (ibdsav < 0) xnew (-ibdsav) = sl (-ibdsav)
      If (ibdsav > 0) xnew (ibdsav) = su (ibdsav)
!
!     Prepare for the iterative method that assembles the constrained Cauchy
!     step in W. The sum of squares of the fixed components of W is formed in
!     WFIXSQ, and the free components of W are set to BIGSTP.
!
      bigstp = adelt + adelt
      iflag = 0
100   wfixsq = zero
      ggfree = zero
      Do 110 i = 1, n
         w (i) = zero
         tempa = dmin1 (xopt(i)-sl(i), glag(i))
         tempb = dmax1 (xopt(i)-su(i), glag(i))
         If (tempa > zero .Or. tempb < zero) Then
            w (i) = bigstp
            ggfree = ggfree + glag (i) ** 2
         End If
110   Continue
      If (ggfree == zero) Then
         cauchy = zero
         Go To 200
      End If
!
!     Investigate whether more components of W can be fixed.
!
120   temp = adelt * adelt - wfixsq
      If (temp > zero) Then
         wsqsav = wfixsq
         step = dsqrt (temp/ggfree)
         ggfree = zero
         Do 130 i = 1, n
            If (w(i) == bigstp) Then
               temp = xopt (i) - step * glag (i)
               If (temp <= sl(i)) Then
                  w (i) = sl (i) - xopt (i)
                  wfixsq = wfixsq + w (i) ** 2
               Else If (temp >= su(i)) Then
                  w (i) = su (i) - xopt (i)
                  wfixsq = wfixsq + w (i) ** 2
               Else
                  ggfree = ggfree + glag (i) ** 2
               End If
            End If
130      Continue
         If (wfixsq > wsqsav .And. ggfree > zero) Go To 120
      End If
!
!     Set the remaining free components of W and all components of XALT,
!     except that W may be scaled later.
!
      gw = zero
      Do 140 i = 1, n
         If (w(i) == bigstp) Then
            w (i) = - step * glag (i)
            xalt (i) = dmax1 (sl(i), dmin1(su(i), xopt(i)+w(i)))
         Else If (w(i) == zero) Then
            xalt (i) = xopt (i)
         Else If (glag(i) > zero) Then
            xalt (i) = sl (i)
         Else
            xalt (i) = su (i)
         End If
140   gw = gw + glag (i) * w (i)
!
!     Set CURV to the curvature of the KNEW-th Lagrange function along W.
!     Scale W by a factor less than one if that can reduce the modulus of
!     the Lagrange function at XOPT+W. Set CAUCHY to the final value of
!     the square of this function.
!
      curv = zero
      Do 160 k = 1, npt
         temp = zero
         Do 150 j = 1, n
150      temp = temp + xpt (k, j) * w (j)
160   curv = curv + hcol (k) * temp * temp
      If (iflag == 1) curv = - curv
      If (curv >-gw .And. curv <-const*gw) Then
         scale = - gw / curv
         Do 170 i = 1, n
            temp = xopt (i) + scale * w (i)
170      xalt (i) = dmax1 (sl(i), dmin1(su(i), temp))
         cauchy = (half*gw*scale) ** 2
      Else
         cauchy = (gw+half*curv) ** 2
      End If
!
!     If IFLAG is zero, then XALT is calculated as before after reversing
!     the sign of GLAG. Thus two XALT vectors become available. The one that
!     is chosen is the one that gives the larger value of CAUCHY.
!
      If (iflag == 0) Then
         Do 180 i = 1, n
            glag (i) = - glag (i)
180      w (n+i) = xalt (i)
         csave = cauchy
         iflag = 1
         Go To 100
      End If
      If (csave > cauchy) Then
         Do 190 i = 1, n
190      xalt (i) = w (n+i)
         cauchy = csave
      End If
200   Return
End Subroutine altmov

Subroutine prelim (n, npt, x, xl, xu, rhobeg, iprint, maxfun, xbase, &
   xpt, fval, gopt, hq, pq, bmat, zmat, ndim, sl, su, nf, kopt, calfun)
   Implicit real * 8 (a-h, o-z)
   Dimension x (*), xl (*), xu (*), xbase (*), xpt (npt,*), fval (*), &
  & gopt (*), hq (*), pq (*), bmat (ndim,*), zmat (npt,*), sl (*), su &
  & (*)
    procedure(func) :: calfun

!
!     The arguments N, NPT, X, XL, XU, RHOBEG, IPRINT and MAXFUN are the
!       same as the corresponding arguments in SUBROUTINE BOBYQA.
!     The arguments XBASE, XPT, FVAL, HQ, PQ, BMAT, ZMAT, NDIM, SL and SU
!       are the same as the corresponding arguments in BOBYQB, the elements
!       of SL and SU being set in BOBYQA.
!     GOPT is usually the gradient of the quadratic model at XOPT+XBASE, but
!       it is set by PRELIM to the gradient of the quadratic model at XBASE.
!       If XOPT is nonzero, BOBYQB will change it to its usual value later.
!     NF is maintaned as the number of calls of CALFUN so far.
!     KOPT will be such that the least calculated value of F so far is at
!       the point XPT(KOPT,.)+XBASE in the space of the variables.
!
!     SUBROUTINE PRELIM sets the elements of XBASE, XPT, FVAL, GOPT, HQ, PQ,
!     BMAT and ZMAT for the first iteration, and it maintains the values of
!     NF and KOPT. The vector X is also changed by PRELIM.
!
!     Set some constants.
!
   half = 0.5d0
   one = 1.0d0
   two = 2.0d0
   zero = 0.0d0
   rhosq = rhobeg * rhobeg
   recip = one / rhosq
   np = n + 1
!
!     Set XBASE to the initial vector of variables, and set the initial
!     elements of XPT, BMAT, HQ, PQ and ZMAT to zero.
!
   Do 20 j = 1, n
      xbase (j) = x (j)
      Do 10 k = 1, npt
10    xpt (k, j) = zero
      Do 20 i = 1, ndim
20 bmat (i, j) = zero
   Do 30 ih = 1, (n*np) / 2
30 hq (ih) = zero
   Do 40 k = 1, npt
      pq (k) = zero
      Do 40 j = 1, npt - np
40 zmat (k, j) = zero
!
!     Begin the initialization procedure. NF becomes one more than the number
!     of function values so far. The coordinates of the displacement of the
!     next initial interpolation point from XBASE are set in XPT(NF+1,.).
!
   nf = 0
50 nfm = nf
   nfx = nf - n
   nf = nf + 1
   If (nfm <= 2*n) Then
      If (nfm >= 1 .And. nfm <= n) Then
         stepa = rhobeg
         If (su(nfm) == zero) stepa = - stepa
         xpt (nf, nfm) = stepa
      Else If (nfm > n) Then
         stepa = xpt (nf-n, nfx)
         stepb = - rhobeg
         If (sl(nfx) == zero) stepb = dmin1 (two*rhobeg, su(nfx))
         If (su(nfx) == zero) stepb = dmax1 (-two*rhobeg, sl(nfx))
         xpt (nf, nfx) = stepb
      End If
   Else
      itemp = (nfm-np) / n
      jpt = nfm - itemp * n - n
      ipt = jpt + itemp
      If (ipt > n) Then
         itemp = jpt
         jpt = ipt - n
         ipt = itemp
      End If
      xpt (nf, ipt) = xpt (ipt+1, ipt)
      xpt (nf, jpt) = xpt (jpt+1, jpt)
   End If
!
!     Calculate the next value of F. The least function value so far and
!     its index are required.
!
   Do 60 j = 1, n
      x (j) = dmin1 (dmax1(xl(j), xbase(j)+xpt(nf, j)), xu(j))
      If (xpt(nf, j) == sl(j)) x (j) = xl (j)
      If (xpt(nf, j) == su(j)) x (j) = xu (j)
60 Continue
   Call calfun (n, x, f)
   If (iprint == 3) Then
      Print 70, nf, f, (x(i), i=1, n)
70    Format (/ 4 x, 'Function number', i6, '    F =', 1 pd18.10,&
      '    The corresponding X is:' / (2 x, 5d15.6))
   End If
   fval (nf) = f
   If (nf == 1) Then
      fbeg = f
      kopt = 1
   Else If (f < fval(kopt)) Then
      kopt = nf
   End If
!
!     Set the nonzero initial elements of BMAT and the quadratic model in the
!     cases when NF is at most 2*N+1. If NF exceeds N+1, then the positions
!     of the NF-th and (NF-N)-th interpolation points may be switched, in
!     order that the function value at the first of them contributes to the
!     off-diagonal second derivative terms of the initial quadratic model.
!
   If (nf <= 2*n+1) Then
      If (nf >= 2 .And. nf <= n+1) Then
         gopt (nfm) = (f-fbeg) / stepa
         If (npt < nf+n) Then
            bmat (1, nfm) = - one / stepa
            bmat (nf, nfm) = one / stepa
            bmat (npt+nfm, nfm) = - half * rhosq
         End If
      Else If (nf >= n+2) Then
         ih = (nfx*(nfx+1)) / 2
         temp = (f-fbeg) / stepb
         diff = stepb - stepa
         hq (ih) = two * (temp-gopt(nfx)) / diff
         gopt (nfx) = (gopt(nfx)*stepb-temp*stepa) / diff
         If (stepa*stepb < zero) Then
            If (f < fval(nf-n)) Then
               fval (nf) = fval (nf-n)
               fval (nf-n) = f
               If (kopt == nf) kopt = nf - n
               xpt (nf-n, nfx) = stepb
               xpt (nf, nfx) = stepa
            End If
         End If
         bmat (1, nfx) = - (stepa+stepb) / (stepa*stepb)
         bmat (nf, nfx) = - half / xpt (nf-n, nfx)
         bmat (nf-n, nfx) = - bmat (1, nfx) - bmat (nf, nfx)
         zmat (1, nfx) = dsqrt (two) / (stepa*stepb)
         zmat (nf, nfx) = dsqrt (half) / rhosq
         zmat (nf-n, nfx) = - zmat (1, nfx) - zmat (nf, nfx)
      End If
!
!     Set the off-diagonal second derivatives of the Lagrange functions and
!     the initial quadratic model.
!
   Else
      ih = (ipt*(ipt-1)) / 2 + jpt
      zmat (1, nfx) = recip
      zmat (nf, nfx) = recip
      zmat (ipt+1, nfx) = - recip
      zmat (jpt+1, nfx) = - recip
      temp = xpt (nf, ipt) * xpt (nf, jpt)
      hq (ih) = (fbeg-fval(ipt+1)-fval(jpt+1)+f) / temp
   End If
   If (nf < npt .And. nf < maxfun) Go To 50
   Return
End Subroutine prelim

Subroutine rescue (n, npt, xl, xu, iprint, maxfun, xbase, xpt, fval, &
  xopt, gopt, hq, pq, bmat, zmat, ndim, sl, su, nf, delta, kopt, vlag, &
  ptsaux, ptsid, w, calfun)
   Implicit real * 8 (a-h, o-z)
   Dimension xl (*), xu (*), xbase (*), xpt (npt,*), fval (*), xopt &
  & (*), gopt (*), hq (*), pq (*), bmat (ndim,*), zmat (npt,*), sl (*), &
  & su (*), vlag (*), ptsaux (2,*), ptsid (*), w (*)
    procedure(func) :: calfun

!
!     The arguments N, NPT, XL, XU, IPRINT, MAXFUN, XBASE, XPT, FVAL, XOPT,
!       GOPT, HQ, PQ, BMAT, ZMAT, NDIM, SL and SU have the same meanings as
!       the corresponding arguments of BOBYQB on the entry to RESCUE.
!     NF is maintained as the number of calls of CALFUN so far, except that
!       NF is set to -1 if the value of MAXFUN prevents further progress.
!     KOPT is maintained so that FVAL(KOPT) is the least calculated function
!       value. Its correct value must be given on entry. It is updated if a
!       new least function value is found, but the corresponding changes to
!       XOPT and GOPT have to be made later by the calling program.
!     DELTA is the current trust region radius.
!     VLAG is a working space vector that will be used for the values of the
!       provisional Lagrange functions at each of the interpolation points.
!       They are part of a product that requires VLAG to be of length NDIM.
!     PTSAUX is also a working space array. For J=1,2,...,N, PTSAUX(1,J) and
!       PTSAUX(2,J) specify the two positions of provisional interpolation
!       points when a nonzero step is taken along e_J (the J-th coordinate
!       direction) through XBASE+XOPT, as specified below. Usually these
!       steps have length DELTA, but other lengths are chosen if necessary
!       in order to satisfy the given bounds on the variables.
!     PTSID is also a working space array. It has NPT components that denote
!       provisional new positions of the original interpolation points, in
!       case changes are needed to restore the linear independence of the
!       interpolation conditions. The K-th point is a candidate for change
!       if and only if PTSID(K) is nonzero. In this case let p and q be the
!       integer parts of PTSID(K) and (PTSID(K)-p) multiplied by N+1. If p
!       and q are both positive, the step from XBASE+XOPT to the new K-th
!       interpolation point is PTSAUX(1,p)*e_p + PTSAUX(1,q)*e_q. Otherwise
!       the step is PTSAUX(1,p)*e_p or PTSAUX(2,q)*e_q in the cases q=0 or
!       p=0, respectively.
!     The first NDIM+NPT elements of the array W are used for working space.
!     The final elements of BMAT and ZMAT are set in a well-conditioned way
!       to the values that are appropriate for the new interpolation points.
!     The elements of GOPT, HQ and PQ are also revised to the values that are
!       appropriate to the final quadratic model.
!
!     Set some constants.
!
   half = 0.5d0
   one = 1.0d0
   zero = 0.0d0
   np = n + 1
   sfrac = half / dfloat (np)
   nptm = npt - np
!
!     Shift the interpolation points so that XOPT becomes the origin, and set
!     the elements of ZMAT to zero. The value of SUMPQ is required in the
!     updating of HQ below. The squares of the distances from XOPT to the
!     other interpolation points are set at the end of W. Increments of WINC
!     may be added later to these squares to balance the consideration of
!     the choice of point that is going to become current.
!
   sumpq = zero
   winc = zero
   Do 20 k = 1, npt
      distsq = zero
      Do 10 j = 1, n
         xpt (k, j) = xpt (k, j) - xopt (j)
10    distsq = distsq + xpt (k, j) ** 2
      sumpq = sumpq + pq (k)
      w (ndim+k) = distsq
      winc = dmax1 (winc, distsq)
      Do 20 j = 1, nptm
20 zmat (k, j) = zero
!
!     Update HQ so that HQ and PQ define the second derivatives of the model
!     after XBASE has been shifted to the trust region centre.
!
   ih = 0
   Do 40 j = 1, n
      w (j) = half * sumpq * xopt (j)
      Do 30 k = 1, npt
30    w (j) = w (j) + pq (k) * xpt (k, j)
      Do 40 i = 1, j
         ih = ih + 1
40 hq (ih) = hq (ih) + w (i) * xopt (j) + w (j) * xopt (i)
!
!     Shift XBASE, SL, SU and XOPT. Set the elements of BMAT to zero, and
!     also set the elements of PTSAUX.
!
   Do 50 j = 1, n
      xbase (j) = xbase (j) + xopt (j)
      sl (j) = sl (j) - xopt (j)
      su (j) = su (j) - xopt (j)
      xopt (j) = zero
      ptsaux (1, j) = dmin1 (delta, su(j))
      ptsaux (2, j) = dmax1 (-delta, sl(j))
      If (ptsaux(1, j)+ptsaux(2, j) < zero) Then
         temp = ptsaux (1, j)
         ptsaux (1, j) = ptsaux (2, j)
         ptsaux (2, j) = temp
      End If
      If (dabs(ptsaux(2, j)) < half*dabs(ptsaux(1, j))) Then
         ptsaux (2, j) = half * ptsaux (1, j)
      End If
      Do 50 i = 1, ndim
50 bmat (i, j) = zero
   fbase = fval (kopt)
!
!     Set the identifiers of the artificial interpolation points that are
!     along a coordinate direction from XOPT, and set the corresponding
!     nonzero elements of BMAT and ZMAT.
!
   ptsid (1) = sfrac
   Do 60 j = 1, n
      jp = j + 1
      jpn = jp + n
      ptsid (jp) = dfloat (j) + sfrac
      If (jpn <= npt) Then
         ptsid (jpn) = dfloat (j) / dfloat (np) + sfrac
         temp = one / (ptsaux(1, j)-ptsaux(2, j))
         bmat (jp, j) = - temp + one / ptsaux (1, j)
         bmat (jpn, j) = temp + one / ptsaux (2, j)
         bmat (1, j) = - bmat (jp, j) - bmat (jpn, j)
         zmat (1, j) = dsqrt (2.0d0) / dabs (ptsaux(1, j)*ptsaux(2, j))
         zmat (jp, j) = zmat (1, j) * ptsaux (2, j) * temp
         zmat (jpn, j) = - zmat (1, j) * ptsaux (1, j) * temp
      Else
         bmat (1, j) = - one / ptsaux (1, j)
         bmat (jp, j) = one / ptsaux (1, j)
         bmat (j+npt, j) = - half * ptsaux (1, j) ** 2
      End If
60 Continue
!
!     Set any remaining identifiers with their nonzero elements of ZMAT.
!
   If (npt >= n+np) Then
      Do 70 k = 2 * np, npt
         iw = (dfloat(k-np)-half) / dfloat (n)
         ip = k - np - iw * n
         iq = ip + iw
         If (iq > n) iq = iq - n
         ptsid (k) = dfloat (ip) + dfloat (iq) / dfloat (np) + sfrac
         temp = one / (ptsaux(1, ip)*ptsaux(1, iq))
         zmat (1, k-np) = temp
         zmat (ip+1, k-np) = - temp
         zmat (iq+1, k-np) = - temp
70    zmat (k, k-np) = temp
   End If
   nrem = npt
   kold = 1
   knew = kopt
!
!     Reorder the provisional points in the way that exchanges PTSID(KOLD)
!     with PTSID(KNEW).
!
80 Do 90 j = 1, n
      temp = bmat (kold, j)
      bmat (kold, j) = bmat (knew, j)
90 bmat (knew, j) = temp
   Do 100 j = 1, nptm
      temp = zmat (kold, j)
      zmat (kold, j) = zmat (knew, j)
100 zmat (knew, j) = temp
   ptsid (kold) = ptsid (knew)
   ptsid (knew) = zero
   w (ndim+knew) = zero
   nrem = nrem - 1
   If (knew /= kopt) Then
      temp = vlag (kold)
      vlag (kold) = vlag (knew)
      vlag (knew) = temp
!
!     Update the BMAT and ZMAT matrices so that the status of the KNEW-th
!     interpolation point can be changed from provisional to original. The
!     branch to label 350 occurs if all the original points are reinstated.
!     The nonnegative values of W(NDIM+K) are required in the search below.
!
      Call update (n, npt, bmat, zmat, ndim, vlag, beta, denom, knew, &
     & w)
      If (nrem == 0) Go To 350
      Do 110 k = 1, npt
110   w (ndim+k) = dabs (w(ndim+k))
   End If
!
!     Pick the index KNEW of an original interpolation point that has not
!     yet replaced one of the provisional interpolation points, giving
!     attention to the closeness to XOPT and to previous tries with KNEW.
!
120 dsqmin = zero
   Do 130 k = 1, npt
      If (w(ndim+k) > zero) Then
         If (dsqmin == zero .Or. w(ndim+k) < dsqmin) Then
            knew = k
            dsqmin = w (ndim+k)
         End If
      End If
130 Continue
   If (dsqmin == zero) Go To 260
!
!     Form the W-vector of the chosen original interpolation point.
!
   Do 140 j = 1, n
140 w (npt+j) = xpt (knew, j)
   Do 160 k = 1, npt
      sum = zero
      If (k == kopt) Then
         Continue
      Else If (ptsid(k) == zero) Then
         Do 150 j = 1, n
150      sum = sum + w (npt+j) * xpt (k, j)
      Else
         ip = ptsid (k)
         If (ip > 0) sum = w (npt+ip) * ptsaux (1, ip)
         iq = dfloat (np) * ptsid (k) - dfloat (ip*np)
         If (iq > 0) Then
            iw = 1
            If (ip == 0) iw = 2
            sum = sum + w (npt+iq) * ptsaux (iw, iq)
         End If
      End If
160 w (k) = half * sum * sum
!
!     Calculate VLAG and BETA for the required updating of the H matrix if
!     XPT(KNEW,.) is reinstated in the set of interpolation points.
!
   Do 180 k = 1, npt
      sum = zero
      Do 170 j = 1, n
170   sum = sum + bmat (k, j) * w (npt+j)
180 vlag (k) = sum
   beta = zero
   Do 200 j = 1, nptm
      sum = zero
      Do 190 k = 1, npt
190   sum = sum + zmat (k, j) * w (k)
      beta = beta - sum * sum
      Do 200 k = 1, npt
200 vlag (k) = vlag (k) + sum * zmat (k, j)
   bsum = zero
   distsq = zero
   Do 230 j = 1, n
      sum = zero
      Do 210 k = 1, npt
210   sum = sum + bmat (k, j) * w (k)
      jp = j + npt
      bsum = bsum + sum * w (jp)
      Do 220 ip = npt + 1, ndim
220   sum = sum + bmat (ip, j) * w (ip)
      bsum = bsum + sum * w (jp)
      vlag (jp) = sum
230 distsq = distsq + xpt (knew, j) ** 2
   beta = half * distsq * distsq + beta - bsum
   vlag (kopt) = vlag (kopt) + one
!
!     KOLD is set to the index of the provisional interpolation point that is
!     going to be deleted to make way for the KNEW-th original interpolation
!     point. The choice of KOLD is governed by the avoidance of a small value
!     of the denominator in the updating calculation of UPDATE.
!
   denom = zero
   vlmxsq = zero
   Do 250 k = 1, npt
      If (ptsid(k) /= zero) Then
         hdiag = zero
         Do 240 j = 1, nptm
240      hdiag = hdiag + zmat (k, j) ** 2
         den = beta * hdiag + vlag (k) ** 2
         If (den > denom) Then
            kold = k
            denom = den
         End If
      End If
250 vlmxsq = dmax1 (vlmxsq, vlag(k)**2)
   If (denom <= 1.0d-2*vlmxsq) Then
      w (ndim+knew) = - w (ndim+knew) - winc
      Go To 120
   End If
   Go To 80
!
!     When label 260 is reached, all the final positions of the interpolation
!     points have been chosen although any changes have not been included yet
!     in XPT. Also the final BMAT and ZMAT matrices are complete, but, apart
!     from the shift of XBASE, the updating of the quadratic model remains to
!     be done. The following cycle through the new interpolation points begins
!     by putting the new point in XPT(KPT,.) and by setting PQ(KPT) to zero,
!     except that a RETURN occurs if MAXFUN prohibits another value of F.
!
260 Do 340 kpt = 1, npt
      If (ptsid(kpt) == zero) Go To 340
      If (nf >= maxfun) Then
         nf = - 1
         Go To 350
      End If
      ih = 0
      Do 270 j = 1, n
         w (j) = xpt (kpt, j)
         xpt (kpt, j) = zero
         temp = pq (kpt) * w (j)
         Do 270 i = 1, j
            ih = ih + 1
270   hq (ih) = hq (ih) + temp * w (i)
      pq (kpt) = zero
      ip = ptsid (kpt)
      iq = dfloat (np) * ptsid (kpt) - dfloat (ip*np)
      If (ip > 0) Then
         xp = ptsaux (1, ip)
         xpt (kpt, ip) = xp
      End If
      If (iq > 0) Then
         xq = ptsaux (1, iq)
         If (ip == 0) xq = ptsaux (2, iq)
         xpt (kpt, iq) = xq
      End If
!
!     Set VQUAD to the value of the current model at the new point.
!
      vquad = fbase
      If (ip > 0) Then
         ihp = (ip+ip*ip) / 2
         vquad = vquad + xp * (gopt(ip)+half*xp*hq(ihp))
      End If
      If (iq > 0) Then
         ihq = (iq+iq*iq) / 2
         vquad = vquad + xq * (gopt(iq)+half*xq*hq(ihq))
         If (ip > 0) Then
            iw = max0 (ihp, ihq) - iabs (ip-iq)
            vquad = vquad + xp * xq * hq (iw)
         End If
      End If
      Do 280 k = 1, npt
         temp = zero
         If (ip > 0) temp = temp + xp * xpt (k, ip)
         If (iq > 0) temp = temp + xq * xpt (k, iq)
280   vquad = vquad + half * pq (k) * temp * temp
!
!     Calculate F at the new interpolation point, and set DIFF to the factor
!     that is going to multiply the KPT-th Lagrange function when the model
!     is updated to provide interpolation to the new function value.
!
      Do 290 i = 1, n
         w (i) = dmin1 (dmax1(xl(i), xbase(i)+xpt(kpt, i)), xu(i))
         If (xpt(kpt, i) == sl(i)) w (i) = xl (i)
         If (xpt(kpt, i) == su(i)) w (i) = xu (i)
290   Continue
      nf = nf + 1
      Call calfun (n, w, f)
      If (iprint == 3) Then
         Print 300, nf, f, (w(i), i=1, n)
300      Format (/ 4 x, 'Function number', i6, '    F =', 1 pd18.10, &
         '    The corresponding X is:' / (2 x, 5d15.6))
      End If
      fval (kpt) = f
      If (f < fval(kopt)) kopt = kpt
      diff = f - vquad
!
!     Update the quadratic model. The RETURN from the subroutine occurs when
!     all the new interpolation points are included in the model.
!
      Do 310 i = 1, n
310   gopt (i) = gopt (i) + diff * bmat (kpt, i)
      Do 330 k = 1, npt
         sum = zero
         Do 320 j = 1, nptm
320      sum = sum + zmat (k, j) * zmat (kpt, j)
         temp = diff * sum
         If (ptsid(k) == zero) Then
            pq (k) = pq (k) + temp
         Else
            ip = ptsid (k)
            iq = dfloat (np) * ptsid (k) - dfloat (ip*np)
            ihq = (iq*iq+iq) / 2
            If (ip == 0) Then
               hq (ihq) = hq (ihq) + temp * ptsaux (2, iq) ** 2
            Else
               ihp = (ip*ip+ip) / 2
               hq (ihp) = hq (ihp) + temp * ptsaux (1, ip) ** 2
               If (iq > 0) Then
                  hq (ihq) = hq (ihq) + temp * ptsaux (1, iq) ** 2
                  iw = max0 (ihp, ihq) - iabs (iq-ip)
                  hq (iw) = hq (iw) + temp * ptsaux (1, ip) * ptsaux &
                 & (1, iq)
               End If
            End If
         End If
330   Continue
      ptsid (kpt) = zero
340 Continue
350 Return
End Subroutine rescue

Subroutine trsbox (n, npt, xpt, xopt, gopt, hq, pq, sl, su, delta, &
   xnew, d, gnew, xbdi, s, hs, hred, dsq, crvmin)
   Implicit real * 8 (a-h, o-z)
   Dimension xpt (npt,*), xopt (*), gopt (*), hq (*), pq (*), sl (*), &
  & su (*), xnew (*), d (*), gnew (*), xbdi (*), s (*), hs (*), hred &
  & (*)
!
!     The arguments N, NPT, XPT, XOPT, GOPT, HQ, PQ, SL and SU have the same
!       meanings as the corresponding arguments of BOBYQB.
!     DELTA is the trust region radius for the present calculation, which
!       seeks a small value of the quadratic model within distance DELTA of
!       XOPT subject to the bounds on the variables.
!     XNEW will be set to a new vector of variables that is approximately
!       the one that minimizes the quadratic model within the trust region
!       subject to the SL and SU constraints on the variables. It satisfies
!       as equations the bounds that become active during the calculation.
!     D is the calculated trial step from XOPT, generated iteratively from an
!       initial value of zero. Thus XNEW is XOPT+D after the final iteration.
!     GNEW holds the gradient of the quadratic model at XOPT+D. It is updated
!       when D is updated.
!     XBDI is a working space vector. For I=1,2,...,N, the element XBDI(I) is
!       set to -1.0, 0.0, or 1.0, the value being nonzero if and only if the
!       I-th variable has become fixed at a bound, the bound being SL(I) or
!       SU(I) in the case XBDI(I)=-1.0 or XBDI(I)=1.0, respectively. This
!       information is accumulated during the construction of XNEW.
!     The arrays S, HS and HRED are also used for working space. They hold the
!       current search direction, and the changes in the gradient of Q along S
!       and the reduced D, respectively, where the reduced D is the same as D,
!       except that the components of the fixed variables are zero.
!     DSQ will be set to the square of the length of XNEW-XOPT.
!     CRVMIN is set to zero if D reaches the trust region boundary. Otherwise
!       it is set to the least curvature of H that occurs in the conjugate
!       gradient searches that are not restricted by any constraints. The
!       value CRVMIN=-1.0D0 is set, however, if all of these searches are
!       constrained.
!
!     A version of the truncated conjugate gradient is applied. If a line
!     search is restricted by a constraint, then the procedure is restarted,
!     the values of the variables that are at their bounds being fixed. If
!     the trust region boundary is reached, then further changes may be made
!     to D, each one being in the two dimensional space that is spanned
!     by the current D and the gradient of Q at XOPT+D, staying on the trust
!     region boundary. Termination occurs when the reduction in Q seems to
!     be close to the greatest reduction that can be achieved.
!
!     Set some constants.
!
   half = 0.5d0
   one = 1.0d0
   onemin = - 1.0d0
   zero = 0.0d0
!
!     The sign of GOPT(I) gives the sign of the change to the I-th variable
!     that will reduce Q from its value at XOPT. Thus XBDI(I) shows whether
!     or not to fix the I-th variable at one of its bounds initially, with
!     NACT being set to the number of fixed variables. D and GNEW are also
!     set for the first iteration. DELSQ is the upper bound on the sum of
!     squares of the free variables. QRED is the reduction in Q so far.
!
   iterc = 0
   nact = 0
   sqstp = zero
   Do 10 i = 1, n
      xbdi (i) = zero
      If (xopt(i) <= sl(i)) Then
         If (gopt(i) >= zero) xbdi (i) = onemin
      Else If (xopt(i) >= su(i)) Then
         If (gopt(i) <= zero) xbdi (i) = one
      End If
      If (xbdi(i) /= zero) nact = nact + 1
      d (i) = zero
10 gnew (i) = gopt (i)
   delsq = delta * delta
   qred = zero
   crvmin = onemin
!
!     Set the next search direction of the conjugate gradient method. It is
!     the steepest descent direction initially and when the iterations are
!     restarted because a variable has just been fixed by a bound, and of
!     course the components of the fixed variables are zero. ITERMAX is an
!     upper bound on the indices of the conjugate gradient iterations.
!
20 beta = zero
30 stepsq = zero
   Do 40 i = 1, n
      If (xbdi(i) /= zero) Then
         s (i) = zero
      Else If (beta == zero) Then
         s (i) = - gnew (i)
      Else
         s (i) = beta * s (i) - gnew (i)
      End If
40 stepsq = stepsq + s (i) ** 2
   If (stepsq == zero) Go To 190
   If (beta == zero) Then
      gredsq = stepsq
      itermax = iterc + n - nact
   End If
   If (gredsq*delsq <= 1.0d-4*qred*qred) Go To 190
!
!     Multiply the search direction by the second derivative matrix of Q and
!     calculate some scalars for the choice of steplength. Then set BLEN to
!     the length of the the step to the trust region boundary and STPLEN to
!     the steplength, ignoring the simple bounds.
!
   Go To 210
50 resid = delsq
   ds = zero
   shs = zero
   Do 60 i = 1, n
      If (xbdi(i) == zero) Then
         resid = resid - d (i) ** 2
         ds = ds + s (i) * d (i)
         shs = shs + s (i) * hs (i)
      End If
60 Continue
   If (resid <= zero) Go To 90
   temp = dsqrt (stepsq*resid+ds*ds)
   If (ds < zero) Then
      blen = (temp-ds) / stepsq
   Else
      blen = resid / (temp+ds)
   End If
   stplen = blen
   If (shs > zero) Then
      stplen = dmin1 (blen, gredsq/shs)
   End If
!
!
!     Reduce STPLEN if necessary in order to preserve the simple bounds,
!     letting IACT be the index of the new constrained variable.
!
   iact = 0
   Do 70 i = 1, n
      If (s(i) /= zero) Then
         xsum = xopt (i) + d (i)
         If (s(i) > zero) Then
            temp = (su(i)-xsum) / s (i)
         Else
            temp = (sl(i)-xsum) / s (i)
         End If
         If (temp < stplen) Then
            stplen = temp
            iact = i
         End If
      End If
70 Continue
!
!     Update CRVMIN, GNEW and D. Set SDEC to the decrease that occurs in Q.
!
   sdec = zero
   If (stplen > zero) Then
      iterc = iterc + 1
      temp = shs / stepsq
      If (iact == 0 .And. temp > zero) Then
         crvmin = dmin1 (crvmin, temp)
         If (crvmin == onemin) crvmin = temp
      End If
      ggsav = gredsq
      gredsq = zero
      Do 80 i = 1, n
         gnew (i) = gnew (i) + stplen * hs (i)
         If (xbdi(i) == zero) gredsq = gredsq + gnew (i) ** 2
80    d (i) = d (i) + stplen * s (i)
      sdec = dmax1 (stplen*(ggsav-half*stplen*shs), zero)
      qred = qred + sdec
   End If
!
!     Restart the conjugate gradient method if it has hit a new bound.
!
   If (iact > 0) Then
      nact = nact + 1
      xbdi (iact) = one
      If (s(iact) < zero) xbdi (iact) = onemin
      delsq = delsq - d (iact) ** 2
      If (delsq <= zero) Go To 90
      Go To 20
   End If
!
!     If STPLEN is less than BLEN, then either apply another conjugate
!     gradient iteration or RETURN.
!
   If (stplen < blen) Then
      If (iterc == itermax) Go To 190
      If (sdec <= 0.01d0*qred) Go To 190
      beta = gredsq / ggsav
      Go To 30
   End If
90 crvmin = zero
!
!     Prepare for the alternative iteration by calculating some scalars and
!     by multiplying the reduced D by the second derivative matrix of Q.
!
100 If (nact >= n-1) Go To 190
   dredsq = zero
   dredg = zero
   gredsq = zero
   Do 110 i = 1, n
      If (xbdi(i) == zero) Then
         dredsq = dredsq + d (i) ** 2
         dredg = dredg + d (i) * gnew (i)
         gredsq = gredsq + gnew (i) ** 2
         s (i) = d (i)
      Else
         s (i) = zero
      End If
110 Continue
   itcsav = iterc
   Go To 210
!
!     Let the search direction S be a linear combination of the reduced D
!     and the reduced G that is orthogonal to the reduced D.
!
120 iterc = iterc + 1
   temp = gredsq * dredsq - dredg * dredg
   If (temp <= 1.0d-4*qred*qred) Go To 190
   temp = dsqrt (temp)
   Do 130 i = 1, n
      If (xbdi(i) == zero) Then
         s (i) = (dredg*d(i)-dredsq*gnew(i)) / temp
      Else
         s (i) = zero
      End If
130 Continue
   sredg = - temp
!
!     By considering the simple bounds on the variables, calculate an upper
!     bound on the tangent of half the angle of the alternative iteration,
!     namely ANGBD, except that, if already a free variable has reached a
!     bound, there is a branch back to label 100 after fixing that variable.
!
   angbd = one
   iact = 0
   Do 140 i = 1, n
      If (xbdi(i) == zero) Then
         tempa = xopt (i) + d (i) - sl (i)
         tempb = su (i) - xopt (i) - d (i)
         If (tempa <= zero) Then
            nact = nact + 1
            xbdi (i) = onemin
            Go To 100
         Else If (tempb <= zero) Then
            nact = nact + 1
            xbdi (i) = one
            Go To 100
         End If
         ratio = one
         ssq = d (i) ** 2 + s (i) ** 2
         temp = ssq - (xopt(i)-sl(i)) ** 2
         If (temp > zero) Then
            temp = dsqrt (temp) - s (i)
            If (angbd*temp > tempa) Then
               angbd = tempa / temp
               iact = i
               xsav = onemin
            End If
         End If
         temp = ssq - (su(i)-xopt(i)) ** 2
         If (temp > zero) Then
            temp = dsqrt (temp) + s (i)
            If (angbd*temp > tempb) Then
               angbd = tempb / temp
               iact = i
               xsav = one
            End If
         End If
      End If
140 Continue
!
!     Calculate HHD and some curvatures for the alternative iteration.
!
   Go To 210
150 shs = zero
   dhs = zero
   dhd = zero
   Do 160 i = 1, n
      If (xbdi(i) == zero) Then
         shs = shs + s (i) * hs (i)
         dhs = dhs + d (i) * hs (i)
         dhd = dhd + d (i) * hred (i)
      End If
160 Continue
!
!     Seek the greatest reduction in Q for a range of equally spaced values
!     of ANGT in [0,ANGBD], where ANGT is the tangent of half the angle of
!     the alternative iteration.
!
   redmax = zero
   isav = 0
   redsav = zero
   iu = 17.0d0 * angbd + 3.1d0
   Do 170 i = 1, iu
      angt = angbd * dfloat (i) / dfloat (iu)
      sth = (angt+angt) / (one+angt*angt)
      temp = shs + angt * (angt*dhd-dhs-dhs)
      rednew = sth * (angt*dredg-sredg-half*sth*temp)
      If (rednew > redmax) Then
         redmax = rednew
         isav = i
         rdprev = redsav
      Else If (i == isav+1) Then
         rdnext = rednew
      End If
170 redsav = rednew
!
!     Return if the reduction is zero. Otherwise, set the sine and cosine
!     of the angle of the alternative iteration, and calculate SDEC.
!
   If (isav == 0) Go To 190
   If (isav < iu) Then
      temp = (rdnext-rdprev) / (redmax+redmax-rdprev-rdnext)
      angt = angbd * (dfloat(isav)+half*temp) / dfloat (iu)
   End If
   cth = (one-angt*angt) / (one+angt*angt)
   sth = (angt+angt) / (one+angt*angt)
   temp = shs + angt * (angt*dhd-dhs-dhs)
   sdec = sth * (angt*dredg-sredg-half*sth*temp)
   If (sdec <= zero) Go To 190
!
!     Update GNEW, D and HRED. If the angle of the alternative iteration
!     is restricted by a bound on a free variable, that variable is fixed
!     at the bound.
!
   dredg = zero
   gredsq = zero
   Do 180 i = 1, n
      gnew (i) = gnew (i) + (cth-one) * hred (i) + sth * hs (i)
      If (xbdi(i) == zero) Then
         d (i) = cth * d (i) + sth * s (i)
         dredg = dredg + d (i) * gnew (i)
         gredsq = gredsq + gnew (i) ** 2
      End If
180 hred (i) = cth * hred (i) + sth * hs (i)
   qred = qred + sdec
   If (iact > 0 .And. isav == iu) Then
      nact = nact + 1
      xbdi (iact) = xsav
      Go To 100
   End If
!
!     If SDEC is sufficiently small, then RETURN after setting XNEW to
!     XOPT+D, giving careful attention to the bounds.
!
   If (sdec > 0.01d0*qred) Go To 120
190 dsq = zero
   Do 200 i = 1, n
      xnew (i) = dmax1 (dmin1(xopt(i)+d(i), su(i)), sl(i))
      If (xbdi(i) == onemin) xnew (i) = sl (i)
      If (xbdi(i) == one) xnew (i) = su (i)
      d (i) = xnew (i) - xopt (i)
200 dsq = dsq + d (i) ** 2
   Return
!
!     The following instructions multiply the current S-vector by the second
!     derivative matrix of the quadratic model, putting the product in HS.
!     They are reached from three different parts of the software above and
!     they can be regarded as an external subroutine.
!
210 ih = 0
   Do 220 j = 1, n
      hs (j) = zero
      Do 220 i = 1, j
         ih = ih + 1
         If (i < j) hs (j) = hs (j) + hq (ih) * s (i)
220 hs (i) = hs (i) + hq (ih) * s (j)
   Do 250 k = 1, npt
      If (pq(k) /= zero) Then
         temp = zero
         Do 230 j = 1, n
230      temp = temp + xpt (k, j) * s (j)
         temp = temp * pq (k)
         Do 240 i = 1, n
240      hs (i) = hs (i) + temp * xpt (k, i)
      End If
250 Continue
   If (crvmin /= zero) Go To 50
   If (iterc > itcsav) Go To 150
   Do 260 i = 1, n
260 hred (i) = hs (i)
   Go To 120
End Subroutine trsbox

Subroutine update (n, npt, bmat, zmat, ndim, vlag, beta, denom, knew,w)

   Implicit real * 8 (a-h, o-z)
   Dimension bmat (ndim,*), zmat (npt,*), vlag (*), w (*)
!
!     The arrays BMAT and ZMAT are updated, as required by the new position
!     of the interpolation point that has the index KNEW. The vector VLAG has
!     N+NPT components, set on entry to the first NPT and last N components
!     of the product Hw in equation (4.11) of the Powell (2006) paper on
!     NEWUOA. Further, BETA is set on entry to the value of the parameter
!     with that name, and DENOM is set to the denominator of the updating
!     formula. Elements of ZMAT may be treated as zero if their moduli are
!     at most ZTEST. The first NDIM elements of W are used for working space.
!
!     Set some constants.
!
   one = 1.0d0
   zero = 0.0d0
   nptm = npt - n - 1
   ztest = zero
   Do 10 k = 1, npt
      Do 10 j = 1, nptm
10 ztest = dmax1 (ztest, dabs(zmat(k, j)))
   ztest = 1.0d-20 * ztest
!
!     Apply the rotations that put zeros in the KNEW-th row of ZMAT.
!
   jl = 1
   Do 30 j = 2, nptm
      If (dabs(zmat(knew, j)) > ztest) Then
         temp = dsqrt (zmat(knew, 1)**2+zmat(knew, j)**2)
         tempa = zmat (knew, 1) / temp
         tempb = zmat (knew, j) / temp
         Do 20 i = 1, npt
            temp = tempa * zmat (i, 1) + tempb * zmat (i, j)
            zmat (i, j) = tempa * zmat (i, j) - tempb * zmat (i, 1)
20       zmat (i, 1) = temp
      End If
      zmat (knew, j) = zero
30 Continue
!
!     Put the first NPT components of the KNEW-th column of HLAG into W,
!     and calculate the parameters of the updating formula.
!
   Do 40 i = 1, npt
      w (i) = zmat (knew, 1) * zmat (i, 1)
40 Continue
   alpha = w (knew)
   tau = vlag (knew)
   vlag (knew) = vlag (knew) - one
!
!     Complete the updating of ZMAT.
!
   temp = dsqrt (denom)
   tempb = zmat (knew, 1) / temp
   tempa = tau / temp
   Do 50 i = 1, npt
50 zmat (i, 1) = tempa * zmat (i, 1) - tempb * vlag (i)
!
!     Finally, update the matrix BMAT.
!
   Do 60 j = 1, n
      jp = npt + j
      w (jp) = bmat (knew, j)
      tempa = (alpha*vlag(jp)-tau*w(jp)) / denom
      tempb = (-beta*w(jp)-tau*vlag(jp)) / denom
      Do 60 i = 1, jp
         bmat (i, j) = bmat (i, j) + tempa * vlag (i) + tempb * w (i)
         If (i > npt) bmat (jp, i-npt) = bmat (i, j)
60 Continue
   Return
End Subroutine update

subroutine bobyqa_test()
!
!     Test problem for BOBYQA, the objective function being the sum of
!     the reciprocals of all pairwise distances between the points P_I,
!     I=1,2,...,M in two dimensions, where M=N/2 and where the components
!     of P_I are X(2*I-1) and X(2*I). Thus each vector X of N variables
!     defines the M points P_I. The initial X gives equally spaced points
!     on a circle. Four different choices of the pairs (N,NPT) are tried,
!     namely (10,16), (10,21), (20,26) and (20,41). Convergence to a local
!     minimum that is not global occurs in both the N=10 cases. The details
!     of the results are highly sensitive to computer rounding errors. The
!     choice IPRINT=2 provides the current X and optimal F so far whenever
!     RHO is reduced. The bound constraints of the problem require every
!     component of X to be in the interval [-1,1].
!
    Implicit real * 8 (a-h, o-z)
    Dimension x (100), xl (100), xu (100), w (500000)
    twopi = 8.0d0 * datan (1.0d0)
    bdl = - 1.0d0
    bdu = 1.0d0
    iprint = 2
    maxfun = 500000
    rhobeg = 1.0d-1
    rhoend = 1.0d-6
    m = 5
    10 n = 2 * m
    Do 20 i = 1, n
       xl (i) = bdl
    20 xu (i) = bdu
    Do 50 jcase = 1, 2
       npt = n + 6
       If (jcase == 2) npt = 2 * n + 1
       Print 30, m, n, npt
    30 Format (/ / 5 x, '2D output with M =', i4, ',  N =', i4,&
       '  and  NPT =', i4)
       Do 40 j = 1, m
          temp = dfloat (j) * twopi / dfloat (m)
          x (2*j-1) = dcos (temp)
    40 x (2*j) = dsin (temp)
       Call bobyqa (n, npt, x, xl, xu, rhobeg, rhoend, iprint, maxfun, w, calfun)
    50 Continue
    m = m + m
    If (m <= 10) Go To 10

contains

    Subroutine calfun (n, x, f)
          Implicit real * 8 (a-h, o-z)
          Dimension x (*)
          f = 0.0d0
          Do 10 i = 4, n, 2
             Do 10 j = 2, i - 2, 2
                temp = (x(i-1)-x(j-1)) ** 2 + (x(i)-x(j)) ** 2
                temp = dmax1 (temp, 1.0d-6)
    10    f = f + 1.0d0 / dsqrt (temp)
          Return
    End Subroutine calfun

End subroutine bobyqa_test

end module bobyqa_module