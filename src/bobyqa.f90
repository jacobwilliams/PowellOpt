 
module bobyqa_module
 
    use kind_module, only: wp
 
    private
 
    abstract interface
    subroutine func (n, x, f)!! calfun interface
        import :: wp
        implicit none
        integer :: n
        real (wp) :: x (*)
        real (wp) :: f
    end subroutine func
    end interface
 
    public :: bobyqa
    public :: bobyqa_test
 
contains
 
    subroutine bobyqa (n, npt, x, xl, xu, rhobeg, rhoend, iprint, maxfun, w, calfun)
        implicit real (wp) (a-h, o-z)
        dimension x (*), xl (*), xu (*), w (*)
        procedure (func) :: calfun
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
        if (npt < n+2 .or. npt > ((n+2)*np)/2) then
            write(*,'(/4X,A)') &
                'Return from BOBYQA because NPT is not in the required interval'
            return
        end if
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
        zero = 0.0_wp
        do j = 1, n
            temp = xu (j) - xl (j)
            if (temp < rhobeg+rhobeg) then
                write(*,'(/4X,A)') &
                    'Return from BOBYQA because one of the differences XU(I)-XL(I) is less than 2*RHOBEG.'
                return
            end if
            jsl = isl + j - 1
            jsu = jsl + n
            w (jsl) = xl (j) - x (j)
            w (jsu) = xu (j) - x (j)
            if (w(jsl) >=-rhobeg) then
                if (w(jsl) >= zero) then
                    x (j) = xl (j)
                    w (jsl) = zero
                    w (jsu) = temp
                else
                    x (j) = xl (j) + rhobeg
                    w (jsl) = - rhobeg
                    w (jsu) = max (xu(j)-x(j), rhobeg)
                end if
            else if (w(jsu) <= rhobeg) then
                if (w(jsu) <= zero) then
                    x (j) = xu (j)
                    w (jsl) = - temp
                    w (jsu) = zero
                else
                    x (j) = xu (j) - rhobeg
                    w (jsl) = min (xl(j)-x(j),-rhobeg)
                    w (jsu) = rhobeg
                end if
            end if
        end do
!
!     Make the call of BOBYQB.
!
        call bobyqb (n, npt, x, xl, xu, rhobeg, rhoend, iprint, maxfun, w(ixb), w(ixp), &
       & w(ifv), w(ixo), w(igo), w(ihq), w(ipq), w(ibmat), w(izmat), ndim, w(isl), &
       & w(isu), w(ixn), w(ixa), w(id), w(ivl), w(iw), calfun)
 
    end subroutine bobyqa
 
    subroutine bobyqb (n, npt, x, xl, xu, rhobeg, rhoend, iprint, maxfun, xbase, xpt, &
   & fval, xopt, gopt, hq, pq, bmat, zmat, ndim, sl, su, xnew, xalt, d, vlag, w, calfun)
        implicit real (wp) (a-h, o-z)
        dimension x (*), xl (*), xu (*), xbase (*), xpt (npt,*), fval (*), xopt (*), gopt &
       & (*), hq (*), pq (*), bmat (ndim,*), zmat (npt,*), sl (*), su (*), xnew (*), xalt &
       & (*), d (*), vlag (*), w (*)
        procedure (func) :: calfun
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
        half = 0.5_wp
        one = 1.0_wp
        ten = 10.0_wp
        tenth = 0.1_wp
        two = 2.0_wp
        zero = 0.0_wp
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
        call prelim (n, npt, x, xl, xu, rhobeg, iprint, maxfun, xbase, xpt, fval, gopt, &
       & hq, pq, bmat, zmat, ndim, sl, su, nf, kopt, calfun)
        xoptsq = zero
        do 10 i = 1, n
            xopt (i) = xpt (kopt, i)
10      xoptsq = xoptsq + xopt (i) ** 2
        fsave = fval (1)
        if (nf < npt) then
            if (iprint > 0) write(*,'(/4X,A)') &
                'Return from BOBYQA because CALFUN has been called MAXFUN times.'
            go to 720
        end if
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
20      if (kopt /= kbase) then
            ih = 0
            do 30 j = 1, n
                do 30 i = 1, j
                    ih = ih + 1
                    if (i < j) gopt (j) = gopt (j) + hq (ih) * xopt (i)
30          gopt (i) = gopt (i) + hq (ih) * xopt (j)
            if (nf > npt) then
                do 50 k = 1, npt
                    temp = zero
                    do 40 j = 1, n
40                  temp = temp + xpt (k, j) * xopt (j)
                    temp = pq (k) * temp
                    do 50 i = 1, n
50              gopt (i) = gopt (i) + temp * xpt (k, i)
            end if
        end if
!
!     Generate the next point in the trust region that provides a small value
!     of the quadratic model subject to the constraints on the variables.
!     The integer NTRITS is set to the number "trust region" iterations that
!     have occurred since the last "alternative" iteration. If the length
!     of XNEW-XOPT is less than HALF*RHO, however, then there is a branch to
!     label 650 or 680 with NTRITS=-1, instead of calculating F at XNEW.
!
60      call trsbox (n, npt, xpt, xopt, gopt, hq, pq, sl, su, delta, xnew, d, w, w(np), &
       & w(np+n), w(np+2*n), w(np+3*n), dsq, crvmin)
        dnorm = min (delta, sqrt(dsq))
        if (dnorm < half*rho) then
            ntrits = - 1
            distsq = (ten*rho) ** 2
            if (nf <= nfsav+2) go to 650
!
!     The following choice between labels 650 and 680 depends on whether or
!     not our work with the current RHO seems to be complete. Either RHO is
!     decreased or termination occurs if the errors in the quadratic model at
!     the last three interpolation points compare favourably with predictions
!     of likely improvements to the model within distance HALF*RHO of XOPT.
!
            errbig = max (diffa, diffb, diffc)
            frhosq = 0.125_wp * rho * rho
            if (crvmin > zero .and. errbig > frhosq*crvmin) go to 650
            bdtol = errbig / rho
            do j = 1, n
                bdtest = bdtol
                if (xnew(j) == sl(j)) bdtest = w (j)
                if (xnew(j) == su(j)) bdtest = - w (j)
                if (bdtest < bdtol) then
                    curv = hq ((j+j*j)/2)
                    do k = 1, npt
                        curv = curv + pq (k) * xpt (k, j) ** 2
                    end do
                    bdtest = bdtest + half * curv * rho
                    if (bdtest < bdtol) go to 650
                end if
            end do
            go to 680
        end if
        ntrits = ntrits + 1
!
!     Severe cancellation is likely to occur if XOPT is too far from XBASE.
!     If the following test holds, then XBASE is shifted so that XOPT becomes
!     zero. The appropriate changes are made to BMAT and to the second
!     derivatives of the current model, beginning with the changes to BMAT
!     that do not depend on ZMAT. VLAG is used temporarily for working space.
!
90      if (dsq <= 1.0e-3_wp*xoptsq) then
            fracsq = 0.25_wp * xoptsq
            sumpq = zero
            do 110 k = 1, npt
                sumpq = sumpq + pq (k)
                sum = - half * xoptsq
                do 100 i = 1, n
100             sum = sum + xpt (k, i) * xopt (i)
                w (npt+k) = sum
                temp = fracsq - half * sum
                do 110 i = 1, n
                    w (i) = bmat (k, i)
                    vlag (i) = sum * xpt (k, i) + temp * xopt (i)
                    ip = npt + i
                    do 110 j = 1, i
110         bmat (ip, j) = bmat (ip, j) + w (i) * vlag (j) + vlag (i) * w (j)
!
!     Then the revisions of BMAT that depend on ZMAT are calculated.
!
            do 150 jj = 1, nptm
                sumz = zero
                sumw = zero
                do 120 k = 1, npt
                    sumz = sumz + zmat (k, jj)
                    vlag (k) = w (npt+k) * zmat (k, jj)
120             sumw = sumw + vlag (k)
                do 140 j = 1, n
                    sum = (fracsq*sumz-half*sumw) * xopt (j)
                    do 130 k = 1, npt
130                 sum = sum + vlag (k) * xpt (k, j)
                    w (j) = sum
                    do 140 k = 1, npt
140             bmat (k, j) = bmat (k, j) + sum * zmat (k, jj)
                do 150 i = 1, n
                    ip = i + npt
                    temp = w (i)
                    do 150 j = 1, i
150         bmat (ip, j) = bmat (ip, j) + temp * w (j)
!
!     The following instructions complete the shift, including the changes
!     to the second derivative parameters of the quadratic model.
!
            ih = 0
            do 170 j = 1, n
                w (j) = - half * sumpq * xopt (j)
                do 160 k = 1, npt
                    w (j) = w (j) + pq (k) * xpt (k, j)
160             xpt (k, j) = xpt (k, j) - xopt (j)
                do 170 i = 1, j
                    ih = ih + 1
                    hq (ih) = hq (ih) + w (i) * xopt (j) + xopt (i) * w (j)
170         bmat (npt+i, j) = bmat (npt+j, i)
            do 180 i = 1, n
                xbase (i) = xbase (i) + xopt (i)
                xnew (i) = xnew (i) - xopt (i)
                sl (i) = sl (i) - xopt (i)
                su (i) = su (i) - xopt (i)
180         xopt (i) = zero
            xoptsq = zero
        end if
        if (ntrits == 0) go to 210
        go to 230
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
190     nfsav = nf
        kbase = kopt
        call rescue (n, npt, xl, xu, iprint, maxfun, xbase, xpt, fval, xopt, gopt, hq, &
       & pq, bmat, zmat, ndim, sl, su, nf, delta, kopt, vlag, w, w(n+np), w(ndim+np), &
       & calfun)
!
!     XOPT is updated now in case the branch below to label 720 is taken.
!     Any updating of GOPT occurs after the branch below to label 20, which
!     leads to a trust region iteration as does the branch to label 60.
!
        xoptsq = zero
        if (kopt /= kbase) then
            do 200 i = 1, n
                xopt (i) = xpt (kopt, i)
200         xoptsq = xoptsq + xopt (i) ** 2
        end if
        if (nf < 0) then
            nf = maxfun
            if (iprint > 0) write(*,'(/4X,A)') &
                'Return from BOBYQA because CALFUN has been called MAXFUN times.'
            go to 720
        end if
        nresc = nf
        if (nfsav < nf) then
            nfsav = nf
            go to 20
        end if
        if (ntrits > 0) go to 60
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
210     call altmov (n, npt, xpt, xopt, bmat, zmat, ndim, sl, su, kopt, knew, adelt, &
       & xnew, xalt, alpha, cauchy, w, w(np), w(ndim+1))
        do 220 i = 1, n
220     d (i) = xnew (i) - xopt (i)
!
!     Calculate VLAG and BETA for the current choice of D. The scalar
!     product of D with XPT(K,.) is going to be held in W(NPT+K) for
!     use when VQUAD is calculated.
!
230     do 250 k = 1, npt
            suma = zero
            sumb = zero
            sum = zero
            do 240 j = 1, n
                suma = suma + xpt (k, j) * d (j)
                sumb = sumb + xpt (k, j) * xopt (j)
240         sum = sum + bmat (k, j) * d (j)
            w (k) = suma * (half*suma+sumb)
            vlag (k) = sum
250     w (npt+k) = suma
        beta = zero
        do 270 jj = 1, nptm
            sum = zero
            do 260 k = 1, npt
260         sum = sum + zmat (k, jj) * w (k)
            beta = beta - sum * sum
            do 270 k = 1, npt
270     vlag (k) = vlag (k) + sum * zmat (k, jj)
        dsq = zero
        bsum = zero
        dx = zero
        do 300 j = 1, n
            dsq = dsq + d (j) ** 2
            sum = zero
            do 280 k = 1, npt
280         sum = sum + w (k) * bmat (k, j)
            bsum = bsum + sum * d (j)
            jp = npt + j
            do 290 i = 1, n
290         sum = sum + bmat (jp, i) * d (i)
            vlag (jp) = sum
            bsum = bsum + sum * d (j)
300     dx = dx + d (j) * xopt (j)
        beta = dx * dx + dsq * (xoptsq+dx+dx+half*dsq) + beta - bsum
        vlag (kopt) = vlag (kopt) + one
!
!     If NTRITS is zero, the denominator may be increased by replacing
!     the step D of ALTMOV by a Cauchy step. Then RESCUE may be called if
!     rounding errors have damaged the chosen denominator.
!
        if (ntrits == 0) then
            denom = vlag (knew) ** 2 + alpha * beta
            if (denom < cauchy .and. cauchy > zero) then
                do 310 i = 1, n
                    xnew (i) = xalt (i)
310             d (i) = xnew (i) - xopt (i)
                cauchy = zero
                go to 230
            end if
            if (denom <= half*vlag(knew)**2) then
                if (nf > nresc) go to 190
                if (iprint > 0) write(*,'(/5X,A)') &
                    'Return from BOBYQA because of much cancellation in a denominator.'
                go to 720
            end if
!
!     Alternatively, if NTRITS is positive, then set KNEW to the index of
!     the next interpolation point to be deleted to make room for a trust
!     region step. Again RESCUE may be called if rounding errors have damaged
!     the chosen denominator, which is the reason for attempting to select
!     KNEW before calculating the next value of the objective function.
!
        else
            delsq = delta * delta
            scaden = zero
            biglsq = zero
            knew = 0
            do 350 k = 1, npt
                if (k == kopt) go to 350
                hdiag = zero
                do 330 jj = 1, nptm
330             hdiag = hdiag + zmat (k, jj) ** 2
                den = beta * hdiag + vlag (k) ** 2
                distsq = zero
                do 340 j = 1, n
340             distsq = distsq + (xpt(k, j)-xopt(j)) ** 2
                temp = max (one, (distsq/delsq)**2)
                if (temp*den > scaden) then
                    scaden = temp * den
                    knew = k
                    denom = den
                end if
                biglsq = max (biglsq, temp*vlag(k)**2)
350         continue
            if (scaden <= half*biglsq) then
                if (nf > nresc) go to 190
                if (iprint > 0) write(*,'(/5X,A)') &
                    'Return from BOBYQA because of much cancellation in a denominator.'
                go to 720
            end if
        end if
!
!     Put the variables for the next calculation of the objective function
!       in XNEW, with any adjustments for the bounds.
!
!
!     Calculate the value of the objective function at XBASE+XNEW, unless
!       the limit on the number of calculations of F has been reached.
!
360     do 380 i = 1, n
            x (i) = min (max(xl(i), xbase(i)+xnew(i)), xu(i))
            if (xnew(i) == sl(i)) x (i) = xl (i)
            if (xnew(i) == su(i)) x (i) = xu (i)
380     continue
        if (nf >= maxfun) then
            if (iprint > 0) write(*,'(/4X,A)') &
                'Return from BOBYQA because CALFUN has been called MAXFUN times.'
            go to 720
        end if
        nf = nf + 1
        call calfun (n, x, f)
        if (iprint == 3) then
            print 400, nf, f, (x(i), i=1, n)
400         format (/ 4 x, 'Function number', i6, '    F =', 1 pd18.10,&
            '   The corresponding X is:' / (2 x, 5d15.6))
        end if
        if (ntrits ==-1) then
            fsave = f
            go to 720
        end if
!
!     Use the quadratic model to predict the change in F due to the step D,
!       and set DIFF to the error of this prediction.
!
        fopt = fval (kopt)
        vquad = zero
        ih = 0
        do 410 j = 1, n
            vquad = vquad + d (j) * gopt (j)
            do 410 i = 1, j
                ih = ih + 1
                temp = d (i) * d (j)
                if (i == j) temp = half * temp
410     vquad = vquad + hq (ih) * temp
        do 420 k = 1, npt
420     vquad = vquad + half * pq (k) * w (npt+k) ** 2
        diff = f - fopt - vquad
        diffc = diffb
        diffb = diffa
        diffa = abs (diff)
        if (dnorm > rho) nfsav = nf
!
!     Pick the next value of DELTA after a trust region step.
!
        if (ntrits > 0) then
            if (vquad >= zero) then
                if (iprint > 0) print 430
430             format (/ 4 x, 'Return from BOBYQA because a trust',&
                ' region step has failed to reduce Q.')
                go to 720
            end if
            ratio = (f-fopt) / vquad
            if (ratio <= tenth) then
                delta = min (half*delta, dnorm)
            else if (ratio <= 0.7_wp) then
                delta = max (half*delta, dnorm)
            else
                delta = max (half*delta, dnorm+dnorm)
            end if
            if (delta <= 1.5_wp*rho) delta = rho
!
!     Recalculate KNEW and DENOM if the new F is less than FOPT.
!
            if (f < fopt) then
                ksav = knew
                densav = denom
                delsq = delta * delta
                scaden = zero
                biglsq = zero
                knew = 0
                do 460 k = 1, npt
                    hdiag = zero
                    do 440 jj = 1, nptm
440                 hdiag = hdiag + zmat (k, jj) ** 2
                    den = beta * hdiag + vlag (k) ** 2
                    distsq = zero
                    do 450 j = 1, n
450                 distsq = distsq + (xpt(k, j)-xnew(j)) ** 2
                    temp = max (one, (distsq/delsq)**2)
                    if (temp*den > scaden) then
                        scaden = temp * den
                        knew = k
                        denom = den
                    end if
460             biglsq = max (biglsq, temp*vlag(k)**2)
                if (scaden <= half*biglsq) then
                    knew = ksav
                    denom = densav
                end if
            end if
        end if
!
!     Update BMAT and ZMAT, so that the KNEW-th interpolation point can be
!     moved. Also update the second derivative terms of the model.
!
        call update (n, npt, bmat, zmat, ndim, vlag, beta, denom, knew, w)
        ih = 0
        pqold = pq (knew)
        pq (knew) = zero
        do 470 i = 1, n
            temp = pqold * xpt (knew, i)
            do 470 j = 1, i
                ih = ih + 1
470     hq (ih) = hq (ih) + temp * xpt (knew, j)
        do 480 jj = 1, nptm
            temp = diff * zmat (knew, jj)
            do 480 k = 1, npt
480     pq (k) = pq (k) + temp * zmat (k, jj)
!
!     Include the new interpolation point, and make the changes to GOPT at
!     the old XOPT that are caused by the updating of the quadratic model.
!
        fval (knew) = f
        do 490 i = 1, n
            xpt (knew, i) = xnew (i)
490     w (i) = bmat (knew, i)
        do 520 k = 1, npt
            suma = zero
            do 500 jj = 1, nptm
500         suma = suma + zmat (knew, jj) * zmat (k, jj)
            sumb = zero
            do 510 j = 1, n
510         sumb = sumb + xpt (k, j) * xopt (j)
            temp = suma * sumb
            do 520 i = 1, n
520     w (i) = w (i) + temp * xpt (k, i)
        do 530 i = 1, n
530     gopt (i) = gopt (i) + diff * w (i)
!
!     Update XOPT, GOPT and KOPT if the new calculated F is less than FOPT.
!
        if (f < fopt) then
            kopt = knew
            xoptsq = zero
            ih = 0
            do 540 j = 1, n
                xopt (j) = xnew (j)
                xoptsq = xoptsq + xopt (j) ** 2
                do 540 i = 1, j
                    ih = ih + 1
                    if (i < j) gopt (j) = gopt (j) + hq (ih) * d (i)
540         gopt (i) = gopt (i) + hq (ih) * d (j)
            do 560 k = 1, npt
                temp = zero
                do 550 j = 1, n
550             temp = temp + xpt (k, j) * d (j)
                temp = pq (k) * temp
                do 560 i = 1, n
560         gopt (i) = gopt (i) + temp * xpt (k, i)
        end if
!
!     Calculate the parameters of the least Frobenius norm interpolant to
!     the current data, the gradient of this interpolant at XOPT being put
!     into VLAG(NPT+I), I=1,2,...,N.
!
        if (ntrits > 0) then
            do 570 k = 1, npt
                vlag (k) = fval (k) - fval (kopt)
570         w (k) = zero
            do 590 j = 1, nptm
                sum = zero
                do 580 k = 1, npt
580             sum = sum + zmat (k, j) * vlag (k)
                do 590 k = 1, npt
590         w (k) = w (k) + sum * zmat (k, j)
            do 610 k = 1, npt
                sum = zero
                do 600 j = 1, n
600             sum = sum + xpt (k, j) * xopt (j)
                w (k+npt) = w (k)
610         w (k) = sum * w (k)
            gqsq = zero
            gisq = zero
            do 630 i = 1, n
                sum = zero
                do 620 k = 1, npt
620             sum = sum + bmat (k, i) * vlag (k) + xpt (k, i) * w (k)
                if (xopt(i) == sl(i)) then
                    gqsq = gqsq + min (zero, gopt(i)) ** 2
                    gisq = gisq + min (zero, sum) ** 2
                else if (xopt(i) == su(i)) then
                    gqsq = gqsq + max (zero, gopt(i)) ** 2
                    gisq = gisq + max (zero, sum) ** 2
                else
                    gqsq = gqsq + gopt (i) ** 2
                    gisq = gisq + sum * sum
                end if
630         vlag (npt+i) = sum
!
!     Test whether to replace the new quadratic model by the least Frobenius
!     norm interpolant, making the replacement if the test is satisfied.
!
            itest = itest + 1
            if (gqsq < ten*gisq) itest = 0
            if (itest >= 3) then
                do 640 i = 1, max0 (npt, nh)
                    if (i <= n) gopt (i) = vlag (npt+i)
                    if (i <= npt) pq (i) = w (npt+i)
                    if (i <= nh) hq (i) = zero
                    itest = 0
640             continue
            end if
        end if
!
!     If a trust region step has provided a sufficient decrease in F, then
!     branch for another trust region calculation. The case NTRITS=0 occurs
!     when the new interpolation point was reached by an alternative step.
!
        if (ntrits == 0) go to 60
        if (f <= fopt+tenth*vquad) go to 60
!
!     Alternatively, find out if the interpolation points are close enough
!       to the best point so far.
!
        distsq = max ((two*delta)**2, (ten*rho)**2)
650     knew = 0
        do 670 k = 1, npt
            sum = zero
            do 660 j = 1, n
660         sum = sum + (xpt(k, j)-xopt(j)) ** 2
            if (sum > distsq) then
                knew = k
                distsq = sum
            end if
670     continue
!
!     If KNEW is positive, then ALTMOV finds alternative new positions for
!     the KNEW-th interpolation point within distance ADELT of XOPT. It is
!     reached via label 90. Otherwise, there is a branch to label 60 for
!     another trust region iteration, unless the calculations with the
!     current RHO are complete.
!
        if (knew > 0) then
            dist = sqrt (distsq)
            if (ntrits ==-1) then
                delta = min (tenth*delta, half*dist)
                if (delta <= 1.5_wp*rho) delta = rho
            end if
            ntrits = 0
            adelt = max (min(tenth*dist, delta), rho)
            dsq = adelt * adelt
            go to 90
        end if
        if (ntrits ==-1) go to 680
        if (ratio > zero) go to 60
        if (max(delta, dnorm) > rho) go to 60
!
!     The calculations with the current value of RHO are complete. Pick the
!       next values of RHO and DELTA.
!
680     if (rho > rhoend) then
            delta = half * rho
            ratio = rho / rhoend
            if (ratio <= 16.0_wp) then
                rho = rhoend
            else if (ratio <= 250.0_wp) then
                rho = sqrt (ratio) * rhoend
            else
                rho = tenth * rho
            end if
            delta = max (delta, rho)
            if (iprint >= 2) then
                if (iprint >= 3) print 690
690             format (5 x)
                print 700, rho, nf
700             format (/ 4 x, 'New RHO =', 1 pd11.4, 5 x, 'Number of',&
                ' function values =', i6)
                print 710, fval (kopt), (xbase(i)+xopt(i), i=1, n)
710             format (4 x, 'Least value of F =', 1 pd23.15, 9 x,&
                'The corresponding X is:'/(2 x, 5d15.6))
            end if
            ntrits = 0
            nfsav = nf
            go to 60
        end if
!
!     Return from the calculation, after another Newton-Raphson step, if
!       it is too short to have been tried before.
!
        if (ntrits ==-1) go to 360
720     if (fval(kopt) <= fsave) then
            do 730 i = 1, n
                x (i) = min (max(xl(i), xbase(i)+xopt(i)), xu(i))
                if (xopt(i) == sl(i)) x (i) = xl (i)
                if (xopt(i) == su(i)) x (i) = xu (i)
730         continue
            f = fval (kopt)
        end if
        if (iprint >= 1) then
            print 740, nf
740         format (/ 4 x, 'At the return from BOBYQA', 5 x,&
            'Number of function values =', i6)
            print 710, f, (x(i), i=1, n)
        end if
        return
    end subroutine bobyqb
 
    subroutine altmov (n, npt, xpt, xopt, bmat, zmat, ndim, sl, su, kopt, knew, adelt, &
   & xnew, xalt, alpha, cauchy, glag, hcol, w)
 
        implicit real (wp) (a-h, o-z)
        dimension xpt (npt,*), xopt (*), bmat (ndim,*), zmat (npt,*), sl (*), su (*), &
       & xnew (*), xalt (*), glag (*), hcol (*), w (*)
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
        half = 0.5_wp
        one = 1.0_wp
        zero = 0.0_wp
        const = one + sqrt (2.0_wp)
        do 10 k = 1, npt
10      hcol (k) = zero
        do 20 j = 1, npt - n - 1
            temp = zmat (knew, j)
            do 20 k = 1, npt
20      hcol (k) = hcol (k) + temp * zmat (k, j)
        alpha = hcol (knew)
        ha = half * alpha
!
!     Calculate the gradient of the KNEW-th Lagrange function at XOPT.
!
        do 30 i = 1, n
30      glag (i) = bmat (knew, i)
        do 50 k = 1, npt
            temp = zero
            do 40 j = 1, n
40          temp = temp + xpt (k, j) * xopt (j)
            temp = hcol (k) * temp
            do 50 i = 1, n
50      glag (i) = glag (i) + temp * xpt (k, i)
!
!     Search for a large denominator along the straight lines through XOPT
!     and another interpolation point. SLBD and SUBD will be lower and upper
!     bounds on the step along each of these lines in turn. PREDSQ will be
!     set to the square of the predicted denominator for each line. PRESAV
!     will be set to the largest admissible value of PREDSQ that occurs.
!
        presav = zero
        do 80 k = 1, npt
            if (k == kopt) go to 80
            dderiv = zero
            distsq = zero
            do 60 i = 1, n
                temp = xpt (k, i) - xopt (i)
                dderiv = dderiv + glag (i) * temp
60          distsq = distsq + temp * temp
            subd = adelt / sqrt (distsq)
            slbd = - subd
            ilbd = 0
            iubd = 0
            sumin = min (one, subd)
!
!     Revise SLBD and SUBD if necessary because of the bounds in SL and SU.
!
            do 70 i = 1, n
                temp = xpt (k, i) - xopt (i)
                if (temp > zero) then
                    if (slbd*temp < sl(i)-xopt(i)) then
                        slbd = (sl(i)-xopt(i)) / temp
                        ilbd = - i
                    end if
                    if (subd*temp > su(i)-xopt(i)) then
                        subd = max (sumin, (su(i)-xopt(i))/temp)
                        iubd = i
                    end if
                else if (temp < zero) then
                    if (slbd*temp > su(i)-xopt(i)) then
                        slbd = (su(i)-xopt(i)) / temp
                        ilbd = i
                    end if
                    if (subd*temp < sl(i)-xopt(i)) then
                        subd = max (sumin, (sl(i)-xopt(i))/temp)
                        iubd = - i
                    end if
                end if
70          continue
!
!     Seek a large modulus of the KNEW-th Lagrange function when the index
!     of the other interpolation point on the line through XOPT is KNEW.
!
            if (k == knew) then
                diff = dderiv - one
                step = slbd
                vlag = slbd * (dderiv-slbd*diff)
                isbd = ilbd
                temp = subd * (dderiv-subd*diff)
                if (abs(temp) > abs(vlag)) then
                    step = subd
                    vlag = temp
                    isbd = iubd
                end if
                tempd = half * dderiv
                tempa = tempd - diff * slbd
                tempb = tempd - diff * subd
                if (tempa*tempb < zero) then
                    temp = tempd * tempd / diff
                    if (abs(temp) > abs(vlag)) then
                        step = tempd / diff
                        vlag = temp
                        isbd = 0
                    end if
                end if
!
!     Search along each of the other lines through XOPT and another point.
!
            else
                step = slbd
                vlag = slbd * (one-slbd)
                isbd = ilbd
                temp = subd * (one-subd)
                if (abs(temp) > abs(vlag)) then
                    step = subd
                    vlag = temp
                    isbd = iubd
                end if
                if (subd > half) then
                    if (abs(vlag) < 0.25_wp) then
                        step = half
                        vlag = 0.25_wp
                        isbd = 0
                    end if
                end if
                vlag = vlag * dderiv
            end if
!
!     Calculate PREDSQ for the current line search and maintain PRESAV.
!
            temp = step * (one-step) * distsq
            predsq = vlag * vlag * (vlag*vlag+ha*temp*temp)
            if (predsq > presav) then
                presav = predsq
                ksav = k
                stpsav = step
                ibdsav = isbd
            end if
80      continue
!
!     Construct XNEW in a way that satisfies the bound constraints exactly.
!
        do 90 i = 1, n
            temp = xopt (i) + stpsav * (xpt(ksav, i)-xopt(i))
90      xnew (i) = max (sl(i), min(su(i), temp))
        if (ibdsav < 0) xnew (-ibdsav) = sl (-ibdsav)
        if (ibdsav > 0) xnew (ibdsav) = su (ibdsav)
!
!     Prepare for the iterative method that assembles the constrained Cauchy
!     step in W. The sum of squares of the fixed components of W is formed in
!     WFIXSQ, and the free components of W are set to BIGSTP.
!
        bigstp = adelt + adelt
        iflag = 0
100     wfixsq = zero
        ggfree = zero
        do 110 i = 1, n
            w (i) = zero
            tempa = min (xopt(i)-sl(i), glag(i))
            tempb = max (xopt(i)-su(i), glag(i))
            if (tempa > zero .or. tempb < zero) then
                w (i) = bigstp
                ggfree = ggfree + glag (i) ** 2
            end if
110     continue
        if (ggfree == zero) then
            cauchy = zero
            go to 200
        end if
!
!     Investigate whether more components of W can be fixed.
!
120     temp = adelt * adelt - wfixsq
        if (temp > zero) then
            wsqsav = wfixsq
            step = sqrt (temp/ggfree)
            ggfree = zero
            do 130 i = 1, n
                if (w(i) == bigstp) then
                    temp = xopt (i) - step * glag (i)
                    if (temp <= sl(i)) then
                        w (i) = sl (i) - xopt (i)
                        wfixsq = wfixsq + w (i) ** 2
                    else if (temp >= su(i)) then
                        w (i) = su (i) - xopt (i)
                        wfixsq = wfixsq + w (i) ** 2
                    else
                        ggfree = ggfree + glag (i) ** 2
                    end if
                end if
130         continue
            if (wfixsq > wsqsav .and. ggfree > zero) go to 120
        end if
!
!     Set the remaining free components of W and all components of XALT,
!     except that W may be scaled later.
!
        gw = zero
        do 140 i = 1, n
            if (w(i) == bigstp) then
                w (i) = - step * glag (i)
                xalt (i) = max (sl(i), min(su(i), xopt(i)+w(i)))
            else if (w(i) == zero) then
                xalt (i) = xopt (i)
            else if (glag(i) > zero) then
                xalt (i) = sl (i)
            else
                xalt (i) = su (i)
            end if
140     gw = gw + glag (i) * w (i)
!
!     Set CURV to the curvature of the KNEW-th Lagrange function along W.
!     Scale W by a factor less than one if that can reduce the modulus of
!     the Lagrange function at XOPT+W. Set CAUCHY to the final value of
!     the square of this function.
!
        curv = zero
        do 160 k = 1, npt
            temp = zero
            do 150 j = 1, n
150         temp = temp + xpt (k, j) * w (j)
160     curv = curv + hcol (k) * temp * temp
        if (iflag == 1) curv = - curv
        if (curv >-gw .and. curv <-const*gw) then
            scale = - gw / curv
            do 170 i = 1, n
                temp = xopt (i) + scale * w (i)
170         xalt (i) = max (sl(i), min(su(i), temp))
            cauchy = (half*gw*scale) ** 2
        else
            cauchy = (gw+half*curv) ** 2
        end if
!
!     If IFLAG is zero, then XALT is calculated as before after reversing
!     the sign of GLAG. Thus two XALT vectors become available. The one that
!     is chosen is the one that gives the larger value of CAUCHY.
!
        if (iflag == 0) then
            do 180 i = 1, n
                glag (i) = - glag (i)
180         w (n+i) = xalt (i)
            csave = cauchy
            iflag = 1
            go to 100
        end if
        if (csave > cauchy) then
            do 190 i = 1, n
190         xalt (i) = w (n+i)
            cauchy = csave
        end if
200     return
    end subroutine altmov
 
    subroutine prelim (n, npt, x, xl, xu, rhobeg, iprint, maxfun, xbase, xpt, fval, gopt, &
   & hq, pq, bmat, zmat, ndim, sl, su, nf, kopt, calfun)
        implicit real (wp) (a-h, o-z)
        dimension x (*), xl (*), xu (*), xbase (*), xpt (npt,*), fval (*), gopt (*), hq &
       & (*), pq (*), bmat (ndim,*), zmat (npt,*), sl (*), su (*)
        procedure (func) :: calfun
 
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
        half = 0.5_wp
        one = 1.0_wp
        two = 2.0_wp
        zero = 0.0_wp
        rhosq = rhobeg * rhobeg
        recip = one / rhosq
        np = n + 1
!
!     Set XBASE to the initial vector of variables, and set the initial
!     elements of XPT, BMAT, HQ, PQ and ZMAT to zero.
!
        do 20 j = 1, n
            xbase (j) = x (j)
            do 10 k = 1, npt
10          xpt (k, j) = zero
            do 20 i = 1, ndim
20      bmat (i, j) = zero
        do 30 ih = 1, (n*np) / 2
30      hq (ih) = zero
        do 40 k = 1, npt
            pq (k) = zero
            do 40 j = 1, npt - np
40      zmat (k, j) = zero
!
!     Begin the initialization procedure. NF becomes one more than the number
!     of function values so far. The coordinates of the displacement of the
!     next initial interpolation point from XBASE are set in XPT(NF+1,.).
!
        nf = 0
50      nfm = nf
        nfx = nf - n
        nf = nf + 1
        if (nfm <= 2*n) then
            if (nfm >= 1 .and. nfm <= n) then
                stepa = rhobeg
                if (su(nfm) == zero) stepa = - stepa
                xpt (nf, nfm) = stepa
            else if (nfm > n) then
                stepa = xpt (nf-n, nfx)
                stepb = - rhobeg
                if (sl(nfx) == zero) stepb = min (two*rhobeg, su(nfx))
                if (su(nfx) == zero) stepb = max (-two*rhobeg, sl(nfx))
                xpt (nf, nfx) = stepb
            end if
        else
            itemp = (nfm-np) / n
            jpt = nfm - itemp * n - n
            ipt = jpt + itemp
            if (ipt > n) then
                itemp = jpt
                jpt = ipt - n
                ipt = itemp
            end if
            xpt (nf, ipt) = xpt (ipt+1, ipt)
            xpt (nf, jpt) = xpt (jpt+1, jpt)
        end if
!
!     Calculate the next value of F. The least function value so far and
!     its index are required.
!
        do 60 j = 1, n
            x (j) = min (max(xl(j), xbase(j)+xpt(nf, j)), xu(j))
            if (xpt(nf, j) == sl(j)) x (j) = xl (j)
            if (xpt(nf, j) == su(j)) x (j) = xu (j)
60      continue
        call calfun (n, x, f)
        if (iprint == 3) then
            print 70, nf, f, (x(i), i=1, n)
70          format (/ 4 x, 'Function number', i6, '    F =', 1 pd18.10,&
            '    The corresponding X is:' / (2 x, 5d15.6))
        end if
        fval (nf) = f
        if (nf == 1) then
            fbeg = f
            kopt = 1
        else if (f < fval(kopt)) then
            kopt = nf
        end if
!
!     Set the nonzero initial elements of BMAT and the quadratic model in the
!     cases when NF is at most 2*N+1. If NF exceeds N+1, then the positions
!     of the NF-th and (NF-N)-th interpolation points may be switched, in
!     order that the function value at the first of them contributes to the
!     off-diagonal second derivative terms of the initial quadratic model.
!
        if (nf <= 2*n+1) then
            if (nf >= 2 .and. nf <= n+1) then
                gopt (nfm) = (f-fbeg) / stepa
                if (npt < nf+n) then
                    bmat (1, nfm) = - one / stepa
                    bmat (nf, nfm) = one / stepa
                    bmat (npt+nfm, nfm) = - half * rhosq
                end if
            else if (nf >= n+2) then
                ih = (nfx*(nfx+1)) / 2
                temp = (f-fbeg) / stepb
                diff = stepb - stepa
                hq (ih) = two * (temp-gopt(nfx)) / diff
                gopt (nfx) = (gopt(nfx)*stepb-temp*stepa) / diff
                if (stepa*stepb < zero) then
                    if (f < fval(nf-n)) then
                        fval (nf) = fval (nf-n)
                        fval (nf-n) = f
                        if (kopt == nf) kopt = nf - n
                        xpt (nf-n, nfx) = stepb
                        xpt (nf, nfx) = stepa
                    end if
                end if
                bmat (1, nfx) = - (stepa+stepb) / (stepa*stepb)
                bmat (nf, nfx) = - half / xpt (nf-n, nfx)
                bmat (nf-n, nfx) = - bmat (1, nfx) - bmat (nf, nfx)
                zmat (1, nfx) = sqrt (two) / (stepa*stepb)
                zmat (nf, nfx) = sqrt (half) / rhosq
                zmat (nf-n, nfx) = - zmat (1, nfx) - zmat (nf, nfx)
            end if
!
!     Set the off-diagonal second derivatives of the Lagrange functions and
!     the initial quadratic model.
!
        else
            ih = (ipt*(ipt-1)) / 2 + jpt
            zmat (1, nfx) = recip
            zmat (nf, nfx) = recip
            zmat (ipt+1, nfx) = - recip
            zmat (jpt+1, nfx) = - recip
            temp = xpt (nf, ipt) * xpt (nf, jpt)
            hq (ih) = (fbeg-fval(ipt+1)-fval(jpt+1)+f) / temp
        end if
        if (nf < npt .and. nf < maxfun) go to 50
        return
    end subroutine prelim
 
    subroutine rescue (n, npt, xl, xu, iprint, maxfun, xbase, xpt, fval, xopt, gopt, hq, &
   & pq, bmat, zmat, ndim, sl, su, nf, delta, kopt, vlag, ptsaux, ptsid, w, calfun)
        implicit real (wp) (a-h, o-z)
        dimension xl (*), xu (*), xbase (*), xpt (npt,*), fval (*), xopt (*), gopt (*), &
       & hq (*), pq (*), bmat (ndim,*), zmat (npt,*), sl (*), su (*), vlag (*), ptsaux &
       & (2,*), ptsid (*), w (*)
        procedure (func) :: calfun
 
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
        half = 0.5_wp
        one = 1.0_wp
        zero = 0.0_wp
        np = n + 1
        sfrac = half / real (np, wp)
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
        do 20 k = 1, npt
            distsq = zero
            do 10 j = 1, n
                xpt (k, j) = xpt (k, j) - xopt (j)
10          distsq = distsq + xpt (k, j) ** 2
            sumpq = sumpq + pq (k)
            w (ndim+k) = distsq
            winc = max (winc, distsq)
            do 20 j = 1, nptm
20      zmat (k, j) = zero
!
!     Update HQ so that HQ and PQ define the second derivatives of the model
!     after XBASE has been shifted to the trust region centre.
!
        ih = 0
        do 40 j = 1, n
            w (j) = half * sumpq * xopt (j)
            do 30 k = 1, npt
30          w (j) = w (j) + pq (k) * xpt (k, j)
            do 40 i = 1, j
                ih = ih + 1
40      hq (ih) = hq (ih) + w (i) * xopt (j) + w (j) * xopt (i)
!
!     Shift XBASE, SL, SU and XOPT. Set the elements of BMAT to zero, and
!     also set the elements of PTSAUX.
!
        do 50 j = 1, n
            xbase (j) = xbase (j) + xopt (j)
            sl (j) = sl (j) - xopt (j)
            su (j) = su (j) - xopt (j)
            xopt (j) = zero
            ptsaux (1, j) = min (delta, su(j))
            ptsaux (2, j) = max (-delta, sl(j))
            if (ptsaux(1, j)+ptsaux(2, j) < zero) then
                temp = ptsaux (1, j)
                ptsaux (1, j) = ptsaux (2, j)
                ptsaux (2, j) = temp
            end if
            if (abs(ptsaux(2, j)) < half*abs(ptsaux(1, j))) then
                ptsaux (2, j) = half * ptsaux (1, j)
            end if
            do 50 i = 1, ndim
50      bmat (i, j) = zero
        fbase = fval (kopt)
!
!     Set the identifiers of the artificial interpolation points that are
!     along a coordinate direction from XOPT, and set the corresponding
!     nonzero elements of BMAT and ZMAT.
!
        ptsid (1) = sfrac
        do 60 j = 1, n
            jp = j + 1
            jpn = jp + n
            ptsid (jp) = real (j, wp) + sfrac
            if (jpn <= npt) then
                ptsid (jpn) = real (j, wp) / real (np, wp) + sfrac
                temp = one / (ptsaux(1, j)-ptsaux(2, j))
                bmat (jp, j) = - temp + one / ptsaux (1, j)
                bmat (jpn, j) = temp + one / ptsaux (2, j)
                bmat (1, j) = - bmat (jp, j) - bmat (jpn, j)
                zmat (1, j) = sqrt (2.0_wp) / abs (ptsaux(1, j)*ptsaux(2, j))
                zmat (jp, j) = zmat (1, j) * ptsaux (2, j) * temp
                zmat (jpn, j) = - zmat (1, j) * ptsaux (1, j) * temp
            else
                bmat (1, j) = - one / ptsaux (1, j)
                bmat (jp, j) = one / ptsaux (1, j)
                bmat (j+npt, j) = - half * ptsaux (1, j) ** 2
            end if
60      continue
!
!     Set any remaining identifiers with their nonzero elements of ZMAT.
!
        if (npt >= n+np) then
            do 70 k = 2 * np, npt
                iw = (real(k-np, wp)-half) / real (n, wp)
                ip = k - np - iw * n
                iq = ip + iw
                if (iq > n) iq = iq - n
                ptsid (k) = real (ip, wp) + real (iq, wp) / real (np, wp) + sfrac
                temp = one / (ptsaux(1, ip)*ptsaux(1, iq))
                zmat (1, k-np) = temp
                zmat (ip+1, k-np) = - temp
                zmat (iq+1, k-np) = - temp
70          zmat (k, k-np) = temp
        end if
        nrem = npt
        kold = 1
        knew = kopt
!
!     Reorder the provisional points in the way that exchanges PTSID(KOLD)
!     with PTSID(KNEW).
!
80      do 90 j = 1, n
            temp = bmat (kold, j)
            bmat (kold, j) = bmat (knew, j)
90      bmat (knew, j) = temp
        do 100 j = 1, nptm
            temp = zmat (kold, j)
            zmat (kold, j) = zmat (knew, j)
100     zmat (knew, j) = temp
        ptsid (kold) = ptsid (knew)
        ptsid (knew) = zero
        w (ndim+knew) = zero
        nrem = nrem - 1
        if (knew /= kopt) then
            temp = vlag (kold)
            vlag (kold) = vlag (knew)
            vlag (knew) = temp
!
!     Update the BMAT and ZMAT matrices so that the status of the KNEW-th
!     interpolation point can be changed from provisional to original. The
!     branch to label 350 occurs if all the original points are reinstated.
!     The nonnegative values of W(NDIM+K) are required in the search below.
!
            call update (n, npt, bmat, zmat, ndim, vlag, beta, denom, knew, w)
            if (nrem == 0) go to 350
            do 110 k = 1, npt
110         w (ndim+k) = abs (w(ndim+k))
        end if
!
!     Pick the index KNEW of an original interpolation point that has not
!     yet replaced one of the provisional interpolation points, giving
!     attention to the closeness to XOPT and to previous tries with KNEW.
!
120     dsqmin = zero
        do 130 k = 1, npt
            if (w(ndim+k) > zero) then
                if (dsqmin == zero .or. w(ndim+k) < dsqmin) then
                    knew = k
                    dsqmin = w (ndim+k)
                end if
            end if
130     continue
        if (dsqmin == zero) go to 260
!
!     Form the W-vector of the chosen original interpolation point.
!
        do 140 j = 1, n
140     w (npt+j) = xpt (knew, j)
        do 160 k = 1, npt
            sum = zero
            if (k == kopt) then
                continue
            else if (ptsid(k) == zero) then
                do 150 j = 1, n
150             sum = sum + w (npt+j) * xpt (k, j)
            else
                ip = ptsid (k)
                if (ip > 0) sum = w (npt+ip) * ptsaux (1, ip)
                iq = real (np, wp) * ptsid (k) - real (ip*np, wp)
                if (iq > 0) then
                    iw = 1
                    if (ip == 0) iw = 2
                    sum = sum + w (npt+iq) * ptsaux (iw, iq)
                end if
            end if
160     w (k) = half * sum * sum
!
!     Calculate VLAG and BETA for the required updating of the H matrix if
!     XPT(KNEW,.) is reinstated in the set of interpolation points.
!
        do 180 k = 1, npt
            sum = zero
            do 170 j = 1, n
170         sum = sum + bmat (k, j) * w (npt+j)
180     vlag (k) = sum
        beta = zero
        do 200 j = 1, nptm
            sum = zero
            do 190 k = 1, npt
190         sum = sum + zmat (k, j) * w (k)
            beta = beta - sum * sum
            do 200 k = 1, npt
200     vlag (k) = vlag (k) + sum * zmat (k, j)
        bsum = zero
        distsq = zero
        do 230 j = 1, n
            sum = zero
            do 210 k = 1, npt
210         sum = sum + bmat (k, j) * w (k)
            jp = j + npt
            bsum = bsum + sum * w (jp)
            do 220 ip = npt + 1, ndim
220         sum = sum + bmat (ip, j) * w (ip)
            bsum = bsum + sum * w (jp)
            vlag (jp) = sum
230     distsq = distsq + xpt (knew, j) ** 2
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
        do 250 k = 1, npt
            if (ptsid(k) /= zero) then
                hdiag = zero
                do 240 j = 1, nptm
240             hdiag = hdiag + zmat (k, j) ** 2
                den = beta * hdiag + vlag (k) ** 2
                if (den > denom) then
                    kold = k
                    denom = den
                end if
            end if
250     vlmxsq = max (vlmxsq, vlag(k)**2)
        if (denom <= 1.0e-2_wp*vlmxsq) then
            w (ndim+knew) = - w (ndim+knew) - winc
            go to 120
        end if
        go to 80
!
!     When label 260 is reached, all the final positions of the interpolation
!     points have been chosen although any changes have not been included yet
!     in XPT. Also the final BMAT and ZMAT matrices are complete, but, apart
!     from the shift of XBASE, the updating of the quadratic model remains to
!     be done. The following cycle through the new interpolation points begins
!     by putting the new point in XPT(KPT,.) and by setting PQ(KPT) to zero,
!     except that a RETURN occurs if MAXFUN prohibits another value of F.
!
260     do 340 kpt = 1, npt
            if (ptsid(kpt) == zero) go to 340
            if (nf >= maxfun) then
                nf = - 1
                go to 350
            end if
            ih = 0
            do 270 j = 1, n
                w (j) = xpt (kpt, j)
                xpt (kpt, j) = zero
                temp = pq (kpt) * w (j)
                do 270 i = 1, j
                    ih = ih + 1
270         hq (ih) = hq (ih) + temp * w (i)
            pq (kpt) = zero
            ip = ptsid (kpt)
            iq = real (np, wp) * ptsid (kpt) - real (ip*np, wp)
            if (ip > 0) then
                xp = ptsaux (1, ip)
                xpt (kpt, ip) = xp
            end if
            if (iq > 0) then
                xq = ptsaux (1, iq)
                if (ip == 0) xq = ptsaux (2, iq)
                xpt (kpt, iq) = xq
            end if
!
!     Set VQUAD to the value of the current model at the new point.
!
            vquad = fbase
            if (ip > 0) then
                ihp = (ip+ip*ip) / 2
                vquad = vquad + xp * (gopt(ip)+half*xp*hq(ihp))
            end if
            if (iq > 0) then
                ihq = (iq+iq*iq) / 2
                vquad = vquad + xq * (gopt(iq)+half*xq*hq(ihq))
                if (ip > 0) then
                    iw = max0 (ihp, ihq) - iabs (ip-iq)
                    vquad = vquad + xp * xq * hq (iw)
                end if
            end if
            do 280 k = 1, npt
                temp = zero
                if (ip > 0) temp = temp + xp * xpt (k, ip)
                if (iq > 0) temp = temp + xq * xpt (k, iq)
280         vquad = vquad + half * pq (k) * temp * temp
!
!     Calculate F at the new interpolation point, and set DIFF to the factor
!     that is going to multiply the KPT-th Lagrange function when the model
!     is updated to provide interpolation to the new function value.
!
            do 290 i = 1, n
                w (i) = min (max(xl(i), xbase(i)+xpt(kpt, i)), xu(i))
                if (xpt(kpt, i) == sl(i)) w (i) = xl (i)
                if (xpt(kpt, i) == su(i)) w (i) = xu (i)
290         continue
            nf = nf + 1
            call calfun (n, w, f)
            if (iprint == 3) then
                print 300, nf, f, (w(i), i=1, n)
300             format (/ 4 x, 'Function number', i6, '    F =', 1 pd18.10,&
                '    The corresponding X is:' / (2 x, 5d15.6))
            end if
            fval (kpt) = f
            if (f < fval(kopt)) kopt = kpt
            diff = f - vquad
!
!     Update the quadratic model. The RETURN from the subroutine occurs when
!     all the new interpolation points are included in the model.
!
            do 310 i = 1, n
310         gopt (i) = gopt (i) + diff * bmat (kpt, i)
            do 330 k = 1, npt
                sum = zero
                do 320 j = 1, nptm
320             sum = sum + zmat (k, j) * zmat (kpt, j)
                temp = diff * sum
                if (ptsid(k) == zero) then
                    pq (k) = pq (k) + temp
                else
                    ip = ptsid (k)
                    iq = real (np, wp) * ptsid (k) - real (ip*np, wp)
                    ihq = (iq*iq+iq) / 2
                    if (ip == 0) then
                        hq (ihq) = hq (ihq) + temp * ptsaux (2, iq) ** 2
                    else
                        ihp = (ip*ip+ip) / 2
                        hq (ihp) = hq (ihp) + temp * ptsaux (1, ip) ** 2
                        if (iq > 0) then
                            hq (ihq) = hq (ihq) + temp * ptsaux (1, iq) ** 2
                            iw = max0 (ihp, ihq) - iabs (iq-ip)
                            hq (iw) = hq (iw) + temp * ptsaux (1, ip) * ptsaux (1, iq)
                        end if
                    end if
                end if
330         continue
            ptsid (kpt) = zero
340     continue
350     return
    end subroutine rescue
 
    subroutine trsbox (n, npt, xpt, xopt, gopt, hq, pq, sl, su, delta, xnew, d, gnew, &
   & xbdi, s, hs, hred, dsq, crvmin)
        implicit real (wp) (a-h, o-z)
        dimension xpt (npt,*), xopt (*), gopt (*), hq (*), pq (*), sl (*), su (*), xnew &
       & (*), d (*), gnew (*), xbdi (*), s (*), hs (*), hred (*)
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
        half = 0.5_wp
        one = 1.0_wp
        onemin = - 1.0_wp
        zero = 0.0_wp
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
        do 10 i = 1, n
            xbdi (i) = zero
            if (xopt(i) <= sl(i)) then
                if (gopt(i) >= zero) xbdi (i) = onemin
            else if (xopt(i) >= su(i)) then
                if (gopt(i) <= zero) xbdi (i) = one
            end if
            if (xbdi(i) /= zero) nact = nact + 1
            d (i) = zero
10      gnew (i) = gopt (i)
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
20      beta = zero
30      stepsq = zero
        do 40 i = 1, n
            if (xbdi(i) /= zero) then
                s (i) = zero
            else if (beta == zero) then
                s (i) = - gnew (i)
            else
                s (i) = beta * s (i) - gnew (i)
            end if
40      stepsq = stepsq + s (i) ** 2
        if (stepsq == zero) go to 190
        if (beta == zero) then
            gredsq = stepsq
            itermax = iterc + n - nact
        end if
        if (gredsq*delsq <= 1.0e-4_wp*qred*qred) go to 190
!
!     Multiply the search direction by the second derivative matrix of Q and
!     calculate some scalars for the choice of steplength. Then set BLEN to
!     the length of the the step to the trust region boundary and STPLEN to
!     the steplength, ignoring the simple bounds.
!
        go to 210
50      resid = delsq
        ds = zero
        shs = zero
        do 60 i = 1, n
            if (xbdi(i) == zero) then
                resid = resid - d (i) ** 2
                ds = ds + s (i) * d (i)
                shs = shs + s (i) * hs (i)
            end if
60      continue
        if (resid <= zero) go to 90
        temp = sqrt (stepsq*resid+ds*ds)
        if (ds < zero) then
            blen = (temp-ds) / stepsq
        else
            blen = resid / (temp+ds)
        end if
        stplen = blen
        if (shs > zero) then
            stplen = min (blen, gredsq/shs)
        end if
!
!
!     Reduce STPLEN if necessary in order to preserve the simple bounds,
!     letting IACT be the index of the new constrained variable.
!
        iact = 0
        do 70 i = 1, n
            if (s(i) /= zero) then
                xsum = xopt (i) + d (i)
                if (s(i) > zero) then
                    temp = (su(i)-xsum) / s (i)
                else
                    temp = (sl(i)-xsum) / s (i)
                end if
                if (temp < stplen) then
                    stplen = temp
                    iact = i
                end if
            end if
70      continue
!
!     Update CRVMIN, GNEW and D. Set SDEC to the decrease that occurs in Q.
!
        sdec = zero
        if (stplen > zero) then
            iterc = iterc + 1
            temp = shs / stepsq
            if (iact == 0 .and. temp > zero) then
                crvmin = min (crvmin, temp)
                if (crvmin == onemin) crvmin = temp
            end if
            ggsav = gredsq
            gredsq = zero
            do 80 i = 1, n
                gnew (i) = gnew (i) + stplen * hs (i)
                if (xbdi(i) == zero) gredsq = gredsq + gnew (i) ** 2
80          d (i) = d (i) + stplen * s (i)
            sdec = max (stplen*(ggsav-half*stplen*shs), zero)
            qred = qred + sdec
        end if
!
!     Restart the conjugate gradient method if it has hit a new bound.
!
        if (iact > 0) then
            nact = nact + 1
            xbdi (iact) = one
            if (s(iact) < zero) xbdi (iact) = onemin
            delsq = delsq - d (iact) ** 2
            if (delsq <= zero) go to 90
            go to 20
        end if
!
!     If STPLEN is less than BLEN, then either apply another conjugate
!     gradient iteration or RETURN.
!
        if (stplen < blen) then
            if (iterc == itermax) go to 190
            if (sdec <= 0.01_wp*qred) go to 190
            beta = gredsq / ggsav
            go to 30
        end if
90      crvmin = zero
!
!     Prepare for the alternative iteration by calculating some scalars and
!     by multiplying the reduced D by the second derivative matrix of Q.
!
100     if (nact >= n-1) go to 190
        dredsq = zero
        dredg = zero
        gredsq = zero
        do 110 i = 1, n
            if (xbdi(i) == zero) then
                dredsq = dredsq + d (i) ** 2
                dredg = dredg + d (i) * gnew (i)
                gredsq = gredsq + gnew (i) ** 2
                s (i) = d (i)
            else
                s (i) = zero
            end if
110     continue
        itcsav = iterc
        go to 210
!
!     Let the search direction S be a linear combination of the reduced D
!     and the reduced G that is orthogonal to the reduced D.
!
120     iterc = iterc + 1
        temp = gredsq * dredsq - dredg * dredg
        if (temp <= 1.0e-4_wp*qred*qred) go to 190
        temp = sqrt (temp)
        do 130 i = 1, n
            if (xbdi(i) == zero) then
                s (i) = (dredg*d(i)-dredsq*gnew(i)) / temp
            else
                s (i) = zero
            end if
130     continue
        sredg = - temp
!
!     By considering the simple bounds on the variables, calculate an upper
!     bound on the tangent of half the angle of the alternative iteration,
!     namely ANGBD, except that, if already a free variable has reached a
!     bound, there is a branch back to label 100 after fixing that variable.
!
        angbd = one
        iact = 0
        do 140 i = 1, n
            if (xbdi(i) == zero) then
                tempa = xopt (i) + d (i) - sl (i)
                tempb = su (i) - xopt (i) - d (i)
                if (tempa <= zero) then
                    nact = nact + 1
                    xbdi (i) = onemin
                    go to 100
                else if (tempb <= zero) then
                    nact = nact + 1
                    xbdi (i) = one
                    go to 100
                end if
                ratio = one
                ssq = d (i) ** 2 + s (i) ** 2
                temp = ssq - (xopt(i)-sl(i)) ** 2
                if (temp > zero) then
                    temp = sqrt (temp) - s (i)
                    if (angbd*temp > tempa) then
                        angbd = tempa / temp
                        iact = i
                        xsav = onemin
                    end if
                end if
                temp = ssq - (su(i)-xopt(i)) ** 2
                if (temp > zero) then
                    temp = sqrt (temp) + s (i)
                    if (angbd*temp > tempb) then
                        angbd = tempb / temp
                        iact = i
                        xsav = one
                    end if
                end if
            end if
140     continue
!
!     Calculate HHD and some curvatures for the alternative iteration.
!
        go to 210
150     shs = zero
        dhs = zero
        dhd = zero
        do 160 i = 1, n
            if (xbdi(i) == zero) then
                shs = shs + s (i) * hs (i)
                dhs = dhs + d (i) * hs (i)
                dhd = dhd + d (i) * hred (i)
            end if
160     continue
!
!     Seek the greatest reduction in Q for a range of equally spaced values
!     of ANGT in [0,ANGBD], where ANGT is the tangent of half the angle of
!     the alternative iteration.
!
        redmax = zero
        isav = 0
        redsav = zero
        iu = 17.0_wp * angbd + 3.1_wp
        do 170 i = 1, iu
            angt = angbd * real (i, wp) / real (iu, wp)
            sth = (angt+angt) / (one+angt*angt)
            temp = shs + angt * (angt*dhd-dhs-dhs)
            rednew = sth * (angt*dredg-sredg-half*sth*temp)
            if (rednew > redmax) then
                redmax = rednew
                isav = i
                rdprev = redsav
            else if (i == isav+1) then
                rdnext = rednew
            end if
170     redsav = rednew
!
!     Return if the reduction is zero. Otherwise, set the sine and cosine
!     of the angle of the alternative iteration, and calculate SDEC.
!
        if (isav == 0) go to 190
        if (isav < iu) then
            temp = (rdnext-rdprev) / (redmax+redmax-rdprev-rdnext)
            angt = angbd * (real(isav, wp)+half*temp) / real (iu, wp)
        end if
        cth = (one-angt*angt) / (one+angt*angt)
        sth = (angt+angt) / (one+angt*angt)
        temp = shs + angt * (angt*dhd-dhs-dhs)
        sdec = sth * (angt*dredg-sredg-half*sth*temp)
        if (sdec <= zero) go to 190
!
!     Update GNEW, D and HRED. If the angle of the alternative iteration
!     is restricted by a bound on a free variable, that variable is fixed
!     at the bound.
!
        dredg = zero
        gredsq = zero
        do 180 i = 1, n
            gnew (i) = gnew (i) + (cth-one) * hred (i) + sth * hs (i)
            if (xbdi(i) == zero) then
                d (i) = cth * d (i) + sth * s (i)
                dredg = dredg + d (i) * gnew (i)
                gredsq = gredsq + gnew (i) ** 2
            end if
180     hred (i) = cth * hred (i) + sth * hs (i)
        qred = qred + sdec
        if (iact > 0 .and. isav == iu) then
            nact = nact + 1
            xbdi (iact) = xsav
            go to 100
        end if
!
!     If SDEC is sufficiently small, then RETURN after setting XNEW to
!     XOPT+D, giving careful attention to the bounds.
!
        if (sdec > 0.01_wp*qred) go to 120
190     dsq = zero
        do 200 i = 1, n
            xnew (i) = max (min(xopt(i)+d(i), su(i)), sl(i))
            if (xbdi(i) == onemin) xnew (i) = sl (i)
            if (xbdi(i) == one) xnew (i) = su (i)
            d (i) = xnew (i) - xopt (i)
200     dsq = dsq + d (i) ** 2
        return
!
!     The following instructions multiply the current S-vector by the second
!     derivative matrix of the quadratic model, putting the product in HS.
!     They are reached from three different parts of the software above and
!     they can be regarded as an external subroutine.
!
210     ih = 0
        do 220 j = 1, n
            hs (j) = zero
            do 220 i = 1, j
                ih = ih + 1
                if (i < j) hs (j) = hs (j) + hq (ih) * s (i)
220     hs (i) = hs (i) + hq (ih) * s (j)
        do 250 k = 1, npt
            if (pq(k) /= zero) then
                temp = zero
                do 230 j = 1, n
230             temp = temp + xpt (k, j) * s (j)
                temp = temp * pq (k)
                do 240 i = 1, n
240             hs (i) = hs (i) + temp * xpt (k, i)
            end if
250     continue
        if (crvmin /= zero) go to 50
        if (iterc > itcsav) go to 150
        do 260 i = 1, n
260     hred (i) = hs (i)
        go to 120
    end subroutine trsbox
 
    subroutine update (n, npt, bmat, zmat, ndim, vlag, beta, denom, knew, w)
 
        implicit real (wp) (a-h, o-z)
        dimension bmat (ndim,*), zmat (npt,*), vlag (*), w (*)
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
        one = 1.0_wp
        zero = 0.0_wp
        nptm = npt - n - 1
        ztest = zero
        do 10 k = 1, npt
            do 10 j = 1, nptm
10      ztest = max (ztest, abs(zmat(k, j)))
        ztest = 1.0e-20_wp * ztest
!
!     Apply the rotations that put zeros in the KNEW-th row of ZMAT.
!
        jl = 1
        do 30 j = 2, nptm
            if (abs(zmat(knew, j)) > ztest) then
                temp = sqrt (zmat(knew, 1)**2+zmat(knew, j)**2)
                tempa = zmat (knew, 1) / temp
                tempb = zmat (knew, j) / temp
                do 20 i = 1, npt
                    temp = tempa * zmat (i, 1) + tempb * zmat (i, j)
                    zmat (i, j) = tempa * zmat (i, j) - tempb * zmat (i, 1)
20              zmat (i, 1) = temp
            end if
            zmat (knew, j) = zero
30      continue
!
!     Put the first NPT components of the KNEW-th column of HLAG into W,
!     and calculate the parameters of the updating formula.
!
        do 40 i = 1, npt
            w (i) = zmat (knew, 1) * zmat (i, 1)
40      continue
        alpha = w (knew)
        tau = vlag (knew)
        vlag (knew) = vlag (knew) - one
!
!     Complete the updating of ZMAT.
!
        temp = sqrt (denom)
        tempb = zmat (knew, 1) / temp
        tempa = tau / temp
        do 50 i = 1, npt
50      zmat (i, 1) = tempa * zmat (i, 1) - tempb * vlag (i)
!
!     Finally, update the matrix BMAT.
!
        do 60 j = 1, n
            jp = npt + j
            w (jp) = bmat (knew, j)
            tempa = (alpha*vlag(jp)-tau*w(jp)) / denom
            tempb = (-beta*w(jp)-tau*vlag(jp)) / denom
            do 60 i = 1, jp
                bmat (i, j) = bmat (i, j) + tempa * vlag (i) + tempb * w (i)
                if (i > npt) bmat (jp, i-npt) = bmat (i, j)
60      continue
        return
    end subroutine update
 
    subroutine bobyqa_test ()
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
        implicit real (wp) (a-h, o-z)
        dimension x (100), xl (100), xu (100), w (500000)
 
        twopi = 8.0_wp * atan (1.0_wp)
        bdl = - 1.0_wp
        bdu = 1.0_wp
        iprint = 2
        maxfun = 500000
        rhobeg = 1.0e-1_wp
        rhoend = 1.0e-6_wp
        m = 5
10      n = 2 * m
        do 20 i = 1, n
            xl (i) = bdl
20      xu (i) = bdu
        do 50 jcase = 1, 2
            npt = n + 6
            if (jcase == 2) npt = 2 * n + 1
            print 30, m, n, npt
30          format (/ / 5 x, '2D output with M =', i4, ',  N =', i4, '  and  NPT =', i4)
            do 40 j = 1, m
                temp = real (j, wp) * twopi / real (m, wp)
                x (2*j-1) = cos (temp)
40          x (2*j) = sin (temp)
            call bobyqa (n, npt, x, xl, xu, rhobeg, rhoend, iprint, maxfun, w, calfun)
50      continue
        m = m + m
        if (m <= 10) go to 10
 
    contains
 
        subroutine calfun (n, x, f)
            implicit real (wp) (a-h, o-z)
            dimension x (*)
            f = 0.0_wp
            do i = 4, n, 2
                do j = 2, i - 2, 2
                    temp = (x(i-1)-x(j-1)) ** 2 + (x(i)-x(j)) ** 2
                    temp = max (temp, 1.0e-6_wp)
                    f = f + 1.0_wp / sqrt (temp)
                end do
            end do
        end subroutine calfun
 
    end subroutine bobyqa_test
 
end module bobyqa_module
