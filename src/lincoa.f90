!*****************************************************************************************
!>
!  LINCOA: **LIN**early **C**onstrained **O**ptimization **A**lgorithm
!
!  The purpose of LINCOA is to seek the least value of a function F of several variables
!  subject to general linear inequality constraints on the variables,
!  when derivatives of F are not available.
!
!# History
!  * M.J.D. Powell, December 6th, 2013 : There are no
!    restrictions on or charges for the use of the software. I hope that the time
!    and effort I have spent on developing the package will be helpful to much
!    research and to many applications.
!  * Jacob Williams, July 2015 : refactoring of the code into modern Fortran.

    module lincoa_module

    use kind_module, only: wp

    private

    abstract interface
        subroutine func (n, x, f)  !! calfun interface
            import :: wp
            implicit none
            integer :: n
            real(wp) :: x (*)
            real(wp) :: f
        end subroutine func
    end interface

    public :: lincoa
    public :: lincoa_test

    contains

!*****************************************************************************************
!>
!  This subroutine seeks the least value of a function of many variables,
!  subject to general linear inequality constraints, by a trust region
!  method that forms quadratic models by interpolation.
!
!  LINCOA solves the following optimization problem:
!```
!   Minimize F(X(1),X(2),...X(N)) subject to:
!   A * X <= B
!```
!
!  Usually there
!  is much freedom in each new model after satisfying the interpolation
!  conditions, which is taken up by minimizing the Frobenius norm of
!  the change to the second derivative matrix of the model. One new
!  function value is calculated on each iteration, usually at a point
!  where the current model predicts a reduction in the least value so
!  far of the objective function subject to the linear constraints.
!  Alternatively, a new vector of variables may be chosen to replace
!  an interpolation point that may be too far away for reliability, and
!  then the new point does not have to satisfy the linear constraints.

    subroutine lincoa (n, npt, m, a, ia, b, x, rhobeg, rhoend, iprint, maxfun, calfun)

        implicit none

        integer,intent(in)                  :: n       !! the number of variables. must be at least 2.
        integer,intent(in)                  :: npt     !! the number of interpolation conditions, which is
                                                       !! required to be in the interval [N+2,(N+1)(N+2)/2]. Typical choices
                                                       !! of the author are NPT=N+6 and NPT=2*N+1. Larger values tend to be
                                                       !! highly inefficent when the number of variables is substantial, due
                                                       !! to the amount of work and extra difficulty of adjusting more points.
        integer,intent(in)                  :: m       !! the number of linear inequality constraints.
        integer,intent(in)                  :: ia      !! the first dimension of the array A, which must be at least N.
        real(wp),dimension(ia,*),intent(in) :: a       !! a matrix whose columns are the constraint gradients, which are
                                                       !! required to be nonzero.
        real(wp),dimension(*),intent(in)    :: b       !! the vector of right hand sides of the constraints, the J-th
                                                       !! constraint being that the scalar product of A(.,J) with X(.) is at
                                                       !! most B(J). The initial vector X(.) is made feasible by increasing
                                                       !! the value of B(J) if necessary.
        real(wp),dimension(*),intent(inout) :: x       !! the vector of variables. Initial values of X(1),X(2),...,X(N)
                                                       !! must be supplied. If they do not satisfy the constraints, then B
                                                       !! is increased as mentioned above. X contains on return the variables
                                                       !! that have given the least calculated F subject to the constraints.
        real(wp),intent(in)                 :: rhobeg  !! RHOBEG and RHOEND must be set to the initial and final values of a
                                                       !! trust region radius, so both must be positive with RHOEND<=RHOBEG.
                                                       !! Typically, RHOBEG should be about one tenth of the greatest expected
                                                       !! change to a variable, and RHOEND should indicate the accuracy that
                                                       !! is required in the final values of the variables.
        real(wp),intent(in)                 :: rhoend  !! RHOBEG and RHOEND must be set to the initial and final values of a
                                                       !! trust region radius, so both must be positive with RHOEND<=RHOBEG.
                                                       !! Typically, RHOBEG should be about one tenth of the greatest expected
                                                       !! change to a variable, and RHOEND should indicate the accuracy that
                                                       !! is required in the final values of the variables.
        integer,intent(in)                  :: iprint  !! The value of IPRINT should be set to 0, 1, 2 or 3, which controls the
                                                       !! amount of printing. Specifically, there is no output if IPRINT=0 and
                                                       !! there is output only at the return if IPRINT=1. Otherwise, the best
                                                       !! feasible vector of variables so far and the corresponding value of
                                                       !! the objective function are printed whenever RHO is reduced, where
                                                       !! RHO is the current lower bound on the trust region radius. Further,
                                                       !! each new value of F with its variables are output if IPRINT=3.
        integer,intent(in)                  :: maxfun  !! an upper bound on the number of calls of CALFUN,
                                                       !! its value being at least NPT+1.
        procedure (func)                    :: calfun  !! It must set
                                                       !! F to the value of the objective function for the variables X(1),
                                                       !! X(2),...,X(N). The value of the argument F is positive when CALFUN
                                                       !! is called if and only if the current X satisfies the constraints

        real(wp),parameter :: zero = 0.0_wp

        real(wp),dimension(:),allocatable :: w
        integer, dimension(n) :: iact !to avoid type mismatch error - JW
        real(wp) :: smallx,sum,temp
        integer :: np,nptm,iamat,ib,iflag,i,iac,ibmat,ifv,igo,ihq,ipq,ipqw,&
                   iqf,irc,isp,istp,iw,ixb,ixn,ixo,ixp,ixs,irf,izmat,j,ndim

!     W is an array used for working space. Its length must be at least
!       M*(2+N) + NPT*(4+N+NPT) + N*(9+3*N) + MAX [ M+3*N, 2*M+N, 2*NPT ].
!       On return, W(1) is set to the final value of F, and W(2) is set to
!       the total number of function evaluations plus 0.5.
        allocate(w(M*(2+N) + NPT*(4+N+NPT) + N*(9+3*N) + MAX(M+3*N, 2*M+N, 2*NPT)))

!
!     Check that N, NPT and MAXFUN are acceptable.
!
        smallx = 1.0e-6_wp * rhoend
        np = n + 1
        nptm = npt - np
        if (n <= 1) then
            print 10
10          format (/ 4 x, 'Return from LINCOA because N is less than 2.')
            return
        end if
        if (npt < n+2 .or. npt > ((n+2)*np)/2) then
            print 20
20          format (/ 4 x, 'Return from LINCOA because NPT is not in',&
            ' the required interval.')
            return
        end if
        if (maxfun <= npt) then
            print 30
30          format (/ 4 x, 'Return from LINCOA because MAXFUN is less', ' than NPT+1.')
            return
        end if
!
!     Normalize the constraints, and copy the resultant constraint matrix
!       and right hand sides into working space, after increasing the right
!       hand sides if necessary so that the starting point is feasible.
!
        iamat = max (m+3*n, 2*m+n, 2*npt) + 1
        ib = iamat + m * n
        iflag = 0
        if (m > 0) then
            iw = iamat - 1
            do j = 1, m
                sum = zero
                temp = zero
                do i = 1, n
                    sum = sum + a (i, j) * x (i)
                    temp = temp + a (i, j) ** 2
                end do
                if (temp == zero) then
                    print 50
50                  format (/ 4 x, 'Return from LINCOA because the gradient of',&
                    ' a constraint is zero.')
                    return
                end if
                temp = sqrt (temp)
                if (sum-b(j) > smallx*temp) iflag = 1
                w (ib+j-1) = max (b(j), sum) / temp
                do i = 1, n
                    iw = iw + 1
                    w (iw) = a (i, j) / temp
                end do
            end do
        end if
        if (iflag == 1) then
            if (iprint > 0) print 70
70          format (/ 4 x, 'LINCOA has made the initial X feasible by',&
            ' increasing part(s) of B.')
        end if
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
        call lincob (n, npt, m, w(iamat), w(ib), x, rhobeg, rhoend, iprint, maxfun, &
                     w(ixb), w(ixp), w(ifv), w(ixs), w(ixo), w(igo), w(ihq), w(ipq), &
                     w(ibmat), w(izmat), ndim, w(istp), w(isp), w(ixn), iact, w(irc), &
                     w(iqf), w(irf), w(ipqw), w, calfun)

        deallocate(w)

    end subroutine lincoa
!*****************************************************************************************

    subroutine lincob (n, npt, m, amat, b, x, rhobeg, rhoend, iprint, maxfun, xbase, xpt, &
                       fval, xsav, xopt, gopt, hq, pq, bmat, zmat, ndim, step, sp, xnew, &
                       iact, rescon, qfac, rfac, pqw, w, calfun)

        implicit real (wp) (a-h, o-z)

        dimension amat (n,*), b (*), x (*), xbase (*), xpt (npt,*), fval (*), xsav (*), &
                  xopt (*), gopt (*), hq (*), pq (*), bmat (ndim,*), zmat (npt,*), &
                  step (*), sp (*), xnew (*), iact (*), rescon (*), qfac (n,*), rfac (*), &
                  pqw (*), w (*)
        procedure (func) :: calfun
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
        real(wp),parameter :: half  = 0.5_wp
        real(wp),parameter :: one   = 1.0_wp
        real(wp),parameter :: tenth = 0.1_wp
        real(wp),parameter :: zero  = 0.0_wp

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
        call prelim (n, npt, m, amat, b, x, rhobeg, iprint, xbase, xpt, fval, xsav, xopt, &
                     gopt, kopt, hq, pq, bmat, zmat, idz, ndim, sp, rescon, step, pqw, w, &
                     calfun)
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
10      knew = 0
        nvala = 0
        nvalb = 0
!
!     Shift XBASE if XOPT may be too far from XBASE. First make the changes
!       to BMAT that do not depend on ZMAT.
!
20      fsave = fopt
        xoptsq = zero
        do i = 1, n
            xoptsq = xoptsq + xopt (i) ** 2
        end do
        if (xoptsq >= 1.0e4_wp*delta*delta) then
            qoptsq = 0.25_wp * xoptsq
            do k = 1, npt
                sum = zero
                do i = 1, n
                    sum = sum + xpt (k, i) * xopt (i)
                end do
                sum = sum - half * xoptsq
                w (npt+k) = sum
                sp (k) = zero
                do i = 1, n
                    xpt (k, i) = xpt (k, i) - half * xopt (i)
                    step (i) = bmat (k, i)
                    w (i) = sum * xpt (k, i) + qoptsq * xopt (i)
                    ip = npt + i
                    do j = 1, i
                        bmat (ip, j) = bmat (ip, j) + step (i) * w (j) + w (i) * step (j)
                    end do
                end do
            end do
!
!     Then the revisions of BMAT that depend on ZMAT are calculated.
!
            do k = 1, nptm
                sumz = zero
                do i = 1, npt
                    sumz = sumz + zmat (i, k)
                    w (i) = w (npt+i) * zmat (i, k)
                end do
                do j = 1, n
                    sum = qoptsq * sumz * xopt (j)
                    do i = 1, npt
                        sum = sum + w (i) * xpt (i, j)
                    end do
                    step (j) = sum
                    if (k < idz) sum = - sum
                    do i = 1, npt
                        bmat (i, j) = bmat (i, j) + sum * zmat (i, k)
                    end do
                end do
                do i = 1, n
                    ip = i + npt
                    temp = step (i)
                    if (k < idz) temp = - temp
                    do j = 1, i
                        bmat (ip, j) = bmat (ip, j) + temp * step (j)
                    end do
                end do
            end do
!
!     Update the right hand sides of the constraints.
!
            if (m > 0) then
                do j = 1, m
                    temp = zero
                    do i = 1, n
                        temp = temp + amat (i, j) * xopt (i)
                    end do
                    b (j) = b (j) - temp
                end do
            end if
!
!     The following instructions complete the shift of XBASE, including the
!       changes to the parameters of the quadratic model.
!
            ih = 0
            do j = 1, n
                w (j) = zero
                do k = 1, npt
                    w (j) = w (j) + pq (k) * xpt (k, j)
                    xpt (k, j) = xpt (k, j) - half * xopt (j)
                end do
                do i = 1, j
                    ih = ih + 1
                    hq (ih) = hq (ih) + w (i) * xopt (j) + xopt (i) * w (j)
                    bmat (npt+i, j) = bmat (npt+j, i)
                end do
            end do
            do j = 1, n
                xbase (j) = xbase (j) + xopt (j)
                xopt (j) = zero
                xpt (kopt, j) = zero
            end do
        end if
!
!     In the case KNEW=0, generate the next trust region step by calling
!       TRSTEP, where SNORM is the current trust region radius initially.
!       The final value of SNORM is the length of the calculated step,
!       except that SNORM is zero on return if the projected gradient is
!       unsuitable for starting the conjugate gradient iterations.
!
        delsav = delta
        ksave = knew
        if (knew == 0) then
            snorm = delta
            do i = 1, n
                xnew (i) = gopt (i)
            end do
            call trstep (n, npt, m, amat, b, xpt, hq, pq, nact, iact, rescon, qfac, rfac, &
                         snorm, step, xnew, w, w(m+1), pqw, pqw(np), w(m+np))
!
!     A trust region step is applied whenever its length, namely SNORM, is at
!       least HALF*DELTA. It is also applied if its length is at least 0.1999
!       times DELTA and if a line search of TRSTEP has caused a change to the
!       active set. Otherwise there is a branch below to label 530 or 560.
!
            temp = half * delta
            if (xnew(1) >= half) temp = 0.1999_wp * delta
            if (snorm <= temp) then
                delta = half * delta
                if (delta <= 1.4_wp*rho) delta = rho
                nvala = nvala + 1
                nvalb = nvalb + 1
                temp = snorm / rho
                if (delsav > rho) temp = one
                if (temp >= half) nvala = zero
                if (temp >= tenth) nvalb = zero
                if (delsav > rho) go to 530
                if (nvala < 5 .and. nvalb < 3) go to 530
                if (snorm > zero) ksave = - 1
                go to 560
            end if
            nvala = zero
            nvalb = zero
!
!     Alternatively, KNEW is positive. Then the model step is calculated
!       within a trust region of radius DEL, after setting the gradient at
!       XBASE and the second derivative parameters of the KNEW-th Lagrange
!       function in W(1) to W(N) and in PQW(1) to PQW(NPT), respectively.
!
        else
            del = max (tenth*delta, rho)
            do i = 1, n
                w (i) = bmat (knew, i)
            end do
            do k = 1, npt
                pqw (k) = zero
            end do
            do j = 1, nptm
                temp = zmat (knew, j)
                if (j < idz) temp = - temp
                do k = 1, npt
                    pqw (k) = pqw (k) + temp * zmat (k, j)
                end do
            end do
            call qmstep (n, npt, m, amat, b, xpt, xopt, nact, iact, rescon, qfac, kopt, &
                         knew, del, step, w, pqw, w(np), w(np+m), ifeas)
        end if
!
!     Set VQUAD to the change to the quadratic model when the move STEP is
!       made from XOPT. If STEP is a trust region step, then VQUAD should be
!       negative. If it is nonnegative due to rounding errors in this case,
!       there is a branch to label 530 to try to improve the model.
!
        vquad = zero
        ih = 0
        do j = 1, n
            vquad = vquad + step (j) * gopt (j)
            do i = 1, j
                ih = ih + 1
                temp = step (i) * step (j)
                if (i == j) temp = half * temp
                vquad = vquad + temp * hq (ih)
            end do
        end do
        do k = 1, npt
            temp = zero
            do j = 1, n
                temp = temp + xpt (k, j) * step (j)
                sp (npt+k) = temp
            end do
            vquad = vquad + half * pq (k) * temp * temp
        end do
        if (ksave == 0 .and. vquad >= zero) go to 530
!
!     Calculate the next value of the objective function. The difference
!       between the actual new value of F and the value predicted by the
!       model is recorded in DIFF.
!
220     nf = nf + 1
        if (nf > maxfun) then
            nf = nf - 1
            if (iprint > 0) print 230
230         format (/ 4 x, 'Return from LINCOA because CALFUN has been',&
            ' called MAXFUN times.')
            go to 600
        end if
        xdiff = zero
        do i = 1, n
            xnew (i) = xopt (i) + step (i)
            x (i) = xbase (i) + xnew (i)
            xdiff = xdiff + (x(i)-xsav(i)) ** 2
        end do
        xdiff = sqrt (xdiff)
        if (ksave ==-1) xdiff = rho
        if (xdiff <= tenth*rho .or. xdiff >= delta+delta) then
            ifeas = 0
            if (iprint > 0) print 250
250         format (/ 4 x, 'Return from LINCOA because rounding errors',&
            ' prevent reasonable changes to X.')
            go to 600
        end if
        if (ksave <= 0) ifeas = 1
        f = real (ifeas, wp)
        call calfun (n, x, f)
        if (iprint == 3) then
            print 260, nf, f, (x(i), i=1, n)
260         format (/ 4 x, 'Function number', i6, '    F =', 1 pd18.10,&
            '    The corresponding X is:' / (2 x, 5d15.6))
        end if
        if (ksave ==-1) go to 600
        diff = f - fopt - vquad
!
!     If X is feasible, then set DFFALT to the difference between the new
!       value of F and the value predicted by the alternative model.
!
        if (ifeas == 1 .and. itest < 3) then
            do k = 1, npt
                pqw (k) = zero
                w (k) = fval (k) - fval (kopt)
            end do
            do j = 1, nptm
                sum = zero
                do i = 1, npt
                    sum = sum + w (i) * zmat (i, j)
                end do
                if (j < idz) sum = - sum
                do k = 1, npt
                    pqw (k) = pqw (k) + sum * zmat (k, j)
                end do
            end do
            vqalt = zero
            do k = 1, npt
                sum = zero
                do j = 1, n
                    sum = sum + bmat (k, j) * step (j)
                end do
                vqalt = vqalt + sum * w (k)
                vqalt = vqalt + pqw (k) * sp (npt+k) * (half*sp(npt+k)+sp(k))
            end do
            dffalt = f - fopt - vqalt
        end if
        if (itest == 3) then
            dffalt = diff
            itest = 0
        end if
!
!     Pick the next value of DELTA after a trust region step.
!
        if (ksave == 0) then
            ratio = (f-fopt) / vquad
            if (ratio <= tenth) then
                delta = half * delta
            else if (ratio <= 0.7_wp) then
                delta = max (half*delta, snorm)
            else
                temp = sqrt (2.0_wp) * delta
                delta = max (half*delta, snorm+snorm)
                delta = min (delta, temp)
            end if
            if (delta <= 1.4_wp*rho) delta = rho
        end if
!
!     Update BMAT, ZMAT and IDZ, so that the KNEW-th interpolation point
!       can be moved. If STEP is a trust region step, then KNEW is zero at
!       present, but a positive value is picked by subroutine UPDATE.
!
        call update (n, npt, xpt, bmat, zmat, idz, ndim, sp, step, kopt, knew, pqw, w)
        if (knew == 0) then
            if (iprint > 0) print 320
320         format (/ 4 x, &
            'Return from LINCOA because the denominator of the updating formula is zero.')
            go to 600
        end if
!
!     If ITEST is increased to 3, then the next quadratic model is the
!       one whose second derivative matrix is least subject to the new
!       interpolation conditions. Otherwise the new model is constructed
!       by the symmetric Broyden method in the usual way.
!
        if (ifeas == 1) then
            itest = itest + 1
            if (abs(dffalt) >= tenth*abs(diff)) itest = 0
        end if
!
!     Update the second derivatives of the model by the symmetric Broyden
!       method, using PQW for the second derivative parameters of the new
!       KNEW-th Lagrange function. The contribution from the old parameter
!       PQ(KNEW) is included in the second derivative matrix HQ. W is used
!       later for the gradient of the new KNEW-th Lagrange function.
!
        if (itest < 3) then
            do k = 1, npt
                pqw (k) = zero
            end do
            do j = 1, nptm
                temp = zmat (knew, j)
                if (temp /= zero) then
                    if (j < idz) temp = - temp
                    do k = 1, npt
                        pqw (k) = pqw (k) + temp * zmat (k, j)
                    end do
                end if
            end do
            ih = 0
            do i = 1, n
                w (i) = bmat (knew, i)
                temp = pq (knew) * xpt (knew, i)
                do j = 1, i
                    ih = ih + 1
                    hq (ih) = hq (ih) + temp * xpt (knew, j)
                end do
            end do
            pq (knew) = zero
            do k = 1, npt
                pq (k) = pq (k) + diff * pqw (k)
            end do
        end if
!
!     Include the new interpolation point with the corresponding updates of
!       SP. Also make the changes of the symmetric Broyden method to GOPT at
!       the old XOPT if ITEST is less than 3.
!
        fval (knew) = f
        sp (knew) = sp (kopt) + sp (npt+kopt)
        ssq = zero
        do i = 1, n
            xpt (knew, i) = xnew (i)
            ssq = ssq + step (i) ** 2
        end do
        sp (npt+knew) = sp (npt+kopt) + ssq
        if (itest < 3) then
            do k = 1, npt
                temp = pqw (k) * sp (k)
                do i = 1, n
                    w (i) = w (i) + temp * xpt (k, i)
                end do
            end do
            do i = 1, n
                gopt (i) = gopt (i) + diff * w (i)
            end do
        end if
!
!     Update FOPT, XSAV, XOPT, KOPT, RESCON and SP if the new F is the
!       least calculated value so far with a feasible vector of variables.
!
        if (f < fopt .and. ifeas == 1) then
            fopt = f
            do j = 1, n
                xsav (j) = x (j)
                xopt (j) = xnew (j)
            end do
            kopt = knew
            snorm = sqrt (ssq)
            do j = 1, m
                if (rescon(j) >= delta+snorm) then
                    rescon (j) = snorm - rescon (j)
                else
                    rescon (j) = rescon (j) + snorm
                    if (rescon(j)+delta > zero) then
                        temp = b (j)
                        do i = 1, n
                            temp = temp - xopt (i) * amat (i, j)
                        end do
                        temp = max (temp, zero)
                        if (temp >= delta) temp = - temp
                        rescon (j) = temp
                    end if
                end if
            end do
            do k = 1, npt
                sp (k) = sp (k) + sp (npt+k)
            end do
!
!     Also revise GOPT when symmetric Broyden updating is applied.
!
            if (itest < 3) then
                ih = 0
                do j = 1, n
                    do i = 1, j
                        ih = ih + 1
                        if (i < j) gopt (j) = gopt (j) + hq (ih) * step (i)
                        gopt (i) = gopt (i) + hq (ih) * step (j)
                    end do
                end do
                do k = 1, npt
                    temp = pq (k) * sp (npt+k)
                    do i = 1, n
                        gopt (i) = gopt (i) + temp * xpt (k, i)
                    end do
                end do
            end if
        end if
!
!     Replace the current model by the least Frobenius norm interpolant if
!       this interpolant gives substantial reductions in the predictions
!       of values of F at feasible points.
!
        if (itest == 3) then
            do k = 1, npt
                pq (k) = zero
                w (k) = fval (k) - fval (kopt)
            end do
            do j = 1, nptm
                sum = zero
                do i = 1, npt
                    sum = sum + w (i) * zmat (i, j)
                end do
                if (j < idz) sum = - sum
                do k = 1, npt
                    pq (k) = pq (k) + sum * zmat (k, j)
                end do
            end do
            do j = 1, n
                gopt (j) = zero
                do i = 1, npt
                    gopt (j) = gopt (j) + w (i) * bmat (i, j)
                end do
            end do
            do k = 1, npt
                temp = pq (k) * sp (k)
                do i = 1, n
                    gopt (i) = gopt (i) + temp * xpt (k, i)
                end do
            end do
            do ih = 1, nh
                hq (ih) = zero
            end do
        end if
!
!     If a trust region step has provided a sufficient decrease in F, then
!       branch for another trust region calculation. Every iteration that
!       takes a model step is followed by an attempt to take a trust region
!       step.
!
        knew = 0
        if (ksave > 0) go to 20
        if (ratio >= tenth) go to 20
!
!     Alternatively, find out if the interpolation points are close enough
!       to the best point so far.
!
530     distsq = max (delta*delta, 4.0_wp*rho*rho)
        do k = 1, npt
            sum = zero
            do j = 1, n
                sum = sum + (xpt(k, j)-xopt(j)) ** 2
            end do
            if (sum > distsq) then
                knew = k
                distsq = sum
            end if
        end do
!
!     If KNEW is positive, then branch back for the next iteration, which
!       will generate a "model step". Otherwise, if the current iteration
!       has reduced F, or if DELTA was above its lower bound when the last
!       trust region step was calculated, then try a "trust region" step
!       instead.
!
        if (knew > 0) go to 20
        knew = 0
        if (fopt < fsave) go to 20
        if (delsav > rho) go to 20
!
!     The calculations with the current value of RHO are complete.
!       Pick the next value of RHO.
!
560     if (rho > rhoend) then
            delta = half * rho
            if (rho > 250.0_wp*rhoend) then
                rho = tenth * rho
            else if (rho <= 16.0_wp*rhoend) then
                rho = rhoend
            else
                rho = sqrt (rho*rhoend)
            end if
            delta = max (delta, rho)
            if (iprint >= 2) then
                if (iprint >= 3) print 570
570             format (5 x)
                print 580, rho, nf
580             format (/ 4 x, 'New RHO =', 1 pd11.4, 5 x, 'Number of',&
                ' function values =', i6)
                print 590, fopt, (xbase(i)+xopt(i), i=1, n)
590             format (4 x, 'Least value of F =', 1 pd23.15, 9 x,&
                'The corresponding X is:'/(2 x, 5d15.6))
            end if
            go to 10
        end if
!
!     Return from the calculation, after branching to label 220 for another
!       Newton-Raphson step if it has not been tried before.
!
        if (ksave ==-1) go to 220
600     if (fopt <= f .or. ifeas == 0) then
            do i = 1, n
                x (i) = xsav (i)
            end do
            f = fopt
        end if
        if (iprint >= 1) then
            print 620, nf
620         format (/ 4 x, 'At the return from LINCOA', 5 x,&
            'Number of function values =', i6)
            print 590, f, (x(i), i=1, n)
        end if
        w (1) = f
        w (2) = real (nf, wp) + half

    end subroutine lincob

    subroutine getact (n, m, amat, b, nact, iact, qfac, rfac, snorm, resnew, resact, g, &
                       dw, vlam, w)

        implicit real (wp) (a-h, o-z)

        dimension amat (n,*), b (*), iact (*), qfac (n,*), rfac (*), resnew (*), &
                  resact (*), g (*), dw (*), vlam (*), w (*)
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
        real(wp),parameter :: one  = 1.0_wp
        real(wp),parameter :: tiny = 1.0e-60_wp
        real(wp),parameter :: zero = 0.0_wp

        tdel = 0.2_wp * snorm
        ddsav = zero
        do i = 1, n
            ddsav = ddsav + g (i) ** 2
            vlam (i) = zero
        end do
        ddsav = ddsav + ddsav
!
!     Set the initial QFAC to the identity matrix in the case NACT=0.
!
        if (nact == 0) then
            do i = 1, n
                do j = 1, n
                    qfac (i, j) = zero
                end do
                qfac (i, i) = one
            end do
            go to 100
        end if
!
!     Remove any constraints from the initial active set whose residuals
!       exceed TDEL.
!
        iflag = 1
        ic = nact
40      if (resact(ic) > tdel) go to 800
50      ic = ic - 1
        if (ic > 0) go to 40
!
!     Remove any constraints from the initial active set whose Lagrange
!       multipliers are nonnegative, and set the surviving multipliers.
!
        iflag = 2
60      if (nact == 0) go to 100
        ic = nact
70      temp = zero
        do i = 1, n
            temp = temp + qfac (i, ic) * g (i)
        end do
        idiag = (ic*ic+ic) / 2
        if (ic < nact) then
            jw = idiag + ic
            do j = ic + 1, nact
                temp = temp - rfac (jw) * vlam (j)
                jw = jw + j
            end do
        end if
        if (temp >= zero) go to 800
        vlam (ic) = temp / rfac (idiag)
        ic = ic - 1
        if (ic > 0) go to 70
!
!     Set the new search direction D. Terminate if the 2-norm of D is zero
!       or does not decrease, or if NACT=N holds. The situation NACT=N
!       occurs for sufficiently large SNORM if the origin is in the convex
!       hull of the constraint gradients.
!
100     if (nact == n) go to 290
        do j = nact + 1, n
            w (j) = zero
            do i = 1, n
                w (j) = w (j) + qfac (i, j) * g (i)
            end do
        end do
        dd = zero
        do i = 1, n
            dw (i) = zero
            do j = nact + 1, n
                dw (i) = dw (i) - w (j) * qfac (i, j)
            end do
            dd = dd + dw (i) ** 2
        end do
        if (dd >= ddsav) go to 290
        if (dd == zero) go to 300
        ddsav = dd
        dnorm = sqrt (dd)
!
!     Pick the next integer L or terminate, a positive value of L being
!       the index of the most violated constraint. The purpose of CTOL
!       below is to estimate whether a positive value of VIOLMX may be
!       due to computer rounding errors.
!
        l = 0
        if (m > 0) then
            test = dnorm / snorm
            violmx = zero
            do j = 1, m
                if (resnew(j) > zero .and. resnew(j) <= tdel) then
                    sum = zero
                    do i = 1, n
                        sum = sum + amat (i, j) * dw (i)
                    end do
                    if (sum > test*resnew(j)) then
                        if (sum > violmx) then
                            l = j
                            violmx = sum
                        end if
                    end if
                end if
            end do
            ctol = zero
            temp = 0.01_wp * dnorm
            if (violmx > zero .and. violmx < temp) then
                if (nact > 0) then
                    do k = 1, nact
                        j = iact (k)
                        sum = zero
                        do i = 1, n
                            sum = sum + dw (i) * amat (i, j)
                        end do
                        ctol = max (ctol, abs(sum))
                     end do
                end if
            end if
        end if
        w (1) = one
        if (l == 0) go to 300
        if (violmx <= 10.0_wp*ctol) go to 300
!
!     Apply Givens rotations to the last (N-NACT) columns of QFAC so that
!       the first (NACT+1) columns of QFAC are the ones required for the
!       addition of the L-th constraint, and add the appropriate column
!       to RFAC.
!
        nactp = nact + 1
        idiag = (nactp*nactp-nactp) / 2
        rdiag = zero
        do j = n, 1, - 1
            sprod = zero
            do i = 1, n
                sprod = sprod + qfac (i, j) * amat (i, l)
            end do
            if (j <= nact) then
                rfac (idiag+j) = sprod
            else
                if (abs(rdiag) <= 1.0e-20_wp*abs(sprod)) then
                    rdiag = sprod
                else
                    temp = sqrt (sprod*sprod+rdiag*rdiag)
                    cosv = sprod / temp
                    sinv = rdiag / temp
                    rdiag = temp
                    do i = 1, n
                        temp = cosv * qfac (i, j) + sinv * qfac (i, j+1)
                        qfac (i, j+1) = - sinv * qfac (i, j) + cosv * qfac (i, j+1)
                        qfac (i, j) = temp
                    end do
                end if
            end if
        end do
        if (rdiag < zero) then
            do i = 1, n
                qfac (i, nactp) = - qfac (i, nactp)
            end do
        end if
        rfac (idiag+nactp) = abs (rdiag)
        nact = nactp
        iact (nact) = l
        resact (nact) = resnew (l)
        vlam (nact) = zero
        resnew (l) = zero
!
!     Set the components of the vector VMU in W.
!
220     w (nact) = one / rfac ((nact*nact+nact)/2) ** 2
        if (nact > 1) then
            do i = nact - 1, 1, - 1
                idiag = (i*i+i) / 2
                jw = idiag + i
                sum = zero
                do j = i + 1, nact
                    sum = sum - rfac (jw) * w (j)
                    jw = jw + j
                end do
                w (i) = sum / rfac (idiag)
             end do
        end if
!
!     Calculate the multiple of VMU to subtract from VLAM, and update VLAM.
!
        vmult = violmx
        ic = 0
        j = 1
250     if (j < nact) then
            if (vlam(j) >= vmult*w(j)) then
                ic = j
                vmult = vlam (j) / w (j)
            end if
            j = j + 1
            go to 250
        end if
        do j = 1, nact
            vlam (j) = vlam (j) - vmult * w (j)
        end do
        if (ic > 0) vlam (ic) = zero
        violmx = max (violmx-vmult, zero)
        if (ic == 0) violmx = zero
!
!     Reduce the active set if necessary, so that all components of the
!       new VLAM are negative, with resetting of the residuals of the
!       constraints that become inactive.
!
        iflag = 3
        ic = nact
270     if (vlam(ic) < zero) go to 280
        resnew (iact(ic)) = max (resact(ic), tiny)
        go to 800
280     ic = ic - 1
        if (ic > 0) go to 270
!
!     Calculate the next VMU if VIOLMX is positive. Return if NACT=N holds,
!       as then the active constraints imply D=0. Otherwise, go to label
!       100, to calculate the new D and to test for termination.
!
        if (violmx > zero) go to 220
        if (nact < n) go to 100
290     dd = zero
300     w (1) = dd
        return
!
!     These instructions rearrange the active constraints so that the new
!       value of IACT(NACT) is the old value of IACT(IC). A sequence of
!       Givens rotations is applied to the current QFAC and RFAC. Then NACT
!       is reduced by one.
!
800     resnew (iact(ic)) = max (resact(ic), tiny)
        jc = ic
810     if (jc < nact) then
            jcp = jc + 1
            idiag = jc * jcp / 2
            jw = idiag + jcp
            temp = sqrt (rfac(jw-1)**2+rfac(jw)**2)
            cval = rfac (jw) / temp
            sval = rfac (jw-1) / temp
            rfac (jw-1) = sval * rfac (idiag)
            rfac (jw) = cval * rfac (idiag)
            rfac (idiag) = temp
            if (jcp < nact) then
                do j = jcp + 1, nact
                    temp = sval * rfac (jw+jc) + cval * rfac (jw+jcp)
                    rfac (jw+jcp) = cval * rfac (jw+jc) - sval * rfac (jw+jcp)
                    rfac (jw+jc) = temp
                    jw = jw + j
                end do
            end if
            jdiag = idiag - jc
            do i = 1, n
                if (i < jc) then
                    temp = rfac (idiag+i)
                    rfac (idiag+i) = rfac (jdiag+i)
                    rfac (jdiag+i) = temp
                end if
                temp = sval * qfac (i, jc) + cval * qfac (i, jcp)
                qfac (i, jcp) = cval * qfac (i, jc) - sval * qfac (i, jcp)
                qfac (i, jc) = temp
            end do
            iact (jc) = iact (jcp)
            resact (jc) = resact (jcp)
            vlam (jc) = vlam (jcp)
            jc = jcp
            go to 810
        end if
        nact = nact - 1
        go to (50, 60, 280), iflag

    end subroutine getact

    subroutine prelim (n, npt, m, amat, b, x, rhobeg, iprint, xbase, xpt, fval, xsav, &
     xopt, gopt, kopt, hq, pq, bmat, zmat, idz, ndim, sp, rescon, step, pqw, w, calfun)

        implicit real (wp) (a-h, o-z)

        dimension amat(n,*),b(*),x(*),xbase(*),xpt(npt,*),fval(*),xsav(*),&
                  xopt(*),gopt(*),hq(*),pq(*),bmat(ndim,*),zmat(npt,*),sp(*),rescon(*),&
                  step(*),pqw(*),w(*)
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
        real(wp),parameter :: half = 0.5_wp
        real(wp),parameter :: one  = 1.0_wp
        real(wp),parameter :: zero = 0.0_wp

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
        do j = 1, n
            xbase (j) = x (j)
            do k = 1, npt
                xpt (k, j) = zero
            end do
            do i = 1, ndim
                bmat (i, j) = zero
            end do
        end do
        do k = 1, npt
            sp (k) = zero
            do j = 1, npt - n - 1
                zmat (k, j) = zero
            end do
        end do
!
!     Set the nonzero coordinates of XPT(K,.), K=1,2,...,min[2*N+1,NPT],
!       but they may be altered later to make a constraint violation
!       sufficiently large. The initial nonzero elements of BMAT and of
!       the first min[N,NPT-N-1] columns of ZMAT are set also.
!
        do j = 1, n
            xpt (j+1, j) = rhobeg
            if (j < npt-n) then
                jp = n + j + 1
                xpt (jp, j) = - rhobeg
                bmat (j+1, j) = half / rhobeg
                bmat (jp, j) = - half / rhobeg
                zmat (1, j) = - reciq - reciq
                zmat (j+1, j) = reciq
                zmat (jp, j) = reciq
            else
                bmat (1, j) = - one / rhobeg
                bmat (j+1, j) = one / rhobeg
                bmat (npt+j, j) = - half * rhosq
            end if
        end do
!
!     Set the remaining initial nonzero elements of XPT and ZMAT when the
!       number of interpolation points exceeds 2*N+1.
!
        if (npt > 2*n+1) then
            do k = n + 1, npt - n - 1
                itemp = (k-1) / n
                ipt = k - itemp * n
                jpt = ipt + itemp
                if (jpt > n) jpt = jpt - n
                xpt (n+k+1, ipt) = rhobeg
                xpt (n+k+1, jpt) = rhobeg
                zmat (1, k) = recip
                zmat (ipt+1, k) = - recip
                zmat (jpt+1, k) = - recip
                zmat (n+k+1, k) = recip
            end do
        end if
!
!     Update the constraint right hand sides to allow for the shift XBASE.
!
        if (m > 0) then
            do j = 1, m
                temp = zero
                do i = 1, n
                    temp = temp + amat (i, j) * xbase (i)
                end do
                b (j) = b (j) - temp
            end do
        end if
!
!     Go through the initial points, shifting every infeasible point if
!       necessary so that its constraint violation is at least 0.2*RHOBEG.
!
        do nf = 1, npt
            feas = one
            bigv = zero
            j = 0
80          j = j + 1
            if (j <= m .and. nf >= 2) then
                resid = - b (j)
                do i = 1, n
                    resid = resid + xpt (nf, i) * amat (i, j)
                end do
                if (resid <= bigv) go to 80
                bigv = resid
                jsav = j
                if (resid <= test) then
                    feas = - one
                    go to 80
                end if
                feas = zero
            end if
            if (feas < zero) then
                do i = 1, n
                    step (i) = xpt (nf, i) + (test-bigv) * amat (i, jsav)
                end do
                do k = 1, npt
                    sp (npt+k) = zero
                    do j = 1, n
                        sp (npt+k) = sp (npt+k) + xpt (k, j) * step (j)
                    end do
                end do
                call update (n,npt,xpt,bmat,zmat,idz,ndim,sp,step,kbase,nf,pqw,w)
                do i = 1, n
                    xpt (nf, i) = step (i)
                end do
            end if
!
!     Calculate the objective function at the current interpolation point,
!       and set KOPT to the index of the first trust region centre.
!
            do j = 1, n
                x (j) = xbase (j) + xpt (nf, j)
            end do
            f = feas
            call calfun (n, x, f)
            if (iprint == 3) then
                print 140, nf, f, (x(i), i=1, n)
140             format (/ 4 x, 'Function number', i6, '    F =', 1 pd18.10,&
                '    The corresponding X is:' / (2 x, 5d15.6))
            end if
            if (nf == 1) then
                kopt = 1
            else if (f < fval(kopt) .and. feas > zero) then
                kopt = nf
            end if
            fval (nf) = f
        end do
!
!     Set PQ for the first quadratic model.
!
        do j = 1, nptm
            w (j) = zero
            do k = 1, npt
                w (j) = w (j) + zmat (k, j) * fval (k)
            end do
        end do
        do k = 1, npt
            pq (k) = zero
            do j = 1, nptm
                pq (k) = pq (k) + zmat (k, j) * w (j)
            end do
        end do
!
!     Set XOPT, SP, GOPT and HQ for the first quadratic model.
!
        do j = 1, n
            xopt (j) = xpt (kopt, j)
            xsav (j) = xbase (j) + xopt (j)
            gopt (j) = zero
        end do
        do k = 1, npt
            sp (k) = zero
            do j = 1, n
                sp (k) = sp (k) + xpt (k, j) * xopt (j)
            end do
            temp = pq (k) * sp (k)
            do j = 1, n
                gopt (j) = gopt (j) + fval (k) * bmat (k, j) + temp * xpt (k, j)
            end do
        end do
        do i = 1, (n*n+n) / 2
            hq (i) = zero
        end do
!
!     Set the initial elements of RESCON.
!
        do j = 1, m
            temp = b (j)
            do i = 1, n
                temp = temp - xopt (i) * amat (i, j)
            end do
            temp = max (temp, zero)
            if (temp >= rhobeg) temp = - temp
            rescon (j) = temp
        end do

    end subroutine prelim

    subroutine qmstep (n, npt, m, amat, b, xpt, xopt, nact, iact, rescon, qfac, kopt, &
                       knew, del, step, gl, pqw, rstat, w, ifeas)

        implicit real (wp) (a-h, o-z)

        dimension amat (n,*), b (*), xpt (npt,*), xopt (*), iact (*), rescon (*), &
                  qfac (n,*), step (*), gl (*), pqw (*), rstat (*), w (*)
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
        real(wp),parameter :: half  = 0.5_wp
        real(wp),parameter :: one   = 1.0_wp
        real(wp),parameter :: tenth = 0.1_wp
        real(wp),parameter :: zero  = 0.0_wp

        test = 0.2_wp * del
!
!     Replace GL by the gradient of LFUNC at the trust region centre, and
!       set the elements of RSTAT.
!
        do k = 1, npt
            temp = zero
            do j = 1, n
                temp = temp + xpt (k, j) * xopt (j)
            end do
            temp = pqw (k) * temp
            do i = 1, n
                gl (i) = gl (i) + temp * xpt (k, i)
            end do
        end do
        if (m > 0) then
            do j = 1, m
                rstat (j) = one
                if (abs(rescon(j)) >= del) rstat (j) = - one
            end do
            do k = 1, nact
                rstat (iact(k)) = zero
            end do
        end if
!
!     Find the greatest modulus of LFUNC on a line through XOPT and
!       another interpolation point within the trust region.
!
        iflag = 0
        vbig = zero
        do k = 1, npt
            if (k == kopt) cycle
            ss = zero
            sp = zero
            do i = 1, n
                temp = xpt (k, i) - xopt (i)
                ss = ss + temp * temp
                sp = sp + gl (i) * temp
            end do
            stp = - del / sqrt (ss)
            if (k == knew) then
                if (sp*(sp-one) < zero) stp = - stp
                vlag = abs (stp*sp) + stp * stp * abs (sp-one)
            else
                vlag = abs (stp*(one-stp)*sp)
            end if
            if (vlag > vbig) then
                ksav = k
                stpsav = stp
                vbig = vlag
            end if
       end do
!
!     Set STEP to the move that gives the greatest modulus calculated above.
!       This move may be replaced by a steepest ascent step from XOPT.
!
        gg = zero
        do i = 1, n
            gg = gg + gl (i) ** 2
            step (i) = stpsav * (xpt(ksav, i)-xopt(i))
        end do
        vgrad = del * sqrt (gg)
        if (vgrad <= tenth*vbig) go to 220
!
!     Make the replacement if it provides a larger value of VBIG.
!
        ghg = zero
        do k = 1, npt
            temp = zero
            do j = 1, n
                temp = temp + xpt (k, j) * gl (j)
            end do
            ghg = ghg + pqw (k) * temp * temp
        end do
        vnew = vgrad + abs (half*del*del*ghg/gg)
        if (vnew > vbig) then
            vbig = vnew
            stp = del / sqrt (gg)
            if (ghg < zero) stp = - stp
            do i = 1, n
                step (i) = stp * gl (i)
            end do
        end if
        if (nact == 0 .or. nact == n) go to 220
!
!     Overwrite GL by its projection. Then set VNEW to the greatest
!       value of |LFUNC| on the projected gradient from XOPT subject to
!       the trust region bound. If VNEW is sufficiently large, then STEP
!       may be changed to a move along the projected gradient.
!
        do k = nact + 1, n
            w (k) = zero
            do i = 1, n
                w (k) = w (k) + gl (i) * qfac (i, k)
            end do
        end do
        gg = zero
        do i = 1, n
            gl (i) = zero
            do k = nact + 1, n
                gl (i) = gl (i) + qfac (i, k) * w (k)
            end do
            gg = gg + gl (i) ** 2
        end do
        vgrad = del * sqrt (gg)
        if (vgrad <= tenth*vbig) go to 220
        ghg = zero
        do k = 1, npt
            temp = zero
            do j = 1, n
                temp = temp + xpt (k, j) * gl (j)
            end do
            ghg = ghg + pqw (k) * temp * temp
        end do
        vnew = vgrad + abs (half*del*del*ghg/gg)
!
!     Set W to the possible move along the projected gradient.
!
        stp = del / sqrt (gg)
        if (ghg < zero) stp = - stp
        ww = zero
        do i = 1, n
            w (i) = stp * gl (i)
            ww = ww + w (i) ** 2
        end do
!
!     Set STEP to W if W gives a sufficiently large value of the modulus
!       of the Lagrange function, and if W either preserves feasibility
!       or gives a constraint violation of at least 0.2*DEL. The purpose
!       of CTOL below is to provide a check on feasibility that includes
!       a tolerance for contributions from computer rounding errors.
!
        if (vnew/vbig >= 0.2_wp) then
            ifeas = 1
            bigv = zero
            j = 0
170         j = j + 1
            if (j <= m) then
                if (rstat(j) == one) then
                    temp = - rescon (j)
                    do i = 1, n
                        temp = temp + w (i) * amat (i, j)
                    end do
                    bigv = max (bigv, temp)
                end if
                if (bigv < test) go to 170
                ifeas = 0
            end if
            ctol = zero
            temp = 0.01_wp * sqrt (ww)
            if (bigv > zero .and. bigv < temp) then
                do k = 1, nact
                    j = iact (k)
                    sum = zero
                    do i = 1, n
                        sum = sum + w (i) * amat (i, j)
                    end do
                    ctol = max (ctol, abs(sum))
                end do
            end if
            if (bigv <= 10.0_wp*ctol .or. bigv >= test) then
                do i = 1, n
                    step (i) = w (i)
                end do
                return
            end if
        end if
!
!     Calculate the greatest constraint violation at XOPT+STEP with STEP at
!       its original value. Modify STEP if this violation is unacceptable.
!
220     ifeas = 1
        bigv = zero
        resmax = zero
        j = 0
230     j = j + 1
        if (j <= m) then
            if (rstat(j) < zero) go to 230
            temp = - rescon (j)
            do i = 1, n
                temp = temp + step (i) * amat (i, j)
            end do
            resmax = max (resmax, temp)
            if (temp < test) then
                if (temp <= bigv) go to 230
                bigv = temp
                jsav = j
                ifeas = - 1
                go to 230
            end if
            ifeas = 0
        end if
        if (ifeas ==-1) then
            do i = 1, n
                step (i) = step (i) + (test-bigv) * amat (i, jsav)
            end do
            ifeas = 0
        end if
!
!     Return the calculated STEP and the value of IFEAS.
!

    end subroutine qmstep

    subroutine trstep (n, npt, m, amat, b, xpt, hq, pq, nact, iact, rescon, qfac, rfac, &
                       snorm, step, g, resnew, resact, d, dw, w)

        implicit real (wp) (a-h, o-z)

        dimension amat (n,*), b (*), xpt (npt,*), hq (*), pq (*), iact (*), rescon (*), &
         qfac (n,*), rfac (*), step (*), g (*), resnew (*), resact (*), d (*), dw (*), w &
         (*)
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
        real(wp),parameter :: half  = 0.5_wp
        real(wp),parameter :: one   = 1.0_wp
        real(wp),parameter :: tiny  = 1.0e-60_wp
        real(wp),parameter :: zero  = 0.0_wp
        real(wp),parameter :: ctest = 0.01_wp

        snsq = snorm * snorm
!
!     Set the initial elements of RESNEW, RESACT and STEP.
!
        if (m > 0) then
            do j = 1, m
                resnew (j) = rescon (j)
                if (rescon(j) >= snorm) then
                    resnew (j) = - one
                else if (rescon(j) >= zero) then
                    resnew (j) = max (resnew(j), tiny)
                end if
            end do
            if (nact > 0) then
                do k = 1, nact
                    resact (k) = rescon (iact(k))
                    resnew (iact(k)) = zero
                end do
            end if
        end if
        do i = 1, n
            step (i) = zero
        end do
        ss = zero
        reduct = zero
        ncall = 0
!
!     GETACT picks the active set for the current STEP. It also sets DW to
!       the vector closest to -G that is orthogonal to the normals of the
!       active constraints. DW is scaled to have length 0.2*SNORM, as then
!       a move of DW from STEP is allowed by the linear constraints.
!
40      ncall = ncall + 1
        call getact (n, m, amat, b, nact, iact, qfac, rfac, snorm, resnew, resact, g, dw, &
                     w, w(n+1))
        if (w(n+1) == zero) go to 320
        scale = 0.2_wp * snorm / sqrt (w(n+1))
        do i = 1, n
            dw (i) = scale * dw (i)
        end do
!
!     If the modulus of the residual of an active constraint is substantial,
!       then set D to the shortest move from STEP to the boundaries of the
!       active constraints.
!
        resmax = zero
        if (nact > 0) then
            do k = 1, nact
                resmax = max (resmax, resact(k))
            end do
        end if
        gamma = zero
        if (resmax > 1.0e-4_wp*snorm) then
            ir = 0
            do k = 1, nact
                temp = resact (k)
                if (k >= 2) then
                    do i = 1, k - 1
                        ir = ir + 1
                        temp = temp - rfac (ir) * w (i)
                    end do
                end if
                ir = ir + 1
                w (k) = temp / rfac (ir)
            end do
            do i = 1, n
                d (i) = zero
                do k = 1, nact
                    d (i) = d (i) + w (k) * qfac (i, k)
                end do
            end do
!
!     The vector D that has just been calculated is also the shortest move
!       from STEP+DW to the boundaries of the active constraints. Set GAMMA
!       to the greatest steplength of this move that satisfies the trust
!       region bound.
!
            rhs = snsq
            ds = zero
            dd = zero
            do i = 1, n
                sum = step (i) + dw (i)
                rhs = rhs - sum * sum
                ds = ds + d (i) * sum
                dd = dd + d (i) ** 2
            end do
            if (rhs > zero) then
                temp = sqrt (ds*ds+dd*rhs)
                if (ds <= zero) then
                    gamma = (temp-ds) / dd
                else
                    gamma = rhs / (temp+ds)
                end if
            end if
!
!     Reduce the steplength GAMMA if necessary so that the move along D
!       also satisfies the linear constraints.
!
            j = 0
110         if (gamma > zero) then
                j = j + 1
                if (resnew(j) > zero) then
                    ad = zero
                    adw = zero
                    do i = 1, n
                        ad = ad + amat (i, j) * d (i)
                        adw = adw + amat (i, j) * dw (i)
                    end do
                    if (ad > zero) then
                        temp = max ((resnew(j)-adw)/ad, zero)
                        gamma = min (gamma, temp)
                    end if
                end if
                if (j < m) go to 110
            end if
            gamma = min (gamma, one)
        end if
!
!     Set the next direction for seeking a reduction in the model function
!       subject to the trust region bound and the linear constraints.
!
        if (gamma <= zero) then
            do i = 1, n
                d (i) = dw (i)
            end do
            icount = nact
        else
            do i = 1, n
                d (i) = dw (i) + gamma * d (i)
            end do
            icount = nact - 1
        end if
        alpbd = one
!
!     Set ALPHA to the steplength from STEP along D to the trust region
!       boundary. Return if the first derivative term of this step is
!       sufficiently small or if no further progress is possible.
!
150     icount = icount + 1
        rhs = snsq - ss
        if (rhs <= zero) go to 320
        dg = zero
        ds = zero
        dd = zero
        do i = 1, n
            dg = dg + d (i) * g (i)
            ds = ds + d (i) * step (i)
            dd = dd + d (i) ** 2
        end do
        if (dg >= zero) go to 320
        temp = sqrt (rhs*dd+ds*ds)
        if (ds <= zero) then
            alpha = (temp-ds) / dd
        else
            alpha = rhs / (temp+ds)
        end if
        if (-alpha*dg <= ctest*reduct) go to 320
!
!     Set DW to the change in gradient along D.
!
        ih = 0
        do j = 1, n
            dw (j) = zero
            do i = 1, j
                ih = ih + 1
                if (i < j) dw (j) = dw (j) + hq (ih) * d (i)
                dw (i) = dw (i) + hq (ih) * d (j)
            end do
        end do
        do k = 1, npt
            temp = zero
            do j = 1, n
                temp = temp + xpt (k, j) * d (j)
            end do
            temp = pq (k) * temp
            do i = 1, n
                dw (i) = dw (i) + temp * xpt (k, i)
            end do
        end do
!
!     Set DGD to the curvature of the model along D. Then reduce ALPHA if
!       necessary to the value that minimizes the model.
!
        dgd = zero
        do i = 1, n
            dgd = dgd + d (i) * dw (i)
        end do
        alpht = alpha
        if (dg+alpha*dgd > zero) then
            alpha = - dg / dgd
        end if
!
!     Make a further reduction in ALPHA if necessary to preserve feasibility,
!       and put some scalar products of D with constraint gradients in W.
!
        alphm = alpha
        jsav = 0
        if (m > 0) then
            do j = 1, m
                ad = zero
                if (resnew(j) > zero) then
                    do i = 1, n
                        ad = ad + amat (i, j) * d (i)
                    end do
                    if (alpha*ad > resnew(j)) then
                        alpha = resnew (j) / ad
                        jsav = j
                    end if
                end if
                w (j) = ad
            end do
        end if
        alpha = max (alpha, alpbd)
        alpha = min (alpha, alphm)
        if (icount == nact) alpha = min (alpha, one)
!
!     Update STEP, G, RESNEW, RESACT and REDUCT.
!
        ss = zero
        do i = 1, n
            step (i) = step (i) + alpha * d (i)
            ss = ss + step (i) ** 2
            g (i) = g (i) + alpha * dw (i)
        end do
        if (m > 0) then
            do j = 1, m
                if (resnew(j) > zero) then
                    resnew (j) = max (resnew(j)-alpha*w(j), tiny)
                end if
            end do
        end if
        if (icount == nact .and. nact > 0) then
            do k = 1, nact
                resact (k) = (one-gamma) * resact (k)
            end do
        end if
        reduct = reduct - alpha * (dg+half*alpha*dgd)
!
!     Test for termination. Branch to label 40 if there is a new active
!       constraint and if the distance from STEP to the trust region
!       boundary is at least 0.2*SNORM.
!
        if (alpha == alpht) go to 320
        temp = - alphm * (dg+half*alphm*dgd)
        if (temp <= ctest*reduct) go to 320
        if (jsav > 0) then
            if (ss <= 0.64_wp*snsq) go to 40
            go to 320
        end if
        if (icount == n) go to 320
!
!     Calculate the next search direction, which is conjugate to the
!       previous one except in the case ICOUNT=NACT.
!
        if (nact > 0) then
            do j = nact + 1, n
                w (j) = zero
                do i = 1, n
                    w (j) = w (j) + g (i) * qfac (i, j)
                end do
            end do
            do i = 1, n
                temp = zero
                do j = nact + 1, n
                    temp = temp + qfac (i, j) * w (j)
                end do
                w (n+i) = temp
            end do
        else
            do i = 1, n
                w (n+i) = g (i)
            end do
        end if
        if (icount == nact) then
            beta = zero
        else
            wgd = zero
            do i = 1, n
                wgd = wgd + w (n+i) * dw (i)
            end do
            beta = wgd / dgd
        end if
        do i = 1, n
            d (i) = - w (n+i) + beta * d (i)
        end do
        alpbd = zero
        go to 150
!
!     Return from the subroutine.
!
320     snorm = zero
        if (reduct > zero) snorm = sqrt (ss)
        g (1) = zero
        if (ncall > 1) g (1) = one

    end subroutine trstep

    subroutine update (n, npt, xpt, bmat, zmat, idz, ndim, sp, step, kopt, knew, vlag, w)

        implicit real (wp) (a-h, o-z)

        dimension xpt (npt,*), bmat (ndim,*), zmat (npt,*), sp (*), step (*), vlag (*), &
                  w (*)

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
        real(wp),parameter :: half = 0.5_wp
        real(wp),parameter :: one  = 1.0_wp
        real(wp),parameter :: zero = 0.0_wp

        nptm = npt - n - 1
!
!     Calculate VLAG and BETA for the current choice of STEP. The first NPT
!       elements of VLAG are set to the values of the Lagrange functions at
!       XPT(KOPT,.)+STEP(.). The first NPT components of W_check are held
!       in W, where W_check is defined in a paper on the updating method.
!
        do k = 1, npt
            w (k) = sp (npt+k) * (half*sp(npt+k)+sp(k))
            sum = zero
            do j = 1, n
                sum = sum + bmat (k, j) * step (j)
            end do
            vlag (k) = sum
        end do
        beta = zero
        do k = 1, nptm
            sum = zero
            do i = 1, npt
                sum = sum + zmat (i, k) * w (i)
            end do
            if (k < idz) then
                beta = beta + sum * sum
                sum = - sum
            else
                beta = beta - sum * sum
            end if
            do i = 1, npt
                vlag (i) = vlag (i) + sum * zmat (i, k)
            end do
        end do
        bsum = zero
        dx = zero
        ssq = zero
        do j = 1, n
            sum = zero
            do i = 1, npt
                sum = sum + w (i) * bmat (i, j)
            end do
            bsum = bsum + sum * step (j)
            jp = npt + j
            do k = 1, n
                sum = sum + bmat (jp, k) * step (k)
            end do
            vlag (jp) = sum
            bsum = bsum + sum * step (j)
            dx = dx + step (j) * xpt (kopt, j)
            ssq = ssq + step (j) ** 2
        end do
        beta = dx * dx + ssq * (sp(kopt)+dx+dx+half*ssq) + beta - bsum
        vlag (kopt) = vlag (kopt) + one
!
!     If KNEW is zero initially, then pick the index of the interpolation
!       point to be deleted, by maximizing the absolute value of the
!       denominator of the updating formula times a weighting factor.
!
        if (knew == 0) then
            denmax = zero
            do k = 1, npt
                hdiag = zero
                do j = 1, nptm
                    temp = one
                    if (j < idz) temp = - one
                    hdiag = hdiag + temp * zmat (k, j) ** 2
                end do
                denabs = abs (beta*hdiag+vlag(k)**2)
                distsq = zero
                do j = 1, n
                    distsq = distsq + (xpt(k, j)-xpt(kopt, j)) ** 2
                end do
                temp = denabs * distsq * distsq
                if (temp > denmax) then
                    denmax = temp
                    knew = k
                end if
            end do
        end if
!
!     Apply the rotations that put zeros in the KNEW-th row of ZMAT.
!
        jl = 1
        if (nptm >= 2) then
            do j = 2, nptm
                if (j == idz) then
                    jl = idz
                else if (zmat(knew, j) /= zero) then
                    temp = sqrt (zmat(knew, jl)**2+zmat(knew, j)**2)
                    tempa = zmat (knew, jl) / temp
                    tempb = zmat (knew, j) / temp
                    do i = 1, npt
                        temp = tempa * zmat (i, jl) + tempb * zmat (i, j)
                        zmat (i, j) = tempa * zmat (i, j) - tempb * zmat (i, jl)
                        zmat (i, jl) = temp
                    end do
                    zmat (knew, j) = zero
                end if
            end do
        end if
!
!     Put the first NPT components of the KNEW-th column of the Z Z^T matrix
!       into W, and calculate the parameters of the updating formula.
!
        tempa = zmat (knew, 1)
        if (idz >= 2) tempa = - tempa
        if (jl > 1) tempb = zmat (knew, jl)
        do i = 1, npt
            w (i) = tempa * zmat (i, 1)
            if (jl > 1) w (i) = w (i) + tempb * zmat (i, jl)
        end do
        alpha = w (knew)
        tau = vlag (knew)
        tausq = tau * tau
        denom = alpha * beta + tausq
        vlag (knew) = vlag (knew) - one
        if (denom == zero) then
            knew = 0
            return
        end if
        sqrtdn = sqrt (abs(denom))
!
!     Complete the updating of ZMAT when there is only one nonzero element
!       in the KNEW-th row of the new matrix ZMAT. IFLAG is set to one when
!       the value of IDZ is going to be reduced.
!
        iflag = 0
        if (jl == 1) then
            tempa = tau / sqrtdn
            tempb = zmat (knew, 1) / sqrtdn
            do i = 1, npt
                zmat (i, 1) = tempa * zmat (i, 1) - tempb * vlag (i)
            end do
            if (denom < zero) then
                if (idz == 1) then
                    idz = 2
                else
                    iflag = 1
                end if
            end if
        else
!
!     Complete the updating of ZMAT in the alternative case.
!
            ja = 1
            if (beta >= zero) ja = jl
            jb = jl + 1 - ja
            temp = zmat (knew, jb) / denom
            tempa = temp * beta
            tempb = temp * tau
            temp = zmat (knew, ja)
            scala = one / sqrt (abs(beta)*temp*temp+tausq)
            scalb = scala * sqrtdn
            do i = 1, npt
                zmat (i, ja) = scala * (tau*zmat(i, ja)-temp*vlag(i))
                zmat (i, jb) = scalb * (zmat(i, jb)-tempa*w(i)-tempb*vlag(i))
            end do
            if (denom <= zero) then
                if (beta < zero) then
                    idz = idz + 1
                else
                    iflag = 1
                end if
            end if
        end if
!
!     Reduce IDZ when the diagonal part of the ZMAT times Diag(DZ) times
!       ZMAT^T factorization gains another positive element. Then exchange
!       the first and IDZ-th columns of ZMAT.
!
        if (iflag == 1) then
            idz = idz - 1
            do i = 1, npt
                temp = zmat (i, 1)
                zmat (i, 1) = zmat (i, idz)
                zmat (i, idz) = temp
            end do
        end if
!
!     Finally, update the matrix BMAT.
!
        do j = 1, n
            jp = npt + j
            w (jp) = bmat (knew, j)
            tempa = (alpha*vlag(jp)-tau*w(jp)) / denom
            tempb = (-beta*w(jp)-tau*vlag(jp)) / denom
            do i = 1, jp
                bmat (i, j) = bmat (i, j) + tempa * vlag (i) + tempb * w (i)
                if (i > npt) bmat (jp, i-npt) = bmat (i, j)
            end do
        end do

    end subroutine update

!*****************************************************************************************
!>
!  Test problem for [[lincoa]].
!
!  Calculate the tetrahedron of least volume that encloses the points
!  `(XP(J),YP(J),ZP(J)), J=1,2,...,NP`. Our method requires the origin
!  to be strictly inside the convex hull of these points. There are
!  twelve variables that define the four faces of each tetrahedron
!  that is considered. Each face has the form `ALPHA*X+BETA*Y+GAMMA*Z=1`,
!  the variables `X(3K-2)`, `X(3K-1)` and `X(3K)` being the values of `ALPHA`,
!  `BETA` and `GAMMA` for the K-th face, K=1,2,3,4. Let the set T contain
!  all points in three dimensions that can be reached from the origin
!  without crossing a face. Because the volume of T may be infinite,
!  the objective function is the smaller of FMAX and the volume of T,
!  where FMAX is set to an upper bound on the final volume initially.
!  There are 4*NP linear constraints on the variables, namely that each
!  of the given points `(XP(J),YP(J),ZP(J))` shall be in T. Let `XS = min
!  XP(J)`, `YS = min YP(J)`, `ZS = min ZP(J)` and `SS = max XP(J)+YP(J)+ZP(J)`,
!  where J runs from 1 to NP. The initial values of the variables are
!  `X(1)=1/XS`, `X(5)=1/YS`, `X(9)=1/ZS`, `X(2)=X(3)=X(4)=X(6)=X(7)=X(8)=0`
!  and `X(10)=X(11)=X(12)=1/SS`, which satisfy the linear constraints,
!  and which provide the bound `FMAX=(SS-XS-YS-ZS)**3/6`. Other details
!  of the test calculation are given below, including the choice of
!  the data points `(XP(J),YP(J),ZP(J)), J=1,2,...,NP`. The smaller final
!  value of the objective function in the case NPT=35 shows that the
!  problem has local minima.

    subroutine lincoa_test ()

        implicit none

        real(wp) :: xp (50), yp (50), zp (50), a (12, 200), b (200), x (12)
        integer :: ia,n,np,j,iw,iprint,jcase,k,i,maxfun,npt,m
        real(wp) :: sumx,sumy,sumz,theta,fmax,rhobeg,rhoend,ss,xs,ys,zs
    !
    !     Set some constants.
    !
        real(wp),parameter :: one  = 1.0_wp
        real(wp),parameter :: two  = 2.0_wp
        real(wp),parameter :: zero = 0.0_wp
        real(wp),parameter :: pi   = 4.0_wp * atan(one)

        ia = 12
        n = 12
    !
    !     Set the data points.
    !
        np = 50
        sumx = zero
        sumy = zero
        sumz = zero
        do j = 1, np
            theta = real (j-1, wp) * pi / real (np-1, wp)
            xp (j) = cos (theta) * cos (two*theta)
            sumx = sumx + xp (j)
            yp (j) = sin (theta) * cos (two*theta)
            sumy = sumy + yp (j)
            zp (j) = sin (two*theta)
            sumz = sumz + zp (j)
        end do
        sumx = sumx / real (np, wp)
        sumy = sumy / real (np, wp)
        sumz = sumz / real (np, wp)
        do j = 1, np
            xp (j) = xp (j) - sumx
            yp (j) = yp (j) - sumy
            zp (j) = zp (j) - sumz
        end do
    !
    !     Set the linear constraints.
    !
        m = 4 * np
        do k = 1, m
            b (k) = one
            do i = 1, n
                a (i, k) = zero
            end do
        end do
        do j = 1, np
            do i = 1, 4
                k = 4 * j + i - 4
                iw = 3 * i
                a (iw-2, k) = xp (j)
                a (iw-1, k) = yp (j)
                a (iw, k) = zp (j)
            end do
        end do
    !
    !     Set the initial vector of variables. The JCASE=1,6 loop gives six
    !       different choices of NPT when LINCOA is called.
    !
        xs = zero
        ys = zero
        zs = zero
        ss = zero
        do j = 1, np
            xs = min (xs, xp(j))
            ys = min (ys, yp(j))
            zs = min (zs, zp(j))
            ss = max (ss, xp(j)+yp(j)+zp(j))
        end do
        fmax = (ss-xs-ys-zs) ** 3 / 6.0_wp
        do jcase = 1, 6
            do i = 2, 8
                x (i) = zero
            end do
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
            print 70, npt, rhoend
70          format (/ / 4 x, 'Output from LINCOA with  NPT =', i4, '  and  RHOEND =', 1 &
           & pd12.4)
            call lincoa(n,npt,m,a,ia,b,x,rhobeg,rhoend,iprint,maxfun,calfun)
        end do

    contains

        subroutine calfun (n, x, f)

            implicit none

            integer :: n
            real(wp) :: x (*)
            real(wp) :: f

            real(wp) :: v12,v13,v14,v23,v24,v34,del1,del2,del3,del4,temp

            f = fmax
            v12 = x (1) * x (5) - x (4) * x (2)
            v13 = x (1) * x (8) - x (7) * x (2)
            v14 = x (1) * x (11) - x (10) * x (2)
            v23 = x (4) * x (8) - x (7) * x (5)
            v24 = x (4) * x (11) - x (10) * x (5)
            v34 = x (7) * x (11) - x (10) * x (8)
            del1 = v23 * x (12) - v24 * x (9) + v34 * x (6)
            if (del1 <= zero) return
            del2 = - v34 * x (3) - v13 * x (12) + v14 * x (9)
            if (del2 <= zero) return
            del3 = - v14 * x (6) + v24 * x (3) + v12 * x (12)
            if (del3 <= zero) return
            del4 = - v12 * x (9) + v13 * x (6) - v23 * x (3)
            if (del4 <= zero) return
            temp = (del1+del2+del3+del4) ** 3 / (del1*del2*del3*del4)
            f = min (temp/6.0_wp, fmax)

        end subroutine calfun

    end subroutine lincoa_test
!*****************************************************************************************

end module lincoa_module
