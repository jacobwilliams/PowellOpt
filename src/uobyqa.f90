!*****************************************************************************************
!>
!  UOBYQA: **U**nconstrained **O**ptimization **BY** **Q**uadratic **A**pproximation
!
!  The purpose of UOBYQA is to seek the least value of a function F of several variables,
!  when derivatives are not available.
!  It uses a trust region method that forms quadratic models by interpolation.
!
!# References
!
!  * "UOBYQA: unconstrained optimization by
!    quadratic approximation" by M.J.D. Powell, Report DAMTP 2000/NA14,
!    University of Cambridge.
!  * "[UOBYQA: unconstrained optimization by quadratic
!    approximation](http://link.springer.com/article/10.1007%2Fs101070100290)" by
!    M.J.D. Powell, Mathematical Programming Series B, Volume
!    92, pages 555-582 (2002).
!
!# History
!  * M.J.D. Powell : It is hoped that the software will
!    be helpful to much future research and to many applications.
!    There are no restrictions on or charges for its use.
!  * Jacob Williams, July 2015 : refactoring of the code into modern Fortran.

    module uobyqa_module

    use kind_module, only: wp

    private

    abstract interface
    subroutine func (n, x, f)  !! calfun interface
        import :: wp
        implicit none
        integer :: n
        real (wp) :: x (*)
        real (wp) :: f
    end subroutine func
    end interface

    public :: uobyqa
    public :: uobyqa_test

    contains
!*****************************************************************************************

!*****************************************************************************************
!>
!  This subroutine seeks the least value of a function of many variables,
!  by a trust region method that forms quadratic models by interpolation.

    subroutine uobyqa (n, x, rhobeg, rhoend, iprint, maxfun, calfun)

        implicit none

        integer,intent(in)    :: n               !! the number of variables and must be at least two
        real(wp),intent(inout),dimension(*) :: x !! Initial values of the variables must be set in X(1),X(2),...,X(N). They
                                                 !! will be changed to the values that give the least calculated F.
        real(wp),intent(in)   :: rhobeg          !! RHOBEG and RHOEND must be set to the initial and final values of a trust
                                                 !! region radius, so both must be positive with RHOEND<=RHOBEG. Typically
                                                 !! RHOBEG should be about one tenth of the greatest expected change to a
                                                 !! variable, and RHOEND should indicate the accuracy that is required in
                                                 !! the final values of the variables.
        real(wp),intent(in)   :: rhoend          !! RHOBEG and RHOEND must be set to the initial and final values of a trust
                                                 !! region radius, so both must be positive with RHOEND<=RHOBEG. Typically
                                                 !! RHOBEG should be about one tenth of the greatest expected change to a
                                                 !! variable, and RHOEND should indicate the accuracy that is required in
                                                 !! the final values of the variables.
        integer,intent(in)    :: iprint          !! The value of IPRINT should be set to 0, 1, 2 or 3, which controls the
                                                 !! amount of printing. Specifically, there is no output if IPRINT=0 and
                                                 !! there is output only at the return if IPRINT=1. Otherwise, each new
                                                 !! value of RHO is printed, with the best vector of variables so far and
                                                 !! the corresponding value of the objective function. Further, each new
                                                 !! value of F with its variables are output if IPRINT=3.
        integer,intent(in)   :: maxfun           !! upper bound on the number of calls of CALFUN.
        procedure (func)     :: calfun           !! It must set F to the value of the objective
                                                 !! function for the variables X(1),X(2),...,X(N).

        real(wp),dimension(:),allocatable :: w
        integer :: npt,ixb,ixo,ixn,ixp,ipq,ipl,ih,ig,id,ivl,iw

        ! Partition the working space array, so that different parts of it can be
        ! treated separately by the subroutine that performs the main calculation.

        npt = (n*n+3*n+2) / 2
        ixb = 1
        ixo = ixb + n
        ixn = ixo + n
        ixp = ixn + n
        ipq = ixp + n * npt
        ipl = ipq + npt - 1
        ih = ipl + (npt-1) * npt
        ig = ih + n * n
        id = ig + n
        ivl = ih
        iw = id + n

        ! The array W will be used for working space
        allocate(w((N**4+8*N**3+23*N**2+42*N+max(2*N**2+4,18*N))/4))

        call uobyqb (n, x, rhobeg, rhoend, iprint, maxfun, npt, w(ixb), w(ixo), w(ixn), &
                     w(ixp), w(ipq), w(ipl), w(ih), w(ig), w(id), w(ivl), w(iw), calfun)

        deallocate(w)

    end subroutine uobyqa
!*****************************************************************************************

    subroutine uobyqb (n, x, rhobeg, rhoend, iprint, maxfun, npt, xbase, xopt, xnew, xpt, &
                       pq, pl, h, g, d, vlag, w, calfun)

        implicit real (wp) (a-h, o-z)

        dimension x (*), xbase (*), xopt (*), xnew (*), xpt (npt,*), pq (*), pl (npt,*), &
                  h (n,*), g (*), d (*), vlag (*), w (*)
        procedure (func) :: calfun
!
!     The arguments N, X, RHOBEG, RHOEND, IPRINT and MAXFUN are identical to
!       the corresponding arguments in SUBROUTINE UOBYQA.
!     NPT is set by UOBYQA to (N*N+3*N+2)/2 for the above dimension statement.
!     XBASE will contain a shift of origin that reduces the contributions from
!       rounding errors to values of the model and Lagrange functions.
!     XOPT will be set to the displacement from XBASE of the vector of
!       variables that provides the least calculated F so far.
!     XNEW will be set to the displacement from XBASE of the vector of
!       variables for the current calculation of F.
!     XPT will contain the interpolation point coordinates relative to XBASE.
!     PQ will contain the parameters of the quadratic model.
!     PL will contain the parameters of the Lagrange functions.
!     H will provide the second derivatives that TRSTEP and LAGMAX require.
!     G will provide the first derivatives that TRSTEP and LAGMAX require.
!     D is reserved for trial steps from XOPT, except that it will contain
!       diagonal second derivatives during the initialization procedure.
!     VLAG will contain the values of the Lagrange functions at a new point X.
!     The array W will be used for working space. Its length must be at least
!     max [ 6*N, ( N**2 + 3*N + 2 ) / 2 ].
!
!     Set some constants.
!
        one = 1.0_wp
        two = 2.0_wp
        zero = 0.0_wp
        half = 0.5_wp
        tol = 0.01_wp
        nnp = n + n + 1
        nptm = npt - 1
        nftest = max (maxfun, 1)
!
!     Initialization. NF is the number of function calculations so far.
!
        rho = rhobeg
        rhosq = rho * rho
        nf = 0
        do i = 1, n
            xbase (i) = x (i)
            do k = 1, npt
                xpt (k, i) = zero
            end do
        end do
        do k = 1, npt
            do j = 1, nptm
                pl (k, j) = zero
            end do
        end do
!
!     The branch to label 120 obtains a new value of the objective function
!     and then there is a branch back to label 50, because the new function
!     value is needed to form the initial quadratic model. The least function
!     value so far and its index are noted below.
!
30      do i = 1, n
            x (i) = xbase (i) + xpt (nf+1, i)
        end do
        go to 120
50      if (nf == 1) then
            fopt = f
            kopt = nf
            fbase = f
            j = 0
            jswitch = - 1
            ih = n
        else
            if (f < fopt) then
                fopt = f
                kopt = nf
            end if
        end if
!
!     Form the gradient and diagonal second derivatives of the initial
!     quadratic model and Lagrange functions.
!
        if (nf <= nnp) then
            jswitch = - jswitch
            if (jswitch > 0) then
                if (j >= 1) then
                    ih = ih + j
                    if (w(j) < zero) then
                        d (j) = (fsave+f-two*fbase) / rhosq
                        pq (j) = (fsave-f) / (two*rho)
                        pl (1, ih) = - two / rhosq
                        pl (nf-1, j) = half / rho
                        pl (nf-1, ih) = one / rhosq
                    else
                        pq (j) = (4.0_wp*fsave-3.0_wp*fbase-f) / (two*rho)
                        d (j) = (fbase+f-two*fsave) / rhosq
                        pl (1, j) = - 1.5_wp / rho
                        pl (1, ih) = one / rhosq
                        pl (nf-1, j) = two / rho
                        pl (nf-1, ih) = - two / rhosq
                    end if
                    pq (ih) = d (j)
                    pl (nf, j) = - half / rho
                    pl (nf, ih) = one / rhosq
                end if
!
!     Pick the shift from XBASE to the next initial interpolation point
!     that provides diagonal second derivatives.
!
                if (j < n) then
                    j = j + 1
                    xpt (nf+1, j) = rho
                end if
            else
                fsave = f
                if (f < fbase) then
                    w (j) = rho
                    xpt (nf+1, j) = two * rho
                else
                    w (j) = - rho
                    xpt (nf+1, j) = - rho
                end if
            end if
            if (nf < nnp) go to 30
!
!     Form the off-diagonal second derivatives of the initial quadratic model.
!
            ih = n
            ip = 1
            iq = 2
        end if
        ih = ih + 1
        if (nf > nnp) then
            temp = one / (w(ip)*w(iq))
            tempa = f - fbase - w (ip) * pq (ip) - w (iq) * pq (iq)
            pq (ih) = (tempa-half*rhosq*(d(ip)+d(iq))) * temp
            pl (1, ih) = temp
            iw = ip + ip
            if (w(ip) < zero) iw = iw + 1
            pl (iw, ih) = - temp
            iw = iq + iq
            if (w(iq) < zero) iw = iw + 1
            pl (iw, ih) = - temp
            pl (nf, ih) = temp
!
!     Pick the shift from XBASE to the next initial interpolation point
!     that provides off-diagonal second derivatives.
!
            ip = ip + 1
        end if
        if (ip == iq) then
            ih = ih + 1
            ip = 1
            iq = iq + 1
        end if
        if (nf < npt) then
            xpt (nf+1, ip) = w (ip)
            xpt (nf+1, iq) = w (iq)
            go to 30
        end if
!
!     Set parameters to begin the iterations for the current RHO.
!
        sixthm = zero
        delta = rho
60      tworsq = (two*rho) ** 2
        rhosq = rho * rho
!
!     Form the gradient of the quadratic model at the trust region centre.
!
70      knew = 0
        ih = n
        do j = 1, n
            xopt (j) = xpt (kopt, j)
            g (j) = pq (j)
            do i = 1, j
                ih = ih + 1
                g (i) = g (i) + pq (ih) * xopt (j)
                if (i < j) g (j) = g (j) + pq (ih) * xopt (i)
                h (i, j) = pq (ih)
            end do
        end do
!
!     Generate the next trust region step and test its length. Set KNEW
!     to -1 if the purpose of the next F will be to improve conditioning,
!     and also calculate a lower bound on the Hessian term of the model Q.
!
        call trstep (n, g, h, delta, tol, d, w(1), w(n+1), w(2*n+1), w(3*n+1), w(4*n+1), &
                     w(5*n+1), evalue)
        temp = zero
        do i = 1, n
            temp = temp + d (i) ** 2
        end do
        dnorm = min (delta, sqrt(temp))
        errtol = - one
        if (dnorm < half*rho) then
            knew = - 1
            errtol = half * evalue * rho * rho
            if (nf <= npt+9) errtol = zero
            go to 290
        end if
!
!     Calculate the next value of the objective function.
!
100     do i = 1, n
            xnew (i) = xopt (i) + d (i)
            x (i) = xbase (i) + xnew (i)
        end do
120     if (nf >= nftest) then
            if (iprint > 0) print 130
130         format (/ 4 x, 'Return from UOBYQA because CALFUN has been',&
            ' called MAXFUN times')
            go to 420
        end if
        nf = nf + 1
        call calfun (n, x, f)
        if (iprint == 3) then
            print 140, nf, f, (x(i), i=1, n)
140         format (/ 4 x, 'Function number', i6, '    F =', 1 pd18.10,&
            '    The corresponding X is:' / (2 x, 5d15.6))
        end if
        if (nf <= npt) go to 50
        if (knew ==-1) go to 420
!
!     Use the quadratic model to predict the change in F due to the step D,
!     and find the values of the Lagrange functions at the new point.
!
        vquad = zero
        ih = n
        do j = 1, n
            w (j) = d (j)
            vquad = vquad + w (j) * pq (j)
            do i = 1, j
                ih = ih + 1
                w (ih) = d (i) * xnew (j) + d (j) * xopt (i)
                if (i == j) w (ih) = half * w (ih)
                vquad = vquad + w (ih) * pq (ih)
            end do
        end do
        do k = 1, npt
            temp = zero
            do j = 1, nptm
                temp = temp + w (j) * pl (k, j)
            end do
            vlag (k) = temp
        end do
        vlag (kopt) = vlag (kopt) + one
!
!     Update SIXTHM, which is a lower bound on one sixth of the greatest
!     third derivative of F.
!
        diff = f - fopt - vquad
        sum = zero
        do k = 1, npt
            temp = zero
            do i = 1, n
                temp = temp + (xpt(k, i)-xnew(i)) ** 2
            end do
            temp = sqrt (temp)
            sum = sum + abs (temp*temp*temp*vlag(k))
        end do
        sixthm = max (sixthm, abs(diff)/sum)
!
!     Update FOPT and XOPT if the new F is the least value of the objective
!     function so far. Then branch if D is not a trust region step.
!
        fsave = fopt
        if (f < fopt) then
            fopt = f
            do i = 1, n
                xopt (i) = xnew (i)
            end do
        end if
        ksave = knew
        if (knew > 0) go to 240
!
!     Pick the next value of DELTA after a trust region step.
!
        if (vquad >= zero) then
            if (iprint > 0) print 210
210         format (/ 4 x, 'Return from UOBYQA because a trust',&
            ' region step has failed to reduce Q')
            go to 420
        end if
        ratio = (f-fsave) / vquad
        if (ratio <= 0.1_wp) then
            delta = half * dnorm
        else if (ratio <= 0.7_wp) then
            delta = max (half*delta, dnorm)
        else
            delta = max (delta, 1.25_wp*dnorm, dnorm+rho)
        end if
        if (delta <= 1.5_wp*rho) delta = rho
!
!     Set KNEW to the index of the next interpolation point to be deleted.
!
        ktemp = 0
        detrat = zero
        if (f >= fsave) then
            ktemp = kopt
            detrat = one
        end if
        do k = 1, npt
            sum = zero
            do i = 1, n
                sum = sum + (xpt(k, i)-xopt(i)) ** 2
            end do
            temp = abs (vlag(k))
            if (sum > rhosq) temp = temp * (sum/rhosq) ** 1.5_wp
            if (temp > detrat .and. k /= ktemp) then
                detrat = temp
                ddknew = sum
                knew = k
            end if
        end do
        if (knew == 0) go to 290
!
!     Replace the interpolation point that has index KNEW by the point XNEW,
!     and also update the Lagrange functions and the quadratic model.
!
240     do i = 1, n
            xpt (knew, i) = xnew (i)
        end do
        temp = one / vlag (knew)
        do j = 1, nptm
            pl (knew, j) = temp * pl (knew, j)
            pq (j) = pq (j) + diff * pl (knew, j)
        end do
        do k = 1, npt
            if (k /= knew) then
                temp = vlag (k)
                do j = 1, nptm
                    pl (k, j) = pl (k, j) - temp * pl (knew, j)
                end do
            end if
        end do
!
!     Update KOPT if F is the least calculated value of the objective
!     function. Then branch for another trust region calculation. The
!     case KSAVE>0 indicates that a model step has just been taken.
!
        if (f < fsave) then
            kopt = knew
            go to 70
        end if
        if (ksave > 0) go to 70
        if (dnorm > two*rho) go to 70
        if (ddknew > tworsq) go to 70
!
!     Alternatively, find out if the interpolation points are close
!     enough to the best point so far.
!
290     do k = 1, npt
            w (k) = zero
            do i = 1, n
                w (k) = w (k) + (xpt(k, i)-xopt(i)) ** 2
            end do
        end do
310     knew = - 1
        distest = tworsq
        do k = 1, npt
            if (w(k) > distest) then
                knew = k
                distest = w (k)
            end if
        end do
!
!     If a point is sufficiently far away, then set the gradient and Hessian
!     of its Lagrange function at the centre of the trust region, and find
!     half the sum of squares of components of the Hessian.
!
        if (knew > 0) then
            ih = n
            sumh = zero
            do j = 1, n
                g (j) = pl (knew, j)
                do i = 1, j
                    ih = ih + 1
                    temp = pl (knew, ih)
                    g (j) = g (j) + temp * xopt (i)
                    if (i < j) then
                        g (i) = g (i) + temp * xopt (j)
                        sumh = sumh + temp * temp
                    end if
                    h (i, j) = temp
                end do
                sumh = sumh + half * temp * temp
            end do
!
!     If ERRTOL is positive, test whether to replace the interpolation point
!     with index KNEW, using a bound on the maximum modulus of its Lagrange
!     function in the trust region.
!
            if (errtol > zero) then
                w (knew) = zero
                sumg = zero
                do i = 1, n
                    sumg = sumg + g (i) ** 2
                end do
                estim = rho * (sqrt(sumg)+rho*sqrt(half*sumh))
                wmult = sixthm * distest ** 1.5_wp
                if (wmult*estim <= errtol) go to 310
            end if
!
!     If the KNEW-th point may be replaced, then pick a D that gives a large
!     value of the modulus of its Lagrange function within the trust region.
!     Here the vector XNEW is used as temporary working space.
!
            call lagmax (n, g, h, rho, d, xnew, vmax)
            if (errtol > zero) then
                if (wmult*vmax <= errtol) go to 310
            end if
            go to 100
        end if
        if (dnorm > rho) go to 70
!
!     Prepare to reduce RHO by shifting XBASE to the best point so far,
!     and make the corresponding changes to the gradients of the Lagrange
!     functions and the quadratic model.
!
        if (rho > rhoend) then
            ih = n
            do j = 1, n
                xbase (j) = xbase (j) + xopt (j)
                do k = 1, npt
                    xpt (k, j) = xpt (k, j) - xopt (j)
                end do
                do i = 1, j
                    ih = ih + 1
                    pq (i) = pq (i) + pq (ih) * xopt (j)
                    if (i < j) then
                        pq (j) = pq (j) + pq (ih) * xopt (i)
                        do k = 1, npt
                            pl (k, j) = pl (k, j) + pl (k, ih) * xopt (i)
                        end do
                    end if
                    do k = 1, npt
                        pl (k, i) = pl (k, i) + pl (k, ih) * xopt (j)
                    end do
                end do
            end do
!
!     Pick the next values of RHO and DELTA.
!
            delta = half * rho
            ratio = rho / rhoend
            if (ratio <= 16.0_wp) then
                rho = rhoend
            else if (ratio <= 250.0_wp) then
                rho = sqrt (ratio) * rhoend
            else
                rho = 0.1_wp * rho
            end if
            delta = max (delta, rho)
            if (iprint >= 2) then
                if (iprint >= 3) print 390
390             format (5 x)
                print 400, rho, nf
400             format (/ 4 x, 'New RHO =', 1 pd11.4, 5 x, 'Number of',&
                        ' function values =', i6)
                print 410, fopt, (xbase(i), i=1, n)
410             format (4 x, 'Least value of F =', 1 pd23.15, 9 x,&
                'The corresponding X is:'/(2 x, 5d15.6))
            end if
            go to 60
        end if
!
!     Return from the calculation, after another Newton-Raphson step, if
!     it is too short to have been tried before.
!
        if (errtol >= zero) go to 100
420     if (fopt <= f) then
            do i = 1, n
                x (i) = xbase (i) + xopt (i)
            end do
            f = fopt
        end if
        if (iprint >= 1) then
            print 440, nf
440         format (/ 4 x, 'At the return from UOBYQA', 5 x,&
            'Number of function values =', i6)
            print 410, f, (x(i), i=1, n)
        end if

    end subroutine uobyqb

    subroutine lagmax (n, g, h, rho, d, v, vmax)

        implicit real (wp) (a-h, o-z)

        dimension g (*), h (n,*), d (*), v (*)
!
!     N is the number of variables of a quadratic objective function, Q say.
!     G is the gradient of Q at the origin.
!     H is the symmetric Hessian matrix of Q. Only the upper triangular and
!       diagonal parts need be set.
!     RHO is the trust region radius, and has to be positive.
!     D will be set to the calculated vector of variables.
!     The array V will be used for working space.
!     VMAX will be set to |Q(0)-Q(D)|.
!
!     Calculating the D that maximizes |Q(0)-Q(D)| subject to ||D|| .LEQ. RHO
!     requires of order N**3 operations, but sometimes it is adequate if
!     |Q(0)-Q(D)| is within about 0.9 of its greatest possible value. This
!     subroutine provides such a solution in only of order N**2 operations,
!     where the claim of accuracy has been tested by numerical experiments.
!
!     Preliminary calculations.
!
        half = 0.5_wp
        halfrt = sqrt (half)
        one = 1.0_wp
        zero = 0.0_wp
!
!     Pick V such that ||HV|| / ||V|| is large.
!
        hmax = zero
        do i = 1, n
            sum = zero
            do j = 1, n
                h (j, i) = h (i, j)
                sum = sum + h (i, j) ** 2
            end do
            if (sum > hmax) then
                hmax = sum
                k = i
            end if
        end do
        do j = 1, n
            v (j) = h (k, j)
        end do
!
!     Set D to a vector in the subspace spanned by V and HV that maximizes
!     |(D,HD)|/(D,D), except that we set D=HV if V and HV are nearly parallel.
!     The vector that has the name D at label 60 used to be the vector W.
!
        vsq = zero
        vhv = zero
        dsq = zero
        do i = 1, n
            vsq = vsq + v (i) ** 2
            d (i) = zero
            do j = 1, n
                d (i) = d (i) + h (i, j) * v (j)
            end do
            vhv = vhv + v (i) * d (i)
            dsq = dsq + d (i) ** 2
        end do
        if (vhv*vhv <= 0.9999_wp*dsq*vsq) then
            temp = vhv / vsq
            wsq = zero
            do i = 1, n
                d (i) = d (i) - temp * v (i)
                wsq = wsq + d (i) ** 2
            end do
            whw = zero
            ratio = sqrt (wsq/vsq)
            do i = 1, n
                temp = zero
                do j = 1, n
                    temp = temp + h (i, j) * d (j)
                end do
                whw = whw + temp * d (i)
                v (i) = ratio * v (i)
            end do
            vhv = ratio * ratio * vhv
            vhw = ratio * wsq
            temp = half * (whw-vhv)
            temp = temp + sign (sqrt(temp**2+vhw**2), whw+vhv)
            do i = 1, n
                d (i) = vhw * v (i) + temp * d (i)
            end do
        end if
!
!     We now turn our attention to the subspace spanned by G and D. A multiple
!     of the current D is returned if that choice seems to be adequate.
!
        gg = zero
        gd = zero
        dd = zero
        dhd = zero
        do i = 1, n
            gg = gg + g (i) ** 2
            gd = gd + g (i) * d (i)
            dd = dd + d (i) ** 2
            sum = zero
            do j = 1, n
                sum = sum + h (i, j) * d (j)
            end do
            dhd = dhd + sum * d (i)
        end do
        temp = gd / gg
        vv = zero
        scale = sign (rho/sqrt(dd), gd*dhd)
        do i = 1, n
            v (i) = d (i) - temp * g (i)
            vv = vv + v (i) ** 2
            d (i) = scale * d (i)
        end do
        gnorm = sqrt (gg)
        if (gnorm*dd <= 0.5e-2_wp*rho*abs(dhd) .or. vv/dd <= 1.0e-4_wp) then
            vmax = abs (scale*(gd+half*scale*dhd))
            return
        end if
!
!     G and V are now orthogonal in the subspace spanned by G and D. Hence
!     we generate an orthonormal basis of this subspace such that (D,HV) is
!     negligible or zero, where D and V will be the basis vectors.
!
        ghg = zero
        vhg = zero
        vhv = zero
        do i = 1, n
            sum = zero
            sumv = zero
            do j = 1, n
                sum = sum + h (i, j) * g (j)
                sumv = sumv + h (i, j) * v (j)
            end do
            ghg = ghg + sum * g (i)
            vhg = vhg + sumv * g (i)
            vhv = vhv + sumv * v (i)
        end do
        vnorm = sqrt (vv)
        ghg = ghg / gg
        vhg = vhg / (vnorm*gnorm)
        vhv = vhv / vv
        if (abs(vhg) <= 0.01_wp*max(abs(ghg), abs(vhv))) then
            vmu = ghg - vhv
            wcos = one
            wsin = zero
        else
            temp = half * (ghg-vhv)
            vmu = temp + sign (sqrt(temp**2+vhg**2), temp)
            temp = sqrt (vmu**2+vhg**2)
            wcos = vmu / temp
            wsin = vhg / temp
        end if
        tempa = wcos / gnorm
        tempb = wsin / vnorm
        tempc = wcos / vnorm
        tempd = wsin / gnorm
        do i = 1, n
            d (i) = tempa * g (i) + tempb * v (i)
            v (i) = tempc * v (i) - tempd * g (i)
        end do
!
!     The final D is a multiple of the current D, V, D+V or D-V. We make the
!     choice from these possibilities that is optimal.
!
        dlin = wcos * gnorm / rho
        vlin = - wsin * gnorm / rho
        tempa = abs (dlin) + half * abs (vmu+vhv)
        tempb = abs (vlin) + half * abs (ghg-vmu)
        tempc = halfrt * (abs(dlin)+abs(vlin)) + 0.25_wp * abs (ghg+vhv)
        if (tempa >= tempb .and. tempa >= tempc) then
            tempd = sign (rho, dlin*(vmu+vhv))
            tempv = zero
        else if (tempb >= tempc) then
            tempd = zero
            tempv = sign (rho, vlin*(ghg-vmu))
        else
            tempd = sign (halfrt*rho, dlin*(ghg+vhv))
            tempv = sign (halfrt*rho, vlin*(ghg+vhv))
        end if
        do i = 1, n
            d (i) = tempd * d (i) + tempv * v (i)
        end do
        vmax = rho * rho * max (tempa, tempb, tempc)

    end subroutine lagmax

    subroutine trstep (n, g, h, delta, tol, d, gg, td, tn, w, piv, z, evalue)

        implicit real (wp) (a-h, o-z)

        dimension g (*), h (n,*), d (*), gg (*), td (*), tn (*), w (*), piv (*), z (*)
!
!     N is the number of variables of a quadratic objective function, Q say.
!     G is the gradient of Q at the origin.
!     H is the Hessian matrix of Q. Only the upper triangular and diagonal
!       parts need be set. The lower triangular part is used to store the
!       elements of a Householder similarity transformation.
!     DELTA is the trust region radius, and has to be positive.
!     TOL is the value of a tolerance from the open interval (0,1).
!     D will be set to the calculated vector of variables.
!     The arrays GG, TD, TN, W, PIV and Z will be used for working space.
!     EVALUE will be set to the least eigenvalue of H if and only if D is a
!     Newton-Raphson step. Then EVALUE will be positive, but otherwise it
!     will be set to zero.
!
!     Let MAXRED be the maximum of Q(0)-Q(D) subject to ||D|| .LEQ. DELTA,
!     and let ACTRED be the value of Q(0)-Q(D) that is actually calculated.
!     We take the view that any D is acceptable if it has the properties
!
!             ||D|| .LEQ. DELTA  and  ACTRED .LEQ. (1-TOL)*MAXRED.
!
!     The calculation of D is done by the method of Section 2 of the paper
!     by MJDP in the 1997 Dundee Numerical Analysis Conference Proceedings,
!     after transforming H to tridiagonal form.
!
!     Initialization.
!
        one = 1.0_wp
        two = 2.0_wp
        zero = 0.0_wp
        delsq = delta * delta
        evalue = zero
        nm = n - 1
        do i = 1, n
            d (i) = zero
            td (i) = h (i, i)
            do j = 1, i
                h (i, j) = h (j, i)
            end do
        end do
!
!     Apply Householder transformations to obtain a tridiagonal matrix that
!     is similar to H, and put the elements of the Householder vectors in
!     the lower triangular part of H. Further, TD and TN will contain the
!     diagonal and other nonzero elements of the tridiagonal matrix.
!
        do k = 1, nm
            kp = k + 1
            sum = zero
            if (kp < n) then
                kpp = kp + 1
                do i = kpp, n
                    sum = sum + h (i, k) ** 2
                end do
            end if
            if (sum == zero) then
                tn (k) = h (kp, k)
                h (kp, k) = zero
            else
                temp = h (kp, k)
                tn (k) = sign (sqrt(sum+temp*temp), temp)
                h (kp, k) = - sum / (temp+tn(k))
                temp = sqrt (two/(sum+h(kp, k)**2))
                do i = kp, n
                    w (i) = temp * h (i, k)
                    h (i, k) = w (i)
                    z (i) = td (i) * w (i)
                end do
                wz = zero
                do j = kp, nm
                    jp = j + 1
                    do i = jp, n
                        z (i) = z (i) + h (i, j) * w (j)
                        z (j) = z (j) + h (i, j) * w (i)
                    end do
                    wz = wz + w (j) * z (j)
                end do
                wz = wz + w (n) * z (n)
                do j = kp, n
                    td (j) = td (j) + w (j) * (wz*w(j)-two*z(j))
                    if (j < n) then
                        jp = j + 1
                        do i = jp, n
                            h (i, j) = h (i, j) - w (i) * z (j) - w (j) * (z(i)-wz*w(i))
                        end do
                    end if
                end do
            end if
        end do
!
!     Form GG by applying the similarity transformation to G.
!
        gsq = zero
        do i = 1, n
            gg (i) = g (i)
            gsq = gsq + g (i) ** 2
        end do
        gnorm = sqrt (gsq)
        do k = 1, nm
            kp = k + 1
            sum = zero
            do i = kp, n
                sum = sum + gg (i) * h (i, k)
            end do
            do i = kp, n
                gg (i) = gg (i) - sum * h (i, k)
            end do
        end do
!
!     Begin the trust region calculation with a tridiagonal matrix by
!     calculating the norm of H. Then treat the case when H is zero.
!
        hnorm = abs (td(1)) + abs (tn(1))
        tdmin = td (1)
        tn (n) = zero
        do i = 2, n
            temp = abs (tn(i-1)) + abs (td(i)) + abs (tn(i))
            hnorm = max (hnorm, temp)
            tdmin = min (tdmin, td(i))
        end do
        if (hnorm == zero) then
            if (gnorm == zero) return
            scale = delta / gnorm
            do i = 1, n
                d (i) = - scale * gg (i)
            end do
            go to 370
        end if
!
!     Set the initial values of PAR and its bounds.
!
        parl = max (zero,-tdmin, gnorm/delta-hnorm)
        parlest = parl
        par = parl
        paru = zero
        paruest = zero
        posdef = zero
        iterc = 0
!
!     Calculate the pivots of the Cholesky factorization of (H+PAR*I).
!
140     iterc = iterc + 1
        ksav = 0
        piv (1) = td (1) + par
        k = 1
150     if (piv(k) > zero) then
            piv (k+1) = td (k+1) + par - tn (k) ** 2 / piv (k)
        else
            if (piv(k) < zero .or. tn(k) /= zero) go to 160
            ksav = k
            piv (k+1) = td (k+1) + par
        end if
        k = k + 1
        if (k < n) go to 150
        if (piv(k) < zero) go to 160
        if (piv(k) == zero) ksav = k
!
!     Branch if all the pivots are positive, allowing for the case when
!     G is zero.
!
        if (ksav == 0 .and. gsq > zero) go to 230
        if (gsq == zero) then
            if (par == zero) go to 370
            paru = par
            paruest = par
            if (ksav == 0) go to 190
        end if
        k = ksav
!
!     Set D to a direction of nonpositive curvature of the given tridiagonal
!     matrix, and thus revise PARLEST.
!
160     d (k) = one
        if (abs(tn(k)) <= abs(piv(k))) then
            dsq = one
            dhd = piv (k)
        else
            temp = td (k+1) + par
            if (temp <= abs(piv(k))) then
                d (k+1) = sign (one,-tn(k))
                dhd = piv (k) + temp - two * abs (tn(k))
            else
                d (k+1) = - tn (k) / temp
                dhd = piv (k) + tn (k) * d (k+1)
            end if
            dsq = one + d (k+1) ** 2
        end if
170     if (k > 1) then
            k = k - 1
            if (tn(k) /= zero) then
                d (k) = - tn (k) * d (k+1) / piv (k)
                dsq = dsq + d (k) ** 2
                go to 170
            end if
            do i = 1, k
                d (i) = zero
            end do
        end if
        parl = par
        parlest = par - dhd / dsq
!
!     Terminate with D set to a multiple of the current D if the following
!     test suggests that it suitable to do so.
!
190     temp = paruest
        if (gsq == zero) temp = temp * (one-tol)
        if (paruest > zero .and. parlest >= temp) then
            dtg = zero
            do i = 1, n
                dtg = dtg + d (i) * gg (i)
            end do
            scale = - sign (delta/sqrt(dsq), dtg)
            do i = 1, n
                d (i) = scale * d (i)
            end do
            go to 370
        end if
!
!     Pick the value of PAR for the next iteration.
!
220     if (paru == zero) then
            par = two * parlest + gnorm / delta
        else
            par = 0.5_wp * (parl+paru)
            par = max (par, parlest)
        end if
        if (paruest > zero) par = min (par, paruest)
        go to 140
!
!     Calculate D for the current PAR in the positive definite case.
!
230     w (1) = - gg (1) / piv (1)
        do i = 2, n
            w (i) = (-gg(i)-tn(i-1)*w(i-1)) / piv (i)
        end do
        d (n) = w (n)
        do i = nm, 1, - 1
            d (i) = w (i) - tn (i) * d (i+1) / piv (i)
        end do
!
!     Branch if a Newton-Raphson step is acceptable.
!
        dsq = zero
        wsq = zero
        do i = 1, n
            dsq = dsq + d (i) ** 2
            wsq = wsq + piv (i) * w (i) ** 2
        end do
        if (par == zero .and. dsq <= delsq) go to 320
!
!     Make the usual test for acceptability of a full trust region step.
!
        dnorm = sqrt (dsq)
        phi = one / dnorm - one / delta
        temp = tol * (one+par*dsq/wsq) - dsq * phi * phi
        if (temp >= zero) then
            scale = delta / dnorm
            do i = 1, n
                d (i) = scale * d (i)
            end do
            go to 370
        end if
        if (iterc >= 2 .and. par <= parl) go to 370
        if (paru > zero .and. par >= paru) go to 370
!
!     Complete the iteration when PHI is negative.
!
        if (phi < zero) then
            parlest = par
            if (posdef == one) then
                if (phi <= phil) go to 370
                slope = (phi-phil) / (par-parl)
                parlest = par - phi / slope
            end if
            slope = one / gnorm
            if (paru > zero) slope = (phiu-phi) / (paru-par)
            temp = par - phi / slope
            if (paruest > zero) temp = min (temp, paruest)
            paruest = temp
            posdef = one
            parl = par
            phil = phi
            go to 220
        end if
!
!     If required, calculate Z for the alternative test for convergence.
!
        if (posdef == zero) then
            w (1) = one / piv (1)
            do i = 2, n
                temp = - tn (i-1) * w (i-1)
                w (i) = (sign(one, temp)+temp) / piv (i)
            end do
            z (n) = w (n)
            do i = nm, 1, - 1
                z (i) = w (i) - tn (i) * z (i+1) / piv (i)
            end do
            wwsq = zero
            zsq = zero
            dtz = zero
            do i = 1, n
                wwsq = wwsq + piv (i) * w (i) ** 2
                zsq = zsq + z (i) ** 2
                dtz = dtz + d (i) * z (i)
            end do
!
!     Apply the alternative test for convergence.
!
            tempa = abs (delsq-dsq)
            tempb = sqrt (dtz*dtz+tempa*zsq)
            gam = tempa / (sign(tempb, dtz)+dtz)
            temp = tol * (wsq+par*delsq) - gam * gam * wwsq
            if (temp >= zero) then
                do i = 1, n
                    d (i) = d (i) + gam * z (i)
                end do
                go to 370
            end if
            parlest = max (parlest, par-wwsq/zsq)
        end if
!
!     Complete the iteration when PHI is positive.
!
        slope = one / gnorm
        if (paru > zero) then
            if (phi >= phiu) go to 370
            slope = (phiu-phi) / (paru-par)
        end if
        parlest = max (parlest, par-phi/slope)
        paruest = par
        if (posdef == one) then
            slope = (phi-phil) / (par-parl)
            paruest = par - phi / slope
        end if
        paru = par
        phiu = phi
        go to 220
!
!     Set EVALUE to the least eigenvalue of the second derivative matrix if
!     D is a Newton-Raphson step. SHFMAX will be an upper bound on EVALUE.
!
320     shfmin = zero
        pivot = td (1)
        shfmax = pivot
        do k = 2, n
            pivot = td (k) - tn (k-1) ** 2 / pivot
            shfmax = min (shfmax, pivot)
        end do
!
!     Find EVALUE by a bisection method, but occasionally SHFMAX may be
!     adjusted by the rule of false position.
!
        ksave = 0
340     shift = 0.5_wp * (shfmin+shfmax)
        k = 1
        temp = td (1) - shift
350     if (temp > zero) then
            piv (k) = temp
            if (k < n) then
                temp = td (k+1) - shift - tn (k) ** 2 / temp
                k = k + 1
                go to 350
            end if
            shfmin = shift
        else
            if (k < ksave) go to 360
            if (k == ksave) then
                if (pivksv == zero) go to 360
                if (piv(k)-temp < temp-pivksv) then
                    pivksv = temp
                    shfmax = shift
                else
                    pivksv = zero
                    shfmax = (shift*piv(k)-shfmin*temp) / (piv(k)-temp)
                end if
            else
                ksave = k
                pivksv = temp
                shfmax = shift
            end if
        end if
        if (shfmin <= 0.99_wp*shfmax) go to 340
360     evalue = shfmin
!
!     Apply the inverse Householder transformations to D.
!
370     nm = n - 1
        do k = nm, 1, - 1
            kp = k + 1
            sum = zero
            do i = kp, n
                sum = sum + d (i) * h (i, k)
            end do
            do i = kp, n
                d (i) = d (i) - sum * h (i, k)
            end do
        end do

    end subroutine trstep

!*****************************************************************************************
!>
!  The Chebyquad test problem (Fletcher, 1965) for N = 2,4,6,8.

    subroutine uobyqa_test ()

        implicit none

        real(wp) :: x (10)
        integer :: iprint,maxfun,i,n
        real(wp) :: rhoend,rhobeg

        iprint = 2
        maxfun = 5000
        rhoend = 1.0e-8_wp
        do n = 2, 8, 2
            do i = 1, n
                x (i) = real (i, wp) / real (n+1, wp)
            end do
            rhobeg = 0.2_wp * x (1)
            print 20, n
20          format (/ / 5 x, '******************' / 5 x, 'Results with N =', i2, / 5 x,&
                    '******************')
            call uobyqa (n, x, rhobeg, rhoend, iprint, maxfun, calfun)
        end do

    contains

        subroutine calfun (n, x, f)

            implicit none

            integer :: n
            real (wp) :: x (*)
            real (wp) :: f

            real(wp) :: y (10, 10)
            real(wp) :: sum
            integer :: j,i,iw,np

            do j = 1, n
                y (1, j) = 1.0_wp
                y (2, j) = 2.0_wp * x (j) - 1.0_wp
            end do
            do i = 2, n
                do j = 1, n
                    y (i+1, j) = 2.0_wp * y (2, j) * y (i, j) - y (i-1, j)
                end do
            end do
            f = 0.0_wp
            np = n + 1
            iw = 1
            do i = 1, np
                sum = 0.0_wp
                do j = 1, n
                    sum = sum + y (i, j)
                end do
                sum = sum / real (n, wp)
                if (iw > 0) sum = sum + 1.0_wp / real (i*i-2*i, wp)
                iw = - iw
                f = f + sum * sum
            end do

        end subroutine calfun

    end subroutine uobyqa_test
!*****************************************************************************************

!*****************************************************************************************
    end module uobyqa_module
!*****************************************************************************************