Subroutine calfun (n, x, f)
      Implicit real * 8 (a-h, o-z)
      Dimension x (*), y (10, 10)
      Do 10 j = 1, n
         y (1, j) = 1.0d0
10    y (2, j) = 2.0d0 * x (j) - 1.0d0
      Do 20 i = 2, n
         Do 20 j = 1, n
20    y (i+1, j) = 2.0d0 * y (2, j) * y (i, j) - y (i-1, j)
      f = 0.0d0
      np = n + 1
      iw = 1
      Do 40 i = 1, np
         sum = 0.0d0
         Do 30 j = 1, n
30       sum = sum + y (i, j)
         sum = sum / dfloat (n)
         If (iw > 0) sum = sum + 1.0 / dfloat (i*i-2*i)
         iw = - iw
40    f = f + sum * sum
      Return
End
Subroutine lagmax (n, g, h, rho, d, v, vmax)
      Implicit real * 8 (a-h, o-z)
      Dimension g (*), h (n,*), d (*), v (*)
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
      half = 0.5d0
      halfrt = dsqrt (half)
      one = 1.0d0
      zero = 0.0d0
!
!     Pick V such that ||HV|| / ||V|| is large.
!
      hmax = zero
      Do 20 i = 1, n
         sum = zero
         Do 10 j = 1, n
            h (j, i) = h (i, j)
10       sum = sum + h (i, j) ** 2
         If (sum > hmax) Then
            hmax = sum
            k = i
         End If
20    Continue
      Do 30 j = 1, n
30    v (j) = h (k, j)
!
!     Set D to a vector in the subspace spanned by V and HV that maximizes
!     |(D,HD)|/(D,D), except that we set D=HV if V and HV are nearly parallel.
!     The vector that has the name D at label 60 used to be the vector W.
!
      vsq = zero
      vhv = zero
      dsq = zero
      Do 50 i = 1, n
         vsq = vsq + v (i) ** 2
         d (i) = zero
         Do 40 j = 1, n
40       d (i) = d (i) + h (i, j) * v (j)
         vhv = vhv + v (i) * d (i)
50    dsq = dsq + d (i) ** 2
      If (vhv*vhv <= 0.9999d0*dsq*vsq) Then
         temp = vhv / vsq
         wsq = zero
         Do 60 i = 1, n
            d (i) = d (i) - temp * v (i)
60       wsq = wsq + d (i) ** 2
         whw = zero
         ratio = dsqrt (wsq/vsq)
         Do 80 i = 1, n
            temp = zero
            Do 70 j = 1, n
70          temp = temp + h (i, j) * d (j)
            whw = whw + temp * d (i)
80       v (i) = ratio * v (i)
         vhv = ratio * ratio * vhv
         vhw = ratio * wsq
         temp = half * (whw-vhv)
         temp = temp + dsign (dsqrt(temp**2+vhw**2), whw+vhv)
         Do 90 i = 1, n
90       d (i) = vhw * v (i) + temp * d (i)
      End If
!
!     We now turn our attention to the subspace spanned by G and D. A multiple
!     of the current D is returned if that choice seems to be adequate.
!
      gg = zero
      gd = zero
      dd = zero
      dhd = zero
      Do 110 i = 1, n
         gg = gg + g (i) ** 2
         gd = gd + g (i) * d (i)
         dd = dd + d (i) ** 2
         sum = zero
         Do 100 j = 1, n
100      sum = sum + h (i, j) * d (j)
110   dhd = dhd + sum * d (i)
      temp = gd / gg
      vv = zero
      scale = dsign (rho/dsqrt(dd), gd*dhd)
      Do 120 i = 1, n
         v (i) = d (i) - temp * g (i)
         vv = vv + v (i) ** 2
120   d (i) = scale * d (i)
      gnorm = dsqrt (gg)
      If (gnorm*dd <= 0.5d-2*rho*dabs(dhd) .Or. vv/dd <= 1.0d-4) &
     & Then
         vmax = dabs (scale*(gd+half*scale*dhd))
         Go To 170
      End If
!
!     G and V are now orthogonal in the subspace spanned by G and D. Hence
!     we generate an orthonormal basis of this subspace such that (D,HV) is
!     negligible or zero, where D and V will be the basis vectors.
!
      ghg = zero
      vhg = zero
      vhv = zero
      Do 140 i = 1, n
         sum = zero
         sumv = zero
         Do 130 j = 1, n
            sum = sum + h (i, j) * g (j)
130      sumv = sumv + h (i, j) * v (j)
         ghg = ghg + sum * g (i)
         vhg = vhg + sumv * g (i)
140   vhv = vhv + sumv * v (i)
      vnorm = dsqrt (vv)
      ghg = ghg / gg
      vhg = vhg / (vnorm*gnorm)
      vhv = vhv / vv
      If (dabs(vhg) <= 0.01d0*dmax1(dabs(ghg), dabs(vhv))) Then
         vmu = ghg - vhv
         wcos = one
         wsin = zero
      Else
         temp = half * (ghg-vhv)
         vmu = temp + dsign (dsqrt(temp**2+vhg**2), temp)
         temp = dsqrt (vmu**2+vhg**2)
         wcos = vmu / temp
         wsin = vhg / temp
      End If
      tempa = wcos / gnorm
      tempb = wsin / vnorm
      tempc = wcos / vnorm
      tempd = wsin / gnorm
      Do 150 i = 1, n
         d (i) = tempa * g (i) + tempb * v (i)
150   v (i) = tempc * v (i) - tempd * g (i)
!
!     The final D is a multiple of the current D, V, D+V or D-V. We make the
!     choice from these possibilities that is optimal.
!
      dlin = wcos * gnorm / rho
      vlin = - wsin * gnorm / rho
      tempa = dabs (dlin) + half * dabs (vmu+vhv)
      tempb = dabs (vlin) + half * dabs (ghg-vmu)
      tempc = halfrt * (dabs(dlin)+dabs(vlin)) + 0.25d0 * dabs &
     & (ghg+vhv)
      If (tempa >= tempb .And. tempa >= tempc) Then
         tempd = dsign (rho, dlin*(vmu+vhv))
         tempv = zero
      Else If (tempb >= tempc) Then
         tempd = zero
         tempv = dsign (rho, vlin*(ghg-vmu))
      Else
         tempd = dsign (halfrt*rho, dlin*(ghg+vhv))
         tempv = dsign (halfrt*rho, vlin*(ghg+vhv))
      End If
      Do 160 i = 1, n
160   d (i) = tempd * d (i) + tempv * v (i)
      vmax = rho * rho * dmax1 (tempa, tempb, tempc)
170   Return
End
!
!     The Chebyquad test problem (Fletcher, 1965) for N = 2,4,6,8.
!
Implicit real * 8 (a-h, o-z)
Dimension x (10), w (10000)
iprint = 2
maxfun = 5000
rhoend = 1.0d-8
Do 30 n = 2, 8, 2
   Do 10 i = 1, n
10 x (i) = dfloat (i) / dfloat (n+1)
   rhobeg = 0.2d0 * x (1)
   Print 20, n
20 Format (/ / 5 x, '******************' / 5 x, 'Results with N =', i2, &
  & / 5 x, '******************')
   Call uobyqa (n, x, rhobeg, rhoend, iprint, maxfun, w)
30 Continue
Stop
End
Subroutine trstep (n, g, h, delta, tol, d, gg, td, tn, w, piv, z, &
& evalue)
   Implicit real * 8 (a-h, o-z)
   Dimension g (*), h (n,*), d (*), gg (*), td (*), tn (*), w (*), piv &
  & (*), z (*)
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
   one = 1.0d0
   two = 2.0d0
   zero = 0.0d0
   delsq = delta * delta
   evalue = zero
   nm = n - 1
   Do 10 i = 1, n
      d (i) = zero
      td (i) = h (i, i)
      Do 10 j = 1, i
10 h (i, j) = h (j, i)
!
!     Apply Householder transformations to obtain a tridiagonal matrix that
!     is similar to H, and put the elements of the Householder vectors in
!     the lower triangular part of H. Further, TD and TN will contain the
!     diagonal and other nonzero elements of the tridiagonal matrix.
!
   Do 80 k = 1, nm
      kp = k + 1
      sum = zero
      If (kp < n) Then
         kpp = kp + 1
         Do 20 i = kpp, n
20       sum = sum + h (i, k) ** 2
      End If
      If (sum == zero) Then
         tn (k) = h (kp, k)
         h (kp, k) = zero
      Else
         temp = h (kp, k)
         tn (k) = dsign (dsqrt(sum+temp*temp), temp)
         h (kp, k) = - sum / (temp+tn(k))
         temp = dsqrt (two/(sum+h(kp, k)**2))
         Do 30 i = kp, n
            w (i) = temp * h (i, k)
            h (i, k) = w (i)
30       z (i) = td (i) * w (i)
         wz = zero
         Do 50 j = kp, nm
            jp = j + 1
            Do 40 i = jp, n
               z (i) = z (i) + h (i, j) * w (j)
40          z (j) = z (j) + h (i, j) * w (i)
50       wz = wz + w (j) * z (j)
         wz = wz + w (n) * z (n)
         Do 70 j = kp, n
            td (j) = td (j) + w (j) * (wz*w(j)-two*z(j))
            If (j < n) Then
               jp = j + 1
               Do 60 i = jp, n
60             h (i, j) = h (i, j) - w (i) * z (j) - w (j) * &
              & (z(i)-wz*w(i))
            End If
70       Continue
      End If
80 Continue
!
!     Form GG by applying the similarity transformation to G.
!
   gsq = zero
   Do 90 i = 1, n
      gg (i) = g (i)
90 gsq = gsq + g (i) ** 2
   gnorm = dsqrt (gsq)
   Do 110 k = 1, nm
      kp = k + 1
      sum = zero
      Do 100 i = kp, n
100   sum = sum + gg (i) * h (i, k)
      Do 110 i = kp, n
110 gg (i) = gg (i) - sum * h (i, k)
!
!     Begin the trust region calculation with a tridiagonal matrix by
!     calculating the norm of H. Then treat the case when H is zero.
!
   hnorm = dabs (td(1)) + dabs (tn(1))
   tdmin = td (1)
   tn (n) = zero
   Do 120 i = 2, n
      temp = dabs (tn(i-1)) + dabs (td(i)) + dabs (tn(i))
      hnorm = dmax1 (hnorm, temp)
120 tdmin = dmin1 (tdmin, td(i))
   If (hnorm == zero) Then
      If (gnorm == zero) Go To 400
      scale = delta / gnorm
      Do 130 i = 1, n
130   d (i) = - scale * gg (i)
      Go To 370
   End If
!
!     Set the initial values of PAR and its bounds.
!
   parl = dmax1 (zero,-tdmin, gnorm/delta-hnorm)
   parlest = parl
   par = parl
   paru = zero
   paruest = zero
   posdef = zero
   iterc = 0
!
!     Calculate the pivots of the Cholesky factorization of (H+PAR*I).
!
140 iterc = iterc + 1
   ksav = 0
   piv (1) = td (1) + par
   k = 1
150 If (piv(k) > zero) Then
      piv (k+1) = td (k+1) + par - tn (k) ** 2 / piv (k)
   Else
      If (piv(k) < zero .Or. tn(k) /= zero) Go To 160
      ksav = k
      piv (k+1) = td (k+1) + par
   End If
   k = k + 1
   If (k < n) Go To 150
   If (piv(k) < zero) Go To 160
   If (piv(k) == zero) ksav = k
!
!     Branch if all the pivots are positive, allowing for the case when
!     G is zero.
!
   If (ksav == 0 .And. gsq > zero) Go To 230
   If (gsq == zero) Then
      If (par == zero) Go To 370
      paru = par
      paruest = par
      If (ksav == 0) Go To 190
   End If
   k = ksav
!
!     Set D to a direction of nonpositive curvature of the given tridiagonal
!     matrix, and thus revise PARLEST.
!
160 d (k) = one
   If (dabs(tn(k)) <= dabs(piv(k))) Then
      dsq = one
      dhd = piv (k)
   Else
      temp = td (k+1) + par
      If (temp <= dabs(piv(k))) Then
         d (k+1) = dsign (one,-tn(k))
         dhd = piv (k) + temp - two * dabs (tn(k))
      Else
         d (k+1) = - tn (k) / temp
         dhd = piv (k) + tn (k) * d (k+1)
      End If
      dsq = one + d (k+1) ** 2
   End If
170 If (k > 1) Then
      k = k - 1
      If (tn(k) /= zero) Then
         d (k) = - tn (k) * d (k+1) / piv (k)
         dsq = dsq + d (k) ** 2
         Go To 170
      End If
      Do 180 i = 1, k
180   d (i) = zero
   End If
   parl = par
   parlest = par - dhd / dsq
!
!     Terminate with D set to a multiple of the current D if the following
!     test suggests that it suitable to do so.
!
190 temp = paruest
   If (gsq == zero) temp = temp * (one-tol)
   If (paruest > zero .And. parlest >= temp) Then
      dtg = zero
      Do 200 i = 1, n
200   dtg = dtg + d (i) * gg (i)
      scale = - dsign (delta/dsqrt(dsq), dtg)
      Do 210 i = 1, n
210   d (i) = scale * d (i)
      Go To 370
   End If
!
!     Pick the value of PAR for the next iteration.
!
220 If (paru == zero) Then
      par = two * parlest + gnorm / delta
   Else
      par = 0.5d0 * (parl+paru)
      par = dmax1 (par, parlest)
   End If
   If (paruest > zero) par = dmin1 (par, paruest)
   Go To 140
!
!     Calculate D for the current PAR in the positive definite case.
!
230 w (1) = - gg (1) / piv (1)
   Do 240 i = 2, n
240 w (i) = (-gg(i)-tn(i-1)*w(i-1)) / piv (i)
   d (n) = w (n)
   Do 250 i = nm, 1, - 1
250 d (i) = w (i) - tn (i) * d (i+1) / piv (i)
!
!     Branch if a Newton-Raphson step is acceptable.
!
   dsq = zero
   wsq = zero
   Do 260 i = 1, n
      dsq = dsq + d (i) ** 2
260 wsq = wsq + piv (i) * w (i) ** 2
   If (par == zero .And. dsq <= delsq) Go To 320
!
!     Make the usual test for acceptability of a full trust region step.
!
   dnorm = dsqrt (dsq)
   phi = one / dnorm - one / delta
   temp = tol * (one+par*dsq/wsq) - dsq * phi * phi
   If (temp >= zero) Then
      scale = delta / dnorm
      Do 270 i = 1, n
270   d (i) = scale * d (i)
      Go To 370
   End If
   If (iterc >= 2 .And. par <= parl) Go To 370
   If (paru > zero .And. par >= paru) Go To 370
!
!     Complete the iteration when PHI is negative.
!
   If (phi < zero) Then
      parlest = par
      If (posdef == one) Then
         If (phi <= phil) Go To 370
         slope = (phi-phil) / (par-parl)
         parlest = par - phi / slope
      End If
      slope = one / gnorm
      If (paru > zero) slope = (phiu-phi) / (paru-par)
      temp = par - phi / slope
      If (paruest > zero) temp = dmin1 (temp, paruest)
      paruest = temp
      posdef = one
      parl = par
      phil = phi
      Go To 220
   End If
!
!     If required, calculate Z for the alternative test for convergence.
!
   If (posdef == zero) Then
      w (1) = one / piv (1)
      Do 280 i = 2, n
         temp = - tn (i-1) * w (i-1)
280   w (i) = (dsign(one, temp)+temp) / piv (i)
      z (n) = w (n)
      Do 290 i = nm, 1, - 1
290   z (i) = w (i) - tn (i) * z (i+1) / piv (i)
      wwsq = zero
      zsq = zero
      dtz = zero
      Do 300 i = 1, n
         wwsq = wwsq + piv (i) * w (i) ** 2
         zsq = zsq + z (i) ** 2
300   dtz = dtz + d (i) * z (i)
!
!     Apply the alternative test for convergence.
!
      tempa = dabs (delsq-dsq)
      tempb = dsqrt (dtz*dtz+tempa*zsq)
      gam = tempa / (dsign(tempb, dtz)+dtz)
      temp = tol * (wsq+par*delsq) - gam * gam * wwsq
      If (temp >= zero) Then
         Do 310 i = 1, n
310      d (i) = d (i) + gam * z (i)
         Go To 370
      End If
      parlest = dmax1 (parlest, par-wwsq/zsq)
   End If
!
!     Complete the iteration when PHI is positive.
!
   slope = one / gnorm
   If (paru > zero) Then
      If (phi >= phiu) Go To 370
      slope = (phiu-phi) / (paru-par)
   End If
   parlest = dmax1 (parlest, par-phi/slope)
   paruest = par
   If (posdef == one) Then
      slope = (phi-phil) / (par-parl)
      paruest = par - phi / slope
   End If
   paru = par
   phiu = phi
   Go To 220
!
!     Set EVALUE to the least eigenvalue of the second derivative matrix if
!     D is a Newton-Raphson step. SHFMAX will be an upper bound on EVALUE.
!
320 shfmin = zero
   pivot = td (1)
   shfmax = pivot
   Do 330 k = 2, n
      pivot = td (k) - tn (k-1) ** 2 / pivot
330 shfmax = dmin1 (shfmax, pivot)
!
!     Find EVALUE by a bisection method, but occasionally SHFMAX may be
!     adjusted by the rule of false position.
!
   ksave = 0
340 shift = 0.5d0 * (shfmin+shfmax)
   k = 1
   temp = td (1) - shift
350 If (temp > zero) Then
      piv (k) = temp
      If (k < n) Then
         temp = td (k+1) - shift - tn (k) ** 2 / temp
         k = k + 1
         Go To 350
      End If
      shfmin = shift
   Else
      If (k < ksave) Go To 360
      If (k == ksave) Then
         If (pivksv == zero) Go To 360
         If (piv(k)-temp < temp-pivksv) Then
            pivksv = temp
            shfmax = shift
         Else
            pivksv = zero
            shfmax = (shift*piv(k)-shfmin*temp) / (piv(k)-temp)
         End If
      Else
         ksave = k
         pivksv = temp
         shfmax = shift
      End If
   End If
   If (shfmin <= 0.99d0*shfmax) Go To 340
360 evalue = shfmin
!
!     Apply the inverse Householder transformations to D.
!
370 nm = n - 1
   Do 390 k = nm, 1, - 1
      kp = k + 1
      sum = zero
      Do 380 i = kp, n
380   sum = sum + d (i) * h (i, k)
      Do 390 i = kp, n
390 d (i) = d (i) - sum * h (i, k)
!
!     Return from the subroutine.
!
400 Return
End
Subroutine uobyqa (n, x, rhobeg, rhoend, iprint, maxfun, w)
   Implicit real * 8 (a-h, o-z)
   Dimension x (*), w (*)
!
!     This subroutine seeks the least value of a function of many variables,
!     by a trust region method that forms quadratic models by interpolation.
!     The algorithm is described in "UOBYQA: unconstrained optimization by
!     quadratic approximation" by M.J.D. Powell, Report DAMTP 2000/NA14,
!     University of Cambridge. The arguments of the subroutine are as follows.
!
!     N must be set to the number of variables and must be at least two.
!     Initial values of the variables must be set in X(1),X(2),...,X(N). They
!       will be changed to the values that give the least calculated F.
!     RHOBEG and RHOEND must be set to the initial and final values of a trust
!       region radius, so both must be positive with RHOEND<=RHOBEG. Typically
!       RHOBEG should be about one tenth of the greatest expected change to a
!       variable, and RHOEND should indicate the accuracy that is required in
!       the final values of the variables.
!     The value of IPRINT should be set to 0, 1, 2 or 3, which controls the
!       amount of printing. Specifically, there is no output if IPRINT=0 and
!       there is output only at the return if IPRINT=1. Otherwise, each new
!       value of RHO is printed, with the best vector of variables so far and
!       the corresponding value of the objective function. Further, each new
!       value of F with its variables are output if IPRINT=3.
!     MAXFUN must be set to an upper bound on the number of calls of CALFUN.
!     The array W will be used for working space. Its length must be at least
!       ( N**4 + 8*N**3 + 23*N**2 + 42*N + max [ 2*N**2 + 4, 18*N ] ) / 4.
!
!     SUBROUTINE CALFUN (N,X,F) must be provided by the user. It must set F to
!     the value of the objective function for the variables X(1),X(2),...,X(N).
!
!     Partition the working space array, so that different parts of it can be
!     treated separately by the subroutine that performs the main calculation.
!
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
   Call uobyqb (n, x, rhobeg, rhoend, iprint, maxfun, npt, w(ixb), &
  & w(ixo), w(ixn), w(ixp), w(ipq), w(ipl), w(ih), w(ig), w(id), &
  & w(ivl), w(iw))
   Return
End
Subroutine uobyqb (n, x, rhobeg, rhoend, iprint, maxfun, npt, xbase, &
& xopt, xnew, xpt, pq, pl, h, g, d, vlag, w)
   Implicit real * 8 (a-h, o-z)
   Dimension x (*), xbase (*), xopt (*), xnew (*), xpt (npt,*), pq (*), &
  & pl (npt,*), h (n,*), g (*), d (*), vlag (*), w (*)
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
   one = 1.0d0
   two = 2.0d0
   zero = 0.0d0
   half = 0.5d0
   tol = 0.01d0
   nnp = n + n + 1
   nptm = npt - 1
   nftest = max0 (maxfun, 1)
!
!     Initialization. NF is the number of function calculations so far.
!
   rho = rhobeg
   rhosq = rho * rho
   nf = 0
   Do 10 i = 1, n
      xbase (i) = x (i)
      Do 10 k = 1, npt
10 xpt (k, i) = zero
   Do 20 k = 1, npt
      Do 20 j = 1, nptm
20 pl (k, j) = zero
!
!     The branch to label 120 obtains a new value of the objective function
!     and then there is a branch back to label 50, because the new function
!     value is needed to form the initial quadratic model. The least function
!     value so far and its index are noted below.
!
30 Do 40 i = 1, n
40 x (i) = xbase (i) + xpt (nf+1, i)
   Go To 120
50 If (nf == 1) Then
      fopt = f
      kopt = nf
      fbase = f
      j = 0
      jswitch = - 1
      ih = n
   Else
      If (f < fopt) Then
         fopt = f
         kopt = nf
      End If
   End If
!
!     Form the gradient and diagonal second derivatives of the initial
!     quadratic model and Lagrange functions.
!
   If (nf <= nnp) Then
      jswitch = - jswitch
      If (jswitch > 0) Then
         If (j >= 1) Then
            ih = ih + j
            If (w(j) < zero) Then
               d (j) = (fsave+f-two*fbase) / rhosq
               pq (j) = (fsave-f) / (two*rho)
               pl (1, ih) = - two / rhosq
               pl (nf-1, j) = half / rho
               pl (nf-1, ih) = one / rhosq
            Else
               pq (j) = (4.0d0*fsave-3.0d0*fbase-f) / (two*rho)
               d (j) = (fbase+f-two*fsave) / rhosq
               pl (1, j) = - 1.5d0 / rho
               pl (1, ih) = one / rhosq
               pl (nf-1, j) = two / rho
               pl (nf-1, ih) = - two / rhosq
            End If
            pq (ih) = d (j)
            pl (nf, j) = - half / rho
            pl (nf, ih) = one / rhosq
         End If
!
!     Pick the shift from XBASE to the next initial interpolation point
!     that provides diagonal second derivatives.
!
         If (j < n) Then
            j = j + 1
            xpt (nf+1, j) = rho
         End If
      Else
         fsave = f
         If (f < fbase) Then
            w (j) = rho
            xpt (nf+1, j) = two * rho
         Else
            w (j) = - rho
            xpt (nf+1, j) = - rho
         End If
      End If
      If (nf < nnp) Go To 30
!
!     Form the off-diagonal second derivatives of the initial quadratic model.
!
      ih = n
      ip = 1
      iq = 2
   End If
   ih = ih + 1
   If (nf > nnp) Then
      temp = one / (w(ip)*w(iq))
      tempa = f - fbase - w (ip) * pq (ip) - w (iq) * pq (iq)
      pq (ih) = (tempa-half*rhosq*(d(ip)+d(iq))) * temp
      pl (1, ih) = temp
      iw = ip + ip
      If (w(ip) < zero) iw = iw + 1
      pl (iw, ih) = - temp
      iw = iq + iq
      If (w(iq) < zero) iw = iw + 1
      pl (iw, ih) = - temp
      pl (nf, ih) = temp
!
!     Pick the shift from XBASE to the next initial interpolation point
!     that provides off-diagonal second derivatives.
!
      ip = ip + 1
   End If
   If (ip == iq) Then
      ih = ih + 1
      ip = 1
      iq = iq + 1
   End If
   If (nf < npt) Then
      xpt (nf+1, ip) = w (ip)
      xpt (nf+1, iq) = w (iq)
      Go To 30
   End If
!
!     Set parameters to begin the iterations for the current RHO.
!
   sixthm = zero
   delta = rho
60 tworsq = (two*rho) ** 2
   rhosq = rho * rho
!
!     Form the gradient of the quadratic model at the trust region centre.
!
70 knew = 0
   ih = n
   Do 80 j = 1, n
      xopt (j) = xpt (kopt, j)
      g (j) = pq (j)
      Do 80 i = 1, j
         ih = ih + 1
         g (i) = g (i) + pq (ih) * xopt (j)
         If (i < j) g (j) = g (j) + pq (ih) * xopt (i)
80 h (i, j) = pq (ih)
!
!     Generate the next trust region step and test its length. Set KNEW
!     to -1 if the purpose of the next F will be to improve conditioning,
!     and also calculate a lower bound on the Hessian term of the model Q.
!
   Call trstep (n, g, h, delta, tol, d, w(1), w(n+1), w(2*n+1), &
  & w(3*n+1), w(4*n+1), w(5*n+1), evalue)
   temp = zero
   Do 90 i = 1, n
90 temp = temp + d (i) ** 2
   dnorm = dmin1 (delta, dsqrt(temp))
   errtol = - one
   If (dnorm < half*rho) Then
      knew = - 1
      errtol = half * evalue * rho * rho
      If (nf <= npt+9) errtol = zero
      Go To 290
   End If
!
!     Calculate the next value of the objective function.
!
100 Do 110 i = 1, n
      xnew (i) = xopt (i) + d (i)
110 x (i) = xbase (i) + xnew (i)
120 If (nf >= nftest) Then
      If (iprint > 0) Print 130
130   Format (/ 4 x, 'Return from UOBYQA because CALFUN has been', ' ca&
     &lled MAXFUN times')
      Go To 420
   End If
   nf = nf + 1
   Call calfun (n, x, f)
   If (iprint == 3) Then
      Print 140, nf, f, (x(i), i=1, n)
140   Format (/ 4 x, 'Function number', i6, '    F =', 1 pd18.10, '    &
     &The corresponding X is:' / (2 x, 5d15.6))
   End If
   If (nf <= npt) Go To 50
   If (knew ==-1) Go To 420
!
!     Use the quadratic model to predict the change in F due to the step D,
!     and find the values of the Lagrange functions at the new point.
!
   vquad = zero
   ih = n
   Do 150 j = 1, n
      w (j) = d (j)
      vquad = vquad + w (j) * pq (j)
      Do 150 i = 1, j
         ih = ih + 1
         w (ih) = d (i) * xnew (j) + d (j) * xopt (i)
         If (i == j) w (ih) = half * w (ih)
150 vquad = vquad + w (ih) * pq (ih)
   Do 170 k = 1, npt
      temp = zero
      Do 160 j = 1, nptm
160   temp = temp + w (j) * pl (k, j)
170 vlag (k) = temp
   vlag (kopt) = vlag (kopt) + one
!
!     Update SIXTHM, which is a lower bound on one sixth of the greatest
!     third derivative of F.
!
   diff = f - fopt - vquad
   sum = zero
   Do 190 k = 1, npt
      temp = zero
      Do 180 i = 1, n
180   temp = temp + (xpt(k, i)-xnew(i)) ** 2
      temp = dsqrt (temp)
190 sum = sum + dabs (temp*temp*temp*vlag(k))
   sixthm = dmax1 (sixthm, dabs(diff)/sum)
!
!     Update FOPT and XOPT if the new F is the least value of the objective
!     function so far. Then branch if D is not a trust region step.
!
   fsave = fopt
   If (f < fopt) Then
      fopt = f
      Do 200 i = 1, n
200   xopt (i) = xnew (i)
   End If
   ksave = knew
   If (knew > 0) Go To 240
!
!     Pick the next value of DELTA after a trust region step.
!
   If (vquad >= zero) Then
      If (iprint > 0) Print 210
210   Format (/ 4 x, 'Return from UOBYQA because a trust', ' region ste&
     &p has failed to reduce Q')
      Go To 420
   End If
   ratio = (f-fsave) / vquad
   If (ratio <= 0.1d0) Then
      delta = half * dnorm
   Else If (ratio <= 0.7d0) Then
      delta = dmax1 (half*delta, dnorm)
   Else
      delta = dmax1 (delta, 1.25d0*dnorm, dnorm+rho)
   End If
   If (delta <= 1.5d0*rho) delta = rho
!
!     Set KNEW to the index of the next interpolation point to be deleted.
!
   ktemp = 0
   detrat = zero
   If (f >= fsave) Then
      ktemp = kopt
      detrat = one
   End If
   Do 230 k = 1, npt
      sum = zero
      Do 220 i = 1, n
220   sum = sum + (xpt(k, i)-xopt(i)) ** 2
      temp = dabs (vlag(k))
      If (sum > rhosq) temp = temp * (sum/rhosq) ** 1.5d0
      If (temp > detrat .And. k /= ktemp) Then
         detrat = temp
         ddknew = sum
         knew = k
      End If
230 Continue
   If (knew == 0) Go To 290
!
!     Replace the interpolation point that has index KNEW by the point XNEW,
!     and also update the Lagrange functions and the quadratic model.
!
240 Do 250 i = 1, n
250 xpt (knew, i) = xnew (i)
   temp = one / vlag (knew)
   Do 260 j = 1, nptm
      pl (knew, j) = temp * pl (knew, j)
260 pq (j) = pq (j) + diff * pl (knew, j)
   Do 280 k = 1, npt
      If (k /= knew) Then
         temp = vlag (k)
         Do 270 j = 1, nptm
270      pl (k, j) = pl (k, j) - temp * pl (knew, j)
      End If
280 Continue
!
!     Update KOPT if F is the least calculated value of the objective
!     function. Then branch for another trust region calculation. The
!     case KSAVE>0 indicates that a model step has just been taken.
!
   If (f < fsave) Then
      kopt = knew
      Go To 70
   End If
   If (ksave > 0) Go To 70
   If (dnorm > two*rho) Go To 70
   If (ddknew > tworsq) Go To 70
!
!     Alternatively, find out if the interpolation points are close
!     enough to the best point so far.
!
290 Do 300 k = 1, npt
      w (k) = zero
      Do 300 i = 1, n
300 w (k) = w (k) + (xpt(k, i)-xopt(i)) ** 2
310 knew = - 1
   distest = tworsq
   Do 320 k = 1, npt
      If (w(k) > distest) Then
         knew = k
         distest = w (k)
      End If
320 Continue
!
!     If a point is sufficiently far away, then set the gradient and Hessian
!     of its Lagrange function at the centre of the trust region, and find
!     half the sum of squares of components of the Hessian.
!
   If (knew > 0) Then
      ih = n
      sumh = zero
      Do 340 j = 1, n
         g (j) = pl (knew, j)
         Do 330 i = 1, j
            ih = ih + 1
            temp = pl (knew, ih)
            g (j) = g (j) + temp * xopt (i)
            If (i < j) Then
               g (i) = g (i) + temp * xopt (j)
               sumh = sumh + temp * temp
            End If
330      h (i, j) = temp
340   sumh = sumh + half * temp * temp
!
!     If ERRTOL is positive, test whether to replace the interpolation point
!     with index KNEW, using a bound on the maximum modulus of its Lagrange
!     function in the trust region.
!
      If (errtol > zero) Then
         w (knew) = zero
         sumg = zero
         Do 350 i = 1, n
350      sumg = sumg + g (i) ** 2
         estim = rho * (dsqrt(sumg)+rho*dsqrt(half*sumh))
         wmult = sixthm * distest ** 1.5d0
         If (wmult*estim <= errtol) Go To 310
      End If
!
!     If the KNEW-th point may be replaced, then pick a D that gives a large
!     value of the modulus of its Lagrange function within the trust region.
!     Here the vector XNEW is used as temporary working space.
!
      Call lagmax (n, g, h, rho, d, xnew, vmax)
      If (errtol > zero) Then
         If (wmult*vmax <= errtol) Go To 310
      End If
      Go To 100
   End If
   If (dnorm > rho) Go To 70
!
!     Prepare to reduce RHO by shifting XBASE to the best point so far,
!     and make the corresponding changes to the gradients of the Lagrange
!     functions and the quadratic model.
!
   If (rho > rhoend) Then
      ih = n
      Do 380 j = 1, n
         xbase (j) = xbase (j) + xopt (j)
         Do 360 k = 1, npt
360      xpt (k, j) = xpt (k, j) - xopt (j)
         Do 380 i = 1, j
            ih = ih + 1
            pq (i) = pq (i) + pq (ih) * xopt (j)
            If (i < j) Then
               pq (j) = pq (j) + pq (ih) * xopt (i)
               Do 370 k = 1, npt
370            pl (k, j) = pl (k, j) + pl (k, ih) * xopt (i)
            End If
            Do 380 k = 1, npt
380   pl (k, i) = pl (k, i) + pl (k, ih) * xopt (j)
!
!     Pick the next values of RHO and DELTA.
!
      delta = half * rho
      ratio = rho / rhoend
      If (ratio <= 16.0d0) Then
         rho = rhoend
      Else If (ratio <= 250.0d0) Then
         rho = dsqrt (ratio) * rhoend
      Else
         rho = 0.1d0 * rho
      End If
      delta = dmax1 (delta, rho)
      If (iprint >= 2) Then
         If (iprint >= 3) Print 390
390      Format (5 x)
         Print 400, rho, nf
400      Format (/ 4 x, 'New RHO =', 1 pd11.4, 5 x, 'Number of', ' func&
        &tion values =', i6)
         Print 410, fopt, (xbase(i), i=1, n)
410      Format (4 x, 'Least value of F =', 1 pd23.15, 9 x, 'The corres&
        &ponding X is:'/(2 x, 5d15.6))
      End If
      Go To 60
   End If
!
!     Return from the calculation, after another Newton-Raphson step, if
!     it is too short to have been tried before.
!
   If (errtol >= zero) Go To 100
420 If (fopt <= f) Then
      Do 430 i = 1, n
430   x (i) = xbase (i) + xopt (i)
      f = fopt
   End If
   If (iprint >= 1) Then
      Print 440, nf
440   Format (/ 4 x, 'At the return from UOBYQA', 5 x, 'Number of funct&
     &ion values =', i6)
      Print 410, f, (x(i), i=1, n)
   End If
   Return
End
