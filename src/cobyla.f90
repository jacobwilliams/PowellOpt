    module cobyla_module

    use kind_module, only: wp
	
    implicit real(wp) (a-h, o-z)

    private

    abstract interface
        subroutine func (N,M,X,F,CON)  !! CALCFC interface
        import :: wp
        implicit none
        integer :: n
        integer :: m
        real(wp) :: x(*)
        real(wp) :: f
        real(wp) :: con(*)
        end subroutine func
    end interface

    public :: cobyla
    public :: cobyla_test

    contains

Subroutine cobyla (n, m, x, rhobeg, rhoend, iprint, maxfun, w, iact, calcfc)
      Dimension x (*), w (*), iact (*)
      procedure(func) :: calcfc
!
!     This subroutine minimizes an objective function F(X) subject to M
!     inequality constraints on X, where X is a vector of variables that has
!     N components. The algorithm employs linear approximations to the
!     objective and constraint functions, the approximations being formed by
!     linear interpolation at N+1 points in the space of the variables.
!     We regard these interpolation points as vertices of a simplex. The
!     parameter RHO controls the size of the simplex and it is reduced
!     automatically from RHOBEG to RHOEND. For each RHO the subroutine tries
!     to achieve a good vector of variables for the current size, and then
!     RHO is reduced until the value RHOEND is reached. Therefore RHOBEG and
!     RHOEND should be set to reasonable initial changes to and the required
!     accuracy in the variables respectively, but this accuracy should be
!     viewed as a subject for experimentation because it is not guaranteed.
!     The subroutine has an advantage over many of its competitors, however,
!     which is that it treats each constraint individually when calculating
!     a change to the variables, instead of lumping the constraints together
!     into a single penalty function. The name of the subroutine is derived
!     from the phrase Constrained Optimization BY Linear Approximations.
!
!     The user must set the values of N, M, RHOBEG and RHOEND, and must
!     provide an initial vector of variables in X. Further, the value of
!     IPRINT should be set to 0, 1, 2 or 3, which controls the amount of
!     printing during the calculation. Specifically, there is no output if
!     IPRINT=0 and there is output only at the end of the calculation if
!     IPRINT=1. Otherwise each new value of RHO and SIGMA is printed.
!     Further, the vector of variables and some function information are
!     given either when RHO is reduced or when each new value of F(X) is
!     computed in the cases IPRINT=2 or IPRINT=3 respectively. Here SIGMA
!     is a penalty parameter, it being assumed that a change to X is an
!     improvement if it reduces the merit function
!                F(X)+SIGMA*MAX(0.0,-C1(X),-C2(X),...,-CM(X)),
!     where C1,C2,...,CM denote the constraint functions that should become
!     nonnegative eventually, at least to the precision of RHOEND. In the
!     printed output the displayed term that is multiplied by SIGMA is
!     called MAXCV, which stands for 'MAXimum Constraint Violation'. The
!     argument MAXFUN is an integer variable that must be set by the user to a
!     limit on the number of calls of CALCFC, the purpose of this routine being
!     given below. The value of MAXFUN will be altered to the number of calls
!     of CALCFC that are made. The arguments W and IACT provide real and
!     integer arrays that are used as working space. Their lengths must be at
!     least N*(3*N+2*M+11)+4*M+6 and M+1 respectively.
!
!     In order to define the objective and constraint functions, we require
!     a subroutine that has the name and arguments
!                SUBROUTINE CALCFC (N,M,X,F,CON)
!                DIMENSION X(*),CON(*)  .
!     The values of N and M are fixed and have been defined already, while
!     X is now the current vector of variables. The subroutine should return
!     the objective and constraint functions at X in F and CON(1),CON(2),
!     ...,CON(M). Note that we are trying to adjust X so that F(X) is as
!     small as possible subject to the constraint functions being nonnegative.
!
!     Partition the working space array W to provide the storage that is needed
!     for the main calculation.
!
      mpp = m + 2
      icon = 1
      isim = icon + mpp
      isimi = isim + n * n + n
      idatm = isimi + n * n
      ia = idatm + n * mpp + mpp
      ivsig = ia + m * n + n
      iveta = ivsig + n
      isigb = iveta + n
      idx = isigb + n
      iwork = idx + n
      
      Call cobylb (n, m, mpp, x, rhobeg, rhoend, iprint, maxfun, &
       w(icon), w(isim), w(isimi), w(idatm), w(ia), w(ivsig), w(iveta), &
       w(isigb), w(idx), w(iwork), iact, calcfc)

End Subroutine cobyla

Subroutine cobylb (n, m, mpp, x, rhobeg, rhoend, iprint, maxfun, con, &
   sim, simi, datmat, a, vsig, veta, sigbar, dx, w, iact, calcfc)
      Dimension x (*), con (*), sim (n,*), simi (n,*), datmat (mpp,*), &
       a (n,*), vsig (*), veta (*), sigbar (*), dx (*), w (*), iact (*)
      procedure(func) :: calcfc
!
!     Set the initial values of some parameters. The last column of SIM holds
!     the optimal vertex of the current simplex, and the preceding N columns
!     hold the displacements from the optimal vertex to the other vertices.
!     Further, SIMI holds the inverse of the matrix that is contained in the
!     first N columns of SIM.
!
      iptem = min (n, 5)
      iptemp = iptem + 1
      np = n + 1
      mp = m + 1
      alpha = 0.25_wp
      beta = 2.1_wp
      gamma = 0.5_wp
      delta = 1.1_wp
      rho = rhobeg
      parmu = 0.0_wp
      If (iprint >= 2) Print 10, rho
10    Format (/ 3 x, 'The initial value of RHO is', 1 pe13.6, 2 x,&
      'and PARMU is set to zero.')
      nfvals = 0
      temp = 1.0_wp / rho
      Do 30 i = 1, n
         sim (i, np) = x (i)
         Do 20 j = 1, n
            sim (i, j) = 0.0_wp
20       simi (i, j) = 0.0_wp
         sim (i, i) = rho
30    simi (i, i) = temp
      jdrop = np
      ibrnch = 0
!
!     Make the next call of the user-supplied subroutine CALCFC. These
!     instructions are also used for calling CALCFC during the iterations of
!     the algorithm.
!
40    If (nfvals >= maxfun .And. nfvals > 0) Then
         If (iprint >= 1) Print 50
50       Format (/ 3 x, 'Return from subroutine COBYLA because the ',&
         'MAXFUN limit has been reached.')
         Go To 600
      End If
      nfvals = nfvals + 1
      Call calcfc (n, m, x, f, con)
      resmax = 0.0_wp
      If (m > 0) Then
         Do 60 k = 1, m
60       resmax = max (resmax,-con(k))
      End If
      If (nfvals == iprint-1 .Or. iprint == 3) Then
         Print 70, nfvals, f, resmax, (x(i), i=1, iptem)
70       Format (/ 3 x, 'NFVALS =', i5, 3 x, 'F =', 1 pe13.6, 4 x,&
         'MAXCV =', 1 pe13.6 / 3 x, 'X =', 1 pe13.6, 1 p4e15.6)
         If (iptem < n) Print 80, (x(i), i=iptemp, n)
80       Format (1 pe19.6, 1 p4e15.6)
      End If
      con (mp) = f
      con (mpp) = resmax
      If (ibrnch == 1) Go To 440
!
!     Set the recently calculated function values in a column of DATMAT. This
!     array has a column for each vertex of the current simplex, the entries of
!     each column being the values of the constraint functions (if any)
!     followed by the objective function and the greatest constraint violation
!     at the vertex.
!
      Do 90 k = 1, mpp
90    datmat (k, jdrop) = con (k)
      If (nfvals > np) Go To 130
!
!     Exchange the new vertex of the initial simplex with the optimal vertex if
!     necessary. Then, if the initial simplex is not complete, pick its next
!     vertex and calculate the function values there.
!
      If (jdrop <= n) Then
         If (datmat(mp, np) <= f) Then
            x (jdrop) = sim (jdrop, np)
         Else
            sim (jdrop, np) = x (jdrop)
            Do 100 k = 1, mpp
               datmat (k, jdrop) = datmat (k, np)
100         datmat (k, np) = con (k)
            Do 120 k = 1, jdrop
               sim (jdrop, k) = - rho
               temp = 0.0_wp
               Do 110 i = k, jdrop
110            temp = temp - simi (i, k)
120         simi (jdrop, k) = temp
         End If
      End If
      If (nfvals <= n) Then
         jdrop = nfvals
         x (jdrop) = x (jdrop) + rho
         Go To 40
      End If
130   ibrnch = 1
!
!     Identify the optimal vertex of the current simplex.
!
140   phimin = datmat (mp, np) + parmu * datmat (mpp, np)
      nbest = np
      Do 150 j = 1, n
         temp = datmat (mp, j) + parmu * datmat (mpp, j)
         If (temp < phimin) Then
            nbest = j
            phimin = temp
         Else If (temp == phimin .And. parmu == 0.0_wp) Then
            If (datmat(mpp, j) < datmat(mpp, nbest)) nbest = j
         End If
150   Continue
!
!     Switch the best vertex into pole position if it is not there already,
!     and also update SIM, SIMI and DATMAT.
!
      If (nbest <= n) Then
         Do 160 i = 1, mpp
            temp = datmat (i, np)
            datmat (i, np) = datmat (i, nbest)
160      datmat (i, nbest) = temp
         Do 180 i = 1, n
            temp = sim (i, nbest)
            sim (i, nbest) = 0.0_wp
            sim (i, np) = sim (i, np) + temp
            tempa = 0.0_wp
            Do 170 k = 1, n
               sim (i, k) = sim (i, k) - temp
170         tempa = tempa - simi (k, i)
180      simi (nbest, i) = tempa
      End If
!
!     Make an error return if SIGI is a poor approximation to the inverse of
!     the leading N by N submatrix of SIG.
!
      error = 0.0_wp
      Do 200 i = 1, n
         Do 200 j = 1, n
            temp = 0.0_wp
            If (i == j) temp = temp - 1.0_wp
            Do 190 k = 1, n
190         temp = temp + simi (i, k) * sim (k, j)
200   error = max (error, Abs(temp))
      If (error > 0.1_wp) Then
         If (iprint >= 1) Print 210
210      Format (/ 3 x, 'Return from subroutine COBYLA because ',&
         'rounding errors are becoming damaging.')
         Go To 600
      End If
!
!     Calculate the coefficients of the linear approximations to the objective
!     and constraint functions, placing minus the objective function gradient
!     after the constraint gradients in the array A. The vector W is used for
!     working space.
!
      Do 240 k = 1, mp
         con (k) = - datmat (k, np)
         Do 220 j = 1, n
220      w (j) = datmat (k, j) + con (k)
         Do 240 i = 1, n
            temp = 0.0_wp
            Do 230 j = 1, n
230         temp = temp + w (j) * simi (j, i)
            If (k == mp) temp = - temp
240   a (i, k) = temp
!
!     Calculate the values of sigma and eta, and set IFLAG=0 if the current
!     simplex is not acceptable.
!
      iflag = 1
      parsig = alpha * rho
      pareta = beta * rho
      Do 260 j = 1, n
         wsig = 0.0_wp
         weta = 0.0_wp
         Do 250 i = 1, n
            wsig = wsig + simi (j, i) ** 2
250      weta = weta + sim (i, j) ** 2
         vsig (j) = 1.0_wp / Sqrt (wsig)
         veta (j) = Sqrt (weta)
         If (vsig(j) < parsig .Or. veta(j) > pareta) iflag = 0
260   Continue
!
!     If a new vertex is needed to improve acceptability, then decide which
!     vertex to drop from the simplex.
!
      If (ibrnch == 1 .Or. iflag == 1) Go To 370
      jdrop = 0
      temp = pareta
      Do 270 j = 1, n
         If (veta(j) > temp) Then
            jdrop = j
            temp = veta (j)
         End If
270   Continue
      If (jdrop == 0) Then
         Do 280 j = 1, n
            If (vsig(j) < temp) Then
               jdrop = j
               temp = vsig (j)
            End If
280      Continue
      End If
!
!     Calculate the step to the new vertex and its sign.
!
      temp = gamma * rho * vsig (jdrop)
      Do 290 i = 1, n
290   dx (i) = temp * simi (jdrop, i)
      cvmaxp = 0.0_wp
      cvmaxm = 0.0_wp
      Do 310 k = 1, mp
         sum = 0.0_wp
         Do 300 i = 1, n
300      sum = sum + a (i, k) * dx (i)
         If (k < mp) Then
            temp = datmat (k, np)
            cvmaxp = max (cvmaxp,-sum-temp)
            cvmaxm = max (cvmaxm, sum-temp)
         End If
310   Continue
      dxsign = 1.0_wp
      If (parmu*(cvmaxp-cvmaxm) > sum+sum) dxsign = - 1.0_wp
!
!     Update the elements of SIM and SIMI, and set the next X.
!
      temp = 0.0_wp
      Do 320 i = 1, n
         dx (i) = dxsign * dx (i)
         sim (i, jdrop) = dx (i)
320   temp = temp + simi (jdrop, i) * dx (i)
      Do 330 i = 1, n
330   simi (jdrop, i) = simi (jdrop, i) / temp
      Do 360 j = 1, n
         If (j /= jdrop) Then
            temp = 0.0_wp
            Do 340 i = 1, n
340         temp = temp + simi (j, i) * dx (i)
            Do 350 i = 1, n
350         simi (j, i) = simi (j, i) - temp * simi (jdrop, i)
         End If
360   x (j) = sim (j, np) + dx (j)
      Go To 40
!
!     Calculate DX=x(*)-x(0). Branch if the length of DX is less than 0.5*RHO.
!
370   iz = 1
      izdota = iz + n * n
      ivmc = izdota + n
      isdirn = ivmc + mp
      idxnew = isdirn + n
      ivmd = idxnew + n
      Call trstlp (n, m, a, con, rho, dx, ifull, iact, w(iz), &
       w(izdota), w(ivmc), w(isdirn), w(idxnew), w(ivmd))
      If (ifull == 0) Then
         temp = 0.0_wp
         Do 380 i = 1, n
380      temp = temp + dx (i) ** 2
         If (temp < 0.25*rho*rho) Then
            ibrnch = 1
            Go To 550
         End If
      End If
!
!     Predict the change to F and the new maximum constraint violation if the
!     variables are altered from x(0) to x(0)+DX.
!
      resnew = 0.0_wp
      con (mp) = 0.0_wp
      Do 400 k = 1, mp
         sum = con (k)
         Do 390 i = 1, n
390      sum = sum - a (i, k) * dx (i)
         If (k < mp) resnew = max (resnew, sum)
400   Continue
!
!     Increase PARMU if necessary and branch back if this change alters the
!     optimal vertex. Otherwise PREREM and PREREC will be set to the predicted
!     reductions in the merit function and the maximum constraint violation
!     respectively.
!
      barmu = 0.0_wp
      prerec = datmat (mpp, np) - resnew
      If (prerec > 0.0_wp) barmu = sum / prerec
      If (parmu < 1.5_wp*barmu) Then
         parmu = 2.0_wp * barmu
         If (iprint >= 2) Print 410, parmu
410      Format (/ 3 x, 'Increase in PARMU to', 1 pe13.6)
         phi = datmat (mp, np) + parmu * datmat (mpp, np)
         Do 420 j = 1, n
            temp = datmat (mp, j) + parmu * datmat (mpp, j)
            If (temp < phi) Go To 140
            If (temp == phi .And. parmu == 0.0_wp) Then
               If (datmat(mpp, j) < datmat(mpp, np)) Go To 140
            End If
420      Continue
      End If
      prerem = parmu * prerec - sum
!
!     Calculate the constraint and objective functions at x(*). Then find the
!     actual reduction in the merit function.
!
      Do 430 i = 1, n
430   x (i) = sim (i, np) + dx (i)
      ibrnch = 1
      Go To 40
440   vmold = datmat (mp, np) + parmu * datmat (mpp, np)
      vmnew = f + parmu * resmax
      trured = vmold - vmnew
      If (parmu == 0.0_wp .And. f == datmat(mp, np)) Then
         prerem = prerec
         trured = datmat (mpp, np) - resmax
      End If
!
!     Begin the operations that decide whether x(*) should replace one of the
!     vertices of the current simplex, the change being mandatory if TRURED is
!     positive. Firstly, JDROP is set to the index of the vertex that is to be
!     replaced.
!
      ratio = 0.0_wp
      If (trured <= 0.0_wp) ratio = 1.0_wp
      jdrop = 0
      Do 460 j = 1, n
         temp = 0.0_wp
         Do 450 i = 1, n
450      temp = temp + simi (j, i) * dx (i)
         temp = Abs (temp)
         If (temp > ratio) Then
            jdrop = j
            ratio = temp
         End If
460   sigbar (j) = temp * vsig (j)
!
!     Calculate the value of ell.
!
      edgmax = delta * rho
      l = 0
      Do 480 j = 1, n
         If (sigbar(j) >= parsig .Or. sigbar(j) >= vsig(j)) Then
            temp = veta (j)
            If (trured > 0.0_wp) Then
               temp = 0.0_wp
               Do 470 i = 1, n
470            temp = temp + (dx(i)-sim(i, j)) ** 2
               temp = Sqrt (temp)
            End If
            If (temp > edgmax) Then
               l = j
               edgmax = temp
            End If
         End If
480   Continue
      If (l > 0) jdrop = l
      If (jdrop == 0) Go To 550
!
!     Revise the simplex by updating the elements of SIM, SIMI and DATMAT.
!
      temp = 0.0_wp
      Do 490 i = 1, n
         sim (i, jdrop) = dx (i)
490   temp = temp + simi (jdrop, i) * dx (i)
      Do 500 i = 1, n
500   simi (jdrop, i) = simi (jdrop, i) / temp
      Do 530 j = 1, n
         If (j /= jdrop) Then
            temp = 0.0_wp
            Do 510 i = 1, n
510         temp = temp + simi (j, i) * dx (i)
            Do 520 i = 1, n
520         simi (j, i) = simi (j, i) - temp * simi (jdrop, i)
         End If
530   Continue
      Do 540 k = 1, mpp
540   datmat (k, jdrop) = con (k)
!
!     Branch back for further iterations with the current RHO.
!
      If (trured > 0.0_wp .And. trured >= 0.1_wp*prerem) Go To 140
550   If (iflag == 0) Then
         ibrnch = 0
         Go To 140
      End If
!
!     Otherwise reduce RHO if it is not at its least value and reset PARMU.
!
      If (rho > rhoend) Then
         rho = 0.5_wp * rho
         If (rho <= 1.5_wp*rhoend) rho = rhoend
         If (parmu > 0.0_wp) Then
            denom = 0.0_wp
            Do 570 k = 1, mp
               cmin = datmat (k, np)
               cmax = cmin
               Do 560 i = 1, n
                  cmin = min (cmin, datmat(k, i))
560            cmax = max (cmax, datmat(k, i))
               If (k <= m .And. cmin < 0.5_wp*cmax) Then
                  temp = max (cmax, 0.0_wp) - cmin
                  If (denom <= 0.0_wp) Then
                     denom = temp
                  Else
                     denom = min (denom, temp)
                  End If
               End If
570         Continue
            If (denom == 0.0_wp) Then
               parmu = 0.0_wp
            Else If (cmax-cmin < parmu*denom) Then
               parmu = (cmax-cmin) / denom
            End If
         End If
         If (iprint >= 2) Print 580, rho, parmu
580      Format (/ 3 x, 'Reduction in RHO to', 1 pe13.6, '  and PARMU =', 1 pe13.6)
         If (iprint == 2) Then
            Print 70, nfvals, datmat (mp, np), datmat (mpp, np), &
             (sim(i, np), i=1, iptem)
            If (iptem < n) Print 80, (x(i), i=iptemp, n)
         End If
         Go To 140
      End If
!
!     Return the best calculated values of the variables.
!
      If (iprint >= 1) Print 590
590   Format (/ 3 x, 'Normal return from subroutine COBYLA')
      If (ifull == 1) Go To 620
600   Do 610 i = 1, n
610   x (i) = sim (i, np)
      f = datmat (mp, np)
      resmax = datmat (mpp, np)
620   If (iprint >= 1) Then
         Print 70, nfvals, f, resmax, (x(i), i=1, iptem)
         If (iptem < n) Print 80, (x(i), i=iptemp, n)
      End If
      maxfun = nfvals

End Subroutine cobylb

Subroutine trstlp (n, m, a, b, rho, dx, ifull, iact, z, zdota, vmultc, &
  sdirn, dxnew, vmultd)
   Dimension a (n,*), b (*), dx (*), iact (*), z (n,*), zdota (*), &
    vmultc (*), sdirn (*), dxnew (*), vmultd (*)
!
!     This subroutine calculates an N-component vector DX by applying the
!     following two stages. In the first stage, DX is set to the shortest
!     vector that minimizes the greatest violation of the constraints
!       A(1,K)*DX(1)+A(2,K)*DX(2)+...+A(N,K)*DX(N) .GE. B(K), K=2,3,...,M,
!     subject to the Euclidean length of DX being at most RHO. If its length is
!     strictly less than RHO, then we use the resultant freedom in DX to
!     minimize the objective function
!              -A(1,M+1)*DX(1)-A(2,M+1)*DX(2)-...-A(N,M+1)*DX(N)
!     subject to no increase in any greatest constraint violation. This
!     notation allows the gradient of the objective function to be regarded as
!     the gradient of a constraint. Therefore the two stages are distinguished
!     by MCON .EQ. M and MCON .GT. M respectively. It is possible that a
!     degeneracy may prevent DX from attaining the target length RHO. Then the
!     value IFULL=0 would be set, but usually IFULL=1 on return.
!
!     In general NACT is the number of constraints in the active set and
!     IACT(1),...,IACT(NACT) are their indices, while the remainder of IACT
!     contains a permutation of the remaining constraint indices. Further, Z is
!     an orthogonal matrix whose first NACT columns can be regarded as the
!     result of Gram-Schmidt applied to the active constraint gradients. For
!     J=1,2,...,NACT, the number ZDOTA(J) is the scalar product of the J-th
!     column of Z with the gradient of the J-th active constraint. DX is the
!     current vector of variables and here the residuals of the active
!     constraints should be zero. Further, the active constraints have
!     nonnegative Lagrange multipliers that are held at the beginning of
!     VMULTC. The remainder of this vector holds the residuals of the inactive
!     constraints at DX, the ordering of the components of VMULTC being in
!     agreement with the permutation of the indices of the constraints that is
!     in IACT. All these residuals are nonnegative, which is achieved by the
!     shift RESMAX that makes the least residual zero.
!
!     Initialize Z and some other variables. The value of RESMAX will be
!     appropriate to DX=0, while ICON will be the index of a most violated
!     constraint if RESMAX is positive. Usually during the first stage the
!     vector SDIRN gives a search direction that reduces all the active
!     constraint violations by one simultaneously.
!
   ifull = 1
   mcon = m
   nact = 0
   resmax = 0.0_wp
   Do 20 i = 1, n
      Do 10 j = 1, n
10    z (i, j) = 0.0_wp
      z (i, i) = 1.0_wp
20 dx (i) = 0.0_wp
   If (m >= 1) Then
      Do 30 k = 1, m
         If (b(k) > resmax) Then
            resmax = b (k)
            icon = k
         End If
30    Continue
      Do 40 k = 1, m
         iact (k) = k
40    vmultc (k) = resmax - b (k)
   End If
   If (resmax == 0.0_wp) Go To 480
   Do 50 i = 1, n
50 sdirn (i) = 0.0_wp
!
!     End the current stage of the calculation if 3 consecutive iterations
!     have either failed to reduce the best calculated value of the objective
!     function or to increase the number of active constraints since the best
!     value was calculated. This strategy prevents cycling, but there is a
!     remote possibility that it will cause premature termination.
!
60 optold = 0.0_wp
   icount = 0
70 If (mcon == m) Then
      optnew = resmax
   Else
      optnew = 0.0_wp
      Do 80 i = 1, n
80    optnew = optnew - dx (i) * a (i, mcon)
   End If
   If (icount == 0 .Or. optnew < optold) Then
      optold = optnew
      nactx = nact
      icount = 3
   Else If (nact > nactx) Then
      nactx = nact
      icount = 3
   Else
      icount = icount - 1
      If (icount == 0) Go To 490
   End If
!
!     If ICON exceeds NACT, then we add the constraint with index IACT(ICON) to
!     the active set. Apply Givens rotations so that the last N-NACT-1 columns
!     of Z are orthogonal to the gradient of the new constraint, a scalar
!     product being set to zero if its nonzero value could be due to computer
!     rounding errors. The array DXNEW is used for working space.
!
   If (icon <= nact) Go To 260
   kk = iact (icon)
   Do 90 i = 1, n
90 dxnew (i) = a (i, kk)
   tot = 0.0_wp
   k = n
100 If (k > nact) Then
      sp = 0.0_wp
      spabs = 0.0_wp
      Do 110 i = 1, n
         temp = z (i, k) * dxnew (i)
         sp = sp + temp
110   spabs = spabs + Abs (temp)
      acca = spabs + 0.1_wp * Abs (sp)
      accb = spabs + 0.2_wp * Abs (sp)
      If (spabs >= acca .Or. acca >= accb) sp = 0.0_wp
      If (tot == 0.0_wp) Then
         tot = sp
      Else
         kp = k + 1
         temp = Sqrt (sp*sp+tot*tot)
         alpha = sp / temp
         beta = tot / temp
         tot = temp
         Do 120 i = 1, n
            temp = alpha * z (i, k) + beta * z (i, kp)
            z (i, kp) = alpha * z (i, kp) - beta * z (i, k)
120      z (i, k) = temp
      End If
      k = k - 1
      Go To 100
   End If
!
!     Add the new constraint if this can be done without a deletion from the
!     active set.
!
   If (tot /= 0.0_wp) Then
      nact = nact + 1
      zdota (nact) = tot
      vmultc (icon) = vmultc (nact)
      vmultc (nact) = 0.0_wp
      Go To 210
   End If
!
!     The next instruction is reached if a deletion has to be made from the
!     active set in order to make room for the new active constraint, because
!     the new constraint gradient is a linear combination of the gradients of
!     the old active constraints. Set the elements of VMULTD to the multipliers
!     of the linear combination. Further, set IOUT to the index of the
!     constraint to be deleted, but branch if no suitable index can be found.
!
   ratio = - 1.0_wp
   k = nact
130 zdotv = 0.0_wp
   zdvabs = 0.0_wp
   Do 140 i = 1, n
      temp = z (i, k) * dxnew (i)
      zdotv = zdotv + temp
140 zdvabs = zdvabs + Abs (temp)
   acca = zdvabs + 0.1_wp * Abs (zdotv)
   accb = zdvabs + 0.2_wp * Abs (zdotv)
   If (zdvabs < acca .And. acca < accb) Then
      temp = zdotv / zdota (k)
      If (temp > 0.0_wp .And. iact(k) <= m) Then
         tempa = vmultc (k) / temp
         If (ratio < 0.0_wp .Or. tempa < ratio) Then
            ratio = tempa
            iout = k
         End If
      End If
      If (k >= 2) Then
         kw = iact (k)
         Do 150 i = 1, n
150      dxnew (i) = dxnew (i) - temp * a (i, kw)
      End If
      vmultd (k) = temp
   Else
      vmultd (k) = 0.0_wp
   End If
   k = k - 1
   If (k > 0) Go To 130
   If (ratio < 0.0_wp) Go To 490
!
!     Revise the Lagrange multipliers and reorder the active constraints so
!     that the one to be replaced is at the end of the list. Also calculate the
!     new value of ZDOTA(NACT) and branch if it is not acceptable.
!
   Do 160 k = 1, nact
160 vmultc (k) = max (0.0_wp, vmultc(k)-ratio*vmultd(k))
   If (icon < nact) Then
      isave = iact (icon)
      vsave = vmultc (icon)
      k = icon
170   kp = k + 1
      kw = iact (kp)
      sp = 0.0_wp
      Do 180 i = 1, n
180   sp = sp + z (i, k) * a (i, kw)
      temp = Sqrt (sp*sp+zdota(kp)**2)
      alpha = zdota (kp) / temp
      beta = sp / temp
      zdota (kp) = alpha * zdota (k)
      zdota (k) = temp
      Do 190 i = 1, n
         temp = alpha * z (i, kp) + beta * z (i, k)
         z (i, kp) = alpha * z (i, k) - beta * z (i, kp)
190   z (i, k) = temp
      iact (k) = kw
      vmultc (k) = vmultc (kp)
      k = kp
      If (k < nact) Go To 170
      iact (k) = isave
      vmultc (k) = vsave
   End If
   temp = 0.0_wp
   Do 200 i = 1, n
200 temp = temp + z (i, nact) * a (i, kk)
   If (temp == 0.0_wp) Go To 490
   zdota (nact) = temp
   vmultc (icon) = 0.0_wp
   vmultc (nact) = ratio
!
!     Update IACT and ensure that the objective function continues to be
!     treated as the last active constraint when MCON>M.
!
210 iact (icon) = iact (nact)
   iact (nact) = kk
   If (mcon > m .And. kk /= mcon) Then
      k = nact - 1
      sp = 0.0_wp
      Do 220 i = 1, n
220   sp = sp + z (i, k) * a (i, kk)
      temp = Sqrt (sp*sp+zdota(nact)**2)
      alpha = zdota (nact) / temp
      beta = sp / temp
      zdota (nact) = alpha * zdota (k)
      zdota (k) = temp
      Do 230 i = 1, n
         temp = alpha * z (i, nact) + beta * z (i, k)
         z (i, nact) = alpha * z (i, k) - beta * z (i, nact)
230   z (i, k) = temp
      iact (nact) = iact (k)
      iact (k) = kk
      temp = vmultc (k)
      vmultc (k) = vmultc (nact)
      vmultc (nact) = temp
   End If
!
!     If stage one is in progress, then set SDIRN to the direction of the next
!     change to the current vector of variables.
!
   If (mcon > m) Go To 320
   kk = iact (nact)
   temp = 0.0_wp
   Do 240 i = 1, n
240 temp = temp + sdirn (i) * a (i, kk)
   temp = temp - 1.0_wp
   temp = temp / zdota (nact)
   Do 250 i = 1, n
250 sdirn (i) = sdirn (i) - temp * z (i, nact)
   Go To 340
!
!     Delete the constraint that has the index IACT(ICON) from the active set.
!
260 If (icon < nact) Then
      isave = iact (icon)
      vsave = vmultc (icon)
      k = icon
270   kp = k + 1
      kk = iact (kp)
      sp = 0.0_wp
      Do 280 i = 1, n
280   sp = sp + z (i, k) * a (i, kk)
      temp = Sqrt (sp*sp+zdota(kp)**2)
      alpha = zdota (kp) / temp
      beta = sp / temp
      zdota (kp) = alpha * zdota (k)
      zdota (k) = temp
      Do 290 i = 1, n
         temp = alpha * z (i, kp) + beta * z (i, k)
         z (i, kp) = alpha * z (i, k) - beta * z (i, kp)
290   z (i, k) = temp
      iact (k) = kk
      vmultc (k) = vmultc (kp)
      k = kp
      If (k < nact) Go To 270
      iact (k) = isave
      vmultc (k) = vsave
   End If
   nact = nact - 1
!
!     If stage one is in progress, then set SDIRN to the direction of the next
!     change to the current vector of variables.
!
   If (mcon > m) Go To 320
   temp = 0.0_wp
   Do 300 i = 1, n
300 temp = temp + sdirn (i) * z (i, nact+1)
   Do 310 i = 1, n
310 sdirn (i) = sdirn (i) - temp * z (i, nact+1)
   Go To 340
!
!     Pick the next search direction of stage two.
!
320 temp = 1.0_wp / zdota (nact)
   Do 330 i = 1, n
330 sdirn (i) = temp * z (i, nact)
!
!     Calculate the step to the boundary of the trust region or take the step
!     that reduces RESMAX to zero. The two statements below that include the
!     factor 1.0E-6 prevent some harmless underflows that occurred in a test
!     calculation. Further, we skip the step if it could be zero within a
!     reasonable tolerance for computer rounding errors.
!
340 dd = rho * rho
   sd = 0.0_wp
   ss = 0.0_wp
   Do 350 i = 1, n
      If (Abs(dx(i)) >= 1.0e-6_wp*rho) dd = dd - dx (i) ** 2
      sd = sd + dx (i) * sdirn (i)
350 ss = ss + sdirn (i) ** 2
   If (dd <= 0.0_wp) Go To 490
   temp = Sqrt (ss*dd)
   If (Abs(sd) >= 1.0e-6_wp*temp) temp = Sqrt (ss*dd+sd*sd)
   stpful = dd / (temp+sd)
   step = stpful
   If (mcon == m) Then
      acca = step + 0.1_wp * resmax
      accb = step + 0.2_wp * resmax
      If (step >= acca .Or. acca >= accb) Go To 480
      step = min (step, resmax)
   End If
!
!     Set DXNEW to the new variables if STEP is the steplength, and reduce
!     RESMAX to the corresponding maximum residual if stage one is being done.
!     Because DXNEW will be changed during the calculation of some Lagrange
!     multipliers, it will be restored to the following value later.
!
   Do 360 i = 1, n
360 dxnew (i) = dx (i) + step * sdirn (i)
   If (mcon == m) Then
      resold = resmax
      resmax = 0.0_wp
      Do 380 k = 1, nact
         kk = iact (k)
         temp = b (kk)
         Do 370 i = 1, n
370      temp = temp - a (i, kk) * dxnew (i)
         resmax = max (resmax, temp)
380   Continue
   End If
!
!     Set VMULTD to the VMULTC vector that would occur if DX became DXNEW. A
!     device is included to force VMULTD(K)=0.0 if deviations from this value
!     can be attributed to computer rounding errors. First calculate the new
!     Lagrange multipliers.
!
   k = nact
390 zdotw = 0.0_wp
   zdwabs = 0.0_wp
   Do 400 i = 1, n
      temp = z (i, k) * dxnew (i)
      zdotw = zdotw + temp
400 zdwabs = zdwabs + Abs (temp)
   acca = zdwabs + 0.1_wp * Abs (zdotw)
   accb = zdwabs + 0.2_wp * Abs (zdotw)
   If (zdwabs >= acca .Or. acca >= accb) zdotw = 0.0_wp
   vmultd (k) = zdotw / zdota (k)
   If (k >= 2) Then
      kk = iact (k)
      Do 410 i = 1, n
410   dxnew (i) = dxnew (i) - vmultd (k) * a (i, kk)
      k = k - 1
      Go To 390
   End If
   If (mcon > m) vmultd (nact) = max (0.0_wp, vmultd(nact))
!
!     Complete VMULTC by finding the new constraint residuals.
!
   Do 420 i = 1, n
420 dxnew (i) = dx (i) + step * sdirn (i)
   If (mcon > nact) Then
      kl = nact + 1
      Do 440 k = kl, mcon
         kk = iact (k)
         sum = resmax - b (kk)
         sumabs = resmax + Abs (b(kk))
         Do 430 i = 1, n
            temp = a (i, kk) * dxnew (i)
            sum = sum + temp
430      sumabs = sumabs + Abs (temp)
         acca = sumabs + 0.1_wp * Abs (sum)
         accb = sumabs + 0.2_wp * Abs (sum)
         If (sumabs >= acca .Or. acca >= accb) sum = 0.0_wp
440   vmultd (k) = sum
   End If
!
!     Calculate the fraction of the step from DX to DXNEW that will be taken.
!
   ratio = 1.0_wp
   icon = 0
   Do 450 k = 1, mcon
      If (vmultd(k) < 0.0_wp) Then
         temp = vmultc (k) / (vmultc(k)-vmultd(k))
         If (temp < ratio) Then
            ratio = temp
            icon = k
         End If
      End If
450 Continue
!
!     Update DX, VMULTC and RESMAX.
!
   temp = 1.0_wp - ratio
   Do 460 i = 1, n
460 dx (i) = temp * dx (i) + ratio * dxnew (i)
   Do 470 k = 1, mcon
470 vmultc (k) = max (0.0_wp, temp*vmultc(k)+ratio*vmultd(k))
   If (mcon == m) resmax = resold + ratio * (resmax-resold)
!
!     If the full step is not acceptable then begin another iteration.
!     Otherwise switch to stage two or end the calculation.
!
   If (icon > 0) Go To 70
   If (step == stpful) Go To 500
480 mcon = m + 1
   icon = mcon
   iact (mcon) = mcon
   vmultc (mcon) = 0.0_wp
   Go To 60
!
!     We employ any freedom that may be available to reduce the objective
!     function before returning a DX whose length is less than RHO.
!
490 If (mcon == m) Go To 480
   ifull = 0
500 Return
End Subroutine trstlp

Subroutine cobyla_test()
!------------------------------------------------------------------------------
!     Main program of test problems in Report DAMTP 1992/NA5.
!------------------------------------------------------------------------------
    !Common nprob
    Dimension x (10), xopt (10), w (3000), iact (51)
    Do 180 nprob = 1, 10
       If (nprob == 1) Then
    !
    !     Minimization of a simple quadratic function of two variables.
    !
          Print 10
    10    Format (/ 7 x, 'Output from test problem 1 (Simple quadratic)')
          n = 2
          m = 0
          xopt (1) = - 1.0_wp
          xopt (2) = 0.0_wp
       Else If (nprob == 2) Then
    !
    !     Easy two dimensional minimization in unit circle.
    !
          Print 20
    20    Format (/ 7 x, 'Output from test problem 2 (2D unit circle ',&
          'calculation)')
          n = 2
          m = 1
          xopt (1) = Sqrt (0.5_wp)
          xopt (2) = - xopt (1)
       Else If (nprob == 3) Then
    !
    !     Easy three dimensional minimization in ellipsoid.
    !
          Print 30
    30    Format (/ 7 x, 'Output from test problem 3 (3D ellipsoid ',&
         'calculation)')
          n = 3
          m = 1
          xopt (1) = 1.0_wp / Sqrt (3.0_wp)
          xopt (2) = 1.0_wp / Sqrt (6.0_wp)
          xopt (3) = - 1.0_wp / 3.0_wp
       Else If (nprob == 4) Then
    !
    !     Weak version of Rosenbrock's problem.
    !
          Print 40
    40    Format (/ 7 x, 'Output from test problem 4 (Weak Rosenbrock)')
          n = 2
          m = 0
          xopt (1) = - 1.0_wp
          xopt (2) = 1.0_wp
       Else If (nprob == 5) Then
    !
    !     Intermediate version of Rosenbrock's problem.
    !
          Print 50
    50    Format (/ 7 x, 'Output from test problem 5 (Intermediate ',&
         'Rosenbrock)')
          n = 2
          m = 0
          xopt (1) = - 1.0_wp
          xopt (2) = 1.0_wp
       Else If (nprob == 6) Then
    !
    !     This problem is taken from Fletcher's book Practical Methods of
    !     Optimization and has the equation number (9.1.15).
    !
          Print 60
    60    Format (/ 7 x, 'Output from test problem 6 (Equation ',&
          '(9.1.15) in Fletcher)')
          n = 2
          m = 2
          xopt (1) = Sqrt (0.5_wp)
          xopt (2) = xopt (1)
       Else If (nprob == 7) Then
    !
    !     This problem is taken from Fletcher's book Practical Methods of
    !     Optimization and has the equation number (14.4.2).
    !
          Print 70
    70    Format (/ 7 x, 'Output from test problem 7 (Equation ',&
         '(14.4.2) in Fletcher)')
          n = 3
          m = 3
          xopt (1) = 0.0_wp
          xopt (2) = - 3.0_wp
          xopt (3) = - 3.0_wp
       Else If (nprob == 8) Then
    !
    !     This problem is taken from page 66 of Hock and Schittkowski's book Test
    !     Examples for Nonlinear Programming Codes. It is their test problem Number
    !     43, and has the name Rosen-Suzuki.
    !
          Print 80
    80    Format (/ 7 x, 'Output from test problem 8 (Rosen-Suzuki)')
          n = 4
          m = 3
          xopt (1) = 0.0_wp
          xopt (2) = 1.0_wp
          xopt (3) = 2.0_wp
          xopt (4) = - 1.0_wp
       Else If (nprob == 9) Then
    !
    !     This problem is taken from page 111 of Hock and Schittkowski's
    !     book Test Examples for Nonlinear Programming Codes. It is their
    !     test problem Number 100.
    !
          Print 90
    90    Format (/ 7 x, 'Output from test problem 9 (Hock and ',&
         'Schittkowski 100)')
          n = 7
          m = 4
          xopt (1) = 2.330499_wp
          xopt (2) = 1.951372_wp
          xopt (3) = - 0.4775414_wp
          xopt (4) = 4.365726_wp
          xopt (5) = - 0.624487_wp
          xopt (6) = 1.038131_wp
          xopt (7) = 1.594227_wp
       Else If (nprob == 10) Then
    !
    !     This problem is taken from page 415 of Luenberger's book Applied
    !     Nonlinear Programming. It is to maximize the area of a hexagon of
    !     unit diameter.
    !
          Print 100
    100   Format (/ 7 x, 'Output from test problem 10 (Hexagon area)')
          n = 9
          m = 14
       End If
       Do 160 icase = 1, 2
          Do 120 i = 1, n
    120   x (i) = 1.0_wp
          rhobeg = 0.5_wp
          rhoend = 0.001_wp
          If (icase == 2) rhoend = 0.0001_wp
          iprint = 1
          maxfun = 2000
          Call cobyla (n, m, x, rhobeg, rhoend, iprint, maxfun, w, iact, calcfc)
          If (nprob == 10) Then
             tempa = x (1) + x (3) + x (5) + x (7)
             tempb = x (2) + x (4) + x (6) + x (8)
             tempc = 0.5_wp / Sqrt (tempa*tempa+tempb*tempb)
             tempd = tempc * Sqrt (3.0_wp)
             xopt (1) = tempd * tempa + tempc * tempb
             xopt (2) = tempd * tempb - tempc * tempa
             xopt (3) = tempd * tempa - tempc * tempb
             xopt (4) = tempd * tempb + tempc * tempa
             Do 130 i = 1, 4
    130      xopt (i+4) = xopt (i)
          End If
          temp = 0.0_wp
          Do 140 i = 1, n
    140   temp = temp + (x(i)-xopt(i)) ** 2
          Print 150, Sqrt (temp)
    150   Format (/ 5 x, 'Least squares error in variables =', 1 pe16.6)
    160 Continue
       Print 170
    170 Format (2 x, '----------------------------------------------',&
        '--------------------')
    180 Continue

contains

    Subroutine calcfc (n, m, x, f, con)
          !Common nprob
          Dimension x (*), con (*)
          If (nprob == 1) Then
    !
    !     Test problem 1 (Simple quadratic)
    !
             f = 10.0_wp * (x(1)+1.0_wp) ** 2 + x (2) ** 2
          Else If (nprob == 2) Then
    !
    !    Test problem 2 (2D unit circle calculation)
    !
             f = x (1) * x (2)
             con (1) = 1.0_wp - x (1) ** 2 - x (2) ** 2
          Else If (nprob == 3) Then
    !
    !     Test problem 3 (3D ellipsoid calculation)
    !
             f = x (1) * x (2) * x (3)
             con (1) = 1.0_wp - x (1) ** 2 - 2.0_wp * x (2) ** 2 - 3.0_wp * x (3) **2
          Else If (nprob == 4) Then
    !
    !     Test problem 4 (Weak Rosenbrock)
    !
             f = (x(1)**2-x(2)) ** 2 + (1.0_wp+x(1)) ** 2
          Else If (nprob == 5) Then
    !
    !     Test problem 5 (Intermediate Rosenbrock)
    !
             f = 10.0_wp * (x(1)**2-x(2)) ** 2 + (1.0_wp+x(1)) ** 2
          Else If (nprob == 6) Then
    !
    !     Test problem 6 (Equation (9.1.15) in Fletcher's book)
    !
             f = - x (1) - x (2)
             con (1) = x (2) - x (1) ** 2
             con (2) = 1.0_wp - x (1) ** 2 - x (2) ** 2
          Else If (nprob == 7) Then
    !
    !     Test problem 7 (Equation (14.4.2) in Fletcher's book)
    !
             f = x (3)
             con (1) = 5.0_wp * x (1) - x (2) + x (3)
             con (2) = x (3) - x (1) ** 2 - x (2) ** 2 - 4.0_wp * x (2)
             con (3) = x (3) - 5.0_wp * x (1) - x (2)
          Else If (nprob == 8) Then
    !
    !     Test problem 8 (Rosen-Suzuki)
    !
             f = x (1) ** 2 + x (2) ** 2 + 2.0_wp * x (3) ** 2 + x (4) ** 2 - &
              5.0_wp * x (1) - 5.0_wp * x (2) - 21.0_wp * x (3) + 7.0_wp * x (4)
             con (1) = 8.0_wp - x (1) ** 2 - x (2) ** 2 - x (3) ** 2 - x (4) &
              ** 2 - x (1) + x (2) - x (3) + x (4)
             con (2) = 10.0_wp - x (1) ** 2 - 2.0_wp * x (2) ** 2 - x (3) ** 2 - &
              2.0_wp * x (4) ** 2 + x (1) + x (4)
             con (3) = 5.0_wp - 2.0_wp * x (1) ** 2 - x (2) ** 2 - x (3) ** 2 - &
              2.0_wp * x (1) + x (2) + x (4)
          Else If (nprob == 9) Then
    !
    !     Test problem 9 (Hock and Schittkowski 100)
    !
             f = (x(1)-10.0_wp) ** 2 + 5.0_wp * (x(2)-12.0_wp) ** 2 + x (3) ** 4 + &
              3.0_wp * (x(4)-11.0_wp) ** 2 + 10.0_wp * x (5) ** 6 + 7.0_wp * x (6) ** 2 &
              + x (7) ** 4 - 4.0_wp * x (6) * x (7) - 10.0_wp * x (6) - 8.0_wp * x &
              (7)
             con (1) = 127.0_wp - 2.0_wp * x (1) ** 2 - 3.0_wp * x (2) ** 4 - x (3) &
              - 4.0_wp * x (4) ** 2 - 5.0_wp * x (5)
             con (2) = 282.0_wp - 7.0_wp * x (1) - 3.0_wp * x (2) - 10.0_wp * x (3) ** &
              2 - x (4) + x (5)
             con (3) = 196.0_wp - 23.0_wp * x (1) - x (2) ** 2 - 6.0_wp * x (6) ** 2 &
              + 8.0_wp * x (7)
             con (4) = - 4.0_wp * x (1) ** 2 - x (2) ** 2 + 3.0_wp * x (1) * x &
              (2) - 2.0_wp * x (3) ** 2 - 5.0_wp * x (6) + 11.0_wp * x (7)
          Else If (nprob == 10) Then
    !
    !     Test problem 10 (Hexagon area)
    !
             f = - 0.5_wp * &
              (x(1)*x(4)-x(2)*x(3)+x(3)*x(9)-x(5)*x(9)+x(5)*x(8)-x(6)*x(7))
             con (1) = 1.0_wp - x (3) ** 2 - x (4) ** 2
             con (2) = 1.0_wp - x (9) ** 2
             con (3) = 1.0_wp - x (5) ** 2 - x (6) ** 2
             con (4) = 1.0_wp - x (1) ** 2 - (x(2)-x(9)) ** 2
             con (5) = 1.0_wp - (x(1)-x(5)) ** 2 - (x(2)-x(6)) ** 2
             con (6) = 1.0_wp - (x(1)-x(7)) ** 2 - (x(2)-x(8)) ** 2
             con (7) = 1.0_wp - (x(3)-x(5)) ** 2 - (x(4)-x(6)) ** 2
             con (8) = 1.0_wp - (x(3)-x(7)) ** 2 - (x(4)-x(8)) ** 2
             con (9) = 1.0_wp - x (7) ** 2 - (x(8)-x(9)) ** 2
             con (10) = x (1) * x (4) - x (2) * x (3)
             con (11) = x (3) * x (9)
             con (12) = - x (5) * x (9)
             con (13) = x (5) * x (8) - x (6) * x (7)
             con (14) = x (9)
          End If
    End Subroutine calcfc

End Subroutine cobyla_test

end module cobyla_module