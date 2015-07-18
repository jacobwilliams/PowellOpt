subroutine addcon (n,m,a,ia,iact,nact,z,u,relacc,indxbd,ztc, &
  cgrad)
implicit real*8 (a-h,o-z)
dimension a(ia,*),iact(*),z(*),u(*),ztc(*),cgrad(*)
np=nact+1
icon=iact(indxbd)
iact(indxbd)=iact(np)
iact(np)=icon
!
!     Form ZTC when the new constraint is a bound.
!
if (icon .gt. m) then
    inewbd=icon-m
    if (inewbd .le. n) then
        temp=-1.0
    else
        inewbd=inewbd-n
        temp=1.0
    end if
    iznbd=inewbd*n-n
    do 10 j=1,n
10     ztc(j)=temp*z(iznbd+j)
!
!     Else form ZTC for an ordinary constraint.
!
else
    do 20 i=1,n
20     cgrad(i)=a(i,icon)
    do 30 j=1,n
    ztc(j)=0.0
    iz=j
    do 30 i=1,n
    ztc(j)=ztc(j)+z(iz)*cgrad(i)
30     iz=iz+n
end if
!
!     Find any Givens rotations to apply to the last columns of Z.
!
j=n
40 jp=j
j=j-1
if (j .gt. nact) then
    if (ztc(jp) .eq. 0.0) goto 40
    if (dabs(ztc(jp)) .le. relacc*dabs(ztc(j))) then
        temp=dabs(ztc(j))
    else if (dabs(ztc(j)) .le. relacc*dabs(ztc(jp))) then
        temp=dabs(ztc(jp))
    else
        temp=dabs(ztc(jp))*dsqrt(1.0+(ztc(j)/ztc(jp))**2)
    end if
    wcos=ztc(j)/temp
    wsin=ztc(jp)/temp
    ztc(j)=temp
!
!     Apply the rotation when the new constraint is a bound.
!
    iz=j
    if (icon .gt. m) then
        do 50 i=1,n
        temp=wcos*z(iz+1)-wsin*z(iz)
        z(iz)=wcos*z(iz)+wsin*z(iz+1)
        z(iz+1)=temp
50         iz=iz+n
        z(iznbd+jp)=0.0
!
!     Else apply the rotation for an ordinary constraint.
!
    else
        wpiv=0.0
        do 60 i=1,n
        tempa=wcos*z(iz+1)
        tempb=wsin*z(iz)
        temp=dabs(cgrad(i))*(dabs(tempa)+dabs(tempb))
        if (temp .gt. wpiv) then
            wpiv=temp
            ipiv=i
        end if
        z(iz)=wcos*z(iz)+wsin*z(iz+1)
        z(iz+1)=tempa-tempb
60         iz=iz+n
!
!     Ensure orthogonality of Z(.,JP) to CGRAD.
!
        sum=0.0
        iz=jp
        do 70 i=1,n
        sum=sum+z(iz)*cgrad(i)
70         iz=iz+n
        if (sum .ne. 0.0) then
            iz=ipiv*n-n+jp
            z(iz)=z(iz)-sum/cgrad(ipiv)
        end if
    end if
    go to 40
end if
!
!     Test for linear independence in the proposed new active set.
!
if (ztc(np) .eq. 0.0) goto 90
if (icon .le. m) then
    sum=0.0
    sumabs=0.0
    iz=np
    do 80 i=1,n
    temp=z(iz)*cgrad(i)
    sum=sum+temp
    sumabs=sumabs+dabs(temp)
80     iz=iz+n
    if (dabs(sum) .le. relacc*sumabs) goto 90
end if
!
!     Set the new diagonal element of U and return.
!
u(np)=1.0/ztc(np)
nact=np
90 return
end
subroutine adjtol (n,m,a,ia,b,xl,xu,x,iact,nact,xbig,relacc,tol, &
  meql)
implicit real*8 (a-h,o-z)
dimension a(ia,*),b(*),xl(*),xu(*),x(*),iact(*),xbig(*)
!
!     Set VIOL to the greatest relative constraint residual of the first
!       NACT constraints.
!
viol=0.0
if (nact .gt. meql) then
    kl=meql+1
    do 20 k=kl,nact
    j=iact(k)
    if (j .le. m) then
        res=b(j)
        resabs=dabs(b(j))
        do 10 i=1,n
        res=res-a(i,j)*x(i)
10         resabs=resabs+dabs(a(i,j)*xbig(i))
    else
        jm=j-m
        if (jm .le. n) then
            res=x(jm)-xl(jm)
            resabs=xbig(jm)+dabs(xl(jm))
        else
            jm=jm-n
            res=xu(jm)-x(jm)
            resabs=xbig(jm)+dabs(xu(jm))
        end if
    end if
    if (res .gt. 0.0) viol=dmax1(viol,res/resabs)
20     continue
end if
!
!     Adjust TOL.
!
tol=0.1*dmin1(tol,viol)
if (tol .le. relacc+relacc) then
    tol=relacc
    do 30 i=1,n
30     xbig(i)=dabs(x(i))
end if
return
end
subroutine conres (n,m,a,ia,b,xl,xu,x,iact,nact,par,g,z,u,xbig, &
  bres,d,ztg,relacc,tol,stepcb,sumres,meql,msat,mtot,indxbd, &
  gm,gmnew,parnew,cgrad)
implicit real*8 (a-h,o-z)
dimension a(ia,*),b(*),xl(*),xu(*),x(*),iact(*),par(*),g(*), &
  z(*),u(*),xbig(*),bres(*),d(*),ztg(*),gm(*),gmnew(*),parnew(*), &
  cgrad(*)
idiff=mtot-msat
!
!     Calculate and partition the residuals of the inactive constraints,
!       and set the gradient vector when seeking feasibility.
!
if (idiff .gt. 0.0) then
    do 10 i=1,n
10     g(i)=0.0
    sumres=0.0
end if
msatk=msat
mdeg=nact
msat=nact
kl=meql+1
do 50 k=kl,mtot
j=iact(k)
!
!     Calculate the residual of the current constraint.
!
if (j .le. m) then
    res=b(j)
    resabs=dabs(b(j))
    do 20 i=1,n
    res=res-x(i)*a(i,j)
20     resabs=resabs+dabs(xbig(i)*a(i,j))
else
    jm=j-m
    if (jm .le. n) then
        res=x(jm)-xl(jm)
        resabs=dabs(xbig(jm))+dabs(xl(jm))
    else
        jm=jm-n
        res=xu(jm)-x(jm)
        resabs=dabs(xbig(jm))+dabs(xu(jm))
    end if
end if
bres(j)=res
!
!     Set TEMP to the relative residual.
!
temp=0.0
if (resabs .ne. 0.0) temp=res/resabs
if (k .gt. msatk .and. temp .lt. 0.0) then
    if (temp+relacc .ge. 0.0) then
        if (j .le. m) then
            sum=dabs(b(j))
            do 30 i=1,n
30             sum=sum+dabs(x(i)*a(i,j))
        else
            jm=j-m
            if (jm .le. n) then
                sum=dabs(x(jm))+dabs(xl(jm))
            else
                sum=dabs(x(jm-n))+dabs(xu(jm-n))
            end if
        end if
        if (dabs(res) .le. sum*relacc) temp=0.0
    end if
end if
!
!     Place the residual in the appropriate position.
!
if (k .le. nact) goto 50
if (k .le. msatk .or. temp .ge. 0.0) then
    msat=msat+1
    if (msat .lt. k) then
        iact(k)=iact(msat)
    end if
    if (temp .gt. tol) then
        iact(msat)=j
    else
        mdeg=mdeg+1
        iact(msat)=iact(mdeg)
        iact(mdeg)=j
    end if
!
!     Update the gradient and SUMRES if the constraint is violated when
!       seeking feasibility.
!
else
    if (j .le. m) then
        do 40 i=1,n
40         g(i)=g(i)+a(i,j)
    else
        j=j-m
        if (j .le. n) then
            g(j)=g(j)-1.0
        else
            g(j-n)=g(j-n)+1.0
        end if
    end if
    sumres=sumres+dabs(res)
end if
50 continue
!
!     Seek the next search direction unless CONRES was called from GETFES
!       and feasibility has been achieved.
!
stepcb=0.0
if (idiff .gt. 0 .and. msat .eq. mtot) goto 60
call getd (n,m,a,ia,iact,nact,par,g,z,u,d,ztg,relacc,ddotg,meql, &
  mdeg,gm,gmnew,parnew,cgrad)
!
!     Calculate the (bound on the) step-length due to the constraints.
!
if (ddotg .lt. 0.0) then
    call stepbd (n,m,a,ia,iact,bres,d,stepcb,ddotg,mdeg,msat, &
      mtot,indxbd)
end if
if (idiff .eq. 0) sumres=ddotg
60 return
end
subroutine delcon (n,m,a,ia,iact,nact,z,u,relacc,idrop)
implicit real*8 (a-h,o-z)
dimension a(ia,*),iact(*),z(*),u(*)
nm=nact-1
if (idrop .eq. nact) goto 60
isave=iact(idrop)
!
!     Cycle through the constraint exchanges that are needed.
!
do 50 j=idrop,nm
jp=j+1
icon=iact(jp)
iact(j)=icon
!
!     Calculate the (J,JP) element of R.
!
if (icon .le. m) then
    rjjp=0.0
    iz=j
    do 10 i=1,n
    rjjp=rjjp+z(iz)*a(i,icon)
10     iz=iz+n
else
    ibd=icon-m
    if (ibd .le. n) then
        izbd=ibd*n-n
        rjjp=-z(izbd+j)
    else
        ibd=ibd-n
        izbd=ibd*n-n
        rjjp=z(izbd+j)
    end if
end if
!
!     Calculate the parameters of the next rotation.
!
ujp=u(jp)
temp=rjjp*ujp
denom=dabs(temp)
if (denom*relacc .lt. 1.0) denom=dsqrt(1.0+denom*denom)
wcos=temp/denom
wsin=1.0/denom
!
!     Rotate Z when a bound constraint is promoted.
!
iz=j
if (icon .gt. m) then
    do 20 i=1,n
    temp=wcos*z(iz+1)-wsin*z(iz)
    z(iz)=wcos*z(iz)+wsin*z(iz+1)
    z(iz+1)=temp
20     iz=iz+n
    z(izbd+jp)=0.0
!
!     Rotate Z when an ordinary constraint is promoted.
!
else
    wpiv=0.0
    do 30 i=1,n
    tempa=wcos*z(iz+1)
    tempb=wsin*z(iz)
    temp=dabs(a(i,icon))*(dabs(tempa)+dabs(tempb))
    if (temp .gt. wpiv) then
        wpiv=temp
        ipiv=i
    end if
    z(iz)=wcos*z(iz)+wsin*z(iz+1)
    z(iz+1)=tempa-tempb
30     iz=iz+n
!
!     Ensure orthogonality to promoted constraint.
!
    sum=0.0
    iz=jp
    do 40 i=1,n
    sum=sum+z(iz)*a(i,icon)
40     iz=iz+n
    if (sum .ne. 0.0) then
        iz=ipiv*n-n+jp
        z(iz)=z(iz)-sum/a(ipiv,icon)
    end if
end if
!
!     Set the new diagonal elements of U.
!
u(jp)=-denom*u(j)
u(j)=ujp/denom
50 continue
!
!     Return.
!
iact(nact)=isave
60 nact=nm
return
end
subroutine eqcons (n,m,meq,a,ia,b,xu,iact,meql,info,z,u,relacc, &
  am,cgrad)
implicit real*8 (a-h,o-z)
dimension a(ia,*),b(*),xu(*),iact(*),z(*),u(*),am(*),cgrad(*)
!
!     Try to add the next equality constraint to the active set.
!
do 50 keq=1,meq
if (meql .lt. n) then
    np=meql+1
    iact(np)=keq
    call addcon (n,m,a,ia,iact,meql,z,u,relacc,np,am,cgrad)
    if (meql .eq. np) goto 50
end if
!
!     If linear dependence occurs then find the multipliers of the
!       dependence relation and apply them to the right hand sides.
!
sum=b(keq)
sumabs=dabs(b(keq))
if (meql .gt. 0) then
    do 10 i=1,n
10     am(i)=a(i,keq)
    k=meql
20     vmult=0.0
    iz=k
    do 30 i=1,n
    vmult=vmult+z(iz)*am(i)
30     iz=iz+n
    vmult=vmult*u(k)
    j=iact(k)
    if (j .le. m) then
        do 40 i=1,n
40         am(i)=am(i)-vmult*a(i,j)
        rhs=b(j)
    else
        jm=j-m-n
        am(jm)=am(jm)-vmult
        rhs=xu(jm)
    end if
    sum=sum-rhs*vmult
    sumabs=sumabs+dabs(rhs*vmult)
    k=k-1
    if (k .ge. 1) goto 20
end if
!
!     Error return if the constraints are inconsistent.
!
if (dabs(sum) .gt. relacc*sumabs) then
    info=5
    goto 60
end if
50 continue
60 return
end
subroutine fgcalc (n,x,f,g)
implicit real*8 (a-h,o-z)
dimension x(*),g(*)
!
!     Calculate the objective function and its gradient.
!
wa=(x(1)-x(3))**2+(x(2)-x(4))**2
wb=(x(3)-x(5))**2+(x(4)-x(6))**2
wc=(x(5)-x(1))**2+(x(6)-x(2))**2
f=1.0/(wa**8)+1.0/(wb**8)+1.0/(wc**8)
g(1)=16.0*((x(3)-x(1))/(wa**9)+(x(5)-x(1))/(wc**9))
g(2)=16.0*((x(4)-x(2))/(wa**9)+(x(6)-x(2))/(wc**9))
g(3)=16.0*((x(5)-x(3))/(wb**9)+(x(1)-x(3))/(wa**9))
g(4)=16.0*((x(6)-x(4))/(wb**9)+(x(2)-x(4))/(wa**9))
g(5)=16.0*((x(1)-x(5))/(wc**9)+(x(3)-x(5))/(wb**9))
g(6)=16.0*((x(2)-x(6))/(wc**9)+(x(4)-x(6))/(wb**9))
return
end

subroutine getd (n,m,a,ia,iact,nact,par,g,z,u,d,ztg,relacc, &
  ddotg,meql,mdeg,gm,gmnew,parnew,cgrad)
implicit real*8 (a-h,o-z)
dimension a(ia,*),iact(*),par(*),g(*),z(*),u(*),d(*),ztg(*), &
  gm(*),gmnew(*),parnew(*),cgrad(*)
!
!     Initialize GM and cycle backwards through the active set.
!
10 do 20 i=1,n
20 gm(i)=g(i)
k=nact
30 if (k .gt. 0) then
!
!     Set TEMP to the next multiplier, but reduce the active set if
!       TEMP has an unacceptable sign.
!
    temp=0.0
    iz=k
    do 40 i=1,n
    temp=temp+z(iz)*gm(i)
40     iz=iz+n
    temp=temp*u(k)
    if (k .gt. meql .and. temp .gt. 0.0) then
        call delcon (n,m,a,ia,iact,nact,z,u,relacc,k)
        goto 10
    end if
!
!     Update GM using the multiplier that has just been calculated.
!
    j=iact(k)
    if (j .le. m) then
        do 50 i=1,n
50         gm(i)=gm(i)-temp*a(i,j)
    else
        jm=j-m
        if (jm .le. n) then
            gm(jm)=gm(jm)+temp
        else
            gm(jm-n)=gm(jm-n)-temp
        end if
    end if
    par(k)=temp
    k=k-1
    goto 30
end if
!
!     Calculate the search direction and DDOTG.
!
ddotg=0.0
if (nact .lt. n) then
    call sdegen (n,m,a,ia,iact,nact,par,z,u,d,ztg,gm,relacc, &
      ddotgm,meql,mdeg,gmnew,parnew,cgrad)
    if (ddotgm .lt. 0.0) then
        do 60 i=1,n
60         ddotg=ddotg+d(i)*g(i)
    end if
end if
return
end
subroutine getfes (n,m,a,ia,b,xl,xu,x,iact,nact,par,info,g,z, &
  u,xbig,relacc,tol,meql,msat,mtot,bres,d,ztg,gm,gmnew,parnew, &
  cgrad)
implicit real*8 (a-h,o-z)
dimension a(ia,*),b(*),xl(*),xu(*),x(*),iact(*),par(*),g(*),z(*), &
  u(*),xbig(*),bres(*),d(*),ztg(*),gm(*),gmnew(*),parnew(*), &
  cgrad(*)
!
!     Make the correction to X for the active constraints.
!
info=0
10 call satact (n,m,a,ia,b,xl,xu,x,iact,nact,info,z,u,xbig,relacc, &
  tol,meql)
if (info .gt. 0) msat=nact
if (msat .eq. mtot) goto 60
!
!     Try to correct the infeasibility.
!
20 msatk=msat
sumrsk=0.0
30 call conres (n,m,a,ia,b,xl,xu,x,iact,nact,par,g,z,u,xbig,bres, &
  d,ztg,relacc,tol,stepcb,sumres,meql,msat,mtot,indxbd,gm,gmnew, &
  parnew,cgrad)
!
!     Include the new constraint in the active set.
!
if (stepcb .gt. 0.0) then
    do 40 i=1,n
    x(i)=x(i)+stepcb*d(i)
40     xbig(i)=dmax1(xbig(i),dabs(x(i)))
    call addcon (n,m,a,ia,iact,nact,z,u,relacc,indxbd,gmnew,cgrad)
end if
!
!     Test whether to continue the search for feasibility.
!
if (msat .lt. mtot) then
    if (stepcb .eq. 0.0) goto 50
    if (msatk .lt. msat) goto 20
    if (sumrsk .eq. 0.0 .or. sumres .lt. sumrsk) then
        sumrsk=sumres
        itest=0
    end if
    itest=itest+1
    if (itest .le. 2) goto 30
!
!     Reduce TOL if it may be too large to allow feasibility.
!
50     if (tol .gt. relacc) then
        call adjtol (n,m,a,ia,b,xl,xu,x,iact,nact,xbig,relacc, &
          tol,meql)
        goto 10
    end if
end if
60 return
end
subroutine getmin (n,m,meq,a,ia,b,xl,xu,x,acc,iact,nact,par, &
  iprint,info,w)
implicit real*8 (a-h,o-z)
dimension a(ia,*),b(*),xl(*),xu(*),x(*),iact(*),par(*),w(*)
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
ig=1
ireskt=ig+n
iz=ireskt+n
iu=iz+n*n
ixbig=iu+n
ibres=ixbig+n
id=ibres+m+n+n
iztg=id+n
igm=iztg+n
ixs=igm+n
igs=ixs+n
!
!     Call the optimization package.
!
call minflc (n,m,meq,a,ia,b,xl,xu,x,acc,iact,nact,par,iprint, &
  info,w(ig),w(iz),w(iu),w(ixbig),w(ireskt),w(ibres),w(id), &
  w(iztg),w(igm),w(ixs),w(igs))
return
end
subroutine initzu (n,m,xl,xu,x,iact,meql,info,z,u,xbig,relacc)
implicit real*8 (a-h,o-z)
dimension xl(*),xu(*),x(*),iact(*),z(*),u(*),xbig(*)
!
!     Set RELACC.
!
ztpar=100.0
relacc=1.0
10 relacc=0.5*relacc
tempa=ztpar+0.5*relacc
tempb=ztpar+relacc
if (ztpar .lt. tempa .and. tempa .lt. tempb) goto 10
!
!     Seek bound inconsistencies and bound equality constraints.
!
meql=0
do 20 i=1,n
if (xl(i) .gt. xu(i)) goto 50
if (xl(i) .eq. xu(i)) meql=meql+1
20 continue
!
!     Initialize U, Z and XBIG.
!
jact=0
nn=n*n
do 30 i=1,nn
30 z(i)=0.0
iz=0
do 40 i=1,n
if (xl(i) .eq. xu(i)) then
    x(i)=xu(i)
    jact=jact+1
    u(jact)=1.0
    iact(jact)=i+m+n
    j=jact
else
    j=i+meql-jact
end if
z(iz+j)=1.0
iz=iz+n
40 xbig(i)=dabs(x(i))
info=1
50 return
end
subroutine ktvec (n,m,a,ia,iact,nact,par,g,reskt,z,u,bres,relaxf, &
  meql,ssqkt,parw,resktw)
implicit real*8 (a-h,o-z)
dimension a(ia,*),iact(*),par(*),g(*),reskt(*),z(*),u(*), &
  bres(*),parw(*),resktw(*)
!
!     Calculate the Lagrange parameters and the residual vector.
!
do 10 i=1,n
10 reskt(i)=g(i)
if (nact .gt. 0) then
    icase=0
20     do 50 kk=1,nact
    k=nact+1-kk
    j=iact(k)
    temp=0.0
    iz=k
    do 30 i=1,n
    temp=temp+z(iz)*reskt(i)
30     iz=iz+n
    temp=temp*u(k)
    if (icase .eq. 0) par(k)=0.0
    if (k .le. meql .or. par(k)+temp .lt. 0.0) then
        par(k)=par(k)+temp
    else
        temp=-par(k)
        par(k)=0.0
    end if
    if (temp .ne. 0.0) then
        if (j .le. m) then
            do 40 i=1,n
40             reskt(i)=reskt(i)-temp*a(i,j)
        else
            jm=j-m
            if (jm .le. n) then
                reskt(jm)=reskt(jm)+temp
            else
                reskt(jm-n)=reskt(jm-n)-temp
            end if
        end if
    end if
50     continue
!
!     Calculate the sum of squares of the KT residual vector.
!
    ssqkt=0.0
    if (nact .eq. n) goto 130
    do 60 i=1,n
60     ssqkt=ssqkt+reskt(i)**2
!
!     Apply iterative refinement to the residual vector.
!
    if (icase .eq. 0) then
        icase=1
        do 70 k=1,nact
70         parw(k)=par(k)
        do 80 i=1,n
80         resktw(i)=reskt(i)
        ssqktw=ssqkt
        goto 20
    end if
!
!     Undo the iterative refinement if it does not reduce SSQKT.
!
    if (ssqktw .lt. ssqkt) then
        do 90 k=1,nact
90         par(k)=parw(k)
        do 100 i=1,n
100         reskt(i)=resktw(i)
        ssqkt=ssqktw
    end if
!
!     Calculate SSQKT when there are no active constraints.
!
else
    ssqkt=0.0
    do 110 i=1,n
110     ssqkt=ssqkt+g(i)**2
end if
!
!     Predict the reduction in F if one corrects any positive residuals
!       of active inequality constraints.
!
relaxf=0.0
if (meql .lt. nact) then
    kl=meql+1
    do 120 k=kl,nact
    j=iact(k)
    if (bres(j) .gt. 0.0) then
        relaxf=relaxf-par(k)*bres(j)
    end if
120     continue
end if
130 return
end
subroutine lsrch (n,x,g,d,xs,gs,relacc,stepcb,ddotg,f,step, &
  nfvals,nfmax,gopt)
implicit real*8 (a-h,o-z)
dimension x(*),g(*),d(*),xs(*),gs(*),gopt(*)
!
!     Initialization.
!
relint=0.9
icount=0
ratio=-1.0
do 10 i=1,n
xs(i)=x(i)
gs(i)=g(i)
gopt(i)=g(i)
if (d(i) .ne. 0.0) then
    temp=dabs(x(i)/d(i))
    if (ratio .lt. 0.0 .or. temp .lt. ratio) ratio=temp
end if
10 continue
step=dmin1(1.0d0,stepcb)
!
!     The following number 1.0D-12 is independent of the working
!       accuracy of the computer arithmetic.
!
stpmin=dmax1(relacc*ratio,1.0d-12*step)
step=dmax1(stpmin,step)
sbase=0.0
fbase=f
ddotgb=ddotg
stplow=0.0
flow=f
dglow=ddotg
stphgh=0.0
stpopt=0.0
fopt=f
dgopt=dabs(ddotg)
!
!     Calculate another function and gradient value.
!
20 do 30 i=1,n
30 x(i)=xs(i)+step*d(i)
call fgcalc (n,x,f,g)
icount=icount+1
dgmid=0.0
do 40 i=1,n
40 dgmid=dgmid+d(i)*g(i)
if (f .le. fopt) then
    if (f .lt. fopt .or. dabs(dgmid) .lt. dgopt) then
        stpopt=step
        fopt=f
        do 50 i=1,n
50         gopt(i)=g(i)
        dgopt=dabs(dgmid)
    end if
end if
if (nfvals+icount .eq. nfmax) goto 70
!
!      Modify the bounds on the steplength or convergence.
!
if (f .ge. fbase+0.1*(step-sbase)*ddotgb) then
    if (stphgh .gt. 0.0 .or. f .gt. fbase .or. dgmid .gt. &
      0.5*ddotg) then
        stphgh=step
        fhgh=f
        dghgh=dgmid
        goto 60
    end if
    sbase=step
    fbase=f
    ddotgb=dgmid
end if
if (dgmid .ge. 0.7*ddotgb) goto 70
stplow=step
flow=f
dglow=dgmid
60 if (stphgh .gt. 0.0 .and. stplow .ge. relint*stphgh) goto 70
!
!     Calculate the next step length or end the iterations.
!
if (stphgh .eq. 0.0) then
    if (step .eq. stepcb) goto 70
    temp=10.0
    if (dgmid .gt. 0.9*ddotg) temp=ddotg/(ddotg-dgmid)
    step=dmin1(temp*step,stepcb)
    goto 20
else if (icount .eq. 1 .or. stplow .gt. 0.0) then
    dgknot=2.0*(fhgh-flow)/(stphgh-stplow)-0.5*(dglow+dghgh)
    if (dgknot .ge. 0.0) then
        ratio=dmax1(0.1d0,0.5d0*dglow/(dglow-dgknot))
    else
        ratio=(0.5*dghgh-dgknot)/(dghgh-dgknot)
    end if
    step=stplow+ratio*(stphgh-stplow)
    goto 20
else
    step=0.1*step
    if (step .ge. stpmin) goto 20
end if
!
!     Return from subroutine.
!
70 if (step .ne. stpopt) then
    step=stpopt
    f=fopt
    do 80 i=1,n
    x(i)=xs(i)+step*d(i)
80     g(i)=gopt(i)
end if
nfvals=nfvals+icount
return
end
!
!     The pentagon problem.
!
implicit real*8 (a-h,o-z)
dimension a(10,15),b(15),xl(6),xu(6),x(6),iact(27),par(20), &
  w(1000)
!
!     The two values of ICASE provide two different values of ACC, the latter
!     accuracy being so demanding that a return with INFO=2 occurs. The
!     final values of the objective function in the two cases agree well
!     and constraint violations are negligible, considering the differences
!     between the final values of the variables.
!
iprint=10
ia=10
n=6
do 100 icase=1,2
acc=1.0d-6
if (icase .eq. 2) acc=1.0d-14
!
!     Set the components of XL, XU and X.
!
do 10 i=1,n
xl(i)=-1.0d6
xu(i)=1.0d6
10 x(i)=0.5d0*dfloat(i-3)
x(2)=0.0d0
x(4)=-1.0d0
x(6)=1.0d0
!
!     Set the constraints.
!
m=0
meq=0
pi=4.0d0*datan(1.0d0)
do 30 k=1,5
theta=0.4d0*dfloat(k-1)*pi
cth=dcos(theta)
sth=dsin(theta)
do 30 j=2,n,2
m=m+1
do 20 i=1,n
20 a(i,m)=0.0d0
a(j-1,m)=cth
a(j,m)=sth
30 b(m)=1.0d0
!
!     Call the optimization package.
!
info=0
print 40, acc,iprint
40 format (//5x,'CALL OF GETMIN WITH  ACC =',1pd11.4, &
  '  AND  IPRINT =',i3)
call getmin (n,m,meq,a,ia,b,xl,xu,x,acc,iact,nact,par,iprint, &
  info,w)
print 50, info
50 format (/5x,'RETURN FROM TOLMIN WITH INFO =',i2)
call fgcalc (n,x,f,w)
print 60, f
60 format (/5x,'FINAL VALUE OF OBJECTIVE FUNCTION =',1pd20.12)
print 70, (x(i),i=1,n)
70 format (/5x,'FINAL COMPONENTS OF X ='//(4x,1p3d20.12))
do 80 k=1,m
do 80 i=1,n
80 b(k)=b(k)-a(i,k)*x(i)
print 90, (b(k),k=1,m)
90 format (/5x,'FINAL CONSTRAINT RESIDUALS ='//(3x,1p6d12.4))
100 continue
stop
end
subroutine minflc (n,m,meq,a,ia,b,xl,xu,x,acc,iact,nact,par, &
  iprint,info,g,z,u,xbig,reskt,bres,d,ztg,gm,xs,gs)
implicit real*8 (a-h,o-z)
dimension a(ia,*),b(*),xl(*),xu(*),x(*),iact(*),par(*),g(*), &
  z(*),u(*),xbig(*),reskt(*),bres(*),d(*),ztg(*),gm(*),xs(*), &
  gs(*)
!
!     Initialize ZZNORM, ITERC, NFVALS and NFMAX.
!
zznorm=-1.0
iterc=0
nfvals=0
nfmax=0
if (info .gt. 0) nfmax=info
!
!     Check the bounds on N, M and MEQ.
!
info=4
if (max0(1-n,-m,meq*(meq-m)) .gt. 0) then
    if (iprint .ne. 0) print 1010
1010     format (/5x,'ERROR RETURN FROM GETMIN BECAUSE A CONDITION', &
      ' ON N, M OR MEQ IS VIOLATED')
    goto 40
end if
!
!     Initialize RELACC, Z, U and TOL.
!
call initzu (n,m,xl,xu,x,iact,meql,info,z,u,xbig,relacc)
tol=dmax1(0.01d0,10.0d0*relacc)
if (info .eq. 4) then
    if (iprint .ne. 0) print 1020
1020     format (/5x,'ERROR RETURN FROM GETMIN BECAUSE A LOWER', &
      ' BOUND EXCEEDS AN UPPER BOUND')
    goto 40
end if
!
!     Add any equality constraints to the active set.
!
if (meq .gt. 0) then
    call eqcons (n,m,meq,a,ia,b,xu,iact,meql,info,z,u,relacc,xs, &
      gs)
    if (info .eq. 5) then
        if (iprint .ne. 0) print 1030
1030         format (/5x,'ERROR RETURN FROM GETMIN BECAUSE THE', &
          ' EQUALITY CONSTRAINTS ARE INCONSISTENT')
        goto 40
    end if
end if
nact=meql
msat=meql
!
!     Add the bounds to the list of constraints.
!
mtot=nact
do 10 i=1,n
if (xl(i) .lt. xu(i)) then
    mtot=mtot+2
    iact(mtot-1)=m+i
    iact(mtot)=m+n+i
end if
10 continue
!
!     Try to satisfy the bound constraints.
!
call getfes (n,m,a,ia,b,xl,xu,x,iact,nact,par,info,g,z,u,xbig, &
  relacc,tol,meql,msat,mtot,bres,d,ztg,gm,reskt,xs,gs)
if (msat .lt. mtot) then
    if (iprint .ne. 0) print 1040
1040     format (/5x,'ERROR RETURN FROM GETMIN BECAUSE THE', &
      ' EQUALITIES AND BOUNDS ARE INCONSISTENT')
    info=6
    goto 40
end if
!
!     Add the ordinary inequalities to the list of constraints.
!
if (m .gt. meq) then
    mp=meq+1
    do 20 k=mp,m
    mtot=mtot+1
20     iact(mtot)=k
end if
!
!     Correct any constraint violations.
!
30 call getfes (n,m,a,ia,b,xl,xu,x,iact,nact,par,info,g,z,u,xbig, &
  relacc,tol,meql,msat,mtot,bres,d,ztg,gm,reskt,xs,gs)
if (msat .lt. mtot) then
    if (iprint .ne. 0) print 1050
1050     format (/5x,'ERROR RETURN FROM GETMIN BECAUSE THE', &
      ' CONSTRAINTS ARE INCONSISTENT')
    info=7
    goto 40
else if (meql .eq. n) then
    if (iprint .ne. 0) print 1060
1060     format (/5x,'GETMIN FINDS THAT THE VARIABLES ARE', &
      ' DETERMINED BY THE EQUALITY CONSTRAINTS')
    goto 40
end if
!
!     Minimize the objective function in the case when constraints are
!       treated as degenerate if their residuals are less than TOL.
!
call minfun (n,m,a,ia,b,xl,xu,x,acc,iact,nact,par,iprint,info,g,z, &
  u,xbig,relacc,zznorm,tol,meql,mtot,iterc,nfvals,nfmax,reskt, &
  bres,d,ztg,gm,xs,gs)
!
!     Reduce TOL if necessary.
!
if (tol .gt. relacc .and. nact .gt. 0) then
    if (nfvals .ne. nfmax) then
        call adjtol (n,m,a,ia,b,xl,xu,x,iact,nact,xbig,relacc,tol, &
          meql)
        goto 30
    else
        info=8
    end if
end if
if (iprint .ne. 0) then
    if (info .eq. 1) print 1070
1070     format (/5x,'GETMIN HAS ACHIEVED THE REQUIRED ACCURACY')
    if (info .eq. 2) print 1080
1080     format (/5x,'GETMIN CAN MAKE NO FURTHER PROGRESS BECAUSE', &
      ' OF ROUNDING ERRORS')
    if (info .eq. 3) print 1090
1090     format (/5x,'GETMIN CAN MAKE NO FURTHER PROGRESS BECAUSE', &
      ' F WILL NOT DECREASE ANY MORE')
    if (info .eq. 8) print 1100
1100     format (/5x,'GETMIN HAS REACHED THE GIVEN LIMIT ON THE', &
      ' NUMBER OF CALLS OF FGCALC')
end if
40 return
end
subroutine minfun (n,m,a,ia,b,xl,xu,x,acc,iact,nact,par,iprint, &
  info,g,z,u,xbig,relacc,zznorm,tol,meql,mtot,iterc,nfvals, &
  nfmax,reskt,bres,d,ztg,gm,xs,gs)
implicit real*8 (a-h,o-z)
dimension a(ia,*),b(*),xl(*),xu(*),x(*),iact(*),par(*),g(*),z(*), &
  u(*),xbig(*),reskt(*),bres(*),d(*),ztg(*),gm(*),xs(*),gs(*)
save f
!
!     Initialize the minimization calculation.
!
msat=mtot
iterk=iterc
nfvalk=nfvals
if (nfvals .eq. 0 .or. info .eq. 1) then
    call fgcalc (n,x,f,g)
    nfvals=nfvals+1
end if
fprev=dabs(f+f+1.0)
iterp=-1
if (iprint .ne. 0) then
    print 1000, tol
1000     format (/5x,'NEW VALUE OF TOL =',1pd13.5)
    iterp=iterc+iabs(iprint)
    if (iterc .eq. 0) iterp=0
end if
!
!     Calculate the next search direction.
!
10 call conres (n,m,a,ia,b,xl,xu,x,iact,nact,par,g,z,u,xbig,bres,d, &
  ztg,relacc,tol,stepcb,ddotg,meql,msat,mtot,indxbd,gm,reskt,xs, &
  gs)
!
!     Calculate the Kuhn Tucker residual vector.
!
call ktvec (n,m,a,ia,iact,nact,par,g,reskt,z,u,bres,relaxf,meql, &
  ssqkt,xs,gs)
!
!     Test for convergence.
!
if (ssqkt .le. acc*acc) then
    info=1
    goto 70
end if
if (ddotg .ge. 0.0) then
    info=2
    goto 70
end if
!
!     Test for termination due to no decrease in F.
!
if (f .ge. fprev) then
    if (tol .eq. relacc .or. nact .eq. 0) then
        if (diff .gt. 0.0) goto 20
    end if
    info=3
    goto 70
end if
20 diff=fprev-f
fprev=f
!
!     Test that more calls of FGCALC are allowed.
!
if (nfvals .eq. nfmax) then
    info=8
    goto 70
end if
!
!     Test whether to reduce TOL and to provide printing.
!
if (tol .gt. relacc .and. iterc .gt. iterk .and. &
  0.1*relaxf .ge. dmax1(diff,-0.5d0*ddotg)) goto 70
if (iterp .eq. iterc) goto 80
!
!     Calculate the step along the search direction.
!
40 iterc=iterc+1
call lsrch (n,x,g,d,xs,gs,relacc,stepcb,ddotg,f,step,nfvals, &
  nfmax,bres)
if (step .eq. 0.0) then
    info=3
    sum=0.0
    do 50 i=1,n
50     sum=sum+dabs(d(i)*gs(i))
    if (ddotg+relacc*sum .ge. 0.0) info=2
    goto 70
end if
!
!     Revise XBIG.
!
do 60 i=1,n
60 xbig(i)=dmax1(xbig(i),dabs(x(i)))
!
!     Revise the second derivative approximation.
!
call zbfgs (n,x,nact,g,z,ztg,xs,gs,zznorm)
!
!     Add a constraint to the active set if it restricts the step.
!
if (step .eq. stepcb) then
    k=iact(indxbd)
    if (k .gt. m) then
        k=k-m
        if (k .le. n) then
            x(k)=xl(k)
        else
            x(k-n)=xu(k-n)
        end if
    end if
    call addcon (n,m,a,ia,iact,nact,z,u,relacc,indxbd,xs,gs)
end if
goto 10
!
!     Printing from the subroutine.
!
70 if (iprint .eq. 0) goto 90
iterk=-1
80 print 1010, iterc,nfvals,f
1010 format (/5x,'ITERS =',i4,5x,'F.VALS =',i4,5x,'F =',1pd15.7)
print 1020, (x(i),i=1,n)
1020 format ('  X =',(1p5d14.5))
print 1030, (g(i),i=1,n)
1030 format ('  G =',(1p5d14.5))
if (iprint .lt. 0) then
    if (nact .eq. 0) then
        print 1050
1050         format (5x,'NO ACTIVE CONSTRAINTS')
    else
        print 1060, (iact(i),i=1,nact)
1060         format (' IA =',(14i5))
        print 1070, (par(i),i=1,nact)
1070         format (' LP =',(1p5d14.5))
    end if
    if (nact .eq. n) then
        print 1080
1080         format (5x,'KT RESIDUAL VECTOR IS ZERO')
    else
        print 1090,(reskt(i),i=1,n)
1090         format (' KT =',(1p5d14.5))
    end if
end if
iterp=iterc+iabs(iprint)
if (iterk .ge. 0) goto 40
90 return
end
subroutine newcon (n,m,a,ia,iact,nact,z,u,d,relacc,mdeg,zzdiag, &
  gmnew,cgrad)
implicit real*8 (a-h,o-z)
dimension a(ia,*),iact(*),z(*),u(*),d(*),zzdiag(*),gmnew(*), &
  cgrad(*)
!
!     Initialization.
!
np=nact+1
khigh=mdeg
iz=0
do 20 i=1,n
zzdiag(i)=0.0
do 10 j=np,n
10 zzdiag(i)=zzdiag(i)+z(iz+j)**2
20 iz=iz+n
!
!     Calculate the scalar products of D with its constraints.
!
30 cvmax=0.0
do 50 k=np,khigh
j=iact(k)
if (j .le. m) then
    sum=0.0
    sumabs=0.0
    sumd=0.0
    do 40 i=1,n
    temp=d(i)*a(i,j)
    sum=sum+temp
    sumabs=sumabs+dabs(temp)
40     sumd=sumd+zzdiag(i)*a(i,j)**2
else
    jm=j-m
    if (jm .le. n) then
        sum=-d(jm)
    else
        jm=jm-n
        sum=d(jm)
    end if
    sumabs=dabs(sum)
    sumd=zzdiag(jm)
end if
!
!     Pick out the most violated constraint, or return if the
!       violation is negligible.
!
if (sum .gt. relacc*sumabs) then
    cviol=sum*sum/sumd
    if (cviol .gt. cvmax) then
        cvmax=cviol
        iadd=k
        savsum=sum
        savabs=sumabs
    end if
end if
50 continue
if (cvmax .le. 0.0) goto 140
if (nact .eq. 0) goto 120
!
!     Set GMNEW to the gradient of the most violated constraint.
!
j=iact(iadd)
if (j .le. m) then
    jmv=0
    do 60 i=1,n
60     gmnew(i)=a(i,j)
else
    jmv=j-m
    do 70 i=1,n
70     gmnew(i)=0.0
    if (jmv .le. n) then
        gmnew(jmv)=-1.0
    else
        jmv=jmv-n
        gmnew(jmv)=1.0
    end if
end if
!
!     Modify GMNEW for the next active constraint.
!
k=nact
80 temp=0.0
iz=k
do 90 i=1,n
temp=temp+z(iz)*gmnew(i)
90 iz=iz+n
temp=temp*u(k)
j=iact(k)
if (j .le. m) then
    do 100 i=1,n
100     gmnew(i)=gmnew(i)-temp*a(i,j)
else
    jm=j-m
    if (jm .le. n) then
        gmnew(jm)=gmnew(jm)+temp
    else
        gmnew(jm-n)=gmnew(jm-n)-temp
    end if
end if
!
!     Revise the values of SAVSUM and SAVABS.
!
sum=0.0
sumabs=0.0
do 110 i=1,n
temp=d(i)*gmnew(i)
sum=sum+temp
110 sumabs=sumabs+dabs(temp)
savsum=dmin1(savsum,sum)
savabs=dmax1(savabs,sumabs)
k=k-1
if (k .ge. 1) goto 80
!
!     Add the new constraint to the active set if the constraint
!       violation is still significant.
!
if (jmv .gt. 0) d(jmv)=0.0
if (savsum .le. relacc*savabs) goto 130
120 k=nact
call addcon (n,m,a,ia,iact,nact,z,u,relacc,iadd,gmnew,cgrad)
if (nact .gt. k) goto 140
!
!     Seek another constraint violation.
!
iadd=np
130 if (np .lt. khigh) then
    k=iact(khigh)
    iact(khigh)=iact(iadd)
    iact(iadd)=k
    khigh=khigh-1
    goto 30
end if
140 return
end
subroutine satact (n,m,a,ia,b,xl,xu,x,iact,nact,info,z,u,xbig, &
  relacc,tol,meql)
implicit real*8 (a-h,o-z)
dimension a(ia,*),b(*),xl(*),xu(*),x(*),iact(*),z(*),u(*), &
  xbig(*)
if (nact .eq. 0) goto 50
do 30 k=1,nact
!
!     Calculate the next constraint residual.
!
j=iact(k)
if (j .le. m) then
    res=b(j)
    resabs=dabs(b(j))
    resbig=resabs
    do 10 i=1,n
    tempa=a(i,j)
    temp=tempa*x(i)
    res=res-temp
    resabs=resabs+dabs(temp)
10     resbig=resbig+dabs(tempa)*xbig(i)
else
    jx=j-m
    if (jx .le. n) then
        res=x(jx)-xl(jx)
        resabs=dabs(x(jx))+dabs(xl(jx))
        resbig=xbig(jx)+dabs(xl(jx))
        savex=xl(jx)
    else
        jx=jx-n
        res=xu(jx)-x(jx)
        resabs=dabs(x(jx))+dabs(xu(jx))
        resbig=xbig(jx)+dabs(xu(jx))
        savex=xu(jx)
    end if
end if
!
!     Shift X if necessary.
!
if (res .ne. 0.0) then
    temp=res/resabs
    if (k .le. meql) temp=-dabs(temp)
    if (tol .eq. relacc .or. temp+relacc .lt. 0.0) then
        info=1
        scale=res*u(k)
        iz=k
        do 20 i=1,n
        x(i)=x(i)+scale*z(iz)
        iz=iz+n
20         xbig(i)=dmax1(xbig(i),dabs(x(i)))
        if (j .gt. m) x(jx)=savex
!
!     Else flag a constraint deletion if necessary.
!
    else if (res/resbig .gt. tol) then
        iact(k)=-iact(k)
    end if
end if
30 continue
!
!     Delete any flagged constraints and then return.
!
idrop=nact
40 if (iact(idrop) .lt. 0) then
    iact(idrop)=-iact(idrop)
    call delcon (n,m,a,ia,iact,nact,z,u,relacc,idrop)
end if
idrop=idrop-1
if (idrop .gt. meql) goto 40
50 return
end
subroutine sdegen (n,m,a,ia,iact,nact,par,z,u,d,ztg,gm,relacc, &
  ddotgm,meql,mdeg,gmnew,parnew,cgrad)
implicit real*8 (a-h,o-z)
dimension a(ia,*),iact(*),par(*),z(*),u(*),d(*),ztg(*),gm(*), &
  gmnew(*),parnew(*),cgrad(*)
mp=meql+1
dtest=0.0
!
!     Calculate the search direction and branch if it is not downhill.
!
10 call sdirn (n,nact,z,d,ztg,gm,relacc,ddotgm)
if (ddotgm .eq. 0.0) goto 120
!
!     Branch if there is no need to consider any degenerate constraints.
!     The test gives termination if two consecutive additions to the
!       active set fail to increase the predicted new value of F.
!
if (nact .eq. mdeg) goto 120
np=nact+1
sum=0.0
do 20 j=np,n
20 sum=sum+ztg(j)**2
if (dtest .gt. 0.0 .and. sum .ge. dtest) then
    if (itest .eq. 1) goto 120
    itest=1
else
    dtest=sum
    itest=0
end if
!
!     Add a constraint to the active set if there are any significant
!       violations of degenerate constraints.
!
k=nact
call newcon (n,m,a,ia,iact,nact,z,u,d,relacc,mdeg,gmnew,parnew, &
  cgrad)
if (nact .eq. k) goto 120
par(nact)=0.0
!
!     Calculate the new reduced gradient and Lagrange parameters.
!
30 do 40 i=1,n
40 gmnew(i)=gm(i)
k=nact
50 temp=0.0
iz=k
do 60 i=1,n
temp=temp+z(iz)*gmnew(i)
60 iz=iz+n
temp=temp*u(k)
parnew(k)=par(k)+temp
if (k .eq. nact) parnew(k)=dmin1(parnew(k),0.0d0)
j=iact(k)
if (j .le. m) then
    do 70 i=1,n
70     gmnew(i)=gmnew(i)-temp*a(i,j)
else
    jm=j-m
    if (jm .le. n) then
        gmnew(jm)=gmnew(jm)+temp
    else
        gmnew(jm-n)=gmnew(jm-n)-temp
    end if
end if
k=k-1
if (k .gt. meql) goto 50
!
!     Set RATIO for linear interpolation between PAR and PARNEW.
!
ratio=0.0
if (mp .lt. nact) then
    ku=nact-1
    do 80 k=mp,ku
    if (parnew(k) .gt. 0.0) then
        ratio=parnew(k)/(parnew(k)-par(k))
        idrop=k
    end if
80     continue
end if
!
!     Apply the linear interpolation.
!
theta=1.0-ratio
do 90 k=mp,nact
90 par(k)=dmin1(theta*parnew(k)+ratio*par(k),0.0d0)
do 100 i=1,n
100 gm(i)=theta*gmnew(i)+ratio*gm(i)
!
!     Drop a constraint if RATIO is positive.
!
if (ratio .gt. 0.0) then
    call delcon (n,m,a,ia,iact,nact,z,u,relacc,idrop)
    do 110 k=idrop,nact
110     par(k)=par(k+1)
    goto 30
end if
!
!     Return if there is no freedom for a new search direction.
!
if (nact .lt. n) goto 10
ddotgm=0.0
120 return
end
subroutine sdirn (n,nact,z,d,ztg,gm,relacc,ddotgm)
implicit real*8 (a-h,o-z)
dimension z(*),d(*),ztg(*),gm(*)
ddotgm=0.0
if (nact .ge. n) goto 60
!
!     Premultiply GM by the transpose of Z.
!
np=nact+1
do 20 j=np,n
sum=0.0
sumabs=0.0
iz=j
do 10 i=1,n
temp=z(iz)*gm(i)
sum=sum+temp
sumabs=sumabs+dabs(temp)
10 iz=iz+n
if (dabs(sum) .le. relacc*sumabs) sum=0.0
20 ztg(j)=sum
!
!     Form D by premultiplying ZTG by -Z.
!
iz=0
do 40 i=1,n
sum=0.0
sumabs=0.0
do 30 j=np,n
temp=z(iz+j)*ztg(j)
sum=sum-temp
30 sumabs=sumabs+dabs(temp)
if (dabs(sum) .le. relacc*sumabs) sum=0.0
d(i)=sum
40 iz=iz+n
!
!     Test that the search direction is downhill.
!
sumabs=0.0
do 50 i=1,n
temp=d(i)*gm(i)
ddotgm=ddotgm+temp
50 sumabs=sumabs+dabs(temp)
if (ddotgm+relacc*sumabs .ge. 0.0) ddotgm=0.0
60 return
end
subroutine stepbd (n,m,a,ia,iact,bres,d,stepcb,ddotg,mdeg,msat, &
  mtot,indxbd)
implicit real*8 (a-h,o-z)
dimension a(ia,*),iact(*),bres(*),d(*)
!
!     Set steps to constraint boundaries and find the least positive one.
!
iflag=0
stepcb=0.0
indxbd=0
k=mdeg
10 k=k+1
if (k .gt. mtot) goto 40
!
!     Form the scalar product of D with the current constraint normal.
!
20     j=iact(k)
    if (j .le. m) then
        sp=0.0
        do 30 i=1,n
30         sp=sp+d(i)*a(i,j)
    else
        jm=j-m
        if (jm .le. n) then
            sp=-d(jm)
        else
            sp=d(jm-n)
        end if
    end if
!
!     The next branch is taken if label 20 was reached via label 50.
!
    if (iflag .eq. 1) goto 60
!
!     Set BRES(J) to indicate the status of the j-th constraint.
!
    if (sp*bres(j) .le. 0.0) then
        bres(j)=0.0
    else
        bres(j)=bres(j)/sp
        if (stepcb .eq. 0.0 .or. bres(j) .lt. stepcb) then
            stepcb=bres(j)
            indxbd=k
        end if
    end if
    go to 10
40 continue
!
!     Try to pass through the boundary of a violated constraint.
!
50 if (indxbd .le. msat) goto 80
    iflag=1
    k=indxbd
    goto 20
60     msat=msat+1
    iact(indxbd)=iact(msat)
    iact(msat)=j
    bres(j)=0.0
    indxbd=msat
    ddotg=ddotg-sp
    if (ddotg .lt. 0.0 .and. msat .lt. mtot) then
!
!     Seek the next constraint boundary along the search direction.
!
        temp=0.0
        kl=mdeg+1
        do 70 k=kl,mtot
        j=iact(k)
        if (bres(j) .gt. 0.0) then
            if (temp .eq. 0.0 .or. bres(j) .lt. temp) then
                temp=bres(j)
                indxbd=k
            end if
        end if
70         continue
        if (temp .gt. 0.0) then
            stepcb=temp
            goto 50
        end if
    end if
80 continue
return
end
subroutine zbfgs (n,x,nact,g,z,ztg,xs,gs,zznorm)
implicit real*8 (a-h,o-z)
dimension x(*),g(*),z(*),ztg(*),xs(*),gs(*)
!
!     Test if there is sufficient convexity for the update.
!
dd=0.0
dg=0.0
temp=0.0
do 10 i=1,n
xs(i)=x(i)-xs(i)
dd=dd+xs(i)**2
temp=temp+gs(i)*xs(i)
gs(i)=g(i)-gs(i)
10 dg=dg+gs(i)*xs(i)
if (dg .lt. 0.1*dabs(temp)) goto 90
!
!     Transform the Z matrix.
!
k=n
20 kp=k
k=k-1
if (k .gt. nact) then
    if (ztg(kp) .eq. 0.0) goto 20
    temp=dabs(ztg(kp))*dsqrt(1.0+(ztg(k)/ztg(kp))**2)
    wcos=ztg(k)/temp
    wsin=ztg(kp)/temp
    ztg(k)=temp
    iz=k
    do 30 i=1,n
    temp=wcos*z(iz+1)-wsin*z(iz)
    z(iz)=wcos*z(iz)+wsin*z(iz+1)
    z(iz+1)=temp
30     iz=iz+n
    goto 20
end if
!
!     Update the value of ZZNORM.
!
if (zznorm .lt. 0.0) then
    zznorm=dd/dg
else
    temp=dsqrt(zznorm*dd/dg)
    zznorm=dmin1(zznorm,temp)
    zznorm=dmax1(zznorm,0.1d0*temp)
end if
!
!     Complete the updating of Z.
!
np=nact+1
temp=dsqrt(dg)
iz=np
do 40 i=1,n
z(iz)=xs(i)/temp
40 iz=iz+n
if (np .lt. n) then
    km=np+1
    do 80 k=km,n
    temp=0.0
    iz=k
    do 50 i=1,n
    temp=temp+gs(i)*z(iz)
50     iz=iz+n
    temp=temp/dg
    sum=0.0
    iz=k
    do 60 i=1,n
    z(iz)=z(iz)-temp*xs(i)
    sum=sum+z(iz)**2
60     iz=iz+n
    if (sum .lt. zznorm) then
        temp=dsqrt(zznorm/sum)
        iz=k
        do 70 i=1,n
        z(iz)=temp*z(iz)
70         iz=iz+n
    end if
80     continue
end if
90 return
end
