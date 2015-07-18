subroutine calfun (n,x,f)
implicit real*8 (a-h,o-z)
common fmax
dimension x(*)
zero=0.0d0
f=fmax
v12=x(1)*x(5)-x(4)*x(2)
v13=x(1)*x(8)-x(7)*x(2)
v14=x(1)*x(11)-x(10)*x(2)
v23=x(4)*x(8)-x(7)*x(5)
v24=x(4)*x(11)-x(10)*x(5)
v34=x(7)*x(11)-x(10)*x(8)
del1=v23*x(12)-v24*x(9)+v34*x(6)
if (del1 .le. zero) goto 10
del2=-v34*x(3)-v13*x(12)+v14*x(9)
if (del2 .le. zero) goto 10
del3=-v14*x(6)+v24*x(3)+v12*x(12)
if (del3 .le. zero) goto 10
del4=-v12*x(9)+v13*x(6)-v23*x(3)
if (del4 .le. zero) goto 10
temp=(del1+del2+del3+del4)**3/(del1*del2*del3*del4)
f=dmin1(temp/6.0d0,fmax)
10 continue
return
end
subroutine getact (n,m,amat,b,nact,iact,qfac,rfac,snorm, &
  resnew,resact,g,dw,vlam,w)
implicit real*8 (a-h,o-z)
dimension amat(n,*),b(*),iact(*),qfac(n,*),rfac(*), &
  resnew(*),resact(*),g(*),dw(*),vlam(*),w(*)
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
one=1.0d0
tiny=1.0d-60
zero=0.0d0
tdel=0.2d0*snorm
ddsav=zero
do 10 i=1,n
ddsav=ddsav+g(i)**2
10 vlam(i)=zero
ddsav=ddsav+ddsav
!
!     Set the initial QFAC to the identity matrix in the case NACT=0.
!
if (nact .eq. 0) then
    do 30 i=1,n
    do 20 j=1,n
20     qfac(i,j)=zero
30     qfac(i,i)=one
    goto 100
end if
!
!     Remove any constraints from the initial active set whose residuals
!       exceed TDEL.
!
iflag=1
ic=nact
40 if (resact(ic) .gt. tdel) goto 800
50 ic=ic-1
if (ic .gt. 0) goto 40
!
!     Remove any constraints from the initial active set whose Lagrange
!       multipliers are nonnegative, and set the surviving multipliers.
!
iflag=2
60 if (nact .eq. 0) goto 100
ic=nact
70 temp=zero
do 80 i=1,n
80 temp=temp+qfac(i,ic)*g(i)
idiag=(ic*ic+ic)/2
if (ic .lt. nact) then
    jw=idiag+ic
    do 90 j=ic+1,nact
    temp=temp-rfac(jw)*vlam(j)
90     jw=jw+j
end if
if (temp .ge. zero) goto 800
vlam(ic)=temp/rfac(idiag)
ic=ic-1
if (ic .gt. 0) goto 70
!
!     Set the new search direction D. Terminate if the 2-norm of D is zero
!       or does not decrease, or if NACT=N holds. The situation NACT=N
!       occurs for sufficiently large SNORM if the origin is in the convex
!       hull of the constraint gradients.
!
100 if (nact .eq. n) goto 290
do 110 j=nact+1,n
w(j)=zero
do 110 i=1,n
110 w(j)=w(j)+qfac(i,j)*g(i)
dd=zero
do 130 i=1,n
dw(i)=zero
do 120 j=nact+1,n
120 dw(i)=dw(i)-w(j)*qfac(i,j)
130 dd=dd+dw(i)**2
if (dd .ge. ddsav) goto 290
if (dd .eq. zero) goto 300
ddsav=dd
dnorm=dsqrt(dd)
!
!     Pick the next integer L or terminate, a positive value of L being
!       the index of the most violated constraint. The purpose of CTOL
!       below is to estimate whether a positive value of VIOLMX may be
!       due to computer rounding errors.
!
l=0
if (m .gt. 0) then
    test=dnorm/snorm
    violmx=zero
    do 150 j=1,m
    if (resnew(j) .gt. zero .and. resnew(j) .le. tdel) then
        sum=zero
        do 140 i=1,n
140         sum=sum+amat(i,j)*dw(i)
        if (sum .gt. test*resnew(j)) then
            if (sum .gt. violmx) then
                l=j
                violmx=sum
            end if
        end if
    end if
150     continue
    ctol=zero
    temp=0.01d0*dnorm
    if (violmx .gt. zero .and. violmx .lt. temp) then
        if (nact .gt. 0) then
            do 170 k=1,nact
            j=iact(k)
            sum=zero
            do 160 i=1,n
160             sum=sum+dw(i)*amat(i,j)
170             ctol=dmax1(ctol,dabs(sum))
        end if
    end if
end if
w(1)=one
if (l .eq. 0) goto 300
if (violmx .le. 10.0d0*ctol) goto 300
!
!     Apply Givens rotations to the last (N-NACT) columns of QFAC so that
!       the first (NACT+1) columns of QFAC are the ones required for the
!       addition of the L-th constraint, and add the appropriate column
!       to RFAC.
!
nactp=nact+1
idiag=(nactp*nactp-nactp)/2
rdiag=zero
do 200 j=n,1,-1
sprod=zero
do 180 i=1,n
180 sprod=sprod+qfac(i,j)*amat(i,l)
if (j .le. nact) then
    rfac(idiag+j)=sprod
else
    if (dabs(rdiag) .le. 1.0d-20*dabs(sprod)) then
        rdiag=sprod
    else
        temp=dsqrt(sprod*sprod+rdiag*rdiag)
        cosv=sprod/temp
        sinv=rdiag/temp
        rdiag=temp
        do 190 i=1,n
        temp=cosv*qfac(i,j)+sinv*qfac(i,j+1)
        qfac(i,j+1)=-sinv*qfac(i,j)+cosv*qfac(i,j+1)
190         qfac(i,j)=temp
    end if
end if
200 continue
if (rdiag .lt. zero) then
    do 210 i=1,n
210     qfac(i,nactp)=-qfac(i,nactp)
end if
rfac(idiag+nactp)=dabs(rdiag)
nact=nactp
iact(nact)=l
resact(nact)=resnew(l)
vlam(nact)=zero
resnew(l)=zero
!
!     Set the components of the vector VMU in W.
!
220 w(nact)=one/rfac((nact*nact+nact)/2)**2
if (nact .gt. 1) then
    do 240 i=nact-1,1,-1
    idiag=(i*i+i)/2
    jw=idiag+i
    sum=zero
    do 230 j=i+1,nact
    sum=sum-rfac(jw)*w(j)
230     jw=jw+j
240     w(i)=sum/rfac(idiag)
end if
!
!     Calculate the multiple of VMU to subtract from VLAM, and update VLAM.
!
vmult=violmx
ic=0
j=1
250 if (j .lt. nact) then
    if (vlam(j) .ge. vmult*w(j)) then
        ic=j
        vmult=vlam(j)/w(j)
    end if
    j=j+1
    goto 250
end if
do 260 j=1,nact
260 vlam(j)=vlam(j)-vmult*w(j)
if (ic .gt. 0) vlam(ic)=zero
violmx=dmax1(violmx-vmult,zero)
if (ic .eq. 0) violmx=zero
!
!     Reduce the active set if necessary, so that all components of the
!       new VLAM are negative, with resetting of the residuals of the
!       constraints that become inactive.
!
iflag=3
ic=nact
270 if (vlam(ic) .lt. zero) goto 280
resnew(iact(ic))=dmax1(resact(ic),tiny)
goto 800
280 ic=ic-1
if (ic .gt. 0) goto 270
!
!     Calculate the next VMU if VIOLMX is positive. Return if NACT=N holds,
!       as then the active constraints imply D=0. Otherwise, go to label
!       100, to calculate the new D and to test for termination.
!
if (violmx .gt. zero) goto 220
if (nact .lt. n) goto 100
290 dd=zero
300 w(1)=dd
return
!
!     These instructions rearrange the active constraints so that the new
!       value of IACT(NACT) is the old value of IACT(IC). A sequence of
!       Givens rotations is applied to the current QFAC and RFAC. Then NACT
!       is reduced by one.
!
800 resnew(iact(ic))=dmax1(resact(ic),tiny)
jc=ic
810 if (jc .lt. nact) then
    jcp=jc+1
    idiag=jc*jcp/2
    jw=idiag+jcp
    temp=dsqrt(rfac(jw-1)**2+rfac(jw)**2)
    cval=rfac(jw)/temp
    sval=rfac(jw-1)/temp
    rfac(jw-1)=sval*rfac(idiag)
    rfac(jw)=cval*rfac(idiag)
    rfac(idiag)=temp
    if (jcp .lt. nact) then
        do 820 j=jcp+1,nact
        temp=sval*rfac(jw+jc)+cval*rfac(jw+jcp)
        rfac(jw+jcp)=cval*rfac(jw+jc)-sval*rfac(jw+jcp)
        rfac(jw+jc)=temp
820         jw=jw+j
    end if
    jdiag=idiag-jc
    do 830 i=1,n
    if (i .lt. jc) then
        temp=rfac(idiag+i)
        rfac(idiag+i)=rfac(jdiag+i)
        rfac(jdiag+i)=temp
    end if
    temp=sval*qfac(i,jc)+cval*qfac(i,jcp)
    qfac(i,jcp)=cval*qfac(i,jc)-sval*qfac(i,jcp)
830     qfac(i,jc)=temp
    iact(jc)=iact(jcp)
    resact(jc)=resact(jcp)
    vlam(jc)=vlam(jcp)
    jc=jcp
    goto 810
end if
nact=nact-1
goto (50,60,280),iflag
end
subroutine lincoa (n,npt,m,a,ia,b,x,rhobeg,rhoend,iprint, &
  maxfun,w)
implicit real*8 (a-h,o-z)
dimension a(ia,*),b(*),x(*),w(*)
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
!
!     Check that N, NPT and MAXFUN are acceptable.
!
zero=0.0d0
smallx=1.0d-6*rhoend
np=n+1
nptm=npt-np
if (n .le. 1) then
    print 10
10     format (/4x,'Return from LINCOA because N is less than 2.')
    goto 80
end if
if (npt .lt. n+2 .or. npt .gt. ((n+2)*np)/2) then
    print 20
20     format (/4x,'Return from LINCOA because NPT is not in', &
      ' the required interval.')
    goto 80
end if
if (maxfun .le. npt) then
    print 30
30     format (/4x,'Return from LINCOA because MAXFUN is less', &
      ' than NPT+1.')
    goto 80
end if
!
!     Normalize the constraints, and copy the resultant constraint matrix
!       and right hand sides into working space, after increasing the right
!       hand sides if necessary so that the starting point is feasible.
!
iamat=max0(m+3*n,2*m+n,2*npt)+1
ib=iamat+m*n
iflag=0
if (m .gt. 0) then
    iw=iamat-1
    do 60 j=1,m
    sum=zero
    temp=zero
    do 40 i=1,n
    sum=sum+a(i,j)*x(i)
40     temp=temp+a(i,j)**2
    if (temp .eq. zero) then
        print 50
50         format (/4x,'Return from LINCOA because the gradient of', &
          ' a constraint is zero.')
        goto 80
    end if
    temp=dsqrt(temp)
    if (sum-b(j) .gt. smallx*temp) iflag=1
    w(ib+j-1)=dmax1(b(j),sum)/temp
    do 60 i=1,n
    iw=iw+1
60     w(iw)=a(i,j)/temp
end if
if (iflag .eq. 1) then
    if (iprint .gt. 0) print 70
70     format (/4x,'LINCOA has made the initial X feasible by', &
      ' increasing part(s) of B.')
end if
!
!     Partition the working space array, so that different parts of it can be
!     treated separately by the subroutine that performs the main calculation.
!
ndim=npt+n
ixb=ib+m
ixp=ixb+n
ifv=ixp+n*npt
ixs=ifv+npt
ixo=ixs+n
igo=ixo+n
ihq=igo+n
ipq=ihq+(n*np)/2
ibmat=ipq+npt
izmat=ibmat+ndim*n
istp=izmat+npt*nptm
isp=istp+n
ixn=isp+npt+npt
iac=ixn+n
irc=iac+n
iqf=irc+m
irf=iqf+n*n
ipqw=irf+(n*np)/2
!
!     The above settings provide a partition of W for subroutine LINCOB.
!
call lincob (n,npt,m,w(iamat),w(ib),x,rhobeg,rhoend,iprint, &
  maxfun,w(ixb),w(ixp),w(ifv),w(ixs),w(ixo),w(igo),w(ihq), &
  w(ipq),w(ibmat),w(izmat),ndim,w(istp),w(isp),w(ixn),w(iac), &
  w(irc),w(iqf),w(irf),w(ipqw),w)
80 return
end
subroutine lincob (n,npt,m,amat,b,x,rhobeg,rhoend,iprint, &
  maxfun,xbase,xpt,fval,xsav,xopt,gopt,hq,pq,bmat,zmat,ndim, &
  step,sp,xnew,iact,rescon,qfac,rfac,pqw,w)
implicit real*8 (a-h,o-z)
dimension amat(n,*),b(*),x(*),xbase(*),xpt(npt,*),fval(*), &
  xsav(*),xopt(*),gopt(*),hq(*),pq(*),bmat(ndim,*), &
  zmat(npt,*),step(*),sp(*),xnew(*),iact(*),rescon(*), &
  qfac(n,*),rfac(*),pqw(*),w(*)
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
half=0.5d0
one=1.0d0
tenth=0.1d0
zero=0.0d0
np=n+1
nh=(n*np)/2
nptm=npt-np
!
!     Set the elements of XBASE, XPT, FVAL, XSAV, XOPT, GOPT, HQ, PQ, BMAT,
!       ZMAT and SP for the first iteration. An important feature is that,
!       if the interpolation point XPT(K,.) is not feasible, where K is any
!       integer from [1,NPT], then a change is made to XPT(K,.) if necessary
!       so that the constraint violation is at least 0.2*RHOBEG. Also KOPT
!       is set so that XPT(KOPT,.) is the initial trust region centre.
!
call prelim (n,npt,m,amat,b,x,rhobeg,iprint,xbase,xpt,fval, &
  xsav,xopt,gopt,kopt,hq,pq,bmat,zmat,idz,ndim,sp,rescon, &
  step,pqw,w)
!
!     Begin the iterative procedure.
!
nf=npt
fopt=fval(kopt)
rho=rhobeg
delta=rho
ifeas=0
nact=0
itest=3
10 knew=0
nvala=0
nvalb=0
!
!     Shift XBASE if XOPT may be too far from XBASE. First make the changes
!       to BMAT that do not depend on ZMAT.
!
20 fsave=fopt
xoptsq=zero
do 30 i=1,n
30 xoptsq=xoptsq+xopt(i)**2
if (xoptsq .ge. 1.0d4*delta*delta) then
    qoptsq=0.25d0*xoptsq
    do 50 k=1,npt
    sum=zero
    do 40 i=1,n
40     sum=sum+xpt(k,i)*xopt(i)
    sum=sum-half*xoptsq
    w(npt+k)=sum
    sp(k)=zero
    do 50 i=1,n
    xpt(k,i)=xpt(k,i)-half*xopt(i)
    step(i)=bmat(k,i)
    w(i)=sum*xpt(k,i)+qoptsq*xopt(i)
    ip=npt+i
    do 50 j=1,i
50     bmat(ip,j)=bmat(ip,j)+step(i)*w(j)+w(i)*step(j)
!
!     Then the revisions of BMAT that depend on ZMAT are calculated.
!
    do 90 k=1,nptm
    sumz=zero
    do 60 i=1,npt
    sumz=sumz+zmat(i,k)
60     w(i)=w(npt+i)*zmat(i,k)
    do 80 j=1,n
    sum=qoptsq*sumz*xopt(j)
    do 70 i=1,npt
70     sum=sum+w(i)*xpt(i,j)
    step(j)=sum
    if (k .lt. idz) sum=-sum
    do 80 i=1,npt
80     bmat(i,j)=bmat(i,j)+sum*zmat(i,k)
    do 90 i=1,n
    ip=i+npt
    temp=step(i)
    if (k .lt. idz) temp=-temp
    do 90 j=1,i
90     bmat(ip,j)=bmat(ip,j)+temp*step(j)
!
!     Update the right hand sides of the constraints.
!
    if (m .gt. 0) then
        do 110 j=1,m
        temp=zero
        do 100 i=1,n
100         temp=temp+amat(i,j)*xopt(i)
110         b(j)=b(j)-temp
    end if
!
!     The following instructions complete the shift of XBASE, including the
!       changes to the parameters of the quadratic model.
!
    ih=0
    do 130 j=1,n
    w(j)=zero
    do 120 k=1,npt
    w(j)=w(j)+pq(k)*xpt(k,j)
120     xpt(k,j)=xpt(k,j)-half*xopt(j)
    do 130 i=1,j
    ih=ih+1
    hq(ih)=hq(ih)+w(i)*xopt(j)+xopt(i)*w(j)
130     bmat(npt+i,j)=bmat(npt+j,i)
    do 140 j=1,n
    xbase(j)=xbase(j)+xopt(j)
    xopt(j)=zero
140     xpt(kopt,j)=zero
end if
!
!     In the case KNEW=0, generate the next trust region step by calling
!       TRSTEP, where SNORM is the current trust region radius initially.
!       The final value of SNORM is the length of the calculated step,
!       except that SNORM is zero on return if the projected gradient is
!       unsuitable for starting the conjugate gradient iterations.
!
delsav=delta
ksave=knew
if (knew .eq. 0) then
    snorm=delta
    do 150 i=1,n
150     xnew(i)=gopt(i)
    call trstep (n,npt,m,amat,b,xpt,hq,pq,nact,iact,rescon, &
      qfac,rfac,snorm,step,xnew,w,w(m+1),pqw,pqw(np),w(m+np))
!
!     A trust region step is applied whenever its length, namely SNORM, is at
!       least HALF*DELTA. It is also applied if its length is at least 0.1999
!       times DELTA and if a line search of TRSTEP has caused a change to the
!       active set. Otherwise there is a branch below to label 530 or 560.
!
    temp=half*delta
    if (xnew(1) .ge. half) temp=0.1999d0*delta
    if (snorm .le. temp) then
        delta=half*delta
        if (delta .le. 1.4d0*rho) delta=rho
        nvala=nvala+1
        nvalb=nvalb+1
        temp=snorm/rho
        if (delsav .gt. rho) temp=one
        if (temp .ge. half) nvala=zero
        if (temp .ge. tenth) nvalb=zero
        if (delsav .gt. rho) goto 530
        if (nvala .lt. 5 .and. nvalb .lt. 3) goto 530
        if (snorm .gt. zero) ksave=-1
        goto 560
    end if
    nvala=zero
    nvalb=zero
!
!     Alternatively, KNEW is positive. Then the model step is calculated
!       within a trust region of radius DEL, after setting the gradient at
!       XBASE and the second derivative parameters of the KNEW-th Lagrange
!       function in W(1) to W(N) and in PQW(1) to PQW(NPT), respectively.
!
else
    del=dmax1(tenth*delta,rho)
    do 160 i=1,n
160     w(i)=bmat(knew,i)
    do 170 k=1,npt
170     pqw(k)=zero
    do 180 j=1,nptm
    temp=zmat(knew,j)
    if (j .lt. idz) temp=-temp
    do 180 k=1,npt
180     pqw(k)=pqw(k)+temp*zmat(k,j)
    call qmstep (n,npt,m,amat,b,xpt,xopt,nact,iact,rescon, &
      qfac,kopt,knew,del,step,w,pqw,w(np),w(np+m),ifeas)
end if
!
!     Set VQUAD to the change to the quadratic model when the move STEP is
!       made from XOPT. If STEP is a trust region step, then VQUAD should be
!       negative. If it is nonnegative due to rounding errors in this case,
!       there is a branch to label 530 to try to improve the model.
!
vquad=zero
ih=0
do 190 j=1,n
vquad=vquad+step(j)*gopt(j)
do 190 i=1,j
ih=ih+1
temp=step(i)*step(j)
if (i .eq. j) temp=half*temp
190 vquad=vquad+temp*hq(ih)
do 210 k=1,npt
temp=zero
do 200 j=1,n
temp=temp+xpt(k,j)*step(j)
200 sp(npt+k)=temp
210 vquad=vquad+half*pq(k)*temp*temp
if (ksave .eq. 0 .and. vquad .ge. zero) goto 530
!
!     Calculate the next value of the objective function. The difference
!       between the actual new value of F and the value predicted by the
!       model is recorded in DIFF.
!
220 nf=nf+1
if (nf .gt. maxfun) then
    nf=nf-1
    if (iprint .gt. 0) print 230
230     format (/4x,'Return from LINCOA because CALFUN has been', &
      ' called MAXFUN times.')
    goto 600
end if
xdiff=zero
do 240 i=1,n
xnew(i)=xopt(i)+step(i)
x(i)=xbase(i)+xnew(i)
240 xdiff=xdiff+(x(i)-xsav(i))**2
xdiff=dsqrt(xdiff)
if (ksave .eq. -1) xdiff=rho
if (xdiff .le. tenth*rho .or. xdiff .ge. delta+delta) then
    ifeas=0
    if (iprint .gt. 0) print 250
250     format (/4x,'Return from LINCOA because rounding errors', &
      ' prevent reasonable changes to X.')
    goto 600
end if
if (ksave .le. 0) ifeas=1
f=dfloat(ifeas)
call calfun (n,x,f)
if (iprint .eq. 3) then
    print 260, nf,f,(x(i),i=1,n)
260     format (/4x,'Function number',i6,'    F =',1pd18.10, &
      '    The corresponding X is:'/(2x,5d15.6))
end if
if (ksave .eq. -1) goto 600
diff=f-fopt-vquad
!
!     If X is feasible, then set DFFALT to the difference between the new
!       value of F and the value predicted by the alternative model.
!
if (ifeas .eq. 1 .and. itest .lt. 3) then
    do 270 k=1,npt
    pqw(k)=zero
270     w(k)=fval(k)-fval(kopt)
    do 290 j=1,nptm
    sum=zero
    do 280 i=1,npt
280     sum=sum+w(i)*zmat(i,j)
    if (j .lt. idz) sum=-sum
    do 290 k=1,npt
290     pqw(k)=pqw(k)+sum*zmat(k,j)
    vqalt=zero
    do 310 k=1,npt
    sum=zero
    do 300 j=1,n
300     sum=sum+bmat(k,j)*step(j)
    vqalt=vqalt+sum*w(k)
310     vqalt=vqalt+pqw(k)*sp(npt+k)*(half*sp(npt+k)+sp(k))
    dffalt=f-fopt-vqalt
end if
if (itest .eq. 3) then
    dffalt=diff
    itest=0
end if
!
!     Pick the next value of DELTA after a trust region step.
!
if (ksave .eq. 0) then
    ratio=(f-fopt)/vquad
    if (ratio .le. tenth) then
        delta=half*delta
    else if (ratio .le. 0.7d0) then
        delta=dmax1(half*delta,snorm)
    else 
        temp=dsqrt(2.0d0)*delta
        delta=dmax1(half*delta,snorm+snorm)
        delta=dmin1(delta,temp)
    end if
    if (delta .le. 1.4d0*rho) delta=rho
end if
!
!     Update BMAT, ZMAT and IDZ, so that the KNEW-th interpolation point
!       can be moved. If STEP is a trust region step, then KNEW is zero at
!       present, but a positive value is picked by subroutine UPDATE.
!
call update (n,npt,xpt,bmat,zmat,idz,ndim,sp,step,kopt, &
  knew,pqw,w)
if (knew .eq. 0) then
    if (iprint .gt. 0) print 320
320     format (/4x,'Return from LINCOA because the denominator' &
      ' of the updating formula is zero.')
    goto 600
end if
!
!     If ITEST is increased to 3, then the next quadratic model is the
!       one whose second derivative matrix is least subject to the new
!       interpolation conditions. Otherwise the new model is constructed
!       by the symmetric Broyden method in the usual way.
!
if (ifeas .eq. 1) then
    itest=itest+1
    if (dabs(dffalt) .ge. tenth*dabs(diff)) itest=0
end if
!
!     Update the second derivatives of the model by the symmetric Broyden
!       method, using PQW for the second derivative parameters of the new
!       KNEW-th Lagrange function. The contribution from the old parameter
!       PQ(KNEW) is included in the second derivative matrix HQ. W is used
!       later for the gradient of the new KNEW-th Lagrange function.       
!
if (itest .lt. 3) then
    do 330 k=1,npt
330     pqw(k)=zero
    do 350 j=1,nptm
    temp=zmat(knew,j)
    if (temp .ne. zero) then
        if (j .lt. idz) temp=-temp
        do 340 k=1,npt
340         pqw(k)=pqw(k)+temp*zmat(k,j)
    end if
350     continue
    ih=0
    do 360 i=1,n
    w(i)=bmat(knew,i)
    temp=pq(knew)*xpt(knew,i)
    do 360 j=1,i
    ih=ih+1
360     hq(ih)=hq(ih)+temp*xpt(knew,j)
    pq(knew)=zero
    do 370 k=1,npt
370     pq(k)=pq(k)+diff*pqw(k)
end if
!
!     Include the new interpolation point with the corresponding updates of
!       SP. Also make the changes of the symmetric Broyden method to GOPT at
!       the old XOPT if ITEST is less than 3.
!
fval(knew)=f
sp(knew)=sp(kopt)+sp(npt+kopt)
ssq=zero
do 380 i=1,n
xpt(knew,i)=xnew(i)
380 ssq=ssq+step(i)**2
sp(npt+knew)=sp(npt+kopt)+ssq
if (itest .lt. 3) then
    do 390 k=1,npt
    temp=pqw(k)*sp(k)
    do 390 i=1,n
390     w(i)=w(i)+temp*xpt(k,i)
    do 400 i=1,n
400     gopt(i)=gopt(i)+diff*w(i)
end if
!
!     Update FOPT, XSAV, XOPT, KOPT, RESCON and SP if the new F is the
!       least calculated value so far with a feasible vector of variables.
!
if (f .lt. fopt .and. ifeas .eq. 1) then
    fopt=f
    do 410 j=1,n
    xsav(j)=x(j)
410     xopt(j)=xnew(j)
    kopt=knew
    snorm=dsqrt(ssq)
    do 430 j=1,m
    if (rescon(j) .ge. delta+snorm) then
        rescon(j)=snorm-rescon(j)
    else
        rescon(j)=rescon(j)+snorm
        if (rescon(j)+delta .gt. zero) then
            temp=b(j)
            do 420 i=1,n
420             temp=temp-xopt(i)*amat(i,j)
            temp=dmax1(temp,zero)
            if (temp .ge. delta) temp=-temp
            rescon(j)=temp
        end if
    end if
430     continue
    do 440 k=1,npt
440     sp(k)=sp(k)+sp(npt+k)
!
!     Also revise GOPT when symmetric Broyden updating is applied.
!
    if (itest .lt. 3) then
        ih=0
        do 450 j=1,n
        do 450 i=1,j
        ih=ih+1
        if (i .lt. j) gopt(j)=gopt(j)+hq(ih)*step(i)
450         gopt(i)=gopt(i)+hq(ih)*step(j)
        do 460 k=1,npt
        temp=pq(k)*sp(npt+k)
        do 460 i=1,n
460         gopt(i)=gopt(i)+temp*xpt(k,i)
    end if
end if
!
!     Replace the current model by the least Frobenius norm interpolant if
!       this interpolant gives substantial reductions in the predictions
!       of values of F at feasible points.
!
if (itest .eq. 3) then
    do 470 k=1,npt
    pq(k)=zero
470     w(k)=fval(k)-fval(kopt)
    do 490 j=1,nptm
    sum=zero
    do 480 i=1,npt
480     sum=sum+w(i)*zmat(i,j)
    if (j .lt. idz) sum=-sum
    do 490 k=1,npt
490     pq(k)=pq(k)+sum*zmat(k,j)
    do 500 j=1,n
    gopt(j)=zero
    do 500 i=1,npt
500     gopt(j)=gopt(j)+w(i)*bmat(i,j)
    do 510 k=1,npt
    temp=pq(k)*sp(k)
    do 510 i=1,n
510     gopt(i)=gopt(i)+temp*xpt(k,i)
    do 520 ih=1,nh
520     hq(ih)=zero
end if
!
!     If a trust region step has provided a sufficient decrease in F, then
!       branch for another trust region calculation. Every iteration that
!       takes a model step is followed by an attempt to take a trust region
!       step.
!
knew=0
if (ksave .gt. 0) goto 20
if (ratio .ge. tenth) goto 20
!
!     Alternatively, find out if the interpolation points are close enough
!       to the best point so far.
!
530 distsq=dmax1(delta*delta,4.0d0*rho*rho)
do 550 k=1,npt
sum=zero
do 540 j=1,n
540 sum=sum+(xpt(k,j)-xopt(j))**2
if (sum .gt. distsq) then
    knew=k
    distsq=sum
end if
550 continue
!
!     If KNEW is positive, then branch back for the next iteration, which
!       will generate a "model step". Otherwise, if the current iteration
!       has reduced F, or if DELTA was above its lower bound when the last
!       trust region step was calculated, then try a "trust region" step
!       instead.
!
if (knew .gt. 0) goto 20
knew=0
if (fopt .lt. fsave) goto 20
if (delsav .gt. rho) goto 20
!
!     The calculations with the current value of RHO are complete.
!       Pick the next value of RHO.
!
560 if (rho .gt. rhoend) then
    delta=half*rho
    if (rho .gt. 250.0d0*rhoend) then
        rho=tenth*rho
    else if (rho .le. 16.0d0*rhoend) then
        rho=rhoend
    else
        rho=dsqrt(rho*rhoend)
    end if 
    delta=dmax1(delta,rho)
    if (iprint .ge. 2) then
        if (iprint .ge. 3) print 570
570         format (5x)
        print 580, rho,nf
580         format (/4x,'New RHO =',1pd11.4,5x,'Number of', &
          ' function values =',i6)
        print 590, fopt,(xbase(i)+xopt(i),i=1,n)
590         format (4x,'Least value of F =',1pd23.15,9x, &
          'The corresponding X is:'/(2x,5d15.6))
    end if
    goto 10
end if
!
!     Return from the calculation, after branching to label 220 for another
!       Newton-Raphson step if it has not been tried before.
!
if (ksave .eq. -1) goto 220
600 if (fopt .le. f .or. ifeas .eq. 0) then
    do 610 i=1,n
610     x(i)=xsav(i)
    f=fopt
end if
if (iprint .ge. 1) then
    print 620, nf
620     format (/4x,'At the return from LINCOA',5x, &
      'Number of function values =',i6)
    print 590, f,(x(i),i=1,n)
end if
w(1)=f
w(2)=dfloat(nf)+half
return
end
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
implicit real*8 (a-h,o-z)
common fmax
dimension xp(50),yp(50),zp(50),a(12,200),b(200),x(12),w(500000)
!
!     Set some constants.
!
one=1.0d0
two=2.0d0
zero=0.0d0
pi=4.0d0*datan(one)
ia=12
n=12
!
!     Set the data points.
!
np=50
sumx=zero
sumy=zero
sumz=zero
do 10 j=1,np
theta=dfloat(j-1)*pi/dfloat(np-1)
xp(j)=dcos(theta)*dcos(two*theta)
sumx=sumx+xp(j)
yp(j)=dsin(theta)*dcos(two*theta)
sumy=sumy+yp(j)
zp(j)=dsin(two*theta)
10 sumz=sumz+zp(j)
sumx=sumx/dfloat(np)
sumy=sumy/dfloat(np)
sumz=sumz/dfloat(np)
do 20 j=1,np
xp(j)=xp(j)-sumx
yp(j)=yp(j)-sumy
20 zp(j)=zp(j)-sumz
!
!     Set the linear constraints.
!
m=4*np
do 30 k=1,m
b(k)=one
do 30 i=1,n
30 a(i,k)=zero
do 40 j=1,np
do 40 i=1,4
k=4*j+i-4
iw=3*i
a(iw-2,k)=xp(j)
a(iw-1,k)=yp(j)
40 a(iw,k)=zp(j)
!
!     Set the initial vector of variables. The JCASE=1,6 loop gives six
!       different choices of NPT when LINCOA is called.
!
xs=zero
ys=zero
zs=zero
ss=zero
do 50 j=1,np
xs=dmin1(xs,xp(j))
ys=dmin1(ys,yp(j))
zs=dmin1(zs,zp(j))
50 ss=dmax1(ss,xp(j)+yp(j)+zp(j))
fmax=(ss-xs-ys-zs)**3/6.0d0
do 80 jcase=1,6
do 60 i=2,8
60 x(i)=zero
x(1)=one/xs
x(5)=one/ys
x(9)=one/zs
x(10)=one/ss
x(11)=one/ss
x(12)=one/ss
!
!     Call of LINCOA, which provides the printing given at the end of this
!       note.
!
npt=5*jcase+10
rhobeg=1.0d0
rhoend=1.0d-6
iprint=1
maxfun=10000
print 70, npt,rhoend
70 format (//4x,'Output from LINCOA with  NPT =',i4, &
  '  and  RHOEND =',1pd12.4)
call lincoa (n,npt,m,a,ia,b,x,rhobeg,rhoend,iprint,maxfun,w)
80 continue
stop
end
subroutine prelim (n,npt,m,amat,b,x,rhobeg,iprint,xbase, &
  xpt,fval,xsav,xopt,gopt,kopt,hq,pq,bmat,zmat,idz,ndim, &
  sp,rescon,step,pqw,w)
implicit real*8 (a-h,o-z)
dimension amat(n,*),b(*),x(*),xbase(*),xpt(npt,*),fval(*), &
  xsav(*),xopt(*),gopt(*),hq(*),pq(*),bmat(ndim,*),zmat(npt,*), &
  sp(*),rescon(*),step(*),pqw(*),w(*)
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
half=0.5d0
one=1.0d0
zero=0.0d0
nptm=npt-n-1
rhosq=rhobeg*rhobeg
recip=one/rhosq
reciq=dsqrt(half)/rhosq
test=0.2d0*rhobeg
idz=1
kbase=1
!
!     Set the initial elements of XPT, BMAT, SP and ZMAT to zero. 
!
do 20 j=1,n
xbase(j)=x(j)
do 10 k=1,npt
10 xpt(k,j)=zero
do 20 i=1,ndim
20 bmat(i,j)=zero
do 30 k=1,npt
sp(k)=zero
do 30 j=1,npt-n-1
30 zmat(k,j)=zero
!
!     Set the nonzero coordinates of XPT(K,.), K=1,2,...,min[2*N+1,NPT],
!       but they may be altered later to make a constraint violation
!       sufficiently large. The initial nonzero elements of BMAT and of
!       the first min[N,NPT-N-1] columns of ZMAT are set also.
!
do 40 j=1,n
xpt(j+1,j)=rhobeg
if (j .lt. npt-n) then
    jp=n+j+1
    xpt(jp,j)=-rhobeg
    bmat(j+1,j)=half/rhobeg
    bmat(jp,j)=-half/rhobeg
    zmat(1,j)=-reciq-reciq
    zmat(j+1,j)=reciq
    zmat(jp,j)=reciq
else
    bmat(1,j)=-one/rhobeg
    bmat(j+1,j)=one/rhobeg
    bmat(npt+j,j)=-half*rhosq
end if
40 continue
!
!     Set the remaining initial nonzero elements of XPT and ZMAT when the
!       number of interpolation points exceeds 2*N+1.
!
if (npt .gt. 2*n+1) then
    do 50 k=n+1,npt-n-1
    itemp=(k-1)/n
    ipt=k-itemp*n
    jpt=ipt+itemp
    if (jpt .gt. n) jpt=jpt-n
    xpt(n+k+1,ipt)=rhobeg
    xpt(n+k+1,jpt)=rhobeg
    zmat(1,k)=recip
    zmat(ipt+1,k)=-recip
    zmat(jpt+1,k)=-recip
50     zmat(n+k+1,k)=recip
end if
!
!     Update the constraint right hand sides to allow for the shift XBASE.
!
if (m .gt. 0) then
    do 70 j=1,m
    temp=zero
    do 60 i=1,n
60     temp=temp+amat(i,j)*xbase(i)
70     b(j)=b(j)-temp
end if
!
!     Go through the initial points, shifting every infeasible point if
!       necessary so that its constraint violation is at least 0.2*RHOBEG.
!
do 150 nf=1,npt
feas=one
bigv=zero
j=0
80 j=j+1
if (j .le. m .and. nf .ge. 2) then
    resid=-b(j)
    do 90 i=1,n
90     resid=resid+xpt(nf,i)*amat(i,j)
    if (resid .le. bigv) goto 80
    bigv=resid
    jsav=j
    if (resid .le. test) then
        feas=-one
        goto 80
    end if
    feas=zero
end if
if (feas .lt. zero) then
    do 100 i=1,n
100     step(i)=xpt(nf,i)+(test-bigv)*amat(i,jsav)
    do 110 k=1,npt
    sp(npt+k)=zero
    do 110 j=1,n
110     sp(npt+k)=sp(npt+k)+xpt(k,j)*step(j)
    call update (n,npt,xpt,bmat,zmat,idz,ndim,sp,step, &
      kbase,nf,pqw,w)
    do 120 i=1,n
120     xpt(nf,i)=step(i)
end if
!
!     Calculate the objective function at the current interpolation point,
!       and set KOPT to the index of the first trust region centre.
!
do 130 j=1,n
130 x(j)=xbase(j)+xpt(nf,j)
f=feas
call calfun (n,x,f)
if (iprint .eq. 3) then
    print 140, nf,f,(x(i),i=1,n)
140     format (/4x,'Function number',i6,'    F =',1pd18.10, &
      '    The corresponding X is:'/(2x,5d15.6))
end if
if (nf .eq. 1) then
    kopt=1
else if (f .lt. fval(kopt) .and. feas .gt. zero) then
    kopt=nf
end if
150 fval(nf)=f
!
!     Set PQ for the first quadratic model.
!
do 160 j=1,nptm
w(j)=zero
do 160 k=1,npt
160 w(j)=w(j)+zmat(k,j)*fval(k)
do 170 k=1,npt
pq(k)=zero
do 170 j=1,nptm
170 pq(k)=pq(k)+zmat(k,j)*w(j)
!
!     Set XOPT, SP, GOPT and HQ for the first quadratic model.
!
do 180 j=1,n
xopt(j)=xpt(kopt,j)
xsav(j)=xbase(j)+xopt(j)
180 gopt(j)=zero
do 200 k=1,npt
sp(k)=zero
do 190 j=1,n
190 sp(k)=sp(k)+xpt(k,j)*xopt(j)
temp=pq(k)*sp(k)
do 200 j=1,n
200 gopt(j)=gopt(j)+fval(k)*bmat(k,j)+temp*xpt(k,j)
do 210 i=1,(n*n+n)/2
210 hq(i)=zero
!
!     Set the initial elements of RESCON.
!
do 230 j=1,m
temp=b(j)
do 220 i=1,n
220 temp=temp-xopt(i)*amat(i,j)
temp=dmax1(temp,zero)
if (temp .ge. rhobeg) temp=-temp
230 rescon(j)=temp  
return
end
subroutine qmstep (n,npt,m,amat,b,xpt,xopt,nact,iact, &
  rescon,qfac,kopt,knew,del,step,gl,pqw,rstat,w,ifeas)
implicit real*8 (a-h,o-z)
dimension amat(n,*),b(*),xpt(npt,*),xopt(*),iact(*), &
  rescon(*),qfac(n,*),step(*),gl(*),pqw(*),rstat(*),w(*)
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
half=0.5d0
one=1.0d0
tenth=0.1d0
zero=0.0d0
test=0.2d0*del
!
!     Replace GL by the gradient of LFUNC at the trust region centre, and
!       set the elements of RSTAT.
!
do 20 k=1,npt
temp=zero
do 10 j=1,n
10 temp=temp+xpt(k,j)*xopt(j)
temp=pqw(k)*temp
do 20 i=1,n
20 gl(i)=gl(i)+temp*xpt(k,i)
if (m .gt. 0) then
    do 30 j=1,m
    rstat(j)=one
30     if (dabs(rescon(j)) .ge. del) rstat(j)=-one
    do 40 k=1,nact
40     rstat(iact(k))=zero
end if
!
!     Find the greatest modulus of LFUNC on a line through XOPT and
!       another interpolation point within the trust region.
!
iflag=0
vbig=zero
do 60 k=1,npt
if (k .eq. kopt) goto 60
ss=zero
sp=zero
do 50 i=1,n
temp=xpt(k,i)-xopt(i)
ss=ss+temp*temp
50 sp=sp+gl(i)*temp
stp=-del/dsqrt(ss)
if (k .eq. knew) then
    if (sp*(sp-one) .lt. zero) stp=-stp
    vlag=dabs(stp*sp)+stp*stp*dabs(sp-one)
else
    vlag=dabs(stp*(one-stp)*sp)
end if
if (vlag .gt. vbig) then
    ksav=k
    stpsav=stp
    vbig=vlag
end if
60 continue
!
!     Set STEP to the move that gives the greatest modulus calculated above.
!       This move may be replaced by a steepest ascent step from XOPT.
!
gg=zero
do 70 i=1,n
gg=gg+gl(i)**2
70 step(i)=stpsav*(xpt(ksav,i)-xopt(i))
vgrad=del*dsqrt(gg)
if (vgrad .le. tenth*vbig) goto 220
!
!     Make the replacement if it provides a larger value of VBIG.
!
ghg=zero
do 90 k=1,npt
temp=zero
do 80 j=1,n
80 temp=temp+xpt(k,j)*gl(j)
90 ghg=ghg+pqw(k)*temp*temp
vnew=vgrad+dabs(half*del*del*ghg/gg)
if (vnew .gt. vbig) then
    vbig=vnew
    stp=del/dsqrt(gg)
    if (ghg .lt. zero) stp=-stp
    do 100 i=1,n
100     step(i)=stp*gl(i)
end if
if (nact .eq. 0 .or. nact .eq. n) goto 220
!
!     Overwrite GL by its projection. Then set VNEW to the greatest
!       value of |LFUNC| on the projected gradient from XOPT subject to
!       the trust region bound. If VNEW is sufficiently large, then STEP
!       may be changed to a move along the projected gradient.
!
do 110 k=nact+1,n
w(k)=zero
do 110 i=1,n
110 w(k)=w(k)+gl(i)*qfac(i,k)
gg=zero
do 130 i=1,n
gl(i)=zero
do 120 k=nact+1,n
120 gl(i)=gl(i)+qfac(i,k)*w(k)
130 gg=gg+gl(i)**2
vgrad=del*dsqrt(gg)
if (vgrad .le. tenth*vbig) goto 220
ghg=zero
do 150 k=1,npt
temp=zero
do 140 j=1,n
140 temp=temp+xpt(k,j)*gl(j)
150 ghg=ghg+pqw(k)*temp*temp
vnew=vgrad+dabs(half*del*del*ghg/gg)
!
!     Set W to the possible move along the projected gradient.
!
stp=del/dsqrt(gg)
if (ghg .lt. zero) stp=-stp
ww=zero
do 160 i=1,n
w(i)=stp*gl(i)
160 ww=ww+w(i)**2
!
!     Set STEP to W if W gives a sufficiently large value of the modulus
!       of the Lagrange function, and if W either preserves feasibility
!       or gives a constraint violation of at least 0.2*DEL. The purpose
!       of CTOL below is to provide a check on feasibility that includes
!       a tolerance for contributions from computer rounding errors.
!
if (vnew/vbig .ge. 0.2d0) then
    ifeas=1
    bigv=zero
    j=0
170     j=j+1
    if (j .le. m) then
        if (rstat(j) .eq. one) then
            temp=-rescon(j)
            do 180 i=1,n
180             temp=temp+w(i)*amat(i,j)
            bigv=dmax1(bigv,temp)
        end if
        if (bigv .lt. test) goto 170
        ifeas=0
    end if
    ctol=zero
    temp=0.01d0*dsqrt(ww)
    if (bigv .gt. zero .and. bigv .lt. temp) then
        do 200 k=1,nact
        j=iact(k)
        sum=zero
        do 190 i=1,n
190         sum=sum+w(i)*amat(i,j)
200         ctol=dmax1(ctol,dabs(sum))
    end if
    if (bigv .le. 10.0d0*ctol .or. bigv .ge. test) then
        do 210 i=1,n
210         step(i)=w(i)
        goto 260
    end if
end if
!
!     Calculate the greatest constraint violation at XOPT+STEP with STEP at
!       its original value. Modify STEP if this violation is unacceptable.
!
220 ifeas=1
bigv=zero
resmax=zero
j=0
230 j=j+1
if (j .le. m) then
    if (rstat(j) .lt. zero) goto 230
    temp=-rescon(j)
    do 240 i=1,n
240     temp=temp+step(i)*amat(i,j)
    resmax=dmax1(resmax,temp)
    if (temp .lt. test) then
        if (temp .le. bigv) goto 230
        bigv=temp
        jsav=j
        ifeas=-1
        goto 230
    end if
    ifeas=0
end if
if (ifeas .eq. -1) then
    do 250 i=1,n
250     step(i)=step(i)+(test-bigv)*amat(i,jsav)
    ifeas=0
end if
!
!     Return the calculated STEP and the value of IFEAS.
!
260 return
end
subroutine trstep (n,npt,m,amat,b,xpt,hq,pq,nact,iact,rescon, &
  qfac,rfac,snorm,step,g,resnew,resact,d,dw,w)
implicit real*8 (a-h,o-z)
dimension amat(n,*),b(*),xpt(npt,*),hq(*),pq(*),iact(*), &
  rescon(*),qfac(n,*),rfac(*),step(*),g(*),resnew(*),resact(*), &
  d(*),dw(*),w(*)
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
half=0.5d0
one=1.0d0
tiny=1.0d-60
zero=0.0d0
ctest=0.01d0
snsq=snorm*snorm
!
!     Set the initial elements of RESNEW, RESACT and STEP.
!
if (m .gt. 0) then
    do 10 j=1,m
    resnew(j)=rescon(j)
    if (rescon(j) .ge. snorm) then
        resnew(j)=-one
    else if (rescon(j) .ge. zero) then
        resnew(j)=dmax1(resnew(j),tiny)
    end if
10     continue
    if (nact .gt. 0) then
        do 20 k=1,nact
        resact(k)=rescon(iact(k))
20         resnew(iact(k))=zero
    end if
end if
do 30 i=1,n
30 step(i)=zero
ss=zero
reduct=zero
ncall=0
!
!     GETACT picks the active set for the current STEP. It also sets DW to
!       the vector closest to -G that is orthogonal to the normals of the
!       active constraints. DW is scaled to have length 0.2*SNORM, as then
!       a move of DW from STEP is allowed by the linear constraints.
!
40 ncall=ncall+1
call getact (n,m,amat,b,nact,iact,qfac,rfac,snorm,resnew, &
  resact,g,dw,w,w(n+1))
if (w(n+1) .eq. zero) goto 320
scale=0.2d0*snorm/dsqrt(w(n+1))
do 50 i=1,n
50 dw(i)=scale*dw(i)
!
!     If the modulus of the residual of an active constraint is substantial,
!       then set D to the shortest move from STEP to the boundaries of the
!       active constraints.
!
resmax=zero
if (nact .gt. 0) then
    do 60 k=1,nact
60     resmax=dmax1(resmax,resact(k))
end if
gamma=zero
if (resmax .gt. 1.0d-4*snorm) then
    ir=0
    do 80 k=1,nact
    temp=resact(k)
    if (k .ge. 2) then
        do 70 i=1,k-1
        ir=ir+1
70         temp=temp-rfac(ir)*w(i)
    end if
    ir=ir+1
80     w(k)=temp/rfac(ir)
    do 90 i=1,n
    d(i)=zero
    do 90 k=1,nact
90     d(i)=d(i)+w(k)*qfac(i,k)
!
!     The vector D that has just been calculated is also the shortest move
!       from STEP+DW to the boundaries of the active constraints. Set GAMMA
!       to the greatest steplength of this move that satisfies the trust
!       region bound.
!
    rhs=snsq
    ds=zero
    dd=zero
    do 100 i=1,n
    sum=step(i)+dw(i)
    rhs=rhs-sum*sum
    ds=ds+d(i)*sum
100     dd=dd+d(i)**2
    if (rhs .gt. zero) then
        temp=dsqrt(ds*ds+dd*rhs)
        if (ds .le. zero) then
            gamma=(temp-ds)/dd
        else
            gamma=rhs/(temp+ds)
        end if
    end if
!
!     Reduce the steplength GAMMA if necessary so that the move along D
!       also satisfies the linear constraints.
!
    j=0
110     if (gamma .gt. zero) then
        j=j+1
        if (resnew(j) .gt. zero) then
            ad=zero
            adw=zero
            do 120 i=1,n
            ad=ad+amat(i,j)*d(i)
120             adw=adw+amat(i,j)*dw(i)
            if (ad .gt. zero) then
                temp=dmax1((resnew(j)-adw)/ad,zero)
                gamma=dmin1(gamma,temp)
            end if
        end if
        if (j .lt. m) goto 110
    end if
    gamma=dmin1(gamma,one)
end if
!
!     Set the next direction for seeking a reduction in the model function
!       subject to the trust region bound and the linear constraints.
!
if (gamma .le. zero) then
    do 130 i=1,n
130     d(i)=dw(i)
    icount=nact
else
    do 140 i=1,n
140     d(i)=dw(i)+gamma*d(i)
    icount=nact-1
end if
alpbd=one
!
!     Set ALPHA to the steplength from STEP along D to the trust region
!       boundary. Return if the first derivative term of this step is
!       sufficiently small or if no further progress is possible.
!
150 icount=icount+1
rhs=snsq-ss
if (rhs .le. zero) goto 320
dg=zero
ds=zero
dd=zero
do 160 i=1,n
dg=dg+d(i)*g(i)
ds=ds+d(i)*step(i)
160 dd=dd+d(i)**2
if (dg .ge. zero) goto 320
temp=dsqrt(rhs*dd+ds*ds)
if (ds .le. zero) then
    alpha=(temp-ds)/dd
else
    alpha=rhs/(temp+ds)
end if
if (-alpha*dg .le. ctest*reduct) goto 320
!
!     Set DW to the change in gradient along D.
!
ih=0
do 170 j=1,n
dw(j)=zero
do 170 i=1,j
ih=ih+1
if (i .lt. j) dw(j)=dw(j)+hq(ih)*d(i)
170 dw(i)=dw(i)+hq(ih)*d(j)
do 190 k=1,npt
temp=zero
do 180 j=1,n
180 temp=temp+xpt(k,j)*d(j)
temp=pq(k)*temp
do 190 i=1,n
190 dw(i)=dw(i)+temp*xpt(k,i)
!
!     Set DGD to the curvature of the model along D. Then reduce ALPHA if
!       necessary to the value that minimizes the model.
!
dgd=zero
do 200 i=1,n
200 dgd=dgd+d(i)*dw(i)
alpht=alpha
if (dg+alpha*dgd .gt. zero) then
    alpha=-dg/dgd
end if
!
!     Make a further reduction in ALPHA if necessary to preserve feasibility,
!       and put some scalar products of D with constraint gradients in W.
!
alphm=alpha
jsav=0
if (m .gt. 0) then
    do 220 j=1,m
    ad=zero
    if (resnew(j) .gt. zero) then
        do 210 i=1,n
210         ad=ad+amat(i,j)*d(i)
        if (alpha*ad .gt. resnew(j)) then
            alpha=resnew(j)/ad
            jsav=j
        end if
    end if
220     w(j)=ad
end if
alpha=dmax1(alpha,alpbd)
alpha=dmin1(alpha,alphm)
if (icount .eq. nact) alpha=dmin1(alpha,one)
!
!     Update STEP, G, RESNEW, RESACT and REDUCT.
!
ss=zero
do 230 i=1,n
step(i)=step(i)+alpha*d(i)
ss=ss+step(i)**2
230 g(i)=g(i)+alpha*dw(i)
if (m .gt. 0) then
    do 240 j=1,m
    if (resnew(j) .gt. zero) then
        resnew(j)=dmax1(resnew(j)-alpha*w(j),tiny)
    end if
240     continue
end if
if (icount .eq. nact .and. nact .gt. 0) then
    do 250 k=1,nact
250     resact(k)=(one-gamma)*resact(k)
end if
reduct=reduct-alpha*(dg+half*alpha*dgd)
!
!     Test for termination. Branch to label 40 if there is a new active
!       constraint and if the distance from STEP to the trust region
!       boundary is at least 0.2*SNORM.
!
if (alpha .eq. alpht) goto 320
temp=-alphm*(dg+half*alphm*dgd)
if (temp .le. ctest*reduct) goto 320
if (jsav .gt. 0) then
    if (ss .le. 0.64d0*snsq) goto 40
    goto 320
end if
if (icount .eq. n) goto 320
!
!     Calculate the next search direction, which is conjugate to the
!       previous one except in the case ICOUNT=NACT.
!
if (nact .gt. 0) then
    do 260 j=nact+1,n
    w(j)=zero
    do 260 i=1,n
260     w(j)=w(j)+g(i)*qfac(i,j)
    do 280 i=1,n
    temp=zero
    do 270 j=nact+1,n
270     temp=temp+qfac(i,j)*w(j)
280     w(n+i)=temp
else
    do 290 i=1,n
290     w(n+i)=g(i)
end if
if (icount .eq. nact) then
    beta=zero
else
    wgd=zero
    do 300 i=1,n
300     wgd=wgd+w(n+i)*dw(i)
    beta=wgd/dgd
end if
do 310 i=1,n
310 d(i)=-w(n+i)+beta*d(i)
alpbd=zero
goto 150
!
!     Return from the subroutine.
!
320 snorm=zero
if (reduct .gt. zero) snorm=dsqrt(ss)
g(1)=zero
if (ncall .gt. 1) g(1)=one
return
end
subroutine update (n,npt,xpt,bmat,zmat,idz,ndim,sp,step, &
  kopt,knew,vlag,w)
implicit real*8 (a-h,o-z)
dimension xpt(npt,*),bmat(ndim,*),zmat(npt,*),sp(*),step(*), &
  vlag(*),w(*)
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
half=0.5d0
one=1.0d0
zero=0.0d0
nptm=npt-n-1
!
!     Calculate VLAG and BETA for the current choice of STEP. The first NPT
!       elements of VLAG are set to the values of the Lagrange functions at
!       XPT(KOPT,.)+STEP(.). The first NPT components of W_check are held
!       in W, where W_check is defined in a paper on the updating method.
!
do 20 k=1,npt
w(k)=sp(npt+k)*(half*sp(npt+k)+sp(k))
sum=zero
do 10 j=1,n
10 sum=sum+bmat(k,j)*step(j)
20 vlag(k)=sum
beta=zero
do 40 k=1,nptm
sum=zero
do 30 i=1,npt
30 sum=sum+zmat(i,k)*w(i)
if (k .lt. idz) then
    beta=beta+sum*sum
    sum=-sum
else
    beta=beta-sum*sum
end if
do 40 i=1,npt
40 vlag(i)=vlag(i)+sum*zmat(i,k)
bsum=zero
dx=zero
ssq=zero
do 70 j=1,n
sum=zero
do 50 i=1,npt
50 sum=sum+w(i)*bmat(i,j)
bsum=bsum+sum*step(j)
jp=npt+j
do 60 k=1,n
60 sum=sum+bmat(jp,k)*step(k)
vlag(jp)=sum
bsum=bsum+sum*step(j)
dx=dx+step(j)*xpt(kopt,j)
70 ssq=ssq+step(j)**2
beta=dx*dx+ssq*(sp(kopt)+dx+dx+half*ssq)+beta-bsum
vlag(kopt)=vlag(kopt)+one
!
!     If KNEW is zero initially, then pick the index of the interpolation
!       point to be deleted, by maximizing the absolute value of the
!       denominator of the updating formula times a weighting factor.
!       
!
if (knew .eq. 0) then
    denmax=zero
    do 100 k=1,npt
    hdiag=zero
    do 80 j=1,nptm
    temp=one
    if (j .lt. idz) temp=-one
80     hdiag=hdiag+temp*zmat(k,j)**2
    denabs=dabs(beta*hdiag+vlag(k)**2)
    distsq=zero
    do 90 j=1,n
90     distsq=distsq+(xpt(k,j)-xpt(kopt,j))**2
    temp=denabs*distsq*distsq
    if (temp .gt. denmax) then
        denmax=temp
        knew=k
    end if
100     continue
end if
!
!     Apply the rotations that put zeros in the KNEW-th row of ZMAT.
!
jl=1
if (nptm .ge. 2) then
    do 120 j=2,nptm
    if (j .eq. idz) then
        jl=idz
    else if (zmat(knew,j) .ne. zero) then
        temp=dsqrt(zmat(knew,jl)**2+zmat(knew,j)**2)
        tempa=zmat(knew,jl)/temp
        tempb=zmat(knew,j)/temp
        do 110 i=1,npt
        temp=tempa*zmat(i,jl)+tempb*zmat(i,j)
        zmat(i,j)=tempa*zmat(i,j)-tempb*zmat(i,jl)
110         zmat(i,jl)=temp
        zmat(knew,j)=zero
    end if
120     continue
end if
!
!     Put the first NPT components of the KNEW-th column of the Z Z^T matrix
!       into W, and calculate the parameters of the updating formula.
!
tempa=zmat(knew,1)
if (idz .ge. 2) tempa=-tempa
if (jl .gt. 1) tempb=zmat(knew,jl)
do 130 i=1,npt
w(i)=tempa*zmat(i,1)
if (jl .gt. 1) w(i)=w(i)+tempb*zmat(i,jl)
130 continue
alpha=w(knew)
tau=vlag(knew)
tausq=tau*tau
denom=alpha*beta+tausq
vlag(knew)=vlag(knew)-one
if (denom .eq. zero) then
    knew=0
    goto 180
end if
sqrtdn=dsqrt(dabs(denom))
!
!     Complete the updating of ZMAT when there is only one nonzero element
!       in the KNEW-th row of the new matrix ZMAT. IFLAG is set to one when
!       the value of IDZ is going to be reduced.
!
iflag=0
if (jl .eq. 1) then
    tempa=tau/sqrtdn
    tempb=zmat(knew,1)/sqrtdn
    do 140 i=1,npt
140     zmat(i,1)=tempa*zmat(i,1)-tempb*vlag(i)
    if (denom .lt. zero) then
        if (idz .eq. 1) then
            idz=2
        else
            iflag=1
        end if
    end if
else
!
!     Complete the updating of ZMAT in the alternative case.
!
    ja=1
    if (beta .ge. zero) ja=jl
    jb=jl+1-ja
    temp=zmat(knew,jb)/denom
    tempa=temp*beta
    tempb=temp*tau
    temp=zmat(knew,ja)
    scala=one/dsqrt(dabs(beta)*temp*temp+tausq)
    scalb=scala*sqrtdn
    do 150 i=1,npt
    zmat(i,ja)=scala*(tau*zmat(i,ja)-temp*vlag(i))
150     zmat(i,jb)=scalb*(zmat(i,jb)-tempa*w(i)-tempb*vlag(i))
    if (denom .le. zero) then
        if (beta .lt. zero) then
            idz=idz+1
        else
            iflag=1
        end if
    end if
end if
!
!     Reduce IDZ when the diagonal part of the ZMAT times Diag(DZ) times
!       ZMAT^T factorization gains another positive element. Then exchange
!       the first and IDZ-th columns of ZMAT.
!
if (iflag .eq. 1) then
    idz=idz-1
    do 160 i=1,npt
    temp=zmat(i,1)
    zmat(i,1)=zmat(i,idz)
160     zmat(i,idz)=temp
end if
!
!     Finally, update the matrix BMAT.
!
do 170 j=1,n
jp=npt+j
w(jp)=bmat(knew,j)
tempa=(alpha*vlag(jp)-tau*w(jp))/denom
tempb=(-beta*w(jp)-tau*vlag(jp))/denom
do 170 i=1,jp
bmat(i,j)=bmat(i,j)+tempa*vlag(i)+tempb*w(i)
if (i .gt. npt) bmat(jp,i-npt)=bmat(i,j)
170 continue
180 return
end
