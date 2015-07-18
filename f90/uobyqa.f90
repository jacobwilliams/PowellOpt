subroutine calfun (n,x,f)
implicit real*8 (a-h,o-z)
dimension x(*),y(10,10)
do 10 j=1,n
y(1,j)=1.0d0
10 y(2,j)=2.0d0*x(j)-1.0d0
do 20 i=2,n
do 20 j=1,n
20 y(i+1,j)=2.0d0*y(2,j)*y(i,j)-y(i-1,j)
f=0.0d0
np=n+1
iw=1
do 40 i=1,np
sum=0.0d0
do 30 j=1,n
30 sum=sum+y(i,j)
sum=sum/dfloat(n)
if (iw .gt. 0) sum=sum+1.0/dfloat(i*i-2*i)
iw=-iw
40 f=f+sum*sum
return
end
subroutine lagmax (n,g,h,rho,d,v,vmax)
implicit real*8 (a-h,o-z)
dimension g(*),h(n,*),d(*),v(*)
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
half=0.5d0
halfrt=dsqrt(half)
one=1.0d0
zero=0.0d0
!
!     Pick V such that ||HV|| / ||V|| is large.
!
hmax=zero
do 20 i=1,n
sum=zero
do 10 j=1,n
h(j,i)=h(i,j)
10 sum=sum+h(i,j)**2
if (sum .gt. hmax) then
    hmax=sum
    k=i
end if
20 continue
do 30 j=1,n
30 v(j)=h(k,j)
!
!     Set D to a vector in the subspace spanned by V and HV that maximizes
!     |(D,HD)|/(D,D), except that we set D=HV if V and HV are nearly parallel.
!     The vector that has the name D at label 60 used to be the vector W.
!
vsq=zero
vhv=zero
dsq=zero
do 50 i=1,n
vsq=vsq+v(i)**2
d(i)=zero
do 40 j=1,n
40 d(i)=d(i)+h(i,j)*v(j)
vhv=vhv+v(i)*d(i)
50 dsq=dsq+d(i)**2
if (vhv*vhv .le. 0.9999d0*dsq*vsq) then
    temp=vhv/vsq
    wsq=zero
    do 60 i=1,n
    d(i)=d(i)-temp*v(i)
60     wsq=wsq+d(i)**2
    whw=zero
    ratio=dsqrt(wsq/vsq)
    do 80 i=1,n
    temp=zero
    do 70 j=1,n
70     temp=temp+h(i,j)*d(j)
    whw=whw+temp*d(i)
80     v(i)=ratio*v(i)
    vhv=ratio*ratio*vhv
    vhw=ratio*wsq
    temp=half*(whw-vhv)
    temp=temp+dsign(dsqrt(temp**2+vhw**2),whw+vhv)
    do 90 i=1,n
90     d(i)=vhw*v(i)+temp*d(i)
end if
!
!     We now turn our attention to the subspace spanned by G and D. A multiple
!     of the current D is returned if that choice seems to be adequate.
!
gg=zero
gd=zero
dd=zero
dhd=zero
do 110 i=1,n
gg=gg+g(i)**2
gd=gd+g(i)*d(i)
dd=dd+d(i)**2
sum=zero
do 100 j=1,n
100 sum=sum+h(i,j)*d(j)
110 dhd=dhd+sum*d(i)
temp=gd/gg
vv=zero
scale=dsign(rho/dsqrt(dd),gd*dhd)
do 120 i=1,n
v(i)=d(i)-temp*g(i)
vv=vv+v(i)**2
120 d(i)=scale*d(i)
gnorm=dsqrt(gg)
if (gnorm*dd .le. 0.5d-2*rho*dabs(dhd) .or. &
  vv/dd .le. 1.0d-4) then
    vmax=dabs(scale*(gd+half*scale*dhd))
    goto 170
end if
!
!     G and V are now orthogonal in the subspace spanned by G and D. Hence
!     we generate an orthonormal basis of this subspace such that (D,HV) is
!     negligible or zero, where D and V will be the basis vectors.
!
ghg=zero
vhg=zero
vhv=zero
do 140 i=1,n
sum=zero
sumv=zero
do 130 j=1,n
sum=sum+h(i,j)*g(j)
130 sumv=sumv+h(i,j)*v(j)
ghg=ghg+sum*g(i)
vhg=vhg+sumv*g(i)
140 vhv=vhv+sumv*v(i)
vnorm=dsqrt(vv)
ghg=ghg/gg
vhg=vhg/(vnorm*gnorm)
vhv=vhv/vv
if (dabs(vhg) .le. 0.01d0*dmax1(dabs(ghg),dabs(vhv))) then
    vmu=ghg-vhv
    wcos=one
    wsin=zero
else
    temp=half*(ghg-vhv)
    vmu=temp+dsign(dsqrt(temp**2+vhg**2),temp)
    temp=dsqrt(vmu**2+vhg**2)
    wcos=vmu/temp
    wsin=vhg/temp
end if
tempa=wcos/gnorm
tempb=wsin/vnorm
tempc=wcos/vnorm
tempd=wsin/gnorm
do 150 i=1,n
d(i)=tempa*g(i)+tempb*v(i)
150 v(i)=tempc*v(i)-tempd*g(i)
!
!     The final D is a multiple of the current D, V, D+V or D-V. We make the
!     choice from these possibilities that is optimal.
!
dlin=wcos*gnorm/rho
vlin=-wsin*gnorm/rho
tempa=dabs(dlin)+half*dabs(vmu+vhv)
tempb=dabs(vlin)+half*dabs(ghg-vmu)
tempc=halfrt*(dabs(dlin)+dabs(vlin))+0.25d0*dabs(ghg+vhv)
if (tempa .ge. tempb .and. tempa .ge. tempc) then
    tempd=dsign(rho,dlin*(vmu+vhv))
    tempv=zero
else if (tempb .ge. tempc) then
    tempd=zero
    tempv=dsign(rho,vlin*(ghg-vmu))
else
    tempd=dsign(halfrt*rho,dlin*(ghg+vhv))
    tempv=dsign(halfrt*rho,vlin*(ghg+vhv))
end if
do 160 i=1,n
160 d(i)=tempd*d(i)+tempv*v(i)
vmax=rho*rho*dmax1(tempa,tempb,tempc)
170 return
end
!
!     The Chebyquad test problem (Fletcher, 1965) for N = 2,4,6,8.
!
implicit real*8 (a-h,o-z)
dimension x(10),w(10000)
iprint=2
maxfun=5000
rhoend=1.0d-8
do 30 n=2,8,2
do 10 i=1,n
10 x(i)=dfloat(i)/dfloat(n+1)
rhobeg=0.2d0*x(1)
print 20, n
20 format (//5x,'******************'/5x, &
  'Results with N =',i2,/5x,'******************')
call uobyqa (n,x,rhobeg,rhoend,iprint,maxfun,w)
30 continue
stop
end
subroutine trstep (n,g,h,delta,tol,d,gg,td,tn,w,piv,z,evalue)
implicit real*8 (a-h,o-z)
dimension g(*),h(n,*),d(*),gg(*),td(*),tn(*),w(*),piv(*),z(*)
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
one=1.0d0
two=2.0d0
zero=0.0d0
delsq=delta*delta
evalue=zero
nm=n-1
do 10 i=1,n
d(i)=zero
td(i)=h(i,i)
do 10 j=1,i
10 h(i,j)=h(j,i)
!
!     Apply Householder transformations to obtain a tridiagonal matrix that
!     is similar to H, and put the elements of the Householder vectors in
!     the lower triangular part of H. Further, TD and TN will contain the
!     diagonal and other nonzero elements of the tridiagonal matrix.
!
do 80 k=1,nm
kp=k+1
sum=zero
if (kp .lt. n) then
    kpp=kp+1
    do 20 i=kpp,n
20     sum=sum+h(i,k)**2
end if
if (sum .eq. zero) then
    tn(k)=h(kp,k)
    h(kp,k)=zero
else
    temp=h(kp,k)
    tn(k)=dsign(dsqrt(sum+temp*temp),temp)
    h(kp,k)=-sum/(temp+tn(k))
    temp=dsqrt(two/(sum+h(kp,k)**2))
    do 30 i=kp,n
    w(i)=temp*h(i,k)
    h(i,k)=w(i)
30     z(i)=td(i)*w(i)
    wz=zero
    do 50 j=kp,nm
    jp=j+1
    do 40 i=jp,n
    z(i)=z(i)+h(i,j)*w(j)
40     z(j)=z(j)+h(i,j)*w(i)
50     wz=wz+w(j)*z(j)
    wz=wz+w(n)*z(n)
    do 70 j=kp,n
    td(j)=td(j)+w(j)*(wz*w(j)-two*z(j))
    if (j .lt. n) then
        jp=j+1
        do 60 i=jp,n
60         h(i,j)=h(i,j)-w(i)*z(j)-w(j)*(z(i)-wz*w(i))
    end if
70     continue
end if
80 continue
!
!     Form GG by applying the similarity transformation to G.
!
gsq=zero
do 90 i=1,n
gg(i)=g(i)
90 gsq=gsq+g(i)**2
gnorm=dsqrt(gsq)
do 110 k=1,nm
kp=k+1
sum=zero
do 100 i=kp,n
100 sum=sum+gg(i)*h(i,k)
do 110 i=kp,n
110 gg(i)=gg(i)-sum*h(i,k)
!
!     Begin the trust region calculation with a tridiagonal matrix by
!     calculating the norm of H. Then treat the case when H is zero.
!
hnorm=dabs(td(1))+dabs(tn(1))
tdmin=td(1)
tn(n)=zero
do 120 i=2,n
temp=dabs(tn(i-1))+dabs(td(i))+dabs(tn(i))
hnorm=dmax1(hnorm,temp)
120 tdmin=dmin1(tdmin,td(i))
if (hnorm .eq. zero) then
    if (gnorm .eq. zero) goto 400
    scale=delta/gnorm
    do 130 i=1,n
130     d(i)=-scale*gg(i)
    goto 370
end if
!
!     Set the initial values of PAR and its bounds.
!
parl=dmax1(zero,-tdmin,gnorm/delta-hnorm)
parlest=parl
par=parl
paru=zero
paruest=zero
posdef=zero
iterc=0
!
!     Calculate the pivots of the Cholesky factorization of (H+PAR*I).
!
140 iterc=iterc+1
ksav=0
piv(1)=td(1)+par
k=1
150 if (piv(k) .gt. zero) then
    piv(k+1)=td(k+1)+par-tn(k)**2/piv(k)
else
    if (piv(k) .lt. zero .or. tn(k) .ne. zero) goto 160
    ksav=k
    piv(k+1)=td(k+1)+par
end if
k=k+1
if (k .lt. n) goto 150
if (piv(k) .lt. zero) goto 160
if (piv(k) .eq. zero) ksav=k
!
!     Branch if all the pivots are positive, allowing for the case when
!     G is zero.
!
if (ksav .eq. 0 .and. gsq .gt. zero) goto 230
if (gsq .eq. zero) then
    if (par .eq. zero) goto 370
    paru=par
    paruest=par
    if (ksav .eq. 0) goto 190
end if
k=ksav
!
!     Set D to a direction of nonpositive curvature of the given tridiagonal
!     matrix, and thus revise PARLEST.
!
160 d(k)=one
if (dabs(tn(k)) .le. dabs(piv(k))) then
    dsq=one
    dhd=piv(k)
else
    temp=td(k+1)+par
    if (temp .le. dabs(piv(k))) then
        d(k+1)=dsign(one,-tn(k))
        dhd=piv(k)+temp-two*dabs(tn(k))
    else
        d(k+1)=-tn(k)/temp
        dhd=piv(k)+tn(k)*d(k+1)
    end if
    dsq=one+d(k+1)**2
end if
170 if (k .gt. 1) then
    k=k-1
    if (tn(k) .ne. zero) then
        d(k)=-tn(k)*d(k+1)/piv(k)
        dsq=dsq+d(k)**2
        goto 170
    end if
    do 180 i=1,k
180     d(i)=zero
end if
parl=par
parlest=par-dhd/dsq
!
!     Terminate with D set to a multiple of the current D if the following
!     test suggests that it suitable to do so.
!
190 temp=paruest
if (gsq .eq. zero) temp=temp*(one-tol)
if (paruest .gt. zero .and. parlest .ge. temp) then
    dtg=zero
    do 200 i=1,n
200     dtg=dtg+d(i)*gg(i)
    scale=-dsign(delta/dsqrt(dsq),dtg)
    do 210 i=1,n
210     d(i)=scale*d(i)
    goto 370
end if
!
!     Pick the value of PAR for the next iteration.
!
220 if (paru .eq. zero) then
    par=two*parlest+gnorm/delta
else
    par=0.5d0*(parl+paru)
    par=dmax1(par,parlest)
end if
if (paruest .gt. zero) par=dmin1(par,paruest)
goto 140
!
!     Calculate D for the current PAR in the positive definite case.
!
230 w(1)=-gg(1)/piv(1)
do 240 i=2,n
240 w(i)=(-gg(i)-tn(i-1)*w(i-1))/piv(i)
d(n)=w(n)
do 250 i=nm,1,-1
250 d(i)=w(i)-tn(i)*d(i+1)/piv(i)
!
!     Branch if a Newton-Raphson step is acceptable.
!
dsq=zero
wsq=zero
do 260 i=1,n
dsq=dsq+d(i)**2
260 wsq=wsq+piv(i)*w(i)**2
if (par .eq. zero .and. dsq .le. delsq) goto 320
!
!     Make the usual test for acceptability of a full trust region step.
!
dnorm=dsqrt(dsq)
phi=one/dnorm-one/delta
temp=tol*(one+par*dsq/wsq)-dsq*phi*phi
if (temp .ge. zero) then
    scale=delta/dnorm
    do 270 i=1,n
270     d(i)=scale*d(i)
    goto 370
end if
if (iterc .ge. 2 .and. par .le. parl) goto 370
if (paru .gt. zero .and. par .ge. paru) goto 370
!
!     Complete the iteration when PHI is negative.
!
if (phi .lt. zero) then
    parlest=par
    if (posdef .eq. one) then
        if (phi .le. phil) goto 370
        slope=(phi-phil)/(par-parl)
        parlest=par-phi/slope
    end if
    slope=one/gnorm
    if (paru .gt. zero) slope=(phiu-phi)/(paru-par)
    temp=par-phi/slope
    if (paruest .gt. zero) temp=dmin1(temp,paruest)
    paruest=temp
    posdef=one
    parl=par
    phil=phi
    goto 220
end if
!
!     If required, calculate Z for the alternative test for convergence.
!
if (posdef .eq. zero) then
    w(1)=one/piv(1)
    do 280 i=2,n
    temp=-tn(i-1)*w(i-1)
280     w(i)=(dsign(one,temp)+temp)/piv(i)
    z(n)=w(n)
    do 290 i=nm,1,-1
290     z(i)=w(i)-tn(i)*z(i+1)/piv(i)
    wwsq=zero
    zsq=zero
    dtz=zero
    do 300 i=1,n
    wwsq=wwsq+piv(i)*w(i)**2
    zsq=zsq+z(i)**2
300     dtz=dtz+d(i)*z(i)
!
!     Apply the alternative test for convergence.
!
    tempa=dabs(delsq-dsq)
    tempb=dsqrt(dtz*dtz+tempa*zsq)
    gam=tempa/(dsign(tempb,dtz)+dtz)
    temp=tol*(wsq+par*delsq)-gam*gam*wwsq
    if (temp .ge. zero) then
        do 310 i=1,n
310         d(i)=d(i)+gam*z(i)
        goto 370
    end if
    parlest=dmax1(parlest,par-wwsq/zsq)
end if
!
!     Complete the iteration when PHI is positive.
!
slope=one/gnorm
if (paru .gt. zero) then
    if (phi .ge. phiu) goto 370
    slope=(phiu-phi)/(paru-par)
end if
parlest=dmax1(parlest,par-phi/slope)
paruest=par
if (posdef .eq. one) then
    slope=(phi-phil)/(par-parl)
    paruest=par-phi/slope
end if
paru=par
phiu=phi
goto 220
!
!     Set EVALUE to the least eigenvalue of the second derivative matrix if
!     D is a Newton-Raphson step. SHFMAX will be an upper bound on EVALUE.
!
320 shfmin=zero
pivot=td(1)
shfmax=pivot
do 330 k=2,n
pivot=td(k)-tn(k-1)**2/pivot
330 shfmax=dmin1(shfmax,pivot)
!
!     Find EVALUE by a bisection method, but occasionally SHFMAX may be
!     adjusted by the rule of false position.
!
ksave=0
340 shift=0.5d0*(shfmin+shfmax)
k=1
temp=td(1)-shift
350 if (temp .gt. zero) then
    piv(k)=temp
    if (k .lt. n) then
        temp=td(k+1)-shift-tn(k)**2/temp
        k=k+1
        goto 350
    end if
    shfmin=shift
else
    if (k .lt. ksave) goto 360
    if (k .eq. ksave) then
        if (pivksv .eq. zero) goto 360
        if (piv(k)-temp .lt. temp-pivksv) then
            pivksv=temp
            shfmax=shift
        else
            pivksv=zero
            shfmax=(shift*piv(k)-shfmin*temp)/(piv(k)-temp)
        end if
    else
        ksave=k
        pivksv=temp
        shfmax=shift
    end if
end if
if (shfmin .le. 0.99d0*shfmax) goto 340
360 evalue=shfmin
!
!     Apply the inverse Householder transformations to D.
!
370 nm=n-1
do 390 k=nm,1,-1
kp=k+1
sum=zero
do 380 i=kp,n
380 sum=sum+d(i)*h(i,k)
do 390 i=kp,n
390 d(i)=d(i)-sum*h(i,k)
!
!     Return from the subroutine.
!
400 return
end
subroutine uobyqa (n,x,rhobeg,rhoend,iprint,maxfun,w)
implicit real*8 (a-h,o-z)
dimension x(*),w(*)
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
npt=(n*n+3*n+2)/2
ixb=1
ixo=ixb+n
ixn=ixo+n
ixp=ixn+n
ipq=ixp+n*npt
ipl=ipq+npt-1
ih=ipl+(npt-1)*npt
ig=ih+n*n
id=ig+n
ivl=ih
iw=id+n
call uobyqb (n,x,rhobeg,rhoend,iprint,maxfun,npt,w(ixb),w(ixo), &
  w(ixn),w(ixp),w(ipq),w(ipl),w(ih),w(ig),w(id),w(ivl),w(iw))
return
end
subroutine uobyqb (n,x,rhobeg,rhoend,iprint,maxfun,npt,xbase, &
  xopt,xnew,xpt,pq,pl,h,g,d,vlag,w)
implicit real*8 (a-h,o-z)
dimension x(*),xbase(*),xopt(*),xnew(*),xpt(npt,*),pq(*), &
  pl(npt,*),h(n,*),g(*),d(*),vlag(*),w(*)
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
one=1.0d0
two=2.0d0
zero=0.0d0
half=0.5d0
tol=0.01d0
nnp=n+n+1
nptm=npt-1
nftest=max0(maxfun,1)
!
!     Initialization. NF is the number of function calculations so far.
!
rho=rhobeg
rhosq=rho*rho
nf=0
do 10 i=1,n
xbase(i)=x(i)
do 10 k=1,npt
10 xpt(k,i)=zero
do 20 k=1,npt
do 20 j=1,nptm
20 pl(k,j)=zero
!
!     The branch to label 120 obtains a new value of the objective function
!     and then there is a branch back to label 50, because the new function
!     value is needed to form the initial quadratic model. The least function
!     value so far and its index are noted below.
!
30 do 40 i=1,n
40 x(i)=xbase(i)+xpt(nf+1,i)
goto 120
50 if (nf .eq. 1) then
    fopt=f
    kopt=nf
    fbase=f
    j=0
    jswitch=-1
    ih=n
else
    if (f .lt. fopt) then
        fopt=f
        kopt=nf
    end if
end if
!
!     Form the gradient and diagonal second derivatives of the initial
!     quadratic model and Lagrange functions.
!
if (nf .le. nnp) then
    jswitch=-jswitch
    if (jswitch .gt. 0) then
        if (j .ge. 1) then
            ih=ih+j
            if (w(j) .lt. zero) then
                d(j)=(fsave+f-two*fbase)/rhosq
                pq(j)=(fsave-f)/(two*rho)
                pl(1,ih)=-two/rhosq
                pl(nf-1,j)=half/rho
                pl(nf-1,ih)=one/rhosq
            else
                pq(j)=(4.0d0*fsave-3.0d0*fbase-f)/(two*rho)
                d(j)=(fbase+f-two*fsave)/rhosq
                pl(1,j)=-1.5d0/rho
                pl(1,ih)=one/rhosq
                pl(nf-1,j)=two/rho
                pl(nf-1,ih)=-two/rhosq
            end if
            pq(ih)=d(j)
            pl(nf,j)=-half/rho
            pl(nf,ih)=one/rhosq
        end if
!
!     Pick the shift from XBASE to the next initial interpolation point
!     that provides diagonal second derivatives.
!
        if (j .lt. n) then
            j=j+1
            xpt(nf+1,j)=rho
        end if
    else
        fsave=f
        if (f .lt. fbase) then
            w(j)=rho
            xpt(nf+1,j)=two*rho
        else
            w(j)=-rho
            xpt(nf+1,j)=-rho
        end if
    end if
    if (nf .lt. nnp) goto 30
!
!     Form the off-diagonal second derivatives of the initial quadratic model.
!
    ih=n
    ip=1
    iq=2
end if
ih=ih+1
if (nf .gt. nnp) then
    temp=one/(w(ip)*w(iq))
    tempa=f-fbase-w(ip)*pq(ip)-w(iq)*pq(iq)
    pq(ih)=(tempa-half*rhosq*(d(ip)+d(iq)))*temp
    pl(1,ih)=temp
    iw=ip+ip
    if (w(ip) .lt. zero) iw=iw+1
    pl(iw,ih)=-temp
    iw=iq+iq
    if (w(iq) .lt. zero) iw=iw+1
    pl(iw,ih)=-temp
    pl(nf,ih)=temp
!
!     Pick the shift from XBASE to the next initial interpolation point
!     that provides off-diagonal second derivatives.
!
    ip=ip+1
end if
if (ip .eq. iq) then
    ih=ih+1
    ip=1
    iq=iq+1
end if
if (nf .lt. npt) then
    xpt(nf+1,ip)=w(ip)
    xpt(nf+1,iq)=w(iq)
    goto 30
end if
!
!     Set parameters to begin the iterations for the current RHO.
!
sixthm=zero
delta=rho
60 tworsq=(two*rho)**2
rhosq=rho*rho
!
!     Form the gradient of the quadratic model at the trust region centre.
!
70 knew=0
ih=n
do 80 j=1,n
xopt(j)=xpt(kopt,j)
g(j)=pq(j)
do 80 i=1,j
ih=ih+1
g(i)=g(i)+pq(ih)*xopt(j)
if (i .lt. j) g(j)=g(j)+pq(ih)*xopt(i)
80 h(i,j)=pq(ih)
!
!     Generate the next trust region step and test its length. Set KNEW
!     to -1 if the purpose of the next F will be to improve conditioning,
!     and also calculate a lower bound on the Hessian term of the model Q.
!
call trstep (n,g,h,delta,tol,d,w(1),w(n+1),w(2*n+1),w(3*n+1), &
  w(4*n+1),w(5*n+1),evalue)
temp=zero
do 90 i=1,n
90 temp=temp+d(i)**2
dnorm=dmin1(delta,dsqrt(temp))
errtol=-one
if (dnorm .lt. half*rho) then
    knew=-1
    errtol=half*evalue*rho*rho
    if (nf .le. npt+9) errtol=zero
    goto 290
end if
!
!     Calculate the next value of the objective function.
!
100 do 110 i=1,n
xnew(i)=xopt(i)+d(i)
110 x(i)=xbase(i)+xnew(i)
120 if (nf .ge. nftest) then
    if (iprint .gt. 0) print 130
130     format (/4x,'Return from UOBYQA because CALFUN has been', &
      ' called MAXFUN times')
    goto 420
end if
nf=nf+1
call calfun (n,x,f)
if (iprint .eq. 3) then
    print 140, nf,f,(x(i),i=1,n)
140      format (/4x,'Function number',i6,'    F =',1pd18.10, &
       '    The corresponding X is:'/(2x,5d15.6))
end if
if (nf .le. npt) goto 50
if (knew .eq. -1) goto 420
!
!     Use the quadratic model to predict the change in F due to the step D,
!     and find the values of the Lagrange functions at the new point.
!
vquad=zero
ih=n
do 150 j=1,n
w(j)=d(j)
vquad=vquad+w(j)*pq(j)
do 150 i=1,j
ih=ih+1
w(ih)=d(i)*xnew(j)+d(j)*xopt(i)
if (i .eq. j) w(ih)=half*w(ih)
150 vquad=vquad+w(ih)*pq(ih)
do 170 k=1,npt
temp=zero
do 160 j=1,nptm
160 temp=temp+w(j)*pl(k,j)
170 vlag(k)=temp
vlag(kopt)=vlag(kopt)+one
!
!     Update SIXTHM, which is a lower bound on one sixth of the greatest
!     third derivative of F.
!
diff=f-fopt-vquad
sum=zero
do 190 k=1,npt
temp=zero
do 180 i=1,n
180 temp=temp+(xpt(k,i)-xnew(i))**2
temp=dsqrt(temp)
190 sum=sum+dabs(temp*temp*temp*vlag(k))
sixthm=dmax1(sixthm,dabs(diff)/sum)
!
!     Update FOPT and XOPT if the new F is the least value of the objective
!     function so far. Then branch if D is not a trust region step.
!
fsave=fopt
if (f .lt. fopt) then
    fopt=f
    do 200 i=1,n
200     xopt(i)=xnew(i)
end if
ksave=knew
if (knew .gt. 0) goto 240
!
!     Pick the next value of DELTA after a trust region step.
!
if (vquad .ge. zero) then
    if (iprint .gt. 0) print 210
210     format (/4x,'Return from UOBYQA because a trust', &
      ' region step has failed to reduce Q')
    goto 420
end if
ratio=(f-fsave)/vquad
if (ratio .le. 0.1d0) then
    delta=half*dnorm
else if (ratio .le. 0.7d0) then
    delta=dmax1(half*delta,dnorm)
else
    delta=dmax1(delta,1.25d0*dnorm,dnorm+rho)
end if
if (delta .le. 1.5d0*rho) delta=rho
!
!     Set KNEW to the index of the next interpolation point to be deleted.
!
ktemp=0
detrat=zero
if (f .ge. fsave) then
    ktemp=kopt
    detrat=one
end if
do 230 k=1,npt
sum=zero
do 220 i=1,n
220 sum=sum+(xpt(k,i)-xopt(i))**2
temp=dabs(vlag(k))
if (sum .gt. rhosq) temp=temp*(sum/rhosq)**1.5d0
if (temp .gt. detrat .and. k .ne. ktemp) then
    detrat=temp
    ddknew=sum
    knew=k
end if
230 continue
if (knew .eq. 0) goto 290
!
!     Replace the interpolation point that has index KNEW by the point XNEW,
!     and also update the Lagrange functions and the quadratic model.
!
240 do 250 i=1,n
250 xpt(knew,i)=xnew(i)
temp=one/vlag(knew)
do 260 j=1,nptm
pl(knew,j)=temp*pl(knew,j)
260 pq(j)=pq(j)+diff*pl(knew,j)
do 280 k=1,npt
if (k .ne. knew) then
    temp=vlag(k)
    do 270 j=1,nptm
270     pl(k,j)=pl(k,j)-temp*pl(knew,j)
end if
280 continue
!
!     Update KOPT if F is the least calculated value of the objective
!     function. Then branch for another trust region calculation. The
!     case KSAVE>0 indicates that a model step has just been taken.
!
if (f .lt. fsave) then
    kopt=knew
    goto 70
end if
if (ksave .gt. 0) goto 70
if (dnorm .gt. two*rho) goto 70
if (ddknew .gt. tworsq) goto 70
!
!     Alternatively, find out if the interpolation points are close
!     enough to the best point so far.
!
290 do 300 k=1,npt
w(k)=zero
do 300 i=1,n
300 w(k)=w(k)+(xpt(k,i)-xopt(i))**2
310 knew=-1
distest=tworsq
do 320 k=1,npt
if (w(k) .gt. distest) then
    knew=k
    distest=w(k)
end if
320 continue
!
!     If a point is sufficiently far away, then set the gradient and Hessian
!     of its Lagrange function at the centre of the trust region, and find
!     half the sum of squares of components of the Hessian.
!
if (knew .gt. 0) then
    ih=n
    sumh=zero
    do 340 j=1,n
    g(j)=pl(knew,j)
    do 330 i=1,j
    ih=ih+1
    temp=pl(knew,ih)
    g(j)=g(j)+temp*xopt(i)
    if (i .lt. j) then
        g(i)=g(i)+temp*xopt(j)
        sumh=sumh+temp*temp
    end if
330     h(i,j)=temp
340     sumh=sumh+half*temp*temp
!
!     If ERRTOL is positive, test whether to replace the interpolation point
!     with index KNEW, using a bound on the maximum modulus of its Lagrange
!     function in the trust region.
!
    if (errtol .gt. zero) then
        w(knew)=zero
        sumg=zero
        do 350 i=1,n
350         sumg=sumg+g(i)**2
        estim=rho*(dsqrt(sumg)+rho*dsqrt(half*sumh))
        wmult=sixthm*distest**1.5d0
        if (wmult*estim .le. errtol) goto 310
    end if
!
!     If the KNEW-th point may be replaced, then pick a D that gives a large
!     value of the modulus of its Lagrange function within the trust region.
!     Here the vector XNEW is used as temporary working space.
!
    call lagmax (n,g,h,rho,d,xnew,vmax)
    if (errtol .gt. zero) then
        if (wmult*vmax .le. errtol) goto 310
    end if
    goto 100
end if
if (dnorm .gt. rho) goto 70
!
!     Prepare to reduce RHO by shifting XBASE to the best point so far,
!     and make the corresponding changes to the gradients of the Lagrange
!     functions and the quadratic model.
!
if (rho .gt. rhoend) then
    ih=n
    do 380 j=1,n
    xbase(j)=xbase(j)+xopt(j)
    do 360 k=1,npt
360     xpt(k,j)=xpt(k,j)-xopt(j)
    do 380 i=1,j
    ih=ih+1
    pq(i)=pq(i)+pq(ih)*xopt(j)
    if (i .lt. j) then
        pq(j)=pq(j)+pq(ih)*xopt(i)
        do 370 k=1,npt
370         pl(k,j)=pl(k,j)+pl(k,ih)*xopt(i)
    end if
    do 380 k=1,npt
380     pl(k,i)=pl(k,i)+pl(k,ih)*xopt(j)
!
!     Pick the next values of RHO and DELTA.
!
    delta=half*rho
    ratio=rho/rhoend
    if (ratio .le. 16.0d0) then
        rho=rhoend
    else if (ratio .le. 250.0d0) then
        rho=dsqrt(ratio)*rhoend
    else
        rho=0.1d0*rho
    end if
    delta=dmax1(delta,rho)
    if (iprint .ge. 2) then
        if (iprint .ge. 3) print 390
390         format (5x)
        print 400, rho,nf
400         format (/4x,'New RHO =',1pd11.4,5x,'Number of', &
          ' function values =',i6)
        print 410, fopt,(xbase(i),i=1,n)
410         format (4x,'Least value of F =',1pd23.15,9x, &
          'The corresponding X is:'/(2x,5d15.6))
    end if
    goto 60
end if
!
!     Return from the calculation, after another Newton-Raphson step, if
!     it is too short to have been tried before.
!
if (errtol .ge. zero) goto 100
420 if (fopt .le. f) then
    do 430 i=1,n
430     x(i)=xbase(i)+xopt(i)
    f=fopt
end if
if (iprint .ge. 1) then
    print 440, nf
440     format (/4x,'At the return from UOBYQA',5x, &
      'Number of function values =',i6)
    print 410, f,(x(i),i=1,n)
end if
return
end
