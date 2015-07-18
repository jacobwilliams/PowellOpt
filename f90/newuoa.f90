subroutine bigden (n,npt,xopt,xpt,bmat,zmat,idz,ndim,kopt, &
  knew,d,w,vlag,beta,s,wvec,prod)
implicit real*8 (a-h,o-z)
dimension xopt(*),xpt(npt,*),bmat(ndim,*),zmat(npt,*),d(*), &
  w(*),vlag(*),s(*),wvec(ndim,*),prod(ndim,*)
dimension den(9),denex(9),par(9)
!
!     N is the number of variables.
!     NPT is the number of interpolation equations.
!     XOPT is the best interpolation point so far.
!     XPT contains the coordinates of the current interpolation points.
!     BMAT provides the last N columns of H.
!     ZMAT and IDZ give a factorization of the first NPT by NPT submatrix of H.
!     NDIM is the first dimension of BMAT and has the value NPT+N.
!     KOPT is the index of the optimal interpolation point.
!     KNEW is the index of the interpolation point that is going to be moved.
!     D will be set to the step from XOPT to the new point, and on entry it
!       should be the D that was calculated by the last call of BIGLAG. The
!       length of the initial D provides a trust region bound on the final D.
!     W will be set to Wcheck for the final choice of D.
!     VLAG will be set to Theta*Wcheck+e_b for the final choice of D.
!     BETA will be set to the value that will occur in the updating formula
!       when the KNEW-th interpolation point is moved to its new position.
!     S, WVEC, PROD and the private arrays DEN, DENEX and PAR will be used
!       for working space.
!
!     D is calculated in a way that should provide a denominator with a large
!     modulus in the updating formula when the KNEW-th interpolation point is
!     shifted to the new position XOPT+D.
!
!     Set some constants.
!
half=0.5d0
one=1.0d0
quart=0.25d0
two=2.0d0
zero=0.0d0
twopi=8.0d0*datan(one)
nptm=npt-n-1
!
!     Store the first NPT elements of the KNEW-th column of H in W(N+1)
!     to W(N+NPT).
!
do 10 k=1,npt
10 w(n+k)=zero
do 20 j=1,nptm
temp=zmat(knew,j)
if (j .lt. idz) temp=-temp
do 20 k=1,npt
20 w(n+k)=w(n+k)+temp*zmat(k,j)
alpha=w(n+knew)
!
!     The initial search direction D is taken from the last call of BIGLAG,
!     and the initial S is set below, usually to the direction from X_OPT
!     to X_KNEW, but a different direction to an interpolation point may
!     be chosen, in order to prevent S from being nearly parallel to D.
!
dd=zero
ds=zero
ss=zero
xoptsq=zero
do 30 i=1,n
dd=dd+d(i)**2
s(i)=xpt(knew,i)-xopt(i)
ds=ds+d(i)*s(i)
ss=ss+s(i)**2
30 xoptsq=xoptsq+xopt(i)**2
if (ds*ds .gt. 0.99d0*dd*ss) then
    ksav=knew
    dtest=ds*ds/ss
    do 50 k=1,npt
    if (k .ne. kopt) then
        dstemp=zero
        sstemp=zero
        do 40 i=1,n
        diff=xpt(k,i)-xopt(i)
        dstemp=dstemp+d(i)*diff
40         sstemp=sstemp+diff*diff
        if (dstemp*dstemp/sstemp .lt. dtest) then
            ksav=k
            dtest=dstemp*dstemp/sstemp
            ds=dstemp
            ss=sstemp
        end if
    end if
50     continue
    do 60 i=1,n
60     s(i)=xpt(ksav,i)-xopt(i)
end if
ssden=dd*ss-ds*ds
iterc=0
densav=zero
!
!     Begin the iteration by overwriting S with a vector that has the
!     required length and direction.
!
70 iterc=iterc+1
temp=one/dsqrt(ssden)
xoptd=zero
xopts=zero
do 80 i=1,n
s(i)=temp*(dd*s(i)-ds*d(i))
xoptd=xoptd+xopt(i)*d(i)
80 xopts=xopts+xopt(i)*s(i)
!
!     Set the coefficients of the first two terms of BETA.
!
tempa=half*xoptd*xoptd
tempb=half*xopts*xopts
den(1)=dd*(xoptsq+half*dd)+tempa+tempb
den(2)=two*xoptd*dd
den(3)=two*xopts*dd
den(4)=tempa-tempb
den(5)=xoptd*xopts
do 90 i=6,9
90 den(i)=zero
!
!     Put the coefficients of Wcheck in WVEC.
!
do 110 k=1,npt
tempa=zero
tempb=zero
tempc=zero
do 100 i=1,n
tempa=tempa+xpt(k,i)*d(i)
tempb=tempb+xpt(k,i)*s(i)
100 tempc=tempc+xpt(k,i)*xopt(i)
wvec(k,1)=quart*(tempa*tempa+tempb*tempb)
wvec(k,2)=tempa*tempc
wvec(k,3)=tempb*tempc
wvec(k,4)=quart*(tempa*tempa-tempb*tempb)
110 wvec(k,5)=half*tempa*tempb
do 120 i=1,n
ip=i+npt
wvec(ip,1)=zero
wvec(ip,2)=d(i)
wvec(ip,3)=s(i)
wvec(ip,4)=zero
120 wvec(ip,5)=zero
!
!     Put the coefficents of THETA*Wcheck in PROD.
!
do 190 jc=1,5
nw=npt
if (jc .eq. 2 .or. jc .eq. 3) nw=ndim
do 130 k=1,npt
130 prod(k,jc)=zero
do 150 j=1,nptm
sum=zero
do 140 k=1,npt
140 sum=sum+zmat(k,j)*wvec(k,jc)
if (j .lt. idz) sum=-sum
do 150 k=1,npt
150 prod(k,jc)=prod(k,jc)+sum*zmat(k,j)
if (nw .eq. ndim) then
    do 170 k=1,npt
    sum=zero
    do 160 j=1,n
160     sum=sum+bmat(k,j)*wvec(npt+j,jc)
170     prod(k,jc)=prod(k,jc)+sum
end if
do 190 j=1,n
sum=zero
do 180 i=1,nw
180 sum=sum+bmat(i,j)*wvec(i,jc)
190 prod(npt+j,jc)=sum
!
!     Include in DEN the part of BETA that depends on THETA.
!
do 210 k=1,ndim
sum=zero
do 200 i=1,5
par(i)=half*prod(k,i)*wvec(k,i)
200 sum=sum+par(i)
den(1)=den(1)-par(1)-sum
tempa=prod(k,1)*wvec(k,2)+prod(k,2)*wvec(k,1)
tempb=prod(k,2)*wvec(k,4)+prod(k,4)*wvec(k,2)
tempc=prod(k,3)*wvec(k,5)+prod(k,5)*wvec(k,3)
den(2)=den(2)-tempa-half*(tempb+tempc)
den(6)=den(6)-half*(tempb-tempc)
tempa=prod(k,1)*wvec(k,3)+prod(k,3)*wvec(k,1)
tempb=prod(k,2)*wvec(k,5)+prod(k,5)*wvec(k,2)
tempc=prod(k,3)*wvec(k,4)+prod(k,4)*wvec(k,3)
den(3)=den(3)-tempa-half*(tempb-tempc)
den(7)=den(7)-half*(tempb+tempc)
tempa=prod(k,1)*wvec(k,4)+prod(k,4)*wvec(k,1)
den(4)=den(4)-tempa-par(2)+par(3)
tempa=prod(k,1)*wvec(k,5)+prod(k,5)*wvec(k,1)
tempb=prod(k,2)*wvec(k,3)+prod(k,3)*wvec(k,2)
den(5)=den(5)-tempa-half*tempb
den(8)=den(8)-par(4)+par(5)
tempa=prod(k,4)*wvec(k,5)+prod(k,5)*wvec(k,4)
210 den(9)=den(9)-half*tempa
!
!     Extend DEN so that it holds all the coefficients of DENOM.
!
sum=zero
do 220 i=1,5
par(i)=half*prod(knew,i)**2
220 sum=sum+par(i)
denex(1)=alpha*den(1)+par(1)+sum
tempa=two*prod(knew,1)*prod(knew,2)
tempb=prod(knew,2)*prod(knew,4)
tempc=prod(knew,3)*prod(knew,5)
denex(2)=alpha*den(2)+tempa+tempb+tempc
denex(6)=alpha*den(6)+tempb-tempc
tempa=two*prod(knew,1)*prod(knew,3)
tempb=prod(knew,2)*prod(knew,5)
tempc=prod(knew,3)*prod(knew,4)
denex(3)=alpha*den(3)+tempa+tempb-tempc
denex(7)=alpha*den(7)+tempb+tempc
tempa=two*prod(knew,1)*prod(knew,4)
denex(4)=alpha*den(4)+tempa+par(2)-par(3)
tempa=two*prod(knew,1)*prod(knew,5)
denex(5)=alpha*den(5)+tempa+prod(knew,2)*prod(knew,3)
denex(8)=alpha*den(8)+par(4)-par(5)
denex(9)=alpha*den(9)+prod(knew,4)*prod(knew,5)
!
!     Seek the value of the angle that maximizes the modulus of DENOM.
!
sum=denex(1)+denex(2)+denex(4)+denex(6)+denex(8)
denold=sum
denmax=sum
isave=0
iu=49
temp=twopi/dfloat(iu+1)
par(1)=one
do 250 i=1,iu
angle=dfloat(i)*temp
par(2)=dcos(angle)
par(3)=dsin(angle)
do 230 j=4,8,2
par(j)=par(2)*par(j-2)-par(3)*par(j-1)
230 par(j+1)=par(2)*par(j-1)+par(3)*par(j-2)
sumold=sum
sum=zero
do 240 j=1,9
240 sum=sum+denex(j)*par(j)
if (dabs(sum) .gt. dabs(denmax)) then
    denmax=sum
    isave=i
    tempa=sumold
else if (i .eq. isave+1) then
    tempb=sum
end if
250 continue
if (isave .eq. 0) tempa=sum
if (isave .eq. iu) tempb=denold
step=zero
if (tempa .ne. tempb) then
    tempa=tempa-denmax
    tempb=tempb-denmax
    step=half*(tempa-tempb)/(tempa+tempb)
end if
angle=temp*(dfloat(isave)+step)
!
!     Calculate the new parameters of the denominator, the new VLAG vector
!     and the new D. Then test for convergence.
!
par(2)=dcos(angle)
par(3)=dsin(angle)
do 260 j=4,8,2
par(j)=par(2)*par(j-2)-par(3)*par(j-1)
260 par(j+1)=par(2)*par(j-1)+par(3)*par(j-2)
beta=zero
denmax=zero
do 270 j=1,9
beta=beta+den(j)*par(j)
270 denmax=denmax+denex(j)*par(j)
do 280 k=1,ndim
vlag(k)=zero
do 280 j=1,5
280 vlag(k)=vlag(k)+prod(k,j)*par(j)
tau=vlag(knew)
dd=zero
tempa=zero
tempb=zero
do 290 i=1,n
d(i)=par(2)*d(i)+par(3)*s(i)
w(i)=xopt(i)+d(i)
dd=dd+d(i)**2
tempa=tempa+d(i)*w(i)
290 tempb=tempb+w(i)*w(i)
if (iterc .ge. n) goto 340
if (iterc .gt. 1) densav=dmax1(densav,denold)
if (dabs(denmax) .le. 1.1d0*dabs(densav)) goto 340
densav=denmax
!
!     Set S to half the gradient of the denominator with respect to D.
!     Then branch for the next iteration.
!
do 300 i=1,n
temp=tempa*xopt(i)+tempb*d(i)-vlag(npt+i)
300 s(i)=tau*bmat(knew,i)+alpha*temp
do 320 k=1,npt
sum=zero
do 310 j=1,n
310 sum=sum+xpt(k,j)*w(j)
temp=(tau*w(n+k)-alpha*vlag(k))*sum
do 320 i=1,n
320 s(i)=s(i)+temp*xpt(k,i)
ss=zero
ds=zero
do 330 i=1,n
ss=ss+s(i)**2
330 ds=ds+d(i)*s(i)
ssden=dd*ss-ds*ds
if (ssden .ge. 1.0d-8*dd*ss) goto 70
!
!     Set the vector W before the RETURN from the subroutine.
!
340 do 350 k=1,ndim
w(k)=zero
do 350 j=1,5
350 w(k)=w(k)+wvec(k,j)*par(j)
vlag(kopt)=vlag(kopt)+one
return
end
subroutine biglag (n,npt,xopt,xpt,bmat,zmat,idz,ndim,knew, &
  delta,d,alpha,hcol,gc,gd,s,w)
implicit real*8 (a-h,o-z)
dimension xopt(*),xpt(npt,*),bmat(ndim,*),zmat(npt,*),d(*), &
  hcol(*),gc(*),gd(*),s(*),w(*)
!
!     N is the number of variables.
!     NPT is the number of interpolation equations.
!     XOPT is the best interpolation point so far.
!     XPT contains the coordinates of the current interpolation points.
!     BMAT provides the last N columns of H.
!     ZMAT and IDZ give a factorization of the first NPT by NPT submatrix of H.
!     NDIM is the first dimension of BMAT and has the value NPT+N.
!     KNEW is the index of the interpolation point that is going to be moved.
!     DELTA is the current trust region bound.
!     D will be set to the step from XOPT to the new point.
!     ALPHA will be set to the KNEW-th diagonal element of the H matrix.
!     HCOL, GC, GD, S and W will be used for working space.
!
!     The step D is calculated in a way that attempts to maximize the modulus
!     of LFUNC(XOPT+D), subject to the bound ||D|| .LE. DELTA, where LFUNC is
!     the KNEW-th Lagrange function.
!
!     Set some constants.
!
half=0.5d0
one=1.0d0
zero=0.0d0
twopi=8.0d0*datan(one)
delsq=delta*delta
nptm=npt-n-1
!
!     Set the first NPT components of HCOL to the leading elements of the
!     KNEW-th column of H.
!
iterc=0
do 10 k=1,npt
10 hcol(k)=zero
do 20 j=1,nptm
temp=zmat(knew,j)
if (j .lt. idz) temp=-temp
do 20 k=1,npt
20 hcol(k)=hcol(k)+temp*zmat(k,j)
alpha=hcol(knew)
!
!     Set the unscaled initial direction D. Form the gradient of LFUNC at
!     XOPT, and multiply D by the second derivative matrix of LFUNC.
!
dd=zero
do 30 i=1,n
d(i)=xpt(knew,i)-xopt(i)
gc(i)=bmat(knew,i)
gd(i)=zero
30 dd=dd+d(i)**2
do 50 k=1,npt
temp=zero
sum=zero
do 40 j=1,n
temp=temp+xpt(k,j)*xopt(j)
40 sum=sum+xpt(k,j)*d(j)
temp=hcol(k)*temp
sum=hcol(k)*sum
do 50 i=1,n
gc(i)=gc(i)+temp*xpt(k,i)
50 gd(i)=gd(i)+sum*xpt(k,i)
!
!     Scale D and GD, with a sign change if required. Set S to another
!     vector in the initial two dimensional subspace.
!
gg=zero
sp=zero
dhd=zero
do 60 i=1,n
gg=gg+gc(i)**2
sp=sp+d(i)*gc(i)
60 dhd=dhd+d(i)*gd(i)
scale=delta/dsqrt(dd)
if (sp*dhd .lt. zero) scale=-scale
temp=zero
if (sp*sp .gt. 0.99d0*dd*gg) temp=one
tau=scale*(dabs(sp)+half*scale*dabs(dhd))
if (gg*delsq .lt. 0.01d0*tau*tau) temp=one
do 70 i=1,n
d(i)=scale*d(i)
gd(i)=scale*gd(i)
70 s(i)=gc(i)+temp*gd(i)
!
!     Begin the iteration by overwriting S with a vector that has the
!     required length and direction, except that termination occurs if
!     the given D and S are nearly parallel.
!
80 iterc=iterc+1
dd=zero
sp=zero
ss=zero
do 90 i=1,n
dd=dd+d(i)**2
sp=sp+d(i)*s(i)
90 ss=ss+s(i)**2
temp=dd*ss-sp*sp
if (temp .le. 1.0d-8*dd*ss) goto 160
denom=dsqrt(temp)
do 100 i=1,n
s(i)=(dd*s(i)-sp*d(i))/denom
100 w(i)=zero
!
!     Calculate the coefficients of the objective function on the circle,
!     beginning with the multiplication of S by the second derivative matrix.
!
do 120 k=1,npt
sum=zero
do 110 j=1,n
110 sum=sum+xpt(k,j)*s(j)
sum=hcol(k)*sum
do 120 i=1,n
120 w(i)=w(i)+sum*xpt(k,i)
cf1=zero
cf2=zero
cf3=zero
cf4=zero
cf5=zero
do 130 i=1,n
cf1=cf1+s(i)*w(i)
cf2=cf2+d(i)*gc(i)
cf3=cf3+s(i)*gc(i)
cf4=cf4+d(i)*gd(i)
130 cf5=cf5+s(i)*gd(i)
cf1=half*cf1
cf4=half*cf4-cf1
!
!     Seek the value of the angle that maximizes the modulus of TAU.
!
taubeg=cf1+cf2+cf4
taumax=taubeg
tauold=taubeg
isave=0
iu=49
temp=twopi/dfloat(iu+1)
do 140 i=1,iu
angle=dfloat(i)*temp
cth=dcos(angle)
sth=dsin(angle)
tau=cf1+(cf2+cf4*cth)*cth+(cf3+cf5*cth)*sth
if (dabs(tau) .gt. dabs(taumax)) then
    taumax=tau
    isave=i
    tempa=tauold
else if (i .eq. isave+1) then
    tempb=tau
end if
140 tauold=tau
if (isave .eq. 0) tempa=tau
if (isave .eq. iu) tempb=taubeg
step=zero
if (tempa .ne. tempb) then
    tempa=tempa-taumax
    tempb=tempb-taumax
    step=half*(tempa-tempb)/(tempa+tempb)
end if
angle=temp*(dfloat(isave)+step)
!
!     Calculate the new D and GD. Then test for convergence.
!
cth=dcos(angle)
sth=dsin(angle)
tau=cf1+(cf2+cf4*cth)*cth+(cf3+cf5*cth)*sth
do 150 i=1,n
d(i)=cth*d(i)+sth*s(i)
gd(i)=cth*gd(i)+sth*w(i)
150 s(i)=gc(i)+gd(i)
if (dabs(tau) .le. 1.1d0*dabs(taubeg)) goto 160
if (iterc .lt. n) goto 80
160 return
end
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
if (iw .gt. 0) sum=sum+1.0d0/dfloat(i*i-2*i)
iw=-iw
40 f=f+sum*sum
return
end
!
!     The Chebyquad test problem (Fletcher, 1965) for N = 2,4,6 and 8,
!     with NPT = 2N+1.
!
implicit real*8 (a-h,o-z)
dimension x(10),w(10000)
iprint=2
maxfun=5000
rhoend=1.0d-6
do 30 n=2,8,2
npt=2*n+1
do 10 i=1,n
10 x(i)=dfloat(i)/dfloat(n+1)
rhobeg=0.2d0*x(1)
print 20, n,npt
20 format (//4x,'Results with N =',i2,' and NPT =',i3)
call newuoa (n,npt,x,rhobeg,rhoend,iprint,maxfun,w)
30 continue
stop
end
subroutine newuoa (n,npt,x,rhobeg,rhoend,iprint,maxfun,w)
implicit real*8 (a-h,o-z)
dimension x(*),w(*)
!
!     This subroutine seeks the least value of a function of many variables,
!     by a trust region method that forms quadratic models by interpolation.
!     There can be some freedom in the interpolation conditions, which is
!     taken up by minimizing the Frobenius norm of the change to the second
!     derivative of the quadratic model, beginning with a zero matrix. The
!     arguments of the subroutine are as follows.
!
!     N must be set to the number of variables and must be at least two.
!     NPT is the number of interpolation conditions. Its value must be in the
!       interval [N+2,(N+1)(N+2)/2].
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
!     (NPT+13)*(NPT+N)+3*N*(N+3)/2.
!
!     SUBROUTINE CALFUN (N,X,F) must be provided by the user. It must set F to
!     the value of the objective function for the variables X(1),X(2),...,X(N).
!
!     Partition the working space array, so that different parts of it can be
!     treated separately by the subroutine that performs the main calculation.
!
np=n+1
nptm=npt-np
if (npt .lt. n+2 .or. npt .gt. ((n+2)*np)/2) then
    print 10
10     format (/4x,'Return from NEWUOA because NPT is not in', &
      ' the required interval')
    go to 20
end if
ndim=npt+n
ixb=1
ixo=ixb+n
ixn=ixo+n
ixp=ixn+n
ifv=ixp+n*npt
igq=ifv+npt
ihq=igq+n
ipq=ihq+(n*np)/2
ibmat=ipq+npt
izmat=ibmat+ndim*n
id=izmat+npt*nptm
ivl=id+n
iw=ivl+ndim
!
!     The above settings provide a partition of W for subroutine NEWUOB.
!     The partition requires the first NPT*(NPT+N)+5*N*(N+3)/2 elements of
!     W plus the space that is needed by the last array of NEWUOB.
!
call newuob (n,npt,x,rhobeg,rhoend,iprint,maxfun,w(ixb), &
  w(ixo),w(ixn),w(ixp),w(ifv),w(igq),w(ihq),w(ipq),w(ibmat), &
  w(izmat),ndim,w(id),w(ivl),w(iw))
20 return
end
subroutine newuob (n,npt,x,rhobeg,rhoend,iprint,maxfun,xbase, &
  xopt,xnew,xpt,fval,gq,hq,pq,bmat,zmat,ndim,d,vlag,w)
implicit real*8 (a-h,o-z)
dimension x(*),xbase(*),xopt(*),xnew(*),xpt(npt,*),fval(*), &
  gq(*),hq(*),pq(*),bmat(ndim,*),zmat(npt,*),d(*),vlag(*),w(*)
!
!     The arguments N, NPT, X, RHOBEG, RHOEND, IPRINT and MAXFUN are identical
!       to the corresponding arguments in SUBROUTINE NEWUOA.
!     XBASE will hold a shift of origin that should reduce the contributions
!       from rounding errors to values of the model and Lagrange functions.
!     XOPT will be set to the displacement from XBASE of the vector of
!       variables that provides the least calculated F so far.
!     XNEW will be set to the displacement from XBASE of the vector of
!       variables for the current calculation of F.
!     XPT will contain the interpolation point coordinates relative to XBASE.
!     FVAL will hold the values of F at the interpolation points.
!     GQ will hold the gradient of the quadratic model at XBASE.
!     HQ will hold the explicit second derivatives of the quadratic model.
!     PQ will contain the parameters of the implicit second derivatives of
!       the quadratic model.
!     BMAT will hold the last N columns of H.
!     ZMAT will hold the factorization of the leading NPT by NPT submatrix of
!       H, this factorization being ZMAT times Diag(DZ) times ZMAT^T, where
!       the elements of DZ are plus or minus one, as specified by IDZ.
!     NDIM is the first dimension of BMAT and has the value NPT+N.
!     D is reserved for trial steps from XOPT.
!     VLAG will contain the values of the Lagrange functions at a new point X.
!       They are part of a product that requires VLAG to be of length NDIM.
!     The array W will be used for working space. Its length must be at least
!       10*NDIM = 10*(NPT+N).
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
nftest=max0(maxfun,1)
!
!     Set the initial elements of XPT, BMAT, HQ, PQ and ZMAT to zero.
!
do 20 j=1,n
xbase(j)=x(j)
do 10 k=1,npt
10 xpt(k,j)=zero
do 20 i=1,ndim
20 bmat(i,j)=zero
do 30 ih=1,nh
30 hq(ih)=zero
do 40 k=1,npt
pq(k)=zero
do 40 j=1,nptm
40 zmat(k,j)=zero
!
!     Begin the initialization procedure. NF becomes one more than the number
!     of function values so far. The coordinates of the displacement of the
!     next initial interpolation point from XBASE are set in XPT(NF,.).
!
rhosq=rhobeg*rhobeg
recip=one/rhosq
reciq=dsqrt(half)/rhosq
nf=0
50 nfm=nf
nfmm=nf-n
nf=nf+1
if (nfm .le. 2*n) then
    if (nfm .ge. 1 .and. nfm .le. n) then
        xpt(nf,nfm)=rhobeg
    else if (nfm .gt. n) then
        xpt(nf,nfmm)=-rhobeg
    end if
else
    itemp=(nfmm-1)/n
    jpt=nfm-itemp*n-n
    ipt=jpt+itemp
    if (ipt .gt. n) then
        itemp=jpt
        jpt=ipt-n
        ipt=itemp
    end if
    xipt=rhobeg
    if (fval(ipt+np) .lt. fval(ipt+1)) xipt=-xipt
    xjpt=rhobeg
    if (fval(jpt+np) .lt. fval(jpt+1)) xjpt=-xjpt
    xpt(nf,ipt)=xipt
    xpt(nf,jpt)=xjpt
end if
!
!     Calculate the next value of F, label 70 being reached immediately
!     after this calculation. The least function value so far and its index
!     are required.
!
do 60 j=1,n
60 x(j)=xpt(nf,j)+xbase(j)
goto 310
70 fval(nf)=f
if (nf .eq. 1) then
    fbeg=f
    fopt=f
    kopt=1
else if (f .lt. fopt) then
    fopt=f
    kopt=nf
end if
!
!     Set the nonzero initial elements of BMAT and the quadratic model in
!     the cases when NF is at most 2*N+1.
!
if (nfm .le. 2*n) then
    if (nfm .ge. 1 .and. nfm .le. n) then
        gq(nfm)=(f-fbeg)/rhobeg
        if (npt .lt. nf+n) then
            bmat(1,nfm)=-one/rhobeg
            bmat(nf,nfm)=one/rhobeg
            bmat(npt+nfm,nfm)=-half*rhosq
        end if
    else if (nfm .gt. n) then
        bmat(nf-n,nfmm)=half/rhobeg
        bmat(nf,nfmm)=-half/rhobeg
        zmat(1,nfmm)=-reciq-reciq
        zmat(nf-n,nfmm)=reciq
        zmat(nf,nfmm)=reciq
        ih=(nfmm*(nfmm+1))/2
        temp=(fbeg-f)/rhobeg
        hq(ih)=(gq(nfmm)-temp)/rhobeg
        gq(nfmm)=half*(gq(nfmm)+temp)
    end if
!
!     Set the off-diagonal second derivatives of the Lagrange functions and
!     the initial quadratic model.
!
else
    ih=(ipt*(ipt-1))/2+jpt
    if (xipt .lt. zero) ipt=ipt+n
    if (xjpt .lt. zero) jpt=jpt+n
    zmat(1,nfmm)=recip
    zmat(nf,nfmm)=recip
    zmat(ipt+1,nfmm)=-recip
    zmat(jpt+1,nfmm)=-recip
    hq(ih)=(fbeg-fval(ipt+1)-fval(jpt+1)+f)/(xipt*xjpt)
end if
if (nf .lt. npt) goto 50
!
!     Begin the iterative procedure, because the initial model is complete.
!
rho=rhobeg
delta=rho
idz=1
diffa=zero
diffb=zero
itest=0
xoptsq=zero
do 80 i=1,n
xopt(i)=xpt(kopt,i)
80 xoptsq=xoptsq+xopt(i)**2
90 nfsav=nf
!
!     Generate the next trust region step and test its length. Set KNEW
!     to -1 if the purpose of the next F will be to improve the model.
!
100 knew=0
call trsapp (n,npt,xopt,xpt,gq,hq,pq,delta,d,w,w(np), &
  w(np+n),w(np+2*n),crvmin)
dsq=zero
do 110 i=1,n
110 dsq=dsq+d(i)**2
dnorm=dmin1(delta,dsqrt(dsq))
if (dnorm .lt. half*rho) then
    knew=-1
    delta=tenth*delta
    ratio=-1.0d0
    if (delta .le. 1.5d0*rho) delta=rho
    if (nf .le. nfsav+2) goto 460
    temp=0.125d0*crvmin*rho*rho
    if (temp .le. dmax1(diffa,diffb,diffc)) goto 460
    goto 490
end if
!
!     Shift XBASE if XOPT may be too far from XBASE. First make the changes
!     to BMAT that do not depend on ZMAT.
!
120 if (dsq .le. 1.0d-3*xoptsq) then
    tempq=0.25d0*xoptsq
    do 140 k=1,npt
    sum=zero
    do 130 i=1,n
130     sum=sum+xpt(k,i)*xopt(i)
    temp=pq(k)*sum
    sum=sum-half*xoptsq
    w(npt+k)=sum
    do 140 i=1,n
    gq(i)=gq(i)+temp*xpt(k,i)
    xpt(k,i)=xpt(k,i)-half*xopt(i)
    vlag(i)=bmat(k,i)
    w(i)=sum*xpt(k,i)+tempq*xopt(i)
    ip=npt+i
    do 140 j=1,i
140     bmat(ip,j)=bmat(ip,j)+vlag(i)*w(j)+w(i)*vlag(j)
!
!     Then the revisions of BMAT that depend on ZMAT are calculated.
!
    do 180 k=1,nptm
    sumz=zero
    do 150 i=1,npt
    sumz=sumz+zmat(i,k)
150     w(i)=w(npt+i)*zmat(i,k)
    do 170 j=1,n
    sum=tempq*sumz*xopt(j)
    do 160 i=1,npt
160     sum=sum+w(i)*xpt(i,j)
    vlag(j)=sum
    if (k .lt. idz) sum=-sum
    do 170 i=1,npt
170     bmat(i,j)=bmat(i,j)+sum*zmat(i,k)
    do 180 i=1,n
    ip=i+npt
    temp=vlag(i)
    if (k .lt. idz) temp=-temp
    do 180 j=1,i
180     bmat(ip,j)=bmat(ip,j)+temp*vlag(j)
!
!     The following instructions complete the shift of XBASE, including
!     the changes to the parameters of the quadratic model.
!
    ih=0
    do 200 j=1,n
    w(j)=zero
    do 190 k=1,npt
    w(j)=w(j)+pq(k)*xpt(k,j)
190     xpt(k,j)=xpt(k,j)-half*xopt(j)
    do 200 i=1,j
    ih=ih+1
    if (i .lt. j) gq(j)=gq(j)+hq(ih)*xopt(i)
    gq(i)=gq(i)+hq(ih)*xopt(j)
    hq(ih)=hq(ih)+w(i)*xopt(j)+xopt(i)*w(j)
200     bmat(npt+i,j)=bmat(npt+j,i)
    do 210 j=1,n
    xbase(j)=xbase(j)+xopt(j)
210     xopt(j)=zero
    xoptsq=zero
end if
!
!     Pick the model step if KNEW is positive. A different choice of D
!     may be made later, if the choice of D by BIGLAG causes substantial
!     cancellation in DENOM.
!
if (knew .gt. 0) then
    call biglag (n,npt,xopt,xpt,bmat,zmat,idz,ndim,knew,dstep, &
      d,alpha,vlag,vlag(npt+1),w,w(np),w(np+n))
end if
!
!     Calculate VLAG and BETA for the current choice of D. The first NPT
!     components of W_check will be held in W.
!
do 230 k=1,npt
suma=zero
sumb=zero
sum=zero
do 220 j=1,n
suma=suma+xpt(k,j)*d(j)
sumb=sumb+xpt(k,j)*xopt(j)
220 sum=sum+bmat(k,j)*d(j)
w(k)=suma*(half*suma+sumb)
230 vlag(k)=sum
beta=zero
do 250 k=1,nptm
sum=zero
do 240 i=1,npt
240 sum=sum+zmat(i,k)*w(i)
if (k .lt. idz) then
    beta=beta+sum*sum
    sum=-sum
else
    beta=beta-sum*sum
end if
do 250 i=1,npt
250 vlag(i)=vlag(i)+sum*zmat(i,k)
bsum=zero
dx=zero
do 280 j=1,n
sum=zero
do 260 i=1,npt
260 sum=sum+w(i)*bmat(i,j)
bsum=bsum+sum*d(j)
jp=npt+j
do 270 k=1,n
270 sum=sum+bmat(jp,k)*d(k)
vlag(jp)=sum
bsum=bsum+sum*d(j)
280 dx=dx+d(j)*xopt(j)
beta=dx*dx+dsq*(xoptsq+dx+dx+half*dsq)+beta-bsum
vlag(kopt)=vlag(kopt)+one
!
!     If KNEW is positive and if the cancellation in DENOM is unacceptable,
!     then BIGDEN calculates an alternative model step, XNEW being used for
!     working space.
!
if (knew .gt. 0) then
    temp=one+alpha*beta/vlag(knew)**2
    if (dabs(temp) .le. 0.8d0) then
        call bigden (n,npt,xopt,xpt,bmat,zmat,idz,ndim,kopt, &
          knew,d,w,vlag,beta,xnew,w(ndim+1),w(6*ndim+1))
    end if
end if
!
!     Calculate the next value of the objective function.
!
290 do 300 i=1,n
xnew(i)=xopt(i)+d(i)
300 x(i)=xbase(i)+xnew(i)
nf=nf+1
310 if (nf .gt. nftest) then
    nf=nf-1
    if (iprint .gt. 0) print 320
320     format (/4x,'Return from NEWUOA because CALFUN has been', &
      ' called MAXFUN times.')
    goto 530
end if
call calfun (n,x,f)
if (iprint .eq. 3) then
    print 330, nf,f,(x(i),i=1,n)
330      format (/4x,'Function number',i6,'    F =',1pd18.10, &
       '    The corresponding X is:'/(2x,5d15.6))
end if
if (nf .le. npt) goto 70
if (knew .eq. -1) goto 530
!
!     Use the quadratic model to predict the change in F due to the step D,
!     and set DIFF to the error of this prediction.
!
vquad=zero
ih=0
do 340 j=1,n
vquad=vquad+d(j)*gq(j)
do 340 i=1,j
ih=ih+1
temp=d(i)*xnew(j)+d(j)*xopt(i)
if (i .eq. j) temp=half*temp
340 vquad=vquad+temp*hq(ih)
do 350 k=1,npt
350 vquad=vquad+pq(k)*w(k)
diff=f-fopt-vquad
diffc=diffb
diffb=diffa
diffa=dabs(diff)
if (dnorm .gt. rho) nfsav=nf
!
!     Update FOPT and XOPT if the new F is the least value of the objective
!     function so far. The branch when KNEW is positive occurs if D is not
!     a trust region step.
!
fsave=fopt
if (f .lt. fopt) then
    fopt=f
    xoptsq=zero
    do 360 i=1,n
    xopt(i)=xnew(i)
360     xoptsq=xoptsq+xopt(i)**2
end if
ksave=knew
if (knew .gt. 0) goto 410
!
!     Pick the next value of DELTA after a trust region step.
!
if (vquad .ge. zero) then
    if (iprint .gt. 0) print 370
370     format (/4x,'Return from NEWUOA because a trust', &
      ' region step has failed to reduce Q.')
    goto 530
end if
ratio=(f-fsave)/vquad
if (ratio .le. tenth) then
    delta=half*dnorm
else if (ratio .le. 0.7d0) then
    delta=dmax1(half*delta,dnorm)
else
    delta=dmax1(half*delta,dnorm+dnorm)
end if
if (delta .le. 1.5d0*rho) delta=rho
!
!     Set KNEW to the index of the next interpolation point to be deleted.
!
rhosq=dmax1(tenth*delta,rho)**2
ktemp=0
detrat=zero
if (f .ge. fsave) then
    ktemp=kopt
    detrat=one
end if
do 400 k=1,npt
hdiag=zero
do 380 j=1,nptm
temp=one
if (j .lt. idz) temp=-one
380 hdiag=hdiag+temp*zmat(k,j)**2
temp=dabs(beta*hdiag+vlag(k)**2)
distsq=zero
do 390 j=1,n
390 distsq=distsq+(xpt(k,j)-xopt(j))**2
if (distsq .gt. rhosq) temp=temp*(distsq/rhosq)**3
if (temp .gt. detrat .and. k .ne. ktemp) then
    detrat=temp
    knew=k
end if
400 continue
if (knew .eq. 0) goto 460
!
!     Update BMAT, ZMAT and IDZ, so that the KNEW-th interpolation point
!     can be moved. Begin the updating of the quadratic model, starting
!     with the explicit second derivative term.
!
410 call update (n,npt,bmat,zmat,idz,ndim,vlag,beta,knew,w)
fval(knew)=f
ih=0
do 420 i=1,n
temp=pq(knew)*xpt(knew,i)
do 420 j=1,i
ih=ih+1
420 hq(ih)=hq(ih)+temp*xpt(knew,j)
pq(knew)=zero
!
!     Update the other second derivative parameters, and then the gradient
!     vector of the model. Also include the new interpolation point.
!
do 440 j=1,nptm
temp=diff*zmat(knew,j)
if (j .lt. idz) temp=-temp
do 440 k=1,npt
440 pq(k)=pq(k)+temp*zmat(k,j)
gqsq=zero
do 450 i=1,n
gq(i)=gq(i)+diff*bmat(knew,i)
gqsq=gqsq+gq(i)**2
450 xpt(knew,i)=xnew(i)
!
!     If a trust region step makes a small change to the objective function,
!     then calculate the gradient of the least Frobenius norm interpolant at
!     XBASE, and store it in W, using VLAG for a vector of right hand sides.
!
if (ksave .eq. 0 .and. delta .eq. rho) then
    if (dabs(ratio) .gt. 1.0d-2) then
        itest=0
    else
        do 700 k=1,npt
700         vlag(k)=fval(k)-fval(kopt)
        gisq=zero
        do 720 i=1,n
        sum=zero
        do 710 k=1,npt
710         sum=sum+bmat(k,i)*vlag(k)
        gisq=gisq+sum*sum
720         w(i)=sum
!
!     Test whether to replace the new quadratic model by the least Frobenius
!     norm interpolant, making the replacement if the test is satisfied.
!
        itest=itest+1
        if (gqsq .lt. 1.0d2*gisq) itest=0
        if (itest .ge. 3) then
            do 730 i=1,n
730             gq(i)=w(i)
            do 740 ih=1,nh
740             hq(ih)=zero
            do 760 j=1,nptm
            w(j)=zero
            do 750 k=1,npt
750             w(j)=w(j)+vlag(k)*zmat(k,j)
760             if (j .lt. idz) w(j)=-w(j)
            do 770 k=1,npt
            pq(k)=zero
            do 770 j=1,nptm
770             pq(k)=pq(k)+zmat(k,j)*w(j)
            itest=0
        end if
    end if
end if
if (f .lt. fsave) kopt=knew
!
!     If a trust region step has provided a sufficient decrease in F, then
!     branch for another trust region calculation. The case KSAVE>0 occurs
!     when the new function value was calculated by a model step.
!
if (f .le. fsave+tenth*vquad) goto 100
if (ksave .gt. 0) goto 100
!
!     Alternatively, find out if the interpolation points are close enough
!     to the best point so far.
!
knew=0
460 distsq=4.0d0*delta*delta
do 480 k=1,npt
sum=zero
do 470 j=1,n
470 sum=sum+(xpt(k,j)-xopt(j))**2
if (sum .gt. distsq) then
    knew=k
    distsq=sum
end if
480 continue
!
!     If KNEW is positive, then set DSTEP, and branch back for the next
!     iteration, which will generate a "model step".
!
if (knew .gt. 0) then
    dstep=dmax1(dmin1(tenth*dsqrt(distsq),half*delta),rho)
    dsq=dstep*dstep
    goto 120
end if
if (ratio .gt. zero) goto 100
if (dmax1(delta,dnorm) .gt. rho) goto 100
!
!     The calculations with the current value of RHO are complete. Pick the
!     next values of RHO and DELTA.
!
490 if (rho .gt. rhoend) then
    delta=half*rho
    ratio=rho/rhoend
    if (ratio .le. 16.0d0) then
        rho=rhoend
    else if (ratio .le. 250.0d0) then
        rho=dsqrt(ratio)*rhoend
    else
        rho=tenth*rho
    end if
    delta=dmax1(delta,rho)
    if (iprint .ge. 2) then
        if (iprint .ge. 3) print 500
500         format (5x)
        print 510, rho,nf
510         format (/4x,'New RHO =',1pd11.4,5x,'Number of', &
          ' function values =',i6)
        print 520, fopt,(xbase(i)+xopt(i),i=1,n)
520         format (4x,'Least value of F =',1pd23.15,9x, &
          'The corresponding X is:'/(2x,5d15.6))
    end if
    goto 90
end if
!
!     Return from the calculation, after another Newton-Raphson step, if
!     it is too short to have been tried before.
!
if (knew .eq. -1) goto 290
530 if (fopt .le. f) then
    do 540 i=1,n
540     x(i)=xbase(i)+xopt(i)
    f=fopt
end if
if (iprint .ge. 1) then
    print 550, nf
550     format (/4x,'At the return from NEWUOA',5x, &
      'Number of function values =',i6)
    print 520, f,(x(i),i=1,n)
end if
return
end
subroutine trsapp (n,npt,xopt,xpt,gq,hq,pq,delta,step, &
  d,g,hd,hs,crvmin)
implicit real*8 (a-h,o-z)
dimension xopt(*),xpt(npt,*),gq(*),hq(*),pq(*),step(*), &
  d(*),g(*),hd(*),hs(*)
!
!     N is the number of variables of a quadratic objective function, Q say.
!     The arguments NPT, XOPT, XPT, GQ, HQ and PQ have their usual meanings,
!       in order to define the current quadratic model Q.
!     DELTA is the trust region radius, and has to be positive.
!     STEP will be set to the calculated trial step.
!     The arrays D, G, HD and HS will be used for working space.
!     CRVMIN will be set to the least curvature of H along the conjugate
!       directions that occur, except that it is set to zero if STEP goes
!       all the way to the trust region boundary.
!
!     The calculation of STEP begins with the truncated conjugate gradient
!     method. If the boundary of the trust region is reached, then further
!     changes to STEP may be made, each one being in the 2D space spanned
!     by the current STEP and the corresponding gradient of Q. Thus STEP
!     should provide a substantial reduction to Q within the trust region.
!
!     Initialization, which includes setting HD to H times XOPT.
!
half=0.5d0
zero=0.0d0
twopi=8.0d0*datan(1.0d0)
delsq=delta*delta
iterc=0
itermax=n
itersw=itermax
do 10 i=1,n
10 d(i)=xopt(i)
goto 170
!
!     Prepare for the first line search.
!
20 qred=zero
dd=zero
do 30 i=1,n
step(i)=zero
hs(i)=zero
g(i)=gq(i)+hd(i)
d(i)=-g(i)
30 dd=dd+d(i)**2
crvmin=zero
if (dd .eq. zero) goto 160
ds=zero
ss=zero
gg=dd
ggbeg=gg
!
!     Calculate the step to the trust region boundary and the product HD.
!
40 iterc=iterc+1
temp=delsq-ss
bstep=temp/(ds+dsqrt(ds*ds+dd*temp))
goto 170
50 dhd=zero
do 60 j=1,n
60 dhd=dhd+d(j)*hd(j)
!
!     Update CRVMIN and set the step-length ALPHA.
!
alpha=bstep
if (dhd .gt. zero) then
    temp=dhd/dd
    if (iterc .eq. 1) crvmin=temp
    crvmin=dmin1(crvmin,temp)
    alpha=dmin1(alpha,gg/dhd)
end if
qadd=alpha*(gg-half*alpha*dhd)
qred=qred+qadd
!
!     Update STEP and HS.
!
ggsav=gg
gg=zero
do 70 i=1,n
step(i)=step(i)+alpha*d(i)
hs(i)=hs(i)+alpha*hd(i)
70 gg=gg+(g(i)+hs(i))**2
!
!     Begin another conjugate direction iteration if required.
!
if (alpha .lt. bstep) then
    if (qadd .le. 0.01d0*qred) goto 160
    if (gg .le. 1.0d-4*ggbeg) goto 160
    if (iterc .eq. itermax) goto 160
    temp=gg/ggsav
    dd=zero
    ds=zero
    ss=zero
    do 80 i=1,n
    d(i)=temp*d(i)-g(i)-hs(i)
    dd=dd+d(i)**2
    ds=ds+d(i)*step(i)
80     ss=ss+step(i)**2
    if (ds .le. zero) goto 160
    if (ss .lt. delsq) goto 40
end if
crvmin=zero
itersw=iterc
!
!     Test whether an alternative iteration is required.
!
90 if (gg .le. 1.0d-4*ggbeg) goto 160
sg=zero
shs=zero
do 100 i=1,n
sg=sg+step(i)*g(i)
100 shs=shs+step(i)*hs(i)
sgk=sg+shs
angtest=sgk/dsqrt(gg*delsq)
if (angtest .le. -0.99d0) goto 160
!
!     Begin the alternative iteration by calculating D and HD and some
!     scalar products.
!
iterc=iterc+1
temp=dsqrt(delsq*gg-sgk*sgk)
tempa=delsq/temp
tempb=sgk/temp
do 110 i=1,n
110 d(i)=tempa*(g(i)+hs(i))-tempb*step(i)
goto 170
120 dg=zero
dhd=zero
dhs=zero
do 130 i=1,n
dg=dg+d(i)*g(i)
dhd=dhd+hd(i)*d(i)
130 dhs=dhs+hd(i)*step(i)
!
!     Seek the value of the angle that minimizes Q.
!
cf=half*(shs-dhd)
qbeg=sg+cf
qsav=qbeg
qmin=qbeg
isave=0
iu=49
temp=twopi/dfloat(iu+1)
do 140 i=1,iu
angle=dfloat(i)*temp
cth=dcos(angle)
sth=dsin(angle)
qnew=(sg+cf*cth)*cth+(dg+dhs*cth)*sth
if (qnew .lt. qmin) then
    qmin=qnew
    isave=i
    tempa=qsav
else if (i .eq. isave+1) then
    tempb=qnew
end if
140 qsav=qnew
if (isave .eq. zero) tempa=qnew
if (isave .eq. iu) tempb=qbeg
angle=zero
if (tempa .ne. tempb) then
    tempa=tempa-qmin
    tempb=tempb-qmin
    angle=half*(tempa-tempb)/(tempa+tempb)
end if
angle=temp*(dfloat(isave)+angle)
!
!     Calculate the new STEP and HS. Then test for convergence.
!
cth=dcos(angle)
sth=dsin(angle)
reduc=qbeg-(sg+cf*cth)*cth-(dg+dhs*cth)*sth
gg=zero
do 150 i=1,n
step(i)=cth*step(i)+sth*d(i)
hs(i)=cth*hs(i)+sth*hd(i)
150 gg=gg+(g(i)+hs(i))**2
qred=qred+reduc
ratio=reduc/qred
if (iterc .lt. itermax .and. ratio .gt. 0.01d0) goto 90
160 return
!
!     The following instructions act as a subroutine for setting the vector
!     HD to the vector D multiplied by the second derivative matrix of Q.
!     They are called from three different places, which are distinguished
!     by the value of ITERC.
!
170 do 180 i=1,n
180 hd(i)=zero
do 200 k=1,npt
temp=zero
do 190 j=1,n
190 temp=temp+xpt(k,j)*d(j)
temp=temp*pq(k)
do 200 i=1,n
200 hd(i)=hd(i)+temp*xpt(k,i)
ih=0
do 210 j=1,n
do 210 i=1,j
ih=ih+1
if (i .lt. j) hd(j)=hd(j)+hq(ih)*d(i)
210 hd(i)=hd(i)+hq(ih)*d(j)
if (iterc .eq. 0) goto 20
if (iterc .le. itersw) goto 50
goto 120
end
subroutine update (n,npt,bmat,zmat,idz,ndim,vlag,beta,knew,w)
implicit real*8 (a-h,o-z)
dimension bmat(ndim,*),zmat(npt,*),vlag(*),w(*)
!
!     The arrays BMAT and ZMAT with IDZ are updated, in order to shift the
!     interpolation point that has index KNEW. On entry, VLAG contains the
!     components of the vector Theta*Wcheck+e_b of the updating formula
!     (6.11), and BETA holds the value of the parameter that has this name.
!     The vector W is used for working space.
!
!     Set some constants.
!
one=1.0d0
zero=0.0d0
nptm=npt-n-1
!
!     Apply the rotations that put zeros in the KNEW-th row of ZMAT.
!
jl=1
do 20 j=2,nptm
if (j .eq. idz) then
    jl=idz
else if (zmat(knew,j) .ne. zero) then
    temp=dsqrt(zmat(knew,jl)**2+zmat(knew,j)**2)
    tempa=zmat(knew,jl)/temp
    tempb=zmat(knew,j)/temp
    do 10 i=1,npt
    temp=tempa*zmat(i,jl)+tempb*zmat(i,j)
    zmat(i,j)=tempa*zmat(i,j)-tempb*zmat(i,jl)
10     zmat(i,jl)=temp
    zmat(knew,j)=zero
end if
20 continue
!
!     Put the first NPT components of the KNEW-th column of HLAG into W,
!     and calculate the parameters of the updating formula.
!
tempa=zmat(knew,1)
if (idz .ge. 2) tempa=-tempa
if (jl .gt. 1) tempb=zmat(knew,jl)
do 30 i=1,npt
w(i)=tempa*zmat(i,1)
if (jl .gt. 1) w(i)=w(i)+tempb*zmat(i,jl)
30 continue
alpha=w(knew)
tau=vlag(knew)
tausq=tau*tau
denom=alpha*beta+tausq
vlag(knew)=vlag(knew)-one
!
!     Complete the updating of ZMAT when there is only one nonzero element
!     in the KNEW-th row of the new matrix ZMAT, but, if IFLAG is set to one,
!     then the first column of ZMAT will be exchanged with another one later.
!
iflag=0
if (jl .eq. 1) then
    temp=dsqrt(dabs(denom))
    tempb=tempa/temp
    tempa=tau/temp
    do 40 i=1,npt
40     zmat(i,1)=tempa*zmat(i,1)-tempb*vlag(i)
    if (idz .eq. 1 .and. temp .lt. zero) idz=2
    if (idz .ge. 2 .and. temp .ge. zero) iflag=1
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
    scalb=scala*dsqrt(dabs(denom))
    do 50 i=1,npt
    zmat(i,ja)=scala*(tau*zmat(i,ja)-temp*vlag(i))
50     zmat(i,jb)=scalb*(zmat(i,jb)-tempa*w(i)-tempb*vlag(i))
    if (denom .le. zero) then
        if (beta .lt. zero) idz=idz+1
        if (beta .ge. zero) iflag=1
    end if
end if
!
!     IDZ is reduced in the following case, and usually the first column
!     of ZMAT is exchanged with a later one.
!
if (iflag .eq. 1) then
    idz=idz-1
    do 60 i=1,npt
    temp=zmat(i,1)
    zmat(i,1)=zmat(i,idz)
60     zmat(i,idz)=temp
end if
!
!     Finally, update the matrix BMAT.
!
do 70 j=1,n
jp=npt+j
w(jp)=bmat(knew,j)
tempa=(alpha*vlag(jp)-tau*w(jp))/denom
tempb=(-beta*w(jp)-tau*vlag(jp))/denom
do 70 i=1,jp
bmat(i,j)=bmat(i,j)+tempa*vlag(i)+tempb*w(i)
if (i .gt. npt) bmat(jp,i-npt)=bmat(i,j)
70 continue
return
end
