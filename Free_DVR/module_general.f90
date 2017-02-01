module general
contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!% a better version of random_seed %!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

subroutine init_random_seed()
  implicit none

  integer ii,nn,value(1:8)
  integer,allocatable :: seed(:)
  double precision flagd

  call random_seed(size=nn)
  allocate(seed(nn))
  call date_and_time(values=value)
  seed = value(8)+37*(/(ii-1,ii=1,nn)/)
  call random_seed(put=seed)
  deallocate(seed)

  do ii=1,value(6)*3600+value(7)*60+value(8)
    call random_number(flagd)
  enddo
end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!% make random number more random %!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

subroutine more_random()
  implicit none

  integer ii,value(1:8)
  double precision flagd

  call date_and_time(values=value)
  do ii=1,value(8)/100
    call random_number(flagd)
  enddo
end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!% gaussian random number generator using box-muller method   %!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!% ref: http://en.wikipedia.org/wiki/gaussian_random_variable %!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

function gaussian_random_number(mean,sigma)
  implicit none

  double precision gaussian_random_number,mean,sigma,pi,r1,r2

  pi=4.0d0*datan(1.0d0)
  call more_random()
  call random_number(r1)
  call more_random()
  call random_number(r2)
  gaussian_random_number=mean+sigma*dsqrt(-2.0d0*dlog(r1))*dcos(2.0d0*pi*r2)
end function

function gaussian_random_number_fast(mean,sigma)
  implicit none

  double precision gaussian_random_number_fast,mean,sigma,pi,r1,r2

  pi=4.0d0*datan(1.0d0)
  call random_number(r1)
  call random_number(r2)
  gaussian_random_number_fast=mean+sigma*dsqrt(-2.0d0*dlog(r1))*dcos(2.0d0*pi*r2)
end function

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!% diagonalize a symmetric matrix using rs   %!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!% this is a rearranged version of dsub1.f90 %!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

subroutine rs(nm,n,a,w,matz,z,fv1,fv2,ierr)
  implicit real*8 (a-h,o-z)
  dimension a(nm,n),w(n),z(nm,n),fv1(n),fv2(n)

  if(n.le.nm) goto 10
  ierr=10*n
  goto 50

  10 if(matz.ne.0) goto 20
  call tred1(nm,n,a,w,fv1,fv2)
  call tqlrat(n,w,fv2,ierr)
  goto 50

  20 call tred2(nm,n,a,w,fv1,z)
  call tql2(nm,n,w,fv1,z,ierr)
  50 return
end subroutine

subroutine tred2(nm,n,a,d,e,z)
  implicit real*8 (a-h,o-z) 
  dimension a(nm,n),d(n),e(n),z(nm,n)

  do 100 i=1,n
  do 100 j=1,i
  z(i,j)=a(i,j)
  100 continue
  if(n.eq.1) goto 320
  do 300 ii=2,n
  i=n+2-ii
  l=i-1
  h=0.0d0
  scale=0.0d0
  if(l.lt.2) goto 130
  do 120 k=1,l
  120 scale=scale+dabs(z(i,k))
  if(scale.ne.0.0d0) goto 1120
  130 e(i)=z(i,l)
  goto 290
  1120 do 150 k=1,l
  z(i,k)=z(i,k)/scale
  h=h+z(i,k)*z(i,k)
  150 continue
  f=z(i,l)
  g=-dsign(dsqrt(h),f)
  e(i)=scale*g
  h=h-f*g
  z(i,l)=f-g
  f=0.0d0
  do 240 j=1,l
  z(j,i)=z(i,j)/(scale*h)
  g=0.0d0
  do 180 k=1,j
  180 g=g+z(j,k)*z(i,k)
  jp1=j+1
  if(l.lt.jp1) goto 220
  do 200 k=jp1,l
  200 g=g+z(k,j)*z(i,k)
  220 e(j)=g/h
  f=f+e(j)*z(i,j)
  240 continue
  hh=f/(h+h)
  do 260 j=1,l
  f=z(i,j)
  g=e(j)-hh*f
  e(j)=g
  do 260 k=1,j
  z(j,k)=z(j,k)-f*e(k)-g*z(i,k)
  260 continue
  do 280 k=1,l
  280 z(i,k)=scale*z(i,k)
  290 d(i)=h
  300 continue
  320 d(1)=0.0d0
  e(1)=0.0d0
  do 500 i=1,n
  l=i-1
  if(d(i).eq.0.0d0) goto 380
  do 360 j=1,l
  g=0.0d0
  do 340 k=1,l
  340 g=g+z(i,k)*z(k,j)
  do 360 k=1,l
  z(k,j)=z(k,j)-g*z(k,i)
  360 continue
  380 d(i)=z(i,i)
  z(i,i)=1.0d0
  if(l.lt.1) goto 500
  do 400 j=1,l
  z(i,j)=0.0d0
  z(j,i)=0.0d0
  400 continue
  500 continue
  return
end subroutine

subroutine tql2(nm,n,d,e,z,ierr)
  implicit real*8 (a-h,o-z)
  real*8 machep
  dimension d(n),e(n),z(nm,n)

  machep=1.0d-30
  ierr=0

  one=1.0d0
  two=2.0d0
  zero=0.0d0

  if(n.eq.1) goto 1001
  do 100 i=2,n
  100 e(i-1)=e(i)
  f=zero
  b=zero
  e(n)=zero
  do 240 l=1,n
  j=0
  h=machep*(dabs(d(l))+dabs(e(l)))
  if(b.lt.h)b=h
  do 110 m=l,n
  if(dabs(e(m)).le.b) goto 120
  110 continue
  120 if(m.eq.l) goto 220
  130 if(j.eq.30) goto 1000
  j=j+1
  l1=l+1
  g=d(l)
  p=(d(l1)-g)/(two*e(l))
  r=dsqrt(p*p+one)
  d(l)=e(l)/(p+dsign(r,p))
  h=g-d(l)
  do 1120 i=l1,n
  1120 d(i)=d(i)-h
  f=f+h
  p=d(m)
  c=one
  s=zero
  mml=m-l
  do 200 ii=1,mml
  i=m-ii
  g=c*e(i)
  h=c*p
  if(dabs(p).lt.dabs(e(i))) goto 150
  c=e(i)/p
  r=sqrt(c*c+one)
  e(i+1)=s*p*r
  s=c/r
  c=one/r
  go to 160
  150 c=p/e(i)
  r=dsqrt(c*c+one)
  e(i+1)=s*e(i)*r
  s=one/r
  c=c*s
  160 p=c*d(i)-s*g
  d(i+1)=h+s*(c*g+s*d(i))
  do 180 k=1,n
  h=z(k,i+1)
  z(k,i+1)=s*z(k,i)+c*h
  z(k,i)=c*z(k,i)-s*h
  180 continue
  200 continue
  e(l)=s*p
  d(l)=c*p
  if(dabs(e(l)).gt.b) goto 130
  220 d(l)=d(l)+f
  240 continue
  do 300 ii=2,n
  i=ii-1
  k=i
  p=d(i)
  do 260 j=ii,n
  if(d(j).ge.p) goto 260
  k=j
  p=d(j)
  260 continue
  if(k.eq.i) goto 300
  d(k)=d(i)
  d(i)=p
  do 280 j=1,n
  p=z(j,i)
  z(j,i)=z(j,k)
  z(j,k)=p
  280 continue
  300 continue
  goto 1001
  1000 ierr=l
  1001 return
end subroutine

real*8 function ssum(n,a,nstp)
  real*8 a
  dimension a(n)

  ssum=a(1)
  if(n.le.1) goto 100
  do i=2,n
  ssum=ssum+a(i)
  enddo
  100 continue
end function

subroutine tred1(nm,n,a,d,e,e2)
  implicit real*8 (a-h,o-z) 
  dimension a(nm,n),d(n),e(n),e2(n)

  zero=0.0d0
  one=1.0d0
  two=2.0d0

  do 100 i=1,n
  100 d(i)=a(i,i)

  do 300 ii=1,n
  i=n+1-ii
  l=i-1
  h=zero
  scale=zero
  if(l.lt.1) goto 130
  do 120 k=1,l
  120 scale=scale+dabs(a(i,k))
  if(scale.ne.zero) goto 140
  130 e(i)=zero
  e2(i)=zero
  goto 290
  140 do 150 k=1,l
  a(i,k)=a(i,k)/scale
  150 h=h+a(i,k)*a(i,k)
  e2(i)=scale*scale*h
  f=a(i,l)
  g=-dsign(dsqrt(h),f)
  e(i)=scale*g
  h=h-f*g
  a(i,l)=f-g
  if(l.eq.1) goto 270
  f=zero
  do 240 j=1,l
  g=zero
  do 180 k=1,j
  180 g=g+a(j,k)*a(i,k)
  jp1=j+1
  if(l.lt.jp1) goto 220
  do 200 k=jp1,l
  200 g=g+a(k,j)*a(i,k)
  220 e(j)=g/h
  f=f+e(j)*a(i,j)
  240 continue

  h=f/(h+h)
  do 260 j=1,l
  f=a(i,j)
  g=e(j)-h*f
  e(j)=g
  do 260 k=1,j
  a(j,k)=a(j,k)-f*e(k)-g*a(i,k)
  260 continue
  270 do 280 k=1,l
  280 a(i,k)=scale*a(i,k)
  290 h=d(i)
  d(i)=a(i,i)
  a(i,i)=h
  300 continue
  return
end subroutine

subroutine tqlrat(n,d,e2,ierr)
  implicit real*8 (a-h,o-z)
  dimension d(n),e2(n)

  one=1.0d0
  two=2.0d0
  zero=0.0d0

  machep=1.0d-30
  ierr=0
  if(n.eq.1) goto 1001
  do 100 i=2,n
  100 e2(i-1)=e2(i)
  f=zero
  b=zero
  e2(n)=zero
  do 290 l=1,n
  j=0
  h=machep*(dabs(d(l))+dsqrt(e2(l)))
  if(b.gt.h) goto 105
  b=h
  c=b*b
  105 do 110 m=l,n
  if(e2(m).le.c) goto 120
  110 continue
  120 if(m.eq.l) goto 210
  130 if(j.eq.30) goto 1000
  j=j+1
  l1=l+1
  s=dsqrt(e2(l))
  g=d(l)
  p=(d(l1)-g)/(two*s)
  r=sqrt(p*p+one)
  d(l)=s/(p+dsign(r,p))
  h=g-d(l)
  do 140 i=l1,n
  140 d(i)=d(i)-h
  f=f+h
  g=d(m)
  if(g.eq.zero) g=b
  h=g
  s=zero
  mml=m-l
  do 200 ii=1,mml
  i=m-ii
  p=g*h
  r=p+e2(i)
  e2(i+1)=s*r
  s=e2(i)/r
  d(i+1)=h+s*(h+d(i))
  g=d(i)-e2(i)/g
  if(g.eq.zero)g=b
  h=g*p/r
  200 continue
  e2(l)=s*g
  d(l)=h
  if(h.eq.zero) goto 210
  if(dabs(e2(l)).le.dabs(c/h)) goto 210
  e2(l)=h*e2(l)
  if(e2(l).ne.zero) goto 130
  210 p=d(l)+f

  if(l.eq.1) goto 250
  do 230 ii=2,l
  i=l+2-ii
  if(p.ge.d(i-1)) goto 270
  d(i)=d(i-1)
  230 continue
  250 i=1
  270 d(i)=p
  290 continue
  goto 1001
  1000 ierr=l
  1001 return
end subroutine

subroutine rsp(nm,n,a,w,fv1,z,ierr)
  implicit real*8 (a-h,o-z)
  dimension a(nm,2),w(nm),fv1(nm),z(nm,n)

  zero=0.0d0
  one=1.0d0

  do 100 i=1,n
  do 50 j=1,n
  z(i,j)=zero
  50 continue
  z(i,i)=one
  w(i)=a(i,2)
  fv1(i)=a(i,1)
  100 continue
  call imtql2(nm,n,w,fv1,z,ierr)
end subroutine

subroutine imtql2(nm,n,d,e,z,ierr)
  implicit real*8 (a-h,o-z)
  real*8 machep
  dimension d(n),e(n),z(nm,n)

  zero=0.0d0
  one=1.0d0
  two=2.0d0

  machep=dble(2.**(-37))
  ierr=0
  if(n.eq.1) goto 1001
  do 100 i=2,n
  100 e(i-1)=e(i)
  e(n)=zero
  do 240 l=1,n
  j=0
  105 do 110 m=l,n
  if(m.eq.n) goto 120
  if(dabs(e(m)).le.machep*(dabs(d(m))+dabs(d(m+1)))) goto 120
  110 continue
  120 p=d(l)
  if(m.eq.l) goto 240
  if(j.eq.30) goto 1000
  j=j+1
  g=(d(l+1)-p)/(two*e(l))
  r=sqrt(g*g+one)
  g=d(m)-p+e(l)/(g+dsign(r,g))
  s=one
  c=one
  p=zero
  mml=m-l
  do 200 ii=1,mml
  i=m-ii
  f=s*e(i)
  b=c*e(i)
  if(dabs(f).lt.dabs(g)) goto 150
  c=g/f
  r=dsqrt(c*c+one)
  e(i+1)=f*r
  s=one/r
  c=c*s
  goto 160
  150 s=f/g
  r=dsqrt(s*s+one)
  e(i+1)=g*r
  c=one/r
  s=s*c
  160 g=d(i+1)-p
  r=(d(i)-g)*s+2.0*c*b
  p=s*r
  d(i+1)=g+p
  g=c*r-b
  do 180 k=1,n
  f=z(k,i+1)
  z(k,i+1)=s*z(k,i)+c*f 
  z(k,i)=c*z(k,i)-s*f
  180 continue
  200 continue
  d(l)=d(l)-p
  e(l)=g
  e(m)=zero
  goto 105
  240 continue
  do 300 ii=2,n
  i=ii-1
  k=i
  p=d(i)
  do 260 j=ii,n
  if(d(j).ge.p) goto 260
  k=j
  p=d(j)
  260 continue
  if(k.eq.i) goto 300
  d(k)=d(i)
  d(i)=p
  do 280 j=1,n
  p=z(j,i)
  z(j,i)=z(j,k)
  z(j,k)=p
  280 continue
  300 continue             
  goto 1001
  1000 ierr=l
  1001 return
end subroutine

subroutine rst(nm,n,w,e,z,ierr)
  implicit real*8 (a-h,o-z)
  dimension w(n),e(n),z(nm,n)

  zero=0.0d0
  one=1.0d0
  two=2.0d0

  if(n.le.nm) goto 10
  ierr=10*n
  goto 50
  10 do 40 i=1,n
  do 30 j=1,n
  z(j,i)=zero
  30 continue
  z(i,i)=one
  40 continue
  call imtql2(nm,n,w,e,z,ierr)
  50 return
end subroutine

end module
