module ncp_dvr
use general
implicit none

!! System
real*8 x_left !left boundary of the Hamiltonian  
real*8 x_right !right boundary of the Hamiltonian
real*8 l0
integer npoint !grid points in reactive region
real*8 xs
real*8 dx
real*8 xbar0 !ceter of initial wavepacket
real*8 k0 !momentum of initial wavepacket
real*8 mass !mass of the particle in atomic unit
real*8 dt !timestep in atomic unit
real*8 pop_min !minimum population, criteria of the simulation
real*8 barrier ! height of gaussian like energy barrier
real*8 population
integer,parameter::nlevel=1

!! Loop control, global
integer nsnap
integer nstep
integer isnap
integer istep

!! Loop control, shared but should not interfere each other
integer ipoint
integer jpoint
integer ilevel
integer jlevel

!! Parameters to be seted
integer ncenter

!! Absorbing region
real*8 gama_re !real part of the complex absorbing potential
real*8 gama_im !imaginary part of the complex absorbing potential
complex*16 gama
real*8 nabsorb_left !absorbing region at the left side
real*8 nabsorb_right ! absorbing region at the right side

!! Hamiltonian
complex*16,allocatable :: H_dvr(:,:)
complex*16,allocatable :: H_absorb(:,:)

!! Wavefunction
complex*16,allocatable :: c_dvr(:)
complex*16,allocatable :: c_absorb(:)

!! Transmission and reflection statistics
real*8,allocatable :: reflect(:,:)
real*8,allocatable :: trans(:,:)
real*8,allocatable :: r_snap(:,:)
real*8,allocatable :: t_snap(:,:)

!! Constants
real*8 pi
real*8 pi2

contains

subroutine simulate
  implicit none

  call setup
  call set_wavefunction
  call set_hamiltonian
  call files(0)
 
  isnap=0
  call write_wavefunction

  do isnap=1,nsnap
  do istep=1,nstep
    call rk4_wave
  end do
  call write_wavefunction
  if (population.lt.pop_min) then
    write(6,*) 'isnap=',isnap,'Too few population left, exit'
    exit
  else if (population.gt.1.1d0) then
    write(6,*) 'fuck'
    exit
  end if
  end do

  call files(1) 
end subroutine

subroutine setup
  implicit none

  open(99,file='DVR_NCP.inp')
  read(99,*) x_left
  read(99,*) x_right
  read(99,*) npoint
  read(99,*) nsnap
  read(99,*) nstep
  read(99,*) barrier
  read(99,*) xbar0
  read(99,*) k0
  read(99,*) mass
  read(99,*) dt
  read(99,*) pop_min
  read(99,*) gama_re
  read(99,*) gama_im
  read(99,*) nabsorb_left
  read(99,*) nabsorb_right
  close(99)

  gama=cmplx(gama_re,gama_im)
  l0=x_right-x_left
  xs=10.0d0/k0
  dx=(x_right-x_left)/npoint
  ncenter=int(-x_left/dx)
  allocate(H_dvr(npoint*nlevel,npoint*nlevel))
  allocate(H_absorb(npoint*nlevel,npoint*nlevel))
  allocate(c_dvr(npoint*nlevel))
  allocate(c_absorb(npoint*nlevel))
  allocate(reflect(nlevel,nstep))
  allocate(trans(nlevel,nstep))
  allocate(t_snap(nsnap,nlevel))
  allocate(r_snap(nsnap,nlevel))
 
  H_dvr=0.0d0
  H_absorb=0.0d0
  c_dvr=0.0d0
  c_absorb=0.0d0
  reflect=0.0d0
  trans=0.0d0
  t_snap=0.0d0
  r_snap=0.0d0

  pi=atan(1.0d0)*4.0
  pi2=pi**2
end subroutine

subroutine files(flag)
  implicit none

  integer flag

  if (flag.eq.0) then
    open(unit=11,file='csave.out')
    open(unit=44,file='trsave.out')
  else
    close(11)
    close(22)
    close(33)
  end if

end subroutine

subroutine write_wavefunction
  implicit none
  complex*16 ccc(nlevel)
  complex*16 www(nlevel)
  real*8 r(nlevel)
  real*8 t(nlevel)
  real*8 pop(npoint,nlevel)
  
  population=0.0d0
  r=0.0d0
  t=0.0d0

  do ipoint=1,npoint
    do ilevel=1,nlevel
      ccc(ilevel)=c_absorb(ipoint+(nlevel-1)*npoint)
    end do
    call convert_diabatic_adiabatic
    do ilevel=1,nlevel
      pop(ipoint,ilevel)=abs(www(ilevel))**2
      population=population+pop(ipoint,ilevel)
      if (ipoint.le.ncenter) then
        r(ilevel)=r(ilevel)+pop(ipoint,ilevel)
      else
        t(ilevel)=t(ilevel)+pop(ipoint,ilevel)
      end if
    end do
  end do
  
  do ilevel=1,nlevel
    write(11*ilevel,'(9999f10.5)') (pop(ipoint,ilevel),ipoint=1,npoint)
  end do
  write(6,'(a,i4.4,3(a,f10.6))') 'isnap=',isnap,' pop=',population,' re=',r(1),' tr=',t(1)

!private subroutine of write_wavefunction,used to do representation trabsfirnatuib
contains

subroutine convert_diabatic_adiabatic
  implicit none
  
  real*8 ee(1:nlevel)
  real*8 pp(1:nlevel,1:nlevel)
  real*8 hh(1:nlevel,1:nlevel)
  real*8 fv1(1:nlevel)
  real*8 fv2(1:nlevel)
  integer ierr

  call Potential(dx*(ipoint-1)+x_left,hh)
  call rs(nlevel,nlevel,hh,ee,1,pp,fv1,fv1,ierr)

  www=0.0d0
  www=matmul(pp,ccc)
end subroutine
end subroutine

subroutine Potential(xx,hh)
  implicit none

  real*8 hh(nlevel,nlevel)
  real*8 xx

  hh=barrier*exp(-1.6d0*xx**2)
end subroutine

subroutine set_wavefunction
  implicit none
  real*8 flagd
  real*8 norm
  real*8 xx
  real*8 cc(npoint,nlevel)
 
  do ipoint=1,npoint
    xx=dx*(ipoint-1)+x_left
    flagd=exp(-(xx-xbar0)**2/(4.0d0*xs**2))
    cc(ipoint,1)=cmplx(flagd*cos(k0*(xx-xbar0)),flagd*sin(k0*(xx-xbar0)))
    norm=norm+abs(cc(ipoint,1))**2 
  end do

! There should be a representation adjustment here, but since this is only a 1D
! model, this part will be improved when I have time  
  cc=cc/sqrt(norm)

  do ipoint=1,npoint
    do ilevel=1,nlevel
      c_dvr(ipoint+(ilevel-1)*npoint)=cc(ipoint,ilevel)
    end do
  end do
  c_absorb=c_dvr
  population=1.0d0
end subroutine

subroutine set_Hamiltonian
  implicit none
  complex*16 V_dvr(npoint,nlevel,nlevel)
  complex*16 V_absorb(npoint,nlevel,nlevel)
  real*8 hh(nlevel,nlevel)
  real*8 xx
  integer ilevel

  do ipoint=1,npoint
    xx=dx*(ipoint-1)+x_left
    call Potential(xx,hh)
    V_dvr(ipoint,:,:)=hh
  end do

  V_absorb=V_dvr
  do ipoint=1,nabsorb_left
    do ilevel=1,nlevel
      V_absorb(nabsorb_left-ipoint+1,ilevel,ilevel)=V_absorb(nabsorb_left-ipoint+1,ilevel,ilevel)+gama*real(nabsorb_left-ipoint+1)/real(nabsorb_left)
    end do
  end do
  do ipoint=1,nabsorb_right
    do ilevel=1,nlevel
      V_absorb(npoint-nabsorb_right+ipoint,ilevel,ilevel)=V_absorb(npoint-nabsorb_right+ipoint,ilevel,ilevel)+gama*real(ipoint)/real(nabsorb_right)
    end do
  end do
  
  H_dvr=0.0d0
  H_absorb=0.0d0
  do ipoint=1,npoint
  do ilevel=1,nlevel
  do jpoint=1,npoint
  do jlevel=1,nlevel
    if(ilevel.eq.jlevel) then
      if(ipoint.eq.jpoint) then
        H_dvr(ipoint+(ilevel-1)*npoint,jpoint+(jlevel-1)*npoint)=H_dvr(ipoint+(ilevel-1)*npoint,jpoint+(jlevel-1)*npoint)+0.5d0/dx**2*pi2/3.0d0/mass
        H_absorb(ipoint+(ilevel-1)*npoint,jpoint+(jlevel-1)*npoint)=H_absorb(ipoint+(ilevel-1)*npoint,jpoint+(jlevel-1)*npoint)+0.5d0/dx**2*pi2/3.0d0/mass
      else
        if(mod(ipoint-jpoint,2).eq.0) then
          H_dvr(ipoint+(ilevel-1)*npoint,jpoint+(jlevel-1)*npoint)=H_dvr(ipoint+(ilevel-1)*npoint,jpoint+(jlevel-1)*npoint)+0.5d0/dx**2*2.0d0/(ipoint-jpoint)**2/mass
          H_absorb(ipoint+(ilevel-1)*npoint,jpoint+(jlevel-1)*npoint)=H_absorb(ipoint+(ilevel-1)*npoint,jpoint+(jlevel-1)*npoint)+0.5d0/dx**2*2.0d0/(ipoint-jpoint)**2/mass
        else
          H_dvr(ipoint+(ilevel-1)*npoint,jpoint+(jlevel-1)*npoint)=H_dvr(ipoint+(ilevel-1)*npoint,jpoint+(jlevel-1)*npoint)-0.5d0/dx**2*2.0d0/(ipoint-jpoint)**2/mass
          H_absorb(ipoint+(ilevel-1)*npoint,jpoint+(jlevel-1)*npoint)=H_absorb(ipoint+(ilevel-1)*npoint,jpoint+(jlevel-1)*npoint)-0.5d0/dx**2*2.0d0/(ipoint-jpoint)**2/mass
        endif
      endif
    endif

    if(ipoint.eq.jpoint) then
      H_dvr(ipoint+(ilevel-1)*npoint,jpoint+(jlevel-1)*npoint)=H_dvr(ipoint+(ilevel-1)*npoint,jpoint+(jlevel-1)*npoint)+V_dvr(ipoint,ilevel,jlevel)
      H_absorb(ipoint+(ilevel-1)*npoint,jpoint+(jlevel-1)*npoint)=H_absorb(ipoint+(ilevel-1)*npoint,jpoint+(jlevel-1)*npoint)+V_absorb(ipoint,ilevel,jlevel)
    endif
  enddo
  enddo
  enddo
  enddo
  
  H_dvr=H_dvr*cmplx(0.0d0,-1.0d0)
  H_absorb=H_absorb*cmplx(0.0d0,-1.0d0)

end subroutine

subroutine trstat
  implicit none
  integer ii
  integer flag
  real*8 r_cumulate(nsnap,nlevel)
  real*8 t_cumulate(nsnap,nlevel)

  do ipoint=1,npoint
  do ilevel=1,nlevel
  if (ipoint.le.ncenter) then
    reflect(ilevel,istep)=reflect(ilevel,istep)+abs(c_absorb(ipoint+(ilevel-1)*npoint))**2-abs(c_dvr(ipoint+(ilevel-1)*npoint))**2
  else
    trans(ilevel,istep)=trans(ilevel,istep)+abs(c_absorb(ipoint+(ilevel-1)*npoint))**2-abs(c_dvr(ipoint+(ilevel-1)*npoint))**2
  end if
  end do
  end do

  if (istep+1.eq.nstep) then
    do ii=1,nstep
      r_snap(isnap,:)=r_snap(isnap,:)+reflect(ii,:)
      t_snap(isnap,:)=t_snap(isnap,:)+trans(ii,:)
    end do
    do ii=1,isnap
      r_cumulate(isnap,:)=r_cumulate(isnap,:)+r_snap(ii,:)
      t_cumulate(isnap,:)=t_cumulate(isnap,:)+t_snap(ii,:)
    end do
    write(44,'(9F20.10)') t_cumulate(isnap,:),r_cumulate(isnap,:),t_snap(isnap,:),r_snap(isnap,:) 
    reflect=0.0d0
    trans=0.0d0
  end if    
  

end subroutine

subroutine rk4_wave
  implicit none
  real*8 dt2
  
  complex*16 c_inter(npoint*nlevel)
  complex*16 cc1(npoint*nlevel)
  complex*16 cc2(npoint*nlevel)
  complex*16 cc3(npoint*nlevel)
  complex*16 dc1(npoint*nlevel)
  complex*16 dc2(npoint*nlevel)
  complex*16 dc3(npoint*nlevel)
  complex*16 dc4(npoint*nlevel)

  dt2=dt*0.5d0
  c_inter=c_absorb

  dc1=matmul(H_dvr,c_inter)
  cc1=c_inter+dc1*dt2
  dc2=matmul(H_dvr,cc1)
  cc2=c_inter+dc2*dt2
  dc3=matmul(H_dvr,cc2)
  cc3=c_inter+dc3*dt
  dc4=matmul(H_dvr,cc3)
  c_dvr=c_inter+1.0d0/3.0d0*(dc1*dt2+dc2*dt+dc3*dt+dc4*dt2)

  dc1=matmul(H_absorb,c_inter)
  cc1=c_inter+dc1*dt2
  dc2=matmul(H_absorb,cc1)
  cc2=c_inter+dc2*dt2
  dc3=matmul(H_absorb,cc2)
  cc3=c_inter+dc3*dt
  dc4=matmul(H_absorb,cc3)
  c_absorb=c_inter+1.0d0/3.0d0*(dc1*dt2+dc2*dt+dc3*dt+dc4*dt2)

  call trstat

end subroutine

end module
