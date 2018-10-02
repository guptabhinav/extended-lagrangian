!TO test blue moon on a double well along x and harmonic along y
program mdcode

!USE data_definitions
!USE nose_hoover_chains

implicit none
real*8 :: ke0,rand,gauss,sigma
integer :: i, ii

integer :: steps_max,step, pfrq
real*8  :: dt,u,u1,u2,te,ke,pe
real*8  :: dummy,t0,t_inst, avg_t
real*8  :: umbr_mean,umbr_k
logical :: us,ext, lang
real*8  :: t0_s, t_s_inst,ke_s,pe_s,k_s,avg_t_s, ddt

integer :: j,ndof,ndof_s, ios, ios2, imts, nmts
real*8  :: glang, glang_s, kbt0, kbt0_s

integer :: meta_max,freq_meta,nmtd,meta_cv
real*8  :: width,h0,cvar,cv_ist,alpha,wtdt,gausspot, avmf
logical :: meta,hilladd,restart,mts

real*8 :: lang_params(6), lang_params_s(6)

real*8, allocatable :: x(:),v(:),mass(:)
real*8,allocatable  :: f(:)
real*8,allocatable  :: height(:),cv_history(:)
real*8, allocatable ::x_s(:),v_s(:),mass_s(:),f_harm(:), f_s(:), rnd(:), rnd_s(:), &
                      v_old(:), v_s_old(:)
real*8, parameter :: kb=3.16e-6  !Boltman constant in a.u. K^-1
real*8, parameter :: amu_to_au=1822.d0
real*8, parameter :: pi=4.d0*atan(1.d0)
integer, parameter :: nd=1, ns=1, ncv=1            !number of dimensions
character(len=20) :: myfilestat, myfilepos

mts=.false.
nmts=0
print *, "pi is =", pi
print *, "kb in HK^-1 =", kb
print *, "amu to au=", amu_to_au
!type(nhc_data)     :: nhc_atoms
!type(nhc_data)     :: nhc_cv

!open files
open( unit =20, file = 'input',status='unknown')
!=====================================================


avg_t=0.0d0
avg_t_s=0.0d0

!=========================
!Langevin parameters
!lang=.false.
!===========================


!===========================
!Nose-hoover parameters
!cdata%nhc_atoms=.false.
!cdata%nhnc_atoms=10
!cdata%nsuzuki_atoms=5
!cdata%nmultint_atoms=5
!cdata%nosefrq_atoms=3900.d0 !in cm-1
!
!cvdata%nhc_cv=.false.
!cvdata%nhnc_cv=10
!cvdata%nsuzuki_cv=5
!cvdata%nmultint_cv=5
!cvdata%nosefrq_cv=8000.d0 !in cm-1
!==================================

meta_cv=2         !The dimension on which MTD bias is acting

allocate(x(nd),v(nd))
allocate(x_s(ns))
allocate(v_s(ns))
allocate(mass(nd))
allocate(mass_s(ns))
allocate(f(nd))
allocate(f_s(ns))
allocate(f_harm(ns))
allocate(rnd(nd))
allocate(rnd_s(ns))
allocate(v_s_old(ns))
allocate(v_old(nd))


!==================================
read(*,*)x(1:nd)
read(*,*)mass(1:nd)
read(*,*)dt
read(*,*)t0
read(*,*)t0_s
read(*,*)lang
read(*,*)glang
read(*,*)glang_s
read(*,*)ext
read(*,*)us
read(*,*)meta
read(*,*)k_s !coupling constant of x and x_s
read(*,*)mass_s(1:ns)
read(*,*)umbr_k ! spring constant for umbrella sampling (acting on x_s)
read(*,*)umbr_mean !eq. of umbrella window (for x_s)
read(*,*)steps_max !md step number
read(*,*)freq_meta !mtd time step
read(*,*)h0
read(*,*)width
read(*,*)wtdt !deltaT parameter in energy unit (a.u.)
read(*,*)pfrq
read(*,*)restart
read(*,*)nmts

if(nmts.gt.0)mts=.true.


!............
write(*,*)'x=',x(1:nd)
write(*,*)'mass=',mass(1:nd)
write(*,*)'dt=',dt
write(*,*)'t0=',t0
write(*,*)'t0_s=',t0_s
if(lang)then 
  write(*,*)'lang is true'
else
  write(*,*)'lang is false'
end if
write(*,*)'glang=',glang
write(*,*)'glang_s=',glang_s
if(ext)then
  write(*,*)'ext is true'
else 
  write(*,*)'ext is false'
end if
if(us)then
  write(*,*)'us is true'
else 
  write(*,*)'us is false'
end if
if(meta)then
  write(*,*)'meta is true'
else 
  write(*,*)'meta is false'
end if
write(*,*)'k_s=',k_s
write(*,*)'mass_s=',mass_s(1:ns)
write(*,*)'umbr_k=',umbr_k
write(*,*)'umbr_mean=',umbr_mean
write(*,*)'steps_max=',steps_max
write(*,*)'freq_meta',freq_meta
write(*,*)'h0=',h0
write(*,*)'width=',width
write(*,*)'wtdt (a.u.)=',wtdt !deltaT parameter in energy unit (a.u.)
write(*,*)'wtdt (K)=',wtdt/kb !deltaT parameter in K is printed
write(*,*)'pfrq=',pfrq !printing frequency
if(restart)then 
  write(*,*)'restart is true'
else
  write(*,*)'restart is false'
end if
if(restart)then 
  write(*,*)'mts is true; nmts=',nmts
else
  write(*,*)'mts is false'
end if

if(restart)then
   myfilestat='old'
   myfilepos='append'
else
   myfilestat='unknown'
   myfilepos='asis'
end if
open( unit =21, file = 'TRAJ',status=trim(myfilestat),position=trim(myfilepos))
open( unit =22, file = 'ENERGIES',status=trim(myfilestat),position=trim(myfilepos))
open( unit =23, file = 'CV_VAL',status=trim(myfilestat),position=trim(myfilepos))
open( unit =24, file = 'AVG_TEMP',status=trim(myfilestat),position=trim(myfilepos))
open( unit =25, file = 'MTD',status=trim(myfilestat),position=trim(myfilepos))

open(unit=88,file='meanforce.dat')


mass(1:nd)=mass(1:nd)*amu_to_au !mass of particle     
mass_s(1:ns)=mass_s(1:ns)*amu_to_au
!==================================


!time_max=20000   !max steps md  
step=0          !strting step md


ndof=nd*1
ndof_s=ns*1

kbt0=kb*t0
kbt0_s=kb*t0_s

ke0=0.5d0*kb*t0*dfloat(nd)
ke_s=0.5d0*kb*t0_s*dfloat(ns)

nmtd=0

if(restart)then
  open(55,file='restart',status='old')
  read(55,*)v(1:nd),v_s(1:nd)
  read(55,*)x(1:ns),x_s(1:ns)
  read(55,*)step
  close(55)
  print *, 'finished reading restart file'
  close(55)
  print *, 'finished reading restart file'
!  open(55,file='MTD',status='old',IOSTAT=ios)
  rewind(25)
!  if(ios.eq.0)then
    nmtd=0
    do
      nmtd=nmtd+1
      if(nmtd.le.meta_max)&
       read(25,*,iostat=ios2) i,j,cv_history(nmtd),height(nmtd),gausspot !*alpha
      if(ios2.ne.0)exit
      if(nmtd.gt.meta_max)stop 'something is wrong...nmtd>metamax'
    end do
    nmtd=nmtd-1
    print *, 'finished reading MTD file for restart: #hills=',nmtd
    close(25)
    open( unit =25, file = 'MTD',status='old',position='APPEND')
else
!initial velocities for the physical degrees of freedom
   print *, 'initial velocity to assign 1'
   do i=1,nd
     call random_number(u1)
     call random_number(u2)
     v(i)= dsqrt(-2.d0*dlog(u1))*dcos(2.d0*pi*u2)*dsqrt(2.d0*ke0/mass(i)/dfloat(nd))
   end do
   call v_scal(nd,v,ke0,mass)
   !initial velocities for the auxiliary degrees of freedom
   print *, 'initial velocity to assign 2'
   do i=1,ns
     call random_number(u1)
     call random_number(u2) 
     v_s(i)=dsqrt(-2.d0*dlog(u1))*dcos(2.d0*pi*u2)*dsqrt(2.d0*ke_s/mass_s(i)/dfloat(ns))
   end do
   call v_scal(ncv,v_s,ke_s,mass_s)
   write(*,*) 'initial velocity done'

!setting initial position of auxiliary variables to physical variables
   if(ns.ne.nd)stop 'ns .ne. nd not implemented yet'
   do i=1,ns
     do j=1,nd
       if (i.eq.j) x_s(i)=x(j)
     end do
   end do
end if
print *, 'initialize done '

!get the cv value
!cvar: s_mtd(t)
!cv_ist: S_mtd(t)
call cv_value(ns,meta_cv,x_s,cvar) ! cvar=x_s(2)
call cv_value(nd,meta_cv,x,cv_ist) ! cv_ist=x(2)
print *, 'initial cvvalues done '

meta_max=steps_max/freq_meta+1   !max mtd steps
print *, 'meta_max =', meta_max
if(meta_max.le.0)stop 'error meta_max <=0'

allocate(height(meta_max))
allocate(cv_history(meta_max))

!get initial gaussian position
if(nmtd.eq.0)then 
  cv_history(1)=cvar
  height(1)=h0
  nmtd=1
end if

!========================================================

print *, 'initial force'
f_s=0.d0
pe_s=0.d0
if (ext) call force_ext(nd,ns,x,x_s,pe_s,f_s,k_s,nmtd,height,width,cvar,cv_history, &
                        meta,f_harm,h0,wtdt,gausspot,us,umbr_k,umbr_mean)
call force(nd,x,pe,f,ncv,ext,f_harm,pe_s,ns)

avmf=0.d0

print *, 'starting MD'
md_loop : do
step= step + 1
ke=0.d0
te=0.d0
pe=0.d0

if(lang)then

!  v_old(1:2)=v(1:2)
!  v_s_old(1:2)=v_s(1:2)
!  call lang_v_update1(nd,dt,kbt0,mass,v,f,rnd,lang_params)
!  call lang_v_update1(ns,dt,kbt0_s,mass_s,v_s,f_s,rnd_s,lang_params_s)

call lang_v_update1_new3(nd,dt,mass,kbt0,x,v,f,rnd,glang)
call lang_x_update_new3(nd,dt,mass,kbt0,x,v,f,rnd,glang)


!if(ext.and.mts)then
!  do imts=1,nmts
!    ddt=dt/dfloat(nmts)
!    call lang_x_update_new(ns,ddt,mass_s,kbt0_s,x_s,v_s,f_s,rnd_s,glang_s)
!    call cv_value(ns,meta_cv,x_s,cvar)
!    call lang_v_update1_new(ns,ddt,mass_s,kbt0_s,x_s,v_s,f_s,rnd_s,glang_s)
!    call force_ext(nd,ns,x,x_s,pe_s,f_s,k_s,nmtd,height,width,cvar,&
!                  cv_history,meta,f_harm,h0,wtdt,gausspot,us,umbr_k,umbr_mean)
!    
!    call lang_v_update2_new(ns,ddt,mass_s,kbt0_s,x_s,v_s,f_s,rnd_s,glang_s)
!  end do
!else if(ext) then
   if(ext)then
      call lang_v_update1_new3(ns,dt,mass_s,kbt0_s,x_s,v_s,f_s,rnd_s,glang_s)
      call lang_x_update_new3 (ns,dt,mass_s,kbt0_s,x_s,v_s,f_s,rnd_s,glang_s)
  end if
!end if
!  call lang_x_update(nd,dt,mass,kbt0,x,v_old,f,rnd,lang_params)
!  call lang_x_update(ns,dt,mass_s,kbt0_s,x_s,v_s_old,f_s,rnd_s,lang_params_s)
else
   call verlet_x_update(nd,dt,mass,x,v,f)
   if(ext)call verlet_x_update(ns,dt,mass_s,x_s,v_s,f_s)
   call verlet_v_update(nd,dt,mass,x,v,f)
   if(ext)call verlet_v_update(ns,dt,mass_s,x_s,v_s,f_s)
end if

call cv_value(ns,meta_cv,x_s,cvar)
call cv_value(nd,meta_cv,x,cv_ist)

if(mod(step,freq_meta).eq.0) then
  nmtd=nmtd+1 !nmtd = current metadynamcis step
  cv_history(nmtd)=cvar
  height(nmtd)= h0*dexp(-gausspot/wtdt)
  print *, "updating hills at ", step, " Hill #=", nmtd, " Height =", height(nmtd)
  write(25,'(2I16,3F16.6)') nmtd,step,cv_history(nmtd),height(nmtd),gausspot !*alpha
end if


if (ext)call force_ext(nd,ns,x,x_s,pe_s,f_s,k_s,nmtd,height,width,cvar,&
                        cv_history,meta,f_harm,h0,wtdt,gausspot,us,umbr_k,umbr_mean)
call force(nd,x,pe,f,ncv,ext,f_harm,pe_s,ns)

if(lang)then
!  call lang_v_update2(nd,dt,mass,v,f,lang_params,glang)
!  call lang_v_update2(ns,dt,mass_s,v_s,f_s,lang_params_s,glang_s)
  call lang_v_update2_new3(nd,dt,mass,kbt0,x,v,f,rnd,glang)
  if(ext)call lang_v_update2_new3(ns,dt,mass_s,kbt0_s,x_s,v_s,f_s,rnd_s,glang_s)
else 
   call verlet_v_update(nd,dt,mass,x,v,f)
   if(ext)call verlet_v_update(ns,dt,mass_s,x_s,v_s,f_s)
end if

!   call lang_new_2(ns,dt,glang_s,mass_s,t0_s,x_s,v_s,f_s)
!   call lang_new_2(nd,dt,glang,mass,t0,x,v,f)
!else
!   call vel(ns,dt,v_s,f_s,mass_s)
!   call vel(nd,dt,v,f,mass)
!   if(nhc_atoms%nhchain)then
!     call nose_evolve(nd,ndof,t0,mass,v,nhc_atoms)
!     call nose_energy(nd,ndof,t0,nhc_atoms)
!   end if
!   if(nhc_cv%nhchain)then
!     call nose_evolve(ns,ndof_s,t0_s,mass_s,v_s,nhc_cv)
!     call nose_energy(ns,ndof_s,t0_s,nhc_cv)
!   end if
!end if

!calculate total energy and temperature
ke=0.0d0
do i=1,nd 
   ke=ke+0.5d0*mass(i)*v(i)*v(i)
end do
t_inst=2.d0*ke/(dfloat(nd)*kb)


ke_s=0.0d0
if(ext) then
  do i=1,ns 
    ke_s=ke_s+0.5d0*mass_s(i)*v_s(i)*v_s(i)
  end do
  t_s_inst=2.d0*ke_s/(dfloat(ns)*kb)
  ke=ke+ke_s
end if

te=ke+pe !pe has pe_s included (inside force routine)

avg_t=avg_t+t_inst
avg_t_s=avg_t_s+t_s_inst

!================================================
avmf=avmf-umbr_k*(x_s(1)-umbr_mean)
!================================================

if(mod(step,pfrq).eq.0) then
  write(*,'(I16,6F16.6)')step,t_inst,t_s_inst,pe,ke,te,gausspot
  write(24,*)step, avg_t/real(step), avg_t_s/real(step)
  write(22,'(I16,5F16.6)')step,t_inst,t_s_inst,pe,ke,te
  write(21,'(I16,2F16.8)')step,x(1:nd),v(1:nd) !,v(1:nd)
  write(23,'(I16,1F16.4)')step,x_s(1:ns)!,v_s(1:ns)
!printing mean force (including the exact one)
  write(88,*)step,umbr_k*(x_s(1)-umbr_mean),avmf/dfloat(step), 0.02d0*4.0d0*umbr_mean*(umbr_mean**2-1.d0**2)
  open(55,file='restart')
   write(55,'(2e16.8)')v(1:nd),v_s(1:nd)
   write(55,'(4e16.8)')x(1:ns),x_s(1:ns)
   write(55,*)step
  close(55)
  print *, 'wrote retart file'
end if

!if(step.eq.2)stop 'stopping'

if (step.ge.steps_max) exit md_loop
end do md_loop

  open(55,file='restart',status='old')
   write(55,'(4e16.8)')v(1:2),v_s(1:2)
   write(55,'(4e16.8)')x(1:2),x_s(1:2)
   write(55,*)step
  close(55)
  print *, 'wrote retart file'

end program mdcode

!===============================================================================
subroutine v_scal(nd,v,ke_init,mass)
implicit none
integer :: nd
real*8 :: v(nd),ke_init,mass(nd)
integer :: i
real*8 :: ke,scal

ke=0.0d0
!write(*,*)'ke0',ke_init
do i=1,nd
ke=ke+0.5d0*mass(i)*v(i)**2
end do
scal=dsqrt(ke_init/ke)
do i=1,nd
v(i)=v(i)*scal
end do
return
end subroutine 
!=====================================================================
subroutine pos(nd,dt,x,v,f,mass)
implicit none
integer :: nd
real*8 :: dt,x(nd),v(nd),f(nd),mass(nd)
integer:: j

do j=1,nd
  x(j)=x(j)+dt*v(j)+0.5d0*dt*dt*f(j)/mass(j)
end do
return
end subroutine
!=====================================================================
subroutine vel(nd,dt,v,f,mass)
implicit none
integer :: nd
real*8 :: dt,v(nd),f(nd),mass(nd)
integer:: j

do j=1,nd
  v(j)=v(j)+0.5d0*dt*f(j)/mass(j)
end do
return
end subroutine
!=====================================================================
subroutine force_ext(nd,ns,x,x_s,pe_s,f_s,k_s,nmtd,height,width,cvar,cv_history,meta,f_harm,h0,wtdt,gausspot,us,umbr_k,umbr_mean)
implicit none
integer :: nd,ns
real*8  :: x(nd),f(nd),pe_s,k_s,f_s(ns),x_s(ns), gausspot, umbr_k, umbr_mean
integer :: i, j, it


real*8  :: expn,cv_history(*)
integer :: meta_cv, nmtd

real*8  :: cvar,diff,width,height(*)
logical :: meta,us

real*8  :: f_harm(ns),h0,wtdt,step_gauss
integer :: umb_cv

pe_s=0.d0
f_s=0.d0
do i=1,ns
  do j=1,nd
     if (i.eq.j) then
       f_s(i)=f_s(i)-k_s*(x_s(i)-x(j)) !v=-kx as f=0.5*k*x**2
       pe_s=pe_s+0.5d0*k_s*(x_s(i)-x(j))**2
    end if
   end do
end do

do i=1,ns
  f_harm(i)=f_s(i)
end do

!force contribution from umbrella.
if(us) then
  umb_cv=1
  do i=1,nd
     if (i.eq.umb_cv) then
       f_s(i)=f_s(i)-umbr_k*(x_s(i)-umbr_mean) !v=-kx as f=0.5*k*x**2
       pe_s=pe_s+0.5d0*umbr_k*(x_s(i)-umbr_mean)**2
     end if
  end do
end if

!force contribution from bias.
!meta=.false.
meta_cv=2
gausspot=0.d0
if(meta) then
  do it= 1,nmtd-1
     diff=(cvar-cv_history(it))
     step_gauss=height(it)*dexp(-0.5d0*(diff/width)**2)
     gausspot= gausspot+step_gauss
     do i=1,nd
       if (i.eq.meta_cv) then
         f_s(i)=f_s(i)+(diff/(width)**2.0d0)*step_gauss
       end if 
     end do
  end do
  pe_s=pe_s+gausspot
!TODO
!  height(nmtd)= h0*dexp(-gausspot/wtdt)
end if

end subroutine    
!==========================================================================================
subroutine force(nd,x,pe,f,ncv,ext,f_harm,pe_s,ns)
implicit none
integer :: nd,ncv,ns
real*8  :: x(nd),f(nd),pe,dummy,umbr_mean,umbr_k,b
real*8  :: f_s(ns),pe_s,a1,a2,a3,a4,w1,w2,w3,w4,p1,p2,p3,p4
real*8  :: x1(3),x2(3),x3(3),x4(3)
integer :: i, icv,j,pt,umb_cv,it
logical :: us, ext

!real*8  :: expn,gausspot,cv_history(10000000)!n_meta)
!integer :: meta_cv, nmtd

!real*8  :: cvar,diff,width,height(10000000)!n_meta)
!logical :: meta

real*8  :: f_harm(*)

real*8, parameter :: v0=0.02d0, a=1.d0, v1=0.1d0 !in atomic units

!pt=4 !potential type
!if(pt.eq.1) then !single well potential
!     pe=0.0d0
!     a=0.3d0
!    do i=1,nd
!     pe=pe+(0.5*a*(x(i)-1.0d0)**2) !0.5*k*x^2
!    end do
!
!    do  i=1,nd
!     f(i)=-a*(x(i)-1.0d0)
!    end do
!end if

f(1) = -v0*4.0d0*x(1)*(x(1)**2-a**2)/a**4 
!f(2) = -v0*4.0d0*x(2)*(x(2)**2-a**2)/a**4 
!f(1) = -v1*(x(1)-0.d0)
!f(2) = -v1*(x(2)-0.d0)
pe=v0*(x(1)**2-a**2)**2 !+ v0*(x(2)**2-a**2)**2
!pe=0.5d0*v1*x(1)**2 + 0.5d0*v1*x(2)**2


!force contribution due to extended part
if(ext) then
  do i=1,nd
     do j=1,ns
       if (i.eq.j) then
         f(i)=f(i)-f_harm(j) !v=-kx as f=0.5*k*x**2
       end if
     end do
  end do
  pe=pe+pe_s
end if

!force contribution from bias.
!meta=.false.
!meta_cv=1
!if( meta .eqv. .true. ) then
!gausspot=0.d0
!do it= 1,nmtd-1
!     expn = 0.0d0
!     diff = (cvar-cv_history(it))
!     gausspot = height(it)*dexp(-0.5d0*(diff/width)**2)
!
!     pe=pe+gausspot
!
!     do i=1,nd
!     if (i.eq.meta_cv) then
!     f(i) = f(i) +(diff/(width)**2.0d0)*gausspot
!     end if 
!     end do
!
!end do
!write(26,*) 'gausspot', gausspot
!end if
!
!
end subroutine    

!===========================================================================
subroutine anderson_therm(nd,dt,nu,v,sigma)
implicit none
integer :: nd,i
real*8  :: u,rand,dt,nu,y,v(nd),sigma
y=nu*dt
do i=1,nd
call random_number(u)
if (u.le.y) then
call gauss(rand,sigma)
v(i)=rand
end if
end do
return
end subroutine
!===========================================================================
subroutine gauss(rand,sigma)
implicit none
real*8 :: x1,x2,rand,u1,u2,sigma
real*8, parameter :: pi=4.d0*atan(1.d0)

call random_number(x1)
call random_number(x2)
rand= sigma*sqrt(-2.0*log(x1))*cos(2.d0*pi*x2)
return
end subroutine
!===========================================================================

subroutine lang_new_1(n,dt,glang,amass,T0ions,rxyz,vxyz,fxyz)
implicit none

integer :: n
real*8 :: T0ions,dt, glang, amass(n), rxyz(n), vxyz(n), fxyz(n)

integer :: iat, k
real*8  :: gfric, noise
real(kind=8) :: pi=4.d0*atan(1.d0) 
real*8, parameter :: kb=3.16e-6  !Boltman constant in a.u. K^-1
!au_to_k=315774.664550534774D0 

gfric=(1.d0-0.5d0*dt*glang)

  do k=1,n
    rxyz(k)=rxyz(k)+dt*vxyz(k)+0.5*dt**2/amass(k)*fxyz(k)
    vxyz(k)=vxyz(k)*gfric+0.5*dt/amass(k)*fxyz(k)
  end do
end subroutine lang_new_1
!===============================================================================

subroutine lang_new_2(n,dt,glang,amass,T0ions,rxyz,vxyz,fxyz)
implicit none

integer :: n
real*8 :: T0ions,dt, glang, amass(n), rxyz(n), vxyz(n), fxyz(n)

integer :: iat, k
real*8  :: gfric, noise
real(kind=8) :: pi=4.d0*atan(1.d0)!, au_to_k=315774.664550534774D0 
real*8, parameter :: kb=3.16e-6  !Boltman constant in a.u. K^-1

real*8 :: dd1
real*8 :: dd2
real*8 :: frand(n)

gfric=(1.d0-0.5d0*dt*glang)
do k=1,n
  CALL RANDOM_NUMBER(dd1)
  CALL RANDOM_NUMBER(dd2)
  if(mod(k,2).eq.0)then
     dd1= dsqrt(-2.d0*dlog(dd1))*dcos(2.d0*pi*dd2)
   else
     dd1= dsqrt(-2.d0*dlog(dd1))*dsin(2.d0*pi*dd2)
   end if 
!    dd1=dd1-0.5d0
!   noise=dsqrt(6.d0*glang*T0ions*amass(k)/dt)
   noise=dsqrt(2.d0*glang*kb*T0ions*amass(k)/dt)
   fxyz(k)=fxyz(k)+2.d0*noise*dd1
end do

do k=1,n
  vxyz(k)=vxyz(k)*gfric+0.5*dt/amass(k)*fxyz(k)
end do

END SUBROUTINE lang_new_2
!============================================================================================

subroutine cv_value(nd,meta_cv,x,cvar)
implicit none
integer :: i,icv,nd,meta_cv
real*8  :: cvar,x(*)
do i=1,nd
  if (i.eq.meta_cv) cvar = x(i)
end do
return
end subroutine
!============================================================================================

!call cv_value(ns,meta_cv,x_s,cvar)
!call cv_value(nd,meta_cv,x,cv_ist)

subroutine lang_cons(dt,lgamma,lang_facts)
implicit none
!=======--------------------------------------------------------------==
!     SAGARMOY MANDAL  !SAGAR HACK
!---------------------------------------------------------------------!
!     CHEMICAL PHYSICS 236(1998) 243-252                              !
!                                                                     !
!     VALUES OF THE COEFFICIENTS CAN BE FOUND IN THE ABOVE PAPER      !
!---------------------------------------------------------------------!
!      IMPLICIT NONE
!     Local Variables
      REAL*8 ::      lgamma, GAMMA_DT, dt, lang_facts(6)
      REAL*8 :: c_0, c_1, c_2, sigma_v, sigma_r, c_rv
!     ==--------------------------------------------------------------==
      GAMMA_DT = lGAMMA*DT
      c_0=dexp(-gamma_dt)
      c_1=(1.d0-c_0)/gamma_dt
      c_2=(1.d0-c_1)/gamma_dt
      sigma_v=dsqrt(1.d0-dexp(-2.d0*gamma_dt))
      sigma_r=dsqrt(dt*dt*(2.d0-(3.d0-4.d0*dexp(-gamma_dt)+dexp(-2.d0*gamma_dt))/gamma_dt)/gamma_dt)
      c_rv=dt*(1.d0-dexp(-gamma_dt))**2/gamma_dt/sigma_r/sigma_v
      !
      lang_facts(1)=c_0
      lang_facts(2)=c_1
      lang_facts(3)=c_2
      lang_facts(4)=sigma_v
      lang_facts(5)=sigma_r
      lang_facts(6)=c_rv
!=======--------------------------------------------------------------==
      RETURN
end subroutine lang_cons


  SUBROUTINE lang_v_update1(nd,dt,kbt,mass,velp,fion,rnd,lang_facts)
    implicit none
  !call lang_v_update1(v,f,rnd)
!=======--------------------------------------------------------------==
    integer :: nd
    REAL*8                   :: velp(nd),fion(nd)&
                                     ,rnd(nd), lang_facts(6), mass(nd), &
                                      dt, kbt
!   local variables
    INTEGER                        :: i, ia, is
    REAL*8                  :: fact
    REAL*8                  :: FACT1, gauss_dist1(2)
    REAL*8                  :: FACT2, FACT3, FACT4, FACT5

   real*8 :: c_0, c_1, c_2, sigma_v, sigma_r, c_rv
   real*8 :: z1, z2, u1, u2
   real*8, parameter :: pi=4.d0*atan(1.d0)

!---------------------------------------------------------------------
      c_0=lang_facts(1)
      c_1=lang_facts(2)
      c_2=lang_facts(3)
      sigma_v=lang_facts(4)
      sigma_r=lang_facts(5)
      c_rv=lang_facts(6)
      !
      !
      DO I=1,nd
        call random_number(u1)
        call random_number(u2) 
        z1=dsqrt(-2.d0*dlog(u1))*dcos(2.d0*pi*u2) 
        z2=dsqrt(-2.d0*dlog(u1))*dsin(2.d0*pi*u2) 
        rnd(i)=z1
        fact=sigma_v*dsqrt(kbt/mass(i))
        fact5=fact*(c_rv*z1+dsqrt(1.d0-c_rv*c_rv))*z2 
        velp(i)=C_0*VELP(i)+ C_0*C_2*FION(I)/C_1/MASS(I)*dt &
                     + FACT5
      ENDDO
!=======--------------------------------------------------------------==
    RETURN
  END SUBROUTINE lang_v_update1


subroutine gauss_dist( mu, sigma, dim, gauss_dist1)
implicit none
!-----------------------------------------------------------------------
!
! ... this function generates an array of numbers taken from a
! normal
! ... distribution of mean value \mu and variance \sigma
!
!IMPLICIT NONE
!
REAL*8, INTENT(IN)    :: mu
REAL*8, INTENT(IN)    :: sigma
INTEGER,  INTENT(IN)        :: dim
REAL*8               :: gauss_dist1( dim )
!local variables
REAL*8               :: x1, x2, w, y1 , y2
INTEGER  :: i
!
!
DO i = 1, dim, 2
  !
  gaussian_loop: DO
     !
     call random_number(y1)
     call random_number(y2)
     !
     x1 = 2.0D0 * y1 - 1.0D0
     x2 = 2.0D0 * y2 - 1.0D0
     !
     w = x1 * x1 + x2 * x2
     !
     IF ( w < 1.0D0 ) EXIT gaussian_loop
     !
  END DO gaussian_loop
  !
  w = dSQRT( ( - 2.0D0 * dLOG( w ) ) / w )
  !
  gauss_dist1(i) = x1 * w * sigma
  !
  IF ( i >= dim ) EXIT
  !
  gauss_dist1(i+1) = x2 * w * sigma
  !
END DO
!
gauss_dist1(1:dim) = gauss_dist1(1:dim) + mu
!
RETURN
!
END subroutine gauss_dist



  SUBROUTINE lang_v_update2(nd,dt,mass,velp,fion,lang_facts,gamma)
!=======--------------------------------------------------------------==
    implicit none
    integer :: nd
    REAL*8                     :: velp(nd), fion(nd), mass(nd), dt, lang_facts(6), gamma
!   local variables
    INTEGER                           :: i, ia, is
    REAL*8                     :: fact
    REAL*8                     :: FACT1
!ocl NOALIAS

    !$omp parallel do private(I,IS,IA,FACT) schedule(static)
!-----------------------------------------------------------
   real*8 :: c_0, c_1, c_2, sigma_v, sigma_r, c_rv

!---------------------------------------------------------------------
      c_0=lang_facts(1)
      c_1=lang_facts(2)
      c_2=lang_facts(3)
      sigma_v=lang_facts(4)
      sigma_r=lang_facts(5)
      c_rv=lang_facts(6)
      !
!      FACT1=2.D0*C_2
      FACT1=(1.d0-c_0/c_1)/gamma/dt
      DO I=1,nd
        VELP(i)=VELP(i)+FION(i)*FACT1*dt/mass(i)
      ENDDO
      !
!=======--------------------------------------------------------------==
    RETURN
  END SUBROUTINE lang_v_update2


  SUBROUTINE lang_x_update(nd,dt,mass,kbt,taup,velp,FION,RND,lang_facts)
  implicit none
  !call lang_x_update(nd,dt,mass,kbt0,x,v,f,rnd,lang_params)
! ==--------------------------------------------------------------==
! ==  UPDATE OF THE POSITIONS FOR VELOCITY VERLET                 ==
! ==--------------------------------------------------------------==
    integer :: nd
    REAL*8    :: taup(nd), &
                  velp(nd),fion(nd),rnd(nd), mass(nd), dt, lang_facts(6), kbt
    REAL*8    :: FACT,FACT1,FACT2,FACT3
    INTEGER   :: i, is

   real*8 :: c_0, c_1, c_2, sigma_v, sigma_r, c_rv

!---------------------------------------------------------------------
      c_0=lang_facts(1)
      c_1=lang_facts(2)
      c_2=lang_facts(3)
      sigma_v=lang_facts(4)
      sigma_r=lang_facts(5)
      c_rv=lang_facts(6)

      DO i=1,nd
        fact3=rnd(i)*sigma_r*dsqrt(kbt/mass(i))
        TAUP(i)=TAUP(i)+DT*VELP(I)*C_1 + &
                c_2*fion(i)*dt*dt/mass(i) + fact3
      ENDDO
!=======--------------------------------------------------------------==
    RETURN
  END SUBROUTINE lang_x_update


subroutine verlet_v_update(nd,dt,mass,x,v,f)
implicit none
integer :: nd
real*8 :: dt, mass(nd), x(nd), v(nd), f(nd)
integer :: i
do i=1,nd
  v(i)=v(i)+dt*0.5d0*f(i)/mass(i)
end do
end subroutine verlet_v_update
subroutine verlet_x_update(nd,dt,mass,x,v,f)
implicit none
integer :: nd
real*8 :: dt, mass(nd), x(nd), v(nd), f(nd)
integer :: i
do i=1,nd
  x(i)=x(i)+dt*v(i)+0.5d0*dt*dt*f(i)/mass(i)
end do
end subroutine verlet_x_update


SUBROUTINE lang_x_update_new(nd,dt,mass,kbt,x,v,f,rnd,gamma)
  implicit none

    integer :: nd
    REAL*8    :: x(nd),v(nd),f(nd),rnd(nd),mass(nd), dt, gamma, kbt
    REAL*8    :: FACT,z1,z2,u1,u2,fact2
    INTEGER   :: i
    real*8, parameter :: pi=4.d0*atan(1.d0)
    DO I=1,nd
        call random_number(u1)
        call random_number(u2) 
        z1=dsqrt(-2.d0*dlog(u1))*dcos(2.d0*pi*u2) 
        z2=dsqrt(-2.d0*dlog(u1))*dsin(2.d0*pi*u2) 
        fact=dsqrt(2.d0*kbt*gamma/mass(i))
        fact2=0.5d0*dt*dt*(f(i)-gamma*v(i)) + &
               fact*dt*dsqrt(dt)*(0.5d0*z1+0.5d0*z2/dsqrt(3.d0))
        !print *, 'fact2 =', fact2
        x(i)=x(i)+dt*v(i)+fact2
        rnd(i)=fact*dsqrt(dt)*z1-gamma*fact2
    END DO
END SUBROUTINE lang_x_update_new

SUBROUTINE lang_v_update1_new(nd,dt,mass,kbt,x,v,f,rnd,gamma)
  implicit none

    integer :: nd
    REAL*8    :: x(nd),v(nd),f(nd),rnd(nd),mass(nd), dt, gamma, kbt
    REAL*8    :: FACT,z1,z2,u1,u2,fact2
    INTEGER   :: i
    real*8, parameter :: pi=4.d0*atan(1.d0)
    DO I=1,nd
        v(i)=v(i)+0.5d0*dt*f(i)-dt*gamma*v(i)+rnd(i)
    END DO
END SUBROUTINE lang_v_update1_new

SUBROUTINE lang_v_update2_new(nd,dt,mass,kbt,x,v,f,rnd,gamma)
  implicit none

    integer :: nd
    REAL*8    :: x(nd),v(nd),f(nd),rnd(nd),mass(nd), dt, gamma, kbt
    REAL*8    :: FACT,z1,z2,u1,u2,fact2
    INTEGER   :: i
    real*8, parameter :: pi=4.d0*atan(1.d0)
    DO I=1,nd
        v(i)=v(i)+0.5d0*dt*f(i)
    END DO
END SUBROUTINE lang_v_update2_new

SUBROUTINE lang_x_update_new2(nd,dt,mass,kbt,x,v,f,rnd,gamma)
  implicit none

    integer :: nd
    REAL*8    :: x(nd),v(nd),f(nd),rnd(nd),mass(nd), dt, gamma, kbt
    REAL*8    :: FACT,z1,z2,u1,u2,fact2
    INTEGER   :: i
    real*8, parameter :: pi=4.d0*atan(1.d0)
    DO I=1,nd
        call random_number(u1)
        call random_number(u2) 
        z1=dsqrt(-2.d0*dlog(u1))*dcos(2.d0*pi*u2) 
        z2=dsqrt(-2.d0*dlog(u1))*dsin(2.d0*pi*u2) 
        fact=dsqrt(2.d0*kbt/gamma/mass(i))
        x(i)=x(i)+f(i)*dt/mass(i)/gamma+fact*z1*dsqrt(dt)
    END DO
END SUBROUTINE lang_x_update_new2

SUBROUTINE lang_x_update_new3(nd,dt,mass,kbt,x,v,f,rnd,gamma)
  implicit none

    integer :: nd
    REAL*8    :: x(nd),v(nd),f(nd),rnd(nd),mass(nd), dt, gamma, kbt
    REAL*8    :: FACT,z1,z2,u1,u2,fact2
    INTEGER   :: i
    real*8, parameter :: pi=4.d0*atan(1.d0)
    DO I=1,nd
        x(i)=x(i)+dt*v(i)
    END DO
END SUBROUTINE lang_x_update_new3

SUBROUTINE lang_v_update1_new3(nd,dt,mass,kbt,x,v,f,rnd,gamma)
  implicit none

    integer :: nd
    REAL*8    :: x(nd),v(nd),f(nd),rnd(nd),mass(nd), dt, gamma, kbt
    REAL*8    :: FACT,z1,z2,u1,u2,fact2
    INTEGER   :: i
    real*8, parameter :: pi=4.d0*atan(1.d0)
    DO I=1,nd
        call random_number(u1)
        call random_number(u2) 
        z1=dsqrt(-2.d0*dlog(u1))*dcos(2.d0*pi*u2) 
        z2=dsqrt(-2.d0*dlog(u1))*dsin(2.d0*pi*u2) 
        rnd(i)=z1
         fact=dt/(2.d0*mass(i))
        v(i)=v(i)+fact*f(i)
        v(i)=v(i)-mass(i)*v(i)*gamma &
          *fact+dsqrt(fact*2.d0*mass(i))*0.5d0*dsqrt(2.d0 &
          *kbt/mass(i)*gamma)*z1
    END DO
END SUBROUTINE lang_v_update1_new3

SUBROUTINE lang_v_update2_new3(nd,dt,mass,kbt,x,v,f,rnd,gamma)
  implicit none

    integer :: nd
    REAL*8    :: x(nd),v(nd),f(nd),rnd(nd),mass(nd), dt, gamma, kbt
    REAL*8    :: FACT,z1,z2,u1,u2,fact2
    INTEGER   :: i
    real*8, parameter :: pi=4.d0*atan(1.d0)
    DO I=1,nd
        fact=dt/(2.d0*mass(i))
        v(i)=v(i)+fact*f(i)
        v(i)=v(i)-mass(i)*v(i)*gamma &
          *fact+dsqrt(fact*2.d0*mass(i))*0.5d0*dsqrt(2.d0 &
          *kbt/mass(i)*gamma)*rnd(i)
    END DO
END SUBROUTINE lang_v_update2_new3
