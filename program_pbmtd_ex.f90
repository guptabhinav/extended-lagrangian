program mdcode
implicit none
real*8 :: ke0,rand,gauss,sigma
integer :: i, ii

integer :: steps_max,step, pfrq
real*8  :: dt,u,u1,u2,te,ke,pe,den
real*8  :: dummy,t0,t_inst, avg_t
real*8  :: umbr_mean,umbr_k
logical :: us,ext, lang
real*8  :: t0_s, t_s_inst,ke_s,pe_s,k_s,avg_t_s, ddt

integer :: j,ndof,ndof_s, ios, ios2, imts, nmts
real*8  :: glang, glang_s, kbt0, kbt0_s, pb_bias

integer :: meta_max,freq_meta,nmtd
real*8  :: alpha,wtdt, avmf
logical :: meta,hilladd,restart,mts

real*8 :: lang_params(6), lang_params_s(6)

real*8, allocatable :: x(:),v(:),mass(:), h0(:)
real*8, allocatable :: f(:), width(:), prob(:)
real*8, allocatable :: height(:,:), cv_history(:,:), gausspot(:)
real*8, allocatable :: x_s(:), v_s(:), mass_s(:), f_harm(:), f_s(:), rnd(:), rnd_s(:)
real*8, allocatable :: v_old(:), v_s_old(:), cvar(:), cv_ist(:)
integer, allocatable :: meta_cv(:)

real*8, parameter :: kb=3.16e-6                    !Boltman constant in a.u. K^-1
real*8, parameter :: amu_to_au=1822.d0
real*8, parameter :: pi=4.d0*atan(1.d0)
integer, parameter :: nd=2, ns=2, ncv=1            !number of dimensions, two well potential
character(len=20) :: myfilestat, myfilepos

mts=.false.
nmts=0
print *, "pi is =", pi
print *, "kb in HK^-1 =", kb
print *, "amu to au=", amu_to_au

open( unit =20, file = 'input',status='unknown')
!=====================================================


avg_t=0.0d0
avg_t_s=0.0d0

!=========================
!Langevin parameters
!lang=.false.
!===========================

allocate(meta_cv(ns))

meta_cv(1)=1
meta_cv(2)=2         !The dimension on which MTD bias is acting

allocate(h0(ns))
allocate(x(nd))
allocate(v(nd))
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
allocate(cvar(ns))
allocate(cv_ist(ns))
allocate(width(ns))
allocate(prob(ns))
allocate(gausspot(ns))

!====================inputs===========================================
read(20,*)x(1:nd)
read(20,*)mass(1:nd)
read(20,*)dt
read(20,*)t0
read(20,*)t0_s
read(20,*)lang
read(20,*)glang
read(20,*)glang_s
read(20,*)ext
read(20,*)us
read(20,*)meta
read(20,*)k_s                             !coupling constant of x and x_s
read(20,*)mass_s(1:ns)
read(20,*)umbr_k                          ! spring constant for umbrella sampling (acting on x_s)
read(20,*)umbr_mean                       !eq. of umbrella window (for x_s)
read(20,*)steps_max                       !md step number
read(20,*)freq_meta                       !mtd time step
read(20,*)h0(1:ns)
read(20,*)width(1:ns)
read(20,*)wtdt                            !deltaT parameter in energy unit (a.u.)
read(20,*)pfrq
read(20,*)restart
read(20,*)nmts                            ! restarting step
!=============================================================================================================

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
write(*,*)'glang=',glang                            ! lang coefficient
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
write(*,*)'h0=',h0(1:ns)
write(*,*)'width=',width(1:ns)
write(*,*)'wtdt (a.u.)=',wtdt                          !deltaT parameter in energy unit (a.u.)
write(*,*)'wtdt (K)=',wtdt/kb                          !deltaT parameter in K is printed
write(*,*)'pfrq=',pfrq                                 !printing frequency
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
open(30, file='data', status='replace')
open(unit=88,file='meanforce.dat')

meta_max=steps_max/freq_meta+1                         !max mtd steps
print *, 'meta_max =', meta_max
if(meta_max.le.0)stop 'error meta_max <=0'


mass(1:nd)=mass(1:nd)*amu_to_au                   !mass of particle     
mass_s(1:ns)=mass_s(1:ns)*amu_to_au
!==================================

allocate(height(ns,meta_max))
allocate(cv_history(ns,meta_max))

step=0                                            !starting step md


ndof=nd*1
ndof_s=ns*1

kbt0=kb*t0
kbt0_s=kb*t0_s

ke0=0.5d0*kb*t0*dfloat(nd)
ke_s=0.5d0*kb*t0_s*dfloat(ns)

nmtd=0

if(restart)then
  open(55,file='restart',status='old')
  read(55,*)v(1:nd),v_s(1:ns)
  read(55,*)x(1:nd),x_s(1:ns)
  read(55,*)step
  close(55)
  print *, 'finished reading restart file'
  close(55)
  print *, 'finished reading restart file'
  rewind(25)
    nmtd=0
    do
      nmtd=nmtd+1
      if(nmtd.le.meta_max)&
       read(25,*,iostat=ios2) i,j,cv_history(1:ns,nmtd),height(1:ns,nmtd),gausspot(1:ns)                  !*alpha
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

call cv_value(ns,meta_cv,x_s,cvar)                     ! cvar=x_s(2)
call cv_value(nd,meta_cv,x,cv_ist)                     ! cv_ist=x(2)
print *, 'initial cv values done '


!get initial gaussian position

if(nmtd.eq.0)then 
  cv_history(1:ns,1)=cvar(1:ns)
  height(1:ns,1)=h0(1:ns)
  nmtd=1
end if
print *, 'cv_history(1,1) = cvar(1)', cv_history(1,1)
print *, 'cv_history(2,1) = cvar(2)', cv_history(2,1)
!========================================================

print *, 'initial force'
f_s=0.d0
pe_s=0.d0
if (ext) call force_ext(nd,ns,x,x_s,pe_s,f_s,k_s,nmtd,height,width,cvar,cv_history,meta_cv, &
                        f,f_harm,h0,wtdt,gausspot,us,umbr_k,umbr_mean,kbt0_s,pb_bias,meta,den)
         call force(nd,x,pe,f,ncv,ext,f_harm,pe_s,ns)

avmf=0.d0

!pb_bias=0.d0

!=================== MD starts here ===============================

print *, 'starting MD'
md_loop : do
step= step + 1
ke=0.d0
te=0.d0
pe=0.d0


if(lang)then


call lang_v_update1_new3(nd,dt,mass,kbt0,x,v,f,rnd,glang)
call lang_x_update_new3(nd,dt,mass,kbt0,x,v,f,rnd,glang)


   if(ext)then
      call lang_v_update1_new3(ns,dt,mass_s,kbt0_s,x_s,v_s,f_s,rnd_s,glang_s)
      call lang_x_update_new3 (ns,dt,mass_s,kbt0_s,x_s,v_s,f_s,rnd_s,glang_s)
  end if
else
   call verlet_x_update(nd,dt,mass,x,v,f)
   if(ext)call verlet_x_update(ns,dt,mass_s,x_s,v_s,f_s)
   call verlet_v_update(nd,dt,mass,x,v,f)
   if(ext)call verlet_v_update(ns,dt,mass_s,x_s,v_s,f_s)
end if

! ==========why? ===========================
call cv_value(nd,meta_cv,x,cv_ist)
call cv_value(ns,meta_cv,x_s,cvar)
!===========================================

!===== well tempered metadynamics and pbmtd============================================

if(mod(step,freq_meta).eq.0) then

  nmtd=nmtd+1                !nmtd = current metadynamcis step

  cv_history(1:ns,nmtd)=cvar(1:ns)

  den = dexp(-gausspot(2)/kbt0_s) + dexp(-gausspot(1)/kbt0_s)
  
  prob(1) = dexp(-gausspot(1)/kbt0_s)/den
  prob(2) = dexp(-gausspot(2)/kbt0_s)/den

  write(100,*)  gausspot(1:ns), den


  height(1:ns,nmtd)= h0(1:ns)*dexp(-gausspot(1:ns)/wtdt)*prob(1:ns)                ! wt_mtd height

  print *, "prob of adding bias 1 and 2 =", prob(1:ns),  'height=', height(1:ns,nmtd)

  write(25,'(2I10,6F16.8)') step,nmtd,height(1:ns,nmtd),gausspot(1:ns),prob(1:ns)                  !*alpha

end if

!=============================================================================================

if (ext)call force_ext(nd,ns,x,x_s,pe_s,f_s,k_s,nmtd,height,width,cvar,cv_history,meta_cv, &
                      f,f_harm,h0,wtdt,gausspot,us,umbr_k,umbr_mean,kbt0_s,pb_bias,meta,den)
    
        call force(nd,x,pe,f,ncv,ext,f_harm,pe_s,ns)


if(lang)then
  call lang_v_update2_new3(nd,dt,mass,kbt0,x,v,f,rnd,glang)
  if(ext)call lang_v_update2_new3(ns,dt,mass_s,kbt0_s,x_s,v_s,f_s,rnd_s,glang_s)
else 
   call verlet_v_update(nd,dt,mass,x,v,f)
   if(ext)call verlet_v_update(ns,dt,mass_s,x_s,v_s,f_s)
end if

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

te=ke+pe                                              !pe has pe_s included (inside force routine)

avg_t=avg_t+t_inst
avg_t_s=avg_t_s+t_s_inst

!================================================
avmf=avmf-umbr_k*(x_s(1)-umbr_mean)
!================================================

if(mod(step,pfrq).eq.0) then
 
  write(30,'(I16,7F16.8)') step, x(1:nd), pe, ke, te, gausspot(1:ns)
 
  write(24,'(I16,4F16.8)') step, t_inst, avg_t/real(step), t_s_inst, avg_t_s/real(step)

  write(22,'(I16,5F16.8)') step, t_inst, t_s_inst, pe, ke, te

  write(21,'(I16,6F16.8)') step, x, v, pe , pb_bias                            

  write(23,'(I16,6F16.8)') step, x_s, v_s, pe_s, pb_bias                             

  write(88,*)step,umbr_k*(x_s(1)-umbr_mean),avmf/dfloat(step), 0.02d0*4.0d0*umbr_mean*(umbr_mean**2-1.d0**2)

  open(55,file='restart')
   write(55,'(4e16.8)')v(1:nd),v_s(1:ns)
   write(55,'(4e16.8)')x(1:nd),x_s(1:ns)
   write(55,*)step
  close(55)
!  print *, 'wrote retart file'
end if


if (step.ge.steps_max) exit md_loop
end do md_loop

!==================== MD finishes here ===============================================

  open(55,file='restart',status='old')
   write(55,'(4e16.8)')v(1:nd),v_s(1:ns)
   write(55,'(4e16.8)')x(1:nd),x_s(1:ns)
   write(55,*)step
  close(55)
  print *, 'wrote retart file'
print *, 'end program'


end program mdcode

!===============================================================================




subroutine v_scal(nd,v,ke_init,mass)
implicit none
integer :: nd
real*8 :: v(nd),ke_init,mass(nd)
integer :: i
real*8 :: ke,scal

ke=0.0d0
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
real*8 :: dt, v(*),f(*),mass(*)
integer:: j

do j=1,nd
  v(j)=v(j)+0.5d0*dt*f(j)/mass(j)
end do
return
end subroutine
!=====================================================================
subroutine force_ext(nd,ns,x,x_s,pe_s,f_s,k_s,nmtd,height,width,cvar,cv_history,meta_cv, &
                     f,f_harm,h0,wtdt,gausspot,us,umbr_k,umbr_mean,kbt0_s,pb_bias,meta,den)
implicit none
integer :: nd,ns
real*8 :: x(*), f(*), pe_s, k_s, f_s(ns), x_s(*)
real*8 :: gausspot(ns), umbr_k, umbr_mean, kbt0_s
real*8 :: expn, cv_history(2,*)
real*8 :: cvar(*), height(2,*), width(*)
real*8 :: f_harm(*), h0(*), wtdt, pb_bias,den

real*8, allocatable :: diff(:), step_gauss(:), dvds(:)

integer :: meta_cv(*), nmtd
integer :: i, j, it,is
integer :: umb_cv

logical :: meta, us

allocate(diff(ns))
allocate(step_gauss(ns))
allocate(dvds(ns))

pe_s=0.d0
f_s=0.d0

do i=1,ns
  do j=1,nd
     if (i.eq.j) then
       f_s(i)=f_s(i)-k_s*(x_s(i)-x(j))                     !v=-kx as f=0.5*k*x**2
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
       f_s(i)=f_s(i)-umbr_k*(x_s(i)-umbr_mean)               !v=-kx as f=0.5*k*x**2
       pe_s=pe_s+0.5d0*umbr_k*(x_s(i)-umbr_mean)**2
     end if
  end do
end if

!force contribution from bias.
!meta=.false.

!meta_cv1=1
!meta_cv2=2

gausspot(1:ns)=0.d0
dvds(1:ns)=0.d0

if(meta) then
  do is=1,ns
     do it= 1,nmtd-1
  
        diff(is)=(cvar(is)-cv_history(is,it))
!        diff(is) = diff(is)*diff(is) 
        step_gauss(is)=height(is,it)*dexp(-0.5d0*(diff(is)/width(is))**2)
 
        gausspot(is) = gausspot(is) + step_gauss(is)

        dvds(is) = dvds(is) + (diff(is)/(width(is))**2)*step_gauss(is)

     end do
                           !     gausspot(is) = gausspot(is) + step_gauss(is)
  end do

!======potential energy and forces in PBMTD ===================================

        den = dexp(-gausspot(1)/kbt0_s) + dexp(-gausspot(2)/kbt0_s) !+ dexp(dfloat(ns))

        f_s(1:2) =f_s(1:2) + dvds(1:2)*dexp(-gausspot(1:2)/kbt0_s)/den
    

        pb_bias =  -kbt0_s*dlog(den) + kbt0_s*dlog(dfloat(ns))

!       pe_s = pe_s - kbt0_s*dlog(pb_bias) !+ kbt0_s*dlog(dfloat(ns)) 
 
        pe_s = pe_s + pb_bias   

end if
!return
end subroutine    
!==========================================================================================

subroutine force(nd,x,pe,f,ncv,ext,f_harm,pe_s,ns)
implicit none
integer :: nd,ncv,ns
real*8  :: x(nd),f(nd),pe,dummy,umbr_mean,umbr_k
real*8  :: f_s(ns),pe_s,a1,a2,a3,a4,w1,w2,w3,w4,p1,p2,p3,p4
real*8  :: x1(3),x2(3),x3(3),x4(3)
integer :: i, icv,j,pt,umb_cv,it
logical :: us, ext
real*8  :: f_harm(*)
real*8, parameter :: v0=0.01d0, a=1.d0, b=1.d0                  !in atomic units

f(1) = -v0*4.0d0*x(1)*(x(1)**2-a**2)/a**4 
f(2) = -v0*4.0d0*x(2)*(x(2)**2-b**2)/b**4 

pe = v0*((x(1)**2-a**2)/a**2)**2 + v0*((x(2)**2-b**2)/b**2)**2                   !+ v0*(x(2)**2-a**2)**2


!force contribution due to extended part
if(ext) then
  do i=1,nd
     do j=1,ns
       if (i.eq.j) then
         f(i)=f(i)-f_harm(j)                                     ! v=-kx as f=0.5*k*x**2
       end if
     end do
  end do
  pe=pe+pe_s
end if

end subroutine    
!=======================================================

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



subroutine cv_value(nd,meta_cv,x,cvar)
implicit none
integer :: i,icv,nd,meta_cv(*)
real*8  :: cvar(*),x(*)
do i=1,nd
  if (i.eq.meta_cv(i)) cvar(i) = x(i)
end do
return
end subroutine
!============================================================================================


subroutine lang_cons(dt,lgamma,lang_facts)
implicit none
!     SAGARMOY MANDAL  !SAGAR HACK
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


subroutine gauss_dist( mu, sigma, dim, gauss_dist1)
implicit none
!-----------------------------------------------------------------------
!
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
